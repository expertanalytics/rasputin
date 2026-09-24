"""A legacy VTK writer: arrays and a vocabulary in, bytes out.

Increment 13 (`docs/increments/13-bundled-mesh.md`). One file carries the
surface, the constraint edges and each edge's feature bits, for ParaView. Like
`ply.py` this module is pure: no path, no file, no `_core`. Its one first-party
import is `features`, which imports only hashlib and pydantic.

THE RULINGS THE BYTES ENCODE.

*Ruling 2*: points are `double`, written with `repr` in text and `>f8` packed.

*Ruling 3*: ASCII by default. Binary legacy VTK is big-endian by the format's
definition, so every packed dtype here is `>`-prefixed -- the PLY writer's
opposite.

*Rulings 4 and 5*: every `LINES` cell precedes every `POLYGONS` cell, in the
file and in every cell array, and every cell array covers every cell, with 0
on triangles. VTK orders cells lines-first whatever the file says and matches
cell data by that index, so any other order gives values to the wrong cells in
silence; an array short of the cell count empties the dataset.

*Ruling 6*: the vocabulary travels as a `(bit, name)` table and a fingerprint
in the dataset's `FieldData`, and each property set on at least one edge gets
a 0/1 cell array named after it. Every mask goes through `vocabulary.names()`,
which refuses a bit the vocabulary does not name.

*Ruling 7*: every string goes through one encoder, which refuses control
characters and non-ASCII, and in text escapes `%` then space as VTK does.
"""

from __future__ import annotations

import re
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

from tin_engine.features import EdgeVocabulary

#: The names this module writes into `FieldData` itself.
RESERVED = frozenset({"feature_bits", "feature_names", "feature_vocabulary"})

_FIELD_NAME = re.compile(r"[a-z][a-z0-9_]*")

#: Legacy type names onto big-endian numpy dtypes.
_DTYPE = {"double": ">f8", "int": ">i4", "unsigned_int": ">u4", "unsigned_char": "u1"}


def write_vtk(
    vertices: npt.ArrayLike,
    *,
    triangles: npt.ArrayLike,
    edges: npt.ArrayLike,
    edge_masks: npt.ArrayLike,
    vocabulary: EdgeVocabulary,
    fields: Sequence[tuple[str, str]] = (),
    binary: bool = False,
) -> bytes:
    """Encode one mesh and its constraint edges as legacy VTK 4.2 PolyData.

    Args:
        vertices: `(N, 3)` float64 coordinates. z is the caller's.
        triangles: `(T, 3)` vertex indices.
        edges: `(E, 2)` vertex indices; may be empty.
        edge_masks: `(E,)` feature masks, one per edge.
        vocabulary: what each bit of a mask means.
        fields: extra dataset strings, such as `("crs", "EPSG:25833")`.
        binary: write packed big-endian records instead of text.

    Raises:
        ValueError: if `vertices` is not `(N, 3)`; if `edge_masks` does not
            have one entry per edge; if a mask carries a bit the vocabulary
            does not name; if a field name is outside `^[a-z][a-z0-9_]*$` or
            reserved; or if a string is not ASCII or has a control character.
    """
    points = np.asarray(vertices, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"vertices must have shape (N, 3); got {points.shape}")
    lines = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    polygons = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    masks = np.asarray(edge_masks, dtype=np.uint32).reshape(-1)
    if len(masks) != len(lines):
        raise ValueError(f"edge_masks has {len(masks)} entries for {len(lines)} edges")
    for name, _ in fields:
        if not _FIELD_NAME.fullmatch(name) or name in RESERVED:
            raise ValueError(f"field name {name!r} is reserved or not ^[a-z][a-z0-9_]*$")

    table = sorted((prop.bit, prop.name) for prop in vocabulary.properties)
    used = {name for mask in np.unique(masks) for name in vocabulary.names(int(mask))}
    cell_masks = np.concatenate([masks, np.zeros(len(polygons), dtype=np.uint32)])

    dataset = [
        _numeric("feature_bits", np.array([b for b, _ in table]), "unsigned_int", binary),
        _strings("feature_names", [n for _, n in table], binary),
        _strings("feature_vocabulary", [vocabulary.fingerprint()], binary),
        *(_strings(name, [value], binary) for name, value in fields),
    ]
    features = [
        _numeric(name, (cell_masks >> bit) & 1, "unsigned_char", binary)
        for bit, name in table
        if name in used
    ]

    out = [
        b"# vtk DataFile Version 4.2\nrasputin mesh\n",
        b"BINARY\n" if binary else b"ASCII\n",
        b"DATASET POLYDATA\n",
        f"FIELD FieldData {len(dataset)}\n".encode("ascii"),
        *dataset,
        f"POINTS {len(points)} double\n".encode("ascii"),
        _body(points, "double", binary),
        _cells("LINES", lines, binary),
        _cells("POLYGONS", polygons, binary),
        f"CELL_DATA {len(cell_masks)}\n".encode("ascii"),
        b"SCALARS feature_mask unsigned_int 1\nLOOKUP_TABLE default\n",
        _body(cell_masks.reshape(-1, 1), "unsigned_int", binary),
    ]
    if features:
        out += [f"FIELD features {len(features)}\n".encode("ascii"), *features]
    return b"".join(out)


def _body(rows: npt.NDArray[np.generic], type_name: str, binary: bool) -> bytes:
    """One block of numbers: packed big-endian, or one text line per row."""
    if binary:
        return rows.astype(_DTYPE[type_name]).tobytes() + b"\n"
    # `str` of a Python float is its `repr`: the shortest text that round trips.
    return "".join(" ".join(map(str, row)) + "\n" for row in rows.tolist()).encode("ascii")


def _cells(keyword: str, cells: npt.NDArray[np.int64], binary: bool) -> bytes:
    """A cell list: each cell is its point count, then its point ids."""
    flat = np.column_stack([np.full(len(cells), cells.shape[1]), cells])
    header = f"{keyword} {len(cells)} {flat.size}\n".encode("ascii")
    return header + _body(flat, "int", binary)


def _numeric(name: str, values: npt.NDArray[np.generic], type_name: str, binary: bool) -> bytes:
    header = f"{name} 1 {len(values)} {type_name}\n".encode("ascii")
    return header + _body(values.reshape(-1, 1), type_name, binary)


def _strings(name: str, values: Sequence[str], binary: bool) -> bytes:
    """A string array: text lines with `%XX` escapes, or length-prefixed bytes."""
    encoded = [_checked(value) for value in values]
    header = f"{name} 1 {len(values)} string\n".encode("ascii")
    if binary:
        return header + b"".join(_prefix(len(v)) + v for v in encoded) + b"\n"
    return header + b"".join(v.replace(b"%", b"%25").replace(b" ", b"%20") + b"\n" for v in encoded)


def _checked(value: str) -> bytes:
    """The one string gate: no control character forges a line, and ASCII only."""
    bad = next((ch for ch in value if ch < " " or ch == "\x7f"), None)
    if bad is not None:
        raise ValueError(f"a string may not contain a control character; got {bad!r}")
    if not value.isascii():
        bad = next(ch for ch in value if not ch.isascii())
        raise ValueError(f"a string must be ASCII; got {bad!r} in {value!r}")
    return value.encode("ascii")


def _prefix(length: int) -> bytes:
    """VTK's binary string length: the top two bits say how wide the length is."""
    for width, tag in ((1, 0b11), (2, 0b10), (4, 0b01)):
        if length < 1 << (8 * width - 2):
            return (tag << (8 * width - 2) | length).to_bytes(width, "big")
    return length.to_bytes(8, "big")
