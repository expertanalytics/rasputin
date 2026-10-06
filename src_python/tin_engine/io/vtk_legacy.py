"""A legacy VTK writer: arrays and a vocabulary in, bytes out.

Increment 13 (`docs/increments/13-bundled-mesh.md`). One file carries the
surface, the constraint edges and each edge's feature bits, for ParaView. Like
`ply.py` this module is pure: no path, no file, no `_core`. Its first-party
imports are `features` and `io/mesh_checks.py`, both pure.

THE RULINGS THE BYTES ENCODE.

*Ruling 2*: points are `double`, written with `repr` in text and `>f8` packed.

*Ruling 3*: ASCII by default. Binary legacy VTK is big-endian by the format's
definition, so every packed dtype here is `>`-prefixed -- the PLY writer's
opposite.

*Rulings 4 and 5*: every `LINES` cell precedes every `POLYGONS` cell, in the
file and in every cell array, and every cell array covers every cell, with 0
on triangles. VTK orders cells lines-first whatever the file says and matches
cell data by that index, so any other order gives values to the wrong cells in
silence; an array short of the cell count empties the dataset. With no
constraint edges the `LINES` block is left out (ruling 4 as revised on
2026-09-28): vtk 9.7.1 rejects `LINES 0 0`, and a file without the block loads.

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
from tin_engine.io.mesh_checks import check_int32, checked_ascii

#: The names this module writes into `FieldData` itself.
#: ``elevation`` is the point array of heights (increment 12), so no dataset
#: string may take the name: ParaView's Color By would offer both.
#: ``land_cover_codes`` names the code system of the cell array
#: ``land_cover_code`` (increment 16c, R3).
RESERVED = frozenset(
    {"feature_bits", "feature_names", "feature_vocabulary", "elevation", "land_cover_codes"}
)

#: The cell array of land-cover codes. Not ``land_cover``: that is 16b's bit 7.
LAND_COVER_CODE = "land_cover_code"

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
    triangle_codes: npt.ArrayLike | None = None,
    land_cover_codes: str = "",
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
        triangle_codes: `(T,)` land-cover codes (increment 16c, R3), written
            as the `int` cell array `land_cover_code`, 0 on every line.
        land_cover_codes: what the codes are, written as the dataset string
            `land_cover_codes` when non-empty and `triangle_codes` is given.

    Raises:
        ValueError: if `vertices` is not `(N, 3)`; if `edge_masks` does not
            have one entry per edge; if a mask carries a bit the vocabulary
            does not name; if a field name is outside `^[a-z][a-z0-9_]*$` or
            reserved; if a string is not ASCII or has a control character; or
            if `triangle_codes` has not one entry per triangle, does not fit
            int32, or comes with a vocabulary naming `land_cover_code`.
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

    table = vocabulary.table()
    codes = None
    if triangle_codes is not None:
        codes = np.asarray(triangle_codes).reshape(-1)
        if len(codes) != len(polygons):
            raise ValueError(f"triangle_codes has {len(codes)} entries, {len(polygons)} triangles")
        check_int32(codes, "triangle_codes")
        if any(name == LAND_COVER_CODE for _, name in table):
            raise ValueError(f"the vocabulary names {LAND_COVER_CODE!r}, the codes' array")
        if land_cover_codes:
            fields = (*fields, ("land_cover_codes", land_cover_codes))
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
    if codes is not None:
        # Ruling 5: every cell array covers every cell; the lines get 0.
        cell_codes = np.concatenate([np.zeros(len(lines), dtype=np.int64), codes])
        features.append(_numeric(LAND_COVER_CODE, cell_codes, "int", binary))

    out = [
        b"# vtk DataFile Version 4.2\nrasputin mesh\n",
        b"BINARY\n" if binary else b"ASCII\n",
        b"DATASET POLYDATA\n",
        f"FIELD FieldData {len(dataset)}\n".encode("ascii"),
        *dataset,
        f"POINTS {len(points)} double\n".encode("ascii"),
        _body(points, "double", binary),
        # Ruling 4 as revised: vtk 9.7.1 rejects `LINES 0 0`, so no edges, no block.
        *([_cells("LINES", lines, binary)] if len(lines) else []),
        _cells("POLYGONS", polygons, binary),
        # z again, as a point array: Color By -> elevation then shows the heights.
        f"POINT_DATA {len(points)}\n".encode("ascii"),
        b"SCALARS elevation double 1\nLOOKUP_TABLE default\n",
        _body(points[:, 2:], "double", binary),
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
    encoded = [checked_ascii(value, "string") for value in values]
    header = f"{name} 1 {len(values)} string\n".encode("ascii")
    if binary:
        return header + b"".join(_prefix(len(v)) + v for v in encoded) + b"\n"
    return header + b"".join(v.replace(b"%", b"%25").replace(b" ", b"%20") + b"\n" for v in encoded)


def _prefix(length: int) -> bytes:
    """VTK's binary string length: the top two bits say how wide the length is."""
    for width, tag in ((1, 0b11), (2, 0b10), (4, 0b01)):
        if length < 1 << (8 * width - 2):
            return (tag << (8 * width - 2) | length).to_bytes(width, "big")
    return length.to_bytes(8, "big")
