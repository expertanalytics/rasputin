"""A PLY writer: arrays in, bytes out.

This module is pure. It takes no path, opens no file, imports nothing
first-party but `features`, `io/models.py` and the shared checks of
`io/mesh_checks.py`, and returns `bytes`, which is what keeps the whole format
testable with no filesystem, no raster and no compiled extension: hand it
three arrays and compare the bytes.

FOUR CHOICES ARE ENCODED IN THE BYTES RATHER THAN IN PROSE.

*Text or binary.* Text is this writer's own default, because a person reads
the output; `rasputin mesh` writes binary unless given `--ascii`. Binary is
little-endian, forced by `<`-prefixed dtypes, so the output does not depend on
the host.

*The coordinates are `double`*, and float32 is not an option. At a
UTM 33N easting of 430 000 the float32 step is 2**19 * 2**-24 = 3.1 cm, which
is thirty times coarser than `cli.DEFAULT_SNAP_SPACING`: two nodes the noder
deliberately kept apart would land on the same coordinate. The ASCII body is
written with `repr`, the shortest representation that round trips, for the same
reason -- a fixed number of decimal places loses the millimetre lattice just as
surely as the narrower dtype would.

*`faces` and `edges` are mutually exclusive.* MDAL's caveat is
that a host application expects either a 1D mesh (edges) or a 2D one (faces),
so a single file holding both is readable and may still display as nothing --
a valid file, a silent half-load, and no error message. The signature is the
one place that can be refused, so it is refused here.

*No z is made anywhere below.* The caller supplies an `(N, 3)`
array; a two-column array is a mistake rather than something to fill in,
because a flat mesh produced by a code path that has no idea it is flat is a
wrong answer that looks like a right one.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from itertools import pairwise

import numpy as np
import numpy.typing as npt

from tin_engine.features import EdgeVocabulary
from tin_engine.io.mesh_checks import check_int32, checked_ascii
from tin_engine.io.models import ClassFractions

#: The edge element's third property: the feature mask of `_core.pyi`'s
#: `NodedPslg.edge_properties`, one `uint32` per edge. PLY fixes no name for
#: it -- `x`/`y`/`z`, `vertex_indices`, `vertex1`/`vertex2` are conventional
#: and this is not -- so readers surface it as a per-edge scalar dataset under
#: whatever this says. 0 means *unclassified*, never *wrong*.
FEATURE_MASK = "feature_mask"

#: Vertex indices are `uint32` in `_core`, so they are `uint` on the wire.
_INDEX = "uint"


def write_ply(
    vertices: npt.ArrayLike,
    *,
    faces: npt.ArrayLike | None = None,
    edges: npt.ArrayLike | None = None,
    edge_properties: npt.ArrayLike | None = None,
    ascii: bool = True,
    comments: Sequence[str] = (),
    vocabulary: EdgeVocabulary | None = None,
    face_codes: npt.ArrayLike | None = None,
    face_fractions: ClassFractions | None = None,
) -> bytes:
    """Encode one mesh as a PLY file.

    Args:
        vertices: `(N, 3)` float64 coordinates. z is the caller's.
        faces: `(T, 3)` vertex indices, or None. The 2D mesh.
        edges: `(E, 2)` vertex indices, or None. The 1D mesh.
        edge_properties: `(E,)` feature masks to ride along with `edges`.
        ascii: write the bodies as text instead of packed records.
        comments: header comment lines, in order, each without its keyword.
        vocabulary: what each bit of `edge_properties` means. Written as
            `feature_bit <bit> <name>` comments sorted by bit and a
            `feature_vocabulary <fingerprint>` comment, after `comments`.
        face_codes: `(T,)` land-cover codes, written as the face property
            `int land_cover_code` after `vertex_indices`.
            The comment naming the code system is the caller's, in `comments`.
        face_fractions: the sparse table of class fractions per face, written
            as the face lists `land_cover_classes` (`list
            uchar uchar`) and `land_cover_fractions` (`list uchar float`)
            between `vertex_indices` and `land_cover_code`.

    Raises:
        ValueError: if `vertices` is not `(N, 3)`; if `faces` and `edges` are
            both given or neither is; if `edge_properties` is given
            without `edges`, or does not have one entry per edge; if a comment
            contains a control character, any of which forges a header line in
            a line-oriented format; if a comment is not ASCII, which the
            header's encoding cannot carry; or if a mask carries a bit
            `vocabulary` does not name; if `face_codes` is given with `edges`,
            has not one entry per face, or does not fit int32; if
            `face_fractions` is given with `edges`, has not one row per face,
            or a row of more than 255 classes.
    """
    points = np.ascontiguousarray(vertices, dtype="<f8")
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"vertices must have shape (N, 3); got {points.shape}")
    if (faces is None) == (edges is None):
        which = "both" if faces is not None else "neither"
        raise ValueError(f"exactly one of faces and edges must be given; got {which}")
    if edge_properties is not None and edges is None:
        raise ValueError("edge_properties needs edges; it has no meaning beside faces")
    if face_codes is not None and faces is None:
        raise ValueError("face_codes needs faces; the edge file has none")
    if face_fractions is not None and faces is None:
        raise ValueError("face_fractions needs faces; the edge file has none")
    if vocabulary is not None:
        # A bit the vocabulary does not name is refused, not written: no
        # reader could tell what it means.
        if edge_properties is not None:
            for mask in np.unique(np.asarray(edge_properties)):
                vocabulary.names(int(mask))
        comments = [
            *comments,
            *(f"feature_bit {bit} {name}" for bit, name in vocabulary.table()),
            f"feature_vocabulary {vocabulary.fingerprint()}",
        ]
    # Any control character forges a header line (\r too, for CRLF-tolerant
    # readers), and --crs is free text: io/mesh_checks.py's gate.
    for comment in comments:
        checked_ascii(comment, "comment")

    blocks = [(_vertex_declaration(len(points)), _vertex_body(points, ascii))]
    if faces is not None:
        blocks.append(
            _face_block(np.ascontiguousarray(faces, dtype="<u4"), face_codes, face_fractions, ascii)
        )
    if edges is not None:
        blocks.append(_edge_block(np.ascontiguousarray(edges, dtype="<u4"), edge_properties, ascii))

    lines = ["ply", f"format {'ascii' if ascii else 'binary_little_endian'} 1.0"]
    lines += [f"comment {comment}" for comment in comments]
    lines += [line for declaration, _ in blocks for line in declaration]
    lines.append("end_header")
    header = ("\n".join(lines) + "\n").encode("ascii")
    return header + b"".join(body for _, body in blocks)


def field_comments(fields: Iterable[tuple[str, str]]) -> list[str]:
    """A file's `(name, value)` fields as header comments, `name value`."""
    return [f"{name} {value}" for name, value in fields]


def _vertex_declaration(count: int) -> list[str]:
    return [f"element vertex {count}", *(f"property double {axis}" for axis in "xyz")]


def _vertex_body(points: npt.NDArray[np.float64], ascii: bool) -> bytes:
    if not ascii:
        return points.tobytes()
    # `repr` and not a format spec: it is the shortest text that reads back as
    # the same double, so the ASCII body is lossless and agrees with the binary
    # one bit for bit. `%.6f` would silently discard the snap lattice.
    return _lines(" ".join(repr(float(c)) for c in point) for point in points)


def _face_block(
    faces: npt.NDArray[np.uint32],
    codes: npt.ArrayLike | None,
    fractions: ClassFractions | None,
    ascii: bool,
) -> tuple[list[str], bytes]:
    """The 2D mesh's block: one variable-length index list per triangle, then
    optionally its class and fraction lists and its land-cover code."""
    faces = faces.reshape(-1, 3)
    declaration = [
        f"element face {len(faces)}",
        f"property list uchar {_INDEX} vertex_indices",
    ]
    columns = [np.full(len(faces), 3), faces]
    layout: list[tuple[str, str] | tuple[str, str, int]] = [("n", "u1"), ("v", "<u4", 3)]
    if fractions is not None:
        if fractions.triangles != len(faces):
            raise ValueError(
                f"face_fractions has {fractions.triangles} rows for {len(faces)} faces"
            )
        declaration += [
            "property list uchar uchar land_cover_classes",
            "property list uchar float land_cover_fractions",
        ]
    if codes is not None:
        values = np.asarray(codes).reshape(-1)
        if len(values) != len(faces):
            raise ValueError(f"face_codes has {len(values)} entries for {len(faces)} faces")
        check_int32(values, "face_codes")
        declaration.append("property int land_cover_code")
        columns.append(values)
        layout.append(("c", "<i4"))
    if fractions is not None:
        return declaration, _face_lists(faces, codes, fractions, ascii)
    if ascii:
        table = np.column_stack(columns).astype(np.int64)
        return declaration, _lines(" ".join(map(str, row)) for row in table.tolist())
    records = np.zeros(len(faces), dtype=np.dtype(layout))
    for (name, *_), column in zip(layout, columns, strict=True):
        records[name] = column
    return declaration, records.tobytes()


def _face_lists(
    faces: npt.NDArray[np.uint32], codes: npt.ArrayLike | None, table: ClassFractions, ascii: bool
) -> bytes:
    """The face records with their class and fraction lists: variable length,
    so the binary body is assembled by offset rather than as fixed records."""
    offsets = table.offsets.astype(np.int64)
    counts = np.diff(offsets)
    if counts.max(initial=0) > 255:
        raise ValueError("a face of more than 255 classes does not fit a uchar list count")
    code = None if codes is None else np.asarray(codes).reshape(-1).astype(np.int64)
    if ascii:
        rows = []
        for t, (a, b) in enumerate(pairwise(offsets)):
            words = [3, *faces[t], counts[t], *table.classes[a:b], counts[t]]
            words += [repr(float(f)) for f in table.fractions[a:b]]
            rows.append(" ".join(map(str, words)) + ("" if code is None else f" {code[t]}"))
        return _lines(rows)
    tail = 0 if code is None else 4
    length = 15 + 5 * counts + tail
    start = np.concatenate([[0], np.cumsum(length)[:-1]])
    out = np.zeros(int(length.sum()), dtype=np.uint8)
    out[start] = 3
    out[start[:, None] + 1 + np.arange(12)] = faces.astype("<u4").view(np.uint8).reshape(-1, 12)
    out[start + 13] = counts
    within = np.arange(len(table.classes)) - np.repeat(offsets[:-1], counts)
    out[np.repeat(start + 14, counts) + within] = table.classes
    out[start + 14 + counts] = counts
    at = np.repeat(start + 15 + counts, counts) + 4 * within
    out[at[:, None] + np.arange(4)] = table.fractions.astype("<f4").view(np.uint8).reshape(-1, 4)
    if code is not None:
        at = start + 15 + 5 * counts
        out[at[:, None] + np.arange(4)] = code.astype("<i4").view(np.uint8).reshape(-1, 4)
    return out.tobytes()


def _edge_block(
    edges: npt.NDArray[np.uint32], properties: npt.ArrayLike | None, ascii: bool
) -> tuple[list[str], bytes]:
    """The 1D mesh's block: two endpoints and, optionally, one feature scalar."""
    columns = ["vertex1", "vertex2"]
    values = [edges.reshape(-1, 2)[:, 0], edges.reshape(-1, 2)[:, 1]]
    if properties is not None:
        masks = np.ascontiguousarray(properties, dtype="<u4").reshape(-1)
        if len(masks) != len(values[0]):
            raise ValueError(f"edge_properties has {len(masks)} entries for {len(values[0])} edges")
        columns.append(FEATURE_MASK)
        values.append(masks)
    declaration = [
        f"element edge {len(values[0])}",
        *(f"property {_INDEX} {column}" for column in columns),
    ]
    if ascii:
        return declaration, _lines(
            " ".join(str(int(v)) for v in record) for record in zip(*values, strict=True)
        )
    records = np.zeros(len(values[0]), dtype=np.dtype([(c, "<u4") for c in columns]))
    for column, column_values in zip(columns, values, strict=True):
        records[column] = column_values
    return declaration, records.tobytes()


def _lines(rows: Iterable[str]) -> bytes:
    """One newline-terminated ASCII line per record, and none when there are none."""
    return "".join(f"{row}\n" for row in rows).encode("ascii")
