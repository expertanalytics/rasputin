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
from tin_engine.io.models import ClassFractions, HeightTable

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
    face_filled: npt.ArrayLike | None = None,
    face_heights: HeightTable | None = None,
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
        face_filled: `(T,)` true where a code was filled, written with
            `face_codes` as the face property `uchar land_cover_filled`.
        face_fractions: the sparse table of class fractions per face, written
            as the face lists `land_cover_classes` (`list
            uchar int`) and `land_cover_fractions` (`list uchar float`)
            between `vertex_indices` and `land_cover_code`.
        face_heights: the DEM height shares per face, written last as the
            face properties `dem_height_shares` (`list uchar float`),
            `int dem_height_first` and `float dem_nodata_share`.

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
            or a row of more than 255 classes; or if `face_heights` has not
            one row per face, or a row of more than 255 shares.
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
    if (face_fractions is not None or face_heights is not None) and faces is None:
        raise ValueError("face_fractions and face_heights need faces; the edge file has none")
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
            _face_block(
                np.ascontiguousarray(faces, dtype="<u4"), face_codes, face_fractions, ascii,
                face_filled, face_heights,
            )
        )  # fmt: skip
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


#: A face property: its declaration, and either fixed values `(T, k)` or a
#: list `(offsets, values)`, with the little-endian dtype it is written in.
_Column = tuple[
    str, npt.NDArray[np.generic] | tuple[npt.NDArray[np.int64], npt.NDArray[np.generic]], str
]


def _face_block(
    faces: npt.NDArray[np.uint32],
    codes: npt.ArrayLike | None,
    fractions: ClassFractions | None,
    ascii: bool,
    filled: npt.ArrayLike | None = None,
    heights: HeightTable | None = None,
) -> tuple[list[str], bytes]:
    """The 2D mesh's block: one index list per triangle, then optionally its
    class and fraction lists, its land-cover code and mark, and its DEM height
    shares, first class and NoData share."""
    faces = faces.reshape(-1, 3)
    corners = (np.arange(len(faces) + 1, dtype=np.int64) * 3, faces.reshape(-1))
    columns: list[_Column] = [(f"list uchar {_INDEX} vertex_indices", corners, "<u4")]
    if fractions is not None:
        if fractions.triangles != len(faces):
            raise ValueError(
                f"face_fractions has {fractions.triangles} rows for {len(faces)} faces"
            )
        offsets = fractions.offsets.astype(np.int64)
        columns += [
            ("list uchar int land_cover_classes", (offsets, fractions.classes), "<i4"),
            ("list uchar float land_cover_fractions", (offsets, fractions.fractions), "<f4"),
        ]
    if codes is not None:
        values = np.asarray(codes).reshape(-1)
        if len(values) != len(faces):
            raise ValueError(f"face_codes has {len(values)} entries for {len(faces)} faces")
        check_int32(values, "face_codes")
        columns.append(("int land_cover_code", values, "<i4"))
        if filled is not None:
            marks = np.asarray(filled).reshape(-1).astype(np.uint8)
            columns.append(("uchar land_cover_filled", marks, "u1"))
    if heights is not None:
        if len(heights.first) != len(faces):
            raise ValueError(f"face_heights has {len(heights.first)} rows for {len(faces)} faces")
        offsets = heights.offsets.astype(np.int64)
        columns += [
            ("list uchar float dem_height_shares", (offsets, heights.shares), "<f4"),
            ("int dem_height_first", heights.first, "<i4"),
            ("float dem_nodata_share", heights.nodata, "<f4"),
        ]
    declaration = [f"element face {len(faces)}", *(f"property {d}" for d, _, _ in columns)]
    return declaration, _records(len(faces), columns, ascii)


def _records(count: int, columns: list[_Column], ascii: bool) -> bytes:
    """`count` face records of `columns`: variable length where a column is a
    list, so the binary body is assembled by offset rather than as fixed records."""
    if ascii:
        words: list[list[str]] = [[] for _ in range(count)]
        for _, values, dtype in columns:
            text = repr if dtype.endswith("f4") else str
            if isinstance(values, tuple):
                offsets, flat = values
                if np.diff(offsets).max(initial=0) > 255:
                    raise ValueError(
                        "a face list of more than 255 values does not fit a uchar count"
                    )
                items = [text(v) for v in flat.tolist()]
                for t, (a, b) in enumerate(pairwise(offsets.tolist())):
                    words[t] += [str(b - a), *items[a:b]]
            else:
                for t, v in enumerate(values.tolist()):
                    words[t].append(text(v))
        return _lines(" ".join(w) for w in words)
    sizes = []
    for _, values, dtype in columns:
        size = np.dtype(dtype).itemsize
        if isinstance(values, tuple):
            counts = np.diff(values[0])
            if counts.max(initial=0) > 255:
                raise ValueError("a face list of more than 255 values does not fit a uchar count")
            sizes.append(1 + size * counts)
        else:
            sizes.append(np.full(count, size, dtype=np.int64))
    length = np.sum(sizes, axis=0) if sizes else np.zeros(count, np.int64)
    at = np.concatenate([[0], np.cumsum(length)[:-1]]).astype(np.int64)
    out = np.zeros(int(length.sum()), dtype=np.uint8)
    for (_, values, dtype), step in zip(columns, sizes, strict=True):
        width = np.dtype(dtype).itemsize
        if isinstance(values, tuple):
            offsets, flat = values
            counts = np.diff(offsets)
            out[at] = counts
            within = np.arange(len(flat)) - np.repeat(offsets[:-1], counts)
            place = np.repeat(at + 1, counts) + width * within
        else:
            flat, place = values, at
        raw = np.ascontiguousarray(flat, dtype=dtype).view(np.uint8).reshape(-1, width)
        out[place[:, None] + np.arange(width)] = raw
        at = at + step
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
