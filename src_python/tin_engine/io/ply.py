"""A PLY writer: arrays in, bytes out.

Increment 10, ruling 6. This module is pure. It takes no path, opens no file,
imports nothing first-party but `features` (increment 13), and returns
`bytes`, which is what keeps the whole format testable with no filesystem, no
raster and no compiled extension: hand it three arrays and compare the bytes.

THREE RULINGS ARE ENCODED IN THE BYTES RATHER THAN IN PROSE.

*Ruling 1* made binary little-endian the default; increment 13's U2 (a) made
text the default instead, because the user reads the output. Binary stays
behind the flag, and its endianness is forced by `<`-prefixed dtypes so the
output does not depend on the host.

*Ruling 2* makes the coordinates `double`, and float32 is not an option. At a
UTM 33N easting of 430 000 the float32 step is 2**19 * 2**-24 = 3.1 cm, which
is thirty times coarser than `cli.DEFAULT_SNAP_SPACING`: two nodes the noder
deliberately kept apart would land on the same coordinate. The ASCII body is
written with `repr`, the shortest representation that round trips, for the same
reason -- a fixed number of decimal places loses the millimetre lattice just as
surely as the narrower dtype would.

*Ruling 3* is why `faces` and `edges` are mutually exclusive. MDAL's caveat is
that a host application expects either a 1D mesh (edges) or a 2D one (faces),
so a single file holding both is readable and may still display as nothing --
a valid file, a silent half-load, and no error message. The signature is the
one place that can be refused, so it is refused here.

*Ruling 4* is why there is no z anywhere below. The caller supplies an `(N, 3)`
array; a two-column array is a mistake rather than something to fill in,
because a flat mesh produced by a code path that has no idea it is flat is a
wrong answer that looks like a right one.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import numpy.typing as npt

from tin_engine.features import EdgeVocabulary

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
) -> bytes:
    """Encode one mesh as a PLY file.

    Args:
        vertices: `(N, 3)` float64 coordinates. z is the caller's (ruling 4).
        faces: `(T, 3)` vertex indices, or None. The 2D mesh.
        edges: `(E, 2)` vertex indices, or None. The 1D mesh.
        edge_properties: `(E,)` feature masks to ride along with `edges`.
        ascii: write the bodies as text instead of packed records.
        comments: header comment lines, in order, each without its keyword.
        vocabulary: what each bit of `edge_properties` means. Written as
            `feature_bit <bit> <name>` comments sorted by bit and a
            `feature_vocabulary <fingerprint>` comment (increment 13, ruling
            9), after `comments`.

    Raises:
        ValueError: if `vertices` is not `(N, 3)`; if `faces` and `edges` are
            both given or neither is (ruling 3); if `edge_properties` is given
            without `edges`, or does not have one entry per edge; if a comment
            contains a control character, any of which forges a header line in
            a line-oriented format; if a comment is not ASCII, which the
            header's encoding cannot carry; or if a mask carries a bit
            `vocabulary` does not name.
    """
    points = np.ascontiguousarray(vertices, dtype="<f8")
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"vertices must have shape (N, 3); got {points.shape}")
    if (faces is None) == (edges is None):
        which = "both" if faces is not None else "neither"
        raise ValueError(f"exactly one of faces and edges must be given; got {which}")
    if edge_properties is not None and edges is None:
        raise ValueError("edge_properties needs edges; it has no meaning beside faces")
    if vocabulary is not None:
        # Increment 7's mechanism 3: a bit nobody names is refused, not written.
        if edge_properties is not None:
            for mask in np.unique(np.asarray(edge_properties)):
                vocabulary.names(int(mask))
        table = sorted((prop.bit, prop.name) for prop in vocabulary.properties)
        comments = [
            *comments,
            *(f"feature_bit {bit} {name}" for bit, name in table),
            f"feature_vocabulary {vocabulary.fingerprint()}",
        ]
    for comment in comments:
        # Any control character, not just \n. A PLY header is line-oriented and
        # \r terminates a line for every CRLF-tolerant reader, which is most of
        # them -- so an unguarded \r forges a header line exactly as \n would.
        # Measured before this guard was widened: `--crs "x\rcomment forged"`
        # wrote that second line into the header and exited 0.
        bad = next((ch for ch in comment if ch < " " or ch == "\x7f"), None)
        if bad is not None:
            raise ValueError(
                f"a comment may not contain control characters; got {bad!r}"
            )
        # ASCII is the header's encoding, so a non-ASCII comment cannot be
        # written. --crs is unvalidated free text by ruling 5 and a degree sign
        # in a projection string is ordinary, so this is a refusal a caller
        # meets, not an internal invariant: it must be the documented
        # ValueError and not a UnicodeEncodeError escaping from the encode
        # below.
        if not comment.isascii():
            bad = next(ch for ch in comment if not ch.isascii())
            raise ValueError(
                f"a comment must be ASCII; got {bad!r} in {comment!r}"
            )

    blocks = [(_vertex_declaration(len(points)), _vertex_body(points, ascii))]
    if faces is not None:
        blocks.append(_face_block(np.ascontiguousarray(faces, dtype="<u4"), ascii))
    if edges is not None:
        blocks.append(_edge_block(np.ascontiguousarray(edges, dtype="<u4"), edge_properties, ascii))

    lines = ["ply", f"format {'ascii' if ascii else 'binary_little_endian'} 1.0"]
    lines += [f"comment {comment}" for comment in comments]
    lines += [line for declaration, _ in blocks for line in declaration]
    lines.append("end_header")
    header = ("\n".join(lines) + "\n").encode("ascii")
    return header + b"".join(body for _, body in blocks)


def _vertex_declaration(count: int) -> list[str]:
    return [f"element vertex {count}", *(f"property double {axis}" for axis in "xyz")]


def _vertex_body(points: npt.NDArray[np.float64], ascii: bool) -> bytes:
    if not ascii:
        return points.tobytes()
    # `repr` and not a format spec: it is the shortest text that reads back as
    # the same double, so the ASCII body is lossless and agrees with the binary
    # one bit for bit. `%.6f` would silently discard the snap lattice.
    return _lines(" ".join(repr(float(c)) for c in point) for point in points)


def _face_block(faces: npt.NDArray[np.uint32], ascii: bool) -> tuple[list[str], bytes]:
    """The 2D mesh's block: one variable-length index list per triangle."""
    declaration = [
        f"element face {len(faces)}",
        f"property list uchar {_INDEX} vertex_indices",
    ]
    if ascii:
        return declaration, _lines(
            " ".join(str(int(v)) for v in (len(face), *face)) for face in faces
        )
    records = np.zeros(len(faces), dtype=np.dtype([("n", "u1"), ("v", "<u4", (3,))]))
    records["n"] = 3
    records["v"] = faces.reshape(-1, 3)
    return declaration, records.tobytes()


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
