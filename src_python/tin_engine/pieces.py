"""A cut run's pieces, sliced from its one start triangulation.

`docs/increments/23-basin-scale.md`, *The seam protocol*, steps 1, 3 and 5.
A piece is a component of the start
triangulation across every edge that is not a seam (`landcover.regions`
with only the seams blocking); its id `(j, i, k)` is its cell's row and
column and its rank in that cell by lowest start-triangle index, a function
of the input alone. Its start mesh is the start triangulation's triangles
with its label, in their order, each seam edge fanned at its seam-pass
points: index arithmetic, no geometry.

Pure NumPy over arrays; the seam pass itself is the caller's callable, so
this module never calls the core.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from tin_engine.io.mesh_index import PieceId, SeamRecord
from tin_engine.landcover import regions

if TYPE_CHECKING:  # the seam pass's result, read here; the call is the caller's
    from tin_engine._core import SeamOutcome


def piece_ids(
    triangles: npt.NDArray[np.int64],
    vertices: npt.NDArray[np.float64],
    seams: npt.NDArray[np.int64],
    cell: Callable[[float, float], tuple[int, int]],
) -> tuple[npt.NDArray[np.int64], list[PieceId]]:
    """Each start triangle's piece, as an index into the ids, which are sorted.
    `cell(x, y)` is the `(j, i)` of the cell holding a point; a component is
    placed by its lowest triangle's centroid, strictly inside its cell."""
    roots = regions(triangles, seams)
    firsts = np.unique(roots)  # each component's lowest triangle, ascending
    centres = vertices[triangles[firsts]].mean(axis=1)
    ranks: dict[tuple[int, int], int] = {}
    ids: list[PieceId] = []
    for x, y in centres.tolist():
        j, i = cell(x, y)
        ids.append((j, i, ranks.get((j, i), 0)))
        ranks[(j, i)] = ids[-1][2] + 1
    order = sorted(range(len(ids)), key=ids.__getitem__)
    position = np.empty(len(ids), dtype=np.int64)
    position[order] = np.arange(len(ids))
    return position[np.searchsorted(firsts, roots)], [ids[k] for k in order]


@dataclass(frozen=True, slots=True)
class Slice:
    """A piece's start mesh, as `refine` takes it, and its seams: `chains[r]`
    is seam edge `r`'s local vertex chain from the edge's first end, and `z`
    the seam record's height of each seam vertex that has one."""

    vertices: npt.NDArray[np.float64]
    triangles: npt.NDArray[np.uint32]
    flags: npt.NDArray[np.uint8]
    edges: npt.NDArray[np.uint32]
    masks: npt.NDArray[np.uint32]
    chains: dict[int, list[int]]
    z: dict[int, float]


def start_slice(
    vertices: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.int64],
    flags: npt.NDArray[np.uint8],
    edges: npt.NDArray[np.uint32],
    masks: npt.NDArray[np.uint32],
    keep: npt.NDArray[np.bool_],
    seam: int,
    seam_pass: Callable[[int], SeamOutcome],
) -> Slice:
    """The start triangles `keep` selects, in order, vertices renumbered in
    their global order (so the slice's order is the whole's, restricted), and
    each seam edge `(a, b)` of triangle `(a, b, c)` with seam-pass points
    `p1 .. pk` made `(a, p1, c), .., (pk, b, c)`; the first replaces the
    triangle, the rest go at the end, the seam edges taken in edge order.
    A vertex at the end of several seam edges takes the z of the lowest (N17)."""
    tri, cons = triangles[keep], flags[keep].astype(np.int64)
    used = np.unique(tri)
    n = len(vertices)

    def key(a: npt.NDArray[np.int64], b: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
        return np.minimum(a, b) * n + np.maximum(a, b)

    held = [key(tri[:, e], tri[:, (e + 1) % 3])[(cons >> e) & 1 == 1] for e in range(3)]
    pairs = edges.astype(np.int64)
    rows = np.flatnonzero(np.isin(key(pairs[:, 0], pairs[:, 1]), np.concatenate(held)))
    tris = np.searchsorted(used, tri).tolist()
    bits = cons.tolist()
    xy = [(x, y) for x, y in vertices[used].tolist()]
    where = {frozenset((t[e], t[(e + 1) % 3])): k for k, t in enumerate(tris) for e in range(3)}
    out_edges: list[tuple[int, int]] = []
    out_masks: list[int] = []
    chains: dict[int, list[int]] = {}
    heights: dict[int, float] = {}
    for r in rows.tolist():
        a, b = (int(v) for v in np.searchsorted(used, pairs[r]))
        m = int(masks[r])
        if not m & seam:
            out_edges.append((a, b))
            out_masks.append(m)
            continue
        found = seam_pass(r)
        points, zs = np.asarray(found.points).tolist(), np.asarray(found.z).tolist()
        ends = [found.z_a, found.z_b]
        if tuple(found.a) != xy[a]:
            points, zs, ends = points[::-1], zs[::-1], ends[::-1]
        chain = chains[r] = [a, *range(len(xy), len(xy) + len(points)), b]
        xy += [(x, y) for x, y in points]
        for v, z in zip(chain, [ends[0], *zs, ends[1]], strict=True):
            if z is not None:
                heights.setdefault(v, z)
        out_edges += itertools.pairwise(chain)
        out_masks += [m] * (len(chain) - 1)
        t = where[frozenset((a, b))]
        e = next(e for e in range(3) if {tris[t][e], tris[t][(e + 1) % 3]} == {a, b})
        run = chain if tris[t][e] == a else chain[::-1]
        c, last = tris[t][(e + 2) % 3], len(run) - 2
        qc, cp = (bits[t] >> (e + 1) % 3) & 1, (bits[t] >> (e + 2) % 3) & 1
        subs = [[run[i], run[i + 1], c] for i in range(last + 1)]
        sub_bits = [1 | (qc << 1) * (i == last) | (cp << 2) * (i == 0) for i in range(last + 1)]
        tris[t], bits[t] = subs[0], sub_bits[0]
        tris += subs[1:]
        bits += sub_bits[1:]
        where[frozenset((run[-1], c))] = len(tris) - 1 if last else t
    return Slice(
        vertices=np.array(xy, dtype=np.float64).reshape(-1, 2),
        triangles=np.array(tris, dtype=np.uint32).reshape(-1, 3),
        flags=np.array(bits, dtype=np.uint8),
        edges=np.array(out_edges, dtype=np.uint32).reshape(-1, 2),
        masks=np.array(out_masks, dtype=np.uint32),
        chains=chains,
        z=heights,
    )


def seam_records(
    piece: PieceId,
    piece_slice: Slice,
    vertices: npt.NDArray[np.float64],
    z: npt.NDArray[np.float64],
    valid: npt.NDArray[np.bool_],
) -> list[SeamRecord]:
    """What the piece writes along each seam edge, from the refined arrays
    (start vertices keep their indices); a vertex without data has NaN z."""
    return [
        SeamRecord(
            piece=piece,
            edge=r,
            xyz=tuple(
                (
                    float(vertices[v, 0]),
                    float(vertices[v, 1]),
                    float(z[v]) if valid[v] else math.nan,
                )
                for v in chain
            ),
        )
        for r, chain in piece_slice.chains.items()
    ]
