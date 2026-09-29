"""A land-cover class per triangle: components first, one point per component.

Increment 16c (`docs/increments/16c-landcover-labels.md`, R1). Pure: numpy and
shapely in, arrays and counts out; nothing first-party, never `_core`, never a
path.

The triangles are split into connected components across every edge that is
not a constraint edge (the spread of Shewchuk's Triangle `-A`, blocked by
segments). Each component is then looked up once, at the incentre of its
triangle with the largest inradius: that point is at least `r` from every
constraint edge, so when `r` exceeds `margin` (which bounds how far the noder
moved an input boundary) it lies strictly on one side of every input
boundary. A component whose largest `r` is at most `margin` is still labelled
by its point, and counted `thin`.

Overlaps (Default D2): the polygon of smallest area wins, ties to the smaller
code, so the result does not depend on the polygons' order.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry.base import BaseGeometry


@dataclass(frozen=True)
class CoverLabels:
    """The codes, and the counts of R3's stderr line."""

    codes: npt.NDArray[np.int32]
    regions: int
    outside: int
    overlapped: int
    thin: int


def regions(triangles: npt.ArrayLike, edges: npt.ArrayLike) -> npt.NDArray[np.int64]:
    """Each triangle's component across unconstrained edges, as the smallest
    triangle index in that component."""
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    cut = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    n = int(max(tri.max(initial=-1), cut.max(initial=-1))) + 1
    sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    keys = _keys(sides, n)
    owner = np.tile(np.arange(len(tri), dtype=np.int64), 3)
    free = ~np.isin(keys, _keys(cut, n))
    order = np.argsort(keys[free], kind="stable")
    keys, owner = keys[free][order], owner[free][order]
    # A manifold mesh: an interior edge's key appears exactly twice.
    pair = np.flatnonzero(keys[1:] == keys[:-1])
    u, v = owner[pair], owner[pair + 1]
    parent = np.arange(len(tri), dtype=np.int64)
    while True:
        ru, rv = parent[u], parent[v]
        if np.array_equal(ru, rv):
            return parent
        # Hook the larger root onto the smaller: labels only fall, no cycle.
        np.minimum.at(parent, np.maximum(ru, rv), np.minimum(ru, rv))
        while not np.array_equal(parent, jumped := parent[parent]):
            parent = jumped


def label_triangles(
    vertices: npt.ArrayLike,
    triangles: npt.ArrayLike,
    edges: npt.ArrayLike,
    *,
    polygons: Sequence[tuple[BaseGeometry, int]],
    margin: float,
) -> CoverLabels:
    """A code per triangle from `(polygon, code)` pairs; 0 in no polygon.

    Args:
        vertices: `(N, 2)` or `(N, 3)`; only x and y are used.
        triangles: `(T, 3)` vertex indices.
        edges: `(E, 2)` constraint edges; they block the spread.
        polygons: the coded polygons, in the vertices' CRS.
        margin: how far a constraint edge may lie from its input boundary.
    """
    xy = np.asarray(vertices, dtype=np.float64)[:, :2]
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    ids = regions(tri, edges)
    a, b, c = xy[tri[:, 0]], xy[tri[:, 1]], xy[tri[:, 2]]
    la, lb, lc = (np.hypot(*(q - p).T) for p, q in ((b, c), (c, a), (a, b)))
    perimeter = la + lb + lc
    cross = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    safe = np.where(perimeter > 0, perimeter, 1.0)
    r = np.where(perimeter > 0, np.abs(cross) / safe, 0.0)
    centre = (la[:, None] * a + lb[:, None] * b + lc[:, None] * c) / safe[:, None]
    # Per component, the triangle of largest r, ties to the lowest index.
    order = np.lexsort((np.arange(len(tri)), -r, ids))
    best = order[np.unique(ids[order], return_index=True)[1]]
    comp_codes, hits = _lookup(centre[best], polygons)
    by_id = np.zeros(len(tri), dtype=np.int32)
    by_id[ids[best]] = comp_codes
    codes = by_id[ids]
    return CoverLabels(
        codes=codes,
        regions=len(best),
        outside=int((hits == 0).sum()),
        overlapped=int((hits > 1).sum()),
        thin=int((r[best] <= margin).sum()),
    )


def _keys(pairs: npt.NDArray[np.int64], n: int) -> npt.NDArray[np.int64]:
    """Undirected edge keys, `min * n + max`."""
    return np.minimum(pairs[:, 0], pairs[:, 1]) * n + np.maximum(pairs[:, 0], pairs[:, 1])


def _lookup(
    points: npt.NDArray[np.float64], polygons: Sequence[tuple[BaseGeometry, int]]
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    """Each point's code by D2, and how many polygons it lies in."""
    codes = np.zeros(len(points), dtype=np.int64)
    hits = np.zeros(len(points), dtype=np.int64)
    if not polygons or not len(points):
        return codes, hits
    tree = shapely.STRtree([p for p, _ in polygons])
    found, which = tree.query(shapely.points(points), predicate="intersects")
    rank = np.array([(p.area, code) for p, code in polygons], dtype=[("a", "f8"), ("c", "i8")])
    # Sort the hits by point, then area, then code: each point's first hit wins.
    order = np.lexsort((rank["c"][which], rank["a"][which], found))
    found, which = found[order], which[order]
    np.add.at(hits, found, 1)
    head = np.unique(found, return_index=True)[1]
    codes[found[head]] = rank["c"][which[head]]
    return codes, hits
