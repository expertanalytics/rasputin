"""A land-cover class per triangle: components first, one point per component.

`docs/increments/16c-landcover-labels.md`, *The blueprint: data flow and boundaries*.

Pure: numpy and shapely in, arrays and counts out; never `_core`, never a
path.

The triangles are split into connected components across every edge that is
not a constraint edge (the spread of Shewchuk's Triangle `-A`, blocked by
segments). Each component is then looked up once, at the incentre of its
triangle with the largest inradius: that point is at least `r` from every
constraint edge, so when `r` exceeds `margin` (which bounds how far the noder
moved an input boundary) it lies strictly on one side of every input
boundary. A component whose largest `r` is at most `margin` is still labelled
by its point, and counted `thin`.

Overlaps: the polygon of smallest area wins, ties to the smaller
code, so the result does not depend on the polygons' order.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry.base import BaseGeometry

#: Shared border lengths are compared rounded to this many decimals, so a tie
#: stays a tie (the land-cover refill's rounding).
DECIMALS = 6


@dataclass(frozen=True)
class CoverLabels:
    """The codes, and the counts the CLI prints on stderr."""

    codes: npt.NDArray[np.int32]
    regions: int
    outside: int
    overlapped: int
    thin: int
    rest: int


def regions(triangles: npt.ArrayLike, edges: npt.ArrayLike) -> npt.NDArray[np.int64]:
    """Each triangle's component across unconstrained edges, as the smallest
    triangle index in that component."""
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    cut = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    n = int(max(tri.max(initial=-1), cut.max(initial=-1))) + 1
    sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    keys = _keys(sides, n)
    # Any sort will do: the components and their labels do not depend on the
    # order of the pairs. Sides are three blocks of T, so the owner is i % T.
    order = np.argsort(keys)
    keys, owner = keys[order], order % max(len(tri), 1)
    free = ~_member(keys, np.unique(_keys(cut, n)))
    keys, owner = keys[free], owner[free]
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
    rest: Sequence[tuple[BaseGeometry, int]] = (),
) -> CoverLabels:
    """A code per triangle from `(polygon, code)` pairs; 0 in no polygon.

    A triangle whose region took 0 then takes the code of the `rest` polygon
    holding its own label point, by the same rule; per triangle, since a
    region the borders leave uncoded may hold several classes of the map.

    Args:
        vertices: `(N, 2)` or `(N, 3)`; only x and y are used.
        triangles: `(T, 3)` vertex indices.
        edges: `(E, 2)` constraint edges; they block the spread.
        polygons: the coded polygons, in the vertices' CRS.
        margin: how far a constraint edge may lie from its input boundary.
        rest: the whole map's coded polygons, for the triangles left 0.
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
    best = _best(ids, r)
    comp_codes, hits = _lookup(centre[best], polygons)
    by_id = np.zeros(len(tri), dtype=np.int32)
    by_id[ids[best]] = comp_codes
    codes = by_id[ids]
    zero = np.flatnonzero(codes == 0)
    from_rest = 0
    if len(rest) and len(zero):
        found = _lookup(centre[zero], rest)[0]
        codes[zero] = found
        from_rest = int((found != 0).sum())
    return CoverLabels(
        codes=codes,
        regions=len(best),
        outside=int((hits == 0).sum()),
        overlapped=int((hits > 1).sum()),
        thin=int((r[best] <= margin).sum()),
        rest=from_rest,
    )


def _keys(pairs: npt.NDArray[np.int64], n: int) -> npt.NDArray[np.int64]:
    """Undirected edge keys, `min * n + max`."""
    return np.minimum(pairs[:, 0], pairs[:, 1]) * n + np.maximum(pairs[:, 0], pairs[:, 1])


def _member(values: npt.NDArray[np.int64], table: npt.NDArray[np.int64]) -> npt.NDArray[np.bool_]:
    """True where a value is in `table`, which is sorted and unique."""
    if not len(table):
        return np.zeros(len(values), dtype=bool)
    at = np.minimum(np.searchsorted(table, values), len(table) - 1)
    return np.asarray(table[at] == values, dtype=bool)


def _best(ids: npt.NDArray[np.int64], r: npt.NDArray[np.float64]) -> npt.NDArray[np.int64]:
    """Per component, in increasing id, the triangle of largest `r` (finite),
    ties to the lowest index; `ids[i] == i` marks a component's smallest index."""
    t = len(ids)
    top = np.full(t, -np.inf)
    np.maximum.at(top, ids, r)
    low = np.full(t, t, dtype=np.int64)
    at_top = np.flatnonzero(r == top[ids])
    np.minimum.at(low, ids[at_top], at_top)
    return low[np.flatnonzero(ids == np.arange(t))]


def _lookup(
    points: npt.NDArray[np.float64], polygons: Sequence[tuple[BaseGeometry, int]]
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    """Each point's code (the smallest containing polygon wins, ties to the
    smaller code), and how many polygons it lies in."""
    codes = np.zeros(len(points), dtype=np.int64)
    hits = np.zeros(len(points), dtype=np.int64)
    if not polygons or not len(points):
        return codes, hits
    # The polygons are the query, so GEOS prepares each one for its own query
    # and frees it after; the caller's geometries are left unprepared.
    tree = shapely.STRtree(shapely.points(points))
    which, found = tree.query([p for p, _ in polygons], predicate="intersects")
    rank = np.array([(p.area, code) for p, code in polygons], dtype=[("a", "f8"), ("c", "i8")])
    # Sort the hits by point, then area, then code: each point's first hit wins.
    order = np.lexsort((rank["c"][which], rank["a"][which], found))
    found, which = found[order], which[order]
    np.add.at(hits, found, 1)
    head = np.unique(found, return_index=True)[1]
    codes[found[head]] = rank["c"][which[head]]
    return codes, hits


@dataclass(frozen=True)
class Filled:
    """The codes with no 0 left, where they came from the fill, the number of
    groups filled, and the largest group's area and the centroid of its
    largest triangle (0, nan, nan with nothing filled)."""

    codes: npt.NDArray[np.int32]
    filled: npt.NDArray[np.bool_]
    groups: int
    largest: tuple[float, float, float]


def fill_no_class(
    vertices: npt.ArrayLike,
    triangles: npt.ArrayLike,
    codes: npt.ArrayLike,
    *,
    water: frozenset[int],
    asking: npt.ArrayLike | None = None,
) -> Filled:
    """Each group of `asking` triangles (default: code 0) connected through
    shared edges takes the code of the nearest donor triangle: one not asking
    whose code is neither 0 nor in `water`, else (none such) any not asking
    with a code. The distance is between closed triangles in the vertices'
    plane, exact; a tie goes to the code whose donors share the most edge
    length with the group, then to the lower code. No donor at all is a
    `ValueError`. The result depends on no triangle or vertex order."""
    xy = np.asarray(vertices, dtype=np.float64)[:, :2]
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    out = np.array(codes, dtype=np.int32).reshape(-1)
    ask = out == 0 if asking is None else np.asarray(asking, dtype=bool).reshape(-1)
    if not ask.any():
        return Filled(out, ask & False, 0, (0.0, float("nan"), float("nan")))
    coded = ~ask & (out != 0)
    donor = coded & ~np.isin(out, list(water))
    donor = donor if donor.any() else coded
    if not donor.any():
        raise ValueError(
            "the land cover gives no triangle of the mesh a class, so code 0 cannot be "
            "filled; it does not cover the domain"
        )
    t = len(tri)
    sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    keys = _keys(sides, len(xy))
    order = np.argsort(keys, kind="stable")
    keys, owner = keys[order], order % t
    pair = np.flatnonzero(keys[1:] == keys[:-1])
    u, v = owner[pair], owner[pair + 1]
    once = np.ones(len(keys), dtype=bool)
    once[pair] = once[pair + 1] = False
    border = np.zeros(t, dtype=bool)
    border[owner[once]] = True
    where = np.flatnonzero(ask)
    group = np.full(t, -1, dtype=np.int64)
    group[where] = where[regions(tri[where], np.zeros((0, 2), np.int64))]
    # A group's nearest donor is nearest to one of its rim triangles, and lies
    # on the donors' rim: a vertex shared with a non-donor or on the border.
    rim = ask & border
    np.logical_or.at(rim, u, ask[u] & ~ask[v])
    np.logical_or.at(rim, v, ask[v] & ~ask[u])
    open_vertex = np.zeros(len(xy), dtype=bool)
    open_vertex[tri[~donor].ravel()] = True
    open_vertex[sides[order[once]].ravel()] = True
    edge = np.flatnonzero(donor & open_vertex[tri].any(axis=1))
    queries = np.flatnonzero(rim)
    asked, given = shapely.polygons(xy[tri[queries]]), shapely.polygons(xy[tri[edge]])
    q, found = shapely.STRtree(given).query_nearest(asked, all_matches=True)
    gap = shapely.distance(asked[q], given[found])
    g = group[queries[q]]
    best = np.full(t, np.inf)
    np.minimum.at(best, g, gap)
    near = gap == best[g]
    candidates: dict[int, set[int]] = defaultdict(set)
    for gi, c in zip(g[near], out[edge[found[near]]], strict=True):
        candidates[int(gi)].add(int(c))
    shared: dict[tuple[int, int], float] = defaultdict(float)
    across = (ask[u] & donor[v]) | (ask[v] & donor[u])
    a, b = np.where(ask[u], u, v)[across], np.where(ask[u], v, u)[across]
    ends = sides[order[pair[across]]]
    lengths = np.hypot(*(xy[ends[:, 1]] - xy[ends[:, 0]]).T)
    for gi, c, length in zip(group[a], out[b], lengths, strict=True):
        shared[int(gi), int(c)] += float(length)
    picks = np.zeros(t, dtype=np.int32)
    for gi, options in candidates.items():
        picks[gi] = min(options, key=lambda c: (-round(shared[gi, c], DECIMALS), c))
    out[where] = picks[group[where]]
    area = plan_areas(xy, tri)
    totals = np.zeros(t)
    np.add.at(totals, group[where], area[where])
    top = int(min(candidates, key=lambda gi: (-totals[gi], gi)))
    members = np.flatnonzero(group == top)
    biggest = members[np.argmax(area[members])]
    cx, cy = xy[tri[biggest]].mean(axis=0)
    return Filled(
        out, ask.copy(), len(np.unique(group[where])), (float(totals[top]), float(cx), float(cy))
    )


def plan_areas(xy: npt.NDArray[np.float64], tri: npt.NDArray[np.int64]) -> npt.NDArray[np.float64]:
    """Each triangle's area in the plane of ``xy``."""
    a, b, c = (xy[tri[:, k]] for k in range(3))
    cross = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    return np.asarray(0.5 * np.abs(cross), dtype=np.float64)
