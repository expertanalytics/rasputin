"""Chosen classes taken from preferred polygons, the rest from the base land cover.

`docs/increments/107-prefer-classes.md`, *2. What "prefer" means where the sources disagree*.

Each preferred polygon that meets a triangle takes its class; it is cut out of
every base polygon, and an earlier source's polygons are cut out of a later
one's. The base polygons of a replaced class are dropped: each connected piece
left of them takes the class of the kept base polygons sharing the longest
border with it (a tie: the lower code), and a piece with no such neighbour is
left out, uncovered. Pure: shapely geometry and arrays in, no file, no CRS.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry.base import BaseGeometry

from tin_engine.cells import Float64
from tin_engine.polygon_fractions import subdivide

#: Half the width of the band around a piece's boundary in which a kept
#: polygon's boundary counts as shared, in the CRS's unit: subdivision puts
#: points on the cut lines, so a shared border is only nearly equal on both
#: sides (a chosen width, far above coordinate rounding in metres and far
#: below a land cover's narrowest polygon).
BAND = 0.01
#: Border lengths are compared rounded to this many decimals, so a tie stays a tie.
DECIMALS = 6

Pairs = Sequence[tuple[int, BaseGeometry]]


@dataclass(frozen=True)
class Areas:
    """Inside the triangles: the replaced class's area in the base, the area
    the preferred polygons give it, and the area its leftover pieces gave to
    neighbouring classes."""

    replaced: float
    given: float
    refilled: float


def _array(geometries: Sequence[BaseGeometry]) -> npt.NDArray[Any]:
    array = np.empty(len(geometries), dtype=object)
    array[:] = list(geometries)
    return array


def _cut(polygons: Pairs, by: Sequence[BaseGeometry]) -> list[tuple[int, BaseGeometry]]:
    """Each polygon less the union of those in `by` it meets, empty ones dropped."""
    shapes = _array([g for _, g in polygons])
    cutters = _array(by)
    hits: dict[int, list[int]] = defaultdict(list)
    if len(cutters) and len(shapes):
        for i, j in zip(
            *shapely.STRtree(cutters).query(shapes, predicate="intersects"), strict=True
        ):
            hits[int(i)].append(int(j))
    cut = []
    for i, (code, polygon) in enumerate(polygons):
        if i in hits:
            polygon = polygon.difference(shapely.union_all(cutters[hits[i]]))
        if not polygon.is_empty:
            cut.append((code, polygon))
    return cut


def _neighbours(pieces: Sequence[BaseGeometry], kept: Pairs) -> list[int | None]:
    """Per piece, the class of the kept polygons sharing the longest border with it."""
    bands = shapely.buffer(shapely.boundary(_array(pieces)), BAND)
    shapes = _array([g for _, g in kept])
    near = np.unique(shapely.STRtree(bands).query(shapes, predicate="intersects")[0])
    subs: list[BaseGeometry] = []
    codes: list[int] = []
    for k in near:
        cut = subdivide(shapes[k])
        subs += cut
        codes += [kept[k][0]] * len(cut)
    lengths: list[dict[int, float]] = [defaultdict(float) for _ in pieces]
    if subs:
        sub_array = _array(subs)
        b, s = shapely.STRtree(sub_array).query(bands, predicate="intersects")
        shared = shapely.length(shapely.intersection(bands[b], shapely.boundary(sub_array[s])))
        for i, j, length in zip(b, s, shared, strict=True):
            lengths[i][codes[j]] += float(length)
    return [
        max(per, key=lambda c: (round(per[c], DECIMALS), -c)) if any(per.values()) else None
        for per in lengths
    ]


@dataclass(frozen=True)
class Replaced:
    """The cut and refill: per kept base polygon its index in the base, its
    code and what is left of it; the preferred polygons as (code, geometry);
    per refilled piece the code it replaced, the piece and the class it took."""

    kept: list[tuple[int, int, BaseGeometry]]
    given: list[tuple[int, BaseGeometry]]
    refilled: list[tuple[int, BaseGeometry, int]]

    @property
    def polygons(self) -> list[tuple[int, BaseGeometry]]:
        """Every polygon as (code, geometry): kept, given, refilled."""
        kept = [(c, g) for _, c, g in self.kept]
        return [*kept, *self.given, *((n, p) for _, p, n in self.refilled)]


def replace(
    base: Pairs,
    chosen: Sequence[tuple[int, Sequence[BaseGeometry]]],
    near: Sequence[BaseGeometry] | npt.NDArray[Any],
) -> Replaced:
    """`base` (code, polygon) with each code of `chosen` (replaced code,
    preferred polygons, first wins an overlap, the same code twice unions)
    taken from its preferred polygons; a preferred polygon or a leftover piece
    is kept only where it meets a geometry of `near`."""
    tree = shapely.STRtree(near)

    def meeting(pairs: Pairs) -> list[tuple[int, BaseGeometry]]:
        if not pairs:
            return []
        found = tree.query(_array([g for _, g in pairs]), predicate="intersects")[0]
        return [pairs[i] for i in np.unique(found)]

    preferred: list[tuple[int, BaseGeometry]] = []
    for code, polygons in chosen:
        own = meeting([(code, g) for g in polygons])
        preferred += _cut(own, [g for _, g in preferred])
    replaced = {code for code, _ in chosen}
    cut = _cut([(i, g) for i, (_, g) in enumerate(base)], [g for _, g in preferred])
    kept = [(i, base[i][0], g) for i, g in cut if base[i][0] not in replaced]
    left = meeting(
        [(base[i][0], p) for i, g in cut if base[i][0] in replaced for p in shapely.get_parts(g)]
    )
    classes = _neighbours([p for _, p in left], [(c, g) for _, c, g in kept]) if left else []
    refilled = [(c, p, n) for (c, p), n in zip(left, classes, strict=True) if n is not None]
    return Replaced(kept, preferred, refilled)


def areas_inside(
    base: Pairs, done: Replaced, codes: Sequence[int], corners: Float64
) -> dict[int, Areas]:
    """Per replaced code of ``codes``, its areas inside the triangles
    ``corners`` (T, 3, k >= 2): the base polygons' of that code (``base``,
    valid), the preferred polygons' and the refilled pieces' of `done`."""
    xy = np.asarray(corners, np.float64)[:, :, :2]
    triangles = shapely.polygons(xy) if len(xy) else np.zeros(0, dtype=object)
    tree = shapely.STRtree(triangles)

    def inside(polygons: Sequence[BaseGeometry]) -> float:
        # Each triangle cut by the few-vertex pieces it meets, not by one union
        # of thousands of polygons; the polygons given never overlap.
        pieces = _array([s for g in polygons for s in subdivide(g)])
        if not len(pieces):
            return 0.0
        p, t = tree.query(pieces, predicate="intersects")
        return float(shapely.area(shapely.intersection(triangles[t], pieces[p])).sum())

    return {
        code: Areas(
            replaced=inside([g for c, g in base if c == code]),
            given=inside([g for c, g in done.given if c == code]),
            refilled=inside([p for c, p, _ in done.refilled if c == code]),
        )
        for code in codes
    }
