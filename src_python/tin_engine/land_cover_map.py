"""The one land-cover map: the class-polygon files joined, cut by the preferred sources.

`docs/increments/111-one-land-cover-map.md`, *3. Design*.

The files' polygons, in the order given, are the base: where two files share
area, the first given keeps it, and a later file's polygon is cut by the
earlier files' polygons it shares area with. Each preferred source then
replaces its class (`preferred.replace`, near the read region), and each
refilled piece is unioned with the polygons of the class it took whose
borders meet its own, so no old border runs between them. Each polygon
reaching the box of the read region widened by `MARGIN` is valid, as the
readers refuse any other (a ring that crosses itself makes GEOS's overlays
throw), and is cut to that box. So every polygon of the map's ``areas`` is
valid, and may be overlaid unchecked. A polygon with no area in the box
goes to ``outside`` as read, whole: nothing overlays it,
and only the edges step reads it, by its rings. Pure: shapely geometry in,
no file, no CRS.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry

from tin_engine.feature_input import GROUPS, MARGIN, LandCoverClasses, _polygonal
from tin_engine.preferred import BAND, Replaced, replace

Pairs = Sequence[tuple[int, BaseGeometry]]


@dataclass(frozen=True)
class LandCoverMap:
    """The map's ``(code, polygon)`` areas in the DEM's CRS, each valid; the
    polygons the read box misses, as read; the class file; the ``--prefer``
    cut (None without ``--prefer``); per overlap between two files, the
    index of the file kept and the area cut from the later one."""

    areas: tuple[tuple[int, BaseGeometry], ...]
    outside: tuple[tuple[int, BaseGeometry], ...]
    classes: LandCoverClasses
    replaced: Replaced | None
    overlap: tuple[tuple[int, BaseGeometry], ...]


def build(
    base: Sequence[Pairs],
    preferred: Sequence[tuple[int, Sequence[BaseGeometry]]],
    classes: LandCoverClasses,
    region: Polygon,
    *,
    cut: bool = True,
    outside: Pairs = (),
    threads: int = 1,
) -> LandCoverMap:
    """The map of ``base`` (one sequence of (code, polygon) per file, in the
    order given) with each code of ``preferred`` (replaced code, polygons, as
    `preferred.replace` takes them) taken from its polygons near ``region``.
    The polygons reaching ``region`` are taken as valid, as the readers
    check them. With ``cut``, each file's polygons are first cut near
    ``region``, those with no area there set aside (`_fixed`). Without it
    ``base`` is taken as already so, as the ``areas`` of a map built with
    ``cut`` are, and ``outside`` (that map's) is passed through. The cut
    runs in ``threads`` threads."""
    xy = shapely.get_coordinates(region.exterior)
    fixed = [_fixed(pairs, xy, threads) if cut else (pairs, ()) for pairs in base]
    areas, overlap = _joined([f for f, _ in fixed])
    far = (*outside, *(pair for _, out in fixed for pair in out))
    if not preferred:
        return LandCoverMap(tuple(areas), far, classes, None, tuple(overlap))
    done = replace(areas, preferred, [region])
    return LandCoverMap(tuple(_refilled_joined(done)), far, classes, done, tuple(overlap))


def _fixed(pairs: Pairs, xy: npt.NDArray[np.float64], threads: int) -> tuple[Pairs, Pairs]:
    """``pairs`` reaching the box of ``xy`` widened by `MARGIN`, each valid,
    cut to the box with `shapely.intersection`, kept when its polygonal part
    has area; and the rest as given. One whose bounding box misses the box,
    or lies inside it, is kept as read uncut. Never `shapely.clip_by_rect`:
    it is wrong on a valid polygon with an edge through a corner of the box."""
    lo, hi = xy.min(axis=0) - MARGIN, xy.max(axis=0) + MARGIN
    box = shapely.box(*lo, *hi)
    shapes = np.array([g for _, g in pairs], dtype=object)
    x0, y0, x1, y1 = shapely.bounds(shapes).reshape(-1, 4).T
    near = (x0 <= hi[0]) & (x1 >= lo[0]) & (y0 <= hi[1]) & (y1 >= lo[1])
    within = (x0 >= lo[0]) & (x1 <= hi[0]) & (y0 >= lo[1]) & (y1 <= hi[1])
    cut = shapes.copy()
    ends = np.cumsum(shapely.get_num_coordinates(shapes))
    edges = np.searchsorted(ends, np.linspace(0, ends[-1] if len(ends) else 0, GROUPS + 1)[1:-1])
    groups = [ix[near[ix]] for ix in np.split(np.arange(len(pairs)), edges)]

    def cut_group(ix: npt.NDArray[np.intp]) -> None:
        for i in ix:
            cut[i] = _cut(shapes[i], box if not within[i] else None)

    with ThreadPoolExecutor(max(threads, 1)) as pool:
        list(pool.map(cut_group, [ix for ix in groups if len(ix)]))
    areas: list[tuple[int, BaseGeometry]] = []
    outside: list[tuple[int, BaseGeometry]] = []
    for (code, g), c, reach in zip(pairs, cut, near, strict=True):
        part = c if reach else None
        (outside if part is None else areas).append((code, g if part is None else part))
    return areas, outside


def _cut(geometry: BaseGeometry, box: Polygon | None) -> BaseGeometry | None:
    """The valid ``geometry`` itself if a Polygon inside ``box`` (None), else
    its polygonal part cut to ``box`` part by part, so its parts keep their
    order; one part as a Polygon; None if it has no area."""
    if box is None and isinstance(geometry, Polygon):
        return geometry
    parts = shapely.get_parts(geometry)
    pieces = shapely.get_parts(shapely.get_parts(parts if box is None else parts & box))
    part = _polygonal(shapely.geometrycollections(list(pieces)))
    return part if part is not None and part.area > 0 else None


def _joined(
    base: Sequence[Pairs],
) -> tuple[list[tuple[int, BaseGeometry]], list[tuple[int, BaseGeometry]]]:
    """The files' polygons, a later one less what it shares with earlier
    ones (pairs that only touch share none); and each area cut, with the
    index of the file that kept it."""
    kept: list[tuple[int, BaseGeometry]] = list(base[0]) if base else []
    owner = [0] * len(kept)
    overlap: list[tuple[int, BaseGeometry]] = []
    for j, pairs in enumerate(base[1:], start=1):
        hits: dict[int, list[int]] = defaultdict(list)
        if kept and pairs:
            shapes = np.array([g for _, g in kept], dtype=object)
            later = np.array([g for _, g in pairs], dtype=object)
            found = shapely.STRtree(shapes).query(later, predicate="intersects")
            for i, k in zip(*found, strict=True):
                shared = shapely.intersection(later[i], shapes[k])
                if shared.area > 0:
                    hits[int(i)].append(int(k))
                    overlap.append((owner[k], shared))
        added = []
        for i, (code, polygon) in enumerate(pairs):
            if i in hits:
                polygon = later[i].difference(shapely.union_all([shapes[k] for k in hits[i]]))
                parts = [
                    p
                    for p in shapely.get_parts(polygon)
                    if p.geom_type == "Polygon" and not p.is_empty
                ]
                if not parts:
                    continue
                polygon = parts[0] if len(parts) == 1 else MultiPolygon(parts)
            added.append((code, polygon))
        kept += added
        owner += [j] * len(added)
    return kept, overlap


def _refilled_joined(done: Replaced) -> list[tuple[int, BaseGeometry]]:
    """The kept polygons, each joined with the refilled pieces and kept
    polygons of its class whose borders meet its own (within the neighbour
    search's band, through each other too), the group under its first kept
    polygon; then the preferred polygons."""
    nodes = [(c, g) for _, c, g in done.kept] + [(n, p) for _, p, n in done.refilled]
    parent = list(range(len(nodes)))

    def root(i: int) -> int:
        while parent[i] != i:
            i = parent[i]
        return i

    k = len(done.kept)
    if done.refilled:
        pieces = np.array([p for _, p, _ in done.refilled], dtype=object)
        bands = shapely.buffer(shapely.boundary(pieces), BAND)
        rings = shapely.boundary(np.array([g for _, g in nodes], dtype=object))
        for b, j in zip(*shapely.STRtree(rings).query(bands, predicate="intersects"), strict=True):
            if nodes[j][0] == nodes[k + b][0]:
                a, c = root(k + int(b)), root(int(j))
                parent[max(a, c)] = min(a, c)  # the root is the group's first
    groups: dict[int, list[BaseGeometry]] = {}
    for i, (_, g) in enumerate(nodes):
        groups.setdefault(root(i), []).append(g)
    joined = [
        (nodes[r][0], gs[0] if len(gs) == 1 else shapely.union_all(gs)) for r, gs in groups.items()
    ]
    return joined + list(done.given)
