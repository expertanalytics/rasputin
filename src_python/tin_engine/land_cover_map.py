"""The one land-cover map: the class-polygon files joined, cut by the preferred sources.

`docs/increments/111-one-land-cover-map.md`, *3. Design*.

The files' polygons, in the order given, are the base: where two files share
area, the first given keeps it, and a later file's polygon is cut by the
earlier files' polygons it shares area with. Each preferred source then
replaces its class (`preferred.replace`, near the read region), and each
refilled piece is unioned with the polygons of the class it took whose
borders meet its own, so no old border runs between them. The polygons stay
unchecked: a ring that crosses itself is repaired only where the shares
need it. Pure: shapely geometry in, no file, no CRS.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry

from tin_engine.feature_input import LandCoverClasses, _repaired
from tin_engine.preferred import BAND, Replaced, replace

Pairs = Sequence[tuple[int, BaseGeometry]]


@dataclass(frozen=True)
class LandCoverMap:
    """The map's ``(code, polygon)`` areas in the DEM's CRS, unchecked; the
    class file; the ``--prefer`` cut (None without ``--prefer``); and per
    overlap between two files, the index of the file kept and the area cut
    from the later one."""

    areas: tuple[tuple[int, BaseGeometry], ...]
    classes: LandCoverClasses
    replaced: Replaced | None
    overlap: tuple[tuple[int, BaseGeometry], ...]


def build(
    base: Sequence[Pairs],
    preferred: Sequence[tuple[int, Sequence[BaseGeometry]]],
    classes: LandCoverClasses,
    region: Polygon,
) -> LandCoverMap:
    """The map of ``base`` (one sequence of (code, polygon) per file, in the
    order given) with each code of ``preferred`` (replaced code, polygons, as
    `preferred.replace` takes them) taken from its polygons near ``region``."""
    areas, overlap = _joined(base, shapely.get_coordinates(region.exterior))
    if not preferred:
        return LandCoverMap(tuple(areas), classes, None, tuple(overlap))
    done = replace(areas, preferred, [region])
    return LandCoverMap(tuple(_refilled_joined(done)), classes, done, tuple(overlap))


def _joined(
    base: Sequence[Pairs], xy: npt.NDArray[np.float64]
) -> tuple[list[tuple[int, BaseGeometry]], list[tuple[int, BaseGeometry]]]:
    """The files' polygons, a later one less what it shares with earlier
    ones (pairs that only touch share none); and each area cut, with the
    index of the file that kept it. The polygons of each pair are repaired
    near the points ``xy`` first (`_repaired`), as the shares step repairs."""
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
                for n, g in ((i, later), (k, shapes)):
                    g[n] = _repaired(g[n], xy) or Polygon()
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
        (nodes[r][0], gs[0] if len(gs) == 1 else shapely.union_all(gs))
        for r, gs in groups.items()
        if r < k  # a piece always has a kept neighbour; else it is left uncovered
    ]
    return joined + list(done.given)
