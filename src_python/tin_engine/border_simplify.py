"""Land-cover borders simplified within a band, each polygon keeping its area
(increment 32, ``docs/increments/32-landcover-simplify.md``, section 6).

Polygons to flat rings for :func:`tin_engine._core.simplify_borders` and back,
with the same polygons, parts and holes in the same order. No I/O.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry

from tin_engine import _core
from tin_engine._core import BorderCounts

__all__ = ["BorderResult", "simplify_borders"]

_REFUSED = {
    _core.BorderStatus.InvalidBand: "the band must be finite and at least 0, got {band}",
    _core.BorderStatus.BadRings: "a ring has fewer than three vertices or a coordinate that is "
    "not finite",
    _core.BorderStatus.InvalidClearance: "the clearance must be finite and at least 0, got "
    "{clearance}",
}


@dataclass(frozen=True, slots=True)
class BorderResult:
    """The simplified polygons, same count and order as the input, and the counts."""

    polygons: tuple[BaseGeometry, ...]
    counts: BorderCounts


def simplify_borders(
    polygons: Sequence[BaseGeometry], band_m: float, clearance_m: float = 0.0
) -> BorderResult:
    """Every border of the coverage ``polygons`` moved at most ``band_m`` from
    its source, every part keeping its area, junctions and the outer boundary
    fixed, no new vertex or edge closer than ``clearance_m`` to a vertex or edge
    it does not share an end with. A refused input is a ``ValueError`` saying why."""
    parts = [[p for p in shapely.get_parts(g) if not p.is_empty] for g in polygons]
    rings = [
        shapely.get_coordinates(r)[:-1] for ps in parts for p in ps for r in shapely.get_rings(p)
    ]
    points = np.concatenate(rings) if rings else np.empty((0, 2))
    first = np.cumsum([0, *map(len, rings)], dtype=np.uint64)
    out = _core.simplify_borders(points, first, band_m, clearance=clearance_m)
    if out.status != _core.BorderStatus.Ok:
        why = _REFUSED[out.status].format(band=band_m, clearance=clearance_m)
        raise ValueError(f"simplify_borders: {why}")
    xy, starts = np.asarray(out.points), np.asarray(out.ring_starts).tolist()
    pieces = iter(xy[starts[k] : starts[k + 1]] for k in range(len(rings)))
    simplified: list[BaseGeometry] = []
    for g, ps in zip(polygons, parts, strict=True):
        new = [Polygon(next(pieces), [next(pieces) for _ in p.interiors]) for p in ps]
        simplified.append(g if not new else new[0] if isinstance(g, Polygon) else MultiPolygon(new))
    return BorderResult(tuple(simplified), out.counts)
