"""Sub-catchment units as mesh constraints, and each cell's unit id.

The units are selected by their inner point, repaired into an exact coverage,
their borders simplified within a band that the triangle budget sets, and
their rings given to the mesh as one constraint feature. No I/O.

The band's derivation: `docs/increments/106-wednesday-cells.md`, *What it does*.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt
import shapely
from pyproj import CRS
from pyproj.crs import ProjectedCRS
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
from shapely.geometry import LineString
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient

from tin_engine.border_simplify import simplify_borders
from tin_engine.crs import parse_crs
from tin_engine.domain import DomainPolygon
from tin_engine.features import TerrainFeature

#: Metres: the snapping and gap width of the repair, and the clearance of the
#: simplification, as for land cover.
REPAIR_M = 1.0
#: Metres added to twice the band around the units' box for the frame, so the
#: frame always reaches beyond the band.
FRAME_M = 1000.0
#: The band is this fraction of a triangle's edge at the budget: the
#: simplification leaves segments about eight bands long on natural borders.
EDGES_PER_BAND = 8.0
TOP_ID = 2**31 - 1


def select(units: Sequence[BaseGeometry], domain: BaseGeometry) -> list[int]:
    """The indices, in order, of the units whose inner point lies in `domain`."""
    inner = shapely.point_on_surface(np.asarray(units, dtype=object))
    return [int(k) for k in np.flatnonzero(shapely.covers(domain, inner))]


def repair(units: Sequence[BaseGeometry]) -> list[BaseGeometry]:
    """The units as an exact coverage: gaps and overlaps under 1 m closed."""
    cleaned = shapely.coverage_clean(
        np.asarray(units, dtype=object), gap_width=REPAIR_M, snapping_distance=REPAIR_M
    )
    return list(cleaned)


def simplify(units: Sequence[BaseGeometry], band_m: float) -> list[BaseGeometry]:
    """The units, in order, each border moved at most `band_m` and each unit
    keeping its area. A frame around their union is simplified with them, so
    the union's outline is a border like the others and no gap opens."""
    union = shapely.coverage_union_all(units)
    grow = 2.0 * band_m + FRAME_M
    x0, y0, x1, y1 = union.bounds
    frame = shapely.box(x0 - grow, y0 - grow, x1 + grow, y1 + grow).difference(union)
    out = simplify_borders([*units, frame], band_m, REPAIR_M).polygons
    return list(out[: len(units)])


def band_m(area_m2: float, triangles: int) -> float:
    """An eighth of the edge of an equilateral triangle of area
    `area_m2 / triangles`: a border then gets about one vertex per triangle edge."""
    return math.sqrt(4.0 * area_m2 / (math.sqrt(3.0) * triangles)) / EDGES_PER_BAND


def equal_area(domain: DomainPolygon) -> DomainPolygon:
    """`domain` itself when projected; in degrees, moved to the Lambert
    azimuthal equal-area projection centred on its centroid, on its CRS's
    ellipsoid, so areas are true and distances near it nearly so."""
    crs = parse_crs(domain.crs)
    if not crs.is_geographic:
        return domain
    centre = domain.polygon.centroid
    laea = LambertAzimuthalEqualAreaConversion(centre.y, centre.x)
    return domain.to_crs(CRS(ProjectedCRS(conversion=laea, geodetic_crs=crs.geodetic_crs)))


def simplify_domain(
    domain: DomainPolygon, triangles: int | None, tolerance: float | None
) -> tuple[DomainPolygon, float]:
    """`domain`'s outline and holes, each ring moved at most the band and the
    area kept, in its own CRS, and the band: `tolerance` when given, else
    `band_m` of its area at `triangles`, both in metres of `equal_area`."""
    work = equal_area(domain)
    if tolerance is not None:
        band = tolerance
    else:
        assert triangles is not None
        band = band_m(work.polygon.area, triangles)
    if band > 0:
        work = work.model_copy(update={"polygon": orient(simplify([work.polygon], band)[0])})
    return work.to_crs(domain.crs), band


def check_ids(fids: Sequence[object], name: str) -> list[int]:
    """The fids as ids Shyft's `catchment_id` holds; any other is a ValueError."""
    for fid in fids:
        if isinstance(fid, bool) or not isinstance(fid, int) or not 1 <= fid <= TOP_ID:
            raise ValueError(f"{name}: feature {fid} has no id Shyft's catchment_id can hold")
    return [int(f) for f in fids]  # type: ignore[call-overload]


def cell_ids(
    units: Sequence[tuple[int, BaseGeometry]], centroids: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.int64], int]:
    """Each centroid's unit id, the nearest unit's when it lies in none, and
    how many lay in none."""
    ids = np.array([i for i, _ in units], dtype=np.int64)
    tree = shapely.STRtree([g for _, g in units])
    points = shapely.points(centroids)
    hit, unit = tree.query(points, predicate="within")
    out = np.full(len(points), -1, dtype=np.int64)
    out[hit] = unit
    strays = np.flatnonzero(out < 0)
    if len(strays):
        out[strays] = tree.query_nearest(points[strays], all_matches=False)[1]
    return ids[out], len(strays)


def borders(units: Sequence[BaseGeometry], mask: int) -> TerrainFeature:
    """Every ring of the units as one constraint feature carrying `mask`."""
    rings = shapely.get_rings(shapely.get_parts(np.asarray(units, dtype=object)))
    lines = tuple(LineString(shapely.get_coordinates(r)) for r in rings)
    return TerrainFeature(fid="catchments", mask=mask, lines=lines)
