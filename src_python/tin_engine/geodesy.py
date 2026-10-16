"""Metres near a point, cell sizes and geodesic areas (increment 40).

`docs/increments/40-geographic-catchment.md`, sections 5.2 to 5.4. A DEM in
degrees is flooded on its own grid; what is a distance is measured in a
:class:`LocalFrame`, an affine map to metres at one point, and what is an
area is geodesic, on the ellipsoid of the CRS. On a projected CRS every
function returns what the metre code computed before, bit for bit.

Layer 1: pyproj through `crs.parse_crs`, NumPy and shapely. No `_core`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from shapely.geometry.base import BaseGeometry

from tin_engine.crs import parse_crs
from tin_engine.io.models import RasterMeta


@dataclass(frozen=True, slots=True)
class LocalFrame:
    """An affine map from a CRS's coordinates to metres near one point.
    Identity for a projected CRS. For a CRS in degrees:
    X = (lon - lon0) * kx,  Y = (lat - lat0) * ky, with kx = (pi/180) N cos(lat0)
    and ky = (pi/180) M, N and M the ellipsoid's prime-vertical and meridian
    radii at lat0: the exact metres per degree at lat0."""

    lon0: float
    lat0: float
    kx: float
    ky: float

    @classmethod
    def identity(cls) -> LocalFrame:
        return cls(0.0, 0.0, 1.0, 1.0)

    @classmethod
    def at(cls, crs: str, point: tuple[float, float]) -> LocalFrame:
        """The frame at `point` (in `crs`): the identity unless `crs` is in
        degrees, then the metres per degree of the CRS's own ellipsoid."""
        parsed = parse_crs(crs)
        if not parsed.is_geographic:
            return cls.identity()
        unit = parsed.axis_info[0].unit_name
        if unit != "degree":
            raise ValueError(f"a geographic CRS in {unit}, not in degrees, is not supported")
        geod = parsed.get_geod()
        assert geod is not None  # a geographic CRS has an ellipsoid
        lon0, lat0 = float(point[0]), float(point[1])
        s = math.sin(math.radians(lat0))
        w = 1.0 - geod.es * s * s
        n, m = geod.a / math.sqrt(w), geod.a * (1.0 - geod.es) / w**1.5
        kx = math.pi / 180.0 * n * math.cos(math.radians(lat0))
        return cls(lon0, lat0, kx, math.pi / 180.0 * m)

    @property
    def is_identity(self) -> bool:
        return (self.lon0, self.lat0, self.kx, self.ky) == (0.0, 0.0, 1.0, 1.0)

    def forward(self, xy: npt.ArrayLike) -> npt.NDArray[np.float64]:
        """`(N, 2)` points of the CRS in metres of the frame."""
        p = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        if self.is_identity:
            return p
        return np.column_stack([(p[:, 0] - self.lon0) * self.kx, (p[:, 1] - self.lat0) * self.ky])

    def inverse(self, xy: npt.ArrayLike) -> npt.NDArray[np.float64]:
        """`(N, 2)` points in metres of the frame back in the CRS."""
        p = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        if self.is_identity:
            return p
        return np.column_stack([p[:, 0] / self.kx + self.lon0, p[:, 1] / self.ky + self.lat0])

    def units(self, metres: float) -> tuple[float, float]:
        """`metres` in the CRS's units, east-west and north-south."""
        return metres / self.kx, metres / self.ky

    def cell_m(self, meta: RasterMeta) -> float:
        """A cell of `meta`: its longer side in metres (section 5.3)."""
        return max(meta.delta_x * self.kx, meta.delta_y * self.ky)

    def meta(self, meta: RasterMeta) -> RasterMeta:
        """`meta`'s node grid in the frame, for the burn; `meta` itself on
        the identity. The `crs` is a label, never parsed."""
        if self.is_identity:
            return meta
        update = {
            "x_min": (meta.x_min - self.lon0) * self.kx,
            "y_max": (meta.y_max - self.lat0) * self.ky,
            "delta_x": meta.delta_x * self.kx,
            "delta_y": meta.delta_y * self.ky,
            "geographic": False,
            "epsg": None,
            "crs": f"metres from ({self.lon0}, {self.lat0}) in {meta.crs}",
        }
        return meta.model_copy(update=update)


def cell_areas_by_row(meta: RasterMeta) -> npt.NDArray[np.float64]:
    """m² of one cell in each row: geodesic (the CRS's ellipsoid, one `Geod`
    polygon per row) when the grid is in degrees; dx * dy in every row
    otherwise. A row's cell is centred on its nodes: latitude y_r - dy/2 to
    y_r + dy/2, dx wide."""
    if not meta.geographic:
        return np.full(meta.rows, meta.delta_x * meta.delta_y, dtype=np.float64)
    geod = parse_crs(meta.crs).get_geod()
    assert geod is not None
    lons, half = [0.0, meta.delta_x, meta.delta_x, 0.0], meta.delta_y / 2
    out = np.empty(meta.rows, dtype=np.float64)
    for r in range(meta.rows):
        y = meta.y_max - r * meta.delta_y
        out[r] = abs(geod.polygon_area_perimeter(lons, [y + half, y + half, y - half, y - half])[0])
    return out


def mean_cell_area(mask: npt.NDArray[np.uint8], meta: RasterMeta) -> float:
    """m² per in-node of `mask`: dx * dy exactly on a projected grid; on a
    grid in degrees, the row-weighted sum over the in-nodes over their count."""
    if not meta.geographic:
        return meta.delta_x * meta.delta_y
    counts = np.count_nonzero(mask, axis=1)
    return float(counts @ cell_areas_by_row(meta) / counts.sum())


def area_m2(geometry: BaseGeometry, crs: str) -> float:
    """`geometry.area` on a projected CRS; the absolute geodesic area on the
    CRS's ellipsoid on one in degrees."""
    parsed = parse_crs(crs)
    if not parsed.is_geographic:
        return float(geometry.area)
    geod = parsed.get_geod()
    assert geod is not None
    return abs(geod.geometry_area_perimeter(geometry)[0])


__all__ = ["LocalFrame", "area_m2", "cell_areas_by_row", "mean_cell_area"]
