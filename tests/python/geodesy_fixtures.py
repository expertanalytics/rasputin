"""Ellipsoid oracles and grids in degrees for increment 40's suites.

`docs/increments/40-geographic-catchment.md`, sections 5.2 to 5.4 and 7. Test
support only, sharing no code with `tin_engine.geodesy`: the metres per degree
and the cell areas below are the closed forms written out here, from the
ellipsoid's `a` and `f`, so a test that compares the module with them compares
two independent computations.

- `metres_per_degree(lat)`: `kx = (pi/180) N cos(lat)`, `ky = (pi/180) M`, N
  and M the prime-vertical and meridian radii at `lat` (section 5.2).
- `band_area(lat0, lat1, width)`: the area of the band between two latitudes,
  `width` degrees of longitude wide, by the authalic closed form (R4), not by
  pyproj's `Geod`.

THE TWIN of `gauge_fixtures.valley()` (R7, R8): the same array on a grid in
EPSG:4326 whose spacing is 10 m at its centre latitude (`dx = 10 / kx`,
`dy = 10 / ky`), the centre node (row `gf.RIM_ROW`, column `gf.CC`) put where
the projected valley's centre node lies in EPSG:25833. A lattice position maps
to the twin by `twin_xy`; a point of the projected valley maps to within about
half a metre of the same lattice position on the twin (UTM's scale, 0.9996 at
the central meridian, over 120 rows; the meridian convergence 300 m east of it),
a twentieth of a cell.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import numpy.typing as npt
from pyproj import Transformer

import gauge_fixtures as gf
from tin_engine.io.models import DemTile, RasterMeta

#: WGS 84, as EPSG:4326 names it.
WGS84_A, WGS84_F = 6378137.0, 1.0 / 298.257223563
#: International 1924 (ED50, EPSG:4230): an ellipsoid far enough from WGS 84
#: (a 251 m longer semi-major axis) to show which one a module read.
INTL_A, INTL_F = 6378388.0, 1.0 / 297.0
#: One arc-second, GLO-30's spacing.
ARCSEC = 1.0 / 3600.0


def metres_per_degree(lat: float, a: float = WGS84_A, f: float = WGS84_F) -> tuple[float, float]:
    """`(kx, ky)`: metres per degree of longitude and of latitude at `lat`."""
    e2 = f * (2.0 - f)
    s = math.sin(math.radians(lat))
    w = 1.0 - e2 * s * s
    n, m = a / math.sqrt(w), a * (1.0 - e2) / w**1.5
    return math.pi / 180.0 * n * math.cos(math.radians(lat)), math.pi / 180.0 * m


def _authalic(lat: float, a: float, f: float) -> float:
    """The ellipsoid's area from the equator to `lat`, per radian of longitude."""
    e2 = f * (2.0 - f)
    e, b = math.sqrt(e2), a * (1.0 - f)
    s = math.sin(math.radians(lat))
    return (
        b * b / 2.0 * (s / (1.0 - e2 * s * s) + math.log((1.0 + e * s) / (1.0 - e * s)) / (2.0 * e))
    )


def band_area(
    lat0: float, lat1: float, width: float, a: float = WGS84_A, f: float = WGS84_F
) -> float:
    """m² between latitudes `lat0` < `lat1`, `width` degrees of longitude wide."""
    return math.radians(width) * (_authalic(lat1, a, f) - _authalic(lat0, a, f))


def geo_meta(
    rows: int,
    cols: int,
    *,
    lon0: float,
    lat0: float,
    dx: float,
    dy: float,
    epsg: int = 4326,
    nodata: float | None = None,
) -> RasterMeta:
    """A point-registered grid in degrees, node (0, 0) at `(lon0, lat0)`."""
    return RasterMeta(
        x_min=lon0,
        y_max=lat0,
        delta_x=dx,
        delta_y=dy,
        cols=cols,
        rows=rows,
        epsg=epsg,
        geographic=True,
        nodata=nodata,
        nodata_source="absent" if nodata is None else "tag",
        pixel_is_area=False,
        vertical_unit_assumed=True,
    )


def geo_tile(
    z: npt.NDArray[np.floating],
    *,
    lon0: float,
    lat0: float,
    dx: float,
    dy: float,
    nodata: float | None = None,
) -> DemTile:
    rows, cols = z.shape
    meta = geo_meta(rows, cols, lon0=lon0, lat0=lat0, dx=dx, dy=dy, nodata=nodata)
    return DemTile(meta=meta, array=z)


def to_4326(x: float, y: float, crs: str = gf.EPSG) -> tuple[float, float]:
    lon, lat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(x, y)
    return float(lon), float(lat)


# ---------------------------------------------------------------- the twin

TWIN_CENTRE = to_4326(*gf.lat(gf.CC, gf.RIM_ROW))
TWIN_KX, TWIN_KY = metres_per_degree(TWIN_CENTRE[1])
TWIN_DX, TWIN_DY = gf.D / TWIN_KX, gf.D / TWIN_KY
TWIN_LON0 = TWIN_CENTRE[0] - gf.CC * TWIN_DX
TWIN_LAT0 = TWIN_CENTRE[1] + gf.RIM_ROW * TWIN_DY


def twin_xy(col: float, row: float) -> tuple[float, float]:
    """A lattice position (col, row) on the twin, in EPSG:4326."""
    return (TWIN_LON0 + col * TWIN_DX, TWIN_LAT0 - row * TWIN_DY)


def twin_tile(z: npt.NDArray[np.floating], nodata: float | None = None) -> DemTile:
    return geo_tile(z, lon0=TWIN_LON0, lat0=TWIN_LAT0, dx=TWIN_DX, dy=TWIN_DY, nodata=nodata)


def twin_column_line(col: float, row0: float, row1: float) -> tuple[tuple[float, float], ...]:
    """`gf.column_line` on the twin: a straight line down a column, downstream."""
    return (twin_xy(col, row0), twin_xy(col, row1))


def lattice_of(meta: RasterMeta, xy: Any) -> list[tuple[int, int]]:
    """Points of `meta`'s CRS as the (row, col) of the nodes they are."""
    pts = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
    r, c = meta.index_of(pts[:, 0], pts[:, 1])
    return [(round(float(a)), round(float(b))) for a, b in zip(r, c, strict=True)]
