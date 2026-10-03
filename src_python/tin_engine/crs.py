"""CRS text and the one reprojection site (increment 15b).

`docs/increments/15-dem-mosaic.md` R2, R9 and I8. `_transformer` holds the only
`Transformer.from_crs` in `src_python/`, always with `always_xy=True`: x is
easting or longitude and y northing or latitude, whatever the CRS's own axis
order, which is the convention a GeoTIFF's model space and GeoJSON both use.

Every refusal is a `ValueError`. pyproj's `CRSError` is a `RuntimeError`, and
would escape every `except ValueError` that turns a refusal into a usage error.

Pure: pyproj and numpy. No paths, no `_core`; CRS never crosses into the core.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict
from pyproj import CRS, Proj, Transformer, get_ellps_map
from pyproj.crs import ProjectedCRS
from pyproj.exceptions import CRSError

Xy = npt.NDArray[np.float64]


def parse_crs(text: str | CRS) -> CRS:
    """Anything `CRS.from_user_input` accepts, or a `ValueError` naming it."""
    if isinstance(text, CRS):
        return text
    try:
        return CRS.from_user_input(text)
    except CRSError as exc:
        raise ValueError(f"cannot read the CRS {text!r}: {exc}") from exc


def reprojector(src: str | CRS, dst: str | CRS) -> Callable[[Any], Xy]:
    """A map from an `(N, 2)` array-like of `(x, y)` in `src` to an `(N, 2)`
    float64 array in `dst`: pyproj's `always_xy` transform, bit for bit.
    A point with no image comes back as `inf`; the caller decides."""
    transformer = _transformer(src, dst)

    def apply(xy: Any) -> Xy:
        points = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        x, y = transformer.transform(points[:, 0], points[:, 1])
        return np.column_stack([np.asarray(x, np.float64), np.asarray(y, np.float64)])

    return apply


def transform_bounds(
    src: str | CRS, dst: str | CRS, box: tuple[float, float, float, float]
) -> tuple[float, float, float, float]:
    """`(x_min, y_min, x_max, y_max)` in `src` to the box holding its image in
    `dst`: pyproj's `transform_bounds`, densified, `always_xy` (23a-2)."""
    x0, y0, x1, y1 = _transformer(src, dst).transform_bounds(*box, densify_pts=21)
    return x0, y0, x1, y1


def transform_description(src: str | CRS, dst: str | CRS) -> str:
    """What PROJ picked for `src` to `dst` (a datum shift included), for the record."""
    return str(_transformer(src, dst).description)


def transform_definition(src: str | CRS, dst: str | CRS) -> str:
    """PROJ's pipeline for `src` to `dst`, to see which steps it takes (16b R5)."""
    return str(_transformer(src, dst).definition)


def crs_label(crs: str | CRS) -> str:
    """`EPSG:n` when pyproj finds an exact EPSG code (so OGC's CRS84 is not
    `EPSG:4326`), otherwise ASCII text pyproj parses back to the same CRS."""
    parsed = parse_crs(crs)
    code = parsed.to_epsg(min_confidence=100)
    if code is not None:
        return f"EPSG:{code}"
    return parsed.to_string().encode("ascii", "backslashreplace").decode("ascii")


class CrsSuggestion(BaseModel):
    """A conformal CRS for a box (15c, D8): WKT2 on the DEM's datum, its family's name,
    and its worst point-scale and areal errors over the box."""

    model_config = ConfigDict(frozen=True)

    proj: str
    family: str
    proj4: str = ""  # the same projection as PROJ text, without the datum
    max_scale_error: float
    max_areal_error: float


def suggest_crs(box: tuple[float, float, float, float], geographic_crs: str | CRS) -> CrsSuggestion:
    """D8: the family by the box's `(W, E, S, N)` latitude and shape, on the
    DEM's own ellipsoid, its scale factor balanced over the box (`s = 2 /
    (kmin + kmax)` from a unit-scale evaluation), angles to 0.1 degrees."""
    west, east, south, north = box
    lon_c, lat_c = (west + east) / 2, (south + north) / 2
    el = parse_crs(geographic_crs).ellipsoid
    assert el is not None  # a geographic CRS has one
    a, rf = el.semi_major_metre, el.inverse_flattening
    names = [k for k, v in get_ellps_map().items() if v.get("a") == a and v.get("rf") == rf]
    ellps = f"+ellps={names[0]}" if names else f"+a={a!r} +rf={rf!r}"
    common = f"+lon_0={round(lon_c, 1)}"
    if abs(lat_c) >= 70 or max(abs(south), abs(north)) >= 80:
        family, head = "polar stereographic", f"+proj=stere +lat_0={math.copysign(90.0, lat_c)}"
    elif (east - west) * math.cos(math.radians(lat_c)) <= north - south:
        family, head = "transverse Mercator", "+proj=tmerc +lat_0=0"
    elif abs(lat_c) <= 15:
        family, head = "Mercator", "+proj=merc"
    else:
        sixth = (north - south) / 6
        family = "Lambert conformal conic"
        head = f"+proj=lcc +lat_1={round(south + sixth, 1)} +lat_2={round(north - sixth, 1)}"
        head += f" +lat_0={round(lat_c, 1)}"
    k = "+k" if family == "transverse Mercator" else "+k_0"
    lon, lat = np.meshgrid(np.linspace(west, east, 21), np.linspace(south, north, 21))

    def proj(s: float) -> tuple[str, float, float, float, float]:
        text = f"{head} {common} {k}={s!r} +x_0=0 +y_0=0 {ellps} +units=m +no_defs"
        f = Proj(text).get_factors(lon.ravel(), lat.ravel())
        scales = np.concatenate([f.meridional_scale, f.parallel_scale])
        areal = float(np.abs(np.asarray(f.areal_scale) - 1).max())
        worst = float(np.abs(scales - 1).max())
        return text, float(scales.min()), float(scales.max()), worst, areal

    _, kmin, kmax, _, _ = proj(1.0)
    text, _, _, scale_error, areal_error = proj(round(2 / (kmin + kmax), 6))
    # On the DEM's own datum, not its ellipsoid alone, so no ballpark transform
    # joins the two; WKT2, because no PROJ +datum names SIRGAS 2000.
    conversion = CRS(text).coordinate_operation
    wkt = ProjectedCRS(conversion, geodetic_crs=parse_crs(geographic_crs)).to_wkt()
    return CrsSuggestion(
        proj=wkt,
        proj4=text,
        family=family,
        max_scale_error=scale_error,
        max_areal_error=areal_error,
    )


def _transformer(src: str | CRS, dst: str | CRS) -> Transformer:
    return Transformer.from_crs(parse_crs(src), parse_crs(dst), always_xy=True)
