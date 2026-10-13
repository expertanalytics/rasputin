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
import threading
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict
from pyproj import CRS, Proj, Transformer, get_ellps_map
from pyproj.crs import ProjectedCRS
from pyproj.enums import WktVersion
from pyproj.exceptions import CRSError, ProjError

from tin_engine.io.models import CUSTOM_SRS_ID, UNDEFINED_SRS, SrsEntry
from tin_engine.run_record import escaped_ascii

Xy = npt.NDArray[np.float64]
#: Increment 35, 3.1 (b): an array of `2 * CHUNK_POINTS` points or more is
#: split into chunks of at least `CHUNK_POINTS` (a chunk of work, about 60 ms
#: of PROJ on one machine) and mapped on up to `threads` threads, each with its
#: own `Transformer` (38: the count is the caller's; one worker by default).
CHUNK_POINTS = 250_000


def parse_crs(text: str | CRS) -> CRS:
    """Anything `CRS.from_user_input` accepts, or a `ValueError` naming it."""
    if isinstance(text, CRS):
        return text
    try:
        return CRS.from_user_input(text)
    except CRSError as exc:
        raise ValueError(f"cannot read the CRS {text!r}: {exc}") from exc


def reprojector(src: str | CRS, dst: str | CRS, *, threads: int = 1) -> Callable[[Any], Xy]:
    """A map from an `(N, 2)` array-like of `(x, y)` in `src` to an `(N, 2)`
    float64 array in `dst`: pyproj's `always_xy` transform, bit for bit.
    A point with no image comes back as `inf`; the caller decides. With
    `threads` above 1, a large array is mapped in chunks on up to that many
    threads (`CHUNK_POINTS`), one `Transformer` per thread since pyproj's are
    not thread-safe; a projection is applied point by point, so the bits are
    those of one call."""
    local = threading.local()
    local.transformer = _transformer(src, dst)

    def one(points: Xy) -> Xy:
        if not hasattr(local, "transformer"):
            local.transformer = _transformer(src, dst)
        x, y = local.transformer.transform(points[:, 0], points[:, 1])
        return np.column_stack([np.asarray(x, np.float64), np.asarray(y, np.float64)])

    def apply(xy: Any) -> Xy:
        points = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        if threads < 2 or len(points) < 2 * CHUNK_POINTS:
            return one(points)
        chunks = np.array_split(points, min(threads, len(points) // CHUNK_POINTS))
        with ThreadPoolExecutor(len(chunks)) as pool:
            return np.concatenate(list(pool.map(one, chunks)))

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


def same_crs(a: str | CRS, b: str | CRS) -> bool:
    """Whether coordinates in `a` and in `b` name the same points: PROJ calls
    them equivalent once both are in x-then-y order, and the operation it
    builds between them is its `noop` (`docs/increments/python-audit.md`,
    section 9). No operation is False."""
    t = _equivalent_xy(a, b)
    return t is not None and str(t.definition).startswith("proj=noop")


def epsg_matches(crs: CRS) -> tuple[int, ...]:
    """The EPSG codes PROJ offers for `crs`, at every confidence, that sit on
    `crs`'s own base geographic CRS and that PROJ calls equivalent to `crs`
    once both are in x-then-y order (the first leg of `same_crs`), lowest
    first. Empty when there are none (`docs/increments/geotiff-crs-by-parameters.md`,
    section 4). A candidate whose axes carry a meridian (EPSG's polar CRSs)
    is also taken when the operation to it is PROJ's `noop`."""
    base = _base_code(crs)
    offered = {int(m.code) for m in crs.list_authority(auth_name="EPSG", min_confidence=0)}
    return tuple(
        code
        for code in sorted(offered)
        if (_base_code(candidate := CRS.from_epsg(code)) == base and _equivalent_xy(crs, candidate))
        or _polar_noop(crs, candidate)
    )


def _polar_noop(crs: CRS, candidate: CRS) -> bool:
    """`candidate`'s axes carry a meridian, its base is `crs`'s in any axis
    order, and the operation from `crs` to it is `noop`."""
    axes = candidate.to_json_dict().get("coordinate_system", {}).get("axis", [])
    a, b = crs.geodetic_crs, candidate.geodetic_crs
    if not any("meridian" in axis for axis in axes) or a is None or b is None:
        return False
    if not a.equals(b, ignore_axis_order=True):
        return False
    try:
        return str(_transformer(crs, candidate).definition).startswith("proj=noop")
    except ProjError:
        return False


def metre_twin(crs: CRS) -> tuple[CRS, float]:
    """`crs` with its axes in metres and its linear conversion parameters
    (false easting and northing, easting and northing at an origin) converted
    to metres, named "<name>, in metres"; and the metres per unit of its axes.
    `crs` itself and 1.0 when its axes are already in metres."""
    factor = float(crs.axis_info[0].unit_conversion_factor)
    if factor == 1.0:
        return crs, 1.0
    data = crs.to_json_dict()
    data.pop("id", None)
    data["name"] = f"{crs.name}, in metres"
    data["conversion"].pop("id", None)
    for parameter in data["conversion"]["parameters"]:
        unit = parameter.get("unit")
        if isinstance(unit, dict) and unit.get("type") == "LinearUnit":
            parameter["value"] *= unit["conversion_factor"]
            parameter["unit"] = "metre"
    for axis in data["coordinate_system"]["axis"]:
        axis["unit"] = "metre"
    return CRS.from_json_dict(data), factor


def _base_code(crs: CRS) -> int | None:
    return None if crs.geodetic_crs is None else crs.geodetic_crs.to_epsg(min_confidence=100)


def _equivalent_xy(a: str | CRS, b: str | CRS) -> Transformer | None:
    """The transformer from `a` to `b` when PROJ calls its two x-then-y CRSs
    equivalent (`same_crs`'s first leg), else None; no operation is None."""
    try:
        t = _transformer(a, b)
    except ProjError:
        return None
    src, dst = t.source_crs, t.target_crs
    return t if src is not None and dst is not None and src.equals(dst) else None


def transform_label(src: str | CRS, dst: str | CRS) -> str:
    """'none' when `same_crs(src, dst)`, else `transform_description(src, dst)`."""
    return "none" if same_crs(src, dst) else transform_description(src, dst)


def single_crs(texts: Iterable[str], refusal: type[ValueError] = ValueError) -> str:
    """The one CRS text among `texts` (a DEM's tiles' `meta.crs`), or
    `refusal` naming them all."""
    found = sorted(set(texts))
    if len(found) != 1:
        raise refusal(
            f"the DEM files are in {len(found)} different CRSs ({', '.join(found)}); "
            "all must be in one CRS"
        )
    return found[0]


def crs_label(crs: str | CRS) -> str:
    """`EPSG:n` when pyproj finds an exact EPSG code (so OGC's CRS84 is not
    `EPSG:4326`), otherwise ASCII text pyproj parses back to the same CRS."""
    parsed = parse_crs(crs)
    code = parsed.to_epsg(min_confidence=100)
    if code is not None:
        return f"EPSG:{code}"
    return escaped_ascii(parsed.to_string())


def srs_entry(crs: str | None) -> SrsEntry:
    """The `gpkg_spatial_ref_sys` row for `crs` (increment 36): an EPSG code is
    its own `srs_id` with WKT1_GDAL (what QGIS and GDAL read); a CRS pyproj
    parses but cannot code is 100000 with its WKT2; no CRS, or text pyproj
    refuses, is the reserved undefined cartesian row -1. Never raises."""
    try:
        parsed = parse_crs(crs) if crs else None
    except ValueError:
        parsed = None
    if parsed is None:
        return UNDEFINED_SRS
    code = parsed.to_epsg()
    if code is None:
        wkt2 = parsed.to_wkt()
        return SrsEntry(
            srs_id=CUSTOM_SRS_ID,
            organization="rasputin",
            code=CUSTOM_SRS_ID,
            name=parsed.name,
            definition=wkt2,
        )
    wkt1 = parsed.to_wkt(version=WktVersion.WKT1_GDAL)
    return SrsEntry(srs_id=code, organization="EPSG", code=code, name=parsed.name, definition=wkt1)


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


class Site(BaseModel):
    """A point's latitude, longitude and grid convergence, in degrees."""

    model_config = ConfigDict(frozen=True)

    lat_deg: float
    lon_deg: float
    #: The true azimuth of grid north (pyproj's meridian convergence).
    convergence_deg: float


def site_of(crs: str | CRS, x: float, y: float) -> Site:
    """Where `(x, y)` in the projected `crs` lies on the globe; a geographic
    CRS is a ValueError."""
    c = parse_crs(crs)
    if not c.is_projected:
        raise ValueError("a projected CRS in metres is needed; this one is in degrees")
    lon, lat = _transformer(c, c.geodetic_crs or "EPSG:4326").transform(x, y)
    factors = Proj(c).get_factors(lon, lat)
    return Site(lat_deg=lat, lon_deg=lon, convergence_deg=factors.meridian_convergence)
