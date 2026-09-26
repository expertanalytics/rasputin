"""The domain polygon: what ``--domain`` reads (increment 16, R1).

One ``Polygon``, holes allowed, from GeoJSON (a bare geometry, a ``Feature`` or
a ``FeatureCollection`` of exactly one) or WKT, by suffix. Vertices are used as
given: no densifying, simplifying or snapping. The polygon is oriented outer
counter-clockwise and holes clockwise (increment 3's winding contract).

CRS is required and never transformed. A GeoJSON file's is its ``crs`` member,
and a file without one is EPSG:4326 by RFC 7946; WKT has none, so it comes from
the caller. The must-match rule is this increment's scope only (U1, and the
user's caveat on it), so it lives in one replaceable function, :func:`check_crs`.

Every vertex must lie in the DEM's node rectangle, border included (U4 a):
outside it there is no bilinear z (R0).

Pure: json, shapely, pyproj and ``RasterMeta``. No ``_core``, no typer.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import shapely
import shapely.wkt
from pydantic import BaseModel, ConfigDict
from pyproj import CRS
from pyproj.exceptions import CRSError
from shapely.geometry import Polygon, shape
from shapely.geometry.polygon import orient
from shapely.validation import explain_validity

from tin_engine.io.models import RasterMeta

GEOJSON_SUFFIXES = (".geojson", ".json")
WKT_SUFFIXES = (".wkt",)
GEOJSON_DEFAULT_EPSG = 4326  # RFC 7946: no crs member means WGS 84


class DomainError(ValueError):
    """A domain file refused, in words for the person who wrote it."""


class DomainPolygon(BaseModel):
    """The domain as read: an oriented shapely polygon and its EPSG code."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    polygon: Polygon
    epsg: int


def read_domain(path: Path, meta: RasterMeta, crs: str | None = None) -> DomainPolygon:
    """Read, check and orient the domain in ``path`` against the DEM ``meta``.

    ``crs`` is the ``--domain-crs`` text: required for WKT, and for GeoJSON it
    must agree with the file's own. Raises :class:`DomainError` on any refusal.
    """
    suffix = path.suffix.lower()
    if suffix not in GEOJSON_SUFFIXES + WKT_SUFFIXES:
        raise DomainError(
            f"{path.name}: unknown suffix {suffix or '(none)'}; use "
            f"{', '.join(GEOJSON_SUFFIXES + WKT_SUFFIXES)}"
        )
    try:
        text = path.read_text()
    except (OSError, UnicodeDecodeError) as exc:
        raise DomainError(f"cannot read {path}: {exc}") from exc
    try:
        if suffix in WKT_SUFFIXES:
            if crs is None:
                raise DomainError(f"{path.name}: WKT carries no CRS; give it with --domain-crs")
            geometry, epsg = shapely.wkt.loads(text), _epsg(crs)
        else:
            geometry, epsg = _geojson(json.loads(text))
            if crs is not None and _epsg(crs) != epsg:
                raise DomainError(
                    f"{path.name} is in EPSG:{epsg} but --domain-crs says EPSG:{_epsg(crs)}"
                )
    except DomainError:
        raise
    except (ValueError, TypeError, KeyError, AttributeError, shapely.errors.GEOSException) as exc:
        raise DomainError(f"cannot parse {path.name}: {exc}") from exc

    if not isinstance(geometry, Polygon):
        raise DomainError(f"{path.name}: the domain must be one Polygon, got {geometry.geom_type}")
    if geometry.is_empty:
        raise DomainError(f"{path.name}: the domain polygon is empty")
    if not geometry.is_valid:
        raise DomainError(f"{path.name}: invalid polygon, {explain_validity(geometry)}")
    check_crs(epsg, meta)
    _check_extent(geometry, meta)
    return DomainPolygon(polygon=orient(geometry, sign=1.0), epsg=epsg)


def check_crs(epsg: int, meta: RasterMeta) -> None:
    """U1 (a): the domain's CRS must be the DEM's, since nothing reprojects yet."""
    if epsg != meta.epsg:
        raise DomainError(
            f"the domain is in EPSG:{epsg} and the DEM in EPSG:{meta.epsg}; they must match "
            f"(a GeoJSON file without a crs member is EPSG:{GEOJSON_DEFAULT_EPSG})"
        )


def _epsg(text: str) -> int:
    try:
        code = CRS.from_user_input(text).to_epsg()
    except CRSError as exc:
        raise DomainError(f"cannot read the CRS {text!r}: {exc}") from exc
    if code is None:
        raise DomainError(f"the CRS {text!r} has no EPSG code")
    return int(code)


def _geojson(doc: dict[str, Any]) -> tuple[Any, int]:
    """The geometry of a bare geometry, a Feature or a one-feature collection,
    and the EPSG code of the document's ``crs`` member."""
    member = doc.get("crs")
    epsg = _epsg(member["properties"]["name"]) if member is not None else GEOJSON_DEFAULT_EPSG
    if doc.get("type") == "FeatureCollection":
        features = doc["features"]
        if len(features) != 1:
            raise DomainError(f"need exactly one feature, got {len(features)}")
        doc = features[0]
    if doc.get("type") == "Feature":
        doc = doc["geometry"]
    return shape(doc), epsg


def _check_extent(polygon: Polygon, meta: RasterMeta) -> None:
    """U4 (a): every vertex in the node rectangle, as the core's ``cell_of`` has it."""
    x_max = meta.x_min + (meta.cols - 1) * meta.delta_x
    y_min = meta.y_max - (meta.rows - 1) * meta.delta_y
    for ring in (polygon.exterior, *polygon.interiors):
        for x, y in ring.coords:
            outside = max(meta.x_min - x, x - x_max, y_min - y, y - meta.y_max)
            if outside > 0:
                raise DomainError(
                    f"vertex ({x}, {y}) is {outside:g} m outside the DEM's node rectangle "
                    f"x {meta.x_min} .. {x_max}, y {y_min} .. {meta.y_max}"
                )
