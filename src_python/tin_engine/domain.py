"""The domain polygon: what ``--domain`` reads (increment 16, R1).

One ``Polygon``, holes allowed, from GeoJSON (a bare geometry, a ``Feature`` or
a ``FeatureCollection`` of exactly one) or WKT, by suffix. Vertices are used as
given: no densifying, simplifying or snapping. The polygon is oriented outer
counter-clockwise and holes clockwise (increment 3's winding contract).

CRS is required. A GeoJSON file's is its ``crs`` member, and a file without one
is EPSG:4326 by RFC 7946; WKT has none, so it comes from the caller. The domain
keeps its own CRS, and :meth:`DomainPolygon.to_crs` moves it into the DEM's
(increment 15b, ``15-dem-mosaic.md`` R9), replacing 16's must-match rule.

Every vertex must lie in the DEM's node rectangle, border included (U4 a):
outside it there is no bilinear z (R0). That is :func:`check_extent`, run after
the transform, against the mosaic.

Pure: json, numpy, shapely, pyproj, pydantic, ``tin_engine.crs`` and
``RasterMeta``. No ``_core``, no typer.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import shapely
import shapely.wkt
from pydantic import BaseModel, ConfigDict
from pyproj import CRS
from shapely.geometry import Polygon, shape
from shapely.geometry.polygon import orient
from shapely.validation import explain_validity

from tin_engine.crs import parse_crs, reprojector, same_crs
from tin_engine.io.models import RasterMeta

GEOJSON_SUFFIXES = (".geojson", ".json")
WKT_SUFFIXES = (".wkt",)
GEOJSON_DEFAULT_CRS = "EPSG:4326"  # RFC 7946: no crs member means WGS 84


class DomainError(ValueError):
    """A domain file refused, in words for the person who wrote it."""


class DomainPolygon(BaseModel):
    """The domain: an oriented shapely polygon, and the text of its CRS."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    polygon: Polygon
    crs: str

    def to_crs(self, dst: str | CRS) -> DomainPolygon:
        """The domain in ``dst``: each vertex through pyproj's ``always_xy``
        transform, once, and re-oriented; rings stay straight (R9). The same CRS,
        by ``same_crs``, returns the same polygon labelled ``target.to_string()``.
        """
        source, target = parse_crs(self.crs), parse_crs(dst)
        if same_crs(source, target):
            return self.model_copy(update={"crs": target.to_string()})
        move = reprojector(source, target)
        rings = [np.asarray(r.coords) for r in (self.polygon.exterior, *self.polygon.interiors)]
        moved = [move(r) for r in rings]
        for ring, image in zip(rings, moved, strict=True):
            bad = ~np.isfinite(image).all(axis=1)
            if bad.any():
                x, y = ring[bad][0]
                raise DomainError(
                    f"the domain vertex ({x}, {y}) in {self.crs} has no image in "
                    f"{target.to_string()}"
                )
        polygon = orient(Polygon(moved[0], moved[1:]), sign=1.0)
        return DomainPolygon(polygon=polygon, crs=target.to_string())


def read_domain(path: Path, crs: str | None = None) -> DomainPolygon:
    """Read, check and orient the domain in ``path``, in its own CRS.

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
            geometry, own = shapely.wkt.loads(text), crs
            _parsed(own)
        else:
            geometry, own = _geojson(json.loads(text))
            if crs is not None and not same_crs(_parsed(crs), _parsed(own)):
                raise DomainError(f"{path.name} is in {own} but --domain-crs says {crs}")
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
    return DomainPolygon(polygon=orient(geometry, sign=1.0), crs=own)


def _parsed(text: str) -> CRS:
    try:
        return parse_crs(text)
    except ValueError as exc:
        raise DomainError(str(exc)) from exc


def _geojson(doc: dict[str, Any]) -> tuple[Any, str]:
    """The geometry of a bare geometry, a Feature or a one-feature collection,
    and the CRS text of the document's ``crs`` member, checked."""
    member = doc.get("crs")
    own = str(member["properties"]["name"]) if member is not None else GEOJSON_DEFAULT_CRS
    _parsed(own)
    if doc.get("type") == "FeatureCollection":
        features = doc["features"]
        if len(features) != 1:
            raise DomainError(f"need exactly one feature, got {len(features)}")
        doc = features[0]
    if doc.get("type") == "Feature":
        doc = doc["geometry"]
    return shape(doc), own


def check_extent(domain: DomainPolygon, meta: RasterMeta) -> None:
    """U4 (a): every vertex in ``meta``'s node rectangle, as the core's
    ``cell_of`` has it, with ``domain`` already in the DEM's CRS."""
    polygon = domain.polygon
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
