"""The domain file: what ``--domain`` reads (increment 16, R1; audit PR C).

One ``Polygon``, holes allowed, from GeoJSON (a bare geometry, a ``Feature`` or
a ``FeatureCollection`` of exactly one) or WKT, by suffix. Vertices are used as
given: no densifying, simplifying or snapping. The polygon is oriented outer
counter-clockwise and holes clockwise (increment 3's winding contract).

CRS is required. A GeoJSON file's is its ``crs`` member, read by
:func:`tin_engine.io.geojson.read_collection`, and a file without one is
EPSG:4326 by RFC 7946; WKT has none, so it comes from the caller.

The file is opened by ``io/repository.py``: GeoJSON through ``read_json``, as
``--features`` and the station files are (a UTF-8 byte order mark is
skipped), WKT through ``read_text``, as UTF-8. Moved here from
``tin_engine.domain`` (``docs/increments/python-audit.md``, section 11), which
keeps the domain's types, ``to_crs`` and ``check_extent``.
"""

from __future__ import annotations

from pathlib import Path

import shapely
import shapely.wkt
from pyproj import CRS
from shapely.geometry import Polygon, shape
from shapely.geometry.polygon import orient
from shapely.validation import explain_validity

from tin_engine.crs import parse_crs, same_crs
from tin_engine.domain import DomainError, DomainPolygon

from .geojson import GEOJSON_SUFFIXES, RFC7946_CRS, read_collection
from .repository import read_json, read_text

WKT_SUFFIXES = (".wkt",)


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
        raw = read_text(path) if suffix in WKT_SUFFIXES else read_json(path)
    except (OSError, UnicodeDecodeError) as exc:
        raise DomainError(f"cannot read {path}: {exc}") from exc
    except ValueError as exc:
        raise DomainError(f"cannot parse {path.name}: {exc}") from exc
    try:
        if suffix in WKT_SUFFIXES:
            if crs is None:
                raise DomainError(f"{path.name}: WKT carries no CRS; give it with --domain-crs")
            geometry, own = shapely.wkt.loads(raw), crs
            _parsed(own)
        else:
            features, own = read_collection(raw, default_crs=RFC7946_CRS)
            if len(features) != 1:
                raise DomainError(f"need exactly one feature, got {len(features)}")
            geometry = shape(features[0]["geometry"])
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


__all__ = ["WKT_SUFFIXES", "read_domain"]
