"""GeoJSON: the one reading path and its `crs` rule, and the catchment file's bytes.

One reading path for every GeoJSON input: :func:`read_collection` turns a
parsed document into its features and the text of its CRS, for `--domain`,
`--land-cover`, `--prefer`, `--lines`, `catchment --lakes` and the station, reference, river and NVE
lake files. The `crs` member is
GeoJSON 2008's named CRS; a file without one is RFC 7946's WGS 84 where the
caller allows a default, and a null member is refused whatever the default.
A `Feature` is one feature, a collection its features, and anything else a
geometry wrapped as one feature with no properties.

Both `rasputin catchment` and `rasputin station-catchments` write a catchment
file, which `mesh --domain` reads back; this module builds
its bytes and opens nothing, as `io/ply.py` and `io/vtk_legacy.py` leave the
writing to their caller.
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import shapely
from shapely.geometry import Polygon, shape

from tin_engine.crs import parse_crs

GEOJSON_SUFFIXES = (".geojson", ".json")
RFC7946_CRS = "EPSG:4326"  # RFC 7946: no crs member means WGS 84


def read_collection(doc: object, *, default_crs: str | None) -> tuple[list[dict[str, Any]], str]:
    """The features of one GeoJSON object and the text of its `crs` member.
    Opens nothing; a refusal is a ValueError in words without the file's
    name, which the caller adds."""
    if not isinstance(doc, dict):
        raise ValueError(
            "not a GeoJSON object; the file must hold a FeatureCollection, a Feature or a geometry"
        )
    if "crs" not in doc:
        if default_crs is None:
            raise ValueError("no crs member; the file must name its CRS")
        crs = default_crs
    elif doc["crs"] is None:
        raise ValueError("the crs member is null; the file must name its CRS")
    else:
        props = doc["crs"].get("properties") if isinstance(doc["crs"], dict) else None
        name = props.get("name") if isinstance(props, dict) else None
        if name is None:
            raise ValueError("the crs member has no name; it must name the CRS")
        crs = str(name)
    parse_crs(crs)
    if doc.get("type") == "Feature":
        return [doc], crs
    if doc.get("type") != "FeatureCollection" and "features" not in doc:
        return [{"type": "Feature", "properties": {}, "geometry": doc}], crs
    features = doc.get("features")
    if not isinstance(features, list) or not all(isinstance(f, dict) for f in features):
        raise ValueError("no features list; the file is not a FeatureCollection")
    return list(features), crs


@dataclass(frozen=True, slots=True)
class CollectionSummary:
    """A GeoJSON object's header, as `rasputin sources` lists it: the
    geometry type names present, the feature count, the attribute names in
    first-seen order, the bounds in its CRS (None with no geometry), whether
    it has a `crs` member, and the CRS a reader taking RFC 7946's default
    reads it in."""

    geometry: frozenset[str]
    features: int
    attributes: tuple[str, ...]
    extent: tuple[float, float, float, float] | None
    crs_stated: bool
    crs: str
    rows: list[dict[str, Any]]


def collection_summary(doc: object) -> CollectionSummary:
    """The header of one parsed GeoJSON object; refused as `read_collection` refuses."""
    rows, crs = read_collection(doc, default_crs=RFC7946_CRS)
    geometries = [f["geometry"] for f in rows if isinstance(f.get("geometry"), dict)]
    attributes = dict.fromkeys(k for f in rows for k in (f.get("properties") or {}))
    x0, y0, x1, y1 = map(float, shapely.total_bounds([shape(g) for g in geometries]))
    return CollectionSummary(
        geometry=frozenset(str(g.get("type")) for g in geometries),
        features=len(rows),
        attributes=tuple(attributes),
        extent=None if math.isnan(x0) else (x0, y0, x1, y1),
        crs_stated=isinstance(doc, dict) and "crs" in doc,
        crs=crs,
        rows=rows,
    )


def feature_collection(crs: str, features: list[dict[str, Any]]) -> dict[str, Any]:
    """{"type": "FeatureCollection", "crs": <named crs>, "features": features}."""
    return {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": crs}},
        "features": features,
    }


def catchment_geojson(polygon: Polygon, crs: str, properties: Mapping[str, object]) -> bytes:
    """UTF-8 JSON of a FeatureCollection naming `crs`, with one Feature: the
    polygon's exterior ring (holes are not written) and `properties`."""
    ring = [[list(xy) for xy in polygon.exterior.coords]]
    feature = {
        "type": "Feature",
        "properties": dict(properties),
        "geometry": {"type": "Polygon", "coordinates": ring},
    }
    return json.dumps(feature_collection(crs, [feature])).encode("utf-8")


__all__ = [
    "GEOJSON_SUFFIXES",
    "RFC7946_CRS",
    "CollectionSummary",
    "catchment_geojson",
    "collection_summary",
    "feature_collection",
    "read_collection",
]
