"""GeoJSON: the one reading path and its `crs` rule, and the catchment file's bytes.

Audit PR C (`docs/increments/python-audit.md`, section 11, "The one reading
path"): :func:`read_collection` turns a parsed document into its features
and the text of its CRS, for `--domain`, `--features`, `catchment --lakes`
and the station, reference, river and NVE lake files. The `crs` member is
GeoJSON 2008's named CRS; a file without one is RFC 7946's WGS 84 where the
caller allows a default, and a null member is refused whatever the default.
A `Feature` is one feature, a collection its features, and anything else a
geometry wrapped as one feature with no properties.

`docs/increments/29-nve-reference-catchments.md`, "The batch" ("The GeoJSON
writer moves"). Both `rasputin catchment` and `rasputin station-catchments`
write a catchment file, which `mesh --domain` reads back; this module builds
its bytes and opens nothing, as `io/ply.py` and `io/vtk_legacy.py` leave the
writing to their caller.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

from shapely.geometry import Polygon

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


def feature_collection(crs: str, features: list[dict[str, Any]]) -> dict[str, Any]:
    """{"type": "FeatureCollection", "crs": <named crs>, "features": features}."""
    return {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": crs}},
        "features": features,
    }


def catchment_geojson(polygon: Polygon, crs: str, properties: Mapping[str, object]) -> bytes:
    """UTF-8 JSON of a FeatureCollection naming `crs`, with one Feature: the
    polygon's exterior ring and `properties`, as increment 22 wrote it."""
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
    "catchment_geojson",
    "feature_collection",
    "read_collection",
]
