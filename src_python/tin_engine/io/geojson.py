"""The catchment file's bytes: one polygon as GeoJSON (increments 22 and 29).

`docs/increments/29-nve-reference-catchments.md`, "The batch" ("The GeoJSON
writer moves"). Both `rasputin catchment` and `rasputin station-catchments`
write a catchment file, which `mesh --domain` reads back; this module builds
its bytes and opens nothing, as `io/ply.py` and `io/vtk_legacy.py` leave the
writing to their caller.
"""

from __future__ import annotations

import json
from collections.abc import Mapping

from shapely.geometry import Polygon


def catchment_geojson(polygon: Polygon, crs: str, properties: Mapping[str, object]) -> bytes:
    """UTF-8 JSON of a FeatureCollection naming `crs`, with one Feature: the
    polygon's exterior ring and `properties`, as increment 22 wrote it."""
    doc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": crs}},
        "features": [
            {
                "type": "Feature",
                "properties": dict(properties),
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[list(xy) for xy in polygon.exterior.coords]],
                },
            }
        ],
    }
    return json.dumps(doc).encode("utf-8")


__all__ = ["catchment_geojson"]
