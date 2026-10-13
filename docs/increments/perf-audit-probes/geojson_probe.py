"""Reading a big GeoJSON FeatureCollection: master's way and two vectorised ones.
Usage: python geojson_probe.py <features.geojson>
"""

import json
import sys
import time
from pathlib import Path

import shapely
from shapely.geometry import shape

path = Path(sys.argv[1])


def t(label, f):
    t0 = time.perf_counter()
    r = f()
    print(f"  {label:<58} {time.perf_counter() - t0:7.3f} s")
    return r


raw = t("read bytes", path.read_bytes)
doc = t("json.loads (stdlib)", lambda: json.loads(raw))
feats = doc["features"]
a = t("shape() per feature (master)", lambda: [shape(f["geometry"]) for f in feats])
b = t("from_geojson on the whole text (GEOS reader)", lambda: shapely.from_geojson(raw.decode()))
parts = shapely.get_parts(b)
print("  same geometries:", len(parts) == len(a) and all(shapely.equals_exact(parts, a, 0)))
texts = t("json.dumps per geometry", lambda: [json.dumps(f["geometry"]) for f in feats])
c = t("from_geojson over the texts (vectorised)", lambda: shapely.from_geojson(texts))
print("  same geometries:", bool(shapely.equals_exact(c, a, 0).all()))
vertices = len(shapely.get_coordinates(a))
print(f"  {len(feats)} features, {vertices} vertices, {len(raw) / 1e6:.0f} MB")
