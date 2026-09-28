"""Write the 4326 test catchment for the 16b-1/2 acceptance (Ola's case).

Not a real catchment: a smooth irregular 48-vertex ring in lon/lat, centred on
the corner where DTM10 tiles 6602_1, 6602_2, 6603_3 and 6603_4 meet
(300 000, 6 650 000 in EPSG:25833, 11.4198 E 59.9387 N), so the DEM is a
four-tile mosaic, inside both DTM10_UTM33_20260925 and the Norway CORINE
extract's 254 tiles. Radii are chosen for a few hundred km². Deterministic.
"""
import json, math, sys
from pathlib import Path
from pyproj import Transformer
from shapely.geometry import Polygon
from shapely.ops import transform

LON0, LAT0 = 11.419790, 59.938743
ring = []
for k in range(48):
    a = 2 * math.pi * k / 48
    r = 1.0 + 0.18 * math.sin(3 * a) + 0.10 * math.cos(5 * a + 0.7)
    ring.append((round(LON0 + 0.17 * r * math.cos(a), 6), round(LAT0 + 0.085 * r * math.sin(a), 6)))
ring.append(ring[0])
poly = Polygon(ring)
assert poly.is_valid
utm = transform(Transformer.from_crs("EPSG:4326", "EPSG:25833", always_xy=True).transform, poly)
doc = {"type": "FeatureCollection",
       "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::4326"}},
       "features": [{"type": "Feature", "properties": {"name": "16b12-acceptance catchment",
                     "area_km2_in_25833": round(utm.area / 1e6, 1)},
                     "geometry": {"type": "Polygon", "coordinates": [ring]}}]}
Path(sys.argv[1]).write_text(json.dumps(doc, indent=1) + "\n")
print(f"area {utm.area / 1e6:.1f} km2, bounds 25833 {[round(b) for b in utm.bounds]}")
