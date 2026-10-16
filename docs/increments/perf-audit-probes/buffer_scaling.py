"""How GEOS buffer time scales with the distance, on one staircase outline.
Usage: python buffer_scaling.py <domain.geojson, projected>

(The run that produced the audit's figures also timed synthetic staircase blobs;
they stayed under 2,009 vertices and 0.12 s and said nothing about scaling, so
that part is left out.)
"""

import json
import sys
import time
from pathlib import Path

import shapely
from shapely.geometry import shape

g = shape(json.loads(Path(sys.argv[1]).read_text())["features"][0]["geometry"])
if g.geom_type == "MultiPolygon":
    g = max(g.geoms, key=lambda p: p.area)
print(len(shapely.get_coordinates(g)), "vertices")
for d in (5.0, 10.0, 20.0, 43.8, 100.0):
    t0 = time.perf_counter()
    g.buffer(d, join_style="mitre")
    print(f"  mitre d={d:6.1f}: {time.perf_counter() - t0:7.3f} s")
