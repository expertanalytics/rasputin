"""Edge statistics and buffer time of an outline, moved to a projected CRS first.
Usage: python outline_probe.py <outline.geojson> <its CRS> <projected CRS> <d>
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import shapely
from pyproj import Transformer
from shapely.geometry import shape

path, src, dst, d = Path(sys.argv[1]), sys.argv[2], sys.argv[3], float(sys.argv[4])
doc = json.loads(path.read_text())
feats = doc.get("features", [doc])
g = shapely.union_all([shape(f["geometry"]) for f in feats])
t = Transformer.from_crs(src, dst, always_xy=True)
g = shapely.transform(g, lambda xy: np.column_stack(t.transform(xy[:, 0], xy[:, 1])))
parts = shapely.get_parts(g)
p = max(parts, key=lambda q: q.area)
pxy = shapely.get_coordinates(p.exterior)
v = np.diff(pxy, axis=0)
seg = np.hypot(*v.T)
axis = np.mean((np.abs(v[:, 0]) < 1e-6) | (np.abs(v[:, 1]) < 1e-6))
print(
    f"{path.name}: {len(parts)} parts, {len(shapely.get_coordinates(g))} vertices, "
    f"largest part {len(pxy)} vertices, area {p.area / 1e6:.0f} km2, "
    f"edge median {np.median(seg):.1f} m, axis-parallel edges {100 * axis:.0f} %"
)
t0 = time.perf_counter()
p.buffer(d, join_style="mitre")
print(f"  buffer({d}, mitre) of the largest part: {time.perf_counter() - t0:.3f} s")
