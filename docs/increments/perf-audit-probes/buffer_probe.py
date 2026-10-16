"""Time the domain buffers the pipeline makes, in isolation, and cheaper equivalents.
Usage: python buffer_probe.py <domain.geojson, projected> <grid spacing h in metres>
"""

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import shape

path = Path(sys.argv[1])
g = shape(json.loads(path.read_text())["features"][0]["geometry"])
if g.geom_type == "MultiPolygon":
    g = max(g.geoms, key=lambda p: p.area)
xy = shapely.get_coordinates(g)
seg = np.hypot(*np.diff(xy, axis=0).T)
print(
    f"{path.name}: {len(xy)} vertices, valid {g.is_valid}, area {g.area / 1e6:.0f} km2, "
    f"edge median {np.median(seg):.1f} m, min {seg.min():.3g} m, geos {shapely.geos_version_string}"
)


def t(label, f):
    t0 = time.perf_counter()
    r = f()
    print(f"  {label:<62} {time.perf_counter() - t0:8.3f} s")
    return r


h = 31 if len(sys.argv) < 3 else float(sys.argv[2])
d = math.sqrt(2) * h


def mitre(p):
    return p.buffer(d, join_style="mitre")


a = t(f"buffer({d:.1f}, mitre) [target_grid_for, dem_input:213]", lambda: mitre(g))
b = t("buffer(100) round [feature_input.source_region]", lambda: g.buffer(100.0))
t("buffer(100) round, quad_segs=2", lambda: g.buffer(100.0, quad_segs=2))
t("buffer(100) of exterior ring only, round", lambda: g.exterior.buffer(100.0))
t("convex_hull", lambda: g.convex_hull)
c = t("buffer(100) of convex_hull, round", lambda: g.convex_hull.buffer(100.0))
t(f"buffer({d:.1f}) mitre of simplify(1 m)", lambda: mitre(shapely.simplify(g, 1.0)))
t("simplify(0) (collinear vertices dropped)", lambda: shapely.simplify(g, 0.0))
print("  simplify(0) vertices:", len(shapely.get_coordinates(shapely.simplify(g, 0.0))))
s0 = t(f"buffer({d:.1f}, mitre) of simplify(0)", lambda: mitre(shapely.simplify(g, 0.0)))
print("  identical with simplify(0) first:", s0.equals_exact(a, 0))
# hull of buffer vs buffer of hull (feature_input only needs the hull)
diff = shapely.convex_hull(b).symmetric_difference(c).area
print(f"  hull(buffer(P,100)) vs buffer(hull(P),100): sym diff {diff:.1f} m2")
