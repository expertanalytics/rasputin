"""P grown by d = P united with its boundary grown by d. Grow the boundary in short,
overlapping pieces (GEOS buffer is fast on short lines), then unite. Compared with
GEOS's direct buffer of the polygon.
Usage: python chunked_buffer.py <domain.geojson, projected> <d> <mitre|round>
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import shape

g = shape(json.loads(Path(sys.argv[1]).read_text())["features"][0]["geometry"])
d = float(sys.argv[2])
join = sys.argv[3]  # mitre or round
xy = shapely.get_coordinates(g.exterior)  # closed
n = len(xy) - 1


def chunked(k: int) -> shapely.Geometry:
    starts = np.arange(0, n, k)
    lines = [shapely.linestrings(xy[s : min(s + k + 2, n + 1)]) for s in starts]  # 1 edge overlap
    if starts[-1] + k + 2 > n + 1:  # wrap: the last piece runs past the start vertex
        lines[-1] = shapely.linestrings(np.vstack([xy[starts[-1] : n + 1], xy[1:3]]))
    grown = shapely.buffer(np.array(lines), d, join_style=join, cap_style="flat")
    return shapely.union_all(np.concatenate([[g], grown]))


t0 = time.perf_counter()
exact = g.buffer(d, join_style=join)
print(f"{n} vertices, d = {d}, join {join}: GEOS buffer {time.perf_counter() - t0:.2f} s")
for k in (250, 1000, 4000):
    t0 = time.perf_counter()
    c = chunked(k)
    tc = time.perf_counter() - t0
    diff = c.symmetric_difference(exact).area
    far = c.exterior.hausdorff_distance(exact.exterior)
    print(
        f"  pieces of {k:5d} edges: {tc:6.2f} s; symmetric difference with GEOS's {diff:.3g} m2 "
        f"of {exact.area / 1e6:.0f} km2; Hausdorff {far:.3g} m; "
        f"bounds equal {c.bounds == exact.bounds}"
    )
