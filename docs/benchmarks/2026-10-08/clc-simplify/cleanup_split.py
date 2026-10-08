"""The 20c-3 land-cover clean-up's steps, timed one by one, and how far
`shapely.coverage_simplify` moves each class at each tolerance.

    /Users/skavhaug/projects/rasputin/.venv/bin/python cleanup_split.py > cleanup_split.txt

A re-enactment of `feature_input._clean` (repair 1 m, same-class merge,
coverage_simplify with simplify_boundary=False), not a call into it: the
region here is the fused outline buffered by 1 km, where rasputin clips to its
own read region, so the times are indicative. The outline rule is left out.
Distance moved: the Hausdorff distance between a class's polygon before and
after simplification (boundaries only), in metres; vertices after.
"""

import json
import time

import numpy as np
import shapely
from pyproj import Transformer
from shapely.geometry import shape

OUTLINE = "/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/fused_catchment_glo30.geojson"
CORINE = "/Users/skavhaug/projects/rasputin_data/germany_corine/isar_loisach_clc2018_3035.geojson"

t = time.perf_counter()
dom = shape(json.load(open(OUTLINE))["features"][0]["geometry"])
region = dom.buffer(1000)
to = Transformer.from_crs(3035, 25832, always_xy=True)
feats = json.load(open(CORINE))["features"]
geoms = [shapely.transform(shape(f["geometry"]), lambda xy: np.column_stack(to.transform(xy[:, 0], xy[:, 1]))) for f in feats]
codes = [int(f["properties"]["Code_18"]) for f in feats]
print(f"read + move: {time.perf_counter() - t:.3f} s, {len(geoms)} polygons")


def step(name, fn, *a, **k):
    t0 = time.perf_counter()
    out = fn(*a, **k)
    print(f"{name}: {time.perf_counter() - t0:.3f} s")
    return out


clipped = step("clip to region", lambda: [shapely.intersection(g, region) for g in geoms])
keep = [i for i, g in enumerate(clipped) if not g.is_empty and g.area > 0]
polys = np.array([clipped[i] for i in keep], dtype=object)
codes = [codes[i] for i in keep]
print(f"vertices after clip: {int(shapely.get_num_coordinates(polys).sum())}")
polys = step("coverage_clean 1 m", shapely.coverage_clean, polys, snapping_distance=1.0, gap_width=1.0, merge_strategy="min_area")
groups: dict[int, list[int]] = {}
for i, c in enumerate(codes):
    groups.setdefault(c, []).append(i)
merged = step("merge same class", lambda: np.array([shapely.coverage_union_all(polys[g]) for g in groups.values()], dtype=object))
print(f"vertices after merge: {int(shapely.get_num_coordinates(merged).sum())}")
for ft in (10, 30, 100):
    s = step(f"coverage_simplify {ft} m", shapely.coverage_simplify, merged, ft, simplify_boundary=False)
    hd = [shapely.hausdorff_distance(a.boundary, b.boundary) for a, b in zip(merged, s)]
    da = [100 * (b.area - a.area) / a.area for a, b in zip(merged, s)]
    print(f"  vertices {int(shapely.get_num_coordinates(s).sum())}; largest class move {max(hd):.1f} m;"
          f" median class move {np.median(hd):.1f} m; class area change {min(da):+.2f} % to {max(da):+.2f} %")
