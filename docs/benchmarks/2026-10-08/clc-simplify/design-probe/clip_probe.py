"""Design review round 1 of increment 32: time the domain clip against
today's region clip (the region approximated in EPSG:25832 as the domain's
convex hull plus 100 m), coverage_clean after each, and reduce_ring's time
per collapse on the domain-clipped rings at 50 m. Run with the repo venv;
output in clip_probe.txt."""
import json, time, numpy as np, shapely
from pyproj import Transformer
from shapely.geometry import shape, Polygon
from tin_engine._core import reduce_ring
from tin_engine.feature_input import _polygonal
OUT="/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/fused_catchment_glo30.geojson"
COR="/Users/skavhaug/projects/rasputin_data/germany_corine/isar_loisach_clc2018_3035.geojson"
dom = shape(json.load(open(OUT))["features"][0]["geometry"])
to = Transformer.from_crs(3035, 25832, always_xy=True)
geoms = [shapely.transform(shape(f["geometry"]), lambda xy: np.column_stack(to.transform(xy[:,0], xy[:,1]))) for f in json.load(open(COR))["features"]]
region = shapely.convex_hull(dom).buffer(100.0)
print("features", len(geoms), "domain vertices", shapely.get_num_coordinates(dom), "region vertices", shapely.get_num_coordinates(region))
def clip(target):
    t = time.perf_counter(); c = [_polygonal(shapely.intersection(g, target)) for g in geoms]; return time.perf_counter()-t, c
for name, target in (("region", region), ("domain", dom)):
    ts = []
    for _ in range(3):
        dt, c = clip(target); ts.append(dt)
    polys = np.array([g for g in c if g is not None], dtype=object)
    t = time.perf_counter(); cl = shapely.coverage_clean(polys, snapping_distance=1.0, gap_width=1.0, merge_strategy="min_area"); tc = time.perf_counter()-t
    print(name, "clip s", [round(x,3) for x in ts], "pieces", len(polys), "vertices", int(shapely.get_num_coordinates(polys).sum()), "coverage_clean s", round(tc,3))
# reduce_ring per collapse: every exterior of the domain-clipped pieces, band 50
rings = [np.asarray(p.exterior.coords)[:-1] for g in c if g is not None for p in shapely.get_parts(g)]
tot_t = tot_c = 0; nv = 0
for r in rings:
    if len(r) < 5: continue
    pr = Polygon(r)
    if not pr.exterior.is_ccw: r = r[::-1]
    r = r - r[0]
    t = time.perf_counter(); o = reduce_ring(np.ascontiguousarray(r), 50.0, np.empty((0,2))); tot_t += time.perf_counter()-t
    tot_c += o.collapses; nv += len(r)
print("reduce_ring rings", len(rings), "vertices", nv, "collapses", tot_c, "s", round(tot_t,3), "us/collapse", round(1e6*tot_t/max(tot_c,1),2))
