"""Does keeping the box-clipped polygon (what the check already computes) pay
for the check in the clean-up's first step, intersection with the region?"""
import sys, time, statistics, numpy as np, shapely
d = np.load(sys.argv[1], allow_pickle=True)
shapes = shapely.from_wkb(d["wkb"]); xy = d["xy"]
region = shapely.Polygon(xy); shapely.prepare(region)
box = (*(xy.min(axis=0) - 100.0), *(xy.max(axis=0) + 100.0))
def med(f, n=3):
    ts = []
    for _ in range(n):
        t0 = time.perf_counter(); r = f(); ts.append(time.perf_counter() - t0)
    return statistics.median(ts), r
t_raw, a = med(lambda: [shapely.intersection(g, region) for g in shapes])
clipped = shapely.clip_by_rect(shapes, *box)
t_clip, b = med(lambda: [shapely.intersection(g, region) for g in clipped])
t_c, _ = med(lambda: shapely.clip_by_rect(shapes, *box))
t_v, _ = med(lambda: shapely.is_valid(clipped))
same = sum(abs(x.area - y.area) <= 1e-6 * max(x.area, 1) for x, y in zip(a, b))
print(f"{sys.argv[1].split('/')[-1]}: region vertices {len(xy)}")
print(f"  clean-up intersection, raw polygons:     {t_raw:.3f} s")
print(f"  clean-up intersection, clipped polygons: {t_clip:.3f} s  (areas equal for {same}/{len(a)})")
print(f"  clip to box {t_c:.3f} s, is_valid of clipped {t_v:.3f} s")
