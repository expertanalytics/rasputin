"""Which clipped polygons the validity check spends its time on."""
import sys, time, numpy as np, shapely
d = np.load(sys.argv[1], allow_pickle=True)
shapes = shapely.from_wkb(d["wkb"]); xy = d["xy"]
box = (*(xy.min(axis=0) - 100.0), *(xy.max(axis=0) + 100.0))
c = shapely.clip_by_rect(shapes, *box)
ts = np.empty(len(c))
for i, g in enumerate(c):
    t0 = time.perf_counter(); shapely.is_valid(g); ts[i] = time.perf_counter() - t0
o = np.argsort(ts)[::-1]
nv = shapely.get_num_coordinates(c); nh = shapely.get_num_interior_rings(shapely.get_parts(c)) if False else None
print(f"total {ts.sum():.3f} s; top 5 {ts[o[:5]].sum():.3f} s, top 20 {ts[o[:20]].sum():.3f} s")
for i in o[:6]:
    parts = shapely.get_parts(c[i]); holes = int(shapely.get_num_interior_rings(parts).sum())
    print(f"  {ts[i]:.4f} s  vertices {nv[i]}  parts {len(parts)}  holes {holes}  code {d['codes'][i]}")
# bigger than 10k vertices
big = nv > 10000
print(f"polygons over 10k clipped vertices: {big.sum()}, their time {ts[big].sum():.3f} s")
