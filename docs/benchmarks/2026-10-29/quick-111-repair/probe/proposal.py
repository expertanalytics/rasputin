"""The proposed form against the branch's, on the dumped polygons (medians of 5):
check = clip to the box + is_valid; then the clean-up's intersection with the region."""
import sys, time, statistics, numpy as np, shapely
from concurrent.futures import ThreadPoolExecutor
d = np.load(sys.argv[1], allow_pickle=True)
shapes = shapely.from_wkb(d["wkb"]); xy = d["xy"]
region = shapely.Polygon(xy); shapely.prepare(region)
box = (*(xy.min(axis=0) - 100.0), *(xy.max(axis=0) + 100.0))
pool = ThreadPoolExecutor(8)
nv = shapely.get_num_coordinates(shapes)
def balanced(k=32):  # index groups of about equal vertex count, biggest first
    order = np.argsort(nv)[::-1]; loads = np.zeros(k); groups = [[] for _ in range(k)]
    for i in order:
        j = int(np.argmin(loads)); groups[j].append(i); loads[j] += nv[i]
    return [np.array(g) for g in groups if g]
groups = balanced()
def branch():
    bad = ~shapely.is_valid(shapely.clip_by_rect(shapes, *box))
    return bad, [shapely.intersection(g, region) for g in shapes]
def proposed():
    def one(ix):
        c = shapely.clip_by_rect(shapes[ix], *box); return ix, c, ~shapely.is_valid(c)
    clipped = np.empty(len(shapes), dtype=object); bad = np.zeros(len(shapes), bool)
    for ix, c, b in pool.map(one, groups):
        clipped[ix] = c; bad[ix] = b
    return bad, [shapely.intersection(g, region) for g in clipped]
def med(f, n=5):
    ts = []
    for _ in range(n):
        t0 = time.perf_counter(); r = f(); ts.append(time.perf_counter() - t0)
    return statistics.median(ts), min(ts), max(ts), r
def raw():
    return None, [shapely.intersection(g, region) for g in shapes]
m0 = med(raw); m1 = med(branch); m2 = med(proposed)
eq = all(shapely.equals_exact(a, b, 0) or shapely.equals(a, b) for a, b in zip(m1[3][1], m2[3][1]))
print(f"{sys.argv[1].split('/')[-1]}:")
for name, m in (("master-like: intersection only (no check)", m0),
                ("branch: check + intersection", m1), ("proposed: threaded check, clipped kept", m2)):
    print(f"  {name:44s} {m[0]:.3f} s [{m[1]:.3f}-{m[2]:.3f}]")
print(f"  flagged: branch {int(m1[3][0].sum())}, proposed {int(m2[3][0].sum())}; clean-up inputs equal: {eq}")
