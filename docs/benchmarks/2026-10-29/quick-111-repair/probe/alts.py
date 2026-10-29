"""Time ways to find the invalid polygons near the region (medians of 5)."""
import sys, time, statistics, numpy as np, shapely
from concurrent.futures import ThreadPoolExecutor
MARGIN = 100.0
d = np.load(sys.argv[1], allow_pickle=True)
shapes = shapely.from_wkb(d["wkb"]); xy = d["xy"]
box = (*(xy.min(axis=0) - MARGIN), *(xy.max(axis=0) + MARGIN))
pool = ThreadPoolExecutor(8)
def chunks(a, k=64):
    return np.array_split(a, min(k, len(a)))
def threaded(f, a):
    return np.concatenate(list(pool.map(f, chunks(a))))
def A(): return ~shapely.is_valid(shapely.clip_by_rect(shapes, *box))
def B(): return threaded(lambda c: ~shapely.is_valid(shapely.clip_by_rect(c, *box)), shapes)
def E(): return ~shapely.is_simple(shapely.boundary(shapely.clip_by_rect(shapes, *box)))
def W(): return ~shapely.is_valid(shapes)
def WT(): return threaded(lambda c: ~shapely.is_valid(c), shapes)
def CLIP(): return shapely.clip_by_rect(shapes, *box)
def ISV_ONLY(c=shapely.clip_by_rect(shapes, *box)): return ~shapely.is_valid(c)
def time_of(f, n=5):
    ts = []
    for _ in range(n):
        t0 = time.perf_counter(); r = f(); ts.append(time.perf_counter() - t0)
    return statistics.median(ts), min(ts), max(ts), int(np.sum(r)) if r.dtype == bool else None
nv = shapely.get_num_coordinates(shapes)
print(f"{sys.argv[1].split('/')[-1]}: {len(shapes)} polygons, {nv.sum()} vertices, max {nv.max()}")
for name, f in [("A clip+is_valid (branch)", A), ("  clip only", CLIP), ("  is_valid of clipped only", ISV_ONLY),
                ("B clip+is_valid, 8 threads", B), ("E clip+is_simple(boundary)", E),
                ("W whole is_valid (master shares)", W), ("WT whole is_valid, 8 threads", WT)]:
    n = 3 if name.startswith("W") else 5
    m, lo, hi, bad = time_of(f, n)
    print(f"{name:36s} {m:.3f} s [{lo:.3f}-{hi:.3f}]  flagged {bad}")
