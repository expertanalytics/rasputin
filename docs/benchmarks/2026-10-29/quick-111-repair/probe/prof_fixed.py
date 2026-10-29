"""Run `rasputin mesh` under cProfile; time land_cover_map._fixed's parts."""
import cProfile, pstats, sys, time, io
import numpy as np, shapely
from tin_engine import land_cover_map as lcm
from tin_engine.feature_input import MARGIN
orig = lcm._fixed
stats = []
def timed(pairs, xy):
    t0 = time.perf_counter()
    box = (*(xy.min(axis=0) - MARGIN), *(xy.max(axis=0) + MARGIN))
    shapes = np.array([g for _, g in pairs] or [shapely.Polygon()], dtype=object)
    t1 = time.perf_counter(); clipped = shapely.clip_by_rect(shapes, *box)
    t2 = time.perf_counter(); bad = ~shapely.is_valid(clipped)
    t3 = time.perf_counter(); full = ~shapely.is_valid(shapes)
    t4 = time.perf_counter()
    out = orig(pairs, xy); t5 = time.perf_counter()
    nv = shapely.get_num_coordinates(shapes); ncv = shapely.get_num_coordinates(clipped)
    inside = shapely.contains_properly(shapely.box(*box), shapes)
    stats.append(dict(n=len(pairs), clip=t2-t1, valid_clipped=t3-t2, valid_whole=t4-t3,
        fixed_total=t5-t4, bad_clipped=int(bad.sum()), bad_whole=int(full.sum()),
        vertices=int(nv.sum()), clipped_vertices=int(ncv.sum()), wholly_inside=int(inside.sum()),
        vmax=int(nv.max()), repaired=out[1]))
    return out
lcm._fixed = timed
from tin_engine.cli import app
pr = cProfile.Profile()
try:
    pr.enable(); app(sys.argv[2:], standalone_mode=False)
finally:
    pr.disable()
    pr.dump_stats(sys.argv[1])
    for s in stats: print("FIXED", s, file=sys.stderr)
