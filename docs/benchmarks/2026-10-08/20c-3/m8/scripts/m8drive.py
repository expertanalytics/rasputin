"""M8: run `rasputin mesh` in-process, wrapping shapely.coverage_clean (the repair)
and snap_to_outline (for the domain), and measure the repair as M7 did.
usage: m8drive.py TAG -- <mesh args>"""
import sys, time, numpy as np, shapely
import tin_engine.feature_input as fi
from tin_engine.cli import app

tag = sys.argv[1]; args = sys.argv[sys.argv.index("--") + 1:]
cleans, doms = [], []
orig_clean, orig_snap = shapely.coverage_clean, fi.snap_to_outline
def clean(polys, **kw):
    t = time.perf_counter(); out = orig_clean(polys, **kw); cleans.append((polys, out, kw, time.perf_counter() - t)); return out
def snap(polygons, outline, distance):
    doms.append(outline); return orig_snap(polygons, outline, distance)
shapely.coverage_clean = clean; fi.snap_to_outline = snap
t0 = time.perf_counter()
try: app(["mesh", *args], standalone_mode=False)
except SystemExit: pass
print(f"{tag}: whole run in-process {time.perf_counter() - t0:.2f} s")

def short(C, lim=0.1):
    n = 0
    for g in C:
        for p in shapely.get_parts(g):
            for r in (p.exterior, *p.interiors):
                xy = shapely.get_coordinates(r); n += int((np.hypot(*np.diff(xy, axis=0).T) < lim).sum())
    return n
def gaps(C, dom):
    U = shapely.coverage_union_all(C)
    hs = [shapely.Polygon(r) for p in shapely.get_parts(U) for r in p.interiors]
    hs = [h for h in hs if dom.intersects(h)]
    return [(round(2 * shapely.maximum_inscribed_circle(h).length, 4), round(h.area, 2)) for h in hs]
dom = doms[0]
for polys, out, kw, dt in cleans:
    moved = sum(shapely.symmetric_difference(a, b).area for a, b in zip(polys, out))
    print(f"{tag}: coverage_clean {kw} {dt:.2f} s; polygons {len(polys)}; vertices "
          f"{shapely.get_num_coordinates(polys).sum()} -> {shapely.get_num_coordinates(out).sum()}; "
          f"land cover moved {moved:.2f} m2; ring edges < 10 cm {short(polys)} -> {short(out)}; "
          f"emptied {sum(g.is_empty for g in out)}")
    print(f"{tag}: gaps touching the domain (width m, area m2) before {gaps(polys, dom)} after {gaps(out, dom)}")
    print(f"{tag}: valid after {shapely.coverage_is_valid(out)}")
    np.save(f"{sys.argv[0].rsplit('/', 1)[0]}/m8-{tag}-cover.npy",
            np.array([shapely.to_wkb(polys), shapely.to_wkb(out)], dtype=object), allow_pickle=True)
    np.save(f"{sys.argv[0].rsplit('/', 1)[0]}/m8-{tag}-dom.npy", np.array([shapely.to_wkb(dom)], dtype=object), allow_pickle=True)
