"""R3 (e) follow-up: where the population-3 length with the merge on lies, and which polygons own it.
usage: r3e_where.py OUT_CRS -- <mesh args>"""
import sys, numpy as np, shapely, pyproj
from shapely.ops import transform
import tin_engine.feature_input as fi
from tin_engine.cli import app
crs = sys.argv[1]; args = sys.argv[sys.argv.index("--") + 1:]
calls = []; orig = fi.snap_to_outline
def wrap(p, o, d):
    r = orig(p, o, d); calls.append((r, o)); return r
fi.snap_to_outline = wrap
try: app(["mesh", *args], standalone_mode=False)
except SystemExit: pass
to = pyproj.Transformer.from_crs(crs, "EPSG:3035", always_xy=True).transform
r, D = calls[0]; D35 = transform(to, D); P35 = [transform(to, g) for g in r.polygons]
for i, ls in enumerate(r.lines):
    for line in ls:
        c = shapely.get_coordinates(transform(to, line)); a, b = c[:-1], c[1:]
        for which, k, v in (("N", 1, 3811923.31), ("E", 0, 4585680.19)):
            on = (np.abs(a[:, k] - v) <= 1) & (np.abs(b[:, k] - v) <= 1)
            for s, e in zip(a[on], b[on]):
                seg = shapely.intersection(shapely.LineString([s, e]), D35)
                if seg.length == 0: continue
                mid = shapely.Point(seg.interpolate(0.5, normalized=True)); side = [j for j, g in enumerate(P35) if g.distance(mid) < 1e-6]
                print(f"{which} line: polygon {i}, {seg.length:.3f} m from {s.round(2)} to {e.round(2)}; polygons touching its middle {side}")
