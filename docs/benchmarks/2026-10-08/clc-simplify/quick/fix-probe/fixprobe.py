# Probe of fixes F1 (no cuts on edges lying on one outline segment) and F2 (a ring with one kept input vertex measured by symmetric difference)
import pickle, sys, time, inspect, textwrap, numpy as np, shapely
import tin_engine.feature_input as fi
orig_ring, orig_loops = fi._Outline.ring, fi._loops
def on_outline(self, xy):
    a, b = shapely.points(xy[:-1]), shapely.points(xy[1:])
    i, j = self.tree.query(shapely.linestrings(np.stack([xy[:-1], xy[1:]], 1)), predicate="dwithin", distance=fi.IN_LINE)
    segs = self.tree.geometries[j]
    ok = (shapely.distance(a[i], segs) <= fi.IN_LINE) & (shapely.distance(b[i], segs) <= fi.IN_LINE)
    out = np.zeros(len(xy) - 1, bool); out[i[ok]] = True
    return out
src = textwrap.dedent(inspect.getsource(orig_ring)).replace(
    "near[self.tree.query(edges, predicate=\"dwithin\", distance=self.d)[0]] = True",
    "near[self.tree.query(edges, predicate=\"dwithin\", distance=self.d)[0]] = True\n    cut = near & ~on_outline(self, xy)")
src = src.replace("if near[k]:", "if cut[k]:")
assert src.count("cut[k]") == 1 and "on_outline(self" in src
ns = {**fi.__dict__, "on_outline": on_outline}; exec(src, ns); new_ring = ns["ring"]
lsrc = textwrap.dedent(inspect.getsource(orig_loops)).replace("if not at:", "if len(at) < 2:")
ns2 = dict(fi.__dict__); exec(lsrc, ns2); new_loops = ns2["_loops"]
def run(polys, outline, d, f1, f2):
    fi._Outline.ring = new_ring if f1 else orig_ring
    fi._loops = new_loops if f2 else orig_loops
    t = time.perf_counter(); r = fi.snap_to_outline(polys, outline, d); return r, time.perf_counter() - t
for name in sys.argv[1:]:
    polys, outline, d = pickle.load(open(name, "rb"))
    r0, t0 = run(polys, outline, d, False, False)
    r1, t1 = run(polys, outline, d, True, False)
    r2, t2 = run(polys, outline, d, True, True)
    same_lines = all(len(a) == len(b) and all(x.equals_exact(y, 0) for x, y in zip(a, b)) for a, b in zip(r0.lines, r1.lines))
    same_area = max(abs(a.area - b.area) for a, b in zip(r0.polygons, r1.polygons))
    same_poly = all(a.equals(b) for a, b in zip(r0.polygons, r1.polygons))
    lost = sum(shapely.difference(shapely.intersection(o, outline), shapely.intersection(n, outline)).area for o, n in zip(polys, r2.polygons))
    print(f"{name}: today {t0:.2f}s area {r0.area_changed:.0f} | F1 {t1:.2f}s lines bit-identical {same_lines}, polygons equal {same_poly} (largest area diff {same_area:.2e} m2) | F1+F2 {t2:.2f}s area {r2.area_changed:.0f} (real change: {lost:.0f})")
