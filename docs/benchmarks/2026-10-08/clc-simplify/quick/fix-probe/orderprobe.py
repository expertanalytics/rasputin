# Option R: today's base steps through the outline rule, then the domain clip, the simplifier, and the linework less what lies on the outline
import pickle, time, numpy as np, shapely
import tin_engine.feature_input as fi
from tin_engine.border_simplify import simplify_borders
polys, outline, d = pickle.load(open("band0.pkl", "rb"))   # region-clipped, repaired, merged: base's input to the rule
t = time.perf_counter(); r = fi.snap_to_outline(polys, outline, d); t_rule = time.perf_counter() - t
t = time.perf_counter(); clipped = [fi._polygonal(shapely.intersection(p, outline)) for p in r.polygons]; t_clip = time.perf_counter() - t
keep = [c for c in clipped if c is not None]
print("valid after rule+clip:", all(c.is_valid for c in keep), " coverage valid:", bool(shapely.coverage_is_valid(np.array(keep, dtype=object))))
t = time.perf_counter(); res = simplify_borders(keep, 50.0); t_simp = time.perf_counter() - t
print("simplifier counts:", res.counts)
segs = shapely.linestrings(np.concatenate([np.stack([x[:-1], x[1:]], 1) for x in [shapely.get_coordinates(g) for g in shapely.get_rings(outline)]]))
tree = shapely.STRtree(segs)
t = time.perf_counter(); nlines = 0; nedges = 0; non = 0
for p in res.polygons:
    for q in shapely.get_parts(p):
        for ring in shapely.get_rings(q):
            xy = shapely.get_coordinates(ring)
            a, b = shapely.points(xy[:-1]), shapely.points(xy[1:])
            i, j = tree.query(shapely.linestrings(np.stack([xy[:-1], xy[1:]], 1)), predicate="dwithin", distance=fi.IN_LINE)
            ok = (shapely.distance(a[i], segs[j]) <= fi.IN_LINE) & (shapely.distance(b[i], segs[j]) <= fi.IN_LINE)
            on = np.zeros(len(xy) - 1, bool); on[i[ok]] = True
            nlines += len(fi._chains(xy[:-1], on.tolist())); nedges += len(on); non += int(on.sum())
t_lines = time.perf_counter() - t
V = sum(int(shapely.get_num_coordinates(g)) for g in res.polygons)
print(f"rule {t_rule:.2f}s (area {r.area_changed:.0f} m2), domain clip {t_clip:.2f}s, simplify {t_simp:.2f}s, linework {t_lines:.2f}s; {V} polygon vertices after; {nedges} ring edges, {non} on the outline, {nlines} lines")
pickle.dump((res.polygons, outline), open("optionR.pkl", "wb"))
