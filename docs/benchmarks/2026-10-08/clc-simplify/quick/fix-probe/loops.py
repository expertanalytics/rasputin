import pickle, sys, numpy as np, shapely
import tin_engine.feature_input as fi
polys, outline, d = pickle.load(open(sys.argv[1], "rb"))
path = fi._Outline(outline, d)
rows = []
for pi, p in enumerate(polys):
    for qi, q in enumerate(shapely.get_parts(p)):
        for ri, ring in enumerate(shapely.get_rings(q)):
            xy = shapely.get_coordinates(ring)
            done = path.ring(xy)
            if done is None: continue
            at = [i for i in done[2] if i >= 0]
            for L in fi._loops(xy, *done):
                a = shapely.intersection(L, outline).area
                if a > 0: rows.append((a, pi, qi, ri, len(xy), len(at), Polygon_area := shapely.area(q), L.centroid.coords[0]))
rows.sort(reverse=True)
tot = sum(r[0] for r in rows)
print(sys.argv[1], "loops with area", len(rows), "sum (not unioned)", round(tot))
for r in rows[:12]:
    print(f"  {r[0]:10.0f} m2 poly {r[1]} part {r[2]} ring {r[3]} ring-verts {r[4]} kept-input {r[5]} part-area {r[6]:.0f} at {r[7][0]:.0f},{r[7][1]:.0f}")
print("  share of top 12:", round(sum(r[0] for r in rows[:12]) / tot, 3), " loops > 1000 m2:", sum(r[0] > 1000 for r in rows), "their sum", round(sum(r[0] for r in rows if r[0] > 1000)))
