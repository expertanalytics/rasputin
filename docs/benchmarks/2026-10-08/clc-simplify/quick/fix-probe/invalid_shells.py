exec(open("fixprobe.py").read().split("for name in sys.argv")[0])
import collections
for name in ("band50.pkl", "clip_nosimp.pkl", "band0.pkl"):
    polys, outline, d = pickle.load(open(name, "rb"))
    path = fi._Outline(outline, d); reasons = collections.Counter(); inval_in = 0; ex = []
    for pi, p in enumerate(polys):
        inval_in += not p.is_valid
        for q in shapely.get_parts(p):
            rs = [orig_ring(path, shapely.get_coordinates(r)) for r in shapely.get_rings(q)]
            if all(x is None for x in rs): continue
            cl = [shapely.get_coordinates(r) if x is None else x[0] for r, x in zip(shapely.get_rings(q), rs)]
            if len(cl[0]) < 3: continue
            s = shapely.Polygon(cl[0], [h for h in cl[1:] if len(h) >= 3])
            if not s.is_valid:
                why = shapely.is_valid_reason(s); reasons[why.split("[")[0]] += 1
                if len(ex) < 4: ex.append((pi, len(cl) - 1, why, round(q.area), round(s.area)))
    print(name, "invalid input polygons", inval_in, dict(reasons)); [print("   ", e) for e in ex]
