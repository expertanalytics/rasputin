exec(open("fixprobe.py").read().split("for name in sys.argv")[0])
polys, outline, d = pickle.load(open("band50.pkl", "rb"))
r0, _ = run(polys, outline, d, False, False); r1, _ = run(polys, outline, d, True, False)
for i, (a, b, o) in enumerate(zip(r0.polygons, r1.polygons, polys)):
    if abs(a.area - b.area) > 1:
        print(i, "input", round(o.area), "today", round(a.area), a.is_valid, len(shapely.get_parts(a)), "F1", round(b.area), b.is_valid, len(shapely.get_parts(b)))
        diff = shapely.symmetric_difference(a, b); parts = sorted(shapely.get_parts(diff), key=lambda g: -g.area)[:3]
        for g in parts: print("   diff part", round(g.area), g.centroid)
# validity of the rebuilt shells
path = fi._Outline(outline, d)
bad = [0, 0]
for p in polys:
    for q in shapely.get_parts(p):
        for k, f in enumerate((orig_ring, new_ring)):
            rs = [f(path, shapely.get_coordinates(r)) for r in shapely.get_rings(q)]
            if all(x is None for x in rs): continue
            cl = [shapely.get_coordinates(r) if x is None else x[0] for r, x in zip(shapely.get_rings(q), rs)]
            if len(cl[0]) >= 3 and not shapely.Polygon(cl[0], [h for h in cl[1:] if len(h) >= 3]).is_valid: bad[k] += 1
print("invalid rebuilt shells today / F1:", bad)
