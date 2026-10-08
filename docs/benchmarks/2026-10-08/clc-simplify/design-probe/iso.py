import json, shapely, border_apsc as m
c = json.load(open("overshoot_chain.json"))["chains"][0]
src = [tuple(p) for p in c["source"]]
x0, y0 = src[0]
src = [(round(x - x0, 2), round(y - y0, 2)) for x, y in src]  # local origin, cm
for un in (True, False):
    m.UNANCHORED = un
    ch = m.Chain(list(src), False, False, False)
    m.simplify([ch], 50.0)
    o = m.walk(ch)
    L1, L2 = shapely.LineString(src), shapely.LineString(o)
    f = float(shapely.distance(shapely.points(shapely.get_coordinates(shapely.segmentize(L1, 0.1))), L2).max())
    b = float(shapely.distance(shapely.points(shapely.get_coordinates(shapely.segmentize(L2, 0.1))), L1).max())
    print("unanchored" if un else "anchored", len(o), "vertices; source->simplified", round(f, 2), "simplified->source", round(b, 2))
json.dump({"band_m": 50.0, "note": "German fused case, CORINE border, local origin, rounded to 1 cm", "source": src}, open("overshoot_fixture.json", "w"))
