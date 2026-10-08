"""The overshoot border alone: today's check against the anchored one, on the
committed overshoot_fixture.json (run from this folder)."""
import json, shapely, border_apsc as m
fx = json.load(open("overshoot_fixture.json"))
src, band = [tuple(p) for p in fx["source"]], fx["band_m"]
for un in (True, False):
    m.UNANCHORED = un
    ch = m.Chain(list(src), False, False, False)
    m.simplify([ch], band)
    o = m.walk(ch)
    L1, L2 = shapely.LineString(src), shapely.LineString(o)
    f = float(shapely.distance(shapely.points(shapely.get_coordinates(shapely.segmentize(L1, 0.1))), L2).max())
    b = float(shapely.distance(shapely.points(shapely.get_coordinates(shapely.segmentize(L2, 0.1))), L1).max())
    print("unanchored" if un else "anchored", len(o), "vertices; source->simplified", round(f, 2), "simplified->source", round(b, 2))
