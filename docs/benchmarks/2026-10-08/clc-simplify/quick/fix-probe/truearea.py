import pickle, sys, numpy as np, shapely
import tin_engine.feature_input as fi
for name in sys.argv[1:]:
    polys, outline, d = pickle.load(open(name, "rb"))
    r = fi.snap_to_outline(polys, outline, d)
    old = [shapely.intersection(p, outline) for p in polys]
    new = [shapely.intersection(p, outline) for p in r.polygons]
    # area that changed polygon: each moved piece leaves one polygon and joins another, so half the sum of losses+gains
    lost = sum(shapely.difference(o, n).area for o, n in zip(old, new))
    gained = sum(shapely.difference(n, o).area for o, n in zip(old, new))
    print(f"{name}: reported {r.area_changed:.0f} m2; lost by its polygon {lost:.0f} m2, gained {gained:.0f} m2; class totals inside: old {sum(o.area for o in old):.0f} new {sum(n.area for n in new):.0f}")
    # the largest loop's part: is it still there?
    p = shapely.get_parts(polys[10])[1]
    print("   poly10 part1 area", round(p.area), "still covered by poly10 after:", round(shapely.intersection(p, r.polygons[10]).area))
