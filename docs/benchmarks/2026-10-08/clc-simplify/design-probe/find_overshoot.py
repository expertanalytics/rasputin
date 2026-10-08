"""Find the border whose simplification by increment 22's check overshoots."""
import json, sys
import numpy as np, shapely
import border_apsc as m
src = open("border_apsc.py").read()
# reuse main's preparation by running it with a hook: simplest is to re-run the steps
m.UNANCHORED = True
captured = {}
orig = m.simplify
def hook(chains, band):
    r = orig(chains, band)
    for c in chains:
        if all(c.fixed) or len(c.fine) < 4:
            continue
        o = m.walk(c)
        L1 = shapely.LineString(c.fine + ([c.fine[0]] if c.closed else []))
        L2 = shapely.LineString(o + ([o[0]] if c.closed else []))
        pts = shapely.points(shapely.get_coordinates(shapely.segmentize(L1, 0.5)))
        d = float(shapely.distance(pts, L2).max())
        if d > band + 1e-6:
            captured.setdefault("chains", []).append({"closed": c.closed, "source": c.fine, "simplified": o, "overshoot_m": d})
    return r
m.simplify = hook
sys.argv = ["x", "50", "/dev/null", "unanchored"]
m.main()
json.dump(captured, open(sys.argv_out if hasattr(sys, "argv_out") else "overshoot_chain.json", "w"))
for c in captured["chains"]:
    print(len(c["source"]), "source vertices,", len(c["simplified"]), "simplified, overshoot", round(c["overshoot_m"], 2), "closed", c["closed"])
