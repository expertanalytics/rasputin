"""Run the branch CLI, catching coverage_gap's inputs; then the largest part's
inscribed-circle centre at the code's tolerance, the default, and with no tile cuts."""
import sys, shapely
import tin_engine.land_cover_map as lcm
from tin_engine import cli
seen = {}
orig = lcm.coverage_gap
def spy(domain, covers, *, threads=1):
    seen["domain"], seen["covers"] = domain, list(covers)
    return orig(domain, covers, threads=threads)
lcm.coverage_gap = spy
if hasattr(cli, "coverage_gap"): cli.coverage_gap = spy
from tin_engine.plugins import main
sys.argv = ["rasputin"] + sys.argv[1:]
try: main()
except SystemExit: pass
dom, cov = seen["domain"], seen["covers"]
whole = shapely.difference(dom, shapely.union_all(cov))
tiled = lcm._polygonal(shapely.union_all([g for g in []])) if False else None
parts = sorted(shapely.get_parts(whole), key=lambda p: -p.area)
big = parts[0]
print("untiled parts", len(parts), "areas", [round(p.area) for p in parts], "sum", round(sum(p.area for p in parts), 1))
for tol in (lcm.GAP_NOISE / 10, None):
    c = shapely.maximum_inscribed_circle(big, tol) if tol else shapely.maximum_inscribed_circle(big)
    print(f"untiled, tolerance {tol}: centre {c.coords[0][0]:.2f} {c.coords[0][1]:.2f} radius {c.length:.3f}")
g = orig(dom, cov, threads=8)
print("code's parts:", [(round(a), round(x, 2), round(y, 2)) for a, x, y in g.parts], "sum", round(g.area, 1))
for x, y in ((136216, 6830714), (g.parts[0][1], g.parts[0][2])):
    print(f"distance to the gap's edge from {x:.2f} {y:.2f}: {shapely.distance(shapely.Point(x, y), big.boundary):.3f}")
