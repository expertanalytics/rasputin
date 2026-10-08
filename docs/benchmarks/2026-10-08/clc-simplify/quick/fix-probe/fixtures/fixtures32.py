"""Increment 32, section 15.3: the fixtures of tests 17-20b, measured at clearance 0
with the worktree's _core (border_collapse.hpp as at 633808b9). Each fixture: two
100 m squares, left [(0,0), A, B, C, D, (0,100)] and right [A, (200,0), (200,100), D, C, B],
A = (100,0) and D = (100,100) junctions, one collapse possible on A-B-C-D; optionally a
triangle island in one square (its ring and the hole's). See fx.py for the re-enacted
placements and the quantities of 15.2 (over every vertex and edge, a superset of what the
grid looks at)."""
from island import *
import subprocess, numpy as np
_hpp = "include/terrain/vector_simplify/border_collapse.hpp"
_root = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True).stdout.strip()
print("kernel source:", _hpp, "blob", subprocess.run(["git", "-C", _root, "hash-object", _hpp], capture_output=True, text=True).stdout.strip(), "(the _core built from it is the .venv's)")
# F12: the fold of fix-design review round 6. A and D junctions of three polygons each;
# the C-D placement E = (50, 0) lies beyond A, and E-D passes 0.5 m from A.
FA, FB, FC, FD = (0.0, 0.5), (0.0, -1.25), (-20.0, 0.0), (-100.0, 0.0)
FOLD = (FA, FD, [
    [FA, (0.0, 50.0), (-50.0, 50.0)],
    [FA, (-50.0, 50.0), (-100.0, 50.0), FD, FC, FB],
    [FA, FB, FC, FD, (-100.0, -50.0), (100.0, -50.0), (100.0, 50.0), (0.0, 50.0)],
])
print("fold coverage valid:", bool(shapely.coverage_is_valid(np.array([shapely.Polygon(r) for r in FOLD[2]], dtype=object))),
      "polygons valid:", all(shapely.Polygon(r).is_valid for r in FOLD[2]))
def r2(isl): return [(round(x, 2), round(y, 2)) for x, y in isl]
def ccw(isl): return isl if shapely.Polygon(isl).exterior.is_ccw else isl[::-1]
S4 = ((84.0, 10.0), (84.0, 55.0))
E4, _ = winner(*S4)
def vertex_island(d):
    for side in (1, -1):
        for P, Q in ((A, E4), (E4, D)):
            isl = ccw(r2(make_vertex_island(E4, P, Q, d, side)))
            if island_ok(*S4, E4, isl): return isl
isl19 = ccw(r2([(77.17583701182855, 9.489137146105495), (70.39874960542065, 13.371286845164954), (75.43935523515924, 19.337214676227575)]))
cases = [
    ("17a short edge at a junction", (102.0, 24.0), (99.0, 40.0), None),
    ("17b the other placement (M9)", (104.1, 52.4), (99.5, 4.9), None),
    ("18 vertex 0.5 m from a new edge", *S4, vertex_island(0.5)),
    ("19 E 0.5 m from an edge", *S4, isl19),
    ("20b 2 m to spare, junction A (M10)", *S4, None),
    ("20b(a) vertex 1.5 m from a new edge", *S4, vertex_island(1.5)),
    ("20b(b) 1.5 m at the junction", (116.5, 45.0), (84.5, 55.0), None),
    ("F12 the fold (M12), band 60", FB, FC, None, 60.0, FOLD),
]
for name, B, C, isl, *more in cases:
    band, frame = more if more else (50.0, None)
    fA, fD, frs = frame if frame else (A, D, rings(B, C, isl))
    print("=" * 100)
    pl, win = report(name, B, C, tuple(isl) if isl else None, band, frame)
    if isl: print("   island ok (inside one square, crossing neither chain, no vertex swept):", island_ok(B, C, win[1], isl))
    for lab, e, dev, q in pl:
        if (lab, e) != win[:2]:
            new = shapely.LineString([fA, e, fD]); others = [shapely.LineString([r[k], r[(k+1) % len(r)]]) for r in frs for k in range(len(r))]
            others = [s for s in others if not {tuple(s.coords[0]), tuple(s.coords[1])} & {B, C}]
            bad = [s for s in others if s.crosses(new) or (s.intersects(new) and not ({tuple(s.coords[0]), tuple(s.coords[1])} & {fA, fD}))]
            print(f"   other placement as a collapse: crosses or touches another edge: {bool(bad)}; deviation {dev:.3f} (band {band:g})")
