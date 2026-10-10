import sys, pickle, pathlib, numpy as np, shapely
sys.path.insert(0, "/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tools")
import bench
m = bench.read_vtk_ascii(pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "fix.vtk"))
P = m.points[:, :2]; T = m.triangles
def ang(p, q, r):
    u, v = q - p, r - p
    return np.degrees(np.arctan2(np.abs(u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0]), (u*v).sum(1)))
a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
A3 = np.stack([ang(a, b, c), ang(b, c, a), ang(c, a, b)], 1); mn = A3.min(1)
src, simp, band, kw, dom = pickle.load(open("simp.pkl", "rb"))
bnd = dom.boundary
def table(ps):
    d = {}
    for pi, p in enumerate(ps):
        for q in shapely.get_parts(p):
            for r in shapely.get_rings(q):
                for x in shapely.get_coordinates(r)[:-1]: d.setdefault(tuple(x), set()).add(pi)
    return d
S, R = table(simp), table(src)
Sk = np.array(list(S.keys()))
E = {tuple(sorted(e)) for e in m.edges.tolist()}
print("under 1 deg:", int((mn < 1).sum()), "under 0.4:", int((mn < 0.4).sum()))
for t in np.argsort(mn)[:2]:
    k = int(A3[t].argmin()); tri = T[t]; vi = tri[k]; o = [tri[(k+1) % 3], tri[(k+2) % 3]]
    print(f"{mn[t]:.4f} deg at apex ({P[vi][0]:.3f},{P[vi][1]:.3f}); arms {np.hypot(*(P[o[0]]-P[vi])):.2f} / {np.hypot(*(P[o[1]]-P[vi])):.2f}; opposite {np.hypot(*(P[o[0]]-P[o[1]])):.3f}; constrained arms {[tuple(sorted((vi, x))) in E for x in o]}, opposite {tuple(sorted(o)) in E}")
    for x in tri:
        key = tuple(P[x]); dd = np.hypot(*(Sk - P[x]).T); j = int(dd.argmin())
        kind = "no cover vertex within 1 cm"
        if dd[j] < 0.01:
            ck = tuple(Sk[j]); kind = ("source vertex" if ck in R else "simplifier E") + f", polys {sorted(S[ck])}"
        print(f"     ({P[x][0]:.3f},{P[x][1]:.3f}) line vertex {bool((m.edges == x).any())}; {kind}; {shapely.distance(shapely.Point(P[x]), bnd):.2f} m from outline")
