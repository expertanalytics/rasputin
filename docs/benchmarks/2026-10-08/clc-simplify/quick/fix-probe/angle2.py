import sys, pickle, numpy as np, shapely
sys.path.insert(0, "/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tools")
import bench, pathlib
m = bench.read_vtk_ascii(pathlib.Path("b50.vtk"))
P = m.points[:, :2]; T = m.triangles
a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
def ang(p, q, r):
    u, v = q - p, r - p
    return np.degrees(np.arctan2(np.abs(u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0]), (u*v).sum(1)))
mn = np.stack([ang(a, b, c), ang(b, c, a), ang(c, a, b)], 1).min(1)
simp = pickle.load(open("band50.pkl", "rb"))[0]; src = pickle.load(open("clip_nosimp.pkl", "rb"))[0]
def vset(ps): return {tuple(x) for x in np.concatenate([shapely.get_coordinates(p) for p in ps]).tolist()}
S, R = vset(simp), vset(src)
rings = [r for p in simp for q in shapely.get_parts(p) for r in shapely.get_rings(q)]
tree = shapely.STRtree(rings)
print("under 1 deg:", int((mn < 1).sum()), " under 0.1:", int((mn < 0.1).sum()))
for t in np.argsort(mn)[:8]:
    tri = T[t]; L = [np.hypot(*(P[tri[i]] - P[tri[(i+1)%3]])) for i in range(3)]
    i = int(np.argmin(L)); u, v = tri[i], tri[(i+1)%3]
    def kind(x):
        k = tuple(P[x]); s = "E(new, simplifier)" if k in S and k not in R else "source vertex" if k in R else "mesh/refine point"
        pt = shapely.Point(P[x]); near = tree.query(pt, predicate="dwithin", distance=0.5)
        return f"{s}; on {sum(shapely.distance(pt, rings[j]) < 1e-9 for j in near)} ring(s)"
    print(f"{mn[t]:.4f} deg: short side {L[i]*1000:.1f} mm between [{kind(u)}] and [{kind(v)}]")
