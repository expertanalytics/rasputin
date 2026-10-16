# how close does each simplifier-made vertex (E) come to a vertex or edge of another border, against the source's own closest approach
import pickle, numpy as np, shapely
simp = pickle.load(open("band50.pkl", "rb"))[0]; src = pickle.load(open("clip_nosimp.pkl", "rb"))[0]
def verts(ps): return np.unique(np.concatenate([shapely.get_coordinates(p) for p in ps]), axis=0)
def edges(ps):
    out = []
    for p in ps:
        for q in shapely.get_parts(p):
            for r in shapely.get_rings(q):
                xy = shapely.get_coordinates(r); out.append(np.stack([xy[:-1], xy[1:]], 1))
    e = np.concatenate(out); e = np.sort(e.view([('x', float), ('y', float)]).reshape(-1, 2), axis=1).view(float).reshape(-1, 2, 2)
    return np.unique(e, axis=0)
for name, ps in (("source (clipped, repaired)", src), ("simplified", simp)):
    V = verts(ps); E = edges(ps); L = shapely.linestrings(E); tree = shapely.STRtree(L)
    pts = shapely.points(V)
    i, j = tree.query(pts, predicate="dwithin", distance=1.0)
    d = shapely.distance(pts[i], L[j])
    keep = d > 1e-6  # the vertex is not on this edge (not its end, not lying on it)
    dmin = np.full(len(V), np.inf); np.minimum.at(dmin, i[keep], d[keep])
    print(f"{name}: {len(V)} vertices; vertex-to-other-edge clearance under 1 cm: {int((dmin < 0.01).sum())}, under 10 cm: {int((dmin < 0.1).sum())}, under 1 m: {int((dmin < 1).sum())}; smallest {dmin.min()*1000:.2f} mm")
