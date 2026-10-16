import pickle, numpy as np, shapely
exec(open("clear.py").read().split("for name, ps in")[0])
polysR, outline = pickle.load(open("optionR.pkl", "rb"))
base = pickle.load(open("band0.pkl", "rb"))[0]
import tin_engine.feature_input as fi
ruled = [fi._polygonal(shapely.intersection(p, outline)) for p in fi.snap_to_outline(base, outline, 5.0).polygons]
ruled = [p for p in ruled if p is not None]
bnd = outline.boundary
for name, ps in (("after rule and clip, before simplifying", ruled), ("option R, simplified", polysR)):
    V = verts(ps); E = edges(ps); L = shapely.linestrings(E); tree = shapely.STRtree(L); pts = shapely.points(V)
    i, j = tree.query(pts, predicate="dwithin", distance=5.0); dd = shapely.distance(pts[i], L[j]); k = dd > 1e-6
    dmin = np.full(len(V), np.inf); np.minimum.at(dmin, i[k], dd[k])
    db = shapely.distance(pts, bnd); off = db > 1e-6
    print(f"{name}: {len(V)} vertices; clearance to a non-incident edge under 1 m: {int((dmin < 1).sum())}, under 5 m: {int((dmin < 5).sum())}, smallest {dmin.min():.4f} m; vertices off the outline but within 5 m of it: {int((off & (db < 5)).sum())}")
