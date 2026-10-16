"""land cover, decomposed, on the captured label_triangles inputs; candidate
rewrites checked against the captured codes (array_equal)."""
import sys, time, json, pickle
import numpy as np, shapely
import tin_engine.landcover as lc
c, S = sys.argv[1], sys.argv[2]
z = np.load(f"{S}/lc_{c}.npz"); want = np.load(f"{S}/lc_{c}_codes.npy")
polys = [(shapely.from_wkb(w), code) for w, code in pickle.load(open(f"{S}/lc_{c}_polys.pkl", "rb"))]
V, T, E, margin = z["vertices"], z["triangles"], z["edges"], float(z["margin"])
def med(f, k=3):
    s = []
    for _ in range(k):
        t = time.perf_counter(); r = f(); s.append(time.perf_counter() - t)
    return round(sorted(s)[k // 2], 3), r
out = {"catchment": c}
tri = np.asarray(T, dtype=np.int64).reshape(-1, 3); cut = np.asarray(E, dtype=np.int64).reshape(-1, 2)
pv = shapely.get_num_coordinates([p for p, _ in polys])
out["sizes"] = dict(triangles=len(tri), vertices=len(V), constraint_edges=len(cut), polygons=len(polys),
                    polygon_vertices_total=int(pv.sum()), polygon_vertices_max=int(pv.max()), polygon_vertices_median=int(np.median(pv)))
out["label_triangles_s"], res = med(lambda: lc.label_triangles(V, T, E, polygons=polys, margin=margin))
assert np.array_equal(res.codes, want)
out["sizes"]["components"] = res.regions
# regions() sub-steps
n = int(max(tri.max(), cut.max())) + 1
sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]); keys = lc._keys(sides, n); ck = lc._keys(cut, n)
out["regions_s"], ids_ref = med(lambda: lc.regions(tri, cut))
out["regions: isin_s"], free = med(lambda: ~np.isin(keys, ck))
out["regions: argsort stable_s"], _ = med(lambda: np.argsort(keys[free], kind="stable"))
out["regions: argsort default_s"], _ = med(lambda: np.argsort(keys[free]))
out["sizes"]["sides"] = len(keys)
# candidate: one argsort of all keys, cut keys removed by searchsorted on the sorted array
def regions_B(tri, cut):
    n = int(max(tri.max(initial=-1), cut.max(initial=-1))) + 1
    sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    keys = lc._keys(sides, n); owner = np.tile(np.arange(len(tri), dtype=np.int64), 3)
    order = np.argsort(keys); keys, owner = keys[order], owner[order]
    ck = np.unique(lc._keys(cut, n)); pos = np.searchsorted(ck, keys); pos[pos == len(ck)] = 0
    free = ck[pos] != keys; keys, owner = keys[free], owner[free]
    pair = np.flatnonzero(keys[1:] == keys[:-1]); u, v = owner[pair], owner[pair + 1]
    parent = np.arange(len(tri), dtype=np.int64)
    while True:
        ru, rv = parent[u], parent[v]
        if np.array_equal(ru, rv): return parent
        np.minimum.at(parent, np.maximum(ru, rv), np.minimum(ru, rv))
        while not np.array_equal(parent, jumped := parent[parent]): parent = jumped
out["regions_B (one unstable argsort + searchsorted)_s"], ids_b = med(lambda: regions_B(tri, cut))
out["regions_B identical"] = bool(np.array_equal(ids_b, ids_ref))
# best-per-component selection
xy = np.asarray(V, dtype=np.float64)[:, :2]; ids = ids_ref
a, b, cc = xy[tri[:, 0]], xy[tri[:, 1]], xy[tri[:, 2]]
la, lb, lcc = (np.hypot(*(q - p).T) for p, q in ((b, cc), (cc, a), (a, b)))
per = la + lb + lcc; cross = (b[:, 0]-a[:, 0])*(cc[:, 1]-a[:, 1])-(cc[:, 0]-a[:, 0])*(b[:, 1]-a[:, 1])
safe = np.where(per > 0, per, 1.0); r = np.where(per > 0, np.abs(cross)/safe, 0.0)
centre = (la[:, None]*a + lb[:, None]*b + lcc[:, None]*cc)/safe[:, None]
def best_A():
    order = np.lexsort((np.arange(len(tri)), -r, ids)); return order[np.unique(ids[order], return_index=True)[1]]
def best_B():
    m = np.full(len(tri), -np.inf); np.maximum.at(m, ids, r)
    cand = np.flatnonzero(r == m[ids]); first = np.full(len(tri), len(tri)); np.minimum.at(first, ids[cand], cand)
    roots = np.unique(ids); return first[roots]
out["best: lexsort_s"], bA = med(best_A)
out["best_B (maximum.at/minimum.at)_s"], bB = med(best_B)
out["best_B identical"] = bool(np.array_equal(bA, bB))
pts = centre[bA]
# _lookup variants
out["_lookup_s"], (codesA, hitsA) = med(lambda: lc._lookup(pts, polys))
geoms = [p for p, _ in polys]
tree = shapely.STRtree(geoms)
bbox_pairs = tree.query(shapely.points(pts))
out["sizes"]["lookup_points"] = len(pts); out["sizes"]["bbox_candidate_pairs"] = bbox_pairs.shape[1]
out["sizes"]["vertices_in_candidate_polygons(sum over pairs)"] = int(pv[bbox_pairs[1]].sum())
def lookup_prepared():
    g2 = [shapely.from_wkb(shapely.to_wkb(g)) for g in geoms]  # fresh objects, so preparing is counted
    t = shapely.STRtree(g2); pp = shapely.points(pts)
    i, j = t.query(pp)  # bbox candidates
    shapely.prepare(np.array(g2, dtype=object)[np.unique(j)])
    ok = shapely.intersects(np.array(g2, dtype=object)[j], pp[i])
    return i[ok], j[ok]
def lookup_reverse():
    g2 = np.array([shapely.from_wkb(shapely.to_wkb(g)) for g in geoms], dtype=object)
    pp = shapely.points(pts); t = shapely.STRtree(pp)
    j, i = t.query(g2, predicate="intersects")  # polygon is the input: prepared by the query
    o = np.lexsort((j, i)); return i[o], j[o]
ref = tree.query(shapely.points(pts), predicate="intersects"); ro = np.lexsort((ref[1], ref[0])); ref = (ref[0][ro], ref[1][ro])
for name, f in {"lookup: bbox candidates then prepared intersects": lookup_prepared,
                "lookup: tree of points, query polygons": lookup_reverse}.items():
    s, (i, j) = med(f); o = np.lexsort((j, i))
    out[name + "_s"] = s; out[name + " same pairs"] = bool(np.array_equal(i[o], ref[0]) and np.array_equal(j[o], ref[1]))

def lookup_C(points, polygons):
    codes = np.zeros(len(points), dtype=np.int64); hits = np.zeros(len(points), dtype=np.int64)
    g = np.array([p for p, _ in polygons], dtype=object)
    tree = shapely.STRtree(g); pp = shapely.points(points)
    i, j = tree.query(pp); shapely.prepare(g[np.unique(j)])
    ok = shapely.intersects(g[j], pp[i]); found, which = i[ok], j[ok]
    rank = np.array([(p.area, code) for p, code in polygons], dtype=[("a", "f8"), ("c", "i8")])
    order = np.lexsort((rank["c"][which], rank["a"][which], found)); found, which = found[order], which[order]
    np.add.at(hits, found, 1); head = np.unique(found, return_index=True)[1]
    codes[found[head]] = rank["c"][which[head]]
    return codes, hits
def label_C(vertices, triangles, edges, polygons, margin):
    xy = np.asarray(vertices, dtype=np.float64)[:, :2]; tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    ids = regions_B(tri, np.asarray(edges, dtype=np.int64).reshape(-1, 2))
    a, b, c = xy[tri[:, 0]], xy[tri[:, 1]], xy[tri[:, 2]]
    la, lb, lc_ = (np.hypot(*(q - p).T) for p, q in ((b, c), (c, a), (a, b)))
    per = la + lb + lc_; cross = (b[:, 0]-a[:, 0])*(c[:, 1]-a[:, 1])-(c[:, 0]-a[:, 0])*(b[:, 1]-a[:, 1])
    safe = np.where(per > 0, per, 1.0); r = np.where(per > 0, np.abs(cross)/safe, 0.0)
    centre = (la[:, None]*a + lb[:, None]*b + lc_[:, None]*c)/safe[:, None]
    m = np.full(len(tri), -np.inf); np.maximum.at(m, ids, r)
    cand = np.flatnonzero(r == m[ids]); first = np.full(len(tri), len(tri)); np.minimum.at(first, ids[cand], cand)
    best = first[np.unique(ids)]
    cc, hits = lookup_C(centre[best], polygons)
    by_id = np.zeros(len(tri), dtype=np.int32); by_id[ids[best]] = cc
    return by_id[ids], hits, r[best]
fresh = lambda: [(shapely.from_wkb(shapely.to_wkb(p)), code) for p, code in polys]
secs = []
for _ in range(3):
    pl = fresh(); t = time.perf_counter(); codesC, hitsC, rb = label_C(V, T, E, pl, margin); secs.append(time.perf_counter() - t)
out["combined (regions_B + best_B + prepared lookup)_s"] = round(sorted(secs)[1], 3)
out["combined identical codes"] = bool(np.array_equal(codesC, want))
out["combined identical counts"] = [int((hitsC == 0).sum()), int((hitsC > 1).sum()), int((rb <= margin).sum())] == [res.outside, res.overlapped, res.thin]
print(json.dumps(out, indent=1, default=int))
