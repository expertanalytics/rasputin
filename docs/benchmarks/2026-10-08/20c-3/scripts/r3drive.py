"""R3 (d) and (e): wrap snap_to_outline, run `rasputin mesh` in-process, measure
its polygons as covfar.py does and its lines near M2's population-3 lines.
usage: r3drive.py TAG OUT_CRS -- <mesh args>"""
import sys, numpy as np, shapely, pyproj
from shapely.ops import transform
import tin_engine.feature_input as fi
from tin_engine.cli import app
tag, crs = sys.argv[1], sys.argv[2]; args = sys.argv[sys.argv.index("--") + 1:]
calls = []
orig = fi.snap_to_outline
def wrap(polygons, outline, distance):
    r = orig(polygons, outline, distance); calls.append((r, outline)); return r
fi.snap_to_outline = wrap
try: app(["mesh", *args], standalone_mode=False)
except SystemExit: pass
def covfar(G, D):
    G = np.array([g for g in G if not g.is_empty], dtype=object)
    ins = shapely.STRtree(G).query(D, predicate="intersects")
    E = shapely.coverage_invalid_edges(G[ins])
    return sum(shapely.difference(shapely.intersection(e, D), D.boundary.buffer(1e-3)).length for e in E if e is not None and not e.is_empty)
far = sum(covfar(r.polygons, D) for r, D in calls)
print(f"{tag}: calls {len(calls)}; (d) mismatched border inside, beyond 1 mm of the outline: {far:.3f} m")
# plant: move one vertex shared by two polygons, far inside, 0.5 m in one polygon only
r, D = calls[0]; G = list(r.polygons); tree = shapely.STRtree(G); inner = D.buffer(-50); planted = None
for i, g in enumerate(G):
    if g.is_empty or g.geom_type != "Polygon": continue
    c = shapely.get_coordinates(g.exterior)
    for k in range(1, len(c) - 1):
        p = shapely.Point(c[k])
        if not inner.contains(p): continue
        nb = [j for j in tree.query(p, predicate="intersects") if j != i]
        if nb:
            c2 = c.copy(); c2[k] += (0.5, 0.0)
            if k == 0: c2[-1] = c2[0]
            G[i] = shapely.Polygon(c2, [h.coords for h in g.interiors]); planted = (i, k, nb[0]); break
    if planted: break
print(f"{tag}: plant {planted}; (d) with one shared vertex moved 0.5 m: {covfar(G, D):.3f} m")
# (e): lines in EPSG:3035, segments with both ends within 1 m of N 3811923.31 or E 4585680.19, inside the domain
to = pyproj.Transformer.from_crs(crs, "EPSG:3035", always_xy=True).transform
tot = 0.0
for r, D in calls:
    D35 = transform(to, D)
    for ls in r.lines:
        for line in ls:
            c = shapely.get_coordinates(transform(to, line))
            a, b = c[:-1], c[1:]
            on = ((np.abs(a[:, 1] - 3811923.31) <= 1) & (np.abs(b[:, 1] - 3811923.31) <= 1)) | ((np.abs(a[:, 0] - 4585680.19) <= 1) & (np.abs(b[:, 0] - 4585680.19) <= 1))
            if on.any():
                segs = shapely.linestrings(np.stack([a[on], b[on]], 1))
                tot += float(shapely.length(shapely.intersection(segs, D35)).sum())
print(f"{tag}: (e) population-3 length inside the domain: {tot:.3f} m")
