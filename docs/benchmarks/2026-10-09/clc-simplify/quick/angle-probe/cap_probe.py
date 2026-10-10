# Probe: split every simplified land-cover edge longer than L into equal pieces (collinear points,
# the same on both sides of a shared edge: interpolated from the lexicographically smaller end),
# then mesh Numedalslagen and report triangles and the worst angle. L = sys.argv[1] metres (0 = off).
import sys, math, pathlib, numpy as np, shapely
from shapely.geometry import Polygon, MultiPolygon
import tin_engine.feature_input as fi
from tin_engine import cli
from tin_engine.border_simplify import BorderResult
L = float(sys.argv[1]); real = fi.simplify_borders; added = [0]
def ring(xy):
    out = []
    for p, q in zip(xy[:-1], xy[1:]):
        out.append(tuple(p)); d = math.dist(p, q)
        if L > 0 and d > L:
            n = math.ceil(d / L); lo, hi = sorted((tuple(p), tuple(q)))
            pts = [(lo[0] + (hi[0]-lo[0])*k/n, lo[1] + (hi[1]-lo[1])*k/n) for k in range(1, n)]
            out += pts if lo == tuple(p) else pts[::-1]; added[0] += n - 1
    return out
def poly(g):
    parts = [Polygon(ring(shapely.get_coordinates(q.exterior)), [ring(shapely.get_coordinates(h)) for h in q.interiors]) for q in shapely.get_parts(g)]
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)
def wrap(polys, band, *a, **kw):
    out = real(polys, band, *a, **kw)
    lengths = [math.dist(p, q) for g in out.polygons for r in shapely.get_rings(g) for p, q in zip(shapely.get_coordinates(r)[:-1], shapely.get_coordinates(r)[1:])]
    print(f"simplified edges: {len(lengths)} (both sides counted), longest {max(lengths):.0f} m, over 250 m {sum(x > 250 for x in lengths)}, over 500 m {sum(x > 500 for x in lengths)}", file=sys.stderr)
    return BorderResult(tuple(poly(g) for g in out.polygons), out.counts)
fi.simplify_borders = wrap
D = "/Users/skavhaug/projects/rasputin_data"; out = f"cap{int(L)}.vtk"
try:
    cli.app(["mesh", "--dem", f"{D}/DTM10_UTM33_20260925", "--domain", f"{D}/numedalslagen_outline_nve.geojson",
             "--features", f"{D}/corine2018_dtm10_utm33.gpkg", "--features-layer", "corine2018",
             "--features-map", "corine", "--tolerance", "10", "--ascii", "--out", out])
except SystemExit: pass
sys.path.insert(0, "/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tools")
import bench
m = bench.read_vtk_ascii(pathlib.Path(out)); P = m.points[:, :2]; T = m.triangles
def ang(p, q, r):
    u, v = q - p, r - p
    return np.degrees(np.arctan2(np.abs(u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0]), (u*v).sum(1)))
a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
mn = np.stack([ang(a, b, c), ang(b, c, a), ang(c, a, b)], 1).min(1)
print(f"L={L:g} m: points added {added[0]}, triangles {len(T)}, worst {mn.min():.3f} deg, under 0.4 {int((mn < 0.4).sum())}, under 1 {int((mn < 1).sum())}, five worst {np.round(np.sort(mn)[:5], 3).tolist()}", file=sys.stderr)
