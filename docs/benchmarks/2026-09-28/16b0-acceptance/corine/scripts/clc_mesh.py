"""16b probe: CORINE linework + a square domain through the real CLI engine. Scratch only."""
import sys, time, collections
from pathlib import Path
import numpy as np, shapely
from shapely.geometry import box, MultiLineString, LineString
sys.path.insert(0, str(Path(__file__).parent))
import clc_probe
from tin_engine import cli
from tin_engine.domain import DomainPolygon
from tin_engine.io.repository import TiffDemRepository
from tin_engine._core import ChainRole, node as core_node

LC = 1 << 7  # scratch "land_cover" bit; C++ bits are opaque

def linework(clipped, dom, dedupe):
    """(segments with masks) inside dom: rings as lines, clipped to the domain."""
    segs = collections.defaultdict(int); dup = []
    for oid, code, g in clipped:
        for p in getattr(g, "geoms", [g]):
            if p.geom_type != "Polygon": continue
            for rr in (p.exterior, *p.interiors):
                a = np.asarray(rr.coords)
                for i in range(len(a) - 1):
                    s = (tuple(a[i]), tuple(a[i + 1]))
                    if dedupe: segs[(min(s), max(s))] |= LC
                    else: dup.append(s)
    lines = list(segs) if dedupe else dup
    ml = MultiLineString([list(s) for s in lines])
    inside = ml.intersection(dom)
    return inside

def chains_of(inside, dom):
    pts = {}; chains = []
    def idx(c):
        c = (float(c[0]), float(c[1]))
        if c not in pts: pts[c] = len(pts)
        return pts[c]
    rings = [dom.exterior, *dom.interiors]
    for k, r in enumerate(rings):
        chains.append(([idx(c) for c in list(r.coords)[:-1]], ChainRole.Outer if k == 0 else ChainRole.Hole, 0))
    merged = shapely.line_merge(inside) if inside.geom_type != "LineString" else inside
    parts = list(getattr(merged, "geoms", [merged]))
    nclosed = 0
    for ln in parts:
        if ln.geom_type != "LineString": continue
        cs = list(ln.coords)
        ids = [idx(c) for c in cs]
        if ids[0] == ids[-1]: nclosed += 1
        chains.append((ids, ChainRole.Breakline, LC))
    xy = np.array(list(pts), dtype=np.float64)
    return xy, chains, len(parts), nclosed

def run(tile, sq, tol, dedupe=True, with_features=True, min_angle=25.0, feet=True):
    meta, clipped, rect = clc_probe.probe(Path(tile), sub=sq) if with_features else (None, [], None)
    dom = box(*sq)
    repo = TiffDemRepository([Path(tile)]); t = repo.load(repo.footprints()[0].name)
    d = DomainPolygon(polygon=shapely.geometry.polygon.orient(dom, 1.0), crs="EPSG:25833")
    if with_features:
        inside = linework(clipped, dom, dedupe)
        xy, chains, nparts, nclosed = chains_of(inside, dom)
        print(f"feature chains {nparts} ({nclosed} closed), vertices incl. domain {len(xy)}, dedupe={dedupe}")
        cli._domain_chains = lambda domain, name: (xy, chains, "probe")
    clock = cli.PhaseClock()
    t0 = time.perf_counter()
    r = cli._dem_mesh(t, "dem", None, True, cli.DEFAULT_SNAP_SPACING, tol, clock, domain=d, domain_name="sq", min_angle=min_angle, feet=feet)
    wall = time.perf_counter() - t0
    tr = r.trimmed
    V = np.asarray(tr.vertices); T = np.asarray(tr.triangles)
    P = V[T][:, :, :2]
    def ang(a, b, c):
        u = b - a; v = c - a
        return np.degrees(np.arctan2(np.abs(u[:,0]*v[:,1]-u[:,1]*v[:,0]), (u*v).sum(1)))
    A = np.minimum(np.minimum(ang(P[:,0],P[:,1],P[:,2]), ang(P[:,1],P[:,2],P[:,0])), ang(P[:,2],P[:,0],P[:,1]))
    e = np.asarray(tr.edges); m = np.asarray(tr.edge_masks)
    print(f"tol {tol}: triangles {len(T)}, vertices {len(V)}, constrained edges {len(e)} (land_cover bit on {int(np.sum(m & LC > 0))}), "
          f"min angle median {np.median(A):.1f}, <1deg {np.mean(A<1)*100:.2f} %, worst {A.min():.4f}; wall {wall:.2f} s")
    print("phases:", ", ".join(f"{n} {v:.2f}" for n, v in clock.phases() if v > 0.05))
    return tr

if __name__ == "__main__":
    tile = sys.argv[1]; sq = tuple(map(float, sys.argv[2:6])); tol = float(sys.argv[6])
    mode = sys.argv[7] if len(sys.argv) > 7 else "dedupe"
    ma = float(sys.argv[8]) if len(sys.argv) > 8 else 25.0
    run(tile, sq, tol, dedupe=(mode == "dedupe"), with_features=(mode != "none"), min_angle=ma)
