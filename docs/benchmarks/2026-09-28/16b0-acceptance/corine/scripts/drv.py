"""16b-0 acceptance driver (@perf, 2026-09-28). One case per process; one JSON line out.

  drv.py --pkg PKG node   --layout A|B --side KM [--repeats N]
  drv.py --pkg PKG ladder --n RUNGS [--repeats N]      crossing grid (k = n^2 crossings)
  drv.py --pkg PKG comb   --n LINES [--repeats N]      n long parallel east-west lines, no crossing
  drv.py --pkg PKG mesh   --side KM --tol T --features none|A --min-angle DEG

Run from the repo root (clc_probe reads ../rasputin_data/corine_sql/... relative).
The pkg tree is put first on sys.path and the editable finder dropped (bench.py's child technique).
"""
import argparse, hashlib, json, resource, statistics, sys, time
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--pkg", required=True)
ap.add_argument("what")
ap.add_argument("--layout", default="A")
ap.add_argument("--side", type=float, default=48)
ap.add_argument("--n", type=int, default=300)
ap.add_argument("--repeats", type=int, default=3)
ap.add_argument("--tol", type=float, default=1.0)
ap.add_argument("--features", default="A")
ap.add_argument("--min-angle", type=float, default=25.0)
a = ap.parse_args()

sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, a.pkg)
SP = Path(__file__).resolve().parent.parent
sys.path.insert(1, str(SP))
sys.path.insert(1, str(Path(__file__).resolve().parent))  # committed copy: clc_*.py beside it
import numpy as np, shapely  # noqa: E402
from shapely.geometry import box  # noqa: E402
import tin_engine, tin_engine._core as core  # noqa: E402
from tin_engine import cli  # noqa: E402
from tin_engine._core import ChainRole  # noqa: E402

assert str(tin_engine.__file__).startswith(a.pkg) and str(core.__file__).startswith(a.pkg), (tin_engine.__file__, core.__file__)
TILE = Path("/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6603_4_10m_z33.tif")
CX, CY = 325000.0, 6675000.0  # step 0's centre: its 48 km square is 301000..349000 x 6651000..6699000
out = {"what": a.what, "pkg": a.pkg, "so_sha256": hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest()[:16]}


def square(km):
    h = km * 500.0
    return (CX - h, CY - h, CX + h, CY + h)


def digest(noded):
    v = np.ascontiguousarray(np.asarray(noded.vertices, dtype=np.float64))
    m = sorted(cli._chain_masks(noded).items())
    return hashlib.sha256(v.tobytes() + repr(m).encode()).hexdigest()[:16], len(v), len(m)


def time_node(xy, chains):
    ts = []
    for _ in range(a.repeats):
        c = cli.PhaseClock()
        run = cli._engine(xy, chains, True, cli.DEFAULT_SNAP_SPACING, c)
        ts.append(dict(c.phases())["start mesh: node"])
        assert run.ok, (run.status, run.message)
    d, nv, ne = digest(run.noded)
    out.update(node_s=ts, node_median=statistics.median(ts), noded_vertices=nv, noded_edges=ne,
               noded_digest=d, start_triangles=len(run.mesh.triangles))


if a.what == "node":
    import clc_probe, clc_mesh
    import io, contextlib
    sq = square(a.side)
    with contextlib.redirect_stdout(io.StringIO()):
        meta, clipped, rect = clc_probe.probe(TILE, sub=sq)
    dom = shapely.geometry.polygon.orient(box(*sq), 1.0)
    if a.layout == "A":  # deduplicated, clipped as lines, line-merged (clc_mesh)
        inside = clc_mesh.linework(clipped, dom, True)
        xy, chains, _, _ = clc_mesh.chains_of(inside, dom)
    else:  # B: every clipped ring as a closed Breakline (clc_noder.py)
        pts = {}; chains = []
        def idx(c):
            c = (float(c[0]), float(c[1])); pts.setdefault(c, len(pts)); return pts[c]
        chains.append(([idx(c) for c in list(dom.exterior.coords)[:-1]], ChainRole.Outer, 0))
        for oid, code, g in clipped:
            for p in getattr(g, "geoms", [g]):
                if p.geom_type != "Polygon":
                    continue
                for rr in (p.exterior, *p.interiors):
                    chains.append(([idx(c) for c in rr.coords], ChainRole.Breakline, clc_mesh.LC))
        xy = np.array(list(pts))
    segs = sum(len(ids) - 1 for ids, _, _ in chains) + 1  # +1: the Outer ring closes
    out.update(layout=a.layout, side_km=a.side, square=sq, chains=len(chains), input_vertices=len(xy),
               chain_positions=sum(len(ids) for ids, _, _ in chains), input_segments=segs)
    time_node(xy, chains)

elif a.what in ("ladder", "comb"):
    O = np.array([500000.0, 6600000.0])
    L = 700.0 * max(1.0, a.n / 300.0)  # keep rung spacing >= ~2.2 m
    ring = np.array([[0, 0], [L, 0], [L, L], [0, L]], dtype=float) + O
    pts = [ring]; chains = [([0, 1, 2, 3], ChainRole.Outer, 0)]; nid = 4
    offs = np.linspace(20.0, L - 20.0, a.n)
    for y in offs:
        pts.append(np.array([[10.0, y], [L - 10.0, y]]) + O); chains.append(([nid, nid + 1], ChainRole.Breakline, 1)); nid += 2
    if a.what == "ladder":
        for x in offs + 0.5:
            pts.append(np.array([[x, 10.0], [x, L - 10.0]]) + O); chains.append(([nid, nid + 1], ChainRole.Breakline, 2)); nid += 2
    xy = np.vstack(pts)
    out.update(n=a.n, side_m=L, input_segments=len(chains) - 1 + 4)
    time_node(xy, chains)

elif a.what == "mesh":
    import clc_probe, clc_mesh
    import io, contextlib
    from tin_engine.domain import DomainPolygon
    from tin_engine.io.repository import TiffDemRepository
    sq = square(a.side)
    dom = shapely.geometry.polygon.orient(box(*sq), 1.0)
    if a.features == "A":
        with contextlib.redirect_stdout(io.StringIO()):
            meta, clipped, rect = clc_probe.probe(TILE, sub=sq)
        inside = clc_mesh.linework(clipped, dom, True)
        xy, chains, _, _ = clc_mesh.chains_of(inside, dom)
        cli._domain_chains = lambda domain, name: (xy, chains, "probe")
        out.update(input_vertices=len(xy), chains=len(chains))
    repo = TiffDemRepository([TILE]); t = repo.load(repo.footprints()[0].name)
    d = DomainPolygon(polygon=dom, crs="EPSG:25833")
    clock = cli.PhaseClock(); t0 = time.perf_counter()
    r = cli._dem_mesh(t, "dem", None, True, cli.DEFAULT_SNAP_SPACING, a.tol, clock, domain=d, domain_name="sq",
                      min_angle=a.min_angle, feet=True)
    wall = time.perf_counter() - t0
    tr = r.trimmed
    V = np.asarray(tr.vertices); T = np.asarray(tr.triangles); P = V[T][:, :, :2]
    def ang(p, q, s):
        u = q - p; v = s - p
        return np.degrees(np.arctan2(np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]), (u * v).sum(1)))
    A = np.minimum(np.minimum(ang(P[:, 0], P[:, 1], P[:, 2]), ang(P[:, 1], P[:, 2], P[:, 0])), ang(P[:, 2], P[:, 0], P[:, 1]))
    ref = r.refinement
    out.update(side_km=a.side, tol=a.tol, features=a.features, min_angle=a.min_angle, triangles=len(T), vertices=len(V),
               start_triangles=r.start_triangles, angle_median=float(np.median(A)), under1_pct=float(np.mean(A < 1) * 100),
               worst=float(A.min()), wall_s=wall, phases={k: round(v, 4) for k, v in clock.phases()},
               refinement=repr(ref)[:600] if ref is not None else None)
else:
    raise SystemExit(a.what)

out["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
print("RESULT " + json.dumps(out), flush=True)
