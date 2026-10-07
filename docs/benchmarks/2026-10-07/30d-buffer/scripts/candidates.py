"""The candidate fixes in isolation, and how GEOS buffer time scales on the
staircase outlines. Each timing is the median of three runs (all three kept).

Usage: PKG=<dir> python candidates.py <lagan|ljungan_flasjo|numedalslagen> <pieces|region|scaling>

pieces   GEOS's mitred buffer of the domain against the outline grown in pieces
         of 250, 500 and 1,000 edges (launch.py's piecewise_buffer, the same
         code the patched runs use). Differences: the area of the symmetric
         difference; the largest distance from either boundary to the other,
         sampled at every vertex and every 0.25 m along each edge
         (shapely.segmentize), so a gap narrower than that between samples
         could be missed; how far each pokes out of the other (the same
         sampling, only the points the other does not contain);
         the largest difference of the four bounds; whether
         every vertex of each lies within 1e-6 m of the other's boundary.
region   feature_input.source_region as master computes it (convex hull of the
         domain grown by 100 m, round) against the hull grown by 100 m.
scaling  mitred buffer time against the distance, on the whole polygon; and
         against the vertex count, on open pieces of the outline (flat ends)
         at the pipeline's distance.
"""

import json
import math
import os
import statistics
import sys
import time
from pathlib import Path

sys.meta_path[:] = [f for f in sys.meta_path if "skbc" not in type(f).__module__]
sys.path.insert(0, str(Path(os.environ["PKG"]).resolve()))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import shapely  # noqa: E402

from tin_engine import feature_input as fi  # noqa: E402
from tin_engine.io.domain_file import read_domain  # noqa: E402

DATA = Path("/Users/skavhaug/projects/rasputin_data")
c, what = sys.argv[1], sys.argv[2]
if c == "numedalslagen":
    path = Path("/Users/skavhaug/projects/rasputin_scratch/norway/numedalslagen/numedalslagen_outline_nve.geojson")
    crs, d = "EPSG:25833", math.sqrt(2) * 10  # dem_input.py:250, the DTM10 cell diagonal
else:
    name = "lagan_mouth" if c == "lagan" else c
    path = DATA / f"sweden_smhi_svar/{name}_svar2022_3006.geojson"
    crs, d = "EPSG:3006", math.sqrt(2) * 31  # target_grid.py:107, dem_input.py:213
domain = read_domain(path, None).to_crs(crs)
P = domain.polygon
n = len(P.exterior.coords) - 1


def timed(f, repeat=3):  # type: ignore[no-untyped-def]
    runs, r = [], None
    for _ in range(repeat):
        t0 = time.perf_counter()
        r = f()
        runs.append(time.perf_counter() - t0)
    return r, runs


def emit(**row):  # type: ignore[no-untyped-def]
    print(json.dumps({"catchment": c, "test": what, "edges": n, **row}), flush=True)


def to_boundary(pts, b):  # type: ignore[no-untyped-def]
    """Distance from each point to b's boundary, through an STRtree of its edges."""
    segs = []
    for ring in shapely.get_rings(b) if b.geom_type == "Polygon" else [r for p in b.geoms for r in shapely.get_rings(p)]:
        xy = shapely.get_coordinates(ring)
        segs.append(shapely.linestrings(np.stack([xy[:-1], xy[1:]], axis=1)))
    tree = shapely.STRtree(np.concatenate(segs))
    _, dist = tree.query_nearest(shapely.points(pts), return_distance=True, all_matches=False)
    return dist


def far(a, b, step=0.25):  # type: ignore[no-untyped-def]
    """Largest distance from points sampled on a's boundary to b's boundary."""
    return float(to_boundary(shapely.get_coordinates(shapely.segmentize(a.boundary, step)), b).max())


def outside(a, b, step=0.25):  # type: ignore[no-untyped-def]
    """How far a pokes out of b: the largest distance to b's boundary from the
    points sampled on a's boundary that b does not contain (0 if none)."""
    pts = shapely.get_coordinates(shapely.segmentize(a.boundary, step))
    shapely.prepare(b)
    out = pts[~shapely.contains_xy(b, pts[:, 0], pts[:, 1])]
    return float(to_boundary(out, b).max()) if len(out) else 0.0


def vertices_within(a, b, tol=1e-6):  # type: ignore[no-untyped-def]
    return bool((to_boundary(shapely.get_coordinates(a.boundary), b) <= tol).all())


if what == "pieces":
    from launch_pieces import piecewise_buffer

    exact, runs = timed(lambda: P.buffer(d, join_style="mitre"))
    emit(method="GEOS buffer", distance=d, runs=runs, median=statistics.median(runs),
         vertices_out=int(shapely.get_num_coordinates(exact)), area_m2=exact.area)  # fmt: skip
    for k in (250, 500, 1000):
        got, runs = timed(lambda k=k: piecewise_buffer(P, d, k))
        bounds = float(np.max(np.abs(np.subtract(got.bounds, exact.bounds))))
        emit(method=f"pieces of {k}", distance=d, runs=runs, median=statistics.median(runs),
             vertices_out=int(shapely.get_num_coordinates(got)), geom_type=got.geom_type,
             holes=len(getattr(got, "interiors", [])), exact_holes=len(exact.interiors),
             sym_diff_m2=got.symmetric_difference(exact).area,
             far_got_to_exact_m=far(got, exact), far_exact_to_got_m=far(exact, got),
             got_outside_exact_m=outside(got, exact), exact_outside_got_m=outside(exact, got),
             bounds_max_diff_m=bounds, bounds_equal=got.bounds == exact.bounds,
             vertices_within_1e6=vertices_within(got, exact) and vertices_within(exact, got),
             equals_exact_0=got.equals_exact(exact, 0.0) if got.geom_type == exact.geom_type else False)  # fmt: skip
elif what == "region":
    master, runs = timed(lambda: fi.source_region(domain, crs, crs))
    emit(method="master: hull of the grown domain", runs=runs, median=statistics.median(runs), area_m2=master.area)

    def hull_first():  # type: ignore[no-untyped-def]
        hull = shapely.convex_hull(P)
        ring = shapely.segmentize(hull.buffer(fi.MARGIN).exterior, fi.DENSIFY)
        return shapely.convex_hull(shapely.multipoints(shapely.get_coordinates(ring)))

    cand, runs = timed(hull_first)
    emit(method="candidate: the hull grown", runs=runs, median=statistics.median(runs), area_m2=cand.area,
         sym_diff_m2=cand.symmetric_difference(master).area,
         far_cand_to_master_m=far(cand, master), far_master_to_cand_m=far(master, cand),
         master_outside_cand_m=outside(master, cand),
         bounds_max_diff_m=float(np.max(np.abs(np.subtract(cand.bounds, master.bounds)))),
         cand_covers_domain_grown_99_5=bool(cand.covers(P.buffer(fi.MARGIN - 0.5))))  # fmt: skip
elif what == "scaling":
    for dist in (5.0, 10.0, 20.0, 31.0, d, 62.0, 100.0):
        _, runs = timed(lambda dist=dist: P.buffer(dist, join_style="mitre"))
        emit(shape="polygon", join="mitre", distance=dist, vertices=n, runs=runs, median=statistics.median(runs))
    _, runs = timed(lambda: P.buffer(fi.MARGIN))
    emit(shape="polygon", join="round", distance=fi.MARGIN, vertices=n, runs=runs, median=statistics.median(runs))
    xy = shapely.get_coordinates(P.exterior)
    for m in (250, 500, 1000, 2000, 4000, 8000, 16000, 32000, n):
        if m > n:
            continue
        line = shapely.linestrings(xy[: m + 1])
        _, runs = timed(lambda line=line: line.buffer(d, join_style="mitre", cap_style="flat"))
        emit(shape="open piece", join="mitre", distance=d, vertices=m + 1, runs=runs, median=statistics.median(runs))
