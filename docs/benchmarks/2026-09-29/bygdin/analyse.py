"""Increment 22 acceptance on Bygdin: the numbers in README.md, from run.sh's
logs and outputs. Usage, from the repository root with .venv active:

    python docs/benchmarks/2026-09-29/bygdin/analyse.py $SCRATCH > .../logs/analysis.txt
"""

import json
import re
import sqlite3
import statistics
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import shape

sys.path.insert(0, "tools")
import bench  # noqa: E402  quality() and read_vtk_ascii(), unchanged

from tin_engine.io.geopackage import decode_geometry  # noqa: E402

HERE = Path("docs/benchmarks/2026-09-29/bygdin")
LOG = HERE / "logs"
SCRATCH = Path(sys.argv[1])
NVE_KM2 = 305.54  # delfeltAreal_km2 of delfelt 1187
CELL = 10.0


def polygon(path: Path) -> shapely.Polygon:
    return shape(json.loads(path.read_text())["features"][0]["geometry"])


def hausdorff_1m(a: shapely.LineString, b: shapely.LineString) -> tuple[float, float]:
    """Directed distances a->b and b->a: every 1 m sample of one boundary to
    the other boundary's segments (exact point-segment distance)."""

    def directed(p: shapely.LineString, q: shapely.LineString) -> float:
        pts = shapely.points(shapely.get_coordinates(shapely.segmentize(p, 1.0)))
        c = shapely.get_coordinates(q)
        segs = shapely.linestrings(np.stack([c[:-1], c[1:]], axis=1))
        _, d = shapely.STRtree(segs).query_nearest(pts, return_distance=True, all_matches=False)
        return float(d.max())

    return directed(a, b), directed(b, a)


print("## catchment")
walls, floods, traces, reduces, rss, foot = [], [], [], [], [], []
for i in (1, 2, 3):
    err = (LOG / f"catchment_t20_run{i}.err").read_text()
    floods.append(sum(float(x) for x in re.findall(r"flood ([\d.]+) s", err)))
    traces.append(float(re.search(r"traced in ([\d.]+) s", err).group(1)))
    reduces.append(float(re.search(r"tolerance 20 m, ([\d.]+) s", err).group(1)))
    walls.append(float(re.search(r"([\d.]+) real", err).group(1)))
    rss.append(int(re.search(r"(\d+)\s+maximum resident", err).group(1)))
    foot.append(int(re.search(r"(\d+)\s+peak memory footprint", err).group(1)))
    per_window = re.findall(r"window (\d): .*?flood ([\d.]+) s", err)
    print(f"run {i}: windows {per_window}, flood sum {floods[-1]:.2f}, trace {traces[-1]},"
          f" reduce {reduces[-1]}, wall {walls[-1]}, maxRSS {rss[-1]/1e6:.0f} MB,"
          f" peak footprint {foot[-1]/1e6:.0f} MB")
med = statistics.median
print(f"median: flood {med(floods):.2f} s, trace {med(traces):.2f} s, reduce {med(reduces):.2f} s,"
      f" wall {med(walls):.2f} s, other (read, seed, write, start-up) {med(walls)-med(floods)-med(traces)-med(reduces):.2f} s,"
      f" maxRSS {med(rss)/1e6:.0f} MB, peak footprint {med(foot)/1e6:.0f} MB")

fine = polygon(SCRATCH / "catchment_t0.geojson")
nve = polygon(HERE / "nve_delfelt_1187.geojson")
print(f"\nNVE polygon: {nve.area/1e6:.4f} km2 by shapely, {shapely.get_num_coordinates(nve)} coords,"
      f" valid {nve.is_valid}; published {NVE_KM2} km2")
print("\n| tolerance | vertices | area km2 | vs fine m2 | vs NVE 305.54 | Hausdorff densify=0.05 | 1 m: fine->red, red->fine | valid |")
for t in (0, 10, 20, 50):
    red = polygon(SCRATCH / f"catchment_t{t}.geojson")
    n = shapely.get_num_coordinates(red.exterior) - 1
    h = shapely.hausdorff_distance(fine.exterior, red.exterior, densify=0.05)
    f2r, r2f = hausdorff_1m(fine.exterior, red.exterior)
    print(f"| {t} | {n} | {red.area/1e6:.6f} | {red.area-fine.area:+.3g} |"
          f" {100*(red.area/1e6-NVE_KM2)/NVE_KM2:+.3f} % | {h:.2f} m | {f2r:.2f}, {r2f:.2f} m | {red.is_valid} |")
    props = json.loads((SCRATCH / f"catchment_t{t}.geojson").read_text())["features"][0]["properties"]
    if t == 20:
        print(f"  properties: {props}")

f2n, n2f = hausdorff_1m(fine.exterior, nve.exterior)
print(f"\nfine outline vs NVE boundary, 1 m samples: fine->NVE {f2n:.0f} m, NVE->fine {n2f:.0f} m")

print("\n## node overlap with NVE (DEM nodes at multiples of 10 m)")
x0, y0, x1, y1 = shapely.union(fine, nve).bounds
xs = np.arange(np.floor(x0 / CELL) * CELL, x1 + CELL, CELL)
ys = np.arange(np.floor(y0 / CELL) * CELL, y1 + CELL, CELL)
gx, gy = np.meshgrid(xs, ys)
shapely.prepare(fine)
shapely.prepare(nve)
ours = shapely.contains_xy(fine, gx, gy)
theirs = shapely.contains_xy(nve, gx, gy)
both = ours & theirs
print(f"ours {ours.sum()} nodes (catchment reports 3049095), NVE {theirs.sum()} nodes, both {both.sum()}")
print(f"of NVE's nodes in ours: {100*both.sum()/theirs.sum():.2f} %; of ours in NVE's: {100*both.sum()/ours.sum():.2f} %")
print(f"ours only {(ours & ~theirs).sum()} nodes ({(ours & ~theirs).sum()*1e-4:.2f} km2),"
      f" NVE only {(theirs & ~ours).sum()} nodes ({(theirs & ~ours).sum()*1e-4:.2f} km2)")
diff = shapely.difference(fine, nve)
parts = sorted(getattr(diff, "geoms", [diff]), key=lambda g: -g.area)[:3]
print("largest ours-not-NVE parts:", [(round(p.area / 1e6, 3), tuple(round(v) for v in p.centroid.coords[0])) for p in parts])
diff = shapely.difference(nve, fine)
parts = sorted(getattr(diff, "geoms", [diff]), key=lambda g: -g.area)[:3]
print("largest NVE-not-ours parts:", [(round(p.area / 1e6, 3), tuple(round(v) for v in p.centroid.coords[0])) for p in parts])

print("\n## meshes (median of 3 binary runs; quality from the ASCII run)")
print("| run | triangles | vertices | start vertices | constraint edges | refine s | total s | wall s | maxRSS MB | worst angle | max degree | max error | within tol | Delaunay checked/ambiguous/violations |")
for tag, tol in (("red_t10", 10), ("fine_t10", 10), ("red_t1", 1), ("fine_t1", 1), ("feat_t10", 10)):
    rows = []
    for i in (1, 2, 3):
        md = (LOG / f"mesh_{tag}_run{i}.stats.md").read_text()
        err = (LOG / f"mesh_{tag}_run{i}.err").read_text()
        rows.append((
            float(re.search(r"^\| refine \| ([\d.]+)", md, re.M).group(1)),
            float(re.search(r"^\| \*\*total\*\* \| \*\*([\d.]+)", md, re.M).group(1)),
            float(re.search(r"([\d.]+) real", err).group(1)),
            int(re.search(r"(\d+)\s+maximum resident", err).group(1)),
        ))
    get = lambda k: med(r[k] for r in rows)  # noqa: E731
    md = (LOG / f"mesh_{tag}_ascii.stats.md").read_text()
    err = (LOG / f"mesh_{tag}_ascii.err").read_text()
    item = lambda name: int(re.search(rf"^\| {name} \| (\d+)", md, re.M).group(1))  # noqa: E731
    max_error = float(re.search(r"achieved max error ([\d.]+) m", err).group(1))
    m = bench.read_vtk_ascii(SCRATCH / f"mesh_{tag}_ascii.vtk")
    q = bench.quality(m.points, m.triangles, m.edges, tol, max_error)
    print(f"| {tag} | {item('output triangles')} | {item('output vertices')} | {item('start vertices')} |"
          f" {item('constraint edges')} | {get(0):.3f} | {get(1):.3f} | {get(2):.2f} | {get(3)/1e6:.0f} |"
          f" {q.worst_angle:.4g}° | {q.max_degree} | {max_error:.4f} | {q.within_tolerance} |"
          f" {q.delaunay_checked}/{q.delaunay_ambiguous}/{q.delaunay_violations} |")
    print(f"   runs (refine, total, wall, maxRSS): {rows}")

print("\n## land cover inside the reduced catchment (CORINE 2018, clipped by shapely)")
red = polygon(SCRATCH / "catchment_t20.geojson")
conn = sqlite3.connect("../rasputin_data/corine2018_dtm10_utm33.gpkg")
bx = red.bounds
ids = [r[0] for r in conn.execute(
    "SELECT id FROM rtree_corine2018_geom WHERE maxx >= ? AND minx <= ? AND maxy >= ? AND miny <= ?",
    (bx[0], bx[2], bx[1], bx[3]))]
area: dict[str, float] = {}
for fid in ids:
    blob, code = conn.execute("SELECT geom, code_18 FROM corine2018 WHERE fid = ?", (fid,)).fetchone()
    a = shapely.intersection(decode_geometry(blob)[1], red).area
    if a > 0:
        area[code] = area.get(code, 0.0) + a
total = sum(area.values())
print(f"{len(ids)} candidate rows; covered {total/1e6:.3f} km2 of {red.area/1e6:.3f}")
for code, a in sorted(area.items(), key=lambda kv: -kv[1]):
    print(f"| {code} | {a/1e6:.3f} km2 | {100*a/red.area:.2f} % |")
m = bench.read_vtk_ascii(SCRATCH / "mesh_feat_t10_ascii.vtk")
text = (SCRATCH / "mesh_feat_t10_ascii.vtk").read_text(encoding="latin-1").split("\n")
for name in ("land_cover", "water"):
    i = next(k for k, line in enumerate(text) if line.startswith(f"{name} 1 "))
    n = int(text[i].split()[2])
    vals = np.array(" ".join(text[i + 1 : i + 1 + n]).split()[:n], dtype=int)
    print(f"mesh cells with {name}=1: {int(vals.sum())} (of {n} cells: {len(m.edges)} lines + {len(m.triangles)} triangles)")
