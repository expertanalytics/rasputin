"""Tables for 15f-3-acceptance.md part (2), Bygdin (@perf). python summarize_bygdin.py SCRATCH > tables-bygdin.md
SCRATCH holds bygdin.sh's meshes (<side>_t<tol>_ascii.vtk). Timings: the four --binary runs a side
and tolerance (blocks B N N B); /usr/bin/time -l real and peak RSS; --stats rows and phases."""
import json
import re
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
W = HERE.parents[3]
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path[:0] = [str(W / "build-bench/pkg"), str(W / "tools")]
import bench  # noqa: E402

L = HERE / "bygdin"
S = Path(sys.argv[1])
ROWS = ["line_points_checked", "line_points_on_nodata", "line_points_inserted", "line_check_dem_nodes_inserted",
        "line_points_refused", "line_points_refused_max_error_m", "line_max_error_m", "line_points_duplicate",
        "max_error_m", "refinement_rounds", "points_inserted", "points_snapped_to_lines"]
PHASES = ["refine", "edge strip: generate", "edge strip: scan (parallel)", "edge strip: split + flip (serial)",
          "trim", "land cover", "write: encode", "other"]


def stats(p):
    t = p.read_text()
    named = {m[2]: m[1].strip() for m in re.finditer(r"^\|[^|]*\| ([^|]*) \| `([a-z_0-9]+)` \|$", t, re.M)}
    ph = {m[1]: float(m[2]) for m in re.finditer(r"^\| ([^|]+?) \| \*?\*?([\d.]+)\*?\*? \| \*?\*?[\d.]+ %", t, re.M)}
    sizes = dict(re.findall(r"^\| (output vertices|output triangles|constraint edges) \| (\d+) \|$", t, re.M))
    return named, ph, sizes


def timed(err):
    t = err.read_text()
    return float(re.search(r"([\d.]+) real", t)[1]), int(re.search(r"(\d+)\s+maximum resident", t)[1])


for tol in ("10", "1"):
    print(f"### Bygdin, tolerance {tol} m\n")
    cols = {}
    for side in ("base", "15f3"):
        runs = sorted(L.glob(f"{side}_t{tol}_b*_run*.err"))
        reals, rss = zip(*(timed(r) for r in runs))
        named, ph, sizes = stats(L / f"{side}_t{tol}_ascii.stats.md")
        phm = {k: st.median([stats(Path(str(r)[:-4] + ".stats.md"))[1].get(k, 0.0) for r in runs]) for k in PHASES}
        m = bench.read_vtk_ascii(S / f"{side}_t{tol}_ascii.vtk")
        q = bench.quality(m.points, m.triangles, m.edges, float(tol), float(named["max_error_m"]))
        under10 = None
        cols[side] = dict(real=st.median(reals), real_range=f"{min(reals):.2f}-{max(reals):.2f}", rss=max(rss) / 2**30,
                          named=named, ph=phm, sizes=sizes, q=q, n=len(runs), sha=m.sha256)
    print("| item | base c193cb1 | 15f-3 83c7fd2 |\n|---|---:|---:|")
    b, n = cols["base"], cols["15f3"]
    print(f"| process wall, median of {b['n']} / {n['n']} (s) | {b['real']:.2f} ({b['real_range']}) | {n['real']:.2f} ({n['real_range']}) |")
    print(f"| peak RSS, max (GiB) | {b['rss']:.2f} | {n['rss']:.2f} |")
    for k in ("output vertices", "output triangles", "constraint edges"):
        print(f"| {k} | {b['sizes'][k]} | {n['sizes'][k]} |")
    for k in ROWS:
        print(f"| `{k}` | {b['named'].get(k, '-')} | {n['named'].get(k, '-')} |")
    for k in PHASES:
        print(f"| phase: {k} (s, median) | {b['ph'][k]:.3f} | {n['ph'][k]:.3f} |")
    for k in ("worst_angle", "angle_median", "share_under_1", "max_degree", "within_tolerance", "delaunay_violations", "delaunay_checked"):
        print(f"| quality: {k} | {getattr(b['q'], k)} | {getattr(n['q'], k)} |")
    print(f"| mesh sha256 | `{b['sha'][:16]}` | `{n['sha'][:16]}` |")
    print()
