"""Tolerance sweep of ``rasputin mesh --domain`` over one resampled basin piece.

A measurement script, not production code (ROADMAP basin item 2.1); nothing
imports it. Run from the repository root with the venv's python, after
``tools/bench.py``'s Release build exists in ``build-bench/pkg``:

    python docs/benchmarks/2026-10-01/basin-piece/run_sweep.py PIECE GRID_TIF OUTLINE \
        SOURCE_WINDOW_TIF OUT_DIR SCRATCH [--tolerances 1,2,5,10,20,50] [--repeats 3]

Per tolerance:
- ``--repeats`` timed runs, ``--binary``, through ``tools/bench.py _child``
  (the Release package, the CLI's default thread count, refine timed inside
  the process), each under ``/usr/bin/time -l`` for the max RSS, the parent
  timing the whole process;
- one ``--ascii`` run for quality: ``tools/bench.py``'s ``quality()`` (worst
  angle, max degree, the tolerance check, the constrained Delaunay check);
- the **final check** of Q6's ruling, measured, not built: the mesh, read
  linearly at every source-DEM node inside the domain (projected), against
  the source value. It counts the source nodes that a final check would find
  outside the tolerance. The **control** reads the mesh at every node of the
  grid it was refined against, inside the domain, where refine guarantees
  the tolerance: a non-zero control count means this measurement is wrong;
- ``pmset -g batt`` before and after.

Results: ``OUT_DIR/results.json`` and the per-run logs and ``--stats``
files under ``OUT_DIR/logs``. Meshes go to SCRATCH.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
PKG = ROOT / "build-bench" / "pkg"
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path[:0] = [str(PKG), str(ROOT / "tools"), str(Path(__file__).parent)]

import bench  # noqa: E402  tools/bench.py: quality() and read_vtk_ascii()
import matplotlib.tri as mtri  # noqa: E402
import prep_dem  # noqa: E402
import pyproj  # noqa: E402
import shapely  # noqa: E402
import shapely.ops  # noqa: E402
import tifffile  # noqa: E402

NODATA = prep_dem.NODATA


def _run(args: list[str], log: Path, threads: int = 0) -> dict:
    t0 = time.perf_counter()
    done = subprocess.run(["/usr/bin/time", "-l", sys.executable, str(ROOT / "tools/bench.py"),
                           "_child", "--pkg", str(PKG), "--threads", str(threads), "--", *args],
                          capture_output=True, text=True, cwd=ROOT)  # fmt: skip
    proc_s = time.perf_counter() - t0
    log.write_text(f"$ rasputin {' '.join(args)}\n--- stdout\n{done.stdout}\n--- stderr\n{done.stderr}")
    if done.returncode != 0:
        raise SystemExit(f"run failed ({done.returncode}), see {log}")
    rec = json.loads(re.search(r"^BENCH (.*)$", done.stderr, re.M).group(1))
    rec["proc_s"] = proc_s
    rec["max_rss_bytes"] = int(re.search(r"(\d+)\s+maximum resident set size", done.stderr).group(1))
    return rec


def _sizes(stats_md: Path) -> dict:
    text = stats_md.read_text()
    get = lambda item: int(re.search(rf"^\| {item} \| (\d+)", text, re.M).group(1))  # noqa: E731
    total = float(re.search(r"^\| \*\*total\*\* \| \*\*([\d.]+)\*\*", text, re.M).group(1))
    threads = re.search(r"Threads: (\d+)", text).group(1)
    return {"vertices": get("output vertices"), "triangles": get("output triangles"),
            "start_vertices": get("start vertices"), "dropped": get("vertices without data dropped"),
            "command_total_s": total, "threads": int(threads)}  # fmt: skip


def _errors(tri: mtri.Triangulation, z: np.ndarray, px: np.ndarray, py: np.ndarray,
            ref: np.ndarray, tol: float) -> dict:  # fmt: skip
    interp = mtri.LinearTriInterpolator(tri, z)
    errs, outside = [], 0
    for i in range(0, px.size, 2_000_000):
        v = interp(px[i : i + 2_000_000], py[i : i + 2_000_000])
        m = np.ma.getmaskarray(v)
        outside += int(m.sum())
        errs.append(np.abs(np.asarray(v[~m]) - ref[i : i + 2_000_000][~m]))
    if px.size == 0:
        raise SystemExit("final check: no nodes to check; the frames disagree")
    e = np.concatenate(errs)
    over = e > tol + 1e-6  # refine's own comparison carries float noise of this order
    return {"nodes": int(e.size), "outside_mesh": outside, "over_tolerance": int(over.sum()),
            "over_share": float(over.mean()), "max": float(e.max()),
            "p99": float(np.quantile(e, 0.99)), "median": float(np.median(e))}  # fmt: skip


def final_check(mesh: bench.VtkMesh, grid_tif: Path, window_tif: Path, outline: Path,
                tol: float, shift: float = 0.0) -> dict:  # fmt: skip
    """``shift`` moves the mesh east by that many metres: the plant that the
    control must catch."""
    pts = mesh.points
    tri = mtri.Triangulation(pts[:, 0] + shift, pts[:, 1], mesh.triangles)
    poly, crs = prep_dem._outline(outline)
    grid, gx0, gy0, h, _ = prep_dem.read_geotiff(grid_tif)
    with tifffile.TiffFile(grid_tif) as tif:  # the grid's own CRS, never assumed
        epsg = int(tif.geotiff_metadata["ProjectedCSTypeGeoKey"])
    to_grid = pyproj.Transformer.from_crs(crs, f"EPSG:{epsg}", always_xy=True)
    poly_g = shapely.ops.transform(to_grid.transform, poly)
    rr, cc = np.mgrid[0 : grid.shape[0], 0 : grid.shape[1]]
    gx, gy = gx0 + cc * h, gy0 - rr * h
    inside = shapely.contains_xy(poly_g, gx, gy) & (grid != NODATA)
    control = _errors(tri, pts[:, 2], gx[inside], gy[inside], grid[inside].astype(np.float64), tol)
    src, lonf, latf, d, _ = prep_dem.read_geotiff(window_tif)
    rr, cc = np.mgrid[0 : src.shape[0], 0 : src.shape[1]]
    valid = src != NODATA
    fwd = pyproj.Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    px, py = fwd.transform(lonf + cc[valid] * d, latf - rr[valid] * d)
    sz = src[valid].astype(np.float64)
    inside = shapely.contains_xy(poly_g, px, py)  # in the mesh's own frame, as the mesh was cut
    px, py, sz = px[inside], py[inside], sz[inside]
    # The boundary strip: within one grid spacing of the domain's edge. Refine
    # checks grid nodes only (increment 16: "at every valid node inside the
    # domain"), so a sliver along an edge holding no grid node is never checked;
    # the final check would check source nodes there. Reported apart.
    strip = shapely.contains_xy(poly_g.boundary.buffer(h), px, py)
    return {"control_resampled_grid": control,
            "source_dem": _errors(tri, pts[:, 2], px, py, sz, tol),
            "source_dem_interior": _errors(tri, pts[:, 2], px[~strip], py[~strip], sz[~strip], tol),
            "source_dem_strip": _errors(tri, pts[:, 2], px[strip], py[strip], sz[strip], tol),
            "strip_width_m": h}  # fmt: skip


def _pmset() -> str:
    return subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout


def main() -> None:
    piece, grid_tif, outline, window_tif, out, scratch = sys.argv[1:7]
    rest = sys.argv[7:]
    tols = [float(t) for t in (rest[rest.index("--tolerances") + 1] if "--tolerances" in rest
                               else "1,2,5,10,20,50").split(",")]  # fmt: skip
    repeats = int(rest[rest.index("--repeats") + 1]) if "--repeats" in rest else 3
    out_dir, scr = Path(out), Path(scratch)
    logs = out_dir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    scr.mkdir(parents=True, exist_ok=True)
    results = {"piece": piece, "grid": grid_tif, "outline": outline, "source_window": window_tif,
               "package": str(PKG), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "runs": []}  # fmt: skip
    for tol in tols:
        tag = f"{piece}_t{tol:g}"
        base = ["mesh", "--dem", grid_tif, "--domain", outline, "--tolerance", f"{tol:g}"]
        rec = {"tolerance": tol, "pmset_before": _pmset(), "timed": []}
        for r in range(repeats):
            rec["timed"].append(_run(base + ["--binary", "--out", str(scr / f"{tag}.bin.vtk"),
                                             "--stats", str(logs / f"{tag}_run{r + 1}.stats.md")],
                                            logs / f"{tag}_run{r + 1}.log"))  # fmt: skip
            rec["timed"][-1].update(_sizes(logs / f"{tag}_run{r + 1}.stats.md"))
        rec["pmset_after"] = _pmset()
        vtk = scr / f"{tag}.ascii.vtk"
        q = _run(base + ["--ascii", "--out", str(vtk), "--stats", str(logs / f"{tag}_ascii.stats.md")],
                 logs / f"{tag}_ascii.log")  # fmt: skip
        mesh = bench.read_vtk_ascii(vtk)
        rec["quality"] = bench.quality(mesh.points, mesh.triangles, mesh.edges, tol,
                                       q["max_error"]).model_dump()  # fmt: skip
        rec["quality"]["mesh_sha256"] = mesh.sha256
        rec["final_check"] = final_check(mesh, Path(grid_tif), Path(window_tif), Path(outline), tol)
        results["runs"].append(rec)
        med = lambda k: float(np.median([t[k] for t in rec["timed"]]))  # noqa: E731
        print(f"tol {tol:g}: {rec['timed'][0]['triangles']} triangles, {rec['timed'][0]['vertices']} "
              f"vertices, proc {med('proc_s'):.2f} s, refine {med('refine_s'):.2f} s; "
              f"source nodes over tol {rec['final_check']['source_dem']['over_tolerance']}, "
              f"control over {rec['final_check']['control_resampled_grid']['over_tolerance']}",
              flush=True)  # fmt: skip
        (out_dir / "results.json").write_text(json.dumps(results, indent=1))
        vtk.unlink()  # the ASCII mesh is large; the binary one is kept in SCRATCH


if __name__ == "__main__":
    main()
