"""15c-2 acceptance, part (b): ``rasputin mesh --dem anadem-v1 --out-crs`` on the
Velhas piece, with the basin-piece's independent final check on each mesh.

A measurement script, not production code; nothing imports it. Run from the
repository root with the venv's python, after ``tools/bench.py``'s Release
build exists in ``build-bench/pkg``:

    python docs/benchmarks/2026-10-02/15c-2-acceptance/run_geo.py OUTLINE WINDOW_TIF CACHE \
        OUT_DIR SCRATCH --crs 'EPSG:31983;+proj=...' --tolerances 50,20,10,5 [--repeats 3] \
        [--no-check] [--shift 15]
    python .../run_geo.py check OUT_DIR/results.json OUTLINE WINDOW_TIF [WORKERS [MIN_TOL]]

Per (CRS, tolerance):
- ``--repeats`` timed ``--binary`` runs through ``tools/bench.py _child`` (the
  Release package, the CLI's default thread count), each under
  ``/usr/bin/time -l`` for the max RSS; the ``--stats`` file of each is kept;
- one ``--ascii`` run: ``bench.quality()`` (worst angle, max degree, the
  tolerance check, the constrained Delaunay check) and the independent final
  check: ``basin-piece-anadem/run_sweep.py``'s ``_errors`` (matplotlib's
  linear interpolator, a different reader and locator than rasputin's) at
  every valid ANADEM node of the source window inside the domain, projected
  into the mesh's CRS by pyproj, split into the interior and the strip within
  one grid spacing of the domain's edge. With ``--shift`` the mesh is moved
  east by that many metres first: the control that must fail.
- ``pmset -g batt`` before and after.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
PKG = ROOT / "build-bench" / "pkg"
PIECE = ROOT / "docs/benchmarks/2026-10-02/basin-piece-anadem"
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path[:0] = [str(PKG), str(ROOT / "tools"), str(PIECE)]

import bench  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import prep_dem  # noqa: E402
import pyproj  # noqa: E402
import shapely  # noqa: E402
import shapely.ops  # noqa: E402
import tifffile  # noqa: E402
from run_sweep import _errors, _pmset  # noqa: E402
from run_sweep import _run as _run_once  # noqa: E402

FAILURES: list[dict] = []


def _run(args: list[str], log: Path, tries: int = 5) -> dict:
    """``run_sweep._run``, retried: 15c-2 at d6ac78d fails some runs at random
    (a GEOSException or SIGBUS in the check points' block pool). Every failed
    attempt is recorded in FAILURES with its log kept, never dropped."""
    for k in range(1, tries + 1):
        try:
            return _run_once(args, log)
        except SystemExit:
            kept = log.with_name(f"{log.stem}.fail{k}.log")
            log.rename(kept)
            text = kept.read_text(errors="replace")
            kind = ("GEOSException" if "GEOSException" in text else
                    "terminated abnormally (signal)" if "terminated abnormally" in text else "other")
            FAILURES.append({"log": kept.name, "attempt": k, "kind": kind,
                             "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})  # fmt: skip
    raise SystemExit(f"{tries} attempts failed, see {log}")


def _stats(stats_md: Path) -> dict:
    text = stats_md.read_text()
    get = lambda item: int(re.search(rf"^\| {item} \| (\d+)", text, re.M).group(1))  # noqa: E731
    rows = dict(re.findall(r"^\| ([a-z][^|]*?) \| ([\d.]+) \|", text, re.M))
    total = float(re.search(r"^\| \*\*total\*\* \| \*\*([\d.]+)\*\*", text, re.M).group(1))
    return {"vertices": get("output vertices"), "triangles": get("output triangles"),
            "threads": int(re.search(r"Threads: (\d+)", text).group(1)),
            "command_total_s": total, "rows": rows}  # fmt: skip


def independent_check(vtk: Path, window: Path, outline: Path, crs: str, h: float,
                      tol: float, shift: float = 0.0) -> dict:  # fmt: skip
    mesh = bench.read_vtk_ascii(vtk)
    pts = mesh.points
    tri = mtri.Triangulation(pts[:, 0] + shift, pts[:, 1], mesh.triangles)
    poly, dom_crs = prep_dem._outline(outline)
    to_mesh = pyproj.Transformer.from_crs(dom_crs, crs, always_xy=True)
    poly_m = shapely.ops.transform(to_mesh.transform, poly)
    src, lonf, latf, d, _ = prep_dem.read_geotiff(window)
    with tifffile.TiffFile(window) as tif:  # the window's own CRS, never assumed
        src_epsg = int(tif.geotiff_metadata["GeographicTypeGeoKey"])
    rr, cc = np.mgrid[0 : src.shape[0], 0 : src.shape[1]]
    valid = src != prep_dem.NODATA
    fwd = pyproj.Transformer.from_crs(f"EPSG:{src_epsg}", crs, always_xy=True)
    px, py = fwd.transform(lonf + cc[valid] * d, latf - rr[valid] * d)
    sz = src[valid].astype(np.float64)
    inside = shapely.contains_xy(poly_m, px, py)
    px, py, sz = px[inside], py[inside], sz[inside]
    strip = shapely.contains_xy(poly_m.boundary.buffer(h), px, py)
    return {"shift_m": shift, "strip_width_m": h,
            "interior": _errors(tri, pts[:, 2], px[~strip], py[~strip], sz[~strip], tol),
            "strip": _errors(tri, pts[:, 2], px[strip], py[strip], sz[strip], tol),
            "mesh_sha256": mesh.sha256,
            "quality": None}  # fmt: skip


def main() -> None:
    outline, window, cache, out, scratch = (Path(a) for a in sys.argv[1:6])
    rest = sys.argv[6:]
    opt = lambda k, d: rest[rest.index(k) + 1] if k in rest else d  # noqa: E731
    crss = opt("--crs", "EPSG:31983").split(";")  # a PROJ string may hold commas
    tols = [float(t) for t in opt("--tolerances", "50,20,10,5").split(",")]
    repeats, shift = int(opt("--repeats", "3")), float(opt("--shift", "0"))
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    path = out / "results.json"
    results = json.loads(path.read_text()) if path.exists() else {"runs": []}
    results.update({"outline": str(outline), "window": str(window), "package": str(PKG)})
    for crs in crss:
        for tol in tols:
            name = "tmerc" if crs.startswith("+proj=tmerc") else crs.replace(":", "")
            tag = f"{name}_t{tol:g}"
            base = ["mesh", "--dem", "anadem-v1", "--cache", str(cache), "--domain", str(outline),
                    "--out-crs", crs, "--tolerance", f"{tol:g}"]  # fmt: skip
            rec: dict = {"crs": crs, "tolerance": tol, "pmset_before": _pmset(), "timed": []}
            rec["started_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            for r in range(repeats):
                st = logs / f"{tag}_run{r + 1}.stats.md"
                t = _run(base + ["--binary", "--out", str(scratch / f"{tag}.bin.vtk"),
                                 "--stats", str(st)], logs / f"{tag}_run{r + 1}.log")  # fmt: skip
                t.update(_stats(st))
                rec["timed"].append(t)
            rec["pmset_after"] = _pmset()
            vtk = scratch / f"{tag}.ascii.vtk"
            q = _run(base + ["--ascii", "--out", str(vtk)], logs / f"{tag}_ascii.log")
            head = urllib.parse.unquote(vtk.read_text(errors="replace")[:20000])  # %20 in fields
            h = float(re.search(r"square ([\d.]+) m grid", head).group(1))
            m = re.search(r"checked against (\d+) source nodes: (\d+) inserted in (\d+) rounds, "
                          r"max error ([^ ]+) m", head)  # fmt: skip
            rec["phase2"] = {"source_nodes": int(m.group(1)), "inserted": int(m.group(2)),
                             "rounds": int(m.group(3)), "max_error": float(m.group(4))}  # fmt: skip
            mesh = bench.read_vtk_ascii(vtk)
            rec["quality"] = bench.quality(mesh.points, mesh.triangles, mesh.edges, tol,
                                           q["max_error"]).model_dump()  # fmt: skip
            rec["vtk"], rec["h"] = str(vtk), h
            if "--no-check" not in rest:
                rec["check"] = independent_check(vtk, window, outline, crs, h, tol)
            if shift:
                rec["control"] = independent_check(vtk, window, outline, crs, h, tol, shift)
            results["runs"].append(rec)
            results["failures"] = FAILURES
            path.write_text(json.dumps(results, indent=1))
            med = lambda k: float(np.median([x[k] for x in rec["timed"]]))  # noqa: E731
            if "check" not in rec:
                continue
            c = rec["check"]
            print(f"{crs} tol {tol:g}: {rec['timed'][0]['triangles']} tri, proc {med('proc_s'):.2f} s, "
                  f"rss {med('max_rss_bytes') / 2**30:.2f} GiB, phase2 +{rec['phase2']['inserted']}; "
                  f"over: interior {c['interior']['over_tolerance']}, strip {c['strip']['over_tolerance']}"
                  + (f"; control(+{shift:g} m) {rec['control']['interior']['over_tolerance']}"
                     if shift else ""), flush=True)  # fmt: skip


def _one(job: tuple) -> dict:
    return independent_check(*job)


def check_all(results_json: Path, outline: Path, window: Path, workers: int,
              min_tol: float = 0.0) -> None:  # fmt: skip
    """The independent check for every run in ``results_json`` that lacks
    one, ``workers`` processes at a time (after the timed runs, so it never
    shares the machine with them)."""
    from concurrent.futures import ProcessPoolExecutor

    results = json.loads(results_json.read_text())
    todo = [r for r in results["runs"] if "check" not in r and r["tolerance"] >= min_tol]
    jobs = [(Path(r["vtk"]), window, outline, r["crs"], r["h"], r["tolerance"]) for r in todo]
    with ProcessPoolExecutor(workers) as pool:
        for rec, c in zip(todo, pool.map(_one, jobs), strict=True):
            rec["check"] = c
            results_json.write_text(json.dumps(results, indent=1))
            print(f"{rec['crs'][:20]} tol {rec['tolerance']:g}: interior over "
                  f"{c['interior']['over_tolerance']} of {c['interior']['nodes']}, strip over "
                  f"{c['strip']['over_tolerance']} of {c['strip']['nodes']}", flush=True)  # fmt: skip


if __name__ == "__main__":
    if sys.argv[1] == "check":  # check RESULTS_JSON OUTLINE WINDOW [WORKERS [MIN_TOL]]
        a = sys.argv[2:]
        check_all(Path(a[0]), Path(a[1]), Path(a[2]), int(a[3]) if len(a) > 3 else 4,
                  float(a[4]) if len(a) > 4 else 0.0)  # fmt: skip
    else:
        main()
