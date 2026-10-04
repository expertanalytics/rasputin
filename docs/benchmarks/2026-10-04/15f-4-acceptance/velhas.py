"""15f-4 acceptance: the Velhas piece on ANADEM. Copied from 15f-3-acceptance/velhas.py;
PY points at 15f-3's venv, as 15f-4's worktree has none.

15c-2's run_geo.py, rewritten for the current CLI: its VTK-header regex no
longer matches, because the "checked against N source nodes" sentence is gone
from both trees. Phase 2's figures are read from --stats instead.
A measurement script, not production code; nothing imports it.

    python velhas.py run LABEL TREE OUT SCRATCH [--tolerances 20,10,5] [--repeats 2]
    python velhas.py check OUT [--workers 3] [--shift 30]

`run`: per tolerance, --repeats timed --binary runs through TREE's
tools/bench.py child (TREE/build-bench/pkg, the editable finder dropped, CLI
default threads), each under /usr/bin/time -l, with --stats kept. Then one
--ascii run with --stats, and bench.quality() on it. pmset before and after
each tolerance. Everything goes to OUT/results.json.
`check`: for each ascii mesh, 15c-2's independent source-node check
(run_geo.independent_check: matplotlib's interpolator at every valid ANADEM
node inside the domain, interior and strip), and with --shift the control.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
W = HERE.parents[3]
DATA = Path("/Users/skavhaug/projects/rasputin_data")
SF = DATA / "sao_francisco_piece"
OUTLINE = SF / "bho2017_5k_76949_outline_epsg4674.geojson"
WINDOW = SF / "bho2017_5k_76949_anadem_window_epsg4674.tif"
CRS = "EPSG:31983"
PY = "/Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/.venv/bin/python"  # 15f-4 has no venv


def opt(rest: list[str], k: str, d: str) -> str:
    return rest[rest.index(k) + 1] if k in rest else d


def pmset() -> str:
    return subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout.strip()


def swap() -> str:
    return subprocess.run(["sysctl", "-n", "vm.swapusage"], capture_output=True, text=True).stdout.strip()


def rows(stats_md: Path) -> dict[str, str]:
    """Every two-column row of a --stats file, and the total."""
    text = stats_md.read_text()
    out = dict(re.findall(r"^\| ([a-z][^|]*?) \| ([^|]+?) \|$", text, re.M))
    tot = re.search(r"^\| \*\*total\*\* \| \*\*([\d.]+)\*\*", text, re.M)
    if tot:
        out["total_s"] = tot.group(1)
    return out


def one(tree: Path, args: list[str], log: Path) -> dict:
    argv = ["/usr/bin/time", "-l", PY, str(tree / "tools/bench.py"), "_child", "--pkg",
            str(tree / "build-bench/pkg"), "--threads", "0", "--", *args]  # fmt: skip
    t0 = time.perf_counter()
    p = subprocess.run(argv, capture_output=True, text=True, cwd=tree)
    wall = time.perf_counter() - t0
    log.write_text(p.stdout + "\n--- stderr\n" + p.stderr)
    if p.returncode != 0:
        raise SystemExit(f"{log}: exit {p.returncode}")
    child = json.loads(re.search(r"^BENCH (\{.*\})$", p.stderr, re.M).group(1)) if "BENCH " in p.stderr else {}
    rss = int(re.search(r"(\d+)\s+maximum resident set size", p.stderr).group(1))
    real = float(re.search(r"([\d.]+) real", p.stderr).group(1))
    return {"wall_s": wall, "real_s": real, "max_rss_bytes": rss, "child": child}


def run(label: str, tree: Path, out: Path, scratch: Path, rest: list[str]) -> None:
    sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
    sys.path[:0] = [str(W / "build-bench/pkg"), str(W / "tools")]
    import bench

    tols = [float(t) for t in opt(rest, "--tolerances", "20,10,5").split(",")]
    repeats = int(opt(rest, "--repeats", "2"))
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=tree).stdout.strip()
    res: dict = {"label": label, "tree": str(tree), "commit": commit, "crs": CRS, "runs": []}
    for tol in tols:
        tag = f"t{tol:g}"
        base = ["mesh", "--dem", "anadem-v1", "--cache", str(DATA / "cache"), "--domain", str(OUTLINE),
                "--out-crs", CRS, "--tolerance", f"{tol:g}"]  # fmt: skip
        rec: dict = {"tolerance": tol, "pmset_before": pmset(), "swap_before": swap(), "timed": []}
        rec["started_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        for r in range(repeats):
            st = logs / f"{tag}_run{r + 1}.stats.md"
            t = one(tree, base + ["--binary", "--out", str(scratch / f"{tag}.bin.vtk"), "--stats", str(st)],
                    logs / f"{tag}_run{r + 1}.log")  # fmt: skip
            t["stats"] = rows(st)
            rec["timed"].append(t)
        vtk = scratch / f"{tag}.ascii.vtk"
        st = logs / f"{tag}_ascii.stats.md"
        q = one(tree, base + ["--ascii", "--out", str(vtk), "--stats", str(st)], logs / f"{tag}_ascii.log")
        q["stats"] = rows(st)
        rec["ascii"] = q
        mesh = bench.read_vtk_ascii(vtk)
        max_err = float(q["stats"].get("max_error_m", "nan").split()[0])
        rec["quality"] = bench.quality(mesh.points, mesh.triangles, mesh.edges, tol, max_err).model_dump()
        rec["mesh_sha256"] = mesh.sha256
        rec["vtk"] = str(vtk)
        rec["pmset_after"], rec["swap_after"] = pmset(), swap()
        res["runs"].append(rec)
        (out / "results.json").write_text(json.dumps(res, indent=1))
        print(f"{label} {tag}: real {[round(x['real_s'], 2) for x in rec['timed']]} s, "
              f"rss {max(x['max_rss_bytes'] for x in rec['timed']) / 2**30:.2f} GiB", flush=True)


def check(out: Path, rest: list[str]) -> None:
    sys.path[:0] = [str(W / "docs/benchmarks/2026-10-02/15c-2-acceptance")]
    import run_geo  # sets its own paths; independent_check uses matplotlib, not tin_engine

    shift = float(opt(rest, "--shift", "0"))
    path = out / "results.json"
    res = json.loads(path.read_text())
    for rec in res["runs"]:
        h = 30.0
        rec["check"] = run_geo.independent_check(Path(rec["vtk"]), WINDOW, OUTLINE, CRS, h, rec["tolerance"])
        if shift:
            rec["control"] = run_geo.independent_check(Path(rec["vtk"]), WINDOW, OUTLINE, CRS, h,
                                                       rec["tolerance"], shift)  # fmt: skip
        path.write_text(json.dumps(res, indent=1))
        c = rec["check"]
        print(f"{res['label']} t{rec['tolerance']:g}: interior over {c['interior']['over_tolerance']} of "
              f"{c['interior']['nodes']}, strip over {c['strip']['over_tolerance']} of {c['strip']['nodes']}"
              + (f"; control +{shift:g} m: {rec['control']['interior']['over_tolerance']}" if shift else ""),
              flush=True)  # fmt: skip


if __name__ == "__main__":
    if sys.argv[1] == "run":
        a = sys.argv[2:]
        run(a[0], Path(a[1]), Path(a[2]), Path(a[3]), a[4:])
    else:
        check(Path(sys.argv[2]), sys.argv[3:])
