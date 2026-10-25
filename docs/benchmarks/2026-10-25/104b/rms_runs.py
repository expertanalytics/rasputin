"""The error summary's cost for 104b's acceptance (@perf, 2026-10-25).

Runs the quick check's cases (docs/benchmarks/quick/cases.toml) from one tree's
Release package (`<tree>/build-bench/pkg`, made by `tools/bench.py run`),
each with and without `--rms`, interleaved, `--stats` and `--record` kept.
The child is the venv's Python with the editable finder dropped, as
bench.py's child does. `pmset -g batt` is read before and after each run.

    python rms_runs.py --tree TREE --out DIR [--repeats 3] [--cases tile,velhas]
        [--sweep tile:1,2,4,8]

`--sweep CASE:T1,T2,...` adds `--rms --threads T` runs of one case, for the
pass's thread scaling. Writes DIR/runs.jsonl (one line per run: wall, the
`--stats` total and error summary row, the mesh's sha256, power); the
`--stats` and `--record` files go to DIR/stats/. Meshes go to
/tmp/perf-104b-meshes/ (not kept; rerun to regenerate).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import tomllib
from pathlib import Path
from string import Template

HERE = Path(__file__).resolve().parent
MESHES = Path("/tmp/perf-104b-meshes")
PRELUDE = (
    "import sys; pkg = sys.argv.pop(1);"
    "sys.meta_path[:] = [f for f in sys.meta_path if 'ScikitBuild' not in type(f).__name__];"
    "sys.path.insert(0, pkg); import tin_engine, tin_engine._core;"
    "assert tin_engine.__file__.startswith(pkg) and tin_engine._core.__file__.startswith(pkg);"
    "sys.argv[0] = 'rasputin'; from tin_engine.cli import app; app()"
)
ROW = re.compile(r"^\| (\*\*)?(?P<phase>[^|*]+?)(\*\*)? \| (\*\*)?(?P<s>[0-9.]+)(\*\*)? \|", re.M)


def power() -> str:
    return subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout.strip()


def phases(md: Path) -> dict[str, float]:
    text = md.read_text()
    timings = text[text.index("## Timings") :]
    return {m["phase"]: float(m["s"]) for m in ROW.finditer(timings)}


def run(tree: Path, label: str, args: list[str], out: Path, rep: int) -> None:
    MESHES.mkdir(parents=True, exist_ok=True)
    (out / "stats").mkdir(parents=True, exist_ok=True)
    stem = f"{label}-r{rep}"
    mesh = MESHES / f"{stem}.vtk"
    md, rec = out / "stats" / f"{stem}.md", out / "stats" / f"{stem}.json"
    argv = [sys.executable, "-c", PRELUDE, str(tree / "build-bench/pkg"), *args,
            "--out", str(mesh), "--binary", "--stats", str(md), "--record", str(rec)]  # fmt: skip
    before = power()
    t0 = time.perf_counter()
    proc = subprocess.run(argv, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    row: dict[str, object] = {"label": label, "rep": rep, "args": args, "wall_s": wall,
                              "exit": proc.returncode, "stderr": proc.stderr[-1500:],
                              "power_before": before, "power_after": power()}  # fmt: skip
    if proc.returncode == 0:
        ph = phases(md)
        record = json.loads(rec.read_text())
        row |= {"total_s": ph.get("total"), "error_summary_s": ph.get("error summary"),
                "mesh_sha256": hashlib.sha256(mesh.read_bytes()).hexdigest(),
                "record": {k: record.get(k) for k in ("triangles", "max_error_m", "rms_error_m",
                           "error_nodes_measured", "tolerance_m")}}  # fmt: skip
        mesh.unlink()
    with open(out / "runs.jsonl", "a") as f:
        f.write(json.dumps(row) + "\n")
    print(f"{label} r{rep} exit {proc.returncode} total {row.get('total_s')} "
          f"summary {row.get('error_summary_s')}", flush=True)  # fmt: skip


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tree", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--cases", default="")
    p.add_argument("--sweep", default="")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    data = os.environ.get("RASPUTIN_DATA", str(a.tree.resolve().parents[3] / "rasputin_data"))
    cases = tomllib.loads((a.tree / "docs/benchmarks/quick/cases.toml").read_text())["case"]
    keep = set(filter(None, a.cases.split(",")))
    plan: list[tuple[str, list[str]]] = []
    for case in cases:
        if keep and case["name"] not in keep:
            continue
        args = [Template(x).substitute(TREE=str(a.tree), RASPUTIN_DATA=data) for x in case["args"]]
        for t in case["threads"]:
            th = [] if t == 0 else ["--threads", str(t)]
            name = case["name"] + ("" if t == 0 else f"-t{t}")
            plan += [(f"{name}-plain", args + th), (f"{name}-rms", args + th + ["--rms"])]
    sweeps: list[tuple[str, list[str]]] = []
    if a.sweep:
        name, counts = a.sweep.split(":")
        case = next(c for c in cases if c["name"] == name)
        args = [Template(x).substitute(TREE=str(a.tree), RASPUTIN_DATA=data) for x in case["args"]]
        sweeps = [(f"{name}-rms-t{t}", [*args, "--threads", t, "--rms"]) for t in counts.split(",")]
    for rep in range(1, a.repeats + 1):  # interleaved, so a drift spreads over every label
        for label, args in plan + sweeps:
            run(a.tree, label, args, a.out, rep)


if __name__ == "__main__":
    main()
