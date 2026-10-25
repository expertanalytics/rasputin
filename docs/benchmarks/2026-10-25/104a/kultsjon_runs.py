"""Kultsjön timing pairs for 104a's acceptance (@perf, 2026-10-25).

Each run is one child process of the venv's Python with a tree's Release
package (`<tree>/build-bench/pkg`, made by `tools/bench.py run`) first on
`sys.path` and the editable finder dropped, as bench.py's child does.
Wall time is taken in this parent. `pmset -g batt` is read before each run.

    python kultsjon_runs.py --head TREE --base TREE --out DIR --repeats 5
        [--extra-tolerances 1e6,53.99,...]

Writes DIR/runs.jsonl (one line per run) and keeps each run's --stats and
--record in DIR/<label>-r<i>.{md,json}; meshes go to /tmp/perf-104a-meshes/
(kept out of git; rerun this script to regenerate them).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

HOME = Path.home()
CATCH = HOME / "rasputin-output/opal-se2-2026-10-21/catchments/kultsjon.geojson"
FEATURES = HOME / ".cache/rasputin/clc2018/features.gpkg"
COMMON = [
    "--dem", "glo30", "--domain", str(CATCH), "--domain-crs", "EPSG:3006",
    "--out-crs", "EPSG:3006", "--features", str(FEATURES), "--features-legend", "corine",
    "--features-repair", "1.0", "--features-tolerance", "50.0",
    "--features-outline-snap", "5.0", "--threads", "4", "--memory-budget", "4G",
    "--z-unit", "m",
]  # fmt: skip
PRELUDE = (
    "import sys; pkg = sys.argv.pop(1);"
    "sys.meta_path[:] = [f for f in sys.meta_path if 'ScikitBuild' not in type(f).__name__];"
    "sys.path.insert(0, pkg); import tin_engine, tin_engine._core;"
    "assert tin_engine.__file__.startswith(pkg) and tin_engine._core.__file__.startswith(pkg);"
    "sys.argv[0] = 'rasputin'; from tin_engine.cli import app; app()"
)


def power() -> str:
    return subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout.strip()


MESHES = Path("/tmp/perf-104a-meshes")


def run(label: str, tree: Path, args: list[str], out: Path, suffix: str, rep: int) -> dict:
    meshes = MESHES
    meshes.mkdir(parents=True, exist_ok=True)
    stem = f"{label}-r{rep}"
    argv = [sys.executable, "-c", PRELUDE, str(tree / "build-bench/pkg"), "mesh",
            "--out", str(meshes / f"{stem}{suffix}"), *COMMON, *args,
            "--stats", str(out / f"{stem}.md"), "--record", str(out / f"{stem}.json")]  # fmt: skip
    before = power()
    t0 = time.perf_counter()
    proc = subprocess.run(argv, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    row = {"label": label, "rep": rep, "tree": str(tree), "args": args, "suffix": suffix,
           "wall_s": wall, "exit": proc.returncode, "stderr": proc.stderr[-2000:],
           "power_before": before, "power_after": power()}  # fmt: skip
    with open(out / "runs.jsonl", "a") as f:
        f.write(json.dumps(row) + "\n")
    print(f"{label} r{rep} exit {proc.returncode} {wall:.3f} s", flush=True)
    return row


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--head", type=Path, required=True)
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--extra-tolerances", default="")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    pairs = [
        ("base-t10", a.base, ["--tolerance", "10.0"]),
        ("head-t10", a.head, ["--tolerance", "10.0"]),
        ("head-n65893", a.head, ["--triangles", "65893"]),
    ]
    for rep in range(1, a.repeats + 1):  # interleaved, so a drift spreads over all three
        for label, tree, args in pairs:
            run(label, tree, args, a.out, ".gpkg", rep)
    for label, tree, args in pairs:  # one .vtk each, for byte comparison
        run(label, tree, args, a.out, ".vtk", 0)
    for t in filter(None, a.extra_tolerances.split(",")):
        run(f"head-t{t}", a.head, ["--tolerance", t], a.out, ".vtk", 0)


if __name__ == "__main__":
    main()
