"""Kultsjön checks beside the timing pairs (@perf, 2026-10-25): one run each.

- the budget mesh against `--tolerance <tolerance_m found>` (design test A4 on Kultsjön);
- each budget trial re-run as `--tolerance T`, for its largest error (the
  error-against-count points), with T rebuilt exactly from the floor's error
  (E_F / 4, / 16, / 64) and the found tolerance from the record;
- P1's 2 % cap (`--max-triangles 67210`) and P4 (`--triangles 30000`, `27000`).

    python kultsjon_points.py --head TREE --out DIR
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from kultsjon_runs import run


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--head", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    found = json.loads((a.out / "head-n65893-r0.json").read_text())["tolerance_m"]
    run("head-tfound", a.head, ["--tolerance", repr(found)], a.out, ".vtk", 0)
    run("head-tfloor", a.head, ["--tolerance", "1e6"], a.out, ".vtk", 0)
    ef = json.loads((a.out / "head-tfloor-r0.json").read_text())["max_error_m"]
    for k in (4, 16, 64):
        run(f"head-tEF-over-{k}", a.head, ["--tolerance", repr(ef / k)], a.out, ".vtk", 0)
    run("head-n65893-cap67210", a.head, ["--triangles", "65893", "--max-triangles", "67210"], a.out, ".vtk", 0)
    run("head-n30000", a.head, ["--triangles", "30000"], a.out, ".vtk", 0)
    run("head-n27000", a.head, ["--triangles", "27000"], a.out, ".vtk", 0)


if __name__ == "__main__":
    main()
