"""Refine's thread scaling on the basin piece: the serial share the basin will meet.

A measurement script, not production code; nothing imports it. The refine
call is forced to each thread count by ``tools/bench.py _child`` (its
``--threads``), and samples are interleaved round-robin over thread counts
so a drift spreads over all of them, as ``bench.py run`` does.

    python docs/benchmarks/2026-10-02/basin-piece-anadem/scaling.py GRID_TIF OUTLINE OUT_DIR SCRATCH \
        [--tolerance 1] [--threads 1,2,4,6,8,10] [--repeats 3]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import run_sweep  # noqa: E402

grid, outline, out, scratch = sys.argv[1:5]
rest = sys.argv[5:]
opt = lambda k, d: rest[rest.index(k) + 1] if k in rest else d  # noqa: E731
tol = float(opt("--tolerance", "1"))
threads = [int(t) for t in opt("--threads", "1,2,4,6,8,10").split(",")]
repeats = int(opt("--repeats", "3"))
out_dir = Path(out)
(out_dir / "logs").mkdir(parents=True, exist_ok=True)
rec = {"tolerance": tol, "grid": grid, "outline": outline, "pmset_before": run_sweep._pmset(),
       "samples": []}  # fmt: skip
for r in range(repeats):
    for t in threads:
        tag = f"scaling_t{tol:g}_threads{t}_run{r + 1}"
        s = run_sweep._run(["mesh", "--dem", grid, "--domain", outline, "--tolerance", f"{tol:g}",
                            "--binary", "--out", str(Path(scratch) / "scaling.bin.vtk")],
                           out_dir / "logs" / f"{tag}.log", threads=t)  # fmt: skip
        s.update(threads=t, run=r + 1)
        rec["samples"].append(s)
        print(f"{tag}: refine {s['refine_s']:.3f} s, process {s['proc_s']:.2f} s", flush=True)
rec["pmset_after"] = run_sweep._pmset()
rec["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
(out_dir / "scaling.json").write_text(json.dumps(rec, indent=1))
one = np.median([s["refine_s"] for s in rec["samples"] if s["threads"] == 1])
print("\n| threads | refine (median) | speed-up | process wall (median) |\n|---:|---:|---:|---:|")
for t in threads:
    m = np.median([s["refine_s"] for s in rec["samples"] if s["threads"] == t])
    p = np.median([s["proc_s"] for s in rec["samples"] if s["threads"] == t])
    print(f"| {t} | {m:.2f} s | {one / m:.2f} | {p:.2f} s |")
