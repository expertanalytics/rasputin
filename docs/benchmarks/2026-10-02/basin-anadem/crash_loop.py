"""Crash-rate loop: ``rasputin mesh --dem anadem-v1 --out-crs`` on the Velhas
piece, N process starts per CRS, every exit status recorded.

A measurement script, not production code; nothing imports it.

    python crash_loop.py RASPUTIN CACHE OUTLINE OUT_JSON SCRATCH N CRS [CRS ...]

RASPUTIN is the command that runs the CLI (the venv's ``rasputin``). Each run
is ``mesh --dem anadem-v1 --cache CACHE --domain OUTLINE --out-crs CRS
--tolerance 10 --binary --out SCRATCH/x.vtk``. The runs alternate between the
CRSs. A run fails if it exits non-zero or on a signal; its stderr tail is kept.
New macOS crash reports (``~/Library/Logs/DiagnosticReports/python*.ips``
written after the loop started) are listed at the end.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path


def main() -> None:
    rasputin, cache, outline, out_json, scratch = sys.argv[1:6]
    n, crss = int(sys.argv[6]), sys.argv[7:]
    Path(scratch).mkdir(parents=True, exist_ok=True)
    reports = Path.home() / "Library/Logs/DiagnosticReports"
    t_start = time.time()
    runs: list[dict] = []
    for i in range(n):
        for k, crs in enumerate(crss):
            cmd = [rasputin, "mesh", "--dem", "anadem-v1", "--cache", cache, "--domain", outline,
                   "--out-crs", crs, "--tolerance", "10", "--binary",
                   "--out", str(Path(scratch) / f"x{k}.vtk")]  # fmt: skip
            t0 = time.time()
            p = subprocess.run(cmd, capture_output=True, text=True)
            rec = {"i": i, "crs_index": k, "returncode": p.returncode,
                   "wall_s": round(time.time() - t0, 3)}  # fmt: skip
            if p.returncode != 0:
                rec["stderr_tail"] = p.stderr[-2000:]
            runs.append(rec)
        bad = [r for r in runs if r["returncode"] != 0]
        print(f"{len(runs)} runs, {len(bad)} failed", flush=True)
        Path(out_json).write_text(json.dumps({"crss": crss, "runs": runs}, indent=1))
    new = sorted(p.name for p in reports.glob("python*.ips") if p.stat().st_mtime >= t_start)
    result = {"crss": crss, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t_start)),
              "runs": runs, "failed": sum(r["returncode"] != 0 for r in runs),
              "new_crash_reports": new}  # fmt: skip
    Path(out_json).write_text(json.dumps(result, indent=1))
    print(f"done: {len(runs)} runs, {result['failed']} failed, crash reports {new}")


if __name__ == "__main__":
    main()
