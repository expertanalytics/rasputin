"""The plant for run_sweep.py's final check: can its control fail?

Meshes the piece once (ASCII), then runs ``final_check`` twice: as is, and
with the mesh moved ``SHIFT`` metres east. The control (mesh against the grid
it was refined on) must be 0 over tolerance as is and non-zero when shifted.

    python docs/benchmarks/2026-10-01/basin-piece/plant_check.py GRID_TIF OUTLINE \
        SOURCE_WINDOW_TIF SCRATCH [TOL] [SHIFT]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import run_sweep  # noqa: E402  (sets sys.path to the Release package first)

grid, outline, window, scratch = sys.argv[1:5]
tol = float(sys.argv[5]) if len(sys.argv) > 5 else 20.0
shift = float(sys.argv[6]) if len(sys.argv) > 6 else 15.0
vtk = Path(scratch) / "plant.ascii.vtk"
run_sweep._run(["mesh", "--dem", grid, "--domain", outline, "--tolerance", f"{tol:g}", "--ascii",
                "--out", str(vtk)], Path(scratch) / "plant.log")  # fmt: skip
mesh = run_sweep.bench.read_vtk_ascii(vtk)
for s in (0.0, shift):
    r = run_sweep.final_check(mesh, Path(grid), Path(window), Path(outline), tol, shift=s)
    c = r["control_resampled_grid"]
    print(json.dumps({"tolerance": tol, "shift_m": s, "control_over_tolerance": c["over_tolerance"],
                      "control_max": c["max"], "control_outside_mesh": c["outside_mesh"],
                      "source_over_tolerance": r["source_dem"]["over_tolerance"]}))  # fmt: skip
vtk.unlink()
