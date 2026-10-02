"""ANADEM against GLO-30 on the piece's 30 m grid: registration and roughness.

A measurement script, not production code; nothing imports it. Both grids
were made by ``prep_dem.py resample`` on the same target nodes. If the ANADEM
window were mis-registered (for instance read as node-registered when it is
area-registered, half a cell off), the difference would shrink when one grid
is shifted by a cell; it must be smallest at shift (0, 0).

    python docs/benchmarks/2026-10-02/basin-piece-anadem/align_check.py ANADEM_GRID GLO30_GRID
"""

from __future__ import annotations

import json
import sys

import numpy as np
import tifffile

a = tifffile.imread(sys.argv[1]).astype(np.float64)
g = tifffile.imread(sys.argv[2]).astype(np.float64)
assert a.shape == g.shape, (a.shape, g.shape)
out: dict = {"shift_rows_cols": {}}
for dr in (-1, 0, 1):
    for dc in (-1, 0, 1):
        A = a[1 + dr : a.shape[0] - 1 + dr, 1 + dc : a.shape[1] - 1 + dc]
        G = g[1:-1, 1:-1]
        m = (A != -9999) & (G != -9999)
        d = A[m] - G[m]
        out["shift_rows_cols"][f"{dr},{dc}"] = {"mean": round(float(d.mean()), 3),
                                                "median_abs": round(float(np.median(np.abs(d))), 3),
                                                "rms": round(float(np.sqrt((d**2).mean())), 3)}  # fmt: skip
best = min(out["shift_rows_cols"], key=lambda k: out["shift_rows_cols"][k]["rms"])
out["best_shift"] = best
for name, z in (("anadem", a), ("glo30", g)):
    lap = np.abs(z[1:-1, 2:] + z[1:-1, :-2] + z[2:, 1:-1] + z[:-2, 1:-1] - 4 * z[1:-1, 1:-1])
    out[f"{name}_mean_abs_laplacian_m"] = round(float(lap.mean()), 3)
    out[f"{name}_p99_abs_laplacian_m"] = round(float(np.quantile(lap, 0.99)), 3)
print(json.dumps(out, indent=1))
if best != "0,0":
    raise SystemExit(f"best fit at shift {best}, not 0,0: mis-registered")
