"""THROWAWAY probe, not production code: reduce_ring on a diagonal strip of cells.

Backs the Q4 ruling text in docs/research/raster-to-vector.md. A strip of n
cells (i, i), traced as one counter-clockwise ring: lower staircase out, upper
staircase back, so each pinch vertex (k, k), k = 1..n-1, appears twice. With
eps > 0 the lower pass's pinch moves by (+eps, -eps) and the upper pass's by
(-eps, +eps), opening each pinch. Run with the repository's venv.
"""

import numpy as np

from tin_engine._core import reduce_ring


def area(r: np.ndarray) -> float:
    x, y = r[:, 0], r[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def strip(n: int, eps: float) -> np.ndarray:
    lower = [(0.0, 0.0)]
    for i in range(n):  # (i+1, i), then the pinch (i+1, i+1) or the far corner
        lower.append((i + 1.0, float(i)))
        e = eps if i + 1 < n else 0.0
        lower.append((i + 1.0 + e, i + 1.0 - e))
    upper = []
    for i in range(n - 1, -1, -1):  # (i, i+1), then the pinch (i, i) except at 0
        upper.append((float(i), i + 1.0))
        if i > 0:
            upper.append((i - eps, i + eps))
    return np.array(lower + upper)


for eps in (0.0, 1e-9):
    ring = strip(40, eps)
    print(f"eps={eps:g}: {len(ring)} vertices, area {area(ring):.6f}")
    for tol in (0.2, 0.4, 0.5, 0.6, 0.7, 1.0, 1.5, 2.0):
        out = reduce_ring(ring, tol, np.zeros((0, 2)))
        r = np.asarray(out.ring)
        print(
            f"  tol {tol:>3}: {len(r):>3} vertices, collapses {out.collapses:>3},"
            f" rejected_crossing {out.rejected_crossing:>4}, area {area(r):.6f}"
        )
