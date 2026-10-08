"""The final check against the source DEM (increment 15c-2, D5 and D6).

Phase 1 (`refine` on the resampled grid) is the caller's, unchanged. Phase 2
is here: the source's nodes, already in the target CRS, filed in a
`CheckPoints` on the target grid's frame, then `refine_points` from phase 1's
mesh until every one is within the tolerance (J2). Numbers only cross into
`_core`: the grid's frame in metres, the points, the mesh.
"""

from __future__ import annotations

import time
from collections.abc import Iterable

import numpy as np

from tin_engine._core import (
    CheckPoints,
    ConstraintCheckPoints,
    LineTolerance,
    PointRefineOutcome,
    RefineOutcome,
    refine_points,
)
from tin_engine.stats import PhaseClock
from tin_engine.target_grid import Block, TargetGrid


def run(
    start: RefineOutcome,
    grid: TargetGrid,
    checks: Iterable[Block],
    tolerance: float,
    clock: PhaseClock,
    strip: ConstraintCheckPoints | None = None,
    feet: bool = False,
    field: LineTolerance | None = None,
) -> tuple[PointRefineOutcome, int]:
    """Phase 2 from phase 1's outcome `start`; the outcome and the number of
    check points stored. `clock` gets D7's rows. `strip`, the edge strip on
    the target grid (15f, D6), joins the same loop. `feet` puts points near
    a line onto it first (20c, R4). `field`, if any, sets each triangle's
    tolerance (increment 33)."""
    h = float(grid.spacing)
    store = CheckPoints(
        x_min=grid.col0 * h, y_max=-grid.row0 * h, spacing=h, rows=grid.rows, cols=grid.cols
    )
    t0, adding = time.perf_counter(), 0.0
    for xy, z in checks:
        t1 = time.perf_counter()
        store.add(xy, z)
        adding += time.perf_counter() - t1
    clock.add("check points: project", time.perf_counter() - t0 - adding)
    with clock.phase("check points: store"):
        store.freeze()
    arrays = (start.vertices, start.triangles, start.z, start.valid, start.edges, start.masks)
    out = refine_points(
        store,
        *(np.asarray(a) for a in arrays),
        tolerance=tolerance,
        strip=strip,
        constraint_feet=feet,
        field=field,
    )
    clock.add("check points: store", adding)
    clock.add("final check: scan (parallel)", out.scan_seconds)
    clock.add("final check: split + flip (serial)", out.split_seconds)
    return out, store.size
