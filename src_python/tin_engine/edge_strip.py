"""The edge strip (increment 15f, D5 and D6): check points where constraints
cross grid lines, and between them, after ``refine``.

Two calls into ``_core`` and their ``--stats`` rows; no geometry. ``generate``
runs on both paths while the tile is held; ``run`` is the projected path's
strip run (``refine_strip``). On the reprojected path the strip goes to
``final_check.run`` instead, which joins it to the source nodes in one loop.
"""

from __future__ import annotations

import numpy as np

from tin_engine._core import (
    ConstraintCheckPoints,
    PointRefineOutcome,
    RasterView,
    RefineOutcome,
    constraint_check_points,
    refine_strip,
)
from tin_engine.stats import PhaseClock


def generate(view: RasterView, start: RefineOutcome, clock: PhaseClock) -> ConstraintCheckPoints:
    """The strip on ``start``'s constraint edges, z from ``view``."""
    with clock.phase("edge strip: generate"):
        return constraint_check_points(view, np.asarray(start.vertices), np.asarray(start.edges))


def run(
    view: RasterView,
    strip: ConstraintCheckPoints,
    start: RefineOutcome,
    tolerance: float,
    clock: PhaseClock,
    feet: bool = False,
) -> PointRefineOutcome:
    """``refine_strip`` from ``start``; ``clock`` gets the run's own two times.
    ``feet`` puts points near a line onto it first (20c, R4)."""
    arrays = (start.vertices, start.triangles, start.z, start.valid, start.edges, start.masks)
    out = refine_strip(
        view, strip, *(np.asarray(a) for a in arrays), tolerance=tolerance, constraint_feet=feet
    )
    clock.add("edge strip: scan (parallel)", out.scan_seconds)
    clock.add("edge strip: split + flip (serial)", out.split_seconds)
    return out
