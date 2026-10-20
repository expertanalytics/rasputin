"""The edge strip: check points where constraints cross grid lines, and
between them, checked after ``refine``, since a thin triangle along a
constraint can hold no DEM node, and ``refine`` checks a triangle only at the
DEM nodes inside it.

Two calls into ``_core`` and their ``--stats`` rows; no geometry. ``generate``
runs on both paths while the tile is held; ``run`` is the projected path's
strip run (``refine_strip``). On the reprojected path the strip goes to
``final_check.run`` instead, which joins it to the source nodes in one loop.
"""

from __future__ import annotations

import numpy as np

from tin_engine._core import (
    ConstraintCheckPoints,
    LineTolerance,
    PointRefineOutcome,
    RasterView,
    RefineOutcome,
    SlopeTolerance,
    constraint_check_points,
    refine_strip,
)
from tin_engine.stats import PhaseClock


def generate(
    view: RasterView, start: RefineOutcome, clock: PhaseClock, frozen: int = 0
) -> ConstraintCheckPoints:
    """The strip on ``start``'s constraint edges, z from ``view``; none on an
    edge whose mask meets ``frozen``, which no point may split."""
    edges = np.asarray(start.edges)
    if frozen:
        edges = edges[(np.asarray(start.masks) & frozen) == 0]
    with clock.phase("edge strip: generate"):
        return constraint_check_points(view, np.asarray(start.vertices), edges)


def run(
    view: RasterView,
    strip: ConstraintCheckPoints,
    start: RefineOutcome,
    tolerance: float,
    clock: PhaseClock,
    feet: bool = False,
    field: LineTolerance | None = None,
    slope: SlopeTolerance | None = None,
    threads: int = 0,
    frozen: int = 0,
) -> PointRefineOutcome:
    """``refine_strip`` from ``start``; ``clock`` gets the run's own two times.
    ``feet`` puts points near a line onto it first; ``field``, if any, sets
    each triangle's tolerance; ``slope``, if any, holds each point to its
    slope's tolerance too; ``threads`` goes to the C++ as given (0 is every
    core). No vertex goes on an edge whose mask meets ``frozen``."""
    arrays = (start.vertices, start.triangles, start.z, start.valid, start.edges, start.masks)
    out = refine_strip(
        view,
        strip,
        *(np.asarray(a) for a in arrays),
        tolerance=tolerance,
        threads=threads,
        frozen_mask=frozen,
        constraint_feet=feet,
        field=field,
        slope=slope,
    )
    clock.add("edge strip: scan (parallel)", out.scan_seconds)
    clock.add("edge strip: split + flip (serial)", out.split_seconds)
    return out
