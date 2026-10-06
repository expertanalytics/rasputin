"""The catchment's calls into ``_core`` (python-audit.md, section 11).

Layer 3, as ``edge_strip.py`` is the edge strip's: the flood and the flow
accumulation over a decoded tile, each the core call on ``raster.to_core``'s
view, and the ring reduction re-exported. ``catchment`` imports this module,
never ``_core`` or ``raster``, so it never holds a core raster view.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from tin_engine import _core
from tin_engine._core import AccumulateOutcome, ReduceStatus, UpstreamOutcome, reduce_ring
from tin_engine.io.models import DemTile
from tin_engine.raster import to_core


def upstream(tile: DemTile, seed: npt.NDArray[np.uint8]) -> UpstreamOutcome:
    """The nodes of ``tile`` that drain to the nodes ``seed`` marks."""
    return _core.upstream(to_core(tile), seed)


def accumulate(tile: DemTile) -> AccumulateOutcome:
    """The flow accumulation over ``tile``."""
    return _core.accumulate(to_core(tile))


__all__ = ["ReduceStatus", "UpstreamOutcome", "accumulate", "reduce_ring", "upstream"]
