"""The catchment's calls into ``_core``.

One of the adapters that alone may import ``_core`` (layer 3 in
``tests/python/test_layering.py``), as ``edge_strip.py`` is for the edge
strip: the flood and the flow accumulation over a decoded tile, each the core
call on ``raster.to_core_lattice``'s view, and the ring reduction re-exported.
The view is a unit lattice (spacing 1, no coordinates), since the hydrology
reads only the grid's shape, values and NoData; so a DEM in degrees floods as
one in metres. ``catchment`` imports this module, never ``_core`` or
``raster``, so it never holds a core raster view.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from tin_engine import _core
from tin_engine._core import AccumulateOutcome, ReduceStatus, UpstreamOutcome, reduce_ring
from tin_engine.io.models import DemTile
from tin_engine.raster import to_core_lattice


def upstream(tile: DemTile, seed: npt.NDArray[np.uint8]) -> UpstreamOutcome:
    """The nodes of ``tile`` that drain to the nodes ``seed`` marks."""
    return _core.upstream(to_core_lattice(tile), seed)


def accumulate(tile: DemTile) -> AccumulateOutcome:
    """Each node's drainage over ``tile``, by the Priority-Flood that fills its pits.

    Runs on the node grid alone: spacing and CRS are unused. A node drains to
    the neighbour that flooded it, its lowest filled neighbour, ties first in,
    first out; not D8's steepest descent. ``.count`` (uint32, rows x cols) is 0
    on NoData, else the number of nodes draining through the node, itself
    included. ``.flow_to`` (uint8) is that neighbour,
    ``3*(dr+1)+(dc+1)`` with ``dr``, ``dc`` its row and column minus the
    node's; 255 for an outlet or NoData. Both are read-only. A tile of 2^32
    nodes or more raises ``ValueError``. Releases the GIL. The result's type
    and its other attributes are not promised.
    """
    return _core.accumulate(to_core_lattice(tile))


__all__ = ["ReduceStatus", "UpstreamOutcome", "accumulate", "reduce_ring", "upstream"]
