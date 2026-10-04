"""The partition of a run's window into cells, one piece job per cell (increment 23c).

Pure arithmetic on the window: no DEM, no mesh, no machine. The rule is
`docs/increments/23-basin-scale.md`, "The partition (23c; B2, B3 and B14 as
ruled)", and the memory estimate is "The memory estimate and the defaults";
`docs/increments/23-probes/partition.py` implements the same rule and its
figures are the design's.

Neither the count nor the cut reads the machine (K5): the budget is a
constant the caller passes, never the machine's memory or cores.
`lines` is "The chains": the partition lines clipped to the domain.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry import LineString, Polygon

#: `--memory-budget`'s default, 16 GiB in bytes (B14 as ruled).
DEFAULT_MEMORY_BUDGET = 2**34

#: `b(T)`: bytes per window node at tolerance T (m), from the basin-piece
#: sweep (`docs/benchmarks/2026-10-01/basin-piece/README.md`).
_BYTES_PER_NODE = (
    (0.0, 637),
    (1.0, 267),
    (2.0, 155),
    (5.0, 65),
    (10.0, 37),
    (20.0, 25),
    (50.0, 19),
)


def bytes_per_node(tolerance: float) -> int:
    """`b(T)`: linear between the table's columns, rounded up; 19 above 50 m."""
    for (t0, b0), (t1, b1) in itertools.pairwise(_BYTES_PER_NODE):
        if tolerance <= t1:
            return math.ceil(b0 + (b1 - b0) * (tolerance - t0) / (t1 - t0))
    return _BYTES_PER_NODE[-1][1]


@dataclass(frozen=True, slots=True)
class Partition:
    """`nx` x `ny` cells of `dx` x `dy` lattice nodes; the last ones may be narrower."""

    nx: int
    ny: int
    dx: int
    dy: int


def partition(
    cols: int,
    rows: int,
    tolerance: float,
    *,
    pieces: int = 1,
    memory_budget: int = DEFAULT_MEMORY_BUDGET,
) -> Partition:
    """The cells of a `cols` x `rows` window: the design's steps 1-3.

    The count approximates `max(pieces, ceil(N b / B))` and is not a floor;
    only `dx dy b <= memory_budget` is guaranteed. Refuses an empty window,
    `pieces < 1`, and a budget under one node's bytes, for which no cell fits.
    """
    if cols < 1 or rows < 1:
        raise ValueError(f"the window is empty: cols {cols} and rows {rows} must be at least 1")
    if pieces < 1:
        raise ValueError(f"pieces must be at least 1, got {pieces}")
    b = bytes_per_node(tolerance)
    if memory_budget < b:
        raise ValueError(
            f"the memory budget of {memory_budget} bytes is under one node's {b} bytes"
        )
    wanted = max(pieces, -(-cols * rows * b // memory_budget))
    if wanted <= 1:
        return Partition(1, 1, cols, rows)
    nx = max(1, round(math.sqrt(wanted * cols / rows)))
    ny = max(1, round(wanted / nx))
    dx, dy = -(-cols // nx), -(-rows // ny)
    # Ends: it runs only while a cell is over budget, and the budget is at
    # least one node's bytes (refused above otherwise), so a 1 x 1 cell fits.
    # A step adds one to nx only while dx >= dy and the cell is over, so
    # dx > 1, i.e. nx < cols; likewise ny only while dy > 1, i.e. ny < rows.
    # So nx stays at most cols and ny at most rows: at most cols + rows - 2 steps.
    while dx * dy * b > memory_budget:
        if dx >= dy:
            nx += 1
        else:
            ny += 1
        dx, dy = -(-cols // nx), -(-rows // ny)
    return Partition(-(-cols // dx), -(-rows // dy), dx, dy)


def lines(
    polygon: Polygon, xs: Sequence[float], ys: Sequence[float]
) -> list[npt.NDArray[np.float64]]:
    """The vertical lines `x = xs[i]` and the horizontal `y = ys[j]` clipped to
    `polygon` as linework (16b R6), each part an `(n, 2)` chain. A part's
    constant coordinate is set back to the line's own value, so a crossing
    GEOS computes with the outline stays exactly on the lattice line; a line
    that only touches the polygon gives no part."""
    x0, y0, x1, y1 = polygon.bounds
    out = []
    for axis, values in ((0, xs), (1, ys)):
        for v in values:
            line = LineString([(v, y0), (v, y1)] if axis == 0 else [(x0, v), (x1, v)])
            for part in shapely.get_parts(shapely.intersection(line, polygon)):
                if isinstance(part, LineString) and part.length > 0:
                    xy = shapely.get_coordinates(part)
                    xy[:, axis] = v
                    out.append(xy)
    return out
