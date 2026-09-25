"""The mesh domain for a DEM before refinement exists: a stride subsample.

Increment 12, R1. Every ``stride``-th node of the DEM in both directions,
starting at node ``(0, 0)`` and always including the last row and column, so the
mesh covers the whole tile. The perimeter of that subsample is the outer ring;
every other node is a free vertex. A stand-in: refinement (ROADMAP gap 2)
chooses vertices by error budget instead.

Pure numpy. It imports neither ``_core`` nor the reader: a ``RasterMeta`` is all
it needs.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

from tin_engine.io.models import RasterMeta

#: The default subsample has at most this many nodes on its longer side.
MAX_NODES_PER_SIDE = 256


def default_stride(meta: RasterMeta) -> int:
    """The smallest stride giving at most :data:`MAX_NODES_PER_SIDE` per side."""
    longest = max(meta.rows, meta.cols)
    return max(1, math.ceil((longest - 1) / (MAX_NODES_PER_SIDE - 1)))


def _axis(count: int, stride: int) -> list[int]:
    """Every ``stride``-th index along one axis, plus the last."""
    picked = list(range(0, count, stride))
    return picked if picked[-1] == count - 1 else [*picked, count - 1]


def subsample(meta: RasterMeta, stride: int) -> tuple[npt.NDArray[np.float64], list[int]]:
    """The subsampled nodes as ``(N, 2)`` float64, and the perimeter ring.

    Nodes are row-major over the picked rows and columns. Coordinates use
    ``RasterGeometry::node``'s expression exactly, ``x_min + col * delta_x`` and
    ``y_max - row * delta_y``: any other spelling can land an ulp off the node
    and change which cell ``bilinear`` picks.

    The ring is counter-clockwise and implicitly closed (no repeated first
    index): along the bottom row, up the right column, back along the top row,
    and down the left column.
    """
    if stride < 1:
        raise ValueError(f"stride must be a positive integer, got {stride}")
    rows = _axis(meta.rows, stride)
    cols = _axis(meta.cols, stride)
    r = np.asarray(rows, dtype=np.float64)[:, None]
    c = np.asarray(cols, dtype=np.float64)[None, :]
    x = np.broadcast_to(meta.x_min + c * meta.delta_x, (len(rows), len(cols)))
    y = np.broadcast_to(meta.y_max - r * meta.delta_y, (len(rows), len(cols)))
    xy = np.column_stack([x.ravel(), y.ravel()])

    nr, nc = len(rows), len(cols)

    def at(i: int, j: int) -> int:
        return i * nc + j

    ring = (
        [at(nr - 1, j) for j in range(nc)]
        + [at(i, nc - 1) for i in range(nr - 2, -1, -1)]
        + [at(0, j) for j in range(nc - 2, -1, -1)]
        + [at(i, 0) for i in range(1, nr - 1)]
    )
    return xy, ring
