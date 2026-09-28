"""The fine outline of a node mask: marching squares (increment 22, PR 1).

`docs/increments/22-auto-catchment.md`, "The fine outline". Every vertex is
the midpoint of a lattice edge between an in-node and an out-node; in each
square of four nodes a segment joins two such midpoints with the in-nodes on
its left in the world frame (x = col east, y = -row north). The one ambiguous
square, two in-nodes on a diagonal, keeps them joined: the catchment is
8-connected and the outside 4-connected (Kong and Rosenfeld 1989). The mask is
padded by one row and column of out-nodes, so every ring closes.

Each boundary midpoint starts exactly one segment and ends exactly one, so
the segments close into simple rings that share no point; outer rings are
counter-clockwise in the world frame and holes clockwise.

Numpy only: no `_core`, no paths.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

Ring = npt.NDArray[np.float64]

#: A square's corners counter-clockwise in the world frame, as (row, col)
#: offsets from its north-west node: south-west, south-east, north-east,
#: north-west. Edge k runs from corner k to corner k + 1.
_CORNERS = ((1, 0), (1, 1), (0, 1), (0, 0))


def trace(mask: npt.ArrayLike) -> list[Ring]:
    """Every ring of `mask`'s in-nodes (non-zero), as `(k, 2)` float arrays of
    `(row, col)` in the mask's own indices, each vertex a half on one axis.
    Not closed: the first vertex is not repeated."""
    padded = np.pad(np.asarray(mask) != 0, 1)
    rows, cols = padded.shape
    full = [padded[r : rows - 1 + r, c : cols - 1 + c] for r, c in _CORNERS]
    count: npt.NDArray[np.uint8] = np.sum([f.astype(np.uint8) for f in full], axis=0)
    # Only squares with in- and out-corners carry a segment.
    i, j = np.nonzero((count > 0) & (count < 4))
    corner = [f[i, j] for f in full]
    # A midpoint's key is its doubled padded (row, col), flattened: unique.
    width = 2 * cols
    mid = [
        (2 * i + _CORNERS[k][0] + _CORNERS[(k + 1) % 4][0]) * width
        + 2 * j
        + _CORNERS[k][1]
        + _CORNERS[(k + 1) % 4][1]
        for k in range(4)
    ]
    starts: list[Any] = []
    ends: list[Any] = []
    for k in range(4):
        # Edge k leaves the in-nodes (corner k in, k + 1 out): the segment
        # ends at the next edge counter-clockwise that enters them again.
        leaving = corner[k] & ~corner[(k + 1) % 4]
        after = corner[(k + 2) % 4]
        end = np.where(
            after,
            mid[(k + 1) % 4],
            np.where(corner[(k + 3) % 4], mid[(k + 2) % 4], mid[(k + 3) % 4]),
        )
        starts.append(mid[k][leaving])
        ends.append(end[leaving])
    start, stop = np.concatenate(starts), np.concatenate(ends)
    order = np.argsort(start, kind="stable")
    following = order[np.searchsorted(start, stop, sorter=order)]
    seen = np.zeros(len(start), dtype=bool)
    rings: list[Ring] = []
    for first in order:
        if seen[first]:
            continue
        loop = []
        s = int(first)
        while not seen[s]:
            seen[s] = True
            loop.append(start[s])
            s = int(following[s])
        keys = np.asarray(loop, dtype=np.int64)
        rings.append(np.column_stack([keys // width / 2 - 1, keys % width / 2 - 1]))
    return rings
