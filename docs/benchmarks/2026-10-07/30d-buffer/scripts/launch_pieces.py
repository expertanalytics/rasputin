"""Candidate 3 of docs/increments/perf-audit.md (branch worktree-perf-audit):
a polygon grown by d is the polygon united with its exterior grown by d; the
exterior is grown in pieces of `size` edges, each overlapping the next by one
edge, with flat ends. Used by launch.py's --patch 3 and by candidates.py."""

import numpy as np
import shapely
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry


def piecewise_buffer(poly: Polygon, distance: float, size: int) -> BaseGeometry:
    """The polygon united with its exterior grown in pieces of `size` edges."""
    xy = shapely.get_coordinates(poly.exterior)
    n = len(xy) - 1
    lines = []
    for s in range(0, n, size):
        stop = s + size + 2
        piece = xy[s:stop] if stop <= n + 1 else np.vstack([xy[s:], xy[1 : stop - n]])
        lines.append(shapely.linestrings(piece))
    grown = shapely.buffer(np.array(lines), distance, join_style="mitre", cap_style="flat")
    return shapely.union_all(np.concatenate([[poly], grown]))
