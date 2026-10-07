"""The catchment outline grown by a distance with mitred corners (increment
30d, `docs/increments/30d-outline-buffer-speed.md` section 3.3).

GEOS's buffer is slow on a long raster-traced outline (a staircase) once the
distance spans two or more steps. There the region is built as the polygon
united with its boundary grown in overlapping pieces, P + B = P u (dP + B);
everywhere else it is GEOS's own buffer, bit for bit.
"""

from __future__ import annotations

import numpy as np
import shapely
from shapely.geometry import Polygon

#: Edges per piece of the ring; each piece overlaps the next by one edge.
PIECE_EDGES = 1000
#: The fewest exterior vertices (closing vertex not counted) grown in pieces.
PIECES_FROM = 5000


def uses_pieces(polygon: Polygon, distance: float) -> bool:
    """The gate: no holes, at least `PIECES_FROM` vertices, and a median edge
    of at most half the distance (the distance spans two steps or more)."""
    xy = shapely.get_coordinates(polygon.exterior)
    if polygon.interiors or len(xy) - 1 < PIECES_FROM:
        return False
    return bool(np.median(np.hypot(*np.diff(xy, axis=0).T)) <= distance / 2)


def grow_mitred(polygon: Polygon, distance: float) -> Polygon:
    """`polygon` grown by `distance` with mitred corners: GEOS's buffer off
    the gate, the pieced union on it."""
    if not uses_pieces(polygon, distance):
        grown = polygon.buffer(distance, join_style="mitre")
    else:
        xy = shapely.get_coordinates(polygon.exterior)
        n = len(xy) - 1
        pieces = []
        for start in range(0, n, PIECE_EDGES):
            stop = start + PIECE_EDGES + 2
            wrapped = xy[start:stop] if stop <= n + 1 else np.vstack([xy[start:], xy[1 : stop - n]])
            pieces.append(shapely.linestrings(wrapped))
        lines = shapely.buffer(np.array(pieces), distance, join_style="mitre", cap_style="flat")
        grown = shapely.union_all(np.concatenate([[polygon], lines]))
    assert isinstance(grown, Polygon)
    return grown
