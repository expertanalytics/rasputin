"""Terrains, rivers and a hand burn for increment 29 PR 2's suites.

`docs/increments/29-nve-reference-catchments.md`, "Placing the gauge",
"Following the river in the DEM", "Sensitivity" and "The red suites" (PR 2).
Imports only what earlier increments shipped (PR 1's `accumulate`, PR 3's
river reader, 22's fixtures), so a missing PR 2 module fails the tests that
use it, not the collection of this helper.

Every terrain is on a 10 m lattice, DTM10's, because the design's numbers
(the 30 m floor of `U`, the 500 m end cap, 35 diagonal or 50 straight
extension steps) are DTM10's. The origin is the corner of a UTM 33 tile.

THE VALLEY. `valley()` is a basin with a rim, a valley floor down column
`CC` falling 0.5 m per row to the south, an exit slot through the southern rim
down column `CC` to the data's south edge, and, outside the rim, ground that
falls to the data's edges, so a catchment inside the rim is clear of every
edge. With `dam=True` an embankment crosses the valley at rows `DAM_ROWS`,
crest `DAM_CREST`, and a notch along row `NOTCH_ROW` leaves the valley floor
eastward through the rim to the east edge, far below the crest. Its first node
(column `CC + 1`) is `notch_drop` above the floor (0.5 m by default) and it
falls 0.01 m per column. In the raw DEM the valley above the embankment
drains out through the notch, so the catchment of a node below the
embankment is the lower valley only; burnt along column `CC`, the embankment
no longer holds the water and the whole valley drains past that node. With
`notch_drop=0` the notch's first node is 0.01 m below the floor beside the
chain, and takes the flow of the chain node above it. These claims were
measured on PR 1's flood when the suite was written, and the tests that rely
on one assert it as a premise.

THE HAND BURN. `hand_burn(z, chain)` is the design's descent, step 4, along a
given chain of nodes, in the array's own dtype: `z'[0] = z[0]`, then
`z'[k] = min(z[k], z'[k-1] - 0.001)`. It is the test's oracle for a chain the
test already knows, not a burn: it does not find the chain.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry import LineString

from mosaic_fixtures import whole
from tin_engine.io.models import DemTile

X0 = 500_000.0
Y0 = 6_600_000.0
D = 10.0
EPSG = "EPSG:25833"
CELL_KM2 = D * D / 1e6

ROWS, COLS, CC = 240, 61, 30
#: The rim: an ellipse round (RIM_ROW, CC) with these semi-axes in nodes.
RIM_ROW, RIM_ROWS, RIM_COLS = 120.0, 110.0, 20.0
DAM_ROWS = (99, 100, 101)
DAM_CREST = 290.0
NOTCH_ROW = 70
#: The row where the exit slot leaves the valley floor.
SLOT_ROW = 180

#: `flow_to`'s code for the neighbour (dr, dc): 3*(dr+1)+(dc+1).
DIRECTION = {(dr, dc): 3 * (dr + 1) + (dc + 1) for dr in (-1, 0, 1) for dc in (-1, 0, 1)}
OUTLET = 255


def lat(col: float, row: float) -> tuple[float, float]:
    """A lattice position (col, row) as a point in EPSG:25833."""
    return (X0 + D * col, Y0 - D * row)


def floor_z(row: float) -> float:
    """The valley floor's elevation at `row` (column `CC`)."""
    return 300.0 - 0.5 * row


def valley(
    dam: bool = False, notch_drop: float = 0.5, dtype: type = np.float64
) -> npt.NDArray[np.floating]:
    """THE VALLEY of the module's docstring, `ROWS` x `COLS` nodes."""
    r, c = np.indices((ROWS, COLS)).astype(np.float64)
    b = ((c - CC) / RIM_COLS) ** 2 + ((r - RIM_ROW) / RIM_ROWS) ** 2
    rim = np.where(b <= 1.0, 1000.0 * b**4, 1000.0 * (2.0 - b))
    floor = floor_z(r) + 3.0 * np.abs(c - CC)
    z = np.where(b <= 1.0, np.maximum(rim, floor), rim)
    slot = np.arange(ROWS) >= SLOT_ROW
    z[slot, CC] = floor_z(SLOT_ROW) - 2.0 * (np.arange(ROWS)[slot] - SLOT_ROW)
    if dam:
        inside = b[DAM_ROWS[0] : DAM_ROWS[-1] + 1] <= 1.0
        band = z[DAM_ROWS[0] : DAM_ROWS[-1] + 1]
        band[inside] = np.maximum(band[inside], DAM_CREST)
        east = np.arange(CC, COLS)
        z[NOTCH_ROW, east[1:]] = floor_z(NOTCH_ROW) + notch_drop - 0.01 * (east[1:] - CC)
    return z.astype(dtype)


#: The large terrain of `two_basins()`: the gauge's own valley down column
#: `BIG_CC`, and a tributary basin far to the east joining it at `JOIN_ROW`.
BIG_ROWS, BIG_COLS, BIG_CC, JOIN_ROW = 300, 700, 100, 152


def two_basins() -> npt.NDArray[np.float64]:
    """A valley like `valley()` (no embankment) down column `BIG_CC`, and a
    tributary basin spanning columns 140 to 660 whose floor falls west into a
    slot along row `JOIN_ROW` that joins the valley floor there. Ground
    outside both rims falls to the data's edges. The catchment of the valley
    node at row 150 lies within 2000 m of column `BIG_CC`; that of row 153
    takes in the tributary, more than 5 km east."""
    r, c = np.indices((BIG_ROWS, BIG_COLS)).astype(np.float64)
    edge = np.minimum(np.minimum(r, BIG_ROWS - 1 - r), np.minimum(c, BIG_COLS - 1 - c))
    background = 500.0 + 2.0 * edge
    bm = ((c - BIG_CC) / 20.0) ** 2 + ((r - 150.0) / 110.0) ** 2
    main = np.maximum(1000.0 * bm**4, floor_z(r) + 3.0 * np.abs(c - BIG_CC))
    bt = ((c - 400.0) / 260.0) ** 2 + ((r - JOIN_ROW) / 90.0) ** 2
    trib = np.maximum(1000.0 * bt**4, 240.0 + 0.02 * (c - 140.0) + 3.0 * np.abs(r - JOIN_ROW))
    z = np.minimum(background, np.where(bm <= 1.0, main, np.inf))
    z = np.minimum(z, np.where(bt <= 1.0, trib, np.inf))
    slot = np.arange(BIG_ROWS) >= SLOT_ROW
    z[slot, BIG_CC] = floor_z(SLOT_ROW) - 2.0 * (np.arange(BIG_ROWS)[slot] - SLOT_ROW)
    east = np.arange(BIG_CC + 1, 231)
    z[JOIN_ROW, east] = floor_z(JOIN_ROW) + 1.0 + 0.05 * (east - BIG_CC - 1)
    return z


def channel(
    rows: int,
    cols: int,
    path: Sequence[tuple[float, float]],
    *,
    top: float = 300.0,
    slope: float = 0.5,
    wall: float = 3.0,
    width: float = 0.0,
    dtype: type = np.float64,
) -> npt.NDArray[np.floating]:
    """A channel along `path`, (row, col) vertices in node units: at each node,
    `top - slope * s + wall * max(0, d - width)`, where `s` is the arc length
    (in nodes) of the node's nearest point on the path and `d` its distance
    from the path. `width` widens the floor to a band of that half-width, so
    nodes across it tie in elevation. A path that runs to the data's edge
    gives the channel an outlet there."""
    line = LineString([(c, r) for r, c in path])
    r, c = np.indices((rows, cols)).astype(np.float64)
    points = shapely.points(c.ravel(), r.ravel())
    s = shapely.line_locate_point(line, points).reshape(rows, cols)
    d = shapely.distance(line, points).reshape(rows, cols)
    return (top - slope * s + wall * np.maximum(0.0, d - width)).astype(dtype)


def tile_of(z: npt.NDArray[np.floating], nodata: float | None = None) -> DemTile:
    rows, cols = z.shape
    return whole(rows, cols, array=z, x_min=X0, y_max=Y0, dx=D, dy=D, nodata=nodata)


def column_line(col: float, row0: float, row1: float) -> tuple[tuple[float, float], ...]:
    """A straight mapped line down a column, from row0 to row1 (downstream)."""
    return (lat(col, row0), lat(col, row1))


def hand_burn(
    z: npt.NDArray[np.floating], chain: Sequence[tuple[int, int]]
) -> npt.NDArray[np.floating]:
    """The descent of "Following the river", step 4, along `chain`, in z's dtype."""
    out = np.array(z, copy=True)
    step = np.asarray(0.001, dtype=out.dtype)
    for k in range(1, len(chain)):
        lowered = out[chain[k - 1]] - step
        out[chain[k]] = min(out[chain[k]], lowered)
    return out


def drains_along(flow_to: npt.NDArray[np.uint8], chain: Sequence[tuple[int, int]]) -> int:
    """The index of the first chain node whose `flow_to` is not its successor,
    or len(chain) - 1 when every node but the last drains into the next."""
    for k in range(len(chain) - 1):
        (r0, c0), (r1, c1) = chain[k], chain[k + 1]
        if flow_to[r0, c0] != DIRECTION.get((r1 - r0, c1 - c0), -1):
            return k
    return len(chain) - 1


def is_taut(chain: Sequence[tuple[int, int]]) -> bool:
    """8-connected, and no two nodes that are not consecutive equal or 8-neighbours."""
    nodes = [tuple(int(v) for v in n) for n in chain]
    for k in range(len(nodes) - 1):
        (r0, c0), (r1, c1) = nodes[k], nodes[k + 1]
        if max(abs(r1 - r0), abs(c1 - c0)) != 1:
            return False
    for i in range(len(nodes)):
        for j in range(i + 2, len(nodes)):
            (r0, c0), (r1, c1) = nodes[i], nodes[j]
            if max(abs(r1 - r0), abs(c1 - c0)) <= 1:
                return False
    return True
