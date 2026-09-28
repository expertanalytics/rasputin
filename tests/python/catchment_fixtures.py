"""Terrains, lakes and an in-memory repository for increment 22's suites.

`docs/increments/22-auto-catchment.md`, "The red suites" (PR 1). Imports only
what earlier increments shipped, so a missing 22 module fails the tests that
use it, not the collection of this helper.

Every terrain is on a 100 m lattice, so the window's 2000 m margin
(`WINDOW_MARGIN_M`) is 20 cells and a catchment can outgrow the first window
on a grid of a few hundred nodes.

THE BOWL. An elliptical basin: z = 1000 (1 - |1 - b|) with
b = ((c - cc) / ax)^2 + ((r - rc) / ay)^2, rising from the centre to a rim at
b = 1 and falling outside it to the data's edges. A small tilt breaks the
mirror ties. The lake is a 7 x 7-node box round the centre, and its outlet is
a channel down column cc from the centre to the south edge, strictly
descending and below the lake: the lake is not a closed depression, as a real
lake with a river out of it is not. What drains into the lake is the part of
the basin that reaches it before the channel.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Mapping

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry import Polygon, box

from mosaic_fixtures import whole
from tin_engine.io.models import DemTile
from tin_engine.io.repository import TileFootprint

X0 = 500_000.0
Y0 = 6_600_000.0
D = 100.0
EPSG = "EPSG:25833"
MARGIN_CELLS = 20  # WINDOW_MARGIN_M / D


def lat(col: float, row: float) -> tuple[float, float]:
    """A lattice position (col, row) as a point in EPSG:25833."""
    return (X0 + D * col, Y0 - D * row)


def bowl(
    rows: int = 400,
    cols: int = 200,
    rc: int = 200,
    cc: int = 100,
    ax: float = 8.0,
    ay: float = 55.0,
) -> npt.NDArray[np.float32]:
    r, c = np.indices((rows, cols)).astype(np.float64)
    b = ((c - cc) / ax) ** 2 + ((r - rc) / ay) ** 2
    z = 1000.0 * (1.0 - np.abs(1.0 - b)) + 0.0131 * c + 0.0077 * r
    channel = np.arange(rows) >= rc
    z[channel, cc] = -1.0 - 0.5 * (np.arange(rows)[channel] - rc)
    return z.astype(np.float32)


def tile_of(z: npt.NDArray[np.floating], nodata: float | None = None) -> DemTile:
    rows, cols = z.shape
    return whole(rows, cols, array=z, x_min=X0, y_max=Y0, dx=D, dy=D, nodata=nodata)


def lake_box(rc: int = 200, cc: int = 100, half: float = 3.5) -> Polygon:
    """The lake: a box round (rc, cc), its sides half a cell from any node."""
    x0, y0 = lat(cc - half, rc + half)
    x1, y1 = lat(cc + half, rc - half)
    return box(x0, y0, x1, y1)


class MemoryRepository:
    """A `DemRepository` over tiles held in memory; records every load."""

    def __init__(self, tiles: Mapping[str, DemTile]) -> None:
        self.tiles = dict(tiles)
        self.loads: list[str] = []

    def footprints(self) -> tuple[TileFootprint, ...]:
        return tuple(
            TileFootprint(name=name, meta=tile.meta, dtype=tile.array.dtype)
            for name, tile in sorted(self.tiles.items())
        )

    def load(self, name: str) -> DemTile:
        self.loads.append(name)
        return self.tiles[name]


def node_grid(tile: DemTile) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    m = tile.meta
    r, c = np.indices((m.rows, m.cols)).astype(np.float64)
    return m.x_min + c * m.delta_x, m.y_max - r * m.delta_y


def seeds_in(tile: DemTile, lake: shapely.Geometry) -> npt.NDArray[np.uint8]:
    """The design's seed mask: DEM nodes inside the lake, `shapely.contains_xy`."""
    x, y = node_grid(tile)
    return shapely.contains_xy(lake, x, y).astype(np.uint8)


def full_flood(tile: DemTile, seed: npt.NDArray[np.uint8]) -> npt.NDArray[np.uint8]:
    """One flood over the whole raster: the reference the window loop must equal."""
    from tin_engine._core import upstream
    from tin_engine.raster import to_core

    return np.asarray(upstream(to_core(tile), seed).mask, dtype=np.uint8)


def filled(mask: npt.ArrayLike) -> npt.NDArray[np.uint8]:
    """`mask` with its holes filled: every out-node the padding cannot reach
    through 4-connected out-nodes becomes in (the (8, 4) pairing)."""
    m = np.pad(np.asarray(mask) != 0, 1)
    outside = np.zeros_like(m)
    outside[0, 0] = True
    queue = deque([(0, 0)])
    rows, cols = m.shape
    while queue:
        r, c = queue.popleft()
        for rr, cc in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
            if 0 <= rr < rows and 0 <= cc < cols and not m[rr, cc] and not outside[rr, cc]:
                outside[rr, cc] = True
                queue.append((rr, cc))
    return (~outside[1:-1, 1:-1]).astype(np.uint8)


def placed(
    mask: npt.ArrayLike, window_x_min: float, window_y_max: float, shape: tuple[int, int]
) -> npt.NDArray[np.uint8]:
    """A window's mask on the whole raster's lattice (origin X0, Y0, spacing D)."""
    m = np.asarray(mask, dtype=np.uint8)
    r0 = round((Y0 - window_y_max) / D)
    c0 = round((window_x_min - X0) / D)
    out = np.zeros(shape, dtype=np.uint8)
    out[r0 : r0 + m.shape[0], c0 : c0 + m.shape[1]] = m
    return out
