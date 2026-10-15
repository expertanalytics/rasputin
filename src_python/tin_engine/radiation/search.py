"""The horizon of every receiver, searched on the DEM around the mesh.

The near raster is the DEM within `d_near_m` of the mesh's box at full
resolution, read by the same `open_dem` the mesh is made with; coverage ends
where the files end. The search stops at `d_near_m`: no coarse far raster is
read. In the melt season terrain further away hides only a sun a few degrees
up, too weak to matter; the measurement is in
`docs/increments/46-hbv-snow.md`, *10.4 Which radiation steps, and how much of each*.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np

from tin_engine.dem_input import DemRequest, open_dem
from tin_engine.io.models import Bounds, DemTile
from tin_engine.radiation.frame import ReceiverFrame
from tin_engine.radiation.geometry import RadiationGeometry
from tin_engine.radiation.request import RadiationRequest
from tin_engine.raster import to_core
from tin_engine.resources import Resources


def near_raster(
    sources: Sequence[Path],
    box: tuple[float, float, float, float],
    request: RadiationRequest,
    resources: Resources | None = None,
) -> DemTile:
    """The DEM over `box` (x_min, y_min, x_max, y_max, in the DEM's CRS) grown
    by `request.d_near_m` on every side, as far as the files cover it, in
    float32. Refusals are `open_dem`'s."""
    d = request.d_near_m
    grown = Bounds.of((box[0] - d, box[1] - d, box[2] + d, box[3] + d))
    tile = open_dem(DemRequest(sources=tuple(sources), bounds=grown), resources).tile
    if tile.array.dtype == np.float32:
        return tile
    return DemTile(meta=tile.meta, array=tile.array.astype(np.float32))


def search_horizons(
    frame: ReceiverFrame, near: DemTile, request: RadiationRequest, *, threads: int = 0
) -> RadiationGeometry:
    """The horizons of `frame`'s receivers on `near`, to `request.d_near_m`,
    with the Earth's curvature less refraction at the frame's radius; `threads`
    threads (0: every core)."""
    return RadiationGeometry.searched(frame, to_core(near), request, threads)
