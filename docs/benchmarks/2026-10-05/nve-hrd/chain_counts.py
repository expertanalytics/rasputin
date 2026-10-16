"""Raw and burnt counts along a station's burnt chain (@perf, step 4).

Usage: chain_counts.py DATA DEM HALF_M CORRIDOR STATION [STATION ...]

Places the gauge as the batch does, burns the reach (`burn.burn_reach`,
corridor CORRIDOR m) in a window of P +-HALF_M, runs `_core.accumulate` on
the raw and the burnt window, and prints, for every chain node within 300 m
of the placed node, its arc (m, + downstream), the raw count and the burnt
count in km2 (lower bounds if the window cuts the catchment)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from tin_engine._core import accumulate
from tin_engine.burn import burn_reach
from tin_engine.catchment import _plan
from tin_engine.crs import reprojector
from tin_engine.dem_input import repository_for
from tin_engine.gauge import Gauge, place
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_stations
from tin_engine.mosaic import Bounds, assemble
from tin_engine.raster import to_core

data, dem, half, corridor = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
stations, s_crs = read_stations(data / "stations.geojson")
segments, crs, _ = read_segments(data / "rivers.geojson")
repository, _ = repository_for((dem,))
move = reprojector(s_crs, crs)
by = {s.station: s for s in stations}
for sid in sys.argv[5:]:
    st = by[sid]
    ((x, y),) = move([(st.x, st.y)])
    p = place(Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river), segments)
    assert p is not None
    reach = p.reach.model_copy(update={"corridor": corridor})
    px, py = p.position
    plan = _plan(repository.footprints(), Bounds(x_min=px - half, y_min=py - half, x_max=px + half, y_max=py + half))
    tile = assemble(plan, repository.load).tile
    burnt, path = burn_reach(tile, reach)
    raw = np.asarray(accumulate(to_core(tile)).count)
    brn = np.asarray(accumulate(to_core(burnt)).count)
    cell = tile.meta.delta_x * tile.meta.delta_y / 1e6
    print(f"== {sid} {st.name}: corridor {corridor:g} m, {len(path.chain)} chain nodes, placed {path.placed}")
    for k, (r, c) in enumerate(path.chain.tolist()):
        if abs(path.arc[k]) <= 300.0:
            z0, z1 = float(tile.array[r, c]), float(burnt.array[r, c])
            print(f"  {path.arc[k]:7.1f} m  raw {raw[r, c] * cell:9.3f}  burnt {brn[r, c] * cell:9.3f} km2"
                  f"  z {z0:8.2f} -> {z1:8.2f}{'  <- placed' if k == path.placed else ''}")
