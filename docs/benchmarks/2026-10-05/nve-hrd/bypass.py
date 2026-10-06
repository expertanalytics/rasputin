"""Does the river's flow pass beside the placed node? (@perf, step 4)

Usage: bypass.py DATA DEM HERE HALF_M OUT_CSV [STATION ...]

For each station (default: every `uncertain` row of HERE/results.csv),
burns the reach as the batch does in a window of P +-HALF_M m, runs
`_core.accumulate` on the burnt window, and writes the placed node's count
and the largest count of any node at most 3 nodes (30 m) from it along each axis (km2; lower bounds when
the window cuts the catchment), and the same largest count in the raw DEM."""

from __future__ import annotations

import csv
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

data, dem, here, half, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), float(sys.argv[4]), Path(sys.argv[5])
rows = {r["station"]: r for r in csv.DictReader(open(here / "results.csv", encoding="utf-8"))}
ids = sys.argv[6:] or [s for s, r in rows.items() if r["class"] == "uncertain"]
stations, s_crs = read_stations(data / "stations.geojson")
segments, crs, _ = read_segments(data / "rivers.geojson")
repository, _ = repository_for((dem,))
move = reprojector(s_crs, crs)
by = {s.station: s for s in stations}
with out.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["station", "placed_km2", "near_max_km2", "raw_near_max_km2", "half_m"])
    for sid in ids:
        st = by[sid]
        ((x, y),) = move([(st.x, st.y)])
        p = place(Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river), segments)
        assert p is not None
        px, py = p.position
        plan = _plan(repository.footprints(), Bounds(x_min=px - half, y_min=py - half, x_max=px + half, y_max=py + half))
        tile = assemble(plan, repository.load).tile
        burnt, path = burn_reach(tile, p.reach)
        cell = tile.meta.delta_x * tile.meta.delta_y / 1e6
        b = np.asarray(accumulate(to_core(burnt)).count)
        a = np.asarray(accumulate(to_core(tile)).count)
        r, c = path.chain[path.placed]
        sl = (slice(max(r - 3, 0), r + 4), slice(max(c - 3, 0), c + 4))
        w.writerow([sid, f"{b[r, c] * cell:.4f}", f"{b[sl].max() * cell:.4f}", f"{a[sl].max() * cell:.4f}", half])
        fh.flush()
