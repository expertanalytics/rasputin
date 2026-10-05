"""Does a station's catchment depend on the window it was computed in? (@perf)

Usage: window_check.py DATA DEM HERE HALF_M OUT_CSV [STATION ...]

For each station (default: every river-seeded row of HERE/results.csv that
is not refused), burns the reach as the batch does in a window of P
+-HALF_M m and reads the burnt count at the placed node (`_core.accumulate`).
The batch's own count is `a0`. A count larger here than the batch's, when
the batch's catchment was clear of its window's edge, means the routing
depends on the window. A count here is a lower bound when this window cuts
the catchment."""

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
ids = sys.argv[6:] or [s for s, r in rows.items() if r["seeded_by"] == "river" and r["class"] != "refused"]
stations, s_crs = read_stations(data / "stations.geojson")
segments, crs, _ = read_segments(data / "rivers.geojson")
repository, _ = repository_for((dem,))
move = reprojector(s_crs, crs)
by = {s.station: s for s in stations}
with out.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["station", "class", "batch_a0_km2", "count_km2", "half_m", "edge_cut", "ratio"])
    for sid in ids:
        st = by[sid]
        ((x, y),) = move([(st.x, st.y)])
        p = place(Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river), segments)
        assert p is not None
        px, py = p.position
        try:
            plan = _plan(repository.footprints(), Bounds(x_min=px - half, y_min=py - half, x_max=px + half, y_max=py + half))
            tile = assemble(plan, repository.load).tile
            burnt, path = burn_reach(tile, p.reach)
        except Exception as exc:  # noqa: BLE001  (a probe: say it and go on)
            print(sid, "skipped:", exc, file=sys.stderr)
            continue
        acc = accumulate(to_core(burnt))
        r, c = path.chain[path.placed]
        k = float(np.asarray(acc.count)[r, c]) * tile.meta.delta_x * tile.meta.delta_y / 1e6
        bits = int(np.asarray(acc.reach)[r, c])
        a0 = float(rows[sid]["a0"])
        w.writerow([sid, rows[sid]["class"], f"{a0:.4f}", f"{k:.4f}", half, bits & 1, f"{k / a0:.3f}" if a0 else ""])
        fh.flush()
        print(sid, rows[sid]["class"], f"{a0:.3f}", f"{k:.3f}", bits, file=sys.stderr)
