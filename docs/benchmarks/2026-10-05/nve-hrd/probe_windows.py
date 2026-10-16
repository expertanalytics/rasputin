"""Where a station's flood reached, window by window (@perf, step 4).

Usage: probe_windows.py DATA DEM OUT_DIR STATION

Runs `station-catchments --only STATION` and logs, for every flood, the
window's bounds and the in-nodes' bounding box in metres (EPSG:25833), by
wrapping `catchment._grow`'s flood callable. Nothing else changes."""

from __future__ import annotations

import sys

import tin_engine.catchment as c
from tin_engine.cli import app

data, dem, out, sid = sys.argv[1:5]
grow = c._grow


def logged(footprints, repository, bounds, flood, burnt):  # type: ignore[no-untyped-def]
    def wrapped(tile):  # type: ignore[no-untyped-def]
        o, kept = flood(tile)
        m = tile.meta
        x = (m.x_min + o.col_min * m.delta_x, m.x_min + o.col_max * m.delta_x)
        y = (m.y_max - o.row_max * m.delta_y, m.y_max - o.row_min * m.delta_y)
        print(f"window x {m.x_min:.0f}..{m.x_min + (m.cols - 1) * m.delta_x:.0f} "
              f"y {m.y_max - (m.rows - 1) * m.delta_y:.0f}..{m.y_max:.0f}; in-nodes "
              f"{o.nodes_in} x {x[0]:.0f}..{x[1]:.0f} y {y[0]:.0f}..{y[1]:.0f}; "
              f"edge {o.touches_edge} nodata {o.touches_nodata}", file=sys.stderr)
        return o, kept
    return grow(footprints, repository, bounds, wrapped, burnt)


c._grow = logged
argv = ["station-catchments", "--dem", dem, "--stations", f"{data}/stations.geojson",
        "--rivers", f"{data}/rivers.geojson", "--reference", f"{data}/reference.geojson",
        "--lakes", f"{data}/lakes.geojson", "--out-dir", out, "--only", sid]  # fmt: skip
try:
    app(argv)
except SystemExit as exc:
    if exc.code not in (0, None):
        raise
