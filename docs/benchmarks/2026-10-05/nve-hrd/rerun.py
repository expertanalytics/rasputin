"""Step 4's sensitivity re-runs (@perf): `rasputin station-catchments` on the
stations named, with the burn's corridor and the map radius changed.

Usage: rerun.py CORRIDOR MAP_RADIUS DATA DEM OUT_DIR STATION [STATION ...]

The map radius is the command's own `--map-radius`. The corridor (the burn's
half-width, `Reach.corridor`, 30 m by default) has no option, so this script
wraps `gauge.place` as `catchment_batch` imported it and returns the same
placement with the reach's corridor replaced; nothing else changes. The
command then runs in this process, as `rasputin` would.
"""

from __future__ import annotations

import sys

import tin_engine.catchment_batch as batch
from tin_engine.cli import app

corridor, radius = float(sys.argv[1]), sys.argv[2]
data, dem, out, only = sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:]
original = batch.place


def place(*args, **kwargs):  # type: ignore[no-untyped-def]
    p = original(*args, **kwargs)
    if p is None:
        return None
    return p.model_copy(update={"reach": p.reach.model_copy(update={"corridor": corridor})})


batch.place = place
argv = ["station-catchments", "--dem", dem, "--stations", f"{data}/stations.geojson"]
argv += ["--rivers", f"{data}/rivers.geojson", "--reference", f"{data}/reference.geojson"]
argv += ["--lakes", f"{data}/lakes.geojson", "--map-radius", radius, "--out-dir", out]
for s in only:
    argv += ["--only", s]
try:
    app(argv)
except SystemExit as exc:
    if exc.code not in (0, None):
        raise
