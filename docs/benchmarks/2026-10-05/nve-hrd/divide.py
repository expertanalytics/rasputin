"""Where ours and NVE's differ, and how flat the ground is there (@perf, step 4).

Usage: divide.py DATA DEM HERE STATION [STATION ...]

For the largest piece of (ours minus NVE's) and of (NVE's minus ours): its
area, the length of its edge shared with the other polygon (within 15 m),
and the raw DEM along that shared edge, sampled every 20 m at the nearest
node: the median height and the share of samples whose 3x3 neighbourhood
spans under 0.5 m (flat ground, as on a lake surface or a bog)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from shapely.geometry import shape

from tin_engine.catchment import _plan
from tin_engine.dem_input import repository_for
from tin_engine.io.station_set import read_references
from tin_engine.mosaic import Bounds, assemble

data, dem, here = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
refs, _ = read_references(data / "reference.geojson")
repository, _ = repository_for((dem,))
for sid in sys.argv[4:]:
    f = json.loads((here / "catchments" / f"{sid}.geojson").read_text())
    f = f["features"][0] if "features" in f else f
    ours, ref = shape(f["geometry"]), refs[sid]
    print(f"== {sid}")
    for label, a, b in (("ours only", ours, ref), ("NVE's only", ref, ours)):
        d = a.difference(b)
        if d.is_empty:
            print(f"  {label}: none"); continue
        big = max(getattr(d, "geoms", [d]), key=lambda g: g.area)
        shared = big.boundary.intersection(b.buffer(15.0))
        n = int(shared.length // 20.0)
        if n < 2:
            print(f"  {label}: largest {big.area / 1e6:.2f} km2, shares no edge"); continue
        x0, y0, x1, y1 = big.bounds
        plan = _plan(repository.footprints(), Bounds(x_min=x0 - 50, y_min=y0 - 50, x_max=x1 + 50, y_max=y1 + 50))
        tile = assemble(plan, repository.load).tile
        m, z = tile.meta, np.asarray(tile.array, dtype=np.float64)
        pts = [shared.interpolate(t) for t in np.linspace(0, shared.length, n)]
        hs, flat = [], 0
        for p in pts:
            r, c = round((m.y_max - p.y) / m.delta_y), round((p.x - m.x_min) / m.delta_x)
            win = z[max(r - 1, 0) : r + 2, max(c - 1, 0) : c + 2]
            hs.append(z[r, c]); flat += int(win.max() - win.min() < 0.5)
        print(f"  {label}: largest {big.area / 1e6:.2f} km2 of {d.area / 1e6:.2f}; shared edge {shared.length / 1e3:.1f} km; "
              f"heights there median {np.median(hs):.1f} m, {np.min(hs):.1f}-{np.max(hs):.1f}; flat (3x3 under 0.5 m) {100 * flat / len(pts):.0f} %")
