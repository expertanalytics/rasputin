"""Where the river above a gauge drains in a small window (@perf, step 4).

Usage: trace_exit.py DATA DEM BIG_HALF SMALL_HALF STATION

In a window of P +-BIG_HALF m, follows the main stem up from the placed node
(at each node the upstream neighbour of largest raw count) for 400 nodes.
Then, in the window P +-SMALL_HALF m, starts at the first of those stem
nodes 1500 m or more from P that lies inside, follows `flow_to` down, and
says where the path ends (a window edge, or the placed node) and the
lowest and highest raw heights on it. Raw DEM, no burn."""

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

data, dem = Path(sys.argv[1]), Path(sys.argv[2])
big, small, sid = float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
stations, s_crs = read_stations(data / "stations.geojson")
segments, crs, _ = read_segments(data / "rivers.geojson")
repository, _ = repository_for((dem,))
st = {s.station: s for s in stations}[sid]
((x, y),) = reprojector(s_crs, crs)([(st.x, st.y)])
p = place(Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river), segments)
assert p is not None
px, py = p.position


def window(half):  # type: ignore[no-untyped-def]
    plan = _plan(repository.footprints(), Bounds(x_min=px - half, y_min=py - half, x_max=px + half, y_max=py + half))
    tile = assemble(plan, repository.load).tile
    _, path = burn_reach(tile, p.reach)
    m = tile.meta
    r, c = path.chain[path.placed]
    return tile, accumulate(to_core(tile)), (m.x_min + c * m.delta_x, m.y_max - r * m.delta_y)


N8 = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if (dr, dc) != (0, 0)]
tile, acc, (gx, gy) = window(big)
m = tile.meta
count, flow = np.asarray(acc.count), np.asarray(acc.flow_to)
r, c = round((m.y_max - gy) / m.delta_y), round((gx - m.x_min) / m.delta_x)
stem = []
for _ in range(400):
    best = None
    for dr, dc in N8:
        rr, cc = r + dr, c + dc
        if 0 <= rr < count.shape[0] and 0 <= cc < count.shape[1]:
            if flow[rr, cc] == 3 * (-dr + 1) + (-dc + 1) and (best is None or count[rr, cc] > count[best]):
                best = (rr, cc)
    if best is None:
        break
    r, c = best
    stem.append((m.x_min + c * m.delta_x, m.y_max - r * m.delta_y, count[r, c] * 1e-4))
print(f"{sid}: placed node ({gx:.0f}, {gy:.0f}); in +-{big:g} m its count is "
      f"{count[round((m.y_max - gy) / m.delta_y), round((gx - m.x_min) / m.delta_x)] * 1e-4:.2f} km2")
tile2, acc2, _ = window(small)
m2 = tile2.meta
flow2, z2 = np.asarray(acc2.flow_to), np.asarray(tile2.array)
start = next(s for s in stem if np.hypot(s[0] - gx, s[1] - gy) >= 1500.0
             and m2.x_min < s[0] < m2.x_min + (m2.cols - 1) * m2.delta_x
             and m2.y_max - (m2.rows - 1) * m2.delta_y < s[1] < m2.y_max)
print(f"stem node ({start[0]:.0f}, {start[1]:.0f}), {np.hypot(start[0] - gx, start[1] - gy):.0f} m from it, "
      f"count {start[2]:.2f} km2 in the big window")
r, c = round((m2.y_max - start[1]) / m2.delta_y), round((start[0] - m2.x_min) / m2.delta_x)
gr, gc = round((m2.y_max - gy) / m2.delta_y), round((gx - m2.x_min) / m2.delta_x)
zs, steps = [], 0
while True:
    zs.append(float(z2[r, c]))
    if (r, c) == (gr, gc):
        print(f"  in +-{small:g} m: reaches the placed node after {steps} steps"); break
    code = int(flow2[r, c])
    if steps > 200000:
        print(f"  in +-{small:g} m: no end after {steps} steps (a cycle)"); break
    if code in (4, 255) or not (0 <= r < flow2.shape[0] and 0 <= c < flow2.shape[1]):
        print(f"  in +-{small:g} m: path ends at ({m2.x_min + c * m2.delta_x:.0f}, {m2.y_max - r * m2.delta_y:.0f}) "
              f"after {steps} steps, row {r} of {flow2.shape[0]}, col {c} of {flow2.shape[1]}, flow code {code}"); break
    r, c = r + code // 3 - 1, c + code % 3 - 1
    steps += 1
    if not (0 <= r < flow2.shape[0] and 0 <= c < flow2.shape[1]):
        print(f"  in +-{small:g} m: path leaves the window after {steps} steps"); break
print(f"  heights on the path: first {zs[0]:.2f}, lowest {min(zs):.2f}, highest {max(zs):.2f}, last {zs[-1]:.2f}")
