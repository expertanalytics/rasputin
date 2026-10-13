"""Increment 40's probe: flood the ANADEM fixture on its own grid in degrees.

Run from the worktree root as

    <python with tin_engine and a built _core> native_flood.py <dir holding tin_engine>

The `sys.meta_path` line drops an editable install of another checkout, which
would otherwise answer `import tin_engine`. Nothing here changes production
code: the core is called on a unit lattice (x = column, y = -row), which is
what the design's `raster.to_core_lattice` will build, and the areas are
computed here the way the design's `geodesy` module will.

It prints the fixture's grid, the cell in metres at the centre row, the flood
from the seed of the fixture's reference catchment (`tests/fixtures/velhas/
NOTICE`, section 2: seed 597180 7852380 in EPSG:31983), the flood and
accumulation times, the area by row-weighted node count and by the geodesic
area of the traced outline, and the overlap with the reference catchment that
the resampled route produced.
"""

# ruff: noqa: E402 -- the imports follow the sys.path set-up on purpose

import json
import sys
import time
from pathlib import Path

sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, sys.argv[1])

import numpy as np
from shapely.geometry import Polygon, shape

from tin_engine import _core
from tin_engine.crs import parse_crs, reprojector
from tin_engine.dem_input import repository_for
from tin_engine.outline import trace

repo, _ = repository_for((Path("tests/fixtures/velhas/anadem_velhas.tif"),))
fp = repo.footprints()[0]
m = fp.meta
arr = np.asarray(repo.load(fp.name).array)
print("grid", m.crs, m.rows, "x", m.cols, "spacing", m.delta_x, m.delta_y)
geod = parse_crs(m.crs).get_geod()
lat = m.y_max - m.rows / 2 * m.delta_y
_, _, dx = geod.inv(0, lat, m.delta_x, lat)
_, _, dy = geod.inv(0, lat - m.delta_y / 2, 0, lat + m.delta_y / 2)
print(f"centre latitude {lat:.4f}: cell {dx:.4f} m east-west, {dy:.4f} m north-south")

((x, y),) = reprojector("EPSG:31983", m.crs)([(597180, 7852380)])
r, c = map(round, m.index_of(x, y))
print("seed", x, y, "node", r, c)
view = _core.raster_view(arr, x_min=0.0, y_max=0.0, delta_x=1.0, delta_y=1.0, nodata=m.nodata)
seed = np.zeros(arr.shape, np.uint8)
seed[r, c] = 1
t = time.perf_counter()
out = _core.upstream(view, seed)
t_up = time.perf_counter() - t
t = time.perf_counter()
_core.accumulate(view)
t_acc = time.perf_counter() - t
print(f"upstream {t_up * 1e3:.1f} ms, accumulate {t_acc * 1e3:.1f} ms; nodes {out.nodes_in}, "
      f"edge {out.touches_edge}, nodata {out.touches_nodata}")

mask = np.asarray(out.mask).reshape(arr.shape)


def row_area(row: int) -> float:
    top = m.y_max - row * m.delta_y + m.delta_y / 2
    lons = [0, m.delta_x, m.delta_x, 0]
    a, _ = geod.polygon_area_perimeter(lons, [top, top, top - m.delta_y, top - m.delta_y])
    return abs(a)


areas = np.array([row_area(k) for k in range(m.rows)])
print(f"row-weighted node area {(mask.sum(axis=1) * areas).sum() / 1e6:.4f} km2")
rings = trace(mask)
ring = max(rings, key=len)
lon, la = m.node_xy(ring[:, 0], ring[:, 1])
a, _ = geod.polygon_area_perimeter(lon, la)
print(f"rings {len(rings)}, outline vertices {len(ring)}, "
      f"geodesic outline area {abs(a) / 1e6:.4f} km2")
ref = shape(json.loads(Path("tests/fixtures/velhas/catchment.geojson").read_text()))
ours = Polygon(reprojector(m.crs, "EPSG:31983")(np.column_stack([lon, la])))
inter = ours.intersection(ref).area
print(f"reference {ref.area / 1e6:.4f} km2 (EPSG:31983); "
      f"shared / reference {inter / ref.area:.4f}, "
      f"shared / ours {inter / ours.area:.4f}")
