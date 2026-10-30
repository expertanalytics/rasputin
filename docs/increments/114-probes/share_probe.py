"""Test 5's shares, measured four ways: the box's outline (2000 steps a side)
or its bounds, against the union of the bowl's four tiles' `_lonlat`
outlines (what `sources` holds) or against `COVERAGE` as a plain degree
rectangle (what `box_probe.py` used).

    PYTHONPATH=tests/python:src_python python docs/increments/114-probes/share_probe.py
"""

import numpy as np
import shapely

from catchment_fixtures import bowl, tile_of
from dataset_fixtures import COVERAGE
from mosaic_fixtures import quadrants
from tin_engine import source_report
from tin_engine.crs import reprojector

STEPS = 2000
to_lonlat = reprojector("EPSG:25833", "EPSG:4326")
tiles = quadrants(tile_of(bowl()), row_cut=200, col_cut=100, overlap=1)
union = shapely.union_all(
    [source_report._lonlat(to_lonlat, t.meta.node_box()) for t in tiles.values()]
)
rect = shapely.box(*COVERAGE)


def outline(box: tuple[float, float, float, float]) -> shapely.Polygon:
    x0, y0, x1, y1 = box
    s = np.linspace(0.0, 1.0, STEPS + 1)
    xs = np.concatenate(
        [x0 + (x1 - x0) * s, np.full_like(s, x1), x1 - (x1 - x0) * s, np.full_like(s, x0)]
    )
    ys = np.concatenate(
        [np.full_like(s, y0), y0 + (y1 - y0) * s, np.full_like(s, y1), y1 - (y1 - y0) * s]
    )
    return shapely.Polygon(to_lonlat(np.column_stack([xs, ys])))


for name, box in [
    ("INSIDE", (508522.7, 6576808.3, 511380.3, 6582383.5)),
    ("OUTSIDE", (555607.8, 6651832.7, 561353.6, 6663057.1)),
    ("EDGE", (517060.4, 6573496.4, 522780.6, 6579094.1)),
]:
    ring = outline(box)
    env = shapely.box(*ring.bounds)
    print(
        f"{name}: tiles by area {union.intersection(ring).area / ring.area:.6f}"
        f", tiles by bounds {union.intersection(env).area / env.area:.6f}"
        f", rectangle by area {rect.intersection(ring).area / ring.area:.6f}"
        f", rectangle by bounds {rect.intersection(env).area / env.area:.6f}"
    )
