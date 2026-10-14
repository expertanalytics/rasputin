# ruff: noqa: E402, E501, E702, N812
"""Increment 44 probe 2: today's row_areas against geodesy.cell_areas_by_row, bit for bit,
on the real MapBiomas cache header and the test fixture; GLO-30 nominal box per band."""
import io
import os
import sys

os.environ["PROJ_NETWORK"] = "OFF"
tree, data = sys.argv[1], sys.argv[2]
sys.path[:0] = [tree + "/tests/python"]
import tin_engine

assert tin_engine.__file__.startswith(tree)
from pathlib import Path

import numpy as np

import class_fixtures as K
from tin_engine.geodesy import cell_areas_by_row
from tin_engine.io.geotiff import read_page
from tin_engine.io.models import ClassWindow, IndexWindow
from tin_engine.io.repository import CacheRepository
from tin_engine.land_cover_fractions import row_areas


def compare(meta, label):
    worst, rows_checked = 0, 0
    for r0 in sorted({0, meta.rows // 3, meta.rows // 2, meta.rows - 3000 if meta.rows > 3000 else 0}):
        n = min(3000, meta.rows - r0)
        w = IndexWindow(row0=r0, col0=0, rows=n, cols=1)
        x, y = meta.node_xy(r0, 0)
        cw = ClassWindow(x_min=x - meta.delta_x / 2, y_max=y + meta.delta_y / 2, delta_x=meta.delta_x,
                         delta_y=meta.delta_y, array=np.zeros((n, 1), np.uint8))
        a, b = row_areas(cw), cell_areas_by_row(meta.windowed(w))
        worst = max(worst, int((a != b).sum())); rows_checked += n
        rel = float(np.max(np.abs(a - b) / a))
    print(f"{label}: rows checked {rows_checked}, rows differing in a bit (worst window) {worst}, max rel diff {rel:.3g}")

repo = CacheRepository(Path(data) / "cache", "mapbiomas-c11", nodata=0.0)
fps = repo.footprints()
for fp in fps[:2]:
    m = fp.meta
    print(fp.name, "crs", m.crs, "rows", m.rows, "cols", m.cols, "delta", repr(m.delta_x), repr(m.delta_y), "y_max", repr(m.y_max))
    compare(m, fp.name)
meta, _, _ = read_page(io.BytesIO(K.class_tiff()))
compare(meta, "fixture 1/128 deg")

print("== GLO-30 nominal growth per latitude band (x: 4 dy + 2 dx, y: 6 dy), degrees")
bands = [("0-50", 1), ("50-60", 1.5), ("60-70", 2), ("70-80", 3), ("80-85", 5), ("85-90", 10)]
for name, (dx, dy) in [("today 0.001, 0.001", (0.001, 0.001)), ("design 0.0027778, 0.00027778", (0.0027778, 0.00027778))]:
    gx, gy = 4 * dy + 2 * dx, 6 * dy
    short = [b for b, s in bands if 4 / 3600 + 2 * s / 3600 > gx or gy < 6 / 3600]
    print(f"{name}: grows x {gx:.6f}, y {gy:.6f}; bands whose own growth exceeds it: {short or 'none'}")
