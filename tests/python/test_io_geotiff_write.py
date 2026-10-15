"""`io.geotiff_write.encode_geotiff`, read back by the core's own reader.

`docs/increments/51-open-core.md`, *9.1 Public (P1 unless marked)*, T9. The
encoder stays public while the suites that exercised it move out, so this
round trip is what keeps it covered.
"""

from __future__ import annotations

import io

import numpy as np

from tin_engine.io.geotiff import decode_dem
from tin_engine.io.geotiff_write import encode_geotiff


def test_a_float32_grid_reads_back_with_its_origin_step_and_crs() -> None:
    # 4 x 5, distinct values that float32 holds exactly, one negative.
    grid = (np.arange(20, dtype=np.float32).reshape(4, 5) * 1.5 - 3.0).astype(np.float32)
    x0, y0 = 500_000.0, 6_800_000.0  # a UTM 33N easting and northing in Norway
    tile = decode_dem(io.BytesIO(encode_geotiff(grid, x0, y0, 10.0, 25833)))
    assert tile.array.dtype == np.float32
    np.testing.assert_array_equal(tile.array, grid)
    meta = tile.meta
    assert (meta.rows, meta.cols) == (4, 5)
    assert (meta.x_min, meta.y_max) == (x0, y0)
    assert (meta.delta_x, meta.delta_y) == (10.0, 10.0)
    assert (meta.epsg, meta.crs, meta.geographic) == (25833, "EPSG:25833", False)
    assert meta.pixel_is_area is False
