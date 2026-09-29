"""`_core.upstream`, the flood's binding (increment 22, PR 1).

`docs/increments/22-auto-catchment.md`, "Where it runs" and "Flow and
membership: one flood". The flood's behaviour is the C++ suite's
(`tests/cpp/unit/test_hydrology_upstream.cpp`); this file pins what crosses
the binding.

Interface assumed (names chosen here, stated in the handback; they mirror the
C++ `UpstreamOutcome`):

- `_core.upstream(view: RasterView, seed: NDArray[uint8]) -> UpstreamOutcome`,
  `seed` of the raster's shape `(rows, cols)`, non-zero for a seed. Any other
  shape is a `ValueError`. Releases the GIL.
- `UpstreamOutcome.mask`: a `uint8` array of shape `(rows, cols)`, 1 for an
  in-node and 0 otherwise, owned by the outcome (it outlives the view).
  `nodes_in`, `row_min`, `row_max`, `col_min`, `col_max` (inclusive),
  `touches_edge`, `touches_nodata`.
"""

from __future__ import annotations

import gc
from typing import Any

import numpy as np
import pytest

from tin_engine._core import raster_view

ROWS, COLS, RV, W = 21, 15, 10, 6


def valley(dtype: type = np.float32) -> np.ndarray:
    """The C++ suite's V-valley: falling west at 0.37 a column, sides rising at
    1 a row to ridges at rows 4 and 16, outer slopes at 2.9; column 0 a wall
    at 1000 except the seed (10, 0) at -1."""
    r, c = np.indices((ROWS, COLS))
    off = np.abs(r - RV).astype(dtype)
    side = np.where(off <= W, off, dtype(W) - dtype(2.9) * (off - dtype(W)))
    z = (dtype(0.37) * c.astype(dtype) + side).astype(dtype)
    z[:, 0] = 1000
    z[RV, 0] = -1
    z.flags.writeable = False
    return z


def view_of(z: np.ndarray, nodata: float | None = None) -> Any:
    return raster_view(
        z, x_min=500_000.0, y_max=6_600_000.0, delta_x=10.0, delta_y=5.0, nodata=nodata
    )


def seed_at(*nodes: tuple[int, int]) -> np.ndarray:
    seed = np.zeros((ROWS, COLS), dtype=np.uint8)
    for r, c in nodes:
        seed[r, c] = 1
    return seed


@pytest.fixture(scope="module")
def upstream() -> Any:
    from tin_engine._core import upstream as the_upstream

    return the_upstream


def expected_valley() -> np.ndarray:
    e = np.zeros((ROWS, COLS), dtype=np.uint8)
    e[RV, 0] = 1
    e[RV - 5 : RV + 6, 1 : COLS - 1] = 1
    return e


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_the_valley_through_the_binding(upstream: Any, dtype: type) -> None:
    out = upstream(view_of(valley(dtype)), seed_at((RV, 0)))
    assert isinstance(out.mask, np.ndarray)
    assert out.mask.dtype == np.uint8
    assert out.mask.shape == (ROWS, COLS)
    assert np.array_equal(out.mask, expected_valley())
    assert out.nodes_in == 11 * 13 + 1
    assert (out.row_min, out.row_max, out.col_min, out.col_max) == (5, 15, 0, 13)
    assert out.touches_edge is True
    assert out.touches_nodata is False


def test_the_mask_outlives_the_view(upstream: Any) -> None:
    view = view_of(valley())
    out = upstream(view, seed_at((RV, 0)))
    del view
    gc.collect()
    assert np.array_equal(out.mask, expected_valley())


def test_a_bool_seed_is_accepted_as_well(upstream: Any) -> None:
    out = upstream(view_of(valley()), seed_at((RV, 0)).astype(bool))
    assert np.array_equal(out.mask, expected_valley())


@pytest.mark.parametrize(
    "shape", [(ROWS, COLS - 1), (ROWS + 1, COLS), (COLS, ROWS), (ROWS * COLS,), (0, 0)]
)
def test_a_seed_of_another_shape_is_a_value_error(upstream: Any, shape: tuple[int, ...]) -> None:
    with pytest.raises(ValueError, match=r"(?i)shape|size"):
        upstream(view_of(valley()), np.zeros(shape, dtype=np.uint8))


def test_the_sentinel_is_nodata_and_never_in(upstream: Any) -> None:
    z = np.array(valley(), copy=True)
    z[RV, 5] = -9999.0
    z.flags.writeable = False
    out = upstream(view_of(z, nodata=-9999.0), seed_at((RV, 5), (RV, 6)))
    assert out.mask[RV, 5] == 0
    assert out.mask[RV, 6] == 1
    assert out.touches_nodata is True


def test_no_seed_no_catchment(upstream: Any) -> None:
    out = upstream(view_of(valley()), np.zeros((ROWS, COLS), dtype=np.uint8))
    assert out.nodes_in == 0
    assert not out.mask.any()
    assert out.touches_edge is False and out.touches_nodata is False
