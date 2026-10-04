"""`_core.accumulate`, the accumulation's binding (increment 29, PR 1).

`docs/increments/29-nve-reference-catchments.md`, "Accumulation, from the same
flood (C++)" and "The red suites", PR 1, Python. The accumulation's behaviour
is the C++ suite's (`tests/cpp/unit/test_hydrology_accumulate.cpp`, the
invariant-critical one); this file pins what crosses the binding: shapes,
dtypes, ownership, and the oracle against `_core.upstream` on two DEMs.

Interface assumed (names chosen here, stated in the handback; they mirror the
C++ `AccumulateOutcome` as `upstream`'s binding mirrors `UpstreamOutcome`):

- `_core.accumulate(view: RasterView) -> AccumulateOutcome`. Releases the GIL.
- `AccumulateOutcome.count`: `uint32`, `(rows, cols)`; 0 on NoData, else the
  nodes draining through the node, itself included.
- `AccumulateOutcome.reach`: `uint8`, `(rows, cols)`; bit 0 the catchment
  touches the window's edge, bit 1 it touches NoData.
- `AccumulateOutcome.flow_to`: `uint8`, `(rows, cols)`; `3*(dr+1)+(dc+1)` of
  the neighbour the node drains to, 255 for an outlet and on NoData.
- All three owned by the outcome: they outlive the view and its array.

Not tested here: the 2^32-node refusal's mapping to `ValueError`; a raster
that size cannot be built in a test, and the C++ suite pins the refusal.
"""

from __future__ import annotations

import gc
from typing import Any

import numpy as np
import pytest

from tin_engine._core import raster_view

ROWS, COLS, RV, W = 21, 15, 10, 6
SENTINEL = -9999.0
OUTLET = 255
WEST = 3  # dr = 0, dc = -1


def valley(dtype: type = np.float32) -> np.ndarray:
    """The C++ suites' V-valley (as in test_core_upstream.py): falling west at
    0.37 a column, sides rising at 1 a row to ridges at rows 4 and 16, outer
    slopes at 2.9; column 0 a wall at 1000 except the outlet (10, 0) at -1."""
    r, c = np.indices((ROWS, COLS))
    off = np.abs(r - RV).astype(dtype)
    side = np.where(off <= W, off, dtype(W) - dtype(2.9) * (off - dtype(W)))
    z = (dtype(0.37) * c.astype(dtype) + side).astype(dtype)
    z[:, 0] = 1000
    z[RV, 0] = -1
    z.flags.writeable = False
    return z


def pitted() -> np.ndarray:
    """A 9 x 11 DEM with pits, flats, equal heights and a sentinel hole, drawn
    from five levels with a fixed seed so the oracle has ties to decide."""
    rng = np.random.default_rng(29_101)
    z = rng.integers(0, 5, size=(9, 11)).astype(np.float64)
    z[4, 5] = SENTINEL
    z[2, 7] = SENTINEL
    z.flags.writeable = False
    return z


def view_of(z: np.ndarray, nodata: float | None = None) -> Any:
    return raster_view(
        z, x_min=500_000.0, y_max=6_600_000.0, delta_x=10.0, delta_y=5.0, nodata=nodata
    )


@pytest.fixture(scope="module")
def core() -> Any:
    from tin_engine import _core

    # Fails here, by name, until the binding exists.
    assert hasattr(_core, "accumulate"), "tin_engine._core has no accumulate"
    return _core


def neighbour(r: int, c: int, code: int) -> tuple[int, int]:
    return r + code // 3 - 1, c + code % 3 - 1


def assert_oracle(core: Any, z: np.ndarray, nodata: float | None) -> None:
    """For every node with data: count == upstream's nodes_in, both bits equal
    upstream's flags, and the node flow_to names has this node in its mask."""
    view = view_of(z, nodata)
    out = core.accumulate(view)
    invalid = np.isnan(z) if nodata is None else (np.isnan(z) | (z == nodata))
    rows, cols = z.shape
    for r in range(rows):
        for c in range(cols):
            if invalid[r, c]:
                assert out.count[r, c] == 0
                assert out.reach[r, c] == 0
                assert out.flow_to[r, c] == OUTLET
                continue
            seed = np.zeros(z.shape, dtype=np.uint8)
            seed[r, c] = 1
            u = core.upstream(view, seed)
            node = f"node ({r}, {c})"
            assert int(out.count[r, c]) == u.nodes_in, node
            assert bool(out.reach[r, c] & 1) is u.touches_edge, node
            assert bool(out.reach[r, c] & 2) is u.touches_nodata, node
            code = int(out.flow_to[r, c])
            if code != OUTLET:
                nr, nc = neighbour(r, c, code)
                seed_n = np.zeros(z.shape, dtype=np.uint8)
                seed_n[nr, nc] = 1
                assert core.upstream(view, seed_n).mask[r, c] == 1, node


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_shapes_and_dtypes(core: Any, dtype: type) -> None:
    out = core.accumulate(view_of(valley(dtype)))
    for name, want in (("count", np.uint32), ("reach", np.uint8), ("flow_to", np.uint8)):
        array = getattr(out, name)
        assert isinstance(array, np.ndarray), name
        assert array.dtype == want, name
        assert array.shape == (ROWS, COLS), name


def test_the_valley_through_the_binding(core: Any) -> None:
    out = core.accumulate(view_of(valley()))
    assert out.count[RV, 0] == 11 * 13 + 1
    assert out.flow_to[RV, 0] == OUTLET
    assert out.reach[RV, 0] & 1
    assert np.all(out.flow_to[RV, 1 : COLS - 1] == WEST)


def test_the_arrays_outlive_the_view_and_its_array(core: Any) -> None:
    z = np.array(valley(), copy=True)
    view = view_of(z)
    out = core.accumulate(view)
    expected = (out.count.copy(), out.reach.copy(), out.flow_to.copy())
    del view, z
    gc.collect()
    # Fresh allocations of the same size, to reuse freed memory if the
    # outcome borrowed it.
    _ = [np.full((ROWS, COLS), 7.0, dtype=np.float32) for _ in range(8)]
    assert np.array_equal(out.count, expected[0])
    assert np.array_equal(out.reach, expected[1])
    assert np.array_equal(out.flow_to, expected[2])


def test_the_oracle_on_the_valley(core: Any) -> None:
    assert_oracle(core, valley(), None)


def test_the_oracle_on_pits_flats_and_nodata(core: Any) -> None:
    z = pitted()
    out = core.accumulate(view_of(z, SENTINEL))
    # The oracle must be able to disagree: the hole sets bit 1 somewhere and
    # some count exceeds one.
    assert (out.reach & 2).any()
    assert out.count.max() > 1
    assert_oracle(core, z, SENTINEL)


def test_nan_is_nodata_as_well(core: Any) -> None:
    z = np.array(valley(), copy=True)
    z[RV, 5] = np.nan
    out = core.accumulate(view_of(z))
    assert out.count[RV, 5] == 0
    assert out.reach[RV, 5] == 0
    assert out.flow_to[RV, 5] == OUTLET
    assert out.reach[RV, 6] & 2
