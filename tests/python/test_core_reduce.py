"""`_core.reduce_ring`, the outline reduction through the binding (increment 22, PR 2).

`docs/increments/22-auto-catchment.md`, "The reduction" and "The red suites"
(PR 2, `test_core_reduce.py`): on rings traced from random blobs and from the
PR 1 synthetic catchments, the five guarantees and determinism.

Interface assumed (names chosen here, stated in the handback; they mirror the
C++ `ReduceOutcome` in `tests/cpp/unit/test_area_collapse.cpp`):

- `_core.reduce_ring(ring, tolerance, keep) -> ReduceOutcome`. `ring` is an
  `(N, 2)` float64 array-like of `(x, y)`, open (the first vertex not
  repeated), counter-clockwise; `keep` an `(K, 2)` array-like of points that
  must stay strictly inside. Any other shape is a `ValueError`. Releases the
  GIL.
- `ReduceOutcome.ring`: an `(M, 2)` float64 array, open. `.status`: a
  `_core.ReduceStatus` (`Ok`, `InvalidTolerance`, `NotCounterClockwise`,
  `TooFewVertices`). Counts: `collinear`, `collapses`, `rejected_crossing`,
  `rejected_seed`, `rejected_tolerance`.

The guarantees, as the design words its tests: area `|A_r - A_f| <= 1e-9
A_f`; simple (shapely `is_valid`, `is_simple`); the seed strictly inside
(shapely `contains`); every fine vertex within the tolerance of the reduced
ring and every reduced vertex within it of the fine ring; the Hausdorff
distance (`densify=0.05`) at most the tolerance plus one cell.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
import shapely
from shapely.geometry import LinearRing, Point, Polygon

from catchment_fixtures import EPSG, MemoryRepository, bowl, lake_box, lat, tile_of
from mosaic_fixtures import quadrants
from test_catchment import plus_with_a_hole

CELL = 10.0


@pytest.fixture(scope="module")
def core() -> Any:
    import tin_engine._core as core

    assert hasattr(core, "reduce_ring"), "tin_engine._core has no reduce_ring"
    return core


def open_ring(coords: Any) -> np.ndarray:
    a = np.asarray(coords, dtype=np.float64)
    return a[:-1] if len(a) > 1 and np.array_equal(a[0], a[-1]) else a


def check(
    core: Any, fine: np.ndarray, tolerance: float, keep: np.ndarray, cell: float = CELL
) -> np.ndarray:
    """Reduce `fine` and check every guarantee; the reduced ring."""
    out = core.reduce_ring(fine, tolerance, keep)
    assert out.status == core.ReduceStatus.Ok
    red = np.asarray(out.ring, dtype=np.float64)
    assert red.ndim == 2 and red.shape[1] == 2
    assert np.isfinite(red).all()
    assert len(fine) - out.collinear - out.collapses == len(red)
    f, r = Polygon(fine), Polygon(red)
    assert r.exterior.is_ccw
    assert abs(r.area - f.area) <= 1e-9 * f.area, (r.area, f.area)
    assert LinearRing(red).is_simple and r.is_valid
    for k in keep:
        assert r.contains(Point(k))
    slack = 1e-9 * max(1.0, tolerance)
    assert shapely.distance(shapely.points(fine), r.exterior).max() <= tolerance + slack
    assert shapely.distance(shapely.points(red), f.exterior).max() <= tolerance + slack
    assert shapely.hausdorff_distance(f.exterior, r.exterior, densify=0.05) <= tolerance + cell
    # Deterministic: the same input again gives the same bits.
    again = core.reduce_ring(fine, tolerance, keep)
    assert np.asarray(again.ring).tobytes() == np.asarray(out.ring).tobytes()
    return red


def blob_ring(seed: int) -> tuple[np.ndarray, np.ndarray]:
    """The outer ring round a random blob of 10 m nodes (a union of disks),
    in metres from its window's lower-left corner, and one in-node inside it."""
    from tin_engine.outline import trace

    rng = np.random.default_rng(22_300 + seed)
    n = 48
    r, c = np.indices((n, n))
    mask = np.zeros((n, n), dtype=bool)
    for _ in range(int(rng.integers(3, 9))):
        cr, cc = rng.integers(12, 36, size=2)
        mask |= np.hypot(r - cr, c - cc) <= rng.uniform(2.0, 10.0)
    mask[24, 24] = True
    rings = [open_ring(ring) for ring in trace(mask)]
    world = [np.column_stack([ring[:, 1] * CELL, (n - 1 - ring[:, 0]) * CELL]) for ring in rings]
    outer = max(world, key=lambda w: Polygon(w).area if LinearRing(w).is_ccw else -1.0)
    keep = np.array([[24 * CELL, (n - 1 - 24) * CELL]])
    if not Polygon(outer).contains(Point(keep[0])):
        keep = np.empty((0, 2))
    return outer, keep


@pytest.mark.parametrize("tolerance", [0.0, 5.0, 10.0, 20.0, 50.0])
@pytest.mark.parametrize("seed", range(12))
def test_random_blobs_keep_every_guarantee(core: Any, seed: int, tolerance: float) -> None:
    fine, keep = blob_ring(seed)
    red = check(core, fine, tolerance, keep)
    if tolerance == 0.0:
        # Only the collinear pass: every reduced vertex is a fine vertex.
        assert {tuple(p) for p in red.tolist()} <= {tuple(p) for p in fine.tolist()}


def test_blobs_are_reduced_substantially_at_two_cells(core: Any) -> None:
    kept = []
    for seed in range(12):
        fine, keep = blob_ring(seed)
        zero = core.reduce_ring(fine, 0.0, keep)
        two = core.reduce_ring(fine, 2 * CELL, keep)
        kept.append(len(two.ring) / len(zero.ring))
    # A reducer that stopped after the collinear pass would keep all of them.
    assert np.median(kept) < 0.8


def catchment_ring(fine: Polygon) -> tuple[np.ndarray, np.ndarray]:
    """A PR 1 fine outline shifted to its lower-left corner, as the design has
    Python do before the reduction."""
    x0, y0 = fine.bounds[0], fine.bounds[1]
    return open_ring(np.asarray(fine.exterior.coords)) - (x0, y0), np.array([x0, y0])


@pytest.mark.parametrize("tolerance", [0.0, 100.0, 200.0, 500.0])
def test_the_pr1_bowl_catchment_keeps_every_guarantee(core: Any, tolerance: float) -> None:
    from tin_engine.catchment import CatchmentRequest, delineate

    z = bowl()
    result = delineate(
        CatchmentRequest(seed=lat(100, 200), seed_crs=EPSG, lakes=(lake_box(),), lakes_crs=EPSG),
        MemoryRepository(quadrants(tile_of(z), row_cut=200, col_cut=100, overlap=1)),
    )
    fine, origin = catchment_ring(result.fine)
    seed = np.array([lat(100, 200)]) - origin
    check(core, fine, tolerance, seed, cell=100.0)


def test_the_pr1_plus_catchment_keeps_every_guarantee(core: Any) -> None:
    from tin_engine.catchment import CatchmentRequest, delineate

    z, lake = plus_with_a_hole()
    result = delineate(
        CatchmentRequest(seed=lat(35, 34), seed_crs=EPSG, lakes=(lake,), lakes_crs=EPSG),
        MemoryRepository({"t.tif": tile_of(z)}),
    )
    fine, origin = catchment_ring(result.fine)
    check(core, fine, 200.0, np.array([lat(35, 34)]) - origin, cell=100.0)


def test_the_statuses_cross_the_binding(core: Any) -> None:
    square = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]])
    none = np.empty((0, 2))
    status = core.ReduceStatus
    assert core.reduce_ring(square, 1.0, none).status == status.Ok
    assert core.reduce_ring(square, -1.0, none).status == status.InvalidTolerance
    assert core.reduce_ring(square, float("nan"), none).status == status.InvalidTolerance
    assert core.reduce_ring(square, float("inf"), none).status == status.InvalidTolerance
    assert core.reduce_ring(square[::-1].copy(), 1.0, none).status == status.NotCounterClockwise
    assert core.reduce_ring(square[:3], 1.0, none).status == status.TooFewVertices


@pytest.mark.parametrize("shape", [(4,), (4, 3), (2, 4, 2)])
def test_a_ring_of_another_shape_is_a_value_error(core: Any, shape: tuple[int, ...]) -> None:
    with pytest.raises(ValueError, match=r"(?i)shape|\(N, ?2\)"):
        core.reduce_ring(np.zeros(shape), 1.0, np.empty((0, 2)))


def test_a_keep_of_another_shape_is_a_value_error(core: Any) -> None:
    square = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]])
    with pytest.raises(ValueError, match=r"(?i)shape|\(K, ?2\)|\(N, ?2\)"):
        core.reduce_ring(square, 1.0, np.zeros((3,)))
