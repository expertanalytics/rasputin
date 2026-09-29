"""The fine outline: marching squares on a catchment's node mask (increment 22, PR 1).

`docs/increments/22-auto-catchment.md`, "The fine outline" and "The red
suites" (PR 1, `test_outline.py`). No `_core`: the tracer is numpy only.

Interface assumed (the design fixes the module, the function and the
guarantee; the representation below is chosen here and stated in the
handback):

- `tin_engine.outline.trace(mask) -> list[numpy.ndarray]`. `mask` is a 2-D
  array, bool or uint8, non-zero for an in-node. Each ring is an `(k, 2)`
  float array of `(row, col)` in the mask's own indices (so a vertex on the
  padding side of row 0 has row -0.5). Every vertex is the midpoint of a
  lattice edge between an in-node and an out-node: one coordinate is an
  integer and the other a half. A ring may or may not repeat its first vertex
  at the end; this suite accepts both and counts vertices without the repeat.
- In the world frame (x = col, y = -row: x east, y north, as the tracer's
  `x = x_min + col * dx`, `y = y_max - row * dy` gives), rings at even
  nesting depth (outer) are counter-clockwise and rings at odd depth (holes)
  clockwise.

THE AREA ORACLE is independent of any tracing: marching squares with the
(8, 4) saddle rule gives each square of four nodes a fixed in-area by its
corners (0 in: 0; 1: 1/8; 2 adjacent: 1/2; 2 diagonal, joined: 3/4; 3: 7/8;
4: 1). The signed sum of the rings' areas (outer positive, holes negative)
must equal the sum over every square of the padded mask.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
import shapely
from shapely.geometry import LinearRing, Point, Polygon

from importscan import first_party_imports

Ring = npt.NDArray[np.float64]
Trace = Callable[[Any], list[Any]]


@pytest.fixture(scope="module")
def trace() -> Trace:
    from tin_engine.outline import trace as the_trace

    return the_trace


def open_ring(ring: npt.ArrayLike) -> Ring:
    """The ring's vertices without a repeated closing vertex."""
    a = np.asarray(ring, dtype=np.float64)
    assert a.ndim == 2 and a.shape[1] == 2, a.shape
    if len(a) > 1 and np.array_equal(a[0], a[-1]):
        a = a[:-1]
    return a


def world(ring: Ring) -> list[tuple[float, float]]:
    """(row, col) -> (x, y) = (col, -row)."""
    return [(float(c), float(-r)) for r, c in ring]


def signed_area(ring: Ring) -> float:
    """Shoelace in the world frame: positive counter-clockwise."""
    xy = np.asarray(world(ring))
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def square_area(mask: npt.ArrayLike) -> float:
    """The oracle: the marching-squares in-area of `mask`, square by square."""
    m = np.pad(np.asarray(mask) != 0, 1).astype(int)
    a, b, c, d = m[:-1, :-1], m[:-1, 1:], m[1:, :-1], m[1:, 1:]  # tl, tr, bl, br
    n = a + b + c + d
    diagonal = (n == 2) & (a == d)
    area = np.select(
        [n == 1, (n == 2) & ~diagonal, diagonal, n == 3, n == 4], [1 / 8, 1 / 2, 3 / 4, 7 / 8, 1.0]
    )
    return float(area.sum())


def depth(ring: Ring, others: list[Ring]) -> int:
    """How many other rings enclose this one (tested at its first vertex)."""
    p = Point(world(ring)[0])
    return sum(Polygon(world(o)).contains(p) for o in others)


def check_guarantee(mask: npt.NDArray[np.uint8], rings: list[Ring]) -> None:
    """The design's guarantee, every clause, on one mask."""
    rows, cols = mask.shape
    for ring in rings:
        assert len(ring) >= 4
        frac = np.abs(ring - np.round(ring))
        halves = np.isclose(frac, 0.5)
        # One coordinate an integer, the other a half, on every vertex.
        assert np.all(halves.sum(axis=1) == 1), ring
        assert np.all(np.isclose(frac, 0.0) | halves)
        assert ring[:, 0].min() >= -0.5 and ring[:, 0].max() <= rows - 0.5
        assert ring[:, 1].min() >= -0.5 and ring[:, 1].max() <= cols - 0.5
        lr = LinearRing(world(ring))
        assert lr.is_valid and lr.is_simple, shapely.validation.explain_validity(Polygon(lr))
    # No two rings share a point.
    for a, b in itertools.combinations(rings, 2):
        assert LinearRing(world(a)).disjoint(LinearRing(world(b)))
    # Outer counter-clockwise, holes clockwise, by nesting depth.
    for i, ring in enumerate(rings):
        d = depth(ring, rings[:i] + rings[i + 1 :])
        assert (signed_area(ring) > 0) == (d % 2 == 0), (i, d, signed_area(ring))
    # Area: the signed sum is the square-by-square oracle.
    assert sum(map(signed_area, rings)) == pytest.approx(square_area(mask), abs=1e-12)
    # In-nodes inside an odd number of rings, out-nodes an even number, each at
    # least a quarter of the cell's diagonal from every ring.
    polygons = [Polygon(world(r)) for r in rings]
    boundaries = [LinearRing(world(r)) for r in rings]
    quarter = np.sqrt(2.0) / 4
    for r, c in itertools.product(range(rows), range(cols)):
        p = Point(float(c), float(-r))
        inside = sum(poly.contains(p) for poly in polygons)
        assert inside % 2 == (1 if mask[r, c] else 0), (r, c, inside)
        for b in boundaries:
            assert b.distance(p) >= quarter - 1e-12, (r, c)


def traced(trace: Trace, mask: npt.ArrayLike) -> list[Ring]:
    m = np.asarray(mask, dtype=np.uint8)
    rings = [open_ring(r) for r in trace(m)]
    check_guarantee(m, rings)
    return rings


def test_the_tracer_imports_no_core() -> None:
    import tin_engine.outline as outline

    assert not any(name.startswith("tin_engine._core") for name in first_party_imports(outline))


def test_an_empty_mask_has_no_ring(trace: Trace) -> None:
    assert list(trace(np.zeros((4, 5), dtype=np.uint8))) == []


def test_a_single_node_is_one_diamond_of_half_a_cell(trace: Trace) -> None:
    mask = np.zeros((5, 6), dtype=np.uint8)
    mask[2, 3] = 1
    (ring,) = traced(trace, mask)
    assert len(ring) == 4
    assert {tuple(v) for v in ring.tolist()} == {(1.5, 3.0), (2.0, 3.5), (2.5, 3.0), (2.0, 2.5)}
    assert signed_area(ring) == pytest.approx(0.5)


def test_a_two_by_two_block_is_an_octagon_of_three_and_a_half(trace: Trace) -> None:
    mask = np.zeros((5, 5), dtype=np.uint8)
    mask[1:3, 2:4] = 1
    (ring,) = traced(trace, mask)
    assert len(ring) == 8
    assert signed_area(ring) == pytest.approx(3.5)


def test_an_l_is_one_ring(trace: Trace) -> None:
    mask = np.zeros((6, 6), dtype=np.uint8)
    mask[1:5, 1] = 1
    mask[4, 1:4] = 1
    (ring,) = traced(trace, mask)
    assert signed_area(ring) > 0


def test_a_diagonal_pair_is_one_ring_by_the_saddle_rule(trace: Trace) -> None:
    mask = np.zeros((4, 4), dtype=np.uint8)
    mask[1, 1] = mask[2, 2] = 1
    (ring,) = traced(trace, mask)
    assert signed_area(ring) == pytest.approx(1 / 8 * 6 + 3 / 4)


def test_the_anti_diagonal_pair_is_one_ring_too(trace: Trace) -> None:
    mask = np.zeros((4, 4), dtype=np.uint8)
    mask[1, 2] = mask[2, 1] = 1
    assert len(traced(trace, mask)) == 1


def test_a_pinch_is_one_simple_ring(trace: Trace) -> None:
    # Two 2 x 2 blocks meeting at one diagonal.
    mask = np.zeros((6, 6), dtype=np.uint8)
    mask[1:3, 1:3] = 1
    mask[3:5, 3:5] = 1
    assert len(traced(trace, mask)) == 1


def test_a_ring_of_nodes_is_an_outer_ring_and_a_clockwise_hole(trace: Trace) -> None:
    mask = np.zeros((5, 5), dtype=np.uint8)
    mask[1:4, 1:4] = 1
    mask[2, 2] = 0
    rings = traced(trace, mask)
    assert len(rings) == 2
    outer, hole = sorted(rings, key=signed_area, reverse=True)
    assert signed_area(outer) > 0
    assert signed_area(hole) == pytest.approx(-0.5)  # the diamond around the out-node


def test_a_diagonal_ring_holds_a_hole_the_out_node_cannot_leave(trace: Trace) -> None:
    # A plus of four in-nodes round an out-node: the four are 8-connected, so
    # the centre is enclosed (the (8, 4) pairing), a clockwise hole of half a cell.
    mask = np.zeros((5, 5), dtype=np.uint8)
    for r, c in ((1, 2), (2, 1), (2, 3), (3, 2)):
        mask[r, c] = 1
    rings = traced(trace, mask)
    assert len(rings) == 2
    assert sorted(map(signed_area, rings)) == pytest.approx([-0.5, 4.5])


def test_disjoint_pieces_are_disjoint_rings(trace: Trace) -> None:
    mask = np.zeros((7, 9), dtype=np.uint8)
    mask[1:3, 1:3] = 1
    mask[4, 6] = 1
    mask[1, 7] = 1
    assert len(traced(trace, mask)) == 3


@pytest.mark.parametrize("dtype", [np.uint8, np.bool_])
def test_a_mask_filling_the_window_closes_through_the_padding(trace: Trace, dtype: Any) -> None:
    mask = np.ones((3, 4), dtype=dtype)
    rings = [open_ring(r) for r in trace(mask)]
    check_guarantee(mask.astype(np.uint8), rings)
    (ring,) = rings
    assert ring[:, 0].min() == -0.5 and ring[:, 0].max() == 2.5
    assert ring[:, 1].min() == -0.5 and ring[:, 1].max() == 3.5


def test_a_single_row_and_a_single_column(trace: Trace) -> None:
    traced(trace, np.ones((1, 5)))
    traced(trace, np.ones((5, 1)))
    traced(trace, np.ones((1, 1)))


@pytest.mark.parametrize("seed", range(40))
def test_random_masks_meet_the_guarantee(trace: Trace, seed: int) -> None:
    rng = np.random.default_rng(22_100 + seed)
    rows, cols = (int(v) for v in rng.integers(1, 13, size=2))
    density = [0.2, 0.5, 0.8][seed % 3]
    mask = (rng.random((rows, cols)) < density).astype(np.uint8)
    traced(trace, mask)


def test_a_checkerboard_is_one_ring_with_holes(trace: Trace) -> None:
    # Every in-node touches its diagonal neighbours: one 8-connected piece.
    r, c = np.indices((7, 7))
    mask = ((r + c) % 2 == 0).astype(np.uint8)
    rings = traced(trace, mask)
    outer = [ring for ring in rings if signed_area(ring) > 0]
    assert len(outer) == 1
