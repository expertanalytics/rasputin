"""The refinement binding: increment 14, R8 and T6/T8 through Python.

`docs/increments/14-adaptive-refinement.md`. `_core.refine(view, mesh, edges,
masks, *, tolerance, threads=0)` refines a start mesh against the DEM and
returns one result object. The C++ suites carry the tolerance oracle and the
conformity checks; this file checks what crosses the boundary.

The binding is `refine` in `bindings/core.cpp`, typed in
`src_python/tin_engine/_core.pyi`.

Every new symbol is fetched inside a fixture, so a missing one fails its own
tests and leaves the rest of the session collecting.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from tin_engine import _core
from tin_engine._core import ChainRole
from tin_engine.cli import DEFAULT_SNAP_SPACING, _constraint_arrays, _engine

X_MIN, Y_MAX, DX, DY = 500000.0, 7000000.0, 10.0, 5.0

Refine = Callable[..., Any]


@pytest.fixture
def refine() -> Refine:
    fn: Refine = _core.refine  # type: ignore[attr-defined]
    return fn


def grid_nodes(rows: int, cols: int, stride: int) -> tuple[np.ndarray, list[int]]:
    """Every `stride`-th node plus the last, by RasterGeometry::node's expression,
    and the counter-clockwise perimeter ring."""
    rs = [*range(0, rows, stride)] + ([rows - 1] if (rows - 1) % stride else [])
    cs = [*range(0, cols, stride)] + ([cols - 1] if (cols - 1) % stride else [])
    xy = np.array([(X_MIN + c * DX, Y_MAX - r * DY) for r in rs for c in cs], dtype=np.float64)
    nr, nc = len(rs), len(cs)
    ring = (
        [(nr - 1) * nc + j for j in range(nc)]
        + [i * nc + nc - 1 for i in range(nr - 2, -1, -1)]
        + [j for j in range(nc - 2, -1, -1)]
        + [i * nc for i in range(1, nr - 1)]
    )
    return xy, ring


def start(array: np.ndarray, stride: int) -> tuple[Any, Any, np.ndarray, np.ndarray]:
    """The view, the start mesh and its constraint arrays, as the CLI builds them."""
    rows, cols = array.shape
    xy, ring = grid_nodes(rows, cols, stride)
    run = _engine(xy, [(ring, ChainRole.Outer, 0)], True, DEFAULT_SNAP_SPACING)
    assert run.mesh is not None and run.noded is not None, run.message
    edges, masks = _constraint_arrays(run.mesh, run.noded)
    array.setflags(write=False)
    view = _core.raster_view(array, x_min=X_MIN, y_max=Y_MAX, delta_x=DX, delta_y=DY)
    return view, run.mesh, edges, masks


def rough(rows: int, cols: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).uniform(0.0, 100.0, (rows, cols)).astype(np.float32)


def plane(rows: int, cols: int) -> np.ndarray:
    r, c = np.indices((rows, cols))
    return np.asarray(3 * c - 2 * r + 7, dtype=np.float32)


class TestTheResult:
    def test_a_plane_needs_no_refinement(self, refine: Refine) -> None:
        view, mesh, edges, masks = start(plane(9, 13), 4)
        out = refine(view, mesh, edges, masks, tolerance=0.0)
        assert out.ok(), out.message
        assert out.inserted == 0
        assert out.max_error == 0.0
        assert_array_equal(np.asarray(out.triangles), np.asarray(mesh.triangles))
        assert_array_equal(np.asarray(out.vertices), np.asarray(mesh.vertices))

    def test_shapes_dtypes_and_the_tolerance(self, refine: Refine) -> None:
        array = rough(17, 15, 1)
        view, mesh, edges, masks = start(array, 8)
        out = refine(view, mesh, edges, masks, tolerance=2.0)
        assert out.ok(), out.message
        v, z, valid = np.asarray(out.vertices), np.asarray(out.z), np.asarray(out.valid)
        tris, e, m = np.asarray(out.triangles), np.asarray(out.edges), np.asarray(out.masks)
        assert v.dtype == np.float64 and v.ndim == 2 and v.shape[1] == 2
        assert z.dtype == np.float64 and z.shape == (len(v),)
        assert valid.dtype == np.bool_ and valid.shape == (len(v),)
        assert tris.ndim == 2 and tris.shape[1] == 3
        assert e.ndim == 2 and e.shape[1] == 2 and m.shape == (len(e),)
        assert out.inserted == len(v) - len(np.asarray(mesh.vertices)) > 0
        assert out.rounds >= 2
        assert 0.0 <= out.max_error <= 2.0
        assert out.uncovered == 0
        # Every vertex is a node, carrying that node's value.
        cols = (v[:, 0] - X_MIN) / DX
        rows = (Y_MAX - v[:, 1]) / DY
        assert (cols == np.round(cols)).all() and (rows == np.round(rows)).all()
        assert_array_equal(z, array[rows.astype(int), cols.astype(int)].astype(np.float64))
        assert valid.all()

    def test_tolerance_is_keyword_only(self, refine: Refine) -> None:
        view, mesh, edges, masks = start(plane(5, 5), 4)
        with pytest.raises(TypeError):
            refine(view, mesh, edges, masks, 1.0)


class TestDeterminism:
    """T6, once through the binding."""

    @pytest.mark.parametrize("threads", [2, 7, 0])
    def test_bit_identical_for_any_thread_count(self, refine: Refine, threads: int) -> None:
        view, mesh, edges, masks = start(rough(33, 33, 7), 16)
        ref = refine(view, mesh, edges, masks, tolerance=3.0, threads=1)
        out = refine(view, mesh, edges, masks, tolerance=3.0, threads=threads)
        for name in ("vertices", "z", "valid", "triangles", "edges", "masks"):
            assert_array_equal(np.asarray(getattr(out, name)), np.asarray(getattr(ref, name)), name)
        for name in ("rounds", "inserted", "flips", "max_error", "uncovered"):
            assert getattr(out, name) == getattr(ref, name), name


class TestRefusals:
    """T8: a status, not an exception."""

    @pytest.mark.parametrize("tolerance", [-1.0, math.nan, math.inf])
    def test_a_bad_tolerance(self, refine: Refine, tolerance: float) -> None:
        view, mesh, edges, masks = start(plane(5, 5), 4)
        out = refine(view, mesh, edges, masks, tolerance=tolerance)
        assert not out.ok()
        assert out.status == _core.RefineStatus.InvalidTolerance  # type: ignore[attr-defined]
        assert out.message


class TestOffNodeStart:
    """Increment 16 (R2): an off-node start ring is refined; outside the node
    rectangle is ``OutsideGrid``, which replaces ``OffLattice``."""

    @staticmethod
    def ring(array: np.ndarray, ring_xy: list[tuple[float, float]]) -> tuple[Any, Any, Any, Any]:
        run = _engine(
            np.array(ring_xy),
            [([*range(len(ring_xy))], ChainRole.Outer, 0)],
            True,
            DEFAULT_SNAP_SPACING,
        )
        assert run.mesh is not None and run.noded is not None, run.message
        edges, masks = _constraint_arrays(run.mesh, run.noded)
        array.setflags(write=False)
        view = _core.raster_view(array, x_min=X_MIN, y_max=Y_MAX, delta_x=DX, delta_y=DY)
        return view, run.mesh, edges, masks

    def test_an_off_node_ring_is_refined_and_kept_where_it_is(self, refine: Refine) -> None:
        ring = [
            (X_MIN + 12.3, Y_MAX - 73.3),
            (X_MIN + 137.7, Y_MAX - 72.9),
            (X_MIN + 136.1, Y_MAX - 6.7),
            (X_MIN + 13.9, Y_MAX - 7.1),
        ]
        view, mesh, edges, masks = self.ring(rough(17, 15, 3), ring)
        out = refine(view, mesh, edges, masks, tolerance=1.0)
        assert out.ok(), out.message
        assert out.inserted > 0
        start = np.asarray(mesh.vertices)
        assert_array_equal(np.asarray(out.vertices)[: len(start)], start)

    def test_a_vertex_outside_the_grid_is_outside_grid(self, refine: Refine) -> None:
        ring = [
            (X_MIN - 5.0, Y_MAX - 73.3),
            (X_MIN + 137.7, Y_MAX - 72.9),
            (X_MIN + 136.1, Y_MAX - 6.7),
        ]
        view, mesh, edges, masks = self.ring(rough(17, 15, 3), ring)
        out = refine(view, mesh, edges, masks, tolerance=1.0)
        assert not out.ok()
        assert out.status == _core.RefineStatus.OutsideGrid  # type: ignore[attr-defined]
        assert out.message
        assert not hasattr(_core.RefineStatus, "OffLattice")
