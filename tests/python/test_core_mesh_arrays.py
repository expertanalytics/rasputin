"""A start mesh from arrays, through the binding: what 23c's start slice needs.

`docs/increments/23-basin-scale.md`, "The seam protocol" step 3: a piece's
start mesh is a NumPy slice of the one start triangulation, with each seam
edge split into a fan, and `refine` is then run on it (step 4). Today
`refine` takes an `IndexedMesh2`, which is "Not constructible from Python"
(`_core.pyi`), and the only producer is `triangulate`. The design's PR table
for 23c has no binding row, so the slice cannot reach `refine` as designed.

PINNED HERE, for `@architect` to confirm or replace (handback, choice 1):

    _core.indexed_mesh(vertices, triangles, constrained_edges) -> IndexedMesh2

`vertices` `(N, 2)` float64, `triangles` `(T, 3)` integer indices,
`constrained_edges` `(T,)` with bit `e` set iff edge `(v[e], v[(e + 1) % 3])`
is constrained (the mesh's own convention). It copies; a shape, an index out
of range, a mask above 7 or a non-finite coordinate is a `ValueError`. It
checks no orientation: `refine` already refuses a clockwise triangle
(`NotCounterClockwise`), and that refusal is what the degeneracy policy's
"a fan that is not counter-clockwise" relies on, so it is pinned here too.

HOW THIS FILE GOES RED: `_core.indexed_mesh` does not exist; the `indexed_mesh`
fixture fails every test on `AttributeError`.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from test_core_refine import rough, start
from tin_engine import _core

Fn = Callable[..., Any]


@pytest.fixture
def indexed_mesh() -> Fn:
    fn: Fn = _core.indexed_mesh  # type: ignore[attr-defined]
    return fn


def arrays(mesh: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.array(mesh.vertices, dtype=np.float64),
        np.array(mesh.triangles, dtype=np.uint32),
        np.array(mesh.constrained_edges, dtype=np.uint8),
    )


class TestRoundTrip:
    def test_the_arrays_come_back_as_given(self, indexed_mesh: Fn) -> None:
        _, mesh, _, _ = start(rough(17, 15, 3), 4)
        v, t, c = arrays(mesh)
        rebuilt = indexed_mesh(v, t, c)
        assert isinstance(rebuilt, _core.IndexedMesh2)
        assert_array_equal(np.asarray(rebuilt.vertices), v)
        assert_array_equal(np.asarray(rebuilt.triangles), t)
        assert_array_equal(np.asarray(rebuilt.constrained_edges), c)
        assert rebuilt.triangle_count == len(t)

    def test_it_copies(self, indexed_mesh: Fn) -> None:
        _, mesh, _, _ = start(rough(9, 9, 4), 4)
        v, t, c = arrays(mesh)
        rebuilt = indexed_mesh(v, t, c)
        v[:] = 0.0
        t[:] = 0
        assert not (np.asarray(rebuilt.vertices) == 0.0).all()
        assert np.asarray(rebuilt.triangles).any()

    @pytest.mark.parametrize("tolerance", [0.0, 2.0, 10.0])
    def test_refine_on_the_rebuilt_mesh_is_refine_on_the_original(
        self, indexed_mesh: Fn, tolerance: float
    ) -> None:
        """K1 for the slice path: a whole-domain slice refines bit for bit as
        the triangulation it was cut from."""
        view, mesh, edges, masks = start(rough(17, 15, 5), 4)
        ref = _core.refine(view, mesh, edges, masks, tolerance=tolerance)
        out = _core.refine(view, indexed_mesh(*arrays(mesh)), edges, masks, tolerance=tolerance)
        assert ref.ok() and out.ok(), (ref.message, out.message)
        for name in ("vertices", "z", "valid", "triangles", "edges", "masks"):
            assert_array_equal(np.asarray(getattr(out, name)), np.asarray(getattr(ref, name)), name)
        for name in ("rounds", "inserted", "flips", "max_error"):
            assert getattr(out, name) == getattr(ref, name), name


class TestRefusals:
    def good(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        v = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0], [10.0, 10.0]])
        t = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.uint32)
        c = np.array([0b101, 0b011], dtype=np.uint8)
        return v, t, c

    def test_the_good_case_is_accepted(self, indexed_mesh: Fn) -> None:
        assert indexed_mesh(*self.good()).triangle_count == 2

    @pytest.mark.parametrize(
        "case",
        [
            "vertices (N, 3)",
            "triangles (T, 2)",
            "index out of range",
            "masks of another length",
            "mask above 7",
            "NaN vertex",
            "infinite vertex",
        ],
    )
    def test_refused(self, indexed_mesh: Fn, case: str) -> None:
        v, t, c = self.good()
        if case == "vertices (N, 3)":
            v = np.column_stack([v, np.zeros(len(v))])
        elif case == "triangles (T, 2)":
            t = t[:, :2]
        elif case == "index out of range":
            t[1, 1] = 4
        elif case == "masks of another length":
            c = c[:1]
        elif case == "mask above 7":
            c[0] = 8
        elif case == "NaN vertex":
            v[3, 0] = np.nan
        else:
            v[3, 1] = np.inf
        with pytest.raises(ValueError):
            indexed_mesh(v, t, c)

    def test_a_clockwise_triangle_is_refines_refusal(self, indexed_mesh: Fn) -> None:
        """Not the constructor's: `refine` names it, as for a bad fan."""
        array = rough(3, 3, 6)
        array.setflags(write=False)
        view = _core.raster_view(array, x_min=0.0, y_max=20.0, delta_x=10.0, delta_y=10.0)
        v = np.array([[0.0, 0.0], [20.0, 0.0], [0.0, 20.0], [20.0, 20.0]])
        t = np.array([[0, 2, 1], [1, 2, 3]], dtype=np.uint32)  # the first is clockwise
        c = np.array([0b101, 0b110], dtype=np.uint8)
        edges = np.array([[0, 1], [1, 3], [3, 2], [2, 0]], dtype=np.uint32)
        masks = np.zeros(4, dtype=np.uint32)
        out = _core.refine(view, indexed_mesh(v, t, c), edges, masks, tolerance=1.0)
        assert not out.ok()
        assert out.status == _core.RefineStatus.NotCounterClockwise, out.status
