"""Frozen edges through the binding: increment 23b, FE1, FE2 and FE5.

`docs/increments/23-basin-scale.md`, "The seam protocol" step 4 and K1, K2.
The C++ suites carry the oracles (`tests/cpp/property/prop_refinement_frozen.cpp`);
this file checks what crosses the boundary.

As N13 ("Settled after 23b's red step") rules it:

    _core.refine(..., *, tolerance, threads=0, min_angle_deg=0.0,
                 constraint_feet=False, frozen_mask=0)
    _core.refine_points(..., *, tolerance, threads=0, frozen_mask=0)
    _core.refine_strip(..., *, tolerance, threads=0, frozen_mask=0)   (N18)
    PointRefineOutcome.on_frozen            int
    PointRefineOutcome.on_frozen_max_error  float
    a negative frozen_mask is a TypeError (no uint32 conversion), as pybind11
    refuses any negative for an unsigned parameter.

Every new symbol is fetched inside a fixture, so a missing one fails its own
tests and leaves the rest of the session collecting.
"""

from __future__ import annotations

import ast
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from test_core_edge_strip import Start, world
from test_core_refine import DX, DY, X_MIN, Y_MAX, plane, start
from tin_engine import _core

STUB = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "_core.pyi"

Fn = Callable[..., Any]


@pytest.fixture
def refine() -> Fn:
    fn: Fn = _core.refine  # type: ignore[attr-defined]
    return fn


@pytest.fixture
def refine_points() -> Fn:
    fn: Fn = _core.refine_points  # type: ignore[attr-defined]
    return fn


def spiked_outline() -> tuple[Any, Any, Any, Any, tuple[float, float]]:
    """A plane over 9 x 13 nodes, the outline's vertices every 4th node, and a
    40 m spike on node (row 0, col 2): exactly on the top outline edge between
    (0, 0) and (0, 4), so only a split of that edge can repair it. The CLI gives
    the outline no property bit, so its edges get mask 1 here, standing in for
    the seam bit."""
    array = plane(9, 13)
    array[0, 2] += 40.0
    view, mesh, edges, masks = start(array, 4)
    masks = np.ones_like(np.asarray(masks))
    return view, mesh, edges, masks, (X_MIN + 2 * DX, Y_MAX - 0 * DY)


def has_vertex(out: Any, p: tuple[float, float]) -> bool:
    return bool((np.asarray(out.vertices) == np.array(p)).all(axis=1).any())


class TestRefineFrozenMask:
    def test_a_frozen_outline_is_never_split(self, refine: Fn) -> None:
        view, mesh, edges, masks, spike = spiked_outline()
        outline = 1
        control = refine(view, mesh, edges, masks, tolerance=0.5)
        assert control.ok(), control.message
        assert has_vertex(control, spike)

        out = refine(view, mesh, edges, masks, tolerance=0.5, frozen_mask=outline)
        assert out.ok(), out.message
        assert not has_vertex(out, spike)
        # The outline's edges are the start's, pair for pair, mask for mask.
        start_pairs = {tuple(sorted(e)) for e in np.asarray(edges).tolist()}
        out_pairs = {tuple(sorted(e)) for e in np.asarray(out.edges).tolist()}
        assert out_pairs == start_pairs
        assert_array_equal(np.sort(np.asarray(out.masks)), np.sort(np.asarray(masks)))
        assert out.max_error <= 0.5

    def test_a_mask_that_meets_no_edge_changes_nothing(self, refine: Fn) -> None:
        view, mesh, edges, masks, _ = spiked_outline()
        ref = refine(view, mesh, edges, masks, tolerance=0.5)
        out = refine(view, mesh, edges, masks, tolerance=0.5, frozen_mask=1 << 31)
        for name in ("vertices", "z", "valid", "triangles", "edges", "masks"):
            assert_array_equal(np.asarray(getattr(out, name)), np.asarray(getattr(ref, name)), name)
        for name in ("rounds", "inserted", "flips", "max_error"):
            assert getattr(out, name) == getattr(ref, name), name

    def test_frozen_mask_is_unsigned(self, refine: Fn) -> None:
        view, mesh, edges, masks, _ = spiked_outline()
        assert refine(view, mesh, edges, masks, tolerance=0.5, frozen_mask=0).ok()
        with pytest.raises(TypeError):
            refine(view, mesh, edges, masks, tolerance=0.5, frozen_mask=-1)


class TestRefinePointsFrozenMask:
    """The square [0, 8]^2 in lattice units (30 m cells), cut by the diagonal
    (0, 0)-(8, 8) with mask 4, z 0 at every vertex; one check point exactly on
    the diagonal at (4, 4) and one inside, both 100 m off."""

    H = 30.0

    def args(self) -> tuple[np.ndarray, ...]:
        h = self.H
        vertices = np.array([[0, 0], [8 * h, 0], [8 * h, -8 * h], [0, -8 * h]], dtype=np.float64)
        vertices += np.array([X_MIN, Y_MAX])
        # (col, row): 0 (0, 0), 1 (8, 0), 2 (8, 8), 3 (0, 8); counter-clockwise
        # in world, the diagonal 0-2 shared.
        triangles = np.array([[0, 3, 2], [0, 2, 1]], dtype=np.uint32)
        z = np.zeros(4)
        valid = np.ones(4, dtype=np.uint8)
        edges = np.array([[0, 1], [1, 2], [2, 3], [0, 3], [0, 2]], dtype=np.uint32)
        masks = np.array([1, 1, 1, 1, 4], dtype=np.uint32)
        return vertices, triangles, z, valid, edges, masks

    def store(self) -> Any:
        h = self.H
        cp = _core.CheckPoints(x_min=X_MIN, y_max=Y_MAX, spacing=h, rows=9, cols=9)  # type: ignore[attr-defined]
        xy = np.array([[X_MIN + 4 * h, Y_MAX - 4 * h], [X_MIN + 2 * h, Y_MAX - 5 * h]])
        cp.add(xy, np.array([100.0, 100.0], dtype=np.float32))
        cp.freeze()
        return cp

    def test_a_point_on_a_frozen_edge_is_reported_not_inserted(self, refine_points: Fn) -> None:
        on = (X_MIN + 4 * self.H, Y_MAX - 4 * self.H)
        control = refine_points(self.store(), *self.args(), tolerance=1.0)
        assert control.ok(), control.message
        assert has_vertex(control, on)
        assert control.on_frozen == 0

        out = refine_points(self.store(), *self.args(), tolerance=1.0, frozen_mask=4)
        assert out.ok(), out.message
        assert not has_vertex(out, on)
        assert has_vertex(out, (X_MIN + 2 * self.H, Y_MAX - 5 * self.H))
        assert out.on_frozen == 1
        assert out.on_frozen_max_error == 100.0
        assert isinstance(out.on_frozen, int)
        assert isinstance(out.on_frozen_max_error, float)


class TestRefineStripFrozenMask:
    """N18 ("If 15f-3 lands first"): `refine_strip` takes `frozen_mask`, as
    `refine` and `refine_points` do. The start is 15f-3's (refine's output on a
    ring of off-node vertices over rough terrain); its west side, one edge,
    gets the frozen bit F and the other three sides mask 1. The caller builds
    the strip on non-frozen edges only (N18); a strip on a frozen edge is
    refused (N16)."""

    F = 1 << 5

    @pytest.fixture(scope="class")
    def scene(self) -> tuple[Start, tuple[np.ndarray, ...], np.ndarray]:
        s = Start()
        v, t, z, ok, e, _ = s.args()
        west = np.array([self.on_west(v[a]) and self.on_west(v[b]) for a, b in e])
        assert west.sum() == 1, "the west side is one edge of the start"
        masks = np.where(west, self.F, 1).astype(np.uint32)
        return s, (v, t, z, ok, e, masks), west

    @staticmethod
    def on_west(p: np.ndarray) -> bool:
        a, b = world([Start.RING[0][0], Start.RING[1][0]], [Start.RING[0][1], Start.RING[1][1]])
        d = b - a
        return bool(abs((p - a) @ np.array([-d[1], d[0]])) / np.hypot(*d) < 1e-6)

    def added_on_west(self, out: Any, before: int) -> int:
        return sum(self.on_west(p) for p in np.asarray(out.vertices)[before:])

    def test_a_frozen_side_is_not_split(self, scene: Any) -> None:
        """The control is the whole strip with no mask: it splits the west side.

        This scene cannot see the DEM rescan's frozen skip: with the strip
        built on non-frozen edges and no mask, the run adds 91 vertices and
        none on the west side (measured at 472d91d), so the rescan never
        reaches it here. That skip stays pinned by refine's own frozen tests
        above, since refine_strip's rescan is refine's scan (N16).
        """
        s, args, west = scene
        v, _, _, _, e, masks = args
        everything = _core.constraint_check_points(s.view, v, e)  # type: ignore[attr-defined]
        control = _core.refine_strip(s.view, everything, *args, tolerance=0.5)  # type: ignore[attr-defined]
        assert control.ok(), control.message
        assert self.added_on_west(control, len(v)) > 0, "the control does not split the west side"

        unfrozen = e[(masks & self.F) == 0]
        strip = _core.constraint_check_points(s.view, v, unfrozen)  # type: ignore[attr-defined]
        out = _core.refine_strip(s.view, strip, *args, tolerance=0.5, frozen_mask=self.F)  # type: ignore[attr-defined]
        assert out.ok(), out.message
        assert self.added_on_west(out, len(v)) == 0
        frozen_pairs = {tuple(sorted(map(int, p))) for p in e[west]}
        out_edges, out_masks = np.asarray(out.edges), np.asarray(out.masks)
        kept = {tuple(sorted(map(int, p))) for p in out_edges[(out_masks & self.F) != 0]}
        assert kept == frozen_pairs  # the frozen side is the start's edge, unsplit

    def test_a_strip_on_a_frozen_edge_is_refused(self, scene: Any) -> None:
        """N16: the caller filters; a strip edge that is frozen is a programming error."""
        s, args, _ = scene
        v, _, _, _, e, _ = args
        everything = _core.constraint_check_points(s.view, v, e)  # type: ignore[attr-defined]
        with pytest.raises(RuntimeError, match=r"^refine_strip: "):
            _core.refine_strip(s.view, everything, *args, tolerance=0.5, frozen_mask=self.F)  # type: ignore[attr-defined]

    def test_frozen_mask_0_is_todays_result(self, scene: Any) -> None:
        s, args, _ = scene
        v, _, _, _, e, _ = args
        strip = _core.constraint_check_points(s.view, v, e)  # type: ignore[attr-defined]
        ref = _core.refine_strip(s.view, strip, *args, tolerance=0.5)  # type: ignore[attr-defined]
        out = _core.refine_strip(s.view, strip, *args, tolerance=0.5, frozen_mask=0)  # type: ignore[attr-defined]
        for name in ("vertices", "z", "valid", "triangles", "edges", "masks"):
            assert_array_equal(np.asarray(getattr(out, name)), np.asarray(getattr(ref, name)), name)
        for name in (
            "rounds",
            "inserted",
            "flips",
            "max_error",
            "strip_inserted",
            "nodes_inserted",
        ):
            assert getattr(out, name) == getattr(ref, name), name

    def test_frozen_mask_is_unsigned(self, scene: Any) -> None:
        s, args, _ = scene
        v, _, _, _, e, masks = args
        strip = _core.constraint_check_points(s.view, v, e[(masks & self.F) == 0])  # type: ignore[attr-defined]
        assert _core.refine_strip(s.view, strip, *args, tolerance=0.5, frozen_mask=0).ok()  # type: ignore[attr-defined]
        with pytest.raises(TypeError):
            _core.refine_strip(s.view, strip, *args, tolerance=0.5, frozen_mask=-1)  # type: ignore[attr-defined]


class TestTheStub:
    def test_the_stub_declares_frozen_mask_and_on_frozen(self) -> None:
        tree = ast.parse(STUB.read_text())
        funcs = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
        classes = {n.name: n for n in tree.body if isinstance(n, ast.ClassDef)}
        for name in ("refine", "refine_points", "refine_strip"):
            assert "frozen_mask" in [a.arg for a in funcs[name].args.kwonlyargs], name
        members = {
            n.name for n in classes["PointRefineOutcome"].body if isinstance(n, ast.FunctionDef)
        }
        assert {"on_frozen", "on_frozen_max_error"} <= members
