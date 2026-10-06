"""The edge strip through the binding: increment 15f-3, PY1.

`docs/increments/15f-edge-strip.md`, D2-D5 (`bindings/core.cpp`'s row), E1,
E2, E4, E7, L2, L6, L8 and "Tests for @tester" (PY1). Every binding of the
increment: the store `ConstraintCheckPoints` (read-only `size`, `no_data`,
`duplicates`, `edge_count`), the generator `constraint_check_points(view,
vertices, edges)`, `refine_points(..., strip=None)`, `refine_strip(view,
strip, vertices, triangles, z, valid, edges, masks, *, tolerance, threads=0)`
and `PointRefineOutcome`'s new fields. Each releases the GIL.

CHOSEN HERE, where the design is silent:

- `strip` is keyword-only on `refine_points`, after `threads`: D5 writes it
  after the existing parameters, whose `tolerance` and `threads` already are
  keyword-only, so it cannot be positional there.
- Every refusal of `constraint_check_points`, the binding's shape checks and
  the C++ `std::invalid_argument` alike, is a `ValueError` whose message
  contains `constraint_check_points` (15f-1's settled message, D2); a NaN
  vertex is outside the node rectangle. `refine_strip`'s shape refusals are a
  `ValueError` naming `refine_strip`, as `refine_points`' name it.
- The store is not constructible from Python (only the generator builds it,
  D3), a `TypeError`, and its four properties are read-only.

The property cases carry `tester.md` §3D's oracles from `strip_oracle.py` as
L8 assigns them: on `refine_strip` (the projected path) the strip oracle (E1),
the DEM-node tolerance oracle (E2) and the constrained-Delaunay oracle; on
`refine_points(..., strip)` (the reprojected path) the strip oracle, J2 by
15c's RP3 oracle, and the Delaunay oracle. The start is refine's output on a
ring with off-node vertices, so its constraint edges cross grid lines.

Went red at `4157dab` because none of these names was bound, so every
fixture that fetched one failed with `AttributeError`, and `refine_points`
refused the `strip` keyword with `TypeError`. The oracle self-checks at the
end passed then as now.

Not invariant-critical (ES2, ES3, ES5 are, in C++); no mutation round.
"""

from __future__ import annotations

import ast
from collections.abc import Callable
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from gil_probe import ticks_during
from strip_oracle import (
    Grid,
    delaunay_violations,
    edge_segments,
    exact_vertices,
    node_findings,
    ruled_points,
    strip_findings,
)
from test_core_cdt import RELEASED_TICKS
from test_core_refine_points import scattered, worst_excess
from tin_engine import _core
from tin_engine._core import ChainRole
from tin_engine.cli import DEFAULT_SNAP_SPACING, _constraint_arrays, _engine

X_MIN, Y_MAX, H = 500000.0, 7000000.0, 30.0
STUB = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "_core.pyi"
#: A written vertex on a constraint is off its line by world rounding only.
ON_EDGE = 1e-6

Factory = Callable[..., Any]


@pytest.fixture
def generator() -> Factory:
    fn: Factory = _core.constraint_check_points  # type: ignore[attr-defined]
    return fn


@pytest.fixture
def refine_strip() -> Factory:
    fn: Factory = _core.refine_strip  # type: ignore[attr-defined]
    return fn


@pytest.fixture
def store_class() -> Any:
    return _core.ConstraintCheckPoints  # type: ignore[attr-defined]


def world(col: Any, row: Any) -> np.ndarray:
    return np.column_stack([X_MIN + np.asarray(col) * H, Y_MAX - np.asarray(row) * H])


def view_of(array: np.ndarray, **kwargs: Any) -> Any:
    return _core.raster_view(array, x_min=X_MIN, y_max=Y_MAX, delta_x=H, delta_y=H, **kwargs)


def rough(n: int, seed: int) -> np.ndarray:
    """Relief at about two cells' wavelength plus noise, float32."""
    r, c = np.indices((n, n), dtype=np.float64)
    noise = np.random.default_rng(seed).normal(0.0, 1.0, (n, n))
    z = 100 + 12 * np.sin(c / 2.3) * np.cos(r / 1.9) + 6 * np.sin((r + 2 * c) / 1.7) + noise
    return z.astype(np.float32)


class Start:
    """refine's output from a ring of off-node vertices on a rough DEM: the
    start of a strip run, as the CLI has it after `refine`."""

    #: (col, row), counter-clockwise in world (y up): down the west side first.
    RING: ClassVar = [(2.3, 2.6), (3.1, 21.7), (22.6, 22.3), (21.7, 2.2)]

    def __init__(self, n: int = 25, tolerance: float = 2.0, seed: int = 3) -> None:
        self.n = n
        self.dem = rough(n, seed)
        self.dem.setflags(write=False)
        self.view = view_of(self.dem)
        self.grid = Grid(X_MIN, Y_MAX, H, H, self.dem.astype(np.float64))
        xy = world([c for c, _ in self.RING], [r for _, r in self.RING])
        run = _engine(xy, [([0, 1, 2, 3], ChainRole.Outer, 0)], True, DEFAULT_SNAP_SPACING)
        assert run.mesh is not None and run.noded is not None, run.message
        edges, masks = _constraint_arrays(run.mesh, run.noded)
        out = _core.refine(self.view, run.mesh, edges, masks, tolerance=tolerance)
        assert out.ok(), out.message
        self.out = out

    def args(self) -> tuple[np.ndarray, ...]:
        o = self.out
        return tuple(
            np.asarray(a) for a in (o.vertices, o.triangles, o.z, o.valid, o.edges, o.masks)
        )

    def strip(self, generator: Factory) -> Any:
        return generator(self.view, np.asarray(self.out.vertices), np.asarray(self.out.edges))

    def oracle_points(self) -> tuple[np.ndarray, np.ndarray]:
        """E1's points: crossings and midpoints of the start's constraint edges."""
        return ruled_points(edge_segments(self.args()[0], self.args()[4]), self.grid)

    def store(self, xy: np.ndarray, z: np.ndarray) -> Any:
        cp = _core.CheckPoints(x_min=X_MIN, y_max=Y_MAX, spacing=H, rows=self.n, cols=self.n)
        cp.add(xy, z)
        cp.freeze()
        return cp


def strip_check(start: Start, out: Any, tolerance: float) -> Any:
    xy, z = start.oracle_points()
    return strip_findings(
        xy,
        z,
        np.asarray(out.vertices),
        np.asarray(out.z),
        np.asarray(out.edges),
        tolerance=tolerance,
        on_edge=ON_EDGE,
        slack=0.0,
    )


def delaunay_check(start: Start, out: Any) -> int:
    v = np.asarray(out.vertices)
    exact = exact_vertices(start.grid, v, start.args()[0])
    return delaunay_violations(
        start.grid, v, np.asarray(out.triangles), np.asarray(out.edges), exact
    )


@pytest.fixture(scope="module")
def start() -> Start:
    return Start()


# ------------------------------------------------------------------ the generator and its store


#: One edge, (0.5, 0.5) to (3.5, 1.5) in (col, row): column crossings at
#: K = 1, 2, 3 (t = 1/6, 1/2, 5/6), and the row crossing R = 1 at t = 1/2 is
#: the node (2, 1), merged with the column crossing there (D2, steps 4-5).
#: So 3 crossings kept, 2c + 1 = 7 points, 1 duplicate.
HAND = world([0.5, 3.5], [0.5, 1.5])


def flat(n: int = 5, dtype: Any = np.float32) -> np.ndarray:
    r, c = np.indices((n, n))
    return (3.0 * c + 2.0 * r).astype(dtype)


class TestTheGenerator:
    @pytest.mark.parametrize("dtype", [np.float32, np.float64], ids=["float32", "float64"])
    def test_counts_on_one_edge_by_hand(self, generator: Factory, dtype: Any) -> None:
        strip = generator(view_of(flat(dtype=dtype)), HAND, np.array([[0, 1]], np.uint32))
        assert (strip.size, strip.duplicates, strip.no_data, strip.edge_count) == (7, 1, 0, 1)

    def test_the_counts_do_not_depend_on_the_pairs_order(self, generator: Factory) -> None:
        strip = generator(view_of(flat()), HAND, np.array([[1, 0]], np.uint32))
        assert (strip.size, strip.duplicates, strip.edge_count) == (7, 1, 1)

    def test_two_vertices_at_one_position_give_no_point(self, generator: Factory) -> None:
        """CC8 at the binding: the one midpoint candidate sits on both ends."""
        xy = np.vstack([HAND, HAND[:1]])
        strip = generator(view_of(flat()), xy, np.array([[0, 1], [0, 2]], np.uint32))
        assert (strip.size, strip.duplicates, strip.edge_count) == (7, 2, 2)

    def test_no_edges_gives_an_empty_store(self, generator: Factory) -> None:
        strip = generator(view_of(flat()), HAND, np.zeros((0, 2), np.uint32))
        assert (strip.size, strip.duplicates, strip.no_data, strip.edge_count) == (0, 0, 0, 0)

    @pytest.mark.parametrize("void", ["nan", "sentinel"])
    def test_points_beside_nodata_are_dropped_and_counted(
        self, generator: Factory, void: str
    ) -> None:
        """CC5 at the binding: node (row 1, col 3) is a corner of the cell that
        holds the crossing at K = 3; NoData by NaN or by the sentinel."""
        array = flat()
        array[1, 3] = np.nan if void == "nan" else -9999.0
        kwargs = {"nodata": -9999.0} if void == "sentinel" else {}
        strip = generator(view_of(array, **kwargs), HAND, np.array([[0, 1]], np.uint32))
        assert strip.no_data >= 1
        assert strip.size + strip.no_data == 7
        assert strip.duplicates == 1

    def test_the_store_is_read_only_and_built_only_by_the_generator(
        self, generator: Factory, store_class: Any
    ) -> None:
        with pytest.raises(TypeError):
            store_class()
        strip = generator(view_of(flat()), HAND, np.array([[0, 1]], np.uint32))
        for name in ("size", "no_data", "duplicates", "edge_count"):
            with pytest.raises(AttributeError):
                setattr(strip, name, 0)


class TestGeneratorRefusals:
    @pytest.mark.parametrize(
        ("vertices", "edges"),
        [
            (np.zeros((2, 3)), np.array([[0, 1]], np.uint32)),
            (np.zeros(4), np.array([[0, 1]], np.uint32)),
            (HAND, np.array([[0, 1, 1]], np.uint32)),
            (HAND, np.array([0, 1], np.uint32)),
            (np.array([["a", "b"], ["c", "d"]]), np.array([[0, 1]], np.uint32)),
            (HAND, np.array([["a", "b"]])),
        ],
        ids=["xy-3-columns", "xy-1d", "edges-3-columns", "edges-1d", "xy-text", "edges-text"],
    )
    def test_shapes_and_dtypes(self, generator: Factory, vertices: Any, edges: Any) -> None:
        with pytest.raises(ValueError, match="constraint_check_points"):
            generator(view_of(flat()), vertices, edges)

    @pytest.mark.parametrize(
        ("vertices", "edges"),
        [
            (HAND, np.array([[0, 2]], np.uint32)),
            (HAND, np.array([[1, 1]], np.uint32)),
            (world([0.5, 4.5], [0.5, 1.5]), np.array([[0, 1]], np.uint32)),
            (np.array([[np.nan, Y_MAX - 15.0], HAND[1]]), np.array([[0, 1]], np.uint32)),
            (np.array([[np.inf, Y_MAX - 15.0], HAND[1]]), np.array([[0, 1]], np.uint32)),
        ],
        ids=[
            "index-out-of-range",
            "degenerate-edge",
            "outside-the-rectangle",
            "nan-end",
            "inf-end",
        ],
    )
    def test_refused_inputs(self, generator: Factory, vertices: Any, edges: Any) -> None:
        with pytest.raises(ValueError, match="constraint_check_points"):
            generator(view_of(flat()), vertices, edges)


# ------------------------------------------- refine_strip: the projected path


class TestRefineStrip:
    @pytest.mark.parametrize("tolerance", [0.0, 0.5, 2.0])
    def test_e1_e2_and_delaunay_by_the_oracles(
        self, generator: Factory, refine_strip: Factory, tolerance: float
    ) -> None:
        """E2 holds for a start that is refine's output at the run's own
        tolerance, as the CLI gives it (15f's S1): the triangles the run does
        not write keep refine's guarantee at that tolerance."""
        start = Start(tolerance=tolerance)
        strip = start.strip(generator)
        out = refine_strip(start.view, strip, *start.args(), tolerance=tolerance)
        assert out.ok(), out.message
        v, oz = np.asarray(out.vertices), np.asarray(out.z)
        assert_array_equal(v[: len(start.args()[0])], start.args()[0])
        found = strip_check(start, out, tolerance)
        assert found.unlocated == 0, found
        assert found.over == 0, found
        nodes = node_findings(
            start.grid, v, oz, np.asarray(out.valid), np.asarray(out.triangles), tolerance=tolerance
        )
        assert nodes.nodes > 0 and nodes.over == 0, nodes
        assert delaunay_check(start, out) == 0
        if tolerance > 0:
            assert out.strip_refused == 0

    def test_without_the_strip_run_the_oracle_finds_points_over(self, start: Start) -> None:
        """The control: refine's output alone, at the tolerance it was refined
        to, against the same oracle. The green above is the strip run's doing."""
        found = strip_check(start, start.out, 2.0)
        assert found.unlocated == 0
        assert found.over > 0, found

    def test_the_outcome_fields(
        self, start: Start, generator: Factory, refine_strip: Factory
    ) -> None:
        strip = start.strip(generator)
        out = refine_strip(start.view, strip, *start.args(), tolerance=0.5)
        assert out.ok(), out.message
        assert isinstance(out, _core.PointRefineOutcome)
        assert out.strip_points == strip.size > 0
        assert out.strip_inserted > 0
        assert out.strip_inserted + out.nodes_inserted == out.inserted
        assert 0.0 <= out.strip_max_error <= 0.5
        assert (out.strip_refused, out.strip_refused_max_error) == (0, 0.0)
        assert 0.0 <= out.max_error <= 0.5
        assert out.uncovered == 0
        for name in ("strip_points", "strip_inserted", "strip_refused", "nodes_inserted"):
            assert isinstance(getattr(out, name), int), name
        for name in ("strip_max_error", "strip_refused_max_error"):
            assert isinstance(getattr(out, name), float), name

    def test_e7_nothing_to_do_nothing_done(
        self, start: Start, generator: Factory, refine_strip: Factory
    ) -> None:
        strip = start.strip(generator)
        out = refine_strip(start.view, strip, *start.args(), tolerance=1e6)
        assert out.ok(), out.message
        assert out.inserted == 0
        for name, given in zip(
            ("vertices", "triangles", "z", "valid", "edges", "masks"), start.args(), strict=True
        ):
            assert_array_equal(np.asarray(getattr(out, name)), given, name)

    @pytest.mark.parametrize("threads", [2, 8])
    def test_e4_bit_identical_for_any_thread_count_and_edge_order(
        self, start: Start, generator: Factory, refine_strip: Factory, threads: int
    ) -> None:
        v, t, z, ok, e, k = start.args()
        ref = refine_strip(
            start.view, start.strip(generator), v, t, z, ok, e, k, tolerance=0.5, threads=1
        )
        order = np.random.default_rng(threads).permutation(len(e))
        e2, k2 = e[order][:, ::-1].copy(), k[order].copy()
        strip = generator(start.view, v, e2)
        out = refine_strip(start.view, strip, v, t, z, ok, e2, k2, tolerance=0.5, threads=threads)
        for name in ("vertices", "z", "valid", "triangles"):
            assert_array_equal(np.asarray(getattr(out, name)), np.asarray(getattr(ref, name)), name)
        for name in ("rounds", "inserted", "strip_inserted", "nodes_inserted", "max_error"):
            assert getattr(out, name) == getattr(ref, name), name


class TestRefineStripRefusals:
    @pytest.mark.parametrize("tolerance", [-1.0, float("nan"), float("inf")])
    def test_a_bad_tolerance_is_a_status(
        self, start: Start, generator: Factory, refine_strip: Factory, tolerance: float
    ) -> None:
        """L2 (1): a status, not a throw, with refine's wording."""
        out = refine_strip(start.view, start.strip(generator), *start.args(), tolerance=tolerance)
        assert not out.ok()
        assert out.status == _core.RefineStatus.InvalidTolerance
        assert out.message == "refine_strip: tolerance must be finite and >= 0"

    def test_a_strip_on_another_geometry(
        self, start: Start, generator: Factory, refine_strip: Factory
    ) -> None:
        """L2 (2): `std::logic_error`, a RuntimeError here."""
        moved = _core.raster_view(start.dem, x_min=X_MIN - H, y_max=Y_MAX, delta_x=H, delta_y=H)
        strip = generator(moved, start.args()[0], start.args()[4])
        with pytest.raises(RuntimeError, match="refine_strip"):
            refine_strip(start.view, strip, *start.args(), tolerance=0.5)

    def test_a_strip_edge_that_is_not_a_constraint_edge(
        self, start: Start, generator: Factory, refine_strip: Factory
    ) -> None:
        """L2 (3): a strip built on a triangle's unconstrained edge."""
        v, t, *_ = start.args()
        constrained = {tuple(sorted(map(int, e))) for e in start.args()[4]}
        loose = next(
            (int(a), int(b))
            for a, b in ((t[i, j], t[i, (j + 1) % 3]) for i in range(len(t)) for j in range(3))
            if tuple(sorted((int(a), int(b)))) not in constrained
        )
        strip = generator(start.view, v, np.array([loose], np.uint32))
        with pytest.raises(RuntimeError, match="refine_strip"):
            refine_strip(start.view, strip, *start.args(), tolerance=0.5)

    @pytest.mark.parametrize(
        "index",
        [0, 1, 2, 3, 4, 5],
        ids=[
            "vertices-3-columns",
            "triangles-2-columns",
            "z-short",
            "valid-short",
            "edges-3-columns",
            "masks-short",
        ],
    )
    def test_shapes(
        self, start: Start, generator: Factory, refine_strip: Factory, index: int
    ) -> None:
        args = list(start.args())
        bad = args[index]
        if bad.ndim == 2:
            args[index] = np.hstack([bad, bad[:, :1]]) if index != 1 else bad[:, :2].copy()
        else:
            args[index] = bad[:-1].copy()
        with pytest.raises(ValueError, match="refine_strip"):
            refine_strip(start.view, start.strip(generator), *args, tolerance=0.5)

    def test_a_triangle_index_out_of_range(
        self, start: Start, generator: Factory, refine_strip: Factory
    ) -> None:
        v, t, z, ok, e, k = start.args()
        t = t.copy()
        t[0, 0] = len(v)
        with pytest.raises(ValueError, match="refine_strip"):
            refine_strip(start.view, start.strip(generator), v, t, z, ok, e, k, tolerance=0.5)

    def test_tolerance_and_threads_are_keyword_only(
        self, start: Start, generator: Factory, refine_strip: Factory
    ) -> None:
        with pytest.raises(TypeError):
            refine_strip(start.view, start.strip(generator), *start.args(), 0.5)

    def test_the_strip_must_be_a_constraint_store(
        self, start: Start, refine_strip: Factory
    ) -> None:
        not_a_strip = start.store(*scattered(start.n, 1, seed=1))
        with pytest.raises(TypeError):
            refine_strip(start.view, not_a_strip, *start.args(), tolerance=0.5)


# ----------------------------- refine_points(..., strip): the reprojected path


class TestRefinePointsWithAStrip:
    @pytest.mark.parametrize("tolerance", [0.0, 0.5, 2.0])
    def test_e1_j2_and_delaunay_by_the_oracles(
        self, start: Start, generator: Factory, tolerance: float
    ) -> None:
        xy, z = scattered(start.n, 1, seed=15)
        strip = start.strip(generator)
        out = _core.refine_points(
            start.store(xy, z), *start.args(), tolerance=tolerance, strip=strip
        )
        assert out.ok(), out.message
        v, t = np.asarray(out.vertices), np.asarray(out.triangles)
        oz, ok = np.asarray(out.z), np.asarray(out.valid)
        given = start.args()[0]
        assert_array_equal(v[: len(given)], given)
        found = strip_check(start, out, tolerance)
        assert found.unlocated == 0 and found.over == 0, found
        excess = worst_excess(xy, z, v, t, oz, ok, tolerance, skip=given)
        assert excess <= 1e-9, f"a check point is {excess} m over tolerance"
        assert delaunay_check(start, out) == 0
        assert out.strip_points == strip.size
        assert out.nodes_inserted == 0
        if tolerance > 0:
            assert out.strip_refused == 0

    def test_without_a_strip_the_strip_fields_are_zero(self, start: Start) -> None:
        out = _core.refine_points(
            start.store(*scattered(start.n, 1, seed=15)), *start.args(), tolerance=0.5
        )
        assert out.ok(), out.message
        assert (out.strip_points, out.strip_inserted, out.strip_refused, out.nodes_inserted) == (
            0,
            0,
            0,
            0,
        )
        assert (out.strip_max_error, out.strip_refused_max_error) == (0.0, 0.0)

    def test_a_strip_on_another_geometry(self, start: Start, generator: Factory) -> None:
        moved = _core.raster_view(start.dem, x_min=X_MIN, y_max=Y_MAX + H, delta_x=H, delta_y=H)
        strip = generator(moved, start.args()[0], start.args()[4])
        with pytest.raises(RuntimeError, match="refine_points"):
            _core.refine_points(
                start.store(*scattered(start.n, 1, seed=15)),
                *start.args(),
                tolerance=0.5,
                strip=strip,
            )

    def test_a_strip_edge_that_is_not_a_constraint_edge(
        self, start: Start, generator: Factory
    ) -> None:
        v, t, *_ = start.args()
        constrained = {tuple(sorted(map(int, e))) for e in start.args()[4]}
        loose = next(
            (int(a), int(b))
            for a, b in ((t[i, j], t[i, (j + 1) % 3]) for i in range(len(t)) for j in range(3))
            if tuple(sorted((int(a), int(b)))) not in constrained
        )
        strip = generator(start.view, v, np.array([loose], np.uint32))
        with pytest.raises(RuntimeError, match="refine_points"):
            _core.refine_points(
                start.store(*scattered(start.n, 1, seed=15)),
                *start.args(),
                tolerance=0.5,
                strip=strip,
            )


# ------------------------------------------------------------------ the GIL


class TestTheGilIsReleased:
    def test_by_the_generator(self, generator: Factory) -> None:
        """1,500 edges across a 1,001-node grid: about 4.5 M points, about
        0.12 s at 27 ns a point (measured on the C++ alone, arm64, -O2)."""
        n, k = 1001, 1500
        dem = np.zeros((n, n), np.float32)
        row0 = 0.5 + (n - 2.0) * np.arange(k) / k
        ends = np.empty((2 * k, 2))
        ends[0::2] = world(np.full(k, 0.25), row0)
        ends[1::2] = world(np.full(k, n - 1.25), (n - 1.5) - (row0 - 0.5))
        edges = np.arange(2 * k, dtype=np.uint32).reshape(k, 2)
        view = view_of(dem)
        strip, ticks, elapsed = ticks_during(lambda: generator(view, ends, edges))
        assert strip.size > 1_000_000
        assert elapsed >= 0.05, (
            f"constraint_check_points took only {elapsed:.3f}s -- too fast to measure GIL "
            "release; enlarge the edges rather than lowering this bound"
        )
        assert ticks >= RELEASED_TICKS, (
            f"only {ticks} ticks in {elapsed:.3f}s: constraint_check_points appears to hold the GIL"
        )

    @staticmethod
    def square(n: int, generator: Factory) -> tuple[Any, Any, tuple[np.ndarray, ...]]:
        """A two-triangle start with corners a quarter cell inside a random
        n x n DEM, its four sides constrained, and the strip on them. At
        tolerance 0 the strip run reaches every node: about 0.14 s for n = 301
        (C++ alone, one thread), so n = 401."""
        rng = np.random.default_rng(9)
        dem = (rng.integers(0, 6400, (n, n)) / 64.0).astype(np.float32)
        a, b = 0.25, n - 1.25
        corners = world([a, a, b, b], [a, b, b, a])
        triangles = np.array([[0, 1, 2], [0, 2, 3]], np.uint32)
        edges = np.array([[0, 1], [1, 2], [2, 3], [0, 3]], np.uint32)
        args = (
            corners,
            triangles,
            np.full(4, 10.0),
            np.ones(4, bool),
            edges,
            np.ones(4, np.uint32),
        )
        view = view_of(dem)
        return view, generator(view, corners, edges), args

    def test_by_refine_strip(self, generator: Factory, refine_strip: Factory) -> None:
        view, strip, args = self.square(401, generator)
        out, ticks, elapsed = ticks_during(
            lambda: refine_strip(view, strip, *args, tolerance=0.0, threads=1)
        )
        assert out.ok(), out.message
        assert elapsed >= 0.05, (
            f"refine_strip took only {elapsed:.3f}s -- too fast to measure GIL release; "
            "enlarge the grid rather than lowering this bound"
        )
        assert ticks >= RELEASED_TICKS, (
            f"only {ticks} ticks in {elapsed:.3f}s: refine_strip appears to hold the GIL"
        )

    def test_by_refine_points_with_a_strip(self, generator: Factory) -> None:
        """15c-1's workload (two check points a cell, every one inserted at
        tolerance 0, from a two-triangle start) with the strip on its sides."""
        n = 301
        _, strip, args = self.square(n, generator)
        rng = np.random.default_rng(9)
        k = 2 * (n - 1) ** 2
        xy = world(
            rng.integers(1, (n - 1) * 1024, k) / 1024.0, rng.integers(1, (n - 1) * 1024, k) / 1024.0
        )
        z = rng.integers(0, 6400, k).astype(np.float32) / np.float32(64.0)
        cp = _core.CheckPoints(x_min=X_MIN, y_max=Y_MAX, spacing=H, rows=n, cols=n)
        cp.add(xy, z)
        cp.freeze()
        out, ticks, elapsed = ticks_during(
            lambda: _core.refine_points(cp, *args, tolerance=0.0, threads=1, strip=strip)
        )
        assert out.ok(), out.message
        assert out.strip_points == strip.size > 0
        assert elapsed >= 0.05, f"refine_points took only {elapsed:.3f}s; enlarge the store"
        assert ticks >= RELEASED_TICKS, (
            f"only {ticks} ticks in {elapsed:.3f}s: refine_points with a strip appears to hold "
            "the GIL"
        )


# ------------------------------------------------------------------ the stubs


class TestStubs:
    """`_core.pyi` is hand-written; mypy checks consistency, not completeness."""

    def declared(self) -> dict[str, ast.ClassDef | ast.FunctionDef]:
        tree = ast.parse(STUB.read_text(encoding="utf-8"))
        return {
            node.name: node
            for node in tree.body
            if isinstance(node, ast.ClassDef | ast.FunctionDef)
        }

    @pytest.mark.parametrize(
        "name", ["ConstraintCheckPoints", "constraint_check_points", "refine_strip"]
    )
    def test_the_stub_declares_the_new_surface(self, name: str) -> None:
        assert name in self.declared()

    def test_the_store_has_its_four_properties(self) -> None:
        cls = self.declared().get("ConstraintCheckPoints")
        assert isinstance(cls, ast.ClassDef), "ConstraintCheckPoints is not stubbed"
        members = {node.name for node in cls.body if isinstance(node, ast.FunctionDef)}
        assert {"size", "no_data", "duplicates", "edge_count"} <= members

    def test_the_generators_parameters(self) -> None:
        fn = self.declared().get("constraint_check_points")
        assert isinstance(fn, ast.FunctionDef), "constraint_check_points is not stubbed"
        assert [a.arg for a in fn.args.args] == ["view", "vertices", "edges"]

    def test_refine_strips_parameters(self) -> None:
        fn = self.declared().get("refine_strip")
        assert isinstance(fn, ast.FunctionDef), "refine_strip is not stubbed"
        positional = [a.arg for a in fn.args.args]
        assert positional == [
            "view",
            "strip",
            "vertices",
            "triangles",
            "z",
            "valid",
            "edges",
            "masks",
        ]
        # 23b's N18 adds frozen_mask; N18 gives no order, so it goes last, after
        # threads, as on refine_points (test_core_frozen.py).
        assert [a.arg for a in fn.args.kwonlyargs] == ["tolerance", "threads", "frozen_mask"]

    def test_the_outcome_carries_the_strip_fields(self) -> None:
        cls = self.declared().get("PointRefineOutcome")
        assert isinstance(cls, ast.ClassDef), "PointRefineOutcome is not stubbed"
        members = {node.name for node in cls.body if isinstance(node, ast.FunctionDef)}
        assert {
            "strip_points",
            "strip_inserted",
            "strip_max_error",
            "strip_refused",
            "strip_refused_max_error",
            "nodes_inserted",
        } <= members


# ------------------------------------------------------------------ the oracles can fail


class TestTheOraclesCanFail:
    """Each oracle above, shown failing on a planted defect, so its green
    means something. They need no strip, and passed at the red step too."""

    def test_the_strip_oracle_on_a_planted_shift(self, start: Start) -> None:
        xy, z = start.oracle_points()
        v, _, oz, _, e, _ = start.args()
        found = strip_findings(xy, z, v, oz + 6.0, e, tolerance=2.0, on_edge=ON_EDGE, slack=0.0)
        assert found.over > 0

    def test_the_node_oracle_on_a_planted_shift(self, start: Start) -> None:
        v, t, oz, ok, _, _ = start.args()
        assert node_findings(start.grid, v, oz, ok, t, tolerance=2.0).over == 0
        assert node_findings(start.grid, v, oz + 6.0, ok, t, tolerance=2.0).over > 0

    def test_the_delaunay_oracle_on_a_planted_flip(self) -> None:
        """A convex quad that is not cocircular, triangulated both ways (both
        counter-clockwise): exactly one diagonal is Delaunay."""
        grid = Grid(X_MIN, Y_MAX, H, H, np.zeros((5, 5)))
        v = world([0.0, 2.0, 2.0, 0.2], [0.0, 0.0, 1.0, 0.8])
        good = np.array([[0, 3, 2], [0, 2, 1]])
        bad = np.array([[0, 3, 1], [3, 2, 1]])
        none = np.zeros((0, 2), np.int64)
        flips = [delaunay_violations(grid, v, t, none) for t in (good, bad)]
        assert sorted(flips) == [0, 1], flips
        # A constrained diagonal is exempt.
        assert delaunay_violations(grid, v, bad, np.array([[1, 3]])) == 0
