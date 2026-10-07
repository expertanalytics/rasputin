"""The final check through the binding: increment 15c-1, RP3 and RP9.

`docs/increments/15c-geographic-dem.md`, D4 to D6 and "Tests for @tester".
`_core.CheckPoints(x_min, y_max, spacing, rows, cols)` files check points by
cell; `.add(xy, z)`, `.freeze()`, `.size`, `.duplicates`, `.outside`.
`_core.refine_points(points, vertices, triangles, z, valid, edges, masks, *,
tolerance, threads)` refines phase 1's mesh against them and returns a
`PointRefineOutcome`: `RefineOutcome`'s fields plus `coincident` and
`coincident_max_error`. Both are bound in `bindings/core.cpp` and stubbed in
`src_python/tin_engine/_core.pyi`.

CHOSEN HERE, where D6 is silent: `add` takes exactly float64 ``(N, 2)`` and
float32 ``(N,)`` and refuses anything else with ValueError (D6: "validated
shapes and dtypes"; a float64 z is the producer's to round, D4); `add` after
`freeze` is a RuntimeError naming the freeze.

RP3 is J2 by an independent oracle in NumPy: it locates every given point in
every output triangle itself (barycentric, 1e-12 slack, so a point on an edge
is tested in both) and never reads the store's order or the scan's records.
Its control plants the output's z shifted by twice the tolerance.

Every new symbol is fetched inside a fixture, so a missing one fails its own
tests and leaves the rest of the session collecting.
"""

from __future__ import annotations

import ast
import pickle
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from gil_probe import ticks_during
from test_core_cdt import RELEASED_TICKS
from tin_engine import _core
from tin_engine._core import ChainRole
from tin_engine.cli import DEFAULT_SNAP_SPACING, _constraint_arrays, _engine

X_MIN, Y_MAX, H = 500000.0, 7000000.0, 30
STUB = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "_core.pyi"

Factory = Callable[..., Any]


@pytest.fixture
def check_points() -> Factory:
    cls: Factory = _core.CheckPoints  # type: ignore[attr-defined]
    return cls


@pytest.fixture
def refine_points() -> Factory:
    fn: Factory = _core.refine_points  # type: ignore[attr-defined]
    return fn


def world(col: np.ndarray, row: np.ndarray) -> np.ndarray:
    return np.column_stack([X_MIN + col * H, Y_MAX - row * H])


def surface(col: np.ndarray, row: np.ndarray) -> np.ndarray:
    return np.asarray(50.0 + 20.0 * np.sin(0.3 * col) + 15.0 * np.cos(0.25 * row))


class Phase1:
    """Phase 1 as the CLI runs it: refine on a square grid, here against the
    smooth surface, from the grid's perimeter ring."""

    def __init__(self, n: int, tolerance: float) -> None:
        self.n = n
        r, c = np.indices((n, n), dtype=np.float64)
        dem = surface(c, r).astype(np.float32)
        dem.setflags(write=False)
        ring = (
            [(0, j) for j in range(0, n - 1, 4)]
            + [(i, n - 1) for i in range(0, n - 1, 4)]
            + [(n - 1, j) for j in range(n - 1, 0, -4)]
            + [(i, 0) for i in range(n - 1, 0, -4)]
        )
        # Counter-clockwise in world: down the west side first.
        ring = ring[::-1]
        xy = world(np.array([c for _, c in ring], float), np.array([r for r, _ in ring], float))
        run = _engine(xy, [([*range(len(ring))], ChainRole.Outer, 0)], True, DEFAULT_SNAP_SPACING)
        assert run.mesh is not None and run.noded is not None, run.message
        edges, masks = _constraint_arrays(run.mesh, run.noded)
        view = _core.raster_view(dem, x_min=X_MIN, y_max=Y_MAX, delta_x=float(H), delta_y=float(H))
        out = _core.refine(view, run.mesh, edges, masks, tolerance=tolerance)
        assert out.ok(), out.message
        self.out = out

    def args(self) -> tuple[np.ndarray, ...]:
        o = self.out
        return tuple(
            np.asarray(a) for a in (o.vertices, o.triangles, o.z, o.valid, o.edges, o.masks)
        )

    def store(self, check_points: Factory, xy: np.ndarray, z: np.ndarray) -> Any:
        cp = check_points(x_min=X_MIN, y_max=Y_MAX, spacing=H, rows=self.n, cols=self.n)
        cp.add(xy, z)
        cp.freeze()
        return cp


def scattered(n: int, per_cell: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Dyadic positions strictly inside cells (k / 1024, k in 1..1023), so each
    stored position is the given one; z the surface plus +-4 m of noise,
    rounded to 1/64 m so it is a float32 exactly."""
    rng = np.random.default_rng(seed)
    cells = (n - 1) * (n - 1) * per_cell
    base_r, base_c = np.divmod(np.repeat(np.arange((n - 1) * (n - 1)), per_cell), n - 1)
    col = base_c + rng.integers(1, 1024, cells) / 1024.0
    row = base_r + rng.integers(1, 1024, cells) / 1024.0
    noise = rng.integers(-256, 257, cells) / 64.0
    z = np.round((surface(col, row) + noise) * 64.0) / 64.0
    return world(col, row), z.astype(np.float32)


def worst_excess(
    xy: np.ndarray,
    zp: np.ndarray,
    vertices: np.ndarray,
    triangles: np.ndarray,
    z: np.ndarray,
    valid: np.ndarray,
    tolerance: float,
    skip: np.ndarray,
) -> float:
    """J2 by brute force: the largest |plane - z_p| - tolerance over every
    (point, triangle) pair whose closed triangle holds the point, the plane
    from the output's z. Points in `skip` (start vertices) are left out."""
    keep = ~(xy[:, None, :] == skip[None, :, :]).all(axis=2).any(axis=1)
    px, py = xy[keep, 0] - X_MIN, xy[keep, 1] - Y_MAX
    pz = zp[keep].astype(np.float64)
    tri = triangles[valid[triangles].all(axis=1)]
    v = vertices - np.array([X_MIN, Y_MAX])
    a, b, c = v[tri[:, 0]], v[tri[:, 1]], v[tri[:, 2]]
    two_a = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    assert (two_a > 0).all(), "an output triangle is not counter-clockwise"
    worst = -np.inf
    for lo in range(0, len(px), 256):
        x, y = px[lo : lo + 256, None], py[lo : lo + 256, None]

        def weight(p: np.ndarray, q: np.ndarray, x: np.ndarray = x, y: np.ndarray = y) -> Any:
            return (
                (q[:, 0] - p[:, 0]) * (y - p[:, 1]) - (q[:, 1] - p[:, 1]) * (x - p[:, 0])
            ) / two_a

        wa, wb, wc = weight(b, c), weight(c, a), weight(a, b)
        inside = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
        plane = wa * z[tri[:, 0]] + wb * z[tri[:, 1]] + wc * z[tri[:, 2]]
        err = np.abs(plane - pz[lo : lo + 256, None])
        if inside.any():
            worst = max(worst, float(err[inside].max()) - tolerance)
    return worst


@pytest.fixture(scope="module")
def phase1() -> Phase1:
    return Phase1(n=25, tolerance=4.0)


class TestJ2ByAnIndependentOracle:
    """RP3."""

    @pytest.mark.parametrize("tolerance", [0.0, 0.5, 2.0, 8.0])
    def test_every_check_point_is_within_tolerance(
        self, phase1: Phase1, check_points: Factory, refine_points: Factory, tolerance: float
    ) -> None:
        xy, z = scattered(phase1.n, 2, seed=15)
        cp = phase1.store(check_points, xy, z)
        out = refine_points(cp, *phase1.args(), tolerance=tolerance, threads=0)
        assert out.ok(), out.message
        # Phase 2 has work exactly when the oracle finds a check point over
        # tolerance on phase 1's mesh alone. On this fixture phase 1's worst is
        # about 7.15 m, so 8.0 is the case where nothing may be inserted.
        p1 = phase1.args()
        work = worst_excess(xy, z, p1[0], p1[1], p1[2], p1[3], tolerance, skip=p1[0]) > 1e-9
        assert (out.inserted > 0) == work, (out.inserted, work)
        assert 0.0 <= out.max_error <= tolerance
        v, t = np.asarray(out.vertices), np.asarray(out.triangles)
        oz, ok = np.asarray(out.z), np.asarray(out.valid)
        start = np.asarray(phase1.out.vertices)
        assert_array_equal(v[: len(start)], start)
        if not work:
            assert_array_equal(t, p1[1])
        excess = worst_excess(xy, z, v, t, oz, ok, tolerance, skip=start)
        assert excess <= 1e-9, f"a check point is {excess} m over tolerance"

    @pytest.mark.parametrize("tolerance", [0.0, 0.5, 2.0, 8.0])
    def test_the_oracle_fails_on_a_planted_shift(
        self, phase1: Phase1, check_points: Factory, refine_points: Factory, tolerance: float
    ) -> None:
        xy, z = scattered(phase1.n, 2, seed=15)
        cp = phase1.store(check_points, xy, z)
        out = refine_points(cp, *phase1.args(), tolerance=tolerance, threads=0)
        assert out.ok(), out.message
        shift = 2.0 * tolerance if tolerance > 0 else 1.0
        planted = np.asarray(out.z) + shift
        start = np.asarray(phase1.out.vertices)
        excess = worst_excess(
            xy,
            z,
            np.asarray(out.vertices),
            np.asarray(out.triangles),
            planted,
            np.asarray(out.valid),
            tolerance,
            skip=start,
        )
        assert excess > 1e-9

    def test_without_phase_2_the_oracle_finds_check_points_over(self, phase1: Phase1) -> None:
        """Phase 1's mesh alone, judged against the check points: the noise is
        +-4 m, so at 0.5 m some are over. RP3's green is phase 2's doing."""
        xy, z = scattered(phase1.n, 2, seed=15)
        v, t, oz, ok, _, _ = phase1.args()
        assert worst_excess(xy, z, v, t, oz, ok, 0.5, skip=v) > 1e-9


class TestTheResult:
    def test_fields_types_and_counts(self, check_points: Factory, refine_points: Factory) -> None:
        p1 = Phase1(n=13, tolerance=4.0)
        xy, z = scattered(p1.n, 1, seed=3)
        # One more exactly at a start vertex, 2.5 m off its z: coincident.
        v0 = np.asarray(p1.out.vertices)[0]
        z0 = np.float32(np.asarray(p1.out.z)[0] + 2.5)
        xy, z = np.vstack([xy, v0]), np.append(z, z0).astype(np.float32)
        cp = p1.store(check_points, xy, z)
        assert cp.size == len(xy)
        assert cp.duplicates == 0
        assert cp.outside == 0
        out = refine_points(cp, *p1.args(), tolerance=1.0)
        assert out.ok(), out.message
        assert out.status == _core.RefineStatus.Ok
        assert np.asarray(out.vertices).dtype == np.float64
        assert np.asarray(out.triangles).dtype == np.uint32
        assert np.asarray(out.z).shape == (len(np.asarray(out.vertices)),)
        assert isinstance(out.inserted, int) and out.inserted > 0
        assert out.coincident == 1
        assert out.coincident_max_error == pytest.approx(
            abs(float(z0) - float(np.asarray(p1.out.z)[0])), abs=1e-12
        )
        for name in ("rounds", "flips", "uncovered", "carved", "max_error"):
            getattr(out, name)

    def test_outside_and_nan_are_counted(self, check_points: Factory) -> None:
        cp = check_points(x_min=X_MIN, y_max=Y_MAX, spacing=H, rows=5, cols=5)
        xy = np.array(
            [[X_MIN - 1.0, Y_MAX - 10.0], [np.nan, Y_MAX - 10.0], [X_MIN + 15.0, Y_MAX - 15.0]]
        )
        cp.add(xy, np.array([1.0, 1.0, np.nan], dtype=np.float32))
        cp.freeze()
        assert (cp.size, cp.outside) == (0, 3)

    @pytest.mark.parametrize("threads", [1, 2, 8])
    def test_bit_identical_for_any_thread_count(
        self, check_points: Factory, refine_points: Factory, threads: int
    ) -> None:
        p1 = Phase1(n=17, tolerance=4.0)
        xy, z = scattered(p1.n, 2, seed=8)
        ref = refine_points(p1.store(check_points, xy, z), *p1.args(), tolerance=1.0, threads=1)
        back = p1.store(check_points, xy[::-1].copy(), z[::-1].copy())
        out = refine_points(back, *p1.args(), tolerance=1.0, threads=threads)
        for name in ("vertices", "z", "valid", "triangles", "edges", "masks"):
            assert_array_equal(np.asarray(getattr(out, name)), np.asarray(getattr(ref, name)), name)
        for name in ("rounds", "inserted", "flips", "max_error", "coincident"):
            assert getattr(out, name) == getattr(ref, name), name


class TestBinding:
    """RP9."""

    @pytest.fixture
    def store(self, check_points: Factory) -> Any:
        return check_points(x_min=X_MIN, y_max=Y_MAX, spacing=H, rows=5, cols=5)

    @pytest.mark.parametrize(
        ("xy", "z"),
        [
            (np.zeros((3, 3)), np.zeros(3, np.float32)),
            (np.zeros(6), np.zeros(3, np.float32)),
            (np.zeros((3, 2)), np.zeros((3, 1), np.float32)),
            (np.zeros((3, 2)), np.zeros(4, np.float32)),
            (np.zeros((3, 2), np.float32), np.zeros(3, np.float32)),
            (np.zeros((3, 2)), np.zeros(3, np.float64)),
            (np.zeros((3, 2), np.int64), np.zeros(3, np.float32)),
        ],
        ids=[
            "xy-3-columns",
            "xy-1d",
            "z-2d",
            "length-mismatch",
            "xy-float32",
            "z-float64",
            "xy-int",
        ],
    )
    def test_add_refuses_shapes_and_dtypes(self, store: Any, xy: np.ndarray, z: np.ndarray) -> None:
        with pytest.raises(ValueError, match="add"):
            store.add(xy, z)
        store.freeze()
        assert store.size == 0

    def test_pickled_float64_xy_and_float32_z_are_accepted(self, store: Any) -> None:
        """An array that crossed a process boundary (multiprocessing pickles
        it) has an equal dtype that is not the same object. The check is on
        what the dtype is, not on which object it is."""
        xy = pickle.loads(pickle.dumps(np.array([[X_MIN + 15.0, Y_MAX - 15.0]])))
        z = pickle.loads(pickle.dumps(np.array([7.0], np.float32)))
        assert xy.dtype == np.float64 and z.dtype == np.float32
        store.add(xy, z)
        store.freeze()
        assert (store.size, store.outside) == (1, 0)

    @pytest.mark.parametrize("dtype", [np.float32, np.int64], ids=["float32", "int64"])
    def test_pickled_xy_of_another_dtype_is_still_refused(self, store: Any, dtype: type) -> None:
        xy = pickle.loads(pickle.dumps(np.array([[15, 15]], dtype)))
        z = pickle.loads(pickle.dumps(np.ones(1, np.float32)))
        with pytest.raises(ValueError, match="add"):
            store.add(xy, z)

    def test_refine_points_refuses_a_store_never_frozen(
        self, check_points: Factory, refine_points: Factory
    ) -> None:
        p1 = Phase1(n=9, tolerance=4.0)
        cp = check_points(x_min=X_MIN, y_max=Y_MAX, spacing=H, rows=p1.n, cols=p1.n)
        cp.add(*scattered(p1.n, 1, seed=1))
        with pytest.raises(RuntimeError, match="not frozen"):
            refine_points(cp, *p1.args(), tolerance=1.0)

    def test_add_after_freeze_is_refused(self, store: Any) -> None:
        store.freeze()
        with pytest.raises(RuntimeError, match=r"freeze|frozen"):
            store.add(np.array([[X_MIN + 15.0, Y_MAX - 15.0]]), np.ones(1, np.float32))

    def test_tolerance_and_threads_are_keyword_only(
        self, check_points: Factory, refine_points: Factory
    ) -> None:
        p1 = Phase1(n=9, tolerance=4.0)
        cp = p1.store(check_points, *scattered(p1.n, 1, seed=1))
        with pytest.raises(TypeError):
            refine_points(cp, *p1.args(), 1.0)

    def test_refine_points_releases_the_gil(
        self, check_points: Factory, refine_points: Factory
    ) -> None:
        # From a two-triangle start over a 301-node grid, every check point
        # (two per cell) goes in at tolerance 0: enough serial work to watch.
        n = 301
        corners = world(
            np.array([0.0, 0.0, n - 1.0, n - 1.0]), np.array([0.0, n - 1.0, n - 1.0, 0.0])
        )
        triangles = np.array([[0, 1, 2], [0, 2, 3]], np.uint32)
        zc, valid = np.zeros(4), np.ones(4, bool)
        edges = np.array([[0, 1], [1, 2], [2, 3], [0, 3]], np.uint32)
        masks = np.ones(4, np.uint32)
        rng = np.random.default_rng(9)
        k = 2 * (n - 1) ** 2
        xy = world(
            rng.integers(1, (n - 1) * 1024, k) / 1024.0, rng.integers(1, (n - 1) * 1024, k) / 1024.0
        )
        z = rng.integers(0, 6400, k).astype(np.float32) / np.float32(64.0)
        cp = check_points(x_min=X_MIN, y_max=Y_MAX, spacing=H, rows=n, cols=n)
        cp.add(xy, z)
        cp.freeze()
        out, ticks, elapsed = ticks_during(
            lambda: refine_points(
                cp, corners, triangles, zc, valid, edges, masks, tolerance=0.0, threads=1
            )
        )
        assert out.ok(), out.message
        assert elapsed >= 0.05, (
            f"refine_points took only {elapsed:.3f}s -- too fast to measure GIL release; "
            "enlarge the store rather than lowering this bound"
        )
        assert ticks >= RELEASED_TICKS, (
            f"only {ticks} ticks in {elapsed:.3f}s: refine_points appears to hold the GIL"
        )


class TestStubs:
    """`_core.pyi` is hand-written; mypy checks consistency, not completeness."""

    def declared(self) -> dict[str, ast.ClassDef | ast.FunctionDef]:
        tree = ast.parse(STUB.read_text(encoding="utf-8"))
        return {
            node.name: node
            for node in tree.body
            if isinstance(node, ast.ClassDef | ast.FunctionDef)
        }

    @pytest.mark.parametrize("name", ["CheckPoints", "PointRefineOutcome", "refine_points"])
    def test_the_stub_declares_the_new_surface(self, name: str) -> None:
        assert name in self.declared()

    def test_refine_points_takes_tolerance_threads_strip_and_frozen_mask_by_keyword(self) -> None:
        """`strip` since increment 15f-3 (D5; keyword-only, as chosen in
        `test_core_edge_strip.py`), `frozen_mask` since 23b (N13),
        `constraint_feet` since 20c-1 (pin 10)."""
        fn = self.declared().get("refine_points")
        assert isinstance(fn, ast.FunctionDef), "refine_points is not stubbed"
        # 23b (N13) adds frozen_mask, after 15f-3's strip (merged first);
        # 20c-1 (pin 10) adds constraint_feet, last, after frozen_mask.
        assert [a.arg for a in fn.args.kwonlyargs] == [
            "tolerance",
            "threads",
            "strip",
            "frozen_mask",
            "constraint_feet",
        ]

    def test_point_refine_outcome_carries_the_two_counters(self) -> None:
        cls = self.declared().get("PointRefineOutcome")
        assert isinstance(cls, ast.ClassDef), "PointRefineOutcome is not stubbed"
        members = {node.name for node in cls.body if isinstance(node, ast.FunctionDef)}
        assert {"coincident", "coincident_max_error", "vertices", "max_error"} <= members
