"""The seam pass through the binding: increment 23b, SP1 to SP4.

`docs/increments/23-basin-scale.md`, "The seam protocol" step 2 and "Tests
@tester can write red" (23b, "C++ Catch2 and through the binding"). The C++
half is `tests/cpp/property/prop_refinement_seam.cpp`; this file carries the
two oracles the design writes in Python: SP1's every-point check against a
bilinear surface computed in NumPy, and SP2's check points recomputed in exact
rationals with `fractions`.

As N8, N12 and N13 ("Settled after 23b's red step") rule it:

    _core.refine_seam(view, a, b, *, tolerance) -> SeamOutcome
        a, b       (x, y) world points, the seam edge's ends
    SeamOutcome
        .a, .b              (x, y) tuples, ordered so a < b by (x, y)
        .z_a, .z_b          float, or None where the end has no height
        .points             (K, 2) float64, read-only: the inserted points, from a
        .z                  (K,) float64, read-only: their heights
        .s                  (K,) float64, read-only: their parameters, strictly
                            increasing in (0, 1)
        .check_points       int: the check points the pass measured
        .no_data            int: check points skipped for a NoData stencil
        .max_error          float: the largest check-point error at the end
    a refused input (tolerance negative or not finite, a == b, an end outside
    the node rectangle) is a ValueError naming refine_seam.

Every new symbol is fetched inside a fixture, so a missing one fails its own
tests and leaves the rest of the session collecting.
"""

from __future__ import annotations

import ast
import math
from collections.abc import Callable
from fractions import Fraction
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from tin_engine import _core

# An exact frame: x = 2 col, y = -row, so every dyadic lattice position is a
# float in the world and back. Non-square cells, so a row/col swap cannot pass.
X_MIN, Y_MAX, DX, DY = 0.0, 0.0, 2.0, 1.0
STUB = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "_core.pyi"

Seam = Callable[..., Any]


@pytest.fixture
def refine_seam() -> Seam:
    fn: Seam = _core.refine_seam  # type: ignore[attr-defined]
    return fn


def view(dem: np.ndarray, nodata: float | None = None) -> Any:
    dem.setflags(write=False)
    if nodata is None:
        return _core.raster_view(dem, x_min=X_MIN, y_max=Y_MAX, delta_x=DX, delta_y=DY)
    return _core.raster_view(dem, x_min=X_MIN, y_max=Y_MAX, delta_x=DX, delta_y=DY, nodata=nodata)


def world(col: float, row: float) -> tuple[float, float]:
    return (X_MIN + col * DX, Y_MAX - row * DY)


def rough(n: int, seed: int) -> np.ndarray:
    """Heights in 1/64 m steps, so float32 holds them exactly."""
    rng = np.random.default_rng(seed)
    return (rng.integers(0, 64 * 64, (n, n)) / 64.0).astype(np.float32)


def bilinear(dem: np.ndarray, col: np.ndarray, row: np.ndarray) -> np.ndarray:
    """The DEM's bilinear surface, written here: the cell holding the point
    (the last cell on the far sides), weights from the fractions in it."""
    rows, cols = dem.shape
    c0 = np.minimum(np.floor(col).astype(int), cols - 2)
    r0 = np.minimum(np.floor(row).astype(int), rows - 2)
    tx, ty = col - c0, row - r0
    z = dem.astype(np.float64)
    return (
        z[r0, c0] * (1 - tx) * (1 - ty)
        + z[r0, c0 + 1] * tx * (1 - ty)
        + z[r0 + 1, c0] * (1 - tx) * ty
        + z[r0 + 1, c0 + 1] * tx * ty
    )


def polyline(out: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The pass's polyline in the lattice: (col, row) and z of a, the points, b."""
    xy = np.vstack([np.array([out.a]), np.asarray(out.points).reshape(-1, 2), np.array([out.b])])
    z = np.concatenate([[out.z_a], np.asarray(out.z), [out.z_b]]).astype(np.float64)
    return (xy[:, 0] - X_MIN) / DX, (Y_MAX - xy[:, 1]) / DY, z


class TestSP1EveryPointOnAGridLine:
    """On a grid-line seam the bilinear surface is linear between nodes, so the
    tolerance holds at every point of the seam, here at 1,001 points."""

    @pytest.mark.parametrize("tolerance", [0.0, 0.5, 2.0, 8.0])
    @pytest.mark.parametrize("vertical", [True, False])
    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_every_point_is_within_tolerance(
        self, refine_seam: Seam, tolerance: float, vertical: bool, seed: int
    ) -> None:
        dem = rough(17, seed)
        a, b = (world(8, 0), world(8, 16)) if vertical else (world(0, 8), world(16, 8))
        out = refine_seam(view(dem), a, b, tolerance=tolerance)
        pc, pr, pz = polyline(out)
        s = np.linspace(0.0, 1.0, 1001)
        col, row = pc[0] + s * (pc[-1] - pc[0]), pr[0] + s * (pr[-1] - pr[0])
        along = np.hypot(pc - pc[0], pr - pr[0])
        at = np.hypot(col - pc[0], row - pr[0])
        assert (np.diff(along) > 0).all(), "the points are not ordered from a to b"
        err = np.abs(bilinear(dem, col, row) - np.interp(at, along, pz))
        assert err.max() <= tolerance + 1e-9, f"{err.max()} m off at s = {s[err.argmax()]}"

    def test_nodes_go_in_at_their_own_height(self, refine_seam: Seam) -> None:
        dem = rough(17, 4)
        out = refine_seam(view(dem), world(8, 0), world(8, 16), tolerance=0.0)
        pc, pr, pz = polyline(out)
        nodes = (pc == np.round(pc)) & (pr == np.round(pr))
        assert nodes.sum() >= 3
        assert_array_equal(pz[nodes], dem[pr[nodes].astype(int), pc[nodes].astype(int)])


# --------------------------------------------------------------------- SP2


def exact_check_points(
    a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]
) -> list[tuple[Fraction, Fraction, Fraction]]:
    """15f's check points for the edge a -> b in exact rationals: (s, col, row)
    for every crossing of a column or row line strictly between the ends (one
    point where both families cross at a node), and the midpoint of each two
    neighbours of (a, crossings..., b). a and b in lattice (col, row)."""
    (c0, r0), (c1, r1) = a, b
    cuts: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for k in range(math.floor(min(c0, c1)) + 1, math.ceil(max(c0, c1))):
        t = (k - c0) / (c1 - c0)
        cuts[t] = (Fraction(k), r0 + t * (r1 - r0))
    for r in range(math.floor(min(r0, r1)) + 1, math.ceil(max(r0, r1))):
        t = (r - r0) / (r1 - r0)
        cuts[t] = (c0 + t * (c1 - c0), Fraction(r))
    chain = (
        [(Fraction(0), c0, r0)] + [(t, *cuts[t]) for t in sorted(cuts)] + [(Fraction(1), c1, r1)]
    )
    out = [chain[i] for i in range(1, len(chain) - 1)]
    for (tl, cl, rl), (tr, cr, rr) in pairwise(chain):
        out.append(((tl + tr) / 2, (cl + cr) / 2, (rl + rr) / 2))
    return sorted(out)


def exact_height(dem: np.ndarray, col: Fraction, row: Fraction) -> Fraction:
    rows, cols = dem.shape
    if col.denominator == 1 and row.denominator == 1:
        return Fraction(float(dem[int(row), int(col)]))
    c0, r0 = min(math.floor(col), cols - 2), min(math.floor(row), rows - 2)
    tx, ty = col - c0, row - r0
    z = [[Fraction(float(dem[r0 + i, c0 + j])) for j in (0, 1)] for i in (0, 1)]
    return (
        z[0][0] * (1 - tx) * (1 - ty)
        + z[0][1] * tx * (1 - ty)
        + z[1][0] * (1 - tx) * ty
        + z[1][1] * tx * ty
    )


def lattice(p: tuple[float, float]) -> tuple[Fraction, Fraction]:
    return (
        (Fraction(p[0]) - Fraction(X_MIN)) / Fraction(DX),
        (Fraction(Y_MAX) - Fraction(p[1])) / Fraction(DY),
    )


ENDS = [
    ((1.25, 0.5), (14.75, 15.5)),
    ((15.5, 2.125), (0.375, 9.0)),
    ((2.0, 2.0), (14.0, 14.0)),  # the diagonal, through nodes
    ((4.0, 1.0), (12.0, 15.0)),  # through nodes every other row
    ((7.25, 7.125), (8.875, 7.5)),  # under two cells long
]


class TestSP2ExactCheckPoints:
    @pytest.mark.parametrize("tolerance", [0.0, 0.25, 1.0, 5.0])
    @pytest.mark.parametrize(("ca", "cb"), ENDS)
    def test_every_check_point_is_within_tolerance(
        self,
        refine_seam: Seam,
        tolerance: float,
        ca: tuple[float, float],
        cb: tuple[float, float],
    ) -> None:
        dem = rough(17, 7)
        out = refine_seam(view(dem), world(*ca), world(*cb), tolerance=tolerance)
        a, b = lattice(out.a), lattice(out.b)
        assert out.a <= out.b
        points = exact_check_points(a, b)
        assert out.check_points == len(points)

        # The output's polyline, its parameters by exact projection onto a -> b.
        ux, uy = b[0] - a[0], b[1] - a[1]

        def param(p: tuple[Fraction, Fraction]) -> Fraction:
            return ((p[0] - a[0]) * ux + (p[1] - a[1]) * uy) / (ux * ux + uy * uy)

        knots = [(Fraction(0), Fraction(out.z_a))]
        knots += [
            (param(lattice(tuple(p))), Fraction(float(z)))
            for p, z in zip(out.points, out.z, strict=True)
        ]
        knots += [(Fraction(1), Fraction(out.z_b))]
        assert all(x[0] < y[0] for x, y in pairwise(knots))

        def mesh(s: Fraction) -> Fraction:
            for (sl, zl), (sr, zr) in pairwise(knots):
                if sl <= s <= sr:
                    return zl + (s - sl) / (sr - sl) * (zr - zl)
            raise AssertionError(f"parameter {s} outside [0, 1]")

        worst = max(
            (abs(exact_height(dem, c, r) - mesh(s)) for s, c, r in points), default=Fraction(0)
        )
        assert float(worst) <= tolerance + 1e-9, f"a check point is {float(worst)} m off"
        assert abs(out.max_error - float(worst)) <= 1e-9 * max(1.0, float(worst))

        # Every inserted point is one of the check points, at its exact height
        # to rounding: the pass inserts nothing else.
        for p, z in zip(out.points, out.z, strict=True):
            c, r = lattice(tuple(p))
            near = min(points, key=lambda q: abs(q[1] - c) + abs(q[2] - r))
            assert abs(float(near[1] - c)) + abs(float(near[2] - r)) <= 1e-9
            assert abs(float(exact_height(dem, near[1], near[2])) - float(z)) <= 1e-9

    def test_ties_go_to_the_smallest_parameter(self, refine_seam: Seam) -> None:
        # Row 0, cols 0 to 4, z 0 10 0 10 0: nodes 1 and 3 tie at 10; at
        # tolerance 7 the first insertion settles the other, so the tie
        # decides: col 1, from a = col 0, whichever end is given first.
        dem = np.array([[0, 10, 0, 10, 0]] * 2, dtype=np.float32)
        for a, b in ((world(0, 0), world(4, 0)), (world(4, 0), world(0, 0))):
            out = refine_seam(view(dem.copy()), a, b, tolerance=7.0)
            assert out.a == world(0, 0)
            assert_array_equal(np.asarray(out.points), np.array([world(1, 0)]))
            assert_array_equal(np.asarray(out.z), [10.0])
            assert_array_equal(np.asarray(out.s), [0.25])


class TestSP3Sameness:
    @pytest.mark.parametrize(("ca", "cb"), [((8.0, 0.0), (8.0, 16.0)), *ENDS])
    def test_the_edge_reversed_gives_the_same_output(
        self, refine_seam: Seam, ca: tuple[float, float], cb: tuple[float, float]
    ) -> None:
        v = view(rough(17, 9))
        f = refine_seam(v, world(*ca), world(*cb), tolerance=0.5)
        r = refine_seam(v, world(*cb), world(*ca), tolerance=0.5)
        assert (f.a, f.b, f.z_a, f.z_b) == (r.a, r.b, r.z_a, r.z_b)
        assert_array_equal(np.asarray(f.points), np.asarray(r.points))
        assert_array_equal(np.asarray(f.z), np.asarray(r.z))
        assert_array_equal(np.asarray(f.s), np.asarray(r.s))
        assert (f.check_points, f.max_error) == (r.check_points, r.max_error)


class TestSP4TheBound:
    @pytest.mark.parametrize(("ca", "cb"), [((8.0, 0.0), (8.0, 16.0)), *ENDS])
    def test_insertions_never_exceed_check_points(
        self, refine_seam: Seam, ca: tuple[float, float], cb: tuple[float, float]
    ) -> None:
        out = refine_seam(view(rough(17, 11)), world(*ca), world(*cb), tolerance=0.0)
        assert 0 < len(out.points) <= out.check_points

    def test_nodata_stencils_are_skipped_and_counted(self, refine_seam: Seam) -> None:
        dem = rough(17, 12)
        clean = refine_seam(view(dem.copy()), world(8, 0), world(8, 16), tolerance=0.0)
        dem[8, 9] = -9999.0  # one column off the seam
        out = refine_seam(view(dem, nodata=-9999.0), world(8, 0), world(8, 16), tolerance=0.0)
        assert clean.no_data == 0
        assert out.no_data == 2  # the midpoints at rows 7.5 and 8.5; nodes keep their own value
        assert out.check_points == clean.check_points - 2
        rows = (Y_MAX - np.asarray(out.points)[:, 1]) / DY
        assert not np.isin(rows, [7.5, 8.5]).any()

    def test_an_end_on_nodata_has_no_height(self, refine_seam: Seam) -> None:
        dem = rough(17, 13)
        dem[0, 8] = -9999.0
        out = refine_seam(view(dem, nodata=-9999.0), world(8, 0), world(8, 16), tolerance=1e9)
        assert out.b == world(8, 0)
        assert out.z_b is None
        assert out.z_a == float(dem[16, 8])


class TestTheBinding:
    def test_arrays_are_read_only_float64(self, refine_seam: Seam) -> None:
        out = refine_seam(view(rough(17, 1)), world(1.25, 0.5), world(14.75, 15.5), tolerance=0.0)
        k = len(out.z)
        assert k > 0
        for arr, shape in ((out.points, (k, 2)), (out.z, (k,)), (out.s, (k,))):
            a = np.asarray(arr)
            assert a.dtype == np.float64
            assert a.shape == shape
            with pytest.raises(ValueError, match="read-only"):
                a[0] = 0.0
        assert isinstance(out.check_points, int)
        assert isinstance(out.no_data, int)
        assert isinstance(out.max_error, float)

    def test_no_insertion_gives_empty_arrays(self, refine_seam: Seam) -> None:
        plane = np.add.outer(np.arange(9.0), 2.0 * np.arange(9.0)).astype(np.float32)
        out = refine_seam(view(plane), world(0, 4), world(8, 4), tolerance=0.5)
        assert np.asarray(out.points).shape == (0, 2)
        assert np.asarray(out.z).shape == (0,)
        assert out.max_error == 0.0

    @pytest.mark.parametrize("tolerance", [-1.0, math.nan, math.inf])
    def test_a_bad_tolerance_is_a_value_error(self, refine_seam: Seam, tolerance: float) -> None:
        with pytest.raises(ValueError, match=r"refine_seam.*tolerance"):
            refine_seam(view(rough(9, 1)), world(1, 1), world(5, 3), tolerance=tolerance)

    def test_a_degenerate_edge_is_a_value_error(self, refine_seam: Seam) -> None:
        with pytest.raises(ValueError, match="refine_seam"):
            refine_seam(view(rough(9, 1)), world(1, 1), world(1, 1), tolerance=0.5)

    def test_an_end_outside_the_grid_is_a_value_error(self, refine_seam: Seam) -> None:
        with pytest.raises(ValueError, match=r"refine_seam.*outside"):
            refine_seam(view(rough(9, 1)), world(1, 1), world(9.5, 3), tolerance=0.5)


class TestTheStub:
    def test_the_stub_declares_refine_seam_and_its_outcome(self) -> None:
        tree = ast.parse(STUB.read_text())
        funcs = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
        classes = {n.name: n for n in tree.body if isinstance(n, ast.ClassDef)}
        assert "refine_seam" in funcs
        fn = funcs["refine_seam"]
        assert [a.arg for a in fn.args.args] == ["view", "a", "b"]
        assert [a.arg for a in fn.args.kwonlyargs] == ["tolerance"]
        assert "SeamOutcome" in classes
        members = {n.name for n in classes["SeamOutcome"].body if isinstance(n, ast.FunctionDef)}
        assert {
            "a",
            "b",
            "z_a",
            "z_b",
            "points",
            "z",
            "s",
            "check_points",
            "no_data",
            "max_error",
        } <= members


class TestTheStrictInequality:
    """Insert while an error is strictly over the tolerance: an error equal to
    it is within."""

    def test_a_plane_at_tolerance_zero_needs_no_point(self, refine_seam: Seam) -> None:
        dem = np.array([[0, 1, 2, 3, 4]] * 2, dtype=np.float32)
        out = refine_seam(view(dem), world(0, 0), world(4, 0), tolerance=0.0)
        assert out.check_points == 7
        assert np.asarray(out.points).shape == (0, 2)
        assert out.max_error == 0.0

    def test_an_error_equal_to_the_tolerance_needs_no_point(self, refine_seam: Seam) -> None:
        dem = np.array([[0, 2, 0]] * 2, dtype=np.float32)
        out = refine_seam(view(dem), world(0, 0), world(2, 0), tolerance=2.0)
        assert out.check_points == 3
        assert np.asarray(out.points).shape == (0, 2)
        assert out.max_error == 2.0
