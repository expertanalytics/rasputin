"""The error summary through the binding: tests B1 to B3 at the core.

`docs/increments/104-triangle-budget.md`, *4.3 The error summary* and *10. The
red suite*. `_core.error_summary(view_or_store, vertices, triangles, z, valid,
*, threads=0)` scans a mesh as given, once, and returns a read-only
`ErrorSummary` with `max_error`, `rms_error` and `count`: the largest and the
root-mean-square |z - plane| over the DEM nodes (a `RasterView`) or the check
points (a frozen `CheckPoints`) inside the mesh, each counted once. A node on
an edge shared by two triangles belongs to the lower triangle index; a vertex
is never a point; a NoData node is not a point; a void triangle (a NoData
corner) measures nothing.

The oracle (`summary_oracle`) is independent of the scan: it decides
membership with exact integer orientations on the frame scaled by 1024 (every
position here is a multiple of 1/1024 of a cell), takes the owner of a point
as the lowest-indexed triangle whose closed set holds it, and evaluates the
plane from the barycentric weights itself.

Invariant-critical (each node once). Every new symbol is fetched inside a
fixture, so a missing one fails its own tests and leaves the rest collecting.
"""

from __future__ import annotations

import ast
import math
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_allclose

from gil_probe import ticks_during
from test_core_cdt import RELEASED_TICKS
from test_core_refine import rough, start
from test_core_refine_points import X_MIN as PX
from test_core_refine_points import Y_MAX as PY
from test_core_refine_points import H as PH
from test_core_refine_points import Phase1, scattered
from tin_engine import _core

STUB = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "_core.pyi"

#: A small frame of unit cells: world x = X0 + col, y = Y0 - row.
X0, Y0 = 500000.0, 7000000.0
#: Frame positions are multiples of 1 / SCALE, so the oracle's orientations are integers.
SCALE = 1024

Summary = Callable[..., Any]


@pytest.fixture
def error_summary() -> Summary:
    fn: Summary = _core.error_summary  # type: ignore[attr-defined]
    return fn


def view_of(array: npt.NDArray[np.float32], dx: float = 1.0, dy: float = 1.0) -> Any:
    array.setflags(write=False)
    return _core.raster_view(array, x_min=X0, y_max=Y0, delta_x=dx, delta_y=dy)


def world(col: npt.ArrayLike, row: npt.ArrayLike) -> npt.NDArray[np.float64]:
    c, r = np.asarray(col, np.float64), np.asarray(row, np.float64)
    return np.column_stack([X0 + c, Y0 - r])


@dataclass(frozen=True)
class Expected:
    count: int
    max_error: float
    sum_squares: float

    @property
    def rms(self) -> float:
        return math.sqrt(self.sum_squares / self.count) if self.count else 0.0


def summary_oracle(
    points: npt.NDArray[np.float64],
    zp: npt.NDArray[np.float64],
    vertices: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.integer[Any]],
    z: npt.NDArray[np.float64],
    valid: npt.NDArray[np.bool_],
) -> Expected:
    """The summary by brute force. `points` and `vertices` are in the frame
    (col, row); a point equal to a vertex is skipped; a point in no closed
    triangle is outside; a point's owner is the lowest index among the
    triangles holding it (the only one inside, both on a shared edge, the one
    on a border edge); an owner with an invalid corner measures nothing."""
    p = np.rint(np.asarray(points) * SCALE).astype(np.int64)
    v = np.rint(np.asarray(vertices) * SCALE).astype(np.int64)
    assert_allclose(p / SCALE, points, rtol=0, atol=1e-9, err_msg="a point off the 1/1024 frame")
    assert_allclose(v / SCALE, vertices, rtol=0, atol=1e-9, err_msg="a vertex off the 1/1024 frame")
    tri = np.asarray(triangles, np.int64)
    a, b, c = v[tri[:, 0]], v[tri[:, 1]], v[tri[:, 2]]
    two_a = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    assert (two_a != 0).all(), "a zero-area triangle"
    sign = np.sign(two_a)
    is_vertex = (p[:, None, :] == v[None, :, :]).all(axis=2).any(axis=1)
    void = ~(valid[tri].all(axis=1))
    count, top, squares = 0, 0.0, 0.0
    for k in np.nonzero(~is_vertex)[0]:
        x, y = p[k]

        def orient(s: Any, t: Any, x: Any = x, y: Any = y) -> Any:
            return (
                (t[:, 0] - s[:, 0]) * (y - s[:, 1]) - (t[:, 1] - s[:, 1]) * (x - s[:, 0])
            ) * sign

        o_ab, o_bc, o_ca = orient(a, b), orient(b, c), orient(c, a)
        holding = np.nonzero((o_ab >= 0) & (o_bc >= 0) & (o_ca >= 0))[0]
        if holding.size == 0:
            continue
        t = int(holding.min())
        if void[t]:
            continue
        wa, wb, wc = (float(o[t]) / float(two_a[t] * sign[t]) for o in (o_bc, o_ca, o_ab))
        plane = wa * z[tri[t, 0]] + wb * z[tri[t, 1]] + wc * z[tri[t, 2]]
        error = abs(float(zp[k]) - plane)
        count, top, squares = count + 1, max(top, error), squares + error * error
    return Expected(count, top, squares)


def lattice_oracle(
    array: npt.NDArray[np.float32],
    vertices_world: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.integer[Any]],
    z: npt.NDArray[np.float64],
    valid: npt.NDArray[np.bool_],
    dx: float = 1.0,
    dy: float = 1.0,
) -> Expected:
    """`summary_oracle` over every valid node of `array`, vertices given in world."""
    rows, cols = array.shape
    r, c = np.indices((rows, cols))
    keep = ~np.isnan(array.ravel())
    nodes = np.column_stack([c.ravel(), r.ravel()]).astype(np.float64)[keep]
    zp = array.ravel().astype(np.float64)[keep]
    frame = np.column_stack([(vertices_world[:, 0] - X0) / dx, (Y0 - vertices_world[:, 1]) / dy])
    return summary_oracle(nodes, zp, frame, triangles, z, valid)


def agrees(got: Any, want: Expected) -> None:
    assert got.count == want.count
    assert_allclose(got.max_error, want.max_error, rtol=1e-12, atol=1e-12)
    assert_allclose(got.rms_error, want.rms, rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------- B1: a plane


def plane_tile() -> npt.NDArray[np.float32]:
    """3 x 4, as `geotiff_fixtures.elevations()`'s shape, holding the plane
    z = 0.5 col + 0.25 row + 3 (every value exact in float32)."""
    r, c = np.indices((3, 4))
    return np.asarray(0.5 * c + 0.25 * r + 3.0, dtype=np.float32)


def corner_mesh(
    array: npt.NDArray[np.float32],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.uint32], npt.NDArray[np.float64]]:
    """The two-triangle mesh of the grid's corners, counter-clockwise in world:
    the diagonal from the top-left to the bottom-right corner."""
    rows, cols = array.shape
    cr = [(0, 0), (0, rows - 1), (cols - 1, rows - 1), (cols - 1, 0)]
    vertices = world([c for c, _ in cr], [r for _, r in cr])
    z = np.array([array[r, c] for c, r in cr], np.float64)
    return vertices, np.array([[0, 1, 2], [0, 2, 3]], np.uint32), z


class TestB1Plane:
    @pytest.mark.parametrize("threads", [1, 4, 0])
    def test_a_plane_has_no_error_and_eight_nodes(
        self, error_summary: Summary, threads: int
    ) -> None:
        """The twelve nodes less the four vertices; the diagonal (3 by 2 cells)
        passes through no node, and the eight border and inner nodes count once."""
        array = plane_tile()
        vertices, triangles, z = corner_mesh(array)
        got = error_summary(
            view_of(array), vertices, triangles, z, np.ones(4, bool), threads=threads
        )
        assert got.count == 8
        assert got.max_error == 0.0
        assert got.rms_error == 0.0

    def test_the_fields_are_read_only(self, error_summary: Summary) -> None:
        array = plane_tile()
        vertices, triangles, z = corner_mesh(array)
        got = error_summary(view_of(array), vertices, triangles, z, np.ones(4, bool))
        for name in ("max_error", "rms_error", "count"):
            getattr(got, name)  # present, so the AttributeError below is the setter's
            with pytest.raises(AttributeError):
                setattr(got, name, 0)

    def test_threads_is_keyword_only(self, error_summary: Summary) -> None:
        array = plane_tile()
        vertices, triangles, z = corner_mesh(array)
        with pytest.raises(TypeError):
            error_summary(view_of(array), vertices, triangles, z, np.ones(4, bool), 1)

    def test_a_nodata_node_is_not_a_point(self, error_summary: Summary) -> None:
        array = plane_tile()
        array[1, 1] = np.nan
        vertices, triangles, z = corner_mesh(array)
        got = error_summary(view_of(array), vertices, triangles, z, np.ones(4, bool))
        assert got.count == 7
        assert got.max_error == 0.0 and got.rms_error == 0.0

    def test_a_void_triangle_measures_nothing(self, error_summary: Summary) -> None:
        """The bottom-left corner NoData: the lower triangle is void and its four
        nodes (two on the bottom row, one on the left, one inside) are not
        points; the upper triangle's four still are."""
        array = np.random.default_rng(5).uniform(0.0, 50.0, (3, 4)).astype(np.float32)
        array[2, 0] = np.nan
        vertices, triangles, z = corner_mesh(array)
        valid = ~np.isnan(z)
        got = error_summary(view_of(array), vertices, triangles, z, valid)
        want = lattice_oracle(array, vertices, triangles, z, valid)
        assert want.count == 4
        agrees(got, want)
        assert got.rms_error > 0.0

    def test_every_node_a_vertex_counts_none(self, error_summary: Summary) -> None:
        """The saturated mesh (every node a vertex): nothing to measure; the
        root-mean-square of no nodes is reported as 0."""
        array = np.random.default_rng(2).uniform(0.0, 50.0, (2, 2)).astype(np.float32)
        vertices, triangles, z = corner_mesh(array)
        got = error_summary(view_of(array), vertices, triangles, z, np.ones(4, bool))
        assert (got.count, got.max_error, got.rms_error) == (0, 0.0, 0.0)

    def test_no_triangles_counts_none(self, error_summary: Summary) -> None:
        array = plane_tile()
        vertices, _, z = corner_mesh(array)
        empty = np.zeros((0, 3), np.uint32)
        got = error_summary(view_of(array), vertices, empty, z, np.ones(4, bool))
        assert (got.count, got.max_error, got.rms_error) == (0, 0.0, 0.0)


# ---------------------------------------------------------------- B2: each node once


def fan() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.uint32]]:
    """The four-triangle fan of a 5 x 5 grid from its centre: the corners
    (top-left, bottom-left, bottom-right, top-right) and the centre, every
    triangle counter-clockwise in world. Each diagonal half is an interior
    edge with one node on it; the four sides are border edges."""
    cr = [(0, 0), (0, 4), (4, 4), (4, 0), (2, 2)]
    vertices = world([c for c, _ in cr], [r for _, r in cr])
    return vertices, np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]], np.uint32)


def noise(rows: int, cols: int, seed: int) -> npt.NDArray[np.float32]:
    return np.random.default_rng(seed).uniform(0.0, 100.0, (rows, cols)).astype(np.float32)


def corner_z(array: npt.NDArray[np.float32], vertices: npt.NDArray[np.float64]) -> Any:
    c = np.rint(vertices[:, 0] - X0).astype(int)
    r = np.rint(Y0 - vertices[:, 1]).astype(int)
    return array[r, c].astype(np.float64)


class TestB2EachNodeOnce:
    @pytest.mark.parametrize("threads", [1, 3, 0])
    def test_the_fan_counts_every_node_but_the_five_vertices(
        self, error_summary: Summary, threads: int
    ) -> None:
        """25 nodes less 5 vertices: every other node is inside a triangle, on a
        border edge, or on one of the four diagonal halves, each counted once
        (24 if the diagonal nodes count twice, 16 if never)."""
        array = noise(5, 5, 11)
        vertices, triangles = fan()
        z, valid = corner_z(array, vertices), np.ones(5, bool)
        got = error_summary(view_of(array), vertices, triangles, z, valid, threads=threads)
        assert got.count == 25 - 5
        agrees(got, lattice_oracle(array, vertices, triangles, z, valid))

    def test_a_void_owner_keeps_its_edge_nodes(self, error_summary: Summary) -> None:
        """The top-left corner NoData voids triangles 0 and 3. The node on the
        diagonal half shared by triangles 0 and 1 is triangle 0's (the lower
        index), so it is not measured; the one shared by 1 and 2 is triangle
        1's, and the one shared by 2 and 3 triangle 2's, so both are."""
        array = noise(5, 5, 12)
        array[0, 0] = np.nan
        vertices, triangles = fan()
        z = corner_z(array, vertices)
        valid = ~np.isnan(z)
        got = error_summary(view_of(array), vertices, triangles, z, valid)
        want = lattice_oracle(array, vertices, triangles, z, valid)
        # Triangle 1: three bottom-row nodes, one inside, the diagonal node at
        # (col 3, row 3); triangle 2: three right-column nodes, one inside, (3, 1).
        assert want.count == 10
        agrees(got, want)

    def fan_store(self, extra: list[tuple[float, float]]) -> tuple[Any, int, Any, Any]:
        """A store on the 5 x 5 frame: the five vertices' positions, points
        exactly on the four diagonal halves and on the four sides, and points
        inside; dyadic positions, so each stored position is the given one.
        Returns (store, how many coincide with a vertex, positions, z)."""
        coincident = [(0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0), (2.0, 2.0)]
        diagonals = [(t, t) for t in (0.5, 1.0, 1.5, 2.5, 3.0, 3.5)]
        diagonals += [(t, 4.0 - t) for t in (0.5, 1.0, 1.5, 2.5, 3.0, 3.5)]
        sides = [(0.5, 0.0), (4.0, 1.25), (2.75, 4.0), (0.0, 3.5), (2.0, 0.0), (4.0, 2.0)]
        inside = [(1.0, 2.0), (2.0, 3.0), (2.5, 1.0), (3.25, 2.125), (0.75, 1.5), (1.5, 0.25)]
        inside += [(2.0 + 1 / 1024, 2.0)]  # 1/1024 of a cell from the centre: still a point
        at = np.array(coincident + diagonals + sides + inside + extra, np.float64)
        z = np.random.default_rng(13).integers(0, 6400, len(at)).astype(np.float32) / 64
        store = _core.CheckPoints(x_min=X0, y_max=Y0, spacing=1.0, rows=5, cols=5)  # type: ignore[attr-defined]
        store.add(world(at[:, 0], at[:, 1]), z)
        store.freeze()
        assert store.size == len(at) and store.duplicates == 0 and store.outside == 0
        return store, len(coincident), at, z

    @pytest.mark.parametrize("threads", [1, 3, 0])
    def test_the_store_counts_every_point_but_the_coincident(
        self, error_summary: Summary, threads: int
    ) -> None:
        store, coincident, at, zp = self.fan_store([])
        vertices, triangles = fan()
        z = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        valid = np.ones(5, bool)
        got = error_summary(store, vertices, triangles, z, valid, threads=threads)
        assert got.count == len(at) - coincident
        frame = np.array([(0, 0), (0, 4), (4, 4), (4, 0), (2, 2)], np.float64)
        agrees(got, summary_oracle(at, zp.astype(np.float64), frame, triangles, z, valid))

    def test_a_void_owner_keeps_its_edge_points(self, error_summary: Summary) -> None:
        """As on the lattice: an invalid top-left vertex voids triangles 0 and 3,
        and the points on their diagonal halves follow the lower index."""
        store, _, at, zp = self.fan_store([])
        vertices, triangles = fan()
        z = np.array([np.nan, 20.0, 30.0, 40.0, 50.0])
        valid = np.array([False, True, True, True, True])
        got = error_summary(store, vertices, triangles, z, valid)
        frame = np.array([(0, 0), (0, 4), (4, 4), (4, 0), (2, 2)], np.float64)
        want = summary_oracle(at, zp.astype(np.float64), frame, triangles, z, valid)
        assert 0 < want.count < len(at) - 5
        agrees(got, want)

    def test_an_invalid_vertex_with_a_height_is_still_void(self, error_summary: Summary) -> None:
        """`valid` decides, not the height: the top-left vertex marked invalid
        with a finite z gives what the NaN z gives."""
        store, _, _, _ = self.fan_store([])
        vertices, triangles = fan()
        valid = np.array([False, True, True, True, True])
        nan = error_summary(store, vertices, triangles, np.array([np.nan, 20, 30, 40, 50.0]), valid)
        got = error_summary(store, vertices, triangles, np.array([10.0, 20, 30, 40, 50]), valid)
        assert (got.count, got.max_error, got.rms_error) == (
            nan.count,
            nan.max_error,
            nan.rms_error,
        )

    def test_a_store_never_frozen_is_refused(self, error_summary: Summary) -> None:
        store = _core.CheckPoints(x_min=X0, y_max=Y0, spacing=1.0, rows=5, cols=5)  # type: ignore[attr-defined]
        store.add(world([1.5], [1.5]), np.ones(1, np.float32))
        vertices, triangles = fan()
        with pytest.raises(RuntimeError, match="not frozen"):
            error_summary(store, vertices, triangles, np.zeros(5), np.ones(5, bool))

    def test_an_off_node_corner_skips_no_node_beside_it(self, error_summary: Summary) -> None:
        """The fan with its centre a quarter cell off node (2, 2): the four
        corners are the only vertices on nodes, so 21 nodes count, (2, 2)
        among them; no coincidence radius applies. The centre's height is the
        bilinear one, as `scan` takes it."""
        array = noise(5, 5, 14)
        vertices, triangles = fan()
        vertices[4] = world([2.25], [2.0])[0]
        z = corner_z(array, vertices)
        z[4] = 0.75 * float(array[2, 2]) + 0.25 * float(array[2, 3])
        valid = np.ones(5, bool)
        got = error_summary(view_of(array), vertices, triangles, z, valid)
        want = lattice_oracle(array, vertices, triangles, z, valid)
        assert want.count == 25 - 4
        agrees(got, want)


# ---------------------------------------------------------------- B3: the producer's own error


class TestB3SameAsTheRefinement:
    """`scan` and `scan_points` carry the accumulator, so the summary's largest
    error is the one the refinement reports for the same mesh, bit for bit."""

    def test_the_native_path(self, error_summary: Summary) -> None:
        array = rough(33, 33, 7)
        view, mesh, edges, masks = start(array, 16)
        out = _core.refine(view, mesh, edges, masks, tolerance=3.0)  # type: ignore[attr-defined]
        assert out.ok(), out.message
        v, t = np.asarray(out.vertices), np.asarray(out.triangles)
        z, valid = np.asarray(out.z), np.asarray(out.valid)
        got = error_summary(view, v, t, z, valid)
        assert got.max_error == out.max_error
        # test_core_refine's frame: DX 10, DY 5 from (500000, 7000000).
        frame_v = np.column_stack([(v[:, 0] - 500000.0) / 10.0, (7000000.0 - v[:, 1]) / 5.0])
        r, c = np.indices(array.shape)
        nodes = np.column_stack([c.ravel(), r.ravel()]).astype(np.float64)
        want = summary_oracle(nodes, array.ravel().astype(np.float64), frame_v, t, z, valid)
        agrees(got, want)
        assert want.count == array.size - len(v)

    def test_the_store_path(self, error_summary: Summary) -> None:
        p1 = Phase1(n=9, tolerance=4.0)
        xy, zp = scattered(p1.n, 2, seed=3)
        store = p1.store(_core.CheckPoints, xy, zp)  # type: ignore[attr-defined]
        out = _core.refine_points(store, *p1.args(), tolerance=1.0)  # type: ignore[attr-defined]
        assert out.ok(), out.message
        assert out.coincident == 0
        v, t = np.asarray(out.vertices), np.asarray(out.triangles)
        z, valid = np.asarray(out.z), np.asarray(out.valid)
        got = error_summary(store, v, t, z, valid)
        assert got.max_error == out.max_error
        raw = np.column_stack([(v[:, 0] - PX) / PH, (PY - v[:, 1]) / PH])
        frame_v = np.rint(raw * SCALE) / SCALE  # inserted vertices: stored positions
        assert_allclose(frame_v, raw, rtol=0, atol=1e-6)
        frame_p = np.column_stack([(xy[:, 0] - PX) / PH, (PY - xy[:, 1]) / PH])
        want = summary_oracle(frame_p, zp.astype(np.float64), frame_v, t, z, valid)
        agrees(got, want)
        assert 0 < want.count < len(xy)


# ---------------------------------------------------------------- determinism and refusals


class TestThreads:
    @pytest.mark.parametrize("threads", [2, 3, 7, 0])
    def test_bit_identical_for_any_thread_count_on_the_lattice(
        self, error_summary: Summary, threads: int
    ) -> None:
        view, mesh, edges, masks = start(rough(65, 65, 9), 16)
        out = _core.refine(view, mesh, edges, masks, tolerance=2.0)  # type: ignore[attr-defined]
        arrays = [np.asarray(a) for a in (out.vertices, out.triangles, out.z, out.valid)]
        ref = error_summary(view, *arrays, threads=1)
        got = error_summary(view, *arrays, threads=threads)
        assert (got.count, got.max_error, got.rms_error) == (
            ref.count,
            ref.max_error,
            ref.rms_error,
        )

    @pytest.mark.parametrize("threads", [2, 3, 7, 0])
    def test_bit_identical_for_any_thread_count_in_the_store(
        self, error_summary: Summary, threads: int
    ) -> None:
        p1 = Phase1(n=33, tolerance=4.0)
        store = p1.store(_core.CheckPoints, *scattered(p1.n, 3, seed=4))  # type: ignore[attr-defined]
        out = _core.refine_points(store, *p1.args(), tolerance=2.0)  # type: ignore[attr-defined]
        arrays = [np.asarray(a) for a in (out.vertices, out.triangles, out.z, out.valid)]
        ref = error_summary(store, *arrays, threads=1)
        got = error_summary(store, *arrays, threads=threads)
        assert (got.count, got.max_error, got.rms_error) == (
            ref.count,
            ref.max_error,
            ref.rms_error,
        )

    def test_the_summary_releases_the_gil(self, error_summary: Summary) -> None:
        """One thread over 16 M nodes in two triangles, summarised `calls` times in
        a row: long enough to watch. One call took 0.048 s on a Linux CI runner,
        so `calls = 6` keeps the timed work near 0.3 s there, six times the bound;
        repeating the call costs no memory, unlike a larger grid."""
        n, calls = 4001, 6
        array = np.random.default_rng(1).uniform(0.0, 100.0, (n, n)).astype(np.float32)
        vertices, triangles, z = corner_mesh(array)
        view = view_of(array)
        valid = np.ones(4, bool)
        results, ticks, elapsed = ticks_during(
            lambda: [
                error_summary(view, vertices, triangles, z, valid, threads=1) for _ in range(calls)
            ]
        )
        assert [r.count for r in results] == [n * n - 4] * calls
        assert elapsed >= 0.05, (
            f"{calls} error_summary calls took only {elapsed:.3f}s -- too fast to measure "
            "GIL release; add calls rather than lowering this bound"
        )
        assert ticks >= RELEASED_TICKS, (
            f"only {ticks} ticks in {elapsed:.3f}s: error_summary appears to hold the GIL"
        )


class TestRefusals:
    """Nothing to report a status in, so a vertex the frame cannot hold is a
    ValueError naming it, as the lattice is built from the vertices."""

    @pytest.mark.parametrize(
        "corner", [(X0 - 0.5, Y0), (X0 + 3.0, Y0 - 2.5), (math.nan, Y0), (X0, math.inf)]
    )
    def test_a_vertex_outside_the_node_rectangle(
        self, error_summary: Summary, corner: tuple[float, float]
    ) -> None:
        array = plane_tile()
        vertices, triangles, z = corner_mesh(array)
        vertices[0] = corner
        with pytest.raises(ValueError, match=r"outside"):
            error_summary(view_of(array), vertices, triangles, z, np.ones(4, bool))


class TestStubs:
    """`_core.pyi` is hand-written; mypy checks consistency, not completeness."""

    def declared(self) -> dict[str, ast.ClassDef | ast.FunctionDef]:
        tree = ast.parse(STUB.read_text(encoding="utf-8"))
        return {
            node.name: node
            for node in tree.body
            if isinstance(node, ast.ClassDef | ast.FunctionDef)
        }

    def test_the_stub_declares_the_summary(self) -> None:
        cls = self.declared().get("ErrorSummary")
        assert isinstance(cls, ast.ClassDef), "ErrorSummary is not stubbed"
        members = {node.name for node in cls.body if isinstance(node, ast.FunctionDef)}
        assert {"max_error", "rms_error", "count"} <= members

    def test_error_summary_takes_threads_by_keyword(self) -> None:
        tree = ast.parse(STUB.read_text(encoding="utf-8"))
        found = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and node.name == "error_summary"
        ]
        assert found, "error_summary is not stubbed"
        for fn in found:  # one stub per overload
            assert [a.arg for a in fn.args.kwonlyargs] == ["threads"]
            assert [a.arg for a in fn.args.args][1:] == ["vertices", "triangles", "z", "valid"]
