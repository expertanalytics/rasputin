"""A tolerance that follows the slope: the binding.

Increment 34 (``docs/increments/34-slope-tolerance.md``, sections 4.1, 4.5 and
9). ``_core.SlopeTolerance(view, N, F, START, END, threads)`` computes every
node's slope class of ``view`` (Horn's, in half degrees rounded up, 255 for
NoData) and the tolerance each class allows; ``refine``, ``refine_strip`` and
``refine_points`` take it as ``slope=``, last by keyword, after 33's
``field``.

- test 10: ``histogram()`` sums to the node count and matches this file's own
  classes on V1 (section 3's Horn, written here in NumPy); ``slope=None`` is
  today's call.
- the binding's refusals: a bad ramp names its bound (test 9's words); a slope
  built on another raster geometry is refused (section 4.5).
- test 6 on the check points as ``scattered(65, 4, 7)`` draws them
  (``test_core_refine_points.py``): positions and noise in that draw order, z
  V1's bilinear surface plus the noise (README item 4b: 129 mixed cells, 253
  points nearest a corner below 30 degrees, 136 of those with noise above
  2 m). Every point ends within ``t`` of its cell's largest valid corner
  class; on V1n, the 16 points around the NoData node within N = 2.

CHOSEN HERE, where section 4.5 is silent: ``histogram()`` returns 256 counts
as a sequence of ints; ``RefineOutcome.max_slope_share`` is a float
property.
"""

from __future__ import annotations

import ast
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_array_equal

from test_core_refine_points import STUB
from tin_engine import _core
from tin_engine._core import ChainRole
from tin_engine.cli import DEFAULT_SNAP_SPACING, _constraint_arrays, _engine

X0, Y0, D = 500_000.0, 6_600_000.0, 10.0
N = 65
NODATA = -9999.0
#: V1n's NoData node: row 32, column 30 (x = 300 m), on the wall.
HOLE = (32, 30)

F64 = npt.NDArray[np.float64]


# ---------------------------------------------------------------- fixtures


def valley(rows: int = N, cols: int = N, dx: float = D, dy: float = D) -> npt.NDArray[np.float32]:
    """Section 9's valley, as ``34-probes/fixture_figures.py`` builds it: a
    flat floor (x below 200 m), a 40-degree wall to x = 440 m, a plateau, and
    3 m bumps; computed in float64, stored as float32."""
    r, c = np.indices((rows, cols), dtype=np.float64)
    x, y = c * dx, r * dy
    base = np.tan(np.radians(40.0)) * np.clip(x - 200.0, 0.0, 240.0)
    bumps = 3.0 * np.sin(2 * np.pi * y / 80.0) * np.sin(2 * np.pi * x / 130.0)
    return np.ascontiguousarray((base + bumps).astype(np.float32))


def v1n() -> npt.NDArray[np.float32]:
    z = valley()
    z[HOLE] = NODATA
    return z


def view_of(z: npt.NDArray[np.float32], x0: float = X0) -> Any:
    return _core.raster_view(z, x_min=x0, y_max=Y0, delta_x=D, delta_y=D, nodata=NODATA)


def horn(z: npt.NDArray[np.float32], dx: float = D, dy: float = D) -> F64:
    """Section 3's slope in degrees, NaN at NoData: Horn's 3 by 3 differences,
    a missing neighbour (outside, or NoData) filled by its reflection through
    the node when the opposite one is there, else an edge one by z(node) and a
    corner one by z(row neighbour) + z(column neighbour) - z(node)."""
    zz = z.astype(np.float64)
    ok = zz != NODATA
    rows, cols = zz.shape
    p = np.pad(np.where(ok, zz, np.nan), 1, constant_values=np.nan)

    def at(di: int, dj: int) -> F64:
        return p[1 + di : 1 + di + rows, 1 + dj : 1 + dj + cols]

    def edge(di: int, dj: int) -> F64:
        v, o = at(di, dj), at(-di, -dj)
        return np.where(np.isnan(v), np.where(np.isnan(o), zz, 2.0 * zz - o), v)

    def corner(di: int, dj: int) -> F64:
        v, o = at(di, dj), at(-di, -dj)
        fill = edge(di, 0) + edge(0, dj) - zz
        return np.where(np.isnan(v), np.where(np.isnan(o), fill, 2.0 * zz - o), v)

    a, b, c = corner(-1, -1), edge(-1, 0), corner(-1, 1)
    d, f = edge(0, -1), edge(0, 1)
    g, h, i = corner(1, -1), edge(1, 0), corner(1, 1)
    gx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * dx)
    gy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * dy)
    return np.where(ok, np.degrees(np.arctan(np.hypot(gx, gy))), np.nan)


def classes(z: npt.NDArray[np.float32]) -> npt.NDArray[np.int64]:
    s = horn(z)
    return np.where(np.isnan(s), 255, np.ceil(2.0 * np.nan_to_num(s))).astype(np.int64)


def near_boundary(z: npt.NDArray[np.float32]) -> npt.NDArray[np.bool_]:
    """Within 1e-9 degrees of a class boundary (test 2's window); 0 excluded."""
    s = np.nan_to_num(horn(z))
    near: npt.NDArray[np.bool_] = (s != 0.0) & (np.abs(2.0 * s - np.round(2.0 * s)) <= 2e-9)
    return near


def ramp(cls: npt.NDArray[np.int64], near: float, far: float, start: float, end: float) -> F64:
    """Section 3's t(class / 2): N from END up, F up to START, linear between."""
    s = cls / 2.0
    t = np.where(s >= end, near, far).astype(np.float64)
    if end > start:
        mid = (s > start) & (s < end)
        t[mid] = far + (near - far) * (s[mid] - start) / (end - start)
    return np.where(cls > 180, far, t)


def cell_classes(cls: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
    """Each cell's largest valid corner class, 0 when none is valid."""
    v = np.where(cls == 255, 0, cls)
    top: npt.NDArray[np.int64] = np.maximum.reduce([v[:-1, :-1], v[:-1, 1:], v[1:, :-1], v[1:, 1:]])
    return top


def outline_start(view: Any) -> tuple[Any, npt.NDArray[np.uint32], npt.NDArray[np.uint32]]:
    """The node rectangle's outline as the start mesh, with its constraint arrays."""
    x1, y1 = X0 + (N - 1) * D, Y0 - (N - 1) * D
    xy = np.array([(X0, y1), (x1, y1), (x1, Y0), (X0, Y0)])
    built = _engine(xy, [([0, 1, 2, 3], ChainRole.Outer, 0)], True, DEFAULT_SNAP_SPACING)
    assert built.mesh is not None and built.noded is not None, built.message
    edges, masks = _constraint_arrays(built.mesh, built.noded)
    return built.mesh, edges, masks


def drawn(seed: int = 7, per_cell: int = 4) -> tuple[F64, F64, F64]:
    """``scattered(65, 4, 7)``'s draw, in its order: (col, row, noise)."""
    rng = np.random.default_rng(seed)
    cells = (N - 1) * (N - 1) * per_cell
    base_r, base_c = np.divmod(np.repeat(np.arange((N - 1) * (N - 1)), per_cell), N - 1)
    col = base_c + rng.integers(1, 1024, cells) / 1024.0
    row = base_r + rng.integers(1, 1024, cells) / 1024.0
    noise = rng.integers(-256, 257, cells) / 64.0
    return col, row, noise


def bilinear(z: npt.NDArray[np.float32], col: F64, row: F64) -> F64:
    zz = z.astype(np.float64)
    c0 = np.minimum(np.floor(col).astype(int), N - 2)
    r0 = np.minimum(np.floor(row).astype(int), N - 2)
    tx, ty = col - c0, row - r0
    out: F64 = (
        zz[r0, c0] * (1 - tx) * (1 - ty)
        + zz[r0, c0 + 1] * tx * (1 - ty)
        + zz[r0 + 1, c0] * (1 - tx) * ty
        + zz[r0 + 1, c0 + 1] * tx * ty
    )
    return out


def check_points() -> tuple[F64, npt.NDArray[np.float32], F64, F64]:
    """(xy, z, col, row): z V1's bilinear surface plus the draw's noise,
    rounded to 1/64 m so it is a float32 exactly."""
    col, row, noise = drawn()
    z = np.round((bilinear(valley(), col, row) + noise) * 64.0) / 64.0
    xy = np.column_stack([X0 + col * D, Y0 - row * D])
    return xy, z.astype(np.float32), col, row


def store(xy: F64, z: npt.NDArray[np.float32]) -> Any:
    cp = _core.CheckPoints(x_min=X0, y_max=Y0, spacing=D, rows=N, cols=N)
    cp.add(xy, z)
    cp.freeze()
    return cp


def excess(out: Any, xy: F64, zp: npt.NDArray[np.float32], allowed: F64) -> tuple[float, int, int]:
    """G4 by brute force: the largest ``|plane - z_p| - allowed`` over every
    (point, triangle) pair whose closed triangle (three valid vertices) holds
    the point, the plane from the output's z; points at an output vertex are
    skipped. Returns (worst excess, points compared, points in no triangle)."""
    vertices, tri = np.asarray(out.vertices), np.asarray(out.triangles)
    z, valid = np.asarray(out.z), np.asarray(out.valid).astype(bool)
    at_vertex = (xy[:, None, :] == vertices[None, :, :]).all(axis=2).any(axis=1)
    p = xy[~at_vertex] - [X0, Y0]
    pz, pa = zp[~at_vertex].astype(np.float64), allowed[~at_vertex]
    v = vertices - [X0, Y0]
    worst, held = -np.inf, np.zeros(len(p), bool)
    for a, b, c in tri[valid[tri].all(axis=1)]:
        (ax, ay), (bx, by), (cx, cy) = v[a], v[b], v[c]
        area = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        wa = ((bx - p[:, 0]) * (cy - p[:, 1]) - (by - p[:, 1]) * (cx - p[:, 0])) / area
        wb = ((cx - p[:, 0]) * (ay - p[:, 1]) - (cy - p[:, 1]) * (ax - p[:, 0])) / area
        wc = 1.0 - wa - wb
        inside = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
        if not inside.any():
            continue
        held |= inside
        plane = wa * z[a] + wb * z[b] + wc * z[c]
        worst = max(worst, float((np.abs(plane - pz) - pa)[inside].max()))
    return worst, len(p), int((~held).sum())


def same_mesh(a: Any, b: Any) -> None:
    for name in ("vertices", "triangles", "z", "valid", "edges", "masks"):
        assert_array_equal(np.asarray(getattr(a, name)), np.asarray(getattr(b, name)), name)


# ---------------------------------------------------------------- test 10


class TestHistogram:
    def test_it_sums_to_the_node_count_and_matches_the_classes_here(self) -> None:
        z = valley()
        s = _core.SlopeTolerance(view_of(z), 2.0, 10.0, 30.0, 30.0, 1)
        h = np.asarray(s.histogram(), dtype=np.int64)
        assert h.shape == (256,)
        assert int(h.sum()) == N * N
        assert int(h[255]) == 0
        want = np.bincount(classes(z).ravel(), minlength=256)
        # Equal, except that a node within 1e-9 degrees of a class boundary may
        # sit one class above (test 2's window): never lower.
        loose = np.bincount(classes(z)[near_boundary(z)], minlength=256)
        below_got, below_want = np.cumsum(h), np.cumsum(want)
        assert (below_got <= below_want).all()
        assert (below_want - below_got <= np.cumsum(loose)).all()
        assert int(h[60:181].sum()) == 1496  # section 9: held to 2 m by the step at 30

    def test_nodata_is_counted_in_entry_255(self) -> None:
        s = _core.SlopeTolerance(view_of(v1n()), 2.0, 10.0, 30.0, 30.0, 1)
        h = np.asarray(s.histogram(), dtype=np.int64)
        assert int(h[255]) == 1
        assert int(h.sum()) == N * N
        assert int(h[60:181].sum()) == 1495

    def test_one_and_eight_threads_agree(self) -> None:
        view = view_of(valley())
        one = _core.SlopeTolerance(view, 2.0, 10.0, 25.0, 35.0, 1)
        eight = _core.SlopeTolerance(view, 2.0, 10.0, 25.0, 35.0, 8)
        assert list(one.histogram()) == list(eight.histogram())


class TestCalls:
    def test_slope_none_is_todays_call(self) -> None:
        z = valley()
        view = view_of(z)
        mesh, edges, masks = outline_start(view)
        plain = _core.refine(view, mesh, edges, masks, tolerance=10.0)
        same_mesh(plain, _core.refine(view, mesh, edges, masks, tolerance=10.0, slope=None))

    def test_refine_takes_a_slope_and_reports_its_share(self) -> None:
        z = valley()
        view = view_of(z)
        mesh, edges, masks = outline_start(view)
        s = _core.SlopeTolerance(view, 2.0, 10.0, 30.0, 30.0, 1)
        out = _core.refine(view, mesh, edges, masks, tolerance=10.0, slope=s)
        assert out.ok(), out.message
        assert isinstance(out.max_slope_share, float)
        assert 0.0 < out.max_slope_share <= 1.0 + 1e-12
        plain = _core.refine(view, mesh, edges, masks, tolerance=10.0)
        assert len(np.asarray(out.triangles)) > len(np.asarray(plain.triangles))
        assert plain.max_slope_share == 0.0

    @pytest.mark.parametrize(
        ("ramp_args", "word"),
        [
            ((0.0, 10.0, 25.0, 35.0), "near"),
            ((12.0, 10.0, 25.0, 35.0), "near"),
            ((2.0, float("nan"), 25.0, 35.0), "far"),
            ((2.0, 10.0, -1.0, 35.0), "start"),
            ((2.0, 10.0, 36.0, 35.0), "start"),
            ((2.0, 10.0, 25.0, 90.0), "end"),
        ],
    )
    def test_a_bad_ramp_is_a_value_error_naming_the_bound(
        self, ramp_args: tuple[float, float, float, float], word: str
    ) -> None:
        with pytest.raises(ValueError, match=word):
            _core.SlopeTolerance(view_of(valley()), *ramp_args, 1)


class TestOtherGeometry:
    """Section 4.5: classes are per node of the grid they were built on, so a
    call on another raster geometry is refused, not silently wrong."""

    WORDS = "the slope tolerance was built on another raster geometry"

    def test_refine_on_another_view(self) -> None:
        z = valley()
        view = view_of(z)
        mesh, edges, masks = outline_start(view)
        s = _core.SlopeTolerance(view, 2.0, 10.0, 30.0, 30.0, 1)
        wider = np.ascontiguousarray(np.hstack([z[:, :1], z]))
        other = _core.raster_view(wider, x_min=X0 - D, y_max=Y0, delta_x=D, delta_y=D)
        with pytest.raises(ValueError, match=self.WORDS):
            _core.refine(other, mesh, edges, masks, tolerance=10.0, slope=s)

    def test_refine_points_with_a_store_on_another_geometry(self) -> None:
        z = valley()
        view = view_of(z)
        mesh, edges, masks = outline_start(view)
        phase1 = _core.refine(view, mesh, edges, masks, tolerance=10.0)
        s = _core.SlopeTolerance(view, 2.0, 10.0, 30.0, 30.0, 1)
        xy, zp, _, _ = check_points()
        cp = _core.CheckPoints(x_min=X0 - D, y_max=Y0, spacing=D, rows=N, cols=N + 1)
        cp.add(xy, zp)
        cp.freeze()
        args = (phase1.vertices, phase1.triangles, phase1.z, phase1.valid, phase1.edges)
        with pytest.raises(ValueError, match=self.WORDS):
            _core.refine_points(cp, *args, phase1.masks, tolerance=10.0, slope=s)


# ---------------------------------------------------------------- test 6, scattered(65, 4, 7)


class TestCheckPoints:
    """33's resampled-path setup on V1, through the binding: phase 1 with the
    slope, then the final check with it (N = 2, F = 10, a step at 30)."""

    def run(self, z: npt.NDArray[np.float32], s: Any, slope: bool = True) -> tuple[Any, F64]:
        view = view_of(z)
        mesh, edges, masks = outline_start(view)
        kw = {"slope": s} if slope else {}
        phase1 = _core.refine(view, mesh, edges, masks, tolerance=10.0, **kw)
        assert phase1.ok(), phase1.message
        xy, zp, _, _ = check_points()
        o = phase1
        out = _core.refine_points(
            store(xy, zp), o.vertices, o.triangles, o.z, o.valid, o.edges, o.masks,
            tolerance=10.0, **kw,
        )  # fmt: skip
        assert out.ok(), out.message
        return out, xy

    def test_the_draw_is_the_designs(self) -> None:
        """README item 4b's figures, from this draw: 129 cells straddle 30
        degrees; 253 points in them are nearest a corner below it; 136 of those
        carry noise above 2 m."""
        cls = classes(valley())
        col, row, noise = drawn()
        lo = np.minimum.reduce([cls[:-1, :-1], cls[:-1, 1:], cls[1:, :-1], cls[1:, 1:]])
        hi = np.maximum.reduce([cls[:-1, :-1], cls[:-1, 1:], cls[1:, :-1], cls[1:, 1:]])
        mixed = (hi >= 60) & (lo < 60)
        r0, c0 = np.floor(row).astype(int), np.floor(col).astype(int)
        exposed = mixed[r0, c0] & (cls[np.round(row).astype(int), np.round(col).astype(int)] < 60)
        assert int(mixed.sum()) == 129
        assert int(exposed.sum()) == 253
        assert int((exposed & (np.abs(noise) > 2.0)).sum()) == 136

    def test_every_point_within_its_cells_t(self) -> None:
        z = valley()
        s = _core.SlopeTolerance(view_of(z), 2.0, 10.0, 30.0, 30.0, 1)
        out, xy = self.run(z, s)
        _, zp, col, row = check_points()
        cell = cell_classes(classes(z))
        own = cell[np.floor(row).astype(int), np.floor(col).astype(int)]
        allowed = ramp(own, 2.0, 10.0, 30.0, 30.0)
        worst, compared, unheld = excess(out, xy, zp, allowed)
        assert compared > 0
        assert unheld == 0
        assert worst <= 1e-9 * 250.0  # |z| <= 205 m on V1
        assert 0.0 < out.max_slope_share <= 1.0 + 1e-12
        # The check can fail here: the final check without the slope breaks it.
        plain, _ = self.run(z, s, slope=False)
        assert excess(plain, xy, zp, allowed)[0] > 0.0

    def test_v1n_the_points_around_the_nodata_node_end_within_n(self) -> None:
        """README item 4b: four cells with the NoData node as a corner, largest
        valid corner classes 83, 86, 83, 81; 4 points each, 8 of the 16 with
        noise above 2 m. z is V1's bilinear surface there (the NoData node's
        V1 value standing in) plus the noise."""
        z = v1n()
        s = _core.SlopeTolerance(view_of(z), 2.0, 10.0, 30.0, 30.0, 1)
        out, xy = self.run(z, s)
        _, zp, col, row = check_points()
        r0, c0 = np.floor(row).astype(int), np.floor(col).astype(int)
        around = np.isin(r0, (31, 32)) & np.isin(c0, (29, 30))
        assert int(around.sum()) == 16
        cell = cell_classes(classes(z))
        assert [int(cell[r, c]) for r, c in ((31, 29), (31, 30), (32, 29), (32, 30))] == [
            83, 86, 83, 81,
        ]  # fmt: skip
        worst, compared, unheld = excess(out, xy[around], zp[around], np.full(16, 2.0))
        assert unheld == 0
        assert worst <= 1e-9 * 250.0, (worst, compared)


# ---------------------------------------------------------------- the stub


class TestStub:
    """``_core.pyi`` declares what the binding adds (mypy reads the stub)."""

    @staticmethod
    def declared() -> dict[str, ast.stmt]:
        tree = ast.parse(STUB.read_text(encoding="utf-8"))
        return {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef | ast.ClassDef)}

    def test_slope_tolerance_is_declared_with_histogram(self) -> None:
        cls = self.declared().get("SlopeTolerance")
        assert isinstance(cls, ast.ClassDef)
        assert "histogram" in {n.name for n in cls.body if isinstance(n, ast.FunctionDef)}

    @pytest.mark.parametrize("name", ["refine", "refine_strip", "refine_points"])
    def test_the_entry_points_take_slope_last_after_field(self, name: str) -> None:
        fn = self.declared().get(name)
        assert isinstance(fn, ast.FunctionDef)
        assert [a.arg for a in fn.args.kwonlyargs][-2:] == ["field", "slope"]

    def test_the_outcome_carries_max_slope_share(self) -> None:
        cls = self.declared().get("RefineOutcome")
        assert isinstance(cls, ast.ClassDef)
        assert "max_slope_share" in {n.name for n in cls.body if isinstance(n, ast.FunctionDef)}
