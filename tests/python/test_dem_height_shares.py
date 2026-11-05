"""Each triangle's DEM height shares through `_core.dem_height_shares`.

`docs/increments/137-dem-stats.md`, section 5: tests 4 to 8, 19, the binding
halves of tests 1, 2, 13 and 17, and invariants 2, 8 and 9 at the binding,
with the refusals of 3.4. Every scene is in the DEM's index space: pixel
(r, c) is [c, c + 1] x [r, r + 1] and node (r, c) its centre. The figures
are the design's (an exact `Fraction` oracle); the scenes without a hand
figure are checked against `height_fixtures.oracle_table` (GEOS's clip).
`class_m=1.0, half=10` unless said.

Bound on a share: 1e-7 absolute, the float32 rounding of a value at most 1
(6e-8) with room; a share the design gives as 0 is asserted exactly 0, and
first and count exactly.

INVARIANT-CRITICAL (section 5): the overlap weights, the residual and the
span. The kill record over section 5's mutation targets is
`docs/increments/137-probes/kill-record.md`.

PINNED BY THIS SUITE (the design names the call; the spelling here):
`_core.dem_height_shares(corners, z, dem, nodata, row_area, *, class_m,
half, threads=0) -> (offsets, first, shares, nodata)`; `dem` float32 or
float64 2-D (another dtype is a `TypeError` naming `dem`, as `class_areas`
names `classes`); every refusal of 3.4 is a `ValueError`, the one past the
window saying `outside`, a zero-area one `zero area`, a non-finite one
`finite`, and one of `class_m` or `half` naming it.
"""

from __future__ import annotations

import ast
import math
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from gil_probe import silence_during
from height_fixtures import Table, oracle_table, row_sums, table_of
from test_core_cdt import RELEASED_SILENCE
from tin_engine import _core

STUB = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "_core.pyi"

#: On a share: float32 rounding of a value in [0, 1] is at most 6e-8.
SHARE = 1e-7
HALF = 10


def call(
    corners: Any,
    z: Any,
    dem: Any,
    nodata: float | None = None,
    row_area: Any = None,
    *,
    class_m: float = 1.0,
    half: int = HALF,
    threads: int = 0,
) -> tuple[np.ndarray, ...]:
    grid = dem if isinstance(dem, np.ndarray) else np.asarray(dem, dtype=np.float64)
    weights = np.ones(grid.shape[0]) if row_area is None else np.asarray(row_area, np.float64)
    out = _core.dem_height_shares(
        np.asarray(corners, dtype=np.float64).reshape(-1, 3, 2),
        np.asarray(z, dtype=np.float64).reshape(-1, 3),
        grid,
        nodata,
        weights,
        class_m=class_m,
        half=half,
        threads=threads,
    )
    return tuple(out)


def shares_of(*args: Any, **kw: Any) -> Table:
    offsets, first, shares, nodata = call(*args, **kw)
    assert (offsets.dtype, first.dtype) == (np.uint32, np.int32), (offsets.dtype, first.dtype)
    assert (shares.dtype, nodata.dtype) == (np.float32, np.float32), (shares.dtype, nodata.dtype)
    assert offsets.ndim == first.ndim == shares.ndim == nodata.ndim == 1
    assert len(offsets) == len(first) + 1 == len(nodata) + 1
    assert offsets[0] == 0 and offsets[-1] == len(shares)
    return table_of(offsets, first, shares, nodata)


def expect(
    table: Table, t: int, first: int, count: int, want: dict[int, float], nod: float = 0.0
) -> None:
    """Row `t` against `{k: share}` (k in units of w), the design's figures:
    first and count exactly, listed shares within `SHARE`, the rest exactly 0."""
    assert (int(table.first[t]), int(table.count[t])) == (first, count), (
        table.first[t],
        table.count[t],
    )
    expected = np.zeros(count)
    for k, value in want.items():
        expected[k + HALF - first] = value
    row = table.row(t)
    assert_allclose(row, expected, rtol=0, atol=SHARE)
    assert_array_equal(row[expected == 0.0], 0.0)
    assert table.nodata[t] == pytest.approx(nod, abs=SHARE)
    if nod == 0.0:
        assert table.nodata[t] == 0.0


# ------------------------------------------------- test 1 and 2 at the binding

#: Test 1's DEM: 11 x 11 nodes, 2 m on columns 1 and 2, 0 elsewhere.
STEPS = np.zeros((11, 11))
STEPS[:, 1:3] = 2.0
#: The node rectangle's corners in index space: nw, ne, sw, se.
NW, NE, SW, SE = (0.5, 0.5), (10.5, 0.5), (0.5, 10.5), (10.5, 10.5)


class TestStepsAtTheBinding:
    """Tests 1 and 2 in index space, both diagonals. The corners lie on nodes
    and the diagonal runs through nodes, so the held-node test decides the
    span: the triangle away from the 2 m columns holds their diagonal nodes."""

    @pytest.mark.parametrize(
        ("west", "east"),
        [((NW, SE, SW), (NW, NE, SE)), ((NW, NE, SW), (NE, SE, SW))],
        ids=["nw-se", "ne-sw"],
    )
    def test_1_both_triangles_on_either_diagonal(self, west: Any, east: Any) -> None:
        table = shares_of([west, east], np.zeros((2, 3)), STEPS)
        expect(table, 0, 10, 3, {0: 33 / 50, 2: 17 / 50})
        expect(table, 1, 10, 3, {0: 47 / 50, 2: 3 / 50})

    def test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator(self) -> None:
        dem = STEPS.copy()
        dem[5, 3] = -9999.0
        table = shares_of([(NW, SE, SW), (NW, NE, SE)], np.zeros((2, 3)), dem, -9999.0)
        # 16/25 and 17/50 over all the area; without NoData in it, 32/49.
        expect(table, 0, 10, 3, {0: 16 / 25, 2: 17 / 50}, nod=1 / 50)
        expect(table, 1, 10, 3, {0: 47 / 50, 2: 3 / 50})

    def test_2_a_nan_pixel_is_nodata_without_a_sentinel(self) -> None:
        dem = STEPS.copy()
        dem[5, 3] = np.nan
        table = shares_of([(NW, SE, SW)], np.zeros(3), dem)
        expect(table, 0, 10, 3, {0: 16 / 25, 2: 17 / 50}, nod=1 / 50)

    def test_without_a_sentinel_the_sentinel_value_is_a_height(self) -> None:
        dem = STEPS.copy()
        dem[5, 3] = -9999.0
        table = shares_of([(NW, SE, SW)], np.zeros(3), dem, None)
        # Held, so the span runs to the clamp at -10.
        expect(table, 0, 0, 13, {-10: 1 / 50, 0: 16 / 25, 2: 17 / 50})


# ------------------------------------------------- test 4: class edges and clamp


#: Test 4's DEM: 6 x 6 pixels of 0 with six marked nodes inside the triangle
#: and one outside it.
EDGE_DEM = np.zeros((6, 6))
for (r, c), value in {
    (0, 0): 3.0, (0, 1): 10.5, (0, 2): -10.6, (1, 0): -10.5, (1, 1): -0.5, (2, 0): 0.5,
    (5, 5): 7.0,
}.items():  # fmt: skip
    EDGE_DEM[r, c] = value
#: Area 18; it holds the 21 nodes with r + c <= 5.
EDGE_TRIANGLE = ((0.0, 0.0), (6.0, 0.0), (0.0, 6.0))


class TestClassEdges:
    """Test 4: an edge belongs to the class above it, classes are centred on
    0, residuals past +-(n + 1/2) w are clamped into the end classes."""

    def test_every_edge_and_both_clamps_where_the_design_says(self) -> None:
        table = shares_of([EDGE_TRIANGLE], [0, 0, 0], EDGE_DEM)
        expect(table, 0, 0, 21, {-10: 2 / 18, 0: 13 / 18, 1: 1 / 18, 3: 1 / 18, 10: 1 / 18})

    def test_6_the_plane_runs_through_the_written_z(self) -> None:
        table = shares_of([EDGE_TRIANGLE], [1.0, 1.0, 1.0], EDGE_DEM)
        expect(table, 0, 0, 21, {-10: 2 / 18, -1: 13 / 18, 0: 1 / 18, 2: 1 / 18, 10: 1 / 18})

    def test_the_width_divides_the_residual(self) -> None:
        """Test 20's width at the binding: w = 2.5, n = 4, edges (j + 1/2) 2.5.
        3.0 is class +1 (1.25 <= 3 < 3.75), 10.5 class +4 (clamped from 4.2),
        -10.6 and -10.5 class -4, -0.5 and 0.5 class 0."""
        table = shares_of([EDGE_TRIANGLE], [0, 0, 0], EDGE_DEM, class_m=2.5, half=4)
        assert (table.first[0], table.count[0]) == (0, 9)
        assert_allclose(
            table.row(0), [2 / 18, 0, 0, 0, 14 / 18, 1 / 18, 0, 0, 1 / 18], rtol=0, atol=SHARE
        )

    def test_a_residual_past_n_plus_1_metres_stays_inside_the_classes(self) -> None:
        """The clamp is in units of w: at w = 2.5, n = 4, residuals of +6 and
        -6 m lie past n + 1 = 5 but inside the classes (5 < 6 < 8.75, the end
        classes' edge): classes +2 and -2, not the end classes. Kills a
        quotient taken without the division by w (`class_of`,
        height_shares.hpp:65), whose early clamp then fires at
        floor(6 + 1/2) = 6 > n + 1. That mutant errs on the residuals from
        n + 3/2 m to (n - 1/2) w, either sign, a class not empty once
        w > (n + 3/2) / (n - 1/2): here 5.5 to 8.75 m."""
        dem = np.zeros((6, 6))
        dem[0, 0], dem[0, 1] = 6.0, -6.0
        table = shares_of([EDGE_TRIANGLE], [0, 0, 0], dem, class_m=2.5, half=4)
        assert (table.first[0], table.count[0]) == (2, 5)
        row = table.row(0)
        assert_allclose(row, [1 / 18, 0, 16 / 18, 0, 1 / 18], rtol=0, atol=SHARE)
        assert_array_equal(row[[1, 3]], 0.0)

    def test_a_width_of_1_216_keeps_its_own_edges(self) -> None:
        """ANADEM's and GLO-30's width, n = 4: 0.5 lies in class 0 (below
        0.608), 3.0 in class +2 (1.824 <= 3 < 3.04), 10.5 clamped to +4."""
        table = shares_of([EDGE_TRIANGLE], [0, 0, 0], EDGE_DEM, class_m=1.216, half=4)
        assert (table.first[0], table.count[0]) == (0, 9)
        assert_allclose(
            table.row(0), [2 / 18, 0, 0, 0, 14 / 18, 0, 1 / 18, 0, 1 / 18], rtol=0, atol=SHARE
        )

    def test_a_residual_on_an_edge_whose_quotient_rounds_low(self) -> None:
        """3.4: the class is exact against the stated edges, not the quotient's
        floor. At w = 1.216, (7 + 1/2) w and (-4 + 1/2) w (the lower edges of
        classes +8 and -3, each one product) give r / w + 1/2 just under 8 and
        -3, so the floor alone is one class low; the comparisons raise it."""
        width = 1.216
        dem = np.zeros((6, 6))
        dem[0, 0], dem[0, 1] = (7 + 0.5) * width, (-4 + 0.5) * width
        assert math.floor(dem[0, 0] / width + 0.5) == 7
        assert math.floor(dem[0, 1] / width + 0.5) == -4
        table = shares_of([EDGE_TRIANGLE], [0, 0, 0], dem, class_m=width)
        expect(table, 0, 7, 12, {-3: 1 / 18, 0: 16 / 18, 8: 1 / 18})

    def test_a_pixel_outside_the_triangle_counts_nowhere(self) -> None:
        table = shares_of([EDGE_TRIANGLE], [0, 0, 0], EDGE_DEM)
        assert table.row(0)[HALF + 7] == 0.0


# ------------------------------------------------- test 5: no held node


class TestNoHeldNode:
    """Test 5 and 3.1: a triangle holding no node has one class, that of its
    area-weighted mean residual, with the residual at the pixel centre."""

    def test_5_the_plane_at_the_pixel_centre_not_the_piece_centroid(self) -> None:
        table = shares_of([((1.1, 1.1), (1.4, 1.1), (1.1, 1.4))], [0.0, 3.0, 0.0], np.zeros((3, 3)))
        # Plane 4.0 at (1.5, 1.5): class -4; at the centroid it would be 1.0, class -1.
        expect(table, 0, 6, 1, {-4: 1.0})

    def test_the_mean_residual_not_the_pieces_spread(self) -> None:
        """Pieces of area 0.01, 0.035, 0.035 and 0.045 over pixels of 0, 0, 0
        and 4 m: the mean is 1.44, class +1, though the pieces span 0 to +4."""
        dem = np.zeros((3, 3))
        dem[1, 1] = 4.0
        table = shares_of([((0.9, 0.9), (1.4, 0.9), (0.9, 1.4))], [0, 0, 0], dem)
        expect(table, 0, 11, 1, {1: 1.0})

    def test_the_mean_is_clamped_to_the_end_class(self) -> None:
        dem = np.full((3, 3), 50.0)
        table = shares_of([((1.1, 1.1), (1.4, 1.1), (1.1, 1.4))], [0, 0, 0], dem)
        expect(table, 0, 20, 1, {10: 1.0})

    def test_held_nodes_all_nodata_with_valid_pieces_take_the_mean(self) -> None:
        """Section 6: a held node on NoData joins no span."""
        dem = np.zeros((3, 3))
        dem[0, 1] = dem[1, 0] = 2.0
        dem[1, 1] = np.nan
        # Holds node (1, 1) only, NaN. Area 1.28: pixel (1, 1) 0.68 (NoData),
        # (0, 1) and (1, 0) 0.2 each at 2 m, the rest 0.2 at 0 m. The valid
        # mean is 0.8 / 0.6 = 1.33, class +1, holding 0.6 / 1.28 = 15/32.
        table = shares_of([((0.8, 0.8), (2.4, 0.8), (0.8, 2.4))], [0, 0, 0], dem)
        expect(table, 0, 11, 1, {1: 15 / 32}, nod=17 / 32)


# ------------------------------------------------- 3.2: the held-node margin


class TestHeldNodeMargin:
    """3.2: a node is held when each edge function over its edge's length is
    >= -1e-9. A vertical west edge 1e-10 east of the node column x = 1.5 holds
    those nodes; 1e-8 east of it does not. Node (1, 1) is 5 m, the rest 0;
    pixel (1, 1)'s piece is 0.5 of the triangle's 4.14."""

    @staticmethod
    def triangle(x: float) -> tuple[tuple[float, float], ...]:
        return ((x, 0.2), (3.8, 1.5), (x, 3.8))

    DEM = np.zeros((4, 4))
    DEM[1, 1] = 5.0

    def test_a_node_within_the_margin_is_held(self) -> None:
        table = shares_of([self.triangle(1.5 + 1e-10)], [0, 0, 0], self.DEM)
        expect(table, 0, 10, 6, {0: 1 - 25 / 207, 5: 25 / 207})

    def test_a_node_on_the_edge_is_held(self) -> None:
        table = shares_of([self.triangle(1.5)], [0, 0, 0], self.DEM)
        expect(table, 0, 10, 6, {0: 1 - 25 / 207, 5: 25 / 207})

    def test_a_node_past_the_margin_is_not_held_and_its_piece_folds(self) -> None:
        table = shares_of([self.triangle(1.5 + 1e-8)], [0, 0, 0], self.DEM)
        expect(table, 0, 10, 1, {0: 1.0})


# ------------------------------------------------- test 7: row areas


class TestRowAreas:
    """Test 7: a piece weighs its index area times its row's pixel area."""

    def test_two_triangles_over_rows_of_area_one_and_three(self) -> None:
        dem = np.array([[0.0, 0.0], [2.0, 2.0]])
        corners = [((0, 0), (2, 0), (0, 2)), ((2, 0), (2, 2), (0, 2))]
        table = shares_of(corners, np.zeros((2, 3)), dem, row_area=[1.0, 3.0])
        expect(table, 0, 10, 3, {0: 1 / 2, 2: 1 / 2})
        expect(table, 1, 10, 3, {0: 1 / 10, 2: 9 / 10})


# ------------------------------------------------- test 19: the fold


class TestFold:
    """Test 19 (departure c): a piece of a pixel centred outside the triangle
    folds into the nearer end of the held nodes' span."""

    DEM = np.zeros((5, 5))
    DEM[0, 0], DEM[1, 3] = 2.0, 5.0
    TRIANGLE = ((0.2, 0.2), (4.0, 0.2), (0.2, 4.0))

    def test_19_the_piece_centred_outside_folds_into_class_2(self) -> None:
        table = shares_of([self.TRIANGLE], [0, 0, 0], self.DEM)
        expect(table, 0, 10, 3, {0: 328 / 361, 2: 33 / 361})

    def test_19_a_low_piece_centred_outside_folds_into_the_low_end(self) -> None:
        dem = -self.DEM
        table = shares_of([self.TRIANGLE], [0, 0, 0], dem)
        expect(table, 0, 8, 3, {-2: 33 / 361, 0: 328 / 361})


# ------------------------------------------------- invariant 9 at the binding


class TestOneClassAtHalfZero:
    """Invariant 9 and 3.1's degenerate run: at n = 0 every triangle with
    data has one class, index 0, and still writes its NoData share."""

    def test_test_4s_scene(self) -> None:
        table = shares_of([EDGE_TRIANGLE], [0, 0, 0], EDGE_DEM, half=0)
        assert (table.first[0], table.count[0]) == (0, 1)
        assert table.row(0)[0] == pytest.approx(1.0, abs=SHARE)

    def test_with_nodata(self) -> None:
        dem = STEPS.copy()
        dem[5, 3] = np.nan
        table = shares_of([(NW, SE, SW), (NW, NE, SE)], np.zeros((2, 3)), dem, half=0)
        assert_array_equal(table.first, [0, 0])
        assert_array_equal(table.count, [1, 1])
        assert table.row(0)[0] == pytest.approx(49 / 50, abs=SHARE)
        assert table.nodata[0] == pytest.approx(1 / 50, abs=SHARE)


# ------------------------------------------------- all NoData at the binding


class TestAllNoData:
    """3.3: a triangle whose pixels are all NoData has count 0, first 0 and
    nodata exactly 1."""

    def test_count_0_first_0_nodata_1(self) -> None:
        dem = np.zeros((3, 3))
        dem[1, 1] = -9999.0
        table = shares_of(
            [((1.1, 1.1), (1.4, 1.1), (1.1, 1.4)), ((0.0, 0.0), (3.0, 0.0), (0.0, 3.0))],
            np.zeros((2, 3)),
            dem,
            -9999.0,
        )
        assert (table.first[0], table.count[0], table.nodata[0]) == (0, 0, 1.0)
        assert table.count[1] > 0 and table.offsets[1] == 0


# ------------------------------------------------- test 8: resolution


def wave(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return 3.0 * np.sin(2 * np.pi * x / 400.0) * np.cos(2 * np.pi * y / 300.0)


def wave_scene(h: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """[0, 1000]^2 m at node spacing `h` (node (r, c) at x = h c, y = 1000 - h r),
    the 32 triangles of a 250 m grid split nw-se, corners at the wave."""
    n = round(1000 / h) + 1
    r, c = np.mgrid[0:n, 0:n].astype(np.float64)
    dem = wave(c * h, 1000.0 - r * h)
    g = np.arange(0.0, 1001.0, 250.0)
    world = []
    for x0, x1 in pairwise(g):
        for y0, y1 in pairwise(g):
            nw, ne, sw, se = (x0, y1), (x1, y1), (x0, y0), (x1, y0)
            world += [(nw, se, sw), (nw, ne, se)]
    w = np.array(world)
    corners = np.stack([w[..., 0] / h + 0.5, (1000.0 - w[..., 1]) / h + 0.5], axis=-1)
    return corners, wave(w[..., 0], w[..., 1]), dem, np.full(n, h * h)


def on_classes(table: Table, t: int, half: int) -> dict[int, float]:
    return {int(table.first[t]) + i: float(v) for i, v in enumerate(table.row(t))}


@pytest.fixture(scope="module")
def both() -> tuple[Table, Table]:
    """Test 8's two DEMs of the wave, 10 m and 5 m."""
    coarse, fine = wave_scene(10.0), wave_scene(5.0)
    return (
        shares_of(*coarse[:3], row_area=coarse[3]),
        shares_of(*fine[:3], row_area=fine[3]),
    )


class TestResolution:
    """Test 8: a 10 m and a 5 m DEM of one smooth surface (probe: distance
    median 0.012, largest 0.018)."""

    def test_the_same_count_for_every_triangle(self, both: tuple[Table, Table]) -> None:
        coarse, fine = both
        assert len(coarse.first) == len(fine.first) == 32
        assert_array_equal(coarse.count, fine.count)

    def test_classes_minus_5_to_plus_5_occur(self, both: tuple[Table, Table]) -> None:
        for table in both:
            seen = {k for t in range(32) for k, v in on_classes(table, t, HALF).items() if v > 0}
            assert {k - HALF for k in seen} == set(range(-5, 6))

    def test_the_histograms_lie_within_0_03(self, both: tuple[Table, Table]) -> None:
        coarse, fine = both
        for t in range(32):
            a, b = on_classes(coarse, t, HALF), on_classes(fine, t, HALF)
            keys = set(a) | set(b)
            assert 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys) <= 0.03, t

    @pytest.mark.parametrize("h", [10.0, 5.0])
    def test_against_geos_and_every_row_sums_to_one(self, h: float) -> None:
        corners, z, dem, area = wave_scene(h)
        table = shares_of(corners, z, dem, row_area=area)
        want = oracle_table(corners, z, dem, None, area, class_m=1.0, half=HALF)
        assert_array_equal(table.first, want.first)
        assert_array_equal(table.offsets, want.offsets)
        assert_allclose(table.shares, want.shares, rtol=0, atol=1e-6)
        assert_array_equal(table.nodata, 0.0)
        # Invariant 2: (count + 1) 2^-24.
        assert (np.abs(row_sums(table) - 1.0) <= (table.count + 1) * 2.0**-24).all()


class TestAgainstTheOracle:
    """Random NoData and NaN pixels on the wave, off-node corners, and three
    widths: the producer's table is the oracle's (invariants 2 and 8)."""

    @pytest.mark.parametrize(("class_m", "half"), [(1.0, 10), (1.216, 3), (0.5, 0)])
    def test_the_same_table(self, class_m: float, half: int) -> None:
        rng = np.random.default_rng(1370)
        corners, z, dem, area = wave_scene(10.0)
        corners = corners + rng.uniform(-0.3, 0.3, corners.shape)
        corners = np.clip(corners, 0.0, 101.0)
        dem = dem.copy()
        holes = rng.choice(dem.size, 300, replace=False)
        dem.flat[holes[:150]] = -9999.0
        dem.flat[holes[150:]] = np.nan
        z = z + rng.uniform(-1.0, 1.0, z.shape)
        table = shares_of(corners, z, dem, -9999.0, area, class_m=class_m, half=half)
        want = oracle_table(corners, z, dem, -9999.0, area, class_m=class_m, half=half)
        assert_array_equal(table.first, want.first)
        assert_array_equal(table.offsets, want.offsets)
        assert_allclose(table.shares, want.shares, rtol=0, atol=1e-6)
        assert_allclose(table.nodata, want.nodata, rtol=0, atol=1e-6)
        assert (table.count <= 2 * half + 1).all()
        assert (table.first + table.count <= 2 * half + 1).all()
        assert (np.abs(row_sums(table) - 1.0) <= (table.count + 1) * 2.0**-24).all()


# ------------------------------------------------- test 13: threads


class TestThreads:
    """Test 13 at the binding: rows go to per-block buffers joined in
    triangle order, so the bytes do not depend on the thread count."""

    def test_threads_1_2_7_and_all_give_equal_bytes(self) -> None:
        corners, z, dem, area = wave_scene(5.0)
        corners, z = np.tile(corners, (40, 1, 1)), np.tile(z, (40, 1))
        one = call(corners, z, dem, None, area, threads=1)
        for threads in (2, 7, 0):
            other = call(corners, z, dem, None, area, threads=threads)
            for a, b in zip(one, other, strict=True):
                assert a.tobytes() == b.tobytes(), threads


class TestOrder:
    """3.4: the canonical rotation takes the corners and their z together."""

    def test_any_rotation_or_winding_gives_the_same_bytes(self) -> None:
        corners, z, dem, area = wave_scene(10.0)
        base = call(corners, z, dem, None, area)
        for order in ([1, 2, 0], [2, 0, 1], [0, 2, 1], [2, 1, 0], [1, 0, 2]):
            other = call(corners[:, order], z[:, order], dem, None, area)
            for a, b in zip(base, other, strict=True):
                assert a.tobytes() == b.tobytes(), order


# ------------------------------------------------- refusals (3.4, test 17)


class TestRefusals:
    """Before any work, for the whole batch."""

    INSIDE = ((0.5, 0.5), (1.5, 0.5), (1.0, 1.5))
    DEM = np.zeros((4, 4))

    @pytest.mark.parametrize("class_m", [0.0, -1.0, np.nan, np.inf])
    def test_class_m_not_finite_and_above_0(self, class_m: float) -> None:
        with pytest.raises(ValueError, match="class_m"):
            call([self.INSIDE], [0, 0, 0], self.DEM, class_m=class_m)

    @pytest.mark.parametrize("half", [-1, 128])
    def test_half_outside_0_to_127(self, half: int) -> None:
        with pytest.raises(ValueError, match="half"):
            call([self.INSIDE], [0, 0, 0], self.DEM, half=half)

    def test_half_127_is_taken(self) -> None:
        table = shares_of([self.INSIDE], [0, 0, 0], self.DEM, half=127)
        assert (table.first[0], table.count[0]) == (127, 1)

    def test_a_corner_past_the_window_says_outside(self) -> None:
        with pytest.raises(ValueError, match="outside"):
            call([((3.5, 0.5), (4.5, 0.5), (4.0, 1.5))], [0, 0, 0], self.DEM)
        with pytest.raises(ValueError, match="outside"):
            call([self.INSIDE, ((0.5, -0.1), (1.5, 0.5), (1.0, 1.5))], np.zeros((2, 3)), self.DEM)

    def test_the_window_edge_is_inside(self) -> None:
        table = shares_of([((3.0, 0.5), (4.0, 0.5), (3.5, 1.5))], [0, 0, 0], self.DEM)
        expect(table, 0, 10, 1, {0: 1.0})

    @pytest.mark.parametrize(
        "corners",
        [((0.5, 0.5), (1.5, 1.5), (2.5, 2.5)), ((0.5, 0.5), (0.5, 0.5), (1.0, 1.5))],
        ids=["collinear", "duplicate vertex"],
    )
    def test_a_zero_area_triangle(self, corners: Any) -> None:
        with pytest.raises(ValueError, match="zero area"):
            call([corners], [0, 0, 0], self.DEM)

    @pytest.mark.parametrize("bad", [np.nan, np.inf])
    def test_a_non_finite_corner(self, bad: float) -> None:
        with pytest.raises(ValueError, match="finite"):
            call([((0.5, bad), (1.5, 0.5), (1.0, 1.5))], [0, 0, 0], self.DEM)

    @pytest.mark.parametrize("bad", [np.nan, -np.inf])
    def test_a_non_finite_corner_height(self, bad: float) -> None:
        with pytest.raises(ValueError, match="finite"):
            call([self.INSIDE], [0, bad, 0], self.DEM)

    def test_z_not_three_per_triangle(self) -> None:
        for z in (np.zeros((1, 2)), np.zeros((2, 3))):
            with pytest.raises(ValueError, match="z"):
                _core.dem_height_shares(
                    np.array([self.INSIDE], np.float64),
                    z,
                    self.DEM,
                    None,
                    np.ones(4),
                    class_m=1.0,
                    half=HALF,
                )

    def test_corners_of_the_wrong_shape(self) -> None:
        with pytest.raises(ValueError, match="corners"):
            _core.dem_height_shares(
                np.zeros((1, 3, 3)),
                np.zeros((1, 3)),
                self.DEM,
                None,
                np.ones(4),
                class_m=1.0,
                half=HALF,
            )

    def test_row_area_not_one_per_row(self) -> None:
        with pytest.raises(ValueError, match="row_area"):
            call([self.INSIDE], [0, 0, 0], self.DEM, row_area=np.ones(3))

    def test_a_dem_neither_float32_nor_float64_is_a_type_error(self) -> None:
        with pytest.raises(TypeError, match="dem"):
            _core.dem_height_shares(
                np.array([self.INSIDE], np.float64),
                np.zeros((1, 3)),
                np.zeros((4, 4), np.int16),
                None,
                np.ones(4),
                class_m=1.0,
                half=HALF,
            )

    def test_a_float32_dem_is_taken(self) -> None:
        table = shares_of([EDGE_TRIANGLE], [0, 0, 0], EDGE_DEM.astype(np.float32))
        expect(table, 0, 0, 21, {-10: 2 / 18, 0: 13 / 18, 1: 1 / 18, 3: 1 / 18, 10: 1 / 18})

    def test_no_triangles_give_empty_arrays(self) -> None:
        table = shares_of(np.zeros((0, 3, 2)), np.zeros((0, 3)), self.DEM)
        assert_array_equal(table.offsets, [0])
        assert len(table.first) == len(table.shares) == len(table.nodata) == 0


# ------------------------------------------------- the binding's surface (3.4)


class TestSurface:
    def test_the_gil_is_released(self) -> None:
        corners, z, dem, area = wave_scene(5.0)
        corners, z = np.tile(corners, (200, 1, 1)), np.tile(z, (200, 1))
        _, silence, elapsed = silence_during(
            lambda: [call(corners, z, dem, None, area, threads=1) for _ in range(3)]
        )
        assert elapsed >= 0.05, f"took only {elapsed:.3f}s: add calls rather than lower the bound"
        assert silence <= RELEASED_SILENCE, f"no tick for {silence:.0%}: the GIL appears held"

    def test_the_stub_declares_it_with_keyword_only_class_m_half_and_threads(self) -> None:
        tree = ast.parse(STUB.read_text(encoding="utf-8"))
        found = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and node.name == "dem_height_shares"
        ]
        assert found, "dem_height_shares is not stubbed"
        for fn in found:
            assert [a.arg for a in fn.args.args] == ["corners", "z", "dem", "nodata", "row_area"]
            assert [a.arg for a in fn.args.kwonlyargs] == ["class_m", "half", "threads"]
