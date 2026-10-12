"""The exact overlap through `_core.class_areas`, the cell areas per row and the
`ClassFractions` table (26b).

`docs/increments/26-land-cover.md`, 1.1, 1.2, section 7 F2 and section 10's
"26b": F2 against shapely on random triangulations of a 20 x 20 class grid
(GEOS's clip is the oracle, never ours), the per-cell sums, and threads. The
pinned spellings are in `fraction_fixtures`' docstring.

HOW THIS FILE GOES RED: `_core.class_areas`, `io.models.ClassFractions` and
`land_cover_fractions.row_areas` are reached through `getattr` or an import at
call time (`AttributeError`, `ModuleNotFoundError`), so each test fails on its
own and the rest of `tests/python` collects. `TestWindowCorner` runs on 26a's
code and is green: it pins what 26b's corner arithmetic depends on.
"""

from __future__ import annotations

import importlib
import io
from itertools import pairwise
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from class_fixtures import LAT0, LON0, STEP, class_grid
from cog_fixtures import BytesBlocks
from fraction_fixtures import (
    CELL,
    CLIP,
    cell_cover,
    geod_cell_area,
    geod_row_areas,
    geos_class_areas,
    mixed_grid,
    rows_of,
    shoelace,
    square_cover,
    square_triangulation,
    unique_grid,
)
from geographic_fixtures import ANADEM_STEP, geographic_keys, tiff
from tin_engine.io.models import ClassWindow, IndexWindow

EXTENT = 20
SEEDS = (26, 27, 28)
#: A window at 15 S with the fixture's 1/128 degree cells: about 7.4e5 m2 each.
WINDOW = ClassWindow(
    x_min=LON0, y_max=LAT0, delta_x=STEP, delta_y=STEP, array=mixed_grid(EXTENT, EXTENT)
)


@pytest.fixture(scope="module")
def core() -> ModuleType:
    return importlib.import_module("tin_engine._core")


def class_areas(core: ModuleType, corners: Any, grid: Any, row_area: Any, **kw: Any) -> Any:
    return core.class_areas(
        np.asarray(corners, dtype=np.float64), grid, np.asarray(row_area, dtype=np.float64), **kw
    )


def well_formed(
    offsets: np.ndarray, classes: np.ndarray, areas: np.ndarray, triangles: int
) -> None:
    """1.2's table: T + 1 offsets from 0, classes strictly ascending within a
    row, every area positive; the dtypes the design names."""
    assert (offsets.dtype, classes.dtype, areas.dtype) == (np.uint32, np.uint8, np.float64)
    assert offsets.shape == (triangles + 1,) and offsets[0] == 0
    assert offsets[-1] == len(classes) == len(areas)
    assert (np.diff(offsets.astype(np.int64)) >= 0).all()
    assert (areas > 0).all()
    for a, b in pairwise(offsets):
        assert (np.diff(classes[a:b].astype(np.int64)) > 0).all()


# ------------------------------------------------------------ F2 against GEOS


class TestF2AgainstGeos:
    @pytest.mark.parametrize("seed", SEEDS)
    def test_every_triangles_class_areas_are_geoss_in_cell_units(
        self, core: ModuleType, seed: int
    ) -> None:
        grid = mixed_grid(EXTENT, EXTENT)
        corners = square_triangulation(seed, EXTENT)
        offsets, classes, areas = class_areas(core, corners, grid, np.ones(EXTENT))
        well_formed(offsets, classes, areas, len(corners))
        want = geos_class_areas(corners, grid, np.ones(EXTENT))
        for t, (got, oracle) in enumerate(zip(rows_of(offsets, classes, areas), want, strict=True)):
            assert set(got) == set(oracle), f"triangle {t}"
            for k, area in oracle.items():
                assert got[k] == pytest.approx(area, abs=CLIP), f"triangle {t}, class {k}"

    @pytest.mark.parametrize("seed", SEEDS)
    def test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows(
        self, core: ModuleType, seed: int
    ) -> None:
        grid = mixed_grid(EXTENT, EXTENT)
        corners = square_triangulation(seed, EXTENT)
        weights = geod_row_areas(LAT0, STEP, EXTENT)
        assert weights.min() > 7e5 and np.ptp(weights) > 1.0, "rows must differ"
        offsets, classes, areas = class_areas(core, corners, grid, weights)
        want = geos_class_areas(corners, grid, weights)
        for t, (got, oracle) in enumerate(zip(rows_of(offsets, classes, areas), want, strict=True)):
            assert set(got) == set(oracle), f"triangle {t}"
            for k, area in oracle.items():
                assert got[k] == pytest.approx(area, rel=1e-9), f"triangle {t}, class {k}"

    @pytest.mark.parametrize("seed", SEEDS)
    def test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square(
        self, core: ModuleType, seed: int
    ) -> None:
        """F2: summed over the triangles, each class's area equals its area
        inside the mesh's footprint, to 1e-9 relative. The footprint is the
        square, whose cover per cell is exact by hand."""
        grid = mixed_grid(EXTENT, EXTENT)
        corners = square_triangulation(seed, EXTENT)
        weights = geod_row_areas(LAT0, STEP, EXTENT)
        _, classes, areas = class_areas(core, corners, grid, weights)
        cover = square_cover(EXTENT, EXTENT, EXTENT) * weights[:, None]
        want = {int(k): float(cover[grid == k].sum()) for k in np.unique(grid)}
        want = {k: v for k, v in want.items() if v > 0}
        got: dict[int, float] = {}
        for k, a in zip(classes.tolist(), areas.tolist(), strict=True):
            got[k] = got.get(k, 0.0) + a
        assert set(got) == set(want)
        for k, total in want.items():
            assert got[k] == pytest.approx(total, rel=1e-9), f"class {k}"
        assert 0 in got, "the no-data class is a class like any other (1.1)"

    def test_the_cover_oracle_agrees_with_geos_on_the_square(self) -> None:
        """The oracle can fail: GEOS's cover of the tiling equals the hand cover."""
        corners = square_triangulation(SEEDS[0], EXTENT)
        got = cell_cover(corners, EXTENT, EXTENT)
        assert np.abs(got - square_cover(EXTENT, EXTENT, EXTENT)).max() <= 1e-9


# ------------------------------------------------------------ per-cell sums


class TestPerCellSums:
    @pytest.mark.parametrize("seed", SEEDS)
    def test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area(
        self, core: ModuleType, seed: int
    ) -> None:
        """F2, per cell, to 1e-12 of a cell: with one class per cell the pieces
        of a cell are the entries of its class. Each shared edge's crossings
        are computed once (1.2), so the two sides leave no gap or overlap."""
        side = 15  # 225 cells, each its own class
        grid = unique_grid(side, side)
        corners = square_triangulation(seed, side)
        offsets, classes, areas = class_areas(core, corners, grid, np.ones(side))
        well_formed(offsets, classes, areas, len(corners))
        summed = np.zeros(side * side + 1)
        np.add.at(summed, classes.astype(np.int64), areas)
        cover = square_cover(side, side, side)
        assert np.abs(summed[1:].reshape(side, side) - cover).max() <= CELL
        assert summed[0] == 0.0

    def test_each_triangles_pieces_sum_to_its_area(self, core: ModuleType) -> None:
        grid = mixed_grid(EXTENT, EXTENT)
        corners = square_triangulation(SEEDS[1], EXTENT)
        offsets, _, areas = class_areas(core, corners, grid, np.ones(EXTENT))
        per_triangle = np.add.reduceat(areas, offsets[:-1].astype(np.int64))
        assert np.abs(per_triangle - shoelace(corners)).max() <= CELL


# ------------------------------------------------------------ threads


class TestThreads:
    def test_one_thread_four_threads_and_all_cores_give_the_same_bits(
        self, core: ModuleType
    ) -> None:
        grid = mixed_grid(EXTENT, EXTENT)
        corners = np.concatenate([square_triangulation(s, EXTENT, points=400) for s in SEEDS])
        assert len(corners) > 1000
        weights = geod_row_areas(LAT0, STEP, EXTENT)
        serial = class_areas(core, corners, grid, weights, threads=1)
        for threads in (4, 0):
            parallel = class_areas(core, corners, grid, weights, threads=threads)
            for ours, theirs in zip(serial, parallel, strict=True):
                assert ours.dtype == theirs.dtype
                assert_array_equal(ours, theirs)


# ------------------------------------------------------------ refusals


class TestRefusals:
    GRID = mixed_grid(4, 4)
    INSIDE = (((0.5, 0.5), (1.5, 0.5), (1.0, 1.5)),)

    def test_a_triangle_past_the_window_is_refused_naming_outside(self, core: ModuleType) -> None:
        with pytest.raises(ValueError, match="outside"):
            class_areas(core, [[[3.5, 0.5], [4.5, 0.5], [4.0, 1.5]]], self.GRID, np.ones(4))
        with pytest.raises(ValueError, match="outside"):
            class_areas(core, [[[0.5, -0.1], [1.5, 0.5], [1.0, 1.5]]], self.GRID, np.ones(4))

    def test_the_window_edge_is_inside(self, core: ModuleType) -> None:
        _, _, areas = class_areas(
            core, [[[3.0, 0.5], [4.0, 0.5], [3.5, 1.5]]], self.GRID, np.ones(4)
        )
        assert areas.sum() == pytest.approx(0.5, abs=CELL)

    def test_a_zero_area_triangle_is_refused(self, core: ModuleType) -> None:
        with pytest.raises(ValueError):
            class_areas(core, [[[0.5, 0.5], [1.5, 1.5], [2.5, 2.5]]], self.GRID, np.ones(4))

    def test_a_non_finite_corner_is_refused(self, core: ModuleType) -> None:
        with pytest.raises(ValueError):
            class_areas(core, [[[0.5, np.nan], [1.5, 0.5], [1.0, 1.5]]], self.GRID, np.ones(4))
        with pytest.raises(ValueError):
            class_areas(core, [[[0.5, 0.5], [np.inf, 0.5], [1.0, 1.5]]], self.GRID, np.ones(4))

    def test_another_dtype_or_a_strided_grid_is_a_type_error(self, core: ModuleType) -> None:
        with pytest.raises(TypeError):
            class_areas(core, self.INSIDE, self.GRID.astype(np.uint16), np.ones(4))
        wide = mixed_grid(4, 8)
        with pytest.raises(TypeError):
            class_areas(core, self.INSIDE, wide[:, ::2], np.ones(4))

    def test_wrong_shapes_are_value_errors(self, core: ModuleType) -> None:
        with pytest.raises(ValueError):
            class_areas(core, [[[0.5, 0.5], [1.5, 0.5]]], self.GRID, np.ones(4))
        with pytest.raises(ValueError):
            class_areas(core, self.INSIDE, self.GRID, np.ones(3))
        with pytest.raises(ValueError):
            class_areas(core, self.INSIDE, self.GRID.ravel(), np.ones(4))

    def test_no_triangles_is_one_offset(self, core: ModuleType) -> None:
        offsets, classes, areas = class_areas(core, np.zeros((0, 3, 2)), self.GRID, np.ones(4))
        assert offsets.tolist() == [0] and len(classes) == 0 and len(areas) == 0


# ------------------------------------------------------------ cell areas per row


class TestRowAreas:
    @pytest.fixture(scope="class")
    def fractions(self) -> ModuleType:
        return importlib.import_module("tin_engine.land_cover_fractions")

    def test_each_rows_cell_area_is_the_ellipsoids(self, fractions: ModuleType) -> None:
        got = fractions.row_areas(WINDOW)
        assert got.dtype == np.float64 and got.shape == (EXTENT,)
        want = geod_row_areas(LAT0, STEP, EXTENT)
        assert got == pytest.approx(want, rel=1e-9)

    def test_the_designs_figures_at_thirty_metres(self, fractions: ModuleType) -> None:
        """1.1: 887.5 m2 at 7 S and 836.0 m2 at 21 S for a MapBiomas cell."""
        for lat, want in ((-7.0, 887.5), (-21.0, 836.0)):
            window = ClassWindow(
                x_min=-45.0, y_max=lat, delta_x=ANADEM_STEP, delta_y=ANADEM_STEP,
                array=np.zeros((1, 1), dtype=np.uint8),
            )  # fmt: skip
            assert fractions.row_areas(window)[0] == pytest.approx(want, abs=0.5)
            assert geod_cell_area(lat, ANADEM_STEP) == pytest.approx(want, abs=0.5)

    def test_the_longitude_does_not_matter(self, fractions: ModuleType) -> None:
        east = WINDOW.model_copy(update={"x_min": LON0 + 3.0})
        assert_array_equal(fractions.row_areas(east), fractions.row_areas(WINDOW))


# ------------------------------------------------------------ the table


class TestClassFractions:
    OFFSETS = np.array([0, 2, 3, 5], dtype=np.uint32)
    CLASSES = np.array([3, 15, 33, 0, 24], dtype=np.uint8)
    FRACTIONS = np.array([0.25, 0.75, 1.0, 0.4, 0.6], dtype=np.float32)

    def make(self, **changes: Any) -> Any:
        from tin_engine.io.models import ClassFractions

        fields: dict[str, Any] = {
            "offsets": self.OFFSETS, "classes": self.CLASSES, "fractions": self.FRACTIONS,
        }  # fmt: skip
        return ClassFractions(**(fields | changes))

    def test_a_table_keeps_its_arrays_and_counts_its_triangles(self) -> None:
        table = self.make()
        assert table.triangles == 3
        assert_array_equal(table.offsets, self.OFFSETS)
        assert_array_equal(table.classes, self.CLASSES)
        assert_array_equal(table.fractions, self.FRACTIONS)
        assert (table.offsets.dtype, table.classes.dtype, table.fractions.dtype) == (
            np.uint32, np.uint8, np.float32,
        )  # fmt: skip

    def test_it_is_frozen(self) -> None:
        table = self.make()
        with pytest.raises(ValueError):
            table.triangles_ = 1  # type: ignore[attr-defined]
        with pytest.raises(ValueError):
            table.offsets = self.OFFSETS  # type: ignore[misc]

    @pytest.mark.parametrize(
        "changes",
        [
            {"offsets": np.array([1, 2, 3, 5], dtype=np.uint32)},
            {"offsets": np.array([0, 3, 2, 5], dtype=np.uint32)},
            {"offsets": np.array([0, 2, 3, 4], dtype=np.uint32)},
            {"offsets": np.array([0, 2, 3, 5], dtype=np.int64)},
            {"classes": np.array([3, 15, 33, 0], dtype=np.uint8)},
            {"classes": np.array([3, 15, 33, 0, 24], dtype=np.int32)},
            {"fractions": np.array([0.25, 0.75, 1.0, 0.4, 0.6], dtype=np.float64)},
            {"fractions": np.array([0.25, 0.75, 1.0, 0.4], dtype=np.float32)},
        ],
        ids=[
            "offsets-not-from-0", "offsets-decreasing", "offsets-end-short", "offsets-int64",
            "classes-short", "classes-int32", "fractions-float64", "fractions-short",
        ],
    )  # fmt: skip
    def test_a_malformed_table_is_refused(self, changes: dict[str, Any]) -> None:
        with pytest.raises(ValueError):
            self.make(**changes)

    def test_an_empty_table_has_no_triangles(self) -> None:
        table = self.make(
            offsets=np.array([0], dtype=np.uint32),
            classes=np.zeros(0, dtype=np.uint8),
            fractions=np.zeros(0, dtype=np.float32),
        )
        assert table.triangles == 0


# ------------------------------------------------------------ the window's corner


class TestWindowCorner:
    """26b's corner arithmetic reads `ClassWindow.x_min` and `y_max` as the
    outer corner of the first cell. `decode_class_window` (26a) moves the
    node half a cell out whatever the file's registration, because
    `read_page` puts the node at the cell's centre for both: pinned here, so
    a PixelIsPoint class raster is not half a cell off in silence (26a's code
    review, round 1)."""

    @pytest.mark.parametrize("area", [True, False], ids=["PixelIsArea", "PixelIsPoint"])
    def test_the_windows_corner_is_half_a_cell_outside_the_first_node(self, area: bool) -> None:
        from tin_engine.io.cog import decode_class_window
        from tin_engine.io.geotiff import read_page

        grid = class_grid()[:16, :16]
        data = tiff(
            grid, x0=LON0, y0=LAT0, step_x=STEP, shorts=geographic_keys(area=area), nodata=None,
            tile=(16, 16),
        ).getvalue()  # fmt: skip
        meta = read_page(io.BytesIO(data), nodata=0.0)[0]
        window = decode_class_window(
            BytesBlocks(data), meta, IndexWindow(row0=0, col0=0, rows=16, cols=16)
        )
        node_x, node_y = meta.node_xy(0, 0)
        assert (window.x_min, window.y_max) == (node_x - STEP / 2, node_y + STEP / 2)
        # PixelIsArea: the tie point is the corner itself; PixelIsPoint: the
        # tie point is the node, and the cell around it reaches half a cell out.
        want = (LON0, LAT0) if area else (LON0 - STEP / 2, LAT0 + STEP / 2)
        assert (window.x_min, window.y_max) == want
        assert_array_equal(window.array, grid)
