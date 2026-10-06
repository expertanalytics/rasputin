"""`tin_engine.io.models`: the node arithmetic of a DEM grid, the NoData rule,
and the two value types moved down from layer 1 (audit PR A, `audit-lattice`).

`docs/increments/python-audit-pr-a.md`, "The types" and red test 1:

- `Bounds` (from `mosaic`) and `TileFootprint` (from `io.repository`) live in
  `io.models`; `Bounds.of(box)` takes shapely's `(x_min, y_min, x_max, y_max)`.
- `RasterMeta.node_xy(rows, cols)` is the core's `RasterGeometry::node`,
  `x_min + col * delta_x`, `y_max - row * delta_y`, bit for bit, elementwise,
  with no broadcasting between rows and cols.
- `RasterMeta.index_of(x, y)` is the fractional `(row, col)`, unrounded.
- `RasterMeta.node_box()` is the node rectangle, flat for one row or one node.
- `RasterMeta.windowed(window)` is `io.cog.window_meta`'s expression.
- `valid_mask(values, nodata)`: NoData is NaN or the sentinel; +-inf is data
  (Ola's ruling, `python-audit.md` section 7, ruling 1).

The grid is non-dyadic on purpose (origin 0.3, 0.7, spacing 0.1): a sum of
steps, a reordered expression or a cast moves an ulp there, and every
comparison below is exact (`==`, `tobytes`), so it shows.

HOW THIS FILE GOES RED. `tin_engine.io.models` is reached through a
module-scoped fixture and every new name is read off it, so each missing name
fails its own tests (`AttributeError`) and the rest of `tests/python` still
collects.
"""

from __future__ import annotations

import importlib
import math
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from pydantic import ValidationError

X_MIN, Y_MAX, STEP = 0.3, 0.7, 0.1
ROWS, COLS = 7, 6
SENTINEL = -32767.0


@pytest.fixture(scope="module")
def models() -> ModuleType:
    return importlib.import_module("tin_engine.io.models")


def grid(models: ModuleType, rows: int = ROWS, cols: int = COLS, **fields: Any) -> Any:
    """The non-dyadic grid, `rows x cols` nodes; `fields` override the rest."""
    base: dict[str, Any] = {
        "x_min": X_MIN,
        "y_max": Y_MAX,
        "delta_x": STEP,
        "delta_y": STEP,
        "rows": rows,
        "cols": cols,
        "epsg": 25833,
        "nodata": SENTINEL,
        "nodata_source": "tag",
        "pixel_is_area": False,
        "vertical_unit_assumed": True,
    }
    return models.RasterMeta(**(base | fields))


def same_bits(got: Any, expected: Any) -> bool:
    got, expected = np.asarray(got), np.asarray(expected)
    return (
        got.dtype == expected.dtype
        and got.shape == expected.shape
        and (got.tobytes() == expected.tobytes())
    )


# ------------------------------------------------------------------ the moved types


class TestMovedTypes:
    def test_bounds_and_tile_footprint_live_in_io_models(self, models: ModuleType) -> None:
        mosaic = importlib.import_module("tin_engine.mosaic")
        repository = importlib.import_module("tin_engine.io.repository")
        assert models.Bounds is mosaic.Bounds, "one class, moved, not a copy"
        assert models.TileFootprint is repository.TileFootprint, "one class, moved, not a copy"

    def test_bounds_of_equals_the_keyword_form(self, models: ModuleType) -> None:
        box = (X_MIN, 0.1, 0.8, Y_MAX)
        keywords = models.Bounds(x_min=X_MIN, y_min=0.1, x_max=0.8, y_max=Y_MAX)
        assert models.Bounds.of(box) == keywords

    @pytest.mark.parametrize(
        "box",
        [
            (0.0, 0.0, 0.0, 1.0),
            (0.0, 0.0, 1.0, 0.0),
            (1.0, 0.0, 0.0, 1.0),
            (0.0, 1.0, 1.0, 0.0),
            (math.nan, 0.0, 1.0, 1.0),
            (0.0, 0.0, math.inf, 1.0),
            (0.0, -math.inf, 1.0, 1.0),
        ],
        ids=["flat_x", "flat_y", "inverted_x", "inverted_y", "nan", "+inf", "-inf"],
    )
    def test_bounds_of_refuses_what_the_constructor_refuses(
        self, models: ModuleType, box: tuple[float, float, float, float]
    ) -> None:
        keywords = dict(zip(("x_min", "y_min", "x_max", "y_max"), box, strict=True))
        reason = r"need finite x_min < x_max and y_min < y_max"
        with pytest.raises(ValidationError, match=reason):
            models.Bounds(**keywords)
        with pytest.raises(ValidationError, match=reason):
            models.Bounds.of(box)


# ------------------------------------------------------------------ node_xy and index_of


class TestNodeXY:
    def test_python_ints_bit_for_bit(self, models: ModuleType) -> None:
        m = grid(models)
        for r in range(ROWS):
            for c in range(COLS):
                x, y = m.node_xy(r, c)
                assert (type(x), type(y)) == (float, float)
                assert (x, y) == (X_MIN + c * STEP, Y_MAX - r * STEP), (r, c)

    def test_a_numpy_scalar_stays_one(self, models: ModuleType) -> None:
        """A caller's `np.float64` stays one: a message printing it would change."""
        x, y = grid(models).node_xy(np.int64(3), np.int64(2))
        assert (type(x), type(y)) == (np.float64, np.float64)
        assert (x, y) == (X_MIN + 2 * STEP, Y_MAX - 3 * STEP)

    @pytest.mark.parametrize("dtype", [np.int32, np.int64, np.float64])
    def test_arrays_bit_for_bit(self, models: ModuleType, dtype: Any) -> None:
        r, c = np.indices((ROWS, COLS), dtype=dtype)
        x, y = grid(models).node_xy(r, c)
        assert same_bits(x, X_MIN + c * STEP)
        assert same_bits(y, Y_MAX - r * STEP)

    def test_a_row_band_and_a_column_range_of_different_lengths(self, models: ModuleType) -> None:
        """`x` reads only `cols`, `y` only `rows`: no broadcasting between them."""
        rows, cols = np.arange(2, 7), np.arange(COLS)  # 5 rows, 6 cols
        x, y = grid(models).node_xy(rows, cols)
        assert same_bits(x, X_MIN + cols * STEP)
        assert same_bits(y, Y_MAX - rows * STEP)


class TestIndexOf:
    def test_bit_for_bit_at_nodes_and_between_them(self, models: ModuleType) -> None:
        m = grid(models)
        r, c = np.indices((2 * ROWS - 1, 2 * COLS - 1), dtype=np.float64) / 2.0
        x, y = X_MIN + c * STEP, Y_MAX - r * STEP
        row, col = m.index_of(x, y)
        assert same_bits(row, (Y_MAX - y) / STEP)
        assert same_bits(col, (x - X_MIN) / STEP)

    def test_python_floats(self, models: ModuleType) -> None:
        row, col = grid(models).index_of(0.55, 0.25)
        assert (type(row), type(col)) == (float, float)
        assert (row, col) == ((Y_MAX - 0.25) / STEP, (0.55 - X_MIN) / STEP)

    def test_it_is_unrounded(self, models: ModuleType) -> None:
        """Rounding stays with each caller (three rules, three questions)."""
        row, col = grid(models).index_of(X_MIN + 0.25 * STEP, Y_MAX - 1.75 * STEP)
        assert (row, col) == pytest.approx((1.75, 0.25), abs=1e-12)  # unit-scale grid, < 1

    def test_rounding_it_at_every_node_gives_the_node(self, models: ModuleType) -> None:
        m = grid(models)
        r, c = np.indices((ROWS, COLS))
        row, col = m.index_of(*m.node_xy(r, c))
        assert np.array_equal(np.round(row), r) and np.array_equal(np.round(col), c)


# ------------------------------------------------------------------ node_box and windowed


class TestNodeBox:
    @pytest.mark.parametrize(("rows", "cols"), [(6, 5), (1, 5), (6, 1), (1, 1)])
    def test_the_node_rectangle_bit_for_bit(self, models: ModuleType, rows: int, cols: int) -> None:
        box = grid(models, rows, cols).node_box()
        x_max, y_min = X_MIN + (cols - 1) * STEP, Y_MAX - (rows - 1) * STEP
        assert box == (X_MIN, y_min, x_max, Y_MAX)
        assert all(type(v) is float for v in box)

    def test_a_one_row_or_one_node_grid_is_flat_not_refused(self, models: ModuleType) -> None:
        x_min, y_min, x_max, y_max = grid(models, 1, 5).node_box()
        assert y_min == y_max and x_min < x_max
        x_min, y_min, x_max, y_max = grid(models, 1, 1).node_box()
        assert (x_min, y_min) == (x_max, y_max)


class TestWindowed:
    def test_the_corner_moves_by_whole_cells_and_nothing_else_changes(
        self, models: ModuleType
    ) -> None:
        m = grid(models, crs="EPSG:25833", geographic=False)
        w = models.IndexWindow(row0=3, col0=2, rows=4, cols=3)
        expected = m.model_copy(
            update={
                "x_min": m.x_min + w.col0 * m.delta_x,
                "y_max": m.y_max - w.row0 * m.delta_y,
                "rows": w.rows,
                "cols": w.cols,
            }
        )
        got = m.windowed(w)
        assert got == expected
        assert (got.x_min, got.y_max) == (X_MIN + 2 * STEP, Y_MAX - 3 * STEP)
        assert got.model_dump(exclude={"x_min", "y_max", "rows", "cols"}) == m.model_dump(
            exclude={"x_min", "y_max", "rows", "cols"}
        )


# ------------------------------------------------------------------ the NoData rule


class TestValidMask:
    VALUES = (1.0, math.nan, math.inf, -math.inf, SENTINEL)

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    def test_nan_and_the_sentinel_are_nodata_infinity_is_data(
        self, models: ModuleType, dtype: Any
    ) -> None:
        mask = models.valid_mask(np.array(self.VALUES, dtype=dtype), SENTINEL)
        assert mask.dtype == np.bool_
        assert mask.tolist() == [True, False, True, True, False]

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    def test_without_a_sentinel_only_nan_is_nodata(self, models: ModuleType, dtype: Any) -> None:
        mask = models.valid_mask(np.array(self.VALUES, dtype=dtype), None)
        assert mask.tolist() == [True, False, True, True, True]

    def test_the_shape_is_kept(self, models: ModuleType) -> None:
        values = np.array(self.VALUES * 2, dtype=np.float32).reshape(2, 5)
        assert models.valid_mask(values, SENTINEL).shape == (2, 5)
