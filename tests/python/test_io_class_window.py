"""`ClassWindow`: `uint8` class blocks, empty tiles as zeros, `PixelIsArea` corners,
strips over a box (26a, `tin_engine.io.models`, `io.cog`, `io.repository`).

`docs/increments/26-land-cover.md`, section 3.1 ("Reading") and section 10's
"26a (Python)". The pixel oracle is `class_fixtures.expected`: the grid as
written, each empty tile's rectangle zeroed by hand from tifffile's page
attributes, sharing nothing with the reader. The pinned spellings are in
`class_fixtures`' docstring.

HOW THIS FILE GOES RED: `ClassWindow`, `decode_class_window`,
`load_class_window` and `class_strips` are reached through fixtures or
`getattr` at call time (`AttributeError`), so each test fails on its own and
the rest of `tests/python` collects.
"""

from __future__ import annotations

import importlib
import io
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from class_fixtures import (
    COLS,
    KEY,
    LAT0,
    LON0,
    ROWS,
    STEP,
    class_grid,
    class_tiff,
    expected,
    object_id,
)
from cog_fixtures import AbsentBlocks, BytesBlocks, same_bytes, sliced, write_cache
from geotiff_fixtures import needs_codecs
from tin_engine.io.models import Bounds, IndexWindow

OBJ = object_id(2025)
EMPTY = (4, 8)  # the middle tile and the last (bottom-right) one
#: Degrees; the coordinates are of magnitude 45, checked at 45.
DEGREES = 1e-9


@pytest.fixture(scope="module")
def models() -> ModuleType:
    return importlib.import_module("tin_engine.io.models")


@pytest.fixture(scope="module")
def cog() -> ModuleType:
    return importlib.import_module("tin_engine.io.cog")


@pytest.fixture(scope="module")
def repository() -> ModuleType:
    return importlib.import_module("tin_engine.io.repository")


def meta_of(data: bytes) -> Any:
    from tin_engine.io.geotiff import read_page

    return read_page(io.BytesIO(data), nodata=0.0)[0]


def window(row0: int, col0: int, rows: int, cols: int) -> IndexWindow:
    return IndexWindow(row0=row0, col0=col0, rows=rows, cols=cols)


def decode(cog: ModuleType, source: Any, data: bytes, w: IndexWindow, **kw: Any) -> Any:
    return cog.decode_class_window(source, meta_of(data), w, **kw)


class TestTheModel:
    def test_a_class_window_is_a_read_only_uint8_grid_with_its_cells_corner(
        self, models: ModuleType
    ) -> None:
        grid = class_grid()
        cw = models.ClassWindow(x_min=LON0, y_max=LAT0, delta_x=STEP, delta_y=STEP, array=grid)
        assert cw.array.dtype == np.uint8 and not cw.array.flags.writeable
        assert (cw.rows, cw.cols) == (ROWS, COLS)
        assert same_bytes(np.asarray(cw.array), grid)
        assert cw.corner(0, 0) == (LON0, LAT0)
        x, y = cw.corner(5, 7)
        assert x == pytest.approx(LON0 + 7 * STEP, abs=DEGREES)
        assert y == pytest.approx(LAT0 - 5 * STEP, abs=DEGREES)
        grid[0, 0] = 255  # the caller's buffer is not shared
        assert cw.array[0, 0] != 255

    @pytest.mark.parametrize(
        "array",
        [class_grid().astype(np.float32), class_grid()[0], class_grid().astype(np.int16)],
        ids=["float32", "one_dimension", "int16"],
    )
    def test_anything_but_a_2d_uint8_array_is_refused(
        self, models: ModuleType, array: np.ndarray
    ) -> None:
        with pytest.raises(ValueError, match="uint8"):
            models.ClassWindow(x_min=LON0, y_max=LAT0, delta_x=STEP, delta_y=STEP, array=array)


class TestDecodeClassWindow:
    @pytest.mark.parametrize(
        "bounds",
        [(0, 0, ROWS, COLS), (2, 3, 5, 6), (10, 12, 20, 20), (16, 32, 16, 16), (39, 47, 1, 1)],
        ids=["whole", "inside_one_tile", "across_four", "on_tile_edges", "last_cell"],
    )
    def test_a_window_is_the_grids_slice_in_uint8(
        self, cog: ModuleType, bounds: tuple[int, int, int, int]
    ) -> None:
        data = class_tiff()
        w = window(*bounds)
        cw = decode(cog, BytesBlocks(data), data, w)
        assert cw.array.dtype == np.uint8
        assert same_bytes(np.asarray(cw.array), sliced(expected(data), w))

    def test_empty_tiles_read_as_zeros(self, cog: ModuleType) -> None:
        data = class_tiff(empty=EMPTY)
        oracle = expected(data)
        assert not same_bytes(oracle, class_grid()), "the emptied tiles must change the oracle"
        cw = decode(cog, BytesBlocks(data), data, window(0, 0, ROWS, COLS))
        assert same_bytes(np.asarray(cw.array), oracle)
        assert not np.asarray(cw.array)[16:32, 16:32].any()
        assert np.asarray(cw.array)[:16, :16].all()

    def test_pixel_is_area_the_corner_is_the_cells_corner_not_its_centre(
        self, cog: ModuleType
    ) -> None:
        data = class_tiff()
        meta = meta_of(data)
        assert meta.pixel_is_area and meta.x_min == pytest.approx(LON0 + STEP / 2, abs=DEGREES)
        cw = decode(cog, BytesBlocks(data), data, window(5, 7, 10, 9))
        assert cw.x_min == pytest.approx(LON0 + 7 * STEP, abs=DEGREES)
        assert cw.y_max == pytest.approx(LAT0 - 5 * STEP, abs=DEGREES)
        assert (cw.delta_x, cw.delta_y) == (STEP, STEP)
        assert (cw.rows, cw.cols) == (10, 9)
        whole = decode(cog, BytesBlocks(data), data, window(0, 0, ROWS, COLS))
        assert (whole.x_min, whole.y_max) == pytest.approx((LON0, LAT0), abs=DEGREES)

    @needs_codecs
    def test_lzw_with_the_horizontal_predictor_decodes_through_imagecodecs(
        self, cog: ModuleType
    ) -> None:
        data = class_tiff(empty=(4,), compression="lzw", predictor=2)
        cw = decode(cog, BytesBlocks(data), data, window(0, 0, ROWS, COLS))
        assert same_bytes(np.asarray(cw.array), expected(data))

    def test_one_thread_and_eight_give_the_same_bytes(self, cog: ModuleType) -> None:
        data = class_tiff(empty=EMPTY)
        w = window(3, 5, 30, 40)
        one = decode(cog, BytesBlocks(data), data, w, threads=1)
        eight = decode(cog, BytesBlocks(data), data, w, threads=8)
        assert same_bytes(np.asarray(one.array), np.asarray(eight.array))
        assert same_bytes(np.asarray(one.array), sliced(expected(data), w))

    def test_a_missing_block_is_not_cached_before_any_block_is_read(self, cog: ModuleType) -> None:
        data = class_tiff()
        source = AbsentBlocks(data, absent=(4,))
        with pytest.raises(cog.NotCached) as caught:
            decode(cog, source, data, window(10, 12, 20, 20))  # tiles 0, 1, 4, 5
        assert (caught.value.missing, caught.value.needed) == (1, 4)


@pytest.fixture
def cache(tmp_path: Path) -> tuple[Path, bytes]:
    data = class_tiff(empty=EMPTY)
    write_cache(tmp_path / "cache", KEY, {OBJ: data}, crs="EPSG:4326")
    return tmp_path / "cache", data


def strips(repo: Any, box: Bounds, rows: int) -> list[Any]:
    return list(repo.class_strips(OBJ, box, rows=rows))


def stacked(windows: list[Any]) -> np.ndarray:
    return np.vstack([np.asarray(w.array) for w in windows])


class TestTheRepository:
    def test_load_class_window_is_the_decode_from_the_cache(
        self, cog: ModuleType, repository: ModuleType, cache: tuple[Path, bytes]
    ) -> None:
        root, data = cache
        w = window(10, 12, 20, 20)
        cw = repository.CacheRepository(root, KEY, nodata=0.0).load_class_window(OBJ, w)
        direct = decode(cog, BytesBlocks(data), data, w)
        assert same_bytes(np.asarray(cw.array), np.asarray(direct.array))
        assert (cw.x_min, cw.y_max) == (direct.x_min, direct.y_max)
        assert same_bytes(np.asarray(cw.array), sliced(expected(data), w))

    def test_strips_over_a_box_cover_its_cells_top_to_bottom_and_no_more(
        self, repository: ModuleType, cache: tuple[Path, bytes]
    ) -> None:
        """Cells 10..23 x 4..9 (the box's corners lie inside those cells), in
        strips of at most 5 rows: 5, 5 and 4."""
        root, data = cache
        repo = repository.CacheRepository(root, KEY, nodata=0.0)
        box = Bounds(
            x_min=LON0 + 4.5 * STEP,
            y_min=LAT0 - 23.5 * STEP,
            x_max=LON0 + 9.5 * STEP,
            y_max=LAT0 - 10.5 * STEP,
        )
        got = strips(repo, box, rows=5)
        assert [s.rows for s in got] == [5, 5, 4]
        assert all(s.cols == 6 for s in got)
        for strip in got:
            assert strip.x_min == pytest.approx(LON0 + 4 * STEP, abs=DEGREES)
        y = LAT0 - 10 * STEP
        for strip in got:
            assert strip.y_max == pytest.approx(y, abs=DEGREES)
            y -= strip.rows * STEP
        assert same_bytes(stacked(got), expected(data)[10:24, 4:10])

    def test_a_box_on_cell_edges_takes_the_cells_whose_interior_it_meets(
        self, repository: ModuleType, cache: tuple[Path, bytes]
    ) -> None:
        root, data = cache
        repo = repository.CacheRepository(root, KEY, nodata=0.0)
        box = Bounds(
            x_min=LON0 + 5 * STEP,
            y_min=LAT0 - 20 * STEP,
            x_max=LON0 + 10 * STEP,
            y_max=LAT0 - 10 * STEP,
        )
        (strip,) = strips(repo, box, rows=100)
        assert (strip.rows, strip.cols) == (10, 5)
        assert strip.corner(0, 0) == pytest.approx((LON0 + 5 * STEP, LAT0 - 10 * STEP), abs=DEGREES)
        assert same_bytes(np.asarray(strip.array), expected(data)[10:20, 5:10])

    def test_a_box_past_the_rasters_edge_is_clipped_to_it(
        self, repository: ModuleType, cache: tuple[Path, bytes]
    ) -> None:
        root, data = cache
        repo = repository.CacheRepository(root, KEY, nodata=0.0)
        box = Bounds(x_min=LON0 - 1.0, y_min=LAT0 - 1.0, x_max=LON0 + 2.5 * STEP, y_max=LAT0 + 1.0)
        got = strips(repo, box, rows=16)
        assert [s.rows for s in got] == [16, 16, 8]
        assert same_bytes(stacked(got), expected(data)[:, :3])

    def test_a_box_meeting_no_cell_is_a_value_error(
        self, repository: ModuleType, cache: tuple[Path, bytes]
    ) -> None:
        root, _ = cache
        repo = repository.CacheRepository(root, KEY, nodata=0.0)
        box = Bounds(x_min=LON0 - 2.0, y_min=LAT0 - 2.0, x_max=LON0 - 1.0, y_max=LAT0 - 1.0)
        with pytest.raises(ValueError) as caught:
            strips(repo, box, rows=16)
        assert not isinstance(caught.value, repository.CacheError)

    def test_a_missing_block_is_not_cached(self, repository: ModuleType, tmp_path: Path) -> None:
        data = class_tiff()
        root = tmp_path / "cache"
        write_cache(root, KEY, {OBJ: data}, skip={OBJ: (4,)}, crs="EPSG:4326")
        repo = repository.CacheRepository(root, KEY, nodata=0.0)
        with pytest.raises(repository.NotCached):
            repo.load_class_window(OBJ, window(10, 12, 20, 20))
        box = Bounds(
            x_min=LON0 + 20 * STEP,
            y_min=LAT0 - 30 * STEP,
            x_max=LON0 + 25 * STEP,
            y_max=LAT0 - 20 * STEP,
        )
        with pytest.raises(repository.NotCached):
            strips(repo, box, rows=5)
