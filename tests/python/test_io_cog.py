"""`tin_engine.io.cog`: decoding only the blocks a window needs (23a-1, W1-W3, W5).

`docs/increments/23-basin-scale.md`, "Windowed source reads (23a-1)". The
oracle for pixels is tifffile decoding the whole full-resolution page, then
slicing (`cog_fixtures.whole_page`); the oracle for which blocks a window
needs is a brute-force test of every block's pixel rectangle, computed from
tifffile's page attributes (`cog_fixtures.brute_force_meeting`).

HOW THIS FILE GOES RED. `tin_engine.io.cog`, `read_page` and the moved
`IndexWindow` are reached through module-scoped fixtures, so while they are
missing each test fails on its own (`ModuleNotFoundError`, `AttributeError`)
and the rest of `tests/python` still collects.
"""

from __future__ import annotations

import importlib
import io
import re
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import tifffile

from cog_fixtures import (
    COLS,
    ROWS,
    VARIANTS,
    AbsentBlocks,
    BytesBlocks,
    Variant,
    brute_force_meeting,
    build,
    page_of,
    same_bytes,
    sliced,
    variant_params,
    whole_page,
    with_zero_byte_count,
    write_object,
)
from geotiff_fixtures import needs_codecs
from tin_engine.io.models import GeoTiffError


@pytest.fixture(scope="module")
def cog() -> ModuleType:
    return importlib.import_module("tin_engine.io.cog")


@pytest.fixture(scope="module")
def geotiff() -> ModuleType:
    return importlib.import_module("tin_engine.io.geotiff")


@pytest.fixture(scope="module")
def window() -> Any:
    """`IndexWindow`, from its new home (23a-1 decided 4)."""
    return importlib.import_module("tin_engine.io.models").IndexWindow


@pytest.fixture(scope="module")
def repository() -> ModuleType:
    return importlib.import_module("tin_engine.io.repository")


#: (row0, col0, rows, cols) on the 50 x 70 grid of 16 x 16 tiles / 8-row strips.
WINDOWS = {
    "inside_one_block": (2, 3, 5, 6),
    "across_four_blocks": (10, 12, 10, 10),
    "on_block_edges": (16, 32, 16, 16),
    "reaching_last_row_and_col": (40, 60, 10, 10),
    "one_row": (17, 5, 1, 40),
    "one_col": (3, 47, 30, 1),
    "whole_raster": (0, 0, ROWS, COLS),
}


def decoded(
    cog: ModuleType, geotiff: ModuleType, source: Any, data: bytes, w: Any, **kw: Any
) -> Any:
    meta, dtype, _ = geotiff.read_page(io.BytesIO(data), nodata=None)
    return cog.decode_window(source, meta, dtype, w, **kw)


def local(cog: ModuleType, geotiff: ModuleType, data: bytes) -> Any:
    _, _, page = geotiff.read_page(io.BytesIO(data), nodata=None)
    return cog.LocalTiffBlocks(page, io.BytesIO(data), "fixture.tif")


class TestReadPage:
    """`read_page`: `read_header` plus the page, through the same `_header`."""

    @pytest.mark.parametrize("variant", variant_params())
    def test_meta_and_dtype_are_read_headers(self, geotiff: ModuleType, variant: Variant) -> None:
        data = build(variant)
        meta, dtype, page = geotiff.read_page(io.BytesIO(data), nodata=None)
        assert (meta, dtype) == geotiff.read_header(io.BytesIO(data), nodata=None)
        assert page.shape == (ROWS, COLS)
        assert page.dataoffsets == page_of(data).dataoffsets

    def test_a_header_prefix_is_enough(self, geotiff: ModuleType) -> None:
        """The cache's `header.bin` ends before the first block and the overview's IFD."""
        data = build(VARIANTS[0])
        prefix = data[: min(page_of(data).dataoffsets)]
        assert geotiff.read_page(io.BytesIO(prefix), nodata=None)[:2] == geotiff.read_header(
            io.BytesIO(data), nodata=None
        )

    def test_refusals_are_read_headers(self, geotiff: ModuleType) -> None:
        with pytest.raises(GeoTiffError, match=r"TIFF structure"):
            geotiff.read_page(io.BytesIO(b"not a tiff"), nodata=None)


class TestW1WindowEqualsWholeSliced:
    """W1: every window, every variant, from local and cached blocks alike."""

    @pytest.mark.parametrize("name", WINDOWS)
    @pytest.mark.parametrize("variant", variant_params())
    def test_local_blocks(
        self, cog: ModuleType, geotiff: ModuleType, window: Any, variant: Variant, name: str
    ) -> None:
        data = build(variant)
        w = window(**dict(zip(("row0", "col0", "rows", "cols"), WINDOWS[name], strict=True)))
        tile = decoded(cog, geotiff, local(cog, geotiff, data), data, w)
        assert same_bytes(np.asarray(tile.array), sliced(whole_page(data), w))
        meta = geotiff.read_header(io.BytesIO(data), nodata=None)[0]
        assert tile.meta == meta.windowed(w)
        assert not tile.array.flags.writeable

    @pytest.mark.parametrize("name", WINDOWS)
    @pytest.mark.parametrize("variant", variant_params())
    def test_cached_blocks_are_identical(
        self,
        cog: ModuleType,
        geotiff: ModuleType,
        repository: ModuleType,
        window: Any,
        tmp_path: Path,
        variant: Variant,
        name: str,
    ) -> None:
        data = build(variant)
        w = window(**dict(zip(("row0", "col0", "rows", "cols"), WINDOWS[name], strict=True)))
        directory = write_object(tmp_path / "obj", data)
        header = (directory / "header.bin").read_bytes()
        meta, dtype, page = geotiff.read_page(io.BytesIO(header), nodata=None)
        cached = cog.decode_window(
            repository.CachedBlocks(directory, page, "src/obj"), meta, dtype, w
        )
        from_file = decoded(cog, geotiff, local(cog, geotiff, data), data, w)
        assert cached.meta == from_file.meta
        assert same_bytes(np.asarray(cached.array), np.asarray(from_file.array))
        assert same_bytes(np.asarray(cached.array), sliced(whole_page(data), w))

    def test_int16_is_promoted_to_float32(
        self, cog: ModuleType, geotiff: ModuleType, window: Any
    ) -> None:
        data = build(next(v for v in VARIANTS if v.name == "tiled_int16"))
        tile = decoded(
            cog, geotiff, local(cog, geotiff, data), data, window(row0=0, col0=0, rows=4, cols=4)
        )
        assert tile.array.dtype == np.float32

    def test_windowed_moves_the_corner_by_whole_cells(self, window: Any) -> None:
        """`RasterMeta.windowed`, which replaced `cog.window_meta` (audit PR A)."""
        from mosaic_fixtures import meta as make_meta

        whole = make_meta(rows=ROWS, cols=COLS, nodata=-9999.0)
        w = window(row0=7, col0=11, rows=5, cols=3)
        expected = whole.model_copy(
            update={
                "x_min": whole.x_min + 11 * whole.delta_x,
                "y_max": whole.y_max - 7 * whole.delta_y,
                "rows": 5,
                "cols": 3,
            }
        )
        assert whole.windowed(w) == expected

    @pytest.mark.parametrize(
        "bounds",
        [(45, 0, 10, 10), (0, 65, 5, 10), (-1, 0, 2, 2), (0, -1, 2, 2), (0, 0, 0, 5), (0, 0, 5, 0)],
        ids=[
            "past_last_row",
            "past_last_col",
            "negative_row0",
            "negative_col0",
            "no_rows",
            "no_cols",
        ],
    )
    def test_a_window_outside_or_empty_is_a_value_error(
        self, cog: ModuleType, geotiff: ModuleType, window: Any, bounds: tuple[int, ...]
    ) -> None:
        data = build(VARIANTS[0])
        source = BytesBlocks(data)
        w = window(**dict(zip(("row0", "col0", "rows", "cols"), bounds, strict=True)))
        with pytest.raises(ValueError) as caught:
            decoded(cog, geotiff, source, data, w)
        assert not isinstance(caught.value, GeoTiffError), "a caller's bug, not a file's"
        assert source.calls == []


DTM10_FILE = next(
    (
        p / "rasputin_data" / "DTM10_UTM33_20260925" / "6400_1_10m_z33.tif"
        for p in Path(__file__).resolve().parents
        if (p / "rasputin_data" / "DTM10_UTM33_20260925" / "6400_1_10m_z33.tif").is_file()
    ),
    None,
)


@needs_codecs
@pytest.mark.skipif(DTM10_FILE is None, reason="needs ../rasputin_data/DTM10_UTM33_20260925")
def test_w1_a_dtm10_window_across_a_512_tile_corner(
    cog: ModuleType, geotiff: ModuleType, window: Any
) -> None:
    assert DTM10_FILE is not None
    w = window(row0=500, col0=1010, rows=30, cols=20)  # meets tiles (0,1) (0,2) (1,1) (1,2)
    with DTM10_FILE.open("rb") as stream:
        meta, dtype, page = geotiff.read_page(stream, nodata=None)
        assert page.chunks == (512, 512)
        tile = cog.decode_window(cog.LocalTiffBlocks(page, stream, DTM10_FILE.name), meta, dtype, w)
    with tifffile.TiffFile(DTM10_FILE) as tif:
        oracle = tif.pages.first.asarray()
    assert same_bytes(np.asarray(tile.array), sliced(oracle, w).astype(np.float32))
    assert tile.meta == meta.windowed(w)


class TestW2Determinism:
    @pytest.mark.parametrize("variant", variant_params())
    def test_one_thread_and_eight_give_the_same_bytes(
        self, cog: ModuleType, geotiff: ModuleType, window: Any, variant: Variant
    ) -> None:
        data = build(variant)
        w = window(row0=5, col0=9, rows=40, cols=55)
        one = decoded(cog, geotiff, BytesBlocks(data), data, w, threads=1)
        eight = decoded(cog, geotiff, BytesBlocks(data), data, w, threads=8)
        assert np.asarray(one.array).tobytes() == np.asarray(eight.array).tobytes()
        assert same_bytes(np.asarray(one.array), sliced(whole_page(data), w))


class TestW3OnlyTheNeededBlocks:
    @pytest.mark.parametrize("name", WINDOWS)
    @pytest.mark.parametrize("variant", variant_params())
    def test_blocks_meeting_is_the_brute_force_set(
        self, cog: ModuleType, window: Any, variant: Variant, name: str
    ) -> None:
        page = page_of(build(variant))
        w = window(**dict(zip(("row0", "col0", "rows", "cols"), WINDOWS[name], strict=True)))
        assert cog.blocks_meeting(page, w) == brute_force_meeting(page, w)

    @pytest.mark.parametrize("name", WINDOWS)
    @pytest.mark.parametrize("variant", variant_params())
    def test_each_needed_block_is_read_once_and_no_other(
        self, cog: ModuleType, geotiff: ModuleType, window: Any, variant: Variant, name: str
    ) -> None:
        data = build(variant)
        w = window(**dict(zip(("row0", "col0", "rows", "cols"), WINDOWS[name], strict=True)))
        source = BytesBlocks(data)
        decoded(cog, geotiff, source, data, w)
        assert sorted(source.calls) == list(brute_force_meeting(source.page, w))


class TestW5RefusalsBeforeDecoding:
    def test_a_missing_block_is_not_cached_before_any_block_is_read(
        self, cog: ModuleType, geotiff: ModuleType, window: Any
    ) -> None:
        data = build(VARIANTS[0])
        w = window(row0=10, col0=12, rows=10, cols=10)  # blocks 0, 1, 5, 6
        source = AbsentBlocks(data, absent=(1, 6, 19))
        with pytest.raises(cog.NotCached) as caught:
            decoded(cog, geotiff, source, data, w)
        assert (caught.value.missing, caught.value.needed) == (2, 4)
        assert caught.value.where == "cache/obj"
        assert isinstance(caught.value, cog.CacheError)
        assert isinstance(caught.value, ValueError)

    @pytest.mark.parametrize("variant", [VARIANTS[0], VARIANTS[2]], ids=lambda v: v.name)
    def test_a_sparse_block_in_the_window_is_refused_naming_it(
        self, cog: ModuleType, geotiff: ModuleType, window: Any, variant: Variant
    ) -> None:
        data = with_zero_byte_count(build(variant), 1)
        w = window(row0=10, col0=12, rows=10, cols=10)  # meets block 1 in both layouts
        assert 1 in brute_force_meeting(page_of(data), w)
        with pytest.raises(GeoTiffError, match=r"block 1\b") as caught:
            decoded(cog, geotiff, local(cog, geotiff, data), data, w)
        # "sparse", not just the block: without the check the codec fails on
        # the empty bytes and still names the block (review round 1).
        assert "sparse" in str(caught.value) and "fixture.tif" in str(caught.value)

    def test_a_sparse_block_outside_the_window_is_never_looked_at(
        self, cog: ModuleType, geotiff: ModuleType, window: Any
    ) -> None:
        data = with_zero_byte_count(build(VARIANTS[0]), 19)  # the last tile
        w = window(row0=0, col0=0, rows=20, cols=20)
        tile = decoded(cog, geotiff, local(cog, geotiff, data), data, w)
        assert same_bytes(np.asarray(tile.array), sliced(whole_page(data), w))

    def test_a_corrupt_block_is_a_geotiff_error_naming_where_and_the_index(
        self, cog: ModuleType, geotiff: ModuleType, window: Any
    ) -> None:
        data = build(VARIANTS[0])
        page = page_of(data)
        buffer = bytearray(data)
        offset, count = page.dataoffsets[6], page.databytecounts[6]
        buffer[offset : offset + count] = b"\xff" * count
        w = window(row0=10, col0=12, rows=10, cols=10)
        with pytest.raises(GeoTiffError) as caught:
            decoded(cog, geotiff, BytesBlocks(bytes(buffer), where="bad.tif"), bytes(buffer), w)
        assert "bad.tif" in str(caught.value)
        assert re.search(r"block 6\b", str(caught.value)), caught.value
