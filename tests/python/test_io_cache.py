"""The tile cache's read side and the repositories' windowed loads (23a-1, W4, W5, W9).

`docs/increments/23-basin-scale.md`, "Windowed source reads (23a-1)" and "The
cache". A cache here is written by `cog_fixtures.write_cache` from a file's
own byte ranges, so the cached pixels have a file to be compared with.

HOW THIS FILE GOES RED. The new names (`CacheRepository`, `CachedBlocks`,
`CacheManifest`, `NotCached`, `load_window`, `check`) are reached through
module-scoped fixtures or attributes, so each test fails on its own.
"""

from __future__ import annotations

import importlib
import io
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pydantic
import pytest

from cog_fixtures import (
    VARIANTS,
    block_file,
    brute_force_meeting,
    build,
    manifest_for,
    page_of,
    same_bytes,
    sliced,
    whole_page,
    write_cache,
)
from geotiff_fixtures import (
    EPSG_WGS84,
    GEOGRAPHIC_TYPE,
    PROJECTED_CS_TYPE,
    TIE_X,
    TIE_Y,
    with_keys,
)
from tin_engine.io.models import GeoTiffError

SOURCE = "test-utm33"
TILED = VARIANTS[0]
STRIPPED = VARIANTS[2]


@pytest.fixture(scope="module")
def repo() -> ModuleType:
    return importlib.import_module("tin_engine.io.repository")


@pytest.fixture(scope="module")
def cog() -> ModuleType:
    return importlib.import_module("tin_engine.io.cog")


@pytest.fixture(scope="module")
def geotiff() -> ModuleType:
    return importlib.import_module("tin_engine.io.geotiff")


@pytest.fixture(scope="module")
def window() -> Any:
    return importlib.import_module("tin_engine.io.models").IndexWindow


@pytest.fixture
def data() -> bytes:
    return build(TILED)


@pytest.fixture
def cache(tmp_path: Path, data: bytes) -> Path:
    """The cache root, holding one source of one object, `tile`."""
    write_cache(tmp_path / "cache", SOURCE, {"tile": data})
    return tmp_path / "cache"


def cached_blocks(repo: ModuleType, geotiff: ModuleType, directory: Path) -> Any:
    header = (directory / "header.bin").read_bytes()
    _, _, page = geotiff.read_page(io.BytesIO(header), nodata=None)
    return repo.CachedBlocks(directory, page, f"{SOURCE}/tile")


class TestW4TheCachesReadSide:
    def test_every_block_present_means_none_missing(
        self, repo: ModuleType, geotiff: ModuleType, cache: Path, data: bytes
    ) -> None:
        blocks = cached_blocks(repo, geotiff, cache / SOURCE / "tile")
        everything = tuple(range(len(page_of(data).dataoffsets)))
        assert blocks.missing(everything) == ()

    @pytest.mark.parametrize("change", [b"", b"\x00"], ids=["one_byte_short", "one_byte_long"])
    def test_a_block_file_of_the_wrong_length_is_missing(
        self, repo: ModuleType, geotiff: ModuleType, cache: Path, data: bytes, change: bytes
    ) -> None:
        directory = cache / SOURCE / "tile"
        path = block_file(directory, page_of(data), 6)
        body = path.read_bytes()
        path.write_bytes(body[:-1] if not change else body + change)
        assert cached_blocks(repo, geotiff, directory).missing((5, 6, 7)) == (6,)

    def test_a_part_file_is_ignored(
        self, repo: ModuleType, geotiff: ModuleType, cache: Path, data: bytes
    ) -> None:
        directory = cache / SOURCE / "tile"
        path = block_file(directory, page_of(data), 6)
        whole = path.read_bytes()
        path.unlink()
        path.with_name(path.name + ".part").write_bytes(whole)  # right bytes, wrong name
        assert cached_blocks(repo, geotiff, directory).missing((6,)) == (6,)

    def test_block_reads_the_files_exact_bytes(
        self, repo: ModuleType, geotiff: ModuleType, cache: Path, data: bytes
    ) -> None:
        page = page_of(data)
        blocks = cached_blocks(repo, geotiff, cache / SOURCE / "tile")
        for index in (0, 6, 19):
            offset, count = page.dataoffsets[index], page.databytecounts[index]
            assert blocks.block(index) == data[offset : offset + count]

    def test_the_manifest_round_trips_through_json(self, repo: ModuleType, data: bytes) -> None:
        manifest = manifest_for(SOURCE, {"tile": data, "other": build(STRIPPED)})
        again = repo.CacheManifest.model_validate_json(manifest.model_dump_json())
        assert again == manifest
        with pytest.raises(pydantic.ValidationError, match="frozen"):
            manifest.source = "x"

    def test_a_header_whose_sha256_differs_is_a_cache_error(
        self, repo: ModuleType, cache: Path
    ) -> None:
        header = cache / SOURCE / "tile" / "header.bin"
        body = bytearray(header.read_bytes())
        body[-1] ^= 0xFF
        header.write_bytes(bytes(body))
        with pytest.raises(repo.CacheError, match=r"re-fetch with --refresh") as caught:
            repo.CacheRepository(cache, SOURCE).footprints()
        assert not isinstance(caught.value, repo.NotCached)

    def test_a_manifest_of_another_source_is_a_cache_error(
        self, repo: ModuleType, tmp_path: Path, data: bytes
    ) -> None:
        directory = write_cache(tmp_path / "cache", SOURCE, {"tile": data})
        other = manifest_for("glo30", {"tile": data})
        (directory / "manifest.json").write_text(other.model_dump_json())
        with pytest.raises(repo.CacheError) as caught:
            repo.CacheRepository(tmp_path / "cache", SOURCE).footprints()
        assert not isinstance(caught.value, repo.NotCached)
        assert "glo30" in str(caught.value) and SOURCE in str(caught.value)

    def test_no_manifest_is_not_cached(self, repo: ModuleType, tmp_path: Path) -> None:
        with pytest.raises(repo.NotCached, match="not in the cache") as caught:
            repo.CacheRepository(tmp_path / "cache", SOURCE).footprints()
        assert caught.value.needed == 0

    def test_construction_reads_no_file(self, repo: ModuleType, cache: Path) -> None:
        repo.CacheRepository(cache / "absent", SOURCE)  # no error: nothing is read
        repository = repo.CacheRepository(cache, SOURCE)
        (cache / SOURCE / "manifest.json").unlink()  # after construction
        with pytest.raises(repo.NotCached):
            repository.footprints()

    def test_footprints_are_the_objects_by_id_with_the_files_header(
        self, repo: ModuleType, geotiff: ModuleType, tmp_path: Path
    ) -> None:
        east = build(TILED, x_min=TIE_X + 700.0)
        west = build(TILED)
        write_cache(tmp_path / "cache", SOURCE, {"west": west, "east": east})
        prints = repo.CacheRepository(tmp_path / "cache", SOURCE).footprints()
        assert [p.name for p in prints] == ["east", "west"]
        for footprint, blob in zip(prints, (east, west), strict=True):
            meta, dtype = geotiff.read_header(io.BytesIO(blob), nodata=None)
            assert (footprint.meta, footprint.dtype) == (meta, dtype)


class TestLoadWindow:
    """`load_window` on both repositories gives the file's pixels, sliced."""

    def test_the_cache_repository(
        self, repo: ModuleType, cog: ModuleType, window: Any, cache: Path, data: bytes
    ) -> None:
        repository = repo.CacheRepository(cache, SOURCE)
        (footprint,) = repository.footprints()
        w = window(row0=10, col0=12, rows=30, cols=40)
        tile = repository.load_window("tile", w)
        assert same_bytes(np.asarray(tile.array), sliced(whole_page(data), w))
        assert tile.meta == cog.window_meta(footprint.meta, w)

    def test_the_tiff_repository(
        self, repo: ModuleType, cog: ModuleType, window: Any, tmp_path: Path, data: bytes
    ) -> None:
        path = tmp_path / "tile.tif"
        path.write_bytes(data)
        repository = repo.TiffDemRepository([path])
        (footprint,) = repository.footprints()
        w = window(row0=0, col0=50, rows=50, cols=20)
        tile = repository.load_window("tile.tif", w)
        assert same_bytes(np.asarray(tile.array), sliced(whole_page(data), w))
        assert tile.meta == cog.window_meta(footprint.meta, w)


def two_objects(tmp_path: Path, skip: dict[str, tuple[int, ...]]) -> tuple[Path, bytes, bytes]:
    """West and east halves of a 50 x 139 mosaic, sharing one node column."""
    west, east = build(TILED), build(STRIPPED, x_min=TIE_X + 69 * 10.0)
    write_cache(tmp_path / "cache", SOURCE, {"east": east, "west": west}, skip=skip)
    return tmp_path / "cache", west, east


class TestW5CheckThePlan:
    def test_check_counts_the_plans_missing_and_needed_blocks_once(
        self, repo: ModuleType, tmp_path: Path
    ) -> None:
        from tin_engine.mosaic import Bounds, plan_mosaic

        # Missing: west 6 and 19 (19 is outside the box), east strip 1 and 6.
        root, west, east = two_objects(tmp_path, {"west": (6, 19), "east": (1, 6)})
        repository = repo.CacheRepository(root, SOURCE)
        box = Bounds(
            x_min=TIE_X + 200.0, y_min=TIE_Y - 120.0, x_max=TIE_X + 900.0, y_max=TIE_Y - 20.0
        )
        plan = plan_mosaic(repository.footprints(), box, None)
        needed = {
            p.name: brute_force_meeting(page_of(west if p.name == "west" else east), p.source)
            for p in plan.tiles
        }
        missing = {"west": (6, 19), "east": (1, 6)}
        expected_missing = sum(len(set(needed[n]) & set(missing[n])) for n in needed)
        assert expected_missing >= 2, "the fixture must miss a needed block in each object"
        with pytest.raises(repo.NotCached) as caught:
            repository.check(plan)
        assert caught.value.missing == expected_missing
        assert caught.value.needed == sum(len(v) for v in needed.values())

    def test_check_passes_when_every_needed_block_is_present(
        self, repo: ModuleType, tmp_path: Path
    ) -> None:
        from tin_engine.mosaic import Bounds, plan_mosaic

        root, _, _ = two_objects(tmp_path, {"west": (19,)})  # bottom-right tile only
        repository = repo.CacheRepository(root, SOURCE)
        box = Bounds(x_min=TIE_X, y_min=TIE_Y - 100.0, x_max=TIE_X + 300.0, y_max=TIE_Y)
        assert repository.check(plan_mosaic(repository.footprints(), box, None)) is None

    def test_load_window_repeats_the_check_per_window(
        self, repo: ModuleType, window: Any, tmp_path: Path
    ) -> None:
        root, _, _ = two_objects(tmp_path, {"west": (6,)})
        repository = repo.CacheRepository(root, SOURCE)
        repository.footprints()
        with pytest.raises(repo.NotCached) as caught:
            repository.load_window("west", window(row0=10, col0=12, rows=10, cols=10))
        assert (caught.value.missing, caught.value.needed) == (1, 4)

    def test_the_tiff_repositorys_check_is_a_no_op(self, repo: ModuleType, tmp_path: Path) -> None:
        from tin_engine.mosaic import plan_mosaic

        path = tmp_path / "tile.tif"
        path.write_bytes(build(TILED))
        repository = repo.TiffDemRepository([path])
        assert repository.check(plan_mosaic(repository.footprints(), None, None)) is None


class TestW9GeographicUntil15c2:
    def test_a_cached_geographic_source_gets_the_readers_geokey_refusal(
        self, repo: ModuleType, tmp_path: Path
    ) -> None:
        keys = with_keys({PROJECTED_CS_TYPE: None, GEOGRAPHIC_TYPE: EPSG_WGS84})
        data = build(TILED, geokeys=keys)
        write_cache(tmp_path / "cache", "geo", {"tile": data}, crs="EPSG:4326")
        with pytest.raises(GeoTiffError, match=r"GeographicTypeGeoKey \(2048\) = 4326") as caught:
            repo.CacheRepository(tmp_path / "cache", "geo").footprints()
        assert not isinstance(caught.value, repo.CacheError)
