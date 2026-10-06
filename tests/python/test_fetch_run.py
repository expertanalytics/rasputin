"""F3-F8 and F10: `fetch` against the local range server (23a-2, `tin_engine.fetch.run`).

`docs/increments/23-basin-scale.md`, "The cache", "Downloading" and 23a-2's
F3-F8 and F10. No network: the server is `fetch_fixtures.RangeServer` on
127.0.0.1, retry delays are zero, and the cache is written under `tmp_path`.
The oracles are the served file's own bytes and the server's request log;
nothing here reads the planner's records to decide what should have been
fetched, except where a test says it compares a dry run with a real one.

API pinned in `fetch_fixtures`' docstring. HOW THIS FILE GOES RED: the `api`
fixture imports `tin_engine.fetch.*` and `CacheWriter`, so every test errors
at setup until 23a-2 lands; the rest of `tests/python` collects.
"""

from __future__ import annotations

import datetime
import importlib
import importlib.metadata
import os
from collections.abc import Iterator, Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import pytest

from cog_fixtures import VARIANTS, build, with_zero_byte_count
from fetch_fixtures import (
    CHANGED,
    MIB,
    PROJECTED_CRS,
    RangeServer,
    block_files,
    block_span,
    glo30_name,
    glo30_tile,
    padded,
    page_of,
    projected,
    snapshot,
    wide,
    with_long_header,
)
from geotiff_fixtures import TIE_X, TIE_Y

if TYPE_CHECKING:
    from tin_engine.io.models import Bounds

SOURCE = "test-src"
CITES = ("First cited work, 2024.", "Second cited work, 2025.")


@pytest.fixture(scope="module")
def api() -> SimpleNamespace:
    return SimpleNamespace(
        plan=importlib.import_module("tin_engine.fetch.plan"),
        run=importlib.import_module("tin_engine.fetch.run"),
        http=importlib.import_module("tin_engine.fetch.http"),
        repository=importlib.import_module("tin_engine.io.repository"),
        sources=importlib.import_module("tin_engine.sources"),
    )


@pytest.fixture
def server() -> Iterator[RangeServer]:
    with RangeServer() as running:
        yield running


def entry(api: SimpleNamespace, server: RangeServer, path: str, **changes: Any) -> Any:
    fields: dict[str, Any] = {
        "id": SOURCE,
        "kind": "one-cog",
        "url": server.url(path),
        "crs": PROJECTED_CRS,
        "nodata": -9999.0,
        "credit": "Test credit line",
        "licence_note": "Test licence note",
        "cite": CITES,
    }
    return api.sources.RemoteSource(**(fields | changes))


def served(server: RangeServer, path: str, data: bytes) -> str:
    """`data` served at `path`; returns `path`."""
    server.put(path, data)
    return path


def bounds(**corners: float) -> Bounds:
    """A `Bounds`, from `io.models` (audit PR A moved it there from `mosaic`),
    read at call time so a missing name fails the test, not the collection."""
    box: Bounds = importlib.import_module("tin_engine.io.models").Bounds(**corners)
    return box


def whole_box(data: bytes) -> Bounds:
    """A box a little inside the served file's node rectangle."""
    page = page_of(data)
    rows, cols = int(page.imagelength), int(page.imagewidth)
    scale = page.tags[33550].value
    return bounds(
        x_min=TIE_X + 0.3 * scale[0],
        y_min=TIE_Y - (rows - 1.3) * scale[1],
        x_max=TIE_X + (cols - 1.3) * scale[0],
        y_max=TIE_Y - 0.3 * scale[1],
    )


def part_box(x0: float, x1: float, rows: float = 30.0) -> Bounds:
    """Columns `x0`..`x1` (in cells of 10 m) of the top `rows` cells (5 m)."""
    return bounds(
        x_min=TIE_X + x0 * 10, y_min=TIE_Y - rows * 5, x_max=TIE_X + x1 * 10, y_max=TIE_Y - 0.4
    )


async def fetch(api: SimpleNamespace, source: Any, root: Path, box: Bounds, **request: Any) -> Any:
    client = api.http.RangeClient(delays=(0, 0, 0))
    writer = api.repository.CacheWriter(root, source.id)
    asked = api.plan.FetchRequest(source=source.id, box=box, margin=0, **request)
    return await api.run.fetch(asked, source, client, writer)


def requested_blocks(server: RangeServer, path: str, data: bytes) -> list[int]:
    """Every block whose span lies inside a logged block range, once per range."""
    page = page_of(data)
    spans = [block_span(page, i) for i in range(len(page.dataoffsets))]
    return [
        i
        for start, stop in server.block_ranges(path)
        for i, (a, b) in enumerate(spans)
        if b > a and start <= a and b <= stop
    ]


def present_blocks(object_dir: Path, data: bytes) -> set[int]:
    page = page_of(data)
    across = -(-int(page.imagewidth) // int(page.chunks[1]))
    out = set()
    for name, body in block_files(object_dir).items():
        if name.endswith(".bin"):
            row, col = name.removesuffix(".bin").split("/")
            index = int(row) * across + int(col)
            assert body == data[slice(*block_span(page, index))], name
            out.add(index)
    return out


def manifest(api: SimpleNamespace, root: Path) -> Any:
    return api.repository.CacheManifest.model_validate_json(
        (root / SOURCE / "manifest.json").read_bytes()
    )


def object_dir(root: Path) -> Path:
    """The one object of a one-cog source."""
    (directory,) = [p for p in (root / SOURCE).iterdir() if p.is_dir()]
    return directory


# --------------------------------------------------------------------------
# F3: the bytes
# --------------------------------------------------------------------------


class TestF3TheBytes:
    async def test_every_block_is_the_files_byte_range(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = projected()
        source = entry(api, server, served(server, "p.tif", data))
        report = await fetch(api, source, tmp_path, whole_box(data))
        check_whole_fetch(server, data, tmp_path, report)

    async def test_a_second_run_requests_no_present_block(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = projected()
        source = entry(api, server, served(server, "p.tif", data))
        first = await fetch(api, source, tmp_path, whole_box(data))
        server.reset()
        again = await fetch(api, source, tmp_path, whole_box(data))
        assert server.block_ranges("p.tif") == []
        assert (again.fetched, again.present, again.requests) == (0, first.needed, 0)

    async def test_a_sparse_block_is_an_empty_file_never_requested(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = padded(with_zero_byte_count(build(VARIANTS[0]), 6))
        source = entry(api, server, served(server, "s.tif", data))
        report = await fetch(api, source, tmp_path, whole_box(data))
        check_sparse(server, data, tmp_path, report)


def check_whole_fetch(server: RangeServer, data: bytes, root: Path, report: Any) -> None:
    page = page_of(data)
    directory = object_dir(root)
    assert present_blocks(directory, data) == set(range(len(page.dataoffsets)))
    header = (directory / "header.bin").read_bytes()
    assert len(header) == MIB and header == data[:MIB]
    asked = requested_blocks(server, "p.tif", data)
    assert sorted(asked) == sorted(set(asked)), "a block was requested twice"
    assert len(server.block_ranges("p.tif")) == report.requests
    assert (report.needed, report.fetched, report.present) == (len(page.dataoffsets),) * 2 + (0,)


def check_sparse(server: RangeServer, data: bytes, root: Path, report: Any) -> None:
    sparse = object_dir(root) / "blocks" / "1" / "1.bin"  # index 6, five across
    assert sparse.is_file() and sparse.stat().st_size == 0
    assert 6 not in requested_blocks(server, "s.tif", data)
    assert report.empty == 1
    assert report.fetched == report.needed - 1


# --------------------------------------------------------------------------
# F4: resume
# --------------------------------------------------------------------------


class TestF4Resume:
    async def test_a_failed_run_resumes_to_a_clean_runs_tree(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        """Nodes 10..300 of the top 40 rows: about 20 blocks in each of the
        three block rows, 80 KiB apart in the file, so three ranges. With one
        connection and 500 after the prefix and one range, the run stops
        part way."""
        data = wide()
        source = entry(api, server, served(server, "w.tif", data))
        box = part_box(10.5, 300.5, rows=40)
        broken, clean = tmp_path / "broken", tmp_path / "clean"
        await fetch(api, source, clean, box, connections=1)
        wanted = present_blocks(object_dir(clean), data)
        assert len(server.block_ranges("w.tif")) == 3, "the box must give three ranges"
        server.reset()
        server.faults.fail_after = 2
        with pytest.raises(api.http.FetchError):
            await fetch(api, source, broken, box, connections=1)
        missing = missing_after_failure(broken, data, wanted)
        plant_part(broken)
        server.reset()
        await fetch(api, source, broken, box, connections=1)
        assert sorted(requested_blocks(server, "w.tif", data)) == missing
        assert snapshot(broken / SOURCE) == snapshot(clean / SOURCE)


def missing_after_failure(root: Path, data: bytes, wanted: set[int]) -> list[int]:
    """The blocks a clean run of the same box has that the failed run left
    out; the failure must have come part way, after some blocks and before all."""
    directory = object_dir(root)
    assert not [n for n in block_files(directory) if n.endswith(".part")]
    have = present_blocks(directory, data)
    assert have, "the failure came before any block: k too small"
    assert have < wanted, "the failure came after every block: k too large"
    return sorted(wanted - have)


def plant_part(root: Path) -> None:
    stray = object_dir(root) / "blocks" / "0" / "99.bin.part"
    stray.parent.mkdir(parents=True, exist_ok=True)
    stray.write_bytes(b"half a block")


# --------------------------------------------------------------------------
# F5: incremental
# --------------------------------------------------------------------------


class TestF5Incremental:
    async def test_a_second_domain_requests_only_its_new_blocks(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = wide()
        source = entry(api, server, served(server, "w.tif", data))
        await fetch(api, source, tmp_path, part_box(10.5, 400.5))
        before = present_blocks(object_dir(tmp_path), data)
        server.reset()
        await fetch(api, source, tmp_path, part_box(300.5, 700.5))
        after = present_blocks(object_dir(tmp_path), data)
        asked = set(requested_blocks(server, "w.tif", data))
        assert asked == after - before and not asked & before
        server.reset()
        await fetch(api, source, tmp_path, part_box(10.5, 400.5))
        check_requests_logged_once(api, tmp_path)


def check_requests_logged_once(api: SimpleNamespace, root: Path) -> None:
    requests = manifest(api, root).requests
    assert len(requests) == 2
    assert len({r.region_sha256 for r in requests}) == 2
    assert {r.date for r in requests} == {datetime.date.today()}


# --------------------------------------------------------------------------
# F6: refusals
# --------------------------------------------------------------------------


def nothing_written(root: Path) -> bool:
    directory = root / SOURCE
    return not directory.exists() or all(
        not block_files(p) for p in directory.iterdir() if p.is_dir()
    )


class TestF6Refusals:
    async def test_200_instead_of_206_is_refused_without_reading_the_body(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = padded(projected(), 16 * MIB)
        source = entry(api, server, served(server, "p.tif", data))
        server.faults.ignore_range = True
        with pytest.raises(api.http.FetchError, match=r"p\.tif"):
            await fetch(api, source, tmp_path, whole_box(data))
        server.wait_done()
        assert server.log and all(e.sent < len(data) for e in server.log)
        assert nothing_written(tmp_path)

    @pytest.mark.parametrize("fault", ["short", "wrong_content_range"])
    async def test_a_bad_body_or_content_range_is_refused(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path, fault: str
    ) -> None:
        data = projected()
        source = entry(api, server, served(server, "p.tif", data))
        setattr(server.faults, fault, True)
        with pytest.raises(api.http.FetchError, match=r"p\.tif"):
            await fetch(api, source, tmp_path, whole_box(data))
        assert nothing_written(tmp_path)

    @pytest.mark.parametrize("change", ["total", "last_modified"])
    async def test_a_block_response_from_a_changed_remote(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path, change: str
    ) -> None:
        data = projected()
        source = entry(api, server, served(server, "p.tif", data))
        if change == "total":
            server.faults.block_total = 1
        else:
            server.faults.block_last_modified = CHANGED
        with pytest.raises(api.http.FetchError, match="--refresh"):
            await fetch(api, source, tmp_path, whole_box(data))
        assert nothing_written(tmp_path)

    @pytest.mark.parametrize("change", ["length", "last_modified"])
    async def test_a_known_objects_prefix_from_a_changed_remote(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path, change: str
    ) -> None:
        data = wide()
        source = entry(api, server, served(server, "w.tif", data))
        await fetch(api, source, tmp_path, part_box(10.5, 100.5))
        before = snapshot(tmp_path)
        if change == "length":
            server.put("w.tif", data + bytes(16))
        else:
            server.put("w.tif", data, CHANGED)
        server.reset()
        with pytest.raises(api.http.FetchError, match="--refresh"):
            await fetch(api, source, tmp_path, part_box(500.5, 700.5))
        assert server.block_ranges("w.tif") == []
        assert snapshot(tmp_path) == before

    async def test_4xx_is_tried_once(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        source = entry(api, server, served(server, "p.tif", projected()))
        server.faults.status = 404
        with pytest.raises(api.http.FetchError, match=r"p\.tif"):
            await fetch(api, source, tmp_path, whole_box(projected()))
        assert len(server.log) == 1

    async def test_5xx_is_tried_four_times(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        source = entry(api, server, served(server, "p.tif", projected()))
        server.faults.status = 500
        with pytest.raises(api.http.FetchError, match=r"p\.tif"):
            await fetch(api, source, tmp_path, whole_box(projected()))
        assert len(server.log) == 4

    async def test_a_held_lock_is_refused(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = projected()
        source = entry(api, server, served(server, "p.tif", data))
        with (
            api.repository.CacheWriter(tmp_path, SOURCE),
            pytest.raises(api.repository.CacheError, match=SOURCE),
        ):
            await fetch(api, source, tmp_path, whole_box(data))
        assert nothing_written(tmp_path)

    async def test_a_header_crs_other_than_the_catalogues_is_refused(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = projected()
        source = entry(api, server, served(server, "p.tif", data), crs="EPSG:3035")
        with pytest.raises(api.http.FetchError, match="3035"):
            await fetch(api, source, tmp_path, whole_box(data))
        assert nothing_written(tmp_path)


# --------------------------------------------------------------------------
# F7: the header
# --------------------------------------------------------------------------


class TestF7TheHeader:
    async def test_a_long_header_doubles_the_prefix(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = with_long_header(3 * MIB)
        source = entry(api, server, served(server, "h.tif", data))
        await fetch(api, source, tmp_path, whole_box(data))
        below_first_block = min(page_of(data).dataoffsets)
        prefixes = [b for a, b in server.ranges("h.tif") if b <= below_first_block or a == 0]
        assert sorted(prefixes) == [MIB, 2 * MIB, 4 * MIB]
        assert (object_dir(tmp_path) / "header.bin").read_bytes() == data[: 4 * MIB]

    async def test_a_header_past_64_mib_is_refused(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = with_long_header(65 * MIB)
        source = entry(api, server, served(server, "h.tif", data))
        with pytest.raises(api.http.FetchError, match="64"):
            await fetch(api, source, tmp_path, whole_box(data))
        assert max(b for _, b in server.ranges("h.tif")) == 64 * MIB
        assert nothing_written(tmp_path)


# --------------------------------------------------------------------------
# F8: GLO-30
# --------------------------------------------------------------------------


class TestF8Glo30:
    async def test_an_unlisted_tile_is_no_tile_and_never_requested(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        listed = [glo30_name(59, 10), glo30_name(59, 11)]
        for name, (lat, lon) in zip(listed, [(59, 10), (59, 11)], strict=True):
            server.put(f"{name}/{name}.tif", glo30_tile(lat, lon))
        server.put("tileList.txt", ("\n".join([*listed, glo30_name(10, 10)]) + "\n").encode())
        source = api.sources.RemoteSource(
            id=SOURCE,
            kind="cog-tiles",
            url_template=server.url("{name}/{name}.tif"),
            tile_list_url=server.url("tileList.txt"),
            crs="EPSG:4326",
            nodata=None,
            credit="Test credit",
            licence_note="Test licence",
        )
        box = bounds(x_min=10.8, y_min=59.8, x_max=11.2, y_max=60.2)
        report = await fetch(api, source, tmp_path, box)
        unlisted = sorted([glo30_name(60, 10), glo30_name(60, 11)])
        assert sorted(report.no_tile) == unlisted
        assert report.objects == 2
        assert not [e for e in server.log if "N60" in e.path]
        assert sorted(manifest(api, tmp_path).objects) == sorted(listed)


# --------------------------------------------------------------------------
# F10: the write side
# --------------------------------------------------------------------------


class TestF10TheWriteSide:
    async def test_a_failed_manifest_write_leaves_the_previous_one(
        self,
        api: SimpleNamespace,
        server: RangeServer,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        data = wide()
        source = entry(api, server, served(server, "w.tif", data))
        await fetch(api, source, tmp_path, part_box(10.5, 100.5))
        previous = (tmp_path / SOURCE / "manifest.json").read_bytes()
        plant_failing_manifest_replace(monkeypatch)
        with pytest.raises(OSError, match="planted"):
            await fetch(api, source, tmp_path, part_box(500.5, 700.5))
        assert (tmp_path / SOURCE / "manifest.json").read_bytes() == previous
        assert len(manifest(api, tmp_path).requests) == 1

    async def test_notice_holds_credit_licence_and_every_citation(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = projected()
        source = entry(api, server, served(server, "p.tif", data))
        await fetch(api, source, tmp_path, whole_box(data))
        check_notice(api, source, tmp_path)

    async def test_without_install_metadata_the_fetch_completes_and_records_unknown(
        self,
        api: SimpleNamespace,
        server: RangeServer,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """A source checkout run without `pip install` has no rasputin metadata.

        `rasputin version` already prints "unknown (...)" then; the manifest's
        `rasputin_version` takes the same "unknown" placeholder instead of the
        fetch dying with PackageNotFoundError.
        """
        plant_missing_install_metadata(monkeypatch)
        data = projected()
        source = entry(api, server, served(server, "p.tif", data))
        report = await fetch(api, source, tmp_path, whole_box(data))
        check_whole_fetch(server, data, tmp_path, report)
        assert manifest(api, tmp_path).rasputin_version.startswith("unknown")

    async def test_refresh_empties_the_cache_and_refetches(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = wide()
        source = entry(api, server, served(server, "w.tif", data))
        await fetch(api, source, tmp_path, part_box(10.5, 100.5))
        first = present_blocks(object_dir(tmp_path), data)
        server.reset()
        await fetch(api, source, tmp_path, part_box(500.5, 600.5), refresh=True)
        second = present_blocks(object_dir(tmp_path), data)
        assert second and not first & second
        assert set(requested_blocks(server, "w.tif", data)) == second
        assert len(manifest(api, tmp_path).requests) == 1

    async def test_a_dry_run_writes_nothing_and_reports_the_real_runs_counts(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = wide()
        source = entry(api, server, served(server, "w.tif", data))
        root = tmp_path / "cache"
        root.mkdir()
        dry = await fetch(api, source, root, part_box(10.5, 300.5), dry_run=True)
        assert snapshot(tmp_path) == {"cache/": b""}
        await fetch(api, source, root, part_box(10.5, 100.5))
        before = snapshot(tmp_path)
        dry = await fetch(api, source, root, part_box(10.5, 300.5), dry_run=True)
        assert snapshot(tmp_path) == before
        real = await fetch(api, source, root, part_box(10.5, 300.5))
        assert counts(dry) == counts(real)
        assert dry.present > 0 and dry.fetched > 0


def counts(report: Any) -> tuple[Any, ...]:
    names = ("objects", "needed", "present", "fetched", "empty", "bytes", "requests", "no_tile")
    return tuple(getattr(report, n) for n in names)


def plant_failing_manifest_replace(monkeypatch: pytest.MonkeyPatch) -> None:
    """`os.replace` onto a `manifest.json` raises; every other replace runs."""
    original = os.replace

    def replace(src: Any, dst: Any, *args: Any, **kwargs: Any) -> None:
        if Path(dst).name == "manifest.json":
            raise OSError("planted: the manifest write failed")
        original(src, dst, *args, **kwargs)

    monkeypatch.setattr(os, "replace", replace)


def plant_missing_install_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    """`importlib.metadata.version("rasputin")` raises as it does uninstalled."""
    original = importlib.metadata.version

    def version(name: str) -> str:
        if name == "rasputin":
            raise importlib.metadata.PackageNotFoundError(name)
        return original(name)

    monkeypatch.setattr(importlib.metadata, "version", version)


def check_notice(api: SimpleNamespace, source: Any, root: Path) -> None:
    text = (root / SOURCE / "NOTICE.txt").read_text()
    assert text == api.sources.notice(source)
    for piece in (source.credit, source.licence_note, *CITES):
        assert piece in text, piece


class TestTheCatalogueNotes:
    def test_anadem_cites_laipelt_2024(self, api: SimpleNamespace) -> None:
        cite: Sequence[str] = api.sources.SOURCES["anadem-v1"].cite
        assert any("Laipelt" in c and "2024" in c for c in cite), cite

    def test_cite_defaults_to_none(self, api: SimpleNamespace) -> None:
        bare = api.sources.RemoteSource(
            id="x", kind="one-cog", url="u", crs="EPSG:4326", nodata=None, credit="c",
            licence_note="l",
        )  # fmt: skip
        assert bare.cite == ()

    def test_notice_of_every_catalogue_source(self, api: SimpleNamespace) -> None:
        for source in api.sources.SOURCES.values():
            text = api.sources.notice(source)
            for piece in (source.credit, source.licence_note, *source.cite):
                assert piece in text, (source.id, piece)
