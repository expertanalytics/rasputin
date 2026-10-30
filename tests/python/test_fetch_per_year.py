"""Fetching a one-COG-per-year source of classes (26a, `tin_engine.fetch.run`).

`docs/increments/26-land-cover.md`, section 3.1: the year, the object id per
year (`brazil_coverage-col11_<year>`), empty tiles as no data, the fetch date
per object and the credit rendered with it. No network: the server is
`fetch_fixtures.RangeServer` on 127.0.0.1 serving `class_fixtures.class_tiff`
for each year, the cache is under `tmp_path`, and the oracles are the served
files' bytes and the server's request log. The pinned spellings are in
`class_fixtures`' docstring.

HOW THIS FILE GOES RED: `RemoteSource` refuses `kind="one-cog-per-year"`
(`ValidationError` building the entry), `FetchRequest` has no `year`, and
`CachedObject` has no `fetched`; each test fails on its own and the rest of
`tests/python` collects.
"""

from __future__ import annotations

import datetime
import importlib
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from class_fixtures import (
    COLS,
    KEY,
    PATH_TEMPLATE,
    ROWS,
    expected,
    object_id,
    per_year_entry,
    serve_years,
    served_path,
    whole_box,
)
from cog_fixtures import same_bytes
from crs_fixtures import crs_box
from fetch_fixtures import PROJECTED_CRS, RangeServer, block_files, projected
from test_fetch_run import present_blocks, requested_blocks
from test_fetch_run import whole_box as dem_box
from tin_engine.io.models import IndexWindow

OBJ = object_id(2025)
OLD = object_id(1985)
EMPTY = (4, 8)
EARLIER = datetime.date(2026, 1, 1)


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


def classes(api: SimpleNamespace, server: RangeServer) -> Any:
    return per_year_entry(api.sources, server.url(PATH_TEMPLATE))


async def fetch(api: SimpleNamespace, source: Any, root: Path, **request: Any) -> Any:
    client = api.http.RangeClient(delays=(0, 0, 0))
    writer = api.repository.CacheWriter(root, source.id)
    given = crs_box(source.crs, whole_box())
    asked = api.plan.FetchRequest(source=source.id, box=given, margin=0, **request)
    return await api.run.fetch(asked, source, client, writer)


def manifest(api: SimpleNamespace, root: Path) -> Any:
    return api.repository.CacheManifest.model_validate_json(
        (root / KEY / "manifest.json").read_bytes()
    )


def nothing_written(root: Path) -> bool:
    directory = root / KEY
    return not directory.exists() or all(
        not block_files(p) for p in directory.iterdir() if p.is_dir()
    )


class TestThePerYearObject:
    async def test_the_default_year_is_the_last_and_names_the_object(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        files = serve_years(server, [2025])
        source = classes(api, server)
        report = await fetch(api, source, tmp_path)
        assert report.objects == 1
        assert (tmp_path / KEY / OBJ).is_dir()
        listed = manifest(api, tmp_path).objects
        assert set(listed) == {OBJ}
        assert listed[OBJ].url == source.url_for(2025) == server.url(served_path(2025))
        assert present_blocks(tmp_path / KEY / OBJ, files[2025]) == set(range(9))

    async def test_each_year_is_its_own_object_in_one_manifest(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        files = serve_years(server, [2025, 1985])
        source = classes(api, server)
        await fetch(api, source, tmp_path, year=2025)
        await fetch(api, source, tmp_path, year=1985)
        assert set(manifest(api, tmp_path).objects) == {OBJ, OLD}
        assert present_blocks(tmp_path / KEY / OBJ, files[2025]) == set(range(9))
        assert present_blocks(tmp_path / KEY / OLD, files[1985]) == set(range(9))
        assert files[2025] != files[1985]

    async def test_a_year_outside_the_range_is_refused_before_any_request(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        serve_years(server, [2025])
        with pytest.raises(api.http.FetchError, match=r"1985.*2025") as caught:
            await fetch(api, classes(api, server), tmp_path, year=1984)
        assert "1984" in str(caught.value)
        assert server.log == []
        assert nothing_written(tmp_path)

    async def test_a_year_given_for_a_source_without_years_is_refused(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        data = projected()
        server.put("p.tif", data)
        dem = api.sources.RemoteSource(
            id=KEY, kind="one-cog", url=server.url("p.tif"), crs=PROJECTED_CRS, nodata=-9999.0,
            credit="c", licence_note="l",
        )  # fmt: skip
        client = api.http.RangeClient(delays=(0, 0, 0))
        writer = api.repository.CacheWriter(tmp_path, KEY)
        given = crs_box(PROJECTED_CRS, dem_box(data))
        asked = api.plan.FetchRequest(source=KEY, box=given, margin=0, year=2025)
        with pytest.raises(api.http.FetchError, match="year"):
            await api.run.fetch(asked, dem, client, writer)
        assert server.log == []
        assert nothing_written(tmp_path)


class TestEmptyTiles:
    async def test_empty_tiles_are_never_requested_and_are_empty_files(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        files = serve_years(server, [2025], empty=EMPTY)
        report = await fetch(api, classes(api, server), tmp_path)
        asked = requested_blocks(server, served_path(2025), files[2025])
        assert sorted(asked) == [0, 1, 2, 3, 5, 6, 7]
        blocks = tmp_path / KEY / OBJ / "blocks"
        for name in ("1/1.bin", "2/2.bin"):  # indices 4 and 8, three across
            assert (blocks / name).is_file() and (blocks / name).stat().st_size == 0, name
        assert (report.needed, report.empty, report.fetched) == (9, 2, 7)

    async def test_the_cached_class_window_reads_empty_tiles_as_zeros(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        files = serve_years(server, [2025], empty=EMPTY)
        await fetch(api, classes(api, server), tmp_path)
        repo = api.repository.CacheRepository(tmp_path, KEY, nodata=0.0)
        cw = repo.load_class_window(OBJ, IndexWindow(row0=0, col0=0, rows=ROWS, cols=COLS))
        assert same_bytes(np.asarray(cw.array), expected(files[2025]))


class TestTheFetchDate:
    async def test_each_object_is_dated_the_day_its_header_was_fetched(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        serve_years(server, [2025, 1985])
        source = classes(api, server)
        await fetch(api, source, tmp_path)
        assert manifest(api, tmp_path).objects[OBJ].fetched == datetime.date.today()
        await fetch(api, source, tmp_path, year=1985)
        listed = manifest(api, tmp_path).objects
        assert {o.fetched for o in listed.values()} == {datetime.date.today()}

    async def test_the_first_date_is_kept_on_a_later_run(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        files = serve_years(server, [2025])
        source = classes(api, server)
        await fetch(api, source, tmp_path)
        listed = manifest(api, tmp_path)
        aged = listed.model_copy(
            update={"objects": {OBJ: listed.objects[OBJ].model_copy(update={"fetched": EARLIER})}}
        )
        (tmp_path / KEY / "manifest.json").write_text(aged.model_dump_json())
        server.reset()
        await fetch(api, source, tmp_path)
        assert manifest(api, tmp_path).objects[OBJ].fetched == EARLIER
        assert server.block_ranges(served_path(2025)) == []
        assert present_blocks(tmp_path / KEY / OBJ, files[2025]) == set(range(9))

    async def test_the_notice_carries_each_objects_credit_with_its_date(
        self, api: SimpleNamespace, server: RangeServer, tmp_path: Path
    ) -> None:
        serve_years(server, [2025, 1985])
        source = classes(api, server)
        today = datetime.date.today()
        await fetch(api, source, tmp_path)
        notice = (tmp_path / KEY / "NOTICE.txt").read_text()
        assert source.credit_for(url=source.url_for(2025), date=today) in notice
        assert "{date}" not in notice and "{url}" not in notice
        assert source.licence_note in notice and source.cite[0] in notice
        await fetch(api, source, tmp_path, year=1985)
        notice = (tmp_path / KEY / "NOTICE.txt").read_text()
        assert source.credit_for(url=source.url_for(2025), date=today) in notice
        assert source.credit_for(url=source.url_for(1985), date=today) in notice
