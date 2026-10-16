"""`rasputin fetch --year` and the refusals by `values` through the CLI (26a, `cli.py`).

`docs/increments/26-land-cover.md`, section 3.1 and section 10's "26a
(Python)": `--year` on a per-year source, a year outside the range refused,
`--year` on a source without years refused, and `mesh --dem` refusing a source
of classes. The server is `fetch_fixtures.RangeServer` on 127.0.0.1, its entry
put into `SOURCES` by `monkeypatch` under `class_fixtures.KEY`; the `mesh`
refusal uses the real `mapbiomas-c11` key with the network closed, since the
refusal must come from the catalogue entry before any read.

HOW THIS FILE GOES RED: `fetch` has no `--year` ("No such option", exit 2,
asserted against), `RemoteSource` refuses the per-year kind, and `mesh --dem
mapbiomas-c11` does not name a class source.
"""

from __future__ import annotations

import importlib
import socket
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from class_fixtures import (
    KEY,
    PATH_TEMPLATE,
    REAL_KEY,
    bbox_args,
    object_id,
    per_year_entry,
    serve_years,
)
from cli_helpers import USAGE, invoke
from fetch_fixtures import RangeServer, projected
from test_cli_fetch import KEY as DEM_KEY
from test_cli_fetch import catalogue_with, install, no_network


@pytest.fixture(scope="module")
def sources() -> ModuleType:
    return importlib.import_module("tin_engine.sources")


@pytest.fixture
def server() -> Iterator[RangeServer]:
    with RangeServer() as running:
        yield running


@pytest.fixture
def classes(sources: ModuleType, monkeypatch: pytest.MonkeyPatch, server: RangeServer) -> Any:
    serve_years(server, [2025, 1985])
    return install(sources, monkeypatch, per_year_entry(sources, server.url(PATH_TEMPLATE)))


def fetched(cache: Path, *args: str) -> tuple[int, str]:
    return invoke("fetch", KEY, *args, "--bbox", *bbox_args(), "--cache", str(cache))


class TestFetchYear:
    def test_year_names_the_object_and_writes_its_cache(self, classes: Any, tmp_path: Path) -> None:
        code, output = fetched(tmp_path, "--year", "1985")
        assert code == 0, output
        assert object_id(1985) in output and object_id(2025) not in output, output
        assert (tmp_path / KEY / object_id(1985) / "header.bin").is_file()
        assert not (tmp_path / KEY / object_id(2025)).exists()

    def test_without_year_the_last_year_is_fetched(self, classes: Any, tmp_path: Path) -> None:
        code, output = fetched(tmp_path)
        assert code == 0, output
        assert object_id(2025) in output, output
        assert (tmp_path / KEY / object_id(2025) / "header.bin").is_file()

    def test_a_year_outside_the_range_is_a_usage_error_naming_the_range(
        self, classes: Any, server: RangeServer, tmp_path: Path
    ) -> None:
        code, output = fetched(tmp_path, "--year", "1984")
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "1984" in output and "1985" in output and "2025" in output, output
        assert server.log == []
        assert not (tmp_path / KEY).exists()

    def test_year_with_a_source_without_years_is_a_usage_error(
        self, sources: ModuleType, monkeypatch: pytest.MonkeyPatch, server: RangeServer,
        tmp_path: Path,
    ) -> None:  # fmt: skip
        catalogue_with(sources, monkeypatch, server.put("p.tif", projected()))
        code, output = invoke(
            "fetch", DEM_KEY, "--year", "2025", "--bbox", *bbox_args(), "--cache", str(tmp_path)
        )
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--year" in output and DEM_KEY in output, output
        assert server.log == []


class TestMeshRefusesClasses:
    def test_dem_with_a_class_source_is_refused_naming_classes(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(socket, "socket", no_network)
        monkeypatch.setattr(socket, "create_connection", no_network)
        out = tmp_path / "m.vtk"
        code, output = invoke(
            "mesh", "--dem", REAL_KEY, "--cache", str(tmp_path), "--out", str(out)
        )
        assert code == USAGE, output
        assert "Traceback" not in output, output
        assert REAL_KEY in output and "classes" in output and "--dem" in output, output
        assert "not in the cache" not in output, "refused for the wrong reason"
        assert not out.exists()
