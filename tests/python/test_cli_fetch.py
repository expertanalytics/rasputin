"""`rasputin fetch`, mesh after fetch, the boundary, and B16 (23a-2: F9, F10, F11).

`docs/increments/23-basin-scale.md`, "The CLI", "The I/O boundary", 23a-2's
F9-F11 and Ola's ruling B16 (a): a mesh made from a catalogue source carries
the source's `licence_note` and `cite`, and since increment 25
(`docs/increments/25-plain-output.md`, D2) its credit as `dem_credit` and its
key as `dem_source`, in the `.vtk` and the `.ply` alike. The server
is `fetch_fixtures.RangeServer` on 127.0.0.1, put into `SOURCES` by
`monkeypatch` under the key `test-fetch`; the CLI's retry delays are its own
(1, 2, 4 s), so no test here makes the CLI retry.

F9's oracle for the refusal is the test's own count: the mesh's windows are
planned by the mesh path's planner (`dem_input._domain_plan`, the object the
mesh evaluates) and each window's blocks are counted with 23a-1's
`blocks_meeting` against the block files on disk, not against the fetch's
plan.

HOW THIS FILE GOES RED: there is no `fetch` command yet (Typer answers 2,
"No such command"), no `tin_engine/fetch/http.py`, and the mesh file has no
`licence_note` or `cite` field (`KeyError`).
"""

from __future__ import annotations

import ast
import importlib
import socket
import sys
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from shapely import affinity
from shapely.geometry import Point, Polygon
from typer.testing import CliRunner

import tin_engine
from cog_fixtures import write_cache
from fetch_fixtures import PROJECTED_CRS, RangeServer, page_of, projected, snapshot
from geotiff_fixtures import TIE_X, TIE_Y
from plyread import read_ply
from test_cli_mesh import plain
from test_cli_mesh_domain import geojson
from test_cli_mesh_mosaic import field, same_mesh
from tin_engine.cli import app
from tin_engine.dem_input import _domain_plan
from tin_engine.domain import read_domain
from tin_engine.io.cog import blocks_meeting
from tin_engine.io.repository import CacheRepository
from vtkread import read_vtk

KEY = "test-fetch"
CREDIT = "Test credit line 23a-2"
LICENCE = "Test licence note: no liability"
CITES = ("First cited work, 2024.", "Second cited work, 2025.")
PACKAGE = Path(tin_engine.__file__).resolve().parent
USAGE, REFUSED = 2, 1
TRIANGLE = [(TIE_X + 200, TIE_Y - 140), (TIE_X + 280, TIE_Y - 140), (TIE_X + 240, TIE_Y - 100)]
runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})


def invoke(*args: str) -> tuple[int, str]:
    result = runner.invoke(app, list(args))
    return result.exit_code, plain(result.output)


@pytest.fixture(scope="module")
def sources() -> ModuleType:
    return importlib.import_module("tin_engine.sources")


@pytest.fixture
def server() -> Iterator[RangeServer]:
    with RangeServer() as running:
        yield running


@pytest.fixture
def data() -> bytes:
    return projected()


def catalogue_with(
    sources: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    url: str,
    crs: str = PROJECTED_CRS,
    cite: tuple[str, ...] = CITES,
) -> Any:
    """`SOURCES` plus the test entry, wherever a module bound the name."""
    entry = sources.RemoteSource(
        id=KEY, kind="one-cog", url=url, crs=crs, nodata=-9999.0, credit=CREDIT,
        licence_note=LICENCE, cite=cite,
    )  # fmt: skip
    original = sources.SOURCES
    patched = {**original, KEY: entry}
    for module in list(sys.modules.values()):
        name = getattr(module, "__name__", "") or ""
        if name.startswith("tin_engine") and getattr(module, "SOURCES", None) is original:
            monkeypatch.setattr(module, "SOURCES", patched)
    return entry


@pytest.fixture
def remote(
    sources: ModuleType, monkeypatch: pytest.MonkeyPatch, server: RangeServer, data: bytes
) -> Any:
    return catalogue_with(sources, monkeypatch, server.put("p.tif", data))


def no_network(*args: Any, **kwargs: Any) -> Any:
    raise AssertionError("rasputin mesh opened a socket")


# --------------------------------------------------------------------------
# F9: mesh after fetch (K7)
# --------------------------------------------------------------------------


def own_counts(cache: Path, domain_file: Path) -> tuple[int, int]:
    """(missing, needed) blocks of the mesh's windows against the files on disk."""
    repository = CacheRepository(cache, KEY)
    plan = _domain_plan(repository.footprints(), read_domain(domain_file, None))[0]
    missing = needed = 0
    for placement in plan.tiles:
        page = page_of((cache / KEY / placement.name / "header.bin").read_bytes())
        across = -(-int(page.imagewidth) // int(page.chunks[1]))
        for index in blocks_meeting(page, placement.source):
            path = cache / KEY / placement.name / "blocks" / str(index // across)
            needed += 1
            missing += not (path / f"{index % across}.bin").is_file()
    return missing, needed


class TestF9MeshAfterFetch:
    def test_the_triangles_box_corners_lie_outside_it(self) -> None:
        triangle = Polygon(TRIANGLE)
        x0, y0, x1, y1 = triangle.bounds
        corners = [(x0, y0), (x1, y1), (x0, y1)]
        assert not any(triangle.contains(Point(c)) for c in corners)

    def test_a_mesh_after_fetch_is_offline_and_equals_the_file_by_path(
        self,
        remote: Any,
        data: bytes,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        cache, domain = tmp_path / "cache", geojson(tmp_path / "d.geojson", TRIANGLE)
        code, output = invoke("fetch", KEY, "--domain", str(domain), "--cache", str(cache))
        assert code == 0, output
        tif = tmp_path / "dem.tif"
        tif.write_bytes(data)
        by_path, by_key = tmp_path / "path.vtk", tmp_path / "key.vtk"
        mesh = ("--domain", str(domain), "--tolerance", "1")
        code, output = invoke("mesh", "--dem", str(tif), *mesh, "--out", str(by_path))
        assert code == 0, output
        monkeypatch.setattr(socket, "socket", no_network)
        monkeypatch.setattr(socket, "create_connection", no_network)
        code, output = invoke(
            "mesh", "--dem", KEY, "--cache", str(cache), *mesh, "--out", str(by_key)
        )
        assert code == 0, output
        same_mesh(read_vtk(by_key.read_bytes()), read_vtk(by_path.read_bytes()))

    def test_a_domain_grown_past_the_fetch_is_refused_with_the_tests_counts(
        self, remote: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        cache, domain = tmp_path / "cache", geojson(tmp_path / "d.geojson", TRIANGLE)
        code, output = invoke("fetch", KEY, "--domain", str(domain), "--cache", str(cache))
        assert code == 0, output
        assert own_counts(cache, domain) == (0, own_counts(cache, domain)[1])
        monkeypatch.setattr(socket, "socket", no_network)
        for step in range(1, 13):
            grown = affinity.scale(Polygon(TRIANGLE), 1 + step / 4, 1 + step / 4)
            path = geojson(tmp_path / f"grown{step}.geojson", list(grown.exterior.coords)[:-1])
            missing, needed = own_counts(cache, path)
            if missing:
                break
        else:
            pytest.fail("no growth reached a block the fetch did not take")
        code, output = invoke(
            "mesh", "--dem", KEY, "--cache", str(cache), "--domain", str(path),
            "--tolerance", "1", "--out", str(tmp_path / "m.vtk"),
        )  # fmt: skip
        assert code == USAGE, output
        assert f"{missing:,} of the {needed:,} blocks" in output, output
        assert "".join(f"run: rasputin fetch {KEY} --domain {path}".split()) in "".join(
            output.split()
        )


# --------------------------------------------------------------------------
# F10 through the CLI: the dry run, the report, the exit codes
# --------------------------------------------------------------------------


class TestF10TheCommand:
    BBOX = (str(TIE_X + 101.3), str(TIE_Y - 187.7), str(TIE_X + 333.9), str(TIE_Y - 61.1))

    def test_a_dry_run_writes_nothing(
        self, remote: Any, server: RangeServer, tmp_path: Path
    ) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        code, output = invoke(
            "fetch", KEY, "--bbox", *self.BBOX, "--cache", str(cache), "--dry-run"
        )
        assert code == 0, output
        assert snapshot(tmp_path) == {"cache/": b""}
        assert server.block_ranges("p.tif") == []

    def test_a_run_writes_the_notice_and_names_the_object(
        self, remote: Any, sources: ModuleType, tmp_path: Path
    ) -> None:
        cache = tmp_path / "cache"
        code, output = invoke("fetch", KEY, "--bbox", *self.BBOX, "--cache", str(cache))
        assert code == 0, output
        assert (cache / KEY / "NOTICE.txt").read_text() == sources.notice(remote)

    def test_a_refusal_exits_1(
        self, sources: ModuleType, monkeypatch: pytest.MonkeyPatch, server: RangeServer,
        data: bytes, tmp_path: Path,
    ) -> None:  # fmt: skip
        catalogue_with(sources, monkeypatch, server.put("p.tif", data), crs="EPSG:3035")
        code, output = invoke("fetch", KEY, "--bbox", *self.BBOX, "--cache", str(tmp_path / "c"))
        assert code == REFUSED, output
        assert "3035" in output, output

    def test_domain_and_bbox_together_is_a_usage_error(self, remote: Any, tmp_path: Path) -> None:
        domain = geojson(tmp_path / "d.geojson", TRIANGLE)
        code, output = invoke(
            "fetch", KEY, "--bbox", *self.BBOX, "--domain", str(domain), "--cache", str(tmp_path)
        )
        assert code == USAGE, output
        assert "--bbox" in output and "--domain" in output, output

    def test_an_unknown_source_is_a_usage_error(self, tmp_path: Path) -> None:
        code, output = invoke(
            "fetch", "no-such-source", "--bbox", *self.BBOX, "--cache", str(tmp_path)
        )
        assert code == USAGE, output
        assert "no-such-source" in output and "No such command" not in output, output


# --------------------------------------------------------------------------
# F11: the boundary
# --------------------------------------------------------------------------


def network_imports(tree: ast.AST) -> set[str]:
    """`urllib*` and `http.client` imports anywhere in `tree`, lazy ones included."""
    found: set[str] = set()
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and not node.level:
            base = node.module or ""
            names = [base] + [f"{base}.{alias.name}" for alias in node.names]
        found |= {n for n in names if n.split(".")[0] == "urllib" or n == "http.client"}
    return found


class TestF11TheBoundary:
    def test_the_scan_finds_what_it_forbids(self) -> None:
        planted = "import urllib\nfrom http import client\ndef f():\n    import http.client\n"
        assert network_imports(ast.parse(planted)) == {"urllib", "http.client"}
        assert network_imports(ast.parse("from urllib.parse import urlsplit\n")) == {
            "urllib.parse",
            "urllib.parse.urlsplit",
        }

    def test_only_fetch_http_imports_urllib_or_http_client(self) -> None:
        importers = sorted(
            str(path.relative_to(PACKAGE))
            for path in PACKAGE.rglob("*.py")
            if network_imports(ast.parse(path.read_text()))
        )
        assert importers == ["fetch/http.py"]

    def test_no_cache_root_is_refused_with_b7s_message(
        self, remote: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        code, output = invoke("fetch", KEY, "--bbox", *TestF10TheCommand.BBOX)
        assert code == USAGE, output
        assert "no cache: set RASPUTIN_DATA" in output and "--cache" in output, output


# --------------------------------------------------------------------------
# B16 (a): the mesh file carries the licence note and the citations
# --------------------------------------------------------------------------


class TestB16TheMeshCarriesTheNotes:
    def test_licence_note_cite_and_credit_are_fields(
        self, sources: ModuleType, monkeypatch: pytest.MonkeyPatch, data: bytes, tmp_path: Path
    ) -> None:
        catalogue_with(sources, monkeypatch, "https://example.invalid/p.tif")
        write_cache(tmp_path / "cache", KEY, {"dem": data})
        out = tmp_path / "key.vtk"
        code, output = invoke(
            "mesh", "--dem", KEY, "--cache", str(tmp_path / "cache"), "--out", str(out)
        )
        assert code == 0, output
        vtk = read_vtk(out.read_bytes())
        assert "elevation_source" not in vtk.field_data
        assert field(vtk, "dem_credit") == CREDIT
        assert field(vtk, "dem_source") == KEY
        assert field(vtk, "licence_note") == LICENCE
        cite = field(vtk, "cite")
        assert all(c in cite for c in CITES), cite

    @pytest.mark.parametrize("cite", [CITES, ()], ids=["cited", "uncited"])
    def test_a_ply_carries_them_as_header_comments(
        self,
        sources: ModuleType,
        monkeypatch: pytest.MonkeyPatch,
        data: bytes,
        tmp_path: Path,
        cite: tuple[str, ...],
    ) -> None:
        """PLY has no field data: the notes are `comment` lines, as `crs` is;
        a `cite` line only when the source has citations (D2: the same
        fields as the `.vtk`)."""
        catalogue_with(sources, monkeypatch, "https://example.invalid/p.tif", cite=cite)
        write_cache(tmp_path / "cache", KEY, {"dem": data})
        out = tmp_path / "key.ply"
        code, output = invoke(
            "mesh", "--dem", KEY, "--cache", str(tmp_path / "cache"), "--out", str(out), "--ascii"
        )
        assert code == 0, output
        comments = read_ply(out.read_bytes())[0].comments
        assert not any(c.startswith("elevation ") for c in comments), comments
        assert f"dem_credit {CREDIT}" in comments, comments
        assert f"dem_source {KEY}" in comments, comments
        assert f"licence_note {LICENCE}" in comments, comments
        cited = [c for c in comments if c.startswith("cite ")]
        if cite:
            assert len(cited) == 1 and all(c in cited[0] for c in cite), comments
        else:
            assert cited == [], comments

    def test_a_mesh_by_path_carries_neither(self, data: bytes, tmp_path: Path) -> None:
        tif, out = tmp_path / "dem.tif", tmp_path / "path.vtk"
        tif.write_bytes(data)
        code, output = invoke("mesh", "--dem", str(tif), "--out", str(out))
        assert code == 0, output
        vtk = read_vtk(out.read_bytes())
        assert "licence_note" not in vtk.field_data and "cite" not in vtk.field_data
        assert "dem_credit" not in vtk.field_data
        assert field(vtk, "dem_source") == "dem.tif"
