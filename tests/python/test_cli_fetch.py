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
import re
import shlex
import socket
import sys
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from shapely import affinity
from shapely.geometry import Point, Polygon

import tin_engine
from cli_helpers import geojson, invoke, no_cache_anywhere, runner, says_no_cache
from cog_fixtures import write_cache
from fetch_fixtures import PROJECTED_CRS, RangeServer, block_span, page_of, projected, snapshot
from geotiff_fixtures import TIE_X, TIE_Y
from plyread import read_ply
from test_cli_mesh_mosaic import field, same_mesh
from tin_engine.cli import app
from tin_engine.dem_input import _domain_plan
from tin_engine.io.cog import blocks_meeting
from tin_engine.io.domain_file import read_domain
from tin_engine.io.repository import CacheRepository
from vtkread import read_vtk

KEY = "test-fetch"
CREDIT = "Test credit line 23a-2"
LICENCE = "Test licence note: no liability"
CITES = ("First cited work, 2024.", "Second cited work, 2025.")
PACKAGE = Path(tin_engine.__file__).resolve().parent
USAGE, REFUSED = 2, 1
TRIANGLE = [(TIE_X + 200, TIE_Y - 140), (TIE_X + 280, TIE_Y - 140), (TIE_X + 240, TIE_Y - 100)]


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
    return install(sources, monkeypatch, entry)


def install(sources: ModuleType, monkeypatch: pytest.MonkeyPatch, entry: Any) -> Any:
    """`SOURCES` plus `entry` under its id, wherever a module bound the name."""
    original = sources.SOURCES
    patched = {**original, entry.id: entry}
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
    BOX_ARGS = ("--bbox", *BBOX, "--bbox-crs", PROJECTED_CRS)

    def test_a_dry_run_writes_nothing(
        self, remote: Any, server: RangeServer, tmp_path: Path
    ) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        code, output = invoke("fetch", KEY, *self.BOX_ARGS, "--cache", str(cache), "--dry-run")
        assert code == 0, output
        assert snapshot(tmp_path) == {"cache/": b""}
        assert server.block_ranges("p.tif") == []

    def test_a_run_writes_the_notice_and_names_the_object(
        self, remote: Any, sources: ModuleType, tmp_path: Path
    ) -> None:
        cache = tmp_path / "cache"
        code, output = invoke("fetch", KEY, *self.BOX_ARGS, "--cache", str(cache))
        assert code == 0, output
        assert (cache / KEY / "NOTICE.txt").read_text() == sources.notice(remote)

    def test_without_cache_or_data_root_it_fetches_into_the_xdg_cache_and_says_where(
        self, remote: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The default cache: `$XDG_CACHE_HOME/rasputin`, the source's directory
        under it, and the summary line ends `into <that directory>`. `HOME`
        points elsewhere, so a write under it would be seen."""
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        xdg, home = tmp_path / "xdg", tmp_path / "home"
        monkeypatch.setenv("XDG_CACHE_HOME", str(xdg))
        monkeypatch.setenv("HOME", str(home))
        result = runner.invoke(app, ["fetch", KEY, *self.BOX_ARGS])
        assert result.exit_code == 0, result.output
        where = xdg / "rasputin" / KEY
        assert (where / "NOTICE.txt").is_file(), sorted(map(str, tmp_path.rglob("*")))
        assert not home.exists()
        (summary,) = [line for line in result.output.splitlines() if " objects, " in line]
        assert summary.startswith(f"{KEY}: "), summary
        assert summary.endswith(f" into {where}"), summary

    def test_a_refusal_exits_1(
        self, sources: ModuleType, monkeypatch: pytest.MonkeyPatch, server: RangeServer,
        data: bytes, tmp_path: Path,
    ) -> None:  # fmt: skip
        catalogue_with(sources, monkeypatch, server.put("p.tif", data), crs="EPSG:3035")
        code, output = invoke(
            "fetch", KEY, "--bbox", *self.BBOX, "--bbox-crs", "EPSG:3035",
            "--cache", str(tmp_path / "c"),
        )  # fmt: skip
        assert code == REFUSED, output
        assert "3035" in output, output

    def test_domain_and_bbox_together_is_a_usage_error(self, remote: Any, tmp_path: Path) -> None:
        domain = geojson(tmp_path / "d.geojson", TRIANGLE)
        code, output = invoke(
            "fetch", KEY, *self.BOX_ARGS, "--domain", str(domain), "--cache", str(tmp_path)
        )
        assert code == USAGE, output
        assert "--bbox" in output and "--domain" in output, output

    def test_an_unknown_source_is_a_usage_error(self, tmp_path: Path) -> None:
        code, output = invoke("fetch", "no-such-source", *self.BOX_ARGS, "--cache", str(tmp_path))
        assert code == USAGE, output
        assert "no-such-source" in output and "No such command" not in output, output


# --------------------------------------------------------------------------
# `docs/increments/114-bbox-crs.md`, tests 6 and 8: the box's CRS
# --------------------------------------------------------------------------

#: `TestF10TheCommand.BBOX` moved into EPSG:4326 (`114-probes/box_probe.txt`).
#: Moved back, its edges move by at most 13.5 mm; the grown fetch box's edges
#: lie at least 1.3 m from a node line, so the same blocks are fetched.
BBOX_DEGREES = ("15.001791152", "59.536663469", "15.005904106", "59.537800446")


def dry_run_lines(cache: Path, *box: str) -> list[str]:
    """The plan lines of a dry run, without the summary's `into <cache>` tail."""
    result = runner.invoke(app, ["fetch", KEY, *box, "--cache", str(cache), "--dry-run"])
    assert result.exit_code == 0, result.output
    return [line.partition(" into ")[0] for line in result.stdout.splitlines()]


def untimed(lines: list[str]) -> list[str]:
    """`lines` without the summary's trailing elapsed time, `, <n>.<n> s`."""
    return [re.sub(r", \d+\.\d s$", "", line) for line in lines]


def refusal(*args: str) -> str:
    """The message of a refused run, unwrapped by Rich."""
    result = runner.invoke(app, list(args), standalone_mode=False)
    message = getattr(result.exception, "message", None)
    assert isinstance(message, str), (result.exception, result.output)
    return message


def hinted(message: str) -> list[str]:
    """The `rasputin fetch` the refusal names, split as a shell would, without `rasputin`."""
    assert "; run: rasputin fetch " in message, message
    words = shlex.split(message.rpartition("; run: ")[2])
    assert words[:2] == ["rasputin", "fetch"], words
    return words[1:]


class TestBoxCrs:
    def test_a_dry_run_in_the_sources_crs_plans_as_before(
        self, remote: Any, data: bytes, tmp_path: Path
    ) -> None:
        """Test 6: measured before `--bbox-crs`, with the box read in the source's CRS.
        The summary's elapsed time (about 0.06 s a run) is cut: it is the host's.
        The plan is blocks 0 to 12 (rows 0 to 2, columns 0 to 2 of 5 across) in
        one request; its byte count is that span of the fixture, whose Deflate
        sizes depend on the host's Deflate build (6,274 bytes on macOS, 6,312 on CI's Linux)."""
        lines = untimed(dry_run_lines(tmp_path / "c", *TestF10TheCommand.BOX_ARGS))
        page = page_of(data)
        size = f"{block_span(page, 12)[1] - block_span(page, 0)[0]:,}"
        assert lines == [
            f"p: 9 blocks, {size} bytes to fetch",
            "test-fetch: 1 objects, 9 blocks needed, 0 present, would fetch 9 (0 empty)"
            f" in 1 requests, {size} bytes",
        ], lines

    def test_a_dry_run_in_degrees_plans_the_same(self, remote: Any, tmp_path: Path) -> None:
        """Elapsed time is cut from each line: it is the host's (0.1 s and 0.2 s on CI)."""
        utm = dry_run_lines(tmp_path / "c", *TestF10TheCommand.BOX_ARGS)
        degrees = dry_run_lines(tmp_path / "c", "--bbox", *BBOX_DEGREES, "--bbox-crs", "EPSG:4326")
        assert untimed(degrees) == untimed(utm)

    def test_the_mesh_hint_names_the_boxs_crs_and_runs_as_printed(
        self, remote: Any, data: bytes, tmp_path: Path
    ) -> None:
        """Test 8, `_fetch_run`: a box in degrees on an empty cache. The hint
        repeats the box and its CRS; run as printed it exits 0, and the mesh
        then runs, offline, as the same box in the DEM's CRS meshes the file."""
        cache = tmp_path / "cache"
        box = ("--bbox", *BBOX_DEGREES, "--bbox-crs", "EPSG:4326")
        mesh = ("mesh", "--dem", KEY, "--cache", str(cache), *box, "--out")
        fetch = hinted(refusal(*mesh, str(tmp_path / "m.vtk")))
        assert fetch == ["fetch", KEY, *box, "--cache", str(cache)], fetch
        code, output = invoke(*fetch)
        assert code == 0, output
        code, output = invoke(*mesh, str(tmp_path / "m.vtk"))
        assert code == 0, output
        tif = tmp_path / "dem.tif"
        tif.write_bytes(data)
        by_path = tmp_path / "path.vtk"
        code, output = invoke(
            "mesh", "--dem", str(tif), *TestF10TheCommand.BOX_ARGS, "--out", str(by_path)
        )
        assert code == 0, output
        same_mesh(read_vtk((tmp_path / "m.vtk").read_bytes()), read_vtk(by_path.read_bytes()))

    def test_the_window_hint_names_its_crs_and_runs_as_printed(
        self, remote: Any, tmp_path: Path
    ) -> None:
        """Test 8, `_bbox_fetch`: one block fetched, then a domain whose
        window needs more. The refusal says the window's CRS (`NotCached`),
        the hint passes it as `--bbox-crs`; run as printed it exits 0, and
        the mesh then runs."""
        cache = tmp_path / "cache"
        corner = ("--bbox", str(TIE_X + 10), str(TIE_Y - 30), str(TIE_X + 20), str(TIE_Y - 20))
        code, output = invoke(
            "fetch", KEY, *corner, "--bbox-crs", PROJECTED_CRS, "--cache", str(cache)
        )
        assert code == 0, output
        domain = geojson(tmp_path / "d.geojson", TRIANGLE)
        mesh = (
            "mesh", "--dem", KEY, "--cache", str(cache), "--domain", str(domain),
            "--tolerance", "1", "--domain-max-shift", "5", "--simplify-domain",
            "--out", str(tmp_path / "m.vtk"),
        )  # fmt: skip
        message = refusal(*mesh)
        window = "500200 6599860 500280 6599900"
        assert f"for the window {window} in {PROJECTED_CRS};" in message, message
        fetch = hinted(message)
        assert fetch == [
            "fetch", KEY, "--bbox", *window.split(), "--bbox-crs", PROJECTED_CRS,
            "--cache", str(cache),
        ], fetch  # fmt: skip
        code, output = invoke(*fetch)
        assert code == 0, output
        code, output = invoke(*mesh)
        assert code == 0, output


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

    def test_no_cache_directory_is_refused_naming_the_option_and_both_variables(
        self, remote: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """No `--cache`, `RASPUTIN_DATA`, `XDG_CACHE_HOME` or `HOME`: nowhere to
        write, so a usage error that names every way to give a place."""
        no_cache_anywhere(monkeypatch)
        code, output = invoke("fetch", KEY, *TestF10TheCommand.BOX_ARGS)
        assert code == USAGE, output
        says_no_cache(output)


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
