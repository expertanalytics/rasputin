"""`rasputin mesh --dem <catalogue key>`: meshing from the tile cache, offline (23a-1, W5, W7, W8).

`docs/increments/23-basin-scale.md`, "Windowed source reads (23a-1)", "The
CLI" and "The I/O boundary". A catalogue entry pointing at a projected test
file is put into `SOURCES` by `monkeypatch` (decided 10: geographic sources
still refuse until 15c-2), and its cache is written by `cog_fixtures` from the
file's own byte ranges, so the same file meshed by path is the oracle.

HOW THIS FILE GOES RED. `tin_engine.sources`, `cli.cache_root` and the cache
writer's `CacheManifest` are reached through fixtures, so each test errors on
its own at setup (`ModuleNotFoundError`, `AttributeError`, `ImportError`) and
the rest of `tests/python` still collects.
The two source scans (no `tin_engine.fetch` or `urllib.request` off `fetch/`,
`RASPUTIN_DATA` only in `cli.py`) are the guard 23a-2 must keep: `fetch/` is
skipped, and `cli.py` may import it only inside a function.
"""

from __future__ import annotations

import ast
import importlib
import socket
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

import tin_engine
from cli_driver import USAGE, geojson, invoke, squashed
from cog_fixtures import VARIANTS, build, write_cache
from geotiff_fixtures import TIE_X, TIE_Y
from test_cli_mesh_mosaic import field, same_mesh, terrain
from test_cli_mesh_refine import stats_row
from tin_engine.dem_input import DemRequest, open_dem
from tin_engine.io.geotiff import decode_dem
from vtkread import read_vtk

KEY = "test-utm33"
CREDIT = "Test elevation, credit line 23a-1"
PACKAGE = Path(tin_engine.__file__).resolve().parent
TILED = VARIANTS[0]


@pytest.fixture(scope="module")
def sources() -> ModuleType:
    return importlib.import_module("tin_engine.sources")


@pytest.fixture
def catalogue(sources: ModuleType, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """`SOURCES` plus a projected test entry, wherever a module bound the name."""
    entry = sources.RemoteSource(
        id=KEY,
        kind="one-cog",
        url="https://example.invalid/test-utm33.tif",
        crs="EPSG:25833",
        nodata=-9999.0,
        credit=CREDIT,
        licence_note="test fixture only",
    )
    original = sources.SOURCES
    patched = {**original, KEY: entry}
    for module in list(sys.modules.values()):
        name = getattr(module, "__name__", "") or ""
        if name.startswith("tin_engine") and getattr(module, "SOURCES", None) is original:
            monkeypatch.setattr(module, "SOURCES", patched)
    return patched


@pytest.fixture
def data() -> bytes:
    """A tiled Deflate float32 DEM, terrain enough for `--tolerance` to refine."""
    return build(TILED, array=terrain(50, 70))


@pytest.fixture
def tif(tmp_path: Path, data: bytes) -> Path:
    path = tmp_path / "dem.tif"
    path.write_bytes(data)
    return path


@pytest.fixture
def cache(tmp_path: Path, data: bytes) -> Path:
    write_cache(tmp_path / "cache", KEY, {"dem": data})
    return tmp_path / "cache"


@pytest.fixture
def no_data_root(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)


def no_network(*args: Any, **kwargs: Any) -> Any:
    raise AssertionError("rasputin mesh opened a socket")


# --------------------------------------------------------------------------
# W7: offline
# --------------------------------------------------------------------------


def forbidden_imports(tree: ast.AST, package: str, *, top_level_only: bool = False) -> set[str]:
    """`tin_engine.fetch*` and `urllib.request` imports in `tree`, relative ones resolved.

    With `top_level_only`, imports inside a function body are not looked at.
    """
    found: set[str] = set()

    def visit(node: ast.AST) -> None:
        if top_level_only and isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef):
            return
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                parent = package.rsplit(".", node.level - 1)[0] if node.level > 1 else package
                base = f"{parent}.{base}" if base else parent
            names = [base] + [f"{base}.{alias.name}" for alias in node.names]
        for name in names:
            if (
                name == "urllib.request"
                or name == "tin_engine.fetch"
                or name.startswith("tin_engine.fetch.")
            ):
                found.add(name)
        for child in ast.iter_child_nodes(node):
            visit(child)

    visit(tree)
    return found


def package_modules() -> list[tuple[str, Path]]:
    """Every module of `tin_engine` outside `fetch/`, as (dotted package, path)."""
    out = []
    for path in sorted(PACKAGE.rglob("*.py")):
        relative = path.relative_to(PACKAGE.parent)
        if relative.parts[:2] == ("tin_engine", "fetch"):
            continue
        parts = relative.with_suffix("").parts
        package = ".".join(parts if path.name == "__init__.py" else parts[:-1])
        out.append((package, path))
    return out


class TestW7Offline:
    def test_the_scan_finds_what_it_forbids(self) -> None:
        """The scan can fail: each forbidden form, planted, is found."""
        planted = (
            "import urllib.request\n"
            "from urllib import request\n"
            "from ..fetch import run\n"
            "from tin_engine.fetch.http import get\n"
            "def lazy():\n    from tin_engine import fetch\n"
        )
        found = forbidden_imports(ast.parse(planted), "tin_engine.io")
        assert {"urllib.request", "tin_engine.fetch.run", "tin_engine.fetch.http"} <= found
        assert "tin_engine.fetch" in found
        top = forbidden_imports(ast.parse(planted), "tin_engine.io", top_level_only=True)
        assert "tin_engine.fetch" in top  # via `from ..fetch import run`
        lazy = ast.parse("def f():\n    import urllib.request\n")
        assert forbidden_imports(lazy, "x", top_level_only=True) == set()

    def test_no_module_off_fetch_imports_fetch_or_urllib_request(self) -> None:
        modules = package_modules()
        assert any(path.name == "cli.py" for _, path in modules)
        offenders = {}
        for package, path in modules:
            if path.name == "cli.py" and path.parent == PACKAGE:
                continue  # its lazy import is the next test's
            found = forbidden_imports(ast.parse(path.read_text()), package)
            if found:
                offenders[str(path.relative_to(PACKAGE))] = found
        assert offenders == {}

    def test_cli_imports_fetch_only_inside_a_function(self) -> None:
        tree = ast.parse((PACKAGE / "cli.py").read_text())
        assert forbidden_imports(tree, "tin_engine", top_level_only=True) == set()
        assert "urllib.request" not in forbidden_imports(tree, "tin_engine")

    def test_a_mesh_from_the_cache_is_offline_and_equals_the_file_by_path(
        self,
        catalogue: dict[str, Any],
        cache: Path,
        tif: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        by_path, by_key = tmp_path / "path.vtk", tmp_path / "key.vtk"
        code, output = invoke("mesh", "--dem", str(tif), "--tolerance", "1", "--out", str(by_path))
        assert code == 0, output
        monkeypatch.setattr(socket, "socket", no_network)
        monkeypatch.setattr(socket, "create_connection", no_network)
        code, output = invoke(
            "mesh", "--dem", KEY, "--cache", str(cache), "--tolerance", "1", "--out", str(by_key)
        )
        assert code == 0, output
        same_mesh(read_vtk(by_key.read_bytes()), read_vtk(by_path.read_bytes()))

    def test_the_fields_name_the_key_and_its_credit(
        self, catalogue: dict[str, Any], cache: Path, tmp_path: Path
    ) -> None:
        """Increment 25, D2: a downloaded DEM's file names the dataset in
        ``dem_source`` and carries its credit and licence; the cache block
        names are the ``--stats`` row ``dem_tiles``."""
        out, md = tmp_path / "key.vtk", tmp_path / "key.md"
        code, output = invoke(
            "mesh", "--dem", KEY, "--cache", str(cache), "--out", str(out), "--stats", str(md)
        )
        assert code == 0, output
        vtk = read_vtk(out.read_bytes())
        assert "elevation_source" not in vtk.field_data
        assert field(vtk, "dem_source") == KEY
        assert field(vtk, "dem_credit") == CREDIT
        assert field(vtk, "licence_note") == "test fixture only"
        assert "cite" not in vtk.field_data  # the entry asks for none
        assert "dem_tiles" not in vtk.field_data
        assert stats_row(md.read_text(encoding="utf-8"), "dem_tiles") == "dem"


class TestOpenDemFromTheCache:
    """`DemRequest.cached` and `open_dem` below the CLI."""

    def test_open_dem_gives_the_files_pixels(
        self, catalogue: dict[str, Any], cache: Path, tif: Path
    ) -> None:
        from tin_engine.dem_input import CachedSource

        opened = open_dem(DemRequest(cached=CachedSource(source=KEY, cache=cache)))
        with tif.open("rb") as stream:
            expected = decode_dem(stream)
        assert opened.tile.meta == expected.meta
        assert np.asarray(opened.tile.array).tobytes() == np.asarray(expected.array).tobytes()
        assert KEY in opened.label

    def test_exactly_one_of_sources_and_cached(self, cache: Path, tif: Path) -> None:
        from tin_engine.dem_input import CachedSource

        cached = CachedSource(source=KEY, cache=cache)
        with pytest.raises(ValueError, match="sources"):
            DemRequest(sources=(tif,), cached=cached)
        with pytest.raises(ValueError):
            DemRequest(sources=())


# --------------------------------------------------------------------------
# W5: the CLI's NotCached message
# --------------------------------------------------------------------------


class TestW5TheFetchCommandInTheMessage:
    BBOX = ("500100", "6599900", "500400", "6599980")

    def test_with_bbox(
        self, catalogue: dict[str, Any], tmp_path: Path, data: bytes, no_data_root: None
    ) -> None:
        write_cache(tmp_path / "cache", KEY, {"dem": data}, skip={"dem": (6,)})
        code, output = invoke(
            "mesh",
            "--dem", KEY, "--cache", str(tmp_path / "cache"), "--bbox", *self.BBOX,
            "--out", str(tmp_path / "m.vtk"),
        )  # fmt: skip
        assert code == USAGE, output
        wanted = (
            f"run: rasputin fetch {KEY} --bbox {' '.join(self.BBOX)} --cache {tmp_path / 'cache'}"
        )
        assert squashed(wanted) in squashed(output), output
        assert not (tmp_path / "m.vtk").exists()

    def test_with_domain(
        self, catalogue: dict[str, Any], tmp_path: Path, data: bytes, no_data_root: None
    ) -> None:
        write_cache(tmp_path / "cache", KEY, {"dem": data}, skip={"dem": (6,)})
        ring = [(TIE_X + 150, TIE_Y - 120), (TIE_X + 400, TIE_Y - 120), (TIE_X + 300, TIE_Y - 40)]
        domain = geojson(tmp_path / "d.geojson", ring)
        code, output = invoke(
            "mesh",
            "--dem", KEY, "--cache", str(tmp_path / "cache"), "--domain", str(domain),
            "--tolerance", "1", "--out", str(tmp_path / "m.vtk"),
        )  # fmt: skip
        assert code == USAGE, output
        wanted = f"run: rasputin fetch {KEY} --domain {domain} --cache {tmp_path / 'cache'}"
        assert squashed(wanted) in squashed(output), output

    def test_the_count_of_missing_and_needed_blocks(
        self, catalogue: dict[str, Any], tmp_path: Path, data: bytes, no_data_root: None
    ) -> None:
        write_cache(tmp_path / "cache", KEY, {"dem": data}, skip={"dem": (6, 19)})
        code, output = invoke(
            "mesh",
            "--dem",
            KEY,
            "--cache",
            str(tmp_path / "cache"),
            "--out",
            str(tmp_path / "m.vtk"),
        )
        assert code == USAGE, output
        assert "2 of the 20 blocks" in output, output


# --------------------------------------------------------------------------
# W8: the cache root and the key (B7)
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cache_root() -> Any:
    return importlib.import_module("tin_engine.cli").cache_root


class TestW8CacheRoot:
    def test_the_option_wins(self, cache_root: Any) -> None:
        assert cache_root(Path("/given"), {"RASPUTIN_DATA": "/data"}) == Path("/given")

    def test_the_variable_alone_is_its_cache_directory(self, cache_root: Any) -> None:
        assert cache_root(None, {"RASPUTIN_DATA": "/data"}) == Path("/data/cache")

    def test_neither_is_none(self, cache_root: Any) -> None:
        assert cache_root(None, {}) is None

    def test_through_the_cli_the_option_wins_over_the_variable(
        self,
        catalogue: dict[str, Any],
        tmp_path: Path,
        data: bytes,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        write_cache(tmp_path / "root" / "cache", KEY, {"dem": data})
        monkeypatch.setenv("RASPUTIN_DATA", str(tmp_path / "root"))
        empty = tmp_path / "empty"
        empty.mkdir()
        code, output = invoke(
            "mesh", "--dem", KEY, "--cache", str(empty), "--out", str(tmp_path / "m.vtk")
        )
        assert code == USAGE, output
        assert "not in the cache" in output and squashed(f"--cache {empty}") in squashed(output)

    def test_through_the_cli_the_variable_alone_finds_the_cache(
        self,
        catalogue: dict[str, Any],
        tmp_path: Path,
        data: bytes,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        write_cache(tmp_path / "root" / "cache", KEY, {"dem": data})
        monkeypatch.setenv("RASPUTIN_DATA", str(tmp_path / "root"))
        code, output = invoke("mesh", "--dem", KEY, "--out", str(tmp_path / "m.vtk"))
        assert code == 0, output

    def test_neither_refuses_a_catalogue_key_naming_both(
        self, catalogue: dict[str, Any], tmp_path: Path, no_data_root: None
    ) -> None:
        code, output = invoke("mesh", "--dem", KEY, "--out", str(tmp_path / "m.vtk"))
        assert code == USAGE, output
        assert "no cache: set RASPUTIN_DATA" in output and "--cache" in output, output

    def test_a_path_needs_neither(self, tif: Path, tmp_path: Path, no_data_root: None) -> None:
        code, output = invoke("mesh", "--dem", str(tif), "--out", str(tmp_path / "m.vtk"))
        assert code == 0, output

    def test_dot_slash_glo30_is_a_path_and_glo30_a_key(
        self,
        sources: ModuleType,
        tmp_path: Path,
        data: bytes,
        no_data_root: None,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        assert "glo30" in sources.SOURCES
        monkeypatch.chdir(tmp_path)
        (tmp_path / "glo30").write_bytes(data)
        code, output = invoke("mesh", "--dem", "./glo30", "--out", str(tmp_path / "path.vtk"))
        assert code == 0, output
        code, output = invoke("mesh", "--dem", "glo30", "--out", str(tmp_path / "key.vtk"))
        assert code == USAGE, output
        assert "no cache: set RASPUTIN_DATA" in output, output

    def test_a_key_and_a_path_together_are_refused(
        self, catalogue: dict[str, Any], cache: Path, tif: Path, tmp_path: Path
    ) -> None:
        code, output = invoke(
            "mesh",
            "--dem",
            KEY,
            "--dem",
            str(tif),
            "--cache",
            str(cache),
            "--out",
            str(tmp_path / "m.vtk"),
        )
        assert code == USAGE, output
        assert KEY in output and "--dem" in output, output
        assert not (tmp_path / "m.vtk").exists()

    def test_only_cli_reads_rasputin_data(self) -> None:
        readers = sorted(
            str(path.relative_to(PACKAGE))
            for path in PACKAGE.rglob("*.py")
            if "RASPUTIN_DATA" in path.read_text()
        )
        assert readers == ["cli.py"]


class TestTheCatalogue:
    def test_the_two_ruled_sources(self, sources: ModuleType) -> None:
        assert {"anadem-v1", "glo30"} <= set(sources.SOURCES)
        for key, entry in sources.SOURCES.items():
            assert entry.id == key and entry.credit
        assert sources.SOURCES["anadem-v1"].kind == "one-cog"
        assert sources.SOURCES["glo30"].kind == "cog-tiles"

    def test_sources_imports_pydantic_only(self, sources: ModuleType) -> None:
        tree = ast.parse(Path(sources.__file__).read_text())
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.add((node.module or "").split(".")[0])
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        assert imported - {"__future__", "collections", "typing", "types"} <= {"pydantic"}
