"""The HTTP server reads its sources file again when it changes, from the next call.

`docs/increments/90-mcp-reload-sources.md`, *7. Tests (`@tester`, red first)*.
The unit cases build `tin_engine.mcp.http.SharedRoots` on a temporary tree with
an explicit environment, so the user's real sources file is never read; a
reload is counted by its one ``"read roots now"`` INFO record of logger
``tin_engine.mcp.http``. The over-HTTP cases run the loopback server with the
fake child: the read check that refuses a path outside every root runs before
the child starts, so the child's run count says whether a path was a root.
"""

from __future__ import annotations

import logging
import os
import threading
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anyio
import httpx2
import pytest

pytest.importorskip("mcp")

from dataset_fixtures import local_entry, write_sources
from mcp_http_driver import HttpServer, call, error_text, rasputin_server
from tin_engine.mcp import http as mcp_http
from tin_engine.sources import catalogue, plugin_source_files

LOGGER = "tin_engine.mcp.http"
RELOADED = "read roots now"


@dataclass(frozen=True)
class Tree:
    """Two stub rasters, a sources file under ``data``, and the server's output root."""

    root: Path
    a: Path
    b: Path
    sources: Path
    out_root: Path

    @property
    def cache(self) -> Path:
        return self.root / "data" / "cache"

    @property
    def environ(self) -> dict[str, str]:
        return {
            "RASPUTIN_DATA": str(self.root / "data"),
            "XDG_CONFIG_HOME": str(self.root / "cfg"),
            "HOME": str(self.root),
        }

    def write(self, *names: str, **more: Path) -> None:
        """The sources file naming the rasters in ``names``, and any ``id=path`` given."""
        paths = {"a": self.a, "b": self.b, **more}
        chosen = {n: paths[n] for n in (*names, *more)}
        write_sources(self.sources, {f"src_{n}": local_entry(p) for n, p in chosen.items()})


@pytest.fixture
def tree(tmp_path: Path) -> Tree:
    root = tmp_path.resolve()
    a, b = root / "datasets" / "a.tif", root / "datasets" / "b.tif"
    a.parent.mkdir()
    for raster in (a, b):
        raster.write_bytes(b"II*\x00")
    made = Tree(root, a, b, root / "data" / "sources.toml", root / "out")
    made.write("a")
    return made


def shared(tree: Tree) -> Any:
    return mcp_http.SharedRoots(tree.environ, tree.out_root)


def reloads(caplog: pytest.LogCaptureFixture) -> int:
    return sum(r.name == LOGGER and RELOADED in r.getMessage() for r in caplog.records)


def warnings(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [
        r.getMessage() for r in caplog.records if r.name == LOGGER and r.levelno == logging.WARNING
    ]


@pytest.fixture
def logged(caplog: pytest.LogCaptureFixture) -> pytest.LogCaptureFixture:
    caplog.set_level(logging.INFO, logger=LOGGER)
    return caplog


# ------------------------------------------------- without a server


class TestReload:
    async def test_an_added_source_is_a_root_from_the_next_call(
        self, tree: Tree, logged: pytest.LogCaptureFixture
    ) -> None:
        roots = shared(tree)
        assert tree.a in roots.roots and tree.b not in roots.roots
        logged.clear()
        tree.write("a", "b")
        found = await roots.current()
        assert tree.a in found and tree.b in found, found
        assert roots.roots == found
        assert reloads(logged) == 1

    async def test_a_removed_source_is_no_longer_a_root(self, tree: Tree) -> None:
        roots = shared(tree)
        tree.write("a", "b")
        assert tree.a in await roots.current()
        tree.write("b")
        found = await roots.current()
        assert tree.a not in found and tree.b in found, found

    async def test_a_broken_file_keeps_the_last_good_roots_logs_once_and_heals(
        self, tree: Tree, logged: pytest.LogCaptureFixture
    ) -> None:
        roots = shared(tree)
        logged.clear()
        tree.write("a", "b")
        good = await roots.current()
        with tree.sources.open("a", encoding="utf-8") as f:
            f.write("[source.broken\n")
        assert await roots.current() == good
        assert await roots.current() == good
        (warning,) = warnings(logged)
        assert str(tree.sources) in warning, warning
        assert "Expected ']'" in warning, warning
        tree.write("b")
        healed = await roots.current()
        assert tree.a not in healed and tree.b in healed, healed
        assert reloads(logged) == 2

    async def test_a_source_inside_the_output_root_is_refused_and_the_roots_kept(
        self, tree: Tree, logged: pytest.LogCaptureFixture
    ) -> None:
        roots = shared(tree)
        before = roots.roots
        logged.clear()
        inside = tree.out_root / "ola" / "work" / "x.tif"
        tree.write("a", inside=inside)
        assert await roots.current() == before
        (warning,) = warnings(logged)
        assert "overlap" in warning, warning
        assert str(inside) in warning, warning
        assert reloads(logged) == 0

    async def test_unchanged_content_reloads_nothing_even_when_touched(
        self, tree: Tree, logged: pytest.LogCaptureFixture
    ) -> None:
        roots = shared(tree)
        before = roots.roots
        logged.clear()
        assert await roots.current() == before
        assert await roots.current() == before
        stamp = tree.sources.stat().st_mtime + 10
        os.utime(tree.sources, (stamp, stamp))
        assert await roots.current() == before
        assert reloads(logged) == 0

    async def test_concurrent_arrivals_after_a_change_reload_once(
        self, tree: Tree, logged: pytest.LogCaptureFixture
    ) -> None:
        roots = shared(tree)
        tree.write("a", "b")
        logged.clear()
        found: list[tuple[Path, ...]] = []

        async def arrive() -> None:
            found.append(await roots.current())

        async with anyio.create_task_group() as group:
            for _ in range(10):
                group.start_soon(arrive)
        assert len(found) == 10
        assert all(tree.b in f for f in found), found
        assert reloads(logged) == 1

    async def test_a_missing_file_leaves_the_cache_and_the_plugins_sources(
        self, tree: Tree
    ) -> None:
        """No user sources, as the CLI reads a missing file; the plugins' local sources
        stay, so the expectation holds once a source plugin is installed."""
        roots = shared(tree)
        tree.sources.unlink()
        plugins = catalogue(None, plugin_source_files()).sources.values()
        local = [s.path for s in plugins if s.kind == "local" and s.path is not None]
        expected = tuple(p.expanduser().resolve() for p in [tree.cache, *local])
        assert await roots.current() == expected

    async def test_a_cancelled_reload_is_retried_by_the_next_call(
        self, tree: Tree, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The first reload after a change is cancelled while it reads the catalogue; the
        content is not taken as judged, so the next call reloads and sees the change."""
        roots = shared(tree)
        entered, release = threading.Event(), threading.Event()
        read = mcp_http._shared_roots

        def held(environ: Any) -> tuple[Path, ...]:
            entered.set()
            release.wait(timeout=10)
            return read(environ)

        monkeypatch.setattr(mcp_http, "_shared_roots", held)
        tree.write("a", "b")
        with anyio.fail_after(10):
            async with anyio.create_task_group() as group:
                group.start_soon(roots.current)
                assert await anyio.to_thread.run_sync(entered.wait, 10)
                group.cancel_scope.cancel()
        release.set()
        monkeypatch.setattr(mcp_http, "_shared_roots", read)
        assert tree.b in await roots.current()


# ------------------------------------------------- over HTTP


@pytest.fixture
def server(tree: Tree) -> Iterator[HttpServer]:
    with rasputin_server(tree.root / "server", data=tree.root / "data", jobs=1) as running:
        running.work("ola").mkdir(parents=True, exist_ok=True)
        yield running


def mesh_reading(raster: Path) -> dict[str, Any]:
    return {"out": "mesh.vtk", "dem": [str(raster)]}


class TestOverHttp:
    async def test_a_changed_sources_file_moves_the_read_roots(
        self, server: HttpServer, tree: Tree
    ) -> None:
        tree.write("a", "b")
        added = tree.root / "added-count"
        server.set_spec(count_file=str(added))
        result = await call(server, "ola", "mesh", mesh_reading(tree.b))
        assert not result.is_error, result
        assert added.read_text().count("\n") == 1
        tree.write("b")
        removed = tree.root / "removed-count"
        server.set_spec(count_file=str(removed))
        result = await call(server, "ola", "mesh", mesh_reading(tree.a))
        assert "outside your folder" in error_text(result)
        assert not removed.exists(), "the child ran on a source the file no longer names"

    async def test_a_link_made_before_a_reload_still_downloads(
        self, server: HttpServer, tree: Tree
    ) -> None:
        data = b"kept across a reload"
        (server.work("ola") / "kept.gpkg").write_bytes(data)
        made = await call(server, "ola", "download_link", {"path": "kept.gpkg"})
        assert not made.is_error, made
        content = made.structured_content
        assert isinstance(content, dict), made
        url = str(content["downloads"][0]["url"] if "downloads" in content else content["url"])
        tree.write("a", "b")
        server.set_spec()
        assert not (await call(server, "ola", "version", {})).is_error
        response = httpx2.get(url, timeout=30)
        assert response.status_code == 200, response
        assert response.content == data
