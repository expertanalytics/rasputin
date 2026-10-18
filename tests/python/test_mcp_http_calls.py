"""Increment 69, T5, T6, T9, T11 and T12: whole calls over HTTP.

`docs/increments/69-mcp-http.md`, sections 5.4, 5.6, 6.5, 7 and 9. One
loopback server per module with the fake child (`mcp_fake_child`), which
copies a prepared file into the caller's ``work/`` and prints its name, as a
command prints the files it writes; T9 runs the real ``python -m tin_engine``.

Bounds: the structured content stays under 64 KiB for a 2 MiB output, the
client's event limit being 1 MiB (a test bound, checked at 2 MiB); a cancelled
child is gone within 2 s, as on stdio.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import anyio
import httpx2
import pytest

pytest.importorskip("mcp")

from geotiff_fixtures import micro_tiff
from mcp_driver import Progress, wait_until_dead
from mcp_http_driver import (
    SESSION_REVISION,
    TOKENS,
    HttpServer,
    Link,
    auth,
    call,
    downloads,
    initialize,
    payload,
    post,
    rasputin_server,
    sha256,
    until_exists,
)
from tin_engine.mcp import server as mcp_server
from vtkread import read_vtk

MIB = 1024 * 1024
BIG = payload(2 * MIB, seed=5)


@pytest.fixture(scope="module")
def server(tmp_path_factory: pytest.TempPathFactory) -> Iterator[HttpServer]:
    with rasputin_server(tmp_path_factory.mktemp("calls")) as running:
        for user in TOKENS:
            running.work(user).mkdir(parents=True, exist_ok=True)
        yield running


@pytest.fixture
def big(tmp_path: Path) -> Path:
    source = tmp_path / "big.gpkg"
    source.write_bytes(BIG)
    return source


def resource_links(result: Any) -> list[Any]:
    return [b for b in result.content if b.type == "resource_link"]


def fetch(url: str) -> bytes:
    response = httpx2.get(url, timeout=60)
    assert response.status_code == 200, response
    return response.content


# ---------------------------------------------------------------- T5


class TestACallEndToEnd:
    async def test_a_2_mib_output_comes_back_as_a_link(self, server: HttpServer, big: Path) -> None:
        server.set_spec(copy=[[str(big), "a.gpkg"]], stdout=["a.gpkg"])
        result = await call(server, "ola", "version", {})
        assert not result.is_error, result
        body = result.structured_content
        assert isinstance(body, dict)
        assert body["outputs"] == ["a.gpkg"]
        (download,) = downloads(result)
        assert download["path"] == "a.gpkg"
        assert download["sha256"] == sha256(BIG)
        assert download["bytes"] == len(BIG)
        assert str(download["url"]).startswith(server.base + "/files/ola/")
        assert str(download["expires"]).endswith(("Z", "+00:00")), download["expires"]
        (link,) = resource_links(result)
        assert str(link.uri) == download["url"]
        assert link.name == "a.gpkg"
        assert link.mime_type == "application/geopackage+sqlite3"
        assert link.size == len(BIG)
        assert fetch(str(link.uri)) == BIG
        assert len(json.dumps(body).encode()) < 64 * 1024
        assert (server.work("ola") / "a.gpkg").read_bytes() == BIG
        assert (server.store("ola") / f"{sha256(BIG)}.gpkg").read_bytes() == BIG

    async def test_a_second_identical_call_adds_no_stored_file(
        self, server: HttpServer, big: Path
    ) -> None:
        server.set_spec(copy=[[str(big), "a.gpkg"]], stdout=["a.gpkg"])
        assert not (await call(server, "ola", "version", {})).is_error
        before = sorted(p.name for p in server.store("ola").iterdir())
        assert not (await call(server, "ola", "version", {})).is_error
        assert sorted(p.name for p in server.store("ola").iterdir()) == before

    async def test_only_regular_files_inside_work_are_outputs(
        self, server: HttpServer, big: Path
    ) -> None:
        """``/etc/hosts``, Kari's file and a directory in ``work/`` stay in stdout."""
        kari_file = server.work("kari") / "x"
        kari_file.write_text("kari's")
        kari_line = "../../kari/work/x"
        assert (server.work("ola") / kari_line).resolve() == kari_file.resolve()
        lines = ["a.gpkg", "/etc/hosts", kari_line, "adir"]
        server.set_spec(copy=[[str(big), "a.gpkg"]], mkdir=["adir"], stdout=lines)
        result = await call(server, "ola", "version", {})
        assert not result.is_error, result
        body = result.structured_content
        assert isinstance(body, dict)
        assert body["outputs"] == ["a.gpkg"]
        assert [d["path"] for d in downloads(result)] == ["a.gpkg"]
        assert [b.name for b in resource_links(result)] == ["a.gpkg"]
        assert body["stdout"] == lines

    async def test_a_failed_call_stores_nothing(self, server: HttpServer, tmp_path: Path) -> None:
        source = tmp_path / "fail.gpkg"
        source.write_bytes(payload(1000, seed=9))
        server.set_spec(copy=[[str(source), "fail.gpkg"]], stdout=["fail.gpkg"], exit=1)
        result = await call(server, "ola", "version", {})
        assert result.is_error
        assert not (server.store("ola") / f"{sha256(source.read_bytes())}.gpkg").exists()


# ---------------------------------------------------------------- T6


class TestTwoUsersApart:
    async def test_calls_at_once_write_to_their_own_work(
        self, server: HttpServer, tmp_path: Path
    ) -> None:
        source = tmp_path / "t6.gpkg"
        source.write_bytes(payload(5000, seed=6))
        server.set_spec(copy=[[str(source), "t6.gpkg"]], stdout=["t6.gpkg"])
        results: dict[str, Any] = {}

        async def one(user: str) -> None:
            results[user] = await call(server, user, "version", {})

        async with anyio.create_task_group() as group:
            for user in ("ola", "kari"):
                group.start_soon(one, user)
        for user in ("ola", "kari"):
            assert not results[user].is_error, results[user]
            assert (server.work(user) / "t6.gpkg").read_bytes() == source.read_bytes()
            (download,) = downloads(results[user])
            assert Link.parse(str(download["url"])).user == user

    async def test_kari_cannot_link_olas_file(self, server: HttpServer, tmp_path: Path) -> None:
        mine = server.work("ola") / "mine.gpkg"
        mine.write_bytes(payload(100, seed=7))
        for given in ("../../ola/work/mine.gpkg", str(mine)):
            result = await call(server, "kari", "download_link", {"path": given})
            assert result.is_error, f"Kari got a link to {given}"
        assert not resource_links(result)

    async def test_download_link_gives_a_fresh_link_to_ones_own_file(
        self, server: HttpServer
    ) -> None:
        data = payload(4321, seed=8)
        (server.work("ola") / "again.gpkg").write_bytes(data)
        result = await call(server, "ola", "download_link", {"path": "again.gpkg"})
        assert not result.is_error, result
        (link,) = resource_links(result)
        assert link.name == "again.gpkg"
        assert fetch(str(link.uri)) == data
        body = result.structured_content
        assert isinstance(body, dict)
        found = body["downloads"][0] if "downloads" in body else body
        assert found["sha256"] == sha256(data)

    @pytest.mark.parametrize("given", ["missing.gpkg", "adir2"], ids=["missing", "directory"])
    async def test_download_link_needs_a_regular_file(self, server: HttpServer, given: str) -> None:
        (server.work("ola") / "adir2").mkdir(exist_ok=True)
        result = await call(server, "ola", "download_link", {"path": given})
        assert result.is_error

    def test_olas_session_with_karis_token_is_refused(self, server: HttpServer) -> None:
        opened = post(server, initialize(SESSION_REVISION), auth("ola"))
        assert opened.status_code == 200, opened.text
        session = opened.headers["mcp-session-id"]
        headers = {"mcp-session-id": session, "mcp-protocol-version": SESSION_REVISION}
        note = {"jsonrpc": "2.0", "method": "notifications/initialized"}
        assert post(server, note, {**auth("ola"), **headers}).status_code == 202
        listing = {"jsonrpc": "2.0", "id": 2, "method": "tools/list"}
        assert post(server, listing, {**auth("ola"), **headers}).status_code == 200
        assert post(server, listing, {**auth("kari"), **headers}).status_code == 404


# ---------------------------------------------------------------- T9


@pytest.fixture(scope="module")
def real_server(tmp_path_factory: pytest.TempPathFactory) -> Iterator[HttpServer]:
    root = tmp_path_factory.mktemp("real")
    with rasputin_server(root, child=list(mcp_server.CHILD)) as running:
        running.work("ola").mkdir(parents=True, exist_ok=True)
        yield running


async def test_a_real_mesh_over_http_returns_a_vtk_link(real_server: HttpServer) -> None:
    (real_server.work("ola") / "micro.tif").write_bytes(micro_tiff().getvalue())
    result = await call(real_server, "ola", "mesh", {"dem": ["micro.tif"], "out": "m.vtk"})
    assert not result.is_error, result
    vtk = [d for d in downloads(result) if str(d["path"]).endswith(".vtk")]
    assert [d["path"] for d in vtk] == ["m.vtk"]
    parsed = read_vtk(fetch(str(vtk[0]["url"])))
    assert parsed is not None


# ---------------------------------------------------------------- T11

UVICORN_LOGGERS = ("uvicorn", "uvicorn.error", "uvicorn.access")


async def test_no_log_record_holds_a_token_a_signature_or_a_link(
    server: HttpServer, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.DEBUG)
    for name in UVICORN_LOGGERS:  # uvicorn's own loggers need not propagate to the root
        logging.getLogger(name).addHandler(caplog.handler)
    try:
        source = tmp_path / "log.gpkg"
        source.write_bytes(payload(2000, seed=11))
        server.set_spec(copy=[[str(source), "log.gpkg"]], stdout=["log.gpkg"])
        result = await call(server, "ola", "version", {})
        (download,) = downloads(result)
        link = Link.parse(str(download["url"]))
        assert fetch(str(download["url"])) == source.read_bytes()
    finally:
        for name in UVICORN_LOGGERS:
            logging.getLogger(name).removeHandler(caplog.handler)
    texts = [record.getMessage() for record in caplog.records]
    for text in texts:
        for token in TOKENS.values():
            assert token not in text, text
        assert link.sig not in text, text
        assert "?fn=" not in text, text
    ours = [r.getMessage() for r in caplog.records if r.name == "tin_engine.mcp.http"]
    assert any("ola" in t and "version" in t for t in ours), ours


def test_the_cli_serves_with_no_access_log(monkeypatch: pytest.MonkeyPatch) -> None:
    """uvicorn's access line would hold the query string, and so the signature."""
    import uvicorn
    from typer.testing import CliRunner

    from mcp_http_driver import TOKENS_VALUE
    from tin_engine.cli import app

    seen: dict[str, Any] = {}

    class UvicornReachedError(Exception):
        pass

    def config(self: Any, app: Any, **options: Any) -> None:
        seen.update(options)
        raise UvicornReachedError

    monkeypatch.setattr(uvicorn.Config, "__init__", config)
    monkeypatch.setenv("RASPUTIN_MCP_TOKENS", TOKENS_VALUE)
    result = CliRunner().invoke(app, ["mcp", "--http"])
    assert isinstance(result.exception, UvicornReachedError), (result.output, result.exception)
    assert seen.get("access_log") is False, seen
    assert seen.get("host") == "127.0.0.1" and seen.get("port") == 8765, seen


# ---------------------------------------------------------------- T12


async def test_progress_lines_arrive_before_the_result(
    server: HttpServer, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(mcp_server, "HEARTBEAT_SECONDS", 60.0)  # only the child's lines
    lines = ["first line", "second line", "third line"]
    server.set_spec(stderr=[[0.1, line] for line in lines])
    progress = Progress()
    result = await call(server, "ola", "version", {}, progress_callback=progress)
    assert not result.is_error, result
    assert [m for m in progress.messages if m in lines] == lines
    progress.assert_strictly_increasing()


async def test_closing_the_client_mid_run_kills_the_child(
    server: HttpServer, tmp_path: Path
) -> None:
    pid_file = tmp_path / "pid"
    server.set_spec(pid_file=str(pid_file), sleep=30.0)
    with anyio.fail_after(20):
        async with anyio.create_task_group() as group:
            group.start_soon(call, server, "ola", "version", {})
            await until_exists(pid_file)
            group.cancel_scope.cancel()
    await wait_until_dead(pid_file, within=2.0)
