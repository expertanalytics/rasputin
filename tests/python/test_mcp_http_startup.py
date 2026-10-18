"""Increment 69, T8: what stops ``rasputin mcp --http`` from starting, and Host and Origin.

`docs/increments/69-mcp-http.md`, sections 5.2, 5.6, 5.8, 6.1, 6.3 and 9. The
command runs in process through Typer's runner, with `uvicorn.Config` replaced
by one that stops the run: a refusal must come before uvicorn is configured,
and a missing refusal ends the test instead of serving forever. Host and
Origin go to a loopback server whose public URL is ``https://mcp.example.org``.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
import typer

pytest.importorskip("mcp")

from typer.testing import CliRunner, Result

from mcp_http_driver import (
    TOKENS,
    TOKENS_VALUE,
    HttpServer,
    auth,
    call,
    initialize,
    post,
    rasputin_server,
)
from tin_engine.cli import app

VARIABLE = "RASPUTIN_MCP_TOKENS"


class UvicornReachedError(Exception):
    """Raised where uvicorn would have been configured: the command got past its checks."""


@pytest.fixture
def run(monkeypatch: pytest.MonkeyPatch) -> Any:
    import uvicorn

    def stop(self: Any, *args: Any, **options: Any) -> None:
        raise UvicornReachedError

    monkeypatch.setattr(uvicorn.Config, "__init__", stop)

    def _run(*options: str, tokens: str | None = TOKENS_VALUE) -> Result:
        if tokens is None:
            monkeypatch.delenv(VARIABLE, raising=False)
        else:
            monkeypatch.setenv(VARIABLE, tokens)
        return CliRunner().invoke(app, ["mcp", "--http", *options])

    return _run


def refused(result: Result, *words: str) -> None:
    assert not isinstance(result.exception, UvicornReachedError), "the server started"
    assert result.exit_code == 2, (result.exit_code, result.output, result.exception)
    for word in words:
        assert word in result.output, result.output


class TestStartUpRefusals:
    def test_no_tokens_exits_2_naming_the_variable(self, run: Any) -> None:
        refused(run(tokens=None), VARIABLE)

    def test_a_bad_token_exits_2_naming_the_user_not_the_token(self, run: Any) -> None:
        result = run(tokens=f"ola={'s' * 31}")
        refused(result, VARIABLE, "ola")
        assert "s" * 31 not in result.output

    def test_a_host_off_loopback_needs_a_public_url(self, run: Any) -> None:
        refused(run("--host", "0.0.0.0"), "--public-url")

    def test_a_public_url_off_loopback_must_be_https(self, run: Any) -> None:
        refused(run("--public-url", "http://example.org"), "https")

    def test_a_shared_root_inside_the_out_root_stops_start(
        self, run: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The tile cache under ``--out-root`` would open other users' folders."""
        monkeypatch.setenv("RASPUTIN_DATA", str(tmp_path / "out" / "data"))
        refused(run("--out-root", str(tmp_path / "out")), str(tmp_path / "out"))

    def test_an_out_root_inside_a_shared_root_stops_start(
        self, run: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("RASPUTIN_DATA", str(tmp_path / "data"))
        refused(run("--out-root", str(tmp_path / "data" / "cache" / "out")), str(tmp_path / "data"))

    def test_good_settings_reach_uvicorn(self, run: Any) -> None:
        result = run("--public-url", "https://mcp.example.org", "--quota", "1M")
        assert isinstance(result.exception, UvicornReachedError), (result.output, result.exception)


def test_a_command_named_download_link_stops_start(tmp_path: Path) -> None:
    from tin_engine.mcp import http as mcp_http

    toy = typer.Typer()

    @toy.command("download-link")
    def download_link(path: str) -> None:
        """A plugin command whose tool name the server needs for itself."""

    @toy.command("other")
    def other() -> None:
        """Another command, so the app keeps its subcommands."""

    settings = mcp_http.HttpSettings(
        host="127.0.0.1",
        port=8765,
        public_url="http://127.0.0.1:8765",
        out_root=tmp_path / "out",
        quota=1024,
    )
    tokens = mcp_http.tokens_from_env({VARIABLE: TOKENS_VALUE})
    with pytest.raises(ValueError, match="download-link"):
        mcp_http.build_app(toy, settings, tokens)


# ---------------------------------------------------------------- Host, Origin, environment


@pytest.fixture(scope="module")
def server(tmp_path_factory: pytest.TempPathFactory) -> Iterator[HttpServer]:
    root = tmp_path_factory.mktemp("startup")
    with rasputin_server(root, public_url="https://mcp.example.org") as running:
        yield running


class TestHostAndOrigin:
    def test_the_public_host_passes(self, server: HttpServer) -> None:
        response = post(server, initialize(), {**auth("ola"), "Host": "mcp.example.org"})
        assert response.status_code == 200, response.text

    def test_the_public_origin_passes(self, server: HttpServer) -> None:
        headers = {**auth("ola"), "Origin": "https://mcp.example.org"}
        assert post(server, initialize(), headers).status_code == 200

    def test_another_origin_is_403(self, server: HttpServer) -> None:
        headers = {**auth("ola"), "Origin": "https://evil.example"}
        assert post(server, initialize(), headers).status_code == 403

    def test_another_host_is_refused(self, server: HttpServer) -> None:
        headers = {**auth("ola"), "Host": "evil.example"}
        assert post(server, initialize(), headers).status_code == 421


async def test_the_child_never_sees_the_tokens(
    server: HttpServer, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv(VARIABLE, TOKENS_VALUE)
    env_file = tmp_path / "env.json"
    server.set_spec(env_file=str(env_file))
    result = await call(server, "ola", "version", {})
    assert not result.is_error, result
    env = json.loads(env_file.read_text())
    assert VARIABLE not in env
    for token in TOKENS.values():
        assert all(token not in value for value in env.values())
