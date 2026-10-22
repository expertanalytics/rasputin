"""One way to drive ``rasputin mcp --http`` in the tests: a real uvicorn on a loopback port.

`docs/increments/69-mcp-http.md`, section 9. Test support only; nothing in
`src_python` uses it. The server under test is the Starlette app a caller's
``build(port)`` returns (the suites pass `tin_engine.mcp.http.build_app`),
served by `uvicorn.Server` in a daemon thread on a socket this module binds
first, so the port is known before the app is built (links are made on the
public URL) and no other process can take it in between.

Start-up waits on uvicorn's own ``started`` flag, not on a sleep; shutdown sets
``should_exit`` and joins the thread. The server runs in this process, so a
test may patch a module the server reads at call time (the heartbeat interval)
and read the server's log records with ``caplog``.

The fake child is `mcp_fake_child`, with one spec file per server
(`HttpServer.spec`): a test writes the spec it needs before its call, and the
calls of one module run one after another, so no two tests share a spec at
once.
"""

from __future__ import annotations

import hashlib
import http.client
import json
import socket
import sys
import threading
import time
from collections.abc import AsyncIterator, Callable, Iterator
from contextlib import asynccontextmanager, contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlencode, urlsplit

import httpx2
import uvicorn
from mcp.client import Client
from mcp.client.streamable_http import streamable_http_client
from mcp_types import CallToolResult

from mcp_driver import FAKE_CHILD

#: Two users and their tokens: 40 characters each, above the 32-character floor.
TOKENS = {"ola": "o" * 20 + "la-token-0123456789", "kari": "k" * 20 + "ari-token-012345678"}
#: The value of ``RASPUTIN_MCP_TOKENS`` that configures `TOKENS`.
TOKENS_VALUE = " ".join(f"{user}={token}" for user, token in TOKENS.items())
#: A fixed signing key: the suites patch `secrets.token_bytes` while the app is
#: built, so a test can sign a link the server will accept (a past expiry, say).
KEY = bytes(range(32))
#: The headers an MCP POST needs besides auth.
ACCEPT = {"Accept": "application/json, text/event-stream", "Content-Type": "application/json"}
#: A protocol revision that still opens sessions (``Mcp-Session-Id``).
SESSION_REVISION = "2025-11-25"
#: Seconds the server is given to start; uvicorn on loopback takes about 1.5 s.
START_SECONDS = 20.0


@dataclass
class HttpServer:
    """A running server: its loopback base, the MCP endpoint, its output root and spec."""

    base: str
    out_root: Path
    spec: Path
    port: int

    @property
    def mcp(self) -> str:
        return f"{self.base}/mcp"

    def child(self) -> list[str]:
        """The ``child`` prefix: `mcp_fake_child` reading this server's spec file."""
        return [sys.executable, str(FAKE_CHILD), str(self.spec)]

    def set_spec(self, **spec: Any) -> None:
        """What the next call's fake child does (`mcp_fake_child`'s docstring lists the keys)."""
        self.spec.write_text(json.dumps(spec))

    def work(self, user: str) -> Path:
        return self.out_root / user / "work"

    def store(self, user: str) -> Path:
        return self.out_root / user / "store"


def bound_socket() -> socket.socket:
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind(("127.0.0.1", 0))
    return sock


@contextmanager
def serving(build: Callable[[int], Any], out_root: Path, spec: Path) -> Iterator[HttpServer]:
    """``build(port)`` served on loopback until the block ends."""
    sock = bound_socket()
    port = sock.getsockname()[1]
    spec.write_text("{}")
    app = build(port)
    config = uvicorn.Config(app, log_config=None, access_log=False, lifespan="on")
    server = uvicorn.Server(config)
    thread = threading.Thread(target=server.run, kwargs={"sockets": [sock]}, daemon=True)
    thread.start()
    deadline = time.monotonic() + START_SECONDS
    while not server.started:
        assert thread.is_alive(), "uvicorn stopped before it started"
        assert time.monotonic() < deadline, f"uvicorn not started within {START_SECONDS} s"
        time.sleep(0.02)
    try:
        yield HttpServer(f"http://127.0.0.1:{port}", out_root, spec, port)
    finally:
        server.should_exit = True
        thread.join(timeout=10)
        sock.close()


def auth(user: str) -> dict[str, str]:
    return {"Authorization": f"Bearer {TOKENS[user]}"}


@asynccontextmanager
async def client(server: HttpServer, user: str) -> AsyncIterator[Client]:
    """An SDK client over Streamable HTTP with ``user``'s token."""
    http = httpx2.AsyncClient(headers=auth(user), timeout=60.0)
    async with http, Client(streamable_http_client(server.mcp, http_client=http)) as c:
        yield c


async def call(
    server: HttpServer, user: str, name: str, arguments: dict[str, Any], **options: Any
) -> CallToolResult:
    """One tool call by ``user`` on a fresh connection."""
    async with client(server, user) as c:
        return await c.call_tool(name, arguments, **options)


def post(
    server: HttpServer, body: dict[str, Any], headers: dict[str, str] | None = None
) -> httpx2.Response:
    """One raw JSON-RPC POST to the MCP endpoint."""
    return httpx2.post(server.mcp, json=body, headers={**ACCEPT, **(headers or {})}, timeout=30)


def initialize(revision: str = SESSION_REVISION, request_id: int = 1) -> dict[str, Any]:
    return {
        "jsonrpc": "2.0",
        "id": request_id,
        "method": "initialize",
        "params": {
            "protocolVersion": revision,
            "capabilities": {},
            "clientInfo": {"name": "test", "version": "0"},
        },
    }


def raw_get(server: HttpServer, target: str, headers: dict[str, str] | None = None) -> int:
    """GET ``target`` (path and query) exactly as written, through `http.client`, which
    does not normalise ``..`` as httpx2 does; the status code."""
    connection = http.client.HTTPConnection("127.0.0.1", server.port, timeout=30)
    try:
        connection.request("GET", target, headers=headers or {})
        response = connection.getresponse()
        response.read()
        return response.status
    finally:
        connection.close()


@dataclass(frozen=True)
class Link:
    """A download link taken apart: ``/files/<user>/<name>?fn=&exp=&sig=``."""

    origin: str
    user: str
    name: str
    fn: str
    exp: str
    sig: str

    @classmethod
    def parse(cls, url: str) -> Link:
        parts = urlsplit(url)
        segments = parts.path.split("/")
        assert segments[-3] == "files", url
        query = {k: v[0] for k, v in parse_qs(parts.query).items()}
        assert sorted(query) == ["exp", "fn", "sig"], url
        prefix = "/".join(segments[:-3])
        return cls(
            f"{parts.scheme}://{parts.netloc}{prefix}",
            segments[-2],
            segments[-1],
            query["fn"],
            query["exp"],
            query["sig"],
        )

    def target(self, **changes: str) -> str:
        """The path and query on the server, with any part changed."""
        parts = {**self.__dict__, **changes}
        query = urlencode({"fn": parts["fn"], "exp": parts["exp"], "sig": parts["sig"]})
        return f"/files/{parts['user']}/{parts['name']}?{query}"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def downloads(result: CallToolResult) -> list[dict[str, Any]]:
    content = result.structured_content
    assert isinstance(content, dict), result
    found = content.get("downloads")
    assert isinstance(found, list), content
    return found


def error_text(result: CallToolResult) -> str:
    assert result.is_error, result
    content = result.structured_content
    assert isinstance(content, dict), result
    return str(content["error"])


def payload(size: int, seed: int = 0) -> bytes:
    """``size`` deterministic bytes that differ by ``seed``."""
    block = hashlib.sha256(seed.to_bytes(8, "big")).digest() * 2048  # 64 KiB
    return (block * (size // len(block) + 1))[:size]


def _fixed_key(original: Callable[[int | None], bytes]) -> Callable[[int | None], bytes]:
    def token_bytes(nbytes: int | None = None) -> bytes:
        return KEY if nbytes == len(KEY) else original(nbytes)

    return token_bytes


@contextmanager
def rasputin_server(
    root: Path,
    *,
    quota: int = 20 * 1024**3,
    public_url: str | None = None,
    jobs: int = 1,
    child: list[str] | None = None,
    data: Path | None = None,
) -> Iterator[HttpServer]:
    """`tin_engine.mcp.http.build_app` over the real CLI, served on loopback.

    The output root is ``root/out``; the spec file ``root/spec.json``; the fake
    child unless ``child`` is given. ``data``, when given, is ``RASPUTIN_DATA``
    while the app is built; the app copies the environment then, so it holds for
    the server's life, and the server's read roots are computed from it (its
    ``cache`` folder and its ``sources.toml``, read again when that file
    changes). The signing key is `KEY`, by
    `secrets.token_bytes` patched while the app is built.
    """
    import secrets

    import pytest

    from tin_engine.cli import app as cli_app
    from tin_engine.mcp import http as mcp_http

    out_root, spec = root / "out", root / "spec.json"
    out_root.mkdir(parents=True, exist_ok=True)
    tokens = mcp_http.tokens_from_env({"RASPUTIN_MCP_TOKENS": TOKENS_VALUE})

    def build(port: int) -> Any:
        settings = mcp_http.HttpSettings(
            host="127.0.0.1",
            port=port,
            public_url=public_url or f"http://127.0.0.1:{port}",
            out_root=out_root,
            quota=quota,
        )
        prefix = child or [sys.executable, str(FAKE_CHILD), str(spec)]
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(secrets, "token_bytes", _fixed_key(secrets.token_bytes))
            if data is not None:
                patch.setenv("RASPUTIN_DATA", str(data))
            return mcp_http.build_app(cli_app, settings, tokens, jobs=jobs, child=prefix)

    with serving(build, out_root, spec) as server:
        yield server


async def until_exists(path: Path, within: float = 20.0) -> None:
    """``path`` exists within ``within`` seconds of now, polled every 20 ms."""
    import anyio

    deadline = time.monotonic() + within
    while not await anyio.Path(path).exists():
        assert time.monotonic() < deadline, f"{path} not made within {within} s"
        await anyio.sleep(0.02)
