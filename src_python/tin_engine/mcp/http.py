"""``rasputin mcp --http``: the same tools over Streamable HTTP, results as signed download links.

The only module that imports Starlette or uvicorn. Each user has one static
bearer token, from `TOKENS_VARIABLE`; only its sha256 is kept, and every
stored digest is compared in constant time. A call runs in the user's
``work`` folder, writes only there, and reads only there or from the
server's shared data (the tile cache and the local catalogue sources). Each
output is copied to the user's ``store`` under its sha256 and handed back as a
link this process serves, signed with a key made at start and held only in
memory, so a restart voids every link. Nothing logged holds a token, a
signature or a link, and uvicorn's access log, which would hold the query
string, is off.

`docs/increments/69-mcp-http.md`, *6. Security rules (each pinned by a test, section 9)*
"""

from __future__ import annotations

import asyncio
import hashlib
import hmac
import logging
import os
import re
import secrets
import stat
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Self

import typer
import uvicorn
from mcp.server.auth.middleware.auth_context import get_access_token
from mcp.server.auth.provider import AccessToken
from mcp.server.auth.settings import AuthSettings
from mcp.server.context import ServerRequestContext
from mcp.server.transport_security import TransportSecuritySettings
from pydantic import BaseModel, ConfigDict, model_validator
from starlette.applications import Starlette
from starlette.datastructures import URL
from starlette.requests import Request
from starlette.responses import FileResponse, Response
from starlette.routing import Route

from tin_engine.mcp import links
from tin_engine.mcp.schema import Download, derive
from tin_engine.mcp.server import (
    CHILD,
    DOWNLOAD_LINK,
    INSTRUCTIONS,
    LOG,
    MIME_TYPES,
    TOKENS_VARIABLE,
    build_server,
)
from tin_engine.mcp.store import Store
from tin_engine.mcp.workspace import HTTP_REFUSED, CallScope, overlap
from tin_engine.sources import catalogue, data_root, plugin_source_files, user_source_file, xdg_base

__all__ = ["HttpSettings", "StaticTokens", "build_app", "serve_http", "tokens_from_env"]

#: Seconds a link is valid, and seconds a stored or written file is kept (both chosen;
#: `DOWNLOAD_LINK`'s description states them).
LINK_SECONDS = 3600
KEEP_SECONDS = 7 * 86_400
#: The shortest token accepted, a chosen floor: ``secrets.token_urlsafe(32)`` gives 43.
TOKEN_FLOOR = 32
USER = re.compile(r"[a-z][a-z0-9_-]{0,31}")
NAME = re.compile(r"[0-9a-f]{64}(\.[a-z0-9]{1,8})?")
FILENAME = re.compile(r"[A-Za-z0-9._-]{1,128}")
EXPIRES = re.compile(r"[0-9]{1,12}")
LOOPBACK = ("127.0.0.1", "localhost", "::1")
HTTP_INSTRUCTIONS = INSTRUCTIONS + (
    "\nOver HTTP: each call runs in your own folder; give paths relative to it, as a result's "
    "outputs are. A call reads only your folder and the server's data sources. Each output "
    "comes back as a download link valid 1 hour; download_link makes a fresh one. Files in "
    "your folder are kept 7 days. A change to the server's sources takes a restart."
)


class HttpSettings(BaseModel):
    """Where the server listens, the address clients use, the output root, the quota in bytes."""

    model_config = ConfigDict(frozen=True)

    host: str
    port: int
    public_url: str
    out_root: Path
    quota: int

    @model_validator(mode="after")
    def _public(self) -> Self:
        if URL(self.public_url).scheme != "https" and not _loopback(self.public_url):
            raise ValueError(
                f"--public-url {self.public_url} must be https unless its host is loopback; "
                "a --host off loopback needs a --public-url"
            )
        return self


def _loopback(url: str) -> bool:
    return URL(url).hostname in LOOPBACK


def tokens_from_env(environ: Mapping[str, str]) -> dict[str, bytes]:
    """``user=token`` pairs from `TOKENS_VARIABLE`, as user to the token's sha256.

    A `ValueError` naming the variable and the user, never the token, when the
    variable is missing or empty, or a pair has a bad name, a repeated name, a
    repeated token, or a token shorter than `TOKEN_FLOOR`.
    """
    found: dict[str, bytes] = {}
    pairs = environ.get(TOKENS_VARIABLE, "").split()
    if not pairs:
        raise ValueError(f"{TOKENS_VARIABLE} is missing or empty; set it to user=token pairs")
    for pair in pairs:
        user, equals, token = pair.partition("=")
        if not equals:
            raise ValueError(f"{TOKENS_VARIABLE}: an entry is not user=token")
        if not USER.fullmatch(user):
            raise ValueError(
                f"{TOKENS_VARIABLE}: user {user} is not a lowercase letter then up to 31 of "
                "a-z, 0-9, _ and -"
            )
        digest = hashlib.sha256(token.encode()).digest()
        if user in found or digest in found.values():
            raise ValueError(f"{TOKENS_VARIABLE}: user {user} repeats a name or a token")
        if len(token) < TOKEN_FLOOR:
            raise ValueError(
                f"{TOKENS_VARIABLE}: the token of user {user} is under {TOKEN_FLOOR} characters"
            )
        found[user] = digest
    return found


class StaticTokens:
    """The bearer-token check: the presented token's sha256 against every stored digest."""

    def __init__(self, tokens: Mapping[str, bytes]) -> None:
        self.tokens = dict(tokens)

    async def verify_token(self, token: str) -> AccessToken | None:
        """The user, as both ``client_id`` and ``token``, or None; never stops at a match,
        so the time taken does not depend on which user matched."""
        digest = hashlib.sha256(token.encode()).digest()
        user = None
        for name, known in self.tokens.items():
            if hmac.compare_digest(digest, known):
                user = name
        return None if user is None else AccessToken(token=user, client_id=user, scopes=[])


def _shared_roots(environ: Mapping[str, str]) -> tuple[Path, ...]:
    """The tile cache (``$RASPUTIN_DATA/cache``, else the XDG cache) and the path of every
    local catalogue source, resolved."""
    data, base = data_root(environ), xdg_base(environ, "XDG_CACHE_HOME", ".cache")
    cache = data / "cache" if data else (base / "rasputin" if base else None)
    sources = catalogue(user_source_file(environ), plugin_source_files()).sources.values()
    paths = [s.path for s in sources if s.kind == "local" and s.path is not None]
    return tuple(p.expanduser().resolve() for p in [*([cache] if cache else []), *paths])


@dataclass(frozen=True)
class _Delivery:
    """The quota and the sweep before a call; the store and the signed links after."""

    settings: HttpSettings
    key: bytes

    async def admit(self, scope: CallScope) -> str | None:
        store = Store(scope.write_root.parent)
        await store.sweep(time.time(), KEEP_SECONDS)
        if (used := await store.usage()) >= self.settings.quota:
            return (
                f"your folder holds {used} bytes, at or over the quota of "
                f"{self.settings.quota} bytes; remove files from it, or wait for the sweep "
                "of files older than 7 days at each call's start"
            )
        return None

    async def deliver(self, scope: CallScope, outputs: list[Path]) -> list[Download]:
        store, user = Store(scope.write_root.parent), str(scope.user)
        expires = int(time.time()) + LINK_SECONDS
        found = []
        for place in outputs:
            stored = await store.put(place)
            name, fn = stored.path.name, re.sub(r"[^A-Za-z0-9._-]", "_", place.name)[:128]
            sig = links.sign(self.key, user, name, fn, expires)
            found.append(
                Download(
                    path=str(place.relative_to(scope.write_root)),
                    sha256=stored.sha256,
                    bytes=stored.size,
                    url=links.url(self.settings.public_url, user, name, fn, expires, sig),
                    expires=datetime.fromtimestamp(expires, UTC).isoformat(),
                )
            )
        return found


def build_app(
    app: typer.Typer,
    settings: HttpSettings,
    tokens: Mapping[str, bytes],
    *,
    jobs: int = 1,
    child: Sequence[str] = CHILD,
) -> Starlette:
    """The Starlette app: ``/mcp`` behind the bearer check, and ``/files``.

    A `ValueError` when a command is named ``download-link`` (the server's own
    tool), or when a shared read root and the output root overlap.
    """
    if any(d.tool_name == DOWNLOAD_LINK for d in derive(app)):
        raise ValueError("a command is named download-link, the name of the server's own tool")
    out_root = settings.out_root.expanduser().resolve()
    shared = _shared_roots(os.environ)
    if (words := overlap(out_root, shared)) is not None:
        raise ValueError(words)
    for user in tokens:
        Store(out_root / user).sweep_now(time.time(), KEEP_SECONDS)
    key = secrets.token_bytes(32)
    # uvicorn writes its access line whenever that logger has a handler, access_log or not
    logging.getLogger("uvicorn.access").disabled = True

    def scope_for(ctx: ServerRequestContext[Any]) -> CallScope:
        access = get_access_token()
        if access is None or access.client_id not in tokens:
            raise RuntimeError("a call without a user")
        work = out_root / access.client_id / "work"
        work.mkdir(parents=True, exist_ok=True)
        return CallScope(
            cwd=work,
            write_root=work,
            read_roots=(work, *shared),
            user=access.client_id,
            refused=frozenset(HTTP_REFUSED),
        )

    async def download(request: Request) -> Response:
        user, name = request.path_params["user"], request.path_params["name"]
        fn, exp, sig = (request.query_params.get(k, "") for k in ("fn", "exp", "sig"))
        if user not in tokens or not all(
            p.fullmatch(v) for p, v in ((NAME, name), (FILENAME, fn), (EXPIRES, exp))
        ):
            return Response(status_code=404)
        status, path = 200, out_root / user / "store" / name
        if not links.valid(key, user, name, fn, int(exp), sig, now=int(exp) - 1):
            status = 403
        elif not links.valid(key, user, name, fn, int(exp), sig, now=time.time()):
            status = 410
        else:
            try:
                info = await asyncio.to_thread(path.lstat)
                status = 200 if stat.S_ISREG(info.st_mode) else 404
            except OSError:
                status = 404
        LOG.info("download user=%s file=%s status=%d", user, name[:12], status)
        if status != 200:
            return Response(status_code=status)
        return FileResponse(
            path,
            media_type=MIME_TYPES.get(path.suffix, "application/octet-stream"),
            filename=fn,
            headers={"Cache-Control": "private, no-store", "Referrer-Policy": "no-referrer"},
        )

    public = URL(settings.public_url)
    security = TransportSecuritySettings(
        enable_dns_rebinding_protection=True,
        allowed_hosts=[public.netloc, "127.0.0.1:*", "localhost:*", "[::1]:*"],
        allowed_origins=[
            f"{public.scheme}://{public.netloc}",
            "http://127.0.0.1:*",
            "http://localhost:*",
            "http://[::1]:*",
        ],
    )
    server = build_server(
        app,
        jobs=jobs,
        child=child,
        scope_for=scope_for,
        delivery=_Delivery(settings, key),
        instructions=HTTP_INSTRUCTIONS,
    )
    LOG.info("serving users %s at %s, reading %s", sorted(tokens), settings.public_url, shared)
    return server.streamable_http_app(
        transport_security=security,
        host=settings.host,
        auth=AuthSettings(issuer_url=settings.public_url, resource_server_url=None),
        token_verifier=StaticTokens(tokens),
        custom_starlette_routes=[Route("/files/{user}/{name}", download, methods=["GET"])],
    )


def serve_http(
    app: typer.Typer, settings: HttpSettings, tokens: Mapping[str, bytes], *, jobs: int = 1
) -> None:
    """`build_app` served by uvicorn with no access log, this module logging to stderr in UTC;
    `build_app`'s `ValueError` comes before anything is served."""
    if not LOG.handlers:
        formatter = logging.Formatter("%(asctime)sZ %(name)s %(message)s")
        formatter.converter = time.gmtime
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        LOG.addHandler(handler)
        LOG.setLevel(logging.INFO)
    served = build_app(app, settings, tokens, jobs=jobs)
    config = uvicorn.Config(served, host=settings.host, port=settings.port, access_log=False)
    uvicorn.Server(config).run()
