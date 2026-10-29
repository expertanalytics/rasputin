"""``rasputin mcp``: the commands served over the Model Context Protocol, each call a child process.

The package's only I/O and only async code: the child ``python -m tin_engine
<argv>``, its stdout and stderr read as data, the record and stats sidecars
beside ``out`` (removed after the read), progress notifications carrying
elapsed seconds, a jobs cap, and the kill on cancel. Nothing here imports
`_core`: the kernel runs in the child, behind the CLI, which keeps every rule
it has. A child per call also means a cancelled call can be killed mid-C++,
two calls never share memory or stdout, a crash in the kernel ends only the
child, and a large run's memory is returned when it exits. This process's
stdout is the wire, so nothing here prints.

`docs/increments/37-mcp-server.md`, *6.4 Execution: a child process per call*
"""

from __future__ import annotations

import asyncio
import base64
import json
import logging
import os
import re
import secrets
import shlex
import signal
import stat
import sys
import time
from collections.abc import AsyncIterator, Awaitable, Callable, Sequence
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Protocol

import anyio
import mcp_types as types
import typer
from anyio import Path as AsyncPath
from mcp.server.context import ServerRequestContext
from mcp.server.lowlevel import Server
from mcp.server.stdio import stdio_server
from mcp.shared.exceptions import MCPError
from pydantic import BaseModel, ValidationError

from tin_engine import installed_version
from tin_engine.mcp.schema import DerivedCommand, Download, RunResult, derive
from tin_engine.mcp.workspace import CallScope, linked, refusal, uri
from tin_engine.sources import OUT_ROOT_VARIABLE

__all__ = [
    "DEFAULT_OUT_ROOT",
    "HEARTBEAT_SECONDS",
    "Delivery",
    "RunResult",
    "build_server",
    "serve",
]

#: Seconds between heartbeats while the child prints nothing; read at every tick.
HEARTBEAT_SECONDS = 5.0
#: The production child: the same interpreter and installed `tin_engine` as the server.
CHILD = (sys.executable, "-m", "tin_engine")
#: Typer's plain refusals and no colour in the child, so its stderr parses as text.
CHILD_ENV = {"TYPER_USE_RICH": "0", "NO_COLOR": "1", "TERM": "dumb"}
#: The output root under the home folder when none is given: a chosen name.
DEFAULT_OUT_ROOT = "rasputin-output"
#: The users' bearer tokens over HTTP; never passed to a child.
TOKENS_VARIABLE = "RASPUTIN_MCP_TOKENS"
#: The tool, served over HTTP only, that stores a file of the caller's and links it.
DOWNLOAD_LINK = "download_link"
#: The tool, served over HTTP only, that stores a file sent as a tool argument and links it;
#: its limit is in bytes after decoding, a chosen scale (a catchment polygon), and the
#: longest base64 that can decode to at most that many bytes.
UPLOAD = "upload"
UPLOAD_BYTES = 2 * 1024 * 1024
UPLOAD_BASE64 = (UPLOAD_BYTES + 2) // 3 * 4
UPLOAD_NAME = re.compile(r"[A-Za-z0-9_-][A-Za-z0-9._-]{0,127}")
#: The HTTP transport's log: one line per call; stdio logs nothing.
LOG = logging.getLogger("tin_engine.mcp.http")
MIME_TYPES = {
    ".gpkg": "application/geopackage+sqlite3",
    ".vtk": "application/octet-stream",
    ".ply": "application/octet-stream",
    ".geojson": "application/geo+json",
    ".json": "application/json",
    ".md": "text/markdown",
}
#: Shown to the model at initialize.
INSTRUCTIONS = (
    "Call sources first: it lists the elevation and land-cover sources this machine knows, "
    "local ones read in place and remote ones read from a cache, each with its coverage as a "
    "longitude/latitude box, cell size in metres, CRS, height unit, tile count, description "
    "and credit. Give it bbox (west south east north, degrees) to see what share of your area "
    "each source covers.\nChoose a source by its coverage of your area, its cell size and "
    "what its description says the heights represent.\nPass the chosen source's id as dem to "
    "mesh, catchment and station_catchments. A local source needs no fetch. A remote source "
    "is read from the cache: fetch the area first; a run that lacks data says what to fetch."
    "\nsources also says which sources file and cache it read, and where a local dataset "
    "would be registered.\n"
    "sources also lists vector datasets (polygons, lines, points) with their path, layer, "
    "geometry, feature count, attributes and, for land cover, the class values. Each source's "
    "fits names the tools and arguments it can be given to, and for land cover the class file "
    "to pass as land_cover_classes; what the features are is in its description. Pass its "
    "path to those arguments, and its layer as land_cover_layer or lakes_layer. To mesh: "
    "class polygons go to land_cover with the fit's land_cover_classes, lakes and glaciers "
    "from another source to prefer, rivers and roads to lines; edges names the borders and "
    "lines the mesh follows (all, water, a category, a class code, a lines key). land_cover "
    "alone gives each class's share of every triangle and the cells; with edges naming a "
    "class, add shares for them.\n"
    "Where a polygon source's description says each unit names its downstream unit, "
    "network_catchment joins the units upstream of a unit, or of a seed at a unit's outlet; "
    "check_dem compares the result with the terrain route of catchment.\n"
    "catchment and network_catchment write a polygon that mesh reads as domain; with "
    "simplify_domain its outline follows the triangle budget.\n"
    "Output folders are created as needed, inside the output root that sources reports as "
    "out_root; give output paths in full: a path outside it, through .. or a link too, is "
    "refused.\n"
    """\
Paths are on the machine this server runs on, as for the rasputin command line.
mesh writes a mesh file at out and returns its run record; the mesh itself never travels here.
Write a .gpkg to inspect a mesh with GIS tools; .vtk is for ParaView.
A refusal's words say what to change; the exit code is the command's.
A cancelled mesh may leave an unreadable .vtk or .ply at out; rerun before reading it."""
)


class _Progress:
    """One call's clock and sender: each notification carries elapsed seconds, increasing."""

    def __init__(self, ctx: ServerRequestContext[Any]) -> None:
        self.began = time.monotonic()
        self.token = ctx.meta.get("progress_token") if ctx.meta else None
        # over HTTP a notification reaches the client only on its request's stream
        self.related = None if ctx.request_id is None else str(ctx.request_id)
        self.session = ctx.session
        self.last = -1.0
        self.lock = anyio.Lock()

    def elapsed(self) -> float:
        return time.monotonic() - self.began

    async def send(self, message: str) -> None:
        if self.token is None:
            return
        async with self.lock:  # one sender per call, so values leave in order
            self.last = max(round(self.elapsed(), 3), self.last + 0.001)
            await self.session.send_progress_notification(
                self.token, self.last, None, message, related_request_id=self.related
            )

    @asynccontextmanager
    async def heartbeat(self, phrase: str) -> AsyncIterator[None]:
        """``"<phrase>, N s"`` every `HEARTBEAT_SECONDS` while the body runs; nothing without a
        progress token."""

        async def tick() -> None:
            while True:
                await asyncio.sleep(HEARTBEAT_SECONDS)
                await self.send(f"{phrase}, {self.elapsed():.1f} s")

        task = asyncio.create_task(tick()) if self.token is not None else None
        try:
            yield
        finally:
            if task is not None:
                task.cancel()


class Delivery(Protocol):
    """How a transport admits a call and hands back its outputs (stdio has none: file URIs)."""

    async def admit(self, scope: CallScope) -> str | None:
        """The refusal of a call before it runs (a quota, say), None when it may run."""

    async def deliver(self, scope: CallScope, outputs: list[Path]) -> list[Download]:
        """One `Download` per output, in order."""

    async def receive(self, scope: CallScope, name: str, data: bytes) -> Path:
        """``data`` written as the upload ``name`` in the caller's folder; its place."""


async def _run_child(
    argv: Sequence[str], progress: _Progress, scope: CallScope
) -> tuple[int, list[str], list[str]]:
    """Run ``argv`` in ``scope.cwd`` to the end, or kill it when cancelled; each stderr line
    is a notification. The child's environment names the write root, which its ``sources``
    listing reports, and holds no token."""
    env = {k: v for k, v in os.environ.items() if k != TOKENS_VARIABLE}
    process = await asyncio.create_subprocess_exec(
        *argv,
        stdin=asyncio.subprocess.DEVNULL,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=scope.cwd,
        env={**env, **CHILD_ENV, OUT_ROOT_VARIABLE: str(scope.write_root)},
    )
    assert process.stdout is not None and process.stderr is not None  # both are PIPE
    stderr: list[str] = []

    async def read_stderr(pipe: asyncio.StreamReader) -> None:
        async for raw in pipe:
            stderr.append(raw.decode(errors="replace").rstrip("\r\n"))
            await progress.send(stderr[-1])

    try:
        stdout, _ = await asyncio.gather(process.stdout.read(), read_stderr(process.stderr))
        code = await process.wait()
    finally:
        if process.returncode is None:
            process.kill()
            with anyio.CancelScope(shield=True):
                await process.wait()
    return code, stdout.decode(errors="replace").splitlines(), stderr


def _error(
    exit_code: int | None, error: str, stderr: list[str] | None = None
) -> types.CallToolResult:
    body: dict[str, Any] = {"exit_code": exit_code, "error": error}
    if stderr is not None:
        body["stderr"] = stderr
    return types.CallToolResult(
        content=[types.TextContent(type="text", text=error)], structured_content=body, is_error=True
    )


def _refusal(code: int, stderr: list[str]) -> str:
    """The CLI's own words: its ``Error:`` line(s) without the usage reminder."""
    if code < 0:
        return f"killed by signal {-code} ({signal.Signals(-code).name})"
    words = [line for line in stderr if line and not line.startswith(("Usage:", "Try '"))]
    errors = [line.removeprefix("Error: ") for line in words if line.startswith("Error: ")]
    return f"exit {code}: {' '.join(errors or words)}"


async def _read(path: Path | None) -> str | None:
    return await AsyncPath(path).read_text() if path and await AsyncPath(path).is_file() else None


def _at(scope: CallScope, path: Path) -> Path:
    return path if scope.cwd is None else scope.cwd / path


async def _outputs(
    stdout: list[str], hidden: list[Path], scope: CallScope
) -> list[tuple[str, Path]]:
    """Each output as (path shown, place): the stdout lines resolved from ``scope.cwd``, less
    the sidecars, naming anything that exists; with a read limit, only a regular file
    inside the write root, shown relative to it."""
    taken = [await AsyncPath(p).resolve() for p in hidden]
    outputs = []
    for line in stdout:
        try:  # a line that cannot name a file (too long, say) raises before 3.14
            if not line.strip():
                continue
            place = Path(await AsyncPath(_at(scope, Path(line))).resolve())
            if scope.read_roots is None:
                kept = await AsyncPath(line).exists() and place not in taken
            else:
                inside = place.is_relative_to(scope.write_root)
                kept = inside and await AsyncPath(place).is_file() and place not in taken
            if kept:
                shown = (
                    line if scope.read_roots is None else str(place.relative_to(scope.write_root))
                )
                outputs.append((shown, place))
        except OSError:
            continue
    return outputs


def _link(path: Path, uri: str, size: int | None = None) -> types.ResourceLink:
    return types.ResourceLink(
        type="resource_link",
        name=path.name,
        uri=uri,
        mime_type=MIME_TYPES.get(path.suffix),
        size=size,
    )


async def _call(
    tool: DerivedCommand,
    arguments: BaseModel,
    prefix: Sequence[str],
    progress: _Progress,
    scope: CallScope,
    delivery: Delivery | None,
) -> types.CallToolResult:
    """One validated call: the sidecars, the child, the shaped result."""
    sidecars: list[Path] = []
    record_file = stats_file = None
    update: dict[str, str] = {}
    if tool.tool_name == "mesh":
        out, token = _at(scope, Path(getattr(arguments, "out"))), secrets.token_hex(4)  # noqa: B009
        if "record_path" in arguments.model_fields_set:
            record_file = _at(scope, Path(getattr(arguments, "record_path")))  # noqa: B009
        else:
            record_file = out.parent / f".{out.name}.record-{token}.json"
            update["record_path"] = str(record_file)
            sidecars.append(record_file)
        if getattr(arguments, "stats", None) == "-":
            stats_file = out.parent / f".{out.name}.stats-{token}.md"
            update["stats"] = str(stats_file)
            sidecars.append(stats_file)
    argv = tool.argv(arguments.model_copy(update=update))
    try:
        began = time.monotonic()
        code, stdout, stderr = await _run_child([*prefix, *argv], progress, scope)
        seconds = round(time.monotonic() - began, 3)
        if code != 0:
            return _error(code, _refusal(code, stderr), stderr)
        record = await _read(record_file)
        outputs = await _outputs(stdout, sidecars, scope)
        places, shown = [p for _, p in outputs], [Path(s) for s, _ in outputs]
        downloads = None if delivery is None else await delivery.deliver(scope, places)
        result = RunResult(
            command=shlex.join(["rasputin", *argv]),
            exit_code=code,
            seconds=seconds,
            outputs=[s for s, _ in outputs],
            stdout=stdout,
            stderr=stderr,
            record=json.loads(record) if record else None,
            stats=await _read(stats_file),
            downloads=downloads,
        )
        if downloads is None:
            links = [_link(s, p.as_uri()) for s, p in zip(shown, places, strict=True)]
        else:
            links = [_link(s, d.url, d.bytes) for s, d in zip(shown, downloads, strict=True)]
    finally:
        with anyio.CancelScope(shield=True):
            for sidecar in sidecars:
                await AsyncPath(sidecar).unlink(missing_ok=True)
    body = result.model_dump(exclude={"downloads"} if delivery is None else None)
    return types.CallToolResult(
        content=[types.TextContent(type="text", text=json.dumps(body)), *links],
        structured_content=body,
    )


#: `DOWNLOAD_LINK` as a tool: one path in, one `Download` out.
DOWNLOAD_TOOL = types.Tool(
    name=DOWNLOAD_LINK,
    description="A fresh download link to a file in your folder (a path as in a result's "
    "outputs). Links last 1 hour; stored files are kept 7 days.",
    input_schema={
        "type": "object",
        "properties": {"path": {"type": "string", "description": "A file in your folder."}},
        "required": ["path"],
        "additionalProperties": False,
    },
    output_schema=Download.model_json_schema(),
)


async def _inputs(scope: CallScope, tool: DerivedCommand, arguments: BaseModel) -> BaseModel | str:
    """Over HTTP, the arguments with each link of the caller's own replaced by the stored
    file's relative path, the file touched so it is kept; or the refusal in words: a URI as
    an output, a URI not the caller's link, a stored file no longer kept."""
    if scope.link_prefix is None:
        return arguments
    for argument, given in tool.writes(arguments):
        if uri(given):
            return (
                f"{argument}: a link names a stored file, which no call writes; give a path "
                "in your folder."
            )
    mapping = {}
    for argument, given in tool.reads(arguments):
        if isinstance(found := linked(scope, argument, given), str):
            return found
        if found is not None:
            try:  # a link or a swept file is not kept; a kept one is touched, so stays kept
                kept = stat.S_ISREG((await AsyncPath(found.place).lstat()).st_mode)
                if kept:
                    await anyio.to_thread.run_sync(os.utime, found.place)
            except OSError:
                kept = False
            if not kept:
                return (
                    f"{argument}: the file this link names is no longer kept (files are kept "
                    "7 days after their last use); upload it again, or rerun the call that made it."
                )
            mapping[given] = found.text
    return tool.substitute(arguments, mapping) if mapping else arguments


def _linked_result(place: Path, download: Download) -> types.CallToolResult:
    body = download.model_dump()
    return types.CallToolResult(
        content=[
            types.TextContent(type="text", text=json.dumps(body)),
            _link(place, download.url, download.bytes),
        ],
        structured_content=body,
    )


async def _download_link(
    scope: CallScope, arguments: dict[str, Any], delivery: Delivery
) -> types.CallToolResult:
    """A regular file inside the caller's write root, stored and linked afresh."""
    given = arguments.get("path")
    if not isinstance(given, str) or set(arguments) != {"path"}:
        return _error(None, "path: give one path, of a file in your folder")
    home = Path(await AsyncPath(given).expanduser())
    place = Path(await AsyncPath(_at(scope, home)).resolve())
    if not (place.is_relative_to(scope.write_root) and await AsyncPath(place).is_file()):
        return _error(None, f"path: {given} is not a file inside your folder {scope.write_root}")
    if (words := await delivery.admit(scope)) is not None:
        return _error(None, words)
    (download,) = await delivery.deliver(scope, [place])
    return _linked_result(place, download)


#: `UPLOAD` as a tool: a name and the content, one `Download` out.
UPLOAD_TOOL = types.Tool(
    name=UPLOAD,
    description="Store a file of your own, sent as text or as base64, up to 2 MiB, in "
    "uploads/ of your folder; returns its download link, which any tool takes as an input.",
    input_schema={
        "type": "object",
        "properties": {
            "name": {"type": "string", "description": "The file name; its suffix is the format."},
            "text": {"type": "string", "description": "The content as text (UTF-8)."},
            "base64": {"type": "string", "description": "The content as padded base64."},
        },
        "required": ["name"],
        "additionalProperties": False,
    },
    output_schema=Download.model_json_schema(),
    annotations=types.ToolAnnotations(
        read_only_hint=False, destructive_hint=False, idempotent_hint=True, open_world_hint=False
    ),
)


async def _upload(
    scope: CallScope, arguments: dict[str, Any], delivery: Delivery
) -> types.CallToolResult:
    """The content checked (name, one of text or base64, size before and after decoding),
    admitted, written to ``uploads/<name>``, stored and linked; nothing written when refused."""
    name, given = arguments.get("name"), arguments.get("text", arguments.get("base64"))
    if not isinstance(name, str) or not UPLOAD_NAME.fullmatch(name):
        return _error(
            None,
            "name: give a file name of letters, digits, '.', '_' and '-', not starting with "
            "'.', up to 128 characters.",
        )
    if not isinstance(given, str) or set(arguments) not in ({"name", "text"}, {"name", "base64"}):
        return _error(None, "give the file as text or as base64, one of the two.")
    data: bytes | None = None  # base64 too long to decode: refused by its size, near enough
    if "text" in arguments:
        try:
            data = given.encode()
        except UnicodeEncodeError:  # a lone surrogate has no UTF-8 form
            return _error(None, "text: holds a character with no UTF-8 form; send UTF-8 text.")
    elif len(given) <= UPLOAD_BASE64:
        try:
            data = base64.b64decode(given, validate=True)
        except ValueError:
            return _error(None, "base64: not valid base64 (A-Z, a-z, 0-9, + and /, padded with =).")
    size = len(given) * 3 // 4 if data is None else len(data)
    if data is None or size > UPLOAD_BYTES:
        return _error(
            None,
            f"the file is {size} bytes, over the upload limit of {UPLOAD_BYTES} (2 MiB); a "
            "larger input must be one of the server's data sources.",
        )
    if (words := await delivery.admit(scope)) is not None:
        return _error(None, words)
    try:
        place = await delivery.receive(scope, name, data)
    except ValueError as exc:  # its text names server paths: for the log, not the client
        LOG.warning("upload user=%s refused by the store: %s", scope.user, exc)
        return _error(None, "the server could not store the upload; ask its operator.")
    (download,) = await delivery.deliver(scope, [place])
    LOG.info("upload user=%s file=%s bytes=%d", scope.user, download.sha256[:12], download.bytes)
    return _linked_result(place, download)


def build_server(
    app: typer.Typer,
    *,
    jobs: int = 1,
    child: Sequence[str] = CHILD,
    out_root: Path | None = None,
    scope_for: Callable[[ServerRequestContext[Any]], Awaitable[CallScope]] | None = None,
    delivery: Delivery | None = None,
    instructions: str = INSTRUCTIONS,
) -> Server[Any]:
    """The server over ``app``'s commands; ``child`` is the argv prefix a command is run with.

    ``scope_for`` gives each call its scope; by default a call may write only inside
    ``out_root`` (default ``~/rasputin-output``), resolved once here, neither made nor
    checked until a call writes, and reads anything. With a ``delivery``, a call must be
    admitted by it before and after taking a job slot, its outputs come back as its
    downloads, and `DOWNLOAD_LINK` is a tool.
    """
    root = (out_root or Path.home() / DEFAULT_OUT_ROOT).expanduser().resolve()
    stdio = CallScope(cwd=None, write_root=root, read_roots=None, user=None)

    async def stdio_scope(ctx: ServerRequestContext[Any]) -> CallScope:
        return stdio

    scope_of = scope_for or stdio_scope
    tools = {d.tool_name: d for d in derive(app)}
    semaphore, prefix = anyio.Semaphore(jobs), list(child)
    extra = [] if delivery is None else [DOWNLOAD_TOOL, UPLOAD_TOOL]

    async def list_tools(ctx: ServerRequestContext[Any], params: Any) -> types.ListToolsResult:
        return types.ListToolsResult(tools=[*(d.tool() for d in tools.values()), *extra])

    async def call_tool(
        ctx: ServerRequestContext[Any], params: types.CallToolRequestParams
    ) -> types.CallToolResult:
        if delivery is not None and params.name == DOWNLOAD_LINK:
            return await _download_link(await scope_of(ctx), params.arguments or {}, delivery)
        if delivery is not None and params.name == UPLOAD:
            return await _upload(await scope_of(ctx), params.arguments or {}, delivery)
        tool = tools.get(params.name)
        if tool is None:
            raise MCPError(types.INVALID_PARAMS, f"Unknown tool: {params.name}")
        try:
            arguments = tool.model.model_validate(params.arguments or {})
        except ValidationError as exc:
            return _error(None, "; ".join(f"{e['loc'][0]}: {e['msg']}" for e in exc.errors()))
        scope = await scope_of(ctx)
        if isinstance(inputs := await _inputs(scope, tool, arguments), str):
            return _error(None, inputs)
        arguments = inputs

        async def check() -> str | None:
            given, writes, reads = (
                arguments.model_dump(),
                tool.writes(arguments),
                tool.reads(arguments),
            )
            words = await refusal(scope, tool.tool_name, given, writes, reads)
            return words if words or delivery is None else await delivery.admit(scope)

        if (words := await check()) is not None:
            return _error(None, words)
        progress = _Progress(ctx)
        try:
            async with progress.heartbeat("waiting for a running job"):
                await semaphore.acquire()
            try:
                if (words := await check()) is not None:  # the tree may change while waiting
                    return _error(None, words)
                async with progress.heartbeat("running"):
                    result = await _call(tool, arguments, prefix, progress, scope, delivery)
            finally:
                semaphore.release()
        except Exception as exc:  # the server's own fault, in words
            return _error(None, f"{type(exc).__name__}: {exc}")
        if scope.user is not None:
            body = result.structured_content or {}
            size = sum(d["bytes"] for d in body.get("downloads") or [])
            LOG.info(
                "call user=%s tool=%s exit=%s seconds=%s bytes=%d",
                *(scope.user, tool.tool_name, body.get("exit_code"), body.get("seconds"), size),
            )
        return result

    return Server(
        "rasputin",
        version=installed_version(),
        instructions=instructions,
        on_list_tools=list_tools,
        on_call_tool=call_tool,
    )


async def serve(app: typer.Typer, jobs: int = 1, out_root: Path | None = None) -> None:
    """``rasputin mcp``: the server on this process's stdin and stdout."""
    server = build_server(app, jobs=jobs, out_root=out_root)
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())
