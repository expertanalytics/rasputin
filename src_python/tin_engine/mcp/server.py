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
import json
import os
import secrets
import shlex
import signal
import sys
import time
from collections.abc import AsyncIterator, Sequence
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

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
from tin_engine.mcp.schema import DerivedCommand, RunResult, derive

__all__ = ["HEARTBEAT_SECONDS", "RunResult", "build_server", "serve"]

#: Seconds between heartbeats while the child prints nothing; read at every tick.
HEARTBEAT_SECONDS = 5.0
#: The production child: the same interpreter and installed `tin_engine` as the server.
CHILD = (sys.executable, "-m", "tin_engine")
#: Typer's plain refusals and no colour in the child, so its stderr parses as text.
CHILD_ENV = {"TYPER_USE_RICH": "0", "NO_COLOR": "1", "TERM": "dumb"}
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
            await self.session.send_progress_notification(self.token, self.last, None, message)

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


async def _run_child(argv: Sequence[str], progress: _Progress) -> tuple[int, list[str], list[str]]:
    """Run ``argv`` to the end, or kill it when cancelled; each stderr line is a notification."""
    process = await asyncio.create_subprocess_exec(
        *argv,
        stdin=asyncio.subprocess.DEVNULL,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        env={**os.environ, **CHILD_ENV},
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


async def _outputs(stdout: list[str], hidden: list[Path]) -> list[str]:
    """The stdout lines naming something that exists after the run, less the sidecars."""
    taken = [await AsyncPath(p).resolve() for p in hidden]
    outputs = []
    for line in stdout:
        if (
            line.strip()
            and await AsyncPath(line).exists()
            and (await AsyncPath(line).resolve() not in taken)
        ):
            outputs.append(line)
    return outputs


async def _call(
    tool: DerivedCommand, arguments: BaseModel, prefix: Sequence[str], progress: _Progress
) -> types.CallToolResult:
    """One validated call: the sidecars, the child, the shaped result."""
    argv = tool.argv(arguments)
    sidecars: list[Path] = []
    record_file = stats_file = None
    if tool.tool_name == "mesh":
        out, token = Path(getattr(arguments, "out")), secrets.token_hex(4)  # noqa: B009
        if "record_path" in arguments.model_fields_set:
            record_file = Path(getattr(arguments, "record_path"))  # noqa: B009
        else:
            record_file = out.parent / f".{out.name}.record-{token}.json"
            argv += ["--record", str(record_file)]
            sidecars.append(record_file)
        if getattr(arguments, "stats", None) == "-":
            stats_file = out.parent / f".{out.name}.stats-{token}.md"
            argv[argv.index("--stats") + 1] = str(stats_file)
            sidecars.append(stats_file)
    try:
        began = time.monotonic()
        code, stdout, stderr = await _run_child([*prefix, *argv], progress)
        seconds = round(time.monotonic() - began, 3)
        if code != 0:
            return _error(code, _refusal(code, stderr), stderr)
        record = await _read(record_file)
        result = RunResult(
            command=shlex.join(["rasputin", *argv]),
            exit_code=code,
            seconds=seconds,
            outputs=await _outputs(stdout, sidecars),
            stdout=stdout,
            stderr=stderr,
            record=json.loads(record) if record else None,
            stats=await _read(stats_file),
        )
        links = [
            types.ResourceLink(
                type="resource_link",
                name=Path(p).name,
                uri=(await AsyncPath(p).resolve()).as_uri(),
                mime_type=MIME_TYPES.get(Path(p).suffix),
            )
            for p in result.outputs
        ]
    finally:
        with anyio.CancelScope(shield=True):
            for sidecar in sidecars:
                await AsyncPath(sidecar).unlink(missing_ok=True)
    body = result.model_dump()
    return types.CallToolResult(
        content=[types.TextContent(type="text", text=json.dumps(body)), *links],
        structured_content=body,
    )


def build_server(app: typer.Typer, *, jobs: int = 1, child: Sequence[str] = CHILD) -> Server[Any]:
    """The server over ``app``'s commands; ``child`` is the argv prefix a command is run with."""
    tools = {d.tool_name: d for d in derive(app)}
    semaphore, prefix = anyio.Semaphore(jobs), list(child)

    async def list_tools(ctx: ServerRequestContext[Any], params: Any) -> types.ListToolsResult:
        return types.ListToolsResult(tools=[d.tool() for d in tools.values()])

    async def call_tool(
        ctx: ServerRequestContext[Any], params: types.CallToolRequestParams
    ) -> types.CallToolResult:
        tool = tools.get(params.name)
        if tool is None:
            raise MCPError(types.INVALID_PARAMS, f"Unknown tool: {params.name}")
        try:
            arguments = tool.model.model_validate(params.arguments or {})
        except ValidationError as exc:
            return _error(None, "; ".join(f"{e['loc'][0]}: {e['msg']}" for e in exc.errors()))
        progress = _Progress(ctx)
        try:
            async with progress.heartbeat("waiting for a running job"):
                await semaphore.acquire()
            try:
                async with progress.heartbeat("running"):
                    return await _call(tool, arguments, prefix, progress)
            finally:
                semaphore.release()
        except Exception as exc:  # the server's own fault, in words
            return _error(None, f"{type(exc).__name__}: {exc}")

    return Server(
        "rasputin",
        version=installed_version(),
        instructions=INSTRUCTIONS,
        on_list_tools=list_tools,
        on_call_tool=call_tool,
    )


async def serve(app: typer.Typer, jobs: int = 1) -> None:
    """``rasputin mcp``: the server on this process's stdin and stdout."""
    server = build_server(app, jobs=jobs)
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())
