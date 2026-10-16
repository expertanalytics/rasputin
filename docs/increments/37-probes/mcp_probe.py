"""Increment 37 probe: what the `mcp` SDK makes of the rasputin commands, and
what each execution route costs. Design evidence, not a test: run it by hand,
`python mcp_probe.py`, in a venv that has mcp>=2.3, imagecodecs and
`tin_engine` importable with a current `_core`. `--serve` makes this file the
stdio server the probe spawns, so the stdio round trip goes through the SDK's
own client.

The server here is the design in miniature: the SDK's low-level `Server`, one
tool per Typer command, each tool's input schema derived from the command's
click parameters (one source of truth) with `additionalProperties: false`,
arguments validated by the derived Pydantic model, the command run as a child
process with `TYPER_USE_RICH=0`, and the exit code, output lines, record and
plain-word refusal returned as data.

Parts (each prints lines starting with its letter):
  A  the derived schemas, listed back through the SDK's client; which of the
     parameters needed a note;
  B  error mapping: a refused command (exit 2), an argument of the wrong type,
     an unknown argument, an unknown tool;
  C  a tool awaiting a child process, cancelled by the client's timeout, in
     process and over stdio: does the child die?
  D  two concurrent calls;
  E  timings: interpreter start, `import tin_engine.cli`, `rasputin version`
     and the 1 m benchmark `mesh` as subprocesses; in-process `version`;
     in-process and stdio tool round trips; the mesh through the tool, with
     progress notifications counted.

The production server (increment 37) must not copy two shortcuts taken here:
the handler reads and unlinks temp files with blocking `Path` calls (ASYNC240),
and the measurement lines run past 100 columns so each stays one grep-able line.
"""
# ruff: noqa: E501, ASYNC240

from __future__ import annotations

import asyncio
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Literal

import anyio
import click
import mcp_types as types
import typer
from mcp.client import Client
from mcp.client.stdio import StdioServerParameters
from mcp.server.context import ServerRequestContext
from mcp.server.lowlevel import Server
from mcp.server.stdio import stdio_server
from mcp.shared.exceptions import MCPError
from pydantic import BaseModel, ConfigDict, Field, ValidationError, create_model

HERE = Path(__file__).resolve()
REPO = HERE.parents[3]
DEM = REPO / "tests/fixtures/dem_archive/7908_3_10m_z33.tif"
CLI = [sys.executable, "-c", "from tin_engine.cli import app; app()"]
CHILD_ENV = {**os.environ, "TYPER_USE_RICH": "0", "NO_COLOR": "1", "TERM": "dumb"}

SCALAR = {"str": str, "text": str, "int": int, "integer": int, "float": float, "bool": bool,
          "boolean": bool, "path": str}  # fmt: skip


def py_type(t: click.ParamType) -> tuple[Any, dict[str, Any], str | None]:
    """A click ParamType as an annotation, Field constraints, and a note or None."""
    if "int" in t.name and hasattr(t, "min") and hasattr(t, "max"):  # click's and Typer's IntRange
        return int, {"ge": t.min, "le": t.max}, None
    if isinstance(t, click.Choice):
        return Literal[tuple(t.choices)], {}, None
    inner = getattr(t, "types", None)
    if inner is not None:  # Typer's and click's tuple types both carry .types
        parts = [py_type(x) for x in inner]
        return tuple[tuple(p[0] for p in parts)], {}, next((p[2] for p in parts if p[2]), None)
    if t.name in SCALAR:
        return SCALAR[t.name], {}, None
    return str, {}, f"unmapped click type {t.name!r}"


class Derived:
    """A Typer command as a tool: its argument model and the way back to argv."""

    def __init__(self, cmd: click.Command, tool_name: str) -> None:
        self.cmd, self.tool_name, self.notes = cmd, tool_name, []
        fields: dict[str, Any] = {}
        for p in cmd.params:
            if p.name in ("help", "install_completion", "show_completion"):
                continue
            typ, constraints, note = py_type(p.type)
            if note:
                self.notes.append(f"{p.name}: {note}")
            if getattr(p, "callback", None) is not None:
                self.notes.append(f"{p.name}: callback {p.callback.__name__} stays in the child")
            if p.multiple:
                typ = list[typ]  # type: ignore[valid-type]
            help_text = getattr(p, "help", None) or ""
            if p.required:
                fields[p.name] = (typ, Field(description=help_text, **constraints))
            else:
                default = p.default() if callable(p.default) else p.default
                if default is None or (p.multiple and default == ()):
                    default, typ = None, typ | None  # type: ignore[operator]
                fields[p.name] = (typ, Field(default, description=help_text, **constraints))
        self.model: type[BaseModel] = create_model(
            f"{tool_name}_arguments", __config__=ConfigDict(extra="forbid"), **fields
        )

    def tool(self) -> types.Tool:
        schema = self.model.model_json_schema()
        schema.pop("title", None)
        return types.Tool(
            name=self.tool_name,
            description=(self.cmd.help or "").split("\n\n")[0],
            input_schema=schema,
        )

    def argv(self, args: BaseModel) -> list[str]:
        """The command line the arguments mean; only the keys the caller set."""
        argv = [self.cmd.name or self.tool_name]
        given = args.model_fields_set
        for p in self.cmd.params:
            if p.name not in given:
                continue
            value = getattr(args, p.name)
            if p.param_type_name == "argument":  # Typer 0.27's TyperArgument is no click.Argument
                argv.extend(str(v) for v in (value if p.nargs != 1 else [value]))
                continue
            opt = p.opts[0]
            if getattr(p, "is_flag", False):
                if value:
                    argv.append(opt)
                elif p.secondary_opts:
                    argv.append(p.secondary_opts[0])
                continue
            items = value if p.multiple else [value]
            for item in items:
                argv.append(opt)
                argv.extend(str(v) for v in (item if isinstance(item, tuple) else [item]))
        return argv


def derive_all() -> dict[str, Derived]:
    from tin_engine.cli import app

    root = typer.main.get_command(app)
    return {n.replace("-", "_"): Derived(c, n.replace("-", "_")) for n, c in root.commands.items()}


def error_result(exit_code: int | None, message: str) -> types.CallToolResult:
    return types.CallToolResult(
        content=[types.TextContent(type="text", text=message)],
        structured_content={"exit_code": exit_code, "error": message},
        is_error=True,
    )


def refusal_words(stderr: str) -> str:
    """Click's own `Error:` line(s), without the usage reminder."""
    lines = [ln for ln in stderr.splitlines() if ln and not ln.startswith(("Usage:", "Try '"))]
    return " ".join(lines).removeprefix("Error: ").strip()


CHILDREN: list[asyncio.subprocess.Process] = []


async def run_child(
    ctx: ServerRequestContext[Any], argv: list[str]
) -> tuple[int, list[str], list[str]]:
    """Run `rasputin argv` as a child; a stderr line is a progress tick; die with the caller."""
    proc = await asyncio.create_subprocess_exec(
        *CLI, *argv, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE, env=CHILD_ENV
    )
    CHILDREN.append(proc)
    token = (ctx.meta or {}).get("progress_token")
    err_lines: list[str] = []

    async def read_stderr() -> None:
        assert proc.stderr is not None
        async for raw in proc.stderr:
            err_lines.append(raw.decode(errors="replace").rstrip())
            if token is not None:
                await ctx.session.send_progress_notification(
                    token, len(err_lines), None, err_lines[-1][:200]
                )

    try:
        assert proc.stdout is not None
        async with asyncio.TaskGroup() as tg:
            tg.create_task(read_stderr())
            out_task = tg.create_task(proc.stdout.read())
        await proc.wait()
    finally:
        if proc.returncode is None:
            proc.kill()
            await proc.wait()
    return proc.returncode or 0, out_task.result().decode(errors="replace").splitlines(), err_lines


def build_server() -> tuple[Server[Any], dict[str, Derived]]:
    derived = derive_all()

    async def list_tools(ctx: Any, params: Any) -> types.ListToolsResult:
        slow = types.Tool(
            name="slow",
            description="Await a child `sleep`; kill it if cancelled.",
            input_schema={
                "type": "object",
                "properties": {"seconds": {"type": "number"}, "pid_file": {"type": "string"}},
                "required": ["seconds", "pid_file"],
                "additionalProperties": False,
            },
        )
        return types.ListToolsResult(tools=[d.tool() for d in derived.values()] + [slow])

    async def call_tool(
        ctx: ServerRequestContext[Any], params: types.CallToolRequestParams
    ) -> types.CallToolResult:
        args = params.arguments or {}
        if params.name == "slow":
            proc = await asyncio.create_subprocess_exec("sleep", str(args["seconds"]))
            Path(args["pid_file"]).write_text(str(proc.pid))
            try:
                await proc.wait()
            finally:
                if proc.returncode is None:
                    proc.kill()
                    await proc.wait()
            return types.CallToolResult(content=[], structured_content={"exit_code": proc.returncode})
        d = derived.get(params.name)
        if d is None:
            raise MCPError(code=types.INVALID_PARAMS, message=f"Unknown tool: {params.name}")
        try:
            model = d.model.model_validate(args)
        except ValidationError as exc:
            words = "; ".join(f"{'.'.join(map(str, e['loc']))}: {e['msg']}" for e in exc.errors())
            return error_result(None, words)
        argv = d.argv(model)
        record_tmp: Path | None = None
        if params.name == "mesh" and "record_path" not in model.model_fields_set:
            record_tmp = Path(tempfile.mkstemp(suffix=".json")[1])
            argv += ["--record", str(record_tmp)]
        t0 = time.perf_counter()
        code, out, err = await run_child(ctx, argv)
        result: dict[str, Any] = {
            "exit_code": code,
            "seconds": round(time.perf_counter() - t0, 3),
            "outputs": [ln for ln in out if ln and Path(ln).exists()] if params.name != "palette" else [],
            "stdout": out,
            "stderr": err,
        }
        if record_tmp is not None:
            if record_tmp.exists() and record_tmp.stat().st_size:
                result["record"] = json.loads(record_tmp.read_text())
                same = record_tmp.resolve()  # the CLI prints the resolved path (/private/var on macOS)
                result["stdout"] = [ln for ln in out if Path(ln).resolve() != same]
                result["outputs"] = [p for p in result["outputs"] if Path(p).resolve() != same]
            record_tmp.unlink(missing_ok=True)
        if code != 0:
            return error_result(code, f"exit {code}: {refusal_words(os.linesep.join(err))}")
        return types.CallToolResult(
            content=[types.TextContent(type="text", text=json.dumps(result)[:4000])],
            structured_content=result,
        )

    server: Server[Any] = Server("rasputin-probe", on_list_tools=list_tools, on_call_tool=call_tool)
    return server, derived


async def amedian(fn: Any, n: int) -> tuple[float, float]:
    walls = []
    for _ in range(n):
        t0 = time.perf_counter()
        await fn()
        walls.append(time.perf_counter() - t0)
    return statistics.median(walls), max(walls)


def median_wall(fn: Any, n: int) -> tuple[float, float]:
    walls = []
    for _ in range(n):
        t0 = time.perf_counter()
        fn()
        walls.append(time.perf_counter() - t0)
    return statistics.median(walls), max(walls)


def alive(pid_file: Path) -> str:
    pid = int(pid_file.read_text())
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return f"pid {pid} dead"
    done = subprocess.run(["ps", "-o", "stat=", "-p", str(pid)], capture_output=True, text=True)
    return f"pid {pid} ALIVE (ps stat {done.stdout.strip()!r})" if done.stdout.strip() else f"pid {pid} dead"


async def cancel_check(client: Client, label: str, tmp: Path) -> None:
    pid_file = tmp / f"{label}.pid"
    t0 = time.perf_counter()
    try:
        await client.call_tool("slow", {"seconds": 30, "pid_file": str(pid_file)}, read_timeout_seconds=1.0)
        print(f"C {label}: timeout did not raise")
    except Exception as exc:
        print(f"C {label}: client raised {type(exc).__name__} after {time.perf_counter() - t0:.2f} s")
    await anyio.sleep(0.5)
    print(f"C {label}: child sleep {alive(pid_file)}")


async def probe(tmp: Path) -> None:
    server, derived = build_server()
    async with Client(server) as client:
        listed = await client.list_tools()
        tools = {t.name: t for t in listed.tools}
        print(f"A tools listed: {sorted(tools)}; protocol {client.protocol_version}")
        print(f"A palette input_schema:\n{json.dumps(tools['palette'].input_schema, indent=1)}")
        mesh = tools["mesh"].input_schema
        print(
            f"A mesh: {len(mesh['properties'])} properties, required {mesh.get('required')}, "
            f"additionalProperties {mesh.get('additionalProperties')}, {len(json.dumps(mesh))} bytes"
        )
        for key in ("bbox", "tolerance_near", "binary", "features", "snap_spacing"):
            print(f"A mesh.{key}: {json.dumps(mesh['properties'][key])}")
        print(f"A fetch.connections: {json.dumps(tools['fetch'].input_schema['properties']['connections'])}")
        print(f"A all schemas: {sum(len(json.dumps(t.input_schema)) for t in tools.values())} bytes")
        for name, d in derived.items():
            for note in d.notes:
                print(f"A note {name}: {note}")

        r = await client.call_tool("mesh", {"dem": [str(DEM)], "bbox": [1, 2, 3, 4]})
        print(f"B refused (no out): is_error={r.is_error} structured={r.structured_content}")
        r = await client.call_tool("mesh", {"out": str(tmp / "x.vtk"), "dem": ["/nonexistent.tif"]})
        print(f"B refused (bad dem): is_error={r.is_error} structured={r.structured_content}")
        r = await client.call_tool("palette", {"name": 3})
        print(f"B wrong type: is_error={r.is_error} text={r.content[0].text!r}")
        r = await client.call_tool("palette", {"name": "corine", "bogus": 1})
        print(f"B unknown argument: is_error={r.is_error} text={r.content[0].text!r}")
        r = await client.call_tool("mesh", {"out": str(tmp / "x.vtk"), "tolerance": "ten"})
        print(f"B wrong type nested: is_error={r.is_error} text={r.content[0].text!r}")
        try:
            r = await client.call_tool("nosuch", {})
            print(f"B unknown tool: is_error={r.is_error} text={r.content[0].text!r}")
        except Exception as exc:
            print(f"B unknown tool: client raised {type(exc).__name__}: {exc}")
        r = await client.call_tool("palette", {"name": "corine"})
        sc = r.structured_content or {}
        print(
            f"B palette corine: is_error={r.is_error}, exit {sc.get('exit_code')}, "
            f"{len(sc.get('stdout', []))} stdout lines, {sc.get('seconds')} s, error {sc.get('error')!r}"
        )

        await cancel_check(client, "in-process", tmp)

        t0 = time.perf_counter()
        a, b = await asyncio.gather(client.call_tool("version", {}), client.call_tool("version", {}))
        print(
            f"D two concurrent version calls: {time.perf_counter() - t0:.2f} s, exits "
            f"{(a.structured_content or {}).get('exit_code')}, {(b.structured_content or {}).get('exit_code')}"
        )

        med, mx = await amedian(lambda: client.call_tool("version", {}), 5)
        print(f"E in-process round trip, version tool: median {med:.3f} s max {mx:.3f}")

        ticks: list[str] = []

        async def on_progress(progress: float, total: float | None, message: str | None) -> None:
            ticks.append(message or "")

        out = tmp / "tool.vtk"
        r = await client.call_tool(
            "mesh", {"dem": [str(DEM)], "tolerance": 1.0, "out": str(out)}, progress_callback=on_progress
        )
        sc = r.structured_content or {}
        print(
            f"E mesh through the tool: is_error={r.is_error}, child {sc.get('seconds')} s, "
            f"outputs {sc.get('outputs')}, record keys {len(sc.get('record', {}))}, "
            f"stderr lines {len(sc.get('stderr', []))}, progress ticks {len(ticks)}, "
            f"text block {len(r.content[0].text)} chars"
        )
        print(f"E mesh stderr lines: {sc.get('stderr')}")
        print(f"E mesh record sample: {json.dumps(dict(list(sc.get('record', {}).items())[:6]))}")

    t0 = time.perf_counter()
    params = StdioServerParameters(command=sys.executable, args=[str(HERE), "--serve"])
    async with Client(params) as client:
        connect = time.perf_counter() - t0
        listed = await client.list_tools()
        med, mx = await amedian(lambda: client.call_tool("version", {}), 5)
        print(
            f"E stdio: connect+initialize {connect:.3f} s, {len(listed.tools)} tools, "
            f"version round trip median {med:.3f} s max {mx:.3f}"
        )
        await cancel_check(client, "stdio", tmp)
        ticks = []
        r = await client.call_tool(
            "mesh",
            {"dem": [str(DEM)], "tolerance": 1.0, "out": str(tmp / "stdio.vtk")},
            progress_callback=on_progress,
        )
        sc = r.structured_content or {}
        print(
            f"E stdio mesh through the tool: is_error={r.is_error}, child {sc.get('seconds')} s, "
            f"progress ticks {len(ticks)} {ticks[:2]}, outputs {len(sc.get('outputs', []))}, "
            f"record keys {len(sc.get('record', {}))}"
        )


def subprocess_timings(tmp: Path) -> None:
    def sub(args: list[str]) -> Any:
        return lambda: subprocess.run(args, capture_output=True, check=True)

    med, mx = median_wall(sub([sys.executable, "-c", "pass"]), 5)
    print(f"E python -c pass: median {med:.3f} s max {mx:.3f}")
    med, mx = median_wall(sub([sys.executable, "-c", "import tin_engine.cli"]), 5)
    print(f"E import tin_engine.cli: median {med:.3f} s max {mx:.3f}")
    med, mx = median_wall(sub([sys.executable, "-c", "import mcp.server.lowlevel, mcp.client"]), 5)
    print(f"E import mcp (server+client): median {med:.3f} s max {mx:.3f}")
    med, mx = median_wall(sub([*CLI, "version"]), 5)
    print(f"E rasputin version (subprocess): median {med:.3f} s max {mx:.3f}")

    from typer.testing import CliRunner

    from tin_engine.cli import app

    runner = CliRunner()
    med, mx = median_wall(lambda: runner.invoke(app, ["version"]), 5)
    print(f"E rasputin version (in-process CliRunner): median {med:.4f} s max {mx:.4f}")

    out, rec, walls = tmp / "m.vtk", tmp / "r.json", []
    for _ in range(3):
        out.unlink(missing_ok=True)
        t0 = time.perf_counter()
        done = subprocess.run(
            [*CLI, "mesh", "--dem", str(DEM), "--tolerance", "1", "--out", str(out), "--record", str(rec)],
            capture_output=True, text=True, env=CHILD_ENV,
        )  # fmt: skip
        walls.append(time.perf_counter() - t0)
        if done.returncode != 0:
            print(f"E mesh FAILED: {done.stderr[-400:]}")
            return
    record = json.loads(rec.read_text())
    print(
        f"E rasputin mesh 7908_3 --tolerance 1 (subprocess): median {statistics.median(walls):.2f} s "
        f"of {[round(w, 2) for w in walls]}; stdout {done.stdout.splitlines()}; "
        f"stderr {len(done.stderr.splitlines())} lines; record {rec.stat().st_size} bytes, "
        f"{len(record)} keys; mesh {out.stat().st_size / 1e6:.1f} MB"
    )


def main() -> None:
    if "--serve" in sys.argv:
        server, _ = build_server()

        async def serve() -> None:
            async with stdio_server() as (read, write):
                await server.run(read, write, server.create_initialization_options())

        anyio.run(serve)
        return
    import importlib.metadata as m

    print(f"probe {HERE}\npython {sys.version.split()[0]}, mcp {m.version('mcp')}, typer "
          f"{m.version('typer')}, pydantic {m.version('pydantic')}")  # fmt: skip
    print(subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout.splitlines()[0])
    with tempfile.TemporaryDirectory() as tmp:
        anyio.run(probe, Path(tmp))
        subprocess_timings(Path(tmp))


if __name__ == "__main__":
    main()
