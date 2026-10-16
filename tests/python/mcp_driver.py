"""One way to drive the rasputin MCP server in the tests (increment 37).

`docs/increments/37-mcp-server.md`, section 9. Test support only; nothing in
`src_python` uses it. The suites take these from here rather than from each
other, so no test module is a helper library for another.

The server under test is `tin_engine.mcp.server.build_server(app, *, jobs,
child)`, where ``child`` is the argv prefix the server runs a command with
(production: ``[sys.executable, "-m", "tin_engine"]``).

In process, the SDK's `Client(server)` connects without a transport; it
carries no progress token (section 2, measured), so progress is asserted over
stdio only, through `stdio_server`, which spawns `mcp_stdio_launcher`.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import pairwise
from pathlib import Path
from typing import Any

import anyio
from mcp.client import Client
from mcp.client.stdio import StdioServerParameters
from mcp_types import CallToolResult

HERE = Path(__file__).resolve().parent
FAKE_CHILD = HERE / "mcp_fake_child.py"
LAUNCHER = HERE / "mcp_stdio_launcher.py"

#: The tool names in the CLI's order (section 5), read from the app at test
#: time so the list and the CLI cannot drift; the server's own command is not
#: a tool.
SERVER_COMMAND = "mcp"

#: The sidecar names the server chooses beside ``out`` (section 6.4): a dot,
#: the mesh's file name, ``record-`` or ``stats-``, eight hex digits.
RECORD_SIDECAR = re.compile(r"^\.(?P<out>.+)\.record-[0-9a-f]{8}\.json$")
STATS_SIDECAR = re.compile(r"^\.(?P<out>.+)\.stats-[0-9a-f]{8}\.md$")

#: Every notification of a call carries elapsed seconds, never a count
#: (section 6.5); the heartbeat says so in words.
RUNNING = re.compile(r"^running, [0-9.]+ s$")
WAITING = re.compile(r"^waiting for a running job, [0-9.]+ s$")


def tool_names(app: Any) -> list[str]:
    """The expected tool list: the CLI's shown commands, in its order, ``-`` as ``_``."""
    import typer

    commands = typer.main.get_command(app).commands
    return [
        name.replace("-", "_")
        for name, command in commands.items()
        if name != SERVER_COMMAND and not command.hidden
    ]


def fake_child(tmp_path: Path, **spec: Any) -> list[str]:
    """A ``child`` prefix running `mcp_fake_child` with ``spec`` (its docstring lists the keys)."""
    path = tmp_path / f"fake-child-{len(list(tmp_path.glob('fake-child-*.json')))}.json"
    path.write_text(json.dumps(spec))
    return [sys.executable, str(FAKE_CHILD), str(path)]


def stdio_server(tmp_path: Path, **config: Any) -> StdioServerParameters:
    """The launcher as a stdio server (`mcp_stdio_launcher`: ``child``, ``jobs``, ``heartbeat``).

    The environment is passed whole: the SDK's default is a short allow-list,
    and the server's child must find the same interpreter and data as the test.
    """
    path = tmp_path / "launcher.json"
    path.write_text(json.dumps(config))
    return StdioServerParameters(
        command=sys.executable, args=[str(LAUNCHER), str(path)], env=dict(os.environ)
    )


def production_server() -> StdioServerParameters:
    """``python -m tin_engine mcp``: the entry an agent's configuration names (section 6.8)."""
    return StdioServerParameters(
        command=sys.executable, args=["-m", "tin_engine", "mcp"], env=dict(os.environ)
    )


async def call(server: Any, name: str, arguments: dict[str, Any], **options: Any) -> CallToolResult:
    """One tool call on a fresh connection; ``server`` is a `Server` or stdio parameters."""
    async with Client(server) as client:
        return await client.call_tool(name, arguments, **options)


def structured(result: CallToolResult) -> dict[str, Any]:
    assert isinstance(result.structured_content, dict), result
    return result.structured_content


def text_block(result: CallToolResult) -> str:
    """The one text block of a result."""
    texts = [b.text for b in result.content if b.type == "text"]
    assert len(texts) == 1, result.content
    return texts[0]


def resource_links(result: CallToolResult) -> list[tuple[str, str | None]]:
    """``(uri, mime_type)`` of each ``resource_link`` block, in order."""
    return [(b.uri, b.mime_type) for b in result.content if b.type == "resource_link"]


def dead(pid: int) -> bool:
    """Gone, or a zombie nobody has reaped yet (a kill without a wait)."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return True
    done = subprocess.run(["ps", "-o", "stat=", "-p", str(pid)], capture_output=True, text=True)
    state = done.stdout.strip()
    return state == "" or state.startswith("Z")


def pid_of(pid_file: Path) -> int:
    return int(pid_file.read_text())


async def wait_until_dead(pid_file: Path, within: float) -> None:
    """The pid in ``pid_file`` is dead within ``within`` seconds of now, polled at 50 ms."""
    pid = pid_of(pid_file)
    deadline = time.monotonic() + within
    while not dead(pid):
        assert time.monotonic() < deadline, f"child {pid} still alive {within} s after the cancel"
        await anyio.sleep(0.05)


@dataclass
class Progress:
    """Every progress notification of one call, with the client's own clock at receipt."""

    began: float = field(default_factory=time.monotonic)
    ticks: list[tuple[float, float | None, str | None, float]] = field(default_factory=list)

    async def __call__(self, progress: float, total: float | None, message: str | None) -> None:
        self.ticks.append((progress, total, message, time.monotonic() - self.began))

    @property
    def values(self) -> list[float]:
        return [t[0] for t in self.ticks]

    @property
    def messages(self) -> list[str | None]:
        return [t[2] for t in self.ticks]

    def assert_strictly_increasing(self) -> None:
        values = self.values
        assert all(b > a for a, b in pairwise(values)), values

    def assert_elapsed_seconds(self) -> None:
        """Each value is seconds since the handler began: never more than the client's own
        elapsed time at receipt (the handler began after the call was sent), and ``total``
        is absent (the run's length is unknown)."""
        for progress, total, _message, received in self.ticks:
            assert total is None, self.ticks
            assert 0.0 <= progress <= received, (progress, received, self.ticks)


def option_value(argv: list[str], flag: str) -> str:
    """The value after ``flag``, which appears exactly once in ``argv``."""
    assert argv.count(flag) == 1, (flag, argv)
    return argv[argv.index(flag) + 1]


def intervals(path: Path) -> list[tuple[float, float]]:
    """`mcp_fake_child`'s ``interval_file`` lines as ``(start, end)``, by start."""
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    return sorted((r["start"], r["end"]) for r in rows)


def overlap(a: tuple[float, float], b: tuple[float, float]) -> bool:
    return a[0] < b[1] and b[0] < a[1]


def only_files(directory: Path, names: Sequence[str]) -> None:
    """``directory`` holds exactly ``names`` (so no sidecar was left behind)."""
    assert sorted(p.name for p in directory.iterdir()) == sorted(names), list(directory.iterdir())
