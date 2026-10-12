"""Increment 37, T7 to T9 over stdio: the transport an agent's configuration uses.

`docs/increments/37-mcp-server.md`, sections 6.5, 6.6, 6.8 and 9. The SDK's
client spawns the server as a process and speaks JSON-RPC on its pipes, so
this is where a stray write to the server's stdout would show, and the only
transport that carries a progress token (section 2, measured). T7 spawns the
production entry, ``python -m tin_engine mcp``; T8 and T9 spawn
`mcp_stdio_launcher` with `mcp_fake_child` as the child.

Clocks: a progress value is seconds since the handler began. It is bounded on
both sides without a tolerance: at least the sleeps the child has done before
the line (0.2 s each), and at most the client's own elapsed time at receipt,
since the handler began after the call was sent. A line count would fail the
upper bound at the first line.
"""

from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

pytest.importorskip("mcp")

import anyio
from mcp.client import Client
from mcp.shared.exceptions import MCPError

from mcp_driver import (
    RUNNING,
    WAITING,
    Progress,
    fake_child,
    production_server,
    stdio_server,
    structured,
    tool_names,
    wait_until_dead,
)
from tin_engine.cli import app

# ---------------------------------------------------------------- T7


async def test_python_m_tin_engine_mcp_lists_the_tools_and_runs_version() -> None:
    # connect and initialize measured at 0.465 s, a version call at 0.264 s (section 3)
    with anyio.fail_after(30):
        async with Client(production_server()) as client:
            listed = await client.list_tools()
            result = await client.call_tool("version", {})
    assert [tool.name for tool in listed.tools] == tool_names(app)
    assert not result.is_error, result.content
    body = structured(result)
    assert body["exit_code"] == 0
    assert len(body["stdout"]) == 2, body["stdout"]
    assert body["stdout"][1].startswith("bounds checks: ")
    assert body["outputs"] == []


# ---------------------------------------------------------------- T8


async def test_a_client_timeout_kills_the_child_over_stdio(tmp_path: Path) -> None:
    pid_file = tmp_path / "child.pid"
    params = stdio_server(tmp_path, child=fake_child(tmp_path, pid_file=str(pid_file), sleep=30))
    async with Client(params) as client:
        with pytest.raises(MCPError):
            await client.call_tool("version", {}, read_timeout_seconds=0.5)
        # measured: dead 0.5 s after a 1.0 s timeout (section 3, probe C)
        await wait_until_dead(pid_file, within=1.0)


# ---------------------------------------------------------------- T9


class TestProgress:
    async def test_each_stderr_line_is_a_notification_carrying_the_line(
        self, tmp_path: Path
    ) -> None:
        lines = ["one", "two", "three"]
        params = stdio_server(
            tmp_path, child=fake_child(tmp_path, stderr=[[0.2, s] for s in lines])
        )
        async with Client(params) as client:
            progress = Progress()
            result = await client.call_tool("version", {}, progress_callback=progress)
        assert structured(result)["stderr"] == lines
        assert progress.messages == lines
        progress.assert_strictly_increasing()
        progress.assert_elapsed_seconds()
        for i, value in enumerate(progress.values, start=1):
            assert value >= 0.2 * i, (i, value)  # i sleeps of 0.2 s came before line i

    async def test_heartbeats_say_running_with_elapsed_seconds(self, tmp_path: Path) -> None:
        params = stdio_server(
            tmp_path,
            heartbeat=0.1,
            child=fake_child(tmp_path, stderr=[[0.25, "phase"]], sleep=0.35),
        )
        async with Client(params) as client:
            progress = Progress()
            result = await client.call_tool("version", {}, progress_callback=progress)
        assert structured(result)["exit_code"] == 0
        messages = [m or "" for m in progress.messages]
        assert all(m == "phase" or RUNNING.fullmatch(m) for m in messages), messages
        assert messages.count("phase") == 1
        line = messages.index("phase")
        # heartbeats every 0.1 s over a 0.6 s run: five expected, at least three asserted,
        # one before the line at 0.25 s and one after (a loaded machine delays, never adds)
        assert len(messages) - 1 >= 3, messages
        assert 0 < line < len(messages) - 1, messages
        progress.assert_strictly_increasing()
        progress.assert_elapsed_seconds()

    async def test_a_queued_call_reports_waiting_then_running_on_one_increasing_clock(
        self, tmp_path: Path
    ) -> None:
        params = stdio_server(
            tmp_path, jobs=1, heartbeat=0.1, child=fake_child(tmp_path, sleep=0.5)
        )
        first, second = Progress(), Progress()
        async with Client(params) as client:
            results = await asyncio.gather(
                client.call_tool("version", {}, progress_callback=first),
                client.call_tool("version", {}, progress_callback=second),
            )
        assert all(structured(r)["exit_code"] == 0 for r in results)
        waited = [p for p in (first, second) if any(WAITING.fullmatch(m or "") for m in p.messages)]
        assert len(waited) == 1, "of two concurrent calls, exactly one queues behind the other"
        (queued,) = waited
        kinds = ["waiting" if WAITING.fullmatch(m or "") else "running" for m in queued.messages]
        assert all(RUNNING.fullmatch(m or "") for m, k in zip(queued.messages, kinds, strict=True)
                   if k == "running"), queued.messages  # fmt: skip
        # the first call holds the semaphore for 0.5 s plus start-up: at least three
        # waiting heartbeats at 0.1 s, then the running ones, never interleaved
        assert kinds.count("waiting") >= 3 and "running" in kinds, kinds
        assert kinds == sorted(kinds, key=("waiting", "running").index), kinds
        queued.assert_strictly_increasing()
        queued.assert_elapsed_seconds()
        for progress in (first, second):
            progress.assert_strictly_increasing()
