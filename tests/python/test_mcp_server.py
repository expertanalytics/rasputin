"""Increment 37, T4 to T6, T8 and T10 in process: the server runs the CLI as a child.

`docs/increments/37-mcp-server.md`, sections 5, 6.4 to 6.6, 7 and 9. The
SDK's `Client(server)` connects in process; the child is `mcp_fake_child`
where the behaviour under test is the server's (validation, refusal mapping,
sidecar removal, cancellation, the jobs cap) and the real ``python -m
tin_engine`` where it is the whole round trip (T6, the micro tile and the
Kartverket tile).

Clocks: every bound here is on the 0.5 s fake child and the 0.72 s probe run
(section 3), with the slack each comment states; nothing waits on the wall
clock without a bound.
"""

from __future__ import annotations

import asyncio
import json
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("mcp")

from mcp.client import Client
from mcp.shared.exceptions import MCPError

from geotiff_fixtures import KARTVERKET, micro_tiff, needs_codecs
from mcp_driver import (
    RECORD_SIDECAR,
    STATS_SIDECAR,
    Progress,
    call,
    fake_child,
    intervals,
    only_files,
    option_value,
    overlap,
    resource_links,
    structured,
    text_block,
    wait_until_dead,
)
from tin_engine.cli import app
from tin_engine.mcp import server as mcp_server
from tin_engine.run_record import WORDING
from vtkread import read_vtk

Build = Callable[..., Any]


@pytest.fixture
def build(tmp_path: Path) -> Build:
    """`build_server` on the real app; ``child`` and ``jobs`` as the test says, and the
    output root the test's own directory unless it says otherwise, so a write the test
    names is inside it."""

    def _build(**options: Any) -> Any:
        options.setdefault("out_root", tmp_path)
        return mcp_server.build_server(app, **options)

    return _build


@pytest.fixture
def out_dir(tmp_path: Path) -> Path:
    """Where ``out`` goes: its own directory, so a leftover sidecar is visible."""
    directory = tmp_path / "out"
    directory.mkdir()
    return directory


@pytest.fixture
def micro_dem(tmp_path: Path) -> Path:
    """The 3 x 4 node tile at 10 m x 5 m (`geotiff_fixtures`), no codec needed."""
    path = tmp_path / "micro.tif"
    path.write_bytes(micro_tiff().getvalue())
    return path


def mesh_arguments(dem: Path, out: Path | str, **more: Any) -> dict[str, Any]:
    return {"dem": [str(dem)], "out": str(out), **more}


# ---------------------------------------------------------------- T4


class TestValidationBeforeTheChild:
    async def test_a_missing_out_is_the_models_refusal_verbatim(
        self, tmp_path: Path, build: Build
    ) -> None:
        argv_file = tmp_path / "argv.json"
        server = build(child=fake_child(tmp_path, argv_file=str(argv_file)))
        result = await call(server, "mesh", {"dem": ["d.tif"]})
        assert result.is_error
        assert structured(result) == {"exit_code": None, "error": "out: Field required"}
        assert "out: Field required" in text_block(result)
        assert not argv_file.exists(), "the child ran on arguments the model should have refused"

    @pytest.mark.parametrize(
        ("arguments", "error_prefix"),
        [
            pytest.param(
                {"dem": ["d.tif"], "out": "a.vtk", "tolerance": "ten"},
                "tolerance: Input should be a valid number",
                id="wrong type",
            ),
            pytest.param(
                {"dem": ["d.tif"], "out": "a.vtk", "bogus": 1},
                "bogus: Extra inputs are not permitted",
                id="unknown key",
            ),
            pytest.param(
                {"dem": ["d.tif"], "out": "a.vtk", "bbox": [1.0, 2.0, 3.0]},
                "bbox: ",
                id="bbox of three",
            ),
        ],
    )
    async def test_refused_naming_the_field_and_the_child_never_runs(
        self, tmp_path: Path, build: Build, arguments: dict[str, Any], error_prefix: str
    ) -> None:
        argv_file = tmp_path / "argv.json"
        server = build(child=fake_child(tmp_path, argv_file=str(argv_file)))
        result = await call(server, "mesh", arguments)
        assert result.is_error
        body = structured(result)
        assert body["exit_code"] is None
        assert body["error"].startswith(error_prefix), body
        assert not argv_file.exists()

    async def test_stats_dash_becomes_a_stats_sidecar_beside_out(
        self, tmp_path: Path, out_dir: Path, build: Build
    ) -> None:
        argv_file = tmp_path / "argv.json"
        server = build(child=fake_child(tmp_path, argv_file=str(argv_file)))
        await call(server, "mesh", mesh_arguments(Path("d.tif"), out_dir / "a.vtk", stats="-"))
        argv = json.loads(argv_file.read_text())
        stats = Path(option_value(argv, "--stats"))
        assert stats.parent == out_dir
        match = STATS_SIDECAR.fullmatch(stats.name)
        assert match is not None and match["out"] == "a.vtk", stats.name
        assert RECORD_SIDECAR.fullmatch(Path(option_value(argv, "--record")).name)

    async def test_a_given_record_path_is_the_callers_and_no_sidecar_is_added(
        self, tmp_path: Path, out_dir: Path, build: Build
    ) -> None:
        argv_file = tmp_path / "argv.json"
        server = build(child=fake_child(tmp_path, argv_file=str(argv_file)))
        keep = out_dir / "keep.json"
        await call(
            server, "mesh", mesh_arguments(Path("d.tif"), out_dir / "a.vtk", record_path=str(keep))
        )
        argv = json.loads(argv_file.read_text())
        assert argv.count("--record") == 1
        assert option_value(argv, "--record") == str(keep)


# ---------------------------------------------------------------- T5

REFUSAL = [
    "Usage: rasputin mesh [OPTIONS]",
    "Try 'rasputin mesh --help' for help.",
    "",
    "Error: Invalid value for --dem: cannot read /nonexistent.tif: No such file or directory",
]


class TestRefusalMapping:
    async def test_exit_2_gives_the_error_lines_words_and_the_whole_stderr(
        self, tmp_path: Path, out_dir: Path, build: Build
    ) -> None:
        server = build(child=fake_child(tmp_path, stderr=[[0, line] for line in REFUSAL], exit=2))
        result = await call(
            server, "mesh", {"out": str(out_dir / "a.vtk"), "dem": ["/nonexistent.tif"]}
        )
        assert result.is_error
        body = structured(result)
        assert body["exit_code"] == 2
        assert body["error"] == f"exit 2: {REFUSAL[-1].removeprefix('Error: ')}"
        assert body["stderr"] == REFUSAL
        assert "Invalid value for --dem" in text_block(result)

    async def test_exit_1_likewise(self, tmp_path: Path, build: Build) -> None:
        words = "Error: the catalogue is unreachable"
        server = build(child=fake_child(tmp_path, stderr=[[0, words]], exit=1))
        result = await call(server, "version", {})
        assert result.is_error
        body = structured(result)
        assert (body["exit_code"], body["error"]) == (1, "exit 1: the catalogue is unreachable")

    async def test_progress_lines_before_the_error_line_stay_out_of_error(
        self, tmp_path: Path, build: Build
    ) -> None:
        """A `fetch` failing mid-run prints its progress first (`cli.py`, the ``bytes`` lines,
        then the ``Error:`` line); `error` carries the ``Error:`` line's words, `stderr` all."""
        lines = ["1,000 of 2,000 bytes", "2,000 of 2,000 bytes", "Error: no route to the catalogue"]
        server = build(child=fake_child(tmp_path, stderr=[[0, line] for line in lines], exit=1))
        result = await call(server, "fetch", {"source": "glo30"})
        assert result.is_error
        body = structured(result)
        assert (body["exit_code"], body["error"]) == (1, "exit 1: no route to the catalogue")
        assert body["stderr"] == lines

    async def test_without_an_error_line_every_line_but_the_usage_reminder_is_kept(
        self, tmp_path: Path, build: Build
    ) -> None:
        """The fallback: a child that dies without an ``Error:`` line (a traceback, an
        ``Aborted!``) still has its words in `error`, less ``Usage:`` and ``Try ...``."""
        lines = [
            "Usage: rasputin fetch [OPTIONS] {source}",
            "Try 'rasputin fetch --help' for help.",
            "",
            "1,000 of 2,000 bytes",
            "Aborted!",
        ]
        server = build(child=fake_child(tmp_path, stderr=[[0, line] for line in lines], exit=1))
        result = await call(server, "fetch", {"source": "glo30"})
        assert result.is_error
        body = structured(result)
        assert (body["exit_code"], body["error"]) == (1, "exit 1: 1,000 of 2,000 bytes Aborted!")
        assert body["stderr"] == lines

    async def test_a_signal_is_named(self, tmp_path: Path, build: Build) -> None:
        server = build(child=fake_child(tmp_path, signal=15))
        result = await call(server, "version", {})
        assert result.is_error
        body = structured(result)
        assert (body["exit_code"], body["error"]) == (-15, "killed by signal 15 (SIGTERM)")

    @pytest.mark.parametrize("create_record", [False, True], ids=["no sidecar", "sidecar created"])
    async def test_no_record_sidecar_remains_after_a_refusal(
        self, tmp_path: Path, out_dir: Path, build: Build, create_record: bool
    ) -> None:
        server = build(
            child=fake_child(
                tmp_path, stderr=[[0, REFUSAL[-1]]], exit=2, create_record=create_record
            )
        )
        result = await call(server, "mesh", {"out": str(out_dir / "a.vtk")})
        assert result.is_error
        only_files(out_dir, [])


# ---------------------------------------------------------------- T6


def assert_micro_mesh(result: Any, out: Path) -> dict[str, Any]:
    """A successful `mesh` of the micro tile to ``out``: the result, as section 5 shapes it."""
    assert not result.is_error, result.content
    body = structured(result)
    assert body["exit_code"] == 0
    assert body["command"].startswith("rasputin mesh "), body["command"]
    assert body["seconds"] > 0
    assert [Path(p).resolve() for p in body["outputs"]] == [out.resolve()]
    record = body["record"]
    assert isinstance(record, dict)
    assert list(record)[:2] == ["rasputin_version", "command"]
    assert list(record)[2:] == [k for k in WORDING if k in record], "the record's own order"
    assert body["stderr"][0].startswith("12 triangles"), body["stderr"]
    vtk = read_vtk(out.read_bytes())
    assert len(vtk.polygons) == 12  # 3 x 4 nodes: 2 x 3 cells, two triangles each
    assert json.loads(text_block(result)) == body
    assert resource_links(result) == [(out.resolve().as_uri(), "application/octet-stream")]
    return body


class TestRealMesh:
    async def test_mesh_writes_the_vtk_returns_the_record_and_leaves_no_sidecar(
        self, out_dir: Path, micro_dem: Path, build: Build
    ) -> None:
        out = out_dir / "a.vtk"
        result = await call(build(), "mesh", mesh_arguments(micro_dem, out))
        body = assert_micro_mesh(result, out)
        assert body["stats"] is None
        only_files(out_dir, ["a.vtk"])

    async def test_a_relative_out_lands_in_the_servers_working_directory(
        self, out_dir: Path, micro_dem: Path, build: Build, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.chdir(out_dir)
        result = await call(build(), "mesh", mesh_arguments(micro_dem, "rel.vtk"))
        assert_micro_mesh(result, out_dir / "rel.vtk")
        only_files(out_dir, ["rel.vtk"])

    async def test_stats_dash_returns_the_report_as_data_not_on_stdout(
        self, out_dir: Path, micro_dem: Path, build: Build
    ) -> None:
        out = out_dir / "a.vtk"
        result = await call(build(), "mesh", mesh_arguments(micro_dem, out, stats="-"))
        body = assert_micro_mesh(result, out)
        assert body["stats"].startswith("# rasputin mesh — statistics"), body["stats"][:80]
        assert "## Sizes" in body["stats"]
        markdown = [line for line in body["stdout"] if line.startswith(("#", "|", "`"))]
        assert markdown == [], body["stdout"]
        only_files(out_dir, ["a.vtk"])

    @needs_codecs
    async def test_the_kartverket_tile_at_one_metre(
        self, out_dir: Path, build: Build, capsys: pytest.CaptureFixture[str]
    ) -> None:
        out = out_dir / "k.vtk"
        began = time.perf_counter()
        result = await call(build(), "mesh", mesh_arguments(KARTVERKET, out, tolerance=1.0))
        seconds = time.perf_counter() - began
        assert not result.is_error, result.content
        body = structured(result)
        record = body["record"]
        assert record["tolerance_m"] == 1
        assert record["max_error_m"] <= 1
        triangles = int(body["stderr"][0].split()[0])
        # 464,290 measured on the probe's machine (section 3); the existing real-tile
        # test pins no exact count, so a 12 % window either side does for this one.
        assert 410_000 <= triangles <= 520_000, body["stderr"]
        only_files(out_dir, ["k.vtk"])
        with capsys.disabled():
            print(f"\nT6 Kartverket through the tool: {triangles} triangles, {seconds:.2f} s")


# ---------------------------------------------------------------- T8, T9 in process


class TestCancellation:
    async def test_a_client_timeout_kills_the_child(self, tmp_path: Path, build: Build) -> None:
        pid_file = tmp_path / "child.pid"
        server = build(child=fake_child(tmp_path, pid_file=str(pid_file), sleep=30))
        async with Client(server) as client:
            with pytest.raises(MCPError):
                await client.call_tool("version", {}, read_timeout_seconds=0.5)
            # measured: dead 0.5 s after a 1.0 s timeout (section 3, probe C)
            await wait_until_dead(pid_file, within=1.0)


async def test_in_process_a_call_with_a_progress_callback_completes(
    tmp_path: Path, build: Build
) -> None:
    # The in-process transport carries no progress token (section 2, measured), so
    # the notifications are asserted over stdio (test_mcp_stdio); here only that the
    # run completes with a callback attached.
    server = build(child=fake_child(tmp_path, stderr=[[0, "a line"]]))
    result = await call(server, "version", {}, progress_callback=Progress())
    assert structured(result)["exit_code"] == 0


# ---------------------------------------------------------------- T10


class TestJobsCap:
    async def test_jobs_1_runs_two_calls_one_after_the_other(
        self, tmp_path: Path, build: Build
    ) -> None:
        log = tmp_path / "intervals.jsonl"
        server = build(jobs=1, child=fake_child(tmp_path, sleep=0.5, interval_file=str(log)))
        async with Client(server) as client:
            began = time.monotonic()
            first, second = await asyncio.gather(
                client.call_tool("version", {}), client.call_tool("version", {})
            )
            wall = time.monotonic() - began
        assert not first.is_error and not second.is_error
        runs = intervals(log)
        assert len(runs) == 2 and not overlap(*runs), runs
        assert wall >= 1.0, wall  # two 0.5 s sleeps in sequence; start-up only adds

    async def test_jobs_2_runs_them_side_by_side(self, tmp_path: Path, build: Build) -> None:
        log = tmp_path / "intervals.jsonl"
        server = build(jobs=2, child=fake_child(tmp_path, sleep=0.5, interval_file=str(log)))
        async with Client(server) as client:
            began = time.monotonic()
            await asyncio.gather(client.call_tool("version", {}), client.call_tool("version", {}))
            wall = time.monotonic() - began
        runs = intervals(log)
        assert len(runs) == 2 and overlap(*runs), runs
        # two 0.5 s sleeps side by side: about 0.55 s with the interpreter's start-up;
        # the bound leaves 0.45 s for a loaded machine and still refutes a sequence (1.0 s)
        assert wall < 1.0, wall
