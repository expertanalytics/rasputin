"""The `rasputin.commands` plugin hook: discovery, merging, refusals, `version`, MCP.

`docs/increments/51-open-core.md`, *4.1 The plugin hook (Python)* and *9.1
Public (P1 unless marked)*, T1 to T6. Every plugin here is a real
distribution on the path (`plugin_fixtures`), found by ``importlib.metadata``
itself. In process, the directory goes on ``sys.path``; in a child, on
``PYTHONPATH``. The module under test is imported in a fixture, so each test
fails on its own while it is missing.

Clocks: the only bounds are hang guards on child processes and the MCP
connection (``python -m tin_engine --help`` measured about 0.3 s; a server
start and a ``version`` call under 1 s in `test_mcp_stdio`), each 30 s or more.
"""

from __future__ import annotations

import re
import subprocess
import sys
import tomllib
from collections.abc import Callable, Iterator
from importlib import import_module
from importlib.metadata import entry_points
from pathlib import Path
from types import ModuleType
from typing import Any

import anyio
import pytest
import typer
from typer.testing import CliRunner

from plugin_fixtures import (
    IMPORT_ERROR,
    MESH,
    NOT_TYPER,
    RUNTIME_ERROR,
    SCRIBBLE,
    child_env,
    hello_from,
    lines_naming,
    plain,
    run_rasputin,
    write_plugin,
)
from tin_engine.cli import app as core_app

ROOT = Path(__file__).resolve().parents[2]
PluginPath = Callable[..., Path]


def hook() -> ModuleType:
    """The module under test, imported in each test so each fails on its own."""
    return import_module("tin_engine.plugins")


@pytest.fixture
def on_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Iterator[PluginPath]:
    """Write a plugin (`write_plugin`'s keywords) into a directory of its own
    and put that directory first on ``sys.path``; the plugin modules are
    dropped from ``sys.modules`` afterwards, since variants share a name."""
    modules: list[str] = []

    def _install(**spec: Any) -> Path:
        root = write_plugin(tmp_path / f"site-{len(modules)}", **spec)
        modules.append(spec.get("module", "fakeplug"))
        monkeypatch.syspath_prepend(str(root))
        return root

    yield _install
    for name in modules:
        sys.modules.pop(name, None)


def commands(app: typer.Typer) -> dict[str, Any]:
    """The top-level commands of ``app``, by name, as Click sees them."""
    return dict(typer.main.get_command(app).commands)  # type: ignore[attr-defined]


def core_names() -> list[str]:
    return list(commands(core_app))


def shown_core_names() -> list[str]:
    """The core's commands that `rasputin --help` lists: the hidden ones left out."""
    return [name for name, command in commands(core_app).items() if not command.hidden]


def lists_command(help_text: str, name: str) -> bool:
    """Whether ``--help`` lists ``name`` as a command (a whole word on a line)."""
    return re.search(rf"(?m)^\W*{re.escape(name)}\b", plain(help_text)) is not None


# ---------------------------------------------------------------- the fixture


def test_the_fixture_is_found_and_loaded_by_importlib_metadata(on_path: PluginPath) -> None:
    """The control for every test below: the plugin is a real install as far
    as ``importlib.metadata`` can tell, and its entry point loads a Typer."""
    on_path()
    (found,) = [
        ep
        for ep in entry_points(group="rasputin.commands")
        if ep.dist and ep.dist.name == "fakeplug"
    ]
    assert (found.name, found.value, found.dist.version) == ("hello", "fakeplug:app", "0.1")
    loaded = found.load()
    assert isinstance(loaded, typer.Typer)
    assert CliRunner().invoke(loaded, []).output == "hello\n"


# ---------------------------------------------------------------- T1


def test_application_lists_the_plugin_command_and_runs_it(on_path: PluginPath) -> None:
    plugins = hook()
    on_path()
    app = plugins.application()
    assert "hello" in commands(app)
    assert set(core_names()) <= set(commands(app))
    result = CliRunner().invoke(app, ["hello"])
    assert result.exit_code == 0, result.output
    assert result.output == "hello\n"


def test_the_core_app_is_never_changed(on_path: PluginPath) -> None:
    plugins = hook()
    before = core_names()
    on_path()
    plugins.application()
    assert core_names() == before
    assert "hello" not in commands(core_app)


def test_the_same_plugins_give_the_same_app(on_path: PluginPath) -> None:
    plugins = hook()
    on_path()
    found = plugins.discover()
    assert isinstance(found, tuple)
    first, second = plugins.application(found), plugins.application(found)
    assert first is not second
    assert list(commands(first)) == list(commands(second))
    assert "hello" in commands(first)


def test_discover_finds_one_record_per_entry_point(on_path: PluginPath) -> None:
    plugins = hook()
    baseline = len(plugins.discover())
    on_path()
    assert len(plugins.discover()) == baseline + 1


def test_python_m_tin_engine_lists_and_runs_the_plugin_command(tmp_path: Path) -> None:
    root = write_plugin(tmp_path / "site")
    listed = run_rasputin("--help", path=[root])
    assert listed.returncode == 0, listed.stderr
    assert lists_command(listed.stdout, "hello"), listed.stdout
    ran = run_rasputin("hello", path=[root])
    assert (ran.returncode, ran.stdout) == (0, "hello\n"), ran.stderr


def test_without_a_plugin_there_is_no_hello() -> None:
    """The control for the test above: ``hello`` comes from the plugin alone."""
    listed = run_rasputin("--help")
    assert listed.returncode == 0, listed.stderr
    assert not lists_command(listed.stdout, "hello")


def test_plugins_main_runs_the_application(tmp_path: Path) -> None:
    """The console script's target: ``main()`` runs ``application()``."""
    root = write_plugin(tmp_path / "site")
    code = "from tin_engine.plugins import main; main()"
    ran = subprocess.run(
        [sys.executable, "-c", code, "hello"],
        capture_output=True,
        text=True,
        env=child_env(root),
        timeout=60,
        check=False,
    )
    assert (ran.returncode, ran.stdout) == (0, "hello\n"), ran.stderr


def test_the_console_script_and_the_cmake_root_entry_point_are_declared() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert project["scripts"]["rasputin"] == "tin_engine.plugins:main"
    assert project["entry-points"]["cmake.root"] == {"rasputin": "tin_engine"}


# ---------------------------------------------------------------- T2


def test_a_plugin_command_named_like_a_core_command_is_skipped(tmp_path: Path) -> None:
    root = write_plugin(tmp_path / "site", body=MESH, entry="mesh")
    core = run_rasputin("mesh", "--help")
    with_plugin = run_rasputin("mesh", "--help", path=[root])
    assert with_plugin.returncode == 0, with_plugin.stderr
    assert with_plugin.stdout == core.stdout
    named = lines_naming(with_plugin.stderr, "fakeplug")
    assert len(named) == 1, with_plugin.stderr
    assert "mesh" in named[0]


def test_the_core_mesh_is_kept_in_process(on_path: PluginPath) -> None:
    plugins = hook()
    on_path(body=MESH, entry="mesh")
    merged = commands(plugins.application())["mesh"]
    core = commands(core_app)["mesh"]
    # Typer wraps the registered function afresh on every get_command, so the
    # Click callbacks differ even for cli.app twice; the function it wraps does not.
    assert merged.callback.__wrapped__ is core.callback.__wrapped__
    assert merged.help == core.help


# ---------------------------------------------------------------- T3

BROKEN = {
    "import-error": (IMPORT_ERROR, "fakeplug is broken on purpose"),
    "runtime-error": (RUNTIME_ERROR, "fakeplug fails at import on purpose"),
    "not-a-typer": (NOT_TYPER, "Typer"),
}


@pytest.mark.parametrize("variant", sorted(BROKEN))
def test_a_broken_plugin_is_skipped_with_one_line(variant: str, tmp_path: Path) -> None:
    body, reason = BROKEN[variant]
    root = write_plugin(tmp_path / "site", body=body)
    listed = run_rasputin("--help", path=[root])
    assert listed.returncode == 0, listed.stderr
    for name in shown_core_names():
        assert lists_command(listed.stdout, name), (name, listed.stdout)
    named = lines_naming(listed.stderr, "fakeplug")
    assert len(named) == 1, listed.stderr
    assert reason in named[0]
    version = run_rasputin("version", path=[root])
    assert version.returncode == 0, version.stderr
    assert lines_naming(version.stdout, "fakeplug") == []


@pytest.mark.parametrize("variant", sorted(BROKEN))
def test_discover_never_raises_on_a_broken_plugin(variant: str, on_path: PluginPath) -> None:
    plugins = hook()
    baseline = len(plugins.discover())
    on_path(body=BROKEN[variant][0])
    found = plugins.discover()
    assert len(found) == baseline + 1
    app = plugins.application(found)
    assert list(commands(app)) == core_names()


# ---------------------------------------------------------------- T4


def test_two_plugins_with_one_command_name_keep_the_first_by_distribution_name(
    tmp_path: Path,
) -> None:
    # bplug is earlier on the path, so discovery meets it first; the rule is
    # sorted distribution names, so aplug must win all the same.
    b = write_plugin(tmp_path / "b", dist="bplug", module="bplug", body=hello_from("bplug"))
    a = write_plugin(tmp_path / "a", dist="aplug", module="aplug", body=hello_from("aplug"))
    ran = run_rasputin("hello", path=[b, a])
    assert (ran.returncode, ran.stdout) == (0, "hello from aplug\n"), ran.stderr
    named = lines_naming(ran.stderr, "plug")
    assert len(named) == 1, ran.stderr
    assert "bplug" in named[0]
    assert "hello" in named[0]


def test_the_first_by_distribution_name_wins_in_process(on_path: PluginPath) -> None:
    plugins = hook()
    on_path(dist="aplug", module="aplug", body=hello_from("aplug"))
    on_path(dist="bplug", module="bplug", body=hello_from("bplug"))  # first on sys.path
    result = CliRunner().invoke(plugins.application(), ["hello"])
    assert (result.exit_code, result.output) == (0, "hello from aplug\n")


# ---------------------------------------------------------------- T5


def test_version_lists_each_loaded_plugin(on_path: PluginPath) -> None:
    plugins = hook()
    on_path()
    result = CliRunner().invoke(plugins.application(), ["version"])
    assert result.exit_code == 0, result.output
    lines = result.output.splitlines()
    assert lines[1].startswith("bounds checks: ")
    assert "fakeplug 0.1" in lines


def test_python_m_tin_engine_version_lists_the_plugin(tmp_path: Path) -> None:
    root = write_plugin(tmp_path / "site", version="0.1")
    ran = run_rasputin("version", path=[root])
    assert ran.returncode == 0, ran.stderr
    assert "fakeplug 0.1" in ran.stdout.splitlines()


# ---------------------------------------------------------------- T6


async def test_build_server_derives_a_tool_from_the_plugin(
    on_path: PluginPath, monkeypatch: pytest.MonkeyPatch
) -> None:
    plugins = hook()
    pytest.importorskip("mcp")
    from mcp.client import Client

    from mcp_driver import structured
    from tin_engine.mcp.server import build_server

    root = on_path()
    # The server runs each call as ``python -m tin_engine``; the child finds
    # the plugin the same way, through its path.
    monkeypatch.setenv("PYTHONPATH", child_env(root)["PYTHONPATH"])
    server = build_server(plugins.application(plugins.discover()))
    with anyio.fail_after(30):
        async with Client(server) as client:
            listed = await client.list_tools()
            result = await client.call_tool("hello", {})
    assert "hello" in [tool.name for tool in listed.tools]
    assert not result.is_error, result.content
    body = structured(result)
    assert (body["exit_code"], body["stdout"]) == (0, ["hello"])


class TestT8dAPluginsOutputIsFenced:
    """`docs/increments/57-vector-sources.md`, section 4.6, T8d: a plugin's
    ``writable=True`` option is fenced by the server's output root as a core
    command's is. The plugin writes by itself, so the root is made first."""

    async def scribble(
        self, on_path: PluginPath, monkeypatch: pytest.MonkeyPatch, root: Path, out: Path
    ) -> Any:
        plugins = hook()
        pytest.importorskip("mcp")
        from mcp.client import Client

        from tin_engine.mcp.server import build_server

        site = on_path(body=SCRIBBLE, entry="scribble")
        monkeypatch.setenv("PYTHONPATH", child_env(site)["PYTHONPATH"])
        server = build_server(plugins.application(plugins.discover()), out_root=root)
        with anyio.fail_after(30):
            async with Client(server) as client:
                return await client.call_tool("scribble", {"out": str(out)})

    async def test_inside_the_root_it_writes(
        self, on_path: PluginPath, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        from mcp_driver import structured

        root = tmp_path / "root"
        root.mkdir()
        out = root / "w.txt"
        result = await self.scribble(on_path, monkeypatch, root, out)
        assert not result.is_error, result.content
        assert structured(result)["exit_code"] == 0, structured(result)
        assert out.read_text() == "a line\n"

    async def test_outside_the_root_it_is_refused_and_never_runs(
        self, on_path: PluginPath, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        from mcp_driver import structured

        root = tmp_path / "root"
        root.mkdir()
        out = tmp_path / "outside.txt"
        result = await self.scribble(on_path, monkeypatch, root, out)
        assert result.is_error, result.content
        body = structured(result)
        assert body["exit_code"] is None, body
        assert str(body["error"]).startswith("out: "), body
        assert f"outside the output root {root.resolve()}" in body["error"], body
        assert not out.exists()  # the child would have written it


async def test_rasputin_mcp_serves_the_plugin_command(tmp_path: Path) -> None:
    """``python -m tin_engine mcp`` builds its server from ``application()``."""
    pytest.importorskip("mcp")
    from mcp import StdioServerParameters
    from mcp.client import Client

    root = write_plugin(tmp_path / "site")
    params = StdioServerParameters(
        command=sys.executable, args=["-m", "tin_engine", "mcp"], env=child_env(root)
    )
    with anyio.fail_after(30):
        async with Client(params) as client:
            listed = await client.list_tools()
    assert "hello" in [tool.name for tool in listed.tools]


@pytest.mark.parametrize(
    ("body", "entry", "reason"),
    [(IMPORT_ERROR, "hello", "fakeplug is broken on purpose"), (MESH, "mesh", "mesh")],
    ids=["broken", "named-like-core"],
)
def test_rasputin_mcp_reports_each_refusal_once(
    body: str, entry: str, reason: str, tmp_path: Path
) -> None:
    """Section 4.1's one stderr line holds for ``mcp`` too, which builds its
    own ``application()`` inside a CLI that already built one. Stdin is empty,
    so the server sees end of input and exits (measured about 0.7 s)."""
    pytest.importorskip("mcp")
    root = write_plugin(tmp_path / "site", body=body, entry=entry)
    served = subprocess.run(
        [sys.executable, "-m", "tin_engine", "mcp"],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        env=child_env(root),
        timeout=60,
        check=False,
    )
    assert served.returncode == 0, served.stderr
    named = lines_naming(served.stderr, "fakeplug")
    assert len(named) == 1, served.stderr
    assert reason in named[0]
