"""What a newcomer reads first: `rasputin --help`, each command's help, and the plugin page.

`docs/increments/54-open-core-cleanup.md`, *5. Tests*, T1 to T6. The help is
rendered by the CLI itself at a fixed width, so what is checked is what a user
sees, including text a `help=` builds from a constant, which the history gate
does not read. The history patterns are imported from the gate by path, so
there is one list.

Not invariant-critical, no mutation round.
"""

from __future__ import annotations

import importlib.util
import re
import subprocess
import sys
import tomllib
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import typer
from typer.testing import CliRunner

from cli_helpers import USAGE, invoke, plain, rough_dem
from plugin_fixtures import HELLO
from tin_engine import cli, plugin_api
from tin_engine.plugins import Plugin, application

ROOT = Path(__file__).resolve().parents[2]
CLI_PY = "src_python/tin_engine/cli.py"
PLUGINS_PAGE = ROOT / "docs" / "plugins.md"

#: Wide enough that Rich wraps no command name or option name.
COLUMNS = 200

#: `rasputin --help`'s command panels and their commands, in order.
PANELS = [
    ("Meshing", ["mesh", "palette"]),
    ("Catchments", ["catchment", "station-catchments"]),
    ("Data", ["fetch", "fetch-stations"]),
    ("Other", ["version", "mcp"]),
]

PANEL_TOP = re.compile(r"^╭─ (?P<title>.+?) ─+╮?$")
COMMAND_ROW = re.compile(r"^│ (?P<name>[a-z][a-z0-9-]*)(?:\s|$)")
ANSI = re.compile(r"\x1b\[[0-9;]*m")

bumpy = rough_dem(16)


def history_gate() -> ModuleType:
    """`tools/check_code_history.py`, imported by path."""
    spec = importlib.util.spec_from_file_location(
        "check_code_history", ROOT / "tools" / "check_code_history.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rendered(app: typer.Typer, *args: str) -> str:
    """``app``'s help for ``args``, as a wide terminal shows it, colour codes removed.

    `TERM` is not `dumb`, as `cli_helpers.help_text` explains: Rich sizes a
    forced terminal whose `TERM` is `dumb` at 80 columns whatever `COLUMNS` says.
    """
    env = {"NO_COLOR": "1", "TERM": "xterm", "COLUMNS": str(COLUMNS)}
    result = CliRunner(env=env).invoke(app, [*args, "--help"])
    assert result.exit_code == 0, result.output
    return ANSI.sub("", result.output)


def panels(text: str) -> list[tuple[str, list[str]]]:
    """The help's panels that list commands, in order, each with its command names."""
    found: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        line = line.rstrip()
        if top := PANEL_TOP.match(line):
            found.append((top["title"], []))
        elif found and (row := COMMAND_ROW.match(line)):
            found[-1][1].append(row["name"])
    return [(title, names) for title, names in found if names and title != "Options"]


def command_names() -> list[str]:
    """Every command of `cli.app`, hidden ones included."""
    return list(typer.main.get_command(cli.app).commands)  # type: ignore[attr-defined]


def first_paragraph(name: str) -> str:
    """The first paragraph of a command's help: the MCP tool's description too."""
    command: Any = typer.main.get_command(cli.app).commands[name]  # type: ignore[attr-defined]
    return " ".join((command.help or "").split("\n\n")[0].split())


def history_hits(text: str) -> list[tuple[str, str]]:
    """Every match of the gate's rules in ``text``, and of the other three bans."""
    gate = history_gate()
    hits = [
        (rule, m.group(0))
        for rule, pattern in gate.RULES.items()
        for m in pattern.finditer(text)
        if not gate.SCIENTIFIC.fullmatch(m.group(0))
    ]
    hits += [("double backtick", "``")] if "``" in text else []
    hits += [("module path", m.group(0)) for m in re.finditer(r"tin_engine\.\w+", text)]
    hits += [("design file", m.group(0)) for m in gate.MD.finditer(text)]
    return hits


def hello_app(panel: str | None = None) -> typer.Typer:
    """The hello plugin of `plugin_fixtures`, run as its module would be;
    with ``panel``, its command names that help panel."""
    body = HELLO
    if panel is not None:
        body = body.replace("@app.command()", f"@app.command(rich_help_panel={panel!r})")
        assert body != HELLO
    namespace: dict[str, Any] = {}
    exec(compile(body, "fakeplug.py", "exec"), namespace)
    app = namespace["app"]
    assert isinstance(app, typer.Typer)
    return app


@pytest.fixture
def mcp_schema() -> Iterator[ModuleType]:
    pytest.importorskip("mcp")
    from tin_engine.mcp import schema

    yield schema


# ---------------------------------------------------------------- T1


class TestTheGateReadsCli:
    """T1: `cli.py` is in the history gate's scope, and passes it."""

    def test_cli_py_is_in_the_gates_paths(self) -> None:
        config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        assert CLI_PY in config["tool"]["check_code_history"]["paths"]

    def test_the_gate_finds_nothing_in_cli_py(self) -> None:
        result = subprocess.run(
            [sys.executable, "tools/check_code_history.py", CLI_PY],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,  # a hang guard: the gate reads one file in well under a second
            check=False,
        )
        assert (result.returncode, result.stdout, result.stderr) == (0, "", "")


# ---------------------------------------------------------------- T2


class TestTheRenderedHelpCarriesNoHistory:
    """T2: no history, markup, module path or design file name in what is shown."""

    def test_the_root_help(self) -> None:
        text = plain(rendered(cli.app))
        assert history_hits(text) == [], text

    @pytest.mark.parametrize("name", command_names())
    def test_each_commands_help(self, name: str) -> None:
        text = plain(rendered(cli.app, name))
        assert history_hits(text) == [], text

    def test_each_mcp_tool_description(self, mcp_schema: ModuleType) -> None:
        described = {d.tool_name: d.tool().description for d in mcp_schema.derive(cli.app)}
        hits = {name: history_hits(text) for name, text in described.items()}
        assert {name: found for name, found in hits.items() if found} == {}, described

    def test_the_bans_see_what_they_ban(self) -> None:
        """The checks above can fail: each ban catches its own example."""
        examples = (
            "increment 16c",
            "R3",
            "Ola",
            "``--out``",
            "tin_engine.mosaic",
            "13-bundled-mesh.md",
        )
        assert [text for text in examples if not history_hits(text)] == []


# ---------------------------------------------------------------- T3


class TestMeshWithoutTheGallery:
    """T3: `--dem` is required; the gallery, `--flat` and `--crs` are gone from `mesh`."""

    def test_no_dem_is_clicks_missing_option(self, tmp_path: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke("mesh", "--out", str(target))
        assert code == USAGE, output
        assert "Missing option" in output and "--dem" in output, output
        assert "gallery" not in output and "fixture" not in output, output
        assert not target.exists()

    def test_a_fixture_name_is_an_extra_argument(self, tmp_path: Path, bumpy: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke("mesh", "catchment", "--dem", str(bumpy), "--out", str(target))
        assert code == USAGE, output
        assert "unexpected extra argument" in output and "catchment" in output, output
        assert not target.exists()

    @pytest.mark.parametrize("option", [("--flat",), ("--crs", "EPSG:25833")])
    def test_flat_and_crs_are_unknown_options(
        self, tmp_path: Path, bumpy: Path, option: tuple[str, ...]
    ) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke("mesh", "--dem", str(bumpy), *option, "--out", str(target))
        assert code == USAGE, output
        assert f"No such option: {option[0]}" in output, output
        assert not target.exists()

    def test_the_help_names_no_gallery(self) -> None:
        text = plain(rendered(cli.app, "mesh"))
        for word in ("gallery", "fixture", "--flat", "--crs"):
            assert word not in text, word

    def test_the_first_paragraph_names_the_three_outputs(self) -> None:
        first = first_paragraph("mesh")
        for suffix in (".vtk", ".ply", ".gpkg"):
            assert suffix in first, first


# ---------------------------------------------------------------- T4


class TestGroupingAndOrder:
    """T4: the root help groups the commands in panels, in the table's order."""

    def test_the_panels_and_their_commands(self) -> None:
        assert panels(rendered(cli.app)) == PANELS

    def test_draw_is_hidden_and_still_runs(self) -> None:
        assert "draw" not in {n for _, names in panels(rendered(cli.app)) for n in names}
        code, output = invoke("draw", "--list")
        assert code == 0, output
        assert "catchment" in output, "the gallery is listed"

    def test_no_tool_is_offered_for_draw(self, mcp_schema: ModuleType) -> None:
        tools = [d.tool_name for d in mcp_schema.derive(cli.app)]
        assert "draw" not in tools and "mesh" in tools, tools

    def test_a_plugin_without_a_panel_goes_under_plugins_after_the_core(self) -> None:
        plugin = hello_app()
        (info,) = plugin.registered_commands
        unset = info.rich_help_panel  # Typer's placeholder for "no panel named"
        app = application([Plugin("fakeplug", "0.1", "hello", plugin)])
        assert panels(rendered(app)) == [*PANELS, ("Plugins", ["hello"])]
        assert plugin.registered_commands == [info], "the plugin's own app is changed"
        assert info.rich_help_panel is unset, "the plugin's own command is changed"

    def test_a_plugin_that_names_a_panel_keeps_it(self) -> None:
        plugin = hello_app(panel="Snow")
        app = application([Plugin("fakeplug", "0.1", "hello", plugin)])
        shown = dict(panels(rendered(app)))
        assert shown.get("Snow") == ["hello"], shown
        assert "Plugins" not in shown, shown


# ---------------------------------------------------------------- T5


class TestOneDescription:
    """T5: the root help shows `cli.DESCRIPTION`, and `pyproject.toml` says the same."""

    def test_the_root_help_shows_it(self) -> None:
        assert plain(cli.DESCRIPTION) in plain(rendered(cli.app))

    def test_pyproject_says_the_same(self) -> None:
        config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        assert config["project"]["description"] == cli.DESCRIPTION


# ---------------------------------------------------------------- T6


class TestThePluginsPage:
    """T6: `docs/plugins.md` names the entry-point group and every promised name."""

    @staticmethod
    def page() -> str:
        assert PLUGINS_PAGE.is_file(), f"{PLUGINS_PAGE} does not exist"
        return PLUGINS_PAGE.read_text(encoding="utf-8")

    def test_it_names_the_group(self) -> None:
        assert "rasputin.commands" in self.page()

    def test_it_names_every_name_of_the_interface(self) -> None:
        page = self.page()
        missing = [
            name for name in plugin_api.__all__ if not re.search(rf"\b{re.escape(name)}\b", page)
        ]
        assert missing == []
