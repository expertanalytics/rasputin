"""Increment 37, T1 to T3: each tool's schema is derived from its Typer command.

`docs/increments/37-mcp-server.md`, sections 5, 6.3 and 9. `derive(app)`
reads the click parameters of every command and gives, per command, a tool
name, a Pydantic model (extra forbidden) and the way back to argv; the schema
is the model's JSON Schema. T1 and T2 use a toy app with no kernel behind it,
one parameter per row of section 5's type table; T3 derives the real app and
compares each schema with the command's own parameters at test time, so the
schema and the CLI cannot drift apart unnoticed.

Per-property ``title`` (Pydantic's, from the field name) is ignored in every
comparison: it carries nothing the name does not.
"""

from __future__ import annotations

import enum
import json
from pathlib import Path
from typing import Annotated, Any

import pytest

pytest.importorskip("mcp")

import typer
from pydantic import BaseModel, ValidationError

from mcp_driver import tool_names
from tin_engine.cli import app
from tin_engine.mcp import schema as mcp_schema

# ---------------------------------------------------------------- the toy app


class Mode(enum.StrEnum):
    fast = "fast"
    slow = "slow"


def toy_app() -> typer.Typer:
    """Two commands covering section 5's type table, and the server's own command."""
    toy = typer.Typer()

    @toy.command()
    def make_thing(
        name: Annotated[str, typer.Argument(help="What to make.")],
        count: Annotated[int, typer.Option("--count", help="How many.")] = 3,
        level: Annotated[int, typer.Option("--level", min=1, max=5, help="One to five.")] = 2,
        ratio: Annotated[float, typer.Option("--ratio", help="A share.")] = 0.5,
        binary: Annotated[bool, typer.Option("--binary/--ascii", help="Binary or text.")] = True,
        box: Annotated[
            tuple[float, float, float, float] | None,
            typer.Option("--box", help="XMIN YMIN XMAX YMAX."),
        ] = None,
        near: Annotated[tuple[Path, float] | None, typer.Option("--near", help="FILE N.")] = None,
        tag: Annotated[list[str] | None, typer.Option("--tag", help="Repeatable.")] = None,
        mode: Annotated[Mode, typer.Option("--mode", help="A choice.")] = Mode.fast,
        note: Annotated[str | None, typer.Option("--note", help="Optional text.")] = None,
        where: Annotated[Path, typer.Option("--where", help="A path.")] = Path("here"),
    ) -> None:
        """Make a thing.

        The second paragraph is not the description.
        """

    @toy.command()
    def plain() -> None:
        """No parameters at all."""

    @toy.command()
    def mcp() -> None:
        """The server's own command: not a tool."""

    return toy


NUMBER = {"type": "number"}
STRING = {"type": "string"}
#: Section 5's type table, one property per row, as JSON Schema without ``title``.
EXPECTED_PROPERTIES: dict[str, dict[str, Any]] = {
    "name": {**STRING, "description": "What to make."},
    "count": {"type": "integer", "default": 3, "description": "How many."},
    "level": {
        "type": "integer", "minimum": 1, "maximum": 5, "default": 2, "description": "One to five.",
    },
    "ratio": {**NUMBER, "default": 0.5, "description": "A share."},
    "binary": {"type": "boolean", "default": True, "description": "Binary or text."},
    "box": {
        "anyOf": [
            {"type": "array", "prefixItems": [NUMBER] * 4, "minItems": 4, "maxItems": 4},
            {"type": "null"},
        ],
        "default": None,
        "description": "XMIN YMIN XMAX YMAX.",
    },
    "near": {
        "anyOf": [
            {"type": "array", "prefixItems": [STRING, NUMBER], "minItems": 2, "maxItems": 2},
            {"type": "null"},
        ],
        "default": None,
        "description": "FILE N.",
    },
    "tag": {
        "anyOf": [{"type": "array", "items": STRING}, {"type": "null"}],
        "default": None,
        "description": "Repeatable.",
    },
    "mode": {
        "enum": ["fast", "slow"], "type": "string", "default": "fast", "description": "A choice.",
    },
    "note": {"anyOf": [STRING, {"type": "null"}], "default": None, "description": "Optional text."},
    "where": {**STRING, "default": "here", "description": "A path."},
}  # fmt: skip


def without_titles(properties: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        name: {k: v for k, v in prop.items() if k != "title"} for name, prop in properties.items()
    }


def derived(application: typer.Typer) -> dict[str, Any]:
    """`derive(app)` by tool name."""
    return {d.tool_name: d for d in mcp_schema.derive(application)}


def click_command(application: typer.Typer, tool_name: str) -> Any:
    return typer.main.get_command(application).commands[tool_name.replace("_", "-")]


def z_unit_refused(model: type[BaseModel], value: str | None) -> bool:
    """Whether validation names `z_unit`; the other fields are left out, so they may fail too."""
    try:
        model.model_validate({"z_unit": value})
    except ValidationError as error:
        return any(e["loc"] == ("z_unit",) for e in error.errors())
    return False


# ---------------------------------------------------------------- T1


class TestSchemaFromAToyApp:
    def test_one_tool_per_command_named_with_underscores_and_the_server_excluded(self) -> None:
        tools = derived(toy_app())
        assert list(tools) == ["make_thing", "plain"]
        assert tools["make_thing"].tool().name == "make_thing"

    def test_every_row_of_the_type_table(self) -> None:
        schema = derived(toy_app())["make_thing"].tool().input_schema
        assert schema["type"] == "object"
        assert "title" not in schema
        assert schema["additionalProperties"] is False
        assert schema["required"] == ["name"]
        assert without_titles(schema["properties"]) == EXPECTED_PROPERTIES

    def test_a_command_without_parameters(self) -> None:
        schema = derived(toy_app())["plain"].tool().input_schema
        assert schema["properties"] == {}
        assert schema.get("required", []) == []
        assert schema["additionalProperties"] is False

    def test_the_description_is_the_docstrings_first_paragraph(self) -> None:
        tools = derived(toy_app())
        assert tools["make_thing"].tool().description == "Make a thing."
        assert tools["plain"].tool().description == "No parameters at all."

    @pytest.mark.parametrize(
        "arguments",
        [
            pytest.param({}, id="missing required argument"),
            pytest.param({"name": "a", "box": [1.0, 2.0, 3.0]}, id="tuple of the wrong arity"),
            pytest.param({"name": "a", "count": "many"}, id="wrong type"),
            pytest.param({"name": "a", "bogus": 1}, id="unknown key"),
            pytest.param({"name": "a", "mode": "medium"}, id="outside the choice"),
            pytest.param({"name": "a", "level": 6}, id="above the range"),
        ],
    )
    def test_the_model_refuses_what_the_cli_would(self, arguments: dict[str, Any]) -> None:
        model = derived(toy_app())["make_thing"].model
        with pytest.raises(ValidationError):
            model.model_validate(arguments)


# ---------------------------------------------------------------- T2


class TestArgvRoundTrip:
    @pytest.mark.parametrize(
        ("arguments", "argv"),
        [
            pytest.param({"name": "a"}, ["make-thing", "a"], id="only set keys appear"),
            pytest.param(
                {"name": "a", "binary": False}, ["make-thing", "a", "--ascii"], id="false flag"
            ),
            pytest.param(
                {"name": "a", "binary": True}, ["make-thing", "a", "--binary"], id="true flag"
            ),
            pytest.param(
                {"name": "a", "box": [1.5, 2.5, 3.5, 4.5]},
                ["make-thing", "a", "--box", "1.5", "2.5", "3.5", "4.5"],
                id="a tuple expands",
            ),
            pytest.param(
                {"name": "a", "near": ["lines.geojson", 2.5]},
                ["make-thing", "a", "--near", "lines.geojson", "2.5"],
                id="a path-and-number tuple",
            ),
            pytest.param(
                {"name": "a", "tag": ["x", "y"]},
                ["make-thing", "a", "--tag", "x", "--tag", "y"],
                id="a list repeats the option",
            ),
            pytest.param(
                {"mode": "slow", "count": 7, "name": "a"},
                ["make-thing", "a", "--count", "7", "--mode", "slow"],
                id="the command's parameter order, not the caller's",
            ),
            pytest.param(
                {"name": "a", "where": "there/x.vtk"},
                ["make-thing", "a", "--where", "there/x.vtk"],
                id="a path is passed as given",
            ),
        ],
    )
    def test_argv_is_the_command_line_the_arguments_mean(
        self, arguments: dict[str, Any], argv: list[str]
    ) -> None:
        tool = derived(toy_app())["make_thing"]
        assert tool.argv(tool.model.model_validate(arguments)) == argv

    def test_a_command_without_parameters_is_its_name(self) -> None:
        tool = derived(toy_app())["plain"]
        assert tool.argv(tool.model.model_validate({})) == ["plain"]


# ---------------------------------------------------------------- T3

#: Section 1's eight commands as tool names. Their parameter counts are not
#: pinned: the test below holds each schema's properties equal to the command's
#: own click parameters, which implies the count, and a pinned count would turn
#: red on every master merge that adds a CLI option (as `mesh` 34 -> 36 did when
#: `--land-cover` and `--land-cover-year` landed) without saying anything about
#: the schema.
TOOL_NAMES = {
    "mesh", "draw", "fetch", "catchment", "station_catchments",
    "fetch_stations", "palette", "version",
}  # fmt: skip

#: Section 3's snapshots of five `mesh` properties, without ``title``.
MESH_SNAPSHOTS: dict[str, dict[str, Any]] = {
    "bbox": {
        "anyOf": [
            {"maxItems": 4, "minItems": 4, "prefixItems": [NUMBER] * 4, "type": "array"},
            {"type": "null"},
        ],
        "default": None,
        "description": "With --dem, mesh only XMIN YMIN XMAX YMAX in the DEM's CRS, snapped "
        "outward to its nodes.",
    },
    "tolerance_near": {
        "anyOf": [
            {"maxItems": 2, "minItems": 2, "prefixItems": [STRING, NUMBER], "type": "array"},
            {"type": "null"},
        ],
        "default": None,
        "description": "Hold the mesh to N metres on these lines (.geojson, .gpkg or .gml). "
        "Needs --tolerance-ramp.",
    },
    "binary": {
        "default": True,
        "description": "Binary files, small and fast; --ascii writes text you can read.",
        "type": "boolean",
    },
    "features": {
        "anyOf": [{"items": STRING, "type": "array"}, {"type": "null"}],
        "default": None,
        "description": "With --domain, polygons and lines as constraints (.geojson, .json, "
        ".gpkg or .gml), clipped to the domain, each edge carrying its class map's bits. "
        "Repeatable: each --features is one source, paired by position with the Nth of "
        "--features-crs, --features-layer and --features-map.",
    },
    "snap_spacing": {
        "default": 0.001,
        "description": "Snap grid spacing for the noder, in the input's units.",
        "type": "number",
    },
}

RUN_RESULT_FIELDS = {
    "command",
    "exit_code",
    "seconds",
    "outputs",
    "stdout",
    "stderr",
    "record",
    "stats",
}


class TestTheRealApp:
    def test_the_eight_commands_in_the_clis_order(self) -> None:
        tools = derived(app)
        assert list(tools) == tool_names(app)
        assert set(tools) == TOOL_NAMES

    def test_no_drift_from_the_click_parameters(self) -> None:
        tools = derived(app)
        for name, d in tools.items():
            params = {p.name for p in click_command(app, name).params}
            assert set(d.tool().input_schema["properties"]) == params, name

    def test_mesh_requires_out_and_matches_the_snapshots(self) -> None:
        schema = derived(app)["mesh"].tool().input_schema
        assert schema["required"] == ["out"]
        assert schema["additionalProperties"] is False
        properties = without_titles(schema["properties"])
        assert {k: properties[k] for k in MESH_SNAPSHOTS} == MESH_SNAPSHOTS

    def test_an_integer_range_carries_its_minimum(self) -> None:
        connections = derived(app)["fetch"].tool().input_schema["properties"]["connections"]
        assert {k: connections[k] for k in ("type", "minimum", "default")} == {
            "type": "integer", "minimum": 1, "default": 8,
        }  # fmt: skip

    def test_t9_threads_and_memory_budget_derive_from_the_cli(self) -> None:
        """Increment 38, T9: the two options reach the `mesh`, `catchment` and
        `station_catchments` tools through the derived schema; `argv` puts
        them on the command line as the CLI's words. The `minimum: 0` bound
        is kept: a bound of zero is a bound (`schema.py` dropped every falsy
        bound before 38, the one `mcp/` line the increment changes)."""
        tools = derived(app)
        for name in ("mesh", "catchment", "station_catchments"):
            properties = tools[name].tool().input_schema["properties"]
            threads, budget = properties["threads"], properties["memory_budget"]
            assert threads["default"] is None, name
            assert {"type": "integer", "minimum": 0} in threads["anyOf"], name
            assert {"type": "null"} in threads["anyOf"], name
            assert budget["default"] is None and {"type": "string"} in budget["anyOf"], name
        mesh = tools["mesh"]
        with_threads = mesh.model.model_validate({"out": "x.vtk", "threads": 3})
        assert mesh.argv(with_threads)[-2:] == ["--threads", "3"]
        with_budget = mesh.model.model_validate({"out": "x.vtk", "memory_budget": "8G"})
        assert mesh.argv(with_budget)[-2:] == ["--memory-budget", "8G"]
        with pytest.raises(ValidationError):
            mesh.model.model_validate({"out": "x.vtk", "threads": -1})

    @pytest.mark.parametrize("name", ["mesh", "catchment", "station_catchments"])
    def test_z_unit_admits_the_three_spellings_and_null(self, name: str) -> None:
        """Increment 48, T3: the schema lists the spellings, and the model holds to them."""
        tool = derived(app)[name]
        branches = tool.tool().input_schema["properties"]["z_unit"].get("anyOf", [])
        strings = [b for b in branches if b.get("type") == "string"]
        listed = {v for b in strings for v in b.get("enum", [b.get("const")])}
        assert listed == {"m", "ft", "us-ft"}, branches
        assert {"type": "null"} in branches, branches
        assert [s for s in ("m", "ft", "us-ft", None) if z_unit_refused(tool.model, s)] == []
        assert all(z_unit_refused(tool.model, s) for s in ("furlong", "M", "us_ft"))

    def test_the_whole_tool_list_stays_small(self) -> None:
        # 20,141 bytes measured 2026-10-13, after --land-cover-map (increment 44);
        # the bound keeps the original ~20 % headroom (16,575 against 20,000) for help text.
        size = sum(len(json.dumps(d.tool().input_schema)) for d in mcp_schema.derive(app))
        assert size < 24_000, size

    def test_descriptions_and_annotations(self) -> None:
        tools = {name: d.tool() for name, d in derived(app).items()}
        assert tools["version"].description == (
            "Print the installed rasputin version, then whether it is bounds-checked."
        )
        version = tools["version"].annotations
        assert version is not None
        assert (version.read_only_hint, version.idempotent_hint) == (True, True)
        for name in ("fetch", "fetch_stations"):
            annotations = tools[name].annotations
            assert annotations is not None and annotations.open_world_hint is True, name
        for name, tool in tools.items():
            assert tool.annotations is not None and tool.annotations.destructive_hint is False, name

    def test_every_tool_declares_the_run_result_as_its_output(self) -> None:
        for d in mcp_schema.derive(app):
            output = d.tool().output_schema
            assert output is not None, d.tool_name
            assert set(output["properties"]) >= RUN_RESULT_FIELDS, d.tool_name
