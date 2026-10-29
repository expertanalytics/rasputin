"""The three border shifts, named by what the number is (Ola, 2026-10-29).

`--edges-max-shift` (land-cover borders), `--catchments-max-shift`
(sub-catchment borders) and `--domain-max-shift` (the `--simplify-domain`
outline) are each the farthest, in metres, a border may move from where it
was when simplified. Each help says that, what 0 means and the default, so an
agent can use the option from its help alone; the MCP mesh tool, built from
the CLI options, carries the three and none of the old names.
"""

from __future__ import annotations

import pytest
import typer.main

from cli_helpers import USAGE, invoke
from tin_engine.cli import app

SHIFTS = ("edges_max_shift", "catchments_max_shift", "domain_max_shift")
OLD = ("edges_tolerance", "features_tolerance", "features_band", "catchments_band", "domain_band")


def mesh_helps() -> dict[str, str]:
    command = typer.main.get_command(app).commands["mesh"]  # type: ignore[attr-defined]
    return {p.name: getattr(p, "help", None) or "" for p in command.params if p.name}


@pytest.mark.parametrize("name", SHIFTS)
def test_each_help_says_what_the_number_is_what_0_means_and_the_default(name: str) -> None:
    text = mesh_helps()[name]
    assert text.startswith("Metres: "), text
    assert "the farthest" in text and "may move from where it was" in text, text
    assert "0: no simplification" in text, text
    assert "Default: " in text, text


def test_the_old_names_are_gone_from_the_cli() -> None:
    helps = mesh_helps()
    assert not set(OLD) & set(helps), sorted(helps)
    for text in helps.values():
        for old in OLD:
            assert old.replace("_", "-") not in text, text


def test_the_mcp_mesh_tool_carries_the_three_shifts_and_no_old_name() -> None:
    pytest.importorskip("mcp")
    from tin_engine.mcp import schema as mcp_schema

    mesh_tool = {d.tool_name: d for d in mcp_schema.derive(app)}["mesh"]
    properties = mesh_tool.tool().input_schema["properties"]
    assert set(SHIFTS) <= set(properties), sorted(properties)
    assert not set(OLD) & set(properties), sorted(properties)
    for name in SHIFTS:
        assert "the farthest" in properties[name].get("description", ""), properties[name]


@pytest.mark.parametrize("old", ["--edges-tolerance", "--catchments-band", "--domain-band"])
def test_the_old_options_are_unknown(old: str) -> None:
    code, output = invoke("mesh", old, "5")
    assert code == USAGE, output
    assert f"No such option: {old}" in output, output
