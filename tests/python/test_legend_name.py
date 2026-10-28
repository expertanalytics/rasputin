"""Old names of the meshing rules are refused, not aliased.

The "class map" names, and the "legend" names of what is now the meshing
rules (`docs/increments/108-legend-split.md`, *4. The design*, the Names
table): each old option and sources-file key is refused, the newer ones
with a message naming the new name.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from cli_helpers import USAGE, invoke
from test_catalogue_data import PER_YEAR, USER_ID, changed, refusal, write
from tin_engine import sources
from tin_engine.cli import app
from tin_engine.mcp import schema as mcp_schema
from tin_engine.mcp.server import INSTRUCTIONS

COLOURS_ONLY = "a legend is the colour scheme only"


def test_the_class_map_options_are_refused() -> None:
    for option, value in (("--features-map", "corine"), ("--land-cover-map", "mapbiomas-c11")):
        code, output = invoke("mesh", option, value)
        assert code == USAGE, output
        assert f"No such option: {option}" in output, output


@pytest.mark.parametrize(
    ("old", "new", "value"),
    [
        ("--features-legend", "--features-rules", "corine"),
        ("--land-cover-legend", "--land-cover-classes", "mapbiomas-c11"),
    ],
)
def test_the_legend_options_are_refused_naming_the_new(old: str, new: str, value: str) -> None:
    code, output = invoke("mesh", old, value)
    assert code == USAGE, output
    assert f"{old} is now {new}; {COLOURS_ONLY}" in " ".join(output.split()), output
    assert output.count(old) == 1, output


def test_a_sources_file_with_the_class_map_key_is_refused(tmp_path: Path) -> None:
    user = write(
        tmp_path / "my-sources.toml", {f"source.{USER_ID}": changed(PER_YEAR, class_map="corine")}
    )
    refusal(sources, user, "class_map")


def test_a_sources_file_with_the_legend_key_is_refused_naming_rules(tmp_path: Path) -> None:
    user = write(
        tmp_path / "my-sources.toml", {f"source.{USER_ID}": changed(PER_YEAR, legend="corine")}
    )
    refusal(sources, user, "legend is now rules (a source names its meshing rules)")


def test_the_mcp_mesh_tool_takes_only_the_new_names() -> None:
    (mesh,) = (d for d in mcp_schema.derive(app) if d.tool_name == "mesh")
    schema = mesh.tool().input_schema
    properties = set(schema["properties"])
    assert {"features_rules", "land_cover_rules"} <= properties, sorted(properties)
    old = {"features_map", "land_cover_map", "features_legend", "land_cover_legend"}
    assert properties.isdisjoint(old), sorted(properties)
    assert schema.get("additionalProperties") is False


def test_the_mcp_instructions_name_the_new_arguments() -> None:
    assert "features_rules" in INSTRUCTIONS
    assert "features_legend" not in INSTRUCTIONS and "land_cover_legend" not in INSTRUCTIONS
