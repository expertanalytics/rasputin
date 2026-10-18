"""The old "class map" names of a land-cover legend are refused, not aliased."""

from __future__ import annotations

from pathlib import Path

from cli_helpers import USAGE, invoke
from test_catalogue_data import PER_YEAR, USER_ID, changed, refusal, write
from tin_engine import sources
from tin_engine.cli import app
from tin_engine.mcp import schema as mcp_schema


def test_the_old_mesh_options_are_refused() -> None:
    for option, value in (("--features-map", "corine"), ("--land-cover-map", "mapbiomas-c11")):
        code, output = invoke("mesh", option, value)
        assert code == USAGE, output
        assert f"No such option: {option}" in output, output


def test_a_sources_file_with_the_old_key_is_refused(tmp_path: Path) -> None:
    user = write(
        tmp_path / "my-sources.toml", {f"source.{USER_ID}": changed(PER_YEAR, class_map="corine")}
    )
    refusal(sources, user, "class_map")


def test_the_mcp_mesh_tool_takes_only_the_new_names() -> None:
    (mesh,) = (d for d in mcp_schema.derive(app) if d.tool_name == "mesh")
    properties = set(mesh.tool().input_schema["properties"])
    assert {"features_legend", "land_cover_legend"} <= properties, sorted(properties)
    assert properties.isdisjoint({"features_map", "land_cover_map"}), sorted(properties)
