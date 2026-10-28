"""Old names of the land-cover class file are refused, not aliased.

The "class map" names, the "legend" names (`docs/increments/108-legend-split.md`,
*4. The design*, the Names table) and the meshing-rules and `--features`
names (`docs/increments/111-one-land-cover-map.md`, section 2: `RENAMED` and
`_RenamedOptions` go with their targets, so every old option is simply
unknown, and an old sources-file key is an unknown key like any other).

PINNED HERE, where the design leaves it open (listed for `@architect`): the
MCP instructions name the arguments `land_cover_classes`, `prefer`, `lines`,
`edges` and `shares`, and no `features_*` argument.
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

OLD_OPTIONS = (
    ("--features-map", "corine"),
    ("--land-cover-map", "mapbiomas-c11"),
    ("--features-legend", "corine"),
    ("--land-cover-legend", "mapbiomas-c11"),
    ("--features-rules", "corine"),
    ("--land-cover-rules", "corine"),
    ("--features", "x.geojson"),
    ("--features-layer", "x"),
    ("--features-crs", "EPSG:25833"),
)


@pytest.mark.parametrize(("option", "value"), OLD_OPTIONS, ids=[o for o, _ in OLD_OPTIONS])
def test_an_old_option_is_unknown(option: str, value: str) -> None:
    code, output = invoke("mesh", option, value)
    assert code == USAGE, output
    assert f"No such option: {option}" in output, output


def test_a_sources_file_with_the_class_map_key_is_refused(tmp_path: Path) -> None:
    user = write(
        tmp_path / "my-sources.toml", {f"source.{USER_ID}": changed(PER_YEAR, class_map="corine")}
    )
    refusal(sources, user, "class_map")


@pytest.mark.parametrize("key", ["legend", "rules"])
def test_a_sources_file_with_an_old_key_is_refused_as_unknown(tmp_path: Path, key: str) -> None:
    """`legend` (once named `rules`) and `rules` (now `classes`) are unknown keys."""
    user = write(
        tmp_path / "my-sources.toml", {f"source.{USER_ID}": changed(PER_YEAR, **{key: "corine"})}
    )
    message = refusal(sources, user, key)
    assert "is now" not in message, message


def test_the_mcp_mesh_tool_takes_only_the_new_names() -> None:
    (mesh,) = (d for d in mcp_schema.derive(app) if d.tool_name == "mesh")
    schema = mesh.tool().input_schema
    properties = set(schema["properties"])
    assert {"land_cover_classes", "edges", "lines", "shares"} <= properties, sorted(properties)
    old = {
        *("features_map", "land_cover_map", "features_legend", "land_cover_legend"),
        *("features_rules", "land_cover_rules", "features", "features_layer", "features_crs"),
    }
    assert properties.isdisjoint(old), sorted(properties)
    assert schema.get("additionalProperties") is False


def test_the_mcp_instructions_name_the_new_arguments() -> None:
    for name in ("land_cover_classes", "prefer", "lines", "edges", "shares"):
        assert name in INSTRUCTIONS, name
    for old in ("features_rules", "features_layer", "features_legend", "land_cover_legend"):
        assert old not in INSTRUCTIONS, old
