"""An agent with only the MCP tools finds the DEM, the lakes and the land cover, and uses them.

`docs/increments/57-vector-sources.md`, sections 4.6 and 4.9, tests T7 (the
acceptance run over MCP only) and T10 (the instructions).

The server is the production one, `build_server(cli.app)` with the real
``python -m tin_engine`` child, in process. The home is a temporary directory
whose `~/.config/rasputin/sources.toml` registers fixture L, fixture V's lake
file and its land-cover GeoPackage under ids that name nothing real;
`RASPUTIN_DATA` and `XDG_CONFIG_HOME` are unset, as for a client that sets no
`env`. After that setup the test acts only through tool calls: every path,
layer, class map and id it passes is read off a tool result. Only the oracle
at the end reads the files written.

The output folders do not exist before the calls: the tools make them.

Bounds: each call is one child run, measured under 0.11 s of work on these
fixtures plus the interpreter's start.
"""

from __future__ import annotations

import json
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("mcp")

from mcp.client import Client

from catchment_fixtures import EPSG
from dataset_fixtures import INSIDE, LOCAL, SEED, local_entry, write_local, write_sources
from mcp_driver import call, structured
from tin_engine.cli import app
from tin_engine.mcp import server as mcp_server
from vector_fixtures import (
    CATCHMENT_NODES,
    COVER_DESCRIPTION,
    COVER_ID,
    LAKE_DESCRIPTION,
    LAKES_ID,
    LANDCOVER,
    MESH_TRIANGLES,
    TRIANGLES_BY_CODE,
    fixture_v,
    vector_entry,
)


@pytest.fixture
def home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "home"
    lakes, cover = fixture_v(tmp_path / "data" / "vectors")
    write_sources(
        root / ".config" / "rasputin" / "sources.toml",
        {
            LOCAL: local_entry(write_local(tmp_path / "data" / "bowl")),
            LAKES_ID: vector_entry(lakes, LAKE_DESCRIPTION),
            COVER_ID: vector_entry(cover, COVER_DESCRIPTION, layer=LANDCOVER, class_map="corine"),
        },
    )
    monkeypatch.setenv("HOME", str(root))
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    return root


def server() -> Any:
    return mcp_server.build_server(app)


def body_of(result: Any) -> dict[str, Any]:
    assert not result.is_error, result.content
    body = structured(result)
    assert body["exit_code"] == 0, body
    return body


def fitting(document: dict[str, Any], tool: str, argument: str) -> list[dict[str, Any]]:
    """The listed sources one of whose fits names `tool` and `argument`."""
    return [
        s
        for s in document["sources"]
        if any((f["tool"], f["argument"]) == (tool, argument) for f in s["fits"])
    ]


def the_fit(source: dict[str, Any], tool: str, argument: str) -> dict[str, Any]:
    (fit,) = [f for f in source["fits"] if (f["tool"], f["argument"]) == (tool, argument)]
    return fit


# ------------------------------------------------------------ T7


async def test_t7_find_then_catchment_then_mesh_over_mcp_only(home: Path, tmp_path: Path) -> None:
    listed = body_of(await call(server(), "sources", {"bbox": list(INSIDE)}))
    (line,) = listed["stdout"]
    document = json.loads(line)

    dems = [s for s in fitting(document, "catchment", "dem") if s["box_covered"] == 1.0]
    (dem,) = dems
    # Land cover fits `lakes` too (it holds polygons); what the features are
    # is in the description, which the agent reads.
    polygons = fitting(document, "catchment", "lakes")
    (lakes,) = [s for s in polygons if "lake" in s["description"].lower()]
    (cover,) = fitting(document, "mesh", "features")
    features_map = the_fit(cover, "mesh", "features")["class_map"]

    outline = tmp_path / "new" / "deeper" / "c.geojson"
    assert not outline.parent.exists()
    catchment = body_of(await call(server(), "catchment", {
        "dem": [dem["id"]], "seed": list(SEED), "seed_crs": EPSG, "lakes": lakes["path"],
        "out": str(outline),
    }))  # fmt: skip
    (written,) = [o for o in catchment["outputs"] if o.endswith(".geojson")]

    out = tmp_path / "meshes" / "run1" / "m.gpkg"
    assert not out.parent.parent.exists()
    mesh = body_of(await call(server(), "mesh", {
        "dem": [dem["id"]], "domain": written, "tolerance": 5.0,
        "features": [cover["path"]], "features_layer": [cover["layer"]],
        "features_map": [features_map], "out": str(out),
    }))  # fmt: skip
    record = mesh["record"]
    assert isinstance(record, dict), mesh
    said = json.dumps(record["features"])
    assert "layer landcover" in said and "class map corine" in said, said

    # The oracle: the files written, as measured.
    assert (dem["id"], lakes["id"], cover["id"]) == (LOCAL, LAKES_ID, COVER_ID)
    (feature,) = json.loads(outline.read_text())["features"]
    assert feature["properties"]["nodes"] == CATCHMENT_NODES
    with closing(sqlite3.connect(out)) as con:
        (count,) = con.execute("SELECT count(*) FROM triangles").fetchone()
        codes = dict(con.execute("SELECT code, count(*) FROM triangles GROUP BY code").fetchall())
    assert count == MESH_TRIANGLES
    assert codes == TRIANGLES_BY_CODE


# ------------------------------------------------------------ T10


class TestT10Instructions:
    async def test_the_vector_guidance(self) -> None:
        async with Client(server()) as client:
            instructions = client.instructions
        assert instructions is not None
        for words in ("fits", "Pass its path", "Output folders are created as needed"):
            assert words in instructions, (words, instructions)

    def test_no_source_id_is_named(self) -> None:
        from tin_engine.sources import SOURCES

        text = mcp_server.INSTRUCTIONS
        assert "Output folders are created as needed" in text  # the new guidance is there
        for key in SOURCES:
            assert key not in text, key
        for named in ("corine", "lakes.geojson", "norway"):
            assert named not in text.lower(), named
