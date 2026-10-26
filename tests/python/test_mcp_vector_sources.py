"""An agent with only the MCP tools finds the DEM, the lakes and the land cover, and uses them.

`docs/increments/57-vector-sources.md`, sections 4.6 and 4.9, tests T7 (the
acceptance run over MCP only), T8c (the output root's fence), T8e (the root
reported) and T10 (the instructions).

The server is the production one, `build_server(cli.app)` with the real
``python -m tin_engine`` child, in process. The home is a temporary directory
whose `~/.config/rasputin/sources.toml` registers fixture L, fixture V's lake
file and its land-cover GeoPackage under ids that name nothing real;
`RASPUTIN_DATA` and `XDG_CONFIG_HOME` are unset, as for a client that sets no
`env`. After that setup the test acts only through tool calls: every path,
layer, legend and id it passes, and the output root it writes in, is read
off a tool result. Only the oracle at the end reads the files written.

The server is built with no output root, so the default `<HOME>/rasputin-output`
applies; neither it nor the output folders exist before the calls: the tools
make them.

T8c runs the server with `mcp_fake_child`, which records its argv: a refused
call is an error result with no exit code, the child never runs (its argv
file stays absent) and nothing is made (`rglob` before and after).

Bounds: each call is one child run, measured under 0.11 s of work on these
fixtures plus the interpreter's start.
"""

from __future__ import annotations

import json
import sqlite3
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anyio
import pytest

pytest.importorskip("mcp")

from mcp.client import Client

from catchment_fixtures import EPSG
from cli_helpers import USAGE, invoke, runner
from dataset_fixtures import INSIDE, LOCAL, SEED, local_entry, write_local, write_sources
from mcp_driver import call, fake_child, option_value, production_server, structured
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
            COVER_ID: vector_entry(cover, COVER_DESCRIPTION, layer=LANDCOVER, legend="corine"),
        },
    )
    monkeypatch.setenv("HOME", str(root))
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    return root


def server() -> Any:
    """The production server with no output root: the default applies."""
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


async def test_t7_find_then_catchment_then_mesh_over_mcp_only(home: Path) -> None:
    default_root = (home / "rasputin-output").resolve()
    assert not default_root.exists()
    listed = body_of(await call(server(), "sources", {"bbox": list(INSIDE)}))
    (line,) = listed["stdout"]
    document = json.loads(line)
    assert document["out_root"] == str(default_root)
    out_root = default_root  # the listing's out_root, asserted equal above

    dems = [s for s in fitting(document, "catchment", "dem") if s["box_covered"] == 1.0]
    (dem,) = dems
    # Land cover fits `lakes` too (it holds polygons); what the features are
    # is in the description, which the agent reads.
    polygons = fitting(document, "catchment", "lakes")
    (lakes,) = [s for s in polygons if "lake" in s["description"].lower()]
    (cover,) = fitting(document, "mesh", "features")
    features_legend = the_fit(cover, "mesh", "features")["legend"]

    outline = out_root / "new" / "deeper" / "c.geojson"
    assert not out_root.exists()
    catchment = body_of(await call(server(), "catchment", {
        "dem": [dem["id"]], "seed": list(SEED), "seed_crs": EPSG, "lakes": lakes["path"],
        "out": str(outline),
    }))  # fmt: skip
    (written,) = [o for o in catchment["outputs"] if o.endswith(".geojson")]

    out = out_root / "meshes" / "run1" / "m.gpkg"
    assert not out.parent.parent.exists()
    mesh = body_of(await call(server(), "mesh", {
        "dem": [dem["id"]], "domain": written, "tolerance": 5.0,
        "features": [cover["path"]], "features_layer": [cover["layer"]],
        "features_legend": [features_legend], "out": str(out),
    }))  # fmt: skip
    record = mesh["record"]
    assert isinstance(record, dict), mesh
    said = json.dumps(record["features"])
    assert "layer landcover" in said and "legend corine" in said, said

    # The oracle: the files written, as measured.
    assert (dem["id"], lakes["id"], cover["id"]) == (LOCAL, LAKES_ID, COVER_ID)
    (feature,) = json.loads(outline.read_text())["features"]
    assert feature["properties"]["nodes"] == CATCHMENT_NODES
    with closing(sqlite3.connect(out)) as con:
        (count,) = con.execute("SELECT count(*) FROM triangles").fetchone()
        codes = dict(con.execute("SELECT code, count(*) FROM triangles GROUP BY code").fetchall())
    assert count == MESH_TRIANGLES
    assert codes == TRIANGLES_BY_CODE


# ------------------------------------------------------------ T8c


@dataclass
class Fence:
    """A server over a fake child, fenced at ``root``; ``argv_file`` appears only if it ran."""

    root: Path
    argv_file: Path
    server: Any


def fenced(tmp_path: Path, root: Path, tag: str = "argv") -> Fence:
    argv_file = tmp_path / f"{tag}.json"
    child = fake_child(tmp_path, argv_file=str(argv_file))
    return Fence(root, argv_file, mcp_server.build_server(app, out_root=root, child=child))


def everything_under(directory: Path) -> list[Path]:
    return sorted(directory.rglob("*"))


def assert_refused(result: Any, argument: str, root: Path) -> str:
    """An error result before the child: no exit code, the argument and the root named."""
    assert result.is_error, result.content
    body = structured(result)
    assert body["exit_code"] is None, body
    error = str(body["error"])
    assert error.startswith(f"{argument}: "), error
    assert f"outside the output root {root.resolve()}" in error, error
    assert "rasputin mcp --out-root" in error, error
    return error


def mesh_into(out: Path | str, **more: Any) -> dict[str, Any]:
    return {"dem": ["d.tif"], "out": str(out), **more}


class TestT8cTheFence:
    async def test_inside_the_absent_root_the_child_runs(self, tmp_path: Path) -> None:
        fence = fenced(tmp_path, tmp_path / "root")
        out = fence.root / "a" / "m.vtk"
        result = await call(fence.server, "mesh", mesh_into(out))
        assert structured(result)["exit_code"] == 0, structured(result)
        assert option_value(json.loads(fence.argv_file.read_text()), "--out") == str(out)

    async def test_dotdot_out_of_the_root_is_refused(self, tmp_path: Path) -> None:
        fence = fenced(tmp_path, tmp_path / "root")
        given = fence.root / ".." / "escape" / "m.vtk"
        before = everything_under(tmp_path)
        result = await call(fence.server, "mesh", mesh_into(given))
        error = assert_refused(result, "out", fence.root)
        assert f"{given} resolves to {(tmp_path / 'escape' / 'm.vtk').resolve()}" in error, error
        assert not fence.argv_file.exists()
        assert everything_under(tmp_path) == before

    async def test_a_symlinked_folder_pointing_outside_is_refused(self, tmp_path: Path) -> None:
        fence = fenced(tmp_path, tmp_path / "root")
        (tmp_path / "outside").mkdir()
        fence.root.mkdir()
        (fence.root / "link").symlink_to(tmp_path / "outside")
        before = everything_under(tmp_path)
        result = await call(fence.server, "mesh", mesh_into(fence.root / "link" / "m.vtk"))
        assert_refused(result, "out", fence.root)
        assert not fence.argv_file.exists()
        assert everything_under(tmp_path) == before

    async def test_a_symlinked_file_pointing_outside_is_refused(self, tmp_path: Path) -> None:
        fence = fenced(tmp_path, tmp_path / "root")
        (tmp_path / "outside").mkdir()
        fence.root.mkdir()
        (fence.root / "l.vtk").symlink_to(tmp_path / "outside" / "real.vtk")
        before = everything_under(tmp_path)
        result = await call(fence.server, "mesh", mesh_into(fence.root / "l.vtk"))
        assert_refused(result, "out", fence.root)
        assert not fence.argv_file.exists()
        assert everything_under(tmp_path) == before

    @pytest.mark.parametrize(
        ("argument", "name"),
        [("record_path", "r.json"), ("out_edges", "e.ply"), ("stats", "s.md")],
    )
    async def test_each_other_output_outside_is_refused_naming_it(
        self, tmp_path: Path, argument: str, name: str
    ) -> None:
        fence = fenced(tmp_path, tmp_path / "root")
        outside = tmp_path / "outside" / name
        before = everything_under(tmp_path)
        arguments = mesh_into(fence.root / "m.vtk", **{argument: str(outside)})
        assert_refused(await call(fence.server, "mesh", arguments), argument, fence.root)
        assert not fence.argv_file.exists()
        assert everything_under(tmp_path) == before

    async def test_stats_dash_passes(self, tmp_path: Path) -> None:
        fence = fenced(tmp_path, tmp_path / "root")
        result = await call(fence.server, "mesh", mesh_into(fence.root / "m.vtk", stats="-"))
        assert structured(result)["exit_code"] == 0, structured(result)
        assert fence.argv_file.exists()

    async def test_a_relative_out_inside_the_working_directory_passes(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        fence = fenced(tmp_path, tmp_path / "root")
        fence.root.mkdir()
        monkeypatch.chdir(fence.root)
        result = await call(fence.server, "mesh", mesh_into("m.vtk"))
        assert structured(result)["exit_code"] == 0, structured(result)
        assert fence.argv_file.exists()

    @pytest.mark.parametrize("spelling", ["link", "real"])
    async def test_a_symlinked_root_passes_through_either_name(
        self, tmp_path: Path, spelling: str
    ) -> None:
        (tmp_path / "real").mkdir()
        (tmp_path / "link").symlink_to(tmp_path / "real")
        fence = fenced(tmp_path, tmp_path / "link")
        result = await call(fence.server, "mesh", mesh_into(tmp_path / spelling / "m.vtk"))
        assert structured(result)["exit_code"] == 0, structured(result)
        assert fence.argv_file.exists()

    async def test_a_read_only_cache_outside_the_root_is_not_fenced(self, tmp_path: Path) -> None:
        """`mesh --cache` reads tiles; only output paths are fenced."""
        fence = fenced(tmp_path, tmp_path / "root")
        arguments = mesh_into(fence.root / "m.vtk", cache=str(tmp_path / "tiles"))
        result = await call(fence.server, "mesh", arguments)
        assert structured(result)["exit_code"] == 0, structured(result)
        assert fence.argv_file.exists()

    async def test_held_fetch_cache_outside_the_root_is_refused(self, tmp_path: Path) -> None:
        """The design's default on a question open to Ola: `fetch` writes tiles into
        its cache, so the cache is held inside the root. If Ola rules the cache free,
        this test flips: the child runs."""
        fence = fenced(tmp_path, tmp_path / "root")
        before = everything_under(tmp_path)
        arguments = {"source": "any", "cache": str(tmp_path / "outside" / "tiles")}
        assert_refused(await call(fence.server, "fetch", arguments), "cache", fence.root)
        assert not fence.argv_file.exists()
        assert everything_under(tmp_path) == before


# ------------------------------------------------------------ T8e


class TestT8eTheRootReported:
    async def test_sources_over_the_server_lists_the_resolved_root(self, tmp_path: Path) -> None:
        root = tmp_path / "some" / ".." / "root"
        listed = body_of(await call(mcp_server.build_server(app, out_root=root), "sources", {}))
        (line,) = listed["stdout"]
        assert json.loads(line)["out_root"] == str((tmp_path / "root").resolve())

    def test_sources_from_a_shell_lists_no_root(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("RASPUTIN_MCP_OUT_ROOT", raising=False)
        result = runner.invoke(app, ["sources"])
        assert result.exit_code == 0, result.output
        (line,) = result.stdout.splitlines()
        document = json.loads(line)
        assert "out_root" in document and document["out_root"] is None, document

    async def test_rasputin_mcp_out_root_reports_it(self, tmp_path: Path) -> None:
        given = tmp_path / "given"
        with anyio.fail_after(60):  # a server start and one child run, each under 1 s
            result = await call(production_server("--out-root", str(given)), "sources", {})
        (line,) = body_of(result)["stdout"]
        assert json.loads(line)["out_root"] == str(given.resolve())

    def test_out_parent_is_no_longer_an_option(self) -> None:
        code, output = invoke("mesh", "--dem", "d.tif", "--out", "x.vtk", "--out-parent", "x")
        assert code == USAGE, output
        assert "No such option" in output and "--out-parent" in output, output


# ------------------------------------------------------------ T10


class TestT10Instructions:
    async def test_the_vector_guidance(self) -> None:
        async with Client(server()) as client:
            instructions = client.instructions
        assert instructions is not None
        for words in ("fits", "Pass its path", "Output folders are created as needed", "out_root"):
            assert words in instructions, (words, instructions)

    def test_no_source_id_is_named(self) -> None:
        from tin_engine.sources import SOURCES

        text = mcp_server.INSTRUCTIONS
        assert "Output folders are created as needed" in text  # the new guidance is there
        for key in SOURCES:
            assert key not in text, key
        for named in ("corine", "lakes.geojson", "norway"):
            assert named not in text.lower(), named

    def test_the_route_from_an_outline_to_a_mesh_is_named(self) -> None:
        """106e, test 8 (`docs/increments/106-wednesday-cells.md`, section 8):
        the sentence names `simplify_domain`, an argument the mesh tool has."""
        from tin_engine.mcp.schema import derive

        text = " ".join(mcp_server.INSTRUCTIONS.split())
        assert "simplify_domain" in text, text
        assert "mesh reads as domain" in text, text
        (tool,) = [d for d in derive(app) if d.tool_name == "mesh"]
        assert "simplify_domain" in tool.tool().input_schema["properties"]
