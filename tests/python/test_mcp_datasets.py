"""An agent finds the local elevation data, chooses it, and passes it by name, over MCP.

`docs/increments/55-dataset-discovery.md`, sections 4.3 to 4.6, tests T3
(discovery), T4 (choice from the listing), T5 (selection then use), T8 (the
instructions) and T9 (the server started with no environment but `HOME` and
`PATH`, the condition Codex starts it in).

The server is the production one: `build_server(cli.app)` with the real
``python -m tin_engine`` child in process (T3 to T5, T8), and ``python -m
tin_engine mcp`` over stdio (T9). No test sets `RASPUTIN_DATA` or
`XDG_CONFIG_HOME`: the user's sources file is found under `HOME`, a temporary
directory, as it is found for a client whose configuration names only the
command. The child inherits the server's environment, so the `monkeypatch`
of `HOME` reaches it.

Bounds: every call is one child run, measured under 0.1 s of work on these
fixtures plus about 0.5 s of interpreter start; the stdio test's 60 s is a
hang guard, not a timing claim.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("mcp")

import anyio
from mcp.client import Client
from mcp.client.stdio import StdioServerParameters

from catchment_fixtures import EPSG
from dataset_fixtures import (
    COARSE,
    COVERAGE,
    CREDIT,
    DESCRIPTION,
    INSIDE,
    LICENCE,
    LOCAL,
    NODES_WITH_LAKE,
    OUTSIDE,
    PARTLY,
    SEED,
    SHARE_PARTLY,
    SHIPPED,
    by_id,
    lake_file,
    local_entry,
    write_coarse,
    write_local,
    write_sources,
)
from mcp_driver import call, structured, text_block
from tin_engine.cli import app
from tin_engine.mcp import server as mcp_server

#: Degrees: the coverage is rounded to 1e-6 by the listing, so the measured
#: box holds to within one rounding step.
DEGREES = 1e-6


@pytest.fixture
def home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A home whose `~/.config/rasputin/sources.toml` registers fixture L and its 200 m copy;
    `RASPUTIN_DATA` and `XDG_CONFIG_HOME` unset, as in a client that sets no `env`."""
    root = tmp_path / "home"
    write_sources(
        root / ".config" / "rasputin" / "sources.toml",
        {
            LOCAL: local_entry(write_local(tmp_path / "data" / "bowl")),
            COARSE: local_entry(write_coarse(tmp_path / "data" / "coarse" / "bowl-200.tif")),
        },
    )
    monkeypatch.setenv("HOME", str(root))
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    return root


def server(out_root: Path | None = None) -> Any:
    """The production server; a test that writes gives its own directory as the output root."""
    return mcp_server.build_server(app, out_root=out_root)


def body_of(result: Any) -> dict[str, Any]:
    assert not result.is_error, result.content
    body = structured(result)
    assert body["exit_code"] == 0, body
    return body


def document_of(result: Any) -> dict[str, Any]:
    """The listing: one JSON document, the one line of the tool result's stdout."""
    stdout = body_of(result)["stdout"]
    assert len(stdout) == 1, stdout
    document = json.loads(stdout[0])
    assert isinstance(document, dict), document
    return document


async def listed(**arguments: Any) -> dict[str, Any]:
    return document_of(await call(server(), "sources", arguments))


# ------------------------------------------------------------ T3


class TestT3Discovery:
    async def test_the_local_entry_and_its_figures(self, home: Path) -> None:
        local = by_id(await listed())[LOCAL]
        assert local["kind"] == "local"
        assert local["values"] == "heights"
        assert local["present"] is True
        assert Path(local["path"]) == home.parent / "data" / "bowl"
        assert local["tiles"] == 4
        assert local["crs"] == EPSG
        assert "UTM zone 33N" in local["crs_name"], local["crs_name"]
        assert local["cell_m"] == [100.0, 100.0]
        assert (local["z_unit"], local["z_unit_from"]) == ("metre", "catalogue")
        assert local["coverage"] == pytest.approx(list(COVERAGE), abs=DEGREES)
        assert local["box_covered"] is None  # no bbox given
        assert (local["description"], local["credit"]) == (DESCRIPTION, CREDIT)
        assert local["licence_note"] == LICENCE
        assert local["problem"] is None

    async def test_the_places_read(self, home: Path) -> None:
        document = await listed()
        expected = home / ".config" / "rasputin" / "sources.toml"
        assert Path(document["user_file"]) == expected
        assert [Path(p) for p in document["looked_at"]] == [expected]
        assert Path(document["register_at"]) == expected
        assert document["user_file_problem"] is None

    async def test_the_shipped_sources_are_listed_uncached_with_descriptions(
        self, home: Path
    ) -> None:
        found = by_id(await listed())
        for key in SHIPPED:
            entry = found[key]
            assert entry["present"] is False, entry
            assert entry["description"].strip(), entry
            assert entry["coverage"] is None and entry["cell_m"] is None, entry

    async def test_package_sources_come_first_then_the_users_in_file_order(
        self, home: Path
    ) -> None:
        ids = [s["id"] for s in (await listed())["sources"]]
        assert ids[-2:] == [LOCAL, COARSE], ids
        assert set(ids[:-2]) == set(SHIPPED), ids

    def test_the_tool_is_read_only(self) -> None:
        tools = {d.tool_name: d.tool() for d in mcp_server.derive(app)}
        annotations = tools["sources"].annotations
        assert annotations is not None
        assert annotations.read_only_hint is True
        assert annotations.destructive_hint is False


# ------------------------------------------------------------ T4


class TestT4Choice:
    @pytest.mark.parametrize(
        ("box", "share"),
        [(INSIDE, 1.0), (PARTLY, SHARE_PARTLY), (OUTSIDE, 0.0)],
        ids=["inside", "partly", "outside"],
    )
    async def test_the_share_of_the_box_covered(
        self, home: Path, box: tuple[float, ...], share: float
    ) -> None:
        # Rounded to 1e-4 by the listing; the measured share is stated to that step.
        found = by_id(await listed(bbox=list(box)))
        assert found[LOCAL]["box_covered"] == pytest.approx(share, abs=1e-4)
        for key in SHIPPED:
            assert found[key]["box_covered"] is None, found[key]

    async def test_two_local_sources_told_apart_by_cell_size_alone(self, home: Path) -> None:
        found = by_id(await listed(bbox=list(INSIDE)))
        assert found[LOCAL]["cell_m"] == [100.0, 100.0]
        assert found[COARSE]["cell_m"] == [200.0, 200.0]
        assert found[COARSE]["tiles"] == 1
        assert found[COARSE]["box_covered"] == 1.0
        assert found[COARSE]["description"] == found[LOCAL]["description"]

    async def test_a_point_box(self, home: Path) -> None:
        """A point is a box with west equal to east and south equal to north."""
        lon, lat = INSIDE[0], INSIDE[1]
        found = by_id(await listed(bbox=[lon, lat, lon, lat]))
        assert found[LOCAL]["box_covered"] == 1.0
        far = by_id(await listed(bbox=[OUTSIDE[0], OUTSIDE[1], OUTSIDE[0], OUTSIDE[1]]))
        assert far[LOCAL]["box_covered"] == 0.0


# ------------------------------------------------------------ T5


def chosen(document: dict[str, Any]) -> str:
    """What an agent reads off the listing: of the sources covering the whole
    box, the one with the smallest cell."""
    covering = [s for s in document["sources"] if s["box_covered"] == 1.0]
    assert covering, document
    return str(min(covering, key=lambda s: s["cell_m"][1])["id"])


class TestT5SelectionThenUse:
    async def test_choose_then_catchment_then_mesh(self, home: Path, tmp_path: Path) -> None:
        key = chosen(await listed(bbox=list(INSIDE)))
        assert key == LOCAL
        lakes = lake_file(tmp_path / "lakes.geojson")
        outline = tmp_path / "catchment.geojson"
        catchment = await call(
            server(tmp_path),
            "catchment",
            {"out": str(outline), "dem": [key], "seed": list(SEED), "seed_crs": EPSG,
             "lakes": str(lakes)},
        )  # fmt: skip
        body_of(catchment)
        (feature,) = json.loads(outline.read_text())["features"]
        assert feature["properties"]["nodes"] == NODES_WITH_LAKE
        out = tmp_path / "mesh.vtk"
        mesh = await call(
            server(tmp_path), "mesh", {"out": str(out), "dem": [key], "domain": str(outline),
                               "tolerance": 5.0},
        )  # fmt: skip
        record = body_of(mesh)["record"]
        assert isinstance(record, dict), record
        assert record["dem_credit"] == CREDIT
        assert record["licence_note"] == LICENCE
        assert out.is_file()

    async def test_fetch_of_a_local_source_is_refused(self, home: Path) -> None:
        result = await call(server(), "fetch", {"source": LOCAL, "bbox": list(INSIDE)})
        assert result.is_error, result.content
        assert "needs no fetch" in text_block(result), text_block(result)


# ------------------------------------------------------------ T8


class TestT8Instructions:
    async def test_the_agent_is_sent_to_the_listing_not_to_a_download(self) -> None:
        async with Client(server()) as client:
            instructions = client.instructions
        assert instructions is not None
        assert "Call sources first" in instructions, instructions
        assert "A local source needs no fetch" in instructions, instructions
        assert "must be fetched first" not in instructions, instructions

    def test_no_source_country_or_resolution_is_named(self) -> None:
        """General guidance only: the agent reasons on the listing's figures."""
        words = mcp_server.INSTRUCTIONS.lower()
        for named in ("norway", "dtm10", "glo30", "glo-30", "anadem", "10 m", "30 m"):
            assert named not in words, named


# ------------------------------------------------------------ T9


async def test_a_server_with_only_home_and_path_finds_the_local_source(home: Path) -> None:
    """The Codex condition: the client's configuration names the command and
    nothing else, so the server sees the SDK's default environment, here with
    `HOME` the test's home; no `RASPUTIN_DATA`, no `XDG_CONFIG_HOME`."""
    parameters = StdioServerParameters(
        command=sys.executable,
        args=["-m", "tin_engine", "mcp"],
        env={"HOME": str(home), "PATH": os.environ["PATH"]},
    )
    with anyio.fail_after(60):
        async with Client(parameters) as client:
            result = await client.call_tool("sources", {})
    document = document_of(result)
    local = by_id(document)[LOCAL]
    assert (local["present"], local["tiles"], local["crs"]) == (True, 4, EPSG)
    assert Path(document["user_file"]) == home / ".config" / "rasputin" / "sources.toml"
