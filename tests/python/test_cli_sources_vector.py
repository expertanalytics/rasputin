"""`rasputin sources` over vector sources, `register_at`, and quoted fetch hints, on the CLI.

`docs/increments/57-vector-sources.md`, sections 4.3, 4.4 and 4.8, tests T6
(the listing of fixture V and fixture L) and T9 (the review ideas:
`register_at` when the file read is `$RASPUTIN_DATA`'s, and fetch hints that
parse back with `shlex.split`).

The user's sources file is found as an agent's server finds it, in
`$XDG_CONFIG_HOME/rasputin/sources.toml` (`conftest.py` gives each test its
own empty `XDG_CONFIG_HOME`), except where T9 sets `RASPUTIN_DATA`.

Bounds: coverage is rounded by the listing to 1e-6 degrees and compared to
that step; a share is rounded to 1e-4 and compared to that step. Both were
measured on these fixtures, whose largest extent is 20 km by 40 km.
"""

from __future__ import annotations

import json
import os
import shlex
from pathlib import Path
from typing import Any

import pytest

from catchment_fixtures import EPSG
from cli_helpers import runner
from cog_fixtures import write_cache
from dataset_fixtures import (
    CACHED,
    INSIDE,
    LOCAL,
    OUTSIDE,
    PARTLY,
    SEED,
    SEED_BLOCK_ROW,
    SHIPPED,
    block_row,
    by_id,
    cached_bytes,
    cached_entry,
    lake_file,
    local_entry,
    write_coarse,
    write_local,
    write_sources,
)
from tin_engine.cli import app
from vector_fixtures import (
    CODE,
    COVER_COVERAGE,
    COVER_DESCRIPTION,
    COVER_ID,
    COVER_SHARES,
    LAKE_COVERAGE,
    LAKE_DESCRIPTION,
    LAKE_SHARE_INSIDE,
    LAKE_SHARE_PARTLY,
    LAKES_ID,
    LAND_CODES,
    LANDCOVER,
    OTHER,
    fits_of,
    fixture_v,
    vector_entry,
)

DEGREES = 1e-6
SHARE = 1e-4

DEMS = {("mesh", "--dem", None), ("catchment", "--dem", None),
        ("station-catchments", "--dem", None)}  # fmt: skip


@pytest.fixture
def config(monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)
    return Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"


@pytest.fixture
def registered(config: Path, tmp_path: Path) -> Path:
    """Fixture L, the lake file and the land-cover GeoPackage (its layer named) in the user file."""
    lakes, cover = fixture_v(tmp_path / "v")
    return write_sources(
        config,
        {
            LOCAL: local_entry(write_local(tmp_path / "tiles")),
            LAKES_ID: vector_entry(lakes, LAKE_DESCRIPTION),
            COVER_ID: vector_entry(cover, COVER_DESCRIPTION, layer=LANDCOVER),
        },
    )


def listing(*args: str) -> dict[str, Any]:
    """`rasputin sources ARGS`: exit 0 and one JSON document on one stdout line."""
    result = runner.invoke(app, ["sources", *args])
    assert result.exit_code == 0, result.output
    lines = result.stdout.splitlines()
    assert len(lines) == 1, result.stdout
    document = json.loads(lines[0])
    assert isinstance(document, dict), document
    return document


def boxed(box: tuple[float, ...]) -> dict[str, Any]:
    return listing("--bbox", *map(str, box))


# ------------------------------------------------------------ T6


class TestT6VectorListing:
    def test_the_land_cover_geopackage(self, registered: Path, tmp_path: Path) -> None:
        cover = by_id(boxed(INSIDE))[COVER_ID]
        assert cover["values"] == "features"
        assert cover["present"] is True
        assert cover["problem"] is None
        assert Path(cover["path"]) == tmp_path / "v" / "landcover.gpkg"
        assert cover["layer"] == LANDCOVER
        assert cover["geometry"] == ["Geometry"]
        assert cover["features"] == 2
        assert cover["attributes"] == [CODE]
        assert cover["coverage"] == pytest.approx(list(COVER_COVERAGE), abs=DEGREES)
        assert cover["box_covered"] == pytest.approx(COVER_SHARES[0], abs=SHARE)
        assert cover["classes"] == [
            {"attribute": CODE, "values": list(LAND_CODES), "distinct": len(LAND_CODES)}
        ]
        assert (cover["description"], cover["origin"]) == (COVER_DESCRIPTION, "user")

    def test_the_land_cover_fits(self, registered: Path) -> None:
        cover = by_id(listing())[COVER_ID]
        assert fits_of(cover) == {
            ("mesh", "--features", "corine"),
            ("mesh", "--features", "corine-water"),
            ("catchment", "--lakes", None),
            ("network-catchment", "--network", None),
        }
        spelt = {(f["tool"], f["argument"]) for f in cover["fits"]}
        assert spelt == {
            ("mesh", "features"),
            ("catchment", "lakes"),
            ("network_catchment", "network"),
        }

    def test_the_lake_file(self, registered: Path, tmp_path: Path) -> None:
        lakes = by_id(boxed(INSIDE))[LAKES_ID]
        assert (lakes["values"], lakes["present"], lakes["problem"]) == ("features", True, None)
        assert lakes["layer"] is None
        assert lakes["geometry"] == ["Polygon"]
        assert lakes["features"] == 1
        assert lakes["attributes"] == ["name"]
        assert lakes["classes"] == []
        assert lakes["coverage"] == pytest.approx(list(LAKE_COVERAGE), abs=DEGREES)
        assert lakes["box_covered"] == pytest.approx(LAKE_SHARE_INSIDE, abs=SHARE)
        assert fits_of(lakes) == {("catchment", "--lakes", None),
                                  ("station-catchments", "--lakes", None),
                                  ("network-catchment", "--network", None)}  # fmt: skip

    def test_the_raster_fields_stay_empty_for_a_vector_source(self, registered: Path) -> None:
        found = by_id(listing())
        for key in (LAKES_ID, COVER_ID):
            assert (found[key]["cell_m"], found[key]["z_unit"]) == (None, None), found[key]

    @pytest.mark.parametrize(
        ("box", "cover_share", "lake_share"),
        [(PARTLY, COVER_SHARES[1], LAKE_SHARE_PARTLY), (OUTSIDE, COVER_SHARES[2], 0.0)],
        ids=["partly", "outside"],
    )
    def test_the_share_of_other_boxes(
        self, registered: Path, box: tuple[float, ...], cover_share: float, lake_share: float
    ) -> None:
        found = by_id(boxed(box))
        assert found[COVER_ID]["box_covered"] == pytest.approx(cover_share, abs=SHARE)
        assert found[LAKES_ID]["box_covered"] == pytest.approx(lake_share, abs=SHARE)

    def test_the_raster_source_fits_dem_on_every_command(self, registered: Path) -> None:
        local = by_id(listing())[LOCAL]
        assert fits_of(local) == DEMS
        assert local["geometry"] == [] and local["features"] is None
        assert local["origin"] == "user"

    def test_the_shipped_sources_are_the_packages(self, registered: Path) -> None:
        found = by_id(listing())
        for key in SHIPPED:
            assert found[key]["origin"] == "package", key
        assert found["mapbiomas-c11"]["values"] == "classes"
        assert fits_of(found["mapbiomas-c11"]) == {("mesh", "--land-cover", None)}
        # Not cached, it still fits --dem: the run says what to fetch.
        glo30 = fits_of(found["glo30"])
        assert {("mesh", "--dem", None), ("catchment", "--dem", None)} <= glo30 <= DEMS, glo30

    def test_the_note_says_a_vector_coverage_is_its_extent(self, registered: Path) -> None:
        note = listing()["note"]
        assert isinstance(note, str) and "extent" in note, note

    def test_no_note_on_extents_without_a_vector_source_read(
        self, config: Path, tmp_path: Path
    ) -> None:
        """The package's own vector entries (`nve-lakes`, `clc2018`) are listed
        but not fetched, so no listed coverage is an extent box: no note."""
        write_sources(config, {LOCAL: local_entry(write_local(tmp_path / "tiles"))})
        note = listing()["note"] or ""
        assert "extent" not in note, note

    def test_a_geopackage_with_two_layers_and_none_named(
        self, config: Path, tmp_path: Path
    ) -> None:
        lakes, cover = fixture_v(tmp_path / "v")
        write_sources(config, {
            LAKES_ID: vector_entry(lakes, LAKE_DESCRIPTION),
            COVER_ID: vector_entry(cover, COVER_DESCRIPTION),
        })  # fmt: skip
        document = listing()
        assert document["user_file_problem"] is None
        found = by_id(document)
        problem = found[COVER_ID]["problem"]
        assert isinstance(problem, str), found[COVER_ID]
        assert LANDCOVER in problem and OTHER in problem, problem
        assert found[COVER_ID]["fits"] == []  # a file not read fits no vector option
        assert found[LAKES_ID]["problem"] is None
        assert found[LAKES_ID]["features"] == 1

    def test_a_missing_vector_file_fits_nothing(self, config: Path, tmp_path: Path) -> None:
        write_sources(config, {LAKES_ID: vector_entry(tmp_path / "gone.geojson", "x")})
        gone = by_id(listing())[LAKES_ID]
        assert gone["present"] is False
        assert gone["fits"] == []

    def test_the_cover_with_its_legend_fits_that_map_alone(
        self, config: Path, tmp_path: Path
    ) -> None:
        _lakes, cover = fixture_v(tmp_path / "v")
        write_sources(config, {
            COVER_ID: vector_entry(cover, COVER_DESCRIPTION, layer=LANDCOVER, legend="corine")
        })  # fmt: skip
        found = by_id(listing())[COVER_ID]
        assert fits_of(found) == {
            ("mesh", "--features", "corine"),
            ("catchment", "--lakes", None),
            ("network-catchment", "--network", None),
        }


# ------------------------------------------------------------ T9


class TestT9RegisterAt:
    def test_with_rasputin_data_read_register_at_is_the_xdg_file(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        root = tmp_path / "data"
        read = write_sources(root / "sources.toml", {LOCAL: local_entry(write_local(root / "t"))})
        monkeypatch.setenv("RASPUTIN_DATA", str(root))
        xdg = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
        document = listing()
        assert document["user_file"] == str(read)
        assert document["register_at"] == str(xdg)
        note = document["note"]
        assert isinstance(note, str), document
        assert "RASPUTIN_DATA" in note and str(xdg) in note, note

    def test_with_only_home_register_at_is_under_home(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        root = tmp_path / "data"
        write_sources(root / "sources.toml", {LOCAL: local_entry(write_local(root / "t"))})
        monkeypatch.setenv("RASPUTIN_DATA", str(root))
        monkeypatch.delenv("XDG_CONFIG_HOME")
        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        document = listing()
        wanted = tmp_path / "home" / ".config" / "rasputin" / "sources.toml"
        assert document["register_at"] == str(wanted)

    def test_reading_the_xdg_file_needs_no_such_sentence(self, registered: Path) -> None:
        document = listing()
        assert document["register_at"] == str(registered)
        assert "RASPUTIN_DATA" not in (document["note"] or ""), document["note"]


def hint(*args: str) -> list[str]:
    """The fetch hint of a run refused for want of cached data, split as a shell would."""
    result = runner.invoke(app, list(args), standalone_mode=False)
    message = getattr(result.exception, "message", None)
    assert isinstance(message, str) and "; run: " in message, (result.exception, result.output)
    return shlex.split(message.rpartition("; run: ")[2])


class TestT9QuotedHints:
    @pytest.fixture
    def spaced(self, config: Path, tmp_path: Path) -> Path:
        write_sources(config, {CACHED: cached_entry()})
        directory = tmp_path / "with space"
        directory.mkdir()
        return directory

    def test_mesh_with_domain_and_cache(self, spaced: Path) -> None:
        domain = lake_file(spaced / "my domain.geojson")
        cache = spaced / "my cache"
        got = hint("mesh", "--dem", CACHED, "--domain", str(domain), "--cache", str(cache),
                   "--tolerance", "5", "--out", str(spaced / "m.vtk"))  # fmt: skip
        assert got == ["rasputin", "fetch", CACHED, "--domain", str(domain), "--cache", str(cache)]

    def test_mesh_land_cover(self, spaced: Path) -> None:
        domain = lake_file(spaced / "my domain.geojson")
        cache = spaced / "my cache"
        dem = write_coarse(spaced / "dem" / "bowl.tif")
        got = hint("mesh", "--dem", str(dem), "--domain", str(domain), "--tolerance", "5",
                   "--land-cover", "mapbiomas-c11", "--cache", str(cache),
                   "--out", str(spaced / "m.vtk"))  # fmt: skip
        assert got[:3] == ["rasputin", "fetch", "mapbiomas-c11"], got
        assert got[got.index("--domain") + 1] == str(domain), got
        assert got[got.index("--cache") + 1] == str(cache), got

    def test_catchment_with_a_cache_holding_a_space(self, spaced: Path) -> None:
        cache = spaced / "my cache"
        write_cache(cache, CACHED, {"dem": cached_bytes()}, skip={"dem": block_row(SEED_BLOCK_ROW)})
        got = hint("catchment", "--dem", CACHED, "--cache", str(cache),
                   "--seed", repr(SEED[0]), repr(SEED[1]), "--seed-crs", EPSG,
                   "--out", str(spaced / "c.geojson"))  # fmt: skip
        assert got[:4] == ["rasputin", "fetch", CACHED, "--bbox"], got
        assert got[-2:] == ["--cache", str(cache)], got
        assert len(got) == 10, got  # four numbers between --bbox and --cache


# ------------------------------------------------- the AR50 legend, test 14


class TestAr50Listing:
    """Test 14 of `docs/increments/73-gml-vector-source.md`: fixture A, converted
    with `gml_to_geopackage` and registered with `legend = "ar50"`, lists
    `arealtype`'s classes (fixture A holds 30 once and 60 five times) and fits
    `--features` with `ar50` alone."""

    @pytest.fixture
    def ar50_entry(self, config: Path, tmp_path: Path) -> Path:
        import importlib
        from contextlib import closing

        from gml32_fixtures import AR50_SAMPLE
        from tin_engine.io.repository import create_sqlite

        gml_import = importlib.import_module("tin_engine.gml_import")
        target = tmp_path / "ar50.gpkg"
        with AR50_SAMPLE.open("rb") as stream, closing(create_sqlite(target)) as conn:
            gml_import.gml_to_geopackage(stream, conn)
        write_sources(config, {"ar50-sample": vector_entry(target, "AR50", legend="ar50")})
        return target

    def test_14_the_classes_and_fits(self, ar50_entry: Path) -> None:
        found = by_id(listing())["ar50-sample"]
        assert found["problem"] is None, found
        assert found["classes"] == [
            {"attribute": "arealtype", "values": ["30", "60"], "distinct": 2}
        ]
        assert fits_of(found) == {
            ("mesh", "--features", "ar50"),
            ("catchment", "--lakes", None),
            ("network-catchment", "--network", None),
        }
