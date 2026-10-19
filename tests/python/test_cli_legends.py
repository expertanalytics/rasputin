"""Legends as data through `rasputin mesh` and `rasputin palette` (41).

`docs/increments/41-class-maps-as-data.md`, section 7: the CLI halves of T3
(the `rasputin palette corine` bytes), T5 (name or path, for `mesh` and
`palette`), T7 (the derived clean-up distances, the overrides and the record),
T8 (text codes on triangles), T9 (vocabulary extras in the `.vtk`), T10
(MapBiomas: `--land-cover`'s `.gpkg` names and colours, the `.vtk` unchanged,
the notice) and T11 (the licence in the record). The scene is the CLI suites'
`rough_dem(16)`: 17 by 21 nodes, cells 10 m by 5 m, so the DEM cell `h` of
section 5.2 is 10 m, the larger spacing; the smaller, 5 m, would give other
numbers in T7's fourth case, which is how that case tells them apart.

Numeric bounds: the record's distances are compared exactly. Each is one
division of two small integers-as-floats (section 5.1), so the expected
literals (0.4, 0.1, ...) are the same doubles; checked up to a 1 km cell in
`test_legends.TestDerivedCleanup`.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The record field for a map's licence is `features_licence`; absent when
  every map's licence is empty.
- `rasputin palette PATH` writes the file's `preset` as the preset's `Name`.
- The three help texts end "Default: from the legend's positional
  accuracy and the DEM cell (CORINE: N)." with N 1, 50 and 5.
- A `--features-legend` file error, an unknown name, an extra repeating a
  default name, and more than 22 extras in one run are usage errors (exit 2).
- A `--land-cover` catalogue entry's `legend` is honoured from the entry
  in `SOURCES`, so a test entry carrying `legend="mapbiomas-c11"` gets
  MapBiomas's names and colours.

RED at the commit that adds this file: no `tin_engine.legends` and no
`data/legends/`; `--features-legend` takes only the four built-in names
(a path is "unknown map"); the clean-up defaults are fixed at 1, 50 and 5;
`palette` knows only `corine`; `--land-cover` writes no class names; no
extra names, no text codes, no licence row.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any

import pytest
import typer
from shapely.geometry import LineString

import tin_engine.cli as cli
from class_fixtures import PATH_TEMPLATE, per_year_entry
from cli_helpers import SQUARE, USAGE, geojson, invoke, rough_dem, runner
from feature_fixtures import Feat, write_geojson
from fraction_fixtures import write_class_cache, write_dem, write_domain
from legend_fixtures import basenames, coded_map, corine_like, uncoded_map, write
from legend_snapshot import PRESET_BYTES
from recordread import file_field
from test_cli_features_cleanup import Spy, run, switches
from test_cli_fetch import install
from test_cli_mesh_cache import no_network
from test_cli_mesh_features import edges_at
from test_cli_mesh_landcover import EAST, WEST, coded, rel
from tin_engine.cli import app
from tin_engine.features import DEFAULT_VOCABULARY, EdgeProperty, EdgeVocabulary
from vtkread import read_vtk

SHIPPED = ("corine", "corine-water", "clc18_kode", "mapbiomas-c11", "ar50")
LAND_COVER, WATER = DEFAULT_VOCABULARY.mask("land_cover"), DEFAULT_VOCABULARY.mask("water")
CORINE_DISTANCES = (1.0, 50.0, 5.0)
OPTIONS = ("features_repair", "features_tolerance", "features_outline_snap")
DERIVED = "Default: from the legend's positional accuracy and the DEM cell (CORINE: {})."
LEVEE = LineString([rel(30.1, -40.3), rel(170.2, -40.3)])
ROAD = LineString([rel(30.1, -15.2), rel(170.2, -15.2)])

bumpy = rough_dem(16)


@pytest.fixture
def square(tmp_path: Path) -> Path:
    return geojson(tmp_path / "square.geojson", SQUARE)


@pytest.fixture
def corine(tmp_path: Path) -> Path:
    """`Code_18` forest (311) west, water (512) east, as the clean-up suite's."""
    return write_geojson(
        tmp_path / "corine.geojson", [coded(1, WEST, "311"), coded(2, EAST, "512")]
    )


@pytest.fixture
def gallery(tmp_path: Path) -> Path:
    from test_cli_mesh_features import GALLERY

    return write_geojson(tmp_path / "gallery.geojson", GALLERY)


@pytest.fixture
def levees(tmp_path: Path) -> Path:
    return write_geojson(
        tmp_path / "levees.geojson",
        [Feat("l", LEVEE, {"kind": "levee"}), Feat("r", ROAD, {"kind": "road"})],
    )


def recorded(record: dict[str, Any]) -> tuple[float, float, float]:
    return (
        record["features_repair_m"],
        record["features_tolerance_m"],
        record["features_outline_snap_m"],
    )


def mesh_params() -> dict[str, Any]:
    return {p.name: p for p in typer.main.get_command(cli.app).commands["mesh"].params}


# ------------------------------------------------------------- T3, T5: palette


class TestPalette:
    def test_t3_palette_corine_writes_masters_bytes(self) -> None:
        result = runner.invoke(app, ["palette", "corine"])
        assert result.exit_code == 0, result.output
        assert result.stdout_bytes == PRESET_BYTES

    def test_t3_palette_corine_out_writes_masters_bytes(self, tmp_path: Path) -> None:
        target = tmp_path / "p.json"
        code, said = invoke("palette", "corine", "--out", str(target))
        output = basenames(said)
        assert code == 0, output
        assert target.read_bytes() == PRESET_BYTES

    def test_t5_palette_of_a_path(self, tmp_path: Path) -> None:
        path = write(tmp_path / "mine.toml", coded_map())
        result = runner.invoke(app, ["palette", str(path)])
        assert result.exit_code == 0, result.output
        (preset,) = json.loads(result.stdout)
        assert preset["Name"] == "rasputin mine"
        assert preset["Annotations"] == [
            *("0", "0 No polygon; constraint lines"),
            *("1", "1 Forest"),
            *("2", "2 Lake"),
        ]

    def test_t5_palette_of_the_shipped_mapbiomas_map(self) -> None:
        result = runner.invoke(app, ["palette", "mapbiomas-c11"])
        assert result.exit_code == 0, result.output
        (preset,) = json.loads(result.stdout)
        maps = importlib.import_module("tin_engine.legends")
        assert preset["Name"] == maps.load("mapbiomas-c11").preset
        assert "33 River, Lake and Ocean" in preset["Annotations"]
        assert len(preset["Annotations"]) == 2 * 34

    def test_t5_palette_refuses_an_uncoded_map(self, tmp_path: Path) -> None:
        path = write(tmp_path / "levees.toml", uncoded_map())
        code, said = invoke("palette", str(path))
        output = basenames(said)
        assert code == USAGE, output
        assert "has no classes to colour" in output
        assert "levees.toml" in output

    def test_t5_palette_refuses_an_unknown_name_listing_the_shipped(self) -> None:
        code, said = invoke("palette", "corin")
        output = basenames(said)
        assert code == USAGE, output
        for name in SHIPPED:
            assert name in output, name


# ------------------------------------------------------------- T5: mesh


class TestFeaturesMapOption:
    def test_t5_the_help_lists_the_shipped_names(self) -> None:
        text = mesh_params()["features_legend"].help or ""
        for name in (*SHIPPED, "property"):
            assert name in text, name

    def test_t5_an_unknown_name_is_refused_listing_the_shipped(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        code, said = invoke(
            "mesh",
            *("--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"),
            *("--out", str(tmp_path / "x.vtk"), "--features", str(corine)),
            *("--features-legend", "corin"),
        )
        output = basenames(said)
        assert code == USAGE, output
        for name in ("corin", *SHIPPED):
            assert name in output, name

    def test_t5_a_missing_path_is_refused_naming_it_before_the_dem_is_read(
        self, tmp_path: Path, square: Path, corine: Path
    ) -> None:
        """Section 4, "Maps load once, early": every file error comes before
        any reading of the DEM, so a DEM that is not there is never reached."""
        code, said = invoke(
            "mesh",
            *("--dem", str(tmp_path / "absent.tif"), "--domain", str(square)),
            *("--tolerance", "1", "--out", str(tmp_path / "x.vtk")),
            *("--features", str(corine), "--features-legend", str(tmp_path / "nope.toml")),
        )
        output = basenames(said)
        assert code == USAGE, output
        assert "nope.toml" in output
        assert "absent.tif" not in output

    def test_t5_a_bad_file_is_refused_naming_it(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        path = write(tmp_path / "bad.toml", corine_like("bad", "40").replace('"#4f8f3f"', '"x"'))
        code, said = invoke(
            "mesh",
            *("--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"),
            *("--out", str(tmp_path / "x.vtk"), "--features", str(corine)),
            *("--features-legend", str(path)),
        )
        output = basenames(said)
        assert code == USAGE, output
        assert "bad.toml" in output and "colour" in output


# ------------------------------------------------------------- T7


class TestDerivedDistances:
    """T7 on `rough_dem` (h = 10 m): the record carries the values used, and
    `open_features` is handed the same."""

    def test_t7_corine_with_no_option_is_todays_1_50_5(
        self,
        tmp_path: Path,
        bumpy: Path,
        square: Path,
        corine: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        spy = Spy(monkeypatch)
        _, record, _ = run(
            tmp_path, bumpy, square, "--features", str(corine), "--features-legend", "corine"
        )
        assert recorded(record) == CORINE_DISTANCES
        (request,) = spy.requests
        assert switches(request) == (1.0, True, 50.0, 5.0)

    def test_t7_one_override_replaces_only_its_own(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        _, record, _ = run(
            tmp_path,
            bumpy,
            square,
            *("--features", str(corine), "--features-legend", "corine"),
            *("--features-tolerance", "20"),
        )
        assert recorded(record) == (1.0, 20.0, 5.0)

    @pytest.mark.parametrize(
        ("accuracy", "want"),
        [
            ("40", (0.4, 20.0, 2.0)),
            # The cell floor: s = max(4, 10) with h the larger spacing, 10 m;
            # the smaller (5 m) would give 0.05, 2.5, 0.25.
            ("4", (0.1, 5.0, 0.5)),
        ],
    )
    def test_t7_a_user_maps_accuracy(
        self,
        tmp_path: Path,
        bumpy: Path,
        square: Path,
        corine: Path,
        monkeypatch: pytest.MonkeyPatch,
        accuracy: str,
        want: tuple[float, float, float],
    ) -> None:
        spy = Spy(monkeypatch)
        path = write(tmp_path / "mine.toml", corine_like("mine", accuracy))
        report, record, _ = run(
            tmp_path, bumpy, square, "--features", str(corine), "--features-legend", str(path)
        )
        assert recorded(record) == want
        (request,) = spy.requests
        assert (request.repair_m, request.tolerance_m, request.outline_snap_m) == want
        assert "legend mine" in record["features"]
        assert report

    def test_t7_a_user_map_with_one_override(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        path = write(tmp_path / "mine.toml", corine_like("mine", "40"))
        _, record, _ = run(
            tmp_path,
            bumpy,
            square,
            *("--features", str(corine), "--features-legend", str(path)),
            *("--features-repair", "0.7"),
        )
        assert recorded(record) == (0.7, 20.0, 2.0)

    def test_t7_two_coded_sources_take_the_larger_accuracy(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        path = write(tmp_path / "mine.toml", corine_like("mine", "40"))
        other = write_geojson(tmp_path / "other.geojson", [coded(3, WEST, "311")])
        _, record, _ = run(
            tmp_path,
            bumpy,
            square,
            *("--features", str(corine), "--features-legend", str(path)),
            *("--features", str(other), "--features-legend", "corine"),
        )
        assert recorded(record) == CORINE_DISTANCES

    def test_t7_an_uncoded_only_run_records_1_50_5(
        self, tmp_path: Path, bumpy: Path, square: Path, gallery: Path
    ) -> None:
        _, record, _ = run(tmp_path, bumpy, square, "--features", str(gallery))
        assert recorded(record) == CORINE_DISTANCES

    @pytest.mark.parametrize(
        ("param", "corine_value"),
        [("features_repair", "1"), ("features_tolerance", "50"), ("features_outline_snap", "5")],
    )
    def test_t7_the_help_says_derived(self, param: str, corine_value: str) -> None:
        help_text = mesh_params()[param].help or ""
        assert help_text.endswith(DERIVED.format(corine_value)), help_text

    @pytest.mark.parametrize("param", OPTIONS)
    def test_t7_the_options_default_to_none(self, param: str) -> None:
        assert mesh_params()[param].default is None


# ------------------------------------------------------------- T8


class TestTextCodesOnTriangles:
    def test_t8_triangles_carry_1_and_2_and_the_lake_border_water(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        cover = write_geojson(
            tmp_path / "cover.geojson",
            [Feat("w", WEST, {"cover": "Forest"}), Feat("e", EAST, {"cover": "Lake"})],
        )
        path = write(tmp_path / "cover.toml", coded_map(codes="Forest = 1\nLake = 2\n"))
        run(tmp_path, bumpy, square, "--features", str(cover), "--features-legend", str(path))
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        codes = vtk.cell_array("land_cover_code").values
        lines = len(vtk.lines)
        assert {1, 2} <= {int(c) for c in codes[lines:]}
        assert set(edges_at(vtk, rel(160.9, -60.2))) == {LAND_COVER | WATER}
        assert set(edges_at(vtk, rel(40.3, -60.2))) == {LAND_COVER}


# ------------------------------------------------------------- T9


class TestVocabularyExtrasInTheMesh:
    def test_t9_levee_on_bit_10_in_the_table_and_on_its_line(
        self, tmp_path: Path, bumpy: Path, square: Path, levees: Path
    ) -> None:
        path = write(tmp_path / "levees.toml", uncoded_map())
        run(tmp_path, bumpy, square, "--features", str(levees), "--features-legend", str(path))
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        bits = [int(b) for b in vtk.field_data["feature_bits"].values]
        names = [str(n) for n in vtk.field_data["feature_names"].values]
        assert (10, "levee") in set(zip(bits, names, strict=True))
        assert 9 not in bits
        run_vocabulary = EdgeVocabulary(
            properties=(*DEFAULT_VOCABULARY.properties, EdgeProperty(name="levee", bit=10))
        )
        assert file_field(vtk, "feature_vocabulary") == run_vocabulary.fingerprint()
        assert set(edges_at(vtk, rel(30.1, -40.3))) == {1 << 10}
        assert set(edges_at(vtk, rel(30.1, -15.2))) == {DEFAULT_VOCABULARY.mask("road")}

    def test_t9_without_extras_the_fingerprint_is_the_defaults(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        path = write(tmp_path / "mine.toml", corine_like("mine", "40"))
        run(tmp_path, bumpy, square, "--features", str(corine), "--features-legend", str(path))
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        assert file_field(vtk, "feature_vocabulary") == DEFAULT_VOCABULARY.fingerprint()

    def test_t9_an_extra_named_river_is_refused(
        self, tmp_path: Path, bumpy: Path, square: Path, levees: Path
    ) -> None:
        path = write(tmp_path / "rivers.toml", uncoded_map(extra='"levee", "river"'))
        code, said = invoke(
            "mesh",
            *("--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"),
            *("--out", str(tmp_path / "x.vtk"), "--features", str(levees)),
            *("--features-legend", str(path)),
        )
        output = basenames(said)
        assert code == USAGE, output
        assert "rivers.toml" in output and "river" in output.replace("rivers.toml", "")
        assert "unknown map" not in output

    def test_t9_twenty_three_extras_over_two_maps_are_refused(
        self, tmp_path: Path, bumpy: Path, square: Path, levees: Path
    ) -> None:
        """Section 6: the run's vocabulary is the default plus every map's
        extras in order of first appearance (``levee`` once); bits 10 to 31
        hold 22, and these are 23. Each map alone fits."""
        first = ", ".join(['"levee"', *(f'"a{i}"' for i in range(11))])
        second = ", ".join(['"levee"', *(f'"b{i}"' for i in range(11))])
        one = write(tmp_path / "one.toml", uncoded_map("one", extra=first))
        two = write(tmp_path / "two.toml", uncoded_map("two", extra=second))
        code, said = invoke(
            "mesh",
            *("--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"),
            *("--out", str(tmp_path / "x.vtk")),
            *("--features", str(levees), "--features-legend", str(one)),
            *("--features", str(levees), "--features-legend", str(two)),
        )
        output = basenames(said)
        assert code == USAGE, output
        assert "22" in output


# ------------------------------------------------------------- T10, T11


class TestMapBiomasThroughTheCli:
    @pytest.fixture
    def scene(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
        import socket

        monkeypatch.setattr(socket, "socket", no_network)
        monkeypatch.setattr(socket, "create_connection", no_network)
        return {
            "cache": write_class_cache(tmp_path / "cache"),
            "dem": write_dem(tmp_path / "dem.tif"),
            "domain": write_domain(tmp_path / "domain.geojson"),
        }

    def land_cover_mesh(
        self,
        tmp_path: Path,
        scene: dict[str, Path],
        monkeypatch: pytest.MonkeyPatch,
        out: str,
        **entry: Any,
    ) -> Path:
        """One `--land-cover` run of the test class source, its entry carrying
        ``entry`` (`legend`, say), written to ``out``."""
        sources = importlib.import_module("tin_engine.sources")
        template = "https://example.invalid/" + PATH_TEMPLATE
        key = install(sources, monkeypatch, per_year_entry(sources, template, **entry)).id
        target = tmp_path / out
        code, said = invoke(
            "mesh",
            *("--dem", str(scene["dem"]), "--domain", str(scene["domain"]), "--tolerance", "3"),
            *("--cache", str(scene["cache"]), "--land-cover", key, "--out", str(target)),
        )
        output = basenames(said)
        assert code == 0, output
        return target

    def test_t10_the_gpkg_triangles_carry_mapbiomas_names_and_colours(
        self, tmp_path: Path, scene: dict[str, Path], monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import sqlite3
        from contextlib import closing

        target = self.land_cover_mesh(
            tmp_path, scene, monkeypatch, "named.gpkg", legend="mapbiomas-c11"
        )
        palette = importlib.import_module("tin_engine.legends").load("mapbiomas-c11").palette
        with closing(sqlite3.connect(target)) as conn:
            rows = conn.execute("SELECT code, land_cover, colour FROM triangles").fetchall()
        assert len({r[0] for r in rows}) > 2
        for code, name, colour in rows:
            assert (name, colour) == palette[code], code

    def test_t10_the_vtk_carries_no_new_field(
        self, tmp_path: Path, scene: dict[str, Path], monkeypatch: pytest.MonkeyPatch
    ) -> None:
        plain = self.land_cover_mesh(tmp_path, scene, monkeypatch, "plain.vtk")
        named = self.land_cover_mesh(
            tmp_path, scene, monkeypatch, "named.vtk", legend="mapbiomas-c11"
        )
        before, after = read_vtk(plain.read_bytes()), read_vtk(named.read_bytes())
        assert sorted(after.field_data) == sorted(before.field_data)

    def test_t10_a_features_run_writes_the_mapbiomas_notice(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        cover = write_geojson(
            tmp_path / "mb.geojson",
            [Feat(1, WEST, {"class": 3}), Feat(2, EAST, {"class": 33})],
        )
        _, record, _ = run(
            tmp_path, bumpy, square, "--features", str(cover), "--features-legend", "mapbiomas-c11"
        )
        notice = importlib.import_module("tin_engine.legends").load("mapbiomas-c11").notice
        assert record["features_notice"] == notice
        assert "MapBiomas Project - Collection 11" in notice
        # MapBiomas at 30 m over a 10 m cell: s = 30.
        assert recorded(record) == (0.3, 15.0, 1.5)


class TestLicenceInTheRecord:
    def test_t11_a_licence_is_recorded(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        text = corine_like("mine", "40").replace(
            'unlisted = "land"\n', 'unlisted = "land"\nlicence = "CC BY 4.0"\n'
        )
        path = write(tmp_path / "mine.toml", text)
        _, record, _ = run(
            tmp_path, bumpy, square, "--features", str(corine), "--features-legend", str(path)
        )
        assert record["features_licence"] == "CC BY 4.0"

    def test_t11_corine_writes_no_licence(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        _, record, _ = run(
            tmp_path, bumpy, square, "--features", str(corine), "--features-legend", "corine"
        )
        assert "features_licence" not in record


# ------------------------------------------------- the AR50 legend, tests 12, 13


class TestAr50ThroughTheCli:
    """Tests 12 and 13 of `docs/increments/73-gml-vector-source.md`: `rasputin
    palette ar50`, and a mesh coded by `ar50`: forest (30) west, fresh water
    (81) east. AR50's 2 m accuracy over the 10 m cell: s = max(2, 10) = 10."""

    @pytest.fixture
    def ar50(self, tmp_path: Path) -> Path:
        return write_geojson(
            tmp_path / "ar50.geojson",
            [Feat(1, WEST, {"arealtype": 30}), Feat(2, EAST, {"arealtype": 81})],
        )

    def test_12_palette_ar50(self) -> None:
        result = runner.invoke(app, ["palette", "ar50"])
        assert result.exit_code == 0, result.output
        (preset,) = json.loads(result.stdout)
        assert preset["Name"] == "rasputin AR50 NIBIO"
        assert len(preset["Annotations"]) == 2 * 10
        assert "81 Freshwater" in preset["Annotations"]

    def test_13_the_triangles_borders_and_record(
        self, tmp_path: Path, bumpy: Path, square: Path, ar50: Path
    ) -> None:
        _, record, _ = run(
            tmp_path, bumpy, square, "--features", str(ar50), "--features-legend", "ar50"
        )
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        codes = vtk.cell_array("land_cover_code").values
        assert {30, 81} <= {int(c) for c in codes[len(vtk.lines) :]}
        assert set(edges_at(vtk, rel(160.9, -60.2))) == {LAND_COVER | WATER}
        assert set(edges_at(vtk, rel(40.3, -60.2))) == {LAND_COVER}
        legend = importlib.import_module("tin_engine.legends").load("ar50")
        assert record["features_notice"] == legend.notice
        assert record["features_licence"] == legend.licence
        assert recorded(record) == (0.1, 5.0, 0.5)

    def test_13_the_gpkg_triangles_carry_ar50_names_and_colours(
        self, tmp_path: Path, bumpy: Path, square: Path, ar50: Path
    ) -> None:
        import sqlite3
        from contextlib import closing

        target = tmp_path / "x.gpkg"
        code, said = invoke(
            "mesh",
            *("--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"),
            *("--out", str(target), "--features", str(ar50), "--features-legend", "ar50"),
        )
        assert code == 0, basenames(said)
        palette = importlib.import_module("tin_engine.legends").load("ar50").palette
        with closing(sqlite3.connect(target)) as conn:
            rows = conn.execute("SELECT code, land_cover, colour FROM triangles").fetchall()
        assert {30, 81} <= {r[0] for r in rows}
        for found, name, colour in rows:
            assert (name, colour) == palette[found], found
