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

- The record field for a map's licence is `land_cover_licence`; absent when
  every map's licence is empty.
- `rasputin palette PATH` writes the file's `preset` as the preset's `Name`.
- The three help texts end "Default: from the land-cover class file's positional
  accuracy and the DEM cell (CORINE: N)." with N 1, 50 and 5.
- A `--land-cover-classes` file error and an unknown name are usage errors
  (exit 2).
- A `--land-cover` catalogue entry's `classes` is honoured from the entry
  in `SOURCES`, so a test entry carrying `classes="mapbiomas-c11"` gets
  MapBiomas's names and colours.
- `docs/increments/111-one-land-cover-map.md`: a class-polygon run is
  `--land-cover FILE --land-cover-classes NAME --edges all`; the record's
  `edges` entry names the class file (`classes NAME`), and `WORDING["edges"]`
  is "Each source of edges: file, layer, land-cover class file or line kind,
  features, lines, vertices" (the design gives the row, not this wording).

RED at the commit that adds this file: no `tin_engine.legends` and no
`data/legends/`; `--features-rules` takes only the four built-in names
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

import tin_engine.cli as cli
from class_fixtures import PATH_TEMPLATE, per_year_entry
from cli_helpers import SQUARE, USAGE, geojson, invoke, rough_dem, runner
from feature_fixtures import Feat, write_geojson
from fraction_fixtures import write_class_cache, write_dem, write_domain
from legend_fixtures import basenames, coded_map, corine_like, legend_text, write
from legend_snapshot import PRESET_BYTES
from recordread import file_field
from test_cli_features_cleanup import Spy, run, switches
from test_cli_fetch import install
from test_cli_mesh_cache import no_network
from test_cli_mesh_features import edges_at
from test_cli_mesh_landcover import EAST, WEST, coded, rel
from tin_engine.cli import app
from tin_engine.features import DEFAULT_VOCABULARY
from vtkread import read_vtk

SHIPPED = ("corine", "clc18_kode", "mapbiomas-c11", "ar50")
LAND_COVER, WATER = DEFAULT_VOCABULARY.mask("land_cover"), DEFAULT_VOCABULARY.mask("water")
CORINE_DISTANCES = (1.0, 50.0, 5.0)
OPTIONS = ("edges_repair", "edges_max_shift", "edges_outline_snap")
DERIVED = (
    "Default: from the land-cover class file's positional accuracy and the DEM cell (CORINE: {})."
)

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


def cover(path: Path, classes: str = "corine") -> tuple[str, ...]:
    """``path`` as the land cover under ``classes``, every class's borders chosen."""
    return ("--land-cover", str(path), "--land-cover-classes", classes, "--edges", "all")


def recorded(record: dict[str, Any]) -> tuple[float, float, float]:
    return (
        record["edges_repair_m"],
        record["edges_max_shift_m"],
        record["edges_outline_snap_m"],
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
        listed = [(1, "Forest", "#1f8d49"), (2, "Lake", "#2532e4")]
        path = write(tmp_path / "mine.toml", legend_text("mine", "rasputin mine", listed))
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
        assert preset["Name"] == maps.load_legend("mapbiomas-c11").preset
        assert "33 River, Lake and Ocean" in preset["Annotations"]
        assert len(preset["Annotations"]) == 2 * 34

    def test_t5_palette_refuses_a_land_cover_class_file(self, tmp_path: Path) -> None:
        path = write(tmp_path / "mine.toml", coded_map())
        code, said = invoke("palette", str(path))
        output = basenames(said)
        assert code == USAGE, output
        assert "mine.toml" in output and "attribute" in output, output

    def test_t5_palette_refuses_an_unknown_name_listing_the_shipped(self) -> None:
        code, said = invoke("palette", "corin")
        output = basenames(said)
        assert code == USAGE, output
        assert "unknown legend corin; use ar50, corine, mapbiomas-c11" in output, output
        for name in ("corine-water", "clc18_kode"):
            assert name not in output, name


# ------------------------------------------------------------- T5: mesh


class TestClassesOption:
    def test_t5_the_help_lists_the_shipped_names(self) -> None:
        text = mesh_params()["land_cover_classes"].help or ""
        for name in SHIPPED:
            assert name in text, name
        for gone in ("corine-water", "property"):
            assert gone not in text, gone

    def test_t5_an_unknown_name_is_refused_listing_the_shipped(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        code, said = invoke(
            "mesh",
            *("--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"),
            *("--out", str(tmp_path / "x.vtk"), *cover(corine, "corin")),
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
            *cover(corine, str(tmp_path / "nope.toml")),
        )
        output = basenames(said)
        assert code == USAGE, output
        assert "nope.toml" in output
        assert "absent.tif" not in output

    def test_t5_a_bad_file_is_refused_naming_it(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        path = write(
            tmp_path / "bad.toml", corine_like("bad", "40").replace("code = 311", "code = 0")
        )
        code, said = invoke(
            "mesh",
            *("--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"),
            *("--out", str(tmp_path / "x.vtk"), *cover(corine, str(path))),
        )
        output = basenames(said)
        assert code == USAGE, output
        assert "bad.toml" in output and "code" in output


# ------------------------------------------------------------- T7


class TestDerivedDistances:
    """T7 on `rough_dem` (h = 10 m): the record carries the values used, and
    `edges_of` is handed the same."""

    def test_t7_corine_with_no_option_is_todays_1_50_5(
        self,
        tmp_path: Path,
        bumpy: Path,
        square: Path,
        corine: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        spy = Spy(monkeypatch)
        _, record, _ = run(tmp_path, bumpy, square, *cover(corine))
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
            *cover(corine),
            *("--edges-max-shift", "20"),
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
        report, record, _ = run(tmp_path, bumpy, square, *cover(corine, str(path)))
        assert recorded(record) == want
        (request,) = spy.requests
        assert (request.repair_m, request.max_shift_m, request.outline_snap_m) == want
        assert "classes mine" in record["edges"]
        assert report

    def test_t7_a_user_map_with_one_override(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        path = write(tmp_path / "mine.toml", corine_like("mine", "40"))
        _, record, _ = run(
            tmp_path,
            bumpy,
            square,
            *cover(corine, str(path)),
            *("--edges-repair", "0.7"),
        )
        assert recorded(record) == (0.7, 20.0, 2.0)

    @pytest.mark.parametrize(
        ("param", "corine_value"),
        [("edges_repair", "1"), ("edges_max_shift", "50"), ("edges_outline_snap", "5")],
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
        data = write_geojson(
            tmp_path / "cover.geojson",
            [Feat("w", WEST, {"cover": "Forest"}), Feat("e", EAST, {"cover": "Lake"})],
        )
        path = write(tmp_path / "cover.toml", coded_map(codes="Forest = 1\nLake = 2\n"))
        run(tmp_path, bumpy, square, *cover(data, str(path)))
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        codes = vtk.cell_array("land_cover_code").values
        lines = len(vtk.lines)
        assert {1, 2} <= {int(c) for c in codes[lines:]}
        assert set(edges_at(vtk, rel(160.9, -60.2))) == {LAND_COVER | WATER}
        assert set(edges_at(vtk, rel(40.3, -60.2))) == {LAND_COVER}


# ------------------------------------------------------------- T9


class TestTheVocabularyInTheMesh:
    def test_t9_the_fingerprint_is_the_defaults(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        """A land-cover class file adds no edge name: the run's vocabulary is
        the default one."""
        path = write(tmp_path / "mine.toml", corine_like("mine", "40"))
        run(tmp_path, bumpy, square, *cover(corine, str(path)))
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        assert file_field(vtk, "feature_vocabulary") == DEFAULT_VOCABULARY.fingerprint()


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
        ``entry`` (`classes`, say), written to ``out``."""
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
            tmp_path, scene, monkeypatch, "named.gpkg", classes="mapbiomas-c11"
        )
        palette = (
            importlib.import_module("tin_engine.land_cover_classes")
            .load("mapbiomas-c11")
            .legend.colours
        )
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
            tmp_path, scene, monkeypatch, "named.vtk", classes="mapbiomas-c11"
        )
        before, after = read_vtk(plain.read_bytes()), read_vtk(named.read_bytes())
        assert sorted(after.field_data) == sorted(before.field_data)

    def test_t10_a_class_polygon_run_writes_the_mapbiomas_notice(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        data = write_geojson(
            tmp_path / "mb.geojson",
            [Feat(1, WEST, {"class": 3}), Feat(2, EAST, {"class": 33})],
        )
        _, record, _ = run(tmp_path, bumpy, square, *cover(data, "mapbiomas-c11"))
        notice = (
            importlib.import_module("tin_engine.land_cover_classes").load("mapbiomas-c11").notice
        )
        assert record["land_cover_notice"] == notice
        assert "MapBiomas Project - Collection 11" in notice
        # MapBiomas at 30 m over a 10 m cell: s = 30.
        assert recorded(record) == (0.3, 15.0, 1.5)


class TestLicenceInTheRecord:
    def test_t11_a_licence_is_recorded(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        text = corine_like("mine", "40").replace(
            'unlisted = "keep"\n', 'unlisted = "keep"\nlicence = "CC BY 4.0"\n'
        )
        path = write(tmp_path / "mine.toml", text)
        _, record, _ = run(tmp_path, bumpy, square, *cover(corine, str(path)))
        assert record["land_cover_licence"] == "CC BY 4.0"

    def test_t11_corine_writes_no_licence(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        _, record, _ = run(tmp_path, bumpy, square, *cover(corine))
        assert "land_cover_licence" not in record


class TestTheNamesInTheRecord:
    """`docs/increments/108-legend-split.md`, the Names table: the record's
    `edges` entry names the land-cover class file; the mesh text
    `land_cover_codes` keeps the word "legend" and takes the class file's
    name. `--edges water` (what `corine-water` chose) changes neither."""

    def test_edges_water(self, tmp_path: Path, bumpy: Path, square: Path, corine: Path) -> None:
        _, record, _ = run(
            tmp_path,
            bumpy,
            square,
            *("--land-cover", str(corine), "--land-cover-classes", "corine", "--edges", "water"),
        )
        assert ", classes corine: " in record["edges"], record["edges"]
        assert ", legend corine; " in record["land_cover_codes"], record["land_cover_codes"]

    def test_the_wording_of_the_edges_field(self) -> None:
        wording = importlib.import_module("tin_engine.run_record").WORDING
        assert wording["edges"] == (
            "Each source of edges: file, layer, land-cover class file or line kind, "
            "features, lines, vertices"
        )
        assert "features" not in wording


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
        assert len(preset["Annotations"]) == 2 * 19
        assert "81 Freshwater" in preset["Annotations"]
        assert (
            "5052 Open firm ground, barren to non-continuous vegetation cover"
            in (preset["Annotations"])
        )
        assert (
            "3031 Forest, coniferous: at least 50 % of the forested area covered by conifers"
            in (preset["Annotations"])
        )

    def test_13_the_triangles_borders_and_record(
        self, tmp_path: Path, bumpy: Path, square: Path, ar50: Path
    ) -> None:
        _, record, _ = run(tmp_path, bumpy, square, *cover(ar50, "ar50"))
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        codes = vtk.cell_array("land_cover_code").values
        assert {30, 81} <= {int(c) for c in codes[len(vtk.lines) :]}
        assert set(edges_at(vtk, rel(160.9, -60.2))) == {LAND_COVER | WATER}
        assert set(edges_at(vtk, rel(40.3, -60.2))) == {LAND_COVER}
        legend = importlib.import_module("tin_engine.land_cover_classes").load("ar50")
        assert record["land_cover_notice"] == legend.notice
        assert record["land_cover_licence"] == legend.licence
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
            *("--out", str(target), *cover(ar50, "ar50")),
        )
        assert code == 0, basenames(said)
        palette = (
            importlib.import_module("tin_engine.land_cover_classes").load("ar50").legend.colours
        )
        with closing(sqlite3.connect(target)) as conn:
            rows = conn.execute("SELECT code, land_cover, colour FROM triangles").fetchall()
        assert {30, 81} <= {r[0] for r in rows}
        for found, name, colour in rows:
            assert (name, colour) == palette[found], found


class TestAr50SplitsThroughTheCli:
    """`docs/increments/91-ar50-vegetation-cover.md`, section 8, test 7:
    coniferous forest (30, treslag 31) west, open firm ground with barren to
    non-continuous vegetation (50, vegetasjonsdekke 52) east."""

    SPLIT_TEXT = (
        "class 30 split by treslag into 3031, 3032, 3033, 3099, "
        "class 50 split by vegetasjonsdekke into 5051, 5052, 5054, 5055, 5099"
    )

    def cover(self, tmp_path: Path, *, second: bool = True) -> Path:
        west: dict[str, Any] = {"arealtype": 30}
        east: dict[str, Any] = {"arealtype": 50}
        if second:
            west["treslag"], east["vegetasjonsdekke"] = 31, 52
        return write_geojson(tmp_path / "split.geojson", [Feat(1, WEST, west), Feat(2, EAST, east)])

    def test_the_triangles_border_and_record(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        _, record, _ = run(tmp_path, bumpy, square, *cover(self.cover(tmp_path), "ar50"))
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        codes = {int(c) for c in vtk.cell_array("land_cover_code").values[len(vtk.lines) :]}
        assert {3031, 5052} <= codes
        assert not {30, 50} & codes
        # The shared side's south end: every edge there is a land-cover border.
        assert set(edges_at(vtk, rel(100.7, -60.2))) == {LAND_COVER}
        assert self.SPLIT_TEXT in record["land_cover_codes"]
        assert record["land_cover_codes"].startswith(
            f"AR50 arealtype code, attribute arealtype, {self.SPLIT_TEXT}, legend ar50; "
        )

    def test_the_gpkg_triangles_carry_the_sub_class_names_and_colours(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        import sqlite3
        from contextlib import closing

        target = tmp_path / "x.gpkg"
        code, said = invoke(
            "mesh",
            *("--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"),
            *("--out", str(target), *cover(self.cover(tmp_path), "ar50")),
        )
        assert code == 0, basenames(said)
        palette = (
            importlib.import_module("tin_engine.land_cover_classes").load("ar50").legend.colours
        )
        with closing(sqlite3.connect(target)) as conn:
            rows = conn.execute("SELECT code, land_cover, colour FROM triangles").fetchall()
        assert {3031, 5052} <= {r[0] for r in rows}
        for found, name, colour in rows:
            assert (name, colour) == palette[found], found

    def test_without_the_second_attributes_the_parents_codes(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        """Section 6: an export carrying only `arealtype` meshes as before."""
        run(tmp_path, bumpy, square, *cover(self.cover(tmp_path, second=False), "ar50"))
        vtk = read_vtk((tmp_path / "x.vtk").read_bytes())
        codes = {int(c) for c in vtk.cell_array("land_cover_code").values[len(vtk.lines) :]}
        assert {30, 50} <= codes
        assert not {3031, 5052} & codes
