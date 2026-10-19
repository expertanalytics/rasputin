"""`tin_engine.legends`: land-cover legends and palettes as data files (41).

`docs/increments/41-class-maps-as-data.md`, section 4 (the file, the loader,
the shipped maps, MapBiomas, the option), section 5 (the derived clean-up
distances), section 6 (the vocabulary) and section 7 (T1 to T12; the CLI
halves are in `test_cli_legends.py`). Fixtures by hand
(`legend_fixtures`); master's CORINE objects from the snapshot
(`legend_snapshot`), so the shipped CORINE files are checked against what
the code held before they existed, not against themselves.

The MapBiomas oracle is the Collection 11 legend CSV,
`tests/python/data/mapbiomas_c11_legend.csv`, fetched by `@tester` on
2026-10-12 from
`https://brasil.mapbiomas.org/wp-content/uploads/sites/3/2026/08/legend_code_mapbiomas_brazil_collection_11.csv`
(linked from `brasil.mapbiomas.org/en/codigos-de-legenda/`), unchanged: 33
classes, `class_id,class_name_pt_br,class_name_en,hex_code`.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `legends.load(value: str) -> Legend`: a shipped name, `property`, or a
  path (section 4's rule: ends in `.toml` or holds a path separator). Every
  refusal is a `ValueError` (a `FeatureError` is one) whose message names the
  file (its name, at least) and, for a bad key or value, the key.
- `load("property")` is today's `property` map over `DEFAULT_VOCABULARY`.
- `derived_cleanup(accuracy_m, cell)` returns a frozen model with fields
  `repair_m`, `tolerance_m` and `outline_snap_m` (`FeatureRequest`'s names).
- A MapBiomas class's label is the CSV's `class_name_en` verbatim (typos
  included, as the design review found them); colours as in the CSV.
- An uncoded map needs no `[edges]`, `notice` or `preset`; every map here
  states `unlisted`.
- `Legend.extra_names` is a tuple of the map's `[vocabulary].extra`;
  `Legend.accuracy_m` a float; `Legend.licence` `""` when absent.
- `feature_input._mask(the_legend, vocabulary, value, what)` keeps its signature.

RED at the commit that adds this file: `tin_engine.legends` does not exist
(`ModuleNotFoundError` in the `cm` fixture and in `feature_fixtures.legend`),
`tin_engine/data/legends/` does not exist, and the MapBiomas catalogue
entry has no `legend`.
"""

from __future__ import annotations

import csv
import importlib
import math
import tomllib
from importlib import resources
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from pydantic import ValidationError

import feature_fixtures as ff
from feature_fixtures import Feat, domain_of, square, write_geojson
from legend_fixtures import basenames, coded_map, uncoded_map, write
from legend_snapshot import CORINE_MAPS, CORINE_NATURAL, DATA, PRESET_NAME, SNAPSHOT
from tin_engine.features import DEFAULT_VOCABULARY

SHIPPED = ("corine", "corine-water", "clc18_kode", "mapbiomas-c11", "ar50")
#: T2's values, as the snapshot lists them.
VALUES: list[Any] = [value for value, _ in SNAPSHOT["maps"]["corine"]["masks"]]
LAND_COVER, WATER = DEFAULT_VOCABULARY.mask("land_cover"), DEFAULT_VOCABULARY.mask("water")
#: The snapshot's domain and polygon for `opened` (see `legend_snapshot`).
DOMAIN = domain_of(square(0, 0, 100, 100))
INSIDE = square(10, 10, 50, 50)
WEST, EAST = square(10, 10, 50, 50), square(50, 10, 90, 50)
MAPBIOMAS_NOTICE = (
    "Contains data from the MapBiomas Project - Collection 11 of the Annual Series of Land "
    "Use and Land Cover Maps of Brazil (CC BY-SA 4.0)"
)


@pytest.fixture(scope="module")
def cm() -> ModuleType:
    return importlib.import_module("tin_engine.legends")


def shipped_dir() -> Any:
    return resources.files("tin_engine") / "data" / "legends"


def shipped_toml(name: str) -> dict[str, Any]:
    return tomllib.loads((shipped_dir() / f"{name}.toml").read_text(encoding="utf-8"))


def mapbiomas_csv() -> dict[int, tuple[str, str]]:
    """The Collection 11 legend CSV: code -> (English label, colour)."""
    with (DATA / "mapbiomas_c11_legend.csv").open(encoding="utf-8", newline="") as f:
        return {int(r["class_id"]): (r["class_name_en"], r["hex_code"]) for r in csv.DictReader(f)}


def refused(cm: ModuleType, path: Path, *words: str) -> str:
    """`load(path)` refused naming the file and every one of ``words``."""
    with pytest.raises(ValueError) as caught:
        cm.load(str(path))
    message = basenames(str(caught.value))
    for word in (path.name, *words):
        assert word in message, f"{word!r} not in: {message}"
    return message


# ------------------------------------------------------------- T1, T12


class TestShippedMaps:
    def test_t12_the_shipped_files_are_package_data(self) -> None:
        """T1 and T12: exactly the shipped ones, found through `importlib.resources`."""
        names = {p.name for p in shipped_dir().iterdir() if p.name.endswith(".toml")}
        assert names == {f"{name}.toml" for name in SHIPPED}

    @pytest.mark.parametrize("name", SHIPPED)
    def test_t1_each_loads_by_name_and_its_name_is_its_stem(
        self, cm: ModuleType, name: str
    ) -> None:
        assert shipped_toml(name)["name"] == name
        assert cm.load(name).name == name

    def test_t1_the_three_corine_class_tables_are_equal(self) -> None:
        corine, water, kode = (shipped_toml(name) for name in CORINE_MAPS)
        assert corine["class"] == water["class"] == kode["class"]
        assert corine["no_class"] == water["no_class"] == kode["no_class"]
        assert len(corine["class"]) == 44


# ------------------------------------------------------------------- T2


class TestCorineBehaviourUnchanged:
    """T2: masks and codes equal master's for every CORINE map and value."""

    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t2_name_attribute_code_system_and_notice(self, cm: ModuleType, name: str) -> None:
        want = SNAPSHOT["maps"][name]
        got = cm.load(name)
        assert (got.name, got.attribute, got.codes, got.notice) == (
            want["name"],
            want["attribute"],
            want["codes"],
            want["notice"],
        )

    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t2_the_masks_of_the_designs_values(self, cm: ModuleType, name: str) -> None:
        fi, the_legend = ff.feature_input(), cm.load(name)
        for value, mask in SNAPSHOT["maps"][name]["masks"]:
            assert fi._mask(the_legend, DEFAULT_VOCABULARY, value, "t2") == mask, value

    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t2_the_mask_of_every_level_3_class(self, cm: ModuleType, name: str) -> None:
        fi, the_legend = ff.feature_input(), cm.load(name)
        got = {
            code: fi._mask(the_legend, DEFAULT_VOCABULARY, code, "t2")
            for code in SNAPSHOT["maps"][name]["every_class"]
        }
        assert got == SNAPSHOT["maps"][name]["every_class"]

    def test_t2_water_and_land_as_the_design_states_them(self) -> None:
        """The snapshot's own numbers, as section 7 words them: water 384,
        other codes 128, `corine-water`'s land codes dropped."""
        assert (LAND_COVER | WATER, LAND_COVER) == (384, 128)
        masks = SNAPSHOT["maps"]["corine"]["every_class"]
        assert {c for c, m in masks.items() if m == 384} == {"511", "512", "521", "522", "523"}
        dropped = SNAPSHOT["maps"]["corine-water"]["every_class"]
        assert {c for c, m in dropped.items() if m is not None} == {
            "511",
            "512",
            "521",
            "522",
            "523",
        }

    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t2_codes_and_masks_through_open_features(
        self, cm: ModuleType, tmp_path: Path, name: str
    ) -> None:
        fi, the_legend = ff.feature_input(), cm.load(name)
        for i, (value, want) in enumerate(SNAPSHOT["maps"][name]["opened"]):
            path = write_geojson(
                tmp_path / f"{i}.geojson", [Feat(1, INSIDE, {the_legend.attribute: value})]
            )
            request = fi.FeatureRequest(sources=(fi.FeatureSource(path=path, legend=the_legend),))
            if want == "refused":
                with pytest.raises(fi.FeatureError):
                    fi.open_features(request, DOMAIN, ff.UTM33)
                continue
            fs = fi.open_features(request, DOMAIN, ff.UTM33)
            assert [[f.code, f.mask] for f in fs.features] == want, value


# ------------------------------------------------------------------- T3


class TestCorinePaletteUnchanged:
    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t3_the_palette_is_masters_in_masters_order(self, cm: ModuleType, name: str) -> None:
        palette = cm.load(name).palette
        assert list(palette.items()) == list(CORINE_NATURAL.items())
        assert next(iter(palette)) == 0

    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t3_the_preset_title(self, cm: ModuleType, name: str) -> None:
        assert cm.load(name).preset == PRESET_NAME == "rasputin CORINE natural"


# ------------------------------------------------------------------- T4


class TestTheBaseFilesLoad:
    """The refusals below each change one line of these; both load as written."""

    def test_the_coded_base(self, cm: ModuleType, tmp_path: Path) -> None:
        got = cm.load(str(write(tmp_path / "base.toml", coded_map())))
        assert (got.name, got.attribute, got.codes, got.accuracy_m) == (
            "mine",
            "cover",
            "Test cover code",
            40.0,
        )
        assert dict(got.palette) == {
            0: ("No polygon; constraint lines", "#3c3c3c"),
            1: ("Forest", "#1f8d49"),
            2: ("Lake", "#2532e4"),
        }
        assert dict(got.classes) == {"1": ("land_cover",), "2": ("land_cover", "water")}
        assert got.otherwise == ("land_cover",)
        assert (got.preset, got.notice, got.licence) == ("rasputin mine", "Test data, made up", "")

    def test_the_uncoded_base(self, cm: ModuleType, tmp_path: Path) -> None:
        got = cm.load(str(write(tmp_path / "levees.toml", uncoded_map())))
        assert (got.name, got.attribute, got.codes, got.palette) == ("levees", "kind", "", None)
        assert dict(got.classes) == {"levee": ("levee",), "road": ("road",)}
        assert got.otherwise == "refuse"
        assert got.extra_names == ("levee",)

    @pytest.mark.parametrize(
        ("unlisted", "otherwise"),
        [("land", ("land_cover",)), ("drop", "drop"), ("refuse", "refuse")],
    )
    def test_unlisted_becomes_otherwise(
        self, cm: ModuleType, tmp_path: Path, unlisted: str, otherwise: Any
    ) -> None:
        got = cm.load(str(write(tmp_path / "u.toml", coded_map(unlisted=unlisted))))
        assert got.otherwise == otherwise

    def test_a_class_edges_override_and_drop(self, cm: ModuleType, tmp_path: Path) -> None:
        """Section 4's table: a class's own `edges` wins; `"drop"` is None."""
        classes = (
            '[[class]]\ncode = 1\nlabel = "Forest"\ncolour = "#1f8d49"\nedges = ["wall"]\n\n'
            '[[class]]\ncode = 2\nlabel = "Lake"\ncolour = "#2532e4"\nedges = "drop"\n'
        )
        got = cm.load(str(write(tmp_path / "o.toml", coded_map(classes=classes))))
        assert dict(got.classes) == {"1": ("wall",), "2": None}


ACC = "positional_accuracy_m"
#: (file stem, the line changed, its replacement, words the refusal names).
#: The stem is the file's name, which every refusal must carry.
ONE_LINE = [
    ("no-attribute", 'attribute = "cover"\n', "", ("attribute",)),
    ("colour-word", 'colour = "#1f8d49"', 'colour = "green"', ("colour",)),
    ("colour-upper", 'colour = "#1f8d49"', 'colour = "#1F8D49"', ("colour",)),
    ("colour-short", 'colour = "#1f8d49"', 'colour = "#1f8d4"', ("colour",)),
    ("duplicate-code", "code = 2\n", "code = 1\n", ("code",)),
    ("code-zero", "code = 1\n", "code = 0\n", ("code",)),
    ("code-two-31", "code = 1\n", "code = 2147483648\n", ("code",)),
    ("unknown-key", 'unlisted = "land"\n', 'unlisted = "land"\nlegend = "x"\n', ("legend",)),
    ("unknown-class-key", 'label = "Forest"', 'label = "Forest"\nlabl = "x"', ("labl",)),
    ("unlisted-word", 'unlisted = "land"', 'unlisted = "keep"', ("unlisted",)),
    ("edges-name", 'land = ["land_cover"]', 'land = ["land_cover", "levee"]', ("levee",)),
    ("no-accuracy", "positional_accuracy_m = 40\n", "", ("positional_accuracy_m",)),
    ("accuracy-zero", f"{ACC} = 40", f"{ACC} = 0", (ACC,)),
    ("accuracy-negative", f"{ACC} = 40", f"{ACC} = -5", (ACC,)),
    ("accuracy-nan", f"{ACC} = 40", f"{ACC} = nan", (ACC,)),
    ("accuracy-inf", f"{ACC} = 40", f"{ACC} = inf", (ACC,)),
    ("coded-with-value", 'label = "Forest"', 'label = "Forest"\nvalue = "forest"', ("value",)),
]  # fmt: skip


class TestFileRefusals:
    """T4: each refusal names the file, and the key where there is one."""

    @pytest.mark.parametrize(
        ("stem", "old", "new", "words"), ONE_LINE, ids=[c[0] for c in ONE_LINE]
    )
    def test_t4_one_line_changed(
        self, cm: ModuleType, tmp_path: Path, stem: str, old: str, new: str, words: tuple[str, ...]
    ) -> None:
        text = coded_map()
        assert old in text, "the fixture changed under this test"
        refused(cm, write(tmp_path / f"{stem}.toml", text.replace(old, new, 1)), *words)

    @pytest.mark.parametrize(
        ("stem", "codes", "word"),
        [
            ("codes-zero", "Forest = 0\n", "Forest"),
            ("codes-two-31", "Forest = 2147483648\n", "Forest"),
            ("codes-no-class", "Forest = 7\n", "Forest"),
        ],
    )
    def test_t4_a_codes_value_outside_the_range_or_naming_no_class(
        self, cm: ModuleType, tmp_path: Path, stem: str, codes: str, word: str
    ) -> None:
        refused(cm, write(tmp_path / f"{stem}.toml", coded_map(codes=codes)), word)

    def test_t4_an_uncoded_class_with_a_colour(self, cm: ModuleType, tmp_path: Path) -> None:
        text = uncoded_map().replace('edges = ["road"]', 'edges = ["road"]\ncolour = "#6e6e6e"')
        refused(cm, write(tmp_path / "uncoded-colour.toml", text), "colour")

    def test_t4_an_uncoded_class_with_a_code(self, cm: ModuleType, tmp_path: Path) -> None:
        text = uncoded_map().replace('edges = ["road"]', 'edges = ["road"]\ncode = 3')
        refused(cm, write(tmp_path / "uncoded-code.toml", text), "code")

    def test_t4_a_toml_syntax_error_names_its_line(self, cm: ModuleType, tmp_path: Path) -> None:
        text = coded_map().replace('attribute = "cover"', 'attribute = "cover')
        assert text.splitlines()[1] == 'attribute = "cover'
        refused(cm, write(tmp_path / "broken.toml", text), "line 2")


# ------------------------------------------------------------------- T5


class TestNameOrPath:
    def test_t5_a_name_is_the_shipped_file(self, cm: ModuleType) -> None:
        with resources.as_file(shipped_dir() / "corine.toml") as path:
            assert cm.load(str(path)) == cm.load("corine")

    def test_t5_a_toml_suffix_or_a_separator_is_a_path(
        self, cm: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        write(tmp_path / "x.toml", coded_map("x"))
        write(tmp_path / "dir" / "x", coded_map("y"))
        write(tmp_path / "corine.toml", coded_map("not-the-shipped-one"))
        monkeypatch.chdir(tmp_path)
        assert cm.load("./x.toml").name == "x"
        assert cm.load("x.toml").name == "x"
        assert cm.load("dir/x").name == "y"
        assert cm.load("corine.toml").name == "not-the-shipped-one"
        assert cm.load("corine").name == "corine"

    def test_t5_an_unknown_name_is_refused_listing_the_shipped_names(self, cm: ModuleType) -> None:
        with pytest.raises(ValueError) as caught:
            cm.load("corin")
        for word in ("corin", *SHIPPED, "property"):
            assert word in str(caught.value), word

    @pytest.mark.parametrize("value", ["./nope.toml", "nowhere/nope"])
    def test_t5_a_missing_path_is_refused_naming_it(
        self, cm: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, value: str
    ) -> None:
        monkeypatch.chdir(tmp_path)
        with pytest.raises(ValueError, match="nope"):
            cm.load(value)

    def test_t5_property_is_todays_identity_map(self, cm: ModuleType) -> None:
        got = cm.load("property")
        assert (got.name, got.attribute, got.codes, got.notice) == ("property", "property", "", "")
        assert dict(got.classes) == {p.name: (p.name,) for p in DEFAULT_VOCABULARY.properties}
        assert got.otherwise == "refuse"
        assert got.palette is None


# ------------------------------------------------------------------- T6


class TestDerivedCleanup:
    """T6, section 5.1: s = max(a, h); repair s/100, simplify s/2, snap
    min(s/20, 50). Exact equality: division gives these floats exactly
    (`100.0/100 == 1.0`, `30/100 == 0.3`, and so on, in CPython 3.12+)."""

    @pytest.mark.parametrize(
        ("accuracy", "cell", "want"),
        [
            (100.0, 10.0, (1.0, 50.0, 5.0)),
            (100.0, 31.0, (1.0, 50.0, 5.0)),
            (100.0, 100.0, (1.0, 50.0, 5.0)),
            (100.0, 200.0, (2.0, 100.0, 10.0)),
            (30.0, 30.0, (0.3, 15.0, 1.5)),
            (4.0, 10.0, (0.1, 5.0, 0.5)),
            (40.0, 10.0, (0.4, 20.0, 2.0)),
            (100.0, 1000.0, (10.0, 500.0, 50.0)),
            (100.0, 1500.0, (15.0, 750.0, 50.0)),
            (2.5, 1.0, (0.025, 1.25, 0.125)),
        ],
    )
    def test_t6_the_three_distances(
        self, cm: ModuleType, accuracy: float, cell: float, want: tuple[float, float, float]
    ) -> None:
        got = cm.derived_cleanup(accuracy, cell)
        assert (got.repair_m, got.tolerance_m, got.outline_snap_m) == want

    def test_t6_the_snap_stays_under_the_margin(self, cm: ModuleType) -> None:
        """Section 5.1's cap: half `feature_input.MARGIN`, 100 m; checked at
        cells up to 2 km, where s / 20 would reach the margin itself."""
        for cell in (999.0, 1000.0, 1001.0, 2000.0):
            snap = cm.derived_cleanup(100.0, cell).outline_snap_m
            assert snap <= 50.0 < ff.feature_input().MARGIN
        assert cm.derived_cleanup(100.0, 999.0).outline_snap_m == 999.0 / 20

    def test_t6_the_result_is_frozen(self, cm: ModuleType) -> None:
        got = cm.derived_cleanup(100.0, 10.0)
        with pytest.raises((TypeError, ValueError, ValidationError)):
            got.repair_m = 3.0
        assert all(math.isfinite(v) for v in (got.repair_m, got.tolerance_m, got.outline_snap_m))


# ------------------------------------------------------------------- T8


def cover_source(tmp_path: Path, unlisted: str, *, bog: bool = True) -> Any:
    fi = ff.feature_input()
    the_legend = cm_load(write(tmp_path / f"cover-{unlisted}.toml", coded_map(
        unlisted=unlisted, codes='Forest = 1\nLake = 2\n'
    )))  # fmt: skip
    features = [Feat("west", WEST, {"cover": "Forest"}), Feat("east", EAST, {"cover": "Lake"})]
    if bog:
        features.append(Feat("bog-7", square(20, 60, 40, 80), {"cover": "Bog"}))
    path = write_geojson(tmp_path / f"cover-{unlisted}.geojson", features)
    return fi.FeatureRequest(sources=(fi.FeatureSource(path=path, legend=the_legend),))


def cm_load(path: Path) -> Any:
    return importlib.import_module("tin_engine.legends").load(str(path))


class TestTextCodes:
    """T8 (B6): `[codes]` turns attribute text into class codes."""

    def test_t8_the_codes_table_is_kept(self, tmp_path: Path) -> None:
        the_legend = cm_load(write(tmp_path / "c.toml", coded_map(codes="Forest = 1\nLake = 2\n")))
        assert dict(the_legend.string_codes) == {"Forest": 1, "Lake": 2}
        assert dict(the_legend.classes)["Lake"] == ("land_cover", "water")
        assert dict(the_legend.classes)["Forest"] == ("land_cover",)

    def test_t8_text_values_carry_their_codes_and_the_lake_water(self, tmp_path: Path) -> None:
        fi = ff.feature_input()
        fs = fi.open_features(cover_source(tmp_path, "drop", bog=False), DOMAIN, ff.UTM33)
        got = {f.fid: (f.code, f.mask) for f in fs.features}
        assert got == {"west": (1, LAND_COVER), "east": (2, LAND_COVER | WATER)}

    @pytest.mark.parametrize("unlisted", ["refuse", "land"])
    def test_t8_an_unlisted_text_is_refused_naming_file_feature_and_value(
        self, tmp_path: Path, unlisted: str
    ) -> None:
        fi = ff.feature_input()
        with pytest.raises(fi.FeatureError) as caught:
            fi.open_features(cover_source(tmp_path, unlisted), DOMAIN, ff.UTM33)
        for word in (f"cover-{unlisted}.geojson", "bog-7", "Bog"):
            assert word in str(caught.value), word

    def test_t8_an_unlisted_text_is_dropped_with_drop(self, tmp_path: Path) -> None:
        fi = ff.feature_input()
        fs = fi.open_features(cover_source(tmp_path, "drop"), DOMAIN, ff.UTM33)
        assert sorted(f.fid for f in fs.features) == ["east", "west"]

    def test_t8_without_a_codes_table_a_text_is_refused_as_today(self, tmp_path: Path) -> None:
        fi = ff.feature_input()
        the_legend = cm_load(write(tmp_path / "plain.toml", coded_map(unlisted="land")))
        path = write_geojson(tmp_path / "plain.geojson", [Feat("f", WEST, {"cover": "Forest"})])
        request = fi.FeatureRequest(sources=(fi.FeatureSource(path=path, legend=the_legend),))
        with pytest.raises(fi.FeatureError, match="Forest"):
            fi.open_features(request, DOMAIN, ff.UTM33)


# ------------------------------------------------------------------- T9


class TestVocabularyExtras:
    def test_t9_an_extra_repeating_a_default_name_is_refused(
        self, cm: ModuleType, tmp_path: Path
    ) -> None:
        refused(cm, write(tmp_path / "river.toml", uncoded_map(extra='"levee", "river"')), "river")

    def test_t9_twenty_two_extras_load_and_twenty_three_are_refused(
        self, cm: ModuleType, tmp_path: Path
    ) -> None:
        """Bits 10 to 31 hold 22 names (section 6)."""
        names = ["levee", *(f"extra{i}" for i in range(1, 23))]
        quoted = ", ".join(f'"{n}"' for n in names)
        many = write(tmp_path / "many.toml", uncoded_map(extra=quoted))
        refused(cm, many, "23")
        fits = ", ".join(f'"{n}"' for n in names[:22])
        assert (
            len(cm.load(str(write(tmp_path / "fits.toml", uncoded_map(extra=fits)))).extra_names)
            == 22
        )

    def test_t9_an_edges_name_from_the_extras_is_accepted(
        self, cm: ModuleType, tmp_path: Path
    ) -> None:
        got = cm.load(str(write(tmp_path / "ok.toml", uncoded_map())))
        assert dict(got.classes)["levee"] == ("levee",)


# ------------------------------------------------------------ T10, T11


class TestMapBiomas:
    def test_t10_the_design_anchors(self, cm: ModuleType) -> None:
        palette = cm.load("mapbiomas-c11").palette
        assert palette[3] == ("Forest Formation", "#1f8d49")
        assert palette[4] == ("Savanna Formation", "#7dc975")
        assert palette[5] == ("Mangrove", "#04381d")
        assert palette[33] == ("River, Lake and Ocean", "#2532e4")
        assert palette[31] == ("Aquaculture", "#091077")

    def test_t10_every_class_is_the_legend_csvs(self, cm: ModuleType) -> None:
        palette = dict(cm.load("mapbiomas-c11").palette)
        assert len(mapbiomas_csv()) == 33
        assert {code: v for code, v in palette.items() if code} == mapbiomas_csv()

    def test_t10_no_class_water_and_the_notice(self, cm: ModuleType) -> None:
        got = cm.load("mapbiomas-c11")
        assert got.palette[0] == ("No class; constraint lines", "#3c3c3c")
        water = {key for key, edges in got.classes.items() if edges and "water" in edges}
        assert water == {"33", "31"}
        assert "26" not in got.classes
        assert got.classes["33"] == ("land_cover", "water")
        assert got.classes["3"] == ("land_cover",)
        assert got.otherwise == ("land_cover",)
        assert got.notice == MAPBIOMAS_NOTICE
        assert (got.attribute, got.accuracy_m) == ("class", 30.0)
        assert got.codes
        assert got.licence == ""

    def test_t10_the_catalogue_entry_names_the_map(self) -> None:
        sources = importlib.import_module("tin_engine.sources")
        assert sources.SOURCES["mapbiomas-c11"].legend == "mapbiomas-c11"


class TestLicence:
    def test_t11_a_licence_is_kept_and_corines_is_empty(
        self, cm: ModuleType, tmp_path: Path
    ) -> None:
        path = write(tmp_path / "lic.toml", coded_map(licence="CC BY 4.0"))
        assert cm.load(str(path)).licence == "CC BY 4.0"
        assert all(cm.load(name).licence == "" for name in CORINE_MAPS)


# ------------------------------------------------------- the AR50 legend


#: `docs/increments/73-gml-vector-source.md`, section 4.5's table: code,
#: label, colour, and NIBIO's RGB (None where NIBIO leaves the class unfilled).
AR50_TABLE = (
    (10, "Developed area and traffic", "#e86474", (232, 100, 116)),
    (20, "Agriculture", "#ffd16e", (255, 209, 110)),
    (30, "Forest", "#9ecc73", (158, 204, 115)),
    (50, "Open firm ground", "#d9d9d9", (217, 217, 217)),
    (60, "Marsh", "#73dfff", (115, 223, 255)),
    (70, "Snow and glacier", "#ffffff", (255, 255, 255)),
    (81, "Freshwater", "#91e7ff", (145, 231, 255)),
    (82, "Ocean", "#d2ffff", (210, 255, 255)),
    (99, "Not surveyed", "#8c8c8c", None),
)
AR50_NOTICE = (
    "Contains data under the Norwegian licence for Open Government Data (NLOD) "
    "distributed by NIBIO: AR50 land-resource map 1:50 000, re-encoded"
)
AR50_LICENCE = (
    "Norwegian licence for Open Government Data (NLOD) 1.0, https://data.norge.no/nlod/en/1.0"
)


class TestAr50:
    """Test 11 of `docs/increments/73-gml-vector-source.md`: the shipped
    `ar50` legend over AR50's land type, `arealtype`."""

    def test_the_table_colours_are_nibios_rgb(self) -> None:
        """The table's own check: each colour is its RGB in hex, so a slip
        between the two columns fails here, not only in the file."""
        for code, _, colour, rgb in AR50_TABLE:
            if rgb is not None:
                assert colour == "#{:02x}{:02x}{:02x}".format(*rgb), code

    def test_the_header(self, cm: ModuleType) -> None:
        got = cm.load("ar50")
        assert (got.name, got.attribute, got.codes) == ("ar50", "arealtype", "AR50 arealtype code")
        assert got.accuracy_m == 2.0
        assert got.notice == AR50_NOTICE
        assert got.licence == AR50_LICENCE
        assert got.preset == "rasputin AR50 NIBIO"

    def test_the_palette_is_exactly_the_tables(self, cm: ModuleType) -> None:
        want = {0: ("No polygon; constraint lines", "#3c3c3c")}
        want |= {code: (label, colour) for code, label, colour, _ in AR50_TABLE}
        assert dict(cm.load("ar50").palette) == want

    def test_fresh_water_and_ocean_carry_water_and_the_rest_land_cover(
        self, cm: ModuleType
    ) -> None:
        got = cm.load("ar50")
        land = ("land_cover",)
        want = {str(code): land for code, *_ in AR50_TABLE}
        want |= {"81": ("land_cover", "water"), "82": ("land_cover", "water")}
        assert dict(got.classes) == want
        assert got.otherwise == land
