"""`tin_engine.land_cover_classes`: the land-cover class files, what a class
polygon's attribute means.

`docs/increments/41-class-maps-as-data.md`, section 4 (the file, the loader,
the shipped maps, MapBiomas, the option), section 5 (the derived clean-up
distances), section 6 (the vocabulary) and section 7 (T1 to T12; the CLI
halves are in `test_cli_legends.py`); `docs/increments/108-legend-split.md`,
section 5, tests 1, 3 and 4 (the rules file apart from its legend, `[cells]`,
the colours and cell mappings unchanged, AR50's known codes);
`docs/increments/111-one-land-cover-map.md`, sections 2 and 8 (meshing
rules become land-cover class files: `[edges]`, `class.edges`, `class.value`
and `[vocabulary]` go, `unlisted = "land"` becomes `"keep"`, `code_system` is
required, `property` and `corine-water` go). Fixtures by
hand (`legend_fixtures`); master's CORINE objects from the snapshots
(`legend_snapshot`, `tests/python/data/legends_5de1144d.json`), so the
shipped files are checked against what the code held before, not against
themselves.

The MapBiomas oracle is the Collection 11 legend CSV,
`tests/python/data/mapbiomas_c11_legend.csv`, fetched by `@tester` on
2026-10-12 from
`https://brasil.mapbiomas.org/wp-content/uploads/sites/3/2026/08/legend_code_mapbiomas_brazil_collection_11.csv`
(linked from `brasil.mapbiomas.org/en/codigos-de-legenda/`), unchanged: 33
classes, `class_id,class_name_pt_br,class_name_en,hex_code`.

`legends_5de1144d.json`, written by `@tester` from `legends.load(name)` at
`5de1144d` for each shipped name: `preset`; `colours` as `[code, label,
colour]` in the palette's order; `categories` (code as text to category);
`prefer` (category to code); `water_codes`, sorted.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `land_cover_classes.load(value: str) -> LandCoverClasses`: a shipped name
  or a path (ends in `.toml` or holds a path separator). Every refusal is a
  `ValueError` whose message names the file (its name, at least) and, for a
  bad key or value, the key. The shipped files are
  `tin_engine/data/land_cover_classes/<name>.toml`.
- `LandCoverClasses.classes` holds the listed values: each class's code as
  text, and each `[codes]` text.
- A border's edge names follow from its class: `land_cover`, and `water` for
  a water class, read through `edges_of` with every class chosen. A value
  `"0511"` is code 511, a water class, so its border carries `water` too
  (master keyed the names on the text and gave `land_cover` only).
- `land_cover_classes.known_codes(rules) -> frozenset[int]`: the codes a class
  polygon may carry; 0 among them only when the rules file says
  `accept_no_class = true` (only with `code_system`, and the legend must
  colour code 0); the legend's `[no_class]` decides nothing. Every shipped
  rules file says it.
- A refusal of a legend-only key in a rules file names the key and the word
  "legend"; one of a `[cells]` misplacement names "[cells]".
- `derived_cleanup(accuracy_m, cell)` lives in `land_cover_classes`, unchanged.
- A MapBiomas class's label is the CSV's `class_name_en` verbatim; colours as
  in the CSV.
"""

from __future__ import annotations

import csv
import importlib
import json
import math
import re
import tomllib
from importlib import resources
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from pydantic import ValidationError

import feature_fixtures as ff
from feature_fixtures import Feat, domain_of, square, write_geojson
from legend_fixtures import basenames, coded_map, corine_like, legend_text, write
from legend_snapshot import CORINE_MAPS, CORINE_NATURAL, DATA, PRESET_NAME, SNAPSHOT
from tin_engine.features import DEFAULT_VOCABULARY

SHIPPED = ("corine", "clc18_kode", "mapbiomas-c11", "ar50")
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
    return importlib.import_module("tin_engine.land_cover_classes")


def shipped_dir() -> Any:
    return resources.files("tin_engine") / "data" / "land_cover_classes"


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


def edges_of_values(tmp_path: Path, name: str, values: list[Any]) -> list[list[Any]]:
    """Each value on one small square inside `DOMAIN`, read as `--land-cover`
    with the class file ``name`` and `--edges all`: `[code, mask]` per
    feature, in the values' order."""
    attribute = ff.classes(name).attribute
    squares = [square(2 + 2 * (i % 40), 2 + 2 * (i // 40), 3 + 2 * (i % 40), 3 + 2 * (i // 40))
               for i in range(len(values))]  # fmt: skip
    path = write_geojson(
        tmp_path / f"{name}.geojson",
        [Feat(i, g, {attribute: v}) for i, (g, v) in enumerate(zip(squares, values, strict=True))],
    )
    fs = ff.open_one(path, DOMAIN, name)
    assert [f.fid for f in fs.features] == list(range(len(values)))
    return [[f.code, f.mask] for f in fs.features]


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

    def test_t1_the_two_corine_class_tables_are_equal(self) -> None:
        corine, kode = (shipped_toml(name) for name in CORINE_MAPS)
        assert corine["class"] == kode["class"]
        assert corine["legend"] == kode["legend"] == "corine"
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
    def test_t2_the_mask_of_every_level_3_class(
        self, cm: ModuleType, tmp_path: Path, name: str
    ) -> None:
        """Through `edges_of`, every class chosen: master's masks."""
        want = SNAPSHOT["maps"][name]["every_class"]
        got = edges_of_values(tmp_path, name, list(want))
        assert {str(code): mask for code, mask in got} == want

    def test_t2_water_and_land_as_the_design_states_them(self) -> None:
        """The snapshot's own numbers, as section 7 words them: water 384,
        other codes 128."""
        assert (LAND_COVER | WATER, LAND_COVER) == (384, 128)
        masks = SNAPSHOT["maps"]["corine"]["every_class"]
        assert {c for c, m in masks.items() if m == 384} == {"511", "512", "521", "522", "523"}

    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t2_codes_and_masks_through_edges_of(
        self, cm: ModuleType, tmp_path: Path, name: str
    ) -> None:
        """Master's codes and refusals; the mask follows the code, so
        `"0511"` (code 511) carries `water` where master's did not."""
        fi, the_legend = ff.feature_input(), cm.load(name)
        water = {511, 512, 521, 522, 523}
        for i, (value, want) in enumerate(SNAPSHOT["maps"][name]["opened"]):
            path = write_geojson(
                tmp_path / f"{i}.geojson", [Feat(1, INSIDE, {the_legend.attribute: value})]
            )
            if want == "refused":
                with pytest.raises(fi.FeatureError):
                    ff.open_one(path, DOMAIN, name)
                continue
            fs = ff.open_one(path, DOMAIN, name)
            ((code, mask),) = want
            assert mask == (LAND_COVER | WATER if code in water else LAND_COVER) or value == "0511"
            by_code = LAND_COVER | WATER if code in water else LAND_COVER
            assert [[f.code, f.mask] for f in fs.features] == [[code, by_code]], value


# ------------------------------------------------------------------- T3


class TestCorinePaletteUnchanged:
    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t3_the_palette_is_masters_in_masters_order(self, cm: ModuleType, name: str) -> None:
        palette = cm.load(name).legend.colours
        assert list(palette.items()) == list(CORINE_NATURAL.items())
        assert next(iter(palette)) == 0

    @pytest.mark.parametrize("name", CORINE_MAPS)
    def test_t3_the_preset_title(self, cm: ModuleType, name: str) -> None:
        assert cm.load(name).legend.preset == PRESET_NAME == "rasputin CORINE natural"


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
        assert dict(got.legend.colours) == {
            0: ("No polygon; constraint lines", "#3c3c3c"),
            1: ("Forest", "#1f8d49"),
            2: ("Lake", "#2532e4"),
        }
        assert got.classes == frozenset({"1", "2"})
        assert got.water_codes == frozenset({2})
        assert got.otherwise == "keep"
        assert (got.legend.preset, got.notice, got.licence) == (
            "rasputin mine",
            "Test data, made up",
            "",
        )

    @pytest.mark.parametrize("unlisted", ["keep", "drop", "refuse"])
    def test_unlisted_becomes_otherwise(
        self, cm: ModuleType, tmp_path: Path, unlisted: str
    ) -> None:
        got = cm.load(str(write(tmp_path / "u.toml", coded_map(unlisted=unlisted))))
        assert got.otherwise == unlisted


ACC = "positional_accuracy_m"
#: (file stem, the line changed, its replacement, words the refusal names).
#: The stem is the file's name, which every refusal must carry.
ONE_LINE = [
    ("no-attribute", 'attribute = "cover"\n', "", ("attribute",)),
    ("duplicate-code", "code = 2\n", "code = 1\n", ("code",)),
    ("code-zero", "code = 1\n", "code = 0\n", ("code",)),
    ("code-two-31", "code = 1\n", "code = 2147483648\n", ("code",)),
    ("unknown-key", 'unlisted = "keep"\n', 'unlisted = "keep"\nlayer = "x"\n', ("layer",)),
    ("unknown-class-key", "code = 1\n", 'code = 1\nlabl = "x"\n', ("labl",)),
    ("unlisted-land", 'unlisted = "keep"', 'unlisted = "land"', ("unlisted",)),
    ("class-edges", "code = 1\n", 'code = 1\nedges = ["land_cover"]\n', ("edges",)),
    ("no-code-system", 'code_system = "Test cover code"\n', "", ("code_system",)),
    ("no-accuracy", "positional_accuracy_m = 40\n", "", ("positional_accuracy_m",)),
    ("accuracy-zero", f"{ACC} = 40", f"{ACC} = 0", (ACC,)),
    ("accuracy-negative", f"{ACC} = 40", f"{ACC} = -5", (ACC,)),
    ("accuracy-nan", f"{ACC} = 40", f"{ACC} = nan", (ACC,)),
    ("accuracy-inf", f"{ACC} = 40", f"{ACC} = inf", (ACC,)),
    ("coded-with-value", "code = 1\n", 'code = 1\nvalue = "forest"\n', ("value",)),
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

    @pytest.mark.parametrize(
        ("stem", "table", "word"),
        [
            ("edges-table", '[edges]\nwater = ["land_cover", "water"]\nland = ["land_cover"]\n',
             "edges"),
            ("vocabulary-table", '[vocabulary]\nextra = ["levee"]\n', "vocabulary"),
        ],
    )  # fmt: skip
    def test_t4_a_deleted_table_is_refused(
        self, cm: ModuleType, tmp_path: Path, stem: str, table: str, word: str
    ) -> None:
        """`[edges]` and `[vocabulary]` go: the edge names follow from the class."""
        refused(cm, write(tmp_path / f"{stem}.toml", coded_map() + "\n" + table), word)

    @pytest.mark.parametrize("name", ["all", "water", "river"])
    def test_t4_a_category_named_as_a_selector_is_refused(
        self, cm: ModuleType, tmp_path: Path, name: str
    ) -> None:
        """So an `--edges` selector means one thing: `all`, `water` and the
        line keys are no category names."""
        classes = f'[[class]]\ncode = 1\nlabel = "F"\ncolour = "#1f8d49"\ncategory = "{name}"\n'
        refused(cm, write(tmp_path / f"category-{name}.toml", coded_map(classes=classes)), name)

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
        for word in ("corin", *SHIPPED):
            assert word in str(caught.value), word

    @pytest.mark.parametrize("value", ["./nope.toml", "nowhere/nope"])
    def test_t5_a_missing_path_is_refused_naming_it(
        self, cm: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, value: str
    ) -> None:
        monkeypatch.chdir(tmp_path)
        with pytest.raises(ValueError, match="nope"):
            cm.load(value)

    @pytest.mark.parametrize("name", ["property", "corine-water"])
    def test_t5_property_and_corine_water_are_unknown_names(
        self, cm: ModuleType, name: str
    ) -> None:
        with pytest.raises(ValueError) as caught:
            cm.load(name)
        assert name in str(caught.value)
        for word in SHIPPED:
            assert word in str(caught.value), word


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


def cover_source(tmp_path: Path, unlisted: str, *, bog: bool = True) -> tuple[Path, str]:
    """Forest and Lake by text, and a Bog the `[codes]` table lacks; the
    data file and the class file's path."""
    classes = write(tmp_path / f"cover-{unlisted}.toml", coded_map(
        unlisted=unlisted, codes='Forest = 1\nLake = 2\n'
    ))  # fmt: skip
    features = [Feat("west", WEST, {"cover": "Forest"}), Feat("east", EAST, {"cover": "Lake"})]
    if bog:
        features.append(Feat("bog-7", square(20, 60, 40, 80), {"cover": "Bog"}))
    return write_geojson(tmp_path / f"cover-{unlisted}.geojson", features), str(classes)


def cm_load(path: Path) -> Any:
    return importlib.import_module("tin_engine.land_cover_classes").load(str(path))


class TestTextCodes:
    """T8 (B6): `[codes]` turns attribute text into class codes."""

    def test_t8_the_codes_table_is_kept(self, tmp_path: Path) -> None:
        the_legend = cm_load(write(tmp_path / "c.toml", coded_map(codes="Forest = 1\nLake = 2\n")))
        assert dict(the_legend.string_codes) == {"Forest": 1, "Lake": 2}
        assert the_legend.classes == frozenset({"1", "2", "Forest", "Lake"})
        assert the_legend.water_codes == frozenset({2})

    def test_t8_text_values_carry_their_codes_and_the_lake_water(self, tmp_path: Path) -> None:
        path, classes = cover_source(tmp_path, "drop", bog=False)
        fs = ff.open_one(path, DOMAIN, classes)
        got = {f.fid: (f.code, f.mask) for f in fs.features}
        assert got == {0: (1, LAND_COVER), 1: (2, LAND_COVER | WATER)}

    @pytest.mark.parametrize("unlisted", ["refuse", "keep"])
    def test_t8_an_unlisted_text_is_refused_naming_file_feature_and_value(
        self, tmp_path: Path, unlisted: str
    ) -> None:
        fi = ff.feature_input()
        path, classes = cover_source(tmp_path, unlisted)
        with pytest.raises(fi.FeatureError) as caught:
            ff.open_one(path, DOMAIN, classes)
        for word in (f"cover-{unlisted}.geojson", "bog-7", "Bog"):
            assert word in str(caught.value), word

    def test_t8_an_unlisted_text_is_dropped_with_drop(self, tmp_path: Path) -> None:
        path, classes = cover_source(tmp_path, "drop")
        fs = ff.open_one(path, DOMAIN, classes)
        assert sorted(f.code for f in fs.features) == [1, 2]

    def test_t8_without_a_codes_table_a_text_is_refused_as_today(self, tmp_path: Path) -> None:
        fi = ff.feature_input()
        classes = write(tmp_path / "plain.toml", coded_map(unlisted="keep"))
        path = write_geojson(tmp_path / "plain.geojson", [Feat("f", WEST, {"cover": "Forest"})])
        with pytest.raises(fi.FeatureError, match="Forest"):
            ff.open_one(path, DOMAIN, str(classes))


# ------------------------------------------------------------ T10, T11


class TestMapBiomas:
    def test_t10_the_design_anchors(self, cm: ModuleType) -> None:
        palette = cm.load("mapbiomas-c11").legend.colours
        assert palette[3] == ("Forest Formation", "#1f8d49")
        assert palette[4] == ("Savanna Formation", "#7dc975")
        assert palette[5] == ("Mangrove", "#04381d")
        assert palette[33] == ("River, Lake and Ocean", "#2532e4")
        assert palette[31] == ("Aquaculture", "#091077")

    def test_t10_every_class_is_the_legend_csvs(self, cm: ModuleType) -> None:
        palette = dict(cm.load("mapbiomas-c11").legend.colours)
        assert len(mapbiomas_csv()) == 33
        assert {code: v for code, v in palette.items() if code} == mapbiomas_csv()

    def test_t10_no_class_water_and_the_notice(self, cm: ModuleType) -> None:
        got = cm.load("mapbiomas-c11")
        assert got.legend.colours[0] == ("No class; constraint lines", "#3c3c3c")
        assert got.water_codes == frozenset({33, 31})
        assert "26" not in got.classes
        assert {"33", "3"} <= got.classes
        assert got.otherwise == "keep"
        assert got.notice == MAPBIOMAS_NOTICE
        assert (got.attribute, got.accuracy_m) == ("class", 30.0)
        assert got.codes
        assert got.licence == ""

    def test_t10_the_catalogue_entry_names_the_map(self) -> None:
        sources = importlib.import_module("tin_engine.sources")
        assert sources.SOURCES["mapbiomas-c11"].classes == "mapbiomas-c11"


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
#: `docs/increments/91-ar50-vegetation-cover.md`, section 3's two tables:
#: parent, value, code, label, colour, and NIBIO's RGB (None: NIBIO gives none).
AR50_SPLIT_TABLE = (
    (
        30,
        31,
        3031,
        "Forest, coniferous: at least 50 % of the forested area covered by conifers",
        "#895944",
        (137, 89, 68),
    ),
    (
        30,
        32,
        3032,
        "Forest, deciduous: less than 20 % of the forested area covered by conifers",
        "#8b786f",
        (139, 120, 111),
    ),
    (
        30,
        33,
        3033,
        "Forest, mixed: 20 to 50 % of the forested area covered by conifers",
        "#8b786e",
        (139, 120, 110),
    ),
    (30, 99, 3099, "Forest, tree type not registered", "#9ecc73", None),
    (50, 51, 5051, "Open firm ground, no vegetation", "#c0c0c0", (192, 192, 192)),
    (
        50,
        52,
        5052,
        "Open firm ground, barren to non-continuous vegetation cover",
        "#ffe8c2",
        (255, 232, 194),
    ),
    (50, 54, 5054, "Open firm ground, barren to good vegetation cover", "#e1c790", (225, 199, 144)),
    (
        50,
        55,
        5055,
        "Open firm ground, good to lush non-continuous vegetation cover",
        "#c7a971",
        (199, 169, 113),
    ),
    (50, 99, 5099, "Open firm ground, vegetation cover not registered", "#f8f7ba", (248, 247, 186)),
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
        assert got.legend.preset == "rasputin AR50 NIBIO"

    def test_the_palette_is_exactly_the_tables(self, cm: ModuleType) -> None:
        """The nine classes, then the four forest and five open-firm-ground
        sub-classes of `docs/increments/91-ar50-vegetation-cover.md`, section 3."""
        want = {0: ("No polygon; constraint lines", "#3c3c3c")}
        want |= {code: (label, colour) for code, label, colour, _ in AR50_TABLE}
        want |= {code: (label, colour) for _, _, code, label, colour, _ in AR50_SPLIT_TABLE}
        assert dict(cm.load("ar50").legend.colours) == want
        assert len(want) == 19

    def test_the_sub_class_colours_are_nibios_rgb(self) -> None:
        """Each sub-class colour is its RGB in hex; forest 99, which NIBIO
        does not colour, takes forest's own."""
        for parent, value, code, _, colour, rgb in AR50_SPLIT_TABLE:
            if rgb is None:
                assert (code, colour) == (3099, "#9ecc73")
            else:
                assert colour == "#{:02x}{:02x}{:02x}".format(*rgb), code
            assert code == 100 * parent + value, code

    def test_the_splits(self, cm: ModuleType) -> None:
        """Forest by tree type, open firm ground by vegetation cover, in code order."""
        got = [(s.code, s.attribute, dict(s.codes)) for s in cm.load("ar50").splits]
        assert got == [
            (30, "treslag", {"31": 3031, "32": 3032, "33": 3033, "99": 3099}),
            (
                50,
                "vegetasjonsdekke",
                {"51": 5051, "52": 5052, "54": 5054, "55": 5055, "99": 5099},
            ),
        ]

    def test_fresh_water_and_ocean_carry_water_and_the_rest_land_cover(
        self, cm: ModuleType
    ) -> None:
        got = cm.load("ar50")
        assert got.classes == frozenset(str(code) for code, *_ in AR50_TABLE)
        assert got.water_codes == frozenset({81, 82})
        assert got.otherwise == "keep"


# ------------------------------------------- splits by a second attribute


def split_classes(*splits: str) -> str:
    """Classes 10, 30 and 50, then ``splits`` (each one or more TOML tables)."""
    classes = "".join(
        f'[[class]]\ncode = {code}\nlabel = "{label}"\ncolour = "{colour}"\n\n'
        for code, label, colour in (
            (10, "Built", "#e86474"),
            (30, "Forest", "#9ecc73"),
            (50, "Open", "#d9d9d9"),
        )
    )
    return classes + "\n".join(splits)


def split(code: int, attribute: str, *subs: tuple[Any, int]) -> str:
    """One `[[split]]` of ``code`` by ``attribute``; ``subs`` are (value, code),
    each labelled and coloured from its code."""
    text = f'[[split]]\ncode = {code}\nattribute = "{attribute}"\n\n'
    for value, sub in subs:
        shown = f'"{value}"' if isinstance(value, str) else str(value)
        text += (
            f"[[split.class]]\nvalue = {shown}\ncode = {sub}\n"
            f'label = "Sub {sub}"\ncolour = "#0000{sub % 100:02d}"\n\n'
        )
    return text


#: Section 8, test 1: 50 by `veg` listed before 30 by `tre`, 99 in both.
VEG = split(50, "veg", (51, 5051), (99, 5099))
TRE = split(30, "tre", (31, 3031), (99, 3099))


class TestSplitLoads:
    """`docs/increments/91-ar50-vegetation-cover.md`, section 8, test 1."""

    def test_the_splits_in_code_order_and_their_codes(self, cm: ModuleType, tmp_path: Path) -> None:
        text = coded_map(classes=split_classes(VEG, TRE))
        got = cm.load(str(write(tmp_path / "two.toml", text)))
        assert [(s.code, s.attribute, dict(s.codes)) for s in got.splits] == [
            (30, "tre", {"31": 3031, "99": 3099}),
            (50, "veg", {"51": 5051, "99": 5099}),
        ]

    def test_the_sub_classes_join_the_legend_not_the_classes(
        self, cm: ModuleType, tmp_path: Path
    ) -> None:
        got = cm.load(str(write(tmp_path / "two.toml", coded_map(classes=split_classes(VEG, TRE)))))
        assert dict(got.legend.colours) == {
            0: ("No polygon; constraint lines", "#3c3c3c"),
            10: ("Built", "#e86474"),
            30: ("Forest", "#9ecc73"),
            50: ("Open", "#d9d9d9"),
            3031: ("Sub 3031", "#000031"),
            3099: ("Sub 3099", "#000099"),
            5051: ("Sub 5051", "#000051"),
            5099: ("Sub 5099", "#000099"),
        }
        assert set(got.classes) == {"10", "30", "50"}

    def test_a_legend_without_splits_has_none(self, cm: ModuleType, tmp_path: Path) -> None:
        assert cm.load(str(write(tmp_path / "base.toml", coded_map()))).splits == ()

    @pytest.mark.parametrize("value", ["52", 52, "52.0", " 52 "])
    def test_a_value_is_keyed_as_a_whole_number(
        self, cm: ModuleType, tmp_path: Path, value: Any
    ) -> None:
        """Section 3, the key of a value: the loader keys a sub-class's
        `value` as a feature's value is keyed at run time."""
        text = coded_map(classes=split_classes(split(50, "veg", (value, 5052))))
        (got,) = cm.load(str(write(tmp_path / "key.toml", text))).splits
        assert dict(got.codes) == {"52": 5052}

    def test_a_text_value_with_letters_is_kept(self, cm: ModuleType, tmp_path: Path) -> None:
        text = coded_map(classes=split_classes(split(50, "veg", ("5a", 5052))))
        (got,) = cm.load(str(write(tmp_path / "letters.toml", text))).splits
        assert dict(got.codes) == {"5a": 5052}


def bare(text: str) -> str:
    """``text`` less its label and colour lines, which a rules file refuses."""
    return re.sub(r'^(label|colour) = ".*"\n', "", text, flags=re.MULTILINE)


def coded_with(*splits: str) -> str:
    """`coded_map` over `split_classes` with ``splits``."""
    return coded_map(classes=split_classes(*splits))


#: `VEG` with its first sub-class given a key a sub-class may not carry.
WET = VEG.replace("code = 5051\n", "code = 5051\nwater = true\n")
LINED = VEG.replace("code = 5051\n", 'code = 5051\nedges = ["levee"]\n')
#: (file stem, the legend text, words the refusal names besides the file).
SPLIT_REFUSALS = [
    ("unlisted-class", coded_with(split(40, "veg", (51, 4051))), ("40",)),
    ("class-twice", coded_with(TRE, split(30, "veg", (51, 3051))), ("30",)),
    ("first-attribute", coded_with(split(50, "COVER", (51, 5051))), ("COVER",)),
    ("sub-is-wet", coded_with(WET), ("water",)),
    ("sub-lines", coded_with(LINED), ("edges",)),
    ("one-key", coded_with(split(50, "veg", (51, 5051), ("51.0", 5052))), ("51",)),
    ("listed-code", coded_with(split(50, "veg", (51, 10))), ("10",)),
    ("code-twice", coded_with(split(50, "veg", (51, 5051), (52, 5051))), ("5051",)),
    ("code-in-two", coded_with(TRE, split(50, "veg", (51, 3031))), ("3031",)),
]


class TestSplitRefusals:
    """Section 8, test 2: section 3's rules 1 to 5, each refused naming the file."""

    @pytest.mark.parametrize(
        ("stem", "text", "words"), SPLIT_REFUSALS, ids=[c[0] for c in SPLIT_REFUSALS]
    )
    def test_refused(
        self, cm: ModuleType, tmp_path: Path, stem: str, text: str, words: tuple[str, ...]
    ) -> None:
        refused(cm, write(tmp_path / f"{stem}.toml", text), *words)

    def test_two_splits_may_share_an_attribute_and_a_value(
        self, cm: ModuleType, tmp_path: Path
    ) -> None:
        """Rules 3 and 5 allow it: the refusals above are their changes alone."""
        text = coded_map(
            classes=split_classes(split(50, "veg", (99, 5099)), split(30, "veg", (99, 3099)))
        )
        got = cm.load(str(write(tmp_path / "shared.toml", text)))
        assert [(s.code, s.attribute) for s in got.splits] == [(30, "veg"), (50, "veg")]


# ------------------------------------- the rules file apart from its legend


#: 108's snapshot of the five shipped names at `5de1144d` (module docstring).
BEFORE: dict[str, Any] = json.loads((DATA / "legends_5de1144d.json").read_text("utf-8"))
#: (file stem, the line changed in `coded_map()`, its replacement, words named).
RULES_ONE_LINE = [
    ("label-in-rules", "code = 1\n", 'code = 1\nlabel = "Forest"\n', ("label", "legend")),
    ("colour-in-rules", "code = 1\n", 'code = 1\ncolour = "#1f8d49"\n', ("colour", "legend")),
    ("preset-in-rules", 'unlisted = "keep"\n', 'unlisted = "keep"\npreset = "p"\n',
     ("preset", "legend")),
    ("category-in-class", "code = 1\n", 'code = 1\ncategory = "forest"\n', ("category", "[cells]")),
]  # fmt: skip


def without_legend(text: str) -> str:
    """``text`` less its `legend` line."""
    return re.sub(r'^legend = ".*"\n', "", text, flags=re.MULTILINE)


class TestRulesRefusals:
    """108, test 1: each refusal names the rules file and what is wrong."""

    @pytest.mark.parametrize(
        ("stem", "old", "new", "words"), RULES_ONE_LINE, ids=[c[0] for c in RULES_ONE_LINE]
    )
    def test_a_legend_key_or_a_cells_key_in_its_old_place(
        self, cm: ModuleType, tmp_path: Path, stem: str, old: str, new: str, words: tuple[str, ...]
    ) -> None:
        text = coded_map()
        assert old in text, "the fixture changed under this test"
        refused(cm, write(tmp_path / f"{stem}.toml", text.replace(old, new, 1)), *words)

    def test_a_no_class_table_in_rules(self, cm: ModuleType, tmp_path: Path) -> None:
        text = coded_map() + '\n[no_class]\nlabel = "None"\ncolour = "#3c3c3c"\n'
        refused(cm, write(tmp_path / "no-class.toml", text), "no_class", "legend")

    def test_a_coded_file_without_a_legend(self, cm: ModuleType, tmp_path: Path) -> None:
        refused(cm, write(tmp_path / "bare.toml", without_legend(coded_map())), "legend")

    def test_a_class_code_missing_from_the_legend(self, cm: ModuleType, tmp_path: Path) -> None:
        colours = write(tmp_path / "forest-only.toml", legend_text("f", "p", [(1, "F", "#1f8d49")]))
        text = re.sub(r'^legend = ".*"$', f'legend = "{colours.name}"', coded_map(), flags=re.M)
        message = refused(cm, write(tmp_path / "lake-missing.toml", text))
        assert re.search(r"\b2\b", message), message

    def test_a_sub_class_code_missing_from_the_legend(self, cm: ModuleType, tmp_path: Path) -> None:
        listed = [(10, "B", "#e86474"), (30, "F", "#9ecc73"), (50, "O", "#d9d9d9")]
        colours = write(
            tmp_path / "parents.toml", legend_text("p", "p", [*listed, (5051, "S", "#000051")])
        )
        text = coded_with(VEG)
        text = re.sub(r'^legend = ".*"$', f'legend = "{colours.name}"', text, flags=re.M)
        message = refused(cm, write(tmp_path / "sub-missing.toml", text))
        assert "5099" in message and "5051" not in message, message

    def test_a_top_level_prefer(self, cm: ModuleType, tmp_path: Path) -> None:
        text = coded_map() + "\n[prefer]\nlake = 2\n"
        refused(cm, write(tmp_path / "prefer.toml", text), "prefer", "[cells]")

    def test_a_code_in_two_categories(self, cm: ModuleType, tmp_path: Path) -> None:
        classes = '[[class]]\ncode = 1\nlabel = "F"\ncolour = "#1f8d49"\ncategory = "forest"\n'
        text = coded_map(classes=classes).replace("forest = [1]", "forest = [1]\nlake = [1]")
        message = refused(cm, write(tmp_path / "twice.toml", text))
        assert re.search(r"\b1\b", message), message

    def test_a_cells_code_not_listed(self, cm: ModuleType, tmp_path: Path) -> None:
        classes = '[[class]]\ncode = 1\nlabel = "F"\ncolour = "#1f8d49"\ncategory = "forest"\n'
        text = coded_map(classes=classes).replace("forest = [1]", "forest = [1, 9]")
        message = refused(cm, write(tmp_path / "unlisted-cell.toml", text))
        assert re.search(r"\b9\b", message), message


class TestRulesNameTheirLegend:
    def test_a_shipped_legend_by_name(self, cm: ModuleType, tmp_path: Path) -> None:
        text = re.sub(r'^legend = ".*"$', 'legend = "corine"', corine_like("w", "40"), flags=re.M)
        got = cm.load(str(write(tmp_path / "w.toml", text)))
        assert dict(got.legend.colours) == dict(cm.load("corine").legend.colours)
        assert got.legend.name == "corine"

    def test_a_legend_path_is_relative_to_the_rules_file(
        self, cm: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        colours = write(
            tmp_path / "colours" / "mine.toml",
            legend_text(
                "mine-colours", "rasputin mine", [(1, "F", "#1f8d49"), (2, "L", "#2532e4")]
            ),
        )
        text = re.sub(
            r'^legend = ".*"$', 'legend = "../colours/mine.toml"', coded_map(), flags=re.M
        )
        rules = write(tmp_path / "rules" / "mine.toml", text)
        monkeypatch.chdir(tmp_path / "colours")
        got = cm.load(str(rules))
        assert got.legend.name == "mine-colours"
        assert dict(got.legend.colours)[2] == ("L", "#2532e4")
        assert colours.exists()

    def test_cells_categories_and_prefer(self, cm: ModuleType, tmp_path: Path) -> None:
        """`[cells]` holds the categories and `prefer`; a split's sub-classes
        take their parent's category."""
        classes = (
            split_classes(VEG).replace(
                'code = 50\nlabel = "Open"', 'code = 50\nlabel = "Open"\ncategory = "forest"'
            )
            + '\n[[class]]\ncode = 2\nlabel = "Lake"\ncolour = "#2532e4"\nwater = true\n'
            + 'category = "lake"\n'
        )
        got = cm.load(
            str(write(tmp_path / "cells.toml", coded_map(classes=classes, prefer="lake = 2\n")))
        )
        assert dict(got.categories) == {50: "forest", 5051: "forest", 5099: "forest", 2: "lake"}
        assert dict(got.prefer) == {"lake": 2}
        assert got.water_codes == frozenset({2})


class TestSameColoursAndCells:
    """108, test 3: each shipped name's legend and cell mappings equal `5de1144d`'s."""

    def test_the_snapshot_holds_the_five(self) -> None:
        """`corine-water` goes; the other four are checked against it."""
        assert set(BEFORE) == {*SHIPPED, "corine-water"}

    @pytest.mark.parametrize("name", SHIPPED)
    def test_colours_preset_categories_prefer_and_water(self, cm: ModuleType, name: str) -> None:
        want, got = BEFORE[name], cm.load(name)
        assert [[c, *v] for c, v in got.legend.colours.items()] == want["colours"]
        assert got.legend.preset == want["preset"]
        assert {str(c): str(v) for c, v in got.categories.items()} == want["categories"]
        assert {str(k): v for k, v in got.prefer.items()} == want["prefer"]
        assert sorted(got.water_codes) == want["water_codes"]

    @pytest.mark.parametrize("name", SHIPPED)
    def test_the_rules_file_holds_no_colour(self, name: str) -> None:
        text = (shipped_dir() / f"{name}.toml").read_text(encoding="utf-8")
        for key in ("label", "colour", "preset", "no_class", "category"):
            assert not re.search(rf"^\[?{key}\b", text, flags=re.MULTILINE), key


class TestKnownCodes:
    """108, test 4: the codes a class polygon may carry are the rules' own."""

    def test_ar50_knows_its_sub_classes(self, cm: ModuleType) -> None:
        got = cm.known_codes(cm.load("ar50"))
        assert {3031, 3032, 3033, 3099, 5051, 5052, 5054, 5055, 5099} <= got

    @pytest.mark.parametrize("name", SHIPPED)
    def test_known_codes_are_the_old_palette_keys(self, cm: ModuleType, name: str) -> None:
        want = {code for code, *_ in BEFORE[name]["colours"]}
        assert cm.known_codes(cm.load(name)) == frozenset(want)

    def test_zero_is_the_rules_setting_not_the_legends(
        self, cm: ModuleType, tmp_path: Path
    ) -> None:
        """Class code 0 (no polygon, constraint lines only) is accepted when
        the rules say ``accept_no_class = true``; the legend's `[no_class]`
        is its colour and label only, and decides nothing."""
        off = cm.load(str(write(tmp_path / "r.toml", coded_map())))
        (named,) = re.findall(r'^legend = "(.*)"$', coded_map(), flags=re.M)
        assert "[no_class]" in (tmp_path / named).read_text(encoding="utf-8")
        on = cm.load(str(write(tmp_path / "r0.toml", coded_map(accept_no_class=True))))
        assert cm.known_codes(off) == {1, 2}
        assert cm.known_codes(on) == {0, 1, 2}
        assert 0 in off.legend.colours and 0 in on.legend.colours

    @pytest.mark.parametrize("name", SHIPPED)
    def test_every_shipped_rules_file_accepts_code_0(self, name: str) -> None:
        assert shipped_toml(name).get("accept_no_class") is True, name

    def test_accepting_zero_needs_its_colour_in_the_legend(
        self, cm: ModuleType, tmp_path: Path
    ) -> None:
        listed = [(1, "F", "#1f8d49"), (2, "L", "#2532e4")]
        colours = write(tmp_path / "plain.toml", legend_text("p", "p", listed, no_class=False))
        text = re.sub(
            r'^legend = ".*"$', f'legend = "{colours.name}"', coded_map(accept_no_class=True),
            flags=re.M,
        )  # fmt: skip
        message = refused(cm, write(tmp_path / "zero-uncoloured.toml", text), "legend")
        assert re.search(r"\b0\b", message), message
