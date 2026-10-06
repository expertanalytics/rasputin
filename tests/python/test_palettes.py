"""`tin_engine.palettes` and `rasputin palette`: natural colours for CORINE (16c, R4).

`docs/increments/16c-landcover-labels.md`, R4 and "Tests for @tester". The
colours are Ola's to change (D6), so they are pinned by family, as hue,
saturation and lightness ranges, never as exact triples: forest green, heath a
lighter green than conifers, bare rock grey, peat bog brown, every water class
blue, glaciers near white and not warm.

Pinned beyond the design's text:

- An `Annotations` value is compared as `int(value)`, so the preset may carry
  codes as strings (ParaView's own presets do) or as numbers.
- Each annotation's label contains the class's `LABEL3` name.
- `NanColor` is present and magenta (R4: "A code the table lacks draws in
  `NanColor`, magenta").

Committed red at `196147e`: `tin_engine.palettes` did not exist yet and
`rasputin palette` was no command, so every test failed on
`ModuleNotFoundError` or on the CLI's usage error for an unknown command. Both
landed in `0487ed0` and the suite has been green since.
"""

from __future__ import annotations

import colorsys
import json
import re
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any

import pytest

from cli_helpers import runner
from gpkg_fixtures import OLA_NORWAY
from tin_engine import palettes
from tin_engine.cli import app

#: R4's "in the Norway extract" column: 34 codes.
NORWAY = frozenset(
    {
        *(111, 112, 121, 122, 123, 124, 131, 132, 133, 141, 142),
        *(211, 222, 231, 242, 243),
        *(311, 312, 313, 321, 322, 324, 331, 332, 333, 334, 335),
        *(411, 412, 423),
        *(511, 512, 522, 523),
    }
)

#: The 44 CLC level-3 classes (Kosztra et al. 2019).
CLC = frozenset(
    {
        *(111, 112, 121, 122, 123, 124, 131, 132, 133, 141, 142),
        *(211, 212, 213, 221, 222, 223, 231, 241, 242, 243, 244),
        *(311, 312, 313, 321, 322, 323, 324, 331, 332, 333, 334, 335),
        *(411, 412, 421, 422, 423),
        *(511, 512, 521, 522, 523),
    }
)

PRESET_NAME = "rasputin CORINE natural"
HEX = re.compile(r"#[0-9a-fA-F]{6}")


def table() -> dict[int, tuple[str, str]]:
    return dict(palettes.CORINE_NATURAL)


def rgb(code: int) -> tuple[float, float, float]:
    """The table's colour for `code`, as r, g, b in 0 .. 1."""
    text = table()[code][1]
    assert HEX.fullmatch(text), text
    r, g, b = (int(text[i : i + 2], 16) / 255 for i in (1, 3, 5))
    return r, g, b


def hls(code: int) -> tuple[float, float, float]:
    """Hue in degrees, lightness and saturation (HLS), each 0 .. 1 but hue."""
    h, lightness, s = colorsys.rgb_to_hls(*rgb(code))
    return h * 360, lightness, s


# --------------------------------------------------------------- the table


class TestTheTable:
    def test_every_clc_class_and_zero_has_a_colour(self) -> None:
        assert set(table()) == CLC | {0}

    def test_every_norway_code_has_a_colour(self) -> None:
        assert len(NORWAY) == 34 and NORWAY <= CLC and len(CLC) == 44
        assert NORWAY - set(table()) == set()

    def test_entries_are_a_name_and_a_hex_colour(self) -> None:
        for code, (name, colour) in table().items():
            assert isinstance(name, str) and name, code
            assert HEX.fullmatch(colour), (code, colour)

    def test_the_names_are_label3(self) -> None:
        # Spot checks from `Legend/CLC_legend.csv` as R4 quotes it.
        assert table()[312][0] == "Coniferous forest"
        assert table()[412][0] == "Peat bogs"
        assert table()[335][0] == "Glaciers and perpetual snow"
        assert table()[512][0] == "Water bodies"

    def test_the_colours_are_distinct(self) -> None:
        colours = [c.lower() for _, c in table().values()]
        assert len(set(colours)) == len(colours)


class TestFamilies:
    """D6: the families, by hue ranges, not by exact triple."""

    @pytest.mark.parametrize("code", [311, 312, 313])
    def test_forests_are_green(self, code: int) -> None:
        hue, _, saturation = hls(code)
        assert 75 <= hue <= 165 and saturation >= 0.2, hls(code)

    def test_heath_is_green_and_lighter_than_conifers(self) -> None:
        hue, lightness, saturation = hls(322)
        assert 60 <= hue <= 150 and saturation >= 0.15, hls(322)
        assert lightness > hls(312)[1]

    def test_bare_rock_is_grey(self) -> None:
        _, lightness, saturation = hls(332)
        assert saturation <= 0.08 and 0.25 <= lightness <= 0.75, hls(332)

    def test_peat_bogs_are_brown(self) -> None:
        hue, lightness, saturation = hls(412)
        assert 15 <= hue <= 50 and saturation >= 0.2 and lightness <= 0.55, hls(412)

    @pytest.mark.parametrize("code", sorted(c for c in CLC if c // 100 == 5))
    def test_water_is_blue(self, code: int) -> None:
        hue, _, saturation = hls(code)
        assert 190 <= hue <= 235 and saturation >= 0.25, hls(code)

    def test_glaciers_are_near_white_and_not_warm(self) -> None:
        r, _, b = rgb(335)
        _, lightness, _ = hls(335)
        assert lightness >= 0.9 and b >= r, rgb(335)

    def test_the_official_legend_is_not_used_for_bogs_or_the_sea(self) -> None:
        """Direction 2: the official legend paints peat bogs blue (077-077-255)
        and the sea almost white (230-242-255)."""
        assert table()[412][1].lower() != "#4d4dff"
        assert table()[523][1].lower() != "#e6f2ff"


# ---------------------------------------------------------- the preset


class TestPreset:
    def preset(self) -> dict[str, Any]:
        out = palettes.paraview_preset(palettes.CORINE_NATURAL, "a name")
        assert isinstance(out, list) and len(out) == 1
        assert isinstance(out[0], dict)
        return out[0]

    def test_it_carries_the_name(self) -> None:
        preset = self.preset()
        assert preset["Name"] == "a name"

    def test_annotations_and_colours_have_one_entry_per_code(self) -> None:
        preset = self.preset()
        entries = len(table())
        assert len(preset["Annotations"]) == 2 * entries
        assert len(preset["IndexedColors"]) == 3 * entries
        assert all(0.0 <= float(c) <= 1.0 for c in preset["IndexedColors"])

    def test_they_are_in_the_same_order(self) -> None:
        preset = self.preset()
        values = [int(v) for v in preset["Annotations"][0::2]]
        labels = [str(v) for v in preset["Annotations"][1::2]]
        colours = preset["IndexedColors"]
        assert sorted(values) == sorted(table())
        for k, code in enumerate(values):
            assert table()[code][0] in labels[k], (code, labels[k])
            got = [float(c) for c in colours[3 * k : 3 * k + 3]]
            assert got == pytest.approx(rgb(code), abs=1e-6), code

    def test_the_nan_colour_is_magenta(self) -> None:
        preset = self.preset()
        r, g, b = (float(c) for c in preset["NanColor"])
        assert r >= 0.8 and b >= 0.8 and g <= 0.3, preset["NanColor"]

    def test_it_is_json(self) -> None:
        preset = self.preset()
        assert json.loads(json.dumps(preset)) == preset


# ------------------------------------------------------------- the command


class TestCommand:
    def test_stdout_is_the_preset(self) -> None:
        result = runner.invoke(app, ["palette", "corine"])
        assert result.exit_code == 0, result.output
        expected = palettes.paraview_preset(palettes.CORINE_NATURAL, PRESET_NAME)
        assert json.loads(result.stdout) == expected

    def test_an_unknown_name_is_a_usage_error_listing_the_known(self) -> None:
        result = runner.invoke(app, ["palette", "nosuch"])
        assert result.exit_code == 2, result.output
        assert "nosuch" in result.output and "corine" in result.output
        assert "No such command" not in result.output


# --------------------------------------------------------- Ola's local data


@pytest.mark.skipif(not OLA_NORWAY.exists(), reason=f"{OLA_NORWAY} is not here")
def test_every_code_in_olas_norway_extract_has_a_colour() -> None:
    """R4's query, rerun: every `code_18` in Ola's Norway GeoPackage."""
    uri = Path(OLA_NORWAY).resolve().as_uri() + "?mode=ro"
    with closing(sqlite3.connect(uri, uri=True)) as con:
        codes = {int(c) for (c,) in con.execute("SELECT DISTINCT code_18 FROM corine2018")}
    assert codes == NORWAY
    assert codes <= set(table())
