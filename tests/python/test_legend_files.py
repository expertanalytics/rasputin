"""`tin_engine.legends`: a legend file is a source's colour scheme, and only that.

`docs/increments/108-legend-split.md`, section 4 (the legend file, the shipped
legends, `rasputin palette`) and section 5, test 2. The meshing-rules half is
in `test_legends.py`; the bytes of `rasputin palette corine` in
`test_cli_legends.TestPalette`.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `legends.load_legend(value)` takes a shipped legend name or a path (ends in
  `.toml` or holds a path separator), and returns the frozen `Legend` with
  `name`, `preset` and `colours` (code to (label, colour), code 0 first when
  `[no_class]` is given, then the classes in file order). Every refusal is a
  `ValueError` naming the file and, for a bad key or value, the key.
- The shipped legends are exactly `corine`, `mapbiomas-c11` and `ar50`, in
  `tin_engine/data/legends/`.
- `rasputin palette` given a meshing-rules file is refused (exit 2) naming the
  file.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import pytest
from pydantic import ValidationError

from cli_helpers import USAGE, invoke
from legend_fixtures import basenames, coded_map, legend_text, write
from tin_engine import legends

LEGENDS = ("ar50", "corine", "mapbiomas-c11")
BASE = [(1, "Forest", "#1f8d49"), (2, "Lake", "#2532e4")]


def refused(path: Path, *words: str) -> None:
    """`load_legend(path)` refused naming the file and every one of ``words``."""
    with pytest.raises(ValueError) as caught:
        legends.load_legend(str(path))
    message = basenames(str(caught.value))
    for word in (path.name, *words):
        assert word in message, f"{word!r} not in: {message}"


class TestTheLegendFile:
    def test_the_base_loads(self, tmp_path: Path) -> None:
        got = legends.load_legend(str(write(tmp_path / "c.toml", legend_text("c", "p c", BASE))))
        assert (got.name, got.preset) == ("c", "p c")
        assert list(got.colours.items()) == [
            (0, ("No polygon; constraint lines", "#3c3c3c")),
            (1, ("Forest", "#1f8d49")),
            (2, ("Lake", "#2532e4")),
        ]

    def test_no_class_is_optional(self, tmp_path: Path) -> None:
        text = legend_text("c", "p", BASE, no_class=False)
        assert 0 not in legends.load_legend(str(write(tmp_path / "c.toml", text))).colours

    def test_the_legend_is_frozen(self, tmp_path: Path) -> None:
        got = legends.load_legend(str(write(tmp_path / "c.toml", legend_text("c", "p", BASE))))
        with pytest.raises((TypeError, ValueError, ValidationError)):
            got.name = "other"


#: (file stem, the line changed in the base legend, its replacement, words named).
ONE_LINE = [
    ("extra-key", 'preset = "p"\n', 'preset = "p"\nattribute = "x"\n', ("attribute",)),
    ("extra-class-key", "code = 1\n", "code = 1\nwater = true\n", ("water",)),
    ("no-preset", 'preset = "p"\n', "", ("preset",)),
    ("no-name", 'name = "c"\n', "", ("name",)),
    ("no-label", 'label = "Forest"\n', "", ("label",)),
    ("repeated-code", "code = 2\n", "code = 1\n", ("code",)),
    ("code-zero", "code = 1\n", "code = 0\n", ("code",)),
    ("colour-word", 'colour = "#1f8d49"', 'colour = "green"', ("colour",)),
    ("colour-upper", 'colour = "#1f8d49"', 'colour = "#1F8D49"', ("colour",)),
    ("colour-short", 'colour = "#1f8d49"', 'colour = "#1f8d4"', ("colour",)),
]


class TestLegendRefusals:
    """108, test 2: every key required, nothing else allowed."""

    @pytest.mark.parametrize(
        ("stem", "old", "new", "words"), ONE_LINE, ids=[c[0] for c in ONE_LINE]
    )
    def test_one_line_changed(
        self, tmp_path: Path, stem: str, old: str, new: str, words: tuple[str, ...]
    ) -> None:
        text = legend_text("c", "p", BASE)
        assert old in text, "the fixture changed under this test"
        refused(write(tmp_path / f"{stem}.toml", text.replace(old, new, 1)), *words)

    def test_a_rules_file_is_no_legend(self, tmp_path: Path) -> None:
        refused(write(tmp_path / "rules.toml", coded_map()), "attribute")

    def test_an_unknown_name_lists_the_shipped_legends(self) -> None:
        with pytest.raises(ValueError) as caught:
            legends.load_legend("corine-water")
        for word in ("corine-water", *LEGENDS):
            assert word in str(caught.value), word


class TestShippedLegends:
    def test_exactly_the_three(self) -> None:
        folder = resources.files("tin_engine") / "data" / "legends"
        names = {p.name for p in folder.iterdir() if p.name.endswith(".toml")}
        assert names == {f"{name}.toml" for name in LEGENDS}

    @pytest.mark.parametrize("name", LEGENDS)
    def test_each_loads_by_name(self, name: str) -> None:
        assert legends.load_legend(name).name == name


class TestPaletteTakesLegendsOnly:
    def test_a_rules_name_is_refused_naming_the_legends(self) -> None:
        code, output = invoke("palette", "corine-water")
        assert code == USAGE, output
        assert (
            "unknown legend corine-water; use ar50, corine, mapbiomas-c11, or a .toml file"
            in " ".join(output.split())
        ), output

    def test_a_legend_file(self, tmp_path: Path) -> None:
        path = write(tmp_path / "c.toml", legend_text("c", "rasputin c", BASE))
        code, output = invoke("palette", str(path))
        assert code == 0, output
        assert '"rasputin c"' in output

    def test_a_rules_file_is_refused_naming_it(self, tmp_path: Path) -> None:
        path = write(tmp_path / "rules.toml", coded_map())
        code, said = invoke("palette", str(path))
        output = basenames(said)
        assert code == USAGE, output
        assert "rules.toml" in output, output
