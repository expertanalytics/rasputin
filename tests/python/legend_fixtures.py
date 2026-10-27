"""Meshing-rules and legend files written by hand for the suites.

`docs/increments/108-legend-split.md`, *4. The design*: a legend file is the
colour scheme only (`name`, `preset`, `[no_class]`, `[[class]]` code, label,
colour); a meshing-rules file holds the rest and names its legend. The
fixtures keep the one-text shape the suites were written against: a class
given here carries its `label` and `colour` (and, optionally, `category`),
and `coded_map` sorts them into the two files. Labels and colours go to a
legend file named after a digest of its text, `category` lines to the rules
file's `[cells.categories]`, and the rules file names the legend by that
relative name. `write` writes the legend beside any rules file that names
one made here, so a suite writes one file as before.

`coded_map` is the smallest coded rules file plus the keys the shipped ones
use (`notice`, `unlisted`); `uncoded_map` the shape of the `property` rules.
Every refusal test changes one line of one of these, and
`test_legends.TestTheBaseFilesLoad` shows both load as written, so a
refusal is the change's and nothing else's.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

#: CORINE's code system, verbatim, so a user map can meet `corine` in one mesh.
CORINE_SYSTEM = "CORINE Land Cover level-3 code"
NO_CLASS = ("No polygon; constraint lines", "#3c3c3c")
#: Legend texts made by `coded_map`, by their file name.
LEGENDS: dict[str, str] = {}
_KEY = re.compile(r'^(label|colour|category) = "(.*)"$')
_LEGEND_LINE = re.compile(r'^legend = "(legend-[0-9a-f]{12}\.toml)"$', re.MULTILINE)


def legend_text(
    name: str, preset: str, colours: list[tuple[int, str, str]], no_class: bool = True
) -> str:
    """A legend file: ``colours`` as (code, label, colour), and `[no_class]`."""
    text = f'name = "{name}"\npreset = "{preset}"\n'
    if no_class:
        text += f'\n[no_class]\nlabel = "{NO_CLASS[0]}"\ncolour = "{NO_CLASS[1]}"\n'
    for code, label, colour in colours:
        text += f'\n[[class]]\ncode = {code}\nlabel = "{label}"\ncolour = "{colour}"\n'
    return text


def _sorted_out(classes: str) -> tuple[str, list[tuple[int, str, str]], dict[str, list[int]]]:
    """``classes`` less its label, colour and category lines; the colours per
    code, in order; the codes per category. A code line precedes its label."""
    kept: list[str] = []
    colours: dict[int, dict[str, str]] = {}
    categories: dict[str, list[int]] = {}
    code = 0
    for line in classes.splitlines():
        if line.startswith("[["):
            code = 0
        if line.startswith("code = "):
            code = int(line.removeprefix("code = "))
        found = _KEY.match(line)
        if found is None:
            kept.append(line)
        elif found[1] == "category":
            categories.setdefault(found[2], []).append(code)
        else:
            colours.setdefault(code, {})[found[1]] = found[2]
    listed = [(c, v["label"], v["colour"]) for c, v in colours.items() if len(v) == 2]
    return "\n".join(kept) + "\n", listed, categories


def coded_map(
    name: str = "mine",
    *,
    attribute: str = "cover",
    system: str = "Test cover code",
    accuracy: str = "40",
    unlisted: str = "land",
    classes: str = "",
    codes: str = "",
    extra: str = "",
    licence: str = "",
    prefer: str = "",
) -> str:
    """A coded rules file: Forest (code 1) and Lake (code 2, water), or
    ``classes``; its legend (preset "rasputin NAME") is kept for `write`."""
    classes = classes or (
        '[[class]]\ncode = 1\nlabel = "Forest"\ncolour = "#1f8d49"\n\n'
        '[[class]]\ncode = 2\nlabel = "Lake"\ncolour = "#2532e4"\nwater = true\n'
    )
    classes, colours, categories = _sorted_out(classes)
    legend = legend_text(f"{name}-colours", f"rasputin {name}", colours)
    file = f"legend-{hashlib.sha256(legend.encode()).hexdigest()[:12]}.toml"
    LEGENDS[file] = legend
    head = (
        f'name = "{name}"\n'
        f'attribute = "{attribute}"\n'
        f'code_system = "{system}"\n'
        f"positional_accuracy_m = {accuracy}\n"
        'notice = "Test data, made up"\n'
        f'unlisted = "{unlisted}"\n'
        f'legend = "{file}"\n'
    )
    if licence:
        head += f'licence = "{licence}"\n'
    body = '\n[edges]\nwater = ["land_cover", "water"]\nland = ["land_cover"]\n\n'
    tail = f"\n[codes]\n{codes}" if codes else ""
    tail += f"\n[vocabulary]\nextra = [{extra}]\n" if extra else ""
    if categories:
        tail += "\n[cells.categories]\n" + "".join(
            f"{category} = {listed}\n" for category, listed in categories.items()
        )
    tail += f"\n[cells.prefer]\n{prefer}" if prefer else ""
    return head + body + classes + tail


def corine_like(name: str, accuracy: str) -> str:
    """A coded map in CORINE's code system on `Code_18`, forest and water."""
    classes = (
        '[[class]]\ncode = 311\nlabel = "Broad-leaved forest"\ncolour = "#4f8f3f"\n\n'
        '[[class]]\ncode = 512\nlabel = "Water bodies"\ncolour = "#3e7bb6"\nwater = true\n'
    )
    return coded_map(
        name, attribute="Code_18", system=CORINE_SYSTEM, accuracy=accuracy, classes=classes
    )


def uncoded_map(name: str = "levees", *, extra: str = '"levee"', unlisted: str = "refuse") -> str:
    """An uncoded map on ``kind``: ``levee`` lines on the extra name, ``road`` on road."""
    vocabulary = f"\n[vocabulary]\nextra = [{extra}]\n" if extra else ""
    return (
        f'name = "{name}"\n'
        'attribute = "kind"\n'
        f'unlisted = "{unlisted}"\n'
        f"{vocabulary}\n"
        '[[class]]\nvalue = "levee"\nedges = ["levee"]\n\n'
        '[[class]]\nvalue = "road"\nedges = ["road"]\n'
    )


def write(path: Path, text: str) -> Path:
    """Write ``text`` at ``path``, and beside it the legend it names, when
    `coded_map` made that legend."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    for file in _LEGEND_LINE.findall(text):
        if file in LEGENDS:
            (path.parent / file).write_text(LEGENDS[file], encoding="utf-8")
    return path


def basenames(text: str) -> str:
    """``text`` with every path cut to its last part. A test's ``tmp_path``
    carries the test's name ("..._named_river_i0"), so a word looked for in a
    message is looked for outside the directories, where only the producer
    can have put it. Works on `plain` output, which can split a long path."""
    return re.sub(r"\S*/(\S*)", r"\1", text)
