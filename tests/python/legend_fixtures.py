"""Legend files written by hand for increment 41's suites.

`docs/increments/41-class-maps-as-data.md`, section 4 (the file) and section 7
(the red suite). `coded_map` is the smallest coded map section 4 allows plus
the keys its example shows (`notice`, `preset`, `unlisted`); `uncoded_map` the
shape of today's `property` map. Every refusal test changes one line of one
of these, and `test_legends.TestTheBaseFilesLoad` shows both load as
written, so a refusal is the change's and nothing else's.
"""

from __future__ import annotations

import re
from pathlib import Path

#: CORINE's code system, verbatim, so a user map can meet `corine` in one mesh.
CORINE_SYSTEM = "CORINE Land Cover level-3 code"


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
) -> str:
    """A coded map: Forest (code 1) and Lake (code 2, water), or ``classes``."""
    classes = classes or (
        '[[class]]\ncode = 1\nlabel = "Forest"\ncolour = "#1f8d49"\n\n'
        '[[class]]\ncode = 2\nlabel = "Lake"\ncolour = "#2532e4"\nwater = true\n'
    )
    head = (
        f'name = "{name}"\n'
        f'attribute = "{attribute}"\n'
        f'code_system = "{system}"\n'
        f"positional_accuracy_m = {accuracy}\n"
        'notice = "Test data, made up"\n'
        f'preset = "rasputin {name}"\n'
        f'unlisted = "{unlisted}"\n'
    )
    if licence:
        head += f'licence = "{licence}"\n'
    body = (
        '\n[edges]\nwater = ["land_cover", "water"]\nland = ["land_cover"]\n\n'
        '[no_class]\nlabel = "No polygon; constraint lines"\ncolour = "#3c3c3c"\n\n'
    )
    tail = f"\n[codes]\n{codes}" if codes else ""
    tail += f"\n[vocabulary]\nextra = [{extra}]\n" if extra else ""
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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def basenames(text: str) -> str:
    """``text`` with every path cut to its last part. A test's ``tmp_path``
    carries the test's name ("..._named_river_i0"), so a word looked for in a
    message is looked for outside the directories, where only the producer
    can have put it. Works on `plain` output, which can split a long path."""
    return re.sub(r"\S*/(\S*)", r"\1", text)
