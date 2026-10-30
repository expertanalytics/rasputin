"""Land-cover class files as TOML data, and the clean-up distances they imply.

`docs/increments/41-class-maps-as-data.md`, *4. The class-map file*.

A land-cover class file says what a class polygon's attribute means: the
code system its values are integer codes of, the water classes, the file's
positional accuracy, and which cell category each class counts as
(``[cells]``). It is read with `tomllib`, checked by a strict Pydantic model,
and turned into the frozen
:class:`~tin_engine.feature_input.LandCoverClasses` the readers take. A file
names its legend (`tin_engine.legends`), the colour scheme of its codes: a
shipped legend name or a path relative to the class file. Class code 0 (an
area no polygon covers) is accepted only when the file says
``accept_no_class = true``, and its legend must then colour it. The shipped
files live in ``tin_engine/data/land_cover_classes/``.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from tin_engine.feature_input import Category, LandCoverClasses, Split, split_key
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.legends import Code, Legend, is_path, load_legend, read
from tin_engine.toml_data import parse

#: Keys a class file refuses, and the place each belongs instead.
_LEGEND, _CELLS = "the legend, not in the land-cover class file", "the [cells] table"
_MOVED = {"label": _LEGEND, "colour": _LEGEND, "preset": _LEGEND, "no_class": _LEGEND,
          "category": f"{_CELLS}, under categories", "prefer": _CELLS}  # fmt: skip
#: Names an ``--edges`` selector means already, so no category may take them.
SELECTORS = frozenset({"all", "water", *(p.name for p in DEFAULT_VOCABULARY.properties)})


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class _Class(_Strict):
    code: Code
    water: bool = False


class _SubClass(_Strict):
    value: int | str
    code: Code


class _Split(_Strict):
    code: Code
    attribute: str
    classes: list[_SubClass] = Field(alias="class")


class _Cells(_Strict):
    categories: dict[Category, list[Code]] = {}
    prefer: dict[Category, Code] = {}


def _moved(where: str, table: Any) -> None:
    """Refuse a key that belongs in the legend or in `[cells]`, naming its place."""
    for key in _MOVED.keys() & (table.keys() if isinstance(table, dict) else set()):
        raise ValueError(f"{where}{key} belongs in {_MOVED[key]}")


class _File(_Strict):
    name: str
    attribute: str
    code_system: str
    positional_accuracy_m: Annotated[float, Field(gt=0, allow_inf_nan=False)]
    notice: str = ""
    licence: str = ""
    legend: str
    unlisted: Literal["keep", "drop", "refuse"]
    classes: list[_Class] = Field(alias="class")
    codes: dict[str, Code] = {}
    split: list[_Split] = []
    cells: _Cells | None = None
    accept_no_class: bool = False

    @model_validator(mode="before")
    @classmethod
    def _nothing_moved(cls, data: Any) -> Any:
        if isinstance(data, dict):
            _moved("", data)
            cells = data.get("cells")
            named = cells.get("categories") if isinstance(cells, dict) else None
            for key in SELECTORS & (named.keys() if isinstance(named, dict) else set()):
                raise ValueError(f"[cells] categories: {key} is an --edges selector, no category")
            for i, c in enumerate(data.get("class", [])):
                _moved(f"class {i + 1}: ", c)
            for s in data.get("split", []):
                for c in s.get("class", []) if isinstance(s, dict) else ():
                    _moved(f"split of class {s.get('code')}: ", c)
        return data

    @model_validator(mode="after")
    def _consistent(self) -> _File:
        listed = [c.code for c in self.classes]
        split = [s.code for s in self.split]
        for s in self.split:
            if s.code not in listed or split.count(s.code) > 1:
                raise ValueError(f"split of class {s.code}: not a listed class, or split twice")
            if s.attribute.lower() == self.attribute.lower():
                raise ValueError(f"split of class {s.code}: {s.attribute} is the first attribute")
            keys = [split_key(c.value) for c in s.classes]
            for key, sub in zip(keys, s.classes, strict=True):
                if key is None or keys.count(key) > 1:
                    raise ValueError(
                        f"split of class {s.code}: value {sub.value!r} has no key or recurs"
                    )
        every = [*listed, *(sub.code for s in self.split for sub in s.classes)]
        for code in every:
            if every.count(code) > 1:
                raise ValueError(f"class code {code} is listed twice")
        cells = self.cells or _Cells()
        for text, code in [*self.codes.items(), *cells.prefer.items()]:
            if code not in listed:
                raise ValueError(f"codes or prefer: {text} = {code} names no listed class")
        counted = [code for codes in cells.categories.values() for code in codes]
        for code in counted:
            if code not in every or counted.count(code) > 1:
                raise ValueError(
                    f"[cells] categories: class {code} is not listed, or counted twice"
                )
        return self

    def classes_of(self, legend: Legend) -> LandCoverClasses:
        cells = self.cells or _Cells()
        given = {code: cat for cat, codes in cells.categories.items() for code in codes}
        categories: dict[int, Category] = {}
        for s in self.split:
            if s.code in given:
                categories |= {sub.code: given[s.code] for sub in s.classes}
        return LandCoverClasses(
            name=self.name,
            attribute=self.attribute,
            codes=self.code_system,
            legend=legend,
            classes=frozenset([*(str(c.code) for c in self.classes), *self.codes]),
            otherwise=self.unlisted,
            notice=self.notice,
            licence=self.licence,
            accuracy_m=self.positional_accuracy_m,
            string_codes=self.codes,
            splits=tuple(
                Split(
                    code=s.code,
                    attribute=s.attribute,
                    codes={str(split_key(c.value)): c.code for c in s.classes},
                )
                for s in sorted(self.split, key=lambda s: s.code)
            ),
            categories=categories | given,
            water_codes=frozenset(c.code for c in self.classes if c.water),
            prefer=cells.prefer,
            accept_no_class=self.accept_no_class,
        )


class CleanUpDistances(_Strict):
    """The land-cover clean-up's three distances, in `EdgeCleanup`'s names."""

    repair_m: float
    max_shift_m: float
    outline_snap_m: float


def derived_cleanup(accuracy_m: float, cell: float) -> CleanUpDistances:
    """Section 5.1: s = max(accuracy, cell); repair s/100, simplify s/2,
    snap min(s/20, 50), half the read region's 100 m margin."""
    s = max(accuracy_m, cell)
    return CleanUpDistances(repair_m=s / 100, max_shift_m=s / 2, outline_snap_m=min(s / 20, 50.0))


def _shipped() -> Any:
    return resources.files("tin_engine") / "data" / "land_cover_classes"


def shipped_names() -> list[str]:
    """The land-cover class files in the package, by name."""
    return sorted(p.name[:-5] for p in _shipped().iterdir() if p.name.endswith(".toml"))


def load(value: str) -> LandCoverClasses:
    """A ``--land-cover-classes`` value: a shipped name or a path (it ends in
    ``.toml`` or holds a path separator). The legend the file names is loaded
    with it, and every class and sub-class code must be in it. Every refusal
    is a ``ValueError`` naming the file."""
    text, shown = read(value, _shipped(), shipped_names(), "land-cover class file")
    label = f"land-cover class file {shown}"
    found = parse(text, _File, label)
    named = found.legend
    if is_path(named) and is_path(value):
        named = str(Path(value).parent / named)
    try:
        legend = load_legend(named)
    except ValueError as exc:
        raise ValueError(f"{label}: {exc}") from exc
    codes = (
        [0] * found.accept_no_class
        + [c.code for c in found.classes]
        + [s.code for t in found.split for s in t.classes]
    )
    missing = [str(code) for code in codes if code not in legend.colours]
    if missing:
        raise ValueError(f"{label}: class {', '.join(missing)} not in legend {legend.name}")
    return found.classes_of(legend)


def known_codes(classes: LandCoverClasses) -> frozenset[int]:
    """The codes a class polygon may carry: every listed class and sub-class
    code, and 0 when the file sets `accept_no_class`."""
    listed = {int(key) for key in classes.classes if key not in classes.string_codes}
    listed |= {code for s in classes.splits for code in s.codes.values()}
    return frozenset(listed | ({0} if classes.accept_no_class else set()))
