"""Meshing rules as TOML data files, and the clean-up distances they imply.

`docs/increments/41-class-maps-as-data.md`, *4. The class-map file*.

A meshing-rules file says which of a source's features become which
constraint edges, at what accuracy, and, for a coded source, which cell
category each class counts as (``[cells]``). It is read with `tomllib`,
checked by a strict Pydantic model, and turned into the frozen
:class:`~tin_engine.feature_input.MeshingRules` the feature reader takes. A
coded file names its legend (`tin_engine.legends`), the colour scheme of its
codes: a shipped legend name or a path relative to the rules file. Class code 0
(an area no polygon covers) is accepted only when the file says
``accept_no_class = true``, and its legend must then colour it. The
shipped rules live in ``tin_engine/data/meshing_rules/``; ``property``, the
identity rules over the vocabulary's names, is built here, not shipped.
"""

from __future__ import annotations

from collections.abc import Iterable
from importlib import resources
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from tin_engine.feature_input import Category, MeshingRules, Split, split_key
from tin_engine.features import DEFAULT_VOCABULARY, EdgeProperty, EdgeVocabulary
from tin_engine.legends import Code, Legend, is_path, load_legend, read
from tin_engine.toml_data import parse

Edges = list[str] | Literal["drop"]
#: Bits 10 to 31 hold the rules' extra edge names; bit 9 is the cut's `seam`.
FIRST_EXTRA_BIT, MAX_EXTRAS = 10, 22
_CODED = {"code", "water"}
DEFAULT_NAMES = tuple(p.name for p in DEFAULT_VOCABULARY.properties)
#: Keys a rules file refuses, and the place each belongs instead.
_LEGEND, _CELLS = "the legend, not in the meshing rules", "the [cells] table"
_MOVED = {"label": _LEGEND, "colour": _LEGEND, "preset": _LEGEND, "no_class": _LEGEND,
          "category": f"{_CELLS}, under categories", "prefer": _CELLS}  # fmt: skip
NO_EDGES = (
    'unlisted = "land" needs an [edges] table: '
    "its land list names the edges an unlisted feature gets"
)


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class _Class(_Strict):
    code: Code | None = None
    water: bool = False
    value: str | None = None
    edges: Edges | None = None


class _EdgesTable(_Strict):
    water: Edges
    land: Edges


class _SubClass(_Strict):
    value: int | str
    code: Code


class _Split(_Strict):
    code: Code
    attribute: str
    classes: list[_SubClass] = Field(alias="class")


class _Vocabulary(_Strict):
    extra: list[str] = []


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
    code_system: str = ""
    positional_accuracy_m: Annotated[float, Field(gt=0, allow_inf_nan=False)] | None = None
    notice: str = ""
    licence: str = ""
    legend: str = ""
    unlisted: Literal["land", "drop", "refuse"]
    edges: _EdgesTable | None = None
    classes: list[_Class] = Field(alias="class")
    codes: dict[str, Code] = {}
    vocabulary: _Vocabulary = _Vocabulary()
    split: list[_Split] = []
    cells: _Cells | None = None
    accept_no_class: bool = False

    @model_validator(mode="before")
    @classmethod
    def _nothing_moved(cls, data: Any) -> Any:
        if isinstance(data, dict):
            _moved("", data)
            for i, c in enumerate(data.get("class", [])):
                _moved(f"class {i + 1}: ", c)
            for s in data.get("split", []):
                for c in s.get("class", []) if isinstance(s, dict) else ():
                    _moved(f"split of class {s.get('code')}: ", c)
        return data

    @model_validator(mode="after")
    def _consistent(self) -> _File:
        coded, extras = bool(self.code_system), self.vocabulary.extra
        if coded and self.positional_accuracy_m is None:
            raise ValueError("code_system needs positional_accuracy_m")
        if coded and self.edges is None:
            raise ValueError("code_system needs [edges]")
        if coded != bool(self.legend):
            raise ValueError(
                "a legend needs code_system" if self.legend else "code_system needs a legend"
            )
        if self.accept_no_class and not coded:
            raise ValueError("accept_no_class needs code_system")
        if self.cells is not None and not coded:
            raise ValueError("[cells] needs code_system")
        if self.unlisted == "land" and self.edges is None:
            raise ValueError(NO_EDGES)
        for i, c in enumerate(self.classes):
            given = c.model_fields_set
            wrong = given & ({"value"} if coded else _CODED)
            needed = {"code"} if coded else {"value", "edges"}
            if wrong or needed - given:
                bad = ", ".join(sorted(wrong)) or ", ".join(sorted(needed - given))
                kind = "coded" if coded else "uncoded"
                raise ValueError(f"class {i + 1}: {bad} wrong or missing in {kind} rules")
        listed = [c.code for c in self.classes if c.code is not None]
        if self.split and not coded:
            raise ValueError("[[split]] needs code_system")
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
        repeated = [n for n in extras if n in DEFAULT_NAMES]
        if repeated:
            raise ValueError(f"vocabulary extra {', '.join(repeated)} is a default name")
        if len(extras) > MAX_EXTRAS:
            raise ValueError(f"{len(extras)} vocabulary extras; bits 10 to 31 hold {MAX_EXTRAS}")
        known = {*DEFAULT_NAMES, *extras}
        lists = [c.edges for c in self.classes]
        lists += [self.edges.water, self.edges.land] if self.edges else []
        for names in lists:
            for name in names if isinstance(names, list) else ():
                if name not in known:
                    raise ValueError(f"edges name {name!r} is not in the vocabulary or its extras")
        return self

    def rules(self, legend: Legend | None) -> MeshingRules:
        def edges(names: Edges) -> tuple[str, ...] | None:
            return None if names == "drop" else tuple(names)

        classes: dict[str, tuple[str, ...] | None] = {}
        for c in self.classes:
            chosen = c.edges
            if chosen is None and self.edges is not None:
                chosen = self.edges.water if c.water else self.edges.land
            classes[str(c.code) if c.value is None else c.value] = edges(chosen or [])
        for text, code in self.codes.items():
            classes[text] = classes[str(code)]
        splits = sorted(self.split, key=lambda s: s.code)
        cells = self.cells or _Cells()
        given = {code: cat for cat, codes in cells.categories.items() for code in codes}
        categories: dict[int, Category] = {}
        for s in self.split:
            if s.code in given:
                categories |= {sub.code: given[s.code] for sub in s.classes}
        land = edges(self.edges.land) if self.edges else None
        otherwise = {
            "land": land if land is not None else "drop",
            "drop": "drop",
            "refuse": "refuse",
        }[self.unlisted]
        return MeshingRules(
            name=self.name,
            attribute=self.attribute,
            classes=classes,
            otherwise=otherwise,
            notice=self.notice,
            codes=self.code_system,
            legend=legend,
            licence=self.licence,
            accuracy_m=self.positional_accuracy_m,
            string_codes=self.codes,
            extra_names=tuple(self.vocabulary.extra),
            splits=tuple(
                Split(
                    code=s.code,
                    attribute=s.attribute,
                    codes={str(split_key(c.value)): c.code for c in s.classes},
                )
                for s in splits
            ),
            categories=categories | given,
            water_codes=frozenset(c.code for c in self.classes if c.code and c.water),
            prefer=cells.prefer,
            accept_no_class=self.accept_no_class,
        )


class CleanUpDistances(_Strict):
    """The land-cover clean-up's three distances, in `FeatureRequest`'s names."""

    repair_m: float
    tolerance_m: float
    outline_snap_m: float


def derived_cleanup(accuracy_m: float, cell: float) -> CleanUpDistances:
    """Section 5.1: s = max(accuracy, cell); repair s/100, simplify s/2,
    snap min(s/20, 50), half the read region's 100 m margin."""
    s = max(accuracy_m, cell)
    return CleanUpDistances(repair_m=s / 100, tolerance_m=s / 2, outline_snap_m=min(s / 20, 50.0))


def _shipped() -> Any:
    return resources.files("tin_engine") / "data" / "meshing_rules"


def shipped_names() -> list[str]:
    """The meshing rules in the package, by name."""
    return sorted(p.name[:-5] for p in _shipped().iterdir() if p.name.endswith(".toml"))


def load(value: str) -> MeshingRules:
    """A ``--features-rules`` value: a shipped name, ``property``, or a path (it
    ends in ``.toml`` or holds a path separator). The legend a coded file names
    is loaded with it, and every class and sub-class code must be in it. Every
    refusal is a ``ValueError`` naming the file."""
    if value == "property":
        classes = {n: (n,) for n in DEFAULT_NAMES}
        return MeshingRules(name="property", attribute="property", classes=classes, cleaned=False)
    text, shown = read(value, _shipped(), [*shipped_names(), "property"], "meshing rules")
    label = f"meshing rules file {shown}"
    found = parse(text, _File, label)
    if not found.legend:
        return found.rules(None)
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
    return found.rules(legend)


def known_codes(rules: MeshingRules) -> frozenset[int]:
    """The codes a class polygon may carry: every listed class and sub-class
    code, and 0 when the rules set `accept_no_class`."""
    listed = {int(key) for key in rules.classes if key not in rules.string_codes}
    listed |= {code for s in rules.splits for code in s.codes.values()}
    return frozenset(listed | ({0} if rules.accept_no_class else set()))


def run_vocabulary(rules: Iterable[MeshingRules], more: tuple[str, ...] = ()) -> EdgeVocabulary:
    """Section 6: the default vocabulary plus every rules file's extras and
    then `more`, in order of first appearance, on bits 10 upward; the default
    itself with none."""
    extras = list(dict.fromkeys([*(n for m in rules for n in m.extra_names), *more]))
    if not extras:
        return DEFAULT_VOCABULARY
    if len(extras) > MAX_EXTRAS:
        raise ValueError(f"{len(extras)} extra edge names in one run; bits 10 to 31 hold 22")
    added = (EdgeProperty(name=n, bit=FIRST_EXTRA_BIT + i) for i, n in enumerate(extras))
    return EdgeVocabulary(properties=(*DEFAULT_VOCABULARY.properties, *added))
