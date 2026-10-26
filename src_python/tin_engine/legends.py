"""Land-cover legends as TOML data files, and the clean-up distances they imply.

`docs/increments/41-class-maps-as-data.md`, *4. The class-map file*.

A legend file is read with `tomllib`, checked by a strict Pydantic model, and
turned into the frozen :class:`~tin_engine.feature_input.Legend` the feature
reader takes. The
shipped legends live in ``tin_engine/data/legends/``; ``property``, the
identity legend over the vocabulary's names, is built here, not shipped.
"""

from __future__ import annotations

import os
from collections.abc import Iterable
from importlib import resources
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from tin_engine.feature_input import Category, Legend, Split, split_key
from tin_engine.features import DEFAULT_VOCABULARY, EdgeProperty, EdgeVocabulary
from tin_engine.toml_data import parse

Code = Annotated[int, Field(ge=1, lt=2**31)]
Edges = list[str] | Literal["drop"]
#: Bits 10 to 31 hold a legend's extra edge names; bit 9 is the cut's `seam`.
FIRST_EXTRA_BIT, MAX_EXTRAS = 10, 22
_CODED = {"code", "label", "colour", "water", "category"}
DEFAULT_NAMES = tuple(p.name for p in DEFAULT_VOCABULARY.properties)
Colour = Annotated[str, Field(pattern=r"^#[0-9a-f]{6}$")]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class _Class(_Strict):
    code: Code | None = None
    label: str | None = None
    colour: Annotated[str, Field(pattern=r"^#[0-9a-f]{6}$")] | None = None
    water: bool = False
    category: Category | None = None
    value: str | None = None
    edges: Edges | None = None


class _EdgesTable(_Strict):
    water: Edges
    land: Edges


class _NoClass(_Strict):
    label: str
    colour: Annotated[str, Field(pattern=r"^#[0-9a-f]{6}$")]


class _SubClass(_Strict):
    value: int | str
    code: Code
    label: str
    colour: Colour


class _Split(_Strict):
    code: Code
    attribute: str
    classes: list[_SubClass] = Field(alias="class")


class _Vocabulary(_Strict):
    extra: list[str] = []


class _File(_Strict):
    name: str
    attribute: str
    code_system: str = ""
    positional_accuracy_m: Annotated[float, Field(gt=0, allow_inf_nan=False)] | None = None
    notice: str = ""
    licence: str = ""
    preset: str = ""
    unlisted: Literal["land", "drop", "refuse"]
    edges: _EdgesTable | None = None
    no_class: _NoClass | None = None
    classes: list[_Class] = Field(alias="class")
    codes: dict[str, Code] = {}
    vocabulary: _Vocabulary = _Vocabulary()
    split: list[_Split] = []

    @model_validator(mode="after")
    def _consistent(self) -> _File:
        coded, extras = bool(self.code_system), self.vocabulary.extra
        if coded and self.positional_accuracy_m is None:
            raise ValueError("code_system needs positional_accuracy_m")
        if coded and self.edges is None:
            raise ValueError("code_system needs [edges]")
        for i, c in enumerate(self.classes):
            given = c.model_fields_set
            wrong = given & ({"value"} if coded else _CODED)
            needed = {"code", "label", "colour"} if coded else {"value", "edges"}
            if wrong or needed - given:
                bad = ", ".join(sorted(wrong)) or ", ".join(sorted(needed - given))
                kind = "coded" if coded else "uncoded"
                raise ValueError(f"class {i + 1}: {bad} wrong or missing in a {kind} legend")
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
        for text, code in self.codes.items():
            if code not in listed:
                raise ValueError(f"codes: {text} = {code} names no listed class")
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

    def legend(self) -> Legend:
        def edges(names: Edges) -> tuple[str, ...] | None:
            return None if names == "drop" else tuple(names)

        classes: dict[str, tuple[str, ...] | None] = {}
        palette: dict[int, tuple[str, str]] = {}
        if self.no_class is not None:
            palette[0] = (self.no_class.label, self.no_class.colour)
        for c in self.classes:
            chosen = c.edges
            if chosen is None and self.edges is not None:
                chosen = self.edges.water if c.water else self.edges.land
            classes[str(c.code) if c.value is None else c.value] = edges(chosen or [])
            if c.code is not None and c.label is not None and c.colour is not None:
                palette[c.code] = (c.label, c.colour)
        for text, code in self.codes.items():
            classes[text] = classes[str(code)]
        for s in self.split:
            palette |= {sub.code: (sub.label, sub.colour) for sub in s.classes}
        splits = sorted(self.split, key=lambda s: s.code)
        categories = {c.code: c.category for c in self.classes if c.code and c.category}
        for s in self.split:
            if s.code in categories:
                categories |= {sub.code: categories[s.code] for sub in s.classes}
        land = edges(self.edges.land) if self.edges else None
        otherwise = {
            "land": land if land is not None else "drop",
            "drop": "drop",
            "refuse": "refuse",
        }[self.unlisted]
        return Legend(
            name=self.name,
            attribute=self.attribute,
            classes=classes,
            otherwise=otherwise,
            notice=self.notice,
            codes=self.code_system,
            palette=palette if self.code_system else None,
            preset=self.preset,
            licence=self.licence,
            accuracy_m=self.positional_accuracy_m or 0.0,
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
            categories=categories,
            water_codes=frozenset(c.code for c in self.classes if c.code and c.water),
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
    return resources.files("tin_engine") / "data" / "legends"


def shipped_names() -> list[str]:
    """The legends in the package, by name."""
    return sorted(p.name[:-5] for p in _shipped().iterdir() if p.name.endswith(".toml"))


def load(value: str) -> Legend:
    """A ``--features-legend`` value: a shipped name, ``property``, or a path (it
    ends in ``.toml`` or holds a path separator). Every refusal is a
    ``ValueError`` naming the file."""
    if value == "property":
        classes = {n: (n,) for n in DEFAULT_NAMES}
        return Legend(name="property", attribute="property", classes=classes)
    if value.endswith(".toml") or "/" in value or os.sep in value:
        path, shown = Path(value), value
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeError) as exc:
            raise ValueError(f"legend file {shown} cannot be read: {exc}") from exc
    else:
        shown = f"{value}.toml"
        if value not in shipped_names():
            known = ", ".join([*shipped_names(), "property"])
            raise ValueError(f"unknown legend {value}; use {known}, or a .toml file")
        text = (_shipped() / shown).read_text(encoding="utf-8")
    return parse(text, _File, f"legend file {shown}").legend()


def run_vocabulary(legends: Iterable[Legend], more: tuple[str, ...] = ()) -> EdgeVocabulary:
    """Section 6: the default vocabulary plus every legend's extras and then
    `more`, in order of first appearance, on bits 10 upward; the default itself
    with none."""
    extras = list(dict.fromkeys([*(n for m in legends for n in m.extra_names), *more]))
    if not extras:
        return DEFAULT_VOCABULARY
    if len(extras) > MAX_EXTRAS:
        raise ValueError(f"{len(extras)} extra edge names in one run; bits 10 to 31 hold 22")
    added = (EdgeProperty(name=n, bit=FIRST_EXTRA_BIT + i) for i, n in enumerate(extras))
    return EdgeVocabulary(properties=(*DEFAULT_VOCABULARY.properties, *added))
