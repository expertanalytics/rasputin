"""Legends: the colour scheme of a source's classes, as TOML data files.

`docs/increments/108-legend-split.md`, *4. The design*.

A legend file holds ``name``, ``preset`` (the ParaView preset's name), an
optional ``[no_class]`` (the label and colour of code 0) and ``[[class]]``
tables of ``code``, ``label`` and ``colour``; every key is required and no
other is allowed. It says nothing of what a source's values mean: that is
the land-cover class file's (`tin_engine.land_cover_classes`), which names
its legend. The shipped legends live in ``tin_engine/data/legends/``.
"""

from __future__ import annotations

import os
from collections.abc import Mapping
from importlib import resources
from pathlib import Path
from typing import Annotated, Any

from pydantic import BaseModel, ConfigDict, Field, model_validator

from tin_engine.toml_data import parse

Code = Annotated[int, Field(ge=1, lt=2**31)]
Colour = Annotated[str, Field(pattern=r"^#[0-9a-f]{6}$")]


class Legend(BaseModel):
    """A legend's name, its ParaView preset's name, and each code's label and
    colour: code 0 first when the file gives `[no_class]`, then file order."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    preset: str
    colours: Mapping[int, tuple[str, str]]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class _NoClass(_Strict):
    label: str
    colour: Colour


class _Class(_Strict):
    code: Code
    label: str
    colour: Colour


class _File(_Strict):
    name: str
    preset: str
    no_class: _NoClass | None = None
    classes: list[_Class] = Field(alias="class")

    @model_validator(mode="after")
    def _once(self) -> _File:
        codes = [c.code for c in self.classes]
        for code in codes:
            if codes.count(code) > 1:
                raise ValueError(f"class code {code} is listed twice")
        return self

    def legend(self) -> Legend:
        colours = {} if self.no_class is None else {0: (self.no_class.label, self.no_class.colour)}
        colours |= {c.code: (c.label, c.colour) for c in self.classes}
        return Legend(name=self.name, preset=self.preset, colours=colours)


def _shipped() -> Any:
    return resources.files("tin_engine") / "data" / "legends"


def shipped_names() -> list[str]:
    """The legends in the package, by name."""
    return sorted(p.name[:-5] for p in _shipped().iterdir() if p.name.endswith(".toml"))


def is_path(value: str) -> bool:
    """Whether a name-or-path value is a path: it ends in ``.toml`` or holds a
    path separator."""
    return value.endswith(".toml") or "/" in value or os.sep in value


def read(value: str, shipped: Any, shipped_names: list[str], kind: str) -> tuple[str, str]:
    """The text of a shipped ``kind`` file by name, or of a file by path, and
    the name its refusals give it."""
    if is_path(value):
        try:
            return Path(value).read_text(encoding="utf-8"), value
        except (OSError, UnicodeError) as exc:
            raise ValueError(f"{kind} file {value} cannot be read: {exc}") from exc
    if value not in shipped_names:
        raise ValueError(f"unknown {kind} {value}; use {', '.join(shipped_names)}, or a .toml file")
    return (shipped / f"{value}.toml").read_text(encoding="utf-8"), f"{value}.toml"


def load_legend(value: str) -> Legend:
    """A legend by shipped name or by path. Every refusal is a ``ValueError``
    naming the file."""
    text, shown = read(value, _shipped(), shipped_names(), "legend")
    return parse(text, _File, f"legend file {shown}").legend()
