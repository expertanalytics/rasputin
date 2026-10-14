"""A TOML data file read into a strict Pydantic model, its refusals naming the file.

Every refusal is a `ValueError` reading `<label>: <key path>: <message>`; a
TOML syntax error keeps the parser's words, which give the line.
"""

from __future__ import annotations

import tomllib
from collections.abc import Sequence
from typing import Any

from pydantic import BaseModel, ValidationError


def parse[M: BaseModel](text: str, model: type[M], label: str) -> M:
    """`text` read as TOML and validated as `model`."""
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        raise ValueError(f"{label}: {exc}") from exc
    return validate(data, model, label)


def validate[M: BaseModel](data: Any, model: type[M], label: str, at: Sequence[str] = ()) -> M:
    """`data` validated as `model`; `at` is the key path of `data` in the file."""
    try:
        return model.model_validate(data)
    except ValidationError as exc:
        why = "; ".join(
            f"{'.'.join(map(str, [*at, *e['loc']])) or 'file'}: "
            f"{e['msg'].removeprefix('Value error, ')}"
            for e in exc.errors()
        )
        raise ValueError(f"{label}: {why}") from exc
