"""What one window's refine takes and gives, cut or not.

`cli.py` refines every window through these, an uncut run's one and each
piece of a cut run, so the uncut run needs them without the cut path's
modules.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from tin_engine._core import CheckPoints, RefineOutcome
from tin_engine.budget import Window
from tin_engine.elevation import Trimmed
from tin_engine.run_record import Value
from tin_engine.target_grid import Block, TargetGrid


class RunError(ValueError):
    """A run the engine refuses, in its own words; no file is written."""


class Once[T]:
    """A value handed on once: `take` returns it and forgets it,
    so the holder keeps nothing alive after; a second `take` is a bug."""

    def __init__(self, value: T) -> None:
        self._value: T | None = value

    def take(self) -> T:
        value, self._value = self._value, None
        if value is None:  # explicit: `python -O` strips an assert
            raise RuntimeError("Once.take called twice")
        return value


@dataclass(frozen=True, slots=True)
class How:
    """How a start mesh is refined (with ``grid`` and ``checks``, the final
    check against a resampled grid; ``checks`` as blocks are filed into a
    store by the refine that takes them). With ``budget``, ``tolerance`` is
    set per trial by the search for the window's count."""

    dem: str
    tolerance: float
    min_angle: float
    feet: bool
    grid: TargetGrid | None
    checks: Iterable[Block] | CheckPoints | None
    budget: Window | None = None


@dataclass(frozen=True, slots=True)
class Refined:
    out: RefineOutcome
    final: RefineOutcome
    z: npt.NDArray[np.float64]
    trimmed: Trimmed
    values: Mapping[str, Value]  # the run record's entries
    on_frozen: int
