"""The triangle budget: the tolerance whose mesh has about the triangles asked for.

`docs/increments/104-triangle-budget.md`, *3. The rule*. Pure: the search
reads only what a `trial(T)` callable returns, the written triangle count and
the largest error of the mesh at tolerance T, and never holds a mesh. It runs
the floor (the mesh no tolerance can shrink), then tolerances coarse to fine
from the floor's own largest error, then a secant on log(count - floor)
against log T, kept inside the bracket of tolerances already known too small
and too large. It stops at the first count in the window, at a mesh with
every node on it, or after `max_trials` trials.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

#: Above any height difference on Earth and any NoData sentinel's magnitude,
#: so the floor trial inserts nothing.
FLOOR_TOLERANCE = 1e6
#: Every trial counts; nearly twice the most a search was measured to take.
MAX_TRIALS = 12


@dataclass(frozen=True, slots=True)
class Window:
    """About `asked` triangles and never more than `cap`: every count from
    `low` = 2 asked - cap to `cap`, as far under the ask as the cap is over."""

    asked: int
    cap: int

    def __post_init__(self) -> None:
        if self.asked < 1:
            raise ValueError(f"must be at least 1, got {self.asked}")
        if self.cap <= self.asked:
            raise ValueError(f"the cap {self.cap} must be over the {self.asked} asked for")

    @property
    def low(self) -> int:
        return 2 * self.asked - self.cap

    def holds(self, count: int) -> bool:
        return self.low <= count <= self.cap


@dataclass(frozen=True, slots=True)
class Trial:
    """One mesh: the tolerance it was refined to, its triangles, its largest error."""

    tolerance: float
    triangles: int
    max_error: float


@dataclass(frozen=True, slots=True)
class Found:
    """The trial chosen, every trial in the order run, and whether the
    chosen count is outside the window (the most under the cap instead)."""

    chosen: Trial
    trials: tuple[Trial, ...]
    missed: bool


class UnderFloor(ValueError):  # noqa: N818 -- the name the suite pins
    """The cap is under the floor's count: no tolerance gives so few triangles."""

    def __init__(self, floor: Trial) -> None:
        super().__init__(f"the floor has {floor.triangles} triangles")
        self.floor = floor


def search(
    trial: Callable[[float], Trial],
    window: Window,
    *,
    floor_tolerance: float = FLOOR_TOLERANCE,
    max_trials: int = MAX_TRIALS,
) -> Found:
    """The tolerance whose `trial` count is in `window`. A floor answer is
    reported at the floor's largest error, the tolerance at and above which
    the mesh is the floor. Raises `UnderFloor` when the floor is over the cap.
    A trial with largest error 0 has every node on the mesh and ends the
    search; so does the `max_trials`-th trial. Without a count in the window
    the chosen trial is the one with the most triangles not over the cap."""
    trials: list[Trial] = []
    too_small: float | None = None  # the largest T whose count was over the cap
    too_large: float | None = None  # the smallest T whose count was under the window

    def run(tolerance: float) -> Trial:
        nonlocal too_small, too_large
        t = trial(tolerance)
        trials.append(t)
        if t.triangles > window.cap:
            too_small = tolerance if too_small is None else max(too_small, tolerance)
        elif t.triangles < window.low:
            too_large = tolerance if too_large is None else min(too_large, tolerance)
        return t

    def finish(missed: bool) -> Found:
        if not missed:
            return Found(trials[-1], tuple(trials), False)
        under = [t for t in trials if t.triangles <= window.cap]
        return Found(max(under, key=lambda t: t.triangles), tuple(trials), True)

    def going(t: Trial) -> bool:
        return not window.holds(t.triangles) and t.max_error > 0 and len(trials) < max_trials

    floor = run(floor_tolerance)
    f, e_f = floor.triangles, floor.max_error
    if f > window.cap:
        raise UnderFloor(floor)
    if window.holds(f) or e_f == 0:
        return Found(Trial(e_f, f, e_f), tuple(trials), not window.holds(f))
    t = run(e_f / 4)
    while t.triangles < 1.5 * f and going(t):
        t = run(t.tolerance / 4)
    if going(t):  # a second trial well above the floor
        t = run(t.tolerance / 4)
    while going(t):
        a, b = trials[-2], trials[-1]
        ca, cb = max(a.triangles - f, 1), max(b.triangles - f, 1)
        guess = math.nan
        if a.tolerance != b.tolerance and ca != cb:
            k = math.log(ca / cb) / math.log(a.tolerance / b.tolerance)
            if k < 0:
                guess = b.tolerance * (max(window.asked - f, 1) / cb) ** (1 / k)
        inside = (
            math.isfinite(guess)
            and guess > 0
            and (too_small is None or guess > too_small)
            and (too_large is None or guess < too_large)
        )
        if inside:
            t = run(guess)
        elif too_small is None:
            t = run((too_large if too_large is not None else b.tolerance) / 4)
        elif too_large is None:
            t = run(too_small * 4)
        else:
            t = run(math.sqrt(too_small * too_large))
    return finish(not window.holds(t.triangles))
