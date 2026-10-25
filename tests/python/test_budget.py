"""The triangle budget's search, pure: tests A6 to A9 and A11 (its `search` half).

`docs/increments/104-triangle-budget.md`, *3. The rule*. `budget.search` is
given a `trial(T)` callable and a `Window`; here every `trial` is a synthetic
count function of the tolerance, so no kernel runs. The count functions are
the design's `pure_search.py` ones: c(T) = F for T >= E_F, else
F + round(A ((E_F / T)^k - 1)), with largest error min(T, E_F). The trial
counts pinned below are the ones that rule measured (`pure-search.log`).
"""

from __future__ import annotations

import math
from collections.abc import Callable

import pytest
from tin_engine.budget import (
    FLOOR_TOLERANCE,
    MAX_TRIALS,
    Found,
    Trial,
    UnderFloor,
    Window,
    search,
)

CountFunction = Callable[[float], Trial]


def synthetic(
    floor: int,
    floor_error: float,
    k: float,
    a: float,
    *,
    rise_at: float | None = None,
    gap: tuple[int, int] | None = None,
) -> CountFunction:
    """The design's synthetic count function. `rise_at`: every count below
    that tolerance raised by 3 % (A8). `gap`: every count in the closed range
    mapped to one above it, so no tolerance gives a count inside (A9)."""

    def trial(tolerance: float) -> Trial:
        if tolerance >= floor_error:
            count = floor
        else:
            count = floor + round(a * ((floor_error / tolerance) ** k - 1))
            if rise_at is not None and tolerance < rise_at:
                count = round(count * 1.03)
            if gap is not None and gap[0] <= count <= gap[1]:
                count = gap[1] + 1
        return Trial(tolerance, count, min(tolerance, floor_error))

    return trial


class Recorded:
    """A `trial` that remembers every tolerance it was asked for."""

    def __init__(self, inner: CountFunction) -> None:
        self.inner = inner
        self.asked: list[float] = []

    def __call__(self, tolerance: float) -> Trial:
        self.asked.append(tolerance)
        return self.inner(tolerance)


def default_cap(asked: int, share: float) -> int:
    """The cap the grid uses: ⌊share · N⌋, at least N + 1."""
    return max(math.floor(share * asked), asked + 1)


def landed(found: Found, window: Window) -> None:
    """In the window, not missed, the chosen trial one of those run."""
    assert not found.missed, found
    assert window.holds(found.chosen.triangles), found
    assert found.chosen in found.trials


# ---------------------------------------------------------------- the window


class TestWindow:
    def test_low_is_as_far_under_the_ask_as_the_cap_is_over(self) -> None:
        window = Window(asked=7386, cap=7755)
        assert window.low == 7017

    @pytest.mark.parametrize(
        ("count", "inside"),
        [(7016, False), (7017, True), (7386, True), (7755, True), (7756, False)],
    )
    def test_holds_both_ends_inclusive(self, count: int, inside: bool) -> None:
        assert Window(asked=7386, cap=7755).holds(count) is inside

    @pytest.mark.parametrize(("asked", "cap"), [(100, 100), (100, 99), (0, 5), (-3, 5)])
    def test_refuses_a_cap_not_over_the_ask_and_an_ask_under_one(
        self, asked: int, cap: int
    ) -> None:
        with pytest.raises(ValueError):
            Window(asked=asked, cap=cap)

    def test_the_constants(self) -> None:
        assert FLOOR_TOLERANCE == 1e6
        assert MAX_TRIALS == 12


# ---------------------------------------------------------------- A6: the synthetic grid


def grid() -> list[tuple[CountFunction, float, Window]]:
    """Section 7's 360 searches: (count function, E_F, window)."""
    cases = []
    for floor in (0, 50, 28561):
        floor_error = 216.0 if floor else 50.0
        for k in (0.6, 1.1, 2.0):
            a = 3e5 / (1000**k - 1)
            trial = synthetic(floor, floor_error, k, a)
            lo = max(trial(floor_error / 8).triangles, int(1.3 * floor) + 1, 20)
            hi = trial(floor_error / 2000).triangles
            for i in range(20):
                asked = round(lo * (hi / lo) ** (i / 19))
                for share in (1.05, 1.02):
                    cases.append((trial, floor_error, Window(asked, default_cap(asked, share))))
    return cases


GRID = grid()


class TestA6Grid:
    def test_the_grid_has_360_searches(self) -> None:
        assert len(GRID) == 360

    @pytest.mark.parametrize(("trial", "floor_error", "window"), GRID)
    def test_lands_in_the_window_in_at_most_seven_trials(
        self, trial: CountFunction, floor_error: float, window: Window
    ) -> None:
        recorded = Recorded(trial)
        found = search(recorded, window)
        landed(found, window)
        assert len(found.trials) <= 7, found.trials
        assert recorded.asked[0] == FLOOR_TOLERANCE
        assert recorded.asked[1] == floor_error / 4
        assert [t.tolerance for t in found.trials] == recorded.asked
        assert found.chosen.triangles <= window.cap

    def test_the_same_trials_on_a_second_call(self) -> None:
        for trial, _, window in GRID[::17]:
            assert search(trial, window) == search(trial, window)

    def test_the_first_trial_is_the_floor_tolerance_given(self) -> None:
        recorded = Recorded(synthetic(50, 216.0, 1.1, 3e5 / (1000**1.1 - 1)))
        search(recorded, Window(5000, 5250), floor_tolerance=1e4)
        assert recorded.asked[0] == 1e4


# ---------------------------------------------------------------- A7: the floor


VELHAS_LIKE = synthetic(49, 76.598, 1.5, 20.0)


class TestA7Floor:
    def test_a_window_under_the_floor_raises_with_the_floor_trial(self) -> None:
        with pytest.raises(UnderFloor) as caught:
            search(VELHAS_LIKE, Window(40, 42))
        assert caught.value.floor == Trial(FLOOR_TOLERANCE, 49, 76.598)

    def test_a_window_holding_the_floor_returns_it_at_its_largest_error(self) -> None:
        recorded = Recorded(VELHAS_LIKE)
        found = search(recorded, Window(48, 50))
        assert recorded.asked == [FLOOR_TOLERANCE]
        assert len(found.trials) == 1
        assert not found.missed
        assert found.chosen.triangles == 49
        assert found.chosen.tolerance == 76.598
        assert found.chosen.max_error == 76.598


# ---------------------------------------------------------------- A8: a rise


class TestA8Rise:
    @pytest.mark.parametrize("asked", [40_000, 65_893, 150_000, 300_000, 600_000])
    def test_a_three_percent_rise_below_three_metres_still_lands_in_five(self, asked: int) -> None:
        window = Window(asked, math.floor(1.05 * asked))
        found = search(synthetic(28561, 216.0, 2.0, 80.3, rise_at=3.0), window)
        landed(found, window)
        assert len(found.trials) == 5, found.trials
        assert found.chosen.triangles <= window.cap


# ---------------------------------------------------------------- A9: the trial cap

GAP = synthetic(28561, 216.0, 2.0, 80.3, gap=(60_000, 72_000))
GAP_WINDOW = Window(65_893, 69_187)


class TestA9TrialCap:
    def test_a_gap_wider_than_the_window_ends_missed_after_twelve(self) -> None:
        recorded = Recorded(GAP)
        found = search(recorded, GAP_WINDOW)
        assert found.missed
        assert len(found.trials) == 12
        assert len(recorded.asked) == 12
        assert found.chosen.triangles == 59_949

    def test_the_largest_count_under_the_cap_is_chosen(self) -> None:
        found = search(GAP, GAP_WINDOW)
        under = [t.triangles for t in found.trials if t.triangles <= GAP_WINDOW.cap]
        assert found.chosen.triangles == max(under)
        assert any(t.triangles > GAP_WINDOW.cap for t in found.trials)

    @pytest.mark.parametrize("cap", [3, 4, 6, 9])
    def test_every_trial_counts_toward_the_cap_given(self, cap: int) -> None:
        recorded = Recorded(GAP)
        found = search(recorded, GAP_WINDOW, max_trials=cap)
        assert len(recorded.asked) == cap
        assert len(found.trials) == cap
        assert found.missed
        assert found.chosen.triangles <= GAP_WINDOW.cap

    def test_a_count_stuck_just_above_the_floor_ends_at_the_cap(self) -> None:
        """The coarse-to-fine quartering never reaches 1.5 F and the error
        never reaches 0: only the trial cap ends it."""

        def stuck(tolerance: float) -> Trial:
            count = 1000 if tolerance >= 50.0 else 1001
            return Trial(tolerance, count, min(tolerance, 50.0))

        recorded = Recorded(stuck)
        found = search(recorded, Window(1500, 1575))
        assert len(recorded.asked) == MAX_TRIALS
        assert found.missed
        assert found.chosen.triangles == 1001

    def test_a_plateau_under_one_and_a_half_floors_ends_saturated(self) -> None:
        """Every node in at 1.4 F: the quartering loop stops at the
        saturated trial instead of running on."""

        def plateau(tolerance: float) -> Trial:
            if tolerance >= 50.0:
                return Trial(tolerance, 1000, 50.0)
            count = min(1000 + round(10 * ((50.0 / tolerance) ** 0.5 - 1)), 1400)
            return Trial(tolerance, count, 0.0 if count == 1400 else min(tolerance, 50.0))

        found = search(plateau, Window(1500, 1575))
        assert found.missed
        assert found.chosen.triangles == 1400
        assert found.trials[-1].max_error == 0.0
        assert len(found.trials) <= MAX_TRIALS


# ---------------------------------------------------------------- A11: saturation


def everything(tolerance: float) -> Trial:
    """Every DEM node a vertex at every tolerance: 640 triangles, error 0."""
    return Trial(tolerance, 640, 0.0)


def saturating(tolerance: float) -> Trial:
    """A6's form (F 12, E_F 36.938, k 0.9, A 60), counts capped at 640
    with error 0 there."""
    t = synthetic(12, 36.938, 0.9, 60.0)(tolerance)
    count = min(t.triangles, 640)
    return Trial(tolerance, count, 0.0 if count == 640 else t.max_error)


class TestA11Saturation:
    def test_a_saturated_floor_in_the_window_is_the_answer_at_tolerance_zero(self) -> None:
        recorded = Recorded(everything)
        found = search(recorded, Window(639, 670))
        assert len(recorded.asked) == 1
        assert not found.missed
        assert found.chosen.triangles == 640
        assert found.chosen.tolerance == 0.0

    def test_a_saturated_floor_under_the_window_ends_missed_after_one(self) -> None:
        recorded = Recorded(everything)
        found = search(recorded, Window(1000, 1050))
        assert len(recorded.asked) == 1
        assert found.missed
        assert found.chosen.triangles == 640

    def test_saturation_reached_mid_search_ends_missed_after_three(self) -> None:
        recorded = Recorded(saturating)
        found = search(recorded, Window(700, 735))
        assert len(recorded.asked) == 3
        assert found.missed
        assert found.chosen.triangles == 640
        assert found.chosen.max_error == 0.0
