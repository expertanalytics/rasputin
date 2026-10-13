"""The snow run's job: `SyntheticSpring`, `FrameSpec` and `SnowJob` in `snow/params.py`.

`docs/increments/46-hbv-snow.md`, *4.1 The blueprint* (flags -> `SnowJob`, a
frozen Pydantic model holding `HbvSnowParameters`, `SnowDistribution`,
`SunCorrection`, `SyntheticSpring`, `FrameSpec`, `Resources`), *4.5 The
synthetic case* (2027-04-01T00:00Z, 91 days, hourly steps, S = 4, frames
every 24 hours, 600 mm initial SWE) and *4.8 The Python surface* (times are
timezone-aware UTC; a naive one is refused). `snow.json` records the job as
the model dumped, so the dump must load back to an equal job.

Interface assumed (names chosen here and listed in the handback):

- `SyntheticSpring(start, days, step_minutes, samples_per_step)`; `days`,
  `step_minutes` and `samples_per_step` at least 1.
- `FrameSpec(every_hours=24)`, at least 1.
- `SnowJob(hbv, distribution, correction, spring, frames, resources,
  initial_swe_mm=600.0)`; `hbv`, `spring` and `frames` default to their
  models' defaults, `distribution` to `SnowDistribution.lognormal(0.4, 10)`;
  `correction` and `resources` are required. A frame interval that is not a
  whole number of steps is refused, located on `frames`.
"""

from __future__ import annotations

import datetime as dt
import importlib
from collections.abc import Callable
from types import ModuleType
from typing import Any

import pydantic
import pytest

from tin_engine.resources import Resources

RESOURCES = Resources(threads=2, memory_budget=1 << 30)


@pytest.fixture(scope="module")
def params() -> ModuleType:
    return importlib.import_module("tin_engine.snow.params")


def refused_on(loc: tuple[str, ...], make: Callable[[], object]) -> pydantic.ValidationError:
    with pytest.raises(pydantic.ValidationError) as caught:
        make()
    assert [e["loc"][: len(loc)] for e in caught.value.errors()] == [loc]
    return caught.value


class TestSyntheticSpring:
    def test_defaults_are_section_4_5s(self, params: ModuleType) -> None:
        s = params.SyntheticSpring()
        assert s.start == dt.datetime(2027, 4, 1, tzinfo=dt.UTC)
        assert (s.days, s.step_minutes, s.samples_per_step) == (91, 60, 4)

    def test_a_start_with_another_offset_is_the_same_instant(self, params: ModuleType) -> None:
        oslo = dt.timezone(dt.timedelta(hours=2))
        s = params.SyntheticSpring(start=dt.datetime(2027, 4, 1, 2, 0, tzinfo=oslo))
        assert s.start == dt.datetime(2027, 4, 1, tzinfo=dt.UTC)

    def test_a_naive_start_is_refused_saying_the_time_zone_is_missing(
        self, params: ModuleType
    ) -> None:
        error = refused_on(
            ("start",), lambda: params.SyntheticSpring(start=dt.datetime(2027, 4, 1))
        )
        assert "time zone" in str(error)

    @pytest.mark.parametrize("field", ["days", "step_minutes", "samples_per_step"])
    def test_counts_below_1_are_refused(self, params: ModuleType, field: str) -> None:
        refused_on((field,), lambda: params.SyntheticSpring(**{field: 0}))

    def test_it_is_frozen(self, params: ModuleType) -> None:
        with pytest.raises(pydantic.ValidationError):
            params.SyntheticSpring().days = 3


class TestSnowJob:
    def job(self, params: ModuleType, **kw: object) -> Any:
        return params.SnowJob(
            correction=params.SunCorrection(mode="melt_factor"), resources=RESOURCES, **kw
        )

    def test_defaults_are_section_4_5s(self, params: ModuleType) -> None:
        job = self.job(params)
        assert job.hbv == params.HbvSnowParameters()
        assert job.distribution == params.SnowDistribution.lognormal(0.4, 10)
        assert job.spring == params.SyntheticSpring()
        assert job.frames == params.FrameSpec()
        assert job.frames.every_hours == 24
        assert job.initial_swe_mm == 600.0

    def test_the_json_dump_loads_back_to_the_same_job(self, params: ModuleType) -> None:
        job = self.job(params, initial_swe_mm=250.0, frames=params.FrameSpec(every_hours=6))
        assert params.SnowJob.model_validate_json(job.model_dump_json()) == job

    def test_a_frame_interval_between_steps_is_refused(self, params: ModuleType) -> None:
        spring = params.SyntheticSpring(step_minutes=90)
        refused_on(("frames",), lambda: self.job(params, spring=spring))

    def test_a_frame_interval_of_whole_steps_is_accepted(self, params: ModuleType) -> None:
        spring = params.SyntheticSpring(step_minutes=90)
        self.job(params, spring=spring, frames=params.FrameSpec(every_hours=3))

    @pytest.mark.parametrize("swe", [-1.0, float("nan"), float("inf")])
    def test_a_bad_initial_swe_is_refused(self, params: ModuleType, swe: float) -> None:
        refused_on(("initial_swe_mm",), lambda: self.job(params, initial_swe_mm=swe))

    def test_a_frame_interval_below_1_hour_is_refused(self, params: ModuleType) -> None:
        refused_on(("every_hours",), lambda: params.FrameSpec(every_hours=0))
