"""The snow run's job: `SyntheticSpring`, `FrameSpec` and `SnowJob` in `snow/params.py`.

`docs/increments/46-hbv-snow.md`, *4.1 The blueprint* (flags -> `SnowJob`, a
frozen Pydantic model holding `HbvSnowParameters`, `SnowDistribution`,
`SunCorrection`, `SyntheticSpring`, `FrameSpec`, `Resources`), *4.5 The
synthetic case* (2027-04-01T00:00Z, 91 days, hourly steps, S = 4, frames
every 24 hours, 600 mm initial SWE) and *4.8 The Python surface* (times are
timezone-aware UTC; a naive one is refused). `snow.json` records the job as
the model dumped, so the dump must load back to an equal job.

Interface:

- `SyntheticSpring(start, days, step_minutes, samples_per_step)`; `days`,
  `step_minutes` and `samples_per_step` at least 1.
- `FrameSpec(every_hours=24)`, at least 1.
- `SnowJob(hbv, distribution, correction, spring, frames, resources,
  initial_swe_mm=600.0)`; `hbv`, `spring` and `frames` default to their
  models' defaults, `distribution` to `SnowDistribution.lognormal(0.4, 10)`;
  `correction` and `resources` are required. A frame interval that is not a
  whole number of steps is refused, located on `frames`.
- Also refused on `frames` (code review round 1 of PR 3): a spring that is
  not a whole number of frames, and a run whose last frame is numbered past
  9999. Frames are numbered from 0, one at the start and one every
  `every_hours` to the end, so a run has days * 24 / every_hours + 1 of them;
  the case file's `swe****.var` pattern has four digits, so at most 10,000
  (numbers 0000 to 9999).
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
        # 2 h is 120 minutes, 1 1/3 steps of 90; the default 24 h is exactly 16.
        spring = params.SyntheticSpring(step_minutes=90)
        frames = params.FrameSpec(every_hours=2)
        caught = refused_on(("frames",), lambda: self.job(params, spring=spring, frames=frames))
        assert "a frame every 2 h is not a whole number of 90-minute steps" in str(caught)

    def test_the_default_daily_frame_on_90_minute_steps_is_accepted(
        self, params: ModuleType
    ) -> None:
        self.job(params, spring=params.SyntheticSpring(step_minutes=90))

    def test_a_frame_interval_of_whole_steps_is_accepted(self, params: ModuleType) -> None:
        spring = params.SyntheticSpring(step_minutes=90)
        self.job(params, spring=spring, frames=params.FrameSpec(every_hours=3))

    def test_a_spring_that_is_not_whole_frames_is_refused(self, params: ModuleType) -> None:
        # 91 days is 2184 h, 45 1/2 frames of 48 h: the last half frame was dropped.
        spring = params.SyntheticSpring(days=91)
        frames = params.FrameSpec(every_hours=48)
        caught = refused_on(("frames",), lambda: self.job(params, spring=spring, frames=frames))
        assert "91 days is not a whole number of 48-hour frames" in str(caught)

    def test_a_day_that_is_not_whole_frames_is_refused(self, params: ModuleType) -> None:
        # 24 h is 4 4/5 frames of 5 h; 5 h is 5 whole hourly steps, so only the
        # new rule refuses it.
        spring = params.SyntheticSpring(days=1)
        frames = params.FrameSpec(every_hours=5)
        refused_on(("frames",), lambda: self.job(params, spring=spring, frames=frames))

    def test_a_spring_of_whole_frames_longer_than_a_day_is_accepted(
        self, params: ModuleType
    ) -> None:
        self.job(
            params, spring=params.SyntheticSpring(days=90), frames=params.FrameSpec(every_hours=48)
        )

    def test_more_frames_than_four_digits_can_number_are_refused(self, params: ModuleType) -> None:
        # 417 days of hourly frames: frames 0 to 10,008, 10,009 in all.
        spring = params.SyntheticSpring(days=417)
        frames = params.FrameSpec(every_hours=1)
        caught = refused_on(("frames",), lambda: self.job(params, spring=spring, frames=frames))
        assert (
            "417 days with a frame every 1 h is 10,009 frames, "
            "more than the 10,000 the animation can number" in str(caught)
        )

    def test_a_last_frame_numbered_10000_is_refused(self, params: ModuleType) -> None:
        # 5000 days of 12-hour frames: frames 0 to 10,000, one past the pattern.
        spring = params.SyntheticSpring(days=5000)
        frames = params.FrameSpec(every_hours=12)
        refused_on(("frames",), lambda: self.job(params, spring=spring, frames=frames))

    def test_a_last_frame_numbered_9999_is_accepted(self, params: ModuleType) -> None:
        # 3333 days of 8-hour frames: frames 0 to 9999, exactly 10,000.
        spring = params.SyntheticSpring(days=3333)
        self.job(params, spring=spring, frames=params.FrameSpec(every_hours=8))

    @pytest.mark.parametrize("swe", [-1.0, float("nan"), float("inf")])
    def test_a_bad_initial_swe_is_refused(self, params: ModuleType, swe: float) -> None:
        refused_on(("initial_swe_mm",), lambda: self.job(params, initial_swe_mm=swe))

    def test_a_frame_interval_below_1_hour_is_refused(self, params: ModuleType) -> None:
        refused_on(("every_hours",), lambda: params.FrameSpec(every_hours=0))


class TestShadows:
    """`SnowJob.shadows` (section 10.5): the radiation design's frozen
    `RadiationRequest`, not a copy of its fields; None (the default) runs the
    sun correction under an open sky."""

    def job(self, params: ModuleType, **kw: Any) -> Any:
        return params.SnowJob(
            correction=params.SunCorrection(mode="melt_factor"), resources=RESOURCES, **kw
        )

    def test_the_default_is_no_shadows(self, params: ModuleType) -> None:
        assert self.job(params).shadows is None

    def test_it_holds_the_radiation_request(self, params: ModuleType) -> None:
        from tin_engine.radiation.request import RadiationRequest

        request = RadiationRequest(limb="upper", d_near_m=5000.0, d_far_m=5000.0)
        job = self.job(params, shadows=request)
        assert job.shadows == request
        assert isinstance(job.shadows, RadiationRequest)

    def test_the_dump_loads_back_to_an_equal_job(self, params: ModuleType) -> None:
        from tin_engine.radiation.request import RadiationRequest

        job = self.job(params, shadows=RadiationRequest())
        dumped = job.model_dump(mode="json")
        assert dumped["shadows"]["n_directions"] == 64
        assert params.SnowJob.model_validate(dumped) == job
        assert self.job(params).model_dump(mode="json")["shadows"] is None

    def test_a_bad_request_is_refused_on_shadows(self, params: ModuleType) -> None:
        refused_on(("shadows",), lambda: self.job(params, shadows={"n_directions": 48}))
