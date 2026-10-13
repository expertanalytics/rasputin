"""The HBV snow kernel through the binding: `_core.SnowDistribution`, `SnowState`.

`docs/increments/46-hbv-snow.md`, *7. Red suites for `@tester`* (PR 1): the
hand-computed steps H1-H12 that the C++ suite also runs, the invariants I1-I5
as Hypothesis properties P1-P5, the refusals V1 as `ValueError`s naming the
status, the renormalisation V2, and the GIL released during a step.

Interface assumed (names chosen here and listed in the handback; they mirror
the C++ in section 4.2):

- `_core.SunCorrection`: `none`, `melt_factor`, `threshold`.
- `_core.SnowDistribution(factors)`: a 1-D float sequence; `.factors`, a
  float64 array, renormalised to mean 1. Refusals are `ValueError`s carrying
  the C++ message.
- `_core.SnowState.uniform(triangles, distribution, initial_swe_mm)`;
  read-only properties `swe`, `sca`, `outflow` (float64, shape `(n,)`) and
  `sp`, `sw` (float64, shape `(n, N)`).
- `state.step(*, tx, ts, cx, lw, cfr, mode, r, d, temperature, precipitation,
  radiation_anomaly, dt_seconds, threads)`: returns None; a refused step
  raises `ValueError` whose message names the status (`invalid_parameters`,
  ...) and, for a non-finite value, the first triangle (`triangle 2`), and
  leaves the state as it was. `radiation_anomaly` is empty under `none`.
  Releases the GIL.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gil_probe import ticks_during

DAY = 86400.0
HBV = {"tx": 0.0, "ts": 0.0, "cx": 3.0, "lw": 0.1, "cfr": 0.05}  # section 5's defaults
SUN = {"r": 0.015, "d": 0.02}
FOUR = (0.5, 0.8, 1.2, 1.5)
TEN = (0.4974, 0.6287, 0.7228, 0.8080, 0.8930, 0.9838, 1.0873, 1.2154, 1.3973, 1.7664)
RELEASED_TICKS = 20  # as test_core_noding.py: a released call lets the ticker run

Array = npt.NDArray[np.float64]


@pytest.fixture(scope="module")
def core() -> Any:
    import tin_engine._core as core

    for name in ("SnowDistribution", "SnowState", "SunCorrection"):
        assert hasattr(core, name), f"tin_engine._core has no {name}"
    return core


def ten(core: Any) -> Any:
    f = np.array(TEN)
    return core.SnowDistribution(f / f.mean())


def step(
    core: Any,
    state: Any,
    t: float | Array,
    p: float | Array,
    di: float | Array = 0.0,
    *,
    mode: str = "none",
    dt: float = DAY,
    threads: int = 1,
    **hbv: float,
) -> None:
    """One step; scalars are spread over every triangle."""
    n = state.swe.shape[0]
    on = mode != "none"
    state.step(
        **{**HBV, **hbv},
        mode=getattr(core.SunCorrection, mode),
        **SUN,
        temperature=np.broadcast_to(np.asarray(t, dtype=np.float64), (n,)).copy(),
        precipitation=np.broadcast_to(np.asarray(p, dtype=np.float64), (n,)).copy(),
        radiation_anomaly=(
            np.broadcast_to(np.asarray(di, dtype=np.float64), (n,)).copy() if on else np.empty(0)
        ),
        dt_seconds=dt,
        threads=threads,
    )


def one_bin(core: Any, sp: float, sw: float = 0.0) -> Any:
    """One bin of sp, with sw of water from a priming rain at T = tx = ts = 0."""
    state = core.SnowState.uniform(1, core.SnowDistribution([1.0]), sp)
    if sw > 0.0:
        step(core, state, 0.0, sw)
    assert state.sp[0, 0] == pytest.approx(sp, abs=1e-12)
    assert state.sw[0, 0] == pytest.approx(sw, abs=1e-12)
    return state


def expect(state: Any, sp: float, sw: float, outflow: float, sca: float) -> None:
    assert state.sp[0, 0] == pytest.approx(sp, abs=1e-12)
    assert state.sw[0, 0] == pytest.approx(sw, abs=1e-12)
    assert state.outflow[0] == pytest.approx(outflow, abs=1e-12)
    assert state.swe[0] == pytest.approx(sp + sw, abs=1e-12)
    assert state.sca[0] == sca


# ---------------------------------------------------------------- H1-H12


class TestHandSteps:
    def test_h1_melt_then_h2_refreeze(self, core: Any) -> None:
        s = one_bin(core, 100.0)
        step(core, s, 5.0, 0.0)
        expect(s, 85.0, 8.5, 6.5, 1.0)
        step(core, s, -2.0, 0.0)  # refreeze 0.05 * 3 * 2 = 0.3
        expect(s, 85.3, 8.2, 0.0, 1.0)

    def test_h3_snowfall_by_the_factors(self, core: Any) -> None:
        s = core.SnowState.uniform(1, core.SnowDistribution(FOUR), 100.0)
        before = s.sp.copy()
        step(core, s, -1.0, 10.0)
        np.testing.assert_allclose(s.sp[0] - before[0], [5.0, 8.0, 12.0, 15.0], atol=1e-12)
        assert s.swe[0] == pytest.approx(110.0, abs=1e-12)

    def test_h4_rain_on_snow(self, core: Any) -> None:
        s = one_bin(core, 100.0)
        step(core, s, 2.0, 10.0)
        expect(s, 94.0, 9.4, 6.6, 1.0)

    def test_h5_rain_on_bare_ground(self, core: Any) -> None:
        s = core.SnowState.uniform(1, core.SnowDistribution([1.0]), 0.0)
        step(core, s, 2.0, 10.0)
        expect(s, 0.0, 0.0, 10.0, 0.0)

    def test_h5b_rain_on_a_bare_bin_below_ts(self, core: Any) -> None:
        s = core.SnowState.uniform(1, core.SnowDistribution([1.0]), 0.0)
        step(core, s, -0.5, 10.0, tx=-1.0)
        expect(s, 0.0, 0.0, 10.0, 0.0)

    def test_h6_melt_out_within_a_step(self, core: Any) -> None:
        s = one_bin(core, 10.0)
        step(core, s, 10.0, 0.0)
        expect(s, 0.0, 0.0, 10.0, 0.0)

    def test_h7_four_bins(self, core: Any) -> None:
        s = core.SnowState.uniform(1, core.SnowDistribution(FOUR), 100.0)
        step(core, s, 20.0, 0.0)
        np.testing.assert_allclose(s.sp[0], [0.0, 20.0, 60.0, 90.0], atol=1e-12)
        np.testing.assert_allclose(s.sw[0], [0.0, 2.0, 6.0, 9.0], atol=1e-12)
        assert s.sca[0] == 0.75
        assert s.swe[0] == pytest.approx(46.75, abs=1e-12)
        assert s.outflow[0] == pytest.approx(53.25, abs=1e-12)

    @pytest.mark.parametrize(
        ("mode", "melt"), [("melt_factor", 18.0), ("threshold", 18.0), ("none", 12.0)]
    )
    def test_h8_break_even(self, core: Any, mode: str, melt: float) -> None:
        s = one_bin(core, 100.0)
        step(core, s, 4.0, 0.0, 100.0, mode=mode)
        assert s.sp[0, 0] == pytest.approx(100.0 - melt, abs=1e-12)

    @pytest.mark.parametrize("mode", ["melt_factor", "threshold"])
    def test_h8b_shaded_slope_with_water(self, core: Any, mode: str) -> None:
        s = one_bin(core, 100.0, 5.0)
        step(core, s, 4.0, 0.0, -300.0, mode=mode)
        expect(s, 100.0, 5.0, 0.0, 1.0)

    def test_h9_sub_day_step(self, core: Any) -> None:
        s = one_bin(core, 100.0)
        step(core, s, 5.0, 0.0, dt=3600.0)
        assert s.sp[0, 0] == pytest.approx(99.375, abs=1e-12)

    def test_h10_rain_at_tx(self, core: Any) -> None:
        s = one_bin(core, 100.0)
        step(core, s, 0.0, 10.0)
        expect(s, 100.0, 10.0, 0.0, 1.0)

    def test_h11_refreeze_capped_by_the_water(self, core: Any) -> None:
        s = one_bin(core, 100.0, 1.0)
        step(core, s, -10.0, 0.0, ts=0.0)  # potential 1.5 > sw = 1
        expect(s, 101.0, 0.0, 0.0, 1.0)

    @pytest.mark.parametrize(
        ("mode", "sp", "sw"), [("threshold", 97.0, 3.0), ("melt_factor", 100.0, 0.0)]
    )
    def test_h12_melt_below_ts(self, core: Any, mode: str, sp: float, sw: float) -> None:
        s = one_bin(core, 100.0)
        step(core, s, -1.0, 0.0, 100.0, mode=mode)
        expect(s, sp, sw, 0.0, 1.0)

    def test_uniform_shapes_and_values(self, core: Any) -> None:
        s = core.SnowState.uniform(3, core.SnowDistribution(FOUR), 100.0)
        assert s.sp.shape == s.sw.shape == (3, 4)
        assert s.swe.shape == s.sca.shape == s.outflow.shape == (3,)
        assert all(a.dtype == np.float64 for a in (s.sp, s.sw, s.swe, s.sca, s.outflow))
        np.testing.assert_allclose(s.sp, np.tile([50.0, 80.0, 120.0, 150.0], (3, 1)), atol=1e-12)
        assert (s.sw == 0.0).all() and (s.sca == 1.0).all() and (s.outflow == 0.0).all()

    def test_zero_triangles(self, core: Any) -> None:
        s = core.SnowState.uniform(0, ten(core), 600.0)
        step(core, s, np.empty(0), np.empty(0), np.empty(0), mode="threshold")
        assert s.swe.shape == (0,)


# ---------------------------------------------------------------- P1-P5

N_TRIANGLES = 6
forcing = st.lists(
    st.tuples(
        st.lists(st.floats(-30.0, 30.0), min_size=N_TRIANGLES, max_size=N_TRIANGLES),
        st.lists(
            st.one_of(st.just(0.0), st.floats(0.0, 100.0)),
            min_size=N_TRIANGLES,
            max_size=N_TRIANGLES,
        ),
        st.lists(st.floats(-600.0, 600.0), min_size=N_TRIANGLES, max_size=N_TRIANGLES),
    ),
    min_size=1,
    max_size=40,
)
PROPERTY = settings(
    max_examples=60,
    derandomize=True,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)


def snapshot(s: Any) -> tuple[bytes, ...]:
    return tuple(a.tobytes() for a in (s.sp, s.sw, s.swe, s.sca, s.outflow))


class TestProperties:
    @PROPERTY
    @given(
        steps=forcing,
        initial=st.floats(0.0, 1e4),
        mode=st.sampled_from(["none", "melt_factor", "threshold"]),
    )
    def test_p1_p2_mass_and_bounds(
        self, core: Any, steps: list[Any], initial: float, mode: str
    ) -> None:
        """Scale: swe up to 1e4 mm, P up to 100 mm per step (section 5)."""
        s = core.SnowState.uniform(N_TRIANGLES, ten(core), initial)
        for t, p, di in steps:
            before = s.swe.copy()
            precip = np.array(p)
            step(core, s, np.array(t), precip, np.array(di), mode=mode)
            tol = 1e-9 * np.maximum(1.0, before + precip)
            assert (np.abs(s.swe - (before + precip - s.outflow)) <= tol).all()
            assert (s.sp >= 0.0).all() and (s.sw >= 0.0).all()
            assert (s.sw <= HBV["lw"] * s.sp + 1e-12).all()
            assert ((s.sca >= 0.0) & (s.sca <= 1.0)).all() and (s.outflow >= 0.0).all()

    @PROPERTY
    @given(steps=forcing)
    def test_p3_flat_is_plain(self, core: Any, steps: list[Any]) -> None:
        states = {
            m: core.SnowState.uniform(N_TRIANGLES, ten(core), 600.0)
            for m in ("none", "melt_factor", "threshold")
        }
        for t, p, _ in steps:
            for mode, s in states.items():
                step(core, s, np.array(t), np.array(p), 0.0, mode=mode)
            assert snapshot(states["none"]) == snapshot(states["melt_factor"])
            assert snapshot(states["none"]) == snapshot(states["threshold"])

    @PROPERTY
    @given(steps=forcing, hours=st.integers(1, 24))
    def test_p4_one_bin_is_hbv96(self, core: Any, steps: list[Any], hours: int) -> None:
        s = core.SnowState.uniform(N_TRIANGLES, core.SnowDistribution([1.0]), 600.0)
        for t, p, _ in steps:
            sp, sw = s.sp[:, 0].copy(), s.sw[:, 0].copy()
            step(core, s, np.array(t), np.array(p), dt=3600.0 * hours)
            for i in range(N_TRIANGLES):
                want = hbv96(sp[i], sw[i], t[i], p[i], hours / 24.0)
                got = (s.sp[i, 0], s.sw[i, 0], s.outflow[i])
                tol = 1e-12 * max(1.0, sp[i] + sw[i] + p[i])
                assert got == pytest.approx(want, abs=tol)
                assert s.sca[i] == (1.0 if want[0] > 0.0 else 0.0)

    @PROPERTY
    @given(steps=forcing, mode=st.sampled_from(["none", "melt_factor", "threshold"]))
    def test_p5_threads_do_not_change_a_bit(self, core: Any, steps: list[Any], mode: str) -> None:
        # Each step's six values tiled over 5 000 triangles: enough blocks that
        # eight workers start.
        reps = 5000 // N_TRIANGLES + 1
        one, eight = (
            core.SnowState.uniform(N_TRIANGLES * reps, ten(core), 600.0) for _ in range(2)
        )
        for t, p, di in steps:
            args = (np.tile(t, reps), np.tile(p, reps), np.tile(di, reps))
            step(core, one, *args, mode=mode, threads=1)
            step(core, eight, *args, mode=mode, threads=8)
            assert snapshot(one) == snapshot(eight)


def hbv96(sp: float, sw: float, t: float, p: float, days: float) -> tuple[float, float, float]:
    """Plain single-point HBV-96 (section 4.2's steps with one bin, mode none)."""
    snow = p if t < HBV["tx"] else 0.0
    rain, out = p - snow, 0.0
    sp += snow
    if sp > 0.0:
        sw += rain
    else:
        out += rain
    if t > HBV["ts"]:
        melt = min(HBV["cx"] * days * (t - HBV["ts"]), sp)
        sp, sw = sp - melt, sw + melt
    elif t < HBV["ts"]:
        refreeze = min(HBV["cfr"] * HBV["cx"] * days * (HBV["ts"] - t), sw)
        sp, sw = sp + refreeze, sw - refreeze
    out += max(0.0, sw - HBV["lw"] * sp)
    return sp, min(sw, HBV["lw"] * sp), out


# ---------------------------------------------------------------- V1, V2

NAN, INF = float("nan"), float("inf")


class TestRefusals:
    @pytest.mark.parametrize(
        ("factors", "reason"),
        [
            ([], "between 1 and 64"),
            ([1.0] * 65, "between 1 and 64"),
            ([0.5, NAN, 1.2, 1.5], "finite"),
            ([-0.1, 0.5, 1.2, 2.4], "negative"),
            ([1.5, 1.2, 0.8, 0.5], "ascending"),
            ([0.5, 0.8, 1.2, 1.5 + 4e-6], "mean"),
        ],
    )
    def test_distribution(self, core: Any, factors: list[float], reason: str) -> None:
        with pytest.raises(ValueError, match=reason):
            core.SnowDistribution(factors)

    @pytest.mark.parametrize("swe", [-1.0, NAN])
    def test_uniform(self, core: Any, swe: float) -> None:
        with pytest.raises(ValueError, match="initial"):
            core.SnowState.uniform(3, core.SnowDistribution(FOUR), swe)

    @pytest.fixture
    def used(self, core: Any) -> Any:
        """Three triangles, four bins, one ok step behind them."""
        s = core.SnowState.uniform(3, core.SnowDistribution(FOUR), 100.0)
        step(core, s, 1.0, 5.0)
        assert s.outflow[2] > 0.0 and s.sw[2, 3] > 0.0
        return s

    @pytest.mark.parametrize(
        ("change", "status"),
        [
            ({"cx": -0.1}, "invalid_parameters"),
            ({"lw": -0.1}, "invalid_parameters"),
            ({"cfr": -0.1}, "invalid_parameters"),
            *[({k: NAN}, "invalid_parameters") for k in ("tx", "ts", "cx", "lw", "cfr", "r", "d")],
            ({"temperature": [1.0, 1.0]}, "size_mismatch"),
            ({"temperature": [1.0] * 4}, "size_mismatch"),
            ({"precipitation": [2.0, 2.0]}, "size_mismatch"),
            ({"precipitation": [2.0] * 4}, "size_mismatch"),
            ({"mode": "none"}, "size_mismatch"),  # dI of length n under none
            ({"radiation_anomaly": []}, "size_mismatch"),
            ({"mode": "threshold", "radiation_anomaly": []}, "size_mismatch"),
            *[
                ({k: [v0, v0, bad]}, "non_finite_forcing")
                for k, v0 in (
                    ("temperature", 1.0),
                    ("precipitation", 2.0),
                    ("radiation_anomaly", 50.0),
                )
                for bad in (NAN, INF, -INF)
            ],
            ({"precipitation": [2.0, 2.0, -1e-12]}, "negative_precipitation"),
            *[({"dt_seconds": dt}, "bad_step_length") for dt in (0.0, -3600.0, NAN)],
        ],
    )
    def test_step(self, core: Any, used: Any, change: dict[str, Any], status: str) -> None:
        kwargs: dict[str, Any] = {
            **HBV,
            **SUN,
            "mode": "melt_factor",
            "temperature": [1.0, 1.0, 1.0],
            "precipitation": [2.0, 2.0, 2.0],
            "radiation_anomaly": [50.0, 50.0, 50.0],
            "dt_seconds": DAY,
            "threads": 1,
        }
        kwargs.update(change)
        kwargs["mode"] = getattr(core.SunCorrection, kwargs["mode"])
        for k in ("temperature", "precipitation", "radiation_anomaly"):
            kwargs[k] = np.array(kwargs[k], dtype=np.float64)
        before = snapshot(used)
        with pytest.raises(ValueError, match=status) as caught:
            used.step(**kwargs)
        if status == "non_finite_forcing":
            assert "triangle 2" in str(caught.value)
        assert snapshot(used) == before

    def test_v2_renormalised(self, core: Any) -> None:
        d = core.SnowDistribution([0.5, 0.8, 1.2, 1.5 + 4e-10])
        assert isinstance(d.factors, np.ndarray) and d.factors.dtype == np.float64
        assert abs(d.factors.mean() - 1.0) <= 1e-15


def test_step_releases_the_gil(core: Any) -> None:
    n = 400_000
    s = core.SnowState.uniform(n, ten(core), 600.0)
    t, p, di = np.full(n, 3.0), np.full(n, 1.0), np.full(n, 80.0)

    def run() -> None:
        for _ in range(5):
            step(core, s, t, p, di, mode="melt_factor", dt=3600.0, threads=1)

    _, ticks, elapsed = ticks_during(run)
    assert elapsed >= 0.05, f"took only {elapsed:.3f}s: too fast to see; enlarge the state"
    assert ticks >= RELEASED_TICKS, f"only {ticks} ticks in {elapsed:.3f}s: it holds the GIL"
