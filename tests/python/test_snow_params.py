"""The snow parameters in Pydantic: `tin_engine.snow.params`.

`docs/increments/46-hbv-snow.md`, *7. Red suites for `@tester`* (PR 1, D1),
with sections 4.2 ("`params.py` refuses the same values earlier"), 4.8 and 5.

Interface under test (names chosen here and listed in the handback):

- `HbvSnowParameters(tx, ts, cx, lw, cfr)`, frozen, section 5's defaults
  (0, 0, 3, 0.1, 0.05); cx, lw and cfr refused below 0; any non-finite value
  refused.
- `SnowDistribution(factors, deep=())`, frozen, `.factors` and `.deep`
  tuples of floats; `deep` empty (flat) or one deep-edge factor per bin with
  f_k <= deep_k <= 2 f_k, refused otherwise as the kernel refuses it
  (section 12.3, V3).
- `SnowDistribution.lognormal(cv, bins)`: for each of N equal-probability
  bands of the log-normal of mean 1, the band's mean as the factor and
  mean + sqrt(3) sd as the deep edge (section 12.4: D1, D3), divided
  by the factors' mean; CV 0 gives the one factor 1.0 and no deep; bins
  outside 1..64, a negative or non-finite CV, and a CV past the largest with
  no negative shallow edge refused, naming cv and the limit.
- `SunCorrection(mode, r, d)`, frozen: mode one of "none", "melt_factor",
  "threshold"; r and d default to 0.015 and 0.02; non-finite r or d refused.

A refused model field is a Pydantic `ValidationError` located on that field;
a refused `lognormal` argument a `ValueError` naming it.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from types import ModuleType
from typing import Any

import pydantic
import pytest

# Section 12.4's probe lines (docs/increments/46-probes/probe_sca.out), cv 0.4, ten bins.
D1 = (0.4778, 0.6216, 0.7157, 0.8004, 0.8848, 0.9750, 1.0780, 1.2060, 1.3902, 1.8505)
DEEP = (0.5964, 0.6736, 0.7592, 0.8421, 0.9280, 1.0226, 1.1343, 1.2796, 1.5074, 2.4245)
#: Section 12.4's largest cv with no negative shallow edge, by bins (probe, to 4 decimals).
LIMITS = {1: 0.5774, 2: 0.8282, 4: 1.0130, 10: 1.2412, 64: 1.7121}
BETWEEN = "a bin's deep edge must lie between its mean and twice its mean"


def refused_on(field: str, make: Callable[[], object]) -> None:
    with pytest.raises(pydantic.ValidationError) as caught:
        make()
    assert [e["loc"] for e in caught.value.errors()] == [(field,)]


@pytest.fixture(scope="module")
def params() -> ModuleType:
    import tin_engine.snow.params as params

    return params


class TestHbvSnowParameters:
    def test_defaults_are_section_5s(self, params: ModuleType) -> None:
        p = params.HbvSnowParameters()
        assert (p.tx, p.ts, p.cx, p.lw, p.cfr) == (0.0, 0.0, 3.0, 0.1, 0.05)

    @pytest.mark.parametrize("field", ["cx", "lw", "cfr"])
    def test_negative_refused(self, params: ModuleType, field: str) -> None:
        refused_on(field, lambda: params.HbvSnowParameters(**{field: -0.1}))

    @pytest.mark.parametrize("field", ["tx", "ts", "cx", "lw", "cfr"])
    @pytest.mark.parametrize("value", [math.nan, math.inf])
    def test_non_finite_refused(self, params: ModuleType, field: str, value: float) -> None:
        refused_on(field, lambda: params.HbvSnowParameters(**{field: value}))

    def test_zero_is_allowed(self, params: ModuleType) -> None:
        p = params.HbvSnowParameters(cx=0.0, lw=0.0, cfr=0.0)
        assert (p.cx, p.lw, p.cfr) == (0.0, 0.0, 0.0)

    def test_frozen(self, params: ModuleType) -> None:
        p = params.HbvSnowParameters()
        with pytest.raises(pydantic.ValidationError):
            p.cx = 4.0


class TestSnowDistribution:
    def test_d1_lognormal_factors(self, params: ModuleType) -> None:
        # The band means (section 12.4), not the midpoint quantiles.
        f = params.SnowDistribution.lognormal(0.4, 10).factors
        assert isinstance(f, tuple) and len(f) == 10
        assert f == pytest.approx(D1, abs=1e-4)
        assert sum(f) / len(f) == pytest.approx(1.0, abs=1e-12)

    @pytest.mark.parametrize("bins", [1, 10, 64])
    def test_cv_zero_is_one_factor(self, params: ModuleType, bins: int) -> None:
        d = params.SnowDistribution.lognormal(0.0, bins)
        assert (d.factors, d.deep) == ((1.0,), ())

    def test_d3_ten_bins_deep_edges(self, params: ModuleType) -> None:
        d = params.SnowDistribution.lognormal(0.4, 10)
        assert isinstance(d.deep, tuple) and len(d.deep) == 10
        assert d.deep == pytest.approx(DEEP, abs=1e-4)
        assert all(f <= h <= 2.0 * f for f, h in zip(d.factors, d.deep, strict=True))

    def test_d3_one_bin_is_uniform_of_the_given_cv(self, params: ModuleType) -> None:
        d = params.SnowDistribution.lognormal(0.4, 1)
        assert d.factors == (1.0,)
        assert d.deep == pytest.approx((1.0 + math.sqrt(3.0) * 0.4,), abs=1e-12)
        assert d != params.SnowDistribution.lognormal(0.0, 1)

    def test_d3_a_cv_too_large_for_one_bin(self, params: ModuleType) -> None:
        with pytest.raises(ValueError, match="cv") as caught:
            params.SnowDistribution.lognormal(0.58, 1)
        assert "0.577" in str(caught.value)

    @pytest.mark.parametrize(("bins", "limit"), sorted(LIMITS.items()))
    def test_d3_the_largest_cv_by_bins(self, params: ModuleType, bins: int, limit: float) -> None:
        d = params.SnowDistribution.lognormal(limit - 0.005, bins)
        assert min(2.0 * f - h for f, h in zip(d.factors, d.deep, strict=True)) >= 0.0
        with pytest.raises(ValueError, match="cv") as caught:
            params.SnowDistribution.lognormal(limit + 0.005, bins)
        assert f"{limit:.3f}" in str(caught.value)

    @pytest.mark.parametrize("bins", [0, 65])
    def test_bins_refused(self, params: ModuleType, bins: int) -> None:
        with pytest.raises(ValueError, match="bins"):
            params.SnowDistribution.lognormal(0.4, bins)

    @pytest.mark.parametrize("cv", [-0.1, math.nan, math.inf])
    def test_cv_refused(self, params: ModuleType, cv: float) -> None:
        with pytest.raises(ValueError, match="cv"):
            params.SnowDistribution.lognormal(cv, 10)

    # 1.7 is under 64 bins' limit of 1.7121 (section 12.4); the deep edges go
    # to the kernel with the factors.
    @pytest.mark.parametrize(("cv", "bins"), [(0.4, 10), (1.7, 64), (0.05, 2), (0.4, 1)])
    def test_the_kernel_accepts_them(self, params: ModuleType, cv: float, bins: int) -> None:
        import tin_engine._core as core

        d = params.SnowDistribution.lognormal(cv, bins)
        kernel = core.SnowDistribution(d.factors, d.deep)
        assert list(kernel.factors) == pytest.approx(d.factors, abs=1e-12)
        assert list(kernel.deep) == pytest.approx(d.deep, abs=1e-12)

    @pytest.mark.parametrize(
        "deep",
        [
            (0.5, 0.8, 1.2),
            (0.5, 0.8, 1.2, 1.5, 3.0),
            (0.5, 0.8, 0.99 * 1.2, 1.5),
            (0.5, 0.8, 1.2, 2.01 * 1.5),
            (0.5, math.nan, 1.2, 1.5),
            (0.5, 0.8, 1.2, math.inf),
        ],
    )
    def test_v3_a_bad_deep_is_refused(self, params: ModuleType, deep: tuple[float, ...]) -> None:
        with pytest.raises(pydantic.ValidationError, match="deep"):
            params.SnowDistribution(factors=(0.5, 0.8, 1.2, 1.5), deep=deep)

    def test_v3_a_deep_edge_off_its_bounds_says_so(self, params: ModuleType) -> None:
        with pytest.raises(pydantic.ValidationError, match=BETWEEN):
            params.SnowDistribution(factors=(0.5, 1.5), deep=(0.5, 3.1))

    def test_v3_the_mean_and_twice_the_mean_are_accepted(self, params: ModuleType) -> None:
        d = params.SnowDistribution(factors=(0.0, 0.5, 1.0, 2.5), deep=(0.0, 0.5, 2.0, 5.0))
        assert d.deep == (0.0, 0.5, 2.0, 5.0)
        assert params.SnowDistribution(factors=(1.0,)).deep == ()


class TestSunCorrection:
    def test_defaults_are_section_5s(self, params: ModuleType) -> None:
        s = params.SunCorrection(mode="melt_factor")
        assert (s.mode, s.r, s.d) == ("melt_factor", 0.015, 0.02)

    @pytest.mark.parametrize("mode", ["none", "melt_factor", "threshold"])
    def test_modes(self, params: ModuleType, mode: str) -> None:
        assert params.SunCorrection(mode=mode).mode == mode

    def test_unknown_mode_refused(self, params: ModuleType) -> None:
        refused_on("mode", lambda: params.SunCorrection(mode="melt-factor"))

    @pytest.mark.parametrize("field", ["r", "d"])
    def test_non_finite_refused(self, params: ModuleType, field: str) -> None:
        kwargs: dict[str, Any] = {"mode": "threshold", field: math.nan}
        refused_on(field, lambda: params.SunCorrection(**kwargs))
