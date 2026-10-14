"""The snow parameters in Pydantic: `tin_engine.snow.params`.

`docs/increments/46-hbv-snow.md`, *7. Red suites for `@tester`* (PR 1, D1),
with sections 4.2 ("`params.py` refuses the same values earlier"), 4.8 and 5.

Interface assumed (names chosen here and listed in the handback):

- `HbvSnowParameters(tx, ts, cx, lw, cfr)`, frozen, section 5's defaults
  (0, 0, 3, 0.1, 0.05); cx, lw and cfr refused below 0; any non-finite value
  refused.
- `SnowDistribution.lognormal(cv, bins)`, frozen, `.factors` a tuple of
  floats: the log-normal's quantiles at the bin midpoints (k + 1/2)/N,
  divided by their mean; CV 0 gives the one factor 1.0; bins outside 1..64 and
  a negative or non-finite CV refused.
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

D1 = (0.4974, 0.6287, 0.7228, 0.8080, 0.8930, 0.9838, 1.0873, 1.2154, 1.3973, 1.7664)


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
        f = params.SnowDistribution.lognormal(0.4, 10).factors
        assert isinstance(f, tuple) and len(f) == 10
        assert f == pytest.approx(D1, abs=1e-4)
        assert sum(f) / len(f) == pytest.approx(1.0, abs=1e-12)

    @pytest.mark.parametrize("bins", [1, 10, 64])
    def test_cv_zero_is_one_factor(self, params: ModuleType, bins: int) -> None:
        assert params.SnowDistribution.lognormal(0.0, bins).factors == (1.0,)

    @pytest.mark.parametrize("bins", [0, 65])
    def test_bins_refused(self, params: ModuleType, bins: int) -> None:
        with pytest.raises(ValueError, match="bins"):
            params.SnowDistribution.lognormal(0.4, bins)

    @pytest.mark.parametrize("cv", [-0.1, math.nan, math.inf])
    def test_cv_refused(self, params: ModuleType, cv: float) -> None:
        with pytest.raises(ValueError, match="cv"):
            params.SnowDistribution.lognormal(cv, 10)

    @pytest.mark.parametrize(("cv", "bins"), [(0.4, 10), (2.0, 64), (0.05, 2)])
    def test_the_kernel_accepts_them(self, params: ModuleType, cv: float, bins: int) -> None:
        import tin_engine._core as core

        f = params.SnowDistribution.lognormal(cv, bins).factors
        assert list(core.SnowDistribution(f).factors) == pytest.approx(f, abs=1e-12)


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
