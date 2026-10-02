"""G10, the suggested CRS (increment 15c-2, D8).

`docs/increments/15c-geographic-dem.md`, D8 and "Tests for @tester" G10.
`crs.suggest_crs(box, geographic_crs) -> CrsSuggestion`, frozen Pydantic with
`proj: str`, `family: str`, `max_scale_error: float`, `max_areal_error: float`.

PINNED HERE, where D8 is silent: `box` is a tuple `(W, E, S, N)` in degrees,
the order D8 and its table write it; `geographic_crs` anything `parse_crs`
reads. The family is read from `+proj=` (stere, tmerc, merc, lcc), not from
`family`'s wording, which D8 leaves free.

THE ORACLE shares nothing with `get_factors`: finite differences of the
suggestion's forward transform (pyproj `Transformer`, `always_xy`) over short
geodesic steps (pyproj `Geod`, central, 1 m each way) north and east of each
of D8's 21 x 21 points, edges included. Meridional scale `h` and parallel
scale `k` are projected length over ground length; the areal scale is the
projected parallelogram's area over the ground one's. It is shown able to
fail: the same suggestion with its scale factor replaced by 1 differs by
about half the unit-scale figure on the basin box (2.5e-3).

HOW THIS FILE GOES RED: `suggest_crs` does not exist; every test fails at
the fixture with `AttributeError`.
"""

from __future__ import annotations

import math
import re
import warnings
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
from pyproj import CRS, Geod, Transformer

import tin_engine.crs as crs_module

Box = tuple[float, float, float, float]  # (W, E, S, N)

#: D8's table. The prototype's balanced `k` is pinned where it was recorded.
VELHAS: Box = (-44.674, -43.459, -20.456, -18.457)
BASIN: Box = (-47.7, -36.3, -21.2, -7.0)
NORWAY: Box = (4.5, 31.2, 57.9, 71.2)
EQUATORIAL: Box = (-70.0, -50.0, -5.0, 5.0)
SVALBARD: Box = (10.0, 34.0, 76.4, 80.9)
EUROPE: Box = (-10.0, 30.0, 40.0, 60.0)


@pytest.fixture
def suggest() -> Callable[..., Any]:
    fn: Callable[..., Any] = crs_module.suggest_crs  # type: ignore[attr-defined]
    return fn


def proj4(proj: str) -> str:
    """The suggestion as PROJ.4 text, whatever form `proj` takes (a PROJ
    string, or WKT2 when the DEM's datum cannot be named in one)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # "loses information": the datum
        return str(CRS.from_user_input(proj).to_proj4())


def params(proj: str) -> dict[str, str]:
    """`+key=value` pairs of the suggestion's PROJ.4 text; a bare `+flag` maps to ''."""
    found = dict(re.findall(r"\+(\w+)(?:=(\S+))?", proj4(proj)))
    return {k: v for k, v in found.items()}


def scale_of(p: dict[str, str]) -> float:
    return float(p["k"] if "k" in p else p["k_0"])


def with_scale(proj: str, value: float) -> str:
    return re.sub(r"\+(k|k_0)=\S+", lambda m: f"+{m.group(1)}={value}", proj4(proj))


def oracle(proj: str, box: Box) -> tuple[float, float, float, float, float]:
    """`(kmin, kmax, max_scale_error, max_areal_error, max_angle_error)` over
    D8's 21 x 21 points; the last is the largest departure, in radians, of
    the angle between the projected north and east steps from a right angle."""
    west, east, south, north = box
    lon, lat = np.meshgrid(np.linspace(west, east, 21), np.linspace(south, north, 21))
    lon, lat = lon.ravel(), lat.ravel()
    geod: Geod = CRS.from_user_input(proj).get_geod()
    forward = Transformer.from_crs(CRS.from_user_input(proj).geodetic_crs, proj, always_xy=True)
    d = 1.0
    steps: list[np.ndarray] = []
    for azimuth in (0.0, 180.0, 90.0, 270.0):
        x, y, _ = geod.fwd(lon, lat, np.full(lon.size, azimuth), np.full(lon.size, d))
        px, py = forward.transform(x, y)
        steps.append(np.column_stack([px, py]))
    north_v = (steps[0] - steps[1]) / (2 * d)
    east_v = (steps[2] - steps[3]) / (2 * d)
    h = np.hypot(north_v[:, 0], north_v[:, 1])
    k = np.hypot(east_v[:, 0], east_v[:, 1])
    areal = np.abs(east_v[:, 0] * north_v[:, 1] - east_v[:, 1] * north_v[:, 0])
    scales = np.concatenate([h, k])
    scale_error = float(np.abs(scales - 1).max())
    angle = np.arcsin(np.clip(areal / (h * k), -1.0, 1.0))
    return (
        float(scales.min()),
        float(scales.max()),
        scale_error,
        float(np.abs(areal - 1).max()),
        float(np.abs(angle - np.pi / 2).max()),
    )


def centre(box: Box) -> tuple[float, float]:
    west, east, south, north = box
    return (west + east) / 2, (south + north) / 2


# ------------------------------------------------------------------ families


def taller_or_wider(lat_c: float, ground_width: float) -> Box:
    """A box 2° tall centred on `lat_c`, `ground_width` degrees wide on the ground."""
    width = ground_width / math.cos(math.radians(lat_c))
    return (10.0, 10.0 + width, lat_c - 1.0, lat_c + 1.0)


FAMILIES = [
    pytest.param((0.0, 20.0, 68.0, 72.0), "stere", id="polar_at_70"),
    pytest.param((0.0, 20.0, 67.9, 71.9), "lcc", id="wide_at_69.9"),
    pytest.param((0.0, 20.0, -72.0, -68.0), "stere", id="south_polar_at_minus_70"),
    pytest.param((0.0, 40.0, 58.0, 80.0), "stere", id="reaches_80"),
    pytest.param((0.0, 40.0, 58.0, 79.9), "tmerc", id="reaches_79.9"),
    pytest.param(taller_or_wider(45.0, 1.99), "tmerc", id="taller_by_0.01"),
    pytest.param(taller_or_wider(45.0, 2.01), "lcc", id="wider_by_0.01"),
    pytest.param((0.0, 40.0, 10.0, 20.0), "merc", id="wide_at_15"),
    pytest.param((0.0, 40.0, 10.1, 20.1), "lcc", id="wide_at_15.1"),
    pytest.param((0.0, 40.0, -20.0, -10.0), "merc", id="wide_at_minus_15"),
]


@pytest.mark.parametrize(("box", "family"), FAMILIES)
def test_the_family_and_its_rounded_parameters(
    suggest: Callable[..., Any], box: Box, family: str
) -> None:
    got = suggest(box, "EPSG:4326")
    p = params(got.proj)
    assert p["proj"] == family, got.proj
    assert got.family, "the family's name is printed in the refusal"
    lon_c, lat_c = centre(box)
    assert float(p["lon_0"]) == round(lon_c, 1)
    s = scale_of(p)
    assert s == round(s, 6) and 0.9 < s <= 1.0
    _, _, south, north = box
    if family == "stere":
        assert float(p["lat_0"]) == math.copysign(90.0, lat_c)
    elif family == "tmerc":
        assert float(p["lat_0"]) == 0.0
    elif family == "lcc":
        sixth = (north - south) / 6
        assert float(p["lat_1"]) == round(south + sixth, 1)
        assert float(p["lat_2"]) == round(north - sixth, 1)
        assert float(p["lat_0"]) == round(lat_c, 1)


@pytest.mark.parametrize("epsg", ["EPSG:4674", "EPSG:4326"])
def test_the_datum_is_the_dems_own(suggest: Callable[..., Any], epsg: str) -> None:
    """D8: "The datum is the DEM's own". Amended after `@perf`'s 15c-2
    acceptance: an ellipsoid alone (`+ellps=`) names no datum, so pyproj
    joined the two by a "Ballpark geographic offset". The suggestion's
    geodetic datum is the DEM's, and the DEM-to-suggestion transformation
    is no ballpark one."""
    got = suggest(VELHAS, epsg)
    suggested, dem = CRS(got.proj), CRS(epsg)
    assert suggested.is_projected
    assert suggested.ellipsoid.name == dem.ellipsoid.name
    assert suggested.datum is not None and dem.datum is not None
    assert suggested.datum.name == dem.datum.name, (suggested.datum.name, dem.datum.name)
    described = Transformer.from_crs(dem, suggested, always_xy=True).description
    assert "ballpark" not in described.lower(), described


@pytest.mark.parametrize(
    ("box", "lon_0", "k"),
    [(VELHAS, -44.1, 0.999972), (BASIN, -42.0, 0.997542), (NORWAY, 17.9, 0.996173)],
    ids=["velhas", "basin", "norway"],
)
def test_the_tables_transverse_mercator_rows(
    suggest: Callable[..., Any], box: Box, lon_0: float, k: float
) -> None:
    p = params(suggest(box, "EPSG:4326").proj)
    assert (p["proj"], float(p["lon_0"]), scale_of(p)) == ("tmerc", lon_0, k)


# ------------------------------------------------------------------ the oracle

BOXES = [
    pytest.param(VELHAS, id="velhas"),
    pytest.param(BASIN, id="basin"),
    pytest.param(NORWAY, id="norway"),
    pytest.param(EQUATORIAL, id="equatorial"),
    pytest.param(SVALBARD, id="svalbard"),
    pytest.param(EUROPE, id="europe"),
]

#: D8 asks 1e-6; central differences over 1 m reach it with margin.
AGREE = 1e-6


@pytest.mark.parametrize("box", BOXES)
def test_the_fields_agree_with_finite_differences(suggest: Callable[..., Any], box: Box) -> None:
    got = suggest(box, "EPSG:4326")
    _, _, scale_error, areal_error, _ = oracle(got.proj, box)
    assert got.max_scale_error == pytest.approx(scale_error, abs=AGREE)
    assert got.max_areal_error == pytest.approx(areal_error, abs=AGREE)


@pytest.mark.parametrize("box", BOXES)
def test_the_scale_is_balanced(suggest: Callable[..., Any], box: Box) -> None:
    """The largest and smallest point scales sit at 1 +- max_scale_error, to
    within the rounding of `s` to 6 decimals (5e-7, so 1e-6 between the two)."""
    got = suggest(box, "EPSG:4326")
    kmin, kmax, _, _, _ = oracle(got.proj, box)
    assert kmax - 1 == pytest.approx(got.max_scale_error, abs=1e-6 + AGREE)
    assert 1 - kmin == pytest.approx(got.max_scale_error, abs=1e-6 + AGREE)


def test_the_oracle_can_fail(suggest: Callable[..., Any]) -> None:
    """The basin's suggestion with `s` replaced by 1: the oracle's figure moves
    by about half the unit-scale one (0.49 % against 0.25 %)."""
    got = suggest(BASIN, "EPSG:4326")
    _, _, unit_error, _, _ = oracle(with_scale(got.proj, 1.0), BASIN)
    assert abs(unit_error - got.max_scale_error) > 2e-3
    assert got.max_scale_error == pytest.approx(0.0025, abs=2e-4)
    assert got.max_areal_error == pytest.approx(0.0049, abs=3e-4)


def test_every_suggestion_is_conformal(suggest: Callable[..., Any]) -> None:
    """Angular distortion under 3e-6 in every row of D8's table: the angle
    between the projected north and east steps stays 90 degrees."""
    for box in (VELHAS, BASIN, NORWAY, EQUATORIAL, SVALBARD, EUROPE):
        *_, angle_error = oracle(suggest(box, "EPSG:4326").proj, box)
        assert angle_error < 3e-6, box
