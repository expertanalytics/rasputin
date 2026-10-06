"""A projected CRS given by parameters (3072 = 32767), read as the EPSG code it is.

`docs/increments/geotiff-crs-by-parameters.md`; section and row numbers below
are its. Section 8 lists the red tests; each test here names its number.

The fixtures are micro-TIFFs from `geotiff_fixtures.py`: `AUSTRIA_KEYS` is the
eighteen GeoKeys of Ola's openDEM Austria file (section 1) on the usual 3 x 4
baseline, and each refusal fixture is that set with exactly one change.
`test_geotiff_fixtures.py` shows with tifffile alone that each fixture carries
its change. The 2 GB file itself is read only by red test 12, header only,
and only where it is present.

HOW THIS FILE GOES RED. As in `test_io_geotiff.py`, production modules are
imported inside fixtures, so each test fails on its own. At the design's base
every 3072 = 32767 file is refused with "ProjectedCSTypeGeoKey (3072) = 32767
is not a resolvable EPSG code": the accepting tests fail on that refusal, the
refusal tests fail on `startswith(P)` or a token their row names, and the
cache and substitution tests fail on the names `_epsg_matches_cached` and
`crs.epsg_matches`, which did not exist at the design's base.

Red tests 13 and 14 were added at code review round 1, red against the green
step `c8984e8`: there a NaN or infinite parameter escapes as pyproj's
`CRSError`, a two-double parameter as `TypeError`, a NaN or infinite 2057 is
refused as "disagrees" rather than "not one finite number", and a file in
grads is refused for its prime meridian (row 6) rather than its unit (row 8).

WHAT A REFUSAL MUST SAY. Every message of section 4's thirteen rows begins
with `P` exactly (case-sensitive `startswith`), then names its key by number
and name and the file's value; those tokens are matched case-insensitively,
so the suite pins what is named and not the sentence around it.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterator, Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pyproj
import pytest

from geotiff_fixtures import (
    EPSG_AUSTRIA_LAMBERT,
    EPSG_ETRS89,
    EPSG_MGI,
    EPSG_WGS84,
    EPSG_XIAN_1980,
    GEOGRAPHIC_TYPE,
    GRADS_WITH_PARIS_MERIDIAN,
    NOT_ONE_NUMBER_DEFECTS,
    PARAMETRIC_DEFECTS,
    PROJ_FALSE_ORIGIN_EASTING,
    PROJ_STD_PARALLEL_1,
    PROJ_STD_PARALLEL_2,
    PROJECTED_CS_TYPE,
    USER_DEFINED,
    GeoKeyValue,
    ParametricDefect,
    austria_keys,
    austria_tiff,
    elevations,
    lambert_on,
    micro_tiff,
    transverse_mercator_on,
    with_keys,
)

#: Section 4: the words every refusal of this path begins with.
P = "ProjectedCSTypeGeoKey (3072) = 32767 (a CRS given by parameters)"

#: Ola's file (section 1). Found by walking up from this file, as
#: `test_io_cog.py` finds DTM10, so the test runs from a worktree too.
_AUSTRIA_RELATIVE = Path("rasputin_data") / "austria_dgm10" / "dhm_at_lamb_10m_2018.tif"
AUSTRIA_FILE = next(
    (
        p / _AUSTRIA_RELATIVE
        for p in Path(__file__).resolve().parents
        if (p / _AUSTRIA_RELATIVE).is_file()
    ),
    None,
)


def _names_numbers(message: str, *numbers: int) -> None:
    """Each number appears in `message` as a whole number, not inside another."""
    for number in numbers:
        assert re.search(rf"(?<![\d.]){number}(?!\d)", message), f"{number} not in {message!r}"


@pytest.fixture
def geotiff() -> ModuleType:
    import tin_engine.io.geotiff as module

    return module


@pytest.fixture
def geotiff_error() -> type[Exception]:
    from tin_engine.io.models import GeoTiffError

    return GeoTiffError


@pytest.fixture
def epsg_of(geotiff: ModuleType) -> Callable[[Mapping[int, GeoKeyValue]], int | None]:
    """The EPSG code `read_meta` gives a baseline tile carrying these GeoKeys."""

    def read(geokeys: Mapping[int, GeoKeyValue]) -> int | None:
        meta = geotiff.read_meta(austria_tiff(geokeys))
        assert meta.crs == f"EPSG:{meta.epsg}"
        assert meta.geographic is False
        return meta.epsg  # type: ignore[no-any-return]

    return read


@pytest.fixture
def refused(geotiff: ModuleType, geotiff_error: type[Exception]) -> Callable[..., str]:
    """`read_meta` must raise `GeoTiffError` whose message starts with `P` and
    names each token (case-insensitive substring)."""

    def check(geokeys: Mapping[int, GeoKeyValue], *names: str) -> str:
        with pytest.raises(geotiff_error) as info:
            geotiff.read_meta(austria_tiff(geokeys))
        message = str(info.value)
        assert message.startswith(P), f"does not start with P: {message!r}"
        for token in names:
            assert token.lower() in message.lower(), f"{token!r} not named in {message!r}"
        return message

    return check


@pytest.fixture
def fresh_cache(geotiff: ModuleType) -> Iterator[Any]:
    """The cached matcher, cleared before the test and after it, so neither an
    earlier read nor a substituted matcher's result can leak across tests."""
    cached = geotiff._epsg_matches_cached
    cached.cache_clear()
    yield cached
    cached.cache_clear()


# ---------------------------------------------------------------------------
# Red tests 1 and 2: the Austrian GeoKeys read as EPSG:31287
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("entry", ["read_meta", "read_header", "read_page", "decode_dem"])
def test_austrian_keys_read_as_austria_lambert(geotiff: ModuleType, entry: str) -> None:
    """Red test 1: every entry point reads the file as if it carried 3072 = 31287."""
    result = getattr(geotiff, entry)(austria_tiff())
    meta = result if entry == "read_meta" else result.meta if entry == "decode_dem" else result[0]
    assert meta.epsg == EPSG_AUSTRIA_LAMBERT
    assert meta.crs == "EPSG:31287"
    assert meta.geographic is False
    if entry == "decode_dem":
        np.testing.assert_array_equal(result.array, elevations())


@pytest.mark.parametrize(("first", "second"), [(46.0, 49.0), (49.0, 46.0)], ids=["46_49", "49_46"])
def test_either_order_of_the_standard_parallels(
    epsg_of: Callable[..., int | None], first: float, second: float
) -> None:
    """Red test 2: PROJ's candidate search takes Lambert 2SP's two standard
    parallels in either order (section 2), so both orders read as 31287."""
    keys = austria_keys({PROJ_STD_PARALLEL_1: first, PROJ_STD_PARALLEL_2: second})
    assert epsg_of(keys) == EPSG_AUSTRIA_LAMBERT


# ---------------------------------------------------------------------------
# Red tests 3, 4 and 6: transverse Mercator, one CRS under two codes, the datum filter
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("keys", "expected"),
    [
        (transverse_mercator_on(EPSG_ETRS89, 15.0), 25833),
        (transverse_mercator_on(EPSG_WGS84, 15.0), 32633),
        (
            transverse_mercator_on(
                EPSG_MGI, 13.333333333333334, scale=1.0, false_easting=450000.0,
                false_northing=-5000000.0,
            ),
            31258,
        ),
    ],
    ids=["utm33_etrs89", "utm33_wgs84", "mgi_gauss_krueger_m31"],
)  # fmt: skip
def test_transverse_mercator_by_parameters(
    epsg_of: Callable[..., int | None], keys: dict[int, GeoKeyValue], expected: int
) -> None:
    """Red test 3: 15 E, 0.9996, 500 000, 0 on ETRS89 is 25833 and on WGS 84
    is 32633; MGI Gauss-Krueger M31 on 4312 is 31258."""
    assert epsg_of(keys) == expected


def test_one_crs_under_two_codes_reads_as_the_lower(epsg_of: Callable[..., int | None]) -> None:
    """Red test 4: Xian 1980 / Gauss-Kruger CM 75E is EPSG:2338 and 2370
    (section 6, `same_crs` to each other). Its parameters, taken from
    `CRS.from_epsg(2338)`, on 4610 read as the lower code."""
    params = {int(p.code): p.value for p in pyproj.CRS.from_epsg(2338).coordinate_operation.params}
    keys = transverse_mercator_on(
        EPSG_XIAN_1980,
        params[8802],
        latitude=params[8801],
        scale=params[8805],
        false_easting=params[8806],
        false_northing=params[8807],
    )
    assert epsg_of(keys) == 2338


def test_the_datum_filter_keeps_only_codes_on_the_files_datum(
    epsg_of: Callable[..., int | None],
) -> None:
    """Red test 6: UTM 35N by parameters on ETRS89 (4258). PROJ also calls
    EPSG:3067 (on EUREF-FIN) equivalent; unfiltered, the two would be
    ambiguous (they are not `same_crs`). Filtered, it is 25835."""
    assert epsg_of(transverse_mercator_on(EPSG_ETRS89, 27.0)) == 25835


# ---------------------------------------------------------------------------
# Red test 5: matching nothing stays refused, and PROJ's tolerance both sides
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("keys", "base"),
    [
        (lambert_on(EPSG_WGS84), EPSG_WGS84),
        (transverse_mercator_on(EPSG_ETRS89, 15.0, scale=1.0), EPSG_ETRS89),
        (austria_keys({PROJ_FALSE_ORIGIN_EASTING: 400000.001}), EPSG_MGI),
    ],
    ids=["austrian_on_wgs84", "tm15_scale_1_on_etrs89", "false_easting_plus_1mm"],
)
def test_a_crs_matching_no_code_stays_refused(
    refused: Callable[..., str], keys: dict[int, GeoKeyValue], base: int
) -> None:
    """Red test 5 (row 12). `tm15_scale_1_on_etrs89` is a case where PROJ offers
    six candidates and none is equivalent: candidates are judged, not trusted.
    `false_easting_plus_1mm` is outside PROJ's 1e-10 relative at 400 000 m
    (0.04 mm; section 5)."""
    message = refused(keys, "no EPSG projected CRS", f"EPSG:{base}")
    _names_numbers(message, 3072, 32767, base)


def test_a_false_easting_within_projs_tolerance_still_matches(
    epsg_of: Callable[..., int | None],
) -> None:
    """Red test 5's pair: + 1e-5 m on 400 000 m is 2.5e-11 relative, inside
    PROJ's 1e-10 (section 5), so it reads as 31287. Checked at this one
    magnitude, the Austrian false easting."""
    assert epsg_of(austria_keys({PROJ_FALSE_ORIGIN_EASTING: 400000.00001})) == 31287


# ---------------------------------------------------------------------------
# Red test 7: ambiguous stays refused (the one substitution of production code)
# ---------------------------------------------------------------------------


def test_matches_that_are_not_one_crs_are_refused(
    refused: Callable[..., str],
    geotiff: ModuleType,
    fresh_cache: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Red test 7 (row 13). No EPSG definition reaches this row at PROJ 9.8.1
    (section 6), so `crs.epsg_matches` is substituted to return EPSG:3067 and
    25835, which are not `same_crs`. Substituted wherever `io.geotiff` can see
    it: on `tin_engine.crs`, and on `io.geotiff` too if it imported the name."""
    import tin_engine.crs

    def two_crss(_crs: Any) -> tuple[int, ...]:
        return (3067, 25835)

    monkeypatch.setattr(tin_engine.crs, "epsg_matches", two_crss)
    if hasattr(geotiff, "epsg_matches"):
        monkeypatch.setattr(geotiff, "epsg_matches", two_crss)
    message = refused(transverse_mercator_on(EPSG_ETRS89, 27.0), "not the same CRS")
    _names_numbers(message, 3067, 25835)


# ---------------------------------------------------------------------------
# Red test 8: each check of section 4 on its one defect (rows 1 to 11)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("defect", PARAMETRIC_DEFECTS, ids=[d.name for d in PARAMETRIC_DEFECTS])
def test_each_check_fires_on_its_one_defect(
    refused: Callable[..., str], defect: ParametricDefect
) -> None:
    """Red test 8: the message starts with `P` and names the key's number, its
    name and the file's value (or that it is absent)."""
    refused(austria_keys(defect.changes), *defect.must_name)


def test_a_disagreeing_inverse_flattening_names_the_files_value(
    refused: Callable[..., str],
) -> None:
    """Row 6 for 2059 = 299.0: the file's value is named, not only EPSG:4312's
    own 299.1528128. Matched as a number (`299` or `299.0`, not followed by
    other digits), so the expected value cannot stand in for it."""
    (defect,) = (d for d in PARAMETRIC_DEFECTS if d.name == "inverse_flattening_disagrees")
    message = refused(austria_keys(defect.changes), *defect.must_name)
    assert re.search(r"(?<![\d.])299(\.0+)?(?![\d.])", message), message


def test_every_row_from_1_to_11_has_a_defect() -> None:
    """Red test 8 covers rows 1 to 11; rows 12 and 13 are red tests 5 and 7."""
    assert {d.row for d in PARAMETRIC_DEFECTS} == set(range(1, 12))


# ---------------------------------------------------------------------------
# Red test 9: the two existing user-defined cases now reach check 2
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "changes",
    [{PROJECTED_CS_TYPE: USER_DEFINED}, {PROJECTED_CS_TYPE: USER_DEFINED, GEOGRAPHIC_TYPE: 4326}],
    ids=["user_defined", "user_defined_beside_geographic"],
)
def test_existing_user_defined_cases_reach_the_method_check(
    geotiff: ModuleType, geotiff_error: type[Exception], changes: dict[int, int | None]
) -> None:
    """Red test 9: `test_refuses_missing_crs`'s two 3072 = 32767 cases (in
    `test_io_geotiff.py`, unchanged) keep passing because their tokens are in
    `P`. Here: on the baseline keys, with no 3074 and no 3075, it is check 2
    that fires."""
    with pytest.raises(geotiff_error) as info:
        geotiff.read_meta(micro_tiff(geokeys=with_keys(changes)))
    message = str(info.value)
    assert message.startswith(P), message
    assert "3075" in message and "ProjCoordTransGeoKey" in message, message


# ---------------------------------------------------------------------------
# Red test 10: the match is cached
# ---------------------------------------------------------------------------


def test_a_second_read_of_the_same_keys_does_not_identify_again(
    geotiff: ModuleType, fresh_cache: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Red test 10: a mosaic reads every tile's header. After `cache_clear()`
    (by `fresh_cache`), the first read asks PROJ to identify (calls
    `pyproj.CRS.list_authority`), and the second read of the same keys asks
    nothing more."""
    calls: list[int] = []
    original = pyproj.CRS.list_authority

    def counted(self: pyproj.CRS, *args: Any, **kwargs: Any) -> Any:
        calls.append(1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(pyproj.CRS, "list_authority", counted)
    assert geotiff.read_meta(austria_tiff()).epsg == EPSG_AUSTRIA_LAMBERT
    first = len(calls)
    assert first >= 1, "the first read did not identify, so the count measures nothing"
    assert geotiff.read_meta(austria_tiff()).epsg == EPSG_AUSTRIA_LAMBERT
    assert len(calls) == first


# ---------------------------------------------------------------------------
# Red test 12: Ola's file, header only, where it is present
# ---------------------------------------------------------------------------


@pytest.mark.skipif(AUSTRIA_FILE is None, reason="needs ../rasputin_data/austria_dgm10")
def test_the_real_austrian_file_reads_as_austria_lambert(geotiff: ModuleType) -> None:
    """Red test 12: 58 061 x 31 793 cells of 10 m, tie point (108 875, 586 555),
    pixel is area, so the first node is half a cell in."""
    assert AUSTRIA_FILE is not None
    with AUSTRIA_FILE.open("rb") as stream:
        meta = geotiff.read_meta(stream)
    assert meta.crs == "EPSG:31287"
    assert meta.epsg == EPSG_AUSTRIA_LAMBERT
    assert (meta.cols, meta.rows) == (58061, 31793)
    assert (meta.delta_x, meta.delta_y) == (10.0, 10.0)
    assert (meta.x_min, meta.y_max) == (108880.0, 586550.0)
    assert meta.pixel_is_area is True


# ---------------------------------------------------------------------------
# Red test 13: a parameter that is not one finite number (rows 6 and 11, `N`)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "defect", NOT_ONE_NUMBER_DEFECTS, ids=[d.name for d in NOT_ONE_NUMBER_DEFECTS]
)
def test_a_parameter_that_is_not_one_finite_number_is_refused(
    refused: Callable[..., str], defect: ParametricDefect
) -> None:
    """Red test 13: `refused` takes `GeoTiffError` only, so a `TypeError` or
    pyproj's `CRSError` escaping fails here. The message starts with `P` and
    names the key's number and name, `not one finite number`, and the file's
    value as Python writes it (`nan`, `inf`, or the tuple), matched as a whole
    token so that, say, the `inf` inside another word cannot stand in for it."""
    message = refused(austria_keys(defect.changes), *defect.must_name)
    (value,) = defect.changes.values()
    assert re.search(rf"(?<![\w.]){re.escape(repr(value))}(?![\w.])", message), message


# ---------------------------------------------------------------------------
# Red test 14: a file in grads is refused for its unit (row 8 before row 6)
# ---------------------------------------------------------------------------


def test_a_file_in_grads_is_refused_for_its_unit_not_its_meridian(
    refused: Callable[..., str],
) -> None:
    """Red test 14: 2054 = 9105 (grad) with 2061 = 2.5969213, Paris in grads.
    2061 is in 2054's unit, so it is the unit that is refused (row 8), and the
    meridian is never compared against EPSG:4312's in degrees (row 6)."""
    message = refused(GRADS_WITH_PARIS_MERIDIAN, "GeogAngularUnitsGeoKey", "2054", "9105")
    assert "disagrees" not in message.lower(), message
