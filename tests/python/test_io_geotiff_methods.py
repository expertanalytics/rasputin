"""Projection methods given by parameters: increment 43, section 6, tests 10 and 11.

`docs/increments/43-units-and-projections.md`, section 4.5. One file per row
of the method table: the parameters of an EPSG CRS (pyproj 3.8.0 / PROJ
9.8.1), written to the GeoKeys that table names, read back as that CRS's
code. Each file with its longitude parameter moved 1 degree is refused as
"cannot be named". The shapes the table refuses are refused naming their keys.

The refusal contract is the earlier design's
(`docs/increments/geotiff-crs-by-parameters.md`, section 4): every message
starts with `P`, then names its keys by number.

HOW THIS FILE GOES RED. At the design's base only 3075 = 1 and 8 are read:
every other file is refused naming 3075, and the polar files would match no
EPSG code even with a method row, because their EPSG CRSs label their axes by
meridian (section 4.5, "Matching widened for polar axes"). Test 11 is green
at the base and must stay so.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any, NamedTuple

import pyproj
import pytest

from geotiff_fixtures import (
    EPSG_AUSTRIA_LAMBERT,
    EPSG_ETRS89,
    GEOGRAPHIC_TYPE,
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    PIXEL_IS_AREA,
    PROJ_COORD_TRANS,
    PROJ_FALSE_EASTING,
    PROJ_FALSE_NORTHING,
    PROJ_LINEAR_UNITS,
    PROJ_NAT_ORIGIN_LAT,
    PROJ_NAT_ORIGIN_LONG,
    PROJ_SCALE_AT_NAT_ORIGIN,
    PROJ_STD_PARALLEL_1,
    PROJ_STD_PARALLEL_2,
    PROJECTED_CS_TYPE,
    PROJECTION,
    USER_DEFINED,
    GeoKeyValue,
    austria_tiff,
    micro_tiff,
    transverse_mercator_on,
)
from tin_engine import crs as crs_rules
from tin_engine.io import geotiff
from tin_engine.io.models import GeoTiffError

P = "ProjectedCSTypeGeoKey (3072) = 32767 (a CRS given by parameters)"
METRE = 9001

# GeoKeys the earlier fixtures do not name (GeoTIFF 1.0 names).
PROJ_CENTER_LONG = 3088
PROJ_CENTER_LAT = 3089
PROJ_SCALE_AT_CENTER = 3093
PROJ_AZIMUTH_ANGLE = 3094
PROJ_STRAIGHT_VERT_POLE_LONG = 3095
PROJ_RECTIFIED_GRID_ANGLE = 3096

# ProjCoordTransGeoKey (3075) codes.
CT_HOTINE_A, CT_MERCATOR, CT_LAEA, CT_ALBERS, CT_POLAR = 3, 7, 10, 11, 15
CT_HOTINE_B = 9815
CT_AZIMUTHAL_EQUIDISTANT = 12

WGS84, NAD83, GDM2000, CH1903_PLUS = 4326, 4269, 4742, 4150


class Row(NamedTuple):
    """One file: the base geographic CRS (2048), 3075, the parameters by
    GeoKey, the longitude key moved for the refusal, and the code it reads as."""

    base: int
    method: int
    parameters: Mapping[int, float]
    longitude: int
    expected: int


#: Section 6, test 10. Parameters copied from each EPSG entry (pyproj 3.8.0).
ROWS: dict[str, Row] = {
    "polar_b_3031": Row(WGS84, CT_POLAR, {
        PROJ_NAT_ORIGIN_LAT: -71.0, PROJ_STRAIGHT_VERT_POLE_LONG: 0.0,
        PROJ_FALSE_EASTING: 0.0, PROJ_FALSE_NORTHING: 0.0,
    }, PROJ_STRAIGHT_VERT_POLE_LONG, 3031),
    "polar_a_5041": Row(WGS84, CT_POLAR, {
        PROJ_NAT_ORIGIN_LAT: 90.0, PROJ_STRAIGHT_VERT_POLE_LONG: 0.0,
        PROJ_SCALE_AT_NAT_ORIGIN: 0.994, PROJ_FALSE_EASTING: 2_000_000.0,
        PROJ_FALSE_NORTHING: 2_000_000.0,
    }, PROJ_STRAIGHT_VERT_POLE_LONG, 5041),
    # EPSG:32661 (UPS North (N,E)) has 5041's parameters; one CRS, two axis
    # orders, and the lowest code is taken.
    "polar_a_32661_as_5041": Row(WGS84, CT_POLAR, {
        PROJ_NAT_ORIGIN_LAT: 90.0, PROJ_STRAIGHT_VERT_POLE_LONG: 0.0,
        PROJ_SCALE_AT_NAT_ORIGIN: 0.994, PROJ_FALSE_EASTING: 2_000_000.0,
        PROJ_FALSE_NORTHING: 2_000_000.0,
    }, PROJ_STRAIGHT_VERT_POLE_LONG, 5041),
    "laea_3035": Row(EPSG_ETRS89, CT_LAEA, {
        PROJ_CENTER_LAT: 52.0, PROJ_CENTER_LONG: 10.0,
        PROJ_FALSE_EASTING: 4_321_000.0, PROJ_FALSE_NORTHING: 3_210_000.0,
    }, PROJ_CENTER_LONG, 3035),
    "laea_6931": Row(WGS84, CT_LAEA, {
        PROJ_CENTER_LAT: 90.0, PROJ_CENTER_LONG: 0.0,
        PROJ_FALSE_EASTING: 0.0, PROJ_FALSE_NORTHING: 0.0,
    }, PROJ_CENTER_LONG, 6931),
    "mercator_a_3395": Row(WGS84, CT_MERCATOR, {
        PROJ_NAT_ORIGIN_LAT: 0.0, PROJ_NAT_ORIGIN_LONG: 0.0, PROJ_SCALE_AT_NAT_ORIGIN: 1.0,
        PROJ_FALSE_EASTING: 0.0, PROJ_FALSE_NORTHING: 0.0,
    }, PROJ_NAT_ORIGIN_LONG, 3395),
    "mercator_b_3994": Row(WGS84, CT_MERCATOR, {
        PROJ_STD_PARALLEL_1: -41.0, PROJ_NAT_ORIGIN_LONG: 100.0,
        PROJ_FALSE_EASTING: 0.0, PROJ_FALSE_NORTHING: 0.0,
    }, PROJ_NAT_ORIGIN_LONG, 3994),
    "albers_3338": Row(NAD83, CT_ALBERS, {
        PROJ_NAT_ORIGIN_LAT: 50.0, PROJ_NAT_ORIGIN_LONG: -154.0,
        PROJ_STD_PARALLEL_1: 55.0, PROJ_STD_PARALLEL_2: 65.0,
        PROJ_FALSE_EASTING: 0.0, PROJ_FALSE_NORTHING: 0.0,
    }, PROJ_NAT_ORIGIN_LONG, 3338),
    "hotine_a_3375": Row(GDM2000, CT_HOTINE_A, {
        PROJ_CENTER_LAT: 4.0, PROJ_CENTER_LONG: 102.25,
        PROJ_AZIMUTH_ANGLE: 323.025796466667, PROJ_RECTIFIED_GRID_ANGLE: 323.130102361111,
        PROJ_SCALE_AT_CENTER: 0.99984, PROJ_FALSE_EASTING: 804_671.0, PROJ_FALSE_NORTHING: 0.0,
    }, PROJ_CENTER_LONG, 3375),
    "hotine_b_2056": Row(CH1903_PLUS, CT_HOTINE_B, {
        PROJ_CENTER_LAT: 46.9524055555556, PROJ_CENTER_LONG: 7.43958333333333,
        PROJ_AZIMUTH_ANGLE: 90.0, PROJ_RECTIFIED_GRID_ANGLE: 90.0, PROJ_SCALE_AT_CENTER: 1.0,
        PROJ_FALSE_EASTING: 2_600_000.0, PROJ_FALSE_NORTHING: 1_200_000.0,
    }, PROJ_CENTER_LONG, 2056),
}  # fmt: skip


def parametric(
    base: int, method: int, parameters: Mapping[int, GeoKeyValue | None]
) -> dict[int, GeoKeyValue]:
    """A metre CRS given by parameters on `base`; a `None` parameter is left out."""
    keys: dict[int, GeoKeyValue] = {
        GT_MODEL_TYPE: 1,
        GT_RASTER_TYPE: PIXEL_IS_AREA,
        GEOGRAPHIC_TYPE: base,
        PROJECTED_CS_TYPE: USER_DEFINED,
        PROJECTION: USER_DEFINED,
        PROJ_COORD_TRANS: method,
        PROJ_LINEAR_UNITS: METRE,
    }
    keys |= {k: v for k, v in parameters.items() if v is not None}
    return keys


def read(keys: Mapping[int, GeoKeyValue]) -> Any:
    return geotiff.read_meta(micro_tiff(geokeys=keys))


def refused(keys: Mapping[int, GeoKeyValue], *names: str) -> str:
    """Refused with a message starting with `P` and naming each of `names`
    (case-insensitive substring)."""
    with pytest.raises(GeoTiffError) as info:
        read(keys)
    message = str(info.value)
    assert message.startswith(P), f"does not start with P: {message!r}"
    for name in names:
        assert name.lower() in message.lower(), f"{name!r} not in {message!r}"
    return message


# ---------------------------------------------------------------------------
# Test 10: one file per row
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("row", ROWS.values(), ids=ROWS.keys())
def test_each_row_reads_as_its_code(row: Row) -> None:
    meta = read(parametric(row.base, row.method, row.parameters))
    assert (meta.epsg, meta.crs, meta.geographic) == (row.expected, f"EPSG:{row.expected}", False)


@pytest.mark.parametrize("row", ROWS.values(), ids=ROWS.keys())
def test_each_row_moved_one_degree_cannot_be_named(row: Row) -> None:
    moved = {**row.parameters, row.longitude: row.parameters[row.longitude] + 1.0}
    refused(parametric(row.base, row.method, moved), "cannot be named")


@pytest.mark.parametrize("row", ROWS.values(), ids=ROWS.keys())
def test_each_row_refuses_its_longitude_key_absent(row: Row) -> None:
    """Section 4.5: "any listed key absent" is refused, naming the key."""
    absent = {**row.parameters, row.longitude: None}
    refused(parametric(row.base, row.method, absent), str(row.longitude), "absent")


MERCATOR_B = ROWS["mercator_b_3994"].parameters
POLAR_B = ROWS["polar_b_3031"].parameters


@pytest.mark.parametrize(
    ("method", "parameters", "names"),
    [
        (
            CT_MERCATOR,
            {**MERCATOR_B, PROJ_SCALE_AT_NAT_ORIGIN: 1.0},
            ("3078", "3092"),
        ),
        (
            CT_MERCATOR,
            {**MERCATOR_B, PROJ_STD_PARALLEL_1: None},
            ("3078", "3092"),
        ),
        (
            CT_POLAR,
            {**POLAR_B, PROJ_SCALE_AT_NAT_ORIGIN: 0.97},
            ("3092", "0.97"),
        ),
        (
            CT_MERCATOR,
            {**MERCATOR_B, PROJ_NAT_ORIGIN_LAT: 10.0},
            ("3081", "10"),
        ),
    ],
    ids=[
        "mercator_with_3078_and_3092",
        "mercator_with_neither",
        "polar_off_the_pole_with_a_scale",
        "mercator_b_with_a_latitude_of_origin",
    ],
)
def test_a_shape_the_table_refuses(
    method: int, parameters: Mapping[int, float | None], names: tuple[str, ...]
) -> None:
    refused(parametric(WGS84, method, parameters), *names)


def test_an_unread_method_lists_the_readable_ones() -> None:
    """3075 = 12 (azimuthal equidistant): the message names 3075 and its value,
    and lists the methods that are read."""
    message = refused(parametric(WGS84, CT_AZIMUTHAL_EQUIDISTANT, {}), "3075")
    assert re.search(r"(?<!\d)12(?!\d)", message), message
    flat = " ".join(message.lower().split())
    for method in (
        "transverse mercator",
        "lambert conic conformal",
        "oblique mercator",
        "lambert azimuthal equal area",
        "albers",
        "polar stereographic",
    ):
        assert method in flat, f"{method!r} not listed in {message!r}"
    assert flat.count("mercator") >= 3, "Mercator itself is not listed"


# ---------------------------------------------------------------------------
# Test 11: `epsg_matches` unchanged off the poles
# ---------------------------------------------------------------------------


def built_as_the_reader_builds(code: int) -> pyproj.CRS:
    """EPSG:`code`'s method and parameters on its base, axes east and north in
    metres: the shape of CRS the reader matches."""
    source = pyproj.CRS.from_epsg(code)
    conversion = source.to_json_dict()["conversion"]
    assert source.geodetic_crs is not None
    base = pyproj.CRS.from_epsg(source.geodetic_crs.to_epsg())
    axes = [
        {"name": n, "abbreviation": n[0], "direction": d, "unit": "metre"}
        for n, d in (("Easting", "east"), ("Northing", "north"))
    ]
    return pyproj.CRS.from_json_dict({
        "type": "ProjectedCRS", "name": "unknown", "base_crs": base.to_json_dict(),
        "conversion": {"name": "unknown", "method": conversion["method"],
                       "parameters": conversion["parameters"]},
        "coordinate_system": {"subtype": "Cartesian", "axis": axes},
    })  # fmt: skip


class TestOffThePolesUnchanged:
    def test_the_austrian_keys_still_read_as_31287(self) -> None:
        assert geotiff.read_meta(austria_tiff()).epsg == EPSG_AUSTRIA_LAMBERT

    def test_utm_35n_on_etrs89_still_reads_as_25835(self) -> None:
        keys = transverse_mercator_on(EPSG_ETRS89, 27.0)
        assert geotiff.read_meta(austria_tiff(keys)).epsg == 25835

    @pytest.mark.parametrize("code", [25835, EPSG_AUSTRIA_LAMBERT])
    def test_the_match_is_that_code_alone(self, code: int) -> None:
        assert crs_rules.epsg_matches(built_as_the_reader_builds(code)) == (code,)
