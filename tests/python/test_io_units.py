"""Height units and feet coordinates: increment 43, section 6, tests 1 to 9 and 12.

`docs/increments/43-units-and-projections.md`. Tests 1 to 8 are the
invariant-critical suite (section 6): each factor of `UNITS`, the file's unit
over the caller's, the contradiction refusal, the NoData mask of `to_metres`,
the factor on `delta_*` as well as on the corner, and the metre case that
writes nothing.

Fixtures are 3 x 3 micro-GeoTIFFs written with tifffile (`geotiff_fixtures`).
Every expected height is the section's own figure: the cell times the unit's
metres in float64, cast once to float32 (section 4.3). A float32 multiply
gives other values in four cells of tests 1 and 2, so the comparison is bit
for bit, not approximate. Coordinates are float64 products of the file value
and the factor, compared exactly too.

Red in b1ec3a52: at the design's base 4099 other than 9001 was refused, 4096
was never read, a projected CRS in feet was refused, and the readers took no
`z_unit=`. The kill record of tests 1 to 8 is in the commit that wrote this
sentence.
"""

from __future__ import annotations

import io
import re
from collections.abc import Mapping
from typing import Any

import numpy as np
import pyproj
import pytest
import tifffile

from geotiff_fixtures import (
    GEOGRAPHIC_TYPE,
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    PIXEL_IS_AREA,
    PIXEL_IS_POINT,
    PROJ_COORD_TRANS,
    PROJ_FALSE_ORIGIN_EASTING,
    PROJ_FALSE_ORIGIN_LAT,
    PROJ_FALSE_ORIGIN_LONG,
    PROJ_FALSE_ORIGIN_NORTHING,
    PROJ_LINEAR_UNITS,
    PROJ_STD_PARALLEL_1,
    PROJ_STD_PARALLEL_2,
    PROJECTED_CS_TYPE,
    PROJECTION,
    USER_DEFINED,
    VERTICAL_UNITS,
    GeoKeyValue,
    micro_tiff,
    with_keys,
)
from tin_engine import mosaic
from tin_engine.io import cog, geotiff
from tin_engine.io.models import GeoTiffError, IndexWindow, TileFootprint

VERTICAL_CS = 4096
METRE, FOOT, US_FOOT, CLARKE_FOOT = 9001, 9002, 9003, 9005
#: Metres per unit: EPSG's definitions (section 5). Not pyproj's
#: 0.30480060960121924, which differs from 1200/3937 in the last bit.
FT, FT_US = 0.3048, 1200 / 3937

#: Test 1's cells; the last is the NoData sentinel.
HEIGHTS = np.array([[100, 200, 300], [400, 500, 600], [700, 800, -9999]], dtype=np.int16)
SENTINEL = "-9999"
#: Section 6, test 1: `HEIGHTS` in US survey feet, in metres.
IN_US_FEET = np.array(
    [
        [30.480060577392578, 60.960121154785156, 91.440185546875],
        [121.92024230957031, 152.40029907226562, 182.88037109375],
        [213.3604278564453, 243.84048461914062, -9999.0],
    ],
    dtype=np.float32,
)
#: Section 6, test 2: `HEIGHTS` in feet, in metres.
IN_FEET = np.array(
    [
        [30.479999542236328, 60.959999084472656, 91.44000244140625],
        [121.91999816894531, 152.39999389648438, 182.8800048828125],
        [213.36000061035156, 243.83999633789062, -9999.0],
    ],
    dtype=np.float32,
)
AS_READ = HEIGHTS.astype(np.float32)

#: Tests 6 and 8: EPSG:2229 (NAD83 / California zone 5 (ftUS)), area-registered,
#: tie (0, 0) at (6 500 000, 2 100 000) ftUS, 3 x 3 ftUS cells.
FEET_TIE = (0.0, 0.0, 0.0, 6_500_000.0, 2_100_000.0, 0.0)
FEET_SCALE = (3.0, 3.0, 0.0)
CALIFORNIA_5_FT, CALIFORNIA_5_M = 2229, 26945
FEET_KEYS: Mapping[int, GeoKeyValue] = {
    GT_MODEL_TYPE: 1,
    GT_RASTER_TYPE: PIXEL_IS_AREA,
    PROJECTED_CS_TYPE: CALIFORNIA_5_FT,
    VERTICAL_UNITS: US_FOOT,
}
#: Section 6, test 6: each `(file value) * 1200/3937`, the corner moved half a cell first.
FEET_X_MIN, FEET_Y_MAX, FEET_DELTA = 1981204.419608839, 640080.8229616459, 0.9144018288036575


def same(a: np.ndarray, b: np.ndarray) -> bool:
    """Bit for bit: dtype, shape and bytes."""
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def heights(changes: Mapping[int, GeoKeyValue | None], array: np.ndarray = HEIGHTS) -> io.BytesIO:
    """The usual EPSG:25833 baseline holding `array`, NoData -9999, with `changes`."""
    return micro_tiff(array, geokeys=with_keys(changes), nodata=SENTINEL)


def feet(
    changes: Mapping[int, GeoKeyValue | None] | None = None,
    tiepoint: tuple[float, ...] = FEET_TIE,
) -> io.BytesIO:
    """Test 6's file, with `changes` to its GeoKeys and another tie point if given."""
    keys = with_keys(changes, base=FEET_KEYS)
    return micro_tiff(HEIGHTS, tiepoint=tiepoint, scale=FEET_SCALE, geokeys=keys, nodata=SENTINEL)


def windowed(stream: io.BytesIO, **kw: Any) -> np.ndarray:
    """The whole page through `read_page` and `cog.decode_window`, as a mesh
    from the cache reads it."""
    meta, dtype, page = geotiff.read_page(stream, **kw)
    window = IndexWindow(row0=0, col0=0, rows=meta.rows, cols=meta.cols)
    blocks = cog.LocalTiffBlocks(page, stream, "units.tif")
    return np.asarray(cog.decode_window(blocks, meta, dtype, window).array)


def refusal(stream: io.BytesIO, *names: str, **kw: Any) -> str:
    """`decode_dem` refuses with `GeoTiffError` naming each of `names`."""
    with pytest.raises(GeoTiffError) as info:
        geotiff.decode_dem(stream, **kw)
    message = str(info.value)
    for name in names:
        assert name in message, f"{name!r} not in {message!r}"
    return message


# ---------------------------------------------------------------------------
# Tests 1 and 2: 4099 in feet or US survey feet
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("unit", "expected", "name"),
    [(US_FOOT, IN_US_FEET, "US survey foot"), (FOOT, IN_FEET, "foot")],
    ids=["us_survey_foot", "foot"],
)
class TestHeightsInFeet:
    """Tests 1 and 2: 4099 = 9003 or 9002 on a metre CRS; NoData untouched."""

    def test_decode_dem_gives_metres(self, unit: int, expected: np.ndarray, name: str) -> None:
        tile = geotiff.decode_dem(heights({VERTICAL_UNITS: unit}))
        assert same(tile.array, expected), tile.array.tolist()
        assert tile.meta.z_unit == name
        assert tile.meta.vertical_unit_assumed is False
        assert tile.meta.nodata == -9999.0

    def test_decode_window_gives_metres(self, unit: int, expected: np.ndarray, name: str) -> None:
        assert same(windowed(heights({VERTICAL_UNITS: unit})), expected)

    def test_the_header_names_the_file_unit(
        self, unit: int, expected: np.ndarray, name: str
    ) -> None:
        meta = geotiff.read_meta(heights({VERTICAL_UNITS: unit}))
        assert (meta.z_unit, meta.vertical_unit_assumed) == (name, False)


class TestAConvertedHeightEqualToNoData:
    """Test 2: -10000 ft is -3048.0 m in float64 and in float32, which is the
    sentinel -3048: the cell would turn into NoData, so the file is refused."""

    COLLIDES = np.array([[-10000, 1], [2, 3]], dtype=np.int16)

    def stream(self) -> io.BytesIO:
        keys = with_keys({VERTICAL_UNITS: FOOT})
        return micro_tiff(self.COLLIDES, geokeys=keys, nodata="-3048")

    def test_the_figure(self) -> None:
        assert -10000 * FT == -3048.0 and np.float32(-10000 * FT) == np.float32(-3048.0)

    def test_decode_dem_refuses(self) -> None:
        refusal(self.stream(), "-3048")

    def test_decode_window_refuses(self) -> None:
        with pytest.raises(GeoTiffError, match="-3048"):
            windowed(self.stream())


# ---------------------------------------------------------------------------
# Test 3: VerticalGeoKey (4096)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("changes", "expected", "name"),
    [
        ({VERTICAL_CS: 6360}, IN_US_FEET, "US survey foot"),  # NAVD88 height (ftUS)
        ({VERTICAL_CS: 8228}, IN_FEET, "foot"),  # NAVD88 height (ft)
        ({VERTICAL_CS: 5703}, AS_READ, "metre"),  # NAVD88 height
        ({VERTICAL_CS: 6360, VERTICAL_UNITS: US_FOOT}, IN_US_FEET, "US survey foot"),
        ({VERTICAL_CS: USER_DEFINED, VERTICAL_UNITS: FOOT}, IN_FEET, "foot"),
    ],
    ids=["6360_alone", "8228_alone", "5703_alone", "6360_and_4099_agree", "32767_and_4099"],
)
def test_the_vertical_crs_gives_the_unit(
    changes: Mapping[int, GeoKeyValue], expected: np.ndarray, name: str
) -> None:
    tile = geotiff.decode_dem(heights(changes))
    assert same(tile.array, expected), tile.array.tolist()
    assert (tile.meta.z_unit, tile.meta.vertical_unit_assumed) == (name, False)


@pytest.mark.parametrize(
    ("changes", "names"),
    [
        ({VERTICAL_CS: 6360, VERTICAL_UNITS: METRE}, ("4096", "6360", "4099", "9001")),
        ({VERTICAL_CS: 6358}, ("4096", "6358")),  # NAVD88 depth (ftUS): axis down
        ({VERTICAL_CS: 4326}, ("4096", "4326")),  # WGS 84: not vertical
        ({VERTICAL_CS: 1}, ("4096",)),  # no such EPSG code
        ({VERTICAL_UNITS: CLARKE_FOOT}, ("4099", "9005", "9001", "9002", "9003")),
    ],
    ids=["4096_and_4099_disagree", "a_depth", "not_vertical", "no_such_code", "unit_not_read"],
)
def test_the_vertical_keys_refused(
    changes: Mapping[int, GeoKeyValue], names: tuple[str, ...]
) -> None:
    """Each message names the key by number, and what it holds."""
    refusal(heights(changes), *names)
    with pytest.raises(GeoTiffError):  # the header alone refuses it too
        geotiff.read_meta(heights(changes))


# ---------------------------------------------------------------------------
# Test 4: the caller's unit, `z_unit=`
# ---------------------------------------------------------------------------


class TestTheCallersUnit:
    def test_a_silent_file_takes_the_callers_unit(self) -> None:
        tile = geotiff.decode_dem(heights({}), z_unit="US survey foot")
        assert same(tile.array, IN_US_FEET), tile.array.tolist()
        assert (tile.meta.z_unit, tile.meta.vertical_unit_assumed) == ("US survey foot", True)

    def test_through_decode_window(self) -> None:
        assert same(windowed(heights({}), z_unit="foot"), IN_FEET)

    @pytest.mark.parametrize("entry", ["read_meta", "read_header", "read_page", "decode_dem"])
    def test_every_entry_point_takes_it(self, entry: str) -> None:
        result = getattr(geotiff, entry)(heights({}), z_unit="foot")
        meta = (
            result if entry == "read_meta" else result.meta if entry == "decode_dem" else result[0]
        )
        assert (meta.z_unit, meta.vertical_unit_assumed) == ("foot", True)

    def test_the_file_wins_when_both_say_the_same(self) -> None:
        """Not assumed: the file says it (the file over the caller)."""
        tile = geotiff.decode_dem(heights({VERTICAL_UNITS: US_FOOT}), z_unit="US survey foot")
        assert same(tile.array, IN_US_FEET)
        assert (tile.meta.z_unit, tile.meta.vertical_unit_assumed) == ("US survey foot", False)

    @pytest.mark.parametrize(
        ("unit", "caller", "patterns"),
        [
            (METRE, "foot", (r"metre", r"f(oo|ee)t")),
            (US_FOOT, "foot", (r"survey", r"f(oo|ee)t.*f(oo|ee)t")),
        ],
        ids=["metre_file_foot_caller", "us_foot_file_foot_caller"],
    )
    def test_a_contradiction_is_refused_naming_both(
        self, unit: int, caller: str, patterns: tuple[str, ...]
    ) -> None:
        """Both units named, in either number ("foot" or "feet")."""
        message = refusal(heights({VERTICAL_UNITS: unit}), z_unit=caller)
        for pattern in patterns:
            assert re.search(pattern, message, re.I | re.S), f"{pattern!r} not in {message!r}"


# ---------------------------------------------------------------------------
# Test 5: a silent metre file is read as today, and nothing is written
# ---------------------------------------------------------------------------


class TestASilentMetreFile:
    def test_bytes_identical_to_the_files_own(self) -> None:
        stream = heights({})
        own = tifffile.imread(io.BytesIO(stream.getvalue())).astype(np.float32)
        tile = geotiff.decode_dem(stream)
        assert same(tile.array, own)
        assert (tile.meta.z_unit, tile.meta.vertical_unit_assumed) == ("metre", True)

    def test_bytes_identical_through_decode_window(self) -> None:
        assert same(windowed(heights({})), AS_READ)

    @pytest.mark.parametrize("changes", [{}, {VERTICAL_UNITS: METRE}], ids=["silent", "stated"])
    def test_to_metres_does_not_write_a_metre_array(
        self, changes: Mapping[int, GeoKeyValue]
    ) -> None:
        """A read-only array: any write, even of the same values, raises."""
        meta = geotiff.read_meta(heights(changes))
        array = AS_READ.copy()
        array.flags.writeable = False
        geotiff.to_metres(array, meta)
        assert same(array, AS_READ)

    def test_to_metres_writes_a_feet_array_in_place(self) -> None:
        meta = geotiff.read_meta(heights({VERTICAL_UNITS: US_FOOT}))
        array = AS_READ.copy()
        assert geotiff.to_metres(array, meta) is None
        assert same(array, IN_US_FEET)


# ---------------------------------------------------------------------------
# Tests 6 to 8: coordinates in feet, read as the CRS's metre twin
# ---------------------------------------------------------------------------


class TestFeetCoordinatesNamedTwin:
    """Test 6: EPSG:2229 is read as EPSG:26945, every coordinate times 1200/3937."""

    def test_the_grid_in_metres(self) -> None:
        meta = geotiff.read_meta(feet())
        assert (meta.crs, meta.epsg, meta.geographic) == (f"EPSG:{CALIFORNIA_5_M}", 26945, False)
        assert meta.file_crs == f"EPSG:{CALIFORNIA_5_FT}"
        assert meta.x_min == FEET_X_MIN
        assert meta.y_max == FEET_Y_MAX
        assert meta.delta_x == FEET_DELTA
        assert meta.delta_y == FEET_DELTA

    def test_the_figures(self) -> None:
        assert FEET_X_MIN == (6_500_000.0 + 1.5) * FT_US
        assert FEET_Y_MAX == (2_100_000.0 - 1.5) * FT_US
        assert FEET_DELTA == 3.0 * FT_US

    def test_heights_in_metres_too(self) -> None:
        tile = geotiff.decode_dem(feet())
        assert same(tile.array, IN_US_FEET)

    def test_a_metre_file_has_no_file_crs(self) -> None:
        assert geotiff.read_meta(heights({})).file_crs is None


FLORIDA_EAST_FT = 2881


class TestFeetCoordinatesUnnamedTwin:
    """Test 7: EPSG:2881's twin matches no EPSG code (its false easting is
    5e-10 relative from EPSG:2777's), so it is carried as its own text."""

    POINT = (700_000.0, 600_000.0)

    def meta(self) -> Any:
        keys = {PROJECTED_CS_TYPE: FLORIDA_EAST_FT, GT_RASTER_TYPE: PIXEL_IS_POINT}
        tie = (0.0, 0.0, 0.0, *self.POINT, 0.0)
        return geotiff.read_meta(feet(keys, tiepoint=tie))

    def test_no_epsg_code_and_the_twins_name(self) -> None:
        meta = self.meta()
        assert meta.epsg is None
        assert pyproj.CRS.from_user_input(meta.crs).name == (
            "NAD83(HARN) / Florida East (ftUS), in metres"
        )
        assert meta.file_crs == f"EPSG:{FLORIDA_EAST_FT}"
        assert (meta.x_min, meta.y_max) == (self.POINT[0] * FT_US, self.POINT[1] * FT_US)

    def test_a_point_and_its_metre_image_land_on_the_same_place(self) -> None:
        twin = pyproj.CRS.from_user_input(self.meta().crs)
        x, y = self.POINT
        by_feet = pyproj.Transformer.from_crs(FLORIDA_EAST_FT, 4326, always_xy=True)
        by_twin = pyproj.Transformer.from_crs(twin, 4326, always_xy=True)
        lon_f, lat_f = by_feet.transform(x, y)
        lon_m, lat_m = by_twin.transform(x * FT_US, y * FT_US)
        assert (lon_m - lon_f, lat_m - lat_f) == (0.0, 0.0)


class TestFeetCoordinatesSilentHeights:
    """Test 8: a feet DEM whose height unit is unstated is never guessed."""

    def test_refused_naming_the_option(self) -> None:
        refusal(feet({VERTICAL_UNITS: None}), "--z-unit")

    def test_the_header_refuses_too(self) -> None:
        with pytest.raises(GeoTiffError, match="--z-unit"):
            geotiff.read_meta(feet({VERTICAL_UNITS: None}))

    def test_read_with_the_callers_unit(self) -> None:
        tile = geotiff.decode_dem(feet({VERTICAL_UNITS: None}), z_unit="US survey foot")
        assert tile.meta.crs == f"EPSG:{CALIFORNIA_5_M}"
        assert (tile.meta.z_unit, tile.meta.vertical_unit_assumed) == ("US survey foot", True)
        assert same(tile.array, IN_US_FEET)


# ---------------------------------------------------------------------------
# Test 9: a feet CRS given by its parameters
# ---------------------------------------------------------------------------

#: EPSG:2229's parameters as GeoKeys: Lambert 2SP (3075 = 8) on NAD83 (4269),
#: 3076 = 9003, the false origin in US survey feet.
PARAMETRIC_FEET: Mapping[int, GeoKeyValue] = {
    GT_MODEL_TYPE: 1,
    GT_RASTER_TYPE: PIXEL_IS_AREA,
    GEOGRAPHIC_TYPE: 4269,
    PROJECTED_CS_TYPE: USER_DEFINED,
    PROJECTION: USER_DEFINED,
    PROJ_COORD_TRANS: 8,
    PROJ_LINEAR_UNITS: US_FOOT,
    PROJ_FALSE_ORIGIN_LAT: 33.5,
    PROJ_FALSE_ORIGIN_LONG: -118.0,
    PROJ_STD_PARALLEL_1: 35.4666666666667,
    PROJ_STD_PARALLEL_2: 34.0333333333333,
    PROJ_FALSE_ORIGIN_EASTING: 6561666.667,
    PROJ_FALSE_ORIGIN_NORTHING: 1640416.667,
    VERTICAL_UNITS: US_FOOT,
}


def test_parametric_feet_read_as_the_metre_code() -> None:
    stream = micro_tiff(
        HEIGHTS, tiepoint=FEET_TIE, scale=FEET_SCALE, geokeys=PARAMETRIC_FEET, nodata=SENTINEL
    )
    meta = geotiff.read_meta(stream)
    assert (meta.crs, meta.epsg) == (f"EPSG:{CALIFORNIA_5_M}", CALIFORNIA_5_M)
    assert (meta.x_min, meta.y_max) == (FEET_X_MIN, FEET_Y_MAX)
    assert (meta.delta_x, meta.delta_y) == (FEET_DELTA, FEET_DELTA)


# ---------------------------------------------------------------------------
# Test 12: the mosaic
# ---------------------------------------------------------------------------


def footprint(name: str, stream: io.BytesIO) -> TileFootprint:
    meta, dtype = geotiff.read_header(stream)
    return TileFootprint(name=name, meta=meta, dtype=dtype)


class TestMosaic:
    def test_tiles_in_two_height_units_are_refused_naming_both(self) -> None:
        """Two point-registered tiles sharing a column: aligned, both selected."""
        west = heights({VERTICAL_UNITS: METRE})
        east = micro_tiff(
            HEIGHTS,
            tiepoint=(0.0, 0.0, 0.0, 500_020.0, 6_600_000.0, 0.0),
            geokeys=with_keys({VERTICAL_UNITS: FOOT}),
            nodata=SENTINEL,
        )
        tiles = [footprint("west.tif", west), footprint("east.tif", east)]
        with pytest.raises(mosaic.MosaicError) as info:
            mosaic.plan_mosaic(tiles)
        message = str(info.value)
        assert "west.tif" in message and "east.tif" in message, message

    def test_one_feet_tile_is_returned_as_loaded(self) -> None:
        """As a metre tile is (I7): the plan's meta equals the tile's, `z_unit`
        and `file_crs` included, so no canvas is made."""
        tile = geotiff.decode_dem(feet())
        plan = mosaic.plan_mosaic([footprint("feet.tif", feet())])
        assert plan.meta == tile.meta
        result = mosaic.assemble(plan, lambda name: tile)
        assert result.tile is tile
