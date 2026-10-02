"""G1, the reader: a geographic 2D GeoTIFF in degrees (increment 15c-2, D6).

`docs/increments/15c-geographic-dem.md`, D6 and "Tests for @tester" G1.
`io/geotiff.py` accepts a geographic 2D CRS from GeographicTypeGeoKey (2048),
with degree units (2054 = 9102, or absent), on every entry point: `decode_dem`,
`read_meta`, `read_header` and `read_page` (23a-2's `geographic=` flag was the
fetch path's alone; 15c-2 lifts the mesh path's refusal, 23a-1 W9). The meta
says so: `RasterMeta.geographic` true and `crs` `EPSG:n`. Geocentric, 3D,
user-defined and non-degree geographic CRSs stay refused, each naming the key.

`RasterMeta` gains `crs: str` (filled from `epsg` as `EPSG:n` when not given)
and `geographic: bool = False`; `epsg` becomes `int | None`, so a target CRS
without an EPSG code (D8's PROJ string) has a meta.

HOW THIS FILE GOES RED. Before 15c-2 the reader refuses 2048 on the mesh
path, and `RasterMeta` ignores the two new fields (Pydantic's default
`extra="ignore"`), so the accepting tests fail on the refusal or on a
missing attribute. The refusal tests are guards: the 3D (4979) and
non-degree cases are new behaviour for the fetch path, but the mesh path
refuses them today for the wrong reason, so each also pins the key it names.
"""

from __future__ import annotations

import io
from collections.abc import Callable
from typing import Any, ClassVar

import numpy as np
import pytest
import tifffile

from geographic_fixtures import (
    ANADEM_STEP,
    DEGREE,
    GEOG_CITATION,
    GEOG_INV_FLATTENING,
    GEOG_SEMI_MAJOR_AXIS,
    GRAD,
    GRS80_A,
    GRS80_RF,
    LAT0,
    LON0,
    NODATA,
    RADIAN,
    geographic_keys,
    projected_tiff,
    rough,
    tiff,
)
from tin_engine.io import geotiff
from tin_engine.io.models import GeoTiffError, RasterMeta

Reader = Callable[[io.BytesIO], RasterMeta]

#: Every entry point, each reduced to the meta it returns.
READERS: dict[str, Reader] = {
    "decode_dem": lambda s: geotiff.decode_dem(s).meta,
    "read_meta": lambda s: geotiff.read_meta(s),
    "read_header": lambda s: geotiff.read_header(s)[0],
    "read_page": lambda s: geotiff.read_page(s)[0],
}


def anadem_like(epsg: int = 4326, **changes: int | None) -> io.BytesIO:
    """ANADEM's header numbers (15-dem-mosaic.md B2): float32, Deflate, tiled,
    area-registered, NoData -9999 in 42113, the 0.000269...° step."""
    return tiff(
        rough(32, 48),
        x0=LON0,
        y0=LAT0,
        step_x=ANADEM_STEP,
        shorts=geographic_keys(epsg, area=True, **changes),
        compression="deflate",
        tile=(16, 16),
    )


def cog_key_set() -> io.BytesIO:
    """The OpenTopography COG's key set ("Measured by @architect"): 1024 = 2,
    1025 = 1, 2048 = 4674 with its citation (2049), 2054 = 9102, and the
    ellipsoid's semi-major axis (2057) and inverse flattening (2059)."""
    return tiff(
        rough(32, 48),
        x0=-82.51654705336689,
        y0=14.079475111062084,
        step_x=ANADEM_STEP,
        shorts=geographic_keys(4674, area=True),
        doubles={GEOG_SEMI_MAJOR_AXIS: GRS80_A, GEOG_INV_FLATTENING: GRS80_RF},
        texts={GEOG_CITATION: "SIRGAS 2000"},
        compression="deflate",
        tile=(16, 16),
    )


def geokeys_of(stream: io.BytesIO) -> dict[str, Any]:
    stream.seek(0)
    with tifffile.TiffFile(stream) as tif:
        keys: dict[str, Any] = tif.geotiff_metadata or {}
    stream.seek(0)
    return keys


class TestTheFixturesCarryTheirKeys:
    """tifffile alone, no production code: each builder writes what it says."""

    def test_anadem_like(self) -> None:
        keys = geokeys_of(anadem_like())
        assert int(keys["GeographicTypeGeoKey"]) == 4326
        assert int(keys["GeogAngularUnitsGeoKey"]) == DEGREE
        assert int(keys["GTRasterTypeGeoKey"]) == 1
        assert "ProjectedCSTypeGeoKey" not in keys

    def test_the_cog_key_set(self) -> None:
        keys = geokeys_of(cog_key_set())
        assert int(keys["GeographicTypeGeoKey"]) == 4674
        assert int(keys["GeogAngularUnitsGeoKey"]) == DEGREE
        assert str(keys["GeogCitationGeoKey"]) == "SIRGAS 2000"
        assert float(keys["GeogSemiMajorAxisGeoKey"]) == GRS80_A
        assert float(keys["GeogInvFlatteningGeoKey"]) == GRS80_RF


@pytest.mark.parametrize("reader", READERS.values(), ids=READERS.keys())
class TestAccepted:
    def test_anadems_header_numbers_read_as_geographic_4326(self, reader: Reader) -> None:
        meta = reader(anadem_like())
        assert meta.geographic is True  # type: ignore[attr-defined]
        assert meta.crs == "EPSG:4326"  # type: ignore[attr-defined]
        assert meta.epsg == 4326
        assert (meta.delta_x, meta.delta_y) == (ANADEM_STEP, ANADEM_STEP)
        # Area-registered: the first node is half a cell in from the tie point.
        assert meta.x_min == LON0 + ANADEM_STEP / 2
        assert meta.y_max == LAT0 - ANADEM_STEP / 2
        assert (meta.rows, meta.cols, meta.nodata) == (32, 48, NODATA)

    def test_the_cog_key_set_reads_as_4674(self, reader: Reader) -> None:
        meta = reader(cog_key_set())
        assert meta.geographic is True  # type: ignore[attr-defined]
        assert meta.crs == "EPSG:4674"  # type: ignore[attr-defined]
        assert meta.epsg == 4674

    def test_an_absent_2054_is_degrees(self, reader: Reader) -> None:
        meta = reader(anadem_like(GEOG_ANGULAR_UNITS=None))
        assert meta.crs == "EPSG:4326"  # type: ignore[attr-defined]

    def test_a_projected_file_is_not_geographic(self, reader: Reader) -> None:
        stream = projected_tiff(rough(4, 5), x0=500_000.0, y0=6_600_000.0, step=10.0, epsg=25833)
        meta = reader(stream)
        assert meta.geographic is False  # type: ignore[attr-defined]
        assert meta.crs == "EPSG:25833"  # type: ignore[attr-defined]


def test_the_pixels_of_a_geographic_file_decode_unchanged() -> None:
    tile = geotiff.decode_dem(anadem_like())
    np.testing.assert_array_equal(tile.array, rough(32, 48))


@pytest.mark.parametrize(
    ("epsg", "changes", "names"),
    [
        (4979, {}, ("2048", "GeographicTypeGeoKey", "4979")),  # WGS 84 3D
        (32767, {}, ("2048", "GeographicTypeGeoKey", "32767")),  # user-defined
        (4978, {}, ("2048", "GeographicTypeGeoKey", "4978", "Geocentric")),
        (4326, {"GEOG_ANGULAR_UNITS": RADIAN}, ("2054", "GeogAngularUnitsGeoKey", "9101")),
        (4326, {"GEOG_ANGULAR_UNITS": GRAD}, ("2054", "GeogAngularUnitsGeoKey", "9105")),
        (4807, {}, ("2048", "GeographicTypeGeoKey", "4807")),  # NTF (Paris): axes in grads
    ],
    ids=["geographic_3d", "user_defined", "geocentric", "radians", "grads", "crs_axes_in_grads"],
)
@pytest.mark.parametrize("reader", READERS.values(), ids=READERS.keys())
def test_refused_naming_the_key(
    reader: Reader, epsg: int, changes: dict[str, int], names: tuple[str, ...]
) -> None:
    with pytest.raises(GeoTiffError) as caught:
        reader(anadem_like(epsg, **changes))
    message = str(caught.value)
    for name in names:
        assert name.lower() in message.lower(), f"{name!r} not in {message!r}"


class TestRasterMeta:
    BASE: ClassVar[dict[str, Any]] = dict(
        x_min=0.0,
        y_max=0.0,
        delta_x=30.0,
        delta_y=30.0,
        cols=2,
        rows=2,
        nodata=None,
        nodata_source="absent",
        pixel_is_area=False,
        vertical_unit_assumed=True,
    )

    def test_crs_is_filled_from_epsg(self) -> None:
        meta = RasterMeta(**self.BASE, epsg=31983)
        assert meta.crs == "EPSG:31983"  # type: ignore[attr-defined]
        assert meta.geographic is False  # type: ignore[attr-defined]

    def test_a_crs_without_an_epsg_code(self) -> None:
        proj = "+proj=tmerc +lat_0=0 +lon_0=-44.1 +k=0.999972 +x_0=0 +y_0=0 +ellps=GRS80"
        meta = RasterMeta(**self.BASE, epsg=None, crs=proj)
        assert meta.epsg is None
        assert meta.crs == proj  # type: ignore[attr-defined]
