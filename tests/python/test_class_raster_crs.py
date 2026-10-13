"""A class raster's own CRS and no-data: cell areas, fractions and decoding.

`docs/increments/44-catalogue-as-data.md`, sections 4.5 to 4.7 and section 7's
T8 (a projected class raster, in metres and in US survey feet), T9 (the
ellipsoid from the raster's CRS) and T10 (no-data in a class raster). The
fixtures are the design's probes P2, P3 and P9 (`docs/increments/44-probes/`),
written by `geotiff_fixtures.micro_tiff`.

The oracles: for areas, `pyproj.Geod` built here from the ellipsoid's own
numbers (or `ellps="WGS84"`, which is what every existing mesh used) and the
grid's `delta_x * delta_y` for a projected raster; for fractions, counting
cells by hand (a triangle whose legs lie on cell edges covers whole cells and
half cells, so its class shares are ratios of small integers); for decoding,
the written grid with the no-data cells and each empty tile zeroed.

Numeric bounds: the projected fractions are compared exactly in float32: the
EPSG:3035 triangle's corners are on cell corners (scale 1e6 m, exact in
binary), and the feet raster's float64 ratio is 0.7499999999995 (P9), 5e-13
from 0.75, far below float32's half-ulp at 0.75 (3e-8). The EPSG:4230 cell
area has relative tolerance 1e-9 against the probe's six printed decimals of
a 3.8e5 m2 value (P4); the WGS 84 rows are compared bit for bit, since a
changed bit changes a MapBiomas mesh's bytes.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The refusal of a stored 0 beside a no-data of 255 is a `GeoTiffError`
  whose message holds the block source's name, `255` and the word `0` (the
  design's sentence is not pinned word for word).
- `TiffDemRepository.load_class_window` applies the file's `GDAL_NODATA`
  exactly as `decode_class_window` does through `LocalTiffBlocks`.
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from pyproj import CRS, Geod

from class_fixtures import AREA_KEYS, LAT0, LON0, ROWS, STEP, class_tiff, expected
from cog_fixtures import with_zero_byte_count
from geotiff_fixtures import (
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    PIXEL_IS_AREA,
    PROJECTED_CS_TYPE,
    micro_tiff,
)
from tin_engine.fetch.plan import parse_prefix
from tin_engine.io.cog import LocalTiffBlocks, decode_class_window
from tin_engine.io.geotiff import read_page
from tin_engine.io.models import ClassWindow, GeoTiffError, IndexWindow
from tin_engine.io.repository import TiffDemRepository
from tin_engine.land_cover_fractions import assign, row_areas
from tin_engine.sources import RemoteSource

#: P3: 8 x 8 cells of 100 m in EPSG:3035, upper-left cell corner (X0, Y0).
X0, Y0 = 4_000_000.0, 3_000_000.0
#: P9: 8 x 8 cells of 100 US survey feet in EPSG:2229, tie at (FT_X0, FT_Y0).
FT_X0, FT_Y0 = 6_500_000.0, 2_100_000.0
#: P9: 100 US survey feet in metres, the metre twin's spacing.
FT_DELTA = 30.480060960121918
#: P9: FT_DELTA squared, the twin's cell area in m2.
FT_CELL = 929.0341161327483
#: The MapBiomas Collection 11 header's spacing and top, and its row count.
MB_STEP, MB_TOP, MB_ROWS = 0.00026949458523585647, 5.422904791408522, 146_501
VERTICAL_UNITS, US_SURVEY_FOOT = 4099, 9003


def halves() -> np.ndarray:
    """8 x 8 cells: columns 0-3 class 10, 4-7 class 20."""
    grid = np.full((8, 8), 10, np.uint8)
    grid[:, 4:] = 20
    return grid


def projected_classes(epsg: int, x0: float, y0: float) -> bytes:
    keys = {GT_MODEL_TYPE: 1, GT_RASTER_TYPE: PIXEL_IS_AREA, PROJECTED_CS_TYPE: epsg}
    return micro_tiff(
        halves(), tiepoint=(0, 0, 0, x0, y0, 0), scale=(100.0, 100.0, 0), geokeys=keys,
        tile=(16, 16),
    ).getvalue()  # fmt: skip


def feet_classes(*, vertical: int | None = None) -> bytes:
    keys = {GT_MODEL_TYPE: 1, GT_RASTER_TYPE: PIXEL_IS_AREA, PROJECTED_CS_TYPE: 2229}
    if vertical is not None:
        keys[VERTICAL_UNITS] = vertical
    return micro_tiff(
        halves(), tiepoint=(0, 0, 0, FT_X0, FT_Y0, 0), scale=(100.0, 100.0, 0), geokeys=keys,
        tile=(16, 16),
    ).getvalue()  # fmt: skip


def window(x_min: float, y_max: float, delta: float, rows: int, cols: int = 1) -> ClassWindow:
    return ClassWindow(
        x_min=x_min, y_max=y_max, delta_x=delta, delta_y=delta,
        array=np.zeros((rows, cols), dtype=np.uint8),
    )  # fmt: skip


def geod_rows(geod: Geod, w: ClassWindow, rows: range) -> np.ndarray:
    """The oracle: `geod`'s area of one cell at longitude 0 per row, the row's
    top `y_max - r * delta_y` and bottom one `delta_y` below it."""
    lons = [0.0, w.delta_x, w.delta_x, 0.0]
    out = []
    for r in rows:
        top = w.y_max - r * w.delta_y
        area, _ = geod.polygon_area_perimeter(lons, [top, top, top - w.delta_y, top - w.delta_y])
        out.append(abs(float(area)))
    return np.array(out)


# ------------------------------------------------------------ T8, metres


class TestT8ProjectedInMetres:
    @pytest.fixture
    def path(self, tmp_path: Path) -> Path:
        target = tmp_path / "laea_classes.tif"
        target.write_bytes(projected_classes(3035, X0, Y0))
        return target

    def test_every_row_is_the_grid_cells_area(self) -> None:
        got = row_areas(window(X0, Y0, 100.0, 8, 8), "EPSG:3035")
        assert got.dtype == np.float64
        assert_array_equal(got, np.full(8, 10_000.0))

    def test_the_triangle_takes_three_quarters_of_class_10(self, path: Path) -> None:
        """Legs of 8 cells on the raster's top and left edges: class 10
        (columns 0-3) covers 24 of the 32 cells' area, class 20 the other 8."""
        vertices = np.array([[X0, Y0], [X0 + 800.0, Y0], [X0, Y0 - 800.0]])
        got = assign(
            vertices, np.array([[0, 2, 1]]), mesh_crs="EPSG:3035",
            repository=TiffDemRepository([path]), name=path.name,
        )  # fmt: skip
        assert got.table.classes.tolist() == [10, 20]
        assert got.table.fractions.tolist() == [0.75, 0.25]
        assert got.codes.tolist() == [10]


# ------------------------------------------------------------ T8, US survey feet


class TestT8ProjectedInFeet:
    @pytest.fixture
    def path(self, tmp_path: Path) -> Path:
        target = tmp_path / "feet_classes.tif"
        target.write_bytes(feet_classes())
        return target

    def test_without_a_unit_the_reader_asks_for_one(self, path: Path) -> None:
        """Increment 43's height-unit check, unchanged: why both land-cover
        routes pass `z_unit="metre"`."""
        with pytest.raises(GeoTiffError, match="z-unit"):
            TiffDemRepository([path]).footprints()

    def test_with_metre_the_meta_is_the_metre_twins(self, path: Path) -> None:
        (footprint,) = TiffDemRepository([path], z_unit="metre").footprints()
        meta = footprint.meta
        assert (meta.crs, meta.file_crs) == ("EPSG:26945", "EPSG:2229")
        assert meta.delta_x == meta.delta_y == FT_DELTA

    def test_every_row_is_the_twins_cell_area(self, path: Path) -> None:
        (footprint,) = TiffDemRepository([path], z_unit="metre").footprints()
        meta = footprint.meta
        corner = window(meta.x_min - meta.delta_x / 2, meta.y_max + meta.delta_y / 2, FT_DELTA, 8)
        got = row_areas(corner, meta.crs)
        assert_array_equal(got, np.full(8, meta.delta_x * meta.delta_y))
        assert got[0] == pytest.approx(FT_CELL, rel=1e-15)

    def test_the_triangle_takes_three_quarters_of_class_10(self, path: Path) -> None:
        """Index-space corners (col, row) (2, 0), (6, 0), (2, 4): class 10
        (columns 2-3) covers 6 of the 8 cells' area."""
        repository = TiffDemRepository([path], z_unit="metre")
        (footprint,) = repository.footprints()
        meta = footprint.meta
        cx, cy = meta.x_min - meta.delta_x / 2, meta.y_max + meta.delta_y / 2
        vertices = np.array(
            [
                [cx + 2 * FT_DELTA, cy],
                [cx + 6 * FT_DELTA, cy],
                [cx + 2 * FT_DELTA, cy - 4 * FT_DELTA],
            ]
        )
        got = assign(
            vertices, np.array([[0, 2, 1]]), mesh_crs=meta.crs, repository=repository,
            name=path.name,
        )  # fmt: skip
        assert got.table.classes.tolist() == [10, 20]
        assert got.table.fractions.tolist() == [0.75, 0.25]
        assert got.codes.tolist() == [10]

    def test_the_fetcher_reads_a_classes_entry_in_metres(self) -> None:
        entry = RemoteSource(
            id="feet-classes", kind="one-cog", url="https://example.invalid/f.tif",
            values="classes", crs="EPSG:2229", credit="c", licence_note="l",
        )  # fmt: skip
        parsed = parse_prefix(feet_classes(), nodata=None, z_unit=entry.read_unit)
        assert parsed is not None
        meta, dtype, _ = parsed
        # parse_prefix returns the dtype the readers decode to (geotiff.PROMOTION):
        # a uint8 file decodes to float32, as on master.
        assert (meta.crs, dtype) == ("EPSG:26945", np.dtype(np.float32))


# ------------------------------------------------------------ T9


class TestT9TheEllipsoid:
    def test_wgs84_rows_are_todays_bit_for_bit(self) -> None:
        w = window(LON0, LAT0, STEP, ROWS, 48)
        assert_array_equal(
            row_areas(w, "EPSG:4326"), geod_rows(Geod(ellps="WGS84"), w, range(ROWS))
        )

    def test_mapbiomas_rows_are_todays_bit_for_bit(self) -> None:
        """P7: a `Geod` from `CRS.get_geod()` differs in 8,862 of these 20,929
        rows; built from the ellipsoid's `a` and `rf`, in none."""
        w = window(-74.0, MB_TOP, MB_STEP, MB_ROWS)
        sample = range(0, MB_ROWS, 7)
        got = row_areas(w, "EPSG:4326")
        assert len(got) == MB_ROWS
        assert_array_equal(got[::7], geod_rows(Geod(ellps="WGS84"), w, sample))

    def test_international_1924_is_its_own_ellipsoid(self) -> None:
        """P4: the top row of a 1/128-degree grid at 60 N in EPSG:4230."""
        got = row_areas(window(0.0, 60.0, STEP, 11), "EPSG:4230")
        assert got[0] == pytest.approx(379_523.462526, rel=1e-9)
        assert got[0] != pytest.approx(379_488.217823, rel=1e-9), "that is WGS 84's"
        assert got[10] == pytest.approx(380_413.219975, rel=1e-9)

    def test_a_sphere_has_flattening_0(self) -> None:
        sphere = "+proj=longlat +R=6371000 +no_defs"
        assert CRS(sphere).ellipsoid.inverse_flattening == 0.0
        w = window(0.0, 60.0, STEP, 4)
        assert_array_equal(row_areas(w, sphere), geod_rows(Geod(a=6_371_000.0, f=0.0), w, range(4)))


# ------------------------------------------------------------ T10


def nodata_grid() -> np.ndarray:
    """P2: 32 x 32, classes 7 (columns 0-15) and 9 (16-31), the top-left
    4 x 4 cells 255."""
    grid = np.full((32, 32), 7, np.uint8)
    grid[:, 16:] = 9
    grid[0:4, 0:4] = 255
    return grid


def nodata_tiff(grid: np.ndarray, *, empty: int | None = 3) -> bytes:
    data = micro_tiff(
        grid, tiepoint=(0, 0, 0, LON0, LAT0, 0), scale=(STEP, STEP, 0), geokeys=AREA_KEYS,
        tile=(16, 16), nodata="255",
    ).getvalue()  # fmt: skip
    return data if empty is None else with_zero_byte_count(data, empty)


def decoded(data: bytes, name: str, rows: int = 32) -> np.ndarray:
    stream = io.BytesIO(data)
    meta, _, page = read_page(stream)
    w = IndexWindow(row0=0, col0=0, rows=rows, cols=32)
    return decode_class_window(LocalTiffBlocks(page, stream, name), meta, w).array


def nodata_oracle() -> np.ndarray:
    """The written grid, its 255 cells and the empty tile 3 (rows and
    columns 16-31) as 0."""
    want = nodata_grid()
    want[want == 255] = 0
    want[16:, 16:] = 0
    return want


class TestT10NoData:
    def test_the_tagged_no_data_becomes_0_beside_the_empty_tile(self) -> None:
        got = decoded(nodata_tiff(nodata_grid()), "nd.tif")
        assert sorted(set(got.ravel().tolist())) == [0, 7, 9]
        assert int((got == 0).sum()) == 272, "16 from 255 and 256 from the empty tile"
        assert_array_equal(got, nodata_oracle())

    def test_the_path_repository_applies_the_tag_too(self, tmp_path: Path) -> None:
        path = tmp_path / "nd.tif"
        path.write_bytes(nodata_tiff(nodata_grid()))
        repository = TiffDemRepository([path], z_unit="metre")
        w = IndexWindow(row0=0, col0=0, rows=32, cols=32)
        assert_array_equal(repository.load_class_window(path.name, w).array, nodata_oracle())

    def test_a_stored_0_beside_no_data_255_is_refused(self) -> None:
        grid = nodata_grid()
        grid[20, 2] = 0
        with pytest.raises(GeoTiffError) as caught:
            decoded(nodata_tiff(grid), "nd0.tif")
        message = str(caught.value)
        assert "nd0.tif" in message and "255" in message, message
        assert re.search(r"\b0\b", message), message

    def test_a_stored_0_outside_the_window_read_is_not_seen(self) -> None:
        """4.6: the check sees only the window read; rows 0-15 miss (20, 2)."""
        grid = nodata_grid()
        grid[20, 2] = 0
        got = decoded(nodata_tiff(grid), "nd0.tif", rows=16)
        assert_array_equal(got, nodata_oracle()[:16])

    def test_without_a_tag_nothing_changes(self) -> None:
        """P1: the 26a fixture, no `GDAL_NODATA`: its oracle, 356 zeros."""
        data = class_tiff(empty=(0,))
        stream = io.BytesIO(data)
        meta, _, page = read_page(stream)
        assert meta.nodata is None
        w = IndexWindow(row0=0, col0=0, rows=ROWS, cols=48)
        got = decode_class_window(LocalTiffBlocks(page, stream, "c.tif"), meta, w).array
        assert_array_equal(got, expected(data))
        assert int((got == 0).sum()) == 356
