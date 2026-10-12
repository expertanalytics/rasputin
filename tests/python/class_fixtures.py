"""Synthetic MapBiomas-like class rasters and a per-year catalogue entry (26a).

`docs/increments/26-land-cover.md`, section 3.1 ("The catalogue entry") and
section 10's "26a (Python)". No real MapBiomas bytes: a class raster here is a
40 x 48 `uint8` grid of MapBiomas class numbers in EPSG:4326, `PixelIsArea`
(the tie point is the upper-left cell's outer corner), 16 x 16 tiles (3 x 3
blocks), written by tifffile through `geotiff_fixtures.micro_tiff`. Cell
`(r, c)` of year `y` holds `CLASSES[(7 r + 3 c + y) % 16]`, so two years
differ cell by cell, and the lower-left 10 x 10 cells hold 0 (no data by
MapBiomas's convention; no `GDAL_NODATA` tag, as the real file). An empty
tile (byte count 0, outside Brazil in the real file) is made with
`cog_fixtures.with_zero_byte_count`; `expected` is the pixel oracle, the
written grid with each empty tile's rectangle zeroed, sharing nothing with
the reader.

PINNED BY THESE SUITES (the design names the parts, not their spelling):
- `tin_engine.sources.RemoteSource`: `kind` gains `"one-cog-per-year"`
  (`url_template` with `{year}`); `years: tuple[int, ...] = ()`, ascending,
  non-empty for that kind; `values: Literal["heights", "classes"] =
  "heights"`; `url_for(year=None) -> str`, the last year by default, a
  `ValueError` naming the first and last year for one outside `years`;
  `credit_for(*, url, date) -> str`, the credit with `{date}` (ISO
  `YYYY-MM-DD`) and `{url}` filled in, unchanged for a credit without them.
  `check_values(source, wanted, option) -> None`, a `ValueError` naming the
  source, its `values` and `option` when they differ. The entry is
  `SOURCES["mapbiomas-c11"]`.
- `tin_engine.fetch.plan.FetchRequest.year: int | None = None`; `fetch.run`
  takes the object id from `url_for(year)` as for `one-cog`
  (`brazil_coverage-col11_<year>`), refuses a year outside `years` and a year
  given for a source without years (`FetchError`, before any request), and
  writes the rendered credit of every object fetched into `NOTICE.txt`.
- `tin_engine.io.repository.CachedObject.fetched: datetime.date | None =
  None`: the day the header was first fetched, kept on later runs.
- `tin_engine.io.models.ClassWindow(x_min, y_max, delta_x, delta_y, array)`:
  frozen; `array` a read-only 2-D `uint8` copy; `rows`, `cols`;
  `corner(row, col) -> (x, y)`, the upper-left corner of that cell.
- `tin_engine.io.cog.decode_class_window(source, meta, window, *,
  threads=None) -> ClassWindow`: `meta` as `read_page` gives it (node
  corner shifted inward half a cell), the window's corner moved back out;
  an empty block is zeros; a block the source lacks is `NotCached` before any
  block is read.
- `tin_engine.io.repository.CacheRepository.load_class_window(name,
  window) -> ClassWindow` and `class_strips(name, box, *, rows) ->
  Iterator[ClassWindow]`: the cells whose interior meets `box` (open cells,
  clipped to the raster; none is a `ValueError`), top to bottom in strips of
  at most `rows` rows.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from cog_fixtures import block_rectangles, page_of, with_zero_byte_count
from fetch_fixtures import padded
from geotiff_fixtures import (
    EPSG_WGS84,
    GEOGRAPHIC_TYPE,
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    PIXEL_IS_AREA,
    micro_tiff,
)
from tin_engine.io.models import Bounds

KEY = "test-classes"
REAL_KEY = "mapbiomas-c11"
YEARS = tuple(range(1985, 2026))
LAST_YEAR = 2025
#: The classes seen over the basin ("What was measured", "The MapBiomas file").
CLASSES = (3, 4, 9, 11, 12, 15, 21, 24, 25, 30, 33, 39, 41, 48, 62, 75)
ROWS, COLS = 40, 48
TILE = (16, 16)
#: 1/128 of a degree: every cell corner is exact in binary, so a box edge put
#: on a cell edge lies on it, not one ulp to either side.
LON0, LAT0, STEP = -45.0, -15.0, 0.0078125
AREA_KEYS = {GT_MODEL_TYPE: 2, GT_RASTER_TYPE: PIXEL_IS_AREA, GEOGRAPHIC_TYPE: EPSG_WGS84}
#: The design's credit format (section 3.1), with the test's own words.
CREDIT = "Test classes - Collection 11, accessed on {date} through the link: {url}"
LICENCE = "CC BY-SA 4.0 (test)"
CITE = ("Test citation, 2020.",)


#: The served path's template, as the bucket's: the year in the file name.
PATH_TEMPLATE = "mb/brazil_coverage-col11_{year}.tif"


def object_id(year: int) -> str:
    return f"brazil_coverage-col11_{year}"


def served_path(year: int) -> str:
    return PATH_TEMPLATE.format(year=year)


def class_grid(year: int = LAST_YEAR) -> np.ndarray:
    r, c = np.indices((ROWS, COLS))
    grid = np.asarray(CLASSES, dtype=np.uint8)[(7 * r + 3 * c + year) % len(CLASSES)]
    grid[30:, :10] = 0
    return np.ascontiguousarray(grid)


def class_tiff(
    year: int = LAST_YEAR,
    *,
    empty: Sequence[int] = (),
    compression: str = "deflate",
    predictor: int | None = None,
) -> bytes:
    """One year's raster as GeoTIFF bytes, tiles `empty` with byte count 0."""
    extra: dict[str, Any] = {} if predictor is None else {"predictor": predictor}
    stream = micro_tiff(
        class_grid(year),
        tiepoint=(0.0, 0.0, 0.0, LON0, LAT0, 0.0),
        scale=(STEP, STEP, 0.0),
        geokeys=AREA_KEYS,
        tile=TILE,
        compression=compression,
        **extra,
    )
    data = stream.getvalue()
    for index in empty:
        data = with_zero_byte_count(data, index)
    return data


def expected(data: bytes, year: int = LAST_YEAR) -> np.ndarray:
    """The pixel oracle: the written grid, each empty tile's rectangle zeroed."""
    page = page_of(data)
    grid = class_grid(year)
    for index, (r0, r1, c0, c1) in enumerate(block_rectangles(page)):
        if page.databytecounts[index] == 0:
            grid[r0:r1, c0:c1] = 0
    return grid


def per_year_entry(sources: Any, url_template: str, **changes: Any) -> Any:
    """A `one-cog-per-year` entry of classes, as the design's MapBiomas one."""
    fields: dict[str, Any] = {
        "id": KEY,
        "kind": "one-cog-per-year",
        "url_template": url_template,
        "years": YEARS,
        "values": "classes",
        "crs": "EPSG:4326",
        "nodata": 0.0,
        "credit": CREDIT,
        "licence_note": LICENCE,
        "cite": CITE,
    }
    return sources.RemoteSource(**(fields | changes))


def serve_years(
    server: Any, years: Sequence[int], *, empty: Sequence[int] = ()
) -> dict[int, bytes]:
    """Each year's raster served at `served_path(year)`; the bytes by year."""
    files = {year: class_tiff(year, empty=empty) for year in years}
    for year, data in files.items():
        server.put(served_path(year), padded(data))
    return files


def whole_box() -> Bounds:
    """A `Bounds` a little inside the raster's cells, in degrees."""
    return Bounds(
        x_min=LON0 + 0.3 * STEP,
        y_min=LAT0 - (ROWS - 0.3) * STEP,
        x_max=LON0 + (COLS - 0.3) * STEP,
        y_max=LAT0 - 0.3 * STEP,
    )


def bbox_args() -> tuple[str, str, str, str]:
    """`whole_box` as `--bbox` words."""
    return (
        str(LON0 + 0.3 * STEP),
        str(LAT0 - (ROWS - 0.3) * STEP),
        str(LON0 + (COLS - 0.3) * STEP),
        str(LAT0 - 0.3 * STEP),
    )
