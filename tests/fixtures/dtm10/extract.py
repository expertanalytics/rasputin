"""Cut increment 15a's real multi-tile fixtures out of Ola's DTM10 archive.

    python tests/fixtures/dtm10/extract.py ../rasputin_data/DTM10_UTM33_20220924

Test support, not production code (`15-dem-mosaic.md`, "Test data"). Every
window is cut from the **archive**, never from the committed benchmark tile,
which is a different release (N5). Source: Kartverket's DTM10, UTM 33, the
archive's 2022-09-24 release, © Kartverket, CC BY 4.0.

Two fixture directories, each usable as ``--dem DIR``:

- ``seam/``: 6400_4 | 6400_1, N3's first edge pair. East-west neighbours whose
  last and first 51 columns are the same nodes. Rows 3072-3327 of both, and
  the 307 columns nearest the seam of each (256 own plus the 51-column
  overlap). Every overlapping node is valid in both and equal (asserted here).
- ``lattices/``: 7707_1, one of N2's eight tiles half a cell east-west off the
  main lattice, over 7707_2, its aligned southern neighbour. 7707_1's last 96
  rows and 7707_2's first 96, so they share 51 node rows (in y), in columns
  that meet in x but are 5 m apart.

Written with the tiles' own georeferencing (tie point shifted to the window,
area-registered, EPSG:25833, NoData -32767 in tag 42113), float32, Deflate, so
reading them needs no `imagecodecs`. The file names are the source tiles'.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import tifffile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "python"))

from geotiff_fixtures import (  # noqa: E402
    GEOGRAPHIC_TYPE,
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    METRE,
    PIXEL_IS_AREA,
    PROJ_LINEAR_UNITS,
    PROJECTED_CS_TYPE,
    micro_tiff,
)

EPSG = 25833
ETRS89 = 4258
NODATA = -32767
KEYS = {
    GT_MODEL_TYPE: 1,
    GT_RASTER_TYPE: PIXEL_IS_AREA,
    GEOGRAPHIC_TYPE: ETRS89,
    PROJECTED_CS_TYPE: EPSG,
    PROJ_LINEAR_UNITS: METRE,
}


def cut(source: Path, rows: slice, cols: slice, target: Path) -> np.ndarray:
    """Write `source[rows, cols]` to `target` with the window's own tie point."""
    with tifffile.TiffFile(source) as tif:
        page = tif.pages.first
        tie = page.tags[33922].value
        scale = page.tags[33550].value
        assert tif.geotiff_metadata["ProjectedCSTypeGeoKey"] == EPSG
        assert int(tif.geotiff_metadata["GTRasterTypeGeoKey"]) == PIXEL_IS_AREA
        assert page.tags[42113].value.strip() == str(NODATA)
        array = page.asarray()
    window = np.ascontiguousarray(array[rows, cols])
    x0 = tie[3] + cols.start * scale[0]
    y0 = tie[4] - rows.start * scale[1]
    stream = micro_tiff(
        window,
        tiepoint=(0.0, 0.0, 0.0, x0, y0, 0.0),
        scale=tuple(scale),
        geokeys=KEYS,
        nodata=str(NODATA),
        compression="deflate",
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(stream.getvalue())
    return window


def main(archive: Path) -> None:
    name = "{}_10m_z33.tif".format
    rows = slice(3072, 3328)
    west = cut(
        archive / name("6400_4"), rows, slice(5051 - 307, 5051), HERE / "seam" / name("6400_4")
    )
    east = cut(archive / name("6400_1"), rows, slice(0, 307), HERE / "seam" / name("6400_1"))
    overlap_west, overlap_east = west[:, -51:], east[:, :51]
    assert (overlap_west != NODATA).all() and (overlap_east != NODATA).all()
    assert np.array_equal(overlap_west, overlap_east), "the seam must agree bit for bit"
    # The control: shifted by one column, the overlap must disagree somewhere.
    assert not np.array_equal(west[:, -52:-1], overlap_east)

    cols = slice(2000, 2096)
    north = cut(
        archive / name("7707_1"), slice(5053 - 96, 5053), cols, HERE / "lattices" / name("7707_1")
    )
    south = cut(archive / name("7707_2"), slice(0, 96), cols, HERE / "lattices" / name("7707_2"))
    assert (north != NODATA).mean() > 0.9 and (south != NODATA).mean() > 0.9
    for path in sorted(HERE.glob("*/*.tif")):
        print(f"{path.relative_to(HERE)}: {path.stat().st_size} bytes")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
