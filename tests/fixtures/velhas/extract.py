"""Cut increment 15c-2's ANADEM fixture out of the tile cache (Q17).

    rasputin fetch anadem-v1 --bbox -44.14 -19.49 -44.04 -19.39   # degrees, EPSG:4674
    python tests/fixtures/velhas/extract.py $RASPUTIN_DATA/cache

(The cache it was cut from was fetched for the whole Velhas piece; the
`fetch` line, which covers the window, was not run for this commit.)

Test support, not production code (`15c-geographic-dem.md`, "Fixtures"). The
window is 320 x 320 nodes of OpenTopography's ANADEM COG
(`anadem_v1_compressed_COG.tif`, DOI 10.5069/G9736P4G), its upper-left node at
COG row ROW0, column COL0. It is written with the COG's own key set (1024 = 2,
1025 = 1 area-registered, 2048 = 4674 with the citation 2049, 2054 = 9102,
2057 and 2059 for GRS80), its spacing, the tie point moved to the window's
first cell, NoData -9999 in tag 42113, float32, Deflate, 128 x 128 tiles, so
reading it needs no `imagecodecs`. The values are the COG's, unchanged.

`catchment.geojson` beside it is derived from this file (see `NOTICE`).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from geographic_fixtures import (
    GEOG_CITATION,
    GEOG_INV_FLATTENING,
    GEOG_SEMI_MAJOR_AXIS,
    GRS80_A,
    GRS80_RF,
    geographic_keys,
    tiff,
)
from tin_engine.io.models import IndexWindow
from tin_engine.io.repository import CacheRepository

ROW0, COL0, SIZE = 124_232, 142_421, 320
OUT = Path(__file__).resolve().parent / "anadem_velhas.tif"


def main(cache: Path) -> None:
    repository = CacheRepository(cache, "anadem-v1")
    (footprint,) = repository.footprints()
    m = footprint.meta
    window = IndexWindow(row0=ROW0, col0=COL0, rows=SIZE, cols=SIZE)
    array = np.asarray(repository.load_window(footprint.name, window).array, dtype=np.float32)
    # The window's first cell's corner: its first node, half a cell out.
    x0 = m.x_min + COL0 * m.delta_x - m.delta_x / 2
    y0 = m.y_max - ROW0 * m.delta_y + m.delta_y / 2
    stream = tiff(
        array,
        x0=x0,
        y0=y0,
        step_x=m.delta_x,
        step_y=m.delta_y,
        shorts=geographic_keys(4674, area=True),
        doubles={GEOG_SEMI_MAJOR_AXIS: GRS80_A, GEOG_INV_FLATTENING: GRS80_RF},
        texts={GEOG_CITATION: "SIRGAS 2000"},
        compression="deflate",
        compressionargs={"level": 9},
        tile=(128, 128),
    )
    OUT.write_bytes(stream.getvalue())
    print(OUT, array.shape, float(array.min()), float(array.max()))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
