"""Peak NumPy memory of one `resample` block on the São Francisco basin's
target grid (@architect, 2026-10-02, docs/research/basin-memory-options.md).

Usage, from the worktree root:
    python blockprobe.py <dir holding tin_engine with a built _core> <basin wkt file> <basin geojson>

The grid is `target_grid_for(basin, <wkt>, 30)`. The source is a fake
geographic EPSG:4674 raster at ANADEM's 0.000269494585236 degree spacing over
the basin's lon/lat box grown by 0.5 degree; its `window` is a zero-byte
`np.broadcast_to` view, so the source allocates nothing and the trace holds
only resample's own arrays. One 256-row block (the default `block_rows`) at
the grid's middle row, threads=1, measured by tracemalloc (NumPy allocations;
pyproj's internal C buffers are not traced).
"""

import json
import sys
import time
import tracemalloc

pkg, wkt_file, geojson = sys.argv[1:4]
# An editable install of tin_engine would shadow `pkg`: drop its finder.
sys.meta_path = [f for f in sys.meta_path if "editable" not in repr(f).lower()]
sys.path.insert(0, pkg)

import numpy as np  # noqa: E402
import shapely  # noqa: E402

import tin_engine.target_grid as tg  # noqa: E402
from tin_engine.crs import crs_label, reprojector  # noqa: E402
from tin_engine.domain import DomainPolygon  # noqa: E402
from tin_engine.io.models import RasterMeta  # noqa: E402

target = crs_label(open(wkt_file).read().strip())
geom = shapely.geometry.shape(json.load(open(geojson))["features"][0]["geometry"])
poly = shapely.transform(geom, reprojector("EPSG:4674", target))
grid = tg.target_grid_for(DomainPolygon(polygon=poly, crs=target), target, 30)
print(f"grid {grid.rows} rows x {grid.cols} cols = {grid.rows * grid.cols / 1e9:.3f} G nodes")

x0, y0, x1, y1 = geom.bounds
d = 0.000269494585236
meta = RasterMeta(
    x_min=x0 - 0.5, y_max=y1 + 0.5, delta_x=d, delta_y=d,
    cols=int((x1 - x0 + 1) / d), rows=int((y1 - y0 + 1) / d),
    epsg=4674, crs="EPSG:4674", geographic=True, nodata=-9999.0,
    nodata_source="tag", pixel_is_area=False, vertical_unit_assumed=False,
)  # fmt: skip


class Source:
    meta = meta

    def window(self, r0: int, r1: int, c0: int, c1: int) -> np.ndarray:
        return np.broadcast_to(np.float32(500.0), (r1 - r0, c1 - c0))


block = grid.model_copy(update={"row0": grid.row0 + grid.rows // 2, "rows": 256})
tracemalloc.start()
t = time.time()
out = tg.resample(block, Source(), threads=1, block_rows=256)
_, peak = tracemalloc.get_traced_memory()
print(
    f"one 256-row block, 1 thread: tracemalloc peak {peak / 1e9:.3f} GB, "
    f"{peak / (256 * grid.cols):.1f} B per node; its canvas {out.array.nbytes / 1e9:.3f} GB; "
    f"wall {time.time() - t:.2f} s"
)
