"""What is alive after each step of the reprojected open, and what a garbage
collection frees. A measurement script, not production code.

    python probe_alive.py DOMAIN OUT_CRS_WKT_FILE CACHE OUT_JSON

Repeats ``dem_input._open_reprojected``'s steps one by one (same calls, same
arguments) and reads this process's ``phys_footprint`` after each, then after
``gc.collect()``, then after ``malloc_zone_pressure_relief`` (malloc hands
its freed but kept pages back): the mosaic's shape and dtype, the target grid's, and whether
the drop seen during phase 1 of the `run_phases.py` runs is memory
malloc had freed but kept.
"""

from __future__ import annotations

import ctypes
import gc
import json
import math
import os
import sys
from pathlib import Path

from tin_engine.crs import crs_label, reprojector
from tin_engine.domain import DomainPolygon
from tin_engine.io.repository import CacheRepository
from tin_engine.mosaic import assemble, plan_mosaic
from tin_engine.target_grid import TileWindows, default_spacing, resample, source_region, target_grid_for

_libc, _buf = ctypes.CDLL("/usr/lib/libSystem.B.dylib"), ctypes.create_string_buffer(512)


def fp() -> int:
    _libc.proc_pid_rusage(os.getpid(), 4, _buf)
    return int.from_bytes(_buf.raw[72:80], "little")


def main() -> None:
    domain_path, wkt_path, cache, out = (Path(a) for a in sys.argv[1:5])
    rec: dict = {}

    def note(key: str) -> None:
        rec[key] = fp()
        gc.collect()
        rec[key + " +gc"] = fp()
        _libc.malloc_zone_pressure_relief(None, 0)  # malloc returns its free pages
        rec[key + " +relief"] = fp()
        print(f"{key:12s} {rec[key] / 1e9:6.2f} GB, after gc {rec[key + ' +gc'] / 1e9:6.2f} GB,"
              f" after malloc relief {rec[key + ' +relief'] / 1e9:6.2f} GB", flush=True)

    target = crs_label(wkt_path.read_text().strip())
    geo = json.loads(domain_path.read_text())["features"][0]["geometry"]
    import shapely.geometry

    domain = DomainPolygon(polygon=shapely.geometry.shape(geo), crs="EPSG:4674").to_crs(target)
    repo = CacheRepository(cache, "anadem-v1")
    footprints = repo.footprints()
    meta = footprints[0].meta
    c = domain.polygon.centroid
    ((ax, ay),) = reprojector(target, meta.crs)([(c.x, c.y)])
    grid = target_grid_for(domain, target, default_spacing(meta, (ax, ay)))
    grown = domain.polygon.buffer(math.sqrt(2) * grid.spacing, join_style="mitre")
    box, needed = source_region(grid, meta, grown)
    plan = plan_mosaic(footprints, box, needed)
    note("planned")
    mosaic = assemble(plan, repo.load, needed, load_window=repo.load_window)
    a = mosaic.tile.array
    rec.update(mosaic_shape=list(a.shape), mosaic_dtype=str(a.dtype), mosaic_nbytes=a.nbytes,
               mosaic_base_nbytes=a.base.nbytes if a.base is not None else None,
               grid_rows=grid.rows, grid_cols=grid.cols, plan_tiles=len(plan.tiles))  # fmt: skip
    note("assembled")
    tile = resample(grid, TileWindows(mosaic.tile), os.cpu_count() or 1)
    rec.update(tile_nbytes=tile.array.nbytes)
    note("resampled")
    out.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
