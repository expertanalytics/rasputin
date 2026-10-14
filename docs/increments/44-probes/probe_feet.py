# ruff: noqa: E402, E501, E702, N812, B023
"""Probe for increment 44, P9: a class raster in a feet CRS, through the
height-unit check as merged. Run as <venv python> probe_feet.py <tree> <scratch dir> [<data root>]."""
import io
import os
import sys
from pathlib import Path

os.environ["PROJ_NETWORK"] = "OFF"
tree, scratch = sys.argv[1], Path(sys.argv[2])
sys.path[:0] = [tree + "/tests/python"]
import tin_engine

assert tin_engine.__file__.startswith(tree), tin_engine.__file__
import numpy as np

from geotiff_fixtures import (
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    PIXEL_IS_AREA,
    PROJECTED_CS_TYPE,
    VERTICAL_UNITS,
    micro_tiff,
)
from tin_engine import _core
from tin_engine import land_cover_fractions as F
from tin_engine.io.cog import LocalTiffBlocks, decode_class_window
from tin_engine.io.geotiff import read_page
from tin_engine.io.models import GeoTiffError, IndexWindow
from tin_engine.io.repository import TiffDemRepository

# 8 x 8 cells of 100 US survey feet in EPSG:2229 (NAD83 / California zone 5
# (ftUS)), area-registered, tie (0, 0) at (6,500,000, 2,100,000) ftUS,
# columns 0-3 class 10 and 4-7 class 20; no vertical keys.
p = np.full((8, 8), 10, np.uint8); p[:, 4:] = 20
KEYS = {GT_MODEL_TYPE: 1, GT_RASTER_TYPE: PIXEL_IS_AREA, PROJECTED_CS_TYPE: 2229}


def tiff(keys):
    return micro_tiff(p, tiepoint=(0, 0, 0, 6_500_000.0, 2_100_000.0, 0), scale=(100.0, 100.0, 0),
                      geokeys=keys, tile=(16, 16)).getvalue()


data = tiff(KEYS)
scratch.mkdir(parents=True, exist_ok=True)
path = scratch / "feet_classes.tif"
path.write_bytes(data)
print(f"fixture {path.name}: {len(data)} bytes")


def attempt(label, call):
    try:
        got = call()
        print(f"{label}: read; {got}")
    except GeoTiffError as exc:
        print(f"{label}: GeoTiffError {exc}")


def meta_line(meta):
    return (f"crs {meta.crs} file_crs {meta.file_crs} delta {meta.delta_x!r} "
            f"x_min {meta.x_min!r} y_max {meta.y_max!r} z_unit {meta.z_unit} assumed {meta.vertical_unit_assumed}")


print("== P9 the reader as merged")
attempt("read_page, no z_unit (today's route for a class raster)", lambda: meta_line(read_page(io.BytesIO(data))[0]))
attempt("TiffDemRepository([file]).footprints(), no z_unit", lambda: meta_line(TiffDemRepository([path]).footprints()[0].meta))
attempt("TiffDemRepository([file], z_unit='metre').footprints()", lambda: meta_line(TiffDemRepository([path], z_unit="metre").footprints()[0].meta))
for unit, name in ((9001, "metre"), (9003, "US survey foot")):
    attempt(f"with VerticalUnitsGeoKey 4099 = {unit} ({name}), z_unit='metre'",
            lambda: meta_line(read_page(io.BytesIO(tiff({**KEYS, VERTICAL_UNITS: unit})), z_unit="metre")[0]))

print("== P9 fractions on the metre twin, with z_unit='metre'")
meta, _, page = read_page(io.BytesIO(data), z_unit="metre")
s = io.BytesIO(data); _, _, page = read_page(s, z_unit="metre")
w = decode_class_window(LocalTiffBlocks(page, s, path.name), meta, IndexWindow(row0=0, col0=0, rows=8, cols=8))
cx, cy = meta.x_min - meta.delta_x / 2, meta.y_max + meta.delta_y / 2
# Cells (2, 0), (6, 0), (2, 4): inside the raster, since the metre twin's
# corner is not exact; class 10 holds 6 of its 8 cells of area.
V = np.array([[cx + 2 * meta.delta_x, cy], [cx + 6 * meta.delta_x, cy], [cx + 2 * meta.delta_x, cy - 4 * meta.delta_y]], float)
T = np.array([[0, 2, 1]])
print("triangle (metres, mesh CRS = meta.crs):", V.tolist())
rows = np.full(8, meta.delta_x * meta.delta_y)
local = F.index_space(V, meta.crs, meta)[T]
print("corners in index space:", local.tolist())
o, c, ar = _core.class_areas(local, w.array, rows, threads=1)
print("rows of delta_x * delta_y m2:", rows[0], "| classes", c.tolist(), "fractions", (ar / ar.sum()).tolist(),
      "areas m2", ar.tolist())
print("today's row_areas on it:", F.row_areas(w)[:2].tolist())
path.unlink()

if len(sys.argv) > 3:
    print("== P9b the cached MapBiomas header (argument 3: the data root)")
    import tifffile
    header = (Path(sys.argv[3]) / "cache/mapbiomas-c11/brazil_coverage-col11_2025/header.bin").read_bytes()
    keys = tifffile.TiffFile(io.BytesIO(header)).geotiff_metadata or {}
    print("vertical GeoKeys:", sorted(k for k in keys if "Vertical" in k) or "none")
    same = read_page(io.BytesIO(header))[0] == read_page(io.BytesIO(header), z_unit="metre")[0]
    print("meta with z_unit='metre' equals meta with none:", same)
