# ruff: noqa: E402, E501, E702, N812, B905
"""Probe for increment 44: run as <venv python> probe_catalogue.py <tree>."""
import io
import os
import sys

os.environ["PROJ_NETWORK"] = "OFF"
tree = sys.argv[1]
sys.path[:0] = [tree + "/tests/python"]
import tin_engine

assert tin_engine.__file__.startswith(tree), tin_engine.__file__
import numpy as np
from pyproj import CRS

import class_fixtures as K
from geotiff_fixtures import (
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    PIXEL_IS_AREA,
    PROJECTED_CS_TYPE,
    micro_tiff,
)
from tin_engine import land_cover_fractions as F
from tin_engine.fetch.run import TILE_SPACING, _tile_name
from tin_engine.io.cog import LocalTiffBlocks, decode_class_window
from tin_engine.io.geotiff import read_page
from tin_engine.io.models import IndexWindow, TileFootprint


class Local:
    """A stand-in for the path repository the design adds: one file's class windows."""
    def __init__(self, data, name="x"):
        self.data, self.name = data, name
        meta, _, self.page = read_page(io.BytesIO(data))
        self.meta = meta
    def footprints(self):
        return (TileFootprint(name=self.name, meta=self.meta, dtype=np.dtype(np.uint8)),)
    def load_class_window(self, name, window):
        s = io.BytesIO(self.data)
        meta, _, page = read_page(s)
        return decode_class_window(LocalTiffBlocks(page, s, name), meta, window)

print("== P1 a class raster read by path (the test fixture, EPSG:4326, 40 x 48, tile 16)")
data = K.class_tiff(empty=(0,))
r = Local(data)
w = r.load_class_window("x", IndexWindow(row0=0, col0=0, rows=K.ROWS, cols=K.COLS))
print("meta crs", r.meta.crs, "geographic", r.meta.geographic, "nodata", r.meta.nodata, r.meta.nodata_source)
print("equals oracle:", bool((w.array == K.expected(data)).all()), "zeros:", int((w.array == 0).sum()))

print("== P2 GDAL_NODATA 255 on a uint8 class raster; an empty tile")
g = np.full((32, 32), 7, np.uint8); g[:, 16:] = 9; g[0:4, 0:4] = 255
d255 = micro_tiff(g, tiepoint=(0, 0, 0, K.LON0, K.LAT0, 0), scale=(K.STEP, K.STEP, 0),
                  geokeys=K.AREA_KEYS, tile=(16, 16), nodata="255").getvalue()
from cog_fixtures import with_zero_byte_count

d255e = with_zero_byte_count(d255, 3)
r2 = Local(d255e)
a = r2.load_class_window("x", IndexWindow(row0=0, col0=0, rows=32, cols=32)).array
print("meta nodata", r2.meta.nodata, r2.meta.nodata_source, "| cells 255:", int((a == 255).sum()),
      "cells 0 (empty tile today):", int((a == 0).sum()), "classes:", sorted(set(a.ravel().tolist())))
try:
    micro_tiff(g, tiepoint=(0, 0, 0, K.LON0, K.LAT0, 0), scale=(K.STEP, K.STEP, 0), geokeys=K.AREA_KEYS,
               nodata="-9999").getvalue()
    read_page(io.BytesIO(micro_tiff(g, tiepoint=(0, 0, 0, K.LON0, K.LAT0, 0), scale=(K.STEP, K.STEP, 0),
               geokeys=K.AREA_KEYS, nodata="-9999").getvalue()))
    print("nodata -9999 on uint8: accepted")
except Exception as e:
    print("nodata -9999 on uint8:", type(e).__name__, str(e)[:100])

print("== P3 a projected class raster (EPSG:3035, 8 x 8 cells of 100 m, classes 10 | 20)")
p = np.full((8, 8), 10, np.uint8); p[:, 4:] = 20
X0, Y0 = 4_000_000.0, 3_000_000.0
keys = {GT_MODEL_TYPE: 1, GT_RASTER_TYPE: PIXEL_IS_AREA, PROJECTED_CS_TYPE: 3035}
dp = micro_tiff(p, tiepoint=(0, 0, 0, X0, Y0, 0), scale=(100.0, 100.0, 0), geokeys=keys, tile=(16, 16)).getvalue()
rp = Local(dp)
print("meta crs", rp.meta.crs, "epsg", rp.meta.epsg, "geographic", rp.meta.geographic, "delta", rp.meta.delta_x)
wp = rp.load_class_window("x", IndexWindow(row0=0, col0=0, rows=8, cols=8))
print("row_areas today on it:", F.row_areas(wp)[:3].tolist())
# One triangle: corner (0,0),(8,0),(0,8) in cells -> class 10 holds 24 of 32 = 0.75.
V = np.array([[X0, Y0], [X0 + 800, Y0], [X0, Y0 - 800]], float)
T = np.array([[0, 2, 1]])
try:
    got = F.assign(V, T, mesh_crs="EPSG:3035", repository=rp, name="x")
    print("assign today:", got.table.classes.tolist(), got.table.fractions.tolist(), got.codes.tolist())
except Exception as e:
    print("assign today:", type(e).__name__, str(e)[:120])
const = np.full(8, 1e4)
from tin_engine import _core

local = F.index_space(V, "EPSG:3035", rp.meta)[T]
o, c, ar = _core.class_areas(local, wp.array, const, threads=1)
print("with constant 1e4 m2 rows:", c.tolist(), (ar / ar.sum()).tolist(), "areas m2", ar.tolist())

print("== P4 the ellipsoid of a geographic class raster: rows of 1/128 deg at 60 N")
for code in ("EPSG:4326", "EPSG:4674", "EPSG:4230", "EPSG:4277"):
    geod = CRS(code).get_geod()
    lons = [0.0, K.STEP, K.STEP, 0.0]
    a0 = abs(geod.polygon_area_perimeter(lons, [60.0, 60.0, 60 - K.STEP, 60 - K.STEP])[0])
    a1 = abs(geod.polygon_area_perimeter(lons, [60.0 - 10 * K.STEP] * 2 + [60 - 11 * K.STEP] * 2)[0])
    print(f"{code} {CRS(code).ellipsoid.name:28s} a={geod.a} 1/f={1/geod.f if geod.f else 0:.6f} "
          f"cell {a0:.6f} m2, 10 rows down {a1:.6f}, ratio {a1/a0:.12f}")

print("== P5 GLO-30 tile names from a template")
tpl = "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM"
ok = True
for lat, lon in [(60, 10), (-15, -45), (0, 0), (-1, -1), (89, 179), (-90, -180)]:
    t = tpl.format(ns="N" if lat >= 0 else "S", ew="E" if lon >= 0 else "W", lat=abs(lat), lon=abs(lon))
    ok &= t == _tile_name(lat, lon)
print("template equals _tile_name on 6 corners:", ok, "| TILE_SPACING today", TILE_SPACING)


print("== P2b a stored class 0 beside GDAL_NODATA 255")
z = g.copy(); z[20, 2] = 0
dz = micro_tiff(z, tiepoint=(0, 0, 0, K.LON0, K.LAT0, 0), scale=(K.STEP, K.STEP, 0),
                geokeys=K.AREA_KEYS, tile=(16, 16), nodata="255").getvalue()
az = Local(dz).load_class_window("x", IndexWindow(row0=0, col0=0, rows=32, cols=32)).array
print("cells 0:", int((az == 0).sum()), "at", [tuple(map(int, v)) for v in zip(*(az == 0).nonzero())])
