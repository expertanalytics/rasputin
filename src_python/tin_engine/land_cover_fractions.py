"""Exact land-cover class fractions per triangle (increment 26b).

`docs/increments/26-land-cover.md`, 1.1 "What is computed, and in which
coordinates", 1.2 and section 6. Every mesh vertex is moved once to the class
raster's longitude and latitude (`crs.reprojector`, the one place a transformer
is built) and then to the raster's index space, where every cell is a unit
square; the triangles are cut into chunks, each chunk's class window read from
the cache and handed with the triangle corners to `_core.class_areas`, whose
exact areas become fractions. The cutoff and the ledger are increment 26c's;
this module writes the exact fractions.

Not `land_cover.py`: 16c's `landcover.py` labels triangles from polygons and is
not replaced. Blocking; the C++ call releases the GIL, so an API worker runs
`assign` in `asyncio.to_thread`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from pyproj import Geod

from tin_engine import _core
from tin_engine.crs import reprojector
from tin_engine.io.models import ClassFractions, ClassWindow, IndexWindow, RasterMeta
from tin_engine.io.repository import CacheRepository

#: Triangles per chunk, each chunk with its own class window (1.3).
CHUNK = 2**18

_GEOD = Geod(ellps="WGS84")


@dataclass(frozen=True, slots=True)
class Assigned:
    """The fractions table and, per triangle, the largest class (ties to the
    smaller class number), which is the file's `land_cover_code`."""

    table: ClassFractions
    codes: npt.NDArray[np.int32]


def row_areas(window: ClassWindow) -> npt.NDArray[np.float64]:
    """The WGS 84 ellipsoidal area in m2 of one cell of each row of `window`
    (1.1: computed once per row by `pyproj.Geod`). A cell's area does not
    depend on its longitude, so every row is measured at longitude 0."""
    lons = [0.0, window.delta_x, window.delta_x, 0.0]
    out = np.empty(window.rows, dtype=np.float64)
    for r in range(window.rows):
        top = window.y_max - r * window.delta_y
        area, _ = _GEOD.polygon_area_perimeter(
            lons, [top, top, top - window.delta_y, top - window.delta_y]
        )
        out[r] = abs(float(area))
    return out


def index_space(
    vertices: npt.ArrayLike, mesh_crs: str, meta: RasterMeta
) -> npt.NDArray[np.float64]:
    """Every vertex's `(col, row)` in the class raster `meta` describes: moved
    once to the raster's CRS, then `col = (lon - lon0) / s`, `row = (lat0 - lat) / s`
    with `(lon0, lat0)` the raster's upper-left cell corner (1.1)."""
    xy = np.asarray(vertices, dtype=np.float64)[:, :2]
    lonlat = reprojector(mesh_crs, meta.crs)(xy)
    corner_x, corner_y = meta.x_min - meta.delta_x / 2, meta.y_max + meta.delta_y / 2
    cols = (lonlat[:, 0] - corner_x) / meta.delta_x
    rows = (corner_y - lonlat[:, 1]) / meta.delta_y
    return np.column_stack([cols, rows])


def window_of(corners: npt.NDArray[np.float64], meta: RasterMeta) -> IndexWindow:
    """The cells holding every corner plus one cell on each side, clipped to
    the raster; a mesh reaching past the raster is a `ValueError`."""
    low, high = corners.reshape(-1, 2).min(axis=0), corners.reshape(-1, 2).max(axis=0)
    if not (low.min() >= 0.0 and high[0] <= meta.cols and high[1] <= meta.rows):
        raise ValueError("the mesh reaches outside the land-cover raster")
    lo, hi = np.floor(low).astype(np.int64) - 1, np.ceil(high).astype(np.int64) + 1
    c0, r0 = max(int(lo[0]), 0), max(int(lo[1]), 0)
    c1, r1 = min(int(hi[0]), meta.cols), min(int(hi[1]), meta.rows)
    return IndexWindow(row0=r0, col0=c0, rows=r1 - r0, cols=c1 - c0)


def assign(
    vertices: npt.ArrayLike,
    triangles: npt.ArrayLike,
    *,
    mesh_crs: str,
    repository: CacheRepository,
    name: str,
    chunk: int = CHUNK,
    threads: int = 0,
) -> Assigned:
    """The exact class fractions of every triangle from the cached class raster
    `name` of `repository`, in chunks of `chunk` triangles in row order."""
    metas = {f.name: f.meta for f in repository.footprints()}
    if name not in metas:
        raise ValueError(f"{name} is not in the cache")
    meta = metas[name]
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    corners = index_space(vertices, mesh_crs, meta)[tri]
    parts: list[tuple[npt.NDArray[np.uint32], npt.NDArray[np.uint8], npt.NDArray[np.float64]]] = []
    for start in range(0, len(tri), chunk):
        piece = corners[start : start + chunk]
        window = window_of(piece, meta)
        classes = repository.load_class_window(name, window)
        local = piece - np.array([window.col0, window.row0], dtype=np.float64)
        parts.append(_core.class_areas(local, classes.array, row_areas(classes), threads=threads))
    offsets = np.zeros(len(tri) + 1, dtype=np.uint32)
    at = 0
    for t, (o, _, _) in enumerate(parts):
        if at + int(o[-1]) > np.iinfo(np.uint32).max:
            raise ValueError(
                "more than 2^32 - 1 class entries; the file's offsets cannot hold them"
            )
        offsets[t * chunk + 1 : t * chunk + len(o)] = o[1:] + np.uint32(at)
        at += int(o[-1])
    classes_all = np.concatenate([p[1] for p in parts] or [np.zeros(0, dtype=np.uint8)])
    areas = np.concatenate([p[2] for p in parts] or [np.zeros(0)])
    return Assigned(
        fractions_of(offsets, classes_all, areas), dominant(offsets, classes_all, areas)
    )


def fractions_of(
    offsets: npt.NDArray[np.uint32], classes: npt.NDArray[np.uint8], areas: npt.NDArray[np.float64]
) -> ClassFractions:
    """`f[t, k] = a[t, k] / sum_k a[t, k]` (1.1), as the file's float32."""
    counts = offsets[1:].astype(np.int64) - offsets[:-1].astype(np.int64)
    if (counts == 0).any():
        raise ValueError("a triangle has no class area; it is below the arithmetic's resolution")
    totals = np.add.reduceat(areas, offsets[:-1].astype(np.int64)) if len(areas) else areas
    fractions = (areas / np.repeat(totals, counts)).astype(np.float32)
    return ClassFractions(offsets=offsets, classes=classes, fractions=fractions)


def dominant(
    offsets: npt.NDArray[np.uint32], classes: npt.NDArray[np.uint8], areas: npt.NDArray[np.float64]
) -> npt.NDArray[np.int32]:
    """Per triangle, the class of the largest area, ties to the smaller class
    (1.4, step 7); classes are ascending within a row, so the first maximum."""
    starts = offsets[:-1].astype(np.int64)
    if len(areas) == 0:
        return np.zeros(len(starts), dtype=np.int32)
    counts = offsets[1:].astype(np.int64) - offsets[:-1].astype(np.int64)
    top = np.repeat(np.maximum.reduceat(areas, starts), counts)
    index = np.where(areas == top, np.arange(len(areas)), len(areas))
    return classes[np.minimum.reduceat(index, starts)].astype(np.int32)
