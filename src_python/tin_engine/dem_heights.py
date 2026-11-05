"""Each triangle's DEM height shares: its area per height class about its own
plane, and its share over NoData, from the DEM's own pixels.

The classes are one width `class_m` (the DEM's relative vertical accuracy),
centred on 0, `2 * half + 1` of them, `half` the least whole number whose
outer edge `(half + 1/2) * class_m` lies beyond the tolerance, so every DEM
node a refinement checked lies inside them. The mesh's corners are moved once
into the DEM's index space, the triangles cut into chunks, and each chunk's
window of pixels handed with its corners to `_core.dem_height_shares`. Pure:
no file, no CLI. Blocking; the C++ call releases the GIL.
"""

from __future__ import annotations

import math
from typing import Final

import numpy as np
import numpy.typing as npt

from tin_engine import _core
from tin_engine.io.models import HeightTable
from tin_engine.land_cover_fractions import grid_row_areas, index_space
from tin_engine.target_grid import SourceWindows

#: The class width when the DEM's source states no relative vertical accuracy,
#: metres: Kartverket's DTM10 standard deviation, the finest DEM meshed.
DEFAULT_CLASS_M: Final = 1.0
#: The most classes on either side of 0: a `.ply` list holds at most 255.
MAX_HALF: Final = 127
#: Triangles per chunk, each with its own window of pixels.
CHUNK: Final = 2**18


def half_count(tolerance_m: float, class_m: float) -> int:
    """The least whole `n >= 0` with `(n + 1/2) * class_m > tolerance_m`. A
    ValueError unless both are finite and above 0, or when `n > MAX_HALF`."""
    if not all(math.isfinite(v) and v > 0.0 for v in (tolerance_m, class_m)):
        raise ValueError(
            f"the tolerance ({tolerance_m}) and the class width ({class_m}) must be finite "
            "and above 0"
        )
    n = max(0, math.floor(tolerance_m / class_m + 0.5))
    while n > 0 and (n - 0.5) * class_m > tolerance_m:
        n -= 1
    while (n + 0.5) * class_m <= tolerance_m:
        n += 1
    if n > MAX_HALF:
        raise ValueError(
            f"--tolerance {tolerance_m:g} at class width {class_m:g} m needs {2 * n + 1} height "
            f"classes; at most {2 * MAX_HALF + 1}"
        )
    return n


def edges_m(half: int, class_m: float) -> npt.NDArray[np.float64]:
    """The `2 * half + 2` class edges `(j + 1/2) * class_m`, `j = -half - 1 ... half`."""
    return np.asarray((np.arange(-half - 1, half + 1) + 0.5) * class_m, dtype=np.float64)


def height_shares(
    windows: SourceWindows,
    vertices: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.uint32],
    *,
    tolerance_m: float,
    class_m: float,
    mesh_crs: str | None,
    threads: int,
    chunk: int = CHUNK,
) -> tuple[HeightTable, int]:
    """The height table of the mesh `vertices` `(n, 3)`, `triangles` `(T, 3)`
    over `windows`, and the number of triangles with no data. `mesh_crs`
    None: the vertices are in the DEM's CRS; else they are moved from it. A
    vertex past the DEM, or with no image in its CRS, is a ValueError."""
    half = half_count(tolerance_m, class_m)
    meta = windows.meta
    xyz = np.asarray(vertices, dtype=np.float64)
    if mesh_crs is None:
        index = np.column_stack(
            [(xyz[:, 0] - meta.x_min) / meta.delta_x, (meta.y_max - xyz[:, 1]) / meta.delta_y]
        )
        index += 0.5
    else:
        index = index_space(xyz, mesh_crs, meta, threads)
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    corners, heights = index[tri], xyz[:, 2][tri]
    if not np.isfinite(corners).all():
        raise ValueError("a mesh vertex has no position in the DEM's coordinate system")
    flat = corners.reshape(-1, 2)
    if len(flat) and (
        flat.min() < 0.0 or flat[:, 0].max() > meta.cols or flat[:, 1].max() > meta.rows
    ):
        raise ValueError("the mesh reaches outside the DEM")
    top = meta.y_max + meta.delta_y / 2
    parts = []
    for start in range(0, len(tri), chunk):
        piece = corners[start : start + chunk]
        lo = np.floor(piece.reshape(-1, 2).min(axis=0)).astype(np.int64)
        hi = np.ceil(piece.reshape(-1, 2).max(axis=0)).astype(np.int64)
        c0, r0 = int(lo[0]), int(lo[1])
        c1, r1 = min(int(hi[0]), meta.cols), min(int(hi[1]), meta.rows)
        dem = np.ascontiguousarray(windows.window(r0, r1, c0, c1))
        areas = grid_row_areas(
            top - r0 * meta.delta_y, meta.delta_x, meta.delta_y, r1 - r0, meta.crs
        )
        local = piece - np.array([c0, r0], dtype=np.float64)
        parts.append(
            _core.dem_height_shares(
                local,
                heights[start : start + chunk],
                dem,
                meta.nodata,
                areas,
                class_m=class_m,
                half=half,
                threads=threads,
            )
        )
    offsets = np.zeros(len(tri) + 1, dtype=np.uint32)
    at = 0
    for k, (o, *_) in enumerate(parts):
        if at + int(o[-1]) > np.iinfo(np.uint32).max:
            raise ValueError(
                "more than 2^32 - 1 height shares; the file's offsets cannot hold them"
            )
        offsets[k * chunk + 1 : k * chunk + len(o)] = o[1:] + np.uint32(at)
        at += int(o[-1])
    joined = [np.concatenate([p[i] for p in parts]) for i in (1, 2, 3)] if parts else []
    first, shares, nodata = joined or (
        np.zeros(0, np.int32), np.zeros(0, np.float32), np.zeros(0, np.float32)
    )  # fmt: skip
    table = HeightTable(
        offsets=offsets, first=first, shares=shares, nodata=nodata,
        edges_m=edges_m(half, class_m), class_m=class_m,
    )  # fmt: skip
    return table, int((table.count == 0).sum())
