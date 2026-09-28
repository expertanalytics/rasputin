"""Tiles for increment 15a's suites, built from `RasterMeta` and `DemTile` directly.

`15-dem-mosaic.md`, "Tests for @tester": planning and assembly tests need no
file, and a dict is the repository. Everything here imports only what
increment 11 already shipped (`io.models`), so a missing 15a module fails the
tests that use it, not the collection of this helper.

The whole grids are asymmetric (rows != cols, `delta_x` != `delta_y`, every
node a distinct value), so a transposed index or a swapped axis cannot pass.
Values are small integers plus halves, exact in float32.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from typing import Any

import numpy as np

from tin_engine.io.models import DemTile, RasterMeta

X0 = 500_000.0
Y0 = 6_600_000.0
DX = 10.0
DY = 5.0
EPSG = 25833
SENTINEL = -32767.0


def meta(
    *,
    rows: int,
    cols: int,
    x_min: float = X0,
    y_max: float = Y0,
    dx: float = DX,
    dy: float = DY,
    epsg: int = EPSG,
    nodata: float | None = None,
    area: bool = False,
    vertical_unit_assumed: bool = True,
    nodata_source: str | None = None,
) -> RasterMeta:
    source = nodata_source or ("absent" if nodata is None else "tag")
    return RasterMeta(
        x_min=x_min,
        y_max=y_max,
        delta_x=dx,
        delta_y=dy,
        cols=cols,
        rows=rows,
        epsg=epsg,
        nodata=nodata,
        nodata_source=source,
        pixel_is_area=area,
        vertical_unit_assumed=vertical_unit_assumed,
    )


def values(rows: int, cols: int, dtype: Any = np.float32) -> np.ndarray:
    """Node (r, c) holds `100 r + c + 0.5`: distinct, and exact in float32."""
    r, c = np.indices((rows, cols))
    return np.asarray(100 * r + c + 0.5).astype(dtype)


def whole(rows: int = 9, cols: int = 13, **kwargs: Any) -> DemTile:
    """One tile to cut up. Default 9 x 13, point-registered, dx 10, dy 5."""
    dtype = kwargs.pop("dtype", np.float32)
    array = kwargs.pop("array", None)
    grid = values(rows, cols, dtype) if array is None else array
    return DemTile(meta=meta(rows=rows, cols=cols, **kwargs), array=grid)


def piece(source: DemTile, r0: int, r1: int, c0: int, c1: int, **changes: Any) -> DemTile:
    """`source[r0:r1, c0:c1]` as its own tile, placed where a producer would put it.

    `x_min` is computed as `x_min + c0 * dx`, the way a producer's tie point
    would give it; `changes` override `RasterMeta` fields (to plant a defect),
    and `array=` replaces the values.
    """
    m = source.meta
    array = changes.pop("array", None)
    fields = m.model_dump()
    fields.update(
        x_min=m.x_min + c0 * m.delta_x,
        y_max=m.y_max - r0 * m.delta_y,
        rows=r1 - r0,
        cols=c1 - c0,
    )
    fields.update(changes)
    grid = source.array[r0:r1, c0:c1] if array is None else array
    return DemTile(meta=RasterMeta(**fields), array=grid)


def quadrants(source: DemTile, row_cut: int, col_cut: int, overlap: int) -> dict[str, DemTile]:
    """Four tiles, named by compass, sharing `overlap` node lines at each cut.

    `overlap` 0 is area-registered neighbours that abut (disjoint node sets),
    1 is point-registered neighbours sharing a row and a column, more is a
    true overlap. North and west tiles extend `overlap` lines past the cut.
    """
    rows, cols = source.array.shape
    north, south = (0, row_cut + overlap), (row_cut, rows)
    west, east = (0, col_cut + overlap), (col_cut, cols)
    return {
        "nw.tif": piece(source, *north, *west),
        "ne.tif": piece(source, *north, *east),
        "sw.tif": piece(source, *south, *west),
        "se.tif": piece(source, *south, *east),
    }


def blocks(
    source: DemTile, rows: int, cols: int, skip: Iterable[tuple[int, int]] = ()
) -> dict[str, DemTile]:
    """`source` cut into abutting `rows x cols`-node blocks, named `b<i><j>.tif`.

    Blocks at the (block row, block column) pairs in `skip` are left out, to
    make a hole in the tile set.
    """
    out: dict[str, DemTile] = {}
    omitted = set(skip)
    n_rows, n_cols = source.array.shape
    for i in range(n_rows // rows):
        for j in range(n_cols // cols):
            if (i, j) not in omitted:
                out[f"b{i}{j}.tif"] = piece(
                    source, i * rows, (i + 1) * rows, j * cols, (j + 1) * cols
                )
    return out


class Loads:
    """A `load` for `assemble`: a dict lookup that records every call."""

    def __init__(self, tiles: Mapping[str, DemTile]) -> None:
        self.tiles = dict(tiles)
        self.calls: list[str] = []

    def __call__(self, name: str) -> DemTile:
        self.calls.append(name)
        return self.tiles[name]


def never(name: str) -> DemTile:
    """A `load` that must not be reached (I6: header refusals fire before pixels)."""
    raise AssertionError(f"load({name!r}) was called")


def footprints(footprint: Callable[..., Any], tiles: Mapping[str, DemTile]) -> list[Any]:
    """`TileFootprint`s for `tiles`, built with the class the test fixture imported."""
    return [footprint(name=name, meta=tile.meta) for name, tile in tiles.items()]


def same_array(a: np.ndarray, b: np.ndarray) -> bool:
    """Equal bit for bit: dtype, shape and bytes (so NaN matches NaN, -0 != 0)."""
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


# ---------------------------------------------------------------------------
# Ola's Q1 revised (2026-09-28): the oracle for overlaps and seams
# ---------------------------------------------------------------------------


def own_depth(rows: int, cols: int) -> np.ndarray:
    """Each node's distance, in nodes, to its own tile's nearest border.

    `min(r, c, rows - 1 - r, cols - 1 - c)`: 0 on the border. The tile's own
    border is its whole grid (`meta.rows x meta.cols`), never the window a
    request uses of it and never the mosaic's.
    """
    r, c = np.indices((rows, cols))
    return np.minimum.reduce([r, c, rows - 1 - r, cols - 1 - c])


def on_canvas(
    tiles: Mapping[str, DemTile], grid: RasterMeta
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Each tile's valid values (float64) and own depths on `grid`'s nodes.

    NaN and -1 where the tile does not cover the node or holds NoData there
    (NaN, or its sentinel). Tiles are placed by their coordinates, so this
    shares no index arithmetic with `plan_mosaic`.
    """
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, tile in tiles.items():
        m = tile.meta
        row0 = round((grid.y_max - m.y_max) / grid.delta_y)
        col0 = round((m.x_min - grid.x_min) / grid.delta_x)
        value = np.full((grid.rows, grid.cols), np.nan)
        depth = np.full((grid.rows, grid.cols), -1)
        array = np.asarray(tile.array, dtype=np.float64)
        valid = ~np.isnan(array)
        if m.nodata is not None:
            valid &= array != m.nodata
        tile_depth = own_depth(m.rows, m.cols)
        for i, j in zip(*np.nonzero(valid), strict=True):
            r, c = row0 + int(i), col0 + int(j)
            if 0 <= r < grid.rows and 0 <= c < grid.cols:
                value[r, c], depth[r, c] = array[i, j], tile_depth[i, j]
        out[name] = (value, depth)
    return out


def deepest_interior(tiles: Mapping[str, DemTile], grid: RasterMeta) -> np.ndarray:
    """The Q1-revised mosaic, node by node: of the tiles holding a valid value,
    the one the node lies deepest in; ties to the name that sorts first.
    NaN where no tile holds a valid value."""
    placed = on_canvas(tiles, grid)
    expected = np.full((grid.rows, grid.cols), np.nan)
    for r in range(grid.rows):
        for c in range(grid.cols):
            held = [(-d[r, c], name) for name, (v, d) in placed.items() if d[r, c] >= 0]
            if held:
                expected[r, c] = placed[min(held)[1]][0][r, c]
    return expected


def winners(tiles: Mapping[str, DemTile], grid: RasterMeta) -> list[list[str]]:
    """The name of the tile each node's value is taken from, '' for none."""
    placed = on_canvas(tiles, grid)
    return [
        [
            min(((-d[r, c], n) for n, (_, d) in placed.items() if d[r, c] >= 0), default=(0, ""))[1]
            for c in range(grid.cols)
        ]
        for r in range(grid.rows)
    ]


#: Ola, 2026-09-28: "Ignore below 1mm". A seam counts a node only where
#: |a - b| >= 1 mm, 0.001 in the DEM's (metre) units, compared in float64.
SEAM_THRESHOLD = 0.001


def seams_of(
    tiles: Mapping[str, DemTile], grid: RasterMeta
) -> list[tuple[str, str, int, float, float]]:
    """The seam report the rule asks for, pair by pair, sorted by name:
    `(first, second, nodes, largest, median)` over the mosaic's nodes where
    both tiles hold a valid value and `|a - b| >= SEAM_THRESHOLD`;
    `|difference|` in float64, the median of an even count the mean of the
    middle two. Pairs with no such node are left out."""
    placed = on_canvas(tiles, grid)
    out = []
    names = sorted(placed)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            va, vb = placed[a][0], placed[b][0]
            both = ~np.isnan(va) & ~np.isnan(vb)
            differ = both & (np.abs(np.where(both, va - vb, 0.0)) >= SEAM_THRESHOLD)
            if differ.any():
                gaps = np.abs(va[differ] - vb[differ])
                out.append((a, b, int(differ.sum()), float(gaps.max()), float(np.median(gaps))))
    return out


def seams(result: Any) -> list[tuple[str, str, int, float, float]]:
    """A `Mosaic`'s (or `DemInput`'s) seam report as plain tuples."""
    return [(s.first, s.second, s.nodes, s.largest, s.median) for s in result.seams]


def shifted_by(tile: DemTile, offset: float) -> DemTile:
    """`tile` with `offset` added to every value, NoData kept: a tile that
    disagrees with its neighbours everywhere, by exactly `offset`."""
    array = np.asarray(tile.array).copy()
    valid = ~np.isnan(array)
    if tile.meta.nodata is not None:
        valid &= array != tile.meta.nodata
    array[valid] += np.asarray(offset, dtype=array.dtype)
    return DemTile(meta=tile.meta, array=array)
