"""The target grid a reprojected DEM is resampled onto, and its check points
(increment 15c-2, `docs/increments/15c-geographic-dem.md` D2-D4).

A `TargetGrid` is a sub-rectangle of one global lattice in the target CRS:
node `(R, K)` sits at `(K h, -R h)` bit for bit (J6), `h` whole metres.
`resample` fills it bilinearly from a `SourceWindows`; `check_point_blocks`
yields the source's own nodes, moved into the target CRS, for the final
check against the source (D5). Every transform goes through
`crs.reprojector` (J5), one transformer per block, so a thread never shares
one. CRS stays here: what reaches `_core` is metres on the grid.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from typing import Annotated, Any, Protocol

import numpy as np
import numpy.typing as npt
import shapely
from pydantic import BaseModel, ConfigDict, Field, StrictInt
from shapely.prepared import prep

from tin_engine.crs import parse_crs, reprojector
from tin_engine.domain import DomainPolygon
from tin_engine.io.models import Bounds, DemTile, RasterMeta, valid_mask

#: Source nodes per check-point block side (D4).
BLOCK = 512

#: Target nodes per resample block (increment 15e, fix 1): bounds one block's
#: transient, ~135 B per node, whatever the grid's width.
BLOCK_NODES = 1 << 20

Block = tuple[npt.NDArray[np.float64], npt.NDArray[np.float32]]


class TargetGrid(BaseModel):
    """`rows x cols` nodes of the global lattice of spacing `spacing` in `crs`,
    the upper-left one global node `(row0, col0)` (D2)."""

    model_config = ConfigDict(frozen=True)

    crs: str
    spacing: Annotated[StrictInt, Field(ge=1, le=1000)]
    row0: int
    col0: int
    rows: int
    cols: int

    def xy(self, r0: int, r1: int) -> npt.NDArray[np.float64]:
        """Rows `r0` to `r1` (exclusive) of nodes, row-major, as `(N, 2)`."""
        r, c = np.indices((r1 - r0, self.cols), dtype=np.float64)
        h = float(self.spacing)
        x, y = (self.col0 + c) * h, -(self.row0 + r0 + r) * h
        return np.column_stack([x.ravel(), y.ravel()])

    def node_box(self) -> tuple[float, float, float, float]:
        """The node rectangle `(x_min, y_min, x_max, y_max)`."""
        h = float(self.spacing)
        x0, y1 = self.col0 * h, -self.row0 * h
        return x0, y1 - (self.rows - 1) * h, x0 + (self.cols - 1) * h, y1


class SourceWindows(Protocol):
    """The source DEM, read by window, never as a whole (D1)."""

    @property
    def meta(self) -> RasterMeta: ...

    def window(self, r0: int, r1: int, c0: int, c1: int) -> npt.NDArray[Any]: ...


class TileWindows:
    """`SourceWindows` over an assembled tile: slices, no copy."""

    def __init__(self, tile: DemTile) -> None:
        self._tile = tile

    @property
    def meta(self) -> RasterMeta:
        return self._tile.meta

    def window(self, r0: int, r1: int, c0: int, c1: int) -> npt.NDArray[Any]:
        return self._tile.array[r0:r1, c0:c1]


def default_spacing(meta: RasterMeta, at: tuple[float, float]) -> int:
    """The source's north-south node spacing in metres at `at` (in the DEM's
    CRS), rounded to whole metres (Q16 (a)); a projected source's own."""
    if not meta.geographic:
        return max(1, round(meta.delta_y))
    lon, lat = at
    half = meta.delta_y / 2
    geod = parse_crs(meta.crs).get_geod()
    assert geod is not None  # a geographic CRS has an ellipsoid
    _, _, metres = geod.inv(lon, lat - half, lon, lat + half)
    return max(1, round(float(metres)))


def target_grid_for(domain: DomainPolygon, target: str, spacing: int) -> TargetGrid:
    """The domain (in `target`) grown by the cell diagonal with mitred
    corners, its bounds snapped outward to the lattice (D2)."""
    h = spacing
    grown = domain.polygon.buffer(math.sqrt(2) * h, join_style="mitre")
    x0, y0, x1, y1 = grown.bounds
    col0, col1 = math.floor(x0 / h), math.ceil(x1 / h)
    row0, row1 = math.floor(-y1 / h), math.ceil(-y0 / h)
    return TargetGrid(
        crs=target, spacing=h, row0=row0, col0=col0, rows=row1 - row0 + 1, cols=col1 - col0 + 1
    )


def source_region(grid: TargetGrid, meta: RasterMeta, grown: Any) -> tuple[Bounds, Any]:
    """The box the grid's nodes read in the DEM's CRS, grown by two source
    cells, and the needed region: `grown` (in the target CRS) moved there and
    grown by a source cell's diagonal. A box across +-180 degrees longitude, or
    reaching a pole, is refused (D2)."""
    n = max(grid.rows, grid.cols)
    x0, y0, x1, y1 = grid.node_box()
    edge = np.linspace(0.0, 1.0, n + 1)
    ring = np.concatenate(
        [
            np.column_stack([x0 + (x1 - x0) * edge, np.full(n + 1, y1)]),
            np.column_stack([np.full(n + 1, x1), y1 - (y1 - y0) * edge]),
            np.column_stack([x1 - (x1 - x0) * edge, np.full(n + 1, y0)]),
            np.column_stack([np.full(n + 1, x0), y0 + (y1 - y0) * edge]),
        ]
    )
    back = reprojector(grid.crs, meta.crs)
    image = back(ring)
    if not np.isfinite(image).all():
        raise ValueError(f"the target grid in {grid.crs} has no image in the DEM's {meta.crs}")
    if meta.geographic and (np.ptp(image[:, 0]) > 180 or np.abs(image[:, 1]).max() >= 90):
        raise ValueError(
            f"the source region in {meta.crs} crosses the antimeridian (180 degrees "
            "longitude) or reaches a pole; neither is read"
        )
    gx, gy = 2 * meta.delta_x, 2 * meta.delta_y
    (bx0, by0), (bx1, by1) = image.min(axis=0), image.max(axis=0)
    box = Bounds(x_min=bx0 - gx, y_min=by0 - gy, x_max=bx1 + gx, y_max=by1 + gy)
    moved = shapely.transform(grown, back)
    needed = moved.buffer(math.hypot(meta.delta_x, meta.delta_y), join_style="mitre")
    return box, needed


def _pool[T](threads: int, work: Callable[[T], Any], items: Iterable[T]) -> Iterator[Any]:
    """`work` over `items` on `threads` workers, in order, `2 threads` at a time."""
    batch, chunk = list(items), max(1, 2 * threads)
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        for lo in range(0, len(batch), chunk):
            yield from pool.map(work, batch[lo : lo + chunk])


def rows_per_block(cols: int) -> int:
    """The most rows of `cols` nodes within `BLOCK_NODES`, and at least one."""
    return max(1, BLOCK_NODES // cols)


def resample(
    grid: TargetGrid, source: SourceWindows, threads: int, *, block_rows: int | None = None
) -> DemTile:
    """The grid's nodes bilinear from the source in its own index space (D3).
    A node whose 2 x 2 stencil touches NoData or leaves the source is NoData:
    the source's sentinel, or NaN without one. Each node from its own
    coordinates alone, so block size and thread count change no value (J3).
    `block_rows` defaults to `rows_per_block(grid.cols)`."""
    rows_per = rows_per_block(grid.cols) if block_rows is None else block_rows
    m = source.meta
    fill = np.nan if m.nodata is None else m.nodata
    canvas = np.empty((grid.rows, grid.cols), dtype=source.window(0, 1, 0, 1).dtype)

    def block(r0: int) -> None:
        r1 = min(r0 + rows_per, grid.rows)
        lonlat = reprojector(grid.crs, m.crs)(grid.xy(r0, r1))
        row, col = m.index_of(lonlat[:, 0], lonlat[:, 1])
        inside = (col >= 0) & (col <= m.cols - 1) & (row >= 0) & (row <= m.rows - 1)
        out = np.full(col.size, fill, dtype=np.float64)
        if inside.any():
            col, row = col[inside], row[inside]
            c0 = np.minimum(np.floor(col), m.cols - 2).astype(np.intp)
            q0 = np.minimum(np.floor(row), m.rows - 2).astype(np.intp)
            wr0, wc0 = int(q0.min()), int(c0.min())
            win = np.asarray(source.window(wr0, int(q0.max()) + 2, wc0, int(c0.max()) + 2))
            i, j = q0 - wr0, c0 - wc0
            a, b = win[i, j].astype(np.float64), win[i, j + 1].astype(np.float64)
            c, d = win[i + 1, j].astype(np.float64), win[i + 1, j + 1].astype(np.float64)
            fx, fy = col - c0, row - q0
            z = (1 - fy) * ((1 - fx) * a + fx * b) + fy * ((1 - fx) * c + fx * d)
            ok = valid_mask(a, m.nodata) & valid_mask(b, m.nodata) & valid_mask(c, m.nodata)
            out[inside] = np.where(ok & valid_mask(d, m.nodata), z, fill)
        canvas[r0:r1] = out.reshape(r1 - r0, grid.cols)

    for _ in _pool(threads, block, range(0, grid.rows, rows_per)):
        pass
    epsg = parse_crs(grid.crs).to_epsg(min_confidence=100)
    h = float(grid.spacing)
    meta = RasterMeta(
        x_min=grid.col0 * h,
        y_max=-grid.row0 * h,
        delta_x=h,
        delta_y=h,
        cols=grid.cols,
        rows=grid.rows,
        epsg=epsg,
        crs=grid.crs,
        nodata=m.nodata,
        nodata_source=m.nodata_source,
        pixel_is_area=False,
        vertical_unit_assumed=m.vertical_unit_assumed,
    )
    # The canvas is adopted, not copied (15e, fix 2): it is C-contiguous of
    # the meta's shape and never handed out writable (R7's condition).
    return DemTile._adopt(meta, canvas)


def check_point_blocks(
    grid: TargetGrid, source: SourceWindows, domain: DomainPolygon, threads: int
) -> Iterator[Block]:
    """The source's valid nodes in the target CRS, one block per 512 x 512
    source nodes, in row-major block order (D4). A block whose outline's image
    misses the domain is skipped whole; non-finite images and points outside
    the grid's node rectangle are dropped. `domain` is in the target CRS."""
    m, h = source.meta, float(grid.spacing)
    # GEOS prepared geometries are not thread-safe, and `prep` prepares its
    # argument in place: each block prepares its own copy (@perf's crash).
    wkb = shapely.to_wkb(domain.polygon)
    x0, y0, x1, y1 = grid.node_box()

    def block(origin: tuple[int, int]) -> Block | None:
        r0, c0 = origin
        r1, c1 = min(r0 + BLOCK, m.rows), min(c0 + BLOCK, m.cols)
        move = reprojector(m.crs, grid.crs)
        rows, cols = np.arange(r0, r1, dtype=np.float64), np.arange(c0, c1, dtype=np.float64)
        top, bottom = np.full(cols.size, rows[0]), np.full(cols.size, rows[-1])
        left, right = np.full(rows.size, cols[0]), np.full(rows.size, cols[-1])
        rr = np.concatenate([top, rows, bottom[::-1], rows[::-1]])
        cc = np.concatenate([cols, right, cols[::-1], left])
        image = move(np.column_stack(m.node_xy(rr, cc)))
        image = image[np.isfinite(image).all(axis=1)]
        hull = shapely.MultiPoint(image).convex_hull.buffer(h) if len(image) else None
        if hull is None or not prep(shapely.from_wkb(wkb)).intersects(hull):
            return None
        z = np.asarray(source.window(r0, r1, c0, c1))
        r, c = np.nonzero(valid_mask(z, m.nodata))
        xy = move(np.column_stack(m.node_xy(r0 + r, c0 + c)))
        keep = np.isfinite(xy).all(axis=1)
        keep &= (xy[:, 0] >= x0) & (xy[:, 0] <= x1) & (xy[:, 1] >= y0) & (xy[:, 1] <= y1)
        return xy[keep], z[r[keep], c[keep]].astype(np.float32)

    origins = [(r, c) for r in range(0, m.rows, BLOCK) for c in range(0, m.cols, BLOCK)]
    for found in _pool(threads, block, origins):
        if found is not None:
            yield found
