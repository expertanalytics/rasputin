"""What a fetch takes, decided from headers alone (increment 23a-2). Pure.

`docs/increments/23-basin-scale.md`, "Planning (`fetch/plan.py`, pure)". The
blocks fetched are those meeting the domain's *box* (decided 1), because every
window the mesh decodes is a box: the box in the frame (`out_crs`, else the
source CRS), grown by `margin` times the source's north-south spacing in
metres, moved into the source CRS, then grown by two source cells as 15c-2's
`source_region` grows its own. The missing blocks become byte ranges, runs
with gaps of at most 64 KiB coalesced up to 8 MiB.
"""

from __future__ import annotations

import io
import math
from collections.abc import Collection, Sequence
from typing import Any

import numpy as np
import tifffile
from pydantic import BaseModel, ConfigDict, model_validator

from tin_engine.crs import parse_crs, same_crs, transform_bounds
from tin_engine.domain import DomainPolygon
from tin_engine.fetch.http import FetchError
from tin_engine.io.cog import block_grid, blocks_meeting
from tin_engine.io.geotiff import read_page
from tin_engine.io.models import Bounds, GeoTiffError, IndexWindow, RasterMeta

GAP = 64 * 1024
MAX_RANGE = 8 * 1024 * 1024
#: Metres per degree of latitude, rounded up: a geographic source's spacing.
METRES_PER_DEGREE = 111_320.0


class FetchRequest(BaseModel):
    """One `rasputin fetch`: a source id and exactly one of `domain` and `box`
    (`box` in `out_crs` if given, else in the source CRS)."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    source: str
    domain: DomainPolygon | None = None
    box: Bounds | None = None
    out_crs: str | None = None
    margin: int = 4
    connections: int = 8
    dry_run: bool = False
    refresh: bool = False

    @model_validator(mode="after")
    def _one_region(self) -> FetchRequest:
        if (self.domain is None) == (self.box is None):
            raise ValueError("give exactly one of a domain and a box")
        return self


class ObjectPlan(BaseModel):
    """One remote file: every block meeting the box (ascending), the byte
    ranges of the missing non-sparse ones, and their total length."""

    model_config = ConfigDict(frozen=True)

    object_id: str
    url: str
    blocks: tuple[int, ...]
    ranges: tuple[tuple[int, int], ...]
    bytes: int


class _Prefix(io.BytesIO):
    """A header prefix that raises on any read past its end (F7). It claims to
    be a very long file, so tifffile reads a tag value lying past the prefix
    instead of dropping the tag as out of the file."""

    def seek(self, offset: int, whence: int = 0, /) -> int:
        return super().seek((1 << 50) + offset, 0) if whence == 2 else super().seek(offset, whence)

    def read(self, size: int | None = -1, /) -> bytes:
        if size is None or size < 0 or self.tell() + size > len(self.getbuffer()):
            raise EOFError("a read past the prefix")
        return super().read(size)


def parse_prefix(
    prefix: bytes, *, nodata: float | None
) -> tuple[RasterMeta, np.dtype[Any], tifffile.TiffPage] | None:
    """The header from a prefix of the file, or None while the prefix lacks
    any of it: a read past its end, or an offset or byte count of a block of
    the full-resolution page (tifffile logs a short offset array and goes on,
    so the count is the check). Geographic CRSs are read (decided 6)."""
    try:
        meta, dtype, page = read_page(_Prefix(prefix), nodata=nodata)
        rows, _, across = block_grid(page)
        count = -(-int(page.imagelength) // rows) * across
        if len(page.dataoffsets) != count or len(page.databytecounts) != count:
            return None
    except EOFError:
        return None
    except GeoTiffError as exc:
        if isinstance(exc.__cause__, EOFError):
            return None
        raise
    return meta, dtype, page


def source_box(request: FetchRequest, meta: RasterMeta) -> Bounds:
    """The box whose blocks are fetched, in the source CRS (`meta.crs`)."""
    source = parse_crs(meta.crs)
    frame = parse_crs(request.out_crs) if request.out_crs else source
    if request.domain is not None:
        x0, y0, x1, y1 = request.domain.to_crs(frame).polygon.bounds
    else:
        assert request.box is not None
        b = request.box
        x0, y0, x1, y1 = b.x_min, b.y_min, b.x_max, b.y_max
    if same_crs(frame, source):
        grow = request.margin * meta.delta_y
    else:
        grow = request.margin * meta.delta_y * (METRES_PER_DEGREE if source.is_geographic else 1)
        grown = (x0 - grow, y0 - grow, x1 + grow, y1 + grow)
        (x0, y0, x1, y1), grow = transform_bounds(frame, source, grown), 0.0
    box = (
        x0 - grow - 2 * meta.delta_x,
        y0 - grow - 2 * meta.delta_y,
        x1 + grow + 2 * meta.delta_x,
        y1 + grow + 2 * meta.delta_y,
    )
    if not all(map(math.isfinite, box)) or (
        source.is_geographic and (box[0] < -180 or box[2] > 180 or box[1] < -90 or box[3] > 90)
    ):
        raise FetchError(f"the box {box} crosses ±180° or a pole, or has no image in the source")
    return Bounds(x_min=box[0], y_min=box[1], x_max=box[2], y_max=box[3])


def plan_object(
    object_id: str,
    url: str,
    page: tifffile.TiffPage,
    meta: RasterMeta,
    box: Bounds,
    present: Collection[int],
) -> ObjectPlan:
    """The blocks with a node strictly within one cell of `box` (the box
    snapped outward to the lattice, as 23a-1's windows are), and the ranges of
    those neither `present` nor sparse (a sparse block is never requested,
    decided 7)."""
    xs, ys = meta.node_xy(np.arange(meta.rows), np.arange(meta.cols))
    cols = np.flatnonzero((xs > box.x_min - meta.delta_x) & (xs < box.x_max + meta.delta_x))
    rows = np.flatnonzero((ys > box.y_min - meta.delta_y) & (ys < box.y_max + meta.delta_y))
    if not len(cols) or not len(rows):
        raise FetchError(f"{object_id}: the box {box} meets no block of {url}")
    window = IndexWindow(
        row0=int(rows[0]),
        col0=int(cols[0]),
        rows=int(rows[-1] - rows[0] + 1),
        cols=int(cols[-1] - cols[0] + 1),
    )
    blocks = blocks_meeting(page, window)
    spans = [
        (int(page.dataoffsets[i]), int(page.dataoffsets[i]) + int(page.databytecounts[i]))
        for i in blocks
        if i not in present and page.databytecounts[i] > 0
    ]
    ranges = coalesce(spans)
    return ObjectPlan(
        object_id=object_id,
        url=url,
        blocks=blocks,
        ranges=ranges,
        bytes=sum(b - a for a, b in ranges),
    )


def coalesce(spans: Sequence[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Block spans sorted by offset, runs whose gaps are at most 64 KiB joined
    into ranges of at most 8 MiB; a larger single block is its own range."""
    out: list[tuple[int, int]] = []
    for start, stop in sorted(spans):
        if out and start - out[-1][1] <= GAP and stop - out[-1][0] <= MAX_RANGE:
            out[-1] = (out[-1][0], max(out[-1][1], stop))  # a span may lie inside another
        else:
            out.append((start, stop))
    return tuple(out)


__all__ = ["FetchRequest", "ObjectPlan", "coalesce", "parse_prefix", "plan_object", "source_box"]
