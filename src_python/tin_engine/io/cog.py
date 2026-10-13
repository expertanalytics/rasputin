"""Decode only the blocks a window needs.

A GeoTIFF's full-resolution page is stored as blocks, tiles or strips; a window
of it needs only the blocks it meets. Where those blocks' bytes come from is a
`BlockSource`: an open file (`LocalTiffBlocks`, here) or the tile cache
(`io/repository.py`'s `CachedBlocks`). This module never opens a file.

Every block is decoded into its own slice of one preallocated array, so the
thread count cannot change a value. A block a source lacks is refused before
any block is decoded, and a sparse block (byte count 0, which tifffile would
fill with 0, a valid elevation) is refused: a DEM's NoData must be stored, not
implied.

`decode_class_window` reads a raster of class numbers
(`uint8`, `PixelIsArea`): the same blocks, but an empty block is a run of 0,
no data by MapBiomas's convention, and the window's corner is the cell's.
"""

from __future__ import annotations

import threading
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any, BinaryIO, Protocol

import numpy as np
import tifffile

from .geotiff import _stage, to_metres
from .models import ClassWindow, DemTile, GeoTiffError, IndexWindow, RasterMeta


class BlockSource(Protocol):
    """The full-resolution page (from its header) and its blocks' stored bytes.

    `where` names the source in messages: a file name, or `<source>/<object>`.
    """

    page: tifffile.TiffPage
    where: str

    def block(self, index: int) -> bytes:
        """Block `index`'s stored bytes, exactly."""
        ...

    def missing(self, indices: Sequence[int]) -> tuple[int, ...]:
        """Those of `indices` this source cannot give, ascending."""
        ...


class CacheError(ValueError):
    """A cache that is not the copy it claims to be: no manifest, another
    source's, or a header whose sha256 the manifest does not list."""


class NotCached(CacheError):  # noqa: N818 -- the record's name, pinned by the suite
    """`missing` of the `needed` blocks are not in the cache at `where`;
    `needed` is 0 when the source itself is not there."""

    def __init__(self, missing: int, needed: int, where: str) -> None:
        self.missing, self.needed, self.where = missing, needed, where
        if needed == 0:
            text = f"{where} is not in the cache"
        else:
            text = f"{where}: {missing:,} of the {needed:,} blocks needed are not in the cache"
        super().__init__(text)


class LocalTiffBlocks:
    """Blocks read from an open GeoTIFF stream; it lacks none.

    `seek` and `read` share one stream, so they run under a lock:
    a read is short against a decode, and `os.pread` would exclude `BytesIO`.
    """

    def __init__(self, page: tifffile.TiffPage, stream: BinaryIO, name: str) -> None:
        self.page, self.where = page, name
        self._stream, self._lock = stream, threading.Lock()

    def block(self, index: int) -> bytes:
        offset, count = self.page.dataoffsets[index], self.page.databytecounts[index]
        with self._lock:
            self._stream.seek(offset)
            return self._stream.read(count)

    def missing(self, indices: Sequence[int]) -> tuple[int, ...]:
        return ()


def block_grid(page: tifffile.TiffPage) -> tuple[int, int, int]:
    """A block's (rows, cols) and the blocks per block row (1 for strips)."""
    rows, cols = (int(n) for n in page.chunks[:2])
    return rows, cols, -(-int(page.imagewidth) // cols)


def blocks_meeting(page: tifffile.TiffPage, window: IndexWindow) -> tuple[int, ...]:
    """The indices of the blocks `window` meets, ascending. Pure."""
    rows, cols, across = block_grid(page)
    last_row, last_col = window.row0 + window.rows - 1, window.col0 + window.cols - 1
    return tuple(
        r * across + c
        for r in range(window.row0 // rows, last_row // rows + 1)
        for c in range(window.col0 // cols, last_col // cols + 1)
    )


def decode_window(
    source: BlockSource,
    meta: RasterMeta,
    dtype: np.dtype[Any],
    window: IndexWindow,
    *,
    threads: int = 1,
) -> DemTile:
    """`window` of the page `source` holds, decoded and cast to `dtype`.

    `meta` and `dtype` are the page's, as `read_page` gives them. A window
    outside the raster, or empty, is a `ValueError` (a caller's bug). A block
    `source` lacks is `NotCached`, before any block is read; a sparse or
    undecodable block is a `GeoTiffError` naming `source.where` and the block.

    `threads` decoding workers, one by default (the count is the caller's,
    and no module but `resources.py` reads the machine); the count changes no
    value.
    """
    w, page = window, source.page
    ends = (meta.rows - w.row0 - w.rows, meta.cols - w.col0 - w.cols)
    if min(w.row0, w.col0, *ends) < 0 or min(w.rows, w.cols) < 1:
        raise ValueError(f"{w} is not a non-empty window of {meta.rows} x {meta.cols} nodes")
    indices = blocks_meeting(page, w)
    if absent := source.missing(indices):
        raise NotCached(len(absent), len(indices), source.where)
    for index in indices:
        if page.databytecounts[index] == 0:
            raise GeoTiffError(
                f"{source.where}: block {index} is sparse (byte count 0); "
                "a DEM's NoData must be stored, not implied"
            )
    out: np.ndarray[Any, np.dtype[Any]] = np.empty((w.rows, w.cols), dtype=dtype)
    _decode_into(source, w, indices, out, threads)
    to_metres(out, meta)  # heights in metres
    # The public constructor, which copies `out`: the no-copy route
    # (`DemTile._adopt`) is kept to `assemble` and `resample`, and a test in
    # test_mosaic.py greps the package for its callers. The copy holds a
    # window's heights twice while it is made.
    return DemTile(meta=meta.windowed(w), array=out)


def decode_class_window(
    source: BlockSource, meta: RasterMeta, window: IndexWindow, *, threads: int = 1
) -> ClassWindow:
    """`window` of the class raster `source` holds, as `uint8`.

    `meta` is the page's, as `read_page` gives it, with the upper-left node
    half a cell inside the raster; the window's corner is moved back out to
    the cell's. A block `source` lacks is `NotCached` before any block is
    read; an empty block (byte count 0) is cells of 0, no data.
    """
    w, page = window, source.page
    ends = (meta.rows - w.row0 - w.rows, meta.cols - w.col0 - w.cols)
    if min(w.row0, w.col0, *ends) < 0 or min(w.rows, w.cols) < 1:
        raise ValueError(f"{w} is not a non-empty window of {meta.rows} x {meta.cols} cells")
    if page.dtype != np.uint8:
        raise GeoTiffError(f"{source.where}: classes must be uint8, the file holds {page.dtype}")
    indices = blocks_meeting(page, w)
    if absent := source.missing(indices):
        raise NotCached(len(absent), len(indices), source.where)
    out = np.zeros((w.rows, w.cols), dtype=np.uint8)
    stored = [i for i in indices if page.databytecounts[i] > 0]
    _decode_into(source, w, stored, out, threads)
    x, y = meta.node_xy(w.row0, w.col0)
    return ClassWindow(
        x_min=x - meta.delta_x / 2,
        y_max=y + meta.delta_y / 2,
        delta_x=meta.delta_x,
        delta_y=meta.delta_y,
        array=out,
    )


def _decode_into(
    source: BlockSource,
    w: IndexWindow,
    indices: Sequence[int],
    out: np.ndarray[Any, np.dtype[Any]],
    threads: int,
) -> None:
    """Decode each of `indices` into its own slice of `out`, `w`'s canvas, on
    `threads` workers; none is sparse."""
    page = source.page

    def put(index: int) -> None:
        data = source.block(index)
        with _stage(f"{source.where}, block {index}"):
            segment, (_, _, r0, c0, _), _ = page.decode(data, index)
        assert segment is not None  # a sparse block was refused or skipped by the caller
        rows = slice(max(r0, w.row0), min(r0 + segment.shape[1], w.row0 + w.rows))
        cols = slice(max(c0, w.col0), min(c0 + segment.shape[2], w.col0 + w.cols))
        block = segment[0, rows.start - r0 : rows.stop - r0, cols.start - c0 : cols.stop - c0, 0]
        out[rows.start - w.row0 : rows.stop - w.row0, cols.start - w.col0 : cols.stop - w.col0] = (
            block
        )

    with ThreadPoolExecutor(threads) as pool:
        list(pool.map(put, indices))


__all__ = [
    "BlockSource",
    "CacheError",
    "LocalTiffBlocks",
    "NotCached",
    "block_grid",
    "blocks_meeting",
    "decode_class_window",
    "decode_window",
]
