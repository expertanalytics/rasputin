"""Block-structured micro-GeoTIFFs and a hand-made tile cache, for 23a-1's suites.

`docs/increments/23-basin-scale.md`, "Windowed source reads (23a-1)" and its
test list W1-W9. A source is a 50 x 70 grid written by `tifffile` through
`geotiff_fixtures.micro_tiff`, with one reduced-resolution page after it (an
overview, as DTM10 and the COGs carry), NoData cells and NaNs scattered. The
variants are the block layouts and codecs `decode_window` must handle alike:
16 x 16 tiles with padded edge tiles, 8-row strips with a short last strip,
Deflate, Deflate with the floating-point predictor, LZW, and int16 (promoted
to float32).

The pixel oracle is `whole_page`: tifffile decoding the full-resolution page
in one call, cast to the dtype the reader promises, then sliced. Nothing in it
shares code with `tin_engine.io.cog`. `BytesBlocks` is a `BlockSource` double
built from tifffile's own offsets, so W2, W3 and W5 can count calls without
trusting `LocalTiffBlocks`.

The cache layout is the record's ("The cache"): `<root>/<source>/manifest.json`,
`<root>/<source>/<object>/header.bin` (the file up to its first block) and
`<root>/<source>/<object>/blocks/<row>/<col>.bin` (each block's own byte range).
The manifest is built with `CacheManifest`, imported where it is used, so this
module collects before 23a-1 lands.

PINNED BY THIS SUITE (the record names the fields but not their spelling):
`CacheManifest(source, crs, rasputin_version, objects, requests)`, with
`requests` a sequence of `{domain_sha256, date}`; `CachedObject(url,
content_length, last_modified, header_sha256, header_bytes, block_shape)`,
`block_shape` the (rows, cols) of one block. `RemoteSource(id, kind, url,
crs, nodata, credit, licence_note)` with `url_template` and `tile_list_url`
optional. `CacheRepository(root, source_id)` reads `<root>/<source_id>/`, and
`CachedBlocks(directory, page, where)` takes the object's directory.
"""

from __future__ import annotations

import datetime
import hashlib
import io
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import tifffile

from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff, needs_codecs

ROWS = 50
COLS = 70
TILE = (16, 16)
ROWS_PER_STRIP = 8
NODATA = -9999
NODATA_TEXT = "-9999"

#: The full-resolution page's file dtype to what a reader returns (11 §7),
#: restated here so the oracle does not import the reader's table.
PROMOTED = {np.dtype(np.int16): np.dtype(np.float32), np.dtype(np.float32): np.dtype(np.float32)}

TILE_BYTE_COUNTS = 325
STRIP_BYTE_COUNTS = 279


def grid(dtype: Any = np.float32) -> np.ndarray:
    """Node (r, c) holds `100 r + c` (+ 0.25 for floats): distinct and exact.

    A fixed scatter of cells holds the NoData sentinel, and for floats a
    second scatter holds NaN, so a window's NoData and NaN must come through
    untouched (W1 compares bytes).
    """
    r, c = np.indices((ROWS, COLS))
    values = 100 * r + c + (0.25 if np.dtype(dtype).kind == "f" else 0)
    out = np.asarray(values).astype(dtype)
    rng = np.random.default_rng(23)
    out.flat[rng.choice(out.size, 40, replace=False)] = NODATA
    if out.dtype.kind == "f":
        out.flat[rng.choice(out.size, 15, replace=False)] = np.nan
    return out


@dataclass(frozen=True)
class Variant:
    """One block layout and codec: the dtype written and `tifffile.write` kwargs."""

    name: str
    dtype: Any
    kwargs: Mapping[str, Any] = field(default_factory=dict)
    codecs: bool = False


VARIANTS: tuple[Variant, ...] = (
    Variant("tiled_deflate", np.float32, {"tile": TILE, "compression": "deflate"}),
    # The record's baseline, as GDAL writes float DEMs; tifffile needs
    # imagecodecs to write predictor 3, and the reader needs it to decode.
    Variant(
        "tiled_deflate_fp_predictor",
        np.float32,
        {"tile": TILE, "compression": "deflate", "predictor": 3},
        codecs=True,
    ),
    Variant("stripped", np.float32, {"rowsperstrip": ROWS_PER_STRIP, "compression": "deflate"}),
    Variant("tiled_lzw", np.float32, {"tile": TILE, "compression": "lzw"}, codecs=True),
    Variant("tiled_int16", np.int16, {"tile": TILE, "compression": "deflate"}),
)


def variant_params() -> list[Any]:
    """`VARIANTS` as pytest params, the codec ones marked `needs_codecs`."""
    return [pytest.param(v, id=v.name, marks=(needs_codecs,) if v.codecs else ()) for v in VARIANTS]


def build(
    variant: Variant,
    *,
    array: np.ndarray | None = None,
    geokeys: Mapping[int, int] | None = None,
    x_min: float = TIE_X,
    y_max: float = TIE_Y,
) -> bytes:
    """`variant` as GeoTIFF bytes: page 0 the grid, page 1 its 2x overview."""
    data = grid(variant.dtype) if array is None else array
    overview = np.ascontiguousarray(data[::2, ::2])
    stream = micro_tiff(
        data,
        tiepoint=(0.0, 0.0, 0.0, x_min, y_max, 0.0),
        geokeys=geokeys,
        nodata=NODATA_TEXT,
        extra_pages=[(overview, 1)],
        **variant.kwargs,
    )
    return stream.getvalue()


def page_of(data: bytes) -> tifffile.TiffPage:
    """The full-resolution page, parsed by tifffile alone."""
    page = tifffile.TiffFile(io.BytesIO(data)).pages.first
    assert isinstance(page, tifffile.TiffPage)
    return page


def whole_page(data: bytes) -> np.ndarray:
    """The oracle: tifffile decodes page 0 whole; cast as the reader promises."""
    with tifffile.TiffFile(io.BytesIO(data)) as tif:
        array = tif.pages.first.asarray()
    return array.astype(PROMOTED[array.dtype])


def sliced(array: np.ndarray, window: Any) -> np.ndarray:
    return array[window.row0 : window.row0 + window.rows, window.col0 : window.col0 + window.cols]


def same_bytes(a: np.ndarray, b: np.ndarray) -> bool:
    """dtype, shape and bytes equal, so NaN matches NaN."""
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


# --------------------------------------------------------------------------
# Block geometry from tifffile's page attributes (the W3 brute force)
# --------------------------------------------------------------------------


def blocks_across(page: tifffile.TiffPage) -> int:
    """Blocks per block row: a tiled page's tile columns, 1 for strips."""
    return int(page.chunked[1])


def block_rectangles(page: tifffile.TiffPage) -> list[tuple[int, int, int, int]]:
    """Every block's pixel rectangle `(r0, r1, c0, c1)`, clipped to the raster."""
    rows, cols = int(page.imagelength), int(page.imagewidth)
    if page.is_tiled:
        height, width = int(page.tilelength), int(page.tilewidth)
    else:
        height, width = int(page.rowsperstrip), cols
    out = []
    for index in range(len(page.dataoffsets)):
        r0 = (index // blocks_across(page)) * height
        c0 = (index % blocks_across(page)) * width
        out.append((r0, min(r0 + height, rows), c0, min(c0 + width, cols)))
    return out


def brute_force_meeting(page: tifffile.TiffPage, window: Any) -> tuple[int, ...]:
    """The indices of every block whose pixel rectangle meets `window`."""
    w_r1, w_c1 = window.row0 + window.rows, window.col0 + window.cols
    return tuple(
        i
        for i, (r0, r1, c0, c1) in enumerate(block_rectangles(page))
        if r0 < w_r1 and window.row0 < r1 and c0 < w_c1 and window.col0 < c1
    )


class BytesBlocks:
    """A `BlockSource` double over a file's bytes, recording each `block(i)`."""

    def __init__(self, data: bytes, where: str = "double.tif") -> None:
        self.data = data
        self.page = page_of(data)
        self.where = where
        self.calls: list[int] = []

    def block(self, index: int) -> bytes:
        self.calls.append(index)
        offset, count = self.page.dataoffsets[index], self.page.databytecounts[index]
        return self.data[offset : offset + count]

    def missing(self, indices: Sequence[int]) -> tuple[int, ...]:
        return ()


class AbsentBlocks(BytesBlocks):
    """Reports `absent` as missing; any `block()` call fails the test (W5)."""

    def __init__(self, data: bytes, absent: Sequence[int]) -> None:
        super().__init__(data, where="cache/obj")
        self.absent = tuple(absent)

    def block(self, index: int) -> bytes:
        raise AssertionError(f"block({index}) was called before the missing check")

    def missing(self, indices: Sequence[int]) -> tuple[int, ...]:
        return tuple(i for i in indices if i in self.absent)


def with_zero_byte_count(data: bytes, index: int) -> bytes:
    """The file with block `index`'s TileByteCounts/StripByteCounts entry set to 0."""
    page = page_of(data)
    tag = page.tags.get(TILE_BYTE_COUNTS) or page.tags[STRIP_BYTE_COUNTS]
    size = {3: 2, 4: 4, 16: 8}[tag.dtype]
    inline = tag.count * size <= 4
    start = (tag.valueoffset if not inline else tag.offset + 8) + index * size
    buffer = bytearray(data)
    assert int.from_bytes(buffer[start : start + size], "little") == page.databytecounts[index]
    buffer[start : start + size] = bytes(size)
    patched = bytes(buffer)
    assert page_of(patched).databytecounts[index] == 0
    return patched


# --------------------------------------------------------------------------
# The cache, written from the file's own byte ranges
# --------------------------------------------------------------------------


def header_prefix(data: bytes) -> bytes:
    """`header.bin`: the file up to its first block (the record, "The cache")."""
    return data[: min(page_of(data).dataoffsets)]


def block_file(object_dir: Path, page: tifffile.TiffPage, index: int) -> Path:
    across = blocks_across(page)
    return object_dir / "blocks" / str(index // across) / f"{index % across}.bin"


def write_object(directory: Path, data: bytes, skip: Sequence[int] = ()) -> Path:
    """One object's directory: `header.bin` and every block but `skip`."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "header.bin").write_bytes(header_prefix(data))
    page = page_of(data)
    for index, (offset, count) in enumerate(
        zip(page.dataoffsets, page.databytecounts, strict=True)
    ):
        if index in skip:
            continue
        path = block_file(directory, page, index)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data[offset : offset + count])
    return directory


def manifest_for(source_id: str, objects: Mapping[str, bytes], crs: str = "EPSG:25833") -> Any:
    """The `CacheManifest` a fetch of `objects` would have written."""
    from tin_engine.io.repository import CachedObject, CacheManifest

    entries = {}
    for object_id, data in objects.items():
        page = page_of(data)
        shape = (
            (int(page.tilelength), int(page.tilewidth))
            if page.is_tiled
            else (int(page.rowsperstrip), int(page.imagewidth))
        )
        header = header_prefix(data)
        entries[object_id] = CachedObject(
            url=f"https://example.invalid/{source_id}/{object_id}.tif",
            content_length=len(data),
            last_modified="Thu, 01 Oct 2026 00:00:00 GMT",
            header_sha256=hashlib.sha256(header).hexdigest(),
            header_bytes=len(header),
            block_shape=shape,
        )
    return CacheManifest(
        source=source_id,
        crs=crs,
        rasputin_version="0.2.0.dev0",
        objects=entries,
        requests=[{"domain_sha256": "0" * 64, "date": datetime.date(2026, 10, 1)}],
    )


def write_cache(
    root: Path,
    source_id: str,
    objects: Mapping[str, bytes],
    *,
    skip: Mapping[str, Sequence[int]] | None = None,
    crs: str = "EPSG:25833",
) -> Path:
    """`<root>/<source_id>/` as a fetch would leave it; returns that directory."""
    directory = root / source_id
    if directory.exists():
        shutil.rmtree(directory)
    for object_id, data in objects.items():
        write_object(directory / object_id, data, (skip or {}).get(object_id, ()))
    manifest = manifest_for(source_id, objects, crs)
    (directory / "manifest.json").write_text(manifest.model_dump_json())
    return directory
