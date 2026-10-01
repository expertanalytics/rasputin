"""A DEM held as tiles on disk: headers first, pixels on request (increment 15a).

`docs/increments/15-dem-mosaic.md` R1-R3 and Ola's Q3 ruling: this is the one
module in `io/` that opens files, and it opens them read-only (`"rb"`). It
only stores: which tiles there are (`footprints`, header-only) and one tile's
pixels (`load`). Choosing tiles and stitching them is `tin_engine.mosaic`'s,
the same grid arithmetic for any storage.

Increment 23a-1 adds windowed loads (`load_window`, only the blocks a window
meets) and the tile cache's read side (`CacheRepository`): a source fetched
by `rasputin fetch` (23a-2) is read from `<cache>/<source>/` alone, offline.

Increment 16b (R2, R3) adds :func:`open_geopackage`: SQLite cannot read from a
Python stream, so a GeoPackage's "stream" is a read-only connection opened
here, and `io/geopackage.py` decodes through it.
"""

from __future__ import annotations

import datetime
import hashlib
import io
import sqlite3
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol

import numpy as np
from pydantic import BaseModel, ConfigDict

from .cog import CacheError, LocalTiffBlocks, NotCached, block_grid, blocks_meeting, decode_window
from .geotiff import decode_dem, read_header, read_page
from .models import DemTile, GeoTiffError, IndexWindow, RasterMeta

if TYPE_CHECKING:
    import tifffile

    from tin_engine.mosaic import MosaicPlan

#: `from_directory` lists these suffixes, in any case (R2).
TILE_SUFFIXES = frozenset({".tif", ".tiff"})


class TileFootprint(BaseModel):
    """One tile as its header describes it: the file's name, its node grid, and
    the dtype it decodes to (float32 unless given; S2)."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    name: str
    meta: RasterMeta
    dtype: np.dtype[Any] = np.dtype(np.float32)


class DemRepository(Protocol):
    """Storage only (R1). Sync: `load` is blocking I/O, and an async caller
    wraps it in `asyncio.to_thread`."""

    def footprints(self) -> tuple[TileFootprint, ...]:
        """Every tile's header, sorted by name."""
        ...

    def load(self, name: str) -> DemTile:
        """The whole tile `name`, decoded. `KeyError` for a name not listed."""
        ...


class TiffDemRepository:
    """GeoTIFF tiles given as paths. Construction reads no file.

    Paths are resolved and deduplicated, so two spellings of one file are one
    tile. A footprint's name is the file name, so two files sharing one are
    refused. `nodata` is the caller's sentinel for every tile (`decode_dem`).
    """

    def __init__(self, paths: Iterable[Path], *, nodata: float | None = None) -> None:
        self._paths: dict[str, Path] = {}
        for path in sorted({Path(p).resolve() for p in paths}, key=str):
            if path.name in self._paths:
                raise ValueError(
                    f"two tiles named {path.name}: {self._paths[path.name]} and {path}"
                )
            self._paths[path.name] = path
        self._nodata = nodata
        self._footprints: tuple[TileFootprint, ...] | None = None

    @classmethod
    def from_directory(cls, directory: Path, *, nodata: float | None = None) -> TiffDemRepository:
        """The `*.tif` and `*.tiff` files directly in `directory`, not recursively.

        Side files (`.tfw`, `.aux.xml`) carry nothing the GeoTIFF tags do not.
        """
        paths = [
            p for p in directory.iterdir() if p.suffix.lower() in TILE_SUFFIXES and p.is_file()
        ]
        if not paths:
            raise ValueError(f"no .tif or .tiff file in the directory {directory}")
        return cls(paths, nodata=nodata)

    def footprints(self) -> tuple[TileFootprint, ...]:
        """Read every header once, on the first call, and keep them (R2)."""
        if self._footprints is None:
            headers = {name: self._read(name, read_header) for name in sorted(self._paths)}
            self._footprints = tuple(
                TileFootprint(name=name, meta=meta, dtype=dtype)
                for name, (meta, dtype) in headers.items()
            )
        return self._footprints

    def load(self, name: str) -> DemTile:
        return self._read(name, decode_dem)

    def load_window(self, name: str, window: IndexWindow) -> DemTile:
        """`window` of tile `name`, decoding only the blocks it meets (23a-1)."""

        def read(stream: Any, *, nodata: float | None) -> DemTile:
            meta, dtype, page = read_page(stream, nodata=nodata)
            return decode_window(LocalTiffBlocks(page, stream, name), meta, dtype, window)

        return self._read(name, read)

    def check(self, plan: MosaicPlan) -> None:
        """A file has every block; nothing to check."""

    def _read[T](self, name: str, reader: Callable[..., T]) -> T:
        """Open tile `name` read-only and run `reader`; a refusal names the file."""
        path = self._paths[name]
        with path.open("rb") as stream:
            try:
                return reader(stream, nodata=self._nodata)
            except GeoTiffError as exc:
                raise GeoTiffError(f"{path.name}: {exc}") from exc


class CachedRequest(BaseModel):
    """One fetch asked of the cache: the domain's sha256 and the date."""

    model_config = ConfigDict(frozen=True)

    domain_sha256: str
    date: datetime.date


class CachedObject(BaseModel):
    """One remote file as fetched: its identity and its header prefix's."""

    model_config = ConfigDict(frozen=True)

    url: str
    content_length: int
    last_modified: str
    header_sha256: str
    header_bytes: int
    block_shape: tuple[int, int]


class CacheManifest(BaseModel):
    """`<cache>/<source>/manifest.json`: what the cache is a copy of (23a-2
    writes it). The blocks present are the directory's, not the manifest's."""

    model_config = ConfigDict(frozen=True)

    source: str
    crs: str
    rasputin_version: str
    objects: dict[str, CachedObject]
    requests: tuple[CachedRequest, ...] = ()


class CachedBlocks:
    """One cached object's blocks: `blocks/<row>/<col>.bin` under `directory`.

    A block is present when its file has exactly the header's byte count; a
    `.part` file has another name and is never read.
    """

    def __init__(self, directory: Path, page: tifffile.TiffPage, where: str) -> None:
        self.directory, self.page, self.where = directory, page, where

    def _path(self, index: int) -> Path:
        across = block_grid(self.page)[2]
        return self.directory / "blocks" / str(index // across) / f"{index % across}.bin"

    def block(self, index: int) -> bytes:
        return self._path(index).read_bytes()

    def missing(self, indices: Sequence[int]) -> tuple[int, ...]:
        def present(index: int) -> bool:
            path = self._path(index)
            return path.is_file() and path.stat().st_size == self.page.databytecounts[index]

        return tuple(i for i in indices if not present(i))


class CacheRepository:
    """A cached source, `<root>/<source>/`, as tiles: one per manifest object,
    named by its id. Construction reads no file; `footprints` reads the
    manifest and each object's `header.bin`, and refuses a cache that is not
    the copy it claims to be (`CacheError`) or is absent (`NotCached`)."""

    def __init__(self, root: Path, source: str, *, nodata: float | None = None) -> None:
        self._directory, self._source, self._nodata = Path(root) / source, source, nodata
        self._pages: dict[str, tuple[RasterMeta, np.dtype[Any], tifffile.TiffPage]] = {}

    def footprints(self) -> tuple[TileFootprint, ...]:
        if not self._pages:
            path = self._directory / "manifest.json"
            if not path.is_file():
                raise NotCached(0, 0, f"{self._source} ({self._directory})")
            manifest = CacheManifest.model_validate_json(path.read_bytes())
            if manifest.source != self._source:
                raise CacheError(f"{path} is the manifest of {manifest.source}, not {self._source}")
            for name, entry in sorted(manifest.objects.items()):
                header = (self._directory / name / "header.bin").read_bytes()
                if hashlib.sha256(header).hexdigest() != entry.header_sha256:
                    raise CacheError(
                        f"{self._source}/{name}: header.bin is not the one the manifest "
                        "lists; re-fetch with --refresh"
                    )
                try:
                    self._pages[name] = read_page(io.BytesIO(header), nodata=self._nodata)
                except GeoTiffError as exc:
                    raise GeoTiffError(f"{self._source}/{name}: {exc}") from exc
        return tuple(
            TileFootprint(name=name, meta=meta, dtype=dtype)
            for name, (meta, dtype, _) in self._pages.items()
        )

    def load(self, name: str) -> DemTile:
        meta = self._page(name)[0]
        return self.load_window(name, IndexWindow(row0=0, col0=0, rows=meta.rows, cols=meta.cols))

    def load_window(self, name: str, window: IndexWindow) -> DemTile:
        meta, dtype, _ = self._page(name)
        return decode_window(self._blocks(name), meta, dtype, window)

    def check(self, plan: MosaicPlan) -> None:
        """One `NotCached` for every block the plan needs and the cache lacks,
        before any is decoded (decided 6)."""
        missing = needed = 0
        for placement in plan.tiles:
            blocks = self._blocks(placement.name)
            indices = blocks_meeting(blocks.page, placement.source)
            needed, missing = needed + len(indices), missing + len(blocks.missing(indices))
        if missing:
            raise NotCached(missing, needed, f"{self._source} ({self._directory})")

    def _page(self, name: str) -> tuple[RasterMeta, np.dtype[Any], tifffile.TiffPage]:
        self.footprints()
        return self._pages[name]

    def _blocks(self, name: str) -> CachedBlocks:
        page = self._page(name)[2]
        return CachedBlocks(self._directory / name, page, f"{self._source}/{name}")


def open_geopackage(path: Path) -> sqlite3.Connection:
    """``path`` as a read-only SQLite connection (16b R3). The caller closes it.

    A URI with ``mode=ro``: a missing file is refused, not created, and the
    file's SpatiaLite triggers never fire. ``as_uri`` percent-encodes a ``?``,
    ``#`` or space in the name, so they stay the file's.
    """
    return sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)


__all__ = [
    "CacheError",
    "CacheManifest",
    "CacheRepository",
    "CachedBlocks",
    "CachedObject",
    "DemRepository",
    "NotCached",
    "TiffDemRepository",
    "TileFootprint",
    "open_geopackage",
]
