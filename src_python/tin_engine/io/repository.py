"""A DEM held as tiles on disk: headers first, pixels on request.

The one module in `io/` that opens files; the readers and writers elsewhere in
`io/` take bytes, a stream or a connection from here. Reads open files
read-only (`"rb"`); the two writers are `CacheWriter` and `create_sqlite`.

`TiffDemRepository` stores, it does not choose: which tiles there are
(`footprints`, headers only), one tile's pixels (`load`), or only the blocks a
window meets (`load_window`). Choosing and stitching tiles is
`tin_engine.mosaic`'s, the same grid arithmetic for any storage.

`CacheRepository` reads what `rasputin fetch` copied to `<cache>/<source>/`,
offline: DEM tiles, a land-cover class raster (`load_class_window`,
`class_strips`: the cells a box meets, in strips) and the day each object was
first fetched (`CachedObject.fetched`). A source of map-service features is
cached as one GeoPackage, `features.gpkg`, which the vector readers take by
path. `CacheWriter` is the one writer of that directory; the mesh path only
reads it.

Other files opened here for other modules: `open_geopackage` (SQLite cannot
read a Python stream, so a GeoPackage is read through a read-only connection,
and `io/geopackage.py` decodes through it), `create_sqlite` (the one writable
connection: the GeoPackage file `io/geopackage_write.py` fills), `read_json`
(station, river, domain and feature files: their readers take a path but open
nothing) and `read_text` (a WKT domain).
"""

from __future__ import annotations

import datetime
import fcntl
import hashlib
import io
import json
import math
import os
import shutil
import sqlite3
from collections.abc import Callable, Iterable, Iterator, Sequence
from contextlib import closing
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Protocol

import numpy as np
from pydantic import BaseModel, ConfigDict

from tin_engine.crs import crs_label

from .cog import (
    CacheError,
    LocalTiffBlocks,
    NotCached,
    block_grid,
    blocks_meeting,
    decode_class_window,
    decode_window,
)
from .geopackage_write import create_features_table, create_geopackage
from .geotiff import decode_dem, read_header, read_page
from .models import (
    Bounds,
    ClassWindow,
    DemTile,
    GeoTiffError,
    IndexWindow,
    RasterMeta,
    SrsEntry,
    TileFootprint,
    ZUnit,
)

if TYPE_CHECKING:
    import tifffile

    from tin_engine.mosaic import MosaicPlan

#: `from_directory` lists these suffixes, in any case.
TILE_SUFFIXES = frozenset({".tif", ".tiff"})


class DemRepository(Protocol):
    """Storage only: choosing tiles is the caller's. Sync: `load` is blocking
    I/O, and an async caller wraps it in `asyncio.to_thread`."""

    def footprints(self) -> tuple[TileFootprint, ...]:
        """Every tile's header, sorted by name."""
        ...

    def load(self, name: str) -> DemTile:
        """The whole tile `name`, decoded. `KeyError` for a name not listed."""
        ...

    def load_window(self, name: str, window: IndexWindow) -> DemTile:
        """`window` of tile `name`, decoded."""
        ...

    def check(self, plan: MosaicPlan) -> None:
        """Refuse `plan` if this storage cannot give every tile it reads."""
        ...


class TiffDemRepository:
    """GeoTIFF tiles given as paths. Construction reads no file.

    Paths are resolved and deduplicated, so two spellings of one file are one
    tile. A footprint's name is the file name, so two files sharing one are
    refused. `nodata` is the caller's sentinel and `z_unit` its height unit
    for every tile (`decode_dem`).
    """

    def __init__(
        self,
        paths: Iterable[Path],
        *,
        nodata: float | None = None,
        threads: int = 1,
        z_unit: ZUnit | None = None,
    ) -> None:
        self._paths: dict[str, Path] = {}
        for path in sorted({Path(p).resolve() for p in paths}, key=str):
            if path.name in self._paths:
                raise ValueError(
                    f"two tiles named {path.name}: {self._paths[path.name]} and {path}"
                )
            self._paths[path.name] = path
        self._nodata, self._threads, self._z_unit = nodata, threads, z_unit
        self._footprints: tuple[TileFootprint, ...] | None = None

    @classmethod
    def from_directory(
        cls,
        directory: Path,
        *,
        nodata: float | None = None,
        threads: int = 1,
        z_unit: ZUnit | None = None,
    ) -> TiffDemRepository:
        """The `*.tif` and `*.tiff` files directly in `directory`, not recursively.

        Side files (`.tfw`, `.aux.xml`) carry nothing the GeoTIFF tags do not.
        """
        paths = [
            p for p in directory.iterdir() if p.suffix.lower() in TILE_SUFFIXES and p.is_file()
        ]
        if not paths:
            raise ValueError(f"no .tif or .tiff file in the directory {directory}")
        return cls(paths, nodata=nodata, threads=threads, z_unit=z_unit)

    def footprints(self) -> tuple[TileFootprint, ...]:
        """Read every header once, on the first call, and keep them."""
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
        """`window` of tile `name`, decoding only the blocks it meets."""

        def read(stream: Any, **keys: Any) -> DemTile:
            meta, dtype, page = read_page(stream, **keys)
            blocks = LocalTiffBlocks(page, stream, name)
            return decode_window(blocks, meta, dtype, window, threads=self._threads)

        return self._read(name, read)

    def load_class_window(self, name: str, window: IndexWindow) -> ClassWindow:
        """`window` of the class raster `name`, as `decode_class_window` gives it."""

        def read(stream: Any, **keys: Any) -> ClassWindow:
            meta, _, page = read_page(stream, **keys)
            blocks = LocalTiffBlocks(page, stream, name)
            return decode_class_window(blocks, meta, window, threads=self._threads)

        return self._read(name, read)

    def check(self, plan: MosaicPlan) -> None:
        """A file has every block; nothing to check."""

    def _read[T](self, name: str, reader: Callable[..., T]) -> T:
        """Open tile `name` read-only and run `reader`; a refusal names the file."""
        path = self._paths[name]
        with path.open("rb") as stream:
            try:
                return reader(stream, nodata=self._nodata, z_unit=self._z_unit)
            except GeoTiffError as exc:
                raise GeoTiffError(f"{path.name}: {exc}") from exc


class CachedRequest(BaseModel):
    """One fetch asked of the cache: the sha256 of the region as given (the
    domain's WKB and CRS, or the box's numbers and CRS) and the date, so
    equal requests on one day are logged once."""

    model_config = ConfigDict(frozen=True)

    region_sha256: str
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
    #: The day the header was first fetched; kept over later runs.
    fetched: datetime.date | None = None


class CacheManifest(BaseModel):
    """`<cache>/<source>/manifest.json`: what the cache is a copy of
    (`CacheWriter` writes it). The blocks present are the directory's, not the
    manifest's."""

    model_config = ConfigDict(frozen=True)

    source: str
    crs: str
    rasputin_version: str
    objects: dict[str, CachedObject]
    requests: tuple[CachedRequest, ...] = ()


class FeatureManifest(BaseModel):
    """`<cache>/<source>/manifest.json` of a map-service layer: what
    `features.gpkg` is a copy of (its key `id_field`, its columns `fields`)."""

    model_config = ConfigDict(frozen=True)

    source: str
    kind: Literal["arcgis-features"] = "arcgis-features"
    url: str
    crs: str
    id_field: str
    fields: tuple[str, ...]
    rasputin_version: str
    first_fetched: datetime.date
    requests: tuple[CachedRequest, ...] = ()


#: The features table of a cached map-service layer, in `features.gpkg`.
FEATURES = "features"


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

    def __init__(
        self,
        root: Path,
        source: str,
        *,
        nodata: float | None = None,
        threads: int = 1,
        z_unit: ZUnit | None = None,
    ) -> None:
        self._directory, self._source, self._nodata = Path(root) / source, source, nodata
        self._threads, self._z_unit = threads, z_unit
        self._pages: dict[str, tuple[RasterMeta, np.dtype[Any], tifffile.TiffPage]] = {}
        self._fetched: dict[str, datetime.date | None] = {}

    def footprints(self) -> tuple[TileFootprint, ...]:
        if not self._pages:
            path = self._directory / "manifest.json"
            if not path.is_file():
                raise NotCached(0, 0, f"{self._source} ({self._directory})")
            data = path.read_bytes()
            if json.loads(data).get("kind") == "arcgis-features":
                raise CacheError(
                    f"{self._source} holds features; pass its features.gpkg to a vector option"
                )
            manifest = CacheManifest.model_validate_json(data)
            if manifest.source != self._source:
                raise CacheError(f"{path} is the manifest of {manifest.source}, not {self._source}")
            for name, entry in sorted(manifest.objects.items()):
                self._fetched[name] = entry.fetched
                try:
                    header = (self._directory / name / "header.bin").read_bytes()
                except OSError as exc:
                    raise CacheError(
                        f"{self._source}/{name}: no header.bin; re-fetch with --refresh"
                    ) from exc
                if hashlib.sha256(header).hexdigest() != entry.header_sha256:
                    raise CacheError(
                        f"{self._source}/{name}: header.bin is not the one the manifest "
                        "lists; re-fetch with --refresh"
                    )
                try:
                    self._pages[name] = read_page(
                        io.BytesIO(header), nodata=self._nodata, z_unit=self._z_unit
                    )
                except GeoTiffError as exc:
                    raise GeoTiffError(f"{self._source}/{name}: {exc}") from exc
        return tuple(
            TileFootprint(name=name, meta=meta, dtype=dtype)
            for name, (meta, dtype, _) in self._pages.items()
        )

    def fetched(self, name: str) -> datetime.date | None:
        """The day object `name`'s header was first fetched, if recorded."""
        self.footprints()
        return self._fetched[name]

    def load(self, name: str) -> DemTile:
        meta = self._page(name)[0]
        return self.load_window(name, IndexWindow(row0=0, col0=0, rows=meta.rows, cols=meta.cols))

    def load_window(self, name: str, window: IndexWindow) -> DemTile:
        meta, dtype, _ = self._page(name)
        return decode_window(self._blocks(name), meta, dtype, window, threads=self._threads)

    def load_class_window(self, name: str, window: IndexWindow) -> ClassWindow:
        """`window` of the class raster `name`, as `decode_class_window` gives it."""
        meta = self._page(name)[0]
        return decode_class_window(self._blocks(name), meta, window, threads=self._threads)

    def class_strips(self, name: str, box: Bounds, *, rows: int) -> Iterator[ClassWindow]:
        """The cells of `name` whose interior meets `box`, clipped to the
        raster, top to bottom in strips of at most `rows` rows. A box
        meeting no cell is a `ValueError`."""
        meta = self._page(name)[0]
        x0, y0 = meta.x_min - meta.delta_x / 2, meta.y_max + meta.delta_y / 2  # the raster's corner
        c0 = max(math.floor((box.x_min - x0) / meta.delta_x), 0)
        c1 = min(math.ceil((box.x_max - x0) / meta.delta_x), meta.cols)
        r0 = max(math.floor((y0 - box.y_max) / meta.delta_y), 0)
        r1 = min(math.ceil((y0 - box.y_min) / meta.delta_y), meta.rows)
        if c0 >= c1 or r0 >= r1:
            raise ValueError(
                f"{self._source}/{name}: the box {box} in {crs_label(meta.crs)} meets none of "
                "its cells"
            )
        for start in range(r0, r1, rows):
            window = IndexWindow(row0=start, col0=c0, rows=min(rows, r1 - start), cols=c1 - c0)
            yield self.load_class_window(name, window)

    def check(self, plan: MosaicPlan) -> None:
        """One `NotCached` for every block the plan needs and the cache lacks,
        so a run stops before any decode, with the whole count, rather than at
        the first window that lacks a block."""
        missing = needed = 0
        for placement in plan.tiles:
            blocks = self._blocks(placement.name)
            indices = blocks_meeting(blocks.page, placement.source)
            needed, missing = needed + len(indices), missing + len(blocks.missing(indices))
        if missing:
            where = f"{self._source} ({self._directory})"
            box, crs = plan.meta.node_box(), crs_label(plan.meta.crs)
            raise NotCached(missing, needed, where, box=box, crs=crs)

    def _page(self, name: str) -> tuple[RasterMeta, np.dtype[Any], tifffile.TiffPage]:
        self.footprints()
        return self._pages[name]

    def _blocks(self, name: str) -> CachedBlocks:
        page = self._page(name)[2]
        return CachedBlocks(self._directory / name, page, f"{self._source}/{name}")


class CacheWriter:
    """The write side of `<root>/<source>/`. `manifest` and `present` read
    without the lock, for a dry run; every `put_*` and `discard` run inside
    `with writer:`, which holds `.lock` by `flock` (one fetch per source at a
    time: refused if held, and dropped by the OS when the process dies, so
    there is no stale lock) and clears stray `.part` files. Every file is
    written as `<name>.part` and renamed, so a crash never leaves half a file
    under the real name; only the manifest is fsynced, since a fetch writes
    thousands of block files. The exception is `features.gpkg`: created as a
    `.part` and renamed, then written in place by SQLite transactions, whose
    journal guards a crash."""

    def __init__(self, root: Path, source: str) -> None:
        self.directory, self.source = Path(root) / source, source
        self._lock: int | None = None

    def __enter__(self) -> CacheWriter:
        self.directory.mkdir(parents=True, exist_ok=True)
        fd = os.open(self.directory / ".lock", os.O_CREAT | os.O_RDWR, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            os.close(fd)
            lock = self.directory / ".lock"
            raise CacheError(f"{self.source}: another fetch holds {lock}") from exc
        self._lock = fd
        for stray in self.directory.rglob("*.part"):
            stray.unlink()
        return self

    def __exit__(self, *exc: object) -> None:
        if self._lock is not None:
            os.close(self._lock)  # closing the descriptor drops the flock
            self._lock = None

    def manifest(self) -> CacheManifest | None:
        return self._manifest(CacheManifest)

    def feature_manifest(self) -> FeatureManifest | None:
        return self._manifest(FeatureManifest)

    def _manifest[M: BaseModel](self, model: type[M]) -> M | None:
        path = self.directory / "manifest.json"
        return model.model_validate_json(path.read_bytes()) if path.is_file() else None

    def feature_ids(self, pk: str) -> set[int]:
        """The keys `features.gpkg` holds; none when it is absent. Reads without the lock."""
        path = self.directory / "features.gpkg"
        if not path.is_file():
            return set()
        with closing(open_geopackage(path)) as conn:
            return {int(k) for (k,) in conn.execute(f'SELECT "{pk}" FROM {FEATURES}')}

    def features(
        self, srs: SrsEntry, pk: str, columns: Sequence[tuple[str, str]]
    ) -> sqlite3.Connection:
        """A writable connection to `features.gpkg`, created with an empty
        `MULTIPOLYGON` table keyed by `pk` on first use; journalled, so a crash
        keeps every committed transaction. The caller closes it."""
        assert self._lock is not None, "a CacheWriter writes only inside `with`"
        path = self.directory / "features.gpkg"
        if not path.is_file():
            part = path.with_name(path.name + ".part")
            with closing(create_sqlite(part)) as conn:
                create_geopackage(conn, srs)
                create_features_table(conn, FEATURES, srs.srs_id, "MULTIPOLYGON", pk, columns)
                conn.commit()
            os.replace(part, path)
        conn = create_sqlite(path)
        conn.execute("PRAGMA journal_mode=DELETE")
        conn.execute("PRAGMA synchronous=NORMAL")
        return conn

    def present(self, object_id: str, page: tifffile.TiffPage, indices: Sequence[int]) -> set[int]:
        """Those of `indices` whose block file is there with its byte count."""
        blocks = CachedBlocks(self.directory / object_id, page, f"{self.source}/{object_id}")
        return set(indices) - set(blocks.missing(indices))

    def put_header(self, object_id: str, data: bytes) -> None:
        self._put(self.directory / object_id / "header.bin", data)

    def put_block(self, object_id: str, page: tifffile.TiffPage, index: int, data: bytes) -> None:
        directory = self.directory / object_id
        self._put(CachedBlocks(directory, page, object_id)._path(index), data)

    def put_manifest(self, manifest: CacheManifest | FeatureManifest) -> None:
        self._put(self.directory / "manifest.json", manifest.model_dump_json().encode(), sync=True)

    def put_notice(self, text: str) -> None:
        self._put(self.directory / "NOTICE.txt", text.encode())

    def discard(self) -> None:
        """Everything of the source but the lock (`--refresh`)."""
        for child in self.directory.iterdir():
            if child.is_dir():
                shutil.rmtree(child)
            elif child.name != ".lock":
                child.unlink()

    def _put(self, path: Path, data: bytes, *, sync: bool = False) -> None:
        assert self._lock is not None, "a CacheWriter writes only inside `with`"
        part = path.with_name(path.name + ".part")
        part.parent.mkdir(parents=True, exist_ok=True)
        try:
            with part.open("wb") as stream:
                stream.write(data)
                if sync:
                    stream.flush()
                    os.fsync(stream.fileno())
            os.replace(part, path)
        except BaseException:
            part.unlink(missing_ok=True)
            raise


def open_geopackage(path: Path) -> sqlite3.Connection:
    """``path`` as a read-only SQLite connection. The caller closes it.

    A URI with ``mode=ro``: a missing file is refused, not created, and the
    file's SpatiaLite triggers never fire. ``as_uri`` percent-encodes a ``?``,
    ``#`` or space in the name, so they stay the file's.
    """
    return sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)


def create_sqlite(path: Path) -> sqlite3.Connection:
    """A writable SQLite connection to ``path``, created if
    missing. Removing a stale file first, and renaming the finished one, are
    the caller's: this module deletes and writes nothing outside the
    connection it hands back. The caller closes it, in the thread that opened
    it (SQLite connections are thread-bound)."""
    return sqlite3.connect(path)


def read_json(path: Path) -> Any:
    """The JSON document in ``path``, read-only."""
    with Path(path).open("rb") as stream:
        return json.load(stream)


def read_text(path: Path) -> str:
    """The text in ``path``, decoded as UTF-8 (a WKT domain file)."""
    return Path(path).read_text(encoding="utf-8")


__all__ = [
    "CacheError",
    "CacheManifest",
    "CacheRepository",
    "CacheWriter",
    "CachedBlocks",
    "CachedObject",
    "CachedRequest",
    "DemRepository",
    "FeatureManifest",
    "NotCached",
    "TiffDemRepository",
    "open_geopackage",
    "read_json",
    "read_text",
]
