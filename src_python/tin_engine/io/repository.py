"""A DEM held as tiles on disk: headers first, pixels on request (increment 15a).

`docs/increments/15-dem-mosaic.md` R1-R3 and Ola's Q3 ruling: this is the one
module in `io/` that opens files, and it opens them read-only (`"rb"`). It
only stores: which tiles there are (`footprints`, header-only) and one tile's
pixels (`load`). Choosing tiles and stitching them is `tin_engine.mosaic`'s,
the same grid arithmetic for any storage.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Protocol

from pydantic import BaseModel, ConfigDict

from .geotiff import decode_dem, read_meta
from .models import DemTile, GeoTiffError, RasterMeta

#: `from_directory` lists these suffixes, in any case (R2).
TILE_SUFFIXES = frozenset({".tif", ".tiff"})


class TileFootprint(BaseModel):
    """One tile as its header describes it: the file's name and its node grid."""

    model_config = ConfigDict(frozen=True)

    name: str
    meta: RasterMeta


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
            self._footprints = tuple(
                TileFootprint(name=name, meta=self._read(name, read_meta))
                for name in sorted(self._paths)
            )
        return self._footprints

    def load(self, name: str) -> DemTile:
        return self._read(name, decode_dem)

    def _read[T](self, name: str, reader: Callable[..., T]) -> T:
        """Open tile `name` read-only and run `reader`; a refusal names the file."""
        path = self._paths[name]
        with path.open("rb") as stream:
            try:
                return reader(stream, nodata=self._nodata)
            except GeoTiffError as exc:
                raise GeoTiffError(f"{path.name}: {exc}") from exc


__all__ = ["DemRepository", "TiffDemRepository", "TileFootprint"]
