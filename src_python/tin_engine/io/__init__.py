"""File formats. Decoding and encoding both live here.

`io/` knows formats and knows nothing about `_core`, so no module in it needs a
compiled extension in the process.

Files are opened in exactly one module, `repository.py`: a DEM in many tiles
has to be listed and read tile by tile, below `cli.py`. It reads files
read-only, and it also opens a GeoPackage as a read-only SQLite connection,
because SQLite cannot read a Python stream; that connection is
`geopackage.py`'s stream. Its two writers are the tile cache's
(`CacheWriter`) and the one writable SQLite connection the GeoPackage writer
fills. Every other module here takes or returns streams, bytes or an open
connection, so it is testable with no filesystem.
"""

from __future__ import annotations

from .ply import write_ply
from .vtk_legacy import write_vtk

__all__ = ["write_ply", "write_vtk"]
