"""A GeoTIFF encoder for a height grid: array in, bytes out, through tifffile.

Point-registered (PixelIsPoint): node (0, 0), row 0 the northernmost, lies at
(x0, y0); float32, one strip, uncompressed. Pure: no path, no file.
"""

from __future__ import annotations

import io

import numpy as np
import numpy.typing as npt
import tifffile

_TIEPOINT, _SCALE, _GEOKEYS = 33922, 33550, 34735


def encode_geotiff(grid: npt.ArrayLike, x0: float, y0: float, step: float, epsg: int) -> bytes:
    """`grid` as a float32 GeoTIFF in the projected CRS `epsg`, nodes `step` apart."""
    # GeoKey directory: version 1.1.0, three keys: projected model, PixelIsPoint, the EPSG code.
    keys = [1, 1, 0, 3, 1024, 0, 1, 1, 1025, 0, 1, 2, 3072, 0, 1, epsg]
    tags = [
        (_TIEPOINT, "d", 6, (0.0, 0.0, 0.0, x0, y0, 0.0), True),
        (_SCALE, "d", 3, (step, step, 0.0), True),
        (_GEOKEYS, "H", len(keys), keys, True),
    ]
    stream = io.BytesIO()
    with tifffile.TiffWriter(stream) as writer:
        writer.write(np.asarray(grid, dtype=np.float32), extratags=tags)
    return stream.getvalue()
