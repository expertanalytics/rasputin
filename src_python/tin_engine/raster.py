"""The one adapter from a decoded tile into the core raster view.

This is the only module that builds a core raster, so what crosses into C++ is
decided in one place. Exactly the array, the four affine scalars and the
optional sentinel cross. Rows and columns come from the array's
shape on the C++ side; the CRS, ``nodata_source``, ``pixel_is_area`` and
``vertical_unit_assumed`` stay in Python, because nothing in C++ reads them.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from tin_engine._core import RasterView, raster_view
from tin_engine.io.models import DemTile


class SharedArray(np.ndarray[Any, Any]):
    """A read-only array that ``copy.deepcopy`` returns as itself.

    ``dataclasses.asdict`` deep-copies every field, so a plain array would
    cross into C++ as a writeable copy. An array nothing may write needs no
    copy, so the view stays zero-copy over the tile's own buffer.
    """

    def __deepcopy__(self, memo: dict[int, Any] | None = None, /) -> SharedArray:
        return self


@dataclass(frozen=True, slots=True)
class RasterArgs:
    """Exactly what crosses into C++ for a raster view, in ``raster_view``'s
    keyword names: ``raster_view(**dataclasses.asdict(args))``."""

    array: SharedArray
    x_min: float
    y_max: float
    delta_x: float
    delta_y: float
    nodata: float | None


def core_args(tile: DemTile) -> RasterArgs:
    """The tile's array and affine scalars, checked once for any module's ``raster_view``.

    The array must be read-only, which a validated ``DemTile`` guarantees: the
    view reads the buffer after this returns, so a writeable one could change
    under it. A tile in degrees is refused: the core measures in metres.
    """
    if tile.array.flags.writeable:
        raise ValueError("core_args needs a read-only array; build the tile through DemTile")
    meta = tile.meta
    if meta.geographic:  # the core measures in metres; a degree is not one
        raise ValueError(f"a geographic DEM ({meta.crs}, degrees) cannot enter the core")
    shared = tile.array.view(SharedArray)
    shared.flags.writeable = False
    return RasterArgs(shared, meta.x_min, meta.y_max, meta.delta_x, meta.delta_y, meta.nodata)


def to_core(tile: DemTile) -> RasterView:
    """A zero-copy view over ``tile.array``, which it keeps alive; `core_args`'s refusals."""
    return raster_view(**asdict(core_args(tile)))


def to_core_lattice(tile: DemTile) -> RasterView:
    """The tile's array on a unit lattice (x = column, y = -row, both spacings
    1), for the hydrology, which reads only the shape, the values and NoData.
    A tile in degrees is accepted: no coordinate crosses."""
    if tile.array.flags.writeable:
        raise ValueError("to_core_lattice needs a read-only array; build the tile through DemTile")
    nodata = tile.meta.nodata
    return raster_view(tile.array, x_min=0.0, y_max=0.0, delta_x=1.0, delta_y=1.0, nodata=nodata)
