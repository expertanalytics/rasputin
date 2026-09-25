"""The one adapter from a decoded tile into the core raster view.

Increment 12, R5, and ``project_structure.md``'s ``raster`` section: this is the
only module that builds a core raster. Exactly the array, the four affine
scalars and the optional sentinel cross. Rows and columns come from the array's
shape on the C++ side; the CRS, ``nodata_source``, ``pixel_is_area`` and
``vertical_unit_assumed`` stay in Python, because nothing in C++ reads them.
"""

from __future__ import annotations

from tin_engine._core import RasterView, raster_view
from tin_engine.io.models import DemTile


def to_core(tile: DemTile) -> RasterView:
    """A zero-copy view over ``tile.array``, which it keeps alive.

    The array must be read-only, which a validated ``DemTile`` guarantees: the
    view reads the buffer after this returns, so a writeable one could change
    under it.
    """
    if tile.array.flags.writeable:
        raise ValueError("to_core needs a read-only array; build the tile through DemTile")
    meta = tile.meta
    return raster_view(
        tile.array,
        x_min=meta.x_min,
        y_max=meta.y_max,
        delta_x=meta.delta_x,
        delta_y=meta.delta_y,
        nodata=meta.nodata,
    )
