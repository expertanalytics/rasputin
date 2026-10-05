"""The domain polygon's types: :class:`DomainPolygon`, its ``to_crs``, and :func:`check_extent`.

What ``--domain`` reads, and the rules for reading it, are
``tin_engine.io.domain_file``'s (audit PR C, ``docs/increments/python-audit.md``
section 10). The domain keeps its own CRS, and :meth:`DomainPolygon.to_crs`
moves it into the DEM's (increment 15b, ``15-dem-mosaic.md`` R9), replacing
16's must-match rule.

Every vertex must lie in the DEM's node rectangle, border included (U4 a):
outside it there is no bilinear z (R0). That is :func:`check_extent`, run after
the transform, against the mosaic.

Pure: numpy, shapely, pyproj, pydantic, ``tin_engine.crs`` and ``RasterMeta``.
No ``_core``, no typer.
"""

from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict
from pyproj import CRS
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient

from tin_engine.crs import parse_crs, reprojector, same_crs
from tin_engine.io.models import RasterMeta


class DomainError(ValueError):
    """A domain file refused, in words for the person who wrote it."""


class DomainPolygon(BaseModel):
    """The domain: an oriented shapely polygon, and the text of its CRS."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    polygon: Polygon
    crs: str

    def to_crs(self, dst: str | CRS) -> DomainPolygon:
        """The domain in ``dst``: each vertex through pyproj's ``always_xy``
        transform, once, and re-oriented; rings stay straight (R9). The same CRS,
        by ``same_crs``, returns the same polygon labelled ``target.to_string()``.
        """
        source, target = parse_crs(self.crs), parse_crs(dst)
        if same_crs(source, target):
            return self.model_copy(update={"crs": target.to_string()})
        move = reprojector(source, target)
        rings = [np.asarray(r.coords) for r in (self.polygon.exterior, *self.polygon.interiors)]
        moved = [move(r) for r in rings]
        for ring, image in zip(rings, moved, strict=True):
            bad = ~np.isfinite(image).all(axis=1)
            if bad.any():
                x, y = ring[bad][0]
                raise DomainError(
                    f"the domain vertex ({x}, {y}) in {self.crs} has no image in "
                    f"{target.to_string()}"
                )
        polygon = orient(Polygon(moved[0], moved[1:]), sign=1.0)
        return DomainPolygon(polygon=polygon, crs=target.to_string())


def check_extent(domain: DomainPolygon, meta: RasterMeta) -> None:
    """U4 (a): every vertex in ``meta``'s node rectangle, as the core's
    ``cell_of`` has it, with ``domain`` already in the DEM's CRS."""
    polygon = domain.polygon
    x_max = meta.x_min + (meta.cols - 1) * meta.delta_x
    y_min = meta.y_max - (meta.rows - 1) * meta.delta_y
    for ring in (polygon.exterior, *polygon.interiors):
        for x, y in ring.coords:
            outside = max(meta.x_min - x, x - x_max, y_min - y, y - meta.y_max)
            if outside > 0:
                raise DomainError(
                    f"vertex ({x}, {y}) is {outside:g} m outside the DEM's node rectangle "
                    f"x {meta.x_min} .. {x_max}, y {y_min} .. {meta.y_max}"
                )
