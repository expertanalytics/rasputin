"""The domain types :class:`DomainPolygon` and :class:`CrsBox`, and :func:`check_extent`.

A :class:`CrsBox` is a ``--bbox`` with its CRS.

What ``--domain`` reads, and the rules for reading it, are
``tin_engine.io.domain_file``'s. The domain keeps its own CRS, and
:meth:`DomainPolygon.to_crs` moves it into the DEM's, so a domain need not be
drawn in the DEM's CRS.

Every vertex must lie in the DEM's node rectangle, border included: outside it
there is no bilinear z. That is :func:`check_extent`, run after the transform,
against the mosaic.

Pure: numpy, shapely, pyproj, pydantic, ``tin_engine.crs`` and ``RasterMeta``.
No ``_core``, no typer.
"""

from __future__ import annotations

import math
from typing import Self

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator
from pyproj import CRS
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient

from tin_engine.crs import Xy, crs_label, parse_crs, reprojector, same_crs
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
        transform, once, and re-oriented; edges between vertices stay
        straight lines in ``dst`` and are not densified. The same CRS,
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


#: Equal steps along each side of a `CrsBox` when it is moved: chosen for a
#: box up to a 12 x 14 degree basin, whose moved edges then fall inside the
#: true image by at most 1.6 mm; the shortfall grows with the side squared.
EDGE_STEPS = 2000


class BoxError(ValueError):
    """A box whose image in another CRS cannot be one box there."""


class CrsBox(BaseModel):
    """A box `(x_min, y_min, x_max, y_max)` and the text of its CRS, x first
    (easting or longitude) whatever the CRS's axis order. Finite, not
    inverted (equal sides are a line or a point), and within ±90° of
    latitude in a geographic CRS."""

    model_config = ConfigDict(frozen=True)

    crs: str
    box: tuple[float, float, float, float]

    @model_validator(mode="after")
    def _a_box(self) -> Self:
        x0, y0, x1, y1 = self.box
        crs = parse_crs(self.crs)
        if not all(map(math.isfinite, self.box)):
            raise ValueError(f"every value of the box must be finite, got {self.box}")
        if x0 > x1 or y0 > y1:
            raise ValueError("XMIN must not exceed XMAX, nor YMIN YMAX")
        if crs.is_geographic and (y0 < -90 or y1 > 90):
            raise ValueError(f"latitude must lie within ±90° in {crs_label(crs)}")
        return self

    def outline(self, dst: str | CRS) -> Xy:
        """The box's sides, each at `EDGE_STEPS + 1` equal steps, corners
        repeated, moved point by point into `dst`: `(4 * (EDGE_STEPS + 1), 2)`."""
        x0, y0, x1, y1 = self.box
        s = np.linspace(0.0, 1.0, EDGE_STEPS + 1)
        xs = np.concatenate([x0 + (x1 - x0) * s, np.full_like(s, x1), x1 - (x1 - x0) * s,
                             np.full_like(s, x0)])  # fmt: skip
        ys = np.concatenate([np.full_like(s, y0), y0 + (y1 - y0) * s, np.full_like(s, y1),
                             y1 - (y1 - y0) * s])  # fmt: skip
        xy = np.column_stack([xs, ys])
        return xy if same_crs(self.crs, dst) else reprojector(self.crs, dst)(xy)

    def to_crs(self, dst: str | CRS) -> tuple[float, float, float, float]:
        """The box in `dst`: the numbers untouched in the same CRS, else the
        bounds of `outline(dst)`, which lie on or inside the true image's.
        Refuses an outline point with no image, and, into a geographic CRS,
        an image across ±180° or round a pole (neighbours over 180° apart)."""
        if same_crs(self.crs, dst):
            return self.box
        target = parse_crs(dst)
        xy = self.outline(target)
        where = f"the box {self.box} in {crs_label(self.crs)}"
        if not np.isfinite(xy).all():
            raise BoxError(f"{where} has no image in {crs_label(target)}")
        if target.is_geographic and (np.abs(np.diff(xy[:, 0])) > 180).any():
            raise BoxError(
                f"{where} crosses ±180° or holds a pole in {crs_label(target)}; split it"
            )
        (x0, y0), (x1, y1) = xy.min(axis=0), xy.max(axis=0)
        return float(x0), float(y0), float(x1), float(y1)


def check_extent(domain: DomainPolygon, meta: RasterMeta) -> None:
    """Refuses a domain with a vertex outside ``meta``'s node rectangle, as
    the core's ``cell_of`` has it, with ``domain`` already in the DEM's CRS."""
    polygon = domain.polygon
    _, y_min, x_max, _ = meta.node_box()
    for ring in (polygon.exterior, *polygon.interiors):
        for x, y in ring.coords:
            outside = max(meta.x_min - x, x - x_max, y_min - y, y - meta.y_max)
            if outside > 0:
                raise DomainError(
                    f"vertex ({x}, {y}) is {outside:g} m outside the DEM's node rectangle "
                    f"x {meta.x_min} .. {x_max}, y {y_min} .. {meta.y_max}"
                )
