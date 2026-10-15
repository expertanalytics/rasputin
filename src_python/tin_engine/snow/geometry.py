"""Each triangle's normal, slope, aspect, plan area and centroid height.

Normals are the receiver frame's: unit, turned up, in true east-north-up,
each taken through the frame's own Jacobian at the triangle's centroid, so the
snow run's slope and aspect and the shadow evaluator see the same normal. A
triangle the frame finds to have zero area gets the vertical normal and is
counted, never refused.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from pyproj import CRS

from tin_engine.radiation.frame import receivers_and_zero_area

F64 = npt.NDArray[np.float64]


@dataclass(frozen=True)
class TriangleGeometry:
    """Per triangle: `normals` `(T, 3)`, `centroid_z`, plan `area` (m²),
    `slope_deg`, `aspect_deg` (true, clockwise from north, in [0, 360));
    `zero_area`, which triangles have no area; `degenerate`, their count."""

    normals: F64
    centroid_z: F64
    area: F64
    slope_deg: F64
    aspect_deg: F64
    zero_area: npt.NDArray[np.bool_]
    degenerate: int


def triangle_geometry(
    vertices: npt.ArrayLike, triangles: npt.ArrayLike, crs: str | CRS
) -> TriangleGeometry:
    """The geometry of `triangles` (`(T, 3)` indices) over `vertices` (`(n, 3)`,
    metres) in the projected `crs`. The receiver frame's refusals, naming the
    triangle by its index in `triangles`, are a ValueError: a non-finite
    coordinate, a vertex index out of range, a centroid too near a pole."""
    v = np.asarray(vertices, dtype=np.float64).reshape(-1, 3)
    t = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    frame, zero_area = receivers_and_zero_area(v, t, crs)
    n = frame.normal
    a = v[t[:, 0]]
    # Edges from the first corner, so large map coordinates cancel before the product.
    d1, d2 = v[t[:, 1]] - a, v[t[:, 2]] - a
    plan = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    return TriangleGeometry(
        normals=n,
        centroid_z=(a[:, 2] + v[t[:, 1], 2] + v[t[:, 2], 2]) / 3.0,
        area=0.5 * np.abs(plan),
        slope_deg=np.degrees(np.arctan2(np.hypot(n[:, 0], n[:, 1]), n[:, 2])),
        aspect_deg=np.degrees(np.arctan2(n[:, 0], n[:, 1])) % 360.0,
        zero_area=zero_area,
        degenerate=int(zero_area.sum()),
    )
