"""Each triangle's normal, slope, aspect, plan area and centroid height.

Normals are unit, per unit inclined area, turned up, and in true
east-north-up: the grid normal is turned about the vertical by the grid
convergence (true azimuth = grid azimuth + convergence), one value for the
whole mesh. A zero-area triangle gets the vertical normal and is counted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

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
    vertices: npt.ArrayLike, triangles: npt.ArrayLike, convergence_deg: float
) -> TriangleGeometry:
    """The geometry of `triangles` (`(T, 3)` indices) over `vertices` (`(n, 3)`,
    metres). A non-finite coordinate is a ValueError."""
    v = np.asarray(vertices, dtype=np.float64).reshape(-1, 3)
    t = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    if not np.isfinite(v).all():
        raise ValueError("vertex coordinates must be finite")
    a = v[t[:, 0]]
    # Edges from the first corner, so large map coordinates cancel before the product.
    cross = np.cross(v[t[:, 1]] - a, v[t[:, 2]] - a).reshape(-1, 3)
    length = np.linalg.norm(cross, axis=1)
    flat = length == 0.0
    n = np.where(flat[:, None], [0.0, 0.0, 1.0], cross / np.where(flat, 1.0, length)[:, None])
    n[n[:, 2] < 0.0] *= -1.0
    g = np.radians(convergence_deg)
    x, y = n[:, 0].copy(), n[:, 1].copy()
    n[:, 0], n[:, 1] = x * np.cos(g) + y * np.sin(g), y * np.cos(g) - x * np.sin(g)
    return TriangleGeometry(
        normals=n,
        centroid_z=(a[:, 2] + v[t[:, 1], 2] + v[t[:, 2], 2]) / 3.0,
        area=0.5 * np.abs(cross[:, 2]),
        slope_deg=np.degrees(np.arctan2(np.hypot(n[:, 0], n[:, 1]), n[:, 2])),
        aspect_deg=np.degrees(np.arctan2(n[:, 0], n[:, 1])) % 360.0,
        zero_area=flat,
        degenerate=int(flat.sum()),
    )
