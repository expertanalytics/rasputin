"""The receiver frame of a mesh: per triangle, the numbers the core needs.

The mesh lives in a projected computation frame; the sun lives in true
directions on the ellipsoid. Each triangle becomes a receiver at its
centroid, with its geodetic latitude and longitude, a 2x2 matrix J taking a
frame displacement (dx, dy) to local true (east, north) metres, its unit
upward normal in true (east, north, up), and its exclusion radius: the
largest ground distance from the centroid to a vertex.

For a projected CRS, J = (1/m) [[cos g, sin g], [-sin g, cos g]], with g the
meridian convergence and m the point scale: a grid direction's true azimuth is
its grid azimuth plus g, and one frame unit is 1/m ground metres. The local
vertical is the ellipsoid normal; the deflection of the vertical (in Norway
some arcseconds) is ignored.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from pyproj import CRS

from tin_engine.crs import frame_factors, parse_crs, reprojector

F64 = npt.NDArray[np.float64]

#: East and north are undefined at the poles; a centroid beyond this latitude
#: is refused.
MAX_LATITUDE_DEG = 89.9
#: A chosen threshold: a triangle whose ground edges, from its first vertex,
#: span a parallelogram smaller than this fraction of the product of their
#: lengths (an angle under about 6e-8 degrees) has zero area and no normal.
ZERO_AREA_SINE = 1e-9


@dataclass(frozen=True)
class ReceiverFrame:
    """One receiver per triangle, in the triangles' order.

    `x`, `y`, `z_obs`: the centroid in the frame, z the triangle's mean height
    plus the observer offset. `lat`, `lon`: the centroid's geodetic degrees.
    `jacobian` (n, 2, 2): (east, north) = J (dx, dy). `normal` (n, 3): unit,
    upward, in true (east, north, up). `exclusion_radius`: ground metres.
    `earth_radius_m`: the Gaussian mean radius sqrt(M N) of the CRS's
    ellipsoid at the receivers' mean latitude, NaN when there are none.
    """

    x: F64
    y: F64
    z_obs: F64
    lat: F64
    lon: F64
    jacobian: F64
    normal: F64
    exclusion_radius: F64
    earth_radius_m: float


def receivers_from_mesh(
    points: npt.ArrayLike,
    triangles: npt.ArrayLike,
    crs: str | CRS,
    *,
    observer_offset: float = 0.0,
) -> ReceiverFrame:
    """The receiver frame of each triangle of a mesh.

    `points` is (n, 3): x and y in the projected `crs`, z in metres;
    `triangles` is (m, 3) vertex indices. `observer_offset` (metres, at least
    0) raises every observer above the triangle's plane. A `ValueError`, naming
    the first offending triangle (0-based) where there is one: a vertex index
    out of range, a vertex that is not finite, a zero-area triangle, a centroid
    beyond 89.9 degrees of latitude; also a negative offset and a geographic
    CRS.
    """
    if not observer_offset >= 0.0:
        raise ValueError(f"the observer offset must be 0 or more metres, not {observer_offset}")
    parsed = parse_crs(crs)
    if not parsed.is_projected:
        raise ValueError(f"{crs} is not a projected CRS; the mesh must be in a projected CRS")
    pts = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    _refuse_first(np.any((tri < 0) | (tri >= len(pts)), axis=1), "has a vertex index out of range")
    v = pts[tri]  # (m, 3 vertices, 3 coordinates)
    _refuse_first(~np.all(np.isfinite(v), axis=(1, 2)), "has a vertex that is not finite")

    # Edges from the first vertex: differences of nearby coordinates are exact,
    # so the centroid's offsets keep the precision the large coordinates lose.
    d = v[:, 1:] - v[:, :1]  # (m, 2, 3)
    to_centroid = d.sum(axis=1) / 3.0
    centroid = v[:, 0] + to_centroid
    lonlat = reprojector(parsed, "EPSG:4326")(centroid[:, :2])
    lon, lat = lonlat[:, 0], lonlat[:, 1]
    _refuse_first(
        ~(np.abs(lat) <= MAX_LATITUDE_DEG),
        f"has its centroid beyond {MAX_LATITUDE_DEG} degrees of latitude, where east and north "
        "are undefined",
    )

    gamma, scale = frame_factors(parsed, lon, lat)
    c, s = np.cos(np.radians(gamma)) / scale, np.sin(np.radians(gamma)) / scale
    jac = np.stack([np.stack([c, s], axis=-1), np.stack([-s, c], axis=-1)], axis=-2)

    ground = np.concatenate([np.einsum("mij,mkj->mki", jac, d[:, :, :2]), d[:, :, 2:]], axis=2)
    n = np.cross(ground[:, 0], ground[:, 1])
    area = np.linalg.norm(n, axis=1)
    lengths = np.linalg.norm(ground, axis=2).prod(axis=1)
    _refuse_first(~(area > ZERO_AREA_SINE * lengths), "has zero area, so no normal")
    n = np.where(n[:, 2:] < 0.0, -n, n) / area[:, None]

    offsets = np.concatenate([-to_centroid[:, None, :2], d[:, :, :2] - to_centroid[:, None, :2]], 1)
    radius = np.linalg.norm(np.einsum("mij,mkj->mki", jac, offsets), axis=2).max(
        axis=1, initial=0.0
    )

    return ReceiverFrame(
        x=centroid[:, 0],
        y=centroid[:, 1],
        z_obs=centroid[:, 2] + observer_offset,
        lat=lat,
        lon=lon,
        jacobian=jac.reshape(-1, 2, 2),
        normal=n.reshape(-1, 3),
        exclusion_radius=radius,
        earth_radius_m=_gaussian_radius(parsed, lat),
    )


def _refuse_first(bad: npt.NDArray[np.bool_], reason: str) -> None:
    """A `ValueError` naming the first triangle `bad` marks, if any."""
    if bad.any():
        raise ValueError(f"triangle {int(np.argmax(bad))} {reason}")


def _gaussian_radius(crs: CRS, lat: F64) -> float:
    """sqrt(M N) on `crs`'s ellipsoid at the mean of `lat`; NaN for no latitudes."""
    ellipsoid = crs.ellipsoid
    if ellipsoid is None or lat.size == 0:
        return math.nan
    a, b = ellipsoid.semi_major_metre, ellipsoid.semi_minor_metre
    e2 = 1.0 - (b / a) ** 2
    w = 1.0 - e2 * math.sin(math.radians(float(lat.mean()))) ** 2
    return math.sqrt(a * (1.0 - e2) / w**1.5 * a / math.sqrt(w))
