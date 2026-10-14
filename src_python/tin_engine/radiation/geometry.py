"""A radiation geometry: what the evaluator reads for each receiver.

Per receiver, its geodetic latitude and longitude, its unit normal in true
(east, north, up), its horizon in `n_directions` true azimuths (int16 in 0.01
degrees, direction k at k 360 / n_directions) and how far each direction's
search reached on the DEM (uint8 whole km, 255 when it never left it), with
the request it was built by. Positions in the projected frame are left out:
the evaluator never needs them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from tin_engine import _core
from tin_engine.radiation.request import RadiationRequest

if TYPE_CHECKING:
    from tin_engine.radiation.frame import ReceiverFrame

F64 = npt.NDArray[np.float64]

#: The horizon that hides nothing: -90 degrees in 0.01 degrees.
OPEN_HORIZON = -9000


@dataclass(frozen=True)
class RadiationGeometry:
    """`lat`, `lon` (n,) degrees; `normal` (n, 3); `horizon` (n, N) int16;
    `reach` (n, N) uint8; `n_directions` N; `request`."""

    lat: F64
    lon: F64
    normal: F64
    horizon: npt.NDArray[np.int16]
    reach: npt.NDArray[np.uint8]
    n_directions: int
    request: RadiationRequest

    @classmethod
    def open_sky(cls, frame: ReceiverFrame, request: RadiationRequest) -> RadiationGeometry:
        """`frame`'s receivers with no terrain around them: a horizon of -90
        degrees and a reach of 0 in every direction, so only a receiver's own
        plane can shade it."""
        shape = (len(frame.lat), request.n_directions)
        return cls(
            lat=frame.lat, lon=frame.lon, normal=frame.normal,
            horizon=np.full(shape, OPEN_HORIZON, dtype=np.int16),
            reach=np.zeros(shape, dtype=np.uint8),
            n_directions=request.n_directions, request=request,
        )  # fmt: skip

    @classmethod
    def searched(
        cls, frame: ReceiverFrame, near: _core.RasterView, request: RadiationRequest, threads: int
    ) -> RadiationGeometry:
        """`frame`'s receivers with their horizons searched on the float32
        `near` raster to `request.d_near_m` and no further, with the Earth's
        curvature less refraction at the frame's radius."""
        horizon, reach = _core.horizon_search(
            frame.x, frame.y, frame.z_obs, frame.lat, frame.lon, frame.jacobian, frame.normal,
            frame.exclusion_radius, near=near, far=None, n_directions=request.n_directions,
            curvature=(1 - request.refraction_k) / (2 * frame.earth_radius_m),
            d_near=request.d_near_m, d_far=request.d_near_m, threads=threads,
        )  # fmt: skip
        return cls(
            lat=frame.lat, lon=frame.lon, normal=frame.normal, horizon=horizon, reach=reach,
            n_directions=request.n_directions, request=request,
        )  # fmt: skip
