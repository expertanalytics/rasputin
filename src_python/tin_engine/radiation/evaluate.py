"""The instant evaluator on a radiation geometry, called synchronously.

Per sample time and receiver: lit, cos_i, the beam factor b = lit max(0,
cos_i) and the sine of the apparent sun elevation, from the core. The limb
comes from the geometry's request: the sun's centre, or its upper edge.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from tin_engine import _core
from tin_engine.radiation.geometry import RadiationGeometry

#: How far the sun's upper edge stands above its centre, degrees.
UPPER_LIMB_DEG = 0.2667


def evaluate_instant(
    geometry: RadiationGeometry,
    times_utc_ms: npt.ArrayLike,
    *,
    pressure_hpa: float = 1010.0,
    temperature_c: float = 10.0,
    threads: int = 0,
) -> _core.InstantFactors:
    """The factors at each UTC time (int64 milliseconds), each array (times,
    receivers); refraction at `pressure_hpa` and `temperature_c`, `threads`
    threads (0: every core). Releases the GIL while it runs."""
    return _core.evaluate_instant(
        geometry.lat, geometry.lon, geometry.normal, geometry.horizon,
        np.asarray(times_utc_ms, dtype=np.int64),
        limb_deg=UPPER_LIMB_DEG if geometry.request.limb == "upper" else 0.0,
        pressure_hpa=pressure_hpa, temperature_c=temperature_c,
        n_directions=geometry.n_directions, threads=threads,
    )  # fmt: skip
