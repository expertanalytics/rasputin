"""A synthetic mountain to run the snow on: a five-ridged massif, a second
hill and random terrain noise, on a square grid.

The heights are chosen, not measured. Every aspect occurs at every height,
and the noise (a random-phase surface with power spectrum k^-3, no
wavelengths under 40 m, 25 m RMS) gives a mesh the triangle density of
real terrain. The construction is `docs/increments/46-hbv-snow.md`, *4.5 The synthetic case*.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict, Field


class TerrainSpec(BaseModel):
    """The grid: `size` by `size` nodes `step` metres apart, node (0, 0) the
    north-west corner at (`x0`, `y0`) in EPSG `epsg`; `seed` for the noise."""

    model_config = ConfigDict(frozen=True)

    size: int = Field(default=1201, ge=2)
    step: float = Field(default=10.0, gt=0.0)
    x0: float = 494000.0
    y0: float = 6806000.0
    epsg: int = 25833
    seed: int = 46


def mountain(spec: TerrainSpec) -> npt.NDArray[np.float32]:
    """Heights in metres, float32, row 0 the northernmost."""
    n, step, half = spec.size, spec.step, (spec.size - 1) * spec.step / 2
    x, y = np.meshgrid((np.arange(n) * step) - half, half - (np.arange(n) * step))
    r, phi = np.hypot(x, y), np.arctan2(y, x)
    z = 300 + 1600 * np.exp(-((r / 2500) ** 2)) * (
        1 + 0.15 * (1 - np.exp(-((r / 1500) ** 2))) * np.cos(5 * phi)
    )
    z += 250 * np.exp(-((x - 3500) ** 2 + (y + 3000) ** 2) / 900**2)
    rng = np.random.default_rng(spec.seed)
    k = np.hypot(*np.meshgrid(np.fft.fftfreq(n, step), np.fft.fftfreq(n, step)))
    amp = np.zeros_like(k)
    keep = (k > 0) & (k <= 1 / 40.0)
    amp[keep] = k[keep] ** -1.5
    noise = np.real(np.fft.ifft2(amp * np.exp(2j * np.pi * rng.random(k.shape))))
    z += noise * (25.0 / noise.std())
    out: npt.NDArray[np.float32] = z.astype(np.float32)
    return out
