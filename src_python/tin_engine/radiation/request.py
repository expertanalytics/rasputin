"""What a radiation geometry is built with: a frozen request, data only.

`n_directions` horizon directions per receiver, a power of two so that the
directions halve and double cleanly; `refraction_k`, the refraction
coefficient that flattens the Earth's curvature for a ray; `d_near_m` and
`d_far_m`, how far the horizon is searched on the full-resolution DEM and on
the coarse one beyond it; `limb`, whether the sun counts as risen at its
centre or at its upper edge.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator


class RadiationRequest(BaseModel):
    """The settings of a radiation geometry, each with a default."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    n_directions: int = 64
    refraction_k: float = Field(0.13, ge=0.0, lt=1.0, allow_inf_nan=False)
    d_near_m: float = Field(20_000.0, gt=0.0, allow_inf_nan=False)
    d_far_m: float = Field(200_000.0, gt=0.0, allow_inf_nan=False)
    limb: Literal["centre", "upper"] = "centre"

    @field_validator("n_directions")
    @classmethod
    def _power_of_two(cls, n: int) -> int:
        if n < 4 or n & (n - 1):
            raise ValueError(f"must be a power of two, at least 4, not {n}")
        return n

    @field_validator("d_far_m")
    @classmethod
    def _beyond_near(cls, d_far: float, info: ValidationInfo) -> float:
        d_near = info.data.get("d_near_m")
        if d_near is not None and d_far < d_near:
            raise ValueError(f"must be at least d_near_m ({d_near} m), not {d_far} m")
        return d_far
