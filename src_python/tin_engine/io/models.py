"""What a decoded GeoTIFF becomes: `RasterMeta`, `DemTile`, and the one error.

Increment 11, `docs/increments/11-raster-ingestion.md` §7. Both models are
frozen Pydantic V2. `RasterMeta` is a superset of the `_core` boundary
contract: `epsg`, `nodata_source`, `pixel_is_area` and `vertical_unit_assumed`
are provenance and must not be forwarded to `_core` (increment 12's adapter
owns that). This module imports nothing first-party and never `_core`.
"""

from __future__ import annotations

from typing import Any, Literal, Self

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict, Field, StrictInt, field_validator, model_validator


class GeoTiffError(ValueError):
    """A file the reader refuses (§5).

    The message names the tag or GeoKey by number and by name, and the file's
    value. One type on purpose: no caller branches on which refusal fired.
    """


class RasterMeta(BaseModel):
    """The node grid of one tile, in metres of a projected CRS.

    `x_min`/`y_max` are the upper-left *node*, already shifted inward half a
    cell for an area-registered file (§4). `rows`/`cols` are `StrictInt`
    (B1, §14), so the legacy's `(12.0, 16.0)` is refused at the type.
    """

    model_config = ConfigDict(frozen=True)

    x_min: float
    y_max: float
    delta_x: float
    delta_y: float
    cols: StrictInt
    rows: StrictInt
    epsg: int
    # Finite or None (§6, round 2): a NaN sentinel matches no cell under
    # `==`, and a NaN field would also break model equality.
    nodata: float | None = Field(allow_inf_nan=False)
    nodata_source: Literal["tag", "caller", "absent"]
    pixel_is_area: bool
    vertical_unit_assumed: bool

    @model_validator(mode="after")
    def _absent_means_none(self) -> Self:
        # One direction only: a NaN declared by tag or caller is also None.
        if self.nodata_source == "absent" and self.nodata is not None:
            raise ValueError("nodata_source 'absent' requires nodata None")
        return self


class DemTile(BaseModel):
    """One decoded tile: its metadata and a read-only, C-contiguous 2-D array.

    A disagreement between `meta` and `array.shape` is Pydantic's
    `ValidationError`, not `GeoTiffError`: it is a programming error in
    whoever built the model, not a fact about a file (§7, B1).
    """

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    meta: RasterMeta
    array: npt.NDArray[Any]

    @field_validator("array")
    @classmethod
    def _read_only_c_contiguous(cls, value: npt.NDArray[Any]) -> npt.NDArray[Any]:
        if value.ndim != 2 or value.dtype not in (np.float32, np.float64):
            raise ValueError(
                f"need a 2-D float32 or float64 array, got {value.dtype} {value.shape}"
            )
        # A read-only *view*: the caller's array keeps its own flags, and the
        # tile never hands out a writeable handle (§7, the concurrency rule).
        view = np.ascontiguousarray(value).view()
        view.flags.writeable = False
        return view

    @model_validator(mode="after")
    def _shape_agrees(self) -> Self:
        if (self.meta.rows, self.meta.cols) != self.array.shape:
            raise ValueError(
                f"meta (rows, cols) = {(self.meta.rows, self.meta.cols)} "
                f"but array.shape = {self.array.shape}"
            )
        return self
