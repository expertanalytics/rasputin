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


class IndexWindow(BaseModel):
    """`rows x cols` nodes starting at `(row0, col0)` (moved from `mosaic.py` by
    23a-1, because `io/` imports nothing first-party but this module)."""

    model_config = ConfigDict(frozen=True)

    row0: int
    col0: int
    rows: int
    cols: int


class RasterMeta(BaseModel):
    """The node grid of one tile, in metres of a projected CRS, or in degrees
    when `geographic` (15c-2, D6). `crs` is `EPSG:n` when `epsg` is
    given; a CRS without an EPSG code has `epsg` None and its text in `crs`.

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
    epsg: int | None
    crs: str = ""
    geographic: bool = False
    # Finite or None (§6, round 2): a NaN sentinel matches no cell under
    # `==`, and a NaN field would also break model equality.
    nodata: float | None = Field(allow_inf_nan=False)
    nodata_source: Literal["tag", "caller", "absent"]
    pixel_is_area: bool
    vertical_unit_assumed: bool

    @model_validator(mode="before")
    @classmethod
    def _crs_from_epsg(cls, data: Any) -> Any:
        # An EPSG code is the CRS; text stands alone only without one.
        if isinstance(data, dict) and data.get("epsg") is not None:
            data = {**data, "crs": f"EPSG:{data['epsg']}"}
        elif isinstance(data, dict) and not data.get("crs"):
            raise ValueError("a RasterMeta needs an EPSG code or a CRS")
        return data

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
        # An owned read-only copy, stored as a view of it: the caller's buffer
        # is never shared, and the view's flag cannot be set back because its
        # base is read-only too (§7, the concurrency rule, round 3).
        copy = np.array(value, order="C", copy=True)
        copy.flags.writeable = False
        return copy.view()

    @classmethod
    def _adopt(cls, meta: RasterMeta, array: npt.NDArray[Any]) -> DemTile:
        """Take `array` as the tile's own, without the copy (increment 15a, R7).

        Private: its two callers are `mosaic.assemble` and
        `target_grid.resample` (increment 15e), each on a canvas it allocated
        and never hands out writable, so each peaks at one canvas, not two.
        The same checks as the public constructor, plus C-contiguity, which
        the constructor gets from its copy. The **passed** buffer is set
        read-only, so the concurrency rule (§7) holds: nobody keeps a writable
        reference.
        """
        if array.ndim != 2 or array.dtype not in (np.float32, np.float64):
            raise ValueError(
                f"need a 2-D float32 or float64 array, got {array.dtype} {array.shape}"
            )
        if not array.flags.c_contiguous or array.shape != (meta.rows, meta.cols):
            raise ValueError(
                f"need a C-contiguous array of shape {(meta.rows, meta.cols)}, "
                f"got {array.shape}, C-contiguous {array.flags.c_contiguous}"
            )
        array.flags.writeable = False
        return cls.model_construct(meta=meta, array=array)

    @model_validator(mode="after")
    def _shape_agrees(self) -> Self:
        if (self.meta.rows, self.meta.cols) != self.array.shape:
            raise ValueError(
                f"meta (rows, cols) = {(self.meta.rows, self.meta.cols)} "
                f"but array.shape = {self.array.shape}"
            )
        return self
