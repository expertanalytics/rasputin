"""`raster.core_args`: exactly what crosses into C++, checked once, for either module.

`docs/increments/51-open-core.md`, *4.3 The raster across two extension
modules* and *9.1 Public (P1 unless marked)*, T8. ``core_args`` keeps
``to_core``'s refusals word for word, and ``to_core`` becomes
``raster_view(**core_args(tile))``, so the view built either way samples
alike. ``core_args`` is looked up inside each test, so each fails on its own.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np
import pytest

from mosaic_fixtures import DX, DY, X0, Y0, whole
from tin_engine import raster
from tin_engine._core import raster_view, sample
from tin_engine.io.models import DemTile, RasterMeta

#: The fields of ``RasterArgs``, in ``raster_view``'s keyword names.
CROSSING = ("array", "x_min", "y_max", "delta_x", "delta_y", "nodata")


def core_args(tile: DemTile) -> Any:
    return raster.core_args(tile)  # type: ignore[attr-defined]


@pytest.fixture
def tile() -> DemTile:
    """3 x 3 float32, point-registered, dx 10 and dy 5, node (r, c) = 100 r + c + 0.5,
    one NoData node at (1, 2)."""
    values = whole(rows=3, cols=3).array.copy()
    values[1, 2] = -9999.0
    return whole(rows=3, cols=3, array=values, nodata=-9999.0)


def nine_points() -> np.ndarray:
    """Three nodes, three points inside cells, three on cell sides: every
    bilinear path, the NoData node's cells included."""
    xs = X0 + DX * np.array([0.0, 1.0, 2.0, 0.25, 1.5, 0.75, 1.0, 0.5, 2.0])
    ys = Y0 - DY * np.array([0.0, 1.0, 2.0, 0.5, 0.25, 1.75, 0.5, 2.0, 1.5])
    return np.column_stack([xs, ys])


def test_the_record_is_exactly_what_crosses(tile: DemTile) -> None:
    args = core_args(tile)
    assert dataclasses.is_dataclass(args)
    fields = dataclasses.asdict(args)
    assert tuple(fields) == CROSSING
    assert np.shares_memory(fields["array"], tile.array)
    assert not fields["array"].flags.writeable
    meta = tile.meta
    assert (fields["x_min"], fields["y_max"], fields["delta_x"], fields["delta_y"]) == (
        meta.x_min,
        meta.y_max,
        meta.delta_x,
        meta.delta_y,
    )
    assert fields["nodata"] == -9999.0


def test_the_record_is_frozen(tile: DemTile) -> None:
    args = core_args(tile)
    with pytest.raises(dataclasses.FrozenInstanceError):
        args.x_min = 0.0


def test_a_writeable_array_is_refused_with_todays_words(tile: DemTile) -> None:
    # model_construct skips DemTile's validator, which would make the array read-only.
    writeable = DemTile.model_construct(meta=tile.meta, array=tile.array.copy())
    assert writeable.array.flags.writeable
    with pytest.raises(ValueError, match="needs a read-only array; build the tile through DemTile"):
        core_args(writeable)


def test_a_geographic_tile_is_refused_with_todays_words() -> None:
    meta = RasterMeta(
        x_min=10.0,
        y_max=60.0,
        delta_x=0.001,
        delta_y=0.0005,
        cols=3,
        rows=3,
        epsg=4326,
        geographic=True,
        nodata=None,
        nodata_source="absent",
        pixel_is_area=False,
        vertical_unit_assumed=True,
    )
    tile = DemTile(meta=meta, array=np.zeros((3, 3), dtype=np.float32))
    expected = r"a geographic DEM \(EPSG:4326, degrees\) cannot enter the core"
    with pytest.raises(ValueError, match=expected):
        core_args(tile)


def test_to_core_and_raster_view_of_core_args_sample_alike(tile: DemTile) -> None:
    points = nine_points()
    z_core, ok_core = sample(raster.to_core(tile), points)
    z_args, ok_args = sample(raster_view(**dataclasses.asdict(core_args(tile))), points)
    np.testing.assert_array_equal(ok_args, ok_core)
    np.testing.assert_array_equal(z_args, z_core)  # one kernel, one input: bit for bit
    # The probe can fail: the NoData node refuses some points and not others.
    assert ok_core.any() and not ok_core.all()
