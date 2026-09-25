"""The raster view binding and its one Python caller: increment 12, tests 6-10.

`docs/increments/12-dem-to-mesh.md` R4 and R5. `_core.raster_view` builds a
zero-copy view over a C-contiguous float32 or float64 array, keeps the array
alive, and takes the four affine scalars by keyword only. `_core.sample` runs
the batch sampler and returns `(z, valid)`. `tin_engine.raster.to_core` is the
only module that builds one from a `DemTile`.

Tests 7 (the lifetime anchor) and 10 (what crosses the boundary) are
invariant-critical.

Every new symbol is fetched inside a fixture, so a missing one fails its own
tests and leaves the rest of the session collecting.
"""

from __future__ import annotations

import gc
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from geotiff_fixtures import DELTA_X, DELTA_Y, TIE_X, TIE_Y, elevations, micro_tiff
from tin_engine import _core
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import DemTile

SENTINEL = -32767.0
AFFINE: dict[str, float] = {"x_min": TIE_X, "y_max": TIE_Y, "delta_x": DELTA_X, "delta_y": DELTA_Y}

Factory = Callable[..., Any]
Sampler = Callable[[Any, np.ndarray], tuple[np.ndarray, np.ndarray]]


@pytest.fixture
def raster_view() -> Factory:
    factory: Factory = _core.raster_view  # type: ignore[attr-defined]
    return factory


@pytest.fixture
def sample() -> Sampler:
    sampler: Sampler = _core.sample  # type: ignore[attr-defined]
    return sampler


@pytest.fixture
def to_core() -> Callable[[DemTile], Any]:
    from tin_engine.raster import to_core  # type: ignore[import-not-found]

    adapter: Callable[[DemTile], Any] = to_core
    return adapter


def node(row: int, col: int) -> tuple[float, float]:
    """Node (row, col) of the baseline grid, by `RasterGeometry::node`'s expression."""
    return (TIE_X + col * DELTA_X, TIE_Y - row * DELTA_Y)


def points(*xy: tuple[float, float]) -> np.ndarray:
    return np.asarray(xy, dtype=np.float64).reshape(-1, 2)


# Interior nodes only, so every z is exact (the far edges carry rounding, R2).
INTERIOR = [(0, 0), (1, 2), (0, 1), (1, 0)]


class TestAcceptedArrays:
    """Test 6: C-contiguous float32 and float64 are accepted, as is, uncopied."""

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    def test_samples_the_nodes_of_either_float_type(
        self, raster_view: Factory, sample: Sampler, dtype: Any
    ) -> None:
        array = elevations(dtype)
        view = raster_view(array, **AFFINE)
        z, valid = sample(view, points(*(node(r, c) for r, c in INTERIOR)))
        assert valid.all()
        assert_array_equal(z, [10 * r + c for r, c in INTERIOR])

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    def test_never_copies(self, raster_view: Factory, sample: Sampler, dtype: Any) -> None:
        array = elevations(dtype)
        view = raster_view(array, **AFFINE)
        array[1, 2] = 77.0
        z, valid = sample(view, points(node(1, 2)))
        assert valid.tolist() == [True]
        assert z.tolist() == [77.0]

    def test_accepts_a_read_only_array(self, raster_view: Factory, sample: Sampler) -> None:
        array = elevations(np.float32)
        array.flags.writeable = False
        z, _ = sample(raster_view(array, **AFFINE), points(node(1, 1)))
        assert z.tolist() == [11.0]


class TestRefusedArrays:
    """Test 6: anything needing a conversion or a copy is a `TypeError`."""

    @pytest.mark.parametrize(
        "make",
        [
            pytest.param(lambda: np.asfortranarray(elevations(np.float64)), id="fortran"),
            pytest.param(lambda: elevations(np.float64, cols=8)[:, ::2], id="strided"),
            pytest.param(lambda: elevations(np.int16), id="int16"),
            pytest.param(lambda: elevations(np.float64).tolist(), id="list"),
        ],
    )
    def test_refuses(self, raster_view: Factory, make: Callable[[], Any]) -> None:
        candidate = make()
        if isinstance(candidate, np.ndarray):
            # The fixture is what its id says, or the refusal proves nothing.
            assert not candidate.flags.c_contiguous or candidate.dtype == np.int16
        with pytest.raises(TypeError):
            raster_view(candidate, **AFFINE)

    def test_affine_scalars_are_keyword_only(self, raster_view: Factory) -> None:
        """Test 9: a positional call is how the legacy transposed them."""
        with pytest.raises(TypeError):
            raster_view(elevations(np.float64), TIE_X, TIE_Y, DELTA_X, DELTA_Y)


class TestKeepAlive:
    """Test 7 (invariant-critical): the view holds the array, not a dangling pointer."""

    def test_sampling_after_every_python_reference_is_gone(
        self, raster_view: Factory, sample: Sampler
    ) -> None:
        view = raster_view(elevations(np.float64) + 1000.0, **AFFINE)
        gc.collect()
        # Churn the allocator so a freed buffer is likely to be reused.
        _junk = [np.full((3, 4), -1.0) for _ in range(64)]
        z, valid = sample(view, points(*(node(r, c) for r, c in INTERIOR)))
        assert valid.all()
        assert_array_equal(z, [1000.0 + 10 * r + c for r, c in INTERIOR])


class TestSample:
    """Test 8: shapes in and out."""

    def test_returns_z_and_valid_of_shape_n(self, raster_view: Factory, sample: Sampler) -> None:
        view = raster_view(elevations(np.float64), **AFFINE)
        z, valid = sample(view, points(node(0, 0), node(1, 1), (0.0, 0.0)))
        assert z.shape == (3,) and z.dtype == np.float64
        assert valid.shape == (3,) and valid.dtype == np.bool_
        assert valid.tolist() == [True, True, False]

    def test_sentinel_and_nan_corners_are_not_valid(
        self, raster_view: Factory, sample: Sampler
    ) -> None:
        array = elevations(np.float32, rows=3, cols=5)
        array[0, 2] = SENTINEL  # a corner of node (0, 1)'s bilinear cell (0, 1)
        array[2, 4] = np.nan
        view = raster_view(array, **AFFINE, nodata=SENTINEL)
        z, valid = sample(view, points(node(0, 1), node(1, 3), node(0, 3)))
        assert valid.tolist() == [False, False, True]
        assert z[:2].tolist() == [0.0, 0.0]  # invalid z is 0.0, never NaN (R2)
        assert z[2] == 3.0

    def test_empty_points(self, raster_view: Factory, sample: Sampler) -> None:
        z, valid = sample(raster_view(elevations(np.float64), **AFFINE), np.zeros((0, 2)))
        assert z.shape == (0,) and valid.shape == (0,)

    @pytest.mark.parametrize("shape", [(3,), (3, 3), (3, 1), (2, 2, 2)])
    def test_wrong_shaped_points_are_a_value_error(
        self, raster_view: Factory, sample: Sampler, shape: tuple[int, ...]
    ) -> None:
        view = raster_view(elevations(np.float64), **AFFINE)
        with pytest.raises(ValueError):
            sample(view, np.zeros(shape))


class TestToCore:
    """Test 10 (invariant-critical): exactly the array, four scalars and the sentinel cross."""

    def test_samples_a_decoded_tile_at_its_nodes(self, to_core: Any, sample: Sampler) -> None:
        tile = decode_dem(micro_tiff())
        z, valid = sample(to_core(tile), points(*(node(r, c) for r, c in INTERIOR)))
        assert valid.all()
        assert_array_equal(z, [10 * r + c for r, c in INTERIOR])

    def test_the_deltas_are_not_swapped_and_the_origin_is_the_top(
        self, to_core: Any, sample: Sampler
    ) -> None:
        """The baseline's deltas differ (10 and 5) and its values are distinct, so
        a swapped delta, `y_min` for `y_max`, or a transposed index each move
        these off-node points onto different values."""
        tile = decode_dem(micro_tiff())
        x, y = node(0, 0)
        query = points((x + 5.0, y - 2.5), (x + 25.0, y - 7.5))
        z, valid = sample(to_core(tile), query)
        assert valid.all()
        assert z.tolist() == pytest.approx([5.5, 17.5], rel=1e-12)

    def test_forwards_the_sentinel(self, to_core: Any, sample: Sampler) -> None:
        array = elevations(np.float32)
        array[0, 2] = SENTINEL  # a corner of node (0, 1)'s bilinear cell (0, 1)
        tile = decode_dem(micro_tiff(array, nodata="-32767"))
        _, valid = sample(to_core(tile), points(node(0, 1), node(1, 3)))
        assert valid.tolist() == [False, True]

    def test_uses_the_meta_corner_for_an_area_registered_file(
        self, to_core: Any, sample: Sampler
    ) -> None:
        from geotiff_fixtures import GT_RASTER_TYPE, PIXEL_IS_AREA, with_keys

        tile = decode_dem(micro_tiff(geokeys=with_keys({GT_RASTER_TYPE: PIXEL_IS_AREA})))
        meta = tile.meta
        assert meta.x_min != TIE_X  # shifted half a cell by decode_dem
        query = points((meta.x_min + DELTA_X, meta.y_max - DELTA_Y))
        z, valid = sample(to_core(tile), query)
        assert valid.tolist() == [True]
        assert z.tolist() == [11.0]

    def test_the_view_outlives_the_tile(self, to_core: Any, sample: Sampler) -> None:
        view = to_core(decode_dem(micro_tiff()))
        gc.collect()
        z, _ = sample(view, points(node(1, 1)))
        assert z.tolist() == [11.0]

    def test_refuses_a_writeable_array(self, to_core: Any) -> None:
        """R5: `to_core` checks the array is read-only. A validated `DemTile`
        cannot hold a writeable one, so `model_construct` bypasses validation to
        build the case. The exception type is not in the design; ValueError is
        this suite's assumption."""
        tile = decode_dem(micro_tiff())
        writeable = np.array(tile.array, copy=True)
        with pytest.raises(ValueError, match="read-only"):
            to_core(DemTile.model_construct(meta=tile.meta, array=writeable))

    def test_raster_view_takes_no_crs_or_shape(self, raster_view: Factory) -> None:
        """No CRS, EPSG, rows or cols crosses (R4). pybind11 functions have no
        `inspect.signature`, so the bound docstring's first line is the signature."""
        signature = (raster_view.__doc__ or "").splitlines()[0]
        assert signature.startswith("raster_view(")
        for name in ("x_min", "y_max", "delta_x", "delta_y", "nodata"):
            assert f"{name}:" in signature
        for forbidden in ("epsg", "crs", "rows", "cols"):
            assert forbidden not in signature.lower()

    def test_shape_comes_from_the_array(self, raster_view: Factory) -> None:
        with pytest.raises(TypeError):
            raster_view(elevations(np.float64), **AFFINE, rows=3, cols=4)
