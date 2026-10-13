"""The synthetic mountain and its GeoTIFF: `snow/terrain.py`, `io/geotiff_write.py`, test R5.

`docs/increments/46-hbv-snow.md`, *4.5 The synthetic case*: a 12 km square,
10 m grid (1201 x 1201 nodes), EPSG:25833, node (0, 0) at (494000, 6806000),
a five-ridged massif, a second hill and a k^-3 random-phase noise from
`numpy.random.default_rng(46)`, built exactly as the design probe builds it
(`46-probes/probe_spring.py`), which `terrain.py` reproduces bit for bit. The
probe's construction is restated below as the oracle (the suite reads no file
under `docs/`). Section 7's R5: `encode_geotiff` read back by `io/geotiff.py`
gives shape 1201 x 1201, step 10, EPSG:25833 and z from 234.0 to 1924.4 within
0.1 m.

Interface assumed (names chosen here and listed in the handback):

- `snow.terrain.TerrainSpec`, frozen, defaults `size` 1201, `step` 10.0,
  `x0` 494000.0, `y0` 6806000.0 (node (0, 0), the north-west corner), `epsg`
  25833, `seed` 46.
- `snow.terrain.mountain(spec) -> float32 array (size, size)`, row 0 north.
- `io.geotiff_write.encode_geotiff(grid, x0, y0, step, epsg) -> bytes`:
  point-registered (PixelIsPoint), float32, node (0, 0) at (x0, y0).
- `rasputin synthetic-terrain --out PATH` writes the default terrain's GeoTIFF.
"""

from __future__ import annotations

import importlib
import io
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_array_equal

from cli_driver import ran
from tin_engine.io.geotiff import decode_dem


def probe_terrain() -> npt.NDArray[np.float32]:
    """Section 4.5's height field, as the design probe computes it."""
    n_side, step = 1201, 10.0
    xs = (np.arange(n_side) * step) - 6000.0
    ys = 6000.0 - (np.arange(n_side) * step)
    x, y = np.meshgrid(xs, ys)
    r, phi = np.hypot(x, y), np.arctan2(y, x)
    z = 300 + 1600 * np.exp(-((r / 2500) ** 2)) * (
        1 + 0.15 * (1 - np.exp(-((r / 1500) ** 2))) * np.cos(5 * phi)
    )
    z += 250 * np.exp(-((x - 3500) ** 2 + (y + 3000) ** 2) / 900**2)
    rng = np.random.default_rng(46)
    k = np.hypot(*np.meshgrid(np.fft.fftfreq(n_side, step), np.fft.fftfreq(n_side, step)))
    amp = np.zeros_like(k)
    keep = (k > 0) & (k <= 1 / 40.0)
    amp[keep] = k[keep] ** -1.5
    noise = np.real(np.fft.ifft2(amp * np.exp(2j * np.pi * rng.random(k.shape))))
    z += noise * (25.0 / noise.std())
    return z.astype(np.float32)


@pytest.fixture(scope="module")
def terrain() -> ModuleType:
    return importlib.import_module("tin_engine.snow.terrain")


@pytest.fixture(scope="module")
def grid(terrain: ModuleType) -> npt.NDArray[np.float32]:
    out: npt.NDArray[np.float32] = terrain.mountain(terrain.TerrainSpec())
    return out


@pytest.fixture(scope="module")
def encode() -> Any:
    return importlib.import_module("tin_engine.io.geotiff_write").encode_geotiff


class TestMountain:
    def test_the_spec_defaults_are_section_4_5s(self, terrain: ModuleType) -> None:
        spec = terrain.TerrainSpec()
        assert (spec.size, spec.step, spec.x0, spec.y0, spec.epsg, spec.seed) == (
            1201,
            10.0,
            494000.0,
            6806000.0,
            25833,
            46,
        )

    def test_it_is_the_probes_terrain_bit_for_bit(self, grid: npt.NDArray[np.float32]) -> None:
        assert grid.dtype == np.float32
        assert_array_equal(grid, probe_terrain())

    def test_its_heights_span_234_0_to_1924_4_m(self, grid: npt.NDArray[np.float32]) -> None:
        assert float(grid.min()) == pytest.approx(234.0, abs=0.05)
        assert float(grid.max()) == pytest.approx(1924.4, abs=0.05)

    def test_row_0_is_north(self, grid: npt.NDArray[np.float32]) -> None:
        # The second hill sits 3 km south of the centre's row: its rows are below 600.
        hill_row, hill_col = 600 + 300, 600 + 350
        assert grid[hill_row, hill_col] > grid[600 - 300, hill_col] + 100

    def test_it_is_the_same_on_every_call(
        self, terrain: ModuleType, grid: npt.NDArray[np.float32]
    ) -> None:
        assert_array_equal(terrain.mountain(terrain.TerrainSpec()), grid)


class TestGeoTiff:
    """R5: `encode_geotiff` read back through the reader `rasputin mesh` uses."""

    def test_the_default_terrain_reads_back_as_section_4_5_says(
        self, encode: Any, grid: npt.NDArray[np.float32]
    ) -> None:
        data = encode(grid, 494000.0, 6806000.0, 10.0, 25833)
        tile = decode_dem(io.BytesIO(data))
        meta = tile.meta
        assert (meta.rows, meta.cols) == (1201, 1201)
        assert (meta.delta_x, meta.delta_y) == (10.0, 10.0)
        assert meta.epsg == 25833
        assert (meta.x_min, meta.y_max) == (494000.0, 6806000.0)
        assert meta.pixel_is_area is False
        assert float(tile.array.min()) == pytest.approx(234.0, abs=0.1)
        assert float(tile.array.max()) == pytest.approx(1924.4, abs=0.1)
        assert_array_equal(tile.array, grid)

    def test_a_float64_grid_is_stored_as_float32(self, encode: Any) -> None:
        grid = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0 + 1e-9]])
        tile = decode_dem(io.BytesIO(encode(grid, 500000.0, 6600000.0, 5.0, 25833)))
        assert_array_equal(tile.array, grid.astype(np.float32))
        assert (tile.meta.rows, tile.meta.cols, tile.meta.delta_x) == (2, 3, 5.0)


class TestCli:
    def test_synthetic_terrain_writes_the_default_terrain(
        self, tmp_path: Path, grid: npt.NDArray[np.float32]
    ) -> None:
        target = tmp_path / "synth.tif"
        ran("synthetic-terrain", "--out", str(target))
        with target.open("rb") as f:
            tile = decode_dem(f)
        assert tile.meta.epsg == 25833
        assert (tile.meta.x_min, tile.meta.y_max) == (494000.0, 6806000.0)
        assert_array_equal(tile.array, grid)
