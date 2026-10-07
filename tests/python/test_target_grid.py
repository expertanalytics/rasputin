"""G2's gate, G3 the target grid, G4 resampling, G5 check points (increment 15c-2).

`docs/increments/15c-geographic-dem.md`, D2-D4, D6, J3, J6 and "Tests for
@tester" G2-G5. `tin_engine/target_grid.py`:

- `TargetGrid(crs, spacing, row0, col0, rows, cols)`, frozen; node `(R, K)`
  at `(K h, -R h)` in the target CRS (J6); `spacing` whole metres, 1 to 1000.
- `target_grid_for(domain, target, spacing) -> tuple[TargetGrid, Polygon]`:
  the domain grown by the cell diagonal with mitred corners, snapped outward
  to the lattice, and that grown domain (increment 30d, section 3.1: the
  caller no longer grows it a second time).
- `TileWindows(tile)`, a `SourceWindows`: `.meta` and `.window(r0, r1, c0, c1)`.
- `resample(grid, source, threads) -> DemTile`: bilinear in the source's
  index space; a node whose stencil touches NoData or leaves the source is
  NoData (the source's sentinel, or NaN with `nodata=None`).
- `check_point_blocks(grid, source, domain, threads)`: `(xy float64 (N, 2),
  z float32 (N,))` blocks, one per 512 x 512 source block.

PINNED HERE, where D2-D4 are silent:
- `default_spacing(meta, at) -> int`: the source's north-south node spacing
  in metres at `at` (the domain's centroid, in the DEM's CRS), rounded to
  whole metres; for a projected source its `delta_y`, rounded.
- `resample(..., block_rows=256)`: the block of target rows D3 names, as a
  keyword so G4 can vary it.
- `target_grid_for` and `check_point_blocks` take the domain as a
  `DomainPolygon` already in the target CRS (D1 moves it there first).
- The resampled tile's meta: `crs` the target's, `geographic` false,
  `x_min = col0 h`, `y_max = -row0 h`, `delta_x = delta_y = h`.

The oracles are pyproj's own `Transformer` (never `crs.reprojector`), the
source's index arithmetic written out, and `shapely.contains_xy`.

HOW THIS FILE GOES RED: `tin_engine.target_grid` does not exist; it is
imported in a fixture, so each test fails on its own. The `to_core` gate
fails on its assertion (no refusal), because `RasterMeta` ignores
`geographic` before 15c-2.

AMENDED for 30d (`docs/increments/30d-outline-buffer-speed.md`, section 3.1):
the two G3 tests unpack `(grid, grown)`. Red until 30d: `target_grid_for`
returns the grid alone.
"""

from __future__ import annotations

import importlib
import math
import threading
import time
from collections.abc import Callable
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
import shapely
from pydantic import ValidationError
from shapely.geometry import Polygon

from geographic_fixtures import ANADEM_STEP, GLO30_STEP, LAT0, LON0, NODATA, project, rough
from tin_engine.domain import DomainPolygon
from tin_engine.io.models import DemTile, RasterMeta

TARGET = "EPSG:31983"  # SIRGAS 2000 / UTM 23S
H = 30


@pytest.fixture(scope="module")
def tg() -> ModuleType:
    return importlib.import_module("tin_engine.target_grid")


def geographic_tile(
    array: npt.NDArray[Any], *, step: float = ANADEM_STEP, nodata: float | None = None
) -> DemTile:
    meta = RasterMeta(
        x_min=LON0,
        y_max=LAT0,
        delta_x=step,
        delta_y=step,
        cols=array.shape[1],
        rows=array.shape[0],
        epsg=4326,
        crs="EPSG:4326",  # type: ignore[call-arg]
        geographic=True,  # type: ignore[call-arg]
        nodata=nodata,
        nodata_source="absent" if nodata is None else "tag",
        pixel_is_area=False,
        vertical_unit_assumed=True,
    )
    return DemTile(meta=meta, array=array)


def source_nodes(meta: RasterMeta) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    r, c = np.indices((meta.rows, meta.cols), dtype=np.float64)
    return (meta.x_min + c * meta.delta_x).ravel(), (meta.y_max - r * meta.delta_y).ravel()


def covering(tg: ModuleType, meta: RasterMeta, margin: int) -> Any:
    """A target grid over the source's image in TARGET, `margin` nodes past it."""
    lon, lat = source_nodes(meta)
    xy = project("EPSG:4326", TARGET, lon, lat)
    col0 = math.floor(xy[:, 0].min() / H) - margin
    col1 = math.ceil(xy[:, 0].max() / H) + margin
    row0 = math.floor(-xy[:, 1].max() / H) - margin
    row1 = math.ceil(-xy[:, 1].min() / H) + margin
    return tg.TargetGrid(
        crs=TARGET, spacing=H, row0=row0, col0=col0, rows=row1 - row0 + 1, cols=col1 - col0 + 1
    )


def fractional_index(grid: Any, meta: RasterMeta) -> tuple[npt.NDArray[np.float64], ...]:
    """Each target node's fractional (row, col) in the source, by pyproj alone."""
    r, c = np.indices((grid.rows, grid.cols), dtype=np.float64)
    x = (grid.col0 + c) * grid.spacing
    y = -(grid.row0 + r) * grid.spacing
    lonlat = project(TARGET, "EPSG:4326", x.ravel(), y.ravel())
    col = (lonlat[:, 0] - meta.x_min) / meta.delta_x
    row = (meta.y_max - lonlat[:, 1]) / meta.delta_y
    return row.reshape(r.shape), col.reshape(r.shape)


def a_domain(meta: RasterMeta, ring_rc: list[tuple[float, float]]) -> DomainPolygon:
    """A polygon given in the source's (row, col) space, moved into TARGET."""
    lon = np.array([meta.x_min + c * meta.delta_x for _, c in ring_rc])
    lat = np.array([meta.y_max - r * meta.delta_y for r, _ in ring_rc])
    xy = project("EPSG:4326", TARGET, lon, lat)
    return DomainPolygon(polygon=shapely.geometry.polygon.orient(Polygon(xy)), crs=TARGET)


# ------------------------------------------------------------------ G2, the gate


def test_to_core_refuses_a_geographic_tile() -> None:
    """J4: no degrees in `_core`. `raster.to_core` is the one adapter (increment 12)."""
    from tin_engine.raster import to_core

    with pytest.raises(ValueError, match=r"(?i)geographic"):
        to_core(geographic_tile(rough(4, 5)))


# ------------------------------------------------------------------ G3, the grid


class TestTheTargetGrid:
    @pytest.fixture
    def domain(self) -> DomainPolygon:
        meta = geographic_tile(rough(80, 80)).meta
        return a_domain(meta, [(10.3, 12.7), (14.1, 61.9), (66.6, 58.2), (60.2, 9.4)])

    def test_nodes_are_multiples_of_h_bit_for_bit(
        self, tg: ModuleType, domain: DomainPolygon
    ) -> None:
        grid, _ = tg.target_grid_for(domain, TARGET, H)
        assert (grid.crs, grid.spacing) == (TARGET, H)
        source = tg.TileWindows(geographic_tile(rough(80, 80)))
        tile = tg.resample(grid, source, threads=2)
        m = tile.meta
        assert (m.x_min, m.y_max) == (grid.col0 * float(H), -grid.row0 * float(H))
        assert (m.delta_x, m.delta_y, m.rows, m.cols) == (H, H, grid.rows, grid.cols)
        assert (m.crs, m.geographic, m.epsg) == (TARGET, False, 31983)

    def test_the_extent_covers_the_grown_domain_and_no_more(
        self, tg: ModuleType, domain: DomainPolygon
    ) -> None:
        result = tg.target_grid_for(domain, TARGET, H)
        assert isinstance(result, tuple) and len(result) == 2, type(result)
        grid, grown = result
        # 30d, section 3.1: the grown domain it returns is GEOS's mitred buffer
        # (a 4-vertex domain is off the gate of section 3.3), bit for bit.
        expected = domain.polygon.buffer(math.sqrt(2) * H, join_style="mitre")
        assert shapely.to_wkb(grown) == shapely.to_wkb(expected)
        x0, y1 = grid.col0 * H, -grid.row0 * H
        x1, y0 = x0 + (grid.cols - 1) * H, y1 - (grid.rows - 1) * H
        gx0, gy0, gx1, gy1 = grown.bounds
        assert x0 <= gx0 and y0 <= gy0 and x1 >= gx1 and y1 >= gy1
        # Snapped outward by less than a cell on every side.
        assert x0 > gx0 - H and y0 > gy0 - H and x1 < gx1 + H and y1 < gy1 + H

    @pytest.mark.parametrize("spacing", [0, 1001, -30])
    def test_spacing_is_whole_metres_from_1_to_1000(self, tg: ModuleType, spacing: int) -> None:
        with pytest.raises(ValidationError):
            tg.TargetGrid(crs=TARGET, spacing=spacing, row0=0, col0=0, rows=2, cols=2)

    def test_spacing_is_an_integer(self, tg: ModuleType) -> None:
        with pytest.raises(ValidationError):
            tg.TargetGrid(crs=TARGET, spacing=30.5, row0=0, col0=0, rows=2, cols=2)


@pytest.mark.parametrize(
    ("step", "expected"), [(ANADEM_STEP, 30), (GLO30_STEP, 31)], ids=["anadem", "glo30"]
)
def test_the_default_spacing_at_19_south(tg: ModuleType, step: float, expected: int) -> None:
    """Q16 (a): the source's north-south node spacing in metres at the centroid,
    rounded: ANADEM's 29.8 m gives 30, GLO-30's 1 arc-second 31 (D2)."""
    meta = geographic_tile(rough(4, 4), step=step).meta
    assert tg.default_spacing(meta, (-44.0, -19.0)) == expected


def test_the_default_spacing_of_a_projected_source_is_its_own(tg: ModuleType) -> None:
    """Q12 (a): DTM10 gives 10."""
    meta = RasterMeta(
        x_min=500_000.0,
        y_max=6_600_000.0,
        delta_x=10.0,
        delta_y=10.0,
        cols=4,
        rows=4,
        epsg=25833,
        nodata=None,
        nodata_source="absent",
        pixel_is_area=False,
        vertical_unit_assumed=True,
    )
    assert tg.default_spacing(meta, (500_015.0, 6_599_985.0)) == 10


# ------------------------------------------------------------------ G4, resampling


class TestResampling:
    def test_an_affine_source_is_reproduced_at_every_target_node(self, tg: ModuleType) -> None:
        """Bilinear is exact on a function affine in the source's index space,
        so every target node inside the source reads it to 1e-9; every node
        whose stencil leaves the source is NaN (no sentinel in the source)."""
        r, c = np.indices((80, 80), dtype=np.float64)
        tile = geographic_tile(100.0 + 0.5 * c + 0.25 * r)
        grid = covering(tg, tile.meta, margin=3)
        out = tg.resample(grid, tg.TileWindows(tile), threads=4)
        assert out.meta.nodata is None
        row, col = fractional_index(grid, tile.meta)
        inside = (row > 0) & (row < 79) & (col > 0) & (col < 79)
        outside = (row < 0) | (row > 79) | (col < 0) | (col > 79)
        assert inside.sum() > 1000 and outside.sum() > 100, "the grid must straddle the source"
        values = np.asarray(out.array, dtype=np.float64)
        np.testing.assert_allclose(
            values[inside], (100.0 + 0.5 * col + 0.25 * row)[inside], atol=1e-9
        )
        assert np.isnan(values[outside]).all()

    def test_a_nodata_corner_makes_the_node_nodata_whatever_its_weight(
        self, tg: ModuleType
    ) -> None:
        array = rough(80, 80)
        array[40, 40] = NODATA
        tile = geographic_tile(array, nodata=NODATA)
        grid = covering(tg, tile.meta, margin=-5)
        out = tg.resample(grid, tg.TileWindows(tile), threads=4)
        assert out.meta.nodata == NODATA
        row, col = fractional_index(grid, tile.meta)
        integral = (row == np.floor(row)) | (col == np.floor(col))
        touches = (np.floor(row) <= 40) & (np.floor(row) + 1 >= 40)
        touches &= (np.floor(col) <= 40) & (np.floor(col) + 1 >= 40)
        inside = (row > 0) & (row < 79) & (col > 0) & (col < 79) & ~integral
        assert (touches & inside).sum() >= 1
        assert (out.array[touches & inside] == NODATA).all()
        clean = inside & ~touches
        assert np.isfinite(out.array[clean]).all() and (out.array[clean] != NODATA).all()

    def test_identical_for_any_thread_count_and_block_size(self, tg: ModuleType) -> None:
        """J3: each node from its own coordinates alone."""
        tile = geographic_tile(rough(80, 80))
        grid = covering(tg, tile.meta, margin=2)
        source = tg.TileWindows(tile)
        runs = [
            tg.resample(grid, source, threads=threads, block_rows=block)
            for threads, block in ((1, 256), (8, 256), (8, 7), (3, 1))
        ]
        first = runs[0].array.tobytes()  # NaN compares as bytes
        assert all(run.array.tobytes() == first for run in runs[1:])
        assert all(run.meta == runs[0].meta for run in runs[1:])


# ------------------------------------------------------------------ +-inf is data (audit PR A)

# Ola's ruling (`docs/increments/python-audit.md`, section 7, ruling 1;
# `python-audit-pr-a.md`, "The +-infinity ruling" and red test 2): NoData is
# NaN or the sentinel, as in the core's `Raster::is_nodata`; +-inf is data.
# A 6 x 5 source of 10 m cells, one infinity at node (2, 2), resampled in its
# own CRS onto the lattice at 5 m: target node (R, K) reads source (R/2, K/2).
UTM = "EPSG:25833"
INF_AT = (2, 2)
#: Target nodes whose stencil holds the infinity at positive weight (source
#: rows and cols 1.5 to 2.5): +-inf.
POSITIVE = [(r, k) for r in (3, 4, 5) for k in (3, 4, 5)]
#: Target nodes whose stencil holds it at weight 0 (0 x inf): NaN.
ZERO = [(2, 2), (2, 3), (2, 4), (2, 5), (3, 2), (4, 2), (5, 2)]


def plane_tile(dtype: Any, nodata: float | None, at: float | None = None) -> DemTile:
    r, c = np.indices((6, 5))
    z = (100.0 + 3.0 * r + 2.0 * c).astype(dtype)
    if at is not None:
        z[INF_AT] = at
    meta = RasterMeta(
        x_min=500_000.0,
        y_max=6_600_000.0,
        delta_x=10.0,
        delta_y=10.0,
        cols=5,
        rows=6,
        epsg=25833,
        nodata=nodata,
        nodata_source="absent" if nodata is None else "tag",
        pixel_is_area=False,
        vertical_unit_assumed=True,
    )
    return DemTile(meta=meta, array=z)


def half_cells(tg: ModuleType) -> Any:
    return tg.TargetGrid(
        crs=UTM, spacing=5, row0=-6_600_000 // 5, col0=500_000 // 5, rows=11, cols=9
    )


@pytest.mark.parametrize("nodata", [None, -32767.0], ids=["no_sentinel", "sentinel"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("sign", [1.0, -1.0], ids=["+inf", "-inf"])
class TestInfinityIsData:
    def test_resample_carries_it_at_positive_weight_and_is_nan_at_zero_weight(
        self, tg: ModuleType, sign: float, dtype: Any, nodata: float | None
    ) -> None:
        grid = half_cells(tg)
        out = tg.resample(grid, tg.TileWindows(plane_tile(dtype, nodata, sign * math.inf)), 1)
        clean = tg.resample(grid, tg.TileWindows(plane_tile(dtype, nodata)), 1)
        values, reference = np.asarray(out.array), np.asarray(clean.array)
        assert [values[n] for n in POSITIVE] == [sign * math.inf] * len(POSITIVE)
        assert np.isnan([values[n] for n in ZERO]).all(), "0 x inf is NaN, NoData to the core"
        away = np.ones(values.shape, dtype=bool)
        away[tuple(np.transpose(POSITIVE + ZERO))] = False
        assert np.isfinite(reference).all(), "the clean source must give a value everywhere"
        assert values[away].tobytes() == reference[away].tobytes(), "away from it, unchanged"

    def test_check_point_blocks_yields_the_node_with_its_infinity(
        self, tg: ModuleType, sign: float, dtype: Any, nodata: float | None
    ) -> None:
        tile = plane_tile(dtype, nodata, sign * math.inf)
        domain = DomainPolygon(
            polygon=shapely.box(499_990.0, 6_599_940.0, 500_050.0, 6_600_010.0), crs=UTM
        )
        blocks = list(tg.check_point_blocks(half_cells(tg), tg.TileWindows(tile), domain, 1))
        xy = np.concatenate([b[0] for b in blocks])
        z = np.concatenate([b[1] for b in blocks])
        assert xy.shape == (30, 2), "every source node, the infinite one too"
        at = (xy[:, 0] == 500_020.0) & (xy[:, 1] == 6_599_980.0)
        assert at.sum() == 1
        assert z[at][0] == np.float32(sign * math.inf)


# ------------------------------------------------------------------ G5, check points


class TestCheckPoints:
    ROWS = COLS = 600  # four 512 x 512 source blocks, cut at (512, 512)

    @pytest.fixture(scope="class")
    def source(self) -> DemTile:
        array = rough(self.ROWS, self.COLS, seed=5)
        holes = np.random.default_rng(7).integers(0, self.ROWS, (400, 2))
        array[holes[:, 0], holes[:, 1]] = NODATA
        array[[511, 512], [512, 511]] = 650.0  # the corner nodes the arm holds stay valid
        return geographic_tile(array, nodata=NODATA)

    @pytest.fixture(scope="class")
    def domain(self, source: DemTile) -> DomainPolygon:
        """A body in block (0, 0) and a thin arm along the diagonal through the
        blocks' common corner, half a cell and a bit wide: nodes (511, 512)
        and (512, 511) lie inside it, so all four blocks hold inside nodes."""
        meta = source.meta
        body = a_domain(meta, [(100.2, 100.4), (100.1, 400.3), (400.3, 400.2), (400.4, 100.1)])
        w = 0.8 / math.sqrt(2)
        arm = a_domain(
            meta,
            [(390 + w, 390 - w), (560 + w, 560 - w), (560 - w, 560 + w), (390 - w, 390 + w)],
        )
        union = body.polygon.union(arm.polygon)
        assert union.geom_type == "Polygon"
        return DomainPolygon(polygon=shapely.geometry.polygon.orient(union), crs=TARGET)

    def test_every_valid_inside_node_appears_exactly_once(
        self, tg: ModuleType, source: DemTile, domain: DomainPolygon
    ) -> None:
        meta = source.meta
        lon, lat = source_nodes(meta)
        z = source.array.ravel()
        valid = z != NODATA
        xy = project("EPSG:4326", TARGET, lon, lat)
        inside = valid & shapely.contains_xy(domain.polygon, xy[:, 0], xy[:, 1])
        corner = np.zeros((self.ROWS, self.COLS), bool)
        corner[[511, 512], [512, 511]] = True
        assert (inside & corner.ravel()).sum() == 2, "the arm must hold both corner nodes"

        x0, y1 = domain.polygon.bounds[0], domain.polygon.bounds[3]
        x1, y0 = domain.polygon.bounds[2], domain.polygon.bounds[1]
        col0, row0 = math.floor(x0 / H) - 2, math.floor(-y1 / H) - 2
        grid = tg.TargetGrid(
            crs=TARGET,
            spacing=H,
            row0=row0,
            col0=col0,
            rows=math.ceil(-y0 / H) + 2 - row0 + 1,
            cols=math.ceil(x1 / H) + 2 - col0 + 1,
        )
        blocks = list(tg.check_point_blocks(grid, tg.TileWindows(source), domain, threads=4))
        assert all(
            b[0].dtype == np.float64 and b[0].ndim == 2 and b[0].shape[1] == 2 for b in blocks
        )
        assert all(b[1].dtype == np.float32 and b[1].shape == (b[0].shape[0],) for b in blocks)
        got_xy = np.concatenate([b[0] for b in blocks])
        got_z = np.concatenate([b[1] for b in blocks])

        key = {(round(a, 6), round(b, 6)): i for i, (a, b) in enumerate(xy) if valid[i]}
        found = [key.get((round(a, 6), round(b, 6))) for a, b in got_xy]
        assert None not in found, "a check point that is no valid source node"
        index = np.array(found)
        assert np.unique(index).size == index.size, "a source node yielded twice"
        assert (got_z == z[index]).all()
        missing = np.setdiff1d(np.flatnonzero(inside), index)
        assert missing.size == 0, f"{missing.size} inside nodes dropped, e.g. {missing[:5]}"
        # Inside the target grid's node rectangle.
        assert (got_xy[:, 0] >= grid.col0 * H).all()
        assert (got_xy[:, 0] <= (grid.col0 + grid.cols - 1) * H).all()
        assert (got_xy[:, 1] <= -grid.row0 * H).all()
        assert (got_xy[:, 1] >= -(grid.row0 + grid.rows - 1) * H).all()


# ------------------------------------------------------------------ @perf's 15c-2 acceptance


class _UsedFrom:
    """A prepared geometry that records every thread calling into it. Each
    call sleeps 20 ms so blocks overlap in time and the pool's four workers
    all take one; without it, one worker could drain a short queue alone."""

    def __init__(self, inner: Any) -> None:
        self._inner = inner
        self.threads: set[int] = set()

    def __getattr__(self, name: str) -> Any:
        attr = getattr(self._inner, name)
        if not callable(attr):
            return attr

        def call(*args: Any, **kwargs: Any) -> Any:
            self.threads.add(threading.get_ident())
            time.sleep(0.02)
            return attr(*args, **kwargs)

        return call


def test_no_prepared_geometry_is_shared_between_threads(
    tg: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`@perf`'s 15c-2 acceptance: about 15 % of geographic runs crashed
    (SIGBUS, `GEOSException: vector`), each in a pool worker inside GEOS's
    prepared-polygon `intersects`. shapely's prepared geometries are not
    thread-safe, and `check_point_blocks` shared one `prep(domain.polygon)`
    across its workers. Deterministic: every `prep` result is instrumented,
    and none may be called from more than one thread. (A fix that avoids
    `prep` passes trivially; one that calls `shapely.prepare` on the shared
    polygon in place is the same race and is not seen here.)"""
    made: list[_UsedFrom] = []
    real = tg.prep

    def recording(geometry: Any) -> _UsedFrom:
        made.append(_UsedFrom(real(geometry)))
        return made[-1]

    monkeypatch.setattr(tg, "prep", recording)
    monkeypatch.setattr(tg, "BLOCK", 16)  # an 80 x 80 source is 25 blocks
    tile = geographic_tile(rough(80, 80))
    domain = a_domain(tile.meta, [(5.3, 5.7), (6.1, 74.2), (74.6, 73.4), (73.2, 6.8)])
    grid = covering(tg, tile.meta, margin=2)
    blocks = list(tg.check_point_blocks(grid, tg.TileWindows(tile), domain, threads=4))
    assert sum(len(z) for _, z in blocks) > 3000
    shared = [sorted(p.threads) for p in made if len(p.threads) > 1]
    assert not shared, f"a prepared geometry used from {len(shared[0])} threads"


# ------------------------------------------------------------------ 15e, fix 1: node-count blocks
#
# `docs/increments/15e-memory-fixes.md`, fix 1: `BLOCK_NODES = 1 << 20` and
# `rows_per_block(cols) = max(1, BLOCK_NODES // cols)`, read at call time so a
# test can lower it; `resample(..., block_rows=None)` uses it, and an explicit
# `block_rows` still overrides. Went red at 9879805 because neither name existed
# and every block was 256 rows, so one block of an 80 x 80 grid held every node.


class _SpiedReprojector:
    """`crs.reprojector` with every transformer it returns recording the
    number of points it is given. One transformer per block (J5), so the
    record is one entry per block."""

    def __init__(self, real: Callable[[str, str], Callable[[Any], Any]]) -> None:
        self._real = real
        self._lock = threading.Lock()
        self.lengths: list[int] = []

    def __call__(self, source: str, target: str) -> Callable[[Any], Any]:
        transform = self._real(source, target)

        def recording(points: Any) -> Any:
            with self._lock:
                self.lengths.append(len(points))
            return transform(points)

        return recording


@pytest.fixture
def spied(tg: ModuleType, monkeypatch: pytest.MonkeyPatch) -> _SpiedReprojector:
    spy = _SpiedReprojector(tg.reprojector)
    monkeypatch.setattr(tg, "reprojector", spy)
    return spy


class TestBlockSize:
    def test_the_default_is_two_to_the_twenty_nodes(self, tg: ModuleType) -> None:
        assert tg.BLOCK_NODES == 1 << 20

    @pytest.mark.parametrize(
        "cols", [1, 2, 7, 999, 1000, 1001, 1999, 2000, 2001, 10_000], ids=lambda c: f"cols={c}"
    )
    def test_rows_per_block_is_the_most_rows_within_the_budget(
        self, tg: ModuleType, monkeypatch: pytest.MonkeyPatch, cols: int
    ) -> None:
        """`rows x cols <= BLOCK_NODES`, and one more row would not fit; one row
        once a row alone is over the budget (`cols > BLOCK_NODES`)."""
        monkeypatch.setattr(tg, "BLOCK_NODES", 1000, raising=False)
        rows = tg.rows_per_block(cols)
        assert isinstance(rows, int)
        if cols > 1000:
            assert rows == 1
        else:
            assert rows * cols <= 1000 < (rows + 1) * cols

    def test_rows_per_block_reads_the_constant_at_call_time(
        self, tg: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        assert tg.rows_per_block(1024) == 1024
        monkeypatch.setattr(tg, "BLOCK_NODES", 4096)
        assert tg.rows_per_block(1024) == 4

    @pytest.mark.parametrize("budget", ["1000", "5 rows exactly", "under one row"])
    def test_resample_uses_it_by_default(
        self,
        tg: ModuleType,
        monkeypatch: pytest.MonkeyPatch,
        spied: _SpiedReprojector,
        budget: str,
    ) -> None:
        """Each block is one transformer call of `rows_per_block(cols)` rows
        (the last one shorter), so the number of calls is
        `ceil(rows / rows_per_block(cols))` and no call exceeds the budget,
        unless a single row does (then one row per block)."""
        tile = geographic_tile(rough(80, 80))
        grid = covering(tg, tile.meta, margin=2)
        nodes = {"1000": 1000, "5 rows exactly": 5 * grid.cols, "under one row": grid.cols - 1}
        limit = nodes[budget]
        assert grid.cols * 256 > limit and grid.rows < 256, "today's 256 rows must be over it"
        monkeypatch.setattr(tg, "BLOCK_NODES", limit, raising=False)
        tg.resample(grid, tg.TileWindows(tile), threads=4)
        per_block = max(1, limit // grid.cols)  # written out, not `rows_per_block`
        assert len(spied.lengths) == math.ceil(grid.rows / per_block), spied.lengths
        assert max(spied.lengths) == per_block * grid.cols
        assert sum(spied.lengths) == grid.rows * grid.cols
        if budget == "under one row":
            assert set(spied.lengths) == {grid.cols}
        else:
            assert max(spied.lengths) <= limit

    def test_an_explicit_block_rows_still_overrides_it(
        self, tg: ModuleType, monkeypatch: pytest.MonkeyPatch, spied: _SpiedReprojector
    ) -> None:
        tile = geographic_tile(rough(80, 80))
        grid = covering(tg, tile.meta, margin=2)
        monkeypatch.setattr(tg, "BLOCK_NODES", 1000, raising=False)
        tg.resample(grid, tg.TileWindows(tile), threads=4, block_rows=7)
        assert len(spied.lengths) == math.ceil(grid.rows / 7)
        assert max(spied.lengths) == 7 * grid.cols

    def test_the_default_blocks_change_no_value(
        self, tg: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """J3 with the new default: small blocks give today's 256-row bytes."""
        array = rough(80, 80)
        array[40, 40] = NODATA
        tile = geographic_tile(array, nodata=NODATA)
        grid = covering(tg, tile.meta, margin=2)
        source = tg.TileWindows(tile)
        whole = tg.resample(grid, source, threads=1, block_rows=256)
        monkeypatch.setattr(tg, "BLOCK_NODES", 1000, raising=False)
        small = tg.resample(grid, source, threads=4)
        assert small.array.tobytes() == whole.array.tobytes()
        assert small.meta == whole.meta


# ------------------------------------------------------------------ 15e, fix 2: canvas adopted
#
# Fix 2: `resample` returns `DemTile._adopt(meta, canvas)`. Went red at 9879805
# because it called the public constructor, whose array is a read-only view of a copy.


class TestCanvasAdopted:
    def test_the_tile_owns_the_canvas_itself(self, tg: ModuleType) -> None:
        """`base is None`: the array is the `np.empty` canvas, not a view of a copy."""
        tile = geographic_tile(rough(80, 80))
        grid = covering(tg, tile.meta, margin=2)
        out = tg.resample(grid, tg.TileWindows(tile), threads=4)
        held = out.array.base is None
        assert held, "the tile's array is a view of another array (a copy of the canvas)"
        assert out.array.flags.c_contiguous and not out.array.flags.writeable
        assert out.array.shape == (grid.rows, grid.cols)

    def test_adopt_is_called_once_with_the_tiles_array(
        self, tg: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """M15's spy pattern (`test_mosaic.py`, `TestM15Adopt`)."""
        adopted: list[npt.NDArray[Any]] = []
        original = DemTile._adopt

        def spy(meta: Any, array: npt.NDArray[Any]) -> DemTile:
            adopted.append(array)
            return original(meta, array)

        monkeypatch.setattr(DemTile, "_adopt", staticmethod(spy))
        tile = geographic_tile(rough(80, 80))
        out = tg.resample(covering(tg, tile.meta, margin=2), tg.TileWindows(tile), threads=4)
        assert len(adopted) == 1
        assert np.shares_memory(adopted[0], out.array)
