"""23c-3's pieces of the DEM read, one function at a time.

`docs/increments/23-basin-scale.md`, "23c-3: each piece's memory its own",
D1-D5. The end-to-end tests are in `test_cli_mesh_pieces_memory.py` and
`test_cli_mesh_pieces_bytes.py`; these pin the parts they are built from:

- D1: `open_dem` is `assemble_dem` after `plan_dem`, bit for bit.
- D2: `mosaic.restrict` cuts a plan to a window of its canvas on the same
  lattice, and assembling it gives that window of the whole, bit for bit;
  `check_canvas` is the cap's refusal with its wording, and
  `plan_mosaic(..., cap=False)` does not refuse.
- D3: `HeldWindows` reads the run's indices from the part it holds and
  refuses a read outside it; `SourceWindows.dtype` gives `resample` its
  canvas type with no read at the origin, so a sub-grid resampled from a
  held part gets the whole grid's bits there; `grid_meta` is the meta
  `resample` gives; `check_point_blocks(..., within=)` yields exactly the
  whole run's points that lie in `within`, each block judged on its whole
  outline (the clipped-outline case, `BLOCK` patched to 8).
- D5: `basin_run.MARGIN` is 1 and `piece_windows` gives 23c-2's windows;
  `cut_partition` takes the cells and returns their partition.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
import shapely
from shapely.geometry import Polygon

import piece_spies as ps
import pieces_fixtures as pf
from geographic_fixtures import ANADEM_STEP, LAT0, LON0, NODATA, project, rough
from mosaic_fixtures import meta as projected_meta
from mosaic_fixtures import quadrants
from test_cli_mesh_mosaic import write_tiles
from tin_engine.domain import DomainPolygon
from tin_engine.io.models import DemTile, IndexWindow, RasterMeta
from tin_engine.resources import Resources

TARGET = "EPSG:31983"  # SIRGAS 2000 / UTM 23S
H = 30


def same_bits(a: npt.NDArray[Any], b: npt.NDArray[Any]) -> bool:
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


# ------------------------------------------------------------------ D1


def test_open_dem_is_assemble_after_plan(tmp_path: Path) -> None:
    from tin_engine import dem_input
    from tin_engine.io.domain_file import read_domain

    ps.four_tile_scene(tmp_path)
    request = dem_input.DemRequest(
        sources=(tmp_path / "tiles",), domain=read_domain(tmp_path / "square.geojson")
    )
    whole = dem_input.open_dem(request)
    res = Resources.fallback()
    plan = dem_input.plan_dem(request, res)
    assembled = dem_input.assemble_dem(plan, res)
    assert plan.meta == whole.tile.meta
    assert assembled.tile.meta == whole.tile.meta
    assert same_bits(assembled.tile.array, whole.tile.array)
    assert assembled.plan == whole.plan and assembled.seams == whole.seams


# ------------------------------------------------------------------ D2


@pytest.fixture(scope="module")
def four_tiles(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from tin_engine.io.repository import TiffDemRepository

    where = tmp_path_factory.mktemp("tiles")
    whole = DemTile(meta=projected_meta(rows=pf.N, cols=pf.N, dx=pf.H, dy=pf.H), array=pf.hills())
    write_tiles(where, quadrants(whole, row_cut=ps.CUT, col_cut=ps.CUT, overlap=1))
    return TiffDemRepository.from_directory(where)


WINDOWS = {
    "inside nw": IndexWindow(row0=3, col0=4, rows=10, cols=12),
    "across both cuts": IndexWindow(row0=15, col0=17, rows=14, cols=12),
    "on the shared column": IndexWindow(row0=0, col0=23, rows=20, cols=1),
    "the whole canvas": IndexWindow(row0=0, col0=0, rows=pf.N, cols=pf.N),
}


@pytest.mark.parametrize("name", list(WINDOWS))
def test_restrict_is_the_plan_cut_to_a_window(four_tiles: Any, name: str) -> None:
    from tin_engine.mosaic import assemble, plan_mosaic, restrict

    plan = plan_mosaic(four_tiles.footprints())
    w = WINDOWS[name]
    cut = restrict(plan, w)
    assert cut.reference == plan.reference
    assert cut.meta == plan.meta.windowed(w)
    assert cut.window == IndexWindow(
        row0=plan.window.row0 + w.row0, col0=plan.window.col0 + w.col0, rows=w.rows, cols=w.cols
    )
    met = {
        t.name
        for t in plan.tiles
        if t.canvas.row0 < w.row0 + w.rows and w.row0 < t.canvas.row0 + t.canvas.rows
        and t.canvas.col0 < w.col0 + w.cols and w.col0 < t.canvas.col0 + t.canvas.cols
    }  # fmt: skip
    assert {t.name for t in cut.tiles} == met
    for t in cut.tiles:  # inside the window, and the same nodes of the same tile
        c, s = t.canvas, t.source
        assert (
            c.row0 >= 0 and c.row0 + c.rows <= w.rows and c.col0 >= 0 and c.col0 + c.cols <= w.cols
        )
        full = next(p for p in plan.tiles if p.name == t.name)
        assert (s.row0 - full.source.row0, s.col0 - full.source.col0) == (
            c.row0 + w.row0 - full.canvas.row0,
            c.col0 + w.col0 - full.canvas.col0,
        )
    load, by_window = four_tiles.load, four_tiles.load_window
    whole = assemble(plan, load, load_window=by_window).tile
    part = assemble(cut, load, load_window=by_window).tile
    assert part.meta == cut.meta
    rows, cols = slice(w.row0, w.row0 + w.rows), slice(w.col0, w.col0 + w.cols)
    assert same_bits(
        np.ascontiguousarray(part.array), np.ascontiguousarray(whole.array[rows, cols])
    )


def test_check_canvas_refuses_over_the_cap_in_the_caps_words(four_tiles: Any) -> None:
    from tin_engine.mosaic import MosaicError, check_canvas, plan_mosaic

    plan = plan_mosaic(four_tiles.footprints(), cap=False, memory_budget=1)
    size = pf.N * pf.N * 4
    check_canvas(plan, size)  # at the cap: passes
    with pytest.raises(MosaicError) as caught:
        check_canvas(plan, size - 1)
    assert str(caught.value).startswith(
        f"the request is {pf.N} x {pf.N} nodes, {size} bytes at float32, over the memory "
        f"budget of {size - 1} bytes (set with --memory-budget"
    ), str(caught.value)


def test_plan_mosaic_still_caps_by_default(four_tiles: Any) -> None:
    from tin_engine.mosaic import MosaicError, plan_mosaic

    with pytest.raises(MosaicError, match="over the memory budget of 1 bytes"):
        plan_mosaic(four_tiles.footprints(), memory_budget=1)


# ------------------------------------------------------------------ D3


def geographic(array: npt.NDArray[Any], nodata: float | None = None) -> DemTile:
    meta = RasterMeta(
        x_min=LON0, y_max=LAT0, delta_x=ANADEM_STEP, delta_y=ANADEM_STEP,
        cols=array.shape[1], rows=array.shape[0], epsg=4326,
        crs="EPSG:4326",  # type: ignore[call-arg]
        geographic=True,  # type: ignore[call-arg]
        nodata=nodata, nodata_source="absent" if nodata is None else "tag",
        pixel_is_area=False, vertical_unit_assumed=True,
    )  # fmt: skip
    return DemTile(meta=meta, array=array)


def held_part(tile: DemTile, r0: int, c0: int) -> DemTile:
    """Rows and columns from `(r0, c0)` on, as a tile of their own."""
    w = IndexWindow(row0=r0, col0=c0, rows=tile.meta.rows - r0, cols=tile.meta.cols - c0)
    return DemTile(meta=tile.meta.windowed(w), array=tile.array[r0:, c0:])


def covering(meta: RasterMeta) -> Any:
    """A target grid over the source's whole image in TARGET."""
    from tin_engine.target_grid import TargetGrid

    r, c = np.indices((meta.rows, meta.cols), dtype=np.float64)
    lon, lat = meta.x_min + c.ravel() * meta.delta_x, meta.y_max - r.ravel() * meta.delta_y
    xy = project("EPSG:4326", TARGET, lon, lat)
    col0, col1 = math.floor(xy[:, 0].min() / H), math.ceil(xy[:, 0].max() / H)
    row0, row1 = math.floor(-xy[:, 1].max() / H), math.ceil(-xy[:, 1].min() / H)
    return TargetGrid(
        crs=TARGET, spacing=H, row0=row0, col0=col0, rows=row1 - row0 + 1, cols=col1 - col0 + 1
    )


def central_third(grid: Any) -> Any:
    r0, c0 = grid.rows // 3, grid.cols // 3
    return grid.model_copy(
        update={"row0": grid.row0 + r0, "col0": grid.col0 + c0, "rows": r0, "cols": c0}
    )


class TestHeldWindows:
    @pytest.fixture(scope="class")
    def tile(self) -> DemTile:
        return geographic(rough(40, 40, seed=3))

    def test_reads_the_run_indices_from_the_part_it_holds(self, tile: DemTile) -> None:
        from tin_engine.target_grid import HeldWindows

        held = HeldWindows(tile.meta, held_part(tile, 5, 7), 5, 7)
        assert held.meta == tile.meta
        assert held.dtype == np.float32
        assert same_bits(np.asarray(held.window(5, 9, 7, 12)), tile.array[5:9, 7:12])
        assert same_bits(np.asarray(held.window(30, 40, 30, 40)), tile.array[30:40, 30:40])

    @pytest.mark.parametrize("box", [(4, 9, 7, 12), (5, 9, 6, 12), (0, 1, 0, 1), (5, 41, 7, 12)])
    def test_refuses_a_read_outside_what_it_holds(
        self, tile: DemTile, box: tuple[int, int, int, int]
    ) -> None:
        from tin_engine.target_grid import HeldWindows

        held = HeldWindows(tile.meta, held_part(tile, 5, 7), 5, 7)
        with pytest.raises(ValueError):
            held.window(*box)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_tile_windows_dtype_is_its_arrays(dtype: Any) -> None:
    from tin_engine.target_grid import TileWindows

    tile = geographic(rough(12, 12).astype(dtype))
    assert TileWindows(tile).dtype == np.dtype(dtype)


class _Recording:
    """A `SourceWindows` over a tile that records every read."""

    def __init__(self, inner: Any) -> None:
        self.inner, self.reads = inner, []

    @property
    def meta(self) -> RasterMeta:
        return self.inner.meta  # type: ignore[no-any-return]

    @property
    def dtype(self) -> np.dtype[Any]:
        return self.inner.dtype  # type: ignore[no-any-return]

    def window(self, r0: int, r1: int, c0: int, c1: int) -> Any:
        self.reads.append((r0, r1, c0, c1))
        return self.inner.window(r0, r1, c0, c1)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_resample_takes_its_type_from_the_source_not_a_read_at_the_origin(dtype: Any) -> None:
    """The `SourceWindows.dtype` case: a sub-grid well inside the source,
    resampled from a part held from row and column 1 on (where a read at
    (0, 0) raises), gets the whole grid's bits and the source's type."""
    from tin_engine.target_grid import HeldWindows, TileWindows, grid_meta, resample

    tile = geographic(rough(40, 40, seed=3).astype(dtype))
    grid = covering(tile.meta)
    sub = central_third(grid)
    whole = resample(grid, TileWindows(tile), 1)
    recording = _Recording(HeldWindows(tile.meta, held_part(tile, 1, 1), 1, 1))
    part = resample(sub, recording, 1)
    assert (0, 1, 0, 1) not in recording.reads
    assert part.array.dtype == np.dtype(dtype)
    r, c = sub.row0 - grid.row0, sub.col0 - grid.col0
    expect = np.ascontiguousarray(whole.array[r : r + sub.rows, c : c + sub.cols])
    assert np.isfinite(expect).all(), "the sub-grid must lie inside the source's image"
    assert same_bits(np.ascontiguousarray(part.array), expect)
    assert part.meta == grid_meta(sub, tile.meta)
    assert whole.meta == grid_meta(grid, tile.meta)


class TestCheckPointsWithin:
    """`within` keeps the blocks that meet it, each judged on its whole
    outline, values read clipped. Source 40 x 40 nodes at about 30 m with
    NoData holes; `BLOCK` 8, so 25 blocks; the domain a quadrilateral over
    source rows and columns 5 to 22; the target grid over the whole source,
    so no point is dropped for leaving it."""

    @pytest.fixture(scope="class")
    def source(self) -> DemTile:
        array = rough(40, 40, seed=11)
        holes = np.random.default_rng(4).integers(0, 40, (60, 2))
        array[holes[:, 0], holes[:, 1]] = NODATA
        array[23, 23] = 500.0  # the witness node below has data
        return geographic(array, nodata=NODATA)

    @pytest.fixture(scope="class")
    def domain(self, source: DemTile) -> DomainPolygon:
        m = source.meta
        rc = [(5.2, 5.4), (5.1, 22.0), (22.0, 22.0), (22.0, 5.1)]
        lon = np.array([m.x_min + c * m.delta_x for _, c in rc])
        lat = np.array([m.y_max - r * m.delta_y for r, _ in rc])
        xy = project("EPSG:4326", TARGET, lon, lat)
        return DomainPolygon(polygon=shapely.geometry.polygon.orient(Polygon(xy)), crs=TARGET)

    #: (row0, col0, rows, cols) in the source's indices. The last holds of
    #: block (16..23, 16..23), which meets the domain, only node (23, 23),
    #: about 42 m past the domain's corner (22, 22): its own outline's image,
    #: grown by h = 30 m, misses the domain.
    WITHIN = ((0, 0, 40, 40), (3, 11, 14, 18), (9, 9, 9, 9), (23, 23, 10, 10))

    def points(self, blocks: list[Any]) -> dict[tuple[float, float], float]:
        out: dict[tuple[float, float], float] = {}
        for xy, z in blocks:
            for (x, y), v in zip(xy.tolist(), z.tolist(), strict=True):
                assert (x, y) not in out, "a point yielded twice"
                out[(x, y)] = v
        return out

    def node_points(self, m: RasterMeta, w: tuple[int, int, int, int]) -> set[tuple[float, float]]:
        """The target position of every node of `w`, by pyproj alone."""
        return set(self.nodes_at(m, w))

    def nodes_at(
        self, m: RasterMeta, w: tuple[int, int, int, int]
    ) -> dict[tuple[float, float], tuple[int, int]]:
        r0, c0, rows, cols = w
        r, c = np.indices((rows, cols))
        r, c = (r0 + r).ravel(), (c0 + c).ravel()
        lon, lat = m.x_min + c * m.delta_x, m.y_max - r * m.delta_y
        xy = project("EPSG:4326", TARGET, lon, lat).tolist()
        return {
            (float(x), float(y)): (int(i), int(j)) for (x, y), i, j in zip(xy, r, c, strict=True)
        }

    @pytest.mark.parametrize("w", WITHIN)
    def test_yields_the_whole_runs_points_inside_within(
        self, monkeypatch: pytest.MonkeyPatch, source: DemTile, domain: DomainPolygon, w: Any
    ) -> None:
        """Equal to the last bit: both calls move the same nodes by the same
        transform. A point is told to its node by its position to 1 mm
        (nodes are about 30 m apart), from pyproj alone."""
        from tin_engine import target_grid as tg

        monkeypatch.setattr(tg, "BLOCK", 8)
        grid = covering(source.meta)
        everything = self.points(
            list(tg.check_point_blocks(grid, tg.TileWindows(source), domain, 1))
        )
        within = IndexWindow(row0=w[0], col0=w[1], rows=w[2], cols=w[3])
        got = self.points(
            list(tg.check_point_blocks(grid, tg.TileWindows(source), domain, 1, within=within))
        )
        m = source.meta
        node = {
            (round(x, 3), round(y, 3)): rc
            for (x, y), rc in self.nodes_at(m, (0, 0, 40, 40)).items()
        }
        located = {k: node.get((round(k[0], 3), round(k[1], 3))) for k in everything}
        assert None not in located.values(), "a check point at no source node"
        r0, c0, rows, cols = w
        keep = {
            k: v for k, v in everything.items()
            if r0 <= located[k][0] < r0 + rows and c0 <= located[k][1] < c0 + cols  # type: ignore[index]
        }  # fmt: skip
        assert got == keep
        if w == (23, 23, 10, 10):  # the witness is yielded, its block judged whole
            witness = self.node_points(m, (23, 23, 1, 1)).pop()
            assert any(abs(x - witness[0]) < 1e-3 and abs(y - witness[1]) < 1e-3 for x, y in got)

    def test_the_witness_is_one_a_clipped_outline_would_skip(
        self, source: DemTile, domain: DomainPolygon
    ) -> None:
        """The clipped-outline mutant's food: node (23, 23) has data, its
        block's whole outline meets the domain and its own image grown by h
        does not."""
        m = source.meta
        (point,) = self.node_points(m, (23, 23, 1, 1))
        assert not domain.polygon.intersects(shapely.Point(point).buffer(H))
        block = sorted(self.node_points(m, (16, 16, 8, 8)))
        assert shapely.MultiPoint(block).convex_hull.buffer(H).intersects(domain.polygon)
        assert source.array[23, 23] != NODATA


# ------------------------------------------------------------------ D5


def run_meta(n: int = 37) -> RasterMeta:
    return projected_meta(rows=n, cols=n, dx=pf.H, dy=pf.H)


@pytest.mark.parametrize("pieces", [4, 16])
def test_piece_windows_are_23c2s_cells_grown_by_one_node(pieces: int) -> None:
    from tin_engine import basin_run
    from tin_engine.decompose import partition

    assert basin_run.MARGIN == 1
    meta = run_meta()
    cells = partition(meta.cols, meta.rows, 1.0, pieces=pieces)
    ask = basin_run.CutAsk(pieces, 2**34, (pf.TIE_X, pf.TIE_Y), 2, 2)
    cut = basin_run.cut_partition(meta, 1.0, cells, ask)
    assert (cut.record.nx, cut.record.ny, cut.record.dx, cut.record.dy) == (
        cells.nx, cells.ny, cells.dx, cells.dy
    )  # fmt: skip
    found = {k: ps.window_bounds(v) for k, v in basin_run.piece_windows(cut).items()}
    dx, dy, last = cells.dx, cells.dy, meta.rows - 1
    expected = {
        (j, i): (max(j * dy - 1, 0), min((j + 1) * dy + 1, last),
                 max(i * dx - 1, 0), min((i + 1) * dx + 1, last))
        for j in range(cells.ny) for i in range(cells.nx)
    }  # fmt: skip
    assert found == expected
