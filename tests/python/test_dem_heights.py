"""`tin_engine.dem_heights` and `io.models.HeightTable`: the classes of a
run, the mesh's corners into the DEM's index space chunk by chunk, then
`_core.dem_height_shares`.

`docs/increments/137-dem-stats.md`, 3.1, 3.2 and 3.4, tests 3 and 22 of
section 5, invariant 4 at this level, and section 6's degeneracies.
Hand-built tiles in EPSG:25833 (10 m nodes, point-registered, node (0, 0) at
`X0`, `Y0`): a vertex at node (r, c) lies at index (c + 1/2, r + 1/2).

PINNED BY THIS SUITE (3.4 gives the names; the spelling of the rest):
`height_shares` returns `(HeightTable, int)`, the table's `edges_m` and
`class_m` those of the run; a vertex past the tile, or with no image, is a
`ValueError` saying `outside` before the binding; `HeightTable` refusals are
pydantic `ValidationError`s (a `ValueError`) naming the field.
"""

from __future__ import annotations

import io
import math
from typing import Any

import numpy as np
import pyproj
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from height_fixtures import VELHAS, edges, oracle_table, table_of, velhas_source
from tin_engine import dem_heights
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import DemTile, HeightTable, RasterMeta
from tin_engine.target_grid import HeldWindows, TileWindows

X0, Y0, H = 500_000.0, 6_600_000.0, 10.0
SHARE = 1e-7


def tile(array: np.ndarray, nodata: float | None = None) -> DemTile:
    data = np.ascontiguousarray(array)
    data.setflags(write=False)
    rows, cols = data.shape
    meta = RasterMeta(
        x_min=X0,
        y_max=Y0,
        delta_x=H,
        delta_y=H,
        cols=cols,
        rows=rows,
        epsg=25833,
        crs="EPSG:25833",
        nodata=nodata,
        nodata_source="absent" if nodata is None else "tag",
        pixel_is_area=False,
        vertical_unit_assumed=True,
    )
    return DemTile(meta=meta, array=data)


def world(index: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Index-space points with their z, as mesh vertices in the tile's CRS."""
    col, row = index[:, 0], index[:, 1]
    return np.column_stack([X0 + (col - 0.5) * H, Y0 - (row - 0.5) * H, z])


def steps() -> np.ndarray:
    """Test 1's DEM: 11 x 11 nodes, 2 m on columns 1 and 2, 0 elsewhere."""
    a = np.zeros((11, 11), np.float32)
    a[:, 1:3] = 2.0
    return a


#: Test 1's mesh: the node rectangle's corners, split from north-west to
#: south-east; triangle 0 holds both west corners.
STEP_VERTICES = world(np.array([[0.5, 0.5], [10.5, 0.5], [0.5, 10.5], [10.5, 10.5]]), np.zeros(4))
STEP_TRIANGLES = np.array([[0, 3, 2], [0, 1, 3]], dtype=np.uint32)


def run(windows: Any, vertices: np.ndarray, triangles: np.ndarray, **kw: Any) -> Any:
    options: dict[str, Any] = {"tolerance_m": 10.0, "class_m": 1.0, "mesh_crs": None, "threads": 1}
    options.update(kw)
    return dem_heights.height_shares(windows, vertices, triangles, **options)


# ------------------------------------------------- test 22: the classes


class TestHalfCount:
    """3.1: n is the least whole number >= 0 with (n + 1/2) w > tolerance."""

    @pytest.mark.parametrize(
        ("tolerance", "width", "n"),
        [
            (1, 1, 1),
            (0.4, 1, 0),
            (0.5, 1, 1),
            (2.5, 1, 3),
            (10, 1, 10),
            (10, 2.5, 4),
            (5, 1.216, 4),
            (127.4, 1, 127),
            (2, 1, 2),
            (20, 1, 20),
        ],
    )
    def test_n(self, tolerance: float, width: float, n: int) -> None:
        assert dem_heights.half_count(tolerance, width) == n

    def test_a_tolerance_on_an_edge_whose_quotient_rounds_low(self) -> None:
        """The margin is strict against the product, not the quotient: at
        w = 1.216 the tolerance (7 + 1/2) w gives tolerance / w + 1/2 just
        under 8, so the floor alone says 7; (7 + 1/2) w is not above the
        tolerance, so n is 8. Kills a raise step taken while < (non-strict)."""
        width = 1.216
        tolerance = (7 + 0.5) * width
        assert math.floor(tolerance / width + 0.5) == 7
        assert dem_heights.half_count(tolerance, width) == 8

    @pytest.mark.parametrize(
        ("tolerance", "width"),
        [
            (127.5, 1),
            (1, 0),
            (0, 1),
            (math.nan, 1),
            (1, math.nan),
            (-1, 1),
            (1, -1),
            (math.inf, 1),
            (1, math.inf),
        ],
    )
    def test_refused(self, tolerance: float, width: float) -> None:
        with pytest.raises(ValueError):
            dem_heights.half_count(tolerance, width)

    def test_the_constants(self) -> None:
        assert dem_heights.DEFAULT_CLASS_M == 1.0
        assert dem_heights.MAX_HALF == 127


class TestEdges:
    def test_edges_of_n_1_at_1_m_are_exact(self) -> None:
        got = dem_heights.edges_m(1, 1.0)
        assert got.dtype == np.float64
        assert got.tolist() == [-1.5, -0.5, 0.5, 1.5]

    @pytest.mark.parametrize(("half", "width"), [(10, 1.0), (4, 2.5), (4, 1.216), (0, 1.0)])
    def test_the_edges_are_one_product_each(self, half: int, width: float) -> None:
        got = dem_heights.edges_m(half, width)
        assert got.shape == (2 * half + 2,)
        assert got.tobytes() == edges(half, width).tobytes()

    def test_2_5_m_classes_reach_11_25(self) -> None:
        got = dem_heights.edges_m(4, 2.5)
        assert (got[0], got[-1]) == (-11.25, 11.25)


# ------------------------------------------------- height_shares


class TestInTheDemsCrs:
    """`mesh_crs=None`: the vertices are in the tile's CRS."""

    def test_the_two_triangles_of_test_1(self) -> None:
        table, without = run(TileWindows(tile(steps())), STEP_VERTICES, STEP_TRIANGLES)
        assert isinstance(table, HeightTable)
        assert without == 0
        assert table.class_m == 1.0
        assert table.edges_m.tobytes() == edges(10, 1.0).tobytes()
        assert_array_equal(table.offsets, [0, 3, 6])
        assert_array_equal(table.first, [10, 10])
        assert_array_equal(table.count, [3, 3])
        assert_allclose(table.shares, [33 / 50, 0, 17 / 50, 47 / 50, 0, 3 / 50], rtol=0, atol=SHARE)
        assert_array_equal(table.shares[[1, 4]], 0.0)
        assert_array_equal(table.nodata, 0.0)

    def test_the_tolerance_sets_n(self) -> None:
        table, _ = run(TileWindows(tile(steps())), STEP_VERTICES, STEP_TRIANGLES, tolerance_m=2.0)
        assert table.edges_m.tobytes() == edges(2, 1.0).tobytes()
        assert_array_equal(table.first, [2, 2])

    def test_the_width_is_the_callers(self) -> None:
        table, _ = run(TileWindows(tile(steps())), STEP_VERTICES, STEP_TRIANGLES, class_m=2.5)
        assert table.class_m == 2.5
        assert table.edges_m.tobytes() == edges(4, 2.5).tobytes()
        assert_array_equal(table.first, [4, 4])
        assert_allclose(table.shares, [33 / 50, 17 / 50, 47 / 50, 3 / 50], rtol=0, atol=SHARE)

    def test_a_tolerance_past_127_classes_is_refused(self) -> None:
        with pytest.raises(ValueError):
            run(TileWindows(tile(steps())), STEP_VERTICES, STEP_TRIANGLES, tolerance_m=200.0)

    def test_chunks_of_one_triangle_give_the_same_bytes(self) -> None:
        rng = np.random.default_rng(137)
        dem = rng.uniform(-6.0, 6.0, (30, 30)).astype(np.float32)
        grid = np.linspace(0.7, 29.3, 7)
        index = np.array([[c, r] for r in grid for c in grid])
        index += rng.uniform(-0.3, 0.3, index.shape)
        vertices = world(index, rng.uniform(-3.0, 3.0, len(index)))
        quads = [7 * r + c for r in range(6) for c in range(6)]
        triangles = np.array(
            [t for q in quads for t in ([q, q + 1, q + 8], [q, q + 8, q + 7])], dtype=np.uint32
        )
        windows = TileWindows(tile(dem))
        whole, _ = run(windows, vertices, triangles, threads=2)
        one, _ = run(windows, vertices, triangles, threads=2, chunk=1)
        seven, _ = run(windows, vertices, triangles, threads=2, chunk=7)
        for name in ("offsets", "first", "shares", "nodata"):
            assert getattr(whole, name).tobytes() == getattr(one, name).tobytes(), name
            assert getattr(whole, name).tobytes() == getattr(seven, name).tobytes(), name
        want = oracle_table(
            index[triangles],
            vertices[triangles][:, :, 2],
            dem,
            None,
            np.full(30, H * H),
            class_m=1.0,
            half=10,
        )
        assert_array_equal(whole.first, want.first)
        assert_array_equal(whole.offsets, want.offsets)
        assert_allclose(whole.shares, want.shares, rtol=0, atol=1e-6)


class TestAllNoData:
    """Test 3: a triangle inside a NoData pixel; section 6's NaN twin."""

    INDEX = np.array([[1.1, 1.1], [1.4, 1.1], [1.1, 1.4]])

    @pytest.mark.parametrize("sentinel", [-9999.0, None])
    def test_count_0_first_0_nodata_1_and_counted(self, sentinel: float | None) -> None:
        dem = np.zeros((3, 3), np.float32)
        dem[1, 1] = np.nan if sentinel is None else sentinel
        table, without = run(
            TileWindows(tile(dem, nodata=sentinel)),
            world(self.INDEX, np.zeros(3)),
            np.array([[0, 1, 2]], dtype=np.uint32),
        )
        assert without == 1
        assert_array_equal(table.offsets, [0, 0])
        assert_array_equal(table.first, [0])
        assert table.nodata.tolist() == [1.0]
        assert len(table.shares) == 0


class TestOutside:
    """Section 6: a mesh past the DEM is a `ValueError` before the binding."""

    def test_a_vertex_past_the_tile(self) -> None:
        vertices = STEP_VERTICES.copy()
        vertices[3, 0] += 30.0
        with pytest.raises(ValueError, match="outside"):
            run(TileWindows(tile(steps())), vertices, STEP_TRIANGLES)

    def test_a_vertex_with_no_image(self) -> None:
        vertices = STEP_VERTICES.copy()
        vertices[3, 1] = np.inf
        with pytest.raises(ValueError):
            run(TileWindows(tile(steps())), vertices, STEP_TRIANGLES)

    def test_a_window_past_a_held_part_is_its_value_error(self) -> None:
        source = tile(steps())
        held = HeldWindows(source.meta, tile(steps()[:6, :6]), 0, 0)
        with pytest.raises(ValueError, match="outside the part held"):
            run(held, STEP_VERTICES, STEP_TRIANGLES)

    def test_a_held_part_holding_the_mesh_reads_the_sources_indices(self) -> None:
        """The mesh's corners lie within rows and columns 2 to 8; the part
        holds rows and columns 1 to 9 only."""
        index = np.array([[2.5, 2.5], [8.5, 2.5], [2.5, 8.5]])
        vertices, triangles = world(index, np.zeros(3)), np.array([[0, 1, 2]], np.uint32)
        whole, _ = run(TileWindows(tile(steps())), vertices, triangles)
        held = HeldWindows(tile(steps()).meta, tile(steps()[1:10, 1:10]), 1, 1)
        part, _ = run(held, vertices, triangles)
        assert whole.shares.tobytes() == part.shares.tobytes()
        assert whole.first.tobytes() == part.first.tobytes()


class TestReprojected:
    """3.5 and invariant 4 at this level: `mesh_crs` given, the vertices are
    moved to the source's CRS and measured on its own pixels with its own
    ellipsoidal row areas. Velhas (ANADEM, EPSG:4674, area-registered) and a
    fan of triangles in EPSG:31983 inside it."""

    def test_against_the_oracle_on_the_sources_pixels(self) -> None:
        source = velhas_source()
        with (VELHAS / "anadem_velhas.tif").open("rb") as f:
            dem = decode_dem(io.BytesIO(f.read()))
        rng = np.random.default_rng(31983)
        centre = np.array([source.x0 + 160 * source.dx, source.y0 - 160 * source.dy])
        to_utm = pyproj.Transformer.from_crs("EPSG:4674", "EPSG:31983", always_xy=True)
        cx, cy = to_utm.transform(*centre)
        angles = np.sort(rng.uniform(0, 2 * np.pi, 40))
        radius = rng.uniform(300.0, 1500.0, 40)
        ring = np.column_stack([cx + radius * np.cos(angles), cy + radius * np.sin(angles)])
        xy = np.vstack([[cx, cy], ring])
        z = rng.uniform(700.0, 1100.0, len(xy))
        vertices = np.column_stack([xy, z])
        triangles = np.array([[0, i, i % 40 + 1] for i in range(1, 41)], dtype=np.uint32)
        table, without = run(
            TileWindows(dem), vertices, triangles, tolerance_m=200.0 / 2.0, mesh_crs="EPSG:31983"
        )
        assert without == 0
        index = source.index(xy, "EPSG:31983")
        want = oracle_table(
            index[triangles],
            z[triangles],
            source.dem,
            source.nodata,
            source.row_areas(),
            class_m=1.0,
            half=100,
        )
        assert_array_equal(table.first, want.first)
        assert_array_equal(table.offsets, want.offsets)
        assert_allclose(table.shares, want.shares, rtol=0, atol=1e-6)


# ------------------------------------------------- HeightTable (3.4)


def good() -> dict[str, Any]:
    return {
        "offsets": np.array([0, 3, 3, 4], np.uint32),
        "first": np.array([1, 0, 2], np.int32),
        "shares": np.array([0.25, 0.25, 0.5, 1.0], np.float32),
        "nodata": np.array([0.0, 1.0, 0.0], np.float32),
        "edges_m": edges(2, 1.0),
        "class_m": 1.0,
    }


class TestHeightTable:
    def test_a_good_table_and_its_count(self) -> None:
        table = HeightTable(**good())
        assert_array_equal(table.count, [3, 0, 1])
        with pytest.raises(ValueError):
            table.class_m = 2.0

    @pytest.mark.parametrize(
        ("field", "value"),
        [
            ("offsets", np.array([0, 3, 3, 4], np.int64)),
            ("first", np.array([1, 0, 2], np.int64)),
            ("shares", np.array([0.25, 0.25, 0.5, 1.0], np.float64)),
            ("nodata", np.array([0.0, 1.0, 0.0], np.float64)),
            ("edges_m", edges(2, 1.0).astype(np.float32)),
            ("shares", np.array([[0.25, 0.25], [0.5, 1.0]], np.float32)),
            ("offsets", np.array([1, 3, 3, 4], np.uint32)),
            ("offsets", np.array([0, 3, 2, 4], np.uint32)),
            ("offsets", np.array([0, 3, 3, 5], np.uint32)),
            ("first", np.array([1, 0], np.int32)),
            ("nodata", np.array([0.0, 1.0], np.float32)),
            ("edges_m", np.array([0.5])),
            ("edges_m", edges(2, 1.0)[:-1]),
            ("edges_m", edges(2, 1.0)[::-1].copy()),
            ("first", np.array([3, 0, 2], np.int32)),
            ("first", np.array([-1, 0, 2], np.int32)),
            ("first", np.array([1, 0, 5], np.int32)),
        ],
        ids=[
            "offsets int64",
            "first int64",
            "shares float64",
            "nodata float64",
            "edges float32",
            "shares 2-D",
            "offsets from 1",
            "offsets decreasing",
            "offsets past the shares",
            "first short",
            "nodata short",
            "one edge",
            "odd edges",
            "edges decreasing",
            "first + count past the classes",
            "first negative",
            "a last row past the classes",
        ],
    )
    def test_refused_naming_the_field(self, field: str, value: np.ndarray) -> None:
        with pytest.raises(ValueError, match=field):
            HeightTable(**{**good(), field: value})

    def test_the_arrays_are_kept_as_given(self) -> None:
        given = good()
        table = HeightTable(**given)
        assert table.shares is given["shares"]

    def test_a_table_of_no_triangles(self) -> None:
        empty = table_of([0], [], [], [])
        assert len(empty.first) == 0
        table = HeightTable(
            offsets=np.zeros(1, np.uint32),
            first=np.zeros(0, np.int32),
            shares=np.zeros(0, np.float32),
            nodata=np.zeros(0, np.float32),
            edges_m=edges(0, 1.0),
            class_m=1.0,
        )
        assert len(table.count) == 0
