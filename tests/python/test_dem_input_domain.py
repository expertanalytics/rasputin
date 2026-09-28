"""`open_dem` with a domain in its own CRS (increment 15b).

`docs/increments/15-dem-mosaic.md` R1 ("15b adds the domain"), R4 point 5
(the needed region: "the domain polygon grown by one cell"), R6 (the domain's
bounds in the DEM's CRS take `--bbox`'s place, and the two exclude each other)
and R9 (the extent check runs after the transform, against the mosaic's
coverage). Pinned by this suite (see "Pinned by the red suite (15b)"):

- `DemRequest(sources=, bounds=, nodata=, domain=)`: `domain` is a
  `DomainPolygon` as read, in its own CRS, default `None`. A request with both
  `bounds` and `domain` is refused with a `ValueError`.
- `open_dem` moves the domain into the DEM's CRS, plans on its bounds with its
  needed region, checks its extent against the plan, and only then assembles.
  `DemInput.domain` is the domain in the DEM's CRS (`None` without one).
- The plan is the plan of `bounds` equal to the moved domain's bounds, except
  where a domain vertex lies within the 1e-6-cell snap band past a node line:
  there the domain's window is one node line wider on that side. A vertex
  exactly on a node line, or 2e-6 cell past one, widens nothing (review S1).
- "Grown by one cell" is pinned only away from its edge: a missing node 0.73
  cell from a domain vertex refuses the request; missing nodes 2.5 cells or
  more from it are NaN filler. Exactly one cell is not ruled. The cell is the
  chosen plan's, so a tile the domain does not select, and the order names
  sort in, change nothing (review B1).
- Tiles in more than one CRS are refused with a domain, naming the codes
  (review S2).
- Every extent refusal fires before any tile is loaded (I6).
- The seam report (Ola's Q1 revised) survives the domain path: `seams` names
  each disagreeing pair with a domain in the DEM's CRS or in another one, and
  is `()` when the overlaps agree (test amendment after the 15b review).

Synthetic tiles are micro-TIFFs written by `test_dem_input.py`'s helpers;
the real case is the committed DTM10 seam extract (`tests/fixtures/dtm10/`),
with a domain in EPSG:25832, the zone the place lies in.

HOW THIS FILE GOES RED: `tin_engine.dem_input` and `tin_engine.domain` are
imported inside fixtures, and `DemRequest` has no `domain` field before 15b,
so each test fails on its own and collection is unaffected.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any, ClassVar

import numpy as np
import pytest
from pyproj import CRS, Transformer
from shapely.geometry import Polygon, box

from geotiff_fixtures import BASE_KEYS, PROJECTED_CS_TYPE
from mosaic_fixtures import X0, Y0, blocks, piece, quadrants, seams, seams_of, whole
from test_dem_input import SEAM, decoded, tiff_of, write_tiles

Ring = list[tuple[float, float]]
UTM33 = "urn:ogc:def:crs:EPSG::25833"


@pytest.fixture(scope="module")
def di() -> ModuleType:
    return importlib.import_module("tin_engine.dem_input")


@pytest.fixture(scope="module")
def dm() -> ModuleType:
    return importlib.import_module("tin_engine.domain")


@pytest.fixture(scope="module")
def mz() -> ModuleType:
    return importlib.import_module("tin_engine.mosaic")


def to_crs(src: str, dst: str, ring: Ring) -> Ring:
    """pyproj's own `always_xy` transform of `ring`: the oracle."""
    t = Transformer.from_crs(src, dst, always_xy=True)
    xs, ys = t.transform(np.array([p[0] for p in ring]), np.array([p[1] for p in ring]))
    return [(float(x), float(y)) for x, y in zip(xs, ys, strict=True)]


def domain_file(
    path: Path, outer: Ring, holes: tuple[Ring, ...] = (), crs: str | None = UTM33
) -> Path:
    doc: dict[str, Any] = {
        "type": "Polygon",
        "coordinates": [[*r, r[0]] for r in (outer, *holes)],
    }
    if crs is not None:
        doc["crs"] = {"type": "name", "properties": {"name": crs}}
    path.write_text(json.dumps(doc))
    return path


def read(
    dm: ModuleType, tmp_path: Path, utm33: Ring, crs: str, holes: tuple[Ring, ...] = ()
) -> Any:
    """A domain file holding the UTM 33 ring `utm33`, written in `crs`.

    `crs` is "EPSG:25833" (written with its `crs` member), "EPSG:4326" (no
    member, as RFC 7946 has it) or another EPSG code (with its member)."""
    if crs == "EPSG:25833":
        return dm.read_domain(domain_file(tmp_path / "d.geojson", utm33, holes))
    outer = to_crs("EPSG:25833", crs, utm33)
    moved = tuple(to_crs("EPSG:25833", crs, h) for h in holes)
    member = None if crs == "EPSG:4326" else crs
    return dm.read_domain(domain_file(tmp_path / "d.geojson", outer, moved, member))


def request(di: ModuleType, *sources: Path, domain: Any = None, bounds: Any = None) -> Any:
    return di.DemRequest(sources=tuple(sources), bounds=bounds, domain=domain)


@pytest.fixture
def no_load(monkeypatch: pytest.MonkeyPatch) -> None:
    """I6: any tile load from here on fails the test."""
    repository = importlib.import_module("tin_engine.io.repository")

    def refuse(self: Any, name: str, *args: Any, **kwargs: Any) -> Any:
        raise AssertionError(f"load({name!r}) was called")

    monkeypatch.setattr(repository.TiffDemRepository, "load", refuse)


# ------------------------------------------------------------------ tiles

# `whole(9, 13)` in quadrants, point-registered, one shared node line:
# nodes x 500 000 .. 500 120 (dx 10), y 6 600 000 .. 6 599 960 (dy 5); the
# north-west tile holds rows 0-4 and columns 0-6.
IN_NW: Ring = [
    (X0 + 5.3, Y0 - 12.9),
    (X0 + 33.7, Y0 - 12.1),
    (X0 + 32.9, Y0 - 3.1),
    (X0 + 6.1, Y0 - 3.7),
]
ACROSS_ALL: Ring = [
    (X0 + 12.3, Y0 - 36.3),
    (X0 + 107.7, Y0 - 35.9),
    (X0 + 106.1, Y0 - 3.7),
    (X0 + 13.9, Y0 - 4.1),
]


@pytest.fixture
def quad_dir(tmp_path: Path) -> Path:
    write_tiles(tmp_path / "quad", quadrants(whole(9, 13), row_cut=4, col_cut=6, overlap=1))
    return tmp_path / "quad"


# 12 x 12 nodes, dx = dy = 10, cut into four abutting 6 x 6 blocks with the
# south-east one missing: nodes x 500 060 .. 500 110, y 6 599 940 .. 6 599 890
# are in no tile.
@pytest.fixture
def l_dir(tmp_path: Path) -> Path:
    write_tiles(tmp_path / "ell", blocks(whole(12, 12, dy=10.0), 6, 6, skip=[(1, 1)]))
    return tmp_path / "ell"


# A triangle whose corner (53.3, -57.1) has the missing node (60, -60) as a
# bilinear corner: 6.7 m east and 2.9 m south, well within one cell.
NEAR_THE_HOLE: Ring = [(X0 + 5.5, Y0 - 5.5), (X0 + 53.3, Y0 - 57.1), (X0 + 5.5, Y0 - 57.1)]
# An L whose box holds every missing node, and whose nearest point to one,
# (35, -35), is 25 m from each axis of (60, -60): 2.5 cells.
AROUND_THE_HOLE: Ring = [
    (X0 + 5.5, Y0 - 5.5),
    (X0 + 105.5, Y0 - 5.5),
    (X0 + 105.5, Y0 - 35.0),
    (X0 + 35.0, Y0 - 35.0),
    (X0 + 35.0, Y0 - 105.5),
    (X0 + 5.5, Y0 - 105.5),
]


class TestDemRequest:
    def test_domain_defaults_to_none(self, di: ModuleType, quad_dir: Path) -> None:
        assert di.DemRequest(sources=(quad_dir,)).domain is None

    def test_bounds_and_domain_exclude_each_other(
        self, di: ModuleType, dm: ModuleType, mz: ModuleType, quad_dir: Path, tmp_path: Path
    ) -> None:
        """R6 and R11: with a domain, its bounds take `--bbox`'s place."""
        domain = read(dm, tmp_path, IN_NW, "EPSG:25833")
        bounds = mz.Bounds(x_min=X0, y_min=Y0 - 20, x_max=X0 + 50, y_max=Y0)
        with pytest.raises(ValueError):
            request(di, quad_dir, domain=domain, bounds=bounds)

    def test_without_a_domain_the_input_has_none(self, di: ModuleType, quad_dir: Path) -> None:
        assert di.open_dem(request(di, quad_dir)).domain is None


class TestTheDomainChoosesTheTiles:
    @pytest.mark.parametrize("crs", ["EPSG:25833", "EPSG:4326", "EPSG:25832"])
    def test_the_domain_arrives_in_the_dems_crs(
        self, di: ModuleType, dm: ModuleType, quad_dir: Path, tmp_path: Path, crs: str
    ) -> None:
        domain = read(dm, tmp_path, ACROSS_ALL, crs)
        opened = di.open_dem(request(di, quad_dir, domain=domain))
        assert CRS.from_user_input(opened.domain.crs) == CRS.from_epsg(25833)
        given = list(domain.polygon.exterior.coords)
        expected = given if crs == "EPSG:25833" else to_crs(crs, "EPSG:25833", given)
        got = {(float(x), float(y)) for x, y in opened.domain.polygon.exterior.coords}
        assert got == set(expected)

    @pytest.mark.parametrize("crs", ["EPSG:25833", "EPSG:4326"])
    def test_the_plan_is_the_plan_of_the_moved_domains_bounds(
        self,
        di: ModuleType,
        dm: ModuleType,
        mz: ModuleType,
        quad_dir: Path,
        tmp_path: Path,
        crs: str,
    ) -> None:
        opened = di.open_dem(request(di, quad_dir, domain=read(dm, tmp_path, IN_NW, crs)))
        x_min, y_min, x_max, y_max = opened.domain.polygon.bounds
        bounds = mz.Bounds(x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max)
        by_box = di.open_dem(request(di, quad_dir, bounds=bounds))
        assert opened.plan == by_box.plan
        assert [t.name for t in opened.plan.tiles] == ["nw.tif"]

    def test_a_domain_across_every_tile_selects_every_tile(
        self, di: ModuleType, dm: ModuleType, quad_dir: Path, tmp_path: Path
    ) -> None:
        domain = read(dm, tmp_path, ACROSS_ALL, "EPSG:4326")
        opened = di.open_dem(request(di, quad_dir, domain=domain))
        assert [t.name for t in opened.plan.tiles] == ["ne.tif", "nw.tif", "se.tif", "sw.tif"]

    def test_a_single_file_is_cut_to_the_domains_window(
        self, di: ModuleType, dm: ModuleType, quad_dir: Path, tmp_path: Path
    ) -> None:
        """R6 on one file: a window of it, as `--bbox` on one file gives."""
        opened = di.open_dem(
            request(di, quad_dir / "nw.tif", domain=read(dm, tmp_path, IN_NW, "EPSG:4326"))
        )
        # IN_NW spans columns 0..4 and rows 0..3 once snapped outward.
        assert (opened.tile.meta.rows, opened.tile.meta.cols) == (4, 5)
        assert opened.tile.meta.x_min == X0
        assert opened.tile.meta.y_max == Y0


class TestNeededRegion:
    """R4 point 5 with a domain: the polygon grown by one cell, in the DEM's CRS."""

    @pytest.mark.parametrize("crs", ["EPSG:25833", "EPSG:4326"])
    def test_a_missing_node_the_boundary_needs_is_refused(
        self,
        di: ModuleType,
        dm: ModuleType,
        mz: ModuleType,
        l_dir: Path,
        tmp_path: Path,
        no_load: None,
        crs: str,
    ) -> None:
        domain = read(dm, tmp_path, NEAR_THE_HOLE, crs)
        with pytest.raises(mz.MosaicError) as info:
            di.open_dem(request(di, l_dir, domain=domain))
        assert "in no tile" in str(info.value)
        assert "500060" in str(info.value)

    @pytest.mark.parametrize("crs", ["EPSG:25833", "EPSG:4326"])
    def test_missing_nodes_away_from_the_domain_are_filler(
        self, di: ModuleType, dm: ModuleType, l_dir: Path, tmp_path: Path, crs: str
    ) -> None:
        domain = read(dm, tmp_path, AROUND_THE_HOLE, crs)
        opened = di.open_dem(request(di, l_dir, domain=domain))
        tile = opened.tile
        assert (tile.meta.rows, tile.meta.cols) == (12, 12)
        assert np.isnan(tile.array[6:, 6:]).all()
        assert np.isfinite(tile.array[:6, :]).all()
        assert np.isfinite(tile.array[:, :6]).all()

    def test_a_domain_enclosing_a_missing_tile_is_refused(
        self, di: ModuleType, dm: ModuleType, mz: ModuleType, tmp_path: Path, no_load: None
    ) -> None:
        """The interior is needed, not only the boundary."""
        write_tiles(tmp_path / "ring", blocks(whole(18, 18, dy=10.0), 6, 6, skip=[(1, 1)]))
        around: Ring = [
            (X0 + 5.5, Y0 - 5.5),
            (X0 + 164.5, Y0 - 5.5),
            (X0 + 164.5, Y0 - 164.5),
            (X0 + 5.5, Y0 - 164.5),
        ]
        domain = read(dm, tmp_path, around, "EPSG:4326")
        with pytest.raises(mz.MosaicError) as info:
            di.open_dem(request(di, tmp_path / "ring", domain=domain))
        assert "in no tile" in str(info.value)

    def test_a_domain_hole_over_a_missing_tile_is_not_needed(
        self, di: ModuleType, dm: ModuleType, tmp_path: Path
    ) -> None:
        """The domain's own hole, 2.5 cells wider than the missing block on
        every side, is not part of the needed region."""
        write_tiles(tmp_path / "ring", blocks(whole(18, 18, dy=10.0), 6, 6, skip=[(1, 1)]))
        around: Ring = [
            (X0 + 5.5, Y0 - 5.5),
            (X0 + 164.5, Y0 - 5.5),
            (X0 + 164.5, Y0 - 164.5),
            (X0 + 5.5, Y0 - 164.5),
        ]
        # Missing nodes: x 60..110, y -60..-110. The hole: x 35..135, y -35..-135.
        hole: Ring = [
            (X0 + 35.0, Y0 - 35.0),
            (X0 + 35.0, Y0 - 135.0),
            (X0 + 135.0, Y0 - 135.0),
            (X0 + 135.0, Y0 - 35.0),
        ]
        domain = read(dm, tmp_path, around, "EPSG:4326", holes=(hole,))
        opened = di.open_dem(request(di, tmp_path / "ring", domain=domain))
        assert np.isnan(opened.tile.array[6:12, 6:12]).all()


class TestExtentAfterTheTransform:
    """R9: the extent check runs after the transform, against the tiles, and
    fires before any tile is loaded (I6)."""

    def test_a_domain_far_from_every_tile_is_refused(
        self, di: ModuleType, dm: ModuleType, quad_dir: Path, tmp_path: Path, no_load: None
    ) -> None:
        far = [(x - 200_000.0, y + 50_000.0) for x, y in IN_NW]
        with pytest.raises(ValueError):
            di.open_dem(request(di, quad_dir, domain=read(dm, tmp_path, far, "EPSG:4326")))

    @pytest.mark.parametrize("crs", ["EPSG:25833", "EPSG:4326"])
    def test_a_domain_reaching_past_the_tiles_is_refused_as_outside(
        self,
        di: ModuleType,
        dm: ModuleType,
        quad_dir: Path,
        tmp_path: Path,
        no_load: None,
        crs: str,
    ) -> None:
        """The last node column is x = 500 120; one vertex is 3 m past it. The
        window is clamped to the tiles, so the per-node coverage check alone
        would not see it."""
        past = [*ACROSS_ALL[:2], (X0 + 123.0, Y0 - 3.7), ACROSS_ALL[3]]
        with pytest.raises(ValueError) as info:
            di.open_dem(request(di, quad_dir, domain=read(dm, tmp_path, past, crs)))
        assert "outside" in str(info.value)

    def test_a_domain_on_the_last_node_line_is_inside(
        self, di: ModuleType, dm: ModuleType, quad_dir: Path, tmp_path: Path
    ) -> None:
        """16's border rule, unchanged: the node rectangle is closed."""
        corners: Ring = [(X0, Y0 - 40.0), (X0 + 120.0, Y0 - 40.0), (X0 + 120.0, Y0), (X0, Y0)]
        opened = di.open_dem(
            request(di, quad_dir, domain=read(dm, tmp_path, corners, "EPSG:25833"))
        )
        assert (opened.tile.meta.rows, opened.tile.meta.cols) == (9, 13)


class TestRealSeam:
    """The committed DTM10 seam (6400_4 | 6400_1, EPSG:25833, near Flekkefjord,
    which lies in UTM zone 32), with a domain drawn in EPSG:25832."""

    # Across the 51-column overlap at x 49 750 .. 50 250; every vertex off-node.
    ACROSS_THE_SEAM: ClassVar[Ring] = [
        (47_503.3, 6_467_207.7),
        (52_496.1, 6_467_301.9),
        (52_402.7, 6_469_298.3),
        (47_601.9, 6_469_203.1),
    ]

    def test_the_extract_is_where_this_test_thinks(self) -> None:
        node_rect = box(47_190.0, 6_466_980.0, 52_810.0, 6_469_530.0)
        assert node_rect.contains(Polygon(self.ACROSS_THE_SEAM))

    def test_a_utm32_domain_opens_both_tiles(
        self, di: ModuleType, dm: ModuleType, tmp_path: Path
    ) -> None:
        domain = read(dm, tmp_path, self.ACROSS_THE_SEAM, "EPSG:25832")
        opened = di.open_dem(request(di, SEAM, domain=domain))
        assert [t.name for t in opened.plan.tiles] == ["6400_1_10m_z33.tif", "6400_4_10m_z33.tif"]
        given = list(domain.polygon.exterior.coords)
        got = {(float(x), float(y)) for x, y in opened.domain.polygon.exterior.coords}
        assert got == set(to_crs("EPSG:25832", "EPSG:25833", given))


class TestNeededRegionIsGrownByThePlansSpacing:
    """R4 point 5: "grown by one cell" is one cell of the lattice the plan is
    on, so the outcome does not depend on a tile the domain does not select,
    nor on where that tile's name sorts (15b review, B1). Each repository holds
    two spacings in one EPSG: an L of blocks with its south-east block missing,
    and one far-away tile on the other spacing, named to sort first or last."""

    FAR = 5_000.0  # the far tile's west edge, metres east of X0

    @staticmethod
    def far_fine() -> Any:
        """A 1 m tile 5 km east of every coarse node: no domain here meets it."""
        return whole(4, 4, dx=1.0, dy=1.0, x_min=X0 + TestNeededRegionIsGrownByThePlansSpacing.FAR)

    @staticmethod
    def far_coarse() -> Any:
        return whole(
            4, 4, dx=10.0, dy=10.0, x_min=X0 + TestNeededRegionIsGrownByThePlansSpacing.FAR
        )

    @pytest.mark.parametrize("fine", [None, "a_fine.tif", "z_fine.tif"])
    def test_a_node_one_coarse_cell_away_is_needed_whatever_else_is_listed(
        self,
        di: ModuleType,
        dm: ModuleType,
        mz: ModuleType,
        tmp_path: Path,
        no_load: None,
        fine: str | None,
    ) -> None:
        """`NEAR_THE_HOLE` on the 10 m L: the missing node (60, -60) is 7.3 m
        from the domain, inside one 10 m cell and outside one 1 m cell. The
        plan is on the 10 m lattice, so the request is refused, with or
        without a 1 m tile listed first or last."""
        tiles = blocks(whole(12, 12, dy=10.0), 6, 6, skip=[(1, 1)])
        if fine is not None:
            tiles[fine] = self.far_fine()
        write_tiles(tmp_path / "dem", tiles)
        domain = read(dm, tmp_path, NEAR_THE_HOLE, "EPSG:25833")
        with pytest.raises(mz.MosaicError) as info:
            di.open_dem(request(di, tmp_path / "dem", domain=domain))
        assert "1 nodes the request needs are in no tile" in str(info.value)
        assert "500060" in str(info.value)

    @pytest.mark.parametrize("coarse", [None, "a_coarse.tif", "z_coarse.tif"])
    def test_a_node_past_one_fine_cell_is_not_needed_whatever_else_is_listed(
        self, di: ModuleType, dm: ModuleType, tmp_path: Path, coarse: str | None
    ) -> None:
        """The mirror, on a 1 m L (missing nodes x 6..11, y -6..-11) with a far
        10 m tile: `AROUND_THE_HOLE` scaled by a tenth is 2.5 m per axis
        (3.5 m) from the nearest missing node, past one 1 m cell and inside one
        10 m cell. The plan is on the 1 m lattice, so the missing block is NaN
        filler, with or without the 10 m tile listed first or last."""
        tiles = blocks(whole(12, 12, dx=1.0, dy=1.0), 6, 6, skip=[(1, 1)])
        if coarse is not None:
            tiles[coarse] = self.far_coarse()
        write_tiles(tmp_path / "dem", tiles)
        tenth = [(X0 + (x - X0) / 10, Y0 + (y - Y0) / 10) for x, y in AROUND_THE_HOLE]
        opened = di.open_dem(
            request(di, tmp_path / "dem", domain=read(dm, tmp_path, tenth, "EPSG:25833"))
        )
        assert opened.plan.meta.delta_x == 1.0
        assert (opened.tile.meta.rows, opened.tile.meta.cols) == (12, 12)
        assert np.isnan(opened.tile.array[6:, 6:]).all()
        assert np.isfinite(opened.tile.array[:6, :]).all()
        assert np.isfinite(opened.tile.array[:, :6]).all()

    @pytest.mark.parametrize("fine", ["a_fine.tif", "z_fine.tif"])
    def test_an_unselected_tile_does_not_change_the_plan(
        self, di: ModuleType, dm: ModuleType, tmp_path: Path, fine: str
    ) -> None:
        write_tiles(tmp_path / "alone", blocks(whole(12, 12, dy=10.0), 6, 6, skip=[(1, 1)]))
        tiles = blocks(whole(12, 12, dy=10.0), 6, 6, skip=[(1, 1)])
        tiles[fine] = self.far_fine()
        write_tiles(tmp_path / "with", tiles)
        domain = read(dm, tmp_path, AROUND_THE_HOLE, "EPSG:25833")
        alone = di.open_dem(request(di, tmp_path / "alone", domain=domain))
        with_fine = di.open_dem(request(di, tmp_path / "with", domain=domain))
        assert with_fine.plan == alone.plan

    # A rectangle x 5.5..53.3, y -5.5..-57.1 with a notch x 20..30 down to
    # y -20 cut from its north edge. Its south-east corner is `NEAR_THE_HOLE`'s,
    # 7.3 m from the missing 10 m node (60, -60).
    NOTCHED: ClassVar[Ring] = [
        (X0 + 5.5, Y0 - 5.5),
        (X0 + 20.0, Y0 - 5.5),
        (X0 + 20.0, Y0 - 20.0),
        (X0 + 30.0, Y0 - 20.0),
        (X0 + 30.0, Y0 - 5.5),
        (X0 + 53.3, Y0 - 5.5),
        (X0 + 53.3, Y0 - 57.1),
        (X0 + 5.5, Y0 - 57.1),
    ]

    @staticmethod
    def fine(rows: int, cols: int, west: float, north: float) -> Any:
        return whole(rows, cols, dx=1.0, dy=1.0, x_min=X0 + west, y_max=Y0 + north)

    @pytest.mark.parametrize("prefix", ["a_", "z_"])
    def test_the_needed_region_grows_again_when_the_lattice_changes(
        self, di: ModuleType, dm: ModuleType, mz: ModuleType, tmp_path: Path, prefix: str
    ) -> None:
        """The re-plan loop runs until the reach is the chosen lattice's cell.

        The 10 m L plus four 1 m tiles covering `NOTCHED` except a 9 x 15-node
        gap at x 21..29, y -5..-19, which lies 1 m outside the notch. On the
        polygon itself the 1 m lattice is chosen (four tiles against three);
        grown by its 1.41 m diagonal the gap is needed, so the 10 m lattice is
        chosen; grown by that one's 14.1 m diagonal, (60, -60) is needed, so
        the 10 m lattice no longer covers either, and the request is refused
        as mixed-lattice (Q5). The 1 m tiles are selected at every stage:
        selection follows the ungrown box. Growing only
        once, by the first plan's cell, accepts it with (60, -60) as NaN
        filler: the B1 bug by another route."""
        tiles = blocks(whole(12, 12, dy=10.0), 6, 6, skip=[(1, 1)])
        tiles[f"{prefix}left.tif"] = self.fine(54, 16, 5.0, -5.0)
        tiles[f"{prefix}mid.tif"] = self.fine(39, 9, 21.0, -20.0)
        tiles[f"{prefix}rtop.tif"] = self.fine(26, 25, 30.0, -5.0)
        tiles[f"{prefix}rbot.tif"] = self.fine(28, 25, 30.0, -31.0)
        write_tiles(tmp_path / "dem", tiles)
        domain = read(dm, tmp_path, self.NOTCHED, "EPSG:25833")
        with pytest.raises(mz.MosaicError) as info:
            di.open_dem(request(di, tmp_path / "dem", domain=domain))
        assert "the request selects tiles on two lattices" in str(info.value)
        assert f"{prefix}left.tif" in str(info.value)


class TestTheSnapBand:
    """The domain's window against `--bbox`'s, on `quad_dir` (node lines every
    10 m in x, 5 m in y). 15a's window snaps a box edge within 1e-6 cell past
    a node line onto it; a domain vertex there still needs the next line, so
    the domain's window is one line wider on that side (15b, `_past`). A
    vertex exactly on a node line, or 2e-6 cell past one, adds nothing (15b
    review, S1: `<` made `<=` would add a line for every vertex on a node
    line, which a gridded catchment has on every side)."""

    # Every edge on a node line, interior to the 9 x 13 node rectangle:
    # columns 2..6, rows 2..6.
    ON_LINES = (X0 + 20.0, Y0 - 30.0, X0 + 60.0, Y0 - 10.0)
    # Unit vector out of the rectangle for each of x_min, y_min, x_max, y_max,
    # and the (rows, cols) and (x_min, y_max) change one more line there makes.
    OUT = (-1, -1, 1, 1)

    def plans(
        self,
        di: ModuleType,
        dm: ModuleType,
        mz: ModuleType,
        quad_dir: Path,
        tmp_path: Path,
        edges: tuple[float, float, float, float],
    ) -> tuple[Any, Any]:
        """The domain's plan and `--bbox`'s, for the rectangle `edges`."""
        x0, y0, x1, y1 = edges
        ring: Ring = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        domain = read(dm, tmp_path, ring, "EPSG:25833")
        assert domain.polygon.bounds == edges  # bit for bit: same CRS, no transform
        by_domain = di.open_dem(request(di, quad_dir, domain=domain)).plan
        bounds = mz.Bounds(x_min=x0, y_min=y0, x_max=x1, y_max=y1)
        by_box = di.open_dem(request(di, quad_dir, bounds=bounds)).plan
        return by_domain, by_box

    def test_on_node_lines_the_domain_window_is_the_bbox_window(
        self, di: ModuleType, dm: ModuleType, mz: ModuleType, quad_dir: Path, tmp_path: Path
    ) -> None:
        by_domain, by_box = self.plans(di, dm, mz, quad_dir, tmp_path, self.ON_LINES)
        assert (by_box.meta.rows, by_box.meta.cols) == (5, 5)
        assert (by_box.meta.x_min, by_box.meta.y_max) == (X0 + 20.0, Y0 - 10.0)
        assert by_domain == by_box

    @pytest.mark.parametrize("edge", [0, 1, 2, 3], ids=["x_min", "y_min", "x_max", "y_max"])
    @pytest.mark.parametrize("cells", [0.0, 2e-6])
    def test_outside_the_band_the_domain_window_is_the_bbox_window(
        self,
        di: ModuleType,
        dm: ModuleType,
        mz: ModuleType,
        quad_dir: Path,
        tmp_path: Path,
        edge: int,
        cells: float,
    ) -> None:
        edges = list(self.ON_LINES)
        edges[edge] += self.OUT[edge] * cells * (10.0 if edge % 2 == 0 else 5.0)
        by_domain, by_box = self.plans(di, dm, mz, quad_dir, tmp_path, tuple(edges))
        assert by_domain == by_box

    @pytest.mark.parametrize("edge", [0, 1, 2, 3], ids=["x_min", "y_min", "x_max", "y_max"])
    def test_in_the_band_the_domain_window_is_one_line_wider(
        self,
        di: ModuleType,
        dm: ModuleType,
        mz: ModuleType,
        quad_dir: Path,
        tmp_path: Path,
        edge: int,
    ) -> None:
        """A vertex 1e-7 cell past a node line: `--bbox` snaps onto the line,
        the domain's window takes the next one on that side only."""
        edges = list(self.ON_LINES)
        step = 10.0 if edge % 2 == 0 else 5.0
        edges[edge] += self.OUT[edge] * 1e-7 * step
        by_domain, by_box = self.plans(di, dm, mz, quad_dir, tmp_path, tuple(edges))
        exact, _ = self.plans(di, dm, mz, quad_dir, tmp_path, self.ON_LINES)
        assert by_box == exact
        m, e = by_domain.meta, exact.meta
        wider = {
            0: (e.rows, e.cols + 1, e.x_min - 10.0, e.y_max),
            1: (e.rows + 1, e.cols, e.x_min, e.y_max),
            2: (e.rows, e.cols + 1, e.x_min, e.y_max),
            3: (e.rows + 1, e.cols, e.x_min, e.y_max + 5.0),
        }[edge]
        assert (m.rows, m.cols, m.x_min, m.y_max) == wider


class TestOneCrs:
    def test_tiles_in_two_crss_with_a_domain_are_refused_before_any_load(
        self, di: ModuleType, dm: ModuleType, mz: ModuleType, tmp_path: Path, no_load: None
    ) -> None:
        """R6: the domain is moved into the DEM's CRS, so the DEM must have one
        (15b review, S2). Both tiles are in the domain's box."""
        tiles = quadrants(whole(9, 13), row_cut=4, col_cut=6, overlap=1)
        write_tiles(tmp_path / "two", {"nw.tif": tiles["nw.tif"]})
        utm32 = {**BASE_KEYS, PROJECTED_CS_TYPE: 25832}
        (tmp_path / "two" / "ne.tif").write_bytes(
            tiff_of(tiles["ne.tif"], geokeys=utm32).getvalue()
        )
        domain = read(dm, tmp_path, IN_NW, "EPSG:25833")
        with pytest.raises(mz.MosaicError) as info:
            di.open_dem(request(di, tmp_path / "two", domain=domain))
        message = str(info.value)
        assert "the tiles are in 2 CRSs" in message
        assert "EPSG:[25832, 25833]" in message
        assert "a domain needs one" in message


class TestSeamsOnTheDomainPath:
    """Ola's Q1 revised, with a domain: `DemInput.seams` is the mosaic's report,
    as without one. `ne.tif` is planted off `nw.tif` on their shared column
    (global column 6, in both tiles only): +0.5 at row 1, +2.0 at row 2, both
    at least 1 mm, and +0.0005 at row 3, below it. So `ne.tif | nw.tif`,
    nodes 2, max 2.0, median 1.25. Every planted node, (500 060, y 6 599 995
    .. 6 599 985), lies inside `ACROSS_ALL`. A resolution that drops the seams
    whenever a domain is given fails the disagreeing cases (15b review)."""

    @staticmethod
    def disagreeing(tmp_path: Path) -> Path:
        source = whole(9, 13)
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        changed = np.array(tiles["ne.tif"].array)
        for row, by in ((1, 0.5), (2, 2.0), (3, 0.0005)):
            changed[row, 0] += np.float32(by)
        tiles["ne.tif"] = piece(source, 0, 5, 6, 13, array=changed)
        write_tiles(tmp_path / "disagree", tiles)
        return tmp_path / "disagree"

    @pytest.mark.parametrize("crs", ["EPSG:25833", "EPSG:4326", "EPSG:25832"])
    def test_a_disagreeing_pair_is_reported_with_a_domain(
        self, di: ModuleType, dm: ModuleType, tmp_path: Path, crs: str
    ) -> None:
        dem = self.disagreeing(tmp_path)
        domain = read(dm, tmp_path, ACROSS_ALL, crs)
        opened = di.open_dem(request(di, dem, domain=domain))
        assert opened.domain is not None
        assert seams(opened) == [("ne.tif", "nw.tif", 2, 2.0, 1.25)]

    def test_the_report_is_the_one_without_a_domain(
        self, di: ModuleType, dm: ModuleType, mz: ModuleType, tmp_path: Path
    ) -> None:
        """The same report as the node-by-node oracle over the domain's mosaic,
        and as `--bbox` at the moved domain's bounds (the same plan, R6)."""
        dem = self.disagreeing(tmp_path)
        domain = read(dm, tmp_path, ACROSS_ALL, "EPSG:4326")
        opened = di.open_dem(request(di, dem, domain=domain))
        tiles = {p.name: decoded(p) for p in sorted(dem.iterdir())}
        assert seams(opened) == seams_of(tiles, opened.tile.meta)
        x_min, y_min, x_max, y_max = opened.domain.polygon.bounds
        bounds = mz.Bounds(x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max)
        assert seams(opened) == seams(di.open_dem(request(di, dem, bounds=bounds)))

    @pytest.mark.parametrize("crs", ["EPSG:25833", "EPSG:4326"])
    def test_agreeing_overlaps_report_nothing_with_a_domain(
        self, di: ModuleType, dm: ModuleType, quad_dir: Path, tmp_path: Path, crs: str
    ) -> None:
        opened = di.open_dem(request(di, quad_dir, domain=read(dm, tmp_path, ACROSS_ALL, crs)))
        assert [t.name for t in opened.plan.tiles] == ["ne.tif", "nw.tif", "se.tif", "sw.tif"]
        assert opened.seams == ()
