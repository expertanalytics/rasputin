"""DC9, the partition's geometry, end to end; and SP3's 23c case on a mosaic.

`docs/increments/23-basin-scale.md`, "Tests @tester can write red", 23c, DC9:
two pieces in one cell; a partition line through a lake ring (here a hole in
the domain, the ring the partition must not lose); a road along a partition
line (one edge, both bits); a cell corner inside the domain (four pieces at
one node, which the noder must return exactly); a sliver piece with no node.
Then SP3 as 23b's record hands it to 23c: on the projected-mosaic path, two
pieces whose windows meet different tiles around one seam get the same seam
points, read from the files.

Every case also passes both of section 3D's oracles and conformity over the
union (`pieces_fixtures.problems`), because each is a path the run takes.

Where a case needs a partition line's position, it is read from a first cut
run of a domain with the same bounding box (the window, and so the partition,
is a function of the box): nothing here re-derives the window rule.

Piece ids `(j, i, k)`: `j` is the cell row from the top (rows grow down,
`y = -R h`), `i` the column from the left.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Polygon

import pieces_fixtures as pf
from feature_fixtures import Feat, write_geojson
from mosaic_fixtures import meta, quadrants
from pieces_fixtures import TIE_X, TIE_Y, rel
from test_cli_mesh_mosaic import write_tiles
from tin_engine.io.models import DemTile

TOL = "1"


@pytest.fixture(scope="module")
def dem(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, np.ndarray]:
    array = pf.hills()
    return pf.write_dem(tmp_path_factory.mktemp("dem") / "dem.tif", array), array


def args(dem: Path, domain: Path, *extra: str) -> tuple[str, ...]:
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", TOL, *extra)


def learn_lines(where: Path, dem: Path) -> tuple[float, float]:
    """The vertical and horizontal partition line of SQUARE's 2 x 2 cut."""
    first = pf.cut(
        where, *args(dem, pf.write_domain(where / "sq.geojson", pf.SQUARE)), "--pieces", "4"
    )
    xs, ys = pf.seam_lines(first.pieces)
    assert len(xs) == 1 and len(ys) == 1, (xs, ys)
    return xs.pop(), ys.pop()


# ------------------------------------------------------------------ the cases


def test_two_pieces_in_one_cell(tmp_path: Path, dem: tuple[Path, np.ndarray]) -> None:
    """A U whose arms enter the top cell apart: two components, two pieces,
    `k` 0 and 1 in that cell; the base below is one."""
    tif, array = dem
    u_ring = [
        rel(30.3, -370.7),
        rel(370.9, -369.3),
        rel(369.7, -30.1),
        rel(250.3, -29.9),
        rel(249.1, -250.7),
        rel(150.9, -249.3),
        rel(151.3, -31.1),
        rel(31.7, -29.3),
    ]
    domain = pf.write_domain(tmp_path / "u.geojson", u_ring)
    done = pf.cut(tmp_path / "run", *args(tif, domain), "--pieces", "2")
    assert (done.index["partition"]["nx"], done.index["partition"]["ny"]) == (1, 2)
    assert sorted(p.id for p in done.pieces) == [(0, 0, 0), (0, 0, 1), (1, 0, 0)]
    arms = [p for p in done.pieces if p.id[:2] == (0, 0)]
    keys = [{(float(x), float(y)) for x, y, _ in arm.points} for arm in arms]
    assert not keys[0] & keys[1], "the two arms share a vertex"
    assert pf.problems(done, array, Polygon(u_ring), 1.0) == []


def test_a_partition_line_through_a_hole(tmp_path: Path, dem: tuple[Path, np.ndarray]) -> None:
    tif, array = dem
    hole = [rel(180.3, -150.7), rel(239.1, -150.1), rel(240.9, -249.1), rel(181.7, -251.3)]
    domain = pf.write_domain(tmp_path / "holed.geojson", pf.SQUARE, (hole,))
    done = pf.cut(tmp_path / "run", *args(tif, domain), "--pieces", "4")
    xs, ys = pf.seam_lines(done.pieces)
    ring = Polygon(hole)
    assert any(ring.intersects(LineString([(x, TIE_Y), (x, TIE_Y - 400)])) for x in xs)
    assert any(ring.intersects(LineString([(TIE_X, y), (TIE_X + 400, y)])) for y in ys)
    assert sorted(p.id for p in done.pieces) == [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)]
    assert pf.problems(done, array, Polygon(pf.SQUARE, [hole]), 1.0) == []
    # The hole stays a hole: no union triangle meets its interior.
    u = pf.union(done.pieces)
    centroids = u.points[u.triangles, :2].mean(axis=1)
    assert not shapely.contains_xy(ring, centroids[:, 0], centroids[:, 1]).any()


def test_a_road_along_a_partition_line_is_one_edge_with_both_bits(
    tmp_path: Path, dem: tuple[Path, np.ndarray]
) -> None:
    tif, array = dem
    x, _ = learn_lines(tmp_path / "first", tif)
    y0, y1 = TIE_Y - 100.3, TIE_Y - 300.7
    roads = write_geojson(
        tmp_path / "road.geojson", [Feat("r", LineString([(x, y0), (x, y1)]), {"property": "road"})]
    )
    domain = pf.write_domain(tmp_path / "sq.geojson", pf.SQUARE)
    done = pf.cut(tmp_path / "run", *args(tif, domain), "--pieces", "4", "--features", str(roads))
    assert pf.seam_lines(done.pieces)[0] == {x}, "the road moved the partition"
    on_road = 0
    for piece in done.pieces:
        road = pf.bit_named(piece.vtk, "road")
        pairs = [
            tuple(sorted(((float(piece.points[a, 0]), float(piece.points[a, 1])),
                          (float(piece.points[b, 0]), float(piece.points[b, 1])))))
            for a, b in piece.edges
        ]  # fmt: skip
        assert len(set(pairs)) == len(pairs), f"a duplicated edge in {piece.id}"
        for ((xa, ya), (xb, yb)), m in zip(pairs, piece.masks, strict=True):
            if xa == xb == x and y1 <= min(ya, yb) and max(ya, yb) <= y0:
                on_road += 1
                assert m & piece.seam and m & road, (piece.id, (xa, ya), (xb, yb), m)
    assert on_road >= 2, "the road's stretch of the seam is in no piece file"
    assert pf.problems(done, array, Polygon(pf.SQUARE), 1.0) == []


def test_four_pieces_meet_at_one_corner_node(tmp_path: Path, dem: tuple[Path, np.ndarray]) -> None:
    """The cell corner inside the domain is a vertex of all four pieces, at
    the crossing of the two lattice lines exactly, with one z."""
    tif, array = dem
    domain = pf.write_domain(tmp_path / "sq.geojson", pf.SQUARE)
    done = pf.cut(tmp_path / "run", *args(tif, domain), "--pieces", "4")
    xs, ys = pf.seam_lines(done.pieces)
    (x,), (y,) = xs, ys
    zs = []
    for piece in done.pieces:
        at = np.flatnonzero((piece.points[:, 0] == x) & (piece.points[:, 1] == y))
        assert len(at) == 1, (piece.id, at)
        zs.append(np.float64(piece.points[at[0], 2]).tobytes())
    assert len(done.pieces) == 4 and len(set(zs)) == 1
    col, row = (x - TIE_X) / pf.H, (TIE_Y - y) / pf.H
    assert col.is_integer() and row.is_integer()
    assert np.float64(array[int(row), int(col)]).tobytes() == zs[0]
    assert pf.problems(done, array, Polygon(pf.SQUARE), 1.0) == []


def test_a_sliver_piece_with_no_dem_node(tmp_path: Path, dem: tuple[Path, np.ndarray]) -> None:
    """SQUARE's box with its upper-right block cut away but for an L 4.1 m
    wide along the vertical line and 3.3 m along the horizontal one: the top
    right cell's piece holds no DEM node off the seams. The seam pass still
    gives its boundary heights; the run succeeds and the union is sound."""
    tif, array = dem
    x, y = learn_lines(tmp_path / "first", tif)
    tl, bl, br = pf.SQUARE[0], pf.SQUARE[1], pf.SQUARE[2]
    ring = [
        tl,
        bl,
        br,
        (TIE_X + 371.9, y + 3.3),  # the box's right side, as SQUARE's
        (x + 4.1, y + 3.3),
        (x + 4.1, TIE_Y - 30.3),
    ]
    polygon = Polygon(ring)
    assert polygon.bounds == Polygon(pf.SQUARE).bounds
    sliver = polygon.intersection(shapely.box(x, y, TIE_X + 1000, TIE_Y + 1000))
    nx, ny, _ = pf.node_xy(array)
    inside = shapely.contains_xy(sliver, nx, ny)
    assert sliver.area > 0 and not inside.any(), "the sliver holds a DEM node"

    domain = pf.write_domain(tmp_path / "sliver.geojson", ring)
    done = pf.cut(tmp_path / "run", *args(tif, domain), "--pieces", "4")
    assert pf.seam_lines(done.pieces) == ({x}, {y})
    top_right = [p for p in done.pieces if p.id[:2] == (0, 1)]
    assert len(top_right) == 1 and len(top_right[0].triangles) > 0
    assert pf.problems(done, array, polygon, 1.0) == []


# ------------------------------------------------------------------ SP3 in 23c


def test_two_pieces_reading_different_tiles_agree_on_their_seam(tmp_path: Path) -> None:
    """The same DEM as four tiles, cut two nodes past where the partition
    lines fall (column and row 21 of this window), overlapping by one line:
    a piece left of a seam can be served by the west tiles alone, the piece
    right of it needs both. Both still write the same seam."""
    array = pf.hills()
    grid = meta(rows=pf.N, cols=pf.N, dx=pf.H, dy=pf.H)
    whole = DemTile(meta=grid, array=array)
    tiles = tmp_path / "tiles"
    write_tiles(tiles, quadrants(whole, row_cut=23, col_cut=23, overlap=1))
    domain = pf.write_domain(tmp_path / "sq.geojson", pf.SQUARE)
    done = pf.cut(tmp_path / "run", *args(tiles, domain), "--pieces", "4")
    xs, ys = pf.seam_lines(done.pieces)
    assert len(xs) == 1 and len(ys) == 1
    assert pf.problems(done, array, Polygon(pf.SQUARE), 1.0) == []
