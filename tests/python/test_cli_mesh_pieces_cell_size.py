"""A cut on a cell size binary cannot hold (9.7 m, 9.9 m): CS1-CS4, CS7-CS9.

`docs/increments/23-basin-scale.md`, "23c-4: cutting a DEM whose cell size is
not exact in binary", "Red tests". A seam end's place on the lattice comes
from the partition (`basin_run.line_rc`, E1), the seam strip's box from that
place (`basin_run.seam_box`, E2), a corner's height from its node
(`basin_run.corner_heights`, E3), and a seam-pass point on an axis-parallel
seam edge takes the edge's across coordinate (E4). CS5 (the bytes kept) is
in `test_cli_mesh_pieces_cell_size_bytes.py`, CS6 in `test_pieces_slice.py`.

The scenes: `pieces_fixtures.hills()` written at the spacing under the
`SQUARE` domain, `--tolerance 1`. At 9.7 m and 9.9 m with `--pieces 4` the
run window starts at lattice node (2, 2) and the partition is 2 x 2 cells of
19 nodes, so the cut lines are run-window column and row 19 (lattice 21):
`500203.7` / `6599796.3` at 9.7 m, `500207.9` / `6599792.1` at 9.9 m.

Distances: coordinates are UTM metres (x 5e5, y 6.6e6, an ulp about 1e-10
and 1e-9 m); `NEAR` (1e-6 m) is the corner-finding bound, a thousandth of
the 1 mm snap, checked on these 400 m scenes, where a seam end lies 9.3e-10 m
from its line at most (the design's prototype).
"""

from __future__ import annotations

import dataclasses
import math
import re
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from shapely.geometry import Polygon

import pieces_fixtures as pf
from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff
from tin_engine import basin_run
from tin_engine.decompose import Partition
from tin_engine.io.models import DemTile, RasterMeta

SNAP = 1e-3  # m: `--snap-spacing`'s default
NEAR = 1e-6  # m; see the module docstring for the scale


# ------------------------------------------------------------------ the partition, by hand


def window_meta(h: float, n: int, nodata: float | None = None) -> RasterMeta:
    """The run window's meta as the run builds it: lattice node (2, 2) on."""
    return RasterMeta(
        x_min=TIE_X + 2 * h, y_max=TIE_Y - 2 * h, delta_x=h, delta_y=h, cols=n, rows=n,
        epsg=25833, nodata=nodata, nodata_source="tag" if nodata is not None else "absent",
        pixel_is_area=False, vertical_unit_assumed=False,
    )  # fmt: skip


def partitioned(h: float, n: int, cells: int, size: int) -> basin_run.Partitioned:
    """`cells` x `cells` cells of `size` nodes over an `n` x `n` window."""
    ask = basin_run.CutAsk(
        pieces=cells * cells, budget=1 << 34, origin=(TIE_X, TIE_Y), row0=2, col0=2
    )
    return basin_run.cut_partition(
        window_meta(h, n), 1.0, Partition(nx=cells, ny=cells, dx=size, dy=size), ask
    )


@pytest.fixture(scope="module")
def at_97() -> basin_run.Partitioned:
    """The 9.7 m scene's partition with four pieces (the run's own, measured)."""
    return partitioned(9.7, 38, 2, 19)


@pytest.fixture(scope="module")
def at_99() -> basin_run.Partitioned:
    return partitioned(9.9, 37, 2, 19)


def today_box(cut: Any, snap: float, a: Any, b: Any) -> tuple[int, int, int, int]:
    """The strip box as the seam pass builds it before 23c-4: `floor` and
    `ceil` of each end's fractional index in the run window, minus and plus one."""
    m = cut.meta
    (ax, ay), (bx, by) = a, b
    c0, c1 = ((v - m.x_min) / m.delta_x for v in (min(ax, bx), max(ax, bx)))
    r0, r1 = ((m.y_max - v) / m.delta_y for v in (max(ay, by), min(ay, by)))
    return math.floor(r0) - 1, math.ceil(r1) + 1, math.floor(c0) - 1, math.ceil(c1) + 1


def test_the_hand_built_partitions_are_the_runs(
    at_97: basin_run.Partitioned, at_99: basin_run.Partitioned
) -> None:
    """The unit tests below stand on the design's measured numbers: the lines,
    and the divisions that miss the whole number (so the probes can fail)."""
    assert (at_97.xs, at_97.ys) == ([500203.7], [6599796.3])
    assert (at_99.xs, at_99.ys) == ([500207.9], [6599792.1])
    assert (at_97.xs[0] - at_97.meta.x_min) / 9.7 == 18.999999999998803
    assert (at_99.xs[0] - at_99.meta.x_min) / 9.9 == 19.000000000003528


# ------------------------------------------------------------------ E1: line_rc


def test_an_end_at_the_crossing_gets_both_indices(at_97: basin_run.Partitioned) -> None:
    assert basin_run.line_rc(at_97, SNAP, at_97.xs[0], at_97.ys[0]) == (19, 19)


def test_an_end_the_noder_moved_one_ulp_still_gets_the_lines_index(
    at_99: basin_run.Partitioned,
) -> None:
    """At 9.9 m the noder returns the horizontal line as `6599792.100000001`."""
    moved = float(np.nextafter(at_99.ys[0], math.inf))
    assert moved == 6599792.100000001
    assert basin_run.line_rc(at_99, SNAP, at_99.xs[0], moved) == (19, 19)


def test_an_end_on_one_line_gets_that_index_only(at_97: basin_run.Partitioned) -> None:
    south = TIE_Y - 370.668  # on the square's bottom edge, on no horizontal line
    assert basin_run.line_rc(at_97, SNAP, at_97.xs[0], south) == (None, 19)
    west = TIE_X + 30.0
    assert basin_run.line_rc(at_97, SNAP, west, at_97.ys[0]) == (19, None)
    assert basin_run.line_rc(at_97, SNAP, west, south) == (None, None)


def test_a_later_line_gets_its_number_times_the_cell() -> None:
    """Line k (1-based) is run-window column `k dx`: 10, 20, 30 for cells of 10."""
    cut = partitioned(9.7, 38, 4, 10)
    assert len(cut.xs) == 3
    got = [basin_run.line_rc(cut, SNAP, x, y) for x, y in zip(cut.xs, cut.ys, strict=True)]
    assert got == [(10, 10), (20, 20), (30, 30)]


def test_the_bound_is_three_quarters_of_the_snap(at_97: basin_run.Partitioned) -> None:
    """At the 1 mm snap the bound is 0.75 mm: 0.7 mm off is on, 0.8 mm is not."""
    x, south = at_97.xs[0], TIE_Y - 370.668
    assert basin_run.line_rc(at_97, SNAP, x + 0.7e-3, south) == (None, 19)
    assert basin_run.line_rc(at_97, SNAP, x - 0.7e-3, south) == (None, 19)
    assert basin_run.line_rc(at_97, SNAP, x + 0.8e-3, south) == (None, None)


def test_the_bound_shrinks_to_a_quarter_cell_under_a_large_snap(
    at_97: basin_run.Partitioned,
) -> None:
    """`--snap-spacing` has no upper limit, so the bound is `min(0.75 snap,
    h / 4)`: at a 10 m snap on 9.7 m cells it is 2.425 m, not 7.5 m."""
    x, south = at_97.xs[0], TIE_Y - 370.668
    assert basin_run.line_rc(at_97, 10.0, x + 2.4, south) == (None, 19)
    assert basin_run.line_rc(at_97, 10.0, x + 2.5, south) == (None, None)
    assert basin_run.line_rc(at_97, 10.0, x - 5.0, south) == (None, None)


# ------------------------------------------------------------------ E2: seam_box


def test_the_box_of_a_seam_down_from_the_crossing_holds_the_line(
    at_97: basin_run.Partitioned,
) -> None:
    """9.7 m: `floor` of 18.999999999998803 starts the box at column 17, one
    node west of the eastern piece's window (columns 18 on); E2 starts it at
    the line minus one. The south end is on no line: its fractional row is
    (6599980.6 - 6599629.332) / 9.7 = 36.2, so the box ends at row 38."""
    a, b = (at_97.xs[0], at_97.ys[0]), (at_97.xs[0], TIE_Y - 370.668)
    assert today_box(at_97, SNAP, a, b) == (17, 38, 17, 20)
    assert basin_run.seam_box(at_97, SNAP, a, b) == (18, 38, 18, 20)
    assert basin_run.seam_box(at_97, SNAP, b, a) == (18, 38, 18, 20)


def test_the_box_of_a_seam_into_a_moved_crossing_holds_the_line(
    at_99: basin_run.Partitioned,
) -> None:
    """9.9 m: the crossing's column index 19.000000000003528 `ceil`s to 20,
    and plus one is column 21, past piece (0, 0)'s window (to column 20)."""
    corner = (at_99.xs[0], float(np.nextafter(at_99.ys[0], math.inf)))
    west = (TIE_X + 30.0, corner[1])  # column (30.0 - 19.8) / 9.9 = 1.03
    assert today_box(at_99, SNAP, west, corner) == (17, 20, 0, 21)
    assert basin_run.seam_box(at_99, SNAP, west, corner) == (18, 20, 0, 20)


@pytest.mark.parametrize(
    "a, b",
    [
        (pf.rel(210.0, -210.0), pf.rel(210.0, -370.5)),
        (pf.rel(29.9, -210.0), pf.rel(210.0, -210.0)),
        (pf.rel(33.3, -31.0), pf.rel(41.7, -58.2)),
    ],
    ids=["down-from-the-crossing", "into-the-crossing", "on-no-line"],
)
def test_on_a_whole_metre_lattice_the_box_is_todays(a: Any, b: Any) -> None:
    """Where the fractional indices are exact integers the box is unchanged."""
    cut = partitioned(10.0, 38, 2, 19)
    assert (cut.xs, cut.ys) == ([TIE_X + 210.0], [TIE_Y - 210.0])
    assert basin_run.seam_box(cut, SNAP, a, b) == today_box(cut, SNAP, a, b)


# ------------------------------------------------------------------ E3: corner_heights (CS9)


def loaded(array: np.ndarray, meta: RasterMeta, r0: int, c0: int, rows: int, cols: int) -> Any:
    """A piece's DEM as `load` hands it: its tile, from run-window node (r0, c0)."""
    tile_meta = meta.model_copy(update={
        "x_min": meta.x_min + c0 * meta.delta_x, "y_max": meta.y_max - r0 * meta.delta_y,
        "rows": rows, "cols": cols})  # fmt: skip
    tile = DemTile(meta=tile_meta, array=array[r0 : r0 + rows, c0 : c0 + cols])
    return SimpleNamespace(tile=tile, row0=r0, col0=c0, grid=None, checks=None)


#: The ends CS9 asks about: the crossing, the crossing one ulp off in y, an
#: end on the vertical line only, an end on no line.
def ends_at(cut: basin_run.Partitioned) -> np.ndarray:
    x, y = cut.xs[0], cut.ys[0]
    return np.array(
        [(x, y), (x, float(np.nextafter(y, math.inf))), (x, TIE_Y - 370.668), (TIE_X + 30.0, y)]
    )


def run_array(n: int, at: tuple[int, int] | None = None, value: float = 0.0) -> np.ndarray:
    array = (50.0 + np.arange(n * n, dtype=np.float32).reshape(n, n) / 8).astype(np.float32)
    if at is not None:
        array[at] = value
    return array


def test_a_corner_takes_its_nodes_value_in_every_piece(at_97: basin_run.Partitioned) -> None:
    """The corner is run-window node (19, 19): in piece (0, 0)'s tile (from
    node 0) it is element (19, 19), in piece (1, 1)'s (from node 18) (1, 1);
    both give that node's value as a float, the ulp-off end too, and an end
    on one line or none is not a corner."""
    array, xy = run_array(38), ends_at(at_97)
    north_west = loaded(array, at_97.meta, 0, 0, 21, 21)
    south_east = loaded(array, at_97.meta, 18, 18, 20, 20)
    want = {0: float(array[19, 19]), 1: float(array[19, 19])}
    for dem in (north_west, south_east):
        got = basin_run.corner_heights(dem, at_97, SNAP, xy, [0, 1, 2, 3])
        assert got == want
        assert all(type(v) is float for v in got.values()), got


@pytest.mark.parametrize("nodata", [-32767.0, None], ids=["sentinel", "nan"])
def test_a_corner_on_a_nodata_node_has_no_height(
    at_97: basin_run.Partitioned, nodata: float | None
) -> None:
    value = math.nan if nodata is None else nodata
    array = run_array(38, (19, 19), value)
    meta = window_meta(9.7, 38, nodata)
    dem = loaded(array, meta, 18, 18, 20, 20)
    cut = dataclasses.replace(at_97, meta=meta)  # the run's meta and the tile's agree
    got = basin_run.corner_heights(dem, cut, SNAP, ends_at(at_97), [0, 1, 2, 3])
    assert got == {0: None, 1: None}


def test_only_the_ends_asked_for_are_answered(at_97: basin_run.Partitioned) -> None:
    dem = loaded(run_array(38), at_97.meta, 18, 18, 20, 20)
    assert basin_run.corner_heights(dem, at_97, SNAP, ends_at(at_97), [1, 2]).keys() == {1}
    assert basin_run.corner_heights(dem, at_97, SNAP, ends_at(at_97), []) == {}


# ------------------------------------------------------------------ the runs


def scene_args(where: Path, h: float, pieces: str, ring: pf.Ring = pf.SQUARE) -> tuple[str, ...]:
    dem = pf.write_dem(where / "dem.tif", pf.hills(), h=h)
    domain = pf.write_domain(where / "domain.geojson", ring)
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1", "--pieces", pieces)


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> Any:
    """Each square-scene run at `(h, pieces)`, made once for the module."""
    made: dict[tuple[float, str], pf.Cut] = {}

    def get(h: float, pieces: str) -> pf.Cut:
        if (h, pieces) not in made:
            where = tmp_path_factory.mktemp(f"cell-{h}-{pieces}")
            made[h, pieces] = pf.run(where / "run", *scene_args(where, h, pieces))
        return made[h, pieces]

    return get


@pytest.fixture(params=[(9.7, "4"), (9.9, "4"), (9.7, "16")], ids=["9.7m-4", "9.9m-4", "9.7m-16"])
def square_run(request: pytest.FixtureRequest, runs: Any) -> Any:
    """CS1 (9.7 m, 4 pieces: the box undershoots), CS2 (9.9 m, 4: the noder
    moves a line and the box overshoots) and CS3 (9.7 m, 16: nine corners)."""
    h, pieces = request.param
    return SimpleNamespace(h=h, pieces=pieces, done=runs(h, pieces))


def read_cut(done: pf.Cut) -> pf.Cut:
    """`pieces_fixtures.cut`'s reading of a run already made."""
    assert done.code == 0, done.output
    assert done.dir is not None, f"no <out>.pieces directory beside {done.out}"
    return pf.read_index_and_pieces(done)


def test_cs1_cs3_a_cut_on_a_decimal_cell_is_whole_and_correct(square_run: Any) -> None:
    """Exit 0, an index, and every oracle at the run's spacing: conformity,
    one z per (x, y), a valid union, no overlap, constrained Delaunay and the
    tolerance, with seams judged within `NEAR` of a lattice line."""
    done = read_cut(square_run.done)
    assert len(done.pieces) == int(square_run.pieces)
    found = pf.problems(done, pf.hills(), Polygon(pf.SQUARE), 1.0, square_run.h, near=True)
    assert found == []


def holders_near(pieces: list[pf.Piece], x: float, y: float) -> dict[pf.PieceId, set[pf.XYZ]]:
    """Each piece file's vertices within `NEAR` of `(x, y)` in both axes."""
    out: dict[pf.PieceId, set[pf.XYZ]] = defaultdict(set)
    for piece in pieces:
        close = (np.abs(piece.points[:, 0] - x) <= NEAR) & (np.abs(piece.points[:, 1] - y) <= NEAR)
        for p in piece.points[close]:
            out[piece.id].add(pf.xyz(p))
    return dict(out)


def lines(index: dict[str, Any]) -> tuple[list[float], list[float]]:
    """The cut lines' x and y from the index's partition record."""
    part = index["partition"]
    (ox, oy), h, r0, c0 = part["origin"], part["spacing"], part["row0"], part["col0"]
    xs = [ox + (c0 + i * part["dx"]) * h for i in range(1, part["nx"])]
    ys = [oy - (r0 + j * part["dy"]) * h for j in range(1, part["ny"])]
    return xs, ys


def test_cs3_each_corner_is_written_with_one_z_by_its_four_pieces(runs: Any) -> None:
    """E3, at 9.7 m with 16 pieces: at a crossing of two cut lines each of the
    four pieces holds two of its four seam edges; all four write the corner
    with one (x, y, z)."""
    done = read_cut(runs(9.7, "16"))
    xs, ys = lines(done.index)
    assert (len(xs), len(ys)) == (3, 3)
    for x in xs:
        for y in ys:
            held = holders_near(done.pieces, x, y)
            assert len(held) == 4, (x, y, held)
            assert len(set().union(*held.values())) == 1, (x, y, held)


def test_cs4_every_seam_edge_is_straight_bit_for_bit(square_run: Any) -> None:
    """E4: a seam-pass point on an axis-parallel seam edge takes the edge's
    across coordinate, so every seam edge has equal x or equal y at both ends."""
    done = read_cut(square_run.done)
    bent = [
        (a, b)
        for piece in done.pieces
        for a, b in pf.seam_edges(piece)
        if a[0] != b[0] and a[1] != b[1]
    ]
    assert bent == []


#: CS7's ring, in metres from the tie point: the square's extent with two
#: slots cut from the top edge, so a spike points up to the crossing of the
#: two cut lines at 9.7 m, (203.7, -203.7); only the vertical seam below
#: leaves it.
SPIKE: pf.Ring = [pf.rel(x, y) for x, y in [
    (27.1, -28.7), (27.1, -372.3), (371.9, -372.3), (371.9, -28.7), (263.7, -28.7),
    (223.7, -283.7), (203.7, -203.7), (183.7, -283.7), (143.7, -28.7)]]  # fmt: skip


def test_cs7_a_corner_only_one_seam_leaves(tmp_path: Path) -> None:
    """E1 takes every line an end lies on, whether or not a seam leaves the
    end along it: the apex gets the horizontal line's row too, so its box
    stays inside both pieces, and both write it with one z."""
    done = read_cut(pf.run(tmp_path / "run", *scene_args(tmp_path, 9.7, "4", SPIKE)))
    found = pf.problems(done, pf.hills(), Polygon(SPIKE), 1.0, 9.7, near=True)
    assert found == []
    held = holders_near(done.pieces, *pf.rel(203.7, -203.7))
    assert sorted(held) == [(1, 0, 0), (1, 1, 0)], held
    assert len(set().union(*held.values())) == 1, held


def test_cs8_nodata_beside_a_seam_at_9_9_m_still_refuses(tmp_path: Path) -> None:
    """The limit the step leaves: the two pieces' tiles are one rounding apart,
    so a seam vertex beside a NoData node has data in one piece's frame and
    none in its neighbour's; K4 refuses (exit 1, no index) rather than write
    a mesh with zero-area triangles. Flips when the follow-up step lands."""
    array = pf.hills()
    array[10, 22] = -32767.0  # lattice node (row 10, column 22): beside the vertical seam
    dem = tmp_path / "dem.tif"
    dem.write_bytes(micro_tiff(array, scale=(9.9, 9.9, 0.0), nodata="-32767").getvalue())
    domain = pf.write_domain(tmp_path / "domain.geojson", pf.SQUARE)
    args = ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1", "--pieces", "4")
    done = pf.run(tmp_path / "run", *args)
    said = " ".join(done.output.split())
    assert done.code == 1, said
    assert "seam edge 7 differs between pieces (0, 0, 0) and (0, 1, 0)" in said
    assert done.dir is not None and not (done.dir / "index.json").exists()
    assert re.search(r"\bno index\b", said), said
