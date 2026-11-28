"""`rasputin mesh --domain PATH --frozen-holes`: a domain's hole rings kept as given.

`docs/increments/198-frozen-holes.md`, *5. Tests*. With the flag no vertex is
added on a hole ring of ``--domain`` and no hole edge is split, so the mesh can
be joined to another along the hole; the DEM's error along the hole edges is
recorded (``frozen_edges_max_error_m``), not refined. A constraint line that
meets a hole edge between its vertices, and a hole vertex with no DEM under
it, are refused.

The fixture is the design's: a 10 m DEM of 101 x 101 nodes over a local
0 to 1,000 m square, z = (x - 500)^2 / 1,000 m, written with `micro_tiff`'s
tie point, so local (x, y) is the file's (TIE_X + x, TIE_Y - 1,000 + y). The
domain is the whole square with the hole 400 to 600 m, whose bottom and top
edges carry vertices at x = 450.0004, 500.0004 and 550.0004 m (0.4 mm off the
1 mm lattice): 10 hole vertices.

A second domain adds a hole with slanted edges, the diamond with corners
(440, 800), (500, 720), (560, 800) and (500, 880), every corner a DEM node.
Its figures, by hand:

* An edge runs 60 m in x and 80 m in y, 100 m long. The DEM is bilinear
  between nodes and z depends on x alone, so along the edge the DEM is the
  piecewise linear interpolant of the parabola in x, with breaks at multiples
  of 10 m. The ends are nodes, so the chord from end to end is exact there;
  chord minus DEM at a break x is (x - x_a)(x_b - x) / 1,000 m, here 0.5, 0.8,
  0.9, 0.8, 0.5 m at 10 to 50 m from the left end: largest 0.9 m, at the
  edge's middle, (470, 840) or (530, 840), itself a DEM node.
* ``frozen_edges_max_error_m`` is exact for the bilinear DEM on a straight
  edge (design section 3.6), so the figure is 0.9 m (the square's largest is
  0.6 m). The ends' z are at most 3.6 m, so float32 storage moves the figure
  by under 1e-6.
* DEM nodes: local x = 500 + 10 i, y = 800 + 10 j lies on the ring where
  4|i| + 3|j| = 24 (8 nodes: the corners and the four middles) and strictly
  inside where 4|i| + 3|j| < 24 (15 + 2 (13 + 11 + 7 + 5 + 3) = 93 nodes).
  With the square's 80 on its ring and 361 inside, 10,201 - 80 - 361 - 8 - 93
  = 9,659 nodes lie off both rings. Every other node is at least 2 m from a
  diamond edge (a lattice point's distance to 4x + 3y = c, in 10 m units, is
  a multiple of 1/5).

The diamond's corners are on the 1 mm lattice, so they do not move.

A third domain has one hole half a cell off the DEM's columns, corners
(205, 100), (395, 100), (395, 300), (205, 300) (design section 5, *F1 exact
error figure*). The corners' bilinear z is 87.05 and 11.05 m, so the bottom
edge's chord is 87.05 - 0.4 (x - 205); at the node x = 300 the DEM is 40 and
the chord 49.05: error 9.05 m, the largest. Sampling at cell steps from
x = 205 finds 295 and 305 instead, and reads 9.0000005.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import Point, Polygon

from cli_helpers import USAGE, UTM33, Ring, geojson, help_text, invoke, squashed
from feature_fixtures import Feat, write_geojson
from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff
from pieces_fixtures import delaunay_violations
from recordread import file_field
from vtkread import VtkFile, lines_as_array, polygons_as_array, read_vtk

X0, Y0 = TIE_X, TIE_Y - 1000.0  # local (0, 0) in the files' coordinates
OFF = 0.0004  # how far the hole's middle vertices lie off the 1 mm lattice
SNAP = 1e-3  # the default --snap-spacing
#: The farthest the noder's snap moves a vertex: half a lattice diagonal.
SNAP_MOVE = SNAP * math.sqrt(2) / 2
#: How far a written vertex may lie from an input ring and count as on it.
ON_RING = 2e-3
TOLERANCE = 0.1
BIT = 1 << 11  # the internal frozen-hole bit, never written
NODATA = -9999.0
#: The three entries the flag adds to the record.
ENTRIES = ("frozen_holes", "frozen_edges_max_error_m", "hole_vertices_snapped_max_m")


def local(x: float, y: float) -> tuple[float, float]:
    return (X0 + x, Y0 + y)


OUTER: Ring = [local(0, 0), local(1000, 0), local(1000, 1000), local(0, 1000)]
SQUARE_HOLE: Ring = [
    local(400, 400),
    local(450 + OFF, 400),
    local(500 + OFF, 400),
    local(550 + OFF, 400),
    local(600, 400),
    local(600, 600),
    local(550 + OFF, 600),
    local(500 + OFF, 600),
    local(450 + OFF, 600),
    local(400, 600),
]
DIAMOND: Ring = [local(440, 800), local(500, 720), local(560, 800), local(500, 880)]
#: A hole whose x corners lie half a cell off the DEM's columns.
OFF_GRID_HOLE: Ring = [local(205, 100), local(395, 100), local(395, 300), local(205, 300)]
#: The square hole with an 11th vertex 0.2 mm from (450.0004, 400): both snap to (450.000, 400).
MERGING_HOLE: Ring = [SQUARE_HOLE[0], local(450.0002, 400), *SQUARE_HOLE[1:]]


def heights() -> np.ndarray:
    """101 x 101 rows of z = (x - 500)^2 / 1,000 m, x = 10 col: exact at the nodes up to float32."""
    x = np.arange(101) * 10.0
    return np.tile(((x - 500.0) ** 2 / 1000.0).astype(np.float32), (101, 1))


def with_nodata_at_two_hole_corners() -> np.ndarray:
    """`heights` with nodata at local (400, 400) and (400, 600): column 40, rows 60 and 40."""
    z = heights()
    z[60, 40] = z[40, 40] = NODATA
    return z


def write_line(path: Path, x: float) -> Path:
    """A ``--lines`` file: one line at local x from y = 0 to 1,000 m."""
    doc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": UTM33}},
        "features": [
            {
                "type": "Feature",
                "properties": {},
                "geometry": {"type": "LineString", "coordinates": [local(x, 0), local(x, 1000)]},
            }
        ],
    }
    path.write_text(json.dumps(doc))
    return path


# ---------------------------------------------------------------- the runs


class Inputs:
    """Every input file of the suite, written once."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.dem = root / "dem.tif"
        self.dem.write_bytes(micro_tiff(heights(), scale=(10.0, 10.0, 0.0)).getvalue())
        self.nodata_dem = root / "nodata.tif"
        stream = micro_tiff(
            with_nodata_at_two_hole_corners(), scale=(10.0, 10.0, 0.0), nodata="-9999"
        )
        self.nodata_dem.write_bytes(stream.getvalue())
        self.square = geojson(root / "square.geojson", OUTER, (SQUARE_HOLE,))
        self.two_holes = geojson(root / "two-holes.geojson", OUTER, (SQUARE_HOLE, DIAMOND))
        self.no_hole = geojson(root / "no-hole.geojson", OUTER)
        self.off_grid = geojson(root / "off-grid.geojson", OUTER, (OFF_GRID_HOLE,))
        self.merging = geojson(root / "merging.geojson", OUTER, (MERGING_HOLE,))
        unit = Feat(1, Polygon(OUTER), {"name": "u1"})
        self.units = write_geojson(root / "units.geojson", [unit], UTM33)
        self.through_vertices = write_line(root / "through-vertices.geojson", 500.0004)
        self.between_vertices = write_line(root / "between-vertices.geojson", 512.3)


class Run:
    """One `rasputin mesh` run: exit code, output, and the files it wrote."""

    def __init__(self, root: Path, name: str, *args: str, record: bool = True) -> None:
        """With `record` false no ``--record`` is asked: a cut run refuses it on its own."""
        self.out = root / f"{name}.vtk"
        self.record_path = root / f"{name}.json"
        asked = ("--record", str(self.record_path)) if record else ()
        self.code, self.output = invoke("mesh", *args, "--out", str(self.out), *asked)

    def ok(self) -> Run:
        assert self.code == 0, self.output
        return self

    @property
    def vtk(self) -> VtkFile:
        return read_vtk(self.ok().out.read_bytes())

    @property
    def record(self) -> dict[str, object]:
        doc = json.loads(self.ok().record_path.read_text(encoding="utf-8"))
        assert isinstance(doc, dict)
        return doc


@pytest.fixture(scope="module")
def inputs(tmp_path_factory: pytest.TempPathFactory) -> Inputs:
    return Inputs(tmp_path_factory.mktemp("frozen-holes"))


def mesh_args(dem: Path, domain: Path, *extra: str) -> tuple[str, ...]:
    return ("--dem", str(dem), "--domain", str(domain), *extra)


@pytest.fixture(scope="module")
def frozen(inputs: Inputs) -> Run:
    """F1's run: the square hole, --tolerance 0.1, --frozen-holes."""
    args = mesh_args(inputs.dem, inputs.square, "--tolerance", str(TOLERANCE), "--frozen-holes")
    return Run(inputs.root, "frozen", *args)


@pytest.fixture(scope="module")
def plain(inputs: Inputs) -> Run:
    """The same run without the flag."""
    args = mesh_args(inputs.dem, inputs.square, "--tolerance", str(TOLERANCE))
    return Run(inputs.root, "plain", *args)


@pytest.fixture(scope="module")
def frozen_two(inputs: Inputs) -> Run:
    """The square hole and the diamond, --tolerance 0.1, --frozen-holes."""
    args = mesh_args(inputs.dem, inputs.two_holes, "--tolerance", str(TOLERANCE), "--frozen-holes")
    return Run(inputs.root, "frozen-two", *args)


@pytest.fixture(scope="module")
def frozen_budget(inputs: Inputs) -> Run:
    """F1's run with --triangles 3000 in place of --tolerance."""
    args = mesh_args(inputs.dem, inputs.square, "--triangles", "3000", "--frozen-holes")
    return Run(inputs.root, "frozen-budget", *args)


# ---------------------------------------------------------------- geometry, written here


def to_local(xy: np.ndarray) -> np.ndarray:
    return np.asarray(xy, dtype=np.float64)[:, :2] - [X0, Y0]


def ring_distance(xy: np.ndarray, ring: Ring) -> np.ndarray:
    """Each point's distance to the closed `ring`, in local coordinates."""
    p = to_local(xy)
    a = to_local(np.array(ring))
    b = np.roll(a, -1, axis=0)
    d = b - a
    t = np.clip(((p[:, None] - a) * d).sum(-1) / (d * d).sum(-1), 0.0, 1.0)
    return np.asarray(np.linalg.norm(p[:, None] - (a + t[..., None] * d), axis=-1).min(axis=1))


def on_ring(vtk: VtkFile, ring: Ring) -> np.ndarray:
    """The indices of the written vertices on `ring`, within the snap's reach."""
    return np.flatnonzero(ring_distance(vtk.points, ring) < ON_RING)


def matched(vtk: VtkFile, ring: Ring) -> list[int]:
    """For each input vertex of `ring`, the nearest written vertex, required within the snap."""
    p = to_local(vtk.points)
    out = []
    for x, y in to_local(np.array(ring)):
        d = np.hypot(p[:, 0] - x, p[:, 1] - y)
        k = int(np.argmin(d))
        assert d[k] <= SNAP_MOVE + 1e-9, f"no written vertex within the snap of ({x}, {y})"
        out.append(k)
    return out


def assert_ring_kept(vtk: VtkFile, ring: Ring) -> None:
    """I1: the written vertices on `ring` are its input vertices, one each, and
    each input edge is one line cell of the file."""
    found = on_ring(vtk, ring)
    nearest = matched(vtk, ring)
    assert sorted(found.tolist()) == sorted(nearest), (
        f"{len(found)} written vertices on a hole ring of {len(ring)} vertices"
    )
    cells = {tuple(sorted(map(int, e))) for e in lines_as_array(vtk)}
    for i, j in zip(nearest, nearest[1:] + nearest[:1], strict=True):
        assert (min(i, j), max(i, j)) in cells, f"hole edge {i}-{j} is not one line cell"


def dem_nodes() -> tuple[np.ndarray, np.ndarray]:
    """Every node's file (x, y) and z as stored (float32, widened)."""
    r, c = np.indices((101, 101))
    xy = np.column_stack([(X0 + 10.0 * c).ravel(), (Y0 + 1000.0 - 10.0 * r).ravel()])
    return xy, heights().astype(np.float64).ravel()


def node_classes(xy: np.ndarray, holes: tuple[Ring, ...]) -> tuple[np.ndarray, np.ndarray]:
    """(on a hole ring, strictly inside a hole) for each node."""
    on = np.zeros(len(xy), dtype=bool)
    inside = np.zeros(len(xy), dtype=bool)
    p = to_local(xy)
    for ring in holes:
        on |= ring_distance(xy, ring) < 1e-6
        shape = Polygon(to_local(np.array(ring)))
        inside |= np.array([shape.contains(Point(q)) for q in p]) & ~on
    return on, inside


def located_errors(vtk: VtkFile, xy: np.ndarray, z: np.ndarray) -> np.ndarray:
    """``|plane - z|`` at each node against a written triangle holding it (closed,
    1e-12 slack), the plane from the file's own vertices; -inf for a node in no
    triangle. Chunked, so the node-by-triangle table stays small."""
    t = polygons_as_array(vtk)
    origin = vtk.points[:, :2].min(axis=0)
    v = vtk.points[:, :2] - origin
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    za, zb, zc = (vtk.points[t[:, k], 2] for k in range(3))
    two_a = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    out = np.empty(len(z))
    for lo in range(0, len(z), 256):
        qx = xy[lo : lo + 256, 0, None] - origin[0]
        qy = xy[lo : lo + 256, 1, None] - origin[1]

        def weight(
            p: np.ndarray, q: np.ndarray, qx: np.ndarray = qx, qy: np.ndarray = qy
        ) -> np.ndarray:
            return np.asarray(
                ((q[:, 0] - p[:, 0]) * (qy - p[:, 1]) - (q[:, 1] - p[:, 1]) * (qx - p[:, 0]))
                / two_a
            )

        wa, wb, wc = weight(b, c), weight(c, a), weight(a, b)
        holds = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
        error = np.abs(wa * za + wb * zb + wc * zc - z[lo : lo + 256, None])
        out[lo : lo + 256] = np.where(holds, error, -np.inf).max(axis=1)
    return out


def constrained_delaunay_violations(vtk: VtkFile) -> list[tuple[tuple[int, int], int]]:
    """Section 3D's other oracle, exact, in the file's frame: the line cells are the constraints."""
    constrained = {tuple(sorted(map(int, e))) for e in lines_as_array(vtk)}
    return delaunay_violations(vtk.points[:, :2], polygons_as_array(vtk), constrained)


# ---------------------------------------------------------------- F1


HELP = (
    "Keep each hole ring of --domain as given: no vertex added on it and no edge of it split, "
    "so the mesh can be joined to another along the hole. The DEM's error along the hole "
    "edges is recorded, not refined."
)


def test_the_flag_and_its_help_are_in_mesh_help() -> None:
    """Section 3.1's help text, wherever Rich wrapped it."""
    flat = squashed(help_text("mesh"))
    assert squashed("--frozen-holes") in flat
    assert squashed(HELP) in flat


class TestF1Frozen:
    """F1: the hole ring kept exactly, its figures recorded."""

    def test_exactly_the_ten_hole_vertices_lie_on_the_ring(self, frozen: Run) -> None:
        vtk = frozen.vtk
        assert len(on_ring(vtk, SQUARE_HOLE)) == 10
        assert_ring_kept(vtk, SQUARE_HOLE)

    def test_the_vertices_off_the_lattice_moved_0_4_mm_in_x_and_the_rest_not_at_all(
        self, frozen: Run
    ) -> None:
        vtk = frozen.vtk
        for k, (x, y) in zip(matched(vtk, SQUARE_HOLE), SQUARE_HOLE, strict=True):
            px, py = vtk.points[k, :2]
            moved = OFF if abs(x - round(x)) > 1e-6 else 0.0  # the three x = ...0004 per edge
            assert abs(abs(px - x) - moved) <= 1e-9, (x, y, px)
            assert abs(py - y) <= 1e-9, (x, y, py)

    def test_the_record_carries_the_three_entries(self, frozen: Run) -> None:
        record = frozen.record
        assert record["frozen_holes"] == 1
        assert record["hole_vertices_snapped_max_m"] == pytest.approx(OFF, abs=1e-9)
        # 0.6 m, at local x = 420, 430, 470, 480, 520, 530, 570 and 580 on both edges;
        # measured with a stand-in on master: 0.6000000238.
        assert record["frozen_edges_max_error_m"] == pytest.approx(0.6, abs=1e-6)

    def test_the_mesh_file_carries_two_of_them_after_max_error_m(self, frozen: Run) -> None:
        vtk = frozen.vtk
        assert int(file_field(vtk, "frozen_holes")) == 1
        assert float(file_field(vtk, "frozen_edges_max_error_m")) == pytest.approx(0.6, abs=1e-6)
        assert "hole_vertices_snapped_max_m" not in vtk.field_data
        names = list(vtk.field_data)
        at = names.index("max_error_m")
        assert names[at + 1 : at + 3] == ["frozen_holes", "frozen_edges_max_error_m"]

    def test_the_written_vocabulary_is_as_without_the_flag(self, frozen: Run, plain: Run) -> None:
        vtk, before = frozen.vtk, plain.vtk
        for name in ("feature_names", "feature_bits", "feature_vocabulary"):
            assert np.array_equal(
                np.asarray(vtk.field_data[name].values), np.asarray(before.field_data[name].values)
            ), name
        masks = np.asarray(vtk.cell_array("feature_mask").values, dtype=np.int64)
        assert not (masks & BIT).any(), "the internal frozen-hole bit was written"

    def test_a_line_through_hole_vertices_only_is_accepted(self, inputs: Inputs) -> None:
        args = mesh_args(
            inputs.dem,
            inputs.square,
            *("--tolerance", str(TOLERANCE), "--frozen-holes"),
            *("--lines", f"wall={inputs.through_vertices}", "--edges", "wall"),
        )
        vtk = Run(inputs.root, "through", *args).vtk
        assert len(on_ring(vtk, SQUARE_HOLE)) == 10
        assert_ring_kept(vtk, SQUARE_HOLE)

    def test_the_budget_path_keeps_the_same_ten_vertices(self, frozen_budget: Run) -> None:
        vtk = frozen_budget.vtk
        assert len(on_ring(vtk, SQUARE_HOLE)) == 10
        assert_ring_kept(vtk, SQUARE_HOLE)


class TestF1SlantedHole:
    """F1 on a hole with slanted edges: the diamond beside the square."""

    def test_both_rings_are_kept(self, frozen_two: Run) -> None:
        vtk = frozen_two.vtk
        assert len(on_ring(vtk, SQUARE_HOLE)) == 10
        assert len(on_ring(vtk, DIAMOND)) == 4
        assert_ring_kept(vtk, SQUARE_HOLE)
        assert_ring_kept(vtk, DIAMOND)

    def test_the_diamond_corners_on_the_lattice_did_not_move(self, frozen_two: Run) -> None:
        vtk = frozen_two.vtk
        for k, (x, y) in zip(matched(vtk, DIAMOND), DIAMOND, strict=True):
            assert vtk.points[k, 0] == x and vtk.points[k, 1] == y

    def test_the_record_counts_two_holes_and_the_diamonds_larger_error(
        self, frozen_two: Run
    ) -> None:
        record = frozen_two.record
        assert record["frozen_holes"] == 2
        assert record["hole_vertices_snapped_max_m"] == pytest.approx(OFF, abs=1e-9)
        # 0.9 m at a diamond edge's middle (module docstring), above the square's 0.6 m.
        assert record["frozen_edges_max_error_m"] == pytest.approx(0.9, abs=1e-6)
        assert float(file_field(frozen_two.vtk, "frozen_edges_max_error_m")) == pytest.approx(
            0.9, abs=1e-6
        )


class TestF1ExactErrorFigure:
    """F1: the recorded error is exact on the DEM, not read between its nodes."""

    def test_a_hole_half_a_cell_off_the_columns_records_9_05(self, inputs: Inputs) -> None:
        # 9.05 m at the node x = 300 on the bottom and top edges (module docstring);
        # measured on 274a2efe, sampling at cell steps: 9.0000005.
        args = mesh_args(
            inputs.dem, inputs.off_grid, "--tolerance", str(TOLERANCE), "--frozen-holes"
        )
        run = Run(inputs.root, "off-grid", *args)
        assert_ring_kept(run.vtk, OFF_GRID_HOLE)
        assert run.record["frozen_edges_max_error_m"] == pytest.approx(9.05, abs=1e-5)


# ---------------------------------------------------------------- F2


def assert_within_off_the_rings(
    run: Run, holes: tuple[Ring, ...], counts: tuple[int, int, int], bound: float
) -> np.ndarray:
    """Every DEM node off the hole rings and outside the holes is in a written
    triangle within `bound` (plus 1e-9) of its plane. `counts` is (on a ring,
    inside a hole, off), by hand. Returns the errors of the nodes on the rings."""
    vtk = run.vtk
    xy, z = dem_nodes()
    on, inside = node_classes(xy, holes)
    off = ~on & ~inside
    assert (int(on.sum()), int(inside.sum()), int(off.sum())) == counts
    error = located_errors(vtk, xy[off], z[off])
    assert np.isfinite(error).all(), f"{int((~np.isfinite(error)).sum())} nodes in no triangle"
    assert error.max() <= bound + 1e-9, error.max()
    return located_errors(vtk, xy[on], z[on])


class TestF2ToleranceOffTheFrozenEdges:
    """F2: within the tolerance off the frozen edges; the largest error on them recorded."""

    def test_the_9760_nodes_off_the_ring_are_within_the_tolerance(self, frozen: Run) -> None:
        # measured with a stand-in on master: largest 0.0999999829
        on_errors = assert_within_off_the_rings(frozen, (SQUARE_HOLE,), (80, 361, 9760), TOLERANCE)
        recorded = frozen.record["frozen_edges_max_error_m"]
        assert isinstance(recorded, float)
        assert on_errors.max() == pytest.approx(recorded, abs=1e-6)

    def test_the_9659_nodes_off_both_rings_are_within_the_tolerance(self, frozen_two: Run) -> None:
        on_errors = assert_within_off_the_rings(
            frozen_two, (SQUARE_HOLE, DIAMOND), (88, 454, 9659), TOLERANCE
        )
        recorded = frozen_two.record["frozen_edges_max_error_m"]
        assert isinstance(recorded, float)
        assert on_errors.max() == pytest.approx(recorded, abs=1e-6)

    def test_the_budget_path_is_within_its_recorded_max_error_off_the_ring(
        self, frozen_budget: Run
    ) -> None:
        bound = frozen_budget.record["max_error_m"]
        assert isinstance(bound, float)
        assert_within_off_the_rings(frozen_budget, (SQUARE_HOLE,), (80, 361, 9760), bound)

    @pytest.mark.parametrize("name", ["frozen", "frozen_two", "frozen_budget"])
    def test_the_mesh_is_constrained_delaunay(
        self, name: str, request: pytest.FixtureRequest
    ) -> None:
        run: Run = request.getfixturevalue(name)
        assert constrained_delaunay_violations(run.vtk) == []


# ---------------------------------------------------------------- F3


def assert_refused(run: Run, *says: str) -> None:
    """Exit 2 for the stated reason, saying each of `says` (whitespace removed
    on both sides, as Rich may break a word), writing no file."""
    assert run.code == USAGE, run.output
    assert "Traceback" not in run.output
    assert "No such option" not in run.output, "refused for the wrong reason"
    flat = "".join(run.output.split())
    for word in says:
        assert "".join(word.split()) in flat, f"{word!r} not in {run.output!r}"
    assert not run.out.exists()
    assert not run.record_path.exists()
    assert not run.out.with_name(run.out.name + ".pieces").exists()


def file_point(x: float, y: float) -> str:
    fx, fy = local(x, y)
    return f"({fx:.3f}, {fy:.3f})"


class TestF3Refusals:
    """F3: the crossing line, the option refusals and the hole vertex with no DEM."""

    def test_a_line_meeting_hole_edges_between_vertices_is_refused(self, inputs: Inputs) -> None:
        wall = ("--lines", f"wall={inputs.between_vertices}", "--edges", "wall")
        base = mesh_args(inputs.dem, inputs.square, "--tolerance", str(TOLERANCE), *wall)
        Run(inputs.root, "between-plain", *base).ok()
        run = Run(inputs.root, "between", *base, "--frozen-holes")
        assert_refused(
            run,
            "a line meets the edge of hole 1 of --domain between its vertices at",
            "give the line a vertex there or mesh without it",
        )
        # The walk order is not a contract; measured today, (512.3, 600) comes first.
        flat = "".join(run.output.split())
        points = [file_point(512.3, 600.0), file_point(512.3, 400.0)]
        assert any("".join(p.split()) in flat for p in points), run.output

    @pytest.mark.parametrize(
        ("case", "says"),
        [
            ("no-domain", ("--frozen-holes", "needs --domain")),
            (
                "simplify",
                ("--frozen-holes keeps the holes as given; --simplify-domain would move them",),
            ),
            ("pieces", ("--frozen-holes does not run on pieces yet",)),
            ("piece-budget", ("--frozen-holes does not run on pieces yet",)),
        ],
    )
    def test_an_option_it_cannot_run_with_is_refused(
        self, inputs: Inputs, case: str, says: tuple[str, ...]
    ) -> None:
        """Section 3.1's refusals; each run is accepted without the flag (a cut
        run refuses --record on its own, so none is asked)."""
        size = ("--tolerance", str(TOLERANCE))
        args = {
            "no-domain": ("--dem", str(inputs.dem), *size),
            "simplify": mesh_args(
                inputs.dem, inputs.square, *size, "--simplify-domain", "--domain-max-shift", "1"
            ),
            # 4 pieces asked, and a budget the window exceeds (6 pieces measured)
            "pieces": mesh_args(inputs.dem, inputs.square, *size, "--pieces", "4"),
            "piece-budget": mesh_args(inputs.dem, inputs.square, *size, "--piece-budget", "1M"),
        }[case]
        Run(inputs.root, f"{case}-plain", *args, record=False).ok()
        assert_refused(Run(inputs.root, case, *args, "--frozen-holes", record=False), *says)

    def test_catchments_is_refused_before_any_land_cover_is_asked(self, inputs: Inputs) -> None:
        """Section 3.1: --catchments replaces the domain, so its holes would not be
        --domain's. One unit, the outer square, no land cover: today's land-cover
        message would be the wrong reason. No run without the flag: it needs land
        cover and a .gpkg output."""
        args = mesh_args(
            inputs.dem,
            inputs.square,
            *("--tolerance", str(TOLERANCE), "--frozen-holes"),
            *("--catchments", str(inputs.units), "--catchments-max-shift", "0"),
        )
        run = Run(inputs.root, "catchments", *args)
        assert_refused(
            run,
            "--catchments",
            "--frozen-holes keeps the holes of --domain; "
            "--catchments replaces the domain with the union of its units",
        )
        assert "land-cover" not in run.output, "refused for the land-cover reason"

    def test_hole_vertices_merged_by_the_snap_are_refused(self, inputs: Inputs) -> None:
        """Section 3.3's second message: 11 vertices, 10 after the 1 mm snap."""
        base = mesh_args(inputs.dem, inputs.merging, "--tolerance", str(TOLERANCE))
        Run(inputs.root, "merging-plain", *base).ok()  # measured: 5,438 triangles
        assert_refused(
            Run(inputs.root, "merging", *base, "--frozen-holes"),
            "hole 1 of --domain has 11 vertices but 10 after snapping to the 0.001 m grid; "
            "vertices closer than that merged",
        )

    # The first end named: measured with a stand-in on master, (400, 400) on the
    # tolerance run and (400, 600) on the budget run's kept trial; the order is
    # not a contract, so either is accepted.
    @pytest.mark.parametrize(
        "size",
        [("--tolerance", str(TOLERANCE)), ("--triangles", "3000")],
        ids=["tolerance", "budget"],
    )
    def test_a_hole_vertex_with_no_dem_is_refused(
        self, inputs: Inputs, size: tuple[str, str]
    ) -> None:
        name = size[0].strip("-")
        base = mesh_args(inputs.nodata_dem, inputs.square, *size)
        Run(inputs.root, f"nodata-plain-{name}", *base).ok()
        run = Run(inputs.root, f"nodata-{name}", *base, "--frozen-holes")
        assert_refused(
            run,
            "3 of the 10 hole edges of --domain have an end with no DEM under it "
            "(nodata or off the DEM), so the mesh could not be joined along the hole; "
            "give holes that lie on DEM data; first at",
        )
        flat = "".join(run.output.split())
        points = [file_point(400.0, 400.0), file_point(400.0, 600.0)]
        assert any("".join(p.split()) in flat for p in points), run.output


# ---------------------------------------------------------------- F4


class TestF4OffIsIdentity:
    """F4: without the flag nothing changes."""

    def test_without_the_flag_the_hole_edges_are_split_as_today(self, plain: Run) -> None:
        # measured on master: 28 vertices on the hole ring
        assert len(on_ring(plain.vtk, SQUARE_HOLE)) > 10
        record = plain.record
        assert not set(ENTRIES) & set(record)
        assert not set(ENTRIES) & set(plain.vtk.field_data)

    def test_a_domain_without_a_hole_gives_the_same_files(self, inputs: Inputs) -> None:
        base = mesh_args(inputs.dem, inputs.no_hole, "--tolerance", str(TOLERANCE))
        off = Run(inputs.root, "no-hole-off", *base).ok()
        on = Run(inputs.root, "no-hole-on", *base, "--frozen-holes").ok()
        assert on.out.read_bytes() == off.out.read_bytes()
        a, b = dict(off.record), dict(on.record)
        a.pop("command")
        b.pop("command")
        assert a == b
