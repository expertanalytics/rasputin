"""`rasputin mesh --dem PATH --domain PATH --tolerance METRES`: increment 16.

`docs/increments/16-domain-polygon.md`, R1 to R4, with the user's rulings R0
(point heights, bilinear between nodes), no snapping of input vertices, U1 (a)
(the domain's CRS must match the DEM's), U3 (a), U4 (a) (refuse a polygon
outside the node rectangle), U5 (a) (no boundary bit) and U6 (a) (the noder's
1 mm snap is the engine's input precision).

Wording pinned from the design: the ``domain`` field is
``<file name>, 1 ring <h> holes, <n> vertices`` (``hole`` or ``holes``
accepted, since the design gives only the zero-hole example); the sentence says
``start domain boundary, boundary z bilinear`` in place of ``start stride``.
The stderr wording for the start triangle count and the off-node count is not
ruled; this suite pins ``<n> start triangles`` and ``<k> start vertices
off-node``.

The real-tile acceptance (T-real) runs under ``needs_codecs`` only. It builds
the quarter circle from the design's numbers, reports, and asserts only the
tolerance.
"""

from __future__ import annotations

import io
import json
import re
import time
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from shapely.geometry import Point, Polygon

from geotiff_fixtures import KARTVERKET, TIE_X, TIE_Y, micro_tiff, needs_codecs
from plyread import read_ply
from test_cli_mesh_dem import SENTINEL, USAGE, invoke, write_tiff
from test_cli_mesh_refine import NUMBER, field, min_angles_degrees, sentence
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import DemTile
from vtkread import VtkFile, read_vtk

UTM33 = "urn:ogc:def:crs:EPSG::25833"
SNAP = 1e-3  # DEFAULT_SNAP_SPACING, U6 (a)
Ring = list[tuple[float, float]]


def geojson(path: Path, outer: Ring, holes: tuple[Ring, ...] = (), crs: str | None = UTM33) -> Path:
    doc: dict[str, Any] = {
        "type": "Polygon",
        "coordinates": [[*r, r[0]] for r in (outer, *holes)],
    }
    if crs is not None:
        doc["crs"] = {"type": "name", "properties": {"name": crs}}
    path.write_text(json.dumps(doc))
    return path


def mesh(
    tmp_path: Path, tif: Path, domain: Path, *extra: str, out: str = "x.vtk"
) -> tuple[int, str, Path]:
    target = tmp_path / out
    code, output = invoke("--dem", str(tif), "--domain", str(domain), "--out", str(target), *extra)
    return code, output, target


def meshed(tmp_path: Path, tif: Path, domain: Path, *extra: str) -> tuple[VtkFile, str]:
    code, output, target = mesh(tmp_path, tif, domain, "--tolerance", "1", *extra)
    assert code == 0, output
    return read_vtk(target.read_bytes()), output


def bilinear(tile: DemTile, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """R0, written here from the design's text: nodes at their value, bilinear between."""
    m = tile.meta
    col, row = (x - m.x_min) / m.delta_x, (m.y_max - y) / m.delta_y
    c0 = np.minimum(np.floor(col).astype(int), m.cols - 2)
    r0 = np.minimum(np.floor(row).astype(int), m.rows - 2)
    tx, ty = col - c0, row - r0
    a = tile.array.astype(np.float64)
    return (
        a[r0, c0] * (1 - tx) * (1 - ty)
        + a[r0, c0 + 1] * tx * (1 - ty)
        + a[r0 + 1, c0] * (1 - tx) * ty
        + a[r0 + 1, c0 + 1] * tx * ty
    )


def nearest(points: np.ndarray, x: float, y: float) -> int:
    return int(np.argmin(np.hypot(points[:, 0] - x, points[:, 1] - y)))


# ---------------------------------------------------------------- fixtures

# micro_tiff's grid: nodes at x = TIE_X + 10 col, y = TIE_Y - 5 row, EPSG:25833,
# PixelIsPoint. 17 rows x 21 cols: x 500 000 .. 500 200, y 6 599 920 .. 6 600 000.
ROWS, COLS = 17, 21

# A square and a hole, every vertex off-node (x not a multiple of 10, y not of 5).
SQUARE: Ring = [
    (TIE_X + 12.3, TIE_Y - 73.3),
    (TIE_X + 187.7, TIE_Y - 72.9),
    (TIE_X + 186.1, TIE_Y - 6.7),
    (TIE_X + 13.9, TIE_Y - 7.1),
]
HOLE: Ring = [
    (TIE_X + 71.1, TIE_Y - 52.7),
    (TIE_X + 72.3, TIE_Y - 28.3),
    (TIE_X + 121.9, TIE_Y - 28.9),
    (TIE_X + 120.7, TIE_Y - 51.1),
]


@pytest.fixture
def bumpy(tmp_path: Path) -> Path:
    array = np.random.default_rng(16).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))


@pytest.fixture
def square(tmp_path: Path) -> Path:
    return geojson(tmp_path / "square.geojson", SQUARE, (HOLE,))


def island(n: int) -> np.ndarray:
    """T12's island (increment 14b): flat sea at 0, land from 3 m on a circular coast."""
    mid, radius = (n - 1) / 2.0, 0.3 * n
    r, c = np.indices((n, n))
    d = np.hypot(r - mid, c - mid)
    return np.where(d <= radius, 3.0 + 0.5 * (radius - d), 0.0).astype(np.float32)


def cone(n: int) -> np.ndarray:
    mid = (n - 1) / 2.0
    r, c = np.indices((n, n))
    return (100.0 - np.hypot(r - mid, c - mid)).astype(np.float32)


# ---------------------------------------------------------------- tests


class TestDomainOutput:
    """R3 and R4, Z1 and Z3 through the CLI."""

    def test_the_field_the_sentence_and_the_report(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        vtk, output = meshed(tmp_path, bumpy, square)
        (domain,) = vtk.field_data["domain"].values
        assert re.fullmatch(r"square\.geojson, 1 ring 1 holes?, 8 vertices", str(domain)), domain
        text = sentence(vtk)
        assert "start domain boundary, boundary z bilinear" in text
        assert "start stride" not in text
        assert "refined from DEM nodes" in text
        assert field(text, rf"achieved max error {NUMBER} m") <= 1.0
        assert re.search(r"\b\d+ start triangles\b", output), output
        assert re.search(r"\b8 start vertices off-node\b", output), output

    def test_the_ply_carries_the_domain_comment(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        code, output, target = mesh(tmp_path, bumpy, square, "--tolerance", "1", out="x.ply")
        assert code == 0, output
        header, _ = read_ply(target.read_bytes())
        assert any(c.startswith("domain square.geojson, 1 ring 1 hole") for c in header.comments), (
            header.comments
        )

    def test_input_vertices_come_out_where_the_file_put_them_within_the_snap(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        """U6 (a): the noder moves an input vertex by at most half the snap per axis."""
        vtk, _ = meshed(tmp_path, bumpy, square)
        for x, y in SQUARE + HOLE:
            i = nearest(vtk.points, x, y)
            assert abs(vtk.points[i, 0] - x) <= SNAP / 2 + 1e-9
            assert abs(vtk.points[i, 1] - y) <= SNAP / 2 + 1e-9

    def test_an_input_vertex_on_the_snap_grid_comes_out_bit_for_bit(self, tmp_path: Path) -> None:
        """Z3: positions are the input's, never recomputed from the fractional frame.

        With the grid's x_min at 0, x = 12.708 does not survive the frame round
        trip, x_min + ((x - x_min) / dx) * dx; near TIE_X = 500 000 every
        coordinate does, which is why the other fixtures cannot tell. 12.708 is
        on the 1 mm snap grid, so U6's snap leaves it where it is and the output
        must hold it exactly.
        """
        x_min = 0.0
        odd = (12_708 * SNAP, TIE_Y - 73.3)
        assert x_min + ((odd[0] - x_min) / 10.0) * 10.0 != odd[0]
        assert round(odd[0] / SNAP) * SNAP == odd[0]
        array = np.random.default_rng(16).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
        tif = write_tiff(
            tmp_path / "origin.tif", micro_tiff(array, tiepoint=(0.0, 0.0, 0.0, x_min, TIE_Y, 0.0))
        )
        ring = [odd, *((x - TIE_X, y) for x, y in SQUARE[1:])]
        vtk, _ = meshed(tmp_path, tif, geojson(tmp_path / "origin.geojson", ring))
        for x, y in ring:
            snapped = (round(x / SNAP) * SNAP, round(y / SNAP) * SNAP)
            i = nearest(vtk.points, *snapped)
            assert (vtk.points[i, 0], vtk.points[i, 1]) == snapped

    def test_boundary_z_is_bilinear_and_inserted_vertices_are_nodes(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        vtk, _ = meshed(tmp_path, bumpy, square)
        tile = decode_dem(io.BytesIO(bumpy.read_bytes()))
        m = tile.meta
        boundary = {nearest(vtk.points, x, y) for x, y in SQUARE + HOLE}
        idx = np.array(sorted(boundary))
        p = vtk.points[idx]
        assert p[:, 2] == pytest.approx(bilinear(tile, p[:, 0], p[:, 1]), rel=1e-9, abs=1e-9)
        rest = np.delete(vtk.points, idx, axis=0)
        assert len(rest) > 0, "rough terrain at 1 m must insert nodes"
        cols = (rest[:, 0] - m.x_min) / m.delta_x
        rows = (m.y_max - rest[:, 1]) / m.delta_y
        assert (cols == np.round(cols)).all() and (rows == np.round(rows)).all()
        expected = tile.array[rows.astype(int), cols.astype(int)].astype(np.float64)
        assert rest[:, 2] == pytest.approx(expected, rel=1e-9, abs=1e-9)


class TestSyntheticEndToEnd:
    """The T12 cone and island with a square domain holding one hole, all off-node."""

    N = 65  # x 500 000 .. 500 640, y 6 599 680 .. 6 600 000
    OUTER: tuple[tuple[float, float], ...] = (
        (TIE_X + 13.7, TIE_Y - 316.1),
        (TIE_X + 627.3, TIE_Y - 316.9),
        (TIE_X + 626.9, TIE_Y - 3.9),
        (TIE_X + 12.1, TIE_Y - 4.3),
    )
    INNER: tuple[tuple[float, float], ...] = (
        (TIE_X + 250.3, TIE_Y - 199.1),
        (TIE_X + 251.7, TIE_Y - 119.7),
        (TIE_X + 390.7, TIE_Y - 120.3),
        (TIE_X + 389.9, TIE_Y - 198.3),
    )

    @pytest.mark.parametrize("terrain", ["cone", "island"])
    @pytest.mark.parametrize("tolerance", ["1", "0.25"])
    def test_nothing_outside_the_square_or_inside_the_hole(
        self, tmp_path: Path, terrain: str, tolerance: str
    ) -> None:
        array = cone(self.N) if terrain == "cone" else island(self.N)
        tif = write_tiff(tmp_path / f"{terrain}.tif", micro_tiff(array))
        domain = geojson(tmp_path / "d.geojson", list(self.OUTER), (list(self.INNER),))
        code, output, target = mesh(tmp_path, tif, domain, "--tolerance", tolerance)
        assert code == 0, output
        vtk = read_vtk(target.read_bytes())
        assert field(sentence(vtk), rf"achieved max error {NUMBER} m") <= float(tolerance)

        shape = Polygon(self.OUTER, [self.INNER])
        grown, hole = shape.buffer(SNAP), Polygon(self.INNER).buffer(-SNAP)
        xy = vtk.points[:, :2]
        assert all(grown.covers(Point(p)) for p in xy)
        assert not any(hole.contains(Point(p)) for p in xy)
        tris = np.asarray(vtk.polygons, dtype=np.int64)
        centroids = xy[tris].mean(axis=1)
        assert all(shape.contains(Point(c)) for c in centroids)
        u, v = xy[tris[:, 1]] - xy[tris[:, 0]], xy[tris[:, 2]] - xy[tris[:, 0]]
        area = 0.5 * np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
        assert area.sum() == pytest.approx(shape.area, rel=1e-6)


class TestNoDataAndDegenerateInput:
    def test_a_boundary_vertex_in_a_nodata_cell_is_dropped_and_counted(
        self, tmp_path: Path, square: Path
    ) -> None:
        """Z2: SQUARE[0] is at col 1.23, row 14.66; its cell has corner (15, 2)."""
        array = np.random.default_rng(2).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
        array[15, 2] = float(SENTINEL)
        tif = write_tiff(tmp_path / "void.tif", micro_tiff(array, nodata=SENTINEL))
        vtk, _ = meshed(tmp_path, tif, square)
        assert (vtk.points[:, 2] != float(SENTINEL)).all()
        assert field(sentence(vtk), rf"{NUMBER} vertices without data dropped") >= 1
        x, y = SQUARE[0]
        assert np.hypot(vtk.points[:, 0] - x, vtk.points[:, 1] - y).min() > SNAP

    def test_nearly_collinear_and_near_node_vertices(self, tmp_path: Path, bumpy: Path) -> None:
        """A vertex 1 cm off a 170 m edge, and one 1 um from a node (the snap lands it on it)."""
        outer: Ring = [
            (TIE_X + 12.3, TIE_Y - 73.3),
            (TIE_X + 100.0 + 1e-6, TIE_Y - 73.1 - 0.01),
            (TIE_X + 187.7, TIE_Y - 72.9),
            (TIE_X + 186.1, TIE_Y - 6.7),
            (TIE_X + 100.0, TIE_Y - 5.0 + 1e-6),
            (TIE_X + 13.9, TIE_Y - 7.1),
        ]
        domain = geojson(tmp_path / "thin.geojson", outer)
        vtk, _ = meshed(tmp_path, bumpy, domain)
        assert field(sentence(vtk), rf"achieved max error {NUMBER} m") <= 1.0
        for x, y in outer:
            i = nearest(vtk.points, x, y)
            assert abs(vtk.points[i, 0] - x) <= SNAP / 2 + 1e-9
            assert abs(vtk.points[i, 1] - y) <= SNAP / 2 + 1e-9


class TestRefusals:
    """R1, every refusal exit 2 with no file written."""

    @staticmethod
    def refused(
        tmp_path: Path, bumpy: Path, domain: Path, *extra: str, says: tuple[str, ...]
    ) -> None:
        code, output, target = mesh(tmp_path, bumpy, domain, *extra)
        assert code == USAGE, output
        assert "No such option" not in output, "refused for the wrong reason"
        assert "Traceback" not in output
        for word in says:
            assert word in output, output
        assert not target.exists()

    def test_a_multipolygon(self, tmp_path: Path, bumpy: Path) -> None:
        a = [
            [TIE_X + 5, TIE_Y - 50],
            [TIE_X + 50, TIE_Y - 50],
            [TIE_X + 50, TIE_Y - 10],
            [TIE_X + 5, TIE_Y - 50],
        ]
        b = [[x + 100, y] for x, y in a]
        doc = {
            "type": "MultiPolygon",
            "coordinates": [[a], [b]],
            "crs": {"type": "name", "properties": {"name": UTM33}},
        }
        path = tmp_path / "multi.geojson"
        path.write_text(json.dumps(doc))
        self.refused(tmp_path, bumpy, path, "--tolerance", "1", says=("--domain", "MultiPolygon"))

    def test_invalid(self, tmp_path: Path, bumpy: Path) -> None:
        bowtie: Ring = [
            (TIE_X + 10, TIE_Y - 60),
            (TIE_X + 90, TIE_Y - 10),
            (TIE_X + 90, TIE_Y - 60),
            (TIE_X + 10, TIE_Y - 10),
        ]
        path = geojson(tmp_path / "bowtie.geojson", bowtie)
        self.refused(
            tmp_path, bumpy, path, "--tolerance", "1", says=("--domain", "Self-intersection")
        )

    def test_empty(self, tmp_path: Path, bumpy: Path) -> None:
        path = tmp_path / "empty.wkt"
        path.write_text("POLYGON EMPTY")
        self.refused(
            tmp_path,
            bumpy,
            path,
            "--tolerance",
            "1",
            "--domain-crs",
            "EPSG:25833",
            says=("--domain",),
        )

    def test_geojson_without_crs(self, tmp_path: Path, bumpy: Path) -> None:
        path = geojson(tmp_path / "wgs.geojson", SQUARE, crs=None)
        self.refused(tmp_path, bumpy, path, "--tolerance", "1", says=("4326", "25833"))

    def test_a_mismatched_epsg(self, tmp_path: Path, bumpy: Path) -> None:
        path = geojson(tmp_path / "utm32.geojson", SQUARE, crs="EPSG:25832")
        self.refused(tmp_path, bumpy, path, "--tolerance", "1", says=("25832", "25833"))

    def test_wkt_without_domain_crs(self, tmp_path: Path, bumpy: Path) -> None:
        path = tmp_path / "square.wkt"
        path.write_text(Polygon(SQUARE).wkt)
        self.refused(tmp_path, bumpy, path, "--tolerance", "1", says=("--domain-crs",))

    def test_wkt_with_domain_crs_is_accepted(self, tmp_path: Path, bumpy: Path) -> None:
        path = tmp_path / "square.wkt"
        path.write_text(Polygon(SQUARE).wkt)
        code, output, target = mesh(
            tmp_path, bumpy, path, "--tolerance", "1", "--domain-crs", "EPSG:25833"
        )
        assert code == 0, output
        assert target.exists()

    def test_geojson_with_an_agreeing_domain_crs_is_accepted(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        code, output, target = mesh(
            tmp_path, bumpy, square, "--tolerance", "1", "--domain-crs", "EPSG:25833"
        )
        assert code == 0, output
        assert target.exists()

    def test_geojson_with_a_disagreeing_domain_crs(
        self, tmp_path: Path, bumpy: Path, square: Path
    ) -> None:
        self.refused(
            tmp_path,
            bumpy,
            square,
            "--tolerance",
            "1",
            "--domain-crs",
            "EPSG:25832",
            says=("--domain-crs", "25832", "25833"),
        )

    def test_a_vertex_outside_the_node_rectangle(self, tmp_path: Path, bumpy: Path) -> None:
        outer = [*SQUARE[:1], (TIE_X + 200.5, TIE_Y - 72.9), *SQUARE[2:]]
        path = geojson(tmp_path / "wide.geojson", outer)
        self.refused(tmp_path, bumpy, path, "--tolerance", "1", says=("--domain", "outside"))

    def test_domain_with_stride(self, tmp_path: Path, bumpy: Path, square: Path) -> None:
        self.refused(
            tmp_path, bumpy, square, "--tolerance", "1", "--stride", "4", says=("--stride",)
        )

    def test_domain_without_tolerance(self, tmp_path: Path, bumpy: Path, square: Path) -> None:
        self.refused(tmp_path, bumpy, square, says=("--tolerance",))

    def test_domain_without_dem(self, tmp_path: Path, square: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "catchment", "--domain", str(square), "--tolerance", "1", "--out", str(target)
        )
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--domain" in output
        assert not target.exists()

    def test_a_missing_domain_file(self, tmp_path: Path, bumpy: Path) -> None:
        self.refused(
            tmp_path, bumpy, tmp_path / "absent.geojson", "--tolerance", "1", says=("--domain",)
        )


# ---------------------------------------------------------------- T-real

CENTRE = (850_250.0, 7_899_750.0)
RADIUS = 30_000.0
SPACING = 200.0


def quarter_circle() -> Ring:
    """The design's acceptance ring, counter-clockwise: up the east border from
    the centre, along the arc, back along the south border. 536 vertices."""
    cx, cy = CENTRE
    edge = round(RADIUS / SPACING)  # 150 segments per border edge
    arc = int(np.ceil(np.pi / 2 * RADIUS / SPACING))  # 236 segments
    ring: Ring = [(cx, cy + RADIUS * i / edge) for i in range(edge)]
    ring += [
        (cx + RADIUS * np.cos(t), cy + RADIUS * np.sin(t))
        for t in np.linspace(np.pi / 2, np.pi, arc + 1)[:-1]
    ]
    ring += [(cx - RADIUS + RADIUS * i / edge, cy) for i in range(edge)]
    ring[edge] = (cx, cy + RADIUS)  # the arc's ends exactly on the border lines
    ring[edge + arc] = (cx - RADIUS, cy)
    return [(float(x), float(y)) for x, y in ring]


def test_the_quarter_circle_ring_is_the_designs() -> None:
    ring = quarter_circle()
    assert len(ring) == 536
    assert Polygon(ring).is_valid and Polygon(ring).exterior.is_ccw
    steps = np.hypot(*np.diff(np.array([*ring, ring[0]]), axis=0).T)
    assert steps.max() <= SPACING + 1e-6


class TestRealTile:
    """T-real: logged, not thresholded, except the tolerance."""

    @needs_codecs
    def test_meshes_the_quarter_circle(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        domain = geojson(tmp_path / "quarter.geojson", quarter_circle())
        for tolerance in ("10", "1"):
            began = time.perf_counter()
            code, output, target = mesh(tmp_path, KARTVERKET, domain, "--tolerance", tolerance)
            seconds = time.perf_counter() - began
            assert code == 0, output
            vtk = read_vtk(target.read_bytes())
            text = sentence(vtk)
            achieved = field(text, rf"achieved max error {NUMBER} m")
            assert achieved <= float(tolerance)
            tris = np.asarray(vtk.polygons, dtype=np.int64)
            degree = np.bincount(tris.ravel(), minlength=len(vtk.points))
            angles = min_angles_degrees(vtk)
            p99 = np.percentile(degree, 99)
            with capsys.disabled():
                print(
                    f"\nT-real tolerance {tolerance} m: {len(tris)} triangles, "
                    f"{len(vtk.points)} vertices, achieved {achieved} m, {seconds:.1f} s"
                    f"\nT-real degree: median {np.median(degree):.0f}, p99 {p99:.0f}, "
                    f"max {degree.max()}, >= 12: {(degree >= 12).sum()}, "
                    f">= 20: {(degree >= 20).sum()}"
                    f"\nT-real min angle: median {np.median(angles):.2f} deg, "
                    f"{100 * np.mean(angles < 1.0):.2f} % under 1 deg, "
                    f"{100 * np.mean(angles < 10.0):.2f} % under 10 deg, "
                    f"worst {angles.min():.4f} deg"
                    f"\nT-real report: {output.strip()}; {text}"
                )
