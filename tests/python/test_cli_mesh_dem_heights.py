"""`rasputin mesh --dem-heights`: each triangle's share of area per height
class about its own plane, its share over NoData, and the run's class edges,
in the mesh file.

`docs/increments/137-dem-stats.md`, tests 1, 2, 10, 11, 13, 14, 16, 17, 18,
20 and 21 of section 5, invariants 1, 2, 4, 5, 7, 8 and 9, and the stable
names of 3.6. The scenes:

- STEPS (tests 1, 2, 14, 17, 18, 20): `micro_tiff` 11 x 11 nodes 10 m apart,
  2 m on columns 1 and 2 and 0 elsewhere, `--stride 10`: two triangles on
  plane 0 (measured at 98fe80c3: two at `--tolerance` 10 and 2), whichever
  diagonal the mesh takes (read from the file);
- the cut run (tests 10, 13): `piece_spies.square_scene`, `--pieces 4`;
- Velhas (test 11): `tests/fixtures/velhas`, `--out-crs EPSG:31983
  --tolerance 5`, measured on the source's own pixels (EPSG:4674,
  area-registered), never the resampled grid;
- the hills (test 21): `pieces_fixtures.hills(n=201)`, `--stride 200`.

The oracle is `height_fixtures.oracle_table` from the file's own vertices.
Bounds: 1e-7 on a hand figure (float32 rounding), 1e-5 on the cut run and
1e-4 on the reprojected path (section 5); the Velhas figures 0.01 m and
0.005 about the probe's 1.164 m and 0.228 (the resampled grid's 1.121 m and
0.262 lie outside both).

A DEM given by path has no source entry, so w = 1 m (3.1).
"""

from __future__ import annotations

import json
import sqlite3
from contextlib import closing
from itertools import pairwise
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from cli_helpers import help_text, invoke, squashed
from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff
from height_fixtures import (
    CLASS_FROM,
    CLASS_M,
    DEFAULT_FROM,
    NAMES,
    VELHAS,
    WITHOUT_DATA,
    Table,
    edges,
    gpkg_table,
    half_of,
    oracle_table,
    plan_areas,
    ply_table,
    row_sums,
    spread_and_zero,
    velhas_source,
    vtk_table,
)
from piece_spies import square_scene
from pieces_fixtures import Cut, H, cut, hills, write_dem
from plyread import read_ply, vertex_array
from recordread import file_field, ply_fields, stats_row
from tin_engine.io.vtk_read import read_mesh
from vtkread import VtkFile, polygons_as_array, read_vtk

FLAG = "--dem-heights"
HAND = 1e-7
CUT_RUN = 1e-5
REPROJECTED = 1e-4
VELHAS_ARGS = (
    *("--dem", str(VELHAS / "anadem_velhas.tif")),
    *("--domain", str(VELHAS / "catchment.geojson")),
    *("--out-crs", "EPSG:31983", "--tolerance", "5"),
)


def steps(nodata: bool = False) -> np.ndarray:
    a = np.zeros((11, 11), np.float32)
    a[:, 1:3] = 2.0
    if nodata:
        a[5, 3] = -9999.0  # world (30, -50) from the tie point
    return a


def write_steps(where: Path, nodata: bool = False) -> Path:
    where.mkdir(parents=True, exist_ok=True)
    path = where / "steps.tif"
    extra = {"nodata": "-9999"} if nodata else {}
    path.write_bytes(micro_tiff(steps(nodata), scale=(10.0, 10.0, 0.0), **extra).getvalue())
    return path


def run_ok(*args: str) -> str:
    code, output = invoke("mesh", *args)
    assert code == 0, output
    return output


def mesh_steps(
    where: Path, out: str, *extra: str, nodata: bool = False, tolerance: str = "10"
) -> Path:
    tif = write_steps(where, nodata)
    target = where / out
    run_ok(
        "--dem", str(tif), "--stride", "10", "--tolerance", tolerance, "--out", str(target), *extra
    )
    return target


def west_and_east(points: np.ndarray, triangles: np.ndarray) -> tuple[int, int]:
    """The triangle holding both west corners (x at the tie point), then the other."""
    west = [t for t, tri in enumerate(triangles) if (points[tri, 0] == TIE_X).sum() == 2]
    assert len(west) == 1, triangles
    return west[0], 1 - west[0]


def read_any(target: Path) -> tuple[np.ndarray, np.ndarray, Table, np.ndarray]:
    """Points, triangles, the table and the edges of a `.vtk` or `.ply`."""
    blob = target.read_bytes()
    if target.suffix == ".vtk":
        vtk = read_vtk(blob)
        table, edges_m, _ = vtk_table(vtk)
        return np.asarray(vtk.points), polygons_as_array(vtk), table, edges_m
    _, data = read_ply(blob)
    table, edges_m = ply_table(blob)
    faces = np.array([list(f) for f in data["face"]["vertex_indices"]])
    return vertex_array(data), faces, table, edges_m


def assert_row(table: Table, t: int, first: int, want: list[float], nod: float = 0.0) -> None:
    """Row `t`: first and count exactly, each share within `HAND`, zeros exactly."""
    assert (int(table.first[t]), int(table.count[t])) == (first, len(want)), (
        table.first[t],
        table.count[t],
    )
    row = table.row(t)
    assert np.abs(row - np.array(want)).max() <= HAND, row
    assert_array_equal(row[np.array(want) == 0.0], 0.0)
    assert abs(table.nodata[t] - nod) <= HAND


# ------------------------------------------------- tests 1, 2 and 18


@pytest.mark.parametrize("out", ["steps.vtk", "steps.ply"])
class TestSteps:
    """Test 1 (n = 10), test 2 (node (5, 3) NoData) and test 18 (n = 2)."""

    def test_1_both_triangles_whichever_the_diagonal(self, tmp_path: Path, out: str) -> None:
        points, triangles, table, edges_m = read_any(mesh_steps(tmp_path, out, FLAG))
        assert len(triangles) == 2
        west, east = west_and_east(points, triangles)
        assert_row(table, west, 10, [33 / 50, 0.0, 17 / 50])
        assert_row(table, east, 10, [47 / 50, 0.0, 3 / 50])
        assert edges_m.tobytes() == edges(10, 1.0).tobytes()

    def test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share(
        self, tmp_path: Path, out: str
    ) -> None:
        points, triangles, table, _ = read_any(mesh_steps(tmp_path, out, FLAG, nodata=True))
        west, east = west_and_east(points, triangles)
        # The denominator keeps NoData's weight: 16/25, not 32/49.
        assert_row(table, west, 10, [16 / 25, 0.0, 17 / 50], nod=1 / 50)
        assert_row(table, east, 10, [47 / 50, 0.0, 3 / 50])

    def test_18_the_count_follows_the_tolerance(self, tmp_path: Path, out: str) -> None:
        target = mesh_steps(tmp_path, out, FLAG, tolerance="2")
        points, triangles, table, edges_m = read_any(target)
        assert len(triangles) == 2
        assert edges_m.tolist() == [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
        west, east = west_and_east(points, triangles)
        # The 2 m nodes at exactly the tolerance are inside: class +2 is the outermost.
        assert_row(table, west, 2, [33 / 50, 0.0, 17 / 50])
        assert_row(table, east, 2, [47 / 50, 0.0, 3 / 50])


class TestRecord:
    """3.3 and 3.6: the record entries in the file, `--stats` and `--record`."""

    def test_2_no_triangle_is_without_data_and_the_width_is_the_default(
        self, tmp_path: Path
    ) -> None:
        stats = tmp_path / "s.md"
        target = mesh_steps(tmp_path, "steps.vtk", FLAG, "--stats", str(stats), nodata=True)
        vtk = read_vtk(target.read_bytes())
        report = stats.read_text()
        assert int(file_field(vtk, WITHOUT_DATA)) == 0
        assert int(stats_row(report, WITHOUT_DATA)) == 0
        assert float(file_field(vtk, CLASS_M)) == 1.0
        assert float(stats_row(report, CLASS_M)) == 1.0
        assert file_field(vtk, CLASS_FROM) == DEFAULT_FROM
        assert stats_row(report, CLASS_FROM) == DEFAULT_FROM

    def test_the_ply_comments_carry_the_entries(self, tmp_path: Path) -> None:
        header, _ = read_ply(mesh_steps(tmp_path, "steps.ply", FLAG).read_bytes())
        fields = ply_fields(
            tuple(c for c in header.comments if not c.startswith("dem_height_edges_m"))
        )
        assert int(fields[WITHOUT_DATA]) == 0
        assert float(fields[CLASS_M]) == 1.0
        assert fields[CLASS_FROM] == DEFAULT_FROM

    def test_the_json_record_has_the_entries(self, tmp_path: Path) -> None:
        record = tmp_path / "r.json"
        mesh_steps(tmp_path, "steps.vtk", FLAG, "--record", str(record))
        text = json.dumps(json.loads(record.read_text()))
        for name in (WITHOUT_DATA, CLASS_M, CLASS_FROM):
            assert name in text, name


class TestOffWritesNothing:
    """Invariant 1's names: without the flag no array, comment, table, column
    or record entry appears. Runs on every platform (the digests of
    `test_dem_heights_off.py` are keyed by one)."""

    ENTRIES = (WITHOUT_DATA, CLASS_M, CLASS_FROM)

    def test_vtk(self, tmp_path: Path) -> None:
        stats = tmp_path / "s.md"
        vtk = read_vtk(mesh_steps(tmp_path, "steps.vtk", "--stats", str(stats)).read_bytes())
        arrays = set(vtk.field_data) | set(vtk.scalars)
        arrays |= {n for block in vtk.cell_fields.values() for n in block}
        assert not arrays & {*NAMES, *self.ENTRIES}, arrays
        assert not any(name in stats.read_text() for name in self.ENTRIES)

    def test_ply(self, tmp_path: Path) -> None:
        blob = mesh_steps(tmp_path, "steps.ply").read_bytes()
        header = blob[: blob.index(b"end_header")].decode()
        assert not any(name in header for name in (*NAMES, *self.ENTRIES)), header

    def test_gpkg(self, tmp_path: Path) -> None:
        target = mesh_steps(tmp_path, "steps.gpkg")
        with closing(sqlite3.connect(target)) as conn:
            tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master")}
            columns = {r[1] for r in conn.execute("PRAGMA table_info(triangles)")}
            dumped = "\n".join(conn.iterdump())
        assert not tables & {"dem_height_shares", "dem_height_classes"}, tables
        assert "dem_nodata_share" not in columns
        assert not any(name in dumped for name in self.ENTRIES)

    def test_read_mesh_gives_none(self, tmp_path: Path) -> None:
        assert read_mesh(mesh_steps(tmp_path, "steps.vtk").read_bytes()).dem_heights is None


# ------------------------------------------------- test 14: formats agree


class TestFormatsAgree:
    """Test 14: `.vtk`, `.ply` and `.gpkg` of test 1 carry the same numbers."""

    @pytest.fixture
    def three(self, tmp_path: Path) -> dict[str, tuple[Table, np.ndarray]]:
        vtk = read_vtk(mesh_steps(tmp_path / "v", "steps.vtk", FLAG).read_bytes())
        table, edges_m, _ = vtk_table(vtk)
        return {
            "vtk": (table, edges_m),
            "ply": ply_table(mesh_steps(tmp_path / "p", "steps.ply", FLAG).read_bytes()),
            "gpkg": gpkg_table(mesh_steps(tmp_path / "g", "steps.gpkg", FLAG)),
        }

    def test_the_same_table_and_edges(self, three: dict[str, tuple[Table, np.ndarray]]) -> None:
        table, edges_m = three["vtk"]
        assert len(table.first) == 2
        assert len(edges_m) == 22 and (edges_m[0], edges_m[-1]) == (-10.5, 10.5)
        for name in ("ply", "gpkg"):
            other, other_edges = three[name]
            assert_array_equal(other.offsets, table.offsets, err_msg=name)
            assert_array_equal(other.first, table.first, err_msg=name)
            assert_array_equal(other.shares, table.shares, err_msg=name)
            assert_array_equal(other.nodata, table.nodata, err_msg=name)
            assert other_edges.tobytes() == edges_m.tobytes(), name

    def test_the_gpkg_keeps_the_zeros_inside_a_span(self, tmp_path: Path) -> None:
        target = mesh_steps(tmp_path, "steps.gpkg", FLAG)
        with closing(sqlite3.connect(target)) as conn:
            rows = conn.execute("SELECT triangle, class, share FROM dem_height_shares").fetchall()
        assert len(rows) == 6
        assert sorted((t, k) for t, k, s in rows if s == 0.0) == [(1, 11), (2, 11)]

    def test_vtk_line_cells_get_nodata_0_and_no_rows(self, tmp_path: Path) -> None:
        vtk = read_vtk(mesh_steps(tmp_path, "steps.vtk", FLAG).read_bytes())
        assert len(vtk.lines) == 4
        table, _, line_nodata = vtk_table(vtk)
        assert_array_equal(line_nodata, np.zeros(4))
        assert len(table.offsets) == 3

    def test_read_mesh_returns_the_height_table(self, tmp_path: Path) -> None:
        data = mesh_steps(tmp_path, "steps.vtk", FLAG).read_bytes()
        got = read_mesh(data).dem_heights
        table, edges_m, _ = vtk_table(read_vtk(data))
        assert got is not None
        assert_array_equal(got.offsets, table.offsets)
        assert_array_equal(got.first, table.first)
        assert_array_equal(got.shares, table.shares)
        assert_array_equal(got.nodata, table.nodata)
        assert got.edges_m.tobytes() == edges_m.tobytes()

    def test_ascii_files_carry_the_same_numbers(self, tmp_path: Path) -> None:
        binary = vtk_table(read_vtk(mesh_steps(tmp_path / "b", "steps.vtk", FLAG).read_bytes()))
        text = vtk_table(
            read_vtk(mesh_steps(tmp_path / "a", "steps.vtk", FLAG, "--ascii").read_bytes())
        )
        for a, b in zip(binary[0].__dict__.values(), text[0].__dict__.values(), strict=True):
            assert_array_equal(a, b)
        assert_array_equal(binary[1], text[1])
        ply_b = ply_table(mesh_steps(tmp_path / "pb", "steps.ply", FLAG).read_bytes())
        ply_a = ply_table(mesh_steps(tmp_path / "pa", "steps.ply", FLAG, "--ascii").read_bytes())
        for a, b in zip(ply_b[0].__dict__.values(), ply_a[0].__dict__.values(), strict=True):
            assert_array_equal(a, b)
        assert_array_equal(ply_b[1], ply_a[1])


# ------------------------------------------------- test 20: the source's width


class TestTheSourcesWidth:
    """Test 20: test 1's DEM registered as a local source stating
    `relative_vertical_accuracy_m = 2.5`: w = 2.5, n = 4."""

    def test_20(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        tif = write_steps(tmp_path)
        root = tmp_path / "data"
        root.mkdir()
        (root / "sources.toml").write_text(
            "[source.steps-local]\n"
            'kind = "local"\n'
            f"path = {json.dumps(str(tif))}\n"
            'credit = "c"\nlicence_note = "l"\ndescription = "d"\n'
            "relative_vertical_accuracy_m = 2.5\n"
        )
        monkeypatch.setenv("RASPUTIN_DATA", str(root))
        target, stats = tmp_path / "w.vtk", tmp_path / "s.md"
        run_ok(
            "--dem",
            "steps-local",
            "--stride",
            "10",
            "--tolerance",
            "10",
            FLAG,
            "--out",
            str(target),
            "--stats",
            str(stats),
        )
        points, triangles, table, edges_m = read_any(target)
        assert edges_m.tobytes() == edges(4, 2.5).tobytes()
        assert (edges_m[0], edges_m[-1]) == (-11.25, 11.25)
        west, east = west_and_east(points, triangles)
        assert_row(table, west, 4, [33 / 50, 17 / 50])
        assert_row(table, east, 4, [47 / 50, 3 / 50])
        vtk = read_vtk(target.read_bytes())
        assert float(file_field(vtk, CLASS_M)) == 2.5
        assert "steps-local" in file_field(vtk, CLASS_FROM)
        assert "relative_vertical_accuracy_m" in file_field(vtk, CLASS_FROM)
        assert "steps-local" in stats_row(stats.read_text(), CLASS_FROM)


# ------------------------------------------------- test 17: refusals


class TestRefusals:
    """Test 17 at the command: non-zero, before any file, naming --tolerance."""

    @pytest.mark.parametrize(
        ("args", "says"),
        [
            ((), ("--tolerance", "--triangles")),
            (("--tolerance", "0"), ("--tolerance",)),
            (("--tolerance", "200"), ("--tolerance", "401", "255")),
        ],
        ids=["no tolerance", "tolerance 0", "401 classes"],
    )
    def test_refused(self, tmp_path: Path, args: tuple[str, ...], says: tuple[str, ...]) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh",
            "--dem",
            str(write_steps(tmp_path)),
            "--stride",
            "10",
            *args,
            FLAG,
            "--out",
            str(target),
        )
        assert code != 0, output
        assert "Traceback" not in output
        for word in says:
            assert squashed(word) in squashed(output), (word, output)
        assert not target.exists()
        assert not list(tmp_path.glob("x*")), list(tmp_path.iterdir())

    @pytest.mark.parametrize("args", [(), ("--tolerance", "0")], ids=["stride only", "tolerance 0"])
    def test_without_the_flag_they_run_as_today(
        self, tmp_path: Path, args: tuple[str, ...]
    ) -> None:
        run_ok(
            "--dem",
            str(write_steps(tmp_path)),
            "--stride",
            "10",
            *args,
            "--out",
            str(tmp_path / "x.vtk"),
        )


# ------------------------------------------------- test 16: the option


class TestOption:
    def test_the_flag_is_in_mesh_help(self) -> None:
        flat = " ".join(help_text("mesh").split())
        assert FLAG in flat
        for name in NAMES:
            assert name in flat, name

    def test_the_flag_is_an_mcp_argument_off_by_default(self) -> None:
        pytest.importorskip("mcp")
        from tin_engine.cli import app
        from tin_engine.mcp import schema as mcp_schema

        tools = {d.tool_name: d for d in mcp_schema.derive(app)}
        properties = tools["mesh"].tool().input_schema["properties"]
        assert "dem_heights" in properties, sorted(properties)
        assert properties["dem_heights"]["type"] == "boolean"
        assert properties["dem_heights"]["default"] is False


# ------------------------------------------------- the --triangles path


class TestTriangles:
    """3.1 and 3.5: under `--triangles` the classes follow the tolerance the
    search chose, as the file states it (`tolerance_m`; measured at 98fe80c3:
    2.05 m, so n = 2)."""

    def test_the_edges_follow_the_chosen_tolerance(self, tmp_path: Path) -> None:
        dem = write_dem(tmp_path / "hills.tif", hills())
        target = tmp_path / "t.vtk"
        run_ok(
            "--dem", str(dem), "--stride", "10", "--triangles", "300", FLAG, "--out", str(target)
        )
        vtk = read_vtk(target.read_bytes())
        tau = float(file_field(vtk, "tolerance_m"))
        table, edges_m, _ = vtk_table(vtk)
        half = half_of(tau, 1.0)
        assert half == 2, tau
        assert edges_m.tobytes() == edges(half, 1.0).tobytes()
        assert (table.count <= 2 * half + 1).all()
        assert (np.abs(row_sums(table) - 1.0) <= (table.count + 1) * 2.0**-24).all()


# ------------------------------------------------- tests 10 and 13: the cut run


def hills_index(points: np.ndarray) -> np.ndarray:
    """`hills()`'s index space: node (0, 0) at the tie point, 10 m nodes."""
    return np.column_stack([(points[:, 0] - TIE_X) / H + 0.5, (TIE_Y - points[:, 1]) / H + 0.5])


@pytest.fixture(scope="module")
def cut_run(tmp_path_factory: pytest.TempPathFactory) -> Cut:
    where = tmp_path_factory.mktemp("cut")
    return cut(where, *square_scene(where), "--pieces", "4", FLAG)


class TestCutRun:
    """Test 10: each piece file's table, from the one source tile, n = 1."""

    def test_four_piece_files_and_1163_triangles(self, cut_run: Cut) -> None:
        assert len(cut_run.pieces) == 4
        assert sum(len(p.triangles) for p in cut_run.pieces) == 1163

    def test_every_piece_file_states_the_edges(self, cut_run: Cut) -> None:
        for piece in cut_run.pieces:
            _, edges_m, _ = vtk_table(piece.vtk)
            assert edges_m.tolist() == [-1.5, -0.5, 0.5, 1.5], piece.path.name

    def test_every_triangle_of_every_piece_against_the_oracle(self, cut_run: Cut) -> None:
        dem = hills()
        for piece in cut_run.pieces:
            table, _, _ = vtk_table(piece.vtk)
            index = hills_index(piece.points)[piece.triangles]
            z = piece.points[piece.triangles][:, :, 2]
            want = oracle_table(
                index, z, dem, None, np.full(dem.shape[0], H * H), class_m=1.0, half=1
            )
            assert_array_equal(table.first, want.first, err_msg=piece.path.name)
            assert_array_equal(table.offsets, want.offsets, err_msg=piece.path.name)
            assert np.abs(table.shares - want.shares).max() <= CUT_RUN, piece.path.name
            assert np.abs(table.nodata - want.nodata).max() <= CUT_RUN, piece.path.name
            assert (np.abs(row_sums(table) - 1.0) <= (table.count + 1) * 2.0**-24).all()
            assert (table.count <= 3).all()


def height_arrays(vtk: VtkFile) -> dict[str, bytes]:
    out = {
        n: np.asarray(vtk.field_data[n].values).tobytes() for n in NAMES if n != "dem_nodata_share"
    }
    out["dem_nodata_share"] = np.asarray(vtk.cell_array("dem_nodata_share").values).tobytes()
    return out


class TestThreads:
    """Test 13 and invariant 5 at the command: the cut run at 1, 2 and 7
    threads writes the same height arrays."""

    def test_equal_arrays(self, tmp_path: Path) -> None:
        found = []
        for threads in ("1", "2", "7"):
            where = tmp_path / threads
            done = cut(where, *square_scene(where), "--pieces", "4", FLAG, "--threads", threads)
            found.append({p.path.name: height_arrays(p.vtk) for p in done.pieces})
        assert found[0] == found[1] == found[2]


# ------------------------------------------------- test 11: reprojected path


@pytest.fixture(scope="module")
def velhas(tmp_path_factory: pytest.TempPathFactory) -> VtkFile:
    target = tmp_path_factory.mktemp("velhas") / "velhas.vtk"
    run_ok(*VELHAS_ARGS, FLAG, "--out", str(target))
    return read_vtk(target.read_bytes())


@pytest.fixture(scope="module")
def oracle(velhas: VtkFile) -> tuple[Table, np.ndarray]:
    """Test 11's oracle on the source's own pixels, and the plan areas."""
    points, triangles = np.asarray(velhas.points), polygons_as_array(velhas)
    source = velhas_source()
    index = source.index(points[:, :2], "EPSG:31983")
    want = oracle_table(
        index[triangles],
        points[triangles][:, :, 2],
        source.dem,
        source.nodata,
        source.row_areas(),
        class_m=1.0,
        half=5,
    )
    return want, plan_areas(points, triangles)


class TestReprojected:
    """Test 11 and invariant 4: the shares come from the source's pixels."""

    def test_969_triangles_and_n_5(self, velhas: VtkFile) -> None:
        table, edges_m, _ = vtk_table(velhas)
        assert len(table.first) == 969
        assert edges_m.tobytes() == edges(5, 1.0).tobytes()

    def test_the_oracle_on_the_sources_pixels(
        self, velhas: VtkFile, oracle: tuple[Table, np.ndarray]
    ) -> None:
        table, _, _ = vtk_table(velhas)
        want, _ = oracle
        assert_array_equal(table.first, want.first)
        assert_array_equal(table.offsets, want.offsets)
        assert np.abs(table.shares - want.shares).max() <= REPROJECTED
        assert np.abs(table.nodata - want.nodata).max() <= REPROJECTED

    def test_spread_and_class_0_share_are_the_sources(
        self, velhas: VtkFile, oracle: tuple[Table, np.ndarray]
    ) -> None:
        table, _, _ = vtk_table(velhas)
        spread, zero = spread_and_zero(table, 5, 1.0, oracle[1])
        assert abs(spread - 1.164) <= 0.01
        assert abs(zero - 0.228) <= 0.005


# ------------------------------------------------- test 21: convergence


@pytest.fixture(scope="module")
def ladder(tmp_path_factory: pytest.TempPathFactory) -> dict[float, Table]:
    where = tmp_path_factory.mktemp("ladder")
    dem = write_dem(where / "hills.tif", hills(n=201))
    out: dict[float, Table] = {}
    for tolerance in (4.0, 2.0, 1.0, 0.4):
        target = where / f"h{tolerance}.vtk"
        run_ok(
            "--dem",
            str(dem),
            "--stride",
            "200",
            "--tolerance",
            str(tolerance),
            FLAG,
            "--out",
            str(target),
        )
        out[tolerance] = vtk_table(read_vtk(target.read_bytes()))[0]
    return out


class TestConvergence:
    """Test 21 and invariant 9 (probe: 1,546, 11,279, 42,457, 65,379
    triangles; one class 0.6, 20.6, 73.1, 100.0 %; mean count 4.875, 2.492,
    1.283, 1.000)."""

    def test_the_triangle_counts(self, ladder: dict[float, Table]) -> None:
        assert [len(t.first) for t in ladder.values()] == [1546, 11279, 42457, 65379]

    def test_the_mean_count_strictly_falls(self, ladder: dict[float, Table]) -> None:
        means = [float(t.count.mean()) for t in ladder.values()]
        assert all(a > b for a, b in pairwise(means)), means

    def test_the_one_class_share_never_falls(self, ladder: dict[float, Table]) -> None:
        ones = [float((t.count == 1).mean()) for t in ladder.values()]
        assert all(a <= b for a, b in pairwise(ones)), ones

    def test_every_count_is_at_most_2n_plus_1(self, ladder: dict[float, Table]) -> None:
        for tolerance, table in ladder.items():
            assert table.count.max() <= 2 * half_of(tolerance, 1.0) + 1, tolerance

    def test_below_half_a_class_every_triangle_has_one(self, ladder: dict[float, Table]) -> None:
        table = ladder[0.4]
        assert_array_equal(table.count, 1)
        assert_array_equal(table.first, 0)
