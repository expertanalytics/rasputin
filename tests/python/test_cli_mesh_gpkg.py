"""`rasputin mesh --out X.gpkg`: a GeoPackage of the mesh (increment 36).

`docs/increments/36-geopackage-output.md`, sections 4, 5 and 8: the CLI half
of T3 to T9. Two fixtures carry it:

- **F2**, `tests/fixtures/dem_archive/7908_3_10m_z33.tif` at `--tolerance 10`
  (69,571 triangles, 383 constraint edges), meshed once to `.gpkg` with
  `--record` and `--stats` beside it, once to `.vtk`, and twice more to one
  path for determinism and replacement. Needs the `codecs` extra (LZW).
- **F3**, the `bumpy` DEM with two hand-made CORINE rectangles (311 west,
  512 east) under `--land-cover-classes corine --edges all`: the coded path, with the
  `land_cover`, `land_cover_classes` and `layer_styles` tables.

The `.vtk` of a run with the same arguments is the oracle for every figure:
section 4 says `fid - 1` is the triangle's or edge's index in it, and the mesh
is deterministic. Every figure is recomputed from the arrays here, never read
from the file under test.

Pinned beyond the design's text:

- The refusal sentences are matched as "layer of the GeoPackage" with
  "lines" (for `--out-edges`) and "applies to .vtk and .ply" (for `--ascii`),
  each naming its own option.
- `gpkg_contents` of an uncoded run lists exactly `triangles`, `lines` and
  `run`; `layer_styles` is QGIS's table and is not required there.
- The QML's `<category>` elements carry `value` (the code) and `label`
  containing the palette's name; its `<symbol>` holds an
  `<Option name="color" value="r,g,b,255">` for the palette's colour.
- The `run` row for the CRS is named `crs`, as the record names it.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
import struct
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from numpy.testing import assert_allclose, assert_array_equal
from shapely.geometry import Polygon, box

from cli_helpers import SQUARE, USAGE, geojson, invoke, rough_dem
from feature_fixtures import Feat, write_geojson
from geotiff_fixtures import HAS_CODECS, KARTVERKET, TIE_X, TIE_Y, needs_codecs
from legend_snapshot import CORINE_NATURAL
from recordread import stats_rows_named
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.io.geopackage import APPLICATION_ID, GpkgLayer, layer_info, query_features
from tin_engine.io.repository import open_geopackage
from vtkread import VtkFile, lines_as_array, polygons_as_array, read_vtk

UTM33 = 25833
FIXED_TIME = "1970-01-01T00:00:00.000Z"


def rect(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    return box(TIE_X + x0, TIE_Y + y0, TIE_X + x1, TIE_Y + y1)


#: Inside `SQUARE` (x 12.3 .. 187.7, y -73.3 .. -6.7 from the tie point), off
#: every DEM node, as `test_cli_mesh_landcover.py` places them.
WEST = rect(40.3, -60.2, 100.7, -20.4)
EAST = rect(100.7, -60.2, 160.9, -20.4)
F3_POLYGONS = {311: WEST, 512: EAST}

bumpy = rough_dem(16)


@dataclass(frozen=True)
class Written:
    """One mesh written as a GeoPackage, and the `.vtk` of the same arguments."""

    gpkg: Path
    output: str
    vtk: VtkFile
    record: dict[str, Any] | None = None
    stats: str = ""

    @property
    def triangles(self) -> np.ndarray:
        return polygons_as_array(self.vtk)

    @property
    def edges(self) -> np.ndarray:
        return lines_as_array(self.vtk)

    @property
    def points(self) -> np.ndarray:
        return self.vtk.points


def mesh_pair(directory: Path, stem: str, *args: str, extra_gpkg: tuple[str, ...] = ()) -> Written:
    """`rasputin mesh <args> --out <stem>.gpkg` and the same to `.vtk`."""
    gpkg, vtk = directory / f"{stem}.gpkg", directory / f"{stem}.vtk"
    code, output = invoke("mesh", *args, "--out", str(gpkg), *extra_gpkg)
    assert code == 0, output
    code, vtk_output = invoke("mesh", *args, "--out", str(vtk))
    assert code == 0, vtk_output
    record_path = directory / f"{stem}.json"
    stats_path = directory / f"{stem}.md"
    return Written(
        gpkg,
        output,
        read_vtk(vtk.read_bytes()),
        json.loads(record_path.read_text("ascii")) if record_path.exists() else None,
        stats_path.read_text("utf-8") if stats_path.exists() else "",
    )


F2_ARGS = ("--dem", str(KARTVERKET), "--tolerance", "10")


@pytest.fixture(scope="module")
def f2(tmp_path_factory: pytest.TempPathFactory) -> Written:
    if not HAS_CODECS:
        pytest.skip("the tile is LZW; needs the `codecs` extra")
    directory = tmp_path_factory.mktemp("f2")
    extra = ("--record", str(directory / "tile.json"), "--stats", str(directory / "tile.md"))
    return mesh_pair(directory, "tile", *F2_ARGS, extra_gpkg=extra)


@dataclass(frozen=True)
class Twice:
    """F2 written twice to one path: the first file's bytes and rows, the path."""

    path: Path
    first_bytes: bytes
    first_contents: list[tuple[Any, ...]]


@pytest.fixture(scope="module")
def f2_twice(tmp_path_factory: pytest.TempPathFactory) -> Twice:
    if not HAS_CODECS:
        pytest.skip("the tile is LZW; needs the `codecs` extra")
    path = tmp_path_factory.mktemp("twice") / "tile.gpkg"
    code, output = invoke("mesh", *F2_ARGS, "--out", str(path))
    assert code == 0, output
    first = path.read_bytes()
    with closing(open_geopackage(path)) as con:
        rows = con.execute("SELECT * FROM gpkg_contents ORDER BY table_name").fetchall()
    code, output = invoke("mesh", *F2_ARGS, "--out", str(path))
    assert code == 0, output
    return Twice(path, first, rows)


@pytest.fixture
def f3(tmp_path: Path, bumpy: Path) -> Written:
    features = [
        Feat(code, polygon, {"Code_18": str(code)}) for code, polygon in F3_POLYGONS.items()
    ]
    corine = write_geojson(tmp_path / "corine.geojson", features)
    domain = geojson(tmp_path / "square.geojson", SQUARE)
    args = ("--dem", str(bumpy), "--domain", str(domain), "--tolerance", "1")
    args += ("--land-cover", str(corine), "--land-cover-classes", "corine", "--edges", "all")
    return mesh_pair(tmp_path, "cover", *args)


@pytest.fixture
def opened(f2: Written) -> Iterator[sqlite3.Connection]:
    with closing(open_geopackage(f2.gpkg)) as con:
        yield con


@pytest.fixture
def opened3(f3: Written) -> Iterator[sqlite3.Connection]:
    with closing(open_geopackage(f3.gpkg)) as con:
        yield con


def tables(con: sqlite3.Connection) -> set[str]:
    return {row[0] for row in con.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}


def registered(con: sqlite3.Connection) -> dict[str, str]:
    return dict(con.execute("SELECT table_name, data_type FROM gpkg_contents").fetchall())


def table_columns(con: sqlite3.Connection, table: str) -> dict[str, tuple[str, int]]:
    """Column name -> (declared type, primary-key position)."""
    return {row[1]: (row[2], row[5]) for row in con.execute(f"PRAGMA table_info({table})")}


def geometries(con: sqlite3.Connection, table: str, attribute: str) -> list[Any]:
    layer = layer_info(con, table)
    return list(query_features(con, layer, None, attribute))


def header_envelope(con: sqlite3.Connection, layer: GpkgLayer) -> list[tuple[float, ...]]:
    """`minx, maxx, miny, maxy` from bytes 8-40 of each blob, in `fid` order."""
    rows = con.execute(f"SELECT {layer.column} FROM {layer.table} ORDER BY {layer.pk}")
    return [struct.unpack_from("<4d", row[0], 8) for row in rows]


def envelope_of(geometry: Any) -> tuple[float, float, float, float]:
    minx, miny, maxx, maxy = geometry.bounds
    return (minx, maxx, miny, maxy)


def triangle_figures(points: np.ndarray, triangles: np.ndarray) -> dict[str, np.ndarray]:
    """The design's per-triangle figures from the arrays: corner heights,
    plan area with the half, Horn's slope `atan2(|n_xy|, |n_z|)`."""
    xyz = points[triangles]
    e1, e2 = xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]
    normal = np.cross(e1, e2)
    z = xyz[:, :, 2]
    return {
        "z_min": z.min(1),
        "z_max": z.max(1),
        "z_mean": z.mean(1),
        "area_m2": 0.5 * np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]),
        "slope_deg": np.degrees(
            np.arctan2(np.hypot(normal[:, 0], normal[:, 1]), np.abs(normal[:, 2]))
        ),
    }


def column(con: sqlite3.Connection, table: str, name: str) -> np.ndarray:
    return np.array([r[0] for r in con.execute(f"SELECT {name} FROM {table} ORDER BY fid")])


# --------------------------------------------------------------------- T8: CLI


@needs_codecs
class TestCli:
    """T8: the suffix picks the writer; the two options that make no sense
    with it are refused in plain words; the file is replaced, never appended."""

    def test_writes_the_file_and_echoes_its_path(
        self, f2: Written, opened: sqlite3.Connection
    ) -> None:
        assert f2.gpkg.is_file()
        assert str(f2.gpkg) in f2.output
        assert not f2.gpkg.with_suffix(".gpkg.part").exists(), "the .part was renamed onto --out"
        assert opened.execute("PRAGMA application_id").fetchone()[0] == APPLICATION_ID
        assert registered(opened) == {
            "triangles": "features",
            "lines": "features",
            "run": "attributes",
        }

    def test_out_edges_is_refused_because_edges_are_the_lines_layer(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        out, edges = tmp_path / "x.gpkg", tmp_path / "edges.ply"
        code, output = invoke(
            "mesh", "--dem", str(bumpy), "--out", str(out), "--out-edges", str(edges)
        )
        assert code == USAGE, output
        assert "--out-edges" in output
        assert "lines" in output and "layer of the GeoPackage" in output
        assert not out.exists() and not edges.exists()

    def test_ascii_is_refused_because_it_applies_to_vtk_and_ply(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        out = tmp_path / "x.gpkg"
        code, output = invoke("mesh", "--dem", str(bumpy), "--out", str(out), "--ascii")
        assert code == USAGE, output
        assert "--ascii" in output
        assert "applies to .vtk and .ply" in output
        assert not out.exists()

    def test_an_existing_file_is_replaced_not_appended_to(
        self, f2: Written, f2_twice: Twice
    ) -> None:
        with closing(open_geopackage(f2_twice.path)) as con:
            (count,) = con.execute("SELECT count(*) FROM triangles").fetchone()
            assert count == len(f2.triangles) == 69_571
            assert registered(con) == {
                "triangles": "features",
                "lines": "features",
                "run": "attributes",
            }

    def test_record_and_stats_work_beside_it(self, f2: Written) -> None:
        assert f2.record is not None and f2.stats
        assert str(f2.gpkg.with_suffix(".json")) in f2.output
        assert str(f2.gpkg.with_suffix(".md")) in f2.output
        assert f2.record["tolerance_m"] == 10

    def test_fid_minus_one_is_the_triangle_index_in_the_vtk_of_the_same_run(
        self, f2: Written, opened: sqlite3.Connection
    ) -> None:
        found = geometries(opened, "triangles", "code")
        assert [f.fid for f in found] == list(range(1, len(f2.triangles) + 1))
        corners = f2.points[f2.triangles]
        rings = np.concatenate([corners, corners[:, :1]], axis=1)
        got = np.array([np.asarray(f.geometry.exterior.coords) for f in found])
        assert_array_equal(got, rings)
        assert all(f.value is None for f in found), "no coded land cover: NULL codes"


# ------------------------------------------------------------- T4: determinism


@needs_codecs
class TestDeterminism:
    """T4: the same run gives the same bytes; nothing in the file is the clock."""

    def test_two_runs_give_sha256_identical_files(self, f2_twice: Twice) -> None:
        second = f2_twice.path.read_bytes()
        assert len(second) == len(f2_twice.first_bytes)
        assert (
            hashlib.sha256(second).hexdigest() == hashlib.sha256(f2_twice.first_bytes).hexdigest()
        )

    def test_last_change_is_the_fixed_text_on_every_row(self, f2_twice: Twice) -> None:
        with closing(open_geopackage(f2_twice.path)) as con:
            rows = con.execute("SELECT * FROM gpkg_contents ORDER BY table_name").fetchall()
            stamps = {r[0] for r in con.execute("SELECT last_change FROM gpkg_contents")}
        assert rows == f2_twice.first_contents
        assert stamps == {FIXED_TIME}


# --------------------------------------------------- T5: run table and classes


@needs_codecs
class TestRunTable:
    """T5: the record as rows; acceptance (d), the agent's three queries."""

    def test_every_record_entry_is_a_row_with_its_number(
        self, f2: Written, opened: sqlite3.Connection
    ) -> None:
        assert f2.record is not None
        rows = {
            name: (wording, value, number)
            for name, wording, value, number in opened.execute(
                "SELECT name, wording, value, number FROM run ORDER BY fid"
            )
        }
        assert list(rows) == list(f2.record), "record order, rasputin_version and command first"
        for name, expected in f2.record.items():
            wording, value, number = rows[name]
            assert wording, name
            if isinstance(expected, bool) or not isinstance(expected, (int, float)):
                assert (value, number) == (str(expected), None), name
            else:
                assert number == expected, name
                assert float(value) == expected, name
        assert rows["rasputin_version"][2] is None and rows["command"][2] is None
        assert table_columns(opened, "run")["fid"] == ("INTEGER", 1)

    def test_the_three_queries_agree_with_stats_and_the_mesh(
        self, f2: Written, opened: sqlite3.Connection
    ) -> None:
        stats = dict(stats_rows_named(f2.stats))
        (count,) = opened.execute("SELECT count(*) FROM triangles").fetchone()
        assert count == len(f2.triangles)
        (tolerance,) = opened.execute("SELECT value FROM run WHERE name = 'tolerance_m'").fetchone()
        assert tolerance == stats["tolerance_m"] == "10"
        (dem,) = opened.execute("SELECT value FROM run WHERE name = 'dem_source'").fetchone()
        assert dem == stats["dem_source"]
        by_code = opened.execute(
            "SELECT code, sum(area_m2) FROM triangles GROUP BY code"
        ).fetchall()
        area = triangle_figures(f2.points, f2.triangles)["area_m2"].sum()
        assert len(by_code) == 1 and by_code[0][0] is None
        # 1e-9 relative (2.5 m^2) on the 2,546 km^2 tile: two float64 sums
        # of 69,571 terms in different orders. A missing half is 1.27e9 off.
        assert by_code[0][1] == pytest.approx(area, rel=1e-9)


class TestLandCoverClasses:
    """T5 on F3: `land_cover_classes` is the legend and the class shares."""

    def test_classes_sum_to_the_triangles_and_their_area(
        self, f3: Written, opened3: sqlite3.Connection
    ) -> None:
        classes = opened3.execute(
            "SELECT code, land_cover, colour, triangles, area_m2 "
            "FROM land_cover_classes ORDER BY code"
        ).fetchall()
        by_code = opened3.execute(
            "SELECT code, count(*), sum(area_m2) FROM triangles GROUP BY code ORDER BY code"
        ).fetchall()
        assert [c[0] for c in classes] == [b[0] for b in by_code] == [0, 311, 512]
        for (code, name, colour, count, area), (_, n, total) in zip(classes, by_code, strict=True):
            assert (name, colour) == CORINE_NATURAL[code]
            assert count == n
            assert area == pytest.approx(total, rel=1e-12)
        (count,) = opened3.execute("SELECT count(*) FROM triangles").fetchone()
        assert sum(c[3] for c in classes) == count == len(f3.triangles)
        assert table_columns(opened3, "land_cover_classes")["fid"] == ("INTEGER", 1)
        assert registered(opened3)["land_cover_classes"] == "attributes"

    def test_the_plan_areas_sum_to_the_domain(self, opened3: sqlite3.Connection) -> None:
        (total,) = opened3.execute("SELECT sum(area_m2) FROM triangles").fetchone()
        # The noder snaps the ~480 m outline to 1 mm, moving the area by well
        # under 0.5 m^2 on this 11,506 m^2 square; a missing half would be
        # 5,750 m^2 off. Measured 1e-7 m^2 off at F3, about 600 triangles.
        assert total == pytest.approx(Polygon(SQUARE).area, abs=1.0)


# --------------------------------------------------------- T3: attributes (F3)


class TestLandCover:
    """T3 on F3: codes, names and colours from the palette, the figures from
    the arrays, and the `land_cover` layer of the coded polygons."""

    def test_codes_names_and_colours_follow_the_palette(self, opened3: sqlite3.Connection) -> None:
        rows = opened3.execute("SELECT code, land_cover, colour FROM triangles").fetchall()
        assert {r[0] for r in rows} == {0, 311, 512}
        for code, name, colour in rows:
            assert (name, colour) == CORINE_NATURAL[code]

    def test_the_figures_equal_those_computed_from_the_arrays(
        self, f3: Written, opened3: sqlite3.Connection
    ) -> None:
        expected = triangle_figures(f3.points, f3.triangles)
        assert_array_equal(column(opened3, "triangles", "z_min"), expected["z_min"])
        assert_array_equal(column(opened3, "triangles", "z_max"), expected["z_max"])
        assert_allclose(column(opened3, "triangles", "z_mean"), expected["z_mean"], rtol=1e-12)
        # 1e-6 m^2: UTM coordinates near 5e5 m square to 2.5e11, where float64
        # keeps about 3e-5; the cross product of corner differences loses far
        # less. Checked at F3, about 600 triangles of 1 to 100 m^2.
        assert_allclose(column(opened3, "triangles", "area_m2"), expected["area_m2"], atol=1e-6)
        assert_allclose(column(opened3, "triangles", "slope_deg"), expected["slope_deg"], atol=1e-9)
        codes = column(opened3, "triangles", "code")
        assert_array_equal(codes, f3.vtk.cell_array("land_cover_code").values[len(f3.edges) :])

    def test_the_land_cover_layer_holds_each_coded_polygon(
        self, opened3: sqlite3.Connection
    ) -> None:
        layer = layer_info(opened3, "land_cover")
        assert layer.srs_id == UTM33 and layer.rtree is None
        found = geometries(opened3, "land_cover", "code")
        assert sorted(f.value for f in found) == [311, 512]
        for feature in found:
            assert feature.geometry.geom_type == "MultiPolygon"
            assert not feature.geometry.has_z
            # The clean-up may move a vertex to the 1 mm noding grid: a
            # 200 m perimeter times 1 mm is 0.2 m^2 at most.
            assert feature.geometry.area == pytest.approx(F3_POLYGONS[feature.value].area, abs=0.5)
            assert shapely.symmetric_difference(
                feature.geometry, F3_POLYGONS[feature.value]
            ).area == pytest.approx(0.0, abs=0.5)
        rows = opened3.execute(
            "SELECT code, land_cover, colour, area_m2 FROM land_cover ORDER BY fid"
        ).fetchall()
        for (code, name, colour, area), feature in zip(rows, found, strict=True):
            assert (name, colour) == CORINE_NATURAL[code]
            assert area == pytest.approx(feature.geometry.area, rel=1e-12)
        assert header_envelope(opened3, layer) == [envelope_of(f.geometry) for f in found]
        kind = opened3.execute(
            "SELECT geometry_type_name, z FROM gpkg_geometry_columns "
            "WHERE table_name = 'land_cover'"
        ).fetchone()
        assert kind == ("MULTIPOLYGON", 0)

    @needs_codecs
    def test_an_uncoded_run_has_none_of_the_land_cover_tables(
        self, opened: sqlite3.Connection
    ) -> None:
        assert not tables(opened) & {"land_cover", "land_cover_classes", "layer_styles"}
        assert column(opened, "triangles", "code").tolist() == [None] * 69_571

    def test_a_coded_run_with_no_polygon_in_the_domain_writes_without_a_land_cover_layer(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        """A land cover whose only polygon lies outside the domain, with a
        road as `--lines`, leaves no polygon to write: the run succeeds as the
        `.vtk` run does, `land_cover` is absent (nothing to draw), while
        `land_cover_classes` and the codes on `triangles` are still written,
        every triangle in no polygon (code 0). (Was a coded source whose only
        feature was a line; lines go to `--lines`, increment 111.)"""
        road = shapely.LineString([(TIE_X + 30.1, TIE_Y - 40.3), (TIE_X + 170.2, TIE_Y - 40.3)])
        far = shapely.box(TIE_X + 900.3, TIE_Y - 950.1, TIE_X + 950.7, TIE_Y - 900.2)
        corine = write_geojson(tmp_path / "far.geojson", [Feat("f", far, {"Code_18": "311"})])
        roads = write_geojson(tmp_path / "road.geojson", [Feat("r", road, {})])
        domain = geojson(tmp_path / "square.geojson", SQUARE)
        args = ("--dem", str(bumpy), "--domain", str(domain), "--tolerance", "1")
        args += ("--land-cover", str(corine), "--land-cover-classes", "corine")
        args += ("--lines", f"road={roads}", "--edges", "all")
        written = mesh_pair(tmp_path, "road", *args)
        with closing(sqlite3.connect(written.gpkg)) as con:
            assert "land_cover" not in tables(con)
            assert "land_cover" not in registered(con)
            codes = column(con, "triangles", "code")
            assert codes.tolist() == [0] * len(written.triangles)
            assert_array_equal(
                codes, written.vtk.cell_array("land_cover_code").values[len(written.edges) :]
            )
            classes = con.execute(
                "SELECT code, land_cover, colour, triangles FROM land_cover_classes"
            ).fetchall()
            assert classes == [(0, *CORINE_NATURAL[0], len(written.triangles))]


# --------------------------------------------------------------- T6: lines (F3)


class TestLines:
    """T6 on F3: one row per constraint edge, the property columns through the
    vocabulary from the `.vtk`'s `feature_mask`."""

    def test_each_edge_reads_back_with_its_length_and_bits(
        self, f3: Written, opened3: sqlite3.Connection
    ) -> None:
        layer = layer_info(opened3, "lines")
        found = geometries(opened3, "lines", "length_m")
        edges = f3.edges
        assert len(found) == len(edges) == len(f3.vtk.lines)
        ends = f3.points[edges]
        got = np.array([np.asarray(f.geometry.coords) for f in found])
        assert_array_equal(got, ends)
        plan = np.hypot(*(ends[:, 1, :2] - ends[:, 0, :2]).T)
        assert_allclose([f.value for f in found], plan, rtol=1e-12)
        masks = np.asarray(f3.vtk.scalars["feature_mask"].values)[: len(edges)].astype(np.int64)
        used = {name for mask in np.unique(masks) for name in DEFAULT_VOCABULARY.names(int(mask))}
        assert used, "F3 has land-cover boundaries, so at least one property is set"
        assert set(table_columns(opened3, "lines")) - {layer.pk, layer.column} == used | {
            "length_m"
        }
        for name in used:
            bit = DEFAULT_VOCABULARY.mask(name)
            assert_array_equal(column(opened3, "lines", name), (masks & bit) != 0)


# ------------------------------------------------------------------ T7: SRS


class TestSrs:
    """T7 through the CLI: the run's CRS in `gpkg_spatial_ref_sys`."""

    @needs_codecs
    def test_a_dem_run_carries_its_epsg_code(self, opened: sqlite3.Connection) -> None:
        for table in ("triangles", "lines"):
            assert layer_info(opened, table).crs == f"EPSG:{UTM33}"
        org, code, definition = opened.execute(
            "SELECT organization, organization_coordsys_id, definition "
            "FROM gpkg_spatial_ref_sys WHERE srs_id = ?",
            (UTM33,),
        ).fetchone()
        assert (org, code) == ("EPSG", UTM33)
        assert definition.startswith("PROJCS[")
        ids = {r[0] for r in opened.execute("SELECT srs_id FROM gpkg_spatial_ref_sys")}
        assert ids == {-1, 0, 4326, UTM33}


# ---------------------------------------------------------------- T9: style


class TestStyle:
    """T9: QGIS's `layer_styles`, a categorised style on `code`."""

    def test_both_coded_layers_have_a_default_style_on_code(
        self, opened3: sqlite3.Connection
    ) -> None:
        rows = opened3.execute(
            "SELECT f_table_name, useAsDefault, update_time, styleQML FROM layer_styles"
        ).fetchall()
        assert {r[0] for r in rows} == {"triangles", "land_cover"}
        for _, default, stamp, qml in rows:
            assert default == 1
            assert stamp == FIXED_TIME
            root = ET.fromstring(qml)
            renderer = root.find("renderer-v2")
            assert renderer is not None
            assert renderer.get("type") == "categorizedSymbol"
            assert renderer.get("attr") == "code"
            categories = renderer.findall("categories/category")
            assert {int(c.get("value", "")) for c in categories} == set(CORINE_NATURAL)
            symbols = {s.get("name"): s for s in renderer.findall("symbols/symbol")}
            for category in categories:
                code = int(category.get("value", ""))
                name, colour = CORINE_NATURAL[code]
                assert name in category.get("label", "")
                option = symbols[category.get("symbol")].find(".//Option[@name='color']")
                assert option is not None
                rgb = tuple(int(v) for v in option.get("value", "").split(",")[:3])
                assert rgb == tuple(int(colour[i : i + 2], 16) for i in (1, 3, 5))
