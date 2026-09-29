"""`rasputin mesh ... --features-map corine`: a land-cover code per cell (16c).

`docs/increments/16c-landcover-labels.md`: R3 (what the files carry), R5 (what
`feature_input` keeps), R6 (the CLI) and "Tests for @tester", through the
command on the synthetic `bumpy` DEM of `test_cli_mesh_features.py`, with
GeoJSON features carrying `Code_18`. For every fixture: I1 (the spread, from
the file's `LINES` and triangles), I2 (the centroid oracle of
`landcover_fixtures`, from the input polygons), I3 (one value per cell, 0 on
every line), and the codes the design names for it.

Pinned beyond the design's text (the design leaves them open):

- The stderr line is matched as `land cover: <R> regions, <O> outside every
  polygon, <V> in more than one, <S> thinner than the snap` (R3's words).
- The `land_cover_codes` string contains `CORINE Land Cover level-3 code`,
  `Code_18`, `corine` and `0 =`; the PLY comment is `land_cover_codes
  <the same text>`.
- A refused code is a usage error (exit 2) naming the feature and the value,
  and `--features`, as every other `FeatureError` is (16b).
- A `LineString` under a coded map keeps its code and has no polygon.
- `rasputin mesh --help` names `land_cover_code` and `rasputin palette`.

Committed red: no code writes `land_cover_code`, `land_cover_codes` or the
stderr line; `rasputin palette` is no command; `ClassMap` has no `codes` and
`TerrainFeature` no `code` or `polygon`; a covering polygon is dropped as
outside; and `corine` accepts any value. Every test fails on one of those.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal
from shapely.geometry import LineString, Polygon, box
from shapely.geometry.base import BaseGeometry
from typer.testing import CliRunner

import feature_fixtures as ff
from feature_fixtures import Feat, domain_of, write_geojson
from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff
from landcover_fixtures import MARGIN, landcover_oracle, spread_violations, vtk_labels
from plyread import read_ply
from test_cli_mesh import plain
from test_cli_mesh_dem import USAGE, invoke, write_tiff
from test_cli_mesh_domain import COLS, ROWS, SQUARE, geojson
from tin_engine.cli import app
from vtkread import VtkFile, read_vtk

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

LAND_COVER_LINE = re.compile(
    r"land cover: (?P<r>\d+) regions, (?P<o>\d+) outside every polygon, "
    r"(?P<v>\d+) in more than one, (?P<s>\d+) thinner than the snap"
)
FEATURES_LINE = re.compile(
    r"\b(?P<n>\d+) features kept, (?P<o>\d+) dropped outside, (?P<c>\d+) clipped, "
    r"(?P<e>\d+) empty skipped\b"
)
CODES_SYSTEM = "CORINE Land Cover level-3 code"
PRESET_NAME = "rasputin CORINE natural"


def rel(x: float, y: float) -> tuple[float, float]:
    return (TIE_X + x, TIE_Y + y)


def rect(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    (ax, ay), (bx, by) = rel(x0, y0), rel(x1, y1)
    return box(ax, ay, bx, by)


def coded(fid: Any, geometry: BaseGeometry, code: Any) -> Feat:
    return Feat(fid, geometry, {"Code_18": code})


# Inside `SQUARE` (x 12.3 .. 187.7, y -73.3 .. -6.7 from the tie point), off
# every DEM node.
WEST = rect(40.3, -60.2, 100.7, -20.4)
EAST = rect(100.7, -60.2, 160.9, -20.4)
EAST_OFFSET = rect(100.7001, -60.2, 160.9, -20.4)  # 0.1 mm from WEST's east side
LAKE = rect(80.1, -50.2, 120.4, -30.3)
FOREST = rect(30.3, -65.1, 170.2, -15.3)
HOLED_FOREST = Polygon(FOREST.exterior, [LAKE.exterior])
ROAD = LineString([rel(30.1, -40.3), rel(170.2, -40.3)])
COVER = rect(-50.3, -150.1, 300.2, 50.4)
#: The domain with a notch cut down from its north side, x 80 .. 120.
NOTCHED: list[tuple[float, float]] = [
    rel(12.3, -73.3),
    rel(187.7, -72.9),
    rel(186.1, -6.7),
    rel(120.3, -6.9),
    rel(119.9, -40.1),
    rel(80.1, -40.3),
    rel(79.7, -7.0),
    rel(13.9, -7.1),
]
BAND = rect(50.2, -30.1, 150.3, -15.2)


@pytest.fixture
def bumpy(tmp_path: Path) -> Path:
    array = np.random.default_rng(16).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))


@pytest.fixture
def plain_square(tmp_path: Path) -> Path:
    return geojson(tmp_path / "square.geojson", SQUARE)


def run(
    tmp_path: Path, tif: Path, domain: Path, *extra: str, out: str = "x.vtk"
) -> tuple[int, str, Path]:
    target = tmp_path / out
    code, output = invoke(
        "--dem", str(tif), "--domain", str(domain), "--tolerance", "1", "--out", str(target), *extra
    )
    return code, output, target


@dataclass(frozen=True)
class Labelled:
    """A `.vtk` written with land-cover codes, and what I1-I3 need of it."""

    vtk: VtkFile
    output: str
    points: np.ndarray
    lines: np.ndarray
    triangles: np.ndarray
    codes: np.ndarray

    def inside(self, geometry: BaseGeometry) -> np.ndarray:
        centroids = self.points[self.triangles][:, :, :2].mean(axis=1)
        return np.asarray(shapely.contains_xy(geometry, *centroids.T))

    def stderr_counts(self) -> dict[str, int]:
        found = [m.groupdict() for m in LAND_COVER_LINE.finditer(self.output)]
        assert len(found) == 1, self.output
        return {k: int(v) for k, v in found[0].items()}


def labelled_vtk(vtk: VtkFile, output: str, polygons: list[tuple[BaseGeometry, int]]) -> Labelled:
    """I1, I2 and I3 on one file (`landcover_fixtures.vtk_labels`)."""
    lines, triangles, codes = vtk_labels(vtk, polygons, MARGIN)
    return Labelled(vtk, output, vtk.points, lines, triangles, codes)


def polygons_of(features: list[Feat]) -> list[tuple[BaseGeometry, int]]:
    return [
        (f.geometry, int(f.properties["Code_18"]))
        for f in features
        if isinstance(f.geometry, Polygon)
    ]


def mesh_corine(
    tmp_path: Path, tif: Path, features: list[Feat], domain: Path | None = None
) -> Labelled:
    path = write_geojson(tmp_path / "corine.geojson", features)
    domain = domain or geojson(tmp_path / "square.geojson", SQUARE)
    code, output, target = run(
        tmp_path, tif, domain, "--features", str(path), "--features-map", "corine"
    )
    assert code == 0, output
    return labelled_vtk(read_vtk(target.read_bytes()), output, polygons_of(features))


def text_field(vtk: VtkFile, name: str) -> str:
    (value,) = vtk.field_data[name].values
    return str(value)


# ------------------------------------------------------------ the fixtures


class TestFixtures:
    """The design's seven labelling fixtures, end to end."""

    def test_two_squares_side_by_side(self, tmp_path: Path, bumpy: Path) -> None:
        got = mesh_corine(tmp_path, bumpy, [coded(1, WEST, "311"), coded(2, EAST, "512")])
        assert set(got.codes[got.inside(WEST)].tolist()) == {311}
        assert set(got.codes[got.inside(EAST)].tolist()) == {512}
        assert set(got.codes[~got.inside(WEST.union(EAST))].tolist()) == {0}
        # The triangles along the shared side carry both codes.
        seam = LineString([rel(100.7, -55.1), rel(100.7, -25.3)])
        cells = shapely.polygons(got.points[got.triangles][:, :, :2])
        touching = np.asarray(shapely.distance(seam, cells)) <= 1e-9
        assert set(got.codes[touching].tolist()) == {311, 512}
        counts = got.stderr_counts()
        assert counts["o"] >= 1 and counts["v"] == 0

    def test_a_forest_with_an_empty_hole(self, tmp_path: Path, bumpy: Path) -> None:
        got = mesh_corine(tmp_path, bumpy, [coded("f", HOLED_FOREST, "312")])
        assert got.inside(LAKE).any()
        assert set(got.codes[got.inside(LAKE)].tolist()) == {0}
        assert set(got.codes[got.inside(HOLED_FOREST)].tolist()) == {312}

    def test_a_lake_in_the_hole_of_a_holed_forest(self, tmp_path: Path, bumpy: Path) -> None:
        got = mesh_corine(
            tmp_path, bumpy, [coded("f", HOLED_FOREST, "312"), coded("l", LAKE, "512")]
        )
        assert set(got.codes[got.inside(LAKE)].tolist()) == {512}
        assert set(got.codes[got.inside(HOLED_FOREST)].tolist()) == {312}
        assert got.stderr_counts()["v"] == 0

    def test_a_lake_in_an_unholed_forest(self, tmp_path: Path, bumpy: Path) -> None:
        """D2: the smaller polygon wins, and stderr counts the overlap."""
        got = mesh_corine(tmp_path, bumpy, [coded("f", FOREST, "312"), coded("l", LAKE, "512")])
        assert set(got.codes[got.inside(LAKE)].tolist()) == {512}
        assert set(got.codes[got.inside(FOREST.difference(LAKE))].tolist()) == {312}
        assert got.stderr_counts()["v"] > 0

    def test_a_polygon_clipped_by_the_domain_into_two_pieces(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        domain = geojson(tmp_path / "notched.geojson", NOTCHED)
        got = mesh_corine(tmp_path, bumpy, [coded("b", BAND, "324")], domain=domain)
        west = got.inside(BAND.intersection(rect(12, -80, 80.1, 0)))
        east = got.inside(BAND.intersection(rect(119.9, -80, 190, 0)))
        assert west.any() and east.any()
        assert set(got.codes[west | east].tolist()) == {324}

    def test_a_road_across_a_polygon(self, tmp_path: Path, bumpy: Path) -> None:
        got = mesh_corine(tmp_path, bumpy, [coded("f", WEST, "311"), coded("r", ROAD, "122")])
        north = got.inside(WEST.intersection(rect(0, -40.3, 200, 0)))
        south = got.inside(WEST.intersection(rect(0, -80, 200, -40.3)))
        assert north.any() and south.any()
        assert set(got.codes[north | south].tolist()) == {311}

    def test_a_polygon_covering_the_domain(self, tmp_path: Path, bumpy: Path) -> None:
        """R5: kept with no lines, counted kept and not clipped."""
        got = mesh_corine(tmp_path, bumpy, [coded("c", COVER, "333")])
        assert set(got.codes.tolist()) == {333}
        found = [m.groupdict() for m in FEATURES_LINE.finditer(got.output)]
        assert found == [{"n": "1", "o": "0", "c": "0", "e": "0"}], got.output
        assert got.stderr_counts()["o"] == 0

    def test_boundaries_a_tenth_of_a_millimetre_apart(self, tmp_path: Path, bumpy: Path) -> None:
        """Less than the snap apart: the run succeeds and, away from the seam,
        each square carries its own code (I2)."""
        got = mesh_corine(tmp_path, bumpy, [coded(1, WEST, "311"), coded(2, EAST_OFFSET, "512")])
        assert set(got.codes[got.inside(rect(40.3, -60.2, 100.6, -20.4))].tolist()) == {311}
        assert set(got.codes[got.inside(rect(100.8, -60.2, 160.9, -20.4))].tolist()) == {512}

    def test_a_code_given_as_a_json_number(self, tmp_path: Path, bumpy: Path) -> None:
        got = mesh_corine(tmp_path, bumpy, [coded(1, WEST, 311)])
        assert set(got.codes[got.inside(WEST)].tolist()) == {311}


# ------------------------------------------------------------ the record


class TestRecord:
    @pytest.fixture
    def two_squares(self, tmp_path: Path) -> Path:
        return write_geojson(
            tmp_path / "two.geojson", [coded(1, WEST, "311"), coded(2, EAST, "512")]
        )

    def test_the_codes_field(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, two_squares: Path
    ) -> None:
        code, output, target = run(
            tmp_path,
            bumpy,
            plain_square,
            "--features",
            str(two_squares),
            "--features-map",
            "corine",
        )
        assert code == 0, output
        text = text_field(read_vtk(target.read_bytes()), "land_cover_codes")
        for words in (CODES_SYSTEM, "Code_18", "corine", "0 ="):
            assert words in text, text

    @pytest.mark.parametrize("binary", [False, True], ids=["ascii", "binary"])
    def test_every_cell_carries_a_code_in_both_encodings(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, two_squares: Path, binary: bool
    ) -> None:
        extra = ["--binary"] if binary else []
        code, output, target = run(
            tmp_path,
            bumpy,
            plain_square,
            "--features",
            str(two_squares),
            "--features-map",
            "corine",
            *extra,
        )
        assert code == 0, output
        got = labelled_vtk(read_vtk(target.read_bytes()), output, [(WEST, 311), (EAST, 512)])
        assert {311, 512} <= set(got.codes.tolist())
        assert list(got.vtk.scalars) == ["feature_mask"]

    def test_no_coded_map_no_codes(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, two_squares: Path
    ) -> None:
        """D3: `property` writes no labels, nor does a run without features;
        `corine` over the same polygons does (the control)."""
        by_property = write_geojson(
            tmp_path / "p.geojson",
            [Feat(1, WEST, {"property": "land_cover"}), Feat(2, EAST, {"property": "water"})],
        )
        runs = {
            "none": (),
            "property": ("--features", str(by_property)),
            "corine": ("--features", str(two_squares), "--features-map", "corine"),
        }
        files = {}
        for name, extra in runs.items():
            (tmp_path / name).mkdir()
            code, output, target = run(tmp_path / name, bumpy, plain_square, *extra)
            assert code == 0, output
            files[name] = (read_vtk(target.read_bytes()), output)
        for name in ("none", "property"):
            vtk, output = files[name]
            with pytest.raises(KeyError):
                vtk.cell_array("land_cover_code")
            assert "land_cover_codes" not in vtk.field_data
            assert not LAND_COVER_LINE.search(output), output
        vtk, output = files["corine"]
        assert "land_cover_codes" in vtk.field_data
        assert LAND_COVER_LINE.search(output), output

    def test_the_stats_phase(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, two_squares: Path
    ) -> None:
        code, output, _ = run(
            tmp_path,
            bumpy,
            plain_square,
            "--features",
            str(two_squares),
            "--features-map",
            "corine",
            "--stats",
            "-",
        )
        assert code == 0, output
        assert re.search(r"\|\s*land cover\s*\|", output), output

    def test_the_help_names_the_array_and_the_palette(self) -> None:
        result = runner.invoke(app, ["mesh", "--help"])
        text = plain(result.output)
        assert "land_cover_code" in text and "rasputin palette" in text, text


# ------------------------------------------------------------ refusals


class TestRefusals:
    @pytest.mark.parametrize(
        "value",
        ["forest", "0", "-311", "2147483648", "3.5", ["311", "312"]],
        ids=["word", "zero", "negative", "over-int32", "decimal", "list"],
    )
    def test_a_code_that_is_not_a_class_code_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, value: Any
    ) -> None:
        """R5: `int(str(value))` in 1 .. 2**31 - 1, else refused naming the
        feature and the value; a list is refused too."""
        path = write_geojson(tmp_path / "bad.geojson", [coded("f-77", WEST, value)])
        code, output, target = run(
            tmp_path, bumpy, plain_square, "--features", str(path), "--features-map", "corine"
        )
        assert code == USAGE, output
        assert "f-77" in output and "--features" in output, output
        shown = value[0] if isinstance(value, list) else value
        assert str(shown) in output, output
        assert not target.exists()


# ------------------------------------------------------------ .ply


class TestPly:
    def test_the_face_file_carries_the_codes_and_the_edge_file_not(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        path = write_geojson(
            tmp_path / "two.geojson", [coded(1, WEST, "311"), coded(2, EAST, "512")]
        )
        edges_path = tmp_path / "e.ply"
        code, output, target = run(
            tmp_path,
            bumpy,
            plain_square,
            "--features",
            str(path),
            "--features-map",
            "corine",
            "--out-edges",
            str(edges_path),
            out="x.ply",
        )
        assert code == 0, output
        header, data = read_ply(target.read_bytes())
        face = header.element("face")
        assert [p.name for p in face.properties] == ["vertex_indices", "land_cover_code"]
        assert face.properties[1].type_name == "int"
        texts = [c for c in header.comments if c.startswith("land_cover_codes ")]
        assert len(texts) == 1 and CODES_SYSTEM in texts[0], header.comments

        edge_header, edge_data = read_ply(edges_path.read_bytes())
        assert all(p.name != "land_cover_code" for e in edge_header.elements for p in e.properties)
        assert not any(c.startswith("land_cover_codes") for c in edge_header.comments)

        # I1 and I2 on the face file, with the edge file's constraints.
        vertex = data["vertex"]
        points = np.column_stack([vertex["x"], vertex["y"], vertex["z"]])
        triangles = np.stack(list(data["face"]["vertex_indices"])).astype(np.int64)
        codes = np.asarray(data["face"]["land_cover_code"])
        lines = np.column_stack([edge_data["edge"]["vertex1"], edge_data["edge"]["vertex2"]])
        assert spread_violations(triangles, lines, codes) == []
        oracle = landcover_oracle(points, triangles, [(WEST, 311), (EAST, 512)], MARGIN)
        assert oracle.mismatches(codes) == []
        assert {311, 512} <= set(codes.tolist())


# ------------------------------------------------------------ VTK readback


def test_vtk_reads_the_array_back(tmp_path: Path, bumpy: Path, plain_square: Path) -> None:
    """`vtkPolyDataReader` sees `land_cover_code` with one value per cell,
    0 on the lines (which VTK numbers first). Skipped without `vtk`, as
    `test_io_vtk_readback.py` is; CI runs it in the `viewer` extra's step."""
    vtk = pytest.importorskip("vtk")
    path = write_geojson(tmp_path / "two.geojson", [coded(1, WEST, "311"), coded(2, EAST, "512")])
    code, output, target = run(
        tmp_path, bumpy, plain_square, "--features", str(path), "--features-map", "corine"
    )
    assert code == 0, output
    reader = vtk.vtkPolyDataReader()
    errors: list[str] = []
    reader.AddObserver("ErrorEvent", lambda _obj, _event: errors.append("error"))
    reader.SetFileName(str(target))
    reader.ReadAllFieldsOn()
    reader.Update()
    assert not errors
    polydata = reader.GetOutput()
    array = polydata.GetCellData().GetArray("land_cover_code")
    assert array is not None
    values = np.array([array.GetValue(i) for i in range(array.GetNumberOfTuples())])
    assert len(values) == polydata.GetNumberOfCells()
    lines = polydata.GetNumberOfLines()
    assert lines > 0 and not values[:lines].any()
    assert {311, 512} <= set(values[lines:].tolist())
    ours = read_vtk(target.read_bytes()).cell_array("land_cover_code").values
    assert_array_equal(values, np.asarray(ours))


# ------------------------------------------------------------ the palette


def test_palette_out_writes_the_preset(tmp_path: Path) -> None:
    """R4: `rasputin palette corine --out FILE` writes the ParaView preset."""
    from tin_engine.palettes import CORINE_NATURAL, paraview_preset

    target = tmp_path / "corine_natural.json"
    result = runner.invoke(app, ["palette", "corine", "--out", str(target)])
    assert result.exit_code == 0, result.output
    assert json.loads(target.read_text()) == paraview_preset(CORINE_NATURAL, PRESET_NAME)


# ------------------------------------------------- what feature_input keeps (R5)


class TestWhatFeatureInputKeeps:
    DOMAIN = domain_of(Polygon(SQUARE))

    def open(self, tmp_path: Path, features: list[Feat], map_name: str, **kw: Any) -> Any:
        path = write_geojson(tmp_path / f"{map_name}.geojson", features, **kw)
        return ff.open_one(path, self.DOMAIN, map_name)

    @pytest.mark.parametrize(
        ("name", "codes"),
        [
            ("corine", CODES_SYSTEM),
            ("corine-water", CODES_SYSTEM),
            ("clc18_kode", CODES_SYSTEM),
            ("property", ""),
        ],
    )
    def test_the_maps_name_their_code_system(self, name: str, codes: str) -> None:
        assert ff.feature_input().CLASS_MAPS[name].codes == codes

    def test_a_coded_polygon_keeps_its_code_and_polygon(self, tmp_path: Path) -> None:
        fs = self.open(tmp_path, [coded(1, WEST, "311")], "corine")
        (feature,) = fs.features
        assert feature.code == 311
        assert feature.polygon is not None
        assert feature.polygon.equals(WEST)

    def test_a_coded_line_keeps_its_code_and_no_polygon(self, tmp_path: Path) -> None:
        fs = self.open(tmp_path, [coded("r", ROAD, "122")], "corine")
        (feature,) = fs.features
        assert (feature.code, feature.polygon) == (122, None)

    def test_the_property_map_keeps_neither(self, tmp_path: Path) -> None:
        fs = self.open(tmp_path, [Feat(1, WEST, {"property": "land_cover"})], "property")
        (feature,) = fs.features
        assert (feature.code, feature.polygon) == (None, None)

    def test_corine_water_keeps_the_water_only(self, tmp_path: Path) -> None:
        fs = self.open(
            tmp_path, [coded("f", FOREST, "312"), coded("l", LAKE, "512")], "corine-water"
        )
        (feature,) = fs.features
        assert (feature.fid, feature.code) == ("l", 512)

    def test_a_reprojected_polygon_is_moved_into_the_dems_crs(self, tmp_path: Path) -> None:
        lonlat = ff.moved(WEST, "EPSG:25833", "EPSG:4326")
        fs = self.open(tmp_path, [coded(1, lonlat, "311")], "corine", crs=None)
        (feature,) = fs.features
        expected = ff.moved(lonlat, "EPSG:4326")
        assert feature.polygon.hausdorff_distance(expected) <= 1e-6
        assert feature.polygon.hausdorff_distance(WEST) <= 1e-3

    def test_a_covering_polygon_is_kept_and_a_disjoint_one_dropped(self, tmp_path: Path) -> None:
        far = rect(900.3, -900.1, 950.7, -850.2)
        fs = self.open(tmp_path, [coded("c", COVER, "333"), coded("x", far, "311")], "corine")
        (feature,) = fs.features
        assert (feature.fid, feature.code, feature.lines) == ("c", 333, ())
        assert feature.polygon.equals(COVER)
        assert (fs.outside, fs.clipped) == (1, 0)

    def test_an_uncoded_map_still_drops_a_covering_polygon(self, tmp_path: Path) -> None:
        """16b's behaviour for `property`; the coded control keeps it."""
        uncoded = self.open(tmp_path, [Feat("c", COVER, {"property": "land_cover"})], "property")
        assert (len(uncoded.features), uncoded.outside) == (0, 1)
        kept = self.open(tmp_path, [coded("c", COVER, "333")], "corine")
        assert len(kept.features) == 1
