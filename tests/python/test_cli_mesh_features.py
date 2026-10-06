"""`rasputin mesh --dem ... --domain ... --features PATH`: increment 16b-2.

`docs/increments/16b-terrain-polygons.md` R1 (the flags), R10 (what the file
records), "Tests for @tester" (CLI; increment 8's three rows end to end; I3 on
the extract; I7 on the extract) and "Acceptance" (the quarter circle on the
committed tile with the committed extract), with Ola's Q6 (b) (the legacy GML).

Pinned by this suite (see "Pinned by the red suite (16b-1/2)"):

- Flags `--features PATH`, `--features-crs TEXT`, `--features-layer NAME`,
  `--features-map NAME` (default `property`). Each of the last three without
  `--features` is a usage error naming itself; `--features` without
  `--domain` is a usage error naming `--features`; so is any `FeatureError`,
  and no file is written.
- The `--stats` row `features` (increment 25 moved it out of the mesh file,
  `docs/increments/25-plain-output.md` D2 and D4), in the design's words:
  `<file name>[ layer <layer>], class map <map>: <n> features, <c> lines, <v>
  vertices`, where `n` is `len(FeatureSet.features)`, `c` the number of
  feature chains and `v` the number of feature vertices `start_chains` hands
  the engine (its vertices less the domain's). The layer is written for a
  GeoPackage only. The row records the data used, never how the file was
  read: it has no outside-the-domain count, which is 0 behind a spatial index
  and 1 in a full read for the same rows, so an indexed and an unindexed
  GeoPackage give the same row and the same file (the principle of Ola's
  ruling of 2026-09-28, "A table in a run report is output, not data.").
- The `--stats` rows `features_crs` (`crs_label` of the source's CRS) and
  `features_transform` (`none` for the DEM's own CRS, else
  `transform_description(source, dem)`), as `domain_crs` and
  `domain_transform`.
- `features_notice`: the map's notice, still in the mesh file (Ola's ruling
  of 2026-10-03); absent when the map has none.
- `start_mesh` reads `the domain outline and the feature lines`.
- stderr says `lines: <n> vertices read, <m> after joining shared edges and
  adding crossings`; `--stats` has phase rows `features read` and
  `features clip`.
- stderr's features line reads `features: <n> kept (<c> cut at the domain
  outline), <o> outside the domain, <e> empty` (D2's "stderr, reworded"):
  `n` counts the features kept and `c` those kept that crossed the domain
  boundary (R10: "features read, dropped outside, clipped"). How the file was
  read, including the outside count, goes to stderr only. The `--stats` row
  `features read` is a phase name and keeps its name. A GeoPackage layer
  without a spatial index is read all the same, and stderr says `<table>:
  this layer has no spatial index, so every row was read` (R3: "Without an
  index the table is scanned and the report says so").

The field without the dropped-outside count and the stderr line saying
"kept" were committed red as "16b-1/2 red: the file records data counts only;
stderr says kept": the field still carried "(<k> dropped outside)", so every
`FEATURES_FIELD` match failed and the indexed and unindexed files differed,
and stderr still said "features read".

Committed red at `e99c8ea` (amended at `972312c`): the flags did not exist
yet, so every run exited 2 with "No such option", which every usage test
excludes, and the runs failed on their exit code.
`test_the_extract_is_small_and_attributed` checks committed data and passed
already at red. The flags landed in `5079da8` and the suite has been green
since.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Polygon

import feature_fixtures as ff
from cli_driver import SQUARE, USAGE, geojson, invoke, polygon_file, rough_dem
from feature_fixtures import Feat, domain_of, write_geojson
from geotiff_fixtures import KARTVERKET, TIE_X, TIE_Y, needs_codecs
from gpkg_fixtures import (
    DTM10,
    EXTRACT,
    EXTRACT_TABLE,
    LEGACY_GML,
    OLA_EUROPE,
    OLA_NORWAY,
    Layer,
    Row,
    copy_reversed,
    needs_rtree,
    write_gpkg,
)
from landcover_fixtures import triangle_areas, vtk_labels
from plyread import read_ply
from test_cli_mesh_domain import quarter_circle
from test_cli_mesh_refine import file_field, ply_fields, stats_row
from tin_engine.crs import transform_description
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.io.geotiff import decode_dem
from vtkread import VtkFile, lines_as_array, polygons_as_array, read_vtk

V = DEFAULT_VOCABULARY
SNAP = 1e-3
FEATURES_FIELD = re.compile(
    r"(?P<name>.+?)(?: layer (?P<layer>[^,:]+))?, class map (?P<map>[a-z0-9-_]+): "
    r"(?P<n>\d+) features?, (?P<c>\d+) lines?, (?P<v>\d+) vert(?:ex|ices)"
)


def rel(x: float, y: float) -> tuple[float, float]:
    return (TIE_X + x, TIE_Y + y)


# In micro_tiff's grid (x 500 000 .. 500 200, y 6 599 920 .. 6 600 000), inside
# `SQUARE` (no hole): increment 8's three gallery rows, as features.
FOREST = Polygon([rel(30.3, -60.1), rel(60.7, -60.3), rel(60.9, -20.2), rel(30.1, -20.4)])
ROAD_INTO_FOREST = LineString([rel(20.2, -40.3), rel(45.1, -40.3)])
LAKE = Polygon([rel(120.3, -60.2), rel(170.1, -60.4), rel(170.3, -20.1), rel(120.1, -20.3)])
BRIDGE = LineString([rel(100.2, -40.7), rel(180.4, -40.7)])
WALL = LineString([rel(150.3, -15.1), rel(250.2, -15.1)])
GALLERY = [
    Feat("forest", FOREST, {"property": "land_cover"}),
    Feat("road", ROAD_INTO_FOREST, {"property": "road"}),
    Feat("lake", LAKE, {"property": "water"}),
    Feat("bridge", BRIDGE, {"property": "road"}),
    Feat("wall", WALL, {"property": "wall"}),
]


bumpy = rough_dem(16)
plain_square = polygon_file(SQUARE)


@pytest.fixture
def gallery(tmp_path: Path) -> Path:
    return write_geojson(tmp_path / "gallery.geojson", GALLERY)


def mesh(
    tmp_path: Path, tif: Path, domain: Path, *extra: str, out: str = "x.vtk", tolerance: str = "1"
) -> tuple[int, str, Path]:
    """Every run also writes ``--stats`` beside the mesh (``x.md`` for ``x.vtk``),
    unless the caller passes its own ``--stats``; ``report`` reads it."""
    target = tmp_path / out
    stats = [] if "--stats" in extra else ["--stats", str(target.with_suffix(".md"))]
    code, output = invoke(
        "mesh",
        "--dem",
        str(tif),
        "--domain",
        str(domain),
        "--tolerance",
        tolerance,
        "--out",
        str(target),
        *extra,
        *stats,
    )
    return code, output, target


def report(tmp_path: Path, out: str = "x.vtk") -> str:
    """The ``--stats`` report ``mesh`` wrote beside ``out``."""
    return (tmp_path / out).with_suffix(".md").read_text(encoding="utf-8")


def meshed(
    tmp_path: Path, tif: Path, domain: Path, *extra: str, tolerance: str = "1"
) -> tuple[VtkFile, str]:
    code, output, target = mesh(tmp_path, tif, domain, *extra, tolerance=tolerance)
    assert code == 0, output
    return read_vtk(target.read_bytes()), output


def text_field(vtk: VtkFile, name: str) -> str:
    return file_field(vtk, name)


def max_error(vtk: VtkFile) -> float:
    return float(file_field(vtk, "max_error_m"))


def edges(vtk: VtkFile) -> tuple[np.ndarray, np.ndarray]:
    """The constraint edges as `(E, 2, 2)` xy, and their masks."""
    lines = lines_as_array(vtk)
    masks = np.asarray(vtk.cell_array("feature_mask").values)[: len(lines)]
    return vtk.points[lines][:, :, :2], masks.astype(int)


def edges_at(vtk: VtkFile, point: tuple[float, float], near: float = SNAP) -> list[int]:
    xy, masks = edges(vtk)
    d = np.hypot(xy[..., 0] - point[0], xy[..., 1] - point[1]).min(axis=1)
    return sorted(int(m) for m in masks[d <= near])


# ------------------------------------------------------------------ usage


class TestUsage:
    def test_features_without_a_domain_is_refused(
        self, tmp_path: Path, bumpy: Path, gallery: Path
    ) -> None:
        code, output = invoke(
            "mesh",
            "--dem",
            str(bumpy),
            "--features",
            str(gallery),
            "--out",
            str(tmp_path / "x.vtk"),
        )
        assert code == USAGE and "--features" in output, output
        assert "No such option" not in output
        assert not (tmp_path / "x.vtk").exists()

    @pytest.mark.parametrize(
        ("flag", "value"),
        [("--features-crs", "EPSG:25833"), ("--features-layer", "x"), ("--features-map", "corine")],
    )
    def test_a_features_option_without_features_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, flag: str, value: str
    ) -> None:
        code, output, target = mesh(tmp_path, bumpy, plain_square, flag, value)
        assert code == USAGE and flag in output and "No such option" not in output, output
        assert not target.exists()

    def test_a_gallery_fixture_with_features_is_refused(
        self, tmp_path: Path, gallery: Path
    ) -> None:
        code, output = invoke(
            "mesh", "river", "--flat", "--features", str(gallery), "--out", str(tmp_path / "x.vtk")
        )
        assert code == USAGE and "--features" in output and "No such option" not in output, output

    def test_an_unknown_suffix_is_refused_naming_it(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        shp = tmp_path / "roads.shp"
        shp.write_bytes(b"")
        code, output, target = mesh(tmp_path, bumpy, plain_square, "--features", str(shp))
        assert code == USAGE and ".shp" in output and "No such option" not in output, output
        assert not target.exists()

    def test_a_layer_for_geojson_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: Path
    ) -> None:
        code, output, _ = mesh(
            tmp_path, bumpy, plain_square, "--features", str(gallery), "--features-layer", "x"
        )
        assert code == USAGE and "--features-layer" in output and "No such option" not in output, (
            output
        )

    def test_an_unknown_map_is_refused_naming_it(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: Path
    ) -> None:
        code, output, _ = mesh(
            tmp_path, bumpy, plain_square, "--features", str(gallery), "--features-map", "nosuch"
        )
        assert code == USAGE and "nosuch" in output and "No such option" not in output, output

    def test_a_feature_refusal_is_a_usage_error_naming_the_feature(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        path = write_geojson(
            tmp_path / "bad.geojson", [Feat("f-99", FOREST, {"property": "glacier"})]
        )
        code, output, target = mesh(tmp_path, bumpy, plain_square, "--features", str(path))
        assert code == USAGE and "f-99" in output and "--features" in output, output
        assert not target.exists()

    def test_a_disagreeing_features_crs_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: Path
    ) -> None:
        code, output, target = mesh(
            tmp_path, bumpy, plain_square, "--features", str(gallery), "--features-crs", "EPSG:3035"
        )
        assert code == USAGE and "No such option" not in output, output
        assert not target.exists()


# ------------------------------------------------------------ the record


class TestRecord:
    def test_the_fields(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: Path
    ) -> None:
        vtk, _ = meshed(tmp_path, bumpy, plain_square, "--features", str(gallery))
        stats = report(tmp_path)
        match = FEATURES_FIELD.fullmatch(stats_row(stats, "features"))
        assert match, stats_row(stats, "features")
        assert match["name"] == "gallery.geojson" and match["layer"] is None
        assert match["map"] == "property"
        domain = domain_of(Polygon(SQUARE))
        fs = ff.open_one(gallery, domain)
        started = ff.start(domain, fs.features)
        chains = [c for c in started.chains if c[1] == "breakline"]
        assert "outside" not in match.string
        assert int(match["n"]) == len(fs.features) == 5
        assert int(match["c"]) == len(chains)
        assert int(match["v"]) == len(started.vertices) - len(SQUARE)
        assert stats_row(stats, "features_crs") == "EPSG:25833"
        assert stats_row(stats, "features_transform").lower() == "none"
        assert "features_notice" not in vtk.field_data
        for moved in ("features", "features_crs", "features_transform"):
            assert moved not in vtk.field_data, moved

    def test_the_start_mesh_and_the_max_error(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: Path
    ) -> None:
        vtk, _ = meshed(tmp_path, bumpy, plain_square, "--features", str(gallery))
        start = stats_row(report(tmp_path), "start_mesh")
        assert start == "the domain outline and the feature lines"
        assert max_error(vtk) <= 1.0

    def test_the_vocabulary_is_the_writers_fields_and_no_edge_vocabulary(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: Path
    ) -> None:
        """R10 as ruled on 2026-09-28 ("A table in a run report is output, not
        data."): the masks' meaning is the writer's reserved fields, which
        with 16b's vocabulary pair bit 7 with `land_cover` and bit 8 with
        `water`; there is no `edge_vocabulary` field."""
        vtk, _ = meshed(tmp_path, bumpy, plain_square, "--features", str(gallery))
        bits = [int(b) for b in vtk.field_data["feature_bits"].values]
        names = [str(n) for n in vtk.field_data["feature_names"].values]
        assert len(bits) == len(names)
        pairs = set(zip(bits, names, strict=True))
        assert {(7, "land_cover"), (8, "water")} <= pairs
        assert text_field(vtk, "feature_vocabulary") == V.fingerprint()
        assert "edge_vocabulary" not in vtk.field_data

    def test_a_reprojected_source(self, tmp_path: Path, bumpy: Path, plain_square: Path) -> None:
        lonlat = [
            Feat(f.fid, ff.moved(f.geometry, "EPSG:25833", "EPSG:4326"), f.properties)  # type: ignore[arg-type]
            for f in GALLERY
        ]
        path = write_geojson(tmp_path / "wgs.geojson", lonlat, crs=None)
        meshed(tmp_path, bumpy, plain_square, "--features", str(path))
        stats = report(tmp_path)
        assert stats_row(stats, "features_crs") == "EPSG:4326"
        assert stats_row(stats, "features_transform") == transform_description(
            "EPSG:4326", "EPSG:25833"
        )

    def test_the_corine_notice(self, tmp_path: Path, bumpy: Path, plain_square: Path) -> None:
        path = write_geojson(tmp_path / "c.geojson", [Feat(1, FOREST, {"Code_18": "311"})])
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square, "--features", str(path), "--features-map", "corine"
        )
        assert "Copernicus Land Monitoring Service" in text_field(vtk, "features_notice")
        assert FEATURES_FIELD.fullmatch(stats_row(report(tmp_path), "features"))["map"] == "corine"  # type: ignore[index]

    def test_the_ply_comments(self, tmp_path: Path, bumpy: Path, plain_square: Path) -> None:
        """D2: the `.ply` carries what the `.vtk` does: the notice, not `features`."""
        path = write_geojson(tmp_path / "c.geojson", [Feat(1, FOREST, {"Code_18": "311"})])
        args = ("--features", str(path), "--features-map", "corine")
        code, output, target = mesh(tmp_path, bumpy, plain_square, *args, out="x.ply")
        assert code == 0, output
        header, _ = read_ply(target.read_bytes())
        fields = ply_fields(header.comments)
        assert "Copernicus Land Monitoring Service" in fields["features_notice"]
        for moved in ("features", "features_crs", "features_transform"):
            assert moved not in fields, moved

    def test_stderr_and_stats(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: Path
    ) -> None:
        _, output = meshed(
            tmp_path, bumpy, plain_square, "--features", str(gallery), "--stats", "-"
        )
        lines = r"\blines: \d+ vertices read, \d+ after joining shared edges and adding crossings\b"
        assert re.search(lines, output), output
        assert not re.search(r"\bnoded vertices\b", output), output
        # `invoke` collapses the output onto one line (`plain`), so a row is
        # found by its cell, `| features read |`, not by a line start.
        assert re.search(r"\|\s*features read\s*\|", output), output
        assert re.search(r"\|\s*features clip\s*\|", output), output


# ------------------------------------------------------------ the report (R3, R10)

# Two features that cross `SQUARE`'s boundary (a polygon over its west side, the
# wall over its east side), one wholly inside, one far outside: 2 clipped.
WEST_EDGE = Polygon([rel(2.1, -50.3), rel(40.2, -50.1), rel(40.4, -30.2), rel(2.3, -30.4)])
FAR_AWAY = LineString([rel(900.3, -900.1), rel(950.7, -900.2)])
CROSSINGS = [
    Feat("forest", FOREST, {"property": "land_cover"}),
    Feat("west", WEST_EDGE, {"property": "water"}),
    Feat("wall", WALL, {"property": "wall"}),
    Feat("far", FAR_AWAY, {"property": "road"}),
]
CLIPPED = re.compile(r"\((\d+) cut at the domain outline\)")
OUTSIDE = re.compile(r"\b(\d+) outside the domain\b")
FEATURES_LINE = re.compile(
    r"\bfeatures: (?P<n>\d+) kept \((?P<c>\d+) cut at the domain outline\), "
    r"(?P<o>\d+) outside the domain, (?P<e>\d+) empty\b"
)
NO_INDEX = re.compile(r"\bt: this layer has no spatial index, so every row was read\b")


def gpkg_of(path: Path, features: list[Feat], *, rtree: bool) -> Path:
    rows = [Row(i + 1, f.geometry, dict(f.properties)) for i, f in enumerate(features)]
    return write_gpkg(path, [Layer("t", 25833, rows, columns=("property",), rtree=rtree)])


class TestReport:
    def test_stderr_gives_the_clipped_count(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """R10: "features read, dropped outside, clipped". Of `CROSSINGS`,
        `west` and `wall` cross the boundary and are kept, clipped; `far` is
        dropped outside and `forest` is untouched."""
        path = write_geojson(tmp_path / "crossings.geojson", CROSSINGS)
        _, output = meshed(tmp_path, bumpy, plain_square, "--features", str(path))
        assert [int(n) for n in OUTSIDE.findall(output)] == [1], output
        assert [int(n) for n in CLIPPED.findall(output)] == [2], output

    def test_stderr_counts_the_features_kept(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """The line's first count is the features kept, and says so: of
        `CROSSINGS`, three kept, `far` dropped outside, two clipped, none
        empty. "features read" named that count wrongly."""
        path = write_geojson(tmp_path / "crossings.geojson", CROSSINGS)
        _, output = meshed(tmp_path, bumpy, plain_square, "--features", str(path))
        found = [m.groupdict() for m in FEATURES_LINE.finditer(output)]
        assert found == [{"n": "3", "c": "2", "o": "1", "e": "0"}], output
        assert not re.search(r"\bfeatures kept,", output), output

    def test_nothing_crossing_is_zero_clipped(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        path = write_geojson(tmp_path / "inside.geojson", CROSSINGS[:1])
        _, output = meshed(tmp_path, bumpy, plain_square, "--features", str(path))
        assert [int(n) for n in CLIPPED.findall(output)] == [0], output

    def test_a_layer_without_an_index_is_scanned_and_reported(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """R3: "Without an index the table is scanned and the report says so".
        The layer is still read: the features row counts the kept three, and
        the full read's outside `far` is on stderr only."""
        path = gpkg_of(tmp_path / "f.gpkg", CROSSINGS, rtree=False)
        _, output = meshed(tmp_path, bumpy, plain_square, "--features", str(path))
        features = stats_row(report(tmp_path), "features")
        match = FEATURES_FIELD.fullmatch(features)
        assert match is not None and match["n"] == "3", features
        assert [int(n) for n in OUTSIDE.findall(output)] == [1], output
        assert NO_INDEX.search(output), output
        assert not re.search(r"R-tree|scanned", output), output

    @needs_rtree
    def test_an_indexed_layer_gives_the_same_file_and_no_scan_line(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """The control: the same rows behind an R-tree give the same file,
        byte for byte, `features` field included, and nothing about a scan is
        said. The R-tree query never returns `far` where the scan reads and
        drops it; that difference is how the file was read, so it stays on
        stderr and out of the file."""
        files = []
        fields = []
        for rtree in (True, False):
            where = tmp_path / str(rtree)
            where.mkdir()
            path = gpkg_of(where / "f.gpkg", CROSSINGS, rtree=rtree)
            code, output, target = mesh(where, bumpy, plain_square, "--features", str(path))
            assert code == 0, output
            if rtree:
                assert not NO_INDEX.search(output), output
            features = stats_row(report(where), "features")
            match = FEATURES_FIELD.fullmatch(features)
            assert match is not None and match["n"] == "3", features
            fields.append(features)
            files.append(target.read_bytes())
        assert fields[0] == fields[1]
        assert files[0] == files[1]


# ------------------------------------------------- increment 8's rows, end to end


class TestGalleryRows:
    """R6: `road-enters-forest`, `wall-leaves-domain` and `bridge-over-lake`
    as GeoJSON over a synthetic DEM, through the whole command."""

    @pytest.fixture
    def vtk(self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: Path) -> VtkFile:
        return meshed(tmp_path, bumpy, plain_square, "--features", str(gallery))[0]

    def test_the_road_ends_inside_the_forest(self, vtk: VtkFile) -> None:
        west = LineString(FOREST.exterior.coords[3:5])
        cross = shapely.intersection(west, ROAD_INTO_FOREST)
        assert edges_at(vtk, (cross.x, cross.y)) == sorted(
            [V.mask("land_cover"), V.mask("land_cover"), V.mask("road"), V.mask("road")]
        )
        assert edges_at(vtk, ROAD_INTO_FOREST.coords[-1]) == [V.mask("road")]

    def test_the_bridge_is_split_at_both_shores(self, vtk: VtkFile) -> None:
        shore = LineString(LAKE.exterior.coords)
        for point in shapely.get_parts(shapely.intersection(shore, BRIDGE)):
            assert edges_at(vtk, (point.x, point.y)) == sorted(
                [V.mask("water"), V.mask("water"), V.mask("road"), V.mask("road")]
            )

    def test_the_wall_stops_at_the_domain_boundary(self, vtk: VtkFile) -> None:
        xy, masks = edges(vtk)
        wall = xy[masks & V.mask("wall") != 0]
        assert len(wall)
        domain = Polygon(SQUARE)
        points = wall.reshape(-1, 2)
        assert np.asarray(shapely.distance(domain, shapely.points(points))).max() <= SNAP
        assert np.asarray(shapely.distance(domain, shapely.points(vtk.points[:, :2]))).max() <= SNAP
        end = shapely.intersection(Polygon(SQUARE).exterior, WALL)
        assert edges_at(vtk, (end.x, end.y)) == sorted([0, 0, V.mask("wall")])

    def test_the_boundary_carries_no_feature_bit(self, vtk: VtkFile) -> None:
        xy, masks = edges(vtk)
        boundary = Polygon(SQUARE).exterior
        mids = (xy[:, 0] + xy[:, 1]) / 2
        on = np.asarray(shapely.distance(boundary, shapely.points(mids))) <= SNAP
        assert on.any() and set(masks[on].tolist()) == {0}

    def test_both_oracles_hold(self, tmp_path: Path, bumpy: Path, vtk: VtkFile) -> None:
        """Section D's two oracles on the refined output with features: the
        constrained-Delaunay property and the tolerance, from the output."""
        with bumpy.open("rb") as stream:
            tile = decode_dem(stream)
        constrained = {tuple(sorted(map(int, e))) for e in lines_as_array(vtk)}
        triangles = polygons_as_array(vtk)
        assert _delaunay_violations(vtk.points[:, :2], triangles, constrained) == []
        assert _max_error(tile, vtk.points, triangles, Polygon(SQUARE)) <= 1.0 + 1e-9


def _incircle(a: Any, b: Any, c: Any, d: Any) -> int:
    """The sign of the incircle determinant, exactly (Fractions of doubles)."""
    from fractions import Fraction

    rows = []
    for p in (a, b, c):
        dx, dy = (
            Fraction(float(p[0])) - Fraction(float(d[0])),
            Fraction(float(p[1])) - Fraction(float(d[1])),
        )
        rows.append((dx, dy, dx * dx + dy * dy))
    (a1, a2, a3), (b1, b2, b3), (c1, c2, c3) = rows
    det = a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)
    return (det > 0) - (det < 0)


def _orient(a: Any, b: Any, c: Any) -> int:
    from fractions import Fraction

    det = (Fraction(float(b[0])) - Fraction(float(a[0]))) * (
        Fraction(float(c[1])) - Fraction(float(a[1]))
    ) - (Fraction(float(b[1])) - Fraction(float(a[1]))) * (
        Fraction(float(c[0])) - Fraction(float(a[0]))
    )
    return (det > 0) - (det < 0)


def _delaunay_violations(
    xy: np.ndarray, triangles: np.ndarray, constrained: set[tuple[int, ...]]
) -> list[Any]:
    """Every interior non-constraint edge whose opposite apex lies strictly
    inside the other triangle's circumcircle, decided exactly in the file's
    frame (world coordinates as written)."""
    opposite: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for t, tri in enumerate(triangles):
        for e in range(3):
            a, b, c = int(tri[e]), int(tri[(e + 1) % 3]), int(tri[(e + 2) % 3])
            opposite.setdefault((min(a, b), max(a, b)), []).append((t, c))
    bad = []
    for edge, sides in opposite.items():
        if len(sides) != 2 or edge in constrained:
            continue
        (t, _), (_, apex) = sides
        a, b, c = (xy[int(i)] for i in triangles[t])
        if _orient(a, b, c) < 0:
            a, b = b, a
        if _incircle(a, b, c, xy[apex]) > 0:
            bad.append((edge, apex))
    return bad


def _max_error(tile: Any, points: np.ndarray, triangles: np.ndarray, domain: Polygon) -> float:
    """The largest |plane - DEM| over valid DEM nodes inside the domain, each
    in the triangle containing it, recomputed from the output."""
    m = tile.meta
    rows, cols = np.indices((m.rows, m.cols))
    x = m.x_min + cols.ravel() * m.delta_x
    y = m.y_max - rows.ravel() * m.delta_y
    z = tile.array.astype(np.float64).ravel()
    inside = np.asarray(shapely.contains_xy(domain, x, y))
    polys = [Polygon(points[t, :2]) for t in triangles]
    tree = shapely.STRtree(polys)
    worst = 0.0
    for px, py, pz in zip(x[inside], y[inside], z[inside], strict=True):
        hits = [
            int(i)
            for i in tree.query(shapely.Point(px, py))
            if polys[int(i)].covers(shapely.Point(px, py))
        ]
        assert hits, (px, py)
        p0, p1, p2 = points[triangles[hits[0]]]
        det = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])
        l1 = ((px - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (py - p0[1])) / det
        l2 = ((p1[0] - p0[0]) * (py - p0[1]) - (px - p0[0]) * (p1[1] - p0[1])) / det
        plane = p0[2] + l1 * (p1[2] - p0[2]) + l2 * (p2[2] - p0[2])
        worst = max(worst, abs(plane - pz))
    return worst


# --------------------------------------------------------- the committed data


class TestCommittedExtract:
    """Acceptance 16b-1/2 in CI: the quarter circle on the committed tile with
    the committed extract, `--features-map corine`, at 10 m (M5's reference:
    about 54 840 triangles). I7: the tolerance holds, every node covered.

    Increment 16c adds I1, I2, I3 and I5 on the same run (the `test_16c_`
    tests, committed red at `196147e`: the file had no `land_cover_code`;
    green since `0487ed0`); its I4 is
    `test_i3_rows_in_reverse_order_give_the_same_file`, unchanged."""

    @pytest.fixture(scope="class")
    def run(self, tmp_path_factory: pytest.TempPathFactory) -> tuple[VtkFile, str, Path]:
        tmp = tmp_path_factory.mktemp("extract")
        domain = geojson(tmp / "quarter.geojson", quarter_circle())
        code, output, target = mesh(
            tmp,
            KARTVERKET,
            domain,
            "--features",
            str(EXTRACT),
            "--features-map",
            "corine",
            tolerance="10",
        )
        assert code == 0, output
        return read_vtk(target.read_bytes()), output, target

    @needs_codecs
    def test_i7_the_tolerance_holds(self, run: Any) -> None:
        vtk, _, target = run
        assert max_error(vtk) <= 10.0
        assert stats_row(report(target.parent), "dem_nodes_outside_mesh") == "0"

    @needs_codecs
    def test_the_fields(self, run: Any) -> None:
        vtk, _, target = run
        stats = report(target.parent)
        match = FEATURES_FIELD.fullmatch(stats_row(stats, "features"))
        assert match, stats_row(stats, "features")
        assert match["name"] == EXTRACT.name and match["layer"] == "U2018_CLC2018_V2020_20u1"
        assert match["map"] == "corine" and int(match["n"]) > 0
        assert stats_row(stats, "features_crs") == "EPSG:3035"
        assert stats_row(stats, "features_transform") == transform_description(
            "EPSG:3035", "EPSG:25833"
        )
        assert "Copernicus Land Monitoring Service" in text_field(vtk, "features_notice")

    @needs_codecs
    def test_land_cover_everywhere_and_water_on_the_coast(self, run: Any) -> None:
        vtk, _, _ = run
        _, masks = edges(vtk)
        found = set(masks.tolist())
        assert {V.mask("land_cover"), V.mask("land_cover", "water")} <= found
        assert found <= {0, V.mask("land_cover"), V.mask("land_cover", "water")}

    @needs_codecs
    @needs_rtree
    def test_i3_rows_in_reverse_order_give_the_same_file(self, run: Any, tmp_path: Path) -> None:
        """The extract rewritten with its rows and R-tree inserted in reverse
        key order, under the same file name, gives a bit-identical `.vtk`."""
        _, _, first = run
        (tmp_path / "r").mkdir()
        reversed_path = copy_reversed(EXTRACT, tmp_path / "r" / EXTRACT.name)
        domain = geojson(tmp_path / "quarter.geojson", quarter_circle())
        code, output, target = mesh(
            tmp_path,
            KARTVERKET,
            domain,
            "--features",
            str(reversed_path),
            "--features-map",
            "corine",
            out="again.vtk",
            tolerance="10",
        )
        assert code == 0, output
        assert target.read_bytes() == first.read_bytes()

    @needs_codecs
    def test_16c_i1_i2_i3_land_cover_codes(self, run: Any) -> None:
        """Increment 16c on the real mesh: every cell carries a code, 0 on the
        lines (I3); the spread (I1); the centroid oracle against the extract's
        polygons moved to EPSG:25833 (I2). The extract covers the quarter
        circle, so no triangle is 0."""
        vtk, _, _ = run
        _, _, codes = vtk_labels(vtk, _extract_polygons())
        assert 0 not in set(codes.tolist())

    @needs_codecs
    def test_16c_i5_class_areas_match_the_clipped_extract(self, run: Any) -> None:
        """I5: per code, the triangles' xy area equals the area of each
        extract polygon (moved to EPSG:25833) intersected with the quarter
        circle, within 0.01 % of the domain's area."""
        vtk, _, _ = run
        codes = np.asarray(vtk.cell_array("land_cover_code").values)[len(lines_as_array(vtk)) :]
        triangles = polygons_as_array(vtk)
        areas = triangle_areas(vtk.points, triangles)
        domain = Polygon(quarter_circle())
        expected: dict[int, float] = {}
        for polygon, code in _extract_polygons():
            expected[code] = expected.get(code, 0.0) + polygon.intersection(domain).area
        got = {int(c): float(areas[codes == c].sum()) for c in np.unique(codes)}
        assert len(expected) >= 5
        for code in set(expected) | set(got):
            gap = abs(got.get(code, 0.0) - expected.get(code, 0.0))
            assert gap <= 1e-4 * domain.area, (code, got.get(code), expected.get(code))

    @needs_codecs
    def test_the_legacy_gml(self, tmp_path: Path) -> None:
        """Q6 (b): the legacy GML (EPSG:4326) with `clc18_kode`, same domain."""
        domain = geojson(tmp_path / "quarter.geojson", quarter_circle())
        vtk, _ = meshed(
            tmp_path,
            KARTVERKET,
            domain,
            "--features",
            str(LEGACY_GML),
            "--features-map",
            "clc18_kode",
            tolerance="10",
        )
        stats = report(tmp_path)
        assert stats_row(stats, "features_crs") == "EPSG:4326"
        match = FEATURES_FIELD.fullmatch(stats_row(stats, "features"))
        assert match and match["name"] == LEGACY_GML.name and match["layer"] is None
        assert max_error(vtk) <= 10.0
        assert V.mask("land_cover", "water") in set(edges(vtk)[1].tolist())


def _extract_polygons() -> list[tuple[Any, int]]:
    """The committed extract's polygons in the DEM's CRS, with their codes (16c)."""
    return [
        (ff.moved(geometry, "EPSG:3035"), int(code))
        for _, geometry, code in ff.extract_features(EXTRACT, EXTRACT_TABLE)
    ]


# ------------------------------------------------------------ Ola's local data


needs_ola = pytest.mark.skipif(
    not (DTM10.is_dir() and (OLA_NORWAY.exists() or OLA_EUROPE.exists())),
    reason="Ola's rasputin_data (DTM10_UTM33_20260925 and a CORINE GeoPackage) is not here",
)


@needs_ola
@needs_codecs
@pytest.mark.parametrize(
    ("gpkg", "layer"),
    [(OLA_NORWAY, "corine2018"), (OLA_EUROPE, "U2018_CLC2018_V2020_20u1")],
    ids=["norway-25833", "europe-3035"],
)
def test_olas_geopackages_over_the_dtm10_archive(tmp_path: Path, gpkg: Path, layer: str) -> None:
    """Local only: the quarter circle over the archive's tiles with Ola's
    GeoPackages (the Norwegian one has two feature layers, so the layer is
    required; its column is `code_18`, which SQLite matches case-blind)."""
    if not gpkg.exists():
        pytest.skip(f"{gpkg.name} is not here")
    domain = geojson(tmp_path / "quarter.geojson", quarter_circle())
    vtk, _ = meshed(
        tmp_path,
        DTM10,
        domain,
        "--features",
        str(gpkg),
        "--features-layer",
        layer,
        "--features-map",
        "corine",
        tolerance="10",
    )
    match = FEATURES_FIELD.fullmatch(stats_row(report(tmp_path), "features"))
    assert match and match["layer"] == layer and int(match["n"]) > 0
    assert max_error(vtk) <= 10.0
    assert V.mask("land_cover", "water") in set(edges(vtk)[1].tolist())


def test_the_extract_is_small_and_attributed() -> None:
    """ "Test data": a few hundred kB, with a NOTICE naming the source."""
    assert EXTRACT.stat().st_size < 1_000_000
    notice = (EXTRACT.parent / "NOTICE").read_text(encoding="utf-8")
    assert "Copernicus Land Monitoring Service" in notice
    assert "European Environment Agency" in notice
    assert "8e30af4" in notice and LEGACY_GML.name in notice
