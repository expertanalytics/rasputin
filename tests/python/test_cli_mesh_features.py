"""`rasputin mesh --dem ... --domain ... --features PATH`: increment 16b-2.

`docs/increments/16b-terrain-polygons.md` R1 (the flags), R10 (what the file
records), "Tests for @tester" (CLI; increment 8's three rows end to end; I3 on
the extract; I7 on the extract) and "Acceptance" (the quarter circle on the
committed tile with the committed extract), with Ola's Q6 (b) (the legacy GML).

Pinned by this suite (`docs/increments/111-one-land-cover-map.md`, sections 3
and 8, migrate it from `--features`):

- Class polygons are `--land-cover PATH --land-cover-classes NAME`, with
  `--land-cover-layer` and `--land-cover-crs` (each refused without
  `--land-cover`); lines are `--lines KEY=PATH`; `--edges all` chooses every
  class and line key. The gallery's polygons are CORINE forest and lake, its
  lines `--lines road=` and `--lines wall=`; it runs with `--edges-tolerance
  0`, its borders as given, as the uncoded gallery ran.
- `--edges` without `--domain` is a usage error saying "needs --dem, --domain
  and --tolerance" and naming `--edges`; any `FeatureError` is a usage error
  naming the option of its file, and no file is written.
- The `--stats` row `edges`: one entry per source edges were taken from,
  joined by `; `: the map `<file>[ layer <layer>], classes <C>: <n>
  features, <c> lines, <v> vertices`, a line key `<file>, lines <key>: ...`,
  where `n` is that source's count in `FeatureSet.counts`, `c` the number of
  its feature chains and `v` its feature vertices. The layer is written for a
  GeoPackage only. The row records the data used, never how the file was
  read: it has no outside-the-domain count, so an indexed and an unindexed
  GeoPackage give the same row and the same file (Ola's ruling of
  2026-09-28, "A table in a run report is output, not data.").
- The `--stats` rows `land_cover_crs` (`crs_label` of the land-cover file's
  CRS) and `land_cover_transform` (`none` for the DEM's own CRS, else
  `transform_description(source, dem)`), as `domain_crs` and
  `domain_transform`.
- `land_cover_notice`: the class file's notice, still in the mesh file.
- `start_mesh` reads `the domain outline and the feature lines`.
- stderr says `lines: <n> vertices read, <m> after joining shared edges and
  adding crossings`; `--stats` has phase rows `land cover read`, `lines
  read` and `edges`.
- stderr's edges line reads `edges: <n> kept (<c> cut at the domain
  outline), <o> outside the domain, <e> empty`. A GeoPackage layer without a
  spatial index is read all the same, and stderr says `<table>: this layer
  has no spatial index, so every row was read`.

The field without the dropped-outside count and the stderr line saying
"kept" were committed red as "16b-1/2 red: the file records data counts only;
stderr says kept": the field still carried "(<k> dropped outside)", so every
match of the row failed and the indexed and unindexed files differed,
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
from cli_helpers import SQUARE, USAGE, geojson, invoke, polygon_file, rough_dem
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
#: One entry of the `--stats` row `edges`: the map's (``classes``) or a line
#: key's (``key``).
EDGES_ENTRY = re.compile(
    r"(?P<name>.+?)(?: layer (?P<layer>[^,:]+))?, "
    r"(?:classes (?P<classes>[a-z0-9-_]+)|lines (?P<key>[a-z_]+)): "
    r"(?P<n>\d+) features?, (?P<c>\d+) lines?, (?P<v>\d+) vert(?:ex|ices)"
)


def edges_entries(stats: str) -> list[re.Match[str]]:
    """The `edges` row's entries, each matched by `EDGES_ENTRY`."""
    row = stats_row(stats, "edges")
    found = [EDGES_ENTRY.fullmatch(e) for e in row.split("; ")]
    assert all(found), row
    return [m for m in found if m is not None]


def rel(x: float, y: float) -> tuple[float, float]:
    return (TIE_X + x, TIE_Y + y)


# In micro_tiff's grid (x 500 000 .. 500 200, y 6 599 920 .. 6 600 000), inside
# `SQUARE` (no hole): increment 8's three gallery rows, as features.
FOREST = Polygon([rel(30.3, -60.1), rel(60.7, -60.3), rel(60.9, -20.2), rel(30.1, -20.4)])
ROAD_INTO_FOREST = LineString([rel(20.2, -40.3), rel(45.1, -40.3)])
LAKE = Polygon([rel(120.3, -60.2), rel(170.1, -60.4), rel(170.3, -20.1), rel(120.1, -20.3)])
BRIDGE = LineString([rel(100.2, -40.7), rel(180.4, -40.7)])
WALL = LineString([rel(150.3, -15.1), rel(250.2, -15.1)])
#: The gallery: its areas as CORINE classes (forest 311, lake 512), its lines
#: by kind.
GALLERY_COVER = [Feat("forest", FOREST, {"Code_18": "311"}), Feat("lake", LAKE, {"Code_18": "512"})]
GALLERY_LINES = {
    "road": [Feat("road", ROAD_INTO_FOREST, {}), Feat("bridge", BRIDGE, {})],
    "wall": [Feat("wall", WALL, {})],
}


def gallery_files(where: Path, crs: str | None = None) -> tuple[Path, dict[str, Path]]:
    """The gallery's land-cover file and one line file per kind, in
    EPSG:25833, or moved to ``crs`` (written with no `crs` member for
    EPSG:4326)."""

    def written(name: str, features: list[Feat]) -> Path:
        if crs is None:
            return write_geojson(where / name, features)
        lonlat = [Feat(f.fid, ff.moved(f.geometry, ff.UTM33, crs), f.properties) for f in features]  # type: ignore[arg-type]
        return write_geojson(where / name, lonlat, crs=None if crs == "EPSG:4326" else crs)

    cover = written("gallery-cover.geojson", GALLERY_COVER)
    return cover, {key: written(f"gallery-{key}.geojson", f) for key, f in GALLERY_LINES.items()}


def gallery_args(where: Path, crs: str | None = None) -> tuple[str, ...]:
    """The gallery as `rasputin mesh` options: `--land-cover` with `corine`,
    each line kind as `--lines`, `--edges all`, the band off."""
    cover, lines = gallery_files(where, crs)
    keyed = [a for key, path in lines.items() for a in ("--lines", f"{key}={path}")]
    return (
        *("--land-cover", str(cover), "--land-cover-classes", "corine"),
        *keyed,
        *("--edges", "all", "--edges-tolerance", "0"),
    )


bumpy = rough_dem(16)
plain_square = polygon_file(SQUARE)


@pytest.fixture
def gallery(tmp_path: Path) -> tuple[str, ...]:
    return gallery_args(tmp_path)


def cover_args(path: Path, classes: str = "corine", *more: str) -> tuple[str, ...]:
    """One land-cover file with ``classes`` and `--edges all`."""
    return ("--land-cover", str(path), "--land-cover-classes", classes, *more, "--edges", "all")


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
    def test_edges_without_a_domain_is_refused(self, tmp_path: Path, bumpy: Path) -> None:
        cover, _ = gallery_files(tmp_path)
        code, output = invoke(
            "mesh",
            *("--dem", str(bumpy), *cover_args(cover)),
            *("--out", str(tmp_path / "x.vtk")),
        )
        assert code == USAGE and "--edges" in output, output
        assert "needs --dem, --domain and --tolerance" in output, output
        assert "No such option" not in output
        assert not (tmp_path / "x.vtk").exists()

    @pytest.mark.parametrize(
        ("flag", "value"),
        [("--land-cover-crs", "EPSG:25833"), ("--land-cover-layer", "x")],
    )
    def test_a_land_cover_option_without_land_cover_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, flag: str, value: str
    ) -> None:
        code, output, target = mesh(tmp_path, bumpy, plain_square, flag, value)
        assert code == USAGE and flag in output and "No such option" not in output, output
        assert not target.exists()

    @pytest.mark.parametrize("option", ["--land-cover", "--lines"])
    def test_an_unknown_suffix_is_refused_naming_it(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, option: str
    ) -> None:
        shp = tmp_path / "roads.shp"
        shp.write_bytes(b"")
        given = cover_args(shp) if option == "--land-cover" else (
            "--lines", f"road={shp}", "--edges", "road"
        )  # fmt: skip
        code, output, target = mesh(tmp_path, bumpy, plain_square, *given)
        assert code == USAGE and ".shp" in output and "No such option" not in output, output
        assert option in output, output
        assert not target.exists()

    def test_a_layer_for_geojson_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        cover, _ = gallery_files(tmp_path)
        code, output, _ = mesh(
            tmp_path, bumpy, plain_square, *cover_args(cover, "corine", "--land-cover-layer", "x")
        )
        assert code == USAGE and "--land-cover-layer" in output, output
        assert "No such option" not in output

    def test_an_unknown_class_file_is_refused_naming_it(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        cover, _ = gallery_files(tmp_path)
        code, output, _ = mesh(tmp_path, bumpy, plain_square, *cover_args(cover, "nosuch"))
        assert code == USAGE and "nosuch" in output and "No such option" not in output, output

    def test_a_feature_refusal_is_a_usage_error_naming_the_feature(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        path = write_geojson(tmp_path / "bad.geojson", [Feat("f-99", FOREST, {"Code_18": "ice"})])
        code, output, target = mesh(tmp_path, bumpy, plain_square, *cover_args(path))
        assert code == USAGE and "f-99" in output and "--land-cover" in output, output
        assert not target.exists()

    def test_a_disagreeing_land_cover_crs_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        cover, _ = gallery_files(tmp_path)
        code, output, target = mesh(
            tmp_path,
            bumpy,
            plain_square,
            *cover_args(cover, "corine", "--land-cover-crs", "EPSG:3035"),
        )
        assert code == USAGE and "No such option" not in output, output
        assert "--land-cover-crs" in output, output
        assert not target.exists()


# ------------------------------------------------------------ the record


#: The CLI's clean-up for the gallery: CORINE's derived defaults, the band off.
GALLERY_CLEANUP = {"repair_m": 1.0, "merge_same_class": True, "outline_snap_m": 5.0}


class TestRecord:
    def test_the_fields(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: tuple[str, ...]
    ) -> None:
        vtk, _ = meshed(tmp_path, bumpy, plain_square, *gallery)
        stats = report(tmp_path)
        cover, road, wall = edges_entries(stats)
        assert (cover["name"], cover["layer"], cover["classes"]) == (
            "gallery-cover.geojson",
            None,
            "corine",
        )
        assert (road["name"], road["key"]) == ("gallery-road.geojson", "road")
        assert (wall["name"], wall["key"]) == ("gallery-wall.geojson", "wall")
        domain = domain_of(Polygon(SQUARE))
        paths = gallery_files(tmp_path)
        opened = [
            ff.open_one(paths[0], domain, switches=GALLERY_CLEANUP),
            ff.open_lines(paths[1]["road"], domain, "road"),
            ff.open_lines(paths[1]["wall"], domain, "wall"),
        ]
        assert "outside" not in stats_row(stats, "edges")
        every = []
        for match, fs, n in zip((cover, road, wall), opened, (2, 2, 1), strict=True):
            assert int(match["n"]) == len(fs.features) == n
            # The source's own lines and their coordinates (a closed ring
            # repeats its first point); the joined geometry's go to
            # edges_joined.
            own_lines = [line for f in fs.features for line in f.lines]
            assert int(match["c"]) == len(own_lines)
            assert int(match["v"]) == sum(int(shapely.get_num_coordinates(g)) for g in own_lines)
            every += list(fs.features)
        started = ff.start(domain, every)
        chains = [c for c in started.chains if c[1] == "breakline"]
        joined = re.fullmatch(
            r"(\d+) lines?, (\d+) vert(?:ex|ices) after joining shared edges and adding crossings",
            stats_row(stats, "edges_joined"),
        )
        assert joined, stats_row(stats, "edges_joined")
        assert int(joined[1]) == len(chains)
        assert int(joined[2]) == len(started.vertices) - len(SQUARE)
        assert stats_row(stats, "land_cover_crs") == "EPSG:25833"
        assert stats_row(stats, "land_cover_transform").lower() == "none"
        assert "Copernicus Land Monitoring Service" in text_field(vtk, "land_cover_notice")
        for moved in ("edges", "land_cover_crs", "land_cover_transform"):
            assert moved not in vtk.field_data, moved

    def test_the_start_mesh_and_the_max_error(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: tuple[str, ...]
    ) -> None:
        vtk, _ = meshed(tmp_path, bumpy, plain_square, *gallery)
        start = stats_row(report(tmp_path), "start_mesh")
        assert start == "the domain outline and the feature lines"
        assert max_error(vtk) <= 1.0

    def test_the_vocabulary_is_the_writers_fields_and_no_edge_vocabulary(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: tuple[str, ...]
    ) -> None:
        """R10 as ruled on 2026-09-28 ("A table in a run report is output, not
        data."): the masks' meaning is the writer's reserved fields, which
        with 16b's vocabulary pair bit 7 with `land_cover` and bit 8 with
        `water`; there is no `edge_vocabulary` field."""
        vtk, _ = meshed(tmp_path, bumpy, plain_square, *gallery)
        bits = [int(b) for b in vtk.field_data["feature_bits"].values]
        names = [str(n) for n in vtk.field_data["feature_names"].values]
        assert len(bits) == len(names)
        pairs = set(zip(bits, names, strict=True))
        assert {(7, "land_cover"), (8, "water")} <= pairs
        assert text_field(vtk, "feature_vocabulary") == V.fingerprint()
        assert "edge_vocabulary" not in vtk.field_data

    def test_a_reprojected_source(self, tmp_path: Path, bumpy: Path, plain_square: Path) -> None:
        meshed(tmp_path, bumpy, plain_square, *gallery_args(tmp_path, "EPSG:4326"))
        stats = report(tmp_path)
        assert stats_row(stats, "land_cover_crs") == "EPSG:4326"
        assert stats_row(stats, "land_cover_transform") == transform_description(
            "EPSG:4326", "EPSG:25833"
        )

    def test_the_corine_notice(self, tmp_path: Path, bumpy: Path, plain_square: Path) -> None:
        path = write_geojson(tmp_path / "c.geojson", [Feat(1, FOREST, {"Code_18": "311"})])
        vtk, _ = meshed(tmp_path, bumpy, plain_square, *cover_args(path))
        assert "Copernicus Land Monitoring Service" in text_field(vtk, "land_cover_notice")
        (entry,) = edges_entries(report(tmp_path))
        assert entry["classes"] == "corine"

    def test_the_ply_comments(self, tmp_path: Path, bumpy: Path, plain_square: Path) -> None:
        """D2: the `.ply` carries what the `.vtk` does: the notice, not `edges`."""
        path = write_geojson(tmp_path / "c.geojson", [Feat(1, FOREST, {"Code_18": "311"})])
        code, output, target = mesh(tmp_path, bumpy, plain_square, *cover_args(path), out="x.ply")
        assert code == 0, output
        header, _ = read_ply(target.read_bytes())
        fields = ply_fields(header.comments)
        assert "Copernicus Land Monitoring Service" in fields["land_cover_notice"]
        for moved in ("edges", "land_cover_crs", "land_cover_transform"):
            assert moved not in fields, moved

    def test_stderr_and_stats(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: tuple[str, ...]
    ) -> None:
        _, output = meshed(tmp_path, bumpy, plain_square, *gallery, "--stats", "-")
        lines = r"\blines: \d+ vertices read, \d+ after joining shared edges and adding crossings\b"
        assert re.search(lines, output), output
        assert not re.search(r"\bnoded vertices\b", output), output
        # `invoke` collapses the output onto one line (`plain`), so a row is
        # found by its cell, `| land cover read |`, not by a line start.
        for phase in ("land cover read", "lines read", "edges"):
            assert re.search(rf"\|\s*{phase}\s*\|", output), (phase, output)
        for gone in ("features read", "features clip"):
            assert not re.search(rf"\|\s*{gone}\s*\|", output), (gone, output)


# ------------------------------------------------------------ the report (R3, R10)

# Two features that cross `SQUARE`'s boundary (a polygon over its west side, the
# wall over its east side), one wholly inside, one far outside: 2 clipped. The
# polygons are the land cover, the lines one `--lines wall` file.
WEST_EDGE = Polygon([rel(2.1, -50.3), rel(40.2, -50.1), rel(40.4, -30.2), rel(2.3, -30.4)])
FAR_AWAY = LineString([rel(900.3, -900.1), rel(950.7, -900.2)])
CROSSING_COVER = [
    Feat("forest", FOREST, {"Code_18": "311"}),
    Feat("west", WEST_EDGE, {"Code_18": "512"}),
]
CROSSING_LINES = [Feat("wall", WALL, {}), Feat("far", FAR_AWAY, {})]
CLIPPED = re.compile(r"\((\d+) cut at the domain outline\)")
OUTSIDE = re.compile(r"\b(\d+) outside the domain\b")
EDGES_LINE = re.compile(
    r"\bedges: (?P<n>\d+) kept \((?P<c>\d+) cut at the domain outline\), "
    r"(?P<o>\d+) outside the domain, (?P<e>\d+) empty\b"
)
NO_INDEX = re.compile(r"\bt: this layer has no spatial index, so every row was read\b")


def crossing_args(where: Path, cover: list[Feat] | None = None) -> tuple[str, ...]:
    """`CROSSING_COVER` (or ``cover``) as the land cover, `CROSSING_LINES` as
    `--lines wall`, `--edges all`."""
    path = write_geojson(where / "crossings.geojson", CROSSING_COVER if cover is None else cover)
    lines = write_geojson(where / "crossing-lines.geojson", CROSSING_LINES)
    return (*cover_args(path), "--lines", f"wall={lines}")


def gpkg_args(where: Path, *, rtree: bool) -> tuple[str, ...]:
    """The crossings as two GeoPackages, each one table `t`: the land cover
    (`Code_18`) and the lines."""
    rows = [Row(i + 1, f.geometry, dict(f.properties)) for i, f in enumerate(CROSSING_COVER)]
    cover = write_gpkg(where / "f.gpkg", [Layer("t", 25833, rows, rtree=rtree)])
    rows = [Row(i + 1, f.geometry, {}) for i, f in enumerate(CROSSING_LINES)]
    lines = write_gpkg(where / "l.gpkg", [Layer("t", 25833, rows, columns=(), rtree=rtree)])
    return (*cover_args(cover), "--lines", f"wall={lines}")


class TestReport:
    def test_stderr_gives_the_clipped_count(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """R10: "features read, dropped outside, clipped". `west` and `wall`
        cross the boundary and are kept, clipped; `far` is dropped outside
        and `forest` is untouched."""
        _, output = meshed(tmp_path, bumpy, plain_square, *crossing_args(tmp_path))
        assert [int(n) for n in OUTSIDE.findall(output)] == [1], output
        assert [int(n) for n in CLIPPED.findall(output)] == [2], output

    def test_stderr_counts_the_features_kept(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """The line's first count is the features kept, and says so: three
        kept, `far` dropped outside, two clipped, none empty."""
        _, output = meshed(tmp_path, bumpy, plain_square, *crossing_args(tmp_path))
        found = [m.groupdict() for m in EDGES_LINE.finditer(output)]
        assert found == [{"n": "3", "c": "2", "o": "1", "e": "0"}], output
        assert not re.search(r"\bfeatures: \d+ kept\b", output), output

    def test_nothing_crossing_is_zero_clipped(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        path = write_geojson(tmp_path / "inside.geojson", CROSSING_COVER[:1])
        _, output = meshed(tmp_path, bumpy, plain_square, *cover_args(path))
        assert [int(n) for n in CLIPPED.findall(output)] == [0], output

    def test_a_layer_without_an_index_is_scanned_and_reported(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """R3: "Without an index the table is scanned and the report says so".
        The layers are still read: the edges row counts the kept three, and
        the full read's outside `far` is on stderr only."""
        _, output = meshed(tmp_path, bumpy, plain_square, *gpkg_args(tmp_path, rtree=False))
        entries = edges_entries(report(tmp_path))
        assert [int(m["n"]) for m in entries] == [2, 1], stats_row(report(tmp_path), "edges")
        assert [int(n) for n in OUTSIDE.findall(output)] == [1], output
        assert NO_INDEX.search(output), output
        assert not re.search(r"R-tree|scanned", output), output

    @needs_rtree
    def test_an_indexed_layer_gives_the_same_file_and_no_scan_line(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """The control: the same rows behind an R-tree give the same file,
        byte for byte, `edges` row included, and nothing about a scan is
        said. The R-tree query never returns `far` where the scan reads and
        drops it; that difference is how the file was read, so it stays on
        stderr and out of the file."""
        files = []
        fields = []
        for rtree in (True, False):
            where = tmp_path / str(rtree)
            where.mkdir()
            code, output, target = mesh(where, bumpy, plain_square, *gpkg_args(where, rtree=rtree))
            assert code == 0, output
            if rtree:
                assert not NO_INDEX.search(output), output
            entries = edges_entries(report(where))
            assert [int(m["n"]) for m in entries] == [2, 1], stats_row(report(where), "edges")
            fields.append(stats_row(report(where), "edges").replace(str(where), "<dir>"))
            files.append(target.read_bytes())
        assert fields[0] == fields[1]
        assert files[0] == files[1]


# ------------------------------------------------- increment 8's rows, end to end


class TestGalleryRows:
    """R6: `road-enters-forest`, `wall-leaves-domain` and `bridge-over-lake`
    as GeoJSON over a synthetic DEM, through the whole command. The lake is
    a CORINE water class, so its shore carries `land_cover` and `water`."""

    @pytest.fixture
    def vtk(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, gallery: tuple[str, ...]
    ) -> VtkFile:
        return meshed(tmp_path, bumpy, plain_square, *gallery)[0]

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
            shore_mask = V.mask("land_cover", "water")
            assert edges_at(vtk, (point.x, point.y)) == sorted(
                [shore_mask, shore_mask, V.mask("road"), V.mask("road")]
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
    the committed extract, `--land-cover-classes corine --edges all`, at 10 m (M5's reference:
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
            *cover_args(EXTRACT),
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
        (match,) = edges_entries(stats)
        assert match["name"] == EXTRACT.name and match["layer"] == "U2018_CLC2018_V2020_20u1"
        assert match["classes"] == "corine" and int(match["n"]) > 0
        assert stats_row(stats, "land_cover_crs") == "EPSG:3035"
        assert stats_row(stats, "land_cover_transform") == transform_description(
            "EPSG:3035", "EPSG:25833"
        )
        assert "Copernicus Land Monitoring Service" in text_field(vtk, "land_cover_notice")

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
            *cover_args(reversed_path),
            out="again.vtk",
            tolerance="10",
        )
        assert code == 0, output
        assert target.read_bytes() == first.read_bytes()

    @pytest.fixture(scope="class")
    def run_band_0(self, tmp_path_factory: pytest.TempPathFactory) -> VtkFile:
        """The same run with ``--edges-tolerance 0`` (increment 32, test 16:
        the default band moves borders up to 50 m from the extract's, so I2's
        oracle, which labels by the extract's polygons, needs the band off)."""
        tmp = tmp_path_factory.mktemp("extract_band_0")
        domain = geojson(tmp / "quarter.geojson", quarter_circle())
        code, output, target = mesh(
            tmp,
            KARTVERKET,
            domain,
            *cover_args(EXTRACT),
            "--edges-tolerance",
            "0",
            tolerance="10",
        )
        assert code == 0, output
        return read_vtk(target.read_bytes())

    @needs_codecs
    def test_16c_i1_i2_i3_land_cover_codes(self, run_band_0: VtkFile) -> None:
        """Increment 16c on the real mesh: every cell carries a code, 0 on the
        lines (I3); the spread (I1); the centroid oracle against the extract's
        polygons moved to EPSG:25833 (I2). The extract covers the quarter
        circle, so no triangle is 0. On the band-0 run (increment 32)."""
        _, _, codes = vtk_labels(run_band_0, _extract_polygons())
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
            *cover_args(LEGACY_GML, "clc18_kode"),
            tolerance="10",
        )
        stats = report(tmp_path)
        assert stats_row(stats, "land_cover_crs") == "EPSG:4326"
        (match,) = edges_entries(stats)
        assert match["name"] == LEGACY_GML.name and match["layer"] is None
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
        *cover_args(gpkg, "corine", "--land-cover-layer", layer),
        tolerance="10",
    )
    (match,) = edges_entries(report(tmp_path))
    assert match["layer"] == layer and int(match["n"]) > 0
    assert max_error(vtk) <= 10.0
    assert V.mask("land_cover", "water") in set(edges(vtk)[1].tolist())


def test_the_extract_is_small_and_attributed() -> None:
    """ "Test data": a few hundred kB, with a NOTICE naming the source."""
    assert EXTRACT.stat().st_size < 1_000_000
    notice = (EXTRACT.parent / "NOTICE").read_text(encoding="utf-8")
    assert "Copernicus Land Monitoring Service" in notice
    assert "European Environment Agency" in notice
    assert "8e30af4" in notice and LEGACY_GML.name in notice
