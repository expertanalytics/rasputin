"""`rasputin mesh --land-cover POLYGONS [--land-cover-rules NAME|FILE.toml]`:
exact class fractions per triangle from class polygons, in the mesh file.

`docs/increments/105-cell-table.md`, section 4.5 (a vector value of
`--land-cover`: a `.gpkg`, `.geojson` or `.gml` file, or a catalogue vector
source with a legend; read through `feature_input.read_source`, moved to the
mesh CRS, cut per triangle; uncovered area is class 0; a coded
`--features-rules` beside it stays refused) and section 6, T6's CLI case.
The scene is `test_cli_mesh_land_cover_fractions.py`'s: its DEM in SIRGAS
2000 / UTM 23S, its domain, `--tolerance 3`. The class polygons are CORINE
coded (`Code_18`, legend `corine`): 311 and 512 do not fit the `uint8` the
raster route had, so the int32 class numbers are exercised end to end.

THE ORACLE recomputes every triangle's shares from the file's own vertices
with GEOS's `intersection(...).area` of the whole input polygon, moved to the
mesh CRS vertex by vertex with pyproj's own `Transformer` (the producer's
edges are straight in the mesh CRS too), never from the producer's pieces.
Bound 1e-6 on a fraction: float32 rounding is 6e-8, and GEOS's clip at these
coordinates (below 1e7 m, triangles of 100 m and up) loses under 1e-9 of a
triangle's area.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- Rows are compared by value with an absent class read as 0, so a share
  under 1e-6 left out on one side and kept on the other is not a finding.
- A vector file's record strings: `land_cover_source` is the file's name,
  `land_cover_codes` is `<file name> class number; 0 = no data, and every
  constraint line`, and no credit, licence or citation is written; a
  catalogue source is named by its id, as the raster route does.
- A vector file with no legend, an uncoded legend, or with
  `--land-cover-year`, is a usage error (exit 2) naming the option, before
  any phase is timed; a polygon whose code the meshing rules lack is a usage
  error naming the code and `--land-cover`, saying "class <code> is not in
  the meshing rules" even when the legend colours the code.
"""

from __future__ import annotations

import datetime
import os
import re
import sqlite3
from collections.abc import Iterable, Sequence
from contextlib import closing
from pathlib import Path

import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal
from pyproj import Transformer
from shapely.geometry import Polygon, box
from shapely.geometry.base import BaseGeometry

from cli_helpers import USAGE, invoke, squashed
from dataset_fixtures import write_sources
from feature_fixtures import Feat, rules, write_geojson
from fraction_fixtures import (
    MESH_CRS,
    Row,
    dem_origin,
    domain_ring,
    dominant,
    rows_of,
    write_dem,
    write_domain,
)
from gpkg_fixtures import Layer, write_gpkg
from gpkg_fixtures import Row as GpkgRow
from legend_fixtures import basenames, coded_map, legend_text, uncoded_map, write
from plyread import parse_header, read_ply
from recordread import file_field
from tin_engine.io.repository import FeatureManifest
from tin_engine.stats import PhaseClock
from vector_fixtures import vector_entry
from vtkread import VtkFile, read_vtk

TOLERANCE = "3"
ARRAYS = ("land_cover_offsets", "land_cover_classes", "land_cover_fractions")
CREDIT_FIELDS = ("land_cover_credit", "land_cover_licence_note", "land_cover_cite")
#: On a fraction: float32 rounding is 6e-8; GEOS's clip loses far less.
FRACTION = 1e-6
POLYGONS = "cover.geojson"
CATALOGUE_KEY = "cover-test"


def cover_polygons() -> list[tuple[int, BaseGeometry]]:
    """Three CORINE classes over the domain, disjoint, none on a DEM node:
    311 inside it, 512 reaching past its east side, 312 with a hole and
    reaching past its north side. The rest of the domain is uncovered."""
    x0, y0 = dem_origin()
    forest = box(x0 + 300.4, y0 - 2000.2, x0 + 1200.6, y0 - 300.8)
    water = box(x0 + 1400.3, y0 - 1700.5, x0 + 2600.0, y0 - 600.9)
    outer = box(x0 + 1400.7, y0 - 400.3, x0 + 2000.1, y0 + 500.0)
    hole = box(x0 + 1600.2, y0 - 300.6, x0 + 1800.9, y0 - 150.4)
    conifer = Polygon(outer.exterior.coords, [hole.exterior.coords])
    return [(311, forest), (512, water), (312, conifer)]


def to_wgs84(geometry: BaseGeometry) -> BaseGeometry:
    move = Transformer.from_crs(MESH_CRS, "EPSG:4326", always_xy=True)
    return shapely.transform(geometry, lambda xy: np.column_stack(move.transform(*xy.T)))


def to_mesh(geometry: BaseGeometry) -> BaseGeometry:
    move = Transformer.from_crs("EPSG:4326", MESH_CRS, always_xy=True)
    return shapely.transform(geometry, lambda xy: np.column_stack(move.transform(*xy.T)))


def write_cover_geojson(
    path: Path, polygons: Sequence[tuple[int, BaseGeometry]], crs: str | None = MESH_CRS
) -> Path:
    feats = [Feat(k + 1, g, {"Code_18": str(code)}) for k, (code, g) in enumerate(polygons)]
    return write_geojson(path, feats, crs=crs)


def write_cover_gpkg(path: Path, polygons: Sequence[tuple[int, BaseGeometry]]) -> Path:
    rows = [GpkgRow(k + 1, g, {"Code_18": str(code)}) for k, (code, g) in enumerate(polygons)]
    return write_gpkg(path, [Layer("cover", 31983, rows)])


def oracle_rows(
    points: np.ndarray, triangles: np.ndarray, polygons: Iterable[tuple[int, BaseGeometry]]
) -> list[Row]:
    """Per triangle, class to share of plan area, by GEOS; class 0 the rest."""
    polygons = list(polygons)
    rows: list[Row] = []
    for tri in triangles:
        triangle = Polygon(np.asarray(points, dtype=np.float64)[tri, :2])
        row: Row = {}
        for code, polygon in polygons:
            area = float(shapely.intersection(triangle, polygon).area)
            if area > 0.0:
                row[code] = row.get(code, 0.0) + area / triangle.area
        rest = 1.0 - sum(row.values())
        if rest > 0.0:
            row[0] = rest
        rows.append(row)
    return rows


def by_value(got: Sequence[Row], want: Sequence[Row], tolerance: float) -> list[str]:
    """Where two tables differ by more than `tolerance`, an absent class 0."""
    findings: list[str] = []
    for t, (g, w) in enumerate(zip(got, want, strict=True)):
        for k in sorted(set(g) | set(w)):
            if abs(g.get(k, 0.0) - w.get(k, 0.0)) > tolerance:
                findings.append(f"triangle {t} class {k}: {g.get(k)} against {w.get(k)}")
    return findings[:20]


def mesh(dem: Path, domain: Path, target: Path, *extra: str) -> tuple[int, str]:
    return invoke(
        "mesh",
        *("--dem", str(dem), "--domain", str(domain), "--tolerance", TOLERANCE),
        *("--out", str(target)),
        *extra,
    )


def table(vtk: VtkFile) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    for name in ARRAYS:
        assert name in vtk.field_data, f"no {name}; the file has {sorted(vtk.field_data)}"
    offsets, classes, fractions = (np.asarray(vtk.field_data[n].values) for n in ARRAYS)
    return offsets, classes, fractions.astype(np.float32)


def triangles_of(vtk: VtkFile) -> np.ndarray:
    return np.array([np.asarray(p) for p in vtk.polygons], dtype=np.int64)


@pytest.fixture
def dem(tmp_path: Path) -> Path:
    return write_dem(tmp_path / "dem.tif")


@pytest.fixture
def domain(tmp_path: Path) -> Path:
    return write_domain(tmp_path / "domain.geojson")


@pytest.fixture
def cover(tmp_path: Path) -> Path:
    return write_cover_geojson(tmp_path / POLYGONS, cover_polygons())


@pytest.fixture
def phases(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """The names of every phase timed in the run."""
    seen: list[str] = []
    original = PhaseClock.add

    def spy(self: PhaseClock, name: str, seconds: float) -> None:
        seen.append(name)
        original(self, name, seconds)

    monkeypatch.setattr(PhaseClock, "add", spy)
    return seen


def run_vtk(dem: Path, domain: Path, target: Path, polygons: Path, *extra: str) -> VtkFile:
    code, output = mesh(
        dem, domain, target, "--land-cover", str(polygons), "--land-cover-rules", "corine", *extra
    )
    assert code == 0, output
    assert "No such option" not in output, output
    return read_vtk(target.read_bytes())


# ------------------------------------------------------------ the file


class TestTheFile:
    @pytest.fixture
    def run(self, tmp_path: Path, dem: Path, domain: Path, cover: Path) -> VtkFile:
        return run_vtk(dem, domain, tmp_path / "x.vtk", cover)

    def test_the_fractions_are_geos_intersection_areas(self, run: VtkFile) -> None:
        rows = rows_of(*table(run))
        want = oracle_rows(run.points, triangles_of(run), cover_polygons())
        assert by_value(rows, want, FRACTION) == []

    def test_every_row_sums_to_one_with_classes_ascending(self, run: VtkFile) -> None:
        for t, row in enumerate(rows_of(*table(run))):
            assert abs(sum(row.values()) - 1.0) <= FRACTION, f"triangle {t}: {row}"
            assert list(row) == sorted(row), f"triangle {t}: {row}"
            assert min(row.values()) > 0.0, f"triangle {t}: {row}"

    def test_every_class_and_the_uncovered_class_0_are_met(self, run: VtkFile) -> None:
        rows = rows_of(*table(run))
        assert set().union(*rows) == {0, 311, 312, 512}
        assert sum(len(row) > 1 for row in rows) > 10, "mixed triangles are the point"

    def test_the_class_numbers_are_declared_int(self, run: VtkFile) -> None:
        """4.5: 311 and 512 do not fit `unsigned_char`."""
        offsets, classes, fractions = (run.field_data[n] for n in ARRAYS)
        assert (offsets.type_name, classes.type_name, fractions.type_name) == (
            "unsigned_int", "int", "float",
        )  # fmt: skip
        assert {311, 312, 512} <= set(np.asarray(classes.values).tolist())

    def test_the_code_is_the_largest_share_and_0_on_every_line_cell(self, run: VtkFile) -> None:
        codes = np.asarray(run.cell_array("land_cover_code").values)
        lines = len(run.lines)
        assert lines > 0, "the domain's outline is a constraint, so the file has lines"
        assert not codes[:lines].any()
        assert_array_equal(codes[lines:], [dominant(row) for row in rows_of(*table(run))])

    def test_the_record_strings_name_the_file(self, run: VtkFile) -> None:
        assert file_field(run, "land_cover_source") == POLYGONS
        assert file_field(run, "land_cover_codes") == (
            f"{POLYGONS} class number; 0 = no data, and every constraint line"
        )
        assert not set(CREDIT_FIELDS) & set(run.field_data), sorted(run.field_data)

    def test_the_ply_lists_hold_the_same_rows_declared_int(
        self, run: VtkFile, tmp_path: Path, dem: Path, domain: Path, cover: Path
    ) -> None:
        target = tmp_path / "x.ply"
        code, output = mesh(
            dem, domain, target, "--land-cover", str(cover), "--land-cover-rules", "corine"
        )
        assert code == 0, output
        data = target.read_bytes()
        face = parse_header(data).element("face")
        classes = next(p for p in face.properties if p.name == "land_cover_classes")
        assert (classes.count_type, classes.type_name) == ("uchar", "int")
        _, elements = read_ply(data)
        rows = [
            dict(zip(np.asarray(k).tolist(), np.asarray(f, dtype=np.float32).tolist(), strict=True))
            for k, f in zip(
                elements["face"]["land_cover_classes"],
                elements["face"]["land_cover_fractions"],
                strict=True,
            )
        ]
        assert rows == rows_of(*table(run))

    def test_a_gpkg_output_names_and_colours_the_codes_from_the_legend(
        self, tmp_path: Path, dem: Path, domain: Path, cover: Path
    ) -> None:
        target = tmp_path / "named.gpkg"
        code, output = mesh(
            dem, domain, target, "--land-cover", str(cover), "--land-cover-rules", "corine"
        )
        assert code == 0, output
        palette = rules("corine").legend.colours
        with closing(sqlite3.connect(target)) as conn:
            rows = conn.execute("SELECT code, land_cover, colour FROM triangles").fetchall()
        assert {r[0] for r in rows} >= {311, 512}
        for code_, name, colour in rows:
            if code_ != 0:
                assert (name, colour) == palette[code_], code_


class TestAr50SubClasses:
    """`docs/increments/108-legend-split.md`, section 5, test 4: an AR50
    class-polygon run whose forest polygon carries a tree type takes the
    sub-class code, which the rules know, so the run is not refused."""

    def test_a_sub_class_polygon_is_counted(self, tmp_path: Path, dem: Path, domain: Path) -> None:
        (_, forest), (_, water), _ = cover_polygons()
        feats = [
            Feat(1, forest, {"arealtype": 30, "treslag": 31}),
            Feat(2, water, {"arealtype": 81}),
        ]
        path = write_geojson(tmp_path / "ar50.geojson", feats, crs=MESH_CRS)
        target = tmp_path / "x.vtk"
        code, output = mesh(
            dem, domain, target, "--land-cover", str(path), "--land-cover-rules", "ar50"
        )
        assert code == 0, output
        assert {3031, 81} <= set().union(*rows_of(*table(read_vtk(target.read_bytes()))))


# ------------------------------------------------------------ the other readers and CRSs


class TestSources:
    def test_a_geopackage_gives_the_geojsons_table(
        self, tmp_path: Path, dem: Path, domain: Path, cover: Path
    ) -> None:
        gpkg = write_cover_gpkg(tmp_path / "cover.gpkg", cover_polygons())
        want = table(run_vtk(dem, domain, tmp_path / "a.vtk", cover))
        got = table(run_vtk(dem, domain, tmp_path / "b.vtk", gpkg))
        for ours, theirs in zip(got, want, strict=True):
            assert_array_equal(ours, theirs)

    def test_polygons_in_another_crs_are_moved_to_the_mesh_crs(
        self, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        """Section 3C's CRS misalignment: the polygons in WGS 84 (an RFC 7946
        GeoJSON, no `crs` member) against the mesh in UTM 23S. The oracle
        moves each polygon's vertices back with pyproj and clips with GEOS."""
        wgs84 = [(code, to_wgs84(g)) for code, g in cover_polygons()]
        path = write_cover_geojson(tmp_path / "wgs84.geojson", wgs84, crs=None)
        run = run_vtk(dem, domain, tmp_path / "x.vtk", path)
        moved = [(code, to_mesh(g)) for code, g in wgs84]
        want = oracle_rows(run.points, triangles_of(run), moved)
        assert by_value(rows_of(*table(run)), want, FRACTION) == []

    def test_a_catalogue_vector_source_uses_its_legend(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        """4.5: a catalogue vector source with a legend, here a user's local
        entry with `legend = "corine"`; no `--land-cover-rules` is given."""
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        gpkg = write_cover_gpkg(tmp_path / "cover.gpkg", cover_polygons())
        config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
        write_sources(config, {CATALOGUE_KEY: vector_entry(gpkg, "Test cover", rules="corine")})
        want = table(run_vtk(dem, domain, tmp_path / "a.vtk", gpkg))
        target = tmp_path / "b.vtk"
        code, output = mesh(dem, domain, target, "--land-cover", CATALOGUE_KEY)
        assert code == 0, output
        run = read_vtk(target.read_bytes())
        for ours, theirs in zip(table(run), want, strict=True):
            assert_array_equal(ours, theirs)
        assert file_field(run, "land_cover_source") == CATALOGUE_KEY


# ------------------------------------------------------------ a fetched map-service layer

FEATURES_KEY = "cover-layer"
FEATURES_URL = "https://maps.example.org/arcgis/rest/services/Cover/MapServer/0"
FIRST_FETCHED = datetime.date(2020, 1, 2)


class TestAFetchedLayer:
    """4.5: a catalogue vector source of kind `arcgis-features`, as `clc2018`
    is, read from `<cache>/<id>/features.gpkg`. Here a user's entry with
    `legend = "corine"` over the scene's CRS, so no network is needed."""

    @pytest.fixture
    def cache(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
        entry = {
            "kind": "arcgis-features", "url": FEATURES_URL, "crs": MESH_CRS,
            "fields": ["ID", "Code_18"], "rules": "corine",
            "credit": "Test cover, {url} (accessed on {date})", "licence_note": "Test licence",
        }  # fmt: skip
        write_sources(config, {FEATURES_KEY: entry})
        return tmp_path / "cache"

    def test_the_credit_is_the_manifests_url_and_first_fetch_date(
        self, cache: Path, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        """Review round 1, finding 1: as the raster route credits the cached
        object's URL and fetch date, a fetched layer is credited with its
        manifest's `url` and `first_fetched`, not the cache path and today."""
        where = cache / FEATURES_KEY
        where.mkdir(parents=True)
        write_cover_gpkg_layer(where / "features.gpkg", cover_polygons())
        manifest = FeatureManifest(
            source=FEATURES_KEY, url=FEATURES_URL, crs=MESH_CRS, id_field="ID",
            fields=("ID", "Code_18"), rasputin_version="0", first_fetched=FIRST_FETCHED,
        )  # fmt: skip
        (where / "manifest.json").write_text(manifest.model_dump_json(), encoding="utf-8")
        target = tmp_path / "x.vtk"
        code, output = mesh(
            dem, domain, target, "--land-cover", FEATURES_KEY, "--cache", str(cache)
        )
        assert code == 0, output
        run = read_vtk(target.read_bytes())
        assert file_field(run, "land_cover_credit") == (
            f"Test cover, {FEATURES_URL} (accessed on 2020-01-02)"
        )

    def test_an_uncached_layer_is_refused_before_any_phase_naming_the_fetch(
        self, cache: Path, tmp_path: Path, dem: Path, domain: Path, phases: list[str]
    ) -> None:
        """Review round 1, finding 2: refused in the option checks, not after
        meshing with SQLite's "unable to open database file"."""
        target = tmp_path / "x.vtk"
        code, output = mesh(
            dem, domain, target, "--land-cover", FEATURES_KEY, "--cache", str(cache)
        )
        assert code == USAGE, output
        assert phases == [], phases
        assert not target.exists()
        assert "unable to open" not in output, output
        assert "--land-cover" in output, output
        assert f"rasputinfetch{FEATURES_KEY}" in squashed(output), output


def write_cover_gpkg_layer(path: Path, polygons: Sequence[tuple[int, BaseGeometry]]) -> Path:
    """The cover as a fetch writes it: one table, `features`, keyed by `ID`."""
    rows = [GpkgRow(k + 1, g, {"Code_18": str(code)}) for k, (code, g) in enumerate(polygons)]
    return write_gpkg(path, [Layer("features", 31983, rows, pk="ID")])


# ------------------------------------------------------------ refusals


class TestRefusals:
    def refused(
        self, dem: Path, domain: Path, tmp_path: Path, *extra: str, before: list[str] | None = None
    ) -> str:
        target = tmp_path / "x.vtk"
        code, output = mesh(dem, domain, target, *extra)
        assert code == USAGE, output
        assert "Nosuchoption" not in squashed(output), output
        assert "Traceback" not in squashed(output), output
        assert not target.exists()
        if before is not None:
            assert before == [], before
        return output

    def test_a_code_the_legend_lacks_is_refused_naming_it(
        self, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        x0, y0 = dem_origin()
        odd = [*cover_polygons(), (9999, box(x0 + 300.4, y0 - 2200.0, x0 + 900.0, y0 - 2100.0))]
        path = write_cover_geojson(tmp_path / "odd.geojson", odd)
        output = self.refused(
            dem, domain, tmp_path, "--land-cover", str(path), "--land-cover-rules", "corine"
        )
        assert "9999" in squashed(output) and "--land-cover" in output, output

    def test_a_code_the_legend_colours_but_the_rules_lack_names_the_rules(
        self, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        """The legend colours 5 and 6, the rules list only 5 (an unlisted
        class is land for the lines), so the code-6 polygon reaches the
        fractions and is refused there, by the rules."""
        colours = write(
            tmp_path / "five-six.toml",
            legend_text("fs", "p", [(5, "Five", "#1f8d49"), (6, "Six", "#2532e4")]),
        )
        classes = '[[class]]\ncode = 5\nlabel = "Five"\ncolour = "#1f8d49"\n'
        text = coded_map(attribute="Code_18", classes=classes, unlisted="land")
        text = re.sub(r'^legend = ".*"$', f'legend = "{colours.name}"', text, flags=re.M)
        rules_file = write(tmp_path / "five.toml", text)
        (_, forest), (_, water), _ = cover_polygons()
        path = write_cover_geojson(tmp_path / "five-six.geojson", [(5, forest), (6, water)])
        output = self.refused(
            dem, domain, tmp_path, "--land-cover", str(path), "--land-cover-rules", str(rules_file)
        )
        assert "class 6 is not in the meshing rules" in " ".join(output.split()), output
        assert "legend" not in basenames(output.replace("--land-cover-rules", "")), output

    def test_a_file_without_a_legend_is_refused_before_any_phase(
        self, tmp_path: Path, dem: Path, domain: Path, cover: Path, phases: list[str]
    ) -> None:
        output = self.refused(dem, domain, tmp_path, "--land-cover", str(cover), before=phases)
        assert "--land-cover-rules" in output, output

    def test_an_uncoded_legend_is_refused_before_any_phase(
        self, tmp_path: Path, dem: Path, domain: Path, cover: Path, phases: list[str]
    ) -> None:
        toml = tmp_path / "levees.toml"
        toml.write_text(uncoded_map(), encoding="utf-8")
        output = self.refused(
            dem, domain, tmp_path, "--land-cover", str(cover), "--land-cover-rules", str(toml),
            before=phases,
        )  # fmt: skip
        assert "--land-cover-rules" in output and "codes" in output, output

    def test_a_year_with_a_file_is_refused_before_any_phase(
        self, tmp_path: Path, dem: Path, domain: Path, cover: Path, phases: list[str]
    ) -> None:
        output = self.refused(
            dem, domain, tmp_path,
            "--land-cover", str(cover), "--land-cover-rules", "corine",
            "--land-cover-year", "2018", before=phases,
        )  # fmt: skip
        assert "--land-cover-year" in output and "no years" in output, output

    def test_a_coded_features_rules_beside_it_is_refused(
        self, tmp_path: Path, dem: Path, domain: Path, cover: Path
    ) -> None:
        """4.5: one code system per mesh, as for a class raster."""
        inner = Polygon(domain_ring()).buffer(-200.0)
        features = write_geojson(
            tmp_path / "corine.geojson", [Feat(1, inner, {"Code_18": "311"})], crs=MESH_CRS
        )
        output = self.refused(
            dem, domain, tmp_path,
            "--features", str(features), "--features-rules", "corine",
            "--land-cover", str(cover), "--land-cover-rules", "corine",
        )  # fmt: skip
        assert "--land-cover" in output and "corine" in output, output
