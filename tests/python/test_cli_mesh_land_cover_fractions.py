"""`rasputin mesh --land-cover KEY [--land-cover-year Y]`: exact class fractions
per triangle in the mesh file, from the tile cache (26b, `cli.py`,
`land_cover_fractions.py`).

`docs/increments/26-land-cover.md`, 1.1 (the raster's own coordinates), 1.6
(what the file carries), section 6 (the blueprint), section 8 (no
`--land-cover`: nothing; a CORINE map beside it: refused) and section 10's
"26b": the `.vtk` arrays and `.ply` lists read back, `land_cover_code` 0 on
line cells, the CLI options. The scene: a projected DEM in SIRGAS 2000 / UTM
23S over the lower-left quarter of the synthetic class raster of
`class_fixtures` (EPSG:4326, cells of 1/128 degree, with a block of class 0),
the class source installed in `SOURCES` under `class_fixtures.KEY` and two of
its years written into a cache; the network closed.

THE ORACLE (`fraction_fixtures.oracle_fractions`) recomputes every triangle's
fractions from the file's own vertices with pyproj's `Transformer`, GEOS's
clip and `Geod`'s cell areas, sharing nothing with the producer. Bound 1e-6
on a fraction: float32 rounding is 6e-8, and the oracle takes the triangle's
edges straight in index space as the producer does (1.1), so the bend between
the two frames does not enter; it is the producer's claim, not this test's.
"""

from __future__ import annotations

import json
import socket
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from shapely.geometry import Polygon

from class_fixtures import (
    CITE,
    CLASSES,
    KEY,
    LICENCE,
    PATH_TEMPLATE,
    class_grid,
    object_id,
    per_year_entry,
)
from cli_driver import USAGE, invoke, plain, runner
from feature_fixtures import Feat, write_geojson
from fraction_fixtures import (
    FETCHED,
    MESH_CRS,
    Row,
    domain_ring,
    dominant,
    oracle_fractions,
    rows_of,
    same_rows,
    write_class_cache,
    write_dem,
    write_domain,
)
from plyread import read_ply
from recordread import file_field, ply_fields, stats_row
from test_cli_fetch import KEY as DEM_KEY
from test_cli_fetch import catalogue_with, install
from test_cli_mesh_cache import no_network
from tin_engine import sources
from tin_engine.cli import app
from tin_engine.io.repository import CacheRepository
from tin_engine.land_cover_fractions import assign
from vtkread import VtkFile, read_vtk

TOLERANCE = "3"
ARRAYS = ("land_cover_offsets", "land_cover_classes", "land_cover_fractions")
STRINGS = (
    "land_cover_codes",
    "land_cover_source",
    "land_cover_credit",
    "land_cover_licence_note",
    "land_cover_cite",
)
#: On a fraction: float32 rounding is 6e-8; both sides clip the same straight
#: index-space edges, so nothing else is in the difference.
FRACTION = 1e-6
#: F1: the file's fractions sum to one to float32 rounding.
ONE = 1e-6


@pytest.fixture
def classes(monkeypatch: pytest.MonkeyPatch) -> Any:
    """The class source in `SOURCES`; its URL is never opened (the cache serves)."""
    template = "https://example.invalid/" + PATH_TEMPLATE
    return install(sources, monkeypatch, per_year_entry(sources, template))


@pytest.fixture
def heights(monkeypatch: pytest.MonkeyPatch) -> Any:
    return catalogue_with(sources, monkeypatch, "https://example.invalid/p.tif")


@pytest.fixture
def offline(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)


@pytest.fixture
def cache(tmp_path: Path) -> Path:
    return write_class_cache(tmp_path / "cache")


@pytest.fixture
def dem(tmp_path: Path) -> Path:
    return write_dem(tmp_path / "dem.tif")


@pytest.fixture
def domain(tmp_path: Path) -> Path:
    return write_domain(tmp_path / "domain.geojson")


def mesh(
    tmp_path: Path, dem: Path, domain: Path, cache: Path, *extra: str, out: str = "x.vtk"
) -> tuple[int, str, Path]:
    target = tmp_path / out
    code, output = invoke(
        "mesh",
        *("--dem", str(dem), "--domain", str(domain), "--tolerance", TOLERANCE),
        *("--cache", str(cache), "--out", str(target)),
        *extra,
    )
    return code, output, target


@dataclass(frozen=True)
class Covered:
    """A `.vtk` written with land cover, and the table read off it."""

    vtk: VtkFile
    output: str
    offsets: np.ndarray
    classes: np.ndarray
    fractions: np.ndarray

    @property
    def triangles(self) -> np.ndarray:
        return np.array([np.asarray(p) for p in self.vtk.polygons], dtype=np.int64)

    @property
    def rows(self) -> list[Row]:
        return rows_of(self.offsets, self.classes, self.fractions)

    def codes(self) -> np.ndarray:
        return np.asarray(self.vtk.cell_array("land_cover_code").values)

    def string(self, name: str) -> str:
        return file_field(self.vtk, name)


def covered(tmp_path: Path, dem: Path, domain: Path, cache: Path, *extra: str) -> Covered:
    code, output, target = mesh(tmp_path, dem, domain, cache, "--land-cover", KEY, *extra)
    assert code == 0, output
    assert "No such option" not in output, output
    vtk = read_vtk(target.read_bytes())
    for name in ARRAYS:
        assert name in vtk.field_data, f"no {name}; the file has {sorted(vtk.field_data)}"
    offsets, classes, fractions = (np.asarray(vtk.field_data[n].values) for n in ARRAYS)
    return Covered(vtk, output, offsets, classes, fractions.astype(np.float32))


# ------------------------------------------------------------ the file


class TestTheFile:
    @pytest.fixture
    def run(
        self, classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> Covered:
        return covered(tmp_path, dem, domain, cache)

    def test_the_arrays_have_the_designs_types_and_one_row_per_triangle(self, run: Covered) -> None:
        offsets, classes, fractions = (run.vtk.field_data[n] for n in ARRAYS)
        assert (offsets.type_name, classes.type_name, fractions.type_name) == (
            "unsigned_int", "unsigned_char", "float",
        )  # fmt: skip
        assert len(run.offsets) == len(run.vtk.polygons) + 1
        assert run.offsets[0] == 0 and run.offsets[-1] == len(run.classes) == len(run.fractions)
        assert (np.diff(run.offsets.astype(np.int64)) > 0).all(), "a triangle without a class"

    def test_f1_each_triangles_fractions_sum_to_one_with_classes_ascending(
        self, run: Covered
    ) -> None:
        for t, row in enumerate(run.rows):
            assert abs(sum(row.values()) - 1.0) <= ONE, f"triangle {t}: {row}"
            assert list(row) == sorted(row), f"triangle {t}: {row}"
            assert min(row.values()) > 0.0, f"triangle {t}: {row}"

    def test_the_code_is_the_largest_class_and_0_on_every_line_cell(self, run: Covered) -> None:
        codes = run.codes()
        lines = len(run.vtk.lines)
        assert lines > 0, "the domain's outline is a constraint, so the file has lines"
        assert not codes[:lines].any()
        assert_array_equal(codes[lines:], [dominant(row) for row in run.rows])

    def test_the_fractions_are_the_oracles(self, run: Covered) -> None:
        want = oracle_fractions(run.vtk.points, run.triangles, class_grid(2025))
        assert same_rows(run.rows, want, FRACTION) == []

    def test_the_classes_are_the_rasters_and_no_data_is_among_them(self, run: Covered) -> None:
        seen = set(run.classes.tolist())
        assert seen <= set(CLASSES) | {0}
        assert 0 in seen, "the mesh reaches the block of class 0 (1.1: a class like any other)"
        assert sum(len(row) > 1 for row in run.rows) > 10, "mixed triangles are the point"

    def test_the_strings_name_the_source_year_credit_licence_and_citation(
        self, run: Covered
    ) -> None:
        # 1.6 spells both with the source's display name ("MapBiomas Collection
        # 11"); RemoteSource has none, so the catalogue key stands in: it is the
        # name `--land-cover` and `fetch` take, so a reader can act on it.
        assert run.string("land_cover_codes") == (
            f"{KEY} class number; 0 = no data, and every constraint line"
        )
        assert run.string("land_cover_source") == f"{KEY}, 2025"
        credit = run.string("land_cover_credit")
        assert "Test classes" in credit and FETCHED in credit and "example.invalid" in credit
        assert run.string("land_cover_licence_note") == LICENCE
        assert run.string("land_cover_cite") == CITE[0]


# ------------------------------------------------------------ the year


class TestTheYear:
    def test_a_year_is_honoured_and_the_last_year_is_the_default(
        self, classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        latest = covered(tmp_path / "a", dem, domain, cache)
        early = covered(tmp_path / "b", dem, domain, cache, "--land-cover-year", "1985")
        assert "1985" in early.string("land_cover_source")
        assert "2025" in latest.string("land_cover_source")
        want = oracle_fractions(early.vtk.points, early.triangles, class_grid(1985))
        assert same_rows(early.rows, want, FRACTION) == []
        assert early.rows != latest.rows, "the two years differ cell by cell"

    def test_a_year_out_of_range_is_a_usage_error_naming_the_range(
        self, classes: Any, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        code, output, target = mesh(
            tmp_path, dem, domain, cache, "--land-cover", KEY, "--land-cover-year", "1984"
        )
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "1984" in output and "1985" in output and "2025" in output, output
        assert not target.exists()

    def test_a_year_without_land_cover_is_a_usage_error(
        self, classes: Any, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        code, output, target = mesh(tmp_path, dem, domain, cache, "--land-cover-year", "2025")
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--land-cover" in output, output
        assert not target.exists()


# ------------------------------------------------------------ refusals


class TestRefusals:
    def test_a_heights_source_is_refused_naming_classes(
        self, heights: Any, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        code, output, target = mesh(tmp_path, dem, domain, cache, "--land-cover", DEM_KEY)
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert DEM_KEY in output and "classes" in output and "--land-cover" in output, output
        assert not target.exists()

    def test_a_year_not_in_the_cache_is_refused_naming_the_fetch(
        self, classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        code, output, target = mesh(
            tmp_path, dem, domain, cache, "--land-cover", KEY, "--land-cover-year", "2000"
        )
        assert code != 0, output
        assert "Traceback" not in output and "No such option" not in output, output
        assert KEY in output and "fetch" in output and "2000" in output, output
        assert not target.exists()

    def test_a_cache_without_the_source_is_refused_naming_the_fetch(
        self, classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        """A cache root that never saw `fetch KEY`: no footprints at all, the
        same refusal as a missing year."""
        empty = tmp_path / "empty"
        empty.mkdir()
        code, output, target = mesh(tmp_path, dem, domain, empty, "--land-cover", KEY)
        assert code == USAGE, output
        assert "Traceback" not in output and "No such option" not in output, output
        assert KEY in output and "fetch" in output and "2025" in output, output
        assert not target.exists()

    def test_a_coded_features_map_beside_it_is_refused(
        self, classes: Any, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        """Section 8: one code system per mesh (16e's rule)."""
        ring = domain_ring()
        inner = Polygon(ring).buffer(-200.0)
        path = write_geojson(
            tmp_path / "corine.geojson", [Feat(1, inner, {"Code_18": "311"})], crs=MESH_CRS
        )
        code, output, target = mesh(
            tmp_path, dem, domain, cache,
            *("--features", str(path), "--features-map", "corine", "--land-cover", KEY),
        )  # fmt: skip
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--land-cover" in output and "corine" in output, output
        assert not target.exists()

    def test_without_a_dem_it_is_refused(self, classes: Any, tmp_path: Path, cache: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh", "catchment", "--flat", "--land-cover", KEY, "--cache", str(cache),
            "--out", str(target),
        )  # fmt: skip
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--dem" in output and "--land-cover" in output, output
        assert not target.exists()

    def test_an_unknown_key_is_a_usage_error_naming_the_catalogue(
        self, classes: Any, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        code, output, target = mesh(tmp_path, dem, domain, cache, "--land-cover", "no-such-source")
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "no-such-source" in output and KEY in output and "--land-cover" in output, output
        assert not target.exists()

    def test_without_a_cache_root_it_is_refused_naming_the_option(
        self, classes: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        """No `--cache` and no `RASPUTIN_DATA`: there is nowhere to read the
        classes from, and the refusal says where they would come from."""
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh", "--dem", str(dem), "--domain", str(domain), "--tolerance", TOLERANCE,
            "--out", str(target), "--land-cover", KEY,
        )  # fmt: skip
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--cache" in output and "RASPUTIN_DATA" in output, output
        assert not target.exists()

    def test_the_help_names_both_options(self) -> None:
        text = plain(runner.invoke(app, ["mesh", "--help"]).output)
        assert "--land-cover-year" in text and "--land-cover " in text, text


# ------------------------------------------------------------ chunks


class TestChunks:
    def test_small_chunks_give_one_calls_table_bit_for_bit(
        self, classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        """1.2: the rows are concatenated in triangle order, so the table does
        not depend on how the triangles are cut into calls. The CLI's scene is
        far below one chunk, so the stitching of several (each with its own
        window and row areas) is exercised here through `assign` directly; the
        chunk size itself (2^18 in file order, the Hilbert order being 26c's)
        is not pinned."""
        run = covered(tmp_path, dem, domain, cache)
        repository = CacheRepository(cache, KEY, nodata=0.0)
        name = object_id(2025)
        whole = assign(
            run.vtk.points, run.triangles, mesh_crs=MESH_CRS, repository=repository, name=name
        )
        assert len(run.triangles) > 3 * 5, "several chunks, the last one partial"
        pieces = assign(
            run.vtk.points, run.triangles, mesh_crs=MESH_CRS, repository=repository, name=name,
            chunk=5,
        )  # fmt: skip
        for ours, theirs in (
            (whole.table.offsets, pieces.table.offsets),
            (whole.table.classes, pieces.table.classes),
            (whole.table.fractions, pieces.table.fractions),
            (whole.codes, pieces.codes),
        ):
            assert ours.dtype == theirs.dtype
            assert_array_equal(ours, theirs)
        # And one call is what the CLI wrote.
        assert_array_equal(whole.table.offsets, run.offsets)
        assert_array_equal(whole.table.classes, run.classes)
        assert_array_equal(whole.table.fractions, run.fractions)

    def test_a_mesh_past_the_raster_and_an_object_not_cached_are_refused(
        self, classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path, cache: Path
    ) -> None:
        """1.2: a triangle outside the window is a refusal, never a clamp; the
        Python side refuses before the C++ sees it. An object the cache lacks
        is refused by name (the CLI checks first; the API stands on its own)."""
        run = covered(tmp_path, dem, domain, cache)
        repository = CacheRepository(cache, KEY, nodata=0.0)
        east = np.asarray(run.vtk.points, dtype=np.float64) + np.array([5e5, 0.0, 0.0])
        with pytest.raises(ValueError, match="outside"):
            assign(
                east, run.triangles, mesh_crs=MESH_CRS, repository=repository, name=object_id(2025)
            )
        with pytest.raises(ValueError, match="not in the cache"):
            assign(
                run.vtk.points, run.triangles, mesh_crs=MESH_CRS, repository=repository,
                name=object_id(2000),
            )  # fmt: skip


# ------------------------------------------------------------ the control


def test_without_land_cover_the_file_carries_none_of_it(
    classes: Any, tmp_path: Path, dem: Path, domain: Path, cache: Path
) -> None:
    code, output, target = mesh(tmp_path, dem, domain, cache)
    assert code == 0, output
    vtk = read_vtk(target.read_bytes())
    assert not (set(ARRAYS) | set(STRINGS)) & set(vtk.field_data)
    with pytest.raises(KeyError):
        vtk.cell_array("land_cover_code")


# ------------------------------------------------------------ --stats and --record


def test_the_strings_are_stats_rows_and_record_entries(
    classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path, cache: Path
) -> None:
    report, record = tmp_path / "x.md", tmp_path / "x.json"
    run = covered(tmp_path, dem, domain, cache, "--stats", str(report), "--record", str(record))
    text = report.read_text(encoding="utf-8")
    recorded = json.loads(record.read_text(encoding="ascii"))
    for name in STRINGS:
        assert stats_row(text, name) == run.string(name), name
        assert recorded[name] == run.string(name), name


# ------------------------------------------------------------ .ply


def test_the_ply_face_lists_carry_the_same_table(
    classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path, cache: Path
) -> None:
    run = covered(tmp_path, dem, domain, cache)
    (tmp_path / "ply").mkdir()
    code, output, target = mesh(
        tmp_path / "ply", dem, domain, cache, "--land-cover", KEY, out="x.ply"
    )
    assert code == 0, output
    header, data = read_ply(target.read_bytes())
    names = [p.name for p in header.element("face").properties]
    assert names == [
        "vertex_indices",
        "land_cover_classes",
        "land_cover_fractions",
        "land_cover_code",
    ]
    face = data["face"]
    assert len(face["land_cover_classes"]) == len(run.vtk.polygons), "the same mesh twice"
    got = [
        {int(k): float(v) for k, v in zip(c.tolist(), f.tolist(), strict=True)}
        for c, f in zip(face["land_cover_classes"], face["land_cover_fractions"], strict=True)
    ]
    assert same_rows(got, run.rows, 0.0) == []
    assert_array_equal(np.asarray(face["land_cover_code"]), [dominant(row) for row in run.rows])
    fields = ply_fields(header.comments)
    for name in STRINGS:
        assert fields[name] == run.string(name), name


# ------------------------------------------------------------ VTK readback


def test_vtk_reads_the_field_arrays_back(
    classes: Any, offline: None, tmp_path: Path, dem: Path, domain: Path, cache: Path
) -> None:
    """`vtkPolyDataReader` with `ReadAllFieldsOn` hands back the three dataset
    arrays with their lengths, and the cell array 0 on the lines (which VTK
    numbers first). Skipped without `vtk`, as 16c's read-back is; CI runs it
    in the `viewer` extra's step."""
    vtk = pytest.importorskip("vtk")
    code, output, target = mesh(tmp_path, dem, domain, cache, "--land-cover", KEY)
    assert code == 0, output
    reader = vtk.vtkPolyDataReader()
    errors: list[str] = []
    reader.AddObserver("ErrorEvent", lambda _obj, _event: errors.append("error"))
    reader.SetFileName(str(target))
    reader.ReadAllFieldsOn()
    reader.Update()
    assert not errors
    polydata = reader.GetOutput()
    ours = read_vtk(target.read_bytes())
    field = polydata.GetFieldData()
    for name in ARRAYS:
        array = field.GetArray(name)
        assert array is not None, name
        values = np.array([array.GetValue(i) for i in range(array.GetNumberOfTuples())])
        assert_array_equal(values, np.asarray(ours.field_data[name].values))
    assert (
        field.GetArray("land_cover_offsets").GetNumberOfTuples() == polydata.GetNumberOfPolys() + 1
    )
    codes = polydata.GetCellData().GetArray("land_cover_code")
    assert codes is not None
    lines = polydata.GetNumberOfLines()
    assert lines > 0 and not any(codes.GetValue(i) for i in range(lines))
