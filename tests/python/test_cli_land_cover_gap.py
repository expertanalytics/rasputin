"""`rasputin mesh --land-cover ... [--land-cover-partial]`: a part of the domain
no land-cover or `--prefer` polygon covers stops the run; the override fills
it, marked, and records it.

`docs/increments/122-lakes-overrule.md`, section 3.6 (the gap, the refusal,
its words, the hint, the override), section 3.4 (the mark and the record),
section 6, and section 7, tests 10, 10b and 17 to 22. The scene is
`test_cli_land_cover_polygons.py`'s: its DEM (40 x 40 nodes at 60 m) and
domain (4,490,161 m²) in SIRGAS 2000 / UTM 23S (EPSG:31983), `--tolerance 3`.
`partial_cover()` is the scene's three CORINE polygons: they leave
1,905,280.4 m² of the domain uncovered, in two parts, 1.8751 km² around
507344 8315355 and the conifer's hole, 0.0301 km² (measured with shapely
2.2's `maximum_inscribed_circle` at 0.001 m on these fixtures, by
`test_the_gap_figures`). `cover_polygons()` adds heath 322 over the rest of
the DEM's box, so it leaves no gap.

THE ORACLE: the gap figures above, recomputed in `test_the_gap_figures` from
the fixture polygons by GEOS, never from the producer; codes and marks read
back from the `.gpkg`'s `triangles` layer and `run` table. A place is
compared to 1 m (the words give metres); an area to 0.5 m² (the design's
bound on `land_cover_gap_m2`, whose largest input is 1.9e6 m²).

PINNED HERE, where the design leaves it open (listed for `@architect`):

- One part is "(1 part)", more "(N parts)".
- `land_cover_gaps` is the number of parts and `land_cover_gap_m2` a number,
  each written as the `run` table's value; `land_cover_partial` is "on" or
  "off"; with no land cover none of the three is written.
- The `.gpkg`'s `land_cover_filled` is 0 or 1 on every triangle that has a
  `code`.
- The hint's ids in "covered by" are separated by ", ".
"""

from __future__ import annotations

import os
import re
import socket
import sqlite3
from contextlib import closing
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon, box

from class_fixtures import KEY
from cli_helpers import USAGE, help_text, invoke, squashed
from dataset_fixtures import write_sources
from feature_fixtures import Feat, write_geojson
from fraction_fixtures import (
    MESH_CRS,
    dem_origin,
    domain_ring,
    write_class_cache,
    write_dem,
    write_domain,
)
from test_cli_land_cover_polygons import cover_polygons, mesh, partial_cover, write_cover_geojson
from test_cli_mesh_cache import no_network
from test_cli_mesh_cells import corners, scene
from test_cli_mesh_cells import mesh as cells_mesh
from test_cli_mesh_edges import CLEAN_CLC, CLEAN_ID, DTM10, OWN_HOME
from test_cli_mesh_prefer import island_scene, write_lake
from vector_fixtures import vector_entry

__all__ = ["scene"]  # a fixture, used by name

#: On `land_cover_gap_m2`, m²: the design's bound.
GAP_M2 = 0.5
#: The scene's gap, by GEOS on the fixtures (module docstring).
GAP_TOTAL_M2 = 1_905_280.4
LARGEST_KM2, LARGEST_AT = "1.8751", (507344, 8315355)
DOMAIN_KM2 = "4.4902"
#: Section 3.6's refusal, whitespace removed (Rich may break lines anywhere).
REFUSED = (
    "ofthedomaininnopolygonofany--land-coveror--preferfile",
    "missinglandcoverisrefused",
    "Addthefilethatcoversit,orpass--land-cover-partialtofilleachgap"
    "fromthenearestnon-waterclassandmarkit",
)
#: Section 3.2's step-2 refusal.
NO_DONOR = (
    "thelandcovergivesnotriangleofthemeshaclass,socode0cannotbefilled;itdoesnotcoverthedomain"
)
PART = re.compile(r"(\d+\.\d{4}) km2 around (\d+) (\d+) in (EPSG:\d+)")


@pytest.fixture
def dem(tmp_path: Path) -> Path:
    return write_dem(tmp_path / "dem.tif")


@pytest.fixture
def domain(tmp_path: Path) -> Path:
    return write_domain(tmp_path / "domain.geojson")


@pytest.fixture
def gappy(tmp_path: Path) -> Path:
    return write_cover_geojson(tmp_path / "gappy.geojson", partial_cover())


def corine(cover: Path, *more: str) -> tuple[str, ...]:
    return ("--land-cover", str(cover), "--land-cover-classes", "corine", *more)


def record(conn: sqlite3.Connection, name: str) -> str | None:
    row = conn.execute("SELECT value FROM run WHERE name = ?", (name,)).fetchone()
    return None if row is None else str(row[0])


def codes_and_marks(conn: sqlite3.Connection) -> tuple[np.ndarray, np.ndarray]:
    rows = conn.execute("SELECT code, land_cover_filled FROM triangles ORDER BY fid").fetchall()
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


def parts_of(output: str) -> list[tuple[str, int, int, str]]:
    return [(a, int(x), int(y), crs) for a, x, y, crs in PART.findall(output)]


def test_the_gap_figures() -> None:
    """The red premise of 17 and 18, by GEOS: `partial_cover()` leaves two
    parts, 1,905,280.4 m², the larger 1.8751 km² around 507344 8315355."""
    rest = Polygon(domain_ring()).difference(shapely.union_all([g for _, g in partial_cover()]))
    parts = sorted(shapely.get_parts(rest), key=lambda p: -p.area)
    assert len(parts) == 2
    assert rest.area == pytest.approx(GAP_TOTAL_M2, abs=0.05)
    assert f"{parts[0].area / 1e6:.4f}" == LARGEST_KM2
    centre = shapely.get_coordinates(shapely.maximum_inscribed_circle(parts[0], 0.001))[0]
    assert np.abs(np.round(centre) - LARGEST_AT).max() == 0


def test_the_scenes_cover_leaves_no_gap() -> None:
    rest = Polygon(domain_ring()).difference(shapely.union_all([g for _, g in cover_polygons()]))
    assert rest.area == 0.0


def test_the_help_names_the_override() -> None:
    flat = squashed(help_text("mesh"))
    assert "--land-cover-partial" in flat
    assert "Mesheven" in flat and "wherenoland-coverfilecoversthedomain" in flat


# ------------------------------------------------------------ test 17, the refusal


def test_17_a_gap_is_refused_naming_its_parts(
    tmp_path: Path, dem: Path, domain: Path, gappy: Path
) -> None:
    target = tmp_path / "x.gpkg"
    code, output = mesh(dem, domain, target, *corine(gappy))
    assert code == USAGE, output
    assert "Invalid value for --land-cover" in output, output
    flat = squashed(output)
    assert "thelandcoverleaves1.9053km2" in flat, output
    assert "(2parts)" in flat, output
    for words in REFUSED:
        assert words in flat, (words, output)
    parts = parts_of(output)
    assert len(parts) == 2, output
    assert [a for a, *_ in parts] == [LARGEST_KM2, "0.0301"], output
    assert {crs for *_, crs in parts} == {"EPSG:31983"}, output
    assert abs(parts[0][1] - LARGEST_AT[0]) <= 1 and abs(parts[0][2] - LARGEST_AT[1]) <= 1
    assert "covered by" not in output, output
    assert not target.exists()


def test_17_the_refusal_comes_before_the_mesh(
    tmp_path: Path, dem: Path, domain: Path, gappy: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Section 3.6: right after the map is built, before the lines and the
    mesh: no phase after `land cover read` is timed."""
    from tin_engine.stats import PhaseClock

    seen: list[str] = []
    original = PhaseClock.add

    def spy(self: PhaseClock, name: str, seconds: float) -> None:
        seen.append(name)
        original(self, name, seconds)

    monkeypatch.setattr(PhaseClock, "add", spy)
    code, output = mesh(dem, domain, tmp_path / "x.gpkg", *corine(gappy))
    assert code == USAGE, output
    assert "land cover read" in seen, seen
    assert seen[-1] == "land cover read", seen


def test_17_the_pieces_refuse_before_any_piece(
    tmp_path: Path, dem: Path, domain: Path, gappy: Path
) -> None:
    """Section 3.6: on the whole domain once, so `--pieces` writes none."""
    target = tmp_path / "out" / "x.vtk"
    target.parent.mkdir()
    code, output = mesh(dem, domain, target, *corine(gappy), "--edges", "all", "--pieces", "4")
    assert code == USAGE, output
    assert "missinglandcoverisrefused" in squashed(output), output
    assert list(target.parent.iterdir()) == []


# ------------------------------------------------------------ test 10, nothing covers


def outside_cover(tmp_path: Path) -> Path:
    x0, y0 = dem_origin()
    far = box(x0 + 9000.3, y0 + 9000.1, x0 + 9500.7, y0 + 9500.9)
    return write_cover_geojson(tmp_path / "far.geojson", [(311, far)])


def test_10_a_land_cover_outside_the_domain_is_one_part_the_domain(
    tmp_path: Path, dem: Path, domain: Path
) -> None:
    target = tmp_path / "x.gpkg"
    code, output = mesh(dem, domain, target, *corine(outside_cover(tmp_path)))
    assert code == USAGE, output
    assert "Invalid value for --land-cover" in output, output
    flat = squashed(output)
    assert f"thelandcoverleaves{DOMAIN_KM2}km2" in flat, output
    assert "(1part)" in flat, output
    for words in REFUSED:
        assert words in flat, (words, output)
    assert [a for a, *_ in parts_of(output)] == [DOMAIN_KM2], output
    assert not target.exists()


def test_10b_the_same_with_the_override_has_no_class_to_fill_from(
    tmp_path: Path, dem: Path, domain: Path
) -> None:
    target = tmp_path / "x.gpkg"
    code, output = mesh(
        dem, domain, target, *corine(outside_cover(tmp_path), "--land-cover-partial")
    )
    assert code == USAGE, output
    assert "Invalid value for --land-cover" in output, output
    assert NO_DONOR in squashed(output), output
    assert not target.exists()


# ------------------------------------------------------------ test 18, the override


def test_18_the_override_fills_marks_and_records_the_gap(
    tmp_path: Path, dem: Path, domain: Path, gappy: Path
) -> None:
    target = tmp_path / "x.gpkg"
    code, output = mesh(dem, domain, target, *corine(gappy, "--land-cover-partial"))
    assert code == 0, output
    line = "landcover:1.9053km2ofthedomaininnofile(2parts,thelargest1.8751km2around"
    flat = squashed(output)
    assert flat.count(line) == 1, output
    assert flat.count("filledby--land-cover-partial") == 1, output
    with closing(sqlite3.connect(target)) as conn:
        codes, marks = codes_and_marks(conn)
        assert len(codes) and not (codes == 0).any(), "a triangle kept code 0"
        assert set(marks.tolist()) == {0, 1}, set(marks.tolist())
        assert record(conn, "land_cover_partial") == "on"
        assert record(conn, "land_cover_gaps") == "2"
        gap = record(conn, "land_cover_gap_m2")
        assert gap is not None and float(gap) == pytest.approx(GAP_TOTAL_M2, abs=GAP_M2)
        zero = conn.execute("SELECT count(*) FROM land_cover_fractions WHERE code = 0")
        assert zero.fetchone() == (0,)
        filled = record(conn, "land_cover_filled")
        assert filled is not None and int(filled) == int(marks.sum())


def test_18_the_marked_triangles_are_those_reaching_the_gap(
    tmp_path: Path, dem: Path, domain: Path, gappy: Path
) -> None:
    """Section 3.3: a triangle whose 0 share moved is marked; so a triangle is
    marked exactly where it reaches the gap (share over 1e-6), and not where
    it does not (under 1e-9); none lies between."""
    target = tmp_path / "x.gpkg"
    code, output = mesh(dem, domain, target, *corine(gappy, "--land-cover-partial"))
    assert code == 0, output
    covered = shapely.union_all([g for _, g in partial_cover()])
    with closing(sqlite3.connect(target)) as conn:
        _, marks = codes_and_marks(conn)
        triangles = shapely.polygons(corners(conn)[:, :, :2])
    rest = 1.0 - shapely.area(shapely.intersection(triangles, covered)) / shapely.area(triangles)
    assert not ((rest > 1e-9) & (rest <= 1e-6)).any()
    assert ((marks == 1) == (rest > 1e-6)).all()


def test_18_without_land_cover_the_three_keys_are_absent(
    tmp_path: Path, dem: Path, domain: Path
) -> None:
    target = tmp_path / "x.gpkg"
    code, output = mesh(dem, domain, target)
    assert code == 0, output
    with closing(sqlite3.connect(target)) as conn:
        for name in ("land_cover_partial", "land_cover_gaps", "land_cover_gap_m2"):
            assert record(conn, name) is None, name


def test_18_with_full_cover_no_gap_and_the_override_off(
    tmp_path: Path, dem: Path, domain: Path
) -> None:
    cover = write_cover_geojson(tmp_path / "cover.geojson", cover_polygons())
    target = tmp_path / "x.gpkg"
    code, output = mesh(dem, domain, target, *corine(cover))
    assert code == 0, output
    assert "in no file" not in output, output
    with closing(sqlite3.connect(target)) as conn:
        assert record(conn, "land_cover_partial") == "off"
        assert record(conn, "land_cover_gaps") == "0"
        gap = record(conn, "land_cover_gap_m2")
        assert gap is not None and float(gap) == 0.0
        assert record(conn, "land_cover_filled") == "0"
        _, marks = codes_and_marks(conn)
        assert not marks.any()


# ------------------------------------------------------------ test 19, an island


def test_19_an_island_fills_without_the_override(tmp_path: Path, dem: Path, domain: Path) -> None:
    """`--prefer lake=` L less H (39,999.84 m²): the base's water in H borders
    no kept polygon, so it is left out of the map; it is no gap (the base's
    lake covers it), so the run goes on and the fill gives it heath, marked.
    `--edges all` keeps H's border a constraint, so a triangle lies in H or
    out of it."""
    lake, hole, cover = island_scene()
    assert hole.area == pytest.approx(39_999.84, abs=1e-3)
    land = write_cover_geojson(tmp_path / "cover.geojson", cover)
    preferred = write_lake(tmp_path / "lake.geojson", lake.difference(hole))
    target = tmp_path / "x.gpkg"
    args = corine(land, "--prefer", f"lake={preferred}", "--edges", "all")
    code, output = mesh(dem, domain, target, *args)
    assert code == 0, output
    with closing(sqlite3.connect(target)) as conn:
        codes, marks = codes_and_marks(conn)
        centroids = shapely.points(corners(conn)[:, :, :2].mean(axis=1))
        gap = record(conn, "land_cover_gap_m2")
        assert gap is not None and float(gap) == 0.0
        assert record(conn, "land_cover_partial") == "off"
    inside = shapely.within(centroids, hole)
    assert inside.sum() > 0, "no triangle in H"
    assert set(codes[inside].tolist()) == {322}
    assert (marks[inside] == 1).all()
    assert (marks[~inside] == 0).all()
    assert not (codes == 0).any()


# ------------------------------------------------------------ test 20, the hint


def test_20_the_hint_names_the_clean_catalogue_file_covering_each_part(
    tmp_path: Path, dem: Path, domain: Path, gappy: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A user catalogue in `$RASPUTIN_DATA` with corine entries `east` and
    `east-clean-31983` and an ar50 entry `other`, each the scene's gap as one
    211 polygon: each part names `east-clean-31983` (of a raw id and its
    clean sibling, only the sibling), never `east` or `other` (another class
    file)."""
    rest = Polygon(domain_ring()).difference(shapely.union_all([g for _, g in partial_cover()]))
    gap = write_geojson(tmp_path / "gap.geojson", [Feat(1, rest, {"Code_18": "211"})], crs=MESH_CRS)
    data = tmp_path / "data"
    monkeypatch.setenv("RASPUTIN_DATA", str(data))
    write_sources(data / "sources.toml", {
        "east": vector_entry(gap, "East", classes="corine"),
        "east-clean-31983": vector_entry(gap, "East, valid in EPSG:31983", classes="corine"),
        "other": vector_entry(gap, "Other", classes="ar50"),
    })  # fmt: skip
    target = tmp_path / "x.gpkg"
    code, output = mesh(dem, domain, target, *corine(gappy))
    assert code == USAGE, output
    assert len(parts_of(output)) == 2, output
    named = re.findall(r"covered by ([\w\-, ]+?)(?=\s+\d+\.\d{4} km2|\s*$)", output)
    assert named == ["east-clean-31983", "east-clean-31983"], output
    assert not re.search(r"\bother\b", output), output
    assert not re.search(r"\beast\b(?!-clean)", output), output


# ------------------------------------------------------------ test 21, a raster


@pytest.fixture
def blocked(scene: dict[str, Path], tmp_path: Path) -> dict[str, Path]:
    """The CLI scene, its cache written with the block of class 0 (no data)
    that the mesh reaches: `write_class_cache(root, no_data=True)`."""
    return {**scene, "cache": write_class_cache(tmp_path / "blocked", no_data=True)}


def test_21_a_raster_reaching_no_data_is_refused(blocked: dict[str, Path], tmp_path: Path) -> None:
    """The cached class raster with its no-data block, which the mesh
    reaches, and no `--land-cover-partial`."""
    target = tmp_path / "x.gpkg"
    code, output = cells_mesh(
        blocked, target, "--land-cover", KEY, "--land-cover-classes", "mapbiomas-c11"
    )
    assert code == USAGE, output
    assert "Invalid value for --land-cover" in output, output
    flat = squashed(output)
    for words in REFUSED[:2]:
        assert words in flat, (words, output)
    assert not target.exists()


def test_21_the_override_merges_every_zero_share(blocked: dict[str, Path], tmp_path: Path) -> None:
    target = tmp_path / "x.gpkg"
    args = ("--land-cover", KEY, "--land-cover-classes", "mapbiomas-c11", "--land-cover-partial")
    code, output = cells_mesh(blocked, target, *args)
    assert code == 0, output
    with closing(sqlite3.connect(target)) as conn:
        zero = conn.execute("SELECT count(*) FROM land_cover_fractions WHERE code = 0")
        assert zero.fetchone() == (0,)
        codes, marks = codes_and_marks(conn)
        assert not (codes == 0).any()
        assert marks.any()
        assert record(conn, "land_cover_partial") == "on"


# ------------------------------------------------------------ test 22, an edge mode


def test_22_a_class_left_out_of_edges_is_not_filled_from_water(
    tmp_path: Path, dem: Path, domain: Path
) -> None:
    """Section 3.7: `--edges water` over a full cover (test 19's heath 322
    and lake L, no `--prefer`, no gap). The borders' polygons are the
    water's alone, so the land is code 0 on today's labelling, though a
    polygon of the map covers it; the fill's donors would then be water only
    and give it 512. Code 0 means "in no polygon of the map" under every
    `--edges`: the land takes 322 from the map (`label_triangles(...,
    rest=)`), unmarked, the lake 512, and the record's `land_cover_filled`
    is 0."""
    lake, _, cover = island_scene()
    land = write_cover_geojson(tmp_path / "cover.geojson", cover)
    target = tmp_path / "x.gpkg"
    code, output = mesh(dem, domain, target, *corine(land, "--edges", "water"))
    assert code == 0, output
    with closing(sqlite3.connect(target)) as conn:
        codes, marks = codes_and_marks(conn)
        centroids = shapely.points(corners(conn)[:, :, :2].mean(axis=1))
        filled = record(conn, "land_cover_filled")
    wet = shapely.within(centroids, lake)
    assert wet.any() and (~wet).any()
    assert set(codes[wet].tolist()) == {512}
    assert set(codes[~wet].tolist()) == {322}
    assert not marks.any()
    assert filled == "0"


# ------------------------------------------------------------ test 11, Bygdin (Ola's data)

BYGDIN = (
    Path(__file__).resolve().parents[2]
    / "docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson"
)
NVE_LAKES = OWN_HOME / "projects/rasputin_scratch/temp-check-20261021/nve-lakes-25833.geojson"
needs_bygdin = pytest.mark.skipif(
    not (DTM10.is_dir() and CLEAN_CLC.exists() and NVE_LAKES.exists()),
    reason="Ola's DTM10, clean clc2018 and NVE lakes are not here",
)
#: `edges_max_shift_m`'s default: the border clean-up may move a shore so far.
SHIFT = 50.0


@needs_bygdin
def test_11_bygdin_with_nve_lakes_no_code_0_and_the_islands_marked(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Local only: Bygdin with `clc2018-clean-25833` and `--prefer lake=` NVE
    (section 2's run): no triangle code 0; 145 marked, 0.1078 km² (±0.0001);
    every triangle whose centroid lies in an NVE lake farther than 50 m from
    its shore carries 512."""
    import json

    from tin_engine.io.vtk_read import read_mesh

    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
    write_sources(config, {
        CLEAN_ID: vector_entry(CLEAN_CLC, "clc2018, valid in EPSG:25833", classes="corine",
                               layer="features"),
    })  # fmt: skip
    target, rec = tmp_path / "bygdin.vtk", tmp_path / "bygdin.json"
    code, output = invoke(
        "mesh", "--dem", str(DTM10), "--domain", str(BYGDIN), "--tolerance", "10",
        "--edges", "all", "--threads", "8", "--land-cover", CLEAN_ID,
        "--prefer", f"lake={NVE_LAKES}", "--out", str(target), "--record", str(rec),
    )  # fmt: skip
    assert code == 0, output
    record_json = json.loads(rec.read_text(encoding="ascii"))
    assert record_json["land_cover_filled"] == 145
    assert record_json["land_cover_filled_m2"] / 1e6 == pytest.approx(0.1078, abs=1e-4)
    assert record_json["land_cover_gaps"] == 0
    mesh_read = read_mesh(target.read_bytes())
    codes = mesh_read.land_cover_code
    assert codes is not None and not (codes == 0).any()
    assert mesh_read.land_cover_filled is not None
    assert int(np.asarray(mesh_read.land_cover_filled).sum()) == 145
    lakes = shapely.union_all(
        [
            shapely.geometry.shape(f["geometry"])
            for f in json.loads(NVE_LAKES.read_text())["features"]
        ]
    )
    xy = mesh_read.vertices[:, :2][mesh_read.triangles].mean(axis=1)
    centroids = shapely.points(xy)
    deep = shapely.within(centroids, lakes) & (
        shapely.distance(centroids, shapely.boundary(lakes)) > SHIFT
    )
    assert deep.sum() > 1000, deep.sum()
    assert set(codes[deep].tolist()) == {512}
