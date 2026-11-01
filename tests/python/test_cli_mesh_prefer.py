"""`rasputin mesh --land-cover ... --prefer lake=<source>`: the lake class
taken from a preferred polygon source, the rest from the land cover.

`docs/increments/107-prefer-classes.md`, section 2 (the rule, the warning,
the refusals), section 4 (legend `[prefer]` table, CLI, credit) and section 5,
cases 6 and 7. The scene is `test_cli_mesh_cells.py`'s on the polygon route
with `corine` (`cover_polygons`: forest 311 and 312, water 512, the rest
uncovered). The preferred lake rectangle covers the east end of forest 311,
the uncovered land between and the west of water 512. What is left of the
water borders no kept polygon (a 200 m gap on every side), so it is left out
(class 0): the lake share is the preferred lake's alone, and the warning fires
(the preferred lake is smaller than the water it replaces). Since
`docs/increments/122-lakes-overrule.md` the scene's 1.9 km² gap is refused,
so every run here passes `--land-cover-partial` (`corine`), and the
left-out water and the gap are filled from the nearest non-water class,
forest (test 9).

THE ORACLE: per triangle of the `triangles` layer read back, GEOS's
`intersection(...).area` share with the preferred lake (lake), and with the
forest polygons less the lake (forest). Bound `SHARE` = 1e-6, as in
`test_cli_mesh_reservoirs.py`: the cut at these coordinates (below 1e7 m,
triangles of 100 m and up) loses under 1e-9 of a triangle, and the shares are
float32 (6e-8).

PINNED HERE, where the design leaves it open (listed for `@architect`):

- Every refusal is a usage error (exit 2) before any phase is timed, naming
  `--prefer`; the raster and no-land-cover ones also name `--land-cover`, the
  unknown key names the key, the legend refusal names the legend.
- A file's warning names the land cover and the preferred source by file name
  (`cover.geojson`, `lake.geojson`), in the design's sentence.
- The run record holds a row whose name contains `prefer`, whose value names
  the key and carries the three areas (at least three numbers); the credit row
  is `land_cover_prefer_credit`, naming the key and the file.
- `--prefer` with a `.vtk` output is accepted.
- The map's borders (`docs/increments/110-prefer-lines.md`, tests 3 and 5,
  as `docs/increments/111-one-land-cover-map.md` routes them: `--edges all`,
  and a `--lines` river): the lines are the `lines` layer's rows carrying
  the edge name, unioned; with the lines from another source, `--prefer`
  changes no line (exactly equal geometry).
"""

from __future__ import annotations

import re
import sqlite3
from contextlib import closing
from pathlib import Path

import numpy as np
import pytest
import shapely
from numpy.testing import assert_allclose
from shapely.geometry import Polygon, box
from shapely.geometry.base import BaseGeometry

from class_fixtures import KEY
from cli_helpers import USAGE, squashed
from feature_fixtures import Feat, classes, write_geojson
from fraction_fixtures import MESH_CRS, dem_origin, domain_ring
from gpkg_fixtures import Layer, write_gpkg
from gpkg_fixtures import Row as GpkgRow
from legend_fixtures import coded_map, write
from test_cli_land_cover_polygons import cover_polygons, write_cover_geojson
from test_cli_mesh_cells import check_sums, column, corners, mesh, scene, written
from test_cli_mesh_reservoirs import (
    RESERVOIRS_KEY,
    install_reservoirs,
    layer_entry,
    oracle_r,
    phases,
    record,
    reservoir_rectangle,
    write_reservoirs,
)

__all__ = ["phases", "scene"]  # fixtures, used by name

#: On a share: see the module docstring.
SHARE = 1e-6
#: The warning's words, whitespace removed (Rich may break lines anywhere).
TOOK = "tookaneighbouringclass"


def preferred_lake() -> Polygon:
    """Over the east end of forest 311, the land beyond it and the west of
    water 512; no side on a DEM node (the DEM's step is 60 m)."""
    x0, y0 = dem_origin()
    return box(x0 + 1000.2, y0 - 1400.3, x0 + 1700.7, y0 - 800.6)


def big_lake() -> Polygon:
    """Over all of water 512 inside the domain, and more."""
    x0, y0 = dem_origin()
    return box(x0 + 1300.1, y0 - 1800.2, x0 + 2700.3, y0 - 500.4)


def write_lake(path: Path, lake: Polygon) -> Path:
    return write_geojson(path, [Feat(1, lake, {"vatnLnr": "1"})], crs=MESH_CRS)


def share(conn: sqlite3.Connection, shape: shapely.Geometry) -> np.ndarray:
    triangles = shapely.polygons(corners(conn)[:, :, :2])
    return shapely.area(shapely.intersection(triangles, shape)) / shapely.area(triangles)


@pytest.fixture
def cover(tmp_path: Path) -> Path:
    return write_cover_geojson(tmp_path / "cover.geojson", cover_polygons())


@pytest.fixture
def lake(tmp_path: Path) -> Path:
    return write_lake(tmp_path / "lake.geojson", preferred_lake())


def corine(cover: Path, *more: str) -> tuple[str, ...]:
    return (
        "--land-cover-partial",
        "--land-cover",
        str(cover),
        "--land-cover-classes",
        "corine",
        *more,
    )


# ------------------------------------------------------------ the legends


def test_the_legends_name_the_class_each_key_replaces() -> None:
    assert dict(classes("corine").prefer) == {"lake": 512, "glacier": 335}
    assert dict(classes("ar50").prefer) == {"lake": 81, "glacier": 70}


# ------------------------------------------------------------ case 7, end to end


def test_case7_the_lake_share_comes_from_the_preferred_polygons(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    args = corine(cover, "--prefer", f"lake={lake}")
    code, output = mesh(scene, tmp_path / "x.gpkg", *args)
    assert code == 0, output
    with closing(sqlite3.connect(tmp_path / "x.gpkg")) as conn:
        assert_allclose(column(conn, "cells", "lake"), share(conn, preferred_lake()), atol=SHARE)
        # 122, section 3.3: the left-out water and the scene's gap are 0
        # shares that join a non-water class; the scene's only ones are
        # forest, so the forest share is the forest less the lake and the rest.
        forest = shapely.union_all([g for c, g in cover_polygons() if c in (311, 312)])
        expected = share(conn, forest.difference(preferred_lake()))
        got = column(conn, "cells", "forest")
        assert (got >= expected - SHARE).all()
        assert_allclose(got, 1.0 - share(conn, preferred_lake()), atol=SHARE)
        assert not column(conn, "cells", "glacier").any()
        assert_allclose(column(conn, "cells", "unspecified"), 0.0, atol=SHARE)
        check_sums(conn)
        rows = conn.execute("SELECT name, value FROM run WHERE name LIKE '%prefer%'").fetchall()
        credit = record(conn, "land_cover_prefer_credit")
    assert credit is not None and "lake" in credit and "lake.geojson" in credit, credit
    areas = [str(v) for n, v in rows if "credit" not in n and "licence" not in n]
    assert any("lake" in v and len(re.findall(r"\d+(?:\.\d+)?", v)) >= 3 for v in areas), rows
    # less lake than the water it replaced: the warning, in the design's words
    flat = squashed(output)
    assert "lake:lake.geojsongives" in flat and "cover.geojsongave" in flat, output
    assert "ofcover.geojsonwaterhadnolake.geojsonlakeand" + TOOK in flat, output


def test_case7_reservoirs_split_the_preferred_lake_share(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    res = write_reservoirs(tmp_path / "res.geojson", wgs84=False)
    prefer = ("--prefer", f"lake={lake}")
    split = (*prefer, "--reservoirs", str(res))
    with (
        closing(written(scene, tmp_path / "a.gpkg", *corine(cover, *prefer))) as a,
        closing(written(scene, tmp_path / "b.gpkg", *corine(cover, *split))) as b,
    ):
        before = column(a, "cells", "lake")
        lake_b, reservoir = column(b, "cells", "lake"), column(b, "cells", "reservoir")
        r = oracle_r(b, reservoir_rectangle())
        assert_allclose(reservoir, np.minimum(before, r), rtol=0, atol=SHARE)
        assert_allclose(lake_b + reservoir, before, rtol=0, atol=1e-12)
        assert reservoir.max() > 0.5, reservoir.max()
        assert ((r > 0.5) & (before == 0)).any(), "no land under the rectangle"
        check_sums(b)


def test_case7_no_warning_when_the_preferred_lake_gives_more(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    big = write_lake(tmp_path / "big.geojson", big_lake())
    code, output = mesh(scene, tmp_path / "x.gpkg", *corine(cover, "--prefer", f"lake={big}"))
    assert code == 0, output
    assert TOOK not in squashed(output), output


def test_prefer_needs_no_gpkg(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    target = tmp_path / "x.vtk"
    code, output = mesh(scene, target, *corine(cover, "--prefer", f"lake={lake}"))
    assert code == 0, output
    assert target.exists()


# ------------------------------------------------------------ case 6, refusals


def refused(scene: dict[str, Path], tmp_path: Path, phases: list[str], *args: str) -> str:
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, *args)
    assert code == USAGE, output
    assert phases == [], phases
    assert not target.exists()
    flat = squashed(output)
    assert "--prefer" in flat, output
    return flat


def test_case6_an_unknown_key_is_refused_naming_it(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path, phases: list[str]
) -> None:
    flat = refused(scene, tmp_path, phases, *corine(cover, "--prefer", f"river={lake}"))
    assert "river" in flat, flat


def test_case6_a_legend_without_prefer_is_refused_naming_it(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path, phases: list[str]
) -> None:
    classes = (
        '[[class]]\ncode = 311\nlabel = "Forest"\ncolour = "#1f8d49"\ncategory = "forest"\n\n'
        '[[class]]\ncode = 512\nlabel = "Lake"\ncolour = "#2532e4"\nwater = true\n'
        'category = "lake"\n'
    )
    legend = write(tmp_path / "mine.toml", coded_map("mine", attribute="Code_18", classes=classes))
    args = ("--land-cover-partial", "--land-cover", str(cover), "--land-cover-classes", str(legend))
    flat = refused(scene, tmp_path, phases, *args, "--prefer", f"lake={lake}")
    assert "mine" in flat, flat


def test_case6_a_raster_land_cover_is_refused(
    scene: dict[str, Path], tmp_path: Path, lake: Path, phases: list[str]
) -> None:
    args = ("--land-cover-partial", "--land-cover", KEY, "--land-cover-classes", "mapbiomas-c11")
    flat = refused(scene, tmp_path, phases, *args, "--prefer", f"lake={lake}")
    assert "--land-cover" in flat, flat


def test_case6_prefer_alone_is_refused(
    scene: dict[str, Path], tmp_path: Path, lake: Path, phases: list[str]
) -> None:
    flat = refused(scene, tmp_path, phases, "--prefer", f"lake={lake}")
    assert "--land-cover" in flat, flat


def test_case6_a_key_with_no_source_is_refused(
    scene: dict[str, Path], tmp_path: Path, cover: Path, phases: list[str]
) -> None:
    refused(scene, tmp_path, phases, *corine(cover, "--prefer", "lake="))


# ------------------------------------------------------------ the MCP tool


def test_the_mcp_mesh_tool_takes_prefer_as_a_list() -> None:
    pytest.importorskip("mcp")
    from tin_engine.cli import app
    from tin_engine.mcp import schema as mcp_schema

    mesh_tool = {d.tool_name: d for d in mcp_schema.derive(app)}["mesh"]
    prefer = mesh_tool.tool().input_schema["properties"]["prefer"]
    assert "array" in str(prefer), prefer
    arguments = {"out": "x.vtk", "dem": ["d.tif"], "prefer": ["lake=a", "glacier=b"]}
    argv = mesh_tool.argv(mesh_tool.model.model_validate(arguments))
    joined = " ".join(argv)
    assert "--prefer lake=a" in joined and "--prefer glacier=b" in joined, argv


# ------------------------------------------------------------ the map's borders

#: The clean-up's switches off, as the design measured.
AS_READ = ("--edges-max-shift", "0", "--edges-repair", "0", "--edges-outline-snap", "0")
#: On a length in m: the `lines` layer's geometry differs from the polygons'
#: by under 1 mm at these coordinates (checked on the two runs below, 0.6 mm),
#: far below the design's 0.1 m figures.
LENGTH = 0.01


def edge_lines(conn: sqlite3.Connection, name: str) -> BaseGeometry:
    """The union of the `lines` layer's rows that carry the edge name `name`."""
    found = []
    for (blob,) in conn.execute(f'SELECT geom FROM lines WHERE "{name}" = 1'):
        envelope = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[(blob[3] >> 1) & 7]
        found.append(shapely.from_wkb(bytes(blob[8 + envelope :])))
    return shapely.union_all(found)


def test_the_maps_borders_take_the_preferred_lake(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    """The lake's border inside the domain is covered, none of the replaced
    water's border is kept, and the lines are the preferred map's borders:
    9,838.4 m without --prefer, 9,082.7 m with it."""
    common = (*corine(cover), "--edges", "all", "--shares", *AS_READ)
    with (
        closing(written(scene, tmp_path / "a.gpkg", *common)) as a,
        closing(written(scene, tmp_path / "b.gpkg", *common, "--prefer", f"lake={lake}")) as b,
    ):
        before, after = edge_lines(a, "land_cover"), edge_lines(b, "land_cover")
        domain = Polygon(domain_ring())
        polygons = dict(cover_polygons())
        shores = [polygons[311].difference(preferred_lake()), polygons[312], preferred_lake()]
        expected = shapely.union_all([g.boundary for g in shores]).intersection(domain)
        assert before.length == pytest.approx(9838.4, abs=0.05)
        assert expected.length == pytest.approx(9082.7, abs=0.05)  # the premise
        assert after.length == pytest.approx(expected.length, abs=LENGTH)
        shore = preferred_lake().boundary.intersection(domain)
        assert shore.length == pytest.approx(2600.4, abs=0.05)
        covered = shore.intersection(after.buffer(LENGTH)).length
        assert covered == pytest.approx(shore.length, abs=LENGTH), covered
        water = polygons[512].boundary.intersection(domain)
        assert water.length == pytest.approx(2756.4, abs=0.05)
        assert water.intersection(after.buffer(LENGTH)).length <= 10 * LENGTH
        gone = before.difference(after.buffer(1.0)).length
        assert gone == pytest.approx(3350.1, abs=0.05), gone
        # the cell shares stay case 7's, the 0 shares joining forest (122)
        assert_allclose(column(b, "cells", "lake"), share(b, preferred_lake()), atol=SHARE)
        rest = 1.0 - share(b, preferred_lake())
        assert_allclose(column(b, "cells", "forest"), rest, atol=SHARE)
        check_sums(b)


def test_lines_of_another_source_are_unchanged_by_prefer(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    """The preferred lake's ring as a `--lines` river: --prefer changes no line."""
    ring = write_geojson(
        tmp_path / "shore.geojson",
        [Feat(1, shapely.LineString(preferred_lake().exterior.coords), {})],
        crs=MESH_CRS,
    )
    common = (*corine(cover), "--lines", f"river={ring}", "--edges", "river", *AS_READ)
    with (
        closing(written(scene, tmp_path / "a.gpkg", *common)) as a,
        closing(written(scene, tmp_path / "b.gpkg", *common, "--prefer", f"lake={lake}")) as b,
    ):
        plain, preferred = edge_lines(a, "river"), edge_lines(b, "river")
    assert plain.length > 0  # the premise: the lines are there
    assert shapely.equals_exact(shapely.normalize(plain), shapely.normalize(preferred), 0.0)


# ------------------------------------------------------------ a spike under --prefer

#: A spike's sides are this far apart at its base, m: GEOS cannot rebuild
#: it (`docs/increments/111-one-land-cover-map.md`, section 3, step 6).
THIN = 1e-8


def spiked_lake(base: float) -> Polygon:
    """A 150 m by 200 m lake whose top edge, at ``base`` below the DEM's
    origin, runs a spike `THIN` wide at its base up to (991.3, -1987.2) from
    it: invalid, 30,000 m²."""
    x0, y0 = dem_origin()
    b = y0 - base
    return Polygon(
        [
            (x0 + 850, b),
            (x0 + 1000, b + THIN),
            (x0 + 991.3, y0 - 1987.2),
            (x0 + 1000, b),
            (x0 + 1000, b - 200),
            (x0 + 850, b - 200),
        ]
    )


#: Section 5, test 14's lake: its spike runs from 100 m below the box
#: `read_polygons` checks in to 341.9 m into it.
SPIKED_LAKE_BASE = 2429.1


def test_14_the_spiked_lakes_premise() -> None:
    lake = spiked_lake(SPIKED_LAKE_BASE)
    assert not lake.is_valid
    assert lake.area == pytest.approx(30_000.0, abs=1e-3)


def test_14_read_polygons_refuses_a_thin_spike_naming_the_feature(tmp_path: Path) -> None:
    """Test 14, as `docs/increments/119-core-assumes-clean.md` has it: the
    lake is invalid and reaches the box, so the reader refuses it by fid
    with GEOS's reason."""
    from tin_engine.feature_input import InvalidPolygons, read_polygons

    lake = spiked_lake(SPIKED_LAKE_BASE)
    path = write_lake(tmp_path / "lake.geojson", lake)
    with pytest.raises(InvalidPolygons) as caught:
        read_polygons(path, None, MESH_CRS, np.array(domain_ring()), lines="lines go to --lines")
    message = str(caught.value)
    assert message.startswith("lake.geojson: 1 of the 1 polygons"), message
    assert f"feature 1: {shapely.is_valid_reason(lake)}" in message, message


def test_14_prefer_refuses_a_thin_spike_with_no_traceback(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    lake = spiked_lake(SPIKED_LAKE_BASE)
    path = write_lake(tmp_path / "lake.geojson", lake)
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, *corine(cover, "--prefer", f"lake={path}"))
    assert code == USAGE, output
    flat = squashed(output)
    assert "Invalidvaluefor--prefer" in flat, output
    assert squashed(f"feature 1: {shapely.is_valid_reason(lake)}") in flat, output
    assert "Traceback" not in output, output
    assert not target.exists()


# ------------------------------------------------------------ a source with no polygon here


def far_lake() -> Polygon:
    """A lake 20 km south of the scene: no polygon meets the read region."""
    x0, y0 = dem_origin()
    return box(x0 + 1000.2, y0 - 20800.6, x0 + 1700.7, y0 - 20000.3)


def test_a_preferred_file_with_no_polygon_here_is_refused_naming_it_and_the_box(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    """Bygdin, 2026-10-29: a preferred source with no polygon near the domain
    replaced all the water with land and the mesh was written. Refused: a
    usage error naming --prefer, the source, the read region's box and its CRS."""
    lake = write_lake(tmp_path / "lake.geojson", far_lake())
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, *corine(cover, "--prefer", f"lake={lake}"))
    assert code == USAGE, output
    flat = squashed(output)
    assert "--prefer" in flat and "lake.geojson" in flat and "nopolygon" in flat, output
    assert MESH_CRS.replace(" ", "") in flat, output
    assert not target.exists()


def test_a_cached_source_never_fetched_here_is_refused_naming_the_fetch(
    scene: dict[str, Path], tmp_path: Path, cover: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The Bygdin case itself: the cache holds only lakes fetched for another
    region. mesh never fetches, so the refusal names the `rasputin fetch`
    that copies the source for this run's domain."""
    layers = [Layer("features", 31983, [GpkgRow(1, far_lake(), {"magasinNr": "1"})], pk="ID")]
    (scene["cache"] / RESERVOIRS_KEY).mkdir(parents=True)
    write_gpkg(scene["cache"] / RESERVOIRS_KEY / "features.gpkg", layers)
    install_reservoirs(monkeypatch, layer_entry())
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, *corine(cover, "--prefer", f"lake={RESERVOIRS_KEY}"))
    assert code == USAGE, output
    flat = squashed(output)
    assert "--prefer" in flat and "nopolygon" in flat, output
    assert f"rasputinfetch{RESERVOIRS_KEY}--domain" in flat, output
    assert not target.exists()


# ------------------------------------------------------------ 122, test 9: no code 0


def uncovered_share(corners_xy: np.ndarray) -> np.ndarray:
    """Per triangle, its share in no polygon of the map: forest 311 less the
    preferred lake, conifer 312 and the lake; the water left is left out."""
    covered = shapely.union_all(
        [g for c, g in cover_polygons() if c in (311, 312)] + [preferred_lake()]
    )
    triangles = shapely.polygons(corners_xy)
    return 1.0 - shapely.area(shapely.intersection(triangles, covered)) / shapely.area(triangles)


def marked_where_uncovered(marks: np.ndarray, corners_xy: np.ndarray) -> None:
    """`docs/increments/122-lakes-overrule.md`, section 3.3: a triangle whose
    code was filled or whose 0 share moved is marked; so marked exactly where
    the triangle reaches the uncovered area (over 1e-6), never where it does
    not (under 1e-9), and none lies between."""
    rest = uncovered_share(corners_xy)
    assert not ((rest > 1e-9) & (rest <= 1e-6)).any()
    assert (rest > 1e-6).any(), "the scene's left-out water reaches no triangle"
    assert_array_equal_bool(marks == 1, rest > 1e-6)


def assert_array_equal_bool(got: np.ndarray, want: np.ndarray) -> None:
    wrong = np.flatnonzero(got != want)
    assert not len(wrong), f"{len(wrong)} triangles marked wrongly, first {wrong[:10]}"


def test_9_no_triangle_keeps_code_0_in_a_gpkg(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    """122's test 9: the left-out water and the scene's gap (1.9 km², hence
    `--land-cover-partial`) are filled and marked; the record counts them,
    no fractions row is code 0, and the fill's stderr line is printed once."""
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, *corine(cover, "--prefer", f"lake={lake}"))
    assert code == 0, output
    assert squashed(output).count("innopolygontookthenearestnon-waterclass") == 1
    with closing(sqlite3.connect(target)) as conn:
        rows = conn.execute("SELECT code, land_cover_filled FROM triangles ORDER BY fid").fetchall()
        codes = np.array([r[0] for r in rows])
        marks = np.array([r[1] for r in rows])
        assert not (codes == 0).any()
        marked_where_uncovered(marks, corners(conn)[:, :, :2])
        filled = record(conn, "land_cover_filled")
        assert filled is not None and int(filled) == int(marks.sum())
        area = record(conn, "land_cover_filled_m2")
        assert area is not None and float(area) > 0.0
        zero = conn.execute("SELECT count(*) FROM land_cover_fractions WHERE code = 0")
        assert zero.fetchone() == (0,)


def test_9_no_triangle_keeps_code_0_in_a_vtk_or_ply(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    """The same in a `.vtk` (an `int` cell array `land_cover_filled`, 0 on
    every line cell, read back by `read_mesh`) and a `.ply` (a face property
    `uchar land_cover_filled` after `land_cover_code`)."""
    from plyread import read_ply
    from tin_engine.io.vtk_read import read_mesh
    from vtkread import read_vtk

    args = corine(cover, "--prefer", f"lake={lake}")
    vtk_path, ply_path = tmp_path / "x.vtk", tmp_path / "x.ply"
    for target in (vtk_path, ply_path):
        code, output = mesh(scene, target, *args)
        assert code == 0, output
    vtk = read_vtk(vtk_path.read_bytes())
    lines = len(vtk.lines)
    filled = vtk.cell_array("land_cover_filled")
    assert filled.type_name == "int"
    marks = np.asarray(filled.values).astype(int)
    assert not marks[:lines].any(), "a line cell is marked"
    codes = np.asarray(vtk.cell_array("land_cover_code").values)
    assert not (codes[lines:] == 0).any()
    points = np.asarray(vtk.points)
    corners_xy = np.array([points[np.asarray(p)][:, :2] for p in vtk.polygons])
    marked_where_uncovered(marks[lines:], corners_xy)
    arrays = read_mesh(vtk_path.read_bytes())
    assert arrays.land_cover_filled is not None
    assert (np.asarray(arrays.land_cover_filled) == marks[lines:]).all()
    header, elements = read_ply(ply_path.read_bytes())
    names = [p.name for p in header.element("face").properties]
    assert names.index("land_cover_filled") == names.index("land_cover_code") + 1, names
    (kind,) = [
        p.type_name for p in header.element("face").properties if p.name == "land_cover_filled"
    ]
    assert kind == "uchar"
    face = elements["face"]
    assert not (np.asarray(face["land_cover_code"]) == 0).any()
    assert int(np.asarray(face["land_cover_filled"]).sum()) == int(marks.sum())
