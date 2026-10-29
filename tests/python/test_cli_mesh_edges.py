"""`rasputin mesh --land-cover ... --prefer ... --lines KEY=SOURCE --edges SELECTOR
[--shares]`: one land-cover map, and the edges the mesh follows chosen from it.

`docs/increments/111-one-land-cover-map.md`, *5. Tests*, tests 2 to 11 (test
1, `land_cover_map.build`, is `test_land_cover_map.py`'s). The scene is
`test_cli_mesh_prefer.py`'s: `cover_polygons` (forest 311 and 312, water 512)
in `cover.geojson`, the preferred lake in `lake.geojson`, every clean-up
switch 0 (`AS_READ`).

THE ORACLE for a line length: the `lines` layer of the `.gpkg` written, one
row per constraint edge with its plan length and one 0/1 column per edge
name; a class border is a row with `land_cover` set, a river one with
`river` set (the domain outline carries no name, so it is in neither). The
design's figures are to 0.1 m and checked on today's tree with the
equivalent `--features` run to 1e-6 m (9,082.700 m), so bound `LENGTH` =
0.01 m, as the design states. Areas from shapely on the polygons as given,
to 1 m², as the design states. Shares: case 7's oracle (`test_cli_mesh_prefer`),
bound `SHARE` = 1e-6.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- Every refusal is exit 2 with no output file, naming the option the design
  names; no assertion on which phases ran, except the cut-run refusals (no
  DEM block assembled: `cli.assemble_dem` is never called).
- `--shares` without `--land-cover` names both options.
- `land_cover_overlap` states the overlap area in m² (a number in the text
  within 1 of the area) and names the file kept by its file name.
- A class-border run without `--shares` writes no `land_cover_fractions`
  table and no `cells` layer; its `triangles` layer keeps the `code` column.
  A run with no land-cover map writes `code` NULL on every triangle, as the
  GeoPackage writer does for a mesh without codes.
- Test 9's codes per triangle are pinned as a digest of the `code` column in
  `fid` order (int32, little endian) with the counts per code, from today's
  tree (`39595f5f`) with `--features cover.geojson --features-rules corine`
  and the clean-up switches 0, no `--prefer` (today's tree refuses a coded
  `--features` beside `--land-cover`).
- A bowtie in a `--land-cover` file ends the run with exit 0 and no
  exception escaping `rasputin mesh`; one in a pair overlapping far outside
  the read box is counted outside by the edges step, as its partner is.
- `--threads N` reaches the repairing `build` call as `threads=N`; the
  `--prefer` call (`fix=False`) repairs nothing and is not checked.
- A river file is a GeoJSON `FeatureCollection` of one `LineString` with no
  properties, in the mesh CRS.
- Ola's question 4 (a county line with one class on both sides) is taken at
  its default, not followed: `test_10_a_county_line_with_one_class_on_both_sides_is_not_followed`
  is the one test a veto changes.
"""

from __future__ import annotations

import hashlib
import os
import pwd
import re
import socket
import sqlite3
import traceback
from collections.abc import Sequence
from contextlib import closing
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from numpy.testing import assert_allclose, assert_array_equal
from shapely.geometry import LineString, Polygon, box
from shapely.geometry.base import BaseGeometry

import tin_engine.cli as cli
from class_fixtures import KEY
from cli_helpers import USAGE, invoke, runner, squashed
from dataset_fixtures import write_sources
from feature_fixtures import Feat, write_geojson
from fraction_fixtures import MESH_CRS, dem_origin, domain_ring, dominant
from legend_fixtures import basenames, coded_map, write
from test_cli_land_cover_polygons import cover_polygons, write_cover_geojson, write_cover_gpkg
from test_cli_mesh_cache import no_network
from test_cli_mesh_cells import check_sums, column, mesh, scene, written
from test_cli_mesh_prefer import preferred_lake, share
from test_cli_mesh_reservoirs import phases, record, write_reservoirs
from test_legend_name import RENAME_TO_CLASSES
from tin_engine.cli import app
from vector_fixtures import vector_entry
from vtkread import read_vtk

__all__ = ["phases", "scene"]  # fixtures, used by name

#: On a line length, m: the design's figures are to 0.1 m, measured to 1e-6 m.
LENGTH = 0.01
#: The river's length inside the domain, m: the design gives 2,448.4 m, its
#: figure to 0.1 m; shapely on the line gives 2,448.3688 m, so the figure is
#: pinned to the second decimal and held to `LENGTH`, as test 2's are.
RIVER = 2448.37
#: On an area, m²: the design's figures are to 1 m².
AREA = 1.0
#: On a share: case 7's bound (`test_cli_mesh_prefer.py`).
SHARE = 1e-6
#: Every clean-up switch off, as the design measured.
AS_READ = ("--edges-repair", "0", "--edges-tolerance", "0", "--edges-outline-snap", "0")
#: Arable land: a class `--prefer lake=` does not replace, and no water class.
ARABLE = 211


def classes(*files: Path | str) -> tuple[str, ...]:
    given = [a for f in files for a in ("--land-cover", str(f))]
    return (*given, "--land-cover-classes", "corine")


def edges(*selectors: str) -> tuple[str, ...]:
    return tuple(a for s in selectors for a in ("--edges", s))


def river_line() -> LineString:
    x0, y0 = dem_origin()
    return LineString(
        [(x0 + 200.3, y0 - 2500.7), (x0 + 1350.9, y0 - 1100.2), (x0 + 2900.6, y0 - 300.1)]
    )


def planted_box() -> Polygon:
    """`P`: arable land inside forest 311, 379.8 m west of the preferred lake."""
    x0, y0 = dem_origin()
    return box(x0 + 420.3, y0 - 1319.8, x0 + 620.4, y0 - 1119.1)


def domain() -> Polygon:
    return Polygon(domain_ring())


@pytest.fixture
def cover(tmp_path: Path) -> Path:
    return write_cover_geojson(tmp_path / "cover.geojson", cover_polygons())


@pytest.fixture
def lake(tmp_path: Path) -> Path:
    return write_geojson(
        tmp_path / "lake.geojson", [Feat(1, preferred_lake(), {"vatnLnr": "1"})], crs=MESH_CRS
    )


@pytest.fixture
def river(tmp_path: Path) -> Path:
    return write_geojson(tmp_path / "river.geojson", [Feat(1, river_line(), {})], crs=MESH_CRS)


def prefer(lake: Path) -> tuple[str, ...]:
    return ("--prefer", f"lake={lake}")


def lengths(conn: sqlite3.Connection) -> dict[str, float]:
    """Per edge name, the plan length of the `lines` rows that carry it."""
    names = [r[1] for r in conn.execute("PRAGMA table_info(lines)")]
    flags = [n for n in names if n not in ("fid", "geom", "length_m")]
    query = 'SELECT COALESCE(SUM(length_m), 0) FROM lines WHERE "{}" = 1'
    return {n: float(conn.execute(query.format(n)).fetchone()[0]) for n in flags}


def tables(conn: sqlite3.Connection) -> set[str]:
    return {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}


def has_shares(conn: sqlite3.Connection) -> bool:
    return "land_cover_fractions" in tables(conn) or "cells" in tables(conn)


def run_gpkg(scene: dict[str, Path], target: Path, *args: str) -> sqlite3.Connection:
    return written(scene, target, *args)


def check_case7_shares(conn: sqlite3.Connection) -> None:
    """Case 7's oracle: the lake share is the preferred lake's, the forest
    share the forest less the lake, per triangle as written."""
    assert_allclose(column(conn, "cells", "lake"), share(conn, preferred_lake()), atol=SHARE)
    forest = shapely.union_all([g for c, g in cover_polygons() if c in (311, 312)])
    expected = share(conn, forest.difference(preferred_lake()))
    assert_allclose(column(conn, "cells", "forest"), expected, atol=SHARE)
    check_sums(conn)


def refused(scene: dict[str, Path], tmp_path: Path, *args: str) -> str:
    """Exit 2, no file; the `basenames` of the output with no whitespace."""
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, *args)
    assert code == USAGE, output
    assert "No such option" not in output, output  # refused by the run, not the parser
    assert not target.exists()
    return squashed(basenames(output))


# ------------------------------------------------------------ test 2: which borders


@pytest.mark.parametrize(
    ("selector", "expected"),
    [("all", 9082.7), ("lake", 2600.4), ("water", 2600.4), ("forest", 7482.8), ("311", 5600.0)],
)
def test_2_the_chosen_borders_are_the_edges(
    scene: dict[str, Path],
    tmp_path: Path,
    cover: Path,
    lake: Path,
    selector: str,
    expected: float,
) -> None:
    args = (*classes(cover), *prefer(lake), *edges(selector), *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        found = lengths(conn)
    assert found["land_cover"] == pytest.approx(expected, abs=LENGTH), found
    if selector in ("all", "lake", "water"):  # the lake is the only water class left
        assert found["water"] == pytest.approx(2600.4, abs=LENGTH), found
    else:
        assert found.get("water", 0.0) == 0.0, found


def test_2_the_premises_by_shapely() -> None:
    """The design's figures from the polygons alone, inside the domain."""
    dom, polygons = domain(), dict(cover_polygons())
    forest = polygons[311].difference(preferred_lake())
    shores = [forest, polygons[312], preferred_lake()]
    every = shapely.union_all([g.boundary for g in shores]).intersection(dom)
    assert every.length == pytest.approx(9082.7, abs=0.05)
    assert forest.boundary.intersection(dom).length == pytest.approx(5600.0, abs=0.05)
    assert preferred_lake().boundary.intersection(dom).length == pytest.approx(2600.4, abs=0.05)
    trees = shapely.union_all([forest.boundary, polygons[312].boundary]).intersection(dom)
    assert trees.length == pytest.approx(7482.8, abs=0.05)


def test_2_no_edges_without_edges_and_the_run_gives_shares(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    args = (*classes(cover), *prefer(lake), *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        assert lengths(conn).get("land_cover", 0.0) == 0.0
        assert has_shares(conn)
        check_case7_shares(conn)


# ------------------------------------------------------------ test 3: a river beside the map


def test_3_a_river_is_followed_and_the_shares_stay_case_7s(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path, river: Path
) -> None:
    assert river_line().intersection(domain()).length == pytest.approx(RIVER, abs=LENGTH)
    args = (*classes(cover), *prefer(lake), "--lines", f"river={river}", *edges("river"))
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args, "--shares", *AS_READ)) as conn:
        found = lengths(conn)
        assert found["river"] == pytest.approx(RIVER, abs=LENGTH), found
        assert found.get("land_cover", 0.0) == 0.0, found
        check_case7_shares(conn)


def test_3_a_line_no_edges_selector_takes_is_refused(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path, river: Path
) -> None:
    args = (*classes(cover), *prefer(lake), "--lines", f"river={river}", *edges("lake"))
    flat = refused(scene, tmp_path, *args)
    assert "--lines" in flat and "river" in flat and "--edges" in flat, flat


# ------------------------------------------------------------ test 4: case 7 with every border


def test_4_case_7s_shares_with_every_border_and_shares(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    args = (*classes(cover), *prefer(lake), *edges("all"), "--shares", *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        assert lengths(conn)["land_cover"] == pytest.approx(9082.7, abs=LENGTH)
        check_case7_shares(conn)


# ------------------------------------------------------------ test 5: refusals


def test_5_features_is_an_unknown_option(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    code, output = mesh(scene, tmp_path / "x.gpkg", "--features", str(cover))
    assert code == USAGE, output
    assert "No such option: --features" in output, output


@pytest.mark.parametrize(
    "old",
    ["--features-rules", "--features-layer", "--features-crs", "--features-tolerance",
     "--features-repair", "--features-outline-snap", "--land-cover-rules"],
)  # fmt: skip
def test_5_the_old_options_are_unknown(old: str) -> None:
    code, output = invoke("mesh", old, "x")
    assert code == USAGE, output
    assert f"No such option: {old}" in output, output


def test_5_prefer_with_a_line_key_names_lines(
    scene: dict[str, Path], tmp_path: Path, cover: Path, river: Path
) -> None:
    flat = refused(scene, tmp_path, *classes(cover), "--prefer", f"river={river}")
    assert "--prefer" in flat and "river" in flat and "--lines" in flat, flat


def test_5_lines_with_an_area_key_names_prefer(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    args = (*classes(cover), "--lines", f"lake={lake}", *edges("all"))
    flat = refused(scene, tmp_path, *args)
    assert "--lines" in flat and "lake" in flat and "--prefer" in flat, flat


def test_5_a_polygon_source_under_lines_is_refused(
    scene: dict[str, Path], tmp_path: Path, lake: Path
) -> None:
    flat = refused(scene, tmp_path, "--lines", f"river={lake}", *edges("river"))
    assert "--lines" in flat and "lake.geojson" in flat, flat


def test_5_a_line_source_under_prefer_is_refused(
    scene: dict[str, Path], tmp_path: Path, cover: Path, river: Path
) -> None:
    flat = refused(scene, tmp_path, *classes(cover), "--prefer", f"lake={river}")
    assert "--prefer" in flat and "river.geojson" in flat, flat


def test_5_a_second_source_for_one_lines_key_is_refused(
    scene: dict[str, Path], tmp_path: Path, river: Path
) -> None:
    other = write_geojson(tmp_path / "other.geojson", [Feat(1, river_line(), {})], crs=MESH_CRS)
    args = ("--lines", f"river={river}", "--lines", f"river={other}", *edges("river"))
    flat = refused(scene, tmp_path, *args)
    assert "--lines" in flat and "river" in flat, flat


def test_5_an_unknown_selector_lists_the_accepted_ones(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    flat = refused(scene, tmp_path, *classes(cover), *edges("wetlands"))
    assert "--edges" in flat and "wetlands" in flat, flat
    for accepted in ("all", "water", "forest", "lake"):
        assert accepted in flat, (accepted, flat)


def test_5_an_area_selector_without_class_polygons_is_refused(
    scene: dict[str, Path], tmp_path: Path, river: Path
) -> None:
    flat = refused(scene, tmp_path, "--lines", f"river={river}", *edges("river", "lake"))
    assert "--edges" in flat and "lake" in flat, flat


def test_5_an_area_selector_over_a_class_raster_is_refused(
    scene: dict[str, Path], tmp_path: Path
) -> None:
    args = ("--land-cover", KEY, "--land-cover-classes", "mapbiomas-c11", *edges("forest"))
    flat = refused(scene, tmp_path, *args)
    assert "--edges" in flat, flat


def test_5_a_line_key_without_its_lines_is_refused(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    flat = refused(scene, tmp_path, *classes(cover), *edges("river"))
    assert "--edges" in flat and "river" in flat, flat


@pytest.mark.parametrize("category", ["all", "water", "river"])
def test_5_a_class_file_with_a_category_named_as_a_selector_is_refused(
    scene: dict[str, Path], tmp_path: Path, cover: Path, category: str
) -> None:
    body = (
        '[[class]]\ncode = 311\nlabel = "Forest"\ncolour = "#1f8d49"\n'
        f'category = "{category}"\n\n'
        '[[class]]\ncode = 512\nlabel = "Lake"\ncolour = "#2532e4"\nwater = true\n'
        'category = "lake"\n'
    )
    mine = write(tmp_path / "mine.toml", coded_map("mine", attribute="Code_18", classes=body))
    flat = refused(scene, tmp_path, "--land-cover", str(cover), "--land-cover-classes", str(mine))
    assert "categor" in flat and category in flat.replace("mine.toml", ""), flat


def test_5_shares_without_land_cover_is_refused(
    scene: dict[str, Path], tmp_path: Path, river: Path
) -> None:
    flat = refused(scene, tmp_path, "--lines", f"river={river}", *edges("river"), "--shares")
    assert "--shares" in flat and "--land-cover" in flat, flat


@pytest.mark.parametrize("option", ["--land-cover-layer", "--land-cover-crs"])
def test_5_a_land_cover_file_option_without_land_cover_is_refused(
    scene: dict[str, Path], tmp_path: Path, option: str
) -> None:
    value = "cover" if option.endswith("layer") else MESH_CRS
    flat = refused(scene, tmp_path, option, value)
    assert option in flat, flat


@pytest.mark.parametrize("option", ["--land-cover-layer", "--land-cover-crs"])
def test_5_a_land_cover_file_option_with_a_catalogue_key_is_refused(
    scene: dict[str, Path], tmp_path: Path, option: str
) -> None:
    gpkg = write_cover_gpkg(tmp_path / "cover.gpkg", cover_polygons())
    config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
    write_sources(config, {"cover-test": vector_entry(gpkg, "Test cover", classes="corine")})
    value = "cover" if option.endswith("layer") else MESH_CRS
    flat = refused(scene, tmp_path, "--land-cover", "cover-test", option, value, *edges("all"))
    assert option in flat, flat


def test_5_a_layer_for_a_file_that_is_no_geopackage_is_refused(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    flat = refused(scene, tmp_path, *classes(cover), "--land-cover-layer", "cover", *edges("all"))
    assert "--land-cover-layer" in flat, flat


def test_5_a_crs_that_is_not_the_files_own_is_refused(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    args = (*classes(cover), "--land-cover-crs", "EPSG:4326", *edges("all"))
    flat = refused(scene, tmp_path, *args)
    assert "--land-cover-crs" in flat, flat


def test_5_a_sources_file_with_a_rules_key_is_refused(
    scene: dict[str, Path], tmp_path: Path
) -> None:
    gpkg = write_cover_gpkg(tmp_path / "cover.gpkg", cover_polygons())
    config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
    write_sources(config, {"cover-test": vector_entry(gpkg, "Test cover", rules="corine")})
    code, output = mesh(scene, tmp_path / "x.gpkg", "--land-cover", "cover-test", *edges("all"))
    assert code != 0, output
    assert "rules" in output and "cover-test" in output, output
    assert RENAME_TO_CLASSES.search(output), output
    assert not (tmp_path / "x.gpkg").exists()


def test_5_a_catalogue_entry_and_classes_that_differ_are_refused(
    scene: dict[str, Path], tmp_path: Path
) -> None:
    gpkg = write_cover_gpkg(tmp_path / "cover.gpkg", cover_polygons())
    config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
    write_sources(config, {"cover-test": vector_entry(gpkg, "Test cover", classes="corine")})
    args = ("--land-cover", "cover-test", "--land-cover-classes", "ar50", *edges("all"))
    flat = refused(scene, tmp_path, *args)
    assert "--land-cover" in flat, flat


# ------------------------------------------------------------ test 6: the MCP tool


def test_6_the_mcp_mesh_tool_follows_the_cli() -> None:
    pytest.importorskip("mcp")
    from tin_engine.cli import app
    from tin_engine.mcp import schema as mcp_schema
    from tin_engine.mcp.server import INSTRUCTIONS

    mesh_tool = {d.tool_name: d for d in mcp_schema.derive(app)}["mesh"]
    properties = set(mesh_tool.tool().input_schema["properties"])
    assert {"edges", "lines", "shares", "land_cover_classes"} <= properties, sorted(properties)
    assert {"land_cover_layer", "land_cover_crs", "edges_tolerance"} <= properties
    assert not {p for p in properties if p.startswith("features")}, sorted(properties)
    assert "land_cover_rules" not in properties
    for word in ("features_rules", "features_layer"):
        assert word not in INSTRUCTIONS, word
    for word in ("edges", "lines", "shares", "prefer"):
        assert word in INSTRUCTIONS, word
    arguments = {"out": "x.vtk", "dem": ["d.tif"], "edges": ["lake", "river"],
                 "lines": ["river=r.geojson"], "shares": True}  # fmt: skip
    argv = " ".join(mesh_tool.argv(mesh_tool.model.model_validate(arguments)))
    assert "--edges lake" in argv and "--edges river" in argv, argv
    assert "--lines river=r.geojson" in argv and "--shares" in argv, argv


# ------------------------------------------------------------ test 7: Sula, Ola's data

#: The user's own home: the suite's `HOME` is a fresh directory (`conftest.py`).
OWN_HOME = Path(pwd.getpwuid(os.getuid()).pw_dir)
SULA = OWN_HOME / "rasputin-output/mcp/ola/work/sula-73-27-0/sula-official-catchment.geojson"
DTM10 = OWN_HOME / "projects/rasputin_data/DTM10_UTM33_20260925"
CACHE = OWN_HOME / ".cache/rasputin"
needs_sula = pytest.mark.skipif(
    not (SULA.exists() and DTM10.is_dir()
         and (CACHE / "clc2018").is_dir() and (CACHE / "nve-lakes").is_dir()),
    reason="Ola's Sula domain, DTM10 and cached clc2018 and nve-lakes are not here",
)  # fmt: skip


@needs_sula
def test_7_sula_the_direct_command_gives_sols_mesh(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Local only: section 0's command against Sol's record of 28 October
    (3,271 triangles, 30.4411 km², lake 8.935 %), the DEM given by its path
    and the cache by `--cache`, the network closed."""
    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    target = tmp_path / "sula.gpkg"
    code, output = invoke(
        "mesh", "--out", str(target), "--dem", str(DTM10), "--tolerance", "20",
        "--domain", str(SULA), "--out-crs", "EPSG:25833", "--cache", str(CACHE),
        "--land-cover", "clc2018", "--prefer", "lake=nve-lakes", "--edges", "all", "--shares",
        "--edges-repair", "0.1", "--edges-tolerance", "20", "--edges-outline-snap", "0",
        "--threads", "4", "--z-unit", "m",
    )  # fmt: skip
    assert code == 0, output
    with closing(sqlite3.connect(target)) as conn:
        area = column(conn, "cells", "area_m2")
        lake_share = column(conn, "cells", "lake")
    assert len(area) == 3271
    assert area.sum() / 1e6 == pytest.approx(30.4411, abs=5e-5)
    assert 100 * (area * lake_share).sum() / area.sum() == pytest.approx(8.935, abs=5e-4)


# ------------------------------------------------------------ test 8: no land-cover map


@pytest.mark.parametrize("selector", ["river", "all"])
def test_8_lines_alone_give_the_river_and_nothing_else(
    scene: dict[str, Path], tmp_path: Path, river: Path, selector: str
) -> None:
    args = ("--lines", f"river={river}", *edges(selector))
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        found = lengths(conn)
        assert found["river"] == pytest.approx(RIVER, abs=LENGTH), found
        assert set(found) == {"river"}, found
        coded = conn.execute("SELECT count(*) FROM triangles WHERE code IS NOT NULL").fetchone()
        assert coded == (0,), coded
        assert not has_shares(conn)


# ------------------------------------------------------------ test 9: when the shares are computed

#: Today's codes per triangle for `--features cover.geojson --features-rules
#: corine` and the clean-up switches 0 (tree 39595f5f): 1,691 triangles.
TODAYS_CODES = {0: 767, 311: 496, 312: 87, 512: 341}
TODAYS_DIGEST = "ee7aac588b95228c74272e43d8f7ce89b2f62ca957bc9d846655f909b491978e"


def codes_of(conn: sqlite3.Connection) -> np.ndarray:
    rows = conn.execute("SELECT code FROM triangles ORDER BY fid").fetchall()
    return np.array([r[0] for r in rows], dtype="<i4")


def test_9_land_cover_alone_writes_the_shares(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path, phases: list[str]
) -> None:
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *classes(cover), *prefer(lake))) as conn:
        assert "land_cover_fractions" in tables(conn)
        check_case7_shares(conn)
    assert "land cover fractions" in phases, phases


def test_9_class_borders_alone_give_todays_codes_and_no_shares(
    scene: dict[str, Path], tmp_path: Path, cover: Path, phases: list[str]
) -> None:
    stats = tmp_path / "stats.md"
    args = (*classes(cover), *edges("all"), *AS_READ, "--stats", str(stats))
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        assert not has_shares(conn)
        codes = codes_of(conn)
        texts = [str(v) for (v,) in conn.execute("SELECT value FROM run").fetchall()]
    assert any("not computed" in t and "add --shares" in t for t in texts), texts
    assert "land cover fractions" not in phases, phases
    report = stats.read_text(encoding="utf-8")
    assert "| edges |" in report and "land cover fractions" not in report, report
    found, counts = np.unique(codes, return_counts=True)
    assert dict(zip(found.tolist(), counts.tolist(), strict=True)) == TODAYS_CODES
    assert hashlib.sha256(codes.tobytes()).hexdigest() == TODAYS_DIGEST


def test_9_shares_asked_give_the_table_and_codes_from_it(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    args = (*classes(cover), *prefer(lake), *edges("all"), "--shares", *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        assert "land_cover_fractions" in tables(conn)
        rows = conn.execute("SELECT triangle, code, fraction FROM land_cover_fractions").fetchall()
        fids = [r[0] for r in conn.execute("SELECT fid FROM triangles ORDER BY fid").fetchall()]
        codes = codes_of(conn)
    table: dict[int, dict[int, float]] = {}
    for triangle, code, fraction in rows:
        table.setdefault(triangle, {})[code] = fraction
    assert_array_equal(codes, [dominant(table[f]) for f in fids])


def test_9_reservoirs_with_class_borders_need_shares(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    res = write_reservoirs(tmp_path / "res.geojson", wgs84=False)
    flat = refused(scene, tmp_path, *classes(cover), *edges("all"), "--reservoirs", str(res))
    assert "--reservoirs" in flat and "--shares" in flat, flat


def test_9_a_river_alone_gives_no_share_table(
    scene: dict[str, Path], tmp_path: Path, river: Path
) -> None:
    args = ("--lines", f"river={river}", *edges("river"))
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        assert not has_shares(conn)


# ------------------------------------------------------------ test 10: several land-cover files


def split(tmp_path: Path, a: Sequence[tuple[int, BaseGeometry]],
          b: Sequence[tuple[int, BaseGeometry]]) -> tuple[Path, Path]:  # fmt: skip
    return (write_cover_geojson(tmp_path / "a.geojson", a),
            write_cover_geojson(tmp_path / "b.geojson", b))  # fmt: skip


def areas_by_code(conn: sqlite3.Connection) -> dict[int, float]:
    """Triangle area per code, inside the mesh (the domain)."""
    rows = conn.execute("SELECT code, SUM(area_m2) FROM triangles GROUP BY code").fetchall()
    return {int(c): float(a) for c, a in rows}


def overlap_figure(conn: sqlite3.Connection) -> str:
    text = record(conn, "land_cover_overlap")
    assert text is not None, "no land_cover_overlap row"
    return text


def numbers(text: str) -> list[float]:
    return [float(n.replace(",", "")) for n in re.findall(r"\d[\d,]*(?:\.\d+)?", text)]


def test_10_the_planted_box_premises() -> None:
    polygons = dict(cover_polygons())
    p = planted_box()
    assert p.area == pytest.approx(40160.1, abs=0.05)
    assert p.length == pytest.approx(801.6, abs=0.05)
    assert p.distance(preferred_lake()) == pytest.approx(379.8, abs=0.05)
    assert polygons[311].contains(p)
    # Planted check: "overlaps" misses a polygon wholly inside another.
    tree = shapely.STRtree([polygons[311]])
    assert len(tree.query(p, predicate="overlaps")) == 0
    assert len(tree.query(p, predicate="intersects")) == 1


def test_10_a_box_wholly_inside_an_earlier_files_polygon_is_cut_away(
    scene: dict[str, Path], tmp_path: Path, lake: Path
) -> None:
    polygons = dict(cover_polygons())
    a, b = split(tmp_path, [(311, polygons[311]), (312, polygons[312])],
                 [(512, polygons[512]), (ARABLE, planted_box())])  # fmt: skip
    args = (*classes(a, b), *prefer(lake), *edges("all"), "--shares", *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        found = lengths(conn)
        overlap = overlap_figure(conn)
        source = record(conn, "land_cover_source")
        assert ARABLE not in {
            c for (c,) in conn.execute("SELECT DISTINCT code FROM land_cover_fractions")
        }
        check_case7_shares(conn)
    assert found["land_cover"] == pytest.approx(9082.7, abs=LENGTH), found
    assert found["water"] == pytest.approx(2600.4, abs=LENGTH), found
    assert any(abs(n - 40160.1) <= AREA for n in numbers(overlap)), overlap
    assert "a.geojson" in overlap, overlap
    assert source is not None and "a.geojson" in source and "b.geojson" in source, source


def test_10_two_identical_polygons_the_first_given_keeps_the_area(
    scene: dict[str, Path], tmp_path: Path, lake: Path
) -> None:
    """A copy of 312 coded arable land in the second file: the overlap is all
    of 312 inside the domain, and no share is arable."""
    polygons = dict(cover_polygons())
    a, b = split(tmp_path, cover_polygons(), [(ARABLE, polygons[312])])
    args = (*classes(a, b), *prefer(lake), *edges("all"), "--shares", *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        found = lengths(conn)
        overlap = overlap_figure(conn)
        assert ARABLE not in {
            c for (c,) in conn.execute("SELECT DISTINCT code FROM land_cover_fractions")
        }
        check_case7_shares(conn)
    inside = polygons[312].intersection(domain()).area
    assert any(abs(n - inside) <= AREA for n in numbers(overlap)), (inside, overlap)
    assert "a.geojson" in overlap, overlap
    assert found["land_cover"] == pytest.approx(9082.7, abs=LENGTH), found


def test_10_a_bowtie_in_the_later_file_is_repaired_not_a_traceback(
    scene: dict[str, Path], tmp_path: Path, lake: Path
) -> None:
    """Code review round 1: `P`'s corners joined as a bowtie (a ring that
    crosses itself at `P`'s centre), coded arable in the second file, made
    `land_cover_map._joined` throw GEOSException. Repaired as the shares step
    repairs, it is `P`'s west and east triangles, wholly inside forest 311:
    the first file keeps all of it, and the mesh is test 10's box run."""
    p = planted_box()
    x0, y0, x1, y1 = p.bounds
    tie = Polygon([(x0, y0), (x1, y1), (x1, y0), (x0, y1)])
    assert not tie.is_valid and dict(cover_polygons())[311].contains(shapely.make_valid(tie))
    polygons = dict(cover_polygons())
    a, b = split(tmp_path, [(311, polygons[311]), (312, polygons[312])],
                 [(512, polygons[512]), (ARABLE, tie)])  # fmt: skip
    args = (*classes(a, b), *prefer(lake), *edges("all"), "--shares", *AS_READ)
    target = tmp_path / "x.gpkg"
    # `test_cli_mesh_cells.mesh`'s arguments, run here to keep the exception.
    given = ["--dem", str(scene["dem"]), "--domain", str(scene["domain"]), "--tolerance", "3",
             "--cache", str(scene["cache"]), "--out", str(target)]  # fmt: skip
    result = runner.invoke(app, ["mesh", *given, *args])
    # No exception escaped the run, and none was printed as a traceback.
    assert result.exception is None, "".join(traceback.format_exception(result.exception))
    assert result.exit_code == 0, result.output
    assert "Traceback" not in result.output and "GEOSException" not in result.output
    with closing(sqlite3.connect(target)) as conn:
        found = lengths(conn)
        overlap = overlap_figure(conn)
        assert ARABLE not in {
            c for (c,) in conn.execute("SELECT DISTINCT code FROM land_cover_fractions")
        }
        check_case7_shares(conn)
    half = p.area / 2  # 20,080.0 m²
    assert any(abs(n - half) <= AREA for n in numbers(overlap)), (half, overlap)
    assert "a.geojson" in overlap, overlap
    assert found["land_cover"] == pytest.approx(9082.7, abs=LENGTH), found


def lone_bowtie() -> Polygon:
    """A 200 m bowtie coded arable, crossing itself at its centre: its two
    triangles (20,000 m²) lie in forest 311, the east one half on the
    preferred lake. It shares area with no polygon of another file."""
    x0, y0 = dem_origin()
    a, b = x0 + 900.1, y0 - 1300.3
    return Polygon([(a, b), (a + 200, b + 200), (a + 200, b), (a, b + 200)])


def test_10_the_lone_bowtie_premises() -> None:
    tie, polygons = lone_bowtie(), dict(cover_polygons())
    assert not tie.is_valid
    whole = shapely.make_valid(tie)
    assert whole.area == pytest.approx(20000.0, abs=AREA)
    assert polygons[311].intersection(whole).area == pytest.approx(20000.0, abs=AREA)
    assert preferred_lake().intersection(whole).area > 0


def run_kept(scene: dict[str, Path], target: Path, *args: str) -> str:
    """`test_cli_mesh_cells.mesh`'s arguments, run here to keep the exception:
    none escapes, none is printed as a traceback, and the exit is 0. The
    output, stderr included."""
    given = ["--dem", str(scene["dem"]), "--domain", str(scene["domain"]), "--tolerance", "3",
             "--cache", str(scene["cache"]), "--out", str(target)]  # fmt: skip
    result = runner.invoke(app, ["mesh", *given, *args])
    assert result.exception is None, "".join(traceback.format_exception(result.exception))
    assert result.exit_code == 0, result.output
    assert "Traceback" not in result.output and "GEOSException" not in result.output
    return result.output


@pytest.mark.parametrize("with_prefer", [False, True], ids=["alone", "prefer"])
def test_10_a_bowtie_in_a_single_file_is_repaired_once_and_counted(
    scene: dict[str, Path], tmp_path: Path, lake: Path, with_prefer: bool
) -> None:
    """Code review round 2: one `--land-cover` file holding a bowtie in no
    overlapping pair stayed unrepaired, and `--prefer lake=` made
    `preferred.replace` throw GEOSException in `_cut`. Repaired once where
    the map is built, the run ends with exit 0 and the record counts the
    one repair (`land_cover_repaired`), with or without `--prefer`."""
    f = write_cover_geojson(tmp_path / "a.geojson", [*cover_polygons(), (ARABLE, lone_bowtie())])
    target = tmp_path / "x.gpkg"
    run_kept(scene, target, *classes(f), *(prefer(lake) if with_prefer else ()), "--shares")
    with closing(sqlite3.connect(target)) as conn:
        assert record(conn, "land_cover_repaired") == "1"


def test_10_a_bowtie_in_a_single_file_with_edge_cleanup_is_not_a_traceback(
    scene: dict[str, Path], tmp_path: Path
) -> None:
    """Code review round 2: `--edges` with a clean-up switch on the same file
    made the clean-up's clip (`feature_input.py`, `clean`) throw on the
    unrepaired bowtie, on master too. Repaired where the map is built, the
    run ends with exit 0 and follows the class borders."""
    f = write_cover_geojson(tmp_path / "a.geojson", [*cover_polygons(), (ARABLE, lone_bowtie())])
    target = tmp_path / "x.gpkg"
    run_kept(scene, target, *classes(f), *edges("all"), "--edges-repair", "0.1")
    with closing(sqlite3.connect(target)) as conn:
        assert lengths(conn)["land_cover"] > 0


#: `cli._edges`'s stderr line.
EDGES_LINE = re.compile(
    r"\bedges: (?P<n>\d+) kept \((?P<c>\d+) cut at the domain outline\), "
    r"(?P<o>\d+) outside the domain, (?P<e>\d+) empty\b"
)


def far_pair() -> tuple[Polygon, Polygon]:
    """A 200 m square coded 311 and a 200 m bowtie coded arable over it, 50 km
    east of the scene: far outside the read box (`source_region` widened by
    `land_cover_map.MARGIN`), and sharing 20,000 m² with each other."""
    x0, y0 = dem_origin()
    a, b = x0 + 50_000, y0 - 1000
    tie = Polygon([(a, b), (a + 200, b + 200), (a + 200, b), (a, b + 200)])
    return box(a, b, a + 200, b + 200), tie


def test_10_the_far_pair_premises() -> None:
    square, tie = far_pair()
    assert not tie.is_valid
    assert square.intersection(shapely.make_valid(tie)).area == pytest.approx(20000.0, abs=AREA)
    assert square.distance(domain()) > 45_000 and tie.distance(domain()) > 45_000


def outside_count(output: str) -> int:
    found = [int(m["o"]) for m in EDGES_LINE.finditer(output)]
    assert len(found) == 1, output
    return found[0]


@pytest.mark.parametrize("shares", [False, True], ids=["edges", "edges-and-shares"])
def test_10_an_overlapping_pair_with_a_bowtie_far_outside_the_region_is_counted_outside(
    scene: dict[str, Path], tmp_path: Path, shares: bool
) -> None:
    """Code review round 3: two `--land-cover` files whose polygons overlap
    50 km outside the region, one a bowtie. `_fixed` keeps a polygon the read
    box misses whole and unchecked, and `_joined` intersected the pair, so
    the run ended in GEOSException. Left out of the overlap join, the run
    ends with exit 0 and the edges step counts both outside the domain: two
    more than the same files without them."""
    polygons = dict(cover_polygons())
    near_a, near_b = [(311, polygons[311]), (312, polygons[312])], [(512, polygons[512])]
    square, tie = far_pair()
    extra = (*edges("all"), *AS_READ, *(("--shares",) if shares else ()))
    (tmp_path / "far").mkdir()
    (tmp_path / "near").mkdir()
    a, b = split(tmp_path / "far", [*near_a, (311, square)], [*near_b, (ARABLE, tie)])
    with_far = run_kept(scene, tmp_path / "far.gpkg", *classes(a, b), *extra)
    a, b = split(tmp_path / "near", near_a, near_b)
    without = run_kept(scene, tmp_path / "near.gpkg", *classes(a, b), *extra)
    assert outside_count(with_far) == outside_count(without) + 2, (with_far, without)


def test_10_threads_reach_the_map_build(
    scene: dict[str, Path], tmp_path: Path, lake: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Code review round 3: `land_cover_map` checked the cut polygons in a
    fixed 8 threads whatever `--threads` said. `--threads 3` reaches every
    `build` that repairs (each call but `fix=False`'s) as ``threads=3``."""
    calls: list[dict[str, Any]] = []
    real = cli.build

    def spy(*args: Any, **kwargs: Any) -> Any:
        calls.append(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "build", spy)
    f = write_cover_geojson(tmp_path / "a.geojson", cover_polygons())
    run_kept(scene, tmp_path / "x.gpkg", *classes(f), *prefer(lake), "--shares", "--threads", "3")
    fixing = [k for k in calls if k.get("fix", True)]
    assert len(fixing) == 1 and len(calls) == 2, calls
    assert fixing[0].get("threads") == 3, calls


def test_10_the_box_first_given_keeps_its_area_and_border(
    scene: dict[str, Path], tmp_path: Path, lake: Path
) -> None:
    polygons = dict(cover_polygons())
    a, b = split(tmp_path, [(311, polygons[311]), (312, polygons[312])],
                 [(512, polygons[512]), (ARABLE, planted_box())])  # fmt: skip
    args = (*classes(b, a), *prefer(lake), *edges("all"), *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args)) as conn:
        found = lengths(conn)
        overlap = overlap_figure(conn)
    assert found["land_cover"] == pytest.approx(9884.3, abs=LENGTH), found
    assert found["water"] == pytest.approx(2600.4, abs=LENGTH), found
    assert "b.geojson" in overlap, overlap


def halves() -> list[list[tuple[int, BaseGeometry]]]:
    """311 cut at x0 + 750.5: the west half with 312, the east half with 512."""
    x0, _ = dem_origin()
    polygons = dict(cover_polygons())
    cut = x0 + 750.5
    west = polygons[311].intersection(box(cut - 1e4, -1e8, cut, 1e8))
    east = polygons[311].intersection(box(cut, -1e8, cut + 1e4, 1e8))
    return [[(311, west), (312, polygons[312])], [(311, east), (512, polygons[512])]]


def test_10_the_cut_between_two_files_premise() -> None:
    (_, west), (_, east) = halves()[0][0], halves()[1][0]
    line = west.intersection(east)
    assert line.length == pytest.approx(1699.4, abs=0.05)
    assert west.area + east.area == pytest.approx(dict(cover_polygons())[311].area, abs=AREA)


def test_10_a_county_line_without_the_merge_is_an_edge(
    scene: dict[str, Path], tmp_path: Path, lake: Path
) -> None:
    a, b = split(tmp_path, *halves())
    args = (*classes(a, b), *prefer(lake), *edges("all"), *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args, "--no-edges-merge-same-class")) as c:
        found = lengths(c)
    assert found["land_cover"] == pytest.approx(10782.1, abs=LENGTH), found


def test_10_a_county_line_with_one_class_on_both_sides_is_not_followed(
    scene: dict[str, Path], tmp_path: Path, lake: Path
) -> None:
    """Ola's question 4 at its default: the files are one map, so with the
    same-class merge the cut between them is no border. A veto changes this
    test alone."""
    a, b = split(tmp_path, *halves())
    args = (*classes(a, b), *prefer(lake), *edges("all"), *AS_READ)
    with closing(run_gpkg(scene, tmp_path / "x.gpkg", *args, "--edges-merge-same-class")) as c:
        found = lengths(c)
    assert found["land_cover"] == pytest.approx(9082.7, abs=LENGTH), found


def test_10_a_class_raster_with_another_land_cover_is_refused(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    args = ("--land-cover", KEY, "--land-cover", str(cover), "--land-cover-classes", "corine")
    flat = refused(scene, tmp_path, *args)
    assert "--land-cover" in flat, flat


def test_10_two_catalogue_entries_whose_classes_differ_are_refused(
    scene: dict[str, Path], tmp_path: Path
) -> None:
    one = write_cover_gpkg(tmp_path / "one.gpkg", cover_polygons())
    two = write_cover_gpkg(tmp_path / "two.gpkg", cover_polygons())
    config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
    write_sources(config, {
        "cover-one": vector_entry(one, "Test cover one", classes="corine"),
        "cover-two": vector_entry(two, "Test cover two", classes="clc18_kode"),
    })  # fmt: skip
    flat = refused(scene, tmp_path, "--land-cover", "cover-one", "--land-cover", "cover-two")
    assert "--land-cover" in flat, flat


# ------------------------------------------------------------ test 11: a cut run

#: Pinned from today's tree with `--features cover.geojson --features-rules
#: corine`, the clean-up switches 0 and `--pieces 4`.
PIECES = {
    "0-0-0": (447, {0: 171, 311: 276}),
    "0-1-0": (484, {0: 198, 311: 15, 312: 87, 512: 184}),
    "1-0-0": (437, {0: 159, 311: 278}),
    "1-1-0": (438, {0: 258, 311: 14, 512: 166}),
}


def pieces_of(tmp_path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Per piece file, its triangles' corners and codes."""
    (directory,) = [p for p in tmp_path.iterdir() if p.name.endswith(".pieces")]
    out = {}
    for f in sorted(directory.glob("*.vtk")):
        vtk = read_vtk(f.read_bytes())
        codes = np.asarray(vtk.cell_array("land_cover_code").values).astype(int)
        corners = np.array([np.asarray(vtk.points)[c] for c in vtk.polygons])
        out[f.stem] = (corners, codes[len(vtk.lines) :])
    return out


def cut_run(scene: dict[str, Path], tmp_path: Path, *args: str) -> dict[str, Any]:
    code, output = mesh(scene, tmp_path / "x.vtk", *args, *AS_READ, "--pieces", "4")
    assert code == 0, output
    return pieces_of(tmp_path)


def area_by_code(pieces: dict[str, Any]) -> dict[int, float]:
    total: dict[int, float] = {}
    for corners, codes in pieces.values():
        areas = shapely.area(shapely.polygons(corners[:, :, :2]))
        for c in np.unique(codes):
            total[int(c)] = total.get(int(c), 0.0) + float(areas[codes == c].sum())
    return total


def test_11_a_cut_run_follows_the_class_borders_as_today(
    scene: dict[str, Path], tmp_path: Path, cover: Path
) -> None:
    pieces = cut_run(scene, tmp_path, *classes(cover), *edges("all"))
    assert sorted(pieces) == sorted(PIECES)
    for name, (triangles, counts) in PIECES.items():
        corners, codes = pieces[name]
        assert len(corners) == triangles, name
        found, n = np.unique(codes, return_counts=True)
        assert dict(zip(found.tolist(), n.tolist(), strict=True)) == counts, name
    total = area_by_code(pieces)
    for code, area in ((311, 1529799.8), (512, 910920.4), (312, 144160.4)):
        assert total[code] == pytest.approx(area, abs=0.1), (code, total)


def test_11_a_cut_run_takes_the_preferred_map(
    scene: dict[str, Path], tmp_path: Path, cover: Path, lake: Path
) -> None:
    total = area_by_code(cut_run(scene, tmp_path, *classes(cover), *prefer(lake), *edges("all")))
    for code, area in ((311, 1409620.0), (512, 420090.0), (312, 144161.0)):
        assert total[code] == pytest.approx(area, abs=AREA), (code, total)


class Assembled:
    """Every `cli.assemble_dem` call: a DEM block decoded."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.calls = 0
        real = cli.assemble_dem

        def spy(*args: Any, **kwargs: Any) -> Any:
            self.calls += 1
            return real(*args, **kwargs)

        monkeypatch.setattr(cli, "assemble_dem", spy)


@pytest.mark.parametrize(
    ("extra", "named"),
    [((), "--land-cover"), (("--edges", "all", "--shares"), "--shares")],
    ids=["shares-by-default", "shares-asked"],
)
def test_11_a_cut_run_refuses_the_shares(
    scene: dict[str, Path],
    tmp_path: Path,
    cover: Path,
    monkeypatch: pytest.MonkeyPatch,
    extra: tuple[str, ...],
    named: str,
) -> None:
    spy = Assembled(monkeypatch)
    target = tmp_path / "x.vtk"
    code, output = mesh(scene, target, *classes(cover), *extra, *AS_READ, "--pieces", "4")
    assert code == USAGE, output
    assert named in squashed(output), output
    assert "cutrun" in squashed(output), output
    assert spy.calls == 0
    assert not target.exists()
    assert not [p for p in tmp_path.iterdir() if p.name.endswith(".pieces")]
