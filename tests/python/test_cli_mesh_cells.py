"""`rasputin mesh --land-cover ... --out x.gpkg`: the `cells` layer and the
`land_cover_fractions` table, and the legend check that guards them.

`docs/increments/105-cell-table.md`, section 4.1 (the check runs only when the
cells table is written), 4.4 (both tables whenever `--land-cover` meets a
`.gpkg`; no `cells` without a legend) and section 6, T5. The scene is
`test_cli_mesh_land_cover_fractions.py`'s (DEM 40 x 40 at 60 m in
EPSG:31983, the synthetic class cache, `--tolerance 3`). Every expected
column is recomputed from the `triangles` layer and the long-form table read
back, never from the producer's own records.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The record line for a land cover without a legend is the `run` row named
  `cells`, its value saying "no legend".
- The legend check's refusal names the class as `class <code>`.
"""

from __future__ import annotations

import re
import socket
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from class_fixtures import KEY, PATH_TEMPLATE, per_year_entry
from cli_helpers import USAGE, invoke, squashed
from fraction_fixtures import write_class_cache, write_dem, write_domain
from legend_fixtures import coded_map, write
from test_cli_fetch import install
from test_cli_land_cover_polygons import cover_polygons, write_cover_geojson
from test_cli_mesh_cache import no_network
from tin_engine import sources
from tin_engine.io.geopackage import layer_info, query_features
from tin_engine.stats import PhaseClock
from vtkread import read_vtk

#: Section 3: at most 44 classes of float32 rounding (6e-8 each) per row.
FIVE_SUM = 3e-6
#: Measured at `074e6322` (section 6, T5).
TRIANGLES, ENTRIES, AREA_M2 = 1590, 1825, 4_490_161
SCENE_CLASSES = {0, 9, 15, 21, 25, 30, 41, 75}
SHARES = ("glacier", "lake", "reservoir", "forest", "unspecified")


@pytest.fixture
def scene(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> dict[str, Path]:
    """The class source installed and the network closed; DEM, domain, cache."""
    install(
        sources, monkeypatch, per_year_entry(sources, "https://example.invalid/" + PATH_TEMPLATE)
    )
    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    return {
        "dem": write_dem(tmp_path / "dem.tif"),
        "domain": write_domain(tmp_path / "domain.geojson"),
        "cache": write_class_cache(tmp_path / "cache"),
    }


def mesh(scene: dict[str, Path], target: Path, *extra: str) -> tuple[int, str]:
    return invoke(
        "mesh", "--dem", str(scene["dem"]), "--domain", str(scene["domain"]),
        "--tolerance", "3", "--cache", str(scene["cache"]), "--out", str(target), *extra,
    )  # fmt: skip


def written(scene: dict[str, Path], target: Path, *extra: str) -> sqlite3.Connection:
    code, output = mesh(scene, target, *extra)
    assert code == 0, output
    return sqlite3.connect(target)


def column(conn: sqlite3.Connection, table: str, name: str) -> np.ndarray:
    rows = conn.execute(f"SELECT {name} FROM {table} ORDER BY fid").fetchall()
    return np.array([np.nan if r[0] is None else r[0] for r in rows], np.float64)


def shares_of(conn: sqlite3.Connection, codes: set[int]) -> np.ndarray:
    """Per triangle, the summed long-form fractions of `codes`."""
    n = conn.execute("SELECT count(*) FROM triangles").fetchone()[0]
    rows = conn.execute("SELECT triangle, code, fraction FROM land_cover_fractions").fetchall()
    out = np.zeros(n)
    for t, code, fraction in rows:
        out[t - 1] += fraction if code in codes else 0.0
    return out


def corners(conn: sqlite3.Connection) -> np.ndarray:
    found = query_features(conn, layer_info(conn, "triangles"), None, "fid")
    return np.array([np.asarray(f.geometry.exterior.coords)[:3] for f in found if f.geometry])


def check_sums(conn: sqlite3.Connection) -> None:
    five = np.column_stack([column(conn, "cells", c) for c in SHARES])
    assert (five >= 0).all()
    assert np.abs(five.sum(1) - 1).max() <= FIVE_SUM
    assert (column(conn, "cells", "savanna") <= column(conn, "cells", "unspecified")).all()


# ------------------------------------------------------------ a raster with a legend


@pytest.fixture
def mapbiomas(scene: dict[str, Path], tmp_path: Path) -> Any:
    conn = written(
        scene, tmp_path / "x.gpkg", "--land-cover", KEY, "--land-cover-legend", "mapbiomas-c11"
    )
    yield conn
    conn.close()


def test_t5_both_tables_with_the_measured_counts(mapbiomas: sqlite3.Connection) -> None:
    count = mapbiomas.execute
    assert count("SELECT count(*) FROM triangles").fetchone() == (TRIANGLES,)
    assert count("SELECT count(*) FROM cells").fetchone() == (TRIANGLES,)
    assert count("SELECT count(*) FROM land_cover_fractions").fetchone() == (ENTRIES,)
    codes = {r[0] for r in count("SELECT DISTINCT code FROM land_cover_fractions")}
    assert codes == SCENE_CLASSES
    assert {r[0] for r in count("SELECT DISTINCT catchment_id FROM cells")} == {1}


def test_t5_geometry_columns_match_the_triangles(mapbiomas: sqlite3.Connection) -> None:
    area = column(mapbiomas, "cells", "area_m2")
    assert_allclose(area, column(mapbiomas, "triangles", "area_m2"), rtol=1e-9)
    assert abs(area.sum() - AREA_M2) <= 0.5, area.sum()
    assert_allclose(column(mapbiomas, "cells", "z_m"), column(mapbiomas, "triangles", "z_mean"))
    assert_allclose(
        column(mapbiomas, "cells", "slope_deg"), column(mapbiomas, "triangles", "slope_deg"),
        rtol=1e-9, atol=1e-9,
    )  # fmt: skip
    xyz = corners(mapbiomas)
    points = query_features(mapbiomas, layer_info(mapbiomas, "cells"), None, "fid")
    got = np.array([f.geometry.coords[0] for f in points if f.geometry])
    assert_allclose(got, xyz.mean(1), rtol=1e-12)
    n = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    n *= np.sign(n[:, 2:3])
    want = np.degrees(np.arctan2(n[:, 0], n[:, 1])) % 360
    want[(n[:, 0] == 0) & (n[:, 1] == 0)] = np.nan
    gap = (column(mapbiomas, "cells", "aspect_deg") - want + 180) % 360 - 180
    assert np.nanmax(np.abs(gap)) <= 1e-9
    assert_array_equal(np.isnan(gap), np.isnan(want))


def test_t5_shares_from_the_long_form_table(mapbiomas: sqlite3.Connection) -> None:
    assert_allclose(column(mapbiomas, "cells", "forest"), shares_of(mapbiomas, {9}), atol=1e-12)
    for name in ("lake", "glacier", "savanna", "reservoir"):
        assert not column(mapbiomas, "cells", name).any(), name
    check_sums(mapbiomas)


# ------------------------------------------------------------ a raster without a legend


def test_t5_without_a_legend_only_the_fractions_table(
    scene: dict[str, Path], tmp_path: Path
) -> None:
    with closing(written(scene, tmp_path / "x.gpkg", "--land-cover", KEY)) as conn:
        tables = {r[0] for r in conn.execute("SELECT table_name FROM gpkg_contents")}
        assert "cells" not in tables and "land_cover_fractions" in tables, tables
        assert conn.execute("SELECT count(*) FROM land_cover_fractions").fetchone() == (ENTRIES,)
        (value,) = conn.execute("SELECT value FROM run WHERE name = 'cells'").fetchone()
    assert "no legend" in value, value


# ------------------------------------------------------------ the legend check


@pytest.fixture
def lake(tmp_path: Path) -> Path:
    """`coded_map`: Lake (code 2) is water and has no category."""
    return write(tmp_path / "lake.toml", coded_map())


def test_t5_uncategorised_water_refused_for_a_gpkg_before_the_dem(
    scene: dict[str, Path], tmp_path: Path, lake: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    phases: list[str] = []
    original = PhaseClock.add
    monkeypatch.setattr(
        PhaseClock, "add", lambda self, name, s: (phases.append(name), original(self, name, s))
    )
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, "--land-cover", KEY, "--land-cover-legend", str(lake))
    assert code == USAGE, output
    assert "--land-cover-legend" in output and re.search(r"class 2\b", output), output
    assert phases == [], phases
    assert not target.exists() and not (tmp_path / "x.gpkg.part").exists()


def test_t5_the_same_legend_writes_a_vtk(
    scene: dict[str, Path], tmp_path: Path, lake: Path
) -> None:
    target = tmp_path / "x.vtk"
    code, output = mesh(scene, target, "--land-cover", KEY, "--land-cover-legend", str(lake))
    assert code == 0, squashed(output)
    names = set(read_vtk(target.read_bytes()).field_data)
    assert {"land_cover_offsets", "land_cover_classes", "land_cover_fractions"} <= names


# ------------------------------------------------------------ class polygons


def test_t5_polygon_route_writes_cells(scene: dict[str, Path], tmp_path: Path) -> None:
    """Section 4.4: both routes. `cover_polygons`: forest 311 and 312, water 512."""
    cover = write_cover_geojson(tmp_path / "cover.geojson", cover_polygons())
    target = tmp_path / "p.gpkg"
    with closing(
        written(scene, target, "--land-cover", str(cover), "--land-cover-legend", "corine")
    ) as conn:
        assert_allclose(column(conn, "cells", "forest"), shares_of(conn, {311, 312}), atol=1e-12)
        assert_allclose(column(conn, "cells", "lake"), shares_of(conn, {512}), atol=1e-12)
        assert column(conn, "cells", "lake").max() > 0.5, "the water polygon is met"
        check_sums(conn)
