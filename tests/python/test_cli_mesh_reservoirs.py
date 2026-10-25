"""`rasputin mesh --land-cover ... --reservoirs POLYGONS --out x.gpkg`: each
cell's lake share split into lake and reservoir by polygon overlay.

`docs/increments/105-cell-table.md`, section 4.2 (`r_t` the share of the
triangle under the union of the polygons, moved to the mesh CRS;
`reservoir = min(lake, r)`, `lake -= reservoir`; refused without
`--land-cover` or with a `.vtk` output; without the option the run record
says so) and section 6, T7's CLI case. The scene is `test_cli_mesh_cells.py`'s,
on the polygon route with `corine` (`cover_polygons`: forest 311 and 312,
water 512), since the raster scene has no lake class. The reservoir rectangle
covers part of the forest polygon, the uncovered land between, and part of
the water polygon, so `min(lake, r)` takes both branches.

THE ORACLE: `r` per triangle from the `triangles` layer read back, by GEOS's
`intersection(...).area` with the rectangle in the mesh CRS (moved vertex by
vertex with pyproj when the file is in WGS 84, as the producer's edges are
straight in the mesh CRS too). Bound 1e-6 on a share: the cut at these
coordinates (below 1e7 m, triangles of 100 m and up) loses under 1e-9 of a
triangle, and the lake share it is compared with is float32 (6e-8).

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The run record row is named `reservoirs`; with the option, its value names
  the file; without it (on a run that writes the `cells` layer), it says
  `none`. A run without `--land-cover` has no such row.
- Both refusals name `--land-cover` and `.gpkg`, exit 2, before any phase is
  timed. A reservoir file of a suffix rasputin does not read (`.shp`) is a
  usage error naming `--reservoirs` and the suffix.
"""

from __future__ import annotations

import re
import sqlite3
from collections.abc import Iterator
from contextlib import closing
from pathlib import Path

import numpy as np
import pytest
import shapely
from numpy.testing import assert_allclose, assert_array_equal
from shapely.geometry import Polygon, box

from class_fixtures import KEY
from cli_helpers import USAGE, squashed
from feature_fixtures import Feat, write_geojson
from fraction_fixtures import MESH_CRS, dem_origin
from test_cli_land_cover_polygons import cover_polygons, to_mesh, to_wgs84, write_cover_geojson
from test_cli_mesh_cells import check_sums, column, corners, mesh, scene, written
from tin_engine.stats import PhaseClock

__all__ = ["scene"]  # a fixture, used by name

#: On a share: see the module docstring.
SHARE = 1e-6


def reservoir_rectangle() -> Polygon:
    """Over the east end of forest 311, the land beyond it and the west of water
    512; no side on a DEM node (the DEM's step is 60 m)."""
    x0, y0 = dem_origin()
    return box(x0 + 900.2, y0 - 1500.3, x0 + 1900.7, y0 - 900.6)


def write_reservoirs(path: Path, wgs84: bool) -> Path:
    shape = to_wgs84(reservoir_rectangle()) if wgs84 else reservoir_rectangle()
    feats = [Feat(1, shape, {"magasinNr": "1"})]
    return write_geojson(path, feats, crs=None if wgs84 else MESH_CRS)


def record(conn: sqlite3.Connection, name: str) -> str | None:
    row = conn.execute("SELECT value FROM run WHERE name = ?", (name,)).fetchone()
    return None if row is None else str(row[0])


def oracle_r(conn: sqlite3.Connection, reservoir: Polygon) -> np.ndarray:
    triangles = shapely.polygons(corners(conn)[:, :, :2])
    return shapely.area(shapely.intersection(triangles, reservoir)) / shapely.area(triangles)


@pytest.fixture
def cover(tmp_path: Path) -> Path:
    return write_cover_geojson(tmp_path / "cover.geojson", cover_polygons())


@pytest.fixture
def without(scene: dict[str, Path], tmp_path: Path, cover: Path) -> Iterator[sqlite3.Connection]:
    """The polygon route with `corine`, no `--reservoirs`."""
    target = tmp_path / "without.gpkg"
    args = ("--land-cover", str(cover), "--land-cover-legend", "corine")
    with closing(written(scene, target, *args)) as conn:
        yield conn


def test_t7_without_reservoirs_the_record_says_none(without: sqlite3.Connection) -> None:
    value = record(without, "reservoirs")
    assert value is not None and re.search(r"\bnone\b", value), value
    assert not column(without, "cells", "reservoir").any()


@pytest.mark.parametrize("wgs84", [False, True], ids=["mesh-crs", "wgs84"])
def test_t7_reservoir_is_min_of_lake_and_r(
    scene: dict[str, Path], tmp_path: Path, cover: Path, without: sqlite3.Connection, wgs84: bool
) -> None:
    res = write_reservoirs(tmp_path / "res.geojson", wgs84)
    args = ("--land-cover", str(cover), "--land-cover-legend", "corine", "--reservoirs", str(res))
    with closing(written(scene, tmp_path / "with.gpkg", *args)) as conn:
        assert_array_equal(corners(conn), corners(without)), "the reservoirs moved the mesh"
        lake_before = column(without, "cells", "lake")
        lake, reservoir = column(conn, "cells", "lake"), column(conn, "cells", "reservoir")
        rectangle = reservoir_rectangle()
        r = oracle_r(conn, to_mesh(to_wgs84(rectangle)) if wgs84 else rectangle)
        assert_allclose(reservoir, np.minimum(lake_before, r), rtol=0, atol=SHARE)
        assert_allclose(lake + reservoir, lake_before, rtol=0, atol=1e-12)
        assert (reservoir <= lake_before).all()
        # both branches of the min are met: the forest under the rectangle
        # adds nothing, and the water under it is (mostly) regulated
        assert ((r > 0.5) & (lake_before == 0)).any(), "no land under the rectangle"
        assert ((lake_before > 0.5) & (r > 0.5)).any(), "no water under the rectangle"
        assert reservoir.max() > 0.5, reservoir.max()
        for name in ("forest", "glacier", "unspecified", "savanna"):
            assert_array_equal(column(conn, "cells", name), column(without, "cells", name), name)
        check_sums(conn)
        value = record(conn, "reservoirs")
    assert value is not None and "res.geojson" in value, value


# ------------------------------------------------------------ refusals


@pytest.fixture
def phases(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    timed: list[str] = []
    original = PhaseClock.add
    monkeypatch.setattr(
        PhaseClock, "add", lambda self, name, s: (timed.append(name), original(self, name, s))
    )
    return timed


@pytest.mark.parametrize(
    ("suffix", "land_cover"), [(".gpkg", False), (".vtk", True)], ids=["no-land-cover", "vtk"]
)
def test_t7_reservoirs_need_land_cover_and_a_gpkg(
    scene: dict[str, Path], tmp_path: Path, phases: list[str], suffix: str, land_cover: bool
) -> None:
    res = write_reservoirs(tmp_path / "res.geojson", wgs84=False)
    target = tmp_path / f"x{suffix}"
    extra = ("--land-cover", KEY, "--land-cover-legend", "mapbiomas-c11") if land_cover else ()
    code, output = mesh(scene, target, *extra, "--reservoirs", str(res))
    assert code == USAGE, output
    flat = squashed(output)
    assert "--reservoirs" in flat and "--land-cover" in flat and ".gpkg" in flat, output
    assert phases == [], phases
    assert not target.exists()


def test_t7_a_shapefile_is_refused_naming_the_option(
    scene: dict[str, Path], tmp_path: Path
) -> None:
    res = tmp_path / "res.shp"
    res.write_bytes(b"\x00\x00\x27\x0a")
    target = tmp_path / "x.gpkg"
    code, output = mesh(
        scene, target, "--land-cover", KEY, "--land-cover-legend", "mapbiomas-c11",
        "--reservoirs", str(res),
    )  # fmt: skip
    assert code == USAGE, output
    flat = squashed(output)
    assert "--reservoirs" in flat and ".shp" in flat, output
    assert not target.exists()


def test_t7_without_land_cover_no_reservoirs_row(scene: dict[str, Path], tmp_path: Path) -> None:
    with closing(written(scene, tmp_path / "plain.gpkg")) as conn:
        assert record(conn, "reservoirs") is None
