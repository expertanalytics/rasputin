"""`rasputin mesh --catchments UNITS`: the units' borders as constraints and
each cell's `catchment_id`, end to end.

`docs/increments/106-wednesday-cells.md`, section 4, tests 5 to 8. The scene:
a synthetic DEM, 56 x 56 nodes at 60 m in EPSG:31983 (`rough`), node (0, 0)
at `dem_origin()`; the test class raster as land cover with the
`mapbiomas-c11` legend, so a `cells` layer is written. The design's local
frame, [0, 3000] x [0, 3000] m, is placed at `ORIGIN`, 130.5 m inside the
DEM's west and south edges, so no unit corner is a DEM node. Every oracle is
GEOS's, on the `triangles` layer read back, never the producer's records.

The three units have two slanted borders, on no grid line or diagonal:
A (3,900,000 m2), B (4,165,000 m2) and C (935,000 m2) tile the square.
The straddle check can fail: on the run without `--catchments`
(`--triangles 200`), 35 of its 197 triangles have more than `SHARE` of their
area in a second unit (measured at the red step).

PINNED HERE, where the design leaves it open (listed for `@architect`):

- A file's `catchments` record entry is the file's name, as `--reservoirs`
  records a file; the record entries are read from the `.gpkg`'s `run`
  table, `catchments_units` as the text `3`.
- Test 6 compares the `triangles` and `cells` rows of the run without
  `--catchments` with master's, read from
  `data/catchments_rows_e8251ee8.json`: "the same rows as master". The file
  was written on macOS from a run whose rows hashed to the sha256 recorded
  at the design's head `e8251ee8`, before any 106c code
  (`9e9dd934ee4b...1c1ef9265e35e`); geometries are stored as their
  coordinates, x, y, z per vertex.
"""

from __future__ import annotations

import json
import math
import socket
import sqlite3
from collections.abc import Callable, Iterator
from pathlib import Path

import numpy as np
import pytest
import shapely
from numpy.testing import assert_allclose, assert_array_equal
from pyproj import Transformer
from shapely.affinity import translate
from shapely.geometry import Polygon, box
from shapely.ops import transform

from class_fixtures import KEY, PATH_TEMPLATE, per_year_entry
from cli_helpers import USAGE, geojson, invoke, squashed
from feature_fixtures import Feat, write_geojson
from fraction_fixtures import MESH_CRS, dem_origin, write_class_cache
from geographic_fixtures import projected_tiff, rough
from test_cli_fetch import install
from test_cli_mesh_cache import no_network
from test_cli_mesh_cells import column, corners
from tin_engine import sources
from tin_engine.io.geopackage import layer_info, query_features

STEP, NODES = 60.0, 56
#: The local frame's (0, 0) in `MESH_CRS`.
ORIGIN = (dem_origin()[0] + 130.5, dem_origin()[1] - NODES * STEP + STEP + 130.5)
#: A share of a triangle's area; GEOS's clip at coordinates near 8e6 m is
#: good to about 1e-9 m2, and the smallest triangle here is over 1e3 m2.
SHARE = 1e-6
FIDS = (11, 22, 33)


def local(g: Polygon) -> Polygon:
    return translate(g, *ORIGIN)


UNITS = [
    local(Polygon([(0, 0), (3000, 0), (3000, 700), (0, 1900)])),
    local(Polygon([(0, 1900), (3000, 700), (3000, 3000), (1700, 3000)])),
    local(Polygon([(0, 1900), (1700, 3000), (0, 3000)])),
]
SQUARE = local(box(0, 0, 3000, 3000))
RING = [(float(x), float(y)) for x, y in list(SQUARE.exterior.coords)[:-1]]
LAND_COVER = ("--land-cover", KEY, "--land-cover-legend", "mapbiomas-c11")


@pytest.fixture
def scene(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> dict[str, Path]:
    """The class source installed and the network closed; DEM, domain, units, cache."""
    install(
        sources, monkeypatch, per_year_entry(sources, "https://example.invalid/" + PATH_TEMPLATE)
    )
    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    x0, y0 = dem_origin()
    dem = tmp_path / "dem.tif"
    dem.write_bytes(
        projected_tiff(rough(NODES, NODES), x0=x0, y0=y0, step=STEP, epsg=31983).getvalue()
    )
    return {
        "dem": dem,
        "domain": geojson(tmp_path / "domain.geojson", RING, crs=MESH_CRS),
        "units": write_units(tmp_path / "units.geojson", UNITS),
        "cache": write_class_cache(tmp_path / "cache"),
    }  # fmt: skip


def write_units(path: Path, units: list[Polygon], crs: str = MESH_CRS) -> Path:
    pairs = zip(FIDS[: len(units)], units, strict=True)
    return write_geojson(path, [Feat(f, g, {"name": f"u{f}"}) for f, g in pairs], crs)


def mesh(scene: dict[str, Path], target: Path, *extra: str, domain: bool = True) -> tuple[int, str]:
    where = ("--domain", str(scene["domain"])) if domain else ()
    return invoke(
        "mesh", "--dem", str(scene["dem"]), *where, "--cache", str(scene["cache"]),
        "--out", str(target), *extra,
    )  # fmt: skip


def record(conn: sqlite3.Connection, name: str) -> str | None:
    row = conn.execute("SELECT value FROM run WHERE name = ?", (name,)).fetchone()
    return None if row is None else str(row[0])


class Run:
    """One finished run: the file, opened, and what it printed."""

    def __init__(self, conn: sqlite3.Connection, output: str) -> None:
        self.conn, self.output = conn, output


RunOf = Callable[..., Run]


@pytest.fixture
def run_of(scene: dict[str, Path], tmp_path: Path) -> Iterator[RunOf]:
    opened: list[sqlite3.Connection] = []

    def run(name: str, *extra: str) -> Run:
        target = tmp_path / name
        code, output = mesh(scene, target, *LAND_COVER, *extra)
        assert code == 0, output
        opened.append(sqlite3.connect(target))
        return Run(opened[-1], output)

    yield run
    for conn in opened:
        conn.close()


CATCHMENTS_TOLERANCE_0 = ("--triangles", "200", "--features-tolerance", "0")


# ------------------------------------------------------------ 5. end to end


def triangles_of(conn: sqlite3.Connection) -> np.ndarray:
    return shapely.polygons(corners(conn)[:, :, :2])


def test_the_fixture_is_a_coverage_of_the_square() -> None:
    assert [round(u.area) for u in UNITS] == [3_900_000, 4_165_000, 935_000]
    assert bool(shapely.coverage_is_valid(UNITS))
    assert shapely.coverage_union_all(UNITS).symmetric_difference(SQUARE).area < 1e-3
    inner = [(750, 1300), (2250, 1300), (425, 2450)]
    for unit, (x, y) in zip(UNITS, inner, strict=True):
        assert unit.point_on_surface().equals(local(shapely.Point(x, y))), unit.point_on_surface()


def test_no_cell_straddles_and_each_takes_its_units_id(
    run_of: RunOf, scene: dict[str, Path]
) -> None:
    ran = run_of("with.gpkg", "--catchments", str(scene["units"]), *CATCHMENTS_TOLERANCE_0)
    triangles = triangles_of(ran.conn)
    area = shapely.area(triangles)
    shares = np.column_stack(
        [shapely.area(shapely.intersection(triangles, u)) / area for u in UNITS]
    )
    second = np.sort(shares, axis=1)[:, -2]
    assert second.max() <= SHARE, f"{int((second > SHARE).sum())} cells straddle two units"
    assert abs(area.sum() - 9_000_000) <= 1e-6 * 9_000_000, area.sum()
    # First, so a rule taking a cell's first vertex, which lies on a border and so
    # within no unit, fails here as the design names it.
    assert "centroid in no unit" not in ran.output, ran.output
    centroids = shapely.centroid(triangles)
    holder = np.array([[u.contains(c) for u in UNITS] for c in centroids])
    assert (holder.sum(1) == 1).all()
    expected = np.array(FIDS)[holder.argmax(1)]
    assert_array_equal(column(ran.conn, "cells", "catchment_id"), expected)
    assert set(expected) == set(FIDS)
    assert record(ran.conn, "catchments") == "units.geojson"
    assert record(ran.conn, "catchments_units") == "3"
    assert float(record(ran.conn, "catchments_band_m") or "nan") == 0.0


# ------------------------------------------------------------ 6. unchanged without it

#: Master's rows of the run below; see the module docstring.
MASTER_ROWS = json.loads((Path(__file__).parent / "data/catchments_rows_e8251ee8.json").read_text())
#: Platforms differ in the last bits of libm's sin and cos (the DEM) and of
#: PROJ, so floats are compared within an absolute tolerance, everything else
#: exactly: ATOL in metres, square metres and shares, ATOL_DEG for `*_deg`.
#: Scales: x, y up to 8.4e6 m, z about 700 m, areas up to 1e5 m2, shares in
#: [0, 1]. One float32 ulp added to half the DEM's nodes (seeds 1 to 5) moved
#: no x, y or area, z by 6.1e-5 m, slope by 1.3e-4 deg and aspect by
#: 1.2e-3 deg; a vertex moved by 1.5 mm, or its z by 1 cm, fails.
ATOL, ATOL_DEG = 1e-3, 1e-2


def rows_of(conn: sqlite3.Connection, table: str) -> dict[str, list[object]]:
    """The table's columns in fid order; a geometry as its x, y, z per vertex."""
    names = [r[1] for r in conn.execute(f"PRAGMA table_info({table})")]
    found = query_features(conn, layer_info(conn, table), None, "fid")
    rows = conn.execute(f"SELECT * FROM {table} ORDER BY fid").fetchall()
    columns = {n: [r[i] for r in rows] for i, n in enumerate(names)}
    columns["geom"] = [shapely.get_coordinates(f.geometry, include_z=True).tolist() for f in found]
    return columns


def assert_masters(conn: sqlite3.Connection, table: str) -> None:
    got, master = rows_of(conn, table), MASTER_ROWS[table]
    assert list(got) == master["columns"]
    for i, name in enumerate(master["columns"]):
        want = [row[i] for row in master["rows"]]
        where = f"{table}.{name}"
        if name == "geom":
            assert [len(g) for g in got[name]] == [len(g) for g in want], where
            got_xyz, want_xyz = np.concatenate(got[name]), np.concatenate(want)
            assert_allclose(got_xyz, want_xyz, rtol=0, atol=ATOL, err_msg=where)
        elif any(isinstance(v, float) for v in want):
            atol = ATOL_DEG if name.endswith("_deg") else ATOL
            got_v, want_v = np.array(got[name], np.float64), np.array(want, np.float64)
            assert_allclose(got_v, want_v, rtol=0, atol=atol, equal_nan=True, err_msg=where)
        else:
            assert got[name] == want, where


def test_without_catchments_the_rows_are_masters(run_of: RunOf) -> None:
    ran = run_of("without.gpkg", "--triangles", "200")
    assert set(column(ran.conn, "cells", "catchment_id")) == {1}
    assert_masters(ran.conn, "triangles")
    assert_masters(ran.conn, "cells")
    assert record(ran.conn, "catchments") is None


# ------------------------------------------------------------ 8. the band from the budget


def test_the_band_follows_the_budget_end_to_end(run_of: RunOf, scene: dict[str, Path]) -> None:
    """A = 9e6 m2, N = 200: h = 322.371 m, band 40.296 m."""
    ran = run_of("budget.gpkg", "--catchments", str(scene["units"]), "--triangles", "200")
    h = math.sqrt(4 * 9_000_000 / (math.sqrt(3) * 200))
    assert h == pytest.approx(322.371, abs=1e-3)
    band = float(record(ran.conn, "catchments_band_m") or "nan")
    assert band == pytest.approx(40.29637244338282, rel=1e-9)
    assert "a triangle's edge at 200 triangles" in " ".join(ran.output.split()), ran.output


def test_the_given_tolerance_is_named_as_the_band(run_of: RunOf, scene: dict[str, Path]) -> None:
    ran = run_of("given.gpkg", "--catchments", str(scene["units"]), *CATCHMENTS_TOLERANCE_0)
    assert "(--features-tolerance)" in ran.output, ran.output
    assert float(record(ran.conn, "catchments_band_m") or "nan") == 0.0


# ------------------------------------------------------------ 7. refusals


def refused(
    scene: dict[str, Path], tmp_path: Path, *extra: str, says: str, domain: bool = True
) -> None:
    target = tmp_path / "refused.gpkg"
    code, output = mesh(scene, target, *extra, domain=domain)
    assert code == USAGE, output
    assert "Traceback" not in output
    flat = squashed(output)
    assert "Invalidvaluefor--catchments" in flat, output
    assert squashed(says) in flat, output
    assert not target.exists()


def test_refusal_1_needs_a_domain(scene: dict[str, Path], tmp_path: Path) -> None:
    refused(
        scene, tmp_path, *LAND_COVER, "--catchments", str(scene["units"]), *CATCHMENTS_TOLERANCE_0,
        says="needs --domain: the units whose inner point lies in the domain are the "
        "sub-catchments",
        domain=False,
    )  # fmt: skip


@pytest.mark.parametrize("extra", [(), ("--land-cover", KEY)], ids=["no-land-cover", "no-legend"])
def test_refusal_2_needs_a_cells_layer(
    scene: dict[str, Path], tmp_path: Path, extra: tuple[str, ...]
) -> None:
    refused(
        scene, tmp_path, *extra, "--catchments", str(scene["units"]), *CATCHMENTS_TOLERANCE_0,
        says="needs --land-cover with a legend and a .gpkg output",
    )  # fmt: skip


def test_refusal_2_needs_a_gpkg_output(scene: dict[str, Path], tmp_path: Path) -> None:
    target = tmp_path / "refused.vtk"
    code, output = invoke(
        "mesh", "--dem", str(scene["dem"]), "--domain", str(scene["domain"]),
        "--cache", str(scene["cache"]), "--out", str(target), *LAND_COVER,
        "--catchments", str(scene["units"]), *CATCHMENTS_TOLERANCE_0,
    )  # fmt: skip
    assert code == USAGE, output
    flat = squashed(output)
    assert "Invalidvaluefor--catchments" in flat and squashed("a .gpkg output") in flat, output
    assert not target.exists()


def test_refusal_3_units_in_degrees(scene: dict[str, Path], tmp_path: Path) -> None:
    to_wgs84 = Transformer.from_crs(MESH_CRS, "EPSG:4326", always_xy=True).transform
    units = [transform(to_wgs84, u) for u in UNITS]
    path = write_units(tmp_path / "deg.geojson", units, crs="EPSG:4326")
    refused(
        scene, tmp_path, *LAND_COVER, "--catchments", str(path), *CATCHMENTS_TOLERANCE_0,
        says="deg.geojson is in degrees (EPSG:4326); the border band is in metres, so give the "
        "units in a projected CRS",
    )  # fmt: skip


def test_refusal_4_no_inner_point_in_the_domain(scene: dict[str, Path], tmp_path: Path) -> None:
    far = write_units(tmp_path / "far.geojson", [translate(u, 100_000, 0) for u in UNITS])
    refused(
        scene, tmp_path, *LAND_COVER, "--catchments", str(far), *CATCHMENTS_TOLERANCE_0,
        says="no unit of far.geojson has its inner point inside the domain",
    )  # fmt: skip


def test_refusal_6_units_in_parts_that_do_not_touch(scene: dict[str, Path], tmp_path: Path) -> None:
    apart = [local(box(100, 100, 1100, 1100)), local(box(1900, 1900, 2900, 2900))]
    path = write_units(tmp_path / "apart.geojson", apart)
    refused(
        scene, tmp_path, *LAND_COVER, "--catchments", str(path), *CATCHMENTS_TOLERANCE_0,
        says="the units join into 2 parts that do not touch",
    )  # fmt: skip


def test_refusal_7_needs_a_budget_or_a_band(scene: dict[str, Path], tmp_path: Path) -> None:
    refused(
        scene, tmp_path, *LAND_COVER, "--catchments", str(scene["units"]), "--tolerance", "3",
        says="needs --triangles, whose budget sets the border band, or --features-tolerance",
    )  # fmt: skip
