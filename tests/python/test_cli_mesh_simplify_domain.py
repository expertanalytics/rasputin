"""`rasputin mesh --simplify-domain`: the `--domain` outline simplified to the
triangle budget, and a one-part `MultiPolygon` read as its polygon.

`docs/increments/106-wednesday-cells.md`, section 8, tests 1 to 7. The scene
is `test_cli_mesh_catchments.py`'s: 56 x 56 DEM nodes at 60 m in EPSG:31983,
the local frame placed at its `ORIGIN`, the `mapbiomas-c11` land cover, so a
`.gpkg` with a `run` table is written.

`ZZ` is the square [0, 3000]^2 with each side a zigzag of 60 steps of 50 m,
every other vertex pushed 20 m outward: 240 vertices, 9,120,000 m2. `ZZH` is
`ZZ` with a hole, the square [1200, 1800]^2 zigzagged the same way, 12 steps
a side, every other vertex pushed 20 m away from the hole's centre: 288
vertices, 8,736,000 m2. Every zigzag is a whole number of 50 m steps, so the
areas are exact in binary.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The pure simplification is `sub_catchments.simplify_domain(domain,
  triangles, tolerance)`, returning the simplified `DomainPolygon`, in the
  domain's CRS, and the band in metres; `tolerance` wins when given, as
  `_units` reads `--features-tolerance`.
- The stderr line counts a ring's vertices without its closing repeat (ZZ
  is "240", as the design's line says).
- Test 5's digest of the `triangles` and `cells` rows is recorded on 106e's
  design head `ce0680ee`, before any 106e code: "the same rows as 106c's head".
"""

from __future__ import annotations

import hashlib
import json
import math
import sqlite3
from pathlib import Path

import pytest
from pyproj import Geod, Transformer
from shapely.affinity import translate
from shapely.geometry import MultiPolygon, Polygon, mapping
from shapely.geometry.polygon import orient
from shapely.ops import transform

from cli_helpers import USAGE, geojson, invoke, squashed
from fraction_fixtures import MESH_CRS
from test_cli_mesh_catchments import (  # noqa: F401  (fixtures, by name)
    CATCHMENTS_TOLERANCE_0,
    ORIGIN,
    Run,
    RunOf,
    mesh,
    record,
    run_of,
    scene,
)
from tin_engine import sub_catchments
from tin_engine.crs import parse_crs, same_crs
from tin_engine.domain import DomainError, DomainPolygon
from tin_engine.io.domain_file import read_domain

Ring = list[tuple[float, float]]

#: The outward normal of each side of a counter-clockwise square, south first.
NORMALS = ((0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0))


def zigzag(x0: float, x1: float, steps: int, push: float = 20.0) -> Ring:
    """The square [x0, x1]^2, counter-clockwise from (x0, x0), each side cut
    into `steps` steps and every other vertex pushed `push` along the side's
    outward normal; the corners stay."""
    corners = [(x0, x0), (x1, x0), (x1, x1), (x0, x1)]
    ring = []
    for k, (nx, ny) in enumerate(NORMALS):
        (ax, ay), (bx, by) = corners[k], corners[(k + 1) % 4]
        for i in range(steps):
            t, d = i / steps, push if i % 2 else 0.0
            ring.append((ax + t * (bx - ax) + d * nx, ay + t * (by - ay) + d * ny))
    return ring


ZZ = translate(Polygon(zigzag(0.0, 3000.0, 60)), *ORIGIN)
ZZH = translate(Polygon(zigzag(0.0, 3000.0, 60), [zigzag(1200.0, 1800.0, 12)]), *ORIGIN)
#: band_m(9,120,000, 100): sqrt(4 A / (sqrt(3) N)) / 8.
ZZ_BAND_100 = 57.36633624679954


def vertices(p: Polygon) -> tuple[int, ...]:
    """Each ring's vertex count, the closing repeat left out: outer first."""
    return tuple(len(r.coords) - 1 for r in (p.exterior, *p.interiors))


def ring_of(p: Polygon) -> Ring:
    return [(float(x), float(y)) for x, y in list(p.exterior.coords)[:-1]]


def test_the_fixtures_are_as_the_design_states() -> None:
    assert ZZ.is_valid and ZZH.is_valid
    assert vertices(ZZ) == (240,) and sum(vertices(ZZH)) == 288
    assert ZZ.area == pytest.approx(9_120_000, abs=1e-6)
    assert ZZH.area == pytest.approx(8_736_000, abs=1e-6)
    assert sub_catchments.band_m(ZZ.area, 100) == pytest.approx(ZZ_BAND_100, rel=1e-9)


# ------------------------------------------------------------ 1. read


def write_doc(path: Path, geometry: object) -> Path:
    doc = {
        "type": "Feature",
        "properties": {},
        "geometry": mapping(geometry),
        "crs": {"type": "name", "properties": {"name": MESH_CRS}},
    }
    path.write_text(json.dumps(doc))
    return path


def test_1_a_one_part_multipolygon_reads_as_its_polygon_oriented(tmp_path: Path) -> None:
    clockwise = orient(ZZ, sign=-1.0)
    read = read_domain(write_doc(tmp_path / "one.geojson", MultiPolygon([clockwise])))
    assert isinstance(read.polygon, Polygon)
    assert read.polygon.equals_exact(orient(ZZ, sign=1.0), tolerance=0)
    assert read.polygon.exterior.is_ccw


def test_1_a_two_part_multipolygon_is_refused(tmp_path: Path) -> None:
    two = MultiPolygon([ZZ, translate(ZZ, 10_000.0, 0.0)])
    with pytest.raises(DomainError) as caught:
        read_domain(write_doc(tmp_path / "two.geojson", two))
    assert "the domain must be one Polygon, got a MultiPolygon of 2 parts" in str(caught.value)


def test_1_a_polygon_reads_as_today(tmp_path: Path) -> None:
    read = read_domain(write_doc(tmp_path / "poly.geojson", orient(ZZ, sign=-1.0)))
    assert read.polygon.equals_exact(orient(ZZ, sign=1.0), tolerance=0)
    assert same_crs(parse_crs(read.crs), parse_crs(MESH_CRS))


# ------------------------------------------------------------ 2. simplify, projected


def test_2_projected_zz_at_100_triangles() -> None:
    given = DomainPolygon(polygon=ZZ, crs=MESH_CRS)
    out, band = sub_catchments.simplify_domain(given, 100, None)
    assert band == pytest.approx(ZZ_BAND_100, rel=1e-9)
    assert vertices(out.polygon) == (4,)
    assert abs(out.polygon.area / ZZ.area - 1) < 1e-12
    assert out.polygon.is_valid
    assert same_crs(parse_crs(out.crs), parse_crs(MESH_CRS))


def test_2_the_hole_is_simplified_with_the_outline() -> None:
    """Without the frame the hole keeps its 48 vertices and the outer ring its 240."""
    out, band = sub_catchments.simplify_domain(DomainPolygon(polygon=ZZH, crs=MESH_CRS), 100, None)
    assert band == pytest.approx(56.14563624873581, rel=1e-9)
    assert vertices(out.polygon) == (4, 4)
    assert abs(out.polygon.area / ZZH.area - 1) < 1e-12
    assert out.polygon.is_valid


def test_2_the_given_tolerance_is_the_band() -> None:
    out, band = sub_catchments.simplify_domain(DomainPolygon(polygon=ZZ, crs=MESH_CRS), 100, 30.0)
    assert band == 30.0
    assert abs(out.polygon.area / ZZ.area - 1) < 1e-12


# ------------------------------------------------------------ 3. simplify, in degrees

GEOD = Geod(ellps="WGS84")


def geodesic_area(p: Polygon) -> float:
    return abs(GEOD.geometry_area_perimeter(p)[0])


def test_3_in_degrees_the_band_is_in_metres_on_an_equal_area_projection() -> None:
    """The band from the area in square degrees (about 5e-4 m) fails here."""
    to_wgs84 = Transformer.from_crs(MESH_CRS, "EPSG:4326", always_xy=True).transform
    degrees = transform(to_wgs84, ZZ)
    assert degrees.centroid.x == pytest.approx(-44.93, abs=0.01)
    assert degrees.centroid.y == pytest.approx(-15.24, abs=0.01)
    out, band = sub_catchments.simplify_domain(
        DomainPolygon(polygon=degrees, crs="EPSG:4326"), 100, None
    )
    assert same_crs(parse_crs(out.crs), parse_crs("EPSG:4326"))
    assert vertices(out.polygon) == (4,)
    assert abs(geodesic_area(out.polygon) / geodesic_area(degrees) - 1) < 1e-6
    assert band == pytest.approx(57.4, rel=1e-3)


# ------------------------------------------------------------ 4. end to end


def with_domain(given: dict[str, Path], polygon: Polygon, name: str) -> dict[str, Path]:
    return {**given, "domain": geojson(given["dem"].parent / name, ring_of(polygon), crs=MESH_CRS)}


@pytest.fixture
def zz_scene(scene: dict[str, Path]) -> dict[str, Path]:  # noqa: F811
    """The scene with `ZZ` as `--domain`; `run_of` reads `scene` from it."""
    scene.update(with_domain(scene, ZZ, "zz.geojson"))
    return scene


def test_4_without_the_flag_the_refusal_names_it(zz_scene: dict[str, Path], tmp_path: Path) -> None:
    target = tmp_path / "refused.gpkg"
    code, output = mesh(zz_scene, target, "--triangles", "100")
    assert code == USAGE, output
    flat = squashed(output)
    assert squashed("the fewest triangles this DEM gives is 470") in flat, output
    assert squashed("; or give --simplify-domain to simplify the outline to the budget") in flat
    assert not target.exists()


def test_4_with_the_flag_the_budget_is_met(zz_scene: dict[str, Path], run_of: RunOf) -> None:  # noqa: F811
    ran = run_of("flag.gpkg", "--triangles", "100", "--simplify-domain")
    count = ran.conn.execute("SELECT count(*) FROM triangles").fetchone()[0]
    assert count == 101, count
    assert abs(count - 100) <= 5
    line = (
        "domain: outline simplified at 57 m (a triangle's edge at 100 triangles, 459 m, "
        "over 8): 240 to 4 vertices, 9.12 km2"
    )
    assert line in " ".join(ran.output.split()), ran.output
    assert float(record(ran.conn, "domain_band_m") or "nan") == pytest.approx(ZZ_BAND_100, rel=1e-9)


def test_4_the_given_tolerance_sets_the_band(zz_scene: dict[str, Path], run_of: RunOf) -> None:  # noqa: F811
    ran = run_of(
        "given.gpkg", "--triangles", "100", "--simplify-domain", "--features-tolerance", "30"
    )
    assert float(record(ran.conn, "domain_band_m") or "nan") == 30.0
    said = " ".join(ran.output.split())
    assert "domain: outline simplified at 30 m (--features-tolerance): 240 to" in said, said


# ------------------------------------------------------------ 5. unchanged without it

#: sha256 over the `triangles` and `cells` rows, in fid order, of the run
#: below, recorded at `ce0680ee` before any 106e code.
ROWS_DIGEST_400 = "16269330e28c0a2346c67b464bdd759c6eb513a705d6393c4d45565ed2b39ed6"


def rows_digest(conn: sqlite3.Connection) -> str:
    digest = hashlib.sha256()
    for table in ("triangles", "cells"):
        for row in conn.execute(f"SELECT * FROM {table} ORDER BY fid"):
            digest.update(repr(row).encode())
    return digest.hexdigest()


def test_5_a_hand_simplified_domain_meshes_as_today(
    scene: dict[str, Path],  # noqa: F811
    run_of: RunOf,  # noqa: F811
) -> None:
    """`ZZ` simplified as the flag would at 400 triangles, given as the domain."""
    by_hand = sub_catchments.simplify([ZZ], sub_catchments.band_m(ZZ.area, 400))[0]
    assert vertices(by_hand) == (4,)
    scene.update(with_domain(scene, by_hand, "by-hand.geojson"))
    ran = run_of("by-hand.gpkg", "--triangles", "400")
    assert ran.conn.execute("SELECT count(*) FROM triangles").fetchone()[0] == 396
    assert rows_digest(ran.conn) == ROWS_DIGEST_400
    assert record(ran.conn, "domain_band_m") is None
    assert "domain: outline simplified" not in ran.output


# ------------------------------------------------------------ 6. with --catchments


def test_6_with_catchments_the_flag_changes_nothing(
    scene: dict[str, Path],  # noqa: F811
    run_of: RunOf,  # noqa: F811
) -> None:
    units = ("--catchments", str(scene["units"]), *CATCHMENTS_TOLERANCE_0)
    without: Run = run_of("without.gpkg", *units)
    flagged: Run = run_of("flagged.gpkg", *units, "--simplify-domain")
    assert rows_digest(flagged.conn) == rows_digest(without.conn)


# ------------------------------------------------------------ 7. refusals


def refused(given: dict[str, Path], tmp_path: Path, *extra: str, says: str) -> None:
    target = tmp_path / "refused.vtk"
    code, output = invoke(
        "mesh", "--dem", str(given["dem"]), "--out", str(target), "--simplify-domain", *extra
    )
    assert code == USAGE, output
    assert "Traceback" not in output
    flat = squashed(output)
    assert "Invalidvaluefor--simplify-domain" in flat, output
    assert squashed(says) in flat, output
    assert not target.exists()


def test_7_needs_a_domain(scene: dict[str, Path], tmp_path: Path) -> None:  # noqa: F811
    refused(scene, tmp_path, "--triangles", "100", says="needs --domain")


def test_7_needs_a_budget_or_a_band(zz_scene: dict[str, Path], tmp_path: Path) -> None:
    refused(
        zz_scene, tmp_path, "--domain", str(zz_scene["domain"]), "--tolerance", "3",
        says="needs --triangles, whose budget sets the outline band, or --features-tolerance",
    )  # fmt: skip


def test_the_band_is_an_eighth_of_the_edge() -> None:
    """The figures the stderr line of test 4 states: h 459 m at N 100."""
    h = math.sqrt(4 * ZZ.area / (math.sqrt(3) * 100))
    assert round(h) == 459 and round(h / 8) == 57
