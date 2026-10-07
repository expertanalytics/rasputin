"""Growing the catchment outline, made faster (increment 30d).

`docs/increments/30d-outline-buffer-speed.md`, sections 3-5. The seven tests
of its section 5:

1. the gate table (`grow.uses_pieces`);
2. equality with GEOS's mitred buffer on the gated staircase;
3. containment, for both routes, on the staircase and on a smooth outline
   with a sharp notch;
4. off the gate, `grow_mitred` is GEOS's buffer, WKB for WKB;
5. one buffer per reprojected decode (3.1); the tuple `target_grid_for`
   returns is in `test_target_grid.py`;
6. the features region from the convex hull (3.2);
7. byte-identical meshes, master's behaviour against the built one, run in
   the same test (no stored hash, so no platform dependence).

"Staircase" here is a raster-traced outline: every edge one step long and
axis-parallel, a vertex at every step, collinear runs included, as SMHI's
outlines are. "The gate" is `uses_pieces`: true only for a hole-free outline
of at least `PIECES_FROM` vertices whose median edge is at most half the
distance. "Pieces" is the gated route of section 3.3 (the ring cut into
pieces of `PIECE_EDGES` edges, each grown, united with the polygon).

PINNED HERE, where the design is silent:
- `grow_mitred` decides its route by calling `uses_pieces` through the
  module's global, and `uses_pieces` reads `PIECES_FROM` there at call time,
  so a test can force either route with `monkeypatch`.
- The gated route unites with `shapely.union_all`, called through the
  `shapely` module (the design's section 3.3 names it).
- Rule 3's "at most d/2" is inclusive: a median of exactly d/2 is gated.
- `feature_input.source_region` buffers nothing but a convex polygon (the
  domain's hull; section 3.2 replaces the domain's own buffer).

HOW THIS FILE GOES RED: `tin_engine.grow` does not exist. It is imported in a
fixture, so each test that needs it fails on its own with
`ModuleNotFoundError`. Test 6 needs no new name and fails on its assertion
that only a convex polygon is buffered (master buffers the staircase itself).
"""

from __future__ import annotations

import importlib
import math
import sys
from collections.abc import Callable, Iterator
from contextlib import AbstractContextManager, contextmanager
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import shapely
from pyproj import CRS
from shapely.geometry import Polygon

from cli_driver import Ring, geojson, invoke
from geotiff_fixtures import KARTVERKET, needs_codecs
from gpkg_fixtures import EXTRACT
from tin_engine.dem_input import DemRequest, open_dem
from tin_engine.domain import DomainPolygon
from tin_engine.io.domain_file import read_domain

#: Lagan's distance (section 3.3: d = 43.84 m at 10 m steps), and the
#: projected path's cell diagonal at 10 m, under two steps.
LAGAN_D = 43.84
UNDER_TWO_STEPS = 10.0 * math.sqrt(2)
#: Section 4. Containment slack: GEOS's own buffer falls 8e-7 m short of d on
#: Numedalslågen, so 1e-6 is too tight. Equality: Hausdorff and bounds.
#: Scale: coordinates to 7e6 m (UTM northings), d 14-100 m; the largest input
#: here is the 6,000-gon of 30 km radius below.
SLACK = 1e-5
SAME = 1e-6
UTM33 = "EPSG:25833"
VELHAS = Path(__file__).resolve().parents[1] / "fixtures" / "velhas"
VELHAS_CRS = "EPSG:31983"
VELHAS_H = 30  # the velhas target grid's spacing (its NOTICE: default_spacing 30 m)


# ------------------------------------------------------------------ outlines


def _walk(moves: list[tuple[int, int, int]], closed: bool = True) -> list[tuple[int, int]]:
    """Unit-step moves `(dx, dy, count)` from (0, 0), one vertex per step.
    A `closed` walk must end at (0, 0), and that closing vertex is dropped."""
    x = y = 0
    ring = [(0, 0)]
    for dx, dy, count in moves:
        for _ in range(count):
            x, y = x + dx, y + dy
            ring.append((x, y))
    if not closed:
        return ring
    assert ring[-1] == (0, 0), ring[-1]
    return ring[:-1]


def _placed(units: list[tuple[int, int]], step: float, at: tuple[float, float]) -> Polygon:
    """Integer units scaled by `step` and moved to `at`: exact floats when
    `at` and `step` are whole metres."""
    xy = np.asarray(units, dtype=np.float64) * step + np.asarray(at)
    polygon = Polygon(xy)
    assert polygon.is_valid and polygon.exterior.is_ccw
    return polygon


def staircase_disc(radius: float, step: float, at: tuple[float, float]) -> Polygon:
    """A raster-traced disc: columns one step wide, each as tall as the disc
    at its centre, the outline walked one step per vertex (about 8 r / step
    vertices). Counter-clockwise, every edge `step` long."""
    columns = 2 * round(radius / step)
    r = radius / step
    h = [
        max(1, math.floor(math.sqrt(r * r - (k + 0.5 - columns / 2) ** 2))) for k in range(columns)
    ]
    moves = [(0, -1, 2 * h[0])]  # down the west side
    for k in range(columns):  # the south side, west to east
        moves.append((1, 0, 1))
        if k + 1 < columns:
            moves.append((0, int(np.sign(h[k] - h[k + 1])), abs(h[k] - h[k + 1])))
    moves.append((0, 1, 2 * h[-1]))  # up the east side
    for k in reversed(range(columns)):  # the north side, east to west
        moves.append((-1, 0, 1))
        if k > 0:
            moves.append((0, int(np.sign(h[k - 1] - h[k])), abs(h[k - 1] - h[k])))
    units = [(x - columns // 2, y + h[0]) for x, y in _walk(moves)]
    return _placed(units, step, at)


def sawtooth(vertices: int, step: float, at: tuple[float, float]) -> Polygon:
    """A staircase of exactly `vertices` vertices: a strip whose south side
    is teeth one step high, every edge `step` long but the west side, one
    edge (a lattice walk of unit steps has an even count, so one edge is
    long to make any count). Its median edge is `step`."""
    k = 2 + (vertices - 5) % 6  # the east side's steps; the west side is one edge
    teeth = (vertices - 3 - k) // 6
    moves = [(1, 0, 1), (0, 1, 1), (1, 0, 1), (0, -1, 1)] * teeth
    moves += [(1, 0, 1), (0, 1, k), (-1, 0, 2 * teeth + 1)]
    units = _walk(moves, closed=False)  # ends at (0, k): the west side is one edge
    polygon = _placed(units, step, at)
    assert len(polygon.exterior.coords) - 1 == vertices
    return polygon


def circle(radius: float, vertices: int, at: tuple[float, float]) -> list[tuple[float, float]]:
    t = 2 * np.pi * np.arange(vertices) / vertices
    return [(at[0] + radius * math.cos(a), at[1] + radius * math.sin(a)) for a in t]


def notched(radius: float, vertices: int, at: tuple[float, float]) -> Polygon:
    """A smooth outline (a circle, edges about 2 pi r / n) with one sharp
    notch, as Numedalslågen's vertices 7537-7538 (section 3.3): at vertex A
    the outline turns by -162 degrees, runs 60 m back into the interior, turns
    by +149 degrees, runs 60 m on, and rejoins the circle two vertices later."""
    ring = circle(radius, vertices, at)
    a = np.asarray(ring[0])
    tangent = a - np.asarray(ring[-1])  # the edge arriving at A
    heading = math.atan2(tangent[1], tangent[0])
    b = a + 60.0 * np.array(
        [math.cos(heading - math.radians(162)), math.sin(heading - math.radians(162))]
    )
    turned = heading - math.radians(162) + math.radians(149)
    c = b + 60.0 * np.array([math.cos(turned), math.sin(turned)])
    polygon = Polygon([tuple(a), tuple(b), tuple(c), *ring[2:]])
    assert polygon.is_valid and polygon.exterior.is_ccw
    return polygon


def turns(polygon: Polygon) -> np.ndarray:
    """The signed turn at each vertex, in degrees (left positive)."""
    xy = np.asarray(polygon.exterior.coords)[:-1]
    d = np.roll(xy, -1, axis=0) - xy
    heading = np.degrees(np.arctan2(d[:, 1], d[:, 0]))
    return (heading - np.roll(heading, 1) + 180.0) % 360.0 - 180.0


def mitred(polygon: Polygon, distance: float) -> Polygon:
    """GEOS's own buffer, as master grows the outline."""
    grown = polygon.buffer(distance, join_style="mitre")
    assert isinstance(grown, Polygon)
    return grown


def master_source_region(domain: DomainPolygon, dem_crs: Any, source_crs: Any) -> Polygon:
    """`feature_input.source_region` as master has it (`feature_input.py:158-167`
    at `f81b20b7`): the domain itself buffered by 100 m."""
    from tin_engine.crs import reprojector, same_crs
    from tin_engine.feature_input import DENSIFY, MARGIN

    ring = shapely.segmentize(domain.polygon.buffer(MARGIN).exterior, DENSIFY)
    xy = shapely.get_coordinates(ring)
    if not same_crs(source_crs, dem_crs):
        xy = reprojector(dem_crs, source_crs)(xy)
    hull = shapely.convex_hull(shapely.multipoints(xy))
    assert isinstance(hull, Polygon)
    return hull


# ------------------------------------------------------------------ fixtures

#: Somewhere in UTM 33 north, whole metres so the staircase is exact.
NORTH = (400_000.0, 6_600_000.0)


@pytest.fixture
def grow() -> ModuleType:
    return importlib.import_module("tin_engine.grow")


@pytest.fixture(scope="module")
def staircase() -> Polygon:
    """Section 5's gated staircase: radius about 6.5 km, 10 m steps."""
    polygon = staircase_disc(6_500.0, 10.0, NORTH)
    assert len(polygon.exterior.coords) - 1 >= 5_000
    return polygon


@pytest.fixture(scope="module")
def notch() -> Polygon:
    """Smooth, 6,001 vertices, edges about 50 m (Numedalslågen's median is
    48.1 m), with the sharp notch; grown by 14.14 m as Numedalslågen is."""
    polygon = notched(47_750.0, 6_000, NORTH)
    assert len(polygon.exterior.coords) - 1 >= 5_000
    t = turns(polygon)
    assert t[0] == pytest.approx(-162.0) and t[1] == pytest.approx(149.0)
    return polygon


@pytest.fixture
def spied(monkeypatch: pytest.MonkeyPatch) -> list[tuple[Any, Any]]:
    """Every `shapely.buffer` call (which `Polygon.buffer` also goes through):
    the geometry and the distance."""
    calls: list[tuple[Any, Any]] = []
    real = shapely.buffer

    def buffer(geometry: Any, distance: Any, *args: Any, **kwargs: Any) -> Any:
        calls.append((geometry, distance))
        return real(geometry, distance, *args, **kwargs)

    monkeypatch.setattr(shapely, "buffer", buffer)
    return calls


def boundary_gap(a: Polygon, b: Polygon) -> float:
    return float(shapely.hausdorff_distance(a.boundary, b.boundary))


# ------------------------------------------------------------------ 1. the gate


@pytest.mark.parametrize(
    ("outline", "distance", "gated"),
    [
        pytest.param("staircase", LAGAN_D, True, id="staircase, d spans 4 steps"),
        pytest.param("staircase", UNDER_TWO_STEPS, False, id="staircase, d under 2 steps"),
        pytest.param("staircase", 20.0, True, id="staircase, median exactly d/2"),
        pytest.param("staircase", 19.99, False, id="staircase, median just over d/2"),
        pytest.param("sawtooth 5000", LAGAN_D, True, id="5,000 vertices"),
        pytest.param("sawtooth 4999", LAGAN_D, False, id="4,999 vertices"),
        pytest.param("6000-gon", LAGAN_D, False, id="6,000-gon, edges over d/2"),
        pytest.param("staircase with a hole", LAGAN_D, False, id="staircase with a hole"),
    ],
)
def test_1_the_gate(
    grow: ModuleType, staircase: Polygon, outline: str, distance: float, gated: bool
) -> None:
    """Section 3.3's three rules. The vertex count leaves the closing vertex
    out: the 5,000-vertex sawtooth has 5,001 coordinates."""
    outlines: dict[str, Callable[[], Polygon]] = {
        "staircase": lambda: staircase,
        "sawtooth 5000": lambda: sawtooth(5_000, 10.0, NORTH),
        "sawtooth 4999": lambda: sawtooth(4_999, 10.0, NORTH),
        # 30 km radius: edges 31.4 m, over LAGAN_D / 2 = 21.9 m.
        "6000-gon": lambda: Polygon(circle(30_000.0, 6_000, NORTH)),
        "staircase with a hole": lambda: Polygon(
            staircase.exterior.coords,
            [
                shapely.box(
                    NORTH[0] - 500, NORTH[1] - 500, NORTH[0] + 500, NORTH[1] + 500
                ).exterior.coords[::-1]
            ],
        ),
    }
    polygon = outlines[outline]()
    assert polygon.is_valid
    assert grow.uses_pieces(polygon, distance) is gated
    assert (grow.PIECES_FROM, grow.PIECE_EDGES) == (5_000, 1_000)


# ------------------------------------------------------------------ 2. equality


def test_2_the_gated_staircase_grows_as_geos_does(
    grow: ModuleType, staircase: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Section 4's second bullet: boundary Hausdorff distance to GEOS's mitred
    buffer and every bound within 1e-6 m. And it is the pieced route that
    produced it: one `union_all`, of the polygon and more than one piece."""
    assert grow.uses_pieces(staircase, LAGAN_D)
    unions: list[int] = []
    real = shapely.union_all

    def union_all(geometries: Any, *args: Any, **kwargs: Any) -> Any:
        unions.append(len(geometries))
        return real(geometries, *args, **kwargs)

    monkeypatch.setattr(shapely, "union_all", union_all)
    grown = grow.grow_mitred(staircase, LAGAN_D)
    monkeypatch.undo()
    assert len(unions) == 1 and unions[0] > 2, unions
    expected = mitred(staircase, LAGAN_D)
    assert isinstance(grown, Polygon)
    assert boundary_gap(grown, expected) <= SAME
    assert np.max(np.abs(np.subtract(grown.bounds, expected.bounds))) <= SAME


# ------------------------------------------------------------------ 3. containment


@pytest.mark.parametrize("pieces", [True, False], ids=["pieces", "geos"])
@pytest.mark.parametrize(
    ("outline", "distance"),
    [("staircase", LAGAN_D), ("notch", UNDER_TWO_STEPS)],
)
def test_3_either_route_covers_the_round_buffer(
    grow: ModuleType,
    staircase: Polygon,
    notch: Polygon,
    monkeypatch: pytest.MonkeyPatch,
    outline: str,
    distance: float,
    pieces: bool,
) -> None:
    """Section 4's first bullet, on both outlines and both routes (the route
    forced through `uses_pieces`): the region covers the polygon grown by
    d - 1e-5 m with round joins, 64 segments a quarter."""
    polygon = staircase if outline == "staircase" else notch
    monkeypatch.setattr(grow, "uses_pieces", lambda p, d: pieces)
    grown = grow.grow_mitred(polygon, distance)
    assert grown.covers(polygon.buffer(distance - SLACK, quad_segs=64))


# ------------------------------------------------------------------ 4. off the gate


def test_4_off_the_gate_it_is_geos_buffer(grow: ModuleType, notch: Polygon) -> None:
    """On Ola's default for section 9's question: a smooth outline of more
    than 5,000 vertices stays on GEOS, so its region is master's, bit for
    bit. (The notch is where the pieced region is not GEOS's.)"""
    assert not grow.uses_pieces(notch, UNDER_TWO_STEPS)
    grown = grow.grow_mitred(notch, UNDER_TWO_STEPS)
    assert shapely.to_wkb(grown) == shapely.to_wkb(mitred(notch, UNDER_TWO_STEPS))


# ------------------------------------------------------------------ 5. one buffer


def test_5_a_reprojected_decode_grows_the_domain_once(
    grow: ModuleType, monkeypatch: pytest.MonkeyPatch, spied: list[tuple[Any, Any]]
) -> None:
    """Section 3.1: the velhas ANADEM extract (EPSG:4674) onto a 30 m grid in
    EPSG:31983 over its catchment. `grow_mitred` runs once (every module's
    binding of it spied), and the domain is buffered by the target cell's
    diagonal once, not twice."""
    original, calls = grow.grow_mitred, []

    def spy(polygon: Polygon, distance: float) -> Polygon:
        calls.append(distance)
        return original(polygon, distance)

    for name, module in list(sys.modules.items()):
        if name.startswith("tin_engine") and module is not None:
            for attr, value in list(vars(module).items()):
                if value is original:
                    monkeypatch.setattr(module, attr, spy)
    open_dem(
        DemRequest(
            sources=(VELHAS / "anadem_velhas.tif",),
            domain=read_domain(VELHAS / "catchment.geojson"),
            target_crs=VELHAS_CRS,
        )
    )
    diagonal = math.sqrt(2) * VELHAS_H
    assert calls == [pytest.approx(diagonal)]
    assert sum(1 for _, d in spied if np.isscalar(d) and d == pytest.approx(diagonal)) == 1


# ------------------------------------------------------------------ 6. features region


@pytest.mark.parametrize("source_crs", [UTM33, "EPSG:3035"], ids=["own CRS", "LAEA"])
def test_6_the_features_region_from_the_hull(
    staircase: Polygon, spied: list[tuple[Any, Any]], source_crs: str
) -> None:
    """Section 4's third bullet, on the staircase: the region covers the
    domain grown by MARGIN - 0.5 m (round), its bounds within 0.5 m of
    master's region; and section 3.2's change itself, that only a convex
    polygon is buffered. Scale: MARGIN = 100 m, chord error 0.48 m."""
    from tin_engine.crs import reprojector
    from tin_engine.feature_input import MARGIN, source_region

    domain = DomainPolygon(polygon=staircase, crs=UTM33)
    region = source_region(domain, UTM33, source_crs)
    buffered = [g for g, _ in spied]
    spied.clear()
    assert buffered and all(g.equals(g.convex_hull) for g in buffered), [
        len(g.exterior.coords) for g in buffered
    ]
    expected = master_source_region(domain, UTM33, source_crs)
    back = Polygon(reprojector(source_crs, UTM33)(np.asarray(region.exterior.coords)))
    if CRS(source_crs) == CRS(UTM33):
        assert np.max(np.abs(np.subtract(region.bounds, expected.bounds))) <= 0.5
    else:  # the bounds compared where the margin is metres, in the DEM's CRS
        was = Polygon(reprojector(source_crs, UTM33)(np.asarray(expected.exterior.coords)))
        assert np.max(np.abs(np.subtract(back.bounds, was.bounds))) <= 0.5
    assert back.covers(staircase.buffer(MARGIN - 0.5, quad_segs=64))


# ------------------------------------------------------------------ 7. the meshes


@pytest.fixture
def as_master(grow: ModuleType) -> Callable[[], AbstractContextManager[None]]:
    """A context in which the run behaves as master: the gate forced false
    (`PIECES_FROM` above any count) and master's features region."""

    @contextmanager
    def context() -> Iterator[None]:
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(grow, "PIECES_FROM", 10**12)
            patch.setattr("tin_engine.feature_input.source_region", master_source_region)
            yield

    return context


def staircase_ring(polygon: Polygon) -> Ring:
    return [(float(x), float(y)) for x, y in polygon.exterior.coords[:-1]]


def twice(
    tmp_path: Path, as_master: Callable[[], AbstractContextManager[None]], *args: str
) -> tuple[bytes, bytes]:
    """The same `mesh` as master behaves and as built: both files' bytes."""
    with as_master():
        code, output = invoke("mesh", *args, "--out", str(tmp_path / "master.vtk"))
    assert code == 0, output
    code, output = invoke("mesh", *args, "--out", str(tmp_path / "built.vtk"))
    assert code == 0, output
    return (tmp_path / "master.vtk").read_bytes(), (tmp_path / "built.vtk").read_bytes()


def test_7a_reprojected_staircase_mesh_is_byte_identical(
    grow: ModuleType, tmp_path: Path, as_master: Callable[[], AbstractContextManager[None]]
) -> None:
    """The velhas DEM onto EPSG:31983 over a staircase domain of 5 m steps,
    radius 3.5 km, inside the extract's ~9 km: gated at the 30 m grid's
    diagonal, so the built run grows it in pieces."""
    polygon = staircase_disc(3_500.0, 5.0, (595_330.0, 7_849_860.0))
    assert len(polygon.exterior.coords) - 1 >= 5_000
    assert grow.uses_pieces(polygon, math.sqrt(2) * VELHAS_H)
    domain = geojson(tmp_path / "staircase.geojson", staircase_ring(polygon), crs=VELHAS_CRS)
    master, built = twice(
        tmp_path,
        as_master,
        *("--dem", str(VELHAS / "anadem_velhas.tif"), "--domain", str(domain)),
        *("--out-crs", VELHAS_CRS, "--tolerance", "5"),
    )
    assert master == built


@needs_codecs
def test_7b_projected_staircase_with_corine_is_byte_identical(
    grow: ModuleType, tmp_path: Path, as_master: Callable[[], AbstractContextManager[None]]
) -> None:
    """The committed DTM10 tile with the committed CORINE GeoPackage over a
    staircase domain of 5 m steps, radius 3.5 km: the projected path's
    needed region (grown by the 10 m cell's diagonal) is gated, and the
    features region comes from the hull."""
    # 23 CORINE polygons of 6 classes cross this disc; relief 333 m.
    polygon = staircase_disc(3_500.0, 5.0, (827_500.0, 7_905_000.0))
    assert grow.uses_pieces(polygon, UNDER_TWO_STEPS)
    domain = geojson(tmp_path / "staircase.geojson", staircase_ring(polygon), crs=UTM33)
    master, built = twice(
        tmp_path,
        as_master,
        *("--dem", str(KARTVERKET), "--domain", str(domain), "--tolerance", "10"),
        *("--features", str(EXTRACT), "--features-map", "corine"),
    )
    assert master == built
