"""Triangle geometry for the snow run: `snow/geometry.py` and `crs.site_of`.

`docs/increments/46-hbv-snow.md`, *4.9 Degeneracy and refusal policy* and
*10.5 How the lit factor enters the sun correction* (one source of normals):
each triangle's normal, slope and aspect are the receiver frame's per-triangle
true east-north-up normal (`radiation.frame.receivers_from_mesh`, through its
Jacobian J at the triangle's own centroid), never a grid normal turned by one
convergence angle for the whole mesh. A mesh in EPSG:25833 that reaches from
the central meridian (15 E, convergence 0) to Bygdin (8.49 E, convergence
-5.72 degrees) gives each facet its true aspect, which no single turn can do,
and a turn left in or applied twice moves Bygdin's aspect by 5.7 degrees. A
zero-area triangle, by the frame's own test (its ground edges span less than
1e-9 of the product of their lengths, `ZERO_AREA_SINE`), gets the vertical
normal and is counted, never a NaN and never refused; a normal pointing down
is turned up. `site_of` gives the sun's site (and its convergence, which the
normals no longer use) at (506000, 6800000) in EPSG:25833.

Interface under test:

- `snow.geometry.triangle_geometry(vertices, triangles, crs) ->
  TriangleGeometry`, `crs` the mesh's projected CRS, with `normals` (`(T, 3)`
  float64, unit, up >= 0, true east-north-up), `centroid_z`, `area` (plan
  area in the CRS, m^2, as the blueprint names it), `slope_deg`, `aspect_deg`
  (true, from north clockwise, in [0, 360)), `zero_area` (bool `(T,)`) and
  `degenerate` (their count). For every triangle not `zero_area`, `normals`
  equals `receivers_from_mesh(vertices, triangles[~zero_area], crs).normal`
  bit for bit. A non-finite coordinate is a `ValueError` saying "finite".
- `crs.site_of(crs, x, y) -> Site` with `lat_deg`, `lon_deg` and
  `convergence_deg` (pyproj's meridian convergence) at the point; a
  geographic CRS is a `ValueError` saying "degrees".
"""

from __future__ import annotations

import importlib
import math
from types import ModuleType

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pyproj import Proj

from tin_engine.radiation.frame import receivers_from_mesh

F64 = npt.NDArray[np.float64]
CRS = "EPSG:25833"
GAMMA_EAST_EDGE = 0.0984  # deg at (506000, 6800000), EPSG:25833 (probe_sizes.out)
X0, Y0 = 506000.0, 6800000.0
#: Bygdin's centre (61.3694 N, 8.4852 E) in EPSG:25833, convergence -5.724 degrees.
BYGDIN = (152147.55, 6821319.44)
#: On UTM zone 33's central meridian, convergence 0.
CENTRAL = (500000.0, 6800000.0)
#: Slope and aspect agree with the plane built through pyproj's factors to
#: this, in degrees: at 10 m edges J changes by under 1e-8 relative across a
#: triangle, checked at Bygdin, 350 km off the central meridian.
TRUE_DEG = 1e-6


@pytest.fixture(scope="module")
def geometry() -> ModuleType:
    return importlib.import_module("tin_engine.snow.geometry")


def facet(
    slope_deg: float,
    aspect_deg: float,
    *,
    at: tuple[float, float] = (X0, Y0),
    side: float = 10.0,
    clockwise: bool = False,
) -> tuple[F64, npt.NDArray[np.int64]]:
    """One triangle at `at` on the plane of that true slope and true aspect.

    The plane is built in the ground frame and taken to the grid through J =
    (1/m) R(gamma), pyproj's convergence and scale at the centroid: z falls by
    tan(slope) per ground metre along the true aspect."""
    plan = [(0.0, 0.0), (side, 0.0), (0.0, side)]
    if clockwise:
        plan = plan[::-1]
    proj = Proj(CRS)
    lon, lat = proj(at[0] + side / 3, at[1] + side / 3, inverse=True)
    f = proj.get_factors(lon, lat)
    g, m = math.radians(f.meridian_convergence), f.meridional_scale
    j = np.array([[math.cos(g), math.sin(g)], [-math.sin(g), math.cos(g)]]) / m
    s, a = math.radians(slope_deg), math.radians(aspect_deg)
    grad = -math.tan(s) * j.T @ np.array([math.sin(a), math.cos(a)])
    verts = np.array([[at[0] + x, at[1] + y, 1000.0 + grad @ (x, y)] for x, y in plan])
    return verts, np.array([[0, 1, 2]])


def plane_normal(slope_deg: float, aspect_deg: float) -> F64:
    s, a = math.radians(slope_deg), math.radians(aspect_deg)
    return np.array([math.sin(s) * math.sin(a), math.sin(s) * math.cos(a), math.cos(s)])


def joined(*meshes: tuple[F64, npt.NDArray[np.int64]]) -> tuple[F64, npt.NDArray[np.int64]]:
    verts, tris, base = [], [], 0
    for v, t in meshes:
        verts.append(v)
        tris.append(t + base)
        base += len(v)
    return np.concatenate(verts), np.concatenate(tris)


class TestOneSourceOfNormals:
    """Section 10.5: the frame's per-triangle normals, not one turn per mesh."""

    def test_each_facet_keeps_its_true_aspect_far_from_the_central_meridian(
        self, geometry: ModuleType
    ) -> None:
        verts, tris = joined(facet(30.0, 180.0, at=BYGDIN), facet(20.0, 90.0, at=CENTRAL))
        g = geometry.triangle_geometry(verts, tris, CRS)
        assert_allclose(g.slope_deg, [30.0, 20.0], rtol=0, atol=TRUE_DEG)
        assert_allclose(g.aspect_deg, [180.0, 90.0], rtol=0, atol=TRUE_DEG)

    def test_the_normals_are_the_receiver_frames_bit_for_bit(self, geometry: ModuleType) -> None:
        verts, tris = joined(
            facet(30.0, 180.0, at=BYGDIN),
            facet(35.0, 300.0, at=BYGDIN, clockwise=True),
            facet(20.0, 90.0, at=CENTRAL),
        )
        g = geometry.triangle_geometry(verts, tris, CRS)
        frame = receivers_from_mesh(verts, tris, CRS)
        assert g.degenerate == 0
        assert_array_equal(g.normals, frame.normal)
        e, n, u = frame.normal.T
        assert_array_equal(g.slope_deg, np.degrees(np.arctan2(np.hypot(e, n), u)))
        assert_array_equal(g.aspect_deg, np.degrees(np.arctan2(e, n)) % 360.0)

    def test_the_kept_triangles_normals_are_the_frames_of_the_kept_triangles(
        self, geometry: ModuleType
    ) -> None:
        line = np.array([[X0, Y0, 5.0], [X0 + 10.0, Y0 + 10.0, 5.0], [X0 + 20.0, Y0 + 20.0, 5.0]])
        verts, tris = joined(
            facet(30.0, 180.0, at=BYGDIN), (line, np.array([[0, 1, 2]])), facet(20.0, 90.0)
        )
        g = geometry.triangle_geometry(verts, tris, CRS)
        assert_array_equal(g.zero_area, [False, True, False])
        kept = receivers_from_mesh(verts, tris[~g.zero_area], CRS)
        assert_array_equal(g.normals[~g.zero_area], kept.normal)


class TestNormals:
    @pytest.mark.parametrize("aspect", [0.0, 90.0, 180.0, 270.0, 33.0])
    def test_a_30_degree_facet_has_its_planes_normal(
        self, geometry: ModuleType, aspect: float
    ) -> None:
        g = geometry.triangle_geometry(*facet(30.0, aspect), CRS)
        assert_allclose(g.normals[0], plane_normal(30.0, aspect), rtol=0, atol=1e-8)
        assert g.slope_deg[0] == pytest.approx(30.0, abs=TRUE_DEG)
        got = (g.aspect_deg[0] - aspect + 180.0) % 360.0 - 180.0
        assert got == pytest.approx(0.0, abs=TRUE_DEG)

    def test_aspects_are_in_0_to_360(self, geometry: ModuleType) -> None:
        verts, tris = joined(*(facet(30.0, a, at=(X0 + 50.0 * k, Y0)) for k, a in
                               enumerate((0.0, 359.9999, 0.0001, 180.0))))  # fmt: skip
        g = geometry.triangle_geometry(verts, tris, CRS)
        assert ((g.aspect_deg >= 0.0) & (g.aspect_deg < 360.0)).all()

    def test_a_clockwise_triangle_is_turned_up(self, geometry: ModuleType) -> None:
        g = geometry.triangle_geometry(*facet(30.0, 180.0, clockwise=True), CRS)
        assert_allclose(g.normals[0], plane_normal(30.0, 180.0), rtol=0, atol=1e-8)

    def test_normals_are_unit_and_never_point_down(self, geometry: ModuleType) -> None:
        rng = np.random.default_rng(463)
        verts = np.column_stack(
            [X0 + rng.random(60) * 500, Y0 + rng.random(60) * 500, 900 + rng.random(60) * 300]
        )
        tris = rng.integers(0, 60, size=(200, 3))
        tris = tris[
            (tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 0] != tris[:, 2])
        ]
        g = geometry.triangle_geometry(verts, tris, CRS)
        assert_allclose(np.linalg.norm(g.normals, axis=1), 1.0, rtol=0, atol=1e-12)
        assert (g.normals[:, 2] >= 0.0).all()

    def test_a_flat_triangle_points_straight_up(self, geometry: ModuleType) -> None:
        g = geometry.triangle_geometry(*facet(0.0, 0.0), CRS)
        assert_array_equal(g.normals[0], [0.0, 0.0, 1.0])
        assert g.slope_deg[0] == 0.0

    def test_a_vertical_triangle_has_a_level_normal_and_is_not_degenerate(
        self, geometry: ModuleType
    ) -> None:
        verts = np.array([[X0, Y0, 0.0], [X0 + 10.0, Y0, 0.0], [X0, Y0, 10.0]])
        g = geometry.triangle_geometry(verts, np.array([[0, 1, 2]]), CRS)
        assert g.normals[0, 2] == pytest.approx(0.0, abs=1e-12)
        assert g.slope_deg[0] == pytest.approx(90.0)
        assert g.degenerate == 0
        assert g.area[0] == 0.0

    def test_a_millimetre_triangle_at_utm_northings_keeps_its_normal(
        self, geometry: ModuleType
    ) -> None:
        # Scale: 1 mm edges at 6.8e6 m, where a double's step is 9.3e-10 m (1e-6 relative).
        verts, tris = facet(30.0, 180.0, side=1e-3)
        g = geometry.triangle_geometry(verts, tris, CRS)
        assert_allclose(g.normals[0], plane_normal(30.0, 180.0), rtol=0, atol=1e-5)


class TestDegenerate:
    def test_zero_area_triangles_get_the_vertical_normal_and_are_counted(
        self, geometry: ModuleType
    ) -> None:
        verts = np.array(
            [
                [X0, Y0, 100.0],
                [X0 + 10.0, Y0 + 10.0, 110.0],
                [X0 + 20.0, Y0 + 20.0, 120.0],  # collinear with the first two
                [X0 + 5.0, Y0, 100.0],
            ]
        )
        tris = np.array([[0, 1, 2], [0, 0, 3], [0, 3, 1]])  # collinear, repeated corner, fine
        g = geometry.triangle_geometry(verts, tris, CRS)
        assert g.degenerate == 2
        assert_array_equal(g.zero_area, [True, True, False])
        assert_array_equal(g.normals[:2], [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])
        assert np.isfinite(g.normals).all()
        assert np.isfinite(g.slope_deg).all() and np.isfinite(g.aspect_deg).all()
        assert g.normals[2, 2] < 1.0

    def test_a_sliver_under_the_frames_zero_area_limit_is_counted_not_refused(
        self, geometry: ModuleType
    ) -> None:
        # Edges (10, 10) and (20, 20 + 1e-8) m: their cross product is 1e-7 m^2,
        # not 0, but the sine between them is 2.5e-10, under the frame's 1e-9.
        verts = np.array([[X0, Y0, 300.0], [X0 + 10.0, Y0 + 10.0, 300.0],
                          [X0 + 20.0, Y0 + 20.0 + 1e-8, 300.0]])  # fmt: skip
        sliver, fine = joined((verts, np.array([[0, 1, 2]])), facet(30.0, 180.0))
        assert np.cross(sliver[1] - sliver[0], sliver[2] - sliver[0])[2] != 0.0
        g = geometry.triangle_geometry(sliver, fine, CRS)
        assert g.degenerate == 1
        assert_array_equal(g.zero_area, [True, False])
        assert_array_equal(g.normals[0], [0.0, 0.0, 1.0])
        assert g.slope_deg[0] == 0.0

    def test_an_empty_mesh_gives_empty_arrays(self, geometry: ModuleType) -> None:
        g = geometry.triangle_geometry(np.empty((0, 3)), np.empty((0, 3), dtype=np.int64), CRS)
        assert g.normals.shape == (0, 3)
        assert g.area.shape == g.centroid_z.shape == (0,)
        assert g.degenerate == 0

    @pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
    def test_a_non_finite_coordinate_is_refused(self, geometry: ModuleType, bad: float) -> None:
        verts, tris = facet(30.0, 180.0)
        verts[1, 2] = bad
        with pytest.raises(ValueError, match="finite"):
            geometry.triangle_geometry(verts, tris, CRS)


class TestMeasures:
    def test_area_is_the_plan_area_whatever_the_slope(self, geometry: ModuleType) -> None:
        for slope in (0.0, 30.0, 60.0):
            g = geometry.triangle_geometry(*facet(slope, 135.0), CRS)
            assert g.area[0] == pytest.approx(50.0, rel=1e-12)

    def test_centroid_z_is_the_corner_mean(self, geometry: ModuleType) -> None:
        verts, tris = facet(30.0, 180.0)
        g = geometry.triangle_geometry(verts, tris, CRS)
        assert g.centroid_z[0] == pytest.approx(verts[:, 2].mean(), abs=1e-12)


class TestSite:
    """`site_of`: the sun's site, and pyproj's convergence at it."""

    def test_site_of_gives_the_east_edges_convergence(self) -> None:
        site = importlib.import_module("tin_engine.crs").site_of(CRS, X0, Y0)
        assert site.convergence_deg == pytest.approx(GAMMA_EAST_EDGE, abs=5e-5)
        assert site.lat_deg == pytest.approx(61.334, abs=5e-4)
        assert site.lon_deg == pytest.approx(15.112, abs=5e-4)

    def test_the_central_meridian_has_no_convergence(self) -> None:
        site = importlib.import_module("tin_engine.crs").site_of(CRS, 500000.0, Y0)
        assert site.convergence_deg == pytest.approx(0.0, abs=1e-9)

    def test_a_geographic_crs_is_refused_in_degrees(self) -> None:
        with pytest.raises(ValueError, match="degrees"):
            importlib.import_module("tin_engine.crs").site_of("EPSG:4326", 15.0, 61.0)

    @pytest.mark.parametrize(
        ("crs", "unit"), [("EPSG:2263", "US survey foot"), ("EPSG:3359", "foot")]
    )
    def test_a_projected_crs_in_feet_is_refused_naming_its_unit(self, crs: str, unit: str) -> None:
        # Code review round 1 of PR 3: slopes need metres on all three axes.
        with pytest.raises(ValueError, match="metres") as caught:
            importlib.import_module("tin_engine.crs").site_of(crs, 1_000_000.0, 200_000.0)
        assert unit in str(caught.value)
