"""Triangle geometry for the snow run: `snow/geometry.py` and `crs.site_of`, test R3.

`docs/increments/46-hbv-snow.md`, *4.3 The sun correction* and *4.9
Degeneracy and refusal policy*: unit normals per unit inclined area, in true
east-north-up, the grid normal turned by the grid convergence gamma at the
mesh's centre (one gamma per mesh). Section 10.5 pins the sign: true azimuth =
grid azimuth + gamma, so a grid-north vector (0, 1) becomes (sin gamma, cos
gamma). A zero-area triangle gets the vertical normal and is counted, never a
NaN; a normal pointing down is turned up. Section 7's R3: normals turned by
+0.0984 degrees at (506000, 6800000) in EPSG:25833.

Interface:

- `snow.geometry.triangle_geometry(vertices, triangles, convergence_deg) ->
  TriangleGeometry` with `normals` (`(T, 3)` float64, unit, up >= 0),
  `centroid_z`, `area` (plan area, m^2, as the blueprint names it),
  `slope_deg`, `aspect_deg` (true, from north clockwise, in [0, 360)) and
  `degenerate` (the count of zero-area triangles). A non-finite coordinate is
  a `ValueError` saying "finite".
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

F64 = npt.NDArray[np.float64]
GAMMA_EAST_EDGE = 0.0984  # deg at (506000, 6800000), EPSG:25833 (probe_sizes.out)
X0, Y0 = 506000.0, 6800000.0


@pytest.fixture(scope="module")
def geometry() -> ModuleType:
    return importlib.import_module("tin_engine.snow.geometry")


def facet(slope_deg: float, aspect_deg: float, *, side: float = 10.0, clockwise: bool = False):
    """One triangle on the plane of that slope, facing that grid aspect, at (X0, Y0)."""
    s, a = math.radians(slope_deg), math.radians(aspect_deg)
    # Downhill points along the aspect: z falls by tan(s) per metre along it.
    gx, gy = -math.tan(s) * math.sin(a), -math.tan(s) * math.cos(a)
    plan = [(0.0, 0.0), (side, 0.0), (0.0, side)]
    if clockwise:
        plan = plan[::-1]
    verts = np.array([[X0 + x, Y0 + y, 1000.0 + gx * x + gy * y] for x, y in plan])
    return verts, np.array([[0, 1, 2]])


def plane_normal(slope_deg: float, aspect_deg: float) -> F64:
    s, a = math.radians(slope_deg), math.radians(aspect_deg)
    return np.array([math.sin(s) * math.sin(a), math.sin(s) * math.cos(a), math.cos(s)])


class TestNormals:
    @pytest.mark.parametrize("aspect", [0.0, 90.0, 180.0, 270.0, 33.0])
    def test_a_30_degree_facet_has_its_planes_normal(
        self, geometry: ModuleType, aspect: float
    ) -> None:
        g = geometry.triangle_geometry(*facet(30.0, aspect), 0.0)
        assert_allclose(g.normals[0], plane_normal(30.0, aspect), atol=1e-12)
        assert g.slope_deg[0] == pytest.approx(30.0, abs=1e-9)
        assert g.aspect_deg[0] == pytest.approx(aspect, abs=1e-9)

    def test_north_is_aspect_0_not_360(self, geometry: ModuleType) -> None:
        g = geometry.triangle_geometry(*facet(30.0, 0.0), 0.0)
        assert 0.0 <= g.aspect_deg[0] < 1e-9

    def test_a_clockwise_triangle_is_turned_up(self, geometry: ModuleType) -> None:
        g = geometry.triangle_geometry(*facet(30.0, 180.0, clockwise=True), 0.0)
        assert_allclose(g.normals[0], plane_normal(30.0, 180.0), atol=1e-12)

    def test_normals_are_unit_and_never_point_down(self, geometry: ModuleType) -> None:
        rng = np.random.default_rng(463)
        verts = np.column_stack(
            [X0 + rng.random(60) * 500, Y0 + rng.random(60) * 500, 900 + rng.random(60) * 300]
        )
        tris = rng.integers(0, 60, size=(200, 3))
        tris = tris[
            (tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 0] != tris[:, 2])
        ]
        g = geometry.triangle_geometry(verts, tris, 0.0)
        assert_allclose(np.linalg.norm(g.normals, axis=1), 1.0, atol=1e-12)
        assert (g.normals[:, 2] >= 0.0).all()

    def test_a_flat_triangle_points_straight_up(self, geometry: ModuleType) -> None:
        g = geometry.triangle_geometry(*facet(0.0, 0.0), 0.0)
        assert_array_equal(g.normals[0], [0.0, 0.0, 1.0])
        assert g.slope_deg[0] == 0.0

    def test_a_vertical_triangle_has_a_level_normal_and_is_not_degenerate(
        self, geometry: ModuleType
    ) -> None:
        verts = np.array([[X0, Y0, 0.0], [X0 + 10.0, Y0, 0.0], [X0, Y0, 10.0]])
        g = geometry.triangle_geometry(verts, np.array([[0, 1, 2]]), 0.0)
        assert_allclose(np.abs(g.normals[0]), [0.0, 1.0, 0.0], atol=1e-12)
        assert g.slope_deg[0] == pytest.approx(90.0)
        assert g.degenerate == 0
        assert g.area[0] == 0.0

    def test_a_millimetre_triangle_at_utm_northings_keeps_its_normal(
        self, geometry: ModuleType
    ) -> None:
        # Scale: 1 mm edges at 6.8e6 m, where a double's step is 9.3e-10 m (1e-6 relative).
        verts, tris = facet(30.0, 180.0, side=1e-3)
        g = geometry.triangle_geometry(verts, tris, 0.0)
        assert_allclose(g.normals[0], plane_normal(30.0, 180.0), atol=1e-5)


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
        g = geometry.triangle_geometry(verts, tris, 0.0)
        assert g.degenerate == 2
        assert_array_equal(g.normals[:2], [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])
        assert np.isfinite(g.normals).all()
        assert np.isfinite(g.slope_deg).all() and np.isfinite(g.aspect_deg).all()
        assert g.normals[2, 2] < 1.0

    def test_an_empty_mesh_gives_empty_arrays(self, geometry: ModuleType) -> None:
        g = geometry.triangle_geometry(np.empty((0, 3)), np.empty((0, 3), dtype=np.int64), 0.0)
        assert g.normals.shape == (0, 3)
        assert g.area.shape == g.centroid_z.shape == (0,)
        assert g.degenerate == 0

    @pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
    def test_a_non_finite_coordinate_is_refused(self, geometry: ModuleType, bad: float) -> None:
        verts, tris = facet(30.0, 180.0)
        verts[1, 2] = bad
        with pytest.raises(ValueError, match="finite"):
            geometry.triangle_geometry(verts, tris, 0.0)


class TestMeasures:
    def test_area_is_the_plan_area_whatever_the_slope(self, geometry: ModuleType) -> None:
        for slope in (0.0, 30.0, 60.0):
            g = geometry.triangle_geometry(*facet(slope, 135.0), 0.0)
            assert g.area[0] == pytest.approx(50.0, rel=1e-12)

    def test_centroid_z_is_the_corner_mean(self, geometry: ModuleType) -> None:
        verts, tris = facet(30.0, 180.0)
        g = geometry.triangle_geometry(verts, tris, 0.0)
        assert g.centroid_z[0] == pytest.approx(verts[:, 2].mean(), abs=1e-12)


class TestConvergence:
    """R3, and section 10.5's sign: true azimuth = grid azimuth + gamma."""

    def test_site_of_gives_the_east_edges_convergence(self) -> None:
        site = importlib.import_module("tin_engine.crs").site_of("EPSG:25833", X0, Y0)
        assert site.convergence_deg == pytest.approx(GAMMA_EAST_EDGE, abs=5e-5)
        assert site.lat_deg == pytest.approx(61.334, abs=5e-4)
        assert site.lon_deg == pytest.approx(15.112, abs=5e-4)

    def test_the_central_meridian_has_no_convergence(self) -> None:
        site = importlib.import_module("tin_engine.crs").site_of("EPSG:25833", 500000.0, Y0)
        assert site.convergence_deg == pytest.approx(0.0, abs=1e-9)

    def test_a_geographic_crs_is_refused_in_degrees(self) -> None:
        with pytest.raises(ValueError, match="degrees"):
            importlib.import_module("tin_engine.crs").site_of("EPSG:4326", 15.0, 61.0)

    def test_a_grid_north_facet_faces_true_azimuth_gamma(self, geometry: ModuleType) -> None:
        g = geometry.triangle_geometry(*facet(30.0, 0.0), GAMMA_EAST_EDGE)
        gamma = math.radians(GAMMA_EAST_EDGE)
        sin30 = 0.5
        assert_allclose(
            g.normals[0],
            [sin30 * math.sin(gamma), sin30 * math.cos(gamma), math.cos(math.radians(30.0))],
            atol=1e-12,
        )
        assert g.aspect_deg[0] == pytest.approx(GAMMA_EAST_EDGE, abs=1e-9)

    def test_the_turn_is_a_rotation_about_up(self, geometry: ModuleType) -> None:
        for aspect in (90.0, 200.0, 359.0):
            g = geometry.triangle_geometry(*facet(25.0, aspect), -1.5)
            assert g.aspect_deg[0] == pytest.approx((aspect - 1.5) % 360.0, abs=1e-9)
            assert g.slope_deg[0] == pytest.approx(25.0, abs=1e-9)
