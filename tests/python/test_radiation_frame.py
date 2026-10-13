"""The receiver frame: `crs.frame_factors`, `radiation/frame.py`, the request.

`docs/increments/46-hbv-snow.md`, *10.8 Red suites for `@tester`* ("PR 4,
TR2"), and the radiation design (`docs/increments/terrain-radiation.md`, R3,
R10, I3, its degeneracy policy and verification 2): J against a
finite-difference transform at Solund, Oslo, the central meridian, Vardo and
Bygdin; the sign of the convergence through J (true azimuth = grid azimuth +
gamma); a synthetic plane off the central meridian keeps its true aspect;
normals, exclusion radii, observer heights, latitudes and longitudes per
triangle; refusals in plain words; `RadiationRequest`'s defaults and limits.

Interface assumed (names chosen here and listed in the handback):

- `tin_engine.crs.frame_factors(crs, lon, lat) -> (gamma_deg, m)`: two float64
  arrays the length of `lon`, pyproj's `meridian_convergence` (degrees, the
  true azimuth of grid north) and the point scale (`meridional_scale`), from
  `Proj(crs).get_factors`. A geographic CRS is a `ValueError` saying a
  projected CRS is needed.
- `tin_engine.radiation.frame.receivers_from_mesh(points, triangles, crs, *,
  observer_offset=0.0) -> ReceiverFrame`: `points` an (n, 3) float64 array
  (x, y in `crs`, z in metres), `triangles` an (m, 3) integer array. Per
  triangle, at its centroid: `x`, `y`, `z_obs` (mean vertex z plus the
  offset), `lat`, `lon` (degrees, `crs.reprojector` to EPSG:4326),
  `jacobian` (m, 2, 2) with (east, north) = J (dx, dy) and
  J = (1/m) [[cos g, sin g], [-sin g, cos g]], `normal` (m, 3) the unit
  upward normal in true (east, north, up), `exclusion_radius` the largest
  ground distance |J (v - c)| from the centroid to a vertex; and
  `earth_radius_m`, sqrt(M N) on the CRS's ellipsoid at the mean latitude.
  Refused with a `ValueError`: a zero-area triangle, a non-finite vertex, a
  vertex index out of range (each naming "triangle i", 0-based), a negative
  offset, a receiver above 89.9 degrees of latitude, a geographic CRS.
- `tin_engine._core.check_receivers(x, y, z_obs, lat_deg, lon_deg, jacobian,
  normal, exclusion_radius)`: the core's check of a receiver set; `None`, or
  a `ValueError` naming "receiver i" and the reason.
- `tin_engine.radiation.request.RadiationRequest`: frozen Pydantic, extra
  fields forbidden, every field defaulted: `n_directions` 64 (a power of two,
  at least 4), `refraction_k` 0.13 (in [0, 1)), `d_near_m` 20 000 (> 0),
  `d_far_m` 200 000 (>= `d_near_m`), `limb` "centre" or "upper" ("centre").

Reference values come from pyproj built in this file (`get_factors`, and an
always_xy `Transformer` for the finite differences), never from the module
under test; the design's rounded figures (P4, probe B0) are checked as well.

Numeric bounds: J entries are of order 1; the finite difference (1 m steps
at 10^6 m coordinates) agrees with (1/m) Rot(gamma) to 1.6e-9 at all five
sites (probe in the commit message), so 1e-7. Angles to 1e-6 degrees where
both sides are computed from the same factors, 1e-3 against the design's
rounded figures.
"""

from __future__ import annotations

import importlib
import math
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from pyproj import CRS, Proj, Transformer

F64 = npt.NDArray[np.float64]
FRAME = "EPSG:25833"

# GRS80, EPSG:25833's ellipsoid.
A_GRS80 = 6378137.0
E2_GRS80 = (1 / 298.257222101) * (2 - 1 / 298.257222101)

# (lon, lat) and the design's rounded gamma (P4; probe B0 for Bygdin).
SITES = {
    "solund": (4.85, 61.08, -8.906),
    "oslo": (10.75, 59.91, -3.679),
    "central meridian": (15.00, 65.00, 0.000),
    "vardo": (31.11, 70.37, 15.218),
    "bygdin": (8.4852, 61.3694, -5.724),
}

# Offsets of a small triangle's vertices from its centroid, in frame metres;
# they sum to zero, so the centroid is the site itself. Counter-clockwise.
OFFSETS = np.array([[-10.0, -10.0], [10.0, -5.0], [0.0, 15.0]])


@pytest.fixture(scope="module")
def crs_mod() -> ModuleType:
    return importlib.import_module("tin_engine.crs")


@pytest.fixture(scope="module")
def frame() -> ModuleType:
    return importlib.import_module("tin_engine.radiation.frame")


@pytest.fixture(scope="module")
def request_mod() -> ModuleType:
    return importlib.import_module("tin_engine.radiation.request")


@pytest.fixture(scope="module")
def core() -> Any:
    import tin_engine._core as core

    return core


def site_xy(lon: float, lat: float) -> tuple[float, float]:
    x, y = Transformer.from_crs("EPSG:4326", FRAME, always_xy=True).transform(lon, lat)
    return float(x), float(y)


def factors(lon: float, lat: float) -> tuple[float, float]:
    f = Proj(CRS.from_user_input(FRAME)).get_factors(lon, lat)
    return float(f.meridian_convergence), float(f.meridional_scale)


def rot_j(gamma_deg: float, m: float) -> F64:
    c, s = math.cos(math.radians(gamma_deg)), math.sin(math.radians(gamma_deg))
    return np.array([[c, s], [-s, c]]) / m


def radii(lat_deg: float) -> tuple[float, float]:
    """The meridional (M) and prime-vertical (N) radii of curvature on GRS80."""
    s2 = math.sin(math.radians(lat_deg)) ** 2
    w = 1 - E2_GRS80 * s2
    return A_GRS80 * (1 - E2_GRS80) / w**1.5, A_GRS80 / math.sqrt(w)


def fd_jacobian(x: float, y: float) -> F64:
    """J by central differences: 1 m frame steps to lon/lat, then to local
    east and north metres on the ellipsoid."""
    inv = Transformer.from_crs(FRAME, "EPSG:4326", always_xy=True)
    _, lat0 = inv.transform(x, y)
    m_r, n_r = radii(lat0)
    j = np.zeros((2, 2))
    for col, (dx, dy) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        lo1, la1 = inv.transform(x + dx, y + dy)
        lo0, la0 = inv.transform(x - dx, y - dy)
        j[0, col] = n_r * math.cos(math.radians(lat0)) * math.radians(lo1 - lo0) / 2
        j[1, col] = m_r * math.radians(la1 - la0) / 2
    return j


def triangle_at(
    x: float, y: float, z: tuple[float, float, float] = (0.0, 0.0, 0.0)
) -> tuple[F64, npt.NDArray[np.int64]]:
    pts = np.column_stack([OFFSETS + np.array([x, y]), np.asarray(z, dtype=np.float64)])
    return pts, np.array([[0, 1, 2]], dtype=np.int64)


def azimuth(east: float, north: float) -> float:
    return math.degrees(math.atan2(east, north)) % 360.0


def plane_heights(j: F64, slope_deg: float, aspect_deg: float) -> tuple[float, float, float]:
    """Heights at OFFSETS of a plane of the given true slope and aspect (the
    downslope true azimuth): z = g . J (dx, dy), the upward normal being
    proportional to (-g_E, -g_N, 1)."""
    t = math.tan(math.radians(slope_deg))
    a = math.radians(aspect_deg)
    g = np.array([-t * math.sin(a), -t * math.cos(a)])
    z = (OFFSETS @ j.T) @ g
    return float(z[0]), float(z[1]), float(z[2])


# --- crs.frame_factors -----------------------------------------------------


class TestFrameFactors:
    def test_matches_pyproj_and_the_design_at_every_site(self, crs_mod: ModuleType) -> None:
        lon = np.array([s[0] for s in SITES.values()])
        lat = np.array([s[1] for s in SITES.values()])
        gamma, m = crs_mod.frame_factors(FRAME, lon, lat)
        assert np.shape(gamma) == (5,) and np.shape(m) == (5,)
        ref = Proj(CRS.from_user_input(FRAME)).get_factors(lon, lat)
        np.testing.assert_allclose(gamma, ref.meridian_convergence, rtol=0, atol=1e-12)
        np.testing.assert_allclose(m, ref.meridional_scale, rtol=0, atol=1e-12)
        design = np.array([s[2] for s in SITES.values()])
        np.testing.assert_allclose(gamma, design, rtol=0, atol=1e-3)

    def test_takes_a_pyproj_crs(self, crs_mod: ModuleType) -> None:
        gamma, m = crs_mod.frame_factors(CRS.from_epsg(25833), np.array([31.11]), np.array([70.37]))
        assert float(gamma[0]) == pytest.approx(15.218, abs=1e-3)
        assert float(m[0]) == pytest.approx(1.00397, abs=1e-5)

    def test_empty_input_gives_empty_arrays(self, crs_mod: ModuleType) -> None:
        gamma, m = crs_mod.frame_factors(FRAME, np.empty(0), np.empty(0))
        assert np.shape(gamma) == (0,) and np.shape(m) == (0,)

    def test_a_geographic_crs_is_refused(self, crs_mod: ModuleType) -> None:
        with pytest.raises(ValueError, match="projected"):
            crs_mod.frame_factors("EPSG:4326", np.array([8.0]), np.array([61.0]))


# --- J and the sign of the convergence, through receivers_from_mesh -------


class TestJacobian:
    @pytest.mark.parametrize("site", list(SITES))
    def test_j_equals_a_finite_difference_transform(self, frame: ModuleType, site: str) -> None:
        lon, lat, _ = SITES[site]
        x, y = site_xy(lon, lat)
        pts, tri = triangle_at(x, y)
        rf = frame.receivers_from_mesh(pts, tri, FRAME)
        assert rf.jacobian.shape == (1, 2, 2)
        np.testing.assert_allclose(rf.jacobian[0], fd_jacobian(x, y), rtol=0, atol=1e-7)
        np.testing.assert_allclose(rf.jacobian[0], rot_j(*factors(lon, lat)), rtol=0, atol=1e-12)

    @pytest.mark.parametrize("site", list(SITES))
    def test_grid_north_has_true_azimuth_gamma(self, frame: ModuleType, site: str) -> None:
        # True azimuth = grid azimuth + gamma: Solund -8.906, Vardo +15.218,
        # Bygdin -5.724 (the catchment's centre), through J, not through gamma.
        lon, lat, design_gamma = SITES[site]
        x, y = site_xy(lon, lat)
        rf = frame.receivers_from_mesh(*triangle_at(x, y), FRAME)
        east, north = rf.jacobian[0] @ np.array([0.0, 1.0])
        got = (azimuth(east, north) + 180.0) % 360.0 - 180.0
        assert got == pytest.approx(design_gamma, abs=1e-3)

    def test_solund_true_north_is_grid_azimuth_plus_8_906(self, frame: ModuleType) -> None:
        lon, lat, _ = SITES["solund"]
        x, y = site_xy(lon, lat)
        rf = frame.receivers_from_mesh(*triangle_at(x, y), FRAME)
        dx, dy = np.linalg.solve(rf.jacobian[0], np.array([0.0, 1.0]))
        assert azimuth(dx, dy) == pytest.approx(8.906, abs=1e-3)


# --- normals, the synthetic plane ------------------------------------------


class TestNormal:
    @pytest.mark.parametrize(
        ("site", "slope", "aspect"),
        [("solund", 30.0, 180.0), ("vardo", 20.0, 45.0), ("bygdin", 35.0, 300.0)],
    )
    def test_a_plane_off_the_central_meridian_keeps_its_true_aspect(
        self, frame: ModuleType, site: str, slope: float, aspect: float
    ) -> None:
        # Without the convergence the aspect would be off by gamma (8.9, 15.2
        # and 5.7 degrees here), the legacy's error.
        lon, lat, _ = SITES[site]
        x, y = site_xy(lon, lat)
        z = plane_heights(rot_j(*factors(lon, lat)), slope, aspect)
        rf = frame.receivers_from_mesh(*triangle_at(x, y, z), FRAME)
        e, n, u = rf.normal[0]
        assert math.degrees(math.atan2(math.hypot(e, n), u)) == pytest.approx(slope, abs=1e-6)
        assert azimuth(e, n) == pytest.approx(aspect, abs=1e-6)

    def test_the_normal_is_unit_and_upward_for_either_winding(self, frame: ModuleType) -> None:
        x, y = site_xy(*SITES["oslo"][:2])
        pts, _ = triangle_at(x, y, (100.0, 140.0, 90.0))
        tri = np.array([[0, 1, 2], [0, 2, 1]], dtype=np.int64)
        rf = frame.receivers_from_mesh(pts, tri, FRAME)
        assert rf.normal.shape == (2, 3)
        np.testing.assert_allclose(np.linalg.norm(rf.normal, axis=1), 1.0, rtol=0, atol=1e-12)
        assert np.all(rf.normal[:, 2] > 0)
        np.testing.assert_allclose(rf.normal[0], rf.normal[1], rtol=0, atol=1e-15)

    def test_a_flat_triangle_points_straight_up(self, frame: ModuleType) -> None:
        x, y = site_xy(*SITES["vardo"][:2])
        rf = frame.receivers_from_mesh(*triangle_at(x, y, (7.0, 7.0, 7.0)), FRAME)
        np.testing.assert_allclose(rf.normal[0], [0.0, 0.0, 1.0], rtol=0, atol=1e-15)

    def test_a_vertical_triangle_is_taken_with_slope_90(self, frame: ModuleType) -> None:
        # Collinear in plan, not in 3-D: n_U = 0 is allowed.
        x, y = site_xy(*SITES["central meridian"][:2])
        pts = np.array([[x - 10, y, 0.0], [x + 10, y, 0.0], [x, y, 30.0]])
        rf = frame.receivers_from_mesh(pts, np.array([[0, 1, 2]]), FRAME)
        e, n, u = rf.normal[0]
        assert u == pytest.approx(0.0, abs=1e-12)
        assert math.hypot(e, n) == pytest.approx(1.0, abs=1e-12)
        # At the central meridian grid east is true east, so the face looks north or south.
        assert abs(e) == pytest.approx(0.0, abs=1e-12)


# --- position, height, exclusion radius, Earth radius ----------------------


class TestPlacement:
    def test_lat_lon_are_the_centroid_through_pyproj(self, frame: ModuleType) -> None:
        sites = [SITES[s][:2] for s in ("solund", "bygdin", "vardo")]
        pts_l, tri_l = [], []
        for k, (lon, lat) in enumerate(sites):
            p, _ = triangle_at(*site_xy(lon, lat))
            pts_l.append(p)
            tri_l.append([3 * k, 3 * k + 1, 3 * k + 2])
        rf = frame.receivers_from_mesh(np.vstack(pts_l), np.array(tri_l), FRAME)
        cx = np.array([site_xy(lon, lat)[0] for lon, lat in sites])
        cy = np.array([site_xy(lon, lat)[1] for lon, lat in sites])
        np.testing.assert_allclose(rf.x, cx, rtol=0, atol=1e-8)
        np.testing.assert_allclose(rf.y, cy, rtol=0, atol=1e-8)
        lon_ref, lat_ref = Transformer.from_crs(FRAME, "EPSG:4326", always_xy=True).transform(
            rf.x, rf.y
        )
        np.testing.assert_allclose(rf.lon, lon_ref, rtol=0, atol=1e-9)
        np.testing.assert_allclose(rf.lat, lat_ref, rtol=0, atol=1e-9)

    def test_z_obs_is_the_mean_height_plus_the_offset(self, frame: ModuleType) -> None:
        x, y = site_xy(*SITES["bygdin"][:2])
        pts, tri = triangle_at(x, y, (1500.0, 1530.0, 1512.0))
        assert frame.receivers_from_mesh(pts, tri, FRAME).z_obs[0] == pytest.approx(1514.0)
        got = frame.receivers_from_mesh(pts, tri, FRAME, observer_offset=2.5).z_obs[0]
        assert got == pytest.approx(1516.5)

    def test_exclusion_radius_is_the_largest_ground_distance_to_a_vertex(
        self, frame: ModuleType
    ) -> None:
        # The farthest vertex is 15 frame metres north of the centroid; at
        # Vardo the point scale 1.00397 makes that 14.9406 ground metres.
        lon, lat, _ = SITES["vardo"]
        x, y = site_xy(lon, lat)
        rf = frame.receivers_from_mesh(*triangle_at(x, y, (0.0, 50.0, 400.0)), FRAME)
        _, m = factors(lon, lat)
        expected = max(float(np.linalg.norm(rot_j(*factors(lon, lat)) @ d)) for d in OFFSETS)
        assert expected == pytest.approx(15.0 / m, rel=1e-12)
        assert rf.exclusion_radius[0] == pytest.approx(expected, rel=1e-12)

    def test_earth_radius_is_the_gaussian_mean_at_the_mean_latitude(
        self, frame: ModuleType
    ) -> None:
        pts_l, tri_l = [], []
        for k, s in enumerate(("solund", "vardo")):
            p, _ = triangle_at(*site_xy(*SITES[s][:2]))
            pts_l.append(p)
            tri_l.append([3 * k, 3 * k + 1, 3 * k + 2])
        rf = frame.receivers_from_mesh(np.vstack(pts_l), np.array(tri_l), FRAME)
        m_r, n_r = radii(float(np.mean(rf.lat)))
        assert rf.earth_radius_m == pytest.approx(math.sqrt(m_r * n_r), abs=1e-6)

    def test_an_empty_mesh_gives_empty_arrays_of_the_right_shapes(self, frame: ModuleType) -> None:
        rf = frame.receivers_from_mesh(np.empty((0, 3)), np.empty((0, 3), dtype=np.int64), FRAME)
        assert rf.jacobian.shape == (0, 2, 2)
        assert rf.normal.shape == (0, 3)
        for name in ("x", "y", "z_obs", "lat", "lon", "exclusion_radius"):
            assert getattr(rf, name).shape == (0,), name


# --- refusals ----------------------------------------------------------------


class TestRefusals:
    @staticmethod
    def two_triangles(z: float = 0.0) -> tuple[F64, npt.NDArray[np.int64]]:
        x, y = site_xy(*SITES["bygdin"][:2])
        pts = np.array(
            [[x, y, z], [x + 10, y, z], [x, y + 10, z], [x + 10, y + 10, z], [x + 20, y + 20, z]]
        )
        return pts, np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64)

    def test_a_repeated_vertex_is_a_zero_area_triangle(self, frame: ModuleType) -> None:
        pts, tri = self.two_triangles()
        tri[1] = [1, 3, 1]
        with pytest.raises(ValueError, match=r"triangle 1\b.*(zero|area)"):
            frame.receivers_from_mesh(pts, tri, FRAME)

    def test_three_collinear_points_are_a_zero_area_triangle(self, frame: ModuleType) -> None:
        # Collinear in 3-D (so not a vertical triangle): (0,0,0), (10,10,0), (20,20,0).
        pts, tri = self.two_triangles()
        tri[1] = [0, 3, 4]
        with pytest.raises(ValueError, match=r"triangle 1\b.*(zero|area)"):
            frame.receivers_from_mesh(pts, tri, FRAME)

    @pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
    @pytest.mark.parametrize("column", [0, 1, 2])
    def test_a_non_finite_vertex_names_its_triangle(
        self, frame: ModuleType, bad: float, column: int
    ) -> None:
        pts, tri = self.two_triangles()
        pts[3, column] = bad  # vertex 3 is used by triangle 1 only
        with pytest.raises(ValueError, match=r"triangle 1\b.*finite"):
            frame.receivers_from_mesh(pts, tri, FRAME)

    def test_a_vertex_index_out_of_range_names_its_triangle(self, frame: ModuleType) -> None:
        pts, tri = self.two_triangles()
        tri[0, 2] = 5
        with pytest.raises(ValueError, match=r"triangle 0\b"):
            frame.receivers_from_mesh(pts, tri, FRAME)

    def test_a_negative_observer_offset_is_refused(self, frame: ModuleType) -> None:
        with pytest.raises(ValueError, match="offset"):
            frame.receivers_from_mesh(*self.two_triangles(), FRAME, observer_offset=-0.1)

    def test_a_receiver_at_the_pole_is_refused(self, frame: ModuleType) -> None:
        # EPSG:3995 (Arctic polar stereographic) puts the pole at (0, 0).
        pts, tri = triangle_at(0.0, 0.0)
        with pytest.raises(ValueError, match=r"89\.9"):
            frame.receivers_from_mesh(pts, tri, "EPSG:3995")

    def test_a_geographic_frame_is_refused(self, frame: ModuleType) -> None:
        pts = np.array([[8.0, 61.0, 0.0], [8.001, 61.0, 0.0], [8.0, 61.001, 0.0]])
        with pytest.raises(ValueError, match="projected"):
            frame.receivers_from_mesh(pts, np.array([[0, 1, 2]]), "EPSG:4326")


# --- the core's check, through the binding ----------------------------------


class TestCoreCheck:
    @staticmethod
    def arrays(frame: ModuleType) -> dict[str, F64]:
        sites = [SITES[s][:2] for s in ("solund", "oslo", "bygdin")]
        pts_l, tri_l = [], []
        for k, (lon, lat) in enumerate(sites):
            p, _ = triangle_at(*site_xy(lon, lat), (0.0, 5.0, 9.0))
            pts_l.append(p)
            tri_l.append([3 * k, 3 * k + 1, 3 * k + 2])
        rf = frame.receivers_from_mesh(np.vstack(pts_l), np.array(tri_l), FRAME)
        return {
            "x": rf.x,
            "y": rf.y,
            "z_obs": rf.z_obs,
            "lat_deg": rf.lat,
            "lon_deg": rf.lon,
            "jacobian": rf.jacobian,
            "normal": rf.normal,
            "exclusion_radius": rf.exclusion_radius,
        }

    def test_the_frame_passes_the_core_check(self, frame: ModuleType, core: Any) -> None:
        assert core.check_receivers(**self.arrays(frame)) is None

    def test_a_nan_names_its_receiver(self, frame: ModuleType, core: Any) -> None:
        a = self.arrays(frame)
        a["normal"] = a["normal"].copy()
        a["normal"][2, 1] = math.nan
        with pytest.raises(ValueError, match=r"receiver 2\b.*finite"):
            core.check_receivers(**a)

    def test_a_polar_latitude_names_its_receiver(self, frame: ModuleType, core: Any) -> None:
        a = self.arrays(frame)
        a["lat_deg"] = a["lat_deg"].copy()
        a["lat_deg"][1] = 89.95
        with pytest.raises(ValueError, match=r"receiver 1\b.*latitude"):
            core.check_receivers(**a)

    def test_a_negative_exclusion_radius_names_its_receiver(
        self, frame: ModuleType, core: Any
    ) -> None:
        a = self.arrays(frame)
        a["exclusion_radius"] = a["exclusion_radius"].copy()
        a["exclusion_radius"][0] = -1.0
        with pytest.raises(ValueError, match=r"receiver 0\b.*exclusion radius"):
            core.check_receivers(**a)

    def test_arrays_of_unequal_length_are_refused(self, frame: ModuleType, core: Any) -> None:
        a = self.arrays(frame)
        a["z_obs"] = a["z_obs"][:2]
        with pytest.raises(ValueError, match="length"):
            core.check_receivers(**a)


# --- RadiationRequest --------------------------------------------------------


class TestRadiationRequest:
    def test_defaults(self, request_mod: ModuleType) -> None:
        r = request_mod.RadiationRequest()
        assert r.n_directions == 64
        assert r.refraction_k == pytest.approx(0.13)
        assert r.d_near_m == pytest.approx(20_000.0)
        assert r.d_far_m == pytest.approx(200_000.0)
        assert r.limb == "centre"

    def test_is_frozen(self, request_mod: ModuleType) -> None:
        import pydantic

        r = request_mod.RadiationRequest()
        with pytest.raises(pydantic.ValidationError):
            r.n_directions = 32

    def test_the_far_reach_may_equal_the_near_one(self, request_mod: ModuleType) -> None:
        r = request_mod.RadiationRequest(d_near_m=20_000.0, d_far_m=20_000.0, limb="upper")
        assert r.d_far_m == r.d_near_m and r.limb == "upper"

    @pytest.mark.parametrize(
        ("kwargs", "field"),
        [
            ({"n_directions": 48}, "n_directions"),
            ({"n_directions": 2}, "n_directions"),
            ({"n_directions": 0}, "n_directions"),
            ({"refraction_k": 1.0}, "refraction_k"),
            ({"refraction_k": -0.01}, "refraction_k"),
            ({"refraction_k": math.nan}, "refraction_k"),
            ({"d_near_m": 0.0}, "d_near_m"),
            ({"d_near_m": math.inf}, "d_near_m"),
            ({"d_near_m": 30_000.0, "d_far_m": 20_000.0}, "d_far_m"),
            ({"limb": "lower"}, "limb"),
            ({"horizon": 3}, "horizon"),
        ],
    )
    def test_refusals_name_the_field(
        self, request_mod: ModuleType, kwargs: dict[str, Any], field: str
    ) -> None:
        import pydantic

        with pytest.raises(pydantic.ValidationError, match=field):
            request_mod.RadiationRequest(**kwargs)
