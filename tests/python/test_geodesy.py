"""`geodesy`: metres near a point, cell sizes and geodesic areas (increment 40).

`docs/increments/40-geographic-catchment.md`, sections 5.2 to 5.4 and the red
suites R4 and R5. The oracles are `geodesy_fixtures`' closed forms from the
ellipsoid's `a` and `f` (not pyproj, not the module), and pyproj's `Geod`
where the design names it (`area_m2`).

Interface as the design gives it (5.2, 5.4): `LocalFrame(lon0, lat0, kx, ky)`,
frozen, with `identity()`, `at(crs, point)`, `is_identity`, `forward`,
`inverse`, `units(metres) -> (mx, my)`, `cell_m(meta)` and `meta(meta)`;
`cell_areas_by_row(meta)`, `mean_cell_area(mask, meta)` and
`area_m2(geometry, crs)`.

Scales: the figures in metres per degree are checked to 1e-4 m on values of
1e5 (float64 carries about 1e-11 there); degrees to 1e-12 on values below
0.04; areas relative 1e-6 on a 1" cell (480 to 950 m²) and 1e-9 on rings of a
few km². The largest coordinate used is 6.6e6 m (UTM northing) and 60°.
"""

from __future__ import annotations

import dataclasses
from types import ModuleType

import numpy as np
import pytest
from pyproj import Geod
from shapely.geometry import Polygon, box

from geodesy_fixtures import (
    ARCSEC,
    INTL_A,
    INTL_F,
    band_area,
    geo_meta,
    metres_per_degree,
)
from mosaic_fixtures import meta as projected_meta

#: A geographic CRS whose axes are in grads (5.2: refused in plain words).
GRAD_CRS = (
    'GEOGCRS["WGS 84 in grads",DATUM["World Geodetic System 1984",'
    'ELLIPSOID["WGS 84",6378137,298.257223563]],CS[ellipsoidal,2],'
    'AXIS["longitude",east,ANGLEUNIT["grad",0.0157079632679489]],'
    'AXIS["latitude",north,ANGLEUNIT["grad",0.0157079632679489]]]'
)


@pytest.fixture(scope="module")
def geodesy() -> ModuleType:
    import tin_engine.geodesy as module

    return module


def arcsec_grid(lat: float, rows: int = 5, cols: int = 4) -> object:
    """A 1" grid in EPSG:4326 whose middle row lies on `lat`."""
    return geo_meta(rows, cols, lon0=10.0, lat0=lat + (rows // 2) * ARCSEC, dx=ARCSEC, dy=ARCSEC)


# ---------------------------------------------------------------- R5: LocalFrame


def test_the_identity_frame(geodesy: ModuleType) -> None:
    frame = geodesy.LocalFrame.identity()
    assert (frame.lon0, frame.lat0, frame.kx, frame.ky) == (0.0, 0.0, 1.0, 1.0)
    assert frame.is_identity
    with pytest.raises(dataclasses.FrozenInstanceError):
        frame.kx = 2.0  # type: ignore[misc]


def test_a_projected_crs_gives_the_identity_bit_for_bit(geodesy: ModuleType) -> None:
    """Every projected path computes what it does today (5.2's second rule)."""
    frame = geodesy.LocalFrame.at("EPSG:25833", (500_000.0, 6_600_000.0))
    assert frame.is_identity
    assert frame.units(2000.0) == (2000.0, 2000.0)
    assert frame.units(30.0 + 2000.0) == (2030.0, 2030.0)
    xy = np.array([[500_000.1, 6_600_000.3], [512_345.678, 6_543_210.987]])
    assert np.array_equal(frame.forward(xy), xy)
    assert np.array_equal(frame.inverse(xy), xy)
    m = projected_meta(rows=3, cols=4)
    assert frame.meta(m) is m
    assert frame.cell_m(m) == max(m.delta_x, m.delta_y)


@pytest.mark.parametrize(
    ("lat", "kx", "ky", "units"),
    [
        (0.0, 111_319.4908, 110_574.2758, (0.01796631, 0.01808739)),
        (60.0, 55_800.0016, 111_412.2875, (0.03584229, 0.01795134)),
    ],
    ids=["equator", "60N"],
)
def test_metres_per_degree_and_units_at_two_latitudes(
    geodesy: ModuleType, lat: float, kx: float, ky: float, units: tuple[float, float]
) -> None:
    """The design's checked figures (5.2, R5) and the closed form agree."""
    frame = geodesy.LocalFrame.at("EPSG:4326", (10.0, lat))
    assert not frame.is_identity
    assert (frame.lon0, frame.lat0) == (10.0, lat)
    assert frame.kx == pytest.approx(kx, abs=1e-4)
    assert frame.ky == pytest.approx(ky, abs=1e-4)
    oracle = metres_per_degree(lat)
    assert (frame.kx, frame.ky) == pytest.approx(oracle, rel=1e-12)
    mx, my = frame.units(2000.0)
    assert (mx, my) == pytest.approx(units, abs=5e-9)
    assert (mx, my) == pytest.approx((2000.0 / oracle[0], 2000.0 / oracle[1]), rel=1e-12)


def test_the_ellipsoid_is_the_crs_own_not_wgs84(geodesy: ModuleType) -> None:
    """ED50 (International 1924) at the equator: `kx` is its `a` times pi/180,
    111,323.8716 m, not WGS 84's 111,319.49."""
    frame = geodesy.LocalFrame.at("EPSG:4230", (5.0, 0.0))
    assert (frame.kx, frame.ky) == pytest.approx(metres_per_degree(0.0, INTL_A, INTL_F), rel=1e-12)
    assert frame.kx == pytest.approx(111_323.8716, abs=1e-4)


def test_a_geographic_crs_in_grads_is_refused(geodesy: ModuleType) -> None:
    with pytest.raises(ValueError, match=r"(?i)degree"):
        geodesy.LocalFrame.at(GRAD_CRS, (10.0, 60.0))


@pytest.mark.parametrize("lat", [0.0, 60.0, -19.4436])
def test_forward_maps_to_metres_and_inverse_returns(geodesy: ModuleType, lat: float) -> None:
    """X = (lon - lon0) kx, Y = (lat - lat0) ky; back to 1e-9 degrees."""
    frame = geodesy.LocalFrame.at("EPSG:4326", (-44.0, lat))
    lonlat = np.array([[-44.0, lat], [-44.01, lat + 0.02], [-43.97, lat - 0.013]])
    xy = frame.forward(lonlat)
    assert xy.shape == (3, 2) and xy.dtype == np.float64
    expected = np.column_stack([(lonlat[:, 0] + 44.0) * frame.kx, (lonlat[:, 1] - lat) * frame.ky])
    assert np.allclose(xy, expected, rtol=0.0, atol=1e-6)
    assert xy[0].tolist() == [0.0, 0.0]
    assert np.allclose(frame.inverse(xy), lonlat, rtol=0.0, atol=1e-9)
    # A tuple of tuples is array-like too (`reach.line`).
    assert np.allclose(frame.forward(tuple(map(tuple, lonlat))), xy, rtol=0.0, atol=0.0)


@pytest.mark.parametrize(("lat", "cell"), [(0.0, 30.9221), (60.0, 30.9479)], ids=["equator", "60N"])
def test_a_cell_is_its_longer_side_in_metres(geodesy: ModuleType, lat: float, cell: float) -> None:
    """At 0° the east-west side (30.92 m against 30.72) is the longer; at 60°
    the north-south one (30.95 m against 15.50). Section 5.3, question Q4."""
    frame = geodesy.LocalFrame.at("EPSG:4326", (10.0, lat))
    assert frame.cell_m(arcsec_grid(lat)) == pytest.approx(cell, abs=1e-4)


def test_the_anadem_cell_at_the_fixtures_centre(geodesy: ModuleType) -> None:
    """Section 3: 28.2996 m by 29.8324 m at 19.4436° S on SIRGAS 2000."""
    step = 0.00026949458523585647
    m = geo_meta(3, 3, lon0=-44.1, lat0=-19.4436 + step, dx=step, dy=step, epsg=4674)
    frame = geodesy.LocalFrame.at("EPSG:4674", (-44.1, -19.4436))
    assert step * frame.kx == pytest.approx(28.2996, abs=1e-4)
    assert step * frame.ky == pytest.approx(29.8324, abs=1e-4)
    assert frame.cell_m(m) == pytest.approx(29.8324, abs=1e-4)


def test_meta_maps_the_node_grid_into_the_frame(geodesy: ModuleType) -> None:
    """For the burn (5.5): x_min' = (x_min - lon0) kx, delta_x' = delta_x kx,
    likewise in y; not geographic; the shape and NoData unchanged."""
    m = geo_meta(7, 9, lon0=10.0, lat0=60.002, dx=ARCSEC, dy=ARCSEC, nodata=-9999.0)
    frame = geodesy.LocalFrame.at("EPSG:4326", (10.001, 60.0))
    out = frame.meta(m)
    assert out.x_min == pytest.approx((10.0 - 10.001) * frame.kx, rel=1e-12)
    assert out.y_max == pytest.approx((60.002 - 60.0) * frame.ky, rel=1e-12)
    assert out.delta_x == pytest.approx(ARCSEC * frame.kx, rel=1e-12)
    assert out.delta_y == pytest.approx(ARCSEC * frame.ky, rel=1e-12)
    assert (out.rows, out.cols, out.nodata) == (7, 9, -9999.0)
    assert out.geographic is False


# ---------------------------------------------------------------- R4: areas


@pytest.mark.parametrize(
    ("lat", "hand"), [(0.0, 949.7741), (60.0, 479.6918)], ids=["equator", "60N"]
)
def test_a_cell_is_centred_on_its_node(geodesy: ModuleType, lat: float, hand: float) -> None:
    """R4's hand figures: a 1" cell from half a spacing south of its node to
    half a spacing north. The band 60°00'00" to 60°00'01" is 479.6898 m², 4e-6
    off, so a cell with the wrong edges fails here."""
    m = arcsec_grid(lat)
    areas = geodesy.cell_areas_by_row(m)
    assert areas.shape == (m.rows,) and areas.dtype == np.float64
    middle = m.rows // 2
    assert areas[middle] == pytest.approx(hand, rel=1e-6)
    for r in range(m.rows):
        y = m.y_max - r * m.delta_y
        expected = band_area(y - ARCSEC / 2, y + ARCSEC / 2, ARCSEC)
        assert areas[r] == pytest.approx(expected, rel=1e-6)
    if lat == 60.0:  # at the equator the two bands agree to second order
        assert band_area(lat, lat + ARCSEC, ARCSEC) != pytest.approx(areas[middle], rel=1e-6)


def test_cell_areas_read_the_crs_ellipsoid(geodesy: ModuleType) -> None:
    m = geo_meta(3, 3, lon0=5.0, lat0=ARCSEC, dx=ARCSEC, dy=ARCSEC, epsg=4230)
    areas = geodesy.cell_areas_by_row(m)
    expected = band_area(-ARCSEC / 2, ARCSEC / 2, ARCSEC, INTL_A, INTL_F)
    assert areas[1] == pytest.approx(expected, rel=1e-6)


def test_a_projected_grid_has_dx_dy_in_every_row(geodesy: ModuleType) -> None:
    m = projected_meta(rows=6, cols=3, dx=10.0, dy=5.0)
    areas = geodesy.cell_areas_by_row(m)
    assert areas.tolist() == [50.0] * 6


def test_mean_cell_area_on_a_projected_grid_is_dx_dy(geodesy: ModuleType) -> None:
    m = projected_meta(rows=4, cols=4, dx=10.0, dy=5.0)
    mask = np.zeros((4, 4), dtype=np.uint8)
    mask[1:3, 2] = 1
    assert geodesy.mean_cell_area(mask, m) == 50.0


def test_mean_cell_area_in_degrees_weights_the_rows(geodesy: ModuleType) -> None:
    """Three in-nodes on the top row, one on the bottom: (3 a0 + a4) / 4."""
    m = arcsec_grid(60.0)
    mask = np.zeros((m.rows, m.cols), dtype=np.uint8)
    mask[0, :3] = 1
    mask[m.rows - 1, 1] = 1
    rows = [m.y_max - r * m.delta_y for r in (0, m.rows - 1)]
    a0, a4 = (band_area(y - ARCSEC / 2, y + ARCSEC / 2, ARCSEC) for y in rows)
    assert geodesy.mean_cell_area(mask, m) == pytest.approx((3 * a0 + a4) / 4, rel=1e-9)


def test_area_m2_on_a_projected_crs_is_the_planar_area(geodesy: ModuleType) -> None:
    ring = Polygon([(500_000.0, 6_600_000.0), (503_000.0, 6_600_000.0), (500_000.0, 6_601_700.0)])
    assert geodesy.area_m2(ring, "EPSG:25833") == ring.area


@pytest.mark.parametrize(("crs", "ellps"), [("EPSG:4326", "WGS84"), ("EPSG:4230", "intl")])
def test_area_m2_in_degrees_is_geodesic_on_the_crs_ellipsoid(
    geodesy: ModuleType, crs: str, ellps: str
) -> None:
    """|Geod.geometry_area_perimeter|, either orientation of the ring."""
    cw = Polygon([(10.0, 60.0), (10.0, 60.02), (10.05, 60.02), (10.05, 60.0)])
    expected = abs(Geod(ellps=ellps).geometry_area_perimeter(cw)[0])
    assert geodesy.area_m2(cw, crs) == pytest.approx(expected, rel=1e-9)
    ccw = box(10.0, 60.0, 10.05, 60.02)
    assert geodesy.area_m2(ccw, crs) == pytest.approx(expected, rel=1e-9)
    assert expected > 0.0
