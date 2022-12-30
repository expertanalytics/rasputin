import pytest
from numpy import array, cos, sin, linspace, pi, float32, zeros, ndindex, sqrt, meshgrid, int64
from numpy.linalg import norm
from datetime import datetime, timedelta

import pyproj
from shapely.geometry import Polygon, Point

from rasputin.geometry import GeoPolygon
from rasputin.reader import Rasterdata
from rasputin.mesh import Mesh
from rasputin import triangulate_dem


@pytest.fixture
def raster_xm():
    # Place origo in Oslo
    lat, lon = 60, 10
    dx = 100
    dy = 100
    epsgid = 32633
    cart_crs = pyproj.CRS.from_epsg(epsgid)
    geo_crs = pyproj.CRS.from_epsg(4326)
    proj = pyproj.Transformer.from_crs(geo_crs, cart_crs)
    x_0, y_0 = proj.transform(lat, lon)
    m, n = 4, 4
    x = linspace(0, 1, m) - 0.5
    y = linspace(0, 1, n) - 0.5
    x, y = meshgrid(x, y)
    r = sqrt(x**2 + y**2)
    arr = -sin(2*pi*r/n)
    return Rasterdata(shape=(m, n),
                      x_min=x_0,
                      y_max=y_0,
                      delta_x=dx,
                      delta_y=dy,
                      array=arr.astype(float32),
                      coordinate_system=f"+init=epsg:{epsgid}",
                      info={})

@pytest.fixture
def raster():
    m, n = 21, 11
    array = (zeros((m, n))
             + linspace(0, 1, m).reshape(-1,  1)**2
             + linspace(0, 1, n).reshape( 1, -1))
    print(array.shape)
    print(f"delta_x={1/(n - 1)}")
    print(f"delta_y={1/(m - 1)}")
    return Rasterdata(shape=(m, n),
                      x_min=0,
                      y_max=1,
                      delta_x=1/(n - 1),
                      delta_y=1/(m - 1),
                      array=array.astype(float32),
                      coordinate_system="+init=epsg:32633",
                      info={})

@pytest.fixture
def raster_list():
    m0, n0, m1, n1 = 31, 20, 53, 10
    d0, d1 = 0.57, 0.48
    array0 = (zeros((m0, n0))
              + linspace(0, d0, m0).reshape(-1,  1)**2
              + linspace(0,  1, n0).reshape( 1, -1))

    array1 = (zeros((m1, n1))
              + linspace(d1, 1, m1).reshape(-1, 1)**2
              + linspace(0,  1, n1).reshape(1, -1))
    return [Rasterdata(shape=(m0, n0),
                       x_min=0,
                       y_max=1,
                       delta_x=d0/(m0 - 1),
                       delta_y=1/(n0 - 1),
                       array=array0.astype(float32),
                       coordinate_system="+init=epsg:32633",
                       info={}),
            Rasterdata(shape=(m0, n0),
                       x_min=d1,
                       y_max=1,
                       delta_x=(1 - d1)/(m1 - 1),
                       delta_y=1/(n1 - 1),
                       array=array1.astype(float32),
                       coordinate_system="+init=epsg:32633",
                       info={})]

@pytest.fixture
def polygon():
    # Create polygonal approximation to a circle
    N = 117
    x, y = array([0.51, 0.5])
    r = 0.4

    epsg_id = 32633

    polygon = Polygon([(x + r*cos(t), y + r*sin(t))
                      for t in linspace(0, 2*pi, N+1)[:-1]])

    return GeoPolygon(polygon=polygon,
                      crs=pyproj.CRS.from_epsg(epsg_id))

@pytest.fixture
def polygon_w_hole():
    # Create polygonal approximation to a circle
    N1, N2 = 117, 77
    x, y = array([0.51, 0.5])
    r1, r2 = 0.4, 0.25

    epsg_id = 32633

    polygon_1 = Polygon((x + r1*cos(t), y + r1*sin(t))
                        for t in linspace(0, 2*pi, N1 + 1)[:-1])

    polygon_2 = Polygon((x + r2*cos(t), y + r2*sin(t))
                        for t in linspace(0, 2*pi, N2 + 1)[:-1])

    return GeoPolygon(polygon=polygon_1.difference(polygon_2),
                      crs=pyproj.CRS.from_epsg(epsg_id))

def test_mesh(raster, polygon):
    mesh = Mesh.from_raster(data=raster, domain=polygon)

    # Check mesh consistency
    assert len(mesh.points) == mesh.num_points
    assert len(mesh.faces) == mesh.num_faces
    assert mesh.characteristic == 1

    # Check that all points in mesh are inside the polygon
    test_poly = polygon.polygon.buffer(1e-10) # Account for fixed float precision
    for (x, y, _) in mesh.points:
        assert test_poly.contains(Point(x, y))

    # Check that all points in raster not inside polygon are not mesh
    def dist(mesh, pt):
        return norm(mesh.points[:, :2] - pt, axis=0).min()

    for i, j in ndindex(raster.array.shape):
        x = raster.x_min + i*raster.delta_x
        y = raster.y_max - j*raster.delta_y
        assert test_poly.contains(Point(x, y)) or dist(mesh, (x, y)) > 1e-10

    # Some sanity tests
    x, y, z = mesh.points.T
    print(raster.x_max)
    assert raster.x_min <= x.min()
    assert raster.x_max >= x.max()
    assert raster.y_min <= y.min()
    assert raster.y_max >= y.max()
    assert raster.array.min() <= z.min()
    assert raster.array.max() >= z.max()


def test_mesh_sub_mesh(raster, polygon):
    mesh = Mesh.from_raster(data=raster, domain=polygon)
    sub_mesh = mesh.extract_sub_mesh(array([0, 1], dtype=int64))
    assert len(sub_mesh.faces) == 2

def test_mesh_w_hole(raster, polygon_w_hole):
    mesh = Mesh.from_raster(data=raster, domain=polygon_w_hole)

    # Check mesh consistency
    assert len(mesh.points) == mesh.num_points
    assert len(mesh.faces) == mesh.num_faces
    assert mesh.characteristic == 0

    # Check that all points in mesh are inside the polygon
    test_poly = polygon_w_hole.polygon.buffer(1e-10)  # Account for fixed float precision
    for (x, y, _) in mesh.points:
        assert test_poly.contains(Point(x, y))

    # Check that all points in raster not inside polygon are not mesh
    def dist(mesh, pt):
        return norm(mesh.points[:,:2] -pt, axis=0).min()

    for i, j in ndindex(raster.array.shape):
        x = raster.x_min + i*raster.delta_x
        y = raster.y_max - j*raster.delta_y
        assert test_poly.contains(Point(x, y)) or dist(mesh, (x, y)) > 1e-10

    # Some sanity tests
    x, y, z = mesh.points.T
    assert raster.x_min <= x.min()
    assert raster.x_max >= x.max()
    assert raster.y_min <= y.min()
    assert raster.y_max >= y.max()
    assert raster.array.min() <= z.min()
    assert raster.array.max() >= z.max()


def test_mesh_shade(raster_xm):
    mesh = Mesh.from_raster(data=raster_xm)
    tp = datetime(2000, 6, 2, 4) + timedelta(minutes=18)
    shades = triangulate_dem.shade(mesh._cpp, tp.timestamp())
    assert len(shades) == mesh.num_faces
