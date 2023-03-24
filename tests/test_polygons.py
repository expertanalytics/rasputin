from numpy import array, cos, sin, linspace, pi

import pytest

import rasputin.triangulate_dem as td

@pytest.fixture
def rectangle_polygon():
    return td.simple_polygon([[0.48, 0.0],
                              [1.0, 0.0],
                              [1.0, 1.0],
                              [0.48, 1.0]])
@pytest.fixture
def square_polygon():
    return td.simple_polygon([[0.48, 0.0],
                              [0.96, 0.0],
                              [0.96, 0.96],
                              [0.48, 0.96]])

@pytest.fixture
def circle_polygon():
    N = 117
    x, y = array([0.5, 0.5])
    r = 0.4

    return td.simple_polygon([(x + r*cos(t), y + r*sin(t))
                              for t in linspace(0, 2*pi, N+1)[:-1]])

@pytest.mark.skip # Segfaults
def test_intersection_rect_circ(rectangle_polygon, circle_polygon):
    polygon = rectangle_polygon.intersection(circle_polygon)
    assert polygon.num_parts() == 1

@pytest.mark.skip # Segfaults
def test_intersection_square_circ(square_polygon, circle_polygon):
    polygon = square_polygon.intersection(circle_polygon)
    assert len(polygon) == 1

def test_intersection_rect_square(rectangle_polygon, square_polygon):
    polygon = rectangle_polygon.intersection(square_polygon)
    assert isinstance(polygon, td.multi_polygon) and polygon.num_parts() == 1
    assert isinstance(polygon[0], td.polygon)
    assert len(polygon[0].holes()) == 0
    assert isinstance(polygon[0].exterior(), td.simple_polygon)

    polygon = square_polygon.intersection(rectangle_polygon)
    assert isinstance(polygon, td.multi_polygon) and polygon.num_parts() == 1
    assert isinstance(polygon[0], td.polygon)
    assert len(polygon[0].holes()) == 0
    assert isinstance(polygon[0].exterior(), td.simple_polygon)

def test_difference_rect_square(rectangle_polygon, square_polygon):
    polygon = rectangle_polygon.difference(square_polygon)
    assert isinstance(polygon, td.multi_polygon) and polygon.num_parts() == 1
    assert isinstance(polygon[0], td.polygon)
    assert len(polygon[0].holes()) == 0
    assert isinstance(polygon[0].exterior(), td.simple_polygon)

    polygon = square_polygon.difference(rectangle_polygon)
    assert isinstance(polygon, td.multi_polygon) and polygon.num_parts() == 0
