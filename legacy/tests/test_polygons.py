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
def circle_polygon():
    N = 117
    x, y = array([0.5, 0.5])
    r = 0.4

    return td.simple_polygon([(x + r*cos(t), y + r*sin(t))
                              for t in linspace(0, 2*pi, N+1)[:-1]])


@pytest.mark.skip # Segfaults
def test_intersection(rectangle_polygon, circle_polygon):
    polygon = rectangle_polygon.intersection(circle_polygon)
    assert polygon.num_parts() == 1
