"""Bindings for the C++20 geometry primitives in include/terrain/core/point.hpp."""

import math

import pytest

from tin_engine import Point2, Point3, cross, dot


class TestPoint2:
    def test_constructs_and_exposes_coordinates(self) -> None:
        p = Point2(1.5, -2.5)
        assert p.x == 1.5
        assert p.y == -2.5

    def test_defaults_to_origin(self) -> None:
        p = Point2()
        assert p.x == 0.0
        assert p.y == 0.0

    def test_equality_is_by_value(self) -> None:
        assert Point2(1.0, 2.0) == Point2(1.0, 2.0)
        assert Point2(1.0, 2.0) != Point2(1.0, 2.5)

    def test_arithmetic(self) -> None:
        a, b = Point2(1.0, 2.0), Point2(0.5, -1.0)
        assert a + b == Point2(1.5, 1.0)
        assert a - b == Point2(0.5, 3.0)
        assert a * 2.0 == Point2(2.0, 4.0)
        assert 2.0 * a == Point2(2.0, 4.0)
        assert a / 2.0 == Point2(0.5, 1.0)
        assert -a == Point2(-1.0, -2.0)

    def test_dot_and_cross(self) -> None:
        a, b = Point2(1.0, 2.0), Point2(3.0, 4.0)
        assert dot(a, b) == pytest.approx(11.0)
        # 2D cross is the scalar z-component.
        assert cross(a, b) == pytest.approx(-2.0)

    def test_is_hashable_by_value(self) -> None:
        assert len({Point2(1.0, 2.0), Point2(1.0, 2.0)}) == 1
        assert len({Point2(1.0, 2.0), Point2(2.0, 1.0)}) == 2

    def test_repr_round_trips_through_std_format(self) -> None:
        assert repr(Point2(1.0, 2.0)) == "Point2(1, 2)"


class TestPoint3:
    def test_constructs_and_exposes_coordinates(self) -> None:
        p = Point3(1.0, 2.0, 3.0)
        assert (p.x, p.y, p.z) == (1.0, 2.0, 3.0)

    def test_defaults_to_origin(self) -> None:
        assert Point3() == Point3(0.0, 0.0, 0.0)

    def test_arithmetic(self) -> None:
        a, b = Point3(1.0, 2.0, 3.0), Point3(0.5, -1.0, 2.0)
        assert a + b == Point3(1.5, 1.0, 5.0)
        assert a - b == Point3(0.5, 3.0, 1.0)
        assert a * 2.0 == Point3(2.0, 4.0, 6.0)
        assert -a == Point3(-1.0, -2.0, -3.0)

    def test_dot(self) -> None:
        assert dot(Point3(1.0, 2.0, 3.0), Point3(4.0, 5.0, 6.0)) == pytest.approx(32.0)

    def test_cross_is_a_vector_and_orthogonal_to_both(self) -> None:
        a, b = Point3(1.0, 0.0, 0.0), Point3(0.0, 1.0, 0.0)
        c = cross(a, b)
        assert c == Point3(0.0, 0.0, 1.0)
        assert dot(c, a) == pytest.approx(0.0)
        assert dot(c, b) == pytest.approx(0.0)

    def test_is_hashable_by_value(self) -> None:
        assert len({Point3(1.0, 2.0, 3.0), Point3(1.0, 2.0, 3.0)}) == 1

    def test_repr_round_trips_through_std_format(self) -> None:
        assert repr(Point3(1.0, 2.0, 3.0)) == "Point3(1, 2, 3)"


def test_dot_rejects_mixed_dimensions() -> None:
    with pytest.raises(TypeError):
        dot(Point2(1.0, 2.0), Point3(1.0, 2.0, 3.0))


def test_nan_never_compares_equal() -> None:
    nan = float("nan")
    assert Point2(nan, 0.0) != Point2(nan, 0.0)
    assert math.isnan(Point2(nan, 0.0).x)


class TestDocumentation:
    """CLAUDE.md section 4 requires every pybind11-exposed surface to be documented."""

    def test_types_and_functions_carry_docstrings(self) -> None:
        for obj in (Point2, Point3, dot, cross):
            assert obj.__doc__, f"{obj!r} has no docstring"

    def test_properties_carry_docstrings(self) -> None:
        for cls, names in ((Point2, "xy"), (Point3, "xyz")):
            for name in names:
                assert getattr(cls, name).__doc__, f"{cls.__name__}.{name} undocumented"

    def test_cross_documents_its_dimensional_asymmetry(self) -> None:
        # 2D cross returns a scalar, 3D returns a vector. That asymmetry was
        # only ever explained in a C++ comment, invisible from help(cross).
        assert "scalar" in cross.__doc__


class TestImmutability:
    """Point types are hashable, so they must not be mutable.

    Exposing def_readwrite alongside __hash__ lets a caller insert a point into
    a set, mutate a coordinate, and silently corrupt the set -- the element
    becomes unreachable. C++ aggregates get away with this; Python sets do not.
    """

    def test_coordinates_are_read_only(self) -> None:
        p = Point2(1.0, 2.0)
        with pytest.raises(AttributeError):
            p.x = 5.0  # type: ignore[misc]

        q = Point3(1.0, 2.0, 3.0)
        with pytest.raises(AttributeError):
            q.z = 5.0  # type: ignore[misc]

    def test_set_membership_survives_attempted_mutation(self) -> None:
        p = Point2(1.0, 2.0)
        s = {p}
        with pytest.raises(AttributeError):
            p.x = 99.0  # type: ignore[misc]
        assert p in s
