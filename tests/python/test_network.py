"""The network logic of `rasputin network-catchment`: pure, no file.

`docs/increments/85-network-catchment.md`, section 7.1 (the invariant-critical
suite) and the interfaces of section 5. Units are shapely squares of 1 km side
in EPSG:3006 (each 1.000000 km2) unless a test says otherwise, near
x 500 000, y 7 000 000: SWEREF99 TM's central meridian, at 63 degrees north.

The module is imported per test, so its absence fails each test rather than
the collection.

Mutation targets the kill record covers (section 7.1): the upstream walk
(direction, the visited test, the seed included), the cycle refusal, the
duplicate-ID refusal, the network-end count, the outlet radius comparison
(killed in `test_cli_network_catchment.py`, since section 5 gives the
comparison no pure function), the area tolerance comparison, the multipart
refusal, and `unit_id`'s integral-float rule.

Numeric bounds: coordinates are up to 7.0e6 m, where a double's spacing is
about 1e-9 m, so a distance in metres is compared to 1e-6 m and an area in m2
to 1e-9 relative (checked at the largest set here, nine squares). The 0.1 %
bound between EPSG:3006 and EPSG:4326 covers the projection's scale on its
central meridian, 0.9996 (0.04 % in length, 0.08 % in area).
"""

from __future__ import annotations

import importlib
import re
from types import ModuleType
from typing import Any, cast

import pytest
import shapely
from pyproj import Transformer
from shapely.geometry import Point, Polygon, box

CRS = "EPSG:3006"
X0 = 500_000.0
Y0 = 7_000_000.0
KM = 1000.0
KM2 = 1e6
END = "-9999"


def nw() -> ModuleType:
    return importlib.import_module("tin_engine.network")


def square(i: float, j: float, side: float = KM) -> Polygon:
    """The square of `side` metres whose south-west corner is `i`, `j` km from (X0, Y0)."""
    x, y = X0 + i * KM, Y0 + j * KM
    return box(x, y, x + side, y + side)


def link(own: Any, down: str = END, area_m2: float | None = None, name: str | None = None) -> Any:
    return nw().Link(id=own, down=down, upstream_area_m2=area_m2, name=name)


def network(*pairs: tuple[Any, str]) -> Any:
    return nw().Network.build([link(own, down) for own, down in pairs])


def to_4326(geometry: Any) -> Any:
    to = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
    return shapely.transform(geometry, lambda x, y: to.transform(x, y), interleaved=False)


def refusal(call: Any, *args: Any) -> str:
    with pytest.raises(nw().NetworkError) as caught:
        call(*args)
    return str(caught.value)


# ---------------------------------------------------------------------------
# The upstream walk (tests 1 to 3)
# ---------------------------------------------------------------------------


class TestUpstream:
    def test_a_chain_from_each_unit_includes_the_unit_and_only_what_is_above(self) -> None:
        """Test 1: C -> B -> A -> end."""
        chain = network(("A", END), ("B", "A"), ("C", "B"))
        assert chain.upstream("A") == ("A", "B", "C")
        assert chain.upstream("B") == ("B", "C")
        assert chain.upstream("C") == ("C",)

    def test_the_end_of_a_chain_is_counted_by_its_value(self) -> None:
        chain = network(("A", END), ("B", "A"), ("C", "B"))
        assert dict(chain.ends) == {END: 1}

    def test_ends_are_counted_per_value_and_a_linked_down_value_is_no_end(self) -> None:
        net = network(("A", END), ("B", END), ("C", "0"), ("D", "A"), ("E", "D"))
        assert dict(net.ends) == {END: 2, "0": 1}

    def test_a_confluence_collects_every_branch_breadth_first(self) -> None:
        """Test 2: B and C drain to A, D to C."""
        net = network(("A", END), ("B", "A"), ("C", "A"), ("D", "C"))
        found = net.upstream("A")
        assert found[0] == "A"
        assert sorted(found) == ["A", "B", "C", "D"]  # each unit once
        assert found[-1] == "D"  # two links above A, after the units one link above
        assert net.upstream("C") == ("C", "D")

    def test_a_cycle_is_refused_naming_its_units_only(self) -> None:
        """Test 3: unitB -> unitC -> unitB, unitA draining into the cycle."""
        net = network(("unitA", "unitB"), ("unitB", "unitC"), ("unitC", "unitB"), ("unitD", END))
        said = refusal(net.upstream, "unitB")
        assert "unitB" in said and "unitC" in said, said
        assert "unitA" not in said, said  # it drains into the cycle; it is not on it
        assert "unitC" in refusal(net.upstream, "unitC")

    def test_units_off_a_cycle_still_walk(self) -> None:
        net = network(("unitA", "unitB"), ("unitB", "unitC"), ("unitC", "unitB"), ("unitD", END))
        assert net.upstream("unitD") == ("unitD",)
        assert net.upstream("unitA") == ("unitA",)

    def test_a_unit_draining_to_itself_is_a_cycle(self) -> None:
        net = network(("loop7", "loop7"), ("feeder", "loop7"))
        assert "loop7" in refusal(net.upstream, "loop7")


# ---------------------------------------------------------------------------
# Building the network (test 4)
# ---------------------------------------------------------------------------


class TestBuild:
    def test_a_duplicate_own_id_is_refused_naming_it(self) -> None:
        said = refusal(nw().Network.build, [link("dup42"), link("other"), link("dup42", "other")])
        assert "dup42" in said, said

    def test_an_empty_own_id_is_refused(self) -> None:
        refusal(nw().Network.build, [link("A"), link("", "A")])

    def test_a_null_own_id_is_refused(self) -> None:
        refusal(nw().Network.build, [link("A"), link(cast(str, None), "A")])

    def test_links_are_kept_by_own_id(self) -> None:
        built = nw().Network.build([link("A", area_m2=2 * KM2, name="Alpha"), link("B", "A")])
        assert set(built.links) == {"A", "B"}
        assert built.links["A"].upstream_area_m2 == 2 * KM2
        assert built.links["A"].name == "Alpha"
        assert built.links["B"].down == "A"


# ---------------------------------------------------------------------------
# IDs compared as text (test 5)
# ---------------------------------------------------------------------------


class TestUnitId:
    @pytest.mark.parametrize("value", [7, 7.0, "7"])
    def test_an_integral_value_is_written_as_an_integer(self, value: object) -> None:
        assert nw().unit_id(value) == "7"

    def test_a_fraction_is_kept(self) -> None:
        assert nw().unit_id(7.5) == "7.5"

    @pytest.mark.parametrize(
        ("value", "text"),
        [(2120008490.0, "2120008490"), (2120008490, "2120008490"), (-9999.0, "-9999")],
    )
    def test_hydrobasins_ids_and_the_svar_end(self, value: object, text: str) -> None:
        assert nw().unit_id(value) == text

    def test_a_text_id_is_unchanged(self) -> None:
        uuid = "11B6DEC8-0A79-40FE-AF91-7413AE9B4703"
        assert nw().unit_id(uuid) == uuid


# ---------------------------------------------------------------------------
# The join (test 6)
# ---------------------------------------------------------------------------


class TestJoin:
    def test_a_row_of_two_squares(self) -> None:
        joined = nw().join([square(0, 0), square(1, 0)], CRS)
        assert joined.polygon.geom_type == "Polygon"
        assert joined.area_m2 == pytest.approx(2 * KM2, rel=1e-9)
        assert joined.holes == 0
        assert joined.holes_m2 == 0.0

    def test_a_ring_of_eight_round_a_missing_ninth(self) -> None:
        ring = [square(i, j) for i in range(3) for j in range(3) if (i, j) != (1, 1)]
        joined = nw().join(ring, CRS)
        assert joined.polygon.geom_type == "Polygon"
        assert joined.area_m2 == pytest.approx(8 * KM2, rel=1e-9)  # the hole subtracted
        assert joined.holes == 1
        assert joined.holes_m2 == pytest.approx(KM2, rel=1e-9)

    def test_neighbours_that_share_no_vertex_still_join(self) -> None:
        """A T-junction: one 1 x 2 km unit beside two 1 km squares, its side
        holding no vertex at their shared corner (SVAR's neighbours are noded
        this way, which is why the join is `union_all`)."""
        tall = box(X0 + KM, Y0, X0 + 2 * KM, Y0 + 2 * KM)
        joined = nw().join([square(0, 0), square(0, 1), tall], CRS)
        assert joined.polygon.geom_type == "Polygon"
        assert joined.area_m2 == pytest.approx(4 * KM2, rel=1e-9)
        assert joined.holes == 0

    def test_two_squares_one_metre_apart_are_refused_with_the_parts_areas(self) -> None:
        apart = Polygon(
            [(X0 + KM + 1, Y0), (X0 + 2 * KM + 1, Y0), (X0 + 2 * KM + 1, Y0 + KM),
             (X0 + KM + 1, Y0 + KM)]
        )  # fmt: skip
        said = refusal(nw().join, [square(0, 0), apart], CRS)
        assert len(re.findall(r"\b1(?:\.0+)? km2", said)) == 2, said

    def test_a_join_in_degrees_has_the_projected_area(self) -> None:
        """Test 9's degrees half, for areas: geodesic in EPSG:4326."""
        joined = nw().join([to_4326(square(0, 0)), to_4326(square(1, 0))], "EPSG:4326")
        assert joined.area_m2 == pytest.approx(2 * KM2, rel=1e-3)


# ---------------------------------------------------------------------------
# The area check (test 7)
# ---------------------------------------------------------------------------


class TestCheckArea:
    @pytest.mark.parametrize("joined_km2", [2.0, 1.981, 2.019])
    def test_within_one_percent_passes(self, joined_km2: float) -> None:
        assert nw().check_area(joined_km2 * KM2, 2.0 * KM2) is None

    @pytest.mark.parametrize("joined_km2", [1.979, 2.021])
    def test_beyond_one_percent_is_refused_naming_both(self, joined_km2: float) -> None:
        said = refusal(nw().check_area, joined_km2 * KM2, 2.0 * KM2)
        assert f"{joined_km2:.3f}" in said, said
        assert "lacks units" in said, said

    def test_the_tolerance_is_relative_to_the_stated_area(self) -> None:
        """0.5 % of a large basin passes; the same absolute gap on a small one does not."""
        assert nw().check_area(995.0 * KM2, 1000.0 * KM2) is None
        refusal(nw().check_area, 0.5 * KM2, 5.5 * KM2)

    @pytest.mark.parametrize("joined_km2", [99.0, 101.0])
    def test_exactly_one_percent_off_passes(self, joined_km2: float) -> None:
        """Refused only beyond the tolerance. Both gaps and 1 % of 100 km2 are
        exactly 1e6 m2 in binary floating point, so the boundary is met exactly."""
        assert nw().check_area(joined_km2 * KM2, 100.0 * KM2) is None

    def test_the_tolerance_is_not_taken_of_the_joined_area(self) -> None:
        """1.005 km2 over 100: beyond 1 % of the stated area, within 1 % of the
        joined 101.005."""
        refusal(nw().check_area, 101.005 * KM2, 100.0 * KM2)


# ---------------------------------------------------------------------------
# The home unit (test 8)
# ---------------------------------------------------------------------------


class TestHomeUnit:
    def test_the_one_covering_unit(self) -> None:
        point = Point(X0 + 500, Y0 + 500)
        assert nw().home_unit(point, [("SE100", square(0, 0))]) == "SE100"

    def test_a_point_in_no_unit_is_refused(self) -> None:
        refusal(nw().home_unit, Point(X0 - 500, Y0 + 500), [])

    def test_a_point_on_the_shared_edge_of_two_is_refused_naming_both(self) -> None:
        point = Point(X0 + KM, Y0 + 500)
        said = refusal(nw().home_unit, point, [("SE100", square(0, 0)), ("SE300", square(1, 0))])
        assert "SE100" in said and "SE300" in said, said


# ---------------------------------------------------------------------------
# The outlet distance (test 9)
# ---------------------------------------------------------------------------


OUTLET = Point(X0 + 500, Y0)


class TestOutletDistance:
    @pytest.mark.parametrize("metres", [199.9, 200.0, 200.1])
    def test_a_point_north_of_the_outlet(self, metres: float) -> None:
        seed = Point(X0 + 500, Y0 + metres)
        got = nw().outlet_distance_m(seed, OUTLET, CRS)
        assert got == pytest.approx(metres, abs=1e-6)

    def test_a_lake_is_as_near_as_its_nearest_edge(self) -> None:
        lake = box(X0 + 400, Y0 + 85, X0 + 600, Y0 + 300)
        assert nw().outlet_distance_m(lake, OUTLET, CRS) == pytest.approx(85.0, abs=1e-6)

    def test_an_outlet_inside_the_lake_is_at_distance_zero(self) -> None:
        lake = box(X0 + 400, Y0 - 50, X0 + 600, Y0 + 300)
        assert nw().outlet_distance_m(lake, OUTLET, CRS) == 0.0

    @pytest.mark.parametrize("metres", [85.0, 200.0, 13_988.0])
    def test_in_degrees_the_distance_is_in_metres(self, metres: float) -> None:
        seed = to_4326(Point(X0 + 500, Y0 + metres))
        got = nw().outlet_distance_m(seed, to_4326(OUTLET), "EPSG:4326")
        assert got == pytest.approx(metres, rel=1e-3)

    def test_a_lake_in_degrees(self) -> None:
        lake = to_4326(box(X0 + 400, Y0 + 85, X0 + 600, Y0 + 300))
        got = nw().outlet_distance_m(lake, to_4326(OUTLET), "EPSG:4326")
        assert got == pytest.approx(85.0, rel=1e-3)


# ---------------------------------------------------------------------------
# The cross-check's figures (test 10)
# ---------------------------------------------------------------------------


class TestCompare:
    def test_one_unit_shared_one_outside_and_a_strip_of_terrain_only(self) -> None:
        units = {"SE100": square(0, 0), "SE300": square(1, 0)}
        network_polygon = shapely.union_all(list(units.values()))
        strip = box(X0 - 100, Y0, X0, Y0 + KM)  # 0.1 km2 west of SE100
        terrain = shapely.union_all([units["SE100"], strip])
        check = nw().compare(network_polygon, terrain, units, CRS)
        assert check.terrain_area_km2 == pytest.approx(1.1, rel=1e-9)
        assert check.shared_km2 == pytest.approx(1.0, rel=1e-9)
        assert check.network_only_km2 == pytest.approx(1.0, rel=1e-9)
        assert check.terrain_only_km2 == pytest.approx(0.1, rel=1e-9)
        assert len(check.units_mostly_outside) == 1
        assert "SE300" in str(check.units_mostly_outside)
        assert "SE100" not in str(check.units_mostly_outside)

    def test_a_unit_exactly_half_inside_is_not_mostly_outside(self) -> None:
        """Under one half is mostly outside; one half exactly is not."""
        units = {"SE100": square(0, 0), "SE300": square(1, 0)}
        network_polygon = shapely.union_all(list(units.values()))
        terrain = box(X0, Y0, X0 + 1.5 * KM, Y0 + KM)
        check = nw().compare(network_polygon, terrain, units, CRS)
        assert len(check.units_mostly_outside) == 0
        assert check.network_only_km2 == pytest.approx(0.5, rel=1e-9)
