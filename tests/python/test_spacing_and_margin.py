"""Sub-metre resample spacing and the feature read margin in cells: increment 43,
section 6, tests 13 and 14 (`docs/increments/43-units-and-projections.md`,
sections 4.6 and 4.7).

HOW THIS FILE GOES RED. At the design's base `default_spacing` rounds to whole
metres, at least 1, so 0.5 m comes back as 1; `TargetGrid.spacing` is a
strict integer, so 0.5 and 2^-10 are refused; `feature_input.source_region`
takes no `cell` and always grows the domain by 100 m. After the green step
(section 9.5) the floor is 2^-3 m: at 0f67a567 2^-4 and 2^-10 are accepted,
and a 1 cm source gets 2^-6 m.
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError
from shapely.geometry import Point, Polygon

from feature_fixtures import UTM33, at, domain_of, square
from tin_engine import feature_input, target_grid
from tin_engine.io.models import RasterMeta

ARC_SECOND = 1.0 / 3600.0


def geographic(step: float) -> RasterMeta:
    """A WGS 84 grid of `step` degrees whose upper-left node is (10 E, 50.01 N)."""
    return RasterMeta(
        x_min=10.0,
        y_max=50.01,
        delta_x=step,
        delta_y=step,
        cols=100,
        rows=100,
        epsg=4326,
        geographic=True,
        nodata=None,
        nodata_source="absent",
        pixel_is_area=False,
        vertical_unit_assumed=True,
    )


class TestDefaultSpacing:
    """Test 13: below a metre the spacing is the nearest power of two."""

    @pytest.mark.parametrize(
        ("step", "expected"),
        [(0.0000045, 0.5), (ARC_SECOND / 9, 3), (ARC_SECOND, 31)],
        ids=["0.5_m", "one_ninth_arc_second", "one_arc_second"],
    )
    def test_a_geographic_source_at_50_north(self, step: float, expected: float) -> None:
        """0.500531 m, 3.432996 m and 30.896962 m (section 4.6's table)."""
        assert target_grid.default_spacing(geographic(step), (10.0, 50.0)) == expected

    def test_a_centimetre_source_stops_at_the_floor(self) -> None:
        """1e-7 degrees is about 1.1 cm; the floor is 2^-3 m, the finest spacing
        whose nodes the 1 mm snap leaves fixed (section 9.5)."""
        assert target_grid.default_spacing(geographic(1e-7), (10.0, 50.0)) == 0.125

    def test_a_projected_10_m_source_is_its_own(self) -> None:
        meta = geographic(10.0).model_copy(update={"geographic": False, "epsg": 25833})
        assert target_grid.default_spacing(meta, (500_000.0, 6_600_000.0)) == 10


class TestTargetGridSpacing:
    """Test 13: a whole number from 1 to 1000, or 2^-k for k from 1 to 3."""

    @staticmethod
    def grid(spacing: float) -> target_grid.TargetGrid:
        return target_grid.TargetGrid(
            crs="EPSG:25833", spacing=spacing, row0=0, col0=0, rows=2, cols=2
        )

    @pytest.mark.parametrize("spacing", [0.5, 2.0**-3], ids=["half", "2^-3"])
    def test_a_power_of_two_below_a_metre_is_accepted(self, spacing: float) -> None:
        assert self.grid(spacing).spacing == spacing

    @pytest.mark.parametrize(
        "spacing", [0.75, 0, 2.0**-4, 2.0**-10], ids=["0.75", "0", "2^-4", "2^-10"]
    )
    def test_anything_else_below_a_metre_is_refused(self, spacing: float) -> None:
        """2^-4 m and finer move nodes off the noder's 1 mm lattice (section 9.5)."""
        with pytest.raises(ValidationError, match=r"2\^-k m for k 1 to 3; got"):
            self.grid(spacing)


class TestFeatureMargin:
    """Test 14: the region is the domain plus max(100 m, 2 cells)."""

    DOMAIN = domain_of(square(0, 0, 300, 300))
    #: 1.5 km east of the domain's east side, level with its middle.
    OUTSIDE = Point(at(1_800, 150))

    def region(self, cell: float) -> Polygon:
        return feature_input.source_region(self.DOMAIN, UTM33, UTM33, cell=cell)

    def test_a_1_km_cell_reaches_2_km(self) -> None:
        region = self.region(1_000.0)
        assert region.contains(self.OUTSIDE)
        assert not region.contains(Point(at(2_400, 150)))

    def test_a_10_m_cell_keeps_100_m(self) -> None:
        region = self.region(10.0)
        assert not region.contains(self.OUTSIDE)
        assert region.contains(Point(at(390, 150)))
