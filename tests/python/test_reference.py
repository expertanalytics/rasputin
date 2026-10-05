"""`tin_engine.reference`: agreement with NVE's polygon, classes, summary (increment 29, PR 4).

`docs/increments/29-nve-reference-catchments.md`, "Agreement and classes"
and "The red suites", PR 4's `test_reference.py`. Pure: hand-built polygons
on a 10 m node lattice, agreement numbers and rows given by hand; no DEM, no
file.

Interface pinned here (the design names `agreement`, `classify`, `summarise`,
`Summary`, the class names, `match_by` and the cause names; what it leaves
open is chosen here and listed in the handback):

- `agreement(fine, reference, meta) -> Agreement`: our fine outline, NVE's
  polygon (a Polygon or a MultiPolygon, all parts counted), and the final
  window's `RasterMeta`, of which only the origin (`x_min`, `y_max`) and the
  spacing are read: the lattice runs on past the window as far as either
  polygon reaches.
- `Agreement`, a frozen dataclass built by keyword: `ours`, `ref`, `both`
  (node counts), `nve_in_ours` (both / ref), `ours_in_nve` (both / ours),
  `area_ratio` (our area / NVE's), `divide_offset_m` and `cell_m` (the
  lattice spacing, which the offset test's 3 cells are counted in).
- `classify(agreement, gauge, refused=False) -> (class, match_by)`:
  `agreement` None without a reference, `gauge` a `catchment.GaugeResult`
  (its `causes`, the joined list PR 4 moves into `catchment.py`, decide
  `uncertain`) or None. The class is `"refused"`, `"uncertain"`, `"match"`,
  `"close"`, `"miss"`, or None for a well-posed station with no reference;
  `match_by` is `"overlap"` or `"offset"` for a match, else None.
- `summarise(rows) -> Summary`, a pydantic model whose
  `model_dump(mode="json")` has the keys used in `TestSummary` below. A row
  is read by attribute (`station`, `station_class`, `match_by`,
  `refusal_cause`, `causes`, `area_ratio`, `nve_in_ours`, `ours_in_nve`,
  `reference_area_km2`, `fine_area_km2`, `tiles`), the same names
  `catchment_batch.StationResult` carries, so a stand-in is passed here.
- Percentiles are numpy's default (linear interpolation between closest
  ranks). A station's size band is decided by NVE's polygon area, or by our
  fine area where there is no reference; a band's lower bound is inside it.

Numeric bounds: lattice counts are exact integers; ratios are compared at
1e-12 relative; offsets, whose scale is at most 500 m here, at 1e-9 m.
"""

from __future__ import annotations

import importlib
import json
from dataclasses import dataclass
from types import ModuleType
from typing import Any

import pytest
from shapely.geometry import MultiPolygon, Polygon, box

from mosaic_fixtures import meta

X0 = 500_000.0
Y0 = 6_600_000.0
D = 10.0
#: A corner far from the window `META` describes: the lattice runs on.
C0, R0 = 1000, 2000
META = meta(rows=2, cols=2, x_min=X0, y_max=Y0, dx=D, dy=D)

KNOWN_LINE = (
    "2 stations refused because their windows select tiles on two different grids, "
    "which rasputin does not combine, and neither grid covers the window alone "
    "(known refusals, not failures)"
)


@pytest.fixture(scope="module")
def ref() -> ModuleType:
    return importlib.import_module("tin_engine.reference")


@pytest.fixture(scope="module")
def catchment() -> ModuleType:
    return importlib.import_module("tin_engine.catchment")


def cells(col: int, row: int, n: int) -> Polygon:
    """The square of n x n lattice nodes from (row, col), edges halfway
    between nodes: n * n nodes strictly inside, area (n D)^2."""
    return box(
        X0 + D * col - D / 2,
        Y0 - D * (row + n - 1) - D / 2,
        X0 + D * (col + n - 1) + D / 2,
        Y0 - D * row + D / 2,
    )


# ---------------------------------------------------------------------------
# Agreement, on the lattice
# ---------------------------------------------------------------------------


class TestAgreement:
    def test_identical_polygons_agree_fully_with_no_offset(self, ref: ModuleType) -> None:
        square = cells(C0, R0, 100)
        a = ref.agreement(square, square, META)
        assert (a.ours, a.ref, a.both) == (10_000, 10_000, 10_000)
        assert a.nve_in_ours == 1.0 and a.ours_in_nve == 1.0
        assert a.area_ratio == pytest.approx(1.0, rel=1e-12)
        assert a.divide_offset_m == pytest.approx(0.0, abs=1e-9)
        assert a.cell_m == D

    def test_a_square_shifted_one_cell_along_x_scores_half_a_cell(self, ref: ModuleType) -> None:
        """Two of the four sides do not move: (100 + 100) cells of area over
        a 4000 m perimeter is 5 m, half a cell."""
        a = ref.agreement(cells(C0 + 1, R0, 100), cells(C0, R0, 100), META)
        assert (a.ours, a.ref, a.both) == (10_000, 10_000, 9_900)
        assert a.nve_in_ours == pytest.approx(0.99, rel=1e-12)
        assert a.ours_in_nve == pytest.approx(0.99, rel=1e-12)
        assert a.area_ratio == pytest.approx(1.0, rel=1e-12)
        assert a.divide_offset_m == pytest.approx(0.5 * D, abs=1e-9)

    def test_a_square_grown_one_cell_every_side_scores_1_01_cells(self, ref: ModuleType) -> None:
        """404 cells of area between the outlines over 400 cells of NVE's
        perimeter: one cell, plus the corner term."""
        a = ref.agreement(cells(C0 - 1, R0 - 1, 102), cells(C0, R0, 100), META)
        assert (a.ours, a.ref, a.both) == (10_404, 10_000, 10_000)
        assert a.nve_in_ours == 1.0
        assert a.ours_in_nve == pytest.approx(10_000 / 10_404, rel=1e-12)
        assert a.area_ratio == pytest.approx(1.0404, rel=1e-12)
        assert a.divide_offset_m == pytest.approx(1.01 * D, abs=1e-9)

    def test_disjoint_polygons_share_nothing(self, ref: ModuleType) -> None:
        a = ref.agreement(cells(C0 + 200, R0, 100), cells(C0, R0, 100), META)
        assert a.both == 0
        assert a.nve_in_ours == 0.0 and a.ours_in_nve == 0.0
        assert a.divide_offset_m == pytest.approx(20_000 * D * D / 4000.0, abs=1e-9)

    def test_every_part_of_a_multipart_reference_counts(self, ref: ModuleType) -> None:
        """NVE's `19.79.0` Gravå has two parts: both are NVE's catchment.
        The perimeter is both parts' (800 m)."""
        first, second = cells(C0, R0, 10), cells(C0 + 50, R0, 10)
        a = ref.agreement(first, MultiPolygon([first, second]), META)
        assert (a.ours, a.ref, a.both) == (100, 200, 100)
        assert a.nve_in_ours == pytest.approx(0.5, rel=1e-12)
        assert a.ours_in_nve == 1.0
        assert a.area_ratio == pytest.approx(0.5, rel=1e-12)
        assert a.divide_offset_m == pytest.approx(100 * D * D / 800.0, abs=1e-9)

    def test_the_window_does_not_bound_the_count(self, ref: ModuleType) -> None:
        """`META` is a 2 x 2 window 10 km from the squares: only its origin
        and spacing place the lattice."""
        small = meta(rows=2, cols=2, x_min=X0 + D * 7, y_max=Y0 - D * 3, dx=D, dy=D)
        square = cells(C0, R0, 30)
        assert ref.agreement(square, square, small).ours == 900

    def test_the_lattice_follows_the_window_origin(self, ref: ModuleType) -> None:
        """Moved half a cell, the lattice puts the square's edges on nodes,
        and a node on an edge is not strictly inside: 29 x 29, not 30 x 30."""
        half = meta(rows=2, cols=2, x_min=X0 + D / 2, y_max=Y0 + D / 2, dx=D, dy=D)
        square = cells(C0, R0, 30)
        assert ref.agreement(square, square, half).ours == 29 * 29


# ---------------------------------------------------------------------------
# Classes, from the numbers given
# ---------------------------------------------------------------------------


def numbers(
    ref: ModuleType, nve_in_ours: float, ours_in_nve: float, offset_m: float = 100.0
) -> Any:
    """An `Agreement` with the two overlaps and the offset given, on DTM10."""
    return ref.Agreement(
        ours=1000,
        ref=1000,
        both=round(1000 * min(nve_in_ours, ours_in_nve)),
        nve_in_ours=nve_in_ours,
        ours_in_nve=ours_in_nve,
        area_ratio=1.0,
        divide_offset_m=offset_m,
        cell_m=D,
    )


def gauge_result(catchment: ModuleType, causes: tuple[str, ...] = (), swing: float = 0.01) -> Any:
    """A `GaugeResult` whose joined `causes` are `causes`; its sensitivity's
    own causes are the same less `direction`, as `delineate` joins them."""
    own = tuple(c for c in causes if c != "direction")
    s = importlib.import_module("tin_engine.sensitivity").Sensitivity(
        a0=10.0, area_up=9.9, area_down=10.0 * (1 + swing), swing=swing,
        largest_step=0.01, largest_step_at_m=10.0, checked_up_m=30.0, checked_down_m=30.0,
        drains="chain_not_draining" not in own, monotone=not own, causes=own,
        well_posed=not own,
    )  # fmt: skip
    return catchment.GaugeResult(
        node=(X0, Y0), chain=((X0, Y0 + D), (X0, Y0)), node_offset_m=0.0, chain_nodes=2,
        lowered_nodes=0, lowered_max_m=0.0, direction_ok="direction" not in causes,
        end_extended_m=0.0, end_closed=True, downstream_checked="whole", sensitivity=s,
        causes=causes,
    )  # fmt: skip


class TestClassify:
    @pytest.mark.parametrize(
        ("overlaps", "expected"),
        [
            ((0.95, 0.95), ("match", "overlap")),
            ((0.95, 0.9499999), ("close", None)),
            ((0.80, 0.80), ("close", None)),
            ((0.80, 0.7999999), ("miss", None)),
            ((1.0, 0.0), ("miss", None)),
        ],
        ids=["95-both", "just-under-95", "80-both", "just-under-80", "one-way-only"],
    )
    def test_the_overlap_bars_are_inclusive(
        self, ref: ModuleType, catchment: ModuleType, overlaps: Any, expected: Any
    ) -> None:
        a = numbers(ref, *overlaps, offset_m=100.0)
        assert ref.classify(a, gauge_result(catchment)) == expected

    @pytest.mark.parametrize(
        ("offset_m", "expected"),
        [(30.0, ("match", "offset")), (30.000001, ("miss", None))],
        ids=["30-m", "just-over-30-m"],
    )
    def test_the_offset_bar_is_three_cells_inclusive(
        self, ref: ModuleType, catchment: ModuleType, offset_m: float, expected: Any
    ) -> None:
        a = numbers(ref, 0.5, 0.5, offset_m=offset_m)
        assert ref.classify(a, gauge_result(catchment)) == expected

    def test_the_overlap_test_is_tried_first(self, ref: ModuleType, catchment: ModuleType) -> None:
        both = numbers(ref, 0.99, 0.99, offset_m=1.0)
        assert ref.classify(both, gauge_result(catchment)) == ("match", "overlap")

    def test_the_offset_test_rescues_a_close_station(
        self, ref: ModuleType, catchment: ModuleType
    ) -> None:
        a = numbers(ref, 0.90, 0.90, offset_m=20.0)
        assert ref.classify(a, gauge_result(catchment)) == ("match", "offset")

    def test_a_swing_of_exactly_0_05_is_scored_and_just_above_is_uncertain(
        self, ref: ModuleType, catchment: ModuleType
    ) -> None:
        """At the bar the sensitivity raises no cause; just above it raises
        `swing` (`sensitivity.assess`), and the class follows."""
        a = numbers(ref, 0.99, 0.99)
        at_bar = gauge_result(catchment, (), swing=0.05)
        above = gauge_result(catchment, ("swing",), swing=0.05 + 1e-9)
        assert ref.classify(a, at_bar) == ("match", "overlap")
        assert ref.classify(a, above) == ("uncertain", None)

    @pytest.mark.parametrize(
        "causes",
        [
            ("downstream_unread",),
            ("chain_not_draining",),
            ("chain_end_open",),
            ("direction",),
            ("swing", "downstream_unread"),
        ],
    )
    def test_any_cause_is_uncertain_whatever_the_overlaps(
        self, ref: ModuleType, catchment: ModuleType, causes: tuple[str, ...]
    ) -> None:
        perfect = numbers(ref, 1.0, 1.0, offset_m=0.0)
        assert ref.classify(perfect, gauge_result(catchment, causes)) == ("uncertain", None)

    def test_refused_beats_uncertain_beats_the_rest(
        self, ref: ModuleType, catchment: ModuleType
    ) -> None:
        perfect = numbers(ref, 1.0, 1.0, offset_m=0.0)
        shaky = gauge_result(catchment, ("swing",))
        assert ref.classify(perfect, shaky, refused=True) == ("refused", None)
        assert ref.classify(None, None, refused=True) == ("refused", None)
        assert ref.classify(perfect, shaky) == ("uncertain", None)

    def test_without_a_reference_only_refused_and_uncertain_are_classes(
        self, ref: ModuleType, catchment: ModuleType
    ) -> None:
        assert ref.classify(None, gauge_result(catchment)) == (None, None)
        assert ref.classify(None, gauge_result(catchment, ("swing",))) == ("uncertain", None)

    def test_without_a_gauge_the_agreement_alone_decides(self, ref: ModuleType) -> None:
        """The nearest-stream fallback (PR 5) has no sensitivity."""
        assert ref.classify(numbers(ref, 0.99, 0.99), None) == ("match", "overlap")


# ---------------------------------------------------------------------------
# The summary
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Row:
    """The attributes of `StationResult` that `summarise` reads."""

    station: str
    station_class: str | None
    match_by: str | None = None
    refusal_cause: str | None = None
    causes: tuple[str, ...] = ()
    area_ratio: float | None = None
    nve_in_ours: float | None = None
    ours_in_nve: float | None = None
    reference_area_km2: float | None = 50.0
    fine_area_km2: float | None = 50.0
    tiles: int | None = 1


def scored(n: int, value: float, **kw: Any) -> Row:
    return Row(
        f"1.{n}.0", kw.pop("station_class", "match"), match_by=kw.pop("match_by", "overlap"),
        area_ratio=value, nve_in_ours=value, ours_in_nve=value, **kw,
    )  # fmt: skip


PERCENTILE_KEYS = ["min", "p10", "p25", "p50", "p75", "p90", "max"]


class TestSummary:
    def dump(self, ref: ModuleType, rows: list[Row]) -> dict[str, Any]:
        out: dict[str, Any] = ref.summarise(rows).model_dump(mode="json")
        return out

    def test_counts_per_class_and_per_match_by(self, ref: ModuleType) -> None:
        rows = [
            scored(1, 1.0),
            scored(2, 1.0),
            scored(3, 0.9, match_by="offset"),
            scored(4, 0.85, station_class="close", match_by=None),
            scored(5, 0.2, station_class="miss", match_by=None),
            Row("1.6.0", "uncertain", causes=("swing",)),
            Row("1.7.0", "refused", refusal_cause="other"),
        ]
        s = self.dump(ref, rows)
        assert s["stations"] == 7
        assert s["classes"] == {"match": 3, "close": 1, "miss": 1, "uncertain": 1, "refused": 1}
        assert s["match_by"] == {"overlap": 2, "offset": 1}

    def test_the_percentiles_of_the_scored_on_a_known_list(self, ref: ModuleType) -> None:
        """Eleven scored values 0.0 to 1.0, and an uncertain station whose
        numbers (5.0) must not enter."""
        rows = [scored(k, k / 10, station_class="miss", match_by=None) for k in range(11)]
        rows.append(Row("9.9.0", "uncertain", causes=("swing",), area_ratio=5.0,
                        nve_in_ours=5.0, ours_in_nve=5.0))  # fmt: skip
        s = self.dump(ref, rows)
        for key in ("area_ratio", "nve_in_ours", "ours_in_nve"):
            got = s["scored"][key]
            assert list(got) == PERCENTILE_KEYS
            expected = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
            assert [got[k] for k in PERCENTILE_KEYS] == pytest.approx(expected, rel=1e-12)

    def test_size_bands_and_the_share_uncertain_in_each(self, ref: ModuleType) -> None:
        """Bands by NVE's area, lower bound inside: 10 km2 is in 10-100,
        1000 in over 1000. Without a reference, our fine area decides."""
        rows = [
            scored(1, 1.0, reference_area_km2=9.99),
            scored(2, 1.0, reference_area_km2=10.0),
            Row("1.3.0", "uncertain", causes=("swing",), reference_area_km2=50.0),
            Row("1.4.0", "uncertain", causes=("swing",), reference_area_km2=99.0),
            scored(5, 0.2, station_class="miss", match_by=None, reference_area_km2=100.0),
            scored(6, 1.0, reference_area_km2=1000.0),
            scored(7, 1.0, reference_area_km2=None, fine_area_km2=2000.0),
        ]
        bands = self.dump(ref, rows)["by_size"]
        assert list(bands) == ["under 10", "10-100", "100-1000", "over 1000"]
        assert [bands[b]["stations"] for b in bands] == [1, 3, 1, 2]
        assert bands["10-100"]["uncertain_share"] == pytest.approx(2 / 3, rel=1e-12)
        assert bands["under 10"]["uncertain_share"] == 0.0
        assert bands["10-100"]["classes"]["match"] == 1
        assert bands["100-1000"]["area_ratio"]["p50"] == pytest.approx(0.2, rel=1e-12)

    def test_tile_count_groups(self, ref: ModuleType) -> None:
        rows = [scored(k, 1.0, tiles=t) for k, t in enumerate([1, 2, 3, 4, 5, 9])]
        tiles = self.dump(ref, rows)["by_tiles"]
        assert list(tiles) == ["1", "2", "3-4", "5+"]
        assert [tiles[g]["stations"] for g in tiles] == [1, 1, 2, 2]

    def test_each_cause_of_uncertain_is_counted_apart(self, ref: ModuleType) -> None:
        """A station with two causes counts in both."""
        rows = [
            Row("1.1.0", "uncertain", causes=("swing", "downstream_unread")),
            Row("1.2.0", "uncertain", causes=("swing",)),
            Row("1.3.0", "uncertain", causes=("chain_not_draining", "chain_end_open")),
            Row("1.4.0", "uncertain", causes=("direction",)),
        ]
        causes = self.dump(ref, rows)["uncertain_causes"]
        assert causes["swing"] == 2
        assert causes["downstream_unread"] == 1
        assert causes["chain_not_draining"] == 1
        assert causes["chain_end_open"] == 1
        assert causes["direction"] == 1

    def test_mixed_grid_refusals_are_known_refusals_never_scored(self, ref: ModuleType) -> None:
        """Counted by cause, and apart with their line and station numbers in
        row order; a stray number on such a row never reaches the scored."""
        rows = [
            scored(1, 1.0),
            Row("196.11.0", "refused", refusal_cause="mixed_grid", area_ratio=9.0,
                nve_in_ours=9.0, ours_in_nve=9.0),
            Row("2.2.0", "refused", refusal_cause="no_river"),
            Row("156.15.0", "refused", refusal_cause="mixed_grid"),
            Row("3.3.0", "refused", refusal_cause="other"),
        ]  # fmt: skip
        s = self.dump(ref, rows)
        assert s["refusal_causes"] == {"mixed_grid": 2, "no_river": 1, "other": 1}
        known = s["known_refusals"]
        assert known["count"] == 2
        assert known["stations"] == ["196.11.0", "156.15.0"]
        assert known["line"] == KNOWN_LINE
        assert s["classes"]["refused"] == 4
        assert s["scored"]["area_ratio"]["max"] == pytest.approx(1.0, rel=1e-12)

    def test_the_json_is_deterministic(self, ref: ModuleType) -> None:
        rows = [scored(1, 1.0), Row("1.2.0", "uncertain", causes=("swing",))]
        first = ref.summarise(rows).model_dump_json()
        assert first == ref.summarise(list(rows)).model_dump_json()
        assert json.loads(first)["stations"] == 2
