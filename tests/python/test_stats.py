# The report is Markdown with a multiplication sign, ≥ and °, and the golden
# string's lines are the report's, so the ambiguous-character and line-length
# rules do not apply.
# ruff: noqa: RUF001, RUF002, E501
"""`tin_engine.stats`: increment 17, Q1 to Q4, M1 and P1, without the extension.

`docs/increments/17-mesh-stats.md` R1, R3 to R5. Not invariant-critical, so no
mutation round: a wrong statistic is visible on the page.

The module is fetched inside a fixture, so its absence fails these tests and
leaves the rest of the session collecting.

Names the design leaves open, chosen here (the tests are their statement):

- ``PhaseClock(now=...)``: ``now`` is a zero-argument callable returning integer
  nanoseconds, ``time.perf_counter_ns`` by default. ``clock.phases()`` returns
  ``(name, seconds)`` pairs in first-seen order.
- ``Quality`` fields: ``angle_median``, ``angle_under_1``, ``angle_under_10``,
  ``angle_worst`` (degrees; the two shares are fractions in [0, 1]) and
  ``degree_median``, ``degree_p99``, ``degree_max``, ``degree_at_least_12``,
  ``degree_at_least_20`` (ints). The degree median uses the same
  ``inverted_cdf`` percentile as p99, so it is an observed degree.
- ``Refinement(tolerance, max_error, rounds, inserted, flips, uncovered,
  carved=None)``. Ola chose C2 (a); the column sits after ``inserted``,
  because carving splits are a subset of the inserts, and is omitted when
  ``carved`` is None.
- ``Sizes``: ``output_vertices``, ``output_triangles``, ``constraint_edges``,
  ``files`` (``(name, bytes)`` pairs) required; ``dem_nodes`` (rows, cols),
  ``dem_spacing`` (dx, dy), ``domain_vertices``, ``domain_holes``,
  ``start_vertices``, ``start_triangles`` and ``dropped`` default to None, and
  a None row is omitted (R4). Spacing prints once when dx == dy, else
  ``dx × dy``. File sizes are decimal: B, kB, MB, GB, one decimal above bytes.
- ``Report(command, sizes, quality, refinement, phases, total, stats_seconds,
  threads=None)``. ``threads`` is the design's "Threads: 10 (hardware
  concurrency)" sentence, which R1's list of Report fields does not carry;
  printed only when not None. A phase named ``"<parent>: ..."`` where
  ``<parent>`` is an earlier phase is a sub-row: listed, not added to the
  total, and announced by the line under the table.
"""

from __future__ import annotations

import importlib
import math
from types import ModuleType
from typing import Any

import numpy as np
import pytest


@pytest.fixture
def stats() -> ModuleType:
    return importlib.import_module("tin_engine.stats")


def tri_mesh(*triangles: tuple[tuple[float, float], ...]) -> tuple[np.ndarray, np.ndarray]:
    """Separate triangles, each with its own three vertices."""
    xy = np.array([p for t in triangles for p in t], dtype=np.float64)
    tris = np.arange(len(xy), dtype=np.int64).reshape(-1, 3)
    return xy, tris


def fan(centre: tuple[float, float], n: int, offset: int) -> tuple[np.ndarray, np.ndarray]:
    """``n`` triangles around one centre: the centre has degree n, each rim vertex 2."""
    angles = np.linspace(0.0, 2 * math.pi, n, endpoint=False)
    rim = np.column_stack([centre[0] + np.cos(angles), centre[1] + np.sin(angles)])
    xy = np.vstack([np.array([centre]), rim])
    k = np.arange(n)
    tris = np.column_stack([np.zeros(n, dtype=np.int64), 1 + k, 1 + (k + 1) % n]) + offset
    return xy, tris


class TestQuality:
    def test_q1_known_angles(self, stats: ModuleType) -> None:
        xy, tris = tri_mesh(
            ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),  # right isosceles: 45
            ((0.0, 0.0), (math.sqrt(3.0), 0.0), (0.0, 1.0)),  # 30-60-90: 30
        )
        q = stats.quality(xy, tris)
        assert q.angle_median == pytest.approx(37.5, rel=1e-12)
        assert q.angle_worst == pytest.approx(30.0, rel=1e-12)
        assert q.angle_under_1 == 0.0
        assert q.angle_under_10 == 0.0

    def test_q2_a_sliver_keeps_its_precision(self, stats: ModuleType) -> None:
        # acos of the normalised dot product is off by 1.1e-5 relative here
        # (cos is 1 - 2e-12, one ulp from 1 is 1.1e-16); atan2 of cross and
        # dot is exact to the ulp. Measured against both spellings with numpy.
        xy, tris = tri_mesh(((0.0, 0.0), (1.0, 0.0), (0.5, 1e-6)))
        q = stats.quality(xy, tris)
        expected = math.degrees(math.atan(2e-6))
        assert q.angle_worst == pytest.approx(expected, rel=1e-9)
        assert q.angle_median == pytest.approx(expected, rel=1e-9)
        assert q.angle_under_1 == 1.0
        assert q.angle_under_10 == 1.0

    @pytest.mark.parametrize(
        "triangle",
        [
            ((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)),  # collinear
            ((0.0, 0.0), (0.0, 0.0), (1.0, 0.0)),  # a duplicate vertex
            ((3.0, 4.0), (3.0, 4.0), (3.0, 4.0)),  # a point
        ],
        ids=["collinear", "duplicate", "point"],
    )
    def test_a_zero_area_triangle_reports_zero_degrees(
        self, stats: ModuleType, triangle: tuple[tuple[float, float], ...]
    ) -> None:
        xy, tris = tri_mesh(triangle)
        q = stats.quality(xy, tris)
        assert q.angle_worst == 0.0
        assert not math.isnan(q.angle_median)

    def test_the_shares_are_by_triangle(self, stats: ModuleType) -> None:
        xy, tris = tri_mesh(
            ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),  # 45
            ((0.0, 0.0), (1.0, 0.0), (0.5, 0.05)),  # ~5.7: under 10, not under 1
            ((0.0, 0.0), (1.0, 0.0), (0.5, 0.001)),  # ~0.11: under both
            ((0.0, 0.0), (1.0, 0.0), (0.5, 0.5)),  # 45
        )
        q = stats.quality(xy, tris)
        assert q.angle_under_1 == 0.25
        assert q.angle_under_10 == 0.5

    def test_q3_known_degrees(self, stats: ModuleType) -> None:
        xy12, t12 = fan((0.0, 0.0), 12, 0)
        xy20, t20 = fan((10.0, 0.0), 20, len(xy12))
        q = stats.quality(np.vstack([xy12, xy20]), np.vstack([t12, t20]))
        # 34 vertices: 32 rim vertices of degree 2, one 12 and one 20.
        assert q.degree_max == 20
        assert q.degree_at_least_12 == 2
        assert q.degree_at_least_20 == 1
        assert q.degree_median == 2
        assert q.degree_p99 == 20  # inverted_cdf: an observed degree, not 19.x
        for value in (q.degree_median, q.degree_p99, q.degree_max):
            assert isinstance(value, int)

    def test_q4_plan_view_reads_x_and_y_only(self, stats: ModuleType) -> None:
        xy, tris = tri_mesh(
            ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
            ((0.0, 0.0), (math.sqrt(3.0), 0.0), (0.0, 1.0)),
            ((0.0, 0.0), (1.0, 0.0), (0.5, 1e-6)),
        )
        z = np.random.default_rng(17).uniform(-1e4, 1e4, len(xy))
        assert stats.quality(np.column_stack([xy, z]), tris) == stats.quality(xy, tris)


class FakeNow:
    """A monotonic nanosecond source that steps by a fixed amount per call."""

    def __init__(self, step_ns: int) -> None:
        self.t = 0
        self.step = step_ns

    def __call__(self) -> int:
        self.t += self.step
        return self.t


class TestPhaseClock:
    """P1. No test sleeps; every duration comes from the fake source."""

    def test_a_phase_is_the_time_between_enter_and_exit(self, stats: ModuleType) -> None:
        clock = stats.PhaseClock(now=FakeNow(250_000_000))
        with clock.phase("decode"):
            pass
        assert clock.phases() == [("decode", pytest.approx(0.25, abs=1e-12))]

    def test_first_seen_order_and_a_repeated_name_accumulates(self, stats: ModuleType) -> None:
        clock = stats.PhaseClock(now=FakeNow(1_000_000))  # 1 ms per call
        with clock.phase("write: encode"):
            pass
        with clock.phase("write: disk"):
            pass
        with clock.phase("write: encode"):
            pass
        clock.add("refine: scan (parallel)", 0.5)
        clock.add("write: disk", 0.25)
        names = [name for name, _ in clock.phases()]
        assert names == ["write: encode", "write: disk", "refine: scan (parallel)"]
        seconds = dict(clock.phases())
        assert seconds["write: encode"] == pytest.approx(0.002, abs=1e-12)
        assert seconds["write: disk"] == pytest.approx(0.251, abs=1e-12)
        assert seconds["refine: scan (parallel)"] == 0.5

    def test_a_phase_that_raises_is_still_recorded(self, stats: ModuleType) -> None:
        clock = stats.PhaseClock(now=FakeNow(1_000))
        with pytest.raises(ValueError, match="boom"), clock.phase("decode"):
            raise ValueError("boom")
        assert [name for name, _ in clock.phases()] == ["decode"]

    def test_the_default_source_is_monotonic_and_non_negative(self, stats: ModuleType) -> None:
        clock = stats.PhaseClock()
        with clock.phase("x"):
            pass
        ((name, seconds),) = clock.phases()
        assert name == "x" and seconds >= 0.0


# M1. The design's R4 sample, rendered from a fixed Report. Its timings add up:
# the top-level rows sum to 3.999 of a 4.030 total, so "other" is 0.031, and
# the four refine sub-rows sum to refine's 0.412.
GOLDEN = """\
# rasputin mesh — statistics

`rasputin mesh --dem 7908_3_10m_z33.tif --domain quarter.geojson --tolerance 1 --out quarter.vtk --stats quarter.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 5051 × 5051 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 536 |
| start triangles | 534 |
| output vertices | 214210 |
| output triangles | 427779 |
| constraint edges | 536 |
| vertices without data dropped | 0 |
| quarter.vtk | 31.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.03 % | 2.41 % | 0.0117° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 43 | 118 | 9 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered |
|---|---|---|---|---|---|---|
| 1 m | 0.999998 m | 40 | 213674 | 0 | 452067 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| decode | 0.061 | 1.5 % |
| domain read | 0.004 | 0.1 % |
| start mesh: build | 0.001 | 0.0 % |
| start mesh: node | 0.001 | 0.0 % |
| start mesh: triangulate | 0.001 | 0.0 % |
| start mesh: constraint edges | 0.001 | 0.0 % |
| refine | 0.412 | 10.2 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: scan (parallel) | 0.118 | 2.9 % |
| refine: split + flip (serial) | 0.161 | 4.0 % |
| refine: setup + output | 0.133 | 3.3 % |
| trim | 0.028 | 0.7 % |
| write: encode | 3.471 | 86.1 % |
| write: disk | 0.019 | 0.5 % |
| other | 0.031 | 0.8 % |
| **total** | **4.030** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.072 s, not included above.
"""

PHASES = (
    ("decode", 0.061),
    ("domain read", 0.004),
    ("start mesh: build", 0.001),
    ("start mesh: node", 0.001),
    ("start mesh: triangulate", 0.001),
    ("start mesh: constraint edges", 0.001),
    ("refine", 0.412),
    ("refine: legalise start", 0.0),
    ("refine: scan (parallel)", 0.118),
    ("refine: split + flip (serial)", 0.161),
    ("refine: setup + output", 0.133),
    ("trim", 0.028),
    ("write: encode", 3.471),
    ("write: disk", 0.019),
)


@pytest.fixture
def quality(stats: ModuleType) -> Any:
    return stats.Quality(
        angle_median=45.0,
        angle_under_1=0.0003,
        angle_under_10=0.0241,
        angle_worst=0.0117,
        degree_median=6,
        degree_p99=9,
        degree_max=43,
        degree_at_least_12=118,
        degree_at_least_20=9,
    )


@pytest.fixture
def quarter(stats: ModuleType, quality: Any) -> Any:
    return stats.Report(
        command=(
            "rasputin mesh --dem 7908_3_10m_z33.tif --domain quarter.geojson "
            "--tolerance 1 --out quarter.vtk --stats quarter.md"
        ),
        sizes=stats.Sizes(
            dem_nodes=(5051, 5051),
            dem_spacing=(10.0, 10.0),
            domain_vertices=536,
            domain_holes=0,
            start_vertices=536,
            start_triangles=534,
            output_vertices=214210,
            output_triangles=427779,
            constraint_edges=536,
            dropped=0,
            files=(("quarter.vtk", 31_800_000),),
        ),
        quality=quality,
        refinement=stats.Refinement(
            tolerance=1.0,
            max_error=0.999998,
            rounds=40,
            inserted=213674,
            flips=452067,
            uncovered=0,
            carved=0,
        ),
        phases=PHASES,
        total=4.030,
        stats_seconds=0.072,
        threads=10,
    )


@pytest.fixture
def fixture_report(stats: ModuleType, quality: Any) -> Any:
    """A gallery fixture: no DEM, no domain, no refinement, nothing dropped."""
    return stats.Report(
        command="rasputin mesh catchment --flat --out c.ply --stats c.md",
        sizes=stats.Sizes(
            output_vertices=40,
            output_triangles=60,
            constraint_edges=30,
            files=(("c.ply", 2048), ("c-edges.ply", 999)),
        ),
        quality=quality,
        refinement=None,
        phases=(("start mesh: build", 0.002), ("write: encode", 0.001)),
        total=0.004,
        stats_seconds=0.001,
    )


class TestRender:
    """M1."""

    def test_the_golden_report(self, stats: ModuleType, quarter: Any) -> None:
        assert stats.render(quarter) == GOLDEN

    def test_rows_that_do_not_apply_are_omitted(
        self, stats: ModuleType, fixture_report: Any
    ) -> None:
        text = stats.render(fixture_report)
        for absent in (
            "## Refinement",
            "DEM nodes",
            "domain vertices",
            "start vertices",
            "without data dropped",
            "Threads:",
            "Sub-rows",
        ):
            assert absent not in text, absent
        sections = [line for line in text.splitlines() if line.startswith("## ")]
        assert sections == ["## Sizes", "## Quality (plan view, x/y)", "## Timings"]
        assert "| c.ply | 2.0 kB |" in text
        assert "| c-edges.ply | 999 B |" in text

    def test_other_is_total_minus_the_top_level_rows(
        self, stats: ModuleType, fixture_report: Any
    ) -> None:
        # "start mesh: build" has no "start mesh" row, so it is top level.
        text = stats.render(fixture_report)
        assert "| other | 0.001 | 25.0 % |" in text
        assert "| **total** | **0.004** | **100 %** |" in text

    def test_carved_is_omitted_when_not_counted(self, stats: ModuleType, quarter: Any) -> None:
        import dataclasses

        report = dataclasses.replace(
            quarter, refinement=dataclasses.replace(quarter.refinement, carved=None)
        )
        text = stats.render(report)
        assert "| tolerance | achieved max error | rounds | inserted | flips | uncovered |" in text
        assert "carved" not in text

    def test_an_achieved_error_never_prints_above_its_tolerance(
        self, stats: ModuleType, quarter: Any
    ) -> None:
        import dataclasses

        just_under = math.nextafter(1.0, 0.0)
        report = dataclasses.replace(
            quarter, refinement=dataclasses.replace(quarter.refinement, max_error=just_under)
        )
        assert f"| 1 m | {just_under!r} m |" in stats.render(report)

    @pytest.mark.parametrize(
        ("spacing", "holes", "expected"),
        [
            ((10.0, 5.0), 1, ("| DEM nodes | 5051 × 5051 (10 × 5 m) |", "(1 ring, 1 hole)")),
            ((0.5, 0.5), 2, ("| DEM nodes | 5051 × 5051 (0.5 m) |", "(1 ring, 2 holes)")),
        ],
    )
    def test_spacing_and_hole_wording(
        self,
        stats: ModuleType,
        quarter: Any,
        spacing: tuple[float, float],
        holes: int,
        expected: tuple[str, str],
    ) -> None:
        import dataclasses

        sizes = dataclasses.replace(quarter.sizes, dem_spacing=spacing, domain_holes=holes)
        text = stats.render(dataclasses.replace(quarter, sizes=sizes))
        for line in expected:
            assert line in text

    def test_render_is_pure(self, stats: ModuleType, quarter: Any) -> None:
        assert stats.render(quarter) == stats.render(quarter)

