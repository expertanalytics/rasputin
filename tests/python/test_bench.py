"""`tools/bench.py`, the acceptance benchmark (`docs/benchmarks/bench-py.md`).

The design's "Seams for tests/python/test_bench.py" lists what is tested here,
and its "Pinned by the red suite" section fixes every name, field and string
these tests assert on. Nothing here needs the 1 m DEM or a real build: every
subprocess goes through an injected `Runner`, and only the last test starts a
real child, on a synthetic DEM, to prove the `cli.refine` wrap still matches the
real signature.

`tools/` is not a package, so the module is loaded from its path, as
`test_session_state.py` loads its tool. It is loaded in a fixture rather than at
import time, so that while `bench.py` is missing each test fails on its own
instead of one collection error stopping the whole suite.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import math
import re
import subprocess
import sys
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from typer.testing import CliRunner

from geotiff_fixtures import micro_tiff

TOOL = Path(__file__).resolve().parents[2] / "tools" / "bench.py"

cli = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})


@pytest.fixture(scope="module")
def bench() -> Iterator[ModuleType]:
    spec = importlib.util.spec_from_file_location("bench", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Pydantic resolves postponed annotations through sys.modules[__module__].
    sys.modules["bench"] = module
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop("bench", None)


# --------------------------------------------------------------------- pmset

AC_WITH_BATTERY = """\
Now drawing from 'AC Power'
 -InternalBattery-0 (id=7929955)\t100%; charged; 0:00 remaining present: true
"""
AC_DESKTOP = "Now drawing from 'AC Power'\n"
BATTERY = """\
Now drawing from 'Battery Power'
 -InternalBattery-0 (id=7929955)\t88%; discharging; 6:25 remaining present: true
"""
UPS = "Now drawing from 'UPS Power'\n"


@pytest.mark.parametrize(
    ("text", "state", "percent"),
    [
        pytest.param(AC_WITH_BATTERY, "ac", 100, id="ac-laptop"),
        pytest.param(AC_DESKTOP, "ac", None, id="ac-desktop-no-battery-line"),
        pytest.param(BATTERY, "battery", 88, id="battery"),
        pytest.param(UPS, "unknown", None, id="ups"),
        pytest.param("", "unknown", None, id="empty"),
        pytest.param("pmset: command not found\n\x00\x01", "unknown", None, id="garbage"),
    ],
)
def test_parse_pmset(bench: ModuleType, text: str, state: str, percent: int | None) -> None:
    power = bench.parse_pmset(text)
    assert power.state == state
    assert power.percent == percent
    assert power.raw == text


@pytest.mark.parametrize(
    ("before", "after", "state"),
    [
        (AC_WITH_BATTERY, AC_DESKTOP, "ac"),
        (BATTERY, BATTERY, "battery"),
        (AC_WITH_BATTERY, BATTERY, "mixed"),
        (BATTERY, AC_WITH_BATTERY, "mixed"),
        ("", "", "unknown"),
    ],
)
def test_combine_power(bench: ModuleType, before: str, after: str, state: str) -> None:
    combined = bench.combine_power(bench.parse_pmset(before), bench.parse_pmset(after))
    assert combined.state == state


# --------------------------------------------------------------- child output

CHILD_JSON = {
    "refine_s": 0.4812,
    "app_s": 0.9731,
    "max_error": 0.9981,
    "rounds": 7,
    "inserted": 1234,
    "flips": 5678,
}


def bench_line(values: dict[str, Any]) -> str:
    return "BENCH " + json.dumps(values)


def test_parse_child_reads_the_bench_line_among_cli_noise(bench: ModuleType) -> None:
    stderr = "\n".join(
        ["refined from DEM nodes, tolerance 1 m", bench_line(CHILD_JSON), "wrote x.vtk", ""]
    )
    child = bench.parse_child(stderr, run="tile t=4 r=2")
    assert child.refine_s == pytest.approx(0.4812)
    assert child.app_s == pytest.approx(0.9731)
    assert child.max_error == pytest.approx(0.9981)
    assert (child.rounds, child.inserted, child.flips) == (7, 1234, 5678)


@pytest.mark.parametrize(
    "stderr",
    [
        pytest.param("", id="empty"),
        pytest.param("Traceback (most recent call last):\n  boom\n", id="no-bench-line"),
        pytest.param("BENCH {refine_s: 0.5\n", id="malformed-json"),
        pytest.param(
            bench_line({k: v for k, v in CHILD_JSON.items() if k != "refine_s"}),
            id="missing-refine_s",
        ),
        pytest.param(bench_line(CHILD_JSON) + "\n" + bench_line(CHILD_JSON), id="two-lines"),
    ],
)
def test_parse_child_refuses_with_the_run_named(bench: ModuleType, stderr: str) -> None:
    with pytest.raises(bench.ChildError, match=re.escape("quarter t=20 r=3")):
        bench.parse_child(stderr, run="quarter t=20 r=3")


def test_child_error_is_a_value_error(bench: ModuleType) -> None:
    assert issubclass(bench.ChildError, ValueError)


# ------------------------------------------------------- statistics, ceiling


def sample(bench: ModuleType, domain: str, threads: int, repeat: int, refine_s: float) -> Any:
    return bench.Sample(
        domain=domain,
        threads=threads,
        repeat=repeat,
        refine_s=refine_s,
        app_s=refine_s + 0.5,
        proc_s=refine_s + 1.0,
        max_error=0.99,
        rounds=3,
        inserted=10,
        flips=20,
    )


def test_median_stats_groups_by_domain_and_threads(bench: ModuleType) -> None:
    samples = [
        sample(bench, "tile", 4, 0, 3.0),
        sample(bench, "quarter", 1, 0, 0.5),
        sample(bench, "tile", 4, 1, 1.0),
        sample(bench, "tile", 0, 0, 7.0),
        sample(bench, "tile", 4, 2, 2.0),
        sample(bench, "quarter", 1, 1, 0.7),
    ]
    stats = bench.median_stats(samples)
    assert [(s.domain, s.threads) for s in stats] == [("quarter", 1), ("tile", 0), ("tile", 4)]
    quarter, tile0, tile4 = stats
    assert (tile4.n, tile4.median, tile4.min, tile4.max) == (3, 2.0, 1.0, 3.0)
    assert quarter.median == pytest.approx(0.6)  # even count: mean of the middle two
    assert (tile0.n, tile0.median) == (1, 7.0)


def test_ceiling_reports_top_and_best_speedup_ignoring_default(bench: ModuleType) -> None:
    ceiling = bench.ceiling({0: 0.1, 1: 10.0, 2: 6.0, 4: 4.0, 20: 5.0})
    assert ceiling.top_threads == 20
    assert ceiling.at_top == pytest.approx(2.0)
    assert ceiling.best == pytest.approx(2.5)
    assert ceiling.best_threads == 4


# ------------------------------------------------------------------ quality

#: A kite whose long diagonal is the shared edge: the non-Delaunay split. The
#: circumcircle of (-2,0), (2,0), (0,1) is centred at (0,-1.5), radius 2.5, so
#: the apex (0,-1) is strictly inside it.
KITE = np.array([[-2.0, 0.0], [0.0, -1.0], [2.0, 0.0], [0.0, 1.0]])
KITE_LONG_SPLIT = np.array([[0, 2, 3], [0, 1, 2]])
KITE_LONG_EDGE = np.array([[0, 2]])

UTM = (500_000.0, 6_600_000.0)


@pytest.mark.parametrize("offset", [(0.0, 0.0), UTM], ids=["origin", "utm-offset"])
def test_quality_counts_the_non_delaunay_split(
    bench: ModuleType, offset: tuple[float, float]
) -> None:
    q = bench.quality(KITE + offset, KITE_LONG_SPLIT, np.zeros((0, 2), np.int64), 1.0, 0.5)
    assert q.delaunay_checked == 1
    assert q.delaunay_violations == 1


def test_quality_skips_a_constrained_edge(bench: ModuleType) -> None:
    q = bench.quality(KITE, KITE_LONG_SPLIT, KITE_LONG_EDGE, 1.0, 0.5)
    assert q.delaunay_checked == 0
    assert q.delaunay_violations == 0


def test_quality_cocircular_square_is_ambiguous_and_decided_exactly(bench: ModuleType) -> None:
    square = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    q = bench.quality(square, np.array([[0, 1, 2], [0, 2, 3]]), np.zeros((0, 2)), 1.0, 0.5)
    assert q.delaunay_checked == 1
    assert q.delaunay_ambiguous == 1
    assert q.delaunay_violations == 0  # on the circle is not strictly inside


def test_quality_sliver_sets_the_worst_angle(bench: ModuleType) -> None:
    sliver = np.array([[0.0, 0.0], [10.0, 0.0], [5.0, 0.01]])
    q = bench.quality(sliver, np.array([[0, 1, 2]]), np.zeros((0, 2)), 1.0, 0.5)
    assert q.worst_angle == pytest.approx(math.degrees(math.atan2(0.01, 5.0)), rel=1e-12)
    assert q.share_under_1 == 1.0


def test_quality_zero_area_triangle_has_worst_angle_zero(bench: ModuleType) -> None:
    flat = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    q = bench.quality(flat, np.array([[0, 1, 2]]), np.zeros((0, 2)), 1.0, 0.5)
    assert q.worst_angle == 0.0


def fan(spokes: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """A centre vertex and a regular `spokes`-gon: the centre has degree `spokes`."""
    angles = np.linspace(0.0, 2.0 * np.pi, spokes, endpoint=False)
    rim = np.column_stack([np.cos(angles), np.sin(angles)])
    points = np.vstack([[0.0, 0.0], rim])
    triangles = np.array([[0, 1 + k, 1 + (k + 1) % spokes] for k in range(spokes)])
    return points, triangles


def test_quality_fan_sets_max_degree_and_agrees_with_stats(bench: ModuleType) -> None:
    from tin_engine.stats import quality as stats_quality

    points, triangles = fan(8)
    q = bench.quality(points, triangles, np.zeros((0, 2)), 1.0, 0.5)
    reference = stats_quality(np.column_stack([points, np.zeros(len(points))]), triangles)
    assert q.max_degree == 8
    assert q.angle_median == reference.angle_median
    assert q.worst_angle == reference.angle_worst
    assert q.share_under_1 == reference.angle_under_1
    assert q.delaunay_violations == 0


@pytest.mark.parametrize(
    ("max_error", "within"),
    [(0.5, True), (1.0, True), (math.nextafter(1.0, 2.0), False), (2.0, False)],
    ids=["below", "equal", "one-ulp-above", "above"],
)
def test_quality_tolerance_check(bench: ModuleType, max_error: float, within: bool) -> None:
    points, triangles = fan(6)
    q = bench.quality(points, triangles, np.zeros((0, 2)), 1.0, max_error)
    assert q.within_tolerance is within


def write_mesh(path: Path, *, binary: bool = False, fields: Sequence[tuple[str, str]] = ()) -> None:
    from tin_engine.features import DEFAULT_VOCABULARY
    from tin_engine.io.vtk_legacy import write_vtk

    points, triangles = fan(6)
    path.write_bytes(
        write_vtk(
            np.column_stack([points, np.arange(len(points), dtype=np.float64)]),
            triangles=triangles,
            edges=np.array([[1, 2], [2, 3]]),
            edge_masks=np.array([1, 2]),
            vocabulary=DEFAULT_VOCABULARY,
            fields=fields,
            binary=binary,
        )
    )


def test_read_vtk_ascii_reads_what_write_vtk_wrote(bench: ModuleType, tmp_path: Path) -> None:
    write_mesh(tmp_path / "m.vtk")
    mesh = bench.read_vtk_ascii(tmp_path / "m.vtk")
    points, triangles = fan(6)
    np.testing.assert_array_equal(mesh.points[:, :2], points)
    np.testing.assert_array_equal(mesh.points[:, 2], np.arange(len(points)))
    np.testing.assert_array_equal(mesh.triangles, triangles)
    np.testing.assert_array_equal(mesh.edges, [[1, 2], [2, 3]])


def test_mesh_sha256_ignores_the_header(bench: ModuleType, tmp_path: Path) -> None:
    write_mesh(tmp_path / "a.vtk")
    write_mesh(tmp_path / "b.vtk", fields=[("crs", "EPSG:25833")])
    a = bench.read_vtk_ascii(tmp_path / "a.vtk")
    b = bench.read_vtk_ascii(tmp_path / "b.vtk")
    assert re.fullmatch(r"[0-9a-f]{64}", a.sha256)
    assert a.sha256 == b.sha256


def test_read_vtk_ascii_refuses_binary(bench: ModuleType, tmp_path: Path) -> None:
    write_mesh(tmp_path / "bin.vtk", binary=True)
    with pytest.raises(ValueError, match="BINARY"):
        bench.read_vtk_ascii(tmp_path / "bin.vtk")


# ------------------------------------------------------------ run records


def record_dict(**changes: Any) -> dict[str, Any]:
    """A complete, comparable-with-itself `RunRecord` as JSON-able data.

    `changes` are dotted paths (`machine__cpu_brand`) onto new values.
    """
    base: dict[str, Any] = {
        "label": "i21",
        "started": "2026-09-27T10:00:00+00:00",
        "tree": {"commit": "a" * 40, "dirty": False},
        "bench_blob": "b" * 40,
        "child_argv": ["mesh", "--dem", "dem.tif", "--tolerance", "1", "--binary"],
        "build": {
            "no_build": False,
            "type": "Release",
            "cxx_flags_release": "-O3 -DNDEBUG",
            "compiler": "AppleClang 17.0.0",
            "so_sha256": "c" * 64,
        },
        "machine": {
            "cpu_brand": "Apple M4 Pro",
            "p_cores": 10,
            "e_cores": 4,
            "memory_bytes": 48 * 2**30,
            "macos": "26.0",
            "python": "3.14.7",
            "numpy": "2.3.0",
        },
        "power": {"state": "ac", "percent": 100, "raw": AC_WITH_BATTERY},
        "inputs": {
            "dem": "tests/fixtures/dem_archive/7908_3_10m_z33.tif",
            "dem_sha256": "d" * 64,
            "domains": [
                {"name": "tile", "path": None, "sha256": None},
                {"name": "quarter", "path": "quarter.geojson", "sha256": "e" * 64},
            ],
            "tolerance": 1.0,
            "extra_args": [],
        },
        "samples": [
            {
                "domain": "tile",
                "threads": t,
                "repeat": r,
                "refine_s": 1.0 + 0.01 * r,
                "app_s": 2.0,
                "proc_s": 2.5,
                "max_error": 0.99,
                "rounds": 5,
                "inserted": 100,
                "flips": 200,
            }
            for r in range(3)
            for t in (0, 1, 20)
        ],
        "stats": [
            {"domain": d, "threads": t, "n": 3, "median": m, "min": m, "max": m}
            for d in ("tile", "quarter")
            for t, m in ((0, 1.0), (1, 4.0), (20, 1.1))
        ],
        "quality": {
            d: {
                "worst_angle": 20.0,
                "angle_median": 40.0,
                "share_under_1": 0.0,
                "max_degree": 9,
                "within_tolerance": True,
                "delaunay_checked": 1000,
                "delaunay_ambiguous": 3,
                "delaunay_violations": 0,
                "mesh_sha256": "f" * 64,
            }
            for d in ("tile", "quarter")
        },
        "accept_quality": False,
        "threshold_pct": 5.0,
        "verdict": [],
    }
    for dotted, value in changes.items():
        *parents, leaf = dotted.split("__")
        node = base
        for key in parents:
            node = node[key]
        node[leaf] = value
    return base


def make(bench: ModuleType, data: dict[str, Any] | None = None, **changes: Any) -> Any:
    source = copy.deepcopy(data) if data is not None else record_dict()
    for dotted, value in changes.items():
        *parents, leaf = dotted.split("__")
        node = source
        for key in parents:
            node = node[key]
        node[leaf] = value
    return bench.RunRecord.model_validate(source)


def with_median(data: dict[str, Any], domain: str, threads: int, median: float) -> dict[str, Any]:
    out = copy.deepcopy(data)
    for row in out["stats"]:
        if (row["domain"], row["threads"]) == (domain, threads):
            row.update(median=median, min=median, max=median)
    return out


def with_quality(data: dict[str, Any], domain: str, **fields: Any) -> dict[str, Any]:
    out = copy.deepcopy(data)
    out["quality"][domain].update(fields)
    return out


# ------------------------------------------------------------- comparable


def test_identical_records_are_comparable(bench: ModuleType) -> None:
    assert bench.comparable(make(bench), make(bench)) is None


def test_battery_with_battery_is_comparable(bench: ModuleType) -> None:
    battery = {"state": "battery", "percent": 88, "raw": BATTERY}
    assert bench.comparable(make(bench, power=battery), make(bench, power=battery)) is None


def test_machine_facts_outside_the_key_do_not_matter(bench: ModuleType) -> None:
    new = make(bench, machine__macos="26.1", machine__memory_bytes=1, label="other")
    assert bench.comparable(new, make(bench)) is None


@pytest.mark.parametrize(
    ("new_state", "base_state"),
    [
        ("ac", "battery"),
        ("battery", "ac"),
        ("mixed", "mixed"),
        ("unknown", "unknown"),
        ("ac", "mixed"),
        ("unknown", "ac"),
    ],
)
def test_power_mismatch_or_indeterminate_is_not_comparable(
    bench: ModuleType, new_state: str, base_state: str
) -> None:
    new = make(bench, power={"state": new_state, "percent": None, "raw": ""})
    base = make(bench, power={"state": base_state, "percent": None, "raw": ""})
    reason = bench.comparable(new, base)
    assert reason is not None
    assert "power" in reason


@pytest.mark.parametrize(
    ("change", "field"),
    [
        ({"machine__cpu_brand": "Apple M3"}, "cpu_brand"),
        ({"machine__p_cores": 8}, "p_cores"),
        ({"machine__e_cores": 2}, "e_cores"),
        ({"inputs__dem_sha256": "0" * 64}, "dem_sha256"),
        ({"inputs__domains": [{"name": "tile", "path": None, "sha256": None}]}, "domains"),
        ({"inputs__tolerance": 0.5}, "tolerance"),
        ({"inputs__extra_args": ["--start-min-angle", "0"]}, "extra_args"),
    ],
)
def test_mismatch_names_the_field(bench: ModuleType, change: dict[str, Any], field: str) -> None:
    reason = bench.comparable(make(bench, **change), make(bench))
    assert reason is not None
    assert field in reason


# ---------------------------------------------------------------- verdict


def test_no_baseline_when_none_found(bench: ModuleType) -> None:
    v = bench.verdict(make(bench), None)
    assert (v.status, v.exit_code) == ("NO BASELINE", 2)
    assert v.lines[0].startswith("NO BASELINE: ")


def test_ac_against_battery_is_no_baseline_never_a_cross_comparison(bench: ModuleType) -> None:
    battery = make(bench, power={"state": "battery", "percent": 88, "raw": BATTERY})
    slow = make(bench, with_median(record_dict(), "tile", 0, 9.0))
    v = bench.verdict(slow, battery)
    assert (v.status, v.exit_code) == ("NO BASELINE", 2)
    assert v.lines[0].startswith("NO BASELINE: ")
    assert "power" in v.lines[0]


def test_mixed_new_run_is_no_baseline(bench: ModuleType) -> None:
    new = make(bench, power={"state": "mixed", "percent": None, "raw": AC_DESKTOP + BATTERY})
    v = bench.verdict(new, make(bench))
    assert (v.status, v.exit_code) == ("NO BASELINE", 2)
    assert "power" in v.lines[0]


def test_equal_runs_are_accepted(bench: ModuleType) -> None:
    v = bench.verdict(make(bench), make(bench))
    assert (v.status, v.exit_code) == ("ACCEPTED", 0)
    assert not any(line.startswith("REGRESSION") for line in v.lines)


def slower(factor: float) -> dict[str, Any]:
    data = record_dict()
    for row in data["stats"]:
        row.update(median=row["median"] * factor)
    return data


def test_four_percent_slower_is_accepted(bench: ModuleType) -> None:
    v = bench.verdict(make(bench, slower(1.04)), make(bench))
    assert (v.status, v.exit_code) == ("ACCEPTED", 0)


def test_six_percent_slower_is_a_regression_with_its_size(bench: ModuleType) -> None:
    new = make(bench, with_median(record_dict(), "tile", 0, 1.06))
    v = bench.verdict(new, make(bench))
    assert (v.status, v.exit_code) == ("REGRESSION", 1)
    assert len(v.lines) == 1
    assert v.lines[0].startswith("REGRESSION: tile refine_s[t=0] ")
    assert v.lines[0].endswith("(+6.0 %)")


def test_every_regressed_count_gets_a_line(bench: ModuleType) -> None:
    v = bench.verdict(make(bench, slower(1.06)), make(bench))
    assert v.status == "REGRESSION"
    assert len(v.lines) == 6  # two domains x threads 0, 1, 20


def test_threshold_raised_accepts_six_percent(bench: ModuleType) -> None:
    new = make(bench, with_median(record_dict(), "tile", 0, 1.06), threshold_pct=10.0)
    assert bench.verdict(new, make(bench)).status == "ACCEPTED"


def test_threshold_lowered_rejects_four_percent(bench: ModuleType) -> None:
    new = make(bench, with_median(record_dict(), "quarter", 20, 1.1 * 1.04), threshold_pct=3.0)
    v = bench.verdict(new, make(bench))
    assert v.status == "REGRESSION"
    assert v.lines[0].startswith("REGRESSION: quarter refine_s[t=20] ")
    assert v.lines[0].endswith("(+4.0 %)")


@pytest.mark.parametrize(
    ("threshold", "base", "new", "status", "size"),
    [
        (5.0, 20.0, 21.0, "ACCEPTED", None),
        (5.0, 20.0, 21.1, "REGRESSION", "(+5.5 %)"),
        (7.5, 40.0, 43.0, "ACCEPTED", None),
        (7.5, 40.0, 43.2, "REGRESSION", "(+8.0 %)"),
    ],
    ids=["exactly-5-accepted", "5.5-regression", "exactly-7.5-accepted", "8-regression"],
)
def test_time_regression_is_strictly_more_than_the_threshold(
    bench: ModuleType, threshold: float, base: float, new: float, status: str, size: str | None
) -> None:
    """A median exactly `threshold` percent slower is accepted, anything above is not.

    The medians are exact in binary and so are the products `new * 100` and
    `base * (100 + threshold)`, so "exactly 5 %" is exact in the reals. The
    ratio `21.0 / 20.0 - 1` is not: it rounds to 0.050000000000000044, which a
    `pct > threshold` test reads as more than 5 %. The boundary has to be
    decided without that rounding.
    """
    base_record = make(bench, with_median(record_dict(), "tile", 0, base))
    new_record = make(bench, with_median(record_dict(), "tile", 0, new), threshold_pct=threshold)
    v = bench.verdict(new_record, base_record)
    assert v.status == status, v.lines
    if size is None:
        assert not any(line.startswith("REGRESSION") for line in v.lines)
    else:
        assert v.lines == [f"REGRESSION: tile refine_s[t=0] {base:.4f} -> {new:.4f} {size}"]


def test_faster_is_accepted(bench: ModuleType) -> None:
    assert bench.verdict(make(bench, slower(0.5)), make(bench)).status == "ACCEPTED"


def test_one_delaunay_violation_at_equal_time_is_a_regression(bench: ModuleType) -> None:
    new = make(bench, with_quality(record_dict(), "tile", delaunay_violations=1))
    v = bench.verdict(new, make(bench))
    assert (v.status, v.exit_code) == ("REGRESSION", 1)
    assert v.lines == [line for line in v.lines if line.startswith("REGRESSION: ")]
    assert any(line.startswith("REGRESSION: tile delaunay_violations ") for line in v.lines)


def test_tolerance_failure_is_a_regression(bench: ModuleType) -> None:
    new = make(bench, with_quality(record_dict(), "quarter", within_tolerance=False))
    v = bench.verdict(new, make(bench))
    assert v.status == "REGRESSION"
    assert any(line.startswith("REGRESSION: quarter tolerance ") for line in v.lines)


@pytest.mark.parametrize(
    ("fields", "measure"),
    [({"worst_angle": 19.5}, "worst_angle"), ({"max_degree": 10}, "max_degree")],
)
def test_quality_loss_is_a_regression_by_default(
    bench: ModuleType, fields: dict[str, Any], measure: str
) -> None:
    new = make(bench, with_quality(record_dict(), "tile", **fields))
    v = bench.verdict(new, make(bench))
    assert (v.status, v.exit_code) == ("REGRESSION", 1)
    assert any(line.startswith(f"REGRESSION: tile {measure} ") for line in v.lines)


def test_quality_gain_is_accepted(bench: ModuleType) -> None:
    new = make(bench, with_quality(record_dict(), "tile", worst_angle=25.0, max_degree=7))
    assert bench.verdict(new, make(bench)).status == "ACCEPTED"


def test_accept_quality_waives_angle_and_degree_and_still_names_them(bench: ModuleType) -> None:
    data = with_quality(record_dict(), "tile", worst_angle=19.5, max_degree=10)
    v = bench.verdict(make(bench, data, accept_quality=True), make(bench))
    assert (v.status, v.exit_code) == ("ACCEPTED", 0)
    waived = [line for line in v.lines if line.startswith("WAIVED (--accept-quality): ")]
    assert any("tile worst_angle " in line for line in waived)
    assert any("tile max_degree " in line for line in waived)


@pytest.mark.parametrize(
    ("fields", "measure"),
    [
        ({"delaunay_violations": 2}, "delaunay_violations"),
        ({"within_tolerance": False}, "tolerance"),
    ],
)
def test_accept_quality_never_waives_tolerance_or_delaunay(
    bench: ModuleType, fields: dict[str, Any], measure: str
) -> None:
    data = with_quality(record_dict(), "tile", worst_angle=19.5, **fields)
    v = bench.verdict(make(bench, data, accept_quality=True), make(bench))
    assert (v.status, v.exit_code) == ("REGRESSION", 1)
    assert any(line.startswith(f"REGRESSION: tile {measure} ") for line in v.lines)
    assert any(line.startswith("WAIVED (--accept-quality): tile worst_angle ") for line in v.lines)


def test_accept_quality_does_not_waive_time(bench: ModuleType) -> None:
    new = make(bench, with_median(record_dict(), "tile", 0, 1.06), accept_quality=True)
    assert bench.verdict(new, make(bench)).status == "REGRESSION"


# ---------------------------------------------------------- baseline search


def store(root: Path, record: Any, date: str, label: str) -> Path:
    directory = root / date / label
    directory.mkdir(parents=True)
    (directory / "run.json").write_text(record.model_dump_json())
    return directory


def test_find_baseline_takes_the_newest_comparable_ancestor(
    bench: ModuleType, tmp_path: Path
) -> None:
    battery = {"state": "battery", "percent": 80, "raw": BATTERY}

    def at(hour: int, commit: str, **changes: Any) -> Any:
        started = f"2026-09-20T{hour:02d}:00:00+00:00"
        return make(bench, started=started, tree={"commit": commit, "dirty": False}, **changes)

    store(tmp_path, at(1, "old"), "2026-09-20", "oldest")
    chosen = store(tmp_path, at(2, "anc"), "2026-09-20", "chosen")
    store(tmp_path, at(3, "sib"), "2026-09-20", "not-an-ancestor")
    store(tmp_path, at(4, "anc2", power=battery), "2026-09-20", "battery")
    store(tmp_path, at(6, "later"), "2026-09-20", "newer-than-new")
    new = at(5, "new")

    ancestors = {"old", "anc", "anc2", "later"}
    calls: list[tuple[str, str]] = []

    def is_ancestor(candidate: str, descendant: str) -> bool:
        calls.append((candidate, descendant))
        return candidate in ancestors

    found = bench.find_baseline(tmp_path, new, is_ancestor)
    assert found is not None
    directory, record = found
    assert directory == chosen
    assert record.tree.commit == "anc"
    assert all(descendant == "new" for _, descendant in calls)


def test_find_baseline_none_in_an_empty_root(bench: ModuleType, tmp_path: Path) -> None:
    assert bench.find_baseline(tmp_path, make(bench), lambda a, b: True) is None


def test_find_baseline_never_returns_the_record_itself(bench: ModuleType, tmp_path: Path) -> None:
    record = make(bench)
    store(tmp_path, record, "2026-09-27", "i21")
    assert bench.find_baseline(tmp_path, record, lambda a, b: True) is None


# A stored run.json that is not a RunRecord: truncated, empty, and valid JSON
# from a schema that is not this one (the realistic case: RunRecord gains a
# required field and every older record stops validating).
MALFORMED_RUN_JSON = [
    '{"label": "i21", "started": ',
    "",
    json.dumps({"label": "i21", "started": "2026-09-20T02:00:00+00:00"}),
]
MALFORMED_IDS = ["truncated", "empty", "schema-drift"]


def store_malformed(root: Path, text: str, date: str, label: str) -> Path:
    directory = root / date / label
    directory.mkdir(parents=True)
    (directory / "run.json").write_text(text)
    return directory


@pytest.mark.parametrize("text", MALFORMED_RUN_JSON, ids=MALFORMED_IDS)
def test_find_baseline_skips_a_malformed_run_json_with_a_warning_naming_it(
    bench: ModuleType, tmp_path: Path, text: str
) -> None:
    """Pinned in bench-py.md: skip with a warning, never abort the search.

    The file sorts after the valid baseline, so an implementation that stops at
    the first bad file does not find the valid one either.
    """
    valid = store(tmp_path, make(bench, started="2026-09-20T01:00:00+00:00"), "2026-09-20", "good")
    store_malformed(tmp_path, text, "2026-09-20", "zz-broken")
    new = make(bench, started="2026-09-27T10:00:00+00:00")
    with pytest.warns(UserWarning, match="zz-broken"):
        found = bench.find_baseline(tmp_path, new, lambda a, b: True)
    assert found is not None
    assert found[0] == valid


def test_find_baseline_with_only_a_malformed_run_json_is_none(
    bench: ModuleType, tmp_path: Path
) -> None:
    store_malformed(tmp_path, "{", "2026-09-20", "broken")
    with pytest.warns(UserWarning, match="broken"):
        assert bench.find_baseline(tmp_path, make(bench), lambda a, b: True) is None


# --------------------------------------------------------------- evidence


def test_write_evidence_files_and_round_trip(bench: ModuleType, tmp_path: Path) -> None:
    record = make(bench, verdict=["REGRESSION: tile refine_s[t=0] 1.0 -> 1.06 (+6.0 %)"])
    bench.write_evidence(record, tmp_path)
    assert bench.RunRecord.model_validate_json((tmp_path / "run.json").read_text()) == record
    rows = [line.split("\t") for line in (tmp_path / "raw.tsv").read_text().splitlines()]
    assert len(rows) == len(record.samples)
    first = record.samples[0]
    assert rows[0][:3] == [first.domain, str(first.threads), str(first.repeat)]
    assert float(rows[0][3]) == pytest.approx(first.refine_s)
    assert rows[0][4] == "ac"
    readme = (tmp_path / "README.md").read_text()
    assert bench.MARKER in readme.splitlines()
    assert "2.2x" in readme
    assert record.verdict[0] in readme
    assert "--accept-quality" not in readme


def test_accept_quality_is_recorded_in_readme_and_run_json(
    bench: ModuleType, tmp_path: Path
) -> None:
    record = make(
        bench,
        accept_quality=True,
        verdict=["WAIVED (--accept-quality): tile worst_angle 20.0 -> 19.5"],
    )
    bench.write_evidence(record, tmp_path)
    assert json.loads((tmp_path / "run.json").read_text())["accept_quality"] is True
    generated = (tmp_path / "README.md").read_text().split(bench.MARKER)[0]
    assert "--accept-quality" in generated


def test_rerun_keeps_the_prose_below_the_marker(bench: ModuleType, tmp_path: Path) -> None:
    bench.write_evidence(make(bench, verdict=["ACCEPTED"]), tmp_path)
    prose = "\nHand-written by @perf: the fan came on at repeat 3.\n"
    with (tmp_path / "README.md").open("a") as fh:
        fh.write(prose)

    bench.write_evidence(make(bench, verdict=["NO BASELINE: power: ac vs battery"]), tmp_path)
    readme = (tmp_path / "README.md").read_text()
    generated, _, kept = readme.partition(bench.MARKER)
    assert "NO BASELINE: power: ac vs battery" in generated
    assert "ACCEPTED" not in generated
    assert readme.count("Hand-written by @perf") == 1
    assert "Hand-written by @perf: the fan came on at repeat 3." in kept


# ---------------------------------------------------- the Typer app, faked


def completed(bench: ModuleType, stdout: str = "", stderr: str = "", code: int = 0) -> Any:
    return bench.Completed(returncode=code, stdout=stdout, stderr=stderr, wall_s=0.25)


SYSCTL = {
    "machdep.cpu.brand_string": "Apple M4 Pro",
    "hw.perflevel0.physicalcpu": "10",
    "hw.perflevel1.physicalcpu": "4",
    "hw.memsize": str(48 * 2**30),
}


class FakeRunner:
    """Canned output keyed on argv[0]; the child is recognised by `_child`.

    A child whose rasputin argv asks for `--ascii` gets a small mesh written
    to its `--out`, so the parent's quality pass has a file to read.
    """

    def __init__(
        self,
        bench: ModuleType,
        *,
        pmset: Sequence[str] = (AC_WITH_BATTERY, AC_WITH_BATTERY),
        refine_s: float = 1.0,
        child_code: int = 0,
    ) -> None:
        self.bench = bench
        self.pmset = list(pmset)
        self.refine_s = refine_s
        self.child_code = child_code
        self.calls: list[list[str]] = []

    def run(self, argv: Sequence[str], cwd: Path | None = None) -> Any:
        args = [str(a) for a in argv]
        self.calls.append(args)
        program = Path(args[0]).name
        if "_child" in args:
            return self._child(args)
        if program == "pmset":
            return completed(self.bench, self.pmset.pop(0))
        if program == "sysctl":
            return completed(self.bench, SYSCTL.get(args[-1], "0") + "\n")
        if program == "sw_vers":
            return completed(self.bench, "26.0\n")
        if program == "git":
            if "rev-parse" in args:
                return completed(self.bench, "a" * 40 + "\n")
            if "hash-object" in args:
                return completed(self.bench, "b" * 40 + "\n")
            return completed(self.bench)  # status --porcelain, merge-base --is-ancestor
        if program == "cmake":
            return completed(self.bench)
        return completed(self.bench, stderr=f"{program}: not faked", code=127)

    def _child(self, args: list[str]) -> Any:
        rasputin = args[args.index("--") + 1 :]
        if "--ascii" in rasputin:
            write_mesh(Path(rasputin[rasputin.index("--out") + 1]))
        values = {**CHILD_JSON, "refine_s": self.refine_s}
        return completed(self.bench, stderr=bench_line(values) + "\n", code=self.child_code)

    def children(self) -> list[list[str]]:
        return [c for c in self.calls if "_child" in c]


def child_threads(args: list[str]) -> int:
    own = args[: args.index("--")]
    return int(own[own.index("--threads") + 1])


def child_mesh_args(args: list[str]) -> list[str]:
    return args[args.index("--") + 1 :]


@pytest.fixture
def workspace(tmp_path: Path) -> dict[str, Path]:
    dem = tmp_path / "dem.tif"
    dem.write_bytes(b"not decoded: the fake child never reads it")
    places = {
        "dem": dem,
        "tree": tmp_path / "tree",
        "out_root": tmp_path / "evidence",
        "mesh_dir": tmp_path / "meshes",
    }
    for key in ("tree", "out_root", "mesh_dir"):
        places[key].mkdir()
    return places


def run_args(ws: dict[str, Path], label: str, *extra: str) -> list[str]:
    return [
        "run", "--label", label, "--tree", str(ws["tree"]), "--dem", str(ws["dem"]),
        "--domain", "tile", "--threads", "1,2", "--repeats", "2",
        "--out-root", str(ws["out_root"]), "--mesh-dir", str(ws["mesh_dir"]), *extra,
    ]  # fmt: skip


@pytest.fixture
def use_runner(bench: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Callable[[Any], None]:
    def install(fake: Any) -> None:
        monkeypatch.setattr(bench, "make_runner", lambda: fake)

    return install


def run_json(ws: dict[str, Path], label: str) -> dict[str, Any]:
    (path,) = ws["out_root"].glob(f"*/{label}/run.json")
    data: dict[str, Any] = json.loads(path.read_text())
    return data


def test_run_samples_interleaved_with_power_bracketing(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    fake = FakeRunner(bench)
    use_runner(fake)
    result = cli.invoke(bench.app, [*run_args(workspace, "t1"), "--no-build"])
    assert result.exit_code == 2, result.output  # nothing stored yet: NO BASELINE
    assert "NO BASELINE" in result.output

    children = fake.children()
    timing = [c for c in children if "--binary" in child_mesh_args(c)]
    ascii_runs = [c for c in children if "--ascii" in child_mesh_args(c)]
    assert [child_threads(c) for c in timing] == [0, 1, 2, 0, 1, 2]
    assert [child_threads(c) for c in ascii_runs] == [0]
    assert children.index(ascii_runs[0]) == len(children) - 1
    assert all("--pkg" not in c[: c.index("--")] for c in children)  # --no-build: installed
    assert all(child_mesh_args(c)[0] == "mesh" for c in children)

    pmset = [i for i, c in enumerate(fake.calls) if Path(c[0]).name == "pmset"]
    child_at = [i for i, c in enumerate(fake.calls) if "_child" in c]
    assert len(pmset) == 2
    assert pmset[0] < child_at[0] and pmset[1] > child_at[-1]
    assert not any(Path(c[0]).name == "cmake" for c in fake.calls)

    data = run_json(workspace, "t1")
    assert data["build"]["no_build"] is True
    assert data["power"]["state"] == "ac"
    assert len(data["samples"]) == 6
    assert {(s["domain"], s["threads"]) for s in data["stats"]} == {
        ("tile", 0),
        ("tile", 1),
        ("tile", 2),
    }
    assert set(data["quality"]) == {"tile"}
    assert data["quality"]["tile"]["max_degree"] == 6
    assert data["threshold_pct"] == 5.0
    assert data["accept_quality"] is False
    assert data["machine"]["cpu_brand"] == "Apple M4 Pro"
    assert (data["machine"]["p_cores"], data["machine"]["e_cores"]) == (10, 4)
    evidence = next(workspace["out_root"].glob("*/t1"))
    assert {p.name for p in evidence.iterdir()} >= {"run.json", "raw.tsv", "README.md"}


def test_run_passes_extra_mesh_args_and_records_them(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    fake = FakeRunner(bench)
    use_runner(fake)
    args = [*run_args(workspace, "t1"), "--no-build", "--", "--start-min-angle", "0"]
    result = cli.invoke(bench.app, args)
    assert result.exit_code == 2, result.output
    for child in fake.children():
        mesh_args = child_mesh_args(child)
        at = mesh_args.index("--start-min-angle")
        assert mesh_args[at + 1] == "0"
    assert run_json(workspace, "t1")["inputs"]["extra_args"] == ["--start-min-angle", "0"]


def test_run_marks_changed_power_as_mixed(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    use_runner(FakeRunner(bench, pmset=(AC_WITH_BATTERY, BATTERY)))
    result = cli.invoke(bench.app, [*run_args(workspace, "t1"), "--no-build"])
    assert result.exit_code == 2, result.output
    assert "NO BASELINE" in result.output
    assert run_json(workspace, "t1")["power"]["state"] == "mixed"


def test_run_records_threshold_and_accept_quality(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    use_runner(FakeRunner(bench))
    args = [*run_args(workspace, "t1"), "--no-build", "--threshold", "7.5", "--accept-quality"]
    result = cli.invoke(bench.app, args)
    assert result.exit_code == 2, result.output
    data = run_json(workspace, "t1")
    assert data["threshold_pct"] == 7.5
    assert data["accept_quality"] is True
    readme = next(workspace["out_root"].glob("*/t1/README.md")).read_text()
    assert "--accept-quality" in readme


def test_run_refuses_a_debug_cache_before_any_child(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    cache = workspace["tree"] / "build-bench" / "CMakeCache.txt"
    cache.parent.mkdir()
    cache.write_text("CMAKE_BUILD_TYPE:STRING=Debug\n")
    fake = FakeRunner(bench)
    use_runner(fake)
    result = cli.invoke(bench.app, run_args(workspace, "t1"))
    assert result.exit_code == 3, result.output
    assert "Debug" in result.output
    assert fake.children() == []
    assert list(workspace["out_root"].rglob("run.json")) == []


def test_run_stops_on_a_failing_child_naming_it(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    use_runner(FakeRunner(bench, child_code=1))
    result = cli.invoke(bench.app, [*run_args(workspace, "t1"), "--no-build"])
    assert result.exit_code == 3, result.output
    assert "tile" in result.output
    assert list(workspace["out_root"].rglob("run.json")) == []


@pytest.mark.parametrize(
    ("second_refine_s", "exit_code", "status"),
    [(1.04, 0, "ACCEPTED"), (1.06, 1, "REGRESSION")],
)
def test_run_finds_the_previous_run_as_baseline(
    bench: ModuleType,
    workspace: dict[str, Path],
    use_runner: Callable[[Any], None],
    second_refine_s: float,
    exit_code: int,
    status: str,
) -> None:
    use_runner(FakeRunner(bench, refine_s=1.0))
    first = cli.invoke(bench.app, [*run_args(workspace, "base"), "--no-build"])
    assert first.exit_code == 2, first.output

    use_runner(FakeRunner(bench, refine_s=second_refine_s))
    second = cli.invoke(bench.app, [*run_args(workspace, "new"), "--no-build"])
    assert second.exit_code == exit_code, second.output
    assert status in second.output
    assert run_json(workspace, "new")["verdict"]


def test_run_does_not_compare_ac_with_a_stored_battery_run(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    use_runner(FakeRunner(bench, pmset=(BATTERY, BATTERY)))
    cli.invoke(bench.app, [*run_args(workspace, "base"), "--no-build"])
    use_runner(FakeRunner(bench, refine_s=5.0))
    result = cli.invoke(bench.app, [*run_args(workspace, "new"), "--no-build"])
    assert result.exit_code == 2, result.output
    assert "NO BASELINE" in result.output


@pytest.mark.parametrize(
    ("new_data", "exit_code", "status"),
    [
        (record_dict(started="2026-09-27T12:00:00+00:00"), 0, "ACCEPTED"),
        (
            with_median(record_dict(started="2026-09-27T12:00:00+00:00"), "tile", 1, 4.4),
            1,
            "REGRESSION",
        ),
        (
            record_dict(
                started="2026-09-27T12:00:00+00:00",
                power={"state": "battery", "percent": 50, "raw": BATTERY},
            ),
            2,
            "NO BASELINE",
        ),
    ],
    ids=["accepted", "regression", "no-baseline"],
)
def test_compare_rejudges_stored_evidence(
    bench: ModuleType, tmp_path: Path, new_data: dict[str, Any], exit_code: int, status: str
) -> None:
    base_dir, new_dir = tmp_path / "base", tmp_path / "new"
    base_dir.mkdir()
    new_dir.mkdir()
    bench.write_evidence(make(bench), base_dir)
    bench.write_evidence(make(bench, new_data), new_dir)
    result = cli.invoke(bench.app, ["compare", str(new_dir), "--baseline", str(base_dir)])
    assert result.exit_code == exit_code, result.output
    assert status in result.output


# ------------------------------------------------- clean exits, not tracebacks
#
# bench-py.md: a refused or failed run exits 3 and writes no evidence. These pin
# that bad input is refused the same way: exit 3 through `typer.Exit` (so
# `result.exception` is the `SystemExit`, not an escaped `FileNotFoundError` or
# `ValidationError`), a message naming the offending file, and, for `run`,
# refused before any build or child, so a missing input is not found out only
# after the measurement it spoils.


def assert_refused(result: Any, names: str) -> None:
    assert result.exit_code == 3, result.output
    assert isinstance(result.exception, SystemExit), repr(result.exception)
    assert "Traceback" not in result.output
    # The name only: the printed path may be resolved (/private/var vs /var)
    # and a long line may be wrapped.
    assert names in "".join(result.output.split()), result.output


def assert_nothing_started(fake: FakeRunner, ws: dict[str, Path]) -> None:
    assert fake.children() == []
    assert not any(Path(c[0]).name == "cmake" for c in fake.calls)
    assert list(ws["out_root"].rglob("run.json")) == []


def test_run_refuses_a_missing_dem_before_building(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    fake = FakeRunner(bench)
    use_runner(fake)
    workspace["dem"] = workspace["dem"].with_name("missing-dem.tif")
    result = cli.invoke(bench.app, run_args(workspace, "t1"))  # a build is asked for
    assert_refused(result, "missing-dem.tif")
    assert_nothing_started(fake, workspace)


def test_run_refuses_a_missing_domain_before_building(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    fake = FakeRunner(bench)
    use_runner(fake)
    nowhere = workspace["tree"] / "nowhere.geojson"
    args = [*run_args(workspace, "t1"), "--domain", str(nowhere)]  # after --domain tile
    result = cli.invoke(bench.app, args)
    assert_refused(result, "nowhere.geojson")
    assert_nothing_started(fake, workspace)


@pytest.mark.parametrize("text", MALFORMED_RUN_JSON, ids=MALFORMED_IDS)
def test_run_refuses_a_malformed_named_baseline_before_measuring(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None], text: str
) -> None:
    """A `--baseline` the caller named is not skipped: it is refused, and
    before the measurement, which would otherwise be lost with it."""
    fake = FakeRunner(bench)
    use_runner(fake)
    named = store_malformed(workspace["tree"], text, "stored", "bad-baseline")
    result = cli.invoke(bench.app, [*run_args(workspace, "t1"), "--baseline", str(named)])
    assert_refused(result, "bad-baseline")
    assert_nothing_started(fake, workspace)


def test_run_skips_a_malformed_stored_run_and_keeps_its_evidence(
    bench: ModuleType, workspace: dict[str, Path], use_runner: Callable[[Any], None]
) -> None:
    use_runner(FakeRunner(bench, refine_s=1.0))
    first = cli.invoke(bench.app, [*run_args(workspace, "base"), "--no-build"])
    assert first.exit_code == 2, first.output
    (date_dir,) = workspace["out_root"].iterdir()
    store_malformed(workspace["out_root"], "{", date_dir.name, "zz-broken")

    use_runner(FakeRunner(bench, refine_s=1.0))
    with pytest.warns(UserWarning, match="zz-broken"):
        second = cli.invoke(bench.app, [*run_args(workspace, "new"), "--no-build"])
    assert second.exit_code == 0, second.output
    assert "ACCEPTED" in second.output
    assert run_json(workspace, "new")["verdict"][0] == "ACCEPTED"


@pytest.fixture
def stored_pair(bench: ModuleType, tmp_path: Path) -> tuple[Path, Path]:
    base_dir, new_dir = tmp_path / "base", tmp_path / "new"
    base_dir.mkdir()
    new_dir.mkdir()
    bench.write_evidence(make(bench), base_dir)
    bench.write_evidence(make(bench, started="2026-09-27T12:00:00+00:00"), new_dir)
    return base_dir, new_dir


@pytest.mark.parametrize("text", MALFORMED_RUN_JSON, ids=MALFORMED_IDS)
def test_compare_refuses_a_malformed_new_run(
    bench: ModuleType, stored_pair: tuple[Path, Path], text: str
) -> None:
    base_dir, new_dir = stored_pair
    (new_dir / "run.json").write_text(text)
    result = cli.invoke(bench.app, ["compare", str(new_dir), "--baseline", str(base_dir)])
    assert_refused(result, "new/run.json")


def test_compare_refuses_a_directory_without_run_json(
    bench: ModuleType, stored_pair: tuple[Path, Path]
) -> None:
    base_dir, new_dir = stored_pair
    (new_dir / "run.json").unlink()
    result = cli.invoke(bench.app, ["compare", str(new_dir), "--baseline", str(base_dir)])
    assert_refused(result, "new/run.json")


@pytest.mark.parametrize("text", MALFORMED_RUN_JSON, ids=MALFORMED_IDS)
def test_compare_refuses_a_malformed_named_baseline(
    bench: ModuleType, stored_pair: tuple[Path, Path], text: str
) -> None:
    base_dir, new_dir = stored_pair
    (base_dir / "run.json").write_text(text)
    result = cli.invoke(bench.app, ["compare", str(new_dir), "--baseline", str(base_dir)])
    assert_refused(result, "base/run.json")


def test_compare_search_skips_a_malformed_stored_run(
    bench: ModuleType, tmp_path: Path, use_runner: Callable[[Any], None]
) -> None:
    root = tmp_path / "evidence"
    store(root, make(bench), "2026-09-27", "base")
    store_malformed(root, "{", "2026-09-27", "zz-broken")
    new_dir = store(root, make(bench, started="2026-09-27T12:00:00+00:00"), "2026-09-28", "new")
    use_runner(FakeRunner(bench))  # git merge-base --is-ancestor: yes
    with pytest.warns(UserWarning, match="zz-broken"):
        result = cli.invoke(bench.app, ["compare", str(new_dir), "--out-root", str(root)])
    assert result.exit_code == 0, result.output
    assert "ACCEPTED" in result.output


# ------------------------------------------------------- one real child run


@pytest.mark.parametrize("threads", [1, 2])
def test_real_child_wraps_the_real_refine(bench: ModuleType, tmp_path: Path, threads: int) -> None:
    pytest.importorskip("tin_engine._core")
    array = np.random.default_rng(14).uniform(0.0, 50.0, (17, 21)).astype(np.float32)
    dem = tmp_path / "bumpy.tif"
    dem.write_bytes(micro_tiff(array).getvalue())
    out = tmp_path / "m.vtk"
    argv = [
        sys.executable, str(TOOL), "_child", "--threads", str(threads), "--",
        "mesh", "--dem", str(dem), "--tolerance", "1", "--out", str(out), "--binary",
    ]  # fmt: skip
    proc = subprocess.run(argv, capture_output=True, text=True, timeout=300, check=False)
    assert proc.returncode == 0, proc.stderr
    child = bench.parse_child(proc.stderr, run=f"real t={threads}")
    assert child.refine_s > 0.0
    assert child.app_s >= child.refine_s
    assert child.max_error <= 1.0
    assert child.rounds >= 1
    assert out.is_file()
