"""`rasputin mesh --triangles N`: tests A1 to A5, A10, A11 (its CLI half) and A12.

`docs/increments/104-triangle-budget.md`, *10. The red suite*. The budget is a
tolerance found by a search over trial `--tolerance` runs; the mesh written is
the one `--tolerance` writes at the tolerance found. The figures are the
design's, measured on the two fixtures: velhas (the reprojected path, with a
domain) and `bumpy` (the native path, `rough_dem(14)`).
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from typer.testing import Result

from cli_helpers import USAGE, crossing_lines, plain, rough_dem, runner, squashed
from recordread import stats_names
from tin_engine import cli
from tin_engine.cli import app
from vtkread import VtkFile, read_vtk

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
VELHAS = FIXTURES / "velhas"
VELHAS_ARGS = (
    *("--dem", str(VELHAS / "anadem_velhas.tif")),
    *("--domain", str(VELHAS / "catchment.geojson")),
    *("--out-crs", "EPSG:31983"),
)
TIMES = "\N{MULTIPLICATION SIGN}"
#: The record keys a budget run adds; none of them goes into the mesh file.
BUDGET_KEYS = ("triangles_asked", "triangles_cap", "budget_trials", "budget_window_missed")

bumpy = rough_dem(14)
lines = crossing_lines()


@dataclass(frozen=True)
class Run:
    """One successful `rasputin mesh` run: its mesh, record, `--stats` and stderr."""

    vtk: VtkFile
    record: dict[str, Any]
    stats: str
    stderr: str

    @property
    def triangles(self) -> int:
        return len(self.vtk.polygons)


def mesh_run(directory: Path, *args: str, name: str = "x") -> Run:
    out, rec, md = (directory / f"{name}{suffix}" for suffix in (".vtk", ".json", ".md"))
    result = runner.invoke(
        app, ["mesh", *args, "--out", str(out), "--record", str(rec), "--stats", str(md)]
    )
    assert result.exit_code == 0, result.output
    record = json.loads(rec.read_text(encoding="ascii"))
    return Run(read_vtk(out.read_bytes()), record, md.read_text(encoding="utf-8"), result.stderr)


def refused_run(directory: Path, *args: str) -> str:
    """Exit 2 with no traceback; no mesh, record or `--stats` written; the plain output."""
    out, rec, md = (directory / f"refused{suffix}" for suffix in (".vtk", ".json", ".md"))
    result: Result = runner.invoke(
        app, ["mesh", *args, "--out", str(out), "--record", str(rec), "--stats", str(md)]
    )
    output = plain(result.output)
    assert result.exit_code == USAGE, output
    assert "Traceback" not in output
    assert "No such option" not in output, "refused for the wrong reason"
    assert not out.exists() and not rec.exists() and not md.exists()
    return output


def says(output: str, *words: str) -> None:
    """Each word in the output, wherever Rich broke the lines."""
    for word in words:
        assert squashed(word) in squashed(output), f"{word!r} not in {output!r}"


def budget_trials(record: dict[str, Any]) -> list[tuple[float, int]]:
    """`budget_trials` as (tolerance, triangles) pairs: ``1e+06 m: 49; 19.1 m: 157``."""
    pairs = []
    for part in str(record["budget_trials"]).split("; "):
        match = re.fullmatch(r"(\S+) m: (\d+)", part)
        assert match is not None, f"{part!r} in {record['budget_trials']!r}"
        pairs.append((float(match.group(1)), int(match.group(2))))
    return pairs


def says_within_three_figures(run: Run) -> None:
    """The summary's "within ... m" is the tolerance found to three significant
    figures, the figure of "tolerance ... m found in", and its largest difference
    shares that precision."""
    found = f"{run.record['tolerance_m']:.3g}"
    assert repr(run.record["tolerance_m"]) != found, "the tolerance found is already short"
    says(plain(run.stderr), f"tolerance {found} m found in", f"within {found} m of it")
    largest = f"{run.record['max_error_m']:.3g}"
    assert f"{run.record['max_error_m']:.5g}" != largest, "the largest difference is already short"
    says(plain(run.stderr), f"(largest difference {largest} m)")


def timings_note(stats: str) -> str:
    """The `--stats` Timings note: the prose between the heading and its table."""
    after = stats.split("## Timings", 1)[1]
    return squashed(after.split("\n|", 1)[0])


@dataclass
class RefinedCalls:
    count: int = 0


@pytest.fixture
def trials(monkeypatch: pytest.MonkeyPatch) -> RefinedCalls:
    """Counts the trials: every call of `cli._refined`, the one refine of a window."""
    calls = RefinedCalls()
    real: Callable[..., Any] = cli._refined

    def spy(*args: Any, **kwargs: Any) -> Any:
        calls.count += 1
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "_refined", spy)
    return calls


@pytest.fixture(scope="module")
def velhas_7386(tmp_path_factory: pytest.TempPathFactory) -> Run:
    """A1's run, shared by A1 and A4: velhas, `--triangles 7386`, default cap."""
    return mesh_run(tmp_path_factory.mktemp("a1"), *VELHAS_ARGS, "--triangles", "7386")


# ---------------------------------------------------------------- A1: the window


class TestA1Window:
    def test_the_count_lands_in_the_window(self, velhas_7386: Run) -> None:
        assert 7017 <= velhas_7386.triangles <= 7755

    def test_the_record_names_the_tolerance_found_and_the_ask(self, velhas_7386: Run) -> None:
        record = velhas_7386.record
        assert record["tolerance_m"] == pytest.approx(1.03647, abs=0.0005)
        assert record["max_error_m"] <= record["tolerance_m"]
        assert record["triangles_asked"] == 7386
        assert record["triangles_cap"] == 7755
        assert "budget_window_missed" not in record

    def test_budget_trials_names_five_tolerances_from_the_floor(self, velhas_7386: Run) -> None:
        pairs = budget_trials(velhas_7386.record)
        assert len(pairs) == 5
        assert pairs[0] == (1e6, 49)
        assert pairs[-1][1] == velhas_7386.triangles

    def test_the_summary_line_says_the_ask_the_cap_and_the_trials(self, velhas_7386: Run) -> None:
        stderr = plain(velhas_7386.stderr)
        assert f"{velhas_7386.triangles} triangles, about the 7386 asked for (at most 7755)" in (
            stderr
        )
        assert "5 trials" in stderr

    def test_the_summary_says_within_the_tolerance_to_three_figures(self, velhas_7386: Run) -> None:
        """Ola's ruling: a budget run's "within ... m" rounds the tolerance found as
        the first sentence does (`.3g`), not its full `repr`."""
        says_within_three_figures(velhas_7386)

    def test_stats_has_a_row_for_each_budget_key(self, velhas_7386: Run) -> None:
        names = stats_names(velhas_7386.stats)
        assert {"triangles_asked", "triangles_cap", "budget_trials"} <= set(names)

    def test_the_timings_note_says_its_rows_sum_the_five_trials(self, velhas_7386: Run) -> None:
        """Section 4.2: the refine, strip, check and trim rows add every trial's time."""
        note = timings_note(velhas_7386.stats)
        assert squashed("triangle budget") in note, note
        assert squashed("sum") in note and squashed("5 trials") in note, note

    def test_a_tolerance_run_has_no_budget_sentence(self, tmp_path: Path, bumpy: Path) -> None:
        note = timings_note(mesh_run(tmp_path, "--dem", str(bumpy), "--tolerance", "5").stats)
        assert "trial" not in note and squashed("triangle budget") not in note, note

    def test_five_trials_ran_and_none_again(self, tmp_path: Path, trials: RefinedCalls) -> None:
        """The window's trial is the mesh kept: no rerun on a hit."""
        mesh_run(tmp_path, *VELHAS_ARGS, "--triangles", "7386")
        assert trials.count == 5


# ---------------------------------------------------------------- A2: the cap


class TestA2Cap:
    def test_a_two_percent_cap_lands_in_its_window(self, tmp_path: Path) -> None:
        run = mesh_run(tmp_path, *VELHAS_ARGS, "--triangles", "7386", "--max-triangles", "7533")
        assert 7239 <= run.triangles <= 7533
        assert run.record["triangles_cap"] == 7533

    @pytest.mark.parametrize("cap", ["7386", "7000"])
    def test_a_cap_not_over_the_ask_is_refused_naming_both(self, tmp_path: Path, cap: str) -> None:
        output = refused_run(tmp_path, *VELHAS_ARGS, "--triangles", "7386", "--max-triangles", cap)
        says(output, "--max-triangles", "--triangles")

    def test_a_cap_without_an_ask_is_refused(self, tmp_path: Path) -> None:
        output = refused_run(tmp_path, *VELHAS_ARGS, "--max-triangles", "7533")
        says(output, "--max-triangles", "--triangles")

    @pytest.mark.parametrize("asked", ["0", "-5"])
    def test_an_ask_under_one_is_refused(self, tmp_path: Path, asked: str) -> None:
        output = refused_run(tmp_path, *VELHAS_ARGS, "--triangles", asked)
        says(output, "--triangles")


# ---------------------------------------------------------------- A3: the floor


class TestA3Floor:
    def test_a_window_holding_the_floor_writes_it_after_one_trial(self, tmp_path: Path) -> None:
        run = mesh_run(tmp_path, *VELHAS_ARGS, "--triangles", "48")
        assert run.triangles == 49
        assert run.record["tolerance_m"] == pytest.approx(76.598, abs=0.001)
        assert len(budget_trials(run.record)) == 1

    def test_an_ask_under_the_floor_is_refused_naming_it(self, tmp_path: Path) -> None:
        output = refused_run(tmp_path, *VELHAS_ARGS, "--triangles", "40")
        says(output, "49 triangles", "76.6 m", "at least 49")

    def test_the_floor_refusal_is_headed_by_triangles(self, tmp_path: Path) -> None:
        """The refusal names the option asked wrongly, `--triangles`, not `--dem`."""
        output = refused_run(tmp_path, *VELHAS_ARGS, "--triangles", "40")
        says(output, "Invalid value for --triangles")
        assert squashed("Invalid value for --dem") not in squashed(output), output


# ---------------------------------------------------------------- A4: the same mesh


class TestA4SameMesh:
    def test_tolerance_at_the_found_value_writes_the_same_mesh(
        self, tmp_path: Path, velhas_7386: Run
    ) -> None:
        found = velhas_7386.record["tolerance_m"]
        again = mesh_run(tmp_path, *VELHAS_ARGS, "--tolerance", repr(found))
        budget, tolerance = velhas_7386.vtk, again.vtk
        np.testing.assert_array_equal(budget.points, tolerance.points)  # x, y and z
        assert len(budget.polygons) == len(tolerance.polygons)
        np.testing.assert_array_equal(np.array(budget.polygons), np.array(tolerance.polygons))
        assert budget.point_scalars.keys() == tolerance.point_scalars.keys()
        assert budget.field_data.keys() == tolerance.field_data.keys()
        for key, mine in budget.field_data.items():
            theirs = tolerance.field_data[key]
            assert (mine.type_name, mine.components, mine.raw) == (
                theirs.type_name,
                theirs.components,
                theirs.raw,
            ), key
            np.testing.assert_array_equal(mine.values, theirs.values, err_msg=key)
        assert not set(BUDGET_KEYS) & set(budget.field_data)


# ---------------------------------------------------------------- A5: the native path


class TestA5Native:
    def test_bumpy_strided_lands_in_its_window(self, tmp_path: Path, bumpy: Path) -> None:
        run = mesh_run(tmp_path, "--dem", str(bumpy), "--stride", "8", "--triangles", "546")
        assert 519 <= run.triangles <= 573
        assert run.record["triangles_cap"] == 573

    def test_the_summary_says_within_the_tolerance_to_three_figures(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        run = mesh_run(tmp_path, "--dem", str(bumpy), "--stride", "8", "--triangles", "546")
        says_within_three_figures(run)

    def test_a_tolerance_run_says_within_its_tolerance_in_full(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        """Only a budget run rounds: `--tolerance` keeps today's sentence."""
        run = mesh_run(tmp_path, "--dem", str(bumpy), "--tolerance", "4.254655409631253")
        says(plain(run.stderr), "within 4.254655409631253 m of it")

    def test_the_default_cap_is_at_least_one_over_the_ask(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        """⌊1.05 · 10⌋ = 10 is not over 10, so the cap is 11, under the floor of 12."""
        output = refused_run(tmp_path, "--dem", str(bumpy), "--stride", "8", "--triangles", "10")
        says(output, "12 triangles", "at least 12")


# ---------------------------------------------------------------- A10: the options


class TestA10Options:
    def test_triangles_and_tolerance_together_are_refused(self, tmp_path: Path) -> None:
        output = refused_run(tmp_path, *VELHAS_ARGS, "--triangles", "7386", "--tolerance", "1")
        says(output, "one of", "--triangles", "--tolerance")

    def test_triangles_with_pieces_keeps_the_cutting_message(self, tmp_path: Path) -> None:
        output = refused_run(tmp_path, *VELHAS_ARGS, "--triangles", "7386", "--pieces", "2")
        says(output, "cutting needs --tolerance and --domain")

    @pytest.mark.parametrize(
        "flags",
        [
            ("--start-min-angle", "20"),
            ("--start-quality-gain", "1"),
            ("--no-constraint-feet",),
        ],
    )
    def test_refinement_options_need_tolerance_or_triangles(
        self, tmp_path: Path, bumpy: Path, flags: tuple[str, ...]
    ) -> None:
        output = refused_run(tmp_path, "--dem", str(bumpy), *flags)
        says(output, f"{flags[0]} needs --tolerance or --triangles")

    def test_a_domain_needs_tolerance_or_triangles(self, tmp_path: Path) -> None:
        output = refused_run(tmp_path, *VELHAS_ARGS)
        says(output, "--domain needs --tolerance or --triangles")

    def test_refinement_options_are_accepted_with_triangles(self, tmp_path: Path) -> None:
        """Without the start's quality points velhas' floor is 25 triangles
        (measured: `--tolerance 1e6` with these flags), so 24 asks for it."""
        flags = ("--start-min-angle", "0", "--start-quality-gain", "1", "--no-constraint-feet")
        run = mesh_run(tmp_path, *VELHAS_ARGS, "--triangles", "24", *flags)
        assert run.triangles == 25
        assert run.record["start_min_angle_deg"] == 0
        assert run.record["snap_to_lines"] == "off"

    def test_tolerance_near_still_needs_tolerance(
        self, tmp_path: Path, bumpy: Path, lines: Path
    ) -> None:
        output = refused_run(
            tmp_path,
            *("--dem", str(bumpy), "--triangles", "600"),
            *("--tolerance-near", str(lines), "1", "--tolerance-ramp", "0", "40"),
        )
        says(output, "--tolerance-near needs --tolerance")
        assert squashed("--tolerance-near needs --tolerance or") not in squashed(output)

    def test_tolerance_slope_still_needs_tolerance(self, tmp_path: Path, bumpy: Path) -> None:
        output = refused_run(
            tmp_path,
            *("--dem", str(bumpy), "--triangles", "600", "--tolerance-slope", "2", "25", "35"),
        )
        says(output, "--tolerance-slope needs --tolerance")
        assert squashed("--tolerance-slope needs --tolerance or") not in squashed(output)

    def test_the_mcp_mesh_tool_lists_both_options(self) -> None:
        pytest.importorskip("mcp")
        from tin_engine.mcp import schema as mcp_schema

        (tool,) = (d for d in mcp_schema.derive(app) if d.tool_name == "mesh")
        properties = tool.tool().input_schema["properties"]
        assert {"triangles", "max_triangles"} <= set(properties)


# ---------------------------------------------------------------- A11: saturation


class TestA11Saturated:
    def test_a_window_holding_every_node_writes_it_at_tolerance_zero(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        run = mesh_run(tmp_path, "--dem", str(bumpy), "--triangles", "639")
        assert run.triangles == 640
        assert run.record["tolerance_m"] == 0.0
        assert len(budget_trials(run.record)) == 1
        assert "budget_window_missed" not in run.record

    def test_the_timings_note_counts_one_trial(self, tmp_path: Path, bumpy: Path) -> None:
        run = mesh_run(tmp_path, "--dem", str(bumpy), "--triangles", "639")
        note = timings_note(run.stats)
        assert squashed("triangle budget") in note and squashed("1 trial") in note, note

    def test_a_window_above_every_node_writes_every_node_and_says_so(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        run = mesh_run(tmp_path, "--dem", str(bumpy), "--triangles", "1000")
        assert run.triangles == 640
        assert run.record["budget_window_missed"] == "yes"
        stderr = plain(run.stderr)
        assert "640 triangles: every DEM node is a vertex" in stderr, stderr
        assert "about 1000 asked for" in stderr, stderr


# ---------------------------------------------------------------- A12: pieces refused


class TestA12Pieces:
    def test_a_trial_that_would_need_pieces_refuses_before_it_runs(
        self, tmp_path: Path, trials: RefinedCalls
    ) -> None:
        """Velhas' fourth trial, at 1.42 m, would cut the window 2 by 2 under
        a budget of 100 bytes a node; the first three ran."""
        output = refused_run(
            tmp_path, *VELHAS_ARGS, "--triangles", "7386", "--piece-budget", "2054000"
        )
        says(output, "1.42 m", f"2 {TIMES} 2", "2054000 bytes", "--tolerance")
        assert "cutting needs" not in output  # the option checks passed
        assert trials.count == 3

    def test_the_default_piece_budget_lands(self, velhas_7386: Run) -> None:
        assert 7017 <= velhas_7386.triangles <= 7755
