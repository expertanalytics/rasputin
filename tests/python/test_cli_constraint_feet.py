"""Constraint feet through the binding and ``rasputin mesh``: increment 20b, R9.

`docs/increments/20b-min-insertion-distance.md`. Ola chose C1 (a), C2 (a) and
C3 (a). The invariant-critical suite is C++
(`tests/cpp/property/prop_refinement_constraint_feet.cpp`); this file covers
the binding, the CLI and the real tile.

The binding keyword is the design's ``constraint_feet``, default ``False``, and
``RefineOutcome`` gains ``feet`` and ``feet_refused``. The CLI passes
``constraint_feet=True`` unless given ``--no-constraint-feet``.

Since increment 25 (``docs/increments/25-plain-output.md``, D4) the setting
and its counts are ``--stats`` rows, not stderr or a sentence: ``snap_to_lines``
is ``on`` or ``off``, ``points_snapped_to_lines`` and ``snaps_refused`` are the
outcome's ``feet`` and ``feet_refused``. The flag keeps its name (D8).

F5 through the CLI: ``--no-constraint-feet`` gives increment 20's
``RefineOutcome`` bit for bit. ``INCREMENT_20`` was RECORDED FROM 517e0e3's
BUILD (increment 20's production tree, the extension rebuilt and reinstalled
first) with ``test_refine_golden.digest``, before any production change on this
branch. No commit may update it to agree with new code. The 10 m quarter
circle is recorded too: the design's acceptance says the rule never fires
there, so the default must match it as well.

T-real is logged, not thresholded, except the tolerance and, at 1 m, no
triangle under 0.1 degrees.
"""

from __future__ import annotations

import re
import time
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import tin_engine.cli as cli
from geotiff_fixtures import KARTVERKET, micro_tiff, needs_codecs
from test_cli_mesh_dem import USAGE, invoke, write_tiff
from test_cli_mesh_domain import SQUARE, geojson, quarter_circle
from test_cli_mesh_refine import file_field, min_angles_degrees, stats_row
from test_cli_mesh_stats import mesh
from test_refine_golden import GOLDEN, digest, refined
from tin_engine import _core
from vtkread import VtkFile, read_vtk

ROWS, COLS = 17, 21
FEET = re.compile(r"(\d+) constraint feet, (\d+) feet refused")  # the pre-25 stderr report

# Recorded from increment 20's CLI; see the module docstring.
INCREMENT_20 = {
    ("tile", "1"): "a8e8720d37147f3e18d1702362c94418aec87a16cbbd2491d83315f94627060f",
    ("quarter_circle", "1"): "19a0ed5a62a47ef6af1b32b87baebb9b8a42d0c2b445d4d18dc8114600062ed1",
    ("quarter_circle", "10"): "e59df8bfb57f91133a1bc4a9fbefc6cbfc98e200606d331ebc05857dd8eb1a40",
}


@pytest.fixture
def bumpy(tmp_path: Path) -> Path:
    array = np.random.default_rng(20).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))


@pytest.fixture
def box(tmp_path: Path) -> Path:
    return geojson(tmp_path / "box.geojson", SQUARE)


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    """Every keyword set ``cli.refine`` is called with, passed through."""
    seen: list[dict[str, Any]] = []
    real = cli.refine

    def spy(*args: Any, **kwargs: Any) -> Any:
        seen.append(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "refine", spy)
    return seen


def run(tmp_path: Path, *args: str) -> tuple[VtkFile, str, str]:
    """(mesh, ``--stats`` report, stderr)."""
    out, md = tmp_path / "x.vtk", tmp_path / "x.md"
    result = mesh(*args, "--out", str(out), "--stats", str(md))
    return read_vtk(out.read_bytes()), md.read_text(encoding="utf-8"), result.stderr


# ---------------------------------------------------------------- binding


class TestBinding:
    """R9: the keyword, its default, and the two outcome fields."""

    @needs_codecs
    @pytest.mark.parametrize("case", sorted(GOLDEN))
    def test_false_given_explicitly_matches_the_digest(
        self, case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        real = _core.refine

        def with_feet_off(*args: Any, **kwargs: Any) -> _core.RefineOutcome:
            return real(*args, **kwargs, constraint_feet=False)

        monkeypatch.setattr(_core, "refine", with_feet_off)
        out = refined(case, tmp_path, 1, min_angle_deg=0.0)
        assert digest(out) == GOLDEN[case]
        assert (out.feet, out.feet_refused) == (0, 0)

    @needs_codecs
    def test_the_default_is_off(self, tmp_path: Path) -> None:
        out = refined("quarter_circle", tmp_path, 1, min_angle_deg=0.0)
        assert isinstance(out.feet, int)
        assert isinstance(out.feet_refused, int)
        assert (out.feet, out.feet_refused) == (0, 0)

    @needs_codecs
    def test_true_puts_feet_on_the_quarter_circle(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        real = _core.refine

        def with_feet_on(*args: Any, **kwargs: Any) -> _core.RefineOutcome:
            return real(*args, **kwargs, constraint_feet=True)

        monkeypatch.setattr(_core, "refine", with_feet_on)
        out = refined("quarter_circle", tmp_path, 1, min_angle_deg=25.0)
        assert out.feet > 0  # M1: six, all on one boundary segment
        assert out.feet_refused >= 0
        assert out.inserted >= out.feet  # R9: `inserted` counts feet too
        assert out.max_error <= 1.0


# ---------------------------------------------------------------- CLI


class TestTheFlag:
    """R9: on by default in the CLI, ``--no-constraint-feet`` turns it off."""

    def test_the_default_is_on(
        self, tmp_path: Path, bumpy: Path, box: Path, calls: list[dict[str, Any]]
    ) -> None:
        _, report, _ = run(tmp_path, "--dem", str(bumpy), "--domain", str(box), "--tolerance", "1")
        assert [c["constraint_feet"] for c in calls] == [True]
        assert stats_row(report, "snap_to_lines") == "on"

    def test_the_default_is_on_for_a_stride_start(
        self, tmp_path: Path, bumpy: Path, calls: list[dict[str, Any]]
    ) -> None:
        _, report, _ = run(tmp_path, "--dem", str(bumpy), "--tolerance", "1")
        assert [c["constraint_feet"] for c in calls] == [True]
        assert stats_row(report, "snap_to_lines") == "on"

    def test_no_constraint_feet_is_off(
        self, tmp_path: Path, bumpy: Path, box: Path, calls: list[dict[str, Any]]
    ) -> None:
        _, report, stderr = run(
            tmp_path,
            "--dem",
            str(bumpy),
            "--domain",
            str(box),
            "--tolerance",
            "1",
            "--no-constraint-feet",
        )
        assert [c["constraint_feet"] for c in calls] == [False]
        assert stats_row(report, "snap_to_lines") == "off"
        assert stats_row(report, "points_snapped_to_lines") == "0"
        assert stats_row(report, "snaps_refused") == "0"
        assert FEET.search(stderr) is None, stderr  # D7: the counts left stderr

    def test_the_help_names_the_flag(self) -> None:
        code, output = invoke("--help")
        assert code == 0, output
        assert "--no-constraint-feet" in output

    def test_the_setting_is_not_in_the_file(self, tmp_path: Path, bumpy: Path, box: Path) -> None:
        """D2: the file carries what a user of the mesh needs; the setting is --stats'."""
        vtk, _, _ = run(tmp_path, "--dem", str(bumpy), "--domain", str(box), "--tolerance", "1")
        for name in ("snap_to_lines", "constraint_feet", "elevation_source"):
            assert name not in vtk.field_data, name


class TestRefusals:
    """R9: ``--no-constraint-feet`` is refused where it could change nothing."""

    def test_without_tolerance(self, tmp_path: Path, bumpy: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke("--dem", str(bumpy), "--no-constraint-feet", "--out", str(target))
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its use, not unknown
        assert "--no-constraint-feet needs --tolerance" in output
        assert not target.exists()

    def test_without_dem(self, tmp_path: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke("catchment", "--flat", "--no-constraint-feet", "--out", str(target))
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its use, not unknown
        assert "applies only with --dem" in output
        assert "--no-constraint-feet" in output
        assert not target.exists()


class TestReport:
    """R9, as increment 25 words it: ``--stats`` carries the outcome's counts."""

    def test_stats_rows_match_the_outcome(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        seen: list[_core.RefineOutcome] = []
        real = cli.refine

        def spy(*args: Any, **kwargs: Any) -> _core.RefineOutcome:
            out = real(*args, **kwargs)
            seen.append(out)
            return out

        monkeypatch.setattr(cli, "refine", spy)
        _, report, _ = run(tmp_path, "--dem", str(bumpy), "--domain", str(box), "--tolerance", "1")
        (out,) = seen
        assert stats_row(report, "points_snapped_to_lines") == str(out.feet)
        assert stats_row(report, "snaps_refused") == str(out.feet_refused)


# ---------------------------------------------------------------- F5 and T-real


def _cli_outcome(
    case: str, tolerance: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *extra: str
) -> _core.RefineOutcome:
    """The RefineOutcome ``rasputin mesh`` itself computes for ``case``."""
    seen: list[_core.RefineOutcome] = []
    real = cli.refine

    def spy(*args: Any, **kwargs: Any) -> _core.RefineOutcome:
        out = real(*args, **kwargs)
        seen.append(out)
        return out

    monkeypatch.setattr(cli, "refine", spy)
    args = ["--dem", str(KARTVERKET), "--tolerance", tolerance, "--out", str(tmp_path / "x.vtk")]
    if case == "quarter_circle":
        args += ["--domain", str(geojson(tmp_path / "quarter.geojson", quarter_circle()))]
    code, output = invoke(*args, *extra)
    assert code == 0, output
    (out,) = seen
    return out


@needs_codecs
@pytest.mark.parametrize("case", [("tile", "1"), ("quarter_circle", "1")])
def test_no_constraint_feet_matches_increment_20(
    case: tuple[str, str], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _cli_outcome(*case, tmp_path, monkeypatch, "--no-constraint-feet")
    assert digest(out) == INCREMENT_20[case]


@needs_codecs
def test_the_default_changes_the_1m_quarter_circle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The converse, so the test above cannot pass with a rule that never runs."""
    out = _cli_outcome("quarter_circle", "1", tmp_path, monkeypatch)
    assert digest(out) != INCREMENT_20[("quarter_circle", "1")]
    assert out.feet > 0


@needs_codecs
def test_the_default_leaves_the_10m_quarter_circle_alone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Acceptance: at 10 m the rule never fires and the mesh is increment 20's."""
    out = _cli_outcome("quarter_circle", "10", tmp_path, monkeypatch)
    assert out.feet == 0
    assert digest(out) == INCREMENT_20[("quarter_circle", "10")]


class TestRealTile:
    """T-real: the quarter circle at 1 m and 10 m, M1's columns, logged."""

    @needs_codecs
    def test_the_quarter_circle(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
        domain = geojson(tmp_path / "quarter.geojson", quarter_circle())
        for tolerance in ("1", "10"):
            out, md = tmp_path / f"q_{tolerance}.vtk", tmp_path / f"q_{tolerance}.md"
            began = time.perf_counter()
            mesh(
                "--dem",
                str(KARTVERKET),
                "--domain",
                str(domain),
                "--tolerance",
                tolerance,
                "--out",
                str(out),
                "--stats",
                str(md),
            )
            elapsed = time.perf_counter() - began
            vtk = read_vtk(out.read_bytes())
            achieved = float(file_field(vtk, "max_error_m"))
            assert achieved <= float(tolerance)
            angles = min_angles_degrees(vtk)
            under = int((angles < 0.1).sum())
            report = md.read_text(encoding="utf-8")
            feet = tuple(stats_row(report, n) for n in ("points_snapped_to_lines", "snaps_refused"))
            with capsys.disabled():
                print(
                    f"\nT-real feet {tolerance} m: {len(vtk.polygons)} triangles, "
                    f"worst {angles.min():.4f} deg, {under} under 0.1 deg, "
                    f"{int((angles < 1.0).sum())} under 1 deg, achieved {achieved} m, "
                    f"feet {feet}, {elapsed:.1f} s"
                )
            if tolerance == "1":
                assert under == 0  # M1: 3 with the rule off, 0 with it on
