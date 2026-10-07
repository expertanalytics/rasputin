"""The foot rule on every insertion path, through the binding and ``rasputin mesh``.

Increment 20c, PR 20c-1 (``docs/increments/20c-soft-quality.md``, R2 to R5 and
"Tests @tester writes red first", 20c-1, "CLI: the new rows in ``--stats`` and
the record; ``--no-constraint-feet``"). The invariant-critical suite is C++
(``tests/cpp/unit/test_constraint_foot*.cpp``,
``tests/cpp/property/prop_constraint_foot_*.cpp``); this file covers the
binding and the CLI.

R5: no new flag. ``--no-constraint-feet`` turns all three paths off (the
quality start, refinement, the final check); the quality start's feet and the
final check's ``feet`` and ``feet_fallback`` reach ``--stats`` and the record
as rows of their own beside ``points_snapped_to_lines``.

PINNED HERE, where R5 leaves the names open (listed for ``@architect``):

- the binding: ``RefineOutcome.quality_feet``, ``PointRefineOutcome.feet_fallback``
  (``feet`` and ``feet_refused`` it already inherits), and a keyword
  ``constraint_feet`` on ``refine_points`` and ``refine_strip``, default
  ``False``, as ``refine`` has it;
- the three row names: ``start_quality_points_snapped_to_lines`` (the quality
  start's feet), ``final_check_points_snapped_to_lines`` (the final check's
  feet, on either path) and ``final_check_snapped_points_added_anyway`` (its
  fallback: footed points later added where they were);
- the CLI passes ``constraint_feet=True`` to the final check's call unless
  given ``--no-constraint-feet``, on the projected path (``refine_strip``) and
  the reprojected one (``refine_points``).

Every new binding symbol is read inside a test, so a missing one fails its own
test and leaves the rest of the session collecting.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import tin_engine.cli as cli
import tin_engine.edge_strip as edge_strip
import tin_engine.final_check as final_check
from cli_driver import SQUARE, geojson, rough_dem
from recordread import stats_row
from test_cli_mesh_geographic import TARGET, TOLERANCE, domain_4674, geographic_dem
from test_cli_mesh_stats import mesh
from test_core_refine_points import Phase1, scattered, worst_excess
from tin_engine import _core

__all__ = ["domain_4674", "geographic_dem"]  # fixtures, used by name

ROWS = (
    "start_quality_points_snapped_to_lines",
    "final_check_points_snapped_to_lines",
    "final_check_snapped_points_added_anyway",
)

bumpy = rough_dem(20)


@pytest.fixture
def box(tmp_path: Path) -> Path:
    return geojson(tmp_path / "box.geojson", SQUARE)


class Spy:
    """Every call of one ``_core`` function as a module uses it: the keywords
    it was given and what it returned."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch, module: Any, name: str) -> None:
        self.kwargs: list[dict[str, Any]] = []
        self.out: list[Any] = []
        real = getattr(module, name)

        def spy(*args: Any, **kwargs: Any) -> Any:
            self.kwargs.append(kwargs)
            result = real(*args, **kwargs)
            self.out.append(result)
            return result

        monkeypatch.setattr(module, name, spy)


def run(tmp_path: Path, *args: str) -> tuple[str, dict[str, Any]]:
    """(``--stats`` report, record) of one ``rasputin mesh`` run."""
    out, md, rec = tmp_path / "x.vtk", tmp_path / "x.md", tmp_path / "x.json"
    mesh(*args, "--out", str(out), "--stats", str(md), "--record", str(rec))
    return md.read_text(encoding="utf-8"), json.loads(rec.read_text(encoding="ascii"))


# ---------------------------------------------------------------- binding


@pytest.fixture(scope="module")
def phase1() -> Phase1:
    return Phase1(n=25, tolerance=4.0)


class TestBinding:
    def test_refine_reports_the_quality_starts_feet(self, phase1: Phase1) -> None:
        out = phase1.out
        assert isinstance(out.quality_feet, int)
        assert out.quality_feet == 0  # Phase1 runs no quality start

    def test_refine_points_takes_the_switch_and_reports_the_fallback(self, phase1: Phase1) -> None:
        """Scattered points: those in the outermost cells lie within half a
        cell of the ring, so feet go in; J2 holds by the independent oracle."""
        xy, z = scattered(phase1.n, 2, seed=15)
        cp = phase1.store(_core.CheckPoints, xy, z)
        on = _core.refine_points(cp, *phase1.args(), tolerance=0.5, constraint_feet=True)
        assert on.ok(), on.message
        assert isinstance(on.feet_fallback, int)
        assert on.feet > 0
        assert 0 <= on.feet_fallback <= on.feet
        v, t = np.asarray(on.vertices), np.asarray(on.triangles)
        start = np.asarray(phase1.out.vertices)
        excess = worst_excess(xy, z, v, t, np.asarray(on.z), np.asarray(on.valid), 0.5, skip=start)
        assert excess <= 1e-9, f"a check point is {excess} m over tolerance"

    def test_refine_points_is_off_by_default(self, phase1: Phase1) -> None:
        xy, z = scattered(phase1.n, 2, seed=15)
        cp = phase1.store(_core.CheckPoints, xy, z)
        default = _core.refine_points(cp, *phase1.args(), tolerance=0.5)
        off = _core.refine_points(cp, *phase1.args(), tolerance=0.5, constraint_feet=False)
        assert (default.feet, default.feet_fallback) == (0, 0)
        assert np.array_equal(np.asarray(default.vertices), np.asarray(off.vertices))
        assert np.array_equal(np.asarray(default.triangles), np.asarray(off.triangles))

    def test_refine_strip_takes_the_switch(self, phase1: Phase1) -> None:
        view = _core.raster_view(
            np.zeros((phase1.n, phase1.n), np.float32),
            x_min=500000.0,
            y_max=7000000.0,
            delta_x=30.0,
            delta_y=30.0,
        )
        args = phase1.args()
        strip = _core.constraint_check_points(view, args[0], args[4])
        for feet in (False, True):
            out = _core.refine_strip(view, strip, *args, tolerance=1.0, constraint_feet=feet)
            assert out.ok(), out.message
            assert isinstance(out.feet_fallback, int)
            assert out.feet_fallback <= out.feet


# ---------------------------------------------------------------- CLI, projected path


class TestProjectedPath:
    """``--dem`` in the domain's CRS: phase 1 is ``refine``, the final check
    ``refine_strip``."""

    def test_the_default_puts_feet_on_every_path(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        phase = Spy(monkeypatch, cli, "refine")
        final = Spy(monkeypatch, edge_strip, "refine_strip")
        args = ("--dem", str(bumpy), "--domain", str(box), "--tolerance", "1")
        report, record = run(tmp_path, *args)
        assert [k["constraint_feet"] for k in phase.kwargs] == [True]
        assert [k.get("constraint_feet") for k in final.kwargs] == [True]
        (p,), (f,) = phase.out, final.out
        assert stats_row(report, ROWS[0]) == str(p.quality_feet)
        assert stats_row(report, ROWS[1]) == str(f.feet)
        assert stats_row(report, ROWS[2]) == str(f.feet_fallback)
        # Unchanged: refinement's own feet.
        assert stats_row(report, "points_snapped_to_lines") == str(p.feet)
        for name in ROWS:
            assert type(record[name]) is int, name
            assert str(record[name]) == stats_row(report, name)

    def test_no_constraint_feet_turns_all_three_off(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        phase = Spy(monkeypatch, cli, "refine")
        final = Spy(monkeypatch, edge_strip, "refine_strip")
        report, record = run(
            tmp_path,
            *("--dem", str(bumpy), "--domain", str(box), "--tolerance", "1"),
            "--no-constraint-feet",
        )
        assert [k["constraint_feet"] for k in phase.kwargs] == [False]
        assert [k.get("constraint_feet", False) for k in final.kwargs] == [False]
        for name in (*ROWS, "points_snapped_to_lines"):
            assert stats_row(report, name) == "0", name
            assert record[name] == 0, name
        assert stats_row(report, "snap_to_lines") == "off"


# ---------------------------------------------------------------- CLI, reprojected path


class TestReprojectedPath:
    """A geographic DEM with ``--out-crs``: the final check is ``refine_points``
    against the source's nodes."""

    def test_the_default_passes_the_switch_to_refine_points(
        self,
        tmp_path: Path,
        geographic_dem: Path,
        domain_4674: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        final = Spy(monkeypatch, final_check, "refine_points")
        report, _ = run(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
            *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
        )
        assert [k.get("constraint_feet") for k in final.kwargs] == [True]
        (f,) = final.out
        assert stats_row(report, ROWS[1]) == str(f.feet)
        assert stats_row(report, ROWS[2]) == str(f.feet_fallback)

    def test_no_constraint_feet_passes_false(
        self,
        tmp_path: Path,
        geographic_dem: Path,
        domain_4674: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        final = Spy(monkeypatch, final_check, "refine_points")
        report, _ = run(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
            *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
            "--no-constraint-feet",
        )
        assert [k.get("constraint_feet", False) for k in final.kwargs] == [False]
        assert stats_row(report, ROWS[1]) == "0"
        assert stats_row(report, ROWS[2]) == "0"
