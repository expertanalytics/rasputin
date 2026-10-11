"""A tolerance that follows the slope: ``rasputin mesh --tolerance-slope N START END``.

Increment 34 (``docs/increments/34-slope-tolerance.md``, sections 4.5, 5 and
9, tests 11, 12 and 13). DEM nodes on steep ground are held to ``N`` metres:
the general ``--tolerance`` (F) up to START degrees of slope, ``N`` from END
degrees, linear between.

Fixture V2 (section 9): the valley on ``micro_tiff``'s grid, 33 rows by 41
columns, 10 m by 5 m (``test_core_slope_tolerance.valley``); the domain is the
node rectangle, so every node is inside the mesh: 1 353 nodes, 661 of them
held to 2 m by a step at 30 degrees (README item 4).

- test 11: each refusal of section 5, in its words, in its order.
- test 12: the record's four fields; ``max_error_slope_share_of_tolerance <=
  1 + 1e-12``; more vertices on the wall (x from 200 to 440 m) than the
  ``--tolerance F`` mesh; ``N = F`` builds no slope, gives the ``--tolerance
  F`` mesh and says so on stderr; without the flag none of the four fields.
- test 13: on the resampled path (``--out-crs``) the slope is built on the
  target grid and passed to the final check, whose share is the record's.

CHOSEN HERE, where section 4.5 is silent: ``slope_nodes_tightened``'s
percentage is matched as a number to 0.05 of the count's share, not as text.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import tin_engine.cli as cli
import tin_engine.edge_strip as edge_strip
import tin_engine.final_check as final_check
from cli_driver import UTM33, geojson, invoke, refused, write_tiff
from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff
from test_cli_constraint_feet_all_paths import Spy
from test_cli_mesh_geographic import TARGET, domain_4674, geographic_dem
from test_cli_mesh_stats import mesh, seconds
from test_core_slope_tolerance import valley
from vtkread import VtkFile, read_vtk

__all__ = ["domain_4674", "geographic_dem"]  # fixtures, used by name

FIELDS = (
    "tolerance_slope_m",
    "tolerance_slope_deg",
    "slope_nodes_tightened",
    "max_error_slope_share_of_tolerance",
)
ROWS, COLS, DX, DY = 33, 41, 10.0, 5.0
STEEP = ("--tolerance-slope", "2", "30", "30")


@pytest.fixture
def v2(tmp_path: Path) -> Path:
    stream = micro_tiff(valley(ROWS, COLS, DX, DY), scale=(DX, DY, 0.0))
    return write_tiff(tmp_path / "v2.tif", stream)


@pytest.fixture
def box(tmp_path: Path) -> Path:
    """V2's node rectangle, so the mesh holds every node."""
    x1, y1 = TIE_X + (COLS - 1) * DX, TIE_Y - (ROWS - 1) * DY
    ring = [(TIE_X, y1), (x1, y1), (x1, TIE_Y), (TIE_X, TIE_Y)]
    return geojson(tmp_path / "box.geojson", ring, crs=UTM33)


def run(tmp_path: Path, *args: str, out: str = "x") -> tuple[VtkFile, dict[str, Any], str, str]:
    """(mesh, record, output, --stats report) of one ``rasputin mesh`` run."""
    vtk, rec, md = tmp_path / f"{out}.vtk", tmp_path / f"{out}.json", tmp_path / f"{out}.md"
    result = mesh(*args, "--out", str(vtk), "--record", str(rec), "--stats", str(md))
    record = json.loads(rec.read_text(encoding="ascii"))
    return read_vtk(vtk.read_bytes()), record, result.output, md.read_text(encoding="utf-8")


def on_wall(points: np.ndarray) -> int:
    x = np.asarray(points, dtype=np.float64)[:, 0] - TIE_X
    return int(((x >= 200.0) & (x <= 440.0)).sum())


# ---------------------------------------------------------------- test 11: refusals


def steep_args(v2: Path, box: Path, *flags: str) -> tuple[str, ...]:
    return ("--dem", str(v2), "--domain", str(box), *flags)


class TestRefusals:
    """Section 5, in its order; each names --tolerance-slope."""

    def test_without_dem(self, tmp_path: Path) -> None:
        refused(
            tmp_path,
            "mesh",
            *("catchment", "--flat", *STEEP),
            says=("applies only with --dem", "--tolerance-slope"),
        )

    def test_without_tolerance(self, tmp_path: Path, v2: Path) -> None:
        """No --domain, so no other flag asks for --tolerance first."""
        refused(
            tmp_path,
            "mesh",
            *("--dem", str(v2), *STEEP),
            says=("--tolerance-slope needs --tolerance",),
        )

    @pytest.mark.parametrize(
        "values",
        [("nan", "30", "30"), ("inf", "30", "30"), ("2", "nan", "30"), ("2", "30", "inf")],
    )
    def test_a_non_finite_number(
        self, tmp_path: Path, v2: Path, box: Path, values: tuple[str, str, str]
    ) -> None:
        args = steep_args(v2, box, "--tolerance", "10", "--tolerance-slope", *values)
        refused(tmp_path, "mesh", *args, says=("--tolerance-slope", "must be finite, got"))

    @pytest.mark.parametrize("n", ["0", "-1"])
    def test_n_not_above_zero(self, tmp_path: Path, v2: Path, box: Path, n: str) -> None:
        args = steep_args(v2, box, "--tolerance", "10", "--tolerance-slope", n, "30", "30")
        refused(tmp_path, "mesh", *args, says=("--tolerance-slope", f"N must be above 0, got {n}"))

    def test_n_above_f(self, tmp_path: Path, v2: Path, box: Path) -> None:
        args = steep_args(v2, box, "--tolerance", "10", "--tolerance-slope", "12", "30", "30")
        refused(tmp_path, "mesh", *args, says=("--tolerance-slope", "N 12 is above --tolerance 10"))

    def test_start_below_zero(self, tmp_path: Path, v2: Path, box: Path) -> None:
        args = steep_args(v2, box, "--tolerance", "10", "--tolerance-slope", "2", "-5", "30")
        refused(tmp_path, "mesh", *args, says=("--tolerance-slope", "START -5 is below 0"))

    def test_start_above_end(self, tmp_path: Path, v2: Path, box: Path) -> None:
        args = steep_args(v2, box, "--tolerance", "10", "--tolerance-slope", "2", "36", "35")
        refused(tmp_path, "mesh", *args, says=("--tolerance-slope", "START 36 is above END 35"))

    def test_end_not_below_ninety(self, tmp_path: Path, v2: Path, box: Path) -> None:
        args = steep_args(v2, box, "--tolerance", "10", "--tolerance-slope", "2", "30", "90")
        refused(tmp_path, "mesh", *args, says=("--tolerance-slope", "END 90 must be below 90"))

    @pytest.mark.parametrize(
        ("values", "first", "not_said"),
        [
            (("0", "-5", "30"), "N must be above 0", "is below 0"),
            (("12", "-5", "30"), "N 12 is above --tolerance 10", "START -5"),
            (("2", "-5", "95"), "START -5 is below 0", "must be below 90"),
            (("2", "96", "95"), "START 96 is above END 95", "must be below 90"),
        ],
    )
    def test_the_first_broken_bound_in_section_5s_order_is_named(
        self,
        tmp_path: Path,
        v2: Path,
        box: Path,
        values: tuple[str, str, str],
        first: str,
        not_said: str,
    ) -> None:
        args = steep_args(v2, box, "--tolerance", "10", "--tolerance-slope", *values)
        out = refused(tmp_path, "mesh", *args, says=(first,))
        assert not_said not in out, out


# ---------------------------------------------------------------- test 12: runs and record


class TestRecord:
    def test_the_four_fields(self, tmp_path: Path, v2: Path, box: Path) -> None:
        _, record, _, _ = run(tmp_path, *steep_args(v2, box, "--tolerance", "10", *STEEP))
        assert record["tolerance_m"] == 10.0  # still the general value
        assert type(record["tolerance_slope_m"]) is float
        assert record["tolerance_slope_m"] == 2.0
        assert record["tolerance_slope_deg"] == "30 to 30"
        found = re.fullmatch(
            r"661 of 1353 DEM nodes \((\d+(?:\.\d+)?) %\)", record["slope_nodes_tightened"]
        )
        assert found is not None, record["slope_nodes_tightened"]
        assert float(found[1]) == pytest.approx(100.0 * 661 / 1353, abs=0.05)
        share = record["max_error_slope_share_of_tolerance"]
        assert isinstance(share, float)
        assert 0.0 < share <= 1.0 + 1e-12

    def test_a_ramp_says_both_angles(self, tmp_path: Path, v2: Path, box: Path) -> None:
        flags = ("--tolerance", "10", "--tolerance-slope", "2", "25", "35")
        _, record, _, _ = run(tmp_path, *steep_args(v2, box, *flags))
        assert record["tolerance_slope_deg"] == "25 to 35"
        assert record["slope_nodes_tightened"].startswith("676 of 1353 DEM nodes (")

    def test_none_of_them_without_the_flag(self, tmp_path: Path, v2: Path, box: Path) -> None:
        _, record, _, report = run(tmp_path, *steep_args(v2, box, "--tolerance", "10"))
        for name in FIELDS:
            assert name not in record, name
        assert "slope" not in seconds(report)

    def test_the_summary_says_so(self, tmp_path: Path, v2: Path, box: Path) -> None:
        _, _, output, _ = run(tmp_path, *steep_args(v2, box, "--tolerance", "10", *STEEP))
        text = " ".join(output.split())
        assert "Every DEM node inside the mesh is within 10 m of it" in text, text
        steep = "Nodes on steep ground are within the tolerance their slope allows (largest error"
        assert steep in text, text

    def test_the_slope_is_a_stats_phase(self, tmp_path: Path, v2: Path, box: Path) -> None:
        _, _, _, report = run(tmp_path, *steep_args(v2, box, "--tolerance", "10", *STEEP))
        assert "slope" in seconds(report)

    def test_the_slope_reaches_refine_and_the_edge_strip(
        self, tmp_path: Path, v2: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        phase = Spy(monkeypatch, cli, "refine")
        strip = Spy(monkeypatch, edge_strip, "refine_strip")
        _, record, _, _ = run(tmp_path, *steep_args(v2, box, "--tolerance", "10", *STEEP))
        assert [k.get("slope") is not None for k in phase.kwargs] == [True]
        assert [k.get("slope") is not None for k in strip.kwargs] == [True]
        # Section 4.5: on the DEM path the share is the larger of refine's and the strip's.
        assert record["max_error_slope_share_of_tolerance"] == max(
            phase.out[0].max_slope_share, strip.out[0].max_slope_share
        )

    def test_no_slope_without_the_flag(
        self, tmp_path: Path, v2: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        phase = Spy(monkeypatch, cli, "refine")
        strip = Spy(monkeypatch, edge_strip, "refine_strip")
        run(tmp_path, *steep_args(v2, box, "--tolerance", "10"))
        assert [k.get("slope") for k in phase.kwargs + strip.kwargs] == [None, None]


class TestEffect:
    def test_more_vertices_on_the_wall_than_the_general_mesh(
        self, tmp_path: Path, v2: Path, box: Path
    ) -> None:
        steep, _, _, _ = run(tmp_path, *steep_args(v2, box, "--tolerance", "10", *STEEP), out="s")
        plain, _, _, _ = run(tmp_path, *steep_args(v2, box, "--tolerance", "10"), out="p")
        assert on_wall(steep.points) > on_wall(plain.points)

    def test_n_equal_to_f_builds_no_slope_and_is_the_general_mesh(
        self, tmp_path: Path, v2: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        phase = Spy(monkeypatch, cli, "refine")
        flags = ("--tolerance", "10", "--tolerance-slope", "10", "25", "35")
        same, _, output, _ = run(tmp_path, *steep_args(v2, box, *flags), out="same")
        assert [k.get("slope") for k in phase.kwargs] == [None]
        text = " ".join(output.split())
        assert "--tolerance-slope: N equals --tolerance, so slope changes nothing" in text, text
        plain, _, _, _ = run(tmp_path, *steep_args(v2, box, "--tolerance", "10"), out="plain")
        assert_array_equal(same.points, plain.points)
        assert len(same.polygons) == len(plain.polygons)
        for a, b in zip(same.polygons, plain.polygons, strict=True):
            assert_array_equal(a, b)


# ---------------------------------------------------------------- test 13: resampled path


class TestResampled:
    def test_the_slope_is_the_target_grids_and_reaches_the_final_check(
        self,
        tmp_path: Path,
        geographic_dem: Path,
        domain_4674: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """refine runs on the target grid's view, and a slope built on another
        geometry is refused (section 4.5), so a run whose refine and final
        check both got the slope built it on the target grid."""
        phase = Spy(monkeypatch, cli, "refine")
        final = Spy(monkeypatch, final_check, "refine_points")
        args = ("--dem", str(geographic_dem), "--domain", str(domain_4674), "--out-crs", TARGET)
        flags = ("--tolerance", "5", "--tolerance-slope", "1", "10", "20")
        _, record, _, report = run(tmp_path, *args, *flags)
        assert [k.get("slope") is not None for k in phase.kwargs] == [True]
        assert [k.get("slope") is not None for k in final.kwargs] == [True]
        assert record["max_error_slope_share_of_tolerance"] == final.out[-1].max_slope_share
        assert "slope" in seconds(report)


def test_the_flag_is_documented_in_plain_words() -> None:
    code, output = invoke("mesh", "--help")
    assert code == 0, output
    assert "--tolerance-slope" in output
    assert "steep" in " ".join(output.split())
