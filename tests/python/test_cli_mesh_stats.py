"""`rasputin mesh --stats PATH`: increment 17, C1 to C6.

`docs/increments/17-mesh-stats.md` R2, R4 and R5. Ola chose C1 (a): ``--stats -``
writes the report to stdout, after the path line(s). Ola chose C2 (a): the
Refinement table has a ``carved`` column.

Stdout and stderr are read apart (``Result.stdout``), because the report is
Markdown and the stdout contract is exact lines; the shared ``invoke`` helper
merges the two streams and collapses whitespace, which is right for refusals
and wrong here. Refusals go through ``invoke``.

Names and wording as ``test_stats.py`` pins them. Timings are checked only as
parseable, non-negative and consistent (C6): no thresholds, no comparison of
one run with another.

Not invariant-critical, no mutation round.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner, Result

from geotiff_fixtures import micro_tiff
from test_cli_mesh_dem import USAGE, invoke, write_tiff
from test_cli_mesh_domain import HOLE, SQUARE, geojson
from tin_engine.cli import app
from vtkread import read_vtk

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

ReportOf = Callable[..., tuple[str, Result, Path]]

TIMES = "\u00d7"  # the report's multiplication sign
HEADING = "# rasputin mesh — statistics"
ROWS, COLS = 17, 21


@pytest.fixture
def bumpy(tmp_path: Path) -> Path:
    array = np.random.default_rng(17).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))


@pytest.fixture
def square(tmp_path: Path) -> Path:
    return geojson(tmp_path / "square.geojson", SQUARE, (HOLE,))


#: The four kinds of run, as the arguments before ``--out``.
KINDS = ("fixture", "dem", "tolerance", "domain")


def kind_args(kind: str, bumpy: Path, square: Path) -> list[str]:
    return {
        "fixture": ["catchment", "--flat"],
        "dem": ["--dem", str(bumpy)],
        "tolerance": ["--dem", str(bumpy), "--tolerance", "1"],
        "domain": ["--dem", str(bumpy), "--domain", str(square), "--tolerance", "1"],
    }[kind]


def mesh(*args: str) -> Result:
    result = runner.invoke(app, ["mesh", *args])
    assert result.exit_code == 0, result.output
    return result


def same_file(printed: str, path: Path) -> bool:
    return Path(printed).resolve() == path.resolve()


def section(report: str, heading: str) -> str:
    """The text under ``## heading`` up to the next ``## ``."""
    match = re.search(rf"^## {re.escape(heading)}\n(.*?)(?=^## |\Z)", report, re.M | re.S)
    assert match is not None, f"no section {heading!r} in\n{report}"
    return match.group(1)


def table(text: str) -> dict[str, list[str]]:
    """Every table row in ``text`` by its first cell, header and rule excluded."""
    rows: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line.startswith("|") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        rows[cells[0]] = cells[1:]
    return rows


def sections(report: str) -> list[str]:
    return [line[3:] for line in report.splitlines() if line.startswith("## ")]


def seconds(report: str) -> dict[str, float]:
    return {
        name.strip("*"): float(cells[0].strip("*"))
        for name, cells in table(section(report, "Timings")).items()
        if name != "phase"
    }


@pytest.fixture
def report_of(tmp_path: Path, bumpy: Path, square: Path) -> ReportOf:
    """Run one kind with ``--stats`` to a file; return (report, result, out)."""

    def run(kind: str, *extra: str) -> tuple[str, Result, Path]:
        out, md = tmp_path / f"{kind}.vtk", tmp_path / f"{kind}.md"
        args = kind_args(kind, bumpy, square)
        result = mesh(*args, "--out", str(out), "--stats", str(md), *extra)
        return md.read_text(encoding="utf-8"), result, out

    return run


class TestOffByDefault:
    """C1: without ``--stats`` stdout is exactly the path line(s), no report."""

    @pytest.mark.parametrize("kind", KINDS)
    def test_stdout_is_the_path_and_no_markdown_appears(
        self, tmp_path: Path, bumpy: Path, square: Path, kind: str
    ) -> None:
        out = tmp_path / "out" / "x.vtk"
        out.parent.mkdir()
        result = mesh(*kind_args(kind, bumpy, square), "--out", str(out))
        (line,) = result.stdout.splitlines()
        assert same_file(line, out)
        assert result.stdout == f"{line}\n"
        assert HEADING not in result.stdout + result.stderr
        assert list(tmp_path.rglob("*.md")) == []

    def test_ply_with_edges_prints_two_paths(self, tmp_path: Path) -> None:
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        result = mesh("catchment", "--flat", "--out", str(surface), "--out-edges", str(edges))
        lines = result.stdout.splitlines()
        assert len(lines) == 2 and same_file(lines[0], surface) and same_file(lines[1], edges)
        assert list(tmp_path.rglob("*.md")) == []


class TestStatsToAFile:
    """C2."""

    @pytest.mark.parametrize("kind", KINDS)
    def test_the_report_is_written_and_its_path_follows_the_mesh(
        self, tmp_path: Path, bumpy: Path, square: Path, kind: str
    ) -> None:
        plain_out, stats_out = tmp_path / "a" / "x.vtk", tmp_path / "b" / "x.vtk"
        plain_out.parent.mkdir()
        stats_out.parent.mkdir()
        md = tmp_path / "b" / "x.md"
        args = kind_args(kind, bumpy, square)
        mesh(*args, "--out", str(plain_out))
        result = mesh(*args, "--out", str(stats_out), "--stats", str(md))
        lines = result.stdout.splitlines()
        assert len(lines) == 2, result.stdout
        assert same_file(lines[0], stats_out) and same_file(lines[1], md)
        assert md.read_text(encoding="utf-8").startswith(HEADING + "\n")
        assert stats_out.read_bytes() == plain_out.read_bytes()

    def test_ply_with_edges_then_the_report(self, tmp_path: Path) -> None:
        surface, edges, md = tmp_path / "s.ply", tmp_path / "e.ply", tmp_path / "s.md"
        result = mesh(
            "catchment", "--flat", "--out", str(surface), "--out-edges", str(edges),
            "--stats", str(md),
        )
        lines = result.stdout.splitlines()
        assert [same_file(a, b) for a, b in zip(lines, (surface, edges, md), strict=True)] == [
            True, True, True,
        ]
        rows = table(section(md.read_text(encoding="utf-8"), "Sizes"))
        assert "s.ply" in rows and "e.ply" in rows

    def test_out_parent_applies_to_the_report(self, tmp_path: Path) -> None:
        out, md = tmp_path / "x.vtk", tmp_path / "x.md"
        mesh("catchment", "--flat", "--out", str(out), "--stats", str(md),
             "--out-parent", str(tmp_path))
        assert md.is_file()

    def test_the_stderr_summary_is_unchanged(self, tmp_path: Path, bumpy: Path) -> None:
        plain_run = mesh("--dem", str(bumpy), "--tolerance", "1", "--out", str(tmp_path / "a.vtk"))
        stats_run = mesh(
            "--dem", str(bumpy), "--tolerance", "1", "--out", str(tmp_path / "b.vtk"),
            "--stats", str(tmp_path / "b.md"),
        )
        assert stats_run.stderr == plain_run.stderr


class TestStatsToStdout:
    """C3, as Ola chose it: C1 (a), stdout after the path line(s)."""

    @pytest.mark.parametrize("kind", KINDS)
    def test_the_report_follows_the_path_line(
        self, tmp_path: Path, bumpy: Path, square: Path, kind: str,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        monkeypatch.chdir(tmp_path)
        out = tmp_path / "x.vtk"
        result = mesh(*kind_args(kind, bumpy, square), "--out", str(out), "--stats", "-")
        first, rest = result.stdout.split("\n", 1)
        assert same_file(first, out)
        assert rest.startswith(HEADING + "\n")
        assert "## Timings" in rest
        assert HEADING not in result.stderr
        assert not (tmp_path / "-").exists()
        assert list(tmp_path.rglob("*.md")) == []

    def test_after_both_ply_paths(self, tmp_path: Path) -> None:
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        result = mesh(
            "catchment", "--flat", "--out", str(surface), "--out-edges", str(edges),
            "--stats", "-",
        )
        lines = result.stdout.splitlines()
        assert same_file(lines[0], surface) and same_file(lines[1], edges)
        assert lines[2] == HEADING


class TestSectionsPerKind:
    """C4."""

    def test_a_fixture(self, report_of: ReportOf) -> None:
        report, _, out = report_of("fixture")
        assert sections(report) == ["Sizes", "Quality (plan view, x/y)", "Timings"]
        sizes = table(section(report, "Sizes"))
        for absent in ("DEM nodes", "domain vertices", "vertices without data dropped"):
            assert absent not in sizes
        vtk = read_vtk(out.read_bytes())
        assert sizes["output vertices"] == [str(len(vtk.points))]
        assert sizes["output triangles"] == [str(len(vtk.polygons))]
        assert out.name in sizes
        phases = seconds(report)
        for row in (
            "start mesh: build", "start mesh: node", "start mesh: triangulate",
            "start mesh: constraint edges", "write: encode", "write: disk", "other", "total",
        ):
            assert row in phases, row
        for absent in ("decode", "sample", "refine", "trim"):
            assert absent not in phases, absent

    def test_a_dem_without_tolerance(self, report_of: ReportOf) -> None:
        report, _, _ = report_of("dem")
        assert "## Refinement" not in report
        sizes = table(section(report, "Sizes"))
        assert sizes["DEM nodes"] == [f"{ROWS} {TIMES} {COLS} (10 {TIMES} 5 m)"]
        assert "vertices without data dropped" in sizes
        assert "domain vertices" not in sizes
        phases = seconds(report)
        assert "decode" in phases and "sample" in phases and "trim" in phases
        assert not any(name.startswith("refine") for name in phases)

    def test_a_tolerance_run(self, report_of: ReportOf) -> None:
        report, result, _ = report_of("tolerance")
        assert sections(report) == [
            "Sizes", "Quality (plan view, x/y)", "Refinement", "Timings",
        ]
        refinement = table(section(report, "Refinement"))
        header = refinement.pop("tolerance")
        ((tolerance, cells),) = refinement.items()
        row = dict(zip(header, cells, strict=True))
        assert tolerance == "1 m"
        match = re.search(r"(\d+) rounds, (\d+) points inserted, (\d+) flips", result.stderr)
        assert match is not None, result.stderr
        assert (row["rounds"], row["inserted"], row["flips"]) == match.groups()
        assert row["carved"] == "0"  # the bumpy DEM has no NoData
        assert "Threads:" in section(report, "Timings")
        phases = seconds(report)
        for name in (
            "refine", "refine: legalise start", "refine: start quality",
            "refine: scan (parallel)", "refine: split + flip (serial)",
            "refine: setup + output", "trim",
        ):
            assert name in phases, name
        assert "sample" not in phases
        assert "Sub-rows sum to their parent" in report

    def test_a_domain_run(self, report_of: ReportOf) -> None:
        report, _, _ = report_of("domain")
        sizes = table(section(report, "Sizes"))
        assert sizes["domain vertices"] == ["8 (1 ring, 1 hole)"]
        assert "domain read" in seconds(report)


class TestRefusals:
    """C5: refused, exit 2, and neither mesh nor report is written."""

    @staticmethod
    def nothing_written(tmp_path: Path) -> list[Path]:
        return [p for p in tmp_path.rglob("*") if p.suffix in {".vtk", ".ply", ".md"}]

    def test_stats_resolving_to_out(self, tmp_path: Path) -> None:
        out = tmp_path / "x.vtk"
        code, output = invoke("catchment", "--flat", "--out", str(out),
                              "--stats", str(tmp_path / "." / "x.vtk"))
        assert code == USAGE, output
        assert "--stats" in output and "overwrite" in output
        assert self.nothing_written(tmp_path) == []

    def test_stats_resolving_to_out_edges(self, tmp_path: Path) -> None:
        code, output = invoke(
            "catchment", "--flat", "--out", str(tmp_path / "s.ply"),
            "--out-edges", str(tmp_path / "e.ply"), "--stats", str(tmp_path / "e.ply"),
        )
        assert code == USAGE, output
        assert "--stats" in output and "overwrite" in output
        assert self.nothing_written(tmp_path) == []

    def test_stats_outside_out_parent(self, tmp_path: Path) -> None:
        inside = tmp_path / "inside"
        inside.mkdir()
        code, output = invoke(
            "catchment", "--flat", "--out", str(inside / "x.vtk"),
            "--out-parent", str(inside), "--stats", str(tmp_path / "x.md"),
        )
        assert code == USAGE, output
        assert "outside the permitted parent" in output
        assert self.nothing_written(tmp_path) == []

    def test_stats_through_a_symlink(self, tmp_path: Path) -> None:
        target = tmp_path / "real.md"
        target.write_text("keep\n")
        link = tmp_path / "link.md"
        link.symlink_to(target)
        code, output = invoke("catchment", "--flat", "--out", str(tmp_path / "x.vtk"),
                              "--stats", str(link))
        assert code == USAGE, output
        assert "symlink" in output
        assert target.read_text() == "keep\n"
        assert not (tmp_path / "x.vtk").exists()

    def test_a_refused_mesh_run_writes_no_report(self, tmp_path: Path, bumpy: Path) -> None:
        md = tmp_path / "x.md"
        code, output = invoke("--dem", str(bumpy), "--tolerance", "-1",
                              "--out", str(tmp_path / "x.vtk"), "--stats", str(md))
        assert code == USAGE, output
        assert "must be finite and >= 0" in output  # the tolerance, not an unknown --stats
        assert not md.exists() and not (tmp_path / "x.vtk").exists()


class TestTimingsAreSane:
    """C6: parseable, non-negative, consistent. No thresholds."""

    #: Each printed value is rounded to 0.001 s, so a sum of k of them may
    #: differ from its printed total by up to k / 2 ms.
    HALF_MS = 0.0005

    @pytest.mark.parametrize("kind", KINDS)
    def test_every_seconds_and_share_cell(self, report_of: ReportOf, kind: str) -> None:
        report, _, _ = report_of(kind)
        rows = table(section(report, "Timings"))
        rows.pop("phase")
        for name, (secs, share) in rows.items():
            assert float(secs.strip("*")) >= 0.0, name
            assert 0.0 <= float(share.strip("*% ")) <= 100.0, name
        match = re.search(r"Statistics computed in ([0-9.]+) s, not included above\.", report)
        assert match is not None and float(match.group(1)) >= 0.0

    @pytest.mark.parametrize("kind", KINDS)
    def test_top_level_rows_and_other_make_the_total(self, report_of: ReportOf, kind: str) -> None:
        report, _, _ = report_of(kind)
        phases = seconds(report)
        total = phases.pop("total")
        parents = set(phases)
        top = [s for name, s in phases.items() if name.split(": ", 1)[0] not in parents - {name}]
        assert sum(top) == pytest.approx(total, abs=self.HALF_MS * (len(top) + 1))

    @pytest.mark.parametrize("kind", ["tolerance", "domain"])
    def test_the_refine_sub_rows_sum_to_at_most_refine(
        self, report_of: ReportOf, kind: str
    ) -> None:
        report, _, _ = report_of(kind)
        phases = seconds(report)
        subs = [s for name, s in phases.items() if name.startswith("refine: ")]
        assert len(subs) == 5  # increment 20 adds "refine: start quality"
        assert sum(subs) <= phases["refine"] + self.HALF_MS * 6
        four = sum(s for name, s in phases.items()
                   if name.startswith("refine: ") and name != "refine: setup + output")
        assert four <= phases["refine"] + self.HALF_MS * 5
