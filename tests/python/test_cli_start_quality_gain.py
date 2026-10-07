"""The soft criterion and the line split, through the binding and ``rasputin mesh``.

Increment 20c, PR 20c-2 (``docs/increments/20c-soft-quality.md``, R7, R8 and
"Tests @tester writes red first", 20c-2: "T-P3: gain -1 is bit-identical to
20c-1" and "CLI: ``--start-quality-gain`` (refusals: non-finite, above 10),
the new ``--stats`` rows and the record"). The invariant-critical suite is
C++ (``tests/cpp/unit/test_quality_gain.cpp``); this file covers the binding
and the CLI. Not mutation-critical.

R7: ``--start-quality-gain DEG``, default 0; ``-1`` restores the hard rule
and, with it, 20c-1's mesh; refused when non-finite or above 10. R8's split
is on when the gain test is and the constraint feet are (``@architect``'s
"Pins ruled for 20c-2's red step", pin 4: ``--no-constraint-feet`` keeps
every point off the lines). Reported: ``start_quality_points_skipped``
keeps its total, and two new rows, "Tries that would not have improved the
angles" and "Land-cover and outline lines split to improve angles" (the
design's wording, asserted here).

PINNED HERE, where the design leaves it open (listed for ``@architect``):

- the binding: keyword ``min_gain_deg`` on ``_core.refine``, default ``-1.0``
  (off, as ``RefineOptions``), and ``RefineOutcome.quality_no_gain`` and
  ``RefineOutcome.quality_line_splits``;
- the CLI passes ``min_gain_deg`` on every refined run: ``0.0`` by default,
  the given value otherwise, any finite value up to 10 accepted (a negative
  one is the hard rule);
- ``--start-quality-gain`` is refused without ``--tolerance`` and without
  ``--dem``, as ``--start-min-angle`` is;
- the record and ``--stats`` rows: ``start_quality_gain_deg`` (an input,
  beside ``start_min_angle_deg``, printed as ``0``; R7's "``elevation_source``
  adds ``start quality gain 0 deg``" names a field increment 25 replaced by
  such rows), ``start_quality_points_without_gain`` and
  ``start_quality_lines_split``, integers in the record.

T-P3 at the CLI: ``TOPOLOGY_20C1`` was RECORDED FROM 41bda81a (master with
20c-1 merged, no 20c-2 production change) by a scratch program that ran this
file's ``topology`` over the ``RefineOutcome`` ``rasputin mesh`` computes, at
``--tolerance 1`` on the box over ``bumpy`` (alone, and with ``GALLERY``
from ``tests/python/test_cli_mesh_features.py`` as ``--features``) and on
the quarter circle over the Kartverket tile at 1 and 10, each with the
CLI's defaults. It
hashes integers only (triangles, valid, edges, masks, the vertex count and
the counters), so a foot's last bit, which a fused multiply-add can move,
does not. No commit may update it to agree with new code.

Every new binding symbol is read inside a test, so a missing one fails its own
test and leaves the rest of the session collecting.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import typer

import tin_engine.cli as cli
from cli_driver import SQUARE, USAGE, geojson, invoke, rough_dem
from feature_fixtures import write_geojson
from geotiff_fixtures import KARTVERKET, needs_codecs
from recordread import stats_row
from test_cli_mesh_domain import quarter_circle
from test_cli_mesh_features import GALLERY
from test_cli_mesh_stats import mesh
from tin_engine import _core
from vtkread import read_vtk

GAIN = "start_quality_gain_deg"
ROWS = ("start_quality_points_without_gain", "start_quality_lines_split")
WORDING = {
    "start_quality_points_without_gain": "Tries that would not have improved the angles",
    "start_quality_lines_split": "Land-cover and outline lines split to improve angles",
}

#: Recorded from 41bda81a; see the module docstring.
TOPOLOGY_20C1 = {
    "box": "f381fb6ae1d636edbcf285acadbf2607f96b102d938142b01c0a8fd2ebbf8db9",
    "gallery": "3364ca6c823a7ce9eac8c7ab843e685b5081f642c23b798160a921f125f595e8",
    "quarter_circle 1": "307ad79a751cd0a50b819660915343873a325b7923eb93a1550af48c4ad54672",
    "quarter_circle 10": "fbd02f6da4ac149e0fadaed340af0fa8337f9b870c0d98a12cdf5a5de47d0aa5",
}

bumpy = rough_dem(20)


@pytest.fixture
def box(tmp_path: Path) -> Path:
    """SQUARE without its hole: its two start triangles have a 20.7° angle,
    so the quality start has work to do."""
    return geojson(tmp_path / "box.geojson", SQUARE)


class Spy:
    """Every ``cli.refine`` call: the arguments, the keywords and the outcome."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.args: list[tuple[Any, ...]] = []
        self.kwargs: list[dict[str, Any]] = []
        self.out: list[Any] = []
        real = cli.refine

        def spy(*args: Any, **kwargs: Any) -> Any:
            self.args.append(args)
            self.kwargs.append(kwargs)
            result = real(*args, **kwargs)
            self.out.append(result)
            return result

        monkeypatch.setattr(cli, "refine", spy)


def run(tmp_path: Path, *args: str) -> tuple[str, dict[str, Any]]:
    """(``--stats`` report, record) of one ``rasputin mesh`` run."""
    out, md, rec = tmp_path / "x.vtk", tmp_path / "x.md", tmp_path / "x.json"
    mesh(*args, "--out", str(out), "--stats", str(md), "--record", str(rec))
    return md.read_text(encoding="utf-8"), json.loads(rec.read_text(encoding="ascii"))


def wording(report: str, name: str) -> str:
    """The label ``--stats`` prints beside the row ``name`` (its first column)."""
    found = re.findall(rf"^\| (.*?) \| .*? \| `?{re.escape(name)}`? \|$", report, re.M)
    assert len(found) == 1, f"{len(found)} rows named {name!r} in\n{report}"
    return str(found[0])


def topology(out: _core.RefineOutcome) -> str:
    """SHA-256 over the integer part of a refine outcome (see the docstring)."""
    h = hashlib.sha256()
    for a in (out.valid, out.triangles, out.edges, out.masks):
        arr = np.ascontiguousarray(a)
        h.update(f"{arr.dtype.str}{arr.shape}".encode())
        h.update(arr.tobytes())
    h.update(
        f"{len(out.vertices)} {out.rounds} {out.inserted} {out.flips} {out.uncovered} "
        f"{out.carved}".encode()
    )
    return h.hexdigest()


def box_args(bumpy: Path, box: Path) -> tuple[str, ...]:
    return ("--dem", str(bumpy), "--domain", str(box), "--tolerance", "1")


# ---------------------------------------------------------------- binding


class TestBinding:
    """``_core.refine`` takes ``min_gain_deg`` and reports R7's and R8's counts."""

    def test_the_keyword_and_the_counts(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        spy = Spy(monkeypatch)
        run(tmp_path, *box_args(bumpy, box))
        (args,), (kwargs,) = spy.args, spy.kwargs
        given = {k: v for k, v in kwargs.items() if k != "min_gain_deg"}
        default = _core.refine(*args, **given)
        off = _core.refine(*args, **given, min_gain_deg=-1.0)
        on = _core.refine(*args, **given, min_gain_deg=0.0)
        for out in (default, off, on):
            assert out.ok(), out.message
            assert isinstance(out.quality_no_gain, int)
            assert isinstance(out.quality_line_splits, int)
            assert out.quality_skipped >= out.quality_no_gain  # one total, R7's skips in it
        assert (off.quality_no_gain, off.quality_line_splits) == (0, 0)
        assert (default.quality_no_gain, default.quality_line_splits) == (0, 0)
        assert topology(default) == topology(off)  # the default is off (pinned)
        assert np.array_equal(np.asarray(default.vertices), np.asarray(off.vertices))


# ---------------------------------------------------------------- the flag


class TestTheFlag:
    """R7: default 0, forwarded on both start paths, ``-1`` the hard rule."""

    def test_the_default_is_0_on_a_domain_start(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        spy = Spy(monkeypatch)
        report, record = run(tmp_path, *box_args(bumpy, box))
        assert [k["min_gain_deg"] for k in spy.kwargs] == [0.0]
        assert stats_row(report, GAIN) == "0"
        assert record[GAIN] == 0

    def test_the_default_is_0_on_a_stride_start(
        self, tmp_path: Path, bumpy: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        spy = Spy(monkeypatch)
        report, _ = run(tmp_path, "--dem", str(bumpy), "--tolerance", "1")
        assert [k["min_gain_deg"] for k in spy.kwargs] == [0.0]
        assert stats_row(report, GAIN) == "0"

    @pytest.mark.parametrize("value", ["-1", "-0.5", "0", "2.5", "10"])
    def test_a_finite_value_up_to_10_is_passed_and_named(
        self,
        tmp_path: Path,
        bumpy: Path,
        box: Path,
        value: str,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        spy = Spy(monkeypatch)
        report, record = run(tmp_path, *box_args(bumpy, box), "--start-quality-gain", value)
        assert [k["min_gain_deg"] for k in spy.kwargs] == [float(value)]
        assert stats_row(report, GAIN) == value
        assert record[GAIN] == float(value)

    def test_minus_1_turns_both_rules_off(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        spy = Spy(monkeypatch)
        report, record = run(tmp_path, *box_args(bumpy, box), "--start-quality-gain", "-1")
        (out,) = spy.out
        assert (out.quality_no_gain, out.quality_line_splits) == (0, 0)
        for name in ROWS:
            assert stats_row(report, name) == "0", name
            assert record[name] == 0, name

    def test_the_help_names_the_flag(self) -> None:
        code, output = invoke("mesh", "--help")
        assert code == 0, output
        assert "--start-quality-gain" in output

    def test_the_help_ties_the_line_split_to_the_constraint_feet(self) -> None:
        """Pin 4: the line split (R8) runs only with the gain test and the
        constraint feet both on, so the clause of the help that names the split
        names the feet too. The meaning is pinned, not the wording: the clause
        is the ``;``- or ``.``-delimited piece that says "split". Read from the
        option itself, so terminal width cannot wrap the words apart."""
        (option,) = [
            p
            for p in typer.main.get_command(cli.app).commands["mesh"].params
            if "--start-quality-gain" in p.opts
        ]
        clauses = [c for c in re.split(r"[;.]\s", option.help or "") if "split" in c]
        assert clauses, option.help
        for clause in clauses:
            assert re.search(r"constraint[- ]feet", clause), clause


class TestRefusals:
    """R7: non-finite and above 10 are refused, for their value; and, as
    ``--start-min-angle``, without ``--tolerance`` or ``--dem`` (pinned)."""

    @pytest.mark.parametrize("value", ["nan", "inf", "-inf", "10.5", "11"])
    def test_a_bad_value(self, tmp_path: Path, bumpy: Path, box: Path, value: str) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh", *box_args(bumpy, box), "--start-quality-gain", value, "--out", str(target)
        )
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its value, not unknown
        assert "--start-quality-gain" in output
        assert not target.exists()

    def test_without_tolerance(self, tmp_path: Path, bumpy: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh", "--dem", str(bumpy), "--start-quality-gain", "0", "--out", str(target)
        )
        assert code == USAGE, output
        assert "No such option" not in output
        assert "--start-quality-gain" in output
        assert "--tolerance" in output
        assert not target.exists()

    def test_without_dem(self, tmp_path: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh", "catchment", "--flat", "--start-quality-gain", "0", "--out", str(target)
        )
        assert code == USAGE, output
        assert "No such option" not in output
        assert "--start-quality-gain" in output
        assert not target.exists()


# ---------------------------------------------------------------- the report


class TestReport:
    """R7 and R8: the two rows carry the outcome's counts, in the design's
    words; the skipped total is unchanged in meaning."""

    def test_the_rows_match_the_outcome(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        spy = Spy(monkeypatch)
        report, record = run(tmp_path, *box_args(bumpy, box))
        (out,) = spy.out
        assert stats_row(report, ROWS[0]) == str(out.quality_no_gain)
        assert stats_row(report, ROWS[1]) == str(out.quality_line_splits)
        assert stats_row(report, "start_quality_points_skipped") == str(out.quality_skipped)
        for name in ROWS:
            assert type(record[name]) is int, name
            assert str(record[name]) == stats_row(report, name)

    def test_the_rows_say_it_in_the_designs_words(
        self, tmp_path: Path, bumpy: Path, box: Path
    ) -> None:
        report, _ = run(tmp_path, *box_args(bumpy, box))
        for name, words in WORDING.items():
            assert wording(report, name) == words, name

    def test_the_setting_is_not_in_the_file(self, tmp_path: Path, bumpy: Path, box: Path) -> None:
        """D2 of increment 25: the mesh file carries what a user of the mesh
        needs; the setting is ``--stats``' and the record's."""
        out = tmp_path / "x.vtk"
        mesh(*box_args(bumpy, box), "--out", str(out))
        vtk = read_vtk(out.read_bytes())
        for name in (GAIN, *ROWS):
            assert name not in vtk.field_data, name


# ---------------------------------------------------------------- T-P3 at the CLI


def _box_outcome(
    tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch, *extra: str
) -> _core.RefineOutcome:
    """The ``RefineOutcome`` ``rasputin mesh`` computes on the box over ``bumpy``."""
    spy = Spy(monkeypatch)
    code, output = invoke("mesh", *box_args(bumpy, box), "--out", str(tmp_path / "x.vtk"), *extra)
    assert code == 0, output
    (out,) = spy.out
    return out


def _quarter_outcome(
    tolerance: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *extra: str
) -> _core.RefineOutcome:
    """The ``RefineOutcome`` ``rasputin mesh`` computes on the quarter circle."""
    spy = Spy(monkeypatch)
    domain = geojson(tmp_path / "quarter.geojson", quarter_circle())
    args = ["--dem", str(KARTVERKET), "--domain", str(domain), "--tolerance", tolerance]
    code, output = invoke("mesh", *args, "--out", str(tmp_path / "x.vtk"), *extra)
    assert code == 0, output
    (out,) = spy.out
    return out


@pytest.mark.parametrize("scene", ["box", "gallery"])
def test_minus_1_is_20c1s_mesh_on_the_box(
    scene: str, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The box over ``bumpy``, alone and with test_cli_mesh_features' gallery."""
    extra = ["--start-quality-gain", "-1"]
    if scene == "gallery":
        extra += ["--features", str(write_geojson(tmp_path / "gallery.geojson", GALLERY))]
    out = _box_outcome(tmp_path, bumpy, box, monkeypatch, *extra)
    assert (out.quality_no_gain, out.quality_line_splits) == (0, 0)
    assert topology(out) == TOPOLOGY_20C1[scene]


@needs_codecs
@pytest.mark.parametrize("tolerance", ["1", "10"])
def test_minus_1_is_20c1s_mesh_on_the_quarter_circle(
    tolerance: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _quarter_outcome(tolerance, tmp_path, monkeypatch, "--start-quality-gain", "-1")
    assert (out.quality_no_gain, out.quality_line_splits) == (0, 0)
    assert topology(out) == TOPOLOGY_20C1[f"quarter_circle {tolerance}"]


def test_the_default_changes_the_box(
    tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The converse, so the tests above cannot pass with a rule that never
    runs: the box's quality start has candidates R7 refuses."""
    out = _box_outcome(tmp_path, bumpy, box, monkeypatch)
    assert out.quality_no_gain > 0
    assert topology(out) != TOPOLOGY_20C1["box"]
