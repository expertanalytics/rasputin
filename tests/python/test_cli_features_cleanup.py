"""The land-cover clean-up's four flags through ``rasputin mesh`` (20c-3).

``docs/increments/20c-soft-quality.md``, "Design of PR 20c-3" ("Steps",
"Off", "Recorded") and "Tests ``@tester`` writes red first", 20c-3: the CLI
item under the repair ("``--features-repair`` defaults to 0.05 whenever a
land-cover map is used (question 9's default); refused: negative,
non-finite; the input row ``features_repair_m`` and the vertex rows in
``--stats`` and the record; the ``features clean-up`` timing sub-row present
under ``features clip``") and the outline rule's CLI test (5 m whenever
features are given, ruling 6; ``--features-outline-snap 0`` turns it off).
Not mutation-critical. Ruling G3 ("Rulings on 20c-3's green step"): the
merge switch is the text ``on`` or ``off``, the same in ``--stats`` and the
record, as ``snap_to_lines`` is (increment 25, "Values").

**If Ola picks another value for question 9**, only ``REPAIR_DEFAULT`` and
the row text (``"1"``) here change.

PINNED HERE, where the design leaves it open (listed for ``@architect``):

- The CLI hands the four switches to ``open_features`` in its
  ``FeatureRequest`` (``repair_m``, ``merge_same_class``, ``tolerance_m``,
  ``outline_snap_m``; see ``test_feature_repair.py``), on every run with
  ``--features``: defaults 1, on, 0, 5, whatever the map.
- ``--stats`` prints them as ``1``, ``0`` and ``5`` (as
  ``start_quality_gain_deg``); the record holds floats.
- No features, none of the four rows.
- The vertex rows are found by the design's wording, "Land-cover vertices
  before and after clean-up"; their name and value format are not pinned.
- The timing sub-row is the phase ``features clip: clean-up`` (the
  ``"<parent>: ..."`` form ``stats._timings`` lists as a sub-row).
- Refused, for its value: ``--features-repair`` and
  ``--features-tolerance`` negative or non-finite;
  ``--features-outline-snap`` negative, non-finite, or at or above 100.

RED at the commit that adds ``test_the_area_row_says_another_polygon``
(on ``a4d50033``): the row's wording is "... that changed class, m2"
(ruling R4).

RED at the commit that adds this file: none of the four flags exists, so
every run that passes one exits 2 with "No such option" (which every refusal
test excludes), and the default runs hand ``open_features`` a request without
the four fields.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
import typer

import tin_engine.cli as cli
import tin_engine.run_record as rr
from cli_driver import SQUARE, USAGE, geojson, invoke, rough_dem
from feature_fixtures import write_geojson
from recordread import stats_names, stats_row
from test_cli_mesh_features import GALLERY
from test_cli_mesh_landcover import EAST, WEST, coded

REPAIR_DEFAULT = 1.0  # question 9's default, 1 m (Ola, 2026-10-08)
OUTLINE_DEFAULT = 5.0  # ruling 6
ROWS = (
    "features_repair_m",
    "features_merge_same_class",
    "features_tolerance_m",
    "features_outline_snap_m",
)
FLAGS = (
    "--features-repair",
    "--features-merge-same-class",
    "--features-tolerance",
    "--features-outline-snap",
)
VERTEX_WORDING = "Land-cover vertices before and after clean-up"
#: Ruling R4 ("Rulings on 20c-3's code review round 1"): with the merge off
#: the rule can move area between two polygons of one class, so not "changed
#: class".
AREA_WORDING = (
    "Land-cover area inside the outline that the outline rule gave to another polygon, m2"
)

bumpy = rough_dem(16)


@pytest.fixture
def square(tmp_path: Path) -> Path:
    return geojson(tmp_path / "square.geojson", SQUARE)


@pytest.fixture
def corine(tmp_path: Path) -> Path:
    return write_geojson(
        tmp_path / "corine.geojson", [coded(1, WEST, "311"), coded(2, EAST, "512")]
    )


@pytest.fixture
def gallery(tmp_path: Path) -> Path:
    return write_geojson(tmp_path / "gallery.geojson", GALLERY)


class Spy:
    """Every ``cli.open_features`` call's request."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.requests: list[Any] = []
        real = cli.open_features

        def spy(request: Any, *args: Any, **kwargs: Any) -> Any:
            self.requests.append(request)
            return real(request, *args, **kwargs)

        monkeypatch.setattr(cli, "open_features", spy)


def run(tmp_path: Path, bumpy: Path, square: Path, *extra: str) -> tuple[str, dict[str, Any], str]:
    """(``--stats`` report, record, output) of one successful run."""
    out, md, rec = tmp_path / "x.vtk", tmp_path / "x.md", tmp_path / "x.json"
    code, output = invoke(
        "mesh",
        "--dem",
        str(bumpy),
        "--domain",
        str(square),
        "--tolerance",
        "1",
        "--out",
        str(out),
        "--stats",
        str(md),
        "--record",
        str(rec),
        *extra,
    )
    assert code == 0, output
    return md.read_text(encoding="utf-8"), json.loads(rec.read_text(encoding="ascii")), output


def switches(request: Any) -> tuple[float, bool, float, float]:
    return (request.repair_m, request.merge_same_class, request.tolerance_m, request.outline_snap_m)


# ---------------------------------------------------------------- defaults


class TestDefaults:
    def test_a_land_cover_map(
        self,
        tmp_path: Path,
        bumpy: Path,
        square: Path,
        corine: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Question 9's default and ruling 6's, with the merge on (ruling 1 (a))
        and no simplification (ruling 1 (b))."""
        spy = Spy(monkeypatch)
        report, record, _ = run(
            tmp_path, bumpy, square, "--features", str(corine), "--features-map", "corine"
        )
        (request,) = spy.requests
        assert switches(request) == (REPAIR_DEFAULT, True, 0.0, OUTLINE_DEFAULT)
        assert stats_row(report, "features_repair_m") == "1"
        assert stats_row(report, "features_tolerance_m") == "0"
        assert stats_row(report, "features_outline_snap_m") == "5"
        assert record["features_repair_m"] == REPAIR_DEFAULT
        assert record["features_tolerance_m"] == 0.0
        assert record["features_outline_snap_m"] == OUTLINE_DEFAULT
        assert record["features_merge_same_class"] == "on"
        assert VERTEX_WORDING in report

    def test_the_area_row_says_another_polygon(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        """Ruling R4: ``land_cover_area_moved_m2``'s wording, in ``--stats``
        on a land-cover run and in the record's table. Red on ``1e011cec``
        ("... that changed class, m2")."""
        report, _, _ = run(
            tmp_path, bumpy, square, "--features", str(corine), "--features-map", "corine"
        )
        assert "land_cover_area_moved_m2" in stats_names(report)
        assert AREA_WORDING in report
        assert "changed class" not in report
        assert rr.WORDING["land_cover_area_moved_m2"] == AREA_WORDING

    def test_the_outline_rule_whenever_features_are_given(
        self,
        tmp_path: Path,
        bumpy: Path,
        square: Path,
        gallery: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        spy = Spy(monkeypatch)
        report, _, _ = run(tmp_path, bumpy, square, "--features", str(gallery))
        (request,) = spy.requests
        assert request.outline_snap_m == OUTLINE_DEFAULT
        assert stats_row(report, "features_outline_snap_m") == "5"

    def test_no_features_no_rows(self, tmp_path: Path, bumpy: Path, square: Path) -> None:
        report, record, _ = run(tmp_path, bumpy, square)
        for name in ROWS:
            assert name not in stats_names(report), name
            assert name not in record, name


# ---------------------------------------------------------------- given


class TestGiven:
    def test_every_switch_given(
        self,
        tmp_path: Path,
        bumpy: Path,
        square: Path,
        corine: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        spy = Spy(monkeypatch)
        report, record, _ = run(
            tmp_path,
            bumpy,
            square,
            "--features",
            str(corine),
            "--features-map",
            "corine",
            "--features-repair",
            "0.1",
            "--no-features-merge-same-class",
            "--features-tolerance",
            "2",
            "--features-outline-snap",
            "0",
        )
        (request,) = spy.requests
        assert switches(request) == (0.1, False, 2.0, 0.0)
        assert stats_row(report, "features_repair_m") == "0.1"
        assert stats_row(report, "features_tolerance_m") == "2"
        assert stats_row(report, "features_outline_snap_m") == "0"
        assert record["features_merge_same_class"] == "off"

    def test_the_merge_row_says_on_or_off(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        args = ("--features", str(corine), "--features-map", "corine")
        (tmp_path / "on").mkdir()
        (tmp_path / "off").mkdir()
        on, _, _ = run(tmp_path / "on", bumpy, square, *args)
        off, _, _ = run(tmp_path / "off", bumpy, square, *args, "--no-features-merge-same-class")
        assert stats_row(on, "features_merge_same_class") == "on"
        assert stats_row(off, "features_merge_same_class") == "off"

    def test_the_timing_sub_row(
        self, tmp_path: Path, bumpy: Path, square: Path, corine: Path
    ) -> None:
        """``features clean-up`` under ``features clip``, in the timing table
        of the ``--stats`` file ``run`` writes (ruling G1); the rows are found
        by their cells, in order."""
        report, _, _ = run(
            tmp_path, bumpy, square, "--features", str(corine), "--features-map", "corine"
        )
        clip = re.search(r"\|\s*features clip\s*\|", report)
        sub = re.search(r"\|\s*features clip: clean-up\s*\|", report)
        assert clip and sub, report
        assert clip.start() < sub.start()

    def test_the_help_names_the_flags(self) -> None:
        names = {
            opt
            for p in typer.main.get_command(cli.app).commands["mesh"].params
            for opt in (*p.opts, *p.secondary_opts)
        }
        assert {*FLAGS, "--no-features-merge-same-class"} <= names


# ---------------------------------------------------------------- refusals


REFUSED = [
    ("--features-repair", "-0.01"),
    ("--features-repair", "nan"),
    ("--features-repair", "inf"),
    ("--features-tolerance", "-1"),
    ("--features-tolerance", "nan"),
    ("--features-tolerance", "inf"),
    ("--features-outline-snap", "-1"),
    ("--features-outline-snap", "nan"),
    ("--features-outline-snap", "inf"),
    ("--features-outline-snap", "100"),  # the read region's margin (OR7)
    ("--features-outline-snap", "150"),
]


@pytest.mark.parametrize(("flag", "value"), REFUSED)
def test_a_bad_value_is_refused(
    tmp_path: Path, bumpy: Path, square: Path, corine: Path, flag: str, value: str
) -> None:
    target = tmp_path / "x.vtk"
    code, output = invoke(
        "mesh",
        "--dem",
        str(bumpy),
        "--domain",
        str(square),
        "--tolerance",
        "1",
        "--features",
        str(corine),
        "--features-map",
        "corine",
        f"{flag}={value}",
        "--out",
        str(target),
    )
    assert code == USAGE, output
    assert "No such option" not in output  # refused for its value, not unknown
    assert flag in output
    assert not target.exists()
