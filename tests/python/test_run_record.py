"""`tin_engine.run_record`: increment 25's record of one run, without the extension.

`docs/increments/25-plain-output.md`, D1 to D3, D5 and "Tests for @tester"
("The record alone", "The at-vertex rule"). Not invariant-critical, so no
mutation round (README, "Cost constraints").

The module is fetched inside a fixture, so its absence fails these tests and
leaves the rest of the session collecting.

Pinned by the design: ``Entry(name, wording, value, number, in_file)`` and
``RunRecord(entries)``, frozen; ``file_fields(record)`` as ``(name, value)``
pairs; ``stats_rows(record)`` as ``(wording, value, name)``;
``summary(record)``; ``as_json(record, version, command)``; the builders
``refined_record``, ``stride_record`` and ``flat_record``; D2's file fields
and their wording ("The fields, for Ola"); D3's rules; D5's JSON.

PINNED HERE, where D1 is silent (the builders' arguments): every builder takes
keyword arguments only, named as the record's entries where an argument is
one entry's value, plus ``triangles`` (the output's triangle count, which the
summary names, D7). Every argument other than the ones each helper below
passes defaults to "not produced on this path", so the entry is absent (D5:
"absent, not null"). ``refined_record``'s ``max_error_m`` is the measured
figure (refinement's, or the final check's on the reprojected path) before
D2's at-vertex rule raises it. ``snap_to_lines``, the start mesh and the
input descriptions are not passed here; the CLI suites check them. The
warnings of D2 and D3 rule 4 are part of ``summary``'s text, each starting
``Warning:`` (D7's example), since ``summary`` is the one stderr text the
module returns.
"""

from __future__ import annotations

import importlib
import json
import math
import re
from types import ModuleType
from typing import Any

import pytest

#: D3 rule 6, matched case-blind as whole words (``snap`` also in
#: ``snapped``); ``snap_to_lines`` is a name, and names are not searched.
BANNED = re.compile(
    r"(?i)\buncovered\b|\bfeet\b|\bstride\b|\bchains\b|\bcoincident\b|\bcarved\b|"
    r"\bnoded\b|R-tree|\bsnap|\bDEM holes\b|\bvoid\b"
)
NAME = re.compile(r"[a-z][a-z0-9]*(?:_[a-z0-9]+)*")
#: A number a person reads: not part of a word, a file name or ``EPSG:31983``.
NUMBER_TOKEN = re.compile(r"(?<![\w.:])-?\d+(?:\.\d+)?(?:e[-+]?\d+)?(?![\w.:])")

#: D2's table: every mesh-file field and its plain wording ("The fields, for Ola").
FILE_WORDING = {
    "crs": "Coordinate system",
    "tolerance_m": "Tolerance",
    "max_error_m": "Largest height error",
    "dem_source": "DEM",
    "dem_credit": "Credit",
    "licence_note": "Licence",
    "cite": "Please cite",
    "nodata_vertices_removed": "Vertices removed on NoData",
    "heights": "Heights",
}
FLAT_HEIGHTS = "none: every z is 0 (--flat)"
CREDIT = "Agencia Nacional de Aguas (ANA), Brazil; https://doi.org/10.5069/G9736P4G"
LICENCE = "CC BY 4.0 (Creative Commons Attribution 4.0 International)"
CITE = "Laipelt, L., et al. (2024). ANADEM"


@pytest.fixture
def rr() -> ModuleType:
    return importlib.import_module("tin_engine.run_record")


# ---------------------------------------------------------------- builders


def projected(rr: ModuleType, **over: Any) -> Any:
    """A refined run on a projected DEM: the inventory's first run."""
    args: dict[str, Any] = {
        "crs": "EPSG:25833",
        "tolerance_m": 5.0,
        "max_error_m": 4.999692612800061,
        "triangles": 116389,
        "nodata_vertices_removed": 199,
        "dem_nodes_outside_mesh": 0,
        "dem_source": "7908_3_10m_z33.tif",
    }
    return rr.refined_record(**(args | over))


def reprojected(rr: ModuleType, **over: Any) -> Any:
    """The inventory's Velhas run: the final check's figure, phase 1's apart."""
    args: dict[str, Any] = {
        "crs": "EPSG:31983",
        "tolerance_m": 5.0,
        "max_error_m": 4.992618429102549,
        "triangles": 983,
        "nodata_vertices_removed": 0,
        "dem_nodes_outside_mesh": 0,
        "dem_source": "anadem_velhas.tif",
        "resampled_grid_max_error_m": 4.992000034877265,
        "dem_nodes_checked": 21615,
        "dem_nodes_at_vertices": 0,
        "dem_nodes_at_vertices_max_error_m": 0.0,
    }
    return rr.refined_record(**(args | over))


def downloaded(rr: ModuleType, **over: Any) -> Any:
    notes = {"dem_source": "anadem-v1", "dem_credit": CREDIT, "licence_note": LICENCE, "cite": CITE}
    return reprojected(rr, **(notes | over))


def stride(rr: ModuleType, **over: Any) -> Any:
    args: dict[str, Any] = {
        "crs": "EPSG:25833",
        "triangles": 8,
        "nodata_vertices_removed": 6,
        "dem_source": "tile.tif",
    }
    return rr.stride_record(**(args | over))


def flat(rr: ModuleType, **over: Any) -> Any:
    args: dict[str, Any] = {"crs": "EPSG:25833", "triangles": 60}
    return rr.flat_record(**(args | over))


PATHS = ("projected", "reprojected", "downloaded", "stride", "flat")


def record_of(rr: ModuleType, path: str) -> Any:
    return {
        "projected": projected,
        "reprojected": reprojected,
        "downloaded": downloaded,
        "stride": stride,
        "flat": flat,
    }[path](rr)


def fields(rr: ModuleType, record: Any) -> dict[str, str]:
    pairs = list(rr.file_fields(record))
    names = [name for name, _ in pairs]
    assert len(names) == len(set(names)), names
    return dict(pairs)


def rows(rr: ModuleType, record: Any) -> dict[str, tuple[str, str]]:
    """``stats_rows`` by name: (wording, value)."""
    found = list(rr.stats_rows(record))
    names = [name for _, _, name in found]
    assert len(names) == len(set(names)), names
    return {name: (wording, value) for wording, value, name in found}


def numbers(record: Any) -> dict[str, float | int]:
    return {e.name: e.number for e in record.entries if e.number is not None}


# ---------------------------------------------------------------- the types


class TestTheTypes:
    def test_entry_and_record_are_frozen_value_types(self, rr: ModuleType) -> None:
        entry = rr.Entry(
            name="tolerance_m", wording="Tolerance", value="5", number=5.0, in_file=True
        )
        record = rr.RunRecord(entries=(entry,))
        assert record.entries == (entry,)
        with pytest.raises(AttributeError):
            entry.value = "6"  # type: ignore[misc]
        assert rr.file_fields(record) == [("tolerance_m", "5")]
        assert rr.stats_rows(record) == [("Tolerance", "5", "tolerance_m")]

    def test_an_entry_not_in_the_file_is_only_in_stats(self, rr: ModuleType) -> None:
        entry = rr.Entry(
            name="edge_flips", wording="Edge swaps", value="7", number=7, in_file=False
        )
        record = rr.RunRecord(entries=(entry,))
        assert rr.file_fields(record) == []
        assert rr.stats_rows(record) == [("Edge swaps", "7", "edge_flips")]


# ---------------------------------------------------------------- names and values


@pytest.mark.parametrize("path", PATHS)
class TestNamesAndValues:
    def test_names_are_plain_snake_case_and_unique(self, rr: ModuleType, path: str) -> None:
        names = [e.name for e in record_of(rr, path).entries]
        assert len(names) == len(set(names)), names
        for name in names:
            assert NAME.fullmatch(name), name

    def test_a_measured_value_has_its_unit_in_its_name(self, rr: ModuleType, path: str) -> None:
        """D3 rule 2, and D5: a float is a measured value ending ``_m`` or
        ``_deg``; a count is an int and carries no unit."""
        for e in record_of(rr, path).entries:
            if isinstance(e.number, float):
                assert e.name.endswith(("_m", "_deg")), e.name
            elif e.number is not None:
                assert type(e.number) is int, (e.name, e.number)
                assert not e.name.endswith(("_m", "_deg")), e.name
            if e.name.endswith(("_m", "_deg")):
                assert isinstance(e.number, float), (e.name, e.number)

    def test_values_are_ascii_without_control_characters(self, rr: ModuleType, path: str) -> None:
        for e in record_of(rr, path).entries:
            assert e.value.isascii(), (e.name, e.value)
            assert not any(ch < " " or ch == "\x7f" for ch in e.value), (e.name, e.value)

    def test_a_number_prints_as_its_value(self, rr: ModuleType, path: str) -> None:
        """D5: the float ``--record`` writes equals the one the ``--stats``
        text parses to, so the two agree to the bit."""
        for e in record_of(rr, path).entries:
            if e.number is not None:
                assert float(e.value) == e.number, (e.name, e.value, e.number)

    def test_no_banned_word_in_any_value_wording_or_the_summary(
        self, rr: ModuleType, path: str
    ) -> None:
        record = record_of(rr, path)
        for e in record.entries:
            assert not BANNED.search(e.value), (e.name, e.value)
            assert not BANNED.search(e.wording), (e.name, e.wording)
        assert not BANNED.search(rr.summary(record)), rr.summary(record)

    def test_nodata_is_called_nodata(self, rr: ModuleType, path: str) -> None:
        """D3 rule 7."""
        record = record_of(rr, path)
        texts = [e.value for e in record.entries] + [e.wording for e in record.entries]
        for text in [*texts, rr.summary(record)]:
            assert not re.search(r"(?i)\bholes?\b.*\bDEM\b|\bDEM\b.*\bholes?\b", text), text


# ---------------------------------------------------------------- the file fields (D2)


class TestFileFields:
    @pytest.mark.parametrize(
        ("path", "expected"),
        [
            (
                "projected",
                {"crs", "tolerance_m", "max_error_m", "dem_source", "nodata_vertices_removed"},
            ),
            ("reprojected", {"crs", "tolerance_m", "max_error_m", "dem_source"}),
            (
                "downloaded",
                {
                    "crs",
                    "tolerance_m",
                    "max_error_m",
                    "dem_source",
                    "dem_credit",
                    "licence_note",
                    "cite",
                },
            ),
            ("stride", {"crs", "dem_source", "nodata_vertices_removed"}),
            ("flat", {"crs", "heights"}),
        ],
    )
    def test_exactly_d2s_fields_for_each_path(
        self, rr: ModuleType, path: str, expected: set[str]
    ) -> None:
        assert set(fields(rr, record_of(rr, path))) == expected

    @pytest.mark.parametrize("path", PATHS)
    def test_their_wording_is_the_tables(self, rr: ModuleType, path: str) -> None:
        record = record_of(rr, path)
        by_name = rows(rr, record)
        for name in fields(rr, record):
            assert by_name[name][0] == FILE_WORDING[name], name

    def test_the_values(self, rr: ModuleType) -> None:
        got = fields(rr, downloaded(rr))
        assert got["crs"] == "EPSG:31983"
        assert got["tolerance_m"] == "5"
        assert got["max_error_m"] == "4.992618429102549"
        assert got["dem_source"] == "anadem-v1"
        assert (got["dem_credit"], got["licence_note"], got["cite"]) == (CREDIT, LICENCE, CITE)
        assert fields(rr, projected(rr))["nodata_vertices_removed"] == "199"

    def test_flat_says_the_heights_are_not_real(self, rr: ModuleType) -> None:
        assert fields(rr, flat(rr))["heights"] == FLAT_HEIGHTS

    def test_flat_without_a_crs_has_no_crs_field(self, rr: ModuleType) -> None:
        assert fields(rr, flat(rr, crs=None)) == {"heights": FLAT_HEIGHTS}

    def test_a_source_that_asks_no_citation_writes_no_cite(self, rr: ModuleType) -> None:
        assert "cite" not in fields(rr, downloaded(rr, cite=None))


# ---------------------------------------------------------------- the omission rules (D3)


class TestOmission:
    @pytest.mark.parametrize("path", ["projected", "stride"])
    def test_a_zero_count_is_absent_from_the_file_and_present_in_stats(
        self, rr: ModuleType, path: str
    ) -> None:
        record = {"projected": projected, "stride": stride}[path](rr, nodata_vertices_removed=0)
        assert "nodata_vertices_removed" not in fields(rr, record)
        assert rows(rr, record)["nodata_vertices_removed"][1] == "0"

    def test_a_measured_zero_is_written(self, rr: ModuleType) -> None:
        """Rule 3 omits counts only: ``--tolerance 0`` writes both measured values."""
        got = fields(rr, projected(rr, tolerance_m=0.0, max_error_m=0.0))
        assert (got["tolerance_m"], got["max_error_m"]) == ("0", "0")

    def test_every_file_field_is_in_stats_with_the_same_value(self, rr: ModuleType) -> None:
        """D3 rule 1."""
        for path in PATHS:
            record = record_of(rr, path)
            by_name = rows(rr, record)
            for name, value in fields(rr, record).items():
                assert by_name[name][1] == value, (path, name)

    def test_the_at_vertex_entries_are_absent_where_not_measured(self, rr: ModuleType) -> None:
        """D6, the default taken while Ola was away: a projected run (before
        15f-3) has neither entry, rather than a 0 nobody measured; a path that
        measures them has both, 0 included."""
        for record in (projected(rr), stride(rr), flat(rr)):
            names = {e.name for e in record.entries}
            assert "dem_nodes_at_vertices" not in names
            assert "dem_nodes_at_vertices_max_error_m" not in names
        by_name = rows(rr, reprojected(rr))
        assert by_name["dem_nodes_at_vertices"][1] == "0"
        assert by_name["dem_nodes_at_vertices_max_error_m"][1] == "0"

    def test_no_tolerance_means_no_tolerance_fields(self, rr: ModuleType) -> None:
        for record in (stride(rr), flat(rr)):
            names = {e.name for e in record.entries}
            assert not names & {"tolerance_m", "max_error_m", "refinement_rounds"}, names

    def test_the_reprojected_figures_are_stats_rows(self, rr: ModuleType) -> None:
        by_name = rows(rr, reprojected(rr))
        assert by_name["resampled_grid_max_error_m"][1] == "4.992000034877265"
        assert by_name["dem_nodes_checked"][1] == "21615"
        assert "resampled_grid_max_error_m" not in fields(rr, reprojected(rr))


# ---------------------------------------------------------------- max_error_m (D2, D3 rule 2)


class TestMaxError:
    def test_printed_never_above_the_tolerance_it_met(self, rr: ModuleType) -> None:
        just_under = math.nextafter(1.0, 0.0)
        got = fields(rr, projected(rr, tolerance_m=1.0, max_error_m=just_under))
        assert float(got["max_error_m"]) <= float(got["tolerance_m"])
        assert got["max_error_m"] == repr(just_under)

    def test_a_short_value_prints_short(self, rr: ModuleType) -> None:
        got = fields(rr, projected(rr, tolerance_m=5.0, max_error_m=4.5))
        assert (got["tolerance_m"], got["max_error_m"]) == ("5", "4.5")

    def test_the_at_vertex_difference_raises_it(self, rr: ModuleType) -> None:
        """D2: max_error_m = max(the final check's figure, the at-vertex one)."""
        record = reprojected(
            rr, max_error_m=0.9, dem_nodes_at_vertices=3, dem_nodes_at_vertices_max_error_m=0.95
        )
        assert fields(rr, record)["max_error_m"] == "0.95"
        assert numbers(record)["max_error_m"] == 0.95
        assert rows(rr, record)["dem_nodes_at_vertices_max_error_m"][1] == "0.95"
        assert "Warning" not in rr.summary(record)  # 0.95 is within the tolerance of 5

    def test_a_smaller_at_vertex_difference_leaves_it(self, rr: ModuleType) -> None:
        record = reprojected(
            rr, max_error_m=4.9, dem_nodes_at_vertices=2, dem_nodes_at_vertices_max_error_m=0.25
        )
        assert fields(rr, record)["max_error_m"] == "4.9"

    def test_above_the_tolerance_it_reads_above_and_warns(self, rr: ModuleType) -> None:
        """D2: the file does not claim a tolerance the mesh misses, and
        stderr names the count and the difference."""
        record = reprojected(
            rr, max_error_m=4.9, dem_nodes_at_vertices=3, dem_nodes_at_vertices_max_error_m=5.5
        )
        got = fields(rr, record)
        assert got["max_error_m"] == "5.5"
        assert float(got["max_error_m"]) > float(got["tolerance_m"])
        text = rr.summary(record)
        warnings = [line for line in text.splitlines() if line.startswith("Warning:")]
        assert len(warnings) == 1, text
        assert re.search(r"\b3\b", warnings[0]) and "5.5" in warnings[0], warnings[0]

    def test_one_node_above_the_tolerance_is_singular(self, rr: ModuleType) -> None:
        """Plurals are correct English ("Settled after the red step", 8)."""
        record = reprojected(
            rr, max_error_m=4.9, dem_nodes_at_vertices=1, dem_nodes_at_vertices_max_error_m=5.5
        )
        (warning,) = [ln for ln in rr.summary(record).splitlines() if ln.startswith("Warning:")]
        assert warning.startswith("Warning: 1 DEM node on a vertex differs"), warning
        assert "nodes" not in warning, warning

    def test_the_stats_value_is_the_raw_at_vertex_difference(self, rr: ModuleType) -> None:
        record = reprojected(
            rr, max_error_m=4.9, dem_nodes_at_vertices=3, dem_nodes_at_vertices_max_error_m=5.5
        )
        by_name = rows(rr, record)
        assert by_name["dem_nodes_at_vertices"][1] == "3"
        assert by_name["dem_nodes_at_vertices_max_error_m"][1] == "5.5"


# ---------------------------------------------------------------- self-checks (D3 rule 4)


class TestSelfChecks:
    def test_zero_is_silent(self, rr: ModuleType) -> None:
        record = projected(rr)
        assert rows(rr, record)["dem_nodes_outside_mesh"][1] == "0"
        assert "Warning" not in rr.summary(record)
        assert "dem_nodes_outside_mesh" not in fields(rr, record)

    def test_one_node_outside_is_singular(self, rr: ModuleType) -> None:
        text = rr.summary(projected(rr, dem_nodes_outside_mesh=1))
        (warning,) = [ln for ln in text.splitlines() if ln.startswith("Warning:")]
        assert warning.startswith("Warning: 1 DEM node with data lies"), warning
        assert "nodes" not in warning, warning

    def test_a_non_zero_self_check_warns(self, rr: ModuleType) -> None:
        record = projected(rr, dem_nodes_outside_mesh=4)
        assert rows(rr, record)["dem_nodes_outside_mesh"][1] == "4"
        text = rr.summary(record)
        warnings = [line for line in text.splitlines() if line.startswith("Warning:")]
        assert len(warnings) == 1 and re.search(r"\b4\b", warnings[0]), text
        assert not BANNED.search(text), text


# ---------------------------------------------------------------- the summary (D7)


class TestSummary:
    def test_the_projected_summary(self, rr: ModuleType) -> None:
        text = rr.summary(projected(rr))
        assert text.startswith("116389 triangles."), text
        assert "within 5 m" in text and "largest difference 4.9997 m" in text, text
        assert "199 vertices on NoData cells" in text, text

    def test_the_no_tolerance_summary_says_on_or_next_to(self, rr: ModuleType) -> None:
        """Increment 12's one-cell trim, until the sampler fix ("NoData on the
        no-tolerance path"): the stride path's count includes vertices next to
        a NoData cell, and the summary says so; the tolerance path does not."""
        assert "6 vertices on or next to NoData cells were removed" in rr.summary(stride(rr))
        assert "on or next to" not in rr.summary(projected(rr))

    def test_no_nodata_sentence_without_nodata(self, rr: ModuleType) -> None:
        assert "NoData" not in rr.summary(projected(rr, nodata_vertices_removed=0))

    def test_the_reprojected_summary_names_the_original_dem(self, rr: ModuleType) -> None:
        text = rr.summary(reprojected(rr))
        assert text.startswith("983 triangles."), text
        assert "original DEM" in text and "largest difference 4.9926 m" in text, text

    @pytest.mark.parametrize("path", PATHS)
    def test_every_number_is_a_field_to_5_significant_figures(
        self, rr: ModuleType, path: str
    ) -> None:
        record = record_of(rr, path)
        known = {float(f"{v:.5g}") for v in numbers(record).values()}
        known |= {float(record_of_triangles(path)), 0.0}  # "every z is 0" on --flat
        text = rr.summary(record)
        for token in NUMBER_TOKEN.findall(text):
            assert float(token) in known, (token, text)

    @pytest.mark.parametrize("path", PATHS)
    def test_it_is_ascii_and_one_line_without_warnings(self, rr: ModuleType, path: str) -> None:
        text = rr.summary(record_of(rr, path))
        assert text.isascii()
        assert "\n" not in text.rstrip("\n"), text


def record_of_triangles(path: str) -> int:
    return {"projected": 116389, "reprojected": 983, "downloaded": 983, "stride": 8, "flat": 60}[
        path
    ]


# ---------------------------------------------------------------- --record's JSON (D5)


class TestAsJson:
    VERSION = "0.2.0.dev0"
    COMMAND = "rasputin mesh --dem a.tif --tolerance 5 --out a.vtk --record a.json"

    def text(self, rr: ModuleType, record: Any) -> str:
        out = rr.as_json(record, self.VERSION, self.COMMAND)
        assert isinstance(out, str)
        return out

    @pytest.mark.parametrize("path", PATHS)
    def test_the_shape(self, rr: ModuleType, path: str) -> None:
        record = record_of(rr, path)
        text = self.text(rr, record)
        assert text.isascii() and text.endswith("\n") and not text.endswith("\n\n")
        obj = json.loads(text)
        names = [name for _, _, name in rr.stats_rows(record)]
        assert list(obj) == ["rasputin_version", "command", *names]
        assert obj["rasputin_version"] == self.VERSION and obj["command"] == self.COMMAND
        # D5's exact formatting: indent 1, ASCII, one trailing newline.
        assert text == json.dumps(obj, indent=1, ensure_ascii=True) + "\n"

    @pytest.mark.parametrize("path", PATHS)
    def test_the_values_follow_stats(self, rr: ModuleType, path: str) -> None:
        """A count is a JSON integer, a measured value the float its ``--stats``
        text parses to, a text value that text."""
        record = record_of(rr, path)
        obj = json.loads(self.text(rr, record))
        by_name = {e.name: e for e in record.entries}
        for _, value, name in rr.stats_rows(record):
            got, number = obj[name], by_name[name].number
            if isinstance(number, float):
                assert type(got) is float and got == float(value), (name, got, value)
            elif number is not None:
                assert type(got) is int and got == int(value), (name, got, value)
            else:
                assert got == value, (name, got, value)

    def test_a_whole_measured_value_keeps_its_point(self, rr: ModuleType) -> None:
        text = self.text(rr, projected(rr))
        assert '"tolerance_m": 5.0' in text, text
        assert '"nodata_vertices_removed": 199' in text, text

    def test_an_omitted_entry_is_absent_not_null(self, rr: ModuleType) -> None:
        obj = json.loads(self.text(rr, projected(rr)))
        assert None not in obj.values()
        assert "dem_nodes_at_vertices" not in obj and "dem_credit" not in obj

    def test_every_file_field_is_in_it_with_the_same_value(self, rr: ModuleType) -> None:
        for path in PATHS:
            record = record_of(rr, path)
            obj = json.loads(self.text(rr, record))
            for name, value in fields(rr, record).items():
                got = obj[name]
                assert (got if isinstance(got, str) else float(got)) == (
                    value if isinstance(got, str) else float(value)
                ), (path, name, got, value)

    def test_the_same_record_gives_the_same_bytes(self, rr: ModuleType) -> None:
        assert self.text(rr, downloaded(rr)) == self.text(rr, downloaded(rr))
