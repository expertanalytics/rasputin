"""`io/rivers.py`: `RiverSegment`, `read_segments`, the `kind` mapping, exact copies (29, PR 3).

`docs/increments/29-nve-reference-catchments.md`, "Placing the gauge" (the
model and the `kind` rule), step 4 (exact copies), "The data are not clean"
and "The red suites", PR 3's `test_rivers.py`. Hand-written files only.

The 25 `objekttype` values are ELVIS layer 2's, grouped and counted by the
service itself (`outStatistics`, `groupByFieldsForStatistics=objekttype`, run
by @tester on 2026-10-04: exactly 25 groups, the count the design gives; the
blank one is a single space). The design names 18 of them; the other seven
are river spellings (`ElvBekkregulert`, `ElvBekkMitlinje`, `Elvbekk`,
`FiltivElv`, `ElvBekkMidtlije`, `ElvbekkRegulert`, `ElvBekRegulert`).

HOW THIS FILE GOES RED: there is no `tin_engine/io/rivers.py`
(`ModuleNotFoundError` in each test's fixture).
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from pydantic import ValidationError

from nve_fixtures import CRS, collection, point, river, write

LAKE_SPELLINGS = [
    "InnsjøMidtlinje",
    "InnsjoMidtlinje",
    "InnsjøMidtlinjeReg",
    "InnsjøMitlinje",
    "InnsjøMidtloinje",
    "Innsjømidtlinje",
    "InnsjoMidtlin",
    "InnsjøRegulert",
]
RIVER_SPELLINGS = [
    "ElvBekk",
    "ElvBekkMidtlinje",
    "ElvBekkRegulert",
    "ElvelinjeFiktiv",
    "FiktivElv",
    "BreMidtlinje",
    "SK",
    "20.08.2014",
    "ElvBekkregulert",
    "ElvBekkMitlinje",
    "Elvbekk",
    "FiltivElv",
    "ElvBekkMidtlije",
    "ElvbekkRegulert",
    "ElvBekRegulert",
]
NULL, BLANK = None, " "
LINE = [(0.0, 0.0), (100.0, -10.0), (250.0, -30.0)]


@pytest.fixture
def rivers() -> ModuleType:
    import tin_engine.io.rivers as module

    return module


def read(rivers: ModuleType, path: Path, features: Sequence[dict[str, Any]]) -> Any:
    segments, crs = rivers.read_segments(write(path, collection(features)))
    assert crs == CRS
    return segments


def shifted(coords: Sequence[tuple[float, float]], d: float) -> list[tuple[float, float]]:
    return [(x + d, y - d) for x, y in coords]


class TestTheModel:
    def test_the_fields(self, rivers: ModuleType, tmp_path: Path) -> None:
        feature = river(
            8841, LINE, elvid="2-11-1", vassdragsnr="002.DC", elvenavn="Nea", vatnlnr=None
        )
        (s,) = read(rivers, tmp_path / "r.geojson", [feature])
        assert isinstance(s, rivers.RiverSegment)
        assert (s.objectid, s.elvid, s.vassdragsnr, s.name) == (8841, "2-11-1", "002.DC", "Nea")
        assert (s.objekttype, s.kind) == ("ElvBekk", "river")
        assert s.line == tuple(LINE)

    def test_nulls_are_kept_as_none(self, rivers: ModuleType, tmp_path: Path) -> None:
        feature = river(1, LINE, vassdragsnr=None, elvenavn=None, objekttype=None)
        (s,) = read(rivers, tmp_path / "r.geojson", [feature])
        assert (s.vassdragsnr, s.name, s.objekttype) == (None, None, None)

    def test_a_segment_is_frozen(self, rivers: ModuleType, tmp_path: Path) -> None:
        (s,) = read(rivers, tmp_path / "r.geojson", [river(1, LINE)])
        with pytest.raises(ValidationError, match="frozen"):
            s.kind = "lake"

    def test_the_files_own_crs_is_returned(self, rivers: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path / "r.geojson", collection([river(1, LINE)], crs="EPSG:32633"))
        assert rivers.read_segments(path)[1] == "EPSG:32633"


class TestTheRefusals:
    def test_a_file_without_crs(self, rivers: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path / "r.geojson", collection([river(1, LINE)], crs=None))
        with pytest.raises(ValueError, match=r"(?i)\bcrs\b"):
            rivers.read_segments(path)

    def test_a_duplicate_objectid(self, rivers: ModuleType, tmp_path: Path) -> None:
        """Two features under one `objectid` are a broken file, not copies:
        the fetch keeps a segment seen twice once, by `objectid`."""
        doc = collection(
            [river(41, LINE), river(42, shifted(LINE, 500)), river(41, shifted(LINE, 900))]
        )
        with pytest.raises(ValueError, match=r"\b41\b"):
            rivers.read_segments(write(tmp_path / "r.geojson", doc))

    def test_a_point(self, rivers: ModuleType, tmp_path: Path) -> None:
        doc = collection([river(1, LINE), point(0.0, 0.0, objectid=2, elvid="1-1-1")])
        with pytest.raises(ValueError, match="Point"):
            rivers.read_segments(write(tmp_path / "r.geojson", doc))


class TestKind:
    """`kind` is total: `lake` when the casefolded `objekttype` starts with
    `innsj`, or when it is null or blank and `vatnlnr` is set; `river` for
    every other value. The raw value stays in the model, as served."""

    @pytest.mark.parametrize("objekttype", LAKE_SPELLINGS)
    @pytest.mark.parametrize("vatnlnr", [495, None], ids=["vatnlnr", "no-vatnlnr"])
    def test_every_lake_spelling_is_a_lake(
        self, rivers: ModuleType, tmp_path: Path, objekttype: str, vatnlnr: int | None
    ) -> None:
        (s,) = read(
            rivers, tmp_path / "r.geojson", [river(1, LINE, objekttype=objekttype, vatnlnr=vatnlnr)]
        )
        assert (s.kind, s.objekttype) == ("lake", objekttype)

    @pytest.mark.parametrize("objekttype", RIVER_SPELLINGS)
    @pytest.mark.parametrize("vatnlnr", [495, None], ids=["vatnlnr", "no-vatnlnr"])
    def test_every_other_value_is_a_river(
        self, rivers: ModuleType, tmp_path: Path, objekttype: str, vatnlnr: int | None
    ) -> None:
        """`vatnlnr` is set on 496 river features too, so it does not decide alone."""
        (s,) = read(
            rivers, tmp_path / "r.geojson", [river(1, LINE, objekttype=objekttype, vatnlnr=vatnlnr)]
        )
        assert (s.kind, s.objekttype) == ("river", objekttype)

    @pytest.mark.parametrize("objekttype", [NULL, BLANK, ""], ids=["null", "space", "empty"])
    @pytest.mark.parametrize(
        ("vatnlnr", "kind"), [(495, "lake"), (None, "river")], ids=["vatnlnr", "no-vatnlnr"]
    )
    def test_null_and_blank_follow_vatnlnr(
        self,
        rivers: ModuleType,
        tmp_path: Path,
        objekttype: str | None,
        vatnlnr: int | None,
        kind: str,
    ) -> None:
        (s,) = read(
            rivers, tmp_path / "r.geojson", [river(1, LINE, objekttype=objekttype, vatnlnr=vatnlnr)]
        )
        assert (s.kind, s.objekttype) == (kind, objekttype)

    def test_the_25_measured_values(self) -> None:
        assert len(LAKE_SPELLINGS) + len(RIVER_SPELLINGS) + len([NULL, BLANK]) == 25
        assert len(set(LAKE_SPELLINGS) | set(RIVER_SPELLINGS)) == 23

    @pytest.mark.parametrize("objekttype", ["INNSJØMIDTLINJE", "innsjoMIDTLINJE", "Innsj"])
    def test_case_does_not_matter(
        self, rivers: ModuleType, tmp_path: Path, objekttype: str
    ) -> None:
        (s,) = read(rivers, tmp_path / "r.geojson", [river(1, LINE, objekttype=objekttype)])
        assert s.kind == "lake"

    def test_a_value_merely_containing_innsj_is_a_river(
        self, rivers: ModuleType, tmp_path: Path
    ) -> None:
        """ "starts with", not "contains"."""
        (s,) = read(rivers, tmp_path / "r.geojson", [river(1, LINE, objekttype="ElvInnsjø")])
        assert s.kind == "river"


class TestExactCopies:
    """Step 4: within one `elvid`, segments whose vertex lists are equal to
    1 cm are one segment, the smallest `objectid` kept; the count dropped is
    reported (`drop_copies`, which `read_segments` applies)."""

    def ids(
        self, rivers: ModuleType, tmp_path: Path, features: Sequence[dict[str, Any]]
    ) -> list[int]:
        return [s.objectid for s in read(rivers, tmp_path / "r.geojson", features)]

    def test_the_smallest_objectid_is_kept_wherever_it_comes(
        self, rivers: ModuleType, tmp_path: Path
    ) -> None:
        features = [
            river(30, LINE),
            river(10, LINE),
            river(20, LINE),
            river(5, shifted(LINE, 1000)),
        ]
        assert sorted(self.ids(rivers, tmp_path, features)) == [5, 10]

    def test_vertices_within_a_centimetre_are_a_copy(
        self, rivers: ModuleType, tmp_path: Path
    ) -> None:
        """5 mm off in x and in y (7.1 mm apart) at coordinates up to 250 m."""
        features = [river(1, LINE), river(2, shifted(LINE, 0.005))]
        assert self.ids(rivers, tmp_path, features) == [1]

    def test_vertices_five_centimetres_off_are_not(
        self, rivers: ModuleType, tmp_path: Path
    ) -> None:
        features = [river(1, LINE), river(2, shifted(LINE, 0.05))]
        assert sorted(self.ids(rivers, tmp_path, features)) == [1, 2]

    def test_a_copy_is_found_at_utm_scale(self, rivers: ModuleType, tmp_path: Path) -> None:
        """Real coordinates, about 7e6 m north: a 5 mm offset is still a copy
        (float64 spacing there is 9.3e-10 m)."""
        utm = [(318563.63 + x, 6919942.64 + y) for x, y in LINE]
        features = [river(1, utm), river(2, shifted(utm, 0.005))]
        assert self.ids(rivers, tmp_path, features) == [1]

    def test_a_different_vertex_count_is_not_a_copy(
        self, rivers: ModuleType, tmp_path: Path
    ) -> None:
        longer = [LINE[0], (50.0, -5.0), *LINE[1:]]  # the same line, one vertex more
        assert sorted(self.ids(rivers, tmp_path, [river(1, LINE), river(2, longer)])) == [1, 2]

    def test_the_same_geometry_in_two_elvids_is_kept_twice(
        self, rivers: ModuleType, tmp_path: Path
    ) -> None:
        features = [river(1, LINE, elvid="1-1-1"), river(2, LINE, elvid="1-1-2")]
        assert sorted(self.ids(rivers, tmp_path, features)) == [1, 2]

    def test_two_geometries_under_one_strekninglnr_are_both_kept(
        self, rivers: ModuleType, tmp_path: Path
    ) -> None:
        """79.3.0's case: one `strekninglnr`, one `elvid`, two geometries 334 m apart."""
        features = [
            river(1, LINE, strekninglnr="79001", objekttype="InnsjøMidtlinje"),
            river(2, shifted(LINE, 334), strekninglnr="79001", objekttype="InnsjøMidtlinje"),
        ]
        assert sorted(self.ids(rivers, tmp_path, features)) == [1, 2]

    def test_the_count_dropped(self, rivers: ModuleType, tmp_path: Path) -> None:
        """82.4.0's nine pairs and 139.35.0's four groups of five: 25 extra
        features, as measured; plus three segments with no copy."""
        features: list[dict[str, Any]] = []
        oid = 1000
        for k in range(9):
            for _ in range(2):
                features.append(river(oid, shifted(LINE, 300 * k), elvid=f"82-{k}"))
                oid += 1
        for k in range(4):
            for _ in range(5):
                features.append(river(oid, shifted(LINE, 5000 + 300 * k), elvid=f"139-{k}"))
                oid += 1
        for k in range(3):
            features.append(river(oid, shifted(LINE, 9000 + 300 * k), elvid=f"x-{k}"))
            oid += 1
        assert len(features) == 9 * 2 + 4 * 5 + 3
        segments = read(rivers, tmp_path / "r.geojson", features)
        assert len(segments) == len(features) - 25 == 9 + 4 + 3
        assert rivers.drop_copies(segments)[1] == 0  # nothing left to drop

    def test_drop_copies_counts_what_it_drops(self, rivers: ModuleType, tmp_path: Path) -> None:
        """The pure function on segments built from distinct files, so the
        reader's own dropping does not hide the count."""
        segments = []
        for oid, elvid, d in [
            (7, "a", 0.0),
            (3, "a", 0.0),
            (5, "a", 0.0),
            (4, "b", 0.0),
            (9, "a", 50.0),
        ]:
            (s,) = read(
                rivers, tmp_path / f"{oid}.geojson", [river(oid, shifted(LINE, d), elvid=elvid)]
            )
            segments.append(s)
        kept, dropped = rivers.drop_copies(segments)
        assert sorted(s.objectid for s in kept) == [3, 4, 9]
        assert dropped == 2
