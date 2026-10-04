"""`rasputin fetch-stations nve-hrd`: the packaged list, the requests, the files (29, PR 3).

`docs/increments/29-nve-reference-catchments.md`, "The station set", "Data
use" and "The red suites", PR 3's `test_fetch_nve.py`. No network: the one
network call, `RangeClient.get_text`, is replaced on the class by
`nve_fixtures.FakeNve.get_text`, which answers from memory and logs every URL.
The User-Agent is tested on `fetch/http.py` itself, against a server on
127.0.0.1 that records the headers it is sent.

What is tested through the CLI and what through a function: the command is
the only entry the design names for the fetch, so the requests, the files and
the refusals are all observed through `rasputin fetch-stations`; the files are
read back through PR 3's own readers, `io/station_set.py` and `io/rivers.py`.
"""

from __future__ import annotations

import hashlib
import json
import re
import threading
from collections.abc import Iterator
from datetime import UTC, datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from nve_fixtures import (
    ALLOW,
    BLANK_TYPE,
    CHUNK,
    COPIES,
    COPY_HIGH,
    COPY_LOW,
    CRS,
    ENVELOPE_HALF,
    FORBIDDEN,
    LAYER_0_NAME,
    MS_2026_05_28,
    MULTI,
    NEWEST,
    NO_HIERARCHY,
    NULL_TYPE,
    OWNER,
    RENAMED,
    SAME_NUMBER,
    SHARED,
    STATION_NUMBER,
    TIE,
    FakeNve,
    build_fake,
    list_comments,
    list_rows,
)
from test_cli_mesh import plain
from tin_engine.cli import app
from tin_engine.fetch.http import RangeClient

REPO = Path(__file__).resolve().parents[2]
FILES = ("stations.geojson", "reference.geojson", "rivers.geojson", "NOTICE.txt")
GEOJSON = FILES[:3]
REFUSED, USAGE = 1, 2
runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})


def invoke(*args: str) -> tuple[int, str]:
    result = runner.invoke(app, list(args))
    return result.exit_code, plain(result.output)


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> FakeNve:
    """The fake service, installed as `RangeClient.get_text` for the test."""
    service = build_fake()
    monkeypatch.setattr(RangeClient, "get_text", lambda self, url: service.get_text(url))
    return service


def fetch(out: Path, *extra: str) -> tuple[int, str]:
    return invoke("fetch-stations", "nve-hrd", "--out-dir", str(out), *extra)


@pytest.fixture
def fetched(fake: FakeNve, tmp_path: Path) -> Path:
    """One successful fetch into `tmp_path/nve`."""
    out = tmp_path / "nve"
    code, output = fetch(out)
    assert code == 0, output
    assert fake.unexpected == [], fake.unexpected
    return out


def load(path: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return doc


def props_of(path: Path) -> list[dict[str, Any]]:
    return [f["properties"] for f in load(path)["features"]]


def by_station(path: Path) -> dict[str, dict[str, Any]]:
    return {f["properties"]["station"]: f for f in load(path)["features"]}


# --------------------------------------------------------------------------
# The packaged list
# --------------------------------------------------------------------------


class TestThePackagedList:
    """The 140 HRD rows: data, not code, read with `csv` (design: "The packaged list")."""

    def test_140_rows_of_four_columns(self) -> None:
        rows = list_rows()
        assert len(rows) == 140
        assert all(len(r) == 4 and None not in r for r in rows)

    def test_station_numbers_are_unique_and_well_formed(self) -> None:
        numbers = [r["station"] for r in list_rows()]
        assert len(set(numbers)) == len(numbers)
        assert [n for n in numbers if not STATION_NUMBER.match(n)] == []

    def test_series_versions_and_start_years_are_numbers(self) -> None:
        rows = list_rows()
        assert {r["series_version"] for r in rows} <= {"0", "1", "2"}
        assert all(re.fullmatch(r"(18|19|20)\d\d", r["hrd_start_daily"]) for r in rows)

    @pytest.mark.parametrize(
        "row",
        [
            # the design's own example, "The packaged list"
            {
                "station": "2.11.0",
                "series_version": "0",
                "name": "Narsjø",
                "hrd_start_daily": "1931",
            },
            # "2 142 1": the third column is the series version, not part of the number
            {
                "station": "2.142.0",
                "series_version": "1",
                "name": "Knappom",
                "hrd_start_daily": "1917",
            },
            # Table 1's last row
            {
                "station": "313.10.0",
                "series_version": "0",
                "name": "Magnor",
                "hrd_start_daily": "1981",
            },
        ],
        ids=["narsjo", "knappom", "magnor-last"],
    )
    def test_spot_rows(self, row: dict[str, str]) -> None:
        assert row in list_rows()

    def test_magnor_is_the_last_row(self) -> None:
        assert list_rows()[-1]["station"] == "313.10.0"

    def test_the_header_names_the_source_and_credits_nve(self) -> None:
        comments = "\n".join(list_comments())
        assert "norwegian-streamflow-reference-dataset-for-climate-change-studies" in comments
        assert "Kilde: NVE" in comments

    def test_notice_md_credits_nve_under_nlod(self) -> None:
        notice = (REPO / "NOTICE.md").read_text(encoding="utf-8")
        assert "Kilde: NVE" in notice and "NLOD" in notice


class TestTheCatalogue:
    def test_the_nve_hrd_entry(self) -> None:
        from tin_engine.sources import STATION_SOURCES, StationSource

        source = STATION_SOURCES["nve-hrd"]
        assert isinstance(source, StationSource)
        assert source.id == "nve-hrd"
        assert source.list_file == "nve_hrd_2025.csv"
        assert source.crs == "EPSG:25833"
        assert "Kilde: NVE" in source.credit
        assert "NLOD" in source.licence_note


# --------------------------------------------------------------------------
# The client identifies itself (fetch/http.py)
# --------------------------------------------------------------------------


class _Headers(BaseHTTPRequestHandler):
    seen: list[str]

    def log_message(self, format: str, *args: Any) -> None:
        pass

    def do_GET(self) -> None:
        self.seen.append(self.headers.get("User-Agent", ""))
        body = b"0123456789"
        if self.headers.get("Range"):
            self.send_response(206)
            self.send_header("Content-Range", f"bytes 0-9/{len(body)}")
        else:
            self.send_response(200)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


@pytest.fixture
def header_server() -> Iterator[tuple[str, list[str]]]:
    """`http://127.0.0.1:<port>/x` and the User-Agent of every request it gets."""
    seen: list[str] = []
    handler = type("Handler", (_Headers,), {"seen": seen})
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=httpd.serve_forever, kwargs={"poll_interval": 0.02})
    thread.start()
    try:
        yield f"http://127.0.0.1:{httpd.server_address[1]}/x", seen
    finally:
        httpd.shutdown()
        httpd.server_close()
        thread.join()


class TestTheUserAgent:
    """ "The client identifies itself (`User-Agent: rasputin/<version>`)
    instead of the library's default" ("Data use")."""

    def test_get_text_sends_rasputins_user_agent(
        self, header_server: tuple[str, list[str]]
    ) -> None:
        url, seen = header_server
        assert RangeClient(delays=()).get_text(url) == "0123456789"
        assert len(seen) == 1
        assert seen[0].startswith("rasputin/"), seen
        assert "Python-urllib" not in seen[0]

    def test_a_range_request_sends_it_too(self, header_server: tuple[str, list[str]]) -> None:
        url, seen = header_server
        assert RangeClient(delays=()).get(url, 0, 10).data == b"0123456789"
        assert len(seen) == 1 and seen[0].startswith("rasputin/"), seen


# --------------------------------------------------------------------------
# The requests
# --------------------------------------------------------------------------


class TestTheRequests:
    def test_only_the_three_named_layers_are_asked(self, fake: FakeNve, fetched: Path) -> None:
        assert fake.unexpected == []
        assert {c.layer for c in fake.calls} == {0, 2, 38}
        assert not any("hydapi" in c.url.lower() for c in fake.calls)

    @pytest.mark.parametrize("layer", [0, 38, 2])
    def test_out_fields_is_exactly_the_allow_list(
        self, fake: FakeNve, fetched: Path, layer: int
    ) -> None:
        calls = fake.calls_to(layer)
        assert calls
        for call in calls:
            assert "*" not in call.params.get("outFields", "*"), call.url
            assert sorted(call.out_fields) == sorted(ALLOW[layer]), call.url

    def test_no_forbidden_field_is_named_in_any_request(self, fake: FakeNve, fetched: Path) -> None:
        for call in fake.calls:
            assert not [f for f in FORBIDDEN if f in call.url.lower()], call.url

    def test_every_request_asks_for_geojson_in_25833(self, fake: FakeNve, fetched: Path) -> None:
        for call in fake.calls:
            assert call.params.get("outSR") == "25833", call.url
            assert call.params.get("f") == "geojson", call.url

    @pytest.mark.parametrize("layer", [0, 38])
    def test_points_and_polygons_go_in_chunks_of_40_in_list_order(
        self, fake: FakeNve, fetched: Path, layer: int
    ) -> None:
        numbers = [r["station"] for r in list_rows()]
        chunks = [c.stations for c in fake.calls_to(layer)]
        assert [len(c) for c in chunks] == [CHUNK, CHUNK, CHUNK, 140 - 3 * CHUNK]
        assert [n for c in chunks for n in c] == numbers

    def test_one_envelope_query_per_station_2_km_each_way(
        self, fake: FakeNve, fetched: Path
    ) -> None:
        """`map_radius + reach_up + 500 m` round the layer 0 point; the
        coordinates are exact sums of whole metres, so equality is exact."""
        calls = fake.calls_to(2)
        assert len(calls) == 140
        expected = [
            v
            for s in fake.stations
            for v in (
                s.x - ENVELOPE_HALF,
                s.y - ENVELOPE_HALF,
                s.x + ENVELOPE_HALF,
                s.y + ENVELOPE_HALF,
            )
        ]
        got = [v for c in calls for v in c.envelope]
        assert got == pytest.approx(expected, abs=1e-6)  # coordinates up to 7.0e6 m

    def test_requests_go_out_one_at_a_time(self, fake: FakeNve, fetched: Path) -> None:
        assert fake.max_active == 1

    def test_nothing_is_requested_when_the_files_exist(self, fake: FakeNve, fetched: Path) -> None:
        before = {name: (fetched / name).read_bytes() for name in (*FILES, "manifest.json")}
        fake.calls.clear()
        code, output = fetch(fetched)
        assert code == 0, output
        assert fake.calls == []
        assert {name: (fetched / name).read_bytes() for name in before} == before

    def test_refresh_fetches_again_and_writes_the_same_files(
        self, fake: FakeNve, fetched: Path
    ) -> None:
        """Deterministic order: the list file's, then `objectid`; only the
        manifest's fetch time may differ between two fetches."""
        before = {name: (fetched / name).read_bytes() for name in FILES}
        fake.calls.clear()
        code, output = fetch(fetched, "--refresh")
        assert code == 0, output
        assert len(fake.calls) == 4 + 4 + 140
        assert {name: (fetched / name).read_bytes() for name in FILES} == before


# --------------------------------------------------------------------------
# The files
# --------------------------------------------------------------------------


class TestStationsGeojson:
    def test_one_point_per_station_in_list_order_with_a_crs(self, fetched: Path) -> None:
        doc = load(fetched / "stations.geojson")
        assert doc["crs"]["properties"]["name"] == CRS
        assert [f["geometry"]["type"] for f in doc["features"]] == ["Point"] * 140
        assert [f["properties"]["station"] for f in doc["features"]] == [
            r["station"] for r in list_rows()
        ]

    def test_exactly_the_designed_properties(self, fetched: Path) -> None:
        """ "no other layer 0 field is copied"."""
        keys = {
            "station",
            "name",
            "series",
            "nve_area_km2",
            "hrd_start_daily",
            "watercourse",
            "river",
        }
        assert all(set(p) == keys for p in props_of(fetched / "stations.geojson"))

    def test_the_values(self, fake: FakeNve, fetched: Path) -> None:
        features = by_station(fetched / "stations.geojson")
        narsjo, knappom = features[NEWEST], features[COPIES]
        s = fake.station(NEWEST)
        assert narsjo["geometry"]["coordinates"] == pytest.approx([s.x, s.y], abs=1e-6)
        assert narsjo["properties"]["series"] == ["1001.0"]
        assert knappom["properties"]["series"] == ["1001.1"]
        assert narsjo["properties"]["hrd_start_daily"] == 1931
        assert narsjo["properties"]["name"] == "Narsjø"
        assert narsjo["properties"]["nve_area_km2"] == s.area_km2  # layer 0's, not the polygon's
        assert narsjo["properties"]["watercourse"] == s.watercourse

    def test_the_name_is_layer_0s_not_the_lists(self, fake: FakeNve, fetched: Path) -> None:
        """Ola's ruling of 2026-10-04: `name` is layer 0's `stasjonnavn`, as
        served. The fake serves a synthetic name for 311.4.0 that no list row
        holds, so the list's own spelling of it cannot make this pass."""
        assert fake.station(RENAMED).name == LAYER_0_NAME
        features = by_station(fetched / "stations.geojson")
        assert features[RENAMED]["properties"]["name"] == LAYER_0_NAME

    def test_river_is_the_first_name_of_the_hierarchy(self, fake: FakeNve, fetched: Path) -> None:
        features = by_station(fetched / "stations.geojson")
        i = [s.number for s in fake.stations].index(NEWEST)
        assert features[NEWEST]["properties"]["river"] == f"Elv{i}"
        assert features[NO_HIERARCHY]["properties"]["river"] is None


class TestReferenceGeojson:
    def test_one_polygon_per_station_in_list_order_with_a_crs(self, fetched: Path) -> None:
        doc = load(fetched / "reference.geojson")
        assert doc["crs"]["properties"]["name"] == CRS
        assert [f["properties"]["station"] for f in doc["features"]] == [
            r["station"] for r in list_rows()
        ]
        keys = {"station", "reference_area_km2", "reference_updated", "versions"}
        assert all(set(f["properties"]) == keys for f in doc["features"])

    def test_the_newest_version_wins_and_is_recorded(self, fetched: Path) -> None:
        """Served oldest, newest, middle: the newest (2026-05-28, 119.43 km²) wins."""
        props = by_station(fetched / "reference.geojson")[NEWEST]["properties"]
        assert props["reference_area_km2"] == 119.43
        assert props["versions"] == 3
        day = datetime.fromtimestamp(MS_2026_05_28 / 1000, UTC).date().isoformat()
        assert props["reference_updated"] == day == "2026-05-28"

    def test_a_tie_on_date_goes_to_the_larger_objectid(self, fetched: Path) -> None:
        """7102 (51.0 km²) is served before 7101 (50.0 km²), both dated alike."""
        props = by_station(fetched / "reference.geojson")[TIE]["properties"]
        assert props["reference_area_km2"] == 51.0
        assert props["versions"] == 2

    def test_a_single_version_counts_one(self, fetched: Path) -> None:
        assert by_station(fetched / "reference.geojson")[COPIES]["properties"]["versions"] == 1

    def test_a_multipart_polygon_stays_multipart(self, fetched: Path) -> None:
        geometry = by_station(fetched / "reference.geojson")[MULTI]["geometry"]
        assert geometry["type"] == "MultiPolygon" and len(geometry["coordinates"]) == 2


class TestRiversGeojson:
    def test_line_strings_by_objectid_with_a_crs(self, fetched: Path) -> None:
        doc = load(fetched / "rivers.geojson")
        assert doc["crs"]["properties"]["name"] == CRS
        assert {f["geometry"]["type"] for f in doc["features"]} == {"LineString"}
        ids = [f["properties"]["objectid"] for f in doc["features"]]
        assert ids == sorted(ids)
        assert all(set(f["properties"]) == set(ALLOW[2]) for f in doc["features"])

    def test_a_segment_seen_from_two_stations_is_kept_once(self, fetched: Path) -> None:
        ids = [p["objectid"] for p in props_of(fetched / "rivers.geojson")]
        assert len(ids) == len(set(ids))
        assert ids.count(SHARED) == 1

    def test_exact_copies_are_both_written(self, fetched: Path) -> None:
        """Dropping copies is the reader's job, so a user's file is cleaned
        the same way."""
        ids = {p["objectid"] for p in props_of(fetched / "rivers.geojson")}
        assert {COPY_LOW, COPY_HIGH, SAME_NUMBER} <= ids

    def test_objekttype_is_written_as_served_with_vatnlnr(self, fetched: Path) -> None:
        props = {p["objectid"]: p for p in props_of(fetched / "rivers.geojson")}
        assert props[NULL_TYPE]["objekttype"] is None and props[NULL_TYPE]["vatnlnr"] == 495
        assert props[BLANK_TYPE]["objekttype"] == " " and props[BLANK_TYPE]["vatnlnr"] == 0
        assert {p["objekttype"] for p in props.values()} >= {"SK", "InnsjoMidtlin", "ElvBekk"}


class TestNoticeAndManifest:
    def test_notice_credits_nve_under_nlod(self, fetched: Path) -> None:
        notice = (fetched / "NOTICE.txt").read_text(encoding="utf-8")
        assert "Kilde: NVE" in notice and "NLOD" in notice

    def test_no_owner_or_editor_field_reaches_any_file(self, fetched: Path) -> None:
        """The fake leaks `stasjoneier`, `globalid`, `oppdatertav` and a
        discharge normal into every answer; none may be copied."""
        for name in (*FILES, "manifest.json"):
            text = (fetched / name).read_text(encoding="utf-8")
            assert OWNER not in text, name
            assert not [f for f in FORBIDDEN if f in text], name

    def test_the_manifest_holds_each_files_sha256(self, fetched: Path) -> None:
        manifest = (fetched / "manifest.json").read_text(encoding="utf-8")
        for name in FILES:
            digest = hashlib.sha256((fetched / name).read_bytes()).hexdigest()
            assert name in manifest and digest in manifest, name

    def test_a_changed_file_no_longer_matches(self, fetched: Path) -> None:
        """The digest check above can fail: one byte changed, one digest gone."""
        path = fetched / "stations.geojson"
        path.write_bytes(path.read_bytes() + b"\n")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest not in (fetched / "manifest.json").read_text(encoding="utf-8")

    def test_the_manifest_holds_every_url_and_a_utc_time(
        self, fake: FakeNve, fetched: Path
    ) -> None:
        manifest = load(fetched / "manifest.json")
        strings = set(_strings(manifest))
        assert [c.url for c in fake.calls if c.url not in strings] == []
        utc = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(\.\d+)?(Z|\+00:00)$")
        assert any(utc.match(s) for s in strings), manifest


def _strings(value: Any) -> Iterator[str]:
    """Every string anywhere in a JSON value."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for v in value.values():
            yield from _strings(v)
    elif isinstance(value, list):
        for v in value:
            yield from _strings(v)


class TestReadBack:
    """The written files read back through PR 3's own readers."""

    def test_stations(self, fetched: Path) -> None:
        from tin_engine.io.station_set import read_stations

        stations, crs = read_stations(fetched / "stations.geojson")
        assert crs == CRS
        assert [s.station for s in stations] == [r["station"] for r in list_rows()]
        knappom = next(s for s in stations if s.station == COPIES)
        assert knappom.series == ("1001.1",) and knappom.name == "Knappom"
        assert next(s for s in stations if s.station == RENAMED).name == LAYER_0_NAME

    def test_references(self, fetched: Path) -> None:
        from shapely.geometry import MultiPolygon, Polygon

        from tin_engine.io.station_set import read_references

        references, crs = read_references(fetched / "reference.geojson")
        assert crs == CRS
        assert len(references) == 140
        assert isinstance(references[MULTI], MultiPolygon)
        assert isinstance(references[NEWEST], Polygon)

    def test_segments_drop_the_copy_and_keep_the_rest(self, fetched: Path) -> None:
        from tin_engine.io.rivers import read_segments

        segments, crs, dropped = read_segments(fetched / "rivers.geojson")
        assert crs == CRS
        assert dropped == 1  # COPY_HIGH, the one exact copy the fake serves
        ids = [s.objectid for s in segments]
        assert COPY_LOW in ids and COPY_HIGH not in ids and SAME_NUMBER in ids
        assert ids.count(SHARED) == 1
        kinds = {s.objectid: s.kind for s in segments}
        assert kinds[NULL_TYPE] == "lake" and kinds[SAME_NUMBER] == "lake"

    def test_a_blank_type_with_lake_number_0_reads_as_a_river(self, fetched: Path) -> None:
        """NVE sends `vatnlnr` 0 for "no lake"; the fake's blank-type segment
        carries it, as both blank-type features of the real layer do."""
        from tin_engine.io.rivers import read_segments

        segments, _, _ = read_segments(fetched / "rivers.geojson")
        assert {s.objectid: s.kind for s in segments}[BLANK_TYPE] == "river"


# --------------------------------------------------------------------------
# The refusals
# --------------------------------------------------------------------------


class TestTheRefusals:
    """Each refusal exits 1 and names the station and the reason."""

    def test_a_station_without_a_point(self, fake: FakeNve, tmp_path: Path) -> None:
        fake.drop_point.add(COPIES)
        code, output = fetch(tmp_path / "nve")
        assert code == REFUSED, output
        assert COPIES in output and "point" in output.lower(), output

    def test_a_station_without_a_polygon(self, fake: FakeNve, tmp_path: Path) -> None:
        fake.drop_polygon.add(COPIES)
        code, output = fetch(tmp_path / "nve")
        assert code == REFUSED, output
        assert COPIES in output and "polygon" in output.lower(), output

    def test_a_truncated_river_answer(self, fake: FakeNve, tmp_path: Path) -> None:
        """`exceededTransferLimit`: a truncated river never reaches `place`."""
        fake.truncate.add(MULTI)
        code, output = fetch(tmp_path / "nve")
        assert code == REFUSED, output
        assert MULTI in output, output
        assert re.search(r"exceededTransferLimit|limit|truncat", output, re.IGNORECASE), output

    @pytest.mark.parametrize(
        ("layer", "name"),
        [
            (0, "stasjonnr"),
            (38, "stasjonnr"),
            (38, "nedborfeltaareal_km2"),
            (38, "oppdateringsdato"),
            (2, "objectid"),
        ],
    )
    def test_a_field_the_service_left_out_is_named_in_plain_words(
        self, fake: FakeNve, tmp_path: Path, layer: int, name: str
    ) -> None:
        """A field the fetch needs, missing from every answer of its layer:
        the refusal names the field in a sentence, not as the bare `KeyError`
        text `'nedborfeltaareal_km2'`, and writes nothing."""
        fake.omit.add((layer, name))
        code, output = fetch(tmp_path / "nve")
        assert code == REFUSED, output
        assert name in output, output
        assert f"Error: '{name}'" not in output, output
        assert re.search(r"(?i)\b(missing|lacks|without|no)\b", output), output
        assert not (tmp_path / "nve").exists()

    def test_an_unknown_source(self, fake: FakeNve, tmp_path: Path) -> None:
        code, output = invoke("fetch-stations", "no-such-list", "--out-dir", str(tmp_path))
        assert code == USAGE, output
        assert "no-such-list" in output and "No such command" not in output, output
        assert fake.calls == []
