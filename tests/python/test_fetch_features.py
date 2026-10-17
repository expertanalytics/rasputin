"""`fetch/features.py`: polygons of an ArcGIS map-service layer, fetched by box into the cache.

`docs/increments/63-fetch-vectors.md`, sections 5.3 and 5.4, and section 6's
tests 2 to 8. No network: `arcgis_fixtures.FakeLayer.get_text` is the injected
`get_text`, and the cache is `tmp_path`. The layers are the hand-made ones of
`arcgis_fixtures` (ids 3, 7, 11, 12 in `box_all`; id 20 outside it), served
under the shipped `nve-lakes` and `clc2018` entries, whose CRS, id field and
fields they match.

The cache is passed as `CacheWriter(root, source.id)`, not yet entered;
`fetch_features` takes its lock itself, as `fetch/run.fetch` does.
"""

from __future__ import annotations

import datetime
import importlib
import json
import re
import sqlite3
from contextlib import closing
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import shapely

from arcgis_fixtures import CORINE, FAR, IDS, NVE, FakeLayer, LayerData, holds
from gpkg_fixtures import needs_rtree
from tin_engine.crs import transform_bounds
from tin_engine.domain import DomainPolygon
from tin_engine.fetch.http import FetchError
from tin_engine.fetch.plan import FetchRequest
from tin_engine.io.geopackage import decode_geometry, layer_info, layer_summary
from tin_engine.io.models import Bounds
from tin_engine.io.repository import CacheError, CacheRepository, CacheWriter, open_geopackage
from tin_engine.sources import SOURCES, RemoteSource

UTM33 = "EPSG:25833"


@pytest.fixture(scope="module")
def features() -> ModuleType:
    """`tin_engine.fetch.features`, imported at first use: a missing module
    fails each test, not the collection."""
    return importlib.import_module("tin_engine.fetch.features")


def source_of(data: LayerData) -> RemoteSource:
    return SOURCES["nve-lakes" if data.url == NVE.url else "clc2018"]


def bounds(b: tuple[float, float, float, float]) -> Bounds:
    return Bounds(x_min=b[0], y_min=b[1], x_max=b[2], y_max=b[3])


def request(
    b: tuple[float, float, float, float], source: str = "nve-lakes", **kw: Any
) -> FetchRequest:
    return FetchRequest(source=source, box=bounds(b), **kw)


def run(
    features: ModuleType,
    fake: FakeLayer,
    root: Path,
    b: tuple[float, float, float, float] | None = None,
    **kw: Any,
) -> Any:
    source = source_of(fake.data)
    asked = request(b or fake.data.box_all, source.id, **kw)
    return features.fetch_features(asked, source, fake.get_text, CacheWriter(root, source.id))


def stored(root: Path, source: str) -> dict[int, tuple[Any, dict[str, Any]]]:
    """The cached rows by id: geometry and the non-key, non-geometry columns."""
    with closing(open_geopackage(root / source / "features.gpkg")) as conn:
        info = layer_info(conn, None)
        names = [r[1] for r in conn.execute(f'PRAGMA table_info("{info.table}")')]
        rest = [n for n in names if n not in (info.pk, info.column)]
        select = ", ".join(f'"{n}"' for n in (info.pk, info.column, *rest))
        out = {}
        for row in conn.execute(f'SELECT {select} FROM "{info.table}"'):
            out[int(row[0])] = (decode_geometry(row[1])[1], dict(zip(rest, row[2:], strict=True)))
        return out


def stored_ids(root: Path, source: str) -> set[int]:
    path = root / source / "features.gpkg"
    return set(stored(root, source)) if path.is_file() else set()


def message(caught: pytest.ExceptionInfo[BaseException]) -> str:
    return str(caught.value)


def names_id(text: str, i: int) -> bool:
    return re.search(rf"(?<!\d){i}(?!\d)", text) is not None


# ------------------------------------------------------------ 2. pure helpers


class TestBoxInCrs:
    def test_a_box_in_the_source_crs_is_unchanged(self) -> None:
        plan = importlib.import_module("tin_engine.fetch.plan")
        given = (200_000.0, 6_600_000.0, 201_000.0, 6_601_000.0)
        got = plan.box_in_crs(request(given), UTM33)
        assert (got.x_min, got.y_min, got.x_max, got.y_max) == given

    def test_a_box_in_another_frame_is_moved_outward_with_no_growth(self) -> None:
        """A degree box (EPSG:4326) round Numedalslagen, moved into EPSG:25833:
        exactly `crs.transform_bounds` of it (densified, outward), so no growth
        is added; every corner's and edge midpoint's image lies inside."""
        plan = importlib.import_module("tin_engine.fetch.plan")
        given = (8.5, 59.0, 10.0, 60.5)
        got = plan.box_in_crs(request(given, out_crs="EPSG:4326"), UTM33)
        moved = transform_bounds("EPSG:4326", UTM33, given)
        assert (got.x_min, got.y_min, got.x_max, got.y_max) == moved
        from tin_engine.crs import reprojector

        x0, y0, x1, y1 = given
        points = [(x, y) for x in (x0, (x0 + x1) / 2, x1) for y in (y0, (y0 + y1) / 2, y1)]
        for x, y in reprojector("EPSG:4326", UTM33)(points):
            assert got.x_min <= x <= got.x_max and got.y_min <= y <= got.y_max

    def test_a_domain_gives_its_bounds(self) -> None:
        plan = importlib.import_module("tin_engine.fetch.plan")
        ring = [(200_000.0, 6_600_000.0), (201_000.0, 6_600_000.0), (200_500.0, 6_601_500.0)]
        domain = DomainPolygon(polygon=shapely.Polygon(ring), crs=UTM33)
        got = plan.box_in_crs(FetchRequest(source="nve-lakes", domain=domain), UTM33)
        assert (got.x_min, got.y_min, got.x_max, got.y_max) == (
            200_000.0, 6_600_000.0, 201_000.0, 6_601_500.0
        )  # fmt: skip

    def test_a_box_crossing_the_antimeridian_is_refused(self) -> None:
        plan = importlib.import_module("tin_engine.fetch.plan")
        with pytest.raises(FetchError, match="±180°"):
            plan.box_in_crs(request((179.0, 10.0, 181.0, 11.0)), "EPSG:4326")


class TestChunks:
    def test_seven_ids_in_threes(self, features: ModuleType) -> None:
        ids = [1, 2, 3, 5, 8, 13, 21]
        assert features.chunks(ids, 3) == [[1, 2, 3], [5, 8, 13], [21]]

    def test_no_ids_no_chunks(self, features: ModuleType) -> None:
        assert features.chunks([], 3) == []

    def test_one_id_one_chunk(self, features: ModuleType) -> None:
        assert features.chunks([9], 500) == [[9]]

    def test_the_cap_is_500_ids(self, features: ModuleType) -> None:
        assert features.MAX_IDS == 500


# ------------------------------------------------------------ 3. the layer check


class TestLayerCheck:
    def test_the_facts_of_the_nve_layer(self, features: ModuleType) -> None:
        facts = features.layer_facts(source_of(NVE), FakeLayer(NVE).get_text)
        assert facts.geometry == "esriGeometryPolygon"
        assert facts.id_field == "objectid"
        assert facts.max_records == 2000
        assert facts.epsg == 25833
        assert facts.fields["areal_km2"] == "esriFieldTypeDouble"

    def test_wkid_is_taken_when_no_latest_wkid_is_stated(self, features: ModuleType) -> None:
        facts = features.layer_facts(source_of(CORINE), FakeLayer(CORINE, wkid_only=True).get_text)
        assert facts.epsg == 3035

    @pytest.mark.parametrize(
        ("knobs", "words"),
        [
            ({"geometry_type": "esriGeometryPolyline"}, ("only polygon layers are fetched",)),
            ({"stated_epsg": 4326}, ("4326", "25833")),
            ({"missing_fields": ("navn",)}, ("navn",)),
            ({"hashed_ids": True}, ("OIDFieldContainsHashValue",)),
        ],
        ids=["a-line-layer", "a-crs-mismatch", "a-field-the-layer-lacks", "hashed-ids"],
    )
    def test_refused_before_any_query(
        self, features: ModuleType, tmp_path: Path, knobs: dict[str, Any], words: tuple[str, ...]
    ) -> None:
        fake = FakeLayer(NVE, **knobs)
        with pytest.raises(FetchError) as caught:
            run(features, fake, tmp_path)
        for word in words:
            assert word in message(caught), message(caught)
        assert "nve-lakes" in message(caught) and NVE.url in message(caught), message(caught)
        assert len(fake.urls) == 1 and "/query" not in fake.urls[0], fake.urls
        assert stored_ids(tmp_path, "nve-lakes") == set()


# ------------------------------------------------------------ 4. the happy path


@needs_rtree
class TestHappyPath:
    @pytest.fixture
    def fetched(self, features: ModuleType, tmp_path: Path) -> tuple[FakeLayer, Any, Path]:
        fake = FakeLayer(NVE, max_records=2)
        report = run(features, fake, tmp_path)
        return fake, report, tmp_path

    @pytest.fixture
    def days(self) -> set[str]:
        """Today, and yesterday for a run that crossed midnight: the fetch's date is one of them."""
        today = datetime.date.today()
        return {today.isoformat(), (today - datetime.timedelta(days=1)).isoformat()}

    def test_the_requests(self, fetched: tuple[FakeLayer, Any, Path]) -> None:
        """The layer, the ids of the box, then two chunks of two ids (the
        layer's `maxRecordCount` of 2 being below 500), ascending."""
        fake, _, _ = fetched
        assert len(fake.urls) == 4
        assert fake.urls[0].split("?")[0] == NVE.url
        (ids,) = fake.ids_queries()
        assert ids["inSR"] == "25833" and ids["where"] == "1=1"
        assert tuple(float(v) for v in ids["geometry"].split(",")) == NVE.box_all
        assert fake.chunks_asked() == [[3, 7], [11, 12]]
        for q in fake.chunk_queries():
            assert set(q["outFields"].split(",")) == {"objectid", *NVE.asked}
            assert q["outSR"] == "25833"
            assert q["returnGeometry"] == "true" and q["f"] == "geojson"

    def test_the_report(self, fetched: tuple[FakeLayer, Any, Path]) -> None:
        fake, report, _ = fetched
        assert (report.in_box, report.present, report.fetched) == (4, 0, 4)
        assert report.requests == len(fake.urls) == 4
        assert 0 < report.bytes <= fake.served
        assert report.seconds >= 0

    def test_four_rows_keyed_by_id_holding_the_allow_list_only(
        self, fetched: tuple[FakeLayer, Any, Path]
    ) -> None:
        _, _, root = fetched
        rows = stored(root, "nve-lakes")
        assert set(rows) == set(IDS)
        for i, (geometry, values) in rows.items():
            assert set(values) == set(NVE.asked), values
            assert values == {k: NVE.properties[i][k] for k in NVE.asked}, i
            assert geometry.geom_type == "MultiPolygon", i
            served = NVE.shapes[i]
            parts = served.geoms if served.geom_type == "MultiPolygon" else [served]
            assert geometry.equals_exact(shapely.MultiPolygon(list(parts)), 0.0), i

    def test_the_bow_tie_is_stored_as_served_unrepaired(
        self, fetched: tuple[FakeLayer, Any, Path]
    ) -> None:
        _, _, root = fetched
        geometry = stored(root, "nve-lakes")[11][0]
        assert not geometry.is_valid
        ring = list(geometry.geoms[0].exterior.coords)
        assert ring == list(NVE.shapes[11].exterior.coords)

    def test_the_layer_declares_multipolygon_and_its_extents_are_the_polygons(
        self, fetched: tuple[FakeLayer, Any, Path]
    ) -> None:
        _, _, root = fetched
        with closing(open_geopackage(root / "nve-lakes" / "features.gpkg")) as conn:
            info = layer_info(conn, None)
            summary = layer_summary(conn, info)
            assert info.pk == "objectid"
            assert info.crs == UTM33
            assert info.rtree is not None
            assert summary.geometry == "MULTIPOLYGON"
            assert summary.features == 4
            assert set(summary.attributes) == set(NVE.asked)
            assert summary.extent == NVE.extent(IDS)
            contents = conn.execute(
                "SELECT min_x, min_y, max_x, max_y FROM gpkg_contents WHERE table_name = ?",
                (info.table,),
            ).fetchone()
            assert tuple(contents) == NVE.extent(IDS)
            types = {r[1]: r[2].upper() for r in conn.execute(f'PRAGMA table_info("{info.table}")')}
            assert types["vatnlnr"] == "INTEGER"
            assert types["navn"] == "TEXT"
            assert types["areal_km2"] == "REAL"

    def test_the_manifest_is_a_features_manifest(
        self, fetched: tuple[FakeLayer, Any, Path], days: set[str]
    ) -> None:
        _, _, root = fetched
        doc = json.loads((root / "nve-lakes" / "manifest.json").read_text(encoding="utf-8"))
        assert doc["kind"] == "arcgis-features"
        assert doc["source"] == "nve-lakes"
        assert doc["url"] == NVE.url
        assert doc["id_field"] == "objectid"
        assert list(doc["fields"]) == list(NVE.asked)
        assert len(doc["requests"]) == 1
        assert doc["first_fetched"] in days
        repository = importlib.import_module("tin_engine.io.repository")
        manifest = repository.FeatureManifest.model_validate_json(
            (root / "nve-lakes" / "manifest.json").read_bytes()
        )
        assert manifest.source == "nve-lakes"

    def test_the_notice_holds_the_credit_with_the_date(
        self, fetched: tuple[FakeLayer, Any, Path], days: set[str]
    ) -> None:
        _, _, root = fetched
        text = (root / "nve-lakes" / "NOTICE.txt").read_text(encoding="utf-8")
        source = source_of(NVE)
        credits = {
            source.credit_for(url=NVE.url, date=datetime.date.fromisoformat(d)) for d in days
        }
        assert any(c in text for c in credits), text
        assert source.licence_note in text

    def test_the_raster_reader_refuses_the_features_cache_in_plain_words(
        self, fetched: tuple[FakeLayer, Any, Path]
    ) -> None:
        _, _, root = fetched
        with pytest.raises(CacheError) as caught:
            CacheRepository(root, "nve-lakes").footprints()
        assert "nve-lakes holds features; pass its features.gpkg to a vector option" in message(
            caught
        )


@needs_rtree
class TestChunkSize:
    def test_chunks_are_capped_at_500_ids_below_max_record_count(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        """1001 squares of 10 m in a row, a layer whose `maxRecordCount` is
        2000: chunks of 500, 500 and 1."""
        many = {
            i: shapely.box(200_000 + 20 * i, 6_600_000, 200_010 + 20 * i, 6_600_010)
            for i in range(1, 1002)
        }
        props = {i: {"vatnlnr": i, "navn": None, "areal_km2": 0.0001} for i in many}
        data = LayerData(
            NVE.url, 25833, "objectid", NVE.fields, many, props,
            box_all=(199_000, 6_599_000, 230_000, 6_601_000), box_left=NVE.box_left,
            box_right=NVE.box_right, asked=NVE.asked,
        )  # fmt: skip
        fake = FakeLayer(data, max_records=2000)
        report = run(features, fake, tmp_path)
        assert [len(c) for c in fake.chunks_asked()] == [500, 500, 1]
        assert report.fetched == 1001


# ------------------------------------------------------------ 5. resume and overlap


@needs_rtree
class TestResumeAndOverlap:
    def test_an_overlapping_box_asks_only_the_ids_not_cached(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        run(features, FakeLayer(NVE, max_records=2), tmp_path, NVE.box_left)
        assert stored_ids(tmp_path, "nve-lakes") == {3, 7}
        second = FakeLayer(NVE, max_records=2)
        report = run(features, second, tmp_path, NVE.box_right)
        assert second.chunks_asked() == [[11, 12]]
        assert (report.in_box, report.present, report.fetched) == (3, 1, 2)
        assert stored_ids(tmp_path, "nve-lakes") == set(IDS)

    def test_an_interrupted_fetch_keeps_the_first_chunk_and_the_rerun_asks_the_rest(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        """The second chunk (11, 12) is refused down to id 11 alone; the
        first chunk's rows stay, and a rerun asks only what is missing."""
        broken = FakeLayer(NVE, max_records=2, refuse=holds(11))
        with pytest.raises(FetchError):
            run(features, broken, tmp_path)
        kept = stored_ids(tmp_path, "nve-lakes")
        assert {3, 7} <= kept and 11 not in kept
        again = FakeLayer(NVE, max_records=2)
        run(features, again, tmp_path)
        asked = [i for chunk in again.chunks_asked() for i in chunk]
        assert sorted(asked) == sorted(set(IDS) - kept)
        assert stored_ids(tmp_path, "nve-lakes") == set(IDS)

    def test_the_same_box_twice_asks_no_chunk_and_logs_one_request(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        run(features, FakeLayer(NVE), tmp_path)
        second = FakeLayer(NVE)
        report = run(features, second, tmp_path)
        assert second.chunks_asked() == []
        assert (report.in_box, report.present, report.fetched) == (4, 4, 0)
        doc = json.loads((tmp_path / "nve-lakes" / "manifest.json").read_text(encoding="utf-8"))
        assert len(doc["requests"]) == 1

    def test_two_boxes_log_two_requests(self, features: ModuleType, tmp_path: Path) -> None:
        run(features, FakeLayer(NVE), tmp_path, NVE.box_left)
        run(features, FakeLayer(NVE), tmp_path, NVE.box_right)
        doc = json.loads((tmp_path / "nve-lakes" / "manifest.json").read_text(encoding="utf-8"))
        assert len(doc["requests"]) == 2

    def test_refresh_discards_the_cache_and_asks_every_id_again(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        run(features, FakeLayer(NVE), tmp_path, NVE.box_left)
        second = FakeLayer(NVE)
        report = run(features, second, tmp_path, NVE.box_right, refresh=True)
        assert [i for c in second.chunks_asked() for i in c] == [7, 11, 12]
        assert report.present == 0
        assert stored_ids(tmp_path, "nve-lakes") == {7, 11, 12}

    def test_a_box_meeting_nothing_fetches_nothing(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        """The services answer `"objectIds": null` for an empty box."""
        empty = (100_000, 6_000_000, 101_000, 6_001_000)
        fake = FakeLayer(NVE)
        report = run(features, fake, tmp_path, empty)
        assert fake.chunks_asked() == []
        assert (report.in_box, report.present, report.fetched) == (0, 0, 0)
        assert stored_ids(tmp_path, "nve-lakes") == set()


# ------------------------------------------------------------ 6. halving


@needs_rtree
class TestHalving:
    def test_a_refused_chunk_is_halved_until_answered(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        """Any chunk holding both 11 and 12 is refused, as CORINE's two sea
        polygons in one answer were: (3, 7, 11, 12) splits into (3, 7) and
        (11, 12), which splits into (11) and (12)."""
        fake = FakeLayer(CORINE, max_records=4, refuse=holds(11, 12))
        report = run(features, fake, tmp_path)
        assert fake.chunks_asked() == [[3, 7, 11, 12], [3, 7], [11, 12], [11], [12]]
        assert stored_ids(tmp_path, "clc2018") == set(IDS)
        assert report.fetched == 4

    def test_one_id_refused_alone_is_the_fetchs_error_naming_it(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        fake = FakeLayer(CORINE, max_records=4, refuse=holds(12))
        with pytest.raises(FetchError) as caught:
            run(features, fake, tmp_path)
        assert [12] in fake.chunks_asked()
        assert names_id(message(caught), 12), message(caught)
        assert "clc2018" in message(caught), message(caught)


# ------------------------------------------------------------ 7. answer checks


@needs_rtree
class TestAnswerChecks:
    def test_an_id_missing_from_an_answer(self, features: ModuleType, tmp_path: Path) -> None:
        with pytest.raises(FetchError) as caught:
            run(features, FakeLayer(NVE, drop=frozenset({7})), tmp_path)
        assert "the layer changed during the fetch; run again" in message(caught)

    def test_an_id_not_asked_in_an_answer(self, features: ModuleType, tmp_path: Path) -> None:
        with pytest.raises(FetchError) as caught:
            run(features, FakeLayer(NVE, extra=frozenset({FAR})), tmp_path)
        assert "the layer changed during the fetch; run again" in message(caught)

    def test_an_answer_whose_crs_is_not_out_sr(self, features: ModuleType, tmp_path: Path) -> None:
        with pytest.raises(FetchError) as caught:
            run(features, FakeLayer(NVE, crs_name="EPSG:4326"), tmp_path)
        assert "EPSG:4326" in message(caught), message(caught)
        assert stored_ids(tmp_path, "nve-lakes") == set()

    def test_a_feature_with_null_geometry(self, features: ModuleType, tmp_path: Path) -> None:
        with pytest.raises(FetchError) as caught:
            run(features, FakeLayer(NVE, raw_geometry={7: None}), tmp_path)
        assert names_id(message(caught), 7), message(caught)

    def test_a_feature_that_is_not_a_polygon(self, features: ModuleType, tmp_path: Path) -> None:
        line = {"type": "LineString", "coordinates": [[201000.0, 6600000.0], [201100.0, 6600100.0]]}
        with pytest.raises(FetchError) as caught:
            run(features, FakeLayer(NVE, raw_geometry={7: line}), tmp_path)
        assert names_id(message(caught), 7), message(caught)

    @pytest.mark.parametrize("bad", [float("nan"), float("inf")], ids=["nan", "infinity"])
    def test_a_feature_with_a_non_finite_coordinate(
        self, features: ModuleType, tmp_path: Path, bad: float
    ) -> None:
        """Python's `json` reads the `NaN` and `Infinity` tokens a server may write."""
        ring = [
            [201000.0, 6600000.0],
            [bad, 6600000.0],
            [201200.0, 6600200.0],
            [201000.0, 6600000.0],
        ]
        with pytest.raises(FetchError) as caught:
            run(
                features,
                FakeLayer(NVE, raw_geometry={7: {"type": "Polygon", "coordinates": [ring]}}),
                tmp_path,
            )
        assert names_id(message(caught), 7), message(caught)
        assert 7 not in stored_ids(tmp_path, "nve-lakes")

    def test_exceeded_transfer_limit_on_the_ids_answer(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        fake = FakeLayer(NVE, exceeded=True)
        with pytest.raises(FetchError) as caught:
            run(features, fake, tmp_path)
        assert "exceededTransferLimit" in message(caught), message(caught)
        assert fake.chunks_asked() == []

    def test_a_year_is_refused(self, features: ModuleType, tmp_path: Path) -> None:
        fake = FakeLayer(NVE)
        with pytest.raises(FetchError, match="has no years"):
            run(features, fake, tmp_path, year=2020)


# ------------------------------------------------------------ 8. dry run


class TestDryRun:
    def test_only_the_layer_and_ids_are_asked_and_nothing_is_written(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        fake = FakeLayer(NVE, max_records=2)
        report = run(features, fake, tmp_path, dry_run=True)
        assert len(fake.urls) == 2 and fake.chunks_asked() == []
        assert (report.in_box, report.present) == (4, 0)
        assert not (tmp_path / "nve-lakes").exists()

    @needs_rtree
    def test_a_dry_run_counts_the_ids_already_cached(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        run(features, FakeLayer(NVE), tmp_path, NVE.box_left)
        before = sorted(p.name for p in (tmp_path / "nve-lakes").iterdir())
        report = run(features, FakeLayer(NVE), tmp_path, dry_run=True)
        assert (report.in_box, report.present) == (4, 2)
        assert sorted(p.name for p in (tmp_path / "nve-lakes").iterdir()) == before


# ------------------------------------------------------------ the box, end to end


class TestTheBoxSent:
    def test_a_box_in_degrees_is_sent_in_the_source_crs(
        self, features: ModuleType, tmp_path: Path
    ) -> None:
        given = (9.0, 59.4, 9.1, 59.5)
        fake = FakeLayer(NVE)
        run(features, fake, tmp_path, given, out_crs="EPSG:4326", dry_run=True)
        (ids,) = fake.ids_queries()
        sent = tuple(float(v) for v in ids["geometry"].split(","))
        # 1 mm: the URL may write the box to the millimetre. Coordinates in
        # metres up to 6.6e6 (northings), checked at that.
        assert sent == pytest.approx(transform_bounds("EPSG:4326", UTM33, given), abs=1e-3)
        assert ids["inSR"] == "25833"


# ------------------------------------------------------------ the features table


class TestFeaturesTable:
    def test_a_field_type_with_no_column_type_is_refused(self) -> None:
        write = importlib.import_module("tin_engine.io.geopackage_write")
        from tin_engine.crs import srs_entry

        with closing(sqlite3.connect(":memory:")) as conn:
            write.create_geopackage(conn, srs_entry(UTM33))
            with pytest.raises(ValueError, match="esriFieldTypeBlob"):
                write.create_features_table(
                    conn, "features", 25833, "MULTIPOLYGON", "objectid",
                    [("vatnlnr", "esriFieldTypeInteger"), ("raw", "esriFieldTypeBlob")],
                )  # fmt: skip
