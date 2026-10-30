"""The `arcgis-features` kind: the catalogue entries, the listing, the readers, the CLI.

`docs/increments/63-fetch-vectors.md`, sections 5.2, 5.5 and 5.6, and section
6's tests 1 and 9 to 11. No network: `RangeClient.get_text` is replaced on the
class by `arcgis_fixtures.FakeLayer.get_text`, and the cache is `tmp_path`.

Ola's two open questions are taken at their defaults: the two entries ship,
and the CORINE credit names no version (the map service states none).
"""

from __future__ import annotations

import datetime
import importlib
import json
import os
from pathlib import Path
from typing import Any

import pytest
import shapely

from arcgis_fixtures import CODES, CORINE, CORINE_URL, IDS, NVE, NVE_URL, FakeLayer, LayerData
from cli_helpers import USAGE, invoke
from crs_fixtures import crs_box
from gpkg_fixtures import needs_rtree
from tin_engine.feature_input import read_lake_polygons, read_source
from tin_engine.fetch.http import RangeClient
from tin_engine.fetch.plan import FetchRequest
from tin_engine.io.repository import CacheWriter
from tin_engine.sources import SOURCES, RemoteSource, catalogue, notice

CREDIT = "A test credit, {url}, {date}"
LICENCE = "A test licence note"
REFUSED = 1


def arcgis(**changes: Any) -> RemoteSource:
    """An `arcgis-features` entry; a change to None drops that key."""
    table: dict[str, Any] = {
        "id": "a",
        "kind": "arcgis-features",
        "url": NVE_URL,
        "crs": "EPSG:25833",
        "fields": ("vatnlnr", "navn"),
        "credit": CREDIT,
        "licence_note": LICENCE,
    }
    table.update(changes)
    return RemoteSource(**{k: v for k, v in table.items() if v is not None})


def refusal(**changes: Any) -> str:
    with pytest.raises(ValueError) as caught:
        arcgis(**changes)
    return str(caught.value)


def fetch_into(root: Path, data: LayerData) -> None:
    """The layer's four polygons fetched into `root` under its shipped entry."""
    features = importlib.import_module("tin_engine.fetch.features")
    source = SOURCES["nve-lakes" if data.url == NVE.url else "clc2018"]
    x0, y0, x1, y1 = data.box_all
    request = FetchRequest(source=source.id, box=crs_box(source.crs, (x0, y0, x1, y1)))
    features.fetch_features(request, source, FakeLayer(data).get_text, CacheWriter(root, source.id))


def listed(cache: Path) -> dict[str, dict[str, Any]]:
    from tin_engine.source_report import listing

    environ = {"XDG_CONFIG_HOME": os.environ["XDG_CONFIG_HOME"]}  # no user sources file
    report = listing(environ, cache, None)
    return {s["id"]: s for s in report.model_dump(mode="json")["sources"]}


def fits_of(entry: dict[str, Any]) -> set[tuple[str, str, str | None]]:
    return {(f["command"], f["option"], f.get("classes")) for f in entry["fits"]}


# ------------------------------------------------------------ 1. the catalogue


class TestTheShippedEntries:
    def test_nve_lakes(self) -> None:
        entry = SOURCES["nve-lakes"]
        assert entry.kind == "arcgis-features"
        assert entry.values == "features"
        assert entry.url == NVE_URL
        assert entry.crs == "EPSG:25833"
        assert entry.fields == ("vatnlnr", "navn", "areal_km2")
        assert "NLOD" in entry.credit and "NVE" in entry.credit
        assert entry.url_for() == NVE_URL

    def test_clc2018(self) -> None:
        entry = SOURCES["clc2018"]
        assert entry.kind == "arcgis-features"
        assert entry.values == "features"
        assert entry.url == CORINE_URL
        assert entry.crs == "EPSG:3035"
        assert entry.fields == ("ID", "Code_18")
        assert entry.classes == "corine"

    def test_the_corine_credit_names_no_version_and_no_doi(self) -> None:
        """The service states no version; the DOI is the V2020_20u1 release's."""
        text = SOURCES["clc2018"].credit + SOURCES["clc2018"].description
        for word in ("2020_20u1", "V2020", "doi.org", "10.2909"):
            assert word not in text, word

    def test_the_corine_credit_states_modification_and_no_endorsement(self) -> None:
        credit = SOURCES["clc2018"].credit
        assert "Modified" in credit and "Not endorsed by the European Union" in credit

    @pytest.mark.parametrize("key", ["nve-lakes", "clc2018"])
    def test_the_credit_takes_the_url_and_the_date(self, key: str) -> None:
        entry = SOURCES[key]
        day = datetime.date(2026, 10, 17)
        filled = entry.credit_for(url=str(entry.url), date=day)
        assert "2026-10-17" in filled and str(entry.url) in filled
        text = notice(entry, (filled,))
        assert filled in text and entry.licence_note in text

    @pytest.mark.parametrize("key", ["nve-lakes", "clc2018"])
    def test_no_year_applies(self, key: str) -> None:
        with pytest.raises(ValueError, match="has no years"):
            SOURCES[key].url_for(2020)

    def test_the_entries_are_in_the_packaged_toml(self) -> None:
        sources = importlib.import_module("tin_engine.sources")
        text = (Path(sources.__file__).parent / "data" / "sources.toml").read_text(encoding="utf-8")
        assert "[source.nve-lakes]" in text and "[source.clc2018]" in text


class TestTheKindsRules:
    def test_the_kind_gives_features(self) -> None:
        assert arcgis().values == "features"

    def test_the_url_is_the_file(self) -> None:
        assert arcgis().url_for() == NVE_URL

    @pytest.mark.parametrize(
        ("field", "changes"),
        [
            ("fields", {"fields": None}),
            ("fields", {"fields": ()}),
            ("fields", {"fields": ("navn", "navn")}),
            ("url", {"url": None}),
            ("crs", {"crs": None}),
            ("layer", {"layer": "features"}),
            ("z_unit", {"z_unit": "metre"}),
            ("path", {"path": "/data/lakes.gpkg"}),
            ("years", {"years": (2018,)}),
            ("url_template", {"url_template": "https://example.invalid/{year}"}),
            ("nodata", {"nodata": -9999.0}),
            ("tile_list_url", {"tile_list_url": "https://example.invalid/list"}),
            ("values", {"values": "heights"}),
        ],
        ids=lambda v: v if isinstance(v, str) else None,
    )
    def test_refused_naming_the_field(self, field: str, changes: dict[str, Any]) -> None:
        text = refusal(**changes)
        assert f"{field}:" in text, text

    def test_fields_are_refused_for_a_raster_kind(self) -> None:
        with pytest.raises(ValueError) as caught:
            RemoteSource(
                id="r", kind="one-cog", url="https://example.invalid/x.tif", crs="EPSG:25833",
                fields=("a",), credit=CREDIT, licence_note=LICENCE,
            )  # fmt: skip
        assert "fields:" in str(caught.value), str(caught.value)

    def test_fields_are_refused_for_a_local_vector(self) -> None:
        with pytest.raises(ValueError) as caught:
            RemoteSource(
                id="l", kind="local", path=Path("/data/l.gpkg"), fields=("a",),
                credit=CREDIT, licence_note=LICENCE,
            )  # fmt: skip
        assert "fields:" in str(caught.value), str(caught.value)

    def test_a_users_file_may_name_the_kind(self, tmp_path: Path) -> None:
        file = tmp_path / "sources.toml"
        file.write_text(
            "[source.my-lakes]\n"
            'kind = "arcgis-features"\n'
            f'url = "{NVE_URL}"\n'
            'crs = "EPSG:25833"\n'
            'fields = ["vatnlnr"]\n'
            f'credit = "{CREDIT}"\n'
            f'licence_note = "{LICENCE}"\n',
            encoding="utf-8",
        )
        got = catalogue(file).sources["my-lakes"]
        assert (got.kind, got.values, got.fields) == ("arcgis-features", "features", ("vatnlnr",))


# ------------------------------------------------------------ 9. the listing


@needs_rtree
class TestTheListing:
    def test_before_a_fetch_not_present_and_no_fits(self, tmp_path: Path) -> None:
        found = listed(tmp_path)
        for key in ("nve-lakes", "clc2018"):
            assert found[key]["present"] is False, key
            assert found[key]["fits"] == [], key
            assert found[key]["values"] == "features"

    def test_after_a_fetch_clc2018_is_listed_from_its_cache(self, tmp_path: Path) -> None:
        fetch_into(tmp_path, CORINE)
        entry = listed(tmp_path)["clc2018"]
        assert entry["present"] is True
        assert entry["path"] == str(tmp_path / "clc2018" / "features.gpkg")
        assert entry["features"] == 4
        assert entry["geometry"] == ["MultiPolygon"]
        assert set(entry["attributes"]) == {"ID", "Code_18"}
        assert ("mesh", "--land-cover", "corine") in fits_of(entry)
        assert ("catchment", "--lakes", None) in fits_of(entry)
        assert {f[0] for f in fits_of(entry)} <= {"mesh", "catchment", "network-catchment"}
        (coded,) = [c for c in entry["classes"] if c["attribute"] == "Code_18"]
        assert coded["values"] == sorted(CODES[i] for i in IDS)

    def test_after_a_fetch_nve_lakes_fits_catchment_lakes(self, tmp_path: Path) -> None:
        fetch_into(tmp_path, NVE)
        entry = listed(tmp_path)["nve-lakes"]
        assert entry["present"] is True
        assert entry["path"] == str(tmp_path / "nve-lakes" / "features.gpkg")
        assert entry["features"] == 4
        assert ("catchment", "--lakes", None) in fits_of(entry)
        assert "station-catchments" not in {f[0] for f in fits_of(entry)}


# ------------------------------------------------------------ 10. through the readers


@needs_rtree
class TestThroughTheReaders:
    def test_catchment_lakes_reads_the_cached_lakes(self, tmp_path: Path) -> None:
        """A point inside lake 3 (a 100 m square) finds that lake alone."""
        fetch_into(tmp_path, NVE)
        x0, y0, _, _ = NVE.shapes[3].bounds
        lakes, crs = read_lake_polygons(
            tmp_path / "nve-lakes" / "features.gpkg", None, (x0 + 50, y0 + 50), "EPSG:25833"
        )
        assert crs == "EPSG:25833"
        assert len(lakes) == 1
        assert shapely.equals(lakes[0], NVE.shapes[3])

    def test_land_cover_with_corine_reads_the_cached_codes(self, tmp_path: Path) -> None:
        fetch_into(tmp_path, CORINE)
        legend = importlib.import_module("tin_engine.land_cover_classes").load("corine")
        rows = read_source(
            tmp_path / "clc2018" / "features.gpkg", None, legend.attribute,
            lambda _crs: CORINE.extent(IDS),
        ).rows  # fmt: skip
        assert {int(fid): value for fid, _, value in rows} == {i: CODES[i] for i in IDS}
        for fid, geometry, _ in rows:
            assert shapely.equals(geometry, CORINE.shapes[int(fid)]), fid


# ------------------------------------------------------------ 11. the CLI


@pytest.fixture
def served(monkeypatch: pytest.MonkeyPatch) -> dict[str, FakeLayer]:
    """Both fake layers behind `RangeClient.get_text`, routed by URL."""
    fakes = {NVE_URL: FakeLayer(NVE), CORINE_URL: FakeLayer(CORINE)}

    def get_text(self: RangeClient, url: str) -> str:
        (fake,) = [f for base, f in fakes.items() if url.startswith(base)]
        return fake.get_text(url)

    monkeypatch.setattr(RangeClient, "get_text", get_text)
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)  # no user sources file but the test's
    return fakes


def box_args(data: LayerData) -> list[str]:
    """`--bbox` of the layer's four polygons, in its source's CRS."""
    source = SOURCES["nve-lakes" if data.url == NVE.url else "clc2018"]
    return ["--bbox", *map(str, data.box_all), "--bbox-crs", source.crs]


@needs_rtree
class TestTheCli:
    def test_a_dry_run_prints_the_plan_line_and_writes_nothing(
        self, tmp_path: Path, served: dict[str, FakeLayer]
    ) -> None:
        code, output = invoke(
            "fetch", "clc2018", *box_args(CORINE), "--dry-run", "--cache", str(tmp_path)
        )
        assert code == 0, output
        assert "clc2018: 4 features in the box, 0 present" in output, output
        assert served[CORINE_URL].chunks_asked() == []
        assert not (tmp_path / "clc2018").exists()

    def test_a_fetch_prints_the_report_line(
        self, tmp_path: Path, served: dict[str, FakeLayer]
    ) -> None:
        code, output = invoke("fetch", "nve-lakes", *box_args(NVE), "--cache", str(tmp_path))
        assert code == 0, output
        assert "nve-lakes: 4 features in the box, 0 present, fetched 4 in 3 requests" in output, (
            output
        )
        assert f"into {tmp_path / 'nve-lakes'}" in output, output
        assert (tmp_path / "nve-lakes" / "features.gpkg").is_file()

    def test_sources_lists_the_fetched_entry(
        self, tmp_path: Path, served: dict[str, FakeLayer]
    ) -> None:
        code, output = invoke("fetch", "nve-lakes", *box_args(NVE), "--cache", str(tmp_path))
        assert code == 0, output
        code, text = invoke("sources", "--cache", str(tmp_path))
        assert code == 0, text
        found = {s["id"]: s for s in json.loads(text)["sources"]}
        assert found["nve-lakes"]["present"] is True
        assert found["nve-lakes"]["path"] == str(tmp_path / "nve-lakes" / "features.gpkg")

    def test_a_year_is_a_usage_error(self, tmp_path: Path, served: dict[str, FakeLayer]) -> None:
        code, output = invoke(
            "fetch", "nve-lakes", *box_args(NVE), "--year", "2020", "--cache", str(tmp_path)
        )
        assert code == USAGE, output
        assert "has no years" in output, output
        assert served[NVE_URL].urls == []

    def test_a_refusal_is_exit_1_with_its_reason(
        self, tmp_path: Path, served: dict[str, FakeLayer]
    ) -> None:
        served[NVE_URL].geometry_type = "esriGeometryPolyline"
        code, output = invoke("fetch", "nve-lakes", *box_args(NVE), "--cache", str(tmp_path))
        assert code == REFUSED, output
        assert "only polygon layers are fetched" in output, output
        assert "Traceback" not in output
