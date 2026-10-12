"""The MapBiomas catalogue entry and the refusals by `values` (26a, `tin_engine.sources`).

`docs/increments/26-land-cover.md`, section 3.1 ("The catalogue entry"),
section 3.2 and Ola's rulings of 2026-10-04: Collection 11, year 2025 by
default (Q1); CC BY-SA 4.0 in the licence note (Q2 (a)). Pure: no file, no
network. The pinned spellings are in `class_fixtures`' docstring.

HOW THIS FILE GOES RED: `SOURCES` has no `mapbiomas-c11` (`KeyError`),
`RemoteSource` has no `years`, `values`, `url_for` or `credit_for`, and
`sources` has no `check_values` (`ValidationError`, `AttributeError`); each
test fails on its own and the rest of `tests/python` collects.
"""

from __future__ import annotations

import datetime
import importlib
from types import ModuleType
from typing import Any

import pytest

from class_fixtures import REAL_KEY, per_year_entry

URL = (
    "https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection11/"
    "lulc/coverage/brazil_coverage/brazil_coverage-col11_{year}.tif"
)
DOI = "10.3390/rs12172735"
DATE = datetime.date(2026, 10, 11)


@pytest.fixture(scope="module")
def sources() -> ModuleType:
    return importlib.import_module("tin_engine.sources")


@pytest.fixture
def entry(sources: ModuleType) -> Any:
    return sources.SOURCES[REAL_KEY]


class TestTheCatalogueEntry:
    def test_one_cog_per_year_of_classes_in_wgs84_with_0_as_no_data(self, entry: Any) -> None:
        assert entry.id == REAL_KEY
        assert entry.kind == "one-cog-per-year"
        assert entry.values == "classes"
        assert entry.years == tuple(range(1985, 2026))
        assert entry.url is None and entry.url_template == URL
        assert entry.crs == "EPSG:4326"
        assert entry.nodata == 0

    @pytest.mark.parametrize("year", [2025, 1985])
    def test_the_url_for_a_year(self, entry: Any, year: int) -> None:
        assert entry.url_for(year) == URL.format(year=year)

    def test_the_default_year_is_the_last(self, entry: Any) -> None:
        assert entry.url_for() == URL.format(year=2025)

    @pytest.mark.parametrize("year", [1984, 2026])
    def test_a_year_outside_the_range_is_refused_naming_the_range(
        self, entry: Any, year: int
    ) -> None:
        with pytest.raises(ValueError, match=r"1985.*2025") as caught:
            entry.url_for(year)
        assert str(year) in str(caught.value)

    def test_the_credit_licence_and_citation(self, entry: Any) -> None:
        assert "Collection 11" in entry.credit
        assert "{date}" in entry.credit and "{url}" in entry.credit
        assert "CC BY-SA 4.0" in entry.licence_note
        assert len(entry.cite) == 1
        assert DOI in entry.cite[0] and "Souza" in entry.cite[0]

    def test_the_credit_rendered_with_the_fetch_date_and_url(self, entry: Any) -> None:
        url = URL.format(year=2025)
        rendered = entry.credit_for(url=url, date=DATE)
        assert rendered == entry.credit.format(date="2026-10-11", url=url)
        assert f"accessed on 2026-10-11 through the link: {url}" in rendered
        assert "{" not in rendered

    @pytest.mark.parametrize("key", ["anadem-v1", "glo30"])
    def test_the_dem_sources_hold_heights_and_have_no_years(
        self, sources: ModuleType, key: str
    ) -> None:
        dem = sources.SOURCES[key]
        assert dem.values == "heights"
        assert dem.years == ()
        assert dem.credit_for(url="https://example.invalid/x.tif", date=DATE) == dem.credit

    def test_a_per_year_entry_needs_years_and_a_year_in_its_template(
        self, sources: ModuleType
    ) -> None:
        with pytest.raises(ValueError, match="years"):
            per_year_entry(sources, URL, years=())
        with pytest.raises(ValueError, match=r"\{year\}"):
            per_year_entry(sources, "https://example.invalid/one.tif")


class TestRefusalsByValues:
    """`--dem` takes heights only; `--land-cover` and `water --source` (26b,
    26d) take classes only; one check for all of them, named by the option."""

    def test_a_class_source_is_refused_as_heights(self, sources: ModuleType, entry: Any) -> None:
        with pytest.raises(ValueError) as caught:
            sources.check_values(entry, "heights", "--dem")
        message = str(caught.value)
        assert REAL_KEY in message and "classes" in message and "--dem" in message

    def test_a_heights_source_is_refused_as_classes(self, sources: ModuleType) -> None:
        with pytest.raises(ValueError) as caught:
            sources.check_values(sources.SOURCES["anadem-v1"], "classes", "--land-cover")
        message = str(caught.value)
        assert "anadem-v1" in message and "heights" in message and "--land-cover" in message

    def test_a_matching_source_passes(self, sources: ModuleType, entry: Any) -> None:
        assert sources.check_values(entry, "classes", "--land-cover") is None
        assert sources.check_values(sources.SOURCES["anadem-v1"], "heights", "--dem") is None
