"""The catalogue of remote DEM sources `rasputin fetch` copies (increments 23a-1, 23a-2).

`docs/increments/23-basin-scale.md`, "The fetch step and the tile cache" and
B8 as ruled: ANADEM and GLO-30. Data only, importing Pydantic alone, so the
mesh path can name a catalogue key without importing `tin_engine.fetch`
(decided 2): meshing is offline by rule.

Increment 29 adds the catalogue of station lists `rasputin fetch-stations`
copies (`STATION_SOURCES`): NVE's HRD stations, their catchment polygons and
the river lines near them.

Increment 26a adds a source of land-cover classes, MapBiomas, as one COG per
year (`one-cog-per-year`, `years`, `values`; `docs/increments/26-land-cover.md`,
section 3.1): `url_for` picks the year's file, `credit_for` fills the fetch
date and URL into the credit, and `check_values` keeps classes out of `--dem`
and heights out of `--land-cover`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from types import MappingProxyType
from typing import Literal, Protocol, Self

from pydantic import BaseModel, ConfigDict, model_validator


class Dated(Protocol):
    """A `datetime.date`, by its one method used here: this module imports
    only Pydantic (23a-1, decided 2; `test_sources_imports_pydantic_only`)."""

    def isoformat(self) -> str: ...


class RemoteSource(BaseModel):
    """One source: one COG (`url`), one COG per tile (`url_template` and the
    `tile_list_url` saying which tiles exist), or one COG per year
    (`url_template` with `{year}`, `years`). `crs` is the one expected,
    checked against each header; `credit` goes into every mesh made from it."""

    model_config = ConfigDict(frozen=True)

    id: str
    kind: Literal["one-cog", "cog-tiles", "one-cog-per-year"]
    url: str | None = None
    url_template: str | None = None
    tile_list_url: str | None = None
    #: The years a `one-cog-per-year` source has, ascending; the default is the last.
    years: tuple[int, ...] = ()
    #: What the cells hold: heights for a DEM, land-cover class numbers for MapBiomas.
    values: Literal["heights", "classes"] = "heights"
    crs: str
    nodata: float | None
    credit: str
    licence_note: str
    #: The works the distributor asks to cite (23a-2); `notice` and the mesh file carry them.
    cite: tuple[str, ...] = ()

    @model_validator(mode="after")
    def _per_year_has_years(self) -> Self:
        if self.kind == "one-cog-per-year":
            if not self.years:
                raise ValueError(f"{self.id}: a one-cog-per-year source needs its years")
            if "{year}" not in (self.url_template or ""):
                raise ValueError(f"{self.id}: the url_template needs a {{year}} field")
        return self

    def url_for(self, year: int | None = None) -> str:
        """The file to read: `url`, or the year's (the last year when None).
        A year outside `years`, or any year for a source without them, is a
        `ValueError` naming the range."""
        if not self.years:
            if year is not None:
                raise ValueError(f"{self.id} has no years; --year does not apply to it")
            assert self.url is not None, f"{self.id} is not a one-file source"
            return self.url
        year = self.years[-1] if year is None else year
        if year not in self.years:
            raise ValueError(
                f"{self.id} has no year {year}; it covers {self.years[0]} to {self.years[-1]}"
            )
        assert self.url_template is not None
        return self.url_template.format(year=year)

    def credit_for(self, *, url: str, date: Dated) -> str:
        """The credit with `{date}` (ISO) and `{url}` filled in; unchanged without them."""
        return self.credit.format(date=date.isoformat(), url=url)


def check_values(source: RemoteSource, wanted: Literal["heights", "classes"], option: str) -> None:
    """Refuse `source` for `option` unless its cells hold `wanted`."""
    if source.values != wanted:
        raise ValueError(f"{option} takes a source of {wanted}; {source.id} holds {source.values}")


_GLO30 = "https://copernicus-dem-30m.s3.amazonaws.com"

SOURCES: Mapping[str, RemoteSource] = MappingProxyType(
    {
        "anadem-v1": RemoteSource(
            id="anadem-v1",
            kind="one-cog",
            url="https://opentopography.s3.sdsc.edu/raster/ANADEM/ANADEM_be/"
            "anadem_v1_compressed_COG.tif",
            crs="EPSG:4674",
            nodata=-9999.0,
            credit="Agência Nacional de Águas e Saneamento Básico. (2025). ANADEM: A Digital "
            "Terrain Model for South America. Distributed by OpenTopography. "
            "https://doi.org/10.5069/G9736P4G.",
            licence_note="CC BY 4.0 (Creative Commons Attribution 4.0 International)",
            # OpenTopography's acknowledgement asks for the paper as well as the data.
            cite=(
                "Laipelt, L., et al. (2024). ANADEM: A Digital Terrain Model for South America. "
                "Remote Sensing, 16(13), 2321. https://doi.org/10.3390/rs16132321",
            ),
        ),
        "glo30": RemoteSource(
            id="glo30",
            kind="cog-tiles",
            url_template=f"{_GLO30}/{{name}}/{{name}}.tif",
            tile_list_url=f"{_GLO30}/tileList.txt",
            crs="EPSG:4326",
            nodata=None,
            # The licence's Art. 6(b) notice for adapted data, verbatim but "(c)" for the sign.
            credit="produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus "
            "Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union "
            "and ESA; all rights reserved",
            # Art. 6(c): this sentence goes with any distribution of the data, modified or not.
            licence_note="Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "
            '"The organisations in charge of the Copernicus programme by law or by delegation '
            'do not incur any liability for any use of the Copernicus WorldDEM-30"',
        ),
        "mapbiomas-c11": RemoteSource(
            id="mapbiomas-c11",
            kind="one-cog-per-year",
            url_template="https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/"
            "collection11/lulc/coverage/brazil_coverage/brazil_coverage-col11_{year}.tif",
            years=tuple(range(1985, 2026)),
            values="classes",
            crs="EPSG:4326",
            nodata=0.0,
            # MapBiomas's own format (terms of use); the date is the day the file was fetched.
            credit="MapBiomas Project - Collection 11 of the Annual Series of Land Use and Land "
            "Cover Maps of Brazil, accessed on {date} through the link: {url}",
            licence_note="CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 International)",
            cite=(
                "Souza, C. M., et al. (2020). Reconstructing Three Decades of Land Use and Land "
                "Cover Changes in Brazilian Biomes with Landsat Archive and Earth Engine. "
                "Remote Sensing, 12(17), 2735. https://doi.org/10.3390/rs12172735",
            ),
        ),
    }
)


class StationSource(BaseModel):
    """One list of gauging stations: the packaged list `list_file` (in
    `tin_engine/data/`) says which stations; `service_url` serves their points,
    polygons and rivers (29, "The station set")."""

    model_config = ConfigDict(frozen=True)

    id: str
    service_url: str
    list_file: str
    crs: str
    credit: str
    licence_note: str
    cite: tuple[str, ...] = ()


STATION_SOURCES: Mapping[str, StationSource] = MappingProxyType(
    {
        "nve-hrd": StationSource(
            id="nve-hrd",
            service_url="https://kart.nve.no/enterprise/rest/services",
            list_file="nve_hrd_2025.csv",
            crs="EPSG:25833",
            credit="Kilde: NVE. The station list is the 2025 version of NVE's Hydrological "
            "Reference Dataset (Norwegian streamflow reference dataset for climate change "
            "studies); stations, catchments (Totalnedbørfelt til målestasjon) and rivers "
            "(ELVIS) are from NVE's map services.",
            licence_note="Norsk lisens for offentlige data (NLOD). NVE disclaims liability for "
            "errors in the data and their use.",
        ),
    }
)


def notice(source: RemoteSource | StationSource, credits: Sequence[str] = ()) -> str:
    """`<cache>/<source>/NOTICE.txt` (23a-2, decided 8): a rendering of the
    catalogue entry, which stays the one place these words are kept. `credits`
    are the rendered credits of the objects fetched (26a), else the entry's."""
    credit = list(credits) or [source.credit]
    lines = [f"{source.id}", "", "Credit:", *credit, "", "Licence:", source.licence_note]
    if source.cite:
        lines += ["", "Cite:", *source.cite]
    return "\n".join(lines) + "\n"


__all__ = ["SOURCES", "STATION_SOURCES", "RemoteSource", "StationSource", "check_values", "notice"]
