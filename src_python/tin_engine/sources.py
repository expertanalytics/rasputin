"""The catalogue of remote DEM sources `rasputin fetch` copies (increments 23a-1, 23a-2).

`docs/increments/23-basin-scale.md`, "The fetch step and the tile cache" and
B8 as ruled: ANADEM and GLO-30. Data only, importing Pydantic alone, so the
mesh path can name a catalogue key without importing `tin_engine.fetch`
(decided 2): meshing is offline by rule.

Increment 29 adds the catalogue of station lists `rasputin fetch-stations`
copies (`STATION_SOURCES`): NVE's HRD stations, their catchment polygons and
the river lines near them.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Literal

from pydantic import BaseModel, ConfigDict


class RemoteSource(BaseModel):
    """One source: one COG (`url`), or one COG per tile (`url_template` and
    the `tile_list_url` saying which tiles exist). `crs` is the one expected,
    checked against each header; `credit` goes into every mesh made from it."""

    model_config = ConfigDict(frozen=True)

    id: str
    kind: Literal["one-cog", "cog-tiles"]
    url: str | None = None
    url_template: str | None = None
    tile_list_url: str | None = None
    crs: str
    nodata: float | None
    credit: str
    licence_note: str
    #: The works the distributor asks to cite (23a-2); `notice` and the mesh file carry them.
    cite: tuple[str, ...] = ()


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


def notice(source: RemoteSource | StationSource) -> str:
    """`<cache>/<source>/NOTICE.txt` (23a-2, decided 8): a rendering of the
    catalogue entry, which stays the one place these words are kept."""
    lines = [f"{source.id}", "", "Credit:", source.credit, "", "Licence:", source.licence_note]
    if source.cite:
        lines += ["", "Cite:", *source.cite]
    return "\n".join(lines) + "\n"


__all__ = ["SOURCES", "STATION_SOURCES", "RemoteSource", "StationSource", "notice"]
