"""The catalogue of remote DEM sources `rasputin fetch` copies (increment 23a-1).

`docs/increments/23-basin-scale.md`, "The fetch step and the tile cache" and
B8 as ruled: ANADEM and GLO-30. Data only, importing Pydantic alone, so the
mesh path can name a catalogue key without importing `tin_engine.fetch`
(decided 2): meshing is offline by rule.
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
            credit="ANADEM v1, distributed by OpenTopography, doi:10.5069/G9736P4G",
            licence_note="read the licence on the DOI's landing page before redistributing",
        ),
        "glo30": RemoteSource(
            id="glo30",
            kind="cog-tiles",
            url_template=f"{_GLO30}/{{name}}/{{name}}.tif",
            tile_list_url=f"{_GLO30}/tileList.txt",
            crs="EPSG:4326",
            nodata=None,
            credit="Copernicus DEM GLO-30, (c) DLR e.V. 2010-2014 and (c) Airbus Defence "
            "and Space GmbH 2014-2018, provided under COPERNICUS by the European Union and ESA",
            licence_note="Copernicus DEM licence: free use with this credit",
        ),
    }
)

__all__ = ["SOURCES", "RemoteSource"]
