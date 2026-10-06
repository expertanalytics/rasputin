"""The hydrography's value types: river segments, gauging stations, lakes.

Layer 0, no first-party import (python-audit.md, section 11): `gauge` places
on them without importing a codec, and `io/rivers.py` and
`io/station_set.py` read files into them. The fields are the ones increment
29 defines (`docs/increments/29-nve-reference-catchments.md`, "Placing the
gauge", "The station set" and "Lake gauges").
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field
from shapely.geometry import Polygon


class RiverSegment(BaseModel):
    """One mapped line, digitised downstream, in its file's CRS. `objekttype`
    is kept as served; `kind` is what "Placing the gauge" reads."""

    model_config = ConfigDict(frozen=True)

    objectid: int
    elvid: str | None
    vassdragsnr: str | None
    name: str | None
    objekttype: str | None
    kind: Literal["lake", "river"]
    line: tuple[tuple[float, float], ...]


class Station(BaseModel):
    """One gauging station, in its file's CRS. `series` holds the discharge
    series keys (`"1001.0"`), the join to NVE's discharge later."""

    model_config = ConfigDict(frozen=True)

    station: str = Field(pattern=r"^\d+\.\d+\.\d+$")
    name: str | None = None
    x: float
    y: float
    series: tuple[str, ...] = ()
    nve_area_km2: float | None = None
    watercourse: str | None = None
    river: str | None = None


@dataclass(frozen=True, slots=True)
class Lake:
    """One lake polygon part, in its file's CRS (PR 4, "Lake gauges").
    `number` is NVE's `vatnlnr`, None for null or 0; `name` is `navn`."""

    number: int | None
    name: str | None
    polygon: Polygon


__all__ = ["Lake", "RiverSegment", "Station"]
