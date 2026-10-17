"""The catalogue of sources: remote ones `rasputin fetch` and `rasputin fetch-stations`
copy, and local ones read in place.

The entries are data: `tin_engine/data/sources.toml` ships them, read here at
import, and `catalogue` adds the user's own from their sources file, which
`user_source_file` finds (local datasets are registered there).
`SOURCES` holds the rasters `rasputin fetch` copies into the tile cache: two
DEMs (ANADEM and GLO-30) and one land-cover raster (MapBiomas, one COG per
year). This module imports no NumPy, pyproj or `tin_engine.fetch`, so the mesh
path can name a catalogue key without them: meshing runs offline. `url_for`
picks a file (a year's, for MapBiomas), `credit_for` fills the fetch date and
URL into the credit, and `check_values` keeps classes out of `--dem` and
heights out of `--land-cover`.

A local entry may also be a vector file (`.gpkg`, `.geojson`, `.json`): it
holds features, which the vector options take by path. A plugin distribution
adds its own `sources.toml` through the entry point group `rasputin.sources`
(`plugin_source_files`).

`STATION_SOURCES` holds the station lists `rasputin fetch-stations` copies:
NVE's HRD stations, their catchment polygons and the river lines near them.
"""

from __future__ import annotations

import os
import string
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib import resources
from importlib.metadata import entry_points
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal, Protocol, Self

from pydantic import BaseModel, ConfigDict, PositiveFloat, model_validator

from tin_engine.toml_data import parse, validate

ZUnit = Literal["metre", "foot", "US survey foot"]
#: The fields each kind's `url_template` takes, and needs.
_URL_FIELDS = {"one-cog": set(), "cog-tiles": {"name"}, "one-cog-per-year": {"year"}}
#: The suffixes of the files the lake and feature readers both take.
VECTOR_SUFFIXES = (".gpkg", ".geojson", ".json")
#: The entry point group whose value names a package holding a `sources.toml`.
SOURCE_GROUP = "rasputin.sources"
#: The MCP server's output root, set in its children's environment; the listing
#: reports it and nothing is fenced by it.
OUT_ROOT_VARIABLE = "RASPUTIN_MCP_OUT_ROOT"
_REMOTE_FIELDS = ("url", "url_template", "tile_list_url", "tile_name", "tile_spacing", "years")
#: The fields an `arcgis-features` entry refuses: its `url` is the layer, its `fields` the columns.
_LAYER_REFUSED = (*_REMOTE_FIELDS[1:], "path", "layer", "nodata", "z_unit")


class Dated(Protocol):
    """A `datetime.date`, by its one method used here: this module imports no
    third-party package but Pydantic."""

    def isoformat(self) -> str: ...


def _fields(name: str, template: str | None, allowed: set[str], *, needed: bool = False) -> None:
    """Refuse a format string using a field outside `allowed` (a dotted or
    indexed one too, so a file cannot reach an object's attributes), or,
    when `needed`, lacking one of them."""
    used = {f for _, f, _, _ in string.Formatter().parse(template or "") if f is not None}
    if used - allowed or (needed and allowed - used):
        shown = ", ".join(f"{{{f}}}" for f in sorted(used)) or "no field"
        takes = ", ".join(f"{{{f}}}" for f in sorted(allowed)) or "none"
        raise ValueError(f"{name} uses {shown}; it takes {takes}" + (", each" if needed else ""))


class RemoteSource(BaseModel):
    """One source: one COG (`url`), one COG per tile (`url_template` and the
    `tile_list_url` saying which tiles exist), one COG per year
    (`url_template` with `{year}`, `years`), or a local GeoTIFF or directory
    of them (`path`, read in place; its CRS is the files' own), or a local
    vector file, which holds features (`layer`: the GeoPackage's), or the
    polygons of an ArcGIS map-service layer (`url`, `fields`), fetched by box
    into a GeoPackage in the cache. A remote
    source's `crs` is the one expected, checked against each header; `credit`
    goes into every mesh made from it.
    An unknown key is refused, so a misspelt one in a file is not ignored."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    id: str
    kind: Literal["local", "one-cog", "cog-tiles", "one-cog-per-year", "arcgis-features"]
    #: A local source's GeoTIFF, or directory of GeoTIFFs.
    path: Path | None = None
    url: str | None = None
    url_template: str | None = None
    tile_list_url: str | None = None
    #: The object name of the 1° tile whose south-west corner is `{lat}`, `{lon}`
    #: (absolute degrees), `{ns}` N or S and `{ew}` E or W; tiled sources only.
    tile_name: str | None = None
    #: The largest (lon, lat) cell of any tile, in degrees; tiled sources only.
    tile_spacing: tuple[PositiveFloat, PositiveFloat] | None = None
    #: The years a `one-cog-per-year` source has, ascending; the default is the last.
    years: tuple[int, ...] = ()
    #: What the cells hold: heights for a DEM, land-cover class numbers for
    #: MapBiomas; a vector file holds features.
    values: Literal["heights", "classes", "features"] = "heights"
    #: The attributes an `arcgis-features` layer is asked for and kept; its id field is the key.
    fields: tuple[str, ...] = ()
    #: The layer of a GeoPackage with more than one features table.
    layer: str | None = None
    crs: str | None = None
    #: The no-data value the source states; None when it states none.
    nodata: float | None = None
    credit: str
    licence_note: str
    #: The works the distributor asks to cite; `notice` and the mesh file carry them.
    cite: tuple[str, ...] = ()
    #: The legend (`tin_engine.legends`) that names and colours its classes.
    legend: str = ""
    #: The height unit the catalogue states for files that do not.
    z_unit: ZUnit | None = None
    #: What the cells represent, and the resolution and region, as the distributor states them.
    description: str = ""

    @model_validator(mode="before")
    @classmethod
    def _features_by_suffix(cls, data: Any) -> Any:
        """A vector path or a map-service layer stating no `values` holds
        features: the file or the kind says what it is."""
        if not isinstance(data, dict) or "values" in data:
            return data
        if _vector(data.get("path")) or data.get("kind") == "arcgis-features":
            return {**data, "values": "features"}
        return data

    @model_validator(mode="after")
    def _consistent(self) -> Self:
        if self.fields and self.kind != "arcgis-features":
            raise ValueError(f"fields: refused for {self.kind}")
        if self.kind == "arcgis-features":
            return self._layer_consistent()
        suffix = "" if self.path is None else self.path.suffix.lower()
        if self.layer is not None and suffix != ".gpkg":
            raise ValueError("layer: only a .gpkg path has layers")
        if self.values == "features":
            if self.kind != "local":
                raise ValueError(f"values: features need a local source, not {self.kind}")
            if not _vector(self.path):
                raise ValueError("values: features need a .gpkg, .geojson or .json path")
            if self.z_unit is not None:
                raise ValueError("z_unit: features hold no heights")
        elif _vector(self.path):
            raise ValueError(f"values: a {suffix} path holds features")
        if self.kind == "local":
            for field in (*_REMOTE_FIELDS, "crs"):
                if getattr(self, field) not in (None, ()):
                    raise ValueError(f"{field}: refused for a local source, read in place")
            if self.path is None:
                raise ValueError("path: needed for a local source")
            if self.values == "classes":
                raise ValueError("values: a local class raster is given to --land-cover as a file")
            return self
        for field, needed in (("path", False), ("crs", True)):
            if needed != (getattr(self, field) is not None):
                raise ValueError(f"{field}: {'needed' if needed else 'refused'} for {self.kind}")
        if self.kind == "one-cog-per-year" and not self.years:
            raise ValueError("years: a one-cog-per-year source needs its years")
        tiled = self.kind == "cog-tiles"
        for field in ("tile_list_url", "tile_name", "tile_spacing"):
            if tiled != (getattr(self, field) is not None):
                raise ValueError(f"{field}: {'needed' if tiled else 'refused'} for {self.kind}")
        wanted = _URL_FIELDS[self.kind]
        _fields("url_template", self.url_template, wanted, needed=bool(wanted))
        _fields("tile_name", self.tile_name, {"ns", "ew", "lat", "lon"})
        _fields("credit", self.credit, {"date", "url"})
        if self.values == "classes" and self.z_unit is not None:
            raise ValueError("z_unit: class numbers are not heights")
        return self

    def _layer_consistent(self) -> Self:
        for field in ("url", "crs", "fields"):
            if not getattr(self, field):
                raise ValueError(f"{field}: needed for {self.kind}")
        if len(set(self.fields)) != len(self.fields):
            raise ValueError("fields: a field is named twice")
        for field in _LAYER_REFUSED:
            if getattr(self, field) not in (None, ()):
                raise ValueError(f"{field}: refused for {self.kind}")
        if self.values != "features":
            raise ValueError(f"values: {self.kind} holds features")
        _fields("credit", self.credit, {"date", "url"})
        return self

    @property
    def read_unit(self) -> ZUnit | None:
        """The height unit every reader of the source's files is given:
        metres for classes, which hold no heights, else `z_unit`."""
        return "metre" if self.values == "classes" else self.z_unit

    def tile_for(self, lat: int, lon: int) -> str:
        """The name of the tile whose south-west corner is (`lat`, `lon`)."""
        assert self.tile_name is not None, f"{self.id} is not a tiled source"
        ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
        return self.tile_name.format(ns=ns, ew=ew, lat=abs(lat), lon=abs(lon))

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


def _vector(path: object) -> bool:
    """Whether `path` names a vector file, by its suffix."""
    return path is not None and Path(str(path)).suffix.lower() in VECTOR_SUFFIXES


def check_values(source: RemoteSource, wanted: Literal["heights", "classes"], option: str) -> None:
    """Refuse `source` for `option` unless its cells hold `wanted`."""
    if source.values != wanted:
        raise ValueError(f"{option} takes a source of {wanted}; {source.id} holds {source.values}")


class StationSource(BaseModel):
    """One list of gauging stations: the packaged list `list_file` (in
    `tin_engine/data/`) says which stations; `service_url` serves their points,
    polygons and rivers."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    id: str
    service_url: str
    list_file: str
    crs: str
    credit: str
    licence_note: str
    cite: tuple[str, ...] = ()


class _Tables(BaseModel):
    """A catalogue file: `[source.<id>]` and `[station.<id>]` tables."""

    model_config = ConfigDict(extra="forbid")

    source: dict[str, dict[str, Any]] = {}
    station: dict[str, dict[str, Any]] = {}


@dataclass(frozen=True, slots=True)
class Catalogue:
    """The shipped entries, the user's and the plugins', by id; where each
    came from ("package", "user", "plugin <distribution>"), and one line per
    plugin file skipped."""

    sources: Mapping[str, RemoteSource]
    stations: Mapping[str, StationSource]
    origins: Mapping[str, str]
    refused: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class PluginFile:
    """One entry point of `SOURCE_GROUP`: its distribution, and its
    `sources.toml` or why there is none."""

    dist: str
    path: Path | None
    error: str | None = None


def _entries[M: BaseModel](
    tables: dict[str, dict[str, Any]], kind: str, model: type[M], label: str
) -> dict[str, M]:
    """Each table validated as `model`, its key the id; an `id` key is refused."""
    out: dict[str, M] = {}
    for key, table in tables.items():
        if "id" in table:
            raise ValueError(f"{label}: {kind}.{key}.id: the table's key is the id; drop the key")
        out[key] = validate({**table, "id": key}, model, label, (kind, key))
    return out


_SHIPPED = "sources.toml"
_tables = parse(
    (resources.files("tin_engine") / "data" / _SHIPPED).read_text(encoding="utf-8"),
    _Tables,
    _SHIPPED,
)
SOURCES: Mapping[str, RemoteSource] = MappingProxyType(
    _entries(_tables.source, "source", RemoteSource, _SHIPPED)
)
STATION_SOURCES: Mapping[str, StationSource] = MappingProxyType(
    _entries(_tables.station, "station", StationSource, _SHIPPED)
)


def data_root(environ: Mapping[str, str]) -> Path | None:
    """The data root, ``$RASPUTIN_DATA``, or None when it is unset or empty:
    the one reader of the variable."""
    data = environ.get("RASPUTIN_DATA")
    return Path(data) if data else None


def xdg_base(environ: Mapping[str, str], variable: str, under_home: str) -> Path | None:
    """The XDG Base Directory rule: ``$variable``, else ``$HOME/under_home``;
    an empty or relative value is ignored, and None when neither is usable."""
    for value in (environ.get(variable), f"{environ.get('HOME') or '.'}/{under_home}"):
        if (path := Path(value or ".")).is_absolute():
            return path
    return None


def user_source_candidates(environ: Mapping[str, str]) -> tuple[Path, ...]:
    """The places the user's sources file is looked for, in order:
    ``$RASPUTIN_DATA/sources.toml``, then ``rasputin/sources.toml`` in the
    XDG configuration directory. Pure in `environ`."""
    root, config = data_root(environ), xdg_base(environ, "XDG_CONFIG_HOME", ".config")
    places = (root, config / "rasputin" if config else None)
    return tuple(place / "sources.toml" for place in places if place is not None)


def user_source_file(environ: Mapping[str, str]) -> Path | None:
    """The first candidate that is a file: one file only, never merged."""
    return next((p for p in user_source_candidates(environ) if p.is_file()), None)


def plugin_source_files() -> tuple[PluginFile, ...]:
    """The `sources.toml` of each package named in `SOURCE_GROUP`, in
    distribution-name order. Never raises: a value naming an attribute, a
    package that fails to import, or one with no `sources.toml` is recorded
    with the reason."""
    found = []
    for ep in entry_points(group=SOURCE_GROUP):
        dist = ep.dist.name if ep.dist else ep.name
        try:
            if ep.attr:
                raise ValueError(f"{ep.value} names an attribute; name the package")
            path = Path(str(resources.files(ep.module) / _SHIPPED))
            if not path.is_file():
                raise ValueError(f"{ep.module} holds no {_SHIPPED}")
        except Exception as exc:  # a broken plugin is reported, never raised
            found.append(PluginFile(dist, None, f"{type(exc).__name__}: {exc}"))
            continue
        found.append(PluginFile(dist, path))
    return tuple(sorted(found, key=lambda p: p.dist))


def _beside(file: Path, value: str) -> Path:
    """A path a catalogue file names: `~` expanded, a relative one from the file's directory."""
    return file.parent / Path(value).expanduser()


def _read(file: Path) -> dict[str, RemoteSource]:
    """The `[source.<id>]` tables of a catalogue file; a relative `legend`
    or `path` is taken from the file's directory. Refused, as a `ValueError`
    naming the file: anything `RemoteSource` refuses, and a station table
    (fetching a station list follows one service's layout, which is code)."""
    label = str(file)
    try:
        tables = parse(file.read_text(encoding="utf-8"), _Tables, label)
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"{label} cannot be read: {exc}") from exc
    for key in tables.station:
        raise ValueError(f"{label}: station.{key}: station lists are not added from a file")
    for table in tables.source.values():
        given = str(table.get("legend", ""))
        if given.endswith(".toml") or "/" in given or os.sep in given:
            table["legend"] = str(_beside(file, given))
        if "path" in table:
            table["path"] = str(_beside(file, str(table["path"])))
    return _entries(tables.source, "source", RemoteSource, label)


def catalogue(user_file: Path | None, plugins: Sequence[PluginFile] = ()) -> Catalogue:
    """The shipped entries, the `[source.<id>]` tables of `user_file` if it
    exists, then each plugin's. Whether a `path` exists is the reader's to
    find out. The user's file is refused whole, as a `ValueError` naming it,
    when `_read` refuses it or an id is the package's. A plugin's file is
    skipped whole, as one `refused` line, when it cannot be read or one of
    its ids is taken: a plugin never replaces an entry nor stops the core."""
    found, origins, refused = dict(SOURCES), dict.fromkeys(SOURCES, "package"), []
    if user_file is not None and user_file.is_file():
        for key, entry in _read(user_file).items():
            if key in found:
                raise ValueError(
                    f"{user_file}: source.{key}: {key} is a source of the package's "
                    f"{_SHIPPED}; give yours another id"
                )
            found[key], origins[key] = entry, "user"
    for plugin in plugins:
        try:
            if plugin.path is None:
                raise ValueError(plugin.error)
            entries = _read(plugin.path)
            if taken := sorted(entries.keys() & found.keys()):
                raise ValueError(f"source {taken[0]} is taken by the {origins[taken[0]]}")
        except ValueError as exc:
            refused.append(f"plugin {plugin.dist}: its sources skipped: {exc}")
            continue
        found |= entries
        origins |= dict.fromkeys(entries, f"plugin {plugin.dist}")
    return Catalogue(
        MappingProxyType(found), STATION_SOURCES, MappingProxyType(origins), tuple(refused)
    )


def notice(source: RemoteSource | StationSource, credits: Sequence[str] = ()) -> str:
    """The text of a fetch's `NOTICE.txt`: a rendering of the catalogue
    entry, which stays the one place these words are kept. `credits` are the
    rendered credits of the objects fetched (each with its own fetch date),
    else the entry's."""
    credit = list(credits) or [source.credit]
    lines = [f"{source.id}", "", "Credit:", *credit, "", "Licence:", source.licence_note]
    if source.cite:
        lines += ["", "Cite:", *source.cite]
    return "\n".join(lines) + "\n"


__all__ = [
    "OUT_ROOT_VARIABLE",
    "SOURCES",
    "SOURCE_GROUP",
    "STATION_SOURCES",
    "VECTOR_SUFFIXES",
    "Catalogue",
    "PluginFile",
    "RemoteSource",
    "StationSource",
    "catalogue",
    "check_values",
    "data_root",
    "notice",
    "plugin_source_files",
    "user_source_candidates",
    "user_source_file",
    "xdg_base",
]
