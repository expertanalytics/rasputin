"""`tin_engine.plugin_api`: the names a plugin may import, each the core's own object.

`docs/increments/51-open-core.md`, *4.2 The Python API a plugin may import
(`tin_engine.plugin_api`)* and *9.1 Public (P1 unless marked)*, T7. The list
is the design's table, written out here: a name added to or dropped from
``__all__`` is a visible edit of this file. Identity, not equality: a
re-export that wrapped or copied an object would let the two drift apart.
The modules under test are imported inside each test, so each fails on its own.

``SIGNATURES`` records the shape of each promised name, written out by hand:
a promise change is then a visible edit of this file. A function is recorded
as ``str(inspect.signature(...))`` (parameters, kinds, defaults, and the return
annotation as source text, the modules using postponed annotations). A
Pydantic class is recorded as its fields in order, the name alone when
required, ``name=<repr of default>`` otherwise; a dataclass as its field
names in order; an exception class as its bases; a changed field type is not
caught, only names, order, requiredness and defaults. ``THREADS`` is recorded
as its type, the Typer option's flag, ``min`` and ``metavar``.
"""

from __future__ import annotations

import dataclasses
import inspect
import typing
from importlib import import_module

import pytest

#: Section 4.2's table: home module, then the names re-exported from it.
PROMISED: dict[str, tuple[str, ...]] = {
    "crs": ("frame_factors", "parse_crs", "reprojector", "same_crs", "site_of", "srs_entry"),
    "dem_input": ("DemRequest", "open_dem", "repository_for"),
    "io.models": ("Bounds", "DemTile", "RasterMeta"),
    "io.vtk_read": ("read_mesh",),
    "io.geotiff_write": ("encode_geotiff",),
    "raster": ("core_args",),
    "resources": ("Resources",),
    "run_record": ("plural",),
    "cli_kit": ("destination", "cli_resources", "words", "THREADS"),
    "catchment_core": ("accumulate",),
    # 116: feature layers, the GeoPackage writer and the catalogue
    "io.feature_table": ("FeatureTable", "read_feature_table"),
    "feature_input": ("FeatureError",),
    "io.geopackage_write": ("create_geopackage", "create_features_table", "insert_features"),
    "sources": (
        "RemoteSource",
        "cache_root",
        "catalogue",
        "features_file",
        "kept_named",
        "plugin_source_files",
        "user_source_file",
    ),
}

#: The promised shape of each name in ``__all__``. On a deliberate change,
#: update the entry in its own commit and name the change in the PR body.
SIGNATURES: dict[str, str] = {
    "THREADS": "int | None, --threads, min=0, metavar=N",
    "Bounds": "x_min, y_min, x_max, y_max",
    "DemRequest": (
        "sources=(), cached=None, bounds=None, nodata=None, z_unit=None,"
        " domain=None, target_crs=None, box=None"
    ),
    "DemTile": "meta, array",
    "RasterMeta": (
        "x_min, y_max, delta_x, delta_y, cols, rows, epsg, crs='', geographic=False,"
        " nodata, nodata_source, pixel_is_area, vertical_unit_assumed,"
        " z_unit='metre', file_crs=None"
    ),
    "Resources": "threads, memory_budget",
    "accumulate": "(tile: 'DemTile') -> 'AccumulateOutcome'",
    "cli_resources": ("(threads: 'int | None', memory_budget: 'str | None') -> 'Resources'"),
    "core_args": "(tile: 'DemTile') -> 'RasterArgs'",
    "destination": "(out: 'Path', hint: 'str' = '--out') -> 'Path'",
    "encode_geotiff": (
        "(grid: 'npt.ArrayLike', x0: 'float', y0: 'float', step: 'float', epsg: 'int') -> 'bytes'"
    ),
    "frame_factors": ("(crs: 'str | CRS', lon: 'Any', lat: 'Any') -> 'tuple[Xy, Xy]'"),
    "open_dem": ("(request: 'DemRequest', resources: 'Resources | None' = None) -> 'DemInput'"),
    "parse_crs": "(text: 'str | CRS') -> 'CRS'",
    "plural": "(n: 'int', one: 'str', many: 'str') -> 'str'",
    "read_mesh": "(data: 'bytes') -> 'MeshArrays'",
    "repository_for": (
        "(sources: 'tuple[Path, ...]', nodata: 'float | None' = None,"
        " cached: 'CachedSource | None' = None, threads: 'int' = 1,"
        " z_unit: 'ZUnit | None' = None)"
        " -> 'TiffDemRepository | CacheRepository'"
    ),
    "reprojector": (
        "(src: 'str | CRS', dst: 'str | CRS', *, threads: 'int' = 1) -> 'Callable[[Any], Xy]'"
    ),
    "same_crs": "(a: 'str | CRS', b: 'str | CRS') -> 'bool'",
    "site_of": "(crs: 'str | CRS', x: 'float', y: 'float') -> 'Site'",
    "words": "(exc: 'ValueError') -> 'str'",
    "FeatureError": "(ValueError)",
    "FeatureTable": "crs, layer, columns, rows",
    "RemoteSource": (
        "id, kind, path=None, url=None, url_template=None, tile_list_url=None,"
        " tile_name=None, tile_spacing=None, years=(), values='heights', fields=(),"
        " layer=None, crs=None, nodata=None, credit, licence_note, cite=(), classes='',"
        " z_unit=None, description='', named_by=None, position_error_m=None"
    ),
    "cache_root": "(option: 'Path | None', environ: 'Mapping[str, str]') -> 'Path | None'",
    "catalogue": "(user_file: 'Path | None', plugins: 'Sequence[PluginFile]' = ()) -> 'Catalogue'",
    "create_features_table": (
        "(conn: 'sqlite3.Connection', name: 'str', srs_id: 'int', kind: 'str', pk: 'str',"
        " columns: 'Sequence[tuple[str, str]]') -> 'None'"
    ),
    "create_geopackage": "(conn: 'sqlite3.Connection', srs: 'SrsEntry') -> 'None'",
    "features_file": "(source: 'RemoteSource', cache: 'Path | None') -> 'Path'",
    "insert_features": (
        "(conn: 'sqlite3.Connection', name: 'str', srs_id: 'int',"
        " rows: 'Sequence[tuple[int, BaseGeometry, Sequence[Any]]]') -> 'None'"
    ),
    "kept_named": ("(keys: 'Sequence[Any]', names: 'Sequence[Any]', prefix: 'str') -> 'list[int]'"),
    "plugin_source_files": "() -> 'tuple[PluginFile, ...]'",
    "read_feature_table": (
        "(path: 'Path', layer: 'str | None',"
        " box_for: 'Callable[[str], tuple[float, float, float, float]] | None' = None)"
        " -> 'FeatureTable'"
    ),
    "srs_entry": "(crs: 'str | None') -> 'SrsEntry'",
    "user_source_file": "(environ: 'Mapping[str, str]') -> 'Path | None'",
}

NAMES = [(home, name) for home, names in PROMISED.items() for name in names]


def test_all_is_exactly_the_promised_list() -> None:
    api = import_module("tin_engine.plugin_api")
    assert sorted(api.__all__) == sorted(name for _, name in NAMES)
    assert len(api.__all__) == len(set(api.__all__)), "a name listed twice"


@pytest.mark.parametrize(("home", "name"), NAMES, ids=[f"{h}.{n}" for h, n in NAMES])
def test_each_name_is_the_object_in_its_home_module(home: str, name: str) -> None:
    api = import_module("tin_engine.plugin_api")
    assert getattr(api, name) is getattr(import_module(f"tin_engine.{home}"), name)


def test_cli_shares_the_threads_option_with_cli_kit() -> None:
    cli, kit = import_module("tin_engine.cli"), import_module("tin_engine.cli_kit")
    assert cli.THREADS is kit.THREADS


def _shape(name: str, obj: object) -> str:
    """The text ``SIGNATURES`` records for ``obj``, by the rules of the module docstring."""
    if name == "THREADS":
        kind, option = typing.get_args(obj)[:2]
        return f"{kind}, {option.default}, min={option.min}, metavar={option.metavar}"
    if inspect.isclass(obj) and issubclass(obj, BaseException):
        return f"({', '.join(base.__name__ for base in obj.__bases__)})"
    if dataclasses.is_dataclass(obj):
        return ", ".join(field.name for field in dataclasses.fields(obj))
    if inspect.isclass(obj):
        fields = obj.model_fields.items()  # type: ignore[attr-defined]
        return ", ".join(k if f.is_required() else f"{k}={f.default!r}" for k, f in fields)
    return str(inspect.signature(obj))  # type: ignore[arg-type]


def test_signatures_record_exactly_the_promised_list() -> None:
    api = import_module("tin_engine.plugin_api")
    assert sorted(SIGNATURES) == sorted(api.__all__)


@pytest.mark.parametrize("name", sorted(SIGNATURES))
def test_each_promised_name_keeps_its_recorded_signature(name: str) -> None:
    api = import_module("tin_engine.plugin_api")
    current = _shape(name, getattr(api, name))
    assert current == SIGNATURES[name], (
        f"a promised signature changed: plugins (rasputin-lsm among them) break."
        f" If on purpose, update SIGNATURES in its own commit and say so in the"
        f" PR body.\n  {name} recorded: {SIGNATURES[name]}\n  {name} current:  {current}"
    )
