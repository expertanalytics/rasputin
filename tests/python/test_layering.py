"""The dependency map of `tin_engine`: which first-party modules each may import.

`docs/increments/python-audit.md`, sections 5 and 8. Six layers, bottom up:
L0 values, L1 pure algorithms, L2 `io/` codecs, L3 the adapters that alone call
`_core`, L4 pipelines, L5 `cli`. A module imports from its own layer or below,
and only L3 imports `_core`; `UPWARD` names the edges that break this today,
each with the PR that removes it.

Each row is the exact set a module imports, read by `importscan`
(`TYPE_CHECKING` imports count; prose does not). Equality, not a subset: a new
edge is a visible table edit, in a `@tester` commit of the PR that adds it. The
rows replace the per-module firewall tests eight suites used to carry.
"""

from __future__ import annotations

import ast
import inspect
from importlib import import_module
from pathlib import Path

import pytest

import tin_engine
from importscan import first_party_imports

# Grouped by layer, so the map reads as a map. Rows that import nothing share a
# line, and long rows continue under their opening quote; the formatter would
# give each row a line of its own and flush the continuations left (or join a
# two-piece row onto one line where the joined line still fits its limit).
# fmt: off
LAYERS: tuple[dict[str, str], ...] = (
    {  # L0: values
        "io.models": "", "features": "", "run_record": "", "stats": "run_record",
        "sources": "toml_data",  # 44: the catalogue read from TOML
        "palettes": "", "hydrography": "", "tin_engine": "",
        "toml_data": "",  # 44: the shared TOML reader; tomllib and Pydantic only
    },
    {  # L1: pure algorithms
        "crs": "io.models run_record",
        "mosaic": "crs io.models resources",  # 114: a printed box names its CRS
        "target_grid": "crs domain grow io.models",
        "grid_domain": "io.models",
        "domain": "crs io.models",
        "chains": "domain features",
        "elevation": "", "outline": "", "sensitivity": "", "landcover": "",
        "grow": "",  # 30d: the outline grown, shapely and NumPy only
        "resources": "",  # 38: the one module that reads the machine; os and Pydantic only
        "burn": "gauge geodesy io.models",  # 40: the burn in a metre frame
        "geodesy": "crs io.models",  # 40: metres near a point, geodesic areas
        "gauge": "hydrography",
        "reference": "io.models",
        "decompose": "",  # 23c: the partition rule; reads no machine, imports nothing
        "budget": "",  # the triangle budget's search; pure, imports nothing
        "network": "geodesy",  # 85: a catchment joined from a sub-catchment network
    },
    {  # L2: io/ codecs, and the one module that opens a connection
        "io": "io.ply io.vtk_legacy",
        "io.cog": "io.geotiff io.models",
        "io.domain_file": "crs domain io.geojson io.repository",
        "io.geojson": "crs",
        "io.geopackage": "", "io.gml": "", "io.mesh_checks": "",
        "io.geopackage_write": "features io.geopackage io.models",  # 36: the writer
        "io.geotiff": "crs io.models",
        "io.mesh_index": "io.models",
        "io.ply": "features io.mesh_checks io.models",  # 26b: the ClassFractions table
        # 63: the features cache, a GeoPackage
        "io.repository": "crs io.cog io.geotiff io.geopackage_write io.models mosaic",  # 114
        "io.rivers": "hydrography io.station_set",
        "io.station_set": "hydrography io.geojson io.repository",
        "io.vtk_legacy": "features io.mesh_checks io.models",  # 26b: the ClassFractions table
        "fetch.http": "tin_engine",
        "io.geotiff_write": "", "io.vtk_read": "",  # 51: in plugin_api, for plugins
    },
    {  # L3: the only importers of _core
        "raster": "_core io.models",
        "edge_strip": "_core stats",
        "final_check": "_core stats target_grid",
        "catchment_core": "_core io.models raster",
        "border_simplify": "_core",  # 32: land-cover borders within a band
        "pieces": "_core io.mesh_index landcover",  # 23c: a piece's files; _core for types
        "basin_run": "_core decompose elevation features io.mesh_index io.models pieces raster"
                     " refine_job run_record stats target_grid",  # 23c: the cut run
        # 23c-3: what one window's refine takes and gives, cut or not
        "refine_job": "_core budget elevation run_record target_grid",
        "land_cover_fractions": "_core crs io.models io.repository",  # 26b: class_areas
        # 105c: polygons cut per triangle; the table and its helpers from above
        "polygon_fractions": "land_cover_fractions",
        "_core": "",
    },
    {  # L4: pipelines
        "dem_input": "crs domain grow io.models io.repository mosaic resources target_grid",
        # 108: the rules compose their legend; 110: the land cover's lines take the
        # preferred classes (preferred imported inside the function, sideways in L4)
        "feature_input": "border_simplify crs domain features io.geojson io.geopackage io.gml"
                         " io.repository legends preferred",
        "catchment": "burn catchment_core crs gauge geodesy io.models io.repository mosaic"
                     " outline resources sensitivity",
        "catchment_batch": "catchment crs gauge geodesy hydrography io.repository reference"
                           " resources",
        "fetch": "",
        "gml_import": "crs io.geopackage_write io.gml",  # 73: GML 3.2 into a GeoPackage
        "fetch.plan": "crs domain fetch.http io.cog io.geotiff io.models",
        "fetch.run": "crs fetch.http fetch.plan io.models io.repository sources tin_engine",
        "fetch.nve": "fetch.http io.geojson sources",
        "fetch.features": "crs fetch.http fetch.plan io.geopackage_write io.repository sources"
                          " tin_engine",  # 63: an ArcGIS layer's polygons into the cache
        "tolerance_field": "crs feature_input",
        "sub_catchments": "border_simplify crs domain features",  # 106c, 106e: units; the domain
        "cells": "feature_input io.models",  # 105a: the cell table; Legend is feature_input's
        "preferred": "cells polygon_fractions",  # 107: classes from a preferred source
        "legends": "toml_data",  # 108: a legend file is the colour scheme only
        # 108: rules files; 111: renamed land-cover class files
        "land_cover_classes": "feature_input features legends toml_data",
        # 111: the one land-cover map, the --prefer cut on LandCoverClasses
        "land_cover_map": "feature_input preferred",
        # 55: the listing; 57: vector sources; 114: the box with its CRS
        "source_report": "crs dem_input domain feature_input geodesy io.geojson io.geopackage"
                         " land_cover_classes io.models io.repository sources",
        # 85: a network file's links and units; reads through feature_input
        "io.network_file": "crs feature_input io.geojson io.geopackage io.repository network",
        # 116: a feature layer with its columns, for plugins; FeatureError is feature_input's
        "io.feature_table": "feature_input io.geopackage io.repository",
    },
    {  # L5: flags in, files and stderr out
        "cli": "_core basin_run budget catchment catchment_batch cells chains legends cli_kit crs"
               " dem_input domain decompose edge_strip elevation feature_input features"
               " fetch.features fetch.http fetch.nve fetch.plan fetch.run final_check gauge"
               " geodesy gml_import grid_domain hydrography io.cog io.domain_file io.geojson"
               " io.geopackage_write io.mesh_index io.models io.network_file io.ply io.repository"
               " io.rivers io.station_set io.vtk_legacy land_cover_fractions landcover legends"
               " land_cover_classes land_cover_map"
               " mcp.http mcp.server mosaic network palettes plugins polygon_fractions preferred"
               " raster"
               " refine_job resources"
               " run_record source_report sources stats sub_catchments target_grid tin_engine"
               " tolerance_field",
        "__main__": "plugins",  # 37: `python -m tin_engine`, the server's child; 51: with plugins
        "cli_kit": "resources",  # 51: the CLI helpers a plugin's commands share with the core's
        "plugins": "cli",  # 51: the core's commands, then each installed plugin's
        "plugin_api": "catchment_core cli_kit crs dem_input domain feature_input io.feature_table"
                      " io.geopackage_write io.geotiff_write io.models io.vtk_read raster"
                      " resources run_record sources",  # 51: what a plugin may import; 116
        "mcp": "", "mcp.schema": "",  # 37: the MCP server over the commands
        "mcp.server": "mcp.schema mcp.workspace sources tin_engine",  # 57: the output root's
                                                                      # variable; 69: the scope
        # 69: MCP over HTTP; only `mcp.http` imports Starlette or uvicorn
        "mcp.http": "mcp.links mcp.schema mcp.server mcp.store mcp.workspace sources",
        "mcp.links": "", "mcp.store": "", "mcp.workspace": "",
    },
)
# fmt: on

TABLE: dict[str, tuple[int, frozenset[str]]] = {
    module: (layer, frozenset(imports.split()))
    for layer, rows in enumerate(LAYERS)
    for module, imports in rows.items()
}

# F12's edges against the rule, importer -> imported, each with the PR
# (section 6) that removes it.
UPWARD: dict[tuple[str, str], str] = {
    ("cli", "_core"): "H, audit-mesh-run",
}

PACKAGE = Path(tin_engine.__file__).parent
SCANNED = sorted(set(TABLE) - {"_core"})


def dotted(short: str) -> str:
    return short if short == "tin_engine" else f"tin_engine.{short}"


def shortened(name: str) -> str:
    return name.removeprefix("tin_engine.")


def breaks_the_rule(importer: str, imported: str) -> bool:
    """Up a layer, or into `_core` from outside L3."""
    importer_layer, imported_layer = TABLE[importer][0], TABLE[imported][0]
    return imported_layer > importer_layer or (imported == "_core" and importer_layer != 3)


def edges_against_the_rule() -> set[tuple[str, str]]:
    return {
        (importer, imported)
        for importer, (_, row) in TABLE.items()
        for imported in row
        if breaks_the_rule(importer, imported)
    }


def test_the_table_names_every_module() -> None:
    on_disk = {"_core"}
    for path in PACKAGE.rglob("*.py"):
        parts = path.relative_to(PACKAGE).with_suffix("").parts
        on_disk.add(".".join(parts[:-1] if parts[-1] == "__init__" else parts) or "tin_engine")
    assert on_disk - set(TABLE) == set(), "modules with no row"
    assert set(TABLE) - on_disk == set(), "rows for modules that are gone"


@pytest.mark.parametrize("module", SCANNED)
def test_each_module_imports_exactly_its_row(module: str) -> None:
    found = {shortened(n) for n in first_party_imports(import_module(dotted(module)))}
    row = TABLE[module][1]
    assert (sorted(found - row), sorted(row - found)) == ([], []), (
        f"{module}: (imports missing from its row, row edges it no longer imports)"
    )


def test_every_edge_goes_down_or_sideways() -> None:
    assert edges_against_the_rule() - set(UPWARD) == set()


def test_no_upward_exception_is_stale() -> None:
    # An exception for an edge the table no longer has, or for one the rule
    # allows, is deleted by the PR that made it so.
    assert set(UPWARD) - edges_against_the_rule() == set()


def test_no_module_imports_by_name() -> None:
    """No `importlib.import_module` or `__import__`: the imports `ast` cannot see."""
    found = []
    for module in SCANNED:
        tree = ast.parse(inspect.getsource(import_module(dotted(module))))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                func = node.func
                name = getattr(func, "attr", None) or getattr(func, "id", None)
                if name in {"import_module", "__import__"}:
                    found.append(f"{module}:{node.lineno}")
    assert found == []
