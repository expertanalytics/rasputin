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
        "io.models": "", "features": "", "sources": "", "run_record": "", "stats": "",
        "palettes": "", "hydrography": "", "tin_engine": "",
    },
    {  # L1: pure algorithms
        "crs": "run_record",
        "mosaic": "io.models",
        "target_grid": "crs domain io.models",
        "grid_domain": "io.models",
        "domain": "crs io.models",
        "chains": "domain features",
        "elevation": "", "outline": "", "sensitivity": "", "landcover": "", "decompose": "",
        "burn": "gauge io.models",
        "gauge": "hydrography",
        "reference": "io.models",
        "viz": "viz.scene viz.style viz.svg",
        "viz.fixtures": "", "viz.protocols": "", "viz.style": "",
        "viz.scene": "viz.protocols",
        "viz.svg": "viz.scene viz.style",
    },
    {  # L2: io/ codecs, and the one module that opens a connection
        "io": "io.ply io.vtk_legacy",
        "io.cog": "io.geotiff io.models",
        "io.domain_file": "crs domain io.geojson io.repository",
        "io.geojson": "crs",
        "io.geopackage": "", "io.gml": "", "io.mesh_checks": "",
        "io.geotiff": "crs io.models",
        "io.mesh_index": "io.models",
        "io.ply": "features io.mesh_checks",
        "io.repository": "io.cog io.geotiff io.models mosaic",
        "io.rivers": "hydrography io.station_set",
        "io.station_set": "hydrography io.geojson io.repository",
        "io.vtk_legacy": "features io.mesh_checks",
        "fetch.http": "tin_engine",
    },
    {  # L3: the only importers of _core
        "raster": "_core io.models",
        "edge_strip": "_core stats",
        "final_check": "_core stats target_grid",
        "catchment_core": "_core io.models raster",
        "_core": "",
    },
    {  # L4: pipelines
        "dem_input": "crs domain io.models io.repository mosaic target_grid",
        "feature_input": "crs domain features io.geojson io.geopackage io.gml io.repository",
        "catchment": "burn catchment_core crs gauge io.models io.repository mosaic outline"
                     " sensitivity",
        "catchment_batch": "catchment crs gauge hydrography io.repository reference",
        "fetch": "",
        "fetch.plan": "crs domain fetch.http io.cog io.geotiff io.models",
        "fetch.run": "crs fetch.http fetch.plan io.models io.repository sources tin_engine",
        "fetch.nve": "fetch.http io.geojson sources",
    },
    {  # L5: flags in, files and stderr out
        "cli": "_core catchment catchment_batch chains crs dem_input domain edge_strip"
               " elevation feature_input features fetch.http fetch.nve fetch.plan fetch.run"
               " final_check gauge grid_domain hydrography io.cog io.domain_file io.geojson"
               " io.mesh_checks io.models io.ply io.repository io.rivers io.station_set"
               " io.vtk_legacy landcover mosaic palettes raster run_record sources stats"
               " target_grid tin_engine viz.fixtures viz.protocols viz.scene viz.style viz.svg",
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
