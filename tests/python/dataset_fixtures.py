"""Local and cached elevation sources, registered in a user's sources file.

`docs/increments/55-dataset-discovery.md`, *7. Tests*. Two datasets of the
bowl terrain (`catchment_fixtures.bowl`, 400 x 200 nodes at 100 m in
EPSG:25833, corner node (500000, 6600000)):

- Fixture L, a local dataset: the bowl cut into four overlapping GeoTIFFs in
  one directory, as `test_cli_catchment` cuts it. `COARSE` is the same bowl
  as one file at 200 m (every second node), so two local sources differ in
  cell size alone.
- Fixture C, a cached remote-style source: the bowl as one tiled GeoTIFF
  (16 x 16 blocks, Deflate, 25 x 13 = 325 blocks) in the tile-cache layout
  `cog_fixtures.write_cache` writes.

The figures below were measured by `@architect` with the main checkout at
`9a7c312d` (same tree as the branch base), by scratch probes run through the
CLI and the library; the design file lists them.

Test support only: nothing in `src_python` uses it.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from shapely.geometry import mapping

from catchment_fixtures import EPSG, X0, Y0, bowl, lake_box, lat, tile_of
from geotiff_fixtures import micro_tiff
from mosaic_fixtures import quadrants, whole
from test_cli_mesh_mosaic import write_tiles

#: Ids the tests register; none is a shipped id.
LOCAL = "bowl-100m"
COARSE = "bowl-200m"
CACHED = "bowl-cached"
#: The package's entries: three rasters and two vector layers fetched from a map service.
SHIPPED = ("anadem-v1", "glo30", "mapbiomas-c11", "nve-lakes", "nve-glaciers", "clc2018")

CREDIT = "Bowl terrain, a test fixture of rasputin"
LICENCE = "Test fixture only"
DESCRIPTION = "Bare-ground test terrain, 100 m, a bowl with a lake and an outlet channel."

#: Fixture L's coverage, lon/lat W S E N, as measured; compared to 1e-6 degrees.
COVERAGE = (15.0, 59.179575, 15.351879, 59.538349)
#: Boxes (W S E N, degrees) and the share of each that fixture L covers, measured.
INSIDE = (15.15, 59.33, 15.20, 59.38)
PARTLY = (15.2, 59.3, 15.6, 59.4)
OUTSIDE = (16.0, 60.0, 16.1, 60.1)
SHARE_PARTLY = 0.3748

#: The bowl's centre node, inside the lake: lattice (col 100, row 200).
SEED = lat(100, 200)
#: Catchment node counts, measured: from the seed node alone, and with the lake.
NODES_NO_LAKE = 623
NODES_WITH_LAKE = 766

#: Fixture C: blocks per row, block rows, and the block row holding the seed (row 200 // 16).
BLOCKS_ACROSS = 13
BLOCK_ROWS = 25
SEED_BLOCK_ROW = 12


def write_local(directory: Path) -> Path:
    """Fixture L: four tiles in `directory`, which is returned."""
    write_tiles(directory, quadrants(tile_of(bowl()), row_cut=200, col_cut=100, overlap=1))
    return directory


def write_coarse(path: Path) -> Path:
    """The bowl at 200 m: every second node of fixture L, one file."""
    z = bowl()[::2, ::2]
    rows, cols = z.shape
    tile = whole(rows, cols, array=z, x_min=X0, y_max=Y0, dx=200.0, dy=200.0)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tiles(path.parent, {path.name: tile})
    return path


def cached_bytes() -> bytes:
    """Fixture C's one object: the bowl as a 16 x 16 tiled Deflate GeoTIFF."""
    return micro_tiff(
        bowl(),
        tiepoint=(0.0, 0.0, 0.0, X0, Y0, 0.0),
        scale=(100.0, 100.0, 0.0),
        nodata="-9999",
        tile=(16, 16),
        compression="deflate",
    ).getvalue()


def block_row(row: int) -> list[int]:
    """The block indices of fixture C's block row `row`."""
    return list(range(row * BLOCKS_ACROSS, (row + 1) * BLOCKS_ACROSS))


def local_entry(path: Path | None, **changes: Any) -> dict[str, Any]:
    """A `kind = "local"` table for `path` (None: no `path` key), heights in metres.
    A change to None drops that key."""
    entry: dict[str, Any] = {
        "kind": "local",
        "path": None if path is None else str(path),
        "z_unit": "metre",
        "credit": CREDIT,
        "licence_note": LICENCE,
        "description": DESCRIPTION,
    }
    entry.update(changes)
    return {k: v for k, v in entry.items() if v is not None}


def cached_entry(**changes: Any) -> dict[str, Any]:
    """A one-COG remote table in EPSG:25833, as fixture C was fetched from."""
    entry: dict[str, Any] = {
        "kind": "one-cog",
        "url": "https://example.invalid/bowl.tif",
        "crs": EPSG,
        "nodata": -9999.0,
        "z_unit": "metre",
        "credit": CREDIT,
        "licence_note": LICENCE,
        "description": DESCRIPTION,
    }
    entry.update(changes)
    return {k: v for k, v in entry.items() if v is not None}


def _toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int | float):
        return repr(value)
    if isinstance(value, Mapping):  # an inline table
        return "{ " + ", ".join(f"{k} = {_toml_value(v)}" for k, v in value.items()) + " }"
    if isinstance(value, Sequence) and not isinstance(value, str):
        return "[" + ", ".join(_toml_value(v) for v in value) + "]"
    return json.dumps(str(value))  # a TOML basic string for ASCII text


def write_sources(path: Path, entries: Mapping[str, Mapping[str, Any]]) -> Path:
    """A user sources file: one `[source.<id>]` table per entry."""
    lines: list[str] = []
    for key, table in entries.items():
        lines.append(f"[source.{key}]")
        lines += [f"{name} = {_toml_value(value)}" for name, value in table.items()]
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def lake_file(path: Path) -> Path:
    """The bowl's lake (`lake_box`) as GeoJSON with a `crs` member naming EPSG:25833."""
    doc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": EPSG}},
        "features": [
            {"type": "Feature", "properties": {"name": "the lake"},
             "geometry": mapping(lake_box())},
        ],
    }  # fmt: skip
    path.write_text(json.dumps(doc))
    return path


def by_id(document: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """A listing's sources by id; each id once."""
    sources = document["sources"]
    found = {s["id"]: s for s in sources}
    assert len(found) == len(sources), [s["id"] for s in sources]
    return found
