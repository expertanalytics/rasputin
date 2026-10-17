"""Vector sources over the bowl, registered in a user's sources file. Test support only.

`docs/increments/57-vector-sources.md`, *7. Tests*, fixture V, beside 55's
fixture L (`dataset_fixtures.write_local`):

- the lake file, `dataset_fixtures.lake_file`: one Polygon, the bowl's lake
  box, with a `crs` member naming EPSG:25833 and one attribute, `name`;
- the land-cover GeoPackage: two features layers in EPSG:25833, `landcover`
  (row 1 the box x 500000 to 520000, y 6560000 to 6600000 less the lake,
  `Code_18` "311"; row 2 the lake box, "512") and `tiles` (one row), so an
  entry for it needs `layer = "landcover"`. `gpkg_fixtures.write_gpkg`
  declares `GEOMETRY` and leaves `gpkg_contents`' bounds unset, so the extent
  comes from the R-tree.

The figures below were measured by `@architect` with the main checkout at
`a0a1ad45` (the same tree as the branch base) by the scratch probe the design
names, and are rounded by the listing (coverage to 1e-6 degrees, a share to
1e-4).
"""

from __future__ import annotations

import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any

from shapely.geometry import box

from catchment_fixtures import lake_box
from dataset_fixtures import CREDIT, LICENCE, lake_file
from gpkg_fixtures import Layer, Row, write_gpkg

#: The land-cover layer, the second features layer, and the coded attribute.
LANDCOVER = "landcover"
OTHER = "tiles"
CODE = "Code_18"
#: The land-cover box, x0 y0 x1 y1 in EPSG:25833, and its R-tree extent.
LAND_BOX = (500_000.0, 6_560_000.0, 520_000.0, 6_600_000.0)
LAND_CODES = ("311", "512")

#: Ids the tests register; none names anything real.
LAKES_ID = "test-lakes"
COVER_ID = "test-cover"

LAKE_DESCRIPTION = "Lake polygons of the bowl, a test fixture."
COVER_DESCRIPTION = "Land-cover polygons of the bowl, CORINE codes, a test fixture."

#: Coverage, lon/lat W S E N, measured.
LAKE_COVERAGE = (15.169718, 59.355479, 15.182063, 59.361781)
COVER_COVERAGE = (15.0, 59.178672, 15.353647, 59.538349)
#: Shares of 55's boxes (`dataset_fixtures.INSIDE`, `PARTLY`, `OUTSIDE`), measured.
LAKE_SHARE_INSIDE = 0.031
LAKE_SHARE_PARTLY = 0.0
COVER_SHARES = (1.0, 0.3792, 0.0)

#: Measured end to end: the catchment of the bowl's seed with the lake, and
#: the mesh of it at tolerance 5 with the land cover under the map corine.
CATCHMENT_NODES = 766
MESH_TRIANGLES = 1083
TRIANGLES_BY_CODE = {311: 1004, 512: 79}


def landcover_gpkg(path: Path, *, rtree: bool = True) -> Path:
    """The land-cover GeoPackage of fixture V."""
    lake = lake_box()
    rows = [Row(1, box(*LAND_BOX).difference(lake), {CODE: "311"}), Row(2, lake, {CODE: "512"})]
    other = [Row(1, box(500_000.0, 6_560_000.0, 501_000.0, 6_561_000.0), {CODE: "999"})]
    layers = [Layer(LANDCOVER, 25833, rows, rtree=rtree), Layer(OTHER, 25833, other, rtree=rtree)]
    return write_gpkg(path, layers)


def set_contents_bounds(path: Path, table: str, bounds: tuple[float, ...] | None) -> None:
    """`gpkg_contents`' four bounds of `table` set (None: unset)."""
    values = (None,) * 4 if bounds is None else bounds
    with closing(sqlite3.connect(path)) as con:
        con.execute(
            "UPDATE gpkg_contents SET min_x = ?, min_y = ?, max_x = ?, max_y = ? "
            "WHERE table_name = ?",
            (*values, table),
        )
        con.commit()


def vector_entry(path: Path, description: str, **changes: Any) -> dict[str, Any]:
    """A `kind = "local"` table for a vector file, stating no `values`; a
    change to None drops that key."""
    entry: dict[str, Any] = {
        "kind": "local",
        "path": str(path),
        "credit": CREDIT,
        "licence_note": LICENCE,
        "description": description,
    }
    entry.update(changes)
    return {k: v for k, v in entry.items() if v is not None}


def fixture_v(directory: Path) -> tuple[Path, Path]:
    """Fixture V's two files in `directory`: the lake file and the land-cover GeoPackage."""
    directory.mkdir(parents=True, exist_ok=True)
    return lake_file(directory / "lakes.geojson"), landcover_gpkg(directory / "landcover.gpkg")


def fits_of(entry: dict[str, Any]) -> set[tuple[str, str, str | None]]:
    """A listed source's fits as (command, option, legend)."""
    return {(f["command"], f["option"], f.get("legend")) for f in entry["fits"]}
