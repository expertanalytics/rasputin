"""Stations, rivers, references and two tile sets for increment 29 PR 4's suites.

`docs/increments/29-nve-reference-catchments.md`, "The batch" and "The red
suites", PR 4. Imports only what PRs 1 to 3 and 2 shipped (`two_basins`,
`accumulate`, `upstream`, the river and station readers), so a missing PR 4
module fails the tests that use it, not the collection of this helper.

THE TERRAIN is `gauge_fixtures.two_basins()` (300 x 700 nodes, 10 m): a
valley down column `BIG_CC` = 100 and a large tributary basin to the east
joining it along row `JOIN_ROW` = 152. The design names the valley DEM of
`test_cli_catchment.py` for the batch; it has no confluence, so the station
that must be uncertain by a confluence below it could not be built on it, and
this terrain (already the stage B terrain of `test_catchment.py`) is used.
Along column 100 from row 60 down, the floor falls strictly (0.5 m a row,
then 2 m a row in the exit slot from row 180), so a burn along it lowers
nothing: `burnt == raw` is asserted as a premise, and the oracle for each
placed catchment is one flood of the raw terrain.

THE RIVERS (one GeoJSON in EPSG:25833, digitised downstream):

- `A` (objectid 101, `elvid` 2-1-1, watercourse 002.A, "Hovedelva"): down
  column 100 from row 60 to the south edge; and an exact copy, 102, which
  the reader drops (so "4 segments read, 1 exact copy dropped" without D).
- `B` (201, 2-2-1, 002.B, "Kortelva"): down column 100 from row 60 to row
  131, so it ends 10 m below a station at row 130, inside `U` = 30 m.
- `C` (301, 2-3-1, 002.C, "Sideelva"): the tributary, along row 152 from
  column 230 west to column 101.
- `D` (401, 2-4-1, 002.D, "Grenseelva"): along row 152 from column 690 west to
  column 500, for the mixed-grid DEM only.

THE STATIONS, each 16 m east of its line (so `U` = 30 m), in file order:

1. `1.140.0` "Treff" on A at row 140: well posed (swing 0.024 on the raw
   counts); its reference is its own flood: `match` by overlap.
2. `1.135.0` "Bom" on A at row 135: well posed; its reference is its own
   flood moved 3 km east: `miss`.
3. `1.150.0` "Samløp" on A at row 150: the tributary joins 20 m below, so
   the area 30 m down is 13 times the placed one: `uncertain` (swing),
   although its reference is its own flood.
4. `1.130.0` "Slutt" on B at row 130: B ends 10 m below it, the downstream
   side is not read to `U`: `uncertain` (downstream_unread), although its
   reference is its own flood (overlaps 100 %).
5. `1.40.0` "Langt" at column 600, row 40, over 3 km from any line:
   `refused`, `refusal_cause` "no_river".

THE MIXED-GRID DEM: the same terrain in two tiles, `a.tif` (columns 0 to
649) on the main lattice and `b.tif` (columns 640 to 699) moved 5 m east,
half a cell. A station on `D` (`2.600.0` "Grense", column 600) has a first
window (columns 384 to the east edge) that selects both and that neither
covers; the station on A at row 140 grows its window twice, to column 571
(measured on PR 2's code), inside `a.tif` alone.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import shapely
from shapely import affinity
from shapely.geometry import MultiPolygon, Polygon, box

import gauge_fixtures as gf
from catchment_fixtures import filled
from mosaic_fixtures import quadrants, whole
from nve_fixtures import collection, feature, point, river, write
from tin_engine._core import upstream
from tin_engine.io.models import DemTile
from tin_engine.raster import to_core

COL = gf.BIG_CC
#: Metres east of its line each station stands.
EAST_M = 16.0
LINE_TOP = 60


@dataclass(frozen=True)
class Spec:
    station: str
    name: str
    col: float
    row: int
    watercourse: str | None
    river: str | None


TREFF = Spec("1.140.0", "Treff", COL, 140, "002.A", "Hovedelva")
BOM = Spec("1.135.0", "Bom", COL, 135, "002.A", "Hovedelva")
SAMLOP = Spec("1.150.0", "Samløp", COL, 150, "002.A", "Hovedelva")
SLUTT = Spec("1.130.0", "Slutt", COL, 130, "002.B", "Kortelva")
LANGT = Spec("1.40.0", "Langt", 600, 40, None, None)
FIVE = (TREFF, BOM, SAMLOP, SLUTT, LANGT)
GRENSE = Spec("2.600.0", "Grense", 600, 152, "002.D", "Grenseelva")


def xy(spec: Spec) -> tuple[float, float]:
    x, y = gf.lat(spec.col, spec.row)
    if spec is GRENSE:  # on a line along a row: 16 m north of it
        return x, y + EAST_M
    return x + EAST_M, y


def placed_node(spec: Spec) -> tuple[int, int]:
    return (spec.row, round(spec.col))


def basins() -> np.ndarray:
    return gf.two_basins()


def basin_tiles() -> dict[str, DemTile]:
    tile = gf.tile_of(basins())
    return quadrants(tile, row_cut=gf.BIG_ROWS // 2, col_cut=gf.BIG_COLS // 2, overlap=1)


def mixed_tiles() -> dict[str, DemTile]:
    z = basins()
    a = whole(
        gf.BIG_ROWS,
        650,
        array=np.ascontiguousarray(z[:, :650]),
        x_min=gf.X0,
        y_max=gf.Y0,
        dx=gf.D,
        dy=gf.D,
    )
    b = whole(
        gf.BIG_ROWS,
        60,
        array=np.ascontiguousarray(z[:, 640:]),
        x_min=gf.X0 + 640 * gf.D + gf.D / 2,
        y_max=gf.Y0,
        dx=gf.D,
        dy=gf.D,
    )
    return {"a.tif": a, "b.tif": b}


def river_features(mixed: bool = False) -> list[dict[str, Any]]:
    a = list(gf.column_line(COL, LINE_TOP, gf.BIG_ROWS - 1))
    b = list(gf.column_line(COL, LINE_TOP, 131))
    c = [gf.lat(230, gf.JOIN_ROW), gf.lat(101, gf.JOIN_ROW)]
    d = [gf.lat(690, gf.JOIN_ROW), gf.lat(500, gf.JOIN_ROW)]
    features = [
        river(101, a, elvid="2-1-1", vassdragsnr="002.A", elvenavn="Hovedelva"),
        river(102, a, elvid="2-1-1", vassdragsnr="002.A", elvenavn="Hovedelva"),
        river(201, b, elvid="2-2-1", vassdragsnr="002.B", elvenavn="Kortelva"),
        river(301, c, elvid="2-3-1", vassdragsnr="002.C", elvenavn="Sideelva"),
    ]
    if mixed:
        features.append(river(401, d, elvid="2-4-1", vassdragsnr="002.D", elvenavn="Grenseelva"))
    return features


def write_rivers(path: Path, mixed: bool = False, crs: str | None = gf.EPSG) -> Path:
    return write(path, collection(river_features(mixed), crs=crs))


def station_feature(spec: Spec) -> dict[str, Any]:
    props: dict[str, Any] = {
        "station": spec.station,
        "name": spec.name,
        "series": ["1001.0"],
        "nve_area_km2": 12.5,
    }
    if spec.watercourse is not None:
        props["watercourse"] = spec.watercourse
    if spec.river is not None:
        props["river"] = spec.river
    return point(*xy(spec), **props)


def write_stations(path: Path, specs: Sequence[Spec] = FIVE, crs: str | None = gf.EPSG) -> Path:
    return write(path, collection([station_feature(s) for s in specs], crs=crs))


def flood_mask(spec: Spec) -> np.ndarray:
    """The placed node's catchment on the raw terrain (the burn lowers
    nothing on column 100 from row 60 down)."""
    z = basins()
    seed = np.zeros(z.shape, dtype=np.uint8)
    seed[placed_node(spec)] = 1
    return np.asarray(upstream(to_core(gf.tile_of(z)), seed).mask, dtype=np.uint8)


def cells_polygon(mask: np.ndarray) -> Polygon | MultiPolygon:
    """The union of the in-nodes' cells (edges halfway between nodes), built
    from row runs: every in-node strictly inside, every out-node outside."""
    boxes = []
    for r in range(mask.shape[0]):
        cols = np.flatnonzero(mask[r])
        if not cols.size:
            continue
        breaks = np.flatnonzero(np.diff(cols) > 1)
        starts = np.r_[cols[0], cols[breaks + 1]]
        ends = np.r_[cols[breaks], cols[-1]]
        for c0, c1 in zip(starts, ends, strict=True):
            x0, y0 = gf.lat(c0 - 0.5, r + 0.5)
            x1, y1 = gf.lat(c1 + 0.5, r - 0.5)
            boxes.append(box(x0, y0, x1, y1))
    out = shapely.union_all(boxes)
    assert isinstance(out, Polygon | MultiPolygon)
    return out


def reference_of(spec: Spec, shift_m: float = 0.0) -> Polygon | MultiPolygon:
    """The station's own flood, holes filled as our fine outline fills them,
    as cells; moved `shift_m` east."""
    poly = cells_polygon(filled(flood_mask(spec)))
    return affinity.translate(poly, xoff=shift_m) if shift_m else poly


def references() -> dict[str, Polygon | MultiPolygon]:
    """Each of the five stations' reference, as THE STATIONS describes."""
    return {
        TREFF.station: reference_of(TREFF),
        BOM.station: reference_of(BOM, shift_m=3000.0),
        SAMLOP.station: reference_of(SAMLOP),
        SLUTT.station: reference_of(SLUTT),
        LANGT.station: box(*gf.lat(590, 50), *gf.lat(610, 30)),
    }


def write_references(path: Path, crs: str | None = gf.EPSG) -> Path:
    features = [feature(p, station=n, reference_area_km2=p.area / 1e6)
                for n, p in references().items()]  # fmt: skip
    return write(path, collection(features, crs=crs))


def swing_on_raw(spec: Spec, u_rows: int = 3) -> float:
    """The one-sided swing within `u_rows` along column 100, from the raw
    terrain's accumulation (what the sensitivity reads where the burn lowers
    nothing)."""
    from tin_engine._core import accumulate

    count = np.asarray(accumulate(to_core(gf.tile_of(basins()))).count, dtype=np.float64)
    r, c = placed_node(spec)
    a0 = count[r, c]
    return float(max(a0 - count[r - u_rows, c], count[r + u_rows, c] - a0) / a0)


# ---------------------------------------------------------------------------
# PR 4, lake gauges ("Lake gauges: the lake is the seed")
# ---------------------------------------------------------------------------
#
# THE LAKE is a box over the main valley's floor round `SAMLOP` (row 150,
# 16 m east of column 100): nodes of rows 144 to 151 and columns 97 to 103,
# its edges halfway between nodes, so `SAMLOP` lies inside it and the
# tributary's slot (row 152) does not. Measured on PR 2's code by
# `delineate` with this lake (22's path): 5715 nodes, against 5540 for
# `SAMLOP`'s placed node on the river path, so a lake row and a river row
# are told apart by their node counts. `OVERLAPPING` is a second box inside
# it, also round `SAMLOP`. `LAKE_LINE_LAKE` is a box over column 100 at
# `TREFF` (row 140) whose east edge is column 101, so `TREFF`'s `P` (column
# 100) is inside and the station (column 101.6) is 6 m outside.

LAKE_NUMBER, LAKE_NAME = 4110, "Samløpvatnet"


def lake_polygon() -> Polygon:
    return box(*gf.lat(96.5, 151.5), *gf.lat(103.5, 143.5))


def overlapping_polygon() -> Polygon:
    return box(*gf.lat(99.5, 150.5), *gf.lat(102.5, 148.5))


def lake_line_polygon() -> Polygon:
    return box(*gf.lat(97.5, 143.5), *gf.lat(101.0, 136.5))


def lake_feature(
    polygon: Polygon, objectid: int, vatnlnr: int | None = LAKE_NUMBER,
    navn: str | None = LAKE_NAME,
) -> dict[str, Any]:  # fmt: skip
    """One `lakes.geojson` feature with the four fields `fetch-stations` writes."""
    return feature(polygon, objectid=objectid, vatnlnr=vatnlnr, navn=navn,
                   areal_km2=polygon.area / 1e6)  # fmt: skip


def write_lakes(
    path: Path, features: Sequence[dict[str, Any]] | None = None, crs: str | None = gf.EPSG
) -> Path:
    """The lake file: by default THE LAKE alone, objectid 1."""
    chosen = [lake_feature(lake_polygon(), 1)] if features is None else list(features)
    return write(path, collection(chosen, crs=crs))


def lake_line_rivers() -> list[dict[str, Any]]:
    """THE RIVERS with `A` (and its copy) mapped as a lake centreline."""
    out = river_features()
    for f in out[:2]:
        f["properties"] |= {"objekttype": "InnsjøMidtlinje", "vatnlnr": LAKE_NUMBER}
    return out
