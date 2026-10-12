# Increment 40: `catchment` on a DEM in degrees, flooding its own grid

**Status: designed** (`@architect`, 2026-10-12, branch `worktree-geo-catchment`,
from master `f09b9863`). Ola was away when it was written; every open point
is decided on a stated default, listed in section 10 for the morning. Every
line citation below is pinned to `f09b9863`. Increment 38 (threads and memory
as settings, branch `worktree-threads`, approved, not yet merged) changes
`delineate` and `_grow` in `src_python/tin_engine/catchment.py` (a `budget`
argument); the red and green steps start from master after 38 merges, and the
shapes named here survive it.

## 1. What Ola asked

The generality audit's row A1 (`docs/increments/generality-audit.md`, section
2, merged as #237): `rasputin catchment` refuses every DEM in degrees, so
GLO-30 and ANADEM, so the São Francisco basin. Its question Q5 offered
"resample onto the suggested projected grid as `mesh --out-crs` does, or
flood in degrees". The main session put it to Ola as "Shall the design for
geographic catchments flood the DEM's own grid, with the hydrology constants
expressed in cells and geodesic areas, instead of resampling?"

Ola, verbatim, 2026-10-12: "on 5: yes".

So: no resampling. The flood runs on the DEM's own longitude-latitude grid;
distances that mean "a few DEM cells" become cell counts; distances that mean
metres stay metres and are converted where they are used; areas are geodesic.
The audit's row A9 (`--rivers` only in the DEM's CRS) comes along, because a
river file in metres over a DEM in degrees is the common case. Row E4 (the
gauge and burn constants in metres) is settled here for the catchment path.

## 2. Prior art: legacy and literature

*Literature.* The flood is Priority-Flood (Barnes, Lehman and Mulla 2014,
"Priority-Flood: an optimal depression-filling and watershed-labeling
algorithm for digital elevation models", Computers & Geosciences 62, 117-127,
doi:10.1016/j.cageo.2013.04.024), already shipped in increments 22 and 29. It
is a graph algorithm: a priority queue of nodes keyed by filled elevation,
eight neighbours per node, and no distance anywhere. Ours is the same: the
queue key is (level, push counter)
(`include/terrain/hydrology/flood.hpp@f09b9863:81-82`), the neighbours are
the eight index neighbours (`include/terrain/hydrology/flood.hpp@f09b9863:23`),
and the only things read from the raster are its shape, its NoData mask and
its values (`include/terrain/hydrology/flood.hpp@f09b9863:42`,
`include/terrain/hydrology/flood.hpp@f09b9863:53`). A node drains to the
neighbour that reached it first, which on sloping ground is its lowest
already-reached neighbour, not the steepest by slope over distance. That holds
today on square cells too (a diagonal is never weighted by its 1.41 times the
length), so a cell that is 28 m by 30 m changes nothing the method assumes.
This increment departs from nothing in the method.

How other tools treat a grid in degrees, read 2026-10-12:

- **TauDEM** (Tarboton; version 5.3 added geographic coordinates) computes a
  metre width and height for every row: `src/tiffIO.cpp` on the master branch
  of `github.com/dtarb/TauDEM` fills `dxc[j]`, `dyc[j]` per row `j` from the
  row's latitude and the ellipsoid (`geotoLength`), and uses them for slopes
  and contributing areas. This is the per-row approach this design takes for
  areas.
- **WhiteboxTools**: the fetched source of `d8_pointer.rs` and
  `d8_flow_accum.rs` (master of `github.com/jblindsay/whitebox-tools`) shows
  no latitude handling; the diagonal is `sqrt(cell_size_x^2 + cell_size_y^2)`
  in the raster's own units. Its release notes say other tools compute a
  z-factor per row on geographic rasters.
- **GRASS `r.watershed`**: `raster/r.watershed/ram/init_vars.c` (main of
  `github.com/OSGeo/grass`) forms the diagonal from the region's `ew_res` and
  `ns_res` as stored, which in a latitude-longitude location are degrees; the
  manual says nothing on latitude-longitude locations.

So per-row metres for areas is established practice (TauDEM); flooding the
native grid without distances is what Priority-Flood does anywhere. No novelty
is claimed, so no novelty search was needed.

*Legacy.* Nothing to carry. The legacy tree holds no hydrology, and its
geographic code is GeoTIFF key parsing, solar position and a Boost.Geometry
transform for the solar code:

```
$ git grep -l -i -E 'watershed|catchment|priority.?flood|flow.?accum' legacy-archive -- legacy
(no files)
$ git grep -l -i -E 'latitude|geographic|4326|Geod' legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/avalanche.py
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/reader.py
legacy-archive:legacy/rasputin/solar_position.h
legacy-archive:legacy/rasputin/triangulate_dem.h
legacy-archive:legacy/rasputin/web_visualize.py
legacy-archive:legacy/tests/test_gml_repository.py
legacy-archive:legacy/tests/test_land_cover_repository.py
```

## 3. Measured before designing

All on Ola's M1 Max, on AC power, 2026-10-12, with today's core (built from
the same `include/` and `bindings/` as `f09b9863`: the extension was taken
from the `se2-catchments` worktree, whose head is `f74df10c`, and
`git diff --stat f74df10c f09b9863 -- include bindings` is empty). This
worktree may not build C++ in the design step. The probes are in
`docs/increments/40-probes/`; each says how it is run.

**Today's refusal.** `rasputin catchment --dem tests/fixtures/velhas/anadem_velhas.tif
--seed 597180 7852380 --seed-crs EPSG:31983 --out c.geojson` exits 2 with
"Invalid value for --dem: a geographic DEM (EPSG:4674, degrees) cannot enter
the core", raised at `src_python/tin_engine/raster.py@f09b9863:26-27`.

**The flood on the native grid** (`40-probes/native_flood.py`). The fixture is
320 x 320 nodes of ANADEM in EPSG:4674 (SIRGAS 2000), spacing
0.000269494585° both ways; at its centre latitude, 19.4436° S, a cell is
28.2996 m east-west and 29.8324 m north-south. From the seed of the fixture's
reference catchment (`tests/fixtures/velhas/NOTICE`, section 2), node (74, 224):

| | value |
|---|---|
| `upstream` on a unit lattice | 12.0 ms; 12,725 nodes; not cut by the edge or NoData |
| `accumulate` | 12.6 ms |
| area, nodes times each row's geodesic cell area | 10.7445 km² |
| area, geodesic, of the traced outline | 10.7441 km² (one ring, 678 vertices) |
| reference: the resampled route at 30 m in EPSG:31983 (`tests/fixtures/velhas/catchment.geojson`, reduced outline) | 10.6952 km² |
| shared area / reference, shared area / ours | 0.9903, 0.9864 |

So flooding the native grid gives the same catchment as the resampled route,
to within the outline's own reduction and one resampling.

**The whole command, prototyped in a scratch copy** (not committed): with only
the hydrology handed a unit lattice and the window margin set to 0.018°, the
unchanged window loop, mosaic and trace ran on the fixture: three windows
(136 x 135, 263 x 296, 320 x 320 nodes), 12,725 nodes, contained. Its report
showed exactly what this design fixes: areas in square degrees ("0.000000
km2"), the window printed as "x -44--44", a tolerance of 0.000539 "m", and
the start point printed as `(np.float64(-44.07…), …)`, a defect that is
there on a projected DEM too.

**Scaling of the flood** (`40-probes/flood_scaling.py`, the fixture mirrored
and tiled): 0.10 M nodes 0.012 s; 6.55 M nodes 1.14 s; 40.96 M nodes 11.4 s
(`upstream`; `accumulate` within 10 % of it). Throughput falls from 8.6 to
3.6 M nodes/s as the queue grows.

**Metres in the catchment path** (`grep` over `catchment.py`, `burn.py`,
`gauge.py`, `sensitivity.py`, `catchment_batch.py` and the `catchment`
command): `WINDOW_MARGIN_M = 2000`
(`src_python/tin_engine/catchment.py@f09b9863:63`); `JOIN_M = 1.0`,
`U_FLOOR_M = 30.0`, `DOWN_EXTRA_M = 100.0`, `CONFLUENCE_M = 100.0`,
`LAKE_GAP_M = 30.0` (`src_python/tin_engine/gauge.py@f09b9863:29-38`); the
reach's `corridor = 30.0` (`src_python/tin_engine/gauge.py@f09b9863:64`);
`map_radius = 500.0`, `reach_up = 1000.0`
(`src_python/tin_engine/gauge.py@f09b9863:150-151`, and the CLI's defaults at
`src_python/tin_engine/cli.py@f09b9863:2449-2450`); `DROP_M = 0.001`,
`END_CAP_M = 500.0`, `SLACK_M = 1e-6`, `DIRECTION_M = 2.0`
(`src_python/tin_engine/burn.py@f09b9863:25-31`); the unnamed 100 m under
which the direction check is skipped
(`src_python/tin_engine/burn.py@f09b9863:178`); `_SLACK_M = 1e-6`
(`src_python/tin_engine/sensitivity.py@f09b9863:24`). `grow.py`, which the
task named, is the mesh's outline growth (imported only by `target_grid.py`
and `dem_input.py`), not the catchment's.

**Other places that assume metres:** the default outline tolerance, twice
the cell (`src_python/tin_engine/catchment.py@f09b9863:247`); the window margin
added to coordinates (`src_python/tin_engine/catchment.py@f09b9863:232`,
`src_python/tin_engine/catchment.py@f09b9863:321`); the reach's downstream
length (`src_python/tin_engine/catchment.py@f09b9863:345`); the cell area
handed to the sensitivity (`src_python/tin_engine/catchment.py@f09b9863:346`);
the lake, outline and hole areas (`src_python/tin_engine/catchment.py@f09b9863:491-498`);
the outline reduction "in metres from its lower-left corner"
(`src_python/tin_engine/catchment.py@f09b9863:512`); and the command's report
(`src_python/tin_engine/cli.py@f09b9863:2472`,
`src_python/tin_engine/cli.py@f09b9863:2477`,
`src_python/tin_engine/cli.py@f09b9863:2493-2508`).

## 4. The design in one picture

```
 seed, lakes, river file (any CRS)          DEM tiles (one CRS: metres or degrees)
            |                                          |
            v                                          v
   moved into the DEM's CRS (pyproj)        window planned and assembled
            |                                (mosaic; margins converted to
            |                                 the DEM's units by a LocalFrame)
            v                                          |
   gauge placed in metres  ----reach---->   burn, in metres (LocalFrame of the window)
   (LocalFrame of the river file                       |
    when that file is in degrees)                      v
                                            _core.upstream / accumulate on a
                                            UNIT LATTICE: shape, values, NoData.
                                            No degrees, no metres, no CRS.
                                                       |
                                                       v
                                            outline traced on the native grid,
                                            in the DEM's CRS
                                                       |
                         areas: geodesic (pyproj Geod) | reduction: in metres
                         by row for nodes, over the    | (LocalFrame at the
                         ring for outlines             |  outline's centre)
                                                       v
                                  GeoJSON in the DEM's CRS, which mesh --domain reads
```

The C++ core changes not at all. Everything new is Python, in one new module
and small edits to the callers.

## 5. Components and boundaries

### 5.1 The hydrology sees a unit lattice (`raster.py`, `catchment_core.py`)

`raster.py` stays the one module that builds a core raster. It gains:

```python
def to_core_lattice(tile: DemTile) -> RasterView:
    """The tile's array on a unit lattice (x = column, y = -row, both spacings 1),
    for the hydrology, which reads only the shape, the values and NoData. A tile
    in degrees is accepted: no coordinate crosses."""
```

Read-only array required, as `to_core`. `to_core` keeps its refusal of a
geographic tile (`src_python/tin_engine/raster.py@f09b9863:26-27`): the mesh
path still needs metres. `catchment_core.upstream` and
`catchment_core.accumulate` call `to_core_lattice` instead of `to_core`, for
every DEM. Why that is safe on a projected DEM: `upstream`, `accumulate` and
`flood` read `geometry()` only for `rows()`, `cols()` and `size()`
(`include/terrain/hydrology/flood.hpp@f09b9863:42`,
`include/terrain/hydrology/upstream.hpp@f09b9863:36`,
`include/terrain/hydrology/accumulate.hpp@f09b9863:35`), so their results do
not depend on the affine scalars. The core boundary rule holds: no CRS, no
unit and no path crosses.

**C++ needs nothing.** The D8 weight question does not arise: the flood never
weights a neighbour by distance (section 2), and the accumulation counts
nodes; areas are weighted in Python (5.4).

### 5.2 `LocalFrame`: metres near one point (`geodesy.py`, new, layer L1)

```python
@dataclass(frozen=True, slots=True)
class LocalFrame:
    """An affine map from a CRS's coordinates to metres near one point.
    Identity for a projected CRS. For a CRS in degrees:
    X = (lon - lon0) * kx,  Y = (lat - lat0) * ky, with kx = (pi/180) N cos(lat0)
    and ky = (pi/180) M, N and M the ellipsoid's prime-vertical and meridian
    radii at lat0: the exact metres per degree at lat0."""
    lon0: float
    lat0: float
    kx: float
    ky: float

    @classmethod
    def identity(cls) -> LocalFrame: ...                       # (0, 0, 1, 1)
    @classmethod
    def at(cls, crs: str, point: tuple[float, float]) -> LocalFrame: ...
    @property
    def is_identity(self) -> bool: ...
    def forward(self, xy: npt.ArrayLike) -> npt.NDArray[np.float64]: ...   # (N, 2)
    def inverse(self, xy: npt.ArrayLike) -> npt.NDArray[np.float64]: ...
    def units(self, metres: float) -> tuple[float, float]: ... # (metres/kx, metres/ky)
    def cell_m(self, meta: RasterMeta) -> float: ...           # max(dx*kx, dy*ky)
    def meta(self, meta: RasterMeta) -> RasterMeta: ...        # the grid, in the frame
```

Rules:

- `at` reads the ellipsoid with `parse_crs(crs).get_geod()` (its `a` and
  `es`), never a hard-coded WGS 84. A geographic CRS whose axis unit is not
  the degree is a `ValueError` in plain words (GeoTIFF DEMs in grads are not
  known to exist; the check costs a line).
- On the identity frame `forward`, `inverse` and `units` return their inputs'
  values unchanged, and `meta` returns the same `RasterMeta`, so every
  projected path computes bit for bit what it does today.
- `meta(m)` maps the node grid: `x_min' = (x_min - lon0) kx`,
  `delta_x' = delta_x kx`, likewise in y; `geographic=False`, and `crs` a
  label naming the frame, never parsed. It exists for the burn only (5.5).
- Why affine and not per row: an affine map keeps straight lines straight, so
  the outline reduction's guarantees (a simple polygon, the seed inside)
  survive the way back, exactly. Per-row scaling would not. The cost is a
  scale error east-west of `|cos(lat)/cos(lat0) - 1|`, about
  `tan(lat0) * (lat - lat0)` in radians: over a 1.5 km reach at 20° S, 9e-5
  (0.13 m); at 60°, 4e-4 (0.6 m). Over the São Francisco basin's span, about
  7° S to 21° S with the frame at its centre, it reaches 2 % to 4 %, which only
  the outline tolerance sees (5.6).

Checked figures (WGS 84 ellipsoid; SIRGAS 2000's GRS 80 differs in the
inverse flattening by 1.5e-6): `kx` 111,319.4908 and `ky` 110,574.2758 m per
degree at 0°; 105,009.8532 and 110,697.4248 at 19.4436° S; 55,800.0016 and
111,412.2875 at 60°. Times the fixture's spacing at 19.4436° S they give
28.2996 m and 29.8324 m, equal to `Geod.inv` over one cell to the four
decimals printed.

### 5.3 Constants: cells where a cell is meant, metres where a distance is meant

A cell, on a grid whose cells are not square in metres, is its **longer side
in metres** at the latitude in question: `LocalFrame.cell_m`. On DTM10 it is
10 m, so every product below equals today's figure there.

| Constant | Today | Becomes | Why |
|---|---|---|---|
| `U_FLOOR_M` | 30 m | `U_FLOOR_CELLS = 3` cells | Increment 29: "30 m (three cells) is the floor the DEM's resolution gives" (`docs/increments/29-nve-reference-catchments.md@f09b9863:1487`). |
| `Reach.corridor` | 30 m | `place` sets it to 3 cells; the model's default stays 30.0 for hand-built reaches | The burn's half-width across the line: a cross-section of at least seven nodes whatever the DEM. |
| `LAKE_GAP_M` | 30 m | removed; `seed_for` passes `gap = U_FLOOR_CELLS * cell` | Increment 29 defined it as "U's floor". |
| `WINDOW_MARGIN_M` | 2000 m | stays 2000 m, converted with `LocalFrame.units` at the seed | A guess at how far the catchment reaches beyond the seed: a distance, not a resolution. Keeping it in metres keeps every projected run's windows as they are. |
| `DOWN_EXTRA_M`, `CONFLUENCE_M` | 100 m | stay metres | Measured on NVE's mapped rivers (29's 19 of 139 confluences within 100 m of `P`): distances on the map. |
| `END_CAP_M` | 500 m | stays metres | How far a lowered end may run down the valley floor: terrain, not resolution. |
| direction-check skip | 100 m | stays metres, named `DIRECTION_MIN_M` | A chain too short to judge its direction: a length. |
| `JOIN_M` | 1 m | stays metres | The river file's digitising precision. |
| `map_radius`, `reach_up` | 500, 1000 m | stay metres (options) | Distances the person asks for. |
| `DROP_M`, `DIRECTION_M` | 0.001, 2 m | unchanged | Vertical. |
| `SLACK_M`, `_SLACK_M` | 1e-6 m | unchanged | Numerical slack on metre coordinates in the frame. |

All metre figures are applied in metres: in the river file's CRS when it is
projected, in a `LocalFrame` when it is in degrees (5.5), and converted to the
DEM's units with `LocalFrame.units` where they bound a window.

Scale each assumes: the metre constants were set on NVE's 1:50,000 rivers
over DTM10 (increment 29, NVE's 140 reference stations); the
cell constants were checked here at 10 m (DTM10 equality) and on the 28-30 m
ANADEM fixture (section 3). Neither is checked at basin scale; neither
affects the flood itself.

### 5.4 Areas

New in `geodesy.py`:

```python
def cell_areas_by_row(meta: RasterMeta) -> npt.NDArray[np.float64]:
    """m² of one cell in each row: geodesic (the DEM CRS's ellipsoid, one
    Geod polygon per row, as increment 26b's row_areas) when the grid is in
    degrees; dx * dy in every row otherwise."""

def mean_cell_area(mask: npt.NDArray[np.uint8], meta: RasterMeta) -> float:
    """m² per in-node of `mask`: dx * dy exactly on a projected grid; on a grid
    in degrees, the row-weighted sum over the in-nodes divided by their count."""

def area_m2(geometry: BaseGeometry, crs: str) -> float:
    """geometry.area on a projected CRS; |Geod.geometry_area_perimeter| (the
    CRS's ellipsoid) on one in degrees."""
```

Measured: 50,400 rows (14° at 1") in 0.12 s by `Geod`, and a closed form by
authalic latitude agrees to 2.9e-10; the geodesic area of a 100,000-vertex
ring takes 32 ms. `land_cover_fractions.row_areas` (26b, WGS 84 hard-coded,
the audit's row B4) is left as it is; this function does not import it.

Where each is used:

- `Catchment` gains `nodes_area: float` (m², `nodes * mean_cell_area(mask, m)`)
  and `reduced_area: float` (m², `area_m2(reduced, crs)`). `fine_area`,
  `holes_area` and `lake_area` become `area_m2` of the fine outline, of each
  filled hole's ring, and of the lake. On a projected DEM each is the same
  float it is today.
- The sensitivity (`assess`) gets `mean_cell_area(out.mask, m) / 1e6` as its
  `cell_area_km2`, in place of `m.delta_x * m.delta_y / 1e6`. On a projected
  DEM that is the same expression. In degrees it is exact for the area at the
  placed node (that node's catchment is the mask) and, for the few nodes up
  and down the chain, off only by the latitude spread of the area they add or
  lose, which is small beside the 5 % swing bar.
- The command's report reads `nodes_area`, `fine_area`, `reduced_area` and
  `holes_area` from the result; the file's `fine_area_m2` and
  `reduced_area_m2` are those fields. Accumulation stays a node count.

### 5.5 The gauge and the burn

- **River file in any CRS** (row A9). `catchment --rivers` no longer refuses a
  river file whose CRS is not the DEM's: the reach is placed in the river
  file's CRS and moved into the DEM's in `delineate`. `station-catchments`
  keeps the same-CRS rule (`check_reach_crs`), since it compares with
  reference polygons in the river file's CRS.
- **Placement** (`gauge.place`) gains a keyword `cell_m: float`, required: the
  DEM's cell at the station (`LocalFrame.at(dem_crs, station).cell_m(meta)`,
  `meta` the first tile's). `U = min(max(d, U_FLOOR_CELLS * cell_m), map_radius)`
  and `corridor = U_FLOOR_CELLS * cell_m`. `seed_for` (in
  `catchment_batch.py`) takes `cell_m` and passes it on, and when the river
  file's CRS is in degrees it places in a `LocalFrame` at the station: the
  segments and the station moved forward, `place` run, the reach line and
  `position` moved back by `inverse`. Distances it reports stay metres. Lakes
  with a river file in degrees do not occur (only the batch passes lakes, and
  it needs a DEM in metres) and are a `ValueError`.
- **The reach in the DEM's CRS** (`_gauged`). The line is moved from the
  request's CRS into the DEM's (pyproj, refused in plain words when a vertex
  has no image). The burn frame is `LocalFrame.at(dem_crs, centre of the
  line's box)`. When the river file's CRS is not the DEM's, or the DEM is in
  degrees, `reach.at` is measured again in the burn frame: the mapped
  position `P`, interpolated at `at` along the line in the request's CRS (in
  its frame when in degrees), moved into the burn frame and projected onto
  the moved line. Otherwise (today's case) the reach is used as given. The
  downstream length is measured on the line in the burn frame. The window's
  pad is `reach.corridor + WINDOW_MARGIN_M` metres, converted by
  `LocalFrame.units`.
- **The burn** (`burn.burn_reach(window, reach, frame=None)`). With a frame it
  burns on `frame.meta(window.meta)` and `frame.forward(reach.line)`: every
  metre in `burn.py` is then a metre. The returned tile keeps `window.meta`.
  The gap refusal names the node in the DEM's CRS (`window.meta.node_xy`),
  with five decimals in degrees and none in metres. `None` is today's
  behaviour exactly.

### 5.6 The outline and the output

- The trace runs on the native grid as today; rings come out in the DEM's
  CRS. Orientation (the shoelace sign) is the same in degrees.
- **The reduction** runs in metres: `_reduced` maps the fine ring and the
  seed through `LocalFrame.at(dem_crs, centre of the fine outline's box)`,
  takes the origin as the frame ring's minimum, reduces there with
  `_core.reduce_ring` at the tolerance in metres, and maps back with
  `inverse`. The default tolerance is `2 * frame.cell_m(m)` (today's
  `2 * max(dx, dy)` on a projected DEM).
- **The output CRS is the DEM's** (decided; question Q1): the GeoJSON names
  EPSG:4674 for ANADEM, with longitude-latitude pairs, the project's
  `always_xy` order (`src_python/tin_engine/crs.py@f09b9863:4-6`).
  `mesh --domain` reprojects a domain from its own CRS
  (`src_python/tin_engine/domain.py@f09b9863:41-61`), so the file feeds
  `mesh --dem <ANADEM> --out-crs <suggested>` as it is, and QGIS opens it.
- **The report** says degrees in plain words on a DEM in degrees: "window 1:
  longitude -44.13472 to -44.04876, latitude -19.48648 to -19.40051, 320 x
  320 nodes, …" (five decimals); metres as today otherwise. The start point
  prints as two plain numbers (the `np.float64(…)` defect is fixed for both).
  Areas in km² and the tolerance in metres, as today.

### 5.7 Refusals

- `station-catchments` on a DEM in degrees: "station-catchments needs a DEM
  in metres (a projected CRS); for a DEM in degrees, run catchment for each
  station" (question Q3).
- A seed, lake vertex or reach vertex with no image in the DEM's CRS: as
  today, in plain words.
- Everything else (cut by the edge, NoData, memory) as today; the window
  geometry is index-based and unchanged.

### 5.8 What does not change

The flood, the accumulation, the window loop's logic, the mosaic, the trace,
`reduce_ring`, `mesh` and its resampling path. On a projected DEM the
`catchment` file is byte-identical to today's (test R11); only the stderr
start line loses its `np.float64(…)` wrapping.

## 6. Left out, with reasons

- **Whole-basin memory.** The São Francisco basin's box at ANADEM's spacing
  is about 51,900 x 42,700 nodes, 2.2e9: under the core's 2^32 limit, about
  13 GB for an ungauged flood (6 bytes a node) and about 44 GB for a gauged
  one (20 bytes a node). The gauged case is over half of a 32 GB machine.
  Increment 38 makes the budget a setting; a flood in tiles (Barnes 2016,
  "Parallel Priority-Flood depression filling for trillion cell digital
  elevation models", Computers & Geosciences 96) is the route past it. Not
  this increment.
- **`station-catchments` in degrees** (refused, 5.7): its reference
  comparison (`reference.py`) and its fine-area column assume metres.
- **Tiles in two CRSs** (row A2) and **the antimeridian and poles** (row A3):
  unchanged; the window is one grid.
- **26b's `row_areas` on WGS 84** (row B4): its own fix.
- **Weighted accumulation** (area in place of a count, in C++): not needed
  while the sensitivity reads areas near one node (5.4).

## 7. Red suites (`@tester`)

Fixture figures are from section 3 unless stated. The layering table
(`tests/python/test_layering.py`) gains `geodesy` in L1 with
`"crs io.models"`, and `burn` gains `geodesy`; `catchment` and
`catchment_batch` gain `geodesy`.

Tests that encode today's refusals and are reversed or replaced:
`tests/python/test_cli_mesh_geographic.py@f09b9863:463-477` (the catchment
refusal), `tests/python/test_catchment.py@f09b9863:522-536` (a reach in
another CRS), `tests/python/test_cli_catchment.py@f09b9863:715` (a river file
in another CRS); `gauge.LAKE_GAP_M == 30.0`
(`tests/python/test_gauge.py@f09b9863:565`) goes with the constant. Existing
calls of `gauge.place` gain `cell_m=10.0`, which keeps their pinned figures
(30 m floor, 130 m downstream).

1. **R1 The refusal is gone.** `catchment` on `anadem_velhas.tif` with seed
   597180 7852380 in EPSG:31983 exits 0 and writes one feature in EPSG:4674;
   `nodes` 12,725.
2. **R2 The Velhas catchment against the reference.** The written outline,
   moved to EPSG:31983, shares at least 0.98 of its area with
   `tests/fixtures/velhas/catchment.geojson` and the reference shares at
   least 0.98 of its with it (fine outline measured 0.9864 and 0.9903; the
   tester measures the reduced one on the green tree and pins its figure
   with a margin). `fine_area_m2` within 1e-3 relative of 10.7441e6 m².
3. **R3 Same array, two grids, same catchment.** A bowl from
   `catchment_fixtures.py` flooded on its projected meta and on a geographic
   meta (EPSG:4326, 1" spacing, at 0° and at 60° N): identical masks and
   node counts. This is the test that the flood ignores coordinates.
4. **R4 Geodesic area against a hand figure.** On a 1" grid in EPSG:4326, a
   catchment of known rows: `nodes_area` equals the sum over rows of the
   count times the band area from the closed form by authalic latitude,
   written out in the test (not `Geod`), to 1e-6 relative. Hand figures: one
   1" cell is 949.7741 m² at the equator and 479.6898 m² at 60° N.
5. **R5 LocalFrame at two latitudes.** `units(2000.0)` on EPSG:4326 is
   (0.01796631°, 0.01808739°) at 0° and (0.03584229°, 0.01795134°) at 60°;
   `cell_m` on a 1" grid is 30.9221 m at 0° and 30.9479 m at 60° (the
   north-south side, the longer). Identity on a projected CRS: `units(2000.0)`
   is (2000.0, 2000.0), `meta(m) is m`, `forward(xy)` equals `xy` bit for bit.
   `forward` then `inverse` returns the input to 1e-9°.
6. **R6 Cell constants.** `place(..., cell_m=c)` gives `U = 3c` when the
   station is nearer than that, and `corridor = 3c`, at c = 10 and c = 29.8324.
7. **R7 A gauge on a geographic twin.** The 10 m valley of `gauge_fixtures.py`
   rebuilt on a geographic grid whose spacing is 10 m at its latitude
   (`dx = 10/kx`, `dy = 10/ky`), the river line moved likewise: the same chain
   (rows, columns), the same placed node and the same node count as on the
   projected valley; areas within 1e-3.
8. **R8 Rivers in another CRS.** The projected valley with its river file in
   EPSG:4326 gives the projected run's chain and node count; the geographic
   twin with its river file in the projected CRS likewise.
9. **R9 The gap refusal in degrees** names the node with five decimals.
10. **R10 `station-catchments` on a DEM in degrees** is refused in the words
    of 5.7, exit 2, nothing written.
11. **R11 Projected runs are byte-identical.** `catchment` on the DTM10
    fixture (`tests/fixtures/dem_archive/7908_3_10m_z33.tif`, a seed the
    tester picks whose catchment is contained, stated in the test), and
    `catchment --rivers` on the valley of `gauge_fixtures.py` as
    `tests/python/test_cli_catchment.py` runs it: each file's SHA-256 equals
    the figure the tester measures on the red commit's tree, where the code is
    master's.
12. **R12 The file feeds `mesh`.** `mesh --dem anadem_velhas.tif --out-crs
    <the suggested CRS> --domain <R1's file>` exits 0 (the geographic mesh
    suite's own `--out-crs` and tolerance).
13. **R13 The report in degrees** prints "longitude" and "latitude" with five
    decimals, areas in km² near 10.74, no `np.float64`.

No suite is invariant-critical in the sense of `docs/increments/README.md`
(no topology, no predicate): no mutation round is required.

## 8. Estimate, speed

Net production lines (`CLAUDE.md` §2), by file: `geodesy.py` about 80;
`raster.py` 8; `catchment_core.py` 1; `catchment.py` about 45; `burn.py` 8;
`gauge.py` 5; `catchment_batch.py` about 20; `cli.py` about 13. **About 180;
250 at the +39 % worst case. One PR, well under 700.** It touches no C++ file.

**Speed.** Hydrology only; no refine or mesh code, so no `@perf` acceptance
by the rule. The flood is the same core call as today (section 3: 12 ms at
0.1 M nodes, 11.4 s at 41 M). New work per run: one `Geod` call per row of the
final window (0.12 s at 50,400 rows) and geodesic areas of the outlines
(32 ms per 100,000 vertices); the frames are a few multiplications per
vertex. A projected run adds nothing measurable: `mean_cell_area` returns
`dx * dy` without touching the mask. An extrapolation, not a measurement:
at the 3.6 M nodes/s seen at 41 M nodes, a 2.2e9-node basin window would
take at least ten minutes to flood.

## 9. The ideas file

Taken up: *Catchments on a geographic DEM (GLO-30)*. Overlapped: *Product
review (b), generality audit* (its rulings, Q5); *Swedish SE2 river-system
catchments* (its open note that GLO-30 needs this fix). Left out and handed
back: *Neighbouring catchments simplified together*, *Catchment from a lake
seed behind a dike*, *Unbiased catchment height statistics under the
sup-norm*.

## 10. Questions for Ola (defaults taken while away)

- **Q1** Write the catchment in the DEM's CRS (degrees for ANADEM), which
  `mesh --domain` and QGIS read? Default: yes. The other choice was the
  suggested projected CRS, which `mesh --out-crs` would then have to match.
- **Q2** Which constants become cell counts: the gauge's uncertainty floor,
  the burn corridor and the lake gap (three cells each); the window margin,
  confluence and downstream distances, the end cap, map radius and reach
  length stay metres (5.3)? Default: as listed.
- **Q3** `station-catchments` on a DEM in degrees refused in plain words, its
  support a later row? Default: yes.
- **Q4** On a cell that is not square in metres, "a cell" is the longer side?
  Default: yes (the safe side for a floor).
- **Q5** `catchment --rivers` accepts a river file in any CRS, while
  `station-catchments` keeps the same-CRS rule? Default: yes.

## Review

**Design review, round 1, 2026-10-12.** Range `f09b9863..9b433a11`, 0 production lines, @reviewer: CHANGES REQUESTED, three blockers. (1) /Users/skavhaug/projects/rasputin/.claude/worktrees/geo-catchment/docs/increments/40-geographic-catchment.md@9b433a11:64-67: the GRASS bullet is false; `r.watershed/ram/do_cum.c`'s `get_dist` converts a lat/lon grid's spacing to metres with `G_distance`, one average over the window, not per row. (2) /Users/skavhaug/projects/rasputin/.claude/worktrees/geo-catchment/docs/increments/40-geographic-catchment.md@9b433a11:482: 479.6898 m² is the band 60°00'00" to 60°00'01"; a cell centred on a node at 60° N (the probe's `row_area`) is 479.6918 m², a 4e-6 difference that fails R4's 1e-6 tolerance; state the cell's edges or use 479.6918. (3) /Users/skavhaug/projects/rasputin/.claude/worktrees/geo-catchment/docs/increments/40-geographic-catchment.md@9b433a11:157-160 and :297: the margin sites omit `src_python/tin_engine/catchment.py@f09b9863:257`, `:279` and `:290` (the growth loop and its doubling), and the design does not say how the pair `units()` returns passes through `_plan` and `_joined`, which take one float today; name the sites and the margin's shape. Suggestion: :366 should say the station is moved into the DEM's CRS before `LocalFrame.at`. Verified: flood, upstream and accumulate read only rows, cols, size, values and NoData, nothing weights the diagonal; both probes reproduce with master's `_core` (12,725 nodes, 10.7445 and 10.7441 km², overlaps 0.9903 and 0.9864, 11.5 s at 41 M nodes); the `kx`/`ky`, `units` and `cell_m` figures; R11 can fail; no degree, CRS or path reaches C++; the TauDEM and WhiteboxTools claims; increment 38's overlap limited to `_grow` and `_gauged` signatures; check_citations exit 0. Lesson: a claim that a tool's manual says nothing needs a grep over the whole tool.
