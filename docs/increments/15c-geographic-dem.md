# Increment 15c — a geographic DEM, projected onto the TIN's CRS, resampled, and checked against its source

Status: **designed by `@architect`, 2026-10-01; not implemented.** Written
before `@tester`, per `docs/increments/README.md` step 1. It replaces R8 and
R10 of `docs/increments/15-dem-mosaic.md` for 15c, following Ola's Q6 and Q9
rulings of 2026-09-30. Questions Q11-Q17 at the end are Ola's.

## Why a file of its own

`15-dem-mosaic.md` is 1,575 lines and records four sub-increments, two of them
shipped. Its R8 (the lattice frame) and R10 (the output CRS) are superseded by
Q6, and its 15c test list, LOC table and acceptance were written for the
superseded design. Rewriting them in place would leave a reader to work out
which half of a 1,600-line file still holds. So 15c gets this file, and
`15-dem-mosaic.md` keeps the record (rulings, B1-B9, R1-R7, R9's 15b half)
with a pointer here. What 15c still takes from it, by section: B1-B9
(ANADEM's facts), R4-R5 (planning and assembling a mosaic, which the source
side reuses unchanged), R7 (the memory cap), R9 (the one reprojection site),
R12 (what domain decomposition needs). 15d stays designed there.

## What is ruled, and what was measured

**Ruled by Ola, 2026-09-30** (`15-dem-mosaic.md`, "Ruled by Ola"):

- **Q6.** A geographic DEM is projected onto the target TIN CRS and
  resampled onto a square grid there, in parallel. After meshing that grid,
  every source DEM node, projected into the target CRS, is checked against
  the mesh, and a node still outside `--tolerance` is inserted as an
  off-node vertex. The guarantee is against the source DEM, not only the
  resampled grid.
- **Q9.** "32 GB is fine": a dense canvas, no block-sparse raster. Two
  basin-sized canvases do not fit, so the source is read by window into the
  one target canvas.
- **Q10.** ANADEM extracts may be committed, crediting ANADEM and Copernicus.
- **Open:** Q7 (who chooses the target CRS), the basin tolerance, and whether
  commercial use matters for the DEM choice.

**Measured on 2026-10-01** (`docs/benchmarks/2026-10-01/basin-piece/README.md`,
by `@perf`, on Copernicus GLO-30 because ANADEM's host answered 403). The
piece is BHO ottobasin 76949 (upper Rio das Velhas, 11,667.6 km², 13.0 M
source nodes), resampled by hand to a 30 m grid in EPSG:31983 and meshed by
today's pipeline; the basin figures come from 200 random 10 km boxes.
What bears on this design:

- **Surprise 1.** Memory, not time, limits the whole basin: a floor of about
  77 GiB at 1 m, 41 GiB at 2 m, 19 GiB at 5 m and 13 GiB at 10 m in one
  process, 8.6 GiB canvas included. Refine itself is fast (2.4 min for the
  basin at 1 m, extrapolated).
- **Surprise 2.** The final check is a refinement phase, not a touch-up.
  After meshing the resampled grid, the share of interior source nodes over
  tolerance is 7.71 % at 1 m, 2.50 % at 2 m, 0.39 % at 5 m, 0.08 % at 10 m
  over the basin's boxes: 52 M, 17 M, 2.6 M and 0.53 M basin nodes, which is
  44 %, 30 %, 15 % and 8 % of the mesh's vertex count. On the steeper piece,
  22.30 % at 1 m and 0.31 % at 10 m. The resampled grid itself is off the
  source by a median 0.36 m, p99 3.14 m, max 27.8 m.
- **Surprise 3.** A strip along the domain's edge is never checked: a
  triangle whose closed node set is empty has error 0 and converges
  (`include/terrain/refinement/scan.hpp@434c374:56`; `needs_split` at
  `include/terrain/refinement/refine.hpp@434c374:187-189` asks for a node).
  On a 10 km box with four corners, source nodes within 30 m of the edge were
  off by up to 541.6 m at a 20 m tolerance; on the piece's dense BHO outline,
  17.9-38.5 m over the six tolerances. It holds for a projected DEM meshed
  directly too, measured against the DEM's bilinear surface between nodes.

## Prior art: legacy and literature

### Literature

- **Greedy insertion, unchanged in method.** Garland and Heckbert, "Fast
  polygonal approximation of terrains and height fields", CMU-CS-95-181,
  1995, as in increments 14 and 14b. Phase 1 (D1) is exactly that, on the
  resampled grid. Phase 2, the final check, is greedy insertion over
  **scattered** samples rather than a grid, which is the older form of the
  method: De Floriani, Falcidieno and Pienovi, "Delaunay-based representation
  of surfaces defined over arbitrarily shaped domains", *Computer Vision,
  Graphics, and Image Processing* 32:127-140, 1985; Heller, "Triangulation
  algorithms for adaptive terrain modeling", *Proc. 4th Int. Symp. on Spatial
  Data Handling*, 1990. Both recalled, not reread; Garland and Heckbert
  survey them. **What differs:** two sample sets in sequence. The grid is
  meshed first because increments 14-21 make that fast; the scattered source
  nodes then repair what the grid got wrong. The guarantee is the scattered
  method's (every source sample within tolerance), not the grid's.
- **Raster reprojection.** GDAL's `gdalwarp` maps each destination pixel back
  through the inverse transform, approximated by linear interpolation along
  rows within an error threshold (0.125 pixel by default), and resamples with
  a chosen kernel (nearest by default). GDAL documentation, recalled; GDAL is
  prohibited here, so this is reference only. **What differs:** every target
  node gets an exact inverse transform from pyproj, no approximation: the
  basin-piece run did 30.9 M of them, with bilinear lookups, in 0.99 s on 8
  threads (`basin-piece/README.md`, Method step 3). The kernel is bilinear in
  the source's own index space, because increment 16's R0 says DEM values are
  point heights with bilinear between them; any other kernel would invent a
  surface the rest of the pipeline does not assume.
- **Resampling changes terrain.** Usery et al., *J. Geographical Systems*
  6:289-306, 2004; Kienzle, *Transactions in GIS* 8(1):83-111, 2004 (both
  recalled; they study slope and flow). Here it was measured: B6 and
  Surprise 2. The final check exists because of it.
- **Checking a simplified surface against its original samples.** Metro:
  Cignoni, Rocchini and Scopigno, "Metro: measuring error on simplified
  surfaces", *Computer Graphics Forum* 17(2):167-174, 1998 (recalled), which
  measures the distance between two surfaces by sampling. **What differs:**
  the check here is vertical, at the source's own nodes, and repairs as well
  as measures.
- **Bucketing points by a uniform grid** for the per-triangle query of D4:
  Franklin, "Uniform grids: a technique for intersection detection on serial
  and parallel machines", *Proc. Auto-Carto 9*, 1989; Akman, Franklin,
  Kankanhalli and Narayanaswami, "Geometric computing and the uniform grid
  data technique", *Computer-Aided Design* 21(7):410-420, 1989 (both
  recalled). The buckets here are the target grid's own cells, which holds
  about one source node each at the source's resolution.
- **Projections.** Snyder, *Map Projections — A Working Manual*, USGS
  Professional Paper 1395, 1987 (recalled). LCC and UTM are conformal: angles
  are true locally, so a triangle judged Delaunay in the target CRS is judged
  on the ground to first order, which is what Q6 bought over the lattice
  frame's up-to-6.5 % stretch.

**Novelty: none claimed.** The two phases are a composition of known
methods. One sentence this design makes true could later look like a claim:
"every source node of a geographic DEM inside the domain lies within the
tolerance of a TIN built in a projected CRS". It is stated as a property of
this pipeline (J2), not as new. Before anyone calls it new, search Google
Scholar and IEEE Xplore for "TIN reprojected DEM error original samples",
"terrain simplification geographic grid projected mesh guarantee" and
"greedy insertion scattered points after grid refinement". No web search tool
was available to this round.

### Legacy

```sh
$ git grep -l -iE 'resampl|reproject|warp|Transformer' legacy-archive -- legacy
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/avalanche.py
legacy-archive:legacy/rasputin/geometry.py
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/tests/test_gml_repository.py
legacy-archive:legacy/tests/test_mesh.py
$ git grep -n -iE 'resampl|warp' legacy-archive -- legacy
$ git grep -l -iE 'tolerance|max_error|error_bound' legacy-archive -- legacy
```

The last two return nothing: the legacy never resampled a raster and never
measured a mesh against its DEM. The hits of the first are point-by-point
vector transforms (`15-dem-mosaic.md`, Legacy, already reported them) and one
mesh transform: `legacy/rasputin/application.py:101-123` meshed in the
raster's CRS, simplified by an edge-collapse `ratio` (no tolerance), and
transformed the finished vertices into `target_coordinate_system`. That is
the order Q6 rejected (mesh in the DEM's frame, transform out), so **nothing
is carried**. No legacy constant is re-derived, so no `@migration-expert`
pass is needed.

## Scope

After 15c,

```sh
rasputin mesh --dem anadem_23K.tif --domain bho_76949.geojson \
    --out-crs EPSG:31983 --tolerance 5 --out velhas.vtk
```

reads a geographic DEM (one file or a directory of tiles in one geographic
CRS), resamples it onto a square grid in the target CRS, meshes that grid as
today, then checks every source node inside the domain and inserts those still
over 5 m. The mesh is written in the target CRS.

- **In:** geographic DEM tiles (EPSG geographic 2D, degrees) through the
  existing reader and mosaic; the target grid; parallel resampling; the
  source nodes as check points; the final check in C++; the CLI and the
  recorded fields; the same path for a projected DEM whose CRS differs from
  `--out-crs` (Q12).
- **Unchanged:** a projected DEM without `--out-crs`, or with `--out-crs`
  equal to its own CRS, is meshed bit-identically to today (J1).
  `include/terrain/refinement/refine.hpp` is not edited.
- **Out (to 15d):** decoding only the needed windows of a tile, and a
  run on the whole basin. 15c is accepted on the Velhas piece (13.0 M
  source nodes), and refuses a request whose working set passes the memory cap (D6).
- **Out (placed below):** the edge strip for projected DEMs, and constraint
  feet in the final check.

## The design

### D1. Data flow

Two phases, and everything that knows a CRS stays above the line.

```
cli.mesh  --dem --domain --out-crs --tolerance      (parses flags, nothing else)
   |
   v  DemRequest(sources, domain, target_crs)
dem_input.open_dem                                   [Python, orchestration]
   |  footprints (headers only)
   |  DEM CRS == target CRS, or no --out-crs  --> today's 15a/15b path, untouched (J1)
   |  otherwise, the reprojected path:
   |    domain.to_crs(target)                        crs.reprojector, the one site
   |    target_grid_for(domain, target, spacing)  -> TargetGrid        (D2, pure)
   |    source_region(grid, dem crs)              -> Bounds, needed    (D2, pure)
   |    plan_mosaic + assemble                    -> source DemTile    (15a, unchanged)
   |    resample(grid, TileWindows(source))       -> target DemTile    (D3, threads)
   |    check_point_blocks(grid, TileWindows(source), domain)          (D4, a lazy iterator)
   v  DemInput(tile = target tile, ..., checks = that iterator)
cli._dem_mesh -> final_check.run                     [Python, orchestration]
   |  phase 1: refine(to_core(target tile), start, ...)   C++, unchanged
   |  drop the target tile
   |  CheckPoints(grid geometry); add(block) for each block; freeze()
   |  phase 2: refine_points(check points, phase-1 mesh, tolerance)  C++, new (D5)
   v  vertices in the target CRS, metres  -->  .vtk (D7)
====================================== _core boundary ======================================
C++ sees: a square grid in metres (RasterView), a start mesh in metres, and check
points as (x, y) in metres with z. It never sees degrees, a CRS, or a path.
```

Phase 1 is the existing pipeline on a projected grid; nothing in it learns
that the grid was resampled. Phase 2 is new: greedy insertion against the
source nodes, starting from phase 1's mesh. The two phases share the target
grid's lattice frame, so phase 2's Delaunay test is the same one phase 1 used.

The source is read through a `SourceWindows` protocol
(`window(r0, r1, c0, c1) -> array` and `meta`), never as a whole array. In
15c its one implementation, `TileWindows`, slices the assembled source tile;
15d replaces it by decoding only the blocks a window needs. That is how Q9's
"one target canvas" is reached: the design never asks for the source as a
whole, and 15d is a change of provider, not of algorithm.

### D2. The target grid

`TargetGrid` (`target_grid.py`, frozen Pydantic):

| field | meaning |
|---|---|
| `crs: str` | the target CRS, as `crs_label` writes it (`EPSG:n`, or one-line text) |
| `spacing: int` | `h`, whole metres, 1 to 1000 |
| `row0, col0: int` | the global index of the upper-left node |
| `rows, cols: int` | the grid's size in nodes |

- **Global node identity.** Node `(R, K)` sits at `x = K·h`, `y = −R·h` in
  the target CRS. `x_min = col0·h` and `y_max = −row0·h`. With `h` a whole
  number of metres and `|K|, |R| < 2²⁵`, `K·h` is an exact double, so the same
  node has the same coordinates in every run, window and subdomain
  (`15-dem-mosaic.md` R12). 21b's integer incircle needs `dx == dy` and few
  significant bits; whole metres give both.
- **Spacing.** Default: the source's north-south node spacing in metres at
  the domain's centroid, rounded to the nearest whole metre (ANADEM 29.8 m
  gives 30; GLO-30 1″ gives 31). The Python API takes `spacing` explicitly;
  the CLI does not expose it (Q16). B6 measured that a finer grid does not
  remove the resampling error, and phase 2 removes it anyway, so the default
  aims at the source's own density.
- **Extent.** The domain in the target CRS, grown by the cell diagonal
  `√2·h` with mitred corners (15b's needed region, the nodes bilinear z can
  read), its bounds snapped outward to multiples of `h`. With `--bbox`
  instead of `--domain`, the box is read in the target CRS, the frame the mesh
  is built in, and snapped the same way.
- **Source region** (`source_region`): the target box's boundary densified
  to one point per target cell, transformed into the DEM's CRS, bounded, and
  grown by two source cells; the needed polygon is the grown domain
  transformed the same way. Both go to 15a's `plan_mosaic` unchanged, which
  selects the tiles, resolves overlaps (R5) and refuses gaps (I4).
- **Refusals:** a source region crossing ±180° longitude or reaching a pole;
  a domain vertex with no image in the target CRS (pyproj's `inf`); a target
  grid over the memory cap (D6).

### D3. Resampling

`resample(grid, source: SourceWindows, threads) -> DemTile`
(`target_grid.py`):

- One canvas of `grid.rows × grid.cols`, in the source's dtype, filled in
  blocks of 256 target rows on a pool of `threads` workers. Per block: the
  nodes' `(x, y)`, one inverse transform to the DEM's CRS
  (`crs.reprojector(target, source)`), fractional source indices
  `col = (lon − x_min)/dx`, `row = (y_max − lat)/dy`, and bilinear over the
  2 × 2 stencil read from `source.window(...)`. This is the measured
  prototype's shape (`basin-piece/prep_dem.py`, `resample`).
- **NoData rule**, the same as `vertex_z` in `scan.hpp`: a target node is
  NoData when any corner of its stencil is NoData or outside the source,
  whatever its weight. The target tile takes the source's sentinel, or NaN
  with `nodata=None` when the source has none.
- **Deterministic:** each node is computed from its own coordinates alone,
  so block size and thread count cannot change a value.
- The canvas becomes the tile through `DemTile._adopt` (15a R7, no copy). Its
  docstring names `assemble` as the one caller; it gains `resample` as the
  second, under the same rule (the buffer was allocated here and is never
  handed out writable).
- Its values are not trusted for the guarantee: phase 1 is measured against
  them, phase 2 against the source.

### D4. Check points: the source nodes, in the target CRS

**Producer, Python** (`target_grid.py`,
`check_point_blocks(grid, source, domain, threads)`): an iterator of
`(xy: float64 (N, 2), z: float32 (N,))` blocks, one per 512 × 512 block of
source nodes, in source row-major block order.

- A block whose outline, transformed into the target CRS, misses the domain
  is skipped whole (shapely, prepared). Every other block yields all its valid
  nodes transformed into the target CRS, minus non-finite images and points
  outside the target grid's node rectangle. Points outside the domain in a
  boundary block are kept: they are never inside a triangle, so phase 2 never
  scans them (D5), and testing 13 M points against a 7,310-vertex outline
  would cost more than holding the few that are outside.
- Check points come from the **assembled source mosaic**, never from tiles,
  so an overlap of two tiles yields each node once, with R5's value (J7).

**Store, C++** (`include/terrain/refinement/check_points.hpp`,
`class CheckPoints`):

- Built on a `raster::RasterGeometry` (the target grid's `x_min`, `y_max`,
  `h`, `h`, rows, cols): numbers only. `add(xy, z)` converts each point to
  fractional `(col, row)` in that lattice and files it in its cell
  `(floor(row), floor(col))` (the last cell for a point on the far edges).
- **Packing: 16 bytes per point**: the cell's column (`uint32`), the offsets
  inside the cell (two `float`), and z (`float`). Rows are implicit after
  sorting, with one `uint64` start offset per grid row. A `float` offset in
  [0, 1) resolves 2⁻²⁴ of a cell, under 2 µm at 30 m. **The check point is
  the stored position**: membership, error and the inserted vertex all use the
  same reconstructed `double`, so the 2 µm is a fixed displacement of where
  the source node is taken to be, not an inconsistency. z is `float`: exact
  for ANADEM and GLO-30, which are float32; a float64 source is rounded, and
  the run says so.
- `freeze()` sorts by `(cell row, cell col, row offset, col offset, z)`, so
  the store, and so phase 2's tie-breaks, do not depend on the order blocks
  were added. Two points at the same position keep the first and are counted
  as `duplicates` (J7 makes this 0 for a mosaic; it is counted so a change
  shows). `add` after `freeze` is a refusal. Afterwards the store is
  immutable and any number of threads may read it.
- Query: `for_each_in(row, c0, c1, f)`, by binary search on the column within
  the row's slice.

### D5. The final check: `refine_points`

`include/terrain/refinement/refine_points.hpp`, a new header beside
`refine.hpp`, which is not edited (J1):

```cpp
struct PointRefineOptions { double tolerance = 0.0; unsigned threads = 0; };

template <class Store>  // CheckPoints, or a test double with for_each_in
[[nodiscard]] RefineOutcome refine_points(
    const Store& points, const raster::RasterGeometry& grid,
    const IndexedMesh2& start, std::span<const double> z, std::span<const std::uint8_t> valid,
    std::span<const std::array<std::uint32_t, 2>> edges, std::span<const std::uint32_t> masks,
    const PointRefineOptions& options);
```

`RefineOutcome` is reused. `max_error` is over check points; `coincident` and
`coincident_max_error` are added (below). The start is phase 1's output as
numbers: vertices in the target CRS, their z and validity, triangles,
constraint edges and masks.

1. **Lattice.** `detail::to_lattice(grid, start, ...)` as in `refine`: a
   vertex that is a node bit for bit is a node, any other is off-node. A
   vertex z table `zt` holds the given z, NaN where invalid. Phase 2 never
   reads the target grid's values: `grid` is only the frame and the buckets.
2. **Scan** (parallel, read-only, one result per triangle), `scan_points`:
   - for each cell row the triangle meets, the column range the triangle
     covers within that band (from its vertices in the band and its edges'
     crossings of the band's two lines), widened by one cell each side, then
     every stored point in those cells is tested **exactly**: three
     `orient_sign` calls (DefaultKernel on `(col, −row)`, the kernel
     `LatticeMesh` already uses) decide the closed triangle;
   - a point equal to one of the triangle's corners is skipped, as the grid
     scan skips a triangle's own vertices;
   - error `|z − plane(p)|`, the plane from `zt` and the same double
     barycentric expression the off-node grid scan uses; the worst point wins
     by strictly larger error, so ties go to the first in store order,
     whatever the thread count;
   - a triangle with a NaN corner (void) gets 14's rule: the check point
     nearest a void corner, and `uncovered` counts its points.
3. **Split** (serial, triangle-index order), as `refine` does: the worst point
   goes in with `split_inside` (a `MeshVertex` overload, +1 small function in
   `lattice_mesh.hpp`) or, on an edge (one zero orientation), `split_edge`,
   skipped this round if the neighbour across that edge was already touched.
   A constrained edge stays constrained on both halves (as for 20b's feet).
   Its z is the point's z, appended to `zt`. `legalise_around` with
   `lattice_frame(h, h, rows, cols)` follows, and every written slot is
   touched and rescanned next round.
4. **Stop** when no triangle needs a split. **Terminates**: every insertion
   is a stored point not yet a vertex (a vertex lies in a closed triangle only
   as one of its corners, and corners are skipped), and the store is finite.
5. **Coincident check points.** A stored point equal to a *start* vertex is
   never inserted and never scanned. After the loop, each start vertex looks up
   its own cell; a match is counted in `coincident` with its error
   `|z − zt|` in `coincident_max_error`. J2 is stated without these points, and
   the file records both numbers. For a projected source node to land on a
   start vertex to 2 µm is not expected to happen; it is counted so that it
   cannot happen silently.
6. **Output** as `refine`'s, with z from `zt` and `valid` where it is a
   number.

**No constraint feet and no start quality in phase 2.** 20b's foot needs a
height at the foot, and phase 2 has heights only at the check points. Whether
phase 2's insertions make needles near constraints is measured at acceptance
(the worst angle against phase 1's); feet for phase 2, with the point's own z
and 20b R5's re-check, would be a follow-up.

**Why a second loop and not a shared one.** `refine`'s loop and this one have
the same skeleton (parallel scan, serial index-order split, legalise, active
set). A loop templated on its scanner would be the cleaner end state, but it
means editing the performance-tuned loop that 21a and 21b tuned, with a
bit-identity proof and `@perf` acceptance for Norway. A separate header costs
about 60 duplicated lines and leaves Norway's path untouched by construction.
`detail::rebuild_active`, `detail::needs_split` and `detail::to_lattice` are
reused, not copied. The merge is named here as a later refactor.

### D6. Types and the I/O boundary

**What crosses into `_core`** (the boundary rule in `project_structure.md`'s
`raster` section, restated for 15c): every number is in a projected CRS, in
metres, on or relative to a square grid. Concretely: the target grid as a
`RasterView` (unchanged binding), the start mesh, and check points as
`(x, y)` in the target CRS with z. **No degrees, no CRS, no path.**

- `RasterMeta` (`io/models.py`) gains `geographic: bool = False` and
  `crs: str`, filled from `epsg` when not given (`EPSG:n`); `epsg` becomes
  `int | None` so a target CRS without an EPSG code (a fitted LCC, Q11 (b)) has
  a meta. Existing constructions with `epsg=n` are unchanged. Call sites that
  format `EPSG:{meta.epsg}` move to `meta.crs`.
- `io/geotiff.py` accepts a geographic 2D CRS from `GeographicTypeGeoKey`
  (2048) with degree units, as B8 found it would with the one refusal lifted.
  Geocentric, 3D, user-defined and non-degree geographic CRSs stay refused.
- **The gate:** `raster.to_core` refuses a tile with `meta.geographic`. It is
  the one adapter into a core raster (increment 12), so this also keeps a
  geographic DEM out of `rasputin catchment` (increment 22), which is not
  extended in 15c. Tested by planting a geographic tile.
- `crs.reprojector` stays the only `Transformer.from_crs` in `src_python/`
  (I8), used both ways: source to target (check points, the domain, the
  features) and target to source (resampling, the source region).
- `CheckPoints` and `refine_points` are bound in `bindings/core.cpp` only:
  `CheckPoints(x_min, y_max, spacing, rows, cols)`, `.add(xy, z)` (NumPy
  arrays, validated shapes and dtypes), `.freeze()`, `.size`, `.duplicates`;
  `refine_points(points, vertices, triangles, z, valid, edges, masks, *,
  tolerance, threads)`, releasing the GIL like `refine`. `_core.pyi` stubs
  both. `include/terrain/` stays free of pybind11.
- **Memory cap** (15a R7: half of physical memory), checked in `open_dem`
  before any pixel is read: source canvas + target canvas + 16 B per source
  node in the source box. The Velhas piece is about 0.75 GB by that sum
  (7,197 × 4,376 source nodes, 7,347 × 4,208 target nodes). The whole basin
  is refused in 15c by design; 15d brings the source canvas down to a window.
- `final_check.py` (new, ~40 lines): `run(tile, start, checks, tolerance,
  threads) -> (phase-1 outcome, phase-2 outcome)`, the two calls and the
  store's construction, so `cli.py` (1,542 lines) grows by the options and
  fields only. It drops its reference to the target tile before building the
  store; the caller must not hold one either, and `_dem_mesh` is restructured
  so it does not.

### D7. What the file records

- `crs`: the target CRS (`EPSG:n`, else one-line WKT2, ASCII).
- `source_crs`: the DEM's; `source_transform`: pyproj's description of the
  source-to-target transformation, so a datum shift is on record.
- `computation_grid`: e.g. `square 30 m grid in EPSG:31983, node (R, K) at
  (30 K, -30 R), resampled bilinear from EPSG:4326`.
- `elevation_source`, the sentence: phase 1 as today, its achieved max error
  labelled "against the resampled grid", then `checked against N source
  nodes: M inserted in P rounds, max error E m at source nodes, C coincident
  with a start vertex (max F m)`.
- `--stats` rows: `resample`, `check points: project`, `check points: store`,
  `final check: scan (parallel)`, `final check: split + flip (serial)`, and
  the source-node error beside the grid's.

The output is already in the target CRS, so 15-dem-mosaic.md R10's
transform-on-output and its `max_reprojection_z_error_estimate` field are not
built: nothing is transformed after meshing.

## What the final check costs (Surprise 2)

(to be written)

## The edge strip (Surprise 3), placed

(to be written)

## Invariants

(to be written)

## Degeneracy policy

(to be written)

## Not in scope

(to be written)

## PR split and LOC

(to be written)

## Tests for `@tester`

(to be written)

## Acceptance

(to be written)

## Questions for Ola

(to be written)
