# Increment 15c — a geographic DEM, projected onto the TIN's CRS, resampled, and checked against its source

Status: **designed by `@architect`, 2026-10-01; not implemented.** Written
before `@tester`, per `docs/increments/README.md` step 1. It replaces R8 and
R10 of `docs/increments/15-dem-mosaic.md` for 15c, following Ola's Q6 and Q9
rulings of 2026-09-30. **Q11-Q14, Q16 and Q17 were ruled by Ola on
2026-10-01, and the design below is amended to them; Q15 is open, and Ola's
direction on it supersedes this file's dense-canvas and memory reasoning at
basin scale** (marked where it occurs, "Superseded at basin scale").

## Why a file of its own

`15-dem-mosaic.md` is about 1,600 lines and records four sub-increments, two of them
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
  one target canvas. *Revisited by Ola on 2026-10-01 (Q15 below): no dense
  canvas at basin scale.*
- **Q10.** ANADEM extracts may be committed, crediting ANADEM and Copernicus.
- **Open then:** Q7 (who chooses the target CRS; ruled as Q11 below), the
  basin tolerance (Q15, still open), and whether commercial use matters for
  the DEM choice (still open).

**Ruled by Ola, 2026-10-01**, on this file's Q11-Q17 (the questions are kept
as asked at the end; each is marked with its ruling):

- **Q11: (a), a required `--out-crs`, with a suggestion in the refusal.**
  Ola: "perhaps we could automatically find a suggestion CRS based on where
  in the world we are? We could take the bbox and use this as a criterion.
  Something that distorts areas and angles the least, but no very rigorous
  analysis if that is time consuming." So the suggestion is chosen by the
  box's latitude and shape among four conformal projections, and the message
  prints its worst scale and areal distortion over the box (D8).
- **Q12: (a).** A projected DEM with a different `--out-crs` takes the same
  path: resample onto the target grid, then the final check.
- **Q13: (a).** Every source node still over tolerance is inserted.
- **Q14: (a), extended.** The edge strip gets its own increment right after
  15c and before 15d: check points where constraints cross grid lines, z
  bilinear, **plus the midpoint between neighbouring crossings**. The
  guarantee wording is the one proposed (Ola did not reword it), extended by
  the midpoints ("The edge strip", below).
- **Q15: open, reframed.** Memory must not decide the tolerance. Ola: "When
  things are tiled on file, and we refine, why are we talking about memory
  limitation? There must be easy fixes to this." Ola's direction:
  - decompose along constraint edges, which are never flipped: "Only two
    intersecting subdomains sharing the constraint is needed to resolve
    potential corrections on the boundary, as these edges are never
    flipped";
  - read the DEM in windows with no dense canvas (this revisits Q9);
  - a separate fetch step that caches downloaded COG tiles: "we should cache
    the tiles we download, and probably do this as a preprocessing step, or
    our performance will drop while we wait for download".

  These get their own basin-scale design, next, in an increment file of its
  own under `docs/increments/`. It is not designed here. This file marks
  where its dense canvas and its memory reasoning are superseded at basin
  scale; at the Velhas piece, 15c's acceptance size, they stand.
- **Q16: (a).** The grid spacing is automatic; no CLI option.
- **Q17: (a), with ANADEM's new source.** Small GLO-30 and ANADEM extracts
  are committed after a licence check. ANADEM's own host
  (`metadados.snirh.gov.br`) refuses with HTTP 403, and Ola confirmed it.
  OpenTopography hosts ANADEM openly, as one cloud-optimised GeoTIFF:
  `https://opentopography.s3.sdsc.edu/raster/ANADEM/ANADEM_be/anadem_v1_compressed_COG.tif`
  (70.9 GB; range reads answer 206; DOI 10.5069/G9736P4G). Also ruled: the
  BHO-derived outline (the union of 1,163 BHO 2017 5k catchments, Rio das
  Velhas, ottobasin 76949) is committed as a fixture, with ANA's credit and
  the CC BY 4.0 notice. ANA's metadata record for it was not read (its host
  answers 403); it is to be checked when reachable. BHO 2017 5k is marked
  superseded by BHAE in ANA's catalogue.

**Measured on 2026-10-01** (`docs/benchmarks/2026-10-01/basin-piece/README.md`,
by `@perf`, on Copernicus GLO-30 because ANADEM's host answered 403). The
piece is BHO ottobasin 76949 (upper Rio das Velhas, 11,667.6 km², 13.0 M
source nodes), resampled by hand to a 30 m grid in EPSG:31983 and meshed by
today's pipeline; the basin figures come from 200 random 10 km boxes.
What bears on this design:

- **Surprise 1.** Memory, not time, limits the whole basin: a floor of about
  77 GiB at 1 m, 41 GiB at 2 m, 19 GiB at 5 m and 13 GiB at 10 m in one
  process, 8.6 GiB canvas included. Refine itself is fast (2.4 min for the
  basin at 1 m, extrapolated). *Superseded at basin scale:* this is the
  floor of one process holding the whole basin, which Ola's Q15 direction
  (decomposition, windows, no dense canvas) removes. It stays the right
  figure for what 15c alone could do.
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

**Measured by `@architect`, 2026-10-01, for the Q17 ruling:** the header of
OpenTopography's ANADEM COG, by one range read of its first 256 KiB
(`curl -r 0-262143 <url>` answered 206) and a hand parse of the first IFD:

- BigTIFF; **187,853 × 262,864 nodes** (columns × rows), float32, Deflate,
  512 × 512 tiles, NoData `-9999` in tag 42113; further IFDs follow
  (overviews, not inspected).
- Spacing 0.00026949458523585647° in both axes, the same as the SNIRH tiles
  (`15-dem-mosaic.md` B2); tie point (−82.51654705336689°, 14.079475111062084°);
  `GTRasterTypeGeoKey` 1 (area-registered, as B2).
- **The CRS differs from the tiles':** `GTModelTypeGeoKey` 2,
  `GeographicTypeGeoKey` **4674 (SIRGAS 2000)**, not 4326; angular units 9102
  (degrees); and, beside the EPSG code, a citation ("SIRGAS 2000") and the
  ellipsoid's semi-major axis and inverse flattening (keys 2049, 2057, 2059).
  D6's reader takes the CRS from 2048 alone (checking 2054 for degrees) and
  ignores the other three, as today's projected path takes its CRS from 3072
  alone (`io/geotiff.py`, `_projected_epsg`); G1 has a case with this exact
  key set. With a SIRGAS 2000 target such as EPSG:31983, the source and
  target share a datum, so no datum shift enters the resampling.

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
  frame's up-to-6.5 % stretch. The same book is the source of Q11's
  suggestion (D8), as recalled: its guidance to choose a projection by the
  region's latitude and shape (azimuthal near a pole, cylindrical near the
  equator, conic for a mid-latitude east-west extent, transverse cylindrical
  for a north-south one) and its one-sixth rule for a conic's standard
  parallels. **What differs:** only conformal members of each family are
  offered (stereographic, Mercator, LCC, Transverse Mercator), because the
  Delaunay test is about angles; Snyder's equal-area alternatives are not.
  The distortion printed is PROJ's own projection factors (`proj_factors`,
  through pyproj's `Proj.get_factors`), not a derivation of ours.

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
rasputin mesh --dem anadem_velhas.tif --domain bho_76949.geojson \
    --out-crs EPSG:31983 --tolerance 5 --out velhas.vtk
```

reads a geographic DEM (one file or a directory of tiles in one geographic
CRS; here a local window cut from OpenTopography's ANADEM COG), resamples it
onto a square grid in the target CRS, meshes that grid as today, then checks
every source node inside the domain and inserts those still over 5 m. The
mesh is written in the target CRS. Without `--out-crs`, a geographic DEM is
refused, and the refusal suggests a CRS for the box (D8).

- **In:** geographic DEM tiles (EPSG geographic 2D, degrees) through the
  existing reader and mosaic; the target grid; parallel resampling; the
  source nodes as check points; the final check in C++; the CLI and the
  recorded fields; the same path for a projected DEM whose CRS differs from
  `--out-crs` (Q12).
- **Unchanged:** a projected DEM without `--out-crs`, or with `--out-crs`
  equal to its own CRS, is meshed bit-identically to today (J1).
  `include/terrain/refinement/refine.hpp` is not edited.
- **Out (to the basin-scale design that follows, Q15):** decoding only the
  needed windows of a tile, a fetch step that caches downloaded COG tiles,
  decomposition along constraint edges, no dense canvas, and a run on the
  whole basin. That design replaces 15d's (`15-dem-mosaic.md`) where they
  overlap; which parts of 15d survive is its call. 15c is accepted on the
  Velhas piece (13.0 M source nodes), and refuses a request whose working set
  passes the memory cap (D6).
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
   |  geographic DEM, no --out-crs  --> refused; message from crs.suggest_crs(box)  (D8)
   |  projected DEM: CRS == target CRS, or no --out-crs  --> today's 15a/15b path (J1)
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
the basin-scale design (Q15) replaces it by decoding only the blocks a window
needs. That is how Q9's "one target canvas" is reached: the design never asks
for the source as a whole, and the basin-scale design (Q15) is a change of
provider, not of algorithm.

*Superseded at basin scale (Q15, 2026-10-01):* the one dense target canvas is
not the basin's route; Ola's direction is windows with no dense canvas, and
subdomains split along constraint edges. What 15c builds carries over as
seams, not as the plan: `SourceWindows` is where cached, windowed decoding
plugs in; a `TargetGrid` is a sub-rectangle of one global lattice (J6, its
`row0`, `col0`), so a subdomain's grid is just a smaller `TargetGrid`, and
`resample`, `check_point_blocks` and `CheckPoints` work on it unchanged; and
J9, in two halves: phase 2 never flips a constraint edge, the property Ola's
decomposition rests on; and phase 2 does split constraint edges, which a
decomposition along constraints must reconcile, so both subdomains sharing a
constraint agree on its Steiner points. The basin-scale design decides the
rest.

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
  gives 30; GLO-30 1″ gives 31). For a projected source (Q12) it is the
  source's own north-south spacing in its units, which are metres (the reader
  refuses others), rounded the same way: DTM10 gives 10. The Python API takes
  `spacing` explicitly; the CLI does not expose it (Q16, ruled (a)). B6 measured that a finer grid does not
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

- One canvas of `grid.rows × grid.cols` (dense, for the piece; *superseded
  at basin scale* by windows, see D1 and Q15), in the source's dtype, filled in
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
  A point outside the node rectangle, or with a non-finite coordinate or z,
  is dropped and counted in `outside` (the Python producer already drops
  them, so a non-zero count in a run means the producer changed).
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
[[nodiscard]] PointRefineOutcome refine_points(
    const Store& points,  // its geometry() is the grid: the frame and the buckets
    const IndexedMesh2& start, std::span<const double> z, std::span<const std::uint8_t> valid,
    std::span<const std::array<std::uint32_t, 2>> edges, std::span<const std::uint32_t> masks,
    const PointRefineOptions& options);
```

It returns `PointRefineOutcome`, defined in `refine_points.hpp`: a
`RefineOutcome` (with `max_error` over check points) plus `coincident` and
`coincident_max_error` (below). `RefineOutcome` itself is unchanged, so
`refine.hpp` is not edited (J1). The grid is not a separate argument: it is
the store's own `geometry()`, so a store filed on one grid cannot be scanned
against another. The start is phase 1's output as
numbers: vertices in the target CRS, their z and validity, triangles,
constraint edges and masks.

1. **Lattice.** `detail::to_lattice(points.geometry(), start, ...)` as in `refine`: a
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
  `int | None` so a target CRS without an EPSG code has a meta: D8's
  suggestion is a PROJ string, and pasting it into `--out-crs` is the
  expected use. Existing constructions with `epsg=n` are unchanged. Call sites that
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
  arrays, validated shapes and dtypes), `.freeze()`, `.size`, `.duplicates`,
  `.outside`;
  `refine_points(points, vertices, triangles, z, valid, edges, masks, *,
  tolerance, threads)`, releasing the GIL like `refine`. `_core.pyi` stubs
  both. `include/terrain/` stays free of pybind11.
- **Memory cap** (15a R7: half of physical memory), checked in `open_dem`
  before any pixel is read: source canvas + target canvas + 16 B per source
  node in the source box. The Velhas piece is about 0.75 GB by that sum
  (7,197 × 4,376 source nodes, 7,347 × 4,208 target nodes). The whole basin
  is refused in 15c by design. *Superseded at basin scale (Q15):* the cap
  stays in 15c as a guard, but at basin scale memory is not to decide what
  can be meshed; the basin-scale design removes the canvases rather than
  shrinking one of them.
- `final_check.py` (new, ~40 lines): `run(tile, start, checks, tolerance,
  threads) -> (phase-1 outcome, phase-2 outcome)`, the two calls and the
  store's construction, so `cli.py` (about 1,500 lines) grows by the options and
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

### D8. The suggested CRS (Q11)

A geographic DEM without `--out-crs` is refused before any pixel is read
(J10), and the refusal suggests a CRS to copy. `crs.suggest_crs(box,
geographic_crs) -> CrsSuggestion` (`crs.py`, pure; frozen Pydantic: `proj:
str`, `family: str`, `max_scale_error: float`, `max_areal_error: float`).

- **The box:** the domain's bounds in the DEM's geographic CRS (the domain
  transformed once through `crs.reprojector`); with `--bbox` and no
  `--out-crs`, the box cannot be read (it is in the target CRS, D2), so the
  DEM's own footprint from the headers is used. `(W, E, S, N)` in degrees,
  centre `(λc, φc)`; the box's width on the ground is compared with its
  height as `(E − W)·cos φc` against `N − S`.
- **The family, in this order** (all conformal; thresholds settled here),
  each first written with unit scale:
  1. **Polar stereographic** when `|φc| ≥ 70°` or the box reaches `|φ| ≥ 80°`:
     `+proj=stere +lat_0=±90 +lon_0=λc +k_0=s`.
  2. **Transverse Mercator** when the box is at least as tall as it is wide:
     `+proj=tmerc +lat_0=0 +lon_0=λc +k=s`.
  3. **Mercator** when wide and `|φc| ≤ 15°`: `+proj=merc +lon_0=λc +k_0=s`
     (a scale factor in place of a true-scale latitude; PROJ takes either).
  4. **Lambert conformal conic** otherwise (wide, mid-latitude): standard
     parallels at `S + (N − S)/6` and `N − (N − S)/6` (Snyder's one-sixth
     rule), `+lat_0=φc +lon_0=λc +k_0=s`.
- **The scale, balanced, in every family.** One evaluation with `s = 1`
  (below) gives the smallest and largest point scale over the box, `kmin` and
  `kmax`; then `s = 2/(kmin + kmax)`, rounded to 6 decimals, and a second
  evaluation gives the figures reported. For a conformal map this puts the
  worst scale error at `(kmax − kmin)/(kmax + kmin)`, about half the one-sided
  figure, so it is the least a single scale factor can reach for the chosen
  projection. It is one line and one extra evaluation, so it is not the
  time-consuming analysis Ola's "no very rigorous analysis if that is time
  consuming" excused. It replaces three cheaper choices made in the first
  draft: TM with `k = 1` (one-sided error), Mercator true at the centre
  latitude, and stereographic true at the centre latitude. For the LCC it
  trims what the one-sixth rule leaves unbalanced over the box (0.92 % to
  0.80 % on the Europe-wide box). Moving the standard parallels, or the
  stereographic's centre off the pole, would be an optimisation over the box,
  which is the rigorous analysis that was excused; not done.

  Every angle is rounded to 0.1°, so a box moved by a few cells usually gets
  the same suggestion (a centre near a rounding boundary can flip one digit,
  and `s` follows the box). The datum is the DEM's own (`+ellps` of its
  geographic CRS, with no datum shift added): ANADEM's COG is SIRGAS 2000
  (GRS80), GLO-30 is WGS 84.
- **The distortion:** pyproj's `Proj(suggestion).get_factors` on a 21 × 21
  lon-lat grid over the box, edges included; `max_scale_error` is the largest
  `|k − 1|` over the meridional and parallel scales, `max_areal_error` the
  largest `|areal − 1|`. Areal distortion is reported, not minimised: a
  conformal map's areal scale is `k²`, so balancing `k` also roughly halves
  it. `Proj` is not a `Transformer`, so J5's one reprojection site is
  unchanged.
- **The message**, one refusal: the DEM is geographic, `--out-crs` is
  required, then `suggested for this box (<family>): --out-crs "<proj>"`, and
  `worst scale error <a> %, worst areal error <b> % over the box`.

**Measured with a prototype of exactly these rules** (pyproj 3.8.0, a
throwaway script, not committed), the worst scale and areal errors, balanced,
with the unit-scale figure beside them:

| box (W, E, S, N) | suggestion | scale | areal | scale at `s = 1` |
|---|---|---:|---:|---:|
| Velhas piece (−44.674, −43.459, −20.456, −18.457) | TM, `lon_0=-44.1 k=0.999972` | 0.003 % | 0.006 % | 0.006 % |
| São Francisco basin, approx. (−47.7, −36.3, −21.2, −7.0) | TM, `lon_0=-42.0 k=0.997542` | 0.25 % | 0.49 % | 0.49 % |
| Norway, approx. (4.5, 31.2, 57.9, 71.2) | TM, `lon_0=17.9 k=0.996173` | 0.38 % | 0.77 % | 0.77 % |
| Equatorial (−70, −50, −5, 5) | Mercator, `k_0=0.998106` | 0.19 % | 0.38 % | 0.38 % |
| Svalbard, approx. (10, 34, 76.4, 80.9) | polar stereographic, `k_0=0.98983` | 0.39 % | 0.78 % | 1.42 % |
| Europe-wide (−10, 30, 40, 60) | LCC, 43.3°/56.7°, `k_0=0.99881` | 0.80 % | 1.61 % | 0.92 % |

Angular distortion is under 3e-6 in every row: PROJ's factors confirm the
four are conformal. The prototype was run on both WGS 84 and GRS80; every
figure shown is the same on both. For the basin, the suggestion's 0.25 %
compares with 1.11 % (areal 2.24 %) for EPSG:31983 or 31984 over the same box,
by the same 21 × 21 evaluation: its central meridian sits in the box's middle
and its scale is balanced.
The two catchments-in-one-CRS argument for (a) still holds: the suggestion
depends on the box, so two catchments get one CRS only if someone passes the
same `--out-crs` to both.

## What the final check costs (Surprise 2)

**What is measured** is the number of source nodes over tolerance after phase
1 (Surprise 2). That is the first scan of phase 2, not its insertion count:
one insertion can bring its neighbours within tolerance (fewer), and a flip
can push a node that was within tolerance out of it (more). The table uses it
as the proxy, and acceptance measures the real count.

| tolerance | basin: nodes over tol. | as share of phase 1's vertices | piece: share of source nodes over |
|---:|---:|---:|---:|
| 1 m | 52 M | 44 % | 22.30 % |
| 2 m | 17 M | 30 % | |
| 5 m | 2.6 M | 15 % | |
| 10 m | 0.53 M | 8 % | 0.31 % |
| 20 m | 0.12 M | 5 % | |

**Mesh size.** The output grows by about the insertion count: up to roughly
+44 % vertices (and triangles) at 1 m and +8 % at 10 m over the basin. That is
the price of the guarantee being about the source DEM, and Ola accepted it
(Q13 (a), 2026-10-01).

**Time, estimated** from the piece's measured rates (not measured for phase 2
itself):

- *Insertions.* Phase 1 at 1 m spends 4.9 s in serial split and flip for
  about 5.2 M insertions, about 1 µs each. Phase 2's ~2.9 M first-pass
  violators on the piece (22.30 % of 13.0 M) would be ~3 s, plus something for
  the incircle: phase 2's vertices are off-node, so 21b's integer incircle
  does not answer for them; 21b's integer incircle cut refine by about 15 %
  at 8 threads, so phase 2 loses that gain. Say 3-4 s.
- *Scans.* Phase 1's scans total about 1.4 s at 1 m (6.33 s refine less the
  4.9 s split). Phase 2 scans every triangle once over 13.0 M points, then
  only touched ones: about 1-2 s.
- *Fixed cost, every tolerance:* projecting the source box's 31.5 M nodes
  (the prototype resampled 30.9 M nodes in 0.99 s on 8 threads, inverse
  transform and bilinear together) and sorting about 13-20 M points into the
  store (a serial sort, ~1-2 s).
- So on the piece: **about +6-8 s at 1 m** on today's 9.67 s process, and
  about +2-3 s at 10 m on 1.44 s, where the fixed cost dominates.

**At basin scale** (here only to place it; the run is now the basin-scale
design's, Q15): the fixed cost grows
with the source, 734 M ANADEM nodes: roughly 0.5-1 min of projection and a
serial sort of over a minute; the insertions grow with the tolerance (2.6 M
at 5 m: seconds). If the serial sort shows in the basin-scale run's profile, it is
parallelised there.

**Memory.** *Superseded at basin scale (Q15, 2026-10-01).* Everything in this
paragraph after the piece's figure assumes one process holding the whole
basin. Ola ruled that memory must not decide the tolerance, and directed
decomposition along constraint edges, windowed reads and no dense canvas,
designed next in its own increment file. The arithmetic below is kept as the
record of why one process does not reach the basin, not as a limit on it.
The store is 16 B per source node: 208 MB on the piece, 11.7 GB
for the basin's 734 M ANADEM nodes. D1 drops the target canvas (8.6 GB at
basin scale) before the store is built, so phase 2's peak is store plus mesh
(about two copies of it: phase 1's output and phase 2's lattice mesh, which
the basin-scale design measures), not store plus canvas plus mesh. Against Surprise 1's floors (13 GiB at 10 m,
19 GiB at 5 m, canvas included), phase 2 moves the basin's peak to roughly
16 GiB at 10 m and 22 GiB at 5 m: floors plus arithmetic, not measurements.
Those floors count the canvas at 4 B per node; the piece's own intercept was
about 17 B per grid node (Surprise 1), and at that rate the canvas term alone is
about 36 GiB, so the basin fits at no tolerance. So **one process on 32 GB
cannot reach the basin at 1-2 m, with or without the final check, and reaches
it at 5-10 m only if the piece's per-node intercept does not grow with the
canvas**; below that the route is
domain decomposition (ROADMAP item 2.4), which Ola's Q15 direction now
makes the route at every tolerance.

## The edge strip (Surprise 3), placed

**The finding, from the code.** A triangle's error is 0 when its closed node
set is empty (`include/terrain/refinement/scan.hpp@434c374:56`), and a
triangle is split only when its scan names a node
(`include/terrain/refinement/refine.hpp@434c374:187-189`). By Pick's theorem a
triangle whose three corners are nodes and whose closed node set holds no
other node has area half a cell, so the defect needs at least one off-node
corner: a domain or feature vertex (16, 16b) or a constraint foot (20b). It
is therefore not confined to the domain's edge: every constraint, lake shores
and CORINE borders included, can carry such slivers. That follows from the
code; it was measured only along the domain's edge.

**What 15c does about it.** On the reprojected path, phase 2 scans source
nodes in every triangle, off-node corners or not, so source nodes in a strip
are checked and inserted like any other. Surprise 3's cases (541.6 m on a
four-corner box, 17.9-38.5 m on the BHO outline) are exactly such nodes. What
remains is narrower: a sliver thinner than the source spacing can still hold
no source node, and between nodes nothing is checked. Smaller, not gone.

**What 15c does not do: the projected path** (Norway, and any DEM meshed
directly), which has no phase 2. There the `--tolerance` claim of every
domain mesh is false along constraints by up to tens of metres today.

**Placement: a separate increment, right after 15c and before 15d** (ruled by
Ola, Q14 (a), 2026-10-01). Its design in one paragraph, as ruled: check
points at every crossing of a constraint edge with a grid line, z linear
between the two nodes of that cell side (which is the bilinear surface there,
16 R0), **plus the midpoint between each two neighbouring crossings along the
edge**, z bilinear from the one cell that holds it (the piece of a straight
edge between two consecutive grid-line crossings lies in one closed cell),
run through 15c's `refine_points` after `refine`, for any DEM. The guarantee
it adds reads: "at every DEM node, wherever a constraint crosses a grid line,
and at the midpoint between neighbouring crossings". Proposed for that
increment's design, and Ola's to overrule: an edge's two end vertices count
as neighbours too, so the pieces at the ends of the edge get a midpoint as
well. Why not inside 15c:

1. it changes every Norwegian mesh with a domain or features, which ends J1's
   bit-identity on purpose and needs `@perf`'s acceptance on the 1 m
   benchmark;
2. the guarantee it adds is new in kind, at points between nodes (worded
   above, as Ola ruled);
3. its points lie on constraint segments, where a computed crossing can round
   a hair outside the domain and so fall in no triangle: they must be filed
   against their edge, not found by membership. That is 20b's foot geometry,
   and it needs its own red suite;
4. 15c's two PRs are already about 370 and 405 lines.

Estimated 160-210 lines (the midpoints add about 10), one PR, reusing
`CheckPoints` and `refine_points`.
Why not wait longer: until it ships, the tolerance printed on every Norwegian
domain mesh is wrong along its edges. A cheaper stopgap exists (split any
constrained edge of a node-free triangle longer than one cell at its midpoint,
about 40 lines): it narrows the strip and guarantees nothing.

## Invariants

- **J1. Norway untouched.** `refine.hpp` is not edited. A projected DEM
  without `--out-crs`, or with `--out-crs` equal to its CRS, gives a `.vtk`
  bit-identical to master's.
- **J2. The guarantee, against the source.** For every check point p (a valid
  source node inside the domain, at its stored position) that does not
  coincide with a start vertex, and every closed output triangle containing
  p, `|plane(p) − z_p| ≤ tolerance`, the plane from the output vertices' z.
  Phase 1's guarantee at the resampled grid's nodes does **not** survive phase
  2, and is not claimed: the grid is not the truth.
- **J3. Determinism.** Phase 2's output is bit-identical for any thread count
  and any order of `CheckPoints.add` calls; resampling is bit-identical for
  any block size and thread count.
- **J4. No degrees in `_core`.** `to_core` refuses a geographic tile;
  `CheckPoints` and `refine_points` take metres in the target CRS.
- **J5. One reprojection site**, `crs.reprojector`, `always_xy` (15 I8).
- **J6. Global node identity.** Target node `(R, K)` is at `(K·h, −R·h)` bit
  for bit, in every run.
- **J7. One check point per source node**, taken from the assembled mosaic,
  so overlaps resolve once (15 R5); duplicates are counted.
- **J8. Phase 2 terminates**, and every vertex it inserts is a check point
  carrying the point's own z.
- **J9. Constraints hold.** Phase 2 splits constraint edges (both halves keep
  the bit and the mask) and never flips one.
- **J10. Header before pixels.** CRS kinds, the antimeridian and poles, the
  memory cap and coverage are refused before any tile is loaded (15 I6).

## Degeneracy policy

- **A check point on an edge** (one zero orientation): `split_edge`; on a
  constrained edge it becomes a Steiner point of the constraint; on the
  domain's boundary (no neighbour), a 1 → 2 split.
- **On a vertex:** never inserted; counted as coincident if the vertex is a
  start vertex (D5 step 5); if it is a vertex phase 2 inserted, it is that
  vertex's own point.
- **Two points at one stored position:** one kept, `duplicates` counted.
- **Slivers:** membership by exact predicates; where a sliver's doubled area
  rounds to 0 or below, the error is bounded by the largest corner
  difference, as the grid scan does (`scan.hpp`).
- **Void triangles** (a NaN corner): 14's carve rule over check points.
- **NoData:** a NoData or NaN source node is not a check point; a target node
  whose stencil touches NoData is NoData.
- **Unequal source spacings** (`du ≠ dv`, GLO-30 above 50°): allowed; the
  default `h` uses the north-south spacing.
- **No image in the target CRS** (pyproj's `inf`): a domain vertex is
  refused; a source node is dropped (it is outside the domain).
- **A point on the target grid's far edge:** filed in the last cell.
- **An empty store** (a domain smaller than a source cell): one scan, nothing
  inserted, `0 source nodes` recorded.
- **The domain in a third CRS** (BHO in EPSG:4674 over an EPSG:4326 DEM): it
  goes straight to the target CRS, once (15b).
- **Datum shifts** between the source and target CRSs: whatever PROJ picks,
  recorded in `source_transform`.

## Not in scope

- **The basin-scale work** (Q15's direction, its own design next): window
  decoding with no dense canvas, a fetch step caching COG tiles,
  decomposition along constraint edges, the whole basin. 15d's design in
  `15-dem-mosaic.md` is to be re-read against it.
- **The edge strip on the projected path:** its own increment (above, Q14,
  ruled).
- **Constraint feet and start quality in phase 2** (D5).
- **One loop shared by `refine` and `refine_points`** (D5): a later refactor.
- **`rasputin catchment` on a geographic DEM:** refused by the `to_core` gate.
- **Tiles in more than one CRS, and N2's half-cell tiles with their
  neighbours:** the source side is still one 15a mosaic, which refuses both.
- **A CLI option for the grid spacing** (ruled out, Q16 (a)), and vertical
  datums.

**Decided here, which Ola may overrule:** `--bbox` on the reprojected path is
read in the target CRS; z in the store is `float`; check points on a start
vertex are counted, not adopted; phase 2 without feet; a second loop rather
than a shared one; the memory cap sums both canvases and the store. For D8:
only conformal families (a reading of Ola's "distorts areas and angles the
least": zero angle error, area reported, not minimised); the four families'
thresholds (70° and 80° for polar, tall-or-square for TM, 15° for Mercator);
one balanced scale factor per family, with the LCC's parallels and the
stereographic's pole not optimised; angles rounded to 0.1° and `s` to 6
decimals; the DEM's own ellipsoid; the 21 × 21 evaluation grid; and the DEM
footprint as the box when `--bbox` is given without `--out-crs`.

## PR split and LOC

Counted in `CLAUDE.md` §2's unit. Estimates; increment 10 came in 39 % over
and 15a's `mosaic.py` 60 % over, so the worst case applies 39 %, with
60 % beside it, and the largest single module is flagged. Both PRs stay under
700 at either.

| | what | est. | worst |
|---|---|---:|---:|
| **15c-1** | **The engine: check points and the final check, in C++. No CRS anywhere** | | |
| | `include/terrain/refinement/check_points.hpp`: store, `add`, `freeze`, query, `outside` | 77 | |
| | `include/terrain/refinement/refine_points.hpp`: `scan_points` | 70 | |
| | `refine_points.hpp`: the loop, `PointRefineOutcome`, coincident points, output | 110 | |
| | `include/terrain/mesh/lattice_mesh.hpp`: `split_inside(MeshVertex)` | 8 | |
| | `bindings/core.cpp`: `CheckPoints`, `refine_points`, GIL release | 75 | |
| | `_core.pyi`: stubs | 30 | |
| | **15c-1 total** | **370** | **514 (592 at +60 %)** |
| **15c-2** | **The geographic path, in Python, ending with the final check** | | |
| | `io/geotiff.py`: geographic 2D, degree units | 30 | |
| | `io/models.py`: `crs`, `geographic`, `epsg: int \| None` | 15 | |
| | `raster.py`: the gate; `mosaic.py` and `catchment.py` (`catchment.py:134-135`): `EPSG:{...}` through `meta.crs` | 13 | |
| | `target_grid.py`: `TargetGrid`, spacing, extent, `source_region` | 80 | |
| | `target_grid.py`: `SourceWindows`, `TileWindows`, `resample` | 55 | |
| | `target_grid.py`: `check_point_blocks` | 35 | |
| | `dem_input.py`: `target_crs`, the reprojected branch (geographic, and a projected DEM in another CRS, Q12), the memory sum | 45 | |
| | `crs.py`: `suggest_crs`, `CrsSuggestion`: four families, the balanced scale, rounding, `get_factors` over the box (D8) | 29 | |
| | `final_check.py`: the two phases and the store | 40 | |
| | `cli.py`: `--out-crs`, Q11's refusal with D8's message, `_dem_mesh` through `final_check`, fields, stats rows | 63 | |
| | **15c-2 total** | **405** | **563 (648 at +60 %)** |

`target_grid.py` at 170 is the module to watch (15a's `mosaic.py` overran its
estimate by 60 %). If it passes 250, `resample` and `check_point_blocks` move
to `resample.py`.

**Why this split, and this order.** 15c-1 is the risky, novel part, and it
needs no CRS: `@tester` can drive it with scattered points and a hand-built
start mesh, and `@perf` can confirm Norway is untouched. 15c-2 is plumbing on
a prototype that already ran (`prep_dem.py`), and it wires the final check in
the same PR that first lets a geographic DEM through, so **no release ever
writes a geographic mesh without the final check**. Together they are about
775 lines, over the ceiling, so they cannot be one PR. The 2026-10-01
rulings moved 15c-2 from 383 to 405: Q11's suggestion is about 17 lines more
than the single fitted LCC first estimated (the refusal itself moves from
`cli.py` into `crs.py`), and Q12 (a) adds 5. Q13, Q14, Q16 and Q17 add no
production lines to 15c.

**Acceptance class.** 15c-1 adds files under `include/terrain/refinement/`
and touches `include/terrain/mesh/`, so the refine/mesh acceptance rule
applies (`docs/increments/README.md`, "Acceptance"). 15c-2 changes no C++,
but it restructures `_dem_mesh`, which drives refine on Norway's path, so the
rule applies to it too (below).

**Documentation in the same PRs** (not counted): `project_structure.md`, the
`raster` boundary rule as restated in D6 (15c-2); `ROADMAP.md`'s 15 row at
each merge; `NOTICE.md` and a fixture `NOTICE` for the GLO-30 and ANADEM
extracts and the BHO outline (Q17, below). `15-dem-mosaic.md` already points
here (this branch).

## Tests for `@tester`

**Invariant-critical suite, for mutation testing:** 15c-1's J2 oracle (RP3)
and the membership and tie-break tests of `scan_points` (RP2, RP4, RP5).
Everything else is ordinary.

**15c-1, C++ (Catch2) and through the binding:**

- **CP1, the store.** Points are filed by cell; the iteration order after
  `freeze` is the same for every permutation of the `add` calls; two points at
  one position count one duplicate; `add` after `freeze` is refused; a point
  on the far edge is in the last cell; a point outside the node rectangle,
  or with a NaN coordinate or z, is dropped and counted in `outside`.
- **RP1, a plane.** Check points whose z lie on the start mesh's planes:
  nothing inserted, `max_error` 0 to 1e-9.
- **RP2, one bump.** Offsets are exact binary fractions (0.25, 0.5), so the
  stored position equals the given one. One point above tolerance inside a
  triangle is inserted at that position with its own z; one below is not; one exactly at the
  tolerance is not (`>`, as `needs_split`).
- **RP3, J2 by an independent oracle.** Scattered points over a rough
  function, several tolerances (0 included). After `refine_points`, a
  brute-force check in NumPy: for every point, every output triangle whose
  closed area contains it (barycentric, with a 1e-12 slack, so a point on an
  edge is tested in both triangles), error at most tolerance + 1e-9. The
  oracle does its own location, never using the store's order or the scan's
  records (the computational-geometry skill's "producer's relation" rule:
  borrow the predicate, not the records). It must be shown to fail: plant a
  mesh whose z is shifted by twice the tolerance.
- **RP4, determinism.** Threads 1, 2 and 8, and two `add` orders: identical
  outputs, bit for bit.
- **RP5, constraints and edges.** A point exactly on a constrained edge
  splits it; both halves keep the bit and the mask; no constraint is ever
  flipped; a point on the domain's boundary splits 1 → 2.
- **RP6, void.** A start vertex with NaN z: the carve rule, `uncovered` 0 at
  the end.
- **RP7, coincident.** A point exactly at a start vertex with a different z:
  not inserted, `coincident` 1, `coincident_max_error` the difference.
- **RP8, outside.** Points outside every triangle are never inserted and
  never raise `max_error`.
- **RP9, the binding.** Shape and dtype refusals for `add`; the GIL is
  released (a Python thread advances while `refine_points` runs).

**15c-2, Python:**

- **G1, the reader.** A geographic micro-TIFF with ANADEM's header numbers
  (B2) reads with `geographic` true and `crs` `EPSG:4326`; one with the
  OpenTopography COG's key set (2048 = 4674 with 2049, 2054 = 9102, 2057,
  2059; "Measured" above) reads with `crs` `EPSG:4674`; geocentric, 3D,
  user-defined and non-degree geographic keys are refused, each by name.
- **G2, the gate.** `to_core` refuses a geographic tile; `rasputin catchment`
  on one is a usage error.
- **G3, the target grid.** Nodes at multiples of `h`, bit for bit; the extent
  covers the grown domain; default spacing 30 for ANADEM's spacing and 31 for
  GLO-30's, at 19°S.
- **G4, resampling.** A source whose values are affine in its own index space
  is reproduced at every target node to 1e-9 (bilinear is exact there); the
  NoData stencil rule; identical arrays for 1 and 8 threads and two block
  sizes.
- **G5, check points.** Every valid source node inside the domain appears
  exactly once, against a brute-force projection and `shapely.contains_xy`;
  block pruning never drops an inside node (a domain with a thin arm through
  a block corner).
- **G6, end to end.** `rasputin mesh` on a synthetic geographic tile with a
  rough surface, a domain in EPSG:4674 and `--out-crs EPSG:31983`: the `.vtk`
  is in the target CRS, records D7's fields, and an independent final check
  (the shape of `basin-piece/run_sweep.py`'s `final_check`) finds 0 source
  nodes over tolerance, interior and strip. The same run without phase 2
  (through the Python API) finds some: the test can fail.
- **G7, J1.** On the committed projected tile with a domain, `--out-crs` set
  to the DEM's own CRS writes the same bytes as no `--out-crs`.
- **G8, refusals, header before pixels.** A geographic DEM without
  `--out-crs` (Q11 (a)), with D8's suggestion in the message; a source
  region across ±180°; the memory sum over the cap; a domain vertex without
  an image. None loads a tile (a repository double that fails on `load`).
- **G9, a projected DEM in another CRS (Q12).** A metre micro-TIFF in
  EPSG:25833 with `--out-crs EPSG:25832` takes the reprojected path: the
  `.vtk` is in EPSG:25832, records D7's fields, default spacing is the
  source's, and G6's independent final check finds 0 source nodes over
  tolerance.
- **G10, the suggestion (D8).** One box per family, at the thresholds and
  just past them (φc 70° against 69.9°, a box 0.01° taller than wide against one 0.01° wider,
  φc 15° against 15.1°), gives the family and the rounded parameters. The
  scale is balanced: by the oracle below, the largest and smallest point
  scales over the box sit at `1 ± max_scale_error` to within the rounding of
  `s`. `CrsSuggestion.max_scale_error` and `max_areal_error` (the fields, not
  the rounded printed text) are checked by an independent oracle, not by
  calling `get_factors` again: finite differences of the suggestion's forward
  transform over short geodesic steps (pyproj `Geod`) at the 21 × 21 points,
  agreeing to 1e-6 (`@tester` may loosen it if the finite differences need it,
  never past 1e-4). It must be shown able to fail: the oracle run against the
  same suggestion with its scale factor replaced by 1 differs by about half
  the unit-scale figure (2.5e-3 on the basin box).

**Fixtures.** Synthetic micro-TIFFs for G1-G10. Real ones, per Q17 (a), each
committed only after its licence is read and quoted in the fixture `NOTICE`
(the licence check is part of the commit, not of this design):

- about 512 × 512 source nodes of the Velhas piece from **GLO-30**, with the
  Copernicus credit;
- the same window from **ANADEM**, cut by range reads from OpenTopography's
  COG (the URL and DOI 10.5069/G9736P4G in the `NOTICE`, with the window and
  the date it was cut), crediting ANADEM; it is EPSG:4674, so G6's realism
  case also exercises the COG's key set;
- the **BHO-derived Velhas outline**, the union of 1,163 BHO 2017 5k
  catchments of ottobasin 76949 (as `basin-piece/fetch_bho.py` built it),
  with ANA's credit and the CC BY 4.0 notice, and the `NOTICE` saying that
  ANA's metadata record was not read (its host answered 403) and is to be
  checked when reachable, and that BHO 2017 5k is marked superseded by BHAE
  in ANA's catalogue. The outline spans the whole piece and the extracts do
  not, so G6's realism case uses the outline intersected with the extract's
  footprint, computed in the test.

## Acceptance

- **15c-1:** CI green. `@perf`, per `docs/increments/README.md`
  "Acceptance": `tools/bench.py`'s 1 m benchmark and thread-scaling sweep,
  with `pmset -g batt` recorded for each run, against the previous
  increment's run in the same power state (on `NO BASELINE`, master's merge
  commit with `--tree`, back to back). The mesh hash is unchanged and refine
  time is within noise at every thread count. Evidence under
  `docs/benchmarks/<date>/`.
- **15c-2:** `@perf` on the Velhas piece (BHO 76949) through `rasputin mesh
  --out-crs EPSG:31983`, on **ANADEM**, a local window cut from
  OpenTopography's COG by range reads (a one-off cut until the basin-scale
  design's fetch step exists), and once on GLO-30 for comparison with
  Surprise 2, at
  tolerances 1, 2, 5, 10, 20 and 50 m: triangles, phase 1 and phase 2 times,
  the fixed cost (projection, store), peak RSS, and phase 2's insertions
  against Surprise 2's first-pass counts (approximate: Surprise 2 was
  measured on a 30 m grid, and the CLI default gives 31 m on GLO-30). The independent final check of the
  basin-piece run reports **0 source nodes over tolerance, interior and
  strip**, and its control still fails when the mesh is shifted 15 m. Worst
  angle and maximum degree against phase 1's. Phase 2's thread scaling at
  1 m. For information only: phase 2 alone from the domain's start mesh
  (through the Python API), to show what phase 1 buys. And, for J1 on the
  real path: `tools/bench.py`'s 1 m benchmark and thread-scaling sweep
  against 15c-1's merge commit, same power state, with the mesh hash
  unchanged and refine and process time within noise.
- For both: every gate in `CLAUDE.md` §4 green, and CI green.

## Questions for Ola

Numbered on from `15-dem-mosaic.md`'s Q1-Q10. Kept as asked; Ola's rulings
of 2026-10-01 are recorded in full under "What is ruled" and marked on each
question here. Q15 is still open.

**Q11 (was Q7). Who chooses the target CRS for a geographic DEM?** It is now
the CRS the mesh is computed in, not only written in. *Ruled (a), with the
suggestion chosen by the box's location and shape (D8).*
- **(a) Required `--out-crs`; the refusal prints an LCC fitted to the domain
  (Snyder's one-sixth rule) to copy. Recommended.** About 15 lines. The CRS is
  a contract with whatever reads the mesh, and two catchments, or two
  subdomains later, land in the same CRS only if someone chose it.
- (b) Fit an LCC to each domain automatically, recorded as WKT2. About 45
  lines. Two catchments meshed separately end up in different CRSs.
- (c) The UTM zone of the domain's centroid, automatically. About 15 lines,
  an EPSG code. The basin spans zones 23 and 24, so one zone is stretched up
  to 9° past its central meridian: still conformal, scale off by about 1 %.

**Q12. `--out-crs` on a projected DEM in another CRS** (say DTM10 in UTM 33,
mesh wanted in UTM 32). *Ruled (a).*
- **(a) Take the same path: resample onto the target grid, then the final
  check. Recommended.** About 5 lines; the path does not care whether the
  source is geographic.
- (b) Refuse, as today.

**Q13. The final check's price, now measured (Surprise 2).** You ruled it in
before its cost was known. *Ruled (a).*
- **(a) As ruled: insert every source node still over tolerance.
  Recommended.** The guarantee is then about the DEM you gave. Over the basin
  that is up to about +44 % vertices at 1 m, +15 % at 5 m, +8 % at 10 m; on the
  piece about +6-8 s at 1 m.
- (b) Check the source at a looser, second tolerance. Fewer insertions, a
  weaker and two-number guarantee, about 10 more lines.
- (c) Measure and record only, insert nothing. The guarantee is then against
  the resampled grid, which is off the source by p99 3.1 m and max 27.8 m on
  the piece.

**Q14. The edge strip (Surprise 3) where there is no final check** (Norway,
any DEM meshed directly). Today `--tolerance` does not hold in slivers along
constraints; up to tens of metres on a dense outline. *Ruled (a), plus the
midpoint between neighbouring crossings ("The edge strip").*
- **(a) Its own increment right after 15c, before 15d: check points at every
  crossing of a constraint with a grid line, z bilinear there, through 15c's
  `refine_points`. Recommended.** 150-200 lines, one PR; changes every
  Norwegian domain mesh, so it gets `@perf`'s acceptance. The guarantee would
  read "at every DEM node, and wherever a constraint crosses a grid line";
  please confirm that wording or give yours.
- (b) Inside 15c. 15c-2 would approach the ceiling, and Norway's acceptance
  would ride on the basin work.
- (c) Wait. The printed tolerance stays wrong along constraints.
- (d) A stopgap: split long constrained edges of node-free triangles at their
  midpoints. About 40 lines, narrows the strip, guarantees nothing.

**Q15. The basin tolerance** (open since 2026-09-30). *Still open, and
reframed by Ola on 2026-10-01: memory must not decide the tolerance; the
direction (decomposition along constraint edges, windowed reads with no
dense canvas, a cached fetch step) is recorded under "What is ruled" and gets
its own basin-scale design. The options below rest on the one-process memory
reasoning that direction supersedes, and are kept as asked.* With Surprise 1 and
this design, 1-2 m does not fit one process on 32 GB with or without the
final check; 5-10 m may, if the piece's per-node intercept does not grow with the canvas
(at its 17 B per node the basin fits at no tolerance), which 15d measures.
- **(a) Tell us what the hydrology needs; meanwhile 15d's acceptance targets
  10 m, then 5 m. Recommended.**
- (b) 1-2 m is needed: then domain decomposition (ROADMAP item 2.4) comes
  before the whole-basin run, and 15d stops at the piece.

**Q16. The target grid's spacing.** *Ruled (a).*
- **(a) The source's north-south spacing at the domain's centroid, rounded to
  whole metres (30 m for ANADEM, 31 m for GLO-30); no CLI option.
  Recommended.** No extra lines. B6 found a finer grid does not remove the
  resampling error, and the final check removes it anyway.
- (b) Also a `--grid-spacing` option. About 10 lines.

**Q17. Fixtures while ANADEM's host refuses us** (HTTP 403 since 30
September). *Ruled (a) for both GLO-30 and ANADEM, ANADEM from
OpenTopography's COG, and the BHO outline as a fixture ("Fixtures").*
- **(a) Commit a small GLO-30 extract of the Velhas piece now, with the
  Copernicus credit, and an ANADEM one when the host answers. Recommended.**
  GLO-30's licence is recalled to allow redistribution with its notice; it is
  read again before anything is committed.
- (b) Synthetic fixtures only until ANADEM is reachable.

Still open from `15-dem-mosaic.md` and not this design's: whether commercial
use matters for the DEM choice.

## Review

### Round 1, `434c374..9af0b47`: CHANGES REQUESTED (`@reviewer`)

Production LOC 0 (docs only: this file +861, `15-dem-mosaic.md` +8). Prior
art, the I/O boundary, the PR split and Q11-Q17 pass; every number and code
claim checked against the tree. Four blocking prose edits: B1 `RefineOutcome`
gained fields although `refine.hpp` was "not edited" (fixed by
`PointRefineOutcome`); B2 the memory caveat understated (at 17 B per node the
basin fits at no tolerance); B3 the acceptance did not name `tools/bench.py`'s
benchmark and thread-scaling sweep for either PR; B4 a line count this branch
made false.

### Round 2, `9af0b47..529fe08`: APPROVED (`@reviewer`)

Production LOC 0 (+55 / −32 in this file). B1-B4 fixed; suggestions taken.
Estimates: 15c-1 370 (514 at +39 %, 592 at +60 %), 15c-2 383 (532, 613),
each under 700, together about 750, so two PRs. The deviation (`refine_points`
reads the store's `geometry()` instead of taking `grid`) accepted: a mismatch
can no longer be expressed. Gates green. CI not yet run. Non-blocking, for the
red step: name `geometry()` beside `for_each_in` in D4 and the comment at the
`refine_points` signature; "`grid` is only the frame" → "the store's
`geometry()` is only the frame"; add a ±inf case to CP1.

### Round 3, `8b5315c..0a7ca47`: CHANGES REQUESTED (`@reviewer`)

Production LOC 0 (+308 / −64 here, +3 / −1 in `15-dem-mosaic.md`). Ola's
2026-10-01 rulings on Q11-Q17 recorded; the distortion table, the ANADEM COG
header and LOC (370, 403) verified. Two blocking edits: B1 Ola's Q11 words
quoted for keeping TM's scale at 1, though Ola made the excuse conditional
("if that is time consuming") and balancing is one line; B2 three sentences
still gave basin-scale work to 15d.

### Round 4, `0a7ca47..2741dba`: APPROVED (`@reviewer`)

Production LOC 0 (+79 / −53). B1 fixed by balancing the scale factor in
every family (all six rows reproduced, e.g. the basin 0.25 % scale, 0.49 %
area, against 1.11 % for UTM); B2 fixed; suggestions taken. 15c-1 370,
15c-2 405 (563 at +39 %, 648 at +60 %). Gates green. CI runs once #128 targets
master.
