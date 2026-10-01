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
  run on the whole basin. 15c is accepted on the piece (210 MiB-scale), and
  refuses a request whose working set passes the memory cap (D6).
- **Out (placed below):** the edge strip for projected DEMs, and constraint
  feet in the final check.

## The design

### D1. Data flow

### D2. The target grid

### D3. Resampling

### D4. Check points: the source nodes, in the target CRS

### D5. The final check: `refine_points`

### D6. Types and the I/O boundary

### D7. What the file records

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
