# Increment 16b — terrain polygons: interior polygons and polylines as constraints (CORINE over Norway)

Status: **designed, not built.** Questions Q1-Q6 at the end are open. Written
by `@architect` before `@tester`, per `docs/increments/README.md` step 1, on
branch `increment16b-terrain-polygons` off the unmerged 15b branch
(`increment15b-domain-crs`). Nothing in this design depends on 15c or 15d.

**Closes.** `ROADMAP.md`'s 16b row, and increment 16's R6 ("16b's CLI,
designed, not built"). After this increment,

```sh
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20220924/ \
    --domain catchment.geojson --tolerance 1 \
    --features ../rasputin_data/corine_sql/.../U2018_CLC2018_V2020_20u1.gpkg \
    --features-layer U2018_CLC2018_V2020_20u1 --features-map corine \
    --out catchment.vtk
```

meshes the catchment with every CORINE land-cover boundary inside it as a
constraint, each constraint edge carrying the vocabulary bits its classes map
to, read without GDAL. The same works with a GeoJSON `FeatureCollection`
(`--features lakes.geojson`).

**Not closed.** Which triangles lie inside which polygon (a region label per
triangle; proposed as 16c, Q2). Coarsening input geometry (Ola's "later
stage"). Hydro-flattening lakes. A class map read from a file (Q4). The
computation frame for a geographic DEM (15c). 20c's soft quality criterion,
which this increment gives its first real input with interior constraints.

## Ola's direction this design is built on

1. **Input is geometry** (16, direction 1-3): one confining polygon, any
   number of polygons inside it, any number of polylines. "They are not
   discretized, they are already discrete. We will come to input geometry
   coarsening at a later stage." Vertices are used as given.
2. **Input vertices stay where they are, with bilinear z** (16 R2). A crossing
   the noder makes is one more such vertex (16 R6).
3. **16 R6 is the input model**: an interior polygon's exterior and holes are
   `Breakline` chains, closed; a polyline is an open `Breakline`; each carries
   the feature's property bits (increment 7). No new role.
4. **Real data is the working example** (2026-09-27): CORINE Land Cover 2018
   over Norwegian DTM10, read GDAL-free from Ola's local GeoPackage.
5. **Inputs come in their own CRS** and are reprojected into the DEM's, which
   is the computation CRS for a projected DEM (15b, R9). "Feature geometry
   (16b) follows the same path" (ROADMAP, "inputs in their own CRS").
6. **Order of work** (2026-09-27): 15a/15b, then 16b, then auto-catchment.

## What was measured

In the scratchpad, never in the tree: `clc_probe.py` (read, reproject, clip,
statistics), `clc_noder.py` (two chain layouts through the real `_engine`),
`clc_mesh.py` (the whole `cli._dem_mesh` path with feature chains patched in).
The installed `_core` was built from this branch's C++, which has not changed
since `4be156c` (`git log -1 -- include src bindings`). On AC power
(`pmset -g batt`: "AC attached"), Apple clang Release, default threads.

### M1. The GeoPackage, read without GDAL

- One features table for Europe, `U2018_CLC2018_V2020_20u1`: columns
  `OBJECTID`, `Shape` (`MULTIPOLYGON`, EPSG:3035), `Code_18` (`TEXT(3)`, the
  CLC code), `Remark`, `Area_Ha`, `ID`. 2 375 406 rows in its R-tree
  (`rtree_U2018_CLC2018_V2020_20u1_Shape`, columns `id, minx, maxx, miny,
  maxy`). Five more tables for the French overseas departments, each in its
  own CRS. Read with `sqlite3.connect("file:...?mode=ro", uri=True)`.
- Every geometry blob seen (74 + 972 + 859 rows) has flags byte `0x01`:
  little-endian, **no envelope**, not empty, standard WKB after the 8-byte
  header. `shapely.wkb.loads` reads the rest.
- The table's insert/update triggers call `ST_IsEmpty`, a SpatiaLite function
  plain SQLite lacks. Read-only access never fires them.
- The R-tree query costs 11-95 ms per 50 km tile (74-972 candidate rows).
- The largest single feature over Norway is the sea polygon (`523`,
  `OBJECTID` 2372334): a 5.75 MB blob, about 360 000 vertices, 593 × 915 km.
  Any coastal request reads all of it.

### M2. What CORINE looks like over three DTM10 tiles

Features intersecting the tile's node rectangle, reprojected 3035 → 25833
vertex by vertex (15b's `reprojector`), clipped to the rectangle with shapely.

| tile | area | features | classes | raw vertices | clipped ring vertices | parse + reproject |
|---|---|---|---|---|---|---|
| `7908_3` (the benchmark tile, coast) | 71 % covered | 60 | 8 | 201 304 | 10 989 | 142 ms |
| `6602_2` (Oslo, fjord, forest) | 100 % | 875 | 22 | 1 186 709 | 119 715 | 1 742 ms |
| `6603_4` (forest, bogs, lakes) | 100 % | 727 | 17 | 979 519 | 114 900 | 1 460 ms |

`7908_3`'s 11 068 vertices in the ROADMAP row are the same count with each
ring's closing vertex included (79 rings). The uncovered 29 % is open sea
beyond CORINE's `523` polygon.

| tile | segments | min | p1 | median | < 1 m | < 10 m (a cell) | vertex pairs < 1 m apart | closest pair |
|---|---|---|---|---|---|---|---|---|
| `7908_3` | 10 989 | 0.268 m | 19.1 m | 54.5 m | 6 | 61 | 3 | 0.268 m |
| `6602_2` | 119 715 | 0.004 m | 0.71 m | 55.9 m | 5 061 | 8 772 | 2 639 | 0.004 m (52 pairs < 1 cm) |
| `6603_4` | 114 900 | 0.190 m | 0.61 m | 51.5 m | 9 626 | 16 201 | 4 991 | 0.190 m |

**CORINE is an exact partition, bit for bit, after reprojection.** Keyed on
the undirected pair of endpoint coordinates (Python floats, `==`):

| tile | segments shared by exactly 2 features | unshared | of which on the clip border | same-direction duplicates | shared with the same class | vertex within 1 mm of another segment's interior |
|---|---|---|---|---|---|---|
| `7908_3` | 5 114 | 761 | 40 | 0 | 0 | 0 |
| `6602_2` | 59 749 | 217 | 217 | 0 | 0 | 0 |
| `6603_4` | 57 356 | 188 | 188 | 0 | 0 | 0 |

- Every interior boundary arrives **twice, in opposite orientation**, with
  identical endpoint coordinates: pyproj maps equal inputs to equal outputs.
- There are **no T-junctions and no partial overlaps** within 1 mm.
- `7908_3`'s 721 unshared, off-border segments are the outer edge of the sea
  polygon and the coast where nothing lies beyond.
- No two neighbours share a class in these tiles, though CORINE's readme
  warns "some redundant lines between neighbouring polygons with the same
  code are still present" (`Documents/readme_U2018_CLC2018_V2020_20u1.txt`).
- 2 to 4 vertices per tile sit exactly on a DEM node; the rest are off-node.

### M3. The noder on CORINE: correct, and quadratic

A 10 km square and a 48 km square in `6603_4`, the square as `Outer`, the
features' linework inside it as `Breakline`s with one scratch bit, through
`cli._engine` (`build_pslg` → `node` → `triangulate`, 1 mm snap). Two layouts:

- **A, deduplicated in Python:** every segment once, rings clipped as lines
  (`MultiLineString ∩ square`), joined by `shapely.line_merge`.
- **B, rings as given:** every clipped polygon ring as a closed `Breakline`
  (16 R6's letter), so every shared edge is given twice, in opposite
  orientation, and the clip's border segments lie on the `Outer` ring.

| square | layout | chains | input vertices | chain positions | noded vertices | noded edges | start triangles | `node` |
|---|---|---|---|---|---|---|---|---|
| 10 km | A | 108 | 2 607 | — | 2 607 | 2 638 | 5 165 | 69 ms |
| 10 km | B | 63 | 2 607 | 5 295 | 2 607 | 2 638 | 5 165 | 223 ms |
| 48 km | A | 1 415 | 51 395 | — | 51 395 | 51 788 | 102 627 | 24.8 s |
| 48 km | B | 988 | 51 395 | 104 406 | 51 395 | 51 788 | 102 627 | 83.1 s |

- **The noder merges every doubled edge exactly and unions its bits.** A and
  B give the same noded edge set and the same mask on every edge, on both
  squares (compared as sets of coordinate pairs, 0 differences). This is 5b's
  amended guarantee 14(a) (duplicate edges permitted, equal node-id pairs) and
  its node-id dedup (`05b-noder-driver.md`, "The merge happens at the node-id
  edge-key dedup"), holding on real data.
- **The noder is quadratic in edges.** Layout A: 1 764 segments 0.03 s,
  8 568 0.69 s, 19 140 3.39 s, 51 791 24.8 s: each 2.7× more edges costs
  about 7×. Doubling the edges (B) costs 3.4×. Densifying the four 48 km
  domain edges to 100 m changed nothing (28.0 s), so long segments are not
  the cause.
- **The likely cause is 5b's verification pass**, read, not isolated:
  `NodedPslgBuilder::check_guarantee_14` (`noding/noded_pslg_builder.hpp`)
  checks 14(a) over every pair of edges and 14(b) over every (edge, node) pair,
  "NO BROAD PHASE ... The candidates are small by construction -- a noded
  constraint set, not a mesh". CORINE makes that premise false. A `sample`
  profile puts almost all time in one function called from the noder's
  driver; the symbols are stripped, so it names no function. A scratch build
  with the check disabled, which would isolate it, was not run (the sandbox
  refused the build). R8 makes the isolation 16b-0's first step.

### M4. Meshing with CORINE: the whole pipeline

The 48 km square in `6603_4` (2 304 km², forest and bogs), layout A, through
`cli._dem_mesh` (quality start 25°, constraint feet on), against the same
square with no features.

| tol | features | start quality nodes | feet | triangles | min ∠ median | < 1° | worst | `node` | refine | whole |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 m | none | 0 | 0 | 6 444 390 | 36.9° | 0.00 % | 0.029° | — | — | 5.1 s |
| 1 m | CORINE | 154 640 | 30 045 | 6 612 554 (+2.6 %) | 36.9° | 0.03 % | 0.0007° | 24.7 s | — | 30.2 s |
| 10 m | none | 0 | 0 | 132 563 | 34.5° | 0.02 % | 0.24° | — | 0.5 s | 0.55 s |
| 10 m | CORINE, quality off | 0 | 354 | 211 487 (1.6×) | 28.5° | 3.47 % | 0.015° | 24.7 s | 0.22 s | 25.0 s |
| 10 m | CORINE, quality 25° | 154 640 | 341 | 480 127 (3.6×) | 36.9° | 0.25 % | 0.0013° | 24.7 s | 0.37 s | 25.2 s |

- Every run met its tolerance, and no valid node was left uncovered.
- At 1 m the DEM dominates: CORINE adds 2.6 % triangles. At 10 m the input
  geometry dominates: 51 395 forced vertices, median 50 m apart, against a
  mesh that would otherwise need 66 523.
- **The quality start triples the 10 m mesh.** Dense boundaries make many
  start triangles with circumradius above the 14 m floor, and increment 20's
  pass inserts 154 640 DEM nodes to fix them. That is 20c's question (soft
  criterion, C1-C3), now with real data to decide it on.
- **Constraint feet (20b) work on interior constraints**: 30 045 at 1 m, 0
  refused. 20b R8 had this as "not measured".
- The worst angles (0.0007°, 0.0013°) come from input segments shorter than a
  centimetre and vertices closer than a metre (M2), which Ola's model keeps.
- Peak RSS 2.4 GB at 1 m (6.6 M triangles), 0.68 GB at 10 m.

## Prior art: legacy and literature

### Literature

- **Constrained Delaunay triangulation with breaklines**: Chew, "Constrained
  Delaunay triangulations", *Algorithmica* 4:97-108, 1989. Unchanged; 16b
  only supplies more constraints. Recalled, not reread.
- **Snap rounding**, which the noder is (increment 5): Hobby, "Practical
  segment intersection with finite precision output", *Comput. Geom.*
  13:199-214, 1999; Goodrich, Guibas, Hershberger and Tanenbaum, "Snap
  rounding line segments efficiently in two and three dimensions", SoCG 1997;
  cited in `05-noder.md`. 16b changes only how the result is verified (R8).
- **Finding intersecting pairs without checking every pair.** Bentley and
  Ottmann, "Algorithms for reporting and counting geometric intersections",
  *IEEE Trans. Comput.* C-28(9):643-647, 1979 (plane sweep); and the simpler
  one-axis "sort and sweep" over bounding-box intervals used as a broad phase
  in collision detection (Cohen, Lin, Manocha and Ponamgi, "I-COLLIDE", *Symp.
  Interactive 3D Graphics* 1995; Baraff's thesis, Cornell, 1992). Both
  recalled. R8 uses sort and sweep for the verifier. **Difference:** it only
  proposes candidate pairs for the exact `classify<K>` and
  `segment_meets_cell<K>` the verifier already runs, so it decides nothing and
  needs no exactness of its own.
- **Polygon clipping** is GEOS's OverlayNG (Martin Davis, 2020, the engine
  behind `shapely.intersection` since GEOS 3.9; documentation recalled). 16b
  clips linework, not areas (R6), and only in Python.
- **Region attributes on a CDT**: Shewchuk, "Triangle: Engineering a 2D
  quality mesh generator and Delaunay triangulator", *Applied Computational
  Geometry*, LNCS 1148, 1996. Triangle's `-A` spreads a regional attribute
  from a seed point across unconstrained edges until it meets constraints.
  That is the method for 16c's region labels (Q2), recorded here so 16c starts
  from it. Recalled.
- **Small input features limit mesh quality**: Ruppert, "A Delaunay refinement
  algorithm for quality 2-dimensional mesh generation", *J. Algorithms*
  18:548-585, 1995; Shewchuk, "Delaunay refinement algorithms for triangular
  mesh generation", *Comput. Geom.* 22:21-74, 2002. Refinement's size is
  governed by the local feature size, which a 4 mm input segment sets to 4 mm.
  Neither increment 20's quality start nor 20b's feet claims Ruppert's bound;
  M4's worst angles are what that costs with CORINE. Recalled.
- **Land cover and hydrology in TINs**: Jones, Wright and Maidment,
  "Watershed delineation with triangle-based terrain models", *J. Hydraulic
  Eng.* 116(10):1232-1251, 1990; Vivoni, Ivanov, Bras and Entekhabi,
  "Generation of triangulated irregular networks based on hydrological
  similarity", *J. Hydrologic Eng.* 9(4):288-302, 2004 (tRIBS). Both use
  hydrographic and land-cover lines as TIN constraints. Recalled, not reread;
  cited for the practice, not for a method carried here.
- **Coverage simplification**, for the later coarsening increment, not here:
  Visvalingam and Whyatt, "Line generalisation by repeated elimination of
  points", *Cartographic J.* 30(1):46-51, 1993; GEOS's `CoverageSimplifier`
  (2024) applies it to a polygon partition without opening gaps. CORINE's
  exact partition (M2) is what such a method needs.
- **GeoPackage**: OGC 12-128r18, *GeoPackage Encoding Standard* 1.3.1, 2021:
  the geometry blob header (magic `GP`, version, flags, `srs_id`, envelope)
  in its section 2.1.3 and the R-tree extension in Annex F.3. Recalled; R3's
  header layout was checked against the real file (M1), not against the text.

**Novelty: none claimed.** Nothing here is new. If a later increment claims
anything about exact partitions as TIN constraints (for example a region
labelling with a guarantee), it searches first: Google Scholar for "land cover
polygons constrained Delaunay terrain", "coverage constraints TIN
generation", "regional attributes constrained triangulation". No web search
tool was available to this round.

### Legacy

```sh
$ grep -rliE 'corine|land_cover|LandCover|clc|geopackage|gpkg|sqlite' legacy/
legacy/rasputin/land_cover_repository.py
legacy/rasputin/globcov_repository.py
legacy/rasputin/wfs_repository.py
legacy/rasputin/tin_repository.py
legacy/rasputin/application.py
legacy/rasputin/web_visualize.py
legacy/rasputin/gml_repository.py
legacy/tests/test_gml_repository.py
legacy/tests/test_land_cover_repository.py
```

What matters, read directly:

- `legacy/rasputin/gml_repository.py:13-64`: CORINE's 44 classes as an
  `Enum` keyed by the CLC code, read from a Norwegian GML delivery (field
  `clc18_kode`) with `lxml`. `:111-114`: classes above 500 are lakes (a
  "lake material"). **Carried:** the CLC code as the key, and "5xx is water"
  as the one rule the default map uses (R4). **Not carried:** GML, `lxml`,
  and the whole file parsed per call.
- `gml_repository.py:147-150` and `:161-163`: `constraints()` returns `[]`.
  The legacy never made land cover a constraint. `:185-225`: it labelled the
  **finished** mesh by testing each cell centre against every polygon in
  Python (`# TODO: Move to C++ for speed!`), and raised on a centre inside
  none. `legacy/rasputin/application.py:125-151` writes that as the face
  fields `cover_type` and `cover_color`. **Not carried in 16b:** labels are
  16c's (Q2), and a flood fill over unconstrained edges (Triangle's regional
  attributes) replaces the per-centre test, which is ambiguous for a triangle
  thinner than the snap.
- `gml_repository.py:141-145`: the data CRS from the GML's `srsName`, via
  `+init=`. Not carried: CRS comes from the GeoPackage's own
  `gpkg_spatial_ref_sys` (R3) through `crs.parse_crs`.

`@migration-expert` is not needed: nothing numeric is ported, and the one
rule carried (5xx is water) is the CLC nomenclature's own first digit
(`Legend/CLC_legend.csv` in Ola's download).
