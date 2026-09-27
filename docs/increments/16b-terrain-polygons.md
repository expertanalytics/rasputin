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
| 10 m | none | 0 | 0 | 132 563 | 34.5° | 0.02 % | 0.24° | — | — | 0.55 s |
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

### M5. The quarter circle on the committed tile, with CORINE

Increment 16's T-real domain (536 ring vertices, 30 km radius about
(850 250, 7 899 750) in `7908_3`), with the CORINE linework inside it (layout
A), quality start 25°, feet on. This is the case a committed extract makes
testable in CI (see "Test data").

| tol | features | input vertices | triangles | < 1° | worst | `node` |
|---|---|---|---|---|---|---|
| 1 m | none | 536 | 428 217 | 0.00 % | 0.40° | — |
| 1 m | CORINE | 5 627 | 441 863 (+3.2 %) | 0.08 % | 0.040° | 0.33 s |
| 10 m | none | 536 | 31 581 | 0.01 % | 0.29° | — |
| 10 m | CORINE | 5 627 | 54 840 (1.7×) | 0.04 % | 0.040° | 0.33 s |

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

- `legacy/rasputin/gml_repository.py:13-66`: CORINE's 44 classes as an
  `Enum` keyed by the CLC code, read from a Norwegian GML delivery (field
  `clc18_kode`) with `lxml`. `:124-127`: classes above 500 are lakes (a
  "lake material"). **Carried:** the CLC code as the key, and "5xx is water"
  as the one rule the default map uses (R4). **Not carried:** GML, `lxml`,
  and the whole file parsed per call.
- `gml_repository.py:181-182`: `constraints()` returns `[]`.
  The legacy never made land cover a constraint. `:184-226`: it labelled the
  **finished** mesh by testing each cell centre against every polygon in
  Python (`# TODO: Move to C++ for speed!`), and raised on a centre inside
  none. `legacy/rasputin/application.py:126-148` writes that as the face
  fields `cover_type` and `cover_color`. **Not carried in 16b:** labels are
  16c's (Q2), and a flood fill over unconstrained edges (Triangle's regional
  attributes) replaces the per-centre test, which is ambiguous for a triangle
  thinner than the snap.
- `gml_repository.py:155-159`: the data CRS from the GML's `srsName`, via
  `+init=`. Not carried: CRS comes from the GeoPackage's own
  `gpkg_spatial_ref_sys` (R3) through `crs.parse_crs`.

`@migration-expert` is not needed: nothing numeric is ported, and the one
rule carried (5xx is water) is the CLC nomenclature's own first digit
(`Legend/CLC_legend.csv` in Ola's download).

## The blueprint: data flow and boundaries

```
cli.py (composition root)                       flags -> FeatureRequest
  │
  ├─ dem_input.open_dem(DemRequest)              15a/15b, unchanged
  ├─ domain.read_domain(...).to_crs(dem crs)     16/15b, unchanged
  │
  ├─ feature_input.open_features(request, domain, dem_crs) -> FeatureSet     [16b-1]
  │     per source:
  │       GeoJSON: read the file (as domain.py does)
  │       GeoPackage: io.repository.open_geopackage(path)   read-only connection
  │                   io.geopackage.layer_info / query_features(conn, layer, box)
  │     region = domain boundary, densified, in the SOURCE crs, + margin      (R5)
  │     pre-clip each geometry to region (source crs)                          (R5)
  │     ClassMap: attribute value -> vocabulary names, or drop, or refuse      (R4)
  │     crs.reprojector(source, dem crs), vertex by vertex                     (R5)
  │     clip each ring or line to the domain polygon, as LINEWORK              (R6)
  │
  ├─ chains.start_chains(domain, features, vocabulary) -> StartChains       [16b-1]
  │     Outer, Holes (mask 0), then every clipped ring or line as a Breakline
  │     with its mask; closed iff unclipped (index repeated)                   (R7)
  │
  ├─ _engine: build_pslg -> node -> triangulate   noder merges doubled edges  (R7)
  │                                                verifier indexed           (R8, 16b-0)
  └─ refine, trim, write: unchanged code; new fields                           (R10)
```

Boundaries, each one checkable by an import test (`TestModuleIsolation`
style):

- **No path below `cli.py` except** `domain.py` and `feature_input.py` (which
  read GeoJSON text, as `domain.py` already does), `dem_input.py` (which only
  builds a repository) and `io/repository.py` (which opens files). The
  GeoPackage decoder in `io/geopackage.py` takes an open `sqlite3.Connection`,
  never a path: the connection is its stream.
- **No CRS in C++.** Feature geometry reaches `_core` as `float64` `(x, y)` in
  the DEM's CRS, like the domain.
- **No `_core` in** `feature_input.py`, `chains.py`, `io/geopackage.py`: all
  three are testable with no compiled extension. Roles leave `chains.py` as
  the strings `"outer"`, `"hole"`, `"breakline"`, and `cli.py` maps them onto
  `ChainRole` through its existing `ROLES` table.
- **`features.py` keeps importing nothing first-party** (increment 7). The
  class maps live in `feature_input.py`, which imports it.
- **Declarative:** `FeatureRequest` is frozen data. `open_features` is a
  function of the request, the domain and the files, and returns frozen data.
  An async caller runs it in `asyncio.to_thread` (sqlite3 and GEOS block).

## Rulings

### R1. The input, and the flags

- **`--features PATH`**, one source, by suffix: `.geojson`/`.json` is a
  GeoJSON `FeatureCollection`; `.gpkg` is a GeoPackage. Anything else is
  refused naming the suffix. Several sources in one run are the API's
  (`FeatureRequest.sources` is a tuple) and a later CLI change (Q4).
- **Only with `--domain`** (and so only with `--tolerance`, 16 R1). Features
  are clipped to the confining polygon (R6); the stride-grid path has none.
  Refused otherwise, as a usage error.
- **`--features-crs TEXT`**, anything `crs.parse_crs` reads. GeoJSON: the
  file's `crs` member, else EPSG:4326 (RFC 7946), exactly `domain.py`'s rule.
  GeoPackage: the layer's `srs_id` through `gpkg_spatial_ref_sys` (R3). Given
  and disagreeing with the file's own (by `pyproj.CRS` equality) is refused,
  as for `--domain-crs`.
- **`--features-layer NAME`**, GeoPackage only: the features table. Required
  when the file has more than one `features` row in `gpkg_contents` (CORINE
  has six: Europe and five overseas departments); the only one otherwise.
  Refused for GeoJSON.
- **`--features-map NAME`**, one of the built-in class maps (R4): `property`
  (the default), `corine`, `corine-water`. A map read from a file is Q4.
- **Accepted geometry**: `Polygon`, `MultiPolygon`, `LineString`,
  `MultiLineString`. Z and M are dropped (`shapely.force_2d`). A `Point`,
  `MultiPoint` or `GeometryCollection` is refused naming the feature: input is
  geometry, never loose points (Ola's direction 1). An empty geometry, or a
  GeoPackage blob with the empty flag, is skipped and counted.
- **Feature polygons are not validated.** Their rings enter as linework (16
  R6: no role, no winding), and a self-intersecting ring is a breakline the
  noder nodes like any crossing. Only the confining polygon is validated
  (16 R1). This is deliberate: refusing CORINE because one of its 2.4 M
  polygons is invalid would help nobody, and nothing downstream reads a
  feature as an area.

### R2. Homes

| module | holds | imports |
|---|---|---|
| `io/geopackage.py` (new) | `GpkgLayer` (frozen: table, geometry column, `srs_id`, CRS text, R-tree table or `None`); `layer_info(conn, table)`; `decode_geometry(blob) -> BaseGeometry`; `query_features(conn, layer, box, attribute) -> Iterator[RawFeature]` | `sqlite3` (types only), `struct`, shapely, pydantic |
| `io/repository.py` | adds `open_geopackage(path) -> sqlite3.Connection`, read-only | as today, plus `sqlite3` |
| `tin_engine/feature_input.py` (new) | `ClassMap`, the built-in maps; `FeatureSource`, `FeatureRequest`; `TerrainFeature`, `FeatureSet`; `FeatureError`; `open_features(request, domain, dem_crs) -> FeatureSet` | json, shapely, numpy, pydantic, `crs`, `domain`, `features`, `io.repository`, `io.geopackage` |
| `tin_engine/chains.py` (new) | `StartChains`; `start_chains(domain, features, vocabulary)`: 16 R6's table as code. `cli._domain_chains` moves here and is its first half | numpy, shapely, `domain`, `feature_input`, `features` |
| `tin_engine/features.py` | two more vocabulary entries (R4, Q3) | unchanged: nothing first-party |
| `cli.py` | four options, refusals, one call each to `open_features` and `start_chains`, the fields (R10) | adds `feature_input`, `chains` |
| `include/terrain/noding/noded_pslg_builder.hpp` | the verifier's sort-and-sweep (R8, 16b-0) | unchanged |

**15's Q3 stays true as written:** `io/repository.py` is still "the one
module in `io/` that opens files, read-only". SQLite cannot read from a
Python stream, so the GeoPackage decoder's stream is the open connection.
`io/__init__.py`'s rule gains one clause saying so.

### R3. Reading a GeoPackage without GDAL

- **Open read-only**: `sqlite3.connect(path.resolve().as_uri() + "?mode=ro",
  uri=True)`, closed with `contextlib.closing` (sqlite3's own context manager
  commits, it does not close). Read-only also means the file's triggers,
  which need SpatiaLite (M1), never fire.
- **Check it is a GeoPackage**: `PRAGMA application_id` is `0x47504B47`
  (`GPKG`; Ola's file answers `1196444487`, which is that). Otherwise
  refused. `user_version` is recorded, not checked (Ola's is `10200`, 1.2).
- **The layer**: its `gpkg_contents` row (`data_type = 'features'`), its
  `gpkg_geometry_columns` row (column name, `srs_id`), and its
  `gpkg_spatial_ref_sys` row. The CRS text is `EPSG:<organization_coordsys_id>`
  when `organization` is `EPSG` (case-insensitive), else the row's
  `definition` WKT; either goes through `crs.parse_crs`. Ola's file: srs 3035,
  organization `EPSG`.
- **The spatial filter**: when `gpkg_extensions` lists `gpkg_rtree_index` for
  the layer and `rtree_<table>_<column>` exists, the query joins it on the
  box (`maxx >= ? AND minx <= ? AND maxy >= ? AND miny <= ?`, closed). SQLite
  stores R-tree bounds as 32-bit floats rounded outward, so this is a
  superset, which is all a candidate query needs. Without an index the table
  is scanned and the report says so. **Rows are ordered by the primary key**
  (`ORDER BY <pk>`), because an R-tree join returns rows in no stated order
  and the mesh depends on input order (R7).
- **The blob header** (OGC 12-128, 2.1.3), `decode_geometry`:
  - bytes 0-1 `GP`, byte 2 version `0`; else refused;
  - byte 3 flags: bit 0 the header's byte order (for `srs_id` and the
    envelope), bits 1-3 the envelope code (0: none; 1: 32 bytes; 2 and 3:
    48; 4: 64; 5-7 refused), bit 4 empty (skip and count), bit 5 extended
    geometry type (refused: not a standard WKB geometry);
  - bytes 4-7 `srs_id`, which must equal the layer's; else refused, naming the
    row;
  - the rest is WKB for `shapely.wkb.loads`, which reads either byte order.
- **Refusals name the file, the layer and, for a row, its primary key.**
- **The SQLite R-tree module must be compiled in.** It is in the venv's
  Python 3.14 (SQLite 3.53.4: `CREATE VIRTUAL TABLE ... USING rtree` works);
  whether it is in CI's Python 3.12-3.14 on `ubuntu-latest` is not verified.
  If not, the query fails with `no such module: rtree`; the reader turns that
  into a refusal naming the cause rather than falling back silently, and the
  test suite's in-memory GeoPackages detect the module and say so.

### R4. From feature attributes to edge bits: class maps

```python
class ClassMap(BaseModel):                 # frozen, extra="forbid"
    name: str
    attribute: str                         # e.g. "property", "Code_18"
    classes: Mapping[str, tuple[str, ...]] # value -> vocabulary names
    otherwise: Literal["refuse", "drop"] | tuple[str, ...] = "refuse"
    notice: str = ""                       # attribution to carry into the file (R10)
```

- **A value maps to vocabulary names**; the edge mask is
  `vocabulary.mask(*names)` (increment 7: no bare masks). `()` is a
  constraint with no bits. `"drop"` means the feature is not a constraint at
  all. `"refuse"` refuses the run, naming the feature and the value.
- **Every name in a map is checked against the vocabulary when the run
  starts**, through `mask`, which raises on a name it does not know. A map
  can therefore not smuggle a bit number in.
- **Built-in maps:**
  - `property`: attribute `property`; each vocabulary name maps to itself;
    `otherwise = "refuse"`. A feature's `property` may be one name or a list
    of names. This is 16 R6's "each with a `property` naming a vocabulary
    entry".
  - `corine`: attribute `Code_18`; `511`, `512`, `521`, `522`, `523` map to
    `("land_cover", "water")`; `otherwise = ("land_cover",)`. Every CLC
    boundary is a constraint. The rule is the nomenclature's first digit
    (5 = water bodies; `Legend/CLC_legend.csv`), which the legacy used too.
  - `corine-water`: the same five codes map to `("water",)`;
    `otherwise = "drop"`. Only shorelines, river banks and the coast.
- **A shared edge gets the union** of its two features' masks: the noder's
  existing merge (increment 7, "Union, not priority"). A CORINE edge between
  a lake and a forest carries `land_cover | water`; so does one between a
  lake and a river. **"Water on exactly one side" is not an edge bit.** It
  needs to know which side is which, which is a region label (16c, Q2), not a
  property of the edge. Under `corine-water` an edge is kept when at least
  one side is water, because the kept feature's ring carries it.
- **Why not one bit per CLC class.** 44 classes do not fit the 32-bit word,
  and increment 7 ruled land cover "stays face-based". The bits say what
  kind of line an edge is; which class lies on each side is 16c's.
- **Vocabulary additions** (Q3): `land_cover` at bit 7 and `water` at bit 8.
  `coastline` (bit 3) stays for linear coastline data; CORINE's sea
  boundary is `water`, because a lagoon or estuary boundary is water too.
  The fingerprint changes; no stored artifact carries one yet (R10 is the
  first), so nothing is invalidated.

### R5. CRS, and what is read

- **Features go into the DEM's CRS by `crs.reprojector`**, vertex by vertex,
  as `DomainPolygon.to_crs` does (15b R9): edges straight in the DEM's CRS,
  one `from_crs` site, `always_xy`. A vertex with no image (`inf`) is
  refused, naming the feature. `features_crs` and `features_transform` are
  recorded (R10).
- **The region.** The domain (already in the DEM's CRS) has its boundary
  densified to at most 1 km per segment and transformed into the source CRS;
  its convex hull, buffered by 100 m, is the region. The GeoPackage query box
  is the region's bounds. The 100 m is a safety margin, not a tolerance:
  a 1 km chord bends between two conformal or equal-area projections by
  orders of magnitude less (15 B7 measured millimetres per kilometre within a
  UTM zone), and anything the margin lets in is clipped exactly in R6.
- **The pre-clip.** Each candidate geometry is intersected with the region
  in the source CRS before it is reprojected. This bounds the work by the
  domain, not by the feature: CORINE's sea polygon is 360 000 vertices (M1),
  of which a catchment needs a few thousand; and it keeps vertices far from
  the domain out of a transform that may have no image for them. The
  pre-clip only removes geometry the exact clip would remove; I4 pins that
  it changes nothing inside the domain.
- **GeoJSON sources** are read whole (`json.loads`, as `domain.py`), then
  filtered by `intersects(region)` and pre-clipped the same way. Streaming
  large GeoJSON is not needed at this scale (the 2026-09-27 probe's CORINE
  over `7908_3` as GeoJSON, 88 features in EPSG:25833, is 8 MB) and not
  designed.

### R6. Clipping to the confining polygon: linework, not areas

This is the ruling 16 R6 left to 16b ("a feature leaving the confining
polygon is 16b's ruling").

- **Each polygon ring and each line is clipped as a line**:
  `shapely.intersection(ring_as_linestring, domain.polygon)` in the DEM's CRS,
  domain holes included. What is outside the domain is dropped; a feature
  wholly outside is dropped and counted.
- **Not an area clip.** `polygon ∩ domain` would give each clipped polygon a
  stretch of the domain's own boundary as a ring edge, so the domain boundary
  would carry `land_cover` bits wherever a polygon was cut. Clipping rings as
  lines gives open pieces that end on the boundary and adds no edge along it.
- **A piece that lies along the domain boundary** (a feature edge
  coinciding with it, possible when the domain was cut from the same data)
  is kept. The noder merges it with the `Outer` or `Hole` edge and the union
  gives that edge the feature's bits. That is true of the data, so it is
  right.
- **Crossing points** are GEOS's doubles, within rounding of the domain edge.
  The noder snaps them and splits the ring there (05b step 4), as for any
  crossing. No special case.
- **Point pieces** (a ring touching the boundary at one point) are dropped:
  a point is not geometry here.
- **An unclipped ring enters closed; a clipped one as open pieces.** A ring
  wholly inside comes back from GEOS as one closed line (first coordinate
  equals last), which `start_chains` writes as a closed breakline by repeating
  the first index (increment 8, ruling 2). GEOS may start it at a different
  vertex; that is harmless.
- **Extent.** After the clip every vertex is inside the domain's closure up
  to GEOS rounding, and the domain is already inside the DEM's node
  rectangle (16 R1). A crossing point on a domain edge that lies on the
  rectangle's border rounds to the noder's 1 mm grid like every other vertex;
  a DEM whose border is not on that grid can put it up to 0.5 mm outside,
  which `refine` refuses (`OutsideGrid`). That risk already exists for the
  domain's own vertices under 16 U6 and is not made worse; one fixture pins
  the border-on-grid case.
- **Increment 8's gallery rows are the fixtures**, as GeoJSON features over
  a synthetic DEM: `road-enters-forest` (a road ending inside a closed ring),
  `wall-leaves-domain` (a line crossing the outer ring: on this path the
  exterior half is clipped away before the engine, so the mesh has no
  finding), and `bridge-over-lake` (a road crossing a ring twice). The
  gallery's own fixtures are unchanged: they test the engine, not the clip.

### R7. Shared boundaries: the noder merges them

- **Rings go to the engine as given**, clipped (R6), one chain per ring or
  piece: 16 R6's letter. Every CORINE interior edge therefore arrives twice,
  in opposite orientation.
- **The noder already does the right thing** (M3): it merges the two copies
  into one edge by the node-id key and unions their masks (05b guarantee
  14(a), amended; "the merge happens at the node-id edge-key dedup"). Layouts
  A and B gave identical noded edge sets and masks on real data.
- **Rejected: deduplicate in Python first** (layout A). It is exact on CORINE
  (M2: shared endpoints are bit-identical after reprojection) and halves the
  noder's input. But it is a second merge site with its own union to get
  right, about 30 lines, and after R8 it saves a constant factor, not an
  order. It is the first thing to add if @perf finds `node` significant after
  R8.
- **Rejected: `shapely.unary_union` of the linework.** GEOS's noding may snap
  vertices, and it drops per-edge attributes.
- **Same-class neighbours** (CORINE's "redundant lines", readme) stay
  constraints: the union of two equal masks is that mask, and dropping the
  edge would need to know both sides, which R7 deliberately does not compute.
- **Order is the source's**: features by primary key (GeoPackage) or file
  order (GeoJSON); within a polygon, exterior then holes; within a
  `MultiPolygon`, part order. Node ids are a function of the node set
  (05b), but refinement's result depends on vertex order (M3: layouts A and
  B, whose start meshes are identical, refined to 261 377 and 261 375
  triangles), so a stated order is what makes a run reproducible.

### R8. The noder's verifier, indexed [16b-0]

M3's quadratic makes a 48 km CORINE square spend 25 s in `node`, and a
catchment the size of Glomma's (about 42 000 km², 18× that square) would
spend hours. The verifier's premise, "the candidates are small by
construction", is false for land cover.

- **Step 0, before any code: isolate the cost.** Build with
  `check_guarantee_14` returning early, rerun `clc_nodetime.py`'s sizes, and
  confirm the quadratic goes with it. If it does not, the cost is in steps
  1-5, and `@architect` re-scopes 16b-0 before `@tester` starts.
- **The replacement: sort and sweep on x.** Items are the edges (their
  bounding boxes) and the node cells (each node's closed cell box, padded by
  one grid spacing on every side: false positives are free, and the padding
  removes any doubt about how the cell's bounds round). Sort by `(xmin,
  kind, index)`. Keep the active edges whose `xmax >= ` the current `xmin`
  (closed). Test an arriving edge against every active edge, and an arriving
  cell against every active edge, whose closed y-intervals intersect; and an
  arriving edge against every active cell likewise. The tests are the
  existing `classify<K>` and `segment_meets_cell<K>`, unchanged.
- **Completeness in one line each:** two segments that cross or overlap
  share a point, so their closed boxes intersect; a segment that meets a
  closed cell meets its box. Box bounds are `min`/`max` of the endpoints'
  doubles, which is exact.
- **Independent of `BroadPhase`.** No shared code and a different algorithm,
  so 5b's reason for having no index ("a verification pass sharing the index
  would be blind wherever the index is") still holds: the verifier does not
  share the driver's index.
- **What stays the same:** the statuses, and which inputs are refused. **What
  may change:** which violating pair the message names first. No test asserts
  wording (house rule).
- **Worst case** is still quadratic (every x-interval overlapping every
  other, such as many long parallel east-west lines); typical inputs are
  `O(n log n + k)`. The header comment says so, replacing the false premise.
- **The brute force stays, as the test oracle**, in `tests/cpp/support/`.

### R9. Vertex density, z, and refinement

- **Vertices are used as given** (Ola's direction 1): no simplification, no
  minimum segment length, no merging of close vertices beyond the noder's
  1 mm snap (16 U6). Measured: no two CORINE vertices within 1 mm in three
  tiles, 52 pairs within 1 cm in one (M2). A pair within the snap would merge
  into one node; the report shows input and noded vertex counts, so it is
  visible.
- **Segments shorter than a cell** (9 000-16 000 per inland tile, the
  shortest 4 mm) enter unchanged. They cost angles, not correctness: the
  worst angle falls to 0.0007° (M4) because a 4 mm segment is a 4 mm local
  feature (Ruppert 1995). Coarsening is the later increment's; CORINE's exact
  partition (M2) is what a coverage-preserving simplifier needs.
- **z is 16 R2's**: bilinear at an off-node vertex, `value_at` at a node.
  Nothing about a feature changes its height: a lake's shore vertices get the
  DEM's bilinear height, not the lake level. Hydro-flattening is out of
  scope.
- **NoData**: a feature vertex in a cell with a NoData corner is invalid,
  and its triangles are carved (16 R2). Unchanged.
- **Refinement code is unchanged.** Interior constraints are what 14b R5 and
  20b R8 already cover: never flipped, split with bits kept, feet on them
  (M4: 30 045 feet at 1 m, none refused). The tolerance guarantee is
  unchanged: error is measured at DEM nodes only.
- **The quality start (20) is unchanged, and M4 is its input for 20c.** On
  CORINE it inserts 154 640 nodes and triples the 10 m mesh. 16b does not
  change it; Q5 asks whether the acceptance runs should record both settings.

### R10. What the file records, and the report

- **`features`**: source file name, layer (GeoPackage), map name, and counts
  after the clip: e.g. `U2018_CLC2018_V2020_20u1.gpkg:U2018_CLC2018_V2020_20u1,
  map corine, 727 features (12 dropped outside), 988 chains, 51 395
  vertices`. `.vtk` field and `.ply` comment, like `domain`.
- **`features_crs`** and **`features_transform`**, as `domain_crs` and
  `domain_transform` (15b).
- **`edge_vocabulary`**: `DEFAULT_VOCABULARY.fingerprint()` and the names by
  bit (`0 river, 1 road, ..., 7 land_cover, 8 water`). This is increment 7's
  mechanism 2 ("a field of whatever serialized artifact carries a mesh"),
  which no writer has implemented yet: without it the edge masks in the file
  cannot be read. Written whenever the file has an edge block, features or
  not.
- **`features_notice`**: the map's `notice`, when it has one. The CORINE maps
  carry the Copernicus attribution (see "Test data"), because a mesh built
  from CORINE and shared must say so and say it was modified.
- **`elevation_source`**'s start clause becomes
  `start domain boundary and features, vertex z bilinear` when features are
  given.
- **stderr**: features read, dropped outside, clipped, and input versus
  noded vertex counts. `--stats`' phase table gains `features read` and
  `features clip`.

### R11. Scale and memory

- **Per request, not per dataset.** The 8.9 GB GeoPackage is never loaded:
  the R-tree gives candidates in 11-95 ms per tile-sized region (M1), and
  SQLite reads only the pages it touches.
- **Python memory** is the candidates' geometry: about 1 M raw vertices per
  tile-sized region before the pre-clip (16 MB of coordinates plus shapely's
  overhead), about 120 000 after. The probe's peak RSS for a whole tile was
  330-340 MB with nothing else loaded. The DEM dominates a real run (a
  50 km tile is 100 MB of float32; M4's 1 m run peaked at 2.4 GB for 6.6 M
  triangles).
- **The engine** holds the chain positions (about twice the unique vertices,
  R7) and the noded graph: tens of MB for a catchment.
- **A whole-Norway request is not a use case.** A coarse box around Norway
  holds 424 319 CORINE features and 850 MB of blobs. The use case is a
  catchment, which auto-catchment (next in the order of work) will supply.
- **Until R8 lands, `node` is the bottleneck** (M3, M4). After it, @perf
  measures `node` on the 48 km square and on a larger block.

## Invariants

- **I1. Bits, from the input.** For every noded constraint edge, its mask
  contains the union of the masks of every feature whose clipped linework
  contributes geometry to it (05b's guarantee 15, applied to features). On a
  partition, an interior edge's mask is exactly the union of its two
  features' map masks. The oracle is built from the features and the map,
  never from the noder's provenance (05b, "Guarantee 15's oracle").
- **I2. The clip keeps input vertices exact.** Every vertex of the engine's
  feature input is either a reprojected source vertex, bit for bit, or lies
  within 1e-6 m of the domain's boundary (a crossing). None lies outside the
  domain's closure by more than that.
- **I3. Determinism.** The same request gives a bit-identical file whatever
  order SQLite returns candidate rows in.
- **I4. The pre-clip changes nothing inside the domain.** The engine input is
  identical with the pre-clip and without it.
- **I5. Ring shape.** A ring wholly inside the domain enters as one closed
  `Breakline` (first index repeated); a clipped ring as open `Breakline`s; a
  line as open `Breakline`s. No feature chain is `Outer` or `Hole`.
- **I6. The verifier's verdict is the brute force's.** For every input the
  sweep returns the status the brute-force check returns.
- **I7. The tolerance guarantee is unchanged**: every valid DEM node in the
  domain is within tolerance of the output mesh (14 T3's oracle), with
  features.

## Degeneracy policy

| case | outcome |
|---|---|
| a feature edge collinear with a domain edge, fully or partly | kept; the noder merges or splits it and unions bits onto the domain edge (R6) |
| a ring touching the domain boundary at one point | the point piece is dropped; the ring's pieces end there |
| a feature wholly outside the domain, or inside a domain hole | dropped, counted |
| no feature left after the clip | meshed with the domain alone; the report says 0 features (a catchment inside one land-cover polygon is legitimate) |
| two feature vertices within the 1 mm snap | one node (16 U6); input and noded vertex counts in the report |
| a closed feature ring collapsing to two nodes after the snap | the noder's existing handling of a collapsed breakline; one fixture pins the status |
| a self-intersecting feature ring | linework; the noder nodes the crossing (R1) |
| the same edge in two features, same class | one constraint edge, same mask (R7) |
| a partial overlap or T-junction between two sources | the noder splits and merges, as for any input (05b) |
| `Point`, `MultiPoint`, `GeometryCollection` | refused, naming the feature (R1) |
| an empty geometry or empty-flag blob | skipped, counted |
| an unknown attribute value | per the map: refused (naming the feature and value), dropped, or the map's default (R4) |
| a vertex with no image under the transform | refused, naming the feature (R5) |
| a GeoPackage without an R-tree | full scan, reported (R3) |
| SQLite without the R-tree module | refused, naming the cause (R3) |
| a blob with `srs_id` other than the layer's, an extended type, a bad magic or version, envelope code 5-7 | refused, naming the row (R3) |
| not a GeoPackage (`application_id`) | refused (R3) |

## Not in scope

- Region labels per triangle (16c, Q2), and "water on exactly one side".
- Coarsening or simplifying feature geometry.
- Hydro-flattening lakes or rivers.
- Several `--features` in one CLI run, and a class map read from a file (Q4).
- Grade separation (increment 8, ruling 1: never).
- Python-side deduplication of shared edges (R7; kept in pocket).
- Streaming GeoJSON; reading GML, Shapefile, FlatGeobuf or any other format.
- Changing the quality start or feet (20c).

## Sub-increments and LOC

Counted in `CLAUDE.md` §2's unit. Estimates, then the worst overrun seen so
far (+39 %, increment 16) applied.

| | file | what | est. |
|---|---|---|---|
| **16b-0** | `noding/noded_pslg_builder.hpp` | sort-and-sweep for 14(a) and 14(b), cell padding, header comment | ~50 |
| | | **16b-0 total** | **~50 (worst ~70)** |
| **16b-1** | `io/geopackage.py` | `GpkgLayer`, `layer_info`, `decode_geometry`, `query_features` | ~70 |
| | `io/repository.py`, `io/__init__.py` | `open_geopackage` | ~8 |
| | `feature_input.py` | `ClassMap` and three built-in maps, request and result models, GeoJSON reading, region, pre-clip, reprojection, clip, counts, refusals | ~130 |
| | `chains.py` | `start_chains` (the domain half moved from `cli.py`, net ~+25) | ~45 |
| | `features.py` | two vocabulary entries | ~2 |
| **16b-2** | `cli.py` | four options and their refusals, the two calls, `features`, `features_crs`, `features_transform`, `features_notice`, `edge_vocabulary`, the sentence, the report; less `_domain_chains` | ~65 |
| | `stats.py` | two phase rows | ~4 |
| | | **16b-1 + 16b-2 total** | **~325 (worst ~450)** |

**Two PRs, recommended:**

1. **16b-0 alone.** It is C++ in the noder's invariant-critical verifier, it
   stands on its own (every noder input benefits), its test is an
   equivalence against brute force, and it is the one piece whose premise
   (M3's attribution) must be confirmed before anything is written (R8 step
   0). About 50 lines.
2. **16b-1 and 16b-2 together.** 16b-1 alone would ship a reader nobody can
   run; the CLI is what makes the acceptance runnable, the same argument
   increment 16 made against splitting. About 325, worst about 450, under
   700.

If 16b-2 grows past the ceiling in review, the seam is `feature_input.py`
plus `io/geopackage.py` (testable from Python with no CLI) first, `chains.py`
and `cli.py` second.

## Tests for `@tester`

**Invariant-critical (mutation testing required): two suites.**

- **`prop_noding_verifier_sweep` (C++, 16b-0)**, I6. Random PSLGs through
  the builder in both forms, the sweep in production and the brute force
  from `tests/cpp/support/`, asserting equal status. The generators must
  include what the sweep can get wrong: duplicate edges in both
  orientations; collinear partial overlaps; boxes that touch at exactly
  equal `x` (an edge's `xmax` equal to another's `xmin`); vertical and
  horizontal edges (zero-width boxes); a node whose cell touches an edge's
  box only at a corner; long edges spanning every other box. Plus the
  existing 14(a) and 14(b) refusal fixtures, unchanged.
  Mutants: `xmax > xmin` for `>=` in eviction; open y-interval test; cells
  unpadded and cell box computed from the node point instead of its cell;
  arriving edges not tested against active cells; the active list evicted
  one item late or early; sort by `xmax` instead of `xmin`.
- **`test_feature_chains_bits` (Python, 16b-1)**, I1, I2 and I5, with the
  oracle built from the input features and the map. On synthetic partitions
  (a 2×2 grid of squares with one class each, a square with a hole holding a
  second polygon, a lake with a river ending on its shore) and on the
  committed CORINE extract through the real `_engine`.
  Mutants: masks taken from the first contributing feature only (no union);
  a clipped ring closed by repeating its first index; holes' rings dropped;
  area clipping instead of line clipping (the domain boundary then carries
  bits: caught by I1's "exactly" on a partition, and by a direct assertion
  that no `Outer` edge carries bits where no feature edge lies on it);
  `corine`'s 5xx rule off by one class (`52x` treated as land only).

Everything else is unit, property or integration testing:

- **GeoPackage decoding**: blobs with every envelope code 0-4, big- and
  little-endian headers, the empty flag, the extended flag (refused), a bad
  magic, version 1, envelope code 5, an `srs_id` differing from the layer's.
  Built in memory: a test helper writes a minimal GeoPackage (the four
  `gpkg_*` tables, one feature table, its R-tree) with `sqlite3` to
  `tmp_path`, so the reader opens a real file. The helper is test support.
- **Layer resolution**: one table (default), several (refused without
  `--features-layer`), an unknown name, `data_type` other than `features`,
  no R-tree (full scan, reported), a non-EPSG `organization` (WKT path),
  `application_id` wrong.
- **Class maps**: `property` with one name, a list, an unknown name
  (refused), no `property` (refused); `corine` on every 5xx code and on
  `111` and `999`; `corine-water` dropping non-water; a map naming a name the
  vocabulary lacks (refused when the run starts).
- **CRS**: a GeoJSON in EPSG:4326 over a 25833 synthetic DEM; a GeoPackage
  in 3035; `--features-crs` agreeing and disagreeing; a vertex with no image.
- **Region and pre-clip**: I4 on the extract (engine input with and without
  the pre-clip, compared bit for bit).
- **Clip degeneracies** from the table above, one fixture each.
- **Determinism**, I3: the extract's rows returned in a shuffled order (a
  second GeoPackage written with rows inserted in reverse) gives a
  bit-identical `.vtk`.
- **CLI**: `--features` without `--domain`; an unknown suffix;
  `--features-layer` on GeoJSON; the four new fields and the notice in `.vtk`
  and `.ply`; `edge_vocabulary` present on a run without features too.
- **Increment 8's three rows** as GeoJSON over a synthetic DEM, end to end:
  the road ends inside the forest ring; the wall's exterior half is absent
  and the mesh has no edge outside the domain; the bridge road is split at
  two shore crossings with `road` and `water` bits on the right edges.
- **I7 on the extract**: every valid node of the committed tile inside the
  quarter circle within tolerance, at 10 m (fast), with CORINE.

## Test data

- **Synthetic**: GeoJSON documents built in the tests; minimal GeoPackages
  written by the helper above; the micro-DEMs increment 16's suites already
  use.
- **Committed: a CORINE extract over the committed tile `7908_3`.** An
  extraction script under `tests/fixtures/corine/` (test support, as
  `tests/fixtures/dtm10/extract.py` is) reads Ola's GeoPackage, takes the
  features whose R-tree box meets the tile's 3035 image, clips them in 3035
  to that image buffered by 1 km (so the partition is exact over the tile),
  and writes them to a small GeoPackage with the source's table name,
  columns, SRS row and R-tree, blobs in the source's header layout (flags
  `0x01`). Expected size: about 11 000 vertices over the tile (M2), a few
  hundred kB. It makes the quarter circle with CORINE (M5) a CI test.
- **Local only, for @perf**: Ola's GeoPackage and the 48 km square in
  `6603_4` (M4), and a larger block once R8 lands.
- **Attribution.** The GeoPackage's own metadata
  (`Metadata/U2018_CLC2018_V2020_20u1.xml` in Ola's download) grants "full,
  open and free access" under Regulation (EU) No 1159/2013, on four
  conditions: inform the public of the source; do not suggest EU endorsement;
  state it clearly where the data "has been adapted or modified"; and
  acknowledge that the data were produced "with funding by the European
  Union". The extract is modified (clipped), so its `README` says, at least:
  source (CORINE Land Cover 2018, version 2020_20u1, Copernicus Land
  Monitoring Service, European Environment Agency), that it is clipped and
  re-encoded, the funding sentence, and no endorsement. The `corine` maps'
  `notice` carries a one-line form of the same into every mesh built from
  them (R10). **Not verified:** the exact attribution line the Copernicus
  Land Monitoring Service website asks for today (recalled as "© European
  Union, Copernicus Land Monitoring Service 2018, European Environment Agency
  (EEA)"); no web access in this round. Ola or the main session should check
  `land.copernicus.eu` before the extract is committed.
- **Found in passing:** `tests/fixtures/corine/0000_4326_corine2018_4e6064_GML.gml`
  (30 MB, 399 CORINE features over Norway in EPSG:4326, written by OGR,
  committed with the foundation reset `3096ccc`) carries no attribution. Only
  `legacy/tests/test_gml_repository.py` names that directory. 16b's PR adds the
  same attribution for it, or removes it if Ola prefers (Q6), per the rule
  that a documentation defect found in an increment is fixed in its PR.

## Acceptance

- **16b-0**: R8 step 0's isolation recorded; `prop_noding_verifier_sweep`
  green with its mutants killed; `node` on M3's 48 km square in `6603_4`
  (layout B, the one 16b uses) far below its 83 s, with the ratio across
  M3's sizes showing the quadratic gone. Every existing noder and gallery
  test green; the gallery pictures unchanged.
- **16b-1/2**: the quarter circle on the committed tile with the committed
  extract, `--features-map corine`, at 1 m and 10 m, writes a `.vtk`
  ParaView opens, with the fields of R10; constraint edges carry
  `land_cover` and, on the coast, `water`; the tolerance holds; M5 is the
  reference. Then Ola's own case: a catchment in EPSG:4326 over the DTM10
  archive with Ola's GeoPackage.
- **@perf, for both PRs.** README's acceptance rule applies: 16b-0 changes
  what the start mesh costs, and 16b-2 changes what `_dem_mesh`, the code
  that drives `refine`, feeds it. The 1 m benchmark and thread sweep with
  `tools/bench.py` must show no change on the standard inputs (they have no
  features; their start ring is small). In addition, @perf records the first
  **CORINE baseline**: M4's 48 km square at 1 m and 10 m, `node` and
  `refine` phases, triangles, peak memory, power state, under
  `docs/benchmarks/<date>/`. The quality-start tripling (M4) is recorded
  there as 20c's input.
- All gates in `CLAUDE.md` §4, including the TSan job; CI green.

## Questions for Ola

**Q1. Where may the noder's verification spend its time?** M3: 25 s of a
30 s run is the noder on a 48 km square of CORINE, and it grows with the
square of the edge count.
- **(a) Index the verifier (16b-0, R8), keep it on for every run.
  Recommended.** About 50 lines of C++, same statuses, and the brute force
  stays as its test oracle.
- (b) Skip verification in Release builds. Fastest, but 5b made the
  verification the thing that turns a missed T-junction into a refusal
  instead of a wrong mesh; switching it off gives that up.
- (c) Leave it; features on small catchments only. Not viable for the
  auto-catchment work that comes next.

**Q2. Region labels: which triangles lie in which polygon.** 16b makes every
land-cover boundary a constraint edge, but the mesh does not say which class
each triangle has, and "shoreline" (water on one side only) needs that.
- **(a) A separate increment 16c, next after 16b: a label per triangle by
  flood fill across unconstrained edges from one seed per feature, the method
  of Shewchuk's Triangle (`-A`), written as a cell field `land_cover` with the
  CLC code. Recommended.** About 80 lines, pure Python over the trimmed mesh
  and the features; the legacy wrote the same field by testing cell centres.
- (b) Fold it into 16b. Pushes 16b-1/2 to about 400 lines (worst 560),
  still under 700, but adds a second invariant-critical suite to the round.
- (c) Not needed; the edges are enough.

**Q3. Two new vocabulary entries, `land_cover` (bit 7) and `water` (bit 8).**
- **(a) Both, as R4. Recommended.** The CORINE maps need a bit for "a
  land-cover boundary" and one for water; `coastline` stays for linear
  coastline data.
- (b) One bit, `land_cover`, and water left to 16c's labels.
- (c) Reuse `coastline` for every water boundary. Wrong for lakes and rivers.

**Q4. The CLI's reach in 16b.**
- **(a) One `--features` source and the three built-in maps; several sources
  and maps read from a file later. Recommended.** The API already takes
  several sources; the CLI question is how flags pair with sources, which
  deserves its own small design once a second source (roads, rivers) is on
  the table.
- (b) Repeatable `--features` now, each followed by its own
  `--features-map`, paired by position. About 20 more lines, and
  position-paired flags are easy to get wrong.
- (c) A JSON request file (`--features-spec`) listing sources and maps. The
  most general; about 40 more lines.

**Q5. The quality start on CORINE.** It triples the 10 m mesh (M4: 211 487
triangles without it, 480 127 with it), because dense boundaries make many
bad start triangles.
- **(a) Change nothing in 16b; @perf records both settings, and 20c decides
  with this data. Recommended.**
- (b) Turn the quality start off when features are given. Quick, but it
  decides 20c's question by the back door.

**Q6. The 30 MB CORINE GML already in `tests/fixtures/corine/`** (legacy
test data, written by OGR, no attribution).
- **(a) Keep it and add the attribution in 16b's PR. Recommended**: removing
  data from the tree is your call, and the legacy suite names the directory.
- (b) Remove it: nothing outside `legacy/` uses it, and 16b's extract
  replaces it for everything current.
