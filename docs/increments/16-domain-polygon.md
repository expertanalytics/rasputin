# Increment 16 — mesh a domain polygon, starting from its boundary alone

Status: **designed, not started.** The user chose U1 (a), U3 (a), U4 (a),
U5 (a) and U6 (a) on 2026-09-26, with a standing caveat on U1: the domain's CRS
will not always have to match the DEM's (section "Ruled by the user"). Written by `@architect` before `@tester`, per
`docs/increments/README.md` step 1, on branch `increment16-domain-polygon` off
`a58e0ab` (increment 14b). Revised the same day after two user rulings: DEM
values are point heights (R0), and input vertices are not snapped (R2). The
measurements below were taken before the second ruling, with snapped vertices;
see "What was measured" for what they still show.

**Closes.** `ROADMAP.md`'s catchment clip, pulled forward at the user's request
(2026-09-26), and the "stars" the user saw in 14b's output. After this
increment,

```sh
rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif \
    --domain quarter.geojson --tolerance 1 --out quarter.vtk
```

meshes only the inside of the polygon in `quarter.geojson`, with its vertices
where the file puts them. The start mesh is the constrained Delaunay
triangulation of the polygon's rings and nothing else. Refinement then inserts
DEM nodes to `--tolerance` as in 14b.

**Not closed.** Interior polygons and polylines as constraints (16b, R6).
Coarsening input geometry (`parallel_refinement.md` step 2, vector
simplification). Reprojection. A flat-ground size cap (increment 14, U2).
Cell-average DEM semantics (R0). Increment 14's slow carving from a NoData
corner (R5).

## The user's direction this design is built on

Given on 2026-09-26, in the user's words where quoted.

1. **"We don't specify loose nodes. We specify polygons and polylines. Then we
   have a tolerance for resolutions."** Input is geometry. No input is a set of
   free points. The stride grid of start vertices is such a set, so with a
   domain there is no start grid.
2. **The input is one confining polygon (it may have holes), any number of
   other polygons inside it (they may have holes), and any number of
   polylines.** Every ring and line becomes a constraint.
3. **"They are not discretized, they are already discrete. We will come to
   input geometry coarsening at a later stage."** Vertex lists are used as
   given. No arc sampling, no horizontal tolerance.
4. **The only tolerance in the tool is the vertical sup-norm `--tolerance`.**
5. **A DEM value is the exact height at a point, and height between points is
   bilinear** (R0).
6. **"Of course the input vertices will not align to the DEM, and they should
   not need to."** Input vertices stay where they are; their z is bilinear
   (R2).
7. The acceptance domain is a quarter circle, made outside the tool.

## What was measured

In the scratchpad (`dom16.py`), never in the tree. The real tile
`7908_3_10m_z33.tif`: 5051 × 5051 nodes, 10 m, EPSG:25833, nodes at
x = 799 750 … 850 250 and y = 7 899 750 … 7 950 250, so the quarter circle's
centre (850 250, 7 899 750) is the lower-right node exactly. Start meshes go
through the CLI's own `_engine` (`build_pslg` → `node` → `triangulate`, snap
spacing `DEFAULT_SNAP_SPACING`), then `_core.refine` and `elevation.trim`,
as `cli._dem_mesh` does. Apple clang build of 14b, `threads = 0`. "Degree" is
triangles per vertex after the trim. Min angle in world coordinates.

**These runs snapped the ring to DEM nodes, because 14b's `refine` accepts
nothing else.** Under R2 the ring stays off-node. The snapped ring is within
6.84 m of the true one and the interior is refined identically, so M1 and M3
are a fair guide to what R2 will give; they are not a measurement of it.
T-real re-measures after the change.

**M1. The quarter circle, boundary-only start, versus today's tile at stride 40.**
The ring: one vertex every 200 m on the arc (236 segments) and on each 30 km
border edge (150 each), 536 vertices.

| run | tol | start tris | triangles | deg ≥ 12 | deg ≥ 20 | max deg | median min ∠ | < 1° | worst ∠ | refine + CDT |
|---|---|---|---|---|---|---|---|---|---|---|
| today, stride 40 | 1 m | 32 258 | 463 974 | 324 (303 sea) | 196 | 74 | 45.0° | 0.00 % | 0.63° | 0.35 s |
| quarter circle | 1 m | 534 | 426 812 | 237 (217 sea) | 10 | 37 | 45.0° | 0.02 % | 0.31° | 1.00 s |
| today, stride 40 | 10 m | 32 258 | 68 569 | 202 (199 sea) | 196 | 74 | 45.0° | 0.00 % | 1.25° | 0.13 s |
| quarter circle | 10 m | 534 | 30 741 | 55 (44 sea) | 6 | 38 | 33.2° | 0.16 % | 0.24° | 0.61 s |

- The start grid's stars are gone: vertices of degree 20 or more fall from 196
  to 10 (1 m) and 6 (10 m). No boundary vertex has degree 12 or more.
- Triangles fall 8 % at 1 m and 55 % at 10 m. The land is all inside the
  domain, so the saving is sea the start grid used to pin.
- Every run met its tolerance (achieved 1.0000 m and 9.9975 m).
- **Grading at the coast is unbounded, as expected.** Median min angle at 10 m
  drops from 45° to 33°, and 0.16 % of triangles are under 1°: long sea
  triangles between the refined coast and the 200 m ring. The high-degree
  vertices left are inserted vertices, nearly all at z = 0 (217 of 237 at 1 m),
  not start vertices.
- `refine` is slower (1.0 s against 0.35 s at 1 m) because the first rounds
  scan 534 triangles as large as the domain.

**M2. Snapping dense input makes stars; one reason R2 does not snap.** The same
quarter circle given as 20 001 arc vertices (2.4 m apart) snaps to a 10 m
staircase of 5 950 ring vertices. At 10 m: max degree 124, 70 boundary vertices
of degree 12 or more, 2.56 % of triangles under 1°, worst 0.077°. Without
snapping, dense input stays a smooth dense ring; how many triangles it costs is
a question for the coarsening increment.

**M3. Coarse input does not hurt.** The quarter circle with only 67 ring
vertices (arc at ~1.1 km): at 1 m, 429 027 triangles, max degree 21; at 10 m,
31 030 triangles, max degree 19. Refinement splits boundary edges where the DEM
needs it. So the 200 m of the acceptance file is not load-bearing.

**M4. The whole tile, outline only (4 corners, 2 triangles).** 1 m:
431 674 triangles, but **5 054 rounds, 109 s**, max degree 1 770. 10 m: 5 052
rounds, 108 s, 14 % of triangles under 1°. The tile's top row is all NoData and
its corners are NoData vertices. See R5.

## Rulings

### R0. What a DEM value means (the user's reading 1)

- **A DEM value is the exact terrain height at one point**: the node itself for
  `PixelIsPoint`, the cell centre for `PixelIsArea`. `GTRasterTypeGeoKey`
  (1025) decides which. Increment 11 already reads it and moves an
  area-registered grid's origin by half a cell onto the centres
  (`io/geotiff.py`, `_placement`), and refuses a file without it. So "node"
  means the same point under both registrations, and nothing in this increment
  reads the key again.
- **Between nodes, height is the bilinear interpolation of the four surrounding
  values**: increment 12's `raster::bilinear` (`raster/sample.hpp`), including
  its refusal when any of the four is NoData or the point is outside the node
  rectangle.
- **The error measure is unchanged**: `|DEM value − triangle plane|` at each
  DEM node in the triangle's closed node set (14, R2). The tolerance is a
  property of the delivered mesh at every valid node inside the domain.
- **Not chosen: cell-average semantics.** A value as the mean height over its
  cell would call for a cell-average error measure (the triangle's integral over
  each cell against the value) and vertex z fitted to conserve volume rather
  than sampled. A possible later increment, if volume conservation turns out
  to matter (water routing, snow storage).
- **Unverified:** how Kartverket derives its 10 m values (point sample,
  mean or resampled from a finer model). Nobody has checked. Reading 1 is a
  ruling, not a finding about the data.

### R1. The domain input

- **`--domain PATH`**, only with `--dem`, and only with `--tolerance` (a
  boundary with no refinement is a few fan triangles; increment 12's bilinear
  path keeps its stride grid). Refused with `--stride`: a stride grid is loose
  points, which the user's direction excludes.
- **Formats: GeoJSON or WKT, by suffix** (`.geojson`/`.json`, `.wkt`). GeoJSON
  may be a bare geometry, a `Feature`, or a `FeatureCollection` with exactly one
  feature. Parsed with `json` and `shapely.geometry.shape`, or `shapely.wkt`.
  No new dependency; shapely is in `pyproject.toml`.
- **Accepted geometry: one `Polygon`, holes allowed.** Refused: `MultiPolygon`
  (and a collection of several), empty, and not `is_valid` (with shapely's
  `explain_validity` in the message). Orientation is set with
  `shapely.geometry.polygon.orient`: outer counter-clockwise, holes clockwise
  (increment 3's winding contract).
- **Vertices are used as given** (directions 3 and 6). No densifying,
  simplifying or snapping. The noder's own snap grid (`--snap-spacing`,
  default 1 mm) still applies, as to every input the engine has ever taken; see
  U6.
- **Extent: every vertex must lie in the DEM's node rectangle**, border
  included. Outside it there is no bilinear value (R0), so there is no z to
  give the vertex. Refused, naming the vertex and its distance outside. For an
  area-registered DEM this is half a cell inside the image edge. See U4 for
  clipping instead. The quarter circle's straight edges lie on the border
  row and column exactly, so they pass.
- **CRS, required, never transformed.** GeoJSON under RFC 7946 is WGS 84 unless
  it carries the 2008 spec's `crs` member, so: a GeoJSON file's CRS is its
  `crs` member (`urn:ogc:def:crs:EPSG::25833` or `EPSG:25833`), and a file
  without one is EPSG:4326 by the standard and refused as a mismatch. WKT has no
  CRS, so `.wkt` needs `--domain-crs EPSG:n`. The CRS must equal the DEM's EPSG
  code; no `Transformer` is created (increment 15, R6's rule).
- **`GeoPolygon` and `--bbox`.** Increment 15 says `GeoPolygon` lands with the
  clip, as its first caller. This is the clip, but it needs no transform, no
  `intersects` and no `buffer` yet. So this increment lands a smaller
  `DomainPolygon` (shapely polygon plus EPSG code, frozen) in a new
  `tin_engine/domain.py`. When 15 lands, its `plan_mosaic` takes
  `domain.polygon.bounds` in place of `--bbox`; `--domain` and `--bbox` are then
  mutually exclusive. `GeoPolygon` with `transform` grows from `DomainPolygon`
  when reprojection is needed.

### R2. Input vertices stay off the lattice; refinement's stay on it

Increment 14's R2 ("every vertex is a DEM node") becomes: **every vertex that
refinement inserts is a DEM node; a start vertex may be anywhere in the node
rectangle.** What that touches, and what it does not:

**The vertex type.** `LatticeMesh` stores `LatticeVertex{row, col}` in
`uint32`. It becomes `MeshVertex{double col, double row}`: fractional lattice
coordinates, with `is_node()` true when both are integers. For a node the
doubles are those integers exactly, so every node-only computation is the one
14b does today. An off-node vertex's `(col, row)` is computed once, in
`to_lattice`, as `(x − x_min) / dx` and `(y_max − y) / dy`. That rounding
defines where the vertex is for every predicate afterwards (it moves the point
by under a nanometre at UTM scale); its world `(x, y)` is kept alongside, as
given, for the output. `split_inside`, `split_edge` and the scan result keep
taking and returning a node (`LatticeVertex`), because only nodes are ever
inserted.

**Classification.** A start vertex whose world coordinates equal
`RasterGeometry::node` of its rounded `(row, col)` bit for bit is a node, as
today. Every other vertex inside the node rectangle is off-node. Outside the
rectangle `refine` refuses (`OutsideGrid`, replacing `OffLattice`).

**One frame for every predicate.** Orientation and incircle both run on
`(col, −row)` doubles (orientation) and `(col·dx, −row·dy)` (14b's
`LatticeFrame`, incircle), with `DefaultKernel`, which is exact on its inputs.
Node-only triangles keep today's `int64` path, which gives the same signs,
because their coordinates are small integers. A mixed triangle uses the kernel.
There is no second frame, so an inside test near a shared edge cannot disagree
between the two triangles that share it.

**The scan (14 R2) for a triangle with an off-node vertex.**
- The bounding box is `ceil`/`floor` of the fractional extents, clamped to the
  grid.
- Membership of node `p` is three exact kernel orientations in the frame.
  "On an edge" is still `Collinear`, not a tolerance.
- The plane at `p` is evaluated with `double` barycentric weights. It needs no
  exactness. Only the error value comes out of it, and its rounding is far
  below any tolerance.
- Nodes that are vertices are excluded by `is_node()` equality, as today.
  An off-node vertex is not a node, so it is never in a node set.

**"Zero error at vertices."** Still true, and still not needed. The plane
interpolates its three vertices' z, so it passes through them whatever those z
are. The error is measured only at nodes. What changes is **which height a
start vertex carries**:
- a node vertex gets `value_at`, as today;
- an off-node vertex gets `bilinear` (R0).

An input vertex that sits exactly on a node is classified as a node and gets
`value_at`, never `bilinear`. That matters because `bilinear` refuses a point
whose cell has a NoData corner, even when that corner's weight is zero.

The guarantee is unchanged in words, at every valid node inside the domain.
Between an off-node boundary vertex and its nearest nodes, the mesh interpolates
bilinear heights linearly, exactly as it already does between any two nodes.

**The Delaunay frame (14b R3).** `LatticeFrame::at` takes the fractional
coordinates. For nodes the result is bit-identical to today. 14b's termination
argument needs a fixed point set with exact signs. `at` is a deterministic
function of each vertex and the kernel is exact on its output, so the argument
holds unchanged. The Delaunay property is stated, as in 14b, for the points in
that frame.

**NoData near an off-node vertex.** When `bilinear` refuses, the vertex is
invalid, exactly like a NoData node:
- the triangles it belongs to are void and are carved (14 R6);
- the trim drops what is left and counts the vertex.

Carving's "valid node nearest a NoData vertex" uses a squared distance in the
fractional frame (`double`, same tie-break), not in `int64`. On the ground this
means the domain boundary is ragged where it runs within a cell of NoData.
That is the existing NoData ruling, not a new one.

**Termination (14b R4)** holds: every insertion is still a valid DEM node not
yet a vertex, from a finite lattice.

**Degenerate input.** An input vertex a nanometre from a node, or two input
vertices almost coincident, give slivers but no wrong sign. The kernel is exact
and `flip` asserts orientations through it. No special case; one fixture
each (T-deg).

### R3. The start mesh: the rings alone

- Chains: the outer ring as `ChainRole.Outer`, each hole as `ChainRole.Hole`,
  every edge with mask 0 (U5). Into `_engine`, then `refine`. No grid, no
  interior vertex.
- Measured (M1, M3), with a snapped ring: a start of 65 to 5 948 fan triangles
  with no interior vertex refines to tolerance, and 14b's once-only
  `legalise_all` runs on it unchanged.
- The domain need not be convex. The CDT drops outside triangles
  (increment 4), and a node set is the nodes in its triangle, so nodes outside
  the domain are never scanned.
- **The first rounds cost more than M1 shows.** Every start triangle has
  off-node vertices, so the whole domain is scanned once through the filtered
  kernel instead of `int64`: about 7 M nodes × 3 orientations for the quarter
  circle. Unmeasured; the filter's fast path should keep it in the tens of
  milliseconds. T-real reports it.

### R4. What the file records

- A new field **`domain`**: file name, rings (outer plus N holes), vertex
  count, e.g. `quarter.geojson, 1 ring 0 holes, 536 vertices`. Both `.vtk`
  (field) and `.ply` (comment), like `elevation_source`.
- **`elevation_source`** replaces `start stride N` with
  `start domain boundary, boundary z bilinear` on this path. The rest of the
  sentence is unchanged.
- The stderr report adds the start triangle count and how many start vertices
  were off-node.

### R5. NoData, sea, and a known defect in increment 14's carving

- NoData inside the polygon: 14's R6 unchanged, plus R2's off-node case.
- Sea is 0 m, not NoData, so it is meshed and refined only where needed (M1).
- **Known defect, not fixed here.** M4 took 109 s and 5 054 rounds, where
  the stride-grid start takes 0.35 s. The round count equals the tile's row count,
  and the tile's top row is NoData, so the likely cause is carving from a NoData
  corner advancing one node per round. It is a defect of increment 14's carving
  and needs its own look. It bites any domain with a vertex in a large NoData
  region. The quarter circle does not meet it.

### R6. The input model for all three kinds, so 16b is additive

| input | geometry | PSLG role | edge mask |
|---|---|---|---|
| confining polygon, exterior | `Polygon.exterior` | `Outer` | 0 (U5) |
| confining polygon, each hole | `Polygon.interiors` | `Hole` | 0 |
| interior polygon, exterior and each of its holes | rings | `Breakline`, closed (first point = last point) | the feature's property bit (a lake or land-cover class; the default vocabulary has none yet) |
| polyline | `LineString` | `Breakline`, open | the feature's property bit (`road`, `river` …) |

- **No new role for interior polygons.** A role says what a chain does to the
  domain: bound it, cut it, or neither. An interior polygon does neither; what
  it is, is a property (increment 7). Increment 3 already allows a closed
  breakline. Which triangles lie inside a lake is a region attribute derived
  later from the ring's constraint edges (ROADMAP's land-cover partitioning).
- **16b's CLI**, designed, not built: `--features PATH`, a GeoJSON
  `FeatureCollection` of `Polygon`/`MultiPolygon`/`LineString`/`MultiLineString`
  features, each with a `property` naming a vocabulary entry. A feature leaving
  the confining polygon is 16b's ruling, with increment 8's crossing gallery as
  its fixtures.
- **Crossings are now free.** A road over a river meets at a point the noder
  creates, generally not a DEM node. Under R2 it is just another off-node start
  vertex with bilinear z. Nothing in `refine` distinguishes it. This is the
  main reason R2 is better than snapping: snapping could not node crossings at
  all.

## Ruled by the user

- **2026-09-26: R0, reading 1** (point heights, bilinear between).
- **2026-09-26: former U2, option (b)**, no snapping (R2). The snapping
  machinery is dropped from the design, not kept as a fallback.
- **2026-09-26: U1 (a), U3 (a), U4 (a), U5 (a), U6 (a).**
- **2026-09-26, the user on U1: "I don't think the domain CRS should have to
  match the DEM CRS in the future. This must be written down. The questions
  will be in which coordinate system we shall do the math."** So U1 (a)'s
  must-match rule is this increment's scope, not a design principle. A later
  increment lets the domain, feature geometry and DEM each come in their own
  CRS, and must first rule on the **computation CRS**: the one coordinate
  system the noder, CDT, refinement and predicates work in, which today is the
  DEM's projected metric CRS by construction (`project_structure.md`: metres,
  no CRS in C++). Open questions it must answer: whether that stays the DEM's
  CRS (reproject vectors onto the raster grid, raster untouched) or becomes a
  chosen output CRS (which would mean resampling the DEM, and breaks reading
  1's "values are at nodes"); where reprojection happens (Python, pyproj,
  `always_xy`); how snapping to `--snap-spacing` and the DEM node frame relate
  across CRSs; and which CRS the output file records. Nothing in this increment
  may make the must-match rule load-bearing: the CRS check stays one
  replaceable function at the Python boundary.

The options below are kept so the reasons stay on record.

### U1. Where the domain's CRS comes from

- **(a) GeoJSON's `crs` member, or `--domain-crs` for WKT; must equal the DEM's.
  Recommended.** Honest about RFC 7946 (no member means WGS 84) and costs one
  flag.
- (b) Always assume the DEM's CRS. Zero flags, but a WGS 84 GeoJSON then fails
  the extent check, which names the wrong cause.
- (c) Reproject with pyproj now. That is `GeoPolygon.transform`; it waits for a
  caller that needs it.

### U3. The tile path without `--domain`

- **(a) Keep today's start grid on the legacy path, for now. Recommended.**
  `--tolerance` without `--domain` behaves exactly as in 14b, so the two can be
  compared in ParaView. `--stride` is legacy and gains no new behaviour. The
  switch waits for R5's defect.
- (b) Make the tile rectangle the domain, outline only, now. Consistent with
  "no loose nodes", but M4: 109 s instead of 0.35 s on the one real tile.
- (c) Later: the domain defaults to the DEM's valid-data footprint as a
  polygon. Avoids NoData corners; needs a raster-to-polygon step.

### U4. A polygon reaching outside the DEM's node rectangle

- **(a) Refuse (R1). Recommended.** There is no z outside it under R0, and a
  polygon past the tile's edge means the wrong DEM or the wrong CRS. Increment
  15 supplies more tiles.
- (b) Clip to the node rectangle with shapely `intersection`. What the legacy
  did. It can return a `MultiPolygon`, which then needs refusing anyway.

### U5. Does the domain boundary carry a property bit?

- **(a) No, mask 0. Recommended.** The `Outer` or `Hole` role already says it
  is the boundary.
- (b) A `domain_boundary` bit in the default vocabulary, so ParaView can colour
  it.

### U6. The noder's 1 mm snap on input vertices (new)

"Stay exactly where they are" meets the noder: `node` snap-rounds every vertex
to its `--snap-spacing` grid (default 1e-3, anchored at the CRS origin,
increment 5a). An input vertex therefore moves by at most 0.71 mm. Crossings in
16b need the noder, and its robustness rests on that grid.

- **(a) Accept it and say so. Recommended.** Sub-millimetre, recorded as the
  engine's input precision; DEM nodes are unaffected (checked in M1 to M4).
- (b) Bypass the noder when there are no crossings (this increment). Exact
  input positions, but two paths into the CDT, and 16b brings the noder back.

## Conflicts with increments 12, 14 and 14b

- **14 R1** (the start is a stride grid, and the stride is the size cap).
  Superseded on the `--domain` path; no size cap there (14 U2 is the way back).
- **14 R2** (every vertex is a DEM node). Relaxed to "every inserted vertex"
  (R2). The scan's `int64` membership stays for node-only triangles.
- **14 R6** (carving by `int64` lattice distance). Generalised to the
  fractional frame for off-node void vertices. Plus R5's defect.
- **14b R3** (`LatticeFrame` on integers). Takes fractional coordinates;
  bit-identical for nodes.
- **14b R6 and C1** (the start stride). Superseded on the `--domain` path; C1
  survives on the legacy tile path (U3 a).
- **12 R1** (stride subsample without `--tolerance`). Untouched.
- **`refine`'s `OffLattice` refusal** is replaced by `OutsideGrid`. Increment
  14's tests that expect `OffLattice` for an off-node vertex are amended in the
  red commit, with the reason in the message.

## Prior art in `legacy/`

```sh
$ grep -rliE 'polyfile|GeoPolygon|point_inside_polygon|from_polygon_file' legacy/
legacy/rasputin/globcov_repository.py
legacy/rasputin/land_cover_repository.py
legacy/rasputin/wfs_repository.py
legacy/rasputin/triangulate_dem.h
legacy/rasputin/geo_tiff_reader.py
legacy/rasputin/application.py
legacy/rasputin/mesh.py
legacy/rasputin/reader.py
legacy/rasputin/geometry.py
legacy/rasputin/gml_repository.py
legacy/tests/test_gml_repository.py
legacy/tests/test_land_cover_repository.py
legacy/tests/test_mesh.py
legacy/tests/test_raster_repository.py
```

What matters, read directly:

- `legacy/rasputin/application.py:37-39, 89-99`: the domain was `-polyfile`
  (WKT or WKB) or `-x`/`-y` coordinate lists, with the CRS given separately,
  and it was required. Carried across: a file, and a stated CRS (U1).
- `legacy/rasputin/geometry.py:205-259`: `GeoPolygon`, shapely plus pyproj.
  `from_polygon_file` picks the format by suffix, as R1 does. Its `transform`
  waits (U1 c); increment 11 ruling 10 keeps the `not touches` rule for 15.
- `legacy/rasputin/triangulate_dem.h:428-470`: the polygon was intersected with
  the raster extent (U4 b), and each boundary edge not on the raster border was
  densified at raster resolution with **bilinear z at off-node points**, the
  same z rule as R2 (but R2 does not densify). `:474-508`: every DEM point
  inserted, then triangles kept by centroid in the polygon. Not carried: the
  CDT drops outside triangles itself, and there is no point cloud to filter.

`@migration-expert` is not needed: the intent is stated above and nothing
numeric is ported.

## Files

| file | what |
|---|---|
| `include/terrain/mesh/lattice_mesh.hpp` | `MeshVertex` (fractional `col`, `row`, `is_node()`); `build` and `flip` orientation checks through the kernel for mixed triangles, `int64` otherwise |
| `include/terrain/mesh/lawson.hpp` | `LatticeFrame::at` on `MeshVertex` |
| `include/terrain/refinement/scan.hpp` | the mixed-triangle path: fractional box, kernel membership, `double` plane; carving distance in the fractional frame |
| `include/terrain/refinement/refine.hpp` | `to_lattice` classifies node / off-node / outside; `OutsideGrid`; off-node z by `bilinear`, invalid when it refuses; world coordinates of off-node vertices carried to the output as given |
| `bindings/core.cpp`, `_core.pyi` | the status rename |
| `src_python/tin_engine/domain.py` (new) | `DomainPolygon`; `read_domain(path, crs)`: formats, CRS, geometry and extent refusals, orientation, chains. Pure: shapely, numpy, `RasterMeta`; no `_core`, no typer |
| `src_python/tin_engine/cli.py` | `--domain`, `--domain-crs`; R1's refusals; the domain branch in `_dem_mesh`; the `domain` field and the sentence (R4) |
| `project_structure.md` | `domain.py` in the tree |
| `ROADMAP.md` | row 16 |

## Tests for `@tester`

**Invariant-critical (mutation testing required):** two suites.

- **`test_refinement_scan`, extended with mixed triangles.** It decides the
  tolerance guarantee near the boundary.
  - S1: a node exactly on an edge between two off-node vertices is in both
    adjacent triangles' sets, and is `Edge*`, not `Inside`.
  - S2: a node a nanometre outside such an edge is in exactly one.
  - S3: plane error agrees with a brute-force oracle using **the producer's
    relation** (exact kernel on the fractional frame), not world coordinates.
  - S4: node-only triangles give bit-identical results to 14b's.
  - Mutants: the `int64` path used for a mixed triangle (truncated coordinates);
    the box computed with `floor` on both ends; an off-node vertex compared to
    nodes by rounded `(row, col)`, which would exclude a real node.
- **`prop_refinement_refine`, the tolerance oracle (14 T3), with off-node start
  rings.** Random polygons with off-node vertices over synthetic terrain; every
  valid node inside the domain is within tolerance; the output is constrained
  Delaunay in the frame (14b R10). Mutant: off-node z from `value_at` of the
  rounded node instead of `bilinear`.

Everything else is property or integration testing.

- **Z1. Vertex z.** An off-node vertex's z equals `bilinear` at its world point;
  a vertex exactly on a node gets `value_at`, including next to a NoData node
  (where `bilinear` would refuse).
- **Z2. NoData.** An off-node vertex in a cell with a NoData corner is invalid,
  its triangles are carved, and the trim counts it.
- **Z3. Output positions.** Off-node vertices come out at their input world
  coordinates, bit for bit, not recomputed from the frame.
- **T-deg.** A vertex 1e-9 cell from a node; two vertices 1e-9 apart. Refines,
  terminates, Delaunay holds.
- **R1 refusals** through the CLI: `MultiPolygon`, invalid, empty, a GeoJSON
  without `crs`, a mismatched EPSG, `.wkt` without `--domain-crs`, a vertex
  outside the node rectangle, `--domain` with `--stride`, `--domain` without
  `--tolerance`.
- **Synthetic end-to-end:** the T12 cone-and-island raster with a square domain
  holding one hole, all vertices off-node. No output vertex inside the hole or
  outside the square, and the tolerance holds.
- **Real-tile acceptance (T-real), logged, not thresholded:** the quarter
  circle at 1 m and 10 m. Report triangles, rounds, the first round's scan
  time, the degree distribution (median, p99, max, count ≥ 12 and ≥ 20) and min
  angles (median, < 1°, < 10°, worst). Assert the tolerance, that the 536 input
  vertices appear in the output at their input coordinates, and that the land
  is present. M1 is the reference. The input is built by the test from the
  formula, not committed: one vertex every 200 m on the arc of radius 30 km
  about (850 250, 7 899 750) and on both border edges, with the `crs` member for
  EPSG:25833.

## LOC estimate

Counted in `CLAUDE.md` §2's unit.

| file | what | est. |
|---|---|---|
| `mesh/lattice_mesh.hpp` | `MeshVertex`, mixed orientation in `build` and `flip` | ~35 |
| `mesh/lawson.hpp` | frame on `MeshVertex` | ~5 |
| `refinement/scan.hpp` | mixed path, carving distance | ~50 |
| `refinement/refine.hpp` | classification, `OutsideGrid`, bilinear z, world output | ~35 |
| `bindings/core.cpp`, `_core.pyi` | status rename | ~5 |
| `domain.py` | model, reader (2 formats, CRS), refusals, orientation, chains | ~85 |
| `cli.py` | two options, refusals, domain branch, field, sentence | ~45 |
| | **total** | **~260** |

On the worst overrun seen so far (+39 %), about 360. **It fits under 700, so it
is not split.** The one natural seam, if a review round grows it, is C++ first
(off-node start vertices through `refine` and the binding, testable from
Python with a hand-built ring) and the domain reader and CLI second. Only the
second lets the user mesh the quarter circle, so the split would not deliver
the acceptance run in its first half, which is why it is not recommended. 16b
(features, closed and open breaklines) is about 100 more, with no C++ change
now that crossings are off-node vertices.

## Acceptance

- The quarter circle through the CLI writes a `.vtk` ParaView opens, with no
  triangle outside the quarter circle, the input vertices where the file put
  them, the `domain` field, and `elevation_source` saying
  `start domain boundary, boundary z bilinear` and the achieved error.
- Side by side with today's `--tolerance 1` output, the sea stars are gone
  (M1: degree ≥ 20 from 196 to about 10).
- All gates in `CLAUDE.md` §4 green, including the TSan job; CI green.

## Not in scope

- Interior polygons and polylines (16b, R6).
- Coarsening or simplifying input geometry (`parallel_refinement.md` step 2).
- Cell-average DEM semantics (R0).
- Reprojection (U1 c), multi-part domains, clipping (U4 b).
- A size cap for flat ground (14 U2), the tile path's switch to outline only
  (U3), and R5's carving defect.
