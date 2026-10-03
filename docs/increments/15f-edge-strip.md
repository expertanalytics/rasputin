# Increment 15f: the edge strip — check points where constraints cross grid lines

Status: **designed by `@architect`, 2026-10-03, while Ola was away
(unattended).** Written before `@tester`, per `docs/increments/README.md`
step 1. Nothing here is implemented. Choices that would normally go to Ola
were made as defaults; each is marked *default* where it occurs and listed
under "Defaults chosen" at the end.

Item 5 of the order Ola ruled for the basin (`23-basin-scale.md`, "The order
(B12, ruled (a))"). The one-paragraph design it implements is in
`15c-geographic-dem.md`, "The edge strip (Surprise 3), placed", and Q14 there.

**Why "15f".** The edge strip was placed "right after 15c" (Q14), so it takes
a 15 letter. `15d` is taken: it names the window-decoding design in
`15-dem-mosaic.md` (R3 [15d], R7 [15d]), which increment 23a-1 replaced, and
`ROADMAP.md` and `23-basin-scale.md` still say "23a-1 replaces 15d". `15e` is
the memory fixes (branch `worktree-agent-ab601ec06c6508cb9`). So the first free
id is `15f` (*default*). `15e-memory-fixes.md` line 17 on that branch says
"15d is taken by the edge strip". That stale line is noted here for 15e's PR,
which should change it to "15f" before its push. It is not edited from this
branch.

## What is ruled

- **Q14 (a), extended (Ola, 2026-10-01).** Check points at every crossing of a
  constraint edge with a grid line, with z linear between the two nodes of that
  cell side. **Also the midpoint between each two neighbouring crossings**,
  with z bilinear from the one cell that holds it. They are run through 15c's
  `refine_points` after `refine`, for any DEM. The guarantee wording was
  proposed and Ola did not reword it: "at every DEM node, wherever a constraint
  crosses a grid line, and at the midpoint between neighbouring crossings".
- **Proposed in 15c and left for Ola to overrule:** an edge's two end vertices
  count as neighbours, so the pieces at either end of an edge get a midpoint
  as well. *Adopted as the default here.*
- **B12 (a) (Ola, 2026-10-01).** The edge strip comes after 15c-2 and before
  23b.
- **What 23 added** (`23-basin-scale.md`, "What survives of 15c, the edge
  strip and 15d"). The check-point generator is one C++ function,
  `constraint_check_points`, which 23b's seam pass reuses. The strip makes no
  check points on frozen edges. 23b's estimate assumes the generator exists
  (`23-basin-scale.md`, the paragraph after the PR table).

## Three findings that change the one-paragraph design

These came from reading the code the strip builds on. Each one adds lines
beyond 15c's estimate of 160-210.

**F1. Membership cannot find a point on a constraint.** `scan_points`
(`include/terrain/refinement/refine_points.hpp:83`) decides which triangle
holds a check point with three exact orientation tests on the point's stored
position. A crossing computed in floating point lies within rounding of its
edge, on one side or the other. On the domain's outline, a point on the
outer side is in no triangle, so it is never scanned and is silently left
unchecked. On an interior constraint, the point is inserted with
`split_inside`, which puts a vertex a hair away from the constraint and
leaves a near-zero-area triangle that no flip can remove. 15c foresaw this
(its reason 3: "they must be filed against their edge, not found by
membership"). So in this design, strip points are filed **by edge and by
parameter along it**, and scanned per constrained sub-edge, not by triangle
membership (D3, D4).

**F2. On the projected path, an insertion on a constraint can break refine's
guarantee at DEM nodes.** A strip insertion splits a constraint edge and
changes the planes of the two triangles beside it, followed by Lawson flips.
A DEM node in those triangles that was within tolerance of the old plane need
not be within tolerance of the new one. On the reprojected path the grid's
nodes do not matter, because 15c's J2 already gives up the guarantee there. The
source nodes do matter, and they are covered because phase 2 and the strip
share one loop (D1). On the projected path (Norway, any DEM meshed directly) it would turn
"every DEM node within tolerance" false. So the strip run there also
**rescans every triangle it writes against the DEM's nodes**, with refine's own
`scan` (`include/terrain/refinement/scan.hpp:90`), and inserts nodes as refine
would (D4, step 3). Run "after refine" as written, the strip would trade one
guarantee for another.

**F3. Increment 15e drops the target tile before phase 2.** On branch
`worktree-agent-ab601ec06c6508cb9`, `_dem_mesh` runs `del tile` after
`refine` and before `final_check.run`, so the resampled grid is gone when phase
2 runs. Strip points need z from that grid. So they are **generated
right after `refine`, before the tile is dropped**. They are a small, separate
set (D3), and only they are carried into phase 2.

**A remark that is not a finding (F4).** Along a straight segment inside one
cell, the bilinear surface is a quadratic in the parameter along the segment
(the `xy` term makes it so). The mesh along a constraint is linear between
vertices. So the error along one piece between neighbouring crossings, with no
vertex inside it, is a quadratic. Its three checked values (the two crossings
and the midpoint) bound it to **1.25 × tolerance at every point of that piece**:
1.25 is the Lebesgue constant of the three equally spaced nodes {0, ½, 1}
(|L₀| + |L½| + |L₁| peaks at t = ¼ and ¾ with 3/8 + 3/4 + 1/8). The bound
fails on a piece where an inserted midpoint now sits inside, so it is a remark
and not a guarantee. It explains why the midpoint is the right single extra
point: a quadratic departs furthest from its chord there. Q1 asks whether Ola
wants the every-point version.

## Prior art: legacy and literature

### Literature

- **Greedy insertion, unchanged.** Garland and Heckbert, "Fast polygonal
  approximation of terrains and height fields", CMU-CS-95-181, 1995. Over
  scattered samples, it goes back to De Floriani, Falcidieno and Pienovi,
  *CVGIP* 32:127-140, 1985. Both are cited as in 14 and 15c. The strip is that
  method with one more sample set: points on the constraints. **What differs:**
  these samples are filed against the constraint they lie on and inserted on
  it, so the method's guarantee (every sample within tolerance) holds for
  samples that triangle membership would miss (F1).
- **Enumerating where a segment crosses grid lines.** Amanatides and Woo, "A
  fast voxel traversal algorithm for ray tracing", *Eurographics '87*, 1987
  (recalled). For height fields, Musgrave, "Grid tracing: fast ray tracing for
  height fields", Yale research report YALEU/DCS/RR-639, 1988 (recalled).
  These step a ray from one cell boundary to the next. Here the crossings of
  one edge are enumerated directly, by the integers between its ends, in each
  axis. There is no traversal state, so each edge is independent and the
  result does not depend on direction (D2).
- **Draping lines on a surface.** ArcGIS's Interpolate Shape (3D Analyst)
  gives a line z from a surface. On a raster it samples at a fixed distance
  that defaults to the cell size. On a TIN it densifies the line at the
  triangle edges it crosses, its "natural densification"
  ([tool reference](https://pro.arcgis.com/en/pro-app/3.4/tool-reference/3d-analyst/interpolate-shape.htm),
  read 2026-10-03). **What differs:** the grid-line crossings are the raster's
  natural densification. On a bilinear surface they are where the line's
  height profile changes from one quadratic piece to the next (F4). Uniform
  sampling at the cell size misses those breaks. Also, nothing is inserted
  unless the tolerance needs it. The points are checks, not vertices.
- **Breaklines in tolerance-driven TIN simplification.** ArcGIS's Decimate TIN
  Nodes guarantees a z tolerance at the source TIN's data nodes. It copies
  breaklines "without any generalization", outside the node budget
  ([how it works](https://help.arcgis.com/en/arcgisdesktop/10.0/help/00q9/00q900000053000000.htm),
  read 2026-10-03). So that approach keeps the guarantee along breaklines by
  never simplifying them. **What differs:** here constraints are simplified to
  what the tolerance needs, and their heights are checked at the grid's own
  breakpoints.

**Novelty: none claimed.** The strip is greedy insertion with one more set of
samples, and the samples are the textbook crossings. Searched, 2026-10-03:
"TIN simplification breaklines vertical error bound along constraint edges DEM
greedy insertion constrained Delaunay check points on breaklines" and "terrain
profile extraction DEM line intersects grid lines bilinear interpolation cell
boundary crossings". Found: GIS tool documentation (above), SAGA's profile
tools and the GRASS wiki on TINs with breaklines. None of it states a vertical
guarantee at the crossings of constraints with grid lines. F4's 1.25 bound is
elementary interpolation theory and is not claimed as new either. If Ola later
wants to publish the every-point form (Q1), the search should cover
Scholar and IEEE Xplore for "breakline" with "error bound" and "TIN" before
any claim is made.

### Legacy

```sh
$ git grep -l -i -E "breakline|densif|crossing|interpolate_shape|along the edge|edge_points" legacy-archive -- legacy
$ git grep -n -i -E "insert_constraint|Constrained_Delaunay|refine|bilinear" legacy-archive -- legacy
legacy-archive:legacy/rasputin/triangulate_dem.h:11:#include <CGAL/Constrained_Delaunay_triangulation_2.h>
legacy-archive:legacy/rasputin/triangulate_dem.h:49:using ConstrainedDelaunay = Constrained_Delaunay_triangulation_2<Gt, CGAL::Default, CGAL::Exact_predicates_tag>;
legacy-archive:legacy/rasputin/triangulate_dem.h:375:    // Interpolate data using using a bilinear interpolation rule on each cell
legacy-archive:legacy/rasputin/triangulate_dem.h:388:        // Using bilinear interpolation on the celll
legacy-archive:legacy/rasputin/triangulate_dem.h:485:        dtin.insert_constraint(point_sequence.begin(), point_sequence.end(), false);
```

The first grep returns nothing. The second leads to the one legacy routine
on this subject: `interpolate_boundary_points`
(`legacy/rasputin/triangulate_dem.h:427-471` at the `legacy-archive` tag),
read in full. **Its intent:** the domain boundary, clipped to the raster
rectangle, is cut into `max(|Δx|/dx, |Δy|/dy)` equal pieces (at least one; line
452-454). Every piece end gets bilinear z (`get_interpolated_value_at_point`,
line 375-392). All of them are inserted as constraint vertices (line 485).
Edges along the raster's rectangle are skipped (line 443). So the legacy mesh
carried the bilinear height about once per cell along the domain boundary,
unconditionally, with no tolerance. It did not do this for interior
constraints.

**What is carried:** the intent, that a constraint's heights follow the
bilinear surface at the grid's resolution. The strip delivers it with
checks rather than unconditional vertices, and on every constraint, not only
the boundary. **What is not carried:** the uniform subdivision. Its sample
points fall between grid lines, so they miss the breaks of the height profile
(above), and inserting all of them ignores the tolerance. No domain constant
is re-derived: the legacy's only constant is the subdivision count, which is
dropped.

## Scope

In:

- the generator, `constraint_check_points`, in C++, over an explicit list of
  edges (D2);
- its store, `ConstraintCheckPoints`, grouped by edge and sorted along it
  (D3);
- the strip in the refinement loop: `refine_points` takes the strip
  (reprojected path), and a new entry point, `refine_strip`, runs it on the
  projected path with the DEM rescan of F2 (D4);
- the binding, the stubs, a small Python module `edge_strip.py`, and the CLI
  wiring on both paths of `rasputin mesh --tolerance` (D5, D6);
- what the `.vtk` and `--stats` record (D7).

Not in scope:

- **Frozen edges and the seam pass** (23b). The generator takes the edges
  it is given; 23b passes the non-frozen constraint edges to the strip and the
  seam edges to its seam pass. No frozen mask exists before 23b.
- **The every-point form** (F4, Q1): midpoints added recursively when one is
  inserted. Not ruled.
- **Constraint feet and start quality in the strip run.** As in 15c's phase 2
  (its D5). A DEM node the F2 rescan inserts next to a constraint goes in with
  `split_inside`, as phase 2's source nodes do.
- **One loop shared by `refine` and `refine_points`** (15c D5's later
  refactor). This design refactors only `refine_points`' own loop (D4).
  `refine.hpp` gets one extracted helper and no change in behaviour.
- **A CLI switch to turn the strip off** (*default*: none; see "Defaults
  chosen").
- **`rasputin mesh` without `--tolerance`** (sampling only, increment 12):
  there is no guarantee to extend, so there is no strip.

## The design

### D1. Data flow

```
                      projected path (Norway)            reprojected path (15c-2)
start mesh ──► refine(dem) ──► out ──┐                   refine(target grid) ──► out ──┐
                                     ▼                                                  ▼
          constraint_check_points(dem, out.vertices,     constraint_check_points(target grid,
                                  out.edges) = strip         out.vertices, out.edges) = strip
                                     │                   del tile  (15e, fix 3)         │
                                     ▼                                                  ▼
          refine_strip(dem, strip, out ...)              refine_points(source store, out ...,
            strip points + DEM rescan of written           strip): source nodes + strip points,
            triangles, one loop                            one loop
                                     │                                                  │
                                     ▼                                                  ▼
                                   trim                                               trim
```

On both paths one loop handles every point set together, so no set is checked
against a mesh that a later set then changes. The strip is generated from
`refine`'s output constraint edges, the mesh's own sub-edges after phase 1, and
not from the input polygons. Those edges are what the surface follows, and
every one of them has two vertices of the start of the strip run.

Boundaries: Python never builds a point; the generator, the store and the
loop are C++ and see numbers in metres and lattice units only (I/O boundary,
`CLAUDE.md` §2). The C++ never sees a CRS or a path. `edge_strip.py` composes
calls and records times; it holds no geometry.

### D2. The generator: `constraint_check_points`

**Site:** `include/terrain/refinement/constraint_points.hpp` (new).

```cpp
template <raster::RasterSource R>
[[nodiscard]] ConstraintCheckPoints constraint_check_points(
    const R& dem, std::span<const Point2> vertices,
    std::span<const std::array<std::uint32_t, 2>> edges);
```

`dem` is the grid refine read: the DEM itself on the projected path, the
resampled target grid on the reprojected one. `vertices` are world points,
those of the mesh the strip run will start from; `edges` are vertex-index
pairs, normally that mesh's constraint edges (`RefineOutcome::edges`). An end
outside the node rectangle, or an index out of range, is
`std::invalid_argument` (refine already refused such a mesh, so in a run this
is a programming error).

Per edge, in this order:

1. **Lattice ends.** Each end goes to lattice coordinates `(col, row)` by
   **the same function `refine` uses**, extracted from `detail::to_lattice`
   (in `include/terrain/refinement/refine.hpp`) as
   `detail::lattice_position(g, p) -> mesh::MeshVertex`: a node exactly when
   `g.node` of the rounded position gives `p` back bit for bit, otherwise the
   clamped fractional position. `to_lattice` then calls it. The extraction
   changes no behaviour, and it means the generator's ends are, bit for bit,
   the vertices the loop builds.
2. **Canonical direction.** The edge runs from its lower vertex index `P0` to
   its higher `P1`, whatever order the pair was given in. So the points,
   computed in floating point, do not depend on the pair's order (test CC6).
3. **Crossings with column lines.** For every integer `K` with
   `min(c0, c1) < K < max(c0, c1)`: `t = (K − c0) / (c1 − c0)`, the point
   `(K, r0 + t (r1 − r0))`. The column is exactly `K`, so the point is
   on the grid line exactly and is off its edge only by rounding.
   **Crossings with row lines**, likewise for every integer `R` strictly
   between `r0` and `r1`: `(c0 + t (c1 − c0), R)`. A family is empty when the
   edge is parallel to it (`c0 == c1` or `r0 == r1`), so nothing divides by
   zero. Each computed coordinate is clamped to the node rectangle.
4. **Crossings at nodes.** A column crossing whose row rounds to an integer
   `R` with `orient_sign(P0, P1, (K, R)) == 0` (exact) is replaced by the node
   `(K, R)` itself. A row crossing is handled the same way. The two crossings
   of one node then have one position, and step 5 keeps one of them.
5. **Order and duplicates.** Crossings are sorted by `t`, with ties broken by
   position. They are then walked into a list that starts with `P0` at
   `t = 0`. A crossing is dropped, and counted in `duplicates`, when:
   - its `t` is not strictly between the last kept entry's `t` and 1
     (`t ≤ t_last` or `t ≥ 1`); or
   - its position equals the last kept entry's position.

   `P1` is then appended at `t = 1`. The bound `t ≥ 1` is needed because `t`
   can round to exactly 1.0 when `P1` lies within an ulp past a grid line
   (@reviewer's probe, round 1: ends at col 0.8932792255671602 and
   `nextafter(10, 11)`). A crossing's position can never equal an end's: a
   column crossing has `col == K` exactly, and `K` is strictly between the
   ends' columns, even after clamping and node snapping; a row crossing
   likewise has `row == R`. So comparing the first crossing with `P0`, or
   the last with `P1`, by position can never fire. The implementation may
   compare against `P0` or skip it; either way the result is the same. The
   `t` bounds, though, must be checked at both ends.
6. **Midpoints.** Each consecutive pair `(l, r)` of the list, ends included
   (ends count as neighbours, *default*), gives one candidate: the position
   `((l.col + r.col)/2, (l.row + r.row)/2)` with `s = (t_l + t_r)/2`. The
   candidate is dropped, and counted in `duplicates`, when its position equals
   `l`'s or `r`'s, or its `s` equals `t_l` or `t_r`. All of these happen only
   on a piece a few ulps long (@reviewer's probe, round 1:
   `P0.col = nextafter(1, 0)` and a crossing at col 1 give a midpoint at the
   crossing's position, `(1, 1.5)`). The ends themselves are not check
   points: a vertex's z is its own.

**The invariant of `on_edge(k)`**, which steps 5 and 6 make true and which
`@tester` pins:

- **I1.** Every point has `0 < s < 1`.
- **I2.** `s` increases strictly along `on_edge(k)`.
- **I3.** In the sequence `P0, points…, P1` as steps 5 and 6 leave it,
  before step 7's NoData drops, no two neighbours share a position. So
  without NoData, no point is at either end's position, and no two
  consecutive points of `on_edge(k)` share one. Nothing is claimed for
  non-neighbours, or for points that become neighbours only because a NoData
  point between them was dropped. Such a coincidence can only arise at ulp
  scale. The loop (15f-2) refuses that insertion through `foot_fits`, because
  a child triangle would have zero area, and counts it.

I1 is what the loop needs: it scans the open interval `(s_a, s_b)`, so a
point with `s = 0` or `s = 1` would never be checked, silently. I3 keeps it
from inserting a vertex on top of an end.

**A degenerate edge.** A pair with equal indices (`e[0] == e[1]`) is
refused with `std::invalid_argument`, in a message containing
`constraint_check_points`. No mesh has such an edge, so it can only be a
programming error. Two distinct indices at one position give no points:
there is no crossing, and the one midpoint candidate sits on both ends, so
step 6 drops it.
7. **Heights.** Every point's z is `vertex_z(dem, point)`
   (`include/terrain/refinement/scan.hpp:67`), the function refine uses for an
   off-node vertex. At a node it is `value_at`. For a crossing, the cell's
   fraction across the line is exactly 0, so the expression is exactly
   `z₀ (1 − f) + z₁ f` between the two nodes of the cell side, which is Q14's
   "z linear between the two nodes". For a midpoint it is bilinear in the
   cell that holds it. A point for which `vertex_z` refuses (a corner of its
   cell is NoData, whatever its weight, refine's own rule) is not kept, and
   is counted in `no_data`.

The parameter `s` stored with each point is its `t`. Per edge that gives
`2c + 1` points for `c` kept crossings, before NoData and the ulp-scale drops
of step 6. An edge with no crossing gets its one midpoint.

### D3. The store: `ConstraintCheckPoints`

**Site:** `constraint_points.hpp`, beside the generator.

```cpp
struct ConstraintPoint {
    mesh::MeshVertex at;  // lattice (col, row), as the generator computed it
    double z;             // vertex_z there
    double s;             // the parameter along its edge, P0 = 0 to P1 = 1
};
static_assert(sizeof(ConstraintPoint) == 32);

class ConstraintCheckPoints {
public:
    const raster::RasterGeometry& geometry() const noexcept;
    std::size_t edge_count() const noexcept;
    std::array<std::uint32_t, 2> edge(std::size_t k) const;      // (P0, P1), lower index first
    std::span<const ConstraintPoint> on_edge(std::size_t k) const;  // s strictly increasing
    std::size_t size() const noexcept;        // points kept
    std::size_t no_data() const noexcept;     // dropped: NoData in their cell
    std::size_t duplicates() const noexcept;  // dropped: same position or s
private:
    raster::RasterGeometry geometry_;
    std::vector<std::array<std::uint32_t, 2>> edges_;
    std::vector<std::size_t> offsets_;   // edges_.size() + 1, into points_
    std::vector<ConstraintPoint> points_;
};
```

It is built only by the generator and is immutable afterwards, so any number of
threads may read it. `edge(k)` and `on_edge(k)` check `k` (`.at()`) and throw
`std::out_of_range` past `edge_count()`, as built in 15f-1 and kept as the
design. The loop's hot path reads `on_edge` once per sub-edge, so the check
costs nothing measurable. **Why not 15c's `CheckPoints`:** that store keeps a
`float` offset in the cell and finds points by membership (F1). The strip
needs the edge, the parameter, and z in `double`. A `float` z would round a
bilinear value near 2,000 m by up to about 6e-5 m (half of the `float`
spacing there, 2⁻¹³ m), which is visible against the oracle's 1e-9.

**Memory:** 32 bytes per point and 16 per edge in the store, plus about 64
bytes per constrained sub-edge in the loop's map (D4, step 1) while the run
lasts. Estimated for Bygdin with
CORINE at 10 m by `@architect` (shapely over
`docs/benchmarks/2026-09-29/bygdin-landcover/bygdin_reduced_t20.geojson` and
`../rasputin_data/corine2018_dtm10_utm33.gpkg`, summing `|Δx| + |Δy|` over the
cell size): about 17,000 crossings on the outline (128.8 km) and 66,000 on the
CORINE borders inside it (about 577 km). That gives about 170,000 points with
the midpoints, about 5.4 MB. It is small next to the DEM on every path
measured so far. 15e's arena (its fix 4) does not apply: it is a layout for
15c's store, and this store is one contiguous vector sized by the generator.

### D4. The loop

**Site:** `include/terrain/refinement/refine_points.hpp`.

`refine_points`' body (lines 169-292) moves into `detail::point_loop`,
parametrised by what it scans. 15c's entry point calls it with the store
alone and behaves exactly as today: 15c's RP1-RP7 suite is the regression
test. Two callers are added:

```cpp
// The reprojected path: source nodes and the strip, one loop.
template <class Store>
[[nodiscard]] PointRefineOutcome refine_points(
    const Store& points, const IndexedMesh2& start, std::span<const double> z,
    std::span<const std::uint8_t> valid, std::span<const std::array<std::uint32_t, 2>> edges,
    std::span<const std::uint32_t> masks, const PointRefineOptions& options,
    const ConstraintCheckPoints* strip = nullptr);

// The projected path: the strip, and the DEM's nodes in every written triangle.
template <raster::RasterSource R>
[[nodiscard]] PointRefineOutcome refine_strip(
    const R& dem, const ConstraintCheckPoints& strip, const IndexedMesh2& start,
    std::span<const double> z, std::span<const std::uint8_t> valid,
    std::span<const std::array<std::uint32_t, 2>> edges, std::span<const std::uint32_t> masks,
    const PointRefineOptions& options);
```

**Refusals.** A strip whose `geometry()` differs from the store's (or the
DEM's) is `std::logic_error`: a programming error, as 15c's unfrozen store is.
So is a strip edge that is not a constraint edge of `start`, because a strip
built on another mesh would otherwise check nothing, silently. A start
constraint edge with no strip edge is allowed: it has no strip points. That is
how 23b leaves frozen edges out.

**Step 1, sub-edges.** After `to_lattice` and `legalise_all` (unchanged), the
loop builds a map from each constrained sub-edge (its two vertex indices,
unordered) to a record `{k, a, b, s_a, s_b}`: the strip edge it lies on, and
the parameters of its two ends. At the start each strip edge is one sub-edge
with `s` 0 and 1. The parallel scan only reads the map. The serial split phase
updates it (step 4).

**Step 2, the scan** (parallel, read-only, one result per triangle):

- **Source points** (reprojected path): `scan_points` exactly as today.
- **Strip points.** A triangle scans each constrained edge it **owns**: the
  edge goes from the lower vertex index to the higher in this triangle, or no
  triangle lies across it. Each sub-edge therefore has exactly one owner, and
  a boundary edge always has one. For a sub-edge `{k, a, b, s_a, s_b}`, the
  points are those of `on_edge(k)` with `s` strictly between `s_a` and
  `s_b`, found by binary search. A point is skipped if it is marked refused
  (step 3). For each, `σ = (s − s_a)/(s_b − s_a)`, the mesh value is
  `z_a + σ (z_b − z_a)`, with `z` from the loop's vertex table, and the error
  is `|z − mesh value|`. This is the mesh surface on the sub-edge, which both
  triangles beside it share.
  - **A void sub-edge** (an end with no z): no error is computed. The point
    nearest an invalid end (smallest `|s − s_end|`) is named for carving, and
    the sub-edge's points are counted in `uncovered`. This is 14's rule,
    applied along the edge.
- **DEM nodes** (`refine_strip` only): `scan(dem, m, t)`, refine's scan,
  unchanged, in every triangle the loop has written. In round 1 that means the
  triangles `legalise_all` flipped (none on refine's output, normally). After
  round 1 the active set is exactly the written triangles, so every active
  triangle is scanned.
- **Combining.** The first void result in the order source, strip, DEM wins.
  Otherwise the strictly largest error wins, with ties to the earlier set in
  that order, and within the strip to the lower triangle-edge index, then the
  lower `s`. None of this depends on the thread count or on the order of the
  edge list.

**Step 3, the split** (serial, triangle-index order, as today). A source
point or DEM node goes in as today: `split_inside`, or `split_edge` when it is
exactly on an edge, with the existing rule that skips the split this round if
the neighbour across that edge was already touched. A **strip point** goes in
with `split_edge(t, e, at)` on its sub-edge, under the same skip rule.
Because `at` is off the edge by rounding, the split is guarded as 20b's feet
are: `detail::foot_fits(m, t, e, at)` (in `refine.hpp`, reused, not copied)
must find every child strictly counter-clockwise. If it does not, the point is
marked **refused** for the rest of the run and counted. The vertex's z is the
point's z. `legalise_around` follows, as today.

**Step 4, keeping the sub-edges.** Whenever a split cuts a constrained edge
`(a, b)` at a new vertex `q`, whichever set named the point, the record is
replaced by `(a, q)` and `(q, b)`, with `s_q` set as follows:

- for a strip point, its own `s`;
- for any other point (a source node, or a DEM node, exactly on the edge), the
  strip point of that edge **at exactly the same position**, if there is one:
  that point is then consumed. Otherwise `q`'s projection onto the sub-edge,
  mapped into `[s_a, s_b]`.

Node crossings are snapped exactly (D2, step 4), so a DEM node the rescan
inserts on a constraint takes over that crossing's point. Without this, a
second point a hair from the new vertex would stay and, at tolerance 0, could
be inserted again.

**Step 5, termination.** Each insertion is a stored source point that is not
yet a vertex (15c), a DEM node that is not yet a vertex (14), or a strip point
strictly inside its sub-edge. That strip point becomes an end and is never in
an open interval again. Refused points are never named again. All three sets
are finite.

**Step 6, the end.** The output is built as today. Then one pass over the
final sub-edges evaluates every strip point (refused ones included) against
its sub-edge, which sets `strip_max_error` (not refused) and
`strip_refused_max_error` exactly.

**Outcome** (`PointRefineOutcome`, extended; `RefineOutcome` unchanged):

| field | meaning |
|---|---|
| `strip_points` | strip points in the run (`strip.size()`), 0 without a strip |
| `strip_inserted` | inserted strip points, a subset of `inserted` |
| `strip_max_error` | the largest strip-point error at the end, refused excluded |
| `strip_refused`, `strip_refused_max_error` | refused points, and their largest error at the end |
| `nodes_inserted` | `refine_strip` only: DEM nodes the rescan inserted, a subset of `inserted` |

`max_error` keeps 15c's meaning on the reprojected path (over source points).
In `refine_strip` it is the largest DEM-node error over the triangles the run
rescanned. `uncovered` counts every set's points in void triangles or on void
sub-edges; by the stopping rule it is 0 at the end, and it is reported so that
a change shows.

**`detail::PointScan`** (line 73) gets `z` as `double`, which is exact for
15c's `float` z, and records where its point came from (set, strip index,
`s`), so that steps 3 and 4 can act on it.

**Locality, parallelism, robustness** (the computational-geometry skill, §5).
Every change is local to the triangles beside one sub-edge. The scan stays
read-only and per triangle. The serial phase is unchanged in kind. Every
location decision is an exact predicate (`orient_sign`) or an exact
comparison of parameters. The one step that rounding can upset, inserting a
point that is a hair off its edge, is guarded by `foot_fits` and counted when
refused.

### D5. Types and functions, with their sites

| site | what |
|---|---|
| `include/terrain/refinement/constraint_points.hpp` (new) | `ConstraintPoint`, `ConstraintCheckPoints`, `constraint_check_points` (D2, D3) |
| `include/terrain/refinement/refine.hpp`, `detail::to_lattice` | `detail::lattice_position` extracted from it, and called by it; no change in behaviour |
| `include/terrain/refinement/refine_points.hpp` | `detail::point_loop`; the strip scan and ownership; the sub-edge map; `refine_points(..., strip)`; `refine_strip`; the outcome fields; `PointScan` (D4) |
| `bindings/core.cpp` | `ConstraintCheckPoints` (read-only: `size`, `no_data`, `duplicates`, `edge_count`); `constraint_check_points(view, vertices, edges)` over the bound raster variant; `refine_points(..., strip=None)`; `refine_strip(view, strip, vertices, triangles, z, valid, edges, masks, *, tolerance, threads=0)`; the new outcome properties. Every call releases the GIL, as `refine_points` does today (`bindings/core.cpp:1045`) |
| `src_python/tin_engine/_core.pyi` | stubs for the above |
| `src_python/tin_engine/edge_strip.py` (new) | `generate(view, start, clock) -> ConstraintCheckPoints` and `run(view, strip, start, tolerance, clock) -> PointRefineOutcome`: the two calls and their clock rows. No geometry |
| `src_python/tin_engine/final_check.py:22` | `run(..., strip: ConstraintCheckPoints \| None = None)`, passed on to `refine_points` |
| `src_python/tin_engine/cli.py`, `_dem_mesh` (`:1397`) | after `refine` (`:1477`): `strip = edge_strip.generate(...)`; projected path: `final = edge_strip.run(...)`; reprojected path: `final_check.run(..., strip=strip)`; the sentence and the report (D7) |

`edge_strip.py` exists so that `_dem_mesh`, already about 150 lines, grows by
about 15 rather than 40, and so that the orchestration can be tested with a
fake `_core` (the python-development skill's unit level). Both of its
functions are pure apart from the clock they are given.

### D6. Python and the CLI

On the tolerance path of `_dem_mesh`, in this order:

1. `out = refine(...)`, unchanged.
2. `strip = edge_strip.generate(to_core(tile), out, clock)`. It runs while the
   tile is held. On 15e's branch this is before `del tile`.
3. Reprojected (`grid` and `checks` present):
   `final, n = final_check.run(out, grid, checks, tolerance, clock, strip=strip)`.
   Projected: `final = edge_strip.run(to_core(tile), strip, out, tolerance,
   clock)`. On the projected path the tile is the DEM; it is held through the
   run, as it is today through `refine`.
4. `trim(final...)`, as today.

A refusal from either run is a usage error in the engine's words, as today
(`cli.py:1492-1493`).

### D7. What the file and `--stats` record

The `elevation_source` sentence gains one clause on both paths, after today's
final-check clause:

> `; edge strip: {strip_points} check points where constraints cross grid lines and between them ({no_data} without data), {strip_inserted} inserted, max error {strip_max_error} m at them, {strip_refused} refused (max {strip_refused_max_error} m)`

Projected path: "achieved max error" becomes "max error at DEM nodes at most
{max(out.max_error, final.max_error)} m". This is an upper bound (*default*):
refine reports one maximum, not one per triangle, so the exact figure after
the strip run would need one more full scan. The bound is at most the
tolerance whenever the guarantee holds. The stderr report gains
`{strip_inserted} strip points inserted, {nodes_inserted} nodes inserted by
the strip run`.

`--stats` rows, through the `PhaseClock`: `edge strip: generate` on both
paths; on the projected path `edge strip: scan (parallel)` and
`edge strip: split + flip (serial)`. On the reprojected path, the joint loop's
times stay in 15c's `final check:` rows, because one loop serves both sets.

The sentence changes on every tolerance mesh, so whole-file hashes change
everywhere. `@perf` compares geometry (vertices, triangles, z, edges) rather
than bytes where identity is expected ("Acceptance").

## The guarantee, and its wording

**Wording for the projected path** (Ola's, as ruled, with the default on the
ends written in):

> Every valid DEM node inside the domain is within the tolerance, and so is
> every point where a constraint crosses a grid line, and the midpoint between
> each two neighbouring crossings along a constraint, an edge's ends counting
> as neighbours.

**Wording for the reprojected path** (*default*; it says which grid, because
there are two):

> Every valid source DEM node inside the domain is within the tolerance. So is
> every point where a constraint crosses a line of the resampled grid, and the
> midpoint between neighbouring crossings, measured against the resampled
> grid's bilinear surface.

The second sentence is weaker than the first in one respect, and it says so:
between source nodes the truth is the source DEM. The strip there measures
against the grid that phase 1 refined, which is what 23 already assumed for its
seam pass ("The seam pass on the reprojected path reads the resampled target
grid, as the edge strip does").

**As invariants.** "Strip point" means a point `constraint_check_points` keeps
for the start mesh's constraint edges (D2), at its stored position, with its
stored z. "The mesh value at p" means the output's linear interpolation along
the output constraint edge that holds p.

- **E1. The strip.** For every strip point that is not refused, `|z_p − mesh
  value at p| ≤ tolerance`. The run reports `strip_refused` and
  `strip_refused_max_error`; the guarantee is stated without refused points,
  as 15c's J2 is stated without coincident ones.
- **E2. DEM nodes kept (projected path).** After `refine_strip`, every valid
  DEM node in every closed output triangle with three valid vertices is within
  the tolerance of that triangle's plane, as after `refine` (increment 14's
  guarantee). F2 is why this is an invariant and not an assumption.
- **E3. Source nodes kept (reprojected path).** 15c's J2 holds after
  `refine_points(..., strip)` as it does without the strip.
- **E4. Determinism.** The output is bit-identical for any thread count,
  for any order of the edge list, and for either order within each pair
  (D2, step 2; D4, step 2).
- **E5. Constraints hold.** Every strip insertion splits a constrained edge
  (both halves keep the bit and the mask, as 15c's J9). No constraint is ever
  flipped. No strip point is inserted inside a triangle.
- **E6. Termination.** As D4, step 5. Every inserted strip vertex carries its
  point's own z.
- **E7. Nothing to do, nothing done.** If no strip point is over the tolerance
  and the start is already legal, the output's vertices, triangles, z, valid
  flags and constraint edges equal the start's, bit for bit, and `inserted` is
  0.
- **E8. Only constraint edges.** The strip adds vertices only on constraint
  edges, and adds DEM nodes (projected path) only in triangles it has written.
  Everywhere else the mesh is refine's.

**What is not claimed.** Nothing about a point of a constraint other than the
strip points, with two remarks:

- **On a constraint that lies along a grid line**, on the projected path,
  every point of it is within the tolerance. Along a grid line the bilinear
  surface is linear between nodes, the mesh is linear between vertices, and so
  the error is linear between breakpoints. At nodes it is bounded by E1/E2, and
  at vertices it is 0 (a vertex's z is the surface's there). The midpoints
  are then redundant. This is the property 23 states for grid-line seams.
- **Elsewhere**, F4's 1.25 × tolerance holds on every piece between
  neighbouring strip points that no vertex has split. It is not a guarantee
  (Q1).

## Degeneracy policy

- **A constraint edge along a grid line** (`c0 == c1` or `r0 == r1`, on an
  integer). The family parallel to it yields nothing. The perpendicular family
  meets it only at nodes, snapped exactly (D2, step 4). Its midpoints lie on
  cell sides, with z linear between the two nodes. The remark under the
  guarantee applies.
- **An edge parallel to a grid line but between lines** (a vertical edge at
  col 3.4). One family only; no division by zero.
- **A crossing at a node.** One point, the node exactly, with z equal to the
  node's value: the two crossings of a node are made one by the exact test and
  the de-duplication (D2, steps 4-5).
- **A crossing at an end vertex** (an end on a grid line). Not a check point:
  the bounds are strict, and the vertex's own z is the surface's.
- **A crossing within rounding of an end vertex, but not at it.** Kept. Its
  insertion is guarded by `foot_fits`. If refused, it is counted with its
  error at the end, and excluded from E1.
- **An edge shorter than a cell, or one inside one cell.** No crossing; one
  midpoint between its ends (*default*: ends count as neighbours). Under the
  other reading of the ruling it would have no check point at all.
- **Two crossings closer than rounding that do not coincide** (the edge
  passes a hair from a node). Both kept if their positions and parameters
  differ, and one dropped (`duplicates`) if either is equal. At tolerance 0 both
  may be inserted. Termination is not affected.
- **Coincident points from shared edges.** Two constraint edges of one mesh
  meet only at a shared vertex: the start mesh is a triangulation, and the
  noder guarantees no overlap. So their points coincide only at that vertex,
  which is excluded from both. A road on a lake shore is one edge carrying both
  bits, so it has one set of points.
- **A DEM node or a source node exactly on a constrained edge**, inserted by
  its own set. The sub-edge is split there, and a strip point at exactly that
  position is consumed (D4, step 4).
- **A constraint edge with an invalid end** (NoData under a vertex). The carve
  rule along the edge (D4, step 2). `uncovered` is 0 at the end.
- **A strip point whose cell touches NoData.** Not kept, and counted in
  `no_data`, even when the corner that is NoData has weight 0 in the point's
  cell. That is `vertex_z`'s rule, which refine already applies to every
  off-node vertex. A crossing on the right or bottom side of a NoData cell is
  therefore dropped although its side is valid. Accepted, for one rule
  everywhere.
- **A strip point on the node rectangle's far edge.** `vertex_z` uses the
  last cell, as it does for vertices.
- **A computed coordinate a hair outside the node rectangle.** Clamped (D2,
  step 3).
- **Interior and boundary constraint edges.** One owner each (D4, step 2), so
  no sub-edge is scanned twice or not at all. A boundary edge goes 1 → 2 on
  insertion.
- **Slivers along the constraint after a strip insertion** (a crossing a
  thousandth of a cell from a vertex leaves a very short sub-edge). Allowed.
  Their effect on the worst angle and on the share under 10° is measured at
  acceptance, against the previous increment.

## Tests for `@tester`

**Invariant-critical suite, for mutation testing:** ES2 (E1 by an
independent oracle), ES3 (E2) and ES5 (ownership and sub-edge bookkeeping).
The mutants to run are named with each: the rescan of F2 turned off, the
"no triangle across it" clause of ownership dropped, and the sub-edge map not
updated after a split. Everything else is ordinary. No throwaway
implementation is asked for beyond those three mutants.

**C++ (Catch2), the generator (15f-1):**

- **CC1, one edge, by hand.** An edge between off-node ends chosen so that
  crossings are exact binary fractions. Positions, `s` strictly increasing,
  z equal to the hand-computed linear and bilinear values, `2c + 1` points,
  the ends not among them.
- **CC2, along a grid line.** An edge on column 3 from row 0.5 to row 4.5:
  the points are the nodes `(3, 1)` to `(3, 4)`, exact, with their values, and
  the side midpoints, each the mean of its two nodes.
- **CC3, through nodes.** `(0.5, 0.5)` to `(3.5, 3.5)`: each of the three
  nodes on it appears once, at its exact position. `duplicates` counts the
  merged crossings.
- **CC4, the short and the parallel.** An edge inside one cell gives one
  midpoint. A vertical edge at col 2.25 gives row crossings only. An edge with
  an end exactly on a grid line has no point at that end.
- **CC5, NoData.** Points whose cell has a NoData corner are dropped and
  counted in `no_data`; the others are unchanged.
- **CC6, direction and order.** Reversing every pair and permuting the edge
  list gives, per edge, the same points bit for bit.

**Settled after 15f-1's red step (46e69ad).** The design said nothing on two
points; `@tester` chose for both, and the choices stand as part of the
design:

- **The refusal message** of `constraint_check_points` (`std::invalid_argument`,
  D2) contains `constraint_check_points`, in the style of refine's
  `"refine: ..."` messages.
- **`edge_count()` equals the number of edges given.** Each given edge is
  filed once, as its own entry, in the order given and in canonical direction
  (D2, step 2). An edge with no kept points has an empty `on_edge(k)`.

**Added after `@reviewer`'s round 1 on 15f-1** (the ruling is D2, steps 5-6
and I1-I3):

- **CC7, the ends at ulp distance.** Both of the reviewer's probes become
  cases. `P1` at `nextafter(10, 11)` with a crossing at col 10 gives a `t` of
  exactly 1. `P0` at col `nextafter(1, 0)`, row 1.5, to `(4, 3.25)` gives a
  midpoint that rounds onto the crossing. In both, I1-I3 hold, and
  `duplicates` counts what was dropped. Plus a seeded sweep: ends placed one
  to four ulps either side of grid lines, I1-I3 asserted on every edge.
- **CC8, the degenerate edge.** `{i, i}` is refused, with
  `constraint_check_points` in the message. Two distinct vertices at one
  position give an empty `on_edge(k)`, with the midpoint counted in
  `duplicates`.
- **CC9, out of range.** `edge(edge_count())` and `on_edge(edge_count())`
  throw `std::out_of_range`.

**C++ (Catch2), the loop (15f-2):**

- **ES1, the defect, then its repair.** A start mesh with a long constrained
  edge between off-node ends over a DEM with a bump between nodes, and no node
  in the triangles beside it (Surprise 3 reproduced in miniature). After
  `refine` alone, the oracle of ES2 finds a crossing over the tolerance; this
  is the oracle shown to fail. After `refine_strip`, it finds none.
- **ES2, E1 by an independent oracle.** Random DEMs (seeded) and start
  meshes with several constrained chains, some along grid lines, some through
  nodes, and tolerances including 0. The oracle generates crossings and
  midpoints itself in NumPy (or plain C++ in the test), from the start's
  constraint edges. It finds the output constraint edge holding each one by
  its own projection (distance under 1e-9 cells, parameter in [0, 1]). It
  interpolates linearly, and asserts `|error| ≤ tolerance + 1e-9 · max(1, |z|)`.
  It asserts `strip_refused == 0` for these inputs. It never reads the
  store's order or the loop's records (the computational-geometry skill:
  "borrow the producer's *predicate*, never its records"). It must also be
  shown to fail: plant an output whose strip vertices' z are shifted by twice
  the tolerance. **Mutant:** the sub-edge map not updated after a split.
- **ES3, E2 by brute force.** On the same inputs, every valid DEM node in
  every closed output triangle (barycentric with a 1e-12 slack, as 15c's RP3)
  is within tolerance + 1e-9. **Mutant:** the rescan of F2 turned off. The
  suite must kill it on at least one seed. If it does not, add seeds or a
  steeper DEM until it does, and say which.
- **ES4, determinism.** Threads 1, 2 and 8; the edge list permuted and each
  pair reversed: identical outputs, bit for bit.
- **ES5, ownership.** A boundary constrained edge whose only triangle runs it
  from the higher index to the lower is still scanned and repaired. An
  interior constrained edge is repaired once (one vertex per point, no
  duplicate positions in the output). **Mutant:** the "no triangle across it"
  clause dropped.
- **ES6, refused.** If `@tester` can build a start where `foot_fits` refuses
  a strip point (a crossing within rounding of an end, beside a near-flat
  triangle): the point is counted with its error, the run ends, and E1 holds
  for the rest. If it cannot be built, `@tester` says so and the case is not
  faked.
- **ES7, void.** A constraint edge with one invalid end: carving along the
  edge, `uncovered` 0 at the end, termination.
- **ES8, consumption.** A diagonal constrained edge through a DEM node, set up
  so that the rescan inserts that node: no two output vertices share a
  position, and the crossing's point is not inserted a second time (at
  tolerance 0 as well).
- **ES9, the reprojected path.** `refine_points(store, ..., strip)`: J2 by
  15c's RP3 oracle and E1 by ES2's oracle both hold. The RP1-RP7 suite still
  passes with the store alone, unchanged.
- **ES10, refusals.** A strip on another geometry, or with an edge that is
  not a constraint edge of the start: `std::logic_error`, `RuntimeError` in
  Python.
- **ES11, E7.** A start whose strip points are all within the tolerance:
  output equal to the input, `inserted` 0.

**Python (pytest), 15f-2:**

- **PY1, the binding** (every binding of this increment, the store's and the
  generator's included). Shapes and dtypes refused as `ValueError`; the
  properties; the GIL released (a second thread makes progress during a run,
  as 15c-1's binding tests do).
- **PY2, `edge_strip.py` alone,** with a fake `_core`: the calls are made in
  order, and the clock gets its rows.
- **PY3, the projected path end to end.** A synthetic GeoTIFF and a domain
  with off-node corners and one polygon feature. The `.vtk` sentence carries
  the clause. The oracle is built from the **input polygons**, not the
  producer's records: their crossings with the DEM's grid lines, z from the
  DEM array in NumPy, interpolated along the written constraint lines. Every
  one is within the tolerance. The same oracle on `refine`'s output alone
  (through `_core`, in the same test) finds a violation, so the test can fail.
- **PY4, the reprojected path end to end.** 15c-2's geographic fixture with
  `--out-crs`. The final check's independent source-node check is still 0
  over. The outline's crossings with the target grid's lines, against the
  resampled grid (rebuilt in the test through `target_grid`), are within the
  tolerance.
- **PY5, `--stats`** has the new rows on each path.

## PR split and LOC

Counted in `CLAUDE.md` §2's unit. Estimates, with the worst cases at +39 %
(increment 10's overrun) and +60 % (15a's `mosaic.py`), as in 15c and 23.

| PR | what | est. | +39 % | +60 % |
|---|---|---:|---:|---:|
| **15f-1** | **The generator and its store, C++ only** | | | |
| | `constraint_points.hpp`: `ConstraintPoint`, `ConstraintCheckPoints`, `constraint_check_points` | 90 | | |
| | `refine.hpp`: `detail::lattice_position` extracted | 10 | | |
| | **15f-1 total** | **100** | **139** | **160** |
| **15f-2** | **The strip in the loop, the bindings, and the wiring** | | | |
| | `refine_points.hpp`: `point_loop` (the moved body's changed lines), the strip scan and ownership, the sub-edge map, `foot_fits` and refused points, consumption, the F2 rescan, the end pass, two entry points, the outcome, `PointScan` | 203 | | |
| | `bindings/core.cpp`: `ConstraintCheckPoints` and `constraint_check_points` (moved from 15f-1, ruling of 2026-10-03 below) | 30 | | |
| | `bindings/core.cpp`: `refine_points(..., strip)`, `refine_strip`, the outcome properties | 62 | | |
| | `_core.pyi`: all the stubs, the store's and the generator's included | 35 | | |
| | `edge_strip.py` | 30 | | |
| | `final_check.py` | 5 | | |
| | `cli.py`: the calls, the sentence, the report | 25 | | |
| | **15f-2 total** | **390** | **542** | **624** |

**Ruled by `@architect`, 2026-10-03, after 15f-1's red step (46e69ad): option
(a), the bindings move to 15f-2.** `@tester` found that 15f-1, as first split,
would merge the generator's binding and stubs (about 42 lines) with no test
until PY1 in 15f-2. Option (a) moves them to 15f-2, where PY1 tests them. 15f-1
is C++ only, and its red suite (CC1-CC6 and the `lattice_position` cases) covers
all of it. Option (b) would have moved PY1's binding cases into 15f-1. It was
not chosen, for two reasons: the red suite is already committed, and option (b)
would have added a Python round to a PR that nothing in Python calls yet. Both
PRs stay under 700 at both margins.

**Why two PRs** (*default*). As one PR the strip is about 490 lines: 681 at
+39 % and 784 at +60 %, which breaks the convention 15c and 23 kept (under
700 at both). Split this way, 15f-1 is a pure C++ addition that 23b needs
anyway and that `@tester` can drive with nothing but edges and a raster. 15f-2
carries the behaviour change, every binding, and the acceptance. If Ola prefers one PR and
accepts the +39 % margin only, the split can be dropped with no design change.

**Why more than 15c's 160-210.** F1 needs points filed by edge, with the
sub-edge map and the guarded insertion (about 75 lines). F2 needs the DEM
rescan and the second entry point (about 45). Binding two entry points and a
store is about 105. The midpoints themselves cost about 10, as 15c said.

Module to watch: `refine_points.hpp` grows from 294 lines to about 500. If
it passes 550, the strip scan and the sub-edge map move to
`constraint_points.hpp`, beside the store they read.

**Documentation in the same PRs** (not counted): `ROADMAP.md`'s row at
each merge; a pointer from `15c-geographic-dem.md`'s "The edge strip" and from
23's order (this design branch adds both); `project_structure.md`, where it
lists `include/terrain/refinement/` and the Python modules
(`constraint_points.hpp`, 15f-1; `edge_strip.py`, 15f-2). For 15f-1 both are
done on this branch: the ROADMAP row 15f, and `project_structure.md`'s
`refinement/` list. That list also gains 15c's `check_points.hpp` and
`refine_points.hpp`, which it lacked; that is a documentation defect fixed
here, per the README's rule. `project_structure.md` is not a governed file
(`.claude/hooks/guard_governance.py`, `GOVERNED`), so it was edited during
the unattended window.

## Acceptance (`@perf`)

The rule applies: both PRs touch `include/terrain/refinement/`
(`docs/increments/README.md`, "Acceptance").

**15f-1** changes no behaviour. `tools/bench.py run` (the quarter and the
tile, the thread sweep), back to back with `--tree` against the previous
increment's merge commit. Expected: geometry identical, refine time within
noise. Power state recorded.

**15f-2** changes every tolerance mesh with a constraint, so:

1. **`tools/bench.py`**, as above. Expected: the tile's geometry identical
   (its outline lies along grid lines: E7 and the grid-line remark), the
   quarter's not (its domain outline is off-grid). The strip's own phases
   are reported beside refine, and the end-to-end time compared. Refine
   itself is not edited, so its time should be within noise.
2. **Bygdin with CORINE**, 16c's script
   (`docs/benchmarks/2026-09-29/bygdin-landcover/run.sh`), at 10 m and 1 m
   tolerance, master and 15f-2 back to back, same power state. Record: strip
   points, `no_data`, inserted (strip and nodes), refused, rounds, the
   strip's time, peak memory, worst angle and share under 10° (bench.py's
   `quality()`). The **independent check**, a script of `@perf`'s that does
   not import `tin_engine`: the crossings of the input outline and the
   CORINE borders with the DTM10 grid lines, z bilinear from the DEM,
   against the written mesh's constraint lines. On master it measures the
   projected-path Surprise 3 for the first time (count over the tolerance,
   and the largest error). On 15f-2 it must find 0 over, or exactly the
   refused points.
3. **A basin piece.** The DEM-derived Velhas test catchment (15c's
   acceptance domain) at 20, 10 and 5 m, and one BHO level-3 unit, 766 (the
   smallest peak in `docs/benchmarks/2026-10-02/basin-level3/README.md`,
   2.82 GB at 20 m), at 20 and 5 m, ANADEM from the cache. 15c's source-node
   check stays at 0 over. The independent strip check, against the resampled
   grid, is 0 over. Record the strip's time and memory beside the final
   check's.

Evidence under `docs/benchmarks/<date>/15f-2-acceptance/`.

## Which branch to base on

Base the implementation branch on **15e's tip**, `worktree-agent-ab601ec06c6508cb9`
(*default*), unless 15e has merged by then, in which case base it on master.

- 15f-2's Python edits land in exactly the lines 15e rewrites: `_dem_mesh`'s
  signature (`held`, `grid`, `checks` in place of `opened`), and the
  `del tile` before `final_check.run`. Built on master, the merge would
  conflict there, and F3's ordering (generate before `del tile`) would be
  re-derived in the conflict rather than written once.
- 15f-1 touches nothing 15e touches: `check_points.hpp` is not edited here, and
  15e does not edit `refine.hpp` or `refine_points.hpp`. It could start on
  master. Basing both on 15e keeps one line of history.
- Where 15e's store matters: nowhere in the loop's logic. `refine_points`
  reads `CheckPoints` through `geometry()`, `frozen()` and `for_each_in`,
  which 15e keeps. Memory per point differs: 16 B in 15e's arena against 32 B
  in the strip's store. On the basin's level-3 units the strip is a small
  fraction of the source points. Estimated by `@architect` from the BHO
  outlines in `../rasputin_data/sao_francisco_piece/bho2017_50k_level3/`, in
  the basin's output CRS, at 30 m: about 108,000 strip points on unit 766
  against its 46.2 M check points (0.23 %), and 356,000 on 761 against
  270.4 M (0.13 %). The check-point counts are from `basin-level3/README.md`.
  So the strip needs no arena.
- 15e's line "15d is taken by the edge strip" is fixed in 15e's own PR (see
  "Why 15f"), not in this one.

## Defaults chosen

Each would have gone to Ola; each is the recommended option and can be
overruled without redesign.

1. **The id `15f`**, not `15d` (above).
2. **An edge's ends count as neighbours** for midpoints (proposed in 15c).
3. **The reprojected wording** names the resampled grid (the guarantee
   section).
4. **The strip runs inside one loop with the other point sets**, not after
   them: without it, F2 would break the DEM-node guarantee on the projected
   path, and phase 2 could break the strip on the reprojected one.
5. **Refused strip points** (a guarded insertion that would make a triangle
   that is not counter-clockwise) are counted with their error and excluded from E1, like
   15c's coincident points.
6. **Points whose cell touches NoData** are dropped and counted (refine's
   `vertex_z` rule).
7. **"Max error at DEM nodes at most …"** on the projected path: a bound,
   not the exact maximum (D7).
8. **No CLI switch** to turn the strip off. The guarantee is the point; an
   off switch would print a false tolerance again. `@perf` compares against
   the previous merge with `--tree`, which needs no switch.
9. **Two PRs**, 15f-1 and 15f-2.
10. **Base on 15e's branch.**

## Questions for Ola

None blocks `@tester`.

- **Q1. The every-point form** (F4). Today's ruling guarantees the tolerance
  at the strip points. With about 15-25 more lines, a midpoint's insertion
  could add the midpoints of its two halves, recursively. The error would then
  be bounded on every piece, so every point of every constraint would be
  within 1.25 × tolerance, or within tolerance on constraints along grid
  lines. One caveat: at tolerance 0 with a curved piece the recursion does
  not end, so it needs a floor (say 1/64 of a cell) and an exception for
  pieces at the floor. **Recommended: not now**, ship as ruled, and decide
  with 15f-2's acceptance numbers (how many midpoints are inserted at all). If
  Ola wants it, it is best ruled before `@tester` writes ES2, which would then
  sample every piece densely.

## Review

**15f-1, round 1, 2026-10-03.** Range `390b516..00d9239` (design 4bfeb47, red 46e69ad, ruling 08e2e14, green 00d9239). Verdict: CHANGES REQUESTED. LOC: 128 net (137 added, 9 deleted) against an estimate of 100 (+28 %, inside +39 %). Blocking: (1) "s strictly increasing" is false at the P1 end: t can round to exactly 1.0 when the end lies within an ulp past a grid line, so the last crossing and its midpoint with P1 coincide (probe: ends at col 0.8932792255671602 and nextafter(10, 11), last two points both col 10, row 3.25, s 1, duplicates 0); at P0 a crossing and its midpoint can share a position with distinct s; D2 steps 5-6 to be ruled, then a test, then the fix. (2) Red-step scaffolding at tests/cpp/CMakeLists.txt:327-328. (3) Stale citations in this file at :268, :435, :499 (refine.hpp lines moved by the extraction). (4) project_structure.md and the ROADMAP row, which the design assigns to this PR. Also to record: the .at() bounds checks on edge(k)/on_edge(k), and the duplicate comparison against the previous kept crossing. The lattice_position extraction is behaviour-preserving. Not pushed; no CI.
