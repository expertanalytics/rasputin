# Increment 15f: the edge strip — check points where constraints cross grid lines

Status: **designed by `@architect`, 2026-10-03, while Ola was away
(unattended).** Written before `@tester`, per `docs/increments/README.md`
step 1. **15f-1** (the generator and its store) **merged as #148**
(`@reviewer` APPROVED, round 2, "Review" below; `@perf` ACCEPTED,
`docs/benchmarks/2026-10-03/15f-1-acceptance.md`). **15f-2**, the C++ loop:
the C++ red step is `afd2498` on `worktree-15f-2`, and the gaps it found are
ruled under "Settled after 15f-2's red step" (L1-L9). The C++ green step is
`67df94c` (11 of 13 ES cases); its three open questions are ruled under
"Settled after 15f-2's green step" (L10-L13), and the sliver cascade and the
split into three PRs under "Settled after L12 and L13" (L14, L15), and the
coincidence radius scaled with the lattice (L16). L14 and L16 are green at
`99e95d7`; `@reviewer` APPROVED 15f-2 in round 2 ("Review"), and `@perf`
ACCEPTED it (meshes and quality identical to master, refine within noise;
`docs/benchmarks/2026-10-03/15f-2-acceptance.md`); merged as #152.
**15f-3** (the bindings, the Python and the CLI) merged as #161:
green at `7d841f3`, with the green step's questions ruled
under "Settled after 15f-3's green step" (S1-S5), 179 net lines against 187.
`@reviewer` APPROVED it in round 2 ("Review"), and `@perf` ACCEPTED it with a cost finding (`docs/benchmarks/2026-10-04/15f-3-acceptance.md`), ruled under "Settled after 15f-3's acceptance" (A1-A6): 15f-3 ships, and a follow-up, **15f-4**, makes the mesh rebuild cheap. **15f-4** is implemented on `worktree-15f-4` (green `5eeb87a`, 21 net lines against about 30) and ACCEPTED by `@perf` (`41a5ad7`, `docs/benchmarks/2026-10-04/15f-4-acceptance.md`: meshes byte-identical, refine -3.1 to +0.5 %, the empty-strip call 0.355 s against A2's 0.43 s, end to end +9 to +18 % over the base on the projected path); it is in review as #163. Choices that would normally go to Ola
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
(`include/terrain/refinement/refine_points.hpp:113`) decides which triangle
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
`scan` (`scan`, `include/terrain/refinement/scan.hpp:122`), and inserts nodes as refine
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
  `split_inside`, as phase 2's source nodes do, except one within the radius `r(g)` (L16)
  lattice units of a constrained edge while a strip is present: that one goes
  in on the edge (L12).
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
   (`vertex_z`, `include/terrain/refinement/scan.hpp:79`), the function refine uses for an
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
  round 1, every triangle the run has written so far is scanned. That is not
  the same as the active set: a skipped or refused triangle is active without
  having been written. The gate is ruled in L3.
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
must find every child strictly counter-clockwise, and the split's new edges
must be locally Delaunay (`detail::strip_fits`, L10). If it does not, the point is
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
| `bindings/core.cpp` | `ConstraintCheckPoints` (read-only: `size`, `no_data`, `duplicates`, `edge_count`); `constraint_check_points(view, vertices, edges)` over the bound raster variant; `refine_points(..., strip=None)`; `refine_strip(view, strip, vertices, triangles, z, valid, edges, masks, *, tolerance, threads=0)`; the new outcome properties. Every call releases the GIL, as `refine_points` does (the `py::gil_scoped_release` in its binding, `bindings/core.cpp:1182`) |
| `src_python/tin_engine/_core.pyi` | stubs for the above |
| `src_python/tin_engine/edge_strip.py` (new) | `generate(view, start, clock) -> ConstraintCheckPoints` and `run(view, strip, start, tolerance, clock) -> PointRefineOutcome`: the two calls and their clock rows. No geometry |
| `src_python/tin_engine/final_check.py:28` | `run(..., strip: ConstraintCheckPoints \| None = None)`, passed on to `refine_points` |
| `src_python/tin_engine/cli.py`, `_dem_mesh` (`:1466`) | after `refine` (`:1550`): `strip = edge_strip.generate(...)`; projected path: `final = edge_strip.run(...)`; reprojected path: `final_check.run(..., strip=strip)`; the sentence and the report (D7) |

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
(`_dem_mesh`'s `typer.BadParameter(f"{dem}: {out.message}", ...)`, `cli.py:1604`, and the twin after either run at `:1623`, which reads `final.message`).

### D7. What the file and `--stats` record

**Superseded for the file, stderr and `--stats` by increment 25**
(`docs/increments/25-plain-output.md`, D6; Ola's ruling of 2026-10-03 that the
`elevation_source` sentence is replaced by named fields before 15f-3's code
step, and Ola's cut of the same day: the mesh file carries only what a user of
the mesh needs). The strip's counts below go to `--stats` as 25's
`line_points_checked`, `line_max_error_m`, `line_points_refused` and the rest
(25, D6), not to the file; the file's `max_error_m` carries the "at most"
bound. The `--stats` phase rows below stand.

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
  as 15c's J2 is stated without coincident ones. On the reprojected path it
  is also stated without a strip point consumed by a source point (L5): there
  the source's z is the vertex's, and J2 governs.
- **E2. DEM nodes kept (projected path).** For a start that is `refine`'s output at the same tolerance (the only start the CLI gives it; S1), after `refine_strip`, every valid
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
- **ES8, consumption** (as amended by L7). A diagonal constrained edge
  through a DEM node. Projected path: the node is a vertex exactly once,
  whichever set inserts it, and no two output vertices share a position (at
  tolerance 0 as well). Reprojected path: a source point at that node wins the
  tie (D4, "Combining"), its split consumes the strip point there (step 4), and
  the strip point is not inserted a second time.
- **ES9, the reprojected path.** `refine_points(store, ..., strip)`: J2 by
  15c's RP3 oracle and E1 by ES2's oracle both hold. The RP1-RP7 suite still
  passes with the store alone, unchanged.
- **ES10, refusals.** A strip on another geometry, or with an edge that is
  not a constraint edge of the start: `std::logic_error`, `RuntimeError` in
  Python.
- **ES11, E7.** A start whose strip points are all within the tolerance:
  output equal to the input, `inserted` 0.

`@tester` added two cases, and the design adopts them: **ES12**, E1, E2 and
the Delaunay oracle at UTM scale with non-square cells; **ES13**, ends within
ulps of grid lines, with I1-I3 seen at the loop (no coincident vertex, no
triangle that is not strictly counter-clockwise, every point over the tolerance
within 1e-9 cells of an end and counted as refused).

**Settled after 15f-2's red step (afd2498).** `@tester` found places where
D4 says nothing or is wrong; ruled by `@architect`, 2026-10-03 (Ola away, the
recommended option taken in each). They are part of the design, and
`@developer` implements them as written here.

- **L1. A refused strip point leaves its triangle active.** `@tester`'s choice
  (ES6) stands. In the split phase a strip point is handled in this order:
  1. `touched[t]` set: `continue`, as today;
  2. the skip rule (the neighbour across its sub-edge already touched): push
     `t` onto `skipped`, as today;
  3. `foot_fits(m, t, e, at)` false: mark the point refused, add 1 to
     `strip_refused`, push `t` onto `skipped`, and leave `any` true. Nothing is
     inserted, and there is no fallback (refine's refused foot falls back to
     its node; a strip point has none).

  The skip rule comes before `foot_fits` because deferral is retried and
  refusal is permanent: a point is never refused in a round in which it would
  not have been split anyway. The refused marks are one flag per point of the
  strip, written only in the serial phase and read by the next parallel scan,
  so there is no race. Next round `t` is rescanned, its refused point is
  skipped (D4, step 2), and it names its next worst point, of any set. Each
  round of a refused triangle either inserts or refuses a point that was never
  named before, so step 5's termination argument holds.
- **L2. Refusals: order, statuses and texts.** `@tester`'s choices stand. The
  checks run in this order, and the first that fails decides:
  - `refine_strip`: (1) a tolerance that is not finite and ≥ 0 returns
    `RefineStatus::InvalidTolerance` with "refine_strip: tolerance must be
    finite and >= 0", no throw, as `refine` and `refine_points` do; (2)
    `strip.geometry()` different from `dem.geometry()`: `std::logic_error`;
    (3) a strip edge that is not a constraint edge of the start:
    `std::logic_error`; (4) `to_lattice`'s refusals, as `refine`.
  - `refine_points(..., strip)`: (1) the tolerance and (2) the frozen store as
    today; then (3) and (4) as above, against `points.geometry()`.

  Every `std::logic_error` text starts with the entry point's name and a
  colon ("refine_strip: ...", "refine_points: ..."), in the style of
  `"refine_points: the check-point store is not frozen"`; only the name is
  pinned. "Differs" means any of the six fields of `RasterGeometry` (x_min,
  y_max, delta_x, delta_y, cols, rows) compares unequal. `RasterGeometry`
  has no `operator==` today. A defaulted one in
  `include/terrain/raster/geometry.hpp` (one line, as `CellIndex` has) is the
  recommended way. "A constraint edge of the start" means the unordered pair
  appears in `edges`, whatever its mask. It is checked against that span,
  before `to_lattice`.
- **L3. Which triangles the DEM rescan reads (corrects D4, step 2).** In
  `refine_strip` a triangle is scanned for DEM nodes exactly when its slot
  has been **written** in this run. "Written" is one flag per slot that, once
  set, stays set. `legalise_all`'s callback sets it (today that callback is
  `[](std::uint32_t) {}`). So does every slot `touched` in a round, and every
  new slot a split creates. Skipped and refused triangles that were not
  written are scanned for strip points only. This makes E8 ("DEM nodes only
  in triangles it has written") literally true. The other reading, "every
  active triangle after round 1", would scan refine's own triangles when they
  are skipped or refused. On a start that did not come from `refine` at the
  same tolerance, that could insert nodes outside the written triangles.
- **L4. What the outcome reads from the scan.** `PointScan` keeps two
  things apart. One is the winner across sets (point, set, strip index, `s`,
  where), for steps 3 and 4. The other is the non-strip set's own largest
  error (source points on the reprojected path, DEM nodes in `refine_strip`)
  and the `uncovered` counts. At the end, `max_error` is the largest of that
  non-strip error over every slot's last result. A slot with no DEM scan
  counts as 0. `strip_max_error` and `strip_refused_max_error` come from step
  6's pass only, never from scan results. `uncovered` sums every set's count,
  the strip's counted at the owning triangle.
- **L5. Consumption in detail (D4, step 4).**
  - A sub-edge record keeps its ends in the order of `s`, so `s_a < s_b`,
    with `a` the end nearer `P0`.
  - The strip point at exactly `q`'s position is looked up among the points
    with `s` in `(s_a, s_b)`, by exact equality of `at`. A linear walk is
    enough, because the path is taken only for a node or a source point that
    lies exactly on a constrained edge. Binary search on `q`'s projected `s`,
    with a check of both neighbours, is also allowed.
  - Without such a point, `σ = clamp(((q − a) · (b − a)) / |b − a|², 0, 1)`
    in lattice `(col, row)`, and `s_q = s_a + σ (s_b − s_a)`. If `s_q` rounds
    onto `s_a` or `s_b`, one half has an empty interval. That is allowed: it
    holds no point to miss, and a point sitting at `s_q` itself is still
    measured by step 6.
  - A consumed point is not inserted, is not counted in `strip_inserted`, and
    is never named again (it is now an end).
  - Step 6 measures a point whose `s` equals an end's `s` against that end's
    z. On the projected path a consumed point's error is then exactly 0,
    because both z's are `value_at` of the node. On the reprojected path a
    source point that lands exactly on a strip point keeps the source's z,
    because the source is the truth there (J2). The difference goes into
    `strip_max_error` and is reported. E1 is stated without such points (E1
    above). For a source point to land exactly on a crossing or a midpoint
    needs an exact coincidence, which no ES case builds and no real input is
    expected to produce. It is ruled so that the reported figure stays honest.
- **L6. Counting.** `carved` counts the void insertions of every set.
  `strip_inserted` counts every strip insertion, carved ones included, and
  `nodes_inserted` counts every DEM-node insertion of `refine_strip`, carved
  ones included. `inserted` is their sum with the source insertions.
  `strip_refused` counts points, each at most once. An inserted strip vertex
  is output at `(x_min + col dx, y_max − row dy)`, valid, with the point's z,
  as `refine_points` outputs a source point.
- **L7. ES8: the substitution is accepted, and D4 is not changed.** The
  order is already fixed: the serial phase walks triangles in index order, and
  E4 holds. What a test cannot fix is which slot holds which triangle after a
  split, because that is `LatticeMesh`'s numbering, which the design leaves
  open on purpose. On `@tester`'s projected fixture the strip names the node
  in round 1: no triangle has been written, so there is no DEM scan, and ties
  go to the strip anyway. A DEM node can reach consumption only in a later
  round, and only if the triangle that does not own the sub-edge comes first
  in index order and names the node. A test that forced this would pin the
  slot numbering. Step 4 runs one code path for "any other point", source or
  DEM node, so the reprojected section exercises it. The projected section
  keeps "one vertex at the node". Optional for `@tester`, not blocking: add
  `CHECK(out.strip_refused == 0)` to the projected section. At tolerance 0,
  if the DEM route were ever taken and consumption were missing, the
  leftover point would show as a refused insertion on top of the vertex.
- **L8. The oracles (`tester.md` §3D), confirmed.** Every property case on the
  projected path carries three oracles: the constrained-Delaunay oracle
  (§3D), the DEM-node tolerance oracle (§3D; it is E2), and the strip oracle
  (E1, ES2's). On the reprojected path the Delaunay and strip oracles stay.
  The tolerance oracle there is 15c's J2 by RP3's oracle, over the source
  points. This is §3D's oracle read against the DEM that path takes as the
  truth: the source DEM, whose nodes the check points are. Phase 2 gives up
  the resampled grid's nodes on purpose (J2, F2), so no oracle checks them.
  ES2 and ES3 above name only their own oracle. The suite carries all three on
  every projected case (`projected_oracles` in
  `tests/cpp/property/prop_refinement_edge_strip.cpp`), and that is what is
  required.
- **L9. Later, by day: TSan.** Once `prop_refinement_edge_strip` is green,
  it joins the TSan target list in `.github/workflows/main.yaml` (the `tsan`
  job's `--target` list). The strip scan is a new parallel reader of shared
  state: the sub-edge map and the refused flags. That file is not edited from
  the unattended window; it is for 15f-2's PR, with Ola present.

LOC effect: L1, L3, L4 and the `operator==` of L2 add about 15 lines to
`refine_points.hpp`'s 203 and 1 to `geometry.hpp`. 15f-2's estimate becomes
about 406: 564 at +39 % and 650 at +60 %, still under 700.

**Settled after 15f-2's green step (67df94c).** `@developer` raised three
points; ruled by `@architect`, 2026-10-03. None needs Ola. L10 and L12 make
refusals rarer and no claim wider, Ola's wording of the guarantee is
unchanged, and refused points stay excluded from E1 and reported (default 5).

- **L10. `strip_fits` is accepted as a departure from D4, step 3.** A strip
  point is guarded by `detail::strip_fits`: `foot_fits`, and then the two new
  edges from `q` (`q`-`c`, and `q`-`d` across the sub-edge) must be locally
  Delaunay by `must_flip`'s test. The reason: a point that is a hair off its
  sub-edge on the far side can lie outside `t`'s circumcircle. The circle's
  segment beyond the chord is thin near the chord's ends, so this happens
  only there. `legalise_around` tests only the edges opposite `q`, so `q`-`c`
  would stay non-Delaunay, and it cannot be flipped: the flip would cross the
  constrained chain `a`-`q`-`b`. Refusing is the only exact outcome. The
  alternative, inserting anyway and leaving a non-Delaunay edge, breaks
  `tester.md` §3D's Delaunay oracle; without the check ES13 failed on 9
  sweep inputs, and with it on none of 240 (`@developer`, `67df94c`). Its
  ~20 lines were not in the estimate; they are counted below.
  `strip_fits` also guards L12's insertions.
- **L11. ES6: the test changes, and the oracle does not.** `ruled_points`
  merges crossings closer than 1e-12 in parameter. That is right for an
  independent oracle: it cannot reproduce the generator's ulp-level arithmetic
  without copying it, which the computational-geometry skill forbids ("borrow
  the producer's predicate, never its records"). So the refused point `f`
  (`t` ≈ 3.7e-17) is not one of the oracle's points, and
  `REQUIRE(f.over > 0)` (`prop_refinement_edge_strip.cpp:465` at `afd2498`) asserts
  something the oracle cannot see. `@tester` changes the refused section:
  - keep `strip_refused == 1` and `strip_refused_max_error ≈ 10`;
  - replace `REQUIRE(f.over > 0)` and the loop over `over_at` by
    `CHECK(f.over == 0)`: every point the oracle generates is repaired, which
    is "E1 holds for the rest";
  - measure `f` from the test's own knowledge: a hand-built one-point list
    `{OraclePoint{{1.0, 1.5}, 6.0, 0}}` (the plane `3 col + 2 row` there),
    through `strip_findings` at the case's tolerance. It must give `over == 1`
    (its error, about 10, is `strip_refused_max_error`). The position is
    the CC7 probe's, known by hand, so nothing is read from the store.
- **L12. A non-strip point a hair off a constrained edge goes in on that
  edge.** This is option 1 of the three `@developer` listed. In the ES13
  sweep the rescan inserted a DEM node lying within ulps of the constraint
  with `split_inside`. That leaves the sliver F1 describes, the node becomes
  the apex over the sub-edge, and `foot_fits` then refuses strip points in
  the middle of the edge, far from any end. F1's reason for filing strip
  points by edge applies to such a node as well. The rule, in the split phase
  of `point_loop`, **only when a strip is given**:
  - A source point or DEM node that the scan places `Inside` triangle `t` is
    tested against each constrained edge `e` of `t`. It is a candidate when
    its distance to the line of `e` is at most **the coincidence radius `r(g)`** (L16; 1e-10 lattice units as first ruled)
    (Euclidean in `(col, row)`) and its projection falls strictly inside
    `e`. The threshold is a decade inside the oracles' `kOnEdge` (1e-9), so
    the producer and the oracle never disagree at the boundary. It is still
    many ulps at the basin's lattice coordinates (about 3e4).
  - With more than one candidate, take the nearest, and on a tie the lower
    edge index in `t`. The point then goes through L1's order with
    `split_edge(t, e, p)` at **its own position** and its own z: touched,
    then the skip rule, then `strip_fits`. If `strip_fits` refuses, it falls
    back to `split_inside`, as today. Nothing is marked refused, because E2
    and J2 need that point in the mesh. The cut then follows step 4, with
    consumption (L5).
  - Without a strip nothing changes, so the no-strip path stays bit-identical
    (`@developer`'s 72 fixtures).

  The other two options are rejected. An L1 fallback has nothing exact to
  fall back to: a strip point has no other position. Accepting the refusals
  in the test would turn a refusal from a rounding-scale event at a vertex
  into a loss of E1 in the middle of an edge, with an error of any size.

  **What the guarantee then says.** E1 is unchanged: every strip point that
  is not refused is within tolerance. Refused points are counted and their
  largest error is reported. What changes is the expectation stated for
  refusals: they are expected only within rounding (1e-9 lattice units) of a
  vertex of the output, which is either an end of a start edge or a vertex the
  run put on the constraint. One case stays outside that expectation: a
  **start** vertex that lies within rounding of another constraint edge
  without being on it. This can only come from the input. A refusal it causes
  is still counted and reported, and `@perf`'s acceptance records
  `strip_refused` and its largest error on Bygdin and the basin pieces.

  **Who changes what.**
  - `@tester`, first, in one commit with the reason:
    - ES13's assertion on `over_at` becomes "within 1e-9 lattice units of an
      **output** vertex", not only of a start end. A crossing within ulps of
      a node that went in on the edge may be refused at tolerance 0.
    - A new case, **ES14**, promotes the sweep input seed 11, k 0, as an
      explicit start. At tolerance 0.5 it asserts `strip_refused == 0`, the
      strip oracle at 0 over, and the DEM node near (8, 8) an end of an
      output constraint edge (on the chain, so not `stray` and no
      `broken_chain`). At tolerance 0 it asserts every over point within
      1e-9 of an output vertex. ES14 is red until `@developer` lands L12,
      and it is what keeps the widened ES13 from hiding the original failure.
    - Optional, from L7: `CHECK(out.strip_refused == 0)` in ES8's projected
      section.
  - `@developer`, then: L12 in `point_loop`'s split phase, and nothing else in
    behaviour. If ES14 or the ES13 sweep still shows refusals in the middle of
    an edge after L12, report the input back rather than widening the
    threshold.
- **L13. Where the strip machinery goes (corrects "Module to watch").**
  `refine_points.hpp` is 582 physical lines at `67df94c`, which is past the
  550 trigger. The design's 294 baseline was `wc -l` at `390b516`, so the
  trigger is `wc -l` too. It counts 552 non-blank lines and 470 lines that
  are neither blank nor comment; the "534" figure matches none of these
  counts at `67df94c`. The strip machinery (`SubEdge`, `SubEdges`,
  `edge_key`, `along`, `scan_strip`, `cut`, `strip_fits`, and L12's
  candidate test) moves to a new header,
  `include/terrain/refinement/strip_scan.hpp`, and **not** to
  `constraint_points.hpp` as the design first said. That header is the
  generator, which 23b's seam pass reuses without the loop. Moving the loop's
  sub-edge map, the hash map and the Lawson predicates into it would hand
  23b dependencies it does not use. `refine_points.hpp` includes the new
  header. This is `@developer`'s, as its own commit after L12 is green: a
  pure move, with no test file and no behaviour change (the suites and the
  72-fixture identity as the check). `project_structure.md`'s `refinement/`
  list gains the header in the same PR.

LOC effect of L10-L13: `strip_fits` about 20 (already in `67df94c`'s 229
net), L12 about 15, the new header's includes and namespace about 15. The C++
lands near 260 against the design's 218 (+19 %). 15f-2's estimate becomes
about 450: 626 at +39 % and 720 at +60 %. That is over 700 at the wider
margin only, so the convention 15c and 23 kept, under 700 at both margins, no
longer holds for 15f-2 on estimates. The ceiling is on actual lines, and
`@reviewer` counts them. The Python half (bindings, stubs, `edge_strip.py`,
CLI, about 187 estimated) is what is left to land. If the measured C++ plus
that half approaches 700, the bindings and Python move to a 15f-3, which is
the split option (b) that was considered for 15f-1.

**Settled after L12 and L13 (3ff8f1d, 61cdbaa).** Ruled by `@architect`,
2026-10-03.

- **L14. A non-strip point within rounding of a vertex is not inserted.**
  The input: ES13's sweep, seed 11, k 37, at tolerance 0. It fails when the
  build has `-ffp-contract=off`, which is likely how CI's x86-64 runner
  builds. The cascade, as `@developer` traced it:
  1. DEM node (17, 9) lies within ulps of the edge's end `b`. Its projection
     is not strictly inside the sub-edge, so L12 does not apply, and it goes
     in by `split_inside`.
  2. That leaves an unconstrained sliver edge, (14, 7.5000000000000009) to
     (17, 9), lying almost on the constraint.
  3. Node (15, 8), which lies on the constraint's line, now sits in the
     triangle on the far side of the sliver, a triangle with no constrained
     edge. L12 finds nothing and it goes in by `split_inside`.
  4. (15, 8) becomes the apex over the sub-edge, and `foot_fits` refuses the
     strip point (15.5, 8.25) in the middle of the edge.

  The root is step 1: a vertex inserted a few ulps from an existing vertex.
  It makes a sliver that every later local rule has to see through.

  **The rule.** When a strip is given, the scan treats a point of the
  non-strip set (a DEM node in `refine_strip`, a source point in
  `refine_points`) that lies within **the coincidence radius `r(g)`** (L16; 1e-10 lattice units as first ruled) of a corner of
  its triangle as it treats a point equal to that corner. The point is
  skipped: never named, never inserted. This is a geometric test inside the
  read-only scan, with no state and no set of marks. `scan` (in `scan.hpp`)
  and `scan_points` each take a radius that defaults to 0. `refine` passes
  nothing and is bit-identical; so is `refine_points` without a strip. After
  the loop, the end pass finds such points per vertex and reports them in
  `coincident` and `coincident_max_error`:
  - for DEM nodes, it looks at the node nearest each vertex;
  - for sources, it extends 15c's existing pass from exact equality to the
    radius, over every vertex, not only start vertices.

  In `refine_strip` those two fields have no other use. Their meaning there
  is "DEM nodes within `r(g)` (L16) of a vertex they are not; the
  largest |node z − vertex z|".

  With this rule, (17, 9) is never inserted, so no sliver forms. Node
  (15, 8) is not left off the chain, either. As built (`99e95d7`), the strip's
  row-8 crossing goes in on the chain at (14.999999999999998, 8), about
  2e-15 from the node. The node is then within `r(g)` of that vertex, so the
  scan skips it and counts it in `coincident`. (This sentence first said
  L12 puts (15, 8) itself on the edge. Which of the two happens depends on
  which set names the point first. Either way the constraint is represented
  at (15, 8) by a vertex on the chain, which is what the case guards.)

  **The rejected options.**
  - L12 searching for constraints that are not edges of `t`, across slivers:
    this needs a walk of unbounded length and is no longer local
    (computational-geometry skill, §5), and it treats the symptom, not the
    sliver.
  - Accepting the refusal: (15.5, 8.25) is about 0.56 cells from the nearest
    vertex, so this would be the mid-edge loss of E1 that L12 already
    refused to accept.

  **What the guarantees then say.**
  - **E2 gains an exception at rounding scale.** "Every valid DEM node
    inside the domain is within the tolerance", except a node within `r(g)` (L16)
    lattice units of a vertex that it is not. Such a node is counted in
    `coincident`, with its largest difference in `coincident_max_error`.
  - That difference is the bilinear surface's change over at most `r(g)`
    cells, plus the vertex's own rounding. It is far below any tolerance the
    CLI accepts above 0. Only at tolerance 0 can such a node exceed the
    tolerance at all.
  - **E3 (J2) gains the same exception for source points**, in the form 15c
    already states it for exact coincidences. A source point near a vertex
    can differ from the vertex's z by any amount, and `coincident_max_error`
    reports it, as it does today.
  - **E1 is unchanged.** L12's expectation for refusals stands: only within
    1e-9 lattice units of an output vertex, apart from start vertices that
    lie within rounding of another constraint.
  - Any other refusal is counted and reported. An ES case that finds one is
    reported back as a design defect, and the case's tolerance is not
    widened.
  - **ASK OLA (not blocking):** E2's exception changes the wording Ola ruled
    ("every valid DEM node ... within the tolerance") at a scale of `r(g)`, 1e-10 to about 5e-10
    cells. The figure is reported, so nothing is hidden. If Ola wants it
    stated in the `.vtk` sentence, it costs one clause in 15f-3's CLI work.
    *Answered by increment 25* (`docs/increments/25-plain-output.md`, D2;
    default taken while Ola was away): the file's `max_error_m` includes the
    largest difference at these nodes, and their count is in `--stats`.
  - **ASK OLA (not blocking), a separate question:** should the build pin
    `-ffp-contract=off` project-wide? Otherwise arm64 (FMA contraction on)
    and x86-64 can differ in the last bits of any refine output. Bit-identity
    across platforms is not claimed today (E4 is about threads and edge
    order). Pinning it changes every stored hash on arm64, so it is Ola's
    call and outside 15f.

  **Who changes what.**
  - `@tester`, first, red:
    - a case **ES15** with seed 11, k 37's exact `a`, `b` and apex as
      literals, at tolerance 0 and 0.5. It asserts: no strip point over the
      tolerance farther than 1e-9 lattice units from an output vertex; node
      (17, 9) not a vertex; `coincident ≥ 1` and
      `coincident_max_error ≤ 1e-9 · max(1, |z|)`; at least one vertex within
      `r(g)` of node (15, 8), and every such vertex an end of an output
      constraint edge (as amended in `fadcc77`; the node itself is a vertex
      only if L12 inserts it before the strip's crossing goes in); E2 by
      `node_findings`, whose 1e-9 relative
      slack already covers the exception; and the Delaunay oracle;
    - a **unit case for the radius**: a DEM node 1e-11 lattice units from a
      start vertex is not inserted at tolerance 0 and is counted, and one at
      1e-8 is inserted;
    - a **build check**: before handing back, run `prop_refinement_edge_strip`
      in a build configured with `-DCMAKE_CXX_FLAGS=-ffp-contract=off` as
      well as the default one, and name both in the handback.
  - `@developer`, then: the radius parameter on `scan` and `scan_points`
    (default 0), passed as `r(g)` (L16) by `point_loop` when a strip is given; the
    end-pass counting; nothing else. Run both builds, `-ffp-contract=off` and
    the default, before handing back, and confirm that `refine` and the
    no-strip `refine_points` are bit-identical (the 72 fixtures). About 15
    lines.

- **L15. Three PRs: the bindings and Python go to 15f-3.** As L13 foresaw.
  Measured by the line rule of `CLAUDE.md` §2 (blank and `//` lines dropped):
  - 15f-1 is 136 net (`390b516..569b00c`);
  - 15f-2's C++ is 275 net so far (`569b00c..61cdbaa`), about 290 with L14;
  - 15f-2 with the Python half (about 187 estimated) would be about 477, or
    about 590 at +60 % on the Python;
  - the "about 712" figure counts 15f-1 into the same PR (136 + 275 + 187 at
    +60 % ≈ 710), because the branch is stacked on 15f-1, which has not been
    published.

  Either way, the convention 15c and 23 kept no longer holds for one PR.
  Ruled:
  - **15f-1** stays its own PR (approved and accepted), and is published
    and merged first;
  - **15f-2** is the C++ loop alone, as the red step already scoped it: no
    caller in the CLI, so no tolerance mesh changes, and its acceptance is a
    no-change run, as 15f-1's was. It also carries the TSan entry (L9);
  - **15f-3** carries every binding, the stubs, `edge_strip.py`,
    `final_check.py`, the CLI, PY1-PY5, and the full acceptance (Bygdin,
    the basin piece). It is based on 15e's tip, or master once 15e has
    merged ("Which branch to base on"), because its Python edits are the
    ones that meet 15e.

  The PR table and the acceptance below are updated to match.

- **L16. The coincidence radius scales with the lattice: `r(g) = max(1e-10,
  64 · ulp(M))`.** This replaces the fixed 1e-10 of L12 and L14 (Ola asked
  for the check, 2026-10-03). Here `M = max(cols, rows) − 1` is the largest
  lattice coordinate the geometry allows, and `ulp(M)` is
  `std::nextafter(M, inf) − M`. One function,
  `detail::coincidence_radius(const raster::RasterGeometry&)`, is used by
  both L12 (`near_constraint`) and L14 (the scan's corner test).

  **1. How large lattice coordinates get.**
  - Lattice coordinates are measured from the corner of the raster handed to
    C++, not from any global origin:
    `col = (x − g.x_min()) / g.delta_x()` (`detail::lattice_position`,
    `include/terrain/refinement/refine.hpp`).
  - Under 23 that raster is the piece's own target window, snapped to the
    global lattice but with its own `x_min` and `y_max` ("Unchanged by
    23a-1", `23-basin-scale.md`; `window_meta`: `x_min + col0·dx`). So a
    coordinate is at most the extent of one window, never of the basin.
  - Windows measured or designed so far:
    - the Norway DTM10 tile, 5,051 × 5,051 (`7908_3_10m_z33.tif`, the bench
      DEM);
    - 15a's mosaic box, 10,051 × 10,051;
    - the São Francisco level-3 units at 30 m, the largest being 761 at
      19,377 × 27,786 (`docs/benchmarks/2026-10-02/basin-level3/runs/t20/761.stats.md`);
    - the whole basin as a single window, 41,332 × 50,297 at 30 m (23's
      partition table). This is reached only when the memory budget allows
      one piece, and at the default budget even 50 m tolerance cuts it 2 × 2.
  - Under 23's partition at 1 m tolerance the basin is 5 × 7 cells of about
    59.4 M nodes, so pieces are near 8,300 × 7,200 nodes.
  - "1 m" in 23 and in `bench.py` is a tolerance, not a cell size. No 1 m
    DEM is in use. A 1 m DEM over the whole basin, about 1.2 M × 1.5 M nodes,
    is not a window any path could hold. Under 23 its pieces would again be
    windows of about `sqrt(budget / b)` nodes per side.
  - Nothing in the code caps a window's aspect, though. A long, thin
    corridor at 1 m could be 10⁶ nodes long and still fit in memory, so the
    rule must not depend on today's sizes.

  **2. What the rounding is proportional to.** It is proportional to the
  size of the coordinates, not to local differences.
  - Strip points are computed in absolute lattice coordinates: a crossing is
    `(K, r0 + t (r1 − r0))` (D2, step 3). The result has the magnitude of the
    coordinate, and so does its rounding, up to a few `ulp(X)`.
  - `near_constraint` and the scan form differences such as `p − a` first.
    For nearby points those differences are nearly exact, so they add little
    rounding of their own. But the positions they subtract already carry
    their representation error of `ulp(X)`. A vertex inserted a few ulps
    from another one is a few `ulp(X)` away.
  - So the quantity the threshold must exceed grows with `X`, and a fixed
    constant runs out:
    - at the whole-basin window, 1e-10 is about 14 ulps;
    - at 10⁶, one ulp (1.16e-10) is already larger than 1e-10.
  - A second source of offset exists: a world coordinate's own rounding
    (about 9.3e-10 m at a UTM northing near 7-8 × 10⁶ m), divided by the
    cell size. It decides how near an input vertex can come to a node
    without being one. It does not drive the sliver failure: separating two
    points 1e-9 cells apart is easy for the exact predicates. The failure
    comes from lattice arithmetic, so the radius tracks `ulp(X)`.

  **3. The rule and its margins.**
  - 64 ulps is 16 times the few-ulp offsets seen in the ES13 and ES15 inputs.
  - The floor of 1e-10 keeps today's value wherever it was already sound: on
    every lattice up to 8,191 nodes across, `r(g) = 1e-10`. That covers every
    test lattice, the DTM10 tile and the pieces at 1 m tolerance.
  - The extent `M` is used rather than each triangle's own coordinates. It is
    one number per run, uniform and deterministic, and errs on the large side
    near the origin, where 64 `ulp(M)` is still physically nothing (under
    5e-9 m at 10 m cells, up to `M` = 5 × 10⁴).
  - Values:

    | lattice extent `M` | `ulp(M)` | `r(g)` | `r(g)` in ulps |
    |---:|---:|---:|---:|
    | 8 (the ES lattices) | 1.8e-15 | 1e-10 (floor) | about 56,000 |
    | 5,050 (DTM10 tile) | 9.1e-13 | 1e-10 (floor) | 110 |
    | 10,050 (15a box) | 1.8e-12 | 1.16e-10 | 64 |
    | 27,785 (unit 761) | 3.6e-12 | 2.33e-10 | 64 |
    | 50,296 (whole basin, one window) | 7.3e-12 | 4.66e-10 | 64 |
    | 10⁶ (a corridor) | 1.16e-10 | 7.45e-9 | 64 |

  **The oracles.** The tests' "on the edge" distance, `kOnEdge` (1e-9 cells,
  `tests/cpp/support/strip_oracle.hpp`), must stay at least `10 · r(g)` for
  the lattice under test, so that producer and oracle never disagree at the
  boundary. On every lattice in today's suites `r(g)` is 1e-10, so 1e-9 is
  exactly a decade inside and nothing changes. A test on a lattice wider than
  8,191 nodes uses `max(1e-9, 10 · r(g))`. So does @perf's independent strip
  check in 15f-3's acceptance (basin pieces, up to about 2.8 × 10⁴ across,
  where `r(g)` is 2.33e-10).

  **Who changes what.**
  - `@tester`: the ES15 radius case needs no change. Its lattice is 9 × 9,
    where `r(g)` is the 1e-10 floor, so 1e-11 is inside and 1e-8 outside, as
    written. Add, red:
    - a unit case for `coincidence_radius`: `M` = 8 gives exactly 1e-10;
      `M` = 50,296 gives exactly 64 × 2⁻³⁷; `M` = 2²⁰ gives exactly 64 × 2⁻³²;
    - one loop case on a lattice wider than 8,191 (for example 16,385 × 9
      nodes, where `r(g) = 64 · 2⁻³⁸ ≈ 2.33e-10`): a start vertex 2e-10 from
      a node is skipped and counted there, while the same offset on a 9 × 9
      lattice is inserted. This shows the radius scales at the loop, not
      only in the helper.
  - `@developer`: `detail::coincidence_radius(g)` in `strip_scan.hpp`;
    `near_constraint` reads it in place of the literal `1e-10`; L14's radius
    is passed as `coincidence_radius(g)`. About 5 lines. L14 is otherwise
    unchanged and can be built now.

**Settled for 15f-3 after 25.** Increment 25 (#157) replaced the
`elevation_source` sentence with named fields
(`docs/increments/25-plain-output.md`, D2 and D6), and `@tester` amended
15f-3's red suite to match (`86e1074` on `worktree-15f-3`). These are the pins
that amendment relies on, ruled by `@architect`, 2026-10-03. Only P2 changes a
test.

- **P1. Confirmed: on the projected path the at-vertex rows are required.**
  `refine_strip` produces `coincident` and `coincident_max_error` (L14, with
  L16's radius), so `dem_nodes_at_vertices` and
  `dem_nodes_at_vertices_max_error_m` are `--stats` rows there, 0 included
  (25, D6). The file's `max_error_m` is
  `max(refine's max_error, refine_strip's max_error, coincident_max_error)`.
  The first two make the "at most" bound of D7; the third comes from 25's D2.
  So `max_error_m ≥ dem_nodes_at_vertices_max_error_m` holds by
  construction.
- **P2. `line_check_dem_nodes_inserted` is required on the projected path
  and absent on the reprojected path.** It is refine_strip's `nodes_inserted`,
  the work of the DEM rescan (F2), and that rescan does not exist on the
  reprojected path. 25's D6 already rules, for the at-vertex rows, that a
  figure is absent from `--stats` and `--record` on a path that does not
  produce it, "rather than a 0 nobody measured". The same rule applies here.
  `PointRefineOutcome::nodes_inserted` is 0 in `refine_points`, but that 0 is
  a default, not a measurement. **`@tester` changes:** in PY4,
  `assert found.dem_nodes_inserted in (None, 0)` becomes
  `assert found.dem_nodes_inserted is None`. PY3's
  `is not None and >= 0` stands.
- **P3. Confirmed: every other `line_*` row is required on both paths.** The
  strip runs on both: `line_points_checked`, `line_max_error_m`,
  `line_points_on_nodata`, `line_points_refused` and
  `line_points_refused_max_error_m`, `line_points_inserted`, and
  `line_points_duplicate`. Zeros stay in `--stats` (25, D3 rule 3).
- **P4. Confirmed: on the reprojected path the check is
  `at_vertices ≤ max_error_m ≤ max(TOLERANCE, at_vertices)`, not an
  equality.** There `max_error_m = max(the final check's max_error,
  coincident_max_error)` (25, D2). The final check stops only when its own
  figure is at most the tolerance, and that figure is not a `--stats` row, so
  the two-sided bound is the most the test can state without the producer's
  record. Checking `resampled_grid` by its prefix is right. 25's field
  table gives the value's form (`30 m square grid in EPSG:31983, 158 columns
  x 130 rows`), and after the cell size the rest is the CRS and the grid's
  size, which belong to the fixture.
- **P5. Confirmed: the projected path keeps `max_error_m ≤ TOLERANCE`, the
  plain form.** It is safe and it is the stronger test.
  - On that path every vertex a DEM node can lie within `r(g)` of has a z
    from the bilinear surface. Start and strip vertices take `vertex_z`, and
    inserted nodes take their own value. So `coincident_max_error` is at
    most the surface's slope times `r(g)`: about 1e-9 m on these fixtures,
    against `TOLERANCE` = 1 m.
  - The max form is needed only where the tolerance is 0, and
    `test_cli_mesh_refine.py`'s `--tolerance 0` case already uses it (25,
    D6).
  - On the reprojected path the at-vertex difference compares a source z
    with a vertex z and can be any size (15c), which is why P4 takes the max
    form there.
  - If a projected fixture ever fails the plain form, the at-vertex figure
    has outgrown slope × `r(g)`. That is a defect to report, not a reason to
    loosen the test.

Two notes:
- **25 answered L14's question for Ola** (whether the file states the E2
  exception): it does, through `max_error_m`, which includes the nodes at
  vertices, with no extra field (25, D2). L14's first "ASK OLA" line is
  closed by that ruling. The second, pinning `-ffp-contract=off`, stays open
  and outside 15f.
- **Q1 (the every-point form, recursive midpoints) is deferred to 15f-3's
  acceptance.** It is decided with the measured share of midpoints inserted
  on Bygdin and the basin piece. Nothing in 15f-3's design or red suite
  depends on it.

**Settled after 15f-3's green step (7d841f3).** `@developer` left seven
cases in five tests red as disagreeing with the design. Ruled by `@architect`,
2026-10-03. All seven are test-side, so `@tester` changes them in one commit
with the reasons, and `@developer` changes nothing.

- **S1. `test_core_edge_strip.py::TestRefineStrip::test_e1_e2_and_delaunay_by_the_oracles[0.0]`
  and `[0.5]`: fix the fixture.**
  - E2 (DEM nodes kept) holds for a start that is `refine`'s output at the
    same tolerance. That is F2's whole argument: the strip run keeps refine's
    guarantee and repairs only the triangles it writes (L3, E8).
  - The fixture refines `Start` at 2.0 and then runs `refine_strip` at 0.0
    and 0.5. The unwritten triangles hold refine's 2.0 result, so the 100
    and 72 nodes over (worst 1.94 m) are the fixture's, not the run's.
  - The CLI always passes refine's output at the run's own tolerance.
  - E2 above now states this precondition.
  - **`@tester`:** build `Start` at the parametrised tolerance, so that
    `refine` and `refine_strip` both run at 0.0, 0.5 and 2.0. Keep E1, E2
    and Delaunay at each. The control case (refine alone, at 2.0) stays as
    it is.
  - Checking E2 only at 2.0 is rejected, because it would leave E2 untested
    at tolerance 0, where the rescan does the most.
- **S2. `test_cli_mesh_edge_strip.py::test_a_failed_strip_run_writes_nothing`:
  a test bug.** The needle `"planted strip-run failure"` contains spaces,
  but the output is joined with all whitespace removed.
  **`@tester`:** normalise both sides the same way, either
  `" ".join(plain(result.output).split())` or the needle with its spaces
  removed.
- **S3. `test_cli_mesh_plain_output.py::test_max_error_bounds_the_dem_nodes_inside_the_mesh[2]`
  and `[0.5]`: 25's assertion is superseded by P1.** Line 323 asserts
  `dem_nodes_at_vertices` absent on the projected path. That was true before
  15f-3 (25's D6 says "today"), and P1 now requires the rows there.
  **`@tester`:**
  - assert both at-vertex rows present;
  - assert `max_error_m ≥ dem_nodes_at_vertices_max_error_m`;
  - keep `max_error_m ≤ tolerance_m` (P5);
  - update the docstring ("no at-vertex figure before 15f-3").

  This edits a test of increment 25, which has merged. That is right:
  15f-3's behaviour is what changes the answer, and 25's D6 foresaw it.
- **S4. `test_cli_mesh_domain_crs.py::TestTheSameCrs::test_the_mesh_is_16s_bit_for_bit`:
  the test accepts "same up to rounding" for the strip's vertices. No C++
  change.**
  - *Why the meshes differ.* Lattice coordinates are measured from the
    corner of the raster the core is given (L16).
    - The current path cuts the DEM window to the domain. The as-16 path
      opens the whole file. Both are on one lattice, but their origins differ
      by a whole number of cells.
    - Before 15f every inserted vertex was a node, and a node's world
      coordinate `x_min + col·dx` is exact from either origin. So the two
      meshes were bit-identical.
    - A strip point is a computed fraction. Its rounding depends on the
      origin, so 8 of 273 points differ by up to 7.1e-15.
  - *Why not measure from a fixed origin.* It would have to cover the whole
    chain, not only the generator:
    - every off-node position, which today is `(x − x_min)/dx` (16's start
      vertices included);
    - the loop's insertion coordinates;
    - the output's `x_min + col·dx`.

    That is a redesign of the core's frame (15f-1, 15f-2 and 16), for a
    property no user or later increment needs, and recommended against.
  - *Not a defect for 23.* Increment 23 needs bit-identity only where two
    pieces compute the same thing: the seams (K4, the conformity check).
    - The strip makes no points on frozen (seam) edges (23, "What 23
      added").
    - 23b's seam pass runs on "a DEM strip that is a function of the edge
      alone (its bounding box snapped outward to the global lattice and grown
      by one node)" (23, the seam pass). Both neighbours therefore use the
      same origin and get the same bits, whatever their windows.
    - A non-seam constraint edge belongs to one piece, and no other
      computation of it has to agree.
    - **One condition this places on 23b, recorded here for its design:**
      when the seam pass calls `constraint_check_points` (which 23 reuses),
      it passes that edge-determined strip's geometry, never the piece's
      window. 23b's red step tests both sides of a seam bit for bit with
      **different** piece windows.
  - **`@tester`:**
    - keep bit-identity for everything a window cannot change: triangles,
      constraint edges and masks, every start vertex, every vertex that is a
      DEM node, and their z;
    - for the other new vertices (the strip's), assert each coordinate
      within 1e-9 lattice units of its counterpart, and z within
      `1e-9 · max(1, |z|)`;
    - rename the test to say so;
    - keep the quarter-circle variant (`test_the_quarter_circle_is_16s_bit_for_bit`)
      to the same rule if it shows the same difference.

    If connectivity ever differs, a predicate decision flipped on a
    rounding-level difference. That is to be reported, not absorbed.
  - **For Ola (not blocking):** 15b's pin "the same-CRS path gives 16's
    mesh bit for bit" becomes "bit for bit except the edge strip's vertices,
    which agree to rounding". Restoring full bit-identity is the redesign
    above, which is not recommended.
- **S5. `test_cli_mesh_domain.py::TestDomainOutput::test_boundary_z_is_bilinear_and_inserted_vertices_are_nodes`:
  its premise is superseded by E1 and E5.** Its premise, that every
  non-corner vertex is a DEM node, no longer holds: strip insertions are
  crossings and midpoints on the outline and the hole.
  **`@tester`:**
  - a new vertex that is a DEM node keeps the node-value check;
  - every other new vertex must lie on the domain's outline or hole (within
    `ON_INPUT`, 2e-3 m, the noder's 1 mm snap with margin, as in
    `test_cli_mesh_edge_strip.py`) and carry bilinear z, like the
    corners (E6: an inserted strip vertex carries its point's `vertex_z`);
  - rename the test to match.

The citations that moved with `7d841f3` are corrected above: D5's binding
(`bindings/core.cpp:1182`, recomputed when 23c-1 and then increment 29 PR 1 took in master) and `_dem_mesh` lines (`:1466`, `:1550`), `final_check.run` (`:28`), and D6's
`:1561` and `:1580`. The "Review" record of 15f-2's round 2 cites
`cli.py:1417` and others as they were at `17c2d14`. That is history, and it is
left as written.

**Settled after 15f-3's acceptance (9f2f7e5).** `@perf` accepted 15f-3 with a
cost finding (`docs/benchmarks/2026-10-04/15f-3-acceptance.md`, "Where
refine_strip's untimed time goes"): on the projected path, `refine_strip`
spends about 90 % of its time in `detail::to_lattice` rebuilding refine's
output as a `LatticeMesh`. The cost is fixed per call, about 0.8 to 1.4 µs per
triangle of the mesh it receives, and the same with an empty strip. Bygdin at
1 m: +1.5 s (+34 %) process time and +0.29 GiB peak memory. Ruled by
`@architect`, 2026-10-04.

- **A1. 15f-3 is not blocked; the fix is a follow-up PR, 15f-4.**
  - The cost is not new to the code base. On master, the reprojected path's
    final check already calls `to_lattice` on refine's output
    (`refine_points.hpp:248`, reached from `final_check.run`). 15f-3 adds the
    same rebuild to the projected path; on the reprojected path it adds
    nothing of this kind (Velhas: +2.6 to +7.5 % process time, the strip's
    own work).
  - The basin (ANADEM, geographic, meshed in a projected CRS) runs the
    reprojected path. So the per-triangle cost that matters for São
    Francisco is in `to_lattice` itself, and it is there with or without
    15f-3. The fix belongs where both paths get it: in the rebuild.
  - 15f-3 is reviewed (APPROVED, round 2) and accepted on correctness. The
    fix edits `LatticeMesh::build`, which every refine path uses, so its
    acceptance is refine's whole bench (meshes byte-identical), not 15f-3's.
    One concern per PR, and a performance change that can be bisected alone.
  - 15f-4 lands before 23b, whose seam pass rebuilds per tile and would
    multiply the cost.
- **A2. The fix's shape: make the rebuild cheap. Do not fuse the calls, and
  do not pass a mesh handle across the binding.**
  - **`LatticeMesh::build`: no node-based container.** Replace the
    `unordered_map` of directed edges by a flat vertex-bucketed table: count
    each triangle's out-edges per `from` vertex, prefix-sum, fill
    `(to, triangle, slot)` into one `std::vector`, sort each bucket by `to`.
    A duplicate directed edge is two equal `to` values side by side (the
    same refusal as today). The neighbour across `a -> b` is a binary search
    for `a` in `b`'s bucket. Cost is linear plus `O(d log d)` per vertex of
    degree `d`, so a fan of high degree stays cheap (a linear scan per
    lookup would be quadratic in `d`). One allocation per array, freed at
    once.
  - **`to_lattice`'s constraint lookup: no `std::map`.** A sorted
    `std::vector` of `(min, max) -> mask`, searched with `lower_bound`. It
    keeps today's meaning exactly: when an edge is listed twice, the later
    mask wins (`std::map::operator[]` assignment), and entries past
    `masks.size()` are ignored.
  - **`legalise_all` stays** in `refine_strip`. It is 5 to 6 %, it marks the
    triangles it flips as written (L3), and `refine_strip` is a public
    binding that must not assume its start is Delaunay.
  - **Why not a fused C++ call (refine, generate, strip).** It merges three
    steps that `_dem_mesh` composes declaratively; every `--stats` row of the
    three would have to come out of one C++ outcome; it duplicates refine's
    parameter list; and it helps only the projected path, so the
    reprojected path (the basin's) would need a fused twin of its own.
  - **Why not a handle** (refine returns an opaque `LatticeMesh` that Python
    passes to `refine_strip`). It keeps a mutable C++ mesh alive in Python
    between calls, which the binding firewall exists to avoid (Python sees
    arrays; C++ functions are pure over them), and it brings lifetime and
    thread rules to the binding. It also keeps the lattice mesh resident
    beside refine's output arrays through the generator, so peak memory
    does not fall.
  - **When to reopen fusion.** If, after 15f-4, the empty-strip
    `refine_strip` call on Bygdin 1 m (`@perf`'s `prof_strip.py`) still takes
    more than 10 % of the base's process time (0.43 s of 4.32 s), the
    rebuild is still the problem and fusion comes back as a design question.
  - Boundaries are unchanged: no path, no CRS, no new binding.
  - **As built (`5eeb87a`, 21 net lines against about 30).** Two departures
    from the text above, both accepted:
    - The table's entries are `(to, triangle)`, not `(to, triangle, slot)`.
      The slot is not needed: the lookup for triangle `t`'s edge `k` already
      knows `k`, and the neighbour is the triangle found.
    - Every triangle is validated (index range, positive orientation) before
      the table is built, where the old code checked each triangle and
      inserted its edges in one pass. The result is the same: any refusal is
      `nullopt`, and which check fires first is not observable.
- **A3. LOC.** About 30 net (build's adjacency about +20 against the 15 it
  replaces; the constraint table about +8 against 4), 42 at +39 % and 48 at
  +60 %. Folding it into 15f-3 would give about 210, far under 700, so the
  split is for review and acceptance (A1), not for size.
- **A4. Tests first.** The refactor must keep `build`'s and `to_lattice`'s
  behaviour, and some of it is not pinned today (the existing cases are
  "build derives adjacency across the shared edge" and "build refuses a
  clockwise or a zero-area triangle", `test_mesh_lattice_split.cpp`).
  `@tester` adds, against today's code (they pass on master, and must still
  pass after):
  - **B1.** `build` refuses a directed edge used twice (two coincident
    counter-clockwise triangles; three triangles on one edge).
  - **B2.** `build` refuses a vertex index out of range and mismatched array
    lengths.
  - **B3.** `build`'s neighbour table equals a brute-force oracle
    (all pairs of triangles, shared reversed edge) on a grid triangulation
    with shuffled triangle order, boundary edges `kNoNeighbour`, and on a
    fan whose centre has degree at least 1,000.
  - **B4.** `to_lattice` with a constraint edge listed twice with different
    masks gives the later mask, and ignores edges past `masks.size()`.
  - The whole-pipeline identity is `@perf`'s: the bench's tile and quarter
    meshes, Bygdin and Velhas, byte-identical to the base.
- **A5. 15f-4's acceptance (`@perf`).** `tools/bench.py` back to back with
  `--tree`: meshes byte-identical on both domains, refine within noise or
  faster. Bygdin 1 m: process time and peak memory against 15f-3's, and the
  empty-strip call per A2's trigger. Velhas 5 m: the final check's time
  against 15f-3's (it gets the same rebuild).
- **A6. Nothing here needs Ola's ruling.** No output changes and no earlier
  ruling moves. The one choice that is his is when to merge: see Q2 under
  "Questions for Ola".

**Settled after 15f-4's guard tests (2afc2f0).** `@tester`'s B1-B4 pass on
the code as it is (ctest 910/910). Five choices the design left open, ruled
by `@architect`, 2026-10-04:

- **G1. B4 lives in `test_refinement_lattice_position.cpp`, not a new
  target: confirmed.** `to_lattice` is in the header that suite already
  builds (`refine.hpp`), and a target for four cases costs a link for
  nothing. But the file's top comment and the CMake comment above
  `add_terrain_backend_test(test_refinement_lattice_position ...)` still say
  the suite is `lattice_position` only. **`@tester` amends both** to say it
  also pins `to_lattice`'s constraint lookup (15f-4, B4), and tags B4's cases
  `[to_lattice]` instead of `[lattice_position]`, so a filter by tag finds
  them.
- **G2. B3 pins each triangle's slot order and its bits and masks; B4 pins
  that a later mask of 0 still sets the constrained bit, and that masks
  past `edges.size()` are ignored: confirmed.** Each is today's behaviour,
  and each is something A2's rewrite could change silently: refine's output
  depends on triangle order (meshes byte-identical, A5), and a sorted table
  that tested "mask != 0" for presence would drop the bit. The loop's
  `i < edges.size() && i < masks.size()` bounds both arrays, so both
  overhangs are today's meaning.
- **G3. A hand Fisher-Yates over raw `std::mt19937` output: confirmed.**
  `std::mt19937`'s sequence is fixed by the standard; `std::shuffle` and
  `std::uniform_int_distribution` are not, so they would give a different
  shuffle per standard library, and a failure on CI not reproducible on the
  Mac.
- **G4. No NaN or infinite coordinates: confirmed.** A2 does not touch
  `orient_sign` or the orientation check in `build`, and `to_lattice` refuses
  a vertex outside the node rectangle before `build` sees it. Nothing in
  15f-4 moves that behaviour.
- **G5. `n >= kNoNeighbour` not tested: confirmed.** It needs more than
  4 G triangles, which no unit test can allocate. A2 keeps that guard as it
  is; `@reviewer` checks by reading that the line survives the rewrite.

**Python (pytest), 15f-3** (L15; first planned for 15f-2):

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
| **15f-2** | **The strip in the loop, C++ only** (L15) | | | |
| | `refine_points.hpp` and `strip_scan.hpp` (L13): `point_loop`, the strip scan and ownership, the sub-edge map, `strip_fits` (L10) and refused points, consumption, the F2 rescan, L12, the end pass, two entry points, the outcome, `PointScan`; `geometry.hpp`'s `operator==` | 218 | | |
| | the same, **measured** at `61cdbaa` | *275* | | |
| | L14: the radius in `scan` and `scan_points`, the end-pass count | 15 | | |
| | **15f-2 total** (measured 275 + 15) | **290** | **296** | **299** |
| | 15f-2 **measured** at `99e95d7` (L14 and L16 in; `@reviewer` round 1) | *316* | | |
| **15f-3** | **The bindings, the Python and the CLI** (L15) | | | |
| | `bindings/core.cpp`: `ConstraintCheckPoints` and `constraint_check_points` (moved from 15f-1, ruling of 2026-10-03 below) | 30 | | |
| | `bindings/core.cpp`: `refine_points(..., strip)`, `refine_strip`, the outcome properties | 62 | | |
| | `_core.pyi`: all the stubs, the store's and the generator's included | 35 | | |
| | `edge_strip.py` | 30 | | |
| | `final_check.py` | 5 | | |
| | `cli.py`: the calls, the sentence, the report | 25 | | |
| | **15f-3 total** | **187** | **260** | **299** |
| **15f-4** | **The rebuild made cheap** (A1-A3, after 15f-3's acceptance) | | | |
| | `lattice_mesh.hpp`, `build`: flat vertex-bucketed adjacency; `refine.hpp`, `to_lattice`: sorted constraint table | 30 | 42 | 48 |
| | the same, **measured** at `5eeb87a` (`lattice_mesh.hpp` 17, `refine.hpp` 4) | *21* | | |

For 15f-2 the margins apply to the 15 estimated lines only, since the 275 are
measured.

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

*Superseded by L15:* the measured C++ outgrew the two-PR plan, so 15f-2 is the
C++ loop and 15f-3 carries the bindings, the Python and the acceptance.

**Why more than 15c's 160-210.** F1 needs points filed by edge, with the
sub-edge map and the guarded insertion (about 75 lines). F2 needs the DEM
rescan and the second entry point (about 45). Binding two entry points and a
store is about 105. The midpoints themselves cost about 10, as 15c said.

Module to watch: `refine_points.hpp` grows from 294 lines to about 500. If
it passes 550, the strip scan and the sub-edge map move out. It passed (582 at
`67df94c`); L13 rules where they go (`strip_scan.hpp`, not
`constraint_points.hpp`).

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

**15f-2** (L15) adds C++ that no CLI path calls yet, so its acceptance is
15f-1's: `tools/bench.py run` back to back with `--tree` against the previous
merge, with geometry identical and refine time within noise. `refine` itself
is touched only by L14's defaulted radius.

**15f-3** changes every tolerance mesh with a constraint, so:

1. **`tools/bench.py`**, as above. Expected: the tile's geometry identical
   (its outline lies along grid lines: E7 and the grid-line remark), the
   quarter's not (its domain outline is off-grid). The strip's own phases
   are reported beside refine, and the end-to-end time compared. Refine
   itself is not edited, so its time should be within noise.
2. **Bygdin with CORINE**, 16c's script
   (`docs/benchmarks/2026-09-29/bygdin-landcover/run.sh`), at 10 m and 1 m
   tolerance, master and 15f-3 back to back, same power state. Record: strip
   points, `no_data`, inserted (strip and nodes), refused, rounds, the
   strip's time, peak memory, worst angle and share under 10° (bench.py's
   `quality()`). The **independent check**, a script of `@perf`'s that does
   not import `tin_engine`: the crossings of the input outline and the
   CORINE borders with the DTM10 grid lines, z bilinear from the DEM,
   against the written mesh's constraint lines. On master it measures the
   projected-path Surprise 3 for the first time (count over the tolerance,
   and the largest error). On 15f-3 it must find 0 over, or exactly the
   refused points.
3. **A basin piece.** The DEM-derived Velhas test catchment (15c's
   acceptance domain) at 20, 10 and 5 m, and one BHO level-3 unit, 766 (the
   smallest peak in `docs/benchmarks/2026-10-02/basin-level3/README.md`,
   2.82 GB at 20 m), at 20 and 5 m, ANADEM from the cache. 15c's source-node
   check stays at 0 over. The independent strip check, against the resampled
   grid, is 0 over. Record the strip's time and memory beside the final
   check's.

Evidence under `docs/benchmarks/<date>/15f-2-acceptance.md` and its
directory (the no-change run) and `docs/benchmarks/<date>/15f-3-acceptance/` (the rest).

## Which branch to base on

Base the implementation branch on **15e's tip**, `worktree-agent-ab601ec06c6508cb9`
(*default*), unless 15e has merged by then, in which case base it on master.

- 15f-3's Python edits (L15; first planned for 15f-2) land in exactly the lines 15e rewrites: `_dem_mesh`'s
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
  with 15f-3's acceptance numbers (how many midpoints are inserted at all;
  15f-2 changed no mesh, so it had none). **Deferred to 15f-3's acceptance**
  ("Settled for 15f-3 after 25"). If
  Ola wants it, it is best ruled before `@tester` writes ES2, which would then
  sample every piece densely.

- **Q2. When to merge 15f-3** (not blocking; *default*: merge when ready).
  In plain words: "15f-3 makes meshing a Norwegian DEM at 1 m tolerance about
  a third slower end to end (Bygdin: 4.3 s to 5.8 s), because it rebuilds the
  mesh once more. A follow-up of about 30 lines (15f-4) is meant to remove
  most of that cost, for every run. Merge 15f-3 now and the follow-up after, or hold 15f-3
  and merge the two together?"

## Review

**15f-1, round 1, 2026-10-03.** Range `390b516..00d9239` (design 4bfeb47, red 46e69ad, ruling 08e2e14, green 00d9239). Verdict: CHANGES REQUESTED. LOC: 128 net (137 added, 9 deleted) against an estimate of 100 (+28 %, inside +39 %). Blocking: (1) "s strictly increasing" is false at the P1 end: t can round to exactly 1.0 when the end lies within an ulp past a grid line, so the last crossing and its midpoint with P1 coincide (probe: ends at col 0.8932792255671602 and nextafter(10, 11), last two points both col 10, row 3.25, s 1, duplicates 0); at P0 a crossing and its midpoint can share a position with distinct s; D2 steps 5-6 to be ruled, then a test, then the fix. (2) Red-step scaffolding at tests/cpp/CMakeLists.txt@00d9239:327-328. (3) Stale citations in this file at :268, :435, :499 (refine.hpp lines moved by the extraction). (4) project_structure.md and the ROADMAP row, which the design assigns to this PR. Also to record: the .at() bounds checks on edge(k)/on_edge(k), and the duplicate comparison against the previous kept crossing. The lattice_position extraction is behaviour-preserving. Not pushed; no CI.

**15f-1, round 2, 2026-10-03.** Range `390b516..2b3ee11`; this round `00d9239..2b3ee11` (e9583e9 round-1 record, 9b249e1 ruling, a991856 red, 2b3ee11 green). Verdict: APPROVED. LOC: 136 net (145 added: constraint_points.hpp 133, refine.hpp 12; 9 deleted) against an estimate of 100, +36 % net, inside +39 %; the extra lines are the I1-I3 guards, which the design did not cost (expect the same in 15f-2 where the loop meets ulp-scale ends). All four round-1 items closed; both probes rerun against 2b3ee11 (P1: 23 points, duplicates 1, s strictly increasing; P0: 10 points, duplicates 1). At-risk citations re-read: 05b-noder-driver.md:1749 holds; this file's round-1 record is history; tests/python/test_features.py:582 -> project_structure.md:359 was already stale on master (follow-up, not this PR's). Not pushed; no CI.

**15f-2, round 1, 2026-10-03.** Range `9815ea4..fdbf93c` (red afd2498 … test fixes fadcc77, merge fdbf93c). Verdict: CHANGES REQUESTED. LOC: 316 net (379 added, 63 deleted): `refine_points.hpp` 138 net, `strip_scan.hpp` 165, `scan.hpp` 12, `geometry.hpp` 1; estimate about 295 (290 plus L16), +7 %. Code holds against D4 and L1-L16; `refine` and the no-strip `refine_points` unchanged by reading (measured identity is @perf's); fadcc77's ES15/ES16 corrections are fair (ES15 still fails against 3b59934's headers). Default build 897/897, nofma strip suites 47/47, TSan clean locally. Blocking: (1) the L9 TSan entry in `.github/workflows/main.yaml`; (2) red-step scaffolding in `tests/cpp/CMakeLists.txt` and the two test headers ("while the helper is missing", "ES1 to ES11", "see the handback"); (3) doc prose at :1143, :1191 (the (15, 8) behaviour since 99e95d7), the `scan.hpp` citations for `scan` and `vertex_z`, and the status line; (4) the ROADMAP 15f row. Not pushed; no CI.

**15f-2, round 2, 2026-10-03.** Range `beb6505..17c2d14` (4f3e2cb test comments, 5a5486e docs, 17c2d14 TSan list and C++ comments); whole PR `9815ea4..17c2d14`. Verdict: APPROVED. LOC: 316 net (379 added, 63 deleted), unchanged from round 1 (this round's `include/` changes are comments only); about 7 % over the ~295 estimate. All four round-1 blockers closed; no production behaviour change in this range; citations re-read at 17c2d14 (`scan.hpp:96`, `:73`, `refine_points.hpp:101`, `cli.py:1417`, `:1500`, `:1511`, `:1517`, `final_check.py:22`) hold. Not pushed; no CI (the TSan entry's first CI run is pending the push).

**15f-3, code review, round 1, 2026-10-03.** Range `c193cb1..2dda16a` (red 4157dab, merge 8decdbf, amendments 86e1074 and 0a192f0, rulings cfd93b7, green 7d841f3, rulings 2f23a6d, tests 2dda16a). Verdict: CHANGES REQUESTED. LOC: 179 net (222 added) against 187. Code sound: bindings per PY1 (GIL released, read-only store, keyword-only strip, shape checks, refusals), the CLI order per D6, the `--stats`/`--record` rows per P1-P3, `max_error_m` per P1/P4; S4 still fails on real changes; @tester's three pins accepted. 901/901 ctest, 4062 pytest. Blocking: (1) present-tense red-step text in three test files; (2) this file's status line and ROADMAP row 15f say 15f-3 not started (and 15f-2 unmerged); (3) `test_edge_strip.py:174` cites `final_check.py:22` (now :28); (4) `project_structure.md` lacks `edge_strip.py`. Suggestions: a unit test for the refused-points warning; `refined_record`'s docstring on `max_error_m`; a docstring sentence on S4's start-vertex prefix. @perf's acceptance outstanding. Not pushed; no CI.

**15f-3, code review, round 2, 2026-10-03.** Range `ae6bb41..9265916` (e6348b4 tests, 9265916 docs); whole branch `c193cb1..9265916`. Verdict: APPROVED. LOC: 179 net (222 added), unchanged; this round touches no production file. All four round-1 blockers closed; the red-step figures in `test_cli_mesh_edge_strip.py` ("98 points, worst 16 m … 152, worst 47 m") are @tester's record of the `4157dab` run, not rerun. Five touched test files 200 passed; ruff clean; merge-tree against origin/master `1a422b5` clean. Suggestion: fix `tests/python/test_features.py:583`'s citation of `project_structure.md:359` (stale on master since 15f-1). @perf's acceptance (Bygdin, basin piece) outstanding. Not pushed; no CI.

**15f-3, code review, round 3, 2026-10-04.** Range `60e2b45..ac06d39` (83c7fd2 citation, 588e879 and 9f2f7e5 @perf acceptance and profile, 107cdb4 ruling A1-A6, ac06d39 ROADMAP); whole branch `c193cb1..ac06d39`. Verdict: APPROVED. LOC: 179 net (222 added), unchanged; no production file in this range. The ruling's claims hold against code (`include/terrain/refinement/refine_points.hpp@ac06d39:220` calls `to_lattice` on refine's output inside `detail::point_loop`, reached from `final_check.run`; `include/terrain/mesh/lattice_mesh.hpp@ac06d39:131`'s `unordered_map`; `include/terrain/refinement/refine.hpp@ac06d39:154-156`'s `std::map`) and every figure matches the acceptance file; status line and ROADMAP row match the tree; merge-tree against origin/master `2060f14` clean. Not pushed; no CI.

**15f-4, code review, round 1, 2026-10-04.** Range `107cdb4..41a5ad7` (guards 2afc2f0 and 6c5fb31, rulings G1-G5 43c256f, green 5eeb87a, @perf acceptance 41a5ad7). Verdict: CHANGES REQUESTED. LOC: 21 net (41 added; `lattice_mesh.hpp` 17, `refine.hpp` 4) against about 30. Code correct by reading: validation before the table, duplicate refusal by adjacent equal `to`, `lower_bound` neighbour with `kNoNeighbour` default, the guard `n >= kNoNeighbour` kept (`lattice_mesh.hpp:129`, G5); `(to, triangle)` suffices since the slot is the lookup's k; `to_lattice`'s three rules hold (`stable_sort` plus last of `upper_bound`; match by key; loop bound unchanged); ctest 910/910 (not rebuilt); `25-plain-output.md:115`, `:120` follow the moved lines. Blocking: (1) `tests/cpp/unit/test_mesh_lattice_split.cpp@41a5ad7:355-358` and `tests/cpp/unit/test_refinement_lattice_position.cpp@41a5ad7:118-121` still describe the maps as today's; (2) the status line calls 15f-4 a future follow-up; no as-built note (departure from A2's `(to, triangle, slot)`, the validation order, 21 against 30 lines); ROADMAP row 15f says 15f-4 planned. Suggestion: the `build` comment states only why the table is flat, not which callers exist. Not pushed; no CI.

**15f-4, code review, round 2, 2026-10-04.** Range `36b4005..3c72c8b` (0ce1a1a test banners, 4169e4e status line, as-built note and ROADMAP, 3c72c8b `build`'s comment); whole PR `107cdb4..3c72c8b`. Verdict: APPROVED. LOC: 21 net (41 added), unchanged; comment-only code change this round. Both blockers and the suggestion closed; the banners match B1-B4's sections; the as-built note's two departures hold by reading; the guard stays at `lattice_mesh.hpp:129`; citations hold; merge-tree against origin/master clean. Not pushed; no CI.

**15f-4, code review, round 3, 2026-10-04.** Range `3c72c8b..0e58ae6` (d84315c round 2 recorded; 52a56fe merge of origin/master 435aa56, #161; 0e58ae6 merge of origin/master 4f56551, #165-#168 and #162 23b, citations recomputed or pinned). Verdict: CHANGES REQUESTED, closed in the commit that records this round. LOC: 21 net (41 added, 20 removed), unchanged; the merges add none. Outside `docs/` the merged head is master plus exactly 15f-4's diff (same six files, identical added and removed lines against `435aa56..52a56fe`); `lattice_mesh.hpp` and `refine.hpp` changed on both sides in disjoint hunks, and `build` and `to_lattice` do not read `frozen_`. The recomputed citations (`refine.hpp` +8, `lattice_mesh.hpp:129`) quote the same text as on master, and the pins `@ac06d39` and `@00d9239` hold. `check_citations.py --base origin/master` passes. ctest 983/983 on `build-15f-4` (not rebuilt); gates clean. @perf need not re-time: no refine or mesh hunk was resolved by hand, and 15f-4's figures are stated against 15f-3 and the base. Blocking: (1) the status line said 15f-3 is on its branch and 15f-4 not pushed (#161 merged; #163 open); (2) ROADMAP row 15f said the same; (3) round 1's `test_mesh_lattice_split.cpp:355-358` and `test_refinement_lattice_position.cpp:118-121` resolved to the rewritten banners, now pinned `@41a5ad7`. CI on #163 is green at 52a56fe only; 0e58ae6 is not pushed and has no CI.
