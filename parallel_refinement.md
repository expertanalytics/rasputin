# Parallel Local Refinement of Constrained Terrain TINs

Design sketch for the post-CGAL meshing backend in rasputin.

## Pipeline

1. **Feature extraction.** Compute natural catchments from the DEM. Segment rivers, lakes, roads, and land-cover boundaries from other data sources.
2. **Vector simplification.** Simplify the resulting polygons and polylines to a target resolution before any meshing happens (e.g. Visvalingam-Whyatt for area-preserving polygon simplification, Douglas-Peucker for polylines). Topology preservation across features matters here — simplified polygons must not cross each other or the constraint lines.
3. **Constraint noding.** Detect intersections between input constraints, snap-round to a precision grid, split segments at intersections so the result is a conforming planar straight-line graph with no interior crossings.
4. **Initial CDT.** Feed the noded constraints into a one-shot constrained Delaunay triangulation. This is the only place a CDT library is needed.
5. **Adaptive refinement.** Refine triangles until each one approximates the underlying DEM within a per-triangle error tolerance. No simplification after this point.
6. **Final Delaunay-flip pass.** Improve aspect ratios via Lawson edge flips. Constraint edges from step 3 are excluded from flipping; the original constraint topology is preserved end-to-end.

## Constraint noding (PSLG construction)

Input constraints (polygons and polylines) can intersect each other — a road through a forest, a river crossing a land-cover boundary. The CDT requires a non-self-intersecting input, so before triangulation we build a conforming planar straight-line graph (PSLG) by detecting intersections and splitting segments at every crossing.

### Snap rounding

A fixed-spacing integer lattice in the projected CRS. Every input vertex and,
more importantly, every *constructed* intersection point is rounded to the
nearest lattice point.

**Why it exists.** Exact predicates are cheap; exact *constructions* are not.
`orient2d` decides a sign exactly, but the intersection of two segments is
generally irrational and has to be rounded to land in a `double`. The predicates
module settled this deliberately — signs only, no exact constructions — so the
snap grid is what makes that rounding a single deliberate target rather than
whatever the FPU produced.

What snapping buys is a deterministic rounding target and **exact equality for
coincidence**: two vertices that snap to the same index pair are bit-identical
doubles, so dedup has a legal key. It does **not** buy exact arithmetic in the
predicates, and it is worth being precise about that, because the opposite is
easy to assume. Snapped world coordinates are `index * spacing`, and at the
spacings recommended below neither factor is a dyadic rational, so the products
are not integral. Even if they were, `orient2d` forms differences and then
products: exactness needs coordinate magnitudes under roughly 2^26, not 2^53,
and a large projected CRS is far past that.

Snapping in fact *manufactures* the hard cases: the filter falls through and the
exact backend decides on a large fraction of snapped breakline data — 38 % of
grid-collinear triples at 0.1 m. **What it manufactures is fall-through with a
definite sign, not an unresolvable collinear triple.** An earlier version of
this paragraph said the latter, quoting `docs/increments/01-predicates.md`,
which said it too; both are corrected. Because `world(g) = g·spacing` is not
affine at a non-dyadic spacing, three points exactly collinear on the grid are
usually *not* collinear in world coordinates — at 0.1 m only 207 of 2976
general-direction triples come back `Collinear`. The design absorbs this because
`DefaultKernel` is total and decisive, not because the fallback is rare.
Measured in `docs/increments/kernel-sufficiency-audit.md` §5.1.

**A dyadic spacing buys two different properties and only one of them needs
local coordinates.** An earlier version bundled them and dismissed both.

1. **`world` is an exact affine scaling, so grid-collinearity survives into
   world coordinates.** Needs a **negative power of two alone** — 2^-4 m =
   6.25 cm sits inside the range below. `ix · 2^-4` is exact for `|ix| < 2^53`,
   and Web Mercator at that spacing needs only `|ix| ≈ 3.2e8 ≈ 2^28`. Measured
   at CRS-origin anchoring with no local coordinates anywhere: 3007/3007 and
   3016/3016 grid-collinear triples report `Collinear`. **This one is free, and
   it is the property the noder's topology depends on** — `on_segment`, the
   collinear arm of `classify`, the cell corners of the hot-pixel predicate, and
   collinear overlap being detectable at all. Prefer a dyadic spacing.
2. **The determinants are exact in plain double, so the filter passes.** Needs
   coordinate magnitudes under roughly 2^26, hence local coordinates, and that
   *is* in tension with anchoring at the CRS origin. It stays dismissed: it is a
   performance property only, and its absence costs nothing in correctness — at
   2^-4 m the filter falls through 100 % of the time and answers `Collinear`
   100 % of the time.

Property 1 is guidance on choosing the spacing, not a default. Spacing remains a
declared parameter entering at the Python boundary, and no value is safe for all
input (`docs/increments/05-noder.md` risk 2).

**Spacing is not the pixel resolution.** An earlier version of this section said
the natural choice is the raster pixel resolution, on the grounds that
downstream sampling cannot resolve finer features. That conflates two different
resolutions. Sampling resolution is a fact about **elevation**: a 10 m DEM
genuinely cannot tell you z at finer than 10 m. The snap grid is a
**planimetric** decision, and a 10 m DEM says nothing about whether a lake shore
or catchment outline is known to 1 m or to 10 m — that is a property of the
vector input.

Size it for what it is for: fine enough that snapping does not visibly move
input geometry, coarse enough to collapse genuine near-coincidences and keep
coordinates integral. That is a function of input precision and the tolerance
acceptable on a domain boundary — decimetres or centimetres for typical terrain
work — and is not tied to cell size.

**The grid is not raster-aligned.** The bounding polygon is arbitrary and has no
relationship to the DEM lattice, and nothing downstream wants one. Vertex
planimetric positions are arbitrary, and z is obtained by sampling the raster at
whatever position a vertex has (`include/terrain/raster/sample.hpp`) — never by
moving a vertex to where a sample is. That is the whole requirement, and
off-node vertices are the expected path rather than a fallback.

Note the displacement argument belongs to *spacing*, not to the anchor: at
spacing equal to the cell size, snapping moves a domain vertex by up to half a
cell diagonal — about 7 m on a 10 m DEM — whatever the anchor is. Once spacing
is decoupled from cell size, the anchor on its own costs nothing extra; it
shifts the lattice by a sub-spacing offset and the worst case is half a
*snap*-cell diagonal either way.

**Anchor at the origin of the projected CRS.** Not for simplicity — because it
is the only anchor that makes snapping a *pure function of one coordinate*.

Anchor to the domain bounding box and snapping becomes a function of the whole
input: adding one vertex that extends the bbox shifts the lattice and re-snaps
every other vertex, which can change the topology of a mesh already built. It
also forces a global bbox reduction before any vertex can be snapped, which
serialises the front of a module whose design is four parallel-fors. The raster
origin fails the same way for a different reason — it is an input-dependent
quantity, so re-windowing the DEM re-snaps everything.

Anchored at zero, `snap(p)` depends on `p` and the spacing alone. That is what
makes it reproducible across runs and safe to run concurrently.

**Index type: `int64`.** At the spacings above, indices over a large projected
CRS reach 1e8–1e9, which does fit `int32` — but with under 2× headroom, and
intermediate arithmetic during snapping and comparison overflows it immediately.
A finer grid must not force a type change. Note also that "coordinates under
1e7 m" is an order of magnitude and not a bound: Web Mercator reaches
±2.0037e7, and `project_structure.md` admits any projected CRS, including
foot-based ones. Dedup keys on
`GridPoint{std::int64_t ix, iy}` with exact equality, never on
`std::hash<Point2>`, which cannot find a NaN and hashes near-coincident points
apart.

**Where the spacing lives.** Under the old framing it was derivable from
`RasterGeometry::delta_x`. It is now a declared parameter — of input precision
and boundary tolerance — so it enters at the Python boundary as configuration
and must travel with the data. `docs/increments/03-pslg.md` has `NodedPslg`
promising that "coordinates lie on the snap grid", which is uncheckable unless
the type carries the grid. Increment 5 should give it a `SnapGrid{double
spacing}`: one scalar, because anchoring at zero means there is no offset to
carry. That is the concrete payoff of the anchoring decision.

Reference: Goodrich-Guibas-Hershberger-Tanenbaum for classical snap rounding; Halperin & Packer's iterated snap rounding handles the cascading-near-intersection case.

### Algorithm

1. **Broad phase via the raster grid.** Bucket each input segment into the raster cells it touches; only test segment pairs that share a cell. This is a spatial index for candidate pairs and nothing more — it is not the snap grid, and the two are independent. Any uniform bucketing works; the raster's is convenient because its extent and spacing are already to hand.
2. **Pairwise robust intersection** with Shewchuk-style adaptive predicates.
3. **Snap intersection points and near-vertices** to the precision grid.
4. **Split each input segment** at every node whose **cell** the segment meets, in arc-order along the segment. Not "at all intersections lying on it": "lying on" is exact incidence, and snap rounding is defined on hot-pixel *proximity* — a segment is routed through every cell it passes through, not only through cells whose centre it exactly contains. At a non-dyadic spacing the two differ on the great majority of real input, so a design that splits on exact incidence reports no T-junction and never splits the host. Step 1's "pairwise robust intersection" stays as it is; that is where crossings come from.
5. **Deduplicate** endpoints — multiple constraints meeting at the same snapped grid cell collapse to a single node — and deduplicate **edges** on their sorted node-id pairs, merging the edge property sets by union.
6. **Output** a `NodedPslg`: the deduplicated node array, the chains with their roles and windings preserved, the per-edge property array and the snap grid. Not a flat list of `(p0, p1, is_river)` segments plus a vertex array — that cannot drive the CDT, which needs `addOutline`/`addHole` per ring with a winding (`docs/increments/04-cdt.md`), and a flat edge list has thrown the outer/hole distinction away.

### Edge metadata

**An output edge carries a set of properties, not a bit.** The one-bit `is_river`
form that stood here rested on a claim that is false for every linear feature
that is not a river, and the correction is the user's: *at a coarse resolution
the same line segment can be both a road and a river.*

The false claim was that non-river constraint edges "don't need to remember what
feature they came from", because the rest of the semantics is face-based and
queried per cell by point-in-polygon. **That works for area features and cannot
work for linear ones.** Point-in-polygon can tell you a face is forest, because
forest has an interior to be inside; it cannot tell you an edge is a road,
because a road has no interior. Buffering one into a polygon is a tolerance, and
this project has none. So the face-based fallback recovers land cover and
recovers nothing about roads, railways, walls or contours — and
`parallel_refinement.md`'s own step 1 lists roads as an input feature class
alongside rivers and lakes.

- **A property set per output edge**, merged at intersections and coincidences by
  **union**. Union, not logical OR of one bit: the operation is the same `|` on
  the representation, but what it ranges over is a set, and an edge may end up
  carrying two, three or no properties.
- Union is commutative, associative and idempotent, which is what lets the merge
  be a parallel reduce over an unordered set of contributing chains, in any
  order, with no tie-break. That is the same property node ids rest on and it is
  why the merge is not, for example, "the highest-priority contributor wins".
- **If a road runs along a river after snapping, the merged edge is both.** The
  road is *geometrically* forgotten — it added no structure the river was not
  already enforcing, and that half of the old text was right and is kept — but it
  is not forgotten *semantically*. Those are different objects, and conflating
  them is what the one-bit form did.
- Land-cover and other area feature polygons are still fed into the CDT so that
  face boundaries align with feature boundaries, and their **area** semantics
  stay face-based and queried per cell. That half of the original claim is sound
  and is the reason the edge property set is a small vocabulary of *linear*
  features rather than a general attribute channel: the legacy land-cover
  vocabulary alone runs to 23 classes
  (`legacy/rasputin/globcov_repository.py:16-38`), and none of them belongs on an
  edge.
- **The C++ core does not name the properties.** It merges an opaque fixed-width
  bitset (`include/terrain/core/edge_properties.hpp`, 32 bits); the mapping from
  bit position to feature name is Python's, in
  `src_python/tin_engine/features.py`. Ruled and shipped in
  `docs/increments/07-edge-properties.md`, which 5b consumes rather than
  introduces.

### Watch list

- **Self-intersecting polygons** from upstream simplification — the same machinery detects them; flag for resimplification or auto-fix.
- **Three-way+ intersections** at the snap grid — handled implicitly by deduplication.
- **T-junctions** where a polyline endpoint lands on another segment's interior — snap-round and split the host segment.

### Parallelism

- Parallel-for over raster cells for broad phase.
- Parallel-for over candidate pairs for exact tests.
- Parallel-for over segments for splitting.
- Final deduplication: parallel sort + unique.

### Reference implementation to study

**S2Builder** in `s2geometry` (Apache 2.0) — Google's planar/spherical PSLG builder, does snap rounding at scale. Useful reference even if not used as a dependency.

## Refinement algorithm

Embarrassingly parallel, local-only, no hanging nodes by construction.

**Per round:**

```
active = current leaf triangles
parallel-for triangle T in active:
    err_max = 0
    p_max   = null
    for sample s in raster_bbox(T):
        if not inside(T, s): continue
        e = | s.z - plane(T)(s.x, s.y) |
        if e > err_max:
            err_max = e
            p_max   = s
    if err_max > tol:
        mark(T, p_max)
compact marked triangles
for each (T, p) in marked:
    replace T in active with three children (fan from p to T's three vertices)
repeat until no marks
```

**Why this works:**

- **No hanging nodes.** Inserting a point strictly inside `T` and fanning to its three vertices leaves the boundary of `T` untouched. No neighbor sees a new vertex on a shared edge, so the mesh stays conforming with zero cross-triangle communication.
- **Max-error insertion** (rather than centroid) chooses the point that reduces approximation error most — Garland-Heckbert-quality decisions without giving up locality. The scan that decides "refine or not" already touches every sample, so tracking the max costs nothing extra.
- **Cache-aligned inner loop.** Iterating raster samples inside a triangle's bbox in row-major order is sequential reads from the DEM — friendly to both CPU prefetchers and GPU coalesced loads.
- **Convergence.** For any bounded continuous height field, the planar approximation error within a triangle shrinks as the triangle does. A finite-depth refinement reaches tolerance everywhere.

## Data structure

Ternary tree per root triangle from the initial CDT.

```
Triangle {
    v0, v1, v2 : vertex index
    children[3]: triangle index, or NULL if leaf
}
```

- **Vertex array:** append-only, atomic counter for new-vertex insertion.
- **Triangle array:** append-only; children written contiguously per parent.
- **Active list:** rebuilt each round by compacting marked triangles' children.
- **Constraint flag** on edges (or as a per-triangle bitmask of which of its three edges are constrained) — propagated to children for those edges that coincide with the parent's constrained edges.

Final mesh: flatten leaves to a triangle/vertex array.

## Parallelism

- **CPU:** `parallel-for` over the active list with OpenMP or a work-stealing pool. Atomic compaction at round end.
- **GPU:** per-triangle CUDA/HIP kernel for the scan; thrust/CUB compaction between rounds.
- Per-triangle work is bounded by the number of raster samples in its bbox — naturally load-balances as refinement deepens (small triangles touch fewer samples).
- No locks during the scan phase. Only contention is the atomic vertex/triangle counters at insertion time, which is one write per refined triangle per round.

## Library choices

The CDT library is used exactly once, at step 4 — noding is step 3. Everything after is owned in-tree.

Candidates (MIT / BSD only):

- **Detria** (MIT, header-only, C++20) — top pick. Minimal dependency, modern API.
- **poly2tri** (BSD-2) — polygon-with-holes CDT only, with no per-edge constraint
  entry point, so breaklines and therefore rivers cannot be expressed. Not viable
  under this design; see `docs/increments/04-cdt.md`.
- **Geogram** (BSD-3) — overkill under this design since we don't need its remeshing or spatial-search infrastructure.

CGAL is replaced. Vector simplification is step 2 and belongs to Python, where
Shapely covers Douglas-Peucker. Visvalingam-Whyatt has no Shapely equivalent and
would be written in-tree. Nothing in the C++ core simplifies anything.

## Final flip pass

After the refinement loop converges:

```
parallel-for edge e in mesh:
    if constrained(e): continue
    if not locally_delaunay(e): mark
resolve conflicts (two adjacent marked edges can't both flip in the same round)
flip marked edges
repeat until no marks
```

Constraint edges are skipped, so the original polygon and polyline geometry is preserved exactly. Interior-only flipping is local (each flip touches one edge and two triangles) and parallelizable with a small amount of conflict resolution (e.g. graph-coloring the edge-conflict graph or randomized retry).

## Tradeoffs vs. serial Garland-Heckbert

- **Triangle count** to hit a given tolerance: comparable, since we're using max-error insertion.
- **Wall-clock time:** much better at scale; throughput scales with available cores / SMs.
- **Memory:** ternary tree adds a small overhead per internal node (3 child indices) but is trivially flattened at the end.
- **Quality:** initial-CDT topology of un-flipped *constraint* edges persists, so input feature simplification quality matters. Interior quality is recovered by the final flip pass.

## Open question: how is triangle size controlled where terrain is flat?

Raised 2026-09-23. **Not settled.**

Refinement stops on elevation error alone. A flat body therefore stays coarse,
which is usually right — it is the whole reason a TIN beats a regular grid,
since vertex density follows terrain gradient and orographic precipitation
follows terrain gradient too, so the adaptation transfers. Spending cells on
flat ground is the thing this design exists to avoid.

The exception: **a flat body surrounded by steep terrain may need resolution
for water routing**, even though its own elevation residual is
zero. Water collects there. One huge triangle cannot represent where it goes.

Three ways to express that, increasing in what they assume:

1. **A global area cap.** Simplest, and wrong for the reason above: it spends
   cells on isolated flat ground.

2. **Grading** — bound how much adjacent cells may differ in size. Gives the
   wanted behaviour: a lone plain has flat neighbours and stays coarse, a valley
   floor beside refined slopes is pulled finer. **But it costs this design its
   core property.** The ternary tree stores `v0, v1, v2` and `children[3]` and
   no adjacency at all, so grading needs a new structure plus a cross-triangle
   read every round — and "zero cross-triangle communication" is what the
   algorithm is sold on. It also cascades: refining A forces B forces C, a
   propagation needing iteration to converge.

3. **A sizing field sampled on the raster.** Precompute a target-edge-length
   raster from local terrain roughness, **blur it**, and have refinement read it
   exactly as it already reads elevation — one more lookup in a loop already
   scanning the bounding box. The blur is what produces the grading: a valley
   floor inherits a small target from its neighbourhood because the field was
   smoothed across the boundary, not because a triangle asked its neighbour
   anything.

**Recommended: 3.** No adjacency, no cross-triangle reads, no cascade, and the
field is fixed before refinement starts. It also composes — if flow
accumulation exists it becomes another term in the same field, and refinement
neither knows nor cares which inputs built it.

### Why not drive this from flow accumulation directly

Rejected as a *requirement*, not as a mechanism. `auto_catchments.md` already
plans an accumulation raster for catchment delineation, so reusing it would be
nearly free — but **the tool must work when a catchment polygon is supplied
rather than derived**, and then no accumulation exists and computing one would
impose a cost the caller did not ask for. So it cannot be a precondition of
meshing. As an
optional term in the sizing field it remains available and is worth revisiting.

**Noted for later: flow accumulation as a refinement driver.** Where routing
needs cells along flow paths through a large flat interior, a blurred roughness
field will not put them there — grading of any kind only propagates inward from
the edges. Accumulation would. Revisit when the derive-catchments path exists.

### Separation of concerns, which this must not erode

Constrained vertices are fixed before meshing and cannot be removed by
refinement or flipping. So the count of constrained degrees of freedom is
decided upstream, by whatever produced the polygons and polylines — and a
DEM-derived catchment boundary can carry one vertex per pixel edge.

Each stage does one thing and hands on an artefact:

- **auto-catchment** produces a catchment polygon. Nothing else.
- **simplify** reduces its vertex count to the target resolution. Separately,
  and identically, for every other constraining polyline and polygon.
- **mesh** consumes constraints it does not question.

This is step 2 and step 5's "No simplification after this point" already, and
it is restated because the pressure to fold simplification into meshing will
come from whoever finds a catchment with 40 000 vertices. The answer is a better
simplify step, not a mesher that edits its own constraints.

## Open question: does the flip pass leave the tolerance undefined?

Raised 2026-09-23. **Not settled. This section asks a question and does not
answer it.** Whoever writes the refinement increment must rule on it before
`@tester` is briefed, because the answer decides what the suite can assert.

### The question

The refinement loop terminates when every leaf triangle satisfies
`max |s.z - plane(T)(s)| <= tol` over the samples inside it. That guarantee is
about **each triangle's own plane**.

The flip pass then runs. A flip replaces two triangles over a quadrilateral
with two different triangles over the same four vertices. The vertices do not
move, but the two surfaces agree only along the new diagonal — everywhere else
in the quad they differ. Every sample in that quad is now measured against a
plane that did not exist when the tolerance was checked, and nothing checks it
again.

So the guarantee the algorithm delivers is "error was within `tol` at the
moment refinement converged", and the mesh handed to the caller is not that
mesh.

### How large the disagreement can be

Measured on a saddle — four corners of a unit square with alternating heights
0, 1, 0, 1:

```
same four vertices, centre of the quad
  surface with diagonal a-b : z = 0.000
  surface with diagonal c-d : z = 1.000
  the flip moves the surface by 1.000 at that point
```

The full height range of the data, at a single flip, with no vertex moved. A
saddle is the worst case rather than a typical one, and on smooth terrain the
disagreement is bounded by local curvature — but it is not bounded by `tol`,
and nothing in the algorithm bounds it.

### Why this may be worse than a bookkeeping problem

The flip criterion is the circumcircle test — it moves the mesh toward the
Delaunay triangulation. **Delaunay is not optimal for piecewise linear
approximation of a surface.** Dyn, Levin and Rippa showed that triangulations
chosen from the data values outperform Delaunay for exactly this problem (*Data
Dependent Triangulations for Piecewise Linear Interpolation*, IMA Journal of
Numerical Analysis 10(1), 1990).

If that holds here, the flip pass is not neutral with respect to `tol`: it
systematically spends vertical accuracy to buy triangle shape. This document
currently describes it as recovering quality, which is true for shape and may
be false for the thing the tolerance measures.

Note also that the minimum-error triangulation problem is NP-hard and not
approximable within any multiplicative factor unless P = NP, so "flip by error"
is a heuristic, not an optimisation with a known answer.

### One thing the pseudocode omits, and why that is safe

The flip pass tests the circumcircle and does not test whether the quad is
convex, though a flip of a non-convex quad would produce overlapping triangles.
That omission is safe, and the reason is worth writing down because a reader
implementing it will ask.

An edge that fails the circumcircle test always has a convex quad. Measured over
200 000 random configurations: of 30 659 non-convex quads, **zero** were flagged
as not locally Delaunay. So the circumcircle test already excludes every case
where a flip would be illegal, and a separate convexity test would never fire.

Option 5 below does not inherit this. A patch larger than two triangles has no
such guarantee, so patch growth must test the boundary it is building.

### What a ruling has to choose between

1. **Flip only when every sample in the quad stays within `tol` afterwards.**
   Keeps the guarantee; the mesh is less Delaunay than it could be.
2. **Alternate refine and flip until both hold.** Keeps both properties; needs
   an argument that it terminates.
3. **Keep the flip but change its criterion from the circumcircle test to
   approximation error.** Data-dependent, per the reference above; abandons the
   shape guarantee the Delaunay criterion gives.
4. **Accept it and say so.** `tol` becomes a refinement parameter rather than a
   property of the delivered mesh, stated plainly in the API.

5. **Re-triangulate a patch rather than flipping an edge.** Grow a region of
   triangles bounded entirely by unconstrained edges, take the vertices inside
   it, and re-solve that patch — choosing a triangulation that satisfies `tol`,
   inserting points until one does, or both.

   This is strictly more general than options 1 to 3: an edge flip is the
   two-triangle case of it. It keeps the property that makes the flip pass
   parallel — the patch boundary is fixed, so nothing outside changes, and the
   conflict rule is the same one flipping already needs, that patches may not
   overlap. If the patch polygon is convex, any triangulation of its interior
   points is valid and covers it exactly, so the choice is free rather than
   searched under constraints.

   The reason it answers this section's question, where flipping does not: the
   decision stops being "which of two diagonals" and becomes "any triangulation
   of these points that meets the tolerance". `tol` becomes a constraint the
   step solves under rather than a property it happens to preserve.

   Three consequences a design must state:

   - **It cannot live in the ternary tree.** An arbitrary re-triangulation of a
     patch is not a fan subdivision, so this runs after the flatten step, on a
     general mesh. The pipeline already flattens last, so the ordering works —
     but "Data structure" above says the tree *is* the mesh, and it would stop
     being true here.
   - **Patch growth needs a stopping rule**, and it is the design parameter.
     Too small is flipping; too large re-triangulates the mesh and loses the
     locality that made any of this parallel.
   - **Patches are bounded by constraints**, so a region dense with rivers and
     roads gets small patches. That is likely the right behaviour: those are
     the places where the input geometry should dominate the triangulation.

   This has a name. *Higher-order Delaunay triangulation* relaxes the
   circumcircle criterion over a neighbourhood rather than a single edge, and
   is the formalisation of this idea; see "Implementing data-dependent
   triangulations with higher order Delaunay triangulations", ACM SIGSPATIAL
   2016. **Nobody here has read it.** Read it before designing this.

Nothing here rules 4 out. It rules out leaving the document as it is, which
promises a bound the pipeline does not deliver.

Option 5 came from the user and is materially better than 1 to 4, which were
written by an agent. Recorded with that provenance because this project has
already found one rule that acquired the authority of a human decision without
having been one (`project_structure.md`'s CRS paragraph, corrected 2026-09-23).

### One thing to check before designing any of this

This document predates every increment record and has never been through the
protocol. Its "Library choices" section still discusses selecting a CDT, which
increment 4 settled. Read it against the shipped tree before trusting any of
it — the corner-graze work is the precedent for what an unchecked old design
document costs.
