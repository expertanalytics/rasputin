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

Floating-point segment-segment intersection is numerically fragile, so all intersection points and any near-coincident input vertices are snapped to a precision grid. The natural choice for terrain is the raster pixel resolution (or a sub-pixel multiple of it): the downstream sampling can't resolve finer features, and integer coordinates on the snap grid make robust predicates trivial.

Reference: Goodrich-Guibas-Hershberger-Tanenbaum for classical snap rounding; Halperin & Packer's iterated snap rounding handles the cascading-near-intersection case.

### Algorithm

1. **Broad phase via the raster grid.** Bucket each input segment into the raster cells it touches; only test segment pairs that share a cell.
2. **Pairwise robust intersection** with Shewchuk-style adaptive predicates.
3. **Snap intersection points and near-vertices** to the precision grid.
4. **Split each input segment** at all intersections lying on it, in arc-order along the segment.
5. **Deduplicate** endpoints — multiple constraints meeting at the same snapped grid cell collapse to a single node.
6. **Output** a list of `(p0, p1, is_river)` segments plus the deduplicated vertex array.

### Edge metadata

Only one constraint type carries through to the mesh: **river**. The rest of the hydrology and land-cover semantics is face-based (queried per cell via point-in-polygon against the original input polygons), so non-river constraint edges don't need to remember what feature they came from.

- One bit per edge: `is_river`.
- Merge rule at intersections and coincidences: logical OR.
- If a road runs along a river after snapping, the merged edge keeps `is_river = true` and the road is geometrically forgotten — it added no structure the river wasn't already enforcing.
- Land-cover and other feature polygons are still fed into the CDT (so face boundaries align with feature boundaries) but their edges carry no metadata after noding.

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

CGAL is replaced. Boost.Geometry remains useful for the upstream vector simplification step.

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
