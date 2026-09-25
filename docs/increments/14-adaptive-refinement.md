# Increment 14 — adaptive refinement against the DEM, to a sup-norm tolerance

Status: **implemented, in review.** Red `65c163f`; green `1cd8438`, `13c10ec`,
`03642ea`, `9b7cdd8`, `7840ed6`, `8869d5e` (test), `f3e8e18`. Measured 668
non-blank production lines (reviewer's count) against ~570, under the 700
ceiling, so the 14b cut did not fire. Real tile at `--tolerance 1`: 670 554
triangles, achieved max error 0.99996 m, about 2 s in the test and 3 s through
the CLI on a loaded laptop. Written by `@architect` before `@tester`,
per `docs/increments/README.md` step 1. The user chose the recommendation on
all three choices on 2026-09-25: U1 (a) opt-in, U2 (a) no flat-ground size cap
yet, U3 (a) no flip pass (section "Ruled by the user").

**Closes.** `ROADMAP.md` gap 2, first half: the mesh follows the terrain instead
of a regular stride. After this increment,

```sh
rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --out tile.vtk
```

writes a `.vtk` whose triangles are small where the terrain bends and large where
it is planar, and whose `elevation_source` field states the tolerance asked for
and the maximum error achieved. Running the same command without `--tolerance`
writes increment 12's uniform mesh (if U1 goes as recommended), so the two can be
opened side by side in ParaView.

**Not closed.** The flip pass (the output has slivers, see U3), a sizing field for
flat ground (U2), catchment clipping, and breaklines through the refined interior.

## The user's decisions this design is built on

1. **Error measure: the sup-norm.** A triangle's error is the maximum, over the
   DEM nodes it contains, of `|DEM value - the triangle's plane at the node|`.
   Refinement stops when every triangle's error is at most `tolerance`. "Contains"
   is defined exactly in R2.
2. **Concurrency.** Each round's per-triangle scan runs in parallel as a read-only
   pass writing one result slot per triangle. The split phase is single-threaded,
   in triangle order, so the output is identical for every thread count. No locks,
   no races. R7 picks the mechanism.

## What was measured

On this machine (Apple clang 21.0.0), in the scratchpad:

1. **`std::execution::par` does not compile.** `std::for_each(std::execution::par, ...)`
   fails with "no member named 'par' in namespace 'std::execution'". libc++ ships
   the parallel algorithms only behind `-fexperimental-library`, and libstdc++
   needs TBB for them. Either way, a dependency or an experimental flag.
2. **OpenMP is not available.** `c++ -fopenmp` fails with "unsupported option
   '-fopenmp'". Apple clang ships no OpenMP runtime; it would need Homebrew's
   `libomp` on every macOS machine and in CI.
3. **CI** runs the C++ core on `ubuntu-latest` and `macos-latest`, and a separate
   asan+ubsan job on `ubuntu-latest` (`.github/workflows/main.yaml`).

Not measured, and stated as arithmetic rather than a measurement:

4. **How many DEM nodes lie on the start grid's edges.** A stride-`s` grid has, per
   `s x s` block, one horizontal, one vertical and one diagonal edge (the CDT picks
   a diagonal in each cocircular square; either diagonal has slope +-1 in node
   units). Each passes through `s - 1` DEM nodes in its interior. So about
   `3/s` of all nodes lie on an edge: 15 % at stride 20, 2 % at stride 158. R3
   depends on this number being non-zero, not on its size.

Triangle counts and run times on the real tile are not predicted here. The
real-tile test (T10) reports them.

## Rulings

### R1. Refinement starts from a coarse stride grid, built exactly as in increment 12

- The start mesh is increment 12's: `grid_domain.subsample(meta, stride)`, then
  `build_pslg` -> `node` -> `triangulate`. No new start-mesh code.
- **`--stride` means the start grid's spacing in DEM nodes.** It is also the
  largest triangle the output can have: refinement only splits. So it doubles as a
  global size cap, which is the only size control in this increment (U2).
- The start grid is not the tile outline alone. Two triangles over a 50 km tile,
  refined by fans and without flips, give fans of slivers radiating from the
  corners across the whole tile. A coarse grid bounds how thin a fan can get.
- **Default stride with `--tolerance`:** at most 33 nodes on the longer side,
  `max(1, ceil((max(rows, cols) - 1) / 32))`. For the fixture that is 158, so a
  1.58 km start spacing and about 2 000 start triangles. Increment 12's default
  (256 a side) would start at 128 000 triangles, and most would never be refined:
  it would hide the adaptivity the user wants to see. A named constant in
  `grid_domain.py` next to `MAX_NODES_PER_SIDE`.

### R2. What "the nodes a triangle contains" means

Work in **lattice coordinates**: every mesh vertex is a DEM node `(row, col)`.
This holds for the start grid (increment 12, measurement 4: snapping moves no
node at UTM scale) and for every inserted point (R3). So:

- **The node set of triangle `T = (a, b, c)`** is every DEM node `p` in the
  *closed* triangle, minus its three vertices. Membership is three integer
  orientation tests in `int64` on `(col, row)`. They are exact, so "on an edge"
  is `orient == 0`, not a tolerance. The largest value is about `(5051)^2 * 2`,
  far below `2^63`.
- **Nodes on an open edge** belong to both triangles that share the edge. Both
  planes agree there: each is the linear interpolant of the two edge endpoints.
  Rounding can make the two computed values differ in the last bits; the tie is
  settled by R3's deterministic order, never by thread timing.
- **Vertices** have error exactly zero: z at a vertex is the node's own value
  (`value_at`, not `bilinear`), so the plane passes through it. They are excluded
  from the set, and so can never be the maximum.
- **NoData nodes** (the sentinel or NaN, `is_nodata`) are skipped. They add
  nothing to the error.
- **Nodes outside the grid** cannot occur. Every vertex is an in-grid node and the
  triangle is their convex hull, so its bounding box lies in the grid. This is an
  assertion, not a branch.
- **The plane** at `p` is `(o_bc * z_a + o_ca * z_b + o_ab * z_c) / (2A)`, where
  the `o` are the three integer orientations of `p` against the edges and `2A` is
  `orient(a, b, c)`. Integers up to `2^53` convert to `double` exactly. The error
  is `|z_p - plane|` in `double`, in metres.
- **The scan result** for `T` is its maximum error, the node where it occurs, and
  where that node lies: strictly inside, or on edge 0, 1 or 2. Ties go to the
  smallest `(row, col)`, so the result does not depend on the scan order.

A triangle **converges** when its maximum error is at most `tolerance`. A
triangle with no nodes in its set has error 0.

### R3. The inserted point, and why edges must be splittable

The point inserted into `T` is its maximum-error node. It is a DEM node, so its z
is exact.

**The pure fan cannot deliver the sup-norm.** A fan inserts strictly inside `T`
and never touches `T`'s edges. So the interpolant along an existing edge never
changes, and the error at a node on that edge can never go down. By measurement 4,
about `3/s` of all nodes sit on a start-grid edge, and every fan edge between two
nodes whose row and column offsets share a factor passes through more. A ridge
crossing a start-grid edge would keep its error forever, and the loop would never
end. The plan's "no hanging nodes by construction" argument is correct, but it
only covers interior points.

**Ruling: split whatever element holds the node.**

- **Strictly inside `T`:** fan into three children `(a, b, p)`, `(b, c, p)`,
  `(c, a, p)`. `T`'s edges are unchanged, so no neighbour sees a new vertex.
- **On edge `e` of `T`, with a neighbour `U` across `e`:** split `T` into two
  and `U` into two, at `p`. Both triangles that share `e` are split, so no
  hanging node appears. This is the only case that touches a neighbour.
- **On edge `e` of `T`, with no neighbour** (the tile boundary): split `T` into
  two.
- **On a vertex:** impossible (R2).

Every child has positive area. A point strictly inside a triangle, or strictly
inside an edge, cannot make a degenerate child, and the integer tests make
"strictly" exact.

**Termination.** Every insertion adds a DEM node with error above `tolerance`,
which is therefore above zero, which therefore is not already a vertex. The
lattice is finite. So the loop ends, at worst with every valid node a vertex and
every error zero. This holds for `tolerance == 0` too, and needs no continuity
argument.

**Constraint edges.** A split point on a constrained edge is exactly on the
segment (`orient == 0` in integers), so the two halves cover the same segment.
Both halves keep the parent's constrained bit and its property mask. In this
increment the only constraint is the outer ring, mask 0. The same rule serves
breaklines later: it adds Steiner points on constraints but never moves them.

### R4. The data structure: a flat triangle array with adjacency, not a tree

The plan's ternary tree cannot represent R3's edge split: splitting `e` changes
two triangles at once, and a tree stores no neighbours. So:

- `LatticeMesh`: vertices as `(row, col)` pairs, `uint32`; triangles as three
  vertex indices, counter-clockwise; per triangle, three neighbour indices (or
  none) and a 3-bit constrained mask plus three 32-bit property masks (increment
  7's `EdgeProperties` bits). Edge `k` runs from vertex `k` to vertex `k+1`, as in
  `IndexedMesh2`.
- **Splits reuse the parent's slot for the first child and append the others**,
  in a fixed order. The array stays dense, so there is no flatten step: the
  output mesh *is* the array.
- It is built from an `IndexedMesh2` plus the constraint edges and masks. Building
  it checks two things and refuses otherwise: every vertex is bit-equal to
  `RasterGeometry::node` of some `(row, col)`, and every triangle has positive
  integer orientation.
- It depends on `core` only. It knows no raster; `(row, col)` are just integers.
  This matches the dependency diagram: `mesh` below `refinement`.
- Adjacency costs 12 bytes a triangle and is what makes a flip pass cheap later
  (U3).

**Output.** Vertices go back to world `(x, y)` through `RasterGeometry::node`,
with z from `value_at` and `valid = !is_nodata`. Constraint edges are the
constrained triangle edges, each emitted once, with their masks. That is exactly
what `elevation.trim` already takes, so the trim stays unchanged.

### R5. The round

```
active = every triangle, in index order
loop:
    parallel scan: result[i] = scan(active[i])            -- read-only
    marks = active triangles that did not converge, in index order
    if marks is empty: stop
    serial, in order, for each marked T:
        if T was already split this round: skip
        if the split needs neighbour U and U was split this round: skip
        split (R3); mark every resulting triangle as split this round
    active = every triangle split this round, plus every marked triangle
             skipped this round, in index order
```

- A marked triangle is skipped for one of two reasons. Either it was already
  split as someone's neighbour, so its pieces are active; or the neighbour its
  split needs was split this round, and then it is unchanged. *Amended at
  green:* the design first said only split triangles stay active, which dropped
  the second case for good and could leave a triangle above tolerance. Skipped
  triangles stay active, so both cases are scanned again next round and nothing
  is lost. Termination holds because every round with a mark splits at least
  the first marked triangle in index order.
- Triangles that converged and were never touched keep their scan result and are
  not scanned again.
- Everything the split phase does depends only on the marks and the index order.
  Neither depends on the thread count, so the output is bit-identical for any
  thread count. That is the determinism argument; T6 tests it.
- The final achieved maximum error is the largest stored result over triangles
  with three valid vertices.

### R6. NoData

A triangle with a NoData vertex has no plane. Call it a **void triangle**. Start
vertices can be NoData (the fixture's row 0 is entirely NoData, column 0 is 57 %);
inserted vertices never are, because they are chosen from valid nodes.

- **Void triangles are carved, not refined.** If a void triangle's node set holds
  any valid node, it is split at the valid node nearest a NoData vertex (squared
  lattice distance, ties to the smallest `(row, col)`). Its children away from
  the void have all-valid vertices and are then refined normally. This repeats
  until each void triangle holds no valid node.
- Without carving, a start triangle touching the void would be dropped whole.
  With R1's 1.58 km default stride, that loses a 1.58 km strip along the whole top
  edge and most of the left edge of the fixture. Carving loses at most the thin
  triangles between the void and its nearest valid nodes.
- **The trim stays.** After refinement, increment 12's `elevation.trim` drops
  void triangles and counts their NoData vertices, unchanged.
- The scan also counts, per void triangle, valid nodes it still holds when it
  converges. Carving continues until no void triangle's closed node set holds a
  valid node, so the total is always zero when the loop ends. It still goes into
  `elevation_source` as "valid DEM nodes not covered", so a future change that
  breaks carving shows up in the file, not silently.
- Increment 12's one-cell trim around every void (bilinear refuses a cell with
  any NoData corner) does not apply here: vertices are nodes, and their z is
  read directly.

### R7. The concurrency mechanism: `std::jthread` over contiguous chunks

- The scan partitions the active list into `threads` contiguous chunks. Each
  thread scans its chunk and writes `result[i]` for its own `i` only, into a
  vector sized before the threads start. Threads are joined before the split
  phase reads anything.
- **No locks, no atomics, no shared writes.** A deadlock needs a lock and a race
  needs two writers to one location; this has neither. The raster, the vertices
  and the triangles are not written during the scan.
- Threads are created per round and joined at its end. Rounds number in the tens
  to low hundreds, so the creation cost is noise. No pool: a pool is state that
  outlives a call, which is the one thing that makes this hard to reason about.
- Why not the alternatives. `std::execution::par` does not compile on macOS
  without an experimental flag, and needs TBB on Linux (measurement 1). OpenMP
  needs `libomp` on macOS (measurement 2). Both add a dependency or a
  per-platform build switch for a loop that `std::jthread` expresses in a dozen
  lines. `std::jthread` is C++20, in both standard libraries, and needs only
  `Threads::Threads`, which the C++ tests already link.
- `threads == 0` means `std::thread::hardware_concurrency()`, and 1 if that
  returns 0. A chunk is never empty: with fewer triangles than threads, fewer
  threads start.
- The helper is `for_each_chunk(n, threads, fn)` in
  `include/terrain/parallel_util/chunks.hpp`, header-only, knowing nothing of
  meshes, so it is testable alone.
- **The binding releases the GIL** around the whole refinement, as `sample` does.
  No Python object is touched inside.
- **A TSan CI job: yes.** This is the project's first multi-threaded code path, and
  asan+ubsan cannot see a data race. TSan cannot share a build with ASan, so it is
  a separate job on `ubuntu-latest` that builds only the four refinement test
  binaries and runs them directly (ctest names tests by Catch case, so
  `ctest -R refinement` matched only one), with `vm.mmap_rnd_bits=28` so TSan
  starts on that kernel. A new refinement suite must be added to the job's
  target list. That keeps it to one small target, not a second
  full build. A TSan report of zero races on T6 is the independent check of the
  "no shared writes" claim above.

### R8. What crosses the boundary

In: the bound `RasterView`, the start `IndexedMesh2`, the constraint edges
`(E, 2)` and masks `(E,)` exactly as `_constraint_arrays` builds them today,
`tolerance` (float, metres) and `threads` (int, default 0).

Out: one result object holding vertices `(M, 2)` float64, `z (M,)` float64,
`valid (M,)` bool, triangles `(K, 3)`, edges `(F, 2)`, masks `(F,)`, and four
numbers: rounds, points inserted, achieved maximum error, valid nodes not covered.

No path, CRS, lattice index or C++ pointer crosses. The lattice is internal to
C++, which already has `RasterGeometry`.

Refusals come back as a status and message, like the CDT's: an off-lattice
vertex, a clockwise or degenerate start triangle, or a negative or non-finite
tolerance.

### R9. The CLI

- `mesh --dem PATH --tolerance METRES`. With `--tolerance`, the pipeline is:
  decode -> `subsample(stride)` -> build/node/triangulate -> `to_core` ->
  `refine` -> `trim` -> write. `refine` replaces `sample`.
- `--tolerance` must be finite and `>= 0`. With a fixture instead of `--dem`, it is
  a usage error, as `--stride` is.
- `--stride` with `--tolerance` sets the start grid (R1). Without it, the R1
  default applies.
- `elevation_source`, ASCII, for example:
  `refined from DEM nodes, tolerance 1 m, achieved max error 0.9987 m, start
  stride 158, 0 valid DEM nodes not covered, 612 vertices without data dropped`.
  The achieved error is printed with enough digits to show it is not above the
  tolerance. The vertical-unit note from increment 12 is appended when it
  applies.
- The same counts go to stderr, as increment 12's drop count does.
- Without `--tolerance`: see U1.

## Ruled by the user

**The user chose U1 (a), U2 (a) and U3 (a) on 2026-09-25.** The options are
kept so the reasons stay on record.

### U1. Is refinement opt-in in this increment?

- **(a) Opt-in. Recommended.** `--tolerance` turns it on. Without it, increment
  12's uniform mesh is written, unchanged. You can open both in ParaView from
  the same tile, and no existing test or behaviour moves. A later increment makes
  refinement the default once the output has been seen.
- (b) On by default with a default tolerance (1 m). The uniform mode becomes
  `--tolerance` absent plus a new flag, or goes. Closer to the product, but it
  changes increment 12's shipped behaviour in the same PR that introduces the
  thing it is compared against.

### U2. Triangle size where the terrain is flat (the plan's first open question)

Without a size control, a lake stays as a few large triangles. The plan records
why that can be wrong for water routing.

- **(a) Later. Recommended.** In this increment the start stride is the only cap
  (R1): no triangle is larger than one start-grid cell. That is the plan's option
  1, the one it calls wrong for the long run because it spends triangles on
  isolated flat ground. It is free here, visible in ParaView, and lets you judge
  whether flat areas matter before a sizing field is built. The plan's
  recommended option 3, a blurred sizing raster read in the same scan, fits R2
  unchanged: one more lookup per node and one more convergence condition. It is
  a follow-up of perhaps 60 lines.
- (b) Now, as a sizing field. Adds the field's construction (a roughness measure
  and a blur, in C++ or numpy), a second CLI knob, and a second convergence rule
  to test. It would push this increment over the ceiling (see LOC).
- (c) Now, as a plain maximum edge length. Cheap (about 15 lines), but it is the
  global cap again, and `--stride` already gives one.

### U3. The flip pass (the plan's second open question)

A mesh made only of fans and edge splits has slivers: long thin triangles fanning
out from inserted points. They render fine and interpolate correctly, but they
are poor for any solver that runs on the mesh.

- **(a) No flip in this increment. Recommended.** Then the tolerance is a property
  of the mesh delivered: every output triangle's error is at most `tolerance`,
  and the plan's worry, that a flip moves the surface by up to the full data
  range at one saddle, does not arise. The output's quality is stated plainly:
  sliver-prone, not Delaunay. `elevation_source` already says the mesh was
  refined, not flipped. The flip becomes its own increment, where the plan's
  options 1 (flip only if the quad stays within tolerance) and 5 (re-triangulate
  a patch) can be weighed against a real mesh. R4's adjacency is what either
  needs, so nothing here has to be redone.
- (b) A tolerance-preserving flip now (the plan's option 1): Lawson flips that are
  kept only if every node in the new quad stays within tolerance. About 120
  lines, serial. Over the ceiling together with the rest.
- (c) A plain Delaunay flip now (the plan's option 4): tolerance stops being a
  property of the output. Not recommended: it undoes the user's choice of a
  pointwise guarantee.

## Prior art in `legacy/`

```sh
$ grep -rlE 'refine|tolerance|ratio' legacy/
legacy/rasputin/triangulate_dem.h
legacy/bindings.cpp
legacy/rasputin/application.py
legacy/rasputin/reader.py
legacy/rasputin/geo_tiff_reader.py
legacy/rasputin/solar_position.h
legacy/rasputin/mesh.py
legacy/tests/test_gml_repository.py
legacy/rasputin/web_visualize.py
```

The only hit on the subject is the `-ratio` knob (`legacy/rasputin/application.py:41`, "Mesh
coarsening factor"), which calls CGAL's Lindstrom-Turk edge-collapse
simplification with a count-ratio stop predicate (`legacy/bindings.cpp:333`).
That is decimation to a triangle *count*, not refinement to an error
*tolerance*, and it is a CGAL call. **Nothing is carried across.** The other
hits are unrelated uses of the words. `@migration-expert` is not needed.

## Files

| file | what |
|---|---|
| `include/terrain/parallel_util/chunks.hpp` | `for_each_chunk` over `std::jthread` (R7) |
| `include/terrain/mesh/lattice_mesh.hpp` | `LatticeMesh`, build with its two checks, the three splits, output (R4) |
| `include/terrain/refinement/scan.hpp` | the per-triangle scan over a `RasterSource`, including the void rule (R2, R6) |
| `include/terrain/refinement/refine.hpp` | options, result, status, the round loop (R5) |
| `bindings/core.cpp` | `refine(view, mesh, edges, masks, *, tolerance, threads=0)`, GIL released |
| `src_python/tin_engine/_core.pyi` | the stub and the result type |
| `src_python/tin_engine/grid_domain.py` | the refinement default stride (R1) |
| `src_python/tin_engine/cli.py` | `--tolerance`, the refusals, `refine` in place of `sample`, the sentence |
| `.github/workflows/main.yaml` | the TSan job (R7) |
| `project_structure.md` | `mesh/` is a flat array with adjacency, not a ternary tree; `parallel_util/` is a header; `window_for` stays unbuilt, since the lattice makes it unnecessary |
| `parallel_refinement.md` | "Data structure" and "Parallelism" point here; the two open questions record the user's rulings on U2 and U3 |

The scan reads rows through `RasterSource::row`, the first caller of the row
access increment 12 added for this purpose.

## Tests for `@tester`

**Invariant-critical (mutation testing required):** only these two.

- `test_mesh_lattice_split` — the three splits and adjacency. Conformity is
  decided here and nowhere else.
- `test_refinement_scan` — node-set membership (edge, vertex, NoData), the plane,
  the argmax and its tie-break. The tolerance guarantee is decided here.

Everything else is property or integration testing, not mutation-tested.

- **T1. Plane DEM.** `z = 3*col - 2*row + 7`, exact in float32. No insertion, one
  round, output triangles equal input triangles, achieved error 0.
- **T2. Single peak.** Flat DEM with one raised node.
  - Strictly inside a start triangle: one insertion, three children around it.
  - On an interior start edge: one insertion, four triangles replace two.
  - On the tile boundary: one insertion, two triangles replace one, both halves
    constrained with the parent's mask.
- **T3. Tolerance oracle.** An independent brute force (numpy, or a separate C++
  loop that tests every node against every triangle, sharing no code with
  `scan.hpp`): for every output triangle with three valid vertices, the maximum
  error over its node set (R2) is at most `tolerance + 1e-9 * max|z|`. Seeded
  random smooth and rough DEMs, at several tolerances including 0.
- **T4. Conforming.** Every undirected edge is in at most two triangles; those in
  one are exactly on the tile boundary; no vertex lies in the relative interior
  of any edge (exact integer test, all pairs, small meshes); every triangle has
  positive orientation; the sum of `2A` equals the start mesh's.
- **T5. Constraints preserved.** Every output constraint edge lies on an input
  constraint edge, and for each input edge the pieces' lengths (in lattice units)
  sum to its length. Masks are inherited.
- **T6. Determinism.** `threads` = 1, 2, 7 and hardware concurrency give
  bit-identical vertices, z, valid, triangles, edges, masks and the four numbers.
  In C++, and once through the binding. This is also what the TSan job runs.
- **T7. NoData.** A void corner: after refinement no kept triangle has a NoData
  vertex, the reported "valid nodes not covered" equals an oracle count, and
  carving stops. An all-NoData DEM gives an empty mesh and the CLI exits 2 with
  no file.
- **T8. Refusals.** An off-lattice start vertex, a clockwise triangle, and a
  tolerance that is negative, NaN or infinite, each give a status, not a crash.
- **T9. CLI.** `--tolerance 1 --out x.vtk` writes a file whose `elevation_source`
  names the tolerance and an achieved error not above it; `--tolerance` with a
  fixture is exit 2; without `--tolerance` the output is increment 12's, byte for
  byte (if U1 (a)).
- **T10. The real tile** at tolerances 5 m and 1 m, skipped without
  `imagecodecs`. Not marked slow: it runs only in the codecs CI step and takes
  about 2 s locally. Reports triangles, rounds, inserted points, achieved error and
  wall time to the test log. Asserts only that the achieved error is at most the
  tolerance and that there are fewer triangles at 5 m than at 1 m. No timing
  assertion.
- **T11. `for_each_chunk`.** Every index visited exactly once, for `n` below,
  equal to and above `threads`, and for `n == 0`.

## LOC estimate

Counted in `CLAUDE.md` §2's unit. An estimate, not a measurement; the
reconciliation follows the table.

| file | what | est. |
|---|---|---|
| `parallel_util/chunks.hpp` | chunking, `jthread`, join | ~30 |
| `mesh/lattice_mesh.hpp` | types, build and checks, adjacency | ~110 |
| | fan split, edge split (one or two sides), neighbour fix-up | ~110 |
| | output arrays, edges emitted once | ~35 |
| `refinement/scan.hpp` | bbox, integer tests, plane, argmax, void rule | ~80 |
| `refinement/refine.hpp` | options, result, status, the round | ~70 |
| `bindings/core.cpp` | `refine`, result object, variant dispatch | ~60 |
| `_core.pyi` | stubs | ~15 |
| `grid_domain.py` | default start stride | ~8 |
| `cli.py` | option, refusals, call, sentence | ~50 |
| | **total** | **~570** |

Increment 10 overran by 39 %; increment 12 came in on estimate. On the worse
bias this is near 790, over the ceiling. So the increment carries a
pre-agreed cut:

- **If the green branch measures over 700, `for_each_chunk` and the parallel scan
  move to 14b**, with the TSan job and T6's multi-thread cases. 14a runs the scan
  serially through the same `scan` function. That removes about 45 lines and
  changes no output, because R5's split phase never depended on the scan's
  thread count. The user still gets the adaptive `.vtk` from 14a.
- The cut is not the edge split or the void carving. Without the edge split, the
  loop does not terminate (R3). Without carving, the default stride loses a
  1.58 km strip of the fixture (R6).

U2 (b) or (c) and U3 (b) are each outside this estimate. Choosing either makes it
two increments.


### Reconciliation

Measured by `@reviewer` over the added production lines, excluding comments,
docstrings and raw-literal bodies: **668 non-blank** (`lattice_mesh.hpp` 175,
`refine.hpp` 174, `bindings/core.cpp` 93, `scan.hpp` 77, `cli.py` 62,
`_core.pyi` 41, `chunks.hpp` 40, `grid_domain.py` 4, `CMakeLists.txt` 2). That
is about 17 % over the ~570 estimate and under the 700 ceiling, so the agreed 14b cut
did not fire.

## Acceptance

- `rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --out tile.vtk`
  writes a file ParaView opens, with visibly finer triangles on slopes than on
  flat ground, and `elevation_source` stating the tolerance and an achieved error
  not above it.
- All gates in `CLAUDE.md` §4 green, including the TSan job; CI green.

## Not in scope

- The flip pass, and any triangle-quality improvement (U3).
- A sizing field or edge-length cap beyond `--stride` (U2).
- Breaklines or holes inside the tile; clipping to a catchment.
- A GPU scan.
- A `--threads` CLI option. The output does not depend on it; the binding takes it
  for tests.
