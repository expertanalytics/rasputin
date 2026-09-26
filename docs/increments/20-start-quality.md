# Increment 20 — a quality start: minimum-angle Steiner nodes before DEM refinement

Status: **implemented, in review, on provisional defaults.** Red `f7ae381`;
green `affc955`, `a66632c` (ASCII wording), `f5489ab`. Ola asked the main
session to "work until we can check the fruits of our efforts against the
previous run", so the build proceeds on the recommended C1 (a), C2 (a) 25°
and C3 (a). These are the main session's provisional picks, **pending Ola's
review**, not Ola's rulings. C4 is not provisional: Ola sent it to increment
20b, a minimum insertion distance (section "Choices for Ola").
Written by `@architect` before `@tester`, per `docs/increments/README.md`
step 1, on branch `increment20-start-quality` off master `a9a93bc` (increments
14 to 18).

**Closes.** The boundary fans Ola saw in increment 16's output. Ola's idea, in
their words: "At the initial cdt, we could have a set of default quality
criterias. This way, we would not start out with a [bad] triangulation." After
this increment,

```sh
rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif \
    --domain quarter.geojson --tolerance 10 --out quarter.vtk --stats -
```

starts refinement from a mesh whose every triangle has a minimum angle of at
least 25°, unless a stated reason prevents it (R5), and the report says how
many quality nodes were added.

**Not closed.** An angle bound during DEM refinement (R8, C3). Slivers from an
off-node input vertex next to a DEM node (R7, C4). A sizing field or maximum
edge length (increment 14's U2; C2). Splitting input segments for a guaranteed
bound (C1 (b)). Increment 19 (carving) is parked by Ola and not touched.

## The problem, measured

The quarter circle of increment 16's T-real: 536 ring vertices every 200 m,
235 of them off-node. Its start CDT is 534 triangles, all fans from the dense
ring into the empty interior: median minimum angle 0.8°, 64 % under 1°, 97 %
under 10°. DEM refinement never inserts a point in flat sea, so the fans
survive into the output (row "none" in M1).

## What was measured

In the scratchpad (`i20/proto.py`), never in the tree. A **batch** Delaunay
refinement in Python: CDT through the CLI's own `_engine` (boundary chain plus
the Steiner points as loose vertices), find every bad triangle, insert its
Steiner point, re-triangulate, repeat. Candidates are taken worst first, and a
candidate within half a circumradius of one already taken this batch waits for
the next batch. The resulting start mesh then goes through `_core.refine` and
`elevation.trim` as `cli._dem_mesh` does. Apple clang build of master
`a9a93bc`, rebuilt and copied into the venv before the run, 10 threads. The
batch loop is not R4's serial algorithm; it gives the same kind of mesh, not
the same mesh. Its Python time says nothing about the C++ pass (R10).

"Snapped" means the circumcentre is replaced by the nearest DEM node (R3).
"Floor" is the smallest circumradius the pass acts on (R5). "Splits" is how
segment encroachment is handled: none (C1 a), Ruppert's diametral circle, or
Shewchuk's 30° diametral lens (C1 b). Angles in world plan view, by
`stats.quality`, over the trimmed output.

**M1. Start and final meshes.** Baseline row reproduces increment 16's T-real
exactly (427 779 and 30 547 triangles, worst 0.0117° and 0.652°).

| start | splits | start tris (Steiner, ring) | start worst ∠ | tol | triangles | median ∠ | < 1° (count) | < 10° | < 20° | worst ∠ | max deg | deg ≥ 12 | `refine` s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| none (today) | — | 534 (0, 536) | 0.46° | 10 m | 30 547 | 32.9° | 0.31 % (95) | 3.97 % | 20.6 % | 0.652° | 43 | 54 | 0.135 |
| | | | | 1 m | 427 779 | 45.0° | 0.03 % (119) | 0.92 % | 9.7 % | 0.0117° | 43 | 216 | 0.379 |
| 20°, snapped | none | 2 202 (834, 536) | 20.02° | 10 m | 31 613 | 33.7° | 0.01 % (3) | 2.35 % | 17.6 % | 0.616° | 16 | 28 | 0.029 |
| | | | | 1 m | 428 105 | 45.0° | 0.00 % (13) | 0.79 % | 9.6 % | 0.0117° | 21 | 182 | 0.223 |
| **25°, snapped** | **none** | **2 702 (1 084, 536)** | **25.00°** | **10 m** | **32 063** | **33.7°** | **0.01 % (2)** | **2.27 %** | **17.4 %** | **0.880°** | **15** | **26** | **0.028** |
| | | | | **1 m** | **427 639** | **45.0°** | **0.00 % (13)** | **0.80 %** | **9.6 %** | **0.0117°** | **21** | **184** | **0.225** |
| 30°, snapped | none | 4 572 (2 019, 536) | 30.01° | 10 m | 33 547 | 34.7° | 0.00 % (0) | 1.83 % | 15.8 % | 1.33° | 13 | 15 | 0.026 |
| | | | | 1 m | 430 065 | 45.0° | 0.00 % (8) | 0.75 % | 9.5 % | 0.0117° | 18 | 154 | 0.221 |
| 20°, snapped | circle | 7 486 (2 851, 1 786) | 20.02° | 10 m | 36 849 | 33.7° | 0.01 % | 2.03 % | — | 0.880° | 16 | 26 | 0.032 |
| | | | | 1 m | 433 969 | 45.0° | 0.00 % | 0.80 % | — | 0.0280° | 21 | 195 | 0.228 |
| 20°, snapped | lens | 5 427 (2 082, 1 265) | 20.02° | 10 m | 35 128 | 33.7° | 0.01 % | 2.19 % | — | 0.352° | 16 | 32 | 0.031 |
| | | | | 1 m | 431 202 | 45.0° | 0.00 % | 0.80 % | — | 0.0280° | 20 | 186 | 0.226 |
| 20°, **not** snapped | lens | 5 543 (2 129, 1 287) | 20.00° | 10 m | 35 320 | 33.7° | 0.01 % | 2.07 % | — | 0.351° | 15 | 28 | 0.036 |
| | | | | 1 m | 430 712 | 45.0° | 0.00 % | 0.79 % | — | 0.0280° | 21 | 178 | 0.233 |

Every run met its tolerance (achieved 9.9969 to 9.9999 m, and 0.99998 to
1.00000 m). A floor of one cell gave the same mesh as two cells at 20°.

What it shows:

- **The fans go.** At 10 m, triangles under 1° fall from 95 to 2 and the
  maximum vertex degree from 43 to 15, for 5 % more triangles (25°). At 1 m the
  triangle count is unchanged within 0.1 % and the maximum degree halves.
- **`refine` gets faster**: 0.135 to 0.028 s at 10 m, 0.379 to 0.225 s at 1 m.
  The first rounds no longer scan 534 domain-sized triangles.
- **Splitting input segments buys nothing here.** Circle and lens modes add 729
  to 1 250 ring vertices and 10 to 20 % more triangles at 10 m, and the final
  quality is no better. The quarter circle is convex and densely sampled, so
  every bad triangle's circumcentre already lies inside the domain.
- **Snapping costs nothing measurable.** Snapped and unsnapped lens runs agree
  to within 1 %.
- **The worst angle at 1 m does not move** (0.0117°, 3 triangles under 0.1° in
  every snapped no-split run). It is increment 16's off-node arc vertex beside a
  DEM node, and start quality cannot touch it (R7). The split runs make it
  worse (0.0280°): a new off-node midpoint beside a node, the same cause.
- **DEM refinement makes its own thin triangles.** With any start, 16 to 18 %
  of the final triangles at 10 m and 9.5 % at 1 m are under 20°. They are what
  an angle bound during refinement would have to fix (R8).
- 25° against 20°: 250 more Steiner nodes and 450 more final triangles at 10 m,
  and one fewer triangle under 1°. 30°: 935 more nodes than 25° and 1 480 more
  triangles, for none under 1°.

## Rulings

### R1. Detria offers nothing; the pass is ours, on `LatticeMesh`

- `lib/detria/README.md`: only detria's predicates are used, "the triangulator
  itself is not used and will not be". detria has no quality refinement anyway
  (`grep -n -i "steiner\|refine\|quality\|angle" lib/detria/detria.hpp` finds
  only its own Steiner-point input and triangle classification).
- **A new header, `include/terrain/mesh/quality.hpp`**, working on a
  `LatticeMesh` in 14b's `LatticeFrame`. It reuses `split_inside`,
  `split_edge`, `legalise_around` and the exact kernel. It depends on `core`,
  `predicates`, `lattice_mesh.hpp` and `lawson.hpp`. **It knows no raster and
  reads no height.** The grid's `rows` and `cols` come in as two numbers.
- **`refine` calls it once**, after `legalise_all` and before the first scan,
  when `RefineOptions::min_angle_deg > 0`. Rejected: a separate binding that
  returns a new world-coordinate start mesh. It would round-trip every vertex
  through `to_lattice` a second time and duplicate the mesh plumbing, for no
  gain in isolation: `quality.hpp` is testable on a hand-built `LatticeMesh`
  without a DEM.

### R2. The criterion: minimum angle, no size bound (C2)

- A triangle is **bad** when its circumradius-to-shortest-edge ratio exceeds
  `1 / (2 sin θ)`. That is the same as its minimum angle being below θ, and
  needs no trigonometry per triangle. Computed in `double` in the frame. It
  only proposes work; nothing topological depends on it.
- **Default θ = 25°.** Ruppert's proof stops at 20.7°, but that proof is not
  what ends this pass (R6); the lattice is. So θ is a cost choice. M1: 25° is
  where the fans are gone for 5 % more triangles at 10 m; 30° costs 1.9 times
  the Steiner nodes of 25° for one fewer sliver. See C2.
- **No maximum edge length.** M1 fixes the fans by angle alone. A size bound is
  increment 14's U2 (a sizing field for flat ground) and belongs with it.
- **Geometry only.** The pass reads positions, never `z`, never NoData. Which
  nodes exist is the lattice, which is geometry. See R9 for what NoData then
  does.

### R3. Steiner points are DEM nodes

- The Steiner point for a bad triangle is **the DEM node nearest its
  circumcentre**: `col = round(cx / dx)`, `row = round(-cy / dy)` in the frame.
- **Why.** Ola's direction for insertion (ROADMAP, the general point insertion
  policy): every vertex comes from the input geometry as given, or from
  refinement against the DEM, and nothing else. A snapped Steiner point is a
  DEM node, like every point refinement inserts. It is not a loose point: its
  z is `value_at`, exact, and it is excluded from node sets as every vertex is.
  An unsnapped circumcentre would be a third kind of vertex, off-node with
  bilinear z, and a new source of R7's slivers.
- **Exactness is free.** The circumcentre is a `double` proposal. The inserted
  point is an integer node, and every decision after that (location, inside or
  on an edge, legalisation) goes through the existing exact predicates and the
  existing `split_inside` and `split_edge`. No new exact code.
- M1: snapping changes the counts by under 1 %.

### R4. The algorithm, serial

1. **Queue.** Every bad triangle, keyed by `(ratio descending, slot
   ascending)`. Each entry also stores the triangle's three vertex indices.
2. **Pop.** An entry whose slot no longer holds those three vertices is stale
   and dropped. A current entry needs no re-test: the same three vertices give
   the same ratio, so the slot-identity check is the whole staleness test.
3. **Skip rules** (R5), in this order: circumradius below the floor;
   circumcentre outside the node rectangle; the snapped node is already a
   vertex; the walk to it crosses a constrained edge or leaves the mesh.
   A skip is counted and the triangle is not queued again unless rewritten.
4. **Locate.** A visibility walk from the bad triangle towards the node, using
   `orient_sign`: at each triangle, cross the first edge (in edge order) that
   has the node strictly on its far side. It ends in the triangle containing
   the node, or on an edge of it (`Collinear`), or at a vertex (skip). **Bound:**
   `triangle_count()` steps; hitting it is a skip, counted separately and
   asserted zero in the tests. A visibility walk does not cycle in a Delaunay
   triangulation; the bound is there because the proof is for the
   unconstrained case.
5. **Insert** with `split_inside` or `split_edge` (a node exactly on an edge,
   constrained or not; R6), then `legalise_around` with the written slots as
   seeds. Every slot written by the split or a flip is re-tested and queued if
   bad.
6. Stop when the queue is empty.

The floor puts the node within half a circumradius of the circumcentre (R5),
so it is strictly inside the bad triangle's circumcircle. When the walk reaches
it without crossing a constraint, legalisation therefore removes the bad
triangle; when it does not, the triangle is a recorded skip.

### R5. The floor, and what "meets θ" means

- **The pass acts only on triangles with circumradius `R ≥ sqrt(dx² + dy²)`**
  (14 m on a 10 m grid). Then the nearest node is at most `sqrt(dx² + dy²) / 2
  ≤ R / 2` from the circumcentre.
- Below the floor the triangle is DEM refinement's business: it is within
  about a cell of every node it covers.
- **The postcondition**, which is what the tests assert: after the pass, every
  triangle meets θ, or has `R` below the floor, or is recorded as skipped for
  one of R4's reasons. The pass returns the counts per reason. On the quarter
  circle the prototype left none below θ.
- This is not Ruppert's guarantee. Under C1 (a) a bad triangle whose
  circumcentre lies across a constraint stays. It happens at reflex corners and
  along long segments. It does not happen on the quarter circle (M1).

### R6. Constraint segments

Under C1 (a), recommended:

- **The pass never computes a point on a segment.** Input polygons and
  polylines keep exactly their input vertices, as increment 16 R2 has it.
- **A snapped node that lies exactly on a constrained edge** (`Collinear` in
  the exact kernel) is inserted with `split_edge`, which keeps the bit and mask
  on both halves. That is what DEM refinement has done since increment 14. On
  the quarter circle it happens on the two straight edges, which lie on the
  tile's border row and column. The node is on the segment exactly, because
  the predicate that says so is exact on integer coordinates.
- **Tile-border segments** are ordinary constraints to the pass. Today the only
  ones are the domain's own, or the stride grid's outer ring on the legacy
  path. When a mosaic lands (ROADMAP gap 6), internal tile seams are not
  constraints, and the pass does not see them.

Under C1 (b), for the record: split an encroached segment at the midpoint of
its two endpoints' fractional `(col, row)`. That point is exact for a segment
along a lattice row or column (one coordinate is an integer and stays one).
Otherwise it is within half an ulp per coordinate of the true midpoint, and the
constraint becomes a chain within that distance of the input line. That is
about 1e-9 m at UTM scale, six orders below the noder's 1 mm snap (16 U6). It
needs `split_edge` to take a `MeshVertex`, and its output z is bilinear.

### R7. Near-coincident input vertices and DEM nodes: not start quality's job

- **Start quality does not fix the 0.0117° sliver** (M1: unchanged in every
  snapped run). The sliver is born later, when refinement inserts a DEM node
  beside an off-node input vertex. No start mesh prevents that.
- **This increment makes it no worse** under R3. Snapped Steiner points are
  nodes, so they cannot sit beside one. The split variants (C1 b) did make it
  worse (0.0280°).
- What to do about it is Ola's (C4). Measured size: 3 triangles under 0.1° at
  1 m out of 427 639.

### R8. After the start: the criterion is not enforced during DEM refinement (C3)

- 14b's insertion keeps the mesh constrained Delaunay but does not bound
  angles. M1: 17 % of triangles at 10 m and 9.6 % at 1 m end under 20°, with
  any start. Those come from refinement itself, mostly along steep terrain and
  the coast.
- **Recommended: start only.** Ola asked for a good start, and that is what the
  measurements support cheaply.
- Enforcing it throughout means: after each round's serial phase, run R4's
  queue seeded with the touched slots, and mark every slot it writes touched
  (14b R2's rule), so the tolerance stays a property of the delivered mesh.
  About 25 more lines. Termination still holds (R10 plus 14b R4:
  every insertion is a new node). **Cost, unmeasured:** up to a Steiner node
  per triangle under θ, so on the order of 10 to 20 % more triangles at 10 m,
  more rounds, and a larger serial phase. A prototype should measure it before
  Ola chooses it.

### R9. NoData

The pass is geometry only (R2), so it can insert a NoData node. That vertex is
invalid, its triangles are void, and 14 R6's carving handles them, as it does
for an off-node start vertex in a NoData cell (16 R2). The quarter circle has
no NoData. A domain over a large NoData area would feed increment 16 R5's slow
carving, which is increment 19's subject (parked). One test fixture checks the
combination terminates and trims.

### R10. Determinism, termination, cost

- **Determinism.** The pass is serial, runs before the first parallel scan, and
  its order depends only on the mesh: the queue key is `(ratio, slot)`, the
  ratio a fixed `double` expression, ties broken by slot. The output is
  bit-identical for any `threads`. 14b's T6 covers it with quality on.
- **Termination.** Every insertion adds a DEM node that is not yet a vertex,
  from a finite lattice. Legalisation terminates (14b R4). A pop without an
  insertion removes an entry, and entries are only added by writes. So the
  pass ends for any θ and any input angle. **No terminator is needed** under
  C1 (a): there are no segment splits to ping-pong. Small input angles (a road
  meeting a river at 10°) just leave triangles in the wedge that shrink below
  the floor or whose node lies across a constraint; both are recorded skips.
  Under C1 (b), Shewchuk's concentric-shell splitting and "do not split a
  triangle whose short edge joins two segments meeting below 60°" would both be
  needed.
- **Cost.** About 1 000 insertions on the quarter circle, each a walk of a few
  steps, a split and about two flips (14b M5: about 250 ns per flip). Estimated
  well under 10 ms, against `refine`'s measured saving of 0.1 to 0.15 s.
  `@developer` reports the real figure in `--stats`.

### R11. CLI, options, report

- **`--start-min-angle DEG`, default 25.** `0` turns the pass off, and the
  mesh is then bit-identical to increment 18's (the file's provenance
  sentence adds `start quality off`). Refused: negative,
  non-finite, above 35 (in practice Delaunay refinement stops terminating
  somewhere above 33°; with the lattice it still ends, but the cost is no longer a start
  cost). Only meaningful with `--tolerance`; given without it, refused like
  `--domain` without `--tolerance`.
- It applies to both start paths. On the legacy stride grid the full cells are
  right isosceles triangles (45°) when `dx == dy`, so the pass inserts nothing
  there except, possibly, in a partial last row or column.
- **`RefineOptions::min_angle_deg`** (default 0 in C++, so the binding and the
  C++ tests keep 18's behaviour unless they ask). The CLI passes 25.
- **`RefineOutcome`** gains `quality_inserted`, `quality_skipped` (a total;
  per-reason counts are the C++ test surface only) and `quality_seconds`.
- **`elevation_source`** gains `start min angle 25 deg` (or `start quality off`).
  ASCII, not `°`: increment 13's guarantee 7 keeps every file string ASCII,
  found by `@developer` at green.
  The stderr report and `--stats`' refinement table gain the two counts, and
  the timings table a row `refine: start quality`.

## Choices for Ola

*Provisional (2026-09-26), pending Ola:* C1 (a), C2 (a), C3 (a), C4 (a),
picked by the main session so the build can reach a comparison run. Ola
reviews them against the result.

### C1. May the start pass add vertices on input polygons and polylines?

- **(a) No. Recommended.** Input geometry keeps exactly its vertices; only DEM
  nodes that already lie exactly on a segment are inserted there, as
  refinement already does. On the quarter circle this gives the full 25° (M1).
  The cost: the bound is best effort (R5). Bad triangles can remain at reflex
  corners and along long segments, most likely with 16b's interior polylines.
- (b) Yes, Ruppert-style: split encroached segments at their midpoints, off-node
  with bilinear z (R6), with a diametral lens and Shewchuk's terminator. Gives
  the angle bound for any input. On the quarter circle it adds 729 to 1 250
  ring vertices and 10 to 20 % more triangles at 10 m for no better result
  (M1). The new vertices are on the input lines, within an ulp, but they are
  not the input's vertices. About 90 more lines. Revisit with 16b if (a)
  leaves needles along roads.

### C2. The default criterion

- **(a) Minimum angle 25°, no size bound. Recommended** (R2, M1).
- (b) 20°, Ruppert's classical value: 3.5 % more triangles at 10 m instead of
  5 %, three triangles under 1° instead of two.
- (c) 30°: 10 % more triangles at 10 m, none under 1°.
- (d) Any of these plus `--start-max-edge METRES`: about 10 lines. Not
  recommended now; increment 14's U2 sizing field is the place for a size bound.

### C3. Enforce the angle during DEM refinement too?

- **(a) No, the start only. Recommended** (R8).
- (b) Yes, in every round. About 25 more lines; cost unmeasured, estimated 10 to
  20 % more triangles at 10 m. Measure first.

### C4. Slivers from an off-node input vertex beside a DEM node (R7)

- **(a) Leave them for now. Recommended.** 3 triangles under 0.1° in 427 639 at
  1 m. They do not affect the tolerance. Record them in the general point
  insertion policy discussion.
- (b) Refinement does not insert a DEM node within ε (say 1 % of a cell) of an
  existing vertex. That breaks the guarantee "within tolerance at every valid
  node" for those nodes; the guarantee would have to be restated with a slope
  bound. Not recommended.
- (c) An input vertex within ε of a DEM node is moved onto it. It removes the
  cause, but it is snapping, which Ola ruled out in increment 16 (former U2).
  Only if Ola wants to revisit that for ε at the millimetre scale, where it is
  below the noder's own 1 mm snap (16 U6).

## Prior art in `legacy/`

```sh
$ grep -rliE 'ruppert|chew|steiner|min_angle|minimum angle|Delaunay_mesher|Mesh_2|refine_Delaunay|angle bound' legacy/
legacy/rasputin/avalanche.py
```

The one hit is a false positive: "Ava**lancheW**arning" in a REST path
(`legacy/rasputin/avalanche.py:28`). The legacy never refined for quality; it
triangulated every DEM point with CGAL (14b's prior art). **Nothing is carried
across.** `@migration-expert` is not needed.

## Files

| file | what |
|---|---|
| `include/terrain/mesh/quality.hpp` (new) | `QualityOptions{min_angle_deg, rows, cols}`, `QualityOutcome` (inserted, skips per reason, walk-bound hits), `improve<K>(LatticeMesh&, const LatticeFrame&, const QualityOptions&)`; R2 to R5. Core, predicates, `lattice_mesh.hpp`, `lawson.hpp`; no raster |
| `include/terrain/refinement/refine.hpp` | `min_angle_deg` in `RefineOptions`; the call after `legalise_all`; the three outcome fields |
| `bindings/core.cpp`, `_core.pyi` | `min_angle_deg` keyword on `refine`; the outcome fields |
| `src_python/tin_engine/cli.py` | `--start-min-angle`; refusals; the sentence and report (R11) |
| `src_python/tin_engine/stats.py` | two refinement counts, one timing row |
| `tests/cpp/CMakeLists.txt` | the new suite, in the TSan job's list |
| `project_structure.md` | `mesh/quality.hpp` in the tree |
| `ROADMAP.md` | row 20 |

## Tests for `@tester`

**Invariant-critical (mutation testing required): one suite.**

- **`test_mesh_quality`** (new). The postcondition and the constraints are
  decided here.
  - Q1 **postcondition** (R5), on a small analogue of the quarter circle (a
    dense off-node arc and two lattice-line edges, 64 × 64 nodes): every
    triangle meets θ, or has R below the floor, or is counted in a skip reason.
    The reasons are recomputed by the test from the output, not read from the
    pass.
  - Q2 **still constrained Delaunay** in the frame (14b R10's check) and every
    inserted vertex is a node.
  - Q3 **constraints**: the constraint edges as a set of lines are unchanged;
    a constraint edge is split only at a node the exact kernel calls
    `Collinear` with it; bits and masks survive on both halves.
  - Q4 **small input angle**: two breakline segments at 10° inside a square
    domain. The pass ends within a stated insertion bound; the wedge's triangles
    are recorded skips.
  - Q5 **walk bound** never hit, over Q1 to Q4.
  - Mutants to kill: no `legalise_around` after an insert (Q2 fails); the walk
    crossing a constrained edge (Q3); the floor removed (Q4's bound, and a node
    outside the circumcircle); accepting an existing vertex as the node
    (duplicate vertex, `build` of the output fails); the ratio test inverted.

Everything else is property or integration testing.

- **Q6. Frame.** `dx = 1, dy = 3`: a triangle that is 30° in the world but
  not in `(col, -row)` is judged in the world frame.
- **Q7. Off switch.** `min_angle_deg = 0` gives output bit-identical to
  increment 18 on 14b's T12 fixtures and on a domain start.
- **Q8. Stride grid.** With `dx == dy` and a stride dividing the grid, the
  pass inserts nothing on a stride start.
- **Q9. Outside the rectangle.** A boundary triangle whose circumcentre lies
  outside the node rectangle is skipped, not clamped.
- **Q10. NoData** (R9). A domain over a NoData block: the pass may insert NoData
  nodes; `refine` ends and the trim drops them.
- **14b's T3 and T6 re-run with quality on**: tolerance oracle, and determinism
  for threads 1, 2, 7 and hardware concurrency. T6 is what the TSan job runs.
- **CLI.** Default 25 appears in `elevation_source`; `--start-min-angle 0`
  says off; refusals: negative, NaN, 36, and without `--tolerance`; `--stats`
  shows the counts and the timing row.
- **T-real, logged, not thresholded** (except the tolerance): the quarter
  circle at 10 m and 1 m, reporting M1's columns. M1's 25° row is the
  reference; expect the serial pass to differ from the batch prototype in
  detail.

## LOC estimate

Counted in `CLAUDE.md` §2's unit.

| file | what | est. |
|---|---|---|
| `mesh/quality.hpp` | options, outcome, ratio, circumcentre and snap, walk, queue loop | ~100 |
| `refinement/refine.hpp` | option, call, outcome fields, timing | ~15 |
| `bindings/core.cpp`, `_core.pyi` | keyword, fields | ~10 |
| `cli.py` | flag, refusals, sentence, report | ~25 |
| `stats.py` | two counts, one row | ~10 |
| | **total** | **~160** |

On the worst overrun seen so far (+39 %), about 220. **It fits under 700, so it
is not split.** C1 (b) adds about 90, C2 (d) about 10, C3 (b) about 25; all of
them together still fit.

### What landed

The developer counted about 240 production lines against ~160, under 700;
`quality.hpp` is 148 of them (reviewer's count; 130 without blank lines),
mostly the walk and the insert-and-requeue step. The reviewer's total is about
222 (about 205 without blank lines), about 40 % over.
Settled at green: a snapped node the walk ends on counts as "already a vertex";
NaN or <= 0 turns the pass off; the two quality columns come last in the
`--stats` Refinement table; the file says `25 deg`, not `25°` (increment 13's
ASCII rule).

**Comparison on Ola's quarter circle** (`--binary --stats -`, "off" is
`--start-min-angle 0`, whose mesh is identical to increment 18's):

| | 10 m off | 10 m on | 1 m off | 1 m on |
|---|---|---|---|---|
| triangles | 30 547 | 31 581 | 427 779 | 428 225 |
| min angle median | 32.95° | 33.69° | 45.00° | 45.00° |
| share < 1° | 0.31 % | 0.01 % | 0.03 % | 0.00 % |
| share < 10° | 3.97 % | 2.54 % | 0.92 % | 0.83 % |
| worst angle | 0.652° | 0.287° | 0.0117° | 0.0117° |
| degree max | 43 | 16 | 43 | 18 |
| degree >= 20 | 7 | 0 | 9 | 0 |
| quality nodes | 0 | 644 | 0 | 644 |
| refine | 0.134 s | 0.036 s | 0.375 s | 0.240 s |
| total | 0.207 s | 0.106 s | 0.469 s | 0.333 s |

The boundary fans are gone. Two numbers differ from the batch prototype: 644
nodes inserted (prototype 1 084: the serial pass legalises and requeues after
every insert, so one node can clear several bad triangles), and the 10 m worst
angle is 0.287°. That triangle comes from DEM refinement, not from the pass:
`@reviewer` meshed at `--tolerance 100000` (no DEM refinement) and the start
is exactly 25° worst with the pass on (0.461° worst, 64 % under 1°, off). C3 (a)
leaves DEM-refinement triangles alone.
The 1 m worst angle is the near-node sliver that increment 20b addresses (Ola:
a minimum insertion distance scaled by z_tol / |grad z|).

## Acceptance

- The quarter circle at `--tolerance 10`: no fans from the boundary in
  ParaView; `--stats` shows start triangles near 2 700, under 1° near 0.01 %,
  maximum degree near 15, and `refine` faster than increment 18's.
- `--start-min-angle 0` reproduces increment 18's mesh bit for bit; the file's
  provenance sentence adds `start quality off`.
- Known, not fixed: a triangle degenerate in double arithmetic has a NaN ratio
  and is neither queued nor counted (reviewer); rare.
- All gates in `CLAUDE.md` §4 green, including the TSan job; CI green.

## Not in scope

- Segment splitting (C1 b), a size bound (C2 d), the angle during refinement
  (C3 b), and R7's slivers (C4).
- Increment 19's carving, parked by Ola.
