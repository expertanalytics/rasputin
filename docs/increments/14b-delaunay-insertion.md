# Increment 14b — Delaunay insertion in the refinement loop

Status: **implemented, in review.** Red `74d3f67`; green `8987ec7`, `508934d`,
`9c9d36a`; test fix `1aa245b`. Measured 124 non-blank production lines
(reviewer's count) against ~130. Real tile at `--tolerance 1`: 463 974
triangles, median minimum angle 45 deg, 0.0 % under 1 deg, worst 0.63 deg,
achieved max error 0.99998 m. Written by `@architect` before `@tester`, per
`docs/increments/README.md` step 1. The user chose C1 (a), C2 (a) and C3 (a) on
2026-09-26, and asked for a more general point insertion policy to be
discussed later (section "Ruled by the user").
Supersedes increment 14's U3 (a), "no flip pass" (see R8).

**Closes.** Triangle quality of the refined mesh. Increment 14 inserts by fans
and edge splits and never changes an edge, so children keep their parent's long
edges and the 1.58 km start triangles become needles. The user saw it in
ParaView; the main session measured it (median minimum angle 0.9° at 10 m, 50 %
of triangles under 1°). After this increment,

```sh
rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --out tile.vtk
```

writes a **constrained Delaunay** mesh whose every triangle is still within the
tolerance, with about a third fewer triangles than today.

**Not closed.** Quality refinement in the Ruppert sense (Steiner points added
for angle alone), a sizing field for flat ground (increment 14, U2), the
coastline as a constraint (R9), and a parallel insert-and-flip phase (R7).

## The user's direction this design is built on

Garland–Heckbert style Delaunay insertion. After inserting a node, legalise the
edges around it with Lawson flips, never flipping a constraint edge. Then
re-check the error of the triangles that changed. The sup-norm tolerance stays a
property of the delivered mesh.

## What was measured

In the scratchpad, never in the tree: a copy of `lattice_mesh.hpp` with a
`flip`, a copy of `refine.hpp` that legalises after each split, and a driver over
the real tile (`7908_3_10m_z33.tif`, 5051 x 5051, 10 m, dx = dy) with the start
mesh taken from the CLI's own `subsample` -> `_engine` path. Apple clang 21,
`-O2`, 10 cores, `threads = 0`. Triangle counts are before the trim (the trim
drops about 8 000 void triangles). "Min angle" is each triangle's smallest angle
in world coordinates, over triangles with three valid vertices.

**M1. Quality, today versus Delaunay insertion.**

| stride | tol | insertion | triangles | median min angle | < 1° | < 10° | worst |
|---|---|---|---|---|---|---|---|
| 158 | 10 m | fan (today) | 102 143 | 0.93° | 50.5 % | 81.6 % | 0.0015° |
| 158 | 10 m | Delaunay | 47 701 | 28.1° | 19.7 % | 22.0 % | 0.30° |
| 158 | 1 m | fan (today) | 678 522 | 2.11° | 35.8 % | 77.9 % | 0.0012° |
| 158 | 1 m | Delaunay | 444 204 | 45.0° | 1.8 % | 2.6 % | 0.30° |
| 80 | 10 m | Delaunay | 53 212 | 32.9° | 17.2 % | 19.3 % | 0.60° |
| 80 | 1 m | Delaunay | 450 571 | 45.0° | 1.8 % | 2.6 % | 0.27° |
| 40 | 10 m | fan | 98 071 | 14.5° | 14.6 % | 44.9 % | 0.020° |
| 40 | 10 m | Delaunay | 76 427 | 45.0° | 0.0 % | 12.3 % | 1.25° |
| 40 | 1 m | fan | 599 624 | 5.19° | 17.7 % | 64.2 % | 0.019° |
| 40 | 1 m | Delaunay | 473 040 | 45.0° | 0.0 % | 2.3 % | 0.63° |
| 20 | 10 m | Delaunay | 170 159 | 45.0° | 0.0 % | 5.1 % | 2.33° |
| 20 | 1 m | Delaunay | 563 213 | 45.0° | 0.0 % | 1.8 % | 2.49° |

Every Delaunay run: zero edges fail `incircle` (checked over every interior
edge of the output with `DefaultKernel`), and the achieved max error is below
the tolerance (0.99999 m, 9.9994 m).

**M2. Where the remaining slivers are.** At stride 158, all but 13 of the
triangles under 1° have three vertices at z = 0: open sea. 96.4 % of the tile's
valid nodes are exactly 0. A sea start cell is never refined, so a dense
refined coastline faces a start vertex up to 1.58 km away, and Delaunay cannot
avoid a thin triangle there: the apex angle is about coast spacing / distance.
That is a grading problem, not a flip problem. A smaller start stride removes
it (stride 40: 0.0 % under 1°). Stride 80 does not (17 % at 10 m).

**M3. Cocircularity.** In the Delaunay output, 12 % (stride 158, 1 m) to 17 %
(stride 20, 1 m) of interior (edge, opposite apex) pairs are exactly cocircular.
The tie rule (R3) is exercised constantly, not in corner cases.

**M4. Where the time goes, today's shipped CLI** (`mesh --dem ... --tolerance`,
stages timed by wrapping the functions `cli.py` calls):

| tol | total | decode | start CDT | refine | trim | write |
|---|---|---|---|---|---|---|
| 1 m, text `.vtk` | 1.62 s | 0.065 | 0.003 | 0.383 | 0.027 | **1.125** |
| 1 m, `--binary` | 0.50 s | 0.055 | 0.002 | 0.389 | 0.026 | 0.014 |
| 10 m, text `.vtk` | 0.42 s | 0.056 | 0.002 | 0.196 | 0.004 | 0.155 |

At 1 m the text writer is 69 % of the run. Refinement is second.

**M5. Inside `refine`** (prototype timers around the parallel scan and the
serial phase; carving is the void triangles' share of both):

| stride | tol | insertion | refine | scan (parallel) | split + flip (serial) | flips | rounds |
|---|---|---|---|---|---|---|---|
| 158 | 1 m | fan (today) | 0.42 s | 0.295 | 0.052 | 0 | 156 |
| 158 | 1 m | Delaunay | 0.47 s | 0.262 | 0.156 | 454 347 | 205 |
| 40 | 1 m | Delaunay | 0.27 s | 0.079 | 0.153 | 447 956 | 53 |
| 158 | 10 m | fan (today) | 0.20 s | 0.188 | 0.006 | 0 | 156 |
| 40 | 10 m | Delaunay | 0.07 s | 0.049 | 0.011 | 40 111 | 53 |

- The scan's cost is the bounding boxes it walks: 628 M node visits for today's
  1 m run, 111 M for Delaunay at stride 40. Needles have huge boxes; Delaunay
  triangles do not.
- NoData carving: about 8 000 splits and 3 % of node visits in every run. Not a
  cost worth a ruling.
- At stride 40, 1 m, the serial phase is 57 % of `refine` but about 30 % of a
  `--binary` run and under 10 % of a text run. About 2 flips per insertion,
  about 250 ns each including the filtered `incircle`.

**M6. Early exit in the scan** (the user's question, C2). Delaunay, stride 40:

| tol | insertion point | triangles | node visits | scan | serial | refine |
|---|---|---|---|---|---|---|
| 1 m | worst node (a) | 473 040 | 111 M | 0.079 | 0.153 | 0.27 s |
| 1 m | first over tol (b) | 598 470 | 95 M | 0.081 | 0.253 | 0.40 s |
| 1 m | first over 2·tol, else worst (c) | 536 965 | 96 M | 0.076 | 0.203 | 0.33 s |
| 1 m | first over 4·tol, else worst (c) | 501 488 | 99 M | 0.076 | 0.175 | 0.30 s |
| 10 m | worst node (a) | 76 427 | 83 M | 0.049 | 0.011 | 0.07 s |
| 10 m | first over tol (b) | 91 659 | 79 M | 0.049 | 0.018 | 0.08 s |

Early exit saves at most 15 % of node visits, because most visits certify
triangles that are within tolerance, and certification needs the full scan.
It costs 6 to 27 % more triangles. It is slower end to end in every case.

## Rulings

### R1. Legalisation happens inside the serial split phase

- The round of increment 14 R5 is unchanged in shape: a parallel read-only
  scan, then a serial pass over the marks in index order.
- **After each split, legalise around the new vertex `q`** before moving to the
  next mark. Seeds: every triangle the split wrote (the parent's slot, the
  neighbour's slot for an edge split, and the appended slots). For each seed,
  the edge opposite `q` is tested; if it is not constrained, has a neighbour,
  and the neighbour's apex is strictly inside the circle (R3), flip it and push
  both resulting triangles. Both contain `q`, so the test edge is always "the
  edge opposite `q`". A stack, popped last-in first-out, so the order is fixed.
- Everything the serial phase does depends only on the marks, the index order
  and the stack order. None of it depends on the thread count, so the output
  stays bit-identical for any `threads`. The scan stays parallel and read-only.

### R2. What is touched, what is skipped, what is rescanned

- **`touched` is every slot written this round, by a split or by a flip.** Flips
  reuse the two slots they flip, so the array stays dense and there is still no
  flatten step.
- The invariant that makes the scan results trustworthy: **a slot not written
  since its last scan holds the triangle that scan measured.** Error depends
  only on the triangle's three vertices and the DEM, so its stored result is
  still exact.
- The skip rules of R5 stay as they are. A marked `T` whose slot was written
  this round (split as a neighbour, or flipped by an earlier insertion's
  legalisation, which can reach well past `T` and `U`) is skipped: its stored
  result describes a triangle that no longer exists. The neighbour rule for an
  edge split stays too. It is no longer needed for validity (an untouched `T`'s
  neighbour pointer is always current), but dropping it changes 14's round
  structure for no measured gain.
- **Active next round** = every touched slot plus every skipped mark, in index
  order, as today. Nothing else is rescanned, and nothing else needs to be.

### R3. The predicate: `DefaultKernel::incircle` in a scaled lattice frame

- The circle test runs on `Point2{col * dx, -(row * dy)}`, a **`LatticeFrame`**
  built from `RasterGeometry` by `refine` and passed down. `LatticeMesh` stays
  raster-free.
- Why scaled and not raw `(col, -row)`: the circle test is not invariant under
  unequal axis scaling. With `dx != dy` the lattice-frame Delaunay is not the
  world Delaunay, and world is where the user judges the angles. Translation
  does not matter, so `x_min` and `y_max` are left out, which also keeps the
  coordinates small for the filter.
- **Exactness.** `FilteredKernel<DetriaExact>` returns the exact sign for the
  doubles it is given. The doubles are fixed per vertex, so the predicate is a
  consistent function of one point set. That is all the termination argument
  (R4) needs. With `dx == dy` (every tile so far) the points are integers
  times one constant and the answer equals the integer lattice answer.
- Not `int64`. The incircle determinant is degree 4 in coordinate differences;
  it fits `int64` only below about 2^15 nodes a side, and `__int128` is banned
  (`exact.hpp`). The kernel has no such limit.
- **Ties.** Flip only on `Incircle::Inside`. `Cocircular` is not flipped. M3
  shows ties on 12 to 17 % of edges, so this is the common case. The result is
  a Delaunay triangulation, not *the* Delaunay triangulation: which of two
  cocircular diagonals survives is decided by the deterministic order of R1,
  and that is what makes it reproducible.

### R4. Tolerance and termination

**The stop condition** is unchanged in words and now carries more weight: stop
after a round whose scan finds no unconverged triangle among the active ones.
Every inactive triangle is unwritten since a scan that found it converged (R2).
So every triangle in the delivered mesh was scanned in its final form and is
within tolerance. The guarantee is of the mesh delivered, as 14 U3 required.

**Termination, in three parts.**

1. *Each legalisation ends.* A flip on a strictly-inside apex strictly lowers
   the lifted surface (the triangulation lifted onto `z = x² + y²`) over the
   flipped quad and leaves it unchanged elsewhere, so no triangulation of the
   current vertex set can recur. There are finitely many. Constraints only
   remove flips from the choice, so the argument holds for the constrained
   case. It needs the sign to be exact for fixed points (R3): an inexact
   predicate can say "inside" both ways round and cycle. Flipping on a
   cocircular tie would also cycle, which is a second reason for R3's tie rule.
2. *Each round with a mark inserts a point.* The first mark in index order is
   untouched when the serial phase starts, and so is its neighbour, so it
   splits.
3. *Insertions are finite.* Each inserts a valid DEM node that is not a vertex
   (it is in a node set, which excludes vertices). Flips never remove a vertex.
   The lattice is finite.

Holds for `tolerance == 0`, and for void carving, which also inserts only
valid non-vertex nodes.

**Convexity.** A flip needs a strictly convex quad. A strictly-inside apex
guarantees one (the plan's "One thing the pseudocode omits" measures the same
fact); `flip` asserts both new orientations are positive in integers.

### R5. Constraints and the start mesh

- A constrained edge is never flipped. A boundary edge has no neighbour and is
  never a candidate. An edge split on a constraint keeps the constrained bit
  and mask on both halves, as today; the new spokes are unconstrained and
  legalisable. `flip` carries the four outer edges' bits and masks across and
  gives the new diagonal none.
- **The start mesh is legalised once before the first scan**
  (`legalise_all`: every unconstrained interior edge on a stack, each flip
  pushes its four outer edges). On a start mesh from our CDT this flips nothing
  (M1's zero violations were reached without it). It is there because the
  binding accepts any valid mesh, and because only a Delaunay start makes the
  output provably constrained Delaunay: point insertion plus legalisation keeps
  a CDT a CDT, and does not make a non-Delaunay mesh one.
- The output is then **constrained Delaunay** in the frame of R3, with respect
  to the constraint edges as split.

### R6. The start mesh and its stride (C1)

With Delaunay insertion the final mesh is the CDT of the final vertex set, so
start *triangles* stop mattering; start *vertices* still do. They are forced
vertices, and they are the only vertices in unrefined flat ground. M2 shows the
remaining slivers are sea triangles between a refined coast and a far start
vertex. So the default stride now trades triangles for grading, and it is the
user's choice (C1). The R1 argument of increment 14 for a coarse grid (it bounds
how thin a fan can get) no longer applies; the size-cap argument still does.

### R7. Performance, and the parallel path left out

- Measured, not estimated (M4, M5). At the recommended stride, `refine` goes
  from 0.42 s to 0.27 s at 1 m and from 0.20 s to 0.07 s at 10 m. Fatter
  triangles have smaller bounding boxes, which pays for the flips.
- The serial phase is the larger half of `refine` at stride 40 and 1 m (0.153
  of 0.27 s) but not of the run: the text `.vtk` writer is 1.1 of 1.6 s. **So
  the serial phase stays serial in this increment.** The first performance
  target is the writer, which is not this increment's (see "Not in scope").
- **The upgrade path, recorded so nobody re-derives it:** deterministic
  reservations (Blelloch, Fineman, Gibbons, Shun, *Internally deterministic
  parallel algorithms can be fast*, PPoPP 2012). Per sub-round: every candidate
  computes, read-only and in parallel, its cavity (the triangles whose
  circumcircle strictly contains its point and are reachable without crossing a
  constraint, plus the ring of neighbours whose links it would rewrite). Each
  cavity triangle is reserved by the lowest candidate index that wants it. That
  minimum can be taken without atomics by emitting `(triangle, candidate)`
  pairs per chunk and reducing by triangle. A candidate that holds all its
  reservations commits; commits touch disjoint slots and run in parallel;
  losers retry next sub-round. The winner set depends only on indices, so the
  output is independent of the thread count, though **not** identical to the
  serial order's output. Slot numbering for appended triangles needs a prefix
  sum over winners. It is worth building when a profile shows the serial phase
  above half of a `--binary` run on a real tile.
- Cheap serial wins first, if needed: the prototype's per-insertion
  `touched.resize` and the neighbour-index searches are unoptimised.

### R8. Increment 14's U3 is superseded

The user chose 14 U3 (a), no flip pass, because a flip *after* refinement would
change triangles nobody rescans, and the tolerance would stop being a property
of the output (`parallel_refinement.md`, the open question). Delaunay insertion
with rescans removes that reason: flips happen inside the loop, every flipped
triangle is rescanned, and the loop stops only when every triangle is within
tolerance (R4). It is the plan's option 2, "alternate refine and flip until both
hold", with its termination argument supplied. 14's doc records this under U3.

### R9. The coastline as a constraint: later

Not here. The 0 m contour runs between DEM nodes, so its vertices are not
lattice nodes. R2 of increment 14 ("every vertex is a DEM node") is what makes
the scan exact and the node sets well defined, and a contour breaks it. It needs
contour extraction in Python, a snap policy onto the lattice or a relaxed R2,
and the noder: its own increment. M2 says what it would buy: the remaining
slivers sit on the coast.

### R10. The quality property

- **On synthetic terrain the refined mesh is constrained Delaunay**: for every
  interior edge that is not constrained, the neighbour's apex is not
  `Incircle::Inside` the triangle's circle. Checked in world coordinates from
  the outcome, on a geometry with integer `x_min`, `y_max`, `dx`, `dy`, so the
  translation to world is exact and cannot flip a sign.
- **The real tile reports min-angle statistics** to the test log (median, share
  under 1° and 10°, worst), plus flips. No threshold. A threshold would pin M2's
  sea grading, which is a property of the data and the stride, not of the
  algorithm; the Delaunay property is the algorithm's and is asserted.

## Ruled by the user

**The user chose C1 (a), C2 (a) and C3 (a) on 2026-09-26.** The options are kept
so the reasons stay on record. The user also said: "We need to discuss a more
general point insertion policy, at some point." C2 is therefore a ruling for
this increment, not the project's last word on insertion: worst-node insertion
stays until that discussion, which should weigh the options in C2, the sizing
field (increment 14's U2), coastline constraints (R9), and non-node points.

### C1. The default start stride with `--tolerance`

- **(a) At most 129 nodes a side. Recommended.** `ceil((max(rows, cols) - 1) /
  128)`: stride 40 on the fixture. Measured (M1): 0.0 % of triangles under 1° at
  both tolerances, 473 k triangles at 1 m (today 679 k), 76 k at 10 m (today
  102 k), `refine` faster at both. Costs 6 % more triangles than stride 158 at
  1 m and 60 % more at 10 m, nearly all of them flat sea at the start spacing.
- (b) Keep 33 nodes a side (stride 158). Fewest triangles (444 k, 48 k), but
  20 % of triangles under 1° at 10 m, all in the sea (M2).
- (c) Something between. Stride 80 was measured and does not fix the sea
  (17 % under 1° at 10 m).

### C2. Which node is inserted (the user's early-exit question)

- **(a) Keep the worst node. Recommended.** Fewest triangles and fastest
  measured (M6). Deterministic as today.
- (b) The first node over the tolerance, in row-major order within the
  triangle's box. Deterministic too: the scan order is fixed and does not
  depend on threads. But it saves only 15 % of node visits (certifying a
  converged triangle still needs every node) and costs 27 % more triangles at
  1 m, which makes the serial phase and the output bigger. Slower end to end.
  With Delaunay insertion the extra points are placed where the error first
  shows, not where it peaks, so the fatter triangles are spent less well.
- (c) Hybrid: stop at the first node over `k · tolerance`, else insert the
  worst. `k = 2` gives +14 % triangles, `k = 4` +6 %; neither is faster than (a).
  A coarse-subsample prefilter was not measured; it would have the same ceiling,
  because certification dominates.

### C3. Keep the fan mode as an option?

- **(a) No. Recommended.** Delaunay insertion is always on with `--tolerance`.
  Fewer triangles, better shape, same guarantee, faster. A flag would keep a
  mode nobody should pick and double the test surface.
- (b) Keep a `--no-flip` switch, for comparison in ParaView. About 10 lines plus
  tests. The old output is also one `git checkout` away.

## Prior art in `legacy/`

```sh
$ grep -rliE 'incircle|flip|delaunay|legali' legacy/
legacy/rasputin/triangulate_dem.h
```

The hits are CGAL typedefs and calls: `Delaunay_triangulation_2` and
`Constrained_Delaunay_triangulation_2` with `Exact_predicates_tag`
(`legacy/rasputin/triangulate_dem.h:10-11`, `:48-49`), used to triangulate
the DEM whole. There is no incremental insertion or flip logic of its own.
**Nothing is carried across.** `@migration-expert` is not needed.

## Files

| file | what |
|---|---|
| `include/terrain/mesh/lattice_mesh.hpp` | `flip(t, e)`: two slots reused, four outer sides carried, integer orientation asserted (R4, R5) |
| `include/terrain/mesh/lawson.hpp` (new) | `LatticeFrame`; `legalise_around<K>(mesh, q, seeds, frame, on_write)` and `legalise_all<K>(mesh, frame, on_write)`, templated on `pred::GeometryKernel`. Depends on `core`, `predicates` and `lattice_mesh.hpp`; no raster |
| `include/terrain/refinement/refine.hpp` | the frame from `RasterGeometry`; `legalise_all` before round 1; `legalise_around` after each split, marking written slots touched; `flips` in the outcome |
| `bindings/core.cpp`, `_core.pyi` | `flips` on the result |
| `src_python/tin_engine/grid_domain.py` | the default start stride (C1) |
| `src_python/tin_engine/cli.py` | `elevation_source` says "constrained Delaunay"; flips in the stderr report |
| `tests/cpp/CMakeLists.txt` | link `terrain_predicates` into the refinement suites; the new suite in the TSan job's list |
| `project_structure.md` | `mesh/lawson.hpp` in the tree |
| `parallel_refinement.md` | "Final flip pass" and the open question point here |

`mesh/` gains a dependency on `predicates/`. Both sit below `refinement/`, and
`predicates/` depends on `core` only, so no cycle. `lattice_mesh.hpp` itself
stays core-only.

## Tests for `@tester`

**Invariant-critical (mutation testing required):** only these two.

- `test_mesh_lawson` (new) — `flip` and both legalisers. Conformity and the
  Delaunay property are decided here. Mutants to kill: flip on `Cocircular`
  (non-termination, so the test needs a flip-count bound), flip a constrained
  edge, drop one of the four repoints, lose an outer mask, use the raw lattice
  frame instead of the scaled one.
- `prop_refinement_refine`, the tolerance oracle (14's T3) re-run under
  Delaunay insertion. Mutant to kill: **a flipped slot not marked touched**,
  which leaves a stale result and must surface as an oracle failure.

Everything else is property or integration testing.

- **L1. Flip topology.** Two triangles, flip, check the vertex order, all six
  neighbour links in both directions, and that outer masks and constrained bits
  moved with their edges and the diagonal has none.
- **L2. Ties.** A unit square (cocircular): `legalise_all` flips nothing,
  whichever diagonal it starts with. A 17 x 17 lattice with random diagonals
  (seeded): it ends with no strictly-inside apex, and the flip count is bounded.
- **L3. Constraints.** A non-Delaunay quad whose diagonal is constrained is not
  flipped. The same quad unconstrained is.
- **L4. Frame.** A quad that is Delaunay in `(col, -row)` but not with
  `dx = 1, dy = 3`: flipped under the scaled frame, not under the raw one.
- **T12. Constrained Delaunay output (R10).** A cone and an island with a
  step coast (flat sea at 0, land rising from a few metres), about 129 x 129,
  strides 4 and 16, tolerances 5, 1 and 0; also once with `dx != dy`. Every
  unconstrained interior edge passes. Plus one start mesh with an interior
  constraint edge: its pieces are still constrained at the end.
- **T3 and T6 re-run** unchanged in intent: tolerance oracle, determinism for
  threads 1, 2, 7 and hardware concurrency. T6 is what the TSan job runs.
- **T2 amended.** Its counts assume fans ("three children"). Under Delaunay the
  insertion may flip; the test asserts the inserted vertex and that the
  triangles covering the old ones are Delaunay, not a fixed count. An amendment
  in the red commit, with the reason in the message.
- **T10 extended.** Reports min-angle median, share under 1° and 10°, worst,
  flips and rounds. Asserts only the tolerance and that 5 m gives fewer
  triangles than 1 m, as today.

## LOC estimate

Counted in `CLAUDE.md` §2's unit.

| file | what | est. |
|---|---|---|
| `mesh/lattice_mesh.hpp` | `flip` | ~25 |
| `mesh/lawson.hpp` | frame, `legalise_around`, `legalise_all` | ~70 |
| `refinement/refine.hpp` | frame, start pass, call after split, touched, `flips` | ~25 |
| `bindings/core.cpp`, `_core.pyi` | `flips` | ~6 |
| `grid_domain.py` | constant (C1) | ~2 |
| `cli.py` | sentence, report | ~4 |
| | **total** | **~130** |

Increment 14 came in 17 % over. On the worst bias seen here (39 %) this is
about 180, far under 700. C3 (b) would add about 10.

## Acceptance

- The acceptance command above writes a `.vtk` ParaView opens, with no needles
  on land, `elevation_source` stating the tolerance, an achieved error not above
  it, and "constrained Delaunay".
- T10's log shows the real tile's min-angle statistics near M1's for the
  chosen stride.
- All gates in `CLAUDE.md` §4 green, including the TSan job; CI green.

## Not in scope

- Steiner points for angle alone (Ruppert, Chew), and a sizing field (14 U2).
- The coastline, or any contour, as a constraint (R9).
- The parallel insert-and-flip phase (R7).
- The text `.vtk` writer's speed. M4 makes it the largest cost at 1 m, and it
  belongs to increment 13's writer, not to refinement.
