# Increment 20b — a minimum insertion distance from constraints, scaled by the tolerance

Status: **implemented, in review.** Red `5f86eeb`; green `d82de77`, `59d9fa5`;
test fixes `cb7eb96`. Ola chose C1 (a), C2 (a) and C3 (a) on
2026-09-26 (section "Ruled by Ola"); 20 and 20b land together in one PR, and a
retrospective follows. Written by `@architect`
before `@tester`, per `docs/increments/README.md` step 1, on branch
`increment20b-min-insertion-distance` off `increment20-start-quality` (unpushed).

**Closes.** Increment 20's C4, which Ola sent here: the 0.0117° sliver at 1 m on
the quarter circle, which the quality start cannot touch because DEM refinement
makes it. Ola's direction, in their words: "If we have a minimum insertion
distance from existing, constraining points? We would perhaps loose the optimal
insertion point, but we could bilinearly find a better point for our efforts.
This can become a fundamental issue later, so better fix it now." And:
"Epsilon could be related to the z-tolerance. After all, that's what we're
resolving here."

**Not closed.** Thin triangles that DEM refinement makes between nodes alone
(increment 20's C3, the 10 m worst angle of 0.287°). Thin triangles fanning from
an input vertex to a row of nodes nearly in line with it (M2's remaining worst,
0.107°). Increment 19 (carving), parked.

## The problem, measured: it is a segment, not a vertex

The brief read the sliver as an off-node input vertex A a fraction of a metre
from a DEM node N, joined by a tiny edge. **The output has no such edge.** In
the 1 m baseline the shortest edge in the whole mesh is 10.000 m. The three
triangles under 0.1° are needles:

| worst ∠ | vertices | what |
|---|---|---|
| 0.0117° | A1 (820260.632, 7900548.616), A2 (820255.980, 7900348.993), N (820260, 7900520) | N is **3.5 cm from the boundary segment A1-A2** (200 m long), 28.6 m from A1 and 171 m from A2 |
| 0.084°, 0.095° | A2 and consecutive nodes on column x = 820260 | nodes 0.23 to 0.25 m from the same segment |

The segment runs almost along lattice column 2051. Its linear z (0.49 m to
0.0 m) is far from the DEM beside it (3.3 m at N), so refinement must put a
vertex near the segment, and the best node it has lies centimetres from it.
Ola's instinct holds, with one change of object: **what constrains is the
constraint segment, and the input vertex is its special case at the ends.**

The vertex-only rule was prototyped as briefed and changes nothing (M1, row
"vertex disc"). 187 of the 235 off-node vertices have a node within ε, and not
one of those nodes is ever inserted. The reason is the tolerance argument
itself: a node beside an input vertex inherits almost its height through
bilinear z, so once A is a vertex of its triangle the node's error is small.
The same is not true beside the middle of a 200 m segment.

## What was measured

In the scratchpad (`i20b/`), never in the tree. A copy of this branch's
`scan.hpp`, `refine.hpp` and `lattice_mesh.hpp` with the rule behind environment
switches, built as `_core` in Release (Apple clang) and loaded in place of the
installed module by a driver that asserts it did. Each run is the real CLI
(`rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain …
--tolerance T --stats -`, quality start 25°, 10 threads). Angles and edges from
the output `.vtk`.

**The tolerance column is an independent oracle**, not the scan's own number:
every DEM node inside the output mesh, against the linear interpolant of the
output file's vertices and z (matplotlib's `TriFinder`), 7 071 502 nodes for the
quarter circle. It can fail: the same oracle at 0.9 m on the 1 m baseline
reports 13 945 nodes over.

Variants:
- **off**: this branch as it is.
- **vertex disc**: a node within ε of an off-node input vertex is never inserted.
- **grad**: the rule of R2 to R4, ε = clamp(tol / G, floor, half a cell) (R3).
- **fixed**: the same rule with ε = half a cell everywhere.
- **no fallback**: a node that got a foot is never inserted afterwards, even if
  its error stays above tolerance. Shown only to measure what R5 prevents.

"Feet" are points inserted on a segment instead of a node. "Fallback" is a
footed node inserted later anyway because its error stayed above tolerance (R5).

**M1. The quarter circle** (536 ring vertices, 235 off-node).

| run | tol | feet | fallback | triangles | worst ∠ | < 0.1° | < 1° | max deg | shortest edge | oracle max error | nodes over tol | `refine` s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| off | 1 m | — | — | 428 225 | 0.0117° | 3 | 19 | 18 | 10.000 m | 0.99998 | 0 | 0.23 |
| vertex disc | 1 m | — | — | 428 225 | 0.0117° | 3 | 19 | 18 | 10.000 m | 0.99998 | 0 | 0.25 |
| **grad** | **1 m** | **6** | **0** | **428 217** | **0.396°** | **0** | **12** | **18** | **9.997 m** | **0.99998** | **0** | **0.27** |
| fixed | 1 m | 6 | 0 | 428 217 | 0.396° | 0 | 12 | 18 | 9.997 m | 0.99998 | 0 | 0.27 |
| off | 10 m | — | — | 31 581 | 0.287° | 0 | 4 | 16 | 10.000 m | 9.99868 | 0 | 0.036 |
| grad, fixed, vertex disc | 10 m | 0 | 0 | 31 581 | 0.287° | 0 | 4 | 16 | 10.000 m | 9.99868 | 0 | 0.04 |

The baseline reproduces increment 20's comparison exactly (428 225 and 31 581
triangles, 0.0117° and 0.287°). The six feet are all on the segment A1-A2, for
nodes at 0.27 to 1.43 m from it (ε there 0.99 to 1.75 m, set by the slope, not
the cap). At 10 m the rule never fires and the mesh is identical. The remaining
1 m worst, 0.396°, is a node-only triangle in flat sea, 3.5 km long: DEM
refinement's own (increment 20's C3).

`refine` time: three repeats each at 1 m, 0.231 to 0.236 s off and 0.269 s with
the rule. Most likely the prototype's bitmap test on every scanned node (not
isolated), which the design does not have; the rule itself runs on six
insertions. `@developer` reports the real figure.

**M2. A steep adversary.** The quarter circle's arc lies in flat sea, so the
ε scale is not tested there (grad and fixed agree). Two synthetic domains in the
tile's steepest 4 km window (median slope 0.26, centre (823750, 7902250),
radius 1450 m, rotated 0.1234 rad so no edge follows the lattice): a 7-gon
(edges about 1.26 km) and a 60-gon (edges about 150 m). 1 m tolerance.

| domain | run | feet | fallback | triangles | worst ∠ | < 0.1° | < 1° | max deg | shortest edge | oracle max error | nodes over tol |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 7-gon | off | — | — | 17 851 | 0.000169° | 39 | 83 | 21 | 10.000 m | 0.99918 | 0 |
| | **grad** | **64** | **0** | **17 937** | **0.107°** | **0** | **19** | **14** | **10.000 m** | **0.99978** | **0** |
| | fixed, no fallback | 82 | — | 17 751 | 0.330° | 0 | 3 | — | 9.48 m | **1.94** | **3** |
| | fixed | 82 | 3 | 17 757 | 0.330° | 0 | 4 | — | **1.90 m** | 0.99990 | 0 |
| 60-gon | off | — | — | 19 598 | 0.000968° | 9 | 37 | 14 | 5.56 m | 0.99996 | 0 |
| | **grad** | **46** | **0** | **19 480** | **0.400°** | **0** | **6** | **14** | **5.56 m** | **0.99996** | **0** |
| | fixed, no fallback | 68 | — | 19 486 | 1.10° | 0 | 0 | — | 7.36 m | **2.16** | **1** |
| | fixed | 68 | 2 | 19 488 | 1.10° | 0 | 0 | — | 4.77 m | 0.99996 | 0 |

What it shows:

- **The needles go.** Under 0.1°: 3 → 0, 39 → 0, 9 → 0. The worst angle rises
  from 0.0117°, 0.000169° and 0.000968° to 0.40°, 0.11° and 0.40°. Triangle
  counts move by under 1 %.
- **The tolerance holds with the fallback, and only with it.** Without R5, a
  fixed half-cell ε leaves 3 and 1 nodes at 1.94 m and 2.16 m against a 1 m
  tolerance. With R5 every run is within tolerance at every node.
- **The slope scale is what makes the foot enough.** With ε = tol / G the
  fallback fired 0 times in 116 feet, and the largest error left at a footed
  node was 0.32, 0.44 and 0.53 m. With a fixed ε it fired 5 times in 150, and
  each firing puts a node within ε of a foot: a 1.90 m edge in the 7-gon. The
  fallback's edge has no lower bound (the node can be 3.5 cm from its foot), so
  a fixed ε can bring back exactly the defect this increment removes.
- **Fixed ε gives better angles here** (0.33° and 1.10° against 0.11° and
  0.40°): it catches nodes up to 5 m from a segment, where grad stops at about
  1 m on steep ground. That is C1.

## Rulings

### R1. Which vertices count as constraining: the constraint segments, ends included

- **The rule is judged against the constrained edges of the triangle holding
  the node**, as closed segments. A node near a segment's end is near an input
  vertex; that is the special case, handled in R2.
- **Not every existing vertex, and not the input vertices alone.** Measured
  (M1, "vertex disc"): the vertex rule does nothing on the quarter circle.
  Refinement-on-refinement near-duplicates cannot occur: before this increment
  every inserted point is a node, and two nodes are at least `min(dx, dy)`
  apart. Feet are the only new off-node points, and R2 keeps them ε from every
  vertex on their segment.
- **Only the triangle's own constrained edges.** A node within ε of a segment
  lies in the triangle on that segment unless that triangle is itself thinner
  than ε, which is the defect being removed. Scanning further (a spatial index
  over all constraints) is not needed for any measured case; see "Not in scope".

### R2. The rule, in the serial phase

When a marked, non-void triangle T's worst node N is about to be inserted (14b
R1, unchanged up to this point):

1. For each constrained edge `e = (a, b)` of T, in edge order: the foot
   `F = a + s (b − a)`, `s` the projection parameter clamped to `[0, 1]`,
   distances in the world frame `(col·dx, row·dy)`. Skip `e` when
   `dist(N, e) ≥ ε(N)` or when N lies on `e` exactly (`Collinear`: today's edge
   split at N, unchanged).
2. **If F is within ε(N) of `a` or `b`**, the nearest constraint point is an
   existing vertex: insert N as today. (R5 explains why this is safe.)
3. **Otherwise insert F instead of N**, with `split_edge(t, e, F)` on an
   off-node `MeshVertex`, then `legalise_around` with the same seeds as any
   edge split. The neighbour rule of 14 R5 applies: a touched neighbour defers T.
4. **Refusals, each inserting N as today and counted:** the exact orientation
   of any of the two or four children is not positive (F rounds off the segment
   by more than the child's height; measured 0 times); `vertex_z` refuses F
   (a NoData corner in F's cell; N itself is valid, so inserting it is sound).
5. Record N as **footed**. A footed node never gets a second foot (R5).

Only DEM refinement uses the rule. **Carving is unchanged**: a void triangle's
carve node is not an error maximum and the rule does not look at it.

### R3. ε: `tol / G`, clamped to a floor and to half a cell

- **`ε(N) = clamp(tol / G(N), floor, cap)`**, with
  - `G(N)`: the largest bilinear slope bound over the (up to four) cells that
    share N: per cell `hypot(max(|z01 − z00|, |z11 − z10|) / dx,
    max(|z10 − z00|, |z11 − z01|) / dy)`. Cells with a NoData corner are left
    out; with none left, `G = 0`.
  - `cap = min(dx, dy) / 2`. A node inserted further than that from a segment
    stands at least half a cell off it, so its triangle on the segment is not a
    needle. Every node within half a cell of a segment is then a foot
    candidate on flat ground.
  - `floor = min(dx, dy) / 100`: 10 cm on a 10 m grid. It exists only for the
    termination bound (R6), not for exactness: the exact kernel handles F at
    any distance, and R2's orientation check refuses a degenerate child.
    Measured ε never went below 0.77 m, so the floor did not bind.
- **Why `tol / G`.** Bilinear interpolation is Lipschitz with constant `G(N)`
  in each of N's cells, so for the foot F at distance `d < ε ≤ tol / G`,
  `|z_N − z_F| ≤ G·d ≤ tol`. Within ε of the segment the DEM differs from its
  value on the segment by at most the tolerance, which is Ola's "that's what
  we're resolving". Where the terrain is flat, `G → 0` and the cap applies.
- **What the scale does not prove.** The error at N afterwards is
  `|z_N − p(N)|` for the plane p of N's triangle, which is at most
  `|z_N − z_F| + |∇p|·d ≤ tol + |∇p|·ε`. `∇p` is not bounded by `G` in a
  thin triangle. So the scale makes the foot **usually** sufficient (M2: 116
  of 116), not always. R5 closes the gap.
- Computed in `double` in the serial phase from `value_at`; it proposes a point
  and decides nothing topological.

### R4. Increment 14 R2 is amended: inserted points are nodes or feet

Increment 14 R2 ("every inserted point is a DEM node"), as amended by 16 R2,
becomes: **every vertex refinement inserts is a DEM node, or a foot on a
constraint segment**. What that costs:

- **Height.** A foot's z is `vertex_z` (bilinear at its fractional position,
  16 R0), the same function the scan already uses for off-node start vertices.
  The output's world point is `(x_min + col·dx, y_max − row·dy)`; its z is
  `vertex_z` of the stored vertex, not `bilinear` of the world point, so the
  file carries the height the scan used.
- **The frames.** F is computed in `double` fractional `(col, row)` and is then
  a fixed vertex like any off-node start vertex. Orientation in `(col, −row)`
  and incircle in 14b's `LatticeFrame` treat it exactly (16 R2), and 16's
  known limit on the two frames applies to it unchanged. F is within rounding
  of the segment, not on it, unless the segment runs along a lattice row or
  column: the constraint becomes a two-piece chain within about 1e-12 cells of
  the input line, the case increment 20's R6 records for C1 (b). R2's exact
  child-orientation check is what keeps the split valid.
- **The scan is unchanged.** Error is still measured only at DEM nodes. A foot
  is a vertex, never in a node set; `is_node()` is false for it, as for an
  off-node start vertex. Increment 18's row spans already take fractional
  vertices.
- **Input geometry.** The segment keeps its input vertices and gains feet
  strictly inside it, as increment 14's edge splits already add nodes on a
  constraint. The feet carry the parent's constrained bit and mask on both
  halves (`split_edge`, unchanged). This is the first increment that adds
  off-lattice points on input segments; increment 20's C1 (a) (no segment
  splits in the *quality* pass) is not affected.

### R5. The tolerance guarantee is unchanged, because a footed node can still be inserted

- **The stop condition is untouched**: every valid node in every triangle is
  scanned, footed nodes included, and refinement stops only when every error is
  at most tol. So the delivered mesh is within tolerance at every valid node,
  exactly as 14b R4 states it. No restatement, no slope bound.
- **The fallback.** A footed node whose error is still above tolerance is the
  worst node of its triangle in a later round like any other; R2 step 5 stops a
  second foot, so it is inserted as today. That is the only way the rule can
  bring a short edge back, and it is what M2's "no fallback" rows show is
  required: without it, 3 and 1 nodes end at 1.94 m and 2.16 m.
- **With ε = tol / G the fallback did not fire** in 116 feet over three
  domains. Its short edges are the price of an unconditional guarantee; the
  alternative is C2 (b).
- R2 step 2 (F near an end, insert N) is safe for the same reason as the
  vertex-disc measurement: N beside a vertex inherits its height, and M1 shows
  such nodes are never selected in practice. *Correction at review:* step 2 is
  not counted anywhere (the code inserts N and counts nothing); a test now pins
  it with a fixture whose foot falls within ε of a segment end.

### R6. Termination

The old argument ("every insertion is a new DEM node") no longer covers feet.
The new one:

1. **Node insertions are finite**: each inserts a valid node not yet a vertex,
   from a finite lattice (14b R4.3), fallbacks included.
2. **Feet are finite**, two ways, either sufficient:
   - each foot is charged to a distinct node N (R2 step 5: footed once), so
     there are at most as many feet as nodes;
   - each foot is at least `ε ≥ floor` from every vertex on its segment when
     inserted (R2 step 2), so a segment of length L holds at most `L / floor`
     feet.
3. **Legalisation ends** (14b R4.1): a foot is a fixed point with exact signs.
4. **Each round with a mark inserts a point** (14b R4.2): the first mark
   inserts N or F, or a refusal inserts N.

Holds for `tolerance == 0` (then `ε = floor` wherever `G > 0` and the cap where
flat; every foot is a node's single foot).

### R7. Determinism

The rule runs in the serial phase and reads only T, N, the mesh and the DEM,
in 14b R1's order. The footed set is written only there. The output stays
bit-identical for any `threads`; 14b's T6 covers it with the rule on.

### R8. Interactions

- **Increment 20 (quality start).** Its Steiner points are nodes snapped from
  circumcentres of triangles with circumradius at least a cell diagonal, and it
  reads no height, so it has no tolerance to scale ε by. **Unchanged.** Not
  measured separately; no triangle under 0.1° remains in any run with it on.
  If a quality node is ever found beside a segment, the cheap fix is a skip
  reason "within `cap` of a constraint", about 8 lines.
- **Carving (14 R6, 19 parked).** Unchanged (R2).
- **Increment 18's scan.** Unchanged; the rule is after the scan. The
  prototype's per-node bitmap is not part of the design; it is the likely
  source of the prototype's 35 ms (not isolated).
- **16b's interior polylines** (roads, rivers) are constraints with neighbours
  on both sides. The rule applies to them unchanged: `split_edge` splits both
  sides, and R2's four-child orientation check covers it. Not measured: no
  interior constraint exists yet.

### R9. Options, report

- **`RefineOptions::constraint_feet`**, default `false` in C++, so the binding
  and the C++ tests keep 20's behaviour unless they ask. The CLI passes `true`.
- **`--no-constraint-feet`** turns it off; the mesh is then bit-identical to
  increment 20's. `elevation_source` adds `constraint feet on` (or `off`),
  ASCII.
- **`RefineOutcome`** gains `feet` and `feet_refused`; `inserted` counts feet
  too. `--stats`' Refinement table gains `feet`.
- No ε knob. The floor and cap are fixed fractions of the cell; C1 is the only
  open parameter.

## Ruled by Ola

**Ola chose C1 (a) slope-scaled ε, C2 (a) the guarantee kept with the fallback,
and C3 (a) no vertex rule, on 2026-09-26.** The options are kept so the reasons
stay on record.

### C1. ε: slope-scaled or fixed?

- **(a) `clamp(tol / G, cell/100, cell/2)`. Recommended.** Ola's direction.
  The foot is enough by construction in the common case (M2: 0 fallbacks in
  116), so no short edges come back. Worst angles 0.40°, 0.11°, 0.40°.
- (b) Half a cell everywhere. Better angles on steep ground (0.33° and 1.10°)
  because it catches nodes up to 5 m from a segment, but 5 fallbacks in 150
  feet, and a fallback edge can be arbitrarily short (1.90 m measured). About 10
  fewer lines.

### C2. What guarantee?

- **(a) Unchanged: within tolerance at every valid node, with the fallback
  (R5). Recommended.** Measured within tolerance on every run. The price is a
  possible short edge where the fallback fires, never seen with C1 (a).
- (b) No fallback: a footed node is never inserted. No short edge can come
  back, but the guarantee weakens to "within tol, except at footed nodes, where
  the error is at most `tol + |∇p|·ε`" with no bound on `|∇p|`. Measured with a
  fixed ε: 1.94 m and 2.16 m at 1 m. With C1 (a) it was 0.32 to 0.53 m, but
  that is a measurement, not a bound. Not recommended.

### C3. Vertex rule as well?

- **(a) No. Recommended.** Measured to do nothing (M1). R2 step 2 already
  covers the ends of segments.
- (b) Also refuse nodes within ε of any off-node vertex, as first briefed.
  About 15 lines for no measured effect.

## Prior art in `legacy/`

```sh
$ grep -rliE 'min_dist|minimum distance|encroach|project.*segment|foot' legacy/
$
```

No file. The legacy triangulated every DEM point with CGAL and had no input
segments to refine against (14b's prior art). **Nothing is carried across.**
`@migration-expert` is not needed.

## Files

| file | what |
|---|---|
| `include/terrain/refinement/refine.hpp` | `constraint_feet` option; ε (R3); the foot branch (R2); footed set; world point and z of an inserted off-node vertex; `feet`, `feet_refused` |
| `include/terrain/mesh/lattice_mesh.hpp` | `split_edge` and `add_vertex` take `MeshVertex` (a `LatticeVertex` still converts) |
| `bindings/core.cpp`, `_core.pyi` | keyword and two fields |
| `src_python/tin_engine/cli.py` | `--no-constraint-feet`; sentence; report |
| `src_python/tin_engine/stats.py` | one Refinement column |
| `ROADMAP.md` | row 20b |

## Tests for `@tester`

**Invariant-critical (mutation testing required): one suite,
`test_refinement_constraint_feet`** (new, C++).

- F1 **tolerance oracle** (R5). A synthetic DEM with a steep plane plus a bump
  beside a long constraint segment at a small irrational angle to the lattice,
  so nodes sit 0.001 to 0.5 cells from it. The oracle recomputes the error at
  every valid node from the output vertices and z, independent of `scan`.
  Within tolerance everywhere.
- F2 **the needle is gone.** A node at 0.0035 cells from the segment whose
  error exceeds tolerance: the output has a foot on the segment, no triangle
  with that node and both segment ends, and the node is not a vertex (the foot
  sufficed). The foot's z equals `vertex_z` at its stored position.
- F3 **constraints.** The union of constrained edges covers the input segment
  (every foot within 1e-9 cells of the input line, found by the test from the
  output); bits and masks on both halves; input vertices unchanged.
- F4 **fallback and termination.** F2's fixture at `tolerance == 0`, where a
  foot almost never brings its node to zero error: the footed node is inserted
  after its foot, the run ends, feet ≤ footed nodes, and every foot was at
  least `floor` from both ends of its segment piece when inserted.
- F5 **off switch** gives output bit-identical to increment 20 on 14b's T12
  fixtures and a domain start.
- Mutants to kill: no fallback (F1 fails, as measured); a second foot allowed
  (F4's feet bound); step 2 removed (F4's termination bound); z of the foot
  from the nearest node instead of bilinear (F1); constrained bit dropped on a
  half (F3); ε without the cap (a disc reaching a second node; F2's node count).

Everything else is property or integration testing:

- 14b's T3 (tolerance) and T6 (determinism for threads 1, 2, 7 and hardware
  concurrency) re-run with the rule on. T6 is what the TSan job runs.
- CLI: the sentence, `--no-constraint-feet`, the `--stats` column.
- T-real, logged, not thresholded except the tolerance: the quarter circle at
  1 m and 10 m with M1's columns. M1's grad rows are the reference.

## LOC estimate

Counted in `CLAUDE.md` §2's unit.

| file | what | est. |
|---|---|---|
| `refinement/refine.hpp` | option, ε, foot branch with checks, footed set, output point, counters | ~70 |
| `mesh/lattice_mesh.hpp` | two signatures | ~2 |
| `bindings/core.cpp`, `_core.pyi` | keyword, two fields | ~8 |
| `cli.py` | flag, sentence, report | ~15 |
| `stats.py` | one column | ~5 |
| | **total** | **~100** |

On the worst overrun seen so far (increment 20: about +50 %), about 150. **Under
700, not split.** C1 (b) is about 10 fewer; C3 (b) about 15 more.

### What landed

132 production lines against ~100 (+32 %, under 700); most in `refine.hpp`
(about 90 against 70). Settled at green: `--no-constraint-feet` needs
`--tolerance` and `--dem`, like `--start-min-angle`; only the first
constrained edge closer than ε is tried; a refusal is counted only when `N` is
then inserted; when `G` is 0, ε is the cap (otherwise 0/0 at tolerance 0 on
flat ground); the footed set is written only in the serial phase.

Ola's quarter circle, text `.vtk` with `--stats`, in the session's
`paraview/runs/20b-constraint-feet/`:

| run | triangles | worst angle | < 0.1° | < 1° | max degree | feet | refine |
|---|---|---|---|---|---|---|---|
| 1 m, feet | 428 217 | 0.3955° | 0 | 12 | 18 | 6 | 0.241 s |
| 1 m, no feet | 428 225 | 0.0117° | 3 | 19 | 18 | 0 | 0.239 s |
| 10 m, either | 31 581 | 0.2873° | 0 | 4 | 16 | 0 | 0.036 s |

Both 1 m runs reach 0.99998 m. The rule costs about 2 ms. At 10 m it never
fires and the mesh equals increment 20's.

Found at green, fixed by `@tester` in `cb7eb96`: increment 20's Q1 pin
(`inserted == 109`) depended on floating-point contraction (92 under
`-ffp-contract=off`); it is now bounds that hold in both modes and still kill
the stale-check and miscount mutants. The golden CLI test needs both passes
off to reproduce pre-20 output.

## Acceptance

- The quarter circle at `--tolerance 1`: no triangle under 0.1°, worst angle
  near 0.4°, achieved max error at most 1 m; at `--tolerance 10` the mesh is
  identical to increment 20's.
- `--no-constraint-feet` reproduces increment 20's mesh bit for bit.
- All gates in `CLAUDE.md` §4 green, including the TSan job; CI green.

## Not in scope

- Nodes near a constraint that is not an edge of their triangle (R1). Not seen.
- An angle bound during DEM refinement (increment 20's C3); the remaining worst
  angles are refinement's own.
- Fans from an input vertex to nodes nearly in line with it (M2's 0.107°).
- The quality pass's own check against constraints (R8).
