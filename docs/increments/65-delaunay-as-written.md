# Increment 65: the refined mesh is Delaunay where it is written

**Status: designed** (`@architect`, 2026-10-17, on `db12c3e1`). Full path:
it changes C++ refine code, so `@perf`'s acceptance run applies. One public
PR, branch `worktree-refine-world`.

## 1. The defect

`detail::point_loop` (`include/terrain/refinement/refine_points.hpp`), the
last stage of `refine_strip` and `refine_points` and so of every `mesh` run,
decides every edge swap with the circle test on grid coordinates
`(col * dx, -(row * dy))`. It writes an inserted vertex at
`(x_min + col * dx, y_max - row * dy)`, and a start vertex as given; both
differ from the grid point by rounding. Where four corners are nearly on one
circle, the swap is right on the grid and wrong as written. A constraint foot
does it: node N, its foot F on an outline, the outline's crossing A of the
grid line through the node Q above N make the quad N, F, A, Q with right
angles at F and at Q, so it is exactly cocircular, and rounding picks the
side. Found on branch `worktree-23c` (its union test, 4 edges in one piece);
`@reviewer` judged it to be on master too, and it is (fixture W1 below).

Measured on master (the main checkout's `_core`, built 2026-10-16; refine
sources unchanged since `952bd1ee`), with the Python bindings and an exact
rational circle test on the written vertices, over W1's 100 outlines:
19 non-Delaunay edges in 12 outlines with constraint feet on, 0 with them
off; `refine` (phase 1) alone has 0 with feet on (same 100, CDT start).

A second finding, measured on the same run: with feet **off**, outline 82
writes two vertices at one point, (500090.0, 6599960.0): a DEM node `refine`
inserted, and a strip point the loop holds at a grid position a rounding
away from it. Two output triangles are flat as written. Out of scope here
(it needs its own design); this increment reports it (T3), and the main
session carries it as an idea.

## 2. Prior art: legacy and literature

*Literature.* Lawson's flip algorithm (C. L. Lawson, "Software for C1
surface interpolation", in *Mathematical Software III*, 1977): swapping every
illegal edge terminates and gives the Delaunay triangulation; with
constrained edges never swapped, an all-locally-Delaunay triangulation is the
constrained Delaunay one (D. T. Lee and A. K. Lin, "Generalized Delaunay triangulation for planar graphs", *DCG* 1, 1986; L. P. Chew,
"Constrained Delaunay triangulations", *Algorithmica* 4, 1989). Exact
predicates (J. R. Shewchuk, "Adaptive precision floating-point arithmetic and
fast robust geometric predicates", *DCG* 18, 1997) are exact **on the doubles
they are given**; the CGAL kernel design "exact predicates, inexact
constructions" states the same discipline for constructed points such as a
foot: round the point first, then decide on the stored value. The loop
departs from that: it decides on a transform of the stored point. This
increment restores the method as published: the last swap pass decides on
the written doubles. No novelty is claimed.

*Legacy.* Nothing to carry; the legacy tree (CGAL) has no circle test of its
own:

```
$ git grep -l -i -E 'incircle|side_of_oriented_circle' legacy-archive -- legacy
(no output, exit 1)
```

## 3. Design

### Choice: one pass on the written points, inside the loop

Weighed against the alternative, swapping on written coordinates
throughout:

- **Mesh hash (K1 of increment 23: a run that needs no fix is
  bit-identical).** The pass swaps exactly the edges the exact test calls
  illegal as written, which is what `tools/bench.py`'s Delaunay check counts.
  A mesh with no violation gets no swap and the same bytes. Swapping on
  written coordinates throughout changes every decision that touches an
  off-node vertex, so every mesh with an outline changes, correct or not.
- **Cost.** The pass is one adjacency build and one test of every edge, about
  the start's own `legalise_all`, once per run (twice where it swaps). Throughout
  would turn off the integer fast path (`lattice_incircle`, increment 21b)
  wherever `x_min + col * dx` rounds at a node (every degree grid), and
  needs a per-vertex written-point lookup in `mesh/lawson.hpp`, shared with
  `refine.hpp` and `quality.hpp`: a far wider change.

**Why not `0d422717` as it is.** It swaps after the loop has measured the
error, so `max_error`, `max_error_near`, the slope figures and the strip
error describe a triangulation that is not the one written, and a swapped
quad holding a check point can exceed the tolerance unseen. And it skips the
pass silently when a triangle is flat or turned over as written (that case
occurs: outline 82 above). The pass therefore moves into the loop, and the
skip is counted.

### The steps

In `point_loop`, where the loop now breaks (`if (!any) break;`):

1. **Written points.** One lambda gives vertex `i`'s written point (start
   vertex as given, inserted vertex `x_min + col dx, y_max - row dy`); the
   output loop after the loop uses the same one, so the pass and the file
   cannot disagree.
2. **Flat check.** Count the triangles whose exact orientation on the
   written points is not counter-clockwise. If any: store the count in
   `out.flat_as_written`, skip the pass, break (the mesh is as the loop
   made it).
3. **Pass.** Build a `LatticeMesh` on the written points as `(x, -y)`, with
   the loop mesh's triangles, constrained bits and masks per slot, and run
   `legalise_all<DefaultKernel>` under `LatticeFrame{1.0, 1.0}` (integer path
   off, so every test is the exact kernel on the written doubles). Flips add
   to `out.flips` and to `out.flips_as_written`. Zero flips: break.
4. **Back to the grid.** Rebuild the loop mesh from its own grid vertices
   and the pass's triangles, bits and masks (`LatticeMesh::build`, slot `t`
   stays slot `t`), and set the frozen mask again. If `build` refuses (a
   swapped triangle not counter-clockwise on the grid), count its triangles
   in `out.flat_as_written` and break with the mesh before the pass.
5. **Rescan.** Mark every slot whose triangle changed as touched and written,
   `rebuild_active`, and continue the loop: it rescans those triangles and
   inserts where they exceed the tolerance; its next end runs the pass again.

Terminates: a round after a pass either inserts a point not yet a vertex
(a finite set) or ends with a pass that flips nothing (the triangulation is
already Delaunay as written) or breaks at step 2 or 4. The pass's time goes
to `split_seconds` ("edge strip: split + flip (serial)").

Constraints never move (constrained edges are never swapped), so `subs`,
the frozen edges and the feet stay valid. Thread-safe: the pass is serial
in the serial phase, no state outside the call.

### Interfaces

- `PointRefineOutcome` gains `std::size_t flat_as_written` (triangles that
  stopped the pass, steps 2 and 4) and `std::size_t flips_as_written` (the
  passes' flips, a subset of `flips`).
- Bound read-only, both in `_core.pyi`. The CLI records
  `flat_as_written` in the run record as `triangles_flat_as_written`,
  wording "Triangles flat or turned over at the written coordinates; where
  any is, the last edge-swap check was skipped" (`run_record.py`'s
  `WORDING`), from both `final` paths.

Nothing else changes: no new option, no change to `refine`, `lawson.hpp` or
`lattice_mesh.hpp`.

## 4. Tests (red first, `@tester`)

New file `tests/cpp/property/prop_refinement_written_delaunay.cpp`,
registered with `add_terrain_backend_test` in `tests/cpp/CMakeLists.txt`.
Its own oracle, no margin: for each interior edge not in `out.edges`, orient
both triangles counter-clockwise by the exact orientation on `out.vertices`,
and count it when `DefaultKernel::incircle` on `out.vertices` says Inside.
(The two existing Delaunay oracles decide on recovered grid coordinates and
one has a margin, so they cannot see this.)

**Fixture.** Geometry `x_min 500000, y_max 6600000, dx = dy = 10`, 12 x 12
nodes; z at (row r, col c) = `float((r * 7 + c * 13) % 11 - 5)` (integers:
the same on every platform). Outline k = 0..99, four world vertices in
decimetres from the corner, `h(i, lo, hi) = lo + (k * (2i + 7) * 7919 +
i * 104729) % (hi - lo)`:
`(x0 + h(0,50,300)/10, y0 - h(1,700,1050)/10)`,
`(x0 + h(2,800,1050)/10, y0 - h(3,750,1050)/10)`,
`(x0 + h(4,700,1050)/10, y0 - h(5,50,300)/10)`,
`(x0 + h(6,50,400)/10, y0 - h(7,50,300)/10)`;
start = the fan `{0,1,2}, {0,2,3}` on those points (not `quality_fixtures::fan`,
which goes through `world()`), edges `(i, i+1)` mask 1; `refine` with feet as
the test says, tolerance 0.5, threads 1; `constraint_check_points`;
`refine_strip`, same feet and tolerance. Measured on master (Python
bindings, macOS arm64): at most 86 output vertices, 0.3 s for all 100.

- **T1 (red).** Feet on: the oracle counts 0 on every outline. Master: 19
  edges in 12 outlines (first: k = 0, 3, 10, 36, 43, 52, 53, 66). The
  message names k. Linux may round differently; the aggregate must be > 0
  on master on both CI legs, which `@tester` confirms from the red CI run.
- **T2.** Feet off: `flips_as_written == 0` and the oracle counts 0 on every
  outline except those with `flat_as_written > 0` (k = 82 on macOS). This
  pins K1: where nothing was wrong, the pass changes nothing.
- **T3 (red: the field does not exist).** A flat triangle is reported:
  geometry `x_min 0, y_max 30, dx = dy = 3`, 11 x 11 nodes; start vertices
  A (0.5, 0.5), B (0.51, 0.51), C (0.52, 0.52), plus one vertex D on the
  other side of AC from B; triangles (A, B, C) and the one across AC, every
  edge constrained. A, B, C are collinear as written (exact orientation 0)
  and counter-clockwise on the grid (orientation 9.25e-20, computed with
  the loop's `(x - x_min) / dx`). `refine_points` with an empty frozen store:
  `flat_as_written == 1`, triangles as given.
- **T4.** On T1's outlines, every DEM node is within the tolerance of the
  written mesh (`constraint_foot_oracles::node_findings`), and
  `out.max_error` is at least the largest node error found. Guards step 5.

Invariant-critical suite (mutation targets): T1 against the pass removed;
T3 against the count removed; T4 against step 5 removed (break after the
pass). If no outline makes T4 fail under that mutant, `@tester` records it.

## 5. Acceptance (`@perf`)

`tools/bench.py run` on the default DEM, both domains, against the previous
increment's merge commit (back to back with `--tree` on `NO BASELINE`),
power state recorded, evidence in `docs/benchmarks/<date>/65-delaunay-as-written/`.

- **Expected hash: unchanged on both `tile` and `quarter`.** Their last
  recorded quality (2026-10-10, `slope-tol/bench-head`) is 0 violations as
  written, so the pass swaps nothing. `@perf` records each domain's mesh
  sha256 for base and head side by side, and the Delaunay count. A changed
  hash with a base count of 0 is a defect and goes back to `@developer`; a
  changed hash with a base count above 0 is the fix, and is recorded so.
- **Time:** within `bench.py`'s threshold (5 %) at every thread count. Over
  it, the fallback (not built) is to start the pass's stack only from edges
  with an off-node corner, valid only where `x_min + col dx` and
  `y_max - row dy` are exact at every node.

## 6. Size

About 40 net counted lines: the pass, flat count and rebuild in
`refine_points.hpp` (about 32, replacing `0d422717`'s 24), two fields, two
bindings, two stub lines, one CLI line, one wording line. Tests are not
counted.

## 7. Left out

- The coincident written vertices of outline 82 (section 1): reported, not
  fixed.
- `refine` (phase 1) writes the same way but is never the last stage; 0
  violations measured on the 100 outlines.
- Branch `worktree-23c` drops `0d422717` and takes this from master.
