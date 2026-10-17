# Increment 65: the refined mesh is Delaunay where it is written

**Status: designed, revised to design review round 1** (`@architect`,
2026-10-17, on `db12c3e1`). Full path: it changes C++ refine code, so
`@perf`'s acceptance run applies. One public PR, branch
`worktree-refine-world`.

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
19 non-Delaunay edges in 12 outlines with constraint feet on (k = 0, 3, 10,
36, 43, 52, 53, 66, 74, 77, 84, 93; re-run 2026-10-17, same figures), 0 with
them off; `refine` (phase 1) alone has 0 with feet on (same 100, CDT start).

A second finding, measured on the same run: with feet **off**
(`--no-constraint-feet`; the CLI's default is on), outline 82 writes two
vertices at one point, (500090.0, 6599960.0), so two output triangles have
zero area as written. Ola ruled zero-area slivers bugs; the fix is
increment 66 (`docs/increments/66-coincident-written-vertices.md`). This
increment only keeps such a triangle out of its swap pass and counts it.

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
the written doubles. The filter of step 2 (skip a quad whose written corners
are an exact translate of its grid corners) rests only on the circle and
orientation tests being exact and translation-invariant. No novelty is
claimed.

*Legacy.* Nothing to carry; the legacy tree (CGAL) has no circle test of its
own:

```
$ git grep -l -i -E 'incircle|side_of_oriented_circle' legacy-archive -- legacy
(no output, exit 1)
```

*Ideas file.* No entry of `../rasputin_scratch/ideas/ideas.md` is taken up
or overlapped.

## 3. Design

### Choice: one check on the written points, inside the loop

Weighed against the alternative, swapping on written coordinates
throughout:

- **Mesh hash (K1 of increment 23: a run that needs no fix is
  bit-identical).** The check swaps only where the exact test calls an edge
  illegal as written, which is what `tools/bench.py`'s Delaunay check counts.
  A mesh with no such edge builds nothing and keeps its bytes. Swapping on
  written coordinates throughout changes every decision that touches an
  off-node vertex, so every mesh with an outline changes, correct or not.
- **Cost.** Measured on a throwaway prototype (section 5): the check of
  steps 1 to 4 takes at most 0.5 % of any quick-check case's run. A copy of
  the whole mesh on every run, as the first version of this design had it,
  cost 13 to 20 % on three of the seven cases. Throughout would also turn off
  the integer fast path (`lattice_incircle`, increment 21b) and need a
  per-vertex written-point lookup in `mesh/lawson.hpp`, shared with
  `refine.hpp` and `quality.hpp`: a far wider change.

**Why not `0d422717` as it is.** It swaps after the loop has measured the
error, so `max_error`, `max_error_near`, the slope figures and the strip
error describe a triangulation that is not the one written, and a swapped
quad holding a check point can exceed the tolerance unseen. And it skips the
whole pass when one triangle is flat or turned over as written (outline 82
has two). Here the check moves into the loop, and a flat triangle is left
out of it, not a reason to skip it.

### The steps

In `point_loop`, where the loop now breaks (`if (!any) break;`):

1. **Written points.** One lambda gives vertex `i`'s written point `w`
   (start vertex as given, inserted vertex
   `(x_min + col dx, y_max - row dy)`); the output loop after the loop uses
   the same one, so the check and the file cannot disagree. For the exact
   tests the written point is the mesh vertex `MeshVertex{w.x, -w.y}`, so
   `frame()` is `(w.x, w.y)`.
2. **Faithful vertices.** Vertex `i` is faithful when its grid frame point
   `f = frame.at(v)` (the point the loop decides on) is exact,
   `fma(col, dx, -f.x) == 0` and `fma(row, dy, f.y) == 0`, and the written
   point is exactly `f` translated, `x_min + f.x == w.x` and
   `y_max + f.y == w.y` with a zero two-sum error. A triangle with three
   faithful corners is, as written, an exact scaling by `dx, dy > 0` and an
   exact translation of its grid triangle, so it is counter-clockwise as
   written; an edge whose quad has four faithful corners gets, as written,
   the circle answer the loop gave it (the exact tests are invariant under an
   exact translation), and the loop leaves every edge legal in its frame, so
   it is legal. Only a triangle with an unfaithful
   corner is a candidate. Every node of a grid whose `col dx` and
   `x_min + col dx` (and the same for rows) are exact is faithful (the tile: 0 unfaithful of 236,724 vertices).
3. **Flat triangles.** A candidate whose exact orientation on the written
   points is not counter-clockwise is flat. It is counted in
   `out.flat_as_written` (the last check's count) and left out of steps 4
   and 5; it keeps its slot and its corners, and its edges are never
   swapped.
4. **Check.** Every unconstrained edge with a neighbour, of a non-flat
   candidate, whose neighbour is not flat, is tested once (from the lower
   slot when both are candidates) with `detail::quad_flips<DefaultKernel>` on
   the written points under `LatticeFrame{1.0, 1.0}` (integer path off, so
   the exact kernel on the written doubles). None illegal: break. Nothing is
   built and the mesh is the loop's.
5. **Pass** (only when step 4 found an illegal edge). Build a `LatticeMesh`
   on the written points of step 1 (vertex `MeshVertex{x, -y}`, so
   `frame()` is `(x, y)`) from the non-flat triangles in slot
   order, with their constrained bits and masks, keeping each copy slot's
   loop slot; run `legalise_all<DefaultKernel>` under
   `LatticeFrame{1.0, 1.0}`. Flips add to `out.flips` and to
   `out.flips_as_written`.
6. **Back to the grid.** Each copy slot's triangle, bits and masks go back
   into its own loop slot; flat triangles stay as they were. Rebuild the
   loop mesh from its own grid vertices (`LatticeMesh::build`) and set the
   frozen mask again. If `build` refuses (a swapped triangle not
   counter-clockwise on the grid), keep the mesh before the pass, store step
   4's count in `out.illegal_as_written`, and break.
7. **Rescan.** Mark every slot the pass wrote as touched and written,
   `rebuild_active`, and continue the loop: it rescans those triangles and
   inserts where they exceed the tolerance; its next end runs steps 1 to 4
   again.

Terminates: a round after a pass either inserts a point not yet a vertex
(a finite set) or ends with a check that finds no illegal edge, or breaks at
step 6. The check's time goes to `split_seconds` ("edge strip: split + flip
(serial)").

Constraints never move (constrained edges are never swapped), so `subs`,
the frozen edges and the feet stay valid. Thread-safe: the check is serial
in the serial phase, no state outside the call.

### Interfaces

- `PointRefineOutcome` gains three `std::size_t` fields:
  `flat_as_written` (triangles flat or turned over as written, left out of
  the last check), `flips_as_written` (the passes' flips, a subset of
  `flips`) and `illegal_as_written` (edges found illegal as written and left
  so, step 6; 0 otherwise).
- Bound read-only, all three in `_core.pyi`. The CLI records two in the run
  record from both `final` paths (`run_record.py`'s `WORDING`):
  `triangles_flat_as_written`, "Triangles flat or turned over at the written
  coordinates; their edges were left out of the last edge-swap check", and
  `edges_not_delaunay_as_written`, "Edges that fail the empty-circle test at
  the written coordinates and could not be swapped".

Nothing else changes: no new option, no change to `refine`, `lawson.hpp` or
`lattice_mesh.hpp`.

## 4. Tests (red first, `@tester`)

New file `tests/cpp/property/prop_refinement_written_delaunay.cpp`,
registered with `add_terrain_backend_test` in `tests/cpp/CMakeLists.txt`.
Its own oracle, no margin: for each interior edge not in `out.edges` whose
two triangles both have positive exact orientation on `out.vertices`
(an edge of a flat triangle is not counted), orient both counter-clockwise
and count the edge when `DefaultKernel::incircle` on `out.vertices` says
Inside. (The two existing Delaunay oracles decide on recovered grid
coordinates and one has a margin, so they cannot see this.)

**Fixture W1.** Geometry `x_min 500000, y_max 6600000, dx = dy = 10`, 12 x 12
nodes; z at (row r, col c) = `float((r * 7 + c * 13) % 11 - 5)` (integers:
the same on every platform). `x0 = x_min = 500000`, `y0 = y_max = 6600000`.
Outline k = 0..99, four world vertices in decimetres from the corner,
`h(i, lo, hi) = lo + (k * (2i + 7) * 7919 + i * 104729) % (hi - lo)`:
`(x0 + h(0,50,300)/10, y0 - h(1,700,1050)/10)`,
`(x0 + h(2,800,1050)/10, y0 - h(3,750,1050)/10)`,
`(x0 + h(4,700,1050)/10, y0 - h(5,50,300)/10)`,
`(x0 + h(6,50,400)/10, y0 - h(7,50,300)/10)`;
start = the fan `{0,1,2}, {0,2,3}` on those points (not `quality_fixtures::fan`,
which goes through `world()`), edges `(i, i+1)` mask 1; `refine` with feet as
the test says, tolerance 0.5, threads 1; `constraint_check_points`;
`refine_strip`, same feet and tolerance. Measured on master (Python
bindings, macOS arm64): at most 86 output vertices, 0.3 s for all 100.

- **T1 (red).** Feet on: on every outline the oracle counts 0,
  `flat_as_written == 0` and `illegal_as_written == 0`. Master: 19 edges in
  12 outlines (first: k = 0, 3, 10, 36, 43, 52, 53, 66). Measured on the
  prototype: 0 edges and 0 flat triangles on all 100; step 4 finds exactly
  the 19 edges, in the same 12 outlines. The message names k. Linux may
  round differently; the aggregate must be > 0 on master on both CI legs,
  which `@tester` confirms from the red CI run.
- **T2.** Feet off: on every outline `flips_as_written == 0`,
  `illegal_as_written == 0` and the oracle counts 0. Measured on the
  prototype: 0 edges on all 100, `flat_as_written` 2 on k = 82 and 0
  elsewhere (macOS arm64; not asserted, Linux may differ). This pins K1:
  where nothing was wrong, the check changes nothing.
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
  `out.max_error` is at least the largest node error found. Guards step 7.
- **T5 (red).** A flat triangle does not stop the pass: W1's outlines
  k = 0 and k = 10 with feet on, plus a second piece in the cell at the
  grid's corner: A (500000.515625, 6599999.4375), B (500000.53125,
  6599999.421875), C (500000.546875, 6599999.40625) (multiples of 1/64,
  collinear as written, orientation +1.08e-20 on the grid), D (500001.5,
  6599999.75), triangles (A, B, C) and (A, C, D), all five edges constrained,
  mask 1. Measured on the prototype: `flat_as_written == 1`,
  `flips_as_written` 2 and 3, the oracle 0; master: 2 and 3 edges.
  Asserted: `flat_as_written == 1`, `flips_as_written > 0`, the oracle 0.

Invariant-critical suite (mutation targets): T1 against the check removed;
T3 against the count removed; T5 against "skip the check when any triangle
is flat" (`0d422717`'s rule); T2 against step 2's filter removed is not a
target (it changes time, not output); T4 against step 7 removed (break
after the pass). If no outline makes T4 fail under that mutant, `@tester`
records it.

## 5. Acceptance (`@perf`)

**The time estimate, measured** (2026-10-17, Apple M1 Max, AC power,
Release build of `37bea375` plus a throwaway prototype of steps 1 to 6 in a
scratch copy, default flags, every quick-check case
(`docs/benchmarks/quick/cases.toml`), three runs each, the timer around steps
1 to 4 inside `point_loop`):

| case | triangles | unfaithful vertices | tested edges | check (s) | run total (s, baseline) | share |
|---|---|---|---|---|---|---|
| tile | 472,374 | 0 | 0 | 0.0014-0.0015 | 0.42 | 0.4 % |
| quarter | 428,218 | 231 | 557 | 0.0013-0.0015 | 0.36 | 0.4 % |
| numedalslagen | 653,039 | 63,345 | 336,895 | 0.018-0.019 | 8.85 | 0.2 % |
| lagan | 463,342 | 114,945 | 343,765 | 0.017-0.019 | 14.5 | 0.1 % |
| geilo-al-ramp | 64,162 | 106 | 537 | 0.0002 | 0.33 | 0.1 % |
| romsdal-slope | 3,096,202 | 0 | 0 | 0.010 | 3.59 | 0.3 % |
| velhas | 99,236 | 34 | 356 | 0.0003-0.0004 | 0.12 | 0.3 % |

Every case: 0 flat triangles, 0 illegal edges, so no pass is built and the
mesh is unchanged. Run totals are `baseline-apple-m1-max-ac.json`'s medians.
The first version's whole-mesh copy (flat check, `build`,
`legalise_all`), timed in the same runs: 0.076 s on tile (18 %), 0.068 s
on quarter (19 %), 0.70 s on romsdal-slope (20 %), 0.012 s on velhas
(10 %). The prototype is gone; the figures are its stderr lines.

`tools/bench.py run` on the default DEM, both domains, against the previous
increment's merge commit (back to back with `--tree` on `NO BASELINE`),
power state recorded, evidence in `docs/benchmarks/<date>/65-delaunay-as-written/`.

- **Expected hash: unchanged on both `tile` and `quarter`.** Their last
  recorded quality (2026-10-10, `slope-tol/bench-head`) is 0 violations as
  written, and the prototype found no illegal edge on either, so no pass
  runs. `@perf` records each domain's mesh sha256 for base and head side by
  side, and the Delaunay count. A changed hash with a base count of 0 is a
  defect and goes back to `@developer`; a changed hash with a base count
  above 0 is the fix, and is recorded so.
- **Time:** within `bench.py`'s threshold (5 %) at every thread count; the
  estimate above is under 0.5 %.

## 6. Size

About 55 net counted lines: the faithful filter, flat count, check, pass and
rebuild in `refine_points.hpp` (about 45, replacing `0d422717`'s 24), three
fields, three bindings, three stub lines, two CLI lines, two wording lines.
Tests are not counted.

## 7. Left out

- The coincident written vertices of outline 82 (section 1): counted and
  kept out of the check here, fixed in increment 66.
- `refine` (phase 1) writes the same way but is never the last stage; 0
  violations measured on the 100 outlines.
- Branch `worktree-23c` drops `0d422717` and takes this from master.

## Review

- **Design review, round 1, 2026-10-17.** Commit `37bea375`, @reviewer: CHANGES REQUESTED. Asked for a timed estimate of the pass on the quick check's cases, a flat triangle left out rather than stopping the pass, `flat_as_written == 0` in T1, and the fixture's `x0`, `y0` and step 1's vertex defined; answered in this revision.
