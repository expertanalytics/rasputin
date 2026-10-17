# Increment 66: no two vertices written at one point

**Status: designed, design review APPROVED (round 2); R2's test amended
after code review round 1** (`@architect`, 2026-10-17, on `ecdaecd2`, revised
on `a3cb7f40`, R2 amended on `acbf6271`; branch
`worktree-coincident-vertices`, based on increment 65's head; 65 merged as
#278, so master is merged in before the push). Full
path: it changes C++ refine code, so `@perf`'s acceptance run applies. Ola
ruled zero-area slivers bugs (relayed by the main session, 2026-10-17).

**Scope.** Every outline of section 1's probe with a coincident vertex, a
flat triangle or a refused run is in scope, the 4 with flat triangles and
no coincident vertex as written included (section 1, "Flat with no
coincident pair"), by that ruling. Test C2 asserts all 10.

## 1. The defect

Fixture W1 of increment 65 (`docs/increments/65-delaunay-as-written.md`,
section 4; repeated in section 4 below), outline k = 82, constraint feet
**off** (`--no-constraint-feet`; the CLI's default is on), tolerance 0.5:

- The outline's edge from vertex 1 (500094.6, 6599920.9) to vertex 2
  (500088.6, 6599971.9) passes, in exact decimals, through the DEM node
  (col 9, row 4) at (500090.0, 6599960.0); as doubles it misses it by a
  rounding.
- `refine` (phase 1) inserts that node as vertex 7, Inside a triangle, off
  the constraint: phase 1 has no rule for a node a hair off a constrained
  edge.
- `constraint_check_points` on phase 1's edges puts the edge's crossing of
  grid row 4 at col 9.0000000000020552 and its crossing of grid column 9 at
  row 3.9999999999825331 (lattice units), and a midpoint between the two: it
  snaps a crossing to a node only when the node is exactly on the edge
  (`orient_sign == 0`). `refine_strip` inserts the row crossing as vertex 49,
  on the constraint. Vertices 7 and 49 both write as (500090.0, 6599960.0),
  z -3.0 and -2.99999999999589.
- Two output triangles, (49, 7, 11) and (49, 2, 7), have zero area as written
  (`flat_as_written == 2`).

**Why feet on avoids it at this tolerance.** Traced on the probe (section
5): with feet on, phase 1 replaces the node by its foot F on the edge
(`mesh::constraint_foot`), so the edge is split at F (phase 1's edges 1-7 and
7-2, vertex 7 = F, z -2.9999999999949942, the bilinear value at F), and the
strip has no crossing strictly inside a sub-edge near F. So yes: the foot
rule is the reason. It does not hold at tolerance 0: there the node itself
also goes in later (its error at F is about 1e-14, above 0), F and the node
write to one point, and `refine_strip` refuses phase 1's output ("a start
triangle is not counter-clockwise") on k = 82, 1172 and 2922 of 3000
outlines. With the default feet on, the CLI fails on those inputs.

**The wider probe** (section 5): W1 with outlines k = 0..2999, feet on and
off, tolerances 1, 0.5, 0.1, 0.05, 0.01, 0.001, 1e-6 and 0. At every
tolerance above 0, only k = 82 with feet off has coincident vertices or flat
triangles. At tolerance 0 there are more, from both mechanisms, with feet on
and off (section 4, C2).

**Flat with no coincident pair** (k = 1457, 2157 and 2852 at tolerance 0,
feet on and off; printed by the probe's detail mode, section 5). The edge
crosses grid row R and grid column C within about 1e-11 lattice units of
the node (C, R), and neither crossing is that node (the node is not exactly
on the edge). The strip stores the row crossing, at row exactly R, and the
column crossing or the midpoint, at row R plus or minus about 1e-11. As
written, y = y0 - 10 row rounds both to the same double (an ulp of y near
6.6e6 m is 9.3e-10 m, coarser than the 1e-10 m between them), while x stays
2 ulps apart (an ulp of x near 5e5 m is 5.8e-11 m). So the two points are
distinct but share the row's y exactly, and a triangle with a third vertex
on the same grid row is flat. k = 1457: vertices at col 8.99999999998836
and col 9 on row 6, third vertex the node (8, 6). k = 2852: the same on row
7 near col 9, third vertex (8, 7). k = 2157: that on row 7, and three
points on row 3 near col 10 (col 9.99999999998836, its midpoint
9.99999999999418, and 10). R2 snaps every such crossing to the node, so
the near points become one.

## 2. Prior art: legacy and literature

*Literature.* Snap rounding (D. H. Greene and F. F. Yao, "Finite-resolution
computational geometry", FOCS 1986; J. Hershberger, "Stable snap rounding",
*CGTA* 46, 2013) states the discipline: a grid point within the rounding
tolerance of a segment becomes a vertex of that segment. The repository
already applies it in the last stage: L12 of increment 15f
(`docs/increments/15f-edge-strip.md`, L12 and L16: a point Inside but within
the coincidence radius of a constrained edge goes in on it:
`near_constraint`, in `strip_scan.hpp` before this increment, in
`refine.hpp`'s `detail` namespace after it). This increment applies the same rule,
with the same radius, at the two places it is missing: phase 1's split and
the strip generator's crossings. No departure from the method; no novelty
claimed.

*Legacy.* Nothing to carry; the legacy tree (CGAL) has no coincidence or
snapping rule of its own:

```
$ git grep -l -i -E 'coincid|duplicate_vert|snap' legacy-archive -- legacy
(no output, exit 1)
```

*Ideas file.* No entry of `../rasputin_scratch/ideas/ideas.md` is taken up
or overlapped.

## 3. Design

One rule, snap rounding at the coincidence radius, in two places. Neither
alone is enough (probe, section 5). Phase 1's rule alone leaves 4 outlines
with flat triangles at tolerance 0, feet on and off: k = 1457 (flat 1),
2157 (flat 2), 2622 (flat 1, and 1 coincident vertex) and 2852 (flat 1),
all four in C2's list; 1457, 2157 and 2852 are the "flat with no coincident
pair" case of section 1, and in 2622 two strip points write onto node
(8, 6), neither of them the node. The generator's rule alone leaves the 3
refused runs at tolerance 0 with feet on (k = 82, 1172, 2922). Both
together: 0 coincident vertices, 0 flat triangles and 0 refused runs on all
3000 outlines, both feet settings, every tolerance probed (1, 0.5, 0.1,
0.001, 1e-6 and 0 rerun on `a3cb7f40`).

**The radius.** `coincidence_radius(g)` as it is (L16):
max(1e-10, 64 ulp(M)) lattice units, M = max(cols, rows) - 1. Scale assumed:
lattice coordinates up to about 3e4 (a 30 m basin; there 64 ulp is about
2.3e-10). Checked at: W1 (12 x 12 nodes, where it is 1e-10) and the seven
quick-check cases (section 5). The written-point equality of increment 65's
step 2 was the other candidate; rejected, because it depends on the output
frame and the decision is taken on the grid, and because one radius for every
"a hair off a constraint" decision keeps the stages consistent.

**R1. Phase 1 (`refine`, `refine.hpp`'s serial phase).** When the worst node
N of triangle t is Inside t (`r.where == NodeLocation::Inside`) and
`near_constraint(m, t, N, coincidence_radius(g))` names an edge e, and
`strip_fits(m, t, e, N, frame)` holds, N goes in by `split_edge(t, e, N)`
at its own position, with its own z (it is a DEM node). This is tested
**before** the foot search, and when it applies no foot is searched (N is on
the constraint already); this is `point_loop`'s order (L12 before the foot,
`refine_points.hpp`). If `strip_fits` refuses, nothing changes: the foot
search and `split_inside` run as today. The split goes through the existing
`split_edge` path (touched and skip rule, seeds, `legalise_around`).
Frozen edges are excluded by `near_constraint` already. R1 also applies
when t is void (`r.is_void`, a vertex without data): the foot search skips
a void triangle because a foot's z is bilinear over valid cells, but R1
inserts the node itself with its own z, exactly as `split_inside` would,
only on the edge; the hazard it removes is geometric and does not depend on
data. `carved` counts the insertion as today. Not exercised by W1 (no
no-data cells); the prototype counted 0 void firings there.

**R2. The strip generator (`constraint_check_points`,
`constraint_points.hpp`).** A crossing whose rounded node N lies within
`coincidence_radius(g)` of the crossing itself is N: the test is
`std::hypot(at.col - N.col, at.row - N.row) <= radius`, with `at` the
crossing as clamped (Euclidean in (col, row), as `near_constraint`). Today
only a node exactly on the line is (`orient_sign == 0`, which stays as the
first test).

*Why the crossing and not the line (amended after code review round 1).*
The first version tested N's distance to the edge's line. On an edge almost
parallel to a grid line, N can lie within the radius of the line but past
the edge's end, or far along it from the crossing; the snap then moved the
point there (C5 below: an edge from row 3.6875 to 3.9375 got a point at row
4). The two fixes offered were the point test above, or the line test plus
`near_constraint`'s 0 < sigma < 1. The point test is chosen:

- It moves a point by at most the radius. A point's `s` is its crossing's
  t, and the strip's error is taken at `s` (`strip_scan.hpp`'s `along`,
  `seam.hpp`); the sigma test keeps a snap that moves the point by up to half
  a cell along the edge with `s` unchanged (C5 B: the node at s = 0.2 where
  its parameter is 0.5), so its error would be read at the wrong place.
- It is snap rounding's own criterion: N within the radius of a point of
  the segment is within the radius of the segment (a hot pixel, section 2).
  N can lie past an end only when that end is itself within the radius of N,
  the case L14 already handles (a node within ulps of an end); the
  unsnapped crossing is as close to that end.
- It cannot be contracted: one subtraction is zero and the other exact
  (Sterbenz), and `hypot` is a call, so FMA does not change it.
- It snaps nothing the first version left alone, up to rounding. N's
  distance to the line is at most its distance to a point on the line, and
  the stored crossing is on the line to within its rounding (under 4e-12
  lattice units at 3e4), so every crossing it snaps the first version snapped
  too, save a node within that margin above the radius. The prototype's R2
  (section 5) was the line test and fired 0 times on every quick-check case,
  so this one is expected to fire 0 times there too; `@perf`'s acceptance
  hashes check it.
- It still catches every crossing the defect needs snapped. Each such
  crossing writes onto, or into one row with, its node, so it lies within
  ulps of it: section 1's are 2.1e-12 (k = 82, the row crossing), 1.7e-11
  (k = 82, the column crossing) and 1.2e-11 (k = 1457, 2157, 2852) lattice
  units from the node, inside the radius (1e-10 on W1). Assumed, not run
  (no C++ build in this amendment): `@developer`'s green run of C2 and C4
  checks it; C2 failing on any outline goes back to `@architect`.
Two crossings that snap to one node are then one point (the existing
"same position as the entry before" drop), and so is the midpoint between
them. On W1 k = 82's edge, the store goes from 13 points to 11, one more
duplicate, and the point at the node is exactly (9, 4) with z -3.0. R2 also
reaches the frozen seam's check points (`seam.hpp` calls the generator);
same rule, same reason.

**Order of the two.** With R1, the node is an end of two sub-edges, so the
generator has no crossing strictly inside them there; R2 covers the nodes
phase 1 did not insert (the tolerance-0 cases where two crossings near a
node, neither inserted by phase 1, both write at the node).

**Sharing the code.** `coincidence_radius`, `near_constraint` and
`strip_fits` move from `strip_scan.hpp` to `refine.hpp`'s `detail`
namespace (`strip_scan.hpp` and `constraint_points.hpp` both include
`refine.hpp`; the reverse would be a cycle). The distance expression is
not factored out: with R2's point test, `near_constraint` is its only user,
so the `detail::line_distance` of the first version is folded back into
`near_constraint` as the expression it was, and its comment ("one copy, so
near_constraint and the generator cannot be contracted differently") goes
with it. That comment held for clang's default contraction, which fuses
only within one expression, not for GCC's `-ffp-contract=fast`, which may
fuse across statements and differently at each inlined call. No new header.

**Rejected:** dropping a strip point that writes onto an existing vertex.
Vertex 7 would stay a hair off the constraint and lie on it as written, so a
triangle along the constraint would be flat as written instead. Also
rejected: R1 after the foot search (keeps feet-on output byte-identical on
W1), because it leaves the three refused runs at tolerance 0 with feet on.

**Interfaces.** None change: no new option, field, binding or CLI text.
`refine_points` needs nothing: its only caller (`final_check.run`) passes
phase 1's output, which R1 fixes, and a strip, so its own L12 is on; its
strip comes from the generator, which R2 fixes.

**Locality, threads, robustness.** Both rules are local (one triangle and
its neighbour, one crossing). R1 runs in the serial phase, R2 in the
generator; no state outside the call, the output is still independent of the
thread count. The radius is a decade inside the oracles' `kOnEdge` (1e-9),
as L12 ruled.

## 4. Tests (red first, `@tester`)

Added to `tests/cpp/property/prop_refinement_written_delaunay.cpp`, which
holds fixture W1 and its chain (`w1_run`: `refine`, then
`constraint_check_points` on its output, then `refine_strip`; one thread;
feet as the test says on both stages). `w1_run` gains a tolerance parameter
(today fixed at 0.5) and must not `REQUIRE(out.ok())` for C2, which counts a
refused run.

**Names.** That file already holds increment 65's T1-T5, so this
increment's tests are **C1-C5** ("C" for coincident), tagged
`[coincident]`; C5 goes in `tests/cpp/unit/test_refinement_constraint_points.cpp`. Two comments of increment 65 in it become false and
`@tester` rewrites them in the red commit: line 26 ("which is why T1 asserts
per outline and T2 does not assert flat_as_written") and line 287
("flat_as_written is not asserted (2 on k = 82 on macOS arm64; Linux may
round differently)"). Both become: feet-off `flat_as_written` is asserted
by C1 (increment 66); the platform note becomes the measurement below.

**Platforms.** The figures below were measured on macOS arm64 (Apple clang,
`-O2`, FMA contraction on: 253 fused multiply-add instructions in the probe
binary) and again with `-ffp-contract=off` (2 left, the ones in the
predicate library), which is what x86-64 without `-march` compiles to (no
FMA instructions exist there). Every figure of C1-C4 and of the probe was
the same in both builds, on master and with the prototype. Still, `@tester`
confirms from the red CI run that C1, C2 and C4 fail on both CI legs
(ubuntu-latest and macos-latest); a leg where one passes on master goes
back to `@architect` before green.

**Fixture W1** (as increment 65, section 4). Geometry `x_min 500000,
y_max 6600000, dx = dy = 10`, 12 x 12 nodes; z at (row r, col c) =
`float((r * 7 + c * 13) % 11 - 5)`. Outline k, four world vertices in
decimetres from the corner, `h(i, lo, hi) = lo + (k * (2i + 7) * 7919 +
i * 104729) % (hi - lo)`: `(x0 + h(0,50,300)/10, y0 - h(1,700,1050)/10)`,
`(x0 + h(2,800,1050)/10, y0 - h(3,750,1050)/10)`,
`(x0 + h(4,700,1050)/10, y0 - h(5,50,300)/10)`,
`(x0 + h(6,50,400)/10, y0 - h(7,50,300)/10)`; start the fan `{0,1,2},
{0,2,3}`, edges `(i, i+1)` mask 1. Measured (probe, macOS arm64, Apple M1
Max): 3000 outlines, both feet settings, tolerance 0, in 0.94 s at `-O2`.

"Coincident" below: the number of output vertices whose written point
(`out.vertices`) is bit-equal to an earlier vertex's, so a point written
three times counts 2. "Flat": `out.flat_as_written`.

- **C1 (red).** Tolerance 0.5, k = 0..99, feet on and off: 0 coincident and
  `flat_as_written == 0` on every outline. Master: k = 82 feet off, 1
  coincident pair, flat 2. The message names k.
- **C2 (red).** Tolerance 0, feet on and off, k in {82, 1172, 1437, 1457,
  1922, 2137, 2157, 2622, 2852, 2922}: `refine_strip` returns ok, 0
  coincident, `flat_as_written == 0`. Master, as (coincident, flat): feet off
  82 (2, 3), 1172 (1, 3), 1437 (2, 3), 1457 (0, 1), 1922 (2, 3), 2137 (1, 3),
  2157 (0, 2), 2622 (1, 1), 2852 (0, 1), 2922 (1, 3); feet on, refused on
  82, 1172, 2922, and 1437 (2, 3), 1457 (0, 1), 1922 (2, 3), 2137 (1, 3),
  2157 (0, 2), 2622 (1, 1), 2852 (0, 1). The other 2990 outlines of 0..2999
  are 0 on master and stay 0 (probed, not asserted). The message names
  each failing k with its two figures.
- **C3 (red).** R1, the node and not its foot: k = 82, tolerance 0.5, feet
  on. Phase 1's output (`refine` alone) has a vertex at world
  (500090.0, 6599960.0) with `z == -3.0` exactly (the node's value) that is
  an end of two of `out.edges`. Master: that vertex is the foot, z
  -2.9999999999949942.
- **C4 (red).** R2 alone, no refine: `constraint_check_points` on W1's DEM
  with the two vertices (500094.6, 6599920.9), (500088.6, 6599971.9) and the
  edge (0, 1): exactly one point within 1e-6 lattice units of (9, 4), and it
  is `at == (9.0, 4.0)` bit for bit, `z == -3.0`; `size() == 11`,
  `duplicates() == 1`. Master: three points there (col 9.0000000000020552
  row 4; col 9.0000000000010267 row 3.9999999999912665, the midpoint; col 9
  row 3.9999999999825331), `size() == 13`, `duplicates() == 0`.
- **C5 (red against `6441ed57`, the first R2).** A passes on master, B not (master leaves row 4's crossing unsnapped). R2
  snaps only within the radius of the crossing (code review round 1). In
  `test_refinement_constraint_points.cpp`, `unit_frame(16, 8)` and
  `dem_on(g)` as CC7 uses them, so lattice (col, row) equals the given
  doubles; one edge (0, 1); e = 2^-35 (about 2.9e-11). Every figure below is
  exact in binary, so compare with `==`; measured by a Python model of the
  generator (`acbf6271`'s expressions) and checked against FMA: the only
  inexact step, B's t = 0.2, rounds the same fused or not.
  - **A, past the end** (@reviewer's edge, rebuilt from its description,
    since the review recorded no coordinates): P0 `unit(5 - e, 3.6875)`, P1
    `unit(5 + e, 3.9375)`. Premise, asserted with `REQUIRE`: `orient_sign`
    of node (5, 4) is not 0, and its distance to the line, written out in
    the test, is 4.366e-11, below `coincidence_radius(g)` (1e-10); the
    node's sigma is 1.25. One crossing, column 5 at t = 0.5, row 3.8125.
    Expected: `size() == 3`, `duplicates() == 0`, `pts[1].at ==
    MeshVertex{5.0, 3.8125}`, `pts[1].s == 0.5`, `pts[2].at.row == 3.875`,
    and no point's row above 3.9375 (P1's). First R2: `pts[1].at == (5, 4)`
    and `pts[2].at.row == 3.96875`, both past P1.
  - **B, inside the edge, far from the crossing**: P0 `unit(5 - e, 3.6875)`,
    P1 `unit(5 + 4e, 4.3125)`. Column 5 is crossed at t = 0.2, row 3.8125
    (0.1875 from node (5, 4)); row 4 at t = 0.5, col 5 + 1.5e (4.366e-11
    from the node). Expected: `size() == 5`, `duplicates() == 0`, points at
    (5 - e/2, 3.75) s 0.1, (5, 3.8125) s 0.2, (5, 3.90625) s 0.35, (5, 4)
    s 0.5 and (5 + 2e, 4.15625) s 0.75 (the last midpoint is between the
    node and P1). First R2, and the sigma variant:
    `size() == 3`, `duplicates() == 1`, the node at s = 0.2.

Invariant-critical suite (mutation targets): R1 removed (killed by C2 feet
on, C3); R1 moved after the foot search (C2 feet on, C3); R2 removed (C2,
C4); R2 on the line distance, the first version (C5 A and B); R2 on the line
distance with the sigma test, the other fix offered (C5 B). C1 kills neither
single-rule removal (each rule alone fixes k = 82 at 0.5); it pins the
reported defect.

## 5. Acceptance (`@perf`) and the probe

**The probe** (2026-10-17, Apple M1 Max, a throwaway copy of `ecdaecd2` in a
scratch directory, removed afterwards; R1 and R2 behind run-time switches,
the off setting being the tree as is):

- W1, k = 0..2999 (section 4's figures). R1 fires 6 times in the 3000
  outlines at every tolerance from 0 to 0.5, feet on or off (4 times at 1);
  `strip_fits` refused none; none on a void triangle (W1 has none). R2
  fires 14 times at 0.5 and at 0 with R1 on (26 with feet off and 17 with
  feet on, R1 off). Rerun on 2026-10-17 on `a3cb7f40` with a fresh
  prototype: the same counts.
- Every quick-check case (`docs/benchmarks/quick/cases.toml`: tile, quarter,
  numedalslagen, lagan, geilo-al-ramp, romsdal-slope, velhas), feet on and
  off, through the CLI: R1 and R2 each fired 0 times, and every output
  `.vtk` (`--binary`) has the same sha256 with both rules as without. The
  probe's build without the rules reproduces the worktree's own `_core` on
  quarter and velhas (same sha256).

The prototype's R2 was the line test (the first version). The amended
point test snaps a subset of what it snapped, up to the crossing's rounding
(section 3, R2), so the
quick-check figures carry over; W1's carry over where every firing was
within the radius of its crossing, which C2 and C4 check.

**Reproducing it.** The probe driver is kept as
`docs/increments/66-probes/probe.cpp` (build lines in its head); it prints
the coincident, flat and refused figures per outline, and with a fifth
argument each flat triangle's vertices (section 1's "Flat with no
coincident pair"). Built against the tree as is, it gives the master
figures. The shipped code has no switch for either rule, so the
single-rule figures (section 3) are reproducible only by building the
probe against the mutants "R1 removed" and "R2 removed" (section 4); the
firing counts above came from counters in the throwaway prototype and are
not reproducible from the tree.

**Timed estimate** (2026-10-17, Apple M1 Max, AC power before, between and
after; `tools/bench_quick.py run`, default flags, all eight quick-check
cases). The tree `a3cb7f40` was timed first and stored as a temporary
baseline (`--save-baseline`, the stored file restored from git afterwards),
then the prototype (a scratch copy with R1 and R2 on, removed afterwards)
was judged against it: verdict **NO CHANGE**, no case slower beyond the
tool's band, no `MESH CHANGED` line, so all eight outputs are
byte-identical. The tree's medians, total seconds: tile 0.421, tile one
thread 0.733, quarter 0.372, numedalslagen 9.008, lagan 14.742,
geilo-al-ramp 0.340, romsdal-slope 3.627, velhas 0.127. Expected cost of
the added steps: none measurable on these cases.

`tools/bench.py run` on the default DEM, both domains, against the previous
increment's merge commit (back to back with `--tree` on `NO BASELINE`),
power state recorded, evidence in
`docs/benchmarks/<date>/66-coincident-written-vertices/`.

- **Expected hash: unchanged on both `tile` and `quarter`** (the probe:
  neither rule changes either output, feet on or off). A changed hash is a
  defect and goes back to `@developer`.
- **Time:** within `bench.py`'s threshold (5 %) at every thread count. R1
  adds one `near_constraint` test (three edges, a few flops) per Inside
  insertion in the serial phase; R2 one distance per crossing whose rounded
  node is not on the line.

## 6. Size

About 16 net counted lines (the prototype's diff, less its switches and
counters): R1 about 9 in `refine.hpp` (the test, and a flag that R1 fired
guarding the foot search; not `!edge`, which would also skip the foot
search for a node found on an edge, `r.where != Inside`, and change
today's output), R2 about 4 in `constraint_points.hpp` (no `line_distance`: folded back,
section 3), the move of three functions from `strip_scan.hpp` to `refine.hpp`
about 0 net. Tests are not counted.

## 7. Left out

- Feet-on output changes where R1 fires (6 insertions in W1's 3000
  outlines at 0.5, one of them on k = 82): the node goes in instead of its foot. Not pinned
  byte-for-byte; C3 pins the node.
- A counter for R1 or R2 in the outcome or the run record: not needed by any
  test, and Ola's run record has no reader for it.

## Review

- **Design review, round 1, 2026-10-17.** Commit `0ea809cb`, 0 counted lines, @reviewer: CHANGES REQUESTED. No timed estimate (`tools/bench_quick.py` on the quick-check cases); the 4 outlines R1 alone leaves flat are not named, nor how 1457, 2157 and 2852 come out flat with no coincident pair; T1 and T4 shown to fail only on macOS arm64 (FMA rounding); test names T1-T4 collide with increment 65's, and two comments there (lines 26, 287) become false. Suggestions: R1 on a void triangle; the probe figures are reproducible only through the mutants.
- **Design review, round 2, 2026-10-17.** Commit `95c8b19b`, 0 counted lines, @reviewer: APPROVED. All four round-1 items fixed; the R1-fired flag (not `!edge`) keeps today's output; the probe under docs/ builds and matches C1 and C2's master figures. Probe pointer to section 3 fixed by the main session; probe includes and argc left as they are (a docs probe, not built by CI).
- **Code review, round 1, 2026-10-17.** Range `ecdaecd2..6441ed57`, 10 counted lines, @reviewer: CHANGES REQUESTED. R2 tests the rounded node's distance to the edge's line, not that its projection falls inside the edge (`near_constraint`'s 0 < sigma < 1): on an edge nearly parallel to a grid line the snapped point lands past the edge's end (probe: edge from row 3.7 to 3.95, crossing stored at row 4). Fix: snap only within the radius of the crossing, or add the sigma test; @architect amends R2, @tester adds the edge, @developer fixes. Held: the `!edge` flag mutant killed by the pytest hash test is enough; CC7's lost case cannot occur outside the radius; `<=` is right. Wordings: section 2 still places `near_constraint` in `strip_scan.hpp`; the `line_distance` comment holds for clang's default contraction, not GCC's `-ffp-contract=fast`.
- **R2 fix and master merge, 2026-10-17.** `80d3e2ce` (R2 snaps only within the radius of the crossing; C5 red before, green after; ctest 1042/1042 in Release, `-ffp-contract=off` and asan/ubsan). Correction to that commit's message: its pytest was 7843 passed and 50 errors, all from Ola's user `sources.toml` in a format the branch did not yet know (137 passed with an empty config folder). Master merge `d7d7e163` (#277-#282, no conflicts): pytest 8100 passed, ctest 1042/1042. 6 net counted lines.
- **Code review, round 2, 2026-10-17.** Range `acbf6271..edd6f02d`, 6 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session (the subset claim holds up to the crossing's rounding). Round 1's finding closed by C5; no fresh mutant needed ("R2 removed" still killed by C2, C4, CC7; the first R2 and the sigma variant killed by C5). `@perf` next.
