# Increment 20c — fewer slivers: constraint-aware insertion everywhere, and a soft quality criterion

Status: **designed** (`@architect`, 2026-10-06, branch
`worktree-soft-quality` off master `ed12512`), design review rounds 1 to
5 answered and round 6 approved; in three PRs: **20c-1** the foot rule on
every insertion path: **built and approved in code review round 2**
(Ola: "yes, build 20c-1, no push"), 163 net production lines (`python3 tools/count_loc.py c074f900
69f37d1c`). Red tests `a982d531` (`@tester`), their pins ruled (below,
"Pins ruled for 20c-1's red step"); green `8f275220` (`@developer`), its
four open items ruled ("Rulings on 20c-1's green step") and answered
(`b09b84f6` to `34b546c4`); `@perf`'s acceptance `b9e2d482`: every mesh and
quality gate passes, refine time regressed; `@developer`'s exact
performance fix `69f37d1c` (R3); `@perf`'s re-time `c9f62f7c`: the 1 m
benchmark accepted, the split phase still over ruling 5's 2 % on both
catchments; that limit is revised below ("The split-phase limit, judged
after the fix"), 20c-1 meets the revised one, and Ola accepted the
remaining time (question 7, ruling 7). Code review round 1 (`@reviewer`,
on `0346ac63`): changes requested, no mutation record; `@tester`'s
mutation round then ran ("Mutation round for 20c-1": 4 survivors of 16
faults, one killed by the new test `4a64c1ef`), and each survivor is ruled
there; `@tester` then killed the four (`cec2dda4` to `7bc41ebf`, tests
only), so all 16 planted faults are killed; ctest 886 of 886, pytest
5404 passed. Code review round 2 (`@reviewer`, on `74889310`): approved.
**Merged** as PR #207 (master `41bda81a`);
**20c-2** the soft criterion and a split of the
constraint line that blocks the walk to a quality point, after 20c-1
(branch `worktree-soft-quality-2` off master `41bda81a`; Ola: "build
20c-2 and 20c-3"): **built**, 158 net production lines (`python3
tools/count_loc.py 41bda81a 50544838`; "LOC" says why above the ~110
estimate). Red tests `100f5a09` (`@tester`), their pins ruled (below,
"Pins ruled for 20c-2's red step"; R7 gains a slack and R8 needs the feet
on) and amended `5b65941a`; green `d3c939d9` and `96119508`
(`@developer`); mutation-round tests `600b2acb` (`@tester`, three cases
added, tests only); `@perf`'s acceptance `fa136aa7`: every mesh gate
passes, the time limit not met; `@developer`'s speed fix `65ae0792`
(meshes byte-identical); `@perf`'s re-time `5b0afe04`: the whole run and
the 1 m benchmark met, the split phase over 2 % in seconds but not per
split, which Ola ruled is how 20c-2 is judged (ruling 8; "The split-phase
limit for 20c-2, judged per split"); the green step and the fix ruled
("Rulings on 20c-2's green step and speed fix"). Code review round 1
(`@reviewer`, on `8ea0f32d`): changes requested, the `--start-quality-gain`
help promised a line split that `--no-constraint-feet` turns off; fix
round: red `9f54bc95` (`@tester`, the help-text test), green `50544838`
(`@developer`, the help reworded, the atan2-skip comment corrected); the
kill table copied in ("Mutation round for 20c-2"). Code review round 2
(`@reviewer`, on `d4ac395f`): approved; **merged** as PR #214 (master
`4abacf65`);
**20c-3** the repair of land-cover borders that should be one (on at
1 m, ruling 9, Ola's answer to question 9 on 2026-10-08), input coarsening, and Ola's outline rule, on
by default at 5 m (question 6, ruling 6), after 20c-2 (branch
`worktree-soft-quality-3` from 20c-2's head `68cc96b6`; Ola: "build 20c-2
and 20c-3"): **designed, brought up to date 2026-10-08** for ruling 9 and
for what 20c-2 and 30d (PR #215) changed ("Design of PR 20c-3"; M7 the
probe that set the tolerance); **built**: shapely 2.2 `985b4a47`
(`@developer`), red `6459e87c` (`@tester`), green `7334b82b`
(`@developer`; 46 of 49 red tests pass, full suite 3 failures), 331 net
production lines (`python3 tools/count_loc.py 483221d2 7334b82b`, +121 %
on ~150, ruling G5); the green step ruled ("Rulings on 20c-3's green
step": two broken or conflicting tests, the record's switch as text, the
clip's shortened edges, the assumptions); test changes `737e5250`
(`@tester`), green for rulings G3 and G6 (e) `e07d5921` (`@developer`);
mutation round on the outline rule done, test commit `5cc0639c`
(`@tester`, tests only: two survivors of eight faults, M2 and M2b, killed
by two new tests; "Mutation round for 20c-3"), its three pins ruled
("Rulings on 20c-3's mutation round"; OR5's wording corrected);
`@perf`'s timed check and gate runs `16143f4c` and profile `4d52dec0`
(`docs/benchmarks/2026-10-08/20c-3/`): the run +30 % on Lagan and +107 %
on Numedalslagen, the clean-up 8.0 s and 8.8 s against the design's
1.1 s; Numedalslagen's outline-rule gate missed (15 slivers, 11 near the
outline, worst 0.0814°); both ruled ("Rulings on 20c-3's timed check and
gate runs": T1 a defect in ruling G4's cuts, fixed in code; T2 the
rebuild and area measure restricted to what moved; T3 the hotspot,
question 10, ruled 10). Fix round: red `51c24772` (`@tester`, OR9 and the far
part), green T1 `35b7a873` and T2 `8716f062` (`@developer`), 355 net
production lines (`python3 tools/count_loc.py 483221d2 8716f062`);
`@perf`'s re-time `41db7333`: T2 leaves the mesh unchanged, the clean-up
5.92 s and 5.48 s (limits met), the outline-rule rows of the gate met
(the repair rows and three checks of that gate were not run on the fixed
code; "Rulings on 20c-3's code review round 1", R3); the T2 pins
ruled ("Rulings on the fix round"): one changed, the area moved for a
ring with no kept vertex (P1): red `48aca5f7` (`@tester`), green
`1e011cec` (`@developer`), 356 net production lines (`python3
tools/count_loc.py 483221d2 1e011cec`). Code review round 1 (`@reviewer`,
on `1e011cec`): changes requested; ruled ("Rulings on 20c-3's code review
round 1": a hole the rule puts onto the outline is filled, fixed in code;
the record row's wording; the gate's unrun checks, run by `@perf`).
Fix round: red `1f7f8335` (`@tester`, R1 and R4), green `aaee898b`
(`@developer`), 360 net production lines (`python3 tools/count_loc.py
483221d2 aaee898b`); `@perf`'s R3 run `63065a15`: (a) to (d) met, (e)
met once its wording is corrected (R3's outcome, below: the 26.194 m left
are borders between different classes in the source data). Code review
round 2 (`@reviewer`, on `05b0fd7a..07d97ee4`): "Review" below. M8
`d7a3462a` (`@architect`): the repair's tolerance probed at 5 cm, 25 cm
and 1 m on the built code, question 9 asked again with 1 m as its
default; Ola: 1 m (ruling 9's default, 2026-10-08). Red `f6e0f287`
(`@tester`), green `c71ccf57` (`@developer`, three lines changed; 360 net
production lines, unchanged); `@tester`'s docstring line `52f6a5d6`
(tests only). Ola accepted the clean-up's time (ruling 10, question 10,
2026-10-08); the speed-up is its own ROADMAP row. **Built. Next:
`@reviewer`'s code review round 3, on `9a04c3c1..` this record's
commit.** Questions 1 to 4 ruled by Ola, 2026-10-06, and questions 5 to
8, 2026-10-07, ruling 9 (amending 5) 2026-10-07, and its default
(question 9) and question 10 2026-10-08 ("Ola's rulings" below); no
question open. Asked by Ola
("2: yes", 2026-10-06, after calling Lagan's worst angle of 0.000412°
"pretty unacceptable!"). Carries increment 20's C1-C3 rulings
(`docs/increments/20-start-quality.md`, "Ola's rulings").

**In one paragraph.** Most slivers are not made by refinement but by the
two insertion paths that lack 20b's rule (put a point that would land a
hair from a constraint line onto the line instead): the quality start and
the final check against the source DEM. Giving them the rule, and letting
it look one triangle further, removes 58 to 74 % of the slivers for +0.3 to
+0.4 % triangles (Lagan 1 598 to 679, Numedalslagen 1 644 to 433) and
raises the worst angle (Lagan 0.000412° to 0.0024°, Numedalslagen 0.00040°
to 0.00068°). A soft criterion then cuts 13 to 17 % of the triangles; with
it, splitting a constraint line where a quality point lies just beyond it
(Ruppert's response to a segment in the way, with a different trigger)
removes another third of the slivers and lifts
Numedalslagen's worst angle to 0.0086°, still 7 to 11 % below master's
triangle count. What is left is the input's own geometry (millimetre
CORINE segments), which only a tolerance on the input removes: 20c-3's
2 m simplification takes the slivers to 151 and 32, and Ola's outline
rule (land-cover borders within 5 m of the catchment outline put onto it)
to 61 and 4, with Numedalslagen's worst angle at 0.83° (M6). Ola ruled
that a mismatch such as Lagan's 1 cm slit is an error to repair (ruling 9):
20c-3 first joins land-cover borders that come within 1 m of each other
(Ola's default, question 9), which closes the slit, removes the land
cover's millimetre edges, and takes the slivers with the defaults to 26
and 2, the worst angle to 0.73° and 0.83° (M7, M8).

## Prior art: legacy and literature

### Literature

- **Ruppert 1995** (J. Algorithms 18(3), 548-585) and **Chew 1989/1993**:
  Delaunay refinement. A bad triangle gets a Steiner point at its
  circumcentre; a circumcentre that **encroaches** on a constraint segment
  (lies in its diametral circle) is rejected, and the segment is split at its
  midpoint instead. Guarantees a minimum angle of about 20.7° for inputs
  without small angles, with triangle size graded to the local feature size
  (lfs: the distance to the nearest non-incident input feature).
  Increment 20's pass is this with two changes: the Steiner point is the DEM
  node nearest the circumcentre (20 R3), and a rejected point is skipped, not
  turned into a segment split (20 C1 (a)).
- **Shewchuk 2002** ("Delaunay refinement algorithms for triangular mesh
  generation", CGTA 22, 21-74): diametral lenses, concentric-shell splitting
  and the rule for small input angles; and **Triangle's `-Y`/`-YY` switches**
  (https://www.cs.cmu.edu/~quake/triangle.Y.html), which forbid Steiner points
  on boundary or on all segments, with the documented consequence that "the
  resulting mesh may contain triangles of poor aspect ratio". That is exactly
  20 R5's "best effort" under C1 (a), and the reason population 2 below
  exists.
- **Shewchuk 2002, "What is a good linear finite element?"** and
  **Babuška and Aziz 1976** (SIAM J. Numer. Anal. 13): for interpolating a
  surface, the gradient error is driven by the **largest** angle, not the
  smallest. A *needle* (one tiny angle, the others near 90°) interpolates
  slopes well; a *cap* (one angle near 180°) does not. Most of Lagan's
  slivers are caps (M1), so they matter for slope and aspect downstream, not
  only for looks.
- **Üngör 2004 (off-centres)** and **Erten and Üngör 2009** ("Quality
  triangulations with locally optimal Steiner points", SIAM J. Sci. Comput.
  31(3), 2103-2130): the Steiner point is chosen in a region, not fixed at
  the circumcentre, so that the new triangle on the shortest edge just meets
  the target; 30 to 40 % fewer points than circumcentres. Our lattice
  restriction (20 R3) is already a "choose a point in a region" rule; the
  region is the DEM's nodes near the circumcentre.
- **Snap rounding** (Hobby 1999, CGTA 13; Guibas and Marimont 1998) and
  **iterated snap rounding** (Halperin and Packer 2002, CGTA 23(2), 209-225;
  bounded drift: Packer 2008): rounding an arrangement of segments to a grid
  of width h so that, after ISR, every vertex is at least h/2 from every
  non-incident edge. That is a lower bound on the input's local feature size,
  which is what Ruppert's size bound needs. The noder already snap-rounds at
  `--snap-spacing` (increment 5a; default 1 mm, 16 U6); shapely's
  `set_precision` does the same through GEOS.
- **Coverage cleaning** (Davis 2025, "Coverage Cleaning in JTS",
  http://lin-ear-th-inking.blogspot.com/2025/04/coverage-cleaning-in-jts.html;
  ported to GEOS 3.14 as `CoverageCleaner`, PostGIS `ST_CoverageClean`,
  shapely 2.2 `coverage_clean`): a set of polygons meant to be a partition
  is noded with a snapping noder, its faces rebuilt, overlaps given to a
  neighbour, and gaps narrower than a width (the diameter of the largest
  inscribed circle) given to the neighbour with the longest shared border.
  20c-3's repair is this, unchanged, with one tolerance for both the
  snapping and the gap width, and run on the polygons clipped to the read
  region. The snapping noder is not iterated snap rounding (Halperin and
  Packer, above), so it gives no lower bound on the input's feature size,
  and 20c-3 claims none.
- **Garland and Heckbert 1995** (greedy insertion for terrain): the
  tolerance-driven refinement this project builds on (increment 14). It has
  no quality criterion; thin triangles are a known by-product.

**What differs here, and what is claimed.** The design below (a) applies
Ruppert's rejection rule ("a candidate too close to a segment goes on the
segment instead") on the two insertion paths that lack it, with the split at
the candidate's **foot** (orthogonal projection) instead of the segment's
midpoint, as increment 20b already does for refinement. (Shewchuk's 1998
tetrahedral refinement uses an encroaching vertex's orthogonal projection
only to choose which subfacet to split, not as the split point; a midpoint
or circumcentre split is what the guarantees are proved for, and a foot
split gives them up, as 20b's already does.) (b) is input snap rounding at a horizontal tolerance. (c) is a
greedy acceptance test: a quality point goes in only if the worst angle of
the triangles it replaces does not fall; and, with it, Ruppert's response
to an encroached segment (split the segment instead of inserting the
point), again at the point's foot rather than the midpoint, and on a
different trigger: the walk to the point is stopped by the segment, not
the point lying in the segment's diametral circle (R8). **No novelty is claimed.** Searched:
"Delaunay refinement split encroached segment at projection", "Delaunay
refinement cost-benefit Steiner point acceptance", "terrain TIN breaklines
Delaunay refinement minimum angle"; found the works above and nothing that
combines a sup-norm DEM tolerance with these rules, but a novelty claim
would need a deeper search than this design made, so none is made.

### Legacy

```sh
$ git grep -l -i -E "ruppert|encroach|sliver|min_angle|minimum angle|snap_rounding|snap rounding|steiner" legacy-archive -- legacy
$
```

No file (exit status 1). The legacy triangulated every DEM point with CGAL
and had no quality pass (increment 20's prior art). **Nothing is carried
across.**

For 20c-3's repair (ruling 9), the same check on its subject, run
2026-10-08:

```sh
$ git grep -l -i -E "coverage|gap_width|clean_coverage|snapping" legacy-archive -- legacy
$
```

No file (exit status 1). Nothing is carried across.

## The problem, measured

### M1. Which insertion path made each sliver (the investigation's E7)

Lagan (`docs/research/lagan-slivers.md` on branch `worktree-sliver-seam`,
`c83658f1`, corrected at `cd93b2cb`): GLO-30 resampled to a 31 m grid in EPSG:3006, CORINE features,
`--tolerance 10`, 863 897 triangles, 1 598 under 1° (0.185 %), worst
0.000412°. A sliver is a triangle whose smallest plan-view angle is under 1°.

**Method, no C++ change needed.** The vertices of `refine`'s and the final
check's output are in insertion order: the start mesh's `n0` vertices first,
then the quality start's `quality_inserted` (it runs before the first scan),
then refinement's nodes and feet (a foot is the off-node one), then the final
check's points (a point that ends on a constraint edge is a line point,
otherwise a source-DEM node). A driver (`e7.py`, scratchpad, not kept) runs
the real CLI with `refine` and `trim` wrapped to capture those arrays, and
tags every vertex. It reproduces the investigation's figures exactly (863 897
triangles, 1 598 slivers, worst 0.000412°). The sliver's *apex* is the
vertex opposite its longest side.

| apex of the sliver (slivers with a side under 10 cm only in the third row) | slivers | of them, long side a constraint edge |
|---|---|---|
| quality-start node (increment 20) | 639 | 618 |
| input vertex (domain or CORINE) | 516 | 133 |
| any, with a side under 10 cm | 185 | — |
| source-DEM node from the final check | 219 | 196 |
| refinement node (20b's feet active) | 38 | 20 |
| line point (15f's strip) | 1 | — |

- **The quality start is the largest source, not refinement.** Its node is
  the DEM node nearest a circumcentre, up to half a cell diagonal away, and
  nothing stops that node from lying millimetres from a constraint segment.
  The worst triangle is one: a quality node 0.8 mm from a 245 m CORINE edge.
  20b R8 foresaw this ("if a quality node is ever found beside a segment,
  the cheap fix is a skip reason") and left it unmeasured.
- **The final check is second**: source nodes inserted beside a constraint,
  with no foot rule on that path.
- **20b's rule works where it applies**: 38 refinement apexes, from about
  46 000 refinement insertions.
- **Input vertices are the rest**: an input vertex close to another input
  segment or nearly in line with its neighbours (a cap whose long side is a
  constraint: 133; whose long side is free: 383), and the 185 slivers with a
  side shorter than 10 cm (the CORINE seam's millimetre segments,
  population 1).
- **Caps dominate**: 1 191 of the 1 598 have a largest angle over 120°. By
  Babuška-Aziz these are the ones that spoil slope.
- **Two are forced by an input angle**: 2 of 1 598 have their smallest
  angle between two constraint edges meeting at a vertex (the one kind of
  sliver no Steiner point can remove), and both are gaps in the CORINE
  source closed by a short third input edge (review round 3, B3; rerun on
  master with the kept prototype, all switches off, `forcedlist.py`). One
  is M5's 1 cm slit: 0.007080° at EPSG:3006 413 677.295 6 331 664.081, two
  79.13 m sides, a 1 cm closing side; in 20c-3's Lagan mesh it is the
  second-worst triangle, and the worst (0.006697°) stands on the same 1 cm
  side from outside the slit, made by the 2 m simplification (M5). The other is 0.674° at 382 340.000 6 253 375.946, two 5.95 m
  sides and a 7 cm closing side. The other 1 596 are avoidable in
  principle; the input-vertex ones only at a cost (M4).

### M2. Three populations, by where they sit

The investigation named two populations. Its review (main session,
2026-10-06) found a third, and M1's tags place all three. Bands are 15 m
either side of a straight CORINE line in EPSG:3035 (18 m for the easting
pair); "rest" is everything outside the four bands. Script `lines.py`,
scratchpad.

| where | triangles | slivers (rate) | a side under 10 cm | apex: quality / input / final check / refinement |
|---|---|---|---|---|
| **1. CORINE seam**, N 3 772 265.53 | 459 | 117 (25 %) | 117 | 0 / 117 / 0 / 0 |
| 1. easting pair E 4 537 573.71 and 4 537 579.53 | 869 | 23 (2.6 %) | 21 | 0 / 23 / 0 / 0 |
| **3. straight lines without mm segments**, N 3 811 923.31 | 354 | 38 (10.7 %) | 0 | 7 / 1 / 9 / 21 |
| 3. E 4 585 680.19 | 387 | 47 (12.1 %) | 0 | 6 / 5 / 21 / 15 |
| **2. rest of the basin** | 861 829 | 1 373 (0.159 %) | 47 | 626 / 555 / 189 (+1 line point) / 2 |

- **Population 1** (the seam's millimetre vertex pairs): every apex is an
  input vertex and every sliver has a side under 10 cm. Only input
  coarsening removes it (option b).
- **Population 3** (the corrected write-up's F5, `cd93b2cb`): two straight
  cuts through single coniferous-forest polygons, CORINE class 312 on both
  sides, so the border carries no land-cover change. Their segments are
  long (median 811 m and 576 m; on the mesh, 261 constraint edges, median
  223 m, longest 1 336 m; sliver long sides median 448 m), and the apex sits
  a median 2.4 m off the line. Unlike the rest, **refinement nodes are the
  largest apex group (36 of 85)**, even with 20b's feet on. None of those 36
  has a foot within 1 m, so 20b's fallback (a footed node inserted after
  all) did not make them; M3's row "A + R" shows that 20b R1's "only the
  triangle holding the node" made most of them: with the quality start and
  final check footed (M3, row A), 34 refinement apexes remain on the two
  lines, and 6 when the holding triangle's three neighbours are searched too
  (row A+R). The 10 final-check apexes left on E 4 585 680.19 have the same
  cause on that path, whose prototype searched only the holding triangle. The final check's source nodes (30) and the quality start
  (13) are the rest, population 2's mechanism. The write-up measures the
  general law behind it: among triangles with a constraint side, the sliver
  share is 0.05 % at 10-30 m edges and 33 % at 800 m or more.
- **Population 2** is the spread-out bulk: quality start and input vertices
  first, the final check third; refinement made 2 of 1 373. So 36 of
  refinement's 38 sliver apexes in the whole basin sit on the two
  population-3 lines.

### M3. Prototypes on Lagan and Numedalslagen

A scratch copy of `include/` at `ed12512`, with each rule behind an
environment switch, built as `_core` (Release, Apple clang) and loaded in
place of the installed module by the driver, which asserts it did. With
every switch off it reproduces master exactly (863 897 triangles, 1 598
slivers). The scratch copy is not kept. Rules prototyped:

- **Q-skip**: the quality start skips a node closer than half a cell to a
  constrained edge of the triangle it lies in.
- **Q-foot**: instead of skipping, it inserts the node's foot (orthogonal
  projection) on that edge, when the foot is at least half a cell from both
  ends and the split children are counter-clockwise; otherwise it skips.
- **P**: the final check (`refine_points.hpp`) inserts the foot of a source
  or DEM point closer than half a cell to a constrained edge of its
  triangle, once per point; the point itself stays a candidate, as 20b's
  fallback. The prototype gave the foot the point's own z.
- **R**: 20b's refinement foot also searches the three neighbours of the
  triangle holding the node (not across a constrained edge), and splits the
  neighbour's edge.
- **A** = Q-foot + P. **pen**: the quality start accepts a candidate only if
  it does not lower the worst angle among the triangles it replaces (its
  cavity, computed read-only before inserting); "sum" variants used the
  summed angle deficit instead and were worse.
- **post**: after refinement converges, the quality pass once more over the
  whole mesh, then refinement again (increment 20's C3).
- Input variants (M4): `1 cm` and `1 m` grid (shapely `set_precision`, the
  investigation's E2 method), `cov 2 m` (shapely `coverage_simplify`,
  tolerance 2 m, on the 1 cm input).

"Pop. 1", "pop. 3" and "rest" are M2's bands. Every run reports every node
of the original DEM within 10 m (the CLI's own check; `@tester`'s oracle is
the gate, not this). Battery power; no timing is compared.

**Lagan** (GLO-30 via a 31 m grid, `--tolerance 10`):

| run | triangles | slivers | pop. 1 | pop. 3 | rest | caps | < 0.05° | < 10° | worst | median | quality nodes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| master | 863 897 | 1 598 | 140 | 85 | 1 373 | 1 191 | 216 | 5.36 % | 0.000412° | 36.62° | 195 158 |
| Q-skip | 851 797 | 1 003 | 140 | 63 | 800 | 594 | 187 | 4.95 % | 0.002384° | 36.48° | 188 903 |
| Q-foot | 867 652 | 917 | 140 | 78 | 699 | 521 | 186 | 4.71 % | 0.002384° | 36.87° | 197 742 |
| P | 864 205 | 1 390 | 140 | 62 | 1 188 | 983 | 200 | 5.26 % | 0.000412° | 36.64° | 195 158 |
| A | 868 036 | 738 | 140 | 50 | 548 | 340 | 172 | 4.62 % | 0.002384° | 36.87° | 197 742 |
| **A + R** | **867 764** | **709** | 140 | 21 | 548 | 311 | 173 | 4.61 % | 0.002384° | 36.87° | 197 742 |
| **A + R + pen** | **758 308** | **720** | 140 | 27 | 553 | 318 | 172 | 5.03 % | 0.002384° | 34.90° | 130 152 |
| A + R + pen + post | 818 088 | 715 | 141 | 23 | 551 | 314 | 172 | 4.53 % | 0.002384° | 36.67° | 165 278 |
| A + R, 1 cm input | 867 786 | 557 | 7 | 21 | 529 | 303 | 9 | 4.57 % | 0.008589° | 36.87° | 197 623 |
| A + R, 1 m input | 867 693 | 448 | 0 | 13 | 435 | 273 | 2 | 4.34 % | 0.018157° | 36.87° | 197 739 |
| A + R, cov 2 m input | 805 776 | 305 | 0 | 12 | 293 | 172 | 9 | 3.01 % | 0.006608° | 36.87° | 177 620 |
| A + R + pen, cov 2 m input | 704 966 | 312 | 0 | 19 | 293 | 183 | 9 | 3.17 % | 0.006608° | 35.29° | 114 506 |
| *write-up E11: segments ≤ 100 m, 1 cm* | *1 073 613* | *461* | | | *438* | | *95* | | *0.00048°* | | |
| *write-up E12: the cuts merged* | *863 483* | *1 515* | | | *1 373* | | *210* | | *0.00041°* | | |

**Numedalslagen** (DTM10, projected path, CORINE, `--tolerance 10`, the
command in `rasputin_scratch/norway/numedalslagen/..._stats.md`):

| run | triangles | slivers | caps | < 0.05° | < 10° | worst | median | quality nodes |
|---|---|---|---|---|---|---|---|---|
| master | 1 287 334 | 1 644 (0.128 %) | 1 323 | 188 | 2.60 % | 0.000399° | 36.87° | 316 491 |
| **A + R** | **1 290 434** | **481 (0.037 %)** | 174 | 117 | 2.20 % | 0.000680° | 36.87° | 318 861 |
| **A + R + pen** | **1 073 884** | **447 (0.042 %)** | — | 117 | 1.96 % | 0.000398° | 35.54° | 191 968 |

Numedalslagen's master slivers had the quality start as apex in 1 129 of
1 644, input vertices in 301, refinement in 92; after A + R: input 285,
quality 46, refinement 32.

What M3 shows:

- **Constraint-aware insertion on all three paths (A + R) more than halves
  the slivers for under 0.5 % more triangles**: Lagan 1 598 to 709, caps
  1 191 to 311; Numedalslagen 1 644 to 481, caps 1 323 to 174. The
  write-up's E11 leaves fewer (461) but costs 24 % more triangles.
- **A foot beats a skip** in the quality start (917 against 1 003 slivers),
  and the final check's foot (P) matters only together with it.
- **Searching the neighbours (R) is what removes population 3's refinement
  slivers** (M2).
- **The penalty (pen) is a triangle saving, not a sliver fix**: 12.6 % fewer
  triangles on Lagan and 16.8 % on Numedalslagen at the same sliver count,
  for a lower median angle (36.9° to 34.9° and 35.5°) and slightly more
  triangles under 10° on Lagan (4.61 to 5.03 %; Numedalslagen falls, 2.20
  to 1.96 %). The hard 25° pass inserts about 68 000 (Lagan) and 127 000
  (Numedalslagen) nodes that do not raise the worst angle where they go in.
- **Quality during refinement (post) buys little**: 60 000 more triangles
  (+7.9 %) for 5 fewer slivers and 4.53 % instead of 5.03 % under 10°.
- **The worst angle moves, but stays set by the input** (M5 has the
  triangles). A + R raises it from master's 0.000412° to 0.0024° on Lagan
  (master's worst is a cap that the quality start's foot removes) and from
  0.000399° to 0.000680° on Numedalslagen. After A + R the worst triangle on
  both stands on a millimetre CORINE segment: a needle on 3 mm on Lagan, a
  cap on 7 mm on Numedalslagen (M5). A needle's angle is about the
  segment's length over the distance to the opposite corner (in radians). The quality start stops
  on a needle once its circumradius is under a cell diagonal (20 R5), which
  leaves the opposite corner up to about two diagonals away: Lagan's worst,
  3.2 mm over 75 m, has a circumradius of 38 m against the 43.8 m floor.
  Even a node one cell away would give only 0.006° on that segment, so with
  the input as given no on-lattice rule gets near 1°; input coarsening moves
  it further (0.0066° to 0.018°, M4), until a gap between two polygons
  and the simplification beside it set it (M5).
- **What is left after A + R is input geometry**: of Lagan's 709, 185 have
  a side under 10 cm (population 1 and its like), and 494 of the other 524
  have an input vertex as apex (median 0.76 m off the long side, shortest
  side median 1.45 m: CORINE outlines with vertices about a metre apart that
  zigzag). A lattice-restricted pass cannot fix features
  smaller than a cell (20 R5's floor); only off-lattice Steiner points
  (Ruppert proper, graded to a metre) or a coarser input can.

### M4. Input coarsening (option b)

- **The noder's own lever does not work today.** `rasputin mesh ...
  --snap-spacing 0.1` (and `1`) on the Lagan case stops with
  `MalformedOutput. the noder produced output violating its own guarantees;
  this is our bug chain 319 does not tile the buffer`
  (`include/terrain/noding/noded_pslg_builder.hpp@ed125121:304`). A defect
  in the noder at coarse spacings, **outside 20c**; not investigated. It has
  its own unnumbered ROADMAP row (approved by Ola 2026-10-06), with the
  command that reproduces it.
- **In Python, before the noder, it works**: shapely's `set_precision`
  (snap rounding through GEOS) keeps all 7 319 CORINE features valid at 1 cm,
  10 cm and 1 m, each on its own (on Numedalslagen at 1 cm, neighbours
  stop sharing their borders, so 20c-3 no longer uses it: M6 (a)). At 1 cm it removes population 1 (140 to 7) and the worst
  angle rises from 0.0024° to 0.0086° (M3's 1 cm file; the file the kept
  script makes gives 0.0066°, a gap between two polygons, M5); at 1 m to
  0.018°, with 448 slivers.
- **Coverage simplification** (`coverage_simplify`, GEOS 3.12's
  coverage-preserving Visvalingam-Whyatt, shared borders simplified once so
  the partition stays a partition) at 2 m removes 80 132 of 1 744 066
  vertices and takes the slivers to 305, with 7 % fewer triangles. It is the
  only run that touches the input-vertex apexes (494 to 266).
- **E12's merge** (same class on both sides) is a semantic rule, not a
  tolerance: a CORINE border between two polygons of one class marks no
  change of land cover. It removes population 3 by removing its lines; A + R
  already takes population 3 from 85 to 21.

### M5. Re-measured for design review round 1

The scratch copy of M3 was not kept, so the prototype was rebuilt from
`ed12512` (Release, Apple clang, loaded in place of the installed module by
the same driver). With every switch off it reproduces master exactly
(863 897 triangles, 1 598 slivers, 0.000412°). Its rebuild of M3's A + R
gives 867 718 triangles and 695 slivers on Lagan, against M3's 867 764 and
709; the first copy is gone, so the difference is not traced, and **every
gate below is set from this section's figures**, not M3's. New in the
rebuild: the neighbour search on all three paths (R1 as designed; M3 had
it in refinement only), the final-check foot's z as R4.3 says or as the
point's own (to compare), the soft criterion as R7, and **B**: when the
quality start's walk to its node is stopped by a constraint edge, the
node's foot on that edge goes in instead (R8). Counts and angles only; most
runs on AC power, none timed. In every run the CLI's own check holds (every
node of the original DEM within 10 m of the mesh).

**The prototype is kept** (review round 2, S4), outside git:
`../rasputin_scratch/20c-prototype/proto-20c.patch` (it applies to
`ed12512` with `git apply`; each rule behind an environment switch named in
`run2.sh` beside it), with the driver (`e7.py`), the analysis scripts and
the input scripts (`snapfeat.py`, `covsimp.py`). For design review round 2
the patch was rebuilt from the round-1 edits and rerun: it reproduces this
section's 20c-1 and 20c-2 rows on both catchments and the 20c-3 row on
Lagan exactly (triangles, slivers, worst angle).

**Lagan:**

| run | triangles | slivers (share) | of them, a side under 10 cm | caps | < 0.05° | < 10° | worst | median |
|---|---|---|---|---|---|---|---|---|
| master | 863 897 | 1 598 (0.185 %) | 185 | 1 191 | 216 | 5.36 % | 0.000412° | 36.62° |
| A + R, rebuilt (as M3: neighbour search in refinement only, the final check's foot at the point's own z) | 867 718 | 695 (0.080 %) | 185 | 297 | 172 | 4.60 % | 0.002384° | 36.87° |
| 20c-1 as designed, foot at the point's own z | 868 254 | 682 (0.079 %) | 185 | 286 | 172 | 4.53 % | 0.002384° | 36.87° |
| **20c-1 as designed** (R4.3's z) | **867 612** | **679 (0.078 %)** | 185 | 285 | 172 | 4.53 % | 0.002384° | 36.87° |
| 20c-1 + soft criterion | 757 981 | 694 (0.092 %) | 185 | 294 | 172 | 4.97 % | 0.002384° | 34.94° |
| 20c-1 + soft, candidates beside a sliver exempt | 757 955 | 694 (0.092 %) | 185 | 294 | 172 | 4.97 % | 0.002384° | 34.94° |
| 20c-1 + B, hard 25° rule | 994 639 | 436 (0.044 %) | 185 | 68 | 168 | 3.43 % | 0.002384° | 36.87° |
| **20c-2 as designed** (soft + B) | **799 120** | **437 (0.055 %)** | 185 | 69 | 168 | 4.40 % | 0.002384° | 35.61° |
| 20c-1 + soft, 1 cm input | 757 755 | 566 (0.075 %) | 59 | 287 | 43 | 4.95 % | 0.006608° | 34.94° |
| 20c-1 + soft, 2 m coverage input | 704 353 | 292 (0.041 %) | 19 | 162 | 9 | 3.11 % | 0.006608° | 35.33° |
| 20c-2 + 1 cm snap + 2 m coverage input (20c-3 as designed to round 2) | 742 689 | 151 (0.020 %) | 19 | 44 | 5 | 2.61 % | 0.006608° | 35.90° |
| **20c-2 + 2 m coverage input** (20c-3, no 1 cm snap, round 3; simplified in the source's EPSG:3035, then moved to EPSG:3006, as `run5.sh` does; moved first, as 20c-3 is designed, round 4: 742 725 triangles, every other column the same) | **742 717** | **151 (0.020 %)** | 19 | 44 | 5 | 2.61 % | 0.006697° | 35.90° |

**Numedalslagen:**

| run | triangles | slivers (share) | of them, a side under 10 cm | caps | < 0.05° | < 10° | worst | median |
|---|---|---|---|---|---|---|---|---|
| master (M3) | 1 287 334 | 1 644 (0.128 %) | 118 | 1 323 | 188 | 2.60 % | 0.000399° | 36.87° |
| **20c-1 as designed** | **1 290 807** | **433 (0.034 %)** | 118 | 126 | 117 | 2.16 % | 0.000680° | 36.87° |
| 20c-1 + soft criterion | 1 073 767 | 440 (0.041 %) | 118 | 126 | 117 | 1.94 % | 0.000398° | 35.54° |
| 20c-1 + soft, candidates beside a sliver exempt | 1 073 813 | 434 (0.040 %) | 118 | 124 | 117 | 1.94 % | 0.000398° | 35.54° |
| 20c-1 + B, hard 25° rule | 1 945 417 | 297 (0.015 %) | 118 | 37 | 114 | 0.73 % | 0.008623° | 37.87° |
| **20c-2 as designed** (soft + B) | **1 140 774** | **298 (0.026 %)** | 118 | 32 | 114 | 1.46 % | 0.008623° | 36.09° |
| 20c-2 + 1 cm snap + 2 m coverage input (round 3) | 1 014 831 | 101 (0.010 %) | 10 | 87 | 72 | 0.83 % | 0.000006° | 36.03° |
| **20c-2 + 2 m coverage input** (20c-3, no 1 cm snap, round 3) | **1 013 251** | **32 (0.003 %)** | 4 | 21 | 4 | 0.81 % | 0.012446° | 36.03° |

What M5 shows:

- **R4.3, the foot's height along the line, measured (review B4).** On
  Lagan (the reprojected path) all 793 final-check feet took their z from
  the strip, between the two line points that bracket them. That z differs
  from the point's own by 4.2 m on average and 17.3 m at most: the point is
  up to half a cell (15.5 m) off the line, on a slope. With R4.3's z the
  strip has to add 111 line points instead of 166 (master: 266), the final
  check adds 26 708 points instead of 27 029, the mesh has 642 fewer
  triangles and 3 fewer slivers; 19 footed points later went in as
  themselves (the fallback) instead of 10. The line points' largest error
  is 9.999 m either way. **R4.3 stays.** On Numedalslagen (the projected
  path, z from the DEM at the foot) the final check made 5 feet.
- **The 20c-1 gate's margin now covers only the step from prototype to
  production**, not an unmeasured rule: the final-check foot is in the
  measured figure (679, 0.078 %).
- **L12 (below, R4) never fired** in these runs; R4's foot did 793 times on
  Lagan.
- **The quality start's feet** (5 930 on Lagan, 4 233 on Numedalslagen):
  none was a near line with no usable foot (so `skipped_near_line` was 0),
  and after every one of them the bad triangle that asked for it was gone.
- **The soft criterion refused** 32 499 (Lagan) and 63 383
  (Numedalslagen) candidates. Exempting candidates whose cavity holds a
  sliver changes 0 and 6 slivers: not designed in.
- **B borrows Ruppert's response, not his trigger.** Ruppert splits a
  segment when a point would lie in its diametral circle (the circle on the
  segment as diameter); B splits the constraint line that stops the walk
  from the bad triangle to its node, at the node's foot. Under the hard
  25° rule it costs +15 % and +51 % triangles, so it belongs with the soft
  criterion, not in 20c-1. With it, 20c-2 has a third fewer slivers than
  20c-1 (437 and 298), a quarter of its caps (69 and 32), and 7.5 % and 11.4 % fewer
  triangles than master.

**Which triangle sets the worst angle, and why** (review round 1 B2;
round 2 B2 and S5). A *cap* has one angle over 120° (M1's definition); a
*needle* has one tiny angle and two near 90°. In either, the tiny angle is
about the short side times the sine of the angle beside it, over the long
side, in radians; for a needle that sine is about 1.
The quality start's floor (20 R5) is a circumradius of one cell diagonal:
43.8 m on Lagan's 31 m grid, 14.1 m on Numedalslagen's 10 m grid.

| after | Lagan | Numedalslagen |
|---|---|---|
| master | 0.000412°: a cap, its apex a quality node 0.8 mm from a 245 m constraint edge | 0.000399°: a cap, its apex a quality node, long side 152 m |
| 20c-1 | 0.002384°: a needle on a 3.2 mm CORINE segment, apex a quality node 75 m away. Circumradius 38 m, under the floor, so never a candidate | 0.000680°: a cap (angles 0.00068°, 38.3° and 141.7°) on a 7.1 mm CORINE segment, apex a quality node 369 m away (circumradius 298 m). Its circumcentre lies across two constraint edges, so the walk is blocked and the candidate skipped (20 C1 (a)) |
| 20c-1 + soft criterion alone | the same triangle | 0.000398°: the same 7.1 mm segment, apex now a refinement node 631 m away. Still blocked; the soft criterion leaves the mesh around it coarser, and the triangle keeps whatever corner the rest of the mesh gives it (631 / 369 = 1.71 = 0.000680 / 0.000398) |
| 20c-2 (soft + B) | the same triangle, corner for corner (rerun, round 2) | 0.008623°: a needle on a 4.1 mm segment, apex a quality node 27 m away. Circumradius 13.7 m, under the floor |
| 20c-3 (1 cm snap and 2 m coverage input, on 20c-2) | 0.006608°: a gap in the land cover, made only of input lines (below) | 0.000006°: a flat triangle on two near-coincident CORINE borders that the 1 cm snap pulled apart (M6) |
| 20c-3 (2 m coverage input, no snap, round 3) | 0.006697°: the gap's own 1 cm closing edge seen from outside the gap, apex 78.9 m away, made by the simplification (below) | 0.012446°: three input vertices, a 3.6 mm input edge, apex 15.7 m away |

So the remaining thin triangle after 20c-2 is, on both catchments, a
millimetre CORINE segment with a lattice node about a cell or more away,
which the lattice cannot bring closer. After 20c-3 on Lagan it is the
1 cm edge that closes a gap between two CORINE polygons, with its nearest
vertex pushed 78.9 m away by the simplification.

**The gap: in the input as given** (review round 2 B2, checked on the
inputs directly). A wedge of three input vertices (EPSG:3006
413 677.295 6 331 664.081, 413 602.500 6 331 638.254 and
413 602.501 6 331 638.244) whose three sides are all constraint segments:
the border of a forest polygon (CORINE 312) and the border of an arable
polygon (CORINE 231), 79.13 m long each, which share their far vertex and
end 1 cm apart, closed by a 1 cm edge of a third polygon (CORINE 243).
Between the two borders lies a wedge 1 cm wide and 79 m long that belongs
to no polygon (0.40 m² of gap; the source coverage is valid, since a valid
coverage may have gaps). It sits on M2's easting pair (E 4 537 579.53 in
EPSG:3035), population 1. Since the angle between two constraint segments
meeting at a vertex is a ceiling on the worst angle of any mesh of that
input, every mesh of these inputs has it, at 0.007080° in the source
(20c-1's mesh, rerun, hidden there below the worse needles, 0.0024°) and
after the 2 m simplification, which keeps all three corners; the 1 cm snap
of round 2 moved its two near corners by under a centimetre and lowered it
to 0.006608°. Neither the snap nor `coverage_simplify` closes a gap
between polygons; only a gap-closing step would (GEOS's coverage cleaning,
GEOS 3.14; the venv's shapely has GEOS 3.13). Ruling 5 first said no
gap-closing step; Ola reversed it (ruling 9), and 20c-3's repair closes
this gap (M7). It is one of M1's two slivers
forced by an input angle (M1, last point).

**20c-3's worst triangle on Lagan is made by the simplification**
(review round 4, B1; rerun on the kept prototype, round 4). It is the
triangle on the other side of the same 1 cm closing edge: corners
413 602.500 6 331 638.254 and 413 602.501 6 331 638.244 (the edge) and
413 527.395 6 331 662.317 (the apex, 78.9 m away), 0.006697°; one side,
the 1 cm edge, and the 78.9 m side from 413 602.500 to the apex are
constraint segments, the third is free, so it is not forced. In the
source, the border shared by the 243 and 312 polygons has a vertex at
413 596.806 6 331 640.083, 5.98 m from the closing edge, and that is
where the triangle on that edge has its apex: 0.088°. `coverage_simplify`
at 2 m drops that vertex (it lies 4.5 mm off the straight border from
413 527.395 to 413 602.500), so the nearest vertex beside the edge is
then 78.9 m away, and the angle falls by the same factor, to 0.006697°,
just below the gap's own 0.007080°. Round 3 placed it "150 m west of the
slit", measuring from the gap's far vertex; it stands on the gap's near
end.

**A guard against it, measured, and not designed in.** The guard: before
simplifying, mark the two ring neighbours of each end of every source
edge under 10 cm, and after `coverage_simplify` put back every marked
vertex it dropped, in ring order and at the same coordinates in both
polygons that share it, so the result is still a coverage
(`covkeep.py`, kept beside the prototype). Counting an edge or point that
two polygons share once: Lagan has 417 such edges, and the guard puts back 39
points; Numedalslagen (its CORINE as read for the domain's box, 8 017
polygons, not filtered to the domain), 1 362 and 71. Meshes on 20c-2:

| input | Lagan: triangles | slivers | side < 10 cm | worst | Numedalslagen: triangles | slivers | side < 10 cm | worst |
|---|---|---|---|---|---|---|---|---|
| 2 m coverage, moved to the mesh's CRS first (as designed) | 742 725 | 151 | 19 | 0.006697° | 1 013 251 | 32 | 4 | 0.012446° |
| the same, with the guard | 742 831 | 151 | 19 | 0.007080° | 1 013 557 | 32 | 4 | 0.012446° |

The guard lifts Lagan's worst angle only to the gap's own forced
0.007080°, which no step could pass while ruling 5 stood (ruling 9 now
repairs the gap, M7), removes no sliver on
either catchment, and leaves Numedalslagen's worst triangle where it is
(its 3.6 mm edge is not an edge of any land-cover ring, so the guard does
not reach it; not traced further). For about 15 lines it buys 6 % of one
angle on one catchment, so **20c-3 has no guard**, and its gates allow for
the simplified triangle (Lagan's worst-angle gate, 0.005°, sits below
0.006697°). The mechanism to keep in mind: a simplification that drops a
vertex beside a short input edge moves that edge's nearest apex away and
lowers its angle by the ratio of the distances (here 78.9 / 5.98 = 13.2,
and 0.088° / 0.006697° = 13.2).

**M3's 1 cm run (0.008589°) used a different file.** The CLI read 330 529
vertices from it, against 330 559 from the 1 cm file that the kept script
(`snapfeat.py`) makes today, which M5's 1 cm run and the 2 m coverage input
were made from. In M3's file no triangle under 0.1° had its smallest angle
between two constraint edges, so the wedge was not in it. That file is
gone and how it differed is not traced; every 1 cm and 2 m figure in M5
is from the file that can be remade.

### M6. Ola's outline rule, measured on 20c-3's input (design review round 3)

Ola, 2026-10-06, proposed for 20c-3: "if a corine polygon intersects
marginally with the catchment polygon, we should simply adjust it so that
it's excluded. At the same time, the corine boarder (the neighbouring
polygons to the marginally, and now excluded polygon), must be adjusted to
the catchment polygon border to compensate?" Round 2 measured only the
first half (thin pieces handed to a neighbour), and on 20c-2's mesh. Round
3 measures both halves on the mesh 20c-3 makes, as review round 3 (B1)
asked.

All runs: the kept prototype (M5) with 20c-2's switches, rebuilt from
`ed12512` for this round; inputs and scripts in
`../rasputin_scratch/20c-prototype/` (`run5.sh` lists the steps). Counts and
angles only, none timed. The Norwegian CORINE is the GeoPackage rewritten
as GeoJSON (`outline.py` at W = 0); with it 20c-2's Numedalslagen mesh is
reproduced exactly (1 140 774 triangles, 298 slivers, 0.008623°).

**(a) 20c-3's input, corrected: no 1 cm snap.** Round 2's 20c-3 snapped
each CORINE polygon to a 1 cm grid on its own (shapely `set_precision`)
before `coverage_simplify`. On Numedalslagen that breaks the coverage:
two neighbouring forest polygons, snapped apart, no longer share 9.7 km
of border inside the catchment (`covcheck.py`; the source and the 2 m
simplification of the source alone have none), and the mesh gets a flat
triangle of 0.000006°, two constraint edges from one input vertex lying
almost on each other where those polygons meet, with 101 slivers.
Simplifying the source directly gives 32 slivers and 0.0124°, and on
Lagan the same 151 slivers as with the snap (0.006697° instead of
0.006608°, M5). `set_precision` works on one polygon at a time, so
nothing makes two neighbours round their shared border the same way.
**So 20c-3 drops the
1 cm snap** (design below), and 20c-3's baseline is "2 m coverage, no
snap": 151 slivers on Lagan, 32 on Numedalslagen.

**(b) Where 20c-3's slivers sit.** On Lagan, 52 of the 151 have their
centre within 5 m of the outline, 34 of them with a CORINE line within
5 m (24 whose long side lies inside, along a border that runs beside the
outline, 21 of these with two corners on it; 10 with the long side on the
outline; `strip2.py`); 95 lie within 50 m.
On Numedalslagen 24 of 32 lie within 5 m. Round 2's thin-piece form reaches
6 (W = 10 m) and 11 (W = 20 m) of the 151 on round 2's 20c-3 mesh, as
review round 3 found (`strip.py` with `outline.py`'s pieces, rerun): 4 and
7 %.

**(c) The rule, second form** (`snapline.py`, then `tolines.py`). Every
part of a land-cover border within D of the outline is put onto the
outline, and the outline carries it from then on:

1. each border edge that comes within D of the outline is cut into pieces
   of at most D, the same cut for both polygons that share it;
2. each point within D moves to its nearest point on the outline, then to
   an outline vertex within D/2 of that, or else to the nearest multiple
   of D along the outline, so two of these points never lie centimetres
   apart;
3. consecutive moved points are joined along the outline through its own
   vertices; a cut point that did not move is kept only next to one that
   did (the step off the outline), and one in line with its neighbours is
   dropped;
4. the stretches that now lie on the outline are dropped from the
   linework, so each remaining border ends where it meets the outline.

A thin piece narrower than D vanishes, its neighbour reaching the
outline in its place (the first half of Ola's rule); a border running
beside the outline within D goes onto it (the second half). Only vertices
move, and a shared border moves the same way in both polygons: after the
rule no shared border inside the catchment is mismatched (`covfar.py`,
both catchments, D = 5 and 10 m). No overlay.

Step 4 is needed, and shapely's `snap` does not do step 2.
`shapely.snap(polygon, outline, 5)` moves a vertex only onto an outline
*vertex* within 5 m; a border 3 m from a straight outline segment stays
where it is (checked on a toy coverage). Without step 4 (borders written
as polygons lying along the outline) Numedalslagen got 43 slivers, against
10 with it (both with the first version's defect below), and a worst angle
of 0.0027°: the mesh had new input vertices on the outline, millimetres to
centimetres apart, where the rule had put no point (checked at one,
EPSG:25833 94 792.151 6 709 280.591, 50 m from any outline vertex). On
Lagan, whose outline is axis-parallel, step 4 made no difference (61
slivers either way). The first version also lost an outline vertex where a
moved point was rounded onto it, which broke 104 m of shared border on
Numedalslagen; fixed before the runs below.

**(d) Measured.** D = 5 m and 10 m; "area" is the land cover that changes
class inside the catchment (`area.py`).

| input (all on 20c-2) | Lagan: triangles | slivers | side < 10 cm | worst | Numedalslagen: triangles | slivers | side < 10 cm | worst |
|---|---|---|---|---|---|---|---|---|
| 2 m coverage (20c-3 without the rule; Lagan's coverage simplified before the move to EPSG:3006 in all three rows, M5) | 742 717 | 151 (0.020 %) | 19 | 0.006697° | 1 013 251 | 32 (0.003 %) | 4 | 0.012446° |
| **2 m coverage, rule at D = 5 m** | **743 385 (+0.09 %)** | **61 (0.008 %)** | 2 | 0.006697° | **1 015 416 (+0.21 %)** | **4 (0.0004 %)** | 0 | **0.832°** |
| 2 m coverage, rule at D = 10 m | 742 465 (−0.03 %) | 61 (0.008 %) | 2 | 0.006697° | 1 013 906 (+0.06 %) | 4 (0.0004 %) | 0 | 0.832° |
| the input as given (20c-2) | 799 120 | 437 | 185 | 0.002384° | 1 140 774 | 298 | 118 | 0.008623° |
| as given, rule at D = 5 m | 799 897 (+0.10 %) | 345 | 166 | 0.002384° | 1 142 932 (+0.19 %) | 271 | 114 | 0.008623° |

Area changing class at D = 5 m: 0.71 ha on Lagan (6 441 km²) and 0.77 ha on
Numedalslagen (5 548 km²), about 0.0001 % of each; at 10 m, 2.65 and
2.36 ha. CORINE's own minimum mapping unit is 25 ha.

- **On 20c-3's input it removes 60 % and 88 % of the slivers** (151 to 61,
  32 to 4) for under 0.25 % more triangles, and every sliver within 20 m
  of the outline (Lagan: 83 to 0; Numedalslagen: 27 to 0; `strip.py`).
  On Lagan 38 more go, 5 to 50 m inside.
- **Numedalslagen's worst angle rises from 0.0124° to 0.83°**, a triangle
  of two refinement nodes and an input vertex, its shortest side 20 m. Lagan's
  stays at 0.006697°, the triangle the 2 m simplification makes on the
  slit's 1 cm closing edge (M5), which the rule does not change.
- **On the input as given it helps less** (21 % and 9 %), because there
  the millimetre segments elsewhere dominate; the worst angles are those
  segments' and do not move.
- **D = 5 m and 10 m give the same slivers**; 5 m moves about a quarter
  to a third of the area 10 m moves (Lagan 0.71 / 2.65 = 27 %,
  Numedalslagen 0.77 / 2.36 = 33 %). D is in metres in the computation CRS; it assumes DEM cells of
  10 to 31 m (the two catchments, the largest inputs it was measured at)
  and is about half a cell or less there.
- An early version run in Lagan's features' own CRS (EPSG:3035, the
  outline moved there vertex by vertex) gave 301 slivers, 227 with a side
  under 10 cm; not traced further. Every run in the table works in the
  computation CRS, as production will (design below).

**Verdict: it pays off on 20c-3's input, so 20c-3 builds it (its third
step, `--features-outline-snap`), on by default at 5 m (Ola, question 6,
ruling 6).** Its tests, gate and lines below are the same under either
default.

### M7. Gaps in the land cover, and the repair's tolerance (`@architect`, 2026-10-08, for ruling 9)

One short probe, to set the repair's tolerance (S, below). Shapely 2.2.0
with GEOS 3.14.1 in a scratch venv (the repository's venv has 2.1.2 with
GEOS 3.13, which lacks `coverage_clean`); AC power; one run each, no
warm-up. The CORINE polygons that reach the domain, moved vertex by vertex
into the computation CRS (Lagan: the GeoJSON in EPSG:3035 to EPSG:3006;
Numedalslagen: the GeoPackage, already in EPSG:25833). The scripts are kept
in `../rasputin_scratch/20c-prototype/m7/` (`gaps.py`, `snap.py`,
`clip.py`; `clip.py` as kept is the hull version, which `sed` made from the
first). A gap is a hole in the union of those polygons; its width is twice
its largest inscribed circle's radius. Mesh figures were not measured: no
C++ build was allowed in this round, and the meshes are `@perf`'s
measurement of the built code.

| | Lagan | Numedalslagen |
|---|---|---|
| polygons reaching the domain; vertices | 2 790; 1 006 098 | 1 611; 1 056 147 |
| gaps touching the domain | 1, the slit (M5): width 1 cm, 0.40 m² | 0 |
| gaps of any width between 0 and 100 m, anywhere | that one | none |
| finding the gaps (union and widths) | 0.83 s | 0.79 s |
| `coverage_is_valid` (the old design's check) | 15.76 s | 15.84 s |
| ring edges under 10 cm, counted per polygon | 356 | 464 |
| `coverage_clean`, gap width 1 m, snapping GEOS's default | 3.11 s; 0.40 m² moved; edges under 10 cm 347; the slit's two corners still 1.005 cm apart | 4.99 s; 461.88 m² moved; 38 |
| the same, snapping 0 | 2.92 s; 0.40 m²; 357; corners apart | 4.40 s; 0.00 m²; 464 |
| the same, snapping 1 cm | 3.11 s; 0.64 m²; 41; corners apart | 4.95 s; 494.70 m²; 51 |
| the same, snapping 5 cm | 3.10 s; 0.72 m²; 5; **the corners one vertex** | 4.90 s; 1 085.05 m²; 1 |
| clipped to the domain's convex hull grown by 100 m (30d's read region) first: vertices; clip; `coverage_clean` at 5 cm and 5 cm | 328 372; 0.36 s; 0.75 s | 407 753; 0.34 s; 0.76 s |
| the same, clipped to the domain grown by 100 m (mitred) | 310 748; 3.41 s; 0.65 s | 290 567; 1.01 s; 0.54 s |

Every cleaned coverage passed `coverage_is_valid`, and at gap widths of
0.1, 1 and 5 m the results were the same (no gap between 5 cm and 100 m to
fill). The other 5 690 and 2 877 holes in each union are 100 m wide or more
and lie outside the domain: the land cover of polygons that do not reach it.

- **The slit is the only gap on either catchment**, and 5 cm is the
  smallest snapping distance probed that joins its two open corners.
- **GEOS's default snapping distance is not zero**: on Numedalslagen, which
  has no gap, it moved 462 m². 20c-3 passes S explicitly.
- **The check costs more than the repair.** `coverage_is_valid` took 15.8 s,
  `coverage_clean` 3 to 5 s, and after the clip to the read region about
  0.75 s.
- **The hull is the cheaper region to clip to**: a convex region with few
  vertices, 0.35 s, against 1 to 3.4 s for the mitred domain on these
  outlines (Lagan's is a staircase).

### M8. The repair's tolerance on the built code: 5 cm, 25 cm and 1 m (`@architect`, 2026-10-08, for question 9)

Ola, on question 9: try 1 m as well; he thinks 5 cm is still only noise.
One short probe, not a sweep, on the built code at `9a04c3c1` with the
worktree's `.venv` (shapely 2.2.0, GEOS 3.14.1; the extension built after
the last C++ change, `50544838`), AC power, one run each, no warm-up. The
defaults (merge on, simplification off, outline rule at 5 m) plus
`--features-repair S`, with the catchment arguments of
`docs/benchmarks/quick/cases.toml` on branch `worktree-perf-process`
(tolerance 10 m). A driver runs `rasputin mesh` in-process and wraps
`shapely.coverage_clean` and `snap_to_outline` to keep the repair's input
and output and the domain; the mesh checks are `@perf`'s `gate.py`,
`corners.py` and `farcorner.py`
(`docs/benchmarks/2026-10-08/20c-3/scripts/`). Scripts, `--stats` files
and logs: `docs/benchmarks/2026-10-08/20c-3/m8/`. The 5 cm and 1 m runs
were repeated and gave the same `.vtk` byte for byte; times are the first
runs'.

Unlike M7, the land-cover figures are the repair's own input: the polygons
after the clip to the read region (Lagan 3 674 pieces, 380 826 vertices;
Numedalslagen 3 297, 603 977). Area moved is M7's measure, the symmetric
difference of each polygon before and after, summed, so each square metre
given from one polygon to another counts twice.

| S | Lagan | Numedalslagen |
|---|---|---|
| **5 cm** (the default when probed; 1 m since Ola's answer): slivers under 1°; worst angle; triangles | 180 (0.0225 %); 0.1917°; 799 368 | 158 (0.0138 %); 0.6302°; 1 141 983 |
| 25 cm | 176 (0.0220 %); 0.1917°; 799 336 | 158 (0.0138 %); 0.6302°; 1 141 959 |
| **1 m** | **26 (0.0033 %); 0.7293°; 798 554** (−0.10 %) | **2 (0.0002 %); 0.8322°; 1 118 006** (−2.1 %) |
| ring edges under 10 cm (M7's count), before the repair → after, at 5 cm / 25 cm / 1 m | 231 → 0 / 0 / 0 | 393 → 2 / 0 / 0 |
| land cover moved, in all; inside the domain (m², counted twice), at 5 cm / 25 cm / 1 m | 0.6; 0.5 / 370; 83 / 13 342; 6 789 | 548; 151 / 1 128; 151 / 42 164; 17 020 |
| thickest piece of moved land cover, at 5 cm / 25 cm / 1 m | 0.010 / 0.241 / 0.988 m | 0.011 / 0.250 / 0.992 m |
| gaps touching the domain: before; after, at every S | the slit (1 cm wide, 0.40 m²); none | none; none |
| polygons the repair empties, at 1 m | one, 6.35 m² and 0.63 m wide, 1.5 km outside the outline | one, under 0.01 m² and 6 mm wide, at every S (the clip's crumb) |
| `--stats`' clean-up time (the whole land-cover step), at 5 cm / 1 m; `coverage_clean` alone | 5.86 / 5.98 s; 0.90 s at every S | 5.51 / 5.46 s; 1.22 to 1.25 s |
| the whole run, at 5 cm / 1 m (`--stats` total) | 21.80 / 21.62 s | 11.28 / 10.82 s |
| the slit's two open corners both mesh vertices? (`corners.py`), at every S | no: one is, the other is 1.0 cm from the nearest vertex | (not on this catchment) |
| slivers at the slit's far corner (`farcorner.py`), at every S | 0 of the 6 triangles there | |
| slivers with a side under 10 cm; with the centre within 20 m of the outline, at 5 cm / 1 m | 0 / 0; 0 / 0 | 0 / 0; 1 / 0 |

- **What 1 m removes is borders of two polygons that come within 1 m of
  each other**, not mismatched shared borders: the repair's input is a
  valid coverage on both catchments, with no overlap (`m8ov.py`), and
  M7 found no gap between 5 cm and 100 m wide. At 1 m the moved land cover
  is long thin pieces: those thicker than 25 cm are 198 pieces, 45 km of
  border on Lagan (median 93 m long, the longest 7 km, about 0.3 m
  thick on average), and 4 629 pieces, 123 km on Numedalslagen (median
  1.5 m long). They lie spread over the basins, not on M2's CORINE seam
  or easting pair (2 pieces of 198 within 15 m of them on Lagan). No
  vertex moves farther than S (the thickest piece is 0.99 m).
- **That is where the slivers were.** Of the 5 cm meshes' slivers, 157 of
  180 (Lagan) and 156 of 158 (Numedalslagen) lie within 2 m of land cover
  the 1 m repair moved; on the 1 m meshes, 2 of 26 and 0 of 2
  (`m8moved.py`). The check can fail: it counts slivers of either mesh
  against the same moved land cover, and the 1 m meshes' remaining slivers
  mostly lie elsewhere.
- **25 cm gives nothing over 5 cm on the mesh** (4 slivers fewer on Lagan,
  none on Numedalslagen, the same worst angles). The gain is between
  25 cm and 1 m; no S between them, and none above 1 m, was probed.
- **The cost**: land cover moved by under 1 m along about 45 km (Lagan) and
  123 km (Numedalslagen) of border, 6 789 m² and 17 020 m² inside the
  domains counted twice, so about 3 400 m² and 8 500 m² given from one
  polygon to another: 5 × 10⁻⁷ and 1.5 × 10⁻⁶ of their 6 441 and
  5 548 km². CORINE is mapped at 1:100 000, with a 100 m minimum width and
  a positional accuracy the producers state as better than 100 m; a 1 m
  move is 1 % of the smallest feature it may hold. No measurable time.
- **Scale.** S = 1 m assumes a land-cover source whose borders carry no
  meaning below a metre (CORINE: 1:100 000) and DEM cells of 10 to 31 m;
  checked on these two catchments' CORINE, about 380 000 and 600 000
  vertices after the clip, the largest inputs probed. On a source mapped
  to the decimetre it moves real borders; the flag is there for that.
  MapBiomas polygons made by rasputin from a 30 m raster (São Francisco)
  have no two vertices closer than a pixel side unless the conversion
  puts them there; not measured.

## What Ola gets

Measured with `--stats`' own quality table (share under 1°, worst angle) and
the triangle count, on two catchments, with the height tolerance kept
("every DEM node within the tolerance", checked by `@tester`'s oracle) and
the output independent of `threads`:

| after | Lagan: slivers (share) | triangles | worst angle | Numedalslagen: slivers (share) | triangles | worst angle |
|---|---|---|---|---|---|---|
| master today | 1 598 (0.185 %) | 863 897 | 0.000412° | 1 644 (0.128 %) | 1 287 334 | 0.000399° |
| PR 20c-1, constraint-aware insertion (built, measured) | 679 (0.078 %) | +0.43 % | 0.0024° | 433 (0.034 %) | +0.27 % | 0.00068° |
| PR 20c-2, the soft criterion and the line split | about 440 (0.055 %) | −7.5 % | 0.0024° | about 300 (0.026 %) | −11 % | 0.0086° |
| PR 20c-3, input coarsening at 2 m (flag, off by default) | about 150 (0.020 %) | −14 % | 0.0067° | about 30 (0.003 %) | −21 % | 0.012° |
| PR 20c-3 with Ola's outline rule at 5 m (its default, ruling 6) | about 60 (0.008 %) | −14 % | 0.0067° | about 4 (0.0004 %) | −21 % | 0.83° |
| PR 20c-3 with the repair at 1 m too (its default, ruling 9, question 9 as ruled 2026-10-08) | on a mesh with the defaults, which do not coarsen (M8): 26 (0.0033 %); at 5 cm, 180 (0.0225 %); on the land cover the slit is closed and its 1 cm edge gone (M7, M8) | 798 554; at 5 cm 799 368 | 0.73°; at 5 cm 0.19° | on a mesh with the defaults (M8): 2 (0.0002 %); at 5 cm, 158 (0.0138 %) | 1 118 006; at 5 cm 1 141 983 | 0.83°; at 5 cm 0.63° |

20c-1's row is `@perf`'s measurement of the built code, the same
triangle counts, sliver counts and worst angles as M5's prototype
(`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/README.md@c9f62f7c:172-184`).
Its cost in time, on the same runs (lines 203-221 there): the whole run
−0.26 % on Lagan and +0.18 % on Numedalslagen against master; refine,
which is under 1.2 s of those 81 s and 11 s runs, +4.1 % and +1.7 %
(23 ms and 19 ms). The 1 m benchmark, where no constraint is near enough
for a foot, is 0.3 to 2.0 % faster than master at every thread count
(lines 261-306 there). The other rows are from M5, and M6 for 20c-3.
20c-2 has the line split (Ola's yes to question 4); without
it, the soft criterion alone would leave about 690 and 440 slivers with
−12 % and −17 % triangles, and Numedalslagen's worst angle would fall back
to 0.00040° (M5).

**Against the write-up's two input-side baselines** (main session's
default: 20c's options must beat them). E12 alone (the forest cuts merged
away) leaves 1 515 slivers; 20c-1 leaves 679. E11 (every CORINE segment cut
to 100 m or less, 1 cm grid) leaves 461 at 1 073 613 triangles, 24 % more
than master; 20c-2 leaves 437 at 799 120 (26 % fewer triangles than E11),
and 20c-3 at 2 m 151 at 742 717 (61 with the outline rule). So 20c-1 beats E11 per triangle but not in
count; 20c-2 and 20c-3 beat it in both.

The **worst angle** rises with each PR, and stays set by short input
edges: with the input as given, millimetre CORINE segments; after 20c-3
on Lagan, the 1 cm edge that closes a 1 cm wide gap between two CORINE
polygons in the source data, seen from outside the gap, whose nearest
vertex the 2 m simplification drops (0.088° in the source, 0.006697°
after; M5, "Which triangle sets the worst angle"). The gap itself forces
0.007080° on any mesh of the input as given, so no step short of closing
it gets Lagan above about 0.007°; 20c-3's repair now closes it (ruling 9,
M7), measured on the land cover but not yet on a mesh, and keeping that vertex was measured to
gain only about 6 % (0.006697° to 0.007080°; M5). On Numedalslagen, with the outline rule, 0.83°. So the gates hold each PR to the worst angle of the one
before (and 20c-1 to master's), and promise no fixed angle.

What stays: slivers whose apex is an input vertex about a metre from
another part of the same outline (Lagan: about 490 after 20c-1, 250 after
20c-2). Removing
them needs either Steiner points finer than the DEM cell (option c below,
not recommended) or the input tolerance of 20c-3.

**São Francisco.** The rates are what carries over: about 0.03 to 0.08 %
of triangles after 20c-1 instead of 0.13 to 0.19 %, and 0.03 to 0.06 %
after 20c-2 with 7 to 11 % fewer triangles than today, which matters more
there than anywhere. Its land
cover will be rasputin's own polygons from a raster (MapBiomas), whose
vertices rasputin places; 20c-3's tolerance is then part of that conversion
rather than a change to someone else's data.

## The options, with costs

### (a) The foot rule on every insertion path. Recommended, first, own PR.

20b's rule ("a candidate too close to a constraint segment goes onto the
segment, at its foot") on the two paths that lack it, the quality start and
the final check, plus 20b's own blind spot (it looks only at the triangle
holding the node). M5: slivers 1 598 to 679 on Lagan and 1 644 to 433 on
Numedalslagen, caps (the slope-spoiling kind) 1 191 to 285 and 1 323 to 126,
worst angle at least master's on both, for +0.3 to +0.4 % triangles. About 150 lines of C++ and plumbing. No ruling
needed: it is Ola's C1 lean ("allow adding point to the constraints if that
improves the mesh props") applied where it was missing, and 20b's guarantee
argument carries over unchanged. **Worth shipping first on its own**: it is
the largest gain, the cheapest, and 20c-2's measurements stand on it.

### (b) Input coarsening. Ruled by Ola (question 1); its own PR.

Merging CORINE's millimetre vertex pairs (population 1) and, further,
simplifying outlines within a horizontal tolerance. M4: at 1 cm, population
1 goes (140 to 7 slivers) and the worst angle rises twentyfold; coverage
simplification at 2 m on top takes the total to about 305 with 7 % fewer
triangles. And E12's same-class merge removes population 3's lines.

**Against the input model.** Ola's model (16b R9, "vertices are used as
given") and "tolerances set resolution" pull in opposite directions here: a
horizontal tolerance is a resolution set by a tolerance, but it moves the
input's vertices. So it needed a ruling: Ola ruled the same-class merge on
whenever a class map is used, and the distance tolerance a flag, off by
default (question 1, "Ola's rulings"). Where it lives: Python,
before the noder (`feature_input.py`'s side of the I/O boundary), with
shapely's `coverage_simplify`, which keeps a partition a partition
(`set_precision` polygon by polygon does not: M6 (a), so round 3 dropped
it). Not the noder's `--snap-spacing`: it fails above 1 cm today (M4).
About 60 lines of Python, and about 70 more for Ola's outline rule, on by
default at 5 m (question 6, ruling 6; M6). **Worth its own PR**; it does not depend on (a)
in code, but its gate is measured on 20c-2's mesh, so it goes after 20c-2.

### (c) The soft criterion, and splitting constraints in general

Increment 20's rulings sent three things here: a soft criterion instead of
the hard 25° (C2), points on constraints when they improve the mesh (C1),
and quality during refinement, not only at the start (C3). Measured:

- **The soft criterion (penalty) is worth having, for size.** Accepting a
  quality candidate only when it does not lower the worst angle of the
  triangles it replaces gives 12.6 % (Lagan) and 16.8 % (Numedalslagen)
  fewer triangles at about the same sliver count (M5). That is Ola's
  "24.99 should be ok, if not exploited": a node that leaves its
  neighbourhood no better is not paid for. About 70 lines. Its own PR,
  after (a). Ruled on by Ola (question 2: on, gain 0).
- **Splitting the constraint a quality point lies beyond** (B in M5,
  Ruppert's response to an encroached segment, at the foot, on the blocked
  walk rather than his diametral-circle test). With the soft
  criterion it takes the slivers a third below 20c-1 and Numedalslagen's
  worst angle from 0.00068° to 0.0086°; under the hard rule it would cost
  +15 to +51 % triangles. About 25 lines, in 20c-2 (R8); Ola said yes
  (question 4).
- **Points on constraints**: (a) is that, at the foot. Ruppert's full
  version (split encroached segments at midpoints, recursively, with
  off-lattice circumcentres graded down to the input's local feature size)
  is what it would take to fix the remaining input-vertex slivers. Not
  recommended: the features it would grade to are about a metre (M3), far
  below a 10 to 31 m DEM cell, so it buys triangles the DEM cannot fill with
  information; increment 20's M1 already measured 10 to 20 % more triangles
  for no better result where it is not needed. About 90 lines plus a
  terminator (20 R10).
- **Quality during refinement (C3)**: one pass after refinement converges,
  then refinement again, adds 7.9 % triangles for 5 fewer slivers and
  4.53 % instead of 5.03 % of triangles under 10° (M3, "post"). Not a sliver
  fix; a choice about the angle distribution. Ruled by Ola: not in 20c
  (question 3).

## Design of PR 20c-1: the foot rule on every insertion path

### R1. One geometric helper, three policies

- **New header `include/terrain/mesh/constraint_foot.hpp`**, depending on
  `core`, `predicates`, `lattice_mesh.hpp` and `lawson.hpp` (for
  `LatticeFrame`, which itself depends on those three only); no raster, no
  height.

  ```cpp
  enum class FootStatus : std::uint8_t { None, Hit, NearEnd, NotCounterClockwise };
  struct FootSearch {
      FootStatus status = FootStatus::None;
      std::uint32_t owner = 0;  // the triangle whose edge it is (t or a neighbour)
      unsigned edge = 0;
      MeshVertex at{};          // the foot; meaningful for Hit only
  };
  // The first constrained, non-frozen edge closer than delta (world
  // distance in frame f) to p, p not exactly on it: t's edges in edge
  // order, then, for each of t's unconstrained edges in edge order, the
  // neighbour's other two edges. None if there is none; NearEnd if that
  // edge's foot lies within delta of either end of it; NotCounterClockwise
  // if a child of the split would not be strictly counter-clockwise (20b's
  // foot_fits, moved here); otherwise Hit.
  FootSearch constraint_foot(const LatticeMesh&, std::uint32_t t, MeshVertex p,
                             double delta, const LatticeFrame&);
  ```

  20b's `detail::foot_of` and `detail::foot_fits` (`refine.hpp`) become this;
  refinement keeps only its ε (`foot_epsilon`, 20b R3). The status is
  returned, not folded into "no foot", because the three paths answer a
  near line with no usable foot differently: the quality start skips the
  node (R2.1), refinement and the final check insert the point itself
  (R3, R4). (Round 1's `std::optional<FootHit>` could not tell a skip from
  "no line near"; review round 1 did not raise it, the prototype did.)
- **As built** (`include/terrain/mesh/constraint_foot.hpp@69f37d1c`): the
  header also has `foot_reachable(m, t)`, true when `t`, or a neighbour
  across an unconstrained edge of `t`, has a constrained, non-frozen edge:
  exactly the edges `constraint_foot` can answer from, so when it is false
  the search would return `None`. Refinement uses it to skip the search
  (R3, `@developer`'s performance fix `69f37d1c`). `foot_fits` sits in
  `mesh::detail`, where 15f's `strip_fits` calls it.
- **Not across a constrained edge.** A neighbour behind a constraint is on
  the other side of a line the candidate is not near enough to matter for.
- **Determinism**: a fixed search order and "first found" make the hit a
  function of the mesh and p only.

### R2. The quality start (`quality.hpp`, `improve`)

After the node is snapped, judged valid and located, and before insertion:

1. **Near a constraint**: `constraint_foot(m, t, node, δ_q, f)`. A `Hit` is
   inserted with `split_edge(owner, edge, at)` instead of the node
   (`QualityOutcome::feet`), if the validity callable accepts the foot.
   `NearEnd`, `NotCounterClockwise`, or a foot the callable refuses is a
   skip (`skipped_near_line`). M5: 0 such skips on both catchments, against
   5 930 and 4 233 feet; a node half a cell from a line and its foot within
   half a cell of an input vertex would put that vertex nearly inside the
   bad triangle's empty circumcircle, which is why it is rare.
2. **δ_q = min(dx, dy) / 2**, the cap of 20b R3. The pass reads no height, so
   there is no slope term. Scale: relative to the cell; checked at 10 m
   (Numedalslagen, 1.29 M triangles) and 31 m (Lagan, 0.87 M).
3. **The validity callable takes a `MeshVertex`** (it took a
   `LatticeVertex`): `refine` passes "`vertex_z` has a value", which for a
   node is today's "not NoData", so nodes are judged as before.
4. **Whether the bad triangle goes is measured, not guaranteed** (review
   S5). The node is within half a cell diagonal of the circumcentre and the
   foot within half a cell (δ_q) of the node, so the foot is within one
   cell diagonal of the circumcentre, and the circumradius is at least one
   cell diagonal (20 R5): inside or on the circumcircle, not strictly
   inside. Nor does being inside the circle put the bad triangle in the
   foot's cavity when a constraint lies between them. So no guarantee is
   claimed. What holds: if the bad triangle survives, its slot is unchanged,
   so it is not offered again, exactly as after any skip, and termination
   (step 5) is unaffected. Measured (M5): after every quality-start foot on
   both catchments the bad triangle was gone.
5. **Termination.** Nodes: as 20 R10. Feet: each is at least δ_q from both
   ends of the sub-edge it splits, so a segment of length L holds at most
   L / δ_q of them (20b R6's second argument).
6. **Off with the feet.** `QualityOptions::constraint_feet` (false by
   default); `refine` passes `RefineOptions::constraint_feet`, so
   `--no-constraint-feet` turns all three paths off and the mesh is
   bit-identical to **master run with `--no-constraint-feet`**. It is not
   master's default mesh: master's default already has 20b's refinement
   feet on, and this switch turns those off too. Measured so by `@perf` on
   both catchments, arrays and `.vtk` equal
   (`docs/benchmarks/2026-10-07/20c-1/README.md@b9e2d482:40-50` and
   `:206-213`).
7. **No vertex rule.** Skipping a node closer than δ_q to an off-node
   corner of its triangle was prototyped (`QVERT`) and changed no triangle
   on Lagan. Not designed in (as 20b's C3).
8. **As built** (`include/terrain/mesh/quality.hpp@69f37d1c`): as above.
   A foot is counted in `QualityOutcome::feet`, not in `inserted`, and
   reaches `RefineOutcome::quality_feet`; `skipped_near_line` is summed
   into `quality_skipped` (R5). Measured: 5 930 feet on Lagan and 4 233 on
   Numedalslagen, M5's figures
   (`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/raw/lagan-head-r1_stats.md@c9f62f7c:83`,
   `docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/raw/numed-head-r1_stats.md@c9f62f7c:87`).

### R3. Refinement (`refine.hpp`): 20b R1 widened

`foot_of` searches through `constraint_foot`, so the holding triangle's
three neighbours count too (M3 "R"). A `Hit` in a neighbour splits the
neighbour's edge; the deferral rule of 14 R5 applies to the owner and to
the triangle across the split edge (either touched this round: the node
waits). **When the owner is a neighbour, the holding triangle stays
active** (it goes on the round's `skipped` list): its node is still
unconverged, and nothing else would rescan a slot the split did not write.
`NearEnd` inserts the node, as 20b's `foot_of` returning nothing does
today; `NotCounterClockwise` and a foot with no `vertex_z` count in
`feet_refused` and insert the node, as today. Everything else in 20b (ε,
footed once, the fallback) is unchanged.

**As built, with `@developer`'s performance fix**
(`include/terrain/refinement/refine.hpp@69f37d1c:339-350`). The green
code computed ε (`foot_epsilon`, up to 16 DEM reads) and ran the search
for every split; `@perf` measured the split phase +32 % on Lagan, +51 %
on Numedalslagen and +42 % on the 1 m tile against `--no-constraint-feet`,
the tile placing no foot at all
(`docs/benchmarks/2026-10-07/20c-1/README.md@b9e2d482:264-288`). The fix
(`69f37d1c`), under ruling 5 below, does the foot work only when a
constraint is within reach, in three steps:

1. `foot_reachable(m, t)` (R1) false: no search, no ε.
2. Otherwise the search runs first at the cap, min(dx, dy) / 2.
3. Only when that finds an edge are the footed set and ε asked; the
   search is rerun at ε when ε is below the cap.

The answer is the one the green code gave: ε is at most the cap, so
nothing within the cap means nothing within ε, and when ε equals the cap
the first search is the ε search. The meshes are byte-identical to the
green code's on the 1 m tile, the quarter, Lagan and Numedalslagen, feet
on and off
(`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/README.md@c9f62f7c:158-170`;
the quarter's hash at `:259` equals the green code's at
`docs/benchmarks/2026-10-07/20c-1/README.md@b9e2d482:298`).
Measured: 720 feet on Lagan and 4 629 on Numedalslagen, 0 refused
(`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/raw/lagan-head-r1_stats.md@c9f62f7c:84-85`,
`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/raw/numed-head-r1_stats.md@c9f62f7c:88-89`).

### R4. The final check (`refine_points.hpp`, `point_loop`)

For a winner that is a stored source point or a DEM node (not a strip
point, not a void carve point), Inside its triangle, not yet footed:

0. **L12 first** (review B3). With a strip, the final check already puts a
   winner that lies within r(g) of a constrained edge of its triangle onto
   that edge, at its own position, with its own z
   (`include/terrain/refinement/refine_points.hpp@ed125121:329-334`;
   `near_constraint`, `include/terrain/refinement/strip_scan.hpp@ed125121:199`).
   r(g) is L16's coincidence radius, max(1e-10, 64 ulp of the largest
   lattice coordinate) lattice units: a rounding-scale distance, far below
   δ_p. That rule runs unchanged and first; R4 applies only to a winner L12
   left Inside, so at a distance d with r(g) < d < δ_p. A winner on an edge
   (`where` is an edge) is not footed either. M5: L12 fired 0 times on both
   catchments, R4 793 times on Lagan.
1. `constraint_foot(m, t, point, δ_p, frame)`; with a strip, the split must
   also pass `strip_fits` (the guard every split of a strip sub-edge
   already passes, 15f L1), else `foot_fits` alone. Anything but a `Hit`
   (or a `Hit` refused by `strip_fits` or for want of a z) inserts the
   point itself, as today. `feet_refused` counts what it counts in
   refinement (R5).
2. A hit is inserted instead of the point, with `split_edge`, and the
   strip's sub-edges are cut as for any split of a constrained edge
   (15f step 4). The point is recorded as footed (by its exact position)
   and stays a candidate: if its error is still above the tolerance it is
   inserted later as itself, which is 20b R5's fallback. When the owner is
   a neighbour, the holding triangle stays active, as in R3. **The tolerance
   guarantee is the stored points' and is unchanged**: every stored point is
   still scanned, and the loop stops only when every error is within it.
3. **z of the foot**: on the projected path (`refine_strip`, a raster at
   hand), `vertex_z`. On the reprojected path (no raster), linear along the
   line between the two strip points that bracket the foot on its sub-edge
   (they are the source's own values on that line, 15f); with none on one
   side, the sub-edge's end vertex's z on that side; on a constrained edge
   with no strip record, linear between its two ends. A NaN z refuses the
   foot. **Measured (M5, review B4)**: on Lagan's reprojected path the
   strip's z differs from the point's own by 4.2 m on average and 17.3 m at
   most, and with it the strip needs 111 line points instead of 166, the
   mesh has 642 fewer triangles and 3 fewer slivers, and 19 footed points
   fall back instead of 10; the point's own z would put a wrong height on
   the line, which the strip then has to repair.
4. **δ_p = min(dx, dy) / 2 of the run's frame** (the target grid). Scale:
   relative to the cell; checked at 31 m (Lagan) and 10 m (Numedalslagen).
   Fixed rather than slope-scaled: the reprojected path has no raster to
   take a slope from. 20b measured that a fixed ε fires the fallback more
   often (5 in 150); the fallback count is reported (R5) and is part of
   `@perf`'s acceptance.
5. **Counters**: `PointRefineOutcome::feet`, `feet_fallback` (footed points
   later inserted as themselves), `feet_refused`.
6. `PointRefineOptions::constraint_feet`, false by default; the CLI passes
   the same switch as `refine`.
7. **As built** (`include/terrain/refinement/refine_points.hpp@69f37d1c`):
   as above, with three details the design left open. The wait rule is
   R3's: the point waits a round when the owner or the triangle across the
   split edge was touched this round. A foot counts in `inserted` as well
   as in `feet`, and a footed DEM node's foot is not in `nodes_inserted`
   (green-step ruling 5). The final check's `feet_refused` is counted but
   not reported; the record's "Moves onto lines refused" row
   (`snaps_refused`) is refinement's alone. Measured: 793 feet and 19
   fallbacks on Lagan, as in M5, and 5 and 0 on Numedalslagen
   (`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/raw/lagan-head-r1_stats.md@c9f62f7c:86-87`,
   `docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/raw/numed-head-r1_stats.md@c9f62f7c:90-91`).

### R5. Report and options

- No new flag. `--no-constraint-feet` covers all three paths;
  `elevation_source` keeps `constraint feet on/off`.
- `RefineOutcome::quality_feet` (feet from the quality start) and the final
  check's `feet` and `feet_fallback` reach `--stats`' Refinement table and
  the run record (`points_snapped_to_lines` gains the two other paths'
  counts as separate rows; plain words, per "Plain product output").
- The quality start's skips keep one total, `quality_skipped`, as 20 R11.
- **`feet_refused` has one meaning on both paths that have it**
  (refinement and the final check; review round 2, S2): a foot was found
  and not inserted, so the candidate went in as itself. That is a
  `NotCounterClockwise` status, or a `Hit` refused for want of a z or by
  `strip_fits`. `NearEnd` (the foot too close to an end of the edge) is
  not counted on either path, as 20b's `foot_of` does not count it today.
- **As built** (`src_python/tin_engine/run_record.py@69f37d1c:38` and
  `:71-79`, `src_python/tin_engine/cli.py@69f37d1c:1649-1666`): the record
  and `--stats` gain three rows, "Points moved onto lines while improving
  the starting mesh" (`start_quality_points_snapped_to_lines`), "Points the
  final check against the DEM moved onto lines"
  (`final_check_points_snapped_to_lines`) and "Of them, also added where
  they were" (`final_check_snapped_points_added_anyway`); refinement's row
  now reads "Points refinement moved onto lines", and the switch's row
  "Points very close to a line were moved onto it". `--stats` prints the
  record's rows, so `stats.py` did not change.

### R6. Interactions

- **Frozen edges (23b)**: `constraint_foot` never returns a frozen edge,
  as 20b's search does not today (N5).
- **NoData**: a foot whose `vertex_z` is empty is refused on every path.
- **Strip (15f)**: a foot on a strip sub-edge is a split of it, handled by
  the existing `cut`; strip points are not footed (they are on the line).
- **Determinism**: all three rules run in serial phases, in the existing
  order, and read only the mesh, the candidate and (for z) the DEM or the
  strip. Output bit-identical for any `threads`.

## Design of PR 20c-2: the soft criterion and the line split in the quality start

### R7. A candidate pays for itself in worst angle

- θ (`--start-min-angle`, default 25°) stays the **trigger**: only a
  triangle under θ makes a candidate, as today.
- **Acceptance**: before inserting a candidate (node or foot), `improve`
  computes its **cavity** read-only: the triangles whose circumcircle holds
  the candidate, grown from the triangle it lies in (and, for a foot, the
  triangle across the split edge) without crossing a constrained edge, with
  the exact `incircle` in the frame that `legalise_around` uses. The new
  triangles are the candidate joined to the cavity's boundary edges. Let
  `old` be the smallest angle among the cavity's triangles and `new` among
  the new ones, both capped at θ. **Insert only if `new ≥ old + P − s`**,
  with the slack s = 10⁻⁶°; otherwise skip (`skipped_no_gain`).
  **Why the slack** (pins ruling, 20c-2's red step, observation (a)): a
  foot on a constrained edge keeps the angle at that edge's end exactly
  (the same two rays), so where that angle is the cavity's worst, `new`
  equals `old` in exact arithmetic and a bare `≥` is decided by the foot's
  last bit, which a fused multiply-add can move between the CI platforms.
  The slack puts the threshold away from that equality. Scale: the
  difference rounding makes is about a few ulp of the largest coordinate
  over the foot's distance to the end (≥ δ_q, half a cell): four ulp on a
  grid of 2¹⁷ cells a side (the São Francisco basin spans about 1 550 km
  north to south, about 5.2 × 10⁴ cells at 30 m) give
  1.3 × 10⁻⁸°, so s has a margin of about 75 there; derived, not
  measured, and checked on no real input yet (the largest fixture is 33
  cells a side). A candidate that lowers the
  worst angle by less than s is accepted; R7 only refuses, so accepting
  more cannot break termination.
- **P = 0° by default** (Ola's ruling, question 2): a candidate is refused
  only if it would leave its neighbourhood's worst angle lower than before.
  M5: 12.6 % and 16.8 % fewer triangles than 20c-1 at about the same sliver
  count (694 against 679, 440 against 433); 32 499 and 63 383 candidates
  refused. P = 2° was measured too and is
  far too strict (8 437 slivers on Lagan). The objective matters: the summed
  angle deficit instead of the worst angle refused the insertions that
  break up caps (1 393 slivers at P = 0).
- **Why the cavity is exact.** The mesh is constrained Delaunay before each
  insertion (14b R10), so Lawson legalisation after inserting a point gives
  exactly the constrained Bowyer-Watson cavity's triangulation. The read-only
  prediction and what `legalise_around` produces are the same triangle set;
  `@tester` asserts it (T-P1).
- **Termination**: the rule only refuses; 20 R10 and R2 above hold.
- **No exemption for slivers.** Accepting every candidate whose cavity
  holds a sliver, whatever the gain, was measured (M5) and changes 0 and 6
  slivers; not designed in.
- `QualityOptions::min_gain_deg` (P), negative = off (today's hard rule);
  `RefineOptions` and the binding carry it; the CLI passes 0.
- **CLI**: `--start-quality-gain DEG`, default 0; `-1` restores the hard
  rule (bit-identical to 20c-1). Refused: non-finite, above 10.
  The setting is an input row of the run record and `--stats`,
  `start_quality_gain_deg`, beside `start_min_angle_deg` (printed `0` by
  default). (This line first said `elevation_source` gains a clause;
  increment 25 replaced that sentence by such rows, so it is not written.)
- Reported: `start_quality_points_skipped` keeps its total, which holds
  R7's refusals; a new row `start_quality_points_without_gain`, "Tries
  that would not have improved the angles", in `--stats` and the record,
  right after `start_quality_points_skipped`.

### R8. A quality point beyond a constraint line splits the line (Ola's yes, question 4)

- **When**: the walk from the bad triangle to its node (20 R4) stops at a
  constrained edge e of triangle t that has a triangle beyond it (an outline
  edge with nothing beyond is still `skipped_blocked`), and e is not frozen.
  The trigger is this blocked walk, not Ruppert's encroachment test (the
  node lying in e's diametral circle, the circle with e as its diameter):
  a node can lie beyond e yet outside that circle (near one of e's ends),
  or inside the circle on the near side, where the walk does not cross e. The blocked walk is what the
  prototype measured (M5, B), and it needs no extra geometry; no
  diametral-circle test is added.
- **What**: the node's foot on e (its orthogonal projection, in the world
  frame) goes in with `split_edge(t, e, foot)`, if the foot is at least δ_q
  (R2.2) from both ends of e, both children on each side are strictly
  counter-clockwise (`foot_fits`), and the validity callable accepts the
  foot; otherwise the skip stays `skipped_blocked`. Counted in
  `QualityOutcome::line_splits`; the new triangles are offered to the queue
  as after any insertion.
- **Not judged by R7's gain test**, as measured: the split answers the
  line, not the bad triangle alone, and the gain test would see only the
  two triangles beside e.
- **Only with the soft criterion, and only with the feet.** On when R7 is
  (gain ≥ 0) **and** `constraint_feet` is on (pins ruling, 20c-2's red
  step, pin 4): the split puts a point on a line at a foot, and
  `--no-constraint-feet` is the one switch that keeps every point off the
  lines (R2.6, R5). M5's B was measured with the feet on, so the gate is
  unchanged. Gain −1 turns both off and gives 20c-1's mesh bit for bit. Under the hard 25°
  rule the split costs +15 % (Lagan) and +51 % (Numedalslagen) triangles
  (M5), which is why it is not in 20c-1.
- **Termination**: every split is at least δ_q from both ends of the
  sub-edge it cuts, as R2's feet are, so a constraint segment of length L
  takes at most L / δ_q vertices from the two rules together.
- **Measured** (M5): 14 025 splits on Lagan and 20 572 on Numedalslagen;
  slivers 437 and 298 (20c-1: 679 and 433); worst angle 0.0024° (unchanged)
  and 0.0086° (20c-1: 0.00068°); triangles 7.5 % and 11.4 % under master's.
- **The constraint stays the same set of lines**: `split_edge` keeps the
  bit and the mask on both halves (20 Q3's check), as for every foot.
- Reported: a row `start_quality_lines_split`, "Land-cover and outline
  lines split to improve angles", in `--stats` and the run record, right
  after `start_quality_points_snapped_to_lines`.

## Design of PR 20c-3: input repair and coarsening (rulings 1 and 9), and the outline rule (ruling 6)

Brought up to date by `@architect`, 2026-10-08, on branch
`worktree-soft-quality-3` (from 20c-2's approved head `68cc96b6`), for
Ola's ruling 9 (the CORINE slit is an input error to repair; it amends
ruling 5) and for what 20c-2 and increment 30d (PR #215) changed.

**Speed judgment** (Ola's rule for every design): **yes, this changes
the work a run does.** It adds a Python stage to every run with a land-cover
source (the clip and the repair: 1.1 s per catchment in M7's probe, about
+6 % of Lagan's 17 s run (30d's figure) and +17 % of Numedalslagen's 6.6 s),
and the outline rule (on by default) changes the mesh and so refinement's
work; the outline rule itself has not been timed. It also drops the old
design's coverage check, which M7 timed at 15.8 s on each catchment. So
`@perf` runs a short timed check before the review: at most 15 minutes,
Lagan and Numedalslagen once each with the defaults, against master, and
the `features clip` row of `--stats` read; not a sweep.

### Where it sits

Python, in the feature reader (`feature_input.py`), never in `_core`
(the I/O boundary, `CLAUDE.md` §2). It is a **land-cover stage** between
reading a source and clipping its lines to the domain, and it works on one
source at a time: a source is one coverage, and two sources (say CORINE and
a lake file) are not one partition, so a gap between their polygons is not
repaired.

- **Which features.** The polygons of a source whose class map names a
  code system (`ClassMap.codes`, 16c R5: `corine`, `corine-water`,
  `clc18_kode`), after the map's drops. Lines, and polygons under a map
  without codes (`property`), go through today's path unchanged.
- **In the computation CRS.** Each polygon is moved whole, vertex by
  vertex, as `_take` already does for the label polygon (`polygon`), and
  the stage works on those. A land-cover source therefore always takes the
  "moved first, then pre-clipped exactly in the DEM's CRS" route (16b R5's
  second route, `w = 0`, exact), never the geographic pre-clip before the
  move: the repair has to see whole borders in one CRS.
- **Clipped to the read region first.** Each moved polygon is intersected
  with `region` (16b R5's region in the DEM's CRS, which `_take` already
  builds: since 30d, PR #215, merged as `769ed992`, the domain's convex
  hull grown by 100 m). The region is convex and lies at
  least 100 m outside the outline, so the new edges the clip makes are
  outside the domain and farther than any outline-rule D (below) from the
  outline; only polygonal parts are kept. **The clip also shortens edges**
  that cross the region's edge, and such an edge can still reach within D
  of the outline (an input edge longer than about 95 m). The rule therefore
  places its cuts on an edge at multiples of D along the edge's line,
  counted from the foot of the CRS's origin on that line, not from the
  edge's end, so where the clip cut the edge does not move them (as built,
  `7334b82b`; corrected in "Rulings on 20c-3's green step", ruling G4, after
  30d's `test_grow.py::test_7b` caught the first build), **less any
  multiple within D/2 of either end of the edge** (ruling T1: a cut
  centimetres from an edge's own vertex made a centimetre constraint
  edge). M7: this takes Lagan's land cover from 1 006 098
  to 328 372 vertices and the repair from 3.1 s to 0.75 s (the clip 0.36 s).
  The label polygon of each feature is then the clipped, repaired one; the
  region holds the domain, so every triangle's label inside it is found as
  before.
- **Steps, in this order**, the first four on the polygons, the last
  turning them into lines:
  1. **Repair** (`--features-repair METRES`, default 1, ruling 9 and
     question 9): `shapely.coverage_clean(polygons, snapping_distance=S,
     gap_width=S)`.
  2. **Same-class borders dropped** (`--features-merge-same-class`, on;
     `--no-features-merge-same-class` turns it off; ruling 1 (a)):
     the polygons of one class are unioned (`coverage_union_all` per
     class), so a border with no change of class is no constraint. Removes
     population 3's lines (E12). As built: each class becomes one feature
     under the first fid of its class in source order, a MultiPolygon where
     its parts do not touch (ruling G6 (d)).
  3. **Horizontal tolerance** (`--features-tolerance METRES`, default 0,
     off; ruling 1 (b)): `coverage_simplify` at that tolerance,
     `simplify_boundary=False`, as measured (M5, M6), no 1 cm snap (M6 (a)).
     The clip changes the border chains it simplifies only where a border
     crosses the read region's edge, 100 m or more outside the domain;
     the gates' margins are taken to cover that (not measured).
  4. **Ola's outline rule** (`--features-outline-snap METRES`, default 5
     whenever features are given; 0 turns it off; ruling 6): M6 (c)'s four
     steps, unchanged below.
  5. **Lines**: each polygon's rings (or the outline rule's linework) go
     through today's pre-clip and clip to the domain.
- **Off.** With `--features-repair 0 --no-features-merge-same-class
  --features-tolerance 0 --features-outline-snap 0` the stage is skipped
  entirely and the features are byte for byte today's (the off switch the
  gates use). Any one of the four on runs the stage, with the repair at
  whatever S says; S = 0 still runs `coverage_clean` at zero snapping and
  zero gap width, which makes the coverage valid for steps 2 and 3 and on
  a valid coverage moves nothing (M7: 0.00 m² on Numedalslagen).
- **Recorded.** Input rows: `features_repair_m`,
  `features_merge_same_class` (the text `on` or `off`, as `snap_to_lines`;
  ruling G3), `features_tolerance_m`,
  `features_outline_snap_m`. Result rows: the land-cover vertices after the
  clip and after the stage, "Land-cover vertices before and after clean-up";
  for the outline rule, the land-cover area inside the outline that changed
  class (the count of borders moved, designed before, is dropped: ruling
  G6 (g)). A `features clip: clean-up` sub-row under
  `features clip` in the timing table, so the timed check reads the stage's
  own seconds. The vertex, area and timing rows appear only when the stage
  ran. The four flags sit in their own help panel, "Land-cover clean-up"
  (ruling G7).
- **Dependency.** `shapely>=2.2` in `pyproject.toml` (2.2.0 is on PyPI;
  its wheels carry GEOS 3.14, which `coverage_clean` needs). Not a new
  dependency, and not on `CLAUDE.md` §2's list. CI installs with `pip
  install -e .`, which already resolves to the newest shapely; local
  venvs (2.1.2 with GEOS 3.13 today) must be upgraded. The bump is its own
  commit by `@developer` before `@tester`'s red step, with the full pytest
  run on it; any test that changes from the bump alone goes to `@architect`
  before the red step, not into 20c-3's tests.

### Step 1: the repair (ruling 9)

**What it does.** GEOS's coverage cleaner (Davis 2025; prior art below)
nodes the polygons' linework with a snapping noder at distance S, rebuilds
the faces, gives each overlap to a neighbour, and fills each gap narrower
than S (twice its largest inscribed circle's radius) into the neighbour
with the longest shared border. The result is a valid coverage: where two
neighbours' borders came within S of each other, they are now one border.

**The tolerance S: two borders are the same border when they come within S
of each other.** One number for both of the cleaner's tolerances (the
snapping distance and the gap width), because they answer the same
question for a partition that is meant to have no gaps. **Default 1 m**
(question 9 as asked again with M8, Ola's yes 2026-10-08). M7 set the
first default, 5 cm, the smallest distance probed that joins the slit's
corners (the first three points below); M8 then probed 25 cm and 1 m on
the built code with the defaults, and 1 m is the probed distance that
removes the thin triangles (the last two points):

- Lagan's slit has its two open ends 1.005 cm apart. Snapping at 1 cm (with
  a 1 m gap width, M7) fills the gap but leaves its two corners apart, so
  its 1 cm closing edge stays a constraint; snapping at 5 cm makes the
  corners one vertex and the two 79 m borders one.
- Below 100 m wide, the two catchments' land cover has no other gap (M7):
  Lagan one (the slit), Numedalslagen none. So any gap width from 5 cm to
  100 m fills the same gap here. The snapping distance moves more as it
  grows (1 cm: 0.64 m² and 495 m²; 5 cm: 0.72 m² and 1 085 m²; nothing
  larger was probed), so 5 cm is the smallest distance probed that joins
  the slit's corners.
- What else 5 cm changes, measured: ring edges shorter than 10 cm go from
  356 to 5 (Lagan) and from 464 to 1 (Numedalslagen); this is population 1,
  M2's millimetre vertex pairs, the same kind of error (borders that should
  be one, millimetres apart). Land cover moved: 0.72 m² and 1 085 m² of
  5 548 to 6 441 km², about 2 × 10⁻⁷ of the area.
- At 1 m (M8, on the repair's own input after the clip): slivers under 1°
  with the defaults 180 to 26 (Lagan) and 158 to 2 (Numedalslagen) against
  5 cm, worst angle 0.19° to 0.73° and 0.63° to 0.83°; 25 cm gives what
  5 cm gives. What 1 m joins is borders of two polygons that run within a
  metre of each other, 157 of Lagan's 180 slivers and 156 of
  Numedalslagen's 158 at 5 cm lying within 2 m of land cover it moved;
  ring edges under 10 cm go to 0 on both.
- What 1 m moves: about 45 km and 123 km of border, by under 1 m, some
  3 400 m² and 8 500 m² of land cover inside the domains given to a
  neighbour (5 × 10⁻⁷ and 1.5 × 10⁻⁶ of the area), 1 % of CORINE's 100 m
  minimum width; no measurable time.

Scale: S is in metres in the computation CRS; the 1 m default assumes a
land-cover source whose borders carry no meaning below a metre (CORINE:
1:100 000, 25 ha minimum mapping unit, 100 m minimum width) and DEM cells
of 10 to 31 m; checked on the two catchments' CORINE, about 380 000 and
600 000 vertices after the clip (about 1 million each before it), the
largest inputs probed (M8). On a source mapped to the decimetre it moves
real borders; S is a flag for that reason.

**A real gap wider than S stays.** It keeps its two borders as two
constraints, and the land cover there has no class, as today. Step 3 does
not simplify its sides (`simplify_boundary=False` keeps edges on a
coverage's boundary, and a gap's sides are on it; RP6 below checks this).
A gap open to the read region's edge is not a face of the cleaner's, so by
its design (Davis 2025) the gap width does not fill it and only the
snapping acts there; not tested. To be open there, a gap must run at least
100 m past the outline; none was found (M7).

**Against ruling 1 (b) (simplification off by default).** The repair is not
the 2 m simplification: it moves a vertex only onto another border within
S (1 m by default), and it is on because Ola called that kind of mismatch
an error (ruling 9) and set its default (question 9). The
simplification at a distance remains Ola's flag, off.

**How it meets the other steps.**

- **Coarsening (step 3).** Repair first: `coverage_simplify` needs a valid
  coverage (M6 (a): one that is not breaks), and a repaired slit is an
  inner border, simplified like any other. M5's mechanism, the 2 m
  simplification dropping the vertex beside the slit's 1 cm edge (0.088° to
  0.006697°), has no edge left to act on. The old design's coverage check
  before simplifying (`coverage_is_valid`, with its plain error) is
  dropped: `coverage_clean`'s output is a valid coverage by its contract,
  and the check cost 15.8 s per catchment (M7).
- **The outline rule (step 4).** S (1 m by default) is a fifth of D (5 m),
  under D/2, and the repair puts no point on the outline, so the rule's rounding (points it
  places are an outline vertex or more than D/2 apart) is unchanged. Repair
  first, so the rule never sees two copies of one border.
- **The same-class merge (step 2).** `coverage_union` per class needs a
  valid coverage; the repair gives one.
- **The noder.** Unchanged; its `--snap-spacing` (1 mm) stays, and the
  noder's failure at coarser spacing (M4, its own ROADMAP row) is not
  needed.

**Not guaranteed.** The snapping noder is not iterated snap rounding:
nothing bounds the input's local feature size below by S, and a short ring
edge can survive (at 5 cm, on Lagan, 5 of the 356 edges under 10 cm; at
1 m none of the edges under 10 cm on either catchment, M8). The gates below ask for what was
measured, not for a minimum edge length.

### Steps 2 to 4

As designed in rounds 1 to 6, with three changes: the first from the
repair, the other two from 30d:

- **Step 3's coverage check is gone** (above). `--features-tolerance`
  with `--features-repair 0` is allowed: step 1 then runs at zero
  tolerance, which makes the coverage valid.
- **The outline rule never buffers the outline.** 30d found GEOS's buffer
  of Lagan's raster-traced outline (a staircase) to be 98 % of decode. The
  rule finds the border vertices within D of the outline with an
  `STRtree` of the outline's segments (`query(..., predicate="dwithin",
  distance=D)`) and moves each to its nearest point on the nearest
  segment; no `buffer`, no `snap`. D is refused at or above 100 m (the read
  region's margin), so the clip's new edges stay out of the rule's reach;
  the edges it shortens are cut along their line (ruling G4).
- **Against the domain as the mesh gets it** (after increment 22's
  reduction, if any), every ring; 30d grows only the read region, never
  the domain itself, so this is unchanged.

Step 4 in full (unchanged): M6 (c)'s four steps on the polygons' rings, in
the computation CRS. M6 (c)'s rounding (its step 2) then rounds once more:
a multiple of D within D/2 of an outline vertex goes to the vertex (the
prototype did not), so two distinct points the rule puts on the outline are
more than D/2 apart unless both are outline vertices. The result is
linework: rings cut where they leave the outline, the stretches on it
dropped, each piece keeping its polygon's class. D's scale: metres in the
computation CRS, measured at 5 and 10 m on DEM cells of 10 and 31 m (M6
(d)); about half a cell or less is what was measured.

### What 20c-2 built, against what 20c-3 assumed

- **R7, R8 and pin 4.** 20c-3's figures were measured on the prototype with
  20c-2's switches (`run5.sh`: soft criterion at gain 0 and the line split,
  feet on), which is what 20c-2 built (R8 only with the feet on, pin 4).
  The built 20c-2 matches the prototype on Lagan (437 slivers) and is 2
  under it on Numedalslagen (296 against 298). Nothing in 20c-3 depends on
  R7 or R8 beyond the mesh it is measured on; the gates below are measured
  against 20c-2 as built (master after PR #214).
- **Ruling 8** (the split phase judged per split) is a time limit on C++
  refinement; 20c-3 changes no C++ and has no `bench.py` run. Its time is
  the speed judgment above.
- **30d (PR #215, merged as `769ed992`, in this branch's base)**: the
  read region is the domain's convex hull grown by 100 m, so more polygons
  are read than before. The clip to that region (above) is what keeps the
  repair's cost down; the old design's "filter to the polygons that reach
  the domain" (for the coverage check) is gone with the check.

The domain outline itself is not touched (increment 22 has its own
reduction). Measured in M3 to M7 on Lagan and Numedalslagen. The repair and
the merge change the default mesh of every run with a land-cover map, and
step 4 every run with features, so 20c-3's off switch for the gate is all
four off.

## PRs and gates

| PR | what | needs | gate (measured value in brackets, M5) |
|---|---|---|---|
| 20c-1 | R1 to R6: the foot rule on the quality start, refinement's neighbours and the final check | nothing | Lagan: share under 1° ≤ 0.09 % (0.078 %), triangles ≤ +1 % of master's (+0.43 %), worst angle ≥ master's 0.000412° (0.0024°); Numedalslagen: share ≤ 0.04 % (0.034 %), triangles ≤ +1 % (+0.27 %), worst angle ≥ master's 0.000399° (0.00068°); `--no-constraint-feet` bit-identical to master run with `--no-constraint-feet` (master's default has 20b's feet on, so it is not the comparison); tolerance oracle; determinism. **Built: every gate passes, at M5's figures** (`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/README.md@c9f62f7c:186-201`); time: the split-phase limit below |
| 20c-2 | R7 and R8: the soft criterion and the line split | 20c-1 merged | on both catchments: triangles ≤ 0.95 × master's (0.925, 0.886); sliver count ≤ 0.75 × 20c-1's (0.64, 0.69); worst angle ≥ 0.95 × 20c-1's (Lagan 1.000, the same triangle; Numedalslagen 12.7 ×); `--start-quality-gain -1` bit-identical to 20c-1. **Built: every gate passes** (triangles 0.9253 and 0.8865; slivers 437 and 296, 0.64 and 0.68; worst angle 1.000 and 12.68 ×; gain −1 equal to 20c-1 bit for bit; `docs/benchmarks/2026-10-07/20c-2/fix-65ae0792/README.md@5b0afe04:20-28`); time: the split-phase limit, judged per split for 20c-2 (ruling 8) |
| 20c-3 | the repair, input coarsening, and the outline rule (repair on at 1 m, ruling 9 and question 9; outline rule on at 5 m, ruling 6) | 20c-2 merged | measured without the merge and, for the rows M5 and M6 measured, without the repair, every flag given explicitly, so the gate does not depend on a default. **`--features-tolerance 2` alone** (`--features-outline-snap 0 --features-repair 0`): Lagan: slivers with a side under 10 cm ≤ 25 (19, from 185), share under 1° ≤ 0.03 % (0.020 %), worst angle ≥ 0.005° (0.006697°); Numedalslagen: share ≤ 0.006 % (0.003 %), worst angle ≥ 0.008° (0.012446°). **With the outline rule at 5 m on top**: Lagan: sliver count ≤ 90 (61), at most 5 with the centre within 20 m of the outline (0; 83 without the rule), triangles ≤ +1 % of the tolerance-only mesh (+0.09 %), worst angle ≥ 0.95 × the tolerance-only figure (1.000, the same triangle); Numedalslagen: sliver count ≤ 10 (4), at most 5 within 20 m of the outline (0; 27 without), triangles ≤ +1 % (+0.21 %), worst angle ≥ 0.1° (0.832°); no shared border inside the catchment left unmatched by the rule. With the merge on, the two population-3 lines carry no land-cover line between two polygons of the same class (M2's bands; borders between different classes that lie within 1 m of the lines are real and stay, R3's outcome). All four off give 20c-2's mesh bit for bit. **The repair at 5 cm on top of each of those two runs** (5 cm given explicitly, the default when the gate was set; the 1 m default's runs are M8's, below): every threshold of the run without it still met; on Lagan the slit's two open corners (413 602.500 6 331 638.254 and 413 602.501 6 331 638.244) are not both mesh vertices, and no triangle with its smallest angle under 1° has a corner at the slit's far corner (413 677.295 6 331 664.081). The thresholds are not raised to the measured figures (R3: they are the design's figures already) |

**20c-3's gate results so far** (`@perf`, on `745bcb0d`,
`docs/benchmarks/2026-10-08/20c-3/README.md@16143f4c`; "near" is the
sliver's centre within 20 m of the domain outline, the same measure as
the prototype's `strip.py`, checked by reading both scripts). Tolerance
2 m alone: met on both. Outline rule at 5 m on top: Lagan met (65 slivers,
4 near, worst 1.000 × the tolerance-only figure); Numedalslagen **not met**
(15 slivers against ≤ 10, 11 near against ≤ 5, worst 0.0814° against
≥ 0.1°; triangles +0.28 %, met). Repair at 5 cm on top: no threshold that
held without it fails; Lagan's slit corners are not both mesh vertices;
the far-corner check, the shared-border check and the population-3 check
not run. All four off: 20c-2's mesh bit for bit on both. The cause of
Numedalslagen's miss and its fix are ruling T1 below. **After T1 and T2**
(`@perf`, on `8716f062`, `docs/benchmarks/2026-10-08/20c-3/README.md@41db7333:91-130`):
outline rule at 5 m on top, Numedalslagen **met** (4 slivers, 0 near,
worst 0.8322°, triangles +0.27 %), Lagan met (61, 0 near, worst
0.006697°, 1.000 × the tolerance-only figure, +0.09 %): the design's
figures. T2's mesh is T1's byte for byte (`.vtk` SHA-256 Lagan
`ef67e65c…`, Numedalslagen `1f11f291…`). The repair runs above predate
T1, and the shared-border, population-3 and far-corner checks were never
run: what is run and what is ruled out is "Rulings on 20c-3's code review
round 1", R3.

**R3's run** (`@perf`, on `aaee898b`,
`docs/benchmarks/2026-10-08/20c-3/README.md@63065a15:132-183`, raw output
in `stats/r3/`): (a) the default meshes unchanged by R1 (Lagan
`ef67e65c…`, Numedalslagen `1f11f291…`): met. (b) the repair on top of the
outline rule: met on both, Lagan 743 614 triangles (+0.09 %), 59 slivers,
0 near the outline, worst 0.2320°; Numedalslagen 1 016 194 (+0.26 %), 4,
0, 0.8322°; `corners.py` on Lagan: the slit's open corners not both mesh
vertices, met. (c) the slit's far corner: met (6 triangles meet there,
none under 1°). (d) shared borders: 0.000 m on both, met; one shared vertex
moved 0.5 m gives 208.170 m and 337.081 m. (e) met on the corrected
wording (R3's outcome, in "Rulings on 20c-3's code review round 1"):
26.194 m with the merge on, all of it borders between different classes;
161 164.495 m with it off.

Slivers with the defaults (every switch at its default; `@perf`, same
commit), master `483221d2` against the branch:

| catchment | master: triangles | slivers under 1° | worst | branch: triangles | slivers under 1° | worst |
|---|---|---|---|---|---|---|
| Lagan | 799 378 | 437 | 0.002384° | 799 408 | 186 | 0.1917° |
| Numedalslagen | 1 141 207 | 296 | 0.008623° | 1 142 108 | 169 | 0.08145° |

After T1 and T2 (`@perf`, on `8716f062`, same README at `41db7333`),
with the repair's default then at 5 cm: Lagan 799 368 triangles, 437 →
180 slivers under 1°, worst 0.002384° → 0.1917°; Numedalslagen
1 141 983, 296 → 158, worst 0.008623° → 0.6302°.

**With the repair's default at 1 m** (Ola's answer to question 9; M8's
1 m row, `@architect`, on `9a04c3c1`, one run each, AC power,
`docs/benchmarks/2026-10-08/20c-3/m8/stats/m8-s1-lagan.md` and
`m8-s1-num.md`; the code at `c71ccf57` differs only in the default, so
these are its default meshes): Lagan 798 554 triangles, 437 → 26 slivers
under 1°, worst 0.002384° → 0.7293°; Numedalslagen 1 118 006, 296 → 2,
worst 0.008623° → 0.8322°. Not re-timed by `@perf`; M8's single runs
give the land-cover clean-up 5.98 s and 5.46 s and the whole run 21.62 s
and 10.82 s (`--stats` total).

**Why counts, not shares, for 20c-2** (review B1). The soft criterion
removes triangles where the angles are already fine, so even at an equal
sliver count the share under 1° rises by the inverse of the triangle ratio:
with the soft criterion alone, counts are 1.02 × 20c-1's on both catchments
but shares 1.17 × and 1.22 × (M5), and round 1's gate of 1.15 × the share
failed on its own figures. The gate is therefore on the count, with the
triangle count gated separately.

**Why 20c-2's worst-angle gate on Lagan sits at its measured value, and
the margin** (review round 2, B1). Equality is expected: Lagan's worst
triangle after 20c-1 and after 20c-2 is the same triangle, corner for
corner (rerun, round 2: the 3.2 mm CORINE segment and a quality node 75 m
away, M5 "Which triangle sets the worst angle"). Its circumradius (38 m) is
under the quality start's floor (43.8 m), so it is never a candidate, and
neither R7 (which only refuses) nor R8 (which acts only on a blocked walk
from a candidate) can touch it. What can move it is its far corner: the
corner is whatever the rest of the mesh puts there, and a production
difference in which quality nodes go in nearby could change it, as the
soft criterion alone moved Numedalslagen's (M5). The gate therefore allows
5 %, 0.95 × 20c-1's, which still forbids the Numedalslagen kind of fall (a
factor 1.71). **The mechanism's own bound is looser** (review round 3, S1):
the needle stays out of the quality start's reach only while its
circumradius is under the 43.8 m floor, so its far corner stays under
about 87.6 m (twice the floor) and the angle above about 75 / 87.6 ≈ 0.86 ×
0.002384°. The gate's 0.95 is tighter than that: a far corner moved beyond
about 79 m (75 / 0.95) fails the gate although the mechanism allows it, and
that is the case the diagnosis below is for. **If the figure misses, even by a hair**, `@perf` does not
rerun or round: it reports the worst triangle (its corners, side lengths,
and which path put each corner in, with the kept prototype's driver,
`../rasputin_scratch/20c-prototype/e7.py` and `worst.py`) and the PR goes
back to `@developer` for a diagnosis, with `@architect` told; only a
design round changes the gate.

**The fallback for a no to question 4 is retired** (Ola said yes). It
had the tightest margin in this file: a sliver count of 1.022 × 20c-1's
against a bound of 1.05 ×, 2.7 % of headroom.

The share is read from `--stats`' quality table ("< 1°"), the worst angle
from the same table, and the counts from the same file; no new tool.
`@perf` runs both catchments for each PR on AC power and keeps the
`_stats.md` files under `docs/benchmarks/<date>/`, with the commands. The
thresholds are M5's figures with a margin for the production code
differing from the prototype; every rule in a gated PR is in the measured
figure, R4.3's z included.

## Tests `@tester` writes red first

### 20c-1

**Invariant-critical (mutation testing required): `test_constraint_foot`**
(new, C++).

- CF1 **the helper**: on hand-built `LatticeMesh` fixtures, a point at
  0.01, 0.4 and 0.6 cells from a constrained edge of its triangle, and from
  a neighbour's: `Hit` and owner as R1 orders them; `None` at 0.6
  (δ = 0.5); `NearEnd` when the foot is within δ of an end;
  `NotCounterClockwise` when a child would not be; `None` across a
  constrained edge, on a frozen edge, and for a point exactly on the edge;
  the foot is the orthogonal projection (to 1e-12 cells), not the midpoint.
- CF2 **quality start**: a bad triangle whose snapped node lies 0.01 cells
  from a long constrained edge: the output has the foot on the edge and not
  the node; the bad triangle is gone; the constraint edges as a set of lines
  are unchanged; bits and masks on both halves (20 Q3's check). Same with
  the constrained edge on the neighbour.
- CF3 **refinement through a neighbour**: a worst node whose holding
  triangle has no constrained edge but a neighbour does, within ε: the foot
  goes on the neighbour's edge, and the holding triangle is rescanned next
  round (its node still goes in if its error stays above the tolerance);
  with that neighbour touched this round, the node waits a round.
- CF4 **final check**: a stored point 0.02 cells from a constrained edge
  with error above the tolerance: the foot goes in, z as R4.3 (a strip whose
  profile differs from the point's z tells the two apart); at
  `tolerance == 0` the point goes in after its foot (the fallback), the run
  ends, feet ≤ footed points; the tolerance oracle over all stored points
  passes. With a strip, a point within r(g) of the edge goes in by L12 at
  its own position and z, not as a foot (R4.0); a point on an edge is not
  footed.
- Mutants to kill: the end check removed (CF1, and a foot a hair from an
  input vertex); the neighbour search crossing a constrained edge (CF1);
  midpoint instead of foot (CF1); footed-once removed (CF4's bound); the
  fallback removed (CF4's oracle); the foot's z from the point (CF4); the
  holding triangle dropped from the active set (CF3).

Property and integration:

- **Off switch**: `constraint_feet = false` on all three paths gives
  output bit-identical to master's with `constraint_feet = false` (not
  master's default, which has 20b's refinement feet on) on 14b's T12
  fixtures, a domain start and a features start.
- 14b's T3 (tolerance) and T6 (determinism, threads 1, 2, 7 and hardware
  concurrency) re-run with the rule on; T6 is the TSan job's.
- CLI: the new rows in `--stats` and the record; `--no-constraint-feet`.
- Existing pins that change because the quality start and the final check
  now foot (`test_mesh_quality` counts, golden hashes over scenes with
  constraints): each listed with the reason and amended in its own commit,
  as the NoData fix did (20, Q-V3).

#### Pins ruled for 20c-1's red step (`@architect`, 2026-10-07, on `a982d531`)

`@tester`'s red commit fixed some details the design left open. Each is
ruled here; "kept" means `@developer` builds to it as written.

1. `constraint_foot`, `FootStatus` and `FootSearch` live in namespace
   `terrain::mesh`. Kept.
2. "Closer than δ" is strict, and the distance is to the closed segment in
   the world frame (col · dx, row · dy). Kept: it is 20b's `foot_of`
   (`>= eps` is "not near"; the projection clamped to the segment), and
   `LatticeFrame` is that frame.
3. `at` is the orthogonal projection in the world frame; owner, edge and
   `at` are asserted only for `Hit`. Kept; for the other statuses they stay
   unspecified.
4. The `NotCounterClockwise` fixture folds a child for every double within
   64 ulps of the projection in each coordinate. Kept: it makes the case
   independent of how `@developer` rounds the projection.
5. `QualityOutcome::inserted` counts DEM nodes only and quality-start feet
   are in `feet` alone; through `refine`, vertices = start + `quality_inserted`
   + `quality_feet` + `inserted`, and refinement's own feet stay inside
   `inserted` as in 20b. Kept: a foot is not a DEM node, which is what
   `inserted`'s comment promises.
6. CF2's exact counts (`feet == 1`, `inserted == 0`) on its two fixtures,
   taken from running today's `improve` on the mesh the foot leaves. Kept.
7. The final check's feet are counted in `inserted` (so vertices = start +
   `inserted` on every final-check path), and `feet_fallback <= feet`. Kept,
   as refinement counts its feet.
8. R4.3's foot z asserted to 1e-9 m. Kept; scale: z under 100 m on these
   fixtures, so 1e-9 m is about 7e4 ulps of headroom.
9. CF4's `ExactStore` through `refine_points`' `Store` template. Kept:
   `refine_points.hpp` names "a test double with geometry(), frozen() and
   for_each_in" as a Store.
10. Python names: `RefineOutcome.quality_feet`, `PointRefineOutcome.feet_fallback`,
    keyword `constraint_feet=False` on `_core.refine_points` and
    `_core.refine_strip`; the CLI passes the switch to the final check on
    both paths (from `edge_strip.py` and `final_check.py`, which the tests
    spy on by those modules' names). Kept.
11. The three row names `start_quality_points_snapped_to_lines`,
    `final_check_points_snapped_to_lines` and
    `final_check_snapped_points_added_anyway`. Kept: snake_case keys in the
    pattern of `start_quality_points_inserted`. Ola's plain-output rule binds
    the wording `--stats` prints beside each key (`run_record.WORDING`),
    which the tests do not pin, so it is set here:
    - `start_quality_points_snapped_to_lines`: "Points moved onto lines while
      improving the starting mesh", right after `start_quality_points_skipped`;
    - `final_check_points_snapped_to_lines`: "Points the final check against
      the DEM moved onto lines", after `snaps_refused`;
    - `final_check_snapped_points_added_anyway`: "Of them, also added where
      they were", right after it;
    - with three such rows, the existing two are reworded (keys unchanged):
      `points_snapped_to_lines` "Points refinement moved onto lines", and
      `snap_to_lines` "Points very close to a line were moved onto it"
      (it is no longer DEM nodes only).
12. Four C++ targets instead of one. Kept: CF1 starts no threads and stays
    out of the TSan job; the other three do. The design's suite name
    `test_constraint_foot` names all four for mutation testing.
13. R3's "the node waits a round" is not visible in the output; the
    deferral fixture is checked by tolerance, Delaunay and determinism only.
    Kept: the wait is a scheduling rule, not an output property, and
    removing it is not a correctness mutant.

Off-switch digests recorded from `c074f900` (master's production code):
features start, quality on, feet off `0x1253792f2761a97d`; `refine_points`
with a strip, feet off `0xa460db59ad797984`; `refine_strip`, feet off
`0x93e787a9ab2b77fd`. Never re-recorded. The design's "T12 fixtures and a
domain start" half of the off switch is 20b's F5
(`tests/cpp/property/prop_refinement_constraint_feet.cpp`), which now also
covers the quality start, since `refine`'s switch is the pass's (R2.6).

**Open items.**

- (a) **`prop_refinement_constraint_feet.cpp`, before green, by `@tester`,
  its own commit.** No shims: R1 says 20b's helpers "become" the new header.
  - `detail::foot_fits` moves, same signature, to
    `terrain::mesh::detail::foot_fits` in `constraint_foot.hpp` (the helper
    calls it for `NotCounterClockwise`). The "R2 step 4" case changes only
    its qualifier, and keeps its direct kills of the two one-sided mutants,
    which CF1's both-sides fixture cannot tell apart.
  - `detail::foot_of` is removed. Its "R2 step 2" case is rewritten against
    `constraint_foot(m, 0, MeshVertex{col, row}, δ, LatticeFrame{10, 5})`, with
    δ = `foot_epsilon(dem, node, 1.0)` (the cap, 2.5 m, on the flat DEM): the
    far node is a `Hit` on edge 0 at (4, 0.8), the two near-end nodes are
    `NearEnd`, and the control at col 3 is a `Hit`.
  - `foot_epsilon` stays in `refine.hpp`; its cases are unchanged.
  - `classify`: the vertex identity becomes start + `quality_inserted` +
    `quality_feet` + `inserted`, and the feet it finds must number `feet +
    quality_feet`. A quality-start foot passes its other checks as it
    stands: it lies on an input segment, it is the projection of a valid DEM
    node within the cap (δ_q is the cap), and its z is bilinear. The "no two
    feet share a source node" check stays over all feet; if a fixture breaks
    it, `@tester` brings the case to `@architect` and does not weaken it.
  These edits leave the file red by compile until green, like the rest of
  the red step.
- (b) **`test_cli_constraint_feet.py::test_the_default_leaves_the_10m_quarter_circle_alone`,
  after green, by `@tester`, its own commit, only if it fails.** Its claim
  is 20b's: refinement's own rule never fires at 10 m. If the quality start
  or the final check now foots there, the amendment keeps `out.feet == 0`,
  adds `out.feet_refused == 0`, moves the recorded `("quarter_circle", "10")`
  digest into `test_no_constraint_feet_matches_increment_20`'s cases (where
  `--no-constraint-feet` must still match it), and renames the test to say
  what it now checks (refinement puts no feet on the 10 m quarter circle).
  The recorded digest is not re-recorded. Any other existing pin that green
  breaks is handled the same way: listed with its reason, amended by
  `@tester` in its own commit, and a pin that claims "equal to master" is
  moved behind `--no-constraint-feet`, never re-recorded.
- (c) **The TSan job, by `@developer`, in the green commit.**
  `test_constraint_foot_quality`, `test_constraint_foot_refine` and
  `test_constraint_foot_final` join both lists in
  `.github/workflows/main.yaml`'s `tsan` job (the build targets and the run
  list, which must match). `test_constraint_foot` does not; it starts no
  threads. Confirmed.

#### Rulings on 20c-1's green step (`@architect`, 2026-10-07, on `8f275220`)

`@developer`'s green commit (147 net LOC) left four things failing and made
five choices the design did not. Each was checked against the code first.

1. **F2's `feet_refused == 1` becomes 2: correct under R3, `@tester`
   re-pins.** Checked with a scratch build of F2's fixture: the two refused
   nodes are the needle (row 16, col 8) and its neighbour above, (row 15,
   col 8), 0.0158 cells right of the side. That node's foot, (col 7.984,
   row 15.0008), lies in a cell whose column-7 corners are the fixture's
   NaN nodes, so `vertex_z` refuses it, and R3 and R5 say a refused foot is
   counted and the node inserted. 20b never asked: the node's own triangle
   has no constrained edge, only its neighbour does, which is exactly what
   R3 adds. Without the NoData column the same run foots that node (4 feet,
   0 refused). The pin becomes `== 2`, with that reason in the comment, and
   the case also asserts that both nodes are vertices.
2. **CF4's 2 Delaunay violations on macOS: an oracle defect, `@tester`
   fixes the oracle; nothing in production.** Checked: the case fails in a
   default Release build and passes with `-ffp-contract=off`; the library
   (`libterrain_predicates.a`, which holds the exact predicate) built either
   way makes no difference, only the contraction in the test's own
   translation unit, which compiles the header-only producer and its
   `x_min + col * dx` output. The violating quad is vertices 256 (a foot),
   339 (a strip point), 290 and 389 (nodes): two points on the same tilted
   side beside two nodes, cocircular to within rounding (exact incircle +1.77e-8 on terms around
   60). The oracle
   rebuilds (col, row) from world points at y_max = 7e6 m, which carries
   about 1e-9 m of rounding, and asks the exact incircle of that rounded
   copy: the exact predicate on rounded data, not the producer's relation.
   **A small-origin fixture does not fix it**: with y_max = 160 m, where the
   rebuild is good to about 1e-14 m, the same case still fails, with 2
   violations under contraction and 4 without (determinants around 1e-13),
   because the producer routinely makes such near-cocircular quads along a
   footed line. The fix, in `tests/cpp/support/constraint_foot_oracles.hpp`'s
   `delaunay_violations`:
   - an edge is a violation when the exact incircle says Inside **and** the
     apex lies inside the circumcircle by more than
     η = 1e-7 · min(dx, dy) in the producer frame (R − |apex − centre|,
     centre and radius in double). Scale: 5e-7 m on the 10 m × 5 m fixtures,
     about 500 times the rebuild's rounding at y_max = 7e6 m; checked there;
     it holds while half an ulp of the largest world coordinate stays under
     η / 100 (world coordinates below about 4e7 m at 5 m cells). A real
     missed flip is off by a fraction of a cell, not by 5e-7 m.
   - the case reports the depth of every quad it excuses, and `@tester`
     shows the oracle can still fail: with `legalise_around` skipped after a
     final-check foot split (a planted mutant in `refine_points.hpp`,
     restored and `touch`ed afterwards), CF4 must go red with depths far
     above η. That is the evidence the margin does not hide a defect.
   **CI**: `macos-latest` is Apple clang on arm64, which contracts by
   default, so this case would fail there as it does locally. The Linux
   jobs use GCC on x86-64; GCC's default for C++ is to contract, but the
   build passes no `-march`, so the baseline instruction set has no fused
   multiply-add to contract into and the Linux build behaves as
   `-ffp-contract=off` does. That pass is luck, not evidence (the small-origin
   run fails without contraction too). 20b's `delaunay_oracle` in
   `tests/cpp/property/prop_refinement_constraint_feet.cpp` rebuilds the
   same way and is green on both today; it is not changed in 20c-1, and gets
   the same margin if it ever goes red on a platform.
3. **The two stub tests: `@tester` amends them.** Pin 10 puts
   `constraint_feet` on `_core.refine_points` and `_core.refine_strip`, and
   both bindings (`bindings/core.cpp`) and `_core.pyi` take it last, after
   `frozen_mask` (23b's "goes last" rule).
   `tests/python/test_core_edge_strip.py::TestStubs::test_refine_strips_parameters`
   expects `["tolerance", "threads", "frozen_mask", "constraint_feet"]`, and
   `tests/python/test_core_refine_points.py::TestStubs::test_refine_points_takes_tolerance_threads_strip_and_frozen_mask_by_keyword`
   expects `["tolerance", "threads", "strip", "frozen_mask",
   "constraint_feet"]`; each comment says 20c-1 (pin 10) added it.
4. **Stale citations: fixed in this commit.** `refine.hpp` went from 436 to
   405 lines. The six broken citations (`25-plain-output.md` lines 115 and
   120, `27-node-sampling.md` lines 124, 150, 171 and 263) and the at-risk
   ones whose line now says something else (`27-node-sampling.md` 115,
   `15f-edge-strip.md` 64, 1554 and 1904, `23-basin-scale.md` 3267 to 3271)
   are pinned to `ed125121`, where each line was read and says what the text
   claims; `15f-edge-strip.md` 1904 to `17c2d14`, the commit its review
   names, since at `ed125121` line 101 is another field.
   `python3 tools/check_citations.py` exits 0.
5. **The five choices beyond the design:**
   - **The validity callable: changed, to the design.** The green commit
     keeps `improve`'s constraint `std::predicate<const LatticeVertex&>`
     and judges a foot only if the callable also takes a `MeshVertex`; a
     `LatticeVertex`-only callable silently refuses every foot, so feet on
     would quietly do nothing for such a caller. R2.3 says the callable
     takes a `MeshVertex`. First `@tester`, in its own commit: the seven
     callables in `tests/cpp/unit/test_mesh_quality_void.cpp` take
     `const MeshVertex&` (the `asked` vector too; `MeshVertex ==
     LatticeVertex` exists), which compiles and passes against the green
     code as it is. Then `@developer`: the constraint becomes
     `std::predicate<const MeshVertex&>`, the `if constexpr` goes, the foot
     is judged by `valid(s.at)`, and the header comment's sentence about a
     `LatticeVertex`-only callable goes.
   - **`constraint_foot.hpp` includes `lawson.hpp` for `LatticeFrame`:
     kept.** `lawson.hpp` depends on core, `lattice_mesh.hpp` and the
     predicates only, the set R1 allows; R1 now names it.
   - **A footed DEM node's foot is not in `nodes_inserted`; its fallback
     insertion is: kept.** A foot is not a DEM node (pin 5's reason), and
     the field stays a subset of `inserted`.
   - **`feet_fallback` counts any later insertion of a footed position:
     kept.** That is R4.5's definition; a footed point goes in as itself at
     most once (it is then a vertex), so `feet_fallback <= feet` holds.
   - **`foot_epsilon` for every split with feet on: kept for now, measured
     by `@perf`.** 20b computed it only when the triangle had a constrained
     edge; it now reads up to 16 DEM values per split in the serial phase.
     `@perf`'s acceptance reports the split phase's seconds, feet on against
     `--no-constraint-feet`, on both catchments; if the feet-on figure is
     more than 2 % higher, `@developer` computes ε only when the triangle or
     a neighbour across an unconstrained edge has a constrained, non-frozen
     edge (no change to the output). **Done**: `@perf` measured it over
     (`b9e2d482`), and `@developer`'s `69f37d1c` does that and one step
     more (R3, "As built"): no search where `foot_reachable` is false, the
     search to the cap first, ε only when an edge is within the cap. The
     limit itself is judged in the next section.

**Next, in order.** `@tester`, one commit each: F2's re-pin (1); the CF4
oracle margin with its mutant demonstration (2); the two stub tests (3);
`test_mesh_quality_void.cpp`'s callables (5). Then `@developer`: the
`MeshVertex` constraint on `improve` (5). Then the full C++ and Python runs,
`@reviewer`, and `@perf`. (All done: `b09b84f6` to `34b546c4`, then
`@perf`, the fix and the re-time, as the Status says.)

#### The split-phase limit, judged after the fix (`@architect`, 2026-10-07, on `69f37d1c`)

**What `@perf` measured at `69f37d1c`**, median of 6 interleaved runs per
build on AC power, threads 10
(`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/README.md@c9f62f7c:203-241`):

| run | split phase, feet on vs off | per split, on vs off | split phase, on vs master | start quality, on vs off | start quality, on vs master | refine, on vs master | whole run, on vs master |
|---|---|---|---|---|---|---|---|
| 1 m tile, tolerance 1 (no foot placed; one mesh for all three builds) | +0.86 % | +0.86 % | −0.37 % | a 0.6 ms phase | a 0.6 ms phase | −0.64 % | −0.85 % |
| Lagan | +7.1 % | +8.3 % | −3.3 % | +6.9 % | +7.7 % (+25 ms) | +4.1 % (+23 ms) | −0.26 % |
| Numedalslagen | +6.7 % | +8.9 % | −7.7 % | +5.7 % | +6.5 % (+32 ms) | +1.7 % (+19 ms) | +0.18 % |

Master's own spread in the split phase is 2.1 to 2.5 % (same file, lines
235-241). The 1 m benchmark (`tools/bench.py`, both domains, 1 to 20
threads) is 0.3 to 2.0 % faster than master in every cell, the
1-to-20-thread speed-up unchanged (2.56 × against 2.55 × on the tile, 2.54
× against 2.51 × on the quarter; lines 261-321).

**Did the 2 % limit measure the right thing? On the tile yes, on the
catchments no.** It was written to catch work done for nothing: ε's DEM
reads on every split, where no constraint is near. Feet on against feet
off measures that only where both runs do the same work, which is the 1 m
tile: no foot is placed and all three builds give one mesh. There it is
+0.86 %, inside master's own spread, so that work is gone. On the
catchments the two runs build different meshes and do different work.
Feet on places 720 and 4 629 feet in this phase, each a split of a
constrained edge with its own legalisation, which feet off never does; it
also does fewer splits (46 144 against 46 656 on Lagan, 183 368 against
187 233 on Numedalslagen; lines 223-233). So the catchment figure measures
the feature, not overhead, and neither seconds nor seconds per split
compares like with like. Why a split costs more with feet on is **not
known**: nothing has been profiled; it may be the remaining foot search
and ε, or the feet themselves. The start-quality phase, which the limit
never covered, rises about 6 % for the same reason (5 930 and 4 233
feet).

**The limit, revised** (replaces ruling 5's 2 % for 20c-1, and holds for
20c-2, which touches the same paths, except that 20c-2's split phase in
point 2 is judged per split: ruling 8, "The split-phase limit for 20c-2,
judged per split"):

1. **No cost where no foot is placed.** On the 1 m tile, where the meshes
   are byte-identical with feet on and off, the split phase with feet on
   is within 2 % of feet off and of master. 20c-1: +0.86 % and −0.37 %:
   **met**.
2. **No slower than what Ola runs today.** On both catchments, against
   master's default (20b's refinement feet on): the split phase not more
   than 2 % slower, and the whole run within 2 %. 20c-1: split phase
   −3.3 % and −7.7 %, whole run −0.26 % and +0.18 %: **met**.
3. **Reported, not gated**: start quality and refine against master and
   against feet off, as in the table. They pay for 58 to 74 % fewer
   slivers, the product Ola asked for, and no cheaper form has been
   measured.

Scale of the 2 %: relative, a median of at least 6 interleaved runs, since
a single phase's spread on master is itself about 2 %. Checked on the 1 m
tile (464 290 triangles), Lagan (867 612) and Numedalslagen (1 290 807).
At São Francisco's size the serial phases matter more than anywhere; point
2 keeps the refinement split phase no slower than master's, and the
start-quality phase is the one to watch (+6.5 to +7.7 % on master, serial).

**Ola's question** (question 7, now ruling 7: accepted): accept the remaining +7 to 9 %
against feet off, or move the foot search into the parallel scan.
Recommended: accept. The most a move could save is the whole feet-on
excess of the split phase, 2.4 ms on Lagan and 9.7 ms on Numedalslagen,
under 0.1 % of either run, and a scan-time answer would still be checked
again in the serial phase, since an earlier split in the same round can
change a neighbour the search reads.

#### Mutation round for 20c-1 (`@tester`, 2026-10-07, test commit `4a64c1ef`, on code `0346ac63`)

Answers code review round 1's B1 (no kill record) and B2 (the CF4 margin
demonstration's depths). Run on a scratch copy (`tools/scratch_copy.py`),
macOS arm64 Release, with the targets `test_constraint_foot`,
`test_constraint_foot_quality`, `prop_constraint_foot_refine`,
`prop_constraint_foot_final` and 20b's `prop_refinement_constraint_feet`.
Each fault was planted in the production code, the targets rebuilt and
run, and the fault removed. Line numbers are in the test files at
`4a64c1ef` (`test_constraint_foot.cpp` is the CF1 file; `_quality` CF2,
`_refine` CF3, `_final` CF4; "20b" is `prop_refinement_constraint_feet.cpp`).

**The design's seven faults** ("Mutants to kill", 20c-1 above):

| # | fault planted | result | caught by |
|---|---|---|---|
| 1 | the end check removed (`foot_on`: a foot within δ of an end of its edge is still a hit) | killed | CF1 "a foot within delta of either end of its edge is NearEnd", lines 218 and 228 (a hit where NearEnd is expected); 20b "R2 step 2", line 708. The variant that checks only the first end: lines 219 and 708 |
| 2 | the neighbour search crosses a constrained edge | killed | CF1 "no search across a constrained edge", line 255 (constrained) and 262 (frozen) |
| 3 | the edge's midpoint instead of the foot | killed | CF1's position check in 8 of 10 cases (line 121, off by up to 9.4 cells); CF2 line 331; CF3 line 174; CF4 lines 230, 250, 291, 419; 20b line 229 |
| 4 | footed once removed (the final check's `&& !was_footed` dropped, so a point already footed may be footed again) | **survived**; killed by `2bf2acf2` | no fixture reached it; ruled below (b); the kill is in "The survivors, killed" |
| 5 | the fallback removed (a footed point that comes up again is dropped, not inserted as itself) | killed | CF4 "at tolerance 0 every footed point goes in after its foot", line 382 (39, 40 and 40 points with error above 0 m for seeds 1 to 3); also lines 227 (`feet_fallback` 0, not 1) and 253 |
| 6 | the foot's height taken from the point (its own z, not the foot's) | killed | CF4 line 252 (z 80 where 73.5 is expected), line 227, and the projected path's height case, lines 272 and 289 (z off by 0.245 m) |
| 7 | refinement: the holding triangle dropped from the active set when the foot goes on a neighbour's edge | killed | CF3 "the foot goes on the neighbour's edge, and the node still goes in from the rescanned holding triangle", line 177 |

**Further faults in `constraint_foot.hpp`:**

| fault planted | result | caught by |
|---|---|---|
| row and column swapped in the world distance | killed | CF1 "distance and foot are measured in the world frame", lines 293 and 121; 20b line 229 |
| the frozen-edge test dropped | killed | CF1 "a frozen edge is never a foot", lines 272, 277 and 262 |
| the nearest edge instead of the first in search order | killed | CF1 "t's own constrained edges come first, in edge order", lines 118 and 119 |
| "closer than δ" made non-strict (`>= delta` to `> delta`) | **survived**; killed by `cec2dda4` | no case puts a point exactly δ from an edge (CF1 tests 0.4 and 0.6 cells against δ = 0.5); the suite's header comment says it kills this fault, which is false; ruled below (a) |

**Faults against the speed fix `69f37d1c`** (`foot_reachable`, and the search to the cap before ε):

| fault planted | result | caught by |
|---|---|---|
| `foot_reachable` ignores the neighbours (false where only a neighbour has a live edge) | killed | CF3 line 174; 20b F2 line 498 (`feet_refused` 1, not 2) |
| `foot_reachable` ignores the triangle's own edges | killed, by 20b only; now also by `7bc41ebf` | 20b F3 line 514 (no feet at tolerance 0 on slope 0.1287); no 20c test caught it; ruled below (d) |
| `foot_reachable` always true | **survived**; killed by `7bc41ebf` | the output cannot change, only work is added; no test calls `foot_reachable` directly; ruled below (d) |
| ε never computed: the search to the cap alone decides | killed by `4a64c1ef` | before it, only 20b F3's Delaunay oracle reacted (line 365, two violations 1.43e-11 m and 9.77e-11 m inside a circle of radius 6.52 m: rounding, so counted as a survivor). `4a64c1ef` adds CF3 "a node within the cap but beyond eps of the neighbour's edge goes in as itself, not as a foot" (tolerance 0.35 m, slope √2 m/m, so ε = 0.247 cells; the node 0.3 cells from the edge; cap 0.5 cells; the case asserts both facts first). It fails at line 224 on the fault and passes on `0346ac63` (7 CF3 cases, 131 assertions) |

**Also tried**: the final check's holding triangle dropped from the active
set when its foot goes on a neighbour's edge (R4.2): **survived**; no CF4
fixture puts a final-check foot on a neighbour's edge. Ruled below (c);
killed by `8e2fb5fb`.

**The survivors, killed** (`@tester`, 2026-10-07, commits `cec2dda4`,
`2bf2acf2`, `8e2fb5fb`, `7bc41ebf` on `aee60fd2`; scratch copy, macOS
arm64 Release; no production change: `git diff --stat 0346ac63 7bc41ebf
-- include src` is empty). Each new case fails on its planted fault and
passes on the code. Line numbers here are at `7bc41ebf`; the tables above
keep theirs at `4a64c1ef` (CF1 lines after 155 and CF4 lines after 420
have since moved).

| fault planted | killed by | fails at |
|---|---|---|
| "closer than δ" made non-strict (`include/terrain/mesh/constraint_foot.hpp@0346ac63:64`) | `cec2dda4`, CF1 "a point exactly delta from a constrained edge of its triangle is None": `V(5, 0.5)`, after asserting the recomputed distance is exactly 0.5 | `tests/cpp/unit/test_constraint_foot.cpp@7bc41ebf:170` (a Hit, not None), the only failing case; the header claim corrected at `tests/cpp/unit/test_constraint_foot.cpp@7bc41ebf:44` |
| footed once removed (`include/terrain/refinement/refine_points.hpp@0346ac63:386`) | `2bf2acf2`, CF4's two-parallel-lines case (ruling (b)'s fixture) | `tests/cpp/property/prop_constraint_foot_final.cpp@7bc41ebf:514` (a second foot at (5, 0) on A1-B1); CF4's header corrected at `tests/cpp/property/prop_constraint_foot_final.cpp@7bc41ebf:34` (the wedge's bound did not kill this fault) |
| R4.2: the holding triangle dropped from the active set (`include/terrain/refinement/refine_points.hpp@0346ac63:428-429` removed) | `8e2fb5fb`, CF4 "the foot goes on the neighbour's edge, and the point still goes in from the holding triangle" (CF3's geometry; a companion case asserts the holding triangle is unchanged and untouched by the flips) | `tests/cpp/property/prop_constraint_foot_final.cpp@7bc41ebf:602` (N never goes in) |
| `foot_reachable` ignores the triangle's own edges | `7bc41ebf`, CF1's direct `foot_reachable` test | `tests/cpp/unit/test_constraint_foot.cpp@7bc41ebf:309` |
| `foot_reachable` ignores the neighbours | `7bc41ebf` | `tests/cpp/unit/test_constraint_foot.cpp@7bc41ebf:315` |
| `foot_reachable` always true | `7bc41ebf` | `tests/cpp/unit/test_constraint_foot.cpp@7bc41ebf:322` and `:330` |

So every one of the 16 planted faults is now killed by a 20c test. Ruling
(b)'s fixture reaches the guard in the code: CF4 "the two-lines fixture is
what it claims" replays `refine_points`' own steps
(`tests/cpp/property/prop_constraint_foot_final.cpp@7bc41ebf:463-500`):
`legalise_all` makes no flip on the rectangle, the first Hit is on A2-B2 at
x = 5, the split and `legalise_around` flip A1-B2, and in the triangle then
holding the point `constraint_foot` gives a second Hit on A1-B1 at (5, 0).
A first pass over the three `foot_reachable` faults reused a stale binary
(the restore and the planted edit fell in the same second, and filtered
build output hid it); all three were rerun with the object file deleted
before each build, and the lines above are from the rerun. Full runs on
`7bc41ebf`: ctest 886 of 886; `pytest --no-cov` 5404 passed, 17 skipped,
0 failed.

**The CF4 margin demonstration** (green-step ruling 2; code review round
1's B2). Fault: `legalise_around` skipped after a final-check foot split.
Every CF4 case that places a foot then fails in its Delaunay check
(`prop_constraint_foot_final.cpp` line 181). How deep the offending apex
lies inside the circumcircle:

| case | depths | η (the check's allowance) |
|---|---|---|
| strip and no-strip cases | 4.96 m, 3.87 m | 1e-7 m |
| projected path's height case | 6.93 m, 1.29 m | 5e-7 m |
| tolerance 0, seed 1 | 10.33, 0.411, 0.234, 0.0196 m | 1e-7 m |
| tolerance 0, seed 2 | 10.41, 0.288 m | 1e-7 m |
| tolerance 0, seed 3 | 10.33, 0.425 m | 1e-7 m |

η is 1e-7 × the cell side (1 m in these fixtures, 5 m on the projected
path). The shallowest depth, 0.0196 m, is about 2e5 times its η. On
`0346ac63`, without the fault, the
two quads CF4 excuses are 3.85e-10 m and 4.11e-10 m deep, about 1 300
times under η. So η sits between rounding and the smallest real failure
with more than three orders of magnitude on each side.

#### Rulings on the survivors (`@architect`, 2026-10-07, on `4a64c1ef`)

**(a) "Closer than δ" made non-strict.** `@tester` adds the case
proposed: a CF1 fixture point at `V(5, 0.5)`, exactly 0.5 cells (δ) from
the edge in doubles, expecting `None`; the case first asserts that the
distance it recomputes is exactly 0.5. And corrects the header comment
(`tests/cpp/unit/test_constraint_foot.cpp@4a64c1ef:44`), which today
claims the 0.4 and 0.6 cases kill this fault. R1's "closer than δ is
strict" stands; no production change. **Done**, `cec2dda4` ("The
survivors, killed", above).

**(b) Footed once removed. Keep the guard: it is reachable, and the
fixtures so far could not reach it.** `@tester`'s explanation (the cut at
the foot makes the repeated search NearEnd) is right for the line that
was footed, and only for that line:

- *The footed line never gives a second hit.* After the foot F goes in,
  the line is cut at F into constrained edges that end at F, and later
  splits only add vertices on it. F is p's orthogonal projection onto the
  line, and F is a mesh vertex, so it lies in the open interior of no
  edge. For any edge on that line, p's nearest point on the closed edge
  is therefore an end of it (the clamped s is 0 or 1, up to rounding of
  F's stored position, which is far below δ). `foot_on` then returns
  NearEnd (s·len < δ holds at an end) or nothing (too far); never a hit
  (`include/terrain/mesh/constraint_foot.hpp@0346ac63:63-67`).
- *A second line can.* `constraint_foot` returns the first verdict in its
  search order. If, after the foot, p's triangle has an edge of another
  constrained line within δ, and no edge of the footed line comes before
  it, that edge gives a hit. The wedge runs met the footed line's edge
  first every time: two lines meeting at a vertex keep the cut sub-edge
  on p's new triangle.
- *A fixture that does it* (computed with exact predicates on these
  coordinates, not yet run through the code): two parallel constrained
  edges A1 = (0, 0) to B1 = (10, 0) and A2 = (0, 0.6) to B2 = (10, 0.6),
  in a frame with dx = dy = 1 (so δ_p = 0.5 cells), the strip between
  them triangulated (A1, B1, B2), (A1, B2, A2), and a stored point
  p = (5, 0.35) with error above the tolerance (tolerance 0). p lies in
  (A1, B2, A2), whose only edge within δ is A2-B2 (0.25 cells): a hit at
  F = (5, 0.6). After the split, B1 lies inside the circumcircle of
  (A1, B2, F), so Lawson flips A1-B2, and p is now in (A1, B1, F), whose
  edge A1-B1 is 0.35 cells away with its foot (5, 0) 5 cells from either
  end: a second hit. With the guard, p goes in as itself (feet 1,
  fallback 1); without it, a second foot goes in at (5, 0), and p owns
  two feet, which CF4's oracle "no stored point owns two feet" refuses
  (the bound the design named for this fault). CF4's wedge was meant to
  reach this ("a point footed on one line can find the other next
  round"); it does not, for the reason in the second point.
- **Who does what**: `@tester` adds this case to CF4 (or one like it, if
  the hand-built mesh needs more around it) and shows it fails with the
  guard removed. If the code does not reach the guard on it, `@tester`
  says so in the handback, and the guard comes back to me before code
  review round 2. R4's text ("recorded as footed … inserted later as
  itself") already states the guard; it is unchanged, and there is no
  work for `@developer`.
- **Done**, `2bf2acf2`: the code reaches the guard on this fixture
  (replayed, "The survivors, killed", above), so the guard stays and
  nothing came back to me. `@tester` also constrained the strip's two
  short sides, as every CF4 fixture does; they lie 5 cells from p and
  change nothing here.

**(c) The final check's holding triangle dropped (R4.2).** This one is a
correctness fault, not dead code: if the holding triangle t is neither
touched nor kept active when the foot goes on a neighbour's edge, t is
not rescanned and its point is never inserted, which breaks the tolerance
guarantee R4.2 states. `@tester` adds the CF4 analogue of CF3's
neighbour case: a stored point at tolerance 0 whose holding triangle has
no constrained edge but whose neighbour, across an unconstrained edge,
has one within δ_p, laid out so that no flip after the foot touches the
holding triangle (as CF3's case is; a flip would mark it touched and hide
the fault); asserted: the foot is on the neighbour's edge, the
point then goes in as itself, and the tolerance oracle over all stored
points passes. It must fail with
`include/terrain/refinement/refine_points.hpp@0346ac63:428-429` removed.
**Done**, `8e2fb5fb`.

**(d) `foot_reachable`.** Yes, a direct CF1 unit test. It is the contract
the speed fix rests on (where it is false, `constraint_foot` finds
nothing), and today one of its faults is caught only by 20b's F3. On
CF1's fixtures, four cases: true where only the triangle's own edge is
constrained and live; true where only a neighbour's is; false where the
only constrained edge lies beyond a constrained edge of the triangle;
false where the only one is frozen. Each false case also asserts that
`constraint_foot` returns `None` for a point next to that edge. These
kill "ignores own edges" and "ignores neighbours" in a 20c target, and
"always true" too, although that fault changes only the work done
(`@perf`'s split-phase figures cover the work). **Done**, `7bc41ebf`. The
"beyond a constrained edge" case makes t's own constrained edge frozen,
with the neighbour's edge live: the only reading in which that case is
false, since a live constrained edge on t would make it true.

**Next, in order.** `@tester`, in one test commit or one per item: (a),
(b), (c), (d), each shown to fail on its fault and pass on the code; the
handback adds one line per item to the tables above (the commit, and the
line that kills). No production change; no `@developer` step unless (b)'s
fixture does not reach the guard. Then the full C++ and Python runs, and
code review round 2 (`@reviewer`). The kills and the full runs are done
(above); next is code review round 2.

### 20c-2

**Invariant-critical (mutation testing required): `test_quality_gain`**
(new, C++, `tests/cpp/unit/test_quality_gain.cpp`, its own target in
`tests/cpp/CMakeLists.txt`), holding T-P1, T-P2 and LS1:

- **T-P1 (invariant)**: for every accepted candidate in a fixture run, the
  predicted new triangles equal the triangles `legalise_around` wrote, as
  sets of vertex triples.
- **T-P2**: no accepted candidate lowers its cavity's worst angle (capped
  at θ) by more than rounding; a fixture where the only candidate would
  lower it is refused and counted in `skipped_no_gain`.
- **LS1 (R8)**: a bad triangle whose node lies beyond a long constrained
  edge: the foot of the node goes on that edge, not the midpoint; none when
  the foot is within δ_q of an end, on an outline edge with nothing beyond,
  or on a frozen edge; the constraint edges as a set of lines are
  unchanged, bits and masks on both halves; with gain −1, no split; with
  the feet off, no split.
- Mutants to kill: the cavity crossing a constrained edge (T-P1); an
  inexact incircle, a plain double determinant instead of the kernel's
  (T-P1, on a near-cocircular fixture); the foot's second seed missing
  (T-P1); the acceptance test reversed (T-P2); the slack dropped and the
  test made strict, `new > old + P` (T-P2, the equality fixture); the split
  at the midpoint (LS1); the end check removed (LS1); R8 left on at gain −1
  (LS1, T-P3); R8 left on with the feet off (LS1). (The slack dropped with
  `≥` kept is killed only where the platform's rounding puts `new` below
  `old` on CF2's own edge with the feet on; a survivor there is recorded
  with that reason, not chased.)

Not mutation-critical:

- **T-P3**: gain −1 is bit-identical to 20c-1; determinism as T6.
- CLI: `--start-quality-gain` (refusals: non-finite, above 10), the new
  `--stats` rows and the record.

#### Pins ruled for 20c-2's red step (`@architect`, 2026-10-07, on `100f5a09`)

`@tester`'s red commit `100f5a09` (on master `41bda81a`) fixed eleven
details the design left open and made three observations. Each is ruled
here; "kept" means `@developer` builds to it as written. `@tester` built a
throwaway prototype of R7 and R8 beyond its brief; its numbers (pin 10,
observations (a) to (c), the box's 20.87° to 20.66°) are a prototype's, not
measurements of built code, and nothing below rests on them alone.

1. **Names and defaults.** `QualityOptions::min_gain_deg` and
   `RefineOptions::min_gain_deg`, negative = off, default negative, set by
   member; binding keyword `min_gain_deg`, default `-1.0`; outcome fields
   `QualityOutcome::skipped_no_gain` and `line_splits`,
   `RefineOutcome::quality_no_gain` and `quality_line_splits` (Python
   `RefineOutcome.quality_no_gain`, `.quality_line_splits`). Kept: off by
   default keeps every existing C++ and binding caller on 20c-1's output;
   only the CLI turns it on (R7).
2. **The cavity helper.** `terrain::mesh::detail::quality_cavity<K>(m, t,
   on, p, f)` in `quality.hpp`, returning `QualityCavity{removed,
   created}`: `removed` the cavity's slots, `created` the new triangles as
   corner triples, counter-clockwise on (col, −row); `on` 0..2 is t's edge
   (a split, the triangle across seeding too, constrained or not), 3 is
   strictly inside; order free. Kept: it is R7's prediction as a function
   of the mesh alone, which is what lets T-P1 compare it with
   `legalise_around`'s result without the whole pass. `improve` asks it
   for the point it is about to insert (node or R2's foot), after every
   existing skip, and reads `created` only for angles.
3. **Angles and rounding.** In the frame `improve` is given, (col · dx,
   −row · dy), each capped at θ. **Changed:** R7 now accepts `new ≥ old +
   P − s` with s = 10⁻⁶° (R7, "Why the slack"; observation (a) below), and
   T-P2's "by more than rounding" is the same s. Tests that change, by
   `@tester`, its own commit before green: `kRounding` in
   `tests/cpp/unit/test_quality_gain.cpp` becomes 1e-6 with a comment
   citing R7's slack; a new case, CF2's own edge (`own_edge()`) with the
   feet on at gain 0: the first insertion is the foot (10, 0.99), within
   1e-12 cells, it is accepted (`feet ≥ 1`, and `check_gain` holds on
   every step); its premise asserts the foot's two new triangles hold
   18.43° at B and at C (the foot is 10 cells from both B and C). The
   equality fixture's comment names its mutant as "the slack dropped and
   the test made strict". The developer names the slack as a constant in
   `quality.hpp`; the test keeps its own literal.
4. **R8 and `constraint_feet`.** Pinned by the tester as independent of
   the feet. **Changed:** R8 runs when R7 is on **and** `constraint_feet`
   is (R8, "Only with the soft criterion, and only with the feet").
   `--no-constraint-feet` has kept every point off the lines on every
   path since 20c-1 (R2.6, R5); R8 is a fourth such path. Tests that
   change, same commit: in `test_quality_gain.cpp`, every LS1 case that
   expects no split (near an end, outline, frozen, refused foot, gain −1)
   passes `feet = true`, so that what refuses the split is the rule under
   test and not the switch; "splits the line at the node's foot" runs with
   the feet on only, and a new case "LS1: with the feet off, no split" on
   `line_beyond()` at gain 0 asserts `line_splits == 0`,
   `skipped_blocked ≥ 1` and e whole (kills "R8 left on with the feet
   off"); "across the line fixture every insertion ... is the predicted
   cavity" keeps both feet values but asserts `line_splits > 0` only with
   the feet on and `== 0` with them off. The pinned-list comment at the
   head of the file and the Python module docstring ("R8's split is on
   exactly when the gain test is") say so. No other file changes: the
   refine suite's scenes run with the feet on.
5. **An R8 split is not a skip; a declined split counts once in
   `skipped_blocked`.** Kept: that is today's count for a blocked walk, and
   a split is an insertion.
6. **Counts.** `quality_skipped` holds `quality_no_gain` (R7: the total is
   kept); a line split counts in none of `quality_inserted`,
   `quality_feet`, `inserted`; vertices = start + `quality_inserted` +
   `quality_feet` + `quality_line_splits` + `inserted` (and, in `improve`,
   start + `inserted` + `feet` + `line_splits`). Kept: each counter names
   one kind of vertex, as 20c-1's pin 5 has it.
7. **The recorder's hook.** The tests rely on `improve` asking the validity
   callable at least once in every loop turn that inserts, before it
   inserts. Kept as a contract of `improve`, and written here so a later
   change does not break it silently: today it asks about every node
   before the walk; R2 asks about the foot and R8 about its foot. An
   insertion with no call before it would make the recorder see two
   insertions as one, which `recorded`'s `REQUIRE(ok)` fails loudly.
8. **`--start-quality-gain`.** Finite values up to and including 10
   accepted, a negative one is the hard rule; refused: `nan`, `inf`,
   `-inf`, `10.5`, `11`; refused without `--tolerance` and without `--dem`,
   as `--start-min-angle` is. Kept: R7's "refused: non-finite, above 10",
   and the same dependencies as the setting it qualifies.
9. **Report rows.** Input row `start_quality_gain_deg` beside
   `start_min_angle_deg`, printed `0` by default; result rows
   `start_quality_points_without_gain` ("Tries that would not have
   improved the angles") and `start_quality_lines_split` ("Land-cover and
   outline lines split to improve angles"), integers in the record, not in
   the mesh file. Kept. R7's `elevation_source` line was stale (increment
   25 replaced that sentence by rows); R7 and R8 now say the rows and
   where they go (R7's after `start_quality_points_skipped`, R8's after
   `start_quality_points_snapped_to_lines`). The input row's wording,
   which the tests do not pin: "Starting mesh: a point is added only if it
   raises the smallest angle around it by at least this many degrees
   (negative = always added)".
10. **CF2's own-edge start through `refine`, feet off**: `quality_no_gain
    == 1`, `quality_skipped == 1`, nothing inserted by the pass. Kept, and
    not on the prototype's word: the start is one triangle with one bad
    candidate, N = (10, 1); its three new triangles hold a sliver under 1°
    against the old 18.43° (the case's premise in `test_quality_gain.cpp`
    asserts that), far outside any slack; once refused, the queue is empty.
11. **T-P3's digests**, recorded from `41bda81a` (master with 20c-1, the
    red commit's parent, no 20c-2 production change) on the three scenes in
    `tests/cpp/support/quality_gain_fixtures.hpp`, and the CLI digests
    hashing integers only. Kept; never re-recorded. Integers only is right:
    a foot's last bit can differ between the fused multiply-add leg and
    the plain one.

**Observations.**

- (a) **A foot whose cavity's worst angle sits at its edge's end keeps it
  exactly, so a bare `≥` is decided by rounding.** Accepted as a defect of
  R7's text; fixed by the slack (pin 3). In CF2's own-edge case both the
  angle at B (the same two rays) and the angle at C (the foot is 10 cells
  from B and from C, so the triangle is isosceles) equal the old 18.43°.
- (b) **R8's end check is nearly unreachable on grid nodes.** Kept in
  R8 all the same: it is R8's termination argument (at most L / δ_q
  splits per segment), not an optimisation, and real starts are not grid
  nodes (CORINE and outline vertices lie anywhere), as the test's own
  `line_near_end()` fixture, with off-node vertices, reaches it.
- (c) **On outline-only starts the prototype refused nothing at gain 0.**
  Noted; no design change. The whole-run T-P1/T-P2 case requires
  refusals over the Q7 ring and the nine-node mesh together for each
  combination of frame, gain and feet; if at green a gain-0 combination
  finds none, that is a premise for `@tester` to repair (another fixture),
  not a reason to change the code.

**The box tests at `--start-quality-gain -1`.** Kept. The amended tests in
`tests/python/test_cli_start_quality.py` (the `--stats` count, NoData, and
the void tests and their `pass_node` fixture) are about increment 20's pass
when it inserts, and about what NoData does to an inserted node; at gain 0
the box's two candidates are refused, so they would fail or check
nothing. `-1` is the hard rule bit for bit (T-P3), so their claims stand
unchanged; the default on the box is covered by
`test_cli_start_quality_gain.py::test_the_default_changes_the_box`.

**Open items.**

- (a) **`@tester`, before green, its own commit:** pins 3 and 4's test
  changes above.
- (b) **Existing CLI pins at default flags, after green, by `@tester`, its
  own commit, only those that fail.** The CLI's default gain is now 0, so
  a digest recorded before 20c-2 at the default start-quality settings
  may change. Each such test gets `--start-quality-gain -1` (beside
  `--no-constraint-feet` where it has it), never a new digest, and the
  commit lists each with its reason, as 20c-1's open item (b). Likely, by
  reading: `test_cli_constraint_feet.py::test_no_constraint_feet_matches_increment_20`
  and `::test_the_default_leaves_the_10m_quarter_circle_alone`, and
  `test_refine_golden.py::test_the_cli_default_flags_digest_is_unchanged_by_node_sampling`.
- (c) **The TSan job, by `@developer`, in the green commit.**
  `test_quality_gain_refine` joins both lists in
  `.github/workflows/main.yaml`'s `tsan` job (the build targets and the
  run list, which must match). `test_quality_gain` does not; it starts no
  threads.

#### The split-phase limit for 20c-2, judged per split (`@architect`, 2026-10-07, on `5b0afe04`; Ola's ruling 8)

**What `@perf` measured** at `65ae0792` against 20c-1 (master `41bda81a`),
median of 6 interleaved runs after a warm-up, AC power, threads 10
(`docs/benchmarks/2026-10-07/20c-2/fix-65ae0792/README.md@5b0afe04:199-217`,
`:231-257`; the first acceptance, before the fix, is
`docs/benchmarks/2026-10-07/20c-2/README.md@fa136aa7`):

| catchment | split phase, seconds | split phase, per split | refinement's splits | whole run | start quality (reported) |
|---|---|---|---|---|---|
| Lagan | +16.29 % | −0.81 % | 54 097 against 46 144 (+17.2 %) | −0.04 % | +23.5 % (+0.082 s) |
| Numedalslagen | +4.00 % | −2.60 % | 195 799 against 183 368 (+6.8 %) | −0.08 % | +20.1 % (+0.108 s) |

On the 1 m tile, where every build gives one mesh, the split phase is
−4.24 % against feet off and −4.52 % against 20c-1 over 6 runs, −0.34 %
and +1.39 % over 20 (point 1: met). The 1 m benchmark is accepted, −0.8 to
+2.6 % per cell pooled over 3 pairs against a base spread of at most
2.7 %, the speed-up unchanged (same file, `:37-46`).

**Ruling.** For 20c-2, point 2's split phase is judged **per split**: the
median split-phase seconds over refinement's points added (`--stats`),
at most 2 % above 20c-1's. **Met: −0.81 % and −2.60 %.** The rest of the
limit stands as written and is met as measured: point 1 in seconds, point
2's whole run in seconds, point 3 reported.

**Why.** As for 20c-1 ("Did the 2 % limit measure the right thing?"),
the seconds compare two different meshes. R7 refuses 34 600 and 71 590
quality points, so the start adds 29 % and 34 % fewer, and refinement adds
17.2 % and 6.8 % more where the DEM needs them, for a mesh 7.5 % and
11.4 % smaller than master's. The limit is there to catch a dearer split;
a split at gain 0 is not dearer. Per split is not exact like for like
either (the splits land in other places), so the bound stays 2 % and the
whole run stays gated in seconds.

**Reported, not gated.** Per split against feet off: +9.63 % and +14.77 %
(20c-1's, accepted under ruling 7: +8.3 % and +8.9 %); the two runs also
differ in R8, and nothing has been profiled. Start quality, serial, is the
phase to watch at São Francisco's size: +0.082 s and +0.108 s here, 0.1 %
and 1.6 % of the whole runs.

Scale of the 2 %: relative, medians of at least 6 interleaved runs.
Checked on the 1 m tile (464 290 triangles), Lagan (799 372) and
Numedalslagen (1 141 207).

#### Rulings on 20c-2's green step and speed fix (`@architect`, 2026-10-07, on `d3c939d9`, `96119508`, `65ae0792`)

What the developer and `@perf` recorded beyond the design, each ruled;
all kept.

1. **One flip test, `detail::quad_flips`** (`lawson.hpp`): `must_flip`'s
   quad test moved out so the cavity asks the same predicate; R7's "the
   exact `incircle` in the frame that `legalise_around` uses" as one
   function. T-P1 then compares two users of one predicate and cannot see
   a fault inside it; `600b2acb`'s two named-cavity cases (a cocircular
   point on the integer path, and test_mesh_lawson's L4b quad on the
   frame-collinear branch) cover that. The fix forces it inline
   (`gnu::always_inline`), which GCC and Clang, the CI compilers, accept.
2. **`foot_on` takes the reach and the end distance apart**
   (`constraint_foot.hpp`). 20c-1's paths pass δ for both, so they are
   unchanged (gain −1 equals 20c-1 bit for bit); R8 passes no reach limit,
   since the walk has already found the edge, and δ_q (half a cell) for
   the ends, as R8 says.
3. **`GainJudge::capped` skips atan2 when every corner clears θ by
   1e-9 of |cross| + |dot|.** A relative margin, with no length scale: it
   is at most about 8 × 10⁻⁸° of angle, far above the rounding of the test, atan2
   and the degree conversion (about 10⁻¹⁴ relative). The argument holds
   for any θ up to 180°; nothing in the C++ or the binding bounds θ, and
   above 180° (already meaningless as an angle cap) the skip would wrongly
   return θ; the comment says so (`50544838`). A corner within the margin falls back to the exact call, so
   the capped value is unchanged; checked as byte-identical meshes at gain
   0 and −1 on Lagan, Numedalslagen, the tile and the quarter
   (`docs/benchmarks/2026-10-07/20c-2/fix-65ae0792/README.md@5b0afe04:139-163`).
4. **`pays` stops once the answer is known.** The worst angle before only
   falls as more old triangles are read, so the bar `old + P − s` only
   falls; a yes at any point is the final yes. Same answers, as item 3's
   check shows.
5. **The cavity's buffers are reused** inside one `improve` call (the
   judge is a local of the serial pass): no state outlives the call or is
   shared between threads. The value-returning `quality_cavity` stays for
   T-P1.
6. **`@perf`'s method**: the first head run of `bench.py` printed
   `REGRESSION` against the single base run before it; the pooled table
   of 3 pairs counts, as in 20c-1's acceptance. The tile series' slow
   parallel scans from r14 on fell on all three builds in turn and are
   left in the medians; the scan is not under the limit. Why gain 0 leaves
   more to refinement is answered by the counts (ruling 8's "Why"); no
   profile is needed for the ruling.
7. **Open item (b) of the pins ruling needed no commit.** The three CLI
   pins it named pass at the default gain 0, with
   `test_cli_start_quality.py` and `test_cli_start_quality_gain.py` (57
   passed, `@architect`, on this worktree's installed extension), so no
   test needed `--start-quality-gain -1`; `@tester` checks it again on
   the code under review.

#### Mutation round for 20c-2 (`@tester`, 2026-10-07, test commit `600b2acb`, on code `96119508`)

Copied from `@tester`'s hand-back for the round (code review round 1's
S1). The suite is `test_quality_gain`, run on the green code (`d3c939d9`
plus `96119508`) in a scratch copy from `tools/scratch_copy.py 96119508`,
macOS arm64 Release. Before each build the target's object files and
binary were deleted; each log shows the compile line of each test file and
the diff of the planted fault; each binary ran under a timeout; after each
restore the headers were compared byte for byte. Every fault this file
lists was killed except D5b, which the design expected to survive. Three
faults outside the list (QF3, QF4, I1) were not killed by
`test_quality_gain`; `600b2acb` (tests only, +136 lines) kills each with a
new case that fails on its fault and passes on the code (21 cases, 7313
assertions; `test_quality_gain_refine` still passes, 6 cases). Line
numbers in the two tables are in
`tests/cpp/unit/test_quality_gain.cpp@96119508` (T) and
`tests/cpp/property/prop_quality_gain_refine.cpp@96119508` (R); new-case
lines are at `600b2acb`.

**The round ran before the speed fix** (code review round 1's S2).
`65ae0792` rewrote `pays` after the round: the acceptance test
`new_w >= old + gain - s` is now made in two places, inside the loop over
the removed triangles (the early stop, ruling 4 above) and once after it,
so D4 and D5 now sit on two comparisons, not one; and `GainJudge::capped`
(the atan2 skip, ruling 3) did not exist when the round ran, so no fault
was planted in it. Two things carry the round over: the meshes are byte
identical before and after the fix at gain 0 and −1 on Lagan,
Numedalslagen, the tile and the quarter
(`docs/benchmarks/2026-10-07/20c-2/fix-65ae0792/README.md@5b0afe04:139-163`),
and `@reviewer`'s round 1 check of the early stop and the atan2 skip by
argument (the bar only falls, so a yes at any point is final; a corner
inside the margin falls back to the exact call).

**Kill table: the increment file's list**

| # | fault planted (where) | result | killed by test:line, what it checks |
|---|---|---|---|
| D1 | cavity crosses a constrained edge (`quality_cavity`: `is_constrained` test dropped) | killed | T:439, cavity `{0,1}` where `{0}` is expected; T:515, 3 new triangles, not 2 |
| D2 | inexact incircle in the cavity (`quad_flips<pred::FastKernel>` in `quality_cavity` only) | killed | T:478, cavity `{1,0}` where `{1}` is expected (near-cocircular fixture) |
| D3 | foot's second seed missing (the push of `neighbours(t)[on]` removed) | killed | T:499, cavity of 1 slot, not 2; T:177/178, prediction differs from what Lawson wrote (3 slots removed, 5 made), in the line fixture |
| D4 | acceptance reversed (`new_w < old + gain - s`) | killed | T:187, worst angle falls 21.80° to 17.40°; T:590, 642, 680; R:129, R:191 |
| D5 | slack dropped and test made strict (`new_w > old + gain`) | killed | T:641 (the equality fixture gets no insertion); T:680 (`feet` is 0, not ≥ 1) |
| D5b | slack dropped, `>=` kept | **survived**, as the design expected | Only rounding decides it, and only where it puts `new` below `old` on CF2's own edge. Per the increment file, recorded and not chased. |
| D6 | R8 splits at the midpoint, not the foot | killed | T:716/717, 0.861 cells off the foot in col, 0.0154 in row |
| D7 | R8's end check removed (end = 0 in R8's `foot_on` call) | killed | T:768, `line_splits` is 2, not 0 |
| D8 | R8 left on at gain −1 (`split_lines = constraint_feet`) | killed | T:732, 4 splits, not 0; R:110, T-P3 digest differs. 20c-1's foot suites (CF1, CF2, CF3, CF4 and 20b's suite) pass |
| D9 | R8 left on with the feet off (`split_lines = judged`) | killed | T:563 and T:745, 4 splits, not 0 |

**Kill table: further faults, including the new shared pieces**

| # | fault planted | result | killed by |
|---|---|---|---|
| D10 | R8's outline guard dropped (an outline edge with nothing beyond gets split) | killed | T:779, 7 splits, not 0; T:791; R:81, R:129 |
| D11 | the validity callable not asked about R8's foot | killed | T:802, 4 splits, not 0 |
| QF2 | `quad_flips` (shared) uses FastKernel throughout | killed, both sides | Cavity suite: T:476, Lawson's own replay removes 2 slots, not 1; R:110. Lawson side: `test_mesh_lattice_incircle.cpp:894` and `:971`. Line 894 counts `must_flip` decisions (270) that differ from that suite's exact reference copy, so it is a sign disagreement, not a rounding-sized distance. `test_mesh_lawson` survives it. |
| QF3 | `quad_flips`: the branch for a side that is not counter-clockwise in the frame dropped | Lawson side killed; cavity suite **survived**, now killed by 600b2acb | Lawson side: `test_mesh_lawson.cpp:466`. New case "a side collinear in the frame is decided from the triangle across…": fails at 600b2acb line 622, 1 slot removed, not 2 |
| QF4 | `quad_flips`' integer path flips on Cocircular | Lawson side killed; cavity suite **survived** (only R:110 caught it), now killed by 600b2acb | Lawson side: `test_mesh_lattice_incircle.cpp:894`, with 890 failures, then the binary hung past the timeout. I count the hang as no part of the kill; line 894 is the kill. `test_mesh_lawson` survives it. New case "on the integer frame a candidate on the circle of the triangle across does not grow the cavity": fails at line 595, 2 slots removed, not 1 |
| F1 | R8 passes delta (half a cell) for both reach and end | killed | T:561 and T:712 (no splits); R:129–131, R:205 |
| F2 | `constraint_foot` passes infinity for both reach and end | killed | T:547, T:680; R:110, 123, 129; CF1 lines 119 (×9), 170, 310, 316, 348, 402; CF2 272, 328, 407; CF3 174; CF4 230, 252, 272, 416, 480, 513, 567, 599 |
| F3 | inside `foot_on`, the end check reads `reach` (and its mirror, the distance check reading `end`) | killed by the compiler gate | `-Werror=unused-parameter` refuses the build (`end`, or `reach`, unused) |
| I1 | insert routine: the slot across a split edge not seeded or offered | `test_quality_gain` **survived**, now killed by 600b2acb | Already killed elsewhere by `test_mesh_quality.cpp:446/524` and R:129. New case "a line split flips beyond the triangle across the line, as predicted": fails at line 650, the flip of D-E2 is missing from the split |
| I2 | insert routine: written slots never offered to the queue | `test_quality_gain` survived; killed by `test_mesh_quality.cpp:381` (and 438, 441, 478, 527, 546, 599, 705, 720) and R:110/129 | Not added: "every written slot is offered" is increment 20's invariant, which its suite holds; it is not T-P1, T-P2 or LS1 |
| I3 | insert routine: a node on an edge inserted as if strictly inside (`edge = 3`) | `test_quality_gain` survived; killed by `test_mesh_quality.cpp:255` (×4), 473, 720, 727 | Not added, for the same reason as I2: it is the node path's on-edge handling, which increment 20's suite covers |

("I" in the QF4 row is `@tester`.) `edge = on` (without the zero-count
test) is an equivalent mutant: with no zero side, `on` is already 3, so it
was not run. Every kill above is an exact count, a cavity or slot set, a
digest, or a distance far above rounding; none rests on a rounding-sized
oracle violation. "Both suites" for a `quad_flips` fault was read as
Lawson's side (`test_mesh_lawson` with `test_mesh_lattice_incircle`, the
integer path's suite) and the cavity's side (`test_quality_gain`); all
three `quad_flips` faults are caught on both. The QF3 and QF4 cases state
the expected cavity outright, since a fault in the shared `quad_flips`
moves the prediction and `legalise_around` together, and T-P1's
comparison cannot see it.

### 20c-3

Python only, all in `tests/python/`, on hand-made coverages of a few
polygons in a projected CRS (no data files). Before the red step,
`@developer`'s commit raises shapely to 2.2 ("Dependency" in the design);
`@tester` runs on it.

**The repair** (`test_feature_repair.py`, new; not mutation-critical: the
geometry is GEOS's, and what 20c-3 owns is which polygons go in, with which
tolerance, in which order, and each of those is an assertion below):

- **RP1, the slit**: the M5 slit in local coordinates (two polygons that
  share a far vertex, their borders 79 m long and ending 1 cm apart, the
  wedge closed by a 1 cm edge of a third polygon) at S = 0.05: the union of
  the result has no hole; the two long borders are one (the two polygons'
  intersection is a line of length 79 m ± 1 mm); no ring edge shorter than
  0.05 is left; the area that changed polygon (the summed symmetric
  differences) is the wedge's plus the two thin triangles that the moved
  corner sweeps on the third polygon's borders, within 1 mm². At S = 0:
  the wedge is still a hole, and nothing moved. (Checked by `@architect`
  on GEOS 3.14.1 with the third polygon 50 m tall on each side: no hole,
  shortest edge 40.9 m, shared border 79.13 m, 0.4456 m² moved for a
  0.3957 m² wedge; at S = 0 one hole, 0 m² moved. A third polygon only
  1 cm thick collapses to empty, so the fixture must not be one.)
- **RP2, a real gap wider than S**: a strip 2 m wide between two polygons,
  closed at both ends by two more (a strip open at an end is not a hole),
  at S = 0.05: still a hole of the same area (to 1e-9 relative), every
  polygon equal (`shapely.equals`) to the input.
- **RP3, nothing to repair**: a valid coverage with no near miss under S:
  every polygon equal (`shapely.equals`) to the input, and the vertex count
  unchanged.
- **RP4, one source at a time**: two sources, each one polygon, with a
  1 cm gap between them: the gap stays.
- **RP5, land cover only**: a source under the `property` map (no codes)
  whose polygons have the same 1 cm slit gives the same lines as today.
- **RP6, the order**: RP1's coverage with one more vertex on the forest
  side's long border, halfway along and 4.5 mm off the straight line, run
  with `--features-tolerance 2`: that vertex is not in the result. Repair
  first, the two borders are one inner border, which `coverage_simplify`
  simplifies; simplified first, the border is a side of the gap, on the
  coverage's boundary, which `simplify_boundary=False` keeps, and the
  repair then snaps it onto the other border, where it stays. The test
  asserts that premise by calling shapely in the reverse order itself.
  (Checked by `@architect` on GEOS 3.14.1, RP1's fixture with the vertex
  at (39.565, 0.0045): repair first, gone; simplified first, kept.)
- **RP7, the clip**: a polygon much larger than the domain: its label
  polygon after the stage lies within the read region and covers the
  domain; the feature's lines inside the domain are those of today's path.
- **RP8, all off**: with the four switches off (`--features-repair 0
  --no-features-merge-same-class --features-tolerance 0
  --features-outline-snap 0`), `open_features` gives a `FeatureSet` equal
  to today's on the CORINE fixtures already in `tests/python/` (fids,
  masks, codes, lines `equals_exact` at 0, label polygons).
- **CLI**: `--features-repair` defaults to 1 whenever a land-cover map
  is used (question 9's default); refused: negative, non-finite; the input
  row `features_repair_m` and the vertex rows in `--stats` and the record;
  the `features clean-up` timing sub-row present under `features clip`.

**Steps 2 and 3**: a two-polygon coverage of one class merges to one and
keeps the first fid; the partition stays valid and its area changes by less
than tolerance × perimeter; `--features-tolerance 2` with
`--features-repair 0` runs and gives a valid coverage (the coverage check
and its plain error are gone with the old design's step 2).

The outline rule (invariant-critical, mutation testing required: it
rewrites input borders); a CLI test pins its default, 5 m whenever
features are given (ruling 6), and that `--features-outline-snap 0` turns
it off:

- **OR1** a border 3 m inside a straight outline edge, parallel to it for
  100 m, at D = 5: no linework is left within D of the outline except
  where a border crosses the band; the border leaves the outline at a
  point on it. A border 6 m inside is untouched.
- **OR2** two polygons sharing such a border give the same linework from
  either side (the noder sees one border), with the edge given in either
  direction; on a diagonal outline edge too.
- **OR3** a strip 3 m wide between the outline and a neighbour vanishes;
  the neighbour reaches the outline; the area that changes class equals
  the strip's inside the domain.
- **OR4** points the rule puts on the outline are an outline vertex or
  more than D/2 from every other such point; one rounded onto an outline
  vertex keeps that vertex in the path (round 3's prototype defect).
- **OR5** a border crossing the outline at a right angle keeps one
  crossing, within D/2 of where it was, or within D where it goes onto
  an outline vertex (ruling P1); an inlet (the outline's way between
  two moved points longer than twice their distance plus 2 D) is joined
  straight, not round the inlet.
- **OR6** D = 0 is a no-op, byte for byte.
- **OR7** D refused at or above 100 m (the read region's margin), and the
  rule never calls `buffer` on the outline (the rule's own function, called
  directly with `shapely.buffer` patched to raise, on a staircase outline
  of 5 000 or more steps, 30d's case).
- **OR8** (added by ruling G4) the read region does not move the cuts: a
  polygon with one long edge running from 3 m inside a straight outline
  edge, nearly parallel to it, out across the read region's edge, clipped
  by three regions whose margins do not differ by a multiple of D (for
  example 100, 150 and 173.3 m), then `snap_to_outline` at D = 5: the
  linework inside the domain is the same for all three (`equals_exact` at
  1e-9 m). (Checked by `@architect` on `7334b82b` with the polygon
  (3, −537.3), (53.1, −537.3), (53.1, 512.7), (3.4, 512.7) against the
  square 0 to 1000: equal for all three; with the cuts counted from the
  edge's own end instead, the 173.3 m clip differs by 0.62 m.)
- Mutants to kill: the cut not in a fixed order (OR2); the vertex dropped
  at a rounded point (OR4); the stretches on the outline kept (OR1); the
  rounding removed (OR4); the inlet test removed (OR5); the cuts counted
  from the edge's own end, not along its line (OR8).

#### Rulings on 20c-3's green step (`@architect`, 2026-10-08, on `7334b82b`)

Green `7334b82b` on red `6459e87c`: 46 of the 49 red tests pass; the full
suite 6 860 passed, 3 failed, 17 skipped. The three failures and the
other items, ruled; **T** marks a test change for `@tester`, **D** a code
change for `@developer`.

- **G1, `test_cli_features_cleanup.py::TestGiven::test_the_timing_sub_row`
  is a broken test. T.** Its `run()` helper already passes `--stats x.md`;
  the test adds `--stats -`, the last one wins, and `run()` then reads an
  `x.md` that was never written (FileNotFoundError). Drop `--stats -` and
  look for the two rows in the report `run()` returns (its first value):
  the `--stats` file carries the timing table (`stats.py`, `_timings`).
  No code change.
- **G2, `test_cli_mesh_landcover.py::TestFixtures::
  test_a_polygon_clipped_by_the_domain_into_two_pieces`: run it with the
  outline rule off. T.** The test's oracle (`labelled_vtk`) checks every
  triangle against the input polygons as given; the band's borders cross
  the notch's sides, and the rule, on by default, moves them there, as
  designed (OR3: area changes class near the outline). Probed on
  `7334b82b` with the fixture's coordinates: 12.3 m² changes class, all of
  it within 2 D (10 m) of the outline, and each crossing stays one
  crossing. The test is about 16c's labelling of a polygon clipped into
  two pieces, not about the rule, so it passes `--features-outline-snap=0`
  (only that flag; the repair and merge stay at their defaults), with a
  docstring line saying why, as `test_a_lake_in_an_unholed_forest` does
  for the repair. Changing the oracle to the moved polygons would make it
  check the code against itself. No code change.
- **G3, `test_cli_mesh_plain_output.py::TestRecord::test_values_follow_stats
  [features]` wins over `@tester`'s pin. T and D.** Increment 25 already
  settled how a switch is recorded: a text `"on"`/`"off"` equal to its
  `--stats` text (`docs/increments/25-plain-output.md`, "Values", which
  names `snap_to_lines`; `cli.py` writes `"on" if feet else "off"`). A JSON
  boolean would be the record's only one, and `True` passes Python's
  `isinstance(..., int)`, which is how it broke the integer rule. So:
  `@tester` changes `test_cli_features_cleanup.py`'s two asserts to
  `record["features_merge_same_class"] == "on"` and `== "off"`, pins the
  `--stats` text as `on` and `off` in `test_the_merge_row_says_on_or_off`,
  and drops the module docstring's "the merge as a JSON boolean";
  `@developer` writes the row as `"on"`/`"off"` in `cli.py` and removes the
  `isinstance(value, bool)` branch from `run_record._entry`.
  `FeatureRequest.merge_same_class` stays a `bool`: the model is typed, the
  record is text.
- **G4, the read-region clip and the outline rule. Design corrected; T.**
  The design said the clip cannot reach the rule; that held for the new
  edges the clip makes (at least 100 m outside the outline, D under
  100 m), not for the input edges it shortens, which can still run within
  D of the outline. With cuts counted from an edge's end, a different
  region moved them, and 30d's byte-identity test
  (`tests/python/test_grow.py::test_7b_projected_staircase_with_corine_is_byte_identical`)
  failed. `@developer`'s fix, cuts at multiples of D along the edge's
  line counted from the foot of the CRS's origin, in the direction of the
  edge's lexicographically last end, is ruled correct: both polygons
  sharing an edge get the same cuts, and a shortened edge keeps its cuts.
  **To rounding, not bit for bit**: the shortened edge's unit vector and
  offset are recomputed from the new end, so a cut can move in its last
  bits, and a cut that the rule does not move onto the outline carries
  that into the linework. That costs nanometres, and only for an input
  edge longer than about 95 m that both crosses the region's edge and
  comes within D of the outline. The design sentence ("Clipped to the read
  region first") is corrected, and OR8 above pins it with a unit test, so
  the slow real-data test is not the only guard. Scale: D in metres in the
  computation CRS, checked at D = 5 on a 1 km square, coordinates up to
  about 540 m from the origin; the multiples are exact in doubles up to
  2^53 D, far beyond any UTM coordinate.
- **G5, size: 331 counted lines against about 150 (+121 %). Accepted;
  one PR.** `python3 tools/count_loc.py 483221d2 7334b82b`:
  `feature_input.py` 260, `cli.py` 63, `run_record.py` 8. Of
  `feature_input.py`'s, the outline rule (`OutlineSnap` to the end of the
  file) is 188 counted lines against the ~70 estimated; the stage itself
  (`_add`, `_clean`, `_polygonal`) 68, much of `_add` moved from the old
  loop (19 lines removed). The estimate was low because it was not taken
  from the prototype it was measured with: `snapline.py` in
  `../rasputin_scratch/20c-prototype/` is already 146 non-blank lines, and
  the design then asked for more than it does (rounding onto outline
  vertices, cuts the same from either side, the inlet test, a polygon
  rebuilt for the labels and the area changed). Under 700, so no split.
  The design should have named one: the land-cover stage (repair, merge,
  simplification, four flags) and the outline rule were separable, each
  with its own tests and gates, and the rule carried the estimate's risk.
  Splitting now would cost a review round and a CI run for nothing the
  limit asks, so it stays one PR. G3 removes a few lines; no other
  shrinking is asked.
- **G6, `@developer`'s pinned assumptions.**
  (a) *An inlet joined straight drops the parts of that line that lie on
  the outline at either end*: kept; it is the design's "stretches on the
  outline dropped" applied to the straight join.
  (b) *Arc-length rounding*: a placed point goes to the nearest multiple of
  D along its outline ring, counted from the ring's first vertex, then to
  a vertex within D/2: kept. The outline is the same object for every
  polygon, so every polygon rounds alike.
  (c) *"Within D" inclusive, ties to the lowest segment index*: kept;
  deterministic, and GEOS's `dwithin` and `query_nearest(max_distance=)`
  are both inclusive.
  (d) *The merge unites a whole class into one feature under its first
  fid, a MultiPolygon if its parts are apart*: kept; the design's "per
  class" and "the first fid of its class" read that way, and it is
  simplest. Noted: labelling ranks overlapping polygons by area (16c D2,
  the smallest wins), so a polygon of another source now wins an overlap
  against a merged class more often than against the single polygon it
  replaced. That is the useful direction (the finer source wins); within
  one source the repaired coverage has no overlaps, so the ranking does
  not arise there.
  (e) *Dropped polygons count as "outside"*: **changed for one case. T
  and D.** A polygon the clip leaves nothing of lies outside the read
  region, so outside the domain, and stays counted there. A polygon the
  repair gives wholly to its neighbours is not outside the domain, and
  the stderr line would say so; it is counted `empty` instead (one word
  in `_clean`). `@tester` pins it with RP1's degenerate case, a third
  polygon 1 cm thick (which the repair empties, as noted in RP1): `empty`
  rises by one and `outside` does not.
  (f) *The vertex, area and timing rows appear only when the stage ran*:
  kept, as the design's off switch requires (all four off: today's output).
  (g) *The "borders moved" count was not built*: kept out; the design
  drops it. The area row carries the rule's effect, a border count needs
  a definition of its own (vertices? chains?), and the PR is over its
  estimate.
- **G7, the help panel.** The four flags sit in a new help panel,
  "Land-cover clean-up" (`cli.py`, `CLEAN_UP`), because in the main panel
  the long `--[no-]features-merge-same-class` cut off older flags' names
  in an 80-column help table. Noted; the design now says so ("Recorded").

**Next**: `@tester`'s test changes (G1, G2, G3, G4's OR8, G6 (e)), then
`@developer` for G3 and G6 (e), then `@tester`'s mutation round on the
outline rule (invariant-critical; the six targets above), then `@perf`'s
short timed check (at most 15 minutes, "Speed judgment") and the gate runs
of "PRs and gates", then review.

#### Mutation round for 20c-3 (`@tester`, 2026-10-08, test commit `5cc0639c`, on code `e07d5921`)

Copied word for word from `@tester`'s hand-back for the round, as for
20c-2; pins ruled below ("Rulings on 20c-3's mutation round").
It ran on `e07d5921`, before T1 and T2 rewrote `_cut` and `ring`; the
killing tests are unchanged since (code review round 1's R5).

- **Where it ran:** a scratch copy made by `tools/scratch_copy.py e07d5921`, with the worktree's `.venv`. Every mutant was a single text replacement in the copy's `src_python/tin_engine/feature_input.py`. After each run the file was restored and compared byte for byte with the original. The copy is now removed.
- **What ran on each mutant:** `test_feature_outline_rule.py` and `test_feature_repair.py`. No mutant failed a repair test.
- **Line numbers:** they are in `tests/python/test_feature_outline_rule.py@5cc0639c`. I re-ran every mutant against that final file to get them.
- **Code line numbers:** the "where" column cites `feature_input.py@e07d5921`.

| # | fault planted (where) | result | killed by test:line, what it checks |
|---|---|---|---|
| M1 | the cut not in a fixed order: `_cut` never flips to the lexicographically last end (`flip = False`, :767) | killed | :184 `test_or2_one_border_from_either_side[same-straight]` and `[same-diagonal]`. The two polygons' shared segments differ in the last bit (500105.4588256424 against …425). The reversed cases pass, as they should. This kill is rounding-sized: under ruling G4's cuts along the line, direction moves only rounding. But OR2's oracle needs both sides of a border to have bit-identical coordinates ("the noder sees one border"), so it is the right oracle. |
| M2 | the vertex dropped at a rounded point: `_place` loses the second rounding, "a multiple of D within D/2 of an outline vertex goes to the vertex" (`self._vertex(r, m) or` removed, :791) | **survived**, now killed by `5cc0639c` | New test `test_or4_a_rounded_point_next_to_a_vertex_goes_onto_it`, fails at :262: the crossing ends at (60, 0), not at the vertex (61, 0). It survived because `TestTheCorner`'s vertex V = (100, 0) sits at a multiple of D along the outline, so the first rounding already lands on V. `TestRounding` has no point whose own position is more than D/2 from the vertex at 61 while its rounded multiple is within D/2 of it. |
| M2b | same target, second form: a join drops the outline's own vertices between two placed points and goes straight (`_join`'s along-the-outline branch returns `[(q[2], True)]`, :828) | **survived**, now killed by `5cc0639c` | New test `test_or5_a_shallow_notch_is_followed`, fails at :327: the polygon after the rule covers 0 m² of a 2 m by 1 m notch, not 2. It survived because the straight join is still flagged as lying on the outline, so it is dropped from the linework, and every existing test checks linework only. The defect shows in `OutlineSnap.polygons`, which `_Tally._add` receives (`feature_input.py:441-446`), and in `area_changed`. |
| M2c | both roundings onto a vertex removed (`out[i] = (r, m, self._at(r, m))`, :791) | killed | :241 `test_or4_placed_points_are_apart` (the vertex at 61 is never used); :262 |
| M3 | the stretches on the outline kept: `_chains` never cuts (`if True:` for `if not any(on):`, :676) | killed | :135 OR1 (linework left in the band), :207 OR3 (strip length not 0), :281 OR4 corner, :310 OR5 inlet |
| M4 | the rounding removed (`m = s`, :790) | killed | :246 two placed points 2.0 m apart, under D/2 = 2.5; :262 |
| M5 | the inlet test removed (`if True:` for `if min(ahead, length - ahead) <= 2 * gap + 2 * self.d:`, :827) | killed | :309 `test_or5_an_inlet_is_joined_straight`: the mouth's linework is empty, distance NaN, not ≤ 1e-9 m |
| M6 | the cuts counted from the edge's own end, not along its line: `_cut` uses `k = 1 .. ceil(len/D) - 1`, `cuts = lo + k D u` (:770-772) | killed | :432 `test_or8_the_same_lines_inside_the_domain`: the 173.3 m clip gives the vertex (53.1, 6.7) where the 100 m clip gives (53.1, 10), 3.3 m apart. The design measured 0.62 m for its own version of this mutant; mine counts from the lexicographically first end, so the number differs, and both are far above rounding. |

Every kill is an exact coordinate, an area or a length far above rounding, except M1, whose oracle is bit-identity on purpose. I also changed `TestTheCorner`'s docstring: it claimed the test catches "round 3's prototype defect", which it cannot, so it now points to the new test.

**Added in the fix round (`@architect`, from `@tester`'s red step `51c24772`).**
One fault more, ruling T1's target, "`_cut`'s end test removed": **killed** by OR9,
`test_or9_no_edge_shorter_than_half_d` (both directions). Shown in the red
step, not as a separate planted run: the code before T1 (`eb810a90`) is
exactly that mutant, OR9 failed on it (a 2.4 cm edge beside V), and passed
with T1 (`35b7a873`). T2's test (`test_t2_the_far_part_passes_through`) is
not mutation-critical (ruling T2).

#### Rulings on 20c-3's mutation round (`@architect`, 2026-10-08, on `5cc0639c`)

`@tester`'s three pins from the round above, ruled. No code change follows
from any of them.

- **P1, how far a crossing may move: OR5's wording corrected.** Step 4
  (above, "Steps 2 to 4") sends a multiple of D within D/2 of an outline
  vertex to the vertex, and `_place` does that
  (`src_python/tin_engine/feature_input.py@e07d5921:775-791`): a point goes
  first to an outline vertex within D/2 of its own place, else to the
  nearest multiple of D along the outline, then to a vertex within D/2 of
  that multiple. Each of the last two moves is at most D/2 along the
  outline, so a point on the outline moves at most D along it, and no
  farther in a straight line. OR5 said "within D/2 of where it was", which
  only the case without a vertex nearby meets; its fixture is that case,
  so no test disagreed. The design is right and OR5's words were short;
  OR5 now reads "within D/2 of where it was, or within D where it goes
  onto an outline vertex", and `test_or4_a_rounded_point_next_to_a_vertex_goes_onto_it`
  (a crossing moved 2.6 m at D = 5, onto the vertex) is inside it.
- **P2, the rule has two outputs. Noted; the tests now check both.**
  `snap_to_outline` returns `OutlineSnap.lines`, the linework the noder
  gets, and `OutlineSnap.polygons`, the label polygons
  (`_Tally._clean`, `src_python/tin_engine/feature_input.py@e07d5921:435-446`)
  from which `area_changed` is summed. A fault in a stretch that is dropped
  from the linework (M2b, a straight join flagged as on the outline) shows
  only in the polygons, so a test of the rule that reads only the linework
  cannot see it. `test_or5_a_shallow_notch_is_followed` reads the polygon;
  any later test of the rule states which of the two outputs it checks.
- **P3, OR2 is bit for bit; ruling G4's "to rounding" is about another
  object. Confirmed.** Within one run, the two polygons sharing a border
  reach the rule with that border's vertices bit-identical: `_clean` runs
  `coverage_clean` before the rule in every run where the rule runs
  (`src_python/tin_engine/feature_input.py@e07d5921:422-435`; at
  `--features-repair 0` too), and its noder makes shared vertices one.
  `_cut` orders each edge by its lexicographically last end before any
  arithmetic, so both sides compute the same cuts from the same two
  points in the same order, which gives the same bits. That is what OR2's
  "the noder sees one border" asks, and M1's kill (the last bit apart) is
  the oracle working, not a rounding-sized flake. G4's "to rounding, not
  bit for bit" compares one edge clipped by two different read regions,
  so two different runs; it does not weaken OR2.

**Next**: `@perf`'s short timed check (at most 15 minutes, "Speed
judgment") and the gate runs of "PRs and gates", then review.

#### Rulings on 20c-3's timed check and gate runs (`@architect`, 2026-10-08, on `16143f4c` and `4d52dec0`)

`@perf`'s check (`docs/benchmarks/2026-10-08/20c-3/README.md@16143f4c`)
and profile (`profile.md@4d52dec0`, same folder). **D** marks a code
change for `@developer`, **T** a test for `@tester`. The probe behind
these rulings: the outline rule's inputs taken from the CLI at
`4d52dec0` (the four runs below: both catchments, with the defaults and
with the gate's flags), and each candidate change applied by patching
the module from a script outside the tree; AC power, one run each, no
warm-up. The scripts were in the session's scratchpad and are not kept;
the changes are described in full below.

- **T1, Numedalslagen's outline-rule gate: a defect in ruling G4's cuts,
  not a gate on the wrong basis. D and T; the gate stands.**
  *The basis matches.* The design's figures (4 slivers, 0 near, 0.832°;
  M6 (d)) are the prototype's at tolerance 2 m, no merge, no repair, on
  20c-2's switches; `@perf`'s gate runs use the same flags
  (`scripts/gates.sh`), and the tolerance-only mesh they start from
  agrees with the prototype's within the 20c-2 build's own difference
  (29 slivers and 25 near, against 32 and 27; 20c-2 built is 2 slivers
  under the prototype, "What 20c-2 built"). "Near" is the same measure
  in both (`scripts/gate.py` and the prototype's `strip.py`: the
  sliver's centre within 20 m of the domain's boundary).
  *The built rule differs from the prototype in its cuts.* The
  prototype (`snapline.py`, `cut`) split an edge into ⌈L/D⌉ equal pieces,
  each longer than D/2. Ruling G4 put the cuts at multiples of D along
  the edge's line, counted from the foot of the CRS's origin, and such a
  multiple can fall any distance from the edge's own vertex. Where that
  vertex lies just beyond D and the next cut within D, step 3 keeps the
  cut (the step off the outline), and the linework gets an edge of
  centimetres, nearly in line with the input edge; the mesh makes
  slivers on it. In `@perf`'s mesh (`g-tol2snap5-num`), 10 of the 11
  near slivers stand on such an edge, 2.9 to 21 cm long. The worst
  (0.0814°, centre 219 550.8 6 599 776.3) has the input vertex
  219 551.178 6 599 769.462, 5.52 m from the outline, and the cut
  219 551.207 6 599 769.467, 2.9 cm from it. The eleventh (shortest side
  1.1 m, no land-cover line within 30 m) is also gone with the fix.
  *The fix (D):* `_cut` drops a multiple within D/2 of either end of the
  edge. The pieces at an edge's ends are then D/2 to 3D/2 long and those
  inside D; an edge shorter than D gets no cut, as in the prototype. G4
  still holds. A cut within D/2 of a clipped end lies more than 97 m
  from the outline, so it never moves and is never kept, and OR8 stands.
  The test reads both ends alike, so OR2 stands too. Probed: Numedalslagen
  4 slivers, 0 near, worst 0.832158°, 1 016 301 triangles; Lagan 61, 0
  near, 0.006697°, 743 616: M6 (d)'s figures exactly. The rule's test
  files (`test_feature_outline_rule.py`, `test_feature_repair.py`,
  `test_cli_features_cleanup.py`, `test_grow.py`) pass with it patched in
  (76 passed; the patch was shown to take effect by a `_cut` returning no
  cuts, which fails 8).
  *The test (T), OR9:* at D = 5, a straight outline edge and a border
  edge from a vertex V 5.5 m inside it to a point 1 m inside, placed so
  that a multiple of D along the edge's line (counted as `_cut` counts)
  falls within 5 cm of V. No edge of the rule's linework inside the
  domain is shorter than D/2, with the edge given in either direction.
  On `4d52dec0` that edge is under 5 cm, so the test fails. Mutation
  target: the end test removed (OR9 kills it).
  *Whose defect:* ruling G4 ruled the cut correct without comparing its
  piece lengths with the prototype's. `@developer` built what G4 said.
- **T2, the speed: rebuild and measure only what moved. D and T.**
  *What costs the time (profile).* With the merge on, each class is one
  polygon (Numedalslagen: up to 590 parts and 125 000 vertices). For
  every class with a moved ring, the rule rebuilds the whole class with
  `union_all` over all its parts, moved or not (2.17 s), and measures the
  area changed with `symmetric_difference` of the whole class, old
  against new (2.20 s).
  *The design.*
  (1) A part with no moved ring is passed through unchanged, as now, and
  no longer goes into a union. The polygon after the rule is the
  unchanged parts plus the parts of `union_all` over the rebuilt parts
  only, with no union when one part was rebuilt. If that polygon is not
  valid, the rule falls back to today's union of all parts. The guard is
  there because validity is argued, not proven: a rebuilt part changes
  only within reach of the outline, where an unchanged part has no edge.
  It cost about 0.1 to 0.2 s on the changed polygons in the probe. The
  union over the rebuilt parts stays: without it, the probe lost about
  1 000 km² of one Lagan polygon, because `make_valid` output overlapped.
  (2) `_Outline.ring` also returns, for each point it gives back, the
  index of the input vertex it is (−1 for a cut, a moved point or an
  outline vertex). Between two consecutive kept input vertices whose
  stretch changed, the old stretch and the new one form a closed loop
  that bounds the region that changed there. A ring with no kept input
  vertex contributes its old and its new polygon. Each loop is made
  valid, and `area_changed` is the area of the union of all loops of all
  polygons, intersected with the outline. Today it is the union of the
  symmetric differences, intersected with the outline. Neighbours moving
  one border both give the same loop, and the union counts it once, as
  today.
  *Outputs identical (probed on all four inputs).* The lines are
  identical bit for bit, because they come from the same `ring` and
  `_chains`. Every polygon is `shapely.equals` to today's, so the labels
  are the same. `area_changed` is equal to four decimals in m²: 13 215.9606
  and 16 074.3979 with the defaults, 13 467.6887 and 16 185.5619 with the
  gate's flags, for Numedalslagen and Lagan.
  Existing tests that pin this: OR1 to OR5 and OR8 (lines), OR5's notch
  and RP7 (polygons), OR3 (area to 1e-6 m²) and OR6 (area 0).
  *The test (T), not mutation-critical:* a polygon of two parts, one far
  from the outline and one with two separate stretches within D of it.
  The far part comes back `equals_exact` at 0 to its input. The area
  equals the overlay measure that the test computes itself (the input
  and output's symmetric difference, intersected with the outline) to
  1e-6 m².
  *Commits:* T1 and T2 as separate green commits, so that `@perf` can
  check that T2 leaves the mesh the same bit for bit.
  *Target* (the probe, both fixes, defaults, against `@perf`'s figures):
  the rule alone, 5.9 → 2.7 s on Numedalslagen and 4.4 → 2.4 s on Lagan.
  The clean-up 8.78 → 5.51 s and 8.04 → 5.92 s. The whole run
  14.10 → 10.57 s (master 6.81 s, **+55 %**) and 22.87 → 20.59 s (master
  17.64 s, **+17 %**). `@perf`'s limit for the re-time: clean-up at most
  6.0 s on Numedalslagen and 6.5 s on Lagan, one run each.
  *Size:* T1 about 2 counted lines, T2 about +25 (`ring`'s index list,
  the loops, the assembly and its guard, less the overlay). About 360 in
  all, under 700.
- **T3, the hotspot: recorded, and question 10.** After T2,
  Numedalslagen's clean-up is still about 52 % of its run (Lagan's about
  29 %, its `features clip` 49 %). The design's 1.1 s was M7's clip and
  repair alone (0.34 + 0.76 s, scratch probe). It never timed the
  outline rule or the clip of the lines to the domain, and in the built
  stage the clip and the repair alone take about 1.9 s on Numedalslagen
  (0.75 + 1.18 s, profile). So the stage cannot get near 1.1 s with
  the rule on. What remains on Numedalslagen after T2 (the profile, less
  T2):
  - the outline rule, about 2.4 s: its Python loop over ring edges
    0.8 s, `make_valid` of the big rebuilt parts 0.6 s, their union
    0.7 s, the loops 0.45 s;
  - `coverage_clean`, 1.2 s;
  - the clip of the lines to the domain (`_add`), 1.0 s;
  - the clip to the read region, 0.75 s;
  - the merge, 0.1 s.
  Further work, not designed here, could cut that: a vectorised ring
  loop, no whole-part repair and union for the big parts, and no domain
  clip for chains the rule already knows to be inside. It might save
  1.5 to 2 s; that is not measured. The gates do not change. Ola
  accepted the time for 20c-3 (ruling 10); that further work is the
  ROADMAP row "Speed up the land-cover clean-up step".

#### Rulings on the fix round (`@architect`, 2026-10-08, on `8716f062` and `41db7333`)

**Outcome of T1 and T2** (`@perf`'s re-time,
`docs/benchmarks/2026-10-08/20c-3/README.md@41db7333:91-130`). T1: both
outline-rule gates met at the design's figures ("PRs and gates"). T2: the
mesh is T1's byte for byte; the clean-up 5.92 s on Lagan (limit 6.5 s)
and 5.48 s on Numedalslagen (limit 6.0 s), met; the whole run 20.59 s
and 10.51 s, +17 % and +54 % against master (question 10; Ola
accepted it, ruling 10).
Defaults: slivers under 1° 437 → 180 and 296 → 158, worst angle
0.002384° → 0.1917° and 0.008623° → 0.6302°.

`@developer`'s five pins on T2, ruled. One changes code.

- **P1, a ring with no kept input vertex adds its old and its new polygon
  to the area moved. Changed: D and T.** That is T2's wording, and it is
  wrong whenever the two polygons overlap. Probe (a 1 km square outline,
  D = 5, polygon A the 2 × 10 m rectangle (0,100)-(2,110) on the outline,
  B the square less A): `area_changed` is 1 000 000 m² on `8716f062`
  against 20 m² on `35b7a873`. B's eight vertices all lie within D, so
  none is kept, and B's old and new polygons together are the whole
  square. *Fix (D):* such a ring gives the symmetric difference of its
  old and its new polygon, each made valid (a probe patching `_loops` so
  gives 20 m²). It can only lower the figure (the difference lies in the
  union), so `@developer` checks the four real inputs' `area_changed`
  still equal T1's to four decimals. *Test (T):* that fixture,
  `area_changed` 20 m² to 1e-6, and A empty, B the whole square; red on
  `8716f062`. Not mutation-critical. The area is the run record's
  `land_cover_area_moved_m2` only, never the mesh, so no `@perf` rerun.
  Whose defect: T2's design (`@architect`).
- **P2, nested `make_valid` output flattened two levels. Kept.**
  `make_valid` and `union_all` return at most a collection of
  (multi-)geometries, so two levels of `get_parts` reach every polygon;
  `_polygonal` drops the lines.
- **P3, the fallback union is not exercised by any test. Kept.** The
  guard covers a case argued not to arise (a rebuilt part changes only
  near the outline, where an unchanged part has no edge); no fixture that
  reaches it is known. It stays untested, and the review notes it.
- **P4, the holes of a part whose shell collapsed still add loops. Kept,
  a known limit.** The part is gone, so the figure can count the hole's
  area as moved, which overstates `area_changed` by at most the area of
  the holes of parts thinner than about 2D. Related and older (the same
  on `35b7a873`): a band narrower than D along the whole outline comes
  back as the whole domain, not empty, because its shell and hole both
  go onto the outline. Not fixed in 20c-3; the review notes it.
  *The band case is fixed after all* by code review round 1's R1 (below):
  with `make_valid`'s "structure" method the band comes back empty and
  its neighbour as the whole square (probe on `1e011cec`, patched).
- **P5, the citation `tests/python/test_grow.py:185`. Kept, not stale.**
  Its docstring names the commit the lines are read at (`f81b20b7`), and
  `src_python/tin_engine/feature_input.py@f81b20b7:158-167` is master's
  `source_region`, the function the test copies. `check_citations.py`
  lists it only because this branch edits that file. No edit needed.
  `@tester` may rewrite it in the pinned form just used (full path, `@`,
  commit) to take it off that list (one line, a test file); the short
  path with `@` is reported broken.

#### Rulings on 20c-3's code review round 1 (`@architect`, 2026-10-08, on `05b0fd7a`)

The review is the last entry under "Review". Ola was asleep; each ruling
is on the default the main session proposed, and none is a question for
him.

- **R1 (review B1), a hole the outline rule puts onto the outline is
  filled. Fixed: D and T.** `snap_to_outline` rebuilds a changed part as
  `Polygon(shell, holes)` and repairs it with `shapely.make_valid`, whose
  default "linework" method keeps a hole's area once the hole touches the
  shell along a line; with the merge on, the smallest polygon wins the
  label, so a lake inside a thin ring of a rarer class takes the ring's
  class. Reproduced on `1e011cec` with the review's fixture (1 km square
  outline, D = 5 m; water = box(3,200,103,300) ∪ box(500,500,800,700);
  ring = box(2,190,113,310) less the first box; rest = the square less
  both): the ring comes back as 13 560 m² and holds the lake's centre
  (53, 250). *Fix (D):* `shapely.make_valid(shell, method="structure")`
  at that one call; the call in `_loops` (ruling P1) stays, as it repairs
  a single ring with no holes and gives only an area. Patched, the three
  polygons are 70 300, 3 260 and 926 440 m², a partition of the square,
  the ring no longer holds (53, 250), and `area_changed` is 340 m² both
  ways. The same patch fixes P4's band case (a band narrower than D along
  the whole outline, its inside a second polygon): the band comes back
  empty and the inside as the whole square, `area_changed` 11 964 m²,
  against both polygons the whole square today. *Test (T):* the review's
  fixture (the ring's area 3 260 m² to 1e-6 m², the ring not covering
  (53, 250), the three areas summing to 1e6 m²) and the band case (the
  band empty, the inside 1e6 m²); both red on `1e011cec`. Not
  mutation-critical: one keyword. *The mesh:* "structure" and "linework"
  differ only on a rebuilt part that is invalid, where a hole touches or
  crosses its shell (filled before, cut out now) or a ring overlaps
  itself (an overlap counted once, not cut out). The review's patched run
  found no filled hole on Lagan or Numedalslagen, but did not look for
  the second kind, so the meshes should not change and `@perf` checks
  that they do not (R3 (a)). Whose defect: the green step's design
  (`@architect`), which named `make_valid` without its method.
- **R2 (review B2), prose the branch made untrue. Fixed in this commit:**
  PR #214 merged (`4abacf65`) in the status and the ROADMAP row; the
  "Next" line; the LOC table (356 built, `python3 tools/count_loc.py
  483221d2 1e011cec`); the two passages that waited on 30d (merged,
  `769ed992`).
- **R3 (review B3), the gate's unrun parts.** The status said "every
  outline-rule gate met"; what was met is the outline-rule rows. The rest
  of the 20c-3 gate row, each ruled:
  - *Run (`@perf`, on `@developer`'s R1 green commit, one session of at
    most 15 minutes, AC power, the catchment arguments of
    `docs/benchmarks/2026-10-08/20c-3/scripts/gates.sh`):*
    (a) **the mesh unchanged by R1**: default runs on both catchments;
    `.vtk` SHA-256 equal to T2's (Lagan `ef67e65c…`, Numedalslagen
    `1f11f291…`). If either differs, stop and report it to `@architect`.
    (b) **the repair on top of the outline rule, on the fixed code**:
    `--no-features-merge-same-class --features-tolerance 2
    --features-outline-snap 5 --features-repair 0.05` on both catchments,
    against the outline-rule row's thresholds (Lagan: at most 90 slivers,
    at most 5 within 20 m of the outline, triangles at most +1 % of the
    tolerance-only mesh `g-tol2-lagan`, worst angle at least 0.95 ×
    0.006697°; Numedalslagen: at most 10, at most 5, at most +1 %, at
    least 0.1°), and on Lagan `scripts/corners.py` (the slit's two open
    corners not both mesh vertices). The repair on the tolerance-only
    run (`--features-outline-snap 0`) is not rerun: everything changed
    since `@perf`'s runs on `745bcb0d` is inside the outline rule, which
    D = 0 skips (`git diff 745bcb0d 1e011cec -- src_python` touches only
    `snap_to_outline` and `_Outline`).
    (c) **the slit's far corner**, on Lagan's two repair meshes (the new
    one from (b) and `g-tol2rep-lagan.vtk` from the earlier runs): no
    triangle whose smallest angle is under 1° has a corner within 1 mm
    of 413 677.295 6 331 664.081. Restated because the `.vtk` does not
    say which edges are constraints: the gate's wording asked for no
    smallest angle between two constraint edges there; a sliver at that
    corner whatever its edges is the stricter test.
    (d) **no shared border left unmatched**: from the run in (b), the
    polygons `snap_to_outline` returns (a driver that wraps the function
    and runs the mesh command in-process), measured as the prototype's
    `../rasputin_scratch/20c-prototype/covfar.py` does (mismatched
    border inside the domain, more than 1 mm from its outline): 0.000 m
    on both. The probe can fail: the same measure with one vertex of one
    shared border moved 0.5 m gives more than 0.
    (e) **population 3 carries no land-cover line with the merge on**:
    Lagan with the defaults, and again with
    `--no-features-merge-same-class`; from the same driver,
    `snap_to_outline`'s lines moved to EPSG:3035, the total length of
    segments with both ends within 1 m of the line N 3 811 923.31 or of
    the line E 4 585 680.19 (M2's two population-3 lines), inside the
    domain. Merge on: 0 m. Merge off: more than 0 (the probe that can
    fail). Read on the land-cover lines, not the mesh, because the lines
    are the only source of land-cover constraint edges and the `.vtk`
    does not mark constraints.
  - *Ruled out: raising the thresholds to the measured figures.* The
    measured figures after T1 are the design's own (61 and 4 slivers, 0
    near the outline, 0.006697° and 0.8322°, "PRs and gates"), and the
    thresholds were set from those figures with the margin for the
    production code; there is nothing new to raise them to. The repair
    rows keep "every threshold of the run without it" (M7 gave no
    figures of their own).
  - *Outcome of (e) (`@architect`, on `@perf`'s `63065a15`): the check's
    wording was wrong, not the code.* With the merge on, 26.194 m remains:
    13.097 m of border, counted once for each of its two polygons, in four
    stretches of 0.69 to 9.46 m
    (`docs/benchmarks/2026-10-08/20c-3/stats/r3/e-where.txt`), between
    merged polygons 10 and 12 and 10 and 11. With the merge on, each
    polygon is one class (`feature_input.py`'s `_clean`, "one polygon per
    class"), so these are borders between different classes. The source
    data agrees: read from the CORINE GeoPackage in EPSG:3035 at each
    stretch's middle, 0.3 m to either side, three stretches have coniferous
    forest (312, OBJECTID 2000793) against transitional woodland-shrub
    (324, OBJECTIDs 1933839 and 1933737), and the fourth coniferous forest
    (312, 2373011) against mixed forest (313, 1929199). A polygon of another
    class reaches the straight cut there, and its border runs within 1 m of
    the line. Population 3 is the cut through one class (312 on both sides),
    and none of it is left. Corrected wording of (e): with the merge on, no
    land-cover line within 1 m of either line has the same class on both
    sides: 0 m, met; with the merge off, more than 0 (161 164.495 m, the
    probe that can fail). No code change and no new test.
- **R4 (review S2), the record row's wording. Taken: D and T.**
  `run_record.py`'s label for `land_cover_area_moved_m2` says "changed
  class", but with the merge off the outline rule can move area between
  two polygons of one class. *Change (D):* "Land-cover area inside the
  outline that the outline rule gave to another polygon, m2". *Test (T):*
  the label, red on `1e011cec`.
- **R5 (review S1), where the mutation round ran. Taken.** The round ran
  on `e07d5921`, before T1 (`35b7a873`) and T2 (`8716f062`) rewrote
  `_cut` and `ring`. The killing tests are unchanged since: between
  `5cc0639c` and `1e011cec` the test file gained tests and lost only
  docstring lines (`git diff 5cc0639c 1e011cec --
  tests/python/test_feature_outline_rule.py`). The review found that the
  record covers every mutation target plus T1's end test.
- **R6 (review S3), P1's check on the real inputs, recorded.** From
  `@developer`'s hand-back for P1's green step: `area_changed` on the four
  real inputs, to four decimals, Numedalslagen with the defaults
  13 243.1104 m², Lagan with the defaults 16 050.3242 m², Numedalslagen
  with the gate's flags 13 481.8453 m², Lagan with the gate's flags
  16 165.0077 m², the same under T1, T2 and the P1 fix. No ring on those
  inputs reaches P1's branch (a ring with no kept input vertex), so the
  equality says only that the fix does not touch them.

## `@perf`'s acceptance

20c-1 and 20c-2 touch `include/terrain/mesh/` and
`include/terrain/refinement/`, so `docs/increments/README.md`'s acceptance
applies: `tools/bench.py` on the 1 m tile and the thread-scaling sweep,
compared with master's stored run on the same power state. The 1 m tile has
domain polygons (`tools/bench.py`'s `domains`), whose outlines are
constraints, so its mesh hash may change wherever a quality node lies near
an outline; `@perf` reports the triangle count and the share under 1° beside
the timings.
Refine time within noise or better (M3: the whole refine is under 1.4 s of
the 15 to 90 s runs; 20c-2 inserts a third fewer quality nodes). Plus the
gate runs above, and the split-phase limit ("The split-phase limit, judged
after the fix"). 20c-3 is Python only and needs no `bench.py` run; its
gate runs are `@perf`'s, and so is its short timed
check ("Speed judgment", in its design).

**20c-1, done.** First acceptance at `34b546c4`
(`docs/benchmarks/2026-10-07/20c-1/README.md@b9e2d482`): gates pass, refine
time regressed (+6 to +24 % on the 1 m benchmark). Re-time after the fix
at `69f37d1c`
(`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/README.md@c9f62f7c`): the 1 m
benchmark accepted, −0.3 to −2.0 %, meshes byte-identical to master's;
gates pass at the same figures.

**20c-2, done.** First acceptance of `600b2acb`
(`docs/benchmarks/2026-10-07/20c-2/README.md@fa136aa7`): every mesh gate
passes; start quality about +60 %, Numedalslagen's whole run +3.2 %, the
split phase over 2 %. Re-time after the speed fix `65ae0792`
(`docs/benchmarks/2026-10-07/20c-2/fix-65ae0792/README.md@5b0afe04`):
meshes byte-identical to `600b2acb`'s, the 1 m benchmark accepted, the
whole run −0.04 % and −0.08 %, the split phase met per split (ruling 8),
start quality +23 % and +20 % (reported, not gated).

## LOC

Counted in `CLAUDE.md` §2's unit.

| PR | file | est. |
|---|---|---|
| 20c-1 | `mesh/constraint_foot.hpp` (new; 20b's `foot_of`, `foot_fits` moved in, neighbour search) | ~45 |
| | `mesh/quality.hpp` (foot branch, counters, callable type, option) | ~25 |
| | `refinement/refine.hpp` (owner handling, callable; minus the moved helpers) | ~0 net |
| | `refinement/refine_points.hpp` (foot branch, footed set, z along the strip, counters, option) | ~45 |
| | `bindings/core.cpp`, `_core.pyi`, `final_check.py`, `edge_strip.py`, `cli.py`, `stats.py`, `run_record.py` | ~35 |
| | **total** | **~150** |
| 20c-1 built | `python3 tools/count_loc.py c074f900 69f37d1c`: `constraint_foot.hpp` 77, `quality.hpp` 20, `refine.hpp` −21, `refine_points.hpp` 53, `strip_scan.hpp` 0, `bindings/core.cpp` 8, `_core.pyi` 6, `cli.py` 3, `edge_strip.py` 3, `final_check.py` 7, `run_record.py` 7 | **163** |
| 20c-2 | `quality.hpp` (cavity, angles, acceptance) ~55; R8's split ~25; option plumbing, CLI and the two rows ~30 | **~110** |
| 20c-2 built | `python3 tools/count_loc.py 41bda81a 65ae0792`: `quality.hpp` 102, `lawson.hpp` 3, `constraint_foot.hpp` 1, `refine.hpp` 6, `bindings/core.cpp` 5, `_core.pyi` 5, `cli.py` 29, `run_record.py` 6. Over the estimate by 47: the plumbing is 51 against ~30, mostly `cli.py`'s option declaration (the help text alone is five lines), its refusals (non-finite, above 10, without `--tolerance`, without `--dem`: pin 8, three checks the estimate did not see) and the wiring into the run and its rows; the speed fix `65ae0792` adds 27 to `quality.hpp` (the atan2 skip, the early stop, the reused buffer), which no estimate had; the green step's `quality.hpp` was 75 against ~80 | **157** (green 130, `count_loc.py 41bda81a 600b2acb`) |
| 20c-2 review fix | `python3 tools/count_loc.py 41bda81a 50544838`: `cli.py` 30 (+1, the reworded help string); the rest as above | **158** |
| 20c-3 (to 2026-10-07) | `feature_input.py` (merge, tolerance, coverage check) ~35; CLI and record ~25; the outline rule ~70 | ~130 |
| 20c-3 (2026-10-08, with the repair) | `feature_input.py`: the land-cover stage (polygons collected per source, clipped to the read region, repaired, merged, simplified; the coverage check gone) ~45; the outline rule (an `STRtree` of the outline's segments, no buffer) ~70; `cli.py`, `run_record.py`, `stats.py` (four flags, of which `--features-repair` is new, their refusals, the input and vertex rows, the timing sub-row) ~35; `pyproject.toml` is not counted | **~150** |
| 20c-3 green | `python3 tools/count_loc.py 483221d2 7334b82b`: `feature_input.py` 260 (the outline rule 188 against ~70; the stage 68, part of it moved from the old loop), `cli.py` 63, `run_record.py` 8. Over the estimate by 181 (+121 %), almost all in the outline rule, whose estimate was not taken from its 146-line prototype (ruling G5); under 700, one PR | **331** |
| 20c-3 rulings T1, T2 | `feature_input.py`: the cut's end test ~2 (T1); `ring`'s index list, the loops, the assembly and its validity guard, less the whole-polygon overlay ~+25 (T2) | **~360** |
| 20c-3 fix round built | `python3 tools/count_loc.py 483221d2 8716f062`: `feature_input.py` 285, `cli.py` 64, `run_record.py` 6 | **355** |
| 20c-3 P1 built | `python3 tools/count_loc.py 483221d2 1e011cec`: `feature_input.py` 286, `cli.py` 64, `run_record.py` 6 | **356** |
| 20c-3 code review round 1 rulings | R1 a keyword argument on an existing line (ruff may wrap it: 0 to 2 lines), R4 a string: 0 | **~357** |

On the worst overrun seen when estimated (+60 %), 240, 175 and 240. Each under 700. 20c-3 then overran by more than that, +121 % to 331 (ruling G5), still under 700.

## Ola's rulings

Ola, 2026-10-06, on round 1's questions 1 to 3: "defaults on the 20c
questions and add the snap-spacing row". So:

1. **Land-cover outlines before meshing.** (a) Borders between two
   polygons of the same class are dropped: yes, on whenever a class map is
   used. (b) Simplifying outlines by a set distance: available as a flag,
   `--features-tolerance`, off (0) by default. 20c-3 no longer waits.
2. **The soft criterion**: on, with gain P = 0.
3. **Quality during refinement too** (increment 20's C3): no, not in 20c.

Ola, 2026-10-06, on question 4 (may the quality start split a land-cover
or outline line when the point it wants lies just beyond that line?):
"First, yes to Q4, then pass it on." So:

4. **The line split (R8)** is in 20c-2, on whenever the soft criterion is.
   The second part of the question (accept Numedalslagen's worst angle
   falling without it) no longer arises.

Ola, 2026-10-07, on the question whether a later step closes gaps between
neighbouring land-cover polygons (such as Lagan's 1 cm slit, M5), asked
with the default no: "defaults on all four". So:

5. **No gap-closing step.** Gaps between neighbouring land-cover polygons
   stay as the source has them; the thin triangle such a gap forces stays
   in the mesh (M5), and the worst-angle gates do not ask more of 20c.
   **Reversed by ruling 9** (2026-10-07): this was the main session's
   default, taken under "defaults on all four", not a call Ola made.

Ola, 2026-10-07, on questions 7 and 6 (the main session's list put
question 7 first and question 6 second): "1-3 default. 4 must wait". So:

6. **The outline rule is on by default** (question 6, the default yes):
   20c-3's `--features-outline-snap` defaults to 5 m whenever features are
   given; `--features-outline-snap 0` turns it off. 20c-3's code, tests
   and gate are as designed; its CLI test pins 5 m.
7. **20c-1's split-phase time is accepted** (question 7, the default
   accept): the foot search stays in the serial split phase, under the
   revised limit ("The split-phase limit, judged after the fix"); no
   follow-up moves it into the parallel scan.

Ola, 2026-10-07, on question 8, asked by the main session after `@perf`'s
re-time `5b0afe04`: "Yes to both" (the other item was an unrelated push).
So:

8. **20c-2's split phase is judged per split** (the default yes), because
   its meshes differ from 20c-1's; recorded in "The split-phase limit for
   20c-2, judged per split", where it is met.

Question 8 as asked: "20c-2's time limit. The step that splits lines is
over the 2 % limit when measured in seconds: +16 % on Lagan and +4 % on
Ljungan. Measured per split, it's faster: −0.8 % and −2.6 %. The new rule
produces a different mesh, with 17 % and 7 % more splits, and that
accounts for the extra seconds. Should we judge it per split, as we did
for 20c-1? Default: yes, and @architect writes that ruling down before the
code review." "Ljungan" there is a slip for Numedalslagen, the catchment
`@perf` measured; the figures are Numedalslagen's. (For 20c-1 the limit
was revised to point 2 against master, not judged per split; the question
stands on its figures.)

Ola, 2026-10-07 (08:19 UTC), on ruling 5: "I don't recall making the
call that we should not fix errors in the input? I consider the cap in the
CORINE-data an error. The edges should have been shared." And, three
minutes later: "yes to the morning check proposals, build 20c-2 and
20c-3". So:

9. **Errors in the input are repaired** (amends ruling 5). Lagan's 1 cm
   slit (M5), where two land-cover borders that should be one end 1 cm
   apart, is an error in the source, and 20c-3 repairs it: borders of one
   land-cover source that come within a tolerance of each other become one
   border, and a gap narrower than it goes to a neighbour (20c-3's step 1,
   `--features-repair`). The tolerance's default is question 9, ruled
   below.

Ola, 2026-10-08 (06:47 UTC), on question 9 as asked again with M8
(`d7a3462a`), the main session's question to him being "Question: set
the default to 1 m? @architect's default: yes.": "yes, set default to
1 m." So:

9 (default). **The repair's tolerance defaults to 1 m** (question 9's
   default as asked again), on whenever a land-cover map is used;
   `--features-repair 0` turns it off, and a smaller value suits a
   land-cover file drawn to the decimetre. It replaces M7's 5 cm. Red
   `f6e0f287` (`@tester`), green `c71ccf57` (`@developer`); the mesh
   figures are M8's 1 m row ("PRs and gates", "With the repair's default
   at 1 m").

Ola, 2026-10-08 (05:54 UTC), on question 10, as the main session put it
to him: the clean-up step adds about 5.5 to 6 s (+17 % on Lagan, +54 %
on Numedalslagen against master; the step itself takes 5.92 s and 5.48 s,
and the whole run grows by 2.95 s and 3.70 s,
`docs/benchmarks/2026-10-08/20c-3/README.md@41db7333:105-106`), with the options "(a) Default: accept
it now and add a ROADMAP row to speed the step up later" or "(b) Hold
20c-3 for another speed round first". Ola: "Let's accept it. One future
option could be to fix the CORINE data itself, as a preprocessing step,
for areas of interest. Goes into the ideas file." So:

10. **The clean-up's time is accepted for 20c-3** (option (a), the
   default): no further speed round before its code review and push. The
   speed-up is the ROADMAP row "Speed up the land-cover clean-up step"
   (T3's hotspots), which also carries Ola's preprocessing idea as an
   option: repair the CORINE data once per area of interest and keep the
   repaired copy, so a run does not repeat the repair. There is no ideas
   file yet (`docs/ideas.md` does not exist), so the row holds the idea
   until there is one.

Question 9 as asked again with M8 (`d7a3462a`):

9. **How close must two land-cover borders be to count as one?** (asked
   again with M8, after your "try increasing to 1m as well".) 20c-3 joins
   borders of one land-cover file that come within this distance of each
   other, and fills any gap narrower than it. Measured on the built code
   with the defaults, Lagan and Numedalslagen (M8):
   - **5 cm** (today's default): the slit closed; triangles with an angle
     under 1°: 180 and 158; the smallest angle in the mesh 0.19° and 0.63°.
   - **25 cm**: the same as 5 cm (176 and 158).
   - **1 m**: the slit closed as well; **26 and 2** such triangles, the
     smallest angle **0.73° and 0.83°**, and 2 % fewer triangles on
     Numedalslagen. Almost all the triangles 1 m removes sat where two
     land-cover borders ran under a metre apart. The price: borders move
     by under 1 m along about 45 km and 123 km of border, some 3 400 m²
     and 8 500 m² of land cover changing class inside the catchments, out
     of 6 441 and 5 548 km². CORINE's smallest feature is 100 m wide, so
     that is 1 % of it. No measurable time.
   **Default: 1 m, on whenever a land-cover map is used** — it is the
   probed distance that removes the thin triangles, and it moves land cover
   by less than CORINE's own precision. `--features-repair 0` turns the
   joining off, and a smaller value suits a land-cover file drawn to the
   decimetre. Not probed: distances between 25 cm and 1 m, or above 1 m.
   The other choice: keep 5 cm. Either way only the default and the tests
   that pin it change.

Questions 6 and 7 as asked:

- Question 6: **Your outline rule: on by default?** 20c-3 builds it either way, as
   `--features-outline-snap METRES`; this question sets only its default.
   Every part of a
   land-cover border within 5 m of the catchment outline is moved onto the
   outline, and a land-cover piece thinner than that along the outline
   goes to its neighbour (both halves of your proposal, M6). Measured on
   20c-3's input (2 m simplification): slivers 151 to 61 on Lagan and 32 to
   4 on Numedalslagen, none left within 20 m of the outline, under 0.25 %
   more triangles; Numedalslagen's worst angle 0.012° to 0.83°, Lagan's
   unchanged (set elsewhere). On the input as given it removes 21 % and 9 %.
   It changes the class of about 0.7 ha per catchment (0.0001 % of the
   area). **Default: yes, on at 5 m whenever features are given** (0 turns
   it off), since it moves borders only within 5 m of the outline, far
   below CORINE's 25 ha minimum mapping unit, and helps with or without
   `--features-tolerance`. The other choice is off unless asked for (0 by
   default), like `--features-tolerance`. Either way the same code, tests
   and gate go into 20c-3; only the default and the CLI test that pins it
   differ.
- Question 7: **20c-1's refinement step with feet on: accept the time, or move the
   foot search into the parallel part?** After `@developer`'s fix, the
   step of refinement that adds points one at a time is 7 % slower with
   feet on than with `--no-constraint-feet` on Lagan and Numedalslagen
   (8 to 9 % per point added), but 3 % and 8 % faster than master as you
   run it today, and the whole run is within 0.3 % of master. On the 1 m
   benchmark, where no point lands near a line, there is no cost at all.
   Moving the search into the part that runs on all cores could save at
   most 2 to 10 ms per run, and adds code. **Default: accept**, under the
   revised limit ("The split-phase limit, judged after the fix"), which
   20c-1 meets. The other choice: a follow-up PR moves the search and ε
   into the parallel scan, with the serial phase checking the answer
   again.

Question 10 as this file put it (the main session asked it in the
shorter form quoted above ruling 10):

10. **The land-cover clean-up is slower than the design said: accept it
   for now?** You asked to hear about any step that takes a large share of a
   run. The new land-cover step (the clip, the repair of question 9, the
   merge by class and your outline rule) takes 8.0 s on Lagan and 8.8 s
   on Numedalslagen. The design estimated 1.1 s (question 9 as first asked), counting only the clip
   and the repair, and those two alone take about 1.9 s. A fix in this
   PR (ruling T2: rebuild and measure only the land cover that moved,
   with the same output) brought the step to 5.92 s and 5.48 s (measured,
   `@perf`, one run each). The whole run is now +17 % on Lagan (17.64 to
   20.59 s) and +54 % on Numedalslagen (6.81 to 10.51 s) against master. On Numedalslagen that
   step is about half the run. The rest is the outline rule (about
   2.4 s), the repair (1.2 s), and clipping the land cover to the
   catchment (1.0 s) and to the read area (0.75 s). **Default: accept
   that time for 20c-3**, and add a ROADMAP row to make the step faster
   later (perhaps 1.5 to 2 s less; not designed or measured). The other
   choice: hold 20c-3 for a further speed round first. Turning the
   outline rule off by default would save about 2.4 s, but it would give
   up most of the rule's gain near the outline; that was not measured
   with the defaults.

## Questions for Ola

None open.

## Not in scope

- Off-lattice Steiner points graded below a DEM cell (option c's Ruppert
  proper): not recommended (M3, "What is left").
- The noder's `MalformedOutput` at `--snap-spacing` 0.1 m and 1 m (M4): a
  defect of its own, not needed by 20c-3; its own unnumbered ROADMAP row,
  approved by Ola 2026-10-06, not designed.
- The vertex rule in the quality start (R2.7): measured to do nothing.
- The first form of Ola's outline rule (thin clipped pieces handed to a
  neighbour by overlay): reaches 4 to 7 % of 20c-3's slivers and makes
  millimetre segments of its own (M6 as of `a3a2b7ec`, round 2; superseded
  by M6's second form, which covers it).
- Closing gaps between land-cover polygons of two different sources, or
  gaps wider than `--features-repair` (none between 5 cm and 100 m wide
  on either catchment, M7). Gaps narrower than it are 20c-3's step 1 (ruling
  9, which reversed ruling 5).
- Quality during refinement (increment 20's C3): ruled out of 20c by Ola
  (question 3).
- A size bound or sizing field (increment 14's U2).

## Review

20c design review round 1 (@reviewer, ed12512..7fcac1e0, 0 counted LOC, docs only): CHANGES REQUESTED. (B1) the 20c-2 gate "share under 1° ≤ 1.15 × 20c-1's" fails on the design's own Lagan numbers (1.162); (B2) no gate on the worst angle, and the claim that the worst angle "does not improve without 20c-3" contradicts the design's own table; (B3) R4 does not mention the final check's existing near-edge rule L12; (B4) R4.3's interpolated foot height was never measured, and the 20c-1 gate's margin is smaller than what the final-check foot contributes.
20c design review round 2 (@reviewer, 7fcac1e0..3d1df11f, 0 counted LOC, docs only): CHANGES REQUESTED; round 1's B1-B4 answered. (B1) the 20c-2 gate "worst angle ≥ 20c-1's" sits exactly on its measured value on Lagan (0.002384° both), against the file's own "a margin for the production code"; (B2) the 20c-3 worst triangle (0.006608°) is called a spike in the outlines themselves, but M3's 1 cm run reached 0.008589°, impossible if those two 79 m sides met as constraints in the 1 cm input, so the 2 m simplification likely made it; reconcile, and fix the reason (two constraint segments meeting at that angle at a shared vertex); (B3) the 20c-3 gate's "19" slivers with a side under 10 cm has no row in M5.
20c design review round 3 (@reviewer, 3d1df11f..a3a2b7ec, 0 counted LOC, docs only): CHANGES REQUESTED; round 2's B1-B3 answered (B2's 1 cm gap confirmed in the CORINE source, M5's 20c-3 Lagan row reproduced exactly). (B1) M6 measures Ola's outline rule on 20c-2's mesh, but the rule was proposed for 20c-3: on 20c-3's mesh 34 of 151 slivers (22.5 %) lie within 5 m of the outline with a CORINE line nearby, and the thin-piece rule reaches 6 and 11 of them (4 % and 7 %), against M6's 8 % and "at most 2.5 %"; (B2) the file says "Questions for Ola: None open" while M6 sets aside Ola's own proposal; (B3) M1's "None is forced by an input angle" contradicts the new finding that the 1 cm gap's thin triangle is in every mesh of that input.
20c design review round 4 (@reviewer, a3a2b7ec..6ac8ab2d, 0 counted LOC, docs only): CHANGES REQUESTED; round 3's B1-B3, S1 and S2 answered (S1's 0.86 × bound recomputed; both catchments' CORINE pass `coverage_is_valid` in the computation CRS, as 20c-3 step 2 says); 20c-1's part unchanged since round 2 and independent of question 6: ready to build. (B1) 20c-3's Lagan worst triangle (0.006697°) is made by the 2 m simplification, not the input: it stands on the slit's own 1 cm closing edge (EPSG:3006 413 602.500 6 331 638.254 / 413 602.501 6 331 638.244), its apex 413 527.395 6 331 662.317 78.9 m away, because `coverage_simplify` drops the source vertex 413 596.806 6 331 640.083 that stands 5.98 m from that edge (0.088° in the source); so docs/increments/20c-soft-quality.md@6ac8ab2d:164-166, :448, :452-453, :455-456, :480-484, :597 and :643-648 and ROADMAP.md@6ac8ab2d:61 ("stays set by the input", "150 m from the slit") are false, and :474-475's "needs no guard" needs re-judging for it; (B2) question 6's other choice (the rule as a flag, off by default) still builds the rule, but the Status, gate, tests and LOC build and gate it only "if Ola says yes" (docs/increments/20c-soft-quality.md@6ac8ab2d:8-9, :611-613, :1028-1030, :1038, :1173, :1225); say question 6 sets only the default, or add a plain no and what it drops.
20c design review round 5 (@reviewer, 6ac8ab2d..51b4434a, 0 counted LOC, docs only): CHANGES REQUESTED; round 4's B1, B2, S1 and S2 answered (rerun on the source: in the designed order (moved to EPSG:3006 first) `coverage_simplify` at 2 m drops 413 596.806 6 331 640.083, which is 4.47 mm off the border and 5.98 m from the 1 cm edge; angles 0.0882° before and 0.006697° after; the gap's corners 413 602.500 6 331 638.254 / 413 602.501 6 331 638.244 now match the source; 27 % and 33 % correct); 20c-1's part unchanged since round 2: ready to build. (B1) docs/increments/20c-soft-quality.md@51b4434a:739-740 still says the outline rule's ~70 lines come "if he says yes (question 6, M6)": round 4's B2, missed in section (b); (B2) docs/increments/20c-soft-quality.md@51b4434a:502-503 "Lagan has 616 such edges and puts back 78 vertices" counts each shared edge and vertex once per polygon: the source has 417 distinct edges under 10 cm, and the guard puts back 39 distinct points (Numedalslagen's 1 779 and 142 not rerun; restate both in one unit).
20c design review round 6 (@reviewer, 51b4434a..c074f900, 0 counted LOC, docs only): APPROVED; round 5's B1, B2, S1 and S2 answered. (B1) /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/increments/20c-soft-quality.md@c074f900:741-742 now says the outline rule is "built either way (question 6 sets only its default)", and a grep for "says yes|if Ola|if he|question 6" finds no conditional build left outside round 4's quoted record. (B2) :502-505 counts each shared edge and point once: rerun on Numedalslagen's box read (8 017 polygons) gives 1 362 distinct short edges / 71 points put back (1 779 / 142 per polygon), matching the doc. (S1) :1045-1048 rerun: 3 polygons fail `coverage_invalid_edges`, none of them reach the domain, and the 1 611 that do reach it pass `coverage_is_valid`. (S2) :691 "about 6 % (0.006697° to 0.007080°)" checks out: 0.007080 / 0.006697 = 1.057. ROADMAP.md@c074f900:61 matches. 20c-1's part has not changed since round 2.
20c-1 code review round 1 (@reviewer, c074f900..0346ac63, 163 counted LOC): CHANGES REQUESTED. (B1) no mutant kill record for the invariant-critical suite test_constraint_foot reached review, and none is on disk: /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/increments/20c-soft-quality.md@0346ac63:17-18 still names "the mutation round" as next; @tester's handback must cover all seven targets at :1255-1259, including CF3's "holding triangle dropped from the active set", which the brief's list left out; (B2) the CF4 oracle-margin mutant demonstration that ruling 2 asks for (:1430-1434, `legalise_around` skipped, depths far above η) is called done at :1503 but its result is recorded nowhere: one line with the measured depths, or the handback, is needed. Checked and true: LOC 163 against an estimate of ~150 (:1670-1671), all under 700; the green, MeshVertex and fix commits touch no test file; no red-step scaffolding is left; F2's re-pin to 2 refusals adds an assertion that both nodes are vertices, and the stub and void-callable edits are type and keyword changes only, so nothing was weakened; on HEAD, built on macOS arm64 Release (which fuses multiply-adds, as the macOS CI leg does), CF4 excuses 2 quads at 3.85e-10 and 4.11e-10 m, about 1 300 times under η = 5e-7 m; the fix 69f37d1c gives the same answers by reading: /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/include/terrain/mesh/constraint_foot.hpp@0346ac63:78-89 tests the same edge set constraint_foot searches, foot_epsilon clamps to the same foot_cap(g) (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/include/terrain/refinement/refine.hpp@0346ac63:219-241), and nothing found within the cap means nothing found within ε; R1-R5, R2.8, R4.7, the gate row and the split-phase table match /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/benchmarks/2026-10-07/20c-1/README.md@b9e2d482:40-50, :206-213, :264-288, :298 and /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/README.md@c9f62f7c:158-170, :186-241, :259, :261-321, and the raw/*_stats.md lines cited; the revised time limit is argued, not a moved goalpost (ruling 5's 2 % was the trigger for the fix, which is done; the tile at +0.86 % shows no work wasted where no foot is placed), and question 7 stays Ola's. Not pushed, so no CI yet.
20c-1 code review round 2 (@reviewer, 0346ac63..74889310, 163 counted LOC c074f900..74889310, tests and docs only since round 1): APPROVED; round 1's B1 and B2 answered. (B1) /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/increments/20c-soft-quality.md@74889310:1590-1600 covers all seven targets named at :1260-1264, including the "foot a hair from an input vertex" case (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/tests/cpp/unit/test_constraint_foot.cpp@4a64c1ef:237) and CF3's holding triangle (line 177). Each cited kill line was spot-checked and asserts what the table says. :1625-1641 records the four survivors as killed, and each new test asserts its premise first. CF1's exact-δ case checks that the distance is exactly 0.5, which holds in doubles on both the fused multiply-add leg and the plain leg (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/tests/cpp/unit/test_constraint_foot.cpp@7bc41ebf:168-170). Each `foot_reachable` false case also checks that `constraint_foot` returns None (:309-331). CF4's two-lines replay uses the same seeds as `refine_points` for a split edge with no neighbour, {owner, new} (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/include/terrain/refinement/refine_points.hpp@0346ac63:369, :423-426). It asserts no flip on the start, a first Hit on A2-B2, at least one flip, and a second Hit on A1-B1 at x = 5 (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/tests/cpp/property/prop_constraint_foot_final.cpp@7bc41ebf:475-499). So the guard at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/include/terrain/refinement/refine_points.hpp@0346ac63:386 is reached. R4.2's companion case asserts that the holding triangle t is unchanged and that no flip writes to it. The two line corrections hold: CF1 :170 and CF4 :514 at 7bc41ebf. (B2) The depths are recorded at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/increments/20c-soft-quality.md@74889310:1656-1675: the shallowest is 0.0196 m against η = 1e-7 m. Ola's rulings 6 and 7 at :1908-1918 quote "1-3 default. 4 must wait", dated 2026-10-07, with 1 = question 7 and 2 = question 6, matching the transcript. Round 1's record at :1984 is word for word the round 1 handback (2 413 characters, identical). `git diff --stat 0346ac63 74889310 -- include src bindings src_python` is empty, and no red-step scaffolding is left. ctest 886 of 886 passed in build-tester, whose binaries are newer than every source file and list the new cases. `pytest --no-cov` gave 5 404 passed and 17 skipped. mypy, ruff check, ruff format, the prohibited-dependency gate, the detria boundary gate and check_citations --base origin/master are all green, and none of the 110 at-risk citations comes from this round's edits. Not pushed, so no CI yet.

20c-2 code review round 1 (@reviewer, 41bda81a..8ea0f32d, 157 counted LOC by `python3 tools/count_loc.py 41bda81a 8ea0f32d`): CHANGES REQUESTED. (B1) The `--start-quality-gain` help says the run will "split a line it lies beyond" (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/src_python/tin_engine/cli.py@8ea0f32d:694-696). That is false with `--no-constraint-feet`: the line split (R8) only runs when both the gain test and the feet are on (pin 4, /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/docs/increments/20c-soft-quality.md@8ea0f32d:1145-1150; `split_lines = judged && o.constraint_feet`, /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/include/terrain/mesh/quality.hpp@8ea0f32d:256). The binding's doc and `_core.pyi` already say "with constraint_feet"; the CLI help should say the same, for example "and, unless --no-constraint-feet, split a line it lies beyond". Checked and true: LOC; ruling 8's figures; the quoted question; the atan2 skip; the early stop; the mutation record; open item (b); ctest and pytest on a HEAD build; gates (details in the list below). Not pushed, so no CI yet. Suggestions: (S1) copy @tester's kill table into this file as "Mutation round for 20c-2"; (S2) say the round ran on 96119508, before the speed fix rewrote `pays`, and that byte-identical meshes carry it over; (S3) ruling 3 and the comment at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/include/terrain/mesh/quality.hpp@8ea0f32d:199-204 say θ at most 35°, but the argument holds for θ up to 180°; (S4) tell Ola about the two slips (done in chat, 2026-10-07).

20c-2 code review round 2 (@reviewer, 8ea0f32d..d4ac395f, 158 counted LOC by `python3 tools/count_loc.py 41bda81a d4ac395f`, 1 since round 1: the help string): APPROVED; round 1's B1 and S1-S3 answered. (B1) The `--start-quality-gain` help now reads "and, unless --no-constraint-feet, split a line it lies beyond" (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/src_python/tin_engine/cli.py@d4ac395f:694-697). That matches `split_lines = judged && o.constraint_feet` (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/include/terrain/mesh/quality.hpp@d4ac395f:257). The red test reads the option's own help (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/tests/python/test_cli_start_quality_gain.py@9f54bc95:230-244). Applied to the help string at 8ea0f32d, its clause check fails; at HEAD the file's 26 tests pass. Red 9f54bc95 touches only that test file and green 50544838 touches no test. (S1) The 25 table rows at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/docs/increments/20c-soft-quality.md@d4ac395f:2104-2178 are byte-identical to @tester's hand-back (session transcript 6a989a29, line 23672). (S2) :2123-2135 is true: at 96119508, `pays` made one comparison and `GainJudge::capped` did not exist. (S3) The comment at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/include/terrain/mesh/quality.hpp@d4ac395f:199-205 and ruling 3 at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/docs/increments/20c-soft-quality.md@d4ac395f:2072-2078 now give θ ≤ 180°. Nothing in the C++ or the binding sets an upper bound on θ, and `improve` returns early unless θ > 0, so a negative θ never reaches the skip. `git diff 8ea0f32d d4ac395f -- include src bindings` shows only that comment, so round 1's C++ build and ctest still stand. The LOC row (158 = 157 + 1) and the status lines :25-46 match. No red-step scaffolding is left. mypy, ruff check, ruff format, the prohibited-dependency gate, the detria boundary gate and `check_citations --base origin/master` (exit 0) are green, and none of the 67 at-risk citations comes from this round's edits. Not pushed, so no CI yet.

20c-3 code review round 1 (@reviewer, 483221d2..1e011cec, 356 counted LOC by `python3 tools/count_loc.py 483221d2 1e011cec`): CHANGES REQUESTED. (B1) A hole the outline rule moves onto the outline gets filled: `shapely.make_valid(shell)` at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/src_python/tin_engine/feature_input.py@1e011cec:663 uses the default "linework" method, so a part whose hole now touches its shell keeps the hole's area. With the merge on, a lake inside a thin ring of a rarer class then takes the ring's label, because the smallest polygon wins. Probe: 1 km square, D 5 m; water = lakes box(3,200,103,300) and box(500,500,800,700); ring = box(2,190,113,310) less the first lake. After the rule the ring's polygon is 13 560 m² and covers the lake centre, so the ring's class wins over water (70 300 m²). The input labels it water. The bug is there since the green step: e07d5921 and 35b7a873 give the same polygons. `make_valid(..., method="structure")` gives 70 300 / 3 260 / 926 440 m², a partition summing to 1e6, with area_changed unchanged at 340 m². Neither catchment hits it: a patched run found no filled hole on Lagan or Numedalslagen. Either fix it with this fixture as a test, or have @architect rule it a known limit next to P4. (B2) Prose the branch made untrue: "PR #214, in the merge queue" at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/docs/increments/20c-soft-quality.md@1e011cec:46 and /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/ROADMAP.md@1e011cec:67. #214 is merged as 4abacf65, inside this branch's own base. "Next: @tester red for P1, @developer green" at :78-80 and at ROADMAP :67: both are done (48aca5f7, 1e011cec). The LOC table at :3060 says 355 plus "about +3"; the built count is 356. "before 30d merges" at :1317-1318 and "open when this was written ... If 20c-3 merges before 30d" at :1506-1511: 30d is merged (769ed992, in the base). (B3) The 20c-3 gate row at :1525 asks for four things that were neither run nor ruled: the shared-border check, the population-3 check, the slit's far-corner check, and @architect raising the thresholds to the measured figures. @perf records the three checks as not run (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/docs/benchmarks/2026-10-08/20c-3/README.md@41db7333:79-89). The "+ repair 5 cm" runs predate T1. Yet :76 says "every outline-rule gate met". Either run them or have @architect rule them out. Checked and true: 356 lines under 700, size explained, mutation record covers all targets plus T1's end test, P1 reproduced, default meshes equal @perf's T2 (ef67e65c, 1f11f291), full pytest 6 878 passed, real-data files 171 passed from the main checkout, all static and governance gates green, merge-tree clean; not pushed, so no CI yet. Suggestions: (S1) say the mutation round ran on e07d5921, before T1 and T2 rewrote `_cut` and `ring`; (S2) the record row "Land-cover area inside the outline that changed class" counts only what the outline rule moved; "...that the outline rule gave another polygon" would be exact; (S3) P1's four-decimal check on the real inputs is not recorded in any file.

20c-3 code review round 2 (@reviewer, 05b0fd7a..07d97ee4, 360 counted LOC by `python3 tools/count_loc.py 483221d2 aaee898b`; the only production change since round 1 is aaee898b's +4 lines): APPROVED; round 1's B1-B3 and S1-S3 answered. (B1) `shapely.make_valid(shell, method="structure")` at /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/src_python/tin_engine/feature_input.py@aaee898b:663-665. With aaee898b's two code changes reverted in memory and 1f7f8335's tests applied, four tests fail (/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/tests/python/test_feature_outline_rule.py@1f7f8335:593-663, /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/tests/python/test_cli_features_cleanup.py@1f7f8335:175-187); on aaee898b both files pass, 54 of 54. Default meshes rerun: Lagan `ef67e65c…`, Numedalslagen `1f11f291…`, equal to @perf's /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/docs/benchmarks/2026-10-08/20c-3/stats/r3/sha256.txt@63065a15. (B2) /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/docs/increments/20c-soft-quality.md@07d97ee4:46 and /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/ROADMAP.md@07d97ee4:67 say #214 is merged; the "Next" lines name this round; :1519-1523 says 30d is merged; no stale phrase remains outside round 1's quoted record. (B3) :3081-3152 rules (a)-(e) and rules out raising the thresholds; :1561-1573 and the gate row at :1537 match /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/docs/benchmarks/2026-10-08/20c-3/README.md@63065a15:132-183; (e)'s corrected wording checked against the CORINE source (312 against 324 and 313 at each stretch). (S1) :3159-3165 true. (S2) R4 matches /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/src_python/tin_engine/run_record.py@aaee898b:59-61. (S3) R6 agrees with @perf's stats rows. LOC 360 against ~357, under 700. No red-step scaffolding. Full pytest 6 885 passed, 17 skipped; real-data files 417 + 26 passed from the main checkout with this venv; mypy, ruff, prohibited-dependency, detria boundary and check_citations green; merge-tree clean against origin/master 483221d2. Not pushed, so no CI yet. Suggestions: add a "built" 360 row to the LOC table; point @perf's README :178-183 to R3 (e)'s outcome; script the source lookup behind (e).
