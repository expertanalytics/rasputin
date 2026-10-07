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
**After approval, on Ola's yes:** the push, the PR and the enqueue, then a
check that `CI result` is green; whether these have happened is shown by
`gh pr list --head worktree-soft-quality --state all` and `gh pr checks
<pr>`, not by this file;
**20c-2** the soft criterion and a split of the
constraint line that blocks the walk to a quality point, after 20c-1
(branch `worktree-soft-quality-2` off master `41bda81a`; Ola: "build
20c-2 and 20c-3"): red tests `100f5a09` (`@tester`), their pins ruled
(below, "Pins ruled for 20c-2's red step"; R7 gains a slack and R8 needs
the feet on), next `@tester`'s amendment for pins 3 and 4, then
`@developer` (green);
**20c-3** (still to build) input clean-up and coarsening, and Ola's outline rule, on by
default at 5 m (question 6, ruling 6), after
20c-2. Questions 1 to 4 ruled by Ola, 2026-10-06,
and questions 5 to 7, 2026-10-07 ("Ola's rulings" below); none open. Asked by Ola
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
to 61 and 4, with Numedalslagen's worst angle at 0.83° (M6).

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
GEOS 3.14; the venv's shapely has GEOS 3.13), and Ola ruled no
gap-closing step (ruling 5, 2026-10-07). It is one of M1's two slivers
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
0.007080°, which no step can pass under ruling 5, removes no sliver on
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
0.007080° on any mesh, so no step short of closing it (ruled out, ruling
5) gets Lagan above about 0.007°, and keeping that vertex was measured to
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

## Design of PR 20c-3: input coarsening (ruled by Ola, question 1), and the outline rule (ruled by Ola, question 6)

Python, in the feature reader, after the polygons are moved to the
computation CRS (vertex by vertex, as today) and before they are clipped
and the chains are built; never in `_core` (the I/O boundary, `CLAUDE.md`
§2). Three independent steps, in this order, each with its own flag:

1. **Same-class borders dropped** (`--features-merge-same-class`, on
   whenever a class map is used, Ola's ruling; `--no-features-merge-same-class`
   turns it off): under a class map, neighbouring polygons of one class
   are unioned (shapely `coverage_union` per class), so a border with no
   change of class is no constraint. Removes population 3's lines (E12).
2. **A horizontal tolerance** (`--features-tolerance METRES`, default 0,
   off, Ola's ruling): the polygons that reach the domain through
   `coverage_simplify` at that tolerance (`simplify_boundary=False`, as
   measured), with no 1 cm snap first (round 3:
   it breaks the coverage, M6 (a)); polylines are left as given (a line
   network is not a coverage, and nothing here measured them). Before
   simplifying, `coverage_is_valid` over those polygons (the filter to
   the domain is needed: Numedalslagen's CORINE as read for the domain's
   box, 8 017 polygons, fails the check over 3 polygons that do not reach
   the domain, and the 1 611 that do pass it); if they do not
   form a coverage, the run stops with a plain error that says how many
   polygons and how many metres of border do not match, and that
   `--features-tolerance 0` reads them as given (both test catchments
   pass this check). The record and `--stats` say the tolerance and the
   vertex counts before and after.
3. **Ola's outline rule** (`--features-outline-snap METRES`; on by
   default at 5 m whenever features are given, Ola's ruling 6; 0 turns it
   off): M6 (c)'s four
   steps on the polygons' rings, against every ring of the domain polygon
   as the mesh gets it (after increment 22's reduction, if any), in the
   computation CRS. M6 (c)'s rounding (its step 2) then rounds once more: a
   multiple of D within D/2 of an outline vertex goes to the vertex (the
   prototype did not), so two distinct points the rule puts on the
   outline are more than D/2 apart unless both are outline vertices. The
   result is linework: rings cut where they leave the outline, the
   stretches on it dropped, each piece keeping its polygon's class. The
   record and `--stats` say D, the borders moved and the land-cover area
   that changed class. D's scale: metres in the computation CRS, measured
   at 5 and 10 m on DEM cells of 10 and 31 m (M6 (d)); about half a cell
   or less is what was measured.

The domain outline itself is not touched (increment 22 has its own
reduction). About 60 lines for steps 1 and 2 and about 70 for step 3.
Measured in M3 to M6 on Lagan and Numedalslagen (the Norwegian CORINE
rewritten as GeoJSON for round 3, M6). The merge changes the default mesh
of every run with a class map, and step 3 changes every run with
features, since its default is on (ruling 6), so 20c-3's off switch for
the gate is all three flags off, `--features-outline-snap 0` among them.

## PRs and gates

| PR | what | needs | gate (measured value in brackets, M5) |
|---|---|---|---|
| 20c-1 | R1 to R6: the foot rule on the quality start, refinement's neighbours and the final check | nothing | Lagan: share under 1° ≤ 0.09 % (0.078 %), triangles ≤ +1 % of master's (+0.43 %), worst angle ≥ master's 0.000412° (0.0024°); Numedalslagen: share ≤ 0.04 % (0.034 %), triangles ≤ +1 % (+0.27 %), worst angle ≥ master's 0.000399° (0.00068°); `--no-constraint-feet` bit-identical to master run with `--no-constraint-feet` (master's default has 20b's feet on, so it is not the comparison); tolerance oracle; determinism. **Built: every gate passes, at M5's figures** (`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/README.md@c9f62f7c:186-201`); time: the split-phase limit below |
| 20c-2 | R7 and R8: the soft criterion and the line split | 20c-1 merged | on both catchments: triangles ≤ 0.95 × master's (0.925, 0.886); sliver count ≤ 0.75 × 20c-1's (0.64, 0.69); worst angle ≥ 0.95 × 20c-1's (Lagan 1.000, the same triangle; Numedalslagen 12.7 ×); `--start-quality-gain -1` bit-identical to 20c-1 |
| 20c-3 | input coarsening, and the outline rule (on by default at 5 m, ruling 6) | 20c-2 merged | measured without the merge, the rule's flag given explicitly, so the gate does not depend on the default. **`--features-tolerance 2` alone** (`--features-outline-snap 0`): Lagan: slivers with a side under 10 cm ≤ 25 (19, from 185), share under 1° ≤ 0.03 % (0.020 %), worst angle ≥ 0.005° (0.006697°); Numedalslagen: share ≤ 0.006 % (0.003 %), worst angle ≥ 0.008° (0.012446°). **With the outline rule at 5 m on top**: Lagan: sliver count ≤ 90 (61), at most 5 with the centre within 20 m of the outline (0; 83 without the rule), triangles ≤ +1 % of the tolerance-only mesh (+0.09 %), worst angle ≥ 0.95 × the tolerance-only figure (1.000, the same triangle); Numedalslagen: sliver count ≤ 10 (4), at most 5 within 20 m of the outline (0; 27 without), triangles ≤ +1 % (+0.21 %), worst angle ≥ 0.1° (0.832°); no shared border inside the catchment left unmatched by the rule. With the merge on, the two population-3 lines carry no constraint edge (M2's bands). All three flags off give 20c-2's mesh |

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
20c-2, which touches the same paths):

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

### 20c-3

Python only: a two-polygon coverage of one class merges to one; the
partition stays valid and its area changes by less than tolerance ×
perimeter; polygons that do not form a coverage stop the run with the
plain error, and pass at tolerance 0; flags off are a no-op.

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
  crossing, within D/2 of where it was; an inlet (the outline's way between
  two moved points longer than twice their distance plus 2 D) is joined
  straight, not round the inlet.
- **OR6** D = 0 is a no-op, byte for byte.
- Mutants to kill: the cut not in a fixed order (OR2); the vertex dropped
  at a rounded point (OR4); the stretches on the outline kept (OR1); the
  rounding removed (OR4); the inlet test removed (OR5).

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
after the fix"). 20c-3 is Python only and needs no `bench.py` run.

**20c-1, done.** First acceptance at `34b546c4`
(`docs/benchmarks/2026-10-07/20c-1/README.md@b9e2d482`): gates pass, refine
time regressed (+6 to +24 % on the 1 m benchmark). Re-time after the fix
at `69f37d1c`
(`docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/README.md@c9f62f7c`): the 1 m
benchmark accepted, −0.3 to −2.0 %, meshes byte-identical to master's;
gates pass at the same figures.

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
| 20c-3 | `feature_input.py` (merge, tolerance, coverage check) ~35; CLI and record ~25; the outline rule ~70 | **~130** |

On the worst overrun seen (+60 %), 240, 175 and 210. Each under 700.

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

The two questions as asked:

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
- Closing gaps between land-cover polygons (Lagan's 1 cm slit, M5): ruled
  out by Ola (ruling 5, 2026-10-07).
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
