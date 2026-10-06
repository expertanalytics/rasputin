# Increment 20c — fewer slivers: constraint-aware insertion everywhere, and a soft quality criterion

Status: **designed** (`@architect`, 2026-10-06, branch
`worktree-soft-quality` off master `ed12512`), in three PRs: **20c-1** the
foot rule on every insertion path, ready for `@tester`; **20c-2** the soft
criterion, after 20c-1; **20c-3** input coarsening, waiting on Ola
(question 1). Asked by Ola ("2: yes", 2026-10-06, after calling Lagan's
worst angle of 0.000412° "pretty unacceptable!"). Carries increment 20's
C1-C3 rulings (`docs/increments/20-start-quality.md`, "Ola's rulings").

**In one paragraph.** Most slivers are not made by refinement but by the
two insertion paths that lack 20b's rule (put a point that would land a
hair from a constraint line onto the line instead): the quality start and
the final check against the source DEM. Giving them the rule, and letting
20b's rule look one triangle further, removes 56 to 71 % of the slivers
for +0.2 to +0.4 % triangles (Lagan 1 598 to 709, Numedalslagen 1 644 to 481).
A soft criterion then cuts 12 to 17 % of the triangles at the same sliver
count. What is left is the input's own geometry, which only a tolerance on
the input removes, and that is Ola's call.

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
greedy penalty on angle deficit. **No novelty is claimed.** Searched:
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

| apex of the sliver (slivers with a side under 10 cm only in the last row) | slivers | of them, long side a constraint edge |
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
- **None is forced by an input angle**: only 2 of 1 598 have their smallest
  angle between two constraint edges meeting at a vertex (the one kind of
  sliver no Steiner point can remove). Everything else is avoidable in
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
- **Worst angle is the input's.** With the input as given, the worst
  triangle after A + R stands on a 3 mm CORINE segment (two input vertices
  3 mm apart; the third corner a node 75 m away): 0.0024°. No insertion rule
  avoids a triangle on a 3 mm constraint short of grading the mesh down to
  millimetres. Only input coarsening moves the worst angle (0.0086° at 1 cm,
  0.018° at 1 m).
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
  (`include/terrain/noding/noded_pslg_builder.hpp:304`). A defect in the
  noder at coarse spacings, outside 20c; reported, not investigated.
- **In Python, before the noder, it works**: shapely's `set_precision`
  (snap rounding through GEOS) keeps all 7 319 CORINE features valid at 1 cm,
  10 cm and 1 m. At 1 cm it removes population 1 (140 to 7) and the worst
  angle rises from 0.0024° to 0.0086°; at 1 m to 0.018°, with 448 slivers.
- **Coverage simplification** (`coverage_simplify`, GEOS 3.12's
  coverage-preserving Visvalingam-Whyatt, shared borders simplified once so
  the partition stays a partition) at 2 m removes 80 132 of 1 744 066
  vertices and takes the slivers to 305, with 7 % fewer triangles. It is the
  only run that touches the input-vertex apexes (494 to 266).
- **E12's merge** (same class on both sides) is a semantic rule, not a
  tolerance: a CORINE border between two polygons of one class marks no
  change of land cover. It removes population 3 by removing its lines; A + R
  already takes population 3 from 85 to 21.

## What Ola gets

Measured with `--stats`' own quality table (share under 1°, worst angle) and
the triangle count, on two catchments, with the height tolerance kept
("every DEM node within the tolerance", checked by `@tester`'s oracle) and
the output independent of `threads`:

| after | Lagan, slivers (share) | Lagan, triangles | Numedalslagen, slivers (share) | Numedalslagen, triangles |
|---|---|---|---|---|
| master today | 1 598 (0.185 %) | 863 897 | 1 644 (0.128 %) | 1 287 334 |
| PR 20c-1, constraint-aware insertion | about 710 (0.082 %) | +0.4 % | about 480 (0.037 %) | +0.2 % |
| PR 20c-2, the soft criterion | about 720 (0.095 %) | −12 % | about 450 (0.042 %) | −17 % |
| PR 20c-3, input coarsening at 2 m (if Ola says yes) | about 310 (0.044 %) | −18 % | not measured | |

**Against the write-up's two input-side baselines** (main session's
default: 20c's options must beat them). E12 alone (the forest cuts merged
away) leaves 1 515 slivers; 20c-1 leaves 709. E11 (every CORINE segment cut
to 100 m or less, 1 cm grid) leaves 461, fewer than 20c-1, but at 1 073 613
triangles, 24 % more than master; 20c-1 + 20c-2 are at 758 308 (29 % fewer
than E11) with 720, and 20c-3 at 2 m reaches 312 at 704 966. So 20c-1 and
20c-2 beat E11 per triangle but not in count; 20c-3 beats it in both.

The **worst angle** does not improve without 20c-3: with the input as given
it is set by CORINE's own millimetre segments (M3). With 20c-3 it rises
from 0.0004° to about 0.007°. A worst-angle promise is therefore made only
for 20c-3, and only relative to the input tolerance Ola chooses.

What stays: slivers whose apex is an input vertex about a metre from
another part of the same outline (Lagan: about 500 after 20c-1). Removing
them needs either Steiner points finer than the DEM cell (option c below,
not recommended) or the input tolerance of 20c-3.

**São Francisco.** The rates are what carries over: about 0.04 to 0.08 %
of triangles after 20c-1 instead of 0.13 to 0.19 %, and 12 to 17 % fewer
triangles after 20c-2, which matters more there than anywhere. Its land
cover will be rasputin's own polygons from a raster (MapBiomas), whose
vertices rasputin places; 20c-3's tolerance is then part of that conversion
rather than a change to someone else's data.

## The options, with costs

### (a) The foot rule on every insertion path. Recommended, first, own PR.

20b's rule ("a candidate too close to a constraint segment goes onto the
segment, at its foot") on the two paths that lack it, the quality start and
the final check, plus 20b's own blind spot (it looks only at the triangle
holding the node). M3: slivers 1 598 to 709 on Lagan and 1 644 to 481 on
Numedalslagen, caps (the slope-spoiling kind) 1 191 to 311 and 1 323 to 174,
for +0.2 to +0.4 % triangles. About 150 lines of C++ and plumbing. No ruling
needed: it is Ola's C1 lean ("allow adding point to the constraints if that
improves the mesh props") applied where it was missing, and 20b's guarantee
argument carries over unchanged. **Worth shipping first on its own**: it is
the largest gain, the cheapest, and 20c-2's measurements stand on it.

### (b) Input coarsening. Ola's ruling first, then its own PR.

Merging CORINE's millimetre vertex pairs (population 1) and, further,
simplifying outlines within a horizontal tolerance. M4: at 1 cm, population
1 goes (140 to 7 slivers) and the worst angle rises twentyfold; coverage
simplification at 2 m on top takes the total to about 305 with 7 % fewer
triangles. And E12's same-class merge removes population 3's lines.

**Against the input model.** Ola's model (16b R9, "vertices are used as
given") and "tolerances set resolution" pull in opposite directions here: a
horizontal tolerance is a resolution set by a tolerance, but it moves the
input's vertices. So it needs a ruling (question 1). Where it lives: Python,
before the noder (`feature_input.py`'s side of the I/O boundary), with
shapely's `set_precision` and `coverage_simplify`, which keep a partition a
partition. Not the noder's `--snap-spacing`: it fails above 1 cm today (M4).
About 60 lines of Python. **Worth its own PR**, independent of (a), so it
can go in parallel once ruled.

### (c) The soft criterion, and splitting constraints in general

Increment 20's rulings sent three things here: a soft criterion instead of
the hard 25° (C2), points on constraints when they improve the mesh (C1),
and quality during refinement, not only at the start (C3). Measured:

- **The soft criterion (penalty) is worth having, for size.** Accepting a
  quality candidate only when it does not lower the worst angle of the
  triangles it replaces gives 12.6 % (Lagan) and 16.8 % (Numedalslagen)
  fewer triangles at the same sliver count (M3, "pen"). That is Ola's
  "24.99 should be ok, if not exploited": a node that leaves its
  neighbourhood no better is not paid for. About 70 lines. Its own PR,
  after (a).
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
  fix; a choice about the angle distribution (question 3).

## Design of PR 20c-1: the foot rule on every insertion path

### R1. One geometric helper, three policies

- **New header `include/terrain/mesh/constraint_foot.hpp`**, depending on
  `core`, `predicates` and `lattice_mesh.hpp` only; no raster, no height.

  ```cpp
  struct FootHit { std::uint32_t owner; unsigned edge; MeshVertex at; };
  // The first constrained, non-frozen edge closer than delta (world
  // distance in frame f) to p: t's edges in edge order, then, for each of
  // t's unconstrained edges in edge order, the neighbour's other two edges.
  // nullopt if none, or if the first such edge's foot lies within delta of
  // either end of it, or if the split's children are not all strictly
  // counter-clockwise (20b's foot_fits, moved here).
  std::optional<FootHit> constraint_foot(const LatticeMesh&, std::uint32_t t,
                                         MeshVertex p, double delta, const LatticeFrame&);
  ```

  20b's `detail::foot_of` and `detail::foot_fits` (`refine.hpp`) become this;
  refinement keeps only its ε (`foot_epsilon`, 20b R3).
- **Not across a constrained edge.** A neighbour behind a constraint is on
  the other side of a line the candidate is not near enough to matter for.
- **Determinism**: a fixed search order and "first found" make the hit a
  function of the mesh and p only.

### R2. The quality start (`quality.hpp`, `improve`)

After the node is snapped, judged valid and located, and before insertion:

1. **Near a constraint**: `constraint_foot(m, t, node, δ_q, f)`. A hit is
   inserted with `split_edge(owner, edge, at)` instead of the node
   (`QualityOutcome::feet`), if the validity callable accepts the foot.
   A near edge with no usable foot (the foot near an end, a child not
   counter-clockwise, or refused by the callable) is a skip
   (`skipped_near_line`).
2. **δ_q = min(dx, dy) / 2**, the cap of 20b R3. The pass reads no height, so
   there is no slope term. Scale: relative to the cell; checked at 10 m
   (Numedalslagen, 1.29 M triangles) and 31 m (Lagan, 0.87 M).
3. **The validity callable takes a `MeshVertex`** (it took a
   `LatticeVertex`): `refine` passes "`vertex_z` has a value", which for a
   node is today's "not NoData", so nodes are judged as before.
4. **Why the bad triangle still goes.** The node is within half a cell
   diagonal of the circumcentre and the circumradius is at least a cell
   diagonal (20 R5), so the foot, at most half a cell from the node, is
   strictly inside the circumcircle. Legalisation then removes the bad
   triangle as it does for the node.
5. **Termination.** Nodes: as 20 R10. Feet: each is at least δ_q from both
   ends of the sub-edge it splits, so a segment of length L holds at most
   L / δ_q of them (20b R6's second argument).
6. **Off with the feet.** `QualityOptions::constraint_feet` (false by
   default); `refine` passes `RefineOptions::constraint_feet`, so
   `--no-constraint-feet` turns all three paths off and the mesh is
   bit-identical to master's.
7. **No vertex rule.** Skipping a node closer than δ_q to an off-node
   corner of its triangle was prototyped (`QVERT`) and changed no triangle
   on Lagan. Not designed in (as 20b's C3).

### R3. Refinement (`refine.hpp`): 20b R1 widened

`foot_of` searches through `constraint_foot`, so the holding triangle's
three neighbours count too (M3 "R"). A hit in a neighbour splits the
neighbour's edge; the deferral rule of 14 R5 applies to the owner and to
the triangle across the split edge (either touched this round: the node
waits). Everything else in 20b (ε, footed once, the fallback, the refusal
counts) is unchanged.

### R4. The final check (`refine_points.hpp`, `point_loop`)

For a winner that is a stored source point or a DEM node (not a strip
point, not a void carve point), Inside its triangle, not yet footed:

1. `constraint_foot(m, t, point, δ_p, frame)`; with a strip, the split must
   also pass `strip_fits` (the guard every split of a strip sub-edge
   already passes, 15f L1), else `foot_fits` alone.
2. A hit is inserted instead of the point, with `split_edge`, and the
   strip's sub-edges are cut as for any split of a constrained edge
   (15f step 4). The point is recorded as footed (by its exact position)
   and stays a candidate: if its error is still above the tolerance it is
   inserted later as itself, which is 20b R5's fallback. **The tolerance
   guarantee is the stored points' and is unchanged**: every stored point is
   still scanned, and the loop stops only when every error is within it.
3. **z of the foot**: on the projected path (`refine_strip`, a raster at
   hand), `vertex_z`. On the reprojected path (no raster), linear along the
   line between the two strip points that bracket the foot on its sub-edge
   (they are the source's own values on that line, 15f); with none on one
   side, the sub-edge's end vertex's z on that side. The prototype used the
   point's own z; the design's choice keeps the line's profile, which the
   strip's own check (15f step 6) measures.
4. **δ_p = min(dx, dy) / 2 of the run's frame** (the target grid). Scale:
   relative to the cell; checked at 31 m (Lagan) and 10 m (Numedalslagen).
   Fixed rather than slope-scaled: the reprojected path has no raster to
   take a slope from. 20b measured that a fixed ε fires the fallback more
   often (5 in 150); the fallback count is reported (R5) and is part of
   `@perf`'s acceptance.
5. **Counters**: `PointRefineOutcome::feet`, `feet_fallback` (footed points
   later inserted as themselves).
6. `PointRefineOptions::constraint_feet`, false by default; the CLI passes
   the same switch as `refine`.

### R5. Report and options

- No new flag. `--no-constraint-feet` covers all three paths;
  `elevation_source` keeps `constraint feet on/off`.
- `RefineOutcome::quality_feet` (feet from the quality start) and the final
  check's `feet` and `feet_fallback` reach `--stats`' Refinement table and
  the run record (`points_snapped_to_lines` gains the two other paths'
  counts as separate rows; plain words, per "Plain product output").
- The quality start's skips keep one total, `quality_skipped`, as 20 R11.

### R6. Interactions

- **Frozen edges (23b)**: `constraint_foot` never returns a frozen edge,
  as 20b's search does not today (N5).
- **NoData**: a foot whose `vertex_z` is empty is refused on every path.
- **Strip (15f)**: a foot on a strip sub-edge is a split of it, handled by
  the existing `cut`; strip points are not footed (they are on the line).
- **Determinism**: all three rules run in serial phases, in the existing
  order, and read only the mesh, the candidate and (for z) the DEM or the
  strip. Output bit-identical for any `threads`.

## Design of PR 20c-2: the soft criterion in the quality start

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
  the new ones, both capped at θ. **Insert only if `new ≥ old + P`**;
  otherwise skip (`skipped_no_gain`).
- **P = 0° by default**: a candidate is refused only if it would leave its
  neighbourhood's worst angle lower than before. M3: 12.6 % and 16.8 %
  fewer triangles at the same sliver count. P = 2° was measured too and is
  far too strict (8 437 slivers on Lagan). The objective matters: the summed
  angle deficit instead of the worst angle refused the insertions that
  break up caps (1 393 slivers at P = 0).
- **Why the cavity is exact.** The mesh is constrained Delaunay before each
  insertion (14b R10), so Lawson legalisation after inserting a point gives
  exactly the constrained Bowyer-Watson cavity's triangulation. The read-only
  prediction and what `legalise_around` produces are the same triangle set;
  `@tester` asserts it (T-P1).
- **Termination**: the rule only refuses; 20 R10 and R2 above hold.
- `QualityOptions::min_gain_deg` (P), negative = off (today's hard rule);
  `RefineOptions` and the binding carry it; the CLI passes 0.
- **CLI**: `--start-quality-gain DEG`, default 0; `-1` restores the hard
  rule (bit-identical to 20c-1). Refused: non-finite, above 10.
  `elevation_source` adds `start quality gain 0 deg` (ASCII, 20 R11).
- Reported: `start_quality_points_skipped` keeps its total; a new row
  "Tries that would not have improved the angles" in `--stats`.

## Design of PR 20c-3: input coarsening (only if Ola says yes)

Python, in the feature reader, after reprojection to the computation CRS and
before the chains are built; never in `_core` (the I/O boundary, `CLAUDE.md`
§2). Two independent steps, each with its own flag:

1. **Same-class borders dropped** (`--features-merge-same-class`, default
   per question 1): under a class map, neighbouring polygons of one class
   are unioned (shapely `coverage_union` per class), so a border with no
   change of class is no constraint. Removes population 3's lines (E12).
2. **A horizontal tolerance** (`--features-tolerance METRES`, default per
   question 1): polygons of a coverage through `coverage_simplify` at that
   tolerance after a `set_precision` to 1 cm; polylines through
   `set_precision` only (a line network is not a coverage). The record and
   `--stats` say the tolerance and the vertex counts before and after.

The domain outline is not touched (increment 22 has its own reduction).
About 60 lines. Measured in M3 and M4 (Lagan only; the Norwegian CORINE is
a GeoPackage that the scratch scripts did not rewrite).

## PRs and gates

| PR | what | needs | gate (with the input as given) |
|---|---|---|---|
| 20c-1 | R1 to R6: the foot rule on the quality start, refinement's neighbours and the final check | nothing | Lagan: share under 1° ≤ 0.09 %, triangles ≤ +1 % of master's; Numedalslagen: ≤ 0.045 %, ≤ +1 %; `--no-constraint-feet` bit-identical to master; tolerance oracle; determinism |
| 20c-2 | R7: the soft criterion | 20c-1 merged | triangles ≤ 0.90 × 20c-1's on both; share under 1° ≤ 1.15 × 20c-1's; `--start-quality-gain -1` bit-identical to 20c-1 |
| 20c-3 | input coarsening | Ola's answer to question 1 | Lagan at the chosen tolerance: no sliver with a side under 10 cm, worst angle ≥ 0.005°; flags off give 20c-2's mesh |

The share is read from `--stats`' quality table ("< 1°") and the triangle
count from the same file; no new tool. `@perf` runs both catchments for
each PR on AC power and keeps the `_stats.md` files under
`docs/benchmarks/<date>/`, with the commands. The thresholds come from M3
with a margin for the production code differing from the prototype (z of
the final-check foot, R4.3).

## Tests `@tester` writes red first

### 20c-1

**Invariant-critical (mutation testing required): `test_constraint_foot`**
(new, C++).

- CF1 **the helper**: on hand-built `LatticeMesh` fixtures, a point at
  0.01, 0.4 and 0.6 cells from a constrained edge of its triangle, and from
  a neighbour's: hit and owner as R1 orders them; none at 0.6 (δ = 0.5);
  none when the foot is within δ of an end; none across a constrained edge;
  none on a frozen edge; the foot is the orthogonal projection (to 1e-12
  cells), not the midpoint.
- CF2 **quality start**: a bad triangle whose snapped node lies 0.01 cells
  from a long constrained edge: the output has the foot on the edge and not
  the node; the bad triangle is gone; the constraint edges as a set of lines
  are unchanged; bits and masks on both halves (20 Q3's check). Same with
  the constrained edge on the neighbour.
- CF3 **refinement through a neighbour**: a worst node whose holding
  triangle has no constrained edge but a neighbour does, within ε: the foot
  goes on the neighbour's edge; with that neighbour touched this round, the
  node waits a round.
- CF4 **final check**: a stored point 0.02 cells from a constrained edge
  with error above the tolerance: the foot goes in, z as R4.3 (a strip whose
  profile differs from the point's z tells the two apart); at
  `tolerance == 0` the point goes in after its foot (the fallback), the run
  ends, feet ≤ footed points; the tolerance oracle over all stored points
  passes.
- Mutants to kill: the end check removed (CF1, and a foot a hair from an
  input vertex); the neighbour search crossing a constrained edge (CF1);
  midpoint instead of foot (CF1); footed-once removed (CF4's bound); the
  fallback removed (CF4's oracle); the foot's z from the point (CF4).

Property and integration:

- **Off switch**: `constraint_feet = false` on all three paths gives
  output bit-identical to master on 14b's T12 fixtures, a domain start and
  a features start.
- 14b's T3 (tolerance) and T6 (determinism, threads 1, 2, 7 and hardware
  concurrency) re-run with the rule on; T6 is the TSan job's.
- CLI: the new rows in `--stats` and the record; `--no-constraint-feet`.
- Existing pins that change because the quality start and the final check
  now foot (`test_mesh_quality` counts, golden hashes over scenes with
  constraints): each listed with the reason and amended in its own commit,
  as the NoData fix did (20, Q-V3).

### 20c-2

Not mutation-critical beyond T-P1:

- **T-P1 (invariant)**: for every accepted candidate in a fixture run, the
  predicted new triangles equal the triangles `legalise_around` wrote, as
  sets of vertex triples. Mutants: the cavity crossing a constrained edge;
  the inexact incircle; the foot's second seed missing.
- **T-P2**: no accepted candidate lowers its cavity's worst angle (capped
  at θ) by more than rounding.
- **T-P3**: gain −1 is bit-identical to 20c-1; determinism as T6.

### 20c-3

Python only: a two-polygon coverage of one class merges to one; a
millimetre vertex pair is gone at 1 cm; the partition stays valid and its
area changes by less than tolerance × perimeter; flags off are a no-op.

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
gate runs above. 20c-3 is Python only and needs no `bench.py` run.

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
| 20c-2 | `quality.hpp` (cavity, angles, acceptance) ~55; option plumbing and CLI ~25 | **~80** |
| 20c-3 | `feature_input.py` (merge, tolerance) ~35; CLI and record ~25 | **~60** |

On the worst overrun seen (+60 %), 240, 130 and 100. Each under 700.

## Questions for Ola

1. **May rasputin change land-cover outlines before meshing?** Today they
   are used exactly as given (16b R9). Two separate yeses are possible:
   (a) drop borders between two polygons of the same class, which carry no
   land-cover change (removes the two forest cuts at Lagan); (b) simplify
   outlines by up to a set distance, which removes the millimetre seams and
   most of the remaining slivers (305 instead of 709 at 2 m, 7 % fewer
   triangles). **Default: (a) yes, on by default with a class map; (b) yes
   as a flag, `--features-tolerance`, default 0 (off) until you pick a
   default distance.** Until answered, 20c-3 is not built.
2. **The soft criterion's strictness.** With P = 0 a quality point is
   skipped only if it would make its neighbourhood's worst angle worse:
   12 to 17 % fewer triangles at the same number of slivers, but the median
   smallest angle drops from about 37° to 35°. **Default: P = 0, on.**
3. **Quality during refinement too** (increment 20's C3, which you
   suspended to here)? Measured: 8 % more triangles for 5 fewer slivers on
   Lagan. **Default: no, not in 20c.**

## Not in scope

- Off-lattice Steiner points graded below a DEM cell (option c's Ruppert
  proper): not recommended (M3, "What is left").
- The noder's `MalformedOutput` at `--snap-spacing` 0.1 m and 1 m (M4): a
  defect of its own, for whoever owns the noder; not needed by 20c-3.
- The vertex rule in the quality start (R2.7): measured to do nothing.
- A size bound or sizing field (increment 14's U2).
