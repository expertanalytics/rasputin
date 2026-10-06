# Slivers in the Lagan mesh: where they are and why

Status: investigation, 2026-10-06, `@architect`. Diagnosis only: no design and
no fix. Asked by Ola ("I could see slithers in the lagan dataset after the
seams"; "Have the @architect figure out a way to investigate the slivers").
A sliver here is a triangle whose smallest angle, in plan view (x/y), is under
1°. Every number below was measured on the files named, with the scripts
described in "How it was measured"; none is an estimate unless it says so.

## Answer first

- **The visible line of slivers is not the GLO-30 tile edge at 57° N.** It is
  a straight line in the CORINE land-cover data, at northing 3 772 265.53 m
  in EPSG:3035 (CORINE's coordinate system): a seam in CORINE itself, most
  likely where two of its production areas meet (inferred from the geometry,
  not checked against the EEA's documentation). CORINE boundaries that cross
  the line have two vertices on it, at northings 3 772 265.53 and
  3 772 265.5343, a few millimetres apart: 450 segments shorter than 1 cm
  lie on it. In
  SWEREF 99 TM (the mesh's EPSG:3006) the line runs, inside the domain, from
  about 57.006° N at 13.50° E to 56.985° N at about 14.22° E, tilted about 4°
  to the parallels, which is why it looked like the 57° N tile edge: it
  crosses 57° N near 13.73° E. A second CORINE line of the same kind, at
  easting 4 537 579.53, carries 23 more slivers (F2).
- **The DEM tiles play no part.** The two GLO-30 tiles meet on one node
  lattice with no overlap and no step in height at 57° N, and the same mesh
  made without features has 13 slivers in 468 337 triangles, none near 57° N.
- **There is a second, spread-out population**, about 0.16 % of triangles
  everywhere in the basin: a DEM node or grid node placed a fraction of a
  metre to a few metres from a long CORINE edge. It is not a line, so it is
  probably not what Ola saw, but it is most of the count (about 1 365 of
  1 598). It is the case increment 20 already names as open ("a bad triangle
  whose circumcentre lies across a constraint stays ... along long
  segments"), plus two insertion paths that have no constraint feet: the
  quality start and the final check against the source DEM. Increment 20c
  (ROADMAP: a soft quality criterion that may split constraint segments) is
  where it belongs; which path makes how many is experiment E7 (the
  experiment list at the end).
- **A third population, 85 slivers, sits on two other straight CORINE
  lines** (F5): N = 3 811 923.31 and E = 4 585 680.19, with 10 to 12 % of
  the triangles within 15 m of them slivers. They have no mm segments. Each
  line is a straight cut through one coniferous-forest polygon (CORINE class
  312 on both sides of every one of its segments), with a vertex only where
  another boundary meets it, so its segments are among the longest in the
  file (median 811 m and 576 m, against 75 m for the file). The slivers are
  the spread-out population's shape (a node 2 to 3 m from a constraint edge
  hundreds of metres long), concentrated where the edges are longest. E12
  (the polygons on both sides merged, so the cuts are gone) removes 83 of
  the 85 and changes nothing else; E11 (every CORINE segment split to 100 m
  or less) removes the same 83 and two thirds of the spread-out population
  too, at 24 % more triangles. Both are input edits, workarounds, not fixes.
- **The mm-segment population belongs to input coarsening**, which increment
  16b R9 (the ninth numbered ruling in `docs/increments/16b-terrain-polygons.md`)
  left to a later increment. E2 (the experiment list at the end) shows that
  merging CORINE vertices closer than 1 cm removes 105 of the seam's 117
  slivers at no cost in triangles; snapping the features file to a 1 cm grid
  before meshing is a workaround available today, not a fix. It does nothing
  for the third population.
- **No population breaks the height tolerance**: every run below reports
  every source DEM node within the tolerance. Slivers cost angles, and
  whatever consumes the mesh downstream, not accuracy.

## The case

- Mesh: `../rasputin_scratch/sweden/lagan/lagan_clc_glo30_tol10m.vtk`, with
  `_stats.md` and `_record.json` beside it; made at master `44fa7f5` with
  `rasputin mesh --dem glo30 --cache ../rasputin_data/sweden_glo30_cache
  --out-crs EPSG:3006 --domain
  ../rasputin_data/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson
  --features ../rasputin_data/sweden_corine/lagan_clc2018_3035.geojson
  --features-crs EPSG:3035 --features-map corine --tolerance 10 --binary`.
- Re-run at master `ed12512` (the main checkout's `.venv`, its installed
  `_core` built 2026-10-05): the same 863 897 triangles, the same 1 598
  slivers. All runs below are at `ed12512`; the first table's on AC power.
- The main session's quick numbers re-checked: 863 897 triangles, 0.185 %
  under 1° (1 598); 0.73 % within 1 km of 57.0° N against 0.17 % elsewhere.
  Its "1.0 to 1.3 %" was not reproduced exactly (it depends on the strip
  chosen: 0.90 % in the 0.02° band just south of 57.0° N). The excess is
  real but, as below, belongs to the CORINE line, not to the latitude.

## What was found

### F1. Where the slivers are, and what they touch

Every vertex of the output was classed by where it can have come from:
a node of the 31 m resampled grid (x and y whole multiples of 31 m), a node of
the GLO-30 source (latitude a whole multiple of 1″ and longitude of 1.5″, to
1e-3 of a step), or an end of a constraint edge (a `LINES` cell). Every
vertex falls in one of the three: 241 005 grid, 27 000 source, 192 569
constraint.

| slivers (< 1°) | count |
|---|---|
| all | 1 598 |
| with a constraint edge as one side | 1 529 |
| corners: constraint, constraint, grid node | 955 |
| corners: three constraint vertices | 381 |
| corners: constraint, constraint, source node | 200 |
| other mixes | 62 |

So nearly every sliver lies on a constraint edge (domain outline or CORINE
border). One has DEM-made corners only.

### F2. The line is CORINE's, not the DEM's

Sliver centroids moved into EPSG:3035: of the 216 slivers under 0.05°, 125
sit within 50 m of northing 3 772 300, and 33 near easting 4 537 600.

The CORINE file (`lagan_clc2018_3035.geojson`, 7 319 features, 1 744 538
vertices) has many vertices on exactly the same northing or easting, i.e.
straight lines through the data. Counting ring vertices with each ring's
closing vertex included, the eight lines with the most are N = 3 811 923.31
(625), E = 4 519 602.54 (569), E = 4 507 398.15 (541), N = 3 772 265.53
(513), E = 4 585 680.19 (504), N = 3 772 265.5343 (497), E = 4 537 573.71
(433) and E = 4 537 579.529 (349). They are of two kinds:

- **Lines that boundaries cross.** Each boundary that crosses the line has
  a vertex on it, and no segment longer than 3 m lies along it: the seam
  has 7 (0.30 to 0.72 m), E = 4 507 398.15 has one (2.97 m), the other two
  none (distinct segments, both ends on the exact value; checked against
  the CORINE file). N = 3 772 265.53 (the seam)
  and E = 4 537 579.53 are of this kind *and* carry mm pairs: the seam's are
  the 450 segments under 1 cm in the table below, and the easting line has
  its vertices spread over values within 2 mm of 4 537 579.53 (349 at
  4 537 579.529, 53 at 4 537 579.5285, 53 at 4 537 579.53, 24 at
  4 537 579.5293). Within 15 m of the easting line or of E = 4 537 573.71,
  5.8 m west of it, there are 23 slivers in 866 triangles, and 21 of them
  have a side shorter than 10 cm: the seam's mechanism. E = 4 519 602.54 and
  4 507 398.15 are of this kind without mm pairs and carry no slivers (0 in
  the 149 and 33 triangles within 15 m of them inside the domain).
- **Lines that are borders.** N = 3 811 923.31 and E = 4 585 680.19 are
  polygon borders along their whole length. They are F5's.

Reading the rings of the seam shows what it is: in the
rings read, a CORINE boundary that crosses the line has two consecutive
vertices there, one at each value, for example

    4430606.5283 3772265.5343
    4430606.5300 3772265.5300

a segment 4.6 mm long. Along the easting line the pairs are about 6 m apart
instead. The polygons on both sides share these vertices (no gap, no overlap),
which is why a gap-and-overlap check passes. In places the line is itself a
polygon border (66 constraint edges of the mesh lie on it). Over an 80 km
stretch of the line (eastings 4 500 to 4 580 km), the kilometre on each side
has more constraint length and more triangles than the bands further out:
triangles per unit area 86 (south) and 111 (north) against 70 to 89 at 1 to
10 km, constraint length 1.37 and 1.73 against 1.18 to 1.28 (relative
figures: the domain does not fill the whole stretch). That, not the DEM, is
the higher triangle density the main session saw near 57° N.

| segments shorter than | in the CORINE file | of them on N ≈ 3 772 266 | in the mesh (inside the domain) |
|---|---|---|---|
| 1 cm | 543 | 450 | 73 |
| 10 cm | 616 | | 106 |
| 1 m | 2 981 | | 334 |

A 4 mm constraint segment is a 4 mm local feature: any triangle on it with
sides of tens of metres has an angle near 0.005°. Of the 216 slivers under
0.05°, 168 have a corner on a constraint edge shorter than 10 cm. Within
15 m of the line there are 459 triangles and 117 slivers (25 %).

This was already known in general: increment 16b, R9
(`docs/increments/16b-terrain-polygons.md`), ruled that CORINE vertices are
used as given and that "segments shorter than a cell ... cost angles, not
correctness", with coarsening left to a later increment. Lagan shows the cost
concentrated on one line, which makes it visible.

### F3. The 57° N tile edge is clean

- The GLO-30 tiles N56 E013 and N57 E013 are both 1″ by 1.5″ (3600 × 2400
  nodes, `pixel_is_area` false); N57's last row is at 57° 0′ 1″ and N56's first
  at 57° 0′ 0″. `plan_mosaic` puts them on one lattice; `assemble` over 13.40
  to 13.60° E, 56.98 to 57.02° N reports no seam, and the mean absolute
  height step between neighbouring rows is 2.05 m across 57.0° N against
  2.14 m over the window (99th percentile of rows 2.64 m). There is no row
  doubled, missing or stepped.
- Without the CORINE line, the 57° N band is ordinary. Within 1 km of 57.0° N
  and more than 1 km from the CORINE line: 24 slivers in 7 721 triangles
  (0.31 %), against 0.21 to 0.26 % further away; counts that small do not
  separate.
- Near 14° E (the other tile edge in the basin) the rate is 0.19 % within
  1 km, the same as elsewhere.

### F4. The spread-out population

Slivers more than 15 m from any of the straight CORINE lines: 1 365. Their
longest side is a constraint edge in 916; it is long (median 165 m, 10th to
90th percentile 57 to 379 m), and the third corner sits close to it
(median 0.91 m, 90th percentile 2.8 m). The third corner is a grid node in
628, a constraint vertex in 548 and a source node in 189.

For 603 of those with a grid-node corner, the node is less than half a cell
(15.5 m) from the constraint edge and its foot on the edge is more than
15.5 m from both ends. Increment 20b's rule ("constraint feet", shown in the
stats as "points moved onto lines") is meant to replace exactly such a node by
its foot on the edge. That it did not points at a path the rule does not
cover; the candidates, in the code as it stands:

1. the quality start (increment 20) inserts DEM nodes with no feet: 195 158
   of them here, against 46 560 points from DEM refinement;
2. the final check against the source DEM (`refine_points.hpp`) inserts
   points with no feet at all: 27 266 here (`dem_check_points_inserted`),
   which are 27 000 source nodes plus 266 strip points
   (`line_points_inserted`; a strip point is a point on a constraint line,
   where a grid line of the DEM crosses it or halfway between two such
   crossings, checked against the DEM: increment 15f, the edge strip). Reading
   the code, not measured: no other step inserts GLO-30 source nodes, so the
   189 source-node corners above come from it;
3. 20b's rule looks only at the constrained edges of the triangle holding the
   node (20b R1, the first numbered ruling in
   `docs/increments/20b-min-insertion-distance.md`), and a node can be close
   to an edge of a neighbouring triangle;
4. 20b's fallback (20b R5) inserts a footed node later if its error stays
   over the tolerance.

Which of these it is cannot be read from the output file, which does not
record where a vertex came from. E3 and E4 below narrow it; E7 and E8
settle it.

### F5. The third population: two straight cuts through one forest polygon

Within 15 m of N = 3 811 923.31 there are 38 slivers in 354 triangles
(10.7 %), and within 15 m of E = 4 585 680.19 47 in 387 (12.1 %): 85
slivers, inside the domain from 13.81 to 14.70° E at about 57.34° N and from
57.01 to 57.63° N at about 14.39° E. No segment under 10 cm touches either
line, and E2's 1 cm snap leaves them (34 and 61 slivers in E2's mesh).

What the lines are, from the CORINE file:

- **Each is a polygon border along its whole length**: 114 segments
  (166 km) on the northing line and 83 (131 km) on the easting line, every
  one shared by exactly two polygons.
- **The class is the same on both sides of every segment**: 312, coniferous
  forest, on all 197. These are the file's only same-class borders: of the
  37 191 km of borders shared by two polygons, 297.9 km have the same class
  on both sides, and all 297.9 km are on these two lines. A same-class
  border marks no change of land cover; it is a cut through what would
  otherwise be one polygon, most likely a production-area or tile edge or a
  split of an oversized polygon (inferred from the geometry, not checked
  against the EEA's documentation). Merging the
  polygons on either side (E12) joins exactly four features into one.
- **A vertex only where another boundary meets the line**: every one of the
  228 and 166 distinct vertices on them has an edge leaving the line (a
  T-junction). Between junctions the cut is one straight segment: median
  811 m and 576 m long, the longest 6.4 km and 12.0 km, against a median of
  75 m over the file's 1 342 157 segments. 94 of the file's 1 167 segments of
  800 m or more are on these two lines.

What the slivers are: of the 85, 56 have a constraint edge as their longest
side, 22 of these (by a 1 cm test in EPSG:3035) on the line itself and the
rest within 10 cm of it, pieces of the line's segments as split by the
meshing. That edge is long (median 544 m and 445 m) and the third corner is
close to it (median 2.2 m and 2.5 m); it is a grid node in 31, a source node
in 24 and a constraint vertex in 1. All 31 grid nodes are within half a cell
(15.5 m) of the edge with their foot more than 15.5 m from both ends: the
case 20b's feet are meant to catch, as in F4. The shortest sides are tens of
metres (median 42 m and 62 m).

So the third population is F4's mechanism, concentrated where the
constraint segments are longest. Over the whole mesh, the share of slivers
among triangles that have a constraint edge as a side rises with that
edge's length:

| constraint edge in the mesh | triangles on it | slivers among them |
|---|---|---|
| 10 to 30 m | 56 091 | 0.05 % |
| 30 to 100 m | 129 704 | 0.11 % |
| 100 to 200 m | 82 663 | 0.38 % |
| 200 to 400 m | 23 984 | 1.71 % |
| 400 to 800 m | 1 712 | 6.89 % |
| 800 m or more | 98 | 32.7 % |

(Edges under 10 m are the mm and short-segment cases of F2: 96 % of the
triangles on edges under 10 cm are slivers.) In the mesh, the constraint
edges with both ends within 10 cm of the two lines have a median length of
249 m and 280 m, against 56 m for all constraint edges; the figure depends
on that test (a 1 cm test gives 263 m on the easting line).

Two runs settle the cause without C++ (table below): **E11**, every CORINE
segment split into pieces of at most 100 m (`shapely.segmentize`) and then
snapped to 1 cm, and **E12**, the four same-class polygons merged so that
the cuts are gone. Both take the two lines from 85 slivers to 2. E12 changes
nothing else (the seam keeps 117, the rest of the basin 1 373), so the cuts
themselves, not anything about where they lie, make the slivers. E11 also
takes the rest of the basin from 1 365 to 438, which says that long segments
are most of F4 too; it costs 24 % more triangles (863 729 to 1 073 613) and
leaves the worst angle where it was (0.00048°).

Review round 1 counted 116 slivers on the other straight lines, adding the
two easting lines of F2 as 18 and 13; 8 slivers lie within 15 m of both, so
those two hold 23, and they are mm segments (F2), not this population.

What stays open is the same question as F4's: which insertion path puts the
node within metres of the long edge (candidates 1 to 4 of F4). E7 and E8
answer it for both populations; their rows below now name the two lines.

The 20c design answers E7 and E8 (`docs/increments/20c-soft-quality.md` on
branch `worktree-soft-quality`, commit `7fcac1e0`, **not reviewed yet**). It
tagged where each vertex came from and attributed the 1 598 slivers by the
vertex opposite the longest side: quality start (increment 20) 639, input
vertex 516, final check's source-DEM nodes 219, refinement 38, and 185 with
a side under 10 cm.

## Runs done here (no C++ build)

All from the main checkout's `.venv` at `ed12512`, outputs in the session's
scratchpad (not kept). Same command as the case unless stated.

"Seam" is the strip within 15 m of CORINE's line N = 3 772 265.53 (EPSG:3035);
"57° N" the strip within 1 km of 57.0° N.

| run | triangles | < 1° | < 0.05° | < 10° | worst | slivers on a < 10 cm constraint edge | seam: triangles, slivers | 57° N: rate | rest: rate |
|---|---|---|---|---|---|---|---|---|---|
| E0 the case | 863 897 | 1 598 (0.185 %) | 216 | 5.36 % | 0.00041° | 189 | 459, 117 | 0.73 % | 0.17 % |
| E1 no features | 468 337 | 13 (0.003 %) | 0 | 1.43 % | 0.074° | 0 | 63, 0 | 0.00 % | 0.00 % |
| E1 no features, tolerance 2 m | 5 196 462 | 64 (0.001 %) | 3 | 0.96 % | 0.038° | 0 | 875, 0 | 0.00 % | 0.00 % |
| E2 CORINE on a 1 cm grid | 863 729 | 1 480 (0.171 %) | 88 | 5.35 % | 0.00041° | 63 | 331, 12 | 0.35 % | 0.17 % |
| E3 `--start-min-angle 0` | 582 890 | 8 665 (1.487 %) | 185 | 10.58 % | 0.00030° | 198 | 205, 69 | 5.60 % | 1.41 % |
| E4 `--no-constraint-feet` | 866 085 | 2 047 (0.236 %) | 269 | 5.53 % | 0.00039° | 189 | 459, 117 | 0.73 % | 0.22 % |
| E5 tolerance 5 m | 1 716 879 | 2 910 (0.169 %) | 288 | 3.77 % | 0.00041° | 191 | 570, 120 | 0.53 % | 0.16 % |

Read together:

- **E1:** without features the mesh has almost no slivers, at 10 m and at
  2 m, and none near 57° N; the 13 at 10 m are scattered along the domain
  outline. Without features there is also no rise in triangle density at
  57° N (triangles per 0.01° band between 13.5 and 14.2° E stay between
  2 238 and 3 010 at 10 m). The DEM and its tile edge are cleared.
- **E2:** merging CORINE's mm vertex pairs takes the seam from 117 slivers to
  12 and the 57° N strip from 0.73 % to 0.35 %; the rest of the basin does
  not move (0.17 %). The visible line is the mm segments. The worst angle does
  not move: that triangle belongs to the spread-out population (F4), a grid
  node 0.8 mm from a 245 m CORINE edge. In EPSG:3035 its centroid is at about
  (4 558 723, 3 738 444); the grid node is at (4 558 732, 3 738 444), and the
  smallest angle (0.00041°) is at the edge's west end, (4 558 596,
  3 738 441).
- **E3:** without the quality start, slivers rise eightfold, to 1.49 %. The
  quality start removes far more slivers than it could be leaving behind, so
  "the quality start made them" is not the main story; but its own rule, in
  increment 20 R5 (the fifth numbered ruling in
  `docs/increments/20-start-quality.md`), says that under C1 (a) (that
  increment's choice of never putting a point on a constraint segment) "a
  bad triangle whose circumcentre lies across a constraint stays. It happens
  at reflex corners and along long segments." That is the shape of F4's flat
  triangles on long CORINE edges.
- **E4:** increment 20b's feet remove 449 slivers (2 047 to 1 598) and leave
  the seam untouched, as expected. 20b replaces a node that is within a
  distance ε of a constraint segment by its foot on the segment (ε is set by
  the slope, at most half a cell: 20b R3); but it never takes a foot within
  ε of a segment's end, and a 4 mm segment is all ends.
- **E5:** at 5 m the mesh doubles and the slivers nearly double (2 910); the
  rate stays at about 0.17 %, and the seam keeps its 120. The spread-out
  population grows with the mesh; the seam's is fixed by the data.

The four line strips and E11, E12 (E2 re-run with them; the later runs on
battery, which changes the time, not the output). Strips are 15 m wide on
each side, in EPSG:3035; "rest" is everything outside the four.

| run | triangles | < 1° | < 0.05° | worst | seam N 3 772 265.53 | mm line E 4 537 579.53 (with E 4 537 573.71) | cut N 3 811 923.31 | cut E 4 585 680.19 | rest |
|---|---|---|---|---|---|---|---|---|---|
| E0 the case | 863 897 | 1 598 | 216 | 0.00041° | 459, 117 | 866, 23 | 354, 38 | 387, 47 | 1 373 (0.159 %) |
| E2 CORINE on a 1 cm grid | 863 729 | 1 480 | 88 | 0.00041° | 331, 12 | 846, 8 | 345, 34 | 433, 61 | 1 365 (0.158 %) |
| E11 segments split to ≤ 100 m, then 1 cm grid | 1 073 613 | 461 | 95 | 0.00048° | 348, 12 | 849, 9 | 233, 0 | 227, 2 | 438 (0.041 %) |
| E12 the cuts merged away | 863 483 | 1 515 | 210 | 0.00041° | 459, 117 | 866, 23 | 287, 0 | 297, 2 | 1 373 (0.159 %) |

Each strip cell is "triangles, slivers". E11's input has 2 638 970
vertices against 1 744 538; all its features stay valid. Every run reports
every source DEM node within the tolerance.

## How it was measured

- The `.vtk` read directly (binary legacy VTK: `POINTS`, `LINES`,
  `POLYGONS`, big-endian), angles in plan view from the x/y coordinates.
- Coordinates moved with `pyproj` (`always_xy=True`) to EPSG:4326 and
  EPSG:3035 for the latitude bands and the CORINE lines.
- The GLO-30 window assembled with the repository's own `CacheRepository`,
  `plan_mosaic` and `assemble`.
- The CORINE line search: exact repeated northings and eastings among all
  ring vertices, 40 or more vertices on one value.
- E2's input: every feature through `shapely.set_precision(geometry, 0.01)`
  (EPSG:3035 metres). All 7 319 stay valid; segments under 1 cm go from 543
  to 70, under 10 cm from 616 to 144. A grid snap treats both sides of a
  shared border alike, so the partition stays a partition.
- E11's input: `shapely.segmentize(geometry, 100.0)`, then the same
  `set_precision(…, 0.01)` so that the points added on a shared border from
  its two sides coincide. E12's input: the features joined through borders
  with the same `Code_18` on both sides (a union-find over shared segments),
  each group merged by `shapely.union_all`: one group of four, all valid.
- F5's sliver rate by constraint edge length: each triangle side that is a
  `LINES` cell of the mesh, binned by that cell's length.

These are throw-away scripts. The experiments that need them again are
listed below with the persona that would own a kept version.

## Experiment list, cheapest first

Each with what result points where. "Done" rows are in the results table
above.

| # | experiment | cost | who | a result that points where |
|---|---|---|---|---|
| E0 | Classify sliver corners and position (F1, F2) | seconds, Python on the `.vtk` | done | — |
| E1 | The same mesh without `--features` | 90 s run | done | slivers gone: features cause them; still there: DEM or outline |
| E2 | CORINE snapped to a 1 cm grid (`shapely.set_precision`), which merges the mm vertex pairs | 90 s run | done | line gone: the mm segments are the cause of the line |
| E3 | `--start-min-angle 0` (no quality start) | 90 s run | done | spread-out slivers fall: the quality start's nodes are candidate 1 |
| E4 | `--no-constraint-feet` | 90 s run | done | how much 20b removes today; the rate it leaves is the floor without feet |
| E5 | Tolerance 5 m and 2 m, with and without features | 2 to 10 min each | done: 5 m with features, 2 m without; 2 m with features (about 5 M triangles) left to `@perf` | how the spread-out rate scales with the number of DEM insertions near edges |
| E6 | Other CORINE extracts: search for mm vertex pairs on straight lines (as F2) | seconds, Python | done for the other two Swedish files: `ljungan_flasjo_clc2018_3035.geojson` has 183 segments under 1 cm, 125 of them on one line near northing 4 366 070; `acklingen_clc2018_3035.geojson` has 68, at most 19 on one line. Whether those lines fall inside the meshed domains was not checked | a seam in the file inside the domain predicts a line of slivers before any mesh is made |
| E7 | Vertex origin in the output: tag each vertex as start, quality start, refinement node, foot, final-check point or strip point, in a scratch copy of the C++ core | one C++ build plus a 90 s run | `@perf` (scratch copy, not committed); answered by the 20c design, unreviewed (end of F5) | the share of spread-out slivers per origin, and of the slivers on the two cuts of F5, settles candidates 1, 2 and 4 of F4 |
| E8 | For each spread-out sliver whose apex is a refinement node: was a constrained edge within ε in a neighbouring triangle when it went in? Logged in the same scratch copy as E7 | with E7 | `@perf`; answered by the 20c design, unreviewed (end of F5) | yes in most: 20b R1's "own triangle only" is the gap (candidate 3) |
| E9 | The 5 m and 1 m runs over a GLO-30 window that has no features but crosses 57° N and 14° E, rate per km band | 10 min | `@perf` | no excess: the tile edge stays cleared at fine tolerances, where more source nodes go in |
| E11 | Every CORINE segment split to 100 m or less, then the 1 cm grid of E2 (F5) | 90 s run | done | the cuts' slivers gone, and two thirds of the rest: long segments are the cause of F5 and most of F4 |
| E12 | The polygons on each side of a same-class border merged, which removes the two cuts (F5) | 90 s run | done | the cuts' slivers gone and nothing else moves: the cut segments, not their position, cause F5 |
| E13 | Other CORINE extracts: same-class shared borders (as F5) | seconds, Python | done for the other two Swedish files: `ljungan_flasjo_clc2018_3035.geojson` has 11.6 km in 14 segments (median 593 m), `acklingen_clc2018_3035.geojson` none. Whether they fall inside the meshed domain was not checked | a same-class border inside the domain predicts F5's slivers before any mesh is made |
| E10 | Failing tests for what E7 and E8 find: a node close to a long constraint edge, inserted by the quality start, by the final check, and from a neighbouring triangle; a polygon pair with a 4 mm step; a straight border between two polygons of one class with segments hundreds of metres long | a red suite | `@tester`, after a design | — |

What is **not** worth running: more DEM-seam checks for GLO-30 inside one
latitude band. The tiles of a band share one lattice and the mosaic is one
grid before anything is meshed, so there is no tile edge left for the mesh to
see.

## Bearing on huge catchments (São Francisco)

- **DEM tile edges.** Within one GLO-30 latitude band (spacing 1″ × 1″ below
  50°, the band São Francisco lies in) the tiles share one lattice and are
  assembled into one grid before meshing, so their edges leave no mark on the
  mesh (F3). This does not cover heights that disagree from tile to tile in
  the data itself: `dem_seams` checks only where tiles overlap, and GLO-30
  tiles do not overlap. ANADEM is one file, but its vegetation correction was
  computed in 93 tiles of 5° × 5° (`docs/research/anadem-accuracy.md`),
  so a step in height at a 5° line is possible; it would show as a line of
  DEM insertions, not of slivers. A check for it is a separate question from
  this one.
- **GLO-30 across a latitude band edge** (50°, 60°, 70°) changes the
  longitude spacing (1″ → 1.5″ → 2″ → 3″), so the tiles are on two lattices.
  Not São Francisco's case; for northern Sweden it is a question for the
  mosaic (`MixedGridError`), not for slivers.
- **CORINE's production seams** are the transferable finding, but not for
  Brazil, which CORINE does not cover. The general lesson is that land-cover
  data can carry mm-scale segments along its own production seams. São
  Francisco's land cover is planned from MapBiomas, a raster that rasputin
  would turn into polygons itself (`docs/research/sao-francisco-basin.md`,
  `docs/research/raster-to-vector.md`). It then has no seams of CORINE's kind
  unless that conversion works tile by tile; if it does, the E6 search
  belongs in its acceptance.
- **Long straight segments** are the hazard that transfers to any data
  (F5): the sliver share among triangles on a constraint edge rises from
  0.1 % at 30 to 100 m to 7 % at 400 to 800 m. A raster-to-vector conversion
  that works tile by tile and does not merge same-class polygons across
  tile edges would make F5's cuts on purpose, so the E13 search belongs in
  its acceptance too; a simplifier that merges collinear runs into long
  segments would make them as well.
- **The spread-out population** grows with the number of constraint edges and
  insertions, so it is the one that scales with the basin. Its rate here,
  about 0.16 %, would be on the order of ten thousand slivers on a mesh of
  millions of triangles (an estimate, not a measurement).

## Review

- Review round 1 of c83658f1 (docs/research/lagan-slivers.md, docs only, 0 production lines): CHANGES REQUESTED. The claim in F2 that the seam lines are CORINE's "most populated" straight lines is false, 116 slivers sit on other straight CORINE lines that have no mm segments and the doc leaves them unexplained, and labels R1/R5/R9, ε, "strip point" and the first uses of E2/E7 need expanding (@reviewer)
- Review round 2 of c83658f1..cd93b2cb (docs/research/lagan-slivers.md, docs only, 0 production lines): CHANGES REQUESTED. F2's "Lines that boundaries cross. No segment lies along the line" is false: the seam N = 3 772 265.53 has 7 segments lying along it (0.30 to 0.72 m) and E = 4 507 398.15 has one (2.97 m). Every other re-measured figure in F2, F5, the edge-length table, the four-strip table and E13 reproduces (@reviewer)
