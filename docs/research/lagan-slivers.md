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
  SWEREF 99 TM (the mesh's EPSG:3006) the line runs from about 57.006° N at
  13.53° E to 56.986° N at 14.19° E, tilted about 4° to the parallels, which
  is why it looked like the 57° N tile edge: it crosses 57° N near 13.73° E.
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
  where it belongs; which path makes how many is experiment E7.
- **The line's population belongs to input coarsening**, which increment 16b
  (R9) left to a later increment. E2 shows that merging CORINE vertices
  closer than 1 cm removes 105 of the line's 117 slivers at no cost in
  triangles; snapping the features file to a 1 cm grid before meshing is a
  workaround available today, not a fix.
- **Neither population breaks the height tolerance**: every run below reports
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
  slivers. All runs below are at `ed12512`, on AC power.
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
straight lines through the data. The most populated are N = 3 772 265.53
(513 vertices) and N = 3 772 265.5343 (497), E = 4 537 573.71 (433) and
E = 4 537 579.529 (349). Reading the rings shows what they are: in the
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
   source nodes with no feet at all: 27 266 here, and the 189 source-node
   corners above can only come from it;
3. 20b's rule looks only at the constrained edges of the triangle holding the
   node (20b R1), and a node can be close to an edge of a neighbouring
   triangle;
4. 20b's fallback (R5) inserts a footed node later if its error stays over
   the tolerance.

Which of these it is cannot be read from the output file, which does not
record where a vertex came from. E3 and E4 below narrow it; E7 and E8
settle it.

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
  node 0.8 mm from a 245 m CORINE edge, near (4 558 596, 3 738 441) in
  EPSG:3035.
- **E3:** without the quality start, slivers rise eightfold, to 1.49 %. The
  quality start removes far more slivers than it could be leaving behind, so
  "the quality start made them" is not the main story; but its own rule, in
  increment 20 R5, says that under C1 (a) (no point ever put on a constraint
  segment) "a bad triangle whose circumcentre lies across a constraint stays.
  It happens at reflex corners and along long segments." That is the shape
  of F4's flat triangles on long CORINE edges.
- **E4:** increment 20b's feet remove 449 slivers (2 047 to 1 598) and leave
  the seam untouched, as expected: a foot is never taken within ε of a
  segment's end, and a 4 mm segment is all ends.
- **E5:** at 5 m the mesh doubles and the slivers nearly double (2 910); the
  rate stays at about 0.17 %, and the seam keeps its 120. The spread-out
  population grows with the mesh; the seam's is fixed by the data.

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
| E7 | Vertex origin in the output: tag each vertex as start, quality start, refinement node, foot, final-check point or strip point, in a scratch copy of the C++ core | one C++ build plus a 90 s run | `@perf` (scratch copy, not committed) | the share of spread-out slivers per origin settles candidates 1, 2 and 4 of F4 |
| E8 | For each spread-out sliver whose apex is a refinement node: was a constrained edge within ε in a neighbouring triangle when it went in? Logged in the same scratch copy as E7 | with E7 | `@perf` | yes in most: 20b R1's "own triangle only" is the gap (candidate 3) |
| E9 | The 5 m and 1 m runs over a GLO-30 window that has no features but crosses 57° N and 14° E, rate per km band | 10 min | `@perf` | no excess: the tile edge stays cleared at fine tolerances, where more source nodes go in |
| E10 | Failing tests for what E7 and E8 find: a node close to a long constraint edge, inserted by the quality start, by the final check, and from a neighbouring triangle; a polygon pair with a 4 mm step | a red suite | `@tester`, after a design | — |

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
- **The spread-out population** grows with the number of constraint edges and
  insertions, so it is the one that scales with the basin. Its rate here,
  about 0.16 %, would be on the order of ten thousand slivers on a mesh of
  millions of triangles (an estimate, not a measurement).
