# Increment 32: land-cover borders simplified within a band, each class's area kept

**Status:** built. Designed on `0b5e0d4c`, Ola's four rulings (section
12), design review round 2 approved. Red step `44f25968` (`@tester`), green
step `82c0a7d4` (`@developer`, 519 net production lines), mutation round
`30f939b6` (`@tester`, one new test case); `@architect`'s rulings on both
hand-backs are in section 14. Code review round 1 done: changes requested
(claims in comments, docstrings and `project_structure.md` that the change
made false), fixed by `@developer`'s `4d4f0411` and `@architect`'s docs
commit after it; red-step notes put in the past tense (`29fbe35c`). Code
review round 2 asked for one more red note (Review section). Next: that
note, a short round 3, then `@perf`'s short timing check
(`tools/bench_quick.py`, base and branch back to back, section 10). Not refine or mesh code, so no
`bench.py` acceptance run.

## 1. What Ola asked

Ola, 2026-10-08: "So, one thing that we need to implement is a CORINE
simplification, similar to what we do on the outline of the auto-catchments.
We get too many triangles (150k) and still only 777m vtol. We need to research
an algorithm for this."

The probe of step 1 (`0b5e0d4c`,
`docs/benchmarks/2026-10-08/clc-simplify/README.md`) measured today's switch,
`--features-tolerance` (GEOS's coverage simplifier, Visvalingam-Whyatt by
area). It asked two questions; Ola answered "Defaults on all." (2026-10-08),
so:

1. The simplifier keeps every land-cover border within a **distance band** of
   its source border, tied to CORINE's accuracy, starting at **50 m**, and
   **each class keeps its area exactly**.
2. The research also covers the **default start angle** (`--start-min-angle`)
   when land cover is on, because the 25° start pass doubles the 50 m mesh.

What the probe found, in short: the 777 m is the height error of a mesh with
no height refinement and is not the land cover's doing (a height target needs
`--tolerance`); today's tolerance is not a distance bound (borders move up to
3.25 times the value given) and keeps no class area; the clean-up costs 1.1 to
1.4 s, of which the simplification 0.12 s.

## 2. Prior art: legacy and literature

### Literature, read

"Read" below means what was read: the full text where it could be had, the
abstract and the authors' code where the full text is paywalled. Publisher
pages (Taylor & Francis, ACM, SIAM) refused the fetch (HTTP 403); abstracts
came from OpenAlex and Crossref, by DOI.

- **Kronenfeld, Stanislawski, Buttenfield and Brockmeyer 2020**,
  "Simplification of polylines by segment collapse: minimizing areal
  displacement while preserving area", *Int. J. Cartography* 6(1):22-46,
  doi:10.1080/23729333.2019.1631535 (APSC). Read: the abstract, and the first
  author's MIT-licensed Python implementation, `apsc.py` in
  `github.com/geobarry/line-simplify` (June 2020). The abstract: segment
  collapse to Steiner points "under the constraint that the areas of
  adjoining polygons are preserved exactly", self-intersections avoided "by
  testing for intersections with two new line segments". The code: one
  polyline at a time, its two end points never collapsed; the replacement
  point E lies on the line parallel to A-D where A-E-D encloses what A-B-C-D
  did (`__equalAreaLine`), placed where that line meets A-B or C-D, chosen
  by a geometric rule (`__placement_AP_EAmin`); priority is the areal displacement; the stopping
  rule a displacement or a point count; the topology check tests the two new
  segments against the line's own segments and, optionally, against *other*
  lines held fixed (`other_pt_lists`). **What it gives:** the operation this
  increment uses, and the fact that it keeps the area on both sides of a
  shared border, so one collapse keeps both neighbours' areas. **What it does
  not:** no distance bound (its stopping rule is displacement or count), no
  coverage (other lines are fixed obstacles, not simplified together), and
  no test that a whole feature (an island) stays on its side when no segment
  crosses.
- **Buchin, Meulemans, van Renssen and Speckmann 2016**, "Area-preserving
  simplification and schematization of polygonal subdivisions", *ACM TSAS*
  2(1), doi:10.1145/2818373. Read: the abstract (the TU/e author version
  refused the fetch). The edge-move: area and topology kept, no new
  orientations, quadratic time, proven to reduce any non-convex simple
  polygon. **Gives:** the subdivision case with area and topology. **Does not
  (from the abstract):** a distance bound to the source; its aim is
  schematization. An open implementation exists (`jakoblistabarth/schemapify`).
- **de Berg, van Kreveld and Schirra 1998**, "Topologically correct
  subdivision simplification using the bandwidth criterion", *CaGIS*
  25(4):243-257, doi:10.1559/152304098782383007. Read: the abstract, and the
  full text of the conference version, "A new approach to subdivision
  simplification", Auto-Carto 12 (1995), from p. 79 (open at cartogis.org).
  The subdivision is cut into **chains between junctions** (vertices of
  degree three or more) and leaves; "Keep the positions of all leafs and
  junctions fixed"; each chain is simplified so that no point of it is
  further than ε from its simplification, it does not cross itself or other
  chains, and a set of points stays on the same side, where "we temporarily
  add to the set P of points all vertices of other chains of the
  subdivision". O(n(n + m) log n) per chain. **Gives:** exactly this
  increment's frame: chains between fixed junctions, other chains' vertices
  as side points. **Does not:** keep area (vertices are a subset of the
  input, no Steiner points); its ε is one-sided (source to simplification).
- **Saalfeld 1999**, "Topologically consistent line simplification with the
  Douglas-Peucker algorithm", *CaGIS* 26:7-18,
  doi:10.1559/152304099782424901. Read: the abstract. A test added to the
  stopping condition of Douglas-Peucker-like (vertex-subset) algorithms keeps
  the line consistent with itself and its neighbours; a dynamic convex hull
  finds conflicts. **Gives:** the side-of-feature idea for vertex subsets.
  **Does not:** keep area; no Steiner points.
- **Estkowski and Mitchell 2001**, "Simplifying a polygonal subdivision
  while keeping it simple", *Proc. 17th SoCG*, pp. 40-49,
  doi:10.1145/378583.378612. **It exists** (the probe cited it from memory;
  checked by DOI in OpenAlex and Semantic Scholar). Read: the abstract.
  Simplifying a subdivision within an error bound, keeping its topology, with
  no Steiner points, is MIN PB-complete: unless P = NP no polynomial
  algorithm comes within a factor n^(1/5) of the fewest vertices; heuristics
  work well in practice. **Gives:** why this increment is greedy and claims
  no minimum.
- **Guibas, Hershberger, Mitchell and Snoeyink 1993**, "Approximating
  polygons and subdivisions with minimum-link paths", *Int. J. Comput. Geom.
  Appl.* 3(4):383-415 (UBC TR-92-05). Read: the abstract. Fatten the object
  (convolve with disks: a band) and approximate inside it; for subdivisions or
  chains with no self-intersections the best approximation is NP-hard.
- **Mendel 2018**, "Area-preserving subdivision simplification with topology
  constraints: exactly and in practice", *ALENEX 2018*, pp. 117-128,
  doi:10.1137/1.9781611975055.11. Read: the abstract. Removes degree-two
  vertices so that given points stay in their faces, each face's area
  changes by at most a factor, and the lines stay within ε; heuristic for
  continental instances in seconds, ILP for city-size optimum. **Gives:** the
  nearest published combination (area, distance, topology on a subdivision).
  **Does not:** keep area exactly (bounded factor), and adds no Steiner
  points.
- **Haunert and Wolff 2010**, "Area aggregation in map generalisation by
  mixed-integer programming", *IJGIS* 24(12):1871-1897,
  doi:10.1080/13658810903401008. Read: the abstract. Merging areas too small
  for the target scale into neighbours, minimising class change, is NP-hard;
  MIP with heuristics. **Gives:** the method if small pieces are ever to be
  merged (section 5.2 rules against it for now).
- **JTS/GEOS `CoverageSimplifier` (Davis 2023)**, what `--features-tolerance`
  calls today through shapely 2.2.0 on GEOS 3.14.1. Read: the JTS source
  (`CoverageSimplifier.java`, `TPVWSimplifier.java`, `CoverageRingEdges.java`,
  master). Coverage edges are cut at nodes ("inner vertices shared by three
  or more polygons, or boundary vertices shared by two or more") and simplified once
  each; Visvalingam-Whyatt by corner area against the *current* line, so
  removals compound; a corner is removable if no vertex of a nearby edge lies
  in its triangle; "Rings smaller than the area tolerance are removed where
  possible"; the tolerance "equates roughly to the maximum distance" and is
  "the square root of the area tolerance". **Gives:** the edge extraction and
  the vertex-in-corner test this design mirrors. **Does not:** bound the
  distance (the probe measured 3.25 times) or keep any area.
- **Visvalingam and Whyatt 1993**, *Cartographic J.* 30:46-51,
  doi:10.1179/000870493786962263, and **Bose, Cabello, Cheong, Gudmundsson,
  van Kreveld and Speckmann 2006**, *J. Discrete Algorithms* 4(4):554-566,
  doi:10.1016/j.jda.2005.06.008: the baseline and the fewest-vertex
  area-preserving path, as in `docs/increments/22-auto-catchment.md`.
- **CORINE Land Cover 2018** product page (Copernicus Land Monitoring
  Service, `land.copernicus.eu/en/products/corine-land-cover/clc2018`):
  positional accuracy "100 m or better", thematic accuracy "≥ 85%", "25
  ha/100 m" (minimum mapping unit 25 ha, minimum width 100 m).

**What this increment takes and where it departs.** The frame is de Berg,
van Kreveld and Schirra's (chains between fixed junctions, other chains'
vertices as side points); the operation is Kronenfeld et al.'s collapse,
which keeps both neighbours' areas exactly; the band is a two-sided Hausdorff
bound checked by a new, cheap test (section 7). Departures, each with its
reason: (a) the stopping rule is the band, not a count or displacement,
because Ola asked for a band; (b) the priority is the band deviation, as in
`reduce_ring`, not areal displacement, so the two reductions behave alike;
(c) the placement tries both of APSC's lines and keeps the one with the
smaller deviation (APSC picks by a geometric rule); no guarantee is lost,
since any point on the equal-area line keeps the area; (d) all chains are
simplified in one global order with a shared crossing index, not one chain
at a time against fixed others.

### Novelty

Searched (2026-10-08): the papers above and their citing work; web searches
for area-preserving simplification of subdivisions, coverages and shared
boundaries (2019-2026), and for follow-ups to APSC. Found: exact area with
topology (APSC, Buchin et al.), a band with topology and no Steiner points
(de Berg et al., Estkowski and Mitchell), a bounded area factor with a band
and no Steiner points (Mendel). Not found: a coverage simplifier that keeps
every face's area exactly *and* holds a two-sided Hausdorff band, with
topology, using Steiner points. Nor the check of section 7 (anchors on the
source, which turn vertex-to-segment tests into a proved Hausdorff bound at
the cost of today's test). **No novelty is claimed here.** Two limits on the
search: the full texts of Kronenfeld et al. and Buchin et al. were not read,
and Mendel 2018 only by its abstract. A write-up that wants to claim the
combination must read those three first.

### Legacy

Nothing to carry. The archive's simplification is CGAL's 3-D mesh edge
collapse, not vector borders:

```
$ git grep -il -E 'simplif|douglas|visvalingam|coverage|land.?cover|corine' legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/geo_tiff_reader.py
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/rasputin/land_cover_repository.py
legacy-archive:legacy/rasputin/mesh.py
legacy-archive:legacy/rasputin/tin_repository.py
legacy-archive:legacy/rasputin/triangulate_dem.h
legacy-archive:legacy/rasputin/web_visualize.py
legacy-archive:legacy/rasputin/wfs_repository.py
legacy-archive:legacy/tests/test_gml_repository.py
legacy-archive:legacy/tests/test_land_cover_repository.py
$ git grep -n -i -E 'simplif|douglas|visvalingam' legacy-archive -- legacy
```

The second grep's hits are all `CGAL::Surface_mesh_simplification`
(Lindstrom-Turk edge collapse of the 3-D TIN, `legacy/bindings.cpp`,
`legacy/rasputin/triangulate_dem.h`) and the Python `simplify(ratio=...)`
calls on it; the land-cover files match on "land cover", not on
simplification.

## 3. How `reduce_ring` works today, and what carries over

`terrain::vector_simplify::reduce_ring<K>`
(`include/terrain/vector_simplify/area_collapse.hpp`, bound as
`_core.reduce_ring`, used by `catchment.py` for the auto-catchment outline),
read from the code:

1. Input: one closed ring, counter-clockwise, open (first vertex not
   repeated), a tolerance, and keep-points (the catchment seed).
2. A collinear pass drops every vertex exactly between its neighbours
   (`include/terrain/vector_simplify/area_collapse.hpp@b39426c0:177`).
3. The ring is a doubly linked list. Each current edge stands for a range of
   fine (source) vertices: `start` and `count`.
4. For each vertex B, with A before it and C, D after: the candidate point E
   lies on line A-B or on line C-D, where the triangle fan A-E-D has the
   signed area of A-B-C-D (or on the foot of B-C's midpoint when both lines
   are parallel to A-D), and the one with the smaller deviation is kept
   (`include/terrain/vector_simplify/area_collapse.hpp@b39426c0:226-245`). Replacing B, C by E keeps the
   ring's area up to rounding.
5. The deviation (`deviation_of`, `include/terrain/vector_simplify/area_collapse.hpp@b39426c0:131-151`): the
   best split of the fine range between A-E and E-D, the largest distance of
   a fine vertex to its edge; and E's distance to the *nearest fine segment
   anywhere in the range* (`include/terrain/vector_simplify/area_collapse.hpp@b39426c0:146-149`).
6. A global heap ordered by (deviation, id); a candidate over the tolerance
   is never queued; the loop stops at four vertices.
7. At pop: A-E and E-D are tested with the exact `noding::classify` against
   every current edge in the grid cells of the box of A, E, D (touching only
   where they share an end), and every keep-point must have winding number
   zero about A-B-C-D-E and lie on neither new edge
   (`include/terrain/vector_simplify/area_collapse.hpp@b39426c0:275-312`).
8. Apply, then re-evaluate the four candidates whose vertices changed.
   Serial; the output depends only on the input.

**Can it run per shared border, keeping both neighbours' areas?** Yes, with
four changes; the prototype of section 4 made them and measured the result.

- *Both areas.* The collapse keeps the signed area between the old chain
  A-B-C-D and the new A-E-D at zero. That region's area is what one side of
  the border gains and the other loses, so both neighbours keep their area,
  and so does every face, since every face is bounded by borders. Nothing
  new is needed: the same area rule.
- *Fixed ends.* A border runs from junction to junction. A junction (a
  vertex where three or more borders meet, or where a border meets the
  outline) never moves: it may be A or D, never B or C. A border with no
  junction (an island's whole outline) is a closed chain, as today's ring.
- *Topology across the coverage.* One crossing index over all borders. Two
  tests per collapse: (i) A-E and E-D cross or touch no edge of any border,
  except at a shared end (A with its own neighbour and with other borders'
  first edges at junction A; D likewise); (ii) **no vertex of any other
  border, and no vertex of this border outside A..D, lies in the region
  swept** (winding number about A-B-C-D-E not zero) **or on A-E or E-D**.
  Test (ii) is today's keep-point test with every nearby vertex as a
  keep-point; it is de Berg et al.'s "all vertices of other chains" and
  JTS's vertex-in-corner test. Without it, a border could sweep over a small
  island whole, crossing nothing: the island would change face.
- *The band, enforced.* Today's deviation is not a distance bound: E is
  checked against the nearest fine segment anywhere in the range, so the
  fine segment that straddles the split may pass further from the new chain
  than the tolerance. `docs/increments/22-auto-catchment.md` records this
  ("Measured, not guaranteed: the Hausdorff distance", with an open
  question for Ola). The prototype hit it on CORINE: one border of the
  German case at 58.6 m with a 50 m band. Section 7's anchored check closes
  it with a proof, at today's cost.

## 4. The evidence: a prototype on the German fused case

A throwaway Python prototype (not production code):
`docs/benchmarks/2026-10-08/clc-simplify/design-probe/border_apsc.py`. It
re-enacts the land-cover stage (CORINE moved to EPSG:25832, clipped to the
fused outline, `coverage_clean` at 1 m, one polygon per class), cuts the
coverage into borders between junctions (by exact coordinates; an edge used
by one ring is on the outline and stays fixed), runs the collapse of section
3 with the four changes in floating point, rebuilds the polygons and writes
a GeoJSON; then `runs.sh` meshes each with today's `rasputin mesh`
(`--features-tolerance 0`, so the stage only repairs, merges and applies the
outline rule). Mac on AC power. Output: `prototype.txt` and `meshes.txt` in
that folder; meshes in
`rasputin_scratch/germany/isar-loisach/clc_simplify/proto/`.

The coverage: 18 classes in 819 pieces, 151 579 ring vertices, 1 576
junctions, 2 393 borders (481 on the outline, fixed; 27 closed).

| band | polygon vertices | collapses | refused by topology | largest class area change | Hausdorff source to simplified, max | simplified to source, max | ground now in another class |
|---|---|---|---|---|---|---|---|
| 0 | 152 426 | 0 | 0 | 2e-13 % | 0 | 0 | 0 |
| 20 m | 53 082 | 49 672 | 17 | 2.4e-11 % | 20.0 m | 20.0 m | 13.7 km² (1.15 %) |
| **50 m** | **31 878** | 60 274 | 46 | 2.6e-11 % | **50.0 m** | **50.0 m** | 40.2 km² (3.36 %) |
| 100 m | 20 430 | 65 998 | 145 | 5.8e-11 % | 100.0 m | 99.7 m | 77.9 km² (6.52 %) |
| 50 m, today's check (section 3) | 31 866 | 60 280 | 52 | 2.6e-11 % | **58.6 m** (1 border) | 50.0 m | 40.2 km² |

Every output is a valid coverage (`shapely.coverage_is_valid`) of valid
polygons. Hausdorff: each border's source and simplified polyline, both
`segmentize`d to 0.5 m, every point's distance to the other line.

**Ground in another class** is the area where the simplified class polygon
lies outside its source (the probe's "labelled with a different class",
here measured on the polygons before meshing, not on a mesh). It is the
price of any band: border length times mean displacement. At 50 m it is
3.4 % of the catchment, every class's total unchanged, and it lies within
CORINE's own positional accuracy (100 m). For comparison, the probe
measured today's switch on meshes: 0.7 % at FT 30 (borders moved up to
147 m) and 4.9 % at FT 100 (up to 325 m), class totals off by up to 1.9 %.

## 5. Rulings on the four research questions

### 5.1 Per shared border

Yes: section 3's four changes. The pieces are published (section 2); the
anchored band check was written for this design.

### 5.2 Small pieces: not merged

Of the 819 pieces, 169 are under CORINE's 25 ha minimum mapping unit (7.09
km² in all), and **160 of them touch the outline**: they are the cut-off
edges of larger CORINE polygons, real ground of that class, not mapping
noise. (CORINE maps nothing under 25 ha; where the other 9 come from was
not checked.) They hold 3 017 of 151 579 vertices (2 %). Ruling: **no
merge**. Merging would relabel real ground and break the exact class area
Ola ruled, for 2 % of the vertices. The collapse cannot shrink a piece
away (its area is kept), and a closed border keeps at least four vertices.
If a merge is ever wanted, Haunert and Wolff 2010 is the method, as its own
increment.

### 5.3 The band in plan, not over the DEM

Ruling: **2-D**. Reasons:

- The height accuracy is the refinement's job, and it is kept: every DEM
  node inside the mesh is checked against its triangle whatever the borders
  are, and refinement puts points back on a border where the terrain needs
  them (increment 20b inserts on constraint segments). Measured: at 10 m
  vertical tolerance the 50 m band's mesh is 291 932 triangles at 15°
  against 268 175 with no land cover at all (+9 %); the terrain decides.
- A 3-D band could only keep more vertices than the refinement asks for.
- CORINE's accuracy is planar (100 m positional).
- On the mesh without height refinement a 3-D band would buy height, but
  the probe showed that borders are not the fix for that error: the 50 m
  band's minimal mesh has 30 266 triangles and a largest height error of
  779 m, against 150 622 and 777 m unsimplified.
- It keeps the vector kernel free of the raster.

### 5.4 The start angle with land cover: measured

`--start-min-angle` (θ), German case, the 50 m band, 50 m vertical
tolerance (`meshes.txt`):

| θ | triangles | under 1° | under 10° | worst angle |
|---|---|---|---|---|
| 0° | 43 288 | 0.18 % | 6.57 % | 0.069° |
| **15°** | **53 177** | 0.05 % | 2.83 % | 0.473° |
| 20° | 61 364 | 0.04 % | 2.36 % | 0.473° |
| 25° (default) | 72 125 | 0.02 % | 1.87 % | 0.574° |
| 25°, gain 2 / 5 / 10 | 65 141 / 59 194 / 54 358 | 0.09 / 0.11 / 0.11 % | 2.05 / 2.70 / 3.83 % | 0.239° |
| no land cover, 0° / 15° / 25° | 23 684 / 23 726 / 23 912 | 0.00 % | 0.98 / 0.90 / 0.85 % | 1.18° |

Without land cover the start pass costs 1 %; with land cover it is all at
the borders: +67 % at 25°, +23 % at 15°. 15° keeps most of the gain in
slivers (under 10°: 6.57 % to 2.83 %, against 1.87 % at 25°) for a third of
the cost. Raising the quality gain (increment 20c) at 25° leaves more
slivers than lowering θ for about the same triangle count (gain 10: 54 358
triangles, 0.11 % under 1°, 3.83 % under 10°; θ 15°: 53 177, 0.05 %,
2.83 %). At 20 m: 117 127 at 15° against 130 613 at 25°; at
10 m: 291 932 against 298 786.

**Default (ruled, questions 1 and 4): 15° when the land-cover stage ran
with the band above 0, 25° otherwise; an explicit `--start-min-angle`
always wins.** So every mesh without land cover, every mesh with
`--features-tolerance 0`, and every stored benchmark baseline without land
cover stays bit for bit. The start angle only acts with `--tolerance`: the
start pass is part of the refinement, and `--start-min-angle` without
`--tolerance` is refused today (`src_python/tin_engine/cli.py@b39426c0:930-933`).
Measured on this one case: the scale it assumes is a 30 m DEM and CORINE;
the largest input checked is this case.

## 6. The blueprint

### Data flow

```
cli.py mesh
  | --features-tolerance BAND (default 50 m; 0 = off)  -> FeatureRequest.tolerance_m
  v
feature_input._clean(cover, region, ask)                       [Python, shapely]
  1. clip each coded polygon to   ask.tolerance_m > 0 ? the domain polygon
                                                      : the read region (today)
  2. coverage_clean (repair, today)  3. merge by class (today)
  4. ask.tolerance_m > 0: border_simplify.simplify_borders(polys, band)   (new;
     replaces shapely.coverage_simplify: feature_input.py line 434 at b39426c0)
  5. snap_to_outline (today)  -> lines and label polygons, as today
                                    |
border_simplify.simplify_borders(polygons, band) -> BorderResult      [Python]
  polygons -> flat (N, 2) float64 points + ring starts (+ which ring is
  which part's shell or hole, which part which polygon, kept in Python)
                                    |
_core.simplify_borders(points, ring_starts, band) -> BorderOutcome   [binding,
                                    |                                 GIL released]
terrain::vector_simplify::simplify_borders<K>(...)                   [C++]
  a. vertices identified by exact coordinates; edges counted by ring use
  b. junctions; borders (chains) between them; each ring as a list of
     (border, direction)
  c. the collapse loop of section 3 over all borders at once
  d. rings rebuilt from their borders: same rings, same order, same
     orientation
                                    |
  <- points + ring starts; Python rebuilds the same polygons, parts, holes
```

The calls, by name (today's at `src_python/tin_engine/feature_input.py@b39426c0:415-434`):

- **Step 1, band above 0:** `_polygonal(shapely.intersection(g,
  self.domain.polygon))` per coded polygon: the domain polygon `_add` and
  `snap_to_outline` already use, overlay in floating point, **no
  `grid_size`**; `_polygonal` keeps the polygon parts and drops the lines and
  points an intersection leaves where a polygon only touches the outline
  (None if nothing is left, counted in `outside`, as today). **Band 0:**
  today's `shapely.intersection(g, region)`, unchanged.
- **Step 2:** today's `shapely.coverage_clean(polys,
  snapping_distance=ask.repair_m, gap_width=ask.repair_m,
  merge_strategy="min_area")`, where `ask.repair_m` is `--features-repair`
  (default 1 m).
- **Step 3:** today's `shapely.coverage_union_all` per class.

Step 1's domain clip, when the band is on, makes the domain outline the
coverage's outer boundary; its edges are used by one ring, so they are fixed,
and every class's area **inside the domain** is kept exactly. Clipping to
the read region instead would let a collapse straddle the outline and move
area across it. With the band at 0, step 1 is today's and the start angle is
25° (section 5.4), so `--features-tolerance 0` keeps today's mesh bit for bit
(test 14). The outline rule (step 5) still moves borders within 5 m of the
outline onto it, and reports the area it moved (`area_changed`), as today.

### The C++ interface

New header `include/terrain/vector_simplify/border_collapse.hpp`, beside
`area_collapse.hpp`, whose `detail::segment_distance` and `detail::EdgeGrid`
it reuses as they are (all borders' vertices live in one array, so the
grid's vertex ids are global). `reduce_ring` is not changed in this
increment.

```cpp
namespace terrain::vector_simplify {
enum class BorderStatus : std::uint8_t { Ok, InvalidBand, BadRings };
struct BorderCounts {
    std::size_t junctions{}, borders{}, fixed_borders{};
    std::size_t collinear{}, collapses{};
    std::size_t rejected_crossing{}, rejected_side{};
};
struct BorderOutcome {
    std::vector<Point2> points;            // every ring, open, back to back
    std::vector<std::uint64_t> ring_starts; // ring k is [starts[k], starts[k+1])
    BorderStatus status{BorderStatus::Ok};
    BorderCounts counts{};
};
template <pred::GeometryKernel K>
[[nodiscard]] BorderOutcome simplify_borders(std::span<const Point2> points,
                                             std::span<const std::uint64_t> ring_starts,
                                             double band);
}
```

- Input rings are open (first vertex not repeated), any orientation, at
  least three vertices; `ring_starts` has one entry per ring plus the end.
  `BadRings` otherwise (and for a non-finite coordinate); `InvalidBand` for a
  negative or non-finite band. On a refusal, `points` and `ring_starts` are
  empty.
- Output: the same number of rings in the same order, each with the same
  orientation; a ring may start at a different vertex (its first vertex may
  have been collapsed).
- `band == 0` returns the input unchanged, bit for bit.
- Pure: no state outside the call, thread-safe, deterministic (heap ordered
  by deviation, then vertex id).
- The edge rule: an edge used by exactly two rings is a shared border edge;
  an edge used by one ring (the outline) or by three or more (a broken
  coverage) is **fixed**. A junction is a vertex with other than two distinct
  incident edges, or whose two edges are not used by the same rings. So an
  input that is not edge-matched to the bit is simplified less, never wrongly.

### The Python interface

New module `src_python/tin_engine/border_simplify.py`, no I/O:

```python
@dataclass(frozen=True, slots=True)
class BorderResult:
    polygons: tuple[BaseGeometry, ...]   # same count and order as the input
    counts: BorderCounts                 # the binding's counts

def simplify_borders(polygons: Sequence[BaseGeometry], band_m: float) -> BorderResult: ...
```

It raises `ValueError` on a refusal status, with the status in plain words.

### What replaces `--features-tolerance`

Nothing beside it: **the flag keeps its name and changes its meaning**.
Today: "simplify land-cover borders by this much", a GEOS area threshold's
square root, default 0. Now: the band, a distance no border moves past, each
class's area kept, **default 50**. `shapely.coverage_simplify` is no longer
called. The band is measured from the borders as repaired and clipped
(steps 1 to 3), not from the file's: the repair (`--features-repair`, 1 m)
and the outline rule (`--features-outline-snap`, 5 m) move borders on top
of it, as today. The help text becomes "Metres: simplify land-cover borders,
each moved at most this far from its repaired, clipped border and each class
keeping its area; the repair and the outline rule move borders on top of
this. 0 is off. Default: 50."; the record's `features_tolerance_m` sentence
becomes "Land-cover borders simplified, each at most this far from its
repaired, clipped border, each class's area kept (0 = off)". The record and
`--stats` key `features_tolerance_m` stay. Metres in the mesh's projected
CRS, as every other land-cover flag (a geographic DEM already needs
`--out-crs`).

The band applies wherever today's flag did: to every source whose class map
has codes (`corine`, `corine-water`), the polygons the land-cover stage
cleans. 50 m is half CORINE's positional accuracy; it assumes CORINE
(1:100 000) and was checked on the German case (1 195 km², 151 579
vertices). A future source with its own accuracy (NMD, 10 m) would carry its
own default; not this increment.

The start angle: `src_python/tin_engine/cli.py@b39426c0:988` picks `DEFAULT_START_MIN_ANGLE`; it
picks a new `LANDCOVER_START_MIN_ANGLE = 15.0` instead when the land-cover
stage ran (`FeatureSet.cover_vertices is not None`), the band is above 0
(`--features-tolerance`), and no `--start-min-angle` was given. The record's
`start_min_angle_deg` already says which was used. The `--start-min-angle`
help gains "Default: 25, or 15 with land cover simplified
(--features-tolerance above 0)."

### Assessment against the architect's criteria

Data and execution apart: the band is a field of the frozen
`FeatureRequest`; the kernel takes a number. State: the kernel is a pure
function, the binding releases the GIL. Dependencies: none new (shapely and
the kernel already present; CLAUDE.md §2 untouched). Async: one synchronous
call inside feature reading, which already runs off the event loop.

## 7. Guarantees and the checks that enforce them

1. **Every class keeps its area; so does every polygon part**, up to the
   rounding of E (measured: 6e-13 relative at 100 m). Mechanism: section 3.
   Before the outline rule; the outline rule then changes area by what it
   reports.
2. **The band, two-sided.** For every border, the Hausdorff distance between
   its source polyline γ and its simplification is at most the band ε.
   Enforced by **anchors**: every simplified vertex V carries an anchor F, a
   point of γ, with |V − F| ≤ ε; the anchors are in order along γ; an
   original vertex is its own anchor; junctions are fixed and their own
   anchors. An anchor is a position on γ, (segment index i, parameter t in
   [0, 1]), ordered by i, then t. A collapse A-B-C-D to A-E-D is allowed
   only if E's anchor F_E (the point nearest E on one source segment, tried
   for each segment from the one holding F_A to the one holding F_D) has
   **F_A ≤ F_E ≤ F_D in that order**, including on the segments that hold
   F_A and F_D (there, t at least F_A's, or at most F_D's), |E − F_E| ≤ ε,
   every source vertex between F_A and F_E is within ε of segment A-E, and
   every one between F_E and F_D within ε of E-D. A candidate segment whose
   nearest point breaks the order is skipped, as the prototype does
   (`docs/benchmarks/2026-10-08/clc-simplify/design-probe/border_apsc.py@29d7346e:133`).
   The cost per candidate is today's: one pass over the range with prefix
   and suffix maxima.

   *Why that bounds Hausdorff.* By the order condition, the anchors of a
   border's simplified vertices are in order along γ at every step, so
   consecutive anchors cut γ into parts that follow each other and cover
   it. Take an edge S = V-W with anchors F_V ≤ F_W
   and γ' the part of γ from F_V to F_W. (a) γ' lies within ε of S: its
   inner vertices by the check, F_V and F_W because they are within ε of V
   and W, and every point of a segment of γ' because the distance to a
   segment is convex. γ is the union of such parts, so γ is within ε of the
   simplification. (b) Every point P of S is within ε of γ'. Measure x along
   S from V (x = 0) to W (x = L). If γ' has points with x ≤ x_P and with
   x ≥ x_P, it is connected, so some X on it has x_X = x_P; X is within ε of
   S and 0 ≤ x_P ≤ L, so its nearest point on S is P, and |X − P| ≤ ε.
   Otherwise all of γ' lies beyond P, say x > x_P (the other side is the
   same with W): then 0 ≤ x_P < x_{F_V}, and with h the offset of F_V from
   the line, |F_V − P|² = (x_{F_V} − x_P)² + h² ≤ x_{F_V}² + h² =
   |F_V − V|² ≤ ε². ∎
3. **Topology.** The output is a valid coverage with the same polygons,
   parts and holes, every part in the same face relations: each border
   keeps the same two neighbours, each junction the same borders around it,
   each island its enclosing face. Mechanism: tests (i) and (ii) of section
   3. With both, no other edge can reach into the swept region: entering it
   means crossing A-E or E-D (the old edges A-B, B-C, C-D are crossed by
   nothing, the input being planar), or ending inside it, or leaving a
   junction A or D into it and then doing one of those. So the sweep from
   A-B-C-D to A-E-D passes over nothing.
4. **Fixed:** junctions, and the coverage's outer boundary (the domain
   outline when the band is on), at their input coordinates bit for bit,
   except that −0.0 and 0.0 are one vertex, output at its first occurrence's
   bits (section 14, pin 4).
5. **Not guaranteed:** the fewest vertices (greedy; the problem is hard,
   section 2); that the swept regions are small (the band bounds distance,
   not relabelled area, which section 4 measures).

## 8. Degeneracies

- Exactly collinear runs inside a border: dropped first, as today (distance
  0, area unchanged); never at a junction.
- A border of two or three vertices: no collapse is possible (a collapse
  needs four); it stays.
- A closed border (an island): collapses stop at four vertices, as
  `reduce_ring`'s do; the collinear pass alone may leave three (section 14,
  pin 6). A loop through one junction is an open border from the junction
  back to it, so the junction never moves.
- A lens (two borders between the same two junctions): each keeps at least
  three vertices; the area keeps the lens open, tests (i) and (ii) keep the
  two apart.
- E exactly on another edge or vertex: `classify` says Touching or
  Overlapping, the collapse is refused (exact predicates).
- E undefined (A-B or C-D parallel to A-D): that line is skipped; both
  parallel: the midpoint foot, as today.
- Coordinates −0.0 and 0.0 are the same vertex.
- Polygons from `coverage_clean` that are still not edge-matched to the bit:
  their unmatched edges are fixed (section 6).

## 9. Tests `@tester` writes red first

**C++, `tests/cpp/unit/test_border_collapse.cpp` — the invariant-critical
suite of this increment** (README, cost constraints). Mutation targets its
kill record must cover: M1 drop test (ii), the swept-region side test; M2
test crossings against the border's own edges only; M3 let a junction be B
or C; M4 replace the anchored check by `deviation_of`'s (E against any
segment in the range); M5 place E on the wrong side of the equal-area line
(sign of the area flipped); M6 drop the anchor-order condition (accept F_E
before F_A or after F_D).

1. Two squares sharing a zig-zag border: both areas kept (relative 1e-12);
   the border within the band both ways (oracle: points sampled every 1/1000
   of each polyline's length plus every vertex, distance to the other line;
   a lower estimate of the Hausdorff distance, so it cannot be red on a
   correct output); fewer vertices than the input.
2. **The overshoot case** (kills M4): the CORINE border saved as
   `docs/benchmarks/2026-10-08/clc-simplify/design-probe/overshoot_fixture.json`
   (177 vertices, local origin, 1 cm; `iso.py` there shows today's check
   gives 58.7 m against a 50 m band on it alone, the anchored check 48.3 m),
   closed into two polygons by a frame far outside the band; band 50 m;
   both directions within 50 m.
2b. **The anchor order** (kills M6): a border that runs out and back, its
   two arms closer than the band, so that after earlier collapses some E's
   nearest source point lies on the other arm, before F_A or after F_D.
   With the order dropped, the collapse passes the per-collapse check and
   leaves a stretch of γ further than the band from the output; with it,
   both directions within the band (test 1's oracle). `@tester` builds the
   chain and shows the mutant red before relying on it. *As built:* M6
   changes the output on this chain but stays within the band; the test is
   kept as a band test of the out-and-back shape, and M6 is recorded
   unkilled (section 14).
3. Three polygons meeting at a junction, and a junction on the outer
   boundary: junctions unchanged bit for bit (kills M3); outer-boundary
   edges unchanged.
4. An island inside a polygon: the island's ring and the hole's ring come
   out as the same points in opposite order; both areas kept; at least four
   vertices.
5. A small island next to a border with a wide bulge, placed so that the
   bulge's collapse sweeps over the whole island without crossing it: the
   island stays in its face (kills M1).
6. Two borders of different polygon pairs running side by side closer than
   the band: no crossing in the output (kills M2).
7. A junction A where another border leaves into the region the collapse
   would sweep: refused.
8. Area sign: an asymmetric chain where E on the wrong side changes the
   area (kills M5; test 1's area check may already).
9. Refusals: band −1, NaN, ∞ → `InvalidBand`; starts not increasing, a ring
   under three vertices, last start not the point count, a NaN coordinate →
   `BadRings`; empty output on refusal.
10. Band 0: output equals input bit for bit. Same input twice: same bits.
11. Degeneracies of section 8: a collinear run, two- and three-vertex
    borders, a lens, an edge used three times (fixed), two rings that differ
    by 1e-9 on a shared border (their edges fixed, no crash).

**Python, `tests/python/test_border_simplify.py`** (the binding and the
adapter, shapely as the oracle):

12. On a small coverage (hand-made, or a committed CORINE clip if one is
    in the test data): `shapely.coverage_is_valid` after; same number of
    polygons, parts, holes; each part's area within 1e-9 relative; each
    class's boundary Hausdorff to its source (`shapely.hausdorff_distance(a,
    b, densify=0.001)`: every segment cut into 1 000, a lower estimate of
    the true distance, so it cannot be red on a correct output) at most the
    band plus 1e-6 m; the set of class pairs that
    share a border unchanged; every junction point present in the output.
13. `simplify_borders` raises `ValueError` naming the refusal.

**CLI, amendments to `tests/python/test_cli_features_cleanup.py` and the
land-cover CLI suites:**

14. The default: record and `--stats` say `features_tolerance_m` 50.
    **Band 0 is today's mesh bit for bit:** `test_cli_features_cleanup.py`'s
    `run` (the `bumpy` DEM, the `square` domain, `--tolerance 1`) with
    `--features` (the `corine` fixture) `--features-map corine --features-tolerance
    0`; the reference is the SHA-256 of the output's points and triangles
    (read back from the `.vtk`, not the file's bytes, whose header may carry
    record text), computed by `@tester` with the same command at
    `b39426c0` and written into the test as a constant, with that commit
    named beside it; and the record's `start_min_angle_deg` is 25.
15. With land cover, band above 0 and no `--start-min-angle`:
    `start_min_angle_deg` 15; with land cover and `--features-tolerance 0`:
    25; without features: 25; `--start-min-angle 25` with land cover: 25.
16. Existing tests that assume the old default (0) or the old meaning
    (Visvalingam-Whyatt) are updated in the red commit, each with the reason
    in the message.

## 10. Size, split point, speed

Real sizes to scale from: `reduce_ring`'s green commit `ac778e43` counted
403 lines (`python3 tools/count_loc.py ac778e43^ ac778e43`): 286 in
`area_collapse.hpp`, 46 binding, 22 stub, 49 Python. The prototype's
collapse loop is 128 Python lines plus 25 for the anchored check, its
border extraction 54, its rebuild 22.

| part | counted lines (estimate) |
|---|---|
| `border_collapse.hpp`: extraction of junctions and borders, ring rebuild | 80-100 |
| `border_collapse.hpp`: the collapse loop over many borders, anchored band, tests (i) and (ii) | 230-280 |
| binding and `_core.pyi` | 60-75 |
| `border_simplify.py` (polygons to arrays and back) | 30-45 |
| `feature_input._clean`, `cli.py` (flag text, default, start angle), record text | 20-30 |
| **total** | **420-530** |

Under 700 in one PR. **Split point** if the red suite pushes the estimate
past 600: PR A, the kernel, binding, stub and `border_simplify.py` (nothing
the user sees changes); PR B, the wiring, the 50 m default and the 15° start
angle.

**Speed, per step, on the German CORINE case with the default flags** (band
50 m, repair 1 m). Measured by
`docs/benchmarks/2026-10-08/clc-simplify/design-probe/clip_probe.py`
(output `clip_probe.txt` there; Mac on battery, best of three for the
clips), in EPSG:25832, with the read region approximated as the domain's
convex hull plus 100 m:

| step | today | with this design | basis |
|---|---|---|---|
| 1. clip (1 881 CORINE polygons) | 0.12 s, to the read region | 0.12 s, to the domain polygon (516 vertices) | measured |
| 2. `coverage_clean` | 0.35 s (221 147 vertices in) | 0.23 s (152 639 in) | measured |
| 3. merge by class | unchanged | unchanged | |
| 4. simplify | 0.12 s, `coverage_simplify` (step 1's probe) | 0.2-0.6 s, estimated | 60 274 collapses × 3.3 µs, `reduce_ring`'s measured time per collapse (121 177 collapses on the same rings one by one, 0.40 s), times 1 to 3 for test (ii) against every nearby vertex instead of one keep-point |
| 4. array conversion and polygon rebuild | none | under 0.1 s, estimated | 150 000 points through `shapely.get_coordinates` and back |
| 5. outline rule | unchanged | unchanged | |

So the stage costs about what it does today (−0.24 s measured on steps 1,
2 and the dropped simplifier, +0.2 to 0.7 s estimated for the kernel and
conversion), while the mesh it feeds is 5.7 times smaller at 50 m vertical
tolerance (1.6 times at 10 m). The 15° start angle acts only with
`--tolerance` (section 5.4). **Without land cover: 0 s added, nothing
changed** (`_clean` is not called; the start angle stays 25°). With
`--features-tolerance 0`: today's steps exactly. The kernel is O(n log n) in
the border vertices. Not measured at design time (no C++ existed yet). This
is not refine or mesh code, so no `bench.py` acceptance run; the default
changes every land-cover mesh, so the main session should have `@perf` run
`tools/bench_quick.py` once, to time the new stage and see the run get
faster. **As built:** both catchment cases of the quick check
(`docs/benchmarks/quick/cases.toml`: Numedalslagen and Lagan) run with
CORINE, so both meshes are expected to change. No quick-check baseline
exists yet (`docs/benchmarks/quick/` holds only `cases.toml` and
`hotspots.toml`), so `@perf` times the base `b39426c0` and the branch back
to back, on the same power state.

## 11. Expected effect on the German fused case

Triangles (`meshes.txt`, the prototype's borders meshed by today's
`rasputin mesh`; "today" is the probe's FT 0 at the default 25°):

| vertical tolerance | today (25°) | 50 m band, 25° | **50 m band, 15° (ruled)** | 50 m band, 0° | no land cover (25°) |
|---|---|---|---|---|---|
| none (minimal) | 150 697 | 30 266 | 30 266 | 30 266 | 513 |
| 50 m | 303 470 | 72 125 | **53 177** | 43 288 | 23 912 |
| 20 m | 336 943 | 130 613 | **117 127** | 110 665 | 96 189 |
| 10 m | 462 120 | 298 786 | **291 932** | 289 200 | 267 748 |

At 50 m, 5.7 times fewer triangles than today; at 10 m, 1.6 times, the
terrain deciding. The minimal mesh's height error stays (779 m): that is
`--tolerance`'s to fix, not the borders'.

## 12. Questions for Ola

1. **The start angle with land cover.** 15° when land cover is in the mesh,
   25° otherwise (section 5.4: +23 % triangles at 50 m instead of +67 %,
   under-1° triangles 0.05 % instead of 0.02 %). Alternatives: 15° for every
   mesh (measured to cost nothing either way without land cover here, but
   changes every mesh and benchmark), or 25° everywhere. **Default: 15° with
   land cover only.**
2. **The catchment outline gets the same band check.** `reduce_ring`'s band
   can be passed by up to half a source segment
   (`docs/increments/22-auto-catchment.md`, "Measured, not guaranteed: the
   Hausdorff distance", left open for you). The anchored check of section 7
   proves the band at the same cost. Outlines would change slightly.
   **Default: yes, as its own small pull request after this one.**
3. **Lakes-only CORINE (`corine-water`) gets the band too.** **Default: yes,
   50 m, same source accuracy.**

**Ola, 2026-10-08, on the main session's summary of this section** ("Three questions, all default yes: 1. Start angle: use 15° when land cover is present and keep 25° otherwise ... 2. Outline reducer: give `reduce_ring` the same proven distance check, as its own small PR afterwards ... 3. Water-only map: the 50 m limit also applies to the lakes-only CORINE map."): **"Go for it."** *Ruled: all three defaults; build 32 now. Recorded by the main session.*

4. **The start angle when simplification is off** (raised by the design
   review, round 1). With land cover but `--features-tolerance 0`, keep
   today's 25°, so that 0 really is off and gives today's mesh bit for bit.
   **Default: 15° only when the band is above 0.**

**Ola, 2026-10-08, on the main session's question** ("With land cover
present but simplification set to 0, should the mesh keep today's 25° start
angle? Default: yes. Use 15° only when simplification is on, so 0 really
means off."): **"I'm fine with keeping 25 start angle when turning off
simpl."** *Ruled: the default. Relayed to `@architect` by the main session;
the design (sections 5.4 and 6, tests 14 and 15) is written on it.*

## 13. ROADMAP

Row 32: built (green `82c0a7d4`, 519 net production lines; mutation round
`30f939b6`, M6 recorded as not killed, section 14); next `@perf`'s short
check (section 10), then `@reviewer`'s code review.

## 14. As built

### The mutation round (`@tester`'s hand-back, copied unchanged)

#### Mutation round for 32 (`@tester`, 2026-10-08, test commit `30f939b6`, on code `82c0a7d4`)

The suite is `test_border_collapse`, run in a scratch copy from `tools/scratch_copy.py 82c0a7d4`, macOS arm64 Release. Each mutant was one or two text replacements in the copy's `include/terrain/vector_simplify/border_collapse.hpp`. Before each build the suite's object file and binary were deleted, and only `test_border_collapse` was built. Each binary ran under a 300 s timeout. After each run the header was restored, touched, and compared byte for byte with the original. Every mutant was re-run against the final test file (`30f939b6`), so all T line numbers are in `tests/cpp/unit/test_border_collapse.cpp@30f939b6`. H line numbers are in `include/terrain/vector_simplify/border_collapse.hpp@82c0a7d4`. Lines 296 to 412 of T are inside the shared oracle `check()`; the "via" column gives the test's own call.

| # | fault planted (where) | result | killed by test:line, what it checks |
|---|---|---|---|
| M1 | test (ii) dropped entirely: `side` never set, both the swept-region term and the on-new-edge term (H:403-405, wrapped as `w == no_node && (…)`) | killed | Test 5 (via T:692). T:384 ×16: island vertices change winding about rings 0 and 3, so they left their face. T:695 ×4: island vertex winding about P is 0, not 1. T:696: `rejected_side` is 0, not ≥ 1. |
| M1b | only the swept-region term dropped (`winding(loop, p[w]) != 0` removed, H:403) | killed | Same lines as M1: T:384 ×16, T:695 ×4, T:696. |
| M2 | crossings tested against the collapsing border's own edges only (`crosses = owner[u] == owner[b] && (…)`, H:398) | killed | Test 6 (via T:708). T:364 ×2: two output edges are Crossing, not Disjoint. T:709: `rejected_crossing` is 0, not ≥ 1. |
| M3 | a loop through one junction becomes a closed border, so its junction is an inner node that can be B or C (a ring with exactly one junction goes to `border_of(ids, true)`; inserted before H:276) | **survived** at `82c0a7d4`; now killed by `30f939b6` | New case "3. a loop through one junction keeps it" (T:639, via T:655). T:336 ×2: J = (50, 0) is missing from the hole ring and from the island ring. T:402: the per-border cut no longer starts at J. It survived because no case had a loop through one junction: test 4's island has no junction, test 5's islands have none, and the lens has two. |
| M3b | junctions with exactly three edges do not cut borders, so they can be B or C (`junction[id(r)] && degree[id(r)] != 3`, H:264) | killed | Test 3 (via T:618): T:336 ×3 (junctions missing), T:346 ×7 (outline edges gone), T:399. Also tests 1, 1-scale, 2, 2b, 5, 7, 8 and four test-11 cases (12 of 19 cases fail). |
| M4 | E's distance taken to the nearest of any segment in the range, not to its anchor F_E: `std::hypot(e − f)` replaced by the minimum distance to every segment from F_A's to F_D's (H:99-113) | killed | Test 2 (via T:583). T:407 ×2: source to output 58.32 m against the 50 m band. That matches the design's 58.3 m for this probe. |
| M5 | E on the wrong side of the equal-area line (`area = -(cross(vb, vc) + cross(vc, vd))`, H:330) | killed | Test 8 (via T:741). T:329 ×2: area changed by 2.545 m² against an allowed 2.4e-9 m². T:329 also fails in 8 other cases; T:674, T:696, T:709, T:864 fail too (12 of 19 cases fail). |
| M6 | the anchor-order condition dropped (H:110-111 removed: F_E may lie before F_A on F_A's segment, or after F_D on F_D's) | **survived** | Discussed below the table. No test added. |

**M6.** It is not equivalent in output: the output changes, and the suite runs 5227 assertions on it against 5203 on the code. But I found no input where it breaks the band.

- **Fuzz evidence.** A scratch fuzz of out-and-back borders like test 2b used three generator settings with 5 000 + 20 000 + 20 000 seeds: gap 5 to 25 m, offsets up to 0.49 of the gap, band 0.5 to 6.5 times the gap, and long segments in one setting. M6 changed the output on 2 199 of the 45 000. None went past the band; the worst was 0.99981 of the band, the same as the code. The fuzz could fail: with the queue limit planted at 1.5 times the band, it flagged 117 of 2 000.
- **Why one reversal is safe.** Dropping the order lets F_E sit only on F_A's own segment before F_A, or on F_D's own segment after F_D. The piece of γ between F_E and F_A is then a straight piece of one segment with no source vertex in it. It lies within ε of A-E, because its points are interpolations of the two pairs (A, F_A) and (E, F_E), each pair within ε. The piece from F_E to F_D is covered by the suffix check plus convexity. Section 7.2's step (b) needs only a connected part from one anchor to the other, not the order. So each reversed collapse keeps the band.
- **What I could not exclude by argument.** The design's covering argument does rely on the order across later collapses. If three or more output vertices end up anchored in reversed order on one source segment, a later collapse's range can omit a straight piece of that segment that only the removed edges covered. The fuzz never produced such a case.

So my verdict is equivalent for every oracle this suite has and for 45 000 fuzzed chains, but not proven equivalent. If @architect wants it killed, a test would need a long source segment carrying three or more placed vertices. I have not built one.

### Rulings on the mutation round (`@architect`)

- **M3.** M3 and M3b together are the design's M3 (a junction moved into
  a border's middle): M3b moves a three-edge junction, M3 the junction of a
  loop through one junction. Both are killed at `30f939b6`.
- **M6** (the anchor-order condition dropped) stays in the code as defence
  in depth and is recorded as **not killed**, with `@tester`'s reason: the
  output changes, but the band held on all 45 000 fuzzed out-and-back
  chains (worst 0.99981 of the band, the same as the code), and the
  covering argument across later collapses is not proven without the
  order. No new test is asked for. Test 2b stays as a band test of the
  out-and-back shape (section 9).

### Rulings on the green step's pins (`@developer`)

The green step pinned nine behaviours the design left open. All nine are
accepted as the design's:

1. `ring_starts[0]` must be 0; otherwise the call refuses with `BadRings`.
2. An edge used twice by the same ring is fixed (never collapsed).
3. A zero-length edge makes its vertex a junction.
4. Each vertex is output at its first occurrence's coordinates. So −0.0 and
   0.0 in one place are one vertex, and at a band above 0 the output may
   carry the other occurrence's sign bit; at band 0 the input bits are kept
   (section 7, guarantee 4).
5. The fallback placement is used only when both A-B and C-D are parallel
   to A-D.
6. The collinear pass may reduce an open border to two vertices; a closed
   border keeps at least three (section 8).
7. A collapse that fails both test (i) and test (ii) counts under
   `rejected_crossing`.
8. The adapter skips empty parts, and returns a polygon with no parts
   unchanged.
9. `xy_points` names its caller in its refusals.

### As-built facts

- Green step `82c0a7d4`: 519 net production lines by
  `python3 tools/count_loc.py b39426c0 82c0a7d4`. `border_collapse.hpp`
  387 (section 10 estimated 310-380), binding 50, `_core.pyi` 33,
  `border_simplify.py` 35, `cli.py` 9, `feature_input.py` 2,
  `run_record.py` 3. Inside section 10's 420-530 total; no split.
- `ctest` 931 of 931 at `82c0a7d4` (one more case at `30f939b6`, the
  loop-through-one-junction case); ASan and UBSan clean; the suite passes
  built with `-ffp-contract=off`; `pytest` 6 979 passed.
- The CI TSan job (`.github/workflows/main.yaml`) builds and runs
  `test_border_collapse`, which starts its own threads to check purity.
- Speed: not yet measured; `@perf`'s short check is section 10's last
  paragraph.

### Left for later (code review round 1's non-blocking observations)

- A rejected collapse is not queued again unless a later collapse beside it
  re-evaluates its node (`border_collapse.hpp`'s main loop), so a collapse
  blocked only for a while is lost.
- The edge grid's cell side is the larger of the band and the mean edge
  length, and an edge is filed in every cell of its bounding box, so a long
  diagonal fixed edge fills about (length / 50 m)² cells at the default
  band: worth watching at São Francisco scale.

## Review

32 design review, 2026-10-08, @reviewer (b39426c0..a4271af8): CHANGES REQUESTED -- (1) only one of the steps the design adds has a cost on the CORINE case: the new domain clip has no time or estimate, and nothing says what the no-land-cover case costs (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:384-385,408-411,668-676); (2) coverage_clean and the new clip are not named with their method and tolerance (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:384-386, /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@b39426c0:416,422-424), and neither is the densify value of test 12 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:630); (3) the anchored check never says the anchors must stay in order (F_A <= F_E <= F_D), which the proof needs and the prototype enforces (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:525-529 against /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/benchmarks/2026-10-08/clc-simplify/design-probe/border_apsc.py@a4271af8:132-133); (4) "--features-tolerance 0 gives today's mesh bit for bit" is false with --tolerance, because 15° is picked whenever the land-cover stage ran, including at 0 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:411-413,499-502,638-639); (5) the new help and record text, "move borders at most this far", leaves out the 1 m repair and the 5 m outline rule, and the ROADMAP row and Status still say "proposed" and "three questions for Ola" after the rulings (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:3-7,484-487; /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/ROADMAP.md@a4271af8 row 32). Checked and sound: the section 7.2 Hausdorff proof (given ordered anchors), both neighbours' areas kept per collapse, the crossing and sweep argument, tests and mutation targets M1-M5, the size basis (403 at ac778e43), the overshoot fixture (58.7 m vs 48.26 m), citations, rulings, literature limits. Suggestions: iso.py should read the committed overshoot_fixture.json; "the anchored band check is this design's" could read "written for this design"; add mutation target M6 (drop the anchor-order condition).

32 design review round 2, 2026-10-08, @reviewer (247be543..36b0cb90): APPROVED. All five blockers from round 1 are closed. (1) Cost for every step on the German CORINE case: the domain clip 0.12 s against today's 0.12 s; `coverage_clean` 0.35 s to 0.23 s; the simplifier estimated at 0.2-0.6 s from `reduce_ring`'s measured 3.3 µs per collapse; nothing added without land cover (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@36b0cb90:721-745, /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/benchmarks/2026-10-08/clc-simplify/design-probe/clip_probe.txt@36b0cb90:1-5). (2) The clip, `coverage_clean` and `coverage_union_all` calls named with their arguments, matching today's code (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@b39426c0:415-434); test 12's `densify=0.001`. (3) Anchor order F_A ≤ F_E ≤ F_D stated and used by the proof, matching the prototype (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/benchmarks/2026-10-08/clc-simplify/design-probe/border_apsc.py@29d7346e:133); M6 and test 2b cover it. (4) 15° needs a band above 0; test 14 pins band 0 with a SHA-256 reference computed at b39426c0. (5) Help and record text mention the 1 m repair and 5 m outline rule; Status and ROADMAP row current. Ola's ruling recorded word for word (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@36b0cb90:785-795). check_citations exit 0. Suggestion taken in this commit: section 11's header "(proposed)" now "(ruled)".

32 code review round 1, 2026-10-08, @reviewer (d1ebdc6d..552b54e5; 519 counted lines by count_loc.py b39426c0 82c0a7d4, same at HEAD, inside the 420-530 estimate, the 600-line split not triggered): CHANGES REQUESTED -- four comment or docstring claims the change made false: (1) "No ``_core``." (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:20) and "never imports _core" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/project_structure.md@552b54e5:170), but feature_input now imports border_simplify, which calls _core.simplify_borders (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:42,437); (2) "(the CLI's are 1, on, 0, 5)" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:135), but the band default is now 50 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/cli.py@552b54e5:772); (3) cover_vertices counted "after the clip to the read region" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:151-152), but with the band on (the default) the clip is to the domain (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:417); (4) "clipped to the domain as linework, never as areas (R6)" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:6-8), but land-cover polygons are now cut by the domain polygon as areas (:417-418). Checked and sound: the anchored check against the section 7 proof (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@552b54e5:88-118); the area-keeping placement of E (:329-346); a junction never B or C (:326); the crossing and side tests (:384-411); the kill table (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tests/cpp/unit/test_border_collapse.cpp@30f939b6); band 0 equals b39426c0's mesh; the 15/25 degree rule (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/cli.py@552b54e5:983-988); the TSan job (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/.github/workflows/main.yaml@552b54e5:162,178). Runs: ctest 932/932; test_border_collapse under TSan clean; pytest 6979 passed; real-data suites passed (only master's old test_16c failed, by design); ruff, mypy and the governance gates clean; a 300-coverage fuzz (about 220 000 collapses) found no failure. CI: not pushed. @perf's short timing check still to run. Suggestions: list border_collapse.hpp and border_simplify.py in project_structure.md and put the three red notes in the past tense naming 44f25968; rejected collapses are never retried; the edge grid fills (length/50 m)² cells for long diagonal fixed edges (watch at São Francisco scale); Status/ROADMAP order (@perf's check comes after this review).

32 code review round 2, 2026-10-08, @reviewer (67f4f0a3..29fbe35c; count_loc.py b39426c0 HEAD still reports 519 counted lines; `git diff 552b54e5 HEAD -- include bindings` is empty, so no production code changed except comments): CHANGES REQUESTED -- one blocker. A fourth red-step note, added in the red commit 44f25968, is still in the present tense and is now false: "RED at the commit that adds them: the stage calls ``coverage_simplify`` and clips to the read region" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tests/python/test_feature_repair.py@29fbe35c:50-51). Nothing under src_python calls coverage_simplify any more, and the clip is to the domain when the band is above 0 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@29fbe35c:418). Round 1's four blockers are closed, each checked against the code (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@29fbe35c:8, :21, :136, :152-153); round 1's suggestions taken (project_structure.md entries true against /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@29fbe35c:205, :350-351 and /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/bindings/core.cpp@29fbe35c:1408). ruff, check_citations and merge-tree against origin/master clean. CI: not pushed. @perf's short timing check still to run.
