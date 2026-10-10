# Increment 32: land-cover borders simplified within a band, each class's area kept

**Status:** built. Designed on `0b5e0d4c`, Ola's four rulings (section
12), design review round 2 approved. Red step `44f25968` (`@tester`), green
step `82c0a7d4` (`@developer`, 519 net production lines), mutation round
`30f939b6` (`@tester`, one new test case); `@architect`'s rulings on both
hand-backs are in section 14. Code review round 1 done: changes requested
(claims in comments, docstrings and `project_structure.md` that the change
made false), fixed by `@developer`'s `4d4f0411` and `@architect`'s docs
commit after it; red-step notes put in the past tense (`29fbe35c`). Code
review round 2 asked for one more red note (Review section), fixed by
`633808b9`; round 3 approved. `@perf`'s short timing check (`35f220b9`)
found the runs about 70 % slower and the worst angle far lower; cause and
fix in section 15, designed on the defaults of question 5 (section 12),
decided on default while Ola was away, reversible; fix design approved
in fix-design review round 7 (`8a8c2488`). Fix built: red `2fba4349`,
green `05c8f73a` (573 counted lines), test fix `8605f607`; mutation round
M7 to M12 killed (M12b and M13 ruled in section 15.6). `@perf`'s quick
check `b192a2c6` met section 15.5 except Numedalslagen's worst angle,
0.311° against a 0.4° floor; section 15.6 finds that triangle made by the
height refinement, not by the borders' clearance, and restates the floor
(question 6, decided on default while Ola was away, reversible). Test 19b
(`179ad71f`, `@tester`) kills M13. Fix code review round 1 (`e33c0e65`):
production code sound, changes requested in prose and in the test file
(red-step staging left in `test_border_collapse.cpp`, this Status,
sections 13 and 14, the ROADMAP row, `project_structure.md`); fixed in
`fa91c266` (prose) and `2a46bca8` (staging removed). Fix code review round 2
asked for the kill table's copy to be made verbatim (recording commit).
Fix code review round 3 approved. Next: the push on Ola's yes. Not refine
or mesh code, so no `bench.py` acceptance run.

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

*Superseded for the order of the steps by section 15.2* (the outline rule
now runs before the domain clip and the simplifier); the reasons below for
the domain clip still hold.

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
    **As built (`8269e0a8`, after CI):** the digest is architecture-bound
    (arm64 Mac `29963a44…`, x86_64 Linux `9cd16850…`, same shapely and
    numpy), so the constant is one value per platform, with its provenance;
    the Linux value was taken at `d5b2ba4c`, since no Linux run of
    `b39426c0` exists, so there it guards against change from `d5b2ba4c` on.
    A second, platform-independent test shows the simplifier is never
    called at band 0 (a spy on `simplify_borders`; band 2 is the control).
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

5. **The order of the outline rule and the simplifier** (raised by
   `@perf`'s quick check, section 15). Today the land cover is clipped to
   the catchment, simplified, then the 5 m outline rule runs; the rule was
   built for borders that cross the outline, and on borders that end on it
   it takes 4 times as long and reports about 80 times the area it really
   moves. Proposed: run the outline rule first, on the same input as
   before 32, then clip to the catchment and simplify. Effect: the rule's
   time and area exactly as before 32; the 50 m band then counts from the
   borders after the rule, not before it. The other way, keeping the order
   and reworking the outline rule for borders that end on the outline, is
   larger and touched one polygon by about 1 km² in the probe. The price
   of running the rule first: the 50 m simplifier runs after it and can
   bring a border back within 5 m of the outline, which the rule had just
   cleared. On Numedalslagen 7 vertices of 51 349 came back within 5 m,
   measured without a clearance (`probe.txt`, `clearR.py`'s lines); with
   the 1 m clearance of section 15.2 none can come closer than 1 m. Accept those, or keep them out with a second, 5 m clearance
   against the outline only (more code, fewer collapses near the outline).
   **Default: run the outline rule first, and accept the few vertices back
   within 5 m.** Section 15 is designed on it.

   *Decided on default while Ola was away, reversible: both parts (outline
   rule first; the vertices back within 5 m accepted). The main session
   took the defaults, 2026-10-08; recorded by `@architect`.*

6. **The worst-angle floor** (raised by `@perf`'s quick check of the
   fix, section 15.6). I expected the worst angle back at 0.4° or more;
   Numedalslagen has one triangle at 0.311°, the next at 0.487°. It is not
   made by the land-cover borders being too close: it is a height
   refinement point 5.2 m beside a 1 km straight border, which the
   simplifier makes long by design. Splitting long borders every 500 m
   brings the triangles under 1° from 76 to 13, but the worst becomes a
   0.326° triangle elsewhere that touches no border closely. Proposed:
   drop the 0.4° expectation; the borders' promise is the 1 m clearance,
   and the worst single angle is decided by the height refinement, as it
   was before 32 (0.83° there was not a bound either). **Default: drop the
   floor and keep the fix as built; the 500 m split, if wanted, as its own
   small change later.**

   *Decided on default while Ola was away, reversible. The main session
   took the default (drop the 0.4° expectation, keep the fix as built),
   2026-10-09; recorded by `@architect`.*

## 13. ROADMAP

Row 32: built (green `82c0a7d4`, 519 net production lines; mutation round
`30f939b6`, M6 recorded as not killed, section 14); code review approved;
`@perf`'s short check found it slower (section 15). The fix of section 15
is done: red `2fba4349`, green `05c8f73a` (573 counted lines in all),
mutation round M7 to M13 killed (M12b equivalent), quick check `b192a2c6`
(Numedalslagen −11 %, Lagan −1.3 %, outline rule and moved area as base,
triangles −42 %; worst angle restated, question 6); red-step staging
removed (`2a46bca8`); fix code review round 3 approved. Next: the push.

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

#### Mutation round for the fix of section 15 (`@tester`; table rows copied unchanged, the paragraphs around them are the main session's summary)

Summary by the main session of @tester's method (2026-10-09; test file at 8605f607, kernel code at 05c8f73a; scratch copy from tools/scratch_copy.py 8605f607, macOS arm64 Release; per mutant: text replacement in border_collapse.hpp, suite object and binary deleted, only test_border_collapse rebuilt, 300 s limit, JUnit report, header restored and checked byte for byte). Baseline: all 29 cases pass (8459 assertions). T lines are in tests/cpp/unit/test_border_collapse.cpp@8605f607; T:1157 is inside check_clearance, T:491 inside check().

| # | fault planted (where) | result | killed by test:line, what it checks |
|---|---|---|---|
| M7 | vertex-to-new-edge check dropped (the `near = near \|\| (p[w] != A && …A, E…) \|\| (p[w] != D && …E, D…)` line in the grid scan) | killed | Test 18, F18: T:1309 `rejected_clearance` ≥ 1, T:1311 `collapses` 0, T:1312 output unchanged. Also 20b F20a: T:1157 ×2 (clearance oracle at clearance 2), T:1420 (clearance-2 output differs). |
| M8 | E-to-edge check dropped (`near = near \|\| segment_distance(E, p[u], p[v]) < clearance`) | killed | Test 19, F19: T:1373 `rejected_clearance` ≥ 1, T:1374 `collapses` 0, T:1375 output unchanged. |
| M9 | short-edge check at placement dropped (the `consider` condition wrapped as `false && (…)`) | killed | Test 17, F17b: T:1279 `skipped_placements` ≥ 1, T:1281 `collapses` 1, T:1282 one vertex fewer, T:1284 one new vertex. Also F17a T:1244, T:1245 and 18b T:1346. |
| M10 | A and D excluded from the vertex check by node number (`w != a`, `w != d`), not by coordinate | killed | 20b S4: T:1408 `rejected_clearance` 0 at the floor, T:1409 floor output = clearance 0's, T:1414 and T:1415 the same at clearance 1. Also F20a and F20b T:1408, T:1409, T:1414, T:1415, and F17b T:1280, T:1281, T:1282, T:1284. |
| M11 | clearance doubled where read (`clearance *= 2.0` after validation) | killed | 20b F20a: T:1408, T:1409, T:1414, T:1415. 20b F20b: T:1407 `skipped_placements` 0 at the floor, T:1409, T:1413, T:1415. |
| M12 | A's and D's coordinates skipped for both new edges, and the explicit A-to-E-D and D-to-A-E line dropped | killed | Test 18b, F12: T:1347 `rejected_clearance` ≥ 1, T:1348 `collapses` 0, T:1349 output unchanged. |
| M1 | test (ii) dropped entirely: `side` never set (wrapped as `w == no_node && (…)`); re-run because the green step moved the A/D exclusion in this block | killed | Test 5: T:491 ×16 (winding numbers change), T:802 ×4 (island winding about P not 1), T:803 `rejected_side` ≥ 1. |
| M1b | only the swept-region term dropped (`winding(loop, p[w]) != 0`) | killed | Same lines as M1: T:491 ×16, T:802 ×4, T:803. |
| M12b (extra) | only the explicit A-to-E-D and D-to-A-E line dropped (`bool near = false;`) | survived (equivalent, argued) | Survives because the scan already checks A's copies in other borders against E-D, and D's copies against A-E. A junction always has a copy in another border, and that copy's edge starts at A, which is inside the query box, so the grid always returns it. I argued this, not proved it. |
| M13 (extra) | the query box not grown by the clearance | survived | No fixture puts a near vertex in a grid cell that only the grown box reaches. Cells are at least the band wide (50 m here); growing the box by 1 m changes the set of cells only next to a cell boundary. |

Not re-run: M2 to M6 (crossing test, junction rule, anchored check, placement formula unchanged by the green step: git diff 2fba4349 05c8f73a -- include/).
Test 10 fix 8605f607: workers call simplify_borders<DefaultKernel> directly; ASan+UBSan 10/10, TSan halt_on_error=1 10/10 (+3/3 test 10 alone); the pre-fix file reproduces the Catch2 race under TSan (then hangs).

Follow-up (@tester, 179ad71f; summary by the main session, the table row below copied unchanged): test 19b added (S4's collapse, band 76.5, cell boundary x = 76.5 is 0.3 m beyond the box from x = 76.8; island triangle (76.3, 14.5), (70.3, 17.5), (70.3, 11.5), first vertex 0.5 m from A-E and E-D). On 05c8f73a's kernel: 30 cases, 8856 assertions pass; TSan halt_on_error=1 exit 0.

| # | fault planted (where) | result | killed by test:line, what it checks |
|---|---|---|---|
| M13 | the query box not grown by the clearance (both `lo`/`hi` lines in the pop, `q.x ∓ clearance` and `q.y ∓ clearance` made `q.x` and `q.y`), in a scratch copy of HEAD with the new test file | killed | Test 19b: T:1446 `rejected_clearance` ≥ 1, T:1447 `collapses` 0, T:1448 output unchanged (`tests/cpp/unit/test_border_collapse.cpp@179ad71f`). The other 29 cases pass. |

@architect (621d46e2) confirmed M12b equivalent against the built scan (border_collapse.hpp@05c8f73a): the scan reaches A as the end of its own incoming edge P-A, or through another border's or a fixed edge at a junction; D through D-X and its copies.

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
- Speed: `@perf`'s first quick check (`35f220b9`) found the runs about
  70 % slower; section 15 has the cause and the fix. The fix's quick
  check (`b192a2c6`, `docs/benchmarks/2026-10-09/clc-simplify/quick/`,
  base `b39426c0` against `8605f607`): Numedalslagen 9.18 s against
  10.34 s (−11 %), Lagan −1.3 %; the outline rule at base's time, the
  moved area exactly base's; output triangles −42 % on Numedalslagen
  (653 039 against 1 118 006); worst angle 0.311° on Numedalslagen and
  0.49° on Lagan (section 15.6, question 6).

### Left for later (code review round 1's non-blocking observations)

- A rejected collapse is not queued again unless a later collapse beside it
  re-evaluates its node (`border_collapse.hpp`'s main loop), so a collapse
  blocked only for a while is lost.
- The edge grid's cell side is the larger of the band and the mean edge
  length, and an edge is filed in every cell of its bounding box, so a long
  diagonal fixed edge fills about (length / 50 m)² cells at the default
  band: worth watching at São Francisco scale.

## 15. The quick check's finding: cause and fix (`@architect`)

`@perf`'s quick check (`35f220b9`,
`docs/benchmarks/2026-10-08/clc-simplify/quick/README.md`) found both
catchment runs about 70 % slower, the time in the 5 m outline rule
(`snap_to_outline`, 2.76 s to 12.18 s under the profiler on Numedalslagen),
the area the rule reports as moved 80 to 95 times larger, and the worst
angle down from 0.83° to 0.0157° (Numedalslagen) and from 0.73° to 0.00066°
(Lagan). One probe on Numedalslagen, run from the scripts in
`docs/benchmarks/2026-10-08/clc-simplify/quick/fix-probe/` (output, with
the commands, in `probe.txt` there), captured the rule's input three ways:
band off (base `b39426c0`'s input), band on, and band on with the simplifier
replaced by "return the input". Lagan was not probed; the fix below is
expected to act there the same way, and `@perf`'s re-run checks it.

### 15.1 Causes

Three separate causes. The suspicion in the brief holds for the time, not
for the area or the worst angle.

1. **Time: the domain clip, not the simplifier.** With the band on, step 1
   clips land cover to the domain polygon (section 6), so every class ring
   that reaches the outline now runs *along* it: 14 844 ring edges lie on
   the outline, against none at base. The rule treats each as an edge
   within 5 m: it cuts it every 5 m, places and joins every piece in
   Python. And where a border comes within 5 m of the outline, the rule
   now pinches the part against its own outline stretch, so the rebuilt
   shell touches itself: 245 rebuilt shells invalid (202 "Ring
   Self-intersection"), each through `make_valid`, against 3 at base
   (`probe.txt`, steps 2 and 6). The rule on the input clipped but not
   simplified takes 14.0 s, on the simplified one 10.3 s, at base 2.3 s
   (step 2): the simplifier makes it faster, the clip slower.
2. **Moved area: the report is wrong, not the borders.** The area the rule
   really moves, each polygon inside the outline before against after, is
   13 573 m² lost and 14 073 m² gained on the branch, against 13 238 m² at
   base (step 4). The reported 1 097 391 m² comes from `_loops`
   (`src_python/tin_engine/feature_input.py@633808b9:683-703`): for a ring
   with exactly one kept input vertex it builds one loop from the whole
   old ring and the whole new ring; where the two coincide along the
   outline, `make_valid` returns the whole part as the loop. 76 loops over
   1 000 m² make up 1 083 899 m², each the area of a whole part that is
   still in place (the largest, 119 080 m², is 119 043 m² still covered
   after the rule; steps 3 and 4). Simplified rings near the outline have
   few interior vertices, so more of them keep one; the clip alone already
   reports 94 709 m² (7 times base) the same way. It is a defect of the
   rule's bookkeeping that the domain clip exposes; base's input does not
   reach it (its report equals the real change, 13 238 m²).
3. **Worst angle: the simplifier, not the outline rule.** The eight worst
   triangles lie 0.8 to 11.9 km from the outline (step 5). Each has two
   constraint vertices 7 to 32 mm apart, and in the four checked one is a
   junction (three classes meet, fixed) and the other a point E the
   simplifier placed on a border leaving that junction, so the new edge
   A-E is about 1 cm long. Tests (i) and (ii) refuse exact touches only;
   nothing keeps E, or a new edge, any distance from what it does not
   touch. The repaired, clipped source has no vertex within 1 m of an edge
   it is not an end of (the repair's promise, `--features-repair`); the
   simplified output has 12, the closest 5.6 mm (`clear.py`). Same cause
   for Lagan expected, not checked.

### 15.2 Ruling (on the default of question 5)

**Order.** The outline rule moves back to where it was built to act: on the
coverage clipped to the read region, before the domain clip and the
simplifier. `_clean`'s steps become:

```
band 0 (unchanged, today's mesh bit for bit):
  1. clip to the read region  2. coverage_clean  3. merge  4. snap_to_outline -> lines, polygons
band above 0:
  1. clip to the read region  2. coverage_clean  3. merge   (as band 0)
  4. snap_to_outline, as today, on that input        -> polygons, area_changed (lines unused)
  5. clip each polygon to the domain polygon (_polygonal(shapely.intersection(g, domain)))
  6. simplify_borders(polygons, band, clearance=ask.repair_m)
  7. lines: each simplified ring cut where it lies on the outline (_outline_lines)
```

Steps 1 to 4 are then base's, so the rule gets base's input, takes base's
time, moves base's area and reports it correctly, by construction. The
domain clip (step 5) keeps section 6's reason: the outline is the coverage's
outer boundary, fixed in the simplifier, so every class keeps its area
inside the domain. The band is now measured from the borders as repaired,
put to the outline and clipped.

**What the order gives up** (question 5, second part, decided on default).
The simplifier runs after the rule, so it can bring a border back within
5 m of the outline: 0 such vertices after the rule and the clip, 7 of
51 349 after simplifying, without the clearance (`probe.txt`, `clearR.py`'s
lines). With the clearance, none of them can come closer than 1 m to the
outline's edges. They are accepted, and the help text says so (15.4).

**Lines (step 7).** New private helper in `feature_input.py`,
`_outline_lines(polygon, outline) -> tuple[LineString, ...]`: for every ring
of every part, an edge lies on the outline when **one** outline segment is
within `IN_LINE` (1 µm) of both its ends (an `STRtree` of the outline's
segments, query `dwithin IN_LINE`, then the two end distances to that
segment); the ring is then cut with today's `_chains(xy, on)`, which drops
those edges. Same output shape as `OutlineSnap.lines`. A gap border (a class
against no class) is not on the outline and stays a line, as today. Probe
(`orderprobe.py`): 88 876 ring edges, 16 700 on the outline, 2 570 lines,
0.21 s.

**Clearance (the worst angle).** `simplify_borders` gains a clearance:

```cpp
template <pred::GeometryKernel K>
[[nodiscard]] BorderOutcome simplify_borders(std::span<const Point2> points,
                                             std::span<const std::uint64_t> ring_starts,
                                             double band, double clearance = 0.0);
```

- `BorderStatus::InvalidClearance` for a negative or non-finite clearance;
  `BorderCounts::rejected_clearance` counts the refusals at pop, and
  `BorderCounts::skipped_placements` the placements left out at placement
  (tests 17 and 20b read both). A placement is counted where it is
  skipped, inside `consider`, after the finite check and **before** the
  anchored deviation and the band check (the band is checked on the best
  placement only, when it is queued), so a skipped placement counts
  whether or not its deviation is within the band, and every evaluation
  of a candidate counts again.
- **At placement** (`candidate`): a placement with |E − A| or |E − D| under
  the clearance is not considered (the other line's placement still is).
- **At pop**, beside tests (i) and (ii), over the grid cells of the box of
  A, B, C, D, E grown by the clearance: refused if a vertex lies closer
  than the clearance to a new edge it does not end. Edge by edge: against
  **A-E**, every vertex but B and C (by node) and any vertex at A's
  coordinates; against **E-D**, every vertex but B and C and any vertex
  at D's coordinates. So D is checked against A-E and A against E-D.
  A's and D's copies are excluded **by coordinate**, each only from the
  edge it ends, as test (ii) excludes them
  (`include/terrain/vector_simplify/border_collapse.hpp@633808b9:400`):
  the junction's copies in the other borders are other nodes at the same
  point. Also refused if E lies closer than the clearance to any edge
  other than A-B, B-C, C-D (edges ending at A or D included, E ends none
  of them). These cover every new approach: the distance between a new
  edge and an old edge is reached at an end of one of them, and of those
  ends A and D are old (their distances to old edges do not change),
  while E and the old edge's ends are checked. B and C need no copies
  rule: a vertex with a copy in another border has other edges, so it is
  a junction and is never B or C. Distances in
  floating point (`detail::segment_distance`): a refusal criterion, so
  rounding can only refuse a collapse or let one through at the clearance
  to rounding; the exact tests (i) and (ii) stay as they are and still
  decide topology. **A clearance hit does not end the scan**: the scan
  over the cells goes on as today and ends only on a failure of (i) or
  (ii), so a collapse failing (i) is counted under `rejected_crossing`,
  else one failing (ii) under `rejected_side`, and only one failing
  neither but the clearance under `rejected_clearance`. The existing
  counts keep their meaning.
- Clearance 0 gives today's output bit for bit (every existing test).
- Binding: keyword `clearance: float = 0.0`; `border_simplify.simplify_borders(
  polygons, band_m, clearance_m=0.0)`; `_clean` passes `ask.repair_m`.
- Why the repair distance: the repair already promises borders at least
  that far apart (measured: none closer on the source), and the simplifier
  only refuses what would break that promise. Scale assumed: metres, 1 m
  default, CORINE at 1:100 000; checked at Numedalslagen (147 429 vertices
  in). Guarantees 1 to 4 of section 7 are untouched: a refusal only keeps
  vertices. Guarantee added: **no vertex or edge the simplifier creates
  comes closer than the clearance to a vertex or edge it does not share an
  end with.** What it does not promise: vertices the outline rule placed
  (310 under 1 m after the rule and clip, the closest 2.6 mm, base's too,
  `clearR.py`) stay as they are.

**Not in this fix.** `_loops`' one-kept-vertex report (cause 2) and the rule's
self-touching rebuild on borders that end on the outline are defects of the
outline rule (20c-3), out of 32's reach once the rule gets base's input
again. Recorded for a later increment; the first fix tried for the time
(no cuts on edges lying on the outline) changed one polygon by 961 767 m²
through that rebuild (`fixprobe.py`, `rebuild_diff.py`), so the rule is not
to be patched for the clipped input here.

### 15.3 Tests `@tester` writes red first

C++, `tests/cpp/unit/test_border_collapse.cpp`. Every fixture below is
two 100 m squares, left `[(0,0), A, B, C, D, (0,100)]` and right
`[A, (200,0), (200,100), D, C, B]`, A = (100, 0) and D = (100, 100) being
junctions (three distinct edges each), so exactly one collapse, of B and C,
is possible; some add a triangle island (its ring and the hole's ring in
the square that holds it). Each was measured at clearance 0 with the built
`_core` (the kernel as at `633808b9`):
`docs/benchmarks/2026-10-08/clc-simplify/quick/fix-probe/fixtures/fixtures32.py`,
output `fixtures32.txt` there. The two placements (on line A-B, on line
C-D) and their anchored deviations are re-enacted in `fx.py` there, and the
kernel's new vertex equals the re-enacted winner in every case that
collapses. The quantities of 15.2 are taken over every vertex and edge, a
superset of what the grid looks at, so a floor they meet holds for the
kernel. Figures in metres.

| fixture | B, C | island | winner E (deviation) | other E (deviation) | clearance 0: collapses | winner: \|E−A\|, \|E−D\|, other vertex to new edge, A to E-D, D to A-E, E to edge | other: \|E−A\|, \|E−D\| |
|---|---|---|---|---|---|---|---|
| F17a | (102, 24), (99, 40) | none | (100.04, 0.48) (1.969) | (100.04, 102.4) (2.400) | 1 | 0.482, 99.52, 99.96, 0.482, 99.52, 0.480 | 102.4, 2.400 |
| F17b | (104.1, 52.4), (99.5, 4.9) | none | (99.9629, −0.4742) (4.118) | (99.9629, 92.9436) (4.121) | 0 (`rejected_crossing` 1) | 0.476, 100.47, 99.96, 0.037, 100.00, 0.474 | 92.94, 7.057 |
| S4 | (84, 10), (84, 55) | none | (76.8, 14.5) (7.200) | (76.8, 34.75) (7.754) | 1 | 27.36, 88.59, 78.16, 27.36, 88.59, 14.50 | 41.78, 69.25 |
| F18 | as S4 | (95.24, 10.64), (88.45, 14.88), (88.66, 7.67) | as S4 | as S4 | 1 | 27.36, 88.59, **0.494**, 27.36, 88.59, 11.66 | 41.78, 69.25 |
| F19 | as S4 | (75.44, 19.34), (70.4, 13.37), (77.18, 9.49) | as S4 | as S4 | 1 | 27.36, 88.59, 2.580, 27.36, 88.59, **0.497** | 41.78, 69.25 |
| F20a | as S4 | (95.77, 11.49), (88.98, 15.73), (89.19, 8.52) | as S4 | as S4 | 1 | 27.36, 88.59, **1.496**, 27.36, 88.59, 12.21 | 41.78, 69.25 |
| F20b | (116.5, 45), (84.5, 55) | none | (100.55, 1.5) (16.193) | (100.55, 101.597) (16.256) | 1 | **1.598**, 98.50, 99.46, **1.598**, 98.50, **1.500** | 101.6, **1.689** |
| F12 (band 60, own frame, below) | (0, −1.25), (−20, 0) | none | (0, −0.25) (1.000) | (50, 0) (50.000) | 1 | 0.750, 100.0, 49.50, 0.750, 100.0, 0.750 | 50.00, 150.0 |

F12 is the fold of fix-design review round 6, in its own frame: A = (0, 0.5)
and D = (−100, 0), each a junction of three polygons,
`[A, (0,50), (−50,50)]`, `[A, (−50,50), (−100,50), D, C, B]` and
`[A, B, C, D, (−100,−50), (100,−50), (100,50), (0,50)]` (a valid coverage,
checked in `fixtures32.py`). Its other placement E = (50, 0) has deviation
exactly 50.000, so the band is 60. As a collapse it crosses and touches
nothing but at A and D, and its quantities are: other vertices 49.50 from
the new edges, **A 0.500 from E-D**, D 100.0 from A-E, E 50.00 from every
edge. The "A to E-D" and "D to A-E" columns are the two distances section
15.2 now checks; in F17b the other placement, made at clearance 1, has
them at 92.94 and 7.057. `fixtures32.py` prints the git blob of
`border_collapse.hpp` it ran against (`dc879c84`, the blob at `633808b9`,
unchanged since `82c0a7d4`); the `_core` in the `.venv` was built from it.

In F17a the other placement crosses the top edge; in F17b the winner lies
below the bottom edge (so the kernel refuses it by test (i) today) and the
other placement, as a collapse, crosses and touches nothing but at A and D
(checked with shapely in `fixtures32.py`) and is within the band (4.121).
The islands lie inside one square, cross neither chain, and have no vertex
in the swept region (checked there too).

17. **Short edge at a junction.**
    - *F17a.* At clearance 0 the collapse makes E 0.482 m from the junction
      A. At clearance 1: no output edge shorter than 1 m that the input
      did not have; `skipped_placements` ≥ 1 (the winner); the output
      equals the input (the other placement crosses the top edge:
      `rejected_crossing` 1).
    - *F17b* (kills M9). At clearance 0 no collapse (the winner crosses).
      At clearance 1 the winner (|E − A| 0.476) is skipped at placement
      and the collapse is made with the other placement: |E − D| 7.057,
      A 92.94 from E-D, D 7.057 from A-E, E 7.056 from every other edge,
      other vertices 99.999 from the new edges: `skipped_placements` ≥ 1,
      `rejected_clearance` 0, the border has one vertex fewer, and the new
      vertex is (99.9629…, 92.9436…) as the re-enacted placement gives it
      (compared to 1e-9). Under M9 the winner is kept, refused at removal
      by test (i), and the output equals the input: red. The two
      deviations differ by 0.003 m; the winner's identity is fixed by the
      kernel's arithmetic, deterministic, and the test asserts it at
      clearance 0 first (the kernel's `rejected_crossing` 1 there).
18. **A vertex near a new edge** (F18, kills M7). An island vertex 0.494 m
    from the new edge A-E; E is 11.66 m from every edge it does not end
    (at least 1 m, so only the vertex check can refuse). Clearance 1:
    refused (`rejected_clearance` ≥ 1), the output equals the input.
    Clearance 0: made.
18b. **A junction near the far new edge** (F12, band 60, kills M12). At
    clearance 0 the collapse is made with the A-B placement, E = (0, −0.25).
    At clearance 1 that placement is skipped (|E − A| 0.750,
    `skipped_placements` ≥ 1); the C-D placement E = (50, 0) is queued and
    refused at removal, because A lies 0.500 m from E-D
    (`rejected_clearance` ≥ 1); the output equals the input. Under M12 A is
    excluded from both new edges, every other quantity is at least 49.50,
    so the collapse is made with E = (50, 0), E-D passing 0.5 m from the
    junction (a 0.6° corner): red.
19. **E near an edge** (F19, kills M8). The island's edge 0.497 m from E;
    every island vertex at least 2.580 m from A-E and E-D (at least 1 m, so
    only the E-to-edge check can refuse). Clearance 1: refused, output
    equals input. Clearance 0: made.
20b. **Room to spare** (fails a check that refuses too much). Each case
    pins its shape first, with correct code, and asserts the pin before
    the result. *The pin:* in the clearance-0 run, the quantities of 15.2
    (|E − A| and |E − D| at every placement considered; the distance from
    every vertex the check looks at to A-E and E-D, and from E to every
    edge it looks at, at every removal that passes (i) and (ii)) are at
    least a stated floor. *Shown by:* a run at the floor as clearance,
    which must report `skipped_placements` 0 and `rejected_clearance` 0
    (the refusals are strict, so nothing fires there exactly when every
    quantity is at least the floor, and the run then takes the
    clearance-0 run's steps); `@tester` prints both counts and the
    output's equality with the clearance-0 output, bit for bit. The table
    above gives each floor's margin.
    - **2 m to spare, at a junction** (S4, kills M10). Floor 2 m (the
      smallest quantity is 14.50). Result at clearance 1: the output
      equals the clearance-0 output bit for bit. M10 excludes A by node
      number, so the junction's copies in the outer borders sit at 0 from
      A-E: refused, red.
    - **Inside the over-refusal band** (kills M11). F20a: an island vertex
      1.496 m from the new edge A-E. F20b: E 1.500 m from the bottom edge
      leaving the junction A, |E − A| 1.598, the other placement's
      |E − D| 1.689. Floor 1.2 m (in rows F20a and F20b of the table,
      every quantity is at least 1.496). Also a run at clearance 2, whose output must differ from
      the clearance-0 output: at clearance 2 a check fires that changes
      the output (in F20a the removal is refused; in F20b both placements
      are skipped). Result at clearance 1: the output equals the
      clearance-0 output bit for bit, `rejected_clearance` and
      `skipped_placements` 0. Correct code passes by the pin. M11 doubles
      the clearance where it is read, so at clearance 1 every check
      behaves as correct code at clearance 2, whose output differs from
      the clearance-0 output: red.
20. **Refusals and identity.** Clearance −1, NaN, ∞ → `InvalidClearance`,
    empty output; clearance 0 equals the three-argument call bit for bit on
    every existing case.

Mutation targets: M7 drop the vertex-to-new-edge check (test 18); M8
drop the E-to-edge check (test 19); M9 drop the short-edge check at
placement (test 17, F17b); M10 exclude A and D from the vertex check by
node number instead of by coordinate (the junction's copies in the other
borders then sit at distance 0 from A-E and E-D; test 20b, S4); M11 the
clearance doubled where the kernel reads it, so every check uses twice
the value given (test 20b, F20a and F20b); M12 skip A's and D's
coordinates for both new edges instead of each only for its own edge
(test 18b, F12).

Python:

21. `tests/python/test_border_simplify.py`: on both coverages of test 12,
    `hand_made()` and `corine_coverage()` (the committed CORINE extract,
    EPSG:3035), with `clearance_m=1`: every output vertex not in the input
    is at least 1 m less 1e-6 m from every edge it is not an end of. The
    slack is 1e-6, not 1e-9: at EPSG:3035 coordinates near (4.8e6, 5.4e6) m
    one unit of floating-point precision is 9.3e-10 m, so 1e-9 would leave
    about one unit for two distance computations (the kernel's and
    shapely's) that may round differently. (Round 1's 5 % vertex-count
    clause stays dropped: it depended on the coverage's shape.)
22. `tests/python/test_cli_features_cleanup.py` (test 14's `run`, `corine`
    fixture): the record's `land_cover_area_moved_m2` with the default band
    equals the band-0 run's exactly (the rule sees the same input); no
    land-cover line has an edge lying on the domain outline (one outline
    segment within 1 µm of both ends); test 14's band-0 hash unchanged.

### 15.4 Change for `@developer`

`_clean` reordered as in 15.2 (the band-0 path untouched), `_outline_lines`,
the clearance through kernel, binding, stub and adapter, the help and record
text updated. `--features-tolerance` help: "Metres: simplify land-cover
borders, each moved at most this far from its border after the repair and
the outline rule, each class keeping its area; the simplification can bring
a border back within the outline-snap distance of the outline, but not
closer than the repair distance. 0 is off. Default: 50." Record sentence:
"Land-cover borders simplified, each at most this far from its border after
the repair and the outline rule, each class's area kept (0 = off)".

Estimate 60 to 90 counted lines (kernel 30 to 45, binding and stub 8,
Python 25 to 35): 519 becomes about 580 to 610.

**Split or not.** Section 10 set a split at an estimate past 600 (PR A the
kernel, binding and adapter; PR B the wiring), so that the first PR changes
nothing the user sees. This estimate reaches 610. Ruling: **one PR**. The
split point was for code not yet written; the wiring is now built, reviewed
and is what the fix changes, so splitting now separates reviewed code from
its fix and adds a review round without making either PR easier to review.
The margin to the 700 ceiling is 90 lines at the top of the estimate.
**Stop line:** if `count_loc.py b39426c0 <green>` passes 640, `@developer`
hands back before committing more, and the split is then PR A = kernel,
binding, stub, `border_simplify.py` and the clearance; PR B = `_clean`,
`_outline_lines`, `cli.py`, record text.

### 15.5 Expected effect

From the probe's single runs (not `bench.py`, another agent running).
The rows marked "(bench)" take base and branch from `bench_quick.py` and add
the differences the probe measured, so they mix two kinds of run; read them
as estimates. The triangle and angle rows are expectations for the code
**with** the clearance, which no run has had (no C++ was built for the
fix); the probe's 91 132 vertices are without it.

| Numedalslagen | base | branch now | after the fix (expected) |
|---|---|---|---|
| outline rule | 2.27 s | 10.28 s | 2.24 s (same input as base) |
| domain clip, after the rule | - | - | 0.80 s |
| simplifier | - | 0.66 s | 0.63 s, a little more with the clearance |
| lines from the simplified rings | - | - | 0.21 s |
| clean-up phase (bench) | 5.19 s | 13.97 s | about 6.8 s (base + 1.6 s) |
| whole run (bench) | 9.92 s | 16.95 s | about 10.7 s (about +8 %; the refine and land-cover savings, 0.8 s, kept) |
| area moved, reported | 13 238 m² | 1 097 391 m² | 13 238 m², exactly base's |
| worst angle | 0.832° | 0.0157° | back near base's: the eight worst are all simplifier-made pairs under 4 cm; expected at or above 0.4° (the German probe's 0.473°) |
| output triangles | 1 118 006 | 654 920 | about the same as now (91 132 simplified polygon vertices against 87 066) |

Lagan, by the same reasoning: the outline rule back to base's time and
area (16 050 m²), the clean-up about base's 5.65 s plus 1.5 to 2 s. So the
whole run is expected within about 10 % of base, not faster: the design's
"the stage costs about what it does today" (section 10) holds within that,
and the mesh is about 40 % smaller. If `@perf`'s re-run shows the worst
angle under 0.4° on either case, or the area moved different from base's,
it comes back to `@architect`.

### 15.6 As built: the worst angle, and two mutants (`@architect`)

Probe: `docs/benchmarks/2026-10-09/clc-simplify/quick/angle-probe/`
(`probe.txt` has the commands and output), the fix as built (`.venv`'s
`_core` rebuilt from HEAD by `@perf`), Numedalslagen, the quick check's
flags. The mesh is 653 039 triangles, as `@perf`'s.

**The 0.311° triangle.** It is 6.7 km from the outline. Its vertices: a
simplifier point E on the border between polygons 13 and 16; a point P on
that border's constraint 1 038 m from E (a split point the refinement put
on the line, not a cover vertex); and the DEM node (121 300, 6 712 580),
5.209 m from segment E-P. Two of its edges are free; E-P is the
constraint. The refinement puts a DEM node within half a cell (5 m on
this 10 m grid, `delta_p`,
`include/terrain/refinement/refine_points.hpp@05c8f73a:287`) of a
constraint in as its foot on the segment; this node is 0.209 m beyond
that, so it went in as itself, and with nothing else near the 1 km
segment the triangle is a needle: 5.2 m across 1 km. So:

- *Not the simplifier's clearance:* no two cover elements are close here;
  the near point is a DEM node, which the clearance does not see. No
  clearance on the borders can prevent it.
- *Not the outline rule:* 6.7 km inside.
- *Made by:* the height refinement beside a long border. The simplifier
  makes the border long (that is its job); the 1 km edge is within the
  band.
- *The rest of the tail:* the next five (0.487° to 0.520°) are the same
  kind: one cover vertex and DEM nodes or a point on a constraint line.
  None has two cover vertices closer than 1 m.

**Whether a border change would help, measured.** `cap_probe.py` splits
every simplified edge longer than L into equal collinear pieces (the same
points on both sides of a shared edge) and meshes again. At L = 500 m:
8 274 points added, 650 407 triangles (−0.4 %), triangles under 1° from
76 to 13, worst 0.326°: a new needle 5.3 km from the outline, a source
vertex and two diagonal DEM neighbours 14.1 m apart, 200 m away, with no
constraint edge at all. So the worst single angle is the height
refinement's (no quality pass after the start pass), and moves with any
change; a border rule bounds the count under 1°, not the worst one. The
0.4° floor of 15.5 was a wrong expectation: it was drawn from the
millimetre pairs of cause 3, which the clearance removed (none left). The
floor is restated (question 6, default): **the borders' promise is the
1 m clearance, checked by tests 17 to 21; the worst angle is reported,
not bounded.** Lagan, 0.49°, is consistent with this.

**M12b** (only the explicit A-to-E-D and D-to-A-E line dropped):
**equivalent, confirmed.** The grid scan visits every edge in the cells
of the grown box except those starting at A, B or C (by node), and checks
both ends of each edge, A's coordinates only against E-D and D's only
against A-E. A always reaches that scan as an end of a visited edge: if A
is inside its border, as the end of its own incoming edge P-A (it starts
at P, not at A, so it is visited); if A is a junction, as the start of an
edge of another border or of a fixed edge at A (a junction has at least
one, and every live edge is in the grid); in both cases the edge has A in
the box, so its cell is queried. D likewise: its own outgoing edge D-X
starts at D, not A, B or C, so it is visited, and so are the edges of
D's copies. So the explicit line checks nothing the scan does not.
It stays as written, a cheap statement of the rule; no test.

**M13** (the query box not grown by the clearance): **a reachable gap,
worth one case.** The grid's cell side is the larger of the band and the
mean edge length, and its origin is the lowest corner of all points
(`include/terrain/vector_simplify/border_collapse.hpp@05c8f73a`, the
grid's construction). A vertex within the clearance of a new edge, across
a cell boundary that lies within the clearance outside the box of A, B,
C, D, E, whose own edges lie wholly beyond that boundary, is found only
by the grown box. New test **19b** (`@tester`, in the mutation round's
manner: it passes on `05c8f73a` and is red under M13): compute the cell
side and origin from the input by that rule; place an island so that one
of its vertices lies 0.5 m from a new edge across a cell boundary that
is 0.3 m beyond the box, with the island's bounding box wholly beyond
the boundary; measure the fixture at clearance 0 as in 15.3 and print
the cell side, the boundary and the distances. Expected at clearance 1:
refused (`rejected_clearance` ≥ 1), output equals input. Under M13:
made.

## Review

32 design review, 2026-10-08, @reviewer (b39426c0..a4271af8): CHANGES REQUESTED -- (1) only one of the steps the design adds has a cost on the CORINE case: the new domain clip has no time or estimate, and nothing says what the no-land-cover case costs (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:384-385,408-411,668-676); (2) coverage_clean and the new clip are not named with their method and tolerance (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:384-386, /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@b39426c0:416,422-424), and neither is the densify value of test 12 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:630); (3) the anchored check never says the anchors must stay in order (F_A <= F_E <= F_D), which the proof needs and the prototype enforces (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:525-529 against /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/benchmarks/2026-10-08/clc-simplify/design-probe/border_apsc.py@a4271af8:132-133); (4) "--features-tolerance 0 gives today's mesh bit for bit" is false with --tolerance, because 15° is picked whenever the land-cover stage ran, including at 0 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:411-413,499-502,638-639); (5) the new help and record text, "move borders at most this far", leaves out the 1 m repair and the 5 m outline rule, and the ROADMAP row and Status still say "proposed" and "three questions for Ola" after the rulings (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@a4271af8:3-7,484-487; /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/ROADMAP.md@a4271af8 row 32). Checked and sound: the section 7.2 Hausdorff proof (given ordered anchors), both neighbours' areas kept per collapse, the crossing and sweep argument, tests and mutation targets M1-M5, the size basis (403 at ac778e43), the overshoot fixture (58.7 m vs 48.26 m), citations, rulings, literature limits. Suggestions: iso.py should read the committed overshoot_fixture.json; "the anchored band check is this design's" could read "written for this design"; add mutation target M6 (drop the anchor-order condition).

32 design review round 2, 2026-10-08, @reviewer (247be543..36b0cb90): APPROVED. All five blockers from round 1 are closed. (1) Cost for every step on the German CORINE case: the domain clip 0.12 s against today's 0.12 s; `coverage_clean` 0.35 s to 0.23 s; the simplifier estimated at 0.2-0.6 s from `reduce_ring`'s measured 3.3 µs per collapse; nothing added without land cover (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@36b0cb90:721-745, /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/benchmarks/2026-10-08/clc-simplify/design-probe/clip_probe.txt@36b0cb90:1-5). (2) The clip, `coverage_clean` and `coverage_union_all` calls named with their arguments, matching today's code (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@b39426c0:415-434); test 12's `densify=0.001`. (3) Anchor order F_A ≤ F_E ≤ F_D stated and used by the proof, matching the prototype (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/benchmarks/2026-10-08/clc-simplify/design-probe/border_apsc.py@29d7346e:133); M6 and test 2b cover it. (4) 15° needs a band above 0; test 14 pins band 0 with a SHA-256 reference computed at b39426c0. (5) Help and record text mention the 1 m repair and 5 m outline rule; Status and ROADMAP row current. Ola's ruling recorded word for word (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@36b0cb90:785-795). check_citations exit 0. Suggestion taken in this commit: section 11's header "(proposed)" now "(ruled)".

32 code review round 1, 2026-10-08, @reviewer (d1ebdc6d..552b54e5; 519 counted lines by count_loc.py b39426c0 82c0a7d4, same at HEAD, inside the 420-530 estimate, the 600-line split not triggered): CHANGES REQUESTED -- four comment or docstring claims the change made false: (1) "No ``_core``." (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:20) and "never imports _core" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/project_structure.md@552b54e5:170), but feature_input now imports border_simplify, which calls _core.simplify_borders (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:42,437); (2) "(the CLI's are 1, on, 0, 5)" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:135), but the band default is now 50 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/cli.py@552b54e5:772); (3) cover_vertices counted "after the clip to the read region" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:151-152), but with the band on (the default) the clip is to the domain (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:417); (4) "clipped to the domain as linework, never as areas (R6)" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@552b54e5:6-8), but land-cover polygons are now cut by the domain polygon as areas (:417-418). Checked and sound: the anchored check against the section 7 proof (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@552b54e5:88-118); the area-keeping placement of E (:329-346); a junction never B or C (:326); the crossing and side tests (:384-411); the kill table (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tests/cpp/unit/test_border_collapse.cpp@30f939b6); band 0 equals b39426c0's mesh; the 15/25 degree rule (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/cli.py@552b54e5:983-988); the TSan job (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/.github/workflows/main.yaml@552b54e5:162,178). Runs: ctest 932/932; test_border_collapse under TSan clean; pytest 6979 passed; real-data suites passed (only master's old test_16c failed, by design); ruff, mypy and the governance gates clean; a 300-coverage fuzz (about 220 000 collapses) found no failure. CI: not pushed. @perf's short timing check still to run. Suggestions: list border_collapse.hpp and border_simplify.py in project_structure.md and put the three red notes in the past tense naming 44f25968; rejected collapses are never retried; the edge grid fills (length/50 m)² cells for long diagonal fixed edges (watch at São Francisco scale); Status/ROADMAP order (@perf's check comes after this review).

32 code review round 2, 2026-10-08, @reviewer (67f4f0a3..29fbe35c; count_loc.py b39426c0 HEAD still reports 519 counted lines; `git diff 552b54e5 HEAD -- include bindings` is empty, so no production code changed except comments): CHANGES REQUESTED -- one blocker. A fourth red-step note, added in the red commit 44f25968, is still in the present tense and is now false: "RED at the commit that adds them: the stage calls ``coverage_simplify`` and clips to the read region" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tests/python/test_feature_repair.py@29fbe35c:50-51). Nothing under src_python calls coverage_simplify any more, and the clip is to the domain when the band is above 0 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@29fbe35c:418). Round 1's four blockers are closed, each checked against the code (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@29fbe35c:8, :21, :136, :152-153); round 1's suggestions taken (project_structure.md entries true against /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@29fbe35c:205, :350-351 and /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/bindings/core.cpp@29fbe35c:1408). ruff, check_citations and merge-tree against origin/master clean. CI: not pushed. @perf's short timing check still to run.

32 code review round 3, 2026-10-08, @reviewer (29fbe35c..633808b9; no production code changed: `git diff 29fbe35c 633808b9` touches only docs/increments/32-landcover-simplify.md and tests/python/test_feature_repair.py; count_loc.py b39426c0 HEAD still reports 519 counted lines): APPROVED -- round 2's one blocker is closed. The fourth red note is now in the past tense, names the red commit 44f25968, and is true of that commit (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tests/python/test_feature_repair.py@633808b9:50-51, /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@44f25968:434). No present-tense red note added by this increment is left (all five name 44f25968); three older ones on base b39426c0 are not this increment's. check_citations exit 0; ruff clean. CI: not pushed. @perf's short timing check still to run.

32 fix-design review, 2026-10-08, @reviewer (35f220b9..33fc7996; no production code, 519 counted by count_loc.py b39426c0 33fc7996): CHANGES REQUESTED. (1) After the new order the simplifier runs after the 5 m outline rule and can pull borders back within 5 m of the outline (0 such vertices after the rule, 7 after simplifying, /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/benchmarks/2026-10-08/clc-simplify/quick/fix-probe/probe.txt@33fc7996:74-75); neither question 5 nor section 15.2 says so (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@33fc7996:817-828,1040-1046). (2) Tests 17-22 cannot fail a clearance check that refuses too much: add a part where clearance 1 still accepts collapses with room to spare (about test 12's coverage at clearance 0), and a mutant excluding the junction's copies of A and D by node number instead of by coordinate as test (ii) does (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@633808b9:400) (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@33fc7996:1064-1088). (3) The 580-610 counted-line estimate crosses section 10's split trigger ("past 600": PR A kernel, PR B wiring) and section 15.4 says "one PR" without weighing it (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@33fc7996:731-734,1094-1097). Checked and sound: all three causes reproduced by re-running the probes (same figures; times within 0.3 s); band 0 bit for bit; the clearance only refuses or narrows a placement and the exact tests decide topology; clearance 0 is today's output (strict comparisons against 0); speed table labelled as single probe runs; library calls named with tolerances. Suggestions: point section 6's old order to 15.2; exclude A and D by coordinate; a clearance hit must not stop the scan (counts order); orderprobe.py prints BorderCounts fields; say the "after the fix" triangle and angle rows lack the clearance; say the bench rows mix bench_quick and probe differences.

32 fix-design review round 2, 2026-10-08, @reviewer (4b810983..6c8bc16c; no production code; count_loc.py b39426c0 6c8bc16c gives 519): CHANGES REQUESTED. One blocker: M11 (clearance refuses at twice the given distance) is said to be killed by test 20b, but 20b's collapses all have at least 2 m to spare (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@6c8bc16c:1108-1118), so under clearance 1 M11 refuses none and survives; test 21's 5 % margin (:1132-1136) kills it only if the coverage has enough collapses with 1-2 m to spare, unmeasured. Fix: give 20b at least one collapse whose nearest thing is strictly between 1 and 2 m (e.g. 1.5 m) that must still be made at clearance 1, or redefine M11. Round 1's blockers otherwise closed: the 7-vertex fact in question 5 (:832-844), 15.2 (:1023-1028) and the help text (:1147-1152), matching /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/benchmarks/2026-10-08/clc-simplify/quick/fix-probe/probe.txt@6c8bc16c:76-77; "not closer than 1 m" holds; M10 killed by 20b's junction case; "by coordinate" matches /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@633808b9:400; one-PR ruling with the 640-line stop line sound (:1158-1169); six suggestions taken. check_citations exit 0. CI: not pushed.

32 fix-design review round 3, 2026-10-08, @reviewer (0efdc9cd..ad3a251e; no production code; count_loc.py b39426c0 ad3a251e gives 519): CHANGES REQUESTED. One blocker: "Under M11 both are refused, so both go red" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@ad3a251e:1118-1133) does not hold for every allowed shape: a too-close E is skipped at placement and the other line's placement tried (:1051-1052; /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@633808b9:323-348), so if |E − A| or |E − D| is also under 2 m, M11 can make the same collapse with the other E and survive. Fix: also require the clearance-1 output to equal the clearance-0 output bit for bit, or require |E − A| and |E − D| ≥ 2 m for both placements. Checked and sound: the four quantities match 15.2's three checks (:1049-1062); strict refusals make 1.5 m at clearance 1; (b) kills M10; mutation-target line updated (:1140-1144). check_citations exit 0. CI: not pushed.

32 fix-design review round 4, 2026-10-08, @reviewer (0a972b37..0563b7c1; no production code; count_loc.py b39426c0 0563b7c1 gives 519): CHANGES REQUESTED. One blocker: test 20b's bit-for-bit clearance-1 = clearance-0 requirement (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@0563b7c1:1130-1133) holds only if no clearance check fires anywhere in the case: collapses are taken in order of their best placement's deviation (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@633808b9:322-347), so one skipped winning placement or pop refusal elsewhere changes the output and correct code fails. Fix: pin in 20b that in the clearance-0 run the three quantities of 15.2 are at least 1.2 m at every winning placement and every pop except the one 1.5 m approach; @tester prints the minimum first. Round 3's blocker closed. Suggestion: reword "Under M11 both are refused" (:1133) as "both differ from the clearance-0 output". check_citations exit 0. CI: not pushed.

32 fix-design review round 5, 2026-10-08, @reviewer (99a43015..123924a9; no production code; count_loc.py b39426c0 123924a9 gives 519): CHANGES REQUESTED. Round 4's blocker closed (a correct run at clearance f with zero `skipped_placements` and `rejected_clearance` takes the clearance-0 run's steps). Four blockers from a full pass over tests 17-22: (1) 20b's "2 m to spare" case names test 1's two squares, which cannot meet the 2 m floor (measured with the built `_core`: a new vertex 0.14 m from the junction; zero-area zig-zag pieces put placements on A or D); use a hand-built one- or two-collapse shape such as the bulge of tests 5-7 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@123924a9:1126-1130; /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tests/cpp/unit/test_border_collapse.cpp@123924a9:450-456); (2) no test kills M9 (drop the short-edge check at placement): the removal-time E-to-edge check already refuses such an E, so M9 refuses instead of using the other placement; test 17 at clearance 1 must also require the collapse made with the other placement (`skipped_placements` ≥ 1, `rejected_clearance` 0, one vertex fewer), that placement pinned at least 1 m clear (:1099-1102, :1150-1151); (3) test 18 does not keep E at least 1 m from the other border's edges, so the E-to-edge check can refuse in place of the vertex check and M7 survives; state it as test 19 states the reverse (:1103-1109); (4) define M11 as the clearance doubled where it is read, so both checks use twice the value, else "behaves as correct code at clearance 2" is false; and "(so the 1.5 m approach decides something)" overstates what the clearance-2 run shows (:1137-1145, :1155-1157). Sound: the strict-refusal argument step by step; M10 killed by the junction case; M11 (as redefined) killed by the clearance-2 run; dropping test 21's 5 % clause loses no certain kill; tests 20 and 22 by construction; the new count fits the estimate. Wording: "at every removal that passes (i) and (ii)" (:1119-1120); drop "except one stated approach" (:1116-1117); say whether `skipped_placements` counts before or after the band check; test 21 names its coverage (EPSG:3035 coordinates leave about one ulp of the 1e-9 slack). check_citations exit 0. CI: not pushed.

32 fix-design review round 6, 2026-10-08, @reviewer (6a901df1..d6362a15; no production code; count_loc.py b39426c0 d6362a15 gives 519): CHANGES REQUESTED. One blocker, a design gap: the vertex check skips any vertex at A's or D's position for both new edges (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@d6362a15:1061-1068), but A ends A-E only and D ends E-D only, so E-D can pass closer than the clearance to A unchecked, breaking the promise at :1086-1088 and possibly test 21. Example (fx.py's formulas): A=(0,0.5), B=(0,-1.25), C=(-20,0), D=(-100,0) meets the 1 m repair promise; at clearance 1 the A-B placement is skipped (|E-A| 0.75) and the C-D placement E=(50,0) passes every check as designed while E-D passes 0.5 m from A (a 0.6° sliver). Fix: check A against E-D and D against A-E; add the fold as a fixture (band 60 or a shifted point, its deviation is exactly 50.0), mutant M12 "skip A and D for both edges", and the two distances as table columns. Checked and sound: fixtures32.py re-run identical to fixtures32.txt; table matches; floors met; F17b kills M9, F18 M7, F19 M8, S4 M10, F20a/F20b at clearance 2 M11; F17b's winner assertion robust (0.003 m gap vs ~1e-14 m rounding); test 21's 1e-6 m slack sound at EPSG:3035 (spacing 9.3e-10 m); skipped_placements counting matches consider(); check_citations exit 0. Suggestions: :1181 "at least 1.496" applies to rows F20a and F20b only; fixtures32.py prints the hash of border_collapse.hpp; fx.quantities gains the two new distances.

32 fix-design review round 7, 2026-10-08, @reviewer (84381421..8a8c2488; no production code; count_loc.py b39426c0 8a8c2488 gives 519): APPROVED. Round 6's gap closed: the vertex check excludes A's coordinates only for A-E and D's only for E-D (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@8a8c2488:1061-1078); the coverage argument holds (non-crossing segments are closest at an end of one; A and D old; E and old edges' ends checked); no other exclusion over-wide (B and C by node is right: a vertex with a copy in another border is a junction, /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@633808b9:203-206; the E check excludes only the three removed edges). F12 and test 18b hold and kill M12; earlier pins and kills unchanged; fixtures32.py re-run identical, blob dc879c84 is border_collapse.hpp at 633808b9. check_citations exit 0. CI: not pushed. Implementation note for @developer: today's scan skips the edge starting at A (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@633808b9:389), so A against E-D and D against A-E need an explicit check; 18b catches a miss. Suggestions: D in F12 lies in two polygons and on the outline, not three (:1139); M10's "at distance 0 from A-E and E-D" means A's copies from A-E and D's from E-D (:1234); test 21 checks new vertices against edges only, not new edges against old vertices; the Status line was stale (fixed in this recording commit).

32 fix code review, 2026-10-09, @reviewer (0fda2e66..179ad71f; count_loc.py b39426c0 179ad71f gives 573; the fix 54): CHANGES REQUESTED. (1) Red-step staging survives in the invariant-critical suite: concept detection, `clearance_missing`, `staged_clearance` with its dead FAIL branch, the `if constexpr (clearance_built<K>)` guard that would silently switch off test 20's clearance-0 identity check if the four-argument call stopped matching, and the "[[maybe_unused]] … until the clearance is built" notes (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tests/cpp/unit/test_border_collapse.cpp@179ad71f:145-193,198-201,1138); the header still says the file "builds before it exists" and "its kill record covers M7 to M12" (:62-63,70-71); @tester calls the API directly, removes the staging, puts the header in the past tense naming 2fba4349, adds 19b and M13. (2) Prose made false: Status (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@179ad71f:17-23), section 13 (:870-873), section 14's "Speed: not yet measured" (:949-950), ROADMAP row 32 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/ROADMAP.md@179ad71f:58). Sound: per-edge exclusion, explicit A/E-D and D/A-E checks, E against all but the three removed edges, a clearance hit not ending the scan, counts in order (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/include/terrain/vector_simplify/border_collapse.hpp@179ad71f:340-347,389-431); InvalidClearance after InvalidBand (:147-154); box grown by the clearance (:389-392); clearance 0 bit for bit; `_clean` reorder and `_outline_lines` (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/feature_input.py@179ad71f:410-459,714-734); band 0 on base's path; 15.4's wording (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/cli.py@179ad71f:768-771, /Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/src_python/tin_engine/run_record.py@179ad71f:38-41); kill record covers M7-M13, M1, M1b, M12b equivalent; quick check and 15.6 recorded faithfully. LOC 54 against 60-90, 573 against 580-610, under 640. ctest 943/943; pytest on the touched suites 182 passed; gates clean; merge-tree against origin/master clean. CI: not pushed. Suggestions: copy the fix's kill table into this file; name 2fba4349 in the three Python red notes; project_structure.md:526 does not list the clearance.

32 fix code review round 2, 2026-10-09, @reviewer (e33c0e65..2a46bca8; count_loc.py b39426c0 2a46bca8 gives 573; no production code changed since 05c8f73a): CHANGES REQUESTED. One blocker, in prose only: the fix's kill table was labelled "copied unchanged"/"copied verbatim" (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@2a46bca8:912,914,932) but was not: rows M7 and M8 lost the escaped `\|\|` (a bare pipe splits a GitHub table cell, :918-919), M12b's sentence was shortened (:926), and the two paragraphs are summaries; and the "Next:" lines (:26-28, :882-883) name done work. The copy errors were the main session's (its scratch copy); fixed in this recording commit: escapes and the M12b sentence restored, the paragraphs labelled as summaries, the Next lines updated. Round 1's two blockers closed: the staging is gone and `same_at_clearance_0` always runs (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/tests/cpp/unit/test_border_collapse.cpp@2a46bca8:143-172), the header past tense (:62-75); Status, sections 13 and 14, ROADMAP row 32 and project_structure.md true against the code; every CHECK and REQUIRE survives the staging removal. test_border_collapse 30 cases, 8856 assertions under -Werror; check_citations exit 0. CI: not pushed.

32 fix code review round 3, 2026-10-09, @reviewer (2a46bca8..84b89865, one recording commit, prose only; count_loc.py b39426c0 84b89865 gives 573; no production code changed since 05c8f73a): APPROVED. Round 2's blocker closed: all 11 rows of the fix's kill table match @tester's handbacks character for character, by exact line comparison against the session transcript; the `\|\|` escapes restored in M7 and M8 (/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify/docs/increments/32-landcover-simplify.md@84b89865:918-919), M12b's sentence restored (:926); headings and the follow-up line say only the rows are copied, the paragraphs labelled as summaries (:914,916,934); section 14's first-round "copied unchanged" table (:889) also matches its handback line for line; Next lines true (:26-29, :884-885). check_citations exit 0; git diff --check clean. CI: not pushed; must be green before the merge. Not refine or mesh code, so no bench.py run.
