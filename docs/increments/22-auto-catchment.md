# Increment 22 — auto-catchment: the catchment of a lake, from the DEM (Bygdin first)

Status: **designed** (`@architect`, 2026-09-29), on branch
`increment22-autocatchment` off master `b4847d7`. Nothing is built. Ola was
asleep while this was written; every choice he would normally make is marked
"Default (main session / @architect, 2026-09-29), for Ola to confirm", with
the alternative, so the loop can run tonight.

**Closes.** The auto-catchment row of `ROADMAP.md`: a catchment computed from
the DEM, reduced to a polygon a mesh can afford, and handed to `--domain`.
After this increment:

```sh
rasputin catchment --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --seed 8.5425 61.3512 \
    --lakes ../rasputin_data/corine2018_dtm10_utm33.gpkg --lakes-layer corine2018 \
    --out bygdin.geojson
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --domain bygdin.geojson --tolerance 1 --out bygdin.vtk
```

**Not closed.** Stream extraction and Strahler order (`auto_catchments.md`).
Flow accumulation, and with it snapping a pour point to the strongest flow
nearby. A catchment whose area runs off the data (a border with Sweden, the
sea): it is refused, not written. Holes in the catchment are filled, not kept.
Parallel flooding (Barnes 2016). Moving the outline tracer into C++ (below,
"Where it runs"). The depression case under "The window" is a known limit.

## Ola's requirements, quoted

From the ROADMAP row, 2026-09-29:

1. "we will get an extreme amount of points in the catchment polygon. This
   must be taken into account, or the number of triangles will be
   overwhelming." The traced outline is reduced before meshing, to a stated
   horizontal tolerance.
2. "the resulting polygon has approximately the same area as the fine one".
   The reduction keeps the area. Off-node vertices ("We have aleady
   established that the polygons and DEM don't need to match") make exact
   area possible.
3. The reduced polygon is simple (no self-crossings) and contains the seed.

And on what the DEM is: "Isn't the DEM points, really?" The catchment below is
a set of DEM nodes, and its outline is drawn between nodes, never around cell
squares.

Tonight's goal, verbatim: "I would really like you to work while I sleep
tonight, to get a first working example of the autocatchment. If possible I
would like the Norwegian lake Bygdin to be the seed for the calculation."

## What the data says about Bygdin (measured 2026-09-29)

These were measured before the design, because two of the facts in the brief
did not hold, and the design depends on the third.

- **The coordinates in the brief are off.** UTM33 (177934, 6818339), given as
  the dam, is on a slope at 1580 m north of Vinstre, 10 km east of Bygdin and
  outside its catchment. NVE's regulation dam, "BYGDIN REGULERINGSDAM" (dam
  1192), is at (168087, 6815337) in EPSG:25833, about 8.792 E, 61.330 N. The
  lake runs from x 142.3 km to 168.5 km, so its west end is in tile
  `6801_3`, not `6801_2`: the catchment crosses a tile edge, as 15a intended.
- **The lake is not one flat in DTM10.** Of the 399,371 DEM nodes inside
  CORINE's Bygdin polygon, 101,058 are exactly 1057.4 m; the rest scatter
  between 1055.4 and 1059.7 m (1st and 99th percentile). The largest
  connected exactly-flat region is 7.1 km² of a 40 km² lake. A seed defined as
  "the flat under the point" would miss most of the lake.
- **CORINE has the lake.** Feature `fid` 54101 of `corine2018` in
  `corine2018_dtm10_utm33.gpkg`, class 512, 39.94 km², one part, no holes.
  NVE gives 40.03 km² at the highest regulated level.
- **NVE publishes the catchment.** NVE's reservoir catchment ("delfelt") 1187,
  "BYGDIN", has `delfeltAreal_km2` = **305.54**, and its upstream list is
  itself alone, so it is the whole catchment. Source: NVE's map service,
  `https://gis3.nve.no/map/rest/services/Mapservices/VassdragsreguleringVannkraft/MapServer/8/query?where=delfeltNr%3D1187&outFields=*&f=json`,
  fetched 2026-09-29; the same query with `returnGeometry=true&outSR=25833&f=geojson`
  returns the polygon (1376 vertices, 305.539 km² by shapely). The reservoir
  record is layer 9 (`magasinNavn='BYGDIN'`, 40.03 km²). Norwegian Wikipedia's
  Bygdin article gives 305.59 km² citing NVE's REGINE units (2015); the NVE
  service is the number used here.
- **The method below was prototyped** (throwaway Python in the scratchpad, not
  committed) over a 4201 x 3101-node window of 10 m DTM10: 3,049,095 nodes,
  304.91 km², against NVE's 305.54 km² (-0.21 %); 99.1 % of NVE's nodes are
  in ours and 99.3 % of ours in NVE's. The traced outline had 17,804
  vertices. shapely's topology-preserving Douglas-Peucker, used only as a
  rough gauge, took it to 924 vertices at 20 m. The acceptance run measures
  these numbers again with the real code.

## Prior art: legacy and literature

**Legacy.** Nothing to carry over. The ROADMAP's grep still returns no files:

```
$ grep -rliE "watershed|flow.?acc|flow.?dir|pour.?point|catchment" legacy
$ echo $?
1
```

For the reduction:

```
$ grep -rliE "simplif|visvalingam|douglas|peucker" legacy
legacy/rasputin/triangulate_dem.h
legacy/rasputin/geo_tiff_reader.py
legacy/rasputin/mesh.py
legacy/rasputin/application.py
legacy/tests/test_gml_repository.py
legacy/bindings.cpp
```

Every hit is CGAL's surface-mesh edge collapse (Lindstrom-Turk cost), which
simplifies a triangle mesh, not a polygon. Nothing is carried over.

**Literature.** Each reference below was checked against Crossref on
2026-09-29 (title, authors, venue, DOI); abstracts were read where Crossref or
Semantic Scholar had them. No paper was read in full tonight, and where the
design leans on a detail beyond the abstract, it says so.

- **Barnes, Lehman and Mulla 2014a**, "Priority-Flood: An optimal
  depression-filling and watershed-labeling algorithm for digital elevation
  models", *Computers & Geosciences* 62:117-127, doi:10.1016/j.cageo.2013.04.024
  (open copy arXiv:1511.04463). Floods the DEM inward from its edges with a
  priority queue; the abstract says it "can also be adapted to label
  watersheds and determine flow directions". **This increment builds on it.**
  What differs: the flood labels one catchment (the lake's) rather than every
  edge outlet's, and ties are broken first in, first out by a push counter so
  the result is deterministic. Barnes's faster variant (a plain queue inside
  depressions) is an optimisation that can come later.
- **Metz, Mitasova and Harmon 2011**, "Efficient extraction of drainage
  networks from massive, radar-based elevation models with least cost path
  search", *HESS* 15:667-678, doi:10.5194/hess-15-667-2011. GRASS's
  `r.watershed`: routing by least-cost search instead of filling first. The
  flood here routes the same way in effect: each node drains to the node that
  flooded it.
- **O'Callaghan and Mark 1984**, "The extraction of drainage networks from
  digital elevation data", *Computer Vision, Graphics, and Image Processing*,
  doi:10.1016/S0734-189X(84)80047-X. Crossref lists it as volume 27(2), page
  247; it is commonly cited as 28(3):323-344. Cite it by DOI until someone
  resolves which. D8: each node drains to its steepest neighbour, slope
  measured with the diagonal's length. **Departure, and why:** the flood
  drains each node to its *lowest* filled neighbour, not its steepest, so the
  diagonal distance plays no part. The two differ only where a diagonal
  neighbour is lower but less steep, which moves a divide by a node here and
  there. In exchange there is one pass instead of three (fill, flats,
  directions), and no flat resolution at all (next point). D8 comes when
  stream extraction needs accumulation.
- **Garbrecht and Martz 1997**, "The assignment of drainage direction over
  flat surfaces in raster digital elevation models", *J. Hydrology*
  193:204-213, doi:10.1016/S0022-1694(96)03138-1; and **Barnes, Lehman and
  Mulla 2014b**, "An efficient assignment of drainage direction over flat
  surfaces in raster digital elevation models", *Computers & Geosciences*
  62:128-135, doi:10.1016/j.cageo.2013.01.009. Both give flats realistic flow
  paths (towards lower terrain, away from higher). **Not needed here, and
  why:** which catchment a flat belongs to depends only on where it spills,
  and the flood gets that exactly; the path water takes across the flat does
  not change the answer. The lake itself is the seed, so its surface, flat or
  noisy, is never routed at all. The one case the flow path decides is a flat
  that spills over two outlets at the same level: the flood splits it by
  distance in steps from each outlet. That is recorded as a limit, not
  solved.
- **Wang and Liu 2006** (*IJGIS* 20:193-213, doi:10.1080/13658810500433453)
  and **Lindsay 2016** (*Hydrological Processes* 30:846-857,
  doi:10.1002/hyp.10648): filling and breaching. Priority-Flood fills; breach
  is out of scope, as `auto_catchments.md` has it.
- **Kong and Rosenfeld 1989**, "Digital topology: Introduction and survey",
  *CVGIP* 48:357-393, doi:10.1016/0734-189X(89)90147-3. The (8, 4) pairing:
  the catchment is 8-connected (water steps diagonally), so the outside must
  be taken 4-connected for the outline to be a simple curve. The tracer's
  saddle rule below is that pairing. **Lorensen and Cline 1987**, "Marching
  cubes", *SIGGRAPH Computer Graphics* 21(4):163-169, doi:10.1145/37402.37422,
  is the family the tracer belongs to (its 2-D case, marching squares).
- **Kronenfeld, Stanislawski, Buttenfield and Brockmeyer 2020**,
  "Simplification of polylines by segment collapse: minimizing areal
  displacement while preserving area", *International Journal of Cartography*
  6(1):22-46, doi:10.1080/23729333.2019.1631535 (online 2019). APSC. From its
  abstract: segments are collapsed to Steiner points in priority order, with
  placement and displacement functions chosen so that area is preserved
  exactly, and "self-intersections can be avoided by testing for
  intersections with two new line segments associated with each segment
  collapse operation". **The reduction builds on it.** What differs: (a) the
  stopping rule is a horizontal tolerance against the fine outline, because
  that is what Ola asked for, where APSC targets a vertex count or scale;
  (b) the priority is that same deviation, not areal displacement; (c) the
  placement below is ours. The paper is paywalled and its placement function
  was not read. The area guarantee does not depend on it: any point on the
  area-preserving line keeps the area. (d) A containment check for the seed
  is added.
- **Buchin, Meulemans, van Renssen and Speckmann 2016**, "Area-preserving
  simplification and schematization of polygonal subdivisions", *ACM TSAS*
  2(1), doi:10.1145/2818373. The edge-move, which also keeps area and
  topology. The alternative if APSC's placement proves poor on lattice
  outlines.
- **Bose, Cabello, Cheong, Gudmundsson, van Kreveld and Speckmann 2006**,
  "Area-preserving approximations of polygonal paths", *J. Discrete
  Algorithms* 4(4):554-566, doi:10.1016/j.jda.2005.06.008. The optimisation
  version (fewest vertices). This increment is greedy and claims no minimum.
- **Visvalingam and Whyatt 1993** (*Cartographic J.* 30:46-51,
  doi:10.1179/000870493786962263) and **Saalfeld 1999** (*CaGIS* 26:7-18,
  doi:10.1559/152304099782424901): the vertex-removal and topology-consistent
  simplifiers. Neither keeps area; they are why off-node vertices matter.

**Novelty.** None is claimed. Seeding a watershed with a polygon of pour
points is standard in GIS tools; lake-as-seed with Priority-Flood labelling,
a marching-squares outline and APSC with a tolerance band is a combination of
published parts. If a later write-up wants to claim something (for example
the tolerance-band guarantee against the fine outline), the check is still to
be done.

## The design

### Data flow

```
cli.py  catchment
  |  parses flags into a CatchmentRequest (frozen Pydantic), builds the
  |  repository (paths stop here), calls delineate(), writes GeoJSON
  v
catchment.py  delineate(request, repository) -> Catchment    [no paths]
  1. seed point -> DEM CRS (crs.reprojector); lake polygon containing it
     -> DEM CRS             (lakes read by cli.py via feature_input.read_lakes)
  2. window = seed bounds grown by the margin
  3. loop: plan_mosaic + assemble (15a)  -> DemTile, read-only
           seed mask: DEM nodes inside the lake   (shapely.contains_xy)
           _core.upstream(raster_view, seed_mask) -> mask, bbox, flags  [C++]
           contained? done : grow the window and repeat
  4. outline.trace(mask) -> rings, lattice units   (numpy, no _core)
     pick the ring around the seed; drop the rest, count them
  5. _core.reduce_ring(ring, tolerance, keep=[seed])  [C++]
  6. -> Catchment(fine, reduced, counts, areas, windows)   frozen
```

The C++ core sees an elevation array, a byte mask, and a ring of doubles in
metres. No path, file or CRS crosses (CLAUDE.md §2, I/O boundary).

### Where it runs

- **The flood is C++.** It visits every node of the window: 12-13 M nodes for
  Bygdin, far more for Glomma. The Python prototype took 54 s on its
  13 M-node Bygdin window at 10 m; a priority queue in C++ should take a
  few seconds.
  `include/terrain/hydrology/upstream.hpp`, header-only, a template on the
  `RasterSource` concept, bound over the existing `RasterView` variant, GIL
  released.
- **The reduction is C++**, `include/terrain/vector_simplify/area_collapse.hpp`,
  because its crossing tests need the exact kernel (`noding/intersect.hpp`'s
  `classify<K>`, `DefaultKernel` in the binding). Both module names are the
  ones `project_structure.md` already plans.
- **The tracer is Python and numpy**, `src_python/tin_engine/outline.py`,
  never importing `_core`. It is raster topology on a byte mask, O(boundary)
  after one vectorised pass, and the prototype traced Bygdin in 0.09 s. Doing
  it in C++ would add a binding and a C++ suite tonight for no measured gain.
  Default (main session / @architect, 2026-09-29), for Ola to confirm;
  alternative: C++ next to the flood, when a Glomma-size run says so.
- **Planning and I/O are Python**, as 15a's `mosaic.py` and `dem_input.py`
  are. `catchment.py` takes a `DemRepository` (the Protocol in
  `io/repository.py`) rather than paths, so tests pass an in-memory one; the
  five lines in `open_dem` that build a repository from sources move to a
  helper in `dem_input.py` that both commands call.

### The seed

**Default (main session / @architect, 2026-09-29), for Ola to confirm: the
seed is a lake polygon, and the catchment is every DEM node that drains into
any node inside it.** The user gives a point (`--seed X Y`, in `--seed-crs`,
default EPSG:4326 so lon lat), and a polygon source (`--lakes`); the polygon
containing the point is the lake. Its DEM nodes are all seeds. For Bygdin:
`--seed 8.5425 61.3512` (mid-lake, UTM33 about (155000, 6819000), on the
1057.4 m surface, inside CORINE feature 54101).

Why this and not the alternatives:

- *A pour point at the dam, snapped to the strongest flow nearby*: needs flow
  accumulation (a second pass and a stored order), a snap radius that can
  jump to the wrong channel, and flat routing across the lake to reach the
  outlet. More code and more ways to be wrong. The brief's dam coordinate
  shows how easily the point itself is off.
- *The DEM flat under the point*: measured above, DTM10's Bygdin is not one
  flat; the flat would be 7 km² of 40.
- *The lake polygon*: no snapping, no flats, no dependence on where exactly
  the outlet is. Its cost: the polygon comes from outside the DEM. If it
  reaches past the outlet, whatever drains into that stretch of river joins
  the catchment. CORINE's Bygdin polygon ends at x 168460, NVE's catchment at
  x 168680 and the dam at x 168087, so a sliver of river below the dam may be
  included. The acceptance run measures it (nodes in ours and not in NVE's).

The source is any polygon layer: a GeoPackage table (`--lakes-layer`) or a
GeoJSON `FeatureCollection`. A multipolygon contributes the part containing
the point. No class filter: the polygon under the point is the one meant.
Refusals: the point in no polygon, or in two (overlapping input); both name
the point in the source's CRS.

**Who reads the lakes** (revised 2026-09-29, after the PR 1 red step).
`feature_input.py`, the module that already reads 16b's sources, gains
`read_lakes(path, layer, point, point_crs) -> (tuple of shapely geometries,
crs text)`. It reuses 16b's reading by splitting the file-reading half of
`_Tally.source` into a public `read_source(path, layer, attribute, box_for)`,
which `_Tally.source` then calls, so 16b's behaviour is unchanged:

- `.gpkg`: `io/repository.open_geopackage`, `io/geopackage.layer_info`
  (`--lakes-layer`, or the only features table), then
  `io/geopackage.query_features` with the box `box_for(layer CRS)`, here the
  seed point moved into the layer's CRS (a point box; the R-tree widening
  makes it a superset). The attribute column is the layer's primary key,
  since lakes need no class.
- `.geojson` / `.json`: `json` and `shapely.geometry.shape` per feature, the
  CRS from the `crs` member or WGS 84, as 16b does.
- `.gml`: comes free with `read_source`; not advertised.

`read_lakes` returns every polygon or multipolygon the query yields, in the
source's CRS; lines and points are skipped. Which one contains the point, and
the refusals above, are `catchment.py`'s: `CatchmentRequest` carries
`lakes` (a tuple of shapely geometries) and `lakes_crs`, and never a path.
`cli.py` calls `read_lakes` and builds the request. This matches the red
suite as committed: `test_catchment.py` passes geometries, and
`test_cli_catchment.py` passes files through the CLI.

**Without `--lakes`**, the seed is the one DEM node nearest the point, a pour
point with no snapping. It is there because it costs five lines and gives the
synthetic tests a seed without a polygon file. stderr says that a pour point
must lie on the flow line. Snapping is a later increment.

### Flow and membership: one flood

`upstream(z, seed)`:

- Every valid node on the window's edge, and every valid node with a NoData
  8-neighbour, is an outlet: pushed with its own z at the start. NoData nodes
  (the sentinel or NaN, as `RasterView::is_nodata` says) are never pushed and
  never in the catchment.
- Pop the lowest key (level, push counter). For each 8-neighbour not yet
  reached: its level is max(its z, the popped level); its label is *in* if it
  is a seed, else the popped node's label; push it.
- Outlets are labelled *in* if they are seeds, *out* otherwise.

The catchment is the nodes labelled *in*. A node is in it exactly when the
chain of nodes that flooded it, which is its drainage path over the filled
surface, passes through the lake. Levels are `double` (a float32 z is exact
in it). The counter makes equal levels first in, first out, so the result
depends on nothing but the input.

Memory: the elevation array is 15a's, borrowed. The flood's own is one byte
per node (unreached, out, in), which is also the returned mask, plus the queue
(24 bytes per entry; the frontier, not the window, in practice). The seed mask
is one more byte per node, owned by numpy.

Returned: the mask (a numpy `uint8` array the outcome owns), the number of
nodes in, the bounding rows and columns of the in-nodes, and two flags. A seed
mask whose shape is not the raster's is a `ValueError` in the binding.

**The flags** (revised 2026-09-29, after the PR 1 red step; the first wording,
"an in-node on the window's edge", could only ever fire for a seed, because
every edge node and every node beside NoData is an outlet and an outlet is in
only if it is a seed). The flags say where the catchment may continue beyond
what the window knows. An outlet's own drainage is unknown: the window assumed
it leaves, and in the full DEM it may instead run into the catchment. So:

- `touches_edge` is true when some in-node is an edge outlet or is an
  8-neighbour of one. On a grid that is exactly: an in-node in the first or
  last two rows or columns (row <= 1, row >= rows - 2, the same for columns).
- `touches_nodata` is true when some in-node is an outlet beside NoData or is
  an 8-neighbour of one, which puts it within two nodes of NoData.

An outlet can be both kinds; then both flags may be set. No comparison of
heights is made: an edge node lower than its in-neighbour might still, in the
full DEM, fill and spill back, so any contact counts. Conservative by design:
a flag never misses a truncated catchment, and a catchment that merely comes
within one node of an edge it does not cross is flagged too. The window loop
then grows past it, so the price is a growth step, not a wrong answer.

Checked against the red suite as committed (1e1b3bb): the C++
`check_invariants` asserts the flag when an in-node is on the edge or beside
NoData, and its absence when no in-node is within one node of the edge or two
of NoData; this rule sets the flag exactly on the first and never on the
second. The named cases (a seed on the edge, the V-valley, the flat lake far
from every edge, no seed, the degenerate rasters, the binding's seed beside
NoData) agree. **No test changes are needed.**

### The window

A read region grown until the catchment clearly lies inside it, the same
re-plan pattern as 15b's `_domain_plan`:

1. Start from the seed's bounds (the lake polygon's, or the point) grown by
   the margin, `WINDOW_MARGIN_M = 2000` metres. Default (main session /
   @architect, 2026-09-29), for Ola to confirm; no flag tonight.
2. Plan and assemble that box (15a), flood it.
3. and 4., revised 2026-09-29 after the PR 1 red step (the first wording
   could loop for ever near the data's edge, where "the bounds plus the
   margin fit" stays false while the window cannot grow). All comparisons
   are in node indices on the plan's lattice, and all are closed: a box that
   ends exactly on the window's last node line fits.

   Let E be the data's node rectangle on the chosen lattice (the union of its
   tiles' extents, which `plan_mosaic` clamps every window to), and W the
   window just flooded, already inside E.

   a. If `touches_nodata`, refuse: the catchment is truncated by missing
      data, which no growth can fix.
   b. The need N is the in-nodes' bounds grown by the **base** margin (2000
      m, never doubled), clamped to E. A side of W *can grow* when N reaches
      past W on that side; since N is clamped to E, that also means W is not
      yet at E there.
   c. If some side can grow: the step margin doubles (4000 m at the first
      growth, then 8000 m, ...), the next window is W joined with the
      in-nodes' bounds grown by the step margin, clamped to E, and the loop
      repeats from 2. Only the step doubles; the test in b always uses the
      base margin.

   Revised 2026-09-29 at the PR 1 green step (58f6904): the previous
   wording let the test in b use the doubled margin too, so the need grew
   as fast as the window, every step grew, and the loop only stopped at the
   data's edge. @developer's first Bygdin run took 8 minutes and was then
   refused by 15a's coverage check. With the base margin in b, Bygdin gives
   304.91 km² (-0.21 % against NVE) in 6.3 s over 3 windows.
   d. Otherwise decide by the flag. If `touches_edge`, refuse, naming the
      sides (from the bounds: an in-node in the first or last two rows or
      columns), each of which is then at E: the catchment is cut by the
      data's edge. If not, accept.

   **It terminates.** Step c runs only when N reaches past W on some side.
   The next window contains the in-nodes' bounds grown by the step margin,
   which is at least the base margin, so it contains N (both clamped to E)
   and has at least one more row or column than W; every window lies inside
   E, which is finite. So step c runs at most (rows of E + columns of E)
   times. **And it stops early**: after a growth step the window holds the
   last catchment's bounds plus the base margin, so a further step happens
   only if the catchment itself grew past that in the new window. On a
   catchment far from the data's edge the number of windows is one more
   than the number of times the catchment outgrew its window's base margin,
   in practice two or three. The memory cap (step 5) may refuse earlier. Every
   exit is an accept or a refusal from a or d.

   **Accepting means**: the catchment comes no closer than two nodes to any
   edge of the window, and either the window already holds the in-nodes'
   bounds plus the margin, or it is at the data's edge on the sides where it
   does not. The margin is a heuristic against the known limit below; the
   flag is the rule.

   **What a test should pin** (for @tester; described only): on a
   synthetic DEM whose data extends more than four base margins beyond the
   full catchment on every side, and whose catchment is larger than the
   first window (the bowl fixture, placed in a larger raster), the loop
   accepts in at most 3 windows, and the final window reaches the data's
   edge on no side. A second, cheaper pin: a catchment that lies inside the
   first window with the base margin to spare is accepted in exactly one
   window. The first of these fails on the pre-58f6904 loop, which grows
   until it meets the data's edge.

   Refusing a truncated catchment stays the default (main session /
   @architect, 2026-09-29), for Ola to confirm; alternative: write it with a
   warning under an `--allow-truncated` flag.
5. Memory: before each flood, refuse if nodes x (itemsize + 2) exceeds half
   the physical memory (`mosaic.physical_memory`, dtype-aware as 15a's cap
   is), naming the window's size. 15a's own cap on the array still applies
   inside `plan_mosaic`.

### A lake in a closed depression

Added 2026-09-29, after the PR 1 red step. **A lake seed need not have an
outlet, and nothing is refused.** Default (main session / @architect,
2026-09-29), for Ola to confirm.

Within a window every depression drains: Priority-Flood fills it to its
lowest rim and routes it out over that rim, so every lake has an outlet as
far as the flood is concerned, a real one or the spill point of its
depression. The catchment is still "every node whose flooding chain passes
through the lake". For the nodes of the depression itself that is what the
design wants when the lake polygon covers the depression, as it does for a
lake whose DEM surface is its water line: all those nodes are seeds.

Where it falls short: when the filled depression is larger than the lake
polygon (a lake drawn smaller than its DEM basin, or a dry pan beside it
below the spill level), the nodes of the depression outside the polygon are
flooded first in, first out from the spill point, so by distance in steps.
Those reached through lake nodes are in; those reached round the lake from
the spill point are out, although water there would run into the lake. The
error is bounded by the part of the depression outside the polygon and
reached before the lake; on Bygdin, where CORINE's polygon covers the DEM's
water surface, the acceptance's node comparison with NVE's polygon is where
it would show.

The alternative, for later: seed the whole depression. A first flood
computes the filled levels (one value per node, 4 bytes at float32, exact
because every level is some node's z); every node 8-connected to a seed
through nodes whose filled level is above their own z joins the seed; a
second flood labels as now. Twice the flood's time and 4 more bytes per
node, and a behaviour the red suite would need to pin (a seed inside a pit).
Not tonight.

**Known limit.** Not touching the window's edge does not prove the catchment
complete. A closed depression that straddles the window's edge drains out
through the edge in the window, but in the full DEM it may fill and spill into
the catchment. The margin makes this unlikely and each growth step makes it
less likely; nothing here rules it out. The acceptance run's comparison
against NVE is the check on Bygdin.

For Bygdin the first box is the lake's bounds plus 2 km: x 140.3-170.5 km,
y 6811.2-6825.6 km. NVE's catchment reaches y 6832.3 km, so one growth step is
expected, to a window of roughly 4100 x 3000 nodes (12 M, 49 MB at float32).

### The fine outline

`outline.trace(mask) -> list of rings`, each a closed ring of `(row, col)` in
lattice units, as halves (every vertex is the midpoint of a lattice edge
between an in-node and an out-node). This is marching squares on the 0/1
node values: in each square of four nodes, a boundary segment joins the
midpoints of its in/out edges, oriented with the in-nodes on the left. The
one ambiguous square, two in-nodes on one diagonal, is resolved by keeping the
two in-nodes connected (the (8, 4) pairing): the catchment is 8-connected
because water steps diagonally. The mask is padded by one row and column of
out-nodes so every ring closes.

**The guarantee.** Every ring is a simple closed polygon, and no two rings
share a point. Each boundary lattice edge has exactly one midpoint, used by
exactly one incoming and one outgoing segment; two segments in one square
never meet (the saddle case cuts two opposite corners); segments of
different squares meet only at shared midpoints. Every in-node lies strictly
inside an odd number of rings and every out-node inside an even number, at a
distance of at least a quarter of the cell's diagonal. Outer rings are
counter-clockwise in the world frame (x east, y north), holes clockwise.

**Holes and pinches.** A pinch (two parts meeting diagonally) is not a
special case: the saddle rule joins them, so the traced ring is simple. A
spur one node wide becomes a thin, simple ring around it. **Holes are filled**:
the fine outline is the one outer ring that contains the seed point (the
user's point with a lake, the pour node without), and every other ring is
dropped. stderr reports how many outer rings were dropped with their node
counts, and how many holes were filled with their area. Default (main session
/ @architect, 2026-09-29), for Ola to confirm; alternative: keep holes as
domain holes (16 supports them), which then need reducing too. If no ring
contains the seed point, `CatchmentError` (a data oddity: the point in a gap
the flood did not reach).

**Area.** The fine outline's area is close to the in-node count times the
cell area: a straight side sits half a cell outside the in-nodes, and each
convex corner cuts off an eighth of a cell. On Bygdin (prototype) 304.909 km²
against 304.9095 km². The fine area is the reference the reduction keeps.

Why not a polygon through the boundary nodes (the brief's "8-connected,
diagonal steps"): it lies half a cell inside the node area along the whole
perimeter (up to about 0.7 km² on Bygdin: half a cell times the fine outline's 140 km), and it
touches itself at every pinch, which needs a repair pass with its own proofs.
The midpoint ring needs none. Default (main session / @architect,
2026-09-29), for Ola to confirm.

The tracer converts to world coordinates only at the end: `x = x_min + col *
delta_x`, `y = y_max - row * delta_y`; halves of a 10 m lattice on 5 m
multiples are exact in doubles.

### The reduction

`reduce_ring(ring, tolerance, keep) -> ReduceOutcome` in
`vector_simplify/area_collapse.hpp`. Input: a simple counter-clockwise ring of
`Point2` in metres (Python subtracts the window's lower-left corner first, so
coordinates stay below about 10^5 and areas lose nothing), a tolerance in
metres, and points that must stay inside (the seed point). Output: the reduced
ring, a status (`Ok`, `InvalidTolerance` for negative or non-finite,
`NotCounterClockwise`, `TooFewVertices` below 4), and counts (collinear
vertices dropped, collapses made, candidates rejected for crossing, for the
seed, for tolerance).

1. **Collinear pass.** Drop every vertex exactly collinear with its two
   neighbours and between them (exact `orient2d` is zero). The lattice outline
   is full of them. This changes neither area nor shape.
2. **Candidates.** For each edge B-C with neighbours A before and D after,
   the collapse replaces B and C by one new point E. E lies on the line
   parallel to A-D at the signed distance that keeps the area of A-B-C-D
   equal to that of A-E-D (APSC's area rule). On that line, E is where it
   meets line A-B or line C-D, whichever gives the smaller deviation (ties to
   A-B); if both are parallel to it, the foot of B-C's midpoint.
3. **Deviation.** Each current edge carries the range of fine-ring vertices it
   stands for (from the original fine ring, before the collinear pass). For a
   candidate, the range is from A-B's start to C-D's end. Its deviation is the
   larger of: the greatest distance from a fine vertex in that range to the
   chain A-E-D, and E's distance to the fine segments in the range. A
   candidate is admissible when its deviation is at most the tolerance. The
   reference is always the fine ring, so errors never accumulate.
4. **Order.** A heap on (deviation, id of B), where ids are the fine
   indices and then a counter for new points. Least deviation first, so the
   outline moves as little as possible for each vertex it loses.
5. **Checks before a collapse**, all with the exact kernel on the doubles as
   given:
   - *Simple*: A-E and E-D do not meet any current edge except at A (with the
     edge ending at A) and at D (with the edge starting at D), and do not
     overlap those two. A uniform grid of current edges (bucket side the
     tolerance, at least one cell) finds the candidates; it is updated on
     every collapse.
   - *Seed inside*: each keep-point has winding number zero around the closed
     loop A-B-C-D-E-A and lies on neither new edge. (The difference between
     its winding in the old ring and in the new one is exactly that loop's.)
   - A rejected candidate waits until one of its four vertices changes.
6. **Apply** the best admissible candidate, re-evaluate the candidates whose
   four vertices include A, E or D, and repeat until none is admissible or 4
   vertices remain. Stale heap entries are skipped by a per-vertex version.

With tolerance 0 only the collinear pass runs.

**Guarantees**, each with the check the tests make:

- *Area*: equal to the fine ring's up to rounding. Tested as
  `|A_reduced - A_fine| <= 1e-9 * A_fine`; the acceptance run reports the
  difference in m².
- *Simple*: no two non-adjacent edges meet, adjacent edges share only their
  vertex. Tested by shapely `is_valid` and, on small cases, a brute-force
  pairwise test with exact orientation.
- *Seed inside*: the seed point is strictly inside (shapely `contains`).
- *Tolerance*: every fine vertex is within the tolerance of the reduced ring,
  and every reduced vertex within the tolerance of the fine ring. So the
  symmetric Hausdorff distance between the two is at most the tolerance plus
  one cell: tested with shapely `hausdorff_distance(..., densify=0.05)`.
- *Deterministic*: serial, ordered by (deviation, id); the same input gives
  the same bits.
- Not guaranteed: the fewest vertices (it is greedy), or that every
  catchment node is inside (a node within the tolerance of the boundary may
  fall either side).

**The tolerance**: `--outline-tolerance METRES`, default twice the DEM's cell
(20 m on DTM10). Default (main session / @architect, 2026-09-29), for Ola to
confirm. Divides are not known better than a cell or two (ours and NVE's disagree by
about 1.5 % of the area in nodes). The rough gauge above gives about 900
vertices for Bygdin at 20 m, from 17,800.

**Locality and parallelism** (the geometry skill asks): every collapse is
local (four vertices, a grid query); the heap is global and the loop serial.
At about 18,000 vertices that is milliseconds. It parallelises later by
cutting the ring into pieces fixed at their ends, if Glomma needs it.

### The command

```
rasputin catchment --dem PATH [--dem PATH ...] --seed X Y [--seed-crs CRS]
                   [--lakes PATH [--lakes-layer NAME]]
                   [--outline-tolerance METRES] --out FILE.geojson
                   [--out-parent DIR]
```

- `--dem` as `mesh` has it: one directory or several files (15a).
- `--seed X Y`, two numbers, like `--bbox`'s four. `--seed-crs` is anything
  pyproj reads, default `EPSG:4326` (so `LON LAT`), moved into the DEM's CRS
  by `crs.reprojector` (15b's one `from_crs` site, `always_xy`). All tiles in
  one CRS, or refuse, as 15b does.
- `--lakes`: `.gpkg` (with `--lakes-layer` when it has more than one feature
  table, as the CORINE extract does), `.geojson` or `.json`. Its CRS is the
  file's own (GeoPackage `srs_id`, GeoJSON `crs` member or WGS 84).
- `--out`: `.geojson` or `.json`, resolved and checked by the same
  `_destination` as `mesh`. Written: a `FeatureCollection` of one `Feature`,
  the reduced polygon in the DEM's CRS, with a `crs` member naming it
  (`EPSG:25833`), which `domain.read_domain` already reads. Properties: seed
  point and CRS, tolerance, node count, fine and reduced vertex counts, fine
  and reduced areas in m², the windows' sizes. Coordinates written with
  `repr` precision, so the area survives the round trip.
- `--outline-tolerance 0` writes the fine outline (collinear vertices only
  removed), for comparison.
- stderr, one line each: every window (box, nodes, flood seconds, and
  "grown: touches north" or "contained"); the seed (lake polygon with its
  area and seed-node count, or the pour node); the catchment (nodes, node
  area); the fine outline (vertices, area, rings dropped, holes filled); the
  reduced outline (vertices, area, the difference in m² and relative, the
  tolerance, seconds). Areas in km² with enough digits to see the
  difference.
- Every refusal is a non-zero exit and no file, as `mesh`'s are.

`catchment.delineate` is blocking (the C++ calls release the GIL); an async
caller runs it in `asyncio.to_thread`, as 16b's `open_features` is used. The
request and the result are frozen; the CLI is the only place with paths.

### New and changed files

| File | What | Production lines (estimate) |
|---|---|---|
| `include/terrain/hydrology/upstream.hpp` | the flood, `UpstreamOutcome` | 110 |
| `include/terrain/vector_simplify/area_collapse.hpp` | the reduction, `ReduceOutcome`, edge grid | 260 |
| `bindings/core.cpp` | `upstream`, `reduce_ring`, two outcome classes | 80 |
| `src_python/tin_engine/_core.pyi` | their stubs | 30 |
| `src_python/tin_engine/outline.py` | the tracer | 70 |
| `src_python/tin_engine/catchment.py` | request, seed, window loop, result | 170 |
| `src_python/tin_engine/dem_input.py` | repository helper split out | 10 |
| `src_python/tin_engine/feature_input.py` | `read_source` split out of `_Tally.source`, `read_lakes` | 30 |
| `src_python/tin_engine/cli.py` | `catchment` command, report, writer | 100 |
| `project_structure.md` | the two C++ modules and two Python modules | docs |

About 860 lines, over the 700 ceiling (CLAUDE.md §2), so two PRs on this
branch.

### The PR split

- **PR 1, the fine catchment** (about 510 lines): `upstream.hpp` and its
  binding, `outline.py`, `catchment.py`, the repository helper, and the
  `catchment` command writing the fine outline. It already answers "what is
  Bygdin's catchment" and can be compared against NVE. Red, green, review.
- **PR 2, the reduction** (about 350 lines): `area_collapse.hpp` and its
  binding, and the command reducing by default with
  `--outline-tolerance`. Red, green, review, then the Bygdin acceptance run.

Both go on `increment22-autocatchment`, PR 2 stacked on PR 1. Neither touches
refine or mesh code, so the 1 m benchmark and scaling sweep (README, rule 2)
do not apply; the Bygdin run below is this increment's acceptance.

## The red suites

Lean: no throwaway implementations, no mutation round. Default (main session
/ @architect, 2026-09-29), for Ola to confirm; the brute-force oracle and the
invariants on random inputs are the defence. None is named invariant-critical
for mutation testing tonight.

**PR 1, C++ (Catch2), `tests/cpp/unit/hydrology_upstream.cpp`:**

- *Brute-force oracle*. For random DEMs with no interior pit and no ties,
  built as z = 10 x (steps to the nearest edge) + a distinct fraction per
  node, every node drains to its lowest neighbour, and the flood's catchment
  of a random seed set must equal the set of nodes whose descent reaches it,
  node for node. (On such a DEM the flood pops in global z order, so each
  node is flooded by its lowest neighbour; the oracle is exact, not
  approximate.) A few hundred small grids, seeds of one node and of blobs.
- *Depressions*: a closed bowl whose spill leads into the seed's valley is
  wholly in; one that spills elsewhere is wholly out.
- *A flat lake on a plateau*: the lake nodes are the seed; a flat shelf
  beside it that spills into it is in; a shelf that spills away is out.
- *V-valley*: two planes meeting at a valley whose outlet is the seed; the
  catchment is the valley's side slopes up to the ridge lines, exactly.
- *Edge and NoData*: in-nodes on the window's edge set `touches_edge`;
  an in-node beside NoData sets `touches_nodata`; NoData nodes are never in;
  a seed on the edge is allowed; the bounds are the in-nodes' bounds.
- *Invariants on random DEMs with pits and flats*: seeds are in; every
  in-node has an 8-path to a seed through in-nodes; running twice gives the
  same mask.

**PR 1, Python:**

- `test_outline.py` (no `_core`): a single in-node gives one diamond of 4
  vertices and area half a cell; a 2 x 2 block; an L; a diagonal pair is one
  ring (the saddle rule); a ring of nodes gives an outer ring and a clockwise
  hole; random masks: every ring valid and simple, rings pairwise disjoint,
  every in-node inside an odd number of rings and every out-node an even
  number; mask on the window's edge closes (padding).
- `test_catchment.py` with an in-memory repository (as
  `mosaic_fixtures.py` builds): a lake polygon from GeoJSON seeds the nodes
  inside it; the point in no polygon, or in two, is refused; `--seed-crs`
  4326 lands on the same node as the DEM-CRS point; a catchment larger than
  the first window grows and ends equal, node for node, to one flood over the
  whole raster; reaching the data's edge, or NoData, is refused naming the
  side; the memory cap refuses (monkeypatched `physical_memory`); holes are
  filled and extra rings dropped, and the counts say so.
- `test_cli_catchment.py`: a synthetic tiled DEM directory and a GeoJSON lake;
  the output reads back through `read_domain` in the DEM's CRS; `rasputin
  mesh --dem ... --domain out.geojson --tolerance ...` succeeds on it; stderr
  carries the window, fine vertex count and fine area lines; a wrong `--out`
  suffix and `--lakes-layer` without `--lakes` are refused.

**PR 2, C++ and Python:**

- `tests/cpp/unit/area_collapse.cpp`: status for a negative, NaN or infinite
  tolerance, a clockwise ring, fewer than 4 vertices; tolerance 0 removes
  only collinear vertices; a traced rectangle of nodes reduces to at most 8
  vertices at a tolerance of one cell, with its area unchanged; a keep-point
  half a cell inside a notch that a collapse would cut off stays inside; a
  thin corridor (two long sides 1.5 cells apart) at a tolerance of 5 cells
  must not cross itself; the same input twice gives equal bits.
- `test_core_reduce.py`, on rings traced from random blobs and from the
  PR 1 synthetic catchments: the five guarantees above (area, simple, seed
  inside, tolerance at vertices, Hausdorff within tolerance plus a cell), and
  determinism.
- `test_cli_catchment.py` gains: the default reduces; the reduced line
  reports vertices, area and the difference; `--outline-tolerance 0` gives
  the fine ring.

## Acceptance: Bygdin, end to end

Run by `@perf` after PR 2 is green, recorded under
`docs/benchmarks/2026-09-29/bygdin/README.md` with the commands, the commit,
`pmset -g batt`, and the numbers:

1. `rasputin catchment` as at the top. Record every window (box, nodes,
   seconds), node count, fine and reduced vertex counts, fine and reduced
   areas and their difference, the times of flood, trace and reduction.
2. **Sanity on the area**: against NVE's 305.54 km² (delfelt 1187, the URL
   above). Expected within 2 %; the prototype gave -0.21 %. Outside 2 % is a
   finding to explain before merge, not a pass. Also the node overlap with
   NVE's polygon (fetched by the URL, not committed), both ways.
3. `rasputin mesh --dem ... --domain bygdin.geojson --tolerance 1` and
   `--tolerance 10`, each with `--stats`: triangles, vertices, time. And once
   at `--tolerance 10` with the fine outline (`--outline-tolerance 0`), to
   show what the reduction saves, which is Ola's first requirement.
4. A picture is optional; the `.vtk` opens in ParaView.

## What each persona reads

`@tester` and `@developer`: this file, then `docs/increments/15-dem-mosaic.md`
(the plan and the re-plan loop) and `docs/increments/16-domain-polygon.md`
(what `--domain` accepts). `@developer` also reads `noding/intersect.hpp` for
`classify` and `raster/view.hpp`.
