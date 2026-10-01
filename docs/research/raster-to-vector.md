# Research note: a categorical raster as watertight constraint polygons

Status: research note, prior art only, no design. Written by `@architect`,
2026-10-01, at Ola's request. It feeds the increment that turns the São
Francisco basin's land cover (MapBiomas, a raster) into polygons, which comes
after the basin-scale design. How each citation was checked is stated beside
it; the convention is at the head of the Sources section.

## The question

Ola, 2026-10-01: the land-cover step should be "very similar to what we did
for the coarsening of the auto water catchment", but "the auto-catchment
utilizes that the raster is really a point raster and not a partition of the
domain. That does not really hold now."

What changes against increment 22 (`docs/increments/22-auto-catchment.md`):

| | increment 22 (catchment) | land cover (this note) |
|---|---|---|
| what a raster value is | a sample at a lattice node | the class of a whole cell (area-registered) |
| fine boundary | marching squares between nodes, vertices at edge midpoints | the cell edges themselves ("cracks"), a staircase with vertices at cell corners |
| how many regions | one, holes filled | thousands of faces of many classes, a partition of the domain |
| boundary pieces | one closed ring | arcs shared by two faces, meeting at junctions where three or more classes meet |
| simplification | one ring, each edge seen once | each shared arc simplified once, so both neighbours stay watertight; junctions fixed |
| small regions | dropped outer rings, filled holes | patches below a minimum size absorbed into a neighbour |
| CRS | the DEM's projected grid | geographic (EPSG:4326 per Ola), so vertices are reprojected afterwards |

Why it matters: every boundary edge is a constraint the mesh must honour. 16b
measured the cost: on a 48 km DTM10 square at a 10 m vertical tolerance,
CORINE's borders with the 25° quality start made the mesh 3.6x larger
(`docs/increments/16b-terrain-polygons.md`, the size table). CORINE is already
generalised (25 ha minimum mapping unit, 100 m minimum width; its vertices
were a median 54 m apart in 16b's probe). A 30 m MapBiomas staircase, with a
minimum mapping unit of half a hectare, has a vertex at every cell step, so
the reduction has to do much more work than CORINE ever needed.

## Legacy

The legacy tree holds a raster land cover, GlobCover, read but never
vectorised. The grep, against the tag:

```
$ git grep -il -E "polygoni|vectori|sieve|raster.?to.?vec|land.?cover|corine|marching|chain.?code" legacy-archive -- legacy
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/rasputin/land_cover_repository.py
legacy-archive:legacy/rasputin/tin_repository.py
legacy-archive:legacy/rasputin/web_visualize.py
legacy-archive:legacy/tests/test_gml_repository.py
legacy-archive:legacy/tests/test_land_cover_repository.py
```

What the hits hold (read from the tag):

- `globcov_repository.py`: `GlobCovRepository.constraints()` returns `[]`;
  `land_cover()` looks up the GlobCover class of given points by converting
  them to raster indices. `application.py` calls it with the mesh's
  triangle centres (`mesh.cell_centers`) to colour triangles. So the legacy
  answer to a land-cover raster was **sample the class at each triangle
  centre, add no constraints**. That is a real alternative to this note's
  subject (see "What combines", option 0) and the only thing worth carrying
  over: the idea, not the code (Pillow-based pixel reads).
- `gml_repository.py`, `land_cover_repository.py`, the tests: CORINE as
  vector GML, already handled by 16b.
- `tin_repository.py`, `web_visualize.py`: storage and display; nothing on
  vectorisation.

No raster-to-polygon code, no sieve, no simplifier other than CGAL's
surface-mesh collapse (increment 22's grep). Nothing else is carried over.

## Prior art

### 1. Boundary extraction from a labelled raster

**Cracks, not samples.** When a cell is a region, its boundary is made of the
edges between cells (1-cells) and the corners (0-cells). This cell-complex
view is what makes boundaries of a labelled image well defined without the
4- versus 8-connectivity paradox that increment 22 had to settle with the
(8, 4) pairing.

- **Kovalevsky 1989**, "Finite topology as applied to image analysis".
  Treats the image as a cell complex (pixels, cracks, points); boundaries
  are sequences of cracks. *Gives:* the right model for an area-registered
  raster; a boundary between two classes is a set of cracks, shared by
  construction. *Lacks:* any simplification. (Crossref record checked; the
  content is from memory of the paper, not re-read.)
- **Brice and Fennema 1970**, "Scene analysis using regions". Region
  growing over a grid where boundaries are the segments between pixels, and
  regions merge across weak boundaries. Usually credited with the crack-edge
  representation; that attribution is from memory and was not checked in the
  paper. *Gives:* the oldest region-merging-on-cracks reference, relevant to
  small-patch absorption too.
- **Freeman 1961**, "On the encoding of arbitrary geometric
  configurations". Chain codes: a boundary as a sequence of unit steps.
  *Gives:* the compact encoding of a staircase arc (one direction per crack,
  or per run). *Lacks:* topology between regions.
- **Rosenfeld 1970**, "Connectivity in digital pictures", and **Kong and
  Rosenfeld 1989**, "Digital topology: introduction and survey". The
  connectivity pairings. *Gives:* why a diagonal checkerboard corner (A B /
  B A) is the one degenerate case: whether the two A cells are one region or
  two is a choice, not a fact of the data.
- **Suzuki and Abe 1985**, "Topological structural analysis of digitized
  binary images by border following" (the algorithm behind OpenCV's
  `findContours`), and **Chang, Chen and Lu 2004**, "A linear-time
  component-labeling algorithm using contour tracing technique". *Gives:*
  linear-time tracing with the outer/hole hierarchy. *Lacks:* both are
  binary and trace each region separately through pixel centres, so two
  neighbours get two different boundaries; for a partition that is the wrong
  starting point.
- **Braquelaire and Brun 1998**, "Image segmentation with topological maps
  and inter-pixel representation", and **Damiand, Bertrand and Fiorio
  2004**, "Topological model for two-dimensional image representation:
  definition and optimal extraction algorithm". The labelled image as a
  planar map (combinatorial map) on the inter-pixel grid: faces are regions,
  edges are the maximal crack chains between two regions, vertices the
  junctions. The second title claims an optimal extraction algorithm.
  *Gives:* exactly the structure we need (arcs between junctions, each arc
  once, face adjacency). *Lacks:* simplification; these are image-analysis
  data structures. (Crossref records checked; no abstract available; the
  description is from memory and should be read before it is relied on.)
- **Teng, Wang and Liu 2008**, "An efficient algorithm for raster-to-vector
  data conversion" (Two-Arm Chains Edge Tracing, TACET). From its abstract:
  one scan line at a time, all classes in one pass, and it "constructs
  complete area topological relationship by recording the shared edge
  between two polygons only once"; scalable to large images. *Gives:* a
  streaming, row-by-row extraction of the planar map from a classified
  remote-sensing image, which is our case. *Lacks:* simplification.
- **Xu, Chen and Yu 2016**, "Vectorization of classified remote sensing
  raster data to establish topological relations among polygons". Title and
  venue checked; no abstract available, content not read. Listed because
  the title is our question; to be read before the design.
- **Wu and Sullivan 2003**, "Multiple material marching cubes algorithm".
  The node-sample route for many labels (a marching-squares relative of
  increment 22's tracer). *Lacks, for us:* it treats values as point samples,
  which is exactly what Ola says no longer holds.

**GDAL `Polygonize`** (method only; GDAL is prohibited, `CLAUDE.md` §2). Read
from the source on 2026-10-01 (`alg/polygonize.cpp`,
`alg/polygonize_polygonizer.h/.cpp`, master): since a 2023 rewrite it
implements Teng, Wang and Liu's TACET (the header cites the paper's DOI). Two
passes: the first labels connected components (4-connected by default,
`8CONNECTED` optional) row by row with a union-find of polygon ids; the second
walks the grid corners row by row with "two arms" per corner, opening,
extending and closing arcs along cell edges. Output: one polygon feature per
component, with its class value. Coordinates come from integer (row, col)
through the geotransform, so the two copies of a shared boundary are
identical, but each polygon is written separately: **the shared arc is
written twice**, and no simplification or sieving happens.

**GRASS `r.to.vect type=area`** (manual read 2026-10-01). Traces the perimeter
of each area along cell edges; with `-s`, "at each change in direction (i.e.,
each corner), the two midpoints of the corner cell ... are taken, and the line
segment connecting them is used to outline this corner" (corner cutting).
The output is GRASS's topological vector model, where "this border exists once
and is shared between two areas" (GRASS vector overview). *Gives:* the planar
map with shared boundaries. *Lacks:* reduction; `-s` adds no simplification and
does not bound displacement.

### 2. Simplifying a coverage, not a polygon

The requirement: a shared arc is simplified once, junctions do not move, no
arc crosses another, no face disappears or swaps side.

- **Douglas and Peucker 1973** and **Ramer 1972** (the same algorithm,
  found twice), and **Visvalingam and Whyatt 1993** (repeated removal of the
  vertex with the smallest effective triangle area). Single-line algorithms.
  *Gives:* a Hausdorff-type bound (Douglas-Peucker) or an area-based
  priority (Visvalingam-Whyatt). *Lacks:* any topology guarantee: both can
  make a line cross itself or a neighbour, and both keep a subset of the
  input vertices (no Steiner points), so on a staircase they keep lattice
  corners.
- **Saalfeld 1999**, "Topologically consistent line simplification with the
  Douglas-Peucker algorithm". From the abstract: a simple test added to the
  stopping condition of Douglas-Peucker-like algorithms "can guarantee that
  the resulting simplified polyline is topologically consistent with itself
  and with all of its neighboring features"; a dynamically updated convex
  hull finds conflicts. *Gives:* the key observation for us, which is that
  a conflict is not only a crossing. A whole small feature (a neighbour's
  vertex, an island ring) can end up on the other side of a simplified line
  without any segment crossing it, so the test must ask whether anything
  lies in the region swept between the old and the new chain. *Lacks:* area
  preservation, Steiner points.
- **de Berg, van Kreveld and Schirra 1998**, "Topologically correct
  subdivision simplification using the bandwidth criterion". From the
  abstract: a nearly quadratic algorithm that simplifies a chain with a
  maximum error ε, keeps given extra points on the same side, and avoids
  self-intersection; "applied as the main subroutine for subdivision
  simplification" so that "the resulting subdivision is topologically
  correct". *Gives:* the published blueprint for our case: simplify arc by
  arc between fixed junctions, with the other arcs' vertices as the points
  that must stay on their side. *Lacks:* area preservation and Steiner
  points; quadratic per chain.
- **Estkowski and Mitchell 2001**, "Simplifying a polygonal subdivision while
  keeping it simple". From the abstract: minimising the vertex count of a
  topology-preserving subdivision simplification without Steiner points is
  hard even to approximate (no n^(1/5) approximation unless P = NP);
  heuristics work well in practice. *Gives:* the reason to accept a greedy
  method and claim no minimum, as increment 22 already does.
- **Guibas, Hershberger, Mitchell and Snoeyink 1993**, "Approximating
  polygons and subdivisions with minimum-link paths". The fattening view
  (each boundary within a tolerance band, fewest links). *Gives:* the
  theory behind a tolerance band. *Lacks:* a practical subdivision
  algorithm; it is also where the hardness of the subdivision case starts.
- **Mustafa, Krishnan, Varadhan and Venkatasubramanian 2006**, "Dynamic
  simplification and visualization of large maps". From the abstract: every
  point of a simplified chain is within its tolerance of the original, no two
  chains intersect, gigabyte-sized maps, out of core. *Gives:* evidence that
  topology-safe simplification scales to large maps. *Lacks:* area
  preservation; built on graphics hardware.
- **Buchin, Meulemans, van Renssen and Speckmann 2016**, "Area-preserving
  simplification and schematization of polygonal subdivisions". From the
  abstract: an edge-move operation that preserves area and topology,
  quadratic time, and every non-convex simple polygon can be reduced with
  it. *Gives:* **the published area-preserving method for subdivisions**,
  the closest single reference to what we need. *Lacks:* a tolerance band
  against the fine boundary (its target is a vertex count or schematisation).
  An open copy is listed by Semantic Scholar (TU/e repository); not read.
- **Kronenfeld, Stanislawski, Buttenfield and Brockmeyer 2020**, APSC (the
  method of increment 22's reduction). From the abstract: segments collapse
  to Steiner points in priority order, with functions that minimise areal
  displacement "under the constraint that the areas of adjoining polygons are
  preserved exactly", and a test of the two new segments against
  intersections. The paper's demonstration is on single lakes. Paywalled, not
  read in full.

  **Does APSC extend to shared arcs and junctions?** Not in anything found
  (the paper's 21 citing works listed by Semantic Scholar were scanned by
  title on 2026-10-01; none is a coverage extension). The geometry says it
  extends, with two changes:
  - *Area.* A collapse A-B-C-D to A-E-D on an arc shared by faces L and R
    keeps the signed area between the two chains at zero, so whatever L
    gains R loses and vice versa, and the zero sum means both keep their
    area exactly. The abstract's "adjoining polygons" wording points the
    same way. Junctions are arc endpoints; a collapse never moves A or D,
    so junctions stay fixed if only interior vertices (B, C) collapse.
  - *Topology.* The paper's check (the two new segments against
    intersections) is not enough in a coverage, for Saalfeld's reason:
    a small island ring, or another arc's vertex, can lie wholly inside
    the loop A-B-C-D-E-A without touching A-E or E-D, and after the
    collapse it is on the wrong side. Increment 22 handled the one case it
    had (the seed point) with a winding-number test; in a coverage the test
    must cover every vertex of every other arc near the loop, which is
    what de Berg et al.'s "extra points" are.

  These two points are argued here, not taken from a source.
- **TopoJSON** (Bostock, specification and `topojson-simplify`, read
  2026-10-01). The specification stores shared "arcs" once and polygons as
  lists of arc indices. `topojson-simplify` assigns Visvalingam weights per
  arc and keeps "the arc endpoints!" (source comment); the source has **no
  crossing test**, so "topology-preserving" there means "shared arcs stay
  shared", not "arcs do not cross". **Harrower and Bloch 2006**
  (MapShaper) is the same arc-based idea as a web service.
- **GEOS/JTS coverage simplification** (Martin Davis; JTS
  `CoverageSimplifier` documentation and PostGIS `ST_CoverageSimplify`,
  both read 2026-10-01). Visvalingam-Whyatt variant on the edges of a valid
  coverage; nodes ("inner vertices shared by three or more polygons") do not
  move; "if the input is a valid coverage, then so is the result"; edges are
  never removed; inner edges only, if asked; the tolerance is "roughly equal
  to the square root of triangular areas". **It is already in rasputin's
  venv**: shapely 2.1.2 on GEOS 3.13.1 exposes `shapely.coverage_simplify`
  and `shapely.coverage_is_valid` (`pyproject.toml` requires `shapely>=2.0`,
  so relying on it would raise the floor to 2.1). Measured below.

**GRASS `v.generalize`** (manual read 2026-10-01). Simplifies (Douglas-Peucker,
Lang, Reumann-Witkam, vertex reduction) or smooths (Boyle, McMaster sliding
average and distance weighting, Chaikin, Hermite, snakes) each line on its
own. Because GRASS stores a shared boundary once, the result stays shared.
"Boundaries are not translated if they would intersect with themselves or
other boundaries": a whole boundary is left unsimplified rather than
partially simplified, and is written to an error map. *Gives:* the simplest
safe rule (reject and keep the original). *Lacks:* area preservation; the
whole-boundary fallback leaves long staircases untouched exactly where the
geometry is tight. **`v.clean tool=prune`** is the vertex-removal version
that keeps topology and never moves "first and last segment of the
boundary".

### 3. The staircase: smoothing versus simplification

A straight border at an angle becomes a staircase of cell edges. There are
two ways to remove it, and they pull in opposite directions for a mesh:

- **Smoothing** moves or adds vertices to make the line look natural:
  Chaikin 1974 corner cutting, McMaster's averaging in `v.generalize`,
  `r.to.vect -s`. GRASS's manual notes smoothing can raise the vertex count
  by up to 4000 %. *For constraints this is the wrong direction*: each
  vertex is a forced mesh vertex.
- **Simplification with Steiner points** replaces a staircase by the line it
  approximates. A subset-of-vertices method (Douglas-Peucker,
  Visvalingam-Whyatt, GEOS coverage simplification) can only keep lattice
  corners, so a straight diagonal stays a zig-zag of corners at small
  tolerances, and at larger ones it is biased toward one side: on the probe
  below the per-class area changed by up to 50 cells. APSC places new points
  off the lattice, and the line it chooses keeps the area, so a staircase
  becomes its mean line. That is the de-stairing we want.

Digital-geometry results say when a staircase *is* a straight line:

- **Debled-Rennesson and Reveillès 1995**, "A linear algorithm for
  segmentation of digital curves": recognition of digital straight segments
  in linear time (from the abstract: "a definition of digital lines using a
  linear double diophantine inequality").
- **Selinger 2003**, "Potrace: a polygon-based tracing algorithm" (read in
  full, 2026-10-01). Traces paths along pixel edges, resolves ambiguous
  corners by a turn policy, drops paths enclosing fewer than `turdsize`
  pixels ("despeckling"), then fits an optimal polygon where a segment
  "approximates" a path if every path vertex is within 1/2 (max-norm) of
  it, and finally fits Bézier curves. *Gives:* a principled zero-tolerance
  de-stairing criterion (half a cell in max-norm). *Lacks:* it is binary
  (black on white): tracing each class separately breaks the shared arcs.
- **Kopf and Lischinski 2011**, "Depixelizing pixel art". From the abstract:
  resolves the connectedness of diagonal neighbours, reshapes cells, fits
  splines. *Relevant only for* its handling of the diagonal ambiguity; it
  targets sprites of a few hundred pixels and outputs curves. Not a method
  for us.

### 4. Small patches

Every patch is at least one closed constraint ring, so patch count drives mesh
size as much as boundary length.

- **MapBiomas already sieves.** Souza et al. 2020, section 2.3.5 (read
  2026-10-01): contiguous regions of "less than or equal to half a hectare
  (i.e., approximately 5 pixels)" are reclassified to "the predominant class
  value of the neighboring pixels"; that is its minimum mapping unit, and a
  similar filter removes "single pixels or streams of pixels" along class
  borders in transition maps. Collection 3.1 is what the paper describes;
  the filter of later collections was not checked.
- **CORINE's minimum mapping unit** is 25 ha with a 100 m minimum width
  (EEA data hub, CLC2018 V2020_20u1: "Spatial resolution: 25 ha/100 m"). So
  the borders 16b met were coarse because of the MMU as much as anything.
- **GDAL `SieveFilter`** (source read 2026-10-01): raster in, raster out.
  Components by 4- or 8-connectivity; each component's "biggest neighbour"
  by size; a small component is merged into its biggest neighbour, chaining
  through neighbours until one is above the threshold, with a cycle guard.
  Merge rule: **largest neighbour by area**.
- **GRASS `r.reclass.area method=rmarea`** (manual): areas under a threshold
  are "substituted with the value of the respective adjacent area with
  largest shared boundary", via `v.clean tool=rmarea` ("the longest
  boundary with adjacent area is removed"). Merge rule: **longest shared
  boundary**, which also removes the most constraint length per merge.
- **Monmonier 1983**, "Raster-mode area generalization for land use and land
  cover maps". From the abstract: raster mode is better for this; thin or
  small polygons are partitioned and others grown, with class weights and
  priorities. *Gives:* thin features (one or two cells wide, such as
  rivers and roads) survive an area sieve and need a width rule too, as
  CORINE's 100 m has.
- **Haunert and Wolff 2010**, "Area aggregation in map generalisation by
  mixed-integer programming". From the abstract: aggregation of small areas
  in a land-cover partition, minimising class change with a semantic
  distance, compact shapes, hard size constraints; NP-hard, solved by MIP
  with heuristics. *Gives:* class-aware merge rules (forest into savanna
  before forest into water). *Lacks:* the cost is far beyond what a basin
  needs.
- **Saura 2002**, "Effects of minimum mapping unit on land cover data
  spatial configuration and composition". From the abstract: a larger MMU
  misrepresents sparse, fragmented classes and lets dominant classes grow.
  *Gives:* the caution for hydrology: class *areas* move when patches are
  absorbed, independent of how carefully boundaries are simplified.

### 5. Land cover in TINs and hydrological meshes

- **Vivoni, Ivanov, Bras and Entekhabi 2004**, and **Ivanov, Vivoni, Bras and
  Entekhabi 2004** (tRIBS). From the first abstract: TINs that "integrate
  information on the surface topography, hydrographic features and land
  surface characteristics". The mesh is driven by hydrological similarity;
  land cover is not, from the abstract, used as hard constraint edges.
- **Kumar, Bhatt and Duffy 2009** (PIHM's domain decomposition). From the
  abstract: unstructured grids "with user-specified geometrical and physical
  constraints", generated from GIS feature objects. *Gives:* precedent for
  vector layers as mesh constraints in a hydrological model. How land cover
  enters (constraint or attribute) was not read.
- **Marsh, Spiteri, Pomeroy and Wheater 2018** (`mesher`, used by the
  Canadian Hydrological Model, Marsh, Pomeroy and Wheater 2020). Read in its
  documentation, 2026-10-01: a categorical raster is **not vectorised**. It
  is a refinement criterion: with `method: mode` and a `tolerance` of 0.6,
  "each triangle must have 60% of one soil type". Vector constraints are
  shapefiles, which "may be beneficial to simplify ... so-as to avoid the
  creation of many small triangles". *Gives:* the main alternative to
  constraint polygons: triangles that are mostly one class, with no
  constraint edges at all. *Lacks:* class boundaries are not edges, so a
  triangle label is a majority, not exact; and refinement toward a class
  boundary makes small triangles along it anyway.
- **Burkhart et al. 2021** (Shyft v4.8, Ola a co-author): cells carry
  land-type fractions. Listed only as the downstream consumer's view: a
  fraction per cell is a third option beside constraints and majority
  labels.

### 6. A probe: GEOS coverage simplification on a staircase

Run 2026-10-01 in the repository's venv (shapely 2.1.2, GEOS 3.13.1), throwaway
script, not committed. A 60 x 60 raster of four classes (smoothed noise,
quantiles), each class the union of its cells, 5,715 vertices, a valid
coverage (`coverage_is_valid`):

| tolerance (cells) | valid coverage | vertices | worst class area change (cells) | worst Hausdorff, class boundary (cells) |
|---:|---|---:|---:|---:|
| 0.5 | yes | 3,868 | 0.00 | 0.00 |
| 1 | yes | 2,743 | 50 | 1.00 |
| 2 | yes | 2,340 | 39.5 | 2.24 |
| 4 | yes | 2,260 | 119 | 3.77 |

- Watertight in every run: `coverage_is_valid` true, and with
  `simplify_boundary=False` the union keeps the full 3,600 cells with no
  holes. With the default `True`, the frame of the domain is simplified too
  and the union lost 6.5 cells.
- Not area-preserving (the third column), and the tolerance is not a
  distance bound (2.24 cells at tolerance 2).
- A one-cell island inside another class was not removed at tolerances 0.5,
  2 and 5: it shrank to a triangle. JTS's documentation says small rings
  are "removed where possible"; GEOS 3.13.1 did not do it in this case. So a
  sieve before vectorisation is needed whichever simplifier is used.
- The vertex count stops falling: rings keep at least four vertices and
  nodes never move, so the floor is set by the number of faces and
  junctions, which is the sieve's job, not the simplifier's.

## What combines into a method for our case

Pieces, in order, each from a source above; the choice between alternatives
is the design's, and Ola's where marked.

0. **Whether to vectorise at all.** The alternatives are the legacy route
   (label each triangle by the class at its centre, no constraints) and
   `mesher`'s (refine until each triangle is mostly one class). Ola's plan is
   polygons as constraints, and 16c's labels depend on it; but a hybrid
   (constraints for a few classes such as water, majority labels for the
   rest) would cut the constraint load most. Question 1.
1. **Sieve in raster space first**, on the class raster as given (cell
   counts, in its own grid). Relabel components below a minimum area into
   the neighbour with the longest shared boundary (GRASS `rmarea`), which
   removes the most constraint length per merge; GDAL's largest-neighbour
   rule is the simpler alternative. A width rule for one- and two-cell
   threads (Monmonier) is a second, optional pass. Doing this on the raster
   keeps the result a partition by construction. NumPy only.
2. **Extract the planar map on the crack grid**: horizontal and vertical
   cracks where neighbouring classes differ, junctions at cell corners where
   three or more classes meet (or the diagonal two-class pattern), arcs as
   maximal crack chains between junctions, each stored once with its left
   and right class (TACET; Damiand et al.; TopoJSON's arcs). Collinear crack
   runs collapse to one segment at extraction. Coordinates are integers in
   the raster's (row, col) lattice, so shared vertices are exact.
3. **Reproject the fine arcs, then simplify in the computation CRS.**
   Every fine vertex goes through pyproj once; a shared vertex transforms
   identically for both faces, so the map stays watertight, and a crack is
   one cell long, so its bend is negligible. Simplifying in metres makes
   the tolerance a distance and the preserved area a true area. The other
   order (simplify in degrees, reproject the few survivors) is cheaper, but
   the area kept is in square degrees and long edges bend after reprojection
   (16b measured 3-6 m for a 10 km edge between EPSG:4326 and UTM 33 at
   60-71° N; not measured at 7-21° S). Question 3.
4. **Simplify each arc once with APSC, junctions fixed**: increment 22's
   `area_collapse.hpp`, generalised from one closed ring to many open arcs
   with fixed endpoints, a shared edge grid over all arcs, and the
   swept-region test of Saalfeld and de Berg et al. (no other arc's vertex
   inside the loop A-B-C-D-E-A). The tolerance band against the fine arc is
   increment 22's. Area per face is kept exactly by the argument in §2.
   Closed arcs with no junction (an island) need one fixed anchor vertex.
5. **Hand the arcs to 16b as linework.** The noder already merges shared
   edges, but arcs stored once never need it. Clipping to the basin outline
   is 16b's existing step.

The off-the-shelf comparison for step 4 is `shapely.coverage_simplify` with
`simplify_boundary=False`: watertight and already installed, but neither
area-preserving nor distance-bounded (probe, §6). It is the baseline the APSC
generalisation has to beat, and the fallback if it does not.

## Novelty

Searched, 2026-10-01: Crossref bibliographic queries ("Kronenfeld segment
collapse topology polygon coverage simplification", "area preserving segment
collapse shared boundaries polygon coverage", "topology preserving
simplification polygon coverage shared edges", "raster to vector conversion
land cover polygons staircase", "vectorization of classified raster maps
topological consistency", "land cover polygons constraints triangulated
irregular network", "unstructured mesh land cover boundaries hydrological
model constraint"); OpenAlex searches ("area-preserving simplification
polygon coverage shared boundaries", "categorical raster polygonization
simplification", "segment collapse polygonal subdivision area preserving",
"topology-preserving simplification of raster-derived polygons", "land cover
raster vectorization mesh generation constraints", "land cover boundaries
constrained Delaunay triangulation hydrological model"); the works citing
APSC on Semantic Scholar. No general web search engine was available to this
session, so grey literature (theses, software blogs) is under-covered.

Found: every piece separately (planar-map extraction from labelled rasters,
topology-safe subdivision simplification, area-preserving subdivision
simplification by edge-moves, APSC on single polylines, sieving rules, raster
land cover as a mesh refinement criterion). Not found: APSC applied to a
coverage's shared arcs, or any method that takes a categorical raster to
watertight, area-preserving, tolerance-bounded constraint polygons for a TIN.

**No novelty is claimed.** The combination in the previous section is
published parts put together, and Buchin et al. 2016 already have area- and
topology-preserving subdivision simplification. If a later write-up wants to
claim the APSC-on-a-coverage extension or the tolerance-band guarantee per
face, the check to do first is a full read of Buchin et al. 2016, Kronenfeld
et al. 2020 and Xu, Chen and Yu 2016, and a proper web search for grey
literature (USGS work on APSC by Stanislawski and co-authors in particular).

## Open questions for Ola

Answered by Ola on 2026-10-01. Each question is kept as asked, with the ruling
under it. Quotations are Ola's words; the rest is the ruling as relayed.

1. **Constraints for every class, or only some?** Every class boundary as a
   constraint is the most expensive option (16b: 3.6x on CORINE, which is far
   coarser than MapBiomas). The alternatives are a hybrid (water and maybe
   urban as constraints, the rest as a per-triangle majority label or class
   fractions, as `mesher` and the legacy code did), or no land-cover
   constraints at all.

   *Answered: the hybrid.* Water bodies and rivers are constraints; every
   other class is carried as a fraction per triangle.
2. **What minimum patch size, and which merge rule?** MapBiomas's own minimum
   is half a hectare (about 5 cells). For a mesh, the natural scale is the
   horizontal tolerance: for example, absorb patches smaller than a few
   times tolerance². Longest shared boundary (fewest constraints left) or
   largest neighbour (simplest), and is a class-similarity rule wanted (no
   forest absorbed into water)?

   *Answered: a minimum area, no merging.* A water body is a constraint only
   above a minimum area tied to the tolerance (a few times tolerance²); a
   smaller one stays a fraction of the triangles it falls in, and nothing is
   merged into a neighbour. Rivers as constraints come from the BHO drainage
   lines, not from the raster; only wide channels and reservoirs come from
   the raster, as water bodies.
3. **Which area must be kept, and where is the tolerance measured?**
   Simplifying after reprojection keeps each class's area in m² and makes
   the tolerance metres; simplifying in the raster's degrees is cheaper but
   keeps area in square degrees. Is exact per-class area (as for the
   catchment) a requirement here, or is "approximately" enough, which would
   let `shapely.coverage_simplify` do the job?

   *Answered: exact area.* Each lake is traced on cell edges, reprojected to
   metres, and reduced by increment 22's area-preserving segment collapse,
   with a no-crossing check against the other lakes, the BHO river lines and
   the domain outline. Water bodies are disjoint, so the shared-border and
   junction case of sections 2 and "What combines" no longer arises for
   constraints; it remains only if a later step wants class polygons.
4. **The diagonal checkerboard corner** (A B / B A): two A regions touching
   at a point, or one region pinched at it? For constraints, separate faces
   (4-connected, GDAL's default) are the simpler choice; it decides patch
   counts for the sieve too.

   *Answered: 8-connected water, free pinch points.* Water cells touching at
   a corner are one water body; its outline is one ring that passes each
   pinch point twice. The pinch points are not fixed: Ola observed that the
   area-preserving collapse itself turns a diagonal strip into a proper
   simple polygon (width c/√2 for cell size c) when the tolerance is about a
   cell or coarser. So the no-crossing check must let a pinch open but never
   let the two sides cross.
5. **Which MapBiomas collection and year**, and is the class legend to be
   reduced first (MapBiomas has several levels; merging to level 1 or 2 before
   vectorising removes many boundaries for free)?

   *Answered in part.* Ola: "We need to keep high resolution on vegetation
   types and crop farming types in Brazil." The default is MapBiomas's full
   legend, crop types included, with fractions stored sparsely per
   triangle; a class map may coarsen it. **The collection and year are still
   open.**

A further ruling, on the fractions themselves. Ola: "We could even have a
cutoff on the fractions. 0.1% soybean does not carry so much information."
And: "95% corn, 5% soybean _could_ become 100% corn. It's basically for crop
specific transpiration." Ola proposed "a local out-of-balance ledger, trying
to compensate for missing covers" in neighbouring triangles, kept by area,
not by fraction, and ruled that "the general idea, Floyd-Steinberg
dithering, applied to class areas instead of pixel intensities, and related
publications should be used to resolve this." The prior art for that is the
next section.

## Sources

Convention: **Crossref** means the DOI's Crossref record was fetched on
2026-10-01 and its title, authors, venue and year match the citation;
**abstract** means the abstract was read (from OpenAlex, Semantic Scholar or
Crossref); **read** means the text cited was read at the URL given; **not
read** means only the record was checked, and anything said about content is
marked as from memory in the text.

| Citation | DOI or URL | Checked |
|---|---|---|
| Brice, Fennema 1970, *Artificial Intelligence* 1(3-4):205-226 | 10.1016/0004-3702(70)90008-1 | Crossref; not read |
| Braquelaire, Brun 1998, *J. Visual Comm. Image Repr.* 9(1):62-79 | 10.1006/jvci.1998.0374 | Crossref; not read |
| Buchin, Meulemans, van Renssen, Speckmann 2016, *ACM TSAS* 2(1):1-36 | 10.1145/2818373 | Crossref, abstract |
| Burkhart et al. 2021, *GMD* 14(2):821-842 | 10.5194/gmd-14-821-2021 | Crossref |
| Chaikin 1974, *CGIP* 3(4):346-349 | 10.1016/0146-664X(74)90028-8 | Crossref |
| Chang, Chen, Lu 2004, *CVIU* 93(2):206-220 | 10.1016/j.cviu.2003.09.002 | Crossref |
| Damiand, Bertrand, Fiorio 2004, *CVIU* 93(2):111-154 | 10.1016/j.cviu.2003.09.001 | Crossref; not read |
| de Berg, van Kreveld, Schirra 1998, *CaGIS* 25(4):243-257 | 10.1559/152304098782383007 | Crossref, abstract |
| Debled-Rennesson, Reveillès 1995, *IJPRAI* 9(4):635-662 | 10.1142/S0218001495000249 | Crossref, abstract |
| Douglas, Peucker 1973, *Cartographica* 10(2):112-122 | 10.3138/FM57-6770-U75U-7727 | Crossref |
| Estkowski, Mitchell 2001, *SoCG '01*:40-49 | 10.1145/378583.378612 | Crossref, abstract |
| Freeman 1961, *IRE Trans. Electronic Computers* EC-10(2):260-268 | 10.1109/TEC.1961.5219197 | Crossref |
| Guibas, Hershberger, Mitchell, Snoeyink 1993, *IJCGA* 3(4):383-415 | 10.1142/S0218195993000257 | Crossref |
| Harrower, Bloch 2006, *IEEE CG&A* 26(4):22-27 | 10.1109/MCG.2006.85 | Crossref |
| Haunert, Wolff 2010, *IJGIS* 24(12):1871-1897 | 10.1080/13658810903401008 | Crossref, abstract |
| Ivanov, Vivoni, Bras, Entekhabi 2004, *WRR* 40(11) | 10.1029/2004WR003218 | Crossref |
| Kong, Rosenfeld 1989, *CVGIP* 48(3):357-393 | 10.1016/0734-189X(89)90147-3 | Crossref |
| Kopf, Lischinski 2011, *ACM TOG* 30(4) | 10.1145/2010324.1964994 | Crossref, abstract |
| Kovalevsky 1989, *CVGIP* 46(2):141-161 | 10.1016/0734-189X(89)90165-5 | Crossref; not read |
| Kronenfeld, Stanislawski, Buttenfield, Brockmeyer 2020, *Int. J. Cartography* 6(1):22-46 | 10.1080/23729333.2019.1631535 | Crossref, abstract |
| Kumar, Bhatt, Duffy 2009, *IJGIS* 23(12):1569-1596 | 10.1080/13658810802344143 | Crossref, abstract |
| Marsh, Spiteri, Pomeroy, Wheater 2018, *Computers & Geosciences* 119:49-67 | 10.1016/j.cageo.2018.06.009 | Crossref, abstract |
| Marsh, Pomeroy, Wheater 2020, *GMD* 13(1):225-247 (CHM) | 10.5194/gmd-13-225-2020 | Crossref |
| Monmonier 1983, *Cartographica* 20:65-91 | 10.3138/X572-0327-4670-1573 | Crossref, abstract |
| Mustafa, Krishnan, Varadhan, Venkatasubramanian 2006, *IJGIS* 20:273-302 | 10.1080/13658810500390794 | Crossref, abstract |
| Ramer 1972, *CGIP* 1(3):244-256 | 10.1016/S0146-664X(72)80017-0 | Crossref |
| Rosenfeld 1970, *JACM* 17(1):146-160 | 10.1145/321556.321570 | Crossref |
| Saalfeld 1999, *CaGIS* 26(1):7-18 | 10.1559/152304099782424901 | Crossref, abstract |
| Saura 2002, *IJRS* 23(22):4853-4880 | 10.1080/01431160110114493 | Crossref, abstract |
| Selinger 2003, "Potrace: a polygon-based tracing algorithm" | https://potrace.sourceforge.net/potrace.pdf | read |
| Souza et al. 2020 (MapBiomas), *Remote Sensing* 12(17):2735 | 10.3390/rs12172735 | Crossref, abstract; section 2.3.5 read |
| Suzuki, Abe 1985, *CVGIP* 30(1):32-46 | 10.1016/0734-189X(85)90016-7 | Crossref |
| Teng, Wang, Liu 2008, *Annals of GIS* 14:54-62 | 10.1080/10824000809480639 | Crossref, abstract |
| Visvalingam, Whyatt 1993, *Cartographic J.* 30(1):46-51 | 10.1179/000870493786962263 | Crossref |
| Vivoni, Ivanov, Bras, Entekhabi 2004, *J. Hydrologic Eng.* 9(4):288-302 | 10.1061/(ASCE)1084-0699(2004)9:4(288) | Crossref, abstract |
| Wu, Sullivan 2003, *IJNME* 58(2):189-207 | 10.1002/nme.775 | Crossref |
| Xu, Chen, Yu 2016, *Earth Science Informatics* 10:99-113 | 10.1007/s12145-016-0273-3 | Crossref; not read |

Software and data documentation, read 2026-10-01:

- GDAL polygonize and sieve, source on master:
  https://github.com/OSGeo/gdal/blob/master/alg/polygonize.cpp ,
  `alg/polygonize_polygonizer.h`, `alg/polygonize_polygonizer.cpp`,
  `alg/gdalsievefilter.cpp`
- GRASS manuals (8.5): https://grass.osgeo.org/grass-stable/manuals/r.to.vect.html ,
  `v.generalize.html`, `r.reclass.area.html`, `v.clean.html`, `vectorintro.html`
- JTS `CoverageSimplifier`:
  https://locationtech.github.io/jts/javadoc/org/locationtech/jts/coverage/CoverageSimplifier.html ;
  PostGIS: https://postgis.net/docs/ST_CoverageSimplify.html ; shapely 2.1.2
  `coverage_simplify` docstring in the venv
- TopoJSON: https://github.com/topojson/topojson-specification ,
  https://github.com/topojson/topojson-simplify (`src/presimplify.js`)
- mesher: https://mesher-hydro.readthedocs.io/en/latest/configuration.html
- CORINE CLC2018: https://www.eea.europa.eu/en/datahub/datahubitem-view/a5144888-ee5a-4e5d-a7af-ccbf8ba8a8b2
