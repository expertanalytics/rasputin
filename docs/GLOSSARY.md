# Glossary

Terms this project uses in a specific sense. A word here means what this file
says it means, everywhere in the repository.

**Add a term when you first use it in a document, a comment or a commit
message.** A term used and not defined is a term the next reader has to
reverse-engineer.

Each entry says what the term means and points at where it is used canonically.
Definitions, not explanations: the reasoning belongs in the design that needed
it.

---

## Testing

**Oracle.** The independent source of truth a test checks against — what says
the right answer *is*, computed some way other than the code under test computes
it. `prop_noding_no_crossings.cpp`'s `split_sequence` is one: it finds where a
segment should be split by scanning every node directly, where the noder uses a
broad phase and an arc sort.

**Self-confirming oracle.** An oracle that asks the code under test for the
answer it is supposed to be checking, so it passes on any internally consistent
result including a wrong one. Named because it has happened here repeatedly:
`docs/increments/05-noder.md` records four. The rule it produced is that an
oracle is built from the **input**, never from the producer's own record of what
it did.

**Mutant.** A deliberate defect introduced into a working implementation to
check that a suite notices. A mutant that survives means the suite does not
cover what it claims to. A mutant with no killer input is struck rather than
kept, because it is a claim the round would otherwise ship as covered.

**Mutation round.** Running a suite against a set of named mutants. Required
only for an increment's **invariant-critical** suites, named as such in its
design record. `docs/increments/README.md` gives the cost argument.

**Invariant-critical.** Of a suite: the one whose failure would mean a stated
invariant is broken, rather than a value being wrong. These get a mutation
round; the others do not.

**Red step, green step.** The two halves of the TDD loop. The red commit adds a
failing suite and touches no production file; the green commit makes it pass and
touches no test file. `docs/increments/README.md` states the rule and why the
trace matters.

**Finding.** In the renderer: a defect the scene builder detected in the mesh or
the constraint set, drawn as a red dashed overlay. A finding on a picture is the
*correct* output for defective input, not a bug in the drawing.

## Geometry

**PSLG.** Planar straight-line graph: the validated constraint set handed to the
triangulator. Vertices plus chains. `include/terrain/core/pslg.hpp`.

**Chain.** A run of vertex indices with a role and a property set.
`include/terrain/core/pslg.hpp`.

**Role.** What a chain is for: `Outer` (the domain boundary), `Hole` (a region
excluded from the mesh), `Breakline` (a constraint the mesh must respect without
bounding anything). Only `Outer` and `Hole` are closed — `is_closed` in
`pslg.hpp` is the statement.

**Breakline.** A chain the triangulation must have edges along, which encloses
nothing. An **area feature** — a forest, a lake — is a *closed* breakline: its
last index repeats its first, and the terrain inside it is meshed.

**Constraint.** Any chain's edges, considered as something the triangulation may
not cut across.

**Predicate.** An exact geometric test returning a sign or a classification, not
a number: orientation, incircle, segment relation. `include/terrain/predicates/`.

**Kernel.** The bundle of predicates a template is parameterised over, so an
algorithm can be run against a filtered or an exact implementation without
changing. `include/terrain/predicates/kernel.hpp` is the concept.

**Degenerate.** Of an input: one where a predicate returns the boundary case —
collinear points, zero area, coincident vertices. Not an error by itself; each
increment's design states what its degeneracies do.

## Noding and snapping

**Noding.** Splitting constraints at their intersections so that no two
constraint edges cross in their interiors, and every crossing is an explicit
shared vertex. A crossing left unsplit is a *ghost point*: a place the geometry
implies a vertex and the data does not name one.

**Snap grid.** The lattice every noded output coordinate lies on, one scalar
spacing anchored at zero. What it buys is **not** accuracy: it makes coincidence
decidable, because two points snap to the same index if and only if their
snapped coordinates are bit-identical doubles. `include/terrain/core/snap_grid.hpp`.

**Snap rounding.** Constructing intersection points and rounding them to the
snap grid, so the result is representable and coincidence is decidable. The
standard alternative — exact constructions over rationals — left with CGAL.

**Hot pixel.** The closed square cell of side `spacing` centred on a lattice
point: `[(ix-½)s, (ix+½)s]` in each axis, per `snap_grid.hpp`. Borrowed from the
snap-rounding literature; a "pixel" here is a cell of a coordinate lattice and
has nothing to do with rasters or images. The project has three unrelated grids
— the DEM raster, the broad-phase spatial index, and this one — and only this
one affects correctness.

**Broad phase.** The spatial index that proposes candidate segment pairs so the
noder does not test every pair against every other. Its cell size is unrelated
to the snap grid's. False positives are free; false negatives are not.

## Rasters

**Tie point.** The GeoTIFF `ModelTiePointTag` (33922): one raster position
`(col, row)` paired with the world position it sits at. Together with the pixel
scale it is the whole georeferencing of a north-up file.
`docs/increments/11-raster-ingestion-prior-art.md` §3.2.

**Pixel scale.** The GeoTIFF `ModelPixelScaleTag` (33550): the world distance
between neighbouring raster nodes, stored positive in both axes. The fact that
y decreases as the row index grows is implied, not stored.

**Grid-registered (pixel-is-point).** The convention that a raster sample is a
point at a node, so `n` samples span `n - 1` spacings. `GTRasterTypeGeoKey`
(1025) value **2**, `RasterPixelIsPoint`. `include/terrain/raster/geometry.hpp` assumes it.

**Area-registered (pixel-is-area).** The other convention: a sample is a cell,
the tie point names a cell corner rather than a node, and `n` samples span `n`
spacings. `GTRasterTypeGeoKey` value **1**, `RasterPixelIsArea` — the
specification's own default. Reading such a file as grid-registered
shifts everything by half a cell. The reader converts rather than refuses, and
does not lean on that default: an absent key is refused, because half a cell is
the entire quantity in dispute. `docs/increments/11-raster-ingestion.md`
ruling 4.

**Mosaic.** A region assembled from several raster files, each with its own
extent and possibly its own CRS. The legacy walked one in
`legacy/rasputin/reader.py:434-456`.

**Half-cell shift.** The conversion from an area-registered file's declared
corner grid to the node grid `RasterGeometry` needs: `x_min = tie_x + dx/2`,
`y_max = tie_y - dy/2`, with the shape and spacings unchanged. Exact, not a
heuristic.

**NoData sentinel.** The single value that marks an absent sample in a raster.
Discovered in Python from `GDAL_NODATA` (42113); compared in C++ with `==`, so
it must be passed exactly as decoded and in the array's own dtype.
`include/terrain/raster/raster.hpp`'s `is_nodata` also treats NaN as NoData
unconditionally, with no sentinel needed.

**`always_xy`.** The `pyproj.Transformer` argument that forces longitude/x
first, ignoring the authority's declared axis order. Mandatory on every
transformer in `src_python/`. Omitting it on `EPSG:4326` puts a point about
6 000 km away with no exception. `docs/increments/11-raster-ingestion.md`
ruling 8.
