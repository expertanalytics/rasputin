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

## Refinement

**Sup-norm error (of a triangle).** The largest `|z_DEM(n) − P(n)|` over every
valid DEM node `n` in the closed triangle, its own vertices excluded, where `P`
is the plane through the triangle's three vertex heights. Refinement stops when
every triangle's sup-norm error is at most `--tolerance`. Ola's choice, recorded
in `docs/increments/14-adaptive-refinement.md`.

**Worst node.** The node attaining a triangle's sup-norm error; ties go to the
smallest `(row, col)`. Refinement's default insertion point.

**Void triangle.** A triangle with at least one vertex whose DEM height is
NoData. It has no plane, so it has no sup-norm error.

**Carving.** Refinement's treatment of a void triangle: insert the valid DEM
node nearest a NoData vertex (squared lattice distance, ties to the smallest
`(row, col)`), repeatedly, until no void triangle's closed node set holds a
valid node. `docs/increments/14-adaptive-refinement.md` R6.

**Trim.** The step after refinement that drops every triangle with a NoData
vertex, and every vertex then unused, and counts the dropped vertices.
`src_python/tin_engine/elevation.py`.

**Foot (of a node on a segment).** Let `N` be a DEM node and `S` a constraint
segment from `A` to `B`, both in world coordinates. With
`s = clamp(((N − A)·(B − A)) / |B − A|², 0, 1)`, the foot is `F = A + s·(B − A)`,
the point of the closed segment nearest `N`; the distance from `N` to `S` is
`|N − F|`.

In refinement (increment 20b, `include/terrain/refinement/refine.hpp`,
`detail::foot_of`), when the worst node `N` of a non-void triangle `T` has not
already been given a foot, `T`'s three edges are tried in edge order, skipping
any edge that is not constrained and any edge `N` lies exactly on (exact
orientation zero; that case is refinement's ordinary edge split at `N`). On the
**first** remaining edge with `|N − F| < ε`:

- if `F` is within `ε` of either end (`s·|S| < ε` or `(1 − s)·|S| < ε`), no foot
  is used: `N` is inserted, and nothing is counted;
- otherwise `F` is the foot. It is inserted instead of `N`, splitting that edge
  (and the neighbour's side) in two, with height `vertex_z` at `F`: bilinear
  interpolation of the DEM (increment 16, R0). `N` is recorded as footed.

If no edge qualifies, `N` is inserted as usual. A foot is **refused**, `N`
inserted instead and `feet_refused` incremented, in exactly two cases: `F`'s
DEM cell has a NoData corner (`vertex_z` gives no height), or a resulting
triangle on either side of the edge would not be strictly counter-clockwise
under the exact kernel (`detail::foot_fits`). A node is footed at most once. A
split deferred because a neighbour was already written this round (increment
14b's skip rule) is neither counted nor footed; it is retried next round.
`docs/increments/20b-min-insertion-distance.md`.

**ε (constraint-foot distance).** In increment 20b (`detail::foot_epsilon`),
`ε = clamp(tol / G, m / 100, m / 2)` with `m = min(dx, dy)` the smaller DEM
spacing and `tol` the `--tolerance`. `G` is the largest slope bound over the up
to four DEM cells that have `N` as a corner, leaving out any cell outside the
grid or with a NoData corner; a cell's bound is
`hypot(max(|z01 − z00|, |z11 − z10|) / dx, max(|z10 − z00|, |z11 − z01|) / dy)`.
When `G` is 0 (flat, or no valid cell), `ε` is the cap `m / 2`, including at
tolerance 0. Within `ε` of `N`, height changes by at most about `tol`, which is
the scale refinement resolves (Ola).

**Triangle budget.** `rasputin mesh --triangles N [--max-triangles M]`: the
mesh is asked for by its triangle count instead of by `--tolerance`. A
search over the tolerance finds one whose written mesh has between 2N − M
and M triangles (the *window*; M defaults to the larger of ⌊1.05 N⌋ and
N + 1), and writes exactly the mesh `--tolerance` writes at it; the
tolerance found is the record's `tolerance_m`. Never more than M. Each
tolerance tried is a *trial*, at most 12. `src_python/tin_engine/budget.py`;
`docs/increments/104-triangle-budget.md`, section 3.

**Error summary (`--rms`).** One pass over the mesh as written: the
root-mean-square of |DEM height − mesh height| over the DEM nodes (on a
reprojected run, the source's nodes) inside the mesh, each counted once,
those on a vertex skipped. Record keys `rms_error_m` and
`error_nodes_measured`; the summary line adds "root-mean-square X.XX m".
Always on a budget run; on a `--tolerance` run with `--rms`.
`include/terrain/refinement/measure.hpp`;
`docs/increments/104-triangle-budget.md`, section 4.3.

**Floor (of a triangle budget).** The mesh no tolerance can shrink: the
outline, the feature lines and the start-quality points, refined at a
tolerance of 10⁶ m. Its largest error is the tolerance at and above which
every mesh is the floor, so a budget answered by the floor reports that
error as its tolerance. A cap under the floor's count is refused.

**Saturated (trial).** A trial whose largest error is 0: every valid DEM
node is a mesh vertex, so no smaller tolerance changes the mesh, and the
search stops there.

## Land cover

**Legend.** The colour scheme of a code system: a TOML file with each
class's label and colour and the ParaView preset name, nothing more. Shipped
ones are in `src_python/tin_engine/data/legends/` (`corine`,
`mapbiomas-c11`, `ar50`); `rasputin palette` takes a legend's name or file.
A land-cover class file names the legend it colours by. Read by
`tin_engine/legends.py`; `docs/increments/108-legend-split.md`.

**Edges.** The class borders and lines the mesh must follow, so no triangle
crosses them: `--edges` chooses them from what `--land-cover`, `--prefer`
and `--lines` give (`all`, `water`, a category of the land-cover class file,
a class code, a `--lines` kind). Not to be confused with `--out-edges`,
which writes the mesh's edges to a PLY file.
`docs/increments/111-one-land-cover-map.md`.

**Land-cover class file.** A TOML file that says what one land-cover source's
classes mean: the attribute that holds its classes, each class's code, which
classes are water, the source's positional accuracy, its notice, the legend
its codes are coloured by, and, in a `[cells]` table, which cell category and
`--prefer` class each code counts as. Shipped ones are in
`src_python/tin_engine/data/land_cover_classes/`; one is chosen with
`--land-cover-classes`, by name or path, or by a sources-file `classes` key.
A class file may split classes by a second attribute (`[[split]]` tables)
into sub-classes with codes of their own; a value it does not list keeps the
parent's code. Read by `tin_engine/land_cover_classes.py`;
`docs/increments/108-legend-split.md`.
