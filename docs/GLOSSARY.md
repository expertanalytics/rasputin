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
