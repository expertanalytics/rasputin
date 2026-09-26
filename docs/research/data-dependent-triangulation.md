# Research lead: data-dependent triangulation

Status: **to look into.** Ola asked for this on 2026-09-26, after the 1 m
comparison against earlier increments. Written by the main session. The
citations are from memory; each must be checked against the paper before this
document is used in a design or in any publication.

## What it is

Our mesher picks which nodes to insert with a data criterion: the worst
vertical error. Once a node is in, it connects the nodes by geometry alone,
using Delaunay flips (14b). Data-dependent triangulation (DDT) uses the data
for that second choice as well. It keeps the same vertices but flips edges to
lower the error of the piecewise-linear surface, instead of to maximise angles.

- Dyn, Levin and Rippa, "Data dependent triangulations for piecewise linear
  interpolation", IMA J. Numer. Anal. 10 (1990). Cost functions over edges,
  such as the angle between the normals of the two triangles or jumps in the
  gradient. The search is a local optimisation by edge flips.
- Rippa, "Long and thin triangles can be good for linear interpolation", SIAM
  J. Numer. Anal. (1992). For an anisotropic surface, the best triangles are
  long and thin, aligned with the direction of least curvature.
- Garland and Heckbert, "Fast polygonal approximation of terrains and height
  fields", CMU-CS-95-181 (1995). This is the greedy insertion method we
  effectively rebuilt. The report also tests data-dependent triangulation in
  place of Delaunay. Its reported effect on the number of triangles needed for
  a given error has to be read in the paper; don't quote it from memory.

## Why it matters for us

- **Mesh size.** Our output is ~428 k triangles at 1 m on the quarter circle
  (scratchpad `bench1m/REPORT.md`). If DDT reaches the same sup-norm tolerance
  with noticeably fewer triangles, that beats any speedup we can get from the
  scan or the rounds.
- **It conflicts with 20c.** 20c wants a minimum-angle criterion, at the start
  and during refinement. Rippa's result says good triangles for terrain can be
  thin. Ridges, valleys and coastlines are anisotropic. **Before 20c is
  designed, Ola has to decide what angles are for:**
  - *Numerical use downstream.* For example, a solver on the mesh wants shape
    regularity. Then the angle criterion is right, and DDT is at most a
    tie-breaker.
  - *Approximation quality per triangle.* Then an angle floor spends triangles
    on the wrong thing, and DDT (with a weak angle floor against true
    degeneracies) is the better target.
- **It also affects the publication option.** A sup-norm guarantee with
  constraints and data-dependent connectivity may be new. The data-independent
  version (ours today) is 1995 work.

## Questions to answer

1. **Size.** At equal sup-norm tolerance, how many fewer triangles does DDT
   need on our tile at 10 m and 1 m? Measure against today's Delaunay output,
   with the benchmark method (Release, fixed thread count, on AC power).
2. **Cost function.** Which cost fits a sup-norm tolerance? Dyn–Levin–Rippa's
   costs are smoothness measures, not max error. Candidate: flip when it lowers
   the max error over the DEM nodes in the quad. That is exact and local, and
   our row-span scan already computes it.
3. **When to flip.** Options:
   - after each round, over the touched slots;
   - once at the end;
   - both.
   Local optimisation by flips can cycle or stall in a poor local minimum. We
   need a termination argument like the one for Lawson (a potential that
   strictly drops).
4. **Constraints and feet.** Constrained edges never flip, as today. Does
   20b's ε rule still make sense when connectivity follows the data?
5. **Determinism and parallelism.** Can DDT flips run in the parallel rounds
   with the same serial-split discipline? Or does it need its own serial phase?
6. **Angles.** Does DDT produce true slivers (near-zero angles with no
   approximation benefit)? If it does, what is the weakest angle floor that
   removes them?

## A cheap first experiment

Take today's 1 m quarter-circle output. Run a serial post-pass that flips any
unconstrained convex quad whose flip lowers the max error over its DEM nodes,
until none does. Then run greedy insertion again from that connectivity to the
same tolerance, and compare triangle counts. This is a prototype in the
scratchpad, not production code. It tells us whether the idea pays before we
design anything.

## Related

- 20c (soft quality penalty): `docs/increments/20-start-quality.md`, "Ola's
  rulings".
- Also for the retrospective: simplification from dense meshes, meaning
  Garland and Heckbert 1997 (quadric error metrics) and Lindstrom and Turk.
