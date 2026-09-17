---
name: computational-geometry
description: Constrained Delaunay triangulation, terrain meshing and refinement design for the C++ core. Use when working on CDT construction, mesh topology, refinement policies, geometric predicates, breaklines and holes, elevation sampling interfaces, or adversarial geometry (collinearity, cocircularity, extreme scales).
---

# Agent Skill: Computational Geometry & Terrain Meshing

Core domain constraints and architecture for the C++ terrain-meshing engine (CDT and local refinement).

## 1. Domain Pipeline & Data Separation
* **Strict Abstraction:** Isolate geometric topology from terrain sampling, refinement policies, and I/O. Triangulation kernels must not hold DTM or CRS data directly.
* **Pipeline Sequence:**
  1. Accept bounding outer polygon + optional inner holes/breaklines via abstract geometric concepts.
  2. Construct a minimal 2D projected Constrained Delaunay Triangulation (CDT). *Do not pre-densify*.
  3. Associate elevation (Z) lazily or via dynamic sampling interfaces at the mesh boundaries.
  4. Evaluate error against reference DTM and refine locally until error tolerance is satisfied.

## 2. Algorithmic Architecture
* **Locality & Decoupling:** Every geometric mutation must have local consequences. Avoid global priority queues, global retriangulation, and whole-mesh synchronization.
* **Spatial Partitioning:** Design for independent subdomain processing. Interface reconciliation must happen over narrowly defined boundaries without global propagation.
* **Parallel Suitability:** Favor data-parallel layouts (Structure of Arrays where applicable), regular memory access, and branch-minimized paths. Avoid pointer-heavy graph structures and recursive walks.

## 3. Robustness & Numerics
* **Dependencies:** Do not use CGAL.
* **Filtered Arithmetic:** 
  1. Evaluate orientation, incircle, and intersections using fast hardware floats.
  2. Determine numerical ambiguity.
  3. Fall back to exact/arbitrary-precision filtering *only* when the float sign/ordering is unsafe.
* **Degeneracy Handling:** Actively handle non-general positions (duplicate/coincident vertices, long collinear sequences, cocircular points, and acute angles) through explicit test cases.
* **Test in the producer's relation, not the prose one:** an oracle must be
  written in the relation the code under test actually used. The operative test:
  *does this check use the relation the producer used, or the one that reads more
  naturally in prose?* In a snapping module those differ on ~93 % of the input —
  a snapped node is *near* a segment, not *on* it — so an exact-incidence oracle
  over snapped output is red on correct output, and an oracle re-derived from the
  producer's own bookkeeping is a self-confirming invariant that cannot fail when
  the thing it checks is broken. Borrow the producer's *predicate*, never its
  records. Worked example, four occurrences and each one introduced by the fix to
  the previous: `docs/increments/05-noder.md`, guarantees 14 and 15.

## 4. Change & Review Boundaries
* **Change Limit:** see the ceiling in `CLAUDE.md` §2.
* Every major algorithmic change must explicitly document its impact on locality, parallelization, and numerical robustness.

