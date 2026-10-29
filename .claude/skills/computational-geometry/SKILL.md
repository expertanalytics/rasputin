---
name: computational-geometry
description: C++ meshing core: CDT construction, refinement, geometric predicates, breaklines and holes, elevation sampling and adversarial geometry.
---

# Computational Geometry & Terrain Meshing

## 1. Pipeline
Topology stays apart from sampling, refinement policy and I/O; a kernel holds
no DTM or CRS data. Outline, holes and breaklines enter through concepts; a
minimal CDT, not pre-densified; elevation through an interface; local
refinement until the tolerance holds.

## 2. Locality
Every mutation is local: no global queue, retriangulation or
synchronisation. Subdomains run independently and meet at narrow interfaces.
Prefer structure-of-arrays layouts and regular access to pointer graphs and
recursive walks.

## 3. Robustness
* Filtered predicates: floats first, exact arithmetic only when the sign is
  unsafe.
* Handle the cases of `tester.md` §3A, and acute angles.
* Oracles use the producer's predicate, never its records.

## 4. Prior art
Start from the literature (`docs/increments/README.md`, step 1): Garland and
Heckbert 1995 (greedy insertion), Chew and Ruppert (Delaunay refinement),
Üngör (off-centres).

## 5. Review
A major algorithmic change states its effect on locality, parallelism and
robustness.
