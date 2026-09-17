# Testing Strategy

How the rasputin backend stays correct under refactors, across compilers, and at scale.

The goal is "water tight": no module ships without invariant tests, every algorithm has a baseline of synthetic + procedural + real inputs, and CI exercises the suite serial and in parallel under sanitizers.

> **Status.** This document is mostly a target, not a description. Most of the
> modules it specifies (`noding`, `cdt`, `refinement`, `flip`, `vector_simplify`,
> `hydrology`) do not exist yet; the post-CGAL core currently ships geometry
> primitives and raster sampling. Sections below are marked **[live]** where CI
> enforces them today and **[planned]** where they do not. Read a planned
> section as a commitment about what the gate must become before the module it
> governs lands — not as a claim about what runs now. A testing doc that
> overstates its own coverage is worse than none, because it is trusted.
>
> Live today (`.github/workflows/main.yaml`): Catch2 C++ suites via ctest on
> ubuntu and macos, an asan+ubsan Debug build of the same suites, pytest
> across Python 3.11-3.13 with an enforced 85% line coverage floor, mypy
> strict, ruff, and the governance gates in `tools/`.

## Three data tiers [planned]

### Tier 1 — Synthetic primitives (in-code)

No on-disk data. Each test constructs its own input: a single triangle, a hand-built polygon, an analytical surface evaluated at a few points. The expected output is known by construction.

This tier dominates unit tests. It's where edge cases live — three-way intersections, single-pixel pits, near-degenerate triangles, segments coincident after snap rounding — because real data hits these by chance only.

- Fast: thousands of tests per second.
- Deterministic.
- No I/O.
- Bug reproductions land here as new in-code tests.

### Tier 2 — Procedural DEMs

Analytical surfaces sampled to raster:

| Surface              | Used to test                                       |
|----------------------|----------------------------------------------------|
| Tilted plane         | Flow direction is uniform; refinement converges in zero rounds. |
| Paraboloid bowl      | Pit fill removes the single pit; accumulation peaks at center. |
| Gaussian peak        | Catchment delineation around the peak; stream extraction below threshold. |
| Sinusoidal ripple    | Refinement converges with predictable triangle count vs tolerance. |
| Step function (ridge)| Snap-rounding behavior on near-coincident features. |
| Y-shaped valley      | Stream network with known Strahler structure.      |

Parameterized by seed, size, and amplitude. Tiny by default (32×32, 128×128) so tests stay sub-second. The bridge between "I constructed this" and "this resembles a raster."

### Tier 3 — Real DEM fixtures

2–3 small clipped regions (~1 km² each), chosen for diverse characteristics:

| Fixture                          | Tests                                              |
|----------------------------------|----------------------------------------------------|
| Mountain ridge (Norway)          | Sharp features, refinement depth, stream branching |
| Coastal region                   | DEM-boundary handling, truncated catchments        |
| Lake-containing valley           | External lake polygon integration                  |
| Flat agricultural / urban        | Pit-fill behavior on naturally flat areas          |

**Source:** Kartverket (Norwegian Mapping Authority, open data) given the project context. Fall back to SRTM for international reproducibility if needed.

**Storage:** check in if total under ~10 MB; otherwise Git LFS or a download script (`test_data/download.py`) with content-addressed cache. CI fetches once and caches.

Used for: integration tests, performance baselines, visual-regression spot-checks (render the output mesh and diff against a stored PNG within a perceptual tolerance).

## Invariant catalog [partly live]

These are the "water tight" properties — they must hold on every input regardless of tier. Sections marked **[live]** describe modules that exist and are enforced by the suites under `tests/cpp/`; **[planned]** sections describe modules that do not exist yet. Property tests run them against generators that produce random inputs in each module's domain.

### `predicates` [live]

- `orient2d(a,b,c) == reversed(orient2d(a,c,b))` and `orient2d(a,b,c) == orient2d(b,c,a)`, for every triple, unconditionally — including collinear and cocircular input. Two callers asking about the same triple in different vertex orders must not disagree.
- `FilteredKernel<E>` returns exactly what `E` returns, on every finite input. The filter is a performance device and may never change an answer.
- The exact backend is consulted iff the float determinant is inside the error bound: zero fallbacks on well-separated input, at least one on a collinear triple.
- `incircle` is total, precondition-free, and independent of the order of its first three arguments. Collinear input returns `Cocircular` without consulting the backend.
- Degeneracy is exact: three points on an exactly representable line are `Collinear`, four on an exactly representable circle are `Cocircular`.
- Correct at UTM33 magnitudes, where naive determinants have no correct digits.

### `core geometry` [live]

- `point_in_ring` is invariant under cyclic rotation of the vertices and under reversal of the ring.
- Every ring vertex classifies `Boundary`.
- `orientation` negates exactly when the ring is reversed, and is invariant under cyclic rotation.
- `orientation` agrees with `sign(signed_area)` for every ring whose area is well separated from zero — and there exists a sliver where it does not, which is why `signed_area`'s sign may never drive a topology decision.
- `|signed_area|` is invariant under rotation and reversal.
- `bounding_box(ring)` contains every vertex, and no point outside it classifies `Inside` or `Boundary`.
- `Box2::expand` as a fold is order-independent; the empty box is its identity.
- `FastKernel` and `DefaultKernel` agree on `point_in_ring` for well-separated input, and there exists a near-degenerate ring where `FastKernel` is demonstrably wrong under every FMA-contraction mode.
- On a ring built from `RasterGeometry` corner nodes, `contains_strict(p)` implies `point_in_ring(p) == Inside`.

### `core geometry — PSLG` [live]

- A `Pslg` exists only if its build produced an empty diagnostics list; there is no partially valid `Pslg`.
- Every index in `chain_indices()` is `< vertices().size()`; every coordinate in `vertices()` is finite, including unreferenced ones.
- Every `Outer` ring is counterclockwise and every `Hole` ring is clockwise under the validating kernel; neither is collinear.
- No closed chain stores its closure; breaklines are exempt and may be closed polylines.
- `indices_of(i)` sub-spans partition `chain_indices()` contiguously in chain order, with no gap and no overlap.
- The vertex buffer is element-wise equal to the builder's input: no dedup, no reordering, no reversal.
- `ring(c)` never throws for a closed chain.
- Validation is exhaustive: N independently broken chains produce at least N diagnostics. Stage 0 is the one exception and returns early, because past truncation every later diagnostic is noise.
- A `const Pslg` is safe for concurrent read; no accessor mutates or caches.
- Negative, and equally load-bearing: a valid `Pslg` promises **no** simplicity, **no** pairwise disjointness and **no** nesting. Any test asserting one of those is testing the noder and belongs in its catalog.

### `noding` [planned]

- No two output segments intersect in their interior.
- Every input vertex appears in the output vertex set (post-snap).
- Every output vertex is either an input vertex (snapped) or lies on at least two input segments.
- Sum of output segment lengths equals sum of input segment lengths, modulo snap perturbation bounded by the snap-grid spacing.
- `is_river` bit on any output segment is the OR of the bits on the input segments that contributed.

### `cdt`

- Every constraint segment from the noded PSLG appears in the output as a connected chain of triangulation edges.
- The unconstrained sub-triangulation is Delaunay (incircle test passes for every interior edge).
- No degenerate triangle (zero or near-zero area below epsilon).
- Triangle count matches Euler's formula given vertex count and boundary.

### `refinement`

- Every leaf triangle's max sample error ≤ tolerance.
- Every constraint edge from the initial CDT survives intact (a constraint edge is either still present or has been subdivided along the constraint geometry — it has not been crossed or replaced).
- Every inserted vertex lies strictly inside its parent triangle (centroid-fan property).
- Mesh remains conforming after each round (no hanging nodes).

### `flip`

- No constraint edge is flipped.
- After convergence, every non-constraint edge is locally Delaunay.
- Vertex set is unchanged from the input mesh; only edges move.
- Triangle count is unchanged.

### `vector_simplify`

- Visvalingam: simplified polygon's area differs from original by less than a configurable tolerance.
- Douglas-Peucker: every original vertex is within the distance tolerance of the simplified polyline.
- Multi-feature: simplified outputs do not introduce new intersections that the inputs did not have.

### `hydrology` — pit filling

- No pit remains in output (every non-boundary pixel has at least one neighbor of strictly lower elevation, modulo epsilon).
- Output elevation ≥ input elevation everywhere (filling never lowers).
- Boundary pixels are unchanged.

### `hydrology` — flow accumulation

- Mass conservation: Σ accumulation = total contributing pixel count.
- Source pixels (no upstream neighbors) have accumulation 1.
- Boundary outlets have flow direction pointing off-domain.

### `hydrology` — catchment

- Output is connected.
- Every interior pixel's D8-traced path eventually reaches the seed.
- No exterior pixel's path reaches the seed.
- Boundary-touch flag is set iff any catchment pixel sits on a DEM boundary.

### `hydrology` — stream extraction

- Output graph is a tree (or forest, if multiple outlets).
- Strahler order is monotonic non-decreasing from headwater to outlet.
- A confluence of two order-`k` streams produces an order-`k+1` segment downstream.

## Frameworks [partly live]

- **Catch2 v3** for C++ unit and integration tests (in `tests/cpp/`, fetched via `FetchContent`).
- **Catch2 `GENERATE` over a seeded range** for C++ property-based tests.
  rapidcheck was considered and declined; each property suite records why at the
  top of the file, the operative reason being that standing up a second
  framework inside the PR that introduces a module is how neither gets done.
- **pytest** for Python tests (existing).
- **hypothesis** for Python property-based tests.
- **pytest-benchmark** for performance regression tracking on tier-3 fixtures.

Property test generators live alongside the modules they test (e.g. `tests/cpp/property/noding_generators.h` produces random sets of polylines with controllable density of intersections).

## Parallelism testing [planned]

**Run the entire suite at thread counts 1, 2, 8, and a CI-runner-max value.** For deterministic algorithms, results must be bit-identical across thread counts; for non-deterministic ones (e.g. order of insertions affects which equally-valid Delaunay flip is chosen), assert equivalence under a defined norm.

**Sanitizers in CI:**

- **[live]** asan + ubsan Debug build run on every PR (`sanitizers` job in
  `.github/workflows/main.yaml`). Built with `-fno-sanitize-recover=all`,
  which is load-bearing: UBSan's default is to print the diagnostic and
  continue, so the process exits 0 and the job passes green on undefined
  behaviour.
- **[live]** `DefaultKernel::incircle` is exercised against non-counterclockwise
  input in the Debug job, unguarded and in process
  (`tests/cpp/unit/test_predicates_default_kernel.cpp`, plus the property sweep
  in `tests/cpp/property/prop_predicates_detria_agreement.cpp`). This is a
  blocking gate for the predicates module, and it is the one test here that
  passes by *not crashing*. The vendored detria backend asserts its
  counterclockwise precondition under `#ifndef NDEBUG` and its assert handler
  calls `std::raise(SIGTRAP)`, which is not a signal Catch2 installs a handler
  for -- so if `FilteredKernel::incircle`'s normalization ever regresses, the
  asan+ubsan Debug job dies on a signal with no assertion text rather than
  reporting a failure. A forked probe in the same file exists to turn that death
  into a readable report; the unguarded test exists to make sure the real,
  unmediated call path is what is being exercised. Release builds define NDEBUG
  and cannot see this at all, which is why it is listed with the sanitizers.
- **[planned]** tsan build on every PR — catches races in parallel code
  paths. Nothing to race yet; due when the parallel refinement lands.
- **[planned]** msan (MemorySanitizer, clang-only) nightly — catches
  uninitialized reads. Needs an instrumented libc++ to avoid false
  positives, so it is the most expensive of the three to stand up.

**Stress tests:** random thread interleavings via `std::this_thread::yield()` insertions in debug builds, run nightly on tier-2 procedural DEMs. Useful for surfacing races that TSan misses.

## CI matrix [partly live]

One row of this matrix is configured. What runs today is two Release C++
builds (ubuntu, macos), the asan+ubsan Debug job, and three Python versions.
Still absent: the compiler matrix, tsan, the thread-count sweep and everything
nightly.

Target, per PR (marked against what exists):

| Compiler | Build type       | Sanitizer       | Thread count | Status |
|----------|------------------|-----------------|--------------|--------|
| gcc      | Release          | none            | max          | **live** (plus macos/clang Release) |
| gcc      | Debug            | asan + ubsan    | 1            | **live** |
| clang    | Debug            | tsan            | 8            | planned |
| clang    | RelWithDebInfo   | none            | 1            | planned |

Nightly adds:

- clang + msan + debug + serial
- Coverage build (gcov) and report
- Performance benchmark against tier-3 fixtures with regression detection (>5% slowdown fails)

## Coverage targets [partly live]

- **Line coverage ≥ 85%**, enforced project-wide by `--cov-fail-under=85`. Per-module is a review obligation, not machine-checked. Anything lower needs a justification comment in the PR.
- **Invariant coverage**: every documented invariant has at least one test naming it.
- **Edge-case coverage**: every condition handled specially in code (NaN, NoData, empty input, single-element input, boundary intersection, etc.) has a named test.

Python line coverage is enforced today: `--cov-fail-under=85` in pyproject's
`addopts`, so a drop below the floor fails the run rather than being noticed in
review. **C++ coverage is not yet measured** -- `gcovr` is not configured and no
nightly coverage build exists, so the per-module figure above is currently a
Python-only guarantee. Invariant and edge-case coverage are review obligations,
not machine-checked ones.

## Test layout conventions [partly live]

```
tests/cpp/
  unit/
    test_noding_snap_rounding.cpp        # one file per logical concept
    test_noding_pairwise_intersection.cpp
    test_refinement_fan_subdivision.cpp
    ...
  property/
    prop_noding_no_crossings.cpp         # one file per invariant
    prop_refinement_tolerance.cpp
    prop_flip_constraint_preservation.cpp
    ...
  integration/
    int_full_pipeline_paraboloid.cpp     # tier-2 end-to-end
    int_full_pipeline_mountain.cpp       # tier-3 end-to-end
    ...
  data/
    *.tif, *.json                        # tiny fixtures used by C++ tests

tests/
  unit/
    test_mesh_python_api.py
    ...
  integration/
    test_pipeline_end_to_end.py
  fixtures/
    paraboloid_64.tif
    expected_outputs/*.json
```

Test names start with `test_` for example-based and `prop_` for property-based, so a grep for the prefix lists the suite.

## Regression handling [planned]

When a bug is fixed:

1. Add a tier-1 (synthetic) test that fails on the broken code and passes on the fix. Name it after the issue or symptom.
2. If the bug came from real data, also clip the smallest possible region triggering it into `test_data/` as a tier-3 fixture with a comment explaining what it reproduces.
3. If the bug exposed a missing invariant, add the invariant to this doc and write a property test for it.

The point is to never lose a bug twice. Tier-1 tests are cheap and stay in the suite forever.

## What we do not test

- **Performance correctness of third-party libraries** (Detria, RichDEM if used) — we trust their own test suites and only test our wrappers.
- **Third-party TIFF container decoding** — we trust the chosen pure-Python
  reader's own test suite. GeoTIFF *georeferencing* is emphatically ours to
  test: GDAL is prohibited (see `CLAUDE.md` section 2), so tie-point and
  pixel-scale interpretation, GeoKey decoding and CRS construction all have to
  be covered here.
- **Visual aesthetics** of the output mesh beyond the perceptual-diff regression on tier-3 fixtures.

These are out of scope and would inflate the suite without buying confidence in rasputin's own code.
