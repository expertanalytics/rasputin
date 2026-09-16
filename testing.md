# Testing Strategy

How the rasputin backend stays correct under refactors, across compilers, and at scale.

The goal is "water tight": no module ships without invariant tests, every algorithm has a baseline of synthetic + procedural + real inputs, and CI exercises the suite serial and in parallel under sanitizers.

## Three data tiers

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

## Invariant catalog

These are the "water tight" properties — they must hold on every input regardless of tier. Property tests run them against generators that produce random inputs in each module's domain.

### `noding`

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

## Frameworks

- **Catch2 v3** for C++ unit and integration tests (in `tests/cpp/`, fetched via `FetchContent`).
- **rapidcheck** for C++ property-based tests, integrated as a Catch2 extension.
- **pytest** for Python tests (existing).
- **hypothesis** for Python property-based tests.
- **pytest-benchmark** for performance regression tracking on tier-3 fixtures.

Property test generators live alongside the modules they test (e.g. `tests/cpp/property/noding_generators.h` produces random sets of polylines with controllable density of intersections).

## Parallelism testing

**Run the entire suite at thread counts 1, 2, 8, and a CI-runner-max value.** For deterministic algorithms, results must be bit-identical across thread counts; for non-deterministic ones (e.g. order of insertions affects which equally-valid Delaunay flip is chosen), assert equivalence under a defined norm.

**Sanitizers in CI:**

- `RASPUTIN_SANITIZER=asan,ubsan` build run on every PR.
- `RASPUTIN_SANITIZER=tsan` build run on every PR — catches races in parallel code paths.
- `RASPUTIN_SANITIZER=msan` (MemorySanitizer, clang-only) run nightly — catches uninitialized reads.

**Stress tests:** random thread interleavings via `std::this_thread::yield()` insertions in debug builds, run nightly on tier-2 procedural DEMs. Useful for surfacing races that TSan misses.

## CI matrix

Per PR:

| Compiler | Build type       | Sanitizer       | Thread count |
|----------|------------------|-----------------|--------------|
| gcc      | Release          | none            | max          |
| gcc      | Debug            | asan + ubsan    | 1            |
| clang    | Debug            | tsan            | 8            |
| clang    | RelWithDebInfo   | none            | 1            |

Nightly adds:

- clang + msan + debug + serial
- Coverage build (gcov) and report
- Performance benchmark against tier-3 fixtures with regression detection (>5% slowdown fails)

## Coverage targets

- **Line coverage ≥ 85%** per module. Anything lower needs a justification comment in the PR.
- **Invariant coverage**: every documented invariant has at least one test naming it.
- **Edge-case coverage**: every condition handled specially in code (NaN, NoData, empty input, single-element input, boundary intersection, etc.) has a named test.

`gcovr` produces HTML reports; failed coverage thresholds fail CI.

## Test layout conventions

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

## Regression handling

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
