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
> across Python 3.12-3.14 with an enforced 85% line coverage floor, mypy
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
| Step function (ridge)| Refinement behaviour across a sharp elevation discontinuity. |
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

These are the "water tight" properties — they must hold on every input regardless of tier. Sections marked **[live]** describe modules that exist and are enforced by the suites under `tests/cpp/` — or, for the Python `viz` section, under `tests/python/`; **[planned]** sections describe modules that do not exist yet. Property tests run them against generators that produce random inputs in each module's domain.

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

The list below replaces the five bullets that stood here until increment 5b.
Three of them were false, and `docs/increments/05-noder.md` gives the
measurement for each. A fourth — the feature-property bullet — was already
corrected upstream by increment 7, from a one-bit `is_river` OR to a union over
property sets; 5b keeps that correction and sharpens where the merge happens and
what the oracle may be built from. The marker stays `[planned]` until 5c merges: 5b ships a noder that
no production code calls.

- **No two output edges cross in their interiors** — guarantee 14(a). It is *verified* by a second pass over the output, not established by construction, because snap rounding can create a crossing that was not in the input. Two output edges with **equal** node-id pairs are permitted and are what the property union below exists for; a partial collinear overlap is not, and after splitting there is no third case.
- **No node's cell meets an edge it is not an endpoint of** — guarantee 14(b), and the half a reader will assume follows from the one above. It does not: this is the hot-pixel form, `segment_meets_cell`, not exact incidence, and it is strictly the stronger of the two. A T-junction the split pass missed is caught here rather than blessed, which an exact-incidence spelling would not do.
- Every input vertex's **snapped image** is an output node, and `node_of_input_vertex` maps it there. Not "every input vertex appears in the output": after dedup, two input vertices in one cell produce **one** output node.
- Every output vertex is either the snapped image of an input vertex or the snapped image of a **constructed** crossing of two input segments. Not "lies on at least two input segments": commit `e412a43` establishes that a constructed crossing point need not lie on *either* segment after snapping.
- **One-way Hausdorff, not a length sum**: every output edge derived from input segment `e` lies within the closed `h/√2`-neighbourhood of `e`, for grid spacing `h`, and within `k·h/√2` after `k` noding rounds. The sum-of-lengths form that stood here was false twice over — a collinear overlap merge deletes length outright, and a segment split at `m` points accumulates up to `m` displacements, so the bound is not per segment. The `√2` arithmetic is right and is kept: a vertex moves by at most half a cell diagonal.
- The property set on any output edge is the **union** over the input **chains** that contributed to it — guarantee 15. Union rather than a one-bit OR: at a coarse resolution one segment can be both a road and a river, so an edge carries a set (`parallel_refinement.md`, "Edge metadata", and `docs/increments/07-edge-properties.md` for the type). The merge happens at the node-id edge-key dedup, not at a collinear-overlap classification. It is verified against an oracle built from the **input**, never from the dedup's own provenance map, which would only restate the dedup — and asserted as a **subset**, not an equality, because a chain passing near both nodes satisfies the contribution relation without having contributed.
- No two output vertices are equal; no output edge is zero-length; every output coordinate satisfies `grid.snapped(v) == v` bitwise; output chain order and roles match the input chain-for-chain.
- Node ids are independent of input order and of thread count.

### `cdt` [live]

- Every **in-domain** constraint edge appears as an edge of one or two output triangles, and is flagged in the constrained-edge mask of each. A breakline outside the outer ring or inside a hole is legal in a `Pslg` and appears in no interior triangle — measured, and pinned by its own fixture.
- **Delaunay with respect to visibility**, which is the only form true of a CDT: for every interior non-constrained edge `(a,b)` with apexes `c,d`, either `incircle(a,b,c,d) != Inside` **or** the open segment `cd` crosses at least one constraint edge. The plain incircle form this catalog previously stated is false for a constrained triangulation and would fail on correct output.
- Every output triangle is counterclockwise under `DefaultKernel` — which is the strongest true statement, and what the property suite asserts. It implies no zero-area triangle, stated exactly rather than "below epsilon": this project has no epsilon.
- Bit `e` of a triangle's constrained-edge mask is set iff edge `(v[e], v[(e+1)%3])` is a constraint. CGAL's opposite-vertex convention is a rotation of this one, and the two agree on **exactly** the masks `0b000` and `0b111` — so every triangle with one or two constrained edges distinguishes them, and only a triangle with zero or three **fails** to. (An earlier revision of this line had the complement, claiming only a one-bit mask distinguishes them. Re-runnable: `python3 -c "A=lambda m:{(e,(e+1)%3) for e in range(3) if m>>e&1}; B=lambda m:{((e+1)%3,(e+2)%3) for e in range(3) if m>>e&1}; print([m for m in range(8) if A(m)==B(m)])"` prints `[0, 7]`.) A one-bit mask is still the sharpest case and the one `tests/python/test_viz_scene.py` builds its fixture from, because under it the rotated read lands on a *different single* edge rather than on an overlapping pair.
- No outcome carries a mesh with a non-`Ok` status, and none carries an empty mesh with `Ok`. The second is the silent mode: "forgot to add the outline" triangulates successfully to zero interior triangles.
- `triangles == 2n - b - 2 + 2h` for `n` referenced vertices, `b` on any boundary ring and `h` disjoint holes — meaningful only because the wrapper returns in-domain triangles (`forEachTriangle`, not the hole or convex-hull variants). Rings that touch break the formula: a corner-touching hole predicts 6 against a measured 5, so that fixture carries an explicit count.

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

### `viz` [live]

The only Python section of this catalog, and the only one whose module is pure
computation over data handed across the pybind11 boundary rather than C++.
`tests/python/test_viz_scene.py` **[live]** is increment 6b-i's sole
invariant-critical suite (`docs/increments/06-cdt-viewer.md`, "What is worth
testing"); it runs against hand-built fakes of `MeshLike`/`PslgLike` with no
compiled extension in the process, which is what `viz/`'s protocol boundary
buys. `tests/python/test_viz_svg.py` and `tests/python/test_cli_draw.py` are
6b-ii's and now **[live]** — ordinary rather than invariant-critical, because a
renderer's failure mode is a wrong-looking picture and a person catches that
instantly. What they pin is structure and arithmetic: that the document parses
as XML, that element counts equal the scene's primitive counts, that the
viewport flips y and preserves one scale, that a `class` token puts each edge in
the right stroke class, that no failure presentation is a blank page, and that
every gallery fixture is the shape its design row claims. `test_cli_draw.py`
adds the composition root's own obligations — the role mapping, `closed_roles`,
passthrough of the engine's words, the `--labels` refusal and path validation.

Invariants over `build_scene` — all `[live]`:

- Bit `e` of a triangle's mask names the edge `(v[e], v[(e+1)%3])`, matching `indexed_mesh.hpp`. The mutant this suite exists for is the rotation to CGAL's opposite-vertex convention, which draws a plausible picture with every constraint stroke on the wrong edge; see the `cdt` section above for why a one-bit mask is the sharpest fixture for it.
- Every undirected edge is emitted exactly once, with endpoints canonically ordered `a < b`, and in a deterministic order. An interior edge shared by two triangles is neither drawn twice nor double-classified, and a mask bit set on only one of the two incident triangles survives.
- The scene **records** mask-versus-chain disagreement rather than reconciling it: a masked edge matching no chain, and a chain edge in no mask, each produce their own finding and are each still drawn. These are two independently derived answers to "this edge is a constraint", which is the only reason their agreement means anything.
- The role join is the input's verdict and stays separate from the mask's: an edge in two chains takes the first chain's role, and its river bit is the OR over every contributing chain.
- Ring closure is the caller's policy, not the scene's (`closed_roles`); a chain of fewer than three vertices is never closed, and a one-vertex chain is never closed into a self-loop, which would violate the `a < b` ordering above.
- Classification is total: a populated mesh, `Ok` with zero triangles, a missing mesh, and a non-`Ok` status each map to their own kind, and a non-`Ok` status wins over a mesh that — per the `cdt` section — cannot accompany it.
- No drawable mesh means no findings, and the roles are kept regardless: the PSLG-only picture is drawn in role colours, and there is no mask for the chains to disagree with.
- The bbox contains every vertex, pads a zero-extent axis about its own centre while leaving a healthy axis alone, and reports that it padded. Non-finite coordinates are refused with `ValueError` rather than poisoning the bounds silently.
- The vertices and bbox are the **mesh's** when a mesh is drawn, so a backend-introduced vertex is inside the picture; they fall back to the input PSLG's only when there is no mesh.
- `viz.scene` imports no compiled extension and no first-party module but `viz.protocols` — pinned by a test that parses that module, because it is what keeps the suite runnable without a build.

## Frameworks [partly live]

- **Catch2 v3** for C++ unit and integration tests (in `tests/cpp/`, fetched via `FetchContent`).
- **Catch2 `GENERATE` over a seeded range** for C++ property-based tests.
  rapidcheck was considered and declined; each property suite records why at the
  top of the file, the operative reason being that standing up a second
  framework inside the PR that introduces a module is how neither gets done.
- **pytest** for Python tests (existing).
- **hypothesis** for Python property-based tests.
- **pytest-benchmark** for performance regression tracking on tier-3 fixtures.

Property test generators live alongside the modules they test (e.g. `tests/cpp/property/noding_generators.h` produces random sets of polylines with controllable density of intersections). That header arrives with increment 5b, the noder's topology half; increment 5a's property suite generates points, spacings and segment pairs, not polyline sets.

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
    prop_noding_snap_invariants.cpp      # one file per invariant (increment 5a)
    prop_noding_no_crossings.cpp         # increment 5b: guarantees 14 and 15
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

**Exception, and it is load-bearing: a characterisation test of third-party
behaviour we depend on but that is not part of its documented API.** The rule
below says we trust a library's own suite and test only our wrappers. That is
right for behaviour the library promises. It is wrong for behaviour we have
merely observed and built on — detria auto-closing an open polyline is not in
its API docs, was found by reading `createConstrainedEdges`, and lets the CDT
wrapper hand `Pslg::indices_of(c)` straight to `addOutline` with no buffer. A
version bump could take it away silently. The L-shaped-ring test in
`test_cdt_detria_backend.cpp` exists to fail loudly if it does, and must not be
deleted as out of scope.

- **Performance correctness of third-party libraries** (Detria, RichDEM if used) — we trust their own test suites and only test our wrappers.
- **Third-party TIFF container decoding** — we trust the chosen pure-Python
  reader's own test suite. GeoTIFF *georeferencing* is emphatically ours to
  test: GDAL is prohibited (see `CLAUDE.md` section 2), so tie-point and
  pixel-scale interpretation, GeoKey decoding and CRS construction all have to
  be covered here.
- **Visual aesthetics** of the output mesh beyond the perceptual-diff regression on tier-3 fixtures.
- **`viz/svg.py`'s stylesheet.** No colour, stroke width or dash array is
  asserted anywhere, and a **golden-file SVG comparison is rejected outright**
  (`docs/increments/06-cdt-viewer.md`, "What is worth testing"): it would pin the
  stylesheet, fail on every cosmetic improvement, and detect only the class of
  defect a person sees instantly. What *is* asserted is the `class` token on each
  element — which stroke class an edge belongs to is structure, and it is the
  only way to check that role colouring reaches the right edges without naming a
  colour. Written down because filling this gap is the first thing someone will
  propose.

These are out of scope and would inflate the suite without buying confidence in rasputin's own code.
