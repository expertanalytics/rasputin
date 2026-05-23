# Project Structure

Layout for the post-CGAL rasputin backend. Each module corresponds to one design doc and one test directory. C++ implementation lives under `src/rasputin/_core/` (the underscore signals "private — go through the Python API"); the existing Python layer stays at `src/rasputin/*.py` and continues to expose the public API.

## Directory layout

```
src/rasputin/
  __init__.py
  mesh.py                  # public API (existing, will be rewired to new backend)
  geometry.py
  reader.py
  tin_repository.py
  ...                      # other existing Python modules stay put

  _core/                   # C++ implementation, exposed via pybind11
    CMakeLists.txt
    bindings.cpp           # pybind11 module definition (entry point)

    raster_io/             # GeoTIFF read/write, raster<T> container, NoData
    geometry_predicates/   # Shewchuk-style robust 2D predicates
    parallel_util/         # work distribution, atomics, thread pool

    vector_simplify/       # Visvalingam-Whyatt, Douglas-Peucker, topology checks
    hydrology/             # pit fill, flow direction, accumulation, catchments, streams
    noding/                # snap rounding, PSLG construction
    cdt/                   # thin wrapper over Detria (or poly2tri)
    mesh/                  # triangle/vertex data structures, ternary tree, edge tags
    refinement/            # adaptive ternary-tree refinement
    flip/                  # final Lawson edge-flip pass, constraint-respecting

cpp_test/                  # C++ tests (Catch2 + rapidcheck)
  unit/                    # per-module unit tests
  property/                # invariant / property-based tests
  integration/             # cross-module pipeline tests
  data/                    # tiny in-repo fixtures used by C++ tests

tests/                     # Python tests (pytest + hypothesis)
  unit/
  integration/
  fixtures/                # tiny synthetic DEMs, expected outputs
  conftest.py

test_data/                 # larger real-world fixtures
  README.md                # explains LFS / download script
  download.py              # fetches and caches DEM fixtures

lib/                       # bundled header-only third-party (Detria, etc.) if not system-installed

# Existing top-level docs
parallel_refinement.md
auto_catchments.md
project_structure.md       # this file
testing.md
```

## Module responsibilities and dependencies

```
                    raster_io ──────────────┐
                                            ▼
geometry_predicates ──┬─→ vector_simplify   hydrology  (depends on raster_io,
                      │                     │           parallel_util)
                      └─→ noding ←──────────┘   (catchment outline, river polylines)
                              │
                              ▼
                            cdt   (thin wrapper over Detria or poly2tri,
                              │    consumes noded PSLG, produces initial mesh)
                              ▼
                            mesh   (triangle/vertex data; ternary tree; edge tags)
                              │
                              ├──→ refinement   (raster_io for sampling, parallel_util)
                              │
                              └──→ flip         (geometry_predicates, parallel_util)
                                      │
                                      ▼
                                   bindings.cpp ──→ Python API
```

Cardinal rule: dependencies flow downward in this diagram; no upward dependencies (e.g. `mesh` does not depend on `refinement` even though `refinement` produces meshes — the refinement code is built on top of the mesh data structures).

### `raster_io`

GeoTIFF read/write, `Raster<T>` container with projection and NoData, sample-point iteration, raster-cell bbox queries. Foundational; everything that touches DEMs goes through here.

### `geometry_predicates`

Shewchuk's adaptive `orient2d` and `incircle` and friends (public domain — drop-in). Used by `noding`, `cdt`, and `flip`. Header-only.

### `parallel_util`

Wrappers over `<thread>`, `<atomic>`, `std::execution` (or a small work-stealing pool). Compaction primitives. Parallel sort/unique wrappers. No external runtime dependency by default; TBB or OpenMP can be opted in via CMake flag.

### `vector_simplify`

Visvalingam-Whyatt for area-preserving polygon simplification; Douglas-Peucker for polyline length-preserving; multi-feature topology checks (no introduced crossings). Pre-noding step.

### `hydrology`

Implements `auto_catchments.md`: pit filling (priority-flood with epsilon, plus optional Lindsay hybrid), D8 flow direction, parallel flow accumulation (Barnes-Lehman-Mulla), catchment BFS, stream extraction with Strahler ordering, mask polygonization. May depend on RichDEM (MIT) as either reference or external library — decision deferred to first implementation pass.

### `noding`

Implements the constraint-noding section of `parallel_refinement.md`. Raster-grid broad phase, robust pairwise intersection, snap rounding, segment splitting, deduplication. Outputs a clean PSLG with one bit per edge (`is_river`).

### `cdt`

Thin wrapper around the chosen MIT/BSD CDT library (Detria preferred, poly2tri fallback). Translates between rasputin's PSLG representation and the library's API. Used exactly once per mesh build.

### `mesh`

Core data structures: vertex array, triangle array, ternary tree of refinement nodes, per-edge constraint bitmask (currently just `is_river`). Owns the flatten-to-final-mesh step. Used by `refinement`, `flip`, and bindings.

### `refinement`

Implements the refinement loop in `parallel_refinement.md`. Active-list management, per-triangle DEM-sample scan with max-error tracking, fan subdivision at the max-error point, round-by-round compaction.

### `flip`

Final Lawson edge-flip pass. Skips constraint-tagged edges. Parallel with edge-conflict resolution.

### `bindings.cpp`

Single pybind11 module that exposes the C++ API to Python. Built as `rasputin/_core.so` (or similar) and imported by the Python layer.

## Build system

CMake-based, building one shared library + one Python extension. C++20 required (matches existing project). Highlights:

- **Top-level `CMakeLists.txt`** orchestrates the build, finds dependencies, configures sanitizer flags via `RASPUTIN_SANITIZER=asan|ubsan|tsan|none`.
- **`src/rasputin/_core/CMakeLists.txt`** builds each module as an `OBJECT` library, links them into the final pybind11 extension.
- **Detria** (header-only) lives in `lib/detria/` if not system-installed; CMake searches both.
- **External deps under consideration:** RichDEM (optional, MIT), GDAL (for GeoTIFF — already in use indirectly), Eigen (if linear algebra needs grow beyond what we want to hand-roll).
- **No CGAL, no Boost.Geometry** in the new core. Existing Python-layer uses of these are migrated incrementally.

`setup.py` continues to drive the build via the existing `CMakeBuild` class (per `CLAUDE.md`).

## Python API surface

The public Python API (`rasputin.mesh.Mesh`, `rasputin.geometry.Geometry`, etc.) stays stable. Underneath, `Mesh.from_raster()` and friends are rewired to call into `rasputin._core` instead of the old `triangulate_dem` module. Migration can be done one entry point at a time; the old CGAL-backed module can coexist behind a feature flag during the cutover.

## What gets deleted, eventually

After the new backend ships and the Python API is rewired:

- `src/rasputin/triangulate_dem.h` / `bindings.cpp` (the CGAL-based originals)
- CGAL, GMP, MPFR from the CMake dependency list
- Boost.Geometry, if no longer used after vector_simplify is in-tree

These removals are not part of the initial build-out; they're a follow-up once feature parity is reached and tests pass on the new backend.

## Existing files to integrate

The repo already has:

- `src/rasputin/triangulate_dem.h` (CGAL-based, to be replaced)
- `src/rasputin/xml_reader.py`, `xml_reader_new.py` (untracked — pending classification)
- `lib/date/` (existing third-party — stays)
- `src/rasputin/smiley.svg` (presumably for web visualization — stays)

The new `_core/` tree is added alongside the existing C++ files; the old files stay buildable until the new pipeline reaches parity, then are removed in a single cleanup commit.
