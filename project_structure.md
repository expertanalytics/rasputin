# Project Structure

Layout for the post-CGAL rasputin backend. Each module corresponds to one design doc and one test directory. Public C++ headers live under `include/terrain/`, implementation under `src/`, and the pybind11 module under `bindings/`. The Python layer is the `tin_engine` package under `src_python/` (the distribution is still named `rasputin`); the extension is imported as `tin_engine._core`, the underscore signalling "private — go through the Python API".

## Directory layout

Entries marked *(planned)* do not exist yet; the rest are in the tree today.

```
include/terrain/           # public C++ headers, header-only where possible
  core/
    point.hpp              # Point2 / Point3 value types
  raster/
    geometry.hpp           # RasterGeometry, CellIndex: grid <-> world mapping
    raster.hpp             # RasterSource concept, owning Raster<T>, NoData
    sample.hpp             # bilinear interpolation over a RasterSource

src/                       # C++ implementation, one directory per module (planned)
  geometry_predicates/     # Shewchuk-style robust 2D predicates
  parallel_util/           # work distribution, atomics, thread pool
  vector_simplify/         # Visvalingam-Whyatt, Douglas-Peucker, topology checks
  hydrology/               # pit fill, flow direction, accumulation, catchments, streams
  noding/                  # snap rounding, PSLG construction
  cdt/                     # thin wrapper over Detria (or poly2tri)
  mesh/                    # triangle/vertex data structures, ternary tree, edge tags
  refinement/              # adaptive ternary-tree refinement
  flip/                    # final Lawson edge-flip pass, constraint-respecting

bindings/
  core.cpp                 # pybind11 module definition -> tin_engine._core
  legacy_bindings.cpp      # CGAL-based original; not built, kept for reference

src_python/tin_engine/     # public Python API (distribution name: rasputin)
  __init__.py              # re-exports from tin_engine._core
  cli.py                   # Typer entry point declared in pyproject (planned)

tests/
  cpp/                     # C++ tests (Catch2; rapidcheck planned)
  python/                  # Python tests (pytest; hypothesis planned)
  fixtures/                # in-repo test data (DEM, GML, textures, TIN archives)

legacy/                    # archived pre-migration tree, not built
  rasputin/

lib/                       # bundled header-only third-party (Detria, etc.) (planned)

# Existing top-level docs
parallel_refinement.md
auto_catchments.md
project_structure.md       # this file
testing.md
```

## Module responsibilities and dependencies

```
                    raster ────────────────┐
                                            ▼
geometry_predicates ──┬─→ vector_simplify   hydrology  (depends on raster,
                      │                     │           parallel_util)
                      └─→ noding ←──────────┘   (catchment outline, river polylines)
                              │
                              ▼
                            cdt   (thin wrapper over Detria or poly2tri,
                              │    consumes noded PSLG, produces initial mesh)
                              ▼
                            mesh   (triangle/vertex data; ternary tree; edge tags)
                              │
                              ├──→ refinement   (raster for sampling, parallel_util)
                              │
                              └──→ flip         (geometry_predicates, parallel_util)
                                      │
                                      ▼
                                   bindings.cpp ──→ Python API
```

Cardinal rule: dependencies flow downward in this diagram; no upward dependencies (e.g. `mesh` does not depend on `refinement` even though `refinement` produces meshes — the refinement code is built on top of the mesh data structures).

### `raster`

`RasterGeometry` (north-up, grid-registered index/world mapping), the
`RasterSource` concept, an owning `Raster<T>` with NoData, and bilinear
sampling. Foundational; everything that touches DEMs goes through here.

Named `raster`, not `raster_io`: it performs no I/O. Header-only and
templated, so `src/` holds nothing for it. Landed in `7785fea`, which also
fixed two legacy defects — an out-of-bounds bilinear read and a transposed
row/column index.

Still to come here: `window_for` (bbox to index window), deferred until
`refinement` gives it a caller, and `RasterView` over a numpy buffer, which
belongs with the reader.

**Open decision:** GeoTIFF *decoding* is not sited yet. `pyproject.toml` and
`.claude/skills/geospatial-data-formats/SKILL.md` put a pure-Python reader
under `src_python/tin_engine/io/`; this section previously assigned read/write
to C++. Decoding in Python keeps the C++ core testable with no I/O, but the
call belongs to `@architect`.

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

### `bindings/core.cpp`

Single pybind11 module that exposes the C++ API to Python. Built as the `_core` extension, installed into the `tin_engine` package, and re-exported by `src_python/tin_engine/__init__.py`. Currently binds the `Point2`/`Point3` value types; `__repr__` routes through the `std::formatter` specializations in `point.hpp` so the C++ and Python renderings cannot drift.

## Build system

CMake-based, building the header-only core plus one Python extension. C++20 required. Highlights:

- **Top-level `CMakeLists.txt`** orchestrates the build. `terrain_headers` is an `INTERFACE` target carrying `include/`. Two options gate the rest: `RASPUTIN_BUILD_PYTHON` (default `OFF`) adds the pybind11 extension, `RASPUTIN_BUILD_TESTS` (default `ON`) adds the C++ suite.
- **`tests/cpp/CMakeLists.txt`** fetches Catch2 v3 via `FetchContent`, so it needs no manual checkout.
- **Detria** (header-only) lives in `lib/detria/` if not system-installed; CMake searches both.
- **External deps under consideration:** RichDEM (optional, MIT), Eigen (if linear algebra needs grow beyond what we want to hand-roll).
- **No CGAL, no Boost.Geometry, no GDAL** in the new core. Existing Python-layer uses are migrated incrementally.
- **Planned:** per-module `OBJECT` libraries linked into the extension once `src/` is populated, and sanitizer flags via `RASPUTIN_SANITIZER=asan|ubsan|tsan|none`.

Packaging is driven by **scikit-build-core**, declared in `[build-system]` in `pyproject.toml`, which invokes this same CMake build. `pip install .` configures with `RASPUTIN_BUILD_PYTHON=ON` and `RASPUTIN_BUILD_TESTS=OFF` — the C++ tests pull Catch2 over the network and have no business running during an install — and installs `_core` into the `tin_engine` package.

There is no `setup.py`. It was removed with the foundation reset, along with the `CMakeBuild` class it carried.

## Python API surface

The public API is the `tin_engine` package, calling into `tin_engine._core`. The pre-migration `rasputin.*` modules (`mesh.py`, `geometry.py`, `reader.py`, `tin_repository.py` and friends) are archived under `legacy/rasputin/` rather than kept in place, so this is a re-implementation against the new backend rather than a rewiring of stable modules. Porting proceeds one entry point at a time; shapes worth preserving should be read out of `legacy/` before being reintroduced.

## What gets deleted, eventually

After the new backend ships and the Python API is rewired:

- `legacy/rasputin/triangulate_dem.h` and `bindings/legacy_bindings.cpp` (the CGAL-based originals)
- CGAL, GMP, MPFR from the CMake dependency list — already absent from the new `CMakeLists.txt`
- Boost.Geometry, if no longer used after vector_simplify is in-tree

These removals are not part of the initial build-out; they're a follow-up once feature parity is reached and tests pass on the new backend.

## Existing files to integrate

The repo already has:

- `legacy/rasputin/triangulate_dem.h` (CGAL-based, to be replaced)
- `legacy/rasputin/*.py` (the pre-migration Python layer, pending per-module classification)
- `lib/date/` — **removed**. The C++20 `<chrono>` calendar types replaced it; see the `lib/date` note in `CLAUDE.md` section 2.

The new `_core/` tree is added alongside the existing C++ files; the old files stay buildable until the new pipeline reaches parity, then are removed in a single cleanup commit.
