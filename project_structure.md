# Project Structure

Layout for the post-CGAL rasputin backend. Each module corresponds to one design doc and one test directory. Public C++ headers live under `include/terrain/`, implementation under `src/`, and the pybind11 module under `bindings/`. The Python layer is the `tin_engine` package under `src_python/` (the distribution is still named `rasputin`); the extension is imported as `tin_engine._core`, the underscore signalling "private — go through the Python API".

## Directory layout

Entries marked *(planned)* do not exist yet; the rest are in the tree today.

```
include/terrain/           # public C++ headers, header-only where possible
  core/
    point.hpp              # Point2 / Point3 value types
    bbox.hpp               # Box2; empty-box identity, exact closed containment
    segment.hpp            # Segment2 and on_segment<K>
    ring.hpp               # Ring concept, PointRing / IndexedRing, point_in_ring<K>
    pslg.hpp               # Pslg, Chain, ChainRole — validated planar input
    pslg_builder.hpp       # PslgBuilder, PslgDiagnostic, the validator
    indexed_mesh.hpp       # IndexedMesh2: SoA mesh + per-triangle constrained mask
  predicates/
    orientation.hpp        # Orientation / Incircle vocabulary; no includes at all
    exact.hpp              # ExactPredicates concept, filter bounds
    kernel.hpp             # GeometryKernel, FastKernel, FilteredKernel<E>
    detria_exact.hpp       # DetriaExact declaration; does NOT include detria
    default_kernel.hpp     # DefaultKernel = FilteredKernel<DetriaExact>
  cdt/
    result.hpp             # CdtStatus, CdtOptions, CdtOutcome
    triangulate.hpp        # CdtBackend concept and the generic entry point
    constrained_edges.hpp  # ConstraintEdgeSet: sorted keys -> per-triangle mask
    detria_backend.hpp     # DetriaBackend declaration; does NOT include detria
  raster/
    geometry.hpp           # RasterGeometry, CellIndex: grid <-> world mapping
    raster.hpp             # RasterSource concept, owning Raster<T>, NoData
    sample.hpp             # bilinear interpolation over a RasterSource
    view.hpp               # RasterView<T>: non-owning view over a contiguous
                           #   caller-supplied buffer (planned)
    window.hpp             # window_for: bbox -> index window (planned, awaits
                           #   a caller in refinement)

src/                       # C++ implementation, one directory per module
                           #   (only predicates/ and cdt/ exist; rest planned)
  predicates/              # exact orient2d/incircle; namespace terrain::pred
  parallel_util/           # work distribution, atomics, thread pool
  vector_simplify/         # Visvalingam-Whyatt, Douglas-Peucker, topology checks
  hydrology/               # pit fill, flow direction, accumulation, catchments, streams
  noding/                  # snap rounding, PSLG construction
  cdt/                     # thin wrapper over vendored Detria
  mesh/                    # triangle/vertex data structures, ternary tree, edge tags
  refinement/              # adaptive ternary-tree refinement
  flip/                    # final Lawson edge-flip pass, constraint-respecting

bindings/
  core.cpp                 # pybind11 module definition -> tin_engine._core

src_python/tin_engine/     # public Python API (distribution name: rasputin)
  __init__.py              # re-exports from tin_engine._core
  cli.py                   # Typer entry point declared in pyproject
  raster.py                # the ONLY adapter from decoded data into _core (planned)
  _core.pyi                # type stubs for the compiled extension
  viz/                     # CDT -> SVG renderer; never imports _core
    __init__.py            # re-exports Scene, SvgStyle, build_scene, render_svg
    protocols.py           # MeshLike / PslgLike / ChainLike -- typing.Protocol
    scene.py               # build_scene(pslg, mesh=None, ...) -> Scene
    style.py               # SvgStyle: frozen Pydantic V2, canvas geometry only
    svg.py                 # (Scene, SvgStyle) -> str; the stylesheet lives here
    fixtures.py            # the eight-fixture synthetic gallery, declarative
  io/                      # all file decoding lives here (planned)
    __init__.py
    geotiff.py             # TIFF container + GeoKey decoding -> DemTile
    models.py              # Pydantic RasterMeta / DemTile

tests/
  cpp/                     # C++ tests (Catch2; unit/ and property/)
  python/                  # Python tests (pytest; hypothesis planned)
  fixtures/                # in-repo test data (DEM, GML, textures, TIN archives)

legacy/                    # archived pre-migration tree, not built
  rasputin/
  bindings.cpp             # CGAL-based original; not built, kept for reference

lib/                       # vendored third-party
  detria/                  # detria.hpp at a pinned SHA, MIT; see its README

# Existing top-level docs
parallel_refinement.md
auto_catchments.md
project_structure.md       # this file
testing.md
docs/increments/           # per-increment design records; see its README
```

## Module responsibilities and dependencies

```
                    raster ────────────────┐
                                            ▼
predicates ───────────┬─→ vector_simplify   hydrology  (depends on raster,
                      │                     │           parallel_util)
                      └─→ noding ←──────────┘   (catchment outline, river polylines)
                              │
                              ▼
                            cdt   (thin wrapper over vendored Detria,
                                   consumes a PSLG, produces core::IndexedMesh2)

                            mesh   (ternary tree; edge tags — built FROM an
                              │     IndexedMesh2, so cdt and mesh both depend
                              │     downward on core and neither on the other)
                              │
                              ├──→ refinement   (raster for sampling, parallel_util)
                              │
                              └──→ flip         (predicates, parallel_util)
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
`refinement` gives it a caller, and `RasterView<T>` — a non-owning view over a
contiguous caller-supplied buffer. `RasterView` belongs in this module
(`view.hpp`), not in `bindings/`: it is pure C++, and every other zero-copy
source (an mmap'd tile, an HDF5 window, a sub-window of a parent raster)
produces the same type. Only the lifetime anchor sits with the bindings.

When `RasterView` lands it must bring contiguous row access into the
`RasterSource` concept in the same change, not afterwards. A scalar-only
concept forces a row-major scan to recompute `linear_index` per sample and
never to walk a row pointer, which forfeits the entire reason the view exists;
and adding a concept requirement later forces a revisit of every model and
every test double.

**Decided (@architect): GeoTIFF decoding lives in Python**, under
`src_python/tin_engine/io/`. The C++ core never opens a file, never sees a
path, and never links a codec.

The deciding argument is dependency gravity, not testability. GeoTIFF is not
an array format: it is a container plus a GeoKey directory plus a CRS.
Decoding it in C++ pulls CRS interpretation across the firewall, and CRS
interpretation means PROJ — GDAL's own dependency. That is the neighbourhood
this migration exists to leave. Keeping decode in Python also preserves
Tier-1 C++ tests that need no fixtures at all.

`README.md` and `testing.md` already assumed this; only the prose in this
section dissented.

**Boundary contract.** Exactly one Python module, `tin_engine/raster.py`,
constructs a core raster, so dtype, contiguity, writeability and
projected-CRS checks have one place to audit. `io/` produces pure Pydantic
data and never touches `_core`. Across the boundary go one C-contiguous 2-D
`float32` or `float64` array, four keyword-named affine scalars, and an
optional NoData sentinel — nothing else.

- Raster dimensions are derived from `array.shape`, never passed alongside it,
  so a shape/geometry disagreement is unrepresentable rather than validated.
  The legacy passed `(array, x_min, y_max, delta_x, delta_y)` positionally
  (`legacy/rasputin/reader.py:340-345`), which is the shape that let the
  transposed row/column defect live.
- numpy owns the buffer. The bound class holds a `py::object` reference as the
  primary lifetime anchor, with `keep_alive<1,2>` backing it; `.noconvert()`
  forbids a conversion copy that would leave the view pointing at a temporary.
- Integer DEMs are promoted in Python at decode time: 16-bit to `float32`,
  32-bit to `float64`. `int32` does not fit float32's mantissa, and promoting
  it there would silently quantise elevations.
- **CRS never crosses into C++**, now or later. The legacy violated this by
  pushing a proj4 string into the mesh; do not reintroduce it. Python
  reprojects every input into one projected CRS first, and rejects a
  geographic CRS outright — a degrees-based raster interpolates perfectly
  happily and yields a silently distorted mesh.
- NoData discovery is Python's (a container tag); NoData semantics are C++'s
  (already implemented in `raster.hpp` and `sample.hpp`). The sentinel is
  compared with `==`, so it must be passed exactly as decoded, never re-typed
  or rounded.
- Every C++ scan releases the GIL and the adapter sets the array read-only
  first. Refinement runs parallel-for with every thread sampling the same
  raster, so concurrent const access is a requirement, not an aspiration — a
  concept cannot express it, so it is stated here and in `raster.hpp`.

**Obligations this puts on the reader.** `geometry.hpp` declares two
load-bearing conventions — north-up, and grid-registered (pixel-is-point) —
which C++ can no longer verify, so the reader must enforce them. All three of
these are gaps in the legacy, verified:

- `GTRasterTypeGeoKey` (1025) is defined at `legacy/rasputin/reader.py:36` and
  read nowhere. An area-registered TIFF therefore lands half a cell off.
- `ModelTransformationTag` (34264) appears nowhere in the legacy reader, so a
  rotated or sheared transform is silently misread instead of rejected.
  `RasterGeometry` cannot represent one; it must raise.
- Missing georeferencing defaults rather than failing — tie point to zeros,
  pixel scale to `(1.0, 1.0)` — so an un-georeferenced TIFF silently becomes a
  unit-spaced raster at the origin.

**Two unrelated meanings of "window"**, which must not be merged: Python
windows are an I/O concern (which strips or tiles to decode); `window_for` is
an index concern (which cells a triangle's bbox covers). Same word, different
layers, no shared code.

### `predicates`

Exact `orient2d` and `incircle` behind `GeometryKernel`, as `FilteredKernel<E>`: an
inline static filter with an exact fallback. `DefaultKernel` binds it to vendored
detria at a pinned SHA. **Not header-only** — `src/predicates/detria_exact.cpp` is the
only TU that includes `detria.hpp`, enforced by CMake privacy, an `#error` guard and
`tools/check_detria_boundary.py`. Used by `core`, `noding`, `cdt` and `flip`. See
`docs/increments/01-predicates.md`.

### `parallel_util`

Wrappers over `<thread>`, `<atomic>`, `std::execution` (or a small work-stealing pool). Compaction primitives. Parallel sort/unique wrappers. No external runtime dependency by default; TBB or OpenMP can be opted in via CMake flag.

### `vector_simplify`

Visvalingam-Whyatt for area-preserving polygon simplification; Douglas-Peucker for polyline length-preserving; multi-feature topology checks (no introduced crossings). Pre-noding step.

### `hydrology`

Implements `auto_catchments.md`: pit filling (priority-flood with epsilon, plus optional Lindsay hybrid), D8 flow direction, parallel flow accumulation (Barnes-Lehman-Mulla), catchment BFS, stream extraction with Strahler ordering, mask polygonization. May depend on RichDEM (MIT) as either reference or external library — decision deferred to first implementation pass.

### `noding`

Implements the constraint-noding section of `parallel_refinement.md`. Uniform-grid broad phase (the raster's grid is convenient, not required — it is a spatial index, not the snap grid), robust pairwise intersection, snap rounding, segment splitting, deduplication. Outputs a clean PSLG. `is_river` is one bit per **chain**, not per edge — only the
noder can produce an edge whose bit disagrees with its source chain, and it
contributes a sparse per-edge override set then. See `docs/increments/03-pslg.md`.

### `cdt`

Thin wrapper around vendored Detria (`lib/detria/`, pinned SHA, MIT). Translates
between rasputin's `Pslg` and the library's API, and returns `IndexedMesh2` from
`core/`. Used exactly once per mesh build.

**poly2tri is not a fallback.** It has no per-edge constraint entry point — only an
outer polyline plus holes plus Steiner points — so breaklines cannot be expressed,
which is to say rivers cannot be. A Steiner point is not a constraint: the mesh may
triangulate around it however it likes, so a river would become vertices the mesh
happens to contain rather than edges it is required to contain. It also inserts
points, which breaks the per-triangle constrained-edge mask. If Detria is ever
dropped the honest options are a different library or our own CDT over the
`predicates` module, not poly2tri. See `docs/increments/04-cdt.md`.

### `mesh`

Core data structures: vertex array, triangle array, ternary tree of refinement nodes, per-edge constraint bitmask (currently just `is_river`). Owns the flatten-to-final-mesh step. Used by `refinement`, `flip`, and bindings.

### `refinement`

Implements the refinement loop in `parallel_refinement.md`. Active-list management, per-triangle DEM-sample scan with max-error tracking, fan subdivision at the max-error point, round-by-round compaction.

### `flip`

Final Lawson edge-flip pass. Skips constraint-tagged edges. Parallel with edge-conflict resolution.

### `bindings/core.cpp`

Single pybind11 module that exposes the C++ API to Python. Built as the `_core` extension, installed into the `tin_engine` package. The **value types only** are re-exported by `src_python/tin_engine/__init__.py` (`__all__` is `Point2`, `Point3`, `cross`, `dot`); the CDT surface below is reached as `tin_engine._core` and is deliberately not re-exported, because `cli.py` is the sole composition root and `viz/` never imports the extension. Typed from `src_python/tin_engine/_core.pyi`, which is what `mypy --strict` sees.

As of increment 6a (`94f94e2`, `docs/increments/06-cdt-viewer.md`) it binds:

- the `Point2`/`Point3` value types, whose `__repr__` routes through the `std::formatter` specializations in `point.hpp` so the C++ and Python renderings cannot drift;
- the `dot` and `cross` free functions over both point types;
- `Pslg`, `Chain`, `PslgDiagnostic`, `PslgBuildResult`, `IndexedMesh2` and `CdtOutcome`, plus the `ChainRole`, `PslgError` and `CdtStatus` enums and a `describe(CdtStatus)` helper;
- two more free functions: `build_pslg`, which validates and returns a `PslgBuildResult` carrying a diagnostics list rather than raising, and `triangulate`, which wraps the kernel call in `py::gil_scoped_release` -- the only call in the module long enough to be worth the release.

Every array-shaped accessor -- `Pslg.vertices`, `Pslg.chain_indices`, `Pslg.indices_of`, `IndexedMesh2.vertices`, `.triangles`, `.constrained_edges` -- returns a **read-only, zero-copy `py::array_t` whose base object is the owner**, built through the single `readonly_view` helper. Nothing in these six is copied and none outlives its buffer. Scoped deliberately: `build_pslg` *does* copy the coordinates handed to it, and `Pslg.chains` and `PslgBuildResult.diagnostics` rebuild per attribute read.

**`PslgBuilder` is deliberately not exposed, and neither is any kernel template parameter.** The decisive reason is that `PslgBuildResult build() &&` is rvalue-ref-qualified: binding it means either a consuming method on a Python object that remains reachable afterwards, or a lambda that copies. Beyond that, the builder is a mutable accumulator; binding it would put half-built, validation-pending state on the Python side and make `Pslg`'s "validated on construction" guarantee unenforceable from there. Python hands `build_pslg` a finished array of vertices and chains and gets back either a `Pslg` or the reasons it is not one. A new accessor on this surface is a design change that needs a reason in an increment file, not a line in a PR.

## Build system

CMake-based, building the header-only core plus one Python extension. C++20 required. Highlights:

- **Top-level `CMakeLists.txt`** orchestrates the build. `terrain_headers` is an `INTERFACE` target carrying `include/`. Two options gate the rest: `RASPUTIN_BUILD_PYTHON` (default `OFF`) adds the pybind11 extension, `RASPUTIN_BUILD_TESTS` (default `ON`) adds the C++ suite.
- **`tests/cpp/CMakeLists.txt`** fetches Catch2 v3 via `FetchContent`, so it needs no manual checkout.
- **Detria** (header-only) is vendored at a pinned SHA in `lib/detria/`. CMake does
  **not** search for a system copy: version skew in a geometry kernel across machines
  is a reproducibility hazard and vendoring a header costs nothing.
- **External deps under consideration:** RichDEM (optional, MIT), Eigen (if linear algebra needs grow beyond what we want to hand-roll).
- **No CGAL, no Boost.Geometry, no GDAL** in the new core. Existing Python-layer uses are migrated incrementally.
- **Planned:** per-module `OBJECT` libraries linked into the extension once `src/` is populated, and sanitizer flags via `RASPUTIN_SANITIZER=asan|ubsan|tsan|none`.

Packaging is driven by **scikit-build-core**, declared in `[build-system]` in `pyproject.toml`, which invokes this same CMake build. `pip install .` configures with `RASPUTIN_BUILD_PYTHON=ON` and `RASPUTIN_BUILD_TESTS=OFF` — the C++ tests pull Catch2 over the network and have no business running during an install — and installs `_core` into the `tin_engine` package.

There is no `setup.py`. It was removed with the foundation reset, along with the `CMakeBuild` class it carried.

## Python API surface

The public API is the `tin_engine` package, calling into `tin_engine._core`. `tin_engine.viz` is the renderer that turns a triangulation into an SVG a person can look at; it consumes the `typing.Protocol`s in `viz/protocols.py` and **never imports `_core`**, so it is testable with no compiled extension in the process. `cli.py` is the single composition root that joins the two -- the same shape as the rule below that exactly one module adapts decoded raster data into `_core`. The `rasputin draw` command drives it: it looks a fixture up in `viz.fixtures.GALLERY`, maps that fixture's **string** chain roles onto `_core.ChainRole` -- `viz/` may not name the enum, so the mapping is the composition root's -- runs `build_pslg` and `triangulate`, and hands the fixture itself to `build_scene` as the `PslgLike`, which is what lets a fixture the validator rejects still be drawn. It is also the only place a path exists, and it resolves and refuses one before writing.

The pre-migration `rasputin.*` modules (`mesh.py`, `geometry.py`, `reader.py`, `tin_repository.py` and friends) are archived under `legacy/rasputin/` rather than kept in place, so this is a re-implementation against the new backend rather than a rewiring of stable modules. Porting proceeds one entry point at a time; shapes worth preserving should be read out of `legacy/` before being reintroduced.

## What gets deleted, eventually

After the new backend ships and the Python API is rewired:

- `legacy/rasputin/triangulate_dem.h` and `legacy/bindings.cpp` (the CGAL-based originals)
- CGAL, GMP, MPFR from the CMake dependency list — already absent from the new `CMakeLists.txt`
- Boost.Geometry — prohibited in the new core (`CLAUDE.md` §2) and machine-checked;
  the remaining uses are in `legacy/`, which is exempt

These removals are not part of the initial build-out; they're a follow-up once feature parity is reached and tests pass on the new backend.

## Existing files to integrate

The repo already has:

- `legacy/rasputin/triangulate_dem.h` (CGAL-based, to be replaced)
- `legacy/rasputin/*.py` (the pre-migration Python layer, pending per-module classification)
- `lib/date/` — **removed**. The C++20 `<chrono>` calendar types replaced it; `CLAUDE.md` section 2 prohibits external `date` libraries.

The new `_core/` tree is added alongside the existing C++ files; the old files stay buildable until the new pipeline reaches parity, then are removed in a single cleanup commit.
