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
    edge_properties.hpp    # EdgeProperties: a set of 32 opaque feature bits,
                           #   no feature NAME anywhere in terrain::
    pslg.hpp               # Pslg, Chain, ChainRole — validated planar input
    pslg_builder.hpp       # PslgBuilder, PslgDiagnostic, the validator
    snap_grid.hpp          # GridPoint, SnapGrid: the lattice noded output lies on
    noded_pslg.hpp         # NodedPslg: the noder's output, the CDT's input (5b)
    indexed_mesh.hpp       # IndexedMesh2: SoA mesh + per-triangle constrained mask
  predicates/
    orientation.hpp        # Orientation / Incircle vocabulary; no includes at all
    exact.hpp              # ExactPredicates concept, filter bounds
    kernel.hpp             # GeometryKernel, FastKernel, FilteredKernel<E>
    detria_exact.hpp       # DetriaExact declaration; does NOT include detria
    default_kernel.hpp     # DefaultKernel = FilteredKernel<DetriaExact>
  noding/
    intersect.hpp          # SegmentRelation, classify, crossing_point,
                           #   segment_meets_cell (the hot-pixel predicate)
    node_set.hpp           # NodeSet: dedup by sorted grid keys
    broad_phase.hpp        # uniform bucket index, for_each_candidate (5b)
    noded_pslg_builder.hpp # NodeStatus, NodeOutcome, the verifier (5b)
    node.hpp               # NodeOptions and node<K>, the driver (5b)
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
    window.hpp             # window_for: bbox -> index window (planned; unbuilt,
                           #   refinement walks lattice nodes directly)
  parallel_util/
    chunks.hpp             # for_each_chunk: contiguous chunks over std::jthread
  mesh/
    lattice_mesh.hpp       # LatticeMesh: flat triangle array over DEM nodes,
                           #   neighbour links, the three splits (14), flip (14b)
    lawson.hpp             # LatticeFrame, legalise_around / legalise_all:
                           #   Lawson flips on strictly-inside apexes (14b)
  refinement/
    scan.hpp               # per-triangle sup-norm scan, NoData carve point (14)
    refine.hpp             # RefineOptions, RefineOutcome, the round loop (14),
                           #   Delaunay insertion (14b)

src/                       # C++ implementation, one directory per module
                           #   (only predicates/ and cdt/ exist; rest planned)
  predicates/              # exact orient2d/incircle; namespace terrain::pred
  parallel_util/           # (none: header-only, include/terrain/parallel_util/)
  vector_simplify/         # Visvalingam-Whyatt, Douglas-Peucker, topology checks
  hydrology/               # pit fill, flow direction, accumulation, catchments, streams
                           # (no noding/ here, and none planned: the noder is
                           #  header-only, its driver a template on the kernel)
  cdt/                     # thin wrapper over vendored Detria
  mesh/                    # (none: header-only, include/terrain/mesh/)
  refinement/              # (none: header-only, include/terrain/refinement/)
  flip/                    # final Lawson edge-flip pass, constraint-respecting

bindings/
  core.cpp                 # pybind11 module definition -> tin_engine._core

src_python/tin_engine/     # public Python API (distribution name: rasputin)
  __init__.py              # re-exports from tin_engine._core
  cli.py                   # Typer entry point declared in pyproject
  raster.py                # the ONLY adapter from decoded data into _core
  grid_domain.py           # DEM extent -> stride-subsampled nodes + outer ring;
                           #   pure numpy, never imports _core
  elevation.py             # drops mesh vertices the DEM has no data for;
                           #   pure numpy, never imports _core
  features.py              # EdgeVocabulary: which bit means which feature.
                           #   The names C++ refuses to hold. Imports nothing
                           #   first-party and never imports _core
  _core.pyi                # type stubs for the compiled extension
  viz/                     # CDT -> SVG renderer; never imports _core
    __init__.py            # re-exports Scene, SvgStyle, build_scene, render_svg
    protocols.py           # MeshLike / PslgLike / ChainLike -- typing.Protocol
    scene.py               # build_scene(pslg, mesh=None, ...) -> Scene
    style.py               # SvgStyle: frozen Pydantic V2. Canvas geometry,
                           #   plus PropertyStroke — the bit -> CSS-token draw
                           #   precedence, which is structure, not taste, and
                           #   carries positions rather than feature names
    svg.py                 # (Scene, SvgStyle) -> str; the stylesheet lives here
    fixtures.py            # the synthetic gallery, declarative; `rasputin draw
  io/                      # all file decoding AND encoding lives here
    __init__.py
    ply.py                 # arrays -> PLY bytes; takes no path and opens nothing
    vtk_legacy.py          # arrays + EdgeVocabulary -> legacy .vtk bytes, for
                           #   ParaView; takes no path and opens nothing
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

                            mesh   (flat array + adjacency; edge tags — FROM an
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

`RasterView<T>` (`view.hpp`, increment 12) is a non-owning view over a
contiguous caller-supplied buffer. It belongs in this module, not in
`bindings/`: it is pure C++, and every other zero-copy source (an mmap'd tile,
an HDF5 window, a sub-window of a parent raster) produces the same type. Only
the lifetime anchor sits with the bindings. `window_for` (bbox to index
window) stays unbuilt: `refinement` works in lattice coordinates, where a
triangle's bounding box already is its index window (increment 14, R2).

`RasterView` brought contiguous row access into the `RasterSource` concept in
the same change, not afterwards. A scalar-only
concept forces a row-major scan to recompute `linear_index` per sample and
never to walk a row pointer, which forfeits the entire reason the view exists;
and adding a concept requirement later forces a revisit of every model and
every test double.

**Decided (@architect): GeoTIFF decoding lives in Python**, under
`src_python/tin_engine/io/`. The C++ core never opens a file, never sees a
path, and never links a codec.

The deciding argument is dependency gravity, not testability. GeoTIFF is not
an array format: it is a container plus a GeoKey directory plus a CRS.
Decoding it in C++ means a TIFF container library, a codec set and an EPSG
database in the core — three dependencies for a module whose entire job is to
index into an array. Python already has all three and a Pydantic layer to
validate them in. Keeping decode in Python also preserves Tier-1 C++ tests
that need no fixtures at all.

This paragraph used to justify the rule with "CRS interpretation means PROJ —
GDAL's own dependency". That was false and is corrected here: PROJ is a
standalone library, GDAL depends on PROJ rather than the reverse, and
`CLAUDE.md` §2 lists PyProj in the core stack, so this project depends on PROJ
and always did. The legacy used `pyproj` for every transformation and never
linked GDAL. `docs/increments/11-raster-ingestion.md` ruling 12 is the record.

`README.md` and `testing.md` already assumed this; only the prose in this
section dissented.

**Boundary contract.** Exactly one Python module, `tin_engine/raster.py`,
constructs a core raster, so dtype, contiguity, writeability and
projected-CRS checks have one place to audit. `io/` produces pure Pydantic
data and never touches `_core`. Across the boundary go one C-contiguous 2-D
`float32` or `float64` array, four keyword-named affine scalars, and an
optional NoData sentinel — nothing else.

- Across the C++ boundary, raster dimensions are derived from `array.shape`,
  never passed alongside it, so on the C++ side a shape/geometry disagreement
  is unrepresentable rather than validated. The Python side is different:
  `io/models.py`'s `RasterMeta` carries `rows` and `cols` (increment 11,
  choice B1), and `DemTile` *validates* them against `array.shape`.
  `raster.py` must not forward them.
  The legacy passed `(array, x_min, y_max, delta_x, delta_y)` positionally
  (`legacy/rasputin/reader.py:340-345`), which is the shape that let the
  transposed row/column defect live.
- numpy owns the buffer. The bound class holds a `py::object` reference as the
  primary lifetime anchor, with `keep_alive<1,2>` backing it; `.noconvert()`
  forbids a conversion copy that would leave the view pointing at a temporary.
- Integer DEMs are promoted in Python at decode time: 16-bit to `float32`,
  32-bit to `float64`. `int32` does not fit float32's mantissa, and promoting
  it there would silently quantise elevations.
- **No CRS string crosses into C++**, now or later, because nothing there
  reads one. No header in `raster/` has a CRS field and none can use it; a
  field nothing reads is a field that rots. The legacy proved it by pushing a
  proj4 string into the mesh for no consumer. Do not reintroduce it. This is
  an argument about dependency surface and unused state, and nothing more.
- **Every number that crosses is metres in a projected CRS**, because nothing
  on either side can tell otherwise. This is the load-bearing rule, and it is
  separate from the one above: that one is about where a string may live, this
  one is about what the numbers mean. Python must deliver every input in one
  projected, metre CRS, and rejects a geographic CRS outright. Today it does
  this by refusal alone: `io/geotiff.py` refuses any file that is not already
  in one, and nothing reprojects yet. Reprojection arrives no earlier than
  the mosaic increment (`docs/increments/11-raster-ingestion.md` §9, §10).
  Measured, pyproj 3.8.0 / PROJ 9.8.1: a 0.0002777° cell at 60°N is 15.5 m
  east-west and 31.0 m north-south, a 2:1 anisotropy invisible to `sample.hpp`, whose bilinear
  weights would then be computed in degrees and applied to metres. The legacy
  never implemented this rejection — a geographic GeoTIFF happened to raise
  during proj4 assembly, so the defence was a typo. See
  `docs/increments/11-raster-ingestion.md` sections 3 and 5.
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

- `GTRasterTypeGeoKey` (1025) is defined in `legacy/rasputin/reader.py` and
  read nowhere. An area-registered TIFF therefore lands half a cell off. The
  reader converts such a file rather than refusing it — the node grid is the
  declared corner grid shifted inward by half a cell — and refuses only a file
  that does not say which convention it uses. The repository's own DEM fixture
  is area-registered, which is why refusal was the wrong rule;
  `docs/increments/11-raster-ingestion.md` ruling 4 carries the measurement.
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

Header-only. Today one helper, `for_each_chunk(n, threads, fn)` in `chunks.hpp`: contiguous chunks over `std::jthread`, created per call and joined before it returns, no pool. `std::execution::par` and OpenMP were both ruled out in increment 14 (R7): neither builds on macOS without an experimental flag or an extra runtime. Needs only `Threads::Threads`.

### `vector_simplify`

Visvalingam-Whyatt for area-preserving polygon simplification; Douglas-Peucker for polyline length-preserving; multi-feature topology checks (no introduced crossings). Pre-noding step.

### `hydrology`

Implements `auto_catchments.md`: pit filling (priority-flood with epsilon, plus optional Lindsay hybrid), D8 flow direction, parallel flow accumulation (Barnes-Lehman-Mulla), catchment BFS, stream extraction with Strahler ordering, mask polygonization. May depend on RichDEM (MIT) as either reference or external library — decision deferred to first implementation pass.

### `noding`

Implements the constraint-noding section of `parallel_refinement.md`. Uniform-grid broad phase, robust pairwise intersection, snap rounding, segment splitting, deduplication. Outputs a `NodedPslg` — the deduplicated node array, the chains with their roles and windings preserved, the per-edge property array and the snap grid.

**The broad phase is a spatial index and its cell size has no relationship to the snap grid's.** The raster's grid is convenient, not required; the snap grid is not even an option. Here is the number that makes it non-negotiable: at a 5 cm spacing a 100 km domain has a snap lattice of 2e6 × 2e6 = **4e12 cells**, and at a decimetre 1e12. A broad phase bucketed at that spacing would allocate one bucket per cell to hold a few thousand segments. The two grids answer different questions — the snap grid quantises *coordinates* and is sized by input precision; the broad phase filters *candidate pairs* and is sized by segment density.

Feature semantics are **per chain** on `Pslg` and a dense array **per edge** on `NodedPslg`, index-aligned with that type's flat edge enumeration — not a sparse override set, because after noding an edge can descend from two chains at once and so has no single source value to override. The merge is **union**, which is commutative, associative and idempotent and therefore reducible over an unordered set of contributing chains: an edge can be both a road and a river. The one-bit `is_river` form was replaced by a property set in increment 7 (`core/edge_properties.hpp`, up to 32 opaque bits); see `docs/increments/07-edge-properties.md` for the ruling, `docs/increments/03-pslg.md` for the record of the form it replaced, and `docs/increments/05b-noder-driver.md`, "Edge properties", for the shape on `NodedPslg`.

**Header-only: there is no `src/noding/` and none is planned.** The driver is a template on the kernel, as increment 3's validator is.

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

`LatticeMesh` (`lattice_mesh.hpp`, increment 14): a flat triangle array over `(row, col)` DEM-node vertices with, per triangle edge, the neighbour across it, a constrained bit and a 32-bit property mask (increment 7's `EdgeProperties` bits — never one `is_river` flag). **Not a ternary tree:** an edge split changes two triangles at once, which a tree without neighbours cannot represent. A split reuses the parent's slot and appends the rest, so the array is the output as it stands and there is no flatten step. Depends on `core` only; it knows no raster. The adjacency is what a later `flip` pass needs. See `docs/increments/14-adaptive-refinement.md`, R4.

### `refinement`

Refinement against the DEM to a sup-norm tolerance (increment 14). `scan.hpp` measures a triangle's largest `|z - plane|` over the DEM nodes it contains, by exact integer tests; `refine.hpp` runs rounds of parallel scan (`parallel_util`) and serial splits in index order — a fan for an interior node, an edge split on both sides for an edge node — so the output does not depend on the thread count. Triangles with a NoData vertex are carved, not refined. See `docs/increments/14-adaptive-refinement.md`.

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
- **No CGAL and no GDAL** in the new core: prohibited by `CLAUDE.md` §2.
  **No Boost.Geometry** either, which is a scope choice and not a prohibition
  — reading this line as one is what put an unauthored ban in §2 for six
  days. Existing Python-layer uses are migrated incrementally.
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
- Boost.Geometry, whose only uses are in `legacy/triangulate_dem.h` and go when
  it does. Not a prohibited dependency — `CLAUDE.md` §2 is the list, and this
  is not on it

These removals are not part of the initial build-out; they're a follow-up once feature parity is reached and tests pass on the new backend.

## Existing files to integrate

The repo already has:

- `legacy/rasputin/triangulate_dem.h` (CGAL-based, to be replaced)
- `legacy/rasputin/*.py` (the pre-migration Python layer, pending per-module classification)
- `lib/date/` — **removed**. The C++20 `<chrono>` calendar types replaced it; `CLAUDE.md` section 2 prohibits external `date` libraries.

The new `_core/` tree is added alongside the existing C++ files; the old files stay buildable until the new pipeline reaches parity, then are removed in a single cleanup commit.
