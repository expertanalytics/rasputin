# Increment 12 — from a GeoTIFF on disk to an elevated mesh in ParaView

Status: **implemented, in review.** Red `39dba00` (with test fixes `2696196`,
`5a36627`), green `3e2e9ec` (C++) and `af8e199` (Python). Measured: 319 net
production lines (C++ 114, Python 264 added, 59 removed) against ~320.
Unruled points settled in green: `--stride` with a fixture is refused;
`subsample` raises `ValueError` for a stride below 1; `trim`'s `dropped`
counts only vertices without data, since that is what the file says, and a
valid vertex left with no triangle is removed but not counted.

*Amended before merge, at the user's request (2026-09-25).* In ParaView the
Color By menu offered `elevation`, but it was the dataset string describing
the sampling; the heights were only the points' z. The `.vtk` now also writes
`POINT_DATA` with `SCALARS elevation double` equal to z, and the string field
is renamed `elevation_source`. `elevation` is a reserved field name in
`write_vtk`, so no string can shadow the heights again. Red `5750ff8`.

Written by
`@architect` before `@tester`,
per `docs/increments/README.md` step 1. The user chose the recommendation on
all three open choices on 2026-09-24: U1 (a), U2 (a), U3 (a) (section "Ruled by
the user").

**Closes.** `ROADMAP.md` MVP gap 1 (the increment 12 half: the adapter), gap 3
(z per mesh vertex) and the first half of gap 5 (a CLI path from a file on disk
to a mesh on disk). After this increment,

```sh
rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --out tile.vtk
```

writes one `.vtk` file with real elevation that ParaView opens.

**Not closed.** Gap 2 (refinement). The vertices here are a regular subsample
of the DEM's own nodes. That is a stand-in, chosen to be trivially simple and
to be thrown away when refinement lands. Catchment clipping is not here either.

**One increment, not two.** The estimate below is about 320 lines, about 450 on
increment 10's overrun bias. That is under `CLAUDE.md` §2's 700. Splitting the
adapter from the sampling and the CLI would mean the first increment ends with
nothing the user can open, which is what the user asked not to have.

## What was measured

All on `9af8225`, with the installed `_core`. Commands are `uv run python -c`
one-liners over the objects named; `uv run --with imagecodecs` where the real
fixture is decoded.

1. **The CDT keeps interior vertices no chain references.** A square outer ring
   plus three free interior points: `node` keeps all 7, `triangulate` returns 7
   vertices and 8 triangles. A 4 × 4 node grid with the 12 perimeter nodes as
   the outer ring gives 16 vertices and 18 triangles, which is 2 × 3 × 3. So a
   grid subsample needs one outer chain and nothing else.
2. **The real fixture's geometry, after decode:** `x_min=799750.0`,
   `y_max=7950250.0`, `delta 10.0`, `5051 × 5051` nodes, `epsg=25833`,
   `nodata=-32767.0` (from the tag), `pixel_is_area=True` (already shifted in by
   half a cell). Valid z from −1.24 m to 391.7 m.
3. **The fixture has NoData on its outline.** Row 0 is 100 % NoData and column
   0 is 57 % NoData; 0.17 % of all nodes are NoData. At stride 50, 159 of the
   10 404 sampled nodes are NoData. A design that needs z on every outline
   vertex fails on the one real file this project has.
4. **Snapping does not move lattice-aligned grid nodes at UTM scale.** A stride
   50 subsample of the fixture's grid, noded at the CLI's 1e-3 spacing: every
   mesh vertex coordinate is bit-equal to a grid node coordinate, and none lies
   outside `[x_min, x_max] × [y_min, y_max]`.
5. **Cost.** Stride 26 (196 × 196 nodes): 76 050 triangles, 0.18 s for
   `build_pslg` + `node` + `triangulate`. Stride 10 (506 × 506): 510 050
   triangles, 3.0 s.
6. **Without `imagecodecs`, `decode_dem` refuses the fixture** with a
   `GeoTiffError` that names Compression (259) = 5 and points at the `codecs`
   extra (`src_python/tin_engine/io/geotiff.py:169`). Nothing new is needed
   for the no-extra case beyond turning that error into a usage error.

## Rulings

### R1. Vertex selection: a regular stride subsample of the DEM's own nodes

- Take every `stride`-th node in both directions, starting at node `(0, 0)`.
  Always include the last row and the last column, so the mesh covers the
  whole tile extent.
- The outer chain is the perimeter of that subsample, counter-clockwise, with
  property mask 0. Every other subsample node is a free vertex. No breaklines,
  no holes.
- Triangulate with the existing `build_pslg` → `node` → `triangulate` path,
  Delaunay on. Not a hand-rolled two-triangles-per-quad split: the CDT is what
  refinement, clipping and breaklines will go through, and using it now means
  the DEM path and the fixture path hand `write_vtk` the same kind of mesh.
- `--stride N` is optional. The default is the smallest stride that gives at
  most 256 nodes on the longer side:
  `max(1, ceil((max(rows, cols) - 1) / 255))`. For the fixture that is 20, so
  254 × 254 nodes and about 128 000 triangles, well under a second by
  measurement 5.
- Vertex coordinates are computed as `x_min + col * delta_x` and
  `y_max - row * delta_y` in float64, the same expression as
  `RasterGeometry::node` (`include/terrain/raster/geometry.hpp:66`). Any other
  spelling can land one ulp off the node and change which cell `bilinear`
  picks.
- Boundary-only is ruled out: four corners and a flat interior is not terrain.
  Every node is ruled out: 25.5 M points for the fixture.

This lives in a new pure module, `tin_engine/grid_domain.py`. It takes a
`RasterMeta` and a stride and returns the `(N, 2)` vertices and the ring's
indices. It imports neither `_core` nor `io.geotiff`.

### R2. Sampling happens in C++, through the existing `bilinear`, unchanged

- z is sampled at the **mesh's** vertex coordinates, after noding and
  triangulation, because those are the vertices written.
- The sampler is `terrain::raster::bilinear` (`include/terrain/raster/sample.hpp`)
  over the new `RasterView<T>`. Its semantics do not change in this increment.
- A new batch function in `sample.hpp` runs it over a span of points and writes
  `z` and a `valid` flag per point. It lives in the core, not in `bindings/`,
  so C++ can test it with no Python. The binding releases the GIL around it.
- What crosses back is two arrays: `z` float64 `(N,)` and `valid` bool `(N,)`.
  No NaN is used as a "no value" marker, so none can leak into a file.

What `bilinear` does at the edges, read from the code and pinned by the tests
in `tests/cpp/unit/test_raster.cpp`:

- **Outside the grid:** `nullopt`. `cell_of` refuses anything outside
  `[x_min, x_max] × [y_min, y_max]` and any non-finite coordinate.
- **On the last row or column:** in range. `bilinear_cell_of` clamps the cell
  to `rows - 2` / `cols - 2`, so a node on the far edge is the far corner of
  the last cell. There z is the node value up to rounding in `tx`, `ty`
  (about 1e-10 relative), not bit-exact. Tests must use a tolerance there.
- **NoData:** `nullopt` if **any** of the four corners is NoData or NaN, even
  when that corner's weight is zero. So a valid node next to a NoData node also
  gets no z. On a 10 m grid that trims one cell (10 m) around every void. This
  is accepted for now and recorded here; changing the sampler is not this
  increment's job.
- **Area/point registration:** handled before the sampler. `decode_dem` has
  already shifted an area-registered file's corner in by half a cell
  (increment 11, ruling 4), so `RasterMeta` always describes the node grid, and
  R1 places vertices on that grid.

### R3. Vertices with no z are removed, with their triangles and edges

The fixture's outline is NoData (measurement 3), so "refuse the file" would
refuse the only real DEM, and writing NaN into `POINTS` breaks ParaView's
bounds. So, after sampling:

- drop every triangle with a vertex whose `valid` is false;
- drop every constraint edge with such an endpoint;
- drop vertices no remaining triangle uses, and renumber.

The result has a ragged edge or a hole where the DEM has no data. This is a
pure numpy function over arrays, in `tin_engine/elevation.py`, with no `_core`
import. If nothing is left, the command exits non-zero and writes no file. The
count of dropped vertices goes into the file (R5) and to stderr.

Known limit, stated so no one rediscovers it: a DEM whose corner is not on the
1e-3 snap lattice can have boundary vertices snapped outward by up to 0.5 mm.
Those vertices then fall outside the grid and are dropped by this rule, and
the drop count shows it. Measurement 4 shows the real fixture is not affected.
Refinement replaces this code path; do not build a workaround for it here.

### R4. The C++ side: the view, the concept change, the binding

This is increment 12 as `project_structure.md` planned it (the `raster`
section), without changes:

- `include/terrain/raster/view.hpp`: `RasterView<T>`, a non-owning view over a
  contiguous row-major buffer. Constructor takes a `RasterGeometry`, a
  `const T*` and an `optional<T>` NoData. `value_at`, `is_nodata`, `geometry`,
  `row`. Pure C++, no pybind11.
- `raster.hpp`: `RasterSource` gains row access,
  `{ r.row(i) } -> std::same_as<std::span<const typename R::value_type>>`.
  `Raster<T>` and the test double `ConstantRaster` gain `row`. This lands in
  the same change as the view, as `project_structure.md` already rules. There
  is no row-walking caller yet; the user chose to add it now (U2 (a)).
- `bindings/core.cpp`: one bound class, `RasterView`, holding a
  `std::variant<RasterView<float>, RasterView<double>>`, built by
  `raster_view(array, *, x_min, y_max, delta_x, delta_y, nodata=None)`.
  - `array` is `py::array_t<T, c_style>` with `.noconvert()`, tried as float32
    then float64. A non-contiguous or other-dtype array is a `TypeError`, never
    a silent copy.
  - Rows and columns come from `array.shape`, never from an argument.
  - The bound object keeps a `py::object` reference to the array, backed by
    `keep_alive`.
  - The NoData value is converted to `T` exactly once. `decode_dem` already
    refuses a sentinel the cell type cannot hold
    (`src_python/tin_engine/io/geotiff.py:307`), so the conversion is exact.
  - `sample(view, points)` takes a float64 `(N, 2)` array and returns
    `(z, valid)`. It releases the GIL.
- `_core.pyi`: the class, the factory and `sample`.

No CRS string, EPSG code, path or file handle crosses. The four affine scalars
are keyword-only, which is what stops the legacy's positional transposition.

### R5. The Python side, and the CRS and units contract

- `tin_engine/raster.py` is the only module that builds a core raster:
  `to_core(tile: DemTile) -> RasterView`. It passes `tile.array` and
  `meta.x_min`, `meta.y_max`, `meta.delta_x`, `meta.delta_y`, `meta.nodata`.
  It does **not** pass `rows`, `cols`, `epsg`, `nodata_source`,
  `pixel_is_area` or `vertical_unit_assumed`. It checks the array is read-only
  (a `DemTile` guarantees this) and refuses otherwise.
- **Units.** Every number that crosses is metres in a projected CRS. That is
  already guaranteed upstream: `decode_dem` refuses any file whose CRS is not
  projected with metre axes (increment 11, §5). This increment adds no
  reprojection and no new check.
- **CRS in the output: yes.** The `.vtk` records `("crs", f"EPSG:{meta.epsg}")`
  as a dataset field, the same mechanism `--crs` uses today. With `--dem`,
  `--crs` is refused as a usage error, because the file already says what its
  CRS is and a second source could contradict it.
- The `elevation` field changes from `none (z=0, --flat)` to a short ASCII
  sentence, e.g. `bilinear from DEM, stride 20, 612 vertices without data
  dropped`. When `meta.vertical_unit_assumed` is true, the sentence also says
  `vertical unit assumed metres`. The file name is not recorded: `write_vtk`
  refuses non-ASCII strings, and a file name can be non-ASCII.
- `cli.py`'s module docstring says it is "the only Python module that imports
  `tin_engine._core`". That stops being true when `raster.py` lands, and
  `project_structure.md` has said since before increment 11 that `raster.py`
  constructs the core raster. The docstring is corrected in this increment to
  "the only module that has a path". `viz/` still never imports `_core`.

### R6. The CLI: `rasputin mesh --dem PATH`

- `mesh` gains `--dem PATH` and `--stride N`. The fixture name becomes
  optional. Exactly one of the fixture name and `--dem` must be given; both or
  neither is a usage error.
- With `--dem`, `--flat` and `--crs` are refused. `--out`, `--out-edges`,
  `--out-parent`, `--binary/--ascii`, `--delaunay/--no-delaunay` and
  `--snap-spacing` keep their meaning. `.ply` output also works, with the same
  real z, because the writers already take `(N, 3)`.
- The CLI opens the file (it is the module with paths) and passes the stream
  to `decode_dem`. A `GeoTiffError` becomes `typer.BadParameter` carrying the
  reader's own message, so a missing `imagecodecs` shows as exit 2 with the
  existing "the `codecs` extra (imagecodecs)" advice. No traceback.
- A missing or unreadable file is a usage error, not a traceback.
- The engine step is shared with the fixture path. `_triangulated` is split so
  that its core takes `(vertices, chains)` rather than a `Fixture`; the DEM
  path must not import `viz/`, as increment 13 ruling 6 already requires for
  mesh output.
- Order: decode → `grid_domain` → build/node/triangulate → `to_core` →
  `sample` → `elevation` trim → write. Each step is a function of the previous
  step's output. Only the first and the last touch the file system.

### The data flow

```
path ──(cli.py opens)──> stream ──decode_dem──> DemTile ─┬─ meta ──grid_domain──> (xy, ring)
                                                          │                           │
                                                          │          build_pslg/node/triangulate
                                                          │                           │
                                                          └─ to_core ──> RasterView   IndexedMesh2
                                                                             │         │
                                                                    _core.sample(view, mesh.vertices)
                                                                             │
                                                                        (z, valid)
                                                                             │
                                                        elevation.trim(mesh arrays, edges, z, valid)
                                                                             │
                                                     write_vtk / write_ply ──(cli.py writes)──> path
```

C++ sees one array, four scalars, one optional sentinel and one `(N, 2)`
point array. Python sees no C++ pointer.

## Ruled by the user

**The user chose U1 (a), U2 (a) and U3 (a) on 2026-09-24.** The options are
kept so the reasons stay on record.

**U1. What to do with vertices that have no elevation.**
- (a) Drop them with their triangles and edges; record the count. **Recommended.**
  It works on the real fixture and writes no NaN.
- (b) Refuse the run if any vertex has no z. Strictest, but it refuses the
  only real DEM in the repository (measurement 3).
- (c) Write the vertex with z = NaN and let the viewer cope. Keeps the full
  outline, but NaN in `POINTS` breaks ParaView's bounds and spreads into
  anything computed from the file. This project has refused silent NaN
  everywhere else (`include/terrain/raster/raster.hpp:62`).

**U2. Row access in the concept now, or when a row-walking caller exists.**
Two sound principles conflict. `project_structure.md` rules that row access
lands with the view, because adding a concept requirement later means revisiting
every model and test double. On the other side, nothing in this increment walks
a row, and a requirement with no caller is untested in use.
- (a) Add it now, as already ruled. About 10 production lines and one test
  double. **Recommended**: the cost is small and it keeps the written ruling
  true.
- (b) Defer it to refinement and amend `project_structure.md` to say so.

**U3. The CLI shape.**
- (a) `rasputin mesh --dem PATH`, fixture name optional. **Recommended.** One
  command that writes meshes, one set of output options.
- (b) A separate `rasputin dem PATH` command. Cleaner argument parsing, but the
  output options would be duplicated or factored out, which costs lines for no
  user gain.

## Prior art in `legacy/`

```sh
grep -rlni 'bilinear\|interpolate\|point_elevation\|get_interpolated' legacy/
```

returned `legacy/bindings.cpp` and `legacy/rasputin/triangulate_dem.h`.

- `legacy/rasputin/triangulate_dem.h:376` is `get_interpolated_value_at_point`,
  the legacy bilinear sampler. It was already ported in `7785fea` as
  `raster/sample.hpp`, fixing an out-of-bounds read and a row/column
  transposition (`project_structure.md`, `raster` section). Nothing more to
  carry.
- `mesh_from_raster` (`legacy/rasputin/triangulate_dem.h:542`) fed **every**
  raster node to CGAL with its own z, plus interpolated boundary points, and
  then simplified. That is refinement's job here (gap 2). The stride subsample
  of R1 replaces it until then; it is deliberately not a port.

No `@migration-expert` step is needed: the one piece with intent to carry (the
sampler) has already been carried.

## Files

```
include/terrain/raster/view.hpp          # new. RasterView<T>
include/terrain/raster/raster.hpp        # concept gains row(); Raster::row
include/terrain/raster/sample.hpp        # batch bilinear over a span of points
bindings/core.cpp                        # RasterView class, raster_view(), sample()
src_python/tin_engine/_core.pyi          # stubs for the three
src_python/tin_engine/raster.py          # new. to_core(DemTile) -> RasterView
src_python/tin_engine/grid_domain.py     # new. RasterMeta + stride -> (xy, ring). No _core
src_python/tin_engine/elevation.py       # new. drop vertices without z. No _core
src_python/tin_engine/cli.py             # --dem, --stride, refusals, shared engine step, docstring
project_structure.md                     # raster section: view landed; module list gains the three
ROADMAP.md                               # row 12; gaps 1, 3, 5
```

## Tests for `@tester`

C++ (`tests/cpp/unit/test_raster.cpp` or a new `test_raster_view.cpp`):

1. `RasterView<float>` and `RasterView<double>` satisfy `RasterSource`;
   `ConstantRaster` still does after it gains `row`.
2. `RasterView` and `Raster` over the same buffer give the same `value_at`,
   `is_nodata` and `bilinear` at a spread of points, including NaN cells and
   the sentinel.
3. `row(i)` has `cols` elements and its first element is `value_at({i, 0})`,
   for the first and last row.
4. Batch sampling: `valid` is false outside the grid, for non-finite points,
   and next to a NoData corner even at zero weight. `z` equals the node value
   at interior nodes, and matches it to 1e-9 relative on the last row and
   column.
5. The view does not copy: a changed buffer value shows through.

Python binding (`tests/python/test_core_raster.py`):

6. `raster_view` accepts C-contiguous float32 and float64. It raises
   `TypeError` for a Fortran-ordered array, a strided slice, int16 and a list,
   and never copies silently.
7. The view keeps the array alive: drop every Python reference to the array,
   run `gc.collect()`, sample, and get the old values.
8. `sample` returns `(z, valid)` with shapes `(N,)`; wrong-shaped points raise
   `ValueError`.
9. The affine scalars are keyword-only: a positional call raises `TypeError`.

Adapter and pure modules:

10. `to_core` forwards exactly `x_min`, `y_max`, `delta_x`, `delta_y`,
    `nodata` and the array. Pin it by sampling a micro-TIFF tile at known nodes,
    and by a signature check that `raster_view` has no `epsg`/`crs`/`rows`/`cols`
    parameter.
11. `grid_domain`: node counts per side for stride 1, a stride that divides
    `n - 1` and one that does not (the last row and column still appear); the
    ring is counter-clockwise, closed implicitly and lies on the extent; the
    default-stride formula gives at most 256 per side; coordinates are
    bit-equal to `x_min + col * delta_x`.
12. `elevation.trim`: drops the right triangles and edges, renumbers
    consistently, returns empty when all are invalid, and changes nothing when
    all are valid.

CLI (`tests/python/test_cli_mesh.py`), on micro-TIFFs from
`tests/python/geotiff_fixtures.py` written to `tmp_path`:

13. `mesh --dem tile.tif --out x.vtk` exits 0. The file's z values equal the
    tile's values at the chosen nodes, and the file carries
    `crs = EPSG:<code>` and an `elevation` field naming the stride.
14. A micro-TIFF with a NoData row on its edge: exit 0, fewer vertices, the
    count in the `elevation` field.
15. All-NoData tile: non-zero exit, no file.
16. Usage errors, exit 2 and no file: both a fixture name and `--dem`; neither;
    `--dem` with `--flat`; `--dem` with `--crs`; a missing file; a file
    `decode_dem` refuses (its message appears in the output); `--stride 0`.
17. `.ply` output with `--dem` has the same z as the `.vtk`.

Real fixture, marked `needs_codecs` (the marker already in
`tests/python/test_io_geotiff.py:90` at design time, since moved to
`tests/python/geotiff_fixtures.py:37` because a second file needs it):

18. `mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --out x.vtk`
    exits 0. Every z is finite and inside [−1.3, 391.8]. The `crs` field is
    `EPSG:25833`. Some vertices are dropped (the NoData outline). In the
    `viewer` job, VTK's own reader loads it with those bounds.
19. Without `imagecodecs` the same command exits 2 and the output names the
    `codecs` extra. (Runs in the default CI job, where the extra is absent.)

**Invariant-critical suite** (mutation testing required, per
`docs/increments/README.md`): tests 2, 4, 7 and 10 — the view agreeing with the
owning raster, the edge semantics of the batch sampler, the lifetime anchor,
and what crosses the boundary. The CLI and trim tests are not.

## LOC estimate

Counted in `CLAUDE.md` §2's unit. An estimate, not a measurement.

| file | what | est. |
|---|---|---|
| `raster/view.hpp` | the view: ctor, four accessors | ~40 |
| `raster/raster.hpp` | concept line, `Raster::row` | ~10 |
| `raster/sample.hpp` | batch sampler | ~15 |
| `bindings/core.cpp` | variant class, factory, `sample`, GIL release | ~80 |
| `_core.pyi` | three stubs | ~15 |
| `raster.py` | `to_core` and the read-only check | ~20 |
| `grid_domain.py` | stride, node coordinates, ring | ~30 |
| `elevation.py` | trim and renumber | ~25 |
| `cli.py` | options, refusals, engine split, fields | ~85 |
| | **total** | **~320** |

Increment 10 overran its estimate by 39 %. On that bias this lands near 450,
still about 250 under the ceiling. Increment 13 came in under (144 against
155), so the true figure is probably between the two.

## Acceptance

- `rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --out tile.vtk`
  (with the `codecs` extra) writes a file ParaView opens as a terrain surface.
  The relief is about 400 m over a 50 km tile, so it looks flat until Z is
  scaled (ParaView's Transform filter). That is a viewing choice, not a defect.
- All gates in `CLAUDE.md` §4 green, CI green.

## Not in scope

- Refinement, error budgets, any adaptive vertex choice (gap 2).
- Clipping to a catchment polygon.
- Mosaics, reprojection, windows (`window_for`).
- A point-data elevation array in the `.vtk`. ParaView can colour by the Z of
  `Points` without one.
- Changing `bilinear`'s NoData rule (zero-weight corners).
