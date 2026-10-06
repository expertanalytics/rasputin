# Project Structure

Layout for the post-CGAL rasputin backend. Each module corresponds to one design doc and one test directory. Public C++ headers live under `include/terrain/`, implementation under `src/`, and the pybind11 module under `bindings/`. The Python layer is the `tin_engine` package under `src_python/` (the distribution is still named `rasputin`); the extension is imported as `tin_engine._core`, the underscore signalling "private — go through the Python API".

## Directory layout

Entries marked *(planned)* do not exist yet; the rest are in the tree today.

```
include/terrain/           # public C++ headers, header-only where possible
  build_info.hpp           # stdlib_hardening(): the bounds-check mode the standard
                           #   library reports; _core.hardening and the guard (24)
  core/
    point.hpp              # Point2 value type, dot and cross
    bbox.hpp               # Box2; empty-box identity, exact closed containment
    segment.hpp            # Segment2 and on_segment<K>
    ring.hpp               # IndexedRing, orientation<K> of a ring
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
    chunks.hpp             # for_each_block (dynamic blocks) over std::jthread
  mesh/
    lattice_mesh.hpp       # LatticeMesh: flat triangle array over DEM nodes,
                           #   neighbour links, the three splits (14), flip (14b)
    lawson.hpp             # LatticeFrame, legalise_around / legalise_all:
                           #   Lawson flips on strictly-inside apexes (14b)
    quality.hpp            # improve(): minimum-angle Steiner nodes on the start
                           #   mesh, DEM node nearest each circumcentre (20)
  refinement/
    scan.hpp               # per-triangle sup-norm scan, NoData carve point (14)
    refine.hpp             # RefineOptions, RefineOutcome, the round loop (14),
                           #   Delaunay insertion (14b), the quality-start call
                           #   (20), constraint feet (20b)
    check_points.hpp       # CheckPoints: source nodes filed by target-grid
                           #   cell, for the final check (15c)
    refine_points.hpp      # refine_points: greedy insertion against stored
                           #   check points (15c, the final check); with the
                           #   edge strip, and refine_strip (15f)
    strip_scan.hpp         # the edge strip inside the loop: sub-edge map,
                           #   scan by owned sub-edge, guarded insertion (15f)
    constraint_points.hpp  # constraint_check_points, ConstraintCheckPoints:
                           #   where constraint edges cross grid lines, and
                           #   the midpoints between, filed by edge (15f)
  hydrology/
    flood.hpp              # detail::flood: the one Priority-Flood of both (29)
    upstream.hpp           # upstream(z, seed): the seed set's catchment (22)
    accumulate.hpp         # accumulate(z): count, reach, flow_to per node (29)
  vector_simplify/
    area_collapse.hpp      # reduce_ring: area-preserving segment collapse
                           #   (Kronenfeld et al. 2020) to a horizontal
                           #   tolerance, exact crossing tests, keep-points (22)

src/                       # C++ implementation, one directory per module
                           #   (only predicates/ and cdt/ exist; rest planned)
  predicates/              # exact orient2d/incircle; namespace terrain::pred
  parallel_util/           # (none: header-only, include/terrain/parallel_util/)
                           # (no vector_simplify/ or hydrology/ here: both
                           #  are header-only, include/terrain/...)
                           # (no noding/ here, and none planned: the noder is
                           #  header-only, its driver a template on the kernel)
  cdt/                     # thin wrapper over vendored Detria
  mesh/                    # (none: header-only, include/terrain/mesh/)
  refinement/              # (none: header-only, include/terrain/refinement/)
  flip/                    # final Lawson edge-flip pass, constraint-respecting

bindings/
  core.cpp                 # pybind11 module definition -> tin_engine._core

src_python/tin_engine/     # public Python API (distribution name: rasputin)
  __init__.py              # installed_version only; never imports _core
  cli.py                   # Typer entry point declared in pyproject;
                           #   `rasputin palette NAME [--out FILE]` writes a
                           #   ParaView colour preset (16c)
  raster.py                # the ONLY adapter from decoded data into _core
  grid_domain.py           # DEM extent -> stride-subsampled nodes + outer ring;
                           #   pure numpy, never imports _core
  mosaic.py                # plan_mosaic / assemble: select, group by lattice,
                           #   check overlaps and coverage, stitch; no files (15a);
                           #   imports only io.models first-party
  dem_input.py             # --dem/--bbox or a domain -> DemInput(tile, plan,
                           #   label, domain in the DEM's CRS) (15a, 15b);
                           #   takes a DemRepository and TileFootprints, and
                           #   the union of the tiles' node_box
  domain.py                # DomainPolygon (one polygon in its own CRS) and its
                           #   to_crs, check_extent, DomainError; reads no file
                           #   (io/domain_file.py reads --domain); never imports
                           #   _core
  crs.py                   # parse_crs, reprojector: the one Transformer.from_crs
                           #   site, always_xy (15b); pyproj and numpy
  target_grid.py           # TargetGrid (its node_box) on one global lattice,
                           #   source_region, SourceWindows / TileWindows,
                           #   resample (bilinear, threads), check_point_blocks;
                           #   NoData by valid_mask (+-inf is data); never
                           #   imports _core or mosaic (15c-2)
  final_check.py           # run: phase 2, the source nodes filed in a
                           #   CheckPoints, then refine_points from phase 1's
                           #   mesh (15c-2); takes the edge strip, checked in
                           #   the same loop (15f-3)
  edge_strip.py            # generate / run: the edge strip's two _core calls
                           #   (constraint_check_points, refine_strip) and
                           #   their clock rows; holds no geometry (15f-3)
  elevation.py             # drops mesh vertices the DEM has no data for;
                           #   pure numpy, never imports _core
  stats.py                 # --stats: PhaseClock, quality, Report, render to
                           #   Markdown; numpy only, no _core, no typer
  features.py              # EdgeVocabulary: which bit means which feature.
                           #   The names C++ refuses to hold. Imports nothing
                           #   first-party and never imports _core
                           #   PIECE_VOCABULARY: the default plus `seam`
                           #   (bit 9), for a cut run's piece files (23c).
                           #   TerrainFeature: one feature's mask and lines
                           #   (moved from feature_input, audit PR C)
  decompose.py             # partition: the window cut into cells under
                           #   --pieces and --memory-budget, b(T); pure, reads
                           #   no machine (23c)
  feature_input.py         # --features: a GeoJSON, GeoPackage layer or GML
                           #   read, mapped to masks by a ClassMap, pre-clipped,
                           #   moved to the DEM's CRS and clipped to the domain
                           #   as linework -> FeatureSet (16b); opens .gml
                           #   itself, reads GeoJSON through io/repository's
                           #   read_json and io/geojson's read_collection;
                           #   never imports _core. read_source: one file's raw
                           #   rows (16b's reader, shared); read_lake_polygons:
                           #   the polygons of `catchment --lakes` near the
                           #   seed point, in the file's CRS (22)
  outline.py               # trace(mask): marching-squares rings between in-
                           #   and out-nodes, (8, 4) saddle rule; numpy only,
                           #   never imports _core (22)
  catchment.py             # CatchmentRequest -> delineate(request, repo) ->
                           #   Catchment: seed, the window loop over 15a's
                           #   plan, the flood, the fine ring and its
                           #   reduction (via catchment_core, never _core);
                           #   takes a DemRepository and no path (22); with a
                           #   river reach, the pour point is the gauge's node
                           #   on the burnt reach (29)
  catchment_core.py        # upstream, accumulate: the core call on
                           #   raster.to_core's view of a DemTile; reduce_ring
                           #   re-exported; the catchment's _core calls, as
                           #   edge_strip.py is the edge strip's (python-
                           #   audit.md, section 11)
  hydrography.py           # RiverSegment, Station, Lake: the hydrography's
                           #   value types; no first-party import, so gauge
                           #   places on them without a codec, and io/rivers.py
                           #   and io/station_set.py read into them (python-
                           #   audit.md, section 11)
  gauge.py                 # place: a station's foot P on NVE's river lines,
                           #   by tier, and the reach round it (29, PR 2);
                           #   lake_seed -> LakeSeed: the gauges seeded with
                           #   their lake, by containment and one distance
                           #   (PR 4); pure shapely, no DEM, no file, no codec
  burn.py                  # burn_reach: the reach moved onto the window's
                           #   valley floor and burnt in -> (burnt copy,
                           #   GaugePath); numpy only (29, PR 2)
  sensitivity.py           # assess -> Sensitivity: is the area well defined
                           #   along the burnt chain within U; arrays in,
                           #   no DEM, no shapely (29, PR 2)
  reference.py             # agreement, classify, summarise: our catchment
                           #   against NVE's polygon, counted on the DEM's
                           #   node lattice; pure, no file (29, PR 4)
  catchment_batch.py       # seed_for: place, lake_seed, the request (shared
                           #   with catchment --rivers); run_batch: delineate
                           #   in a thread, compare, StationResult rows to a
                           #   BatchSink; no paths (29, PR 4)
  landcover.py             # regions, label_triangles: a land-cover code per
                           #   triangle, components across unconstrained edges,
                           #   one point-in-polygon test per component (16c);
                           #   numpy and shapely, never imports _core
  palettes.py              # CORINE_NATURAL (code -> label, colour) and
                           #   paraview_preset(); data, imports nothing
                           #   first-party (16c)
  chains.py                # start_chains: the domain's rings, then every
                           #   feature line, as (indices, role, mask) (16b);
                           #   never imports _core
  sources.py               # the catalogues `rasputin fetch` and `rasputin
                           #   fetch-stations` copy from: RemoteSource (ANADEM,
                           #   GLO-30; 23a-1, 23a-2), StationSource (NVE's HRD,
                           #   29), notice(); data, imports Pydantic alone, so
                           #   the mesh path names a key without importing fetch/
  fetch/                   # copies remote data to disk; the only package that
                           #   touches the network, and nothing on `rasputin
                           #   mesh`'s path imports it (23a-2)
    __init__.py
    http.py                # RangeClient: byte ranges over HTTP, stdlib only,
                           #   retries; `User-Agent: rasputin/<version>` (29);
                           #   query_url; the one module importing urllib
    plan.py                # FetchRequest -> the blocks meeting the domain's
                           #   box, from headers alone; pure (23a-2)
    run.py                 # fetch(): headers, plan, bounded async download;
                           #   client and cache writer injected (23a-2)
    nve.py                 # fetch_station_set(source, get_text) -> each file's
                           #   bytes (stations, references, rivers, lakes,
                           #   manifest last); field allow-lists; writes
                           #   nothing, cli.py writes the files (29)
  data/                    # package data, read by importlib.resources
    nve_hrd_2025.csv       # NVE's 140 HRD stations, 2025 version: which
                           #   stations, series version, start year; its names
                           #   only check the extraction from the PDF (29)
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
  io/                      # all file decoding AND encoding lives here, with
                           #   known exceptions, encoded in cli.py with the
                           #   standard library: `palette`'s JSON (16c), to
                           #   move here in a follow-up, and station-
                           #   catchments' results.csv and summary.json (29)
    __init__.py
    geojson.py             # read_collection: a parsed GeoJSON object ->
                           #   (features, crs text), the one crs and shape rule
                           #   for every GeoJSON reader (audit PR C);
                           #   feature_collection: the document both writers
                           #   build; catchment_geojson: a polygon, its CRS and
                           #   its properties -> the catchment file's bytes
                           #   (22's writer, moved here in 29); opens nothing
    domain_file.py         # read_domain: --domain's one Polygon (GeoJSON or
                           #   WKT, by suffix) in its own CRS, checked and
                           #   oriented -> DomainPolygon; reads through
                           #   repository's read_json / read_text, opens
                           #   nothing itself (16; moved from domain.py in
                           #   audit PR C)
    ply.py                 # arrays -> PLY bytes; takes no path and opens nothing;
                           #   face_codes= adds the face property
                           #   land_cover_code (16c)
    vtk_legacy.py          # arrays + EdgeVocabulary -> legacy .vtk bytes, for
                           #   ParaView; takes no path and opens nothing;
                           #   triangle_codes=, land_cover_codes= add the cell
                           #   array land_cover_code and its field (16c)
    geotiff.py             # TIFF container + GeoKey decoding -> DemTile
    models.py              # frozen Pydantic Bounds (`of` a shapely-order
                           #   box), IndexWindow, RasterMeta, TileFootprint,
                           #   DemTile; RasterMeta's node arithmetic, the
                           #   core's spelling: node_xy, index_of, node_box,
                           #   windowed; valid_mask, the one NoData rule (NaN
                           #   or the sentinel; +-inf is data); imports
                           #   nothing first-party
    geopackage.py          # GeoPackage layer_info / query_features over an
                           #   open sqlite3.Connection; frozen dataclasses,
                           #   opens nothing, knows no path (16b)
    gml.py                 # read_gml: OGR-written GML2 from a binary stream,
                           #   standard library XML; opens nothing (16b)
    cog.py                 # decode_window: only the blocks a window meets,
                           #   from a BlockSource (an open file or the tile
                           #   cache), its meta `meta.windowed(window)`;
                           #   opens nothing (23a-1)
    station_set.py         # read_stations -> (stations, crs) and
                           #   read_references -> ({station: polygon}, crs):
                           #   fetch-stations' files or a user's points file,
                           #   read by io/geojson's read_collection with no
                           #   default CRS, so the `crs` member is required;
                           #   duplicates and wrong geometry refused with
                           #   ValueError (29); read_nve_lakes -> (Lake parts,
                           #   crs): `station-catchments --lakes`, a
                           #   MultiPolygon split (29, PR 4)
    rivers.py              # kind_of (lake or river, total over NVE's
                           #   objekttype spellings), drop_copies,
                           #   read_segments -> (RiverSegments, crs, copies
                           #   dropped); one LineString per segment (29)
    repository.py          # the ONE io/ module that opens files:
                           #   DemRepository, the Protocol: footprints, load,
                           #   load_window, check; TiffDemRepository lists
                           #   headers, loads tiles (15a); open_geopackage, a read-only SQLite
                           #   connection (16b); CacheRepository and
                           #   CachedBlocks read the tile cache (23a-1);
                           #   CacheWriter writes it under a lock (23a-2);
                           #   read_json, a JSON document, read-only (29);
                           #   read_text, UTF-8 text (audit PR C's WKT domain)
    mesh_index.py          # MeshIndex (a cut run's index.json, frozen,
                           #   unknown keys refused), SeamRecord,
                           #   check_conformity (K4); opens nothing (23c)

tests/
  cpp/                     # C++ tests (Catch2; unit/ and property/)
  python/                  # Python tests (pytest; hypothesis planned)
  fixtures/                # in-repo test data (DEM, CORINE); credits in NOTICE.md

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
predicates ───────────┬─→ noding ──→ vector_simplify   hydrology  (raster only)
                      │   (area_collapse uses noding's classify<K>
                      │    and the predicates' kernel)
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
  projected, metre CRS, and never passes a geographic one on. Since 15c-2
  `io/geotiff.py` reads a geographic DEM, but `raster.to_core` refuses a
  geographic tile: the DEM is resampled in Python onto a square grid in the
  projected, metre `--out-crs` (`target_grid.py`), and only that grid and the
  source nodes moved into it cross. Since increment 15b the domain polygon
  is reprojected in Python too (`crs.py`, the one `from_crs` site); see
  `docs/increments/15c-geographic-dem.md` D6 and
  `docs/increments/11-raster-ingestion.md` §9, §10.
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

Header-only. One helper in `chunks.hpp`, over `std::jthread` created per call and joined before it returns, no pool: `for_each_block(n, threads, BlockSchedule, fn)`, blocks handed out from one atomic counter, run inline below `BlockSchedule::inline_below`, which refine's scan uses (increment 21a). `std::execution::par` and OpenMP were both ruled out in increment 14 (R7): neither builds on macOS without an experimental flag or an extra runtime. Needs only `Threads::Threads`.

### `vector_simplify`

Header-only. `area_collapse.hpp` (increment 22): `reduce_ring`, Kronenfeld, Stanislawski, Buttenfield and Brockmeyer's area-preserving segment collapse (APSC, 2020): each collapse replaces two vertices by one on the line that keeps the area, so the area is kept up to rounding; collapses are taken least deviation first while the fine ring's vertices stay within a horizontal tolerance; each new edge is tested for crossings with `noding::classify<K>` over a uniform grid, and keep-points (the seed) must stay inside. Serial and deterministic. Visvalingam-Whyatt and Douglas-Peucker are not used: neither keeps area. See `docs/increments/22-auto-catchment.md`.

### `hydrology`

Header-only; depends on `raster` only. `flood.hpp` (increment 29) holds the one Priority-Flood (Barnes, Lehman and Mulla 2014), `detail::flood(z, state, on_reach)`, from the window's edge and from the nodes beside NoData, ties first in, first out; it calls `on_reach(j, j)` for each outlet and `on_reach(i, j)` when popped node `i` reaches node `j` first, and returns the outlets beside NoData. Both functions below run it, so they cannot drift. Each node drains to the node that flooded it, its lowest filled neighbour, so there is no separate fill, no flat resolution and no D8 pass.

`upstream.hpp` (increment 22): `upstream(z, seed)`, a template on the `RasterSource` concept, labels a node *in* when it is a seed or was flooded from an *in* node. Returns the mask, the in-nodes' count and bounds, and `touches_edge` / `touches_nodata` (an in-node that is, or neighbours, an outlet of that kind).

`accumulate.hpp` (increment 29): `accumulate(z)`, the accumulation pass. One flood records each node's flooder (`flow_to`, 255 for an outlet and on NoData) and the push order; one sweep in reverse push order adds each node's count and reach bits to its flooder. For every node `c` with data, `count[c]` and the two `reach` bits equal `upstream(z, {c})`'s `nodes_in`, `touches_edge` and `touches_nodata`. Refused with `std::length_error` (a `ValueError` in Python) at 2^32 nodes or more, before any per-node array is allocated, since the count and order are 32-bit.

Serial; the bindings release the GIL. What `auto_catchments.md` also sketches (epsilon filling, D8, streams and Strahler order, RichDEM) is not built; the outline tracer is Python (`outline.py`). See `docs/increments/22-auto-catchment.md` and `docs/increments/29-nve-reference-catchments.md`.

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

Refinement against the DEM to a sup-norm tolerance (increment 14). `scan.hpp` measures a triangle's largest `|z - plane|` over the DEM nodes it contains, by exact integer tests; `refine.hpp` runs rounds of parallel scan (`parallel_util`) and serial splits in index order — a fan for an interior node, an edge split on both sides for an edge node — each followed by Lawson flips (`mesh/lawson.hpp`, increment 14b) that never cross a constrained or boundary edge, with every written slot rescanned, so the output is constrained Delaunay and within tolerance, and does not depend on the thread count. Triangles with a NoData vertex are carved, not refined. Before the first scan, an optional minimum-angle pass (`mesh/quality.hpp`, increment 20) improves the start mesh; during refinement, a worst node within ε of a constraint segment is replaced by its foot on the segment (increment 20b, see `docs/GLOSSARY.md`). See `docs/increments/14-adaptive-refinement.md`.

### `flip`

Final Lawson edge-flip pass. Skips constraint-tagged edges. Parallel with edge-conflict resolution.

### `bindings/core.cpp`

Single pybind11 module that exposes the C++ API to Python. Built as the `_core` extension, installed into the `tin_engine` package. Nothing is re-exported: `src_python/tin_engine/__init__.py` holds only `installed_version` (`__all__`), so the package imports without the extension; the CDT surface below is reached as `tin_engine._core` and is deliberately not re-exported, because `cli.py` is the sole composition root and `viz/` never imports the extension. Typed from `src_python/tin_engine/_core.pyi`, which is what `mypy --strict` sees.

As of increment 6a (`94f94e2`, `docs/increments/06-cdt-viewer.md`) it binds:

- no point value types: coordinates cross as `(N, 2)` float64 arrays (the Python `Point2`, `Point3`, `dot` and `cross` went in PR B of `docs/increments/cpp-audit.md`, section 7);
- no free functions over points: the C++ `dot` and `cross` in `point.hpp` stay inside the core;
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
- **External deps under consideration:** RichDEM (optional, MIT; increment 22 does not use it), Eigen (if linear algebra needs grow beyond what we want to hand-roll).
- **No CGAL and no GDAL** in the new core: prohibited by `CLAUDE.md` §2.
  **No Boost.Geometry** either, which is a scope choice and not a prohibition
  — reading this line as one is what put an unauthored ban in §2 for six
  days. Existing Python-layer uses are migrated incrementally.
- **Planned:** per-module `OBJECT` libraries linked into the extension once `src/` is populated, and sanitizer flags via `RASPUTIN_SANITIZER=asan|ubsan|tsan|none`.

Packaging is driven by **scikit-build-core**, declared in `[build-system]` in `pyproject.toml`, which invokes this same CMake build. `pip install .` configures with `RASPUTIN_BUILD_PYTHON=ON` and `RASPUTIN_BUILD_TESTS=OFF` — the C++ tests pull Catch2 over the network and have no business running during an install — and installs `_core` into the `tin_engine` package.

There is no `setup.py`. It was removed with the foundation reset, along with the `CMakeBuild` class it carried.

## Python API surface

The public API is the `tin_engine` package, calling into `tin_engine._core`. `tin_engine.viz` is the renderer that turns a triangulation into an SVG a person can look at; it consumes the `typing.Protocol`s in `viz/protocols.py` and **never imports `_core`**, so it is testable with no compiled extension in the process. `cli.py` is the single composition root that joins the two -- the same shape as the rule below that exactly one module adapts decoded raster data into `_core`. The `rasputin draw` command drives it: it looks a fixture up in `viz.fixtures.GALLERY`, maps that fixture's **string** chain roles onto `_core.ChainRole` -- `viz/` may not name the enum, so the mapping is the composition root's -- runs `build_pslg` and `triangulate`, and hands the fixture itself to `build_scene` as the `PslgLike`, which is what lets a fixture the validator rejects still be drawn. It is also the only place a path exists, and it resolves and refuses one before writing.

The catchment file's bytes come from `io/geojson.py`'s `catchment_geojson`, which opens nothing; `rasputin catchment` and `rasputin station-catchments` both call it, and `cli.py` keeps only the write.

The pre-migration `rasputin.*` modules (`mesh.py`, `geometry.py`, `reader.py`, `tin_repository.py` and friends) are not in the working tree. They, with the CGAL-based `triangulate_dem.h` and `bindings.cpp`, left it in the release-hygiene PR and are kept in history under the annotated tag `legacy-archive` (`docs/increments/release-hygiene.md`, section 3). A shape worth preserving is read back with `git show legacy-archive:legacy/rasputin/<file>` before it is reintroduced, and an increment's "Legacy" section greps the tag. `tools/check_citations.py` resolves a `legacy/…:N` citation through the tag.

## What was removed

- The CGAL-era pipeline: `legacy/` (30 files, 5,224 lines), with `tools/check_legacy_imports.py` and `@migration-expert`, which existed only for it. History keeps all three (`legacy-archive`).
- CGAL, GMP, MPFR from the CMake dependency list.
- Boost.Geometry, whose only uses were in `legacy/`. It was never a prohibited dependency (`CLAUDE.md` §2 is the list), only an unused one.
- `lib/date/`. The C++20 `<chrono>` calendar types replaced it; `CLAUDE.md` §2 prohibits external `date` libraries.
