# Rasputin roadmap

Rasputin builds terrain TINs: a digital elevation model and a set of linear and
area features go in, a constrained triangulation comes out, coarsened under an
error budget. It is a rewrite of a CGAL-based pipeline onto MIT-licensed
components and predicates of its own.

This table is the index. Each increment's record under `docs/increments/` holds
its design, its rulings and its measurements; `docs/increments/README.md` is the
protocol those records follow.

**Authoring order is not ship order.** Increment 6 is designed sixth and ships
before 5b, because the viewer depends on nothing the noder produces and 5b is the
increment that most needs a picture to check against
(`docs/increments/06-cdt-viewer.md`, "The ordering ruling").

| # | What it is | Status | Record |
|---|---|---|---|
| 1 | Exact geometric predicates | shipped (`d0946bc`, `aeb667a`) | `docs/increments/01-predicates.md` |
| 2 | 2D geometry value types | shipped (`1a60ad4`) | `docs/increments/02-core-geometry.md` |
| 3 | The PSLG and its validator | shipped (`d984ff8`) | `docs/increments/03-pslg.md` |
| 4 | The constrained Delaunay triangulation | shipped (`f22ddd3`) | `docs/increments/04-cdt.md` |
| 5a | The noder's numeric floor: snap grid and hot pixels | shipped (`e090909`, #69) | `docs/increments/05-noder.md` |
| 5b | The noder proper: segment splitting and the noded PSLG | shipped (`93fc772`, #76) | `docs/increments/05b-noder-driver.md` |
| 5c | Wiring the noder through: CDT signature, bindings, CLI | shipped (#77) | `docs/increments/05c-noder-wiring.md` |
| 5d | The corner graze: a single-point cell touch stops counting | designed, **unscheduled** | `docs/increments/05d-corner-graze.md` |
| 6a | The pybind11 CDT surface | shipped (`311459b`, #71) | `docs/increments/06-cdt-viewer.md` |
| 6b-i | `build_scene`: the mesh-and-PSLG-to-geometry mapping | shipped (`ab4681c`, #72) | `docs/increments/06-cdt-viewer.md` |
| 6b-ii | The SVG renderer, the fixture gallery and `rasputin draw` | shipped (`4834568`, #74) | `docs/increments/06-cdt-viewer.md` |
| 7 | Edge property sets, replacing the one-bit `is_river` | shipped (`2e7577c`, #75) | `docs/increments/07-edge-properties.md` |
| 8 | The crossing gallery: roads into forests, bridges over lakes, structures leaving a catchment | shipped (`6d68d9b`, #81) | `docs/increments/08-crossing-gallery.md` |
| 9 | Gallery output: `rasputin gallery` renders all eleven fixtures into a directory the caller names | designed | `docs/increments/09-gallery-output.md` |
| 10 | Mesh output: a PLY writer, binary `double` by default, constraint edges in a second file | shipped (#86) | `docs/increments/10-mesh-output.md` |
| 11 | Raster ingestion, the decode half: GeoTIFF bytes to a validated `DemTile` | shipped (#89) | `docs/increments/11-raster-ingestion.md` |
| 13 | One mesh file for ParaView: legacy `.vtk` with triangles, constraint edges, feature masks and the vocabulary, text by default | designed | `docs/increments/13-bundled-mesh.md` |
| — | `raster/`: grid-to-world geometry and bilinear sampling | shipped (`7785fea`), **no record** | none — predates the protocol |

## What stands between here and an operational MVP

An MVP is one command turning a DEM and a catchment polygon into a terrain TIN
file. Five things are missing, in dependency order. None has an increment record
yet, so each starts at `docs/increments/README.md` step 1.

1. **Raster ingestion, Python side.** The C++ `raster/` module samples; nothing
   decodes a GeoTIFF into it. `project_structure.md` names `raster.py` as the
   only adapter from decoded data into `_core`, and marks it planned. This is
   where CRS stops. Split in two: increment 11 decodes
   (`io/geotiff.py`, bytes to a validated tile, no `_core`); the adapter,
   `RasterView` and the pybind11 buffer surface are increment 12, which does
   not yet have a record.
2. **Refinement.** The largest piece and the actual product: coarsen a dense DEM
   under an error budget instead of triangulating what you are given. Plan in
   `parallel_refinement.md`. The legacy `-ratio 0.4` was this knob.
3. **Elevation assembly.** `IndexedMesh2` is 2D by design and z comes from
   sampling the raster per vertex. Both halves exist; nothing joins them.
4. **Mesh output.** Shipped as increment 10: `rasputin mesh` writes PLY, binary
   `double` by default, with the constraint edges and their feature masks in a
   second file. What remains is real z, which is gap 3.
5. **A CLI that does the job.** `rasputin` has `version`, `draw` and `mesh`,
   and all three take a built-in fixture. There is still no path from a file on
   disk to a mesh on disk: that needs gap 1.

Open and not MVP-blocking: clipping to the catchment polygon
(`auto_catchments.md`, and `legacy/rasputin/geometry.py`'s Shapely intersection
is the prior art), 5d above, and land-cover partitioning, whose foundation is
increment 7's property sets and whose consumer does not exist.
