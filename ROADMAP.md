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
| 12 | DEM to elevated mesh: the zero-copy `RasterView` and its binding, z sampled per mesh vertex, `rasputin mesh --dem file.tif` | shipped (#91) | `docs/increments/12-dem-to-mesh.md` |
| 13 | One mesh file for ParaView: legacy `.vtk` with triangles, constraint edges, feature masks and the vocabulary, text by default | shipped (#90) | `docs/increments/13-bundled-mesh.md` |
| 14 | Adaptive refinement: split triangles at the worst DEM node until every triangle's max error is within `--tolerance`; parallel scan, serial deterministic splits | shipped (#92) | `docs/increments/14-adaptive-refinement.md` |
| 14b | Delaunay insertion: Lawson flips after each refinement insertion, never across a constraint, flipped triangles rescanned; the output is constrained Delaunay and still within `--tolerance` | implemented, in review | `docs/increments/14b-delaunay-insertion.md` |
| 16 | Mesh a domain polygon: `--domain` (GeoJSON or WKT, one polygon with holes, CRS equal to the DEM's), vertices kept where they are with bilinear z (DEM values are point heights), the start mesh is the CDT of its rings alone, then refinement inserts DEM nodes as in 14b. Pulled forward from the catchment-clip entry at the user's request | designed, choices open | `docs/increments/16-domain-polygon.md` |
| — | `raster/`: grid-to-world geometry and bilinear sampling | shipped (`7785fea`), **no record** | none — predates the protocol |

## What stands between here and an operational MVP

An MVP is one command turning a DEM and a catchment polygon into a terrain TIN
file. Six things were missing, in dependency order. Gap 4 has shipped; gaps 1
and 3 and the first half of 5 shipped as increment 12 (#91); gap 2's first half
shipped as increment 14 (#92), and 14b, which fixes its triangle shape, is
implemented and in review.

1. **Raster ingestion, Python side.** The C++ `raster/` module samples; nothing
   decodes a GeoTIFF into it. `project_structure.md` names `raster.py` as the
   only adapter from decoded data into `_core`. This is
   where CRS stops. Split in two: increment 11 decodes
   (`io/geotiff.py`, bytes to a validated tile, no `_core`); the adapter,
   `RasterView` and the pybind11 buffer surface shipped as increment 12.
2. **Refinement.** The largest piece and the actual product: refine a coarse
   triangulation by looking up the DEM until every triangle is within a
   tolerance, instead of triangulating every DEM node. Plan in
   `parallel_refinement.md`. Increment 14 does it to a sup-norm tolerance
   (`mesh --dem --tolerance`); 14b replaces the fans with Delaunay insertion; a flat-ground size cap is later.
3. **Elevation assembly.** `IndexedMesh2` is 2D by design and z comes from
   sampling the raster per vertex. Both halves exist; increment 12 joins
   them.
4. **Mesh output.** Shipped as increment 10 (PLY, for QGIS) and extended by
   increment 13: `rasputin mesh --out x.vtk` writes one text file for ParaView
   holding the triangles, the constraint edges, their feature masks and the
   feature names. Real z is gap 3, in increment 12.
5. **A CLI that does the job.** `rasputin` has `version`, `draw` and `mesh`,
   and all three take a built-in fixture. Increment 12 adds the first path
   from a file on disk to a mesh on disk, `mesh --dem file.tif`, over a regular
   subsample of the DEM's nodes until refinement (gap 2) replaces it.
6. **A DEM in several tiles.** `mesh --dem` takes one `.tif`, but a DEM usually
   arrives as many. Increment 11 set the mosaic aside as its own increment
   (ruling 10) and it never reached this list; the user caught the gap on
   2026-09-25. First cut: a DEM *repository* over a directory (the legacy
   `RasterRepository`'s shape, which the user likes), reading each tile's
   extent from its header, selecting the tiles that cover the area, refusing
   mixed CRS, cell size or grid alignment, and stitching the selection into one
   grid. Aligned tiles are one bigger grid, so everything downstream is
   unchanged. Planned after increment 14; no record yet.

Open and not MVP-blocking: **inputs in their own CRS**. The user (2026-09-26):
"I don't think the domain CRS should have to match the DEM CRS in the future.
This must be written down. The questions will be in which coordinate system we
shall do the math." Increment 16 requires a match for now; a later increment
chooses the computation CRS and reprojects inputs into it (see
`docs/increments/16-domain-polygon.md`, "Ruled by the user"). Also **a general point insertion policy**, which the user
wants to discuss (2026-09-26): which points refinement inserts, beyond
increment 14b's worst DEM node, and how that meets the sizing field, coastline
constraints and points that are not DEM nodes. The user's direction for it
(2026-09-26): input is geometry, never loose points; every vertex comes either
from the input polygons and polylines, used as given, or from refinement against
the DEM to `--tolerance`. Nothing else. Increment 16's U2 (what an input vertex
that is not a DEM node becomes) was the first concrete question in it; the
user ruled that input vertices stay where they are, with bilinear z. Also:
coarsening input geometry (`parallel_refinement.md` step 2), 5d above, and
land-cover partitioning, whose foundation is increment 7's property sets and
whose consumer does not exist.

Known defect: increment 14's NoData carving from a NoData corner appears to
advance one node per round. Meshing the real tile from its outline alone took
5 054 rounds and 109 s against 0.35 s from the stride grid
(`docs/increments/16-domain-polygon.md`, M4 and R5). It needs its own look.

Clipping to the catchment polygon left this list on 2026-09-26: it is increment
16, pulled forward at the user's request. Interior polygons and polylines follow
as 16b, designed in the same record.
