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
| 9 | Gallery output: `rasputin gallery` renders all eleven fixtures into a directory the caller names | designed, **unscheduled** (Ola, 2026-09-27: not important now) | `docs/increments/09-gallery-output.md` |
| 10 | Mesh output: a PLY writer, binary `double` by default, constraint edges in a second file | shipped (#86) | `docs/increments/10-mesh-output.md` |
| 11 | Raster ingestion, the decode half: GeoTIFF bytes to a validated `DemTile` | shipped (#89) | `docs/increments/11-raster-ingestion.md` |
| 12 | DEM to elevated mesh: the zero-copy `RasterView` and its binding, z sampled per mesh vertex, `rasputin mesh --dem file.tif` | shipped (#91) | `docs/increments/12-dem-to-mesh.md` |
| 13 | One mesh file for ParaView: legacy `.vtk` with triangles, constraint edges, feature masks and the vocabulary, text by default | shipped (#90) | `docs/increments/13-bundled-mesh.md` |
| 14 | Adaptive refinement: split triangles at the worst DEM node until every triangle's max error is within `--tolerance`; parallel scan, serial deterministic splits | shipped (#92) | `docs/increments/14-adaptive-refinement.md` |
| 14b | Delaunay insertion: Lawson flips after each refinement insertion, never across a constraint, flipped triangles rescanned; the output is constrained Delaunay and still within `--tolerance` | shipped (#93) | `docs/increments/14b-delaunay-insertion.md` |
| 16 | Mesh a domain polygon: `--domain` (GeoJSON or WKT, one polygon with holes, CRS equal to the DEM's), vertices kept where they are with bilinear z (DEM values are point heights), the start mesh is the CDT of its rings alone, then refinement inserts DEM nodes as in 14b. Pulled forward from the catchment-clip entry at the user's request | shipped (#94) | `docs/increments/16-domain-polygon.md` |
| 17 | `rasputin mesh --stats PATH` (or `-` for stdout): a Markdown report of sizes, plan-view quality (min angle, vertex degree), refinement counters and per-phase timings, with `refine`'s legalise, scan and split-and-flip times from new `RefineOutcome` fields | shipped (#95) | `docs/increments/17-mesh-stats.md` |
| 18 | Row-span scan: refinement visits only the DEM nodes inside each triangle, row by row over contiguous memory (exact integer bounds for node-only triangles, a double estimate corrected by the exact predicate otherwise). Rows are walked as tile segments so a later multi-tile mosaic is additive. Output bit-identical to 17 | shipped (#96) | `docs/increments/18-row-span-scan.md` |
| 20 | Quality start: before DEM refinement, Steiner points at the DEM node nearest each bad triangle's circumcentre until the start mesh has a 25° minimum angle (`--start-min-angle`, 0 is off); geometry only, input segments not split, serial and deterministic. Removes increment 16's boundary fans | landed with 20b as interim (C1-C3 open, to 20c; C4 -> 20b) | `docs/increments/20-start-quality.md` |
| 20b | Minimum insertion distance from constraints: when refinement's worst node lies within ε = clamp(tol / slope, cell/100, cell/2) of a constraint segment, insert the foot on the segment (off-node, bilinear z) instead; a footed node still above tolerance is inserted after all, so the tolerance guarantee is unchanged. Removes the 0.0117° needle at 1 m (increment 20's C4) | shipped with 20 (#97) | `docs/increments/20b-min-insertion-distance.md` |
| — | `tools/bench.py`: the 1 m benchmark and the thread-scaling sweep from one checked-in command, with power state, quality and commit recorded per run; the one-off scripts in `docs/benchmarks/2026-09-26/` are its specification. Rule 2's acceptance run needs it | shipped with branch `tools-bench`'s PR | `docs/benchmarks/bench-py.md` |
| — | The serial phase: profile refine's serial insert-and-flip phase, then parallelise what the profile blames. Scaling tops out at about 2.0-2.2×. Profiled 2026-09-27: serial part about a third of 1-thread refine, mostly Lawson legalisation; the scan stops near 5× from load imbalance | profiled; designed as increment 21 (Ola's rulings 2026-09-27: L1 determinism, one path, at most 2 % more triangles). **21a shipped with branch `increment21a-quick-wins`'s PR**: dynamic scan blocks, active merge, reused flip stack; mesh bit-identical; on AC refine -10 to -15 % at 8 threads, ceiling 2.1x -> 2.3x (`docs/benchmarks/2026-09-27/21a-acceptance.md`). **21b shipped with branch `increment21b-lattice-incircle`'s PR**: an int64 lattice incircle answers 99.97-100 % of refine's incircle tests; mesh bit-identical; on battery refine -15 % at 8 threads, ceiling 2.3x -> 2.5x (`docs/benchmarks/2026-09-27/21b-acceptance.md`). **21c measured** (`docs/benchmarks/2026-09-27/21c/README.md`): option C costs 6.5-10 % more triangles; A1 about 1 %; A0 is bit-identical (on these two inputs, not proven), and with evaluate-once is modelled at only 7-11 % faster refine at 8 threads and slower at 4. **21d is deferred** (Ola, 2026-09-27: "Review and push 21c, then basin-work"), behind the work in "Order of work" below; domain decomposition parallelises scan and split together | `docs/increments/21-parallel-refine.md`, `docs/benchmarks/2026-09-27/serial-profile/README.md` |
| — | Release hardening: measure libc++'s `_LIBCPP_HARDENING_MODE_FAST` (and libstdc++'s assertions on the GCC leg), which bounds-check `std::vector`, `span` and the like in Release, on the 1 m benchmark; switch it on if the cost is small. Today an out-of-range index in shipped code the tests miss is undefined behaviour and crashes Python (Ola, 2026-09-27: "I'm surprised we don't have proper memory control"; CI's ASan/UBSan/TSan cover what the tests reach) | to measure (Ola, 2026-09-27); 21d, which it was placed after, is deferred behind the basin work | none yet |
| 15 | A DEM in several tiles, and the domain in its own CRS: 15a Norway, many tiles in one projected CRS (Ola's 254 DTM10 UTM33 tiles; only the selected tiles must share a lattice, since 8 of them sit half a cell off); 15b the domain polygon reprojected into the DEM's CRS; 15c and 15d the São Francisco basin (geographic DEMs meshed in the DEM's own lattice frame, not resampled; window decoding and memory) | **15a shipped with branch `increment15-dem-mosaic`'s PR**: `--dem DIR` or several files, `--bbox`, tiles selected and stitched on one lattice, Q2-Q5 refusals; overlaps that disagree (real DTM10 tiles exported on different dates, up to 52 m) are split down the middle and reported per seam (Ola's Q1 revised, 2026-09-28). The design's acceptance box, 9 tiles and 10,051² nodes, meshes: 11.05 M triangles in 31 s on battery. **15b built with branch `increment15b-domain-crs`'s PR**: `--domain` in its own CRS (e.g. EPSG:4326) over a UTM33 mosaic, reprojected vertex by vertex. Then 16b (Ola, 2026-09-27: Norway first) | `docs/increments/15-dem-mosaic.md` |
| 16b | Interior polygons and polylines as constraints ("terrain polygons": lakes, land cover, roads, rivers): `--features PATH`, a GeoJSON `FeatureCollection`, each feature naming a vocabulary property; closed or open `Breakline`s, crossings noded, off-node vertices with bilinear z. **Its working example is real data** (Ola, 2026-09-27): CORINE Land Cover 2018 over the benchmark tile `7908_3_10m_z33.tif`. The source is Ola's local copy, `rasputin_data/corine_sql/.../U2018_CLC2018_V2020_20u1.gpkg` (8.2 GB, EPSG:3035, a sibling of this repository, which also holds the 254-tile DTM10 archive for gap 6). It reads without GDAL: sqlite3 over its R-tree, the GeoPackage blob header stripped, `shapely.wkb`, then `pyproj` to 25833, so CRS stops in Python as before. Probed 2026-09-27 in 0.2 s: 60 polygons in 8 classes (heath, bare rock, sparse vegetation, bogs, intertidal flats, water, sea, urban), 11 068 vertices clipped to the tile, median segment 54 m against 10 m cells. The EEA's public ArcGIS service (`image.discomap.eea.europa.eu`, `Corine/CLC2018_WM`) returns the same 11 068 clipped vertices and is the route for anyone without the file. What it forces on 16b's design: neighbouring polygons share their boundaries, so each shared edge arrives twice; the polygons run past the domain and must be clipped; the extract is committed as a fixture with the Copernicus attribution. Placed before 20c because 20c may split constraint segments and should be designed and measured on inputs that have interior ones | **designed** (`@architect`, 2026-09-28): two PRs, 16b-0 (the noder's verifier indexed: it is quadratic, 25 s of a 30 s CORINE run) and 16b-1/2 (GeoJSON and GeoPackage features, class maps, clip as linework, the noder merges shared edges); Q1-Q6 ruled by Ola 2026-09-28 (Q6: the legacy CORINE GML stays and is read). **16b-0 shipped (#107, merged 2026-09-28)**: the verifier's pair search by sort and sweep, always on; `node` on a 48 km CORINE square 24.8 s → 0.067 s. **16b-1/2 built on branch `increment16b12-features`** (2026-09-28): GeoPackage, GeoJSON and the legacy GML (strict reader, fixture repaired) as `--features`; CORINE class maps; the pre-clip keeps whole edges and widens its margin per edge for long ones; the noder merges shared edges; net +642 production lines. Before 20c | `docs/increments/16b-terrain-polygons.md`, `docs/increments/16-domain-polygon.md` (R6) |
| 16c | A land-cover label per triangle: which input polygon (CORINE class) each triangle lies in, by a flood fill over the unconstrained edges. The bits on an edge say what kind of line it is; which class lies on each side is 16c's | planned: ruled by Ola 2026-09-28 (16b's Q2) as a separate increment, next after 16b | `docs/increments/16b-terrain-polygons.md` (Q2, R6) |
| 16d | A record of the input polylines in the mesh. Today a constraint edge keeps only its type bits; the source feature, its attributes (a river's name or order, a road's class) and which edges belong to one polyline are lost, and where edges from two features merge only the union of their bits survives. Proposed: a feature table (one row per input feature: source ID and attributes) and, per constraint edge, the list of features it came from (a list, since a merged edge belongs to several). The noder already tracks each piece's input chain; the open part is the output: a per-edge list in `.vtk`/`.ply`, or a side file | proposed by Ola 2026-09-28; to be designed by `@architect` once Ola has said what it is for; after 16b | - |
| 20c | Soft quality criterion: a penalty that each Steiner node or constraint split must pay for in angle gained, instead of 20's hard 25°; applied at the start and during DEM refinement; may split constraint segments when that improves the mesh. Ola's rulings on 20's C1-C3 | to design after 16b (`@architect` measures cost against 20 first) | `docs/increments/20-start-quality.md` (Ola's rulings) |
| — | Auto-catchment: the watershed upstream of a coordinate, computed from the DEM and handed to `--domain`, so a catchment no longer has to be supplied as a file (Ola, 2026-09-27: "not far into the future"). The textbook route is depression handling (Priority-Flood, Barnes, Lehman and Mulla 2014), D8 flow directions (O'Callaghan and Mark 1984) and accumulation, the pour point snapped to the strongest flow nearby, the upstream cells traced and their outline turned into a polygon; the literature check is `@architect`'s. Open for its design: whether it runs in the C++ core (a 10 m tile is 25 M cells); how a stair-stepped cell outline becomes a domain polygon, which meets input coarsening; **a requirement, not an option (Ola, 2026-09-29): "we will get an extreme amount of points in the catchment polygon. This must be taken into account, or the number of triangles will be overwhelming."** A cell outline has a vertex at every cell step (tens of thousands for a mid-sized catchment at 10 m, far more for Glomma), and every boundary edge is a constraint the mesh must honour; 16b measured what dense constraints cost (CORINE borders make the 10 m mesh 3.6x larger). The design must state how the outline is reduced, to what horizontal tolerance, and what it guarantees (a simple polygon, the pour point inside); the fine outline is defined on the DEM's node lattice (Ola, 2026-09-29: "Isn't the DEM points, really?"), but the reduced polygon's vertices need not be DEM nodes: off-node domain vertices are supported since 16 (Ola: "We have aleady established that the polygons and DEM don't need to match"), so vertices may be placed to keep the area; and the reduced polygon keeps approximately the fine one's area (Ola, 2026-09-29: "the resulting polygon has approximately the same area as the fine one"), since catchment area drives runoff volume. Literature to check: area-preserving polyline simplification (e.g. Bose et al. 2006; Kronenfeld et al. 2020, segment collapse); and that a real catchment crosses tile edges, so it needs gap 6 (a DEM in several tiles) first. Legacy has nothing on it (`grep -rliE "watershed|flow.?acc|flow.?dir|pour.?point|catchment" legacy` returns no files) | **built as increment 22** (2026-09-29, overnight), in two PRs: PR 1, the fine catchment, net 619 lines; PR 2, the outline reduction, net 406 lines. Bygdin: 304.91 km² against NVE's 305.54 km² (−0.21 %), 99.1/99.3 % node overlap with NVE's polygon, 17 812 → 740 outline vertices at 20 m with the area kept, catchment in 6.5 s (`docs/benchmarks/2026-09-29/bygdin.md`). The reduction bounds the distance to the fine outline at its vertices only (measured on Bygdin: 19.92 m at 20 m). Parked for later (Ola, 2026-09-29): "It could be useful to have the option of some tolerance, doing more adaptive coarsening, but write it down and leave it for now": an enforced whole-outline tolerance as an option, which also belongs with 20c's input coarsening. Designed by `@architect`, 2026-09-29, Bygdin first, moved ahead of 20c for Ola's night run: a lake polygon (CORINE) as the seed, a Priority-Flood labelling the lake's catchment, a marching-squares outline, area-preserving segment collapse to `--outline-tolerance` (default twice the cell). Two PRs on one branch; defaults for Ola to confirm are marked in the file | `docs/increments/22-auto-catchment.md` |
| — | `raster/`: grid-to-world geometry and bilinear sampling | shipped (`7785fea`), **no record** | none — predates the protocol |

## What stands between here and an operational MVP

An MVP is one command turning a DEM and a catchment polygon into a terrain TIN
file. Six things were missing, in dependency order. Gap 4 has shipped; gaps 1
and 3 and the first half of 5 shipped as increment 12 (#91); gap 2's first half
shipped as increment 14 (#92), and 14b, which fixes its triangle shape, shipped
as #93.

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
   unchanged. **Designed as increment 15** (`docs/increments/15-dem-mosaic.md`):
   15a and 15b are Norway (many tiles in one projected CRS; the domain in its
   own CRS), 15c and 15d the São Francisco basin (geographic DEMs; memory).

**Order of work, ruled by Ola on 2026-09-27** ("My priorities are to get the
Norwegian cases sorted first"): 15a and 15b (Norwegian multi-tile DEM), then
16b (CORINE terrain polygons over Norway), then auto-catchment. The basin's
15c and 15d follow; 21d and domain decomposition come after.

Open and not MVP-blocking: **inputs in their own CRS**. The user (2026-09-26):
"I don't think the domain CRS should have to match the DEM CRS in the future.
This must be written down. The questions will be in which coordinate system we
shall do the math." **Settled for the domain and a projected DEM by increment
15b**: the domain comes in its own CRS and is reprojected into the DEM's, which
is the computation CRS. Feature geometry (16b) follows the same path; a
geographic DEM's computation frame is 15c (`docs/increments/15-dem-mosaic.md`). Also **a general point insertion policy**, which the user
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

Open, for later (Ola, 2026-09-28): **the surface model, DOM10, alongside the
terrain model.** Kartverket's DOM10 is the top surface (tree crowns, roofs,
bridges) where DTM10 is the bare ground; both are published on hoydedata.no
and, where laser coverage exists, derive from the same NDH laser data (DTM10
falls back to the 2013 contour model elsewhere). Two uses Ola wants to keep open: **shading** (terrain and surface
shadowing, e.g. for solar radiation), which needs the surface, not the ground;
and **canopy and building height, DOM10 minus DTM10**, for vegetation-related
work and land-cover classification. Hydrology keeps meshing the DTM: a DOM
mesh would dam rivers at bridges. Not designed; no increment yet.

Known defect: increment 14's NoData carving from a NoData corner appears to
advance one node per round. Meshing the real tile from its outline alone took
5 054 rounds and 109 s against 0.35 s from the stride grid
(`docs/increments/16-domain-polygon.md`, M4 and R5). It needs its own look.

Clipping to the catchment polygon left this list on 2026-09-26: it is increment
16, pulled forward at the user's request. Interior polygons and polylines follow
as 16b, designed in the same record.
