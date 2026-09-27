# Increment 15 — a DEM in many tiles, inputs in their own CRS, and the computation frame

Status: **Q1-Q5 ruled by Ola; 15a implemented** on branch
`increment15-dem-mosaic` (red `2696bc2`, green `ff7cc8d`), in review. 15b-15d
are designed, not implemented; Q6-Q10 are open. Written by `@architect` before
`@tester`, per `docs/increments/README.md` step 1.

## Ruled by Ola

- **2026-09-27, priority:** "My priorities are to get the Norwegian cases
  sorted first, so let's finish this and turn to the main topics." So the
  increment is split so that **the Norwegian case ships first and on its own**
  (15a, 15b): many aligned tiles in one projected CRS, as in Ola's DTM10
  archive, and a domain or features in their own CRS. The geographic DEM and
  the basin-wide computation frame for the São Francisco basin (15c, 15d) are
  designed here but come **after** the Norwegian sub-increments. Nothing in
  15a or 15b depends on them.
- **2026-09-27, Q1-Q5: "yes to all five, carry on"**, to the recommendations
  in "Questions for Ola":
  - **Q1:** two tiles that disagree at an overlapping node are refused,
    naming both tiles, the node count and the largest difference.
  - **Q2:** `--bbox` stays.
  - **Q3:** the repository lives in `io/repository.py`, the one module in
    `io/` that opens files.
  - **Q4:** tiles with different NoData sentinels are refused.
  - **Q5:** a request mixing the eight half-cell tiles with the main lattice
    is refused; a request inside either lattice is meshed.

  Q6-Q10 (the basin) wait until after the Norwegian sub-increments.
- **2026-09-27, Q5 as read after review ("Yes", to the main session's
  proposal):** @reviewer found that DTM10's 51-node overlaps make a box that
  only reaches into a half-cell tile's overlap strip select that tile and be
  refused as mixed-lattice, although the main lattice covers every node (the
  design's own 15a acceptance box was refused; 98 of 576 20-km boxes around
  the eight tiles). The rule is now, per lattice, whether **its own tiles
  cover every node the request needs**:
  - exactly one lattice covers it: mesh on that lattice, and drop the other
    lattices' tiles from the plan;
  - none covers it: refuse, as before (the Q5 message);
  - several cover it (a box wholly inside an overlap strip): use the lattice
    with **the most tiles in the repository**, ties broken by the name of its
    first tile. Header-only and order-independent.

## Scope

**15a and 15b (Norway, first).** Close `ROADMAP.md` gap 6 (a DEM in several
tiles) for tiles in one projected CRS, and "inputs in their own CRS" for a
projected DEM. After them,

```sh
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20220924/ \
    --domain catchment_wgs84.geojson --tolerance 1 --out catchment.vtk
```

selects the tiles the catchment needs out of 254, stitches them, transforms the
catchment from EPSG:4326 into the DEM's EPSG:25833, and meshes it. `--dem
file.tif` alone behaves exactly as today.

**15c and 15d (the basin, after Norway).** Geographic DEMs (ANADEM, EPSG:4326),
the computation frame, `--out-crs`, and decoding only the needed part of a
2.5 GiB tile. They settle the question Ola put on 2026-09-26, "The questions
will be in which coordinate system we shall do the math"
(`16-domain-polygon.md`, "Ruled by the user").

**Not in this increment.** Domain decomposition (increment 21's option D),
a block-sparse raster in C++, tiles in more than one CRS, and resampling. Each
is named in "Not in scope" with the reason it waits.

## What this builds on

The parked design, `git show f8cedbe:docs/increments/15-dem-repository.md`
(branch `increment15-dem-repository`, never merged), called "the parked design"
below: rulings R1-R8, questions U1-U4. It was written for Ola's DTM10 archive
before anyone had read the archive's headers. They have now been read (N1-N5),
and ANADEM's too (B1-B9).

What is carried across, and what changes:

| parked | here | why |
|---|---|---|
| R1: repository protocol, pure `plan_mosaic`, `assemble` | carried (R1) | still right; a test double is a dict |
| R2: paths only in `io/repository.py` | carried (R2), still a question (Q3) | plus a new `dem_input.py`, so `cli.py` (1101 lines) does not grow |
| R3: `read_meta` split out of `decode_dem` | carried (R3); window decoding added in 15d | an ANADEM tile is 2.5 GiB |
| R4: every tile in the directory on one lattice, or refuse | **changed** (R4): only the **selected** tiles must share a lattice | 8 of Ola's 254 tiles are half a cell off (N2); the parked rule refuses the whole archive |
| R4: canvas origin from the tiles as decoded | **changed** (R4): each lattice has a **reference node**, and every node has a global index | frame coordinates must not depend on the window, for domain decomposition (R12) |
| R5: valid beats NoData, disagreement refused, order-independent | carried (R5), now **with real data behind it** | DTM10 (N3) and ANADEM (B4) overlaps agree bit for bit |
| R5: gaps stay NaN | **changed** (R5): a node the request needs that no tile covers is refused | Ola's ruling, `18-row-span-scan.md` R6: "a missing tile inside the extent is a data error" |
| R6: `--bbox` in the DEM's CRS, no transform | carried (R6); the domain is transformed in 15b (R9) | "inputs in their own CRS" |
| R6: `GeoPolygon` lands with the clip | superseded: 16 shipped `DomainPolygon`; 15b gives it a transform | |
| R7: `MAX_MOSAIC_BYTES = 1 GiB`, the canvas copied twice | **changed** (R7): cap at half of physical memory; one canvas, not two | a Norwegian catchment's box passes 1 GiB (R7) |
| R7: window decoding later | 15d (basin) | DTM10 tiles are 100 MB decoded; ANADEM's are 2.5 GiB |
| R8: `--dem` repeatable, `--bbox`, `dem_tiles` | carried (R11) | |
| U1-U4 | carried to the questions (Q1-Q4), U1 now with measurements | |

## What was measured

On `increment15-dem-mosaic` at `59024e3`: tifffile 2026.9.20, pyproj 3.8.0
(PROJ 9.8.1), numpy 2.5.3, Python 3.14.7, on Ola's Mac (32 GiB,
`sysctl hw.memsize` = 34359738368; 10 cores). The scripts are checked in,
so each number can be re-run:

- `docs/increments/15-probes/dtm10_probe.py headers|seams DIR` (N1-N3), run
  on `../rasputin_data/DTM10_UTM33_20220924`;
- `docs/increments/15-probes/anadem_probe.py headers|seams|resample|frame`
  (B2-B7). It reads ANADEM by HTTP range requests, so no tile is downloaded
  whole, and it needs network access.

They are measurement scripts, not production code, and nothing imports them.

### Norway: the DTM10 archive

- **N1. What the archive is.** 254 `.tif` files (each with a `.tfw` and a
  `.aux.xml` beside it), 19 GB. Every file passes `io/geotiff.py`'s own header
  checks (the probe calls them, so a refusal would be `decode_dem`'s):
  EPSG:25833, 10 m in both axes, area-registered, float32, LZW, tiled, NoData
  −32767 from tag 42113. 243 are 5051 × 5051 nodes; the rest are 5052 or 5053
  on a side, and four are short on one side (7305_3 2881 × 5051, 7405_1
  5051 × 3521, 7405_2 5051 × 3511, 7507_4 4103 × 5052; corrected in review:
  the probe's per-axis tallies had been paired into two tiles that do not
  exist). Header-only read of all
  254: 0.56 s wall.
- **N2. 246 tiles share one lattice; 8 do not.** Taking the north-west-most
  node (x −100250, y 7950250) as reference, 246 tiles sit at integer node
  offsets. **Eight are half a cell (5 m) off east-west**, on integer rows:
  7304_1, 7507_4, 7606_2, 7707_1, 7707_3, 7807_2, 7807_3, 7808_3. Their nodes
  are at x = …745 instead of …750, and they are 5052 wide and mostly 5053 high.
  They look like a second production run. Under the parked R4 (whole directory
  on one lattice) `--dem DIR` would refuse the whole archive.
- **N3. Neighbours overlap by 51 nodes, and the overlaps agree.** A 5051-node
  tile is 50 km plus 510 m, so each edge neighbour shares 51 rows or columns
  (a few share 52 or 53; one pair 53). Over 4 neighbour pairs of the first
  tile, every overlapping cell is equal: 257 601 of 257 601 on each edge
  overlap, 2601 of 2601 on each corner. The control shifts the comparison by
  one node: on the two edge overlaps it finds a maximum difference of 18.0 m
  and 46.3 m, so the probe can fail there. On the two 51 × 51 corner overlaps
  the shifted comparison is also all equal: those corners are flat (sea, most
  likely), so for them the control is **not** able to fail. That is 4 pairs of
  869 overlapping pairs.
- **N4. Size.** The union of all 254 tiles spans 155 051 × 125 051 nodes:
  19.4 G nodes, 72 GiB as float32. The tiles hold 6.45 G nodes (24 GiB), so
  most of the box is sea or outside Norway. A 2 × 2 block of tiles is
  10 051² = 101 M nodes, 404 MB.
- **N5. The committed benchmark tile is a different release from the
  archive's.** `tests/fixtures/dem_archive/7908_3_10m_z33.tif` and the
  archive's file of the same name have the same tie point and shape, but
  910 706 of 25.5 M cells differ, by up to 46.1 m, and 2868 cells are valid in
  the archive and NoData in the fixture. So a real multi-tile fixture must be
  cut from one release, and mixing the committed tile with archive neighbours
  is a mixing of releases.

### The basin: ANADEM

- **B1. Where ANADEM lives.** `https://hge-iph.github.io/anadem/` lists 53
  tiles, `anadem_v1_<zone><band>.tif` under
  `metadados.snirh.gov.br/files/anadem_v1_tiles/`, named by MGRS grid-zone
  designator (17L … 25M). The server accepts range requests. Tiles are 0.58 to
  2.09 GB compressed (`curl -I`). The GitHub repository's `LICENSE` is MIT and
  speaks of "the Software" (see Q10).
- **B2. What an ANADEM tile is.** Every tile near the basin: EPSG:4326
  (`GTModelTypeGeoKey` 2, `GeographicTypeGeoKey` 4326), **area-registered**,
  float32, Deflate, tiled 512 × 512, NoData −9999 in tag 42113, and 4 to 6
  extra pages, all reduced-resolution overviews. Spacing
  **0.00026949458523585647° in both axes**: 0.970″, not 1″. It is 30 m of
  equator in degrees (`d·π·a/180 = 30.000000000000004`), so ANADEM is **not on
  Copernicus's 1″ lattice**. A tile's header costs 64 KiB of reads.
- **B3. The tiles are huge and share one lattice.** They are 6° grid zones,
  not 1° tiles: up to 674 M cells, 2.5 GiB as float32. The basin needs about
  six, not the research note's "~150". Every tile's offset from 23L is an
  integer number of cells to within 2.2e-11 cells, with identical spacing.
  Neighbours overlap by 8 or 9 cells, and 24L/24K by 87 rows. Band edges sit
  near 8.13°S and 16.26°S, not on MGRS's 8° and 16°.
- **B4. Overlapping cells agree bit for bit.** 23L|24L (the 42°W seam, 4096
  cells) and 23L|23K (the 16.26°S seam, 4088 cells): all equal. The control can
  fail: shifted by one column, 34 of 3584 equal and a maximum difference of
  19.9 m; shifted by one row, 129 of 3577 and 3.4 m. In the 87-row 24L|24K
  overlap, 40 776 of 44 544 cells are valid in one tile only, and the 3768
  valid in both are equal. So **valid-beats-NoData is needed on real data**.
  Three samples of three seams, not the archive.
- **B5. Size.** ANADEM's cell at 7°S / 14°S / 21°S is 29.78 / 29.12 / 28.02 m
  east-west and 29.80 / 29.82 / 29.84 m north-south. The basin (636 920 km²,
  OAS) is **734 M nodes**. Its bounding box (48–36°W, 21–7°S, the research
  note's approximation) is **2.31 G nodes, 8.6 GiB as float32**. The basin is
  mostly in 23K, 23L, 24K and 24L, with slivers of 23M and 24M. Whether any of
  it lies west of 48°W (22K, 22L) is to be read from the BHO polygon; not
  checked.
- **B6. What resampling costs.** A real 1536² window of 23K over the Serra do
  Espinhaço (43.86–43.45°W, 19.16–18.74°S, heights 592–1682 m, no NoData) was
  resampled bilinearly onto square grids in the research note's basin LCC. The
  resampled surface was then compared with the source at every inner source
  node:

  | target spacing | max | p99.9 | p99 | median |
  |---:|---:|---:|---:|---:|
  | 30 m | 13.38 m | 5.43 m | 3.07 m | 0.41 m |
  | 20 m | 14.06 m | 3.77 m | 2.13 m | 0.28 m |
  | 10 m | 8.91 m | 1.96 m | 1.11 m | 0.14 m |

  Control: resampling onto the source's own nodes gives 9.7e-10 m. So a mesh
  built to `--tolerance 1` on a resampled grid can be 13 m off the DEM it came
  from, and a finer target grid does not remove that.
- **B7. How far a straight edge bends between frames.** Take an edge that is
  straight in longitude and latitude, and write it in the basin LCC. Its
  midpoint lies off the straight LCC edge by 5 mm at 1 km, 0.14 m at 5 km,
  2.2 m at 20 km and 13.6 m at 50 km (worst direction, at 14°S; within 10 % of
  that at 7°S and 21°S). The error grows with the square of the edge length.
- **B8. What refuses an ANADEM tile today.** `io/geotiff.py`'s private checks
  on 23L's header: `_single_page` passes (the extra pages are reduced),
  `_check_page` gives float32, `_placement` gives the node grid, `_nodata`
  gives (−9999, "tag"), and `_projected_epsg` refuses: "GeographicTypeGeoKey
  (2048) = 4326 is a Geographic 2D CRS; a projected CRS is required". Nothing
  else stands in the way.
- **B9. The core's indices hold at basin size, by grep.**
  `include/terrain/raster/` indexes with `std::size_t` (`geometry.hpp:60`), and
  `RowSpan` carries `uint32` per axis (`mesh/row_spans.hpp:31`). A grep for
  32-bit flat node indices in `include/terrain/` found none. A run on 2.3 G
  nodes is @perf's (15d acceptance).

Not measured: Deflate or LZW decode throughput for a whole tile, the triangle
count at basin scale. The BHO polygon's CRS was checked afterwards by the main
session (web, 2026-09-27): ANA's metadata gives SIRGAS 2000, **EPSG:4674**,
geographic. ANADEM's *data* licence is still not found: no statement beyond the
repository's MIT software licence turned up, so Copernicus GLO-30's terms for
modified data apply at least (its notice, plus a notice of modification).

## Prior art: legacy and literature

### Literature

- **The refinement method is unchanged**: greedy insertion, Garland and
  Heckbert, "Fast polygonal approximation of terrains and height fields",
  CMU-CS-95-181, 1995 (increments 14 and 14b). This increment changes which
  grid refinement runs on and in which frame, not the method.
- **Barycentric coordinates are invariant under affine maps** (textbook; e.g.
  Farin, *Curves and Surfaces for CAGD*, 5th ed., 2002; recalled, not reread).
  R8 rests on it. At a given point, the vertical error of a linear triangle is
  the same in any affine image of the plane. So a sup-norm guarantee computed
  in an affine frame of the DEM's lattice is exactly the guarantee in the DEM's
  own coordinates.
- **Bilinear interpolation error** is O(h²·|f''|) for smooth f (e.g. Ciarlet,
  *The Finite Element Method for Elliptic Problems*, 1978; recalled). Terrain
  is not smooth at 30 m, so the bound says little. B6 measures instead.
- **DEM resampling and reprojection change terrain.** Usery, Finn, Scheidt,
  Ruhl, Beard and Bearden, "Geospatial data resampling and resolution effects
  on watershed modeling", *J. Geographical Systems* 6:289-306, 2004; Kienzle,
  "The effect of DEM raster resolution on first order, second order and
  compound terrain derivatives", *Transactions in GIS* 8(1):83-111, 2004. Both
  recalled, not reread. They study derived attributes (slope, flow). Our
  concern is narrower: resampling changes the surface the sup-norm is measured
  against. **Difference:** we do not resample (R8).
- **Mosaicking and warping practice, for reference only (no GDAL).**
  `gdalbuildvrt` composites overlapping sources in order, later over earlier.
  `gdalwarp` maps every destination pixel back through the inverse transform,
  approximated by linear interpolation along rows within an error threshold
  (default 0.125 pixel). GDAL documentation, recalled. **Differences:** no order
  precedence, since disagreement is refused and valid beats NoData in any
  order (R5); and no warp (R8).
- **Projection choice**: Snyder, *Map Projections — A Working Manual*, USGS
  Professional Paper 1395, 1987, recalled. A Lambert conformal conic suits an
  extent wide east-west at mid latitudes, with the standard parallels about one
  sixth of the latitude range in from each edge. The research note's LCC (10°S
  and 18.5°S over 7–21°S) is close to that rule (9.3°S, 18.7°S). Here it is
  only a suggested **output** CRS (R10).
- **ANADEM**: Laipelt et al., "ANADEM: A Digital Terrain Model for South
  America", *Remote Sensing* 16(13):2321, 2024. The grid facts used here come
  from the tiles' headers (B2), not from the paper.
- **Domain decomposition** prior art (tile-parallel terrain simplification,
  *Remote Sensing* 12(3):437, 2020; Linardakis and Chrisochoides 2006) is cited
  in `21-parallel-refine.md`. Not designed here (R12).

**Novelty: none claimed.** One claim could look tempting later: an exact
sup-norm guarantee against a geographic DEM's own nodes, with the mesh written
in a projected CRS. It is not made. Before anyone makes it, search Google
Scholar and IEEE Xplore for "TIN generation geographic coordinates DEM
reprojection error", "terrain simplification latitude longitude grid
projection" and "greedy insertion DEM geographic lattice". No web search tool
was available to this round; ANADEM's facts were read with `curl` and range
requests.

### Legacy

```sh
$ grep -rlE 'RasterRepository|get_intersections|add_raster|mosaic' legacy/
legacy/bindings.cpp
legacy/rasputin/mesh.py
legacy/rasputin/application.py
legacy/rasputin/reader.py
legacy/rasputin/web_visualize.py
legacy/tests/test_gml_repository.py
legacy/tests/test_land_cover_repository.py
legacy/tests/test_raster_repository.py
$ grep -rliE 'reproject|resampl|Transformer|to_crs|transform\(' legacy/
legacy/rasputin/globcov_repository.py
legacy/rasputin/reader.py
legacy/rasputin/application.py
legacy/rasputin/avalanche.py
legacy/rasputin/geometry.py
legacy/rasputin/gml_repository.py
legacy/tests/test_mesh.py
legacy/tests/test_land_cover_repository.py
legacy/tests/test_gml_repository.py
```

- **The repository** (`legacy/rasputin/reader.py:429-468`) was reported on in
  `11-raster-ingestion-prior-art.md` §3.4, §4.12, §4.13, §5.8 and §8.4, and
  the parked design ruled on it. That ruling is carried unchanged. Carried:
  the repository over a directory, header-only footprints, selection of the
  covering tiles. Not carried: `glob` order, the unchecked CRS, polygon
  subtraction, the `1e-10` stop area, `not touches`.
- **The frame, new here.** `legacy/rasputin/application.py:101-123` transformed
  the domain into the raster's CRS, meshed there, and transformed the finished
  mesh's vertices into a caller-chosen `target_coordinate_system`. That is the
  shape of R8-R10: compute in the DEM's own frame, bring vectors in, send the
  result out. **Carried.** Not carried: every CRS was spelled `+init=` and no
  `Transformer` passed `always_xy` (increment 11 §9); z went through
  `proj.transform(x, y, z)`, which is harmless for a 2D CRS but says nothing;
  the bending of edges (B7) was never measured or recorded; and the output CRS
  was a proj4 string.
- `legacy/rasputin/geometry.py:235-241`, `GeoPolygon.transform`: the domain's
  transform, with the same `+init=` caveat. R9 carries the operation into
  `DomainPolygon` through the one reprojection helper.
- `globcov_repository.py`, `gml_repository.py` and `avalanche.py` transform
  land cover and avalanche inputs point by point. That is 16b's concern.

A new `@migration-expert` pass is not needed: no legacy constant is re-derived.

## Rulings

Sub-increment tags say where each piece lands: **[15a]** and **[15b]** are
Norway, **[15c]** and **[15d]** the basin.

### R1. The interface: a narrow repository, pure planning, assembly beside it [15a]

Carried from the parked R1 with two additions, the lattice key and the
reference node (R4).

```
DemRepository (Protocol)               storage only
  footprints() -> tuple[TileFootprint, ...]       header-only, sorted by name
  load(name: str) -> DemTile                      whole tile, decoded
                                                  [15d: load(name, window)]

plan_mosaic(footprints, bounds, needed=None) -> MosaicPlan    pure, no pixels
assemble(plan, load) -> Mosaic                                pixels, no files
```

- `TileFootprint`: `name: str`, `meta: RasterMeta`.
- `MosaicPlan` (frozen): the mosaic's `RasterMeta`; the lattice's reference
  node and the window's global index offset `(row0, col0)` (R4); per selected
  tile, its name, its expected `RasterMeta`, and two index windows (where it
  sits in the canvas, and which part of it is used). Every refusal that
  headers can decide fires in `plan_mosaic`, before any pixel is read.
- `Mosaic`: the assembled `DemTile` plus the plan, as provenance.
- The repository does not select or assemble. Those steps are the same grid
  arithmetic for every storage backend. The protocol stays at two methods, and
  sync: `load` is blocking I/O, and an async caller wraps it in
  `asyncio.to_thread`.
- **A declarative front, `dem_input.py`.** `DemRequest` (frozen Pydantic:
  `sources: tuple[Path, ...]`, `bounds: Bounds | None`, `nodata: float | None`;
  15b adds the domain) goes in, and `open_dem(request) -> DemInput` returns the
  tile, the plan and the provenance strings. `cli.py` parses flags into a
  `DemRequest` and calls `open_dem`. A GUI backend or an API worker builds the
  same request without Typer. `cli.py` is 1101 lines today, and this keeps the
  mosaic out of it.

### R2. Homes, and where paths live [15a]

| module | holds | imports |
|---|---|---|
| `io/geotiff.py` | `read_meta(stream)`, the header phase (R3) | as today |
| `io/repository.py` | `DemRepository`, `TileFootprint`, `TiffDemRepository` | `io.geotiff`, `io.models` |
| `tin_engine/mosaic.py` | `Bounds`, `IndexWindow`, `MosaicPlan`, `Mosaic`, `MosaicError`, `plan_mosaic`, `assemble` | `io.models`, numpy, shapely |
| `tin_engine/dem_input.py` | `DemRequest`, `DemInput`, `open_dem` | the three above |
| [15b] `tin_engine/crs.py` | `parse_crs`, `reprojector`: the one `Transformer.from_crs` site | pyproj |
| [15c] `tin_engine/frame.py` | `LatticeFrame`, `frame_for` | `io.models`, `crs` |

- Paths live in `io/repository.py` and `dem_input.py` (which only builds the
  repository), and nowhere else below `cli.py`. `read_meta` and `decode_dem`
  still take streams. Where the repository lives is still the parked U3, now
  Q3.
- `TiffDemRepository.from_directory` lists `*.tif` and `*.tiff`
  (case-insensitive), not recursively. That skips DTM10's `.tfw` and
  `.aux.xml` side files, which carry nothing the GeoTIFF tags do not already
  hold. Paths are resolved, deduplicated and sorted by their string form.
- A file in the listing that fails `read_meta` is refused, naming the file,
  and is not skipped.
- Footprints are read once, lazily, and cached. `assemble` checks that each
  loaded tile's `meta` equals the planned one ("changed since it was listed").

### R3. The header phase [15a], and window decoding [15d]

- [15a] Carried from the parked R3: `_header(tif, nodata) -> (RasterMeta,
  dtype)` holds everything `decode_dem` does before `page.asarray()`
  (now `_header` itself, `io/geotiff.py:105-140`). `read_meta(source, *, nodata=None)` calls it (through
  `read_header`, which also returns the decoded dtype, S2) and returns, and `decode_dem` calls it and then decodes. The header phase keeps
  the codec refusal.
- [15d] `decode_dem(source, *, window=None)`: with an `IndexWindow`, decode
  only the TIFF blocks (or strips) that meet it, and return a `DemTile` of the
  window, with its own `x_min`/`y_max`. Both archives are tiled 512 × 512 (N1,
  B2), so the blocks come from `page.dataoffsets` and `page.decode`, decoded on
  a thread pool: zlib and LZW release the GIL. A window's `DemTile` must equal
  the whole tile's `DemTile` sliced to the window, `meta` included (T-window).

### R4. Planning: the selected tiles on one lattice, or refuse [15a]

**What changes from the parked R4.** The parked design required every tile in
the repository to share one lattice. Ola's archive has two lattices (N2), so
that rule refuses the whole archive for a catchment that touches none of the
eight odd tiles. Here:

1. **Footprints are grouped into lattices.** Two tiles are on the same lattice
   when they have the same `epsg`, `delta_x`, `delta_y`, `pixel_is_area` and
   `nodata`, and their node offsets are integers within `ALIGN_TOLERANCE = 1e-6`
   cells. The measured noise is at most 2.2e-11 cells (B3), and exactly 0 for
   DTM10. Grouping looks at headers only and is cheap: 254 headers take 0.56 s.
2. **Each lattice has a reference node**: its north-west-most node, the minimum
   node `x_min` and maximum node `y_max` over its tiles, each taken as decoded.
   A node's **global index** `(R, K)` counts from it. For DTM10's main lattice
   that is (x −100250, y 7950250). The reference depends only on which tiles
   the repository holds. It does not depend on the request or on file order.
3. **Selection** is by the request's box in the DEM's CRS. A tile is selected
   when its node rectangle meets the box, closed on both sides (corrected by
   the red suite: when its node *index range* meets the snapped window; see
   "Pinned by the red suite (15a)"). **All selected
   tiles must be on one lattice** (`refuses_mixed_lattice`). The message names
   a tile from each lattice, and the difference: "7707_1 is 0.5 cell (5 m)
   east-west off 7707_2's lattice". Different CRS, spacing, registration or
   NoData are named as such (the parked `refuses_mixed_crs_mosaic`,
   `refuses_mixed_spacing`, `refuses_mixed_registration`,
   `refuses_mixed_nodata`).
4. **The window.** `c0 = floor((bx_min − X_ref) / dx)`,
   `c1 = ceil((bx_max − X_ref) / dx)`, rows likewise (an edge within 1e-6 cell
   of a node is first snapped onto it, so float noise in a spacing like 0.1 m
   cannot add a node line), clamped to the lattice's
   union. The mosaic's node grid is that index window, snapped outward, so it
   covers the box. From here on everything is integer index arithmetic. The
   mosaic's `x_min` is computed as `X_ref + c0 · dx`, and its `y_max` likewise.
   For DTM10 (integers times 10) that is exact. The parked design's
   measurement 4 (split and re-stitch reproduces the origin bit for bit) is
   pinned again as M9.
5. **Coverage** (changed from the parked design's "gaps stay NaN", after Ola's
   ruling in 18 R6). A node the request **needs** that no selected tile covers
   is refused (`refuses_uncovered_nodes`), naming the uncovered area in the
   DEM's CRS. The needed region is:
   - with `bounds` alone, the whole window;
   - with a domain [15b], the domain polygon grown by one cell. Only it is read
     (bilinear z of a boundary vertex reads the four nodes around it; 16 R2).
     Canvas cells outside the needed region that no tile covers stay NaN, as
     filler, and are never read.

   The check is shapely's `covers(union of tile node rectangles, needed)`, in
   the DEM's CRS. It is exact enough for rectangles on one lattice: every
   corner is a node coordinate. (Corrected by the red suite: the check is per
   node, since abutting area-registered tiles leave a gap between their node
   rectangles; see "Pinned by the red suite (15a)".)
6. **Refused as well:** an empty repository; a box that meets no tile; a window
   under 2 × 2 nodes; a canvas over the memory cap (R7).

`plan_mosaic` never sees a path, a stream or a pixel. Its test double is a list
of `TileFootprint`.

### R5. Assembly: stitching and overlaps [15a]

Carried from the parked R5 almost unchanged.

- The canvas is allocated once, NaN-filled, in `np.result_type` of the
  selected tiles' dtypes (float32 for both archives).
- Tiles are loaded one at a time in plan order. Each tile's used window is
  merged into the canvas, and the tile is dropped.
- **Overlaps**, per node:
  - one value NoData (NaN or the sentinel), the other valid: **valid wins**;
  - both NoData: NoData;
  - both valid and equal (`==`): accepted;
  - both valid and different: **refused** (`refuses_overlap_disagreement`),
    naming both tiles, the number of disagreeing nodes and the largest
    difference. That is the parked U1 (a), now Q1, with measurements: DTM10
    (N3) and ANADEM (B4) overlaps agree bit for bit where sampled.
- The result does not depend on tile order (M10).
- **One tile, no bounds: the loaded `DemTile` is returned as it is.** No
  canvas, no copy, so `--dem file.tif` stays bit-identical to today, memory
  included.
- The mosaic is one regular node grid, so `raster.to_core`, `grid_domain`,
  `refine` and `elevation.trim` run unchanged. Refinement's lattice rule
  (increment 14 R1-R2) holds.

### R6. Which region to mesh [15a]

The union of the selected lattice by default. `--bbox XMIN YMIN XMAX YMAX` in
the DEM's CRS, as in the parked R6 (Q2). With a domain [15b], the domain's
bounds in the DEM's CRS take `--bbox`'s place, and the two flags exclude each
other (16 R1 already said so).

Without `--bbox` or `--domain`, `--dem DIR` on Ola's archive selects every
tile, the half-cell ones included, so the mixed-lattice refusal (Q5) fires
before the memory cap (R7) is reached; `test_no_bounds_selects_both_lattices_
and_is_refused` pins that order. Both refusals are correct, and the lattice
message ends "but a --bbox inside one lattice is meshed". On one lattice alone,
the 72 GiB union box would be refused by the cap, whose message says
"narrow it with --bbox".

### R7. Memory

**[15a] One canvas, capped at half of physical memory.**

- **No second canvas copy.** `DemTile`'s validator copies its array
  (`io/models.py:82`), so the parked design peaked at two canvases. `assemble`
  builds the canvas privately, so it can hand it over without a copy: a
  private constructor in `io/models.py`, `_adopt(meta, array)`. It runs the
  same shape and dtype checks, sets the array read-only, and stores it. Its one
  caller is `assemble`, on a buffer that `assemble` allocated and never
  returns writable. The concurrency rule (increment 11 §7) survives: nobody
  else holds a writable reference.
- **The cap** is `physical memory / 2` (`os.sysconf("SC_PHYS_PAGES") *
  os.sysconf("SC_PAGE_SIZE")`, which works on macOS and Linux), checked in
  `plan_mosaic` before any pixel is read. On Ola's Mac that is 16 GiB, about
  4.3 G float32 nodes. The peak is then the canvas plus one decoded tile (100
  MB for DTM10) plus the mesh. The parked 1 GiB (268 M nodes) would refuse a
  large Norwegian catchment's box outright, and the fixed number did not know
  the machine.
- **What fits.** A 2 × 2 DTM10 block: 404 MB. A 164 km × 164 km box at 10 m:
  1 GiB. All of Norway at 10 m (72 GiB box, 24 GiB of tiles): no. That is DD's
  or a block-sparse raster's problem (R12), not this increment's.

**[15d] The basin.** Its box is 8.6 GiB at float32 (B5), under the 16 GiB cap
on Ola's Mac. So **the basin is meshable in one piece with a dense canvas on
a 32 GiB machine**, with two conditions:

- **Window decoding (R3).** Without it, a request touching 23L decodes all
  2.5 GiB of it, and assembly holds the canvas plus a whole tile. With it, the
  transient is one tile's window.
- **What must be held:** the canvas (8.6 GiB), the mesh, and refine's
  per-triangle state. **What is streamed:** the tiles, one window at a time,
  and within a window, blocks on a thread pool. Nothing streams during refine:
  the scan runs on many threads over a read-only raster (14 R7), and 18 R6
  names lazy decode inside the scan as a hazard. Nothing here decodes lazily.
- On a 16 GiB machine the cap is 8 GiB and the basin box is refused. The
  answers to that are a block-sparse raster that holds only the 512² blocks
  meeting the basin (about 2.7 GiB, since the basin is 32 % of its box), or
  domain decomposition. Neither is designed here (Q9).

### R8. The computation frame [15c; basin, after Norway]

**The question.** The core wants numbers in one Cartesian frame, roughly in
metres, with square cells if 21b's integer incircle is to apply (`dx == dy`,
`21-parallel-refine.md` QW2). ANADEM gives degrees, with cells 2.4 % taller
than wide at 14°S and 6.5 % at 21°S (B5). Three routes:

- **(A) Resample onto a square grid in a projected CRS** (the basin LCC, UTM
  23S or Polyconic), in Python, chunked. Everything downstream is unchanged.
  But the tolerance is then measured against the resampled grid, not the DEM:
  on real terrain, 13 m off at the DEM's own nodes (B6). This also breaks 16
  R0 ("a DEM value is the exact height at a point") for the source, and it
  costs a second full-size array plus 2.3 G inverse transforms.
- **(B) Mesh in a projected CRS and sample the geographic grid.** The DEM's
  nodes are then a curved grid in the mesh's frame. Increment 18's row-span
  scan, and 14's lattice rule that every inserted vertex is a DEM node, both
  assume the nodes form a regular grid in the frame. This route rewrites the
  scan. Rejected.
- **(C) The lattice frame. Recommended (Q6).** Mesh in an affine image of the
  DEM's own lattice, bring vectors in by transform (R9), and send the mesh out
  by transform (R10). Nothing is resampled. The legacy did this
  (`legacy/rasputin/application.py:101-123`).

**Why (C) keeps the guarantee exactly.** The frame is an affine map of
(longitude, latitude). Barycentric coordinates are affine-invariant, so at every
DEM node the vertical error of the mesh in the frame equals its error in
longitude and latitude. The sup-norm tolerance holds at the DEM's own nodes
with no loss. What the frame changes is only horizontal: triangle shape, and
which triangulation is Delaunay, are judged in the frame, not on the ground.

**The frame, precisely.** `LatticeFrame` (frozen) holds the DEM's CRS, the
lattice's reference node in it, the DEM spacing `(du, dv)`, and the frame
spacing `(hx, hy)`:

- **For a projected DEM it is the identity**: `hx = du`, `hy = dv`, and the
  frame's origin is the DEM's origin. Everything Norway does is bit-identical
  to having no frame at all.
- **For a geographic DEM:** `hx = round_to(du · π·a/180, 2⁻¹⁰)`, and `hy` from
  `dv` likewise, where `a` is the CRS ellipsoid's semi-major axis. For ANADEM
  that gives `hx = hy = 30.0` exactly (B2). The node with global index
  `(R, K)` has frame coordinates `(K · hx, −R · hy)`, and `RasterMeta` handed
  to `to_core` carries `x_min = c0 · hx`, `y_max = −r0 · hy`, `delta = h`.
- **Exact.** `h` has at most 15 significant bits and `|K|` stays below 2²¹
  (all of South America is 1.3 M columns). So `K · h` is an exact double,
  whatever the window. The same node has the same frame coordinates in every
  mosaic and every subdomain (R12). 21b's condition "significant bits of dx
  plus `bit_width(max(rows, cols) − 1)` ≤ 53" holds: 4 + 16 for the basin.
- **Square when the DEM's cells are square in degrees**, as ANADEM's are
  (`du == dv`). Then `hx == hy`, and 21b's integer incircle applies. The price
  is anisotropy: a frame unit is `cos φ`-dependent on the ground. At 21°S it
  is 0.934 m east-west against 0.995 m north-south, a ratio of 1.065. That
  moves a 45° angle by at most 1.8°. Scale does not matter: the tolerance is
  vertical.
- **Anisotropy limit** `MAX_FRAME_ANISOTROPY = 1.10` (Q8): the ratio of ground
  metres per frame unit, east-west against north-south, over the request's
  extent. Above it, `frame_for` refuses. The basin passes at 1.065. A
  geographic DEM of Norway (60–71°N, ratio 2 to 3) is refused, and should be:
  the triangles would be judged in a frame stretched two- to threefold.
- **A geographic DEM with `du ≠ dv`** gives `hx ≠ hy`. The kernel then takes
  the filtered incircle, which is correct but slower. It is not refused.

**Who chooses the frame: nobody.** It is a function of the DEM's lattice
alone. It does not depend on the domain, on the request, or on an option. That
is deliberate (R12).

**The boundary rule, restated.** `project_structure.md`'s `raster` section
says "every number that crosses is metres in a projected CRS, because nothing
on either side can tell". Its reason was a 2:1 anisotropy nobody could see
(increment 11 §3). In 15c it becomes: **every number that crosses is in the
computation frame, which is Cartesian, affine to the DEM's lattice, and within
`MAX_FRAME_ANISOTROPY` of square on the ground over the request's extent.** A
projected DEM satisfies it as today. `to_core(tile, frame)` refuses a
geographic tile without a frame: that is the gate, and it is tested (T-gate).

**Resampling is designed out, not forgotten.** It is the only route for tiles
in more than one CRS, and for N2's half-cell tiles if they must be meshed
together with their neighbours. If Ola wants either, it is its own increment,
and its guarantee is against the resampled grid, stated as such.

### R9. Inputs in their own CRS

**[15b] For a projected DEM (Norway).** The domain, and later 16b's features,
come in any CRS pyproj can parse. They are transformed into the DEM's CRS
before anything else sees them.

- **One helper**, `crs.reprojector(src, dst) -> Callable[[xy], xy]`, holds the
  only `Transformer.from_crs` in `src_python/`, always with `always_xy=True`.
  Increment 11 §9's grep test already guards every `from_crs` call site; it
  now has one real site to find. A GeoTIFF's model space is (x = easting or
  longitude, y = northing or latitude) whatever the EPSG axis order, so
  `always_xy` is the convention the tiles themselves use.
- `DomainPolygon` keeps its CRS (`epsg: int` becomes a `crs: str`, the text
  pyproj parsed, since a domain may have no EPSG code). A new
  `DomainPolygon.to_crs(dst) -> DomainPolygon` transforms **vertices only**,
  and rings stay straight in the destination. This is Ola's "already discrete"
  direction (16, direction 3): vertices are used as given, and an edge is the
  straight segment between them in the frame the mesh is built in. Within one
  UTM zone the edges bend by millimetres per kilometre of edge (B7's order of
  magnitude). Not densified.
- `check_crs`'s must-match rule (`domain.py:97`) is replaced by the transform.
  It was kept as "one replaceable function at the Python boundary" for exactly
  this (16, "Ruled by the user"). The extent check (16 R1) now runs after the
  transform, against the mosaic's coverage (R4 point 5), not one tile's
  rectangle.
- A GeoJSON file without a `crs` member is EPSG:4326 (RFC 7946), and is now
  transformed rather than refused.
- **Provenance:** the `.vtk` records `domain_crs` (the input's CRS) and
  `domain_transform` (the pyproj transformer's `description`), so a datum
  shift chosen by PROJ is on record.
- Snapping: the noder's 1 mm grid (16 U6) applies after the transform, in the
  DEM's CRS, as today.

**[15c] For a geographic DEM.** Vectors go through the same helper into the
DEM's geographic CRS (longitude, latitude), then through the frame's affine map
into the frame. A domain given in longitude and latitude, like BHO's is
likely to be, has edges straight in the frame, exactly.

### R10. The output CRS [15c; basin, after Norway]

- **`--out-crs TEXT`** (anything `pyproj.CRS.from_user_input` accepts).
  - Projected DEM: optional, and the default is the DEM's CRS. The output is
    unchanged from today, and nothing is transformed.
  - Geographic DEM: **required** (Q7). The refusal message prints a suggested
    LCC fitted to the domain's extent (Snyder's one-sixth rule, central
    meridian at the middle), so the user can copy it.
- **What is transformed:** vertex (x, y) only, frame → DEM CRS → `--out-crs`.
  z is untouched. Triangles and edges keep their indices.
- **What it costs, recorded rather than hidden.** An edge straight in the
  frame is not straight in the output CRS (B7: 0.14 m on a 5 km edge). The
  vertical effect at a point is about the horizontal displacement times the
  triangle's slope. Per triangle, `δ_T = max over its edges of |T(midpoint in
  frame) − midpoint of T(ends)|`, and `ε_T = δ_T · |∇z_T|` in the output CRS.
  That is one extra transform per edge. The `.vtk` records
  `max_reprojection_z_error_estimate` = max over triangles of `ε_T`. It is an
  estimate (midpoints, not a bound), and it says so. When `--out-crs` is the
  DEM's own CRS, the output map is affine, and the field is exactly 0.
- **Recorded CRS:** `EPSG:n` when pyproj finds an exact EPSG code, otherwise
  single-line WKT2 (ASCII, which `write_vtk` requires).
- `--stats` (17) measures angles and areas in the output coordinates, since
  that is the mesh the user gets.

### R11. The CLI and what the file records [15a, 15b]

Carried from the parked R8:

- **`--dem` becomes repeatable**: exactly one directory, or one or more files.
  A directory mixed with files, or two directories, is a usage error.
- **`--bbox XMIN YMIN XMAX YMAX`**, only with `--dem`; excludes `--domain`
  [15b].
- Refusals from `plan_mosaic` and `assemble` (`MosaicError`, a `ValueError`)
  become `typer.BadParameter(param_hint="--dem")`, like `GeoTiffError` today.
- **Fields:** `crs` as today; for a mosaic, `elevation_source` is prefixed
  with `mosaic of N tiles, R x C nodes; `, and `dem_tiles` lists the file names
  (`; `-separated, escaped with `encode("ascii", "backslashreplace")`). [15b]
  adds `domain_crs` and `domain_transform`. [15c] adds `computation_frame`
  (for example `lattice of EPSG:4326, h = 30 m, reference lon -48.00089 lat
  -8.12998`) and `max_reprojection_z_error_estimate`.
- stderr: one line, `mosaic of N tiles, R x C nodes`.

### R12. What domain decomposition needs from this, and what must not block it

Domain decomposition is not designed here (`21-parallel-refine.md` Q5 moved
it to the large-area work). What 15 settles so that it stays open:

- **One frame for the whole DEM.** It is a function of the lattice alone
  (R8), never of the domain, the window or an option. Every subdomain uses the
  same frame.
- **Global node identity.** `(R, K)` from the lattice's reference node (R4).
  Two subdomains that share a seam agree on every seam node's index, and, by
  R8's exactness, on its frame coordinates bit for bit.
- **The raster for a subdomain is a plan.** `plan_mosaic(footprints, bounds)`
  is pure. A subdomain's box gives its own window, under its own memory cap.
- **Overlaps resolve the same way in every subdomain**, because R5 does not
  depend on order.
- **Vectors are transformed once**, before any partition, so a domain vertex on
  a seam has one frame position.
- **The output CRS is explicit** (R10), not fitted per domain. A fitted
  default would put two subdomains, or two catchments, in different CRSs.

Choices rejected partly because they would block it: resampling per window
(neighbouring windows would disagree at seams); a computation CRS fitted to
each domain; order-dependent overlap resolution; a frame origin at each
mosaic's corner (the same node would get different coordinates); decoding
lazily inside the parallel scan.

## Invariants

- **I1. One lattice per mosaic.** Every node of an assembled mosaic lies on
  the selected lattice. Any selected tile off it is a refusal, never resampled
  or snapped.
- **I2. Order independence.** `plan_mosaic` and `assemble` give the same
  result for every permutation of the footprints.
- **I3. Valid beats NoData; disagreement refuses.** No overlapping node is ever
  silently chosen between two different valid values.
- **I4. Coverage.** Every node the request needs comes from a tile. NaN filler
  exists only outside the needed region.
- **I5. Split and re-stitch is the identity.** A grid cut into tiles
  (any overlap, either registration) and planned without bounds reassembles to
  the original `DemTile`, `meta` and array equal with `==`.
- **I6. Header before pixels.** Every refusal that headers can decide fires
  before any `load`.
- **I7. Single-file path unchanged.** `--dem file.tif` output and memory are
  bit-identical to before 15a.
- **I8. [15b] One reprojection site**, with `always_xy=True`. A vector's
  vertices are transformed exactly once.
- **I9. [15c] The frame is a function of the lattice alone**, and exact: a
  node's frame coordinates are `(K · hx, −R · hy)` bit for bit, in every
  window.
- **I10. [15c] No degrees cross into `_core`.** `to_core` refuses a geographic
  tile without a `LatticeFrame`.
- **I11. [15c] The guarantee is unchanged.** With the lattice frame,
  `--tolerance` holds at every valid DEM node inside the domain, measured in
  the DEM's own coordinates.

## Degeneracy policy

- **Overlap width.** Any width, including none (area-registered neighbours
  that abut: disjoint node sets), one shared line (point-registered), 51 nodes
  (DTM10), 87 rows (ANADEM 24L/24K), or a tile wholly inside another. One rule
  covers all (R5).
- **The same file twice** under two names: every overlap is equal, so it is
  accepted. The recorded tile list shows both.
- **A tile touching the box on one node line** is selected, and contributes
  that line.
- **Float noise in tile offsets** up to `ALIGN_TOLERANCE = 1e-6` cells. Measured
  noise is at most 2.2e-11 (B3).
- **Half-cell offsets** (N2): refused by name, not snapped. Snapping would move
  every value of a tile by 5 m, which is resampling.
- **Spacings not exact in binary** (0.1 m, or ANADEM's degrees): the tolerance
  above, then integer windows (parked M3).
- **A domain vertex on a tile seam, or a bilinear cell straddling a seam:**
  nothing special. The canvas is one array.
- **[15b] Axis order.** Every transform is `always_xy`. A GeoJSON with
  latitude first is a wrong file, not a case, and is refused by the extent
  check.
- **[15b] Datum shifts** (SIRGAS 2000 or ETRS89 to WGS 84): whatever
  transformation PROJ picks, recorded in `domain_transform`.
- **[15c] The antimeridian and the poles.** A lattice whose extent crosses
  ±180° is refused. The poles are refused by the anisotropy limit.
- **[15c] `du ≠ dv`** (Copernicus above 50° uses wider longitude spacing):
  allowed, with the filtered incircle. Tiles of different spacing in one
  request are refused as mixed spacing.

## Not in scope

- **Tiles in more than one CRS, and the eight half-cell tiles meshed with
  their neighbours.** Both need resampling, and R8 explains what that costs.
  A catchment inside the eight's own lattice meshes; one straddling the two
  lattices is refused by name.
- **Domain decomposition**, and all of Norway at 10 m in one mesh (R7, R12).
- **A block-sparse raster in C++** (18 R6's `row_segments` source). Needed
  only to mesh the whole basin on a 16 GiB machine; Q9.
- **Recursive directory walks, remote storage, an async repository API.**
- **Vertical datums.** Heights pass through unchanged, as today.

## Sub-increments and LOC

Counted in `CLAUDE.md` §2's unit (production lines, comments and docstrings
excluded, tests excluded). Estimates. Increment 10 came in 39 % over its
estimate and 14 came in 17 % over; the worst case below applies 39 %.

**Norway first (Ola, 2026-09-27).** 15a and 15b ship on their own and in this
order. 15c and 15d follow later.

| | what | est. | worst |
|---|---|---:|---:|
| **15a** | **Norway: many tiles, one projected CRS** | | |
| | `io/geotiff.py`: `_header`, `read_meta` | 20 | |
| | `io/repository.py`: protocol, footprint, listing, cache, load | 60 | |
| | `mosaic.py` types: `Bounds`, `IndexWindow`, placements, plan, `Mosaic`, error | 40 | |
| | `plan_mosaic`: lattice grouping, reference node, selection, coverage, cap | 105 | |
| | `assemble`: canvas, overlap rule, single-tile fast path, meta check | 55 | |
| | `io/models.py`: `_adopt` | 15 | |
| | `dem_input.py`: `DemRequest`, `open_dem` | 35 | |
| | `cli.py`: `--dem` list, `--bbox`, refusals, fields | 45 | |
| | **15a total** | **375** | **520** |
| **15b** | **Norway: the domain in its own CRS** | | |
| | `crs.py`: `parse_crs`, `reprojector` | 30 | |
| | `domain.py`: `crs: str`, `to_crs`, extent against coverage; `check_crs` removed | 40 | |
| | `dem_input.py`: domain to bounds and the needed region | 20 | |
| | `cli.py`: `--domain` with a mosaic, `domain_crs`, `domain_transform` | 25 | |
| | **15b total** | **115** | **160** |
| **15c** | **Basin, after Norway: geographic DEMs and the frame** | | |
| | `io/geotiff.py`: accept geographic 2D, degree axes | 30 | |
| | `io/models.py`: `RasterMeta` CRS kind | 10 | |
| | `frame.py`: `LatticeFrame`, `frame_for`, anisotropy, affine maps | 95 | |
| | `raster.py`: `to_core(tile, frame)`, the gate | 15 | |
| | output: `--out-crs`, transform, bending estimate, WKT2 field, stats on output | 85 | |
| | **15c total** | **235** | **325** |
| **15d** | **Basin, after Norway: window decoding and basin memory** | | |
| | `io/geotiff.py`: window decode, block selection, thread pool | 55 | |
| | `io/repository.py`, `assemble`: windows through `load` | 20 | |
| | **15d total** | **75** | **105** |

Each is under 700 even at its worst case. 15a and 15b together would be 490,
or 680 at the worst case: too close to the ceiling to merge, and Ola asked for
Norway to ship on its own anyway. No sub-increment touches C++, so the
refine/mesh acceptance rule (README, "Acceptance") does not formally apply. The
basin run in 15d's acceptance is @perf's all the same, because it is the first
run at 30× the benchmark.

**Documentation changed in the same PRs** (not counted):
`project_structure.md` (layout; the `raster` boundary rule in 15c, R8);
`11-raster-ingestion.md` §3 and §10 (15c; `GeoPolygon` superseded);
`16-domain-polygon.md` U1 (15b: the must-match rule replaced); `io/__init__.py`'s
docstring (Q3); `ROADMAP.md` rows 15a-15d and the "inputs in their own CRS"
item.

## Tests for `@tester`

Micro-TIFFs come from `tests/python/geotiff_fixtures.py`. Planning and assembly
tests build `RasterMeta` and `DemTile` directly and use a dict as the
repository. The parked design's lists M1-M12, F1-F4, G1-G2 and C1-C6 are
carried as written there, with the changes below.

**Invariant-critical, and mutation testing required:**

- **15a: `test_mosaic.py`** (planning and assembly). The index arithmetic, the
  lattice grouping and the overlap rule are where a wrong answer is a silently
  wrong terrain.
- **15c: `test_frame.py`.** The frame's exactness and the no-degrees gate are
  where a wrong answer is a silently wrong guarantee.

15b and 15d are ordinary suites. 15b's axis-order cases (below) are the ones
that matter there.

**15a, changed or new against the parked list:**

- M2 becomes lattice grouping: two lattices in one repository are accepted;
  a request inside either one plans; a request straddling both is refused,
  naming a tile of each and the offset (0.5 cell, as in N2).
- New M13, the reference node: the same tile's global index is the same
  whether planned alone, with neighbours, or under any `bounds`. The mosaic's
  `x_min` equals `X_ref + c0 · dx` bit for bit.
- New M14, coverage: a hole in the tile set inside the box is refused, naming
  the area. The same hole outside the needed region is NaN filler and is
  accepted.
- M5 becomes the cap: set from a patched physical-memory function, not a
  constant, and `load` is never called.
- New M15: `_adopt` is reachable only from `assemble`. The mosaic's array is
  read-only, and `np.shares_memory` holds between the canvas and the result
  (no copy), while the single-tile path still returns the loaded tile itself.
- T-real (the parked design's, `needs_codecs`): the committed benchmark tile
  cut into quadrants, 1-pixel overlap, meshes exactly like the whole tile.

**15b:**

- The domain in EPSG:4326 GeoJSON (no `crs` member), in EPSG:25832, and in
  WKT with `--domain-crs EPSG:3035`, each over a 25833 mosaic: the transformed
  vertices equal pyproj's own `always_xy` transform, and the mesh is built.
- **Axis order, able to fail:** a test transform without `always_xy` puts a
  vertex outside the DEM, and the extent check refuses it. Increment 11 §9's
  grep test finds exactly one `from_crs` site.
- `domain_crs` and `domain_transform` are recorded. A domain already in the
  DEM's CRS is not transformed, and the mesh is bit-identical to 16's.

**15c:**

- T-frame: for random lattices (spacing, offset, size), `(K · h, −R · h)`
  equals the planned `RasterMeta`'s node coordinates bit for bit, in every
  window.
- T-invariance: an analytic plane `z = a + b·lon + c·lat` on a geographic
  micro-mosaic meshes to two triangles with zero error, at any tolerance.
- T-gate: `to_core` of a geographic tile without a frame raises.
- T-anisotropy: a geographic extent reaching 20°S passes, and one reaching 30°S
  or 60°N is refused (the ratio is about 1/cos φ: 1.06, 1.15, 2.0).
- T-out: `--out-crs` equal to the DEM's CRS gives
  `max_reprojection_z_error_estimate` exactly 0; an LCC gives a positive value
  that grows with edge length.

**15d:**

- T-window: for every window, including windows on block edges and one-node
  windows, `decode_dem(s, window=w)` equals `decode_dem(s)` sliced to `w`,
  `meta` included. Tiled and stripped micro-TIFFs.
- The mosaic built with windows equals the one built from whole tiles.

### Pinned by the red suite (15a)

Names and behaviours the design left open, fixed by the red commit's suites
(`tests/python/test_mosaic.py`, `test_io_repository.py`, `test_io_read_meta.py`,
`test_dem_input.py`, `test_cli_mesh_mosaic.py`). All were run green against a
scratch implementation that was not committed; `test_mosaic.py` also killed
39 of 39 mutants of it.

**Two corrections to R4.** Both were found by the scratch implementation
failing M1's area-registered case.

- *Selection* (point 3) is by **index window**, not by the tile's node
  rectangle meeting the box. A tile is selected when its node index range
  meets the outward-snapped window, closed on both sides. The rectangle rule
  leaves nodes uncovered: area-registered neighbours at x 0..50 and 60..110
  with `bx_max = 55` snap to column 60, which only the east tile holds, and
  that tile's rectangle does not meet the box.
- *Coverage* (point 5) is **per node**: a node of the window that no selected
  tile covers is uncovered. The design's check was shapely's `covers` over
  the union of node rectangles, and it would refuse every pair of abutting
  area-registered tiles: their node rectangles are one spacing apart, so the
  union has a gap where no node lies.

**`tin_engine/mosaic.py`.**

- `ALIGN_TOLERANCE = 1e-6`; `MosaicError(ValueError)`;
  `physical_memory() -> int`, which is `SC_PHYS_PAGES * SC_PAGE_SIZE`, looked
  up at call time so a test can patch it.
- `Bounds(x_min=, y_min=, x_max=, y_max=)`: a non-finite, inverted or
  zero-extent box raises `ValueError`.
- `IndexWindow(row0, col0, rows, cols)`. `TilePlacement(name, meta, canvas,
  source)`: `canvas` is in canvas indices, `source` in the tile's own.
  Without bounds, `source` is the whole tile.
- `MosaicPlan(meta, reference, window, tiles)` is a frozen **Pydantic** model,
  because the tests reorder `tiles` with `model_copy`. `reference` is
  `(X_ref, Y_ref)`. `window` is the mosaic's first node as a global index,
  plus its shape. `tiles` is sorted by name. A tile's global origin is
  `window.row0 + canvas.row0 - source.row0`, and the same for columns.
- `Mosaic(tile, plan)`. `plan_mosaic(footprints, bounds=None, needed=None)`
  takes exactly those parameters. `needed` is a shapely geometry in the DEM's
  CRS, and it is closed: a node on its boundary is needed. Growing a domain by
  one cell is the caller's job (15b).
- **The cap** counts **the planned decoded dtype**: the itemsize of
  `np.result_type` over the *selected* footprints' `dtype`, and refuses when
  `rows * cols * itemsize > physical_memory() // 2`; the refusal names that
  dtype (`float64`) and `--bbox`. A float64 tile outside the window does not
  count. The cap is on the window, not on the lattice's union. (Amended after
  review, S2; the red suite had pinned 4 bytes a node, which under-counted a
  float64 mosaic by 2×.)
- **The mosaic's `RasterMeta`**: `nodata_source` comes from the first selected
  tile by name, and `vertical_unit_assumed` is true if any selected tile's is.
- **NoData against NoData:** two NaNs give NaN, and two sentinels give the
  sentinel (I5 needs both). NaN against the sentinel is NoData, and which of
  the two is not ruled, only that every order gives the same answer.
- **What messages name** (a case-insensitive substring match):
  - a lattice straddle: both tile names, `0.5 cell`, and `east-west` or
    `north-south`;
  - a CRS difference: both EPSG codes;
  - a spacing difference: `spacing` and the odd value;
  - a registration difference: `registration`;
  - a NoData difference: `nodata` and both values, with `None` for an absent
    sentinel;
  - an overlap disagreement: both names — the two tiles whose *values*
    disagree, not the first tile that merely covers the node (B2) — the count
    as a whole number, and the largest difference;
  - uncovered nodes: the node bounding box of the uncovered nodes, with each
    coordinate written out, not in scientific notation;
  - a changed tile: `changed since it was listed`;
  - the cap: `--bbox`, and the planned dtype.

**Amended after review (tests after green, S1-S3, B1, B2).** Pinned by the
test amendment that follows green `ff7cc8d`:

- **Q5 as read after review (B1),** in `TestB1LatticeByCoverage`
  (`test_mosaic.py`), the real-extract cases in `test_dem_input.py`'s
  `TestRealDtm10`, one CLI case in `test_cli_mesh_mosaic.py`, and the
  acceptance box against Ola's archive in `test_io_repository.py`'s
  `TestB1RealArchive` (headers only; skipped where the archive is absent). When
  one lattice's tiles are selected, nothing changes. When several are, **the
  nodes a lattice must cover** are its nodes in the request's box, snapped
  outward, clamped to the bounding box of every selected tile (any lattice) —
  *not* to that lattice's own union, or a box running past it into the other
  lattice's tiles would be silently cut short. Without a box, that is the
  bounding box of every tile. So a bare `--dem DIR` over two lattices is still
  refused as mixed-lattice (`test_no_bounds_selects_both_lattices_and_is_refused`
  keeps its outcome), as are all the other two-lattice refusals in M2 and C4.
  The chosen lattice's plan equals the plan of its tiles alone. The count is of
  a lattice's tiles in the repository, not of those the box selects; "its first
  tile" is its tile whose name sorts first.
- **The design's acceptance box selects nine tiles, not four.** With the
  51-node overlaps, the index-window selection (the first correction above)
  also selects the main-lattice neighbours 7807_1, 7808_2, 7809_3, 7809_4 and
  7909_3, whose overlap strips reach into the 2 × 2 block. The archive test
  asserts the four are in the plan, 7807_2 is not, and the window is
  10051 × 10051; the Acceptance's "`dem_tiles` listing four files" is not
  pinned. Whether to drop tiles that only duplicate covered nodes is open.
- **The coverage refusal's memory (S1).** Refusing a mostly uncovered window
  with no `needed` peaks, by tracemalloc, below the window's float32 canvas
  (index and coordinate arrays of the uncovered nodes cost 8× it). The
  `needed` path is not pinned.
- **The edge snap (S3)** has its test: at spacing 0.1, the box edges 0.3, 0.7
  and 0.9 add no node line.

**`io/models.py`.** `DemTile._adopt(meta, array)` is a classmethod. It raises
`ValueError` on a shape, dtype, ndim or non-C-contiguous mismatch, sets the
**passed** buffer read-only, and keeps it without a copy. Only `io/models.py`
and `mosaic.py` contain the string `_adopt`.

**`io/repository.py` and `io/geotiff.py`.**

- `read_meta(source, *, nodata=None)`.
- `TileFootprint(name=, meta=, dtype=)`, where `name` is the file's name
  (`path.name`) and `dtype` the numpy dtype the tile decodes to (a `np.dtype`,
  default float32), filled by the repository from the header through
  `geotiff.PROMOTION` (S2).
- `TiffDemRepository(paths, *, nodata=None)` and
  `TiffDemRepository.from_directory(directory, *, nodata=None)`. Construction
  reads no file.
- Refusals:
  - two paths with one file name: a `ValueError` naming the name;
  - an empty directory: a `ValueError` naming the directory;
  - a file that fails `read_meta`: a `GeoTiffError` naming the file;
  - `load` of an unknown name: `KeyError`.
- Every `open` in `repository.py` uses mode `"rb"`.
- No other `io/` module calls a file opener. The guard is an AST scan whose
  own scanner is tested on planted source.
- The docstring of `io/__init__.py` names `repository.py` and no longer says
  "Nothing here opens a file."

**`tin_engine/dem_input.py`.** `DemRequest(sources=, bounds=, nodata=)` is
frozen. A directory mixed with files, or two directories, raises a
`ValueError` whose message contains "director". No sources also raises.
`DemInput` has `tile`, `plan` and `label`. `label` is the directory's name, or
the stem of the first file as given.

**`cli.py`.**

- `--bbox` needs a short `metavar` (the scratch used `BOX`). Typer's default
  `<float float float float>` widens the help's type column until
  `--no-constraint-feet` is cut off at 80 columns, which fails the existing,
  unedited `test_cli_constraint_feet.py::TestTheFlag::test_the_help_names_the_flag`.
- `--bbox` without `--dem`, or with an invalid box, is a usage error naming
  `--bbox`. `--bbox` on a single file meshes a window of that file.
- **`dem_tiles`** holds the selected tiles' names, sorted and `; `-joined. It
  is recorded whenever `--dem` is a directory or several files, even when only
  one tile is selected, because the label is then the directory's and the file
  used must be on record. It is never recorded for a single file.
- The `mosaic of N tiles, R x C nodes; ` prefix and the stderr line appear only
  for N ≥ 2.

**Fixtures.** `tests/fixtures/dtm10/`, cut by `tests/fixtures/dtm10/extract.py`
from Ola's archive (`DTM10_UTM33_20220924`), one release, © Kartverket, CC BY
4.0 (checked by the main session, 2026-09-27, against Geonorge's metadata API
for dataset `dddbb667-1303-4ac5-8640-7ec04c0e3918`: "Åpne data", CC BY 4.0):

- `seam/`: 6400_4 | 6400_1, rows 3072-3327, 307 columns each, with a 51-column
  overlap. It agrees bit for bit, and a one-column shift disagrees; the script
  asserts both.
- `lattices/`: 7707_1 over 7707_2, 96 × 96 each, 5 m east-west apart.

Deflate, 0.6 MB together.

## Test data

**Norway (15a, 15b).**

- *Synthetic:* micro-TIFF mosaics, 2 × 2, overlap 0, 1 and 3, point- and
  area-registered; a two-lattice repository with one tile shifted half a cell.
- *Real, split:* the committed benchmark tile cut into quadrants (above).
- *Real, seam:* an extract of a real DTM10 edge overlap from Ola's archive,
  for example 6400_1 | 6400_4 (N3's first edge pair). Two windows of about
  256 × 307 nodes (256 plus the 51-node overlap), with their real
  georeferencing, Deflate-compressed: well under 1 MB. **Cut both from the
  archive, never one from the committed tile** (N5). Same attribution as the
  committed tile. The extraction script is test support and goes under
  `tests/`, as @tester's.
- *Real, two lattices:* an extract of one of N2's eight tiles and an aligned
  neighbour, for `refuses_mixed_lattice`. @tester picks the pair with
  `dtm10_probe.py headers`.
- *@perf, Norway:* the 2 × 2 block 7908_3, 7908_2, 7808_4, 7808_1 from the
  archive (all on the main lattice): 10 051² = 101 M nodes, 404 MB float32,
  4× the benchmark tile. Local only, not committed.

**The basin (15c, 15d; after Norway).**

- *Synthetic:* a geographic micro-TIFF with ANADEM's exact header numbers
  (EPSG:4326, area-registered, spacing 0.00026949458523585647, tie point
  (−48.00102804928458, −8.129843152810082), NoData −9999). `micro_tiff` needs
  a geographic-CRS parameter, which is test code.
- *Real, seam:* extracts of 23L | 24L (the 42°W seam, bit-equal overlap) and
  24L | 24K (the 87-row overlap with NoData on one side), cut by range reads
  as in `anadem_probe.py seams`. Four 512² blocks, about 2.4 MB compressed.
  Whether they may be committed is Q10.
- *@perf, a real piece of the basin:* **46–44°W × 17.26–15.26°S**, across
  the 23L/23K seam, over the São Francisco near Pirapora–Januária and the
  Paracatu. It is believed to lie inside the basin; @perf checks that against
  BHO. 7422 × 7422 = **55 M nodes, 210 MiB float32**, about 2.2× the
  benchmark. About 240 blocks by range reads, roughly 150 MB of download
  (from the measured 0.5–0.85 MB per block). The Espinhaço window of B6 (2.4 M
  nodes, 45 km, steep) is a quick steep-terrain case.
- *@perf, the whole basin:* the four main tiles plus slivers, about 8 GB of
  download; the box is 8.6 GiB at float32 (B5).

## Acceptance

- **15a:** `rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20220924/ --bbox
  … --tolerance 1` on the 2 × 2 block writes a `.vtk` that ParaView opens,
  with `dem_tiles` listing the four block tiles among the main-lattice tiles
  the index window selects (nine: with 51-node overlaps it also takes 7807_1,
  7808_2, 7809_3, 7809_4 and 7909_3; corrected at the 15a test amendment), and
  not 7807_2, whose lattice does not cover the box (Ola's Q5 reading). A box
  that straddles one of N2's odd tiles and a neighbour so that neither lattice
  covers it is refused, naming both. The quadrant split of the benchmark
  tile meshes identically to the tile. `--dem file.tif` output unchanged.
- **15b:** a catchment polygon in EPSG:4326 over the archive meshes, with
  `domain_crs` and `domain_transform` recorded.
- **15c:** the basin piece above meshes with `--out-crs` set to the basin
  LCC, and records `computation_frame` and the reprojection estimate.
- **15d:** @perf meshes the basin piece and then the whole basin on Ola's Mac,
  recording peak memory, decode time and refine time, under
  `docs/benchmarks/<date>/`, with the power state.
- For all: every gate in `CLAUDE.md` §4 green, and CI green.

## Questions for Ola (Q1-Q5 ruled 2026-09-27, see "Ruled by Ola"; Q6-Q10 open)

Q1-Q5 are about Norway and are needed before 15a starts. Q6-Q10 are about
the basin and can wait until after Norway.

**Q1 (the parked U1). Two tiles give different valid values at the same node.**
Now measured: DTM10 overlaps agree bit for bit on 4 of 869 pairs (N3), and
ANADEM's on 3 seams (B4). But the committed benchmark tile and the archive's
tile of the same name differ at 910 706 nodes (N5), so two releases in one
directory is a real case.
- **(a) Refuse, naming both tiles, the count and the largest difference.
  Recommended.** Valid still beats NoData, and the result does not depend on
  tile order.
- (b) First tile by sorted name wins, silently. That is what the legacy did in
  effect.
- (c) Refuse by default, and add a `--overlap first` flag. About 15 lines more.

**Q2 (the parked U2). `--bbox`.**
- **(a) The union by default, plus an optional `--bbox` in the DEM's CRS.
  Recommended.** Without it, `--dem DIR` on your archive is refused by the cap
  (72 GiB box), and the only way round would be copying files into a smaller
  directory. It is about 20 lines.
- (b) The union only.

**Q3 (the parked U3). Where the repository and its paths live.** `io/__init__.py`
says nothing in `io/` opens a file.
- **(a) `io/repository.py`, with `io/`'s rule narrowed to "decoders and
  encoders take streams; `io/repository.py` is the one module that opens
  files, read-only". Recommended.** Storage access sits next to the formats it
  reads.
- (b) A top-level `tin_engine/dem_repository.py`, leaving `io/`'s rule as it is.
- (c) Paths in `cli.py`. Not recommended: `cli.py` is 1101 lines, and a GUI
  or API worker could not reuse it.

**Q4 (the parked U4). Tiles with different NoData sentinels.** Your ruling
C2 (b) in increment 18 already assumed (a): "differing sentinels are refused".
Both real archives have one sentinel each (N1, B2).
- **(a) Refuse. Recommended; please confirm.**
- (b) Rewrite each tile's sentinel to NaN and give the mosaic no sentinel.
  About 10 lines.

**Q5 (new). The eight half-cell tiles in your archive (N2).** 7304_1,
7507_4, 7606_2, 7707_1, 7707_3, 7807_2, 7807_3 and 7808_3 sit 5 m east-west
off the other 246, and are one or two nodes larger.
- **(a) Refuse a request that mixes them with the main lattice; a request
  inside their own lattice meshes. Recommended.** Nothing is resampled. You
  may know where they came from: a re-delivery of those sheets on a shifted
  grid would explain it, and a fresh download of them might be on the main
  lattice.
- (b) Resample them onto the main lattice (bilinear, a half-cell shift, so
  every value becomes the average of two). That is resampling (R8), and its
  own increment.

**Q6 (the basin). The computation frame.**
- **(a) The lattice frame: mesh in an affine image of the DEM's own grid,
  transform the domain in and the mesh out. Recommended.** The tolerance holds
  exactly at the DEM's own nodes. Triangle shapes are judged in a frame up to
  6.5 % out of square on the ground at the basin's south edge.
- (b) Resample onto a square grid in a projected CRS first. Shapes are judged
  on the ground, but the tolerance then refers to the resampled grid: 13 m off
  the DEM at worst on real Espinhaço terrain (B6). It also costs a second
  array the size of the basin's box.

**Q7 (the basin). The output CRS for a geographic DEM.**
- **(a) Required as `--out-crs`, with a suggested basin-fitted LCC in the
  refusal message. Recommended.** The output CRS is a contract with whatever
  reads the mesh, and a CRS fitted to each domain would put two catchments, or
  two subdomains, in different CRSs.
- (b) Fit an LCC to the domain automatically, and record it as WKT2.
- (c) Write longitude and latitude. Exact, with no bending, but in degrees.

**Q8 (the basin). The anisotropy limit and the square frame.**
- **(a) A square frame (21b's integer incircle applies) and a limit of 10 %:
  the basin passes (6.5 %); geographic DEMs above about 25° latitude are
  refused. Recommended.**
- (b) A frame fitted to the extent's mean latitude: shapes truer (±3 % over
  the basin), cells not square, so the filtered incircle is used and 21b's
  gain (−15 % refine at 8 threads) is lost for geographic DEMs.

**Q9 (the basin). Memory for the whole basin.** Its box at 30 m is 8.6 GiB,
and the design caps the canvas at half of physical memory.
- **(a) A dense canvas, as designed: the basin fits on your 32 GiB Mac.
  Recommended for now.** Is 32 GiB the machine the basin must run on, or must
  it also run on a smaller one?
- (b) A block-sparse raster in C++ holding only the 512² blocks that meet the
  basin (about 2.7 GiB). Its own increment: a new `RasterSource`, the concept
  change 18 R6 describes, and a binding.
- (c) Go to domain decomposition, which also solves memory.

**Q10 (the basin). May an ANADEM extract be committed as a test fixture?**
ANADEM's repository is MIT-licensed, but that licence names "the Software",
and the data derives from Copernicus GLO-30, which carries its own attribution
terms.
- **(a) Commit about 2.4 MB of extracts with both notices, after you or the
  authors confirm the data licence. Recommended.**
- (b) Download on demand in a test marked `network`, skipped in CI.

Still open from the research note, for @perf's basin run: the tolerance
intended for the basin (it sets the triangle count more than anything else),
and whether commercial use matters.
