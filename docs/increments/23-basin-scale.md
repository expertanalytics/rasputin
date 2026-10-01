# Increment 23: basin scale — pieces cut on constraints, windowed DEM reads, a tile cache

Status: **designed by `@architect`, 2026-10-01; not implemented, nothing
ruled yet.** Design only, written before `@tester` per
`docs/increments/README.md` step 1. The questions for Ola are B1-B12 at the
end; each has a recommendation and its cost.

## Why this record, and why its name

Ola asked on 2026-10-01 for a basin-scale design in an increment file of its
own (`15c-geographic-dem.md`, Q15). It gets the next free number, 23 (22 is
auto-catchment). "Basin scale" because what it delivers is the ability to
mesh the São Francisco basin (635,194.5 km², 734 M ANADEM nodes) at any
tolerance on one machine, with memory set by the size of a piece, not of the
basin. It is not a sub-increment of 15: 15 is about DEM input, while this
changes the refinement core (frozen edges, a seam pass), the run (many
pieces, in parallel) and the output (a piece per file and an index). It
absorbs 15d's window decoding and supersedes 15's Q9 ("32 GB is fine", a
dense canvas) at basin scale.

Sub-increments, each its own PR: **23a-1** (windowed decoding and the tile
cache, reading), **23a-2** (the fetch step), **23b** (frozen edges and the
seam pass, in C++), **23c** (the decomposition, run piece by piece), **23d**
(pieces in parallel, resumable runs, `rasputin stitch`). Order and sizes are
under "Order of work".

## Ola's direction, 2026-10-01

Quoted from `15c-geographic-dem.md` (Q15) and the brief for this record:

- **Locality, not memory, decides.** "So, one of our design criterias is
  locality. When things are tiled on file, and we refine, why are we talking
  about memory limitation? There must be easy fixes to this."
- **Cut along constraint edges.** "We could easily cut on constraint-edges,
  right? Only two intersecting subdomains sharing the constraint is needed to
  resolve potential corrections on the boundary, as these edges are never
  flipped."
- **Windowed DEM reads, no dense canvas.** This revisits Q9.
- **A fetch step that caches tiles.** "I think we should cache the tiles we
  download, and probably do this as a preprocessing step, or our performance
  will drop while we wait for download." Meshing never touches the network; a
  missing tile is a refusal that names the fetch command; the cache has a
  manifest, and fetching is resumable and incremental.

The main session's reasoning, which this record was asked to check, is
checked under "The subdomain model" (the union argument) and "The seam
protocol" (where the split exchange goes).

## What was measured for this design

By `@architect`, 2026-10-01, with the two scripts in `docs/increments/23-probes/`
(run as `python docs/increments/23-probes/<script> ../rasputin_data/sao_francisco_piece`;
both were run for this record and print the figures below).

**BHO is an exact coverage** (`bho_coverage.py`). The 1,163 BHO 2017 5k
elementary catchments of ottobasin 76949 (the Velhas piece of
`docs/benchmarks/2026-10-01/basin-piece/README.md`):

- every polygon edge occurs exactly twice (76,459 edges) or once (7,309, the
  outline), bit for bit: neighbouring catchments share their boundary
  vertices as identical doubles, as CORINE's do (16b M2);
- the sum of the areas equals the union's area to 4.4e-16, relative;
- the ottocode prefixes give the Pfafstetter levels inside the piece:

| level | units | km², min / median / max | largest bbox at 30 m | shared boundary |
|---:|---:|---|---:|---:|
| 6 | 9 | 56 / 792 / 5,201 | 11.4 M nodes | 421 km |
| 7 | 73 | 1 / 66 / 1,014 | 2.4 M nodes | 1,768 km |
| 8 | 431 | 0 / 13 / 218 | 0.6 M nodes | 4,577 km |

  For comparison, a grid of 2048-node blocks at 30 m cuts the piece along
  396 km. The outline's 7,309 edges over 785.9 km are 107.5 m apart on
  average, so a BHO boundary used as a cut brings a vertex every ~107 m.
- Pfafstetter units are unbalanced: at level 6 the largest is 93 times the
  smallest and holds 45 % of the piece.

**ANADEM's COG serves as a block cache source** (`anadem_cog.py`; the header
facts in `15c-geographic-dem.md` "Measured by @architect" are not repeated):

- a `HEAD` gives `Content-Length: 70928454424`, `Accept-Ranges: bytes`,
  `Last-Modified: Sun, 07 Jun 2026 02:04:54 GMT`, and **`ETag:
  "00000000000000000000000000000000-1"`, a placeholder**: the ETag cannot
  detect a changed file;
- tifffile parses every one of the 11 pages from the first 8 MiB alone (a
  reader that refuses any read past the prefix did not fire). The
  full-resolution page has 188,638 blocks of 512²; image data starts at byte
  23,716,931,662 (the overviews come first); blocks are in row-major order
  with small gaps between them;
- one block over the Velhas piece (index 90561, 671,506 bytes) was fetched by
  its own range request and decoded alone by `TiffPage.decode`: 512² float32,
  heights 713.4-1,363.7 m. So a block's bytes plus the parsed header are
  enough to decode it; no file is needed;
- **the whole basin is 3,061 blocks, 1.71 GiB compressed** (blocks meeting
  the BHO level-2 outline grown by three source cells; 2.99 GiB as float32);
  the Velhas piece is 82 blocks, 0.05 GiB. The block median is 1,085 bytes
  over the whole COG (ocean and NoData compress to almost nothing) and
  611,639 bytes over the basin.

**GLO-30's tile list exists**: `https://copernicus-dem-30m.s3.amazonaws.com/tileList.txt`
answered 200 and lists names like `Copernicus_DSM_COG_10_N00_00_E006_00_DEM`,
so which 1° tiles exist (land) is known without probing for 404s.

## Prior art: legacy and literature

### Literature

Bibliographic details marked **(verified)** were checked against Crossref on
2026-10-01 (title, authors, venue, volume, pages, DOI); the description of
the method is recalled unless it says otherwise. No full text was reread.

- **Streaming Delaunay.** Isenburg, Liu, Shewchuk and Snoeyink, "Streaming
  computation of Delaunay triangulations", *ACM Trans. Graphics* 25(3):
  1049-1056, 2006, doi:10.1145/1141911.1141992 (verified). Points arrive in
  a stream with *finalization* tags: once a region of space is known to get
  no more points, the triangles whose circumcircles lie in finalized space
  are written out and freed, so memory follows the front, not the data. The
  companion "Generating raster DEM from mass points via TIN streaming",
  GIScience 2006, LNCS 4197:186-198, doi:10.1007/11863939_13 (verified),
  goes the other way (points to raster). **What we take:** locality as the
  memory principle, which is Ola's, and "finish a region, write it out, free
  it": a piece's file is written when the piece is done. **What differs:**
  they triangulate every point, with no constraints and no refinement; we
  choose points by error and must honour constraints. Their memory bound
  comes from spatial finalization of a point stream; ours from pieces bounded
  by constraints, which are independent by construction, so no front exists.
- **I/O-efficient CDT.** Agarwal, Arge and Yi, "I/O-efficient construction
  of constrained Delaunay triangulations", ESA 2005, LNCS (volume recalled
  as 3669) 355-366, doi:10.1007/11561071_33 (verified). An external-memory
  CDT. **Not used:** the vector input (outline, cuts, features) fits in
  memory at basin scale; only the DEM is large, and it is partitioned.
- **Decoupled refinement.** Linardakis and Chrisochoides, "Delaunay
  decoupling method for parallel guaranteed quality planar mesh refinement",
  *SIAM J. Sci. Comput.* 27(4):1394-1423, 2006, doi:10.1137/030602812
  (verified); and the graded version, *SIAM J. Sci. Comput.* 30(4):1875-1891,
  2008, doi:10.1137/060677276 (verified). Separators are refined *before* the
  subdomains, so that each subdomain then refines with no communication and
  the union keeps the Delaunay and quality guarantees. **This is the method
  we build on.** **What differs, and why:** their separators are refined from
  the quality criterion and a sizing function, because what a subdomain will
  need on its boundary is only known through quality bounds. Ours are refined
  by their own one-dimensional tolerance, which is fully known in advance:
  the TIN restricted to a constraint edge is the linear interpolant of the
  edge's two end heights, whatever the triangles on either side are, so the
  error at any point on the edge is a function of the edge's vertices alone
  (the seam pass). **The guarantee we drop:** in their method the separator
  edges end up Delaunay, so the union is Delaunay. Ours stay constraints, so
  the union is constrained Delaunay *with respect to the seams too* (B4).
- **Interfaces first.** Galtier and George, "Prepartitioning as a way to mesh
  subdomains in parallel", 5th International Meshing Roundtable, 1996
  (**unverified**: not in Crossref; the IMR 1996 proceedings carry no DOI).
  Galtier, "Load balancing issues in the prepartitioning method", LNCS,
  1997, pp. 922-936, doi:10.1007/bfb0002835 (verified) is the same line.
  Partition, mesh the interfaces, then mesh the subdomains independently.
  Structurally this is our order (seam pass, then pieces).
- **Exchange on shared segments.** Chernikov and Chrisochoides, "Algorithm
  872: parallel 2D constrained Delaunay mesh generation", *ACM Trans. Math.
  Softw.* 34(1):1-20, 2008, doi:10.1145/1322436.1322442 (verified); and Kot,
  Chernikov and Chrisochoides, "Parallel out-of-core constrained Delaunay
  mesh generation", IDAACS 2005, pp. 183-190, doi:10.1109/idaacs.2005.282967
  (verified). Recalled: when a subdomain splits a segment on its boundary, it
  sends the split to the neighbour sharing that segment, so the two agree;
  the out-of-core version swaps subdomains to disk. **This is the protocol
  Ola describes** ("only two intersecting subdomains sharing the constraint
  is needed"). It is the alternative of B5 (b). Its termination rests on
  Ruppert-type bounds; ours would rest on finite candidate sets (below).
- **Parallel Delaunay by domain decomposition**, for context: Said,
  Weatherill, Morgan and Verhoeven, "Distributed parallel Delaunay mesh
  generation", *Comput. Methods Appl. Mech. Engrg.* 177(1-2):109-125, 1999,
  doi:10.1016/s0045-7825(98)00374-0 (verified); Blelloch, Miller, Hardwick
  and Talmor, "Design and implementation of a practical parallel Delaunay
  algorithm", *Algorithmica* 24(3-4):243-269, 1999, doi:10.1007/pl00008262
  (verified); Wu, Guan and Gong, "ParaStream: a parallel streaming Delaunay
  triangulation algorithm for LiDAR points on multicore architectures",
  *Computers & Geosciences* 37(9):1355-1363, 2011,
  doi:10.1016/j.cageo.2011.01.008 (verified), the terrain case of streaming
  plus parallelism, again for all points. And Chrisochoides, "Parallel mesh
  generation", in *Numerical Solution of PDEs on Parallel Computers*, LNCSE
  51, pp. 237-264, 2006, doi:10.1007/3-540-31619-1_7 (verified). Increment
  21 cites this chapter as "A survey of parallel mesh generation methods";
  Crossref's title is "Parallel Mesh Generation".
- **The union is constrained Delaunay.** Chew, "Constrained Delaunay
  triangulations", *Algorithmica* 4:97-108, 1989, doi:10.1007/bf01553881
  (verified). A triangle is constrained Delaunay when no vertex *visible*
  from its interior lies inside its circumcircle, where constraints block
  visibility. The argument under "The subdomain model" rests on this
  definition.
- **Tiled terrain simplification.** Campos, Quintana, Garcia, Schmitt and
  Spoelstra, "3D simplification methods and large scale terrain tiling",
  *Remote Sensing* 12(3):437, 2020, doi:10.3390/rs12030437 (verified; first
  found by the main session for increment 21): greedy insertion among other
  methods, run tile by tile in parallel, keeping tile-border vertices shared
  between neighbours. Cignoni, Ganovelli, Gobbetti, Marton, Ponchio and
  Scopigno, "BDAM — Batched Dynamic Adaptive Meshes for high performance
  terrain visualization", *Computer Graphics Forum* 22(3), 2003,
  doi:10.1111/1467-8659.00698 (found by OpenAlex; authors and issue
  recalled): TIN patches in a hierarchy, simplified to an error with
  boundaries shared consistently between neighbours. **The closest terrain
  prior art.** **What differs:** both are for rendering, and neither is
  recalled as giving a sup-norm tolerance *on the shared border itself*
  against the DEM between nodes; ours does, by the seam pass. Not reread, so
  this difference is to be confirmed before it is written anywhere public.
- **The seam pass's method.** Douglas and Peucker, "Algorithms for the
  reduction of the number of points required to represent a digitized line
  or its caricature", *Cartographica* 10(2):112-122, 1973,
  doi:10.3138/fm57-6770-u75u-7727 (verified): insert the worst point, split,
  repeat. The seam pass is that in the vertical (height error at check
  points, not lateral distance), which is greedy insertion (Garland and
  Heckbert, CMU-CS-95-181, 1995, as in increment 14) in one dimension.
- **Pfafstetter coding.** Verdin and Verdin, "A topological system for
  delineation and codification of the Earth's river basins", *J. Hydrology*
  218(1-2):1-12, 1999, doi:10.1016/s0022-1694(99)00011-6 (verified). Each
  level divides a basin into nine units (four tributary basins, five
  interbasins), coded by one more digit; BHO's ottocodes are this.
- **Cloud-optimised GeoTIFF.** OGC 21-026, "OGC Cloud Optimized GeoTIFF
  Standard", doi:10.62973/21-026 (found in Crossref, no year given; recalled
  as 2023). Headers and offset arrays at the start, tiles addressable by
  range request. The probe above confirms ANADEM's copy behaves so.

**Novelty: none claimed.** One combination could look new: *a greedy terrain
refinement to an exact sup-norm tolerance, decomposed into pieces along
constraints, where the pieces need no communication because each seam is
refined beforehand by its own one-dimensional error.* Searched on
2026-10-01: Crossref bibliographic queries for every citation above, plus
"parallel construction of triangulated irregular network from DEM",
"tile-based TIN generation large DEM seams" and "out-of-core terrain
simplification TIN large DEM"; OpenAlex full-text search for "parallel
terrain simplification domain decomposition error bound" (4,791 results,
nothing on the subject in the top hits), "TIN generation DEM parallel
partition seam" (51), "greedy insertion terrain triangulation parallel tiles"
(17), "decoupled constrained Delaunay terrain subdomains separator
refinement" (0), "watershed partition parallel mesh generation hydrological"
(854, hydrological models), "out-of-core triangulated irregular network
construction large DEM" (618), "terrain approximation error guarantee tiles
boundary consistency" (402; BDAM and ROAM, rendering) and "Pfafstetter
subbasin parallel mesh" (1, unrelated). Nothing found does this, but Google
Scholar, IEEE Xplore and Scopus were not searched (no access from this
round), and Campos et al. and BDAM were not reread. So no claim is made;
before one is, read those two and run increment 21's search list
(`21-parallel-refine.md`, "Searches to run before any novelty claim").

### Legacy

```sh
$ git grep -l -iE 'subdomain|decompos|partition|seam|stitch' legacy-archive -- legacy
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/triangulate_dem.h
legacy-archive:legacy/rasputin/web_visualize.py
$ git grep -l -iE 'download|urlopen|requests\.|http|cache|fetch' legacy-archive -- legacy
legacy-archive:legacy/rasputin/avalanche.py
legacy-archive:legacy/rasputin/web_visualize.py
legacy-archive:legacy/rasputin/wfs_repository.py
```

Read: the first set is a docstring wish ("partition the whole area into
disjoint topologies", `application.py:27`, repeated in `web_visualize.py`)
and `partition()` in `triangulate_dem.h:833`, which splits a finished face
list in two by a predicate (lakes by slope). Neither decomposes a meshing
problem. The second set fetches over the network *at run time*, with no
cache: `avalanche.py:30` calls `requests.get` on NVE's API, and
`wfs_repository.py` reads a WFS service; `web_visualize.py`'s hits are
`localhost` URLs in help text. That is the pattern Ola's direction rules out
(meshing must never wait on the network), so **nothing is carried**. No
legacy constant is re-derived, so no `@migration-expert` pass is needed.

## The subdomain model

(pending)

## Choosing the cuts

(pending)

## The seam protocol

(pending)

## Tolerance, the final check and the constraint check points per piece

(pending)

## Windowed source reads

(pending)

## The fetch step and the tile cache

(pending)

## Memory and parallelism

(pending)

## Output

(pending)

## What survives of 15c, the edge strip and 15d

(pending)

## Order of work and PR split

(pending)

## Tests @tester can write red

(pending)

## @perf acceptance

(pending)

## Questions for Ola

(pending)
