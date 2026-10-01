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
- tifffile parses every one of the 11 pages from the first 8 MiB alone, with
  every page's offsets numbering its blocks and no read past the prefix
  (tifffile checks each offset against the stream's size before reading and
  logs a refusal rather than raising, so the strict wrapper never sees an
  out-of-range read; the offsets count is the check). Cut to 64 KiB or
  1 MiB, tifffile still returned 11 pages, logged missing tags and gave pages with no offsets (1 MiB: 0 offsets for 188,638
  blocks on the full page), and the wrapper recorded nothing: the offsets
  check is the one that fires, and the probe failed on it. The
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

How each entry was checked. The first pass, on 2026-10-01, had only Crossref
and OpenAlex; the same day the section was redone with web search and full
texts where they could be fetched. The marks:
**Crossref**: the DOI's record matches authors, title, venue, volume and
pages as given. **read**: the full text was read at the URL named, and what is
said about the method comes from it. **code**: the source was read at the
commit named. **secondary**: the work itself was not reached; what is said
comes from the source named. **recalled**: from memory, unchecked.

- **Streaming Delaunay.** Isenburg, Liu, Shewchuk and Snoeyink, "Streaming
  computation of Delaunay triangulations", *ACM Trans. Graphics* 25(3):
  1049-1056, 2006, doi:10.1145/1141911.1141992 (Crossref; read,
  `https://people.eecs.berkeley.edu/~jrs/papers/dtstream.pdf`). Points arrive
  with *finalization* tags; a triangle is written out and freed once its
  circumcircle meets no unfinalized cell, so memory follows the front, not
  the data. The paper says "we have not yet implemented support for
  constrained Delaunay triangulations". The companion "Generating raster DEM
  from mass points via TIN streaming", GIScience 2006, LNCS 4197:186-198,
  doi:10.1007/11863939_13 (Crossref), goes from points to a raster. **What we
  take:** locality as the memory principle, which is Ola's, and "finish a
  region, write it out, free it": a piece's file is written when the piece is
  done. **What differs:** they triangulate every point, with no constraints
  and no refinement; we choose points by error and must honour constraints.
  Their memory bound comes from spatial finalization of a point stream; ours
  from pieces bounded by constraints, which are independent by construction,
  so no front exists.
- **I/O-efficient CDT.** Agarwal, Arge and Yi, "I/O-efficient construction
  of constrained Delaunay triangulations", ESA 2005, LNCS 3669:355-366,
  doi:10.1007/11561071_33 (Crossref; the volume number, recalled in the first
  pass, is from web search results for the book's Springer and DBLP pages). Secondary (Isenburg et al.
  2006, §2): a small random sample splits the points into subproblems that
  Triangle solves in core, over many read passes. **Not used:** the vector
  input (outline, cuts, features) fits in memory at basin scale; only the DEM
  is large, and it is partitioned.
- **Decoupled refinement.** Linardakis and Chrisochoides, "Delaunay
  decoupling method for parallel guaranteed quality planar mesh refinement",
  *SIAM J. Sci. Comput.* 27(4):1394-1423, 2006, doi:10.1137/030602812
  (Crossref); and "Graded Delaunay decoupling method for parallel guaranteed
  quality planar mesh generation", *SIAM J. Sci. Comput.* 30(4):1875-1891,
  2008, doi:10.1137/060677276 (Crossref; read,
  `https://class-pages.pages.cs.odu.edu/crtcpub/assets/publications/pdf/graded-delaunay-decoupling-method-for-parallel-guaranteed-qu-2008.pdf`).
  Read in the 2008 paper: the separators "will be a permanent part of the
  geometry"; a preprocessing step refines them to edge lengths between
  (√3/2)k and 2k, with k from the triangle-area bound and the local feature
  size (Theorem 5.3, for Ruppert's algorithm), so that the unmodified
  sequential mesher (Triangle) run on each subdomain never splits a separator
  edge; the union is then a *conforming Delaunay* triangulation, because the
  open diametral circles of the separator edges stay empty (Proposition 5.2).
  **This is the method we build on.** **What differs, and why:** their
  separator lengths come from the quality bound and a sizing function,
  because what a subdomain needs on its boundary is known only through those
  bounds. Ours come from the seam's own one-dimensional tolerance, which is
  known in advance: the TIN restricted to a constraint edge is the linear
  interpolant of the edge's two end heights, whatever the triangles on either
  side are, so the error at any point on the edge is a function of the
  edge's vertices alone (the seam pass). And they keep the mesher off the
  separators by giving it nothing to split there; we forbid the split
  outright (frozen edges), because a tolerance-driven refine has no size
  bound that would keep it away. **The guarantee we drop:** their union is
  Delaunay; ours is constrained Delaunay with the seams as constraints, which
  is PCDM's guarantee (below) (B4).
- **Interfaces first.** Galtier and George, "Prepartitioning as a way to mesh
  subdomains in parallel", Proc. 5th International Meshing Roundtable, 1996,
  pp. 107-121 or 107-122 (citing works differ). **Not reached:**
  `imr.sandia.gov` answered 403 to every fetch, no other copy was found, and
  the paper has no DOI. Secondary: Linardakis and Chrisochoides 2008 (§2, read
  above) describe it as "a parallel projective Delaunay meshing method which
  guarantees the quality of the elements and eliminates communication, but
  may suffer setbacks in the form of regenerating part of the mesh"; search
  summaries of citing works say the interfaces between subdomains are meshed
  before the subdomains. Galtier, "Load balancing issues in the
  prepartitioning method", LNCS, 1997, pp. 922-936, doi:10.1007/bfb0002835
  (Crossref) is the same line. Structurally this is our order (seam pass,
  then pieces).
- **Exchange on shared segments.** Chernikov and Chrisochoides, "Algorithm
  872: parallel 2D constrained Delaunay mesh generation", *ACM Trans. Math.
  Softw.* 34(1):1-20, 2008, doi:10.1145/1322436.1322442 (Crossref; read,
  `https://class-pages.pages.cs.odu.edu/crtcpub/assets/publications/pdf/algorithm-872-parallel-2d-constrained-delaunay-mesh-generati-2008.pdf`);
  and Kot, Chernikov and Chrisochoides, "Parallel out-of-core constrained
  Delaunay mesh generation", IDAACS 2005, pp. 183-190,
  doi:10.1109/idaacs.2005.282967 (Crossref). Read: subdomains are separated
  by constrained segments, and "if the mesh inside each subdomain is
  Delaunay, then the global mesh is constrained Delaunay"; an encroached
  boundary edge is split at its midpoint and the neighbour is sent a split
  message, `split(p0, p1, p2)`, with care for messages that arrive out of
  order; termination is detected with Dijkstra's algorithm. **This is the
  protocol Ola describes** ("only two intersecting subdomains sharing the
  constraint is needed"). It is the alternative of B5 (b), and its union
  guarantee (constrained Delaunay with respect to the separators) is the one
  B4 (a) gives.
- **Parallel Delaunay and terrain triangulation**, for context: Said,
  Weatherill, Morgan and Verhoeven, "Distributed parallel Delaunay mesh
  generation", *Comput. Methods Appl. Mech. Engrg.* 177(1-2):109-125, 1999,
  doi:10.1016/s0045-7825(98)00374-0 (Crossref); Blelloch, Miller, Hardwick
  and Talmor, "Design and implementation of a practical parallel Delaunay
  algorithm", *Algorithmica* 24(3-4):243-269, 1999, doi:10.1007/pl00008262
  (Crossref); Wu, Guan and Gong, "ParaStream: a parallel streaming Delaunay
  triangulation algorithm for LiDAR points on multicore architectures",
  *Computers & Geosciences* 37(9):1355-1363, 2011,
  doi:10.1016/j.cageo.2011.01.008 (Crossref), streaming plus parallelism for
  all points; Chrisochoides, "Parallel mesh generation", in *Numerical
  Solution of PDEs on Parallel Computers*, LNCSE 51, pp. 237-264, 2006,
  doi:10.1007/3-540-31619-1_7 (Crossref; increment 21 had cited it as "A
  survey of parallel mesh generation methods", corrected on this branch). And
  Puppo, Davis, DeMenthon and Teng, "Parallel terrain triangulation", *IJGIS*
  8(2):105-128, 1994, doi:10.1080/02693799408901989 (Crossref; abstract via
  search summary): greedy selection of grid points into a Delaunay TIN on a
  CM-2, parallel by inserting many points per pass, not by decomposing the
  domain. Kang, Lee, Yang and Park, "A fast digital terrain simplification
  algorithm with a partitioning method", HPC Asia 2000, vol. 2, pp. 613-618,
  doi:10.1109/hpc.2000.843506 (Crossref; abstract, OpenAlex; paywalled, not
  read): greedy insertion run block by block over equal rectangular blocks,
  each insertion looking only at the current block's points and triangles,
  with the block corners inserted first; a serial speed-up (4 to 20 times),
  and the abstract says nothing of how block borders are refined.
- **The union is constrained Delaunay.** Chew, "Constrained Delaunay
  triangulations", *Algorithmica* 4:97-108, 1989, doi:10.1007/bf01553881
  (Crossref). A triangle is constrained Delaunay when no vertex *visible*
  from its interior lies inside its circumcircle, where constraints block
  visibility. The argument under "The subdomain model" rests on this
  definition.
- **Tiled terrain simplification: the closest terrain prior art.** Four
  works, each reread for this pass.
  - Campos, Quintana, Garcia, Schmitt, Spoelstra and Schaap, "3D
    simplification methods and large scale terrain tiling", *Remote Sensing*
    12(3):437, 2020, doi:10.3390/rs12030437 (Crossref; **corrected**: six
    authors, the first pass dropped Schaap). The paper itself refused every
    fetch (MDPI and the Girona repository answered 403 or a bot check), so
    the method was read in the authors' code and its documentation: code,
    `coronis-computing/emodnet_qmgc` at `036aa9c`
    (`src/tin_creation/tin_creation_greedy_insertion_strategy.cpp`,
    `src/tin_creation/tin_creation_simplification_point_set.cpp`,
    `src/base/zoom_tiles_scheduler.h`), and its wiki at `438824c`
    ("General Parameters", "Point Set Simplification Parameters"). Tiles are
    meshed one by one with a choice of methods (greedy insertion after
    Garland and Heckbert with a vertical error, quadric edge collapse,
    point-set simplification). **The first tile built decides a shared
    border; neighbours built later take its border vertices as fixed** ("we
    keep track of tiles that are already built, and we maintain and transfer
    the border vertices to the next tiles to triangulate"). In parallel, a
    tile in progress blocks its eight neighbours, and the wiki says "the
    results are not deterministic" with more than one thread. In the
    point-set route a tile's free borders are first simplified as polylines
    with an error "computed in Z", then the interior. *Gives:* borders that
    agree with no post-pass, and, in the point-set route, a one-dimensional
    vertical-error pass over the border before the interior, which is the
    shape of our seam pass. *Lacks:* a border belongs to whichever tile runs
    first, so the mesh depends on the schedule and the thread count; in the
    greedy route the border is checked only by the first tile's 2D run, and
    the point-set route's interior has no error bound. **Departure:** our
    seam is a function of the seam alone, computed identically by both
    neighbours, so the output does not depend on the order (K5).
  - Cignoni, Ganovelli, Gobbetti, Marton, Ponchio and Scopigno, "BDAM —
    Batched Dynamic Adaptive Meshes for high performance terrain
    visualization", *Computer Graphics Forum* 22(3):505-514, 2003,
    doi:10.1111/1467-8659.00698 (Crossref, **corrected**: authors and issue
    were recalled in the first pass; read,
    `https://vcg.isti.cnr.it/Publications/2003/CGGMPS03a/bdam.pdf`), and the
    same authors' "Planet-sized batched dynamic adaptive meshes (P-BDAM)",
    IEEE Visualization 2003, pp. 147-154, doi:10.1109/VISUAL.2003.1250366
    (Crossref; read, `http://www.crs4.it/vic/data/papers/ieeeviz03-pbdam.pdf`).
    Patches are built bottom up: vertices on the patches' longest edges are
    marked non-modifiable, each square block of four patches is simplified
    by quadric edge collapse to half its vertex count with those vertices
    locked, and the error is measured afterwards as "the maximum vertical
    difference" by rendering both meshes and comparing depth buffers. P-BDAM
    runs the blocks in parallel on five PCs; "synchronization is required
    only at the completion of each bintree level". *Gives:* independent
    simplification of blocks whose shared borders are locked, at planet
    scale. *Lacks:* the target is a vertex count, not a tolerance; the error
    is measured, not bounded, and on a sampled depth buffer; a border is
    locked at the vertices the finer level left, not chosen against the DEM.
  - Hoppe, "Smooth view-dependent level-of-detail control and its
    application to terrain rendering", IEEE Visualization '98, pp. 35-42,
    doi:10.1109/VISUAL.1998.745282 (Crossref; read,
    `https://hhoppe.com/svdlod.pdf`). Blocks are simplified by edge collapse
    until an error threshold, "we constrain ecol's to leave boundary vertices
    untouched", then stitched 2×2 and simplified again; only the last,
    single block simplifies the boundary. The error is the exact L∞ vertical
    deviation from the triangulated grid, found at the vertices of the union
    of the two triangulations: grid points inside faces and grid-line
    crossings inside edges. *Gives:* the exact sup-norm argument our check
    points use (nodes plus row and column crossings, Q14). Hoppe adds that
    for quadtree and bintree subdivisions "all grid line crossings happen to
    fall exactly on grid points", the same reason a grid-line seam's pass is
    exact. *Lacks:* a shared border stays at
    full resolution until both sides are merged (P-BDAM: "some of the borders
    remains not simplified until the whole mesh can be loaded entirely in
    memory"). That is the one alternative to a seam pass, seams at full
    lattice resolution; not taken, because a 2048-node block line would
    carry 2048 vertices per block side at every tolerance.
  - Bertilsson, "Dynamic creation of multi-resolution triangulated irregular
    network", MSc thesis MECS-2015-18, Blekinge Institute of Technology
    (read, `https://www.diva-portal.org/smash/get/diva2:867859/FULLTEXT02.pdf`;
    no DOI). Patches selected independently for rendering; "the border points
    only test their height differences with other border points along the
    same edge. This is required so as to produce identical selection for two
    patches sharing a border", and "computes the same border points multiple
    times to achieve fully independent computation". *Gives:* the closest
    precedent for our seam pass computed by both neighbours. *Lacks:* a
    screen-space heuristic, no error guarantee, no constraints.
  - Zygmunt and Róg, "New approach towards Digital Elevation Model data
    generalisation using the Douglas-Peucker algorithm and Delaunay
    triangulation based on characteristic boundary points", *Measurement*
    260:119849, 2026, doi:10.1016/j.measurement.2025.119849 (Crossref;
    found by `@reviewer`); preprint SSRN 2023, doi:10.2139/ssrn.4639984.
    **Secondary:** the journal is closed (OpenAlex: no open copy),
    ScienceDirect and SSRN both answered 403, and neither record carries an
    abstract, so what follows is from web search summaries of the abstract:
    characteristic points are taken along the dataset borders by
    Douglas-Peucker under a Z-tolerance, so that neighbouring datasets share
    identical boundary points; each dataset is then triangulated by recursive
    Delaunay passes until the Z-tolerance holds; adjacent TINs "can be easily
    joined ... without errors on the boundaries"; sequential and parallel
    runs are reported. *Gives:* on the summary, the seam pass itself (a
    vertical Douglas-Peucker over a shared border, both sides identical,
    before the interior) and a tolerance-driven interior. *Unknown without
    the full text:* whether borders can be constraints or breaklines (seams
    cut along features), whether the tolerance holds between grid nodes or
    only at them, and whether the output is independent of order and thread
    count.
- **Grey literature found by web search** (patents read on their Google
  Patents pages through a summarising fetch; claims not read in full):
  - Starhill et al. (Microsoft), "Maintaining consistent boundaries in
    parallel mesh simplification", US 10,043,309 B2, granted 2018-08-07:
    component meshes of a general 3D model simplified in parallel by edge
    collapse, with a boundary collapse whose result "is independent of the
    data on the interior of the component mesh and enables shared boundaries
    of adjacent component meshes to simplify identically", "without the
    need to synchronize the component meshes". The same principle as our seam pass, for
    general meshes, with no tolerance.
  - Godzaridis and St-Pierre (Bentley Systems), "Multi-resolution tiled 2.5D
    Delaunay triangulation stitching", US 10,255,716 B1, granted 2019-04-09:
    LiDAR terrain tiles triangulated independently, then stitched by removing
    the triangles whose circumcircle crosses a tile boundary and
    retriangulating those points with the neighbours' constraints. Prior art
    for B4 (b), seams removed after the run.
  - HERE `tin-terrain` (code at `b96f3f5`,
    `src/dem2tintiles_workflow.cpp`): groups of web tiles meshed
    independently by a greedy vertical-error method on a crop grown by 100
    cells, then cut to tile bounds; nothing makes two groups agree on their
    shared boundary. Mapbox Martini and Delatin mesh tiles to a maximum
    vertical error independently; their edges do not share vertices, and the
    maintainer's answer is skirts (`mapbox/martini` issue 7, 2019, which
    points to `mapbox/delatin` issue 2). In web terrain the cracks are
    hidden, not prevented.
- **The seam pass's method.** Douglas and Peucker, "Algorithms for the
  reduction of the number of points required to represent a digitized line
  or its caricature", *Cartographica* 10(2):112-122, 1973,
  doi:10.3138/fm57-6770-u75u-7727 (Crossref): insert the worst point, split,
  repeat. The seam pass is that in the vertical (height error at check
  points, not lateral distance), which is greedy insertion (Garland and
  Heckbert, CMU-CS-95-181, 1995, as in increment 14) in one dimension.
- **Pfafstetter coding.** Verdin and Verdin, "A topological system for
  delineation and codification of the Earth's river basins", *J. Hydrology*
  218(1-2):1-12, 1999, doi:10.1016/s0022-1694(99)00011-6 (Crossref). Each
  level divides a basin into nine units (four tributary basins, five
  interbasins), coded by one more digit; BHO's ottocodes are this.
- **Cloud-optimised GeoTIFF.** OGC 21-026, "OGC Cloud Optimized GeoTIFF
  Standard", version 1.0, approved 2023-05-08, published 2023-07-14,
  doi:10.62973/21-026 (read, `https://docs.ogc.org/is/21-026/21-026.html`;
  the year was recalled in the first pass). Headers and offset arrays at the
  start, tiles addressable by range request. The probe above confirms
  ANADEM's copy behaves so.

**Novelty: none claimed.** The combination that could look new: *a greedy
terrain refinement to an exact sup-norm vertical tolerance, split into pieces
along constraints, where each seam is refined beforehand by its own
one-dimensional error and computed identically by both neighbours, so that no
piece waits for another and the output does not depend on the order.* After
the web pass every ingredient has a precedent: separators refined before the
subdomains, with no communication (Linardakis and Chrisochoides); a shared
boundary reduced from its own data alone so both sides agree with no
communication (the Microsoft patent, and Bertilsson's thesis, which computes
it on both sides); tile borders simplified first as polylines by a vertical
error, then the interior (Campos et al., point-set route, first tile owns the
border); the exact sup-norm vertical error over grid points and grid-line
crossings (Hoppe 1998). Closest found: Zygmunt and Róg (2026), which on the
search summary already reduces dataset borders by a vertical Douglas-Peucker
so that neighbours share identical boundary points, then triangulates each
dataset to the Z-tolerance. What remains unconfirmed until its full text is
read: constraints (seams cut along breaklines), the exact sup-norm between
grid nodes, and an output independent of order and thread count. No claim is
made. Before any is made public: the full text of Zygmunt and Róg, a Google
Scholar and Scopus search (web search is not a bibliographic index), the
Campos et al. paper itself (only its code and wiki were read), and Galtier
and George 1996 (not reached).

Searched, 2026-10-01. First pass: Crossref bibliographic queries for every
citation, plus "parallel construction of triangulated irregular network from
DEM", "tile-based TIN generation large DEM seams", "out-of-core terrain
simplification TIN large DEM"; OpenAlex full-text search for "parallel
terrain simplification domain decomposition error bound", "TIN generation DEM
parallel partition seam", "greedy insertion terrain triangulation parallel
tiles", "decoupled constrained Delaunay terrain subdomains separator
refinement", "watershed partition parallel mesh generation hydrological",
"out-of-core triangulated irregular network construction large DEM", "terrain
approximation error guarantee tiles boundary consistency", "Pfafstetter
subbasin parallel mesh". Web pass: one search per closest work (Galtier and
George, Campos et al., BDAM, P-BDAM, Algorithm 872, streaming Delaunay), and
"parallel TIN generation from large DEM tiles seams consistent borders error
bound greedy insertion", "out-of-core terrain simplification TIN tile
boundaries error-bounded massive DEM", "parallel greedy insertion terrain TIN
partition subdomains boundary vertices shared maximum vertical error",
"thesis parallel triangulated irregular network construction large DEM tiles
boundary consistency domain decomposition", "massive DEM TIN construction
parallel blocks block boundary greedy insertion error threshold seamless
merging", "tiled constrained Delaunay terrain mesh tile boundaries fixed
vertices independent tiles TIN LiDAR breaklines out-of-core", "TIN mesh
generation per sub-basin watershed partition parallel hydrological model
shared boundary", "mesh generation each sub-catchment separately shared
boundary polyline", "terrain TIN generation subdomains cut along breaklines
independent refinement no communication vertical tolerance guaranteed on
shared boundary", "Delaunay refinement terrain approximation parallel
subdomain height error separators refined first", "heremaps tin-terrain
zemlya tiles borders", "Delatin OR Martini RTIN tiles max error cracks
skirts", "streaming simplification large meshes processing sequences", "parallel TIN
generation large DEM partition blocks boundary vertices shared error
threshold greedy insertion seamless" (run by `@reviewer`; found Zygmunt and
Róg), and
the patent title above. Found: the works above. The hydrology searches found
the parallel tRIBS of Vivoni, Mascaro, Mniszewski, Fasel, Springer, Ivanov and
Bras, "Real-world hydrologic assessment of a fully-distributed hydrological
model in a parallel computing environment", *J. Hydrology* 409(1-2):483-496,
2011, doi:10.1016/j.jhydrol.2011.08.053 (Crossref; read,
`http://vivoni.asu.edu/pdf/VivoniJH2011.pdf`), which partitions an existing
basin TIN into sub-basins along the channel network to run the simulation, not
to build the mesh; and PIHM pages on mesh constraints, not decomposition.

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
disjoint topologies", `legacy/rasputin/application.py:27`, repeated in `web_visualize.py`)
and `partition()` in `legacy/rasputin/triangulate_dem.h:833`, which splits a finished face
list in two by a predicate (lakes by slope). Neither decomposes a meshing
problem. The second set fetches over the network *at run time*, with no
cache: `legacy/rasputin/avalanche.py:30` calls `requests.get` on NVE's API, and
`wfs_repository.py` reads a WFS service; `web_visualize.py`'s hits are
`localhost` URLs in help text. That is the pattern Ola's direction rules out
(meshing must never wait on the network), so **nothing is carried**. No
legacy constant is re-derived, so no `@migration-expert` pass is needed.

## The subdomain model

### Terms

- **Cut**: a chain carrying the `seam` edge property, a new entry in the
  edge vocabulary (`features.py`). An *artificial* cut is a line the plan
  adds (below, a global grid line); a *natural* cut is a constraint the input
  already has (a BHO unit boundary), which gets the `seam` bit beside its own.
  The core never sees the name: refine is handed `frozen_mask =
  vocabulary.mask("seam")` as a number, as `edge_properties.hpp` requires of
  any policy that treats some lines differently.
- **Piece**: a connected component of the domain's start triangulation across
  every edge that is not a seam. Its boundary is outline and hole edges and
  seam edges. Lakes, roads and other constraints do not divide pieces.
- **Seam edge**: an edge of the start triangulation with the `seam` bit. It
  is shared by two pieces, or by one piece and the outside where a cut runs
  along the outline.
- **Corner**: a vertex where three or more pieces meet, or where a seam meets
  the outline or a feature. Every corner is a vertex of the noded input, so
  every piece around it has it.

### The run, end to end

```
cli.mesh --dem anadem --cache DIR --domain D --out-crs C --tolerance T --block-nodes B
  |
  v  [global, vectors only: no DEM is read here]
plan_basin(request) -> BasinPlan                     (Python, pure; frozen data)
  |  domain, features in the target CRS (15b, 16b)
  |  cuts: grid lines of the global lattice, every B nodes, clipped to the domain  (23c)
  |  build_pslg -> node -> triangulate       the existing engine, run once        (cli._engine)
  |  labels = landcover.regions(triangles, seam edges)   one label per piece      (16c)
  |  per piece: PieceJob(id, start slice, its seam edges, target window, source window)
  v
run_plan(plan, jobs=J, threads=T)  [per piece, in parallel, local]   (async; asyncio.to_thread)
  |  windows: decode_window from cache or local tiles  -> resample (15c D3)   (23a-1, 15c)
  |  seam pass for each of its seam edges: refine_seam(strip, edge, T)        (23b, C++)
  |  start = slice with seam edges split into fans at the seam-pass points    (23c, NumPy)
  |  refine(window, start, frozen_mask=seam)       phase 1                     (14-21, 23b)
  |  refine_points(edge-strip points, frozen_mask)  the edge strip (Q14)       (edge strip)
  |  refine_points(source nodes, frozen_mask)       phase 2, reprojected path  (15c)
  |  write the piece file + its seam record, then free everything
  v
index (conformity check on every seam) -> stitched file, unless --no-stitch   (23c, 23d)
================================ _core boundary ================================
C++ sees: a window RasterView in metres, a start mesh in metres, seam edges as
two points in metres, a tolerance and a frozen mask. No CRS, no path, no URL.
```

### The union argument, checked

The main session's reasoning, item by item.

1. *"A constraint blocks visibility, so the CDT of the union equals the union
   of the subdomain CDTs when the shared constraint's vertices agree."*
   **Correct, with one word changed: "a" CDT, not "the".** By Chew's
   definition, triangle t of piece A is constrained Delaunay if no vertex
   visible from t's interior lies strictly inside its circumcircle. Take a
   vertex w of piece B and a point x inside t. The segment from x to w
   leaves A through its boundary. Either it crosses a seam edge's interior,
   and w is not visible; or it passes through a seam vertex v. In the second
   case, if v is a vertex of t or lies on t's circumcircle, then x is inside
   the circle and v on it, so the ray beyond v is outside the circle and w is
   not inside; if v is strictly inside the circle, v itself is visible from x
   inside A (else w is not visible either) and t was not constrained Delaunay
   in A. So the union of the pieces' CDTs is a CDT of the union, with the
   seams as constraints. On a DEM lattice the CDT is not unique (cocircular
   ties are common, 21 §4: 9 % of incircle tests are exact ties), so the union
   is *a* CDT, not necessarily the one a whole-domain run would pick. This
   design sidesteps the difference: every piece's start mesh is a slice of
   one start triangulation of the whole domain.
2. *"Insertions strictly inside a subdomain never affect the other."*
   **Correct.** `split_inside` writes one triangle's slot and two new ones;
   Lawson flips never cross a constrained edge (14b R5); the quality pass's
   walk stops at one (`include/terrain/mesh/quality.hpp@98da562:143`); a
   triangle's scan reads only the nodes inside it.
3. *"Only splits on a shared constraint need the two neighbours to agree."*
   **Correct, and there are five ways a vertex lands on a constraint in
   today's code and the designed increments:** refine's worst node lying
   exactly on a constrained edge
   (`include/terrain/refinement/refine.hpp@98da562:335`, split at `:356`);
   20b's feet (`include/terrain/refinement/refine.hpp@98da562:339`); the quality pass's node on an edge
   (20 R6); 15c's `refine_points` with a check point on an edge; and the
   edge strip's check points, which lie on constraints by design. This
   design turns off all five on seams (frozen edges) and gives the seam its
   vertices before either piece runs (the seam pass).
4. *"Corners where 3+ subdomains meet are fixed vertices."* **Correct.**
   They are vertices of the noded input, so they are in every start slice
   around them, and a frozen edge is never split, so they stay corners.
5. *"Natural cuts add nothing; artificial cuts add constraints."* **Half
   correct.** A natural cut adds no edge, but freezing it changes what refine
   may do there: no feet and no splits on it except the seam pass's. So the
   mesh along a natural cut still differs from a whole-domain run. An
   artificial cut adds edges, and they stay in the output (B4).
6. *"The cut and the split exchange must be deterministic functions of the
   input."* **Agreed, and one more input must be excluded: the machine.** A
   cut chosen from physical memory or the thread count would make the mesh
   depend on the computer. The cut is a function of the domain, the lattice
   and `--block-nodes`, which is recorded in the file (K5).
7. *"Shared-edge splits are common, given Q14's crossing/midpoint check points
   along constraints."* **Correct under Q14 as designed** (check points
   inserted after refine, by `refine_points`). Here they are moved *before*
   refine on seams, into the seam pass, which is what removes the exchange.

## Choosing the cuts

### Artificial cuts: global grid lines (23c)

- **The lattice.** Node `(R, K)` of the computation grid: 15c's target grid
  on the reprojected path (`x = K·h`, `y = −R·h`, J6), or the mosaic's own
  lattice from 15a's reference node on a projected DEM meshed directly.
- **Blocks.** Block `(i, j)` holds the nodes with `K` in `[iB, (i+1)B]` and
  `R` in `[jB, (j+1)B]`, for `B = --block-nodes` (default 2048, B3). The block
  lines are the same for every request on that lattice: they are tiles on
  file, in Ola's sense, not a property of one domain.
- **When to cut (B2).** If the needed window (the domain grown by the cell
  diagonal, 15b) spans more than `B` nodes in either direction, every block
  line crossing the domain's interior becomes a cut. Otherwise nothing is
  cut, there is one piece and no seam, and the mesh is today's, bit for bit
  (K1).
- **The chains.** Each line, from the window's edge to its edge, is clipped
  to the domain as linework (16b R6) and enters as breakline chains with mask
  `seam`. The noder nodes them like any input (crossings with the outline,
  features and the other lines). A block corner is where two cut lines
  cross, at a node whose coordinates are whole metres when `h` is, so on the
  noder's 1 mm grid; that the noder's crossing arithmetic returns it exactly
  for two axis-parallel segments is for 23c's red suite to pin (DC9), not
  assumed here.
- **Why grid lines** and not, say, a balanced bisection:
  1. *No quality cost.* A node not on a grid-line seam is at least one cell
     from it, and 20b's ε is at most half a cell (`foot_epsilon`'s cap), so a
     foot on a grid-line seam never triggers: freezing it changes nothing a
     foot would have done.
  2. *An exact seam pass.* The nodes on a grid line lie exactly on the seam,
     and the bilinear surface along a grid line is piecewise linear between
     them, so after the seam pass the tolerance holds at **every point** of
     the seam against the DEM's bilinear surface, not only at check points.
  3. *Exact heights.* A seam vertex on a grid line is a node, so its z is the
     node's value in both pieces, bit for bit.
  4. *Locality.* Blocks of a fixed global grid are the same whatever the
     request, so two catchments meshed separately share their block lines.
- **Several pieces in one block.** Where the domain enters a block twice, the
  labelling gives two components and so two pieces.
- **Piece ids.** `(j, i, k)`: block row, block column, and the component's
  rank by its lowest start-triangle index. A function of the input alone.

### Natural cuts: BHO's Pfafstetter levels (later, 23e; B1)

- **The descent.** Start from the domain's own unit (ottobasin 76 for the
  basin). Replace any unit whose window exceeds `B` nodes in either direction
  by its children (one more ottocode digit). Stop at the elementary
  catchments, and cut any still-oversized one by grid lines. Deterministic: a
  function of BHO and `B`. The pieces are then hydrological units, and the
  output keeps that: one file per unit, named by its ottocode.
- **What it needs, measured:** BHO is an exact coverage (every shared edge
  twice, bit for bit), so the noder merges each shared boundary into one
  seam, as it does for CORINE (16b R7).
- **What it costs, from the measurement:** Pfafstetter units are unbalanced
  (56 to 5,201 km² at level 6 in the piece), which the descent handles; and
  BHO's boundaries carry a vertex every ~107 m, used as given (Ola's input
  model). At level 7 the piece's 1,768 km of seams would bring about 16,400
  vertices, against 44,914 mesh vertices at 50 m (+37 %) and 407,376 at 10 m
  (+4 %): arithmetic on the measured lengths and the basin-piece mesh sizes,
  not a run. On general (non-grid) seams frozen feet can leave needles next
  to the seam; to be measured. B10 asks whether shared chains should be
  simplified, once per chain so both sides agree, with increment 22's
  area-preserving collapse.

## The seam protocol

**Recommended (B5 (a)): seams are frozen, after a one-dimensional seam pass
that each neighbour computes identically.** No exchange, no rounds.

1. **One start triangulation, vectors only.** The outline, holes, features
   and cuts go through the existing engine once (`build_pslg` → `node` →
   `triangulate`; `triangulate` takes only a noded PSLG,
   `include/terrain/cdt/triangulate.hpp@98da562:23`). The noder runs once, so
   each seam has one noded geometry, shared by both sides by construction.
   `landcover.regions(triangles, seam_edges)`
   (`src_python/tin_engine/landcover.py@98da562:42`, 16c's components across
   unblocked edges) with only seam edges blocking gives the pieces. This is
   the one global step, and it never reads the DEM (K9).
2. **The seam pass, per seam edge** (`refine_seam`, C++, 23b). For an edge
   `(a, b)`:
   - *Check points* on the open edge: every lattice node exactly on it (exact
     orientation), every crossing with a row or column line (Q14), and the
     midpoint between each two neighbouring crossings, the edge's ends counting
     as neighbours (Q14 as ruled, with the ends as 15c proposed). z is
     bilinear at each; a check point with a NoData stencil is skipped and
     counted. On a grid-line seam the crossings are the nodes and the
     midpoints add nothing (the surface is linear between nodes), so the
     check points are the nodes.
   - *Greedy, one-dimensional.* While some check point `p` has
     `|z_p − lerp(p)| > tolerance`, where `lerp` interpolates between the
     current vertices on either side of `p`, insert the worst one (ties to the
     smallest parameter from `a`), and repeat. Douglas-Peucker in the
     vertical; it ends because each check point is inserted at most once.
   - *Output:* the inserted points in order from `a` to `b`, world `(x, y)`
     and z, plus z for `a` and `b` themselves.
   - *Who computes it: both neighbours, each for its own seam edges.* The
     inputs are the edge's world endpoints (from the one start triangulation),
     the tolerance, and a DEM strip that is a function of the edge alone (its
     bounding box snapped outward to the global lattice and grown by one
     node), with the edge oriented so `a < b` by `(x, y)`. The strip's values
     are a function of each node, not of the window (15c J6 and D3: each
     target node resampled from its own coordinates; 15 R5: overlaps resolved
     independently of order). So both pieces get the same points, bit for
     bit. This keeps every piece job self-contained: no global seam phase
     between planning and pieces. The price is computing each seam twice,
     which is one-dimensional work.
3. **The start slice, with fans** (Python, NumPy). A piece's start mesh is the
   start triangulation's triangles with its label, vertices renumbered in
   first-use order, constraint edges and masks restricted to them. Each seam
   edge `(a, b)` with seam-pass points `p1 … pk`, in triangle `(a, b, c)`,
   becomes `(a, p1, c), (p1, p2, c), …, (pk, b, c)`: index arithmetic, no
   geometry. Each fan triangle is counter-clockwise when every `pi` lies
   between `a` and `b` and `c` is off the edge's line by more than rounding;
   the points lie on the edge to rounding, and no vertex's 1 mm snap cell
   meets an edge it is not an endpoint of (the noder's guarantee 14(b),
   `05-noder.md`), which keeps `c` half a millimetre from the segment, not
   from its line; on a grid-line seam the `pi` are exact nodes on an
   axis-parallel line, so the fan is valid whenever `(a, b, c)` is; on a
   general seam a flat triangle can leave `c` within rounding of the line.
   If it is not valid, refine's own check refuses
   it (`NotCounterClockwise`): loud, never a silent difference between two
   pieces.
4. **Refine with the seam frozen** (23b). `RefineOptions` gains
   `std::uint32_t frozen_mask = 0`; an edge is frozen when its mask meets it.
   `refine` legalises the fans (its first `legalise_all`) and refines as
   today, except that **nothing inserts a vertex on a frozen edge**:
   - the scan does not count a node lying exactly on a frozen edge of the
     triangle (it is the seam pass's); only triangles that have a frozen edge
     take this path, so every other triangle is scanned as today;
   - no foot is taken on a frozen edge; the node itself goes in, as on 20b's
     refused-foot path;
   - the quality pass skips a node that lies on a frozen edge
     (`skipped_frozen`, a new counter);
   - 15c's `refine_points` (with `PointRefineOptions::frozen_mask`) skips a
     check point on a frozen edge and counts it (`on_frozen`, with its error,
     as `coincident` is counted);
   - the edge strip makes no check points on frozen edges;
   - `LatticeMesh::split_edge` asserts that its edge is not frozen, so a
     sixth path added later cannot do it silently in a debug build (the
     assertion is untested, FE6).
5. **Heights on seams.** Each piece writes its seam vertices' z from the seam
   record (the seam pass's output), not from its own bilinear evaluation, so
   both pieces write the same z bit for bit (K4). Inside refine the planes at
   a seam vertex use its own `vertex_z`, which can differ from the record in
   the last bits when the vertex is off-node (a crossing with the outline or a
   feature; never on a grid-line seam's own vertices). The guarantee then
   holds to that rounding, as it already does for any off-node start vertex
   today (`refine` writes `raster::bilinear` at the world point and scans with
   `vertex_z` at the lattice point); the oracles' 1e-9 slack covers it.

**Termination.** The seam pass ends per edge (finite check points, each
inserted once). Refine ends as before: every insertion is a valid node not
yet a vertex, and freezing only removes candidates. There are no rounds
between pieces, so nothing else can fail to end.

**Determinism.** The plan is a function of the input, the lattice and
`--block-nodes`; the noder and the CDT are deterministic; piece ids come from
the input; the seam pass is pure per edge; the slice keeps the start
triangulation's order; each `refine` is deterministic for any thread count
(14 R5; 21's L1 if 21d lands). So the output does not depend on `--jobs`,
`--threads`, which piece runs first, or the machine (K5).

**The alternative (B5 (b)): exchange on shared segments**, the protocol of
Chernikov and Chrisochoides' PCDM and of Ola's description. Pieces refine with
seams splittable; each split a piece makes on a seam (a foot, a node on it, a
phase-2 point, an edge-strip point) is recorded as a request; after a pass,
each seam's requests from both sides are merged (sorted, deduplicated) and
applied to both pieces, which then continue; repeat until a pass makes no
request. It ends, because feet are taken once per node and check points are
finite, and it is deterministic if each pass is. **What it buys:** feet and
splits on natural seams, so no quality loss next to them. **What it costs:**
rounds in which both neighbours must be resumed or rerun (so either both are
live at once, against locality, or a piece is rerun per round), a merge
step, and a second kind of refine run (resume with injected splits). On grid
seams it buys nothing, since no foot triggers there. So (a) now; (b) only if
the natural-cut acceptance shows needles along BHO seams.

## Tolerance, the final check and the constraint check points per piece

- **The guarantee, decomposed** (K3). Every valid DEM node inside the domain
  is within `--tolerance` of the mesh: a node strictly inside a piece by that
  piece's refine, as today; a node on a seam by the seam pass. On a grid-line
  seam the tolerance holds at every point of the seam against the bilinear
  surface; on any other seam, at its check points.
- **The edge strip (Q14), per piece.** Its check points on a piece's
  non-frozen constraints (outline pieces, features) are made and refined in
  that piece, after `refine`, as the edge strip is designed; its points on
  seams are the seam pass's. The wording Ola ruled for its guarantee holds
  over the union.
- **The final check (15c phase 2), per piece.** A piece's check points are
  the source nodes whose projection falls in its target window. Points
  outside the piece's triangles are never scanned (15c D4), so each source
  node is checked by the one piece whose closed triangles hold it, except a
  node exactly on a seam, which both pieces skip and count in `on_frozen`
  with its error. For a geographic source, projected nodes landing exactly on
  a seam are not expected; the count is there so it cannot happen silently.
  J2 holds over the union with that one named exception.
- **The seam pass on the reprojected path** reads the resampled target grid,
  as the edge strip does. The guarantee against the source DEM is phase 2's,
  at source nodes, which are not on seams.
- **The memory cap** (15a R7, half of physical memory) applies to what runs
  at once: the runner starts a piece only while the running pieces' window,
  store and mesh estimates fit under it. It can delay a piece, never change
  one.

## Windowed source reads

- **The seam is 15c's `SourceWindows`** (`window(r0, r1, c0, c1)` and `meta`),
  which 15c built so that windowed decoding plugs in as a change of provider.
- **23a-1 adds a block layer under it** (`io/cog.py`):
  - `BlockSource` (protocol): `page` (the parsed full-resolution TIFF page:
    shape, block shape, dtype, compression, offsets and byte counts) and
    `block(index) -> bytes`. Two implementations: `LocalTiffBlocks`, a local
    tiled GeoTIFF, reading `page.dataoffsets` ranges from the file (15 R3
    [15d], DTM10 and ANADEM's MGRS tiles are 512² tiles); and
    `CachedBlocks`, the cache's block files.
  - `decode_window(blocks, window) -> DemTile`: decodes only the blocks that
    meet the window, with `TiffPage.decode` (shown above to decode a block
    from its bytes alone), on a thread pool (zlib releases the GIL, 15 R3),
    into one array that becomes the tile through `DemTile._adopt` (its third
    caller, under the same rule: allocated here, never handed out writable).
    It equals decoding the whole tile and slicing (15's T-window).
  - `BlockWindows(SourceWindows)`: for a source window, every object (file or
    cached COG) meeting it, decoded by window and assembled by 15a's overlap
    rule (R5, order-independent).
- **No dense canvas at basin scale.** Each piece holds a target window (its
  start slice's bounding box grown by the cell diagonal, snapped outward to
  the lattice) and a source window (15c D2's `source_region` of that target
  window: its image in the source CRS, grown by two source cells). Both are
  about `B²` nodes plus margins. Nothing holds the basin.
- **Header before pixels** (15c J10). Coverage and missing blocks are known
  from the parsed headers and the cache's presence check before any block is
  decoded.
- **A projected DEM meshed directly** (Norway) gets the same: a piece's window
  of the 15a mosaic is `plan_mosaic(footprints, piece bounds)` (pure, 15
  R12), each tile decoded by window.

## The fetch step and the tile cache

### Types (`fetch/sources.py`, frozen Pydantic)

```
RemoteSource      id ("anadem-v1", "glo30"); kind ("one-cog" | "cog-tiles");
                  url, or url_template plus tile_list_url; crs (expected, checked
                  against each header); nodata; credit; licence_note
SOURCES           the catalogue: a Mapping[str, RemoteSource] of data, two entries
FetchRequest      source id, domain (path and CRS), target CRS or None, margin
                  (source cells, default 2), cache root
FetchPlan         source id; objects: (object id, url, block indices, bytes)
CacheManifest     source id, url, content_length, last_modified, header_sha256,
                  header_bytes, crs, block shape, rasputin version, and the
                  requests fetched (domain hash, date)
```

`anadem-v1` is OpenTopography's COG
(`https://opentopography.s3.sdsc.edu/raster/ANADEM/ANADEM_be/anadem_v1_compressed_COG.tif`,
DOI 10.5069/G9736P4G, as ruled in Q17). `glo30` is the AWS bucket, one COG per
1° tile, with `tileList.txt` saying which tiles exist; a tile absent from the
list is sea, recorded as "no tile", not as missing.

### Planning (`fetch/plan.py`, pure)

- **The header** of each object: read a prefix (1 MiB, doubling) until
  every page's offsets number its blocks; tifffile does not raise on a
  short prefix, it logs and returns a page without offsets (measured at
  64 KiB and 1 MiB). ANADEM needs at most 8 MiB (measured). The prefix is
  cached as `header.bin`.
- **The needed region** in the source CRS: the same function 15c's
  `source_region` uses (the target window's image, grown by two source
  cells), applied to the whole domain. Fetch and mesh call this one function,
  so for the same domain and options the blocks meshing needs are a subset of
  the blocks fetch got (K7).
- **Blocks**: the full-resolution page's blocks whose footprint meets the
  region (prepared shapely). Overviews are not fetched.

### The cache (`io/repository.py`, the one module in `io/` that opens files, 15 Q3)

```
<cache>/<source-id>/manifest.json                       identity, written atomically
<cache>/<source-id>/<object-id>/header.bin              the parsed prefix
<cache>/<source-id>/<object-id>/blocks/<row>/<col>.bin  one block, its exact bytes
```

- **The directory is the inventory; the manifest is the identity.** A block is
  present exactly when its file exists with its byte count from the header.
  Blocks are written as `*.part` and renamed, so a crash never leaves a half
  block, and nothing has to be kept in step with the files. The manifest says
  what the cache is a copy of, and what was asked of it.
- **Identity.** Each fetch run sends one `HEAD` and compares
  `Content-Length` and `Last-Modified` with the manifest, and the header
  prefix's sha256; any change refuses ("the remote copy changed;
  `--refresh` re-fetches it"). The ETag is not used: ANADEM's is a
  placeholder (measured). Meshing never checks the remote; it trusts the
  manifest, because it is offline by rule.

### Downloading (`fetch/http.py`, `fetch/run.py`)

- Standard library `urllib.request` with `Range`; no new dependency.
- A response must be 206 with a matching `Content-Range`, and every block's
  length must equal its byte count before the rename.
- Runs of wanted blocks whose byte gaps are at most 64 KiB are coalesced into
  one request of at most 8 MiB (ANADEM's blocks are in row-major order with
  small gaps between them).
- `asyncio` with `asyncio.to_thread` per request, at most `--connections`
  (default 8) at once; 5xx and timeouts retried three times with backoff,
  4xx never.
- **Resumable:** present blocks are never requested; `.part` files are
  removed on start. **Incremental:** a second domain fetches only its missing
  blocks. `--dry-run` prints the plan (objects, blocks present and missing,
  bytes). For the basin on ANADEM that is 3,061 blocks, 1.71 GiB (measured).

### The CLI

```
rasputin fetch anadem-v1 --domain basin.geojson --out-crs EPSG:31983 --cache DIR [--dry-run] [--refresh]
rasputin mesh --dem anadem-v1 --cache DIR --domain basin.geojson --out-crs EPSG:31983 --tolerance 5 --out basin.vtk
```

With `--cache`, `--dem` names a catalogue source; without it, a path, as
today. A block meshing needs and the cache lacks is refused before any
decode, naming the command: "anadem-v1: 37 of the 3,061 blocks this domain
needs are not in DIR; run: rasputin fetch anadem-v1 --domain basin.geojson
--out-crs EPSG:31983 --cache DIR". The cache root comes from `--cache` or
`RASPUTIN_CACHE` (B7).

### The I/O boundary

Network code lives only in `fetch/` (`http.py` is the one importer of
`urllib.request`); cache files are opened only in `io/repository.py`; CRS
stays in Python (the manifest records the header's CRS, checked against the
catalogue's); C++ sees decoded windows as numbers. Two checks make it
testable: no `tin_engine` module on `rasputin mesh`'s path imports
`tin_engine.fetch` or `urllib.request` in its own source (an AST check over
`src_python/tin_engine`; `sys.modules` cannot be the oracle, because pyproj
imports `urllib.request` itself), with `cli.py` importing `tin_engine.fetch`
lazily inside the fetch command; and a mesh run with `socket.socket`
replaced by one that raises (K7). The fetch step records the source's credit
in the manifest, and the mesh file's `elevation_source` carries it.

## Memory and parallelism

**Per piece**, at `B = 2048` (a 61.4 km block at 30 m, 3,775 km²): target
window and source window about 16 MiB each (4 B per node), check points
about 64 MiB (16 B per source node, 15c D4), and the mesh at about 310 B per
triangle (the basin-piece fit, everything in one process included).
Triangles per block, by the basin-piece densities (all arithmetic, not
measured):

| tolerance | basin mean | basin p90 | the steep piece | memory per piece (mean to piece) |
|---:|---:|---:|---:|---:|
| 1 m | 1.41 M | 2.86 M | 3.38 M | 0.5-1.1 GiB |
| 5 m | 0.21 M | 0.53 M | 0.64 M | 0.2-0.3 GiB |
| 10 m | 0.08 M | 0.20 M | 0.26 M | ~0.1 GiB |

So eight pieces at once at 1 m stay under about 9 GiB, and memory follows
`--jobs`, not the basin.

**Two levels of parallelism.** `--jobs J` pieces at once, each refining with
`T` threads (`J × T` at most the core count; the defaults split the cores
evenly). Pieces run in one process on `asyncio.to_thread`: refine releases the
GIL, and so do zlib and most of pyproj and NumPy. If the GIL shows in the
profile, the follow-up is a process pool; pieces are already pure data in and
out. Pieces start largest window first, so the biggest one is not last.

**Does this lift increment 21's serial-phase ceiling? For the basin, yes.**
Today refine reaches 2.05× at 10 threads on the Velhas piece at 1 m, with the
serial split and flip at 77 % (basin-piece README). Pieces run their serial
phases at the same time, so the serial phase is parallel across pieces. On
the Velhas piece (11 blocks at `B = 2048`), if the largest block holds 15-20 %
of the work (not measured), its single-thread refine would be about 2-2.6 s
against 6.32 s for the whole piece at 10 threads today: about 2.5-3×. The
basin's box is 2.31 G nodes, 551 blocks of 4.19 M nodes, and the basin covers
about a third of it, so on the order of 200 pieces: enough to keep every core
busy. The limits then become cores, memory bandwidth and the global vector
step, not the serial phase. Below the block size (one piece, such as the
Norway benchmark with `--block-nodes 0`) nothing changes, and 21d stays the
route for one piece's serial phase. All of this is arithmetic for `@perf` to
measure.

**The basin at 1 m, as arithmetic:** 237 M triangles (215-261 M, the box
sample) at the piece's 12.95 s per 10.43 M triangles on one thread is about
5 min of refine on one core, about 40 s on 8 if it scales; resampling 734 M
source nodes at the prototype's 30.9 M per 0.99 s on 8 threads is about 24 s;
plus phase 2, decoding 3,061 blocks and writing.

## Output

- **Pieces and an index, always, when a run is cut.** Each piece is written
  as soon as it is done, as a normal `.vtk` or `.ply` (whichever `--out`
  names), to `<out>.pieces/<id>.<ext>`, with a seam record beside it: for
  each of its seam edges (by its index in the start triangulation), the
  vertex sequence with `(x, y, z)`. The index, `<out>.pieces/index.json`
  (`MeshIndex`, frozen Pydantic, `io/mesh_index.py`), holds the CRS, the
  tolerance, the decomposition (`block_nodes`, the lattice, the cut kind),
  the source's identity and credit, and per piece its file, sha256, counts
  and window.
- **Conformity is checked when the index is written** (K4): for every seam
  edge, the two pieces' records must be equal bit for bit; a difference fails
  the run naming the seam and the two pieces.
- **One stitched file by default** (B6 (a)): `--out basin.vtk` still means one
  file. The stitcher streams piece by piece (memory: one piece and the seam
  vertices), numbering vertices by first occurrence in piece-id order, so a
  seam vertex gets one number. `--no-stitch` skips it; `rasputin stitch
  <out>.pieces` does it later. The legacy `.vtk` and the PLY both need counts
  in their headers, which the index has before the first byte is written.
- **How a consumer reads it.** The stitched file reads as today. A consumer of
  pieces reads `index.json` and opens the files it wants; each is a complete
  mesh of its piece. Seam edges carry the `seam` bit, named in the file's
  vocabulary, so a consumer can tell a cut from a river.
- **Resumable runs** (23d; decided here, Ola may overrule, B9): each piece's
  job has a hash of its spec (input identities, options, rasputin version),
  stored beside its file; a rerun skips pieces whose hash matches.

## What survives of 15c, the edge strip and 15d

- **15c survives whole**, and at basin scale runs per piece. `TargetGrid` (a
  sub-rectangle of one global lattice, J6) is the piece window; `resample`,
  `check_point_blocks`, `CheckPoints` and `refine_points` are used unchanged
  except for the frozen mask (23b). D8 and the Q11-Q17 rulings stand. Its
  memory cap applies per piece. 15c's acceptance stays on the Velhas piece,
  undecomposed. What 15c already marked "Superseded at basin scale" is
  replaced by this record and nothing else is. If 23a lands before 15c-2,
  15c-2's acceptance reads ANADEM from the cache instead of a one-off cut.
- **The edge strip survives as designed** (15c, "The edge strip", Q14), with
  two additions: its check-point generator is written as one function
  (`constraint_check_points`, in C++) that the seam pass reuses, and it makes
  no check points on frozen edges.
- **Under B12 (a), 15d's window decoding survives** (15 R3 [15d] and its T-window test) and
  moves into 23a-1 as `decode_window` over a `BlockSource`. **15d's basin
  memory plan does not** (15 R7 [15d]: the basin "meshable in one piece with
  a dense canvas", 8.6 GiB): superseded by pieces. 15d's whole-basin run moves
  to this record's acceptance. ROADMAP item 2.3 becomes 23a-1.
- **15 R12** is met as written (one frame, global node identity, a subdomain's
  raster as a plan, order-independent overlaps, vectors transformed once, an
  explicit output CRS). One of its rejections no longer applies: "resampling
  per window (neighbouring windows would disagree at seams)". Under 15c J6 and
  D3 every target node is computed from its own coordinates, so two windows
  agree at every shared node, bit for bit.
- **21d stays deferred.** Pieces give the basin its parallelism (21's option
  D, which Ola's Q5 placed with the large-area work); 21d remains the route
  for the serial phase inside one piece.

## Invariants

- **K1. Small runs untouched.** With no cut (the window within
  `--block-nodes`, or `--block-nodes 0`) and `frozen_mask` 0, every mesh is
  bit-identical to master's: 23b and 23c change nothing on that path.
- **K2. Frozen means frozen.** In a piece, no pass inserts a vertex on a
  frozen edge: refine's split, feet, the quality pass, `refine_points` (phase
  2 and the edge strip). Pinned by FE2-FE5; `split_edge`'s assertion is a
  debug guard against a later path, untested (FE6).
- **K3. The guarantee over the union.** Every valid DEM node inside the
  domain is within `--tolerance`: inside a piece by refine, on a seam by the
  seam pass; on a grid-line seam at every point against the bilinear
  surface. On the reprojected path, every source node too, except those
  exactly on a seam, which are counted (`on_frozen`) with their error.
- **K4. Conformity.** Two pieces sharing a seam edge write the same vertex
  sequence on it, `(x, y, z)` bit for bit; checked when the index is
  written, and a difference fails the run.
- **K5. Determinism.** The output is a function of the inputs and the
  options (`--block-nodes` included), never of `--jobs`, `--threads`, the
  order pieces run in, or the machine's memory.
- **K6. Locality.** A piece's files depend only on its job spec: its start
  slice, its windows and its seam edges. Re-running one piece alone gives the
  same bytes.
- **K7. Meshing is offline.** No network on `rasputin mesh`'s path; blocks
  meshing needs that the cache lacks are refused before any decode, naming
  the fetch command; for the same domain and options, mesh needs a subset of
  what fetch got.
- **K8. The cache is honest.** A block is present exactly when its file has
  its byte count; a fetch refuses a remote whose length, date or header
  changed.
- **K9. One global step, vectors only.** Noding, the start triangulation and
  the labelling run once for the whole domain and never read the DEM; their
  memory follows the vector input.
- **K10. The I/O boundary.** The core sees windows, meshes, edges, a
  tolerance and a mask, all as numbers in metres: no CRS, path or URL.

## Degeneracy policy

- **Nodes on a seam** (every node on a grid-line seam): the seam pass's
  check points, never a piece's scan candidates.
- **A seam along the outline** (a block line on a `--bbox` edge, say): not a
  cut, since it crosses no interior. A natural cut that coincides with an
  outline or feature edge is merged by the noder (masks unioned), frozen, and
  handled by the seam pass like any seam.
- **A seam along a feature edge** (a road on a block line): one edge with both
  bits; the feature keeps its bits, and the seam's rules apply.
- **A piece with no DEM node** (a sliver between a block line and the
  outline): refine has nothing to scan; the seam pass still runs on its seam
  edges, so its boundary still meets the tolerance there.
- **A seam edge shorter than a cell** may have no check point; its two
  vertices are its only heights, as any short constraint edge's are.
- **NoData on a seam:** a check point whose stencil touches NoData is skipped
  and counted; a seam vertex without a height is invalid, and the triangles
  around it are carved by 14's rule in each piece.
- **A fan that is not counter-clockwise** (a vertex within rounding of a seam
  edge's line, which the noder's guarantee 14(b) rules out): refused by
  refine's `NotCounterClockwise`, with the piece and seam named.
- **Several components of the domain in one block:** several pieces.
- **Four pieces at a block corner inside the domain:** the corner is a node
  and a vertex of all four.
- **A source node exactly on a seam** (phase 2): skipped and counted in both
  pieces, `on_frozen`, with its error.
- **The remote changed** between fetches: refused; `--refresh` re-fetches.
- **A GLO-30 tile not in `tileList.txt`:** sea, "no tile"; a domain whose
  needed region falls there gets NoData there, as 15a's coverage rule says
  (refused if 15a would refuse it).

## Not in scope

- **Natural cuts** (BHO Pfafstetter units, 23e) and the simplification of
  their shared chains (B10): designed above in outline, built later with the
  basin's own inputs (ROADMAP item 2.5).
- **Removing artificial seams after the run** (B4 (b)).
- **The exchange protocol** (B5 (b)).
- **Noding in pieces.** The global vector step is fine at basin scale with
  outline, cuts and BHO rivers; a polygonised MapBiomas might not be, and
  would need its own look.
- **A process pool** for pieces, unless the GIL shows in the profile.
- **Fetching other datasets** (BHO, MapBiomas): the catalogue is data, so
  adding a COG source is an entry, but vector sources need their own fetch.
- **Machines other than one:** pieces are pure data in and out, which keeps a
  cluster possible; nothing here builds it.

## Order of work and PR split

### The proposed order (B12)

ROADMAP's "Order of work from 2026-09-30", item 2, today reads: measure a
piece (done), 15c, 15d, parallel refine or domain decomposition, the basin's
own inputs. Proposed in its place:

1. **15c-1**, as designed (`15c-geographic-dem.md`): check points and the
   final check in C++. Needs nothing new; ready for `@tester`.
2. **23a-1**, windowed decoding and the cache's read side. Replaces 15d.
3. **23a-2**, the fetch step. With 23a-1, rasputin fetches and caches its own
   ANADEM, so 15c-2's acceptance no longer needs a one-off cut while
   `metadados.snirh.gov.br` answers 403.
4. **15c-2**, as designed, its acceptance on ANADEM from the cache.
5. **The edge strip** (Q14), as designed, writing `constraint_check_points`.
6. **23b**, frozen edges and the seam pass, in C++.
7. **23c**, the decomposition, piece by piece, with pieces and the index as
   output.
8. **23d**, pieces in parallel, resumable runs, `rasputin stitch`.
9. **The basin run**, `@perf`, no code: 50, 20, 10, 5, 2 and 1 m.
10. **23e and the basin's own inputs**: natural cuts, BHO rivers as
    polylines, MapBiomas land cover.

**Against a ruling.** This reverses part of Q14's placement ruling (the
edge strip right after 15c and before 15d, `15c-geographic-dem.md`, Q14 and
"The edge strip"): under this order 23a-1, which replaces 15d, and 23a-2 land
between 15c-1 and 15c-2, and the edge strip after both. B12 (b) keeps the
ruling.

**Why this order.** 23a has no dependency on the core and gives every later
step real ANADEM data, which is what Ola asked the fetch step for ("or our
performance will drop while we wait for download"). 23b needs `refine_points`
(15c-1) and the edge strip's generator. 23c needs 15c-2's `TargetGrid` and
`resample` per window, and 23b. 23d is what makes the basin fast, but
nothing in it changes a mesh, so it goes after the first correct cut runs.
Release hardening and 21d stay deferred, as ROADMAP has them.

### PR split and LOC

Counted in `CLAUDE.md` §2's unit. Estimates; the worst case applies 39 %
(increment 10's overrun), with 60 % (15a's `mosaic.py`) beside it. Every PR
stays under 700 at both.

| PR | what | est. | +39 % | +60 % |
|---|---|---:|---:|---:|
| **23a-1** | **Windowed decoding and the cache, reading** | | | |
| | `io/cog.py`: `BlockSource`, `LocalTiffBlocks`, `decode_window` on a thread pool | 60 | | |
| | `io/repository.py`: `CacheRepository` read side, `CacheManifest`, presence | 55 | | |
| | `BlockWindows` (15c's `SourceWindows`), the missing-block refusal | 40 | | |
| | `fetch/sources.py`: `RemoteSource`, the catalogue (two entries) | 40 | | |
| | `cli.py`: `--cache`, `RASPUTIN_CACHE`, `--dem <source id>` | 25 | | |
| | **23a-1 total** | **220** | **306** | **352** |
| **23a-2** | **The fetch step** | | | |
| | `fetch/plan.py`: header prefix, needed region, blocks, GLO-30 tile list | 75 | | |
| | `fetch/http.py`: ranged GET, 206 and length checks, coalescing, retries, async bound | 80 | | |
| | `fetch/run.py`: plan, missing, download, verify, put; `--dry-run` | 50 | | |
| | `io/repository.py`: write side, atomic put, manifest, identity check | 35 | | |
| | `cli.py`: `rasputin fetch` | 45 | | |
| | **23a-2 total** | **285** | **396** | **456** |
| **23b** | **Frozen edges and the seam pass, C++** | | | |
| | `lattice_mesh.hpp`: the frozen mask, `is_frozen`, the assertion in `split_edge` | 15 | | |
| | `scan.hpp`: nodes on a frozen edge skipped (only for triangles with one) | 30 | | |
| | `refine.hpp`: `frozen_mask`, no feet on frozen edges | 15 | | |
| | `quality.hpp`: skip on frozen, `skipped_frozen` | 10 | | |
| | `refine_points.hpp`: skip on frozen, `on_frozen` | 15 | | |
| | `seam.hpp`: `refine_seam` (one-dimensional greedy over `constraint_check_points`) | 100 | | |
| | `bindings/core.cpp`, `_core.pyi` | 60 | | |
| | **23b total** | **245** | **341** | **392** |
| **23c** | **The decomposition, piece by piece** | | | |
| | `decompose.py`: when to cut, the grid lines as chains, `BasinPlan` | 75 | | |
| | `features.py`: the `seam` property | 5 | | |
| | `pieces.py`: labels, the start slice, fans | 85 | | |
| | `pieces.py`: `PieceJob`, its windows (`TargetGrid` or mosaic plan) | 65 | | |
| | `basin_run.py`: run pieces in order, async-ready | 45 | | |
| | `io/mesh_index.py`: `MeshIndex`, piece writer, seam records, conformity | 85 | | |
| | `cli.py`: `--block-nodes`, pieces output, fields | 45 | | |
| | **23c total** | **405** | **563** | **648** |
| **23d** | **In parallel, resumable, stitched** | | | |
| | `basin_run.py`: `--jobs`, thread split, largest first, the memory cap | 45 | | |
| | `basin_run.py`: job hashes, skip finished pieces | 35 | | |
| | `stitch.py`: streaming `.vtk` and `.ply`, seam vertices numbered once | 130 | | |
| | `cli.py`: `rasputin stitch`, `--no-stitch` | 25 | | |
| | **23d total** | **235** | **327** | **376** |

Modules to watch: `pieces.py` at 150 (split `fans` out past 220) and
`stitch.py` at 130. 23b assumes the edge strip has written
`constraint_check_points`; if 23b came first it would carry about 40 more.
Natural cuts (23e) are not estimated: they get their own design with the
basin's inputs.

**Documentation in the same PRs** (not counted): `project_structure.md` (the
`fetch/` package and the boundary rule that network code stays in it, 23a-2;
pieces and the index, 23c); `ROADMAP.md`'s rows at each merge, and item 2's
order once Ola rules B12; `15-dem-mosaic.md` and `15c-geographic-dem.md`
pointing here where they hand basin scale on (this branch).

**Acceptance class.** 23a-1 changes how Norway's tiles are decoded, so its
`@perf` run checks the mesh hash and process time (below), though it touches
no refine code. 23b, 23c and 23d touch or drive refine, so
`docs/increments/README.md` "Acceptance" applies in full.

## Tests @tester can write red

**Invariant-critical suites, for mutation testing:** 23b's frozen-edge suite
(FE2-FE5) and the seam pass's oracle (SP1, SP2); 23c's conformity, union and
equality tests (DC2-DC4, DC6). The rest is ordinary.

### 23a-1 (Python)

- **W1, decoding a window** equals decoding the whole file with tifffile and
  slicing, `meta` included, for windows inside one block, across four, on
  block boundaries and at the raster's edge; from a local tiled TIFF and from
  the same blocks in a cache, identical arrays.
- **W2, determinism**: identical arrays for 1 and 8 decode threads.
- **W3, presence**: a block file of the wrong length is absent; a `.part` file
  is ignored; the manifest round-trips.
- **W4, the refusal**: a missing block is refused before any block is
  decoded (a `BlockSource` double that fails on `block()`), and the message
  names the fetch command with the domain and cache given.
- **W5, offline**: no `tin_engine` module on `rasputin mesh`'s path imports
  `tin_engine.fetch` or `urllib.request` in its own source (an AST check over
  `src_python/tin_engine`; pyproj imports `urllib.request`, so `sys.modules`
  cannot be the oracle), and `cli.py` imports `tin_engine.fetch` only inside
  the fetch command; and a mesh from a prepared cache
  with `socket.socket` replaced by one that raises succeeds.
- **W6, overlaps**: `BlockWindows` over two overlapping objects equals 15a's
  `assemble` on the same objects.

### 23a-2 (Python)

- **F1, the fixture**: a local range server (`http.server` in a thread)
  serving a synthetic tiled COG written by tifffile (Deflate, 16² tiles, one
  overview), counting requests, and able to fail after k responses, answer
  500, answer 200 without honouring `Range`, or change `Last-Modified`.
- **F2, the plan**: the blocks equal a brute-force test of every block's
  footprint against the region (shapely), for a domain with a thin arm
  through a block corner.
- **F3, the bytes**: each cached block equals the file's byte range; no
  request asks for a present block; coalesced requests number no more than
  the runs of wanted blocks.
- **F4, resume**: the server fails mid-run; a rerun requests only what is
  missing, and the cache then equals a clean run's, byte for byte.
- **F5, incremental**: a second, overlapping domain requests only its new
  blocks.
- **F6, refusals**: 200 instead of 206; a short body; a changed length or date;
  4xx not retried; 5xx retried three times, then refused.
- **F7, the header**: a header longer than the first prefix makes the reader
  double the prefix; the parse never reads past what was fetched; a prefix
  that cuts the full page's offset array is detected (tifffile alone would
  return the page with no offsets).
- **F8, GLO-30**: a tile missing from the tile list is "no tile", not missing.

### 23b (C++ Catch2 and through the binding)

- **FE1, K1**: every existing refine suite, golden digests included, passes
  unchanged with the mask 0. (This is the check that 23b changes nothing for
  Norway; no new test.)
- **FE2, the scan and the split**: a node exactly on a frozen grid-line edge
  with a large error is never inserted; every other node ends within
  tolerance (a brute-force oracle over all nodes not on frozen edges); with
  the mask off the same node is inserted, so the test can fail.
- **FE3, feet**: a node within ε of a frozen edge goes in itself, not a foot.
- **FE4, the quality pass** never splits a frozen edge; `skipped_frozen`
  counts it.
- **FE5, `refine_points`**: a point exactly on a frozen edge is not inserted;
  `on_frozen` and its error are reported.
- **FE6**: untested, as `EdgeProperties::bit`'s precondition is
  (`tests/cpp/unit/test_edge_properties.cpp`: an assert aborts the process
  and this tree has no death-test harness; Release CI defines `NDEBUG`).
  K2 is pinned by FE2-FE5. A harness is not budgeted in 23b: it would be the
  tree's first, for one assertion that guards a path no caller takes.
- **SP1, a grid-line seam**: the check points are the nodes on it; after
  `refine_seam`, the error is within tolerance at every node, and at 1,000
  points along the line against a bilinear surface computed in NumPy
  independently (the every-point property).
- **SP2, a general seam**: an edge with rational endpoints; the oracle
  recomputes crossings, midpoints and bilinear heights in exact rationals
  (Python `fractions`) and checks every check point within tolerance
  + 1e-9 of the output's piecewise-linear heights; tolerance 0 inserts every
  check point with a nonzero error; ties go to the smallest parameter.
- **SP3, sameness**: the edge given reversed, and the strip window shifted by
  whole blocks, give the same output bit for bit.
  In 23c, on the projected-mosaic path, a case where two pieces' windows
  select different tile sets around one seam: both pieces get the same seam
  points (read from the files, as DC2 does).
- **SP4, the bound**: insertions never exceed check points; NoData stencils
  are skipped and counted.

### 23c (Python, end to end on synthetic rasters)

- **DC1, K1**: a domain within the block size writes the same bytes as
  today's path.
- **DC2, conformity, read from the files**: for every seam edge, both pieces'
  vertex sequences are equal bit for bit. The oracle reads the piece files,
  not the seam records or the index.
- **DC3, a valid union**: in the stitched mesh every edge not on the outline
  is in exactly two triangles, the triangles' area sums to the domain's
  (shapely) to 1e-9 relative, and no two triangles overlap (brute force on a
  small case).
- **DC4, constrained Delaunay with the seams**: an exact incircle oracle
  (Python `fractions`, independent of refine's predicates) finds no violation
  when the seams count as constraints, on a case with cocircular lattice
  points; and finds some when the seams are left out, so the test can fail
  and shows what the union is.
- **DC5, the tolerance over the union**: every DEM node in the domain within
  tolerance of the stitched mesh (brute-force barycentric location, 1e-9
  slack), seams and the edge strip included; the same mesh with z shifted by
  twice the tolerance fails.
- **DC6, the main session's claim, as an equality**: on a node-only case (a
  rectangle on nodes, no features, grid-line seams, no quality start, no
  feet), the decomposed mesh equals, as a set of world-coordinate triangles,
  one `refine` run over the whole start triangulation with the same
  seam-pass points fanned in and the same seams frozen. Why equality can hold: pieces are slices of that triangulation in
  its order, and every operation in a piece is the same as in the whole run,
  restricted to the piece; with only node vertices every coordinate is an
  integer in both frames, so no rounding differs. Control: the whole run
  without freezing differs (`@tester` picks a case with a seam node over
  tolerance).
- **DC7, determinism**: `--jobs` 1 and 3, `--threads` 1 and 4, pieces in
  reversed order: identical piece files and stitched file.
- **DC8, locality**: one piece re-run alone from its job spec writes the same
  bytes.
- **DC9, geometry**: two pieces in one block; a block line through a lake
  ring; a road along a block line (one edge, both bits); a block corner inside
  the domain (four pieces at one node, which the noder must return exactly);
  a sliver piece with no node.
- **DC10, the reprojected path**: a synthetic geographic tile, `--out-crs`,
  cut into four pieces; 15c's independent final check (G6) finds 0 source
  nodes over tolerance on the union; `on_frozen` is 0.

### 23d (Python)

- **PJ1, the stitch**: the stitched file equals the union of the pieces
  (as sets of triangles, seam vertices once), numbered by first occurrence in
  piece-id order; it reads pieces one at a time (a reader double that refuses
  a second open piece).
- **PJ2, resume**: a run stopped after k pieces and rerun computes only the
  rest, and the outputs are identical to a clean run's.
- **PJ3, invalidation**: a changed tolerance changes every piece's hash.

## @perf acceptance

Per `docs/increments/README.md` "Acceptance": `tools/bench.py`'s 1 m benchmark
and thread-scaling sweep, `pmset -g batt` with every run, against the
previous increment's run in the same power state (on `NO BASELINE`, the
previous merge commit with `--tree`, back to back), evidence under
`docs/benchmarks/<date>/`.

- **23a-1:** the benchmark's mesh hash unchanged and process time within
  noise (decoding is now windowed). On ANADEM: decode throughput per block and
  per window, threads 1 to 8.
- **23a-2:** fetch the basin's 3,061 ANADEM blocks: wall time, bytes,
  requests, and a rerun that fetches nothing; an interrupted run resumed.
- **23b:** the README rule in full; the mesh hash unchanged (mask 0), refine
  within noise at every thread count.
- **23c:** the README rule, with `--block-nodes 0` on the benchmark so it
  stays comparable (hash unchanged), plus the benchmark cut by the default
  block size, recorded as a new baseline. The Velhas piece on ANADEM from the
  cache, at 1, 2, 5, 10, 20 and 50 m, cut (11 blocks at `B = 2048`) and uncut
  (15c-2's run): triangles and their difference, worst angle, maximum degree,
  0 constrained-Delaunay violations with seams as constraints, 0 nodes and 0
  source nodes over tolerance by the independent check (its control still
  failing), the seam pass's insertions, time and peak RSS.
- **23d:** scaling on the Velhas piece at 1 m over `--jobs` × `--threads`;
  then **the basin** at 50, 20, 10, 5, 2 and 1 m: triangles against the
  200-box estimate (237 M at 1 m, 215-261 M), wall time per phase, fetch
  bytes, peak RSS, and the stitched file's size. The claim to test: **peak
  RSS under 16 GiB at every tolerance**, i.e. memory no longer decides the
  tolerance. The independent final check over the union (on sampled pieces
  if a full check is too long; `@perf` says which and why).

## Questions for Ola

Numbered B1-B12, so they cannot be confused with 15's and 15c's Q1-Q17.
The literature pass with web search (2026-10-01, "Prior art") changed no
recommendation; it added a note to B4 and an option (c) to B5.

**B1. Which cuts first?**
Ola's direction (2026-10-01) was to cut on constraint edges; (a) does cut on
constraint edges, but ones the plan adds, not ones the input has, and the
pieces are not hydrological units.
- **(a) Global grid lines first (23c); BHO's Pfafstetter units later (23e),
  with the basin's own inputs. Recommended.** Grid lines work for any input,
  Norway included, cost nothing in quality (no foot can trigger next to
  them), give an exact seam pass and exact seam heights, and need no new
  data. Their price is lines in the mesh that mean nothing hydrologically
  (B4), and the mesh near a seam differs from an uncut run's (seam vertices
  placed by the 1D pass; to be measured in 23c's acceptance, cut against
  uncut).
- (b) BHO units first. The pieces are hydrological units, which is what a
  basin model wants as output; but it needs BHO for the whole basin, brings
  a vertex every ~107 m along each seam (about +37 % vertices at 50 m on the
  piece at level 7, arithmetic), and loses feet along seams. Its size is not
  estimated (23e gets its own design).
- (c) Both in 23c: over the ceiling; two PRs anyway.

**B2. When is a run cut?**
- **(a) Only when the needed window is wider or taller than `--block-nodes`;
  then along every block line through the domain. Recommended.** Small meshes
  stay bit-identical to today's (K1). It is one rule, not two code paths: an
  uncut run is the one-piece case.
- (b) Always, along the global grid: one rule with no threshold, but any
  mesh crossing a block line gets a seam, small ones included, and every
  existing mesh that does changes.

**B3. The block size.**
- **(a) 2048 nodes by default (61 km at 30 m, 20 km at 10 m), as
  `--block-nodes`, recorded in the file; 0 means never cut. Recommended.**
  About 1 GiB per piece at 1 m in the steepest terrain measured, so eight at
  once fit a 32 GB machine. The 1 m benchmark tile (25.5 M nodes) would be
  cut by default, so `tools/bench.py` passes `--block-nodes 0` from 23c on.
- (b) 4096: a quarter of the seams, up to about 4 GiB per piece at 1 m, so
  fewer at once.
- (c) 1024: four times the seams, a quarter of the memory.

**B4. The seams in the output.**
- **(a) Kept as constraint edges with a `seam` bit. Recommended.** The mesh is
  constrained Delaunay with the seams as constraints; a consumer can tell a
  seam from a river by its bit. No extra code.
- (b) Removed after the run: for each seam, the strips on both sides
  re-legalised across it and rescanned. The union is then Delaunay across
  artificial seams, but both neighbours' strips must be live at once, and it
  needs its own design (roughly 300 lines).
- *From the literature pass:* (a) is the guarantee Chernikov and
  Chrisochoides' PCDM gives (constrained Delaunay with respect to the
  separators). (b) has a published form, Bentley's stitching patent (remove
  the triangles whose circumcircle crosses a tile boundary, retriangulate
  with the neighbours' points). Linardakis and Chrisochoides get a Delaunay
  union with no post-pass by spacing the separator vertices so that no
  diametral circle is ever entered, which needs a size bound that a
  tolerance-driven refine does not have. Recommendation unchanged.

**B5. How two pieces agree on a seam.**
- **(a) The seam is frozen after a one-dimensional seam pass, which both
  neighbours compute identically. Recommended.** No exchange and no rounds;
  each piece job is self-contained. On grid-line seams no feet are lost (the
  seam vertex count against an uncut run is measured in 23c's acceptance);
  on natural seams it gives up feet next to the seam.
- (b) Exchange on shared segments (PCDM, the protocol Ola described): splits
  requested by either side applied to both, in rounds until none is new.
  Keeps feet on natural seams; needs rounds, a merge step, resumable refine
  runs and both neighbours available. Worth it only if natural seams show
  needles.
- (c) One neighbour owns the seam and the other inherits it, the scheme of
  Campos et al. 2020 (found in the literature pass): the seam is computed
  once, but the mesh then depends on which piece runs first, and their
  documentation says the results are not deterministic with more than one
  thread. Not recommended: it breaks K5. The same pass found (a)'s principle
  (a shared boundary reduced from its own data, so both sides agree without
  talking) in a Microsoft patent and in Bertilsson's 2015 thesis, which
  supports (a).

**B6. The output of a cut run.**
- **(a) Pieces and an index always; one stitched file for `--out x.vtk` by
  default, skipped with `--no-stitch`. Recommended.** `--out` keeps meaning
  one file; pieces are the durable result and what resumes a run; at the
  basin's 1 m the stitched file is several GB and can be skipped.
- (b) Pieces and the index only; `rasputin stitch` on request.

**B7. Where the cache lives.**
- **(a) `--cache DIR` or the `RASPUTIN_CACHE` environment variable, refused
  if neither is set. Recommended.** Explicit, as `--out-crs` is (Q11); one
  more option.
- (b) A default under `../rasputin_data/cache`, next to the repository:
  convenient here, surprising for anyone else.
- (c) The platform's cache directory (`~/Library/Caches/rasputin` on macOS):
  conventional, but 2 GB of DEM blocks hidden there.

**B8. Sources in the first fetch PR.**
- **(a) ANADEM (OpenTopography's COG) and GLO-30 (AWS). Recommended.** GLO-30
  is what the basin-piece measurement used and is global; about 25 of 23a's
  lines.
- (b) ANADEM only.

**B9. Resumable runs** (skip pieces whose job hash matches).
- **(a) In 23d. Recommended.** About 35 lines; a basin run interrupted at 1 m
  resumes where it stopped.
- (b) Not now.

**B10. Natural seams' vertices** (for 23e, asked now so it is not a
surprise). BHO boundaries bring a vertex every ~107 m, used as given under
Ola's input model.
- **(a) Decide when 23e is designed, on measured counts. Recommended.**
- (b) Simplify each shared chain once, between corners, with increment 22's
  area-preserving collapse, so both sides use the same chain.
- (c) Always as given.

**B11. The basin tolerance** (Q15, open). Memory no longer limits it.
- **(a) The basin run measures 50 down to 1 m and Ola chooses from the
  numbers. Recommended.** Nothing in this design depends on the answer.
- (b) Name it now, and the acceptance stops there.

**B12. The order of work** in "Order of work and PR split".
- **(a) As proposed: 15c-1, 23a-1, 23a-2, 15c-2, the edge strip, 23b, 23c,
  23d, the basin run, 23e. Recommended.** The fetch step lands early, so
  every later step runs on ANADEM. This reverses part of Q14's placement
  ruling (the edge strip right after 15c and before 15d): under (a) 23a-1,
  which replaces 15d, and 23a-2 land between 15c-1 and 15c-2, and the edge
  strip after both.
- (b) ROADMAP's current order (15c-1, 15c-2, the edge strip, then 23), with
  15c-2's acceptance on a one-off ANADEM cut. (b) keeps Q14's placement as
  ruled.

**Decided here, which Ola may overrule:** grid lines on the computation
lattice as artificial cuts; piece ids `(j, i, k)`; the seam pass computed by
both neighbours rather than once; seam heights from the seam record; the
ETag ignored; coalescing at 64 KiB gaps and 8 MiB requests; 8 connections and
three retries; the manifest holding identity, not inventory; largest piece
first; one process with threads; the index as JSON; pieces written beside the
stitched file as `<out>.pieces/`.

## Review

### Round 1, `6518336..c6c9241`: CHANGES REQUESTED (`@reviewer`)

Production LOC 0 (docs only; the two probes are measurement scripts nothing
imports). Probes re-run and reproduced; code citations, the union argument,
the arithmetic, determinism, the I/O boundary and LOC held; web spot-checks
matched. Nine blocking edits: a close precedent missed in the novelty
statement (Zygmunt and Róg 2026); a patent misquoted; the ANADEM prefix
evidence assumed tifffile raises on a short read; B12 reversed part of Ola's
Q14 placement without saying so; B1 did not name its departure from Ola's
direction and quoted an unsupported LOC figure; W5's import oracle could not
pass; FE6 had no death-test harness; the fan-validity argument used the
segment where it needed the line; DC6 compared different starts.

### Round 2, `c6c9241..6c3fbc5`: CHANGES REQUESTED (`@reviewer`)

All nine applied. One blocking item: the probe's strict wrapper can never
record a read past the prefix, because tifffile bounds-checks offsets first;
the offsets count is the check.

### Round 3, `6c3fbc5..aa1f7c2`: APPROVED (`@reviewer`)

Whole branch `6518336..aa1f7c2`, 12 commits, +1,731 / −13, production LOC 0.
Probe at 8 MiB reproduces every figure; at 1 MiB it fails on the offsets
assertion. ruff, ruff format, mypy and the governance gates green. CI not yet
run: no PR.
