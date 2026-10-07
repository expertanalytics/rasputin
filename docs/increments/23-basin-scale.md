# Increment 23: basin scale — pieces cut on constraints, windowed DEM reads, a tile cache

Status: **designed by `@architect`, 2026-10-01; B1-B14 ruled by Ola on
2026-10-01 and the design reworked to the rulings; B15 and B16 ruled
2026-10-02. 23a-1 merged (#136); 23a-2 merged (#138); 23b merged (#162,
merge commit `4f56551`; red
`45045dc`, its open points settled under "Settled after 23b's red step
(45045dc)"; green `3c464ec`, N18 in at `185081c`), 249 net production
lines after the master merge (248 before it; the merged sum takes one
more line), `@perf`'s acceptance recorded at `91c7cb5` and again at
`4541e38` (ACCEPTED, `docs/benchmarks/2026-10-04/23b-merged-acceptance.md`);
23c split in two:
23c-1 (about 215 estimated, 156 built) merged (#164, merge commit
`e56a2bf`), no `@perf` run needed, and 23c-2 (about 265) in progress on
`worktree-23c`, not pushed, its state in that branch's copy of this file;
23d onwards not implemented.** Written before `@tester`
per `docs/increments/README.md` step 1. The rulings are under "Ruled by Ola,
2026-10-01", below; the questions are kept as asked at the end, each marked
with its ruling. 23a-1 and 23a-2 are designed in full under "Windowed
source reads (23a-1)" and "The fetch step and the tile cache".

## Why this record, and why its name

Ola asked on 2026-10-01 for a basin-scale design in an increment file of its
own (`15c-geographic-dem.md`, Q15). It gets the next free number, 23 (22 is
auto-catchment). "Basin scale" because what it delivers is the ability to
mesh the São Francisco basin (635,194.5 km², 734 M ANADEM nodes) at any
tolerance on one machine, with memory set by the size of a piece, not of the
basin. It is not a sub-increment of 15: 15 is about DEM input, while this
changes the refinement core (frozen edges, a seam pass, seam removal), the run (many
pieces, in parallel) and the output (a piece per file and an index). It
absorbs 15d's window decoding and supersedes 15's Q9 ("32 GB is fine", a
dense canvas) at basin scale.

Sub-increments, each its own PR: **23a-1** (windowed decoding and the tile
cache, reading), **23a-2** (the fetch step), **23b** (frozen edges and the
seam pass, in C++), **23c** (the partition, run piece by piece), **23d**
(pieces in parallel, resumable runs, `rasputin stitch`), **23f** (seam
removal and thinning, in C++), **23g** (seam cleanup at stitching, in
Python). **23e** (sub-catchments from the DEM) is later and gets its own
design. Order and sizes are under "Order of work".

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

## Ruled by Ola, 2026-10-01

On B1-B14 (the questions as asked are kept at the end). Where the ruling
changes the design, the section named does the rework.

- **B1: (a), artificial cuts; BHO is dropped entirely as a source of
  geometry.** Ola: "I'm not interested in archaic maps". Sub-catchments come
  later, derived from the DEM by extending increment 22 to many outlets, and
  enter the mesh as a code per triangle or as ordinary constraints, never as
  seams (23e, reworded under "Sub-catchments"). BHO is used at most for
  official codes or for validation. The BHO measurements under "What was
  measured" stay as the record B1 was ruled on; nothing in the design uses
  them. Pointer, not a ruling on that record: the raster-to-vector research
  answers its design question 2 with "Rivers as constraints come from the
  BHO drainage lines" (`docs/research/raster-to-vector.md@4b085a5:523-524`, on
  branch `worktree-raster-vector`); that answer is to be revisited, and
  DEM-derived drainage is the likely replacement.
- **B2 and B3: replaced by a partition of the domain's bounding box** into
  Nx × Ny pieces, Nx · Ny about a target count, `--pieces`. Ola: "Compute
  dx, dy, based on concurrency requirements and Nx and Ny, so that Nx*Ny
  approx M*Np". The count is a parameter with a fixed default, **not** the
  detected core count, so that the mesh does not depend on the machine
  (Ola's L1 determinism ruling, increment 21). Piece sides are whole DEM
  spacings, so seams lie on lattice lines (exact seam heights, no feet
  lost). B14 replaced the size cap and the minimum (below). What is lost
  against the global grid of B3 is under "Choosing the cuts".
- **B4: (b), seams removed after the run.** The strips on both sides are
  re-legalised across each seam and rescanned, so the stitched mesh is
  Delaunay across artificial seams; this includes **local seam thinning**:
  seam vertices the tolerance does not need are removed (remove, retriangulate
  the hole, recheck the tolerance in the hole only), per seam, in parallel
  across seams, in a fixed order along a seam, with the corners where four
  pieces meet kept or handled in a final small pass. It happens at
  stitching, and only a band from each neighbour is live. Bentley's
  stitching patent is prior art. Designed under "Seam removal and thinning".
- **B5: (a), the seam frozen after the one-dimensional pass,** computed
  identically by both neighbours, with no communication; B4's cleanup
  repairs the over-density and poor triangles along seams afterwards. The
  acceptance of the stitch PR (23g) compares cut against uncut on the test
  piece after cleanup: vertices near seams, worst angles, triangle count.
  **(b), the exchange, is the named fallback** if the gap is large.
- **B6: (a).** Pieces and the index always; one stitched file by default
  for `--out`, clean (B4's cleanup); `--no-stitch` skips it. Piece files keep
  their seams.
- **B7: a new environment variable, `RASPUTIN_DATA`**, the data root. None
  exists today (the legacy code read `RASPUTIN_DATA_DIR`). The cache lives at
  `$RASPUTIN_DATA/cache`; `--cache` overrides it; a run that needs the cache
  is refused if neither is set.
- **B8: (a)**, ANADEM and GLO-30 in the first fetch PR.
- **B9: (a)**, resumable runs.
- **B10: moot**, BHO being dropped.
- **B11: (a).** `@perf` measures the decomposed basin run from 50 m down to
  1 m, then Ola chooses the tolerance.
- **B12: (a)**, the proposed order, knowingly reversing part of Q14's
  placement (Ola chose it). The order now carries 23f and 23g ("Order of
  work").
- **B13: (c).** The Velhas piece's and the basin's BHO outlines stay as
  measurement domains now; 23e derives both from the DEM, with one comparison
  run of the two outlines; once the DEM-derived ones work, the BHO ones are
  abandoned. Ola: "Go back, we choose option c: a, then b. When b works, we
  can abandon a."
- **B14: one parameter, `--memory-budget`, replaces `N_MIN`, `N_MAX` and
  the per-piece memory cap.** Default 16 GB, fixed, not read from the
  machine. A domain is cut only when its estimated memory at the requested
  tolerance exceeds the budget; pieces are sized to fit it; no hard upper
  limit anywhere; `--pieces` stays, as a request for more pieces. The same
  parameters give the same mesh on any machine, and a machine too small runs
  out of memory. Ola: "we should not limit huge discetisations based on less
  performant hardware. If you have to choose a large number, you should be
  allowed to do so. On limited resources, you should not expect to reproduce
  huge catchments." On the old minimum: "1M DEM points seams way to small"
  (small domains have no runtime or memory problem); on the default: "I'd
  say 16 GB is good." The rule is under "Choosing the cuts".

Also ruled on 2026-10-01, for `15c-geographic-dem.md`: Q17's BHO-outline
fixture is replaced. Ola: "yes, switch to a DEM-derived test catchment". The
test domain is a catchment derived by increment 22 from the ANADEM extract;
no BHO, and no licence question (amended there).

## What was measured for this design

By `@architect`, 2026-10-01, with the scripts in `docs/increments/23-probes/`
(run as `python docs/increments/23-probes/<script> ../rasputin_data/sao_francisco_piece`;
each was run for this record and prints the figures below).

**The partition, on the Velhas piece and the basin** (`partition.py`, run
after the B14 ruling). The rule under "Choosing the cuts", at the default
`--memory-budget` (16 GiB) and no `--pieces` request, on a 30 m lattice in
EPSG:31983. Both outlines are BHO's, kept as measurement domains until 23e
(B13 (c)).

| domain | window | domain covers | tolerance | cells | meet the domain | cell nodes | largest / mean area | seams inside |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Velhas piece | 4,208 × 7,347 = 30.9 M | 42 % | 0.5 m and up | 1 × 1 | 1 | 30.9 M | 1.00 | 0 |
| Velhas piece | | | 0 | 1 × 2 | 2 | 15.5 M | 1.02 | 73 km |
| basin | 41,332 × 50,297 = 2,079 M | 34 % | 50 m | 2 × 2 | 3 | 520 M | 1.24 | 769 km |
| basin | | | 10 m | 2 × 3 | 5 | 346 M | 1.76 | 1,344 km |
| basin | | | 5 m | 3 × 3 | 7 | 231 M | 1.63 | 2,438 km |
| basin | | | 1 m | 5 × 7 | 23 | 59.4 M | 1.93 | 5,110 km |
| basin | | | 0 | 8 × 10 | 45 | 26.0 M | 1.65 | 7,970 km |

The Velhas piece is one piece at every tolerance but 0, so today's mesh;
cut on request, at 1 m `--pieces 16` gives 3 × 5 cells of 2.06 M nodes (12
pieces) and `--pieces 64` 6 × 11 cells of 0.47 M (44 pieces). For the basin
at 1 m, `--pieces 16` changes nothing and `--pieces 64` gives 7 × 9 cells of
33.0 M (37 pieces). The basin's window here is 2.08 G nodes; the 2.31 G of
the basin-piece README (Surprise 1) is that README's canvas, not reconciled
with this one.

**BHO is an exact coverage** (`bho_coverage.py`; the record B1 was ruled on,
not used by the design since BHO was dropped). The 1,163 BHO 2017 5k
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
  bound that would keep it away. **The guarantee we drop, and win back:**
  their union is Delaunay; the union of our pieces is constrained Delaunay
  with the seams as constraints, which is PCDM's guarantee (below). Under
  B4 (b), ruled, the stitched mesh gets Delaunay across the seams back by a
  post-pass (Bentley's patent, below; "Seam removal and thinning"); piece
  files keep the weaker guarantee.
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
  constraint is needed"). It is B5 (b), the named fallback, and its union
  guarantee (constrained Delaunay with respect to the separators) is the one
  our piece files have.
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
    lattice resolution; not taken, because a partition line would carry a
    vertex at every node along it into both pieces at every tolerance, and
    the cleanup's thinning (B4 (b)) would then have to remove nearly all of
    them; the seam pass places only what the seam's own error needs.
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
    for B4 (b), seams removed after the run, which Ola ruled. **What we take:**
    the band (triangles near the boundary, from both tiles) as the only part
    retriangulated. **What differs:** we keep the band's triangles and
    re-legalise them by flips behind a fence of frozen edges, and grow the
    band where a fence edge fails the incircle test, rather than deleting and
    retriangulating; we rescan the band against the tolerance; and we thin
    the seam's vertices, which the patent (as summarised) does not.
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
- **Seam thinning: vertex removal under a tolerance.** Devillers, "On
  deletion in Delaunay triangulations", *IJCGA* 12(3):193-205, 2002,
  doi:10.1142/S0218195902000815 (Crossref): deleting a vertex changes the
  Delaunay triangulation only inside the vertex's star, so retriangulating
  the star's polygon Delaunay-wise restores the whole. Schroeder, Zarge and
  Lorensen, "Decimation of triangle meshes", *ACM SIGGRAPH Computer
  Graphics* 26(2):65-70, 1992, doi:10.1145/142920.134010 (Crossref; content
  recalled): remove a vertex, retriangulate the hole, keep the removal if an
  error criterion holds. Lee, "Comparison of existing methods for building
  triangular irregular network models of terrain from grid digital elevation
  models", *IJGIS* 5(3):267-285, 1991, doi:10.1080/02693799108927855
  (Crossref; content recalled): the *drop heuristic*, removing grid points
  from a full TIN while the vertical error stays within a tolerance.
  Mostafavi, Gold and Dakowicz, "Delete and insert operations in
  Voronoi/Delaunay methods and applications", *Computers & Geosciences*
  29(4):523-530, 2003, doi:10.1016/S0098-3004(03)00017-7 (Crossref). **What
  we take:** the drop heuristic, restricted to the seam pass's vertices and
  rechecked in the hole only (Devillers' locality is what makes "in the hole
  only" exact). **What differs:** the recheck is the exact sup-norm scan at
  grid nodes (and source check points on the reprojected path), not an
  error at the removed point; candidates are taken in a fixed order, not by
  least error, so the result does not depend on a priority queue's ties.
- **Pfafstetter coding.** Verdin and Verdin, "A topological system for
  delineation and codification of the Earth's river basins", *J. Hydrology*
  218(1-2):1-12, 1999, doi:10.1016/s0022-1694(99)00011-6 (Crossref). Each
  level divides a basin into nine units (four tributary basins, five
  interbasins), coded by one more digit; BHO's ottocodes are this. Kept as
  the record B1 was ruled on: BHO is no source of geometry (B1), and its
  codes are at most attached to DEM-derived units later (23e).
- **Cloud-optimised GeoTIFF.** OGC 21-026, "OGC Cloud Optimized GeoTIFF
  Standard", version 1.0, approved 2023-05-08, published 2023-07-14,
  doi:10.62973/21-026 (read, `https://docs.ogc.org/is/21-026/21-026.html`;
  the year was recalled in the first pass). Headers and offset arrays at the
  start, tiles addressable by range request. The probe above confirms
  ANADEM's copy behaves so.
- **HTTP range requests.** RFC 9110 (2022): `Range` §14.2, `Content-Range`
  §14.4, 206 §15.3.7, `If-Range` §13.1.5, which bars a date validator to a
  client holding an entity tag and makes a date strong only under §8.8.2.2
  (read). 23a-2 uses single ranges and checks `Content-Range` and
  `Last-Modified` on each response instead (its "Decided here" 2).
- **GDAL's `/vsicurl/`** is the method 23a-2 follows without the dependency
  (`CLAUDE.md` §2): block reads by range request, consecutive ranges merged
  (`GDAL_HTTP_MERGE_CONSECUTIVE_RANGES`), parallel single-range requests
  because multi-range GETs are "not supported by a majority of servers
  (including AWS S3 or Google GCS)" (`GDAL_HTTP_MULTIRANGE`), retries
  (`GDAL_HTTP_MAX_RETRY`), all in `https://gdal.org/en/stable/user/configoptions.html`
  (read 2026-10-02). It differs in persistence: GDAL's cache is in memory
  per process, rasputin's is a directory of blocks that outlives the run.

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
grid nodes, and an output independent of order and thread count. The seam
removal of B4 (b) adds nothing new either: a band re-legalised across a tile
border is Bentley's patent, and dropping vertices while the vertical error
holds is Lee's drop heuristic, local by Devillers' deletion result. No claim
is made. Before any is made public: the full text of Zygmunt and Róg, a Google
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
After the rulings, for B4 (b): Crossref lookups of Devillers 2002,
Schroeder et al. 1992, Lee 1991 and Mostafavi et al. 2003, and one web
search, "merge parallel Delaunay subdomains remove interface constraint edges
vertex removal decimation terrain tolerance", which found Bentley's patent
again, a 3D divide-and-conquer merge phase and decremental Delaunay patents,
and nothing that thins a removed seam under a vertical tolerance.

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
  edge vocabulary (`features.py`). Every cut is *artificial*: a line the plan
  adds (below, a partition line on the lattice). Natural cuts along input
  constraints (BHO units) were dropped with B1. The core never sees the name: refine is handed `frozen_mask =
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
- **Seam unit**: a maximal chain of seam edges between two corners. It has
  one piece on each side; the cleanup works per unit ("Seam removal and
  thinning").

### The run, end to end

```
cli.mesh --dem anadem-v1 --domain D --out-crs C --tolerance T [--pieces P] [--memory-budget B]
  |
  v  [global, vectors only: no DEM is read here]
plan_basin(request) -> BasinPlan                     (Python, pure; frozen data)
  |  domain, features in the target CRS (15b, 16b)
  |  cuts: the partition of the window into Nx x Ny, lines on the lattice   (23c)
  |  build_pslg -> node -> triangulate       the existing engine, run once        (cli._engine)
  |  labels = landcover.regions(triangles, seam edges)   one label per piece      (16c)
  |  per piece: PieceJob(id, start slice, its seam edges, target window, source window)
  v
run_plan(plan, jobs=J, threads=T)  [per piece, in parallel, local]   (async; asyncio.to_thread)
  |  source window: plan_mosaic + assemble, tiles decoded by window  -> resample (15c D3)   (23a-1, 15c)
  |  seam pass for each of its seam edges: refine_seam(strip, edge, T)        (23b, C++)
  |  start = slice with seam edges split into fans at the seam-pass points    (23c, NumPy)
  |  refine(window, start, frozen_mask=seam)       phase 1                     (14-21, 23b)
  |  refine_points(edge-strip points, frozen_mask)  the edge strip (Q14)       (edge strip)
  |  refine_points(source nodes, frozen_mask)       phase 2, reprojected path  (15c)
  |  write the piece file + its seam record, then free everything
  v
index (conformity check on every seam)                                        (23c)
  v  [stitching, unless --no-stitch; per seam unit in parallel, bands only]
cleanup(unit): two bands -> re-legalise across the seam -> rescan -> thin     (23f, 23g)
final pass (serial): corners, and what a unit could not finish in its band   (23g)
stitched file: pieces' cores + cleaned bands, vertices numbered once         (23d, 23g)
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
   mesh along a natural cut still differs from a whole-domain run. Natural
   cuts were dropped with B1. An artificial cut adds edges; they stay in the
   piece files, and the cleanup removes them from the stitched file (B4).
6. *"The cut and the split exchange must be deterministic functions of the
   input."* **Agreed, and one more input must be excluded: the machine.** A
   cut chosen from physical memory or the thread count would make the mesh
   depend on the computer. The cut is a function of the domain, the lattice,
   the tolerance, `--pieces` and `--memory-budget` (with the `b(T)` table),
   all recorded in the file (K5). Ola's ruling on B2 and B3 says the same: the count is a parameter
   with a fixed default, not the detected core count.
7. *"Shared-edge splits are common, given Q14's crossing/midpoint check points
   along constraints."* **Correct under Q14 as designed** (check points
   inserted after refine, by `refine_points`). Here they are moved *before*
   refine on seams, into the seam pass, which is what removes the exchange.

## Choosing the cuts

### The partition (23c; B2, B3 and B14 as ruled)

- **The lattice.** Node `(R, K)` of the computation grid: 15c's target grid
  on the reprojected path (`x = K·h`, `y = −R·h`, J6), or the mosaic's own
  lattice from 15a's reference node on a projected DEM meshed directly.
- **The window.** The domain's bounding box grown by the cell diagonal (15b)
  and snapped outward to the lattice: `cols × rows` nodes, `N = cols · rows`,
  first node `(R0, K0)`.
- **The rule** (Ola: "Nx*Ny approx M*Np"; B14), integer arithmetic on the
  window, with `P = --pieces` (default 1, no request), `B = --memory-budget`
  in bytes and `b = b(T)` the estimated bytes per window node at the
  tolerance `T` (the table below):
  1. The count: `P' = max(P, ⌈N · b / B⌉)`. If `P' ≤ 1`: **one piece**, no
     cut, today's mesh bit for bit (K1).
  2. Near-square cells: `Nx = max(1, round(√(P' · cols / rows)))`,
     `Ny = max(1, round(P' / Nx))`. The cell count approximates `P'` (Ola's
     "approx") and is not a floor: rounding may give fewer cells than asked
     (`--pieces 16` on the Velhas window gives 3 × 5). Only the budget bound
     of step 3 is guaranteed.
  3. Whole spacings: `dx = ⌈cols / Nx⌉`, `dy = ⌈rows / Ny⌉` nodes. While
     `dx · dy · b > B`, add one to `Nx` if `dx ≥ dy`, else to `Ny`, and
     recompute. Then `Nx = ⌈cols / dx⌉`, `Ny = ⌈rows / dy⌉`, so no column or
     row of cells is empty; the last ones may be narrower.
  4. The partition lines are the lattice columns `K0 + i·dx` (`0 < i < Nx`)
     and rows `R0 + j·dy` (`0 < j < Ny`), each from the window's edge to its
     edge.
  `partition.py` (23-probes) implements this rule as written here, and its
  figures are under "What was measured".
- **The chains.** Each line is clipped to the domain as linework (16b R6) and
  enters as breakline chains with mask `seam`. The noder nodes them like any
  input. A cell corner is where two lines cross, at a node whose coordinates
  are whole metres when `h` is, so on the noder's 1 mm grid; that the
  noder's crossing arithmetic returns it exactly for two axis-parallel
  segments is for 23c's red suite to pin (DC9), not assumed here.
- **Why lattice lines**, as Ola ruled:
  1. *No quality cost.* A node not on a lattice-line seam is at least one
     cell from it, and 20b's ε is at most half a cell (`foot_epsilon`'s cap),
     so a foot on such a seam never triggers: freezing it changes nothing a
     foot would have done.
  2. *An exact seam pass.* The nodes on a lattice line lie exactly on the
     seam, and the bilinear surface along it is piecewise linear between
     them, so after the seam pass the tolerance holds at **every point** of
     the seam against the DEM's bilinear surface, not only at check points.
  3. *Exact heights.* A seam vertex on a lattice line is nearly always a
     node (a cell-side midpoint only at a rounding tie, N14), so its z is
     the node's value; both pieces write the record's z in any case (step 5).
- **Several pieces in one cell.** Where the domain enters a cell twice, the
  labelling gives two components and so two pieces.
- **Piece ids.** `(j, i, k)`: cell row, cell column, and the component's
  rank by its lowest start-triangle index. A function of the input alone.

### The memory estimate and the defaults (B14 as ruled)

- **`--memory-budget`, default 16 GB, read as 16 GiB (2^34 bytes)**, a
  constant: never the machine's memory, so the mesh does not depend on the
  machine (K5). No upper limit on it, on `--pieces` or on a piece: a run whose
  pieces do not fit the machine runs out of memory. 15a R7's refusal at
  half of physical memory goes too (B15, ruled (a)).
- **`b(T)`, bytes per window node at tolerance `T`**, from the basin-piece
  sweep (`docs/benchmarks/2026-10-01/basin-piece/README.md`): its max-RSS fit,
  0.55 GiB + 310 B per triangle, gives about 17 B per grid node (the intercept
  over the 30.9 M-node grid, leaving out the ~66 MiB interpreter) plus 310 B
  per triangle; the
  triangles per node are the Velhas piece's (its triangles at `T` over its
  12,957,257 domain nodes), above the basin's p90 at 1 m (894 against 757
  per km²), so the estimate is high for most terrain. At `T = 0` every node
  is a vertex, two triangles per node. In whole bytes, rounded up:

  | T | 0 | 1 m | 2 m | 5 m | 10 m | 20 m | 50 m |
  |---|---:|---:|---:|---:|---:|---:|---:|
  | triangles per node | 2 | 0.805 | 0.443 | 0.154 | 0.062 | 0.024 | 0.006 |
  | `b(T)`, bytes | 637 | 267 | 155 | 65 | 37 | 25 | 19 |

  Between two columns `b` is linear in `T`, rounded up (the counts fall
  convexly, so the chord is above them); above 50 m it is 19. Every window
  node is costed as a domain node, which overestimates a domain that covers
  part of its window (the Velhas piece's 1 m run: 7.7 GiB estimated, 3.60 GiB
  measured). The table and `B` are recorded in the index.
- **What it gives** (`partition.py`, "What was measured"): the Velhas piece
  is one piece at every tolerance but 0; the basin is cut into 2 × 2 cells
  at 50 m and 5 × 7 at 1 m (23 pieces). The 1 m benchmark tile (25.5 M nodes,
  6.8 GB) and Bygdin at 10 m (3.05 M nodes) stay one piece, so their meshes
  and `tools/bench.py`'s hash are unchanged (K1). Parallelism across pieces
  is asked for with `--pieces`; at the default a domain under the budget
  runs as today, one piece on refine's threads.

### What is lost against the global grid of B3

B3 cut on one global grid, block lines every 2048 nodes of the lattice. The
partition follows the domain instead. Lost:

- **Seams no longer lie at the same coordinates across separate runs.** Two
  catchments meshed separately, or one domain at two values of `--pieces`,
  have their seams in different places. Their stitched files are unaffected
  by the seams themselves (removed), but the mesh near where a seam ran
  differs between the two runs; piece files of separate runs never line up.
- **A small edit moves every seam.** Growing the domain's box by one node
  can move every partition line, so the mesh changes along every seam, not
  only near the edit; and a resumed run (B9) recomputes every piece, since
  every job spec changed. Under the global grid only the touched blocks
  would have rerun.
- **Pieces are not reusable between domains.** A piece of the global grid
  was the same for every domain that covered the whole block; a partition
  cell is not.
- **Seams move with the tolerance**, since `b(T)` sets the count: the basin
  at 50 m and at 1 m is cut differently.
- Kept: seams on lattice lines (exact seam pass, exact heights, no feet
  lost), the mesh independent of the machine, one piece for small domains.

### Sub-catchments (later, 23e; B1 as ruled)

Not cuts. Sub-catchments are derived from the DEM by extending increment 22
from one outlet to many, and enter the mesh as a code per triangle (as 16c's
land-cover labels do) or, where a boundary must be in the mesh, as ordinary
constraints: never as seams, so they do not decide the partition. BHO is used
at most to attach official codes to DEM-derived units, or to validate them.
23e also derives the Velhas piece's and the basin's outlines from the DEM
and runs the one comparison of each against its BHO outline (B13 (c)); the
acceptance runs move to the DEM-derived outlines once they work. 23e gets
its own design with the basin's own inputs.

## The seam protocol

**Ruled (B5 (a)): seams are frozen, after a one-dimensional seam pass that
each neighbour computes identically.** No exchange, no rounds. What the
frozen seam costs (seam vertices the 2D refine would not have placed, and
the thin triangles beside them) is repaired at stitching by B4's cleanup
("Seam removal and thinning").

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
     counted. On a grid-line seam the crossings are the nodes, and 15f's
     generator adds the midpoints of the cell sides between them (`2c + 1`
     check points); the surface is linear between nodes, so the greedy
     inserts nodes, a midpoint only at a rounding tie (N14, "Settled after
     23b's red step").
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

**Determinism.** The plan is a function of the input, the lattice, the
tolerance, `--pieces` and `--memory-budget`; the noder and the CDT are deterministic; piece ids come from
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
finite, and it is deterministic if each pass is. **What it buys:** seam
vertices only where a piece's 2D refine asks for them, so no over-density
along seams to clean up. **What it costs:** rounds in which both neighbours
must be resumed or rerun (so either both are live at once, against locality,
or a piece is rerun per round), a merge step, and a second kind of refine run
(resume with injected splits). On lattice-line seams it buys no feet, since
none triggers there. **Ruled: the named fallback**, taken only if 23g's
comparison of cut against uncut after cleanup shows a large gap (the
thresholds are under "@perf acceptance", 23g).

## Tolerance, the final check and the constraint check points per piece

- **The guarantee, decomposed** (K3). Every valid DEM node inside the domain
  is within `--tolerance` of the mesh: a node strictly inside a piece by that
  piece's refine, as today; a node on a seam by the seam pass. On a grid-line
  seam the tolerance holds at every point of the seam against the bilinear
  surface; on any other seam, at its check points. That is the piece files'
  guarantee. **In the stitched file** the seams are gone, and every node in a
  band the cleanup changed is rescanned, and every node in a thinned hole
  rechecked, so the guarantee is the interior one everywhere: every valid
  node within tolerance; the every-point property along former seams is not
  kept (triangles now cross the line).
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
- **Memory at run time** (B14): the runner (23d) starts a piece only while
  the running pieces' estimates (`dx · dy · b(T)`) fit under
  `--memory-budget`, and always starts one when none is running. It can
  delay a piece, never change or refuse one. 15a R7's refusal at half of
  physical memory is deleted (B15, ruled (a)).

## Windowed source reads (23a-1)

Designed in full for `@tester` on 2026-10-02. **15a's path stays the path**:
`plan_mosaic` from headers, then `assemble` from loads. 23a-1 changes two
things under it: a load decodes only the blocks its placement needs, and a
catalogue source's tiles come from the cache instead of from files. 15c-2's
`SourceWindows` is unchanged (its `TileWindows` slices the assembled source
tile); at basin scale that tile is one piece's source window, so nothing
holds the basin.

### Decided here

1. **No `BlockWindows`.** A windowed `assemble` is the same thing (every
   object meeting the window, decoded by window, overlaps by R5), so a second
   assembler is not written; `SourceWindows` stays 15c-2's to define.
2. **The catalogue is data, outside `fetch/`**: `tin_engine/sources.py`. The
   mesh path must name catalogue keys and must not import `tin_engine.fetch`
   (W5), so the catalogue cannot live there. `CacheManifest` moves to
   `io/repository.py` for the same reason. `fetch/` (23a-2) imports both.
3. **`--dem` becomes `list[str]`**, turned into paths in the CLI. `Path`
   normalises `./glo30` to `glo30`, so the ruled "`./glo30` is a path"
   cannot be honoured on `Path`.
4. **`IndexWindow` moves to `io/models.py`** (re-exported by `mosaic.py`),
   because `io/` imports nothing first-party but `io.models`.
5. **A sparse block (byte count 0) in a needed window is refused.**
   tifffile fills it with 0, a valid elevation; a DEM's NoData must be
   stored, not implied.
6. **The missing-block check runs once per plan, before any decode**
   (`check(plan)`), and `decode_window` repeats it per window as a backstop.
   The CLI adds the fetch command to the message: only the CLI knows the
   flags as given.
7. **Threads:** `decode_window(..., threads=4)`, not a flag; the output does not depend on it (W2).
   Increment 30c replaces the fixed 4 by `os.cpu_count()` workers (`30c-dem-read-speed.md`, 3.1).
8. **Local reads under a lock** (`seek` and `read` on one stream). A read
   is short against a decode; `os.pread` would exclude `BytesIO` fixtures.
9. **`elevation_source` names the source id and the catalogue's `credit`.**
   The manifest does not repeat the credit (one place per fact).
10. **Geographic sources still refuse** until 15c-2 (`read_header`'s GeoKey
    2048 refusal), after the cache checks. 23a-1's offline mesh tests use a
    projected entry added to `SOURCES` by `monkeypatch`.

### Types and functions

`io/cog.py`. Imports `io.models` and `io.geotiff`; never opens a file.

```
class BlockSource(Protocol):
    page: tifffile.TiffPage                      the full-resolution page, from the header
    where: str                                   for messages: a file name, or "<source>/<object>"
    def block(self, index: int) -> bytes: ...    the block's stored bytes
    def missing(self, indices: Sequence[int]) -> tuple[int, ...]: ...
class CacheError(ValueError)                     manifest absent, wrong source, header hash mismatch
class NotCached(CacheError)                      missing: int, needed: int, where: str
class LocalTiffBlocks                            (page, stream, name); missing() is always ()
def blocks_meeting(page, window: IndexWindow) -> tuple[int, ...]       pure, ascending
def window_meta(meta: RasterMeta, window: IndexWindow) -> RasterMeta   x_min + col0*dx, y_max - row0*dy
def decode_window(source: BlockSource, meta: RasterMeta, dtype, window, *, threads=4) -> DemTile
```

`decode_window`:

- A window outside the raster, or with fewer than one row or column, is a
  `ValueError` (a caller's bug, not a file's).
- `source.missing(blocks_meeting(...))` first; any → `NotCached`, before any
  `block()` call.
- Per block, on a `ThreadPoolExecutor(threads)`: `block(i)` →
  `page.decode(data, i)` → the segment cropped to the raster (edge tiles
  come padded, the last strip short) and to the window → cast by
  `PROMOTION` into its own slice of one preallocated output. Slices are
  disjoint, so thread count cannot change a value. A decode failure is a
  `GeoTiffError` naming `where` and the block index (`geotiff._stage`).
- Tiled and stripped pages alike (`page.chunks`). DTM10 is 512² LZW tiles
  with three overviews; ANADEM is a COG.
- The output becomes the tile through the public `DemTile(...)`
  constructor, one extra copy per window. A departure, accepted in review:
  `_adopt` was designed here, but 15a's suite (`test_mosaic.py`, M15)
  reserves `_adopt` for `mosaic.py`. Widen it only if `@perf` shows the copy
  matters.

Checked with tifffile 2026.9.20 before writing this: `page.decode(bytes, i)`
decodes a block from its bytes alone; it does so on a page parsed from a
prefix that ends before the first block; and it still works after its
`TiffFile` is closed.

`io/geotiff.py`: `read_page(source, *, nodata) -> (RasterMeta, dtype,
TiffPage)`, `read_header` plus the page, through the same `_header`, so
every refusal is the same.

`mosaic.py`: `assemble(plan, load, needed=None, *, load_window=None)`. With
`load_window`, each placement is `load_window(name, placement.source)`, its
meta must equal `window_meta(placement.meta, placement.source)` (else the
existing "changed since it was listed" refusal), and it is copied whole.
Without it, nothing changes (`catchment.py` and 15a's suite).

`io/repository.py`:

- `TiffDemRepository.load_window(name, window)`: opens the file read-only,
  `read_page`, `LocalTiffBlocks`, `decode_window`. `check(plan)` is a no-op.
- `CacheManifest` and `CachedObject` (frozen; 23a-2 writes them, 23a-1
  reads them), fields under "The cache" below.
- `CachedBlocks(directory, page, where)`: `block(i)` reads
  `blocks/<i // blocks_across>/<i % blocks_across>.bin`; `missing` lists the
  indices whose file is absent or not exactly the header's byte count. A
  `.part` file has another name and is never read.
- `CacheRepository(cache, source)`. Construction reads no file. The same
  `footprints`, `load`, `load_window` and `check` as `TiffDemRepository`.
  `footprints()` reads `manifest.json` (absent: `NotCached` with
  `needed = 0`, "not in the cache"; another source id: `CacheError`), then
  each listed object's `header.bin`, refusing one whose sha256 differs from
  the manifest's (`CacheError`, "re-fetch with --refresh"). The footprint's
  name is the object id. `check(plan)` sums `missing` and `needed` over
  every placement's `blocks_meeting` and raises one `NotCached` for the
  plan.

`dem_input.py`: `DemRequest` gains `cached: CachedSource | None`
(`source: str`, `cache: Path`), exactly one of `sources` and `cached`;
`repository_for` returns a `CacheRepository` for it, labelled by the source
id; `open_dem` calls `repository.check(plan)` and then
`assemble(..., load_window=repository.load_window)`.

`tin_engine/sources.py`: `RemoteSource` and `SOURCES` (fields under "Types"
below); it imports Pydantic only.

`cli.py` (`mesh`): `--cache DIR`; `RASPUTIN_DATA` read once, in the command;
`cache_root(option, environ) -> Path | None`. `--dem` with exactly one
value that is a key of `SOURCES` is that source; one value is a path
otherwise, and a key mixed with paths is refused. A key with no cache root
is refused (B7's message, under "The CLI"). A `NotCached` becomes the
usage error with `; run: rasputin fetch <key>` plus the run's `--domain` or
`--bbox`, `--out-crs` (once 15c-2 adds it) and `--cache` as given. `dem_tiles` is written for a
catalogue source as for a directory.

### Unchanged by 23a-1

- **No dense canvas at basin scale.** Each piece holds a target window (its
  start slice's bounding box grown by the cell diagonal, snapped outward to
  the lattice) and a source window (15c D2's `source_region` of that target
  window: its image in the source CRS, grown by two source cells). Both are
  about `B²` nodes plus margins.
- **Header before pixels** (15c J10): coverage from the headers, missing
  blocks from `check(plan)`, both before any block is decoded.
- **A projected DEM meshed directly** (Norway): a piece's window is
  `plan_mosaic(footprints, piece bounds)` (pure, 15 R12), each tile decoded by
  window.

## The fetch step and the tile cache

Designed in full for `@tester` on 2026-10-02 (23a-2). `rasputin fetch`
copies what a mesh of the same domain reads into the layout 23a-1 reads;
nothing in it is on `rasputin mesh`'s path (K7).

### Decided here (23a-2)

1. **Fetch the box, not the outline.** 23a-1's `check(plan)` and
   `decode_window` need every block of a placement's window, and every
   window is a box (15c-2's source box, 23c's piece windows inside it). So
   fetch takes the blocks meeting the domain's box, in the run's frame,
   grown and moved into the source CRS ("Planning"). For the basin on
   ANADEM that is about 8,300 blocks against the 3,061 meeting the grown
   outline, and for the Velhas piece 160 against 82 (both counted on the
   BHO outlines grown by three cells, with the probe's tie point and step,
   in ANADEM's CRS; the `--out-crs` box's image is a little larger).
   `@perf`'s piece fetch also took 160 blocks. Making the read side outline-aware (an absent block
   outside the needed region read as NoData) would save the difference once,
   at the cost of a second meaning of "missing"; not done.
2. **No `HEAD`, no `If-Range`.** Every response carries `Content-Range`'s
   total and `Last-Modified`; each is checked against the manifest, so
   identity is checked on every request at no cost. `If-Range` with a date
   is barred to a client that has an entity tag (RFC 9110 §13.1.5), and
   ANADEM's ETag is a placeholder (measured), so it is not used either.
3. **One range per request.** No multi-range GETs: S3 does not serve them
   (GDAL's `GDAL_HTTP_MULTIRANGE` documentation says so of AWS S3 and GCS).
4. **A request is all or nothing.** Its blocks are written only after the
   whole body arrived and every length matched; a failed request writes
   nothing, so resuming has request granularity (at most 8 MiB).
5. **One fetch per source at a time**: `fcntl.flock` on
   `<cache>/<source>/.lock`, non-blocking, refused if held. The OS drops it
   when the process dies, so there is no stale lock; `.part` files are
   cleared under it. Linux and macOS only, as CI is.
6. **The header is parsed with geographic CRSs allowed.** `read_page` gains
   `geographic: bool = False`; only `fetch/` passes `True`. The mesh path's
   refusal (23a-1 W9) stays until 15c-2.
7. **A sparse block (byte count 0) is written as an empty file**, with no
   request. It is then present, and the mesh's own sparse refusal (23a-1,
   decided 5) names it, rather than a misleading `NotCached`.
8. **Source notes: `<cache>/<source>/NOTICE.txt`**, rewritten on every run
   (not `--dry-run`) from the catalogue by `sources.notice(source)`: the
   credit, the licence note (Art. 6(c) for GLO-30) and the works the
   distributor asks to cite (`RemoteSource.cite`; for ANADEM, Laipelt et
   al. 2024, as OpenTopography's acknowledgement asks). The catalogue stays
   the one place; the file is a rendering. What the mesh file carries is
   B16.
9. **The requests log** appends `(region_sha256, date)` once per distinct
   pair: the sha256 of the region as given (the domain's WKB and CRS, or the
   box's four numbers and CRS), so equal requests on one day log once. The
   manifest's field is renamed from `domain_sha256` (23a-1 only reads it).
10. **No fsync per block** (thousands of files); the manifest is written to a
    temporary file, fsynced and `os.replace`d.
11. **`--refresh`** discards the source's cache (manifest, headers, blocks,
    under the lock) and fetches this request anew. Blocks fetched for other
    domains go with it; they are of the old remote copy anyway.

**Departures by `@developer`, accepted in review (23a-2, round 1):**

- `RangeClient` accepts a reply cut short at the end of the file
  (`Content-Range` ending at `total - 1`) only when the range asked for
  ran past it (`stop > total`); any other short reply is refused.
- GLO-30's tiles are chosen from one box grown as if the spacing were a
  fixed 0.001° (`fetch/run.py`'s `TILE_SPACING`, more than any tile's below
  80°), not per tile; each tile's blocks are then planned from its own
  header.
- The box is moved by `crs.transform_bounds`, the one site that wraps
  pyproj's, not by pyproj directly in `fetch/`.
- `FetchReport` also carries the per-object plans (`plans`), and `fetch`
  takes a `progress(done, total)` callback in bytes; the CLI's stderr line
  is that callback.
- A one-file source's object id is the URL's file name without its
  extension (for ANADEM, `anadem_v1_compressed_COG`).
- `rasputin fetch` has `--domain-crs`, with `mesh`'s meaning (the
  `--domain` file's CRS).
- `cite` is written, in `NOTICE.txt` and in the mesh file, only when the
  source has citations.
- B16 (a) is implemented for both formats: `licence_note` and `cite` as
  `.vtk` fields and as `.ply` header comments.

### Types (frozen Pydantic)

```
tin_engine/sources.py (23a-1; data, imports Pydantic only)
RemoteSource      id ("anadem-v1", "glo30"); kind ("one-cog" | "cog-tiles");
                  url, or url_template plus tile_list_url; crs (expected, checked
                  against each header); nodata; credit; licence_note;
                  cite: tuple[str, ...] = () (23a-2)
SOURCES           the catalogue: a Mapping[str, RemoteSource] of data, two entries
notice(source) -> str                                                     (23a-2)
io/repository.py (23a-1 reads, 23a-2 writes)
CacheManifest     source id, crs, rasputin version, objects: {object id: CachedObject},
                  requests: ({region_sha256, date}, ...)
CachedObject      url, content_length, last_modified, header_sha256, header_bytes,
                  block shape (rows, cols)
CacheWriter       (root, source); a context manager holding the lock (decided 5):
                  put_header, put_block, put_manifest, put_notice, discard
fetch/plan.py (23a-2, pure)
FetchRequest      source id; domain (DomainPolygon) or box (Bounds, in out_crs, else
                  the source CRS); out_crs or None; margin (default 4); connections
                  (default 8); dry_run; refresh
ObjectPlan        object id, url, block indices, ranges: ((start, stop), ...), bytes
FetchPlan         source id; objects: (ObjectPlan, ...); no_tile: (name, ...)
fetch/run.py (23a-2)
FetchReport       objects, blocks needed, present, fetched, empty; bytes; requests;
                  no_tile; seconds
```

The date is `datetime.date` (standard library; `CLAUDE.md` §2). The objects
a mesh reads are the manifest's, sorted by id; the blocks it has are the
directory's.

`anadem-v1` is OpenTopography's COG
(`https://opentopography.s3.sdsc.edu/raster/ANADEM/ANADEM_be/anadem_v1_compressed_COG.tif`,
DOI 10.5069/G9736P4G, as ruled in Q17). `glo30` is the AWS bucket, one COG per
1° tile, with `tileList.txt` saying which tiles exist; a tile absent from the
list is sea, reported as "no tile", not as missing.

### Planning (`fetch/plan.py`, pure)

- **The header** of each object: a prefix of 1 MiB, doubled until the
  full-resolution page has an offset and a byte count for every block of
  its grid (`block_grid`), refused past 64 MiB. tifffile does not raise on
  a short prefix, it logs and returns a page without offsets (measured at
  64 KiB and 1 MiB), so the count is the check. The parse reads from a
  stream that raises on any read past the prefix. ANADEM needs 8 MiB
  (measured). The prefix is `header.bin`.
- **The frame** is `out_crs` if given, else the source CRS: the CRS the
  mesh's box is in (15c-2's target grid, or 15a's DEM CRS).
- **The source box**: the domain's bounds in the frame (or the box as
  given), grown by `margin` times the source's north-south spacing in
  metres, moved into the source CRS with pyproj's `transform_bounds`
  (densified, `always_xy`), then grown by two source cells as 15c-2's
  `source_region` grows its box. 15c-2 grows the target box by `√2·h`:
  42 m at the CLI's default `h` on ANADEM (30 m) and 44 m on GLO-30
  (31 m), against the 119 m and 124 m that `margin` 4 gives. When the
  frame is the source CRS, the box is grown by `margin + 2` cells. A box crossing
  ±180° or reaching a pole is refused (15c-2's refusal).
- **The blocks**: the source box as an `IndexWindow` of the page, clipped
  to the raster, then 23a-1's `blocks_meeting`. A box that meets no block
  is refused. For GLO-30 the objects are the listed tiles whose 1° square
  (from the name) meets the source box, each with its own window.
- **The ranges**: the missing blocks sorted by offset, runs whose byte gaps
  are at most 64 KiB coalesced into one `(start, stop)` of at most 8 MiB
  (a single larger block is its own range).
- **K7 by construction**: the mesh's windows lie inside the source box for
  any partition, so a mesh after a fetch of the same domain and options
  finds every block. `check(plan)` stays the authority: a Python-API mesh
  whose spacing `margin` does not cover is refused naming the fetch, and
  `FetchRequest.margin` is the remedy.

### The cache (`io/repository.py`, the one module in `io/` that opens files, 15 Q3)

```
<cache>/<source-id>/manifest.json                       identity, written atomically
<cache>/<source-id>/NOTICE.txt                          credit, licence, citations
<cache>/<source-id>/.lock                               held by a running fetch
<cache>/<source-id>/<object-id>/header.bin              the parsed prefix
<cache>/<source-id>/<object-id>/blocks/<row>/<col>.bin  one block, its exact bytes;
                                                         <row>, <col> in blocks, a strip's col is 0
```

- **The directory is the inventory; the manifest is the identity.** A block is
  present exactly when its file exists with its byte count from the header.
  Blocks are written as `<col>.bin.part` and renamed, so a crash never
  leaves a half block, and nothing has to be kept in step with the files.
- **Identity.** A known object's prefix is re-read each run
  (`header_bytes` long) and its sha256, total length and `Last-Modified`
  compared with the manifest; every block response is checked the same way
  (decided 2). Any change refuses: "the remote copy changed; `--refresh`
  re-fetches it". Meshing never checks the remote; it trusts the manifest,
  because it is offline by rule.
- **Order of writes**: headers and the manifest first (a new object is
  listed before any of its blocks), then blocks, then the request appended
  to the manifest. A crash at any point leaves a cache 23a-1 reads.

### Downloading (`fetch/http.py`, `fetch/run.py`)

- `http.py`, the one importer of `urllib.request`: `RangeClient.get(url,
  start, stop) -> RangeResponse(data, total, last_modified)` and
  `get_text(url)` (GLO-30's tile list, fetched each run, not cached). A
  `Range` response must be 206 with `Content-Range` `bytes start-(stop-1)/total`
  and exactly `stop - start` bytes; a 200 is refused without reading its
  body. Timeouts, connection errors, short bodies and 5xx are retried three
  times with delays `(1, 2, 4)` s (a field, zero in tests); 4xx never. A
  refusal is `FetchError`, naming the URL and the range.
- `run.py`: `async def fetch(request, source, client, writer) -> FetchReport`.
  Each range is `asyncio.to_thread(client.get, ...)` under an
  `asyncio.Semaphore(connections)`; its blocks are cut by the header's
  offsets, length-checked and put (in the same thread). `client` and
  `writer` are parameters, so the run is tested against the local server
  and a writer in `tmp_path`, and an API worker awaits it in its own loop.
- **Resumable:** present blocks are never requested; `.part` files are
  removed on start. **Incremental:** a second domain fetches only its
  missing blocks. **`--dry-run`** reads the headers and prints the plan
  (objects, no-tile names, blocks present and missing, bytes, requests),
  and writes nothing, not even the directory.

### The CLI

```
export RASPUTIN_DATA=../rasputin_data
rasputin fetch anadem-v1 --domain basin.geojson --out-crs EPSG:31983 [--cache DIR] [--dry-run] [--refresh] [--connections 8]
rasputin mesh --dem anadem-v1 --domain basin.geojson --out-crs EPSG:31983 --tolerance 5 --out basin.vtk
```

**The cache root (B7, ruled):** `--cache DIR` if given, else
`$RASPUTIN_DATA/cache`, else refused: "no cache: set RASPUTIN_DATA (the data
root; the cache is $RASPUTIN_DATA/cache) or pass --cache". `RASPUTIN_DATA` is
new; nothing in `src_python` or `tools` reads an environment variable of that
name today, and the legacy code's was `RASPUTIN_DATA_DIR`
(`legacy/rasputin/__init__.py:4`). The test helper
`tests/python/gpkg_fixtures.py` has a module constant `RASPUTIN_DATA`, the
path `../rasputin_data`, not the variable; 23a-1 leaves it alone (a test
reading the variable is `@tester`'s call). The variable is read once, in
`cli.py`, into the request model; nothing below the CLI reads the
environment.

**`--dem`:** a catalogue key (`anadem-v1`, `glo30`) names a cached source;
anything else is a path, as today; a file whose name is a catalogue key is
written as a path (`./glo30`), which is why the option is read as text
(23a-1, "Decided here" 3). Only a catalogue source needs the cache, so a
mesh from local files needs neither variable nor option. A block meshing
needs and the cache lacks is refused before any decode, naming the command:
"anadem-v1: 37 of the 8,300 blocks this domain needs are not in DIR; run:
rasputin fetch anadem-v1 --domain basin.geojson --out-crs EPSG:31983".

**`rasputin fetch SOURCE`** (23a-2): `SOURCE` a catalogue key; `--domain`
or `--bbox`, exactly one, read as `mesh` reads them (the same parsers), so
the command the mesh's refusal prints runs as given; `--out-crs`, the frame
("Planning"); `--cache` and `RASPUTIN_DATA` through 23a-1's `cache_root`;
`--connections` (default 8). It imports `tin_engine.fetch` inside the
command and runs `asyncio.run(fetch(...))`. Output: the report, one line
per object, and a progress line to stderr at each tenth of the bytes. A
refusal (`FetchError`, `CacheError`, a held lock, a changed remote) exits
1 with its message; a usage error exits 2.

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
replaced by one that raises (K7). The mesh file's `elevation_source` names
the source id and the catalogue's `credit` (23a-1).

## Memory and parallelism

**Per piece, within the budget** (B14): a piece's estimate,
`dx · dy · b(T)`, is at most `--memory-budget`, and the runner keeps the
running pieces' sum under it ("Tolerance, the final check and the
constraint check points per piece"), so memory
follows the budget, not the basin. The estimate is high for most terrain
(the Velhas piece at 1 m: 7.7 GiB estimated, 3.60 GiB measured). The cleanup
at stitching holds two bands per seam unit, far less than a piece.

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
the Velhas piece cut with `--pieces 64` (44 pieces, the largest 1.59 times
the mean area, so about 4 % of the area), if the largest piece holds 4-8 % of the work
(steepness varies; not measured), its single-thread refine is about 0.5-1.0 s,
and the total, 12.95 s of one-thread work over 10 cores, about 1.3-1.7 s
with the balance bound (pieces started largest first take at most about
`(W / J) · (1 + r · J / P)` for work `W` on `J` cores, `r` the largest
piece's work over the mean), against 6.32 s for the whole piece at 10 threads
today: about 4-5×, before the cleanup's cost. At 1 m the budget alone cuts
the basin into 35 cells, 23 of which meet it; `--pieces` asks for more where
the cores outnumber them. The limits then become cores, memory bandwidth, the global
vector step and the cleanup, not the serial phase. With one piece (any
domain under the budget, the benchmark included) nothing changes, and 21d stays
the route for one piece's serial phase. All of this is arithmetic for
`@perf` to measure.

**The basin at 1 m, as arithmetic:** 237 M triangles (215-261 M, the box
sample) at the piece's 12.95 s per 10.43 M triangles on one thread is about
5 min of refine on one core, about 40 s on 8 if it scales; resampling 734 M
source nodes at the prototype's 30.9 M per 0.99 s on 8 threads is about 24 s;
plus phase 2, decoding between the 3,061 blocks meeting the outline and the
about 8,300 of its box (windows are boxes), and writing.

## Output

- **Pieces and an index, always, when a run is cut.** Each piece is written
  as soon as it is done, as a normal `.vtk` or `.ply` (whichever `--out`
  names), to `<out>.pieces/<id>.<ext>`, with a seam record beside it: for
  each of its seam edges (by its index in the start triangulation), the
  vertex sequence with `(x, y, z)`. Piece files keep their seams (B6, ruled).
  From 23g on, the piece also writes a band record per seam unit
  (`<id>.band-<unit>.npz`, "Seam removal and thinning"). The index,
  `<out>.pieces/index.json` (`MeshIndex`, frozen Pydantic,
  `io/mesh_index.py`), holds the CRS, the tolerance, the partition
  (`pieces`, `memory_budget`, `b(T)`, the lattice, `Nx`, `Ny`, `dx`, `dy`), the
  source's identity and credit, and per piece its file, sha256, counts and
  window.
- **Conformity is checked when the index is written** (K4): for every seam
  edge, the two pieces' records must be equal bit for bit; a difference fails
  the run naming the seam and the two pieces.
- **One stitched file by default, clean** (B6 (a), ruled): `--out basin.vtk`
  still means one file, with the seams removed and thinned (B4 (b), the next
  section). The stitcher streams piece by piece (memory: one piece and the
  cleaned bands), writing each piece's core (its triangles in no band) and
  the cleaned bands, numbering vertices by first occurrence in piece-id
  order, so a vertex on a former seam gets one number and a removed one none.
  It reads the inputs twice, counting and then writing, because the legacy
  `.vtk` and the PLY both need counts in their headers. `--no-stitch` skips
  it; `rasputin stitch <out>.pieces` does it later, cleanup included. Until
  23g lands, 23d's stitched file keeps the seams.
- **How a consumer reads it.** The stitched file reads as today, with no
  `seam` bit left in it. A consumer of pieces reads `index.json` and opens
  the files it wants; each is a complete mesh of its piece. Seam edges carry
  the `seam` bit, named in the file's vocabulary, so a consumer can tell a cut
  from a river.
- **Resumable runs** (23d; B9, ruled): each piece's job has a hash of its spec
  (input identities, options, rasputin version), stored beside its file; a
  rerun skips pieces whose hash matches. The cleanup of a seam unit is keyed
  by its two pieces' hashes, so it reruns when either does.

## Seam removal and thinning (23f, 23g; B4 (b) as ruled)

**What it does.** At stitching, for each seam unit (a chain of seam edges
between two corners, one piece on each side): take a band of triangles from
each neighbour, clear the seam bit, re-legalise the band so it is Delaunay
across the former seam, rescan it against the tolerance, and remove the seam
vertices the tolerance does not need. Units run in parallel; a final serial
pass takes the corners and whatever a unit could not finish inside its band.
Only the two bands of a unit are live; a piece's core is never read, except
on the rare deferred path (below).

### Data, written by each piece job (23g)

- **Units** come from the plan: the global start triangulation's seam edges,
  split at corners, each with its two pieces and its vertex chain oriented as
  the seam pass orients it (`a < b` by `(x, y)`). A function of the input.
  Units are numbered by their smallest seam-edge index in the start
  triangulation, which the noder and the CDT make a function of the input;
  "lowest-numbered" and "unit id" below mean this number.
- **Zones.** Each triangle of a finished piece is given to exactly one of
  the piece's units: a triangle with an edge on a unit goes to that unit
  (to the lowest-numbered one if it has edges on two, at a corner); any
  other to the unit nearest its centroid (distance to the unit's polyline;
  ties to the lower unit id). The zones of a piece are disjoint and
  cover it, so two units never claim one triangle.
- **Bands.** The band of unit `u` on one side is the triangles of `zone(u)`
  within `R` rings of `u`'s vertices (ring 1: triangles with a vertex on `u`;
  ring k+1: triangles sharing a vertex with ring k). `R = 4`, decided here
  (B14). Its **fence** is every band edge whose other triangle is not in the
  band (an edge of the piece's boundary is a constraint already). For each
  fence edge the record keeps the outside triangle's third vertex, which is
  all the fence check needs, and whether that triangle lies in another
  unit's band (a **shared** fence edge, on a zone border, mostly near
  corners): that unit may change the triangle concurrently, so the recorded
  vertex can go stale.
- **The band record** (`<id>.band-<unit>.npz`, NumPy): vertices `(x, y, z)`
  as float64, triangles, edge masks, the fence edges and their outside
  vertices, and the band's triangle indices in the piece file. Vertices on
  `u` are equal bit for bit in both records (K4), and are matched by exact
  `(x, y)`.

### The unit cleanup (23g drives it; 23f is its C++)

For one unit, from its two band records and a DEM strip covering the bands'
box (decoded and resampled like a piece window; each node computed from its
own coordinates, J6 and D3, so the values are those the pieces used):

1. **Merge** the two bands into one mesh. `u`'s edges lose the `seam` bit;
   an edge left with no bit is free. An edge that is also a feature keeps its
   other bits and stays a constraint. Other units' seam edges in the band
   keep theirs. Fence edges get a run-local `fence` bit.
2. **Re-legalise**: Lawson flips over every free edge (the existing
   `legalise_all`), constraints and fence never flipped.
3. **Rescan**: `refine` over the band with the run's tolerance and options,
   `frozen_mask = seam | fence` (23b), so nothing is inserted on a fence or
   on another unit's seam; on the reprojected path also `refine_points` with
   the band's source check points (15c phase 2), and the edge strip's check
   points on the band's remaining constraints. Former seam edges are free
   now, so feet and splits on them are allowed. A node lying exactly on a
   fence edge is skipped, which loses nothing: a fence edge never changes,
   so the mesh there is what it was when that node was within tolerance.
4. **Thinning**, in a fixed order along `u` from its `a` end: each vertex
   the seam pass placed on `u` (every vertex of `u` but its two corners) is a
   candidate. A candidate whose star leaves the band is deferred. If its star
   lies in the band and every edge at it is free:
   - retriangulate the star's polygon on the side: ears by exact orientation,
     the first valid ear in the polygon's order, then Lawson flips inside the
     hole with exact incircle and the tree's tie rule (Devillers 2002: the
     Delaunay triangulation changes only inside the star, so the result is
     Delaunay);
   - recheck **in the hole only**: every DEM node in the new triangles'
     closed node sets (the scan refine uses), and on the reprojected path
     every source check point inside them;
   - if every one is within tolerance, commit (the star's k triangles become
     k − 2, the vertex is retired); otherwise nothing changes. Nothing has to
     be undone, because nothing is written before the check.
   Removal never moves a vertex and never touches a constraint.
5. **Fence check**, run over the band as it stands after thinning: for each
   free fence edge, the exact incircle test of its inside triangle against
   the recorded outside vertex. A strict violation (an exact tie is not one)
   means the Delaunay repair wants to cross the fence; it is recorded as
   **deferred** (below), never ignored.
6. **Write** the cleaned band (triangles, vertices, retired vertices) and the
   unit's deferred items, sorted.

**The final pass** (serial, 23g), in this order, each sorted by unit id and
position: (a0) **shared fence edges** are rechecked against the current
state, since both sides may have changed; a violation joins (a);
(a) **deferred fence violations**: a region of `R` rings around the
edge, taken from the current state (cleaned bands, and the piece's core when
needed, read one piece at a time), fenced, re-legalised, fence-checked and
grown (`R` doubled) until no violation is left, then rescanned as in step 3;
(b) **deferred thinning candidates**, as step 4 over the current state;
(c) **corners where three or more pieces meet** inside the domain: a lattice
crossing the cut added, so a removal candidate like any seam vertex, its
star taken from the cleaned bands of the units around it. Then, at every
vertex where two or more units meet (kept corners on the outline, a hole or
a feature included), any seam edge that no band held loses its seam bit, and
a region of `R` rings around it is handled as in (a): fenced, re-legalised,
fence-checked, grown until no violation is left, then rescanned as in step
3. A corner where a
seam meets the outline, a hole or a feature is an input-constraint vertex
and is **kept**: removing it would merge two constrained edges that the
noder's 1 mm snap need not leave collinear. So one artificial vertex stays
on the outline per seam that reaches it. The final pass writes patches that
the stitcher applies.

**Why it is correct.**
- *Delaunay across the seams.* After the final pass every free edge of the
  stitched mesh is locally Delaunay: edges in no band were already
  (constrained Delaunay in their piece, and no seam edge among them: every
  edge of `u` is inside band(u) or one of its fence edges, since a triangle
  with an edge on `u` is in `zone(u)`; the one exception, an edge at a
  corner whose two triangles each also have an edge on a lower-numbered
  unit, is taken by the final pass's (c), as an (a) region, rescanned), band edges by steps 2-3, fence
  edges by step 5, run after thinning (against an outside that did not
  change) or by the final pass's (a0) and (a) (where it may have), and holes
  by Devillers' result. Locally Delaunay at every free edge is a
  constrained Delaunay triangulation with respect to the remaining
  constraints (Chew 1989), i.e. the input's constraints only (K11).
- *The tolerance.* A triangle outside every band and hole is unchanged; a
  band triangle is rescanned; a hole's triangles are rechecked before the
  commit. So every valid node is within tolerance in the stitched file
  (K3, stitched).
- *Termination.* Flips end (Lawson); refine ends as before; thinning visits
  each candidate once; the final pass's region grows within a finite mesh.
- *Determinism.* A unit's output is a function of its two band records and
  the DEM; the final pass is serial in a sorted order; so the stitched file
  does not depend on `--jobs` or on which unit runs first (K5).

**Locality and memory.** A unit holds two bands of `R` rings and a strip of
DEM; units are independent, so they run in parallel with `--jobs`. The
deferred path reads one piece's core at a time and is expected rare; 23g's
acceptance counts it.

**Numerical robustness.** Every decision (flip, fence check, ear, hole
legality) uses the tree's exact predicates; the error checks use the scan's
own arithmetic, so the oracle relation is the producer's
(computational-geometry skill, §3).

## What survives of 15c, the edge strip and 15d

- **15c survives whole**, and at basin scale runs per piece. `TargetGrid` (a
  sub-rectangle of one global lattice, J6) is the piece window; `resample`,
  `check_point_blocks`, `CheckPoints` and `refine_points` are used unchanged
  except for the frozen mask (23b). D8 and the Q11-Q17 rulings stand. Its
  memory cap (15a R7) goes (B15, ruled (a)). 15c's acceptance stays on the Velhas piece,
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

- **K1. Small runs untouched.** With one piece (a window whose estimate is
  within `--memory-budget`, with no `--pieces` request) and `frozen_mask` 0,
  every mesh is bit-identical to master's: 23b, 23c and 23g change nothing
  on that path, and `tools/bench.py` (one piece at the default) stays
  comparable with every stored run.
- **K2. Frozen means frozen.** In a piece, no pass inserts a vertex on a
  frozen edge: refine's split, feet, the quality pass, `refine_points` (phase
  2 and the edge strip), and 15f's loop: strip points, `refine_strip`'s DEM
  rescan and L12's on-edge insertion (N7, N16). Pinned by FE2-FE5; `split_edge`'s assertion is a
  debug guard against a later path, untested (FE6).
- **K3. The guarantee over the union.** Every valid DEM node inside the
  domain is within `--tolerance`: inside a piece by refine, on a seam by the
  seam pass; on a grid-line seam at every point against the bilinear
  surface. On the reprojected path, every source node too, except those
  exactly on a seam, or within `r(g)` of one with a strip (N7), which are
  counted (`on_frozen`) with their error. In
  the stitched file, after the cleanup, every valid node and (reprojected)
  every source node is within tolerance, seams included, by the band rescan
  and the hole recheck; the every-point property along former seams is not
  kept.
- **K4. Conformity.** Two pieces sharing a seam edge write the same vertex
  sequence on it, `(x, y, z)` bit for bit; checked when the index is
  written, and a difference fails the run.
- **K5. Determinism.** The output is a function of the inputs and the
  options (`--tolerance`, `--pieces` and `--memory-budget` included; `b(T)`
  and `R` are constants, recorded in the index), never of `--jobs`, `--threads`, the order pieces
  or seam units run in, the core count or the machine's memory.
- **K6. Locality.** A piece's files depend only on its job spec: its start
  slice, its windows and its seam edges. Re-running one piece alone gives the
  same bytes.
- **K7. Meshing is offline.** No network on `rasputin mesh`'s path; blocks
  meshing needs that the cache lacks are refused before any decode, naming
  the fetch command; for the same domain and options, mesh needs a subset of
  what fetch got (the box rule under "Planning").
- **K8. The cache is honest.** A block is present exactly when its file has
  its byte count; a fetch refuses a remote whose length, date or header
  changed.
- **K9. One global step, vectors only.** Noding, the start triangulation and
  the labelling run once for the whole domain and never read the DEM; their
  memory follows the vector input.
- **K10. The I/O boundary.** The core sees windows, meshes, edges, a
  tolerance and a mask, all as numbers in metres: no CRS, path or URL.
- **K11. The stitched file has no seams.** No edge carries the `seam` bit,
  every free edge is locally Delaunay (so the mesh is constrained Delaunay
  with respect to the input's constraints alone), and the vertices on former
  seams are only those the tolerance needs under the fixed thinning order,
  plus one per seam end on the outline, a hole or a feature.
- **K12. Cleanup locality.** A seam unit's cleanup depends only on its two
  band records, the DEM and the options; re-running one unit alone gives the
  same bytes. Only the final pass sees more than one unit.

## Degeneracy policy

- **Nodes on a seam** (every node on a grid-line seam): the seam pass's
  check points, never a piece's scan candidates.
- **A seam along the outline** (a partition line on a `--bbox` edge, say):
  not a cut, since it crosses no interior.
- **A seam along a feature edge** (a road on a partition line): one edge with
  both bits; the feature keeps its bits, and the seam's rules apply. The
  cleanup clears only the seam bit, so the edge stays a constraint and its
  vertices are not thinned.
- **A piece with no DEM node** (a sliver between a partition line and the
  outline): refine has nothing to scan; the seam pass still runs on its seam
  edges, so its boundary still meets the tolerance there.
- **A seam edge shorter than a cell** (no crossing) has one check point,
  its midpoint, and none when the midpoint's cell touches NoData (N15); with
  none, its two vertices are its only heights, as any short constraint
  edge's are.
- **NoData on a seam:** a check point whose stencil touches NoData is skipped
  and counted; a seam vertex without a height is invalid, and the triangles
  around it are carved by 14's rule in each piece.
- **A fan that is not counter-clockwise** (a vertex within rounding of a seam
  edge's line, which the noder's guarantee 14(b) rules out): refused by
  refine's `NotCounterClockwise`, with the piece and seam named.
- **Several components of the domain in one cell:** several pieces.
- **Four pieces at a cell corner inside the domain:** the corner is a node
  and a vertex of all four; the cleanup's final pass may remove it.
- **A lattice-collinear star** (a seam vertex whose star polygon has
  collinear vertices, common on a lattice): an ear with zero area is not
  valid; the first valid ear in the polygon's order is taken. Cocircular
  holes are settled by the tree's incircle tie rule.
- **A seam unit with no vertex between its corners:** nothing to thin; its
  one or more edges are re-legalised like any other.
- **A band that touches the other side's outline** (a sliver piece): the
  outline edges are constraints, so the band stops there with no fence.
- **A deferred fence violation:** never left in place; the final pass grows
  the region until it is gone (K11).
- **A source node exactly on a seam** (phase 2): skipped and counted in both
  pieces, `on_frozen`, with its error.
- **The remote changed** between fetches: refused; `--refresh` re-fetches.
- **A GLO-30 tile not in `tileList.txt`:** sea, "no tile"; a domain whose
  needed region falls there gets NoData there, as 15a's coverage rule says
  (refused if 15a would refuse it).

## Not in scope

- **Sub-catchments from the DEM** (23e): increment 22 extended to many
  outlets, as codes per triangle or ordinary constraints; built later with
  the basin's own inputs (ROADMAP item 2.5). Natural cuts are dropped (B1).
- **The exchange protocol** (B5 (b)): the named fallback, built only if
  23g's acceptance shows a large gap.
- **Noding in pieces.** The global vector step is fine at basin scale with
  the outline and the cuts; DEM-derived rivers or a polygonised MapBiomas
  might not be, and would need their own look.
- **A process pool** for pieces, unless the GIL shows in the profile.
- **Fetching other datasets** (MapBiomas): the catalogue is data, so
  adding a COG source is an entry, but vector sources need their own fetch.
- **Machines other than one:** pieces are pure data in and out, which keeps a
  cluster possible; nothing here builds it.

## Order of work and PR split

### The order (B12, ruled (a))

ROADMAP's "Order of work from 2026-09-30", item 2, read: measure a piece
(done), 15c, 15d, parallel refine or domain decomposition, the basin's own
inputs. Ola ruled B12 (a), the proposed order; it now carries the cleanup:

1. **15c-1**, as designed (`15c-geographic-dem.md`): check points and the
   final check in C++. Needs nothing new; ready for `@tester`.
2. **23a-1**, windowed decoding and the cache's read side. Replaces 15d.
3. **23a-2**, the fetch step. With 23a-1, rasputin fetches and caches its own
   ANADEM, so 15c-2's acceptance no longer needs a one-off cut while
   `metadados.snirh.gov.br` answers 403.
4. **15c-2**, as designed, its acceptance on ANADEM from the cache.
5. **The edge strip** (Q14), as designed, writing `constraint_check_points`
   (`docs/increments/15f-edge-strip.md`, 15f-1 and 15f-2).
6. **23b**, frozen edges and the seam pass, in C++.
7. **23c**, the partition, piece by piece, with pieces and the index as
   output.
8. **23d**, pieces in parallel, resumable runs, `rasputin stitch` (the
   stitched file keeps its seams until 23g).
9. **23f**, seam removal and thinning, in C++.
10. **23g**, the cleanup at stitching; the stitched file is clean from here.
11. **The basin run**, `@perf`, no code: 50, 20, 10, 5, 2 and 1 m (B11);
    then Ola chooses the tolerance.
12. **23e and the basin's own inputs**: sub-catchments from the DEM, the
    DEM-derived Velhas and basin outlines (B13 (c)), rivers
    (DEM-derived drainage likely; the pointer under B1), MapBiomas land
    cover.

**Against a ruling, knowingly.** This reverses part of Q14's placement ruling
(the edge strip right after 15c and before 15d, `15c-geographic-dem.md`, Q14
and "The edge strip"): 23a-1, which replaces 15d, and 23a-2 land between
15c-1 and 15c-2, and the edge strip after both. Ola chose it (B12 (a)).

**Why this order.** 23a has no dependency on the core and gives every later
step real ANADEM data, which is what Ola asked the fetch step for ("or our
performance will drop while we wait for download"). 23b needs `refine_points`
(15c-1) and the edge strip's generator. 23c needs 15c-2's `TargetGrid` and
`resample` per window, and 23b. 23d is what makes the basin fast, but
nothing in it changes a mesh, so it goes after the first correct cut runs.
23f needs 23b's frozen mask and nothing of 23c-23d, so it could move earlier;
it sits next to its driver so that one review reads both. 23g needs 23c's
pieces, 23d's stitcher and 23f. The basin run waits for 23g, so that it
measures what ships (stitched, clean). Release hardening and 21d stay
deferred, as ROADMAP has them.

### PR split and LOC

Counted in `CLAUDE.md` §2's unit. Estimates; the worst case applies 39 %
(increment 10's overrun), with 60 % (15a's `mosaic.py`) beside it. Every PR
stays under 700 at both; the largest, 23c, is 688 at +60 %.

| PR | what | est. | +39 % | +60 % |
|---|---|---:|---:|---:|
| **23a-1** | **Windowed decoding and the cache, reading** | | | |
| | `io/cog.py`: `BlockSource`, `LocalTiffBlocks`, `blocks_meeting`, `window_meta`, `decode_window`, the two errors | 75 | | |
| | `io/geotiff.py`: `read_page` | 8 | | |
| | `io/models.py`: `IndexWindow`, moved from `mosaic.py` | 0 | | |
| | `mosaic.py`: `assemble(..., load_window=)` | 12 | | |
| | `io/repository.py`: `load_window`, `CacheManifest`, `CachedObject`, `CachedBlocks`, `CacheRepository`, `check` | 85 | | |
| | `dem_input.py`: `CachedSource`, the cache repository, `check` before `assemble` | 20 | | |
| | `sources.py`: `RemoteSource`, the catalogue (two entries, B8) | 35 | | |
| | `cli.py`: `--cache`, `RASPUTIN_DATA`, `--dem` as text, the key, the fetch command in the refusal | 35 | | |
| | **23a-1 total** | **270** | **375** | **432** |
| **23a-2** | **The fetch step** | | | |
| | `fetch/plan.py`: header completeness, frame and source box, blocks, GLO-30 tiles, coalescing, identity comparison | 95 | | |
| | `fetch/http.py`: `RangeClient`, 206, `Content-Range` and length checks, retries, `FetchError` | 70 | | |
| | `fetch/run.py`: lock, headers, identity, manifest, bounded async download, report, `--dry-run` | 85 | | |
| | `io/repository.py`: `CacheWriter` (lock, `.part` and rename, atomic manifest, notice, discard) | 50 | | |
| | `io/geotiff.py`: `read_page(..., geographic=)` | 8 | | |
| | `sources.py`: `cite`, `notice` | 15 | | |
| | `cli.py`: `rasputin fetch` | 50 | | |
| | **23a-2 total** | **373** | **518** | **597** |
| **23b** | **Frozen edges and the seam pass, C++** | | | |
| | `lattice_mesh.hpp`: the frozen mask, `is_frozen`, the assertion in `split_edge` | 15 | | |
| | `scan.hpp`: nodes on a frozen edge skipped (only for triangles with one) | 30 | | |
| | `refine.hpp`: `frozen_mask`, no feet on frozen edges | 15 | | |
| | `quality.hpp`: skip on frozen, `skipped_frozen` | 10 | | |
| | `refine_points.hpp`, `strip_scan.hpp`: skip on frozen (N7's radius with a strip), `on_frozen` once per point (N6), a frozen strip edge refused (N16) | 30 | | |
| | `seam.hpp`: `refine_seam` (one-dimensional greedy over `constraint_check_points`, N10's carving) | 115 | | |
| | `bindings/core.cpp`, `_core.pyi` (N13's arrays; `refine_strip`'s `frozen_mask` if 15f-3 is in, N18) | 70 | | |
| | **23b total** (after "Settled after 23b's red step") | **285** | **396** | **456** |
| **23c** | **The partition, piece by piece** | | | |
| | `decompose.py`: the partition rule, `b(T)`, the lines as chains, `BasinPlan` | 95 | | |
| | `features.py`: the `seam` property | 5 | | |
| | `pieces.py`: labels, the start slice, fans | 85 | | |
| | `pieces.py`: `PieceJob`, its windows (`TargetGrid` or mosaic plan) | 65 | | |
| | `basin_run.py`: run pieces in order, async-ready | 45 | | |
| | `io/mesh_index.py`: `MeshIndex`, piece writer, seam records, conformity | 85 | | |
| | `cli.py`: `--pieces`, `--memory-budget`, pieces output, fields | 50 | | |
| | `mosaic.py` and `catchment.py`: the refusals at half of physical memory deleted (B15 (a)) | 0 | | |
| | **23c total** | **430** | **598** | **688** |
| **23d** | **In parallel, resumable, stitched (seams kept)** | | | |
| | `basin_run.py`: `--jobs`, thread split, largest first, admission under the budget | 45 | | |
| | `basin_run.py`: job hashes, skip finished pieces | 35 | | |
| | `stitch.py`: streaming `.vtk` and `.ply`, two passes, vertices numbered once | 130 | | |
| | `cli.py`: `rasputin stitch`, `--no-stitch` | 25 | | |
| | **23d total** | **235** | **327** | **376** |
| **23f** | **Seam removal and thinning, C++** | | | |
| | `lattice_mesh.hpp`: retire a vertex and two triangle slots; compaction on output | 35 | | |
| | `mesh/vertex_removal.hpp`: the star's polygon, ears, Lawson in the hole, on the side | 110 | | |
| | `refinement/thin.hpp`: ordered candidates, star-in-band test, hole recheck, commit | 75 | | |
| | `mesh/fence.hpp`: exact incircle of each fence edge against its outside vertex | 25 | | |
| | `bindings/core.cpp`, `_core.pyi` | 55 | | |
| | **23f total** | **300** | **417** | **480** |
| **23g** | **The cleanup at stitching** | | | |
| | `seams.py`: units from the plan, zones, bands of `R` rings, fences, band records | 100 | | |
| | `seam_cleanup.py`: one unit (merge, bits, legalise, rescan on a DEM strip, fence check, thin, write) | 110 | | |
| | `seam_cleanup.py`: the final pass (shared fences, deferred items, corners, unheld seam edges), patches | 85 | | |
| | `stitch.py`: cores and cleaned bands, retired vertices dropped | 30 | | |
| | `basin_run.py`, `cli.py`: units in parallel, keyed for resume | 20 | | |
| | **23g total** | **345** | **480** | **552** |

Modules to watch: `pieces.py` at 150 (split `fans` out past 220),
`stitch.py` at 160 after 23g, and `vertex_removal.hpp` at 110, whose ear
search stays simple because the hole is a star polygon. 23b assumes the edge
strip has written `constraint_check_points`; if 23b came first it would carry
about 40 more. B4 (b)'s estimate as asked was roughly 300 lines; with
thinning, the final pass and the band records it is about 645 over 23f and
23g, which is why it is two PRs. 23e is not estimated: it gets its own
design with the basin's inputs.

**Documentation in the same PRs** (not counted): `project_structure.md` (the
`fetch/` package and the boundary rule that network code stays in it, 23a-2;
pieces and the index, 23c; the band records and the cleanup, 23g);
`ROADMAP.md`'s rows at each merge, and item 2's order as ruled (B12 (a)),
which this branch's PR can already write; `15-dem-mosaic.md` and
`15c-geographic-dem.md` pointing here where they hand basin scale on (this
branch); the `RASPUTIN_DATA` variable in the README's setup (23a-1).

**Acceptance class.** 23a-1 changes how Norway's tiles are decoded, so its
`@perf` run checks the mesh hash and process time (below), though it touches
no refine code. 23b, 23c, 23d, 23f and 23g touch or drive refine or mesh
code, so
`docs/increments/README.md` "Acceptance" applies in full.

## Tests @tester can write red

**Invariant-critical suites, for mutation testing:** 23b's frozen-edge suite
(FE2-FE5) and the seam pass's oracle (SP1, SP2); 23c's conformity, union and
equality tests (DC2-DC4, DC6); 23f's removal and hole-recheck oracles (VR2,
VR3, TH2); 23g's stitched-Delaunay and tolerance oracles (SC2, SC3, SC8). The
rest is ordinary.

### 23a-1 (Python)

**Fixtures.** A source is `geotiff_fixtures.micro_tiff` with `tile=(16, 16)`
on a 50 × 70 grid (edge tiles padded), Deflate with the floating-point
predictor, one overview page, NoData cells scattered; variants: stripped
(`rowsperstrip=8`, last strip short), LZW, int16 (promoted to float32). A
cache is written by a test helper from the same bytes: `header.bin` is the
file up to its first block, each block the file's own byte range at
`blocks/<row>/<col>.bin`, `manifest.json` from `CacheManifest`. The oracle
for pixels is tifffile's whole-page decode, sliced; nothing in it shares
code with `decode_window`.

- **W1, a window equals the whole, sliced**, `meta` equal to `window_meta`
  of the whole file's: windows inside one block, across four, exactly on
  block edges, reaching the last row and column, one row, one column; every
  variant; from `LocalTiffBlocks` and from `CachedBlocks`, identical. With
  `rasputin_data` present (skipped otherwise): a DTM10 tile, a window across
  a 512² tile corner.
- **W2, determinism**: `threads=1` and `threads=8` give byte-identical
  arrays (`tobytes`, so NaN compares).
- **W3, only the needed blocks**: a recording `BlockSource` double sees
  `block(i)` exactly for `blocks_meeting`'s indices, each once; and
  `blocks_meeting` equals a brute-force test of every block's pixel
  rectangle against the window.
- **W4, the cache's read side**: a block file of the wrong length is
  missing; a `.part` file is ignored; the manifest round-trips through JSON;
  a `header.bin` whose sha256 differs from the manifest's, and a manifest of
  another source, are `CacheError`; no manifest is `NotCached` ("not in the
  cache"). `CacheRepository` construction reads no file.
- **W5, the refusals**: a missing block is `NotCached` before any block is
  decoded (a double whose `block()` fails); `check(plan)` over two objects
  counts the plan's missing and needed blocks, raised once; the CLI's message
  names `rasputin fetch <key>` with the run's `--domain` or `--bbox`,
  `--out-crs` and `--cache` as given. A block of byte count 0 in the window
  (the header's TileByteCounts entry patched) is a `GeoTiffError` naming the
  block; outside the window it is never looked at.
- **W6, windowed assembly**: on 15a's fixtures (quadrants with overlaps,
  mixed dtypes), `assemble(plan, load, load_window=...)` equals
  `assemble(plan, load)`, array and seams; a `load_window` whose meta is not
  `window_meta(placement.meta, placement.source)` gets "changed since it was
  listed". 15a's suite stays green unedited.
- **W7, offline**: no `tin_engine` module on `rasputin mesh`'s path imports
  `tin_engine.fetch` or `urllib.request` in its own source (an AST check over
  `src_python/tin_engine`, with `importscan.py`; pyproj imports
  `urllib.request`, so `sys.modules` cannot be the oracle), and `cli.py`
  imports `tin_engine.fetch` only inside a function; a mesh from a prepared
  cache, with a projected entry put in `SOURCES` by `monkeypatch` and
  `socket.socket` replaced by one that raises, succeeds and gives the same
  vertices and triangles as the same file meshed by path.
- **W8, the cache root and the key (B7)**: `--cache DIR` wins over
  `RASPUTIN_DATA`; with only `RASPUTIN_DATA=R` the cache is `R/cache`; with
  neither, a catalogue `--dem` is refused naming both, and a path `--dem`
  runs as today (`monkeypatch` on the environment); `./glo30` is a path,
  `glo30` a key; a key and a path together are refused; nothing under
  `src_python/tin_engine` but `cli.py` reads `RASPUTIN_DATA` (a source
  scan); `elevation_source` names the key and its credit.
- **W9, geographic, until 15c-2**: a cached geographic source passes the
  cache checks and is then refused with the reader's own GeoKey message, not
  a cache message.

### 23a-2 (Python)

No network: every URL is `http://127.0.0.1:<port>/...`, put in `SOURCES` by
`monkeypatch`; retry delays are zero. Async tests use `pytest-asyncio`.

- **F1, the fixture**: a range server (`http.server` in a thread, its own
  handler) serving files written by `geotiff_fixtures.micro_tiff` (Deflate,
  16² tiles, one overview), one projected (EPSG:32633) and one geographic
  (EPSG:4326), plus a GLO-30-style tile list. It logs every request's range,
  and can fail after k responses, answer 500 or 404, answer 200 ignoring
  `Range`, send a short body, a wrong `Content-Range`, or a changed
  `Last-Modified` or length.
- **F2, the plan**: the blocks equal a brute-force test of every block's
  pixel rectangle against the source box (projected and reprojected
  frames); the source box contains every node of the domain's box in the
  frame, grown by `margin` (sampled densely, moved with pyproj); a box
  outside the raster is refused; ranges cover exactly the missing blocks,
  no range spans a gap over 64 KiB or exceeds 8 MiB unless it is one block.
- **F3, the bytes**: each cached block equals the file's byte range; the
  request log shows no present block requested and one request per range;
  `header.bin` is the file's prefix; a sparse block (TileByteCounts patched
  to 0) is an empty file and was never requested.
- **F4, resume**: the server fails after k responses; the rerun requests
  only the missing blocks, and the cache tree then equals a clean run's,
  byte for byte. A stray `.part` file is gone after any run.
- **F5, incremental**: a second, overlapping domain requests only its new
  blocks; the manifest's requests list both, and the same request twice on
  one day is listed once.
- **F6, refusals**: 200 instead of 206 with a body far larger than the
  socket buffers (several MiB), refused, and the server sent fewer bytes
  than the body; a short body; a wrong `Content-Range`; a changed
  length or `Last-Modified` on a block response, and on a known object's
  prefix, each naming `--refresh`, with no block of that response written;
  4xx tried once; 5xx tried four times, then refused; a held lock refused;
  a header CRS other than the catalogue's refused.
- **F7, the header**: a header longer than 1 MiB makes the reader double
  the prefix (the log shows 1, 2, 4 MiB); a read past the prefix raises in
  the test's stream; a prefix that cuts the full page's offset array is
  detected; past 64 MiB it is refused. A geographic header is read (decided
  6) while `read_page` without the flag still refuses it (23a-1 W9 green).
- **F8, GLO-30**: a tile missing from the list is "no tile" in the plan and
  the report, never requested; two listed tiles become two objects.
- **F9, mesh after fetch (K7)**: fetch the projected fixture for a domain
  whose box corners lie outside it, then `rasputin mesh --dem <key>` with
  `socket.socket` raising succeeds, with the same vertices and triangles as
  the same file meshed by path. Then the domain is grown past the fetched
  box, step by step, until the mesh's window first meets a block fetch did
  not take: that mesh is `NotCached` naming the fetch command, and its
  `missing` and `needed` equal the test's own counts, from `blocks_meeting`
  on the grown windows against the blocks fetched (not from the fetch's
  plan; growing on every side can add more than one block).
- **F10, the write side**: a manifest write that fails before `os.replace`
  (monkeypatched) leaves the previous manifest readable; `NOTICE.txt`
  holds each source's `credit`, `licence_note` and every `cite` entry
  verbatim; `--refresh` empties the source's cache and refetches; dry-run
  leaves `tmp_path` byte-identical (an empty cache stays without a source
  directory) and prints the counts the real run then reports.
- **F11, the boundary**: only `fetch/http.py` imports `urllib` or
  `http.client` (an AST scan over `src_python/tin_engine`); 23a-1's W7 stays green;
  `rasputin fetch` with neither `--cache` nor `RASPUTIN_DATA` is refused
  with B7's message.

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
- **SP1, a grid-line seam**: the check points are the nodes on it and the
  midpoints of the cell sides between them, `2c + 1` for `c` nodes (N14);
  after `refine_seam`, the error is within tolerance at every node, and at 1,000
  points along the line against a bilinear surface computed in NumPy
  independently (the every-point property).
- **SP2, a general seam**: an edge with rational endpoints; the oracle
  recomputes crossings, midpoints and bilinear heights in exact rationals
  (Python `fractions`) and checks every check point within tolerance
  + 1e-9 of the output's piecewise-linear heights; at tolerance 0 every
  check point ends at error 0 (to the slack); ties go to the smallest
  parameter.
- **SP3, sameness**: the edge given reversed gives the same output bit for
  bit; the strip window shifted by whole cells does too where the lattice
  arithmetic is exact in both windows, and agrees to rounding elsewhere
  (N17).
  In 23c, on the projected-mosaic path, a case where two pieces' windows
  select different tile sets around one seam: both pieces get the same seam
  points (read from the files, as DC2 does).
- **SP4, the bound**: insertions never exceed check points; NoData stencils
  are skipped and counted.

### 23c (Python, end to end on synthetic rasters)

- **DC0, the partition rule** (pure, no mesh): `b(T)` at the table's
  columns and rounded up between them (0.5 m gives 452); figures of
  `partition.py` (4,208 × 7,347 nodes at 1 m and `--pieces 64` gives 6 × 11
  cells of 702 × 668; 41,332 × 50,297 at 1 m gives 5 × 7 of 8,267 × 7,186,
  at 50 m 2 × 2; `--pieces 16` on the first gives 3 × 5, 15 cells: the
  count approximates the request and is not a floor); one piece at `N · b = B` and a cut at `N · b = B + 1`;
  a large `--memory-budget` gives one piece above every former cap; for
  windows, tolerances, `--pieces` and budgets drawn at random (a seeded
  generator, a few thousand) no cell's `dx · dy · b` exceeds the budget, no
  row or column of cells is
  empty, `dx` and `dy` are whole nodes and every line is a lattice line; the
  same partition with `os.cpu_count` and `physical_memory` patched to small
  and large values (neither the count nor the cut is the machine's).
- **DC1, K1**: a domain under the budget, and a larger one with a budget
  above its estimate, write the same bytes as today's path.
- **DC11, no machine refusal** (B15, ruled (a)): a mosaic
  over half of a patched `physical_memory` is planned, not refused; 15a's M5
  refusal cases in `test_mosaic.py` and `test_dem_input.py` are inverted
  (pinned behaviour B15 changes); the same for 22's catchment refusal in
  `test_catchment.py`.
- **DC2, conformity, read from the files**: for every seam edge, both pieces'
  vertex sequences are equal bit for bit. The oracle reads the piece files,
  not the seam records or the index.
- **DC3, a valid union**: in the union of the piece files (assembled in the
  test; the stitcher is 23d's) every edge not on the outline
  is in exactly two triangles, the triangles' area sums to the domain's
  (shapely) to 1e-9 relative, and no two triangles overlap (brute force on a
  small case).
- **DC4, constrained Delaunay with the seams**: an exact incircle oracle
  (Python `fractions`, independent of refine's predicates) finds no violation
  when the seams count as constraints, on a case with cocircular lattice
  points; and finds some when the seams are left out, so the test can fail
  and shows what the union is.
- **DC5, the tolerance over the union**: every DEM node in the domain within
  tolerance of the union of the piece files (brute-force barycentric location, 1e-9
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
- **DC7, determinism**: `--threads` 1 and 4, pieces in reversed order:
  identical piece files (`--jobs` is 23d's, PJ4; the reversed order
  moved there, "Settled after 23c's red step" item 12).
- **DC8, locality**: one piece re-run alone from its job spec writes the same
  bytes.
- **DC9, geometry**: two pieces in one cell; a partition line through a
  lake ring; a road along a partition line (one edge, both bits); a cell
  corner inside the domain (four pieces at one node, which the noder must
  return exactly); a sliver piece with no node.
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
- **PJ3, invalidation**: a changed tolerance changes every piece's hash; a
  changed `--pieces` changes the partition and so every hash.
- **PJ4, determinism**: `--jobs` 1 and 3, `--threads` 1 and 4, and pieces
  in reversed order (moved from 23c's DC7, "Settled after 23c's red step"):
  identical piece files and stitched file.

### 23f (C++ Catch2 and through the binding)

- **VR1, retire**: removing a vertex frees its slot and two triangle slots;
  the compacted output has no dead index, and the mesh's adjacency is valid
  (every edge in at most two triangles, consistent neighbours).
- **VR2, the hole is Delaunay**: for random stars (Delaunay meshes of a few
  hundred random and lattice points), removing each interior vertex in turn
  gives exactly the triangles a brute-force Delaunay of the remaining points
  has inside the star, checked by an exact incircle oracle in Python
  `fractions`; cocircular lattice stars accept any Delaunay answer (no
  violation, not equality); the mutant "fan from the first polygon vertex,
  no flips" fails.
- **VR3, collinear links**: a star with three or more collinear polygon
  vertices (a lattice line) is retriangulated with no zero-area triangle and
  no violation.
- **TH1, thinning keeps what it must**: on a ramp-plus-bump raster, a seam
  vertex at the bump stays and the ones on the plane go; with tolerance 0
  none goes; the order is the one given (reversing the candidates changes
  the result only where the test says it may).
- **TH2, the hole recheck is exact**: after thinning, a brute-force oracle
  (every node, barycentric location, 1e-9 slack) finds every node within
  tolerance; the mutant that checks only the removed vertex's own node
  fails on a case built for it (a node in the hole off the removed vertex).
- **TH3, the band limit**: a candidate whose star leaves the band is not
  removed and is reported as deferred.
- **FC1, the fence check**: an exact cocircular fence edge is not a
  violation, one just inside is; the result is the same with the edge given
  reversed.

### 23g (Python, end to end on synthetic rasters)

- **SC1, K1**: one piece, no cleanup, the same bytes as 23c's DC1.
- **SC2, Delaunay across seams** (K11): in the stitched file no edge has the
  `seam` bit, and the exact incircle oracle of DC4, with only the input's
  constraints, finds no violation; the same oracle on 23d's stitched file
  (seams kept) finds some, so the test can fail.
- **SC3, the tolerance after cleanup**: DC5's oracle on the stitched file,
  and on the reprojected path 15c's independent final check (G6), find 0
  over tolerance; the z-shift control fails.
- **SC4, thinning works**: on a case where the seam pass over-densifies (a
  plane crossed by a seam with a bump far from it), the stitched file has
  fewer vertices on the former seam than the piece files, and no vertex the
  tolerance needs is gone (SC3).
- **SC5, determinism and locality** (K5, K12): `--jobs` 1 and 4 and units in
  reversed order give identical stitched files; one unit re-run alone from
  its two band records gives the same bytes.
- **SC6, the final pass**: a case built so that one band's flips reach the
  fence (a long thin triangle across the band edge): the violation is
  deferred, the final pass clears it, and SC2 holds; a four-piece corner is
  removed where the tolerance allows and kept where it does not; a seam's
  end on the outline is kept; a four-piece corner built so that both
  triangles of one seam edge also have an edge on lower-numbered units: the
  stitched file has no `seam` bit there, the edge is locally Delaunay (SC2's
  oracle) and SC3 holds; the same with the corner on a feature (kept), so
  the bit is cleared whether or not the corner is removed.
- **SC7, memory shape**: the cleanup opens band records only, never a piece
  file, on a case with no deferral (a loader double that refuses piece
  files).
- **SC8, zones** (pure, from band records): on DC9's corner cases, every
  triangle with an edge on unit `u` is in zone(u) and band(u), except one
  with edges on two units, which is in the lower-numbered unit's; every edge
  of `u` is in band(u) or one of its fence edges. The mutant "zones by
  nearest centroid only" fails on a case built with a triangle on `u` whose
  centroid is nearer another unit.

## @perf acceptance

Per `docs/increments/README.md` "Acceptance": `tools/bench.py`'s 1 m benchmark
and thread-scaling sweep, `pmset -g batt` with every run, against the
previous increment's run in the same power state (on `NO BASELINE`, the
previous merge commit with `--tree`, back to back), evidence under
`docs/benchmarks/<date>/`.

- **23a-1:** the benchmark's mesh hash unchanged and process time within
  noise (decoding is now windowed). On DTM10 tiles: decode throughput per
  block and per window, threads 1 to 8 (ANADEM's waits for 15c-2, which
  lifts the geographic refusal).
- **23a-2:** fetch the basin's ANADEM box (about 8,300 blocks, "Decided
  here (23a-2)" 1): wall time, bytes, requests, `--dry-run`'s plan against
  the run, a rerun that fetches nothing, and an interrupted run resumed.
- **23b:** the README rule in full; the mesh hash unchanged (mask 0), refine
  within noise at every thread count.
- **23c:** the README rule; the benchmark is one piece at the default, so
  its hash is unchanged. The Velhas piece on ANADEM from the cache, at 1, 2,
  5, 10, 20 and 50 m, cut with `--pieces 64` (6 × 11 cells, 44 pieces) and uncut
  (the default; 15c-2's run):
  triangles and their difference, worst angle, maximum degree, 0
  constrained-Delaunay violations with seams as constraints, 0 nodes and 0
  source nodes over tolerance by the independent check (its control still
  failing), the seam pass's insertions, time and peak RSS. These are the
  before-cleanup figures 23g compares against.
- **23d:** scaling on the Velhas piece at 1 m over `--jobs` × `--threads`,
  cut with `--pieces 64`, against one piece on the same cores; the per-piece
  fixed cost (planning, slicing, windows, seam pass, writing) on that run's
  smallest pieces, which tells how far a `--pieces` request pays.
- **23f:** the README rule in full; the mesh hash unchanged (nothing in
  23f runs on the existing path); vertex removal timed per candidate on a
  synthetic band.
- **23g (B5's comparison, ruled):** the Velhas piece at the same six
  tolerances, cut and cleaned against uncut: **triangle count**; **vertices
  within two cells of a former seam**, against the uncut mesh's in the same
  strips; **worst angle in those strips** and the share of their triangles
  under 5°; 0 violations by the exact incircle oracle with the input's
  constraints only; 0 nodes and 0 source nodes over tolerance; the cleanup's
  time, the units' bands and how many items went to the final pass. **The
  gap is large**, and B5 (b), the exchange, goes back to Ola as the
  fallback, if at any tolerance the cut triangles exceed uncut by more than
  1 %, or the strip vertices by more than 25 %, or the strips' worst angle
  is below half the uncut strips' (thresholds decided here, B14).
- **The basin run** (B11, after 23g): the basin at 50, 20, 10, 5, 2 and 1 m:
  triangles against the 200-box estimate (237 M at 1 m, 215-261 M), wall
  time per phase (cleanup included), fetch bytes, peak RSS, and the stitched
  file's size. The claim to test: **peak RSS under 16 GiB at every
  tolerance**, i.e. memory no longer decides the tolerance. The independent
  final check over the union (on sampled pieces if a full check is too long;
  `@perf` says which and why). Then Ola chooses the tolerance.

## Settled after 23b's red step (45045dc)

Ruled by `@architect`, 2026-10-03, with Ola present. `@tester`'s red step
for 23b pinned 13 choices the design left open (N1-N13) and found six places
where the design was written before 15f-1 and 15f-2 landed and is now out of
date (N14-N19). Each is confirmed or corrected here before `@developer`
starts, so that no red-step choice reaches green unexamined. They are part of
the design, and `@developer` implements them as written here. The pins are
quoted from the headers of `tests/cpp/unit/test_mesh_frozen.cpp`,
`tests/cpp/unit/test_refinement_scan_frozen.cpp`,
`tests/cpp/property/prop_refinement_frozen.cpp`,
`tests/cpp/property/prop_refinement_seam.cpp`,
`tests/python/test_core_frozen.py` and `tests/python/test_core_seam.py` at
`45045dc`.

**Verdict in one line:** all 13 pins stand. N5, N6, N7 and N10 are confirmed
with a precision the tests do not contradict; N16 adds one refusal that
`@tester` must write red (a frozen strip edge); nothing in the committed
suites has to change.

### The 13 pins

- **N1. The mask lives on the mesh. Confirmed.**
  `LatticeMesh::set_frozen_mask(std::uint32_t) noexcept`,
  `frozen_mask() const noexcept` (0 after `build`) and
  `is_frozen(t, e) const noexcept`, true exactly when
  `(mask(t, e) & frozen_mask()) != 0`. Holding it on the mesh is right: the
  scan, the quality pass and `split_edge`'s assertion (FE6) all receive the
  mesh already, so no signature changes and `improve`'s stays as it is.
  `refine`, `refine_points` and `refine_strip` call `set_frozen_mask` once,
  after `to_lattice` and before `legalise_all`.
- **N2. `frozen_mask` is the last member of `RefineOptions` and
  `PointRefineOptions`. Confirmed.** Last, so every existing aggregate and
  designated initialiser compiles unchanged (K1). `refine_strip` takes
  `PointRefineOptions`, so it reads the same member; nothing else is added
  for it.
- **N3. The scan and a node on a frozen edge. Confirmed.** A node exactly on
  a frozen edge of the triangle (exact orientation, inside the closed
  segment, not an end) is never the argmax, never the carve point of a void
  triangle, and not counted in `uncovered`. Only triangles with a frozen edge
  take that path (FS3 pins it). The carve point then is the nearest valid
  node off the frozen edge; a void triangle whose only valid nodes lie on its
  frozen edge carves nothing, and those nodes are the seam pass's (its
  greedy measures them, and N10's carving covers them when a seam end has
  no z). The "not in `uncovered`" half is what keeps refine's
  stopping rule from waiting on nodes it may not insert.
- **N4. `skipped_frozen` is summed into `quality_skipped`. Confirmed.**
  `RefineOutcome::quality_skipped` already reads "every reason summed";
  `QualityOutcome::skipped_frozen` is a new reason. A bad triangle whose
  snapped node lies exactly on a frozen edge of the triangle the walk found
  is counted there and nothing is inserted.
- **N5. Feet and frozen edges. Confirmed, with the rule written out.** In
  `detail::foot_of`, a frozen edge is skipped as an unconstrained one is
  (`continue`), so it is never a foot's edge. If another, non-frozen
  constrained edge of the triangle qualifies under 20b's own rule, its foot
  is taken as today; otherwise there is no foot and the node goes in itself.
  That is not a refused foot (none was computed), so `feet_refused` does not
  count it. FE3's fixture has one candidate edge, so the pin and this rule
  agree.
- **N6. `on_frozen` and `on_frozen_max_error`. Confirmed, with three
  precisions.**
  - It counts **stored check points** (the source set of `refine_points`)
    that lie on a frozen edge, each point once, whatever the number of
    triangles that hold it and whatever the thread count. How it is counted
    once is `@developer`'s choice: an end pass over the store, or the scan
    counting only on frozen edges the triangle owns (15f D4's ownership
    rule: lower to higher vertex index in this triangle, or no triangle
    across), summed over each slot's last result as L4 does for `max_error`.
  - The error is `|z − (z_a + σ (z_b − z_a))|`, with `σ` the projection
    parameter of the point on the frozen edge in lattice `(col, row)` and
    `z_a`, `z_b` the loop's vertex z of its two ends. Frozen edges are never
    split, so the edge in the output is the edge in the start.
  - A point on a frozen edge with an end that has no z is counted and adds
    nothing to `on_frozen_max_error` (there is no linear z to compare with),
    as `coincident_max_error` reads "valid vertices only". It is never a
    carve point (N3's rule, in the point scan).
- **N7. L12 next to a frozen edge: skipped and counted in `on_frozen`.
  Confirmed (`[pinned]` case in `prop_refinement_frozen.cpp`).** With a
  strip given, the point scan treats a point within the coincidence radius
  `r(g)` (15f L16) of a frozen edge of its triangle, projection strictly
  inside it, as N3 treats a point exactly on it: skipped, never named. This
  is a test inside the read-only scan, like L14's corner test, so it needs no
  marks and cannot make a triangle name the same point every round. A frozen
  edge is therefore never an L12 candidate; there is no `split_inside`
  fallback beside a frozen edge, which would leave the sliver F1 and L14
  describe. The order inside the scan: L14's corner test first (the point
  counts in `coincident`), then this frozen test (it counts in `on_frozen`,
  with N6's error at its projection), then the scan as today. Without a
  strip, the frozen test is exact (radius 0), as L12 and L14 are; K1 and
  15c's path are unchanged.
- **N8. `SeamPoint`, `SeamOutcome`, `refine_seam`. Confirmed.**
  `refine_seam(dem, a, b, tolerance)` in
  `include/terrain/refinement/seam.hpp`; `SeamPoint{Point2 at; double z;
  double s;}`; `SeamOutcome{Point2 a, b; std::optional<double> z_a, z_b;
  std::vector<SeamPoint> points; std::size_t check_points, no_data; double
  max_error;}`, with `a < b` by `(x, y)` whichever order the caller gave.
  An inserted point is output at `(x_min + col·dx, y_max − row·dy)` of its
  lattice position with the check point's own z, as L6 outputs a strip point.
  `z_a`, `z_b` are `vertex_z` at the ends' lattice positions, `nullopt` where
  it refuses. The pass passes the generator the vertices `{a, b}` in that
  order, so the generator's canonical `P0` (15f D2, step 2) is `a` and `s`
  runs from `a`.
- **N9. The check points are `constraint_check_points`' for the one edge.
  Confirmed.** That is the reuse 23 designed ("What survives of 15c, the
  edge strip and 15d"). `check_points` and `no_data` are the store's
  `size()` and `no_data()`. The consequences on grid-line seams and short
  seams are N14 and N15.
- **N10. An end on NoData: carve along the edge. Confirmed (`[pinned]` case
  in `prop_refinement_seam.cpp`).** A piece of the polyline with an end that
  has no z has no lerp. While such a piece holds a check point, the one
  nearest its invalid end (smallest `|s − s_end|`, ties to the smaller `s`)
  is inserted. This is 15f D4's void-sub-edge rule, applied in one
  dimension, and it is needed: refine may not carve on a frozen edge (N3),
  so without it the seam's valid nodes next to a NoData end would be covered
  by no one. Because NoData check points are already dropped, the inserted
  point is valid, and each invalid end costs at most one insertion. The set
  inserted does not depend on whether carving runs before or after the
  greedy; `@developer` may do it first.
- **N11. `max_error` leaves out void pieces. Confirmed.** After N10 a void
  piece holds no check point, so there is nothing to leave out at the end;
  the rule matters only for what the greedy compares while it runs.
- **N12. Refusals: `std::invalid_argument` naming `refine_seam`.
  Confirmed.** For a tolerance that is negative or not finite (the text also
  contains "tolerance"), for `a == b` (exact equality of the world points),
  and for an end outside the node rectangle (the text contains "outside").
  `refine_seam` checks these itself before calling the generator, so the
  message names `refine_seam` rather than `constraint_check_points`. It
  throws where `refine` returns a status because it is a pure function with
  no outcome status, as the generator is; pybind11 turns
  `std::invalid_argument` into `ValueError`.
- **N13. The Python surface. Confirmed.**
  `_core.refine_seam(view, a, b, *, tolerance) -> SeamOutcome`; `.a`, `.b`
  tuples; `.z_a`, `.z_b` float or `None`; `.points` `(K, 2)`, `.z` and `.s`
  `(K,)`, all float64 and read-only; `.check_points`, `.no_data`,
  `.max_error`. `frozen_mask=0` as a keyword on `refine` and
  `refine_points`; a negative mask is `TypeError` (pybind11's own refusal of
  a negative for an unsigned parameter, nothing written for it);
  `PointRefineOutcome.on_frozen` and `.on_frozen_max_error`; stubs in
  `_core.pyi`. Added here: `refine_seam` releases the GIL, as every refine
  binding does, because 23d runs pieces on threads (`asyncio.to_thread`)
  and each runs its seam passes.

### What 15f changed in 23b's design

- **N14. A grid-line seam gets `2c + 1` check points, not the nodes
  alone.** The generator adds the midpoint between each two neighbouring
  crossings, ends included; on a grid line those are the midpoints of cell
  sides. **Accepted as is**: filtering them would be a second code path for
  no gain. Along a grid line both the bilinear surface and the polyline are
  linear between nodes, so a midpoint's error is the mean of its two
  neighbours' errors. It is never strictly the worst, and on a tie the node
  before it wins (smaller `s`). So the greedy inserts nodes only, apart from
  a tie broken the other way by rounding, in which case the midpoint goes in
  at its exact half-integer position with `vertex_z`'s z. The every-point
  property (SP1) holds either way. "Exact heights" under "Why lattice lines"
  now reads "nearly always a node"; K4 does not depend on it, since both
  pieces write the record's z (step 5).
- **N15. A seam edge shorter than a cell.** With no crossing it has exactly
  one check point, its midpoint, and none if the cell holding the midpoint
  touches NoData (`vertex_z`'s rule). The degeneracy policy's "may have no
  check point" is reworded to say so (SP5 pins both).
- **N16. The paths onto a seam, after 15f-2.** K2's list was written before
  `refine_strip` existed. The paths now are refine's split, feet, the
  quality pass, `refine_points`' source points, and in 15f's loop the strip
  points, the DEM rescan of `refine_strip`, and L12's on-edge insertion.
  Each is closed:
  - the DEM rescan calls refine's own `scan`, which reads the mesh's mask
    (N1, N3);
  - L12 never takes a frozen edge (N7);
  - strip points: the caller builds the strip on non-frozen edges only, and
    **a strip edge that is frozen is refused**, like L2's other programming
    errors: `std::logic_error`, text starting with the entry point's name
    (`"refine_strip: ..."`, `"refine_points: ..."`), checked after L2's
    check (3) ("not a constraint edge of the start") and before (4).
    Refusing rather than skipping keeps a wrong caller loud; skipping would
    hide a strip that silently checks less than its caller thinks.
  - consumption (15f D4, step 4) acts on strip sub-edges only, and a frozen
    edge has none.

  **`@tester` adds, red:** one case each for `refine_strip` and
  `refine_points(..., strip)`: a strip built over every constraint edge,
  the frozen one included, run with that edge's mask in `frozen_mask`, is
  `std::logic_error` whose text starts with the entry point's name; the same
  strip with `frozen_mask = 0` runs (the control). This is the only change
  to the suites this section asks for.
- **N17. SP3, reworded.** Lattice coordinates are measured from the corner
  of the raster C++ is handed (15f L16), so a crossing computed in the strip
  window and the same crossing computed in a window shifted by whole cells
  can differ in the last bits. The claim is therefore: **the edge reversed
  gives the same output bit for bit; a shifted window gives the same output
  bit for bit where the lattice arithmetic is exact in both windows** (ends
  on nodes, or crossings that are exact binary fractions in both frames),
  and otherwise agrees to rounding. `@tester`'s SP3 case already chooses
  ends that are exact in both frames and says why, so the test stands. What
  the design relies on is not a shifted window but the same window: the
  strip is a function of the edge alone (step 2), so both pieces call
  `refine_seam` on the same raster and get the same bits. On 23c's
  lattice-line seams the arithmetic is exact anyway (N19 below).

  One consequence for 23c, written here so it is not rediscovered there: a
  vertex at the end of two or more seam edges (an outline crossing, a
  feature crossing, a cell corner) gets a `z_a`/`z_b` from each edge's
  record, computed in different strip windows. Off a node these can differ
  in the last bits. **23c's rule: such a vertex takes its z from the record
  of the seam edge with the lowest index in the start triangulation's edge
  list**, a function of the input, the same in every piece, so K4 holds at
  corners too.
- **N18. If 15f-3 lands first.** 15f-3 binds `refine_strip` and
  `refine_points(..., strip)`. Whichever of 15f-3 and 23b lands second adds
  `frozen_mask=0` as a keyword to the `refine_strip` binding and its stub,
  with one binding test (a frozen outline is not split; a negative mask is
  `TypeError`), mirroring `test_core_frozen.py`. 15f-3's
  `edge_strip.generate` takes the edges as given; leaving frozen edges out
  of the strip is the caller's job, and the only caller with a nonzero mask
  is 23c, so the filter (`edges[(masks & frozen) == 0]` in NumPy) is 23c's
  line, not 15f-3's. With 23b's refusal (N16) a missing filter fails loudly.

  **Settled after N18's red step (e26802e).** `@tester`'s four choices:
  1. *Order:* confirmed. `frozen_mask` is the last keyword, after
     `threads`, as on `refine_points`; the binding and the stub follow the
     kw-only pin in `test_core_edge_strip.py`.
  2. *The frozen test:* confirmed, with one amendment. The strip is built on
     non-frozen edges (the caller's job, N18) and the refusal case shows the
     mask reaches the N16 check. That does not show the mask reaches the DEM
     rescan's skip (N1, N3): if the filtered strip with `frozen_mask=0`
     leaves the west side unsplit anyway, the "not split" assertion passes
     with the mask ignored by the rescan. So it must be measured: `@tester`
     runs the filtered strip with `frozen_mask=0` and records how many
     vertices land on the west side. If more than zero, that becomes a
     second control in `test_a_frozen_side_is_not_split`, asserted `> 0`. If
     zero, the scene cannot see the rescan's skip; the test's docstring says
     so, and that skip stays pinned by refine's own frozen tests (`scan`
     is shared, N16), not by a new scene here.
  3. *The refusal case:* confirmed. `std::logic_error` reaches Python as
     `RuntimeError` (as L2's other refusals in `test_core_edge_strip.py`),
     and `^refine_strip: ` is N16's text rule. The control N16 asks for (the
     same full strip with `frozen_mask=0` runs) is
     `test_frozen_mask_0_is_todays_result`.
  4. *Location:* confirmed. `test_core_frozen.py`, beside the other two
     entry points' `frozen_mask` classes.
- **N19. `on_frozen` for `refine_strip`'s DEM nodes: not counted.**
  `on_frozen` counts stored check points only (N6), so it is 0 in
  `refine_strip`, whose stored set is empty. The reasons:
  - a DEM node exactly on a frozen edge is a seam check point (the
    generator snaps node crossings exactly, 15f D2 step 4), so the seam
    pass's guarantee (K3) covers it, and refine does not count it either;
  - L3 limits the DEM rescan to triangles the run has written, so a count
    of DEM nodes near frozen edges would depend on which triangles were
    written, and would not mean anything.

  **The limit this leaves, stated.** On a seam that is not a lattice line
  (a "general seam"), the seam pass's points between two crossings are a
  hair off the original line, so a DEM node that lies exactly on the
  original line but was not inserted can be a hair off the fan sub-edge
  that now holds it. With a strip given, N7's radius skips it, uncounted;
  in refine (no strip, radius 0) it is scanned as an inside node, its error
  is the seam pass's to rounding, and only at a rounding tie with the
  tolerance would it be inserted, a hair from a frozen edge, which K2's
  oracle would report. On a lattice-line seam this cannot happen: the seam
  lies on a column `K` or row `R` that is exact in every window
  (`x = K·h` in whole metres when `h` is, 23c's DC9), every crossing is a
  node, and orientation is exact. **23c makes lattice-line seams only**
  (natural cuts dropped with B1), so 23b's guarantees are exact on every
  seam the run makes. A future general seam (23e, or a cut along an input
  line) reopens this, and the fix then is N7's radius in refine's scan too,
  counted.

### `@tester`'s throwaway implementation

`@tester` built a throwaway implementation in scratch, not committed, to
check that the suites can pass. Against the rules:

- **Allowed, narrowly.** The lean-brief rule (no throwaway, no mutation
  round by default) has one exception: a suite the increment file marks
  invariant-critical. 23b's FE2-FE5 and SP1-SP2 are so marked ("Tests
  `@tester` can write red"). So a throwaway is not outside the rule.
- **But it was used for a different purpose.** The exception is for
  mutation testing, which these suites' own headers place at green
  ("mutation runs at green"). A throwaway built at red to check
  feasibility is not that; it is the cost the rule was written to avoid,
  and its choices can leak into the pins. Here they did not do harm: every
  pin is confirmed above on its own grounds, not on "the throwaway passed".
- **It is not a hand-off.** Scratch is not a channel
  (`.claude/REQUIRED-READING.md`, "Data, scratch and temp folders").
  `@developer` writes 23b from this file and the suites, does not read the
  throwaway, and the main session does not name its path in a brief. Next
  time, the brief should say whether a throwaway is wanted.

### LOC effect

N7 (the frozen test with the radius in the point scan), N10 (carving), N16
(one refusal) and N13's arrays add about 40 lines. The 23b rows of "PR split
and LOC" are updated: about **285**, 396 at +39 % and 456 at +60 %, under
700 at both.

### As built (23b green)

Green is `3c464ec`. `@tester`'s `f00a7b1` closed the two gaps the mutation
pass found and amended 15c's stub test, which pinned `refine_points`'
keyword-only arguments without `frozen_mask` (N13). LOC in `CLAUDE.md` §2's
unit: 258 added and 22 removed, 236 net, against about 285. Choices the
rulings left to `@developer`, and other changes outside the PR table:

- **`on_frozen` once per point (N6), by ownership in the scan.** No end pass
  over the store. A stored point exactly on a frozen edge is counted by the
  triangle that owns the edge (15f D4: lower to higher vertex index in this
  triangle, or no triangle across). A point only within `r(g)` of the edge
  (N7) lies in one triangle, which counts it. Both are summed over each
  slot's last result, as `max_error` is. `PointScan::offer` keeps the count
  when a strip or DEM candidate replaces the scan's winner.
  **A rounding-scale double count, accepted as a report.** A point counted
  for being within `r(g)` of a frozen edge, not on it, is counted by every
  triangle that holds it (closed membership). A point lying exactly on an
  interior edge from a frozen edge's end vertex `v`, at an angle `θ` to that
  frozen edge, is in both triangles beside the interior edge. If each of the
  two has a frozen edge at `v` whose projection test it passes, the point is
  counted twice. That needs a bend in the frozen chain at `v` (two
  collinear frozen edges cannot both hold the projection strictly inside)
  and a distance from `v` between `r(g)` (below it, L14's corner test skips
  the point) and `r(g)/sin θ`. Bends occur at every seam-pass point of a
  seam that is not a lattice line, whose sub-edges are a hair off
  collinear; on lattice-line seams only at seam corners, where the window
  is at most `r(g)` to `r(g)·√2`. `on_frozen` is a report, not an invariant, in the
  same class as `coincident`'s rounding-scale double counts
  (`refine_points.hpp`, the comment above the coincident pass). If only one
  of the two triangles has the frozen edge, the other may name the point and
  split the interior edge there, putting a vertex within `r(g)` of the
  frozen edge, which K2's oracle would report. That also needs exact
  incidence with the interior edge at that distance, so it is rounding-scale
  and is not fixed in 23b; it is the same gap N19 states for general seams.
- **N16 is a second loop.** The frozen-strip-edge refusal runs after L2's
  check (3) has passed for every strip edge, not inside the same loop, so
  a strip that also holds a non-constraint edge is refused for that first,
  as N16 orders.
- **The seam greedy keeps each piece's worst point in a priority queue.** It
  inserts the same set as the naive greedy, since pieces are independent, and
  ties go to the smallest `s` inside a piece. Cost is O(n log n · depth), not
  O(n²).
- **L12's guard is kept though the scans make it unreachable.** Both scans
  skip a point within `r(g)` of a frozen edge, by the same distance and
  projection `near_constraint` computes. But the two copies of that
  expression may be FP-contracted differently, so `near_constraint` also
  excludes frozen edges.
- **CI:** `prop_refinement_frozen` and `prop_refinement_seam` join the TSan
  job in `.github/workflows/main.yaml`.
- **Citations:** three line citations in `15f-edge-strip.md` (:61, :84,
  :355) moved to where `scan_points`, `scan` and `vertex_z` now are.
- **K1:** a scratch program hashes `refine`, `refine_points` (with and without
  a strip) and `refine_strip` over 756 seeded runs. It gives the same hash
  with master's headers and with 23b's when `frozen_mask` is not named or is
  0, both with the default FP contraction and with `-ffp-contract=off`. A
  planted change to the scan changes the hash.

**Mutation runs** (FE2-FE5, SP1, SP2; Release build, one mutant at a time,
the source restored and touched after each). `@developer`'s 27 ran against
`3c464ec`'s suites, and `@tester`'s list against `f00a7b1`'s. Where
`@tester`'s mutant was the same as one of the 27, it is not run twice.
The suites that kill a mutant are `test_mesh_frozen` (mesh),
`test_refinement_scan_frozen` (scan), `prop_refinement_frozen` (frozen) and
`prop_refinement_seam` (seam).

| `@developer`'s mutant | verdict |
|---|---|
| M1 the scan never skips a frozen node | killed (scan, frozen) |
| M2 the frozen skip only in the all-node path | killed (scan, frozen) |
| M3 the frozen skip only off the all-node path | killed (scan, frozen) |
| M4 `frozen_edge_at` does not exclude the edge's ends | survived; equivalent: `for_each_row_span` excludes the triangle's vertices (`include/terrain/mesh/row_spans.hpp`), so the scan never visits a corner, and the point scan's corner test returns before `frozen_edge_at` is called. ("A vertex has error 0" was the first reason given; it is not exact at ulp level, so it is not the reason.) |
| M5 the radius ignored (N7) | killed (frozen) |
| M6 `is_frozen` true for any masked edge once a mask is set | killed (mesh, scan, frozen) |
| M7 `refine` does not set the mask | killed (frozen) |
| M8 `point_loop` does not set the mask | killed (frozen) |
| M9 feet taken on frozen edges | killed (frozen) |
| M10 the quality pass splits frozen edges | killed (mesh, frozen) |
| M11 `skipped_frozen` not summed into `quality_skipped` | killed (frozen) |
| M12 L12 may take a frozen edge | survived; unreachable while the two distance expressions agree (above) |
| M13 `scan_points` has no frozen test | killed (frozen) |
| M14 `on_frozen` counted in both triangles | killed (frozen) |
| M15 the error measured from the wrong end | killed (frozen) |
| M16 a frozen strip edge not refused | killed (frozen) |
| M17 `on_frozen` not summed over slots | killed (frozen) |
| M18 `offer` drops the frozen count | survived at `3c464ec`; killed (frozen) after `f00a7b1` |
| S1 ties between pieces to the largest `s` | survived; equivalent (pieces are independent) |
| S2 inserts at `error >= tolerance` | survived at `3c464ec`, C++ and Python; killed (seam) after `f00a7b1` |
| S3 no carving from the end `b` | killed (seam) |
| S4 lerp by index, not by `s` | killed (seam) |
| S5 ends not ordered | killed (seam) |
| S6 `max_error` reported as 0 | killed (seam) |
| S7 the right piece not re-offered after a split | killed (seam) |
| S8 output x scaled by `dy` | killed (seam) |
| S9 ties inside a piece to the largest `s` | killed (seam) |

| `@tester`'s mutant | verdict |
|---|---|
| M3 the void branch drops the frozen skip | ran: killed (scan, frozen) |
| M4 skipped nodes still counted in `uncovered` | ran: killed (scan, frozen) |
| M5 the scans skip nodes on every constrained edge | ran (both scans, since `frozen_edge_at` is shared): killed (scan, frozen) |
| M6 `is_frozen` returns `is_constrained` | ran: killed (mesh, scan, frozen) |
| M9 `foot_of` refuses the whole triangle when any edge is frozen | survived at `f00a7b1`; killed (frozen) after `c7c225f`. Not equivalent: N5 takes a foot on another qualifying constrained edge, and FE3's first fixture has one candidate edge. `c7c225f` adds an FE3 case with a frozen side and the needle on a second, unfrozen constrained side: the foot must land there; with M9 planted all four of its generator runs fail and no other case does (confirmed by `@tester` and `@reviewer`) |
| M14 `on_frozen` once per triangle, so each point twice | covered by `@developer`'s M14: killed |
| M15 the frozen error measured against 0 | ran: killed (frozen) |
| M17 L12 puts a point onto a frozen edge | covered by `@developer`'s M12: survived, unreachable (above) |
| M18 a frozen mask of 0 freezes everything (K1) | ran: killed (mesh, scan, frozen) |
| S3 interpolation between the two ends only | ran: killed (seam) |
| S4 nodes only, no midpoints | ran: killed (seam) |
| S6 ends ordered as given | covered by `@developer`'s S5: killed |
| S7 no carving from either end | ran (both ends; `@developer`'s S3 removed one): killed (seam) |
| S8 inserted points output with the wrong origin | ran (`x_min` dropped): killed (seam) |

## Settled after 23c's red step (65e3990)

Ruled by `@architect`, 2026-10-04, unattended (Ola asleep). `@tester`'s red
step for 23c (`aacf37c..65e3990`) pinned 15 choices the design left open and
asked three questions. Each is confirmed or corrected here before
`@developer` starts; `@developer` implements them as written here. The pins
are quoted from the headers of `tests/python/test_decompose_partition.py`,
`test_core_mesh_arrays.py`, `test_mesh_index_conformity.py`,
`pieces_fixtures.py` and `test_cli_mesh_pieces*.py` at `65e3990`.

**Verdict in one line:** all 15 pins stand, three with a precision (2, 8,
14); `@tester` writes four more red tests (2, 8, 14 and question B); one
question goes to Ola (B17, the budget's unit), with a default that needs no
test change; 23c is split into two PRs.

### The 15 pins

1. **`_core.indexed_mesh(vertices, triangles, constrained_edges)`.
   Confirmed.** Copies, `ValueError` on a bad shape, index, mask or
   coordinate, no orientation check (refine's `NotCounterClockwise` is the
   check the degeneracy policy names). Edge-property masks stay out of it:
   `refine` already takes `edges` and `masks` beside the mesh, so the seam
   bit reaches the core the way every feature bit does. The binding was
   missing from the PR table; it is a row now (below, about 30 lines).
2. **Pieces in `<out>.pieces/`. Confirmed, with three precisions.**
   - The directory is `--out`'s full name plus `.pieces` (`x.vtk` gives
     `x.vtk.pieces/`), so a `.vtk` and a `.ply` run of one name do not share
     it. Piece files are `<j>-<i>-<k>.<ext>`; the index's `file` is relative
     to the directory, as pinned.
   - **23c writes nothing at `--out` in a cut run** (the stitcher is 23d's),
     and says on stderr where the pieces are. `@tester` asserts it; 23d's
     PJ1 inverts it.
   - **A rerun over an existing directory deletes its `index.json` first and
     writes the new one last** (to a temporary name, then renamed). Without
     this, a rerun that fails part-way leaves the previous run's index
     beside piece files it has overwritten: a valid-looking index of a mesh
     that never existed. Piece files not listed in the index are not the
     run's; 23c deletes nothing else. `@tester` adds one red test: a cut run,
     then a rerun that fails (a truncated DEM, as in the input suite), leaves
     no `index.json`.
3. **Index keys. Confirmed**, and three more fields the design's "Output"
   names, which the pinned tests do not forbid: `source` (the run's
   `elevation_source` sentence, identity and credit), `vocabulary` (the
   piece files' fingerprint, item 4) and per piece `window` (`row0`, `col0`,
   `rows`, `cols` on the lattice). No test change.
4. **`seam` named in the piece files' vocabulary. Confirmed; the bit is 9.**
   See question A.
5. **`MeshIndex` frozen, `extra="forbid"`. Confirmed.** An index is read by
   other programs (23d's stitcher, a consumer); an unknown key is drift, not
   data.
6. **`SeamRecord`, `check_conformity`, `ConformityError`. Confirmed**, NaN
   equal to NaN and 0.0 unequal to -0.0 included (bit for bit means bytes).
   One record passes: a seam edge on the outline has one piece (degeneracy
   policy). The runner knows from the start triangulation how many triangles
   each seam edge has (one or two) and asserts that many records exist
   before calling the check; that is internal and needs no test.
7. **`--memory-budget` in whole bytes, default 2^34; refusals. Confirmed for
   23c**, pending B17. `partition()`'s two `ValueError`s stand: below one
   node's bytes step 3 cannot end.
8. **`--pieces` and `--memory-budget` need `--dem` and `--tolerance`.
   Confirmed. Cutting needs `--domain`: both options without `--domain` are
   a usage error ("needs --domain").** Without a domain the start mesh is
   the stride grid, not a PSLG through `build_pslg`, `node` and
   `triangulate`, so there is nothing for the cuts to enter; cutting it
   would be a second start path. A run without `--domain` (or with `--bbox`)
   is never cut, keeps today's bytes (K1), and the budget does not apply
   to it. `@tester` adds one red test: `--pieces 4` without `--domain` is
   exit 2 naming `--domain`, with no file written.
9. **K1 on the uncut path, no partition field. Confirmed.** The budget and
   `b(T)` are recorded in the index, and only a cut run has one; K5's
   "recorded in the file" reads as "in the index". A run is cut when the
   partition has more than one cell (`nx · ny > 1`), not when the labelling
   finds more than one piece: a domain inside one cell of a multi-cell
   partition writes a pieces directory with one piece.
10. **DC8 as locality under a DEM change in another cell. Confirmed.** It
    tests K6 through what a job may depend on, which is the property. It
    also pins that a piece file's fields are the piece's own (its counts,
    its sentence), never the whole run's.
11. **`counts.on_frozen` per index entry. Confirmed.** `counts` also holds
    the piece's `triangles` and `vertices`.
12. **DC7 with threads forced through a wrapper; the reversed-order half not
    written. Confirmed.** 23c runs pieces in one fixed order and offers no
    hook; the order becomes a parameter with 23d's runner (largest first),
    so the reversed-order half moves to 23d's PJ4.
13. **DC6 on an asymmetric octagon on nodes, the start taken by a spy on
    `triangulate`. Confirmed.** The octagon keeps every crossing on a node,
    which is what the equality needs; the spy makes the oracle's start the
    run's own, so no CDT tie can separate them.
14. **DC9's lake as a hole in the domain. Confirmed, and not enough.** The
    hole is worth keeping (the partition must not lose a ring). But a lake
    as a water polygon (`--features`) is the case the basin is made of: a
    seam crossing a constraint ring, where the crossing is an off-node
    corner, the only seam vertices whose z comes from a bilinear evaluation
    rather than a node (seam protocol step 5). `@tester` adds the water
    variant: a partition line through a water ring given by `--features`,
    the ring not dividing pieces, the crossings vertices of both pieces,
    conformity and both oracles over the union.
15. **Area to 1e-6 for moved domains. Confirmed**, on the measured 4e-7 of
    the uncut mesh.

### The three questions

**A. Where `seam` lives. A piece-file vocabulary, not `DEFAULT_VOCABULARY`.**
`features.py` gains `PIECE_VOCABULARY`: `DEFAULT_VOCABULARY`'s properties
plus `seam` at bit 9, the lowest bit the default leaves free. Cut runs write
their piece files (and later 23d's stitched-with-seams file) with it; the
uncut path keeps `DEFAULT_VOCABULARY`, so K1's bytes, every fingerprint and
`test_features.py` stay as they are. Two more reasons: `feature_input`
builds its classes from `DEFAULT_VOCABULARY`, so a seam there would let a
user tag their own lines as seams; and after 23g's cleanup no `seam` bit is
left in the stitched file (K11), which can go back to the default. The
index records the piece vocabulary's fingerprint (item 3). No test change:
the suites read the bit by name.

**B. A piece whose triangles all fall over NoData. Keep it in the index,
write no file.** Its entry has `"file": null`, `sha256` null, zero
triangles and vertices, and its seam records still go to the conformity
check (the seam pass ran on its edges). Refusing would make a cut run fail
where the same run uncut succeeds; writing a zero-triangle file asks every
reader to handle one; dropping the entry hides that the piece was there.
The run is refused, as uncut ("no data under any triangle; nothing to
write"), only when every piece is empty, and then no index is written.
`@tester` adds two red tests: a cut where one cell is all NoData (exit 0,
that entry as above, the other files present and passing the oracles), and a
cut of an all-NoData DEM (exit 2, no `index.json`). Under NoData a seam edge
can lose its triangle on one side and not the other, so the edges suite's
"exactly two files" holds only where both sides keep theirs; its scene
avoids the case, and no change is asked.

**C. Two claims that cannot run before green.**
- *DC4's control* (violations once seams stop counting as constraints). It
  is evidence that the oracle can fail, and that the frozen seam costs
  Delaunay quality across it, which is what 23g removes. If at green the
  scene shows none, the test is not wrong and the product is not wrong:
  `@tester` makes the control a planted one (one interior edge of the union
  flipped, which the oracle must catch), and "As built" records that the
  scene's seams cost nothing. A scene that can show the cost is a 23g
  concern (its comparison of cut against uncut), not 23c's.
- *DC6's equality.* If it fails at green, `@developer` reduces it to the
  smallest case and names the first differing triangle. If the cause is in
  23c's code (slice order, renumbering, fan order, window origin), it is a
  bug and gets fixed. If the cause is in refine itself (a decision that
  depends on a triangle's index or on the window's origin rather than on
  the triangle), 23c does not change refine: DC6 is reduced to what the
  design needs (DC2-DC5 on the same scene, which K3 and K4 rest on), "As
  built" records the cause, and "The union argument, checked" item 1's
  "sidesteps the difference" is corrected. It then goes to Ola as a
  question, because "a cut run is one refine, restricted" is a claim a
  publication would make.

### PR split and LOC

The rulings add about 50 lines: the binding (~30, `bindings/core.cpp` and
`_core.pyi`), `PIECE_VOCABULARY` (~5), the index's write order and the empty
piece (~10), the `--domain` refusal (~5). That is about **480**: 667 at
+39 % and 768 at +60 %, over the ceiling at the second. **23c is split in
two stacked PRs** along the red step's own files:

| PR | what | est. | +39 % | +60 % |
|---|---|---:|---:|---:|
| **23c-1** | `decompose.py` (rule, `b(T)`), `features.py`'s `PIECE_VOCABULARY`, `io/mesh_index.py` (model, records, conformity), `_core.indexed_mesh` | 215 | 299 | 344 |
| | tests: `test_decompose_partition.py`, `test_core_mesh_arrays.py`, `test_mesh_index_conformity.py` | | | |
| **23c-2** | `decompose.py` (lines as chains, `BasinPlan`), `pieces.py`, `basin_run.py`, the index writer, `cli.py`, the refusals at half of physical memory deleted | 265 | 368 | 424 |
| | tests: `pieces_fixtures.py`, `test_cli_mesh_pieces*.py`, DC11 in `test_mosaic.py`, `test_dem_input.py`, `test_catchment.py` | | | |

DC11 goes with 23c-2, so the memory refusals are deleted in the same PR that
brings the cut: between the two merges a run over half of the machine's
memory is still refused rather than running uncut out of memory. The red
step's commits split along the same line (`edd8a29`, `d994a0e`, `db6a5c0`
to 23c-1; the rest to 23c-2); how the branches are rebuilt is the main
session's.

### Three points from 23c-1's green (3ffae13)

Raised by `@developer`; ruled by `@architect`, 2026-10-04.

1. **"The lattice" in the partition record: named, in 23c-2.** The index
   must say where the partition lines are, or a consumer cannot place a
   seam. Four fields join `PartitionRecord`: `row0` and `col0` (the
   window's first node `(R0, K0)` on the lattice, integers, possibly
   negative), `spacing` (`h`, metres) and `origin` (`[x, y]`, the world
   position of lattice node `(0, 0)` in the index's `crs`). Line `i` is then
   at `x = origin_x + (col0 + i·dx)·h`, line `j` at
   `y = origin_y − (row0 + j·dy)·h`. They belong to 23c-2, because
   `decompose.partition` sees only `cols` and `rows` and the run is what
   knows the lattice; until then no index is written, so adding fields to a
   model that forbids unknown keys breaks nothing. `@tester` (23c-2): a cut
   run's seam lines, read from the piece files, are exactly the lines those
   four fields and `dx`, `dy` give. `@developer` (23c-2): the four fields.
2. **Field types: confirmed, with two changes, both in 23c-1.** `source` and
   `vocabulary` plain strings, `window` as `io.models.IndexWindow` with
   negative `row0`/`col0` allowed, `counts` `{triangles, vertices,
   on_frozen}`: confirmed. Changed: **piece ids are strict integers,
   non-negative** (`["1", "2", "3"]` is refused: an id written as text is
   drift, and lax mode would let it through), and **`file` and `sha256` are
   null together or set together** (a model validator). `@tester` (23c-1)
   adds the two refusals to the schema-drift cases; `@developer` (23c-1)
   the strict type and the validator.
3. **The partition loop: bounded, no change.** The loop runs only while
   `dx · dy · b > B`, and `partition` refuses `B < b`, so inside it
   `dx · dy > 1`. When it grows `nx`, `dx ≥ dy`, so `dx > 1`, so
   `nx < cols` before the step and `nx ≤ cols` after; the same for `ny`
   and `rows`. At `nx = cols`, `ny = rows` the cells are one node and the
   condition fails. So it ends within `cols + rows − 2` steps, each O(1):
   about 90,000 for the basin's window at 1 m. `@developer` puts this
   argument in a comment above the loop; nothing for `@tester` (DC0's
   random draws already reach it).
4. **`partition` refuses an empty window: yes** (`@reviewer`'s suggestion,
   23c-1 round 1). `cols < 1` or `rows < 1` is a `ValueError`, as
   `pieces < 1` and a budget under one node are. The run never passes one
   (the window is the domain's box grown by a cell diagonal, so at least
   2 × 2 nodes), but today `rows == 0` fails as a `ZeroDivisionError` in
   step 2, which says nothing about the cause. `@tester` (23c-1) adds the
   two cases to DC0's refusals; `@developer` the check.

**B17 for Ola** is under "New questions, after the rulings".

## Questions for Ola

Numbered B1-B12, so they cannot be confused with 15's and 15c's Q1-Q17.
The literature pass with web search (2026-10-01, "Prior art") changed no
recommendation; it added a note to B4 and an option (c) to B5. **All twelve
were ruled by Ola on 2026-10-01** ("Ruled by Ola, 2026-10-01", near the
top); they are kept as asked, each marked with its ruling, as are B13 and
B14, asked after them and ruled the same day. B15 and B16 were ruled on
2026-10-02.

**B1. Which cuts first?**
*Ruled (a), and BHO dropped entirely as a geometry source: "I'm not interested in archaic maps".*
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
*Replaced, with B3, by Ola's partition of the bounding box (`--pieces`, a cap, a minimum).*
- **(a) Only when the needed window is wider or taller than `--block-nodes`;
  then along every block line through the domain. Recommended.** Small meshes
  stay bit-identical to today's (K1). It is one rule, not two code paths: an
  uncut run is the one-piece case.
- (b) Always, along the global grid: one rule with no threshold, but any
  mesh crossing a block line gets a seam, small ones included, and every
  existing mesh that does changes.

**B3. The block size.**
*Replaced, with B2 (above).*
- **(a) 2048 nodes by default (61 km at 30 m, 20 km at 10 m), as
  `--block-nodes`, recorded in the file; 0 means never cut. Recommended.**
  About 1 GiB per piece at 1 m in the steepest terrain measured, so eight at
  once fit a 32 GB machine. The 1 m benchmark tile (25.5 M nodes) would be
  cut by default, so `tools/bench.py` passes `--block-nodes 0` from 23c on.
- (b) 4096: a quarter of the seams, up to about 4 GiB per piece at 1 m, so
  fewer at once.
- (c) 1024: four times the seams, a quarter of the memory.

**B4. The seams in the output.**
*Ruled (b), with local seam thinning; designed under "Seam removal and thinning".*
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
*Ruled (a); (b) is the named fallback if 23g's comparison shows a large gap.*
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
*Ruled (a), the stitched file clean; piece files keep their seams.*
- **(a) Pieces and an index always; one stitched file for `--out x.vtk` by
  default, skipped with `--no-stitch`. Recommended.** `--out` keeps meaning
  one file; pieces are the durable result and what resumes a run; at the
  basin's 1 m the stitched file is several GB and can be skipped.
- (b) Pieces and the index only; `rasputin stitch` on request.

**B7. Where the cache lives.**
*Ruled: a new `RASPUTIN_DATA`, the cache at `$RASPUTIN_DATA/cache`, `--cache` overriding, refused if neither.*
- **(a) `--cache DIR` or the `RASPUTIN_CACHE` environment variable, refused
  if neither is set. Recommended.** Explicit, as `--out-crs` is (Q11); one
  more option.
- (b) A default under `../rasputin_data/cache`, next to the repository:
  convenient here, surprising for anyone else.
- (c) The platform's cache directory (`~/Library/Caches/rasputin` on macOS):
  conventional, but 2 GB of DEM blocks hidden there.

**B8. Sources in the first fetch PR.**
*Ruled (a).*
- **(a) ANADEM (OpenTopography's COG) and GLO-30 (AWS). Recommended.** GLO-30
  is what the basin-piece measurement used and is global; about 25 of 23a's
  lines.
- (b) ANADEM only.

**B9. Resumable runs** (skip pieces whose job hash matches).
*Ruled (a).*
- **(a) In 23d. Recommended.** About 35 lines; a basin run interrupted at 1 m
  resumes where it stopped.
- (b) Not now.

**B10. Natural seams' vertices** (for 23e, asked now so it is not a
surprise). *Moot: BHO is dropped (B1).* BHO boundaries bring a vertex every
~107 m, used as given under Ola's input model.
- **(a) Decide when 23e is designed, on measured counts. Recommended.**
- (b) Simplify each shared chain once, between corners, with increment 22's
  area-preserving collapse, so both sides use the same chain.
- (c) Always as given.

**B11. The basin tolerance** (Q15, open). Memory no longer limits it.
*Ruled (a).*
- **(a) The basin run measures 50 down to 1 m and Ola chooses from the
  numbers. Recommended.** Nothing in this design depends on the answer.
- (b) Name it now, and the acceptance stops there.

**B12. The order of work** in "Order of work and PR split".
*Ruled (a), knowingly reversing part of Q14's placement.*
- **(a) As proposed: 15c-1, 23a-1, 23a-2, 15c-2, the edge strip, 23b, 23c,
  23d, the basin run, 23e. Recommended.** The fetch step lands early, so
  every later step runs on ANADEM. This reverses part of Q14's placement
  ruling (the edge strip right after 15c and before 15d): under (a) 23a-1,
  which replaces 15d, and 23a-2 land between 15c-1 and 15c-2, and the edge
  strip after both.
- (b) ROADMAP's current order (15c-1, 15c-2, the edge strip, then 23), with
  15c-2's acceptance on a one-off ANADEM cut. (b) keeps Q14's placement as
  ruled.

### New questions, after the rulings

**B13. The domains of the acceptance runs and of the basin run.**
*Ruled (c): (a) now, (b) in 23e, then the BHO outlines abandoned.* With BHO
dropped as geometry (B1), two domains in this design are still BHO outlines:
the Velhas piece (ottobasin 76949: 15c-2's, 23c's and 23g's acceptance, and
the basin-piece baseline every comparison reads), and the basin itself (BHO
level 2: the basin run, and `partition.py`'s figures).
- **(a) Keep both as measurement domains only, which B1 allows ("BHO at most
  for official codes or validation"). Recommended** for now: every figure
  stays comparable with the basin-piece baseline, and nothing in the code
  reads BHO.
- (b) DEM-derived now: the Velhas catchment by increment 22 on the
  resampled ANADEM grid (a 30.9 M-node window, inside 22's memory cap), and
  the basin's by 22 at basin scale. The basin's is not possible yet: 22
  floods one dense projected window, and the basin's would be the 2.08 G-node
  canvas Ola ruled out; it needs 22 windowed, which is 23e's work. The
  Velhas runs would also lose their baseline.
- (c) (a) now, (b) as part of 23e, with one comparison run of the two
  outlines when it lands.

**B14. The defaults decided here.**
*Ruled: neither; one `--memory-budget`, 16 GB, replaces `N_MIN`, `N_MAX` and
the per-piece cap, and `--pieces` only asks for more ("The memory estimate
and the defaults"). `R = 4` and 23g's thresholds were not addressed and stand
as decided.* As asked: `--pieces 64`, `N_MIN = 2^19` and
`N_MAX = 2^22` nodes (then under "The defaults, by arithmetic"), the bands' `R = 4`
rings, and 23g's gap thresholds (+1 % triangles, +25 % vertices near former
seams, half the worst angle).
- **(a) As decided. Recommended.** The arithmetic is in the record, and 23d
  and 23g measure what it assumed (per-piece cost, balance, deferrals).
  Its cost: every window of 2^20 nodes or more (about 105 km² at 10 m,
  944 km² at 30 m) is cut by default, Bygdin at 10 m and the 1 m benchmark
  tile included, so their default meshes change (cleaned along the seams
  from 23g on); `--pieces 1` gives today's mesh.
- (b) A larger `N_MIN`, `2^21`: domains up to 4.2 M nodes stay whole (Bygdin
  at 10 m most likely among them), but `--pieces` then decides the count only
  over a twofold range of sizes, and the Velhas piece gets 3 × 5 cells instead
  of 60.

**B15. The refusals at half of physical memory.** B14 ruled out a hard
upper limit read from the machine for the partition. Two older refusals do
exactly that: 15a R7 in `plan_mosaic` (`src_python/tin_engine/mosaic.py@5a57793:209`,
a canvas over half of physical memory) and 22's identical one in
`src_python/tin_engine/catchment.py@5a57793:150`. Neither changes a mesh;
each refuses on a small machine what a larger one meshes.
- **(a) Delete both. Recommended.** Matches "a machine too small runs out of
  memory" and "we should not limit huge discetisations [sic] based on less
  performant hardware". 23c deletes them (no lines counted), inverts the
  pinned tests (DC11), and `physical_memory` goes if nothing else reads it.
- (b) Keep both: a clear refusal instead of a crash or swapping, at the cost
  of a machine-dependent limit.
- (c) Delete 15a R7 only (the mosaic canvas, which pieces bound), keep 22's
  until 23e windows the catchment flood.

**Ruled by Ola, 2026-10-02: (a), delete both** ("B15 a").

**B16. What a mesh made from a catalogue source says about its licence.**
Asked 2026-10-02 (23a-2). Today the mesh file's `elevation_source` carries
the catalogue's `credit` only (23a-1). GLO-30's licence, Art. 6(c), asks
that its no-liability sentence accompany any distribution of the data,
modified or not; ANADEM's distributor asks that Laipelt et al. 2024 be
cited. 23a-2 writes all of it to `<cache>/<source>/NOTICE.txt` ("Decided
here (23a-2)" 8), which travels with the cache, not with a mesh.
- **(a) The mesh file also carries `licence_note` and `cite`**, in its
  header beside `elevation_source` (about 5 lines, in 23a-2's `cli.py`).
  Recommended: a mesh handed on keeps the notes the sources ask for.
- (b) `credit` only, as now; the notes stay in `NOTICE.txt` and the docs,
  and whoever distributes a mesh carries them.

**Ruled by Ola, 2026-10-02: (a), the mesh file carries them** ("B16 a"; 23a-2
implements it). Ola also gave the go-ahead to implement 23a-2.

**B17. How `--memory-budget` is written on the command line.** Asked
2026-10-04 (23c's red step). The design says "16 GB, read as 16 GiB"; the
tests take the option as a whole number of bytes, so 16 GiB is typed
`17179869184`.
- **(a) Whole bytes only.** Recommended for 23c, and the default if not
  ruled: no parser, no unit to argue about, and the index records the
  number exactly.
- (b) Also a size with a unit, `16G` or `16GiB`, read as powers of 1024
  (about 10 lines in `cli.py`, and one more red test). The index still
  records bytes. Can follow 23c without changing anything it writes.

**Decided here, which Ola may overrule:** lattice lines on the computation
lattice as artificial cuts; the partition rule's integer details (near-square
cells, the last row and column narrower); `b(T)` from the Velhas piece's
densities, linear between the measured tolerances and rounded up, every
window node costed; 16 GB read as 16 GiB; the runner's admission under the
budget; piece
ids `(j, i, k)`; the seam pass computed by both neighbours rather than once;
seam heights from the seam record; zones by a unit's own edges first (the lowest-numbered at a
corner), else the nearest unit; seam ends on
the outline, holes and features kept; the ETag ignored; coalescing at 64 KiB
gaps and 8 MiB requests; 8 connections and three retries; the manifest
holding identity, not inventory; largest piece first; one process with
threads; the index as JSON; pieces and band records written beside the
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

### Round 4, `3c32d3d..ceb2df2`: CHANGES REQUESTED (`@reviewer`)

Ola's B1-B12 rulings folded in. Partition rule and its arithmetic reproduced
(`partition.py`); prior-art additions match Crossref; LOC under 700 per PR at
+60 %. No rule file changed by the commits whose Bash writes the governance
guard had refused (edits redone with the Edit tool; for Ola to judge). Four
blocking edits: the fence check must run after thinning; stale "not yet
ruled" text in ROADMAP and 15; 15c's DEM-derived fixture cannot exist at the
red step; B14 (a) did not name its cost.

### Round 5, `ceb2df2..1c55a2e`: CHANGES REQUESTED (`@reviewer`)

Round 4 applied; a new zone rule (a triangle with an edge on a unit belongs
to it) found sound and complete. Four gaps: unit ids undefined; the corner
exception path did not rescan; its scope missed kept corners; no test pinned
the rule.

### Round 6, `1c55a2e..0818903`: APPROVED (`@reviewer`)

Whole branch `6518336..0818903`, 19 commits, 8 files, +2,513 / −32,
production LOC 0. 23g 345 (552 at +60 %), 23f + 23g about 645, largest PR 23c
at 656 at +60 %. Citations, ruff, ruff format, mypy and the governance gates
green. CI not yet run: no PR.

### B13 and B14 ruled, round 1, `5a57793..4c182c8`: CHANGES REQUESTED (`@reviewer`)

Production LOC 0; 23c 430 (688 at +60 %). B13 (c) and B14 recorded faithfully; b(T) derives from the basin-piece sweep (17 B per node, 0.805 triangles per node at 1 m); `partition.py` reproduces every row; the 1 m benchmark and Bygdin stay one piece. Blocking: DC0's "at least `--pieces`" could not pass as cells (`--pieces 16` on Velhas gives 3 × 5) and could not fail as P'; deleting 15a R7's refusal was presented as decided under B14 though not asked, and `src_python/tin_engine/catchment.py@96cddea:150` has the same refusal.

### Round 2, `4c182c8..a7c1bfb`: APPROVED (`@reviewer`)

DC0 now pins the cell count, which approximates P' and is not a floor; the refusals at half of physical memory are question B15 (`src_python/tin_engine/mosaic.py@5a57793:209`, `src_python/tin_engine/catchment.py@5a57793:150`), with every dependent place conditional on Ola's answer; units in GiB (7.7 estimated against 3.60 measured). Citations resolve. CI not yet run: no PR.

### 23a-1, round 1, `18316b6..a235716`: CHANGES REQUESTED (`@reviewer`)

LOC 326 net (363 added, 37 removed) against the estimate of 270, which is +21% and inside the +39% worst case. Red came before green, and the green commit touched no test. The local gates were green. Three blocking items:
- The GLO-30 credit was missing the licence's Art. 6(b) notice for adapted data and the Art. 6(c) no-liability sentence.
- The ANADEM credit named no creator, which CC BY 4.0 requires.
- Some prose was made false by the branch: four line citations had drifted, `docs/increments/23-basin-scale.md@0a1e522:1079` still named `_adopt`, and two red-step paragraphs described the tests as still red.

Accepted in round 1: the public `DemTile` constructor instead of `_adopt` (one extra copy per window), NoData taken from the request, and `--cache` ignored for a path `--dem`. Mutation pass: 17 of 19 mutants killed. One survivor was equivalent (thread count). The other showed that the sparse-block test did not pin the refusal message.

### 23a-1, round 2, `a235716..72ec413`: APPROVED (`@reviewer`)

LOC is 336 net (373 added, 37 removed), still under the +39% worst case and far under 700. Both credits now match their sources word for word (GLO-30: the licence's Art. 6(b) notice, Art. 6(c) quoted in `licence_note`; ANADEM: OpenTopography's citation, CC BY 4.0). A missing `header.bin` is a `CacheError`; the sparse-block test fails with the sparse check removed; the four citations resolve to the quoted code; the departure and the `--out-crs` timing are recorded. pytest 3525 passed, 13 skipped in a fresh venv on this worktree's source; ruff, ruff format, mypy, `check_citations` and `check_prohibited_deps` clean. Not blocking: only `credit` reaches the mesh file, so Art. 6(c)'s sentence stays in the catalogue for 23a-2 to carry; `_ascii` writes ANADEM's accented credit as escapes. Remaining: `@perf`'s decode-speed acceptance, then CI after Ola approves the push.

### 23a-2 design, round 1, `99bd723..2120609`: CHANGES REQUESTED (`@reviewer`)

Design only; estimate 373 lines, 518 at +39 % and 597 at +60 %. Correct: block counts recounted from the COG's tie point and step (basin box 8,300 = 100 × 83 against 3,061 for the outline; Velhas 160 = 16 × 10 against 82, matching @perf's fetch); RFC 9110 §13.1.5 says what is claimed; retries, coalescing, lock, `.part`, atomic manifest, `geographic=` for fetch only, `NOTICE.txt` and `cite` consistent; B16 a question. Blocking: F9's "one cell larger is `NotCached`" cannot pass (fetch grows by 6 cells, the mesh window by none); ROADMAP's 285 lines and the 3,061 decoded blocks were stale.

### 23a-2 design, round 2, `2120609..f32a865`: CHANGES REQUESTED (`@reviewer`)

ROADMAP (373), the decode count and F6 fixed. Blocking: F9 asserted which block is missing, but `NotCached` carries only counts. The one-sentence fix (assert `missing` and `needed` against the test's own counts) is in 7f14e7d, checked by the main session against the diff instead of a third round (Ola's two-round cap).

### 23a-2, round 1, `99bd723..226c6dc`: CHANGES REQUESTED (`@reviewer`)

LOC was 598 added and 6 removed, 592 net, against an estimate of 373. That is +60 %, one line over the +60 % column, and under 700. Two packed regions (`FetchRequest(...)` in `cli.py`, `RasterMeta(...)` in `run.py`) are keyword-only constructor calls and still readable. The red step failed its tests (11 failed, 46 errors). The green step touched no test and left only the two tests it reported red. The amendment's reasons held: F4's old box gave a single coalesced range, and Q3 ruled out `CacheWriter`, which the design places in that module. pytest 3489 passed and 16 skipped in a fresh venv; ruff, ruff format, mypy, check_prohibited_deps, check_detria_boundary and check_citations clean. Retries, 206 handling, the lock, `.part` and rename, the atomic manifest, the order of writes and the check for a changed remote all match the design. Only `fetch/http.py` imports networking. `NOTICE.txt` matches the catalogue. All of @developer's departures were sound. Mutation pass: 12 of 14 mutants killed. Blocking:
- B16 (a) was implemented for `.vtk` only. A `.ply` mesh carried neither `licence_note` nor `cite`; a probe confirmed it.
- Prose still called B15 open after Ola ruled it: the Status paragraph, lines ~803 and ~2250, and ROADMAP row 23.

### 23a-2, round 2, `226c6dc..dd76fe8`: APPROVED (`@reviewer`)

- **B16 in `.ply`:** the notes are now written as header comments in both PLY files, and `cite` only when the source has citations. The new PLY cases (cited and uncited) fail against 226c6dc's source and pass at dd76fe8.
- **`coalesce`:** now takes `max(stop, …)`, so a span lying inside another no longer shrinks the range.
- **Progress count:** now updated under a lock.
- **Prose:** the B15 text is correct in the Status paragraph, lines ~803, ~1017, ~1618 and ~2250, the LOC row, DC11 and ROADMAP row 23. The departures are recorded and match the code.
- **Checks:** pytest 3491 passed, 16 skipped in the reviewer's venv (without vtk); ruff, ruff format, mypy and check_citations clean. About 604 production lines, under 700. CI not yet run.

**23b, code review, round 1, 2026-10-03.** Range `595c56a..c7c225f` (red 45045dc, rulings 053e9a2, red 829489c, green 3c464ec, tests f00a7b1, as-built 54578ae, test c7c225f). Verdict: CHANGES REQUESTED. LOC: 235 net (259 added, 24 removed), against about 285. Code correct; K1 holds on the code paths and FE1's digests; M4 and S1 equivalent (M4's reason: row spans exclude vertices), M12/M17 unreachable; the priority-queue greedy matched a naive greedy on 4,000 random cases; N7, N16 and N10 as ruled; bindings fine; TSan entries present. Blocking: (1) `15f-edge-strip.md:562` cites `bindings/core.cpp@44fa7f5:1045` (1057 in the merged tree); (2) the as-built M9 row still says survived (killed after c7c225f); (3) red-step prose in six test headers and the SP1 bullet (2c + 1 check points under N14; N7's on_frozen within r(g) with a strip). Also: merge master (conflict in `tests/cpp/CMakeLists.txt`, keep both; the merged tree passes 954/954) and ROADMAP's 23b row. Not pushed; no CI.

**23b, code review, round 2, 2026-10-03.** Range `c7c225f..e5a6c83` (merge 76de1ea, tests 3882c39, docs e5a6c83); whole PR `b4bcdc3..e5a6c83`. Verdict: APPROVED. LOC: 235 net (259 added, 24 removed), unchanged. All round-1 blockers closed; the merge touches only `tests/cpp/CMakeLists.txt`, both blocks kept; merged HEAD with hardening on: ctest 954/954, pytest frozen/seam/refine_points 108 passed. The single-triangle gap @architect recorded may stay open (rounding-scale, K2 and K4 hold, documented). Suggestions: say the gap also arises on lattice-line seams; name the fix (test proximity against frozen edges incident to the triangle's corners). Not pushed; no CI. Outstanding: @perf's acceptance.

**23b, code review, round 3, 2026-10-04.** Range `c4fb2bf..8da0f2a` (8990504 @perf acceptance: REGRESSION, refine +4.7 % tile / +4.3 % quarter at 1 thread, meshes identical; 8da0f2a the `nodes` loop as one lambda instantiated for frozen and unfrozen, so the unfrozen path drops `on_frozen`); whole PR `b4bcdc3..8da0f2a`. Verdict: APPROVED. LOC: +4 net in this range (16 added, 12 removed, `scan.hpp`), 239 net for the PR. Both paths behave as before (`on_frozen` was `frozen && …`, so the unfrozen instance skips nothing new; the frozen instance skips the same nodes); the comment states the measured cost and claims no gain; the two shifted `scan.hpp` citations in `15f-edge-strip.md` (`:123`, `:80`) re-read and hold; ctest 954/954 in the Release, no-FMA and sanitizer builds (logs, not rebuilt). Before merge: @perf reruns the acceptance at 8da0f2a and appends it, saying the REGRESSION verdict applies to `c4fb2bf` (the timed experiment was an early-return loop, not this lambda); then CI after Ola approves the push.

**23b, code review, round 4, 2026-10-04.** Range `9511563..81c55d6` (8d67c9c @perf: 8da0f2a timed +6.3 % / +6.4 % at 1 thread, REGRESSION; 91c7cb5 @developer: @perf's early-return form, loop body duplicated on purpose; 81c55d6 @perf: ACCEPTED, +1.2 % / +1.8 % at 1 thread, meshes identical, superseding both REGRESSION verdicts); whole PR `b4bcdc3..81c55d6`. Verdict: CHANGES REQUESTED. LOC: +8 net in this range, 247 net for the PR. Both paths behave as before (the `return` leaves the per-segment lambda only); citations `scan.hpp:122`, `:79` hold; the acceptance file's verdict lines are correct. Blocking: the comment in `scan.hpp`'s `nodes` branch says the early return "measured within 1 % of the base", which is true of @perf's experiment, not of 91c7cb5 (+1.2 % / +1.8 %). Comment-only fix, same line count, no re-timing (same binary). Not pushed; no CI.

**23b, code review, round 5, 2026-10-04.** Range `89e8a18..ed5db94` (the `scan.hpp` comment now gives the measured +1.2 % / +1.8 %). Verdict: APPROVED. LOC: 0 this round, 247 net for the PR `b4bcdc3..ed5db94`. Same line count, citations hold; the binary is unchanged, so @perf's ACCEPTED rerun at `91c7cb5` stands. Not pushed; no CI.

**23b, code review, round 6, 2026-10-04.** Range `ed5db94..185081c` (4cd28dc GCC 13 fix in `frozen_oracle.hpp`; e4f07c3 and 24161fe merges of origin/master; e26802e red, 472d91d rulings, f3620b1 red amendment, 185081c green: N18, `refine_strip` takes `frozen_mask`). Verdict: CHANGES REQUESTED. LOC: +1 this round (the `_core.pyi` parameter; the binding's new lines are raw-literal docstring), 248 net for the PR against about 285. N18 is correct: the mask reaches `point_loop`, so N16's refusal shows as `RuntimeError`; f3620b1's "91 added, 0 on the west side" reproduced with mask 0 and with mask 32; pytest frozen, edge strip and refine_points: 102 passed on a rebuilt extension; the frozen and seam property suites and the oracle unit suite pass on clang. N18 alone needs no @perf rerun (bindings and stub only). Blocking: (1) PR #162 conflicts with master #165 (`quality.hpp` comment and fields, `include/terrain/refinement/refine.hpp@ed125121:320-322`), so CI has never run on 4cd28dc; merge, keep both the frozen and the void skips, and @perf reruns the acceptance at the merged head (refine and mesh code resolved; master's 27 and #165 moved the base); (2) citations made false by the merges in this range: `27-node-sampling.md:115,124,148,150,171,172,263` and `25-plain-output.md:115,120` are numbered against master's `refine.hpp` and `scan.hpp`; recompute after the #165 merge; (3) the status line (`23-basin-scale.md:3-7`) and ROADMAP row 23 still say 23a-2 is on its branch (#138 merged) and 23b is not implemented or is 235 lines. Not pushed; no CI.

**23b, code review, round 7, 2026-10-04.** Range `185081c..67ca1ac` (d5aeb82 round 6 recorded; 3403116 merge of origin/master d20126b, #165 void skip and #167; e3a6add citations; 67ca1ac status line and ROADMAP row 23). Verdict: CHANGES REQUESTED. LOC: +1 this round (the merged `quality_skipped` sum), 249 net for the PR against about 285. The merge keeps both skips: `quality.hpp` tests floor, outside, `skipped_void`, the walk, vertex, then `skipped_frozen` (`:130-181`), the order its comment now states; `include/terrain/refinement/refine.hpp@ed125121:320-322` sums all seven. The comment's departure from R4 step 3's list needs no note in `20-start-quality.md`: R4 step 4 already has the walk find the vertex. The round-6 citations in `27-node-sampling.md` and `25-plain-output.md`, and `15f-edge-strip.md:1554,1910`, re-read as quotations at 67ca1ac, hold; status line and ROADMAP row 23 match the tree. `build-23b`, `-nofma` and `-san` current, ctest 974/974 each; touched Python suites 199 passed on the current extension. Blocking: `20-start-quality.md:608-610` says "23b, not yet on master" and calls the conflict "two-line", which the merge of #162 makes false; state the rule instead. Still open: CI's only run on #162 (e4f07c3) is red on GCC 13 (`frozen_oracle.hpp:108`), fixed in 4cd28dc but never run; @perf reruns the acceptance at the merged head against master d20126b (refine and mesh code resolved in the merge). Not pushed; no CI.

**23b, code review, round 8, 2026-10-04.** Range `4541e38..0a3187b` plus `4541e38` (4541e38 round 7 recorded and `20-start-quality.md` states the 23b rule; 0a3187b @perf's acceptance at the merged head). Verdict: CHANGES REQUESTED, closed in the commit that records this round. LOC: 0 this round (no production file in the range), 249 net for the PR against about 285. Round 7's blocker is closed: `20-start-quality.md:608-610` now states the rule, which matches `include/terrain/mesh/quality.hpp@ed125121:64-65,143,182` and `include/terrain/refinement/refine.hpp@ed125121:321-322`. @perf's verdict at `4541e38` against master `d20126b` is ACCEPTED (`docs/benchmarks/2026-10-04/23b-merged-acceptance.md`). The pooled medians, run ranges and changes were recomputed from the eight `raw.tsv` files and match: +0.8 % tile and +0.2 % quarter at 1 thread, -1.7 to +1.7 % over 2 to 20 threads. Each run's commit, module hash, `bench.py` version and AC power state match the doc, and the mesh hashes and quality are identical in all eight runs. `check_citations.py` passes. Blocking: the status line (`23-basin-scale.md:9-10`) and ROADMAP row 23 still said the merged-head acceptance was "owed"; both now say it is recorded. Merge-ready waits on CI: #162's remote head is still `4cd28dc`, which has never had a CI run; not pushed.

**23c-1, code review, round 1, 2026-10-04.** Range `aacf37c..1426245` (red 58fb1f1, 6dd4ffe, 714b0b3; rulings 6ca9a8c; green 3ffae13; rulings 6ce8e57; red b9f154a; green 1426245), stacked on 23b. Verdict: CHANGES REQUESTED. LOC: 154 net (`bindings/core.cpp` 30, `_core.pyi` 3, `decompose.py` 53, `features.py` 3, `io/mesh_index.py` 65) against about 215. No @perf run needed: nothing under `include/` changed, the binding only copies. The loop-bound comment's argument is sound; the binding's GIL and refusals match `refine_points`'s (stricter on negative indices); `project_structure.md` and docstrings match the code. Blocking: (1) "HOW THIS FILE GOES RED" paragraphs in `test_core_mesh_arrays.py:22-23`, `test_decompose_partition.py:27-28`, `test_mesh_index_conformity.py:21-22`; (2) ROADMAP row 23 (still "23c … 430", "23b … in review") and this file's status paragraph ("23b onwards not implemented"). Suggestions: say "integer indices" in `indexed_mesh`'s docstring; a ValueError for `rows == 0` in `partition()`; plain imports for the `dec`/`mi` fixtures. Whichever of 15f-3 and 23c-1 merges second renumbers `tests/python/test_features.py@1426245:583`'s citation. Not pushed; no CI.

**23c-1, code review, round 2, 2026-10-04.** Range `5216cb3..0679950` (c4ec906 test headers, 9633282 status paragraph, ROADMAP row 23 and point 4, 7e5fe0a and e5881a5 the empty-window tests, 0679950 the refusal); whole PR `aacf37c..0679950`. Verdict: APPROVED. LOC: 156 net against about 215. Both blockers closed; the empty-window suggestion taken (point 4, ten cases); citations hold; no @perf run needed. Nit: the status paragraph and ROADMAP say 154 lines built; update to 156 at merge. Not pushed; no CI.

**23c-1, code review, round 3, 2026-10-04.** Range `0679950..b2f6176` (164b83c round 2 recorded; 027e8f9 merge of 23b's head e4f07c3; 9c791d6 merge of origin/master 45acf22, #163; 3b9408f `tests/python/test_features.py@b2f6176:583` cites `project_structure.md:166`; b2f6176 re-cites `bindings/core.cpp@44fa7f5:1174` and `:927-932`, and "156 built"); whole PR `45acf22..b2f6176`. Verdict: CHANGES REQUESTED, on CI alone. LOC: 156 net against about 215. The production and test diffs against master are byte-identical to round 2's; the conflict resolutions keep both sides (core.cpp includes; ROADMAP row 23 is master's text plus the 23c clause; the status line and the 23b rounds 6-8 and 23c-1 rounds 1-2 are all present); the citations hold as quotations; `check_citations.py` passes. Blocking: PR #164's remote head is `027e8f9`, whose CI (run 37177966246) is red on four Linux C++ jobs at `frozen_oracle.hpp:108` (GCC 13), fixed in 4cd28dc, which reaches this branch only through 9c791d6. Push b2f6176 and the push must show every check green; no tree change is needed and no further review round unless the pushed head differs. Remote master is at 72d4608 (#170, #171, docs only), which merges cleanly.

**Row 23 status fix, prose review, round 1, 2026-10-04.** Range `99093af..216513e` (one commit: ROADMAP row 23 and the status line of `23-basin-scale.md`). Verdict: CHANGES REQUESTED. LOC: 0 production lines (prose only; `ROADMAP.md` 1 line changed, `23-basin-scale.md` +7/-5). Checked and true: #162 MERGED at 2026-10-04T09:31:36Z, merge commit `4f56551`; #164 MERGED at 2026-10-04T14:23:36Z, merge commit `e56a2bf`; both are ancestors of HEAD, as are `91c7cb5` and `4541e38`. `check_citations.py --base origin/master` passes. Its 4 at-risk citations were re-read as quotations and hold: `23-basin-scale.md:3267,3271` are dated review records about the status paragraph as it was then, and `27-node-sampling.md:410,412` cite `ROADMAP.md:50`, which is still row 27 (the edit to row 23 adds no lines). Blocking: (1) the commit calls 23c-2 "implemented", but `worktree-23c`'s head `d9a93b6` is a red step that only adds tests (`test_cli_mesh_pieces_settled.py`, +49), and its message says "Red today", so the branch currently has failing tests and the round-3 approval `45e822d` no longer covers its head; (2) the 23c-2 sentence uses `worktree-23c` as its source but leaves out what that branch's own design already says (`b5a73ec`: 23c split in three, 23c-3 about 185 lines, nine PRs, 23c-2 built at 478 net lines), while the same edited sentences on master still say "split in two" and "eight PRs"; (3) `45e822d` and `d9a93b6` are on no remote branch, so a reader of master cannot resolve them. They go stale with the next commit on `worktree-23c`, which breaks the rule "write the rule and the command that resolves it, not the resolved value". Fix: on master, keep the 23b and 23c-1 merges, and for 23c-2 say only that it is in progress on `worktree-23c`, not pushed, and that its state is in that branch's copy of this file (no hashes, no review rounds). Not pushed; no CI run exists for this branch (`gh run list --branch worktree-row23-fix` is empty and there is no remote head).

**Row 23 status fix, prose review, round 2, 2026-10-04.** Range `216513e..d61de61` (d61de61: the 23c-2 sentence in ROADMAP row 23 and in the status line of `23-basin-scale.md`, and round 1 recorded); whole branch `99093af..d61de61`. Verdict: APPROVED. LOC: 0 production lines (prose only; `ROADMAP.md` 1 line changed, `23-basin-scale.md` +8/-4 including the round-1 record). All three blockers from round 1 are closed. (1) 23c-2 is no longer called "implemented". It is now "in progress on `worktree-23c`, not pushed", which matches that branch's head `d9a93b6`, a red step that only adds tests. (2) For 23c-2's state, both texts now point to that branch's copy of this file and no longer restate it. "Split in two" and "eight PRs" stay as master's design until `worktree-23c` (where `b5a73ec` says split in three and nine PRs) merges. (3) The hashes found on no remote branch (`45e822d`, `d9a93b6`) and the review-round count are gone from both texts. Checked and true: #162 MERGED 2026-10-04T09:31:36Z, merge commit `4f56551`; #164 MERGED 2026-10-04T14:23:36Z, merge commit `e56a2bf`; `4f56551`, `e56a2bf`, `91c7cb5`, `4541e38`, `3c464ec`, `45045dc` and `185081c` are ancestors of HEAD; `docs/benchmarks/2026-10-04/23b-merged-acceptance.md` exists; `worktree-23c` is not on origin (`git ls-remote`); the branch contains origin/master `99093af`. The facts in the round-1 record hold (`b5a73ec` is on `worktree-23c` and says split in three, nine PRs and 478 net lines; `d9a93b6` adds only `test_cli_mesh_pieces_settled.py`, +49, and its message says "Red today"). `check_citations.py --base origin/master` passes (exit 0) and lists 6 at-risk citations. Each was re-read as a quotation and holds: `23-basin-scale.md:3267` and `:3271` are still dated 23b review records (rounds 6 and 8); `:3279`, the round-1 record, cites them as such and cites `ROADMAP.md:50` as row 27, which it still is (row 23 is line 51); `27-node-sampling.md:410,412` cite `ROADMAP.md:50` as the sampler row, which is now row 27. Not pushed; no CI.
