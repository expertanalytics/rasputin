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
   20b's feet (`refine.hpp@98da562:339`); the quality pass's node on an edge
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
  features and the other lines); a block corner is the crossing of two axis
  lines, which is exact and lies on the noder's 1 mm grid when `h` is whole
  metres.
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
   the points lie on the edge to rounding, and the noder keeps every vertex at
   least half a snap cell (0.5 mm) from an edge it is not on (05b's hot-pixel
   guarantee), so the fan is valid. If it is not, refine's own check refuses
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
     sixth path added later cannot do it silently.
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
  tifffile parses every page without reading past it. ANADEM needs at most
  8 MiB (measured). The prefix is cached as `header.bin`.
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
testable: an import test that no module on `rasputin mesh`'s path imports
`tin_engine.fetch` or `urllib.request`, and a mesh run with `socket.socket`
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
- **15d's window decoding survives** (15 R3 [15d] and its T-window test) and
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

## Order of work and PR split

(pending)

## Tests @tester can write red

(pending)

## @perf acceptance

(pending)

## Questions for Ola

(pending)
