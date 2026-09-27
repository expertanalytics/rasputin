# Increment 21 — parallel refine: quick wins, then the serial phase

Status: **design only, ruled by Ola 2026-09-27** (section 8, "Ola's rulings"). Written by `@architect`
on 2026-09-27, on branch `serial-profile`, per `docs/increments/README.md`
step 1. No code and no tests exist for it. This file proposes a sequence of
small increments (21a to 21d, section 7); each gets its own Rulings and Tests
sections when Ola has answered section 8 and the measurements in section 6
are in.

**Closes, when all of it lands.** The ROADMAP row "The serial phase": refine
tops out at about 2.0x from 1 to 20 threads. **Not closed:** the scan's own
ceiling near 5x beyond what chunking recovers, E-core scheduling, and very
large areas that do not fit in memory (the virtual mosaic, increment 15
refreshed, `18-row-span-scan.md` R6).

**Non-negotiable, whatever Ola rules on determinism.** Every option below
keeps:

- the sup-norm tolerance as a property of the delivered mesh: every triangle
  is scanned in its final form and is within tolerance
  (`14b-delaunay-insertion.md` R4, stop condition);
- constrained Delaunay output in the lattice frame, constraint edges never
  flipped (14b R3, R5), with a Delaunay oracle on the output (14b R10);
- termination for every tolerance, 0 included (14 R3, 14b R4);
- the I/O boundary and the no-CGAL/no-GDAL mandates (`CLAUDE.md` §2).

## Ola's rulings this design is built on

Recorded verbatim, not paraphrased.

- 2026-09-27, before the profile:
  "Some of the serial parts could be parallellized by multicoloring/dd techniques, but let's wait for the analysis before we get ahead of ourselves."
- 2026-09-27, after the profile, on the determinism contract: **"I think that
  when we are encountering geometries where this is important, ie very large
  areas, we should not rely on bit-identical outputs."**

**What the second ruling settles and what it does not.** It clearly allows
output that differs from today's mesh, which is fixed by triangle-index
order (`14-adaptive-refinement.md` R5, T6). It does **not** say whether the
output may then vary with the thread count, or from run to run. It also
leaves open whether "when we are encountering geometries where this is
important" means a switch (today's path for small inputs, a parallel path
for large ones) or one path for everything. Section 4 lays out the levels;
questions Q1 and Q2 in section 8 ask Ola to choose.

**Ruled by Ola, 2026-09-27: "yes to all"** to the recommendations in
section 8:

- **Q1: L1.** Deterministic and independent of the thread count, in a new
  order; T6 stays, increment 18's golden digests are re-recorded once.
- **Q2: one path** for every input; every mesh changes once when 21d lands.
- **Q3: at most 2 % more triangles** at 1 m on the quarter circle, worst angle
  and max degree no worse, measured in 21c before choosing between C and A1.
- **Q4: 21a and 21b go ahead now**; 21b waits for the tie-classification
  measurement.
- **Q5: domain decomposition is deferred** to the large-area (mosaic) work.

## 1. Prior art: legacy and literature

### Legacy

```sh
$ grep -rliE 'std::thread|pthread|parallel|concurren|openmp|pragma omp|tbb|subdomain|independent.?set|multiprocessing|ThreadPool' legacy/
legacy/rasputin/reader.py
```

The only hit is a false positive: `legacy/rasputin/reader.py:63-64` are the
GeoTIFF key names `ProjStdParallel1GeoKey` and `ProjStdParallel2GeoKey`.
The legacy tree called CGAL's `Delaunay_triangulation_2` and
`Constrained_Delaunay_triangulation_2` whole (14b, "Prior art in `legacy/`")
and ran nothing in parallel. **Nothing is carried across.**
`@migration-expert` is not needed.

### Literature: what was searched

**No literature database or web search was available in this session.** The
list below is from the architect's recall and from the citations already in
this repository (`docs/research/data-dependent-triangulation.md`,
`14b-delaunay-insertion.md` R7, `parallel_refinement.md`). Titles, venues and
years are given as recalled; **each must be checked against the paper before
it is used in a publication**, as `docs/research/data-dependent-triangulation.md`
already requires for its own list. Where the content of a paper is not
recalled with confidence, this file says so rather than quoting it.

**No novelty is claimed in this file.** Section 1's last subsection lists the
searches to run before any claim.

### The method we build on: greedy insertion on a height field

- **De Floriani, Falcidieno and Pienovi, "Delaunay-based representation of
  surfaces defined over arbitrarily shaped domains", Computer Vision,
  Graphics and Image Processing 32 (1985).** Greedy insertion of the
  worst-fitting data point into a Delaunay triangulation until an error
  bound holds. Recalled as the origin of error-driven Delaunay refinement of
  terrain; to be checked.
- **Garland and Heckbert, "Fast polygonal approximation of terrains and
  height fields", CMU-CS-95-181 (1995).** Greedy insertion: one point at a
  time, the globally worst, into a Delaunay triangulation with Lawson flips,
  with each triangle's candidate (its worst point) recomputed only for
  triangles that changed, and a heap over triangles. That is what increments
  14, 14b and 18 rebuilt
  (`docs/retrospectives/2026-09-27-increments-14-to-20b.md`).
- **How we differ from Garland and Heckbert, today.** Our round inserts one
  point *per unconverged triangle*, not the single global worst, in
  triangle-index order, and skips a triangle whose slot an earlier insertion
  of the same round rewrote (14 R5, 14b R2). That is a batched greedy, not
  strict greedy. Whether the report discusses inserting several points per
  pass has to be read in the report and is not quoted here. What the batching
  costs in triangles against strict greedy has not been measured. It matters
  here because every option in section 5 changes the batch.

### Parallel Delaunay insertion and refinement with determinism

- **Blelloch, Fineman, Gibbons and Shun, "Internally deterministic parallel
  algorithms can be fast", PPoPP 2012.** Deterministic reservations: each
  round, every pending item reserves the locations it will touch with a
  priority-min write; an item holding all its reservations commits; losers
  retry. The result depends only on the priorities, not on the thread count.
  Delaunay triangulation and Delaunay refinement are among its benchmarks (the
  PBBS suite). Already cited as the upgrade path in 14b R7.
  **Difference from us:** their refinement inserts circumcentres to fix bad
  angles, with no data; ours inserts DEM nodes chosen by a scan, with a
  tolerance guarantee checked by rescans, and on constraints.
- **Blelloch, Gu, Shun and Sun, "Parallelism in randomized incremental
  algorithms", SPAA 2016; J. ACM 67(5), 2020.** Shows that incremental
  Delaunay triangulation in a *random* insertion order has shallow dependence
  depth (recalled as O(log n) with high probability), and that the parallel
  version gives exactly the sequential result for that order. **What it tells
  us:** priority order matters. A random-like priority (for example a hash of
  the inserted node) gives short dependence chains; a spatially correlated
  order such as our triangle indices may not. Section 6, question 2, asks for
  that depth to be measured.
- **Nguyen, Lenharth and Pingali, "Deterministic Galois: on-demand, portable
  and parameterless", ASPLOS 2014.** Deterministic scheduling of irregular
  programs, Delaunay mesh refinement among them, with output independent of
  the thread count. Same level as Blelloch et al.; a runtime rather than a
  per-algorithm design.
- **Kulkarni, Pingali, Walter, Ramanarayanan, Bala and Chew, "Optimistic
  parallelism requires abstractions", PLDI 2007.** Galois. Delaunay mesh
  refinement is the motivating example: speculative cavity insertion with
  conflict detection and rollback. Nondeterministic.

### Parallel Delaunay refinement by geometric scheduling and decomposition

- **Chernikov and Chrisochoides, "Practical and efficient point insertion
  scheduling method for parallel guaranteed quality Delaunay refinement",
  ICS 2004.** Points whose insertion regions are far enough apart are
  independent; a uniform grid or quadtree of cells is coloured so that cells
  of one colour can be refined concurrently. This is the multicolouring Ola
  named. **Difference:** their independence distance comes from the quality
  bound of Ruppert/Chew refinement (circumradius bounds). Our insertions have
  no such bound: a footprint's extent is data-dependent, so independence has
  to be checked, not derived (section 5, option B).
- **Chernikov and Chrisochoides, "Generalized Delaunay mesh refinement: from
  scalar to parallel", IMR 2006**, and **"Three-dimensional Delaunay
  refinement for multi-core processors", ICS 2008.** Later work in the same
  line. Recalled, not reread.
- **Spielman, Teng and Üngör, "Parallel Delaunay refinement: algorithms and
  analyses", IJCGA 17(1), 2007 (IMR 2002).** Parallel Ruppert-style
  refinement in rounds of independent insertions, with a polylogarithmic
  round bound. Same difference as above: the independence argument rests on
  the quality criterion.
- **Linardakis and Chrisochoides, "Delaunay decoupling method for parallel
  guaranteed quality planar mesh refinement", SIAM J. Sci. Comput. 27(4),
  2006.** Domain decomposition: separators are pre-refined so that
  subdomains refine with no communication, and the union is still Delaunay.
  **Difference:** decoupling relies on knowing ahead how fine the separator
  must be, from the quality criterion and a sizing function. Our refinement
  is driven by DEM error, which is only known by scanning.
- **Galtier and George, "Prepartitioning as a way to mesh subdomains in
  parallel", IMR 1996.** Partition first, then mesh each part; interfaces
  fixed in advance.
- **Chrisochoides, "A survey of parallel mesh generation methods", in
  Numerical Solution of PDEs on Parallel Computers, LNCSE 51, Springer
  (2006).** The survey to read first; recalled as classifying methods into
  tightly coupled (concurrent insertion with synchronisation), partially
  coupled and decoupled (DD).
- **Antonopoulos, Blagojevic, Chernikov, Chrisochoides and Nikolopoulos,
  "Multigrain parallel Delaunay mesh generation: challenges and opportunities
  for multithreaded architectures", ICS 2005.** Fine-grained concurrent point
  insertion within a subdomain. Recalled as finding that the finest grain
  pays only with cheap synchronisation; to be checked, because it bears
  directly on our density (section 2).

### Parallel Lawson flipping

- **Qi, Cao and Tan, "Computing 2D constrained Delaunay triangulation using
  the GPU", I3D 2012 (IEEE TVCG 19(5), 2013).** Recalled as: insert points in
  parallel, at most one per triangle per round, then restore the (constrained)
  Delaunay property by rounds of parallel flips on edges that do not share a
  triangle. **This is structurally close to our round**, which already takes
  one point per triangle. The difference: we legalise after each insertion,
  serially, and skip triangles an earlier insertion rewrote; they insert the
  whole batch, then flip. Section 5, option C.
- **Navarro, Hitschfeld-Kahler and Mateu, "A parallel GPU-based algorithm for
  Delaunay edge-flips", EuroCG 2011.** Parallel flipping of a whole
  triangulation to Delaunay with conflict resolution between edges that share
  a triangle.
- **Lawson, "Software for C1 surface interpolation", in Mathematical
  Software III (1977).** The flip algorithm and its convergence: any sequence
  of flips of locally non-Delaunay edges terminates at the Delaunay
  triangulation, because each flip lowers the lifted surface. That argument
  does not depend on flip order, which is what makes parallel flipping
  terminate (14b R4, part 1).

### Parallel Delaunay with locks (nondeterministic)

- **Kohout, Kolingerová and Žára, "Parallel Delaunay triangulation in E2 and
  E3 for computers with shared memory", Parallel Computing 31 (2005).**
- **Batista, Millman, Pion and Singler, "Parallel geometric algorithms for
  multi-core computers", CGTA 43(8), 2010.** The CGAL parallel triangulation,
  with per-vertex locks. CGAL itself is prohibited here; the method is only
  a reference point.

### Parallel greedy insertion on terrain

**Not found from recall.** No paper on parallel greedy (error-driven)
insertion for height fields is known to the architect, and none is cited in
this repository. This is exactly the gap a novelty claim would sit in, so it
is not filled by guessing.

### Checked by web search, 2026-09-27 (main session)

The design was written from memory. Afterwards the main session checked the
three papers the options rest on, and ran one terrain query:

- **Verified:** Blelloch, Fineman, Gibbons and Shun, "Internally deterministic
  parallel algorithms can be fast", PPoPP 2012
  (doi:10.1145/2145816.2145840). Chernikov and Chrisochoides, "Practical and
  efficient point insertion scheduling method for parallel guaranteed quality
  Delaunay refinement", ICS 2004, pp. 48-57 (doi:10.1145/1006209.1006217); its
  independence condition compares point distance with an upper bound on
  triangle circumradius, as option B assumes. Qi, Cao and Tan, "Computing 2D
  constrained Delaunay triangulation using the GPU", I3D 2012, extended in IEEE
  TVCG 19(5):736-748, 2013; it flips all flippable pairs in parallel, as
  option C assumes.
- **Found, not in the list above:** "3D Simplification Methods and Large Scale
  Terrain Tiling", Remote Sensing 12(3):437, 2020 (mdpi.com/2072-4292/12/3/437).
  It adapts greedy insertion, among other methods, to work tile by tile, in
  parallel, keeping tile-border vertices shared between neighbours. That is
  prior art for option D (domain decomposition by tiles) on terrain. Also "A
  fast digital terrain simplification algorithm with a partitioning method",
  IEEE, 2000 (ieeexplore.ieee.org/document/843506), not read.
- The other citations above are still unverified, and the searches below have
  not been run in full.

### Searches to run before any novelty claim

A claim is not made until these are run in a real database (Google Scholar,
ACM DL, IEEE Xplore, Scopus) and the results recorded here:

- "parallel greedy insertion terrain", "parallel TIN generation DEM",
  "parallel height field approximation", "GPU terrain triangulation error";
- "deterministic parallel Delaunay refinement", "thread-count independent
  mesh generation";
- "parallel Delaunay refinement terrain", "sup-norm" or "L-infinity" with
  "TIN" and "parallel";
- forward citations of Garland and Heckbert 1995 and of Blelloch et al. 2012
  that mention terrain or height fields;
- venues: IMR, SoCG, SPAA, PPoPP, SIGSPATIAL, IJGIS, Computers & Geosciences.

The candidate that might be ours, per the retrospective: a deterministic,
thread-count-independent parallel greedy refinement to an exact sup-norm
tolerance on a DEM lattice, with constraints, whose output is constrained
Delaunay. Blelloch et al. 2012 already makes deterministic parallel Delaunay
refinement not new in itself; only the combination could be, and that is
unchecked.

## 2. What the profile says

The evidence is `docs/benchmarks/2026-09-27/serial-profile/README.md` (commit
`ec1f7b0`, measured at `d6d6beb`, battery, the 1 m benchmark on the quarter
circle). It is not restated here. The design depends on these figures, and
on nothing else from it:

| figure | value | used for |
|---|---|---|
| serial part of 1-thread refine | 34 % (0.159 of 0.462 s) | the Amdahl ceiling, 2.9x |
| split phase at 1 / 8 threads | 0.127 / 0.137 s | what options A-C must parallelise |
| exact incircle fallback | 5.0 % of refine; 9.1 % of incircle tests, all Cocircular | quick win 2 |
| `active` rebuild (collect, sort, unique) | 3.6 % | quick win 3 |
| scan imbalance at 8 threads | 16.5 of 61.8 ms | quick win 1 |
| thread start per round | about 90 µs, 3.8 ms over 41 rounds | why options A-C need a thread team, not a spawn per step |
| insertions per round | 77 % in rounds 9-19, 9,000-18,700 each | parallel slack per round |
| median distance to the nearest same-round insertion, rounds 11-29 | about 3 nodes | independence at the finest grain is rare |
| footprint (slots read or written) | mean 11.0, p99 18, max 36 | the size of a reservation |
| write set | mean 5.5, p99 9 | slots appended and rewritten |
| footprint conflict degree, big rounds 9-19 | mean 3.2-6.4 (7.7 in round 5) | expected winners per sub-round |
| marks deferred because their slot was rewritten | 51 % | the batch semantics every option must match or change on purpose |

Two caveats carried from the evidence. The footprints were measured in the
serial order, on the mesh as earlier insertions left it, so they are not
exactly the footprints a parallel schedule would see. And the profile was
taken at one thread; the split phase is about 8 % slower after a
multi-threaded scan, cause not measured.

## 3. Quick wins that need no new algorithm

All four keep the output **bit-identical to today** at every thread count, so
none needs Ola's determinism ruling. Expected gains are arithmetic on the
profile's figures, not measurements.

### QW1. Work-balanced scan scheduling

- **What.** `for_each_chunk` splits `active` into equal counts
  (`include/terrain/parallel_util/chunks.hpp`). Replace that, for the scan,
  with dynamic scheduling: workers take blocks of `active` from a shared
  atomic counter until it runs out. Block size about `n / (16 * threads)`,
  at least 1. Rounds with few active triangles (the profile's rounds 30-41
  scan fewer than 1,200 each, mostly spawn cost) run inline below a
  threshold to be set from the sweep.
- **Why dynamic and not cost-weighted static chunks.** A static weight (node
  count from each triangle's integer area) needs a serial or prefix-sum pass
  per round, and it cannot see the measured causes that are not work: the
  per-worker slowdown, the E cores, and the late start above 10 threads.
  Dynamic scheduling absorbs all three.
- **Determinism.** Unchanged. The scan is pure and writes one result slot per
  triangle (14 R7, 18 R4), so which thread scans which triangle cannot change
  a result. The atomic is a work counter only. This amends 14 R7's "no
  atomics" in wording, not in substance: R7's point is no shared *result*
  writes, which still holds. TSan sees a correctly synchronised atomic.
- **Exceptions.** `for_each_chunk` promises that the exception of the
  lowest-index chunk that threw is rethrown, "fixed by the chunking, not by
  thread timing" (`include/terrain/parallel_util/chunks.hpp:10-13`, tested at
  `tests/cpp/unit/test_refinement_chunks.cpp:79`). Dynamic scheduling keeps
  that contract: each block records its exception by block index, and the
  lowest block index that threw is rethrown after the join. Every block is
  still run, as today, so which blocks throw does not depend on timing.
- **Expected gain.** Up to the 16.5 ms imbalance at 8 threads, less one
  block's work per thread. Scan roughly 62 -> 47 ms; refine at 8 threads
  roughly 0.229 -> 0.215 s (-6 %). Nothing at 1 thread.
- **Its own small increment:** yes, with QW3 (21a).

### QW2. An integer incircle for quads whose four corners are DEM nodes

- **What the profile shows.** 146,962 incircle tests took the exact path, and
  all 146,962 returned Cocircular. The filter cannot certify a determinant
  that is exactly zero, so every exact tie goes to the adaptive path. That
  costs 5.0 % of refine and 38 % of flip-test time.
- **The inference to check first.** That these ties are cocircular lattice
  quads (for example a grid rectangle's four corners) is inferred, not
  measured (profile README, "What is measured and what is inferred").
  **Measurement, before 21b is designed in detail** (`@perf`, the existing
  instrumentation patch, `scripts/instrument.patch`): at every exact-path
  incircle call in `must_flip`, record
  1. whether all four vertices are nodes (`MeshVertex::is_node`);
  2. the int64 lattice determinant of the four (below), which should be 0;
  3. the shape: axis-aligned rectangle (two distinct rows and two distinct
     columns), or another cocircular lattice quad.

  Also count the filtered-path calls whose four corners are all nodes, since
  the integer path would take those too. Run on the quarter circle at 1 m and
  on the full tile without a domain. The quarter circle has off-node arc
  vertices (`16-domain-polygon.md` R2), and 20b's feet are off-node, so the
  all-node share is not 100 % by construction. The gain scales with it.
- **The design.** A pure function in `mesh/`, beside `lawson.hpp`:
  `lattice_incircle(a, b, c, d) -> std::optional<Incircle>`, called at the
  top of `must_flip`. It answers only when
  - all four vertices are nodes;
  - `dx == dy`, so the frame is the lattice times one positive constant, and
    the sign of the scaled determinant is the sign of the lattice one. The
    integer determinant is taken on `(col, -row)`, as `MeshVertex::frame()`
    orients the frame; on `(col, row)` the sign inverts;
  - every coordinate difference from `d` is at most 2^14 nodes. The
    determinant is then at most 12 · 2^56 < 2^63, and exact in `int64`;
  - **the frame is exact**: `col * dx` and `row * dy` are exactly
    representable for every node of the grid. Sufficient: the significant bits
    of `dx` plus `bit_width(max(rows, cols) - 1)` are at most 53. For the
    benchmark's `dx = 10` (3 significant bits, 13 for 5051) that holds with
    room to spare. It is decided once per refine call.

  Otherwise it returns `nullopt` and today's `FilteredKernel<DetriaExact>`
  path runs. The mesh triangle `a, b, c` is counter-clockwise in integers
  (`LatticeMesh` invariant), and under an exact frame it is so in the frame
  too, so the integer path also skips `must_flip`'s frame `orient2d`.
- **Determinism.** Under the four conditions, `DetriaExact` on the frame
  doubles computes the sign of the same exact number, so every flip decision
  is the same and the output is bit-identical. The exact-frame condition is
  what makes that true; without it (say `dx = 0.1`) the integer answer is the
  *true* lattice answer and the rounded frame's answer can differ on a tie.
  That would probably still terminate (unproven: an integer-strict Inside
  seems to keep a margin far above the rounding), but it would change output,
  so it is excluded and nothing rests on the argument.
- **Expected gain.** Most of the 5.0 %, plus part of the 2.7 % filtered cost
  for all-node quads, less the integer determinant's own cost. Roughly
  0.018 s at every thread count: -4 % at 1 thread, -8 % at 8.
- **Its own small increment:** yes (21b). It adds a predicate path, so its
  suite is invariant-critical (section 7).

### QW3. Rebuild `active` by merging, not sorting

- **What.** `refine.hpp` rebuilds `active` each round by collecting touched
  slots, appending `skipped`, then `std::sort` and `std::unique`. Both inputs
  are already sorted: touched slots are collected by an ascending loop over
  slots, and `skipped` is filled while iterating the sorted `active`. So a
  linear `std::merge` plus `std::unique` gives the same vector. `unique` is
  still needed: a skipped triangle can be touched later in the same round.
- **Determinism.** Same vector, so bit-identical.
- **Expected gain.** Most of 3.6 %: about 0.014 s at every thread count.
- **With QW1 in 21a.**

### QW4 (conditional). Do not rescan a skipped triangle whose slot is still unwritten

- **What.** A mark skipped because its edge neighbour was touched is itself
  unchanged, and if nothing writes its slot later in the round, its stored
  result is still exact (14b R2's invariant). Today it is rescanned anyway.
  It can go to the next round's marks without a scan.
- **Determinism.** The scan is pure, so a rescan of an unwritten slot
  returns the same result. Bit-identical.
- **Gain: unknown.** It saves the node visits of those triangles only (edge
  deferrals are 7 % of mark occurrences). **Measurement:** in the
  instrumented build, split each round's scanned nodes into slots written in
  the previous round and slots skipped-and-unwritten. Build it only if the
  second share is worth a few percent of the scan.

### What the quick wins add up to

Arithmetic, not measurement: 1-thread refine 0.462 -> about 0.43 s; 8
threads 0.229 -> about 0.18 s, which is about 2.4x over the new 1-thread
time. The serial part is then still about three quarters of the 8-thread
time, so the ceiling moves only a little. The quick wins are worth doing
first because they are cheap and bit-identical, and because every option in
section 5 is measured against the mesh they leave.

A smaller item folds into 21a: `legalise_around` allocates its stack per call
(0.7 % of refine); a caller-owned buffer removes it.

### Pinned by the red suite (21a)

Section 3 left these names and signatures open; `@tester` chose them in the
red step, and the suites named below hold them.

- **QW1** (`include/terrain/parallel_util/chunks.hpp`,
  `tests/cpp/unit/test_refinement_chunks_dynamic.cpp`):
  `struct BlockSchedule { std::size_t block = 0; std::size_t inline_below; }`,
  `constexpr std::size_t default_block(std::size_t n, unsigned threads) noexcept`
  returning `max(1, n / (16 * threads))`, and
  `template <class Fn> void for_each_block(std::size_t n, unsigned threads, BlockSchedule, Fn&& fn)`.
  `block == 0` means `default_block(n, threads)` after `threads == 0` is
  resolved as `for_each_chunk` resolves it. Block `k` is
  `[k*b, min(n, (k+1)*b))` and `fn` is called exactly once per block, so the
  set of calls is the same partition for every thread count. With one thread
  or `n < inline_below`, every block runs on the calling thread in ascending
  order. Otherwise at most `min(threads, blocks)` threads call `fn`, and blocks
  run concurrently. Every block runs even when some throw; the exception of the
  lowest block index that threw is rethrown after the join. The default of
  `inline_below` is 21a's to set from the sweep and is not pinned.
  `BlockSchedule{}` is the scan's schedule. `for_each_chunk` and its suite are
  unchanged.
- **QW3** (`include/terrain/refinement/refine.hpp`,
  `tests/cpp/unit/test_refinement_active.cpp`):
  `void detail::rebuild_active(std::span<const char> touched, std::span<const std::uint32_t> skipped, std::vector<std::uint32_t>& active)`.
  It replaces `active`'s contents with what today's collect, sort and unique
  give. Precondition: `skipped` is ascending and below `touched.size()`.
- **The stack** (`include/terrain/mesh/lawson.hpp`,
  `tests/cpp/unit/test_mesh_lawson_stack.cpp`):
  `using FlipStack = std::vector<std::uint32_t>;` and a `legalise_around`
  overload taking `FlipStack& stack` before `on_write`. It must produce the
  same flips, `on_write` sequence and mesh as today's algorithm, which the
  suite copies as its oracle. Allocation counts are not pinned.

At the red commit the three suites were registered only once their header
named `for_each_block`, `rebuild_active` or `FlipStack` (the 18 and 20b
precedent); `065a2bd` removed those guards after green.

### 21a: inline_below

Set to **256** by `@developer` in the green step, from a sweep that is a
developer's quick look, not `@perf`'s acceptance run. Setup: the 1 m
benchmark's tile (`tests/fixtures/dem_archive/7908_3_10m_z33.tif`, tolerance
1, no domain), M1 Max on AC. A throwaway driver wrapped `cli.refine` the way
`tools/bench.py`'s child does and called it 9 times per thread count, taking
the median of `scan_seconds` and of the wall time of refine. `_core` was
rebuilt for each value. Median scan in ms at 8 threads (10 threads within
3 ms of it):

| inline_below | 0 | 64 | 256 | 512 | 1024 | 2048 | 8192 | 32768 |
|---|---|---|---|---|---|---|---|---|
| scan, 8 threads | 58.4-58.8 | 57.2 | 57.3-57.5 | 57.5-57.6 | 58.8 | 65.9 | 83.3 | 119.7 |

From 64 to 1024 the values differ by less than the run-to-run noise; 64-512
are about 1 ms under 0, and 1024 is level with it. From 2048 up, the inline small rounds cost more
than their thread starts save. 256 is in the middle of the flat range. Measured
back to back against the red commit (`7b54a4e`, `for_each_chunk`): scan at
8 threads 65.0 -> 57.7 ms and refine 243 -> 219 ms (-10 %). At 1 thread there
is no difference beyond noise (refine 515 vs 519 ms), as QW1 predicts. The
mesh sha256 that `bench.py` prints for the tile and the quarter circle was the same
before and after.

## 4. Determinism levels

Today's contract (14 R5, 14b R1, tested by 14's T6 and 18's T3 golden
digests) is the strongest level. Ola's ruling allows leaving it; how far is
the question.

| level | what holds | what it buys | what it costs |
|---|---|---|---|
| **L0** bit-identical to today | same mesh as the serial triangle-index order, any thread count | every golden digest and T6 stay as they are | the parallel schedule must reproduce today's greedy order, including which marks are skipped and how appended slots and vertices are numbered (section 5, option A0). Expensive, and its parallelism depends on a dependence depth nobody has measured |
| **L1** deterministic, thread-count-independent, new order | same mesh for 1, 2, 7 or 20 threads, and on every run; different from today's | T6 keeps its meaning. A bug seen on one machine replays on another. A benchmark mesh can be hashed and cited, which the publication option needs. Golden digests re-recorded once | the order must come from the data (priorities, a fixed block grid), never from the thread count. Some bookkeeping to number new slots deterministically |
| **L2** deterministic per thread count | same mesh for the same `threads`; a different mesh for a different one | nothing over L1 for the options below; it would matter only for a DD with one subdomain per thread | the output depends on the machine. Reproducing a user's mesh needs their thread count, which the CLI does not expose (14, "Not in scope") and the output would then have to record. T6 weakens to "same threads, same mesh" |
| **L3** nondeterministic | tolerance and CDT only; the mesh varies run to run | the lock-based designs (Kohout et al.; Galois) become available | failures do not replay. Tests can only be property tests. Cocircular ties (9 % of incircle tests are exact ties) then resolve differently each run, so two runs differ even with no bug. TSan becomes the only race evidence |

**Recommendation: L1.** Every option recommended in section 5 reaches L1
at no extra cost over L2, because the schedule is keyed to the data (a hash
of the inserted node, or a fixed block grid in DEM nodes), not to the thread
count. L2 buys nothing those options need. L3 gives up replayable failures for
a speed-up no measurement says we need. L0 is possible in principle but costs
the most and promises the least (option A0).

**One path, not a size switch.** One reading of the ruling is "keep today's
path for small inputs, a parallel path for very large ones". That keeps two
production refinement paths, which 18's C4 rejected for the scan on the same
grounds: two things to test, and a threshold to tune. Under L1 there is no
need: the 1-thread run of the parallel algorithm is itself deterministic. The
cost is that **every** mesh changes once, small ones included. Question Q2.

## 5. Options for the serial phase

Each option is judged against the measured density: insertions about 3 nodes
apart in rounds 11-29, each conflicting with 3.2-6.4 others in the big
rounds 9-19 (7.7 in round 5, 2,063 insertions), footprints of about 11 slots (p99 18), 9,000-18,700 insertions in a
big round, 51 % of marks deferred by the serial loop.

The ceilings below are arithmetic on the profile, after the quick wins
(1-thread refine about 0.43 s; scan about 47 ms at 8 threads; split
0.109 s at 1 thread, about 0.117 s at 8; rest about 0.015 s). They assume
the named parallel efficiency for the split phase and are not measurements.

### What any option must get right

- **Slot and vertex numbering.** `LatticeMesh` appends triangles and vertices
  with `push_back` in the order insertions happen (`lattice_mesh.hpp`,
  `next_slot`, `add_vertex`). In parallel, appended indices must come from a
  prefix sum over the committed insertions in priority order, with the
  vectors sized before the commit step. Otherwise the numbering, and so the
  next round's order and the output, depend on thread timing.
- **The ring.** A split or flip rewrites the back-pointer in each outer
  neighbour (`repoint`). Those slots are written, so they belong to the
  footprint and must be reserved or owned.
- **The skip rule.** Today a mark whose slot was rewritten earlier in the
  round is skipped and rescanned next round (51 % of marks). An option either
  keeps that semantics under a new order, or changes the batch on purpose and
  measures what it costs in triangles. 14b's C2 measured that less greedy
  choices cost 6-27 % more triangles, so this is not a detail.
- **Constraint feet (20b).** The `footed` set is keyed by node. A node on a
  shared edge can be the worst node of both triangles; the edge split
  reserves both, so slot conflicts already serialise the two. No new rule, but
  the suite must cover it.
- **Thread cost.** At about 90 µs per `for_each_chunk` spawn, any option with
  more than a handful of synchronisation steps per round needs a thread team:
  threads started once per `refine` call, synchronised by `std::barrier`,
  joined at the end of the call. That keeps 14 R7's rule that no state
  outlives a call; it is not a pool.

### Option A. Deterministic reservations (Blelloch et al. 2012)

Per round, after the scan, the marks go through sub-rounds:

1. **Footprint, read-only, in parallel.** Each pending mark computes the
   slots its insertion will read or write: the split triangle (and the edge
   neighbour for an edge split), the cavity (triangles reached across
   unconstrained edges for which `must_flip`'s test holds), and the ring.
2. **Reserve.** Each slot takes the minimum priority of the marks that want
   it. No atomics are needed: emit `(slot, priority)` pairs per chunk and
   reduce by slot; or use an atomic fetch-min, whose result is the same
   whatever the timing.
3. **Commit, in parallel.** A mark holding every slot it asked for commits.
   Slot numbers for its appended triangles and vertex come from a prefix sum
   over committers in priority order. Commits touch disjoint slots.
4. **Skip or retry.** A mark whose own slot was written by a commit is skipped
   and rescanned next round, as today. The rest retry.

- **A0, priority = triangle index (L0).** Reproducing today exactly also
  needs a mark to wait until every lower-index mark whose footprint overlaps
  it has resolved, and slots renumbered at the end of the round into today's
  order. Its parallelism is the dependence depth of the index order within a
  round. Indices are spatially correlated (start slots in grid order,
  appended slots at the end), so chains may be long. **Not recommended**
  unless the depth measurement (section 6, question 2) comes back short.
- **A1, priority = a hash of the inserted node's `(row, col)` (L1).** The
  order is random-like, so dependence is shallow (Blelloch, Gu, Shun and Sun).
  The skip rule is kept, with "earlier" meaning lower priority, so the batch
  per round should stay close to today's 42 % of marks. Mesh quality should
  then be close to today's; that is to be measured, not assumed.
- **Expected parallelism.** With conflict degree d, a mark wins its first
  sub-round with probability about 1/(d+1) under random priorities: 11-18 %
  of 9,000-18,700 marks, that is 1,000-3,400 disjoint commits in the first
  sub-round of a big round, shrinking after. That is ample for 8-16 threads;
  the grain is fine but not too fine. The limit is synchronisation per
  sub-round (two barriers), which a thread team keeps to microseconds.
- **The double-work problem.** Step 1 evaluates the same incircle tests that
  Lawson would evaluate again at commit. The flip test is 13.2 % of refine.
  The fix is to commit from the computed cavity directly (Bowyer-Watson
  style: remove the cavity, fan from the new vertex, fix the ring), so each
  test runs once. That makes the commit a new routine. `legalise_around` then
  becomes its **oracle**: on the same mesh and point, the cavity insert and
  split-then-Lawson must give the same triangulation. That equality is
  expected (both remove exactly the triangles whose circle strictly contains
  the point and that are reachable across unconstrained edges, and neither
  acts on a cocircular tie), but it is a claim to test, not to assume,
  including under `must_flip`'s frame-orientation branch (`lawson.hpp`).
- **Rough ceiling at 8 threads**, split phase at 50-70 % parallel efficiency:
  scan 0.047 + split 0.021-0.029 + rest 0.015 + barriers about 0.005 =
  about 0.09-0.10 s, that is 4.3-4.8x over 1 thread and about 2.3-2.5x
  faster than today's 8-thread time. The scan's own 5x then limits.
- **Complexity:** high. Footprint walk, reservation, cavity commit with
  constraints and edge splits, prefix-sum numbering, thread team. Probably two
  PRs.
- **Level:** L1 (A1), or L0 (A0) at a higher cost.

### Option B. Geometric multicolouring (Chernikov and Chrisochoides 2004)

Ola's "multicolouring". Cut the lattice into fixed B x B-node blocks and
colour them in a 2 x 2 pattern, so blocks of one colour are at least B nodes
apart. Each round runs four colour phases. In a phase, blocks of that colour
are processed in parallel, each block's marks serially in index order, with
today's `split` plus `legalise_around`.

- **Independence is not guaranteed; it must be checked.** Chernikov and
  Chrisochoides derive a safe distance from the quality bound. Ours has none:
  a footprint's extent depends on the data, and early rounds work on start
  triangles up to 40 nodes across (stride 40, 14b C1), larger in flat sea.
  So each insertion needs its footprint computed first (as in option A, step
  1) and must be deferred if it leaves its block's halo. Deferred marks go to
  a serial tail after the four phases, or to the next round.
- **Against the density.** 382 64x64 blocks carry insertions over the run;
  at B = 32 about four times as many, some 380 per colour, which is ample for
  16 threads with dynamic scheduling. Per-block load is uneven: in round 14 a
  64-block takes up to 237 insertions against a mean of about 53 (profile,
  section 3). Dynamic scheduling across hundreds of blocks absorbs that.
  What is **not** known is the share of footprints that leave a halo for a
  given B (section 6, question 6).
- **Numbering.** Each block appends into a reserved range sized by its mark
  count (at most two triangles and one vertex per insertion); ranges are
  compacted at the end of the phase by a prefix sum, and references remapped
  in the slots written that phase.
- **Rough ceiling:** similar to A if few footprints cross; worse if many do,
  because the serial tail grows. Four barriers a round, so synchronisation is
  cheap even with a spawn per phase (though a team is still better).
- **Complexity:** medium. The footprint walk is still needed to decide
  deferral; the commit reuses today's code. Numbering and remap are new.
- **Level:** L1, because the blocks are fixed in DEM nodes and each block's
  order is index order. Not L0: the colour-major order changes which marks
  are skipped.

### Option C. Batch insertion, then parallel Lawson flipping (Qi, Cao and Tan 2012)

Our round already takes one point per triangle, which is the GPU methods'
insertion rule. Change the serial phase to:

1. **Insert the batch.** A fan writes `T` and repoints `T`'s neighbours; an
   edge split writes `T` and `U` and repoints theirs. So two marks collide
   when one's split triangles are the other's split triangles *or its ring*:
   two edge splits of the same edge, and also two fans in adjacent triangles
   (each repoints the other's slot). Colliding marks are resolved by a fixed
   rule (the lower key splits in this step, the other in a second step or
   next round), or the repoints are deferred to a separate pass after all
   splits. Slots are numbered by a prefix sum in index order. Which of these
   is cheaper is a 21d design question.
2. **Flip to CDT in parallel rounds.** Collect the non-Delaunay edges among
   those written. In each flip round, choose a set of edges no two of which
   share a triangle, by a fixed rule (an edge flips if its key is the smallest
   among candidate edges of its two triangles), and flip them in parallel.
   Push the four outer edges of each flip for the next flip round. Repeat
   until none remains. Termination is Lawson's argument, which does not depend
   on order (14b R4, part 1).
3. Everything written is touched and rescanned next round, as today.

- **Against the density.** Parallelism is per edge, not per insertion:
  thousands of candidate edges per flip round, each flip about 250 ns with
  its incircle test (14b M5). The unknowns are how many flip rounds a batch
  needs, and how many more flips arbitrary-order flipping costs than Lawson
  around each point.
- **The big risk: the batch changes.** Today 51 % of marks are skipped
  because an earlier insertion rewrote their triangle; here nearly every mark
  inserts. That is less greedy. By 14b's C2 evidence it may cost several
  percent more triangles, and fewer rounds. It must be measured before being
  chosen (section 6, question 6b). A thinning rule (defer a mark whose
  triangle is adjacent to a lower-key mark's) moves it back towards today's
  batch at little cost.
- **Rough ceiling:** if the batch costs nothing in triangles, the best of the
  three. The predicate runs once per test, the commit is today's `flip`, and
  the grain is finest. Same arithmetic as A with 60-80 % efficiency: about
  0.08-0.09 s at 8 threads.
- **Complexity:** medium. Batch split with numbering, a deterministic
  parallel flip loop, a thread team. `split_*` and `flip` are reused.
- **Level:** L1. The final mesh is the CDT of the round's vertex set, which
  is unique except at cocircular ties. A tie is never flipped, so which
  diagonal survives is decided by the flip schedule, which is fixed by edge
  keys, not by timing. Not L0.

### Option D. Domain decomposition (Linardakis and Chrisochoides 2006; Galtier and George 1996)

Cut the domain into subdomains along artificial constraint edges, refine each
as its own mesh, and merge.

- **The Delaunay property across seams.** Artificial separators are
  constraints, so the union is constrained Delaunay with respect to *them*,
  not only to the input's constraints. Restoring the promise means removing
  the separators after the merge, legalising across them, and rescanning what
  changed: a serial or option-C step along every seam, followed by more
  rounds.
- **Separator refinement.** A worst node on a separator splits that separator,
  which the subdomain on the other side must see. Linardakis and Chrisochoides
  avoid this by pre-refining separators from a sizing function; we have no
  sizing function, only the DEM error, known after scanning.
- **Against the density.** DD's grain is coarse: few synchronisations, good
  cache locality. The density data says finer grains have plenty of work, so
  DD's advantage here is not parallelism.
- **Where DD does belong: very large areas.** Ola's ruling names them. When
  the DEM (a mosaic of the 254-tile archive, `ROADMAP.md` 16b and 18 R6) does
  not fit in memory, subdomains that are refined one or a few at a time are
  the natural shape, and seams are the price. That is the parked mosaic
  increment's problem, not this one's.
- **Level:** L1 if the partition is fixed in DEM nodes; L2 if one subdomain
  per thread.
- **Complexity:** high (seams, merge, re-legalisation), for a gain the other
  options reach more cheaply on one machine.

### Others considered

- **Lock-based concurrent insertion** (Kohout et al.; Batista et al.;
  Galois). L3 only. Not recommended (section 4).
- **Pipelining the next round's scan with this round's splits.** A slot
  written early in the round can be rewritten later in it, so its scan
  cannot start before the round ends without the footprint machinery of
  option A. No cheaper than A.
- **A cheaper serial insert** (Bowyer-Watson in place of split plus Lawson,
  serially). It could cut writes, since a flip rewrites two slots and
  repoints two, but the gain is not estimable from the profile. It is the
  commit routine of option A anyway.

### Recommended order

1. **21a and 21b**, the quick wins. Bit-identical, no ruling needed, and they
   set the baseline every option is measured against.
2. **21c, measurement only** (`@perf`, no production code): the questions in
   section 6, especially footprint extents, the dependence depth of the index
   order, and a scratch simulation of option C's batch (triangles, rounds,
   flip rounds, flips) and of A1's (triangles, rounds, sub-rounds).
3. **21d, the parallel serial phase:** **option C if** 21c shows its batch
   costs at most about 2 % more triangles at 1 m (or the thinning rule brings
   it there); **otherwise option A1.** Option B is a scheduling variant of A
   that still needs A's footprint walk, so it is chosen over A1 only if 21c
   shows almost no footprint crosses a 32-node halo. Option D waits for the
   mosaic.

## 6. @perf's six questions

The six questions are in `@perf`'s handback of the profile round
(2026-09-27), not in the committed README, so they are copied here verbatim.
Each is answered or placed.

1. *"The determinism contract fixes the output by triangle-index order. Would
   a different order be acceptable if it still does not depend on the thread
   count? Any colouring or domain decomposition depends on that answer."*
   **Placed with Ola.** His ruling allows a different order; whether
   thread-count independence is required is Q1 in section 8. The architect
   recommends L1 (section 4).

2. *"Insertions are 3 nodes apart, and each conflicts with 4-8 others in a
   round. What does the parallel Delaunay refinement literature say about the
   granularity of independent sets at that density? Does lattice
   cocircularity change the picture?"*
   **Answered, with one measurement to add.** Deterministic reservations work
   at this grain: with random-like priorities about 1/(d+1) of the pending
   marks win a sub-round, which here is 1,000-3,400 disjoint commits in the
   first sub-round of a big round (section 5, option A). The binding limit is
   synchronisation per sub-round, not the supply of independent work, hence
   the thread team. The geometric colouring literature (Chernikov and
   Chrisochoides) works at a coarser grain and needs a distance bound we do
   not have (option B). Cocircularity does not change the conflict counts
   much (a tie is not flipped, so it stops propagation rather than extending
   it). It matters for determinism: with ties, the triangulation of a vertex
   set is not unique, so any change of order can change which diagonal
   survives. That is why L0 is expensive and why L1 needs a fixed tie rule.
   **Measurement (21c):** the dependence depth of today's index order within
   each round. Build the DAG with an edge j -> i when j < i in index order and
   their footprints share a slot, and report the longest path per round. A
   short depth would make option A0 (L0) worth a second look; a long one
   settles it.

3. *"Every exact-path incircle test is a tie between DEM nodes. Is there a
   cheaper exact decision for quads whose corners are all nodes? That is 5 %
   of refine and about 17 % of the serial phase."*
   **Answered: yes,** an int64 lattice incircle under four stated conditions,
   bit-identical to today (QW2, 21b). The premise is still an inference; the
   classification measurement in QW2 comes first.

4. *"The scan loses 16.5 ms of 61.8 ms to imbalance at 8 threads. Is the
   chunking free to change, given that results are written per slot? And does
   the 5.6x rescan amplification come from the deferral rule?"*
   **Chunking: yes, freely** (QW1). The scan is pure and each result slot has
   one writer, so the partition cannot change a result.
   **Amplification: partly, and mostly not avoidably.** A triangle an
   insertion rewrites is a new triangle and must be scanned; that is most of
   the 5.6x, and the deferral rule only moves it by a round. The avoidable
   part is rescanning a skipped mark whose slot nobody wrote (QW4).
   **Measurement:** per round, scanned nodes split into slots written in the
   previous round and skipped-and-unwritten slots (QW4).

5. *"The split phase slows by about 8 % after a multi-threaded scan. Should
   that be profiled before any design relies on the 1-thread serial figure?"*
   **No separate profile now; use the right figure instead.** The designs in
   section 5 use the 8-thread split time (0.137 s), not the 1-thread one, so
   they do not rely on the slowdown's cause. Profiling it without hardware
   counters (no `xctrace`, profile README) would give another inference.
   **Placed in 21c:** re-measure the split time at 1 and 8 threads after 21a
   and 21b land, since both change what the split phase runs after.

6. *"What partition-boundary data would a domain decomposition need? I could
   measure footprint coordinates as a next step."*
   **Yes, and for colouring as much as for DD.** Measurement for 21c, from the
   instrumentation patch extended with vertex coordinates:
   - (a) per insertion, the bounding box in nodes of all vertices of its
     footprint slots; the distribution of its larger side (p50, p90, p99,
     max) per round;
   - per insertion, whether that box leaves its block's halo, for B = 16, 32,
     64 and 128 and a halo of B/2; the share per round;
   - per round and B, insertions per block (max and mean), for load balance;
   - (b) a scratch simulation (never in the tree) of option C's batch: for
     the 1 m quarter circle, triangles, rounds, flips, flip rounds per refine
     round, and achieved max error, against today's; the same with the
     thinning rule; and option A1's batch (hashed priority, skip rule kept):
     triangles, rounds, and sub-rounds per round.

## 7. Proposed increments, LOC and invariant-critical suites

Counted in `CLAUDE.md` §2's unit. Estimates, not measurements; the worst
overrun recorded so far is +99 % for a whole increment (6a shipped 467 lines
against ~235, `06-cdt-viewer.md:626`) and +116 % for one file (`scene.py`,
`05b-noder-driver.md:379-382`); the table gives each at +66 % (increment 17).
At +99 % the conclusions hold: A1 (~450) comes to ~895 and is split in any
case; C (~320) comes to ~637, under 700. The per-file factor does not apply to
a whole increment, but at +116 % C would be ~691, only just under, so 21d on C
is worth counting early.

| increment | what | est. | at +66 % | determinism | invariant-critical suite (mutation round) |
|---|---|---|---|---|---|
| **21a** | QW1 dynamic scan scheduling with a small-round inline threshold (`parallel_util/`), QW3 merge, `legalise_around`'s caller-owned stack; QW4 only if measured worth it (about +20) | ~60 | ~100 | bit-identical | `test_refinement_chunks_dynamic`, a new suite for the dynamic scheduler (`test_refinement_chunks` is unchanged): every index visited exactly once for every n, thread count and block size. A dropped or doubled block leaves a stale scan result, which breaks the tolerance guarantee silently, so this is where the guarantee is decided |
| **21b** | QW2 `lattice_incircle` in `mesh/`, the exact-frame check, the call in `must_flip` | ~50 | ~85 | bit-identical where it answers | a new `test_mesh_lattice_incircle`: agreement with `DetriaExact` on the frame doubles for random and adversarial node quads (grid rectangles, other cocircular lattice quads such as points on a circle of radius 5, near-overflow differences at 2^14, `dx != dy` and inexact `dx` returning `nullopt`, off-node corners returning `nullopt`). Mutants: bound at 2^15, `dx != dy` not refused, the exact-frame check dropped, a sign flip |
| **21c** | measurement only, `@perf`; evidence under `docs/benchmarks/<date>/` | 0 | 0 | — | none |
| **21d** (option C) | thread team (`std::barrier`, per call), batch split with prefix-sum numbering, deterministic parallel flip rounds, the round loop | ~320 | ~530 | L1 | a new `test_mesh_parallel_lawson`: CDT property, conformity and a flip-count bound on random lattice meshes with cocircular ties and constraints, and identical output for threads 1, 2, 7 and hardware concurrency. Plus `prop_refinement_refine`'s tolerance oracle re-run, with the mutant "a flipped slot not marked touched" |
| **21d** (option A1) | thread team, footprint walk, reservation, cavity commit, prefix-sum numbering, the round loop | ~450 | ~750 | L1 | a new `test_mesh_cavity_insert` with `legalise_around` as its oracle (same triangulation on the same mesh and point, ties included), plus the reservation's disjointness and the same determinism and tolerance suites. **Over the ceiling at the worst bias, so it would be cut in two**: the cavity insert serially first (bit-identical to split-plus-Lawson up to the equality claim), then reservations |

Every 21d option changes the output once. What that does to the existing
suites:

- **14's T6** (`prop_refinement_refine`, "T6: the output is bit-identical for
  1, 2, 7 and all threads", and T16 for off-node rings) compares thread
  counts with each other, not with a stored mesh. So do
  `prop_refinement_quality.cpp:233` and
  `prop_refinement_constraint_feet.cpp:587`, which join the same contract. **Under L1 it stays as it
  is and becomes 21d's determinism test.** Under L2 or L3 it would have to be
  weakened or dropped, which is a concrete cost of those levels.
- **18's T3 golden digests** (`tests/python/test_refine_golden.py`) say that
  no commit may update them to agree with new code. 21d is the deliberate
  exception: `@tester` retires them as "unchanged since increment 17" and
  records new digests from 21d's reviewed output, in their own commit with
  the reason, so they guard against the next unintended change. Q1 in
  section 8 is what authorises that.
- Other tests that pin exact counts of a refined mesh (for example 14's T2
  under 14b's amendment) are for `@tester` to find in the red step.

Every new suite joins the TSan job's list in `.github/workflows/main.yaml`.

**Acceptance for 21a, 21b and 21d** is `@perf`'s run
(`docs/increments/README.md`, "Acceptance"): the 1 m benchmark and the thread
sweep from `tools/bench.py`, one power state against the same one (21a ran on
AC, base and 21a back to back: `docs/benchmarks/2026-09-27/21a-acceptance.md`). For 21a and 21b the mesh
hash must be unchanged. For 21d the mesh hash changes by design, so the
comparison is time, triangle count, worst angle and max degree (Q3).

## 8. Questions for Ola (ruled 2026-09-27; see "Ola's rulings")

**Q1. Which determinism level?** Your ruling allows leaving bit-identical
output. It does not say whether the mesh may then depend on the thread count
or vary between runs. The levels are in section 4.
*Recommendation: L1*, deterministic and independent of the thread count, in a
new order. It costs no speed in the recommended options, keeps failures
replayable, and keeps benchmark meshes citable for a publication.

**Q2. One path for every input, or a switch for large areas?** "When we are
encountering geometries where this is important, ie very large areas" can be
read as keeping today's path for small inputs.
*Recommendation: one path.* Every mesh, small ones included, changes once
when 21d lands; after that there is one algorithm to test, and no threshold
to tune. With L1 the 1-thread run is as reproducible as today's.

**Q3. How many more triangles may the parallel serial phase cost?** Option C
changes which points each round inserts, and that may cost triangles for the
same tolerance (14b's C2 saw 6-27 % for other less greedy rules).
*Recommendation:* at most 2 % more triangles at 1 m on the quarter circle,
with worst angle and max degree no worse, measured by `@perf` in 21c before
21d is chosen. If C misses it, option A1 keeps today's batch rule.

**Q4. Go ahead with 21a and 21b now, before the answers above?** They are
bit-identical and need none of the rulings.
*Recommendation: yes.* 21b after the tie-classification measurement (QW2),
which is a short `@perf` task on the existing instrumentation patch.

**Q5. Domain decomposition deferred to the large-area (mosaic) work?** DD's
advantage is coarse grain and memory, not parallel work, of which the
density data shows plenty at a finer grain.
*Recommendation: yes.* Design it with the refreshed increment 15, where the
DEM no longer fits in memory and subdomains are the natural unit.
