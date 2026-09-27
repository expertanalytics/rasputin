# Increment 21c: measurements for choosing 21d's option (@perf, 2026-09-27)

Status: **done**, with a later section "A0 and the scaling model". Measurement only; no production code changed. The scratch
simulations are `scripts/sim.patch`, which is not applied and never committed
to `include/`, like the serial profile's `instrument.patch`.

Asked by the main session for `docs/increments/21-parallel-refine.md`:
section 5 (options A1, B, C, D), section 6 (the items placed in 21c) and
section 7 item 2. Ola's rule (Q3): **at most 2 % more triangles at 1 m on the
quarter circle, with worst angle and max degree no worse.**

## Findings

1. **Option C misses the 2 % rule by 3 to 5 times.** On the quarter circle,
   batch insertion followed by parallel flipping gives **+10.0 %** triangles.
   With the thinning rule of section 5 it gives +7.5 %. A greedy independent
   set gives +7.5 %, and two-hop thinning +6.5 %. The worst angle is worse in
   every variant (0.18-0.30 deg against 0.3955). Max degree is 18-20 against
   18. The tile shows the same: +5.9 % to +9.2 %. The flipping itself is cheap
   to parallelise: 7 flip rounds per refine round (median, at most 11), with
   about 47,000 independent flips in the first flip round of the largest
   round. So the triangle count is what rules C out, not its parallelism.
2. **Option A1 meets the triangle rule: +0.86 % to +1.27 %** on the quarter
   circle over four hash seeds, and +1.21 % to +1.41 % on the tile. Max degree
   is 17-18 against 18 (tile: 74 against 74). **The worst angle depends on the
   seed.** On the quarter circle, 2 of 4 seeds give a worse worst angle
   (0.2717 deg against 0.3955) and 2 give a better one (0.4458 and 0.4591). On
   the tile all four seeds are better (0.78-1.25 deg against 0.63). In each
   case the worst angle is one interior sliver with three node corners, and
   the count of triangles under 1 deg is 8-13 for A1 against 12 for today.
   Whether "worst angle no worse" should hold for one triangle, or be read
   from that distribution, is for Ola (see "Worst angle").
3. **A1 gave the same triangulation as the serial loop run in hash order**
   ("Hser"): the same triangles, rounds, insertions, flips and quality on all
   8 runs. Two runs were compared triangle by triangle and are the same mesh
   up to numbering. So **A1's +1 % comes from changing the order, not from the
   reservations**. That is measured on these inputs, not proven in general.
   A1 needs a median of 14 sub-rounds per refine round on the quarter circle
   (at most 23, 447 in total). The first sub-round of a big round commits
   5-19 % of the marks.
4. **Option B at B = 32 is not "almost no crossing".** 7.4 % of insertions on
   the quarter circle (7.9 % on the tile) have a footprint that leaves its
   block's 16-node halo. In the big rounds this is 5.2 % (tile: 3.9 %). At
   B = 64 it is 2.8 % (tile 1.2 %), and at B = 128 it is 1.0 % (tile 0.04 %).
   Early rounds cross almost always: 99 % in round 1 at B = 32.
5. **The dependence depth of today's index order is short**: 13-19 in the big
   rounds 9-19 of the quarter circle, which insert 9,300-18,700 each. The
   largest is 24, in round 2. On the tile the depth is 80-101 in rounds 1-3 and
   at most 28 afterwards. Section 6 Q2 said that a short depth "would make
   option A0 (L0) worth a second look". This depth is an upper bound: it
   counts read-read sharing of a slot as a dependence. Whether a parallel
   schedule can reach it is a separate question (see "What this does not
   show").
6. **The split phase at the current tree** (`7f688aa`, battery) takes 92.6 ms at
   1 thread and 96.4 / 98.2 ms at 8 threads (quarter / tile). That is 4-6 %
   slower at 8 threads. Refine takes 431 / 466 ms at 1 thread and 162 / 186 ms
   at 8. The serial share at 1 thread is 26 % on the quarter circle, against
   34 % in the profile.
7. **Side finding, measured: `tools/bench.py`'s tolerance check cannot see a
   stale scan result.** `within_tolerance` is refine's own reported
   `max_error <= tolerance`. With a planted defect (flipped slots not marked
   touched), bench.py reports `within_tolerance: True` while a full rescan of
   the delivered mesh finds 6,012 triangles over tolerance, the worst at
   19.5 m. 21d changes that bookkeeping, so its acceptance needs an
   independent rescan. Adding one is a change to `tools/bench.py` and follows
   the TDD loop.

## Power state

**Battery** for everything: 79 % to 78 %, discharging
(`data/pmset_sims_*.txt`, `data/split_phases.txt.pmset`). The counts and
quality figures do not depend on power; the timings in section 6 do.
Machine: Apple M1 Max (8 P + 2 E), macOS 27.0.

## Method

- **Tree.** `7f688aa` (21a and 21b, branch `increment21c-measurements`). The
  simulations ran in a scratch `git worktree` of `7f688aa` with
  `scripts/sim.patch` applied. The patch adds
  `include/terrain/refinement/sim21c.hpp` and replaces the round loop of
  `refine.hpp`. With `RASPUTIN_SIM` unset, the loop is today's.
  `RASPUTIN_SIM=today` also gives today's mesh, sha256 `1e531976…` on the
  quarter circle. That is the 21b acceptance hash, re-checked after the last
  edit to the patch.
- **Builds.** Both builds are Release, made by `tools/bench.py`'s `build()`
  (`-O3 -DNDEBUG`, AppleClang). The scratch `_core` hashes to `ea1c5449…`
  for the final patch, which made the Hser, plant and re-check runs. The
  other simulation runs used `7d703cda…`, a build of the same patch before
  the `Hser` mode and the plants were added; quarter C, Cthin and A1 rerun on
  `ea1c5449…` give the same mesh sha256 as the stored runs. The unpatched
  `_core` hashes to `2b69e078…`; it was used for the split timings.
- **Benchmark.** DEM `tests/fixtures/dem_archive/7908_3_10m_z33.tif`, tolerance
  1, all other CLI defaults: constraint feet on, and the default start
  quality. Domains: `docs/benchmarks/2026-09-26/quarter.geojson`, and the tile
  (no domain).
- **Driver.** `scripts/sim_driver.py` runs one `rasputin mesh` through
  `tools/bench.py _child` with the scratch pkg, and writes ASCII VTK. It then
  applies bench.py's own `quality()`, the same function the acceptance run
  uses: worst angle, share under 1 deg, max degree, and the constrained
  Delaunay check with an exact fallback. The result is merged with the
  simulation's JSON. `scripts/run_sims.sh` runs every mode on both domains.
  `scripts/analyse21c.py` makes `data/tables.md`, which holds every table
  below plus all per-round rows.
- **Tolerance check, independent.** After the loop, every simulation rescans
  **every** triangle of the final mesh and reports the largest error and the
  number of triangles that still need a split ("full rescan" below). It does
  not rely on the loop's `touched` bookkeeping.

### The modes

| mode | what it simulates |
|---|---|
| `today` | today's loop: index order, skip rule. With `RASPUTIN_SIM_FOOT=1` it also records each insertion's footprint, for questions 2 and 6(a) |
| `C` | option C. **Batch:** every mark of the round is planned on the round-start mesh. A mark's split set is its triangle, plus the edge neighbour for an edge split. A mark inserts when the lowest slot index among the marks wanting each slot of its split set is its own; otherwise it is deferred to the next round. Splits have no legalisation. **Flip rounds:** the candidates are every unconstrained edge of every triangle written so far. Each is tested with today's `must_flip`, on the mesh as it is at the start of the flip round. An edge flips when its key (a hash of its two vertex indices) is the smallest among the flippable edges of both its triangles. The loop repeats over the triangles of every flippable edge until no edge is flippable. Everything written is touched |
| `Cthin` | C, plus section 5's thinning rule: a mark is deferred if an edge-adjacent triangle is marked with a lower slot index |
| `Cmis` | C, plus a greedy independent set: a mark is deferred if an edge-adjacent mark with a lower slot index inserts |
| `Cthin2` | C, with thinning over two edge hops |
| `A1` | option A1, deterministic reservations. The priority is `(splitmix64((row << 32 \| col) ^ seed), slot)` of the worst node. Each sub-round goes through four steps. (1) The skip rule: a mark whose slot was written this round is skipped, and so is an edge split whose neighbour was written; both are rescanned next round, as today. (2) The footprint on the current mesh: the Bowyer-Watson cavity, meaning the base triangle(s) plus the triangles whose circle strictly contains the point, reached across unconstrained edges, plus its ring of edge neighbours. (3) A mark commits when it holds the minimum priority on every footprint slot. (4) Winners commit with today's `split` and `legalise_around`, in priority order; losers retry |
| `Hser` | today's serial loop and skip rule, with the marks taken in A1's priority order: the serial L1 twin of A1 |

### Why the serial simulation is the parallel result

- **C, batch.** Winning split sets are disjoint. A split writes only its own
  split set and appends new slots; it rewrites a neighbour's pointer and
  nothing else of the neighbour. So the final topology does not depend on the
  order of the splits. Committing them in index order numbers the new slots
  as a prefix sum in index order would.
- **C, flips.** Every decision is taken on the mesh at the start of the flip
  round, before any flip of that round. The selected edges share no
  triangle; this was checked on every flip round and counted as
  "consistency" (0 everywhere). A flip reads and writes only its own two
  triangles, plus the pointers of the four outer neighbours. So disjoint flips
  commute.
- **A1.** Winners' footprints are disjoint, and a footprint contains
  everything a commit reads or writes, provided the cavity contains every slot
  Lawson writes. That was checked at every commit: every existing slot that
  `split` or a flip wrote must lie in the cavity computed beforehand. The
  count of misses is "consistency", 0 on all 16 A1 runs. Committing in
  priority order numbers the new slots as a prefix sum in priority order
  would.

### The checks can fail

Each check was run against a planted defect (`RASPUTIN_SIM_PLANT`, in the
patch; `data/plants/`), quarter circle:

| plant | what it breaks | what caught it |
|---|---|---|
| `nocavity` (A1) | the cavity is the base triangles only | consistency 527,101 misses; the mesh then also has 40 Delaunay violations and 3 triangles over tolerance (max 6.83) |
| `notouch` (C) | flipped slots not marked touched | full rescan: 6,012 triangles over tolerance, max 19.53. **bench.py's `within_tolerance` said True**, because refine's reported `max_error` was 0.99999 (finding 7) |
| `noflip` (C) | no flip rounds | Delaunay check: 139,301 violations of 995,557 edges |

`scripts/same_mesh.py` compares two meshes triangle by triangle, as sets of
coordinate triples. It reports SAME for A1 against Hser (quarter seed 0, tile
seed 2), SAME for a mesh against itself, and DIFFERENT for today against A1
(339,917 of today's triangles are not in A1's mesh).

## 1. The 2 % rule, per option

Quarter circle; today is 428,217 triangles, worst angle 0.3955 deg and max
degree 18. Every row is a valid constrained Delaunay mesh within tolerance:
the full rescan's max error is 1.0000 with 0 triangles over tolerance, there
are 0 Delaunay violations, and consistency is 0.

| option | variant | triangles | vs today | 2 % rule | worst angle | max degree | rounds |
|---|---|---:|---:|---|---:|---:|---:|
| C | batch, no thinning | 471,045 | **+10.00 %** | fails | 0.1889 (worse) | 20 (worse) | 20 |
| C | thinning (section 5) | 460,453 | **+7.53 %** | fails | 0.2955 (worse) | 18 | 31 |
| C | greedy independent set | 460,443 | **+7.53 %** | fails | 0.2379 (worse) | 20 (worse) | 24 |
| C | two-hop thinning | 455,911 | **+6.47 %** | fails | 0.1841 (worse) | 20 (worse) | 50 |
| A1 | seed 0 | 433,669 | **+1.27 %** | meets | 0.2717 (worse) | 17 | 36 |
| A1 | seed 1 | 433,609 | **+1.26 %** | meets | 0.4591 | 18 | 37 |
| A1 | seed 2 | 433,079 | **+1.14 %** | meets | 0.4458 | 18 | 38 |
| A1 | seed 3 | 431,879 | **+0.86 %** | meets | 0.2717 (worse) | 18 | 38 |

Tile, against today's 472,374 triangles, worst angle 0.6296 and max degree 74:
C gives +9.20 %, thinning +6.76 %, the independent set +6.98 % and two-hop
thinning +5.92 %. A1 gives +1.21 % to +1.41 %, with worst angle 0.78-1.25 and
max degree 74. The full table, with flips, inserted counts, the share under
1 deg and the Delaunay edge counts, is `data/tables.md` ("Summary").

**Why C costs triangles (measured counts, inferred cause).** The batch
changes: today inserts 213,464 of 509,029 marks (42 %) and skips 51 % because
an earlier insertion rewrote their slot. C inserts 234,878 of 340,966 (69 %)
in half the rounds. More points go in per round before the errors are
re-measured, so the refinement is less greedy. 14b's C2 saw the same effect.
Thinning brings C back only part of the way: +6.5 % at two hops, and the tile
then needs 167 rounds.

### Worst angle

The worst angle is one triangle. The five smallest angles of each quarter
mesh (`scripts/worst_angles.py`) are interior slivers whose three corners lie
on the 10 m grid, with no constraint edge:

- today: 0.3955, 0.4341, 0.4458, 0.6296, 0.6708;
- A1 seed 0 (and Hser): 0.2717, 0.6679, 0.6913, 0.7332, 0.7751.

A1 seed 0 has one triangle below today's worst, and its second worst is
better than today's second worst. That "on the 10 m grid" means "at a DEM
node" is inferred, not checked: it assumes the grid origin is a multiple of
10 m. Under the literal rule, A1 passes or fails on the quarter circle
depending on the hash seed. A seed chosen for this one benchmark would be
tuning to the test.

## 2. Option C: flip rounds

| domain | variant | refine rounds | flip rounds, total | per refine round (median / max) | flips | largest round: candidate edges / flips in its first flip round |
|---|---|---:|---:|---|---:|---|
| quarter | C | 20 | 130 | 7 / 9 | 412,965 | 233,681 / 47,539 (round 9) |
| quarter | thinning | 31 | 192 | 7 / 9 | 422,051 | 154,147 / 30,889 (round 13) |
| tile | C | 45 | 174 | 3.5 / 8 | 410,455 | 278,757 / 57,253 (round 7) |

C does 7 % fewer flips than today, 412,965 against 445,657, for 10 % more
triangles.

## 3. Option A1: sub-rounds

| domain | seed | rounds | sub-rounds, total | per round (median / max) | first sub-round's winners, big rounds (>= 5,000 marks) | sub-rounds until 90 % of a big round's inserts | footprint (cavity + ring): mean, p99, max slots |
|---|---|---:|---:|---|---|---|---|
| quarter | 0 | 36 | 447 | 14 / 21 | 5-19 % of marks | 4-11 | 9.0, 16, 34 |
| quarter | 1-3 | 37-38 | 461-476 | 13-15.5 / 20-23 | 4-19 % | 4-11 | 9.0, 16, 30-34 |
| tile | 0-3 | 50-53 | 376-389 | 3-5 / 19-21 | 5-26 % | 3-11 | 8.9, 16, 28-30 |

In quarter round 14 (44,823 marks, 18,357 inserted), the first sub-rounds
commit 3,821, 2,945, 2,542, 2,136, 1,816… and 20 sub-rounds end the round. A
sub-round count includes a final sub-round that only retires skipped marks.
Per-round lists are in `data/sims/*_A1*.json` (`subround_winners`). Section 5
predicted that about 1/(d+1), 11-18 %, of the marks would win the first
sub-round. The measured share is 5-19 %.

**Inferred, not measured:** with a thread team and two barriers per sub-round,
the quarter circle's 447 sub-rounds are about 900 barriers per refine call.

## 4. Footprint extents: option B and the DD boundary data (section 6, Q6a)

These are measured in today's serial order, on the mesh that earlier
insertions left, as the profile measured slots. A footprint is every slot
the insertion wrote (appended ones included) plus their edge neighbours
afterwards; for the quarter circle this is p50 10, p90 14, p99 18, max 36
slots, which matches the profile's footprint (p50 10, p99 18, max 36). The
extent is the bounding box, in nodes, of every vertex of those slots. The
larger side over all insertions is p50 10, p90 32, p99 142, max 1,548 on the
quarter circle, and p50 10, p90 39, p99 81, max 163 on the tile. A block is
the B × B-node block that holds the inserted point. "Leaves the halo" means
the box reaches more than B/2 nodes outside that block. "Crosses" means the
box touches or crosses the block's boundary.

| B | quarter: leaves the B/2 halo | tile: leaves the B/2 halo | quarter: crosses the boundary | tile: crosses the boundary |
|---:|---:|---:|---:|---:|
| 16 | 18.8 % (big rounds 15.8 %) | 20.8 % (16.1 %) | 81.9 % | 82.7 % |
| 32 | **7.4 %** (5.2 %) | **7.9 %** (3.9 %) | 57.1 % | 58.9 % |
| 64 | 2.8 % (1.8 %) | 1.2 % (0.27 %) | 34.9 % | 37.3 % |
| 128 | 1.0 % (0.55 %) | 0.04 % (0.01 %) | 20.9 % | 23.4 % |
| 256 | – | – | 9.9 % | 12.6 % |
| 512 | – | – | 4.3 % | 7.4 % |
| 1024 | – | – | 1.3 % | 4.6 % |

"Big rounds" means rounds 9-19 on the quarter circle and rounds with at least
5,000 insertions on the tile. The per-round shares, and insertions per block
(max and mean over the blocks hit) for B = 32 and 64, are in
`data/tables.md`; the other values of B are in `data/sims/*_today_foot.json`.
At B = 32 in quarter round 14, the busiest block takes 89 insertions against
a mean of 17.0; at B = 64 it is 237 against 52.9, as in the profile.

**For B (evidence, not a design):** a 32-node halo defers about 7 % of all
insertions to a serial tail, and 56-99 % of those in rounds 1-4. At B = 64 or 128
the tail is 1-3 %, with fewer blocks per colour. **For DD:** a seam every
256 nodes is touched by 10-13 % of insertions, and a seam every 1,024 nodes
by 1.3-4.6 %.

## 5. Dependence depth of today's index order (section 6, Q2)

This is section 6's DAG: an edge from j to i when j comes before i in the
serial order of a round and their footprints (as defined in section 4) share
a slot. The figure is the longest path per round, counted over insertions.
Skipped marks write nothing and do not extend a chain. It is computed exactly
in one pass: `depth(i) = 1 + max(last[s] for s in footprint(i))`.

| quarter round | 1 | 2 | 5 | 9 | 12 | 14 | 16 | 19 | 22 | 25 | 30 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| inserted | 231 | 398 | 2,063 | 9,294 | 16,374 | 18,714 | 17,536 | 11,136 | 4,740 | 1,430 | 108 |
| depth | 18 | 24 | 17 | 19 | 17 | 16 | 16 | 13 | 13 | 9 | 7 |

On the tile the depth is 101, 88 and 80 in rounds 1-3, where the start
triangles are in grid order, then 28 in round 4 and 8-21 in the big rounds.
Every round is in `data/tables.md`.

## 6. The split phase at 1 and 8 threads (section 6, Q5)

Measured at `7f688aa` (unpatched, `_core` `2b69e078…`) with the profile's
`prof_driver.py` (`scripts/split_sweep.sh`): 3 interleaved passes of 5 calls,
15 samples per cell, medians in ms, min-max in brackets. Battery, 78 %.

| domain, threads | refine | scan | split | rest |
|---|---:|---:|---:|---:|
| quarter, 1 | 431.1 | 318.1 | 92.6 (88.8-121.9) | 16.2 |
| quarter, 8 | 161.7 | 47.5 | 96.4 (94.1-99.3) | 16.7 |
| tile, 1 | 465.9 | 343.6 | 92.6 (88.5-107.5) | 26.8 |
| tile, 8 | 186.3 | 60.6 | 98.2 (94.7-108.9) | 27.9 |

The split phase is 4-6 % slower after an 8-thread scan (the profile saw 8 %).
The serial part at 1 thread (refine - scan) is 113 ms, 26 % on the quarter
circle, down from 34 % in the profile. At 8 threads the split phase is 60 %
of refine on the quarter circle and 53 % on the tile. Raw samples:
`data/split_phases.txt`.

## What is measured and what is inferred

**Measured:** every count, quality figure and timing in the tables; that
every simulated mesh passes the full rescan and the Delaunay check; that A1
equals Hser on all 8 runs' counts and on 2 runs triangle by triangle; the
depths; the extents; that each check fails on its planted defect.

**Inferred:**

- that the serial simulations equal the parallel executions (the argument
  above; the disjointness it rests on is checked at every step);
- that C costs triangles because it is less greedy;
- that A1 equals Hser in general (it is measured only here);
- the barrier count and cost for A1;
- that the worst-angle triangles' corners are DEM nodes (the grid-origin
  assumption);
- for option B, that footprints measured in serial order stand in for those a
  colour-phase order would produce.

## What this does not show

- **Parallel speed of any option.** The simulations are serial and none was
  timed as a parallel algorithm.
- **Whether A0 reaches the measured depth.** A0 is option A with priority =
  index, bit-identical to today. The depth assumes each insertion's footprint
  is known. In a reservation scheme it is computed on the mesh at the start
  of a sub-round, and it can change when a lower-index mark commits. Also,
  today's numbering (slots appended in serial order) must be reproduced. A0
  was not simulated here; it was later the same day, in "A0 and the scaling
  model" below.
- **Cthin's round count on the tile** (166 rounds against 53) was not
  investigated.
- **The tile's 20 rounds of about 196 insertions** each (rounds 22-41,
  footprint sides of exactly 40 or 80 nodes) are today's behaviour, seen in
  passing; their cause was not examined.

## Not measured

- **QW4** (section 3, and section 6 Q4: scanned nodes split into slots written
  in the previous round and slots that were skipped and not written). It is
  not placed in 21c by the increment file, and was not measured.
- A C variant whose thinning reproduces today's skip rule exactly. That is
  essentially A1's batch.
- Option B and option D as simulations (only their boundary data, in section
  4).
- An AC run of the split timings. There is no AC figure at `7f688aa`.

## Regenerating

1. `git worktree add --detach <wt> 7f688aa` and `git -C <wt> apply
   docs/benchmarks/2026-09-27/21c/scripts/sim.patch`. Then build with
   bench.py's `build()`, for example
   `python -c "import sys; sys.path.insert(0,'tools'); import bench; from pathlib import Path; print(bench.build(bench.make_runner(), Path('<wt>')))"`.
2. `scripts/run_sims.sh <wt>/build-bench/pkg <out>`. Hser runs:
   `sim_driver.py <pkg> Hser <domain> <out> [--seed N]`. Plants: set
   `RASPUTIN_SIM_PLANT=nocavity|notouch|noflip`.
3. `scripts/analyse21c.py <out>` makes the tables.
4. Split timings: build the unpatched tree the same way, then
   `scripts/split_sweep.sh <repo>/build-bench/pkg <file>`.

The meshes (ASCII VTK, about 22 MB each) were written to the session
scratchpad and deleted by the driver, or kept only for the same-mesh and
worst-angle checks. Steps 1-2 regenerate them; `--keep` keeps them.

## A0 and the scaling model (@perf, 2026-09-27, evening)

Asked by Ola: "simulate A0 first, then we decide", and what is known about the
parallel scaling of A1 against C. Measurement and a model only; no production
code. The simulation is `scripts/sim_a0.patch`, applied **on top of**
`scripts/sim.patch`, never committed to `include/`. Data: `data/a0/`.

### Findings

1. **A0 is bit-identical to today on both domains.** Mesh sha256 as
   bench.py hashes it: quarter `1e531976…` (the 21b acceptance hash), tile
   `11741a81…`, equal to today's. Triangles, rounds, insertions and flips are
   equal (quarter 428,217 / 41 / 213,464 / 445,657; tile 464,290 / 53 /
   219,837 / 445,675), and a hash of the whole lattice mesh (triangles,
   neighbours, vertices) is equal **after every round** (41 and 53 rounds).
   No first divergence exists to report. The rescan, Delaunay and footprint
   checks all pass: full rescan max error 0.99998 with 0 over tolerance, 0
   Delaunay violations, 0 consistency misses, 0 order violations.
2. **Identity needs the end-of-round renumbering.** 212,226 of 213,464
   commits on the quarter circle (99.4 %) happen out of index order. Without
   renumbering (plant `norenumber`) the mesh differs from round 1 on.
3. **Identity is measured on these two inputs, not proven.** A footprint is
   computed at the start of a sub-round. A lower-index mark that loses and
   commits later could, in principle, grow a footprint into a region a
   higher-index mark has already committed. The simulation checks that at
   every commit ("order violations": a footprint slot already stamped by a
   higher-index commit this round). The count was 0 on both domains. A
   production A0 would need this check at runtime with a fallback, or a proof.
4. **A0's sub-round profile is close to A1's on the quarter circle. On the
   tile it is worse in the first rounds.**

   | domain | option | rounds | sub-rounds: median / max / total | first sub-round's winners, big rounds | sub-rounds to 90 % of a big round | largest round: winners per sub-round |
   |---|---|---:|---|---|---|---|
   | quarter | A0 | 41 | 13 / 24 / 469 | 4.6-19.6 % | 4-10 | r14, 45,306 marks: 4045, 3152, 2669, 2316, 1921, 1510, 1184, 827, 523, 295, 142, 82, 29, 12, 6, 1 |
   | quarter | A1 (seed 0) | 36 | 14 / 21 / 447 | 4.6-19.4 % | 4-11 | r14, 44,823: 3821, 2945, 2542, 2136, 1816, 1517, 1225, 852, 639, 388, 234, 117, 63, 29, 16, 6, 4, 3, 3, 1 |
   | tile | A0 | 53 | 2.5 / 101 / 592 | 4.1-26.3 % | 3-50 | r10, 58,626: 4916, 3789, 3453, 3087, 2721, 2258, 1726, 1111, 637, 300, 135, 45, 13, 3, 0 |
   | tile | A1 (seed 0) | 51 | 4.5 / 19 / 376 | 5.3-26.1 % | 3-11 | r11, 56,301: 5268, 3986, 3471, 3056, 2470, 1931, 1390, 910, 542, 249, 113, 44, 12, 5, 2, 1 |

   "Big" means at least 5,000 marks. On the tile, A0's rounds 1-3 need 101, 94
   and 86 sub-rounds, with 2,317-6,351 marks. That is the index order's
   dependence depth measured in section 5 (101, 88, 80), where the start
   triangles are in grid order. A1 needs at most 19 there. A sub-round count
   includes a final sub-round that only retires skipped marks, which can have
   0 winners.
5. **Re-evaluation dominates, for A0 and A1 alike.** A mark that loses is
   evaluated again (skip checks, plan, cavity and ring, reservation) in every
   later sub-round. Evaluations per mark: 4.00 (A0 quarter), 4.13 (A1
   quarter), 4.46 (A0 tile), 4.00 (A1 tile), which is 2.0-2.3 million per
   run. At the measured 279 ns each (quarter) that is 568 ms of serial work,
   about 6 times today's whole split phase (92.6 ms).
6. **Synchronisation on this Mac, measured** (`scripts/sync_bench.cpp`,
   battery): starting and joining 8 fresh `std::jthread`s costs 90-189 µs
   (medians of 3 runs; 96 µs median of medians). `parallel_util::for_each_block`
   with a trivial body costs 83-132 µs (86 µs). The design's 90 µs is in that
   range. With a persistent team, one `std::barrier::arrive_and_wait` of 8
   threads costs 0.82-0.96 µs (0.90 µs); a spinning sense-reversing barrier
   costs 0.72-0.85 µs. At 4 threads: spawn 58 µs, barrier 0.50 µs. At 16
   threads (more than the 10 cores) the barrier costs 20 µs.
7. **The model at 8 threads** (split phase, ms; today's serial split is 92.6
   ms on both domains, and today's refine at 8 threads 161.7 / 186.3 ms):

   | option | quarter, fresh threads per step | quarter, team + barrier | tile, fresh threads | tile, team + barrier | 2 % rule |
   |---|---:|---:|---:|---:|---|
   | A0 | 231 | **93** | 276 | **102** | meets (identical) |
   | A1 | 216 | **88** | 188 | **81** | meets (+1.3 %) |
   | C | 92 | **50** | 114 | **52** | fails (+10 %) |
   | A0, each mark evaluated once (lower bound) | 177 | 39 | 218 | 44 | meets |
   | A1, each mark evaluated once (lower bound) | 159 | 31 | 137 | 30 | meets |

   Refine at 8 threads, adding the measured 8-thread scan and rest (section 6)
   to the team figures: A0 157 / 190 ms, A1 152 / 169 ms, C 114 / 140 ms; the
   lower bounds 104 / 132 (A0) and 95 / 118 (A1). Against today's 161.7 /
   186.3 ms, **as simulated, A0 gains 3 % on the quarter circle and loses 2 %
   on the tile, A1 gains 6-9 %, and C gains 25-30 %.** Most of the gap is re-evaluation. The rows "each mark evaluated
   once" are a bound: they assume a loser keeps its footprint until a slot in
   it is written. That rule was not simulated and its cost is not known.
   **With fresh threads per step, no reservation variant beats today's serial
   split.** 1,100-1,800 synchronisation steps at 96 µs are 108-176 ms.

### The model (all of it is a model, not a measurement)

`scripts/model21c.py data/a0 data/a0/sync_bench.txt` writes `data/a0/model.md`
(4, 8 and 16 threads, both domains). The split phase at T threads is

  sum over steps of ceil(items / T) × per-item cost + syncs × sync cost.

- **Steps and items**, recorded per step by the simulation. A0 and A1: per
  sub-round, an evaluation step (items: marks still pending, `sr_pending`)
  and a commit step (items: winners); A0 also has one renumbering step per
  round. C: per refine round, a batch step (items: marks); per flip round, a
  test step (items: candidate edges, `fr_cand`) and a flip step (items: flips,
  `fr_flips`). The last flip-free test of each refine round is not recorded
  per round; its candidates (total minus recorded) are spread evenly.
- **Syncs per step**, assumed: 3 per sub-round (reserve | check winners and
  count | prefix sum and commit), 1 per renumbering, 2 per C batch, 3 per flip
  round (test | select by minimum key | flip). A step with no items still
  pays its sync. "Fresh threads" charges one spawn and join (96 µs at 8) per
  sync; "team" one barrier (0.90 µs at 8).
- **Per-item costs, measured serially on battery.** Evaluation: 279 ns
  (quarter) and 259 ns (tile), from A0's timer (`t_eval / n_eval`). The same
  routine is used for A1; A1's own timer reads 463-480 ns, because its
  simulation still allocates per candidate and A0's does not. Commit (split
  plus Lawson): today's production split phase per insertion, 92.6 ms /
  213,464 = 434 ns (quarter) and 421 ns (tile). The simulation's own commit
  timer, with its bookkeeping, reads 531-811 ns. A0 renumbering: the
  simulation's whole-mesh copy per round, 1.7 ms (quarter) and 2.1 ms
  (tile) per round, divided by T. A production renumbering would touch only the slots written.
  C: batch 113-118 ns per mark, test 103-104 ns per candidate edge
  (collecting, sorting and deduplicating included), flip 77-87 ns.
- **Not in the model:** load imbalance within a step (every item is charged
  the mean cost), memory bandwidth and cache effects of 8 threads, the E-cores
  (8 P + 2 E), the 4-6 % slower split after an 8-thread scan (section 6), a
  barrier's cost when threads arrive unevenly, and the scan and rest for C and
  A1 (they have 10 % and 1 % more triangles, and fewer rounds; today's
  8-thread figures are used for all). The 16-thread rows run 16 threads on 10
  cores.
- **At 1 thread the model gives 729 (A0), 692 (A1) and 395 (C) ms of split**,
  against today's 92.6. That is the work the parallel schedule adds, before it
  is divided among threads: re-evaluation for A0 and A1, and the test of 3.1
  million candidate edges for C. The C simulation was not written for speed
  (candidate sort per flip round), so C's figures are the simulation's, not a
  design's.

### The checks can fail

Quarter circle, `RASPUTIN_SIM_PLANT` (`data/a0/plants/`):

| plant | what it breaks | what caught it |
|---|---|---|
| `norenumber` | appended slots and vertices keep commit order | sha256 `e7bbaf47…` against `1e531976…`; round hash differs from round 1; 428,117 triangles |
| `noreserve` | in sub-round 1 of each round, the lowest pending mark neither reserves nor wins | order violations: 84; skips caused by a higher-index writer: 11; round hash differs from round 2; 428,209 triangles |
| `nocavity` | the cavity is the base triangles only | consistency: 506,728 misses; order violations: 842; 22 Delaunay violations; full rescan: 1 triangle over tolerance (1.11) |

### Method

- **Tree and build.** A scratch worktree of `7f688aa`, with `scripts/sim.patch`
  and then `scripts/sim_a0.patch` applied, built Release by bench.py's
  `build()` (`-O3 -DNDEBUG`, AppleClang 21). `_core` sha256 `ed06d760…`; a
  second fresh worktree with both patches gives the same sources and the same
  `_core` hash. On this build, quarter and tile A1 (seed 0) and C give the
  same mesh sha256 as the stored 21c runs, so the timers added in
  `sim_a0.patch` did not change what the modes do.
- **What `sim_a0.patch` adds.** Mode `A0`: per sub-round, every pending mark
  in index order applies today's skip rules first (slot written this round:
  skipped, as today, not added to `skipped`; edge split whose neighbour was
  written: added to `skipped`). Then it computes A1's footprint (Bowyer-Watson
  cavity plus ring) on the current mesh and reserves each slot with its index.
  A mark holding every slot wins. Winners commit in index order with today's
  `split` and `legalise_around`. At the end of the round,
  `LatticeMesh::sim_renumber` (scratch) permutes the appended slots and
  vertices into the order a prefix sum in index order gives. It also adds
  per-sub-round and per-flip-round work records, phase timers, and a mesh hash
  per round to `today`, `A0` and `A1`; C's flip selection resets sparsely, and
  `plant()` reads the environment once (the timers were otherwise inflated).
- **Runs.** `scripts/run_a0.sh <pkg> <out>`: today, A0, A1 (seed 0) and C on
  both domains, then the three plants on the quarter circle. The same driver,
  quality checks and full rescan as the rest of 21c.
- **Sync benchmark.** `c++ -std=c++20 -O3 -DNDEBUG -I include
  scripts/sync_bench.cpp`, run for 8, 4 and 16 threads, then 8 twice more
  (`data/a0/sync_bench.txt`): 200 samples for spawn and `for_each_block`, 15
  samples of 2,000 barrier phases each.

### Power state

**Battery** throughout, 72 % to 70 %, discharging (`data/a0/pmset_*.txt`).
Apple M1 Max, 8 P + 2 E cores. The counts and hashes do not depend on power;
the per-item and sync costs, and so the model, do. There is no AC baseline
for any of them.

### Not measured

- A parallel implementation of any option. Every speed above is the model's.
- A rule that stops re-evaluating a loser until its footprint is written, and
  its cost. The "evaluated once" rows only bound it.
- A0 identity on inputs other than these two, and a proof of it.
- A1 seeds 1-3 in this run (their sub-round profiles are in `data/sims/`).
- Load imbalance, and sync costs with real work between barriers.
- An AC run.
