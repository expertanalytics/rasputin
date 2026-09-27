# Serial-phase profile of refine, 2026-09-27

Status: done (@perf). Measurement and analysis only: no fix, no design, no
production code changed. Every figure below is from **battery** power
(`pmset -g batt` beside each timing sweep: `data/*.pmset`). The instrumented runs
(`data/instr_t*_rounds.txt`, and `rounds_t1.md`, `chunks_t8.md` derived from
them) have no committed power record; they ran in the same battery session,
but that is not on disk. Their counts do not depend on power; their timings
do. No AC run was made.

The ask (Ola, 2026-09-27): profile refine's serial phase before anyone designs
a fix. "Some of the serial parts could be parallellized by multicoloring/dd techniques, but let's wait for the analysis before we get ahead of ourselves."

## Findings

1. **About a third of single-thread refine does not parallelise, not half.**
   The serial part is 34 % (0.159 s of 0.462 s): the split + flip phase (27 %)
   plus the between-round work outside both timers (6 %). Amdahl's law with a
   perfectly parallel scan gives a ceiling of 2.9x. Adding the scan's measured
   scaling predicts 2.10x at 8 threads, and 2.02x was measured. The rest of
   the gap is measured too: the split phase runs about 8 % slower after a
   multi-threaded scan.
2. **Inside the serial phase, Lawson legalisation costs more than insertion.**
   The flip test (`must_flip` plus predicates) is 13.2 % of refine. Of that,
   5.0 % is the exact incircle fallback. It fires on 9.1 % of incircle tests,
   and every one of those tests was exactly cocircular. The legalise_around
   loop is 6.3 %, and topology writes (split, flip, put and repoint) are 7.3 %.
   Allocation is 1.5 % and locking is 0 %: the serial phase takes no lock. The
   `active` sort between rounds is 3.6 %.
3. **The scan speeds up about 5x at 8 threads (4.95x), not 8x, and flattens
   from about 8.** Load imbalance across the contiguous chunks costs 16.5 ms of
   the 61.8 ms. Round 1 alone is 6.1 ms of it, about 39 % of that run's 15.8 ms (`data/chunks_t8.md`). Thread start
   costs 3.8 ms, and each worker is about 9 % slower than a lone thread. Above
   10 threads (8 P + 2 E cores) the last chunk starts about 30 ms late. The
   bytes read (2.6 GB/s at 8 threads) show no sign of a memory-bandwidth limit.
   That is inferred from byte counts; no hardware counters were read.
4. **Insertions in the same round are dense and overlap.** 41 rounds; 77 % of
   the 213,464 insertions fall in rounds 9-19, at 9,000-18,700 per round. In
   rounds 5-21 each round touches 69-93 % of the 382 64×64-node blocks that ever receive an insertion (about a fifth of the 64×64-node blocks that intersect `quarter.geojson`: 1,774-1,857 depending on the raster origin, counted by @reviewer with tifffile and shapely; that count is not committed. Refine never inserts into the rest; why, flat ground or NoData, is not checked).
   In rounds 11-29 the median distance to the nearest same-round insertion
   is 3 nodes. An insertion writes 5.5 triangle slots on average (p99 9) and
   reads or writes 11.0 (p99 18). 95 % of insertions share a slot, read or
   written, with another insertion of the same round; 65 % share a written
   slot. The serial loop already defers 51 % of marked triangles because an
   earlier insertion in the same round rewrote their slot.

The 2026-09-26 README said the serial part was "thought to be the serial
insert and flip phase" and "roughly half". The profile confirms the phase and
corrects the share: 34 % at one thread. The low ceiling comes from that 34 %
together with a scan that stops at about 5x.

## Method

- **Commit** `d6d6beb` (branch `serial-profile`), clean apart from this
  directory. The 1 m benchmark: DEM
  `tests/fixtures/dem_archive/7908_3_10m_z33.tif` (5051 × 5051 nodes, 10 m),
  domain `docs/benchmarks/2026-09-26/quarter.geojson`, tolerance 1, all other
  CLI defaults. Every run gave the same output: 41 rounds, 213,464 inserted,
  445,657 flips, 428,217 triangles, and the same mesh sha256 for the plain,
  profiled and instrumented builds (`ff705683…` over the binary VTK).
- **Machine**: Apple M1 Max, 8 P + 2 E cores, 32 GiB, macOS 27.0, AppleClang
  21.0.0, Python 3.14.7. **Power: battery**, 83-86 %, throughout
  (`data/*.pmset`, and `run.json` of the bench.py run). `powermode 0` was
  observed but is not recorded in any committed file.
- **Builds** (all scratch, in gitignored `build-*` directories):
  - `build-prof`: `CMAKE_BUILD_TYPE=Release`, `CMAKE_CXX_FLAGS=-g`, which
    gives `-g -O3 -DNDEBUG … -flto`. Used for the phase sweeps. pybind11 strips
    Release modules, so it has no symbols.
  - `build-prof-lto`: `RelWithDebInfo` with
    `CMAKE_CXX_FLAGS_RELWITHDEBINFO="-O3 -g -DNDEBUG -flto"` and
    `CMAKE_MODULE_LINKER_FLAGS="-flto -Wl,-object_path_lto,<dir>/lto.o"`, then
    `dsymutil`. It is unstripped and has line tables, and runs as fast as the
    Release bench build: 0.456 s against 0.460 s, and 0.455 s against 0.462 s,
    in back-to-back pairs of 5 calls. Used for the profile. A first attempt
    without `-flto` was 16 % slower in the scan, so it was not used.
  - `build-instr`: Release (as bench.py builds it) in a scratch worktree of
    `d6d6beb` with `scripts/instrument.patch` applied. That is a **local patch,
    never committed**, generated by `scripts/instrument.py`. It adds counters
    and per-round records only. Used for questions 2 and 3. Its split-phase
    timings are inflated by its own logging and are not reported.
- **Timing harness**: `scripts/prof_driver.py`, bench.py's child technique
  (the build's `pkg/` first on `sys.path`, and `cli.refine` wrapped with
  threads forced). It calls refine K times with the same arguments and prints
  `RefineOutcome`'s phase seconds (increment 17) per call. `rest` is the
  wrapper's wall time minus the four phase timers: setup, the between-round
  `active` rebuild, output and the pybind return.
  `scripts/sweep.sh` makes 3 interleaved passes over thread counts
  1-8, 10, 12, 16, 20, with 5 calls per process: 15 samples, reported as the
  median (`scripts/summ.py`).
- **Profiler**: `/usr/bin/sample` at 1 ms, attached for 18 s to a driver that
  was making 45 single-thread calls (`xctrace` is not installed: Command Line
  Tools only). `scripts/attribute.py` computes each call-graph node's self
  count and resolves every `_core` address with `atos -i` against the dSYM to
  its inline chain. It files each sample under the innermost frame that
  matches a category. The matching ignores template arguments, which would
  otherwise match on type names. Sampling slowed the calls by about 10 %
  (0.510 s, `data/sample_t1_driver_phases.txt`). The profile's scan share
  (66.0 %) matches the timer's (64 %).
- **bench.py cross-check**: `tools/bench.py run --label serial-profile-bench
  --domain docs/benchmarks/2026-09-26/quarter.geojson --threads
  1,2,4,6,8,10,12,16,20 --repeats 5`. Its evidence is in
  `../serial-profile-bench/`. It is the first stored battery run for the
  quarter at `d6d6beb`. It records `dirty` only because this directory was
  untracked at the time.

## 1. Where single-thread refine time goes

### Phase timers (RefineOutcome, increment 17), battery

Sweep 2 (`data/phases_sweep2_battery.txt`), 15 samples per row, medians:

| threads | refine s | speed-up | scan s | scan speed-up | split s | rest s | serial share |
|---|---|---|---|---|---|---|---|
| 1 | 0.462 | 1.00x | 0.303 | 1.00x | 0.127 | 0.029 | 34 % |
| 2 | 0.332 | 1.39x | 0.174 | 1.74x | 0.126 | 0.029 | 47 % |
| 3 | 0.283 | 1.63x | 0.124 | 2.44x | 0.128 | 0.029 | 56 % |
| 4 | 0.266 | 1.74x | 0.100 | 3.03x | 0.136 | 0.030 | 63 % |
| 5 | 0.251 | 1.84x | 0.083 | 3.68x | 0.136 | 0.030 | 67 % |
| 6 | 0.238 | 1.94x | 0.071 | 4.30x | 0.136 | 0.030 | 70 % |
| 7 | 0.232 | 1.99x | 0.065 | 4.66x | 0.136 | 0.029 | 72 % |
| 8 | 0.229 | 2.02x | 0.061 | 5.00x | 0.137 | 0.030 | 74 % |
| 10 | 0.236 | 1.96x | 0.063 | 4.82x | 0.142 | 0.030 | 73 % |
| 12 | 0.235 | 1.97x | 0.062 | 4.90x | 0.142 | 0.030 | 74 % |
| 16 | 0.227 | 2.03x | 0.054 | 5.63x | 0.143 | 0.030 | 76 % |
| 20 | 0.230 | 2.00x | 0.057 | 5.33x | 0.142 | 0.030 | 75 % |

`legalise` (start mesh) is under 0.1 ms and `quality` about 1 ms at every
count. "Serial share" is (refine - scan) / refine.

Sweep 1 (`data/phases_sweep1_battery.txt`) ran 20 minutes earlier, also on
battery, with the same build. Every 1-thread figure was 14-15 % slower
(refine 0.527 s, scan 0.348 s, split 0.146 s), and the serial share was the
same, 34 %. Its ceiling was 2.25x at 20. The cause of the session-to-session
difference was not found. The bench.py run, made between the two sweeps,
agrees with sweep 2: 0.466 s at 1 thread, best 2.02x at 16. So sweep 2 is the
one reported. Compare figures only within one sweep.

### Amdahl

- Serial part at 1 thread: 0.462 - 0.303 = **0.159 s, 34 %**. With a scan
  of zero cost the ceiling would be 0.462 / 0.159 = **2.9x**.
- With the scan's measured 8-thread time (0.061 s) and the 1-thread serial
  part: 0.159 + 0.061 = 0.220 s, a predicted **2.10x**. Measured: 2.02x
  (0.229 s). The difference is measured: split takes 0.127 s at 1-3 threads
  and 0.136-0.143 s at 4 or more. Why is not measured; it is thought to be
  cache locality, since the scan's results and the mesh lines were last
  touched by other cores.
- The measured ceiling today is **2.03x at 16 threads** (bench.py: 2.02x at
  16), against the 2026-09-26 figure of about 2.2x. That older run was at
  8f47e7e in a different harness and is not a baseline for this one
  (`docs/increments/README.md`).

### Profile of the single-thread call, battery

Taken from `data/attribution_t1.txt` (from `data/sample_t1.txt`): 13,079 self
samples inside `refine`, at 1 ms each. The shares are of the whole refine
call.

| where | samples | share of refine | phase |
|---|---|---|---|
| scan (`scan<>`, row spans, row segments, the chunk lambda) | 8,634 | 66.0 % | parallel |
| legalise_around loop (stack push/pop, neighbour lookups, `touched` marks) | 825 | 6.3 % | split |
| `must_flip` itself (edge lookup in the neighbour, frame points) | 718 | 5.5 % | split |
| exact incircle fallback (`incircleadapt`, expansions) | 653 | 5.0 % | split |
| filtered orient2d/incircle | 349 | 2.7 % | split |
| `LatticeMesh::put` / `repoint` (shared by split and flip) | 339 | 2.6 % | split |
| `LatticeMesh::flip` | 330 | 2.5 % | split |
| `split_inside` / `split_edge` / `add_vertex` | 284 | 2.2 % | split |
| refine body without line info (line 0 in the line table) | 208 | 1.6 % | split, mostly (not resolvable) |
| allocation (malloc/free/memmove) | 198 | 1.5 % | split, mostly |
| split loop bookkeeping (result read, `touched`, `skipped`) | 10 | 0.1 % | split |
| rebuild `active`: collect, `std::sort`, `unique` | 466 | 3.6 % | between rounds (in `rest`) |
| output (vertices, z, triangles, constraint edges) | 39 | 0.3 % | after the loop (in `rest`) |
| pybind conversion inside the call | 18 | 0.1 % | |
| start quality, setup (`to_lattice`, `LatticeMesh::build`) | 8 | 0.1 % | |

Grouped over the split phase, which is about 30 % of refine:

- **Lawson flips dominate insertion.** Flip test 13.2 % (must_flip 5.5, exact
  5.0, filtered 2.7). The legalise_around loop is 6.3 % and flip writes
  2.5 %. Insertion (`split_*`) is 2.2 %. `put`/`repoint` (2.6 %) serves both.
- **Exact predicates on the lattice.** Instrumented counts
  (`data/instr_t1_rounds.txt`, the `S` records) for the split phase: 1,615,895
  incircle tests, of which **146,962 (9.1 %) took the exact path, and all
  146,962 returned Cocircular**. orient2d: 3,232,783 tests, 527 exact. So 9 %
  of the incircle tests cost 38 % of flip-test time (653 of 1,720 samples),
  and every one of them was an exact tie. That fits cocircular DEM nodes (a
  grid rectangle's four corners), but the quads themselves were not
  classified, so that remains an inference.
- **Allocation**: 1.5 % in total. Charged to its owner, 0.7 % is
  under legalise_around, whose only allocation is its per-call `std::vector` stack, and 0.4 % sits under
  must_flip; the remainder is under 0.2 %.
- **Locking**: no sample in any mutex, lock or `psynch` frame. The serial
  phase takes no lock, and at one thread `for_each_chunk` runs inline with no
  thread.
- **Rescans**: the rescan is the next round's scan of every touched slot, so
  it is inside the scan share. Over the run, the scan visits 39,502,901 DEM
  nodes. Round 1 visits 7,078,976, about the whole domain. So each domain node
  is scanned **5.6 times** on average (`data/chunks_t8.md`, nodes column).

## 2. The parallel phase on its own

Scan only, from the instrumented build's per-chunk records
(`data/chunks_sweep.tsv`). Medians of 3 runs, battery. Each figure is summed
over the 41 rounds.

| threads | scan ms | scan speed-up | ideal ms (1-thread / N) | busy thread-ms | sum of mean chunk ms | sum of max chunk ms | imbalance ms (max - mean) | last chunk start ms | join tail ms |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 306.0 | 1.00x | 306.0 | 306 | 306.0 | 306.0 | 0.0 | 0.0 | 0.04 |
| 2 | 177.1 | 1.73x | 153.0 | 315 | 157.3 | 174.7 | 17.4 | 1.6 | 0.82 |
| 4 | 101.2 | 3.02x | 76.5 | 327 | 81.7 | 98.4 | 16.7 | 2.3 | 0.84 |
| 6 | 72.7 | 4.21x | 51.0 | 334 | 55.7 | 69.4 | 13.7 | 3.1 | 0.99 |
| 8 | 61.8 | 4.95x | 38.2 | 334 | 41.8 | 58.3 | 16.5 | 3.8 | 1.05 |
| 10 | 64.2 | 4.77x | 30.6 | 387 | 38.7 | 58.8 | 20.1 | 5.0 | 1.32 |
| 12 | 62.2 | 4.92x | 25.5 | 377 | 31.4 | 53.3 | 21.9 | 29.9 | 1.48 |
| 16 | 54.8 | 5.58x | 19.1 | 372 | 23.2 | 43.7 | 20.5 | 28.8 | 1.69 |
| 20 | 57.8 | 5.29x | 15.3 | 382 | 19.1 | 44.6 | 25.5 | 36.9 | 2.01 |

At 8 threads, the 61.8 ms breaks down into the ideal 38.2 ms plus four
measured losses:

- **Per-thread slowdown: +3.6 ms.** Busy thread-time is 334 ms against
  306 ms alone, so each worker runs 9 % slower. The cause is not measured. It
  is thought to be P-cluster clock or shared-L2 effects.
- **Load imbalance: +16.5 ms, the largest loss.** `for_each_chunk` splits the
  sorted `active` list into equal counts, not equal work. In the single run of `data/chunks_t8.md`, round 1 alone
  accounts for 6.1 ms of it: its slowest chunk took 9.7 ms against a mean of
  3.6 ms (`data/chunks_t8.md`). The imbalance is already 17 ms at 2 threads.
- **Thread start: about 1.8 ms** until the first chunk runs, and 3.8 ms until
  the last one starts (41 spawns of 8 `jthread`s; about 90 µs per round).
- **Join: 1.0 ms.**

Rounds 30-41 each scan fewer than 1,200 triangles, in 0.07-0.17 ms, most of
which is spawn cost.

**Is the scan flat after about 7 threads?** Nearly: 4.95x at 8, 4.8-4.9x at
10-12, 5.3-5.6x at 16-20. From 10 threads on, the busy time rises to
372-387 ms. The two E cores are thought to account for this, but no core
affinity was recorded. From 12 threads on (more threads than cores), the last
chunk starts about 30 ms late. Refine as a whole is flat from about 7 because
the serial part is 70-76 % of the multi-thread time.

**Memory bandwidth**: the scan reads 39.5 M float32 nodes, 158 MB, in
61.8 ms at 8 threads. That is 2.6 GB/s, far below this machine's memory
bandwidth. There is no sign of a bandwidth limit. This is an inference from
byte counts; no hardware counters were read.

## 3. Data bearing on multicolouring or domain decomposition

This section designs neither. The data comes from the local patch
(`scripts/instrument.patch`) at 1 thread. `data/rounds_t1.md` holds all 41
rounds. Selected rounds are shown here, all battery.

Definitions, from the code as instrumented:

- **Write set**: the slots an insertion writes. That is the split triangle
  `t`, the appended slots, the edge neighbour `u` if any, and every slot
  legalise_around's `on_write` reports.
- **Footprint**: the write set plus every slot legalise_around reads, meaning
  each popped triangle and its neighbour across the tested edge. The slots
  `repoint` rewrites are all among these.
- **Conflict**: two insertions of the same round share a slot. The
  footprints are measured in the serial order, so a later insertion's
  footprint is taken on the mesh as earlier ones left it.
- **Nearest-neighbour distance**: between the insertion points of one round,
  in DEM nodes (1 node = 10 m).
- **Blocks**: a 64×64-node partition used only to describe spread. The
  domain's insertions touch 382 of these blocks over the whole run.

| round | active | marked | inserted | deferred (slot touched) | deferred (edge) | flips/ins | footprint mean (p99, max) | write set mean | ins. in a footprint conflict | ins. in a write-write conflict | mean/max conflict degree | median NN dist (nodes) | blocks hit | max ins/block |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1,822 | 480 | 231 | 244 | 5 | 2.18 | 9.9 (15, 16) | 5.2 | 100 % | 96 % | 6.5/15 | 26.9 | 116 | 5 |
| 5 | 5,437 | 4,997 | 2,063 | 2,717 | 217 | 2.42 | 11.0 (18, 24) | 5.5 | 100 % | 89 % | 7.7/21 | 7.8 | 264 | 50 |
| 10 | 44,913 | 29,273 | 11,769 | 15,412 | 2,092 | 2.23 | 11.1 (18, 28) | 5.5 | 99 % | 76 % | 6.1/22 | 3.6 | 344 | 169 |
| 14 | 89,510 | 45,306 | 18,714 | 23,137 | 3,455 | 2.08 | 11.0 (18, 36) | 5.5 | 97 % | 66 % | 4.7/23 | 3.0 | 354 | 237 |
| 20 | 56,930 | 19,662 | 8,675 | 9,570 | 1,417 | 1.89 | 10.8 (18, 26) | 5.4 | 89 % | 50 % | 3.0/17 | 3.0 | 303 | 189 |
| 25 | 11,418 | 3,036 | 1,430 | 1,437 | 169 | 1.85 | 10.7 (18, 28) | 5.4 | 78 % | 38 % | 2.0/13 | 3.2 | 173 | 51 |
| 30 | 1,103 | 217 | 108 | 93 | 16 | 1.95 | 10.8 (22, 22) | 5.4 | 62 % | 27 % | 1.3/7 | 4.1 | 31 | 21 |
| 35 | 59 | 15 | 9 | 5 | 1 | 2.33 | 10.9 (14, 14) | 5.4 | 56 % | 44 % | 1.1/3 | 8.6 | 3 | 7 |

Conflict degree is per footprint conflict: the number of other insertions of
the round an insertion shares a footprint slot with.

Over all rounds:

- **Insertions per round**: 231 in round 1, peaking at 18,714 in round 14,
  then one each in rounds 39-40. Rounds 9-19 hold 164,608 of the 213,464
  insertions (77 %).
- **Spread**: in rounds 5-25 each round touches 173-355 of the 382 blocks
  (45-93 %; 69-93 % in rounds 5-21). In the big rounds a block receives at
  most 154-248 insertions. In rounds 11-29 the median nearest-neighbour
  distance is about 3 nodes (30 m; 2.8-3.6 over those rounds, `data/rounds_t1.md`); in round 1 it is 27 nodes. Insertions are
  spread over all of the domain that refine works on (the 382 blocks above) at once, and they are close together.
- **Footprint size**: write set mean 5.49 slots (p50 5, p90 7, p99 9,
  max 18). Footprint mean 10.96 (p50 10, p90 14, p99 18, max 36). Flips per
  insertion 2.09 (445,657 / 213,464).
- **Overlap within a round**: 202,927 insertions (95 %) share a footprint
  slot with at least one other insertion of the round, in 495,188 pairs.
  138,637 (65 %) share a written slot (111,952 pairs). 76 % overlap an
  insertion earlier in the serial order.
- **The loop's own deferral**: of 509,029 marked (non-converged) triangle
  occurrences, **259,051 (51 %)** were skipped because an earlier insertion
  of the same round had already rewritten their slot. 36,514 (7 %) were
  edge-split deferrals (the neighbour was touched). 213,464 (42 %) were
  inserted. Every deferred triangle is rescanned the next round, which adds
  to the 5.6× scan amplification in section 1.

Not measured here: the geometric extent of a footprint (vertex positions) and
how many footprints would cross a partition boundary of a given size. The
patch records slots, not coordinates.

## What is measured and what is inferred

Measured: every timing and count in the tables; the profile attribution
(which rests on sample's 1 ms sampling and on `atos` line tables, where 1.6 %
had no line); the predicate counts; the chunk timings; the footprint and
conflict counts; and that the split phase takes no lock.

Inferred, and worded that way above:

- that the exact-path incircle ties are cocircular lattice quads;
- the cause of the split phase's slowdown after a parallel scan;
- the cause of the 9 % per-worker slowdown;
- that E cores cause the busy-time rise from 10 threads;
- that there is no bandwidth limit (from byte counts, not counters);
- the cause of sweep 1's 14-15 % slower single-thread figures.

## Regenerating

The scripts are dated one-off evidence, kept verbatim, like `2026-09-26/`.
They hard-code this session's scratchpad path (`S=` / `R=` at the top of
`sweep.sh`; the driver takes absolute `--pkg` and DEM paths). Fix those
paths before rerunning.

1. Phase sweep: configure `build-prof` as above, then assemble
   `build-prof/pkg/tin_engine/` the way bench.py's `build()` does (symlinks
   to `src_python/tin_engine/*`, plus a copy of the `.so`). Then run
   `scripts/sweep.sh <abs pkg> <out> --domain <abs quarter.geojson>` and
   `scripts/summ.py <out>`.
2. Profile: configure `build-prof-lto` as above, run `dsymutil` and assemble
   the pkg. Start `prof_driver.py --threads 1 --repeat 45 --pause 3 -- mesh …`,
   attach `sample <pid> 18 1 -mayDie -file s.txt` once it prints `PID`, then
   run `scripts/attribute.py s.txt <dSYM>`.
3. Instrumented: `git worktree add --detach <dir> d6d6beb`, then
   `python scripts/instrument.py <dir>` (or `git apply
   scripts/instrument.patch`). Build Release in `<dir>/build-instr`, run the
   driver with `RASPUTIN_PROF_OUT=<file>`, then `scripts/analyse.py <file>`
   (rounds) or `--chunks` (per-chunk scan). `scripts/chunksum.py` sums the
   per-chunk records.

Not in the repository:

- the raw instrumented logs with per-insertion `I` records (about 32 MB per
  run);
- the meshes;
- the builds and dSYMs.

All of these were in the session scratchpad
(`/private/tmp/claude-501/…/scratchpad`), which does not survive the session.
Steps 1-3 regenerate them. The `R`, `C` and `S` records of the 1- and 8-thread
instrumented runs are kept in `data/instr_t*_rounds.txt`.
