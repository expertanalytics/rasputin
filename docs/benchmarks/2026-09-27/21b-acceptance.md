# Increment 21b acceptance run (@perf, 2026-09-27)

**Verdict: ACCEPTED.** On battery, measured back to back against the previous
merge (`071e4f3`, 21a): refine at 8 threads is 14.5-14.8 % faster on both
domains, in both pairs. The mesh sha256 is identical on both domains. No
measure regressed.

One thing to act on: the head `a03ae65` has a split phase 4-6 % slower than
21b's green commit `f707322`. The only production change between them is the
non-finite refusal in `lattice_incircle` (see "Green against head").

## Method

- **Why a re-measured base.** The Mac was on battery. The stored `bench.py`
  runs of 21a (`21a*/`) are AC. The only battery `bench.py` run,
  `serial-profile-bench/`, predates 21a. So there was no comparable baseline.
  Master `071e4f3` (the 21a merge, #102) was checked out in a temporary
  `git worktree` and measured with `--tree`, back to back with 21b's head
  `a03ae65`. The worktree was removed afterwards.
- **`--baseline` is explicit.** `071e4f3` is not an ancestor of `a03ae65`: the
  branch forked from 21a's head, not from the merge
  (`git merge-base --is-ancestor 071e4f3 a03ae65` exits 1). `find_baseline`
  would therefore skip it.
- **Two pairs**, run in the order base, 21b, base, 21b between 18:56 and
  19:05 CEST. The verdicts can be reproduced from the stored evidence:

  ```
  .venv/bin/python tools/bench.py compare docs/benchmarks/2026-09-27/21b \
      --baseline docs/benchmarks/2026-09-27/21b-base-071e4f3          # ACCEPTED
  .venv/bin/python tools/bench.py compare docs/benchmarks/2026-09-27/21b-r2 \
      --baseline docs/benchmarks/2026-09-27/21b-base-071e4f3-r2       # ACCEPTED
  .venv/bin/python tools/bench.py compare docs/benchmarks/2026-09-27/21b-base-071e4f3-r2 \
      --baseline docs/benchmarks/2026-09-27/21b-base-071e4f3          # REGRESSION (same build; see Spread)
  .venv/bin/python tools/bench.py compare docs/benchmarks/2026-09-27/21b-r2 \
      --baseline docs/benchmarks/2026-09-27/21b                       # ACCEPTED
  ```

  The two cross pairs (`21b-r2` against `21b-base-071e4f3`, and `21b` against
  `21b-base-071e4f3-r2`) are also ACCEPTED.
- **Command.** All four runs used the main tree's `tools/bench.py`, blob
  `77765b1`, the same as 21a's. Each run did its own Release build into
  `<tree>/build-bench`.

  ```
  .venv/bin/python tools/bench.py run --label <label> [--tree <worktree of 071e4f3>] \
      [--baseline docs/benchmarks/2026-09-27/<base label>] --mesh-dir <scratch> \
      --domain docs/benchmarks/2026-09-26/quarter.geojson --domain tile \
      --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 \
      --repeats 5 --threads 1,2,...,20
  ```
- **Machine and power.** Apple M1 Max (8 P + 2 E), 32 GiB, macOS 27.0. All four
  runs were on **battery**, discharging: 90 % to 86 % over the sequence. The
  `pmset -g batt` text from before and after each run is in its `run.json`.
  `caffeinate -ims` was held for all runs and for the scratch drivers.
- **Builds.** The base `_core` hashes to `654a2442…` in both of its runs. That
  is the same `.so` as 21a's acceptance run of `cade752`. 21b's hashes to
  `dfdeb6ca…` in both runs.
- **The `(dirty)` on `21b` and `21b-r2` is not a source change.** As in 21a,
  the only untracked files were the evidence directories of the runs just
  made. `git status --porcelain --untracked-files=no` was empty.
- **Phase times** (scan, split) are not in `bench.py`'s child record. They were
  measured afterwards, on battery (86 %), by a scratch driver that is not
  checked in. It uses the child's technique: import the tree's
  `build-bench/pkg`, force `threads` into `cli.refine`, and read
  `RefineOutcome.scan_seconds` and `split_seconds`. That was 5 repeats at 1, 8
  and 16 threads, with the two trees interleaved.
- **Meshes** stayed in the session scratchpad and are not kept. To regenerate
  one, rerun `bench.py` as above.

## Results (refine seconds, median of 5; ms)

Pair 1 (`21b-base-071e4f3` then `21b`); pair 2 in brackets. The default thread
count is the CLI's (`threads == 0`, 10 on this machine).

| domain | threads | base 071e4f3 | 21b | change |
|---|---:|---:|---:|---:|
| tile | 1 | 507.1 [519.6] | 488.5 [471.4] | -3.7 % [-9.3 %] |
| tile | 8 | 226.4 [227.4] | 193.0 [194.1] | -14.8 % [-14.6 %] |
| tile | 16 | 223.6 [225.8] | 197.0 [193.9] | -11.9 % [-14.1 %] |
| tile | default (10) | 226.3 [229.5] | 194.1 [193.8] | -14.2 % [-15.5 %] |
| tile | best | 222.8 @12 [225.6 @13] | 191.8 @9 [192.4 @13] | -13.9 % [-14.7 %] |
| quarter | 1 | 447.7 [471.7] | 432.8 [424.8] | -3.3 % [-9.9 %] |
| quarter | 8 | 200.7 [201.2] | 171.1 [172.0] | -14.8 % [-14.5 %] |
| quarter | 16 | 202.9 [204.0] | 173.8 [173.5] | -14.4 % [-15.0 %] |
| quarter | default (10) | 203.1 [203.5] | 171.8 [172.7] | -15.4 % [-15.1 %] |
| quarter | best | 200.4 @20 [201.2 @8] | 170.7 @14 [171.9 @12] | -14.8 % [-14.6 %] |

Every swept thread count from 1 to 20 is faster in both pairs. In pair 1 the
gain is 3.7-15.5 % on the tile and 3.3-16.1 % on the quarter circle. In pair 2
it is 9.3-15.8 % and 9.9-16.5 %. The 1-thread gain is the least stable figure:
it is -3 % to -4 % in pair 1 and -9 % to -10 % in pair 2, because the 1-thread
cell is the noisiest one on battery (see "Spread"). The phase driver's
interleaved samples put it at -6 % (tile) and -8 % (quarter).

**Phase times** (scratch driver, median of 5, min-max in brackets, ms):

| domain, threads | refine | scan | split |
|---|---:|---:|---:|
| quarter, 1 | 453.3 -> 419.1 (-8 %) | 310.7 -> 305.6 (-2 %) | 123.0 (122.2-123.8) -> 95.7 (94.7-98.1) (-22 %) |
| quarter, 8 | 198.5 -> 170.8 (-14 %) | 47.1 -> 47.4 (+1 %) | 131.9 (131.4-134.6) -> 103.8 (102.4-105.0) (-21 %) |
| quarter, 16 | 203.1 -> 174.3 (-14 %) | 45.2 -> 45.3 (0 %) | 138.2 (131.8-138.9) -> 109.4 (102.5-110.0) (-21 %) |
| tile, 1 | 501.2 -> 471.1 (-6 %) | 345.4 -> 345.5 (0 %) | 124.5 (123.4-125.4) -> 96.3 (95.0-96.5) (-23 %) |
| tile, 8 | 225.1 -> 191.7 (-15 %) | 57.7 -> 58.2 (+1 %) | 134.5 (133.3-134.7) -> 102.6 (101.6-103.4) (-24 %) |
| tile, 16 | 227.3 -> 193.9 (-15 %) | 56.1 -> 55.9 (0 %) | 137.7 (134.5-139.0) -> 106.9 (104.8-108.0) (-22 %) |

The whole gain is in the split phase: 27-32 ms at every thread count. The scan
is unchanged. The start mesh's `legalise_all` on the tile also fell from 2.6
to 0.6 ms. It goes through the same `must_flip`, so it takes the integer path
too. Rounds, insertions and flips are the same in both trees: quarter 41 /
213,464 / 445,657, tile 53 / 219,837 / 445,675.

**Ceiling** (speed-up over 1 thread; best, and at 20 threads):

| domain | base 071e4f3 | 21b | 21a acceptance (AC) | 2026-09-26 (battery) |
|---|---|---|---|---|
| tile | 2.28x @12 (2.26x @20) [2.30x, 2.28x] | 2.55x @9 (2.53x @20) [2.45x, 2.40x] | 2.28x | 2.12x @8 (2.10x @20) |
| quarter | 2.23x @20 (2.23x @20) [2.34x, 2.32x] | 2.54x @14 (2.48x @20) [2.47x, 2.46x] | 2.29x | 2.09x @15 (2.07x @20) |

The ceiling rises because the split phase, which does not scale, got shorter.
The 1-thread refine fell by 30 ms and the 8-thread refine by 30-33 ms, so the
ratio rises. The ceiling depends on the noisy 1-thread cell, so it moves by
about 0.1x between the two pairs. From 8 threads on, the curve is flat: every
median from 8 to 20 threads is within 1.4-3.2 % of that run's best. The
2026-09-26 battery figures (medians of `scaling/raw_battery.tsv`, master
`8f47e7e`) came from a different method (`scale.py`, the refine call alone)
and are context here, not a baseline.

**Quality: identical.** Both trees, both pairs:

| domain | worst angle | max degree | within tolerance | Delaunay violations | mesh sha256 (ASCII POINTS on) |
|---|---:|---:|---|---:|---|
| quarter | 0.3955 deg | 18 | yes | 0 of 641,791 | `1e531976f33b38dd883f9cf1945f6303dbac94c92f30b52b8bd7b46db7eab9ee` |
| tile | 0.6296 deg | 74 | yes | 0 of 692,056 | `11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c` |

On both domains, the mesh sha256 of 21b equals the base's, and also the 21a
acceptance run's figures. That is what a bit-identical increment requires. The
scratch driver's binary meshes also hash the same for base, green and head:
quarter `ff705683…`, tile `645919aa…`. These are the hashes @developer
reported at green.

## Spread of the median of 5 (battery)

For each tree, the two runs' medians were compared over all 42 cells (2
domains × thread counts 0-20), as |median(run 1) / median(run 2) - 1|:

| tree | median abs. difference | 90th percentile | max |
|---|---:|---:|---:|
| base 071e4f3 | 0.89 % | 2.74 % | 5.45 % (quarter, 3 threads) |
| 21b | 0.69 % | 1.77 % | 3.62 % (tile, 1 thread) |

Within a single cell, (max - min) / median of the 5 samples had a median of
3.5-4.6 % per run, against 1.1-1.9 % on AC in 21a. There were also single
samples far from the rest. One base tile sample at the default thread count
took 546 ms, against 222-230 ms for the others. One base tile 1-thread sample
took 719 ms, against 499-522 ms. The median absorbed them. Their cause has not
been measured.

Measured the other way round, base run 2 against base run 1, `bench.py`
reports REGRESSION on two quarter cells of the **same build**: +5.4 % at 1
thread and +5.8 % at 3 threads (the `compare` above). This is the second time
the median of 5 of an unchanged build has moved by about 5 % (21a on AC:
4.6 %). On battery it went over the 5 % threshold. A single-cell verdict at
about 5 % is therefore noise until that cell is rerun. 21b's gain at 8 threads
and above is about 15 % in every cell of both pairs, several times this
spread.

## Green against head

Head's split phase is 4-6 % slower than green's. The @developer quick sweep
below reported faster 21b figures than this run measured, so the green commit
`f707322` was built in a temporary worktree, since removed (the same `build()`, `.so` `2e30f2d8…`) and interleaved
with head `a03ae65` by the scratch driver. Battery, 86 %, 5 repeats, ms:

| domain, threads | refine, green -> head | split, green -> head |
|---|---:|---:|
| quarter, 1 | 414.0 -> 421.4 (+1.8 %) | 89.9 (89.5-90.9) -> 95.5 (94.4-118.7) (+6.3 %) |
| quarter, 8 | 174.4 -> 170.2 (-2.4 %) | 99.7 (97.7-106.4) -> 103.6 (102.1-121.4) (+4.0 %) |
| tile, 1 | 468.3 -> 472.5 (+0.9 %) | 91.2 (90.2-98.3) -> 96.6 (95.5-99.1) (+5.9 %) |
| tile, 8 | 185.7 -> 193.0 (+4.0 %) | 97.3 (97.0-98.1) -> 102.8 (102.7-104.4) (+5.6 %) |

The split phase went up by 4-6 ms in all four cells, and the min-max ranges
do not overlap. The refine totals are less clean: green's quarter 8-thread
samples were noisy, with a scan range of 47-61 ms. The only production-code
change from `f707322` to `a03ae65` is the `lawson.hpp` hunk that refuses
non-finite corners in `lattice_incircle` (`git diff --stat f707322 a03ae65 --
include src_cpp src_python`). So the 5 ms is attributed to that hunk by
elimination. How it costs 5 ms has not been profiled. The loss is about a
sixth of 21b's gain, and 21b still passes with a wide margin.

## Against the estimates

| source | measure | estimate or report | measured (pair 1 [pair 2]) |
|---|---|---|---|
| design, QW2 "Expected gain" | refine, every thread count | ~-18 ms: -4 % at 1 thread, -8 % at 8 | 8 threads: tile -33.4 [-33.3] ms (-14.8 %), quarter -29.6 [-29.2] ms (-14.8 %); 1 thread: -6 to -8 % (phase driver) |
| design, "What the quick wins add up to" | 8-thread refine after 21a + 21b | 0.229 -> ~0.18 s (battery profile) | 0.193 / 0.171 s (tile / quarter); 21a's AC base `227522c` was 0.249 / 0.230 s |
| @developer quick sweep (battery, vs 21a head `d3ee2ce`, at green) | quarter, 8 threads | 0.198 -> 0.164 s (-17 %) | 0.201 -> 0.171 s (-14.8 %) [0.201 -> 0.172] |
| same | tile, 8 threads | 0.225 -> 0.187 s (-17 %) | 0.226 -> 0.193 s (-14.8 %) [0.227 -> 0.194] |
| same | split phase, 8 threads | quarter 132 -> 98 ms (-26 %), tile 134 -> 97 ms (-28 %) | quarter 131.9 -> 103.8 ms (-21 %), tile 134.5 -> 102.6 ms (-24 %) |

The base figures agree with @developer's within 1.5 %. The 21b figures are
3-4 % slower than his, and his split figures match green's here (98 and 97 ms
against 99.7 and 97.3 ms). The difference is the green-to-head change above.
The design's estimate was too low by a factor of about 2 at 8 threads. As the
increment file records, the integer path is also cheaper than the filtered
path it replaces for all-node quads, not only than the exact path.

## Addendum: the head slowdown, fixed in 4be156c (@developer, 2026-09-27)

"Green against head" above measured `a03ae65`, whose non-finite refusal (four
`std::isfinite` pairs per `lattice_incircle` call) cost 4-6 % in the split
phase. `4be156c` folds that refusal into the spread test
(`!(std::fabs(x) <= bound)` refuses inf and NaN alike) and drops the explicit
checks; `[nonfinite]` still passes and UBSan stays quiet on it. Timed by
@developer with a scratch driver (not checked in), green `f707322` and head
`4be156c` alternating per repeat, 5 repeats, battery 86 %, discharging, ms,
median (min-max):

| domain, threads | refine, green -> head | split, green -> head |
|---|---:|---:|
| quarter, 1 | 410.3 (409.5-412.9) -> 409.5 (408.8-411.5) | 90.4 (89.1-90.7) -> 89.2 (89.0-90.5) |
| quarter, 8 | 163.9 (163.0-164.0) -> 163.0 (162.1-165.2) | 97.8 (97.4-98.2) -> 97.5 (96.3-99.1) |
| tile, 1 | 462.8 (460.2-466.4) -> 462.9 (460.9-465.5) | 90.6 (89.7-91.2) -> 90.4 (89.8-91.8) |
| tile, 8 | 184.4 (183.8-184.8) -> 185.8 (183.2-186.2) | 96.6 (96.1-96.9) -> 97.5 (95.9-97.6) |

Every range overlaps and every median is within 1.4 %: the gap is closed. The
bench.py verdict above was measured on `a03ae65`, so it understates 21b by
about that gap; it was not re-run on `4be156c`. Meshes equal on both domains.
