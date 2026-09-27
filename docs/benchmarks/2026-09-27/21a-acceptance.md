# Increment 21a acceptance run (@perf, 2026-09-27)

**Verdict: ACCEPTED.** On AC, against the previous merge (`227522c`) measured
back to back: refine at 8 threads is 11 % faster on the tile and 15 % faster on
the quarter circle. The mesh sha256 is identical on both domains. No measure
regressed.

## Method

- **Why a re-measured base.** The only stored `bench.py` run
  (`serial-profile-bench/`) is on battery, so there is no AC baseline. As
  `docs/increments/README.md` ("Acceptance") says to, master `227522c` (the
  previous merge, #101) was checked out in a temporary `git worktree` and
  measured with `--tree`, back to back with 21a's HEAD `cade752`. The worktree
  was removed afterwards.
- **Two pairs.** Each tree was measured twice, in the order base, 21a, base,
  21a, between 12:54 and 13:05 CEST. The second pair (`-r2`) is there to
  measure how much the median of 5 moves between runs (see "Spread" below).
  `bench.py`'s verdicts: `21a` against `21a-base-227522c` (`--baseline`),
  `21a-r2` against `21a-base-227522c-r2` (`--baseline`), and
  `21a-base-227522c-r2` against `21a-base-227522c` (found by the baseline
  search). All three are ACCEPTED.
- **Command** (the main tree's `tools/bench.py`, blob `77765b1`, for all four
  runs; each run did its own Release build into `<tree>/build-bench`):

  ```
  .venv/bin/python tools/bench.py run --label <label> [--tree <worktree of 227522c>] \
      [--baseline docs/benchmarks/2026-09-27/<base label>] --mesh-dir <scratch> \
      --domain docs/benchmarks/2026-09-26/quarter.geojson --domain tile \
      --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 \
      --repeats 5 --threads 1,2,...,20
  ```
- **Machine and power.** Apple M1 Max (8 P + 2 E), 32 GiB, macOS 27.0, on
  **AC** for all four runs (96-98 %, charging; `pmset -g batt` before and after
  each run is in its `run.json`). `caffeinate -ims` was held for the whole
  sequence.
- **The `(dirty)` on `21a` and `21a-r2` is not a source change.** `bench.py`
  sets `dirty` from `git status --porcelain`, which also lists untracked files.
  During the 21a runs the only untracked files were the base runs' evidence
  directories, written minutes earlier. `git status --porcelain
  --untracked-files=no` was empty. The `_core` sha256 is the same in `21a` and
  `21a-r2` (`654a2442…`), and so is the base's in both of its runs
  (`c1d29205…`). Any back-to-back pair where the new tree is the repo will show
  this until the first run's evidence is committed.
- **Scan time** is not in `bench.py`'s child record. It was measured afterwards
  on AC by a scratch driver that uses the child's technique: it imports each
  tree's `build-bench/pkg`, forces `threads` into `cli.refine`, and reads
  `RefineOutcome.scan_seconds`. That was 5 repeats at 1, 8 and 16 threads, with
  the two trees interleaved. The driver is not checked in. Its figures are the
  "scan" rows below.
- **Meshes** stayed in the session scratchpad and are not kept. To regenerate
  one, rerun `bench.py` as above. The quality run is one `--ascii` child per
  domain at the CLI's default thread count.

## Results (refine seconds, median of 5; ms)

Pair 1 (`21a-base-227522c` then `21a`); pair 2 in brackets.

| domain | threads | base 227522c | 21a | change |
|---|---:|---:|---:|---:|
| tile | 1 | 520.8 [518.0] | 496.9 [500.8] | -4.6 % [-3.3 %] |
| tile | 8 | 249.5 [247.8] | 222.4 [222.1] | -10.9 % [-10.4 %] |
| tile | 16 | 245.8 [245.3] | 220.1 [220.1] | -10.5 % [-10.3 %] |
| tile | default (10) | 254.5 [254.5] | 219.4 [219.0] | -13.8 % [-13.9 %] |
| tile | best | 245.8 @16 [245.3 @16] | 218.4 @10 [218.4 @12] | -11.1 % [-11.0 %] |
| quarter | 1 | 464.0 [486.5] | 444.8 [445.0] | -4.1 % [-8.5 %] |
| quarter | 8 | 230.4 [233.9] | 196.5 [197.0] | -14.7 % [-15.8 %] |
| quarter | 16 | 224.4 [227.0] | 194.4 [195.4] | -13.3 % [-13.9 %] |
| quarter | default (10) | 233.5 [232.8] | 196.0 [194.3] | -16.1 % [-16.5 %] |
| quarter | best | 224.4 @16 [227.0 @16] | 194.2 @13 [194.0 @15] | -13.5 % [-14.5 %] |

Every swept thread count, 1 to 20, is faster in pair 1: by 4.6-14.7 % on the
tile and 4.1-16.3 % on the quarter circle.

**Scan at 8 threads** (scratch driver, median of 5, min-max in brackets):
tile 65.3 (64.8-68.1) -> 57.4 (57.0-59.2) ms; quarter 59.9 (59.3-60.5) ->
46.6 (46.5-47.4) ms. At 1 thread the scan is unchanged: tile 342.0 -> 342.8,
quarter 303.7 -> 303.1 ms.

**Ceiling** (speed-up over 1 thread; best, and at 20 threads):

| domain | base 227522c | 21a | 2026-09-26 (AC) |
|---|---|---|---|
| tile | 2.12x @16 (2.09x @20) [2.11x, 2.07x] | 2.28x @10 (2.26x @20) [2.29x, 2.28x] | 2.2x |
| quarter | 2.07x @16 (2.04x @20) [2.14x, 2.11x] | 2.29x @13 (2.27x @20) [2.29x, 2.28x] | 2.2x |

In 21a the curve is flat from 9 threads on: every median from 9 to 20 threads is within 1.3 % of the best. The base is
not monotone: on the tile it rises from 249.5 ms at 8 threads to 257.7 at 9,
then falls slowly to its best at 16. 21a's ceiling is above the 2026-09-26 AC
figure of 2.2x. The base, re-measured today with `bench.py`, is below it
(2.07-2.14x). The 2026-09-26 figure came from a different method (3 repeats,
timings parsed from stdout), so it is context here, not a baseline.

**Quality: identical.** Both trees, both pairs:

| domain | worst angle | max degree | within tolerance | Delaunay violations | mesh sha256 |
|---|---:|---:|---|---:|---|
| quarter | 0.3955 deg | 18 | yes | 0 of 641,791 | `1e531976f33b38dd883f9cf1945f6303dbac94c92f30b52b8bd7b46db7eab9ee` |
| tile | 0.6296 deg | 74 | yes | 0 of 692,056 | `11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c` |

The mesh sha256 of 21a equals the base's on both domains, as the design
requires for a bit-identical increment.

## Spread of the median of 5 (Ola's ruling on the 5 % threshold)

For each tree, the two runs' medians were compared over all 42 cells (2
domains × thread counts 0-20), as |median(run 1) / median(run 2) - 1|:

| tree | median abs. difference | 90th percentile | max |
|---|---:|---:|---:|
| base 227522c | 0.66 % | 1.65 % | 4.62 % (quarter, 1 thread) |
| 21a | 0.24 % | 0.72 % | 0.85 % |

Within a cell, (max - min) / median of the 5 samples had a median of 1.1-1.9 %
per run, and a maximum of 3.8-5.2 % for 21a and 8.5-13.9 % for the base.

The 4.62 % was not one outlier absorbed by the median. In `21a-base-227522c-r2`
all five quarter 1-thread samples were high (0.482-0.550 s, against
0.460-0.467 s in the first run; 4.8 % measured the other way round, against
run 1). The rest of that run's quarter cells were 0-3.4 % slower than run 1,
mostly 1-2 %. The cause has not been measured. So on AC the median of 5 usually moves
by under 2 % between back-to-back runs of the same build, but it has been seen
to move by 4.6 %, just under the 5 % threshold. A single-cell regression of
about 5 % is therefore not conclusive from one run: rerun that cell before
calling it. 21a's gains, 10-16 % at 8 threads and above, are several times that
spread, and they reproduced in both pairs.

## Against the estimates

| source | measure | estimate | measured (pair 1 [pair 2]) |
|---|---|---|---|
| design, "What the quick wins add up to" | 8-thread refine, all quick wins (21a + 21b) | 0.229 -> ~0.18 s | 21a alone: tile 0.249 -> 0.222 s, quarter 0.230 -> 0.197 s |
| design, QW1 + QW3 + stack | 21a alone at 8 threads | ~-10 % (QW1 -6 %, QW3 ~14 ms, stack 0.7 %) | tile -10.9 % [-10.4 %], quarter -14.7 % [-15.8 %] |
| design, QW1 | scan at 8 threads | ~62 -> 47 ms | tile 65.3 -> 57.4 ms, quarter 59.9 -> 46.6 ms |
| design, QW3 + stack | 1-thread refine | ~-17 ms (QW3 14 ms, stack ~3 ms) | tile -24 [-17] ms, quarter -19 [-42] ms |
| @developer quick sweep (vs red `7b54a4e`) | tile, 8 threads: scan / refine | 65.0 -> 57.7 ms / 243 -> 219 ms | 65.3 -> 57.4 ms / 249.5 -> 222.4 ms |

The design's 0.229 s baseline is the battery profile's figure, so it is
compared here by ratio, not by absolute value. 21a alone meets or beats its
share of the estimate. On the tile, the scan gain (-7.9 ms at 8 threads) is
about 20 ms less than the refine gain (-27 ms). That difference is thought to
be the QW3 merge and the caller-owned stack, as it matches the 1-thread gain,
where the scan did not change. No profile has confirmed that split. The
quarter's 1-thread pair-2 figure (-42 ms) is inflated by the base's high cell
described above.
