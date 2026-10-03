# Increment 23b acceptance: frozen edges and the seam pass (@perf, 2026-10-04)

Branch `worktree-agent-a8b29c323dc2ad69d` at `c4fb2bf` (`@reviewer` code
round 2 APPROVED), against its merge base `b4bcdc3` (increment 24's merge),
back to back, AC power throughout. The rule is
`docs/increments/23-basin-scale.md`, "@perf acceptance", 23b: the README
"Acceptance" rule in full, the 1 m benchmark's mesh hash unchanged (mask 0),
and refine within noise at every thread count.

## Verdict

**REGRESSION: refine_s at 1 thread, +4.7 % (tile) and +4.3 % (quarter);
+1 to +4 % at 2 to 20 threads.** The meshes are identical.

1. **Mesh hash unchanged.** All twelve runs give the same tile and quarter
   mesh sha256, the same as increment 24's acceptance, and the same quality.
   Every (domain, thread count) cell has a single (max_error, rounds,
   inserted, flips) tuple over all twelve runs.
2. **Refine is not within noise.** 23b is slower than the base in every
   run at every low thread count. At 1 thread the runs' medians do not
   overlap: tile base 0.5124-0.5232 s against 23b 0.5414-0.5517 s, quarter
   0.4532-0.4648 s against 0.4718-0.4844 s. The pooled change (20 to 25
   samples a side) is under `bench.py`'s 5 % threshold everywhere, but
   `bench.py` itself flagged tile at 1 thread in two of the four back-to-back
   pairs (+5.7 %, +5.5 %). The spread between runs of one binary is about 2 %.
3. **The cause, measured.** 23b adds an `on_frozen` test per DEM node
   inside `scan`'s inner loop (`include/terrain/refinement/scan.hpp`, the
   `nodes` branch). With mask 0 it is always false, but it is still tested
   per node. Experiment E is `c4fb2bf` plus `23b-acceptance/exp-unswitch.diff`
   (12 lines): when the triangle has no frozen edge, that loop runs as in
   the base. E is within -2.4 % to +1.6 % of the base at every cell, and
   +0.5 % / +0.7 % at 1 thread. Its meshes are identical too. The fix is
   for `@developer`. E was not committed to the branch.

## Results: refine_s, pooled median (seconds) and change against the base

B is the base `b4bcdc3` (4 clean runs), N is 23b `c4fb2bf` (5 runs), and E is
the experiment (2 runs). The full table, every thread count with run ranges,
is `23b-acceptance/tables-all.md`.

| domain | threads | B | N | E | N/B % | E/B % |
|---|---:|---:|---:|---:|---:|---:|
| tile | default | 0.2002 | 0.2032 | 0.2003 | +1.5 | +0.0 |
| tile | 1 | 0.5201 | 0.5447 | 0.5226 | +4.7 | +0.5 |
| tile | 2 | 0.3377 | 0.3507 | 0.3405 | +3.9 | +0.8 |
| tile | 4 | 0.2489 | 0.2560 | 0.2500 | +2.8 | +0.4 |
| tile | 8 | 0.2015 | 0.2067 | 0.2026 | +2.6 | +0.5 |
| tile | 20 | 0.2000 | 0.2032 | 0.2002 | +1.6 | +0.1 |
| tile | 2-20, range | | | | +0.6..+3.9 | -2.4..+1.1 |
| quarter | default | 0.1795 | 0.1811 | 0.1795 | +0.9 | -0.0 |
| quarter | 1 | 0.4586 | 0.4782 | 0.4617 | +4.3 | +0.7 |
| quarter | 2 | 0.2967 | 0.3053 | 0.2979 | +2.9 | +0.4 |
| quarter | 4 | 0.2191 | 0.2242 | 0.2194 | +2.3 | +0.1 |
| quarter | 8 | 0.1779 | 0.1810 | 0.1801 | +1.7 | +1.2 |
| quarter | 20 | 0.1797 | 0.1823 | 0.1781 | +1.4 | -0.9 |
| quarter | 2-20, range | | | | +0.3..+2.9 | -1.4..+1.6 |

**Scaling ceiling** (pooled, 1 thread over 20), AC: B tile 2.60x and
quarter 2.55x; N 2.68x and 2.62x (a slower single thread, not better
scaling); E 2.61x and 2.59x. Increment 24's AC hardened runs gave 2.48x to
2.50x.

## Quality (identical in all twelve runs)

| domain | worst angle | max degree | within tolerance | Delaunay violations | mesh sha256 |
|---|---:|---:|---|---:|---|
| tile | 0.6296° | 74 | yes | 0 of 692,056 | `11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c` |
| quarter | 0.3955° | 18 | yes | 0 of 641,791 | `ccebf96a86c6c5e244e4a0281919de4e866fcfe789b66024a290ac2af33771a1` |

## Method

- `tools/bench.py run`, blob `2b7dca1aa1184afea53b33c357ab3b836403d2f3` (the
  branch's), driving every tree through `--tree`, `--hardening on`,
  `--repeats 5`. Defaults otherwise: DEM
  `tests/fixtures/dem_archive/7908_3_10m_z33.tif` (sha256 `aabd0cbc…`),
  tolerance 1 m, domains `tile` and `quarter`, threads 0 (CLI default) and
  1 to 20, interleaved. Release `_core` from `bench.build()` into each
  tree's `build-bench/`: B `0b3362e5…` (the same binary as `24-on-r2`), N
  `91391626…`, E `1aa0f7c4…`.
- Apple M1 Max (8 P + 2 E, 32 GiB), macOS 27.0, Python 3.14.7, numpy 2.5.3,
  AppleClang 21.0.0. `caffeinate -i` was held during each batch. No other
  agent or build ran. `pmset -g batt` was recorded before and after every run
  (in each `run.json` and in the logs): AC, 100 %, charged, every time.
- Three batches, 22:04 to 22:45 UTC (`pairs.sh`, `pairs2.sh`, `pairs3.sh`,
  each with its `.log`): batch 1 ran B N N B, batch 2 ran N B B N, and batch 3
  ran E B N E. The tables come from `summarize.py` (batch 1's pairs,
  `tables.md`) and `summarize_all.py` (all runs, `tables-all.md`).
- **The first base run is left out of the pooled figures.** It ran during
  XProtect scans (load average 19 to 31). At 1 thread it is 4 to 10 % slower than the
  other four base runs of the same binary, and `bench.py` judged it a
  REGRESSION against `24-on-r2`, which is the same binary. Including it
  narrows N/B at 1 thread. In batch 1 alone (`tables.md`) the pooled change
  is -2.1 % to +3.0 %, which would have hidden the effect. That is why
  batches 2 and 3 were run. During batches 2 and 3 `dasd` held about one
  core (load around 7).
- The experiment was sanitize-first. ASan and UBSan builds
  (`-fsanitize=address,undefined`, Debug) of `test_refinement_scan`,
  `_scan_frozen`, `_scan_offnode`, `_refine`, `_frozen_oracle` and
  `_refine_points` all passed on E's tree before its Release build was timed.
- N runs are marked dirty because of the untracked evidence directories.
  No source differed from `c4fb2bf`. E is dirty by its diff.

## Files and clean-up

- Run directories in `bench.py`'s format: `23b-base-b4bcdc3{,-r2..-r5}/`,
  `23b{,-r2..-r5}/` and `23b-exp-unswitch{,-r2}/`.
- `23b-acceptance/`: the drivers, the logs, both summarizers and their
  tables, and `exp-unswitch.diff`.
- The quality meshes (about 690 MB of ASCII VTK), the base worktree and the
  experiment worktree were in this session's scratchpad and have been
  deleted. To regenerate them, `git worktree add --detach <dir> b4bcdc3` and
  the same at `c4fb2bf` with `git apply exp-unswitch.diff`. Build each tree
  with `bench.build(bench.make_runner(), Path(T), "on")`, then run
  `pairs*.sh W B S [E]`. The meshes land in `S/meshes/<label>/`.

## Log
- 2026-10-03T22:04:28Z: both trees built (base _core 0b3362e5..., same as 24-on; 23b _core 91391626...). Batch started: order base, 23b, 23b, base; AC, 100 %, charged.
- 2026-10-03T22:14:19Z: batch 1 done (all AC 100 %). Pooled refine within -2.1..+3.0 %, but same-binary drift up to 9.4 % and per-pair up to +5.7 % (tile, 1 thread, N2/B2); XProtect scans and load average 19-31 during and after the batch. Waiting for the machine to quiet, then batch 2 in the order 23b, base, base, 23b.
- 2026-10-03T22:34:02Z: batch 2 done (N B B N, AC 100 %, load ~7 from dasd). Over all four pairs 23b is slower at every low thread count in every run: tile 1 thread base 0.512-0.523 s (first base run 0.543, disturbed) vs 23b 0.541-0.552 s, pooled +4.1 %; per pair +5.7 % and +5.5 % in two pairs. Same-binary spread among the clean runs is about 2 %, so this is outside noise though under bench.py's 5 % threshold pooled. Mesh hashes identical in all eight runs. Suspect, not yet measured: the per-node `on_frozen` test added to scan.hpp's inner loop. Next: a scratch experiment that unswitches that loop.
- 2026-10-03T22:45Z: batch 3 (E B N E) done, AC 100 %. E within -2.4..+1.6 % of base; N +4.7 / +4.3 % at 1 thread. Verdict REGRESSION.
