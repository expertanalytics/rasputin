# 20c-2 acceptance: the soft quality criterion (R7) and the line split (R8) (@perf, 2026-10-07)

Branch `worktree-soft-quality-2`, code `d3c939d9` and `96119508`, tests up
to `600b2acb` (the commit measured), against master `41bda81a` (20c-1
merged). The rules: `docs/increments/20c-soft-quality.md`, "PRs and gates"
(20c-2's row), "`@perf`'s acceptance", and "The split-phase limit, judged
after the fix"; `docs/increments/README.md`, "Acceptance".

Words used below: **20c-1** is master `41bda81a`; **master** in the
triangle gate is `2764ef71`, the merge just before 20c-1 (its meshes equal
`ed125121`'s, the master the gate's figures came from); **gain 0** is
20c-2's default (`--start-quality-gain 0`); **gain −1** turns R7 and R8
off; **feet off** is `--no-constraint-feet`; **split phase** is the
`--stats` phase `refine: split + flip (serial)`; **start quality** is
`refine: start quality`.

## Verdict

- **Every mesh gate of 20c-2 passes** on both catchments (first generated
  table): triangles 0.9253 × and 0.8865 × master's (limit 0.95 ×); slivers
  under 1° 437 and 296, 0.6436 × and 0.6836 × 20c-1's (limit 0.75 ×);
  worst angle 1.000000 × 20c-1's on Lagan (the same value to every printed
  digit, 0.00238402811°) and 12.68 × on Numedalslagen (0.00862300263°);
  gain −1 gives 20c-1's mesh bit for bit in all 7 repeats on both
  catchments and on the tile; one mesh over 7 repeats at gain 0; largest
  height error under 10 m; no constrained Delaunay violation.
- **1 m benchmark and thread sweep (`tools/bench.py run`): ACCEPTED by
  bench.py** (its 5 % threshold; every run after the first printed `ACCEPTED`,
  each judged against the run stored before it),
  and the meshes equal 20c-1's on both domains. But every one of the 42
  cells is slower, +0.1 % to +3.1 %, against a base-to-base spread of
  median 0.7 %, max 3.5 %. The 1-to-20-thread speed-up is 2.48 × against
  2.49 × (tile) and 2.43 × against 2.50 × (quarter).
- **The split-phase limit (revised for 20c-1, holding for 20c-2): not met.**
  - Point 1 (tile, identical meshes): split phase +1.19 % against feet
    off (met) and **+3.36 % against 20c-1 (not met)** over the 6 timed
    repeats. A longer tile series run afterwards (20 timed runs each,
    alternating) gives +1.29 %; both are reported, and the first is not
    replaced by the second.
  - Point 2 (catchments, against what Ola runs today, 20c-1): **split
    phase +18.96 % (Lagan) and +4.22 % (Numedalslagen), not met; whole run
    −0.75 % (Lagan, met) and +3.23 % (Numedalslagen, not met).** The
    meshes differ: at gain 0 refinement adds 17.2 % and 6.8 % more points
    (54 097 against 46 144; 195 799 against 183 368), and per split the
    phase is +1.47 % and −2.40 %.
  - Point 3 (reported, not gated): start quality +61.7 % (Lagan, +0.219 s)
    and +59.4 % (Numedalslagen, +0.320 s) against 20c-1; refine +39.6 %
    and +29.8 %. At gain −1 every catchment phase is within 1.5 % of 20c-1
    (on the tile up to +2.98 %, the split phase, on identical meshes).
- **Overall: REGRESSION**, on the time limit only: Numedalslagen's whole
  run +3.23 % (limit 2 %), and the split phase in seconds +18.96 % and
  +4.22 % (limit 2 %), with the tile's split phase +3.36 % at 6 repeats.
  Nothing has been profiled. What is measured: the extra split-phase
  seconds on the catchments go with the extra splits (per split +1.5 %
  and −2.4 %); the extra refine time is mostly the start quality phase.
  Why start quality costs about 60 % more while it adds 29 % and 34 %
  fewer points is **not known**; by reading the diff, R7's read-only
  cavity is computed for every candidate and R8 adds splits, but no
  profile says how the time divides. Also by reading only: the one change
  on the split phase's path is `must_flip` now calling the new
  `quad_flips` (`include/terrain/mesh/lawson.hpp`).

## Method

- **Machine and power.** Apple M1 Max, macOS 27.0, 10 hardware threads.
  AC power for every run, battery 100 % charged throughout. `pmset -g batt`
  was read before and after each run (`raw/*.power`, `raw-tile/*.power`;
  `run.json` `power` for the benchmark). The scripts stop before a run off
  AC and set aside a run whose power changed; none was set aside.
  `pmset -g log` shows no sleep, wake or power-source event between 17:55
  and 19:00 local time (runs 15:57Z to 16:59Z).
- **Builds.** `tools/bench.py`'s `build()`: Release, AppleClang, hardening
  on (`libc++ fast` in every run), into `<tree>/build-bench/pkg`. 20c-1:
  a scratch worktree of `41bda81a`, `_core` sha256 `86765ff1…` (the same
  as 20c-1's `69f37d1c` in the re-time). 20c-2: this worktree at
  `600b2acb`, `_core` sha256 `0922d7a4…`. Master `2764ef71`: a scratch
  worktree, `_core` sha256 `19d71363…` (the same as `ed125121`'s), built
  by a one-thread tile run of `bench.py` whose evidence was thrown away.
  Every catchment and tile run prints `tin_engine.__file__` and
  `_core.__file__` (first two lines of each log; table below) and
  `drive.py` refuses to run when either is outside its `pkg` tree; after
  each benchmark run `pairs.sh` prints the same and the `_core` sha256
  (`bench/pairs.log`). Python: this worktree's `.venv` (Python 3.14). The
  20c-2 benchmark runs record `dirty: true` because this directory was
  untracked during the runs; nothing else differed from `600b2acb`.
- **1 m benchmark** (`pairs.sh`; log `bench/pairs.log`): `tools/bench.py
  run` with its defaults: the 7908_3 tile and
  `docs/benchmarks/2026-09-26/quarter.geojson`, tolerance 1, threads
  default and 1 to 20, 5 repeats, hardening on. Three pairs, order 20c-1,
  20c-2, 20c-2, 20c-1, 20c-1, 20c-2. The pooled table from
  `bench_summarize.py` is the one that counts.
- **Catchments and tile** (`run.sh`, `drive.py`; log `runs.log`): Lagan
  (Copernicus GLO-30, SMHI outline, CORINE 2018, tolerance 10 m),
  Numedalslagen (DTM10 UTM33, NVE outline, CORINE 2018, tolerance 10 m)
  and the 1 m tile (tolerance 1), commands in `run.sh`. Four builds timed:
  20c-1, gain 0, gain −1, feet off. Seven repeats (r0 to r6), order
  reversed on every other repeat; r0 is a warm-up, left out of the timing
  medians, so each median is of 6. Master `2764ef71` ran twice per
  catchment (`run.sh pre 0 1`), for the triangle gate only. Threads: the
  CLI default (10). The tile series (`run.sh tile 0 .. 20`, log
  `runs-tile.log`, `raw-tile/`) was added after the 6-run figure for point
  1 missed, to see how much of it is noise.
- **Mesh measures** are exact, from the output arrays (`drive.py`, after
  the run): triangle count, count under 1°, worst angle at full
  precision, `bench.py`'s constrained Delaunay check, and a hash of the
  arrays and of the `.vtk` from `POINTS` on. The height error is
  `max_error_m` from `--stats`; the quality-start counts are from
  `--stats` too.
- **Tables** below are generated by `python3 summarize.py`, which reads
  only `raw/`, `raw-tile/`, `bench/`, and for the cross-check rows
  `../20c-1/fix-69f37d1c/raw/`.
- **Meshes** were in the session scratchpad and are not kept. To
  regenerate them, make worktrees of `41bda81a` and `2764ef71` at the
  paths in the scripts, then run `pairs.sh pre`, `pairs.sh` (it builds
  both trees), `run.sh pre 0 1` and `run.sh 0 1 2 3 4 5 6`.

<!-- tables: generated by summarize.py -->

### The 20c-2 gates (`docs/increments/20c-soft-quality.md`, PRs and gates)

master = 2764ef71 (the merge before 20c-1); 20c-1 = master 41bda81a; 20c-2 = gain 0, the CLI default.

| catchment | gate | threshold | measured | result |
|---|---|---|---|---|
| Lagan | triangles ≤ 0.95 × master's | ≤ 820702.15 (863 897 × 0.95) | 799 372 (0.9253 ×) | pass |
| Lagan | slivers under 1° ≤ 0.75 × 20c-1's | ≤ 509.25 (679 × 0.75) | 437 (0.6436 ×) | pass |
| Lagan | worst angle ≥ 0.95 × 20c-1's | ≥ 0.00226482671° (0.00238402811° × 0.95) | 0.00238402811° (1.000000 ×) | pass |
| Lagan | --start-quality-gain -1 bit-identical to 20c-1 | arrays and .vtk from POINTS on equal, every repeat | 7 vs 7 runs, equal | pass |
| Lagan | determinism (gain 0) | one mesh over every repeat | 7 runs, 1 distinct | pass |
| Lagan | largest height error (gain 0) | ≤ 10 m | 9.999930432 m | pass |
| Lagan | constrained Delaunay (gain 0) | 0 violations | 0 of 1 015 095 checked | pass |
| Numedalslagen | triangles ≤ 0.95 × master's | ≤ 1222967.30 (1 287 334 × 0.95) | 1 141 207 (0.8865 ×) | pass |
| Numedalslagen | slivers under 1° ≤ 0.75 × 20c-1's | ≤ 324.75 (433 × 0.75) | 296 (0.6836 ×) | pass |
| Numedalslagen | worst angle ≥ 0.95 × 20c-1's | ≥ 0.000646271124° (0.000680285394° × 0.95) | 0.00862300263° (12.675566 ×) | pass |
| Numedalslagen | --start-quality-gain -1 bit-identical to 20c-1 | arrays and .vtk from POINTS on equal, every repeat | 7 vs 7 runs, equal | pass |
| Numedalslagen | determinism (gain 0) | one mesh over every repeat | 7 runs, 1 distinct | pass |
| Numedalslagen | largest height error (gain 0) | ≤ 10 m | 9.999935150 m | pass |
| Numedalslagen | constrained Delaunay (gain 0) | 0 violations | 0 of 1 541 174 checked | pass |

### Mesh and quality, exact from the output arrays (every repeat; r0 included)

| catchment | build | triangles | < 1° | share < 1° | worst | max degree | max height error m | Delaunay violations | identical across runs | arrays sha256 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | master 41bda81a (20c-1) | 464 290 | 1 | 0.0002 % | 0.629598608° | 74 | 0.999976 | 0 of 692 056 | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | 20c-2 600b2acb, gain 0 (default) | 464 290 | 1 | 0.0002 % | 0.629598608° | 74 | 0.999976 | 0 of 692 056 | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | 20c-2 600b2acb, --start-quality-gain -1 | 464 290 | 1 | 0.0002 % | 0.629598608° | 74 | 0.999976 | 0 of 692 056 | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | 20c-2 600b2acb, --no-constraint-feet | 464 290 | 1 | 0.0002 % | 0.629598608° | 74 | 0.999976 | 0 of 692 056 | 7 runs, 1 distinct | 5866d649a19a |
| Lagan | master 2764ef71 (before 20c-1) | 863 897 | 1598 | 0.1850 % | 0.000412127124° | 18 | 9.999973 | 0 of 1 129 623 | 2 runs, 1 distinct | 66056f56baf4 |
| Lagan | master 41bda81a (20c-1) | 867 612 | 679 | 0.0783 % | 0.00238402811° | 17 | 9.999754 | 0 of 1 128 715 | 7 runs, 1 distinct | 6e397e1e78c8 |
| Lagan | 20c-2 600b2acb, gain 0 (default) | 799 372 | 437 | 0.0547 % | 0.00238402811° | 15 | 9.999930 | 0 of 1 015 095 | 7 runs, 1 distinct | aaf4136a8720 |
| Lagan | 20c-2 600b2acb, --start-quality-gain -1 | 867 612 | 679 | 0.0783 % | 0.00238402811° | 17 | 9.999754 | 0 of 1 128 715 | 7 runs, 1 distinct | 6e397e1e78c8 |
| Lagan | 20c-2 600b2acb, --no-constraint-feet | 753 619 | 1458 | 0.1935 % | 0.000391667745° | 18 | 9.999930 | 0 of 964 355 | 7 runs, 1 distinct | ce555258b74f |
| Numedalslagen | master 2764ef71 (before 20c-1) | 1 287 334 | 1644 | 0.1277 % | 0.000398708734° | 19 | 9.999972 | 0 of 1 781 776 | 2 runs, 1 distinct | 89d7629522c8 |
| Numedalslagen | master 41bda81a (20c-1) | 1 290 807 | 433 | 0.0335 % | 0.000680285394° | 19 | 9.999972 | 0 of 1 783 018 | 7 runs, 1 distinct | 93f93b681dd2 |
| Numedalslagen | 20c-2 600b2acb, gain 0 (default) | 1 141 207 | 296 | 0.0259 % | 0.00862300263° | 19 | 9.999935 | 0 of 1 541 174 | 7 runs, 1 distinct | dac84b2b640a |
| Numedalslagen | 20c-2 600b2acb, --start-quality-gain -1 | 1 290 807 | 433 | 0.0335 % | 0.000680285394° | 19 | 9.999972 | 0 of 1 783 018 | 7 runs, 1 distinct | 93f93b681dd2 |
| Numedalslagen | 20c-2 600b2acb, --no-constraint-feet | 1 085 653 | 5507 | 0.5073 % | 7.22481777e-06° | 19 | 9.999934 | 0 of 1 479 600 | 7 runs, 1 distinct | 4b432d60dc4e |

### Quality-start and refinement counts (`--stats`; the same in every repeat, else listed)

| catchment | build | start: points added | start: tries skipped (all reasons) | start: of them, without gain (R7) | start: points snapped to lines (feet) | start: lines split (R8) | refinement: points added |
|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | master 41bda81a (20c-1) | 503 | 252 | - | 0 | - | 219 837 |
| 1 m tile (tolerance 1) | 20c-2 600b2acb, gain 0 (default) | 503 | 252 | 0 | 0 | 0 | 219 837 |
| 1 m tile (tolerance 1) | 20c-2 600b2acb, --start-quality-gain -1 | 503 | 252 | 0 | 0 | 0 | 219 837 |
| 1 m tile (tolerance 1) | 20c-2 600b2acb, --no-constraint-feet | 503 | 252 | 0 | 0 | 0 | 219 837 |
| Lagan | master 2764ef71 (before 20c-1) | 195 158 | 110 556 | - | - | - | 46 560 |
| Lagan | master 41bda81a (20c-1) | 192 169 | 104 936 | - | 5 930 | - | 46 144 |
| Lagan | 20c-2 600b2acb, gain 0 (default) | 136 273 | 111 408 | 34 600 | 3 488 | 14 025 | 54 097 |
| Lagan | 20c-2 600b2acb, --start-quality-gain -1 | 192 169 | 104 936 | 0 | 5 930 | 0 | 46 144 |
| Lagan | 20c-2 600b2acb, --no-constraint-feet | 126 404 | 123 780 | 35 613 | 0 | 0 | 56 735 |
| Numedalslagen | master 2764ef71 (before 20c-1) | 316 491 | 79 591 | - | - | - | 183 874 |
| Numedalslagen | master 41bda81a (20c-1) | 314 843 | 75 663 | - | 4 233 | - | 183 368 |
| Numedalslagen | 20c-2 600b2acb, gain 0 (default) | 208 852 | 88 243 | 71 590 | 2 565 | 20 579 | 195 799 |
| Numedalslagen | 20c-2 600b2acb, --start-quality-gain -1 | 314 843 | 75 663 | 0 | 4 233 | 0 | 183 368 |
| Numedalslagen | 20c-2 600b2acb, --no-constraint-feet | 189 036 | 104 679 | 66 173 | 0 | 0 | 205 896 |

### Cross-check with 20c-1's re-time (`../20c-1/fix-69f37d1c/raw`, another session)

| catchment | meshes compared (arrays and .vtk) | result |
|---|---|---|
| 1 m tile (tolerance 1) | master 41bda81a (20c-1) vs 20c-1 69f37d1c | equal |
| Lagan | master 41bda81a (20c-1) vs 20c-1 69f37d1c | equal |
| Lagan | master 2764ef71 (before 20c-1) vs master ed125121 | equal |
| Numedalslagen | master 41bda81a (20c-1) vs 20c-1 69f37d1c | equal |
| Numedalslagen | master 2764ef71 (before 20c-1) vs master ed125121 | equal |

### Phases: median seconds over repeats r1 and on ([6] runs per build), threads 10 (the CLI default)

| catchment | phase | master 41bda81a (20c-1): median (min to max) | 20c-2 600b2acb, gain 0 (default): median (min to max) | 20c-2 600b2acb, --start-quality-gain -1: median (min to max) | 20c-2 600b2acb, --no-constraint-feet: median (min to max) | gain 0 vs -1 | gain 0 vs 20c-1 | gain -1 vs 20c-1 | gain 0 vs feet off |
|---|---|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | refine | 0.1909 (0.1865 to 0.1954) | 0.1953 (0.1925 to 0.1977) | 0.1949 (0.1907 to 0.1970) | 0.1942 (0.1882 to 0.1963) | +0.18 % | +2.31 % | +2.13 % | +0.55 % |
| 1 m tile (tolerance 1) | refine: start quality | 0.0007 (0.0007 to 0.0009) | 0.0010 (0.0010 to 0.0010) | 0.0007 (0.0007 to 0.0008) | 0.0010 (0.0009 to 0.0010) | +41.31 % | +43.75 % | +1.73 % | +3.25 % |
| 1 m tile (tolerance 1) | refine: scan (parallel) | 0.0571 (0.0569 to 0.0590) | 0.0574 (0.0571 to 0.0594) | 0.0577 (0.0568 to 0.0583) | 0.0574 (0.0565 to 0.0577) | -0.50 % | +0.52 % | +1.03 % | +0.02 % |
| 1 m tile (tolerance 1) | refine: split + flip (serial) | 0.1067 (0.1030 to 0.1113) | 0.1103 (0.1057 to 0.1111) | 0.1099 (0.1046 to 0.1116) | 0.1090 (0.1031 to 0.1111) | +0.37 % | +3.36 % | +2.98 % | +1.19 % |
| 1 m tile (tolerance 1) | total | 0.4320 (0.4210 to 0.4330) | 0.4340 (0.4320 to 0.4380) | 0.4335 (0.4280 to 0.4370) | 0.4340 (0.4280 to 0.4390) | +0.12 % | +0.46 % | +0.35 % | +0.00 % |
| Lagan | refine | 0.5759 (0.5680 to 0.6115) | 0.8040 (0.7997 to 0.8211) | 0.5777 (0.5739 to 0.5812) | 0.7299 (0.7260 to 0.7351) | +39.17 % | +39.60 % | +0.31 % | +10.14 % |
| Lagan | refine: start quality | 0.3549 (0.3487 to 0.3832) | 0.5738 (0.5690 to 0.5883) | 0.3547 (0.3519 to 0.3571) | 0.5030 (0.4996 to 0.5079) | +61.77 % | +61.66 % | -0.07 % | +14.07 % |
| Lagan | refine: scan (parallel) | 0.0505 (0.0499 to 0.0520) | 0.0517 (0.0514 to 0.0527) | 0.0505 (0.0503 to 0.0529) | 0.0520 (0.0514 to 0.0527) | +2.25 % | +2.36 % | +0.10 % | -0.59 % |
| Lagan | refine: split + flip (serial) | 0.0382 (0.0371 to 0.0392) | 0.0455 (0.0438 to 0.0461) | 0.0384 (0.0376 to 0.0393) | 0.0415 (0.0410 to 0.0422) | +18.38 % | +18.96 % | +0.49 % | +9.66 % |
| Lagan | total | 72.7235 (71.8020 to 73.6140) | 72.1775 (71.2740 to 73.1020) | 72.4900 (71.7840 to 72.7960) | 72.1510 (71.3620 to 72.9870) | -0.43 % | -0.75 % | -0.32 % | +0.04 % |
| Numedalslagen | refine | 1.1241 (1.1174 to 1.1300) | 1.4589 (1.4499 to 1.4687) | 1.1343 (1.1238 to 1.1705) | 1.3342 (1.3262 to 1.3457) | +28.62 % | +29.78 % | +0.91 % | +9.34 % |
| Numedalslagen | refine: start quality | 0.5389 (0.5354 to 0.5422) | 0.8591 (0.8530 to 0.8670) | 0.5448 (0.5400 to 0.5774) | 0.7472 (0.7377 to 0.7594) | +57.68 % | +59.40 % | +1.09 % | +14.97 % |
| Numedalslagen | refine: scan (parallel) | 0.2867 (0.2849 to 0.2884) | 0.2981 (0.2977 to 0.2988) | 0.2881 (0.2848 to 0.2894) | 0.3020 (0.2999 to 0.3054) | +3.48 % | +3.98 % | +0.48 % | -1.29 % |
| Numedalslagen | refine: split + flip (serial) | 0.1559 (0.1547 to 0.1580) | 0.1625 (0.1607 to 0.1645) | 0.1581 (0.1567 to 0.1604) | 0.1504 (0.1488 to 0.1533) | +2.77 % | +4.22 % | +1.41 % | +8.04 % |
| Numedalslagen | total | 6.6755 (6.6630 to 6.7070) | 6.8910 (6.8640 to 6.9030) | 6.7050 (6.6750 to 6.8310) | 6.7010 (6.6580 to 7.0330) | +2.77 % | +3.23 % | +0.44 % | +2.84 % |

### The split-phase limit (`docs/increments/20c-soft-quality.md`, "The split-phase limit, judged after the fix")

| point | measure | limit | measured (medians) | result |
|---|---|---|---|---|
| 1. tile, no cost where no foot is placed | split phase, gain 0 vs 20c-2 600b2acb, --no-constraint-feet | ≤ +2 % | +1.19 % (meshes identical: True) | met |
| 1. tile, no cost where no foot is placed | split phase, gain 0 vs master 41bda81a (20c-1) | ≤ +2 % | +3.36 % (meshes identical: True) | NOT MET |
| 2. Lagan, no slower than what Ola runs today | refine: split + flip (serial), gain 0 vs 20c-1 | ≤ +2 % | +18.96 % | NOT MET |
| 2. Lagan, no slower than what Ola runs today | total, gain 0 vs 20c-1 | within 2 % | -0.75 % | met |
| 3. Lagan, reported, not gated | refine: start quality, gain 0 vs 20c-1 / vs feet off | - | +61.66 % / +14.07 % | - |
| 3. Lagan, reported, not gated | refine, gain 0 vs 20c-1 / vs feet off | - | +39.60 % / +10.14 % | - |
| 2. Numedalslagen, no slower than what Ola runs today | refine: split + flip (serial), gain 0 vs 20c-1 | ≤ +2 % | +4.22 % | NOT MET |
| 2. Numedalslagen, no slower than what Ola runs today | total, gain 0 vs 20c-1 | within 2 % | +3.23 % | NOT MET |
| 3. Numedalslagen, reported, not gated | refine: start quality, gain 0 vs 20c-1 / vs feet off | - | +59.40 % / +14.97 % | - |
| 3. Numedalslagen, reported, not gated | refine, gain 0 vs 20c-1 / vs feet off | - | +29.78 % / +9.34 % | - |

### The split phase per split (for reading point 2: the meshes differ, so the work differs)

Splits are `points_inserted` from `--stats` (refinement's main loop). µs per split = median split-phase
seconds / splits. Reported beside the limit, which is written in seconds; not a gate.

| catchment | splits, base | splits, head | splits, headm1 | splits, nofeet | µs per split, base | µs per split, head | µs per split, headm1 | µs per split, nofeet | splits, gain 0 vs 20c-1 | per split, gain 0 vs 20c-1 | per split, gain 0 vs feet off |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | 219 837 | 219 837 | 219 837 | 219 837 | 0.4853 | 0.5016 | 0.4998 | 0.4957 | +0.00 % | +3.36 % | +1.19 % |
| Lagan | 46 144 | 54 097 | 46 144 | 56 735 | 0.8285 | 0.8407 | 0.8326 | 0.7310 | +17.24 % | +1.47 % | +15.00 % |
| Numedalslagen | 183 368 | 195 799 | 183 368 | 205 896 | 0.8502 | 0.8298 | 0.8621 | 0.7304 | +6.78 % | -2.40 % | +13.61 % |

### A longer tile series for point 1 (`run.sh tile 0 .. 20`, `raw-tile/`; r0 a warm-up; 20 and 20 timed runs, base and head alternating)

| phase | master 41bda81a (20c-1): median (min to max) | 20c-2 600b2acb, gain 0 (default): median (min to max) | gain 0 vs 20c-1 |
|---|---|---|---|
| refine | 0.1929 (0.1848 to 0.1954) | 0.1949 (0.1878 to 0.2002) | +1.07 % |
| refine: start quality | 0.0007 (0.0007 to 0.0007) | 0.0010 (0.0009 to 0.0011) | +45.15 % |
| refine: scan (parallel) | 0.0576 (0.0568 to 0.0588) | 0.0574 (0.0568 to 0.0602) | -0.27 % |
| refine: split + flip (serial) | 0.1084 (0.1015 to 0.1112) | 0.1098 (0.1033 to 0.1122) | +1.29 % |
| total | 0.4415 (0.4320 to 0.4470) | 0.4330 (0.4260 to 0.4440) | -1.93 % |

Every mesh in the series equal to each other and to `raw/`'s tile meshes: True.

### Spread of 20c-1's own runs (base, (max - min) / median)

| catchment | refine: start quality | refine: split + flip (serial) | refine | total |
|---|---|---|---|---|
| 1 m tile (tolerance 1) | 29.4 % | 7.7 % | 4.7 % | 2.8 % |
| Lagan | 9.7 % | 5.6 % | 7.6 % | 2.5 % |
| Numedalslagen | 1.3 % | 2.1 % | 1.1 % | 0.7 % |

### Power, every run (`pmset -g batt` before and after; 1 m benchmark: `run.json`)

| run | power | battery before -> after | process s |
|---|---|---|---|
| lagan-base-r0 | AC | 100 % charged -> 100 % charged | 73.9 |
| lagan-base-r1 | AC | 100 % charged -> 100 % charged | 74.3 |
| lagan-base-r2 | AC | 100 % charged -> 100 % charged | 73.1 |
| lagan-base-r3 | AC | 100 % charged -> 100 % charged | 73.6 |
| lagan-base-r4 | AC | 100 % charged -> 100 % charged | 74.2 |
| lagan-base-r5 | AC | 100 % charged -> 100 % charged | 73.9 |
| lagan-base-r6 | AC | 100 % charged -> 100 % charged | 75.0 |
| lagan-head-r0 | AC | 100 % charged -> 100 % charged | 73.5 |
| lagan-head-r1 | AC | 100 % charged -> 100 % charged | 72.5 |
| lagan-head-r2 | AC | 100 % charged -> 100 % charged | 73.8 |
| lagan-head-r3 | AC | 100 % charged -> 100 % charged | 72.9 |
| lagan-head-r4 | AC | 100 % charged -> 100 % charged | 73.7 |
| lagan-head-r5 | AC | 100 % charged -> 100 % charged | 74.3 |
| lagan-head-r6 | AC | 100 % charged -> 100 % charged | 73.1 |
| lagan-headm1-r0 | AC | 100 % charged -> 100 % charged | 73.7 |
| lagan-headm1-r1 | AC | 100 % charged -> 100 % charged | 73.9 |
| lagan-headm1-r2 | AC | 100 % charged -> 100 % charged | 73.3 |
| lagan-headm1-r3 | AC | 100 % charged -> 100 % charged | 74.0 |
| lagan-headm1-r4 | AC | 100 % charged -> 100 % charged | 73.1 |
| lagan-headm1-r5 | AC | 100 % charged -> 100 % charged | 73.7 |
| lagan-headm1-r6 | AC | 100 % charged -> 100 % charged | 74.1 |
| lagan-nofeet-r0 | AC | 100 % charged -> 100 % charged | 73.9 |
| lagan-nofeet-r1 | AC | 100 % charged -> 100 % charged | 72.9 |
| lagan-nofeet-r2 | AC | 100 % charged -> 100 % charged | 73.7 |
| lagan-nofeet-r3 | AC | 100 % charged -> 100 % charged | 72.5 |
| lagan-nofeet-r4 | AC | 100 % charged -> 100 % charged | 74.1 |
| lagan-nofeet-r5 | AC | 100 % charged -> 100 % charged | 73.6 |
| lagan-nofeet-r6 | AC | 100 % charged -> 100 % charged | 73.1 |
| lagan-pre-r0 | AC | 100 % charged -> 100 % charged | 74.4 |
| lagan-pre-r1 | AC | 100 % charged -> 100 % charged | 73.7 |
| numed-base-r0 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-base-r1 | AC | 100 % charged -> 100 % charged | 8.6 |
| numed-base-r2 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-base-r3 | AC | 100 % charged -> 100 % charged | 8.6 |
| numed-base-r4 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-base-r5 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-base-r6 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-head-r0 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-head-r1 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-head-r2 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-head-r3 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-head-r4 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-head-r5 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-head-r6 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-headm1-r0 | AC | 100 % charged -> 100 % charged | 8.6 |
| numed-headm1-r1 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-headm1-r2 | AC | 100 % charged -> 100 % charged | 8.7 |
| numed-headm1-r3 | AC | 100 % charged -> 100 % charged | 8.6 |
| numed-headm1-r4 | AC | 100 % charged -> 100 % charged | 8.6 |
| numed-headm1-r5 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-headm1-r6 | AC | 100 % charged -> 100 % charged | 8.6 |
| numed-nofeet-r0 | AC | 100 % charged -> 100 % charged | 8.3 |
| numed-nofeet-r1 | AC | 100 % charged -> 100 % charged | 8.2 |
| numed-nofeet-r2 | AC | 100 % charged -> 100 % charged | 8.2 |
| numed-nofeet-r3 | AC | 100 % charged -> 100 % charged | 8.2 |
| numed-nofeet-r4 | AC | 100 % charged -> 100 % charged | 8.6 |
| numed-nofeet-r5 | AC | 100 % charged -> 100 % charged | 8.2 |
| numed-nofeet-r6 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-pre-r0 | AC | 100 % charged -> 100 % charged | 9.0 |
| numed-pre-r1 | AC | 100 % charged -> 100 % charged | 8.6 |
| tile-base-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r1 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r1 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r6 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-headm1-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r1 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r1 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-nofeet-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| b1m-base-r1 | ac | 100 % charged -> 100 % charged | - |
| b1m-base-r2 | ac | 100 % charged -> 100 % charged | - |
| b1m-base-r3 | ac | 100 % charged -> 100 % charged | - |
| b1m-head-r1 | ac | 100 % charged -> 100 % charged | - |
| b1m-head-r2 | ac | 100 % charged -> 100 % charged | - |
| b1m-head-r3 | ac | 100 % charged -> 100 % charged | - |

### Where each run loaded tin_engine and _core from (first two lines of every `raw/*.log`)

| runs | count | tin_engine.__file__ | _core.__file__ |
|---|---|---|---|
| lagan-base | 7 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/base/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/base/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| lagan-head | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| lagan-headm1 | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| lagan-nofeet | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| lagan-pre | 2 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/pre/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/pre/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| numed-base | 7 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/base/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/base/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| numed-head | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| numed-headm1 | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| numed-nofeet | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| numed-pre | 2 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/pre/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/pre/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| tile-base | 7 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/base/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2/base/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| tile-head | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| tile-headm1 | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| tile-nofeet | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |

### The 1 m benchmark and thread sweep (`tools/bench.py run`, `bench/b1m-*`; `bench_summarize.py bench b1m`)

base = master 41bda81a, head = 20c-2 600b2acb.

base: 3 runs, commits ['41bda81'], power ['ac'], hardening ['libc++ fast']
head: 3 runs, commits ['600b2ac'], power ['ac'], hardening ['libc++ fast']

| domain | mesh sha256 (base) | mesh sha256 (head) | equal |
|---|---|---|---|
| tile | 11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c | 11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c | True |
| quarter | a60fb597b5b738268a3bd27730989392ee86fb030fe060b78fd3a0bcdbb60a34 | a60fb597b5b738268a3bd27730989392ee86fb030fe060b78fd3a0bcdbb60a34 | True |

| domain | threads | refine_s base | refine_s head | change | app_s change | proc_s change |
|---|---:|---:|---:|---:|---:|---:|
| tile | default | 0.1915 | 0.1945 | +1.5 % | +1.2 % | +0.6 % |
| tile | 1 | 0.4861 | 0.4880 | +0.4 % | +0.6 % | +0.1 % |
| tile | 2 | 0.3167 | 0.3188 | +0.7 % | +0.6 % | +0.7 % |
| tile | 3 | 0.2625 | 0.2637 | +0.4 % | +0.4 % | +0.0 % |
| tile | 4 | 0.2370 | 0.2384 | +0.6 % | +0.6 % | +0.3 % |
| tile | 5 | 0.2184 | 0.2204 | +0.9 % | +1.1 % | -0.0 % |
| tile | 6 | 0.2075 | 0.2089 | +0.7 % | +0.3 % | +0.1 % |
| tile | 7 | 0.1973 | 0.1989 | +0.8 % | +0.2 % | +0.3 % |
| tile | 8 | 0.1921 | 0.1950 | +1.5 % | +0.9 % | +0.1 % |
| tile | 9 | 0.1930 | 0.1960 | +1.6 % | +0.8 % | +0.9 % |
| tile | 10 | 0.1932 | 0.1950 | +0.9 % | +0.2 % | +0.4 % |
| tile | 11 | 0.1936 | 0.1951 | +0.8 % | +0.2 % | -0.0 % |
| tile | 12 | 0.1925 | 0.1949 | +1.2 % | +0.5 % | +0.3 % |
| tile | 13 | 0.1892 | 0.1945 | +2.8 % | +1.0 % | +0.4 % |
| tile | 14 | 0.1941 | 0.1948 | +0.3 % | -0.1 % | +0.0 % |
| tile | 15 | 0.1935 | 0.1953 | +1.0 % | +0.8 % | +0.3 % |
| tile | 16 | 0.1906 | 0.1950 | +2.3 % | +1.7 % | +1.1 % |
| tile | 17 | 0.1944 | 0.1945 | +0.1 % | +0.2 % | -0.3 % |
| tile | 18 | 0.1947 | 0.1959 | +0.6 % | +0.3 % | +0.0 % |
| tile | 19 | 0.1907 | 0.1963 | +2.9 % | +0.9 % | +1.0 % |
| tile | 20 | 0.1953 | 0.1965 | +0.6 % | +0.8 % | +0.8 % |
| quarter | default | 0.1745 | 0.1762 | +1.0 % | +0.4 % | +0.4 % |
| quarter | 1 | 0.4298 | 0.4304 | +0.1 % | +0.4 % | +0.2 % |
| quarter | 2 | 0.2794 | 0.2816 | +0.8 % | +0.7 % | +0.3 % |
| quarter | 3 | 0.2290 | 0.2317 | +1.2 % | +0.9 % | +0.5 % |
| quarter | 4 | 0.2087 | 0.2116 | +1.4 % | +1.2 % | +1.0 % |
| quarter | 5 | 0.1932 | 0.1957 | +1.3 % | +0.3 % | -0.2 % |
| quarter | 6 | 0.1824 | 0.1854 | +1.7 % | +0.7 % | +0.1 % |
| quarter | 7 | 0.1753 | 0.1775 | +1.2 % | +0.5 % | +0.4 % |
| quarter | 8 | 0.1701 | 0.1728 | +1.6 % | +0.8 % | +0.1 % |
| quarter | 9 | 0.1739 | 0.1768 | +1.6 % | +0.7 % | +0.3 % |
| quarter | 10 | 0.1735 | 0.1764 | +1.7 % | +1.2 % | +0.5 % |
| quarter | 11 | 0.1712 | 0.1755 | +2.5 % | +1.3 % | +0.7 % |
| quarter | 12 | 0.1734 | 0.1762 | +1.6 % | +0.9 % | +0.5 % |
| quarter | 13 | 0.1732 | 0.1765 | +2.0 % | +1.0 % | +0.6 % |
| quarter | 14 | 0.1716 | 0.1747 | +1.8 % | +0.5 % | +0.3 % |
| quarter | 15 | 0.1722 | 0.1768 | +2.7 % | +1.7 % | +0.5 % |
| quarter | 16 | 0.1745 | 0.1766 | +1.2 % | +0.5 % | +0.2 % |
| quarter | 17 | 0.1726 | 0.1765 | +2.3 % | +1.1 % | +1.0 % |
| quarter | 18 | 0.1734 | 0.1771 | +2.2 % | +1.5 % | +0.6 % |
| quarter | 19 | 0.1742 | 0.1767 | +1.4 % | +0.7 % | +0.8 % |
| quarter | 20 | 0.1716 | 0.1769 | +3.1 % | +1.2 % | +0.5 % |

refine_s pooled change over every cell: +0.1 .. +3.1 %
base-vs-base noise (per cell, max/min of the base runs' medians): median 0.7 %, max 3.5 %

| domain | side | worst angle | median angle | share < 1° | max degree | within tolerance | Delaunay violations (checked) |
|---|---|---|---|---|---|---|---|
| tile | base (1 distinct of 3) | 0.6296° | 45.00° | 0.0002 % | 74 | True | 0 (692056) |
| tile | head (1 distinct of 3) | 0.6296° | 45.00° | 0.0002 % | 74 | True | 0 (692056) |
| quarter | base (1 distinct of 3) | 0.3955° | 45.00° | 0.0028 % | 18 | True | 0 (641792) |
| quarter | head (1 distinct of 3) | 0.3955° | 45.00° | 0.0028 % | 18 | True | 0 (641792) |

| domain | side | refine_s 1 thread | refine_s 20 threads | speed-up |
|---|---|---:|---:|---:|
| tile | base | 0.4861 | 0.1953 | 2.49 x |
| tile | head | 0.4880 | 0.1965 | 2.48 x |
| quarter | base | 0.4298 | 0.1716 | 2.50 x |
| quarter | head | 0.4304 | 0.1769 | 2.43 x |
<!-- end of generated tables -->
