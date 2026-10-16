# 20c-2 re-time after the speed fix 65ae0792 (@perf, 2026-10-07)

Branch `worktree-soft-quality-2` at `65ae0792` (`@developer`'s speed fix in
`include/terrain/mesh/quality.hpp` and `include/terrain/mesh/lawson.hpp`,
said to leave the output unchanged), against master `41bda81a` (20c-1 merged).
The first acceptance, of `600b2acb`, is `../README.md` (commit `fa136aa7`).
The rules: `docs/increments/20c-soft-quality.md`, "The split-phase limit,
judged after the fix" and "`@perf`'s acceptance"; `docs/increments/README.md`,
"Acceptance".

Words used below: **20c-1** is master `41bda81a`; **gain 0** is 20c-2's
default (`--start-quality-gain 0`); **gain −1** turns R7 (the soft quality
criterion) and R8 (the line split) off; **feet off** is
`--no-constraint-feet`; **split phase** is the `--stats` phase
`refine: split + flip (serial)`; **start quality** is
`refine: start quality`; **whole run** is the `--stats` total.

## Verdict

1. **Byte-identity: confirmed.** Every mesh `65ae0792` produced equals the
   first acceptance's mesh from `600b2acb`, in arrays and in the `.vtk` from
   `POINTS` on, in every repeat: gain 0, gain −1 and feet off on Lagan,
   Numedalslagen and the 1 m tile, and 20c-1 again; on the 1 m benchmark,
   the tile and the quarter. So every 20c-2 mesh gate of `../README.md` holds
   unchanged (triangles 0.9253 × and 0.8865 × master's; slivers under 1°
   0.6436 × and 0.6836 × 20c-1's; worst angle 1.000000 × and 12.68 × 20c-1's;
   gain −1 equals 20c-1 bit for bit), and the quality columns below are the
   same figures again.
2. **Time against 20c-1** (median of 6 after a warm-up; the tile also 20
   runs per build, rotating): the fix removed about two thirds of the start
   quality's excess. Start quality is +23.5 % on Lagan (0.431 s against
   0.349 s; before the fix +61.7 %) and +20.1 % on Numedalslagen (0.642 s
   against 0.535 s; before +59.4 %). The split phase barely moved: +16.3 %
   on Lagan (before +19.0 %), +4.0 % on Numedalslagen (before +4.2 %). The
   whole run is −0.04 % on Lagan and −0.08 % on Numedalslagen (before −0.75 %
   and +3.23 %).
3. **1 m benchmark (`tools/bench.py run`, 3 pairs): ACCEPTED** by bench.py
   (every run after the first printed `ACCEPTED` except the head's first,
   which printed `REGRESSION` on two quarter cells, +6.3 % and +5.1 %,
   judged against the single base run stored before it; the pooled table
   is the one that counts). Pooled over 3 runs a side, refine time changes
   by −0.8 % to +2.6 % per cell, against a base-to-base spread of median
   0.8 %, max 2.7 %; the tile cells are all within ±0.8 %; the quarter cells
   are −0.1 % to +2.6 %, all but one at or under +1.6 %. The meshes equal
   20c-1's on both domains, and the 1-to-20-thread speed-up is 2.49 × against
   2.49 × (tile) and 2.44 × against 2.46 × (quarter).
4. **The revised time limit, point by point** (the generated table "4."):
   - Point 1, tile, where no foot is placed: **met.** Split phase −4.24 %
     against feet off and −4.52 % against 20c-1 over 6 runs; −0.34 % and
     +1.39 % over 20 runs. Limit +2 %.
   - Point 2, Lagan: split phase **+16.29 %, not met** (limit +2 %); whole
     run **−0.04 %, met** (limit within 2 %).
   - Point 2, Numedalslagen: split phase **+4.00 %, not met** (limit +2 %);
     whole run **−0.08 %, met** (limit within 2 %).
   - Point 3, reported, not gated: start quality +23.46 % (Lagan) and
     +20.11 % (Numedalslagen) against 20c-1, +14.48 % and +16.96 % against
     feet off; refine +16.05 % and +11.19 % against 20c-1.
   - `@developer`'s indicative figures against these: start quality 0.433 s
     and 0.652 s indicated, 0.431 s and 0.642 s measured; Numedalslagen whole
     run 6.683 s indicated, 6.6305 s measured (20c-1 in this session 6.6360 s);
     tile split phase +0.0 to +0.5 % indicated, −4.52 % (6 runs) and +1.39 %
     (20 runs) measured.
5. **Overall: REGRESSION on point 2's split phase only**, +16.29 % on Lagan
   and +4.00 % on Numedalslagen (limit +2 %). Every other measure is met.
   What is measured about the remaining split phase: the meshes differ, and
   at gain 0 refinement splits 17.2 % and 6.8 % more times (54 097 against
   46 144; 195 799 against 183 368); per split the phase is −0.81 % and
   −2.60 % against 20c-1. So the extra seconds go with the extra splits,
   not with a dearer split. Why gain 0 leaves more for refinement to split
   is not profiled; by the counts, the quality start adds 29 % and 34 % fewer
   points. The split phase's limit is written in seconds, so this verdict
   stands as measured; whether the limit should count work per split for
   20c-2 is a design question, not a measurement.

## Method

- **Machine and power.** Apple M1 Max, macOS 27.0, 10 hardware threads.
  AC power for every run, battery 100 % charged throughout (the power table
  below: 147 catchment and tile runs, 6 benchmark runs). `pmset -g batt` was
  read before and after each run (`raw/*.power`, `raw-tile/*.power`;
  `run.json` `power` for the benchmark). The scripts stop before a run off AC
  and set aside a run whose power changed; none was set aside. `pmset -g log`
  shows no sleep, wake or power-source event between 19:00 and 21:00 local
  time (runs 17:51Z to 18:49Z).
- **Builds.** `tools/bench.py`'s `build()` (Release, AppleClang, hardening on,
  `libc++ fast` in every run) into `<tree>/build-bench/pkg`, done by the
  first benchmark run of each tree. 20c-1: a scratch worktree of `41bda81a`,
  `_core` sha256 `86765ff1…`, the same file as in the first acceptance.
  20c-2: this worktree at `65ae0792`, `_core` sha256 `28f1db86…` (the first
  acceptance's `600b2acb` build was `0922d7a4…`), the benchmark's
  `tree.commit` `65ae079`. Every catchment and tile run prints
  `tin_engine.__file__` and `_core.__file__` (first two lines of each log;
  table below), and `drive.py` refuses to run when either is outside its
  `pkg` tree; after each benchmark run `pairs.sh` prints the same and the
  `_core` sha256 (`bench/pairs.log`). Python: this worktree's `.venv`
  (Python 3.14). The head benchmark runs record `dirty: true` because this
  directory was untracked during the runs; nothing else differed from
  `65ae0792`.
- **1 m benchmark** (`pairs.sh`; log `bench/pairs.log`): `tools/bench.py run`
  with its defaults: the 7908_3 tile and
  `docs/benchmarks/2026-09-26/quarter.geojson`, tolerance 1, threads default
  and 1 to 20, 5 repeats, hardening on. Three pairs, order 20c-1, 20c-2,
  20c-2, 20c-1, 20c-1, 20c-2. A first attempt stopped after its first run
  because `pairs.sh` took bench.py's exit 2 (its `NO BASELINE` verdict) for a
  failure; the check now stops only on exit 3 or higher (an error), that
  run's evidence was deleted, and the series was started again from its
  first run.
- **Catchments and tile** (`run.sh`, `drive.py`; log `runs.log`): Lagan
  (Copernicus GLO-30, SMHI outline, CORINE 2018, tolerance 10 m),
  Numedalslagen (DTM10 UTM33, NVE outline, CORINE 2018, tolerance 10 m) and
  the 1 m tile (tolerance 1), commands in `run.sh`. Four builds timed: 20c-1,
  gain 0, gain −1, feet off. Seven repeats (r0 to r6), order reversed on
  every other repeat; r0 is a warm-up, left out of the timing medians, so
  each median is of 6. Threads: the CLI default (10).
- **Longer tile series** (`run.sh tile 0 .. 20`, log `runs-tile.log`,
  `raw-tile/`): 20c-1, gain 0 and feet off in turn, order reversed on every
  other repeat, r0 a warm-up, so 20 timed runs each. From r14 on, all three
  builds have runs with a slow parallel scan (up to 0.101 s against a median
  of 0.057 s; `raw-tile/*-r14` to `*-r20`), so something else used the CPU
  then; the cause is not measured. Since the builds alternate, each build has
  some of them, and the medians are reported with them in.
- **Mesh measures** are exact, from the output arrays (`drive.py`, after the
  run): triangle count, count under 1°, worst angle at full precision,
  `bench.py`'s constrained Delaunay check, and a hash of the arrays and of
  the `.vtk` from `POINTS` on. The height error is `max_error_m` from
  `--stats`; the quality-start counts are from `--stats` too.
- **Tables** below are generated by `python3 summarize.py`, which reads only
  `raw/`, `raw-tile/`, `bench/`, and for the byte-identity and before-the-fix
  rows the first acceptance's `../raw/`, `../raw-tile/` and `../bench/`. The
  before-the-fix table compares medians from two sessions; its 20c-1 columns
  show the drift between them: within 2 % in every row but the tile's start quality, a 0.6 to 0.7 ms phase, which is −12.5 %.
- **Meshes** were in the session scratchpad and are not kept. To regenerate
  them, make a worktree of `41bda81a` at the path in the scripts, then run
  `pairs.sh` (it builds both trees), `run.sh 0 1 2 3 4 5 6` and
  `run.sh tile 0 1 … 20`.

<!-- tables: generated by summarize.py -->

### 1. Byte-identity with the first acceptance (`../raw`, `../raw-tile`, `../bench`: fa136aa7, measuring 600b2acb)

Arrays (vertices, triangles, constraint edges) and the .vtk from `POINTS` on, every repeat of both sides; the 1 m benchmark's `mesh_sha256` per domain over all its runs.

| domain | meshes compared | runs | sha256 (first 12) | result |
|---|---|---|---|---|
| 1 m tile (tolerance 1) | master 41bda81a (20c-1) vs master 41bda81a (20c-1) | 7 vs 7 | 5866d649a19a | equal |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, gain 0 (default) vs 20c-2 600b2acb, gain 0 | 7 vs 7 | 5866d649a19a | equal |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, --start-quality-gain -1 vs 20c-2 600b2acb, gain -1 | 7 vs 7 | 5866d649a19a | equal |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, --no-constraint-feet vs 20c-2 600b2acb, feet off | 7 vs 7 | 5866d649a19a | equal |
| Lagan | master 41bda81a (20c-1) vs master 41bda81a (20c-1) | 7 vs 7 | 6e397e1e78c8 | equal |
| Lagan | 20c-2 65ae0792, gain 0 (default) vs 20c-2 600b2acb, gain 0 | 7 vs 7 | aaf4136a8720 | equal |
| Lagan | 20c-2 65ae0792, --start-quality-gain -1 vs 20c-2 600b2acb, gain -1 | 7 vs 7 | 6e397e1e78c8 | equal |
| Lagan | 20c-2 65ae0792, --no-constraint-feet vs 20c-2 600b2acb, feet off | 7 vs 7 | ce555258b74f | equal |
| Numedalslagen | master 41bda81a (20c-1) vs master 41bda81a (20c-1) | 7 vs 7 | 93f93b681dd2 | equal |
| Numedalslagen | 20c-2 65ae0792, gain 0 (default) vs 20c-2 600b2acb, gain 0 | 7 vs 7 | dac84b2b640a | equal |
| Numedalslagen | 20c-2 65ae0792, --start-quality-gain -1 vs 20c-2 600b2acb, gain -1 | 7 vs 7 | 93f93b681dd2 | equal |
| Numedalslagen | 20c-2 65ae0792, --no-constraint-feet vs 20c-2 600b2acb, feet off | 7 vs 7 | 4b432d60dc4e | equal |
| 1 m tile (tolerance 1), longer series | master 41bda81a (20c-1) vs master 41bda81a (20c-1) | 21 vs 21 | 5866d649a19a | equal |
| 1 m tile (tolerance 1), longer series | 20c-2 65ae0792, gain 0 (default) vs 20c-2 600b2acb, gain 0 | 21 vs 21 | 5866d649a19a | equal |
| 1 m tile (tolerance 1), longer series | 20c-2 65ae0792, --no-constraint-feet vs 20c-2 600b2acb, feet off | 21 vs 7 | 5866d649a19a | equal |
| 1 m benchmark, tile | base: 41bda81a vs 41bda81a | 3 vs 3 | 11741a81adfa | equal |
| 1 m benchmark, quarter | base: 41bda81a vs 41bda81a | 3 vs 3 | a60fb597b5b7 | equal |
| 1 m benchmark, tile | head: 65ae0792 vs 600b2acb | 3 vs 3 | 11741a81adfa | equal |
| 1 m benchmark, quarter | head: 65ae0792 vs 600b2acb | 3 vs 3 | a60fb597b5b7 | equal |

### Mesh and quality, exact from the output arrays (every repeat; r0 included)

| catchment | build | triangles | < 1° | share < 1° | worst | max degree | max height error m | Delaunay violations | identical across runs | arrays sha256 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | master 41bda81a (20c-1) | 464 290 | 1 | 0.0002 % | 0.629598608° | 74 | 0.999976 | 0 of 692 056 | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, gain 0 (default) | 464 290 | 1 | 0.0002 % | 0.629598608° | 74 | 0.999976 | 0 of 692 056 | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, --start-quality-gain -1 | 464 290 | 1 | 0.0002 % | 0.629598608° | 74 | 0.999976 | 0 of 692 056 | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, --no-constraint-feet | 464 290 | 1 | 0.0002 % | 0.629598608° | 74 | 0.999976 | 0 of 692 056 | 7 runs, 1 distinct | 5866d649a19a |
| Lagan | master 41bda81a (20c-1) | 867 612 | 679 | 0.0783 % | 0.00238402811° | 17 | 9.999754 | 0 of 1 128 715 | 7 runs, 1 distinct | 6e397e1e78c8 |
| Lagan | 20c-2 65ae0792, gain 0 (default) | 799 372 | 437 | 0.0547 % | 0.00238402811° | 15 | 9.999930 | 0 of 1 015 095 | 7 runs, 1 distinct | aaf4136a8720 |
| Lagan | 20c-2 65ae0792, --start-quality-gain -1 | 867 612 | 679 | 0.0783 % | 0.00238402811° | 17 | 9.999754 | 0 of 1 128 715 | 7 runs, 1 distinct | 6e397e1e78c8 |
| Lagan | 20c-2 65ae0792, --no-constraint-feet | 753 619 | 1458 | 0.1935 % | 0.000391667745° | 18 | 9.999930 | 0 of 964 355 | 7 runs, 1 distinct | ce555258b74f |
| Numedalslagen | master 41bda81a (20c-1) | 1 290 807 | 433 | 0.0335 % | 0.000680285394° | 19 | 9.999972 | 0 of 1 783 018 | 7 runs, 1 distinct | 93f93b681dd2 |
| Numedalslagen | 20c-2 65ae0792, gain 0 (default) | 1 141 207 | 296 | 0.0259 % | 0.00862300263° | 19 | 9.999935 | 0 of 1 541 174 | 7 runs, 1 distinct | dac84b2b640a |
| Numedalslagen | 20c-2 65ae0792, --start-quality-gain -1 | 1 290 807 | 433 | 0.0335 % | 0.000680285394° | 19 | 9.999972 | 0 of 1 783 018 | 7 runs, 1 distinct | 93f93b681dd2 |
| Numedalslagen | 20c-2 65ae0792, --no-constraint-feet | 1 085 653 | 5507 | 0.5073 % | 7.22481777e-06° | 19 | 9.999934 | 0 of 1 479 600 | 7 runs, 1 distinct | 4b432d60dc4e |

### Quality-start and refinement counts (`--stats`; the same in every repeat, else listed)

| catchment | build | start: points added | start: tries skipped (all reasons) | start: of them, without gain (R7) | start: points snapped to lines (feet) | start: lines split (R8) | refinement: points added |
|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | master 41bda81a (20c-1) | 503 | 252 | - | 0 | - | 219 837 |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, gain 0 (default) | 503 | 252 | 0 | 0 | 0 | 219 837 |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, --start-quality-gain -1 | 503 | 252 | 0 | 0 | 0 | 219 837 |
| 1 m tile (tolerance 1) | 20c-2 65ae0792, --no-constraint-feet | 503 | 252 | 0 | 0 | 0 | 219 837 |
| Lagan | master 41bda81a (20c-1) | 192 169 | 104 936 | - | 5 930 | - | 46 144 |
| Lagan | 20c-2 65ae0792, gain 0 (default) | 136 273 | 111 408 | 34 600 | 3 488 | 14 025 | 54 097 |
| Lagan | 20c-2 65ae0792, --start-quality-gain -1 | 192 169 | 104 936 | 0 | 5 930 | 0 | 46 144 |
| Lagan | 20c-2 65ae0792, --no-constraint-feet | 126 404 | 123 780 | 35 613 | 0 | 0 | 56 735 |
| Numedalslagen | master 41bda81a (20c-1) | 314 843 | 75 663 | - | 4 233 | - | 183 368 |
| Numedalslagen | 20c-2 65ae0792, gain 0 (default) | 208 852 | 88 243 | 71 590 | 2 565 | 20 579 | 195 799 |
| Numedalslagen | 20c-2 65ae0792, --start-quality-gain -1 | 314 843 | 75 663 | 0 | 4 233 | 0 | 183 368 |
| Numedalslagen | 20c-2 65ae0792, --no-constraint-feet | 189 036 | 104 679 | 66 173 | 0 | 0 | 205 896 |

### 2. Phases: median seconds over repeats r1 and on ([6] runs per build), threads 10 (the CLI default)

| catchment | phase | master 41bda81a (20c-1): median (min to max) | 20c-2 65ae0792, gain 0 (default): median (min to max) | 20c-2 65ae0792, --start-quality-gain -1: median (min to max) | 20c-2 65ae0792, --no-constraint-feet: median (min to max) | gain 0 vs -1 | gain 0 vs 20c-1 | gain -1 vs 20c-1 | gain 0 vs feet off |
|---|---|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | refine | 0.1914 (0.1866 to 0.1934) | 0.1864 (0.1856 to 0.1904) | 0.1905 (0.1873 to 0.1930) | 0.1912 (0.1859 to 0.1925) | -2.18 % | -2.65 % | -0.48 % | -2.54 % |
| 1 m tile (tolerance 1) | refine: start quality | 0.0006 (0.0006 to 0.0007) | 0.0007 (0.0007 to 0.0008) | 0.0006 (0.0006 to 0.0007) | 0.0007 (0.0007 to 0.0008) | +14.87 % | +18.56 % | +3.21 % | +1.56 % |
| 1 m tile (tolerance 1) | refine: scan (parallel) | 0.0569 (0.0566 to 0.0573) | 0.0567 (0.0564 to 0.0575) | 0.0568 (0.0566 to 0.0577) | 0.0569 (0.0567 to 0.0579) | -0.29 % | -0.37 % | -0.08 % | -0.36 % |
| 1 m tile (tolerance 1) | refine: split + flip (serial) | 0.1080 (0.1034 to 0.1093) | 0.1032 (0.1028 to 0.1065) | 0.1070 (0.1041 to 0.1090) | 0.1077 (0.1027 to 0.1078) | -3.54 % | -4.52 % | -1.01 % | -4.24 % |
| 1 m tile (tolerance 1) | total | 0.4240 (0.4200 to 0.4250) | 0.4190 (0.4140 to 0.4220) | 0.4230 (0.4180 to 0.4330) | 0.4235 (0.4180 to 0.4280) | -0.95 % | -1.18 % | -0.24 % | -1.06 % |
| Lagan | refine | 0.5682 (0.5656 to 0.5742) | 0.6593 (0.6569 to 0.6634) | 0.5759 (0.5743 to 0.5967) | 0.6042 (0.6031 to 0.6048) | +14.48 % | +16.05 % | +1.37 % | +9.13 % |
| Lagan | refine: start quality | 0.3491 (0.3453 to 0.3530) | 0.4310 (0.4286 to 0.4369) | 0.3534 (0.3505 to 0.3712) | 0.3765 (0.3745 to 0.3778) | +21.96 % | +23.46 % | +1.23 % | +14.48 % |
| Lagan | refine: scan (parallel) | 0.0501 (0.0497 to 0.0508) | 0.0518 (0.0514 to 0.0521) | 0.0502 (0.0500 to 0.0509) | 0.0513 (0.0507 to 0.0516) | +3.15 % | +3.42 % | +0.26 % | +1.02 % |
| Lagan | refine: split + flip (serial) | 0.0378 (0.0369 to 0.0389) | 0.0440 (0.0431 to 0.0452) | 0.0382 (0.0367 to 0.0391) | 0.0421 (0.0409 to 0.0427) | +15.00 % | +16.29 % | +1.12 % | +4.53 % |
| Lagan | total | 72.3525 (71.9990 to 72.8540) | 72.3210 (72.0380 to 72.5300) | 72.0720 (71.5300 to 72.4190) | 72.0465 (71.5810 to 72.7420) | +0.35 % | -0.04 % | -0.39 % | +0.38 % |
| Numedalslagen | refine | 1.1151 (1.1114 to 1.1204) | 1.2398 (1.2343 to 1.2535) | 1.1179 (1.1147 to 1.1360) | 1.1307 (1.1236 to 1.1403) | +10.91 % | +11.19 % | +0.25 % | +9.65 % |
| Numedalslagen | refine: start quality | 0.5348 (0.5304 to 0.5417) | 0.6423 (0.6415 to 0.6562) | 0.5355 (0.5335 to 0.5379) | 0.5492 (0.5444 to 0.5569) | +19.96 % | +20.11 % | +0.12 % | +16.96 % |
| Numedalslagen | refine: scan (parallel) | 0.2847 (0.2825 to 0.2866) | 0.2982 (0.2976 to 0.2996) | 0.2861 (0.2841 to 0.3012) | 0.2992 (0.2982 to 0.3024) | +4.25 % | +4.73 % | +0.46 % | -0.32 % |
| Numedalslagen | refine: split + flip (serial) | 0.1548 (0.1522 to 0.1575) | 0.1610 (0.1591 to 0.1630) | 0.1554 (0.1545 to 0.1601) | 0.1475 (0.1468 to 0.1494) | +3.60 % | +4.00 % | +0.39 % | +9.14 % |
| Numedalslagen | total | 6.6360 (6.5940 to 6.6590) | 6.6305 (6.5980 to 6.6540) | 6.6510 (6.6060 to 6.7500) | 6.4175 (6.3910 to 6.4760) | -0.31 % | -0.08 % | +0.23 % | +3.32 % |

### 2. The longer tile series (`run.sh tile 0 .. 20`, `raw-tile/`; r0 a warm-up; 20, 20, 20 timed runs, base, head and feet off in turn)

| phase | master 41bda81a (20c-1): median (min to max) | 20c-2 65ae0792, gain 0 (default): median (min to max) | 20c-2 65ae0792, --no-constraint-feet: median (min to max) | gain 0 vs 20c-1 | gain 0 vs feet off |
|---|---|---|---|---|---|
| refine | 0.1915 (0.1866 to 0.2056) | 0.1932 (0.1869 to 0.2045) | 0.1926 (0.1867 to 0.2378) | +0.87 % | +0.30 % |
| refine: start quality | 0.0007 (0.0006 to 0.0008) | 0.0008 (0.0007 to 0.0009) | 0.0008 (0.0007 to 0.0009) | +16.20 % | +1.83 % |
| refine: scan (parallel) | 0.0573 (0.0566 to 0.0702) | 0.0580 (0.0567 to 0.0641) | 0.0574 (0.0567 to 0.1011) | +1.35 % | +1.03 % |
| refine: split + flip (serial) | 0.1058 (0.1004 to 0.1099) | 0.1073 (0.1021 to 0.1144) | 0.1076 (0.0993 to 0.1099) | +1.39 % | -0.34 % |
| total | 0.4285 (0.4150 to 0.4480) | 0.4380 (0.4220 to 0.5080) | 0.4305 (0.4140 to 0.4770) | +2.22 % | +1.74 % |

Every mesh in the series equal to every other: True.

### 4. The split-phase limit (`docs/increments/20c-soft-quality.md`, "The split-phase limit, judged after the fix")

| point | measure | limit | measured (medians) | result |
|---|---|---|---|---|
| 1. tile, no cost where no foot is placed (raw/, 6 runs) | split phase, gain 0 vs 20c-2 65ae0792, --no-constraint-feet | ≤ +2 % | -4.24 % | met |
| 1. tile, no cost where no foot is placed (raw/, 6 runs) | split phase, gain 0 vs master 41bda81a (20c-1) | ≤ +2 % | -4.52 % | met |
| 1. tile, no cost where no foot is placed (raw-tile/, 20 runs) | split phase, gain 0 vs 20c-2 65ae0792, --no-constraint-feet | ≤ +2 % | -0.34 % | met |
| 1. tile, no cost where no foot is placed (raw-tile/, 20 runs) | split phase, gain 0 vs master 41bda81a (20c-1) | ≤ +2 % | +1.39 % | met |
| 2. Lagan, no slower than what Ola runs today | refine: split + flip (serial), gain 0 vs 20c-1 | ≤ +2 % | +16.29 % | NOT MET |
| 2. Lagan, no slower than what Ola runs today | total, gain 0 vs 20c-1 | within 2 % | -0.04 % | met |
| 3. Lagan, reported, not gated | refine: start quality, gain 0 vs 20c-1 / vs feet off | - | +23.46 % / +14.48 % | - |
| 3. Lagan, reported, not gated | refine, gain 0 vs 20c-1 / vs feet off | - | +16.05 % / +9.13 % | - |
| 2. Numedalslagen, no slower than what Ola runs today | refine: split + flip (serial), gain 0 vs 20c-1 | ≤ +2 % | +4.00 % | NOT MET |
| 2. Numedalslagen, no slower than what Ola runs today | total, gain 0 vs 20c-1 | within 2 % | -0.08 % | met |
| 3. Numedalslagen, reported, not gated | refine: start quality, gain 0 vs 20c-1 / vs feet off | - | +20.11 % / +16.96 % | - |
| 3. Numedalslagen, reported, not gated | refine, gain 0 vs 20c-1 / vs feet off | - | +11.19 % / +9.65 % | - |

### The split phase per split (for reading point 2: the meshes differ, so the work differs)

Splits are `points_inserted` from `--stats` (refinement's main loop). µs per split = median split-phase
seconds / splits. Reported beside the limit, which is written in seconds; not a gate.

| catchment | splits, base | splits, head | splits, headm1 | splits, nofeet | µs per split, base | µs per split, head | µs per split, headm1 | µs per split, nofeet | splits, gain 0 vs 20c-1 | per split, gain 0 vs 20c-1 | per split, gain 0 vs feet off |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | 219 837 | 219 837 | 219 837 | 219 837 | 0.4915 | 0.4693 | 0.4865 | 0.4901 | +0.00 % | -4.52 % | -4.24 % |
| Lagan | 46 144 | 54 097 | 46 144 | 56 735 | 0.8194 | 0.8128 | 0.8286 | 0.7414 | +17.24 % | -0.81 % | +9.63 % |
| Numedalslagen | 183 368 | 195 799 | 183 368 | 205 896 | 0.8441 | 0.8222 | 0.8474 | 0.7164 | +6.78 % | -2.60 % | +14.77 % |

### Before and after the fix, side by side (gain 0; `../raw` is another session, so read the 20c-1 column for drift)

| catchment | phase | 20c-1, first session | 600b2acb | 600b2acb vs 20c-1 | 20c-1, this session | 65ae0792 | 65ae0792 vs 20c-1 | 65ae0792 vs 600b2acb |
|---|---|---|---|---|---|---|---|---|
| 1 m tile (tolerance 1) | refine: start quality | 0.0007 | 0.0010 | +43.75 % | 0.0006 | 0.0007 | +18.56 % | -27.82 % |
| 1 m tile (tolerance 1) | refine: split + flip (serial) | 0.1067 | 0.1103 | +3.36 % | 0.1080 | 0.1032 | -4.52 % | -6.44 % |
| 1 m tile (tolerance 1) | refine | 0.1909 | 0.1953 | +2.31 % | 0.1914 | 0.1864 | -2.65 % | -4.56 % |
| 1 m tile (tolerance 1) | total | 0.4320 | 0.4340 | +0.46 % | 0.4240 | 0.4190 | -1.18 % | -3.46 % |
| Lagan | refine: start quality | 0.3549 | 0.5738 | +61.66 % | 0.3491 | 0.4310 | +23.46 % | -24.89 % |
| Lagan | refine: split + flip (serial) | 0.0382 | 0.0455 | +18.96 % | 0.0378 | 0.0440 | +16.29 % | -3.31 % |
| Lagan | refine | 0.5759 | 0.8040 | +39.60 % | 0.5682 | 0.6593 | +16.05 % | -17.99 % |
| Lagan | total | 72.7235 | 72.1775 | -0.75 % | 72.3525 | 72.3210 | -0.04 % | +0.20 % |
| Numedalslagen | refine: start quality | 0.5389 | 0.8591 | +59.40 % | 0.5348 | 0.6423 | +20.11 % | -25.23 % |
| Numedalslagen | refine: split + flip (serial) | 0.1559 | 0.1625 | +4.22 % | 0.1548 | 0.1610 | +4.00 % | -0.92 % |
| Numedalslagen | refine | 1.1241 | 1.4589 | +29.78 % | 1.1151 | 1.2398 | +11.19 % | -15.01 % |
| Numedalslagen | total | 6.6755 | 6.8910 | +3.23 % | 6.6360 | 6.6305 | -0.08 % | -3.78 % |

### Spread of 20c-1's own runs (base, (max - min) / median)

| catchment | refine: start quality | refine: split + flip (serial) | refine | total |
|---|---|---|---|---|
| 1 m tile (tolerance 1) | 10.3 % | 5.5 % | 3.6 % | 1.2 % |
| Lagan | 2.2 % | 5.3 % | 1.5 % | 1.2 % |
| Numedalslagen | 2.1 % | 3.4 % | 0.8 % | 1.0 % |
| 1 m tile (tolerance 1), longer series | 22.9 % | 9.1 % | 9.9 % | 7.7 % |

### Power, every run (`pmset -g batt` before and after; 1 m benchmark: `run.json`)

| run | power | battery before -> after | process s |
|---|---|---|---|
| lagan-base-r0 | AC | 100 % charged -> 100 % charged | 73.5 |
| lagan-base-r1 | AC | 100 % charged -> 100 % charged | 73.4 |
| lagan-base-r2 | AC | 100 % charged -> 100 % charged | 73.6 |
| lagan-base-r3 | AC | 100 % charged -> 100 % charged | 73.8 |
| lagan-base-r4 | AC | 100 % charged -> 100 % charged | 73.3 |
| lagan-base-r5 | AC | 100 % charged -> 100 % charged | 74.2 |
| lagan-base-r6 | AC | 100 % charged -> 100 % charged | 73.7 |
| lagan-head-r0 | AC | 100 % charged -> 100 % charged | 73.7 |
| lagan-head-r1 | AC | 100 % charged -> 100 % charged | 73.7 |
| lagan-head-r2 | AC | 100 % charged -> 100 % charged | 73.3 |
| lagan-head-r3 | AC | 100 % charged -> 100 % charged | 73.6 |
| lagan-head-r4 | AC | 100 % charged -> 100 % charged | 73.3 |
| lagan-head-r5 | AC | 100 % charged -> 100 % charged | 73.6 |
| lagan-head-r6 | AC | 100 % charged -> 100 % charged | 73.5 |
| lagan-headm1-r0 | AC | 100 % charged -> 100 % charged | 72.9 |
| lagan-headm1-r1 | AC | 100 % charged -> 100 % charged | 73.1 |
| lagan-headm1-r2 | AC | 100 % charged -> 100 % charged | 73.5 |
| lagan-headm1-r3 | AC | 100 % charged -> 100 % charged | 72.9 |
| lagan-headm1-r4 | AC | 100 % charged -> 100 % charged | 73.5 |
| lagan-headm1-r5 | AC | 100 % charged -> 100 % charged | 73.3 |
| lagan-headm1-r6 | AC | 100 % charged -> 100 % charged | 73.7 |
| lagan-nofeet-r0 | AC | 100 % charged -> 100 % charged | 73.0 |
| lagan-nofeet-r1 | AC | 100 % charged -> 100 % charged | 73.9 |
| lagan-nofeet-r2 | AC | 100 % charged -> 100 % charged | 73.2 |
| lagan-nofeet-r3 | AC | 100 % charged -> 100 % charged | 73.2 |
| lagan-nofeet-r4 | AC | 100 % charged -> 100 % charged | 73.0 |
| lagan-nofeet-r5 | AC | 100 % charged -> 100 % charged | 73.4 |
| lagan-nofeet-r6 | AC | 100 % charged -> 100 % charged | 72.7 |
| numed-base-r0 | AC | 100 % charged -> 100 % charged | 8.4 |
| numed-base-r1 | AC | 100 % charged -> 100 % charged | 8.4 |
| numed-base-r2 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-base-r3 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-base-r4 | AC | 100 % charged -> 100 % charged | 8.4 |
| numed-base-r5 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-base-r6 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-head-r0 | AC | 100 % charged -> 100 % charged | 8.2 |
| numed-head-r1 | AC | 100 % charged -> 100 % charged | 8.2 |
| numed-head-r2 | AC | 100 % charged -> 100 % charged | 8.3 |
| numed-head-r3 | AC | 100 % charged -> 100 % charged | 8.2 |
| numed-head-r4 | AC | 100 % charged -> 100 % charged | 8.3 |
| numed-head-r5 | AC | 100 % charged -> 100 % charged | 8.2 |
| numed-head-r6 | AC | 100 % charged -> 100 % charged | 8.3 |
| numed-headm1-r0 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-headm1-r1 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-headm1-r2 | AC | 100 % charged -> 100 % charged | 8.6 |
| numed-headm1-r3 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-headm1-r4 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-headm1-r5 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-headm1-r6 | AC | 100 % charged -> 100 % charged | 8.5 |
| numed-nofeet-r0 | AC | 100 % charged -> 100 % charged | 8.3 |
| numed-nofeet-r1 | AC | 100 % charged -> 100 % charged | 7.9 |
| numed-nofeet-r2 | AC | 100 % charged -> 100 % charged | 7.9 |
| numed-nofeet-r3 | AC | 100 % charged -> 100 % charged | 7.9 |
| numed-nofeet-r4 | AC | 100 % charged -> 100 % charged | 8.0 |
| numed-nofeet-r5 | AC | 100 % charged -> 100 % charged | 7.9 |
| numed-nofeet-r6 | AC | 100 % charged -> 100 % charged | 8.0 |
| tile-base-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r1 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-base-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r0 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-head-r1 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r2 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-head-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r4 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-head-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-head-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r1 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-headm1-r2 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-headm1-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-headm1-r4 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-headm1-r5 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-headm1-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r0 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-nofeet-r1 | AC | 100 % charged -> 100 % charged | 2.6 |
| tile-nofeet-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| tile-nofeet-r6 | AC | 100 % charged -> 100 % charged | 2.6 |
| raw-tile/tile-base-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r1 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r5 | AC | 100 % charged -> 100 % charged | 2.6 |
| raw-tile/tile-base-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r7 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r8 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r9 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r10 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r11 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r12 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r13 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r14 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r15 | AC | 100 % charged -> 100 % charged | 2.8 |
| raw-tile/tile-base-r16 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r17 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r18 | AC | 100 % charged -> 100 % charged | 2.9 |
| raw-tile/tile-base-r19 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-base-r20 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r1 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r7 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r8 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r9 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r10 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r11 | AC | 100 % charged -> 100 % charged | 2.6 |
| raw-tile/tile-head-r12 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r13 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r14 | AC | 100 % charged -> 100 % charged | 2.8 |
| raw-tile/tile-head-r15 | AC | 100 % charged -> 100 % charged | 2.8 |
| raw-tile/tile-head-r16 | AC | 100 % charged -> 100 % charged | 2.8 |
| raw-tile/tile-head-r17 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r18 | AC | 100 % charged -> 100 % charged | 2.8 |
| raw-tile/tile-head-r19 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-head-r20 | AC | 100 % charged -> 100 % charged | 2.8 |
| raw-tile/tile-nofeet-r0 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r1 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r2 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r3 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r4 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r5 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r6 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r7 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r8 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r9 | AC | 100 % charged -> 100 % charged | 2.6 |
| raw-tile/tile-nofeet-r10 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r11 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r12 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r13 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r14 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r15 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r16 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r17 | AC | 100 % charged -> 100 % charged | 2.8 |
| raw-tile/tile-nofeet-r18 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r19 | AC | 100 % charged -> 100 % charged | 2.7 |
| raw-tile/tile-nofeet-r20 | AC | 100 % charged -> 100 % charged | 2.7 |
| b1m-base-r1 | ac | 100 % charged -> 100 % charged | - |
| b1m-base-r2 | ac | 100 % charged -> 100 % charged | - |
| b1m-base-r3 | ac | 100 % charged -> 100 % charged | - |
| b1m-head-r1 | ac | 100 % charged -> 100 % charged | - |
| b1m-head-r2 | ac | 100 % charged -> 100 % charged | - |
| b1m-head-r3 | ac | 100 % charged -> 100 % charged | - |

### Where each run loaded tin_engine and _core from (first two lines of every `raw/*.log` and `raw-tile/*.log`)

| runs | count | tin_engine.__file__ | _core.__file__ |
|---|---|---|---|
| raw/lagan-base | 7 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b/base/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b/base/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/lagan-head | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/lagan-headm1 | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/lagan-nofeet | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/numed-base | 7 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b/base/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b/base/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/numed-head | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/numed-headm1 | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/numed-nofeet | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/tile-base | 7 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b/base/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b/base/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/tile-head | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/tile-headm1 | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw/tile-nofeet | 7 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw-tile/tile-base | 21 | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b/base/build-bench/pkg/tin_engine/__init__.py` | `/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b/base/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw-tile/tile-head | 21 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |
| raw-tile/tile-nofeet | 21 | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/__init__.py` | `/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2/build-bench/pkg/tin_engine/_core.cpython-314-darwin.so` |

### 3. The 1 m benchmark and thread sweep (`tools/bench.py run`, `bench/b1m-*`; `bench_summarize.py bench b1m`)

base = master 41bda81a, head = 20c-2 65ae0792.

base: 3 runs, commits ['41bda81'], power ['ac'], hardening ['libc++ fast']
head: 3 runs, commits ['65ae079'], power ['ac'], hardening ['libc++ fast']

| domain | mesh sha256 (base) | mesh sha256 (head) | equal |
|---|---|---|---|
| tile | 11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c | 11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c | True |
| quarter | a60fb597b5b738268a3bd27730989392ee86fb030fe060b78fd3a0bcdbb60a34 | a60fb597b5b738268a3bd27730989392ee86fb030fe060b78fd3a0bcdbb60a34 | True |

| domain | threads | refine_s base | refine_s head | change | app_s change | proc_s change |
|---|---:|---:|---:|---:|---:|---:|
| tile | default | 0.1924 | 0.1917 | -0.3 % | -0.4 % | -0.2 % |
| tile | 1 | 0.4843 | 0.4837 | -0.1 % | -0.0 % | -0.1 % |
| tile | 2 | 0.3151 | 0.3151 | -0.0 % | -0.0 % | -0.1 % |
| tile | 3 | 0.2602 | 0.2612 | +0.4 % | -0.1 % | +0.0 % |
| tile | 4 | 0.2353 | 0.2363 | +0.4 % | -0.6 % | -0.4 % |
| tile | 5 | 0.2177 | 0.2181 | +0.2 % | -0.1 % | +0.0 % |
| tile | 6 | 0.2056 | 0.2063 | +0.4 % | -0.4 % | -0.1 % |
| tile | 7 | 0.1961 | 0.1967 | +0.3 % | +0.1 % | +0.3 % |
| tile | 8 | 0.1913 | 0.1917 | +0.2 % | -0.7 % | -0.3 % |
| tile | 9 | 0.1938 | 0.1935 | -0.2 % | +0.1 % | -0.1 % |
| tile | 10 | 0.1920 | 0.1925 | +0.3 % | +0.1 % | -0.3 % |
| tile | 11 | 0.1926 | 0.1924 | -0.1 % | +0.2 % | +0.2 % |
| tile | 12 | 0.1917 | 0.1931 | +0.8 % | +0.2 % | +0.0 % |
| tile | 13 | 0.1927 | 0.1928 | +0.0 % | +0.6 % | -0.4 % |
| tile | 14 | 0.1926 | 0.1920 | -0.3 % | +0.0 % | +0.2 % |
| tile | 15 | 0.1926 | 0.1910 | -0.8 % | +0.4 % | +0.2 % |
| tile | 16 | 0.1935 | 0.1938 | +0.2 % | -0.0 % | +0.1 % |
| tile | 17 | 0.1929 | 0.1928 | -0.1 % | +0.4 % | +0.4 % |
| tile | 18 | 0.1932 | 0.1927 | -0.3 % | +0.4 % | -0.1 % |
| tile | 19 | 0.1938 | 0.1930 | -0.4 % | -0.0 % | +0.1 % |
| tile | 20 | 0.1941 | 0.1945 | +0.2 % | -0.2 % | +0.0 % |
| quarter | default | 0.1731 | 0.1740 | +0.5 % | +0.3 % | +0.4 % |
| quarter | 1 | 0.4279 | 0.4280 | +0.0 % | -0.2 % | +0.1 % |
| quarter | 2 | 0.2778 | 0.2785 | +0.3 % | +0.7 % | +0.4 % |
| quarter | 3 | 0.2280 | 0.2306 | +1.1 % | +0.8 % | +0.7 % |
| quarter | 4 | 0.2076 | 0.2089 | +0.6 % | +0.3 % | +0.3 % |
| quarter | 5 | 0.1921 | 0.1932 | +0.6 % | +0.7 % | +0.3 % |
| quarter | 6 | 0.1814 | 0.1824 | +0.5 % | +0.7 % | +0.5 % |
| quarter | 7 | 0.1744 | 0.1750 | +0.3 % | -0.1 % | -0.2 % |
| quarter | 8 | 0.1703 | 0.1707 | +0.3 % | +0.5 % | -0.0 % |
| quarter | 9 | 0.1732 | 0.1742 | +0.6 % | +0.3 % | +0.8 % |
| quarter | 10 | 0.1724 | 0.1736 | +0.7 % | +0.6 % | +0.5 % |
| quarter | 11 | 0.1725 | 0.1737 | +0.7 % | +1.0 % | +0.7 % |
| quarter | 12 | 0.1738 | 0.1737 | -0.0 % | -0.0 % | +0.5 % |
| quarter | 13 | 0.1735 | 0.1733 | -0.1 % | +0.5 % | +0.1 % |
| quarter | 14 | 0.1736 | 0.1736 | -0.0 % | +0.6 % | +0.6 % |
| quarter | 15 | 0.1729 | 0.1746 | +1.0 % | +0.7 % | +0.4 % |
| quarter | 16 | 0.1723 | 0.1745 | +1.3 % | +0.9 % | +0.7 % |
| quarter | 17 | 0.1714 | 0.1742 | +1.6 % | +0.7 % | +0.2 % |
| quarter | 18 | 0.1731 | 0.1752 | +1.2 % | +0.9 % | +0.3 % |
| quarter | 19 | 0.1703 | 0.1747 | +2.6 % | +0.7 % | +0.4 % |
| quarter | 20 | 0.1741 | 0.1751 | +0.6 % | +0.7 % | +0.3 % |

refine_s pooled change over every cell: -0.8 .. +2.6 %
base-vs-base noise (per cell, max/min of the base runs' medians): median 0.8 %, max 2.7 %

| domain | side | worst angle | median angle | share < 1° | max degree | within tolerance | Delaunay violations (checked) |
|---|---|---|---|---|---|---|---|
| tile | base (1 distinct of 3) | 0.6296° | 45.00° | 0.0002 % | 74 | True | 0 (692056) |
| tile | head (1 distinct of 3) | 0.6296° | 45.00° | 0.0002 % | 74 | True | 0 (692056) |
| quarter | base (1 distinct of 3) | 0.3955° | 45.00° | 0.0028 % | 18 | True | 0 (641792) |
| quarter | head (1 distinct of 3) | 0.3955° | 45.00° | 0.0028 % | 18 | True | 0 (641792) |

| domain | side | refine_s 1 thread | refine_s 20 threads | speed-up |
|---|---|---:|---:|---:|
| tile | base | 0.4843 | 0.1941 | 2.49 x |
| tile | head | 0.4837 | 0.1945 | 2.49 x |
| quarter | base | 0.4279 | 0.1741 | 2.46 x |
| quarter | head | 0.4280 | 0.1751 | 2.44 x |
<!-- end of generated tables -->
