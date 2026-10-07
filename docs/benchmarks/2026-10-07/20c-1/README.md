# 20c-1 acceptance: the foot rule on every insertion path (@perf, 2026-10-07)

Branch `worktree-soft-quality`, code at `34b546c4`, against master `ed125121`.
The rules: `docs/increments/20c-soft-quality.md`, "PRs and gates" (the 20c-1
row), the paragraph after it on how `@perf` runs the catchments, and "Rulings
on 20c-1's green step" item 5 (`foot_epsilon`: the split phase's seconds,
feet on against `--no-constraint-feet`, with a 2 % limit); and
`docs/increments/README.md`, "Acceptance" (the 1 m benchmark and thread
sweep from `tools/bench.py`).

## Verdict

**REGRESSION: the time measures. Every mesh and quality gate passes.**

- **Quality gates, both catchments: all pass** (table "The 20c-1 gates").
  Lagan and Numedalslagen reproduce the design's M5 figures for 20c-1
  exactly: triangle count, count under 1°, and worst angle.
- **`--no-constraint-feet`**: bit-identical to master run with
  `--no-constraint-feet`, on both catchments. It is **not** identical to
  master's default output. Master's default already has 20b's refinement
  feet, and 20c-1's switch turns those off as well. The gate was read as
  "master with the same switch" (see "Reading of the bit-identity gate").
- **The split phase, feet on against feet off (ruling 5): over the 2 %
  limit on both catchments.** The `refine: split + flip (serial)` row is
  higher with feet on, by the percentages in the table "Feet on against
  feet off". Under ruling 5, `@developer` therefore computes ε only where a
  constrained edge is near. Against master's default, the same phase is
  also higher.
- **1 m benchmark: REGRESSION on refine time, with identical meshes.** The
  tile's and the quarter's meshes are byte-identical to master's (same
  `mesh_sha256`). Refine is slower at every thread count, and the gap grows
  with the thread count. The 1-to-20-thread speed-up falls (table at the
  end). The per-phase run on the 1 m tile (`raw1m/`, rows "1 m tile")
  puts the extra time in `refine: split + flip (serial)`. The parallel scan
  does not change. In that run no foot is placed and the mesh is the same.
  The cause inside the split phase is **thought to be** the `foot_epsilon`
  work that ruling 5 names. That has not been profiled, so it is a
  hypothesis.

## Reading of the bit-identity gate

The design (R2.6) says `--no-constraint-feet` "turns all three paths off and
the mesh is bit-identical to master's". On master, `--no-constraint-feet`
already exists and turns 20b's refinement feet off, and master's default
has them on. 20c-1's switch also turns refinement's feet off. So the only
comparison that can hold is with master under the same switch, which is
the reading the existing test takes:
`tests/python/test_cli_constraint_feet.py::test_no_constraint_feet_matches_increment_20`
compares the switch's output with digests recorded before 20b. Both
readings are in the tables.

## Method

- **Machine and power.** Apple M1 Max, macOS 27.0, 10 hardware threads. AC power
  for every run: `pmset -g batt` was read before and after each run and is
  stored in `raw/*.power`, summarised in "Power, every run". The scripts
  stop before any run that would start off AC, and set aside
  (`*.DISCARD`) any run whose power source changed. No run was set aside.
  The battery stayed at 100 %, charged.
- **Builds.** `tools/bench.py`'s `build()` produced both: Release `-O3
  -DNDEBUG`, AppleClang 21, `RASPUTIN_HARDENING=ON` (`libc++ fast` in every
  run). Each build is in `<tree>/build-bench/pkg`. Base: a scratch worktree
  of `ed125121` (`_core` sha256 `19d71363…`). Head: this worktree
  (`8940a78a…`). Nothing was installed into a venv. Every run loads the
  package from its `pkg` tree and prints `tin_engine.__file__` and
  `_core.__file__` (first lines of each `raw/*.log`), and `drive.py`
  refuses to run when either file comes from elsewhere. The Python used is
  this worktree's `.venv` (Python 3.14) for both trees.
  `pyproject.toml` does not differ between the two commits.
- **Catchments** (`run.sh`, `drive.py`; log `runs.log`). The commands are
  the ones in `../rasputin_scratch/20c-prototype/runs.sh` that set M5's
  figures, with nothing added. Lagan: GLO-30 from the local cache,
  `--out-crs EPSG:3006`, SMHI outline, CORINE 2018 (EPSG:3035),
  `--tolerance 10`. Numedalslagen: DTM10, NVE outline, CORINE 2018 gpkg,
  `--tolerance 10`. Threads: the CLI default (10). There were four builds:
  master, 20c-1, and each of them with `--no-constraint-feet`.
  There were 8 repeats (r0 to r7). The order was reversed on every other
  repeat. r0 is a warm-up and is left out of the timing medians. Exact
  counts and the worst angle come from the trimmed output arrays at full
  precision; the `--stats` file rounds the share to 0.01 %. The phase times
  are the `--stats` phases, unrounded (`PhaseClock.add`).
  The Delaunay check is `tools/bench.py`'s constrained check.
- **1 m tile phases** (`onem.sh`, `raw1m/`): the same four builds on
  `tests/fixtures/dem_archive/7908_3_10m_z33.tif`, tolerance 1, 7 repeats,
  interleaved.
- **1 m benchmark and thread sweep** (`pairs.sh`; log `bench/pairs.log`):
  `tools/bench.py run` with its defaults: the tile and
  `docs/benchmarks/2026-09-26/quarter.geojson`, tolerance 1, threads default
  and 1 to 20, 5 repeats, hardening on. The order was base, head, head,
  base. Each `run.json` records its tree commit, its `.so` hash, the
  `bench.py` blob and the power state. All of these runs were on AC. The
  verdict lines `bench.py` printed (in `pairs.log`) compare each run with
  the previous stored run. The pooled table below, from
  `bench_summarize.py` (copied from `2026-10-06/audit-pr-a/`), is the one
  that counts.
- **Tables** below are generated by `python3 summarize.py`, which reads only
  `raw/`, `raw1m/` and `bench/` and rewrites the block between the markers.
- **Meshes** were in the session scratchpad and are not kept. To regenerate
  them, rerun `run.sh 0 1 2 3 4 5 6 7`, `onem.sh` and `pairs.sh`. The paths
  are inside the scripts. The base tree is a worktree of `ed125121`, built
  with `tools/bench.py`'s `build()`.

<!-- tables: generated by summarize.py -->

### Power, every run

| run | power before and after | battery before -> after | process s |
|---|---|---|---|
| lagan-base-r0 | AC | 100 % charged -> 100 % charged | 82.7 |
| lagan-base-r1 | AC | 100 % charged -> 100 % charged | 82.5 |
| lagan-base-r2 | AC | 100 % charged -> 100 % charged | 82.3 |
| lagan-base-r3 | AC | 100 % charged -> 100 % charged | 82.0 |
| lagan-base-r4 | AC | 100 % charged -> 100 % charged | 83.3 |
| lagan-base-r5 | AC | 100 % charged -> 100 % charged | 83.0 |
| lagan-base-r6 | AC | 100 % charged -> 100 % charged | 82.4 |
| lagan-base-r7 | AC | 100 % charged -> 100 % charged | 83.1 |
| lagan-basenofeet-r0 | AC | 100 % charged -> 100 % charged | 82.3 |
| lagan-basenofeet-r1 | AC | 100 % charged -> 100 % charged | 83.7 |
| lagan-basenofeet-r2 | AC | 100 % charged -> 100 % charged | 82.9 |
| lagan-basenofeet-r3 | AC | 100 % charged -> 100 % charged | 82.4 |
| lagan-basenofeet-r4 | AC | 100 % charged -> 100 % charged | 82.8 |
| lagan-basenofeet-r5 | AC | 100 % charged -> 100 % charged | 81.9 |
| lagan-basenofeet-r6 | AC | 100 % charged -> 100 % charged | 82.8 |
| lagan-basenofeet-r7 | AC | 100 % charged -> 100 % charged | 83.0 |
| lagan-head-r0 | AC | 100 % charged -> 100 % charged | 82.7 |
| lagan-head-r1 | AC | 100 % charged -> 100 % charged | 82.6 |
| lagan-head-r2 | AC | 100 % charged -> 100 % charged | 82.9 |
| lagan-head-r3 | AC | 100 % charged -> 100 % charged | 83.3 |
| lagan-head-r4 | AC | 100 % charged -> 100 % charged | 82.9 |
| lagan-head-r5 | AC | 100 % charged -> 100 % charged | 82.9 |
| lagan-head-r6 | AC | 100 % charged -> 100 % charged | 82.1 |
| lagan-head-r7 | AC | 100 % charged -> 100 % charged | 82.9 |
| lagan-nofeet-r0 | AC | 100 % charged -> 100 % charged | 82.0 |
| lagan-nofeet-r1 | AC | 100 % charged -> 100 % charged | 82.3 |
| lagan-nofeet-r2 | AC | 100 % charged -> 100 % charged | 83.1 |
| lagan-nofeet-r3 | AC | 100 % charged -> 100 % charged | 82.9 |
| lagan-nofeet-r4 | AC | 100 % charged -> 100 % charged | 83.1 |
| lagan-nofeet-r5 | AC | 100 % charged -> 100 % charged | 83.1 |
| lagan-nofeet-r6 | AC | 100 % charged -> 100 % charged | 82.7 |
| lagan-nofeet-r7 | AC | 100 % charged -> 100 % charged | 82.4 |
| numed-base-r0 | AC | 100 % charged -> 100 % charged | 13.0 |
| numed-base-r1 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-base-r2 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-base-r3 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-base-r4 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-base-r5 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-base-r6 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-base-r7 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-basenofeet-r0 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-basenofeet-r1 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-basenofeet-r2 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-basenofeet-r3 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-basenofeet-r4 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-basenofeet-r5 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-basenofeet-r6 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-basenofeet-r7 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-head-r0 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-head-r1 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-head-r2 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-head-r3 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-head-r4 | AC | 100 % charged -> 100 % charged | 13.3 |
| numed-head-r5 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-head-r6 | AC | 100 % charged -> 100 % charged | 13.3 |
| numed-head-r7 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-nofeet-r0 | AC | 100 % charged -> 100 % charged | 13.0 |
| numed-nofeet-r1 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-nofeet-r2 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-nofeet-r3 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-nofeet-r4 | AC | 100 % charged -> 100 % charged | 13.2 |
| numed-nofeet-r5 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-nofeet-r6 | AC | 100 % charged -> 100 % charged | 13.1 |
| numed-nofeet-r7 | AC | 100 % charged -> 100 % charged | 13.1 |
| tile-base-r1 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-base-r2 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-base-r3 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-base-r4 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-base-r5 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-base-r6 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-base-r7 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-basenofeet-r1 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-basenofeet-r2 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-basenofeet-r3 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-basenofeet-r4 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-basenofeet-r5 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-basenofeet-r6 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-basenofeet-r7 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-head-r1 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-head-r2 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-head-r3 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-head-r4 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-head-r5 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-head-r6 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-head-r7 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-nofeet-r1 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-nofeet-r2 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-nofeet-r3 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-nofeet-r4 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-nofeet-r5 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-nofeet-r6 | AC | 100 % charged -> 100 % charged | 0.0 |
| tile-nofeet-r7 | AC | 100 % charged -> 100 % charged | 0.0 |

### Mesh and quality (every repeat; r0 included)

| catchment | build | triangles | < 1° | share < 1° | < 10° | worst | median | max degree | max height error m | Delaunay violations | identical across runs | mesh sha256 (arrays) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Lagan | master ed125121 | 863 897 | 1598 | 0.1850 % | 46269 | 0.000412127° | 36.62° | 18 | 9.999973 | 0 of 1 129 623 (6490 near-tie  exact) | 8 runs, 1 distinct | 66056f56baf4 |
| Lagan | 20c-1 34b546c4 | 867 612 | 679 | 0.0783 % | 39304 | 0.00238403° | 36.87° | 17 | 9.999754 | 0 of 1 128 715 (6452 near-tie  exact) | 8 runs, 1 distinct | 6e397e1e78c8 |
| Lagan | 20c-1, --no-constraint-feet | 866 085 | 2047 | 0.2364 % | 47871 | 0.000391668° | 36.57° | 18 | 9.999973 | 0 of 1 133 054 (6496 near-tie  exact) | 8 runs, 1 distinct | b0adc3504b1c |
| Lagan | master, --no-constraint-feet | 866 085 | 2047 | 0.2364 % | 47871 | 0.000391668° | 36.57° | 18 | 9.999973 | 0 of 1 133 054 (6496 near-tie  exact) | 8 runs, 1 distinct | b0adc3504b1c |
| Numedalslagen | master ed125121 | 1 287 334 | 1644 | 0.1277 % | 33484 | 0.000398709° | 36.87° | 19 | 9.999972 | 0 of 1 781 776 (4737 near-tie  exact) | 8 runs, 1 distinct | 89d7629522c8 |
| Numedalslagen | 20c-1 34b546c4 | 1 290 807 | 433 | 0.0335 % | 27820 | 0.000680285° | 36.87° | 19 | 9.999972 | 0 of 1 783 018 (4692 near-tie  exact) | 8 runs, 1 distinct | 93f93b681dd2 |
| Numedalslagen | 20c-1, --no-constraint-feet | 1 303 047 | 6399 | 0.4911 % | 45443 | 7.22482e-06° | 36.87° | 19 | 9.999972 | 0 of 1 805 691 (4855 near-tie  exact) | 8 runs, 1 distinct | 08a1fb84834d |
| Numedalslagen | master, --no-constraint-feet | 1 303 047 | 6399 | 0.4911 % | 45443 | 7.22482e-06° | 36.87° | 19 | 9.999972 | 0 of 1 805 691 (4855 near-tie  exact) | 8 runs, 1 distinct | 08a1fb84834d |
| 1 m tile (tolerance 1) | master ed125121 | 464 290 | 1 | 0.0002 % | 10758 | 0.629599° | 45.00° | 74 | 0.999976 | 0 of 692 056 (94469 near-tie  exact) | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | 20c-1 34b546c4 | 464 290 | 1 | 0.0002 % | 10758 | 0.629599° | 45.00° | 74 | 0.999976 | 0 of 692 056 (94469 near-tie  exact) | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | 20c-1, --no-constraint-feet | 464 290 | 1 | 0.0002 % | 10758 | 0.629599° | 45.00° | 74 | 0.999976 | 0 of 692 056 (94469 near-tie  exact) | 7 runs, 1 distinct | 5866d649a19a |
| 1 m tile (tolerance 1) | master, --no-constraint-feet | 464 290 | 1 | 0.0002 % | 10758 | 0.629599° | 45.00° | 74 | 0.999976 | 0 of 692 056 (94469 near-tie  exact) | 7 runs, 1 distinct | 5866d649a19a |

### The 20c-1 gates (`docs/increments/20c-soft-quality.md`, PRs and gates)

| catchment | gate | threshold | measured | result |
|---|---|---|---|---|
| Lagan | share under 1° | ≤ 0.09 % | 0.0783 % | pass |
| Lagan | triangles vs master | ≤ +1 % | +0.430 % | pass |
| Lagan | worst angle | ≥ master's 0.000412127° (design: 0.000412°) | 0.00238403° | pass |
| Lagan | --no-constraint-feet bit-identical to master --no-constraint-feet | arrays and .vtk from POINTS on equal | arrays equal, .vtk equal | pass |
| Lagan | (not a gate) --no-constraint-feet vs master's default | - | differ | - |
| Lagan | determinism, 20c-1 default | one mesh over every repeat | 8 runs, 1 distinct | pass |
| Lagan | largest height error (the CLI's own check) | ≤ 10 m | 9.999754 m | pass |
| Numedalslagen | share under 1° | ≤ 0.04 % | 0.0335 % | pass |
| Numedalslagen | triangles vs master | ≤ +1 % | +0.270 % | pass |
| Numedalslagen | worst angle | ≥ master's 0.000398709° (design: 0.000399°) | 0.000680285° | pass |
| Numedalslagen | --no-constraint-feet bit-identical to master --no-constraint-feet | arrays and .vtk from POINTS on equal | arrays equal, .vtk equal | pass |
| Numedalslagen | (not a gate) --no-constraint-feet vs master's default | - | differ | - |
| Numedalslagen | determinism, 20c-1 default | one mesh over every repeat | 8 runs, 1 distinct | pass |
| Numedalslagen | largest height error (the CLI's own check) | ≤ 10 m | 9.999972 m | pass |

### Timings: median seconds over repeats r1 and on (runs per build: [7]), threads 10 (hardware concurrency)

| catchment | phase | master ed125121: median (min to max) | 20c-1 34b546c4: median (min to max) | 20c-1, --no-constraint-feet: median (min to max) | master, --no-constraint-feet: median (min to max) |
|---|---|---|---|---|---|
| Lagan | refine | 0.5510 (0.5452 to 0.5528) | 0.5790 (0.5715 to 0.5833) | 0.5445 (0.5430 to 0.5550) | 0.5439 (0.5413 to 0.5512) |
| Lagan | refine: start quality | 0.3267 (0.3218 to 0.3291) | 0.3470 (0.3425 to 0.3493) | 0.3257 (0.3248 to 0.3335) | 0.3235 (0.3224 to 0.3337) |
| Lagan | refine: scan (parallel) | 0.0514 (0.0504 to 0.0521) | 0.0504 (0.0502 to 0.0516) | 0.0506 (0.0505 to 0.0520) | 0.0508 (0.0504 to 0.0513) |
| Lagan | refine: split + flip (serial) | 0.0391 (0.0380 to 0.0413) | 0.0480 (0.0465 to 0.0500) | 0.0363 (0.0351 to 0.0364) | 0.0360 (0.0340 to 0.0365) |
| Lagan | final check: scan (parallel) | 0.3028 (0.3011 to 0.3059) | 0.3041 (0.3022 to 0.3051) | 0.3108 (0.3073 to 0.3131) | 0.3091 (0.3058 to 0.3109) |
| Lagan | final check: split + flip (serial) | 0.0304 (0.0292 to 0.0316) | 0.0339 (0.0328 to 0.0346) | 0.0319 (0.0314 to 0.0325) | 0.0316 (0.0308 to 0.0319) |
| Lagan | total | 81.2010 (80.6560 to 81.9900) | 81.5520 (80.8010 to 81.9770) | 81.5450 (81.0330 to 81.8120) | 81.5190 (80.5570 to 82.3980) |
| Numedalslagen | refine | 1.0919 (1.0828 to 1.1002) | 1.1826 (1.1685 to 1.2023) | 1.0757 (1.0654 to 1.1001) | 1.0708 (1.0636 to 1.0741) |
| Numedalslagen | refine: start quality | 0.4954 (0.4910 to 0.5050) | 0.5348 (0.5284 to 0.5471) | 0.5029 (0.4968 to 0.5243) | 0.4966 (0.4900 to 0.5009) |
| Numedalslagen | refine: scan (parallel) | 0.2854 (0.2844 to 0.2868) | 0.2854 (0.2835 to 0.2890) | 0.2863 (0.2849 to 0.2919) | 0.2853 (0.2848 to 0.2878) |
| Numedalslagen | refine: split + flip (serial) | 0.1675 (0.1652 to 0.1698) | 0.2209 (0.2162 to 0.2238) | 0.1462 (0.1435 to 0.1467) | 0.1451 (0.1427 to 0.1480) |
| Numedalslagen | final check: scan (parallel) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) |
| Numedalslagen | final check: split + flip (serial) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) |
| Numedalslagen | total | 11.3170 (11.2670 to 11.3320) | 11.3870 (11.3480 to 11.4090) | 11.2520 (11.2390 to 11.3710) | 11.2730 (11.2570 to 11.3130) |
| 1 m tile (tolerance 1) | refine | 0.1940 (0.1874 to 0.1949) | 0.2369 (0.2241 to 0.2443) | 0.1912 (0.1854 to 0.1925) | 0.1928 (0.1865 to 0.1935) |
| 1 m tile (tolerance 1) | refine: start quality | 0.0006 (0.0006 to 0.0007) | 0.0007 (0.0006 to 0.0007) | 0.0006 (0.0006 to 0.0006) | 0.0006 (0.0006 to 0.0007) |
| 1 m tile (tolerance 1) | refine: scan (parallel) | 0.0571 (0.0568 to 0.0577) | 0.0574 (0.0565 to 0.0577) | 0.0570 (0.0566 to 0.0575) | 0.0572 (0.0566 to 0.0587) |
| 1 m tile (tolerance 1) | refine: split + flip (serial) | 0.1090 (0.1031 to 0.1100) | 0.1527 (0.1415 to 0.1598) | 0.1072 (0.1020 to 0.1084) | 0.1080 (0.1019 to 0.1088) |
| 1 m tile (tolerance 1) | final check: scan (parallel) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) |
| 1 m tile (tolerance 1) | final check: split + flip (serial) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) | 0.0000 (0.0000 to 0.0000) |
| 1 m tile (tolerance 1) | total | 0.4130 (0.4060 to 0.4140) | 0.4540 (0.4410 to 0.4630) | 0.4090 (0.4040 to 0.4130) | 0.4120 (0.4050 to 0.4160) |

### Feet on against feet off, and against master (change of the medians)

| catchment | phase | 20c-1 vs 20c-1 feet off | 20c-1 vs master | 20c-1 feet off vs master feet off |
|---|---|---|---|---|
| Lagan | refine | +6.34 % | +5.07 % | +0.11 % |
| Lagan | refine: start quality | +6.53 % | +6.19 % | +0.67 % |
| Lagan | refine: scan (parallel) | -0.41 % | -2.06 % | -0.51 % |
| Lagan | refine: split + flip (serial) | +32.36 % | +22.79 % | +0.78 % |
| Lagan | final check: scan (parallel) | -2.16 % | +0.42 % | +0.55 % |
| Lagan | final check: split + flip (serial) | +6.10 % | +11.44 % | +1.02 % |
| Lagan | total | +0.01 % | +0.43 % | +0.03 % |
| Numedalslagen | refine | +9.94 % | +8.31 % | +0.46 % |
| Numedalslagen | refine: start quality | +6.34 % | +7.95 % | +1.26 % |
| Numedalslagen | refine: scan (parallel) | -0.33 % | -0.01 % | +0.35 % |
| Numedalslagen | refine: split + flip (serial) | +51.10 % | +31.87 % | +0.73 % |
| Numedalslagen | final check: scan (parallel) | - (no such phase) | - (no such phase) | - (no such phase) |
| Numedalslagen | final check: split + flip (serial) | - (no such phase) | - (no such phase) | - (no such phase) |
| Numedalslagen | total | +1.20 % | +0.62 % | -0.19 % |
| 1 m tile (tolerance 1) | refine | +23.93 % | +22.12 % | -0.87 % |
| 1 m tile (tolerance 1) | refine: start quality | +6.85 % | +4.78 % | +0.08 % |
| 1 m tile (tolerance 1) | refine: scan (parallel) | +0.64 % | +0.45 % | -0.20 % |
| 1 m tile (tolerance 1) | refine: split + flip (serial) | +42.39 % | +40.12 % | -0.72 % |
| 1 m tile (tolerance 1) | final check: scan (parallel) | - (no such phase) | - (no such phase) | - (no such phase) |
| 1 m tile (tolerance 1) | final check: split + flip (serial) | - (no such phase) | - (no such phase) | - (no such phase) |
| 1 m tile (tolerance 1) | total | +11.00 % | +9.93 % | -0.73 % |

### The 1 m benchmark and thread sweep (`tools/bench.py run`, `bench/b1m-*`; `bench_summarize.py bench b1m`)

base: 2 runs, commits ['ed12512'], power ['ac'], hardening ['libc++ fast']
head: 2 runs, commits ['34b546c'], power ['ac'], hardening ['libc++ fast']

| domain | mesh sha256 (base) | mesh sha256 (head) | equal |
|---|---|---|---|
| tile | 11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c | 11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c | True |
| quarter | a60fb597b5b738268a3bd27730989392ee86fb030fe060b78fd3a0bcdbb60a34 | a60fb597b5b738268a3bd27730989392ee86fb030fe060b78fd3a0bcdbb60a34 | True |

| domain | threads | refine_s base | refine_s head | change | app_s change | proc_s change |
|---|---:|---:|---:|---:|---:|---:|
| tile | default | 0.1938 | 0.2372 | +22.4 % | +9.1 % | +5.6 % |
| tile | 1 | 0.4818 | 0.5110 | +6.0 % | +3.5 % | +2.6 % |
| tile | 2 | 0.3131 | 0.3424 | +9.4 % | +4.6 % | +3.5 % |
| tile | 3 | 0.2614 | 0.2918 | +11.6 % | +5.3 % | +3.6 % |
| tile | 4 | 0.2362 | 0.2719 | +15.1 % | +6.7 % | +4.0 % |
| tile | 5 | 0.2181 | 0.2547 | +16.8 % | +6.8 % | +4.8 % |
| tile | 6 | 0.2070 | 0.2431 | +17.4 % | +6.8 % | +4.5 % |
| tile | 7 | 0.1969 | 0.2338 | +18.7 % | +7.8 % | +4.8 % |
| tile | 8 | 0.1916 | 0.2285 | +19.3 % | +7.6 % | +5.1 % |
| tile | 9 | 0.1946 | 0.2372 | +21.9 % | +9.0 % | +5.4 % |
| tile | 10 | 0.1936 | 0.2373 | +22.5 % | +9.2 % | +6.1 % |
| tile | 11 | 0.1939 | 0.2374 | +22.5 % | +8.9 % | +5.7 % |
| tile | 12 | 0.1930 | 0.2319 | +20.2 % | +8.4 % | +5.4 % |
| tile | 13 | 0.1938 | 0.2350 | +21.3 % | +8.7 % | +5.9 % |
| tile | 14 | 0.1946 | 0.2365 | +21.5 % | +8.5 % | +5.4 % |
| tile | 15 | 0.1942 | 0.2388 | +23.0 % | +9.7 % | +6.1 % |
| tile | 16 | 0.1951 | 0.2391 | +22.6 % | +9.1 % | +6.1 % |
| tile | 17 | 0.1946 | 0.2348 | +20.6 % | +8.5 % | +6.1 % |
| tile | 18 | 0.1945 | 0.2378 | +22.3 % | +9.1 % | +5.6 % |
| tile | 19 | 0.1940 | 0.2378 | +22.6 % | +9.1 % | +6.0 % |
| tile | 20 | 0.1939 | 0.2393 | +23.4 % | +9.3 % | +5.9 % |
| quarter | default | 0.1758 | 0.2130 | +21.1 % | +9.8 % | +6.1 % |
| quarter | 1 | 0.4250 | 0.4498 | +5.8 % | +4.0 % | +2.6 % |
| quarter | 2 | 0.2773 | 0.3036 | +9.5 % | +5.2 % | +4.0 % |
| quarter | 3 | 0.2291 | 0.2582 | +12.7 % | +6.9 % | +4.1 % |
| quarter | 4 | 0.2086 | 0.2416 | +15.8 % | +8.2 % | +5.1 % |
| quarter | 5 | 0.1941 | 0.2256 | +16.2 % | +8.2 % | +5.3 % |
| quarter | 6 | 0.1839 | 0.2153 | +17.1 % | +8.1 % | +5.0 % |
| quarter | 7 | 0.1768 | 0.2091 | +18.3 % | +8.3 % | +5.1 % |
| quarter | 8 | 0.1716 | 0.2039 | +18.8 % | +8.6 % | +5.4 % |
| quarter | 9 | 0.1759 | 0.2166 | +23.2 % | +10.8 % | +6.4 % |
| quarter | 10 | 0.1742 | 0.2142 | +23.0 % | +10.5 % | +6.0 % |
| quarter | 11 | 0.1751 | 0.2045 | +16.8 % | +8.5 % | +5.4 % |
| quarter | 12 | 0.1757 | 0.2158 | +22.8 % | +10.7 % | +6.4 % |
| quarter | 13 | 0.1757 | 0.2138 | +21.6 % | +10.3 % | +6.2 % |
| quarter | 14 | 0.1721 | 0.2083 | +21.0 % | +9.8 % | +6.3 % |
| quarter | 15 | 0.1759 | 0.2172 | +23.5 % | +10.8 % | +6.3 % |
| quarter | 16 | 0.1755 | 0.2168 | +23.6 % | +10.8 % | +6.6 % |
| quarter | 17 | 0.1752 | 0.2020 | +15.3 % | +7.9 % | +5.2 % |
| quarter | 18 | 0.1763 | 0.2178 | +23.6 % | +11.2 % | +6.8 % |
| quarter | 19 | 0.1768 | 0.2166 | +22.5 % | +10.1 % | +6.1 % |
| quarter | 20 | 0.1763 | 0.2044 | +15.9 % | +8.1 % | +5.5 % |

refine_s pooled change over every cell: +5.8 .. +23.6 %
base-vs-base noise (per cell, max/min of the base runs' medians): median 0.4 %, max 3.3 %

| domain | side | worst angle | median angle | share < 1° | max degree | within tolerance | Delaunay violations (checked) |
|---|---|---|---|---|---|---|---|
| tile | base (1 distinct of 2) | 0.6296° | 45.00° | 0.0002 % | 74 | True | 0 (692056) |
| tile | head (1 distinct of 2) | 0.6296° | 45.00° | 0.0002 % | 74 | True | 0 (692056) |
| quarter | base (1 distinct of 2) | 0.3955° | 45.00° | 0.0028 % | 18 | True | 0 (641792) |
| quarter | head (1 distinct of 2) | 0.3955° | 45.00° | 0.0028 % | 18 | True | 0 (641792) |

| domain | side | refine_s 1 thread | refine_s 20 threads | speed-up |
|---|---|---:|---:|---:|
| tile | base | 0.4818 | 0.1939 | 2.48 x |
| tile | head | 0.5110 | 0.2393 | 2.14 x |
| quarter | base | 0.4250 | 0.1763 | 2.41 x |
| quarter | head | 0.4498 | 0.2044 | 2.20 x |
<!-- end of generated tables -->
