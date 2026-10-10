# Increment 34 acceptance: slope-adaptive vertical tolerance (@perf, 2026-10-10)

Branch `worktree-slope-tol`, head `39e3ae34` (code review round 3 approved,
513 counted lines), against its base, master `9ba38490`. The rules:
`docs/increments/34-slope-tolerance.md`, section 10, "Speed judgment";
`docs/increments/README.md`, "Acceptance".

## Verdict: ACCEPTED

| gate | measured | limit | result |
|---|---|---|---|
| 1 m benchmark and thread sweep, default flags (G1) | 42 cells (2 domains x 21 thread counts): median refine time -2.0 % to +1.9 %, median of the changes +0.1 % | 5 % (bench.py's band) | `ACCEPTED` by `bench.py` |
| mesh SHA-256 unchanged at default flags (G1) | tile `11741a81…` and quarter `a60fb597…`, the same on base and head; the quick check's six default-flag meshes the same as master's baseline | equal | met |
| quality at default flags | worst angle 0.6296° / 0.3955°, max degree 74 / 18, within tolerance, 0 Delaunay violations: the same on both sides | no loss | met |
| quick check, default-flag cases against master's baseline | `NO CHANGE`; totals -0.0 % to +1.2 % | the case's own band, at least 5 % | met |
| Q7: `romsdal-slope` refine time per output triangle against uniform 2 m on tile `6901_3` | **1.039 times** (0.658 µs against 0.633 µs) | at most 1.25 times | met |
| Q7: the `slope` phase (Horn's classes over 2.55e7 nodes) | **0.1065 s** wall (0.1051-0.1066) | under 0.5 s | met |

## 1. The 1 m benchmark and the thread sweep (G1)

`tools/bench.py run`, default flags (tile and quarter domains, tolerance 1,
threads 0 and 1 to 20, 5 repeats each, interleaved), base then head back to
back: `bench-base/` (`--tree` a scratch worktree of `9ba38490`) and
`bench-head/` (`--baseline bench-base`). Both printed `ACCEPTED`. The table
of every cell is `delta.md` (`scripts/delta.py`): 42 cells, the change in
the median refine time from -2.0 % to +1.9 %, median of the changes +0.1 %,
against `bench.py`'s 5 % band. Ceiling (speed-up from 1 thread to 20, best in
brackets): tile 2.54x (2.56x at 11) base, 2.55x (2.58x at 11) head; quarter
2.53x (2.55x at 11) base, 2.58x (2.60x at 13) head.

| domain | threads | base median s | head median s | change |
|---|---:|---:|---:|---:|
| tile | 0 (default) | 0.1864 | 0.1871 | +0.3 % |
| tile | 1 | 0.4768 | 0.4798 | +0.6 % |
| tile | 20 | 0.1876 | 0.1883 | +0.4 % |
| quarter | 0 (default) | 0.1665 | 0.1668 | +0.1 % |
| quarter | 1 | 0.4232 | 0.4311 | +1.9 % |
| quarter | 20 | 0.1672 | 0.1671 | -0.1 % |

Mesh SHA-256 unchanged at default flags: tile `11741a81…`, quarter
`a60fb597…` on both sides (the same hashes as increment 33's acceptance,
`2026-10-09/vtol-field`). Quality the same on both sides: worst angle
0.6296° (tile) and 0.3955° (quarter), max vertex degree 74 and 18, within
tolerance, 0 Delaunay violations (`run.json`, `quality`).

The base `_core` (`a4d38b21…`) is byte for byte the one increment 33's
head produced: no C++ changed between `67af3081` and `9ba38490`. The head
`_core` is `a741d02f…`, the one installed in the venv. `bench-head/run.json`
records the tree as dirty: the uncommitted files were this directory and
`docs/benchmarks/quick/` only, no code (`git status --short` at the time).

## 2. The quick check: the `romsdal-slope` case and the new baseline

`docs/benchmarks/quick/cases.toml` gains `romsdal-slope`, as
`34-slope-tolerance.md` section 10 states it: tile `6901_3` of DTM10
(Romsdalen, Isfjorden, Trollveggen; 5051 by 5051 nodes), `--tolerance 10
--tolerance-slope 2 25 35`, default threads, 1 warm-up and 3 runs, and the
DTM10 directory as its input key. The pin
`test_the_shipped_cases_are_sections_3_table` in
`tests/python/test_bench_quick.py` passes with it
(`pytest tests/python/test_bench_quick.py -k shipped_cases`: 2 passed with
`-k "shipped_cases or hotspots_file"`).

`bench_quick.py run --save-baseline` at head (`RASPUTIN_DATA` set to
`../rasputin_data`; the default resolves relative to the worktree's parent
and does not exist there), judged against the stored baseline master carried
(`quick/base-67af3081.json`, increment 33's), then stored the new one,
`docs/benchmarks/quick/baseline-ac.json` (all six cases, commit `39e3ae34`).
Output: `quick/head-run.log`. Verdict `NO CHANGE`, exit 0; `romsdal-slope`
`NOT JUDGED` (not in the old baseline), as a first record is. Power
`quick/power-before.txt`, `quick/power-after.txt`.

The default-flag cases against master's baseline (`quick/delta.md`,
`scripts/quick_delta.py`; both records on AC, 2026-10-09 and today):

| case | base total s | head total s | change | refine base / head s | mesh same |
|---|---:|---:|---:|---|---|
| tile | 0.420 | 0.420 | -0.0 % | 0.187 / 0.186 | yes |
| tile, 1 thread | 0.711 | 0.717 | +0.8 % | 0.478 / 0.480 | yes |
| quarter | 0.355 | 0.357 | +0.5 % | 0.173 / 0.166 | yes |
| numedalslagen | 9.902 | 9.970 | +0.7 % | 1.209 / 1.194 | yes |
| lagan | 19.456 | 19.557 | +0.5 % | 0.663 / 0.650 | yes |
| geilo-al-ramp | 0.334 | 0.338 | +1.2 % | 0.100 / 0.100 | yes |

The new case: total 3.639 s (3.622-3.641), refine 2.032 s, `slope` 0.106 s,
largest error 9.9999 m of 10, mesh `992e49cb…`. Its phases: refine 56 %
(split and flip, serial, 33 %; scan, parallel, 16 %; setup and output 7 %),
trim 3 %, slope 3 %, decode 2 %.

### 2.1 The two point-scan paths (`PointScan` 104 to 120 bytes)

`PointScan` grew by 16 bytes on every run (as-built pin). The edge strip
(`refine_strip`) runs in every default-flag case; the final check
(`refine_points`) runs in Lagan's `--out-crs` case. Medians, min-max in
brackets, from the two quick records:

| case | phase | base s | head s | change |
|---|---|---|---|---:|
| lagan | check points: store | 0.2755 (0.2754-0.2785) | 0.2833 (0.2806-0.2854) | +7.8 ms, +2.8 % |
| lagan | final check: scan (parallel) | 0.3008 (0.3000-0.3014) | 0.3020 (0.2984-0.3046) | +1.2 ms, +0.4 % |
| lagan | final check: split + flip (serial) | 0.0331 (0.0326-0.0337) | 0.0346 (0.0340-0.0357) | +1.5 ms |
| numedalslagen | edge strip: scan (parallel) | 0.0123 (0.0120-0.0130) | 0.0137 (0.0137-0.0153) | +1.4 ms |
| numedalslagen | edge strip: split + flip (serial) | 0.0080 (0.0080-0.0085) | 0.0088 (0.0086-0.0088) | +0.8 ms |
| quarter | edge strip: scan (parallel) | 0.0009 | 0.0016 | +0.7 ms |
| tile | edge strip: scan (parallel) | 0.001 | 0.001 | +0.1 ms |

What the growth costs: milliseconds. The largest absolute change is Lagan's
`check points: store`, +7.8 ms, outside its min-max spread, which is 0.04 %
of that run's 19.56 s; the edge-strip changes are 0.1 to 1.4 ms. None of
these phases reaches the quick check's judging floor (5 % of the total), and
`refine` itself is 1 to 4 % faster on these cases, within the band. The
cause of the millisecond increases is not measured (no profile; the gates
passed): it is thought to be the larger per-slot record and the slope
branch's checks, not shown.

## 3. Q7: refine time per output triangle on tile `6901_3`

`scripts/q7.py`: one warm-up of each, then 3 runs alternating uniform 2 m
and the slope's ramp, default threads, head's Release build
(`build-bench/pkg`). "Refine" is the `--stats` row `refine`; the classes are
built before it, in the `slope` phase. Raw: `q7/raw.jsonl`; output:
`q7/result.txt`. AC power on every row.

| case | triangles | refine median s (min-max) | µs per triangle | total median s |
|---|---:|---|---:|---:|
| uniform 2 m | 5 916 344 | 3.7465 (3.7425-3.7491) | 0.633 | 7.558 |
| `--tolerance 10 --tolerance-slope 2 25 35` | 3 096 202 | 2.0365 (2.0359-2.0427) | 0.658 | 4.233 |
| **ratio** | | | **1.039** | |

The `slope` phase: 0.1065 s (0.1051-0.1066). The design's estimate for the
tile was 1.4 to 3.5 million triangles, refine about 2.4 to 2.6 s and the
whole run 4 to 5 s; measured 3.10 million, 2.04 s and 4.23 s. The uniform
2 m figure, 0.633 µs per triangle, is the probe's figure to the digit. The
run's record: 7 917 630 of 25 512 601 DEM nodes (31.0 %) held tighter
because of their slope; largest error as a share of its slope's tolerance 1
(the summary line's rounding of the stored share).

## 4. Hotspots (a `--stats` row at 40 % or more of a run)

From the quick check's output (`quick/head-run.log`) and the base record:

| case | phase | share | at base |
|---|---|---:|---|
| romsdal-slope (new) | refine | 55.8 % | no base (the flag is new); split and flip, serial, 33 %, under 40 |
| tile | refine | 44.3 % | 44.7 % |
| tile, 1 thread | refine; refine: scan (parallel) | 66.9 %; 50.1 % | 67.3 %; 50.2 % |
| quarter | refine | 46.6 % | 48.6 % |
| numedalslagen | features clip; its clean-up | 52.6 %; 51.9 % | 52.9 %; 52.3 % |
| geilo-al-ramp | decode | 43.5 % | 42.6 % |

Lagan's `features clip` is 49.6 % (49.9 % at base) against the ruled 40.7 %
in `hotspots.toml`: 8.9 points up, under the 10 that raise it again. Every
hotspot but `romsdal-slope`'s is the same phase at the same share as in
increment 33's acceptance (`2026-10-09/vtol-field/README.md`, section 6),
and none is in `hotspots.toml` with a ruling; `romsdal-slope`'s refine share
is new and goes to Ola with this run (section 10 of the increment file puts
the serial split-and-flip loop at a third of the run on a 3.1-million-triangle
mesh, which is where the refine time goes; not profiled here).

## 5. What went wrong on the way (for the record)

- The first base run was lost after its timing sweep: `bench.py run
  --mesh-dir DIR --label a/b` writes the quality mesh to `DIR/a/b_<domain>.vtk`
  and `mesh` refuses a missing parent directory, so the run exited 3 after
  all 210 timed children. Creating `DIR/slope-tol/` first fixed it. A
  `bench.py` change (make the quality path's parent, or refuse a slash in the
  label before building) needs a failing test from `@tester` first.
- `bench_quick.py`'s default data folder is `<repo>/../rasputin_data`, which
  under `.claude/worktrees/<name>` is `.claude/worktrees/rasputin_data`:
  one run exited 3 on "no such input" before `RASPUTIN_DATA` was set.

## 6. Method

- **Machine and power.** Apple M1 Max, 8 P + 2 E cores, 32 GiB, macOS 27.0,
  Python 3.14.7. AC power for every run, battery 98 to 100 %, charging, read
  before and after each run (`run.json`, the quick-check records, `q7/raw.jsonl`).
  `caffeinate` was running. Nothing else ran beside a timing run.
- **Builds.** `bench.py`'s Release builds (`-O3 -DNDEBUG`, hardening on, every
  child reporting `libc++ fast`): base from a scratch worktree of `9ba38490`,
  head from this worktree's `build-bench`. The worktree's venv `_core` was
  rebuilt from head first (`cmake --build build-pyext --target _core`, copied,
  SHA-256 `a741d02f…`), as `.claude/REQUIRED-READING.md` ("Stale artifacts")
  requires.
- **Inputs.** The 1 m benchmark: `tests/fixtures/dem_archive/7908_3_10m_z33.tif`
  (tile and quarter domains). The quick check and Q7: DTM10 at
  `../rasputin_data/DTM10_UTM33_20260925` (tile `6901_3`), the quick check's
  catchment inputs as `docs/benchmarks/quick/cases.toml` lists them.
- **Meshes** are not committed. They were in the session scratchpad
  (`meshes-base/`, `meshes-head/`, `q7/`), removed after the run; rerun
  `bench.py run` as section 1 says and `scripts/q7.py` to regenerate them.
