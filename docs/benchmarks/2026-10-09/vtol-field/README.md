# Increment 33 acceptance: feature-adaptive vertical tolerance (@perf, 2026-10-09)

Branch `worktree-vtol-field`, head `67af3081` (code review round 3
approved), against its base, master `4cd7e050`. The rules:
`docs/increments/33-feature-tolerance.md`, section 10, "Speed judgment";
`docs/increments/README.md`, "Acceptance".

## Verdict: ACCEPTED

| gate | measured | limit | result |
|---|---|---|---|
| 1 m benchmark and thread sweep, default flags (G1) | 42 cells (2 domains x 21 thread counts): median refine time -2.4 % to +2.4 %, median of the changes -0.2 % | 5 % (bench.py's band) | `ACCEPTED` by `bench.py` |
| mesh SHA-256 unchanged at default flags (G1) | tile `11741a81…` and quarter `a60fb597…`, the same on base and head | equal | met |
| quality at default flags | worst angle 0.6296° / 0.3955°, max degree 74 / 18, within tolerance, 0 Delaunay violations: the same on both sides | no loss | met |
| Geilo-Ål ramp, refine time per output triangle against uniform 1 m on the same section (Q8) | **3.44 times** (1.572 µs against 0.457 µs) | at most 5 times | met |

## 1. The 1 m benchmark and the thread sweep

`tools/bench.py run`, default flags (tile and quarter domains, tolerance 1,
threads 0 and 1 to 20, 5 repeats each, interleaved), base then head back to
back: `bench-base/` (`--tree` a scratch worktree of `4cd7e050`) and
`bench-head/` (`--baseline bench-base`). Both printed `ACCEPTED`. The table
of every cell is `delta.md` (`scripts/delta.py`); the extremes are +2.4 %
(quarter, 13 threads) and -2.4 %. Ceiling (speed-up from 1 thread to 20,
best in brackets): tile 2.47x (2.56x at 12) base, 2.50x (2.54x at 18) head;
quarter 2.52x (2.54x at 14) base, 2.53x (2.55x at 14) head.

| domain | threads | base median s | head median s | change |
|---|---:|---:|---:|---:|
| tile | 0 (default) | 0.1919 | 0.1885 | -1.7 % |
| tile | 1 | 0.4788 | 0.4768 | -0.4 % |
| tile | 20 | 0.1935 | 0.1906 | -1.5 % |
| quarter | 0 (default) | 0.1695 | 0.1683 | -0.7 % |
| quarter | 1 | 0.4232 | 0.4250 | +0.4 % |
| quarter | 20 | 0.1679 | 0.1680 | +0.1 % |

`bench-head/run.json` records the tree as dirty: the uncommitted files were
this directory and `docs/benchmarks/quick/` only, no code
(`git status --short` at the time).

## 2. The quick check: the Geilo-Ål ramp case and the first baseline

`docs/benchmarks/quick/cases.toml` gains `geilo-al-ramp`: DTM10, the
Geilo-Ål section domain, `--tolerance 20 --tolerance-near bergensbanen.geojson 1
--tolerance-ramp 0 3000`, default threads, 1 warm-up and 3 runs. The inputs
are made by `bergen-line/prepare.py` (its README gives the source, the
licence and the selection).

No quick-check baseline existed (`docs/benchmarks/quick/` had none), so
`bench_quick.py run --save-baseline` at head stored the first one,
`docs/benchmarks/quick/baseline-ac.json` (all five cases; it printed
`NO BASELINE: no stored baseline`, exit 2, as a first save does). The
geilo-al-ramp case: 0.334 s in total (median), refine 0.100 s, 64 162
triangles, largest error 19.984 m.

The default-flag cases were also run at the base, with `cases.toml` as on
master, to compare (`quick/base-4cd7e050.json`, `quick/delta.md`):

| case | base total s | head total s | change | mesh same |
|---|---:|---:|---:|---|
| tile | 0.429 | 0.420 | -2.1 % | yes |
| tile, 1 thread | 0.716 | 0.711 | -0.8 % | yes |
| quarter | 0.352 | 0.355 | +1.0 % | yes |
| numedalslagen | 9.890 | 9.902 | +0.1 % | yes |
| lagan | 19.610 | 19.456 | -0.8 % | yes |

## 3. Q8: refine time per output triangle on the Geilo-Ål section

`scripts/q8.py`: one warm-up of each, then 5 runs alternating uniform 1 m and
the ramp, default threads, head's Release build (`build-bench/pkg`).
"Refine" is the `--stats` row `refine`. Raw: `q8/raw.jsonl`; output:
`q8/result.txt`.

| case | triangles | refine median s (min-max) | µs per triangle |
|---|---:|---|---:|
| uniform 1 m | 970 037 | 0.4436 (0.4399-0.4602) | 0.457 |
| ramp, 1 m to 20 m over 0-3 000 m | 64 162 | 0.1008 (0.1005-0.1016) | 1.572 |
| **ratio** | | | **3.44** |

The design estimated 2.2 to 4.6 times; the probe estimated about 53 000
triangles for the ramp, and the run gives 64 162 (+21 %).

## 4. `detail::max_error_near`, timed

A timing patch (`max-error-near/timing.patch`: a steady clock around the
function body, a per-block count of `policy.at` queries, one stderr line per
call) in a scratch copy of `67af3081`. The scratch build ran once first with
`-fsanitize=undefined -fno-sanitize-recover=all` and libc++'s extensive
hardening on the Geilo-Ål ramp case, through `_core`: clean, the same 64 162
triangles. Then Release (`-O3 -DNDEBUG`, hardening fast). Output:
`max-error-near/runs.txt`. Each run makes two calls.

| case | call | seconds | slots | queries |
|---|---|---:|---:|---:|
| Geilo-Ål ramp (3 runs) | first | 0.000498-0.000555 | 64 156 | 964 |
| | second | 0.000141-0.000164 | 64 162 | 1 |
| Hokksund-Bergen corridor (1 run) | first | 0.0220 | 1 587 214 | 30 588 |
| | second | 0.0015 | 1 587 426 | 14 |

On the corridor it is 0.024 s of a 6.3 s run (0.4 %); code review round 1's
estimate for the serial loop it replaced was 4 to 11 s.

## 5. The corridor, once, for the record (not judged)

`scripts/corridor.sh`: Hokksund to Bergen, the three lines, the ramp, head's
Release build, default threads, one run, `/usr/bin/time -l`. Output:
`corridor/` (`--stats` file, stderr, power before and after).

| measure | value |
|---|---|
| output triangles | 1 587 426 (design estimate 1.37 million) |
| total (`--stats`) | 6.369 s; wall 6.94 s (design estimate 4 to 7 s) |
| decode | 1.361 s, 21.4 % |
| refine | 4.367 s, 68.6 % (scan, parallel, 3.723 s, 58.5 %; split and flip, serial, 0.550 s, 8.6 %) |
| refine per output triangle | 2.75 µs |
| largest error / near the lines | 19.9987 m / 0.99996 m |
| peak memory footprint | 2.91 GB |

## 6. Hotspots (a `--stats` row at 40 % or more of a run)

| run | phase | share | at base too |
|---|---|---:|---|
| corridor ramp | refine | 68.6 % | no base (the flag is new) |
| corridor ramp | refine: scan (parallel) | 58.5 % | no base |
| geilo-al-ramp (quick) | decode | 43 % | no base |
| tile (quick) | refine | 45 % | yes, 45 % |
| tile, 1 thread (quick) | refine; refine: scan (parallel) | 67 %; 50 % | yes, 67 %; 50 % |
| quarter (quick) | refine | 49 % | yes, 47 % |
| numedalslagen (quick) | features clip; its clean-up | 53 %; 52 % | yes, the same |

Lagan's `features clip` is 9.71 s of 19.46 s, 49.9 %, against the ruled
40.7 % in `hotspots.toml`: 9.2 points up, under the 10 points that raise it
again; it was 49.7 % at base. The corridor's scan share is the field's
per-triangle queries by the design's account (section 10); no profile was
taken, so that is not measured.

## Method

- **Machine and power.** Apple M1 Max, 8 P + 2 E cores, 32 GiB, macOS 27.0,
  Python 3.14.7, numpy 2.5.3, shapely 2.2.0 (GEOS 3.14.1), pyproj 3.8.0
  (PROJ 9.8.1). AC power for every run, battery 100 %, charged; read before
  and after each run (`run.json`, the quick-check records, `corridor/`,
  `max-error-near/power.txt`). `caffeinate` was running. Nothing else ran.
- **Builds.** `bench.py`'s Release builds (`-O3 -DNDEBUG`, AppleClang 21,
  hardening on, `libc++ fast` in every child): base `_core` SHA-256
  `28f1db86…`, head `a4d38b21…`. The quick check, Q8 and the corridor run
  from the head's `build-bench/pkg`.
- **Inputs.** DTM10 at `../rasputin_data/DTM10_UTM33_20260925`; Bane NOR's
  Banenettverk (NLOD), prepared as `bergen-line/README.md` says.
- **Meshes** are not committed. They were in the session scratchpad
  (`meshes-base/`, `meshes-head/`, `q8/`, `corridor_ramp.vtk`), removed after the run;
  rerun the scripts named above to regenerate them.
