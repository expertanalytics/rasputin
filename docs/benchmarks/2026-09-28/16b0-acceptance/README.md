# Increment 16b-0 acceptance (@perf, 2026-09-28)

**Verdict: ACCEPTED** (battery against a battery baseline measured back to
back). The standard 1 m benchmark is unchanged and the meshes are
bit-identical. `node` on M3's 48 km CORINE square, layout B, fell from
84.1 s to 0.144 s. Short summary: `../16b0-acceptance.md`.

## What was measured, and on what

- **Branch** `increment16b0-verifier-sweep` at `4e83389` (master merged in).
  **Base**: master `dc372a5` (#106), in a temporary detached `git worktree`
  in the session scratchpad. The worktree was removed afterwards.
  `git diff --stat dc372a5 4e83389 -- include src_cpp src_python bindings`
  lists one file, `include/terrain/noding/noded_pslg_builder.hpp`. That header
  is included by `node.hpp`, `core/noded_pslg.hpp` and `bindings/core.cpp`.
  Nothing under `include/terrain/refinement/` includes it.
- **Builds**: `bench.py`'s own Release builds into `<tree>/build-bench`.
  `_core` sha256 prefixes: base `2b69e078db7f`, 16b-0 `190587e7ba0e`, the
  same in all three runs of each. Parts 2-4 import those same `build-bench/pkg`
  trees (bench.py's child technique: drop the editable finder, put the pkg
  first, assert that `tin_engine` and `_core` load from it). So the old
  verifier was measured with the same compiler and flags as the new one. No
  C++ was patched, so there was no scratch build to sanitize.
- **Machine**: Apple M1 Max (8 P + 2 E), 32 GiB, macOS 27.0.
  `caffeinate -ims` was held for every run.
- **Power: battery, discharging, all runs.** The percentages fell from 76 %
  (15:21) to 62 % (15:53). The `pmset -g batt` line is logged before every
  case in `corine/logs/*.log`, and before and after each bench run in its
  `run.json`. An earlier attempt at 12:22 hit a Mac sleep and a bench.py
  mesh-dir error. It produced no evidence, and none of its timings are used
  (the log was kept in the scratchpad as `pairs-void-sleep.log`, not here).
- **Mesh digests.** The bench meshes are hashed by bench.py. For `node`, a
  digest over the noded vertices and the per-edge masks was taken (`drv.py`,
  `digest`).

## 1. The 1 m benchmark and thread sweep (`tools/bench.py`, blob `77765b1`)

Three pairs, run in the order base, 16b-0, base, 16b-0, base, 16b-0, between
15:21 and 15:53 CEST. Every 16b-0 run names its base explicitly
(`--baseline`). Command (`corine/scripts/pairs.sh`, `pair3.sh`):

```
.venv/bin/python tools/bench.py run --label 16b0-acceptance/<label> [--tree <worktree of dc372a5>] \
    [--baseline docs/benchmarks/2026-09-28/16b0-acceptance/base-dc372a5[-rN]] --mesh-dir <scratch> \
    --domain docs/benchmarks/2026-09-26/quarter.geojson --domain tile \
    --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --repeats 5 --threads 1,...,20
```

bench.py's verdicts: pair 1 **ACCEPTED**; pair 2 **REGRESSION** on two cells
(`tile refine_s[t=10] +6.9 %`, `tile refine_s[t=16] +8.0 %`); pair 3
**ACCEPTED**. To reproduce:

```
.venv/bin/python tools/bench.py compare docs/benchmarks/2026-09-28/16b0-acceptance/16b0-r2 \
    --baseline docs/benchmarks/2026-09-28/16b0-acceptance/base-dc372a5-r2
```

**Pair 2's two cells are judged noise, not a regression.** The reasons:

- The same cells were +1.9 % / +0.5 % in pair 1 and +2.1 % / +1.2 % in pair 3.
- The same build moves by more than that between runs. Base run 2 against
  base run 1 exceeds 5 % in 3 of 42 cells, by up to +13.5 %. 16b-0 run 2
  against 16b-0 run 1 exceeds 5 % in 3 of 42 cells, by up to +10.0 %.
- `refine_s` times the `refine` call alone. The only changed header is not
  compiled into refine's code path (see above).
- Pooled over all 42 cells, each cell's median over the three base runs against its median over the three 16b-0 runs gives a median change of -0.28 %
  (min -5.9 %, max +2.8 %), and no cell is above +5 %.

Refine seconds, median of 5 per run (ms). The last column gives the median of
the three pair changes, and the three changes themselves:

| domain | threads | base r1 | 16b-0 r1 | base r2 | 16b-0 r2 | base r3 | 16b-0 r3 | change |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| tile | 1 | 471.7 | 480.9 | 535.6 | 477.1 | 466.7 | 464.1 | -0.6 % (+2.0, -10.9, -0.6) |
| tile | 2 | 312.6 | 313.6 | 346.7 | 309.2 | 303.6 | 302.7 | -0.3 % (+0.3, -10.8, -0.3) |
| tile | 4 | 228.0 | 234.2 | 237.8 | 229.6 | 227.7 | 228.8 | +0.5 % (+2.7, -3.5, +0.5) |
| tile | 8 | 187.4 | 184.4 | 192.0 | 187.2 | 184.9 | 184.6 | -1.6 % (-1.6, -2.5, -0.2) |
| tile | 10 | 184.1 | 187.7 | 186.0 | 198.9 | 180.9 | 184.7 | +2.1 % (+1.9, +6.9, +2.1) |
| tile | 16 | 186.3 | 187.2 | 190.6 | 205.9 | 185.9 | 188.1 | +1.2 % (+0.5, +8.0, +1.2) |
| tile | 20 | 189.2 | 195.8 | 188.0 | 189.2 | 186.8 | 186.8 | +0.6 % (+3.5, +0.6, +0.0) |
| tile | default (10) | 185.9 | 192.2 | 190.1 | 192.3 | 187.6 | 186.0 | +1.2 % (+3.4, +1.2, -0.9) |
| quarter | 1 | 477.2 | 453.9 | 433.1 | 424.1 | 450.5 | 407.7 | -4.9 % (-4.9, -2.1, -9.5) |
| quarter | 2 | 306.3 | 277.5 | 275.3 | 276.5 | 271.4 | 263.7 | -2.8 % (-9.4, +0.4, -2.8) |
| quarter | 4 | 207.8 | 201.9 | 211.4 | 201.3 | 202.6 | 198.7 | -2.9 % (-2.9, -4.8, -1.9) |
| quarter | 8 | 173.0 | 170.8 | 165.6 | 163.4 | 163.4 | 161.6 | -1.3 % (-1.3, -1.3, -1.1) |
| quarter | 10 | 172.4 | 165.2 | 164.6 | 164.4 | 163.4 | 163.6 | -0.1 % (-4.2, -0.1, +0.2) |
| quarter | 16 | 169.9 | 165.5 | 165.8 | 170.3 | 167.0 | 161.7 | -2.6 % (-2.6, +2.7, -3.1) |
| quarter | 20 | 174.3 | 165.7 | 166.3 | 169.1 | 167.1 | 165.5 | -0.9 % (-4.9, +1.7, -0.9) |
| quarter | default (10) | 180.5 | 167.1 | 168.8 | 164.5 | 164.4 | 164.8 | -2.6 % (-7.4, -2.6, +0.3) |

**Ceiling** (speed-up over 1 thread: the best, and at 20 threads), against
the 2.0-2.2x battery baseline:

| run | tile | quarter |
|---|---|---|
| base r1 / r2 / r3 | 2.56x @10 (2.49x) / 2.88x @10 (2.85x) / 2.58x @10 (2.50x) | 2.81x @16 (2.74x) / 2.64x @14 (2.60x) / 2.76x @15 (2.70x) |
| 16b-0 r1 / r2 / r3 | 2.61x @8 (2.46x) / 2.55x @8 (2.52x) / 2.55x @11 (2.48x) | 2.77x @15 (2.74x) / 2.59x @8 (2.51x) / 2.52x @8 (2.46x) |

The ceiling follows the noisy 1-thread cell. Base r2's 2.88x comes from its
slow 1-thread tile sample (535.6 ms). From 8 threads on, the curve is flat in
every run.

**Quality: identical** in all six runs.

| domain | worst angle | max degree | within tolerance | Delaunay violations | mesh sha256 (ASCII, POINTS on) |
|---|---:|---:|---|---:|---|
| quarter | 0.3955 deg | 18 | yes | 0 of 641,791 | `ccebf96a86c6c5e244e4a0281919de4e866fcfe789b66024a290ac2af33771a1` |
| tile | 0.6296 deg | 74 | yes | 0 of 692,056 | `11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c` |

**Found in passing, not 16b-0's.** The three base runs had no `--baseline`,
so bench.py judged master `dc372a5` against the newest comparable stored
ancestor. That is `2026-09-27/21b-r2` (`a03ae65`, battery; identified by its
stats). It reported REGRESSION in quarter at 1-3, 7 and 11 threads (up to
+12.3 %) in r1, in tile at 1, 2 and 13 threads in r2, and in quarter at 1
thread (+6.1 %) in r3. The cells differ from run to run, the 1-thread cell is
the noisiest on battery, and 21b's own evidence has an unchanged build moving
by about 5 %. So no master slowdown is claimed from this. It has not been
re-measured.

The `(dirty)` flag on the 16b-0 runs means only that the untracked evidence
directories existed: `git status --porcelain --untracked-files=no` was empty.

## 2. `node` on the CORINE squares (M3, R8 step 0)

**Inputs.** Squares centred at (325 000, 6 675 000) in `6603_4`, the same
squares as R8's step 0 (the 48 km square is 301 000-349 000 ×
6 651 000-6 699 000). The linework is CORINE 2018 v2020_20u1 from Ola's
European GeoPackage (`../rasputin_data/corine_sql/...gpkg`, 8.86 GB), read by
`clc_probe.py` without GDAL. The DEM's archive does not enter `node`; the
tile is used only for its metadata. It was
`DTM10_UTM33_20260925/6603_4_10m_z33.tif` (sha256 `65b2105a896133e2…`).

- **Layout A**: the linework deduplicated, clipped as lines and line-merged
  (`clc_mesh.py`).
- **Layout B**: every clipped ring as a closed `Breakline`, so every shared
  edge is given twice (`clc_noder.py`'s letter, and the layout 16b uses).

`node` is the `start mesh: node` phase of `cli._engine`, 1 mm snap. 16b-0 is
the median of 3; master is 1 run (`corine/scripts/drv.py`,
`corine_node.sh`). The noded digests of A and B are equal at every size, and
equal between master and 16b-0.

| square | layout | input segments | noded vertices / edges | `node`, master `dc372a5` | `node`, 16b-0 | speed-up | step 0 / brief reference |
|---|---|---:|---|---:|---:|---:|---|
| 12 km | A | 3 767 | 3 728 / 3 763 | 0.134 s | 0.0041 s | 33x | 0.134 / 0.0041 s |
| 24 km | A | 14 843 | 14 698 / 14 839 | 2.06 s | 0.0174 s | 119x | 2.07 / 0.0167 s |
| 36 km | A | 30 514 | 30 238 / 30 510 | 8.69 s | 0.0364 s | 239x | 8.99 / 0.0364 s |
| 48 km | A | 51 792 | 51 395 / 51 788 | 24.85 s | 0.0668 s | 372x | 24.83 / 0.0670 s |
| 12 km | B | 7 486 | 3 728 / 3 763 | 0.449 s | 0.0078 s | 58x | — |
| 24 km | B | 29 586 | 14 698 / 14 839 | 6.88 s | 0.0352 s | 196x | — |
| 36 km | B | 60 864 | 30 238 / 30 510 | not run | 0.0764 s | | — |
| **48 km** | **B** | 103 419 | 51 395 / 51 788 | **84.05 s** | **0.1435 s** | **586x** | M3: 83.1 s |
| 96 km | A | 201 694 | 200 028 / 201 690 | not run (about 6 min extrapolated at n²) | 0.336 s | | — |
| 96 km | B | 403 005 | 200 028 / 201 690 | not run | 0.803 s | | — |
| 144 km | A | 435 049 | 430 886 / 435 045 | not run | 0.877 s | | — |
| 144 km | B | 869 427 | 430 886 / 435 045 | not run | 2.18 s | | — |

**The quadratic is gone.** From 12 to 48 km, the segments grow 13.8x. With
16b-0, `node` grows 16.3x (A) and 18.4x (B); on master it grows 185x and 187x.
From 48 to 144 km, the segments grow 8.4x and `node` grows 13.1x (A) and 15.2x
(B). That is an exponent of about 1.2-1.3 in the segment count: above
n log n, far below 2. Where the excess goes has not been profiled. The 144 km
block (20 736 km², about half of Glomma's catchment by area) nodes in 2.2 s
in layout B. Peak RSS of the process, including Python, shapely and the
GeoPackage read: 1.17 GB (B) and 1.47 GB (A) at 144 km.

## 3. The sweep's admitted worst case

Two synthetic PSLGs (`drv.py ladder|comb`), an `Outer` square with 2-point
`Breakline`s:

- **ladder**: n east-west and n north-south lines crossing, n² crossings.
  This is the GIL test's fixture at n = 300.
- **comb**: n long, parallel east-west lines with no crossing. Every
  x-interval overlaps every other; this is R8's "worst case is still
  quadratic".

The domain side is 700 m, scaled by n/300 above 300 so that the line spacing
stays at about 2.3 m.

| case | n | noded edges | `node`, 16b-0 (median of 3) | step | `node`, master |
|---|---:|---:|---:|---|---:|
| ladder | 300 | 180 604 | 0.262 s | | 210.2 s |
| ladder | 600 | 721 204 | 1.98 s | 7.6x for 4x edges | stopped after ~1 min (at 210 s × 16 it would take ~56 min) |
| ladder | 1 200 | 2 882 404 | 17.3 s | 8.7x for 4x edges | not run |
| comb | 1 000 | 1 004 | 0.0064 s | | 0.0166 s |
| comb | 4 000 | 4 004 | 0.092 s | 14.4x for 4x | 0.246 s |
| comb | 16 000 | 16 004 | 1.37 s | 14.9x for 4x | not run |
| comb | 64 000 | 64 004 | 21.4 s | 15.6x for 4x | not run |

- **The comb is quadratic as admitted**: 4x the lines costs 14.4-15.6x. At
  64 000 parallel lines spanning the domain, `node` takes 21 s. On the ladder
  the sweep is about 800x faster than master at n = 300. The ladder grows as
  about edges^1.5 (n³ in rungs): the n² crossings are output, and at any x
  about n horizontal pieces are active.
- The split between the driver and the verifier has **not been isolated** on
  these cases (no build with the verifier disabled was run). That the
  remaining cost is the sweep's active set is thought to be the case, not
  measured.
- CORINE is far from this case. The 144 km block (870 k segments) nodes in
  2.2 s, against 21 s for 64 k comb lines.

## 4. The first CORINE baseline: M4's 48 km square

The square is 301 000-349 000 × 6 651 000-6 699 000 in `6603_4`
(`DTM10_UTM33_20260925`), layout A. It goes through `cli._dem_mesh` with the
square as the domain, constraint feet on, 1 mm snap, and the default thread
count, using 16b-0's `build-bench/pkg`. There were 3 runs per row, each in
its own process (`corine_mesh.sh` run three times, 15:46-15:48, battery
64 %). Triangles and quality were identical in all 3; the times are the
median. `whole` is the `_dem_mesh` call. It excludes the GeoPackage read and
clip and the tile load (not timed). Peak RSS is `/usr/bin/time -l`
of the whole process.

| tol | features | start quality | triangles | vs none | start quality nodes | feet | min angle median | < 1 deg | worst | `node` | start quality | refine (of which scan / split) | whole | peak RSS |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | none | 25 deg | 6 581 348 | 1.00x | 0 | 0 | 36.9 | 0.00 % | 0.217 | — | — | 5.06 s (1.50 / 3.22) | 5.34 s | 1.93 GB |
| 1 m | none | off | 6 581 348 | 1.00x | 0 | 0 | 36.9 | 0.00 % | 0.217 | — | — | 4.97 s (1.49 / 3.15) | 5.25 s | 1.96 GB |
| 1 m | CORINE | 25 deg | 6 752 517 | 1.03x | 154 640 | 30 311 | 36.9 | 0.03 % | 0.0008 | 0.068 s | 0.196 s | 5.48 s (0.94 / 3.98) | 5.99 s | 2.35 GB |
| 1 m | CORINE | off | 6 655 862 | 1.01x | 0 | 31 990 | 36.9 | 0.04 % | 0.031 | 0.067 s | — | 5.29 s (1.02 / 3.90) | 5.78 s | 2.27 GB |
| 10 m | none | 25 deg | 133 379 | 1.00x | 0 | 0 | 34.5 | 0.03 % | 0.258 | — | — | 0.60 s (0.56 / 0.03) | 0.60 s | 0.41 GB |
| 10 m | none | off | 133 379 | 1.00x | 0 | 0 | 34.5 | 0.03 % | 0.258 | — | — | 0.61 s (0.57 / 0.03) | 0.62 s | 0.41 GB |
| 10 m | CORINE | 25 deg | **480 961** | **3.61x** | 154 640 | 350 | 36.9 | 0.25 % | 0.0013 | 0.067 s | 0.201 s | 0.38 s (0.09 / 0.03) | 0.60 s | 0.69 GB |
| 10 m | CORINE | off | 212 255 | 1.59x | 0 | 359 | 28.6 | 3.45 % | 0.015 | 0.066 s | — | 0.22 s (0.13 / 0.04) | 0.44 s | 0.62 GB |

"min angle median" is the median over triangles of each triangle's smallest
angle, in degrees; "worst" is the smallest in the mesh, in degrees.

- Every run met its tolerance (max error 1.0 m and 9.9995-9.99996 m), with 0
  valid DEM nodes uncovered and 0 feet refused.
- **The quality start triples the 10 m mesh: 3.61x the featureless
  mesh, and 2.27x the same features with the start quality off.** It inserts
  154 640 start-quality nodes at both tolerances. This is 20c's input. At 1 m
  the same nodes add only 1.5 % over quality off.
- **`node` is no longer the cost.** It was 24.7 s of M4's 25-30 s. It is now
  0.07 s, less than the start-quality pass (0.20 s). At 10 m with CORINE, the
  whole `_dem_mesh` call takes 0.44-0.60 s.
- **Against M4** (design time, which DEM archive is not recorded, before 21c):
  the triangle counts differ by 0.2-2.1 % (for example 1 m CORINE 6 612 554 in
  M4 against 6 752 517 here; 10 m quality 480 127 against 480 961), and the
  start-quality nodes are equal (154 640). This baseline is the first one
  measured by @perf with the method written down. M4 is context, not a
  baseline.

## Not measured

- A pure-verifier split (a build with `check_guarantee_14` disabled) on the
  ladder and the comb. The breakdown of what remains is thought to be, not
  measured.
- Master on the 36 km layout B, on 96 km and 144 km, and on ladder 600-1200
  and comb 16 000-64 000. These would take tens of minutes to hours at the
  measured n².
- An AC run of any of the above. There is no AC 16b-0 evidence.

## Files

- `base-dc372a5[-r2|-r3]/`, `16b0[-r2|-r3]/`: bench.py run directories
  (`run.json`, `raw.tsv`, generated `README.md`).
- `corine/scripts/`: `drv.py` (the one-case driver: node, ladder, comb, mesh),
  `corine_node.sh`, `corine_mesh.sh`, `pairs.sh`, `pair3.sh`, and the scratch
  `clc_probe.py` and `clc_mesh.py` that `drv.py` imports (copied from the
  session scratchpad; they read the GeoPackage relative to the repo root).
- `corine/logs/`: raw output with a `pmset` line before every case
  (`node.log`, `mesh-{1,2,3}.log`, `pairs.log`, `pair3.log`).
- The meshes are not kept (they stayed in the session scratchpad). To
  regenerate them, rerun bench.py as above, or `drv.py mesh`. Run from the
  repo root: `.venv/bin/python docs/benchmarks/2026-09-28/16b0-acceptance/corine/scripts/drv.py
  --pkg <tree>/build-bench/pkg mesh --side 48 --tol 10 --features A --min-angle 25`
  (the committed `drv.py` also looks for `clc_*.py` beside itself).
