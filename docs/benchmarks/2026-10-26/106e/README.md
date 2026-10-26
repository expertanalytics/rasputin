# 106e acceptance: `mesh --simplify-domain` (@perf, 2026-10-26)

Head `worktree-demo-fixes-e` at `3736e39a` (106e's green step `2a351e1d`,
its tests `92e39f54`, the kill record `3736e39a`). Base: master `e0d129c0`,
without 106c (without its flags 106c changes no mesh), from a detached
scratch worktree, run back to back with the
head. The rule is the `@perf` paragraph of section 8 of
`docs/increments/106-wednesday-cells.md` (106e).

Apple M1 Max (8 P + 2 E cores), 32 GiB, macOS 27.0.1, Python 3.14.7,
shapely 2.2.0 (GEOS 3.14.1). Every run on AC (`pmset -g batt`: "Now
drawing from 'AC Power'", 100 %), read before and after each run
(`power.txt`). 106e has no C++ (`git diff e0d129c0 3736e39a -- include src
bindings CMakeLists.txt` is empty): the head's and master's Release `_core`
are the same bytes (sha256 `a353ac1e…`), and the same bytes are installed in
the worktree's venv that ran the outline runs. No sanitizer run: no C++
patch. Desktop applications were running; nothing else of mine ran beside a
timing run.

## Verdict

**ACCEPTED.** The 1 m benchmark's meshes are byte-identical to master's, in
the same time; Ljungan's run is within the design's time and its triangle
budget; Kultsjön traced gives the predicted 502. Two of the design's figures
differ, both explained below: Ljungan keeps 198 outline vertices, not 195,
and meshes to 483 triangles, not 503.

## 1. The 1 m benchmark (`tools/bench.py`, `--hardening on`, no flag)

| run | commit | tile hash | quarter hash | tile t0 | tile t1 | tile t12 | tile t20 | tile ceiling | quarter t0 | quarter t1 | quarter t12 | quarter t20 | quarter ceiling |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `106e-head-ac-r1` | 3736e39a | 11741a81 | a60fb597 | 0.2874 | 0.7069 | 0.2914 | 0.3061 | 2.31x | 0.2108 | 0.5784 | 0.2151 | 0.2143 | 2.70x |
| `106e-base-ac-r1` | e0d129c0 | 11741a81 | a60fb597 | 0.2996 | 0.6972 | 0.2940 | 0.2886 | 2.42x | 0.2125 | 0.5699 | 0.2128 | 0.2169 | 2.63x |

Refine medians in seconds, 5 repeats; ceiling is t1 over t20. Quality, both
runs: worst angle 0.6296° / 0.3955°, max degree 74 / 18, within tolerance,
0 Delaunay violations of 692,056 / 641,792 edges. `bench.py compare
106e-head-ac-r1 --baseline 106e-base-ac-r1`: 2 of 42 counts over 5 % (tile
t7 +7.6 %, t20 +6.1 %), from the same `_core` bytes as the base, so noise of
the day. Against the stored `106c-head-ac-r2` (which `bench.py run` picked
itself), both today's head and today's master are 5 to 15 % slower on most
counts, as on 106c's day; thought to be the desktop load, not profiled.

## 2. Outlines with `--simplify-domain --triangles 500`

Command: `rasputin mesh --dem glo30 --out-crs EPSG:3006 --simplify-domain
--triangles 500 --domain D [--land-cover clc2018] --out OUT --stats
NAME.stats.md`, one run each. Ljungan's `D` is
`~/rasputin-output/opal-se2-2026-10-21/catchments-original/ljungan.geojson`
(EPSG:3006, 111,248 vertices); Kultsjön traced is
`~/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon-lake-seeded-check.geojson`
(EPSG:4326, 617 vertices), the outline the design's probes call "Kultsjön,
terrain". The runs with `--land-cover clc2018` wrote `.gpkg`, the others
`.vtk`.

| outline | land cover | band | vertices after (design) | triangles (design) | worst angle | max degree | tolerance found | largest error |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Ljungan, official | clc2018 | 962.57 m | **198** (195) | **483** (503) | 0.214° | 10 | 246.8 m | 243.9 m |
| Ljungan, official | none | 962.57 m | 198 | 483 | 0.214° | 10 | 246.8 m | 243.9 m |
| Kultsjön, traced | clc2018 | 281 m | **122** (122) | **502** (502) | 0.329° | 9 | 152.3 m | 152.2 m |
| Kultsjön, traced | none | 281 m | 122 | 502 | 0.329° | 9 | 152.3 m | 152.2 m |

Every mesh within its tolerance and its budget's 5 % (at most 525). The
land cover does not change the mesh: same counts, angles and errors with
and without it. Ljungan's land cover: 2,547 class shares over 483
triangles; Kultsjön's: 1,361 over 502.

**198, not 195.** The design's probe (`outline_budget.py`) simplified the
ring as the file stores it, clockwise; `read_domain` orients it
counter-clockwise first, and the same band then keeps 198
(`orient_probe.py`: clockwise 195, counter-clockwise 198, the CLI's 198).
**483, not 503.** The design's 503 is a mesh of the probe's 195-vertex
outline; the CLI's 198-vertex outline gives 483, with or without land
cover, inside the budget's 5 %.

**Times** (seconds, the `--stats` phase table):

| | Ljungan, clc2018 | Ljungan, none | Kultsjön traced, clc2018 | Kultsjön traced, none |
|---|---:|---:|---:|---:|
| domain read | **0.100** | 0.095 | 0.002 | 0.002 |
| domain simplify | **2.573** | 2.534 | 0.178 | 0.170 |
| everything else but land cover | 4.601 | 4.219 | 0.795 | 0.800 |
| total less "land cover fractions" | 7.274 | **6.848** | 0.975 | 0.972 |
| land cover fractions (not in the prediction) | 51.299 | | 46.631 | |
| run total | 58.573 | 6.848 | 47.606 | 0.972 |
| wall clock (`time`) | 58.94 | 7.23 | 47.95 | 1.41 |

Prediction for Ljungan, the probe's mesh time measured without land cover:
(4.0 mesh + 2.4 simplify + 0.10 read) x 1.1 = 7.15 s. The like run, without
land cover, takes **6.85 s**, under it. With land cover the rest takes
7.27 s, 1.7 % over, mostly in "decode" (1.92 s against 1.63 s) and
"other" (0.22 s against 0.09 s); one run each, so within the noise of one
sample.

The figure the design left unmeasured, set from these runs: reading
Ljungan's outline, **0.10 s**. The simplification takes **2.57 s** (the
probe 2.43 s; the reviewer's 5.2 s was under load). The land-cover step,
46.6 s for Kultsjön's 502 triangles and 51.3 s for Ljungan's 483, is 106a's
repair of self-crossing land-cover polygons (40 and 48 repaired), measured
on 106c's day (`../106c/README.md`); not 106e's. These times were measured with 106a's old, slow repair; master's corrected one (`4afd75fb`) came in with the merge at `d945ecc6`.

## Files

- `../106e-head-ac-r1/`, `../106e-base-ac-r1/`: `bench.py` output
  (`run.json`, `raw.tsv`, `README.md`), bench.py blob `16ad2ff5`.
- `*.stats.md`, `*.log`: each outline run's `--stats` output and stderr
  (the Ljungan land-cover log also holds `/usr/bin/time -l`).
- `power.txt`: `pmset -g batt` before and after each run.
- `orient_probe.py`: the 195-against-198 check, its command in its
  docstring.

The meshes were written to `rasputin_scratch/perf-106e/` and removed; rerun
the commands above to regenerate them.
