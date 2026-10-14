# Increment 46 PR 6 acceptance: `rasputin snow` with cast shadows

The run that `docs/increments/46-hbv-snow.md` sections 10.7 and 10.9 ask of
PR 6 (the evaluator, the near raster and the snow wiring): Bygdin on DTM10 and
the synthetic 2 m mesh, with shadows, at 1 and 8 threads: wall time,
horizon-search time, evaluator time per step, peak memory; the beam the
shadows remove on 2027-04-15, 05-15 and 06-21 against section 10.6's figures
on DEM nodes; one screenshot of Bygdin. No refine or mesh code is touched, so
the refine benchmark and thread sweep of `docs/increments/README.md` do not
apply. There is no earlier shadowed run to compare with (**no baseline**);
the comparison is with section 10.4's cost estimate.

## The whole run (`rasputin snow`, wall time and memory)

Median of 3 runs, interleaved; range in brackets. Memory is `/usr/bin/time
-l`'s maximum resident set size. "off" is one run each, `--shadows off`.

| case | triangles | days | threads | wall, s | max resident | output |
|---|---:|---:|---:|---:|---:|---:|
| Bygdin, shadows on | 51 527 | 122 | 1 | **176.37** (176.30-177.37) | 590.3 MB | 52 524 563 B |
| Bygdin, shadows on | 51 527 | 122 | 8 | **27.49** (27.47-27.61) | 596.4 MB | same |
| Bygdin, shadows off | 51 527 | 122 | 1 / 8 | 62.18 / 11.14 | 145.7 / 153.7 MB | 52 524 429 B |
| synthetic 2 m, shadows on | 294 770 | 91 | 1 | **507.76** (498.56-509.60) | 481.1 MB (439.7-481.5) | 227 028 196 B |
| synthetic 2 m, shadows on | 294 770 | 91 | 8 | **78.34** (77.65-78.42) | 481.6 MB | same |
| synthetic 2 m, shadows off | 294 770 | 91 | 1 / 8 | 168.38 / 31.46 | 445.4 / 447.0 MB | 227 028 062 B |

- Speed-up 1 to 8 threads with shadows: 6.4 (Bygdin), 6.5 (synthetic).
- **The output does not depend on the thread count**: every file of the 1- and
  8-thread runs is byte-identical except `snow.json` (which records the
  threads), for both cases.
- Shadows add about 440 MB of memory on Bygdin: the near raster is 7314 x
  6195 float32 nodes (181 MB) plus what the search builds on it; the split
  between the two was not measured.

## The parts (`components.py`)

The near raster's read, the horizon search and the evaluator, each timed
alone as `rasputin snow` calls them; the evaluator on all 24 hourly steps of
every 7th day of the spring (4 sun samples per step; 432 steps on Bygdin, 312
on the synthetic). Median of 3 runs; the three agree within 1 %.

| case | threads | near raster read, s | horizon search, s | evaluator, ms per step (95th pct) | evaluator, ns per triangle and sample |
|---|---:|---:|---:|---:|---:|
| Bygdin | 1 | 0.57 | **135.2** | **10.88** (13.9) | 55.3 |
| Bygdin | 8 | 0.57 | **18.7** | **1.82** (2.2) | 9.1 |
| synthetic 2 m | 1 | 0.01 | **330.0** | **62.17** (79.2) | 54.9 |
| synthetic 2 m | 8 | 0.01 | **45.5** | **9.66** (11.9) | 8.5 |

Search speed-up 7.2 (Bygdin) and 7.3 (synthetic) at 8 threads; evaluator 6.0
and 6.4. On Bygdin every ray ran its full 20 km (`reach` 255 everywhere); on
the synthetic DEM, a 12 km square, every ray left the raster first.

**Against section 10.4's estimate** (one thread):

| figure | estimate | measured | ratio |
|---|---:|---:|---:|
| horizon search, Bygdin | about 40 s before the pyramid prunes it (0.8 ms per triangle) | 135.2 s (2.62 ms per triangle) | 3.4 |
| evaluator, per triangle and sample | about 29 ns | 55 ns | 1.9 |
| evaluator, synthetic 2 m over 91 days | about 75 s | 135.8 s (62.17 ms x 2184 steps) | 1.8 |

**The one profile** (`search_reach.py`, output `raw/search_reach.out`), taken
because the search is 3.4 times its estimate: 2000 fixed random Bygdin
receivers, one thread, searched to 2.5, 5, 10 and 20 km on the same near
raster:

| search distance | µs per ray | ms per receiver |
|---:|---:|---:|
| 2.5 km | 17.5 | 1.12 |
| 5 km | 25.1 | 1.61 |
| 10 km | 33.0 | 2.11 |
| 20 km | 42.8 | 2.74 |

Each doubling of the distance adds about 8-10 µs per ray, so the pyramid does
cut the far part of the ray; but the first 2.5 km alone, 500 samples at 5 m,
cost 17.5 µs per ray, about 35 ns per sample if every sample was read, where
the estimate assumed 3 ns. Where those 17.5 µs go is not measured: a `sample`
profile of a one-thread search put almost all its samples in one loop of
`_core` of about 800 bytes, but the extension carries no symbols, so the loop
was not named. Section 10.4's 3 ns per sample is thought to be too low for a
bilinear read plus the pyramid's bookkeeping; a symbolised profile would say.

## What the shadows remove, on the mesh (`shadow_effect.py`)

As the run computes it: one receiver per triangle, 96 sun samples a day (the
midpoints of 15-minute intervals), the beam at the mesh's centre. Weighted by
each triangle's map area, the analogue of section 10.6's evenly spaced DEM
nodes and of test E2. Section 10.6's node figures in brackets.

| site and date | beam removed | mean change in the day's dI, W/m² | 5th percentile, W/m² | share losing more than 20 W/m² |
|---|---:|---:|---:|---:|
| Bygdin, 04-15 | 1.43 % (1.6) | -2.36 (-2.7) | -13.06 (-14.3) | 2.42 % (3.1) |
| Bygdin, 05-15 | 1.18 % (1.3) | -2.78 (-3.0) | -14.71 (-14.9) | 2.96 % (2.9) |
| Bygdin, 06-21 | 1.13 % (1.2) | -3.09 (-3.3) | -15.73 (-16.5) | 3.25 % (3.5) |
| synthetic 2 m, 04-15 | 0.61 % (0.7) | **-0.99** (-1.1) | -5.71 (-5.6) | 0.06 % |
| synthetic 2 m, 05-15 | 0.43 % (0.5) | -1.00 (-1.2) | -5.56 (-5.8) | 0.02 % |
| synthetic 2 m, 06-21 | 0.40 % (0.5) | -1.09 (-1.3) | -5.81 (-6.1) | 0.02 % |

- On the mesh the beam removed and the mean change are 6-17 % smaller than
  on the nodes, on every date; the 5th percentile and the share losing more
  than 20 W/m² are within a quarter of the node figures, either way. The probe
  and the run differ in normals (DEM central differences against triangles at
  10 m tolerance), refraction (none in the probe) and the receivers' places.
  Which of these makes the gap is not measured.
- Counting each triangle once instead (small triangles crowd the steep
  ground), Bygdin's beam removed is 2.51 / 2.22 / 2.16 % and the mean change
  -4.12 / -4.94 / -5.43 W/m² (`raw/shadow_effect-bygdin.json`, "triangle").
- **E2's figure**: the plan-area-weighted domain mean change of the day's dI
  on the 2 m mesh on 2027-04-15 is **-0.99 W/m²**; unweighted, -1.04.
- Over the whole 122-day Bygdin run, the plan-weighted mean of
  `radiation_anomaly_mean` falls by 2.78 W/m² with shadows; no triangle's is
  higher with shadows. On 7.9 % of the map area the melt-out day is later
  with shadows, by up to 27 days. On the synthetic 91-day run: -1.03 W/m²,
  9.7 %, up to 6 days.
- **Much of Bygdin keeps its snow to 1 August** under the synthetic spring:
  32 949 of 51 527 triangles never fall below half cover (32 448 with shadows
  off); the mean SWE is 232 mm of 600 on day 120.

## ParaView

ParaView 6.1.1's `pvpython` (`check_paraview.py`, output `raw/paraview.out`)
read both Bygdin runs: 123 time steps each (0 to 2928 hours), 51 527 cells in
every frame, the six cell arrays.

`bygdin-swe-day110.png`, day 110 (20 July), seen from above, north up, with
the terrain's own shading: left, SWE with shadows on; middle, shadows off
(0 dark blue to 600 mm white: Bygdin lake's valley floor nearly bare, the
peaks and the northern ridge still white); right, their difference, on minus
off, 0 white to 250 mm dark red. The difference is never negative (range 0
to 489 mm); it is near 0 over most of the catchment and shows as red streaks
on the steep slopes of the ridges, most along the northern ridge, where
terrain shades the slopes and keeps snow.

## Machine and power state

Apple M1 Max, 10 cores, 32 GB, macOS 27.0; **on mains power** for every run
(`pmset -g batt` before each, column `power` of `raw/runs.tsv`: "AC Power",
100 %, charged). Nothing else of ours ran beside a timed run; the one-minute
load average at each start is in `runs.tsv` (1.4 to 7.5; the higher figures
come from processes outside this work, which were not identified; the three Bygdin 1-thread runs, at
loads 2.6, 6.4 and 7.5, agree within 0.6 %). Python 3.14.7, NumPy 2.5.4.

## Inputs

- Commit `7cf50a9bb8fbbb84b590733dd9d4966f8b6864049`; its own venv, the
  extension built Release with bounds checks on (`RASPUTIN_HARDENING=ON`, the
  default) by `uv pip install '<worktree>[codecs]' -C
  cmake.build-type=Release`, plus `click` 8.5.0 (see below).
- `$DTM10` = `../rasputin_data/DTM10_UTM33_20260925`.
- `rasputin mesh --dem $DTM10 --domain
  docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson --tolerance 10
  --out bygdin.vtk`: 51 527 triangles, sha1
  `ef5b03bbd976a19b643f9eba41191ddebea77b10`.
- `rasputin synthetic-terrain --out synth.tif` (sha1
  `1ec534188aefa351b2d1a6c61b9d53b77cda8b09`), `rasputin mesh --dem
  synth.tif --tolerance 2 --out mesh2.vtk`: 294 770 triangles, sha1
  `2280b14621fdc411e06d38e8ae6535ed515714d7`.
- Bygdin: `rasputin snow --mesh bygdin.vtk --dem $DTM10 --shadows on|off
  --days 122 --threads N`; synthetic: `rasputin snow --mesh mesh2.vtk --dem
  synth.tif --shadows on --threads N` (off: no `--dem`); every other flag at
  its default (melt-factor correction, from 2027-04-01, hourly steps, 4 sun
  samples per step, daily frames).

## Method and files

- `run_pr6.sh <venv bin> <work dir> $DTM10 raw 3` made every run above, in
  this order: the snow runs interleaved by thread count, the off runs, the
  part timings, the profile. `raw/runs.tsv`, `raw/time-*.txt` (`time -l`),
  `raw/log-*.txt`, `raw/components-*.json`, `raw/search_reach.*`.
- `shadow_effect.py <mesh> <dem> <out.json>`, run after the timings, eight
  threads: `raw/shadow_effect-*.json` (all three weightings) and `.out`.
- `check_paraview.py`: `pvpython check_paraview.py bygdin-shadows-on
  bygdin-shadows-off bygdin-swe-day110.png 110`.

The outputs are not in git. The 8-thread Bygdin runs, for ParaView, are
`../rasputin_scratch/snow-demo/bygdin-shadows-on/` and
`bygdin-shadows-off/` (51 MB each; open `snow.case`). The meshes, `synth.tif`
and the other runs' last outputs are in `../rasputin_scratch/snow-demo/pr6/`.
The commands above regenerate them.

**A fresh install fails without `click`**: `tin_engine/cli.py` imports
`click`, which `pyproject.toml` does not declare and typer 0.27.3 no longer
brings; the venv here needed `click` installed by hand. It predates PR 6
(the import is on `master`).
