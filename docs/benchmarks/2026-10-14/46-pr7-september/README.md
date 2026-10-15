# Increment 46 PR 7: Bygdin to 1 September

The run `docs/increments/46-hbv-snow.md` section 11.6 asks of PR 7: section
10.9's Bygdin commands without `--days`, so at the new defaults (153 days,
1 April to 1 September 2027; the temperature at 300 m rises from -3 to +14 °C
over 122 days, to 1 August, then holds; one bin of snow per triangle), shadows
on and off, 8 threads. No refine or mesh code is touched, so this is a
demonstration, not an acceptance run: there is no verdict. The nearest earlier
figures are PR 6's 122-day, ten-bin run (`../46-pr6-acceptance/`); they differ
in run length and bins, so they are context, not a baseline.

## Wall time and memory

Median of 3 runs, shadows on and off interleaved; range in brackets. Memory is
`/usr/bin/time -l`'s maximum resident set size (MB = 10^6 B).

| case | days | bins | threads | wall, s | max resident, MB | frames | output, B |
|---|---:|---:|---:|---:|---:|---:|---:|
| Bygdin, shadows on | 153 | 1 | 8 | **29.56** (29.09-29.94) | 601.5 (594.5-611.5) | 154 | 65 318 365 |
| Bygdin, shadows off | 153 | 1 | 8 | **9.96** (9.77-10.09) | 149.5 (147.8-149.7) | 154 | 65 318 231 |
| PR 6, shadows on | 122 | 10 | 8 | 27.49 | 596.4 | 123 | 52 524 563 |
| PR 6, shadows off | 122 | 10 | 8 | 11.14 | 153.7 | 123 | 52 524 429 |

Section 11.6 expected about 30 s with shadows, less with one bin, 154 frames
and 63.6 MB of frames: measured 29.56 s, 154 frames, 63 556 416 B of `swe`
and `sca` frames (154 x 2 x 206 352) of the 65.3 MB. With shadows off, 31 more
days run 11 % faster than PR 6's ten-bin run; why (the bins or something else)
was not measured.

## Snow cover

`cover.py`, output `raw/cover.json`. Day 0 is 1 April; day 122 is 1 August,
day 153 is 1 September. One bin, so every triangle's `sca` is 0 or 1. Map
area 304.9 km², 51 527 triangles.

| | shadows on | shadows off |
|---|---:|---:|
| covered on 1 August, share of map area | **0.411** | **0.399** |
| covered on 1 August, share of triangles | 0.649 | 0.640 |
| mean SWE on 1 August, mm (by area) | 84.7 | 81.6 |
| covered on 1 September, share of map area | **0.043** | **0.040** |
| covered on 1 September, share of triangles | 0.157 | 0.149 |
| mean SWE on 1 September, mm (by area) | 6.8 | 6.3 |
| triangles that never melt out | 8 083 | 7 675 |
| last melt-out day | 153 (1 September) | 153 (1 September) |
| median melt-out day, triangles that melt | 126 | 125 |

- **The snow does not all melt by 1 September**: 4.3 % of the map area keeps
  snow with shadows (4.0 % without), up to 590 mm on a triangle; triangles are
  still melting out on the last frame.
- Small triangles crowd the high steep ground, so the share by count is
  higher than by area on both dates.
- Shadows on minus off: SWE never lower with shadows (0 to 513 mm on
  1 August, 0 to 543 mm on 1 September); 10 256 triangles melt out later with
  shadows, by up to 40 days.

## ParaView

ParaView 6.1.1's `pvpython` (`check_paraview.py`, copied from PR 6's with a
date label; output `raw/paraview.out`) read both runs: 154 time steps each
(0 to 3672 hours), 51 527 cells in every frame, the six cell arrays.

`bygdin-swe-day137.png`, day 137 (16 August), from above, north up: left SWE
with shadows on, middle off (0 dark blue to 600 mm white), right on minus off
(0 white to 250 mm dark red; range 0 to 552 mm). The valley and most slopes
are bare; snow is left on the high ground of the northern ridges; the
difference is near 0 almost everywhere, with red patches on steep slopes of
the northern ridges.

## The slow test

`RASPUTIN_SLOW=1 pytest tests/python/test_cli_snow.py::TestBygdin`, once,
same venv, on mains power: 3 passed in 37.5 s, including
`test_153_days_of_daily_frames` (154 daily frames).

## Machine and power state

Apple M1 Max, 10 cores, 32 GB, macOS 27.0; **on mains power** for every run
(`pmset -g batt` before each, column `power` of `raw/runs.tsv`). Nothing else
of ours ran. The one-minute load average at each start rose from 4.5 to 13.2
over the six runs (`runs.tsv`); after them the top processes were a WebKit
page and Spotlight, outside this work. The three runs of each case agree
within 3 %. Python 3.14.7, NumPy 2.5.4 (`raw/pip-freeze.txt`).

## Inputs and method

- Commit `cde378160a0e2e60bdc3adf765dc4525cd664afc` (branch
  `worktree-snow-summer`), installed into a fresh venv with `uv pip install
  '<worktree>[codecs]' -C cmake.build-type=Release click`; bounds checks on
  (the default). The installed Python files were checked equal to the
  branch's.
- The mesh is PR 6's `bygdin.vtk`, reused: sha1
  `ef5b03bbd976a19b643f9eba41191ddebea77b10`, 51 527 triangles, made by
  `rasputin mesh --dem $DTM10 --domain
  docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson --tolerance 10
  --out bygdin.vtk`; `$DTM10` = `../rasputin_data/DTM10_UTM33_20260925`.
- `run_pr7.sh <venv bin> bygdin.vtk $DTM10 <on dir> <off dir> raw 3` ran
  `rasputin snow --mesh bygdin.vtk --dem $DTM10 --shadows on|off --threads 8
  --out <dir>`, every other flag at its default (`raw/snow-on.json` and
  `raw/snow-off.json` record them: 153 days, warming 122, factors `[1.0]`).
- `python cover.py raw/cover.json 0 122 153 -- on=<on dir> off=<off dir>`.
- `pvpython check_paraview.py <on dir> <off dir> bygdin-swe-day137.png 137 "16 August"`.

The outputs are not in git: `../rasputin_scratch/snow-demo/bygdin-sept-shadows-on/`
and `bygdin-sept-shadows-off/` beside the main checkout (65 MB each; open
`snow.case`), from the third repeat. The commands above regenerate them.
