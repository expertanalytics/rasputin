# Increment 46 PR 3 acceptance: `rasputin snow` on the 1 m mesh

The run that `docs/increments/46-hbv-snow.md` sections 8 and 10.7 ask of PR 3
(the snow run and its animation): the default synthetic spring on the 1 m
mesh, at 1 and at 8 threads, with wall time, peak memory and output size, and
the animation checked in ParaView. No refine or mesh code is touched, so the
refine benchmark and thread sweep of `docs/increments/README.md` do not apply,
and there is no earlier `rasputin snow` measurement to compare with: these are
the first figures.

## Result

Median of 3 runs each, interleaved (1, 8, 1, 8, 1, 8 threads):

| threads | wall, s | max resident set | peak memory footprint | output |
|---:|---:|---:|---:|---:|
| 1 | **78.75** (78.61-80.05) | 430.6 MB (410.7 MiB) | 408.7 MB | 498 148 389 B (475.1 MiB), 191 files |
| 8 | **37.32** (37.26-37.46) | 431.1 MB (411.1 MiB) | 409.3 MB | 498 148 389 B (475.1 MiB), 191 files |

- Speed-up 1 to 8 threads: **2.11 times**. CPU time (user) is about the same
  at both counts (82-84 s at 1 thread, 85-86 s at 8), so most of the run does
  not use the extra threads; which part stays serial is not measured here (no
  profile was taken).
- Memory does not grow with threads (the difference is under 0.2 %).
- **The output does not depend on the thread count**: `diff -r` of the 1- and
  8-thread directories differs only in `snow.json`'s recorded `threads`.
- Output size against the design's estimate (section 4.6, 487 796 213 B for
  the geometry and two values over 92 frames): the extra 10 352 176 B are the
  four per-run variables (`slope_deg`, `aspect_deg`, `radiation_anomaly_mean`,
  `melt_out_day`, 2 587 724 B each), `snow.json` and the longer case file.
  In all: geometry 11 654 352 B, 188 variable files of 2 587 724 B,
  `snow.case` 1 039 B, `snow.json` 886 B.
- The program's own line, "92 frames in N s", reads 78.2-79.7 s and
  36.9-37.1 s (`raw/time-*.txt`, first line); the table's wall time adds
  process start and import.

## ParaView

ParaView 6.1.1's `pvpython` (`check_paraview.py`, output `raw/paraview.out`)
opened `snow.case` from the 8-thread run: 92 time steps, 0 to 2184 hours;
cell arrays `swe`, `sca`, `slope_deg`, `aspect_deg`,
`radiation_anomaly_mean`, `melt_out_day`; every one of the 92 frames read
with 646 870 cells. SWE range: frame 0, 600 to 600 mm; frame 45 (day 45,
16 May), 149.4 to 600 mm; frame 91 (day 91), 0 to 592.5 mm. The openvkl
initialisation errors ParaView prints (a volume-rendering module missing from
the app bundle) do not touch the reader or the screenshot.

`swe-day45.png`: frame 45, SWE per triangle, 0 (dark blue) to 600 mm (white),
seen obliquely from above. The massif keeps its snow; the lower ground is
bluer, with less.

## Machine and power state

Apple M1 Max, 10 cores, 32 GB, macOS 27.0; **on mains power** (`pmset -g
batt`: "AC Power", 100 %, charged, before and after). Nothing else of ours
ran; load average 1.79 at the start, 3.76 at the end. Python 3.14.7, NumPy
2.5.4.

## Inputs

- Commit `914e1f574f5d08f90ed601d3aa87982499f03e8b` (branch
  `worktree-snow-anim`), the extension rebuilt Release from it
  (`build-pyext`, `cmake --build build-pyext -j --target _core`) and copied
  into the venv before the runs.
- `rasputin synthetic-terrain --out synth.tif` (defaults; sha1
  `1ec534188aefa351b2d1a6c61b9d53b77cda8b09`), then
  `rasputin mesh --dem synth.tif --tolerance 1 --out mesh1.vtk`: 646 870
  triangles, sha1 `5143d0d3d34d0141599efec0e1a891b5ec80066f` (byte-identical
  to the demo's `mesh1.vtk`).
- `rasputin snow --mesh mesh1.vtk --out out-t<N> --threads <N>`, every other
  flag at its default: melt-factor sun correction, 2027-04-01 for 91 days,
  hourly steps, a frame every 24 hours, 600 mm initial SWE, CV 0.4, 10 bins
  (the full job in `raw/snow-t8.json`).

## Method and files

- `run_snow.sh <rasputin> <work dir> 3` made every run, each under
  `/usr/bin/time -l` (wall time, maximum resident set size, peak memory
  footprint); output size is the sum of the output directory's file sizes.
  It wrote `runs.tsv` with a literal `\t` between fields (`print -r`); the
  committed `runs.tsv` is that file with the `\t` turned into tabs, nothing
  else changed.
- `raw/time-t<N>-r<i>.txt`: each run's program line and `time -l` report.
- `check_paraview.py`: run as
  `pvpython check_paraview.py <dir>/snow.case swe-day45.png 45`.

The outputs (about 475 MiB each) are not in git. They are in
`/Users/skavhaug/projects/rasputin_scratch/snow-demo/acceptance/` (`out-t1/`,
`out-t8/`, with `synth.tif` and `mesh1.vtk`); the three commands above
regenerate them.
