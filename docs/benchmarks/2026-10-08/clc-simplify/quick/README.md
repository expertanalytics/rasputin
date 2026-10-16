# Increment 32, quick check: base `b39426c0` against branch `d6f04b9e`

The short timed check of `docs/increments/32-landcover-simplify.md`, section
10: `tools/bench_quick.py`, base and branch back to back, 2026-10-08.

**Verdict: SLOWER (exit 1).** Both catchment cases take about 70 % longer.
The new simplifier itself costs 0.66 s on Numedalslagen, in line with the
design's estimate. The added time is in the outline rule (`snap_to_outline`,
which moves borders within 5 m of the outline onto it). The branch does not
change that code, but it takes about 4.4 times as long on the simplified
borders. Measured with a profiler on Numedalslagen only; why it is slower
was not measured.

## Method

- Machine: Apple M1 Max (8 performance + 2 efficiency cores), 32 GiB, macOS
  27.0, Python 3.14.7, NumPy 2.5.3, shapely 2.2.0, GEOS 3.14.1, pyproj 3.8.0,
  PROJ 9.8.1. **Power: AC**, the same before and after each run (pmset).
- Tool: `tools/bench_quick.py` and `tools/bench.py` as on master `4cd7e050`
  (blobs `5164c4fe` and `237f1c91`, with the GEOS and PROJ version fields of
  #219), run from a `git archive` copy so the baseline went to a scratch
  folder. Cases: `docs/benchmarks/quick/cases.toml` as at `4cd7e050` (the
  same as on base and branch). Data: `../rasputin_data`.
- Base: a detached worktree of `b39426c0`. Branch: the worktree at
  `d6f04b9e`. Both built Release with bounds checks on (`libc++ fast`) into
  `<tree>/build-bench`, and checked again by the tool at the start of its run.
- Order: base with `--save-baseline` (18:41 to 18:43 CEST), then the branch
  against it (18:43 to 18:47). Runs and warm-ups as in `cases.toml`: tile and
  quarter 5 runs, the two catchments 1 warm-up and 3 runs. Medians, with
  min and max in brackets.
- Triangles and angles: one more `--stats` run per catchment per tree
  (`stats-*.md`). Their totals are within 2.5 % of the quick-check medians.
- Profile: one `cProfile` run of Numedalslagen per tree
  (`profile-numedalslagen.txt`).
- Raw records: `base-b39426c0.json` and `branch-d6f04b9e.json`. The meshes
  were written to a scratch folder and deleted. To make them again, run
  `rasputin mesh` with the case's arguments from `cases.toml` and
  `--stats <file>`.

## Whole run and phases (seconds, median)

The table shows phases at 5 % or more of either run.

| case | measure | base | branch | change | share base | share branch |
|---|---|---|---|---|---|---|
| tile | total | 0.42 | 0.42 | -0.6 % | | |
| tile t=1 | total | 0.72 | 0.72 | -0.4 % | | |
| quarter | total | 0.35 | 0.35 | -0.6 % | | |
| numedalslagen | **total** | 9.92 (9.90-9.93) | **16.95** (16.94-17.01) | **+70.8 %** | 100 % | 100 % |
| numedalslagen | features clip | 5.25 | 14.04 | +167.3 % | 53 % | 83 % |
| numedalslagen | features clip: clean-up | 5.19 | 13.97 | +169.4 % | 52 % | 82 % |
| numedalslagen | refine | 1.20 | 0.67 | -44.7 % | 12 % | 4 % |
| numedalslagen | refine: start quality | 0.62 | 0.04 | -94.1 % | 6 % | 0 % |
| numedalslagen | land cover | 0.58 | 0.30 | -47.5 % | 6 % | 2 % |
| numedalslagen | decode | 1.01 | 1.01 | +0.5 % | 10 % | 6 % |
| lagan | **total** | 19.52 (19.48-19.55) | **34.15** (33.81-34.98) | **+75.0 %** | 100 % | 100 % |
| lagan | features clip | 9.72 | 25.82 | +165.6 % | 50 % | 76 % |
| lagan | features clip: clean-up | 5.65 | 21.67 | +283.4 % | 29 % | 63 % |
| lagan | decode | 4.56 | 4.57 | +0.2 % | 23 % | 13 % |

Tile and quarter have no land cover and are unchanged: same mesh hash, all
phases within the 5 % band. On Lagan, refine also fell from 0.66 s to
0.24 s (`--stats` runs), below the 5 % cut-off of the table.

## Where the clean-up time goes (Numedalslagen, under the profiler)

| step | base | branch |
|---|---|---|
| `_clean` (the clean-up phase) | 5.76 s | 16.09 s |
| of it, `snap_to_outline` (the 5 m outline rule) | 2.76 s | **12.18 s** |
| of it, `simplify_borders` (new) | - | 0.66 s |
| of it, `coverage_clean` | 1.16 s | 0.50 s |

The profiler adds 11 to 15 % to the clean-up time (5.76 s against 5.19 s, 16.09 s against 13.97 s). Inside `snap_to_outline`,
`_loops` is called 503 times instead of 373, and `make_valid` 245 times
instead of 3. Lagan was not profiled.

## Meshes (reported, not judged: CORINE borders change by design)

| | Numedalslagen base | Numedalslagen branch | Lagan base | Lagan branch |
|---|---|---|---|---|
| land-cover vertices after clean-up | 598 595 | 76 220 | 385 640 | 103 572 |
| start triangles | 281 483 | 85 615 | 326 646 | 151 755 |
| **output triangles** | **1 118 006** | **654 920** (-41 %) | **798 554** | **464 119** (-42 %) |
| smallest angle, median | 36.09° | 33.42° | 35.63° | 32.46° |
| triangles under 1° | 0.00 % | 0.01 % | 0.00 % | 0.00 % |
| triangles under 10° | 1.18 % | 3.28 % | 4.31 % | 2.44 % |
| **worst angle** | 0.832° | **0.0157°** | 0.729° | **0.000658°** |
| max vertex degree | 19 | 18 | 14 | 16 |
| vertices of degree 12 or more | 153 | 630 | 85 | 259 |
| largest height error (tolerance 10 m) | 9.99994 | 9.99999 | 9.99999 | 9.99996 |
| area the outline rule moved, m² | 13 238 | 1 097 391 | 16 050 | 1 531 067 |

The tolerance check passed on every case. A Delaunay check was not run:
`bench_quick.py` does not make one.

## Hotspots (a phase at 40 % or more of a run)

- Numedalslagen: features clip 83 %, its clean-up 82 % (base: 53 % and 52 %).
  Not in `hotspots.toml`.
- Lagan: features clip 76 % (ruled at 41 %, now 10 points or more above
  that), its clean-up 63 % (not in `hotspots.toml`).
- In both runs, unchanged: tile refine 45 %, tile t=1 refine 68 % (scan 50 %),
  quarter refine 47 to 48 %. Not in `hotspots.toml`.
