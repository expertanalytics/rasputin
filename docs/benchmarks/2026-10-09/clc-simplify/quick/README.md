# Increment 32 fix, quick check: base `b39426c0` against branch `8605f607`

The short timed check of `docs/increments/32-landcover-simplify.md`, re-run on
the fix of section 15 (green step `05c8f73a`, test fix `8605f607`) with
`tools/bench_quick.py`. Base and branch were run back to back on 2026-10-09.
The earlier check, on `d6f04b9e`, is `../../../2026-10-08/clc-simplify/quick/`.

**Tool verdict: SLOWER (exit 1).** The tool reports SLOWER because the
"features clip" phase grew by more than its 5 % band on both catchments
(+16 % Numedalslagen, +15 % Lagan). Whole runs: Numedalslagen is **11 %
faster**, and Lagan is unchanged (-1.3 %, inside the band). Tile and quarter
are unchanged, with the same mesh hash.

**The fix does what section 15.5 expected for time and moved area. It does not
meet the worst-angle expectation on Numedalslagen: 0.311°, under the 0.4°
floor.** Per section 15.5, that goes back to `@architect`. Lagan's worst
angle, 0.49°, meets the floor.

## Method

- Machine: Apple M1 Max (8 performance + 2 efficiency cores), 32 GiB, macOS
  27.0, Python 3.14.7, NumPy 2.5.3, shapely 2.2.0, GEOS 3.14.1, pyproj 3.8.0,
  PROJ 9.8.1. **Power: AC**, the same before and after every run (pmset,
  recorded in both JSON records).
- Tool: `tools/bench_quick.py` and `tools/bench.py` as on master `4cd7e050`
  (blobs `5164c4fe` and `237f1c91`), run from a `git archive` copy so that the
  base record went to a scratch folder. This is the same method as the
  2026-10-08 check. Cases: `docs/benchmarks/quick/cases.toml` blob `48f20ce6`
  (the same at `4cd7e050` and at `8605f607`). Data: `../rasputin_data`.
- Base: a detached worktree of `b39426c0`. Branch: the `clc-simplify` worktree
  at `8605f607`. The tool built both Release with bounds checks on
  (`libc++ fast`) into `<tree>/build-bench`, and checked them again at the
  start of its run.
- Order: base with `--save-baseline` (00:58 to 01:00 CEST), then the branch
  against it twice (01:00 to 01:03, and 01:03 to 01:05 with `--save-baseline`
  to keep its record). The tables use the second branch run
  (`branch-8605f607.json`). The first run gave Numedalslagen 9.05 s (-12.5 %)
  and the same verdict lines (`run-branch.log`).
- Runs and warm-ups as in `cases.toml`: tile and quarter 5 runs, the two
  catchments 1 warm-up and 3 runs. Medians, with min and max in brackets.
- Triangles and angles: one more `--stats` run per catchment per tree
  (`stats-*.md`), 01:05 to 01:07.
- Profile: one `cProfile` run of Numedalslagen per tree
  (`profile-numedalslagen.txt`).
- The meshes were written to a scratch folder and deleted. To make them
  again, run `bench.py _child --pkg <tree>/build-bench/pkg --threads 0 --` with
  the case's arguments from `cases.toml`, plus `--out <file> --stats <file>`.

## Section 15.5, expected against measured

| measure | base | expected after the fix | measured on `8605f607` | met |
|---|---|---|---|---|
| Numedalslagen whole run | 10.34 s | about +8 % | 9.18 s, **-11.2 %** | yes, better than expected |
| Lagan whole run | 20.10 s | within about 10 % of base | 19.84 s, -1.3 % | yes |
| outline rule (`snap_to_outline`, under the profiler, Numedalslagen) | 2.83 s | base's time | 2.77 s; same call counts as base (373 `_loops`, 3 `make_valid`) | yes |
| clean-up phase, Numedalslagen | 5.29 s | about base + 1.6 s | 6.16 s, base + 0.87 s | yes |
| clean-up phase, Lagan | 5.78 s | base + 1.5 to 2 s | 7.29 s, base + 1.51 s | yes |
| moved area, Numedalslagen | 13 238.0 m² | exactly base's | 13 238.0 m² | yes, to the 0.1 m² the stats print |
| moved area, Lagan | 16 050.3 m² | exactly base's | 16 050.3 m² | yes, to the 0.1 m² the stats print |
| **worst angle, Numedalslagen** | 0.832° | 0.4° or more | **0.311°** | **no** |
| worst angle, Lagan | 0.729° | 0.4° or more | 0.49° | yes |
| output triangles, Numedalslagen | 1 118 006 | about -40 % | 653 039, -41.6 % | yes |
| output triangles, Lagan | 798 554 | about -40 % | 463 342, -42.0 % | yes |

On Numedalslagen, base's decode phase was 0.14 s slower than the branch's
(1.18 s against 1.04 s). The branch does not change decoding, so about 0.14 s
of the 1.16 s gain is probably run-to-run noise (base ran first). Not
measured further.

## Whole run and phases (seconds, median)

The table shows phases at 5 % or more of either run.

| case | measure | base | branch | change | share base | share branch |
|---|---|---|---|---|---|---|
| tile | total | 0.44 | 0.43 | -1.7 % | | |
| tile t=1 | total | 0.75 | 0.72 | -4.7 % | | |
| quarter | total | 0.37 | 0.35 | -3.4 % | | |
| numedalslagen | **total** | 10.34 (10.24-10.41) | **9.18** (9.09-9.29) | **-11.2 %** | 100 % | 100 % |
| numedalslagen | decode | 1.18 | 1.04 | -12.3 % | 11 % | 11 % |
| numedalslagen | features clip | 5.35 | 6.22 | +16.3 % | 52 % | 68 % |
| numedalslagen | features clip: clean-up | 5.29 | 6.16 | +16.5 % | 51 % | 67 % |
| numedalslagen | refine | 1.29 | 0.67 | -48.4 % | 12 % | 7 % |
| numedalslagen | refine: start quality | 0.69 | 0.04 | -94.9 % | 7 % | 0 % |
| numedalslagen | land cover | 0.61 | 0.30 | -50.5 % | 6 % | 3 % |
| lagan | **total** | 20.10 (20.03-20.12) | **19.84** (19.75-20.08) | **-1.3 %** | 100 % | 100 % |
| lagan | decode | 4.65 | 4.67 | +0.5 % | 23 % | 24 % |
| lagan | features clip | 9.96 | 11.47 | +15.2 % | 50 % | 58 % |
| lagan | features clip: clean-up | 5.78 | 7.29 | +26.2 % | 29 % | 37 % |

Tile and quarter have no land cover. Their mesh hashes match base's, and
every phase is inside the 5 % band.

## Where the clean-up time goes (Numedalslagen, under the profiler)

| step | base | branch |
|---|---|---|
| `_clean` (the clean-up phase) | 5.90 s | 6.74 s |
| of it, `snap_to_outline` (the 5 m outline rule) | 2.83 s | 2.77 s |
| of it, `coverage_clean` | 1.18 s | 1.15 s |
| of it, `simplify_borders` (new) | - | 0.68 s |
| of it, shapely calls made from `_clean` itself | 0.77 s | 1.55 s |

## Meshes (reported; CORINE borders change by design)

| | Numedalslagen base | Numedalslagen branch | Lagan base | Lagan branch |
|---|---|---|---|---|
| land-cover vertices after clean-up | 598 595 | 75 014 | 385 640 | 102 441 |
| start triangles | 281 483 | 84 399 | 326 646 | 150 604 |
| **output triangles** | **1 118 006** | **653 039** (-41.6 %) | **798 554** | **463 342** (-42.0 %) |
| smallest angle, median | 36.09° | 33.47° | 35.63° | 32.47° |
| triangles under 1° | 0.00 % | 0.01 % | 0.00 % | 0.00 % |
| triangles under 10° | 1.18 % | 3.19 % | 4.31 % | 2.38 % |
| **worst angle** | 0.832° | **0.311°** | 0.729° | 0.49° |
| max vertex degree | 19 | 17 | 14 | 16 |
| vertices of degree 12 or more | 153 | 635 | 85 | 256 |
| largest height error (tolerance 10 m) | 9.99994 | 10 | 9.99993 | 9.99998 |
| area the outline rule moved, m² | 13 238.0 | 13 238.0 | 16 050.3 | 16 050.3 |

The tolerance check (largest error at most the tolerance) passed on every
case. On Numedalslagen the branch's largest error is printed as exactly 10.
A Delaunay check was not run: `bench_quick.py` does not make one.

## Hotspots (a phase at 40 % or more of a run)

- Numedalslagen: features clip 68 %, its clean-up 67 % (base: 52 % and 51 %).
  Neither is in `hotspots.toml`.
- Lagan: features clip 58 %. It was ruled at 41 % and is now 10 points or more
  above that, so the tool raises it again. Base today: 50 %.
- In both trees, unchanged: tile refine 44 to 46 %, tile t=1 refine 67 to 68 %
  (scan 50 to 51 %), quarter refine 47 to 48 %. None is in `hotspots.toml`.
