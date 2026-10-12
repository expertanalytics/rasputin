# 26b acceptance: exact land-cover class fractions per triangle (@perf, 2026-10-12)

Increment 26b (`docs/increments/26-land-cover.md`, section 14 as built; the
acceptance is section 11, item 1, scoped to what 26b builds: no cutoff, no
ledger, no Hilbert order, no water). Two gates:

1. **Exactness and cost on a basin piece**: the Velhas piece meshed at 10 and
   5 m with `--land-cover mapbiomas-c11 --land-cover-year 2025`.
2. **No change to existing commands**: the quick check (`tools/bench_quick.py`).

Verdicts: **gate 1 ACCEPTED**; **gate 2: no change attributable to 26b**, by a
back-to-back run against the branch's base commit, because the stored
baseline's comparison was not valid today (a foreign load on the machine; the
untouched base commit reads SLOWER against it by the same amounts). Details
in "Gate 2".

## Method

- **Tree**: branch `worktree-mapbiomas-26b` at `4149b50d` (26b green
  `5b8bb3e0`, test fixes `f165282b`). Base for the back-to-back run:
  `5f28292c`, the merge of 26a, checked out as a scratch worktree and
  removed afterwards.
- **Build**: Release `_core` built by the same configure `tools/bench.py`
  uses (`-DCMAKE_BUILD_TYPE=Release -DRASPUTIN_HARDENING=ON`, flags
  `-O3 -DNDEBUG`, AppleClang 21.0.0) into the worktree's `build-bench`, then
  copied into the worktree's venv (`.venv`); sha256 of the `.so`
  `55bb3993519147fe…`. The quick check rebuilt the same directory (a no-op)
  and ran from `build-bench/pkg`.
- **Machine**: Apple M1 Max (8 P + 2 E cores, 32 GiB), macOS 27.0, Python
  3.14.7, numpy 2.5.4, shapely 2.2.0 (GEOS 3.14.1), pyproj 3.8.0 (PROJ
  9.8.1), tifffile 2026.9.20, imagecodecs 2026.10.10 (installed into the
  venv for this run: MapBiomas is LZW with the horizontal predictor).
- **Power**: AC, battery 100 % charged, before and after every block
  (`raw/pmset.log`, `raw/quick_pmset_*.txt`, `raw/quick_b2b_pmset.log`,
  `raw/analysis_pmset.log`).
- **Foreign load, measured**: a Microsoft Teams call ran on the machine for
  the whole session (`Microsoft Teams ModuleHost` 40 to 50 % of a core,
  `WindowServer` 28 to 33 %, `coreaudiod` 12 %, the camera daemon 10 %;
  load average 5.3 to 7.6 on 10 cores; `raw/analysis_pmset.log`,
  `raw/quick_b2b_pmset.log`). I was not blocked by it, but every timing
  below was taken under it, and it is why the stored quick-check baseline
  could not be used (Gate 2).
- **Data**: DEM `rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif`
  (ANADEM resampled to SIRGAS 2000 / UTM 23S at 30 m, the grid the existing
  `velhas76949_t*.bin.vtk` runs of `docs/benchmarks/2026-10-02/basin-piece-anadem/`
  used; their `run_sweep.py` and the meshes' headers say so), domain
  `bho2017_5k_76949_outline_epsg4674.geojson` (one ring, 7,309 vertices,
  11,667.6 km²), tolerances 10 and 5 m, no other option. Land cover:
  MapBiomas Collection 11, year **2025** (the design's default and Ola's Q1
  ruling; the brief's 2023 gives way to the design), fetched with
  `rasputin fetch mapbiomas-c11 --year 2025 --domain <outline> --domain-crs EPSG:4674 --out-crs EPSG:31983 --cache rasputin_data/cache`:
  540 blocks, 5,077,805 bytes, 10.2 s, now in `rasputin_data/cache/mapbiomas-c11/brazil_coverage-col11_2025/`.
- **Runs**: `run_meshes.sh` (this folder): each tolerance with and without
  `--land-cover`, three runs each, under `/usr/bin/time -l` for the peak
  resident set; `--stats` per run in `raw/*.stats.md`, stderr in `raw/*.log`.
  Meshes in `rasputin_scratch/sao_francisco_piece_meshes/26b-fractions/`
  (not committed; `run_meshes.sh` regenerates them).
- **Analysis**: `analyse.py` (this folder; its docstring is the method),
  output `raw/analysis_t10.txt`, `raw/analysis_t5.txt`,
  `raw/analysis_t10.json`, `raw/analysis_t5.json`. The reference is
  independent of the producer: the cached raster read back through the cache
  reader, cell centres tested against the outline with GEOS
  (`shapely.contains_xy`), the cells the outline's boundary crosses found by
  GEOS (`shapely.intersects` of each outline segment with the cells of its
  box), cell areas by `pyproj.Geod` per row on the WGS 84 ellipsoid. The
  outline is moved from EPSG:4674 to the raster's EPSG:4326 by pyproj.
- **Threads**: the CLI's fractions pass runs `_core.class_areas` with
  `threads=0`, all 10 cores (`cli._class_fractions` passes no thread count;
  there is no `--threads` option on `mesh`). The sweep by 1, 2, 4, 8 and 10
  threads calls `land_cover_fractions.assign`, the function the CLI calls,
  on the mesh read back from the file, median of 3; `analyse.py` checks that
  its output equals the file's arrays bit for bit before the numbers count.

## Gate 1: exactness and cost

### The meshes

| tolerance | vertices | triangles | class entries | classes | per triangle, mean | max | min angle worst | vertex degree max | max error | file with / without land cover |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 m | 328,469 | 649,504 | 1,299,387 | 19 | 2.001 | 8 | 0.113° | 14 | 10 m | 35.3 / 23.6 MB |
| 5 m | 789,884 | 1,571,955 | 2,690,575 | 19 | 1.712 | 7 | 0.113° | 14 | 5 m | 82.9 / 56.8 MB |

With and without `--land-cover`, the points, triangles, constraint lines and
feature masks are equal (`np.array_equal` on the files read back, both
tolerances); the quality rows of the `--stats` reports are identical. The
pass changes no geometry. (These meshes differ from the 2026-10-02 ones,
649,845 and 1,576,681 triangles, which were built from `5a57793`; that is
the refine changes of the increments since, not 26b, which touches no refine
code.)

Classes per triangle, 10 m: 1 class 227,784 triangles, 2: 242,707, 3:
136,907, 4: 35,779, 5: 5,660, 6: 622, 7: 40, 8: 5; 5 m: 739,185 / 585,474 /
211,858 / 32,497 / 2,769 / 169 / 3 (`raw/analysis_t*.json`).
The design's section 5 estimate of 1.50 at 10 m and 1.36 at 5 m is *after
the 5 % cutoff and the 1 ha floor* (26c); 2.00 and 1.71 here are before any
cutoff, so not comparable yet.

### Class areas against MapBiomas's own count

Identical for the two meshes to the fourth decimal of a km² (as they must
be: exact fractions of the same domain). `mesh` is the sum over triangles of
the exact class areas; `ref` is the area of the cells whose centres fall
inside the outline (13,819,778 cells); the band is `[interior, touching]`:
the cells the boundary does not cross (`interior`) up to those plus every
cell it crosses (33,945 cells, 28.7 km² in all). The exact area of a class
inside the outline lies in its band; the pass test is `mesh` inside the band
widened by max(0.1 % of `ref`, 1 ha). The outline's ellipsoidal area is
11,667.557 km²; the mesh's exact total is 11,667.557 km².

| class | mesh km² | ref (centres) km² | band [interior, touching] km² | mesh − ref ha | mesh − ref % | outside band | pass |
|---:|---:|---:|---|---:|---:|---:|---|
| 15 | 4198.5320 | 4198.4871 | [4194.6977, 4202.4629] | +4.49 | +0.0011 | 0 | yes |
| 3 | 2823.2657 | 2823.2545 | [2821.4419, 2825.1149] | +1.12 | +0.0004 | 0 | yes |
| 4 | 1191.8110 | 1191.8006 | [1190.5630, 1192.9908] | +1.05 | +0.0009 | 0 | yes |
| 12 | 964.3825 | 964.3493 | [960.9197, 967.8130] | +3.32 | +0.0034 | 0 | yes |
| 24 | 794.7068 | 794.6999 | [793.9541, 795.4425] | +0.68 | +0.0009 | 0 | yes |
| 21 | 778.4885 | 778.4923 | [777.7287, 779.2380] | −0.38 | −0.0005 | 0 | yes |
| 9 | 478.7770 | 478.7854 | [477.4191, 480.1170] | −0.84 | −0.0018 | 0 | yes |
| 30 | 164.2116 | 164.2281 | [163.6365, 164.7627] | −1.66 | −0.0101 | 0 | yes |
| 41 | 85.4943 | 85.4957 | [85.4516, 85.5364] | −0.14 | −0.0016 | 0 | yes |
| 29 | 78.0782 | 78.0894 | [77.6505, 78.4962] | −1.12 | −0.0143 | 0 | yes |
| 33 | 60.8048 | 60.8051 | [60.7831, 60.8262] | −0.04 | −0.0005 | 0 | yes |
| 11 | 20.7332 | 20.7327 | [20.7285, 20.7370] | +0.05 | +0.0024 | 0 | yes |
| 25 | 14.1565 | 14.1578 | [14.1130, 14.1993] | −0.13 | −0.0092 | 0 | yes |
| 39 | 7.7806 | 7.7816 | [7.7757, 7.7858] | −0.10 | −0.0129 | 0 | yes |
| 46, 75, 48, 47, 31 | 3.8421, 1.9019, 0.5246, 0.0530, 0.0127 | the same | [x, x] (no boundary cells) | 0.00 | 0 | 0 | yes |

Every class is inside its band, so the pass is met with the whole tolerance
to spare; the worst difference from the centre count is class 15 (+4.49 ha)
in area and class 29 (−0.0143 %) in share, both far inside their bands (the
bands with boundary cells are 0.84 to 777 ha wide). A second sum,
independent of the index-space construction (the file's float32 fractions
times each triangle's own ellipsoidal area from pyproj on its corners),
agrees: per triangle the two areas differ by at most 8.4e-5 relative (a UTM
straight-edged triangle against the sum of its cell pieces; the
`triangle_area_vs_geod` entry of the JSON files), and the totals are
11,667,557,060 m² (exact), 11,667,557,063 m² (Geod) and 11,667,557,062 m²
(the outline), within 3 m² of one another.

### Per-cell sums on a 1 km square against GEOS

A 1 km square around the mesh's centroid (35 × 36 cells, 1,260, at 10 m;
35 × 37 cells, 1,295, at 5 m, the centroids differ slightly; 62 and 103
triangles meet it). For every cell,
`_core.class_areas` with that cell alone as class 1 over those triangles,
against GEOS's clip of the same triangles with the cell's box times the row's
cell area, and against the cell's own area (the square is interior, so every
cell is wholly covered):

| mesh | worst relative difference to GEOS | worst relative difference to the cell's area | limit |
|---|---:|---:|---:|
| 10 m | 2.7e-16 | 1.3e-16 | 1e-6 |
| 5 m | 2.7e-16 | 2.7e-16 | 1e-6 |

### Misplaced share on squares (bounded, not measured)

Exact fractions misplace nothing; what follows measures the binning, not the
pass: mesh class areas binned by triangle centroid against raster class areas
binned by cell centre, over UTM-aligned squares at least half covered.

| square | squares | 10 m: median / p95 / max | 5 m: median / p95 / max |
|---|---:|---|---|
| 25 km | 18 | 0.06 % / 0.12 % / 0.15 % | 0.03 % / 0.06 % / 0.08 % |
| 5 km | 470 | 0.39 % / 1.22 % / 4.07 % | 0.23 % / 0.66 % / 1.40 % |
| 1 km | 11,655 | 3.2 % / 12.1 % / 100 % | 1.9 % / 6.6 % / 100 % |

The numbers halve from 10 to 5 m with the triangle size, as a binning
artefact does; the design's "5 km 95th percentile within 0.5 %" is a 26c
(cutoff) criterion and is not judged here.

### Time of the fractions pass

`--stats` phase `land cover fractions` (all 10 cores), median of 3, and the
sweep over `assign` (median of 3, the same inputs):

| mesh | triangles | CLI phase s | per M triangles | assign 1 thr | 2 | 4 | 8 | 10 | per M at 8 thr | speed-up 1 → 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 m | 649,504 | 0.493 | 0.76 s | 1.296 | 0.816 | 0.603 | 0.500 | 0.489 | **0.77 s** | 2.6× |
| 5 m | 1,571,955 | 0.984 | 0.63 s | 2.242 | 1.498 | 1.143 | 0.999 | 0.991 | **0.64 s** | 2.2× |

Under 2 s per million triangles at 8 threads (the brief's gate; section 5
of the design gives only the probe's one-thread cost, 0.15 µs per covered
cell plus 2.2 µs per triangle: the 10 m mesh at one thread took 1.30 s
against the probe's 2.2 µs × 649,504 = 1.43 s before the cells, so the
built code is not slower than the probe).

The parts, timed apart (`raw/analysis_*.txt`, "parts"): `_core.class_areas`
alone scales 7.4× (10 m: 0.878 s at 1 thread, 0.118 s at 8) and 7.4× (5 m:
1.444 → 0.195 s). The rest of `assign` is serial Python: the move to index
space 0.06 / 0.12 s, the class-window reads 0.17 / 0.36 s (3 and 6 chunks
of 2^18 triangles, each window nearly the whole raster, 7371 × 4494 cells,
because the triangles are in file order; the Hilbert order is 26c), and
0.15 / 0.33 s not timed separately (`row_areas`, a Python loop of one Geod
call per row per chunk; `fractions_of`; `dominant`; the concatenations).
At 8 threads the C++ is a quarter of the pass (10 m) or a fifth (5 m).

### Peak memory

`/usr/bin/time -l` maximum resident set size, median of 3:

| mesh | with land cover | without | added | per triangle |
|---|---:|---:|---:|---:|
| 10 m | 933.6 MB | 787.2 MB | +146 MB | 225 B |
| 5 m | 1,425.1 MB | 1,383.4 MB | +42 MB | 27 B |

The design's section 5 estimates about 25 bytes per triangle for the pass.
At 5 m the added peak matches; at 10 m it is 6 times that, and the two
disagree because the peak is measured over the whole process: at 10 m the
pass's working set (33 MB class window per chunk, the 650k × 3 × 2 doubles of
corners, 31 MB, pyproj's arrays) stands above refine's peak, at 5 m mostly
below it. Not a gate; recorded for 26c.

### Whole run

Median total of the `--stats` clock, with / without land cover: 10 m
1.494 / 0.949 s; 5 m 3.104 / 2.070 s. Wall clock (`time`): 1.99 / 1.43 s
and 3.83 / 2.82 s.

## Gate 2: no change to existing commands

`RASPUTIN_DATA=…/rasputin_data .venv/bin/python tools/bench_quick.py run --tree <tree>`
against the stored `docs/benchmarks/quick/baseline-ac.json` (commit
`357808b1`, refreshed 2026-10-11, identical on `master` and this branch).
Three runs, in order: the branch (`raw/quick_check.txt`), then the base
commit `5f28292c` and the branch again back to back
(`raw/quick_b2b_base.txt`, `raw/quick_b2b_branch.txt`).

**Meshes**: no `MESH CHANGED` line in any of the three runs, so every case's
mesh SHA-256 equals the baseline's (the script prints the line on any
difference).

**Timings**: every run's verdict is `SLOWER` (exit 1), the branch *and* the
untouched base commit alike, by 7 to 33 % in every case and in phases 26b
does not touch (`decode`, `features read`, `features clip`, `refine`), under
the Teams load measured above. Base against branch, back to back:

| case | baseline s | base 5f28292c s | branch 4149b50d s | branch / base |
|---|---:|---:|---:|---:|
| tile | 0.42 | 0.47 | 0.47 | 1.01 |
| tile t=1 | 0.70 | 0.84 | 0.85 | 1.01 |
| quarter | 0.36 | < 0.37 (in band) | 0.41 | ≤ 1.11 |
| numedalslagen | 8.83 | 9.77 | 9.41 | 0.96 |
| lagan | 14.54 | 16.15 | 16.09 | 1.00 |
| geilo-al-ramp | 0.33 | 0.38 | 0.38 | 0.99 |
| romsdal-slope | 3.58 | 4.36 | 4.31 | 0.99 |

(The base run's `quarter` is not in its output, so it was inside the 5 %
band, under 0.373 s; the branch's 0.41 s is the one case not matched by the
base, 0.05 s on a 0.36 s case, with the first branch run at 0.42 s.) The
quick check against the stored baseline is therefore not a measurement of
26b today; measured back to back, the branch is within 1 % of its base on
every case but `quarter`, where it is at most 11 % slower on a 0.4 s case.
A `mesh` without `--land-cover` runs the same code as before (`cli.py`'s
26b call is behind `if classes is not None`), which the equal meshes and
the equal times agree with. **Before the merge, re-run `tools/bench_quick.py run` on
an idle machine**: a three-minute run, and the formal verdict this gate
asks for.

**HOTSPOT lines** (`hotspots.toml` holds only Lagan `features clip` at
40.7 %): the branch's runs print `numedalslagen features clip 67 %` and
`features clip: clean-up 66 %` (expected: ruled on another unmerged branch
at 68 %), `tile refine 46 %`, `tile t=1 refine 71 %`, `tile t=1 refine:
scan (parallel) 53 %`, `quarter refine 48 %`, `geilo-al-ramp decode 43 %`,
`romsdal-slope refine 58 %`. The base commit prints the same set, and so do
the 2026-10-11 quick checks under `docs/benchmarks/2026-10-11/*/`
(`grep -rh HOTSPOT docs/benchmarks/2026-10-11`); none is new. Lagan
`features clip` is at 6.95 / 16.09 = 43 %, under the 50.7 % that would
raise it again.

## Verdicts

| gate | verdict |
|---|---|
| 1a class areas within the reference band + max(0.1 %, 1 ha) | **ACCEPTED**: all 19 classes inside the band itself, 0 outside, both meshes |
| 1b fractions pass under 2 s per million triangles at 8 threads | **ACCEPTED**: 0.77 s (10 m), 0.64 s (5 m) |
| 1c per-cell sums within 1e-6 relative on a 1 km square against shapely | **ACCEPTED**: 2.7e-16 worst, both meshes |
| 2 existing `mesh` unchanged | meshes equal in all three runs; timings: **no change attributable to 26b** back to back against `5f28292c` (within 1 % on six cases, `quarter` ≤ 11 % on 0.4 s); the stored-baseline verdict is `SLOWER` for base and branch alike under a measured foreign load, so it is to be **re-run idle before the merge** |

## Files

- `run_meshes.sh`: the twelve mesh runs. `analyse.py`: the reference, the
  sweep, the per-cell check, the misplaced share.
- `raw/velhas_t{10,5}_lc{on,off}_r{1,2,3}.{log,stats.md}`: per run.
- `raw/analysis_t{10,5}.{txt,json}`: the analysis output.
- `raw/quick_check.txt`, `raw/quick_b2b_{base,branch}.txt`: the quick checks;
  `raw/*pmset*`: power and load.
