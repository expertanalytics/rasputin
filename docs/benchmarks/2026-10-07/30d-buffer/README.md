# 30d, step 1: the outline-buffer baseline on master `f81b20b7`

`@perf`, 2026-10-07. No production code changed. The tables are written by
`scripts/make_readme.py` from `raw/`; edit `scripts/README.template.md`, not
this file.

## Result

**The one hotspot is GEOS `buffer` of the SMHI staircase outline.** On Lagan
it is 98 % of `decode` and 95–97 % of `features read` (the traced runs). The worst call is
`feature_input.py:161`, which grows the whole domain by 100 m with round
corners: 28 s on Lagan. The GeoPackage path calls it **twice** (once for the
R-tree box, `feature_input.py:294`, and once for the clip region,
`feature_input.py:311`). The GeoJSON path calls it once.

**The best single fix is to grow the features region from the convex hull**
(candidate 2). With the shared buffer (candidate 1) in place, it takes Lagan
with the GeoPackage from 61.6 s to 31.1 s, and `features read` from 30.6 s to
1.1 s. All three candidates together take Lagan from 103.4 s to 16.4–18.1 s
with the GeoPackage, and from 72.7 s to 12.9–13.6 s with the GeoJSON window.
Every `.vtk` was byte-identical to master's, and so was every quality figure.

**Why the GeoPackage path is slower (103.4 s against 72.7 s on Lagan):**
1. The second 100 m buffer adds 28 s to `features read`. The traces show two
   calls, at 28.1 s each.
2. `features clip` is 3.7 s slower. The European GeoPackage hands over whole
   CORINE polygons: 15.9 M vertices on Lagan, against 1.74 M in the GeoJSON
   window (Ljungan: 14.9 M against 0.68 M). That the extra clip time comes
   from those extra vertices is likely but was not measured.
3. Reading costs about the same. The R-tree query takes 1.0 s warm (1.8 s on
   the first query of a process), against 1.5 s to parse and build the
   GeoJSON.
At 1.0 s for an 8.9 GB file, the query does not read the whole file.

**Confirmed without cProfile.** The audit (`docs/increments/perf-audit.md` on
branch `worktree-perf-audit`, PR #208, not on master) measured 35.7 s
`decode` and 30.1 s `features read` with the GeoJSON window. Here, with no
profiler, the median of 3 is 35.9 s and 29.6 s. Ola's afternoon GeoPackage
run (62.7 s `features read`, 112 s in total; one run) is close to the 57.2 s
and 103.4 s measured here.

**Caveat on growing the outline in pieces (candidate 3).** On the staircase
outlines, every piece size gives GEOS's region to within 6.1e-7 m. On
Numedalslågen's smooth NVE outline it does not. There the pieced region
reaches up to 10.8 m beyond GEOS's at one place where the outline
turns back on itself (near vertex 7541), and GEOS's pokes out of the pieced region by 6.9e-6 m, which
fails a 1e-6 m check. Why GEOS's region is smaller there was not measured.
The candidate as patched applies to any mitred polygon of 5,000 vertices or
more, so it would reach Numedalslågen's `dem_input.py:250` buffer. That buffer
takes 5 ms there. No patched mesh run was made on Numedalslågen.

**What was not run.** Ola cut the run: "the question needed minutes, not an
hour". The patched full runs on Numedalslågen were therefore skipped. Its
master runs and probes had already finished, and they serve as the second
control.

## Method

- **Code:** master `f81b20b7`'s `src_python/tin_engine` (`git archive`) in a
  scratch directory, with a `_core` built from the same commit: Release,
  `-O3 -DNDEBUG`, `RASPUTIN_HARDENING=ON` (bounds checks on, libc++ fast).
  The build's SHA-256 is `86765ff1…`, equal to the main checkout's `.venv`
  `_core`. Every run prints `tin_engine.__file__` and `_core.__file__`
  (`raw/runs/*_stderr.txt`).
- **Machine:** Apple M1 Max, 10 cores, 32 GiB, macOS 27.0, Python 3.14.7,
  shapely 2.1.2, GEOS 3.13.1. Threads: 10.
- **Power:** AC for every run. `pmset -g batt` was read before and after each
  run (`raw/**/*_power.txt`), and the generator refuses any run that was not
  on AC both times. Nothing else ran.
- **Inputs:** `--tolerance 10 --binary`.
  - Lagan and Ljungan above Flåsjö: as in Ola's
    `rasputin_scratch/sweden/mesh_sweden.sh` (GLO-30 cache, `--out-crs
    EPSG:3006`, SMHI SVAR outline). Features were either the CORINE GeoJSON
    window `rasputin_data/sweden_corine/<c>_clc2018_3035.geojson` or the
    European GeoPackage `U2018_CLC2018_V2020_20u1.gpkg`.
  - Numedalslågen: DTM10 with the NVE outline. Features were either
    `corine2018_dtm10_utm33.gpkg` (as in the 2026-10-06 profiles) or the
    European GeoPackage.
- **Timed runs:** `scripts/timed_master.sh` and `scripts/timed_candidates.sh`
  (each through `run_mesh.sh` and `launch.py`). Each case had one warm-up
  (`w0`, discarded) and three repeats. Each cell shows the median in bold and
  the three values. The phases are `--stats`'s own rows. The candidate sweep
  was stopped once and resumed; the script skips finished runs.
- **Traces:** one extra run per case with `launch.py --trace`. This is a
  Python wrapper around each `buffer`, `source_region`, `query_features` and
  `read_json` call, not a profiler (`raw/trace/`).
- **Isolated calls:** `scripts/isolated.py` and `scripts/candidates.py`
  (`raw/probes/`), median of 3. Distances between the pieced region and
  GEOS's are sampled every 0.25 m along both boundaries.
- **Candidates:** as in the audit's `run_patched.py`, now in `launch.py` and
  `launch_pieces.py`. (1) Each buffer of the same polygon with the same
  arguments is computed once. (2) `source_region` is grown from the convex
  hull. (3) A mitred polygon of 5,000 vertices or more is grown in pieces of
  K edges and united with the polygon.
- **Meshes** are not committed. They were written to the session scratchpad
  and removed; `run_mesh.sh` regenerates them. Their SHA-256 values are in
  `raw/runs/vtk_sha256.txt`. The GeoJSON window and the GeoPackage give
  slightly different Lagan meshes (867,612 against 867,678 triangles) with
  the same features kept; the cause was not measured.

## Master: phase times, seconds

| catchment | features | decode | features read | features clip | total |
|---|---|---|---|---|---|
| Lagan | GeoJSON window | **35.88** (35.09, 35.88, 35.88) | **29.63** (29.60, 29.74, 29.63) | **3.11** (3.10, 3.11, 3.12) | **72.74** (71.90, 72.86, 72.74) |
| Lagan | European GeoPackage | **35.12** (35.12, 36.05, 35.08) | **57.23** (57.37, 56.98, 57.23) | **6.79** (6.79, 6.78, 6.84) | **103.38** (103.38, 103.91, 103.35) |
| Ljungan above Flåsjö | GeoJSON window | **4.94** (4.75, 4.94, 5.03) | **2.90** (2.85, 2.92, 2.90) | **0.58** (0.56, 0.58, 0.58) | **9.29** (8.93, 9.29, 9.34) |
| Ljungan above Flåsjö | European GeoPackage | **4.89** (4.71, 4.89, 5.10) | **5.56** (5.56, 5.50, 5.66) | **4.36** (4.29, 4.36, 4.42) | **15.57** (15.39, 15.57, 16.01) |
| Numedalslågen | UTM33 GeoPackage | **1.17** (1.18, 1.17, 1.03) | **0.40** (0.40, 0.41, 0.26) | **1.81** (1.81, 1.81, 1.83) | **6.97** (6.97, 6.99, 6.74) |
| Numedalslågen | European GeoPackage | **1.03** (1.02, 1.17, 1.03) | **1.24** (1.25, 1.16, 1.24) | **6.23** (6.22, 6.23, 6.25) | **12.15** (12.15, 12.20, 12.14) |

## Where the time goes: every buffer call (one traced run each)

| catchment | features | call site | calls | distance | join | vertices in | seconds |
|---|---|---|---|---|---|---|---|
| Lagan | GeoJSON window | `feature_input.py:161` | 1 | 100 | round | 56261 | 27.964 |
| Lagan | GeoJSON window | `dem_input.py:213` | 1 | 43.84 | mitre | 56261 | 17.555 |
| Lagan | GeoJSON window | `target_grid.py:107` | 1 | 43.84 | mitre | 56261 | 17.514 |
| Lagan | GeoJSON window | `read_json` (`feature_input.py:419`) | 1 |  |  |  | 0.715 |
| Lagan | GeoJSON window | `target_grid.py:145` | 1 | 0.0005008 | mitre | 40176 | 0.203 |
| Lagan | GeoJSON window | `target_grid.py:243` | 90 | 31 | round | 715/760/980… | 0.035 |
| Lagan | GeoJSON window | *decode* in this run: 36.024 s, of which buffers 35.307 s | | | | | rest 0.717 |
| Lagan | GeoJSON window | *features read* in this run: 29.498 s, of which buffers 27.964 s | | | | | rest 1.534 |
| Lagan | European GeoPackage | `feature_input.py:161` | 2 | 100 | round | 56261 | 56.200 |
| Lagan | European GeoPackage | `target_grid.py:107` | 1 | 43.84 | mitre | 56261 | 17.263 |
| Lagan | European GeoPackage | `dem_input.py:213` | 1 | 43.84 | mitre | 56261 | 16.729 |
| Lagan | European GeoPackage | `query_features` (`feature_input.py:407`) | 1 |  |  | 7866 rows | 1.779 |
| Lagan | European GeoPackage | `target_grid.py:145` | 1 | 0.0005008 | mitre | 40176 | 0.202 |
| Lagan | European GeoPackage | `target_grid.py:243` | 90 | 31 | round | 715/760/980… | 0.048 |
| Lagan | European GeoPackage | *decode* in this run: 34.907 s, of which buffers 34.243 s | | | | | rest 0.664 |
| Lagan | European GeoPackage | *features read* in this run: 58.079 s, of which buffers 56.200 s | | | | | rest 1.879 |
| Ljungan above Flåsjö | GeoJSON window | `feature_input.py:161` | 1 | 100 | round | 23260 | 2.211 |
| Ljungan above Flåsjö | GeoJSON window | `target_grid.py:107` | 1 | 43.84 | mitre | 23260 | 2.090 |
| Ljungan above Flåsjö | GeoJSON window | `dem_input.py:213` | 1 | 43.84 | mitre | 23260 | 2.071 |
| Ljungan above Flåsjö | GeoJSON window | `read_json` (`feature_input.py:419`) | 1 |  |  |  | 0.242 |
| Ljungan above Flåsjö | GeoJSON window | `target_grid.py:145` | 1 | 0.0005008 | mitre | 18451 | 0.082 |
| Ljungan above Flåsjö | GeoJSON window | `target_grid.py:243` | 15 | 31 | round | 710/812/923… | 0.007 |
| Ljungan above Flåsjö | GeoJSON window | *decode* in this run: 4.569 s, of which buffers 4.250 s | | | | | rest 0.319 |
| Ljungan above Flåsjö | GeoJSON window | *features read* in this run: 2.757 s, of which buffers 2.211 s | | | | | rest 0.546 |
| Ljungan above Flåsjö | European GeoPackage | `feature_input.py:161` | 2 | 100 | round | 23260 | 4.372 |
| Ljungan above Flåsjö | European GeoPackage | `target_grid.py:107` | 1 | 43.84 | mitre | 23260 | 2.103 |
| Ljungan above Flåsjö | European GeoPackage | `dem_input.py:213` | 1 | 43.84 | mitre | 23260 | 2.095 |
| Ljungan above Flåsjö | European GeoPackage | `query_features` (`feature_input.py:407`) | 1 |  |  | 1156 rows | 1.020 |
| Ljungan above Flåsjö | European GeoPackage | `target_grid.py:145` | 1 | 0.0005008 | mitre | 18451 | 0.082 |
| Ljungan above Flåsjö | European GeoPackage | `target_grid.py:243` | 15 | 31 | round | 710/812/923… | 0.006 |
| Ljungan above Flåsjö | European GeoPackage | *decode* in this run: 4.601 s, of which buffers 4.287 s | | | | | rest 0.314 |
| Ljungan above Flåsjö | European GeoPackage | *features read* in this run: 5.465 s, of which buffers 4.372 s | | | | | rest 1.093 |
| Numedalslågen | UTM33 GeoPackage | `query_features` (`feature_input.py:407`) | 1 |  |  | 7932 rows | 0.194 |
| Numedalslågen | UTM33 GeoPackage | `feature_input.py:161` | 2 | 100 | round | 14093 | 0.036 |
| Numedalslågen | UTM33 GeoPackage | `dem_input.py:250` | 1 | 14.14 | mitre | 14093 | 0.007 |
| Numedalslågen | UTM33 GeoPackage | *decode* in this run: 1.137 s, of which buffers 0.007 s | | | | | rest 1.130 |
| Numedalslågen | UTM33 GeoPackage | *features read* in this run: 0.291 s, of which buffers 0.036 s | | | | | rest 0.255 |
| Numedalslågen | European GeoPackage | `query_features` (`feature_input.py:407`) | 1 |  |  | 8115 rows | 1.170 |
| Numedalslågen | European GeoPackage | `feature_input.py:161` | 2 | 100 | round | 14093 | 0.036 |
| Numedalslågen | European GeoPackage | `dem_input.py:250` | 1 | 14.14 | mitre | 14093 | 0.005 |
| Numedalslågen | European GeoPackage | *decode* in this run: 1.039 s, of which buffers 0.005 s | | | | | rest 1.034 |
| Numedalslågen | European GeoPackage | *features read* in this run: 1.287 s, of which buffers 0.036 s | | | | | rest 1.251 |

## Each call in isolation

| catchment | step | seconds, median (3 runs) |
|---|---|---|
| Lagan | footprints (tile headers) | **0.000** (0.039, 0.000, 0.000) |
| Lagan | domain.to_crs(target) | **0.000** (0.002, 0.000, 0.000) |
| Lagan | target_grid_for (target_grid.py:107's buffer inside) | **17.341** (17.489, 17.341, 16.889) |
| Lagan | buffer target_grid.py:107 = dem_input.py:213 (the same call) (d = 43.84, 56261 vertices) | **16.860** (16.860, 16.961, 16.858) |
| Lagan | target_grid.source_region (target_grid.py:145's buffer inside) | **0.211** (0.212, 0.211, 0.210) |
| Lagan | buffer target_grid.py:145 (d = 0.0005008, 40176 vertices) | **0.201** (0.201, 0.201, 0.201) |
| Lagan | plan_mosaic | **0.002** (0.003, 0.002, 0.002) |
| Lagan | repository.check | **0.000** (0.001, 0.000, 0.000) |
| Lagan | assemble (decode + seams) | **0.136** (0.152, 0.136, 0.135) |
| Lagan | resample | **0.498** (0.510, 0.494, 0.498) |
| Lagan | buffer feature_input.py:161 (d = 100, 56261 vertices) | **28.152** (28.196, 28.152, 28.097) |
| Lagan | feature_input.source_region(domain, dem, dem) | **28.131** (28.131, 28.007, 28.187) |
| Lagan | read_json (GeoJSON window) | **0.741** (0.693, 0.741, 0.742) |
| Lagan | read_source (GeoJSON window: read_json + shape) | **1.569** (1.654, 1.569, 1.557) |
| Lagan | feature_input.source_region(domain, dem, EPSG:3035) [gpkg box] | **27.997** (28.051, 27.937, 27.997) |
| Lagan | query_features (gpkg, 8.86 GB) | **1.024** (1.082, 1.024, 0.992) |
| Ljungan above Flåsjö | footprints (tile headers) | **0.000** (0.040, 0.000, 0.000) |
| Ljungan above Flåsjö | domain.to_crs(target) | **0.000** (0.002, 0.000, 0.000) |
| Ljungan above Flåsjö | target_grid_for (target_grid.py:107's buffer inside) | **2.141** (2.196, 2.141, 2.115) |
| Ljungan above Flåsjö | buffer target_grid.py:107 = dem_input.py:213 (the same call) (d = 43.84, 23260 vertices) | **2.087** (2.086, 2.087, 2.205) |
| Ljungan above Flåsjö | target_grid.source_region (target_grid.py:145's buffer inside) | **0.090** (0.090, 0.090, 0.089) |
| Ljungan above Flåsjö | buffer target_grid.py:145 (d = 0.0005008, 18451 vertices) | **0.084** (0.084, 0.084, 0.087) |
| Ljungan above Flåsjö | plan_mosaic | **0.000** (0.001, 0.000, 0.000) |
| Ljungan above Flåsjö | repository.check | **0.000** (0.000, 0.000, 0.000) |
| Ljungan above Flåsjö | assemble (decode + seams) | **0.074** (0.079, 0.074, 0.074) |
| Ljungan above Flåsjö | resample | **0.194** (0.211, 0.194, 0.193) |
| Ljungan above Flåsjö | buffer feature_input.py:161 (d = 100, 23260 vertices) | **2.259** (2.313, 2.226, 2.259) |
| Ljungan above Flåsjö | feature_input.source_region(domain, dem, dem) | **2.262** (2.313, 2.262, 2.240) |
| Ljungan above Flåsjö | read_json (GeoJSON window) | **0.256** (0.240, 0.258, 0.256) |
| Ljungan above Flåsjö | read_source (GeoJSON window: read_json + shape) | **0.547** (0.582, 0.547, 0.547) |
| Ljungan above Flåsjö | feature_input.source_region(domain, dem, EPSG:3035) [gpkg box] | **2.229** (2.229, 2.220, 2.272) |
| Ljungan above Flåsjö | query_features (gpkg, 8.86 GB) | **1.009** (1.727, 1.009, 0.946) |
| Numedalslågen | footprints (tile headers) | **0.000** (0.242, 0.000, 0.000) |
| Numedalslågen | _domain_plan (dem_input.py:250's buffer inside) | **0.110** (0.129, 0.110, 0.106) |
| Numedalslågen | buffer dem_input.py:250 (d = 14.14, 14093 vertices) | **0.005** (0.005, 0.005, 0.005) |
| Numedalslågen | repository.check | **0.000** (0.000, 0.000, 0.000) |
| Numedalslågen | assemble (decode + seams) | **0.888** (1.024, 0.888, 0.685) |
| Numedalslågen | buffer feature_input.py:161 (d = 100, 14093 vertices) | **0.017** (0.018, 0.017, 0.017) |
| Numedalslågen | feature_input.source_region(domain, dem, dem) | **0.035** (0.037, 0.035, 0.035) |
| Numedalslågen | feature_input.source_region(domain, dem, EPSG:25833) [gpkg33 box] | **0.035** (0.035, 0.035, 0.035) |
| Numedalslågen | query_features (gpkg33, 0.44 GB) | **0.235** (0.320, 0.235, 0.177) |
| Numedalslågen | feature_input.source_region(domain, dem, EPSG:3035) [gpkg box] | **0.040** (0.043, 0.040, 0.040) |
| Numedalslågen | query_features (gpkg, 8.86 GB) | **1.077** (1.888, 1.077, 1.076) |

What each input hands to the clip:

| catchment | features | rows read | vertices read |
|---|---|---|---|
| Lagan | GeoJSON window | 7319 | 1,744,538 |
| Lagan | European GeoPackage | 7866 | 15,888,613 |
| Ljungan above Flåsjö | GeoJSON window | 851 | 678,103 |
| Ljungan above Flåsjö | European GeoPackage | 1156 | 14,921,039 |
| Numedalslågen | UTM33 GeoPackage | 7932 | 5,417,093 |
| Numedalslågen | European GeoPackage | 8115 | 16,808,734 |

## The candidates in full runs (Lagan, Ljungan), seconds

| catchment | features | variant | decode | features read | features clip | total | `.vtk` vs master | worst angle | max degree | largest error m | DEM nodes outside |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Lagan | GeoJSON window | master | **35.88** (35.09, 35.88, 35.88) | **29.63** (29.60, 29.74, 29.63) | **3.11** (3.10, 3.11, 3.12) | **72.74** (71.90, 72.86, 72.74) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | GeoJSON window | 1: buffer once | **18.81** (18.81, 19.32, 18.13) | **29.98** (29.98, 30.20, 29.81) | **3.14** (3.13, 3.14, 3.14) | **56.18** (56.18, 56.93, 55.41) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | GeoJSON window | 1+2: + region from the hull | **18.88** (18.88, 18.95, 18.77) | **1.52** (1.52, 1.53, 1.50) | **3.15** (3.15, 3.17, 3.12) | **27.82** (27.82, 27.99, 27.68) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | GeoJSON window | 1+2+3, pieces of 1000 | **4.69** (4.69, 4.68, 4.72) | **1.50** (1.50, 1.50, 1.52) | **3.14** (3.19, 3.10, 3.14) | **13.64** (13.64, 13.57, 13.67) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | GeoJSON window | 1+2+3, pieces of 500 | **4.19** (4.18, 4.21, 4.19) | **1.51** (1.51, 1.51, 1.50) | **3.17** (3.15, 3.17, 3.18) | **13.24** (13.11, 13.24, 13.33) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | GeoJSON window | 1+2+3, pieces of 250 | **3.86** (3.83, 3.86, 3.91) | **1.51** (1.51, 1.51, 1.54) | **3.17** (3.16, 3.17, 3.35) | **12.89** (12.82, 12.89, 13.17) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | European GeoPackage | master | **35.12** (35.12, 36.05, 35.08) | **57.23** (57.37, 56.98, 57.23) | **6.79** (6.79, 6.78, 6.84) | **103.38** (103.38, 103.91, 103.35) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | European GeoPackage | 1: buffer once | **19.52** (20.15, 19.52, 18.77) | **30.62** (30.62, 30.80, 29.74) | **6.89** (7.25, 6.89, 6.89) | **61.62** (62.43, 61.62, 59.82) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | European GeoPackage | 1+2: + region from the hull | **18.73** (18.73, 19.62, 18.73) | **1.08** (1.10, 1.08, 1.07) | **6.96** (7.01, 6.96, 6.80) | **31.11** (31.11, 32.26, 30.84) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | European GeoPackage | 1+2+3, pieces of 1000 | **4.80** (4.89, 4.80, 4.64) | **1.85** (1.95, 1.85, 1.05) | **6.95** (7.01, 6.93, 6.95) | **18.06** (18.63, 18.06, 17.00) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | European GeoPackage | 1+2+3, pieces of 500 | **4.32** (4.30, 4.52, 4.32) | **1.90** (1.90, 1.95, 1.32) | **6.85** (6.85, 7.37, 6.82) | **17.49** (17.49, 18.36, 16.78) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Lagan | European GeoPackage | 1+2+3, pieces of 250 | **3.92** (4.03, 3.92, 3.87) | **1.14** (1.79, 1.14, 1.07) | **6.97** (7.13, 6.97, 6.90) | **16.44** (17.39, 16.44, 16.23) | identical | 0.00238° | 17 | 9.99975411631047 | 0 |
| Ljungan above Flåsjö | GeoJSON window | master | **4.94** (4.75, 4.94, 5.03) | **2.90** (2.85, 2.92, 2.90) | **0.58** (0.56, 0.58, 0.58) | **9.29** (8.93, 9.29, 9.34) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | GeoJSON window | 1: buffer once | **2.58** (2.58, 2.58, 2.60) | **2.88** (2.88, 2.88, 2.84) | **0.57** (0.57, 0.57, 0.56) | **6.85** (6.85, 6.89, 6.80) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | GeoJSON window | 1+2: + region from the hull | **2.55** (2.55, 2.66, 2.54) | **0.55** (0.55, 0.55, 0.55) | **0.56** (0.56, 0.56, 0.56) | **4.45** (4.45, 4.55, 4.42) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | GeoJSON window | 1+2+3, pieces of 1000 | **2.11** (2.11, 2.23, 2.11) | **0.55** (0.55, 0.56, 0.55) | **0.56** (0.57, 0.56, 0.56) | **4.02** (4.02, 4.14, 4.02) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | GeoJSON window | 1+2+3, pieces of 500 | **1.87** (1.87, 1.87, 1.91) | **0.57** (0.57, 0.56, 0.57) | **0.57** (0.57, 0.56, 0.57) | **3.81** (3.81, 3.79, 3.87) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | GeoJSON window | 1+2+3, pieces of 250 | **1.85** (1.83, 1.86, 1.85) | **0.55** (0.56, 0.55, 0.55) | **0.56** (0.56, 0.56, 0.57) | **3.76** (3.75, 3.76, 3.77) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | European GeoPackage | master | **4.89** (4.71, 4.89, 5.10) | **5.56** (5.56, 5.50, 5.66) | **4.36** (4.29, 4.36, 4.42) | **15.57** (15.39, 15.57, 16.01) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | European GeoPackage | 1: buffer once | **2.58** (2.56, 2.69, 2.58) | **3.28** (3.26, 3.29, 3.28) | **4.27** (4.36, 4.27, 4.27) | **11.00** (11.00, 11.05, 10.94) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | European GeoPackage | 1+2: + region from the hull | **2.69** (2.69, 2.72, 2.58) | **1.01** (1.01, 1.06, 0.98) | **4.28** (4.28, 4.33, 4.27) | **8.77** (8.77, 8.93, 8.62) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | European GeoPackage | 1+2+3, pieces of 1000 | **2.12** (2.12, 2.10, 2.13) | **1.00** (1.00, 0.99, 1.02) | **4.26** (4.26, 4.26, 4.31) | **8.20** (8.20, 8.17, 8.26) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | European GeoPackage | 1+2+3, pieces of 500 | **1.88** (1.88, 1.87, 1.91) | **1.01** (1.01, 0.97, 1.01) | **4.32** (4.32, 4.27, 4.33) | **8.04** (8.04, 7.92, 8.06) | identical | 0.024° | 20 | 9.999125803873994 | 0 |
| Ljungan above Flåsjö | European GeoPackage | 1+2+3, pieces of 250 | **1.85** (1.86, 1.84, 1.85) | **0.99** (0.98, 0.99, 1.04) | **4.33** (4.33, 4.27, 4.33) | **7.99** (7.99, 7.90, 8.01) | identical | 0.024° | 20 | 9.999125803873994 | 0 |

## Pieces against GEOS's buffer, the pipeline's distance

| catchment | method | seconds, median (3 runs) | sym. difference m² | largest boundary distance m | pokes out of GEOS's m | GEOS's pokes out m | bounds, largest difference m | every vertex within 1e-6 m |
|---|---|---|---|---|---|---|---|---|
| Lagan (56260 edges, d = 43.84 m) | GEOS buffer, mitre | **16.989** (17.163, 16.989, 16.925) | | | | | | |
| Lagan | pieces of 250 | **2.864** (2.867, 2.859, 2.864) | 0.0205 | 6.13e-07 | 6.13e-07 | 0 | 4.34e-07 | yes |
| Lagan | pieces of 500 | **3.201** (3.201, 3.196, 3.209) | 0 | 0 | 0 | 0 | 0 | yes |
| Lagan | pieces of 1000 | **3.696** (3.704, 3.671, 3.696) | 0 | 0 | 0 | 0 | 0 | yes |
| Ljungan above Flåsjö (23259 edges, d = 43.84 m) | GEOS buffer, mitre | **2.188** (2.202, 2.169, 2.188) | | | | | | |
| Ljungan above Flåsjö | pieces of 250 | **1.203** (1.211, 1.191, 1.203) | 0.00441 | 6.13e-07 | 6.13e-07 | 0 | 0 | yes |
| Ljungan above Flåsjö | pieces of 500 | **1.261** (1.261, 1.253, 1.261) | 0.0029 | 6.13e-07 | 6.13e-07 | 0 | 0 | yes |
| Ljungan above Flåsjö | pieces of 1000 | **1.513** (1.509, 1.513, 1.532) | 0 | 0 | 0 | 0 | 0 | yes |
| Numedalslågen (14092 edges, d = 14.14 m) | GEOS buffer, mitre | **0.005** (0.005, 0.005, 0.005) | | | | | | |
| Numedalslågen | pieces of 250 | **0.083** (0.084, 0.083, 0.083) | 101 | 10.8 | 10.8 | 6.92e-06 | 0 | no |
| Numedalslågen | pieces of 500 | **0.069** (0.069, 0.068, 0.069) | 101 | 10.8 | 10.8 | 6.92e-06 | 0 | no |
| Numedalslågen | pieces of 1000 | **0.057** (0.057, 0.058, 0.057) | 101 | 10.8 | 10.8 | 6.92e-06 | 0 | no |

## The features region: from the grown domain against from the hull (100 m)

| catchment | master: hull of the grown domain, s | candidate: the hull grown, s | sym. difference m² | largest boundary distance m | bounds, largest difference m | candidate covers the domain grown by 99.5 m |
|---|---|---|---|---|---|---|
| Lagan | **28.194** (28.194, 28.187, 28.245) | **0.004** (0.004, 0.004, 0.004) | 1.16e+04 | 0.215 | 0 | yes |
| Ljungan above Flåsjö | **2.267** (2.267, 2.287, 2.212) | **0.002** (0.002, 0.002, 0.002) | 5741 | 0.252 | 0 | yes |
| Numedalslågen | **0.037** (0.037, 0.036, 0.037) | **0.001** (0.001, 0.001, 0.001) | 1.302e+04 | 0.232 | 0.0423 | yes |

## How buffer time scales on the staircase

On the staircase outline, time is near zero up to one step (10 m) and jumps
at 20 m. It does not grow steadily with distance: on Lagan, 43.8 m took
16.9 s, 62 m took 11.6 s and 100 m took 55 s. On open pieces it grows faster
than linearly with the vertex count: on Lagan about 0.07 s for 1,000
vertices, 2.0 s for 8,000 and 81 s for the whole ring as one open line. That
is why short pieces are fast. With a smooth outline (Numedalslågen) every
call stays under 20 ms.

| catchment | shape | join | distance m | vertices | seconds, median (3 runs) |
|---|---|---|---|---|---|
| Lagan | polygon | mitre | 5 | 56260 | **0.029** (0.029, 0.028, 0.029) |
| Lagan | polygon | mitre | 10 | 56260 | **0.032** (0.033, 0.032, 0.032) |
| Lagan | polygon | mitre | 20 | 56260 | **1.807** (1.831, 1.807, 1.799) |
| Lagan | polygon | mitre | 31 | 56260 | **2.473** (2.515, 2.473, 2.383) |
| Lagan | polygon | mitre | 43.84 | 56260 | **16.852** (16.954, 16.852, 16.836) |
| Lagan | polygon | mitre | 62 | 56260 | **11.572** (11.572, 12.646, 11.567) |
| Lagan | polygon | mitre | 100 | 56260 | **54.984** (55.929, 53.513, 54.984) |
| Lagan | polygon | round | 100 | 56260 | **27.974** (27.974, 27.963, 28.099) |
| Lagan | open piece | mitre | 43.84 | 251 | **0.013** (0.014, 0.013, 0.013) |
| Lagan | open piece | mitre | 43.84 | 501 | **0.029** (0.029, 0.029, 0.029) |
| Lagan | open piece | mitre | 43.84 | 1001 | **0.070** (0.070, 0.070, 0.070) |
| Lagan | open piece | mitre | 43.84 | 2001 | **0.133** (0.133, 0.132, 0.133) |
| Lagan | open piece | mitre | 43.84 | 4001 | **0.441** (0.440, 0.441, 0.444) |
| Lagan | open piece | mitre | 43.84 | 8001 | **2.024** (2.036, 2.024, 2.012) |
| Lagan | open piece | mitre | 43.84 | 16001 | **11.247** (11.226, 11.269, 11.247) |
| Lagan | open piece | mitre | 43.84 | 32001 | **39.133** (38.857, 39.133, 40.547) |
| Lagan | open piece | mitre | 43.84 | 56261 | **80.981** (80.445, 83.494, 80.981) |
| Ljungan above Flåsjö | polygon | mitre | 5 | 23259 | **0.011** (0.012, 0.011, 0.011) |
| Ljungan above Flåsjö | polygon | mitre | 10 | 23259 | **0.011** (0.011, 0.011, 0.011) |
| Ljungan above Flåsjö | polygon | mitre | 20 | 23259 | **0.482** (0.489, 0.482, 0.478) |
| Ljungan above Flåsjö | polygon | mitre | 31 | 23259 | **0.589** (0.641, 0.589, 0.586) |
| Ljungan above Flåsjö | polygon | mitre | 43.84 | 23259 | **2.156** (2.101, 2.294, 2.156) |
| Ljungan above Flåsjö | polygon | mitre | 62 | 23259 | **2.396** (2.396, 2.504, 2.260) |
| Ljungan above Flåsjö | polygon | mitre | 100 | 23259 | **9.644** (9.753, 9.285, 9.644) |
| Ljungan above Flåsjö | polygon | round | 100 | 23259 | **2.216** (2.216, 2.220, 2.213) |
| Ljungan above Flåsjö | open piece | mitre | 43.84 | 251 | **0.014** (0.015, 0.014, 0.014) |
| Ljungan above Flåsjö | open piece | mitre | 43.84 | 501 | **0.030** (0.030, 0.030, 0.030) |
| Ljungan above Flåsjö | open piece | mitre | 43.84 | 1001 | **0.074** (0.074, 0.074, 0.074) |
| Ljungan above Flåsjö | open piece | mitre | 43.84 | 2001 | **0.152** (0.172, 0.151, 0.152) |
| Ljungan above Flåsjö | open piece | mitre | 43.84 | 4001 | **0.342** (0.367, 0.339, 0.342) |
| Ljungan above Flåsjö | open piece | mitre | 43.84 | 8001 | **0.934** (0.920, 0.948, 0.934) |
| Ljungan above Flåsjö | open piece | mitre | 43.84 | 16001 | **4.382** (4.382, 4.572, 4.292) |
| Ljungan above Flåsjö | open piece | mitre | 43.84 | 23260 | **11.092** (11.523, 10.412, 11.092) |
| Numedalslågen | polygon | mitre | 5 | 14092 | **0.004** (0.004, 0.004, 0.004) |
| Numedalslågen | polygon | mitre | 10 | 14092 | **0.004** (0.004, 0.004, 0.004) |
| Numedalslågen | polygon | mitre | 20 | 14092 | **0.005** (0.006, 0.005, 0.005) |
| Numedalslågen | polygon | mitre | 31 | 14092 | **0.006** (0.006, 0.006, 0.006) |
| Numedalslågen | polygon | mitre | 14.14 | 14092 | **0.005** (0.005, 0.005, 0.005) |
| Numedalslågen | polygon | mitre | 62 | 14092 | **0.010** (0.010, 0.010, 0.010) |
| Numedalslågen | polygon | mitre | 100 | 14092 | **0.015** (0.015, 0.015, 0.015) |
| Numedalslågen | polygon | round | 100 | 14092 | **0.017** (0.017, 0.017, 0.017) |
| Numedalslågen | open piece | mitre | 14.14 | 251 | **0.000** (0.000, 0.000, 0.000) |
| Numedalslågen | open piece | mitre | 14.14 | 501 | **0.000** (0.000, 0.000, 0.000) |
| Numedalslågen | open piece | mitre | 14.14 | 1001 | **0.001** (0.001, 0.001, 0.001) |
| Numedalslågen | open piece | mitre | 14.14 | 2001 | **0.001** (0.001, 0.001, 0.001) |
| Numedalslågen | open piece | mitre | 14.14 | 4001 | **0.003** (0.003, 0.003, 0.003) |
| Numedalslågen | open piece | mitre | 14.14 | 8001 | **0.006** (0.006, 0.006, 0.006) |
| Numedalslågen | open piece | mitre | 14.14 | 14093 | **0.011** (0.011, 0.011, 0.011) |
