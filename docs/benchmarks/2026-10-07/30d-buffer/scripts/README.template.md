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

{{MASTER}}

## Where the time goes: every buffer call (one traced run each)

{{TRACE}}

## Each call in isolation

{{ISOLATED}}

## The candidates in full runs (Lagan, Ljungan), seconds

{{CANDIDATES}}

## Pieces against GEOS's buffer, the pipeline's distance

{{PIECES}}

## The features region: from the grown domain against from the hull (100 m)

{{REGION}}

## How buffer time scales on the staircase

On the staircase outline, time is near zero up to one step (10 m) and jumps
at 20 m. It does not grow steadily with distance: on Lagan, 43.8 m took
16.9 s, 62 m took 11.6 s and 100 m took 55 s. On open pieces it grows faster
than linearly with the vertex count: on Lagan about 0.07 s for 1,000
vertices, 2.0 s for 8,000 and 81 s for the whole ring as one open line. That
is why short pieces are fast. With a smooth outline (Numedalslågen) every
call stays under 20 ms.

{{SCALING}}
