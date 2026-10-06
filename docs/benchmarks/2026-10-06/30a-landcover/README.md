# Increment 30a acceptance: the land cover phase of `rasputin mesh`

@perf, 2026-10-06. Base `a7154ec` (the profile commit the design starts from)
against branch `worktree-landcover-speed` at `3066d60`. The acceptance is
section 9 of `docs/increments/30a-landcover-speed.md`.

**Verdict: ACCEPTED.** The meshes are byte-identical, and the land cover
phase is 4.2 times faster on Numedalslågen and 3.8 times faster on
Skiensvassdraget. The gate was "at most a third of the base".

## Method

- **Machine:** Apple M1 Max (8 performance + 2 efficiency cores), 32 GB,
  macOS 27.0. Python 3.13.15, numpy 2.5.3, shapely 2.1.2 on both sides.
- **Power: AC for every run.** `pmset -g batt` was read before and after
  each run (`raw/stats/*_power.txt`, `raw/probe_power.txt`), and all 26
  readings say "AC Power". The profile in `../bottlenecks/` was made on
  battery, so its seconds are quoted below and not compared.
- **Installs:** the base is `a7154ec`'s non-editable install
  (`worktree-bottlenecks/.venv`). The branch is the worktree's own uv `.venv`,
  an editable install of `3066d60`. No C++ changes between them, and both
  load the same `_core` extension: the sha256 of the two `.so` files is the
  same (`5b0adf95…`). Both are Release builds with bounds checks on (libc++
  fast mode), as each stats file says. So the comparison is like with like
  on build mode and hardening too.
- **Command:** `scripts/stats.sh <side> <catchment> <repeat>`, run by
  `scripts/run_all.sh`. It runs `rasputin mesh` with DEM
  `rasputin_data/DTM10_UTM33_20260925`, domain
  `rasputin_scratch/norway/<c>/<c>_outline_nve.geojson`, `--features
  rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018
  --features-map corine --tolerance 10 --binary --stats`, with 10 threads.
  This is the bottleneck profile's command, and the profile ran it the same way.
- **Repeats:** 3 per catchment per side, with base and branch alternated:
  base first on odd repeats, branch first on even ones. The tables give the
  median. `scripts/summarize.py` turns `raw/stats/` into `raw/summary.txt`.
- **Byte-identity:** the sha256 of every `.vtk` written by the 12 timing runs
  (`raw/stats/vtk_sha256.txt`). Then the design's probe
  (`docs/increments/30a-probes/landcover_bytes.py`), run on the branch after
  the timing runs and never beside them, in `fixtures` mode
  (`raw/probe_fixtures_branch.txt`) and `mesh` mode
  (`raw/probe_mesh_branch.txt`), compared with
  `docs/increments/30a-probes/base_a7154ec.txt`.
- The meshes (58 and 134 MB) are not committed. They were written to the
  session scratchpad. To regenerate them, set `S=` at the top of
  `scripts/stats.sh` and rerun `scripts/run_all.sh`.

## 1. Same output

| check | Numedalslågen | Skiensvassdraget |
|---|---|---|
| `.vtk` sha256, 3 base + 3 branch runs | 1 distinct value: `34f7117e5e528e97…` | 1 distinct value: `ab996019190166f9…` |
| matches the probe's base `vtk=` hash | yes | yes |
| probe `mesh` line equal to the base line (codes, counts, vtk) | yes | yes |
| `land cover:` stderr line, all 6 runs | 1860 areas; 0 / 0 / 0 | 2881 areas; 0 / 0 / 0 |

Probe `fixtures` mode on the branch: 107 passed, 1 skipped. Every one of the
base file's 68 `fixture`/`mesh` lines (66 fixture lines and 2 mesh lines) is
on the branch unchanged. The branch has 6 more `fixture` lines. All of them
come from tests that 30a added (`TestCallerPolygons`, `TestDegenerateMeshes`,
which are absent from `test_landcover.py` at `a7154ec`), so they are new and
not compared.

Quality is the same, because the meshes are byte-identical. From the stats
files: worst angle 0.000399° and 3.15e-05°, max vertex degree 19 and 20,
largest height error 9.99997 m and 9.99999 m (tolerance 10 m), the same on
both sides. The Delaunay check is not part of `mesh --stats`. It was not run
separately, because the mesh bytes are unchanged.

## 2. Time (median of 3, AC, seconds)

| catchment | phase | base `a7154ec` | branch `3066d60` | branch / base |
|---|---|---|---|---|
| Numedalslågen (1.29 M triangles) | **land cover** | 2.847 | **0.674** | **0.237 (4.22× faster)** |
| | total | 13.658 | 11.407 | 0.835 (−2.25 s) |
| Skiensvassdraget (3.03 M triangles) | **land cover** | 6.539 | **1.744** | **0.267 (3.75× faster)** |
| | total | 20.458 | 15.679 | 0.766 (−4.78 s) |

Per-repeat values: land cover, base 2.847 / 2.870 / 2.825 and branch
0.688 / 0.673 / 0.674 (Numedalslågen); base 6.565 / 6.522 / 6.539 and branch
1.744 / 1.800 / 1.741 (Skiensvassdraget). The spread (largest minus smallest,
over the median) is 3.4 % on Skiensvassdraget's branch land cover and 2.2 % or
less on every other row.

**Gate (branch median ≤ a third of the base median):** passed on both,
at 0.237 and 0.267.

The total falls by about the land-cover saving (2.17 s and 4.80 s); the
other phases were not compared row by row. For reference, the
battery-powered profile (`../bottlenecks/`) gave land cover 2.868 s and 6.694 s
at the base, and the design's prototype gave 0.71 s and 1.80 s on battery.
These are a different power state, so they are quoted and not compared.

## Files

- `scripts/stats.sh`, `scripts/run_all.sh`, `scripts/summarize.py`: the
  scripts that produced this, as run.
- `raw/stats/`: each run's `--stats` file, stderr, power readings, and the
  `.vtk` sha256 list.
- `raw/summary.txt`: the output of `summarize.py`.
- `raw/probe_*`: the probe's output on the branch, with power readings.
