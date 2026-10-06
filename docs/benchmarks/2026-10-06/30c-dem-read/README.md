# Increment 30c acceptance: the `decode` phase of `rasputin mesh`

@perf, 2026-10-06. Base: `6c729e97` (30b's approved head; its DEM code equals
master's). Branch: `worktree-dem-read-speed` at `0663efc9` (the green step).
The acceptance is section 9 of `docs/increments/30c-dem-read-speed.md`.

**Verdict: ACCEPTED.** The branch's `decode` median is under the design's
limit of 0.7 of the base's on both catchments (the gate table below). Every
mesh the 12 timing runs wrote is byte-identical to the base's. The design's
probe finds every line of `docs/increments/30c-probes/base_6c729e97.txt`
present and unchanged on the branch. The only new lines come from the red
suite's tests.

Every figure in the tables between the `GENERATED` markers is written there by
`scripts/summarize.py` from `raw/`; none is copied by hand. Rerun it to check
them.

## Method

- **Machine:** Apple M1 Max (8 performance and 2 efficiency cores), 32 GB,
  macOS 27.0. Python 3.14.7, tifffile 2026.9.20, imagecodecs 2026.8.16,
  numpy 2.5.3, shapely 2.1.2 with GEOS 3.13.1, the same on both sides.
  `refine` runs on 10 threads (each stats file says so).
- **Power: AC for every kept run.** `scripts/stats.sh` reads `pmset -g batt`
  before and after each run. It also counts the "Using Batt" lines in
  `pmset -g log` from the run's start, which would catch a switch to battery
  and back inside the run. A run that is not on AC throughout is moved to
  `raw/discarded/` and the set stops. The log filter was checked against an
  earlier window of the same day that was on battery: it found 5 lines there.
  None of the runs was discarded for power.
- **Two runs discarded for another reason.** While the first base run on
  Numedalslågen ran, my own shell was polling in an empty busy loop on one
  core. That run and its paired branch run are in `raw/discarded/`
  (`README.txt` there says so). Both were redone in the same order, base
  first, after the other 10 runs.
- **Installs**, both non-editable, each in its own scratch venv, built with
  `uv pip install ".[codecs]"` from `git archive <commit>`, with pytest,
  pytest-asyncio and hypothesis added for the probe:
  - Before the runs, the installed `tin_engine/*.py` of each venv matched
    `git archive <commit> src_python/tin_engine` byte for byte (`cmp`, 54
    files each, no difference). The two installs' `io/cog.py` differ.
  - The only source difference between the commits is in
    `src_python/tin_engine/io/cog.py` and `src_python/tin_engine/mosaic.py`
    (`git diff --stat 6c729e97 0663efc9 -- src_python include src
    CMakeLists.txt pyproject.toml`).
  - Both venvs load a `_core` extension with the same sha256 (`19d71363…`).
    It is a Release build with bounds checks on (libc++ fast mode), as each
    stats file says.
  - Each run calls its venv's own `python` with `tin_engine.cli.app`, never
    a `bin/rasputin` entry script (the trap in section 8 of the design). Each
    run's stderr starts with the `tin_engine` it imported, and the quality
    table below lists them.
- **Command:** `scripts/stats.sh <side> <catchment> <repeat>`, driven by
  `scripts/run_all.sh`. It runs `rasputin mesh` with DEM
  `rasputin_data/DTM10_UTM33_20260925`, domain
  `rasputin_scratch/norway/<c>/<c>_outline_nve.geojson`, and `--features
  rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018
  --features-map corine --tolerance 10 --binary --stats`. These are the inputs
  and flags of 30a's and 30b's acceptance runs.
- **Repeats:** 3 per catchment per side, base and branch alternated: base
  first on odd repeats, branch first on even ones. Nothing else of mine ran
  during the kept timing runs. The tables give the median.
- **Byte-identity:** the sha256 of every `.vtk` the 12 kept runs wrote
  (`raw/stats/vtk_sha256.txt`). After the timing runs came the design's probe,
  `docs/increments/30c-probes/dem_bytes.py`, in `fixtures` mode and `mesh`
  mode (`scripts/probe.sh`). It was run on both installs from the worktree
  root, so both ran the branch's `tests/python`. Each run was compared with
  `base_6c729e97.txt` by section 6's two `comm` commands
  (`raw/probe_compare.txt`).
  - The base install reproduces the base file: no base line missing. It
    also gives the same 16 new lines as the branch. On the base, 3 tests fail:
    R1's `seven_cores` and `cores_unknown`, and R2. These are the red tests,
    as expected (`FAILED` lines in `raw/probe_fixtures_base.txt`).
  - The comparison was checked once by planting a changed `vtk=` value in a
    copy of the base file: the first `comm` printed that line.
- The meshes (about 58 and 134 MB each, 1.1 GB in all) are not committed.
  They were written to the session scratchpad. To regenerate them, rebuild
  the two venvs as above, set `S=` in `scripts/stats.sh`, and rerun
  `scripts/run_all.sh` (it skips runs whose stats file exists).

## Against the design's prototype

Section 8 of the design measured a prototype the same way: `decode` 1.794 →
1.057 s on Numedalslågen and 1.307 → 0.815 s on Skiensvassdraget. These
medians are within 0.03 s of that on both sides. The whole-run totals here
are 0.2 to 0.35 s below the design's. This is the same machine on AC power,
but in a different session; the cause was not measured.

## Results

<!-- GENERATED by scripts/summarize.py: do not edit by hand -->

Power: 24 of 24 `pmset -g batt` readings (before and after each of 12 kept runs) say AC Power; `pmset -g log` holds 0 "Using Batt" lines inside those runs. Runs discarded for power: 0; for a busy loop of @perf's beside the run: 2 (`raw/discarded/README.txt`).

### Time (median of 3, seconds)

| catchment | phase | base `6c729e97` | branch `0663efc9` | branch / base |
|---|---|---|---|---|
| Numedalslågen | **decode** | 1.788 | **1.051** | **0.588** (-0.737 s) |
|  | total | 7.480 | 6.719 | 0.898 (-0.761 s) |
| Skiensvassdraget | **decode** | 1.323 | **0.842** | **0.636** (-0.481 s) |
|  | total | 11.462 | 10.966 | 0.957 (-0.496 s) |

### Gate: branch `decode` median at most 0.7 of the base's

| catchment | branch / base | verdict |
|---|---|---|
| Numedalslågen | 0.588 | pass |
| Skiensvassdraget | 0.636 | pass |

### Per repeat (seconds, repeats 1 / 2 / 3; spread = (max − min) / median)

| catchment | phase | side | values | spread |
|---|---|---|---|---|
| Numedalslågen | decode | base | 1.822 / 1.788 / 1.776 | 2.6 % |
| Numedalslågen | decode | branch | 1.037 / 1.053 / 1.051 | 1.5 % |
| Numedalslågen | total | base | 7.516 / 7.480 / 7.440 | 1.0 % |
| Numedalslågen | total | branch | 6.695 / 6.762 / 6.719 | 1.0 % |
| Skiensvassdraget | decode | base | 1.325 / 1.323 / 1.306 | 1.4 % |
| Skiensvassdraget | decode | branch | 0.842 / 0.807 / 0.842 | 4.2 % |
| Skiensvassdraget | total | base | 11.508 / 11.462 / 11.393 | 1.0 % |
| Skiensvassdraget | total | branch | 11.016 / 10.925 / 10.966 | 0.8 % |

### Every row of the `--stats` timing table (median of 3, seconds; recorded, not gated)

| catchment | phase | base | branch | branch − base |
|---|---|---|---|---|
| Numedalslågen | domain read | 0.021 | 0.014 | -0.007 |
|  | decode | 1.788 | 1.051 | -0.737 |
|  | features read | 0.265 | 0.263 | -0.002 |
|  | features clip | 1.816 | 1.818 | +0.002 |
|  | start mesh: build | 0.005 | 0.006 | +0.001 |
|  | start mesh: node | 0.482 | 0.482 | +0.000 |
|  | start mesh: triangulate | 0.142 | 0.142 | +0.000 |
|  | start mesh: constraint edges | 0.345 | 0.342 | -0.003 |
|  | refine | 1.109 | 1.100 | -0.009 |
|  | refine: legalise start | 0.008 | 0.008 | +0.000 |
|  | refine: start quality | 0.510 | 0.504 | -0.006 |
|  | refine: scan (parallel) | 0.287 | 0.286 | -0.001 |
|  | refine: split + flip (serial) | 0.170 | 0.168 | -0.002 |
|  | refine: setup + output | 0.136 | 0.135 | -0.001 |
|  | edge strip: generate | 0.113 | 0.113 | +0.000 |
|  | edge strip: scan (parallel) | 0.013 | 0.013 | +0.000 |
|  | edge strip: split + flip (serial) | 0.010 | 0.010 | +0.000 |
|  | trim | 0.051 | 0.051 | +0.000 |
|  | land cover | 0.676 | 0.678 | +0.002 |
|  | write: encode | 0.021 | 0.021 | +0.000 |
|  | write: disk | 0.008 | 0.008 | +0.000 |
|  | other | 0.604 | 0.598 | -0.006 |
|  | total | 7.480 | 6.719 | -0.761 |
| Skiensvassdraget | domain read | 0.016 | 0.023 | +0.007 |
|  | decode | 1.323 | 0.842 | -0.481 |
|  | features read | 0.288 | 0.288 | +0.000 |
|  | features clip | 1.337 | 1.341 | +0.004 |
|  | start mesh: build | 0.009 | 0.009 | +0.000 |
|  | start mesh: node | 1.013 | 1.018 | +0.005 |
|  | start mesh: triangulate | 0.271 | 0.270 | -0.001 |
|  | start mesh: constraint edges | 0.624 | 0.618 | -0.006 |
|  | refine | 2.919 | 2.932 | +0.013 |
|  | refine: legalise start | 0.014 | 0.014 | +0.000 |
|  | refine: start quality | 1.028 | 1.030 | +0.002 |
|  | refine: scan (parallel) | 0.829 | 0.832 | +0.003 |
|  | refine: split + flip (serial) | 0.774 | 0.782 | +0.008 |
|  | refine: setup + output | 0.277 | 0.277 | +0.000 |
|  | edge strip: generate | 0.202 | 0.203 | +0.001 |
|  | edge strip: scan (parallel) | 0.027 | 0.027 | +0.000 |
|  | edge strip: split + flip (serial) | 0.025 | 0.025 | +0.000 |
|  | trim | 0.117 | 0.116 | -0.001 |
|  | land cover | 1.787 | 1.777 | -0.010 |
|  | write: encode | 0.044 | 0.046 | +0.002 |
|  | write: disk | 0.019 | 0.025 | +0.006 |
|  | other | 1.414 | 1.407 | -0.007 |
|  | total | 11.462 | 10.966 | -0.496 |

### Same output over all 6 runs per catchment

| catchment | distinct `.vtk` sha256 | distinct `dem_seams` values | exit |
|---|---|---|---|
| Numedalslågen | 1: `34f7117e5e528e97…` | 1 | exit 0 |
| Skiensvassdraget | 1: `ab996019190166f9…` | 1 | exit 0 |

### Quality and installs (from each `--stats` and stderr file; distinct values over 3 repeats)

| catchment | side | `tin_engine` imported (under the scratchpad) | worst angle | max vertex degree | largest height error, m | threads | bounds checks |
|---|---|---|---|---|---|---|---|
| Numedalslågen | base | `venv-base/lib/python3.14/site-packages/tin_engine/__init__.py` | 0.000399° | 19 | 9.999971563165786 | 10 | on (libc++ fast) |
| Numedalslågen | branch | `venv-branch/lib/python3.14/site-packages/tin_engine/__init__.py` | 0.000399° | 19 | 9.999971563165786 | 10 | on (libc++ fast) |
| Skiensvassdraget | base | `venv-base/lib/python3.14/site-packages/tin_engine/__init__.py` | 3.15e-05° | 20 | 9.999987707205037 | 10 | on (libc++ fast) |
| Skiensvassdraget | branch | `venv-branch/lib/python3.14/site-packages/tin_engine/__init__.py` | 3.15e-05° | 20 | 9.999987707205037 | 10 | on (libc++ fast) |

### Power, each timing run (`pmset -g batt` before and after; "Using Batt" lines in `pmset -g log` during the run)

| catchment | side | repeat | before | after | battery lines |
|---|---|---|---|---|---|
| Numedalslågen | base | 1 | 19:26:22 AC Power, 100%; charged | 19:26:30 AC Power, 100%; charged | 0 |
| Numedalslågen | base | 2 | 19:16:26 AC Power, 100%; charged | 19:16:35 AC Power, 100%; charged | 0 |
| Numedalslågen | base | 3 | 19:17:14 AC Power, 100%; charged | 19:17:22 AC Power, 100%; charged | 0 |
| Numedalslågen | branch | 1 | 19:26:35 AC Power, 100%; charged | 19:26:42 AC Power, 100%; charged | 0 |
| Numedalslågen | branch | 2 | 19:16:14 AC Power, 100%; charged | 19:16:22 AC Power, 100%; charged | 0 |
| Numedalslågen | branch | 3 | 19:17:27 AC Power, 100%; charged | 19:17:34 AC Power, 100%; charged | 0 |
| Skiensvassdraget | base | 1 | 19:15:39 AC Power, 100%; charged | 19:15:52 AC Power, 100%; charged | 0 |
| Skiensvassdraget | base | 2 | 19:16:56 AC Power, 100%; charged | 19:17:09 AC Power, 100%; charged | 0 |
| Skiensvassdraget | base | 3 | 19:17:39 AC Power, 100%; charged | 19:17:52 AC Power, 100%; charged | 0 |
| Skiensvassdraget | branch | 1 | 19:15:57 AC Power, 100%; charged | 19:16:09 AC Power, 100%; charged | 0 |
| Skiensvassdraget | branch | 2 | 19:16:39 AC Power, 100%; charged | 19:16:51 AC Power, 100%; charged | 0 |
| Skiensvassdraget | branch | 3 | 19:17:57 AC Power, 100%; charged | 19:18:09 AC Power, 100%; charged | 0 |

### The byte-identical probe (`raw/probe_compare.txt`, verbatim)

```
== base: tin_engine from venv-base/lib/python3.14/site-packages/tin_engine / tin_engine from venv-base/lib/python3.14/site-packages/tin_engine
== base: base lines missing or changed (must be empty):
== base: lines the base file lacks (allowed on the branch: the red suite's new tests, section 6):
fixture tests/python/test_io_cog.py::TestR1DefaultThreadsAreTheMachines::test_the_pool_is_built_with_the_machines_core_count[cores_unknown] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestR1DefaultThreadsAreTheMachines::test_the_pool_is_built_with_the_machines_core_count[seven_cores] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestR1DefaultThreadsAreTheMachines::test_the_pool_is_built_with_the_machines_core_count[threads_given] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[stripped] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[stripped] #2 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_deflate_fp_predictor] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_deflate_fp_predictor] #2 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_deflate] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_deflate] #2 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_int16] #1 decode_window window=5,9,40,55 dtype=float32 sha=3b62a53c66658977
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_int16] #2 decode_window window=5,9,40,55 dtype=float32 sha=3b62a53c66658977
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_lzw] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_lzw] #2 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_mosaic.py::TestSeamRegionOnCountedNodes::test_p1_edge_cases_under_a_needed_region[even_count] #1 assemble shape=6x9 dtype=float64 sha=6049f69784dc1633 seams=[a.tif/b.tif/4/inf/1.5]
fixture tests/python/test_mosaic.py::TestSeamRegionOnCountedNodes::test_p1_edge_cases_under_a_needed_region[odd_count] #1 assemble shape=6x9 dtype=float64 sha=f310b7c7cf07ba90 seams=[a.tif/b.tif/3/inf/1.0]
fixture tests/python/test_mosaic.py::TestSeamRegionOnCountedNodes::test_r2_the_region_is_tested_only_on_nodes_that_could_count #1 assemble shape=6x9 dtype=float32 sha=8bbe171a32a27678 seams=[a.tif/b.tif/3/5.0/3.0]
== base: counts: base file 1180, this run 1196; fixtures: pytest exit 1, 1194 calls recorded; pytest summary: 3 failed, 717 passed, 4 skipped in 36.46s

== branch: tin_engine from venv-branch/lib/python3.14/site-packages/tin_engine / tin_engine from venv-branch/lib/python3.14/site-packages/tin_engine
== branch: base lines missing or changed (must be empty):
== branch: lines the base file lacks (allowed on the branch: the red suite's new tests, section 6):
fixture tests/python/test_io_cog.py::TestR1DefaultThreadsAreTheMachines::test_the_pool_is_built_with_the_machines_core_count[cores_unknown] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestR1DefaultThreadsAreTheMachines::test_the_pool_is_built_with_the_machines_core_count[seven_cores] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestR1DefaultThreadsAreTheMachines::test_the_pool_is_built_with_the_machines_core_count[threads_given] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[stripped] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[stripped] #2 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_deflate_fp_predictor] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_deflate_fp_predictor] #2 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_deflate] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_deflate] #2 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_int16] #1 decode_window window=5,9,40,55 dtype=float32 sha=3b62a53c66658977
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_int16] #2 decode_window window=5,9,40,55 dtype=float32 sha=3b62a53c66658977
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_lzw] #1 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_io_cog.py::TestW2Determinism::test_p2_the_default_thread_count_changes_no_byte[tiled_lzw] #2 decode_window window=5,9,40,55 dtype=float32 sha=e0c2cdb842b96264
fixture tests/python/test_mosaic.py::TestSeamRegionOnCountedNodes::test_p1_edge_cases_under_a_needed_region[even_count] #1 assemble shape=6x9 dtype=float64 sha=6049f69784dc1633 seams=[a.tif/b.tif/4/inf/1.5]
fixture tests/python/test_mosaic.py::TestSeamRegionOnCountedNodes::test_p1_edge_cases_under_a_needed_region[odd_count] #1 assemble shape=6x9 dtype=float64 sha=f310b7c7cf07ba90 seams=[a.tif/b.tif/3/inf/1.0]
fixture tests/python/test_mosaic.py::TestSeamRegionOnCountedNodes::test_r2_the_region_is_tested_only_on_nodes_that_could_count #1 assemble shape=6x9 dtype=float32 sha=8bbe171a32a27678 seams=[a.tif/b.tif/3/5.0/3.0]
== branch: counts: base file 1180, this run 1196; fixtures: pytest exit 0, 1194 calls recorded; pytest summary: 720 passed, 4 skipped in 30.63s
```

Probe power (`raw/probe_power.txt`):

```
base fixtures before 19:27:13: Now drawing from 'AC Power'  -InternalBattery-0 (id=24117347) 100%; charged; 0:00 remaining present: true 
base fixtures after 19:27:50: Now drawing from 'AC Power'  -InternalBattery-0 (id=24117347) 100%; charged; 0:00 remaining present: true 
base mesh before 19:27:50: Now drawing from 'AC Power'  -InternalBattery-0 (id=24117347) 100%; charged; 0:00 remaining present: true 
base mesh after 19:28:12: Now drawing from 'AC Power'  -InternalBattery-0 (id=24117347) 100%; charged; 0:00 remaining present: true 
branch fixtures before 19:28:12: Now drawing from 'AC Power'  -InternalBattery-0 (id=24117347) 100%; charged; 0:00 remaining present: true 
branch fixtures after 19:28:43: Now drawing from 'AC Power'  -InternalBattery-0 (id=24117347) 100%; charged; 0:00 remaining present: true 
branch mesh before 19:28:43: Now drawing from 'AC Power'  -InternalBattery-0 (id=24117347) 100%; charged; 0:00 remaining present: true 
branch mesh after 19:29:03: Now drawing from 'AC Power'  -InternalBattery-0 (id=24117347) 100%; charged; 0:00 remaining present: true
```

<!-- END GENERATED -->

## The probe at the merge commit `26a5d839` (design section 6, "At the merge")

Run 2026-10-06, 19:44 to 19:46, on AC power throughout (`raw/probe_merge_power.txt`), by
`scripts/probe_merge.sh`: `26a5d839` (this branch with master merged in, #199 and #201) installed
non-editable in a scratch venv from `git archive` (Python 3.14.7, tifffile 2026.9.20, imagecodecs
2026.8.16, numpy 2.5.3, shapely 2.1.2, GEOS 3.13.1, the same as the base's), that venv's own `python`
calling the probe from the worktree root (its `tests/python` and probe equal `26a5d839`'s:
`git diff --stat 26a5d839 ef4abf64 -- tests docs/increments/30c-probes/dem_bytes.py` prints nothing).
The probe's first line names the scratch venv's `site-packages/tin_engine`. Both modes ran twice;
the two runs' `fixture` and `mesh` lines are identical, and run 1 is `raw/probe_merge_26a5d839.txt`
(first line: how it was made).

Against `docs/increments/30c-probes/base_6c729e97.txt`, by section 6's two `comm` commands:

| | lines |
|---|---|
| base `fixture` and `mesh` lines | 1,180 |
| base lines missing or changed in the merge's run | 0 |
| lines the base lacks | 16, the same 16 the branch's run at `bc8d91eb` added (`raw/probe_compare.txt`) |
| `.vtk` hashes (Numedalslågen, Skiensvassdraget) | `34f7117e5e528e97`, `ab996019190166f9`, the base's |

The 16 new lines are the red suite's: `test_io_cog.py` `TestR1DefaultThreadsAreTheMachines` (3) and
`TestW2Determinism::test_p2_…` (10), `test_mosaic.py` `TestSeamRegionOnCountedNodes` P1 (2) and R2 (1).
Master's merge adds none: pytest over the 14 suites gives 726 passed, 4 skipped, against 720 on the
branch at `bc8d91eb`; the 6 more are master's two parametrised tests in `test_io_repository.py`
(`test_a_call_to_the_repositorys_own_reader_is_not_an_opener`,
`test_only_the_repositorys_own_names_are_let_through`, 6 cases by `--collect-only`), which call neither
`assemble` nor `decode_window`. The probe recorded 1,194 calls, as on the branch.
