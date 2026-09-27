# 1 m mesher comparison, increments 14 to HEAD

Kartverket tile `7908_3_10m_z33.tif` (5051 x 5051, 10 m), tolerance 1 m.
Machine: Apple Silicon, 10 cores, **on battery** (`pmset -g batt`: 100 %, discharging).
Threads: 10. This is the default, `hardware_concurrency`. No commit has a thread
option, so the count cannot be set any other way.
Build: Release, `-O3 -DNDEBUG`, same compiler and pybind11 for every commit.

Commits measured. Where an increment landed through a merge, the merge commit is used:
14 `935271c`, 14b `333bcf8`, 16 `1018901`, 17 `2760d36`, 18 `a9a93bc`,
20 `f5489ab`, 20b before the Lawson fix `75a5b6c`, HEAD `08a7711`.
`head_nofeet_q0` is HEAD run with `--start-min-angle 0 --no-constraint-feet`.

Timing columns:
- **refine s** is the `_core.refine` call alone: the median of 3 `--binary` runs, with the min and max in brackets.
- **app s** is the whole CLI command body: decode, refine, trim and the binary write. Interpreter start-up and imports are not included.
- **process s** is wall time for the whole process.

All quality columns are computed by one script (`quality.py`) from the ASCII `.vtk`.
They are plan view (x/y), and degree means incident triangles.
Achieved max error is the value `refine` reports, which the tests prove equals the per-DEM-node oracle.

## Quarter circle (`--domain quarter.geojson`, 16 onward)

| label | refine s (med of 3) | app s | process s | triangles | achieved max err m | rounds | flips | min angle median | % < 1° | n < 1° | n < 0.1° | worst ° | max degree | deg ≥ 12 | deg ≥ 20 | CDT violations (edges) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| i16 | 3.291 (3.268–3.308) | 3.406 | 3.645 | 427779 | 0.999998 | 40 | 452067 | 45.00 | 0.028 | 119 | 3 | 0.0117 | 43 | 216 | 9 | 0 / 641137 |
| i17 | 3.323 (3.310–3.339) | 3.425 | 3.683 | 427779 | 0.999998 | 40 | 452067 | 45.00 | 0.028 | 119 | 3 | 0.0117 | 43 | 216 | 9 | 0 / 641137 |
| i18 | 0.373 (0.373–0.374) | 0.472 | 0.717 | 427779 | 0.999998 | 40 | 452067 | 45.00 | 0.028 | 119 | 3 | 0.0117 | 43 | 216 | 9 | 0 / 641137 |
| i20 | 0.236 (0.236–0.237) | 0.337 | 0.578 | 428225 | 0.999980 | 41 | 445653 | 45.00 | 0.004 | 19 | 3 | 0.0117 | 18 | 185 | 0 | 0 / 641806 |
| i20b | 0.241 (0.240–0.248) | 0.342 | 0.607 | 428217 | 0.999980 | 41 | 445657 | 45.00 | 0.003 | 12 | 0 | 0.3955 | 18 | 184 | 0 | 0 / 641791 |
| head | 0.238 (0.237–0.251) | 0.339 | 0.587 | 428217 | 0.999980 | 41 | 445657 | 45.00 | 0.003 | 12 | 0 | 0.3955 | 18 | 184 | 0 | 0 / 641791 |
| head_nofeet_q0 | 0.378 (0.377–0.381) | 0.475 | 0.721 | 427779 | 0.999998 | 40 | 452067 | 45.00 | 0.028 | 119 | 3 | 0.0117 | 43 | 216 | 9 | 0 / 641137 |

Identical meshes, compared as md5 of everything from `POINTS` on:
- 16 = 17 = 18 = HEAD with both features off.
- 20b (`75a5b6c`) = HEAD. The Lawson fix did not change this mesh.

Cross-check:
- HEAD's `--stats` agrees with `quality.py`: worst 0.396°, 12 under 1°, max degree 18, 184 vertices at degree ≥ 12.
- Agreement with the earlier session numbers:
  - 17 app time is 3.43 s (earlier 3.39 s).
  - 18 app time is 0.47 s (earlier 0.48 s).
  - 20 and 20b app time is 0.34 s (earlier about 0.33 s).
  - HEAD with feet has worst angle 0.396° and 0 triangles under 0.1°.
- No disagreement.

## Whole tile (no domain; 14 and 14b have no domain option)

The start grid is each commit's default `refine_start_stride`:
- 158 at 14.
- 40 from 14b on.

`i14_stride40` reruns 14 with `--stride 40`, so it can be compared like for like with 14b.

| label | refine s (med of 3) | app s | process s | triangles | achieved max err m | rounds | flips | min angle median | % < 1° | n < 1° | n < 0.1° | worst ° | max degree | deg ≥ 12 | deg ≥ 20 | CDT violations (edges) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| i14 | 0.379 (0.370–0.380) | 0.481 | 0.723 | 670554 | 0.999959 | 156 | ? | 2.11 | 35.775 | 239892 | 59524 | 0.0012 | 768 | 23859 | 8577 | 136094 / 1001542 |
| i14_stride40 | 0.161 (0.161–0.162) | 0.319 | 0.545 | 591537 | 0.999992 | 38 | ? | 5.19 | 17.731 | 104886 | 6413 | 0.0188 | 236 | 20584 | 7611 | 124570 / 882925 |
| i14b | 0.268 (0.268–0.276) | 0.426 | 0.667 | 463974 | 0.999985 | 53 | 446326 | 45.00 | 0.000 | 2 | 0 | 0.6296 | 74 | 324 | 196 | 0 / 691582 |
| i16 | 0.282 (0.280–0.286) | 0.437 | 0.680 | 463974 | 0.999985 | 53 | 446326 | 45.00 | 0.000 | 2 | 0 | 0.6296 | 74 | 324 | 196 | 0 / 691582 |
| i17 | 0.288 (0.282–0.288) | 0.442 | 0.693 | 463974 | 0.999985 | 53 | 446326 | 45.00 | 0.000 | 2 | 0 | 0.6296 | 74 | 324 | 196 | 0 / 691582 |
| i18 | 0.261 (0.260–0.263) | 0.417 | 0.661 | 463974 | 0.999985 | 53 | 446326 | 45.00 | 0.000 | 2 | 0 | 0.6296 | 74 | 324 | 196 | 0 / 691582 |
| i20 | 0.264 (0.263–0.271) | 0.421 | 0.681 | 464290 | 0.999976 | 53 | 445675 | 45.00 | 0.000 | 1 | 0 | 0.6296 | 74 | 324 | 197 | 0 / 692056 |
| i20b | 0.266 (0.264–0.266) | 0.416 | 0.662 | 464290 | 0.999976 | 53 | 445675 | 45.00 | 0.000 | 1 | 0 | 0.6296 | 74 | 324 | 197 | 0 / 692056 |
| head | 0.264 (0.263–0.265) | 0.414 | 0.667 | 464290 | 0.999976 | 53 | 445675 | 45.00 | 0.000 | 1 | 0 | 0.6296 | 74 | 324 | 197 | 0 / 692056 |
| head_nofeet_q0 | 0.265 (0.263–0.275) | 0.420 | 0.678 | 463974 | 0.999985 | 53 | 446326 | 45.00 | 0.000 | 2 | 0 | 0.6296 | 74 | 324 | 196 | 0 / 691582 |

Identical meshes:
- 14b = 16 = 17 = 18 = HEAD with both features off.
- 20 = 20b = HEAD.

On the whole tile, the quality start inserts 503 nodes and skips 252. The grid start is all 45° triangles, so these most likely come from the tile's nodata edge. The report does not show this; it is inferred. No constraint feet fire.

## Regressions

A regression here is a step where time rose by more than 10 % or a quality measure got worse. Each step is compared with the previous one, like for like.

- **14 to 14b, whole tile, time +66 %** (0.161 s to 0.268 s refine, both at stride 40). Likely cause: the Lawson flips that 14b adds, 446k of them. This is the intended cost of the feature. In return, CDT violations drop from 124 570 to 0, triangles under 1° drop from 17.7 % to 0.0 %, and max degree drops from 236 to 74.
- **14b to 16 to 17, whole tile, time +5 % then +2 %** (0.268 s, 0.282 s, 0.288 s). Both steps are under the 10 % threshold. Increment 18 then brings it back to 0.261 s. Likely cause: off-node start support in the refine loop (16) and phase timing (17). Not flagged.
- **16 to 17, quarter, +1 %.** Noise. Not flagged.
- **18 to 20, whole tile: deg ≥ 20 rises from 196 to 197.** Every other measure is equal or better, including triangles under 1°, which drop from 2 to 1. Likely cause: the 503 quality-start nodes near the nodata edge change the local topology. Negligible.
- **20 to 20b, quarter, refine +2 %** (0.236 s to 0.241 s). This is within the run-to-run spread. Not flagged.
- **HEAD with features off compared with 18, quarter, +1 %** (0.378 s against 0.373 s) on a byte-identical mesh. The code added in 20 and 20b therefore costs about nothing when it is switched off.

No step made the tolerance worse: every run achieves ≤ 1 m. No run violates constrained Delaunay from 14b on. Increment 14 violates it by design, because it has no flips.

Gains, for context:
- 18 makes the quarter circle 8.8x faster (refine 3.29 s to 0.37 s).
- 20 cuts the quarter-circle refine time by another 37 % (0.373 s to 0.236 s). The quality start leaves fewer bad triangles to split.
- 20 and 20b take the quarter circle's max degree from 43 to 18, degree ≥ 20 from 9 to 0, and the worst angle from 0.0117° to 0.396°.

## Method

The steps were run in this order:
1. `build.sh` creates one `git worktree --detach` per commit under `wt/`. It builds `_core` with
   `cmake -DCMAKE_BUILD_TYPE=Release -DRASPUTIN_BUILD_PYTHON=ON -DRASPUTIN_BUILD_TESTS=OFF -Dpybind11_DIR=<uv pybind11> -DPython_EXECUTABLE=.venv/bin/python`
   and copies the `.so` into the worktree's `src_python/tin_engine/`.
2. `run.py WT args…` removes the venv's scikit-build editable finder from `sys.meta_path`, puts `WT/src_python` first on the path, and asserts that `tin_engine` comes from WT. It also wraps `cli.refine` and `cli.refine_start_stride` to print their timing and stride, and then runs the Typer `app`.
3. `bench.sh DOMAIN LABEL WT [extra…]` runs the following from the repo root:
   - three timed runs:
     `.venv/bin/python run.py WT mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif [--domain …/paraview/inputs/quarter.geojson] --tolerance 1 --binary --out logs/tmp.vtk [extra]`
   - then one run with `--ascii --out runs/<dom>_<label>.vtk`, plus `--stats runs/<dom>_<label>.stats.md` where the commit has `--stats` (17 on).
4. `quality.py` is run on the ASCII file.
5. `summarize.py` builds the tables from `logs/`.

The CDT check in `quality.py` has three parts:
- It visits every interior edge that is not a constraint `LINES` edge, and tests the apex of one side against the circumcircle of the other side's CCW triangle.
- The float determinant is taken on coordinates relative to the apex. When it falls within 4x Shewchuk's error bound, an exact `Fraction` evaluation on the written doubles decides.
- Strictly inside counts as a violation, once per edge. The C++ tests count once per side, so their count can be up to 2x this.

The tests evaluate in the frame (col·dx, −row·dy). This script uses world x/y, which are the same frame translated by the tile origin. For lattice nodes that translation is exact. For off-node vertices (the domain boundary and feet) the written doubles are used as they are.

Worktrees removed after the run. Scripts: `build.sh`, `run.py`, `bench.sh`, `quality.py`, `summarize.py`. Raw logs: `logs/`. Meshes: `runs/<domain>_<label>.vtk` (text) and the `--stats` reports as `runs/*.stats.md`.
