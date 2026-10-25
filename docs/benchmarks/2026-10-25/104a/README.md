# 104a acceptance: the triangle budget (@perf, 2026-10-25)

Head `worktree-tri-budget` at `632cfce5`. Base: `origin/master` at
`6abfbc49`, run from a detached scratch worktree; `git diff 97759b7b
6abfbc49 -- include src_python bindings CMakeLists.txt tools/bench.py` is
empty, so the base's code is the branch's fork point. The rule is
`docs/increments/104-triangle-budget.md`, section 11 (104a): the 1 m
benchmark's mesh hashes equal the baseline's and timings within noise; on
Kultsjön the budget run at most 1.5 s slower than `--tolerance 10`, with
the predictions P1-P4 checked.

Apple M1 Max (8 P + 2 E cores), 32 GiB, macOS 27.0.1. `tools/bench.py`
blob `16ad2ff5`, `--hardening on` (bounds checks `libc++ fast`). Nothing of
mine ran beside a timing run; another Claude session's process (about 10 %
of one core) was running throughout. Power: `pmset -g batt`, in each
`run.json` and in each line of `kultsjon/runs.jsonl`.

## Verdict

**ACCEPTED.** Every mesh a `--tolerance` run writes is byte-identical to
master's, timings are within noise of master measured back to back, and the
budget run on Kultsjön lands exactly as predicted (P1-P4) and costs 0.82 s
over the tolerance run, under the 1.5 s allowed.

## 1. Bench (`tools/bench.py`, the uncut path, 1 m tolerance)

The three battery runs were made at 92-91 %; the charger was plugged in
during `head-r2`, which bench.py marks `mixed` (not comparable). The
comparison was then remade on AC, head-base-head back to back.

| run | power | commit | `_core` sha256 | tile hash | quarter hash | tile t0 | tile t1 | tile t12 | tile t20 | tile ceiling | quarter t0 | quarter t1 | quarter t12 | quarter ceiling | bench.py verdict |
|---|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| head-r1 | battery 92 % | 632cfce5 | 52390c7f | 11741a81 | a60fb597 | 0.2027 | 0.5250 | 0.1981 | 0.2027 | 2.69x | 0.1808 | 0.4706 | 0.1840 | 2.67x | NO BASELINE (no stored battery run) |
| base-r1 | battery 91 % | 6abfbc49 | 52390c7f | 11741a81 | a60fb597 | 0.1997 | 0.5106 | 0.1984 | 0.2002 | 2.61x | 0.1772 | 0.4477 | 0.1782 | 2.57x | ACCEPTED against head-r1 |
| head-r2 | mixed | 632cfce5 | 52390c7f | 11741a81 | a60fb597 | 0.1987 | 0.5121 | 0.1979 | 0.1993 | 2.60x | 0.1754 | 0.4512 | 0.1756 | 2.60x | NO BASELINE (mixed) |
| head-ac-r1 | AC 90 % | 632cfce5 | 52390c7f | 11741a81 | a60fb597 | 0.1984 | 0.5192 | 0.1994 | 0.2082 | 2.61x | 0.1808 | 0.4596 | 0.1790 | 2.64x | REGRESSION against 23c-4 head-r2 (2026-10-19), 21 cells +5 to +9.5 % |
| base-ac-r1 | AC 91 % | 6abfbc49 | 52390c7f | 11741a81 | a60fb597 | 0.2061 | 0.5081 | 0.2008 | 0.2005 | 2.57x | 0.1790 | 0.4535 | 0.1774 | 2.61x | ACCEPTED against head-ac-r1 |
| head-ac-r2 | AC 92 % | 632cfce5 | 52390c7f | 11741a81 | a60fb597 | 0.2027 | 0.5081 | 0.1986 | 0.2023 | 2.57x | 0.1788 | 0.4683 | 0.1816 | 2.66x | ACCEPTED against base-ac-r1 |

Refine medians in seconds, 5 repeats; t0 is the core's default thread
count; the ceiling is t1 over the fastest swept count.

- **Mesh hashes unchanged**: tile `11741a81…`, quarter `a60fb597…` in all
  six runs, equal to 23c-4's accepted runs. Quality identical: worst angle
  0.6296° / 0.3955°, max degree 74 / 18, within tolerance, 0 Delaunay
  violations of 692,056 / 641,792 edges.
- **The `_core` is the same file** on head and base (sha256 `52390c7f…`),
  and the same as 23c-4's: 104a changes no C++, so `refine` cannot differ.
- **head-ac-r1's REGRESSION is machine state, not code**: the base, with the
  same `.so`, measured right after it at the same level (quarter t0 0.1790
  against 23c-4's 0.1668), and head-ac-r2 is ACCEPTED against it. Today's
  AC level is 5-9 % above 2026-10-19's on both trees; thought to be the
  other session's process and the battery charging, not profiled.

## 2. Kultsjön: `--tolerance 10` against `--triangles 65893` (AC)

Section 7's command (GLO-30, catchment `kultsjon.geojson`, CORINE features,
4 threads), to `.gpkg` with `--stats` and `--record`; `kultsjon_runs.py`,
5 repeats interleaved base / head / budget, then one `.vtk` each. Wall is
the whole process (interpreter start included); total is `--stats`' total.
All 18 runs on AC (93 %, charging).

| run | triangles | tolerance (m) | largest error (m) | trials | wall median (min-max) s | `--stats` total median (min-max) s |
|---|---:|---:|---:|---:|---|---|
| master `--tolerance 10` | 65,893 | 10 | 9.9998 | | 3.217 (3.216-3.281) | 2.913 (2.906-2.955) |
| head `--tolerance 10` | 65,893 | 10 | 9.9998 | | 3.197 (3.162-3.237) | 2.896 (2.869-2.915) |
| head `--triangles 65893` | 67,494 | 9.7885 (found) | 9.7885 | 5 | 4.013 (3.977-4.080) | 3.714 (3.666-3.762) |

- **Byte-identical**: the `.vtk` of master and head at 10 m are equal
  (`cmp`; sha256 `a24ae336…`). The GeoPackages differ only in the `run`
  table's `command` row, which holds the output path (`base-t10` against
  `head-t10`); all 64 rows else and every other table are equal.
- **Not slower**: head's tolerance run is 0.02 s faster than master's
  median, within the spread.
- **The budget costs 0.82 s** (total 3.714 against 2.896), under the
  1.5 s allowed and under P3's 5.5 s (wall 4.0 s). The design estimated
  0.6 s for the four extra trials; the 3.37 m trial (291,983 triangles) is
  the dearest. `--stats` says "The triangle budget ran 5 trials: the
  refine, edge strip, final check and trim rows sum all of them", and its
  refine row is 0.259-0.266 s against 0.046 s for one run.
- **Same mesh as `--tolerance` at the tolerance found** (design test A4 on
  Kultsjön): `--tolerance 9.788516166348431` writes a `.vtk` byte-identical
  to the budget run's (sha256 `6c77a035…`).
- Quality, both from `--stats`: worst angle 0.426° both; median smallest
  angle 34.26° (10 m) and 34.32° (budget); max vertex degree 12 both; 0 flat
  triangles and 0 edges failing the empty-circle test as written in both.

### The error-against-count points (the trials)

Each trial tolerance re-run as `--tolerance T` (`kultsjon_points.py`), T
rebuilt from the floor's largest error E_F = 215.974 m; counts equal the
record's `budget_trials` on every row.

| trial | T (m) | triangles | largest error (m) |
|---:|---:|---:|---:|
| 1 (floor) | 10⁶ | 28,561 | 215.974 |
| 2 | 53.994 (E_F / 4) | 29,439 | 53.900 |
| 3 | 13.498 (E_F / 16) | 49,263 | 13.498 |
| 4 | 3.3746 (E_F / 64) | 291,983 | 3.3745 |
| 5, chosen | 9.7885 | 67,494 | 9.7885 |
| the 10 m mesh | 10 | 65,893 | 9.9998 |

### Predictions P1-P4

| prediction | measured | |
|---|---|---|
| P1: 67,494 triangles at 9.79 m in 5 trials (T 9.78-9.80, count 67,400-67,600) | 67,494, T 9.7885, largest error 9.7885 m, 5 trials | met |
| P1, 2 % cap (`--max-triangles 67210`): 65,764 at 10.02 m, 6 trials | 65,764, T 10.0202, largest error 10.0197 m, 6 trials | met |
| P1, RMS 3.00 m | not measured: RMS is 104b | open |
| P2: the trials lie on the ladder's curve | the floor equals the ladder's 28,561; every other trial lies between its ladder neighbours (53.99 m: 29,439 between 28,742 and 30,590; 13.50 m: 49,263 between 38,228 and 62,521; 3.37 m: 291,983 between 165,394 and 447,064; 9.79 m: 67,494 between 65,893 and 69,892) | met |
| P3: under 5.5 s | wall 4.0 s, total 3.7 s | met |
| P4: `--triangles 30000` returns the floor, 28,561 at 215.97 m, 1 trial | 28,561, tolerance 215.974 m, 1 trial | met |
| P4: `--triangles 27000` refused, exit 2, nothing written | exit 2, no mesh, `--stats` or `--record` file; message names 28561 triangles and 216 m and "at least 28561" | met |

Two wordings seen on the way, reported, not judged here:

- With the 2 % cap the summary says "tolerance 10 m found in 6 trials.
  Every node of the original DEM inside the mesh is within 10 m of it
  (largest difference 10.02 m)": the three-figure rounding of 10.0202 reads
  as a largest difference over the tolerance.
- The `--triangles 27000` refusal is headed "Invalid value for --dem".

## 3. The suite, quiet

AC, nothing else of mine running, after rebuilding `build-pyext` and copying
`_core` into the venv: `pytest -n auto`: 10,061 passed, 18 skipped,
1 xfailed, 154 s. `test_core_cdt.py::TestGil::test_the_ticker_distinguishes_a_held_gil_from_a_released_one`
passed in that run and alone.

## Files and regeneration

- `../104a-*/`: bench.py's evidence (`run.json`, `raw.tsv`, `README.md`).
  Its meshes were in `$TMPDIR/rasputin-bench/<label>` and are not kept;
  `python tools/bench.py run --label <label> --tree <tree>` regenerates them.
- `kultsjon_runs.py`, `kultsjon_points.py`: the Kultsjön drivers; they need
  each tree's `build-bench/pkg` from a bench.py run. Output:
  `kultsjon/runs.jsonl` (wall, exit, stderr tail, power before and after)
  and each run's `--stats` and `--record`. Meshes went to
  `/tmp/perf-104a-meshes/` (not kept); rerun the drivers to regenerate.
