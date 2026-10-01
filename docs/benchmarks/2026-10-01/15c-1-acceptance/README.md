# Increment 15c-1 acceptance (@perf, 2026-10-01)

**Verdict: ACCEPTED.** The refine/mesh acceptance rule applies to 15c-1
(`docs/increments/README.md`, "Acceptance"; `docs/increments/15c-geographic-dem.md`,
"Acceptance"). The mesh is byte-identical to the base's on both domains in all
10 runs. Refine time is within noise at every thread count: pooled over 5
back-to-back pairs, no cell moves by more than 5 %. Quality is unchanged. All
runs were on **battery**, and only battery is compared with battery.

## What was measured

- **Branch** `worktree-15c-1` at `6ab7ad5`. The production code is final.
- **Base** is master's merge commit `a130f7c` (#128, the 15c design), in a
  detached `git worktree` in the job's scratch directory.
  `git diff --stat a130f7c 6ab7ad5 -- include src bindings CMakeLists.txt`
  lists four files:
  - `lattice_mesh.hpp`: `split_inside` now takes a `MeshVertex`.
  - `check_points.hpp` and `refine_points.hpp`: new.
  - `bindings/core.cpp`: +93 lines that bind the new headers.

  The `_core` builds therefore differ: sha256 `eb5f06e88a7c…` for the base and
  `7dee5ad918d0…` for the branch. `refine` does not call the new code.
- **Script:** `tools/bench.py`, blob `77765b18`, identical at both commits. It
  built Release (`-O3 -DNDEBUG`, AppleClang 21.0.0) into
  `<tree>/build-bench` for each tree, so neither run used a stale `.so`.
  No C++ was patched for measurement, so there was no scratch build to
  sanitize.
- **Inputs** (bench.py defaults):
  - DEM `tests/fixtures/dem_archive/7908_3_10m_z33.tif`, tolerance 1 m.
  - Domains: `tile` (the whole DEM) and
    `docs/benchmarks/2026-09-26/quarter.geojson`.
  - 5 repeats per cell. Threads 0 (the CLI default) and 1 to 20.
- **Machine:** Apple M1 Max (8 P + 2 E), 32 GiB, macOS 27.0, Python 3.14.7,
  numpy 2.5.3. Every block ran under `caffeinate -ims`.
- **Power: battery, discharging, in every run**, falling from 82 % (21:46) to
  71 % (22:11). Each `run.json` has the `pmset -g batt` text from before and
  after the run, and the logs have one reading before each run.
- **Dirty flag.** The branch runs record `dirty: true`. The only change in the
  tree was this untracked evidence directory (`git status --porcelain` showed
  `?? docs/benchmarks/2026-10-01/15c-1-acceptance/`). The base runs are clean.

## Why the base was re-measured

bench.py's baseline search found no comparable stored run. Every run printed
`NO BASELINE: no comparable stored run`, and the brief says to measure
`a130f7c` with `--tree`, back to back, in that case. The search failed for two
reasons. Both are bench.py defects; they are reported here and not fixed:

1. **`comparable()` compares the domain list in order.** The stored battery
   runs (`2026-09-27/21b-r2`, `2026-09-28/16b12-acceptance/*`) record the
   domains as `[quarter, tile]`. Today's default records them as
   `[tile, quarter]`. The DEM, the hashes and the tolerance are the same, but
   the runs are judged not comparable: `domains: [('tile', None),
   ('quarter', 'b7f7…')] vs [('quarter', 'b7f7…'), ('tile', None)]`. This was
   checked by calling `bench.comparable` on the loaded records.
2. **`find_baseline()` globs only `<out-root>/*/*/run.json`.** Acceptance runs
   stored one level deeper, such as `16b0-acceptance/*`,
   `16b12-acceptance/*` and this directory, are never found.

## Method

- **Five pairs, about 2.3 min per run.**
  - Pairs 1-3 ran base first (`scripts/pairs.sh`). Each branch run was judged
    by `run --baseline` against the base run just before it.
  - Pairs 4-5 ran branch first (`scripts/pairs_reversed.sh`), so that a drift
    over the session (battery level, heat) does not always fall on the
    branch. Each branch run was then judged with
    `bench.py compare … --baseline` against the base run that followed it.
- **Pooled figure.** `scripts/aggregate.py` takes, for each (domain,
  threads) cell, the median of the five per-run medians, for each build. Its
  output is `logs/aggregate.md`.
- **Noise floor.** The same script also compares the same build with itself,
  across its five runs.

## Results

### Mesh and quality: identical

These figures are the same in all 10 runs, base and branch:

| domain | worst angle | max degree | within tolerance | Delaunay checked / ambiguous / violations | mesh sha256 |
|---|---:|---:|---|---|---|
| quarter | 0.3955° | 18 | yes | 641 791 / 80 136 / 0 | `ccebf96a86c6c5e2…` |
| tile | 0.6296° | 74 | yes | 692 056 / 94 469 / 0 | `11741a81adfa17b3…` |

Both mesh hashes also equal those stored by 16b-1/2's acceptance
(`2026-09-28/16b12-acceptance/16b12-r4/run.json`).

### Refine time (s), pooled medians of 5 runs per build

`change` is the branch over the base. A threads value of 0 is the CLI default.

| threads | quarter base | quarter 15c-1 | change | tile base | tile 15c-1 | change |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.1702 | 0.1696 | -0.3 % | 0.1884 | 0.1952 | +3.6 % |
| 1 | 0.4723 | 0.4597 | -2.7 % | 0.5102 | 0.5174 | +1.4 % |
| 2 | 0.2959 | 0.2946 | -0.4 % | 0.3363 | 0.3343 | -0.6 % |
| 3 | 0.2382 | 0.2364 | -0.8 % | 0.2722 | 0.2720 | -0.1 % |
| 4 | 0.2099 | 0.2077 | -1.0 % | 0.2363 | 0.2385 | +0.9 % |
| 5 | 0.1914 | 0.1926 | +0.7 % | 0.2153 | 0.2187 | +1.6 % |
| 6 | 0.1811 | 0.1803 | -0.5 % | 0.2054 | 0.2077 | +1.1 % |
| 7 | 0.1740 | 0.1721 | -1.1 % | 0.1974 | 0.1969 | -0.3 % |
| 8 | 0.1691 | 0.1712 | +1.2 % | 0.1925 | 0.1921 | -0.2 % |
| 9 | 0.1714 | 0.1684 | -1.8 % | 0.1925 | 0.1934 | +0.5 % |
| 10 | 0.1692 | 0.1725 | +2.0 % | 0.1927 | 0.1917 | -0.5 % |
| 11 | 0.1706 | 0.1697 | -0.6 % | 0.1930 | 0.1928 | -0.1 % |
| 12 | 0.1685 | 0.1708 | +1.4 % | 0.1891 | 0.1953 | +3.3 % |
| 13 | 0.1706 | 0.1701 | -0.3 % | 0.1905 | 0.1909 | +0.2 % |
| 14 | 0.1688 | 0.1695 | +0.5 % | 0.1934 | 0.1923 | -0.6 % |
| 15 | 0.1688 | 0.1706 | +1.0 % | 0.1922 | 0.1899 | -1.2 % |
| 16 | 0.1692 | 0.1678 | -0.9 % | 0.1938 | 0.1904 | -1.8 % |
| 17 | 0.1699 | 0.1688 | -0.6 % | 0.1953 | 0.1938 | -0.8 % |
| 18 | 0.1709 | 0.1698 | -0.7 % | 0.1945 | 0.1937 | -0.5 % |
| 19 | 0.1722 | 0.1666 | -3.2 % | 0.1949 | 0.1905 | -2.3 % |
| 20 | 0.1716 | 0.1692 | -1.4 % | 0.1927 | 0.1949 | +1.1 % |

- **All 42 cells:** the median change is -0.39 %, and the range is -3.2 % to
  +3.6 %. No cell is beyond ±5 %. The min-max per cell is in
  `logs/aggregate.md`.
- **Noise floor:** the same build compared with itself, across its own runs,
  moves -13.1 % to +13.1 % per cell (840 cells). 154 of those cells are
  beyond ±5 %.

### Per pair, and bench.py's own verdicts

| pair | order | quarter (median over 21 cells) | tile | bench.py verdict |
|---|---|---:|---:|---|
| 1 | base first | +0.0 % | +1.0 % | REGRESSION on 3 quarter cells (t=2 +5.2 %, t=10 +13.4 %, t=12 +8.3 %) |
| 2 | base first | -2.3 % | +1.0 % | REGRESSION on 4 tile cells (t=1 +14.7 %, t=3 +9.4 %, t=5 +7.4 %, t=12 +5.3 %) |
| 3 | base first | -0.8 % | +3.9 % | REGRESSION on 5 tile cells (t=2, 10, 15, 17, 20; +5.8 % to +10.4 %) |
| 4 | branch first | -0.8 % | -2.7 % | ACCEPTED (`compare`) |
| 5 | branch first | +0.7 % | +1.4 % | ACCEPTED (`compare`) |

- **The flags in pairs 1-3 are read as noise.** Each pair flags different
  cells, and none of them recurs. Pooled over five pairs, none is above
  +3.6 %. The same build moves up to ±13 % between its own runs.
- **Order effect.** In pairs 1-3 the branch always ran second, and the tile's
  median change was +1.0 to +3.9 %. With the branch first (pairs 4-5) it was
  -2.7 % and +1.4 %. This fits a drift over the session landing on the
  second run, but that is an inference: no profile was taken.
- **Not resolved below about 2 %.** The pairs cannot tell a change smaller
  than about 2 % from noise. bench.py's threshold is 5 % on the median.

### Ceiling (speed-up over 1 thread; pooled medians)

| domain | build | 1 → 20 threads | best | at threads |
|---|---|---:|---:|---:|
| quarter | base a130f7c | 2.75x | 2.80x | 12 |
| quarter | 15c-1 | 2.72x | 2.76x | 19 |
| tile | base a130f7c | 2.65x | 2.70x | 12 |
| tile | 15c-1 | 2.65x | 2.72x | 15 |

The 2026-09-26 reference is 2.2x. The ceiling is unchanged.

## Files

- `base-a130f7c-r{1..5}/`, `15c-1-r{1..5}/`: bench.py's evidence
  (`run.json`, `raw.tsv`, a generated `README.md`).
  - bench.py resolves every path, so those files record absolute paths
    (`/Users/…`). bench.py has no option to record relative ones.
  - The verdicts stored in pairs 4-5's branch runs read `NO BASELINE`,
    because they ran first. Their verdicts against their pairs are in
    `logs/pairs_reversed.log`.
- `logs/pairs.log`, `logs/pairs_reversed.log`: the console output, with a
  `pmset` reading before each run. `logs/aggregate.md`: the tables above.
- `scripts/`: the scripts that made the runs and the tables.
- **Meshes** are not committed. The quality meshes (`*_quarter.vtk`,
  `*_tile.vtk`, ASCII) were written to the job's scratch directory, which is
  lost with the session. Rerunning the commands below regenerates them.

## Reproduce

From the 15c-1 worktree at `6ab7ad5`, with a venv that has numpy and typer
and in which `tin_engine.stats` imports:

```sh
git worktree add --detach "$SCRATCH/base-a130f7c" a130f7c
PY=<venv>/bin/python BASE="$SCRATCH/base-a130f7c" MESH="$SCRATCH/meshes" \
  caffeinate -ims sh docs/benchmarks/2026-10-01/15c-1-acceptance/scripts/pairs.sh 1 2 3
PY=<venv>/bin/python BASE="$SCRATCH/base-a130f7c" MESH="$SCRATCH/meshes" \
  caffeinate -ims sh docs/benchmarks/2026-10-01/15c-1-acceptance/scripts/pairs_reversed.sh 4 5
python3 docs/benchmarks/2026-10-01/15c-1-acceptance/scripts/aggregate.py
```

The scripts write under `docs/benchmarks/<today>/15c-1-acceptance/`. Record
`pmset -g batt` before you start: a run on AC power is compared only with an
AC baseline.
