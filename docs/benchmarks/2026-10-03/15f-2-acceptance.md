# 15f-2 acceptance (@perf, 2026-10-03)

Branch `worktree-15f-2` at `5ac2903` (15f-2: the edge strip in the
refinement loop, C++ only; `strip_scan.hpp` new, `refine_points.hpp`
reworked into `detail::point_loop`, `scan.hpp` given a coincidence radius
defaulted to 0; nothing in the CLI calls the strip yet, L15), against master
`29d713b`, master's tip at the time of the run. The PR's range is
`9815ea4..5ac2903`; `29d713b` differs from `9815ea4` only in tests and docs
(`git diff --stat 9815ea4 29d713b`: four docs, three Python test files), so
it is the same code base. The design's acceptance for 15f-2
(`docs/increments/15f-edge-strip.md`, "Acceptance (`@perf`)"): `tools/bench.py
run` back to back with `--tree` against the previous merge, geometry
identical, refine time within noise, power state recorded.

Apple M1 Max (8 P + 2 E, 32 GiB), macOS 27.0, Python 3.14.7, numpy 2.5.3.
**AC power** for all four runs, battery 100 % and charged, `pmset -g batt`
before and after each run (in each `run.json`; `bench.py` records `ac` only
when both agree). Ola was present; no other agent ran and nothing else built
during the runs (14:35-14:40 UTC).

## Verdict

**ACCEPTED.** The `_core` binaries differ (as they must: `scan`'s signature
changed), but the meshes are byte-identical, quality is identical, and refine
times are within noise on both domains at every thread count.

## Method

`tools/bench.py run` (blob `6c5de4c`, unchanged between `29d713b` and
`5ac2903`), Release `_core` built by `bench.py` into each tree's
`build-bench/`: default DEM `tests/fixtures/dem_archive/7908_3_10m_z33.tif`,
tolerance 1, domains `tile` and `quarter`, threads 1, 2, 4, 6, 8, 10, 16, 20
and the CLI default (10), 5 repeats interleaved over thread counts. Four runs,
alternating, in the order A (master), B (15f-2), B2 (15f-2), A2 (master).
Driver: `15f-2-acceptance/pairs.sh` (its paths are this session's
scratchpad), log `15f-2-acceptance/pairs.log`. Evidence in `bench.py`'s
format: `15f-2-base-29d713b/`, `15f-2/`, `15f-2-r2/`,
`15f-2-base-29d713b-r2/`. Tables from `15f-2-acceptance/summarize.py`
(output in `15f-2-acceptance/tables.md`).

`bench.py`'s own verdicts: B and B2 ACCEPTED against A; A2 ACCEPTED against
B2. A, judged against the newest comparable stored run, printed
REGRESSION at 1 and 2 threads (+7.1 % to +8.3 %). `run.json` does not name
the run it picked, but its baseline figures (0.4030, 0.2612, 0.4545,
0.2978) are exactly those of `15f-1-base-390b516-r2` (started 12:34 UTC
today, binary `41a22a6`). That is master against
master across sessions, not this branch: A2, the same `29d713b` binary run
after the branch, is as slow as A (table below), and 15e's own runs of the
same binary (`a6309dc`, `15e/` and `15e-r2/`, 10:24-10:29 UTC) measured
tile 1 thread 0.4822 and 0.4714 s against 0.4880 and 0.4896 s here. The
cause of the session-to-session drift is not measured. The branch's runs
are marked dirty: the untracked evidence directories (`git status` showed
nothing else), no source change.

## Builds and geometry

| run | commit | power | `_core` sha256 | tile proc s (default) | quarter proc s | tile mesh sha256 | quarter mesh sha256 |
|---|---|---|---|---:|---:|---|---|
| 15f-2-base-29d713b (A) | 29d713b | ac | a6309dc2b621 | 0.610 | 0.542 | 11741a81adfa | ccebf96a86c6 |
| 15f-2 (B) | 5ac2903 | ac | 4b5f875da201 | 0.605 | 0.543 | 11741a81adfa | ccebf96a86c6 |
| 15f-2-r2 (B2) | 5ac2903 | ac | 4b5f875da201 | 0.611 | 0.549 | 11741a81adfa | ccebf96a86c6 |
| 15f-2-base-29d713b-r2 (A2) | 29d713b | ac | a6309dc2b621 | 0.604 | 0.552 | 11741a81adfa | ccebf96a86c6 |

- **`_core` changes.** Full sha256: master
  `a6309dc2b621a611b4426fdad59094f451d70fadef6533a062e72fd57edd33b7`,
  branch `4b5f875da20112c9c18efc1602ef387eafc898c51af6361f4cec7f9ea2fc2709`
  (in the `run.json` files and rechecked with `shasum -a 256`). `_core` has
  one translation unit, `bindings/core.cpp`, compiled to LTO bitcode; the
  shipped `.so` keeps no local symbols, so the comparison was made on the
  bitcode (`llvm-dis` of both `core.cpp.o`, function bodies compared with
  type, value and metadata numbering normalised).
  - `refine<RasterView<double|float>>` itself: identical after
    normalisation.
  - Its per-thread block (`parallel_util::for_each_block<... refine ...>`)
    differs in one call: `scan(dem, m, t)` became
    `scan(dem, m, t, 0.0)`, a call to the new four-argument `scan`
    (`...LatticeMeshEjd`). `scan` is not inlined into it in the bitcode, so
    the radius arrives as a runtime 0.0 unless link-time optimisation
    propagates it (not checked on the stripped `.so`).
  - `scan`'s row-span body grew from 484 to 554 IR lines (double DEM): the
    `skip` set setup (skipped at once when `radius == 0`) and the
    `skipped(p)` test, `n_skip != 0 && ...`, in the void branch and the
    off-node branch. The all-node branch, the one with the incremental
    orientation steps, is unchanged in source.
  - `refine_points<CheckPoints>` is now a wrapper over
    `detail::point_loop<CheckPoints, NoSet>`, and `scan_strip`,
    `strip_fits`, `near_constraint` and `cut` are new. `bench.py` does not
    time `refine_points`: its runs have no feature points.
- **The meshes are identical**: `bench.py`'s mesh sha256 (from `POINTS` on)
  agrees in all four runs per domain, and the four ASCII meshes per domain
  were compared byte for byte with `cmp` against A's: identical. The probe
  can tell meshes apart: tile and quarter differ (`cmp`, line 5). The hashes
  equal 15f-1's and 15e's acceptance runs.
- Every sample of a domain, at every thread count, in every run, reports the
  same `max_error`, `rounds`, `inserted` and `flips` (`summarize.py` asserts
  it).

## Quality (identical in all four runs)

| domain | worst angle deg | share < 1 deg | max degree | within tolerance | Delaunay violations | max_error | rounds | inserted | flips |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| tile | 0.6296 | 2.15e-06 | 74 | yes | 0 of 692,056 | 0.9999759 | 53 | 219,837 | 445,675 |
| quarter | 0.3955 | 2.80e-05 | 18 | yes | 0 of 641,791 | 0.9999797 | 41 | 213,464 | 445,657 |

## Refine time (median of 5, seconds)

| domain | threads | A base | B 15f-2 | B2 15f-2 | A2 base | B/A | B2/A2 | (B+B2)/(A+A2) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| quarter | default | 0.1708 | 0.1717 | 0.1689 | 0.1718 | 1.005 | 0.983 | 0.994 |
| quarter | 1 | 0.4329 | 0.4372 | 0.4406 | 0.4407 | 1.010 | 1.000 | 1.005 |
| quarter | 2 | 0.2830 | 0.2820 | 0.2776 | 0.2901 | 0.997 | 0.957 | 0.976 |
| quarter | 4 | 0.2051 | 0.2053 | 0.2037 | 0.2092 | 1.001 | 0.974 | 0.987 |
| quarter | 6 | 0.1784 | 0.1774 | 0.1780 | 0.1790 | 0.995 | 0.994 | 0.994 |
| quarter | 8 | 0.1677 | 0.1677 | 0.1678 | 0.1696 | 1.000 | 0.989 | 0.995 |
| quarter | 10 | 0.1682 | 0.1673 | 0.1677 | 0.1686 | 0.994 | 0.995 | 0.995 |
| quarter | 16 | 0.1698 | 0.1693 | 0.1704 | 0.1686 | 0.997 | 1.011 | 1.004 |
| quarter | 20 | 0.1690 | 0.1711 | 0.1692 | 0.1682 | 1.012 | 1.006 | 1.009 |
| tile | default | 0.1883 | 0.1894 | 0.1897 | 0.1896 | 1.006 | 1.001 | 1.003 |
| tile | 1 | 0.4880 | 0.4897 | 0.4947 | 0.4896 | 1.003 | 1.010 | 1.007 |
| tile | 2 | 0.3190 | 0.3172 | 0.3271 | 0.3155 | 0.995 | 1.037 | 1.015 |
| tile | 4 | 0.2312 | 0.2334 | 0.2335 | 0.2327 | 1.010 | 1.004 | 1.007 |
| tile | 6 | 0.2025 | 0.2030 | 0.2020 | 0.2033 | 1.003 | 0.994 | 0.998 |
| tile | 8 | 0.1883 | 0.1880 | 0.1884 | 0.1899 | 0.998 | 0.992 | 0.995 |
| tile | 10 | 0.1876 | 0.1885 | 0.1887 | 0.1879 | 1.004 | 1.004 | 1.004 |
| tile | 16 | 0.1885 | 0.1882 | 0.1915 | 0.1904 | 0.998 | 1.006 | 1.002 |
| tile | 20 | 0.1918 | 0.1886 | 0.1898 | 0.1908 | 0.983 | 0.995 | 0.989 |

Every per-pair ratio is within 0.957-1.037 and every pooled ratio within
0.976-1.015, inside `bench.py`'s 5 % threshold, and the deviations go both
ways (the largest, at 2 threads, are +3.7 % on the tile and -4.3 % on the
quarter in the same pair). At 1 thread, where a per-node cost in `scan`
would show undiluted, the pooled ratios are 1.007 (tile) and 1.005
(quarter). **Refine did not move beyond noise.** This does not separate a
sub-1 % cost of the extra `skipped(p)` test from noise, and how much of
`scan`'s time falls in the off-node and void branches on these domains was
not profiled. Ceiling (B, `15f-2/README.md`): tile 2.60x at 20 threads over
1, best 2.61x at 8; quarter 2.55x at 20, best 2.61x at 10; against the
2026-09-26 reference of about 2.2x.

## Not done

No micro-timing of `refine_strip` / `refine_points` with a strip (optional
in the brief). It would need a C++ driver, and so a sanitized run before
timing; nothing calls the strip from the CLI yet, and 15f-3's acceptance
reports the strip's own phases where they run. No sanitized scratch build
was needed: no simulation or instrumentation patch was made.

## Meshes and builds

The quality meshes (about 22 MB each, ASCII VTK) were written to this
session's scratchpad and deleted after the comparison; the master worktree
(`git worktree add --detach <dir> 29d713b`) and both `build-bench/`
directories were removed. To regenerate, from the branch's worktree with the
project venv's python as `PY` and `B` a detached worktree of `29d713b`:

```bash
git worktree add --detach $B 29d713b
for T in $B $PWD; do $PY -c "import sys; sys.path.insert(0,'tools'); import bench; \
  from pathlib import Path; print(bench.build(bench.make_runner(), Path('$T')))"; done
TH=1,2,4,6,8,10,16,20; D=docs/benchmarks/2026-10-03
PYTHONPATH=$B/build-bench/pkg $PY tools/bench.py run --label 15f-2-base-29d713b --tree $B --threads $TH --repeats 5
PYTHONPATH=$PWD/build-bench/pkg $PY tools/bench.py run --label 15f-2 --threads $TH --repeats 5 --baseline $D/15f-2-base-29d713b
PYTHONPATH=$PWD/build-bench/pkg $PY tools/bench.py run --label 15f-2-r2 --threads $TH --repeats 5 --baseline $D/15f-2-base-29d713b
PYTHONPATH=$B/build-bench/pkg $PY tools/bench.py run --label 15f-2-base-29d713b-r2 --tree $B --threads $TH --repeats 5 --baseline $D/15f-2-r2
$PY $D/15f-2-acceptance/summarize.py
```

The bitcode comparison: `llvm-dis` (Homebrew LLVM) on each tree's
`build-bench/CMakeFiles/_core.dir/bindings/core.cpp.o`, then compare
`define` bodies with `%name.N`, `!N` and `#N` numbering normalised.
