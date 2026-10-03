# 15f-1 acceptance (@perf, 2026-10-03)

Branch `worktree-agent-a0cb49bea16623b07` at `fe06447` (15f-1: the edge
strip's constraint check-point generator and store, `constraint_points.hpp`,
and `detail::lattice_position` extracted from `refine.hpp`'s `to_lattice`),
against master `390b516`, the previous increment's merge. The design's
acceptance for 15f-1 (`docs/increments/15f-edge-strip.md`, "Acceptance
(`@perf`)"): `tools/bench.py` on the quarter and the tile with the thread
sweep, back to back with `--tree` against master; geometry identical, refine
time within noise, power state recorded.

Apple M1 Max (8 P + 2 E, 32 GiB), macOS 27.0, Python 3.14.7, numpy 2.5.3.
**AC power** for all four runs, battery 100 % and charged, `pmset -g batt`
before and after each run (in each `run.json`; `bench.py` records `ac` only
when both agree). Ola was away; no other agent ran during the runs
(12:31-12:35 UTC), and `caffeinate` held the machine awake.

## Verdict

**ACCEPTED.** The two builds' `_core` are byte-identical, so are the
meshes, and refine times are within noise on both domains at every thread
count.

## Method

`tools/bench.py run` (blob `6c5de4c`, unchanged between `390b516` and
`fe06447`), Release `_core` built by `bench.py` into each tree's
`build-bench/`: default DEM `tests/fixtures/dem_archive/7908_3_10m_z33.tif`,
tolerance 1, domains `tile` and `quarter`, threads 1, 2, 4, 6, 8, 10, 16, 20
and the CLI default (10), 5 repeats interleaved over thread counts. Four runs,
alternating, in the order A (master), B (15f-1), B2 (15f-1), A2 (master), so
that order drift lands on both sides. Driver: `15f-1-acceptance/pairs.sh`
(its paths are this session's scratchpad), log `15f-1-acceptance/pairs.log`.
Evidence in `bench.py`'s format, where `find_baseline` looks:
`15f-1-base-390b516/`, `15f-1/`, `15f-1-r2/`, `15f-1-base-390b516-r2/`.
Tables below from `15f-1-acceptance/summarize.py` (output also in
`15f-1-acceptance/tables.md`).

`bench.py`'s own verdicts: A ACCEPTED against the newest comparable stored
run, `2026-10-02/15c-2-base-e1f6042-r2`; B and B2 ACCEPTED against A; A2
ACCEPTED against B2. The branch's runs are marked dirty: the untracked
evidence directory `docs/benchmarks/2026-10-03/` (`git status` showed nothing
else), no source change.

## Builds and geometry

| run | commit | power | `_core` sha256 | tile proc s (default) | quarter proc s | tile mesh sha256 | quarter mesh sha256 |
|---|---|---|---|---:|---:|---|---|
| 15f-1-base-390b516 (A) | 390b516 | ac | 41a22a6168db | 0.590 | 0.518 | 11741a81adfa | ccebf96a86c6 |
| 15f-1 (B) | fe06447 | ac | 41a22a6168db | 0.583 | 0.520 | 11741a81adfa | ccebf96a86c6 |
| 15f-1-r2 (B2) | fe06447 | ac | 41a22a6168db | 0.581 | 0.516 | 11741a81adfa | ccebf96a86c6 |
| 15f-1-base-390b516-r2 (A2) | 390b516 | ac | 41a22a6168db | 0.578 | 0.520 | 11741a81adfa | ccebf96a86c6 |

- **`_core` is byte-identical** on master and the branch (full sha256
  `41a22a6168dbaa141c18ce4f867e011a13176aec93769a75f6be6c95250485f3`, in all
  four `run.json` and rechecked with `shasum -a 256` on both `build-bench/pkg`
  copies). The identity says something: `bindings/core.cpp` includes
  `refine.hpp` (`build-bench/CMakeFiles/_core.dir/bindings/core.cpp.o.d`
  lists it), so the extraction of `lattice_position` was compiled and inlined
  to the same code. `constraint_points.hpp` is in no `_core` translation unit
  yet (nothing calls it; only its tests include it).
- **The meshes are identical**: `bench.py`'s mesh sha256 (from `POINTS` on)
  agrees in all four runs per domain, and the four ASCII meshes per domain were
  compared byte for byte with `cmp` against A's: identical. The probe can tell
  meshes apart: tile and quarter differ (`cmp`, line 5). The hashes also equal
  the 15c-2 acceptance's (`../2026-10-02/15c-2-acceptance/README.md`).
- Every sample of a domain, at every thread count, in every run, reports the
  same `max_error`, `rounds`, `inserted` and `flips` (`summarize.py` asserts
  it).

## Quality (identical in all four runs)

| domain | worst angle deg | share < 1 deg | max degree | within tolerance | Delaunay violations | max_error | rounds | inserted | flips |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| tile | 0.6296 | 2.15e-06 | 74 | yes | 0 of 692,056 | 0.9999759 | 53 | 219,837 | 445,675 |
| quarter | 0.3955 | 2.80e-05 | 18 | yes | 0 of 641,791 | 0.9999797 | 41 | 213,464 | 445,657 |

## Refine time (median of 5, seconds)

| domain | threads | A base | B 15f-1 | B2 15f-1 | A2 base | B/A | B2/A2 | (B+B2)/(A+A2) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| quarter | default | 0.1628 | 0.1659 | 0.1632 | 0.1651 | 1.019 | 0.988 | 1.004 |
| quarter | 1 | 0.4001 | 0.4002 | 0.4005 | 0.4030 | 1.000 | 0.994 | 0.997 |
| quarter | 2 | 0.2610 | 0.2625 | 0.2613 | 0.2612 | 1.006 | 1.000 | 1.003 |
| quarter | 4 | 0.1971 | 0.1969 | 0.1984 | 0.1962 | 0.999 | 1.011 | 1.005 |
| quarter | 6 | 0.1727 | 0.1732 | 0.1732 | 0.1729 | 1.003 | 1.001 | 1.002 |
| quarter | 8 | 0.1608 | 0.1620 | 0.1616 | 0.1619 | 1.008 | 0.998 | 1.003 |
| quarter | 10 | 0.1648 | 0.1618 | 0.1636 | 0.1655 | 0.982 | 0.989 | 0.985 |
| quarter | 16 | 0.1638 | 0.1657 | 0.1632 | 0.1664 | 1.012 | 0.981 | 0.996 |
| quarter | 20 | 0.1676 | 0.1662 | 0.1669 | 0.1639 | 0.992 | 1.018 | 1.005 |
| tile | default | 0.1859 | 0.1864 | 0.1866 | 0.1844 | 1.003 | 1.012 | 1.007 |
| tile | 1 | 0.4581 | 0.4565 | 0.4556 | 0.4545 | 0.997 | 1.002 | 0.999 |
| tile | 2 | 0.2998 | 0.2985 | 0.2966 | 0.2978 | 0.996 | 0.996 | 0.996 |
| tile | 4 | 0.2241 | 0.2213 | 0.2248 | 0.2235 | 0.987 | 1.006 | 0.997 |
| tile | 6 | 0.1983 | 0.1970 | 0.1972 | 0.1970 | 0.993 | 1.001 | 0.997 |
| tile | 8 | 0.1832 | 0.1838 | 0.1830 | 0.1829 | 1.003 | 1.001 | 1.002 |
| tile | 10 | 0.1859 | 0.1855 | 0.1857 | 0.1856 | 0.998 | 1.001 | 0.999 |
| tile | 16 | 0.1876 | 0.1860 | 0.1813 | 0.1865 | 0.992 | 0.972 | 0.982 |
| tile | 20 | 0.1885 | 0.1879 | 0.1864 | 0.1871 | 0.997 | 0.996 | 0.997 |

Every per-pair ratio is within 0.97-1.02 and every pooled ratio within
0.982-1.007, well inside `bench.py`'s 5 % threshold. With the same binary on
both sides this is the noise floor of the method on this machine, not a
measurement of the branch. Ceiling (B, `15f-1/README.md`): tile 2.43x at 20
threads over 1, best 2.48x at 8; quarter 2.41x at 20, best 2.47x at 10;
against the 2026-09-26 reference of about 2.2x.

## Not done

No micro-timing of `constraint_check_points`: the design's acceptance does
not ask for one for 15f-1, and nothing calls the generator yet. Its time is
measured where it runs, in 15f-2's acceptance ("The strip's own phases are
reported beside refine"). No sanitized scratch build was needed: no
simulation or instrumentation patch was made.

## Meshes and builds

The quality meshes (about 22 MB each, ASCII VTK) were written to this
session's scratchpad and deleted after the comparison; the master worktree
(`git worktree add --detach <dir> 390b516`) and both `build-bench/`
directories were removed. To regenerate, from the branch's worktree with the
project venv's python as `PY` and `B` a detached worktree of `390b516`:

```bash
git worktree add --detach $B 390b516
for T in $B $PWD; do $PY -c "import sys; sys.path.insert(0,'tools'); import bench; \
  from pathlib import Path; print(bench.build(bench.make_runner(), Path('$T')))"; done
TH=1,2,4,6,8,10,16,20; D=docs/benchmarks/2026-10-03
PYTHONPATH=$B/build-bench/pkg $PY tools/bench.py run --label 15f-1-base-390b516 --tree $B --threads $TH --repeats 5
PYTHONPATH=$PWD/build-bench/pkg $PY tools/bench.py run --label 15f-1 --threads $TH --repeats 5 --baseline $D/15f-1-base-390b516
PYTHONPATH=$PWD/build-bench/pkg $PY tools/bench.py run --label 15f-1-r2 --threads $TH --repeats 5 --baseline $D/15f-1-base-390b516
PYTHONPATH=$B/build-bench/pkg $PY tools/bench.py run --label 15f-1-base-390b516-r2 --tree $B --threads $TH --repeats 5 --baseline $D/15f-1-r2
$PY $D/15f-1-acceptance/summarize.py
```

Add `--mesh-dir DIR` to keep the quality meshes somewhere known; each run
takes about a minute.
