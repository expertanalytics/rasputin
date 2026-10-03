# 15c-2 acceptance (@perf, 2026-10-02)

Branch `worktree-15c-2` at `d6ac78d`, against master's merge `e1f6042` (PR 133,
15c-1). Apple M1 Max (8 P + 2 E, 32 GiB), macOS 27.0, Python 3.14.7, **AC
power** for every run (`pmset -g batt` before and after each, in the JSON).
Release `_core` built by `tools/bench.py`; the branch's and master's builds
are **byte-identical** (sha256 `7dee5ad918d0…` in all four `run.json`), as
15c-2 changes no C++.

## Verdict

**REGRESSION**: `rasputin mesh --dem anadem-v1 --out-crs …` **crashes at
random**: 7 of 47 process starts in the sweep below (15 %) died on a signal,
and 2 of 8 in the first probe (`runs/first_failures/`): one SIGBUS, one
`GEOSException: vector`. All 8 macOS crash reports of the session (SIGSEGV 5,
SIGBUS 2, SIGABRT 1; `runs/crash_reports_faulting_threads.txt`) fault in a
Python worker thread inside GEOS's prepared-polygon `intersects` machinery
(the lazy `STRtree` build, or `MonotoneChain` overlaps). The one call on that
path from a worker thread is `check_point_blocks`'s `inside.intersects(...)`
(`src_python/tin_engine/target_grid.py`, the `block` closure run by `_pool` on
a `ThreadPoolExecutor`), where `inside = prep(domain.polygon)` is one prepared
geometry shared by all workers. That concurrent use of it is the cause is
**thought to be** so: it fits every stack, and no single-threaded or
pre-warmed run has been measured.

Everything else measured passes:

- **Norway's path (J1)**: mesh hashes unchanged on both domains, refine and
  process time within noise (below).
- **The geographic path**: the independent final check finds **0 source nodes
  over tolerance, interior and strip**, on every mesh it was run on (table
  below; 50-5 m, both CRSs), and its control fails when the mesh is moved
  15 m. At 1 m it was not run to the end (Not done).

## (a) Norway's path: `tools/bench.py`, 1 m benchmark and thread sweep

`bench.py run` found no comparable stored run newer than 2026-09-27 (refine
has changed since), so master `e1f6042` was measured with `--tree` back to
back, in the order A (master), B (15c-2), B2 (15c-2), A2 (master). Threads
1, 2, 4, 6, 8, 10, 16, 20 and the default (10), 5 repeats, interleaved;
default DEM and both domains (`tile`, `quarter`). Evidence (in
`bench.py`'s format, where `find_baseline` looks): `../15c-2-base-e1f6042/`,
`../15c-2/`, `../15c-2-r2/`, `../15c-2-base-e1f6042-r2/`.

Median refine seconds:

| domain | threads | A base | B 15c-2 | B2 15c-2 | A2 base | B/A | B2/A2 | (B+B2)/(A+A2) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| quarter | default | 0.1666 | 0.1671 | 0.1620 | 0.1655 | 1.003 | 0.979 | 0.991 |
| quarter | 1 | 0.4019 | 0.4338 | 0.4338 | 0.4270 | 1.079 | 1.016 | 1.047 |
| quarter | 2 | 0.2626 | 0.2769 | 0.2829 | 0.2856 | 1.054 | 0.991 | 1.021 |
| quarter | 4 | 0.1975 | 0.2030 | 0.2016 | 0.2039 | 1.028 | 0.989 | 1.008 |
| quarter | 6 | 0.1729 | 0.1776 | 0.1772 | 0.1820 | 1.027 | 0.974 | 1.000 |
| quarter | 8 | 0.1606 | 0.1653 | 0.1663 | 0.1709 | 1.029 | 0.973 | 1.000 |
| quarter | 10 | 0.1665 | 0.1668 | 0.1647 | 0.1632 | 1.002 | 1.009 | 1.006 |
| quarter | 16 | 0.1662 | 0.1653 | 0.1670 | 0.1687 | 0.994 | 0.990 | 0.992 |
| quarter | 20 | 0.1707 | 0.1682 | 0.1646 | 0.1667 | 0.986 | 0.987 | 0.986 |
| tile | default | 0.1858 | 0.1880 | 0.1869 | 0.1993 | 1.012 | 0.938 | 0.974 |
| tile | 1 | 0.4595 | 0.4824 | 0.5155 | 0.5060 | 1.050 | 1.019 | 1.034 |
| tile | 2 | 0.2985 | 0.3135 | 0.3194 | 0.3298 | 1.050 | 0.968 | 1.007 |
| tile | 4 | 0.2249 | 0.2365 | 0.2320 | 0.2378 | 1.052 | 0.976 | 1.013 |
| tile | 6 | 0.1977 | 0.2057 | 0.2008 | 0.2047 | 1.040 | 0.981 | 1.010 |
| tile | 8 | 0.1834 | 0.1881 | 0.1901 | 0.1928 | 1.025 | 0.986 | 1.005 |
| tile | 10 | 0.1868 | 0.1887 | 0.1923 | 0.1853 | 1.010 | 1.038 | 1.024 |
| tile | 16 | 0.1865 | 0.1892 | 0.1989 | 0.1885 | 1.014 | 1.055 | 1.035 |
| tile | 20 | 0.1876 | 0.1912 | 0.1874 | 0.1881 | 1.019 | 0.996 | 1.008 |

| run | power | .so sha256 | tile proc_s (default threads) | quarter proc_s | tile mesh sha256 | quarter mesh sha256 |
|---|---|---|---:|---:|---|---|
| 15c-2-base-e1f6042 | ac | 7dee5ad918d0 | 0.538 | 0.470 | 11741a81adfa | ccebf96a86c6 |
| 15c-2 | ac | 7dee5ad918d0 | 0.556 | 0.484 | 11741a81adfa | ccebf96a86c6 |
| 15c-2-r2 | ac | 7dee5ad918d0 | 0.552 | 0.471 | 11741a81adfa | ccebf96a86c6 |
| 15c-2-base-e1f6042-r2 | ac | 7dee5ad918d0 | 0.588 | 0.479 | 11741a81adfa | ccebf96a86c6 |

`bench.py`'s own verdicts read REGRESSION on a few low thread counts (+5 to
+12 %) for B against A, and also for A2 against B: the same binary, run later,
is slower. The first run of the four (A) was the fastest at 1 to 8 threads on both
domains. With the identical `.so` and the pooled ratio within 0.97-1.05, the
differences are run-order drift, not the branch. Quality equal in all four:
tile worst angle 0.6296°, max degree 74; quarter 0.3955°, 18; 0 Delaunay
violations; within tolerance. `bench.py` marks the branch's runs dirty: the
untracked evidence directory, no source change (the `.so` hashes agree).

## (b) The geographic path: the Velhas piece on ANADEM

`rasputin mesh --dem anadem-v1 --cache ../rasputin_data/cache --domain
…/bho2017_5k_76949_outline_epsg4674.geojson --out-crs X --tolerance T`
(BHO 76949, 11,667.6 km²). The catalogue key is `anadem-v1`; `--dem anadem`
is read as a file path and refused. Without `--out-crs` the CLI refuses and
suggests (`runs/refusal_no_out_crs.out`):

```
+proj=tmerc +lat_0=0 +lon_0=-44.1 +k=0.999972 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs
```

worst scale error 0.00282 %. Both that and EPSG:31983 were run. The default
spacing is 30 m on both. Three timed `--binary` runs per row (median), one
`--ascii` run for quality and the independent check; default threads (10).
Peak RSS from `/usr/bin/time -l`. "refine" is phase 1 (the stats file's
`refine` row); "other" is the stats file's unattributed time (see Surprises).

| out CRS | tol | triangles | vs hand-resampled | phase 2 inserted (rounds) | proc s | refine s | other s | peak RSS GiB | worst angle | max deg | CDT viol. | over tol: interior | strip |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| suggested tmerc | 50 | 76,361 | 1.020 | 680 (8) | 2.79 | 0.11 | 0.19 | 1.87 | 0.0828 | 13 | 0 | 0 of 13,791,979 | 0 of 27,799 |
| suggested tmerc | 20 | 276,076 | 1.042 | 5,712 (10) | 3.13 | 0.24 | 0.44 | 1.82 | 0.0828 | 13 | 0 | 0 of 13,791,979 | 0 of 27,799 |
| suggested tmerc | 10 | 709,134 | 1.091 | 29,508 (10) | 3.92 | 0.42 | 1.08 | 1.83 | 0.0417 | 13 | 0 | 0 of 13,791,979 | 0 of 27,799 |
| suggested tmerc | 5 | 1,860,101 | 1.180 | 143,257 (10) | 7.08 | 1.02 | 3.17 | 1.82 | 0.0487 | 13 | 0 | 0 of 13,791,979 | 0 of 27,799 |
| suggested tmerc | 1 | 13,877,141 | 1.595 | 2,588,049 (11) | 34.18 | 5.88 | 21.06 | 5.72 | 0.0006 | 14 | 0 | not run | not run |
| EPSG:31983 | 50 | 76,179 | 1.018 | 669 (7) | 2.48 | 0.10 | 0.16 | 1.81 | 0.1135 | 14 | 0 | 0 of 13,791,972 | 0 of 27,806 |
| EPSG:31983 | 20 | 276,205 | 1.043 | 5,686 (8) | 2.92 | 0.22 | 0.44 | 1.85 | 0.1135 | 14 | 0 | 0 of 13,791,972 | 0 of 27,806 |
| EPSG:31983 | 10 | 708,581 | 1.090 | 29,368 (11) | 3.88 | 0.43 | 1.15 | 1.87 | 0.0418 | 14 | 0 | 0 of 13,791,972 | 0 of 27,806 |
| EPSG:31983 | 5 | 1,861,967 | 1.181 | 142,643 (12) | 6.64 | 0.97 | 3.02 | 1.94 | 0.0303 | 14 | 0 | 0 of 13,791,972 | 0 of 27,806 |
| EPSG:31983 | 1 | 13,878,353 | 1.595 | 2,589,525 (12) | 32.56 | 5.76 | 20.06 | 6.36 | 0.0004 | 16 | 0 | not run | not run |

Failed attempts (retried): 7
- 2026-10-02T11:06:13Z tmerc_t50_run2.fail1.log: terminated abnormally (signal)
- 2026-10-02T11:06:20Z tmerc_t50_ascii.fail1.log: terminated abnormally (signal)
- 2026-10-02T11:11:13Z EPSG31983_t20_run3.fail1.log: terminated abnormally (signal)
- 2026-10-02T11:11:46Z EPSG31983_t5_run2.fail1.log: terminated abnormally (signal)
- 2026-10-02T11:12:02Z EPSG31983_t5_ascii.fail1.log: terminated abnormally (signal)
- 2026-10-02T11:12:04Z EPSG31983_t5_ascii.fail2.log: terminated abnormally (signal)
- 2026-10-02T11:14:03Z EPSG31983_t1_ascii.fail1.log: terminated abnormally (signal)

"vs hand-resampled" is the triangle count over the 2026-10-02 hand-resampled
ANADEM run (`../basin-piece-anadem/README.md`: 74,841 / 264,833 / 649,845 /
1,576,681 / 8,699,303 at 50 / 20 / 10 / 5 / 1 m). As expected, more, from
phase 2's insertions: at 1 m, 2.59 M inserted vertices make about 5.2 M
triangles, and 8.70 M + 5.18 M = 13.88 M.

**The independent check** is `../basin-piece-anadem/run_sweep.py`'s
(`_errors`: matplotlib's linear interpolator, not rasputin's locator),
re-pointed at the mesh's own CRS: every valid ANADEM node of the source window
(`../rasputin_data/sao_francisco_piece/bho2017_5k_76949_anadem_window_epsg4674.tif`),
projected by pyproj, inside the domain; the strip is the nodes within one grid
spacing (30 m) of the domain's edge. 13,791,979 interior and 27,799 strip
nodes in the suggested CRS (13,791,972 and 27,806 in EPSG:31983), the
interior count the 2026-10-02 README gives as 13.79 M. **Control**
(`runs/control/results.json`): the suggested-CRS 50 m mesh moved 15 m east
finds 1,156 interior nodes over tolerance; unmoved, 0.

## Surprises

1. **The random crash** (Verdict). The script retries a failed start (up to
   5 times), so every figure above is from a completed run; `runs/geo/results.json` `failures` lists each
   failed attempt and `runs/geo/logs/*.fail*.log` keeps its log.
2. **Phase 2's own timing rows always read 0.000.** In all 30 timed stats
   files, `final check: scan (parallel)` and `final check: split + flip
   (serial)` are 0.000 while phase 2 inserted up to 2.59 M vertices.
   `refine_points` never sets `scan_seconds` or `split_seconds`
   (`grep -n "scan_seconds\|split_seconds" include/terrain/refinement/refine_points.hpp`
   finds nothing), so `final_check.run` adds zeros, and phase 2's time lands
   in the stats file's `other` row (20 s of 32.6 s at 1 m). That `other` is
   phase 2 is thought to be so, not measured.
3. **The geographic path costs about 2-4x the hand-resampled run.** 1 m:
   32.6 s and 6.4 GiB, against 8.27 s and 3.11 GiB on the hand-resampled grid
   (phase 1 alone, no final check). The fixed cost at 50 m is 1.8-1.9 GiB
   against 0.62 GiB.
4. **Phase 2 lowers the worst angle.** EPSG:31983, 50 m: 0.1135° (the
   hand-resampled run's 0.113°); 10 m: 0.042°; 5 m: 0.030°; 1 m: 0.0004°. Max
   degree 14 to 16 at 1 m. Point insertion without a quality step; recorded,
   not judged, as master has no phase 2 to compare with.
5. **The suggested CRS has no datum.** The VTK's `source_transform` reads
   "Ballpark geographic offset from SIRGAS 2000 to unknown": `+ellps=GRS80`
   names an ellipsoid only, so pyproj uses a ballpark (null) datum shift.
   SIRGAS 2000 is on GRS80, so the numbers are unaffected here; the label is
   misleading, and on a non-GRS80 source it would not be harmless.

## Not done

Past the 60-minute budget (about 65 minutes), so scaled down:

- **The independent check at 1 m**: four matplotlib checks in parallel did
  not finish the two 1 m meshes (13.9 M triangles) within 30 minutes and were
  stopped; the rows read "not run". rasputin's own phase 2 reports max error
  1.0 m at its 18.8 M source nodes there, which is not an independent check.

- The 2 m tolerance and the GLO-30 comparison run of the increment's full
  acceptance list; phase 2's thread scaling at 1 m; phase 2 alone from the
  domain's start mesh. Ola's brief for this run named 50, 20, 10, 5 and 1 m on
  ANADEM.

## Reproduce

From the worktree root, the venv's python (with numpy, shapely, pyproj,
tifffile, imagecodecs, typer, pydantic, matplotlib) as `PY`,
`D=../rasputin_data` (absolute in the runs), `E=docs/benchmarks/2026-10-02/15c-2-acceptance`,
`B` a worktree of `e1f6042`, `TM` the suggested CRS above, and `SCRATCH`
and `CTRL` directories outside the repository.

```bash
git worktree add --detach $B e1f6042
PYTHONPATH=$B/build-bench/pkg $PY tools/bench.py run --label 15c-2-base-e1f6042 --tree $B \
  --threads 1,2,4,6,8,10,16,20 --repeats 5          # A; builds $B/build-bench first
$PY -c "import sys; sys.path.insert(0,'tools'); import bench; from pathlib import Path; \
print(bench.build(bench.make_runner(), Path('.').resolve()))"   # the branch's pkg, for PYTHONPATH
PYTHONPATH=$PWD/build-bench/pkg $PY tools/bench.py run --label 15c-2 --threads 1,2,4,6,8,10,16,20 \
  --repeats 5 --baseline docs/benchmarks/2026-10-02/15c-2-base-e1f6042   # B; then B2 and A2 the same way
O=$D/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson
W=$D/sao_francisco_piece/bho2017_5k_76949_anadem_window_epsg4674.tif
$PY $E/run_geo.py $O $W $D/cache $E/runs/geo $SCRATCH --crs "$TM;EPSG:31983" \
  --tolerances 50,20,10,5,1 --repeats 3 --no-check            # ~10 min
$PY $E/run_geo.py check $E/runs/geo/results.json $O $W 4      # the independent check
$PY $E/run_geo.py $O $W $D/cache $CTRL $SCRATCH --crs "$TM" --tolerances 50 --repeats 1 --shift 15
$PY $E/summarize.py                                           # the tables above
```

`bench.py`'s parent imports `tin_engine.stats`, so it needs the measured
tree's `build-bench/pkg` on `PYTHONPATH`. Meshes stayed in the session's
scratch directory and are not kept; the commands regenerate them. The ANADEM
window and the cache are described in `../basin-piece-anadem/README.md` and
the 23a-2 acceptance (`../23a-2-acceptance/`).

Data credits as in `../basin-piece-anadem/README.md` (ANADEM, CC BY 4.0;
BHO 2017, ANA).
