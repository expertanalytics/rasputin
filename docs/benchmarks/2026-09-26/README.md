# Benchmarks, 2026-09-26/27: 1 m progression and thread scaling

This is the evidence for `docs/retrospectives/2026-09-27-increments-14-to-20b.md`.
It was copied out of the session scratchpad so that it survives the session.

The scripts are a one-off record. They hard-code scratchpad paths and git
worktrees that no longer exist. A reusable, checked-in tool (`tools/bench.py`,
owned by @perf) will replace them; it is not yet written. To rerun by hand, fix the paths at the top of
`bench1m/bench.sh` and `scaling/scale.py`'s caller.

- DEM: `tests/fixtures/dem_archive/7908_3_10m_z33.tif`. Domain: `quarter.geojson`
  (Ola's quarter circle). Tolerance: 1 m.
- Machine: Apple Silicon, 8 performance + 2 efficiency cores.

## bench1m/

One run per increment, 14 to HEAD (08a7711), in Release worktrees. It ran on
**battery**. `REPORT.md` has the tables, the regressions and the method.
`quality/` has the per-run quality JSON and the `--stats` reports. The meshes
(about 380 MB) are not kept in the repo.

## scaling/

This measures the `_core.refine` call alone, with threads forced by
`scale.py N`, at 1 to 20 threads, taking the median of 5 runs. Code: master
8f47e7e.

- `raw_battery.tsv`: on battery.
- `raw_ac.tsv`: on AC. **The ceiling is confirmed on AC**, so it is not
  throttling.

| threads | quarter AC s | speedup | tile AC s | speedup |
|---|---|---|---|---|
| 1 | 0.499 | 1.00x | 0.539 | 1.00x |
| 2 | 0.360 | 1.39x | 0.379 | 1.42x |
| 4 | 0.280 | 1.78x | 0.303 | 1.78x |
| 8 | 0.247 | 2.02x | 0.255 | 2.11x |
| 10 | 0.243 | 2.05x | 0.262 | 2.06x |
| 16 | 0.230 | 2.17x | 0.254 | 2.12x |
| 20 | 0.231 | 2.16x | 0.256 | 2.10x |

Battery and AC agree to within about 7 % (the largest gap is the quarter circle at 7 threads:
0.258 s on AC against 0.241 s on battery), and they show the same ceiling. Roughly half of the single-thread
refine time does not parallelise. That part is thought to be the serial insert
and flip phase, but it has not been profiled yet.
