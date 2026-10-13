# Quick-check baseline re-saved for increment 45 (@perf, 2026-10-13)

Asked by the main session: increment 45's last step
(`docs/increments/45-linux-defaults.md`, section 4.4, last paragraph). The
increment adds the `velhas` case, stores each case's input list in the record
and names the baseline file after the machine and power state, so the old
`docs/benchmarks/quick/baseline-ac.json` (commit `969c49e2`, 2026-10-12) no
longer applies. It is replaced by `docs/benchmarks/quick/baseline-apple-m1-max-ac.json`
and deleted in the same commit. Not an acceptance run: no refine or mesh code
is touched.

## Method

- Tree: the `linux-defaults` worktree at `87df7c96` (increment 45, code review
  approved); script version `tools/bench_quick.py` and `tools/bench.py` at
  that commit. The tool builds `_core` in Release into the tree's
  `build-bench` (hardening on, the default) and times from there; the venv's
  own extension was also rebuilt and reinstalled first.
- Command: `RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data
  .venv/bin/python tools/bench_quick.py run --save-baseline`, log
  `save-run.log` (with `pmset -g batt`, `uptime` and the time before and
  after).
- Machine: Apple M1 Max (8 P + 2 E cores, 32 GB), macOS 27.0, Python 3.14.7,
  numpy 2.5.4, shapely 2.2.0 / GEOS 3.14.1, pyproj 3.8.0 / PROJ 9.8.1. Power:
  AC, battery 100 %, before and after; no thermal warning recorded. Load
  average 4.4 / 5.9 / 8.1 before (falling after the extension build), 3.1
  after; nothing of ours beside the run. 21:40 to 21:43 local.
- Threads: default, plus `tile t=1`. Runs, DEMs, domains and tolerances: as
  in `docs/benchmarks/quick/cases.toml` at `87df7c96`.

## Result

Verdict printed: `NO BASELINE: no stored baseline` (exit 2), as expected for a
first save under the new file name. Against the old record by hand
(`compare.txt`; medians of `total`, seconds):

| case | inputs | old (`969c49e2`) | new (`87df7c96`) | change | new spread | max error / tolerance | mesh identical |
|---|---|---|---|---|---|---|---|
| tile | 1 | 0.430 | 0.425 | -1.2 % | 9.2 % | 0.99998 / 1 | yes |
| tile t=1 | 1 | 0.711 | 0.725 | +2.0 % | 3.3 % | 0.99998 / 1 | yes |
| quarter | 2 | 0.366 | 0.364 | -0.7 % | 5.4 % | 0.99998 / 1 | yes |
| numedalslagen | 3 | 8.894 | 8.846 | -0.5 % | 0.4 % | 10.0 / 10 | yes |
| lagan | 3 | 14.736 | 14.498 | -1.6 % | 0.2 % | 9.9999 / 10 | yes |
| geilo-al-ramp | 3 | 0.344 | 0.334 | -3.1 % | 2.4 % | 19.98 / 20 | yes |
| romsdal-slope | 1 | 3.624 | 3.588 | -1.0 % | 0.7 % | 9.9999 / 10 | yes |
| velhas | 1 | (new) | 0.125 | new case | 2.0 % | 0.49994 / 0.5 | - |

Every old case is within 3.1 % of the old baseline with a byte-identical mesh,
so the new record is not taken on a loaded machine. All seven cases (eight
runs with `tile t=1`) are present, each with a non-empty `inputs` list; power
`ac`. Hotspot lines (`save-run.log`) are the same seven as the 2026-10-12
refresh, within 1 point each; `velhas` has none.
