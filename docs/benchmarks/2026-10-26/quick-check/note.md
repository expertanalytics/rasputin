# Quick check: is the 1 m benchmark still about 40 % slower than on 2026-10-25?

**No.** On a quiet machine the run matches 2026-10-25's stored run within
noise. The slowdown in this morning's acceptance runs came from the machine,
not from the code.

## What was compared

All three runs are `tools/bench.py run`, Release, bounds checks on (libc++
fast), on AC power, on the same M1 Max. The DEM is `7908_3_10m_z33.tif`,
tolerance 1.0, domains `tile` and `quarter`, 5 repeats, and the figures are
medians of refine time in seconds.

| Domain, threads | 2026-10-25 20:50 (`104b-head-ac-r2`) | 2026-10-26 10:34 (`106c-base-ac-r1`) | now, 16:01 (`quick-check`) | now vs 10-25 |
|---|---|---|---|---|
| tile, 0 (default) | 0.1942 | 0.2758 | 0.1909 | -1.7 % |
| tile, 1 | 0.4770 | 0.6509 | 0.4821 | +1.1 % |
| tile, 4 | 0.2326 | 0.3192 | 0.2349 | +1.0 % |
| tile, 8 | 0.1911 | 0.2735 | 0.1919 | +0.4 % |
| tile, 20 | 0.1947 | 0.2716 | 0.1931 | -0.8 % |
| quarter, 0 (default) | 0.1656 | 0.2044 | 0.1750 | +5.7 % |
| quarter, 1 | 0.4223 | 0.5655 | 0.4229 | +0.1 % |
| quarter, 4 | 0.2057 | 0.2714 | 0.2084 | +1.3 % |
| quarter, 8 | 0.1685 | 0.2071 | 0.1705 | +1.2 % |
| quarter, 20 | 0.1742 | 0.2105 | 0.1762 | +1.2 % |

Over all 21 thread counts the median ratio of now to 2026-10-25 is 1.006
(tile) and 1.012 (quarter). `bench.py` flags only quarter at t=0 (+5.7 %,
just over its 5 % threshold). That is a single cell with no matching move in
the tile, so it is read as noise.

## Why no merged change can be the cause

All three runs loaded the same `_core` build: sha256 `a353ac1e…` in each
`run.json`. No C++ or CMake file changed between `0d751480` (the 10-25 run's
tree) and `fc55c699` (this run's tree). The only `pyproject.toml` change adds
a file to mypy's list. The meshes are identical too: tile `11741a81`, quarter
`a60fb597`, the same as in the 106c runs. So quality is unchanged.

## Machine state

When this run started, the load averages were 8.49, 10.44 and 12.53 (over 1,
5 and 15 minutes); when it finished they were 4.80, 7.83 and 10.99. The load
was falling from earlier work. The 10:34 runs' load was not recorded, so the
cause of the morning's slowdown is thought to be desktop load and was not
profiled.

The raw data are in `run.json` and `raw.tsv` beside this note. The meshes were
not committed; to regenerate them, run
`tools/bench.py run --label quick-check --domain tile --domain docs/benchmarks/2026-09-26/quarter.geojson --baseline docs/benchmarks/2026-10-25/104b-head-ac-r2`
at `fc55c699`.
