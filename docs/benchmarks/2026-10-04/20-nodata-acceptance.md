# Increment 20 fix acceptance: the start-quality pass skips NoData nodes (@perf, 2026-10-04)

Branch `worktree-20-nodata`, HEAD `cce79b0` (code `0d60c69`), against
origin/master `2060f14`, in two reversed batches. AC power for every run:
100 %, charged, `pmset -g batt` before and after each. Swap did not grow
(10,455 → 10,447 MB). The rule is `docs/increments/20-start-quality.md`,
"`@perf`'s acceptance".

## Verdict

**ACCEPTED.**

- **NoData inside the benchmark domains.** `count_nodata.py`
  (`nodata_count.txt`) reads the tile with tifffile and does not import
  tin_engine. The tile has 42,476 NoData nodes (rows 0-2876). The quarter
  domain has 0.
- **Mesh hashes unchanged on both domains:** tile `11741a81…`, quarter
  `ccebf96a…`, the same in all eight runs. Quality is identical (worst
  angle 0.6296° / 0.3955°, max degree 74 / 18, 0 Delaunay violations), and
  so are refine's counters.
- **How many NoData nodes the pass met on the tile: none.**
  - The `--stats` start-quality rows on the tile are identical: 503
    inserted and 252 skipped on both master and the fix
    (`tile_base.stats.md`, `tile_fix.stats.md`; 199 NoData vertices removed
    by trim on both).
  - The fix adds `skipped_void` to the skipped total. With the total, the
    inserts and the whole mesh unchanged, `skipped_void` is 0.
  - So the tile's NoData block lies where the start-quality pass never
    snaps. This is inferred from the counts, not from a debug print.
- **Refine within noise.** -0.5 to +1.3 % over every cell of the sweep
  (threads default and 1 to 20, both domains). At 1 / 20 threads it is
  +0.9 / -0.0 % on the tile and +0.4 / +0.1 % on the quarter. Process time
  is within -0.4 to +0.7 % at the listed cells.

## Method

- **Builds.** Release `_core` from `bench.build()`: master in a detached
  worktree `d63f62c3…`, the fix `63b4cebb…`. Both were built with this
  worktree's `.venv` python. Master's hash differs from earlier nights'
  `0b3362e5…`, which were built with another venv's python.
- **Runs.** `20-nodata-acceptance/pairs.sh` (log `pairs.log`): batch 1
  master, fix, fix, master; batch 2 fix, master, master, fix. 5 repeats,
  threads 0 and 1 to 20, both domains, `caffeinate -i`. It stops off AC or
  on swap growth over 100 MB.
- **Memory at the start.** The compressor held 13.9 GB of RAM (Safari's
  WebKit processes), 66 MB of pages were free, and swap stood at 10.46 of
  11.26 GB. The bench runs peak under 1 GB.
- **Tables.** `summarize.py` → `tables.md` (pooled medians, run ranges,
  `app_s`, `proc_s`, quality and counters).
- Run directories: `20nd-base-r{1..4}/` and `20nd-r{1..4}/`.
- The meshes and the master worktree were in the scratchpad and are
  deleted.

## Log
- 2026-10-04T04:16:46Z: AC 100 %. Memory: compressor 13.9 GB occupied, 66 MB free pages, swap 10.46 of 11.26 GB (Safari's WebKit processes hold the compressed memory). bench runs peak under 1 GB.
- NoData count (`count_nodata.py`): the tile has 42,476 NoData nodes (rows 0-2876); the quarter domain has 0.
