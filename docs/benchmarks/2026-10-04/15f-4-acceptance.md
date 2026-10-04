# Increment 15f-4 acceptance: a cheap mesh rebuild (@perf, 2026-10-04)

Branch `worktree-15f-4` at `5eeb87a` (green, not yet reviewed; stacked on
15f-3). It is measured against 15f-3 at `107cdb4` (its `_core` is
byte-identical to `83c7fd2`'s, `0ee8cd40…`) and against master's merge base
`c193cb1`, back to back. AC power for every run: 100 %, charged, `pmset -g
batt` before and after each. Swap stayed at 10.53 GB and did not grow. The
rule is `docs/increments/15f-edge-strip.md`, A5, and A2's trigger: the
empty-strip `refine_strip` call on Bygdin 1 m must take under 0.43 s.

## Verdict

**ACCEPTED.**

- **Meshes byte-identical to 15f-3 everywhere measured:** the bench's
  tile (`11741a81…`) and quarter (`a60fb597…`), Bygdin 1 m (`661bf16b…`)
  and Velhas 5 m (`c9087951…`).
- **Refine within noise or faster:** -3.1 to +0.5 % against 15f-3 over
  every cell, and -3.0 to +0.3 % against the base.
- **A2's trigger is not met.** The empty-strip `refine_strip` call on
  Bygdin 1 m takes **0.355 s** (median of 10, range 0.351-0.361 s),
  against the 0.43 s threshold. 15f-3 took 1.304 s in the same session, so
  the rebuild is 3.7 times cheaper. Fusion stays closed.

**The morning question, "15f-3 costs a third end to end", answered:**

| process time, median | base c193cb1 | 15f-3 | 15f-4 | 15f-4 vs 15f-3 | 15f-4 vs base |
|---|---:|---:|---:|---:|---:|
| 1 m bench, tile, 1 thread | 0.939 s | 1.289 s | 1.045 s | -19.0 % | +11.3 % |
| 1 m bench, tile, 20 threads | 0.620 s | 0.963 s | 0.725 s | -24.7 % | +17.1 % |
| 1 m bench, quarter, 1 thread | 0.826 s | 1.137 s | 0.931 s | -18.1 % | +12.7 % |
| 1 m bench, quarter, 20 threads | 0.549 s | 0.851 s | 0.648 s | -23.8 % | +18.0 % |
| Bygdin 1 m (6 runs a side) | 4.39 s | 5.89 s | 4.78 s | -18.8 % | +8.9 % |
| Velhas 5 m (4 runs a side) | 6.83 s | 7.31 s | 5.88 s | -19.5 % | -13.8 % |

- **Projected path:** 15f-3's extra cost falls from +34 to +55 % to +9 to
  +18 % over the base.
- **What remains on Bygdin 1 m:** the `other` row is 0.395 s (15f-3:
  1.425 s; base: 0.005 s), which is close to the empty-strip call. The
  strip's timed phases take about 0.04 s.
- **Reprojected path:** Velhas 5 m is now **faster than the base**, by
  13.8 %. That is because the final check rebuilds the same mesh (A5): its
  `other` row falls from 2.03 s (base) to 1.04 s.
- **Peak memory, Bygdin 1 m:** 1.04 GiB against 15f-3's 1.11 and the
  base's 0.82 (max of 6 runs).

## Method

- **Builds and tools.** Release `_core` from `bench.build()` in each tree:
  base `0b3362e5…`, 15f-3 `0ee8cd40…`, 15f-4 `3d2bf22c…`. 15f-4's
  `tools/bench.py` (blob as in each `run.json`) drives every tree through
  `--tree`. The 15f-4 worktree has no `.venv`, so 15f-3's venv runs
  everything.
- **Order.** `15f-4-acceptance/batches.sh` (log `batches.log`) runs batch 1
  as base, 15f-3, 15f-4 and batch 2 as 15f-4, 15f-3, base,
  under `caffeinate -i`. No other agent or build ran.
- **Per slot:**
  - `bench.py run`: tile and quarter, threads 0 and 1 to 20, 5 repeats;
    run directories `15f-4-{base,15f3,15f4}-b{1,2}/`;
  - Bygdin 1 m with CORINE: 3 `--binary` runs under `/usr/bin/time -l`,
    plus one `--ascii` run in batch 1 for the hash (`bygdin/`);
  - Velhas 5 m in EPSG:31983: `velhas.py` (15f-3's, with `PY` pointed at
    15f-3's venv), 2 timed runs and 1 `--ascii` run (`velhas/`).
- **The empty-strip calls.** `prof_strip.py` (15f-3's) with K = 5 on
  Bygdin 1 m, in the order 15f-4, 15f-3, 15f-3, 15f-4, giving 10 full and
  10 empty calls a side (`empty_strip/`). Full strip: 15f-4 0.402 s, 15f-3
  1.428 s. Empty: 0.355 s and 1.304 s. The timed phases were 0.03 s full
  and 0.004 s empty on both.
- **Tables.** Everything comes from `summarize.py` (`tables.md`): pooled
  medians, every thread count, `app_s` and `proc_s`, and the `--stats`
  phases.
- The sanitize-first rule does not apply: no C++ was patched by `@perf`.

## Results in brief (from `tables.md`)

- **Refine (`refine_s`), 15f-4 against 15f-3:** tile -1.9 % at 1 thread
  and -1.9 % at 20; quarter -0.6 % and -1.6 %; -3.1 to +0.5 % over every
  cell.
- **Quality**, identical on all three trees where the mesh is: tile worst
  angle 0.6296°, max degree 74; quarter 0.3955°, 18; 0 Delaunay violations
  in all six runs.
- **Bygdin 1 m `--stats` phases (medians):**

  | phase | base | 15f-3 | 15f-4 |
  |---|---:|---:|---:|
  | refine | 0.600 | 0.608 | 0.597 |
  | strip (generate + scan + split) | - | 0.042 | 0.040 |
  | other | 0.005 | 1.425 | 0.395 |

- **Velhas 5 m `--stats` phases (medians):**

  | phase | base | 15f-3 | 15f-4 |
  |---|---:|---:|---:|
  | final check scan | 0.742 | 0.788 | 0.791 |
  | final check split | 0.168 | 0.175 | 0.180 |
  | other | 2.03 | 2.42 | 1.04 |

## Files and clean-up

- `15f-4-acceptance/`: `batches.sh` and `batches.log`, `summarize.py` and
  `tables.md`, `bygdin/`, `velhas/`, `empty_strip/`, and the copies of
  `velhas.py`, `prof_strip.py` and `strip_check.py` that were used.
- The meshes and the two detached worktrees (`c193cb1`, `107cdb4`) were in
  this session's scratchpad and have been deleted. To regenerate them,
  recreate the worktrees, build each tree with `bench.build`, and run
  `batches.sh W B P S`.

## Log
- 2026-10-04T01:18:32Z: start. AC 100 %, swap 10.53 GB used (flat). The worktree has no .venv; 15f-3's venv (same lock) runs every tree.
- 2026-10-04T01:43:24Z: done. ACCEPTED; empty-strip 0.355 s < 0.43 s.
