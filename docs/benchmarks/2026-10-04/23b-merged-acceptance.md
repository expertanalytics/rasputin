# Increment 23b acceptance at the merged head (@perf, 2026-10-04)

Branch `worktree-agent-a8b29c323dc2ad69d` at `4541e38` (source as of
`3403116`, the merge of master `d20126b` into 23b), against master `d20126b`
back to back. Same inputs and protocol as the `91c7cb5` rerun in
`23b-acceptance.md`.

## Verdict

**ACCEPTED.** At 1 thread refine is +0.8 % (tile) and +0.2 % (quarter)
against master; at 20 threads +1.5 % and +0.9 %; over 2 to 20 threads
-1.7 to +1.7 %. All are under `bench.py`'s 5 % threshold, and `bench.py`
judged all eight runs ACCEPTED. The run ranges of the two sides overlap at
1 thread on both domains. The meshes are identical to master's.

1. **Mesh hashes identical to master** (the unfrozen path, mask 0): tile
   `11741a81…`, quarter `a60fb597…`, the same in all eight runs, master and
   merged head alike. The quarter hash is not the `ccebf96a…` of the
   `91c7cb5` run because master now includes increment 15f-3, which changed
   the quarter mesh (`15f-3-acceptance.md`). Quality is identical, and every
   (domain, thread count) cell has a single (max_error, rounds, inserted,
   flips) tuple over all eight runs.
2. **The 1-thread margin is thinner than at `91c7cb5`**: +0.8 % / +0.2 %
   now, against +1.2 % / +1.8 % then. The merge's conflict resolution in
   `quality.hpp` and `refine.hpp` shows no measurable cost.

## Results: refine_s, pooled median (seconds), merged head against master

B is master `d20126b` (4 runs), F is the merged head `4541e38` (4 runs). The
full table, every thread count with run ranges, is
`23b-merged-acceptance/tables-merged.md`.

| domain | threads | master | merged head | change | master runs | merged-head runs |
|---|---:|---:|---:|---:|---|---|
| tile | default | 0.1973 | 0.1995 | +1.1 % | 0.1967-0.1979 | 0.1979-0.2015 |
| tile | 1 | 0.4833 | 0.4872 | +0.8 % | 0.4754-0.4988 | 0.4818-0.5036 |
| tile | 2 | 0.3203 | 0.3194 | -0.3 % | 0.3137-0.3287 | 0.3165-0.3300 |
| tile | 4 | 0.2386 | 0.2400 | +0.6 % | 0.2370-0.2442 | 0.2376-0.2452 |
| tile | 8 | 0.1956 | 0.1976 | +1.0 % | 0.1947-0.1996 | 0.1960-0.2012 |
| tile | 20 | 0.1991 | 0.2021 | +1.5 % | 0.1987-0.1997 | 0.2009-0.2034 |
| quarter | default | 0.1771 | 0.1791 | +1.1 % | 0.1760-0.1803 | 0.1774-0.1826 |
| quarter | 1 | 0.4270 | 0.4276 | +0.2 % | 0.4212-0.4437 | 0.4253-0.4497 |
| quarter | 2 | 0.2809 | 0.2791 | -0.6 % | 0.2750-0.2904 | 0.2782-0.2914 |
| quarter | 4 | 0.2115 | 0.2115 | +0.0 % | 0.2078-0.2172 | 0.2107-0.2175 |
| quarter | 8 | 0.1734 | 0.1744 | +0.6 % | 0.1712-0.1755 | 0.1732-0.1803 |
| quarter | 20 | 0.1782 | 0.1798 | +0.9 % | 0.1763-0.1788 | 0.1791-0.1812 |
| both | 2-20, range | | | -1.7..+1.7 % | | |

At 9 to 20 threads the merged head is 0.5 to 1.7 % slower in every cell, and
at 20 threads its run range sits just above master's on both domains. That is
under the threshold and about the same-binary spread seen in
`23b-acceptance.md` (about 2 %); its cause was not measured.

**Scaling ceiling** (pooled, 1 thread over 20), AC: master tile 2.43x and
quarter 2.40x; merged head 2.41x and 2.38x. The `91c7cb5` run's base
`b4bcdc3` gave 2.60x and 2.55x: master is faster at 1 thread (tile 0.4833 s
against 0.5158 s) at the same 20-thread time, so the ratio is lower. Which of
master's merges since `b4bcdc3` gives the faster single thread was not
measured.

## Quality (identical in all eight runs)

| domain | worst angle | max degree | within tolerance | Delaunay violations | mesh sha256 |
|---|---:|---:|---|---:|---|
| tile | 0.6296° | 74 | yes | 0 | `11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c` |
| quarter | 0.3955° | 18 | yes | 0 | `a60fb597b5b738268a3bd27730989392ee86fb030fe060b78fd3a0bcdbb60a34` |

## Method

- `tools/bench.py run`, blob `2b7dca1aa1184afea53b33c357ab3b836403d2f3` (the
  branch's, the same blob as the `91c7cb5` run), driving both trees through
  `--tree`, `--hardening on`, `--repeats 5`. Defaults otherwise: DEM
  `tests/fixtures/dem_archive/7908_3_10m_z33.tif`, tolerance 1 m, domains
  `tile` and `quarter`, threads 0 (CLI default) and 1 to 20, interleaved.
  Release `_core` from `bench.build()` into each tree's `build-bench/`:
  master `077f0bfb…`, merged head `318e43f2…`.
- Master was a detached worktree of `d20126b` under `.claude/worktrees/`,
  built in its own `build-bench/`, and removed after the run. The main
  checkout was not built in.
- Apple M1 Max (8 P + 2 E, 32 GiB), macOS 27.0, AppleClang 21.0.0.
  `caffeinate -i` held during the batch. No other agent or build ran; load
  average about 2 before and after. `pmset -g batt` before and after every
  run (`pairs7.log`, and each `run.json`): AC Power, 100 %, charged, every
  time.
- Two balanced batches in opposite order, B F F B then F B B F
  (`23b-merged-acceptance/pairs7.sh`, `pairs7.log`), 07:25 to 07:51 UTC. The
  tables are from `summarize_merged.py` (`tables-merged.md`).
- Merged-head runs are marked dirty because of the untracked evidence
  directories; `git status` showed no other change. Master runs are clean.

## Files and clean-up

- Run directories in `bench.py`'s format: `23bm-base-r1..r4/` (master) and
  `23bm-r1..r4/` (merged head).
- `23b-merged-acceptance/`: the driver, its log, the summarizer and its table.
- The quality meshes were in this session's scratchpad and have been deleted.
  To regenerate: `git worktree add --detach <dir> d20126b`, build each tree
  with `bench.build(bench.make_runner(), Path(T), "on")`, then
  `pairs7.sh W B S`; the meshes land in `S/meshes/<label>/`.

## Log
- 2026-10-04T07:24:35Z: file created; AC 100 % charged; load ~2.
- 2026-10-04T07:25:25Z: both built (master _core 077f0bfb..., merged head _core 318e43f2...); batch 1+2 started, B N N B then N B B N, caffeinate -i.
- 2026-10-04T07:45:55Z: batch 1 (B F F B) done, all AC 100 %.
- 2026-10-04T07:52:26Z: batch 2 (F B B F) done, all AC 100 %. +0.8 / +0.2 % at 1 thread, -1.7..+1.7 % over 2-20; meshes identical. ACCEPTED.
