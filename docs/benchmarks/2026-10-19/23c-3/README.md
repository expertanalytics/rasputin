# 23c-3 acceptance: each piece's memory its own (@perf, 2026-10-19)

`worktree-23c-3` at `09c8aad4` (code review approved on `fefdff59`; the later
commits are prose). Base: 23c-2, `worktree-23c` at `321bac56`, accepted on
2026-10-18 (`docs/benchmarks/2026-10-18/23c-2/`). Master `cef9f3f0` measured
too, for D8. Base and master ran from detached scratch worktrees. The rule is
`docs/increments/23-basin-scale.md`, "23c-3: each piece's memory its own",
the list under "`@perf` applies".

AC power for every run, 100 %, charged (`pmset -g batt` before and after
each: in each `run.json`, `pairs*.log`, `logs/*.pmset`). Swap 6.03 to
6.06 GB used of 7.0 GB from the first run to the last, no growth. Apple M1
Max, 8 P + 2 E cores, 32 GiB, macOS 27.0.1. No other run of mine beside a
timing run; Ola's desktop apps were open (a browser tab at about 50 % of
one core; load average 2.6 to 4.6; `cpu2.log` samples the top users every
15 s during the second bench batch).

## Verdict

**ACCEPTED.** Meshes and piece files unchanged, memory now follows the piece.
Three findings, none of them a gate in the design:

1. `bench.py` reports REGRESSION on the tile (refine +28 %) in all four head
   runs. With the same `_core` binary and the same mesh, this measured as a
   process-layout effect, not the branch (part 1).
2. The pieces' peak exceeds `dx · dy · b(T)` + 72 MiB by 31 to 255 MB, so
   `b(T)` does not cover the reprojected path's extras. The design's "Left
   out" (`b(T)` re-measured with phase 2) is needed.
3. Reading each piece separately costs time: on Velhas at 10 m with 64
   pieces, decode takes 5.00 s against 1.31 s and the whole run 9.0 s
   against 4.9 s.

## 1. Bench (`tools/bench.py`, the uncut path)

- **K1, mesh hashes unchanged** in all ten runs: tile `11741a81…`, quarter
  `a60fb597…`. Quality is the same: worst angle 0.6296° / 0.3955°, max
  degree 74 / 18, within tolerance, 0 Delaunay violations (692,056 /
  641,792 edges checked).
- **The `_core` binary is the same file** on all three trees: sha256
  `52390c7f…` (Release, `--hardening on`). `git diff --stat cef9f3f0
  09c8aad4 -- include src bindings CMakeLists.txt` is empty.
- **Quarter (D8): process time within noise of master's.** Pooled over 4
  head, 4 base and 2 master runs, head against master: -0.8 to +0.6 %
  across default, 1, 8 and 20 threads. Head against 23c-2: -0.6 to -0.3 %.
  Over the 32-layout sweep below: head 0.6261 s, 23c-2 0.6298 s, master
  0.6266 s (default threads). 23c-2's +0.6 to +1.8 % against master is gone.
  - After an uncut quarter run, `sys.modules` holds none of `basin_run`,
    `pieces` and `io.mesh_index` on head and on master; 23c-2 holds all
    three (`import_probe.txt`).
  - Median `import tin_engine.cli` over 8 fresh processes: master 185.0 ms,
    23c-2 190.7 ms, head 188.7 ms.
- **Tile: bench.py's REGRESSION is a layout effect.**
  - Every head run reports the tile's refine 25 to 31 % slower at 2 or more
    threads and 8 % at 1 (`pairs.log`, `pairs2.log`; second batch in
    reversed order, W B B W). The quarter in the same runs is unchanged.
  - The binary and the mesh are the same, and 23c-3 changes no tile-path
    code that refine sees: `tile_probe.py` shows the same array (float32
    5051², C order, a read-only view, page-aligned) on both sides.
  - Outside `bench.py run`, the slow mode follows the process's
    environment size, not the tree. 23c-2 with an extra variable `PAD` of
    0 or 96 bytes: 0.227-0.240 s, 3 of 3 each; with 16 bytes:
    0.188-0.192 s (`env_pad_repeat.txt`). In bench.py, the trees' package
    paths differ in length (the head's is in the repository, the others in
    the scratchpad), so each tree gets its own layout. Why a layout costs
    25 % was not measured; it is thought to be alignment or aliasing in
    memory.
  - **Over 32 layouts** (`layout_sweep.py`: `PAD` 0 to 248 bytes, trees
    interleaved, one sample each), the tile's refine median is head 0.1875
    s, 23c-2 0.1882 s, master 0.1915 s (default threads), and head 0.4773,
    23c-2 0.4854, master 0.4804 s at 1 thread. No difference.
- **Scaling ceiling, tile refine, 1 thread to the default (`--threads` unset, the core's own count)**, from the
  32-layout sweep's medians: head 2.55× (0.4773 / 0.1875), 23c-2 2.58×,
  master 2.51×. bench.py's head figure (2.17×) is in the slow layout.

| domain | threads | measure | M | B | N | N/B | N/M |
|---|---:|---|---:|---:|---:|---:|---:|
| tile | default | refine_s | 0.1914 | 0.1881 | 0.2410 | +28.1 % | +25.9 % |
| tile | 1 | refine_s | 0.4788 | 0.4826 | 0.5218 | +8.1 % | +9.0 % |
| tile | default | proc_s | 0.6937 | 0.6953 | 0.7486 | +7.7 % | +7.9 % |
| quarter | default | refine_s | 0.1669 | 0.1671 | 0.1670 | -0.1 % | +0.1 % |
| quarter | 1 | refine_s | 0.4302 | 0.4230 | 0.4224 | -0.1 % | -1.8 % |
| quarter | default | proc_s | 0.6262 | 0.6329 | 0.6289 | -0.6 % | +0.4 % |
| quarter | 1 | proc_s | 0.8897 | 0.8853 | 0.8826 | -0.3 % | -0.8 % |
| quarter | 20 | proc_s | 0.6348 | 0.6339 | 0.6316 | -0.4 % | -0.5 % |

M = master `cef9f3f0`, B = 23c-2 `321bac56`, N = 23c-3 `09c8aad4`; pooled
medians of bench.py's samples. Every cell and both sweeps are in `tables.md`.

## 2. The Velhas piece: ANADEM from the cache, EPSG:31983

- **K6 on real data: piece files byte-identical to 23c-2's.** `cmp` of every
  file in `velhas.vtk.pieces/`, `index.json` included, between 23c-2 and
  23c-3 at 10 m with `--pieces` 16, 64 and 256 (13, 51 and 154 pieces) and
  at 1 m with 64 (51 pieces): 0 files differ out of 14, 52, 155 and 52. At
  10 m with 64 pieces the start pass inserts 14,951 (median 321, range
  0-566), as in 23c-2's acceptance.
- **Memory follows the piece.** Peak footprint at 10 m:

| `--pieces` | 23c-2 | 23c-3 | largest cell `dx · dy · b(T)` |
|---:|---:|---:|---:|
| 16 | 1,733 MB | 407 MB | 76 MB |
| 64 | 1,689 MB | 155 MB | 17 MB |
| 256 | 1,733 MB | 111 MB | 5 MB |
| uncut | | 1,686 MB | |

  It falls as pieces shrink, where 23c-2's does not. At 1 m with 64 pieces:
  23c-2 1,711 MB, 23c-3 294 MB.
- **The bound is not met.** Peak less `dx · dy · b(T)` less 72 MiB:

| run | `b(T)` | cell | peak | excess |
|---|---:|---|---:|---:|
| 50 m, 64 | 19 | 701 × 668 | 155 MB | 71 MB |
| 20 m, 64 | 25 | 701 × 668 | 155 MB | 68 MB |
| 10 m, 16 | 37 | 1402 × 1470 | 407 MB | 255 MB |
| 10 m, 16, `--threads 1` | 37 | 1402 × 1470 | 249 MB | 97 MB |
| 10 m, 64 | 37 | 701 × 668 | 155 MB | 62 MB |
| 10 m, 256 | 37 | 351 × 350 | 111 MB | 31 MB |
| 5 m, 64 | 65 | 701 × 668 | 156 MB | 50 MB |
| 2 m, 64 | 155 | 701 × 668 | 202 MB | 53 MB |
| 1 m, 64 | 267 | 701 × 668 | 294 MB | 94 MB |

  The excess has to hold the vector step, which is at most 31 MB (the
  256-piece run's whole excess). Beyond it, `b(T)` misses up to 224 MB, at
  10 m with 16 pieces: about 110 B per node of the cell, against `b(T)` =
  37. Not profiled. With 1 thread the 16-piece peak is 249 MB; with 2 or 10
  threads it is 404 or 407 MB (`logs/head_t10_p16_th*`), so the extra does
  not grow with the thread count past 2. Its cause is not measured.
- **Decode costs more.** Phases summed over pieces, 23c-2 against 23c-3:
  at 10 m with 64 pieces, decode 1.31 s against 5.00 s and the check points
  (project and store) 1.49 s against 2.64 s; refine 0.43 s against 0.42 s.
  At 16 pieces, decode 1.27 s against 3.04 s; at 256, 1.26 s against
  4.23 s; at 1 m with 64, 1.27 s against 4.83 s.

| tolerance | cut 64: wall | peak | start inserts (sum; median; range) | uncut: wall | peak | points |
|---:|---:|---:|---|---:|---:|---:|
| 50 m | 8.32 s | 155 MB | 8,871; 140; 0-532 | 2.89 s | 1,731 MB | 40,075 |
| 20 m | 8.52 s | 155 MB | 11,699; 228; 0-546 | 3.28 s | 1,739 MB | 140,547 |
| 10 m | 9.04 s | 155 MB | 14,951; 321; 0-566 | 3.97 s | 1,686 MB | 357,592 |
| 5 m | 9.81 s | 156 MB | 20,554; 459; 0-759 | 5.78 s | 1,707 MB | 933,232 |
| 2 m | 13.38 s | 202 MB | 34,320; 727; 0-1442 | 12.56 s | 2,186 MB | 3,204,612 |
| 1 m | 19.17 s | 294 MB | 48,630; 967; 0-2017 | 24.47 s | 3,937 MB | 6,945,120 |

  An uncut run's start pass inserts 8,036 at every tolerance. Maximum RSS for
  every run, and per-piece peaks, are in `tables.md` and `logs/*.jsonl`.
- **Basin smoke at 50 m: not run.** `rasputin fetch anadem-v1 --domain
  bho2017_level2_76_raw.geojson --out-crs EPSG:31983 --dry-run` reports 100
  of 8,700 blocks to fetch (`basin_dry_run.txt`), so per the design it waits
  for B11's basin run.

## Method

- **Bench:** `pairs.sh` (M B W W B M) and `pairs2.sh` (W B B W); 5 repeats,
  threads default and 1 to 20, domains `tile` and `quarter`, `--hardening
  on`, `caffeinate -i`. Each tree built fresh by bench.py into its own
  `build-bench`; bench.py blob `16ad2ff5…` in all three trees. Runs in
  `2026-10-19/`. `master-r1` reports NO BASELINE (exit 2), being the first.
  The head runs are marked dirty: the only untracked files were this
  directory's. The layout sweeps: `layout_sweep.py`
  (`layout_sweep_t{0,1}.tsv`). `summarize.py` writes `tables.md`.
- **Velhas:** `velhas_all.sh` runs `velhas.sh` for each case: the tree's
  Release package (`build-bench/pkg`) through `velhas_driver.py`, under
  `/usr/bin/time -l`, `MallocLargeCache=0`, default threads, `--binary`,
  `mesh --dem anadem-v1 --cache ~/.cache/rasputin --domain
  bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983
  --tolerance T [--pieces P]`. The driver records each `cli.refine` call
  and the peak footprint before and inside it from `proc_pid_rusage`; the
  process's lifetime peak there equals `/usr/bin/time`'s "peak memory
  footprint" byte for byte (checked on the 50 m run). It also keeps the
  run's phase clock, since `--stats` is refused on a cut run.
- **Not kept:** the piece files and the meshes (in the session scratchpad,
  deleted) and the scratch worktrees. Rerun `velhas_all.sh W B S` and
  `pairs.sh` to regenerate them.
