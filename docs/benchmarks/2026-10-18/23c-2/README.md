# 23c-2 acceptance after the master merge (@perf, 2026-10-18)

The `23c` worktree at `4e140864` (code review round 5 approved; master
`cb713c87` merged in, so `include/`, `src/` and `bindings/` equal master's).
Base: master `cb713c87` in a detached scratch worktree. The rule is
`docs/increments/23-basin-scale.md`, "23c-2 at green", point 8: the bench on
the uncut path (mesh hashes unchanged, K1, "small runs untouched") and one
smoke run of the Velhas piece cut into pieces. This replaces the 2026-10-04
run (`9b7a5b4`, `docs/benchmarks/2026-10-04/23c-2-smoke.md`), which predates
the master merge.

AC power throughout, 100 %, `pmset -g batt` before and after every run
(in each `run.json` and `velhas.pmset`); swap 0.25 MB before and after, no
growth. Apple M1 Max, 8 P + 2 E cores, 32 GiB, macOS 27.0.1. Nothing else
ran beside the timing runs.

## Verdict

**ACCEPTED on both parts.**

### 1. Bench (`tools/bench.py`, the uncut path)

- **Mesh hashes unchanged** on both domains, all four runs: tile
  `11741a81…`, quarter `a60fb597…`. Quality identical: worst angle 0.6296° /
  0.3955°, max degree 74 / 18, within tolerance, 0 Delaunay violations
  (692,056 / 641,792 edges checked). Refine counters at default threads
  identical (rounds, inserted, flips, max error).
- **The `_core` binary is the same file** on both sides: sha256
  `a902c2f4…` (Release, libc++ fast hardening, `--hardening on`). Only the
  Python layer differs.
- **Refine time:** -3.6 to +3.0 % over every cell of the sweep (threads
  default and 1 to 20, both domains). Threshold 5 %.
- **Process time:** tile -0.6 to +1.4 %; quarter +0.6 to +1.8 %, every
  cell positive. The quarter's process time is consistently about 1 % higher
  on the head; the binary is identical, so it lies in Python. Not profiled;
  it is under the 5 % threshold.
- **Scaling ceiling, refine 1 thread to 20 (tile):** base 2.48×, head 2.49×.

| domain | threads | refine B s | refine N s | change | process B s | process N s | change |
|---|---:|---:|---:|---:|---:|---:|---:|
| tile | default | 0.1879 | 0.1907 | +1.5 % | 0.6786 | 0.6844 | +0.9 % |
| tile | 1 | 0.4787 | 0.4738 | -1.0 % | 0.9702 | 0.9648 | -0.6 % |
| tile | 20 | 0.1933 | 0.1905 | -1.5 % | 0.6774 | 0.6849 | +1.1 % |
| quarter | default | 0.1691 | 0.1717 | +1.5 % | 0.6136 | 0.6196 | +1.0 % |
| quarter | 1 | 0.4188 | 0.4183 | -0.1 % | 0.8628 | 0.8702 | +0.9 % |
| quarter | 20 | 0.1700 | 0.1701 | +0.0 % | 0.6130 | 0.6204 | +1.2 % |

B = master `cb713c87`, N = 23c-2 `4e140864`; pooled medians over 10 samples
a cell per side. Every cell, and app time, is in `tables.md`.

### 2. The Velhas smoke: ANADEM from the cache, EPSG:31983, 10 m, `--pieces 64`

- **Completion:** exit 0; **51 piece files** and `index.json` (51 entries)
  from the 6 × 11 partition; every piece's refine status `Ok`.
- **Conformity holds over the piece files:**
  - 3,305 constraint segments appear in exactly two pieces, with identical
    end points and identical z (0 mismatches); none in more than two.
  - The 7,490 segments only one piece has all lie on the domain's outline:
    midpoints within 0.66 mm of the outline projected into EPSG:31983
    (`once_on_outline.txt`).
  - Points the pieces share have identical z (largest difference 0.0).
- **The independent check finds 0 over tolerance** (15c-2's
  `run_geo.independent_check` on the merged pieces): 13,791,972 interior
  source nodes (largest error 9.99997 m) and 27,806 strip nodes (9.99616 m),
  0 outside the mesh. **Its control** (mesh moved 30 m east) still fails:
  847,428 interior and 1,888 strip nodes over.
- **0 constrained-Delaunay violations** with the seams and the outline as
  constraints: 10,795 constraint edges, 1,059,582 edges checked, 7,857 by
  the exact test. Merged mesh: 359,290 points, 711,088 triangles; worst
  angle 0.0283°, max degree 17.
- **Time and peak memory (recorded, not judged):** 6.20 s wall, peak
  footprint 1.64 GB, maximum RSS 1.69 GB (`/usr/bin/time -l`,
  `velhas.err`). The 51 refine calls sum to 0.74 s
  (`velhas_refine_calls.jsonl`).
- **Start-pass inserts:** 14,951 over the 51 pieces; median 321, range 0
  to 566.
- Against 2026-10-04 (different code, not judged): 361,185 points and
  19,866 start inserts then, 359,290 and 14,951 now. The cause is not
  measured here; master changed between the two runs.

## Method

- **Bench:** `pairs.sh` (log `pairs.log`): one balanced batch, base, head,
  head, base; 5 repeats, threads default and 1 to 20, domains `tile` and
  `quarter` (the benchmark's DTM10 tile at 1 m), `--hardening on`,
  `caffeinate -i`, stopping if off AC or if swap grew. Each run is built
  fresh by `bench.py` into its tree's `build-bench`. bench.py blob
  `16ad2ff5…` on both sides. Run directories: `2026-10-18/{base,head}-r{1,2}/`.
  `base-r1` reports NO BASELINE (exit 2), as it is the baseline. The head
  runs are marked dirty: the only untracked files were this directory's.
  `summarize.py` writes `tables.md`.
- **Smoke:** `smoke.sh` (log `smoke.log`). `smoke_driver.py` runs the head's
  Release package (`build-bench/pkg`, built by the head bench runs) with
  `MallocLargeCache=0` under `/usr/bin/time -l` and records every
  `cli.refine` call: `mesh --dem anadem-v1 --cache ~/.cache/rasputin --domain
  bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983 --tolerance
  10 --pieces 64 --ascii`. Checks: `check_pieces.py` (`check_pieces.json`)
  and `once_on_outline.py`, copied from the 2026-10-04 run, the latter's
  `strip_check` import pointed at the in-tree
  `docs/benchmarks/2026-10-04/15f-3-acceptance/`.
- **Interpreter for the checks:** the main checkout's Python 3.14 venv
  (matplotlib, which this worktree's 3.13 venv lacks), with
  `QUALITY_PKG` = this worktree's `src_python`: `tin_engine.stats.quality`
  is numpy only, so no extension is involved. The piece files are this
  branch's.
- **Not kept:** the piece files (ASCII), the merged mesh, the quality
  meshes and the base worktree were in the session scratchpad and are
  deleted. Rerun `smoke.sh` and `pairs.sh` to regenerate them.
