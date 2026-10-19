# 23c-2: bench and a real-data smoke run (@perf, 2026-10-04)

The `23c` worktree, HEAD `0c23255` (code `38ab1a0`). The rule is
`docs/increments/23-basin-scale.md`, "23c-2 at green (38ab1a0)", point 8.
AC power throughout: 100 %, charged, `pmset -g batt` before and after
every run. No run added to swap.

## Verdict

**PASS on both parts.**

1. **Bench (the uncut path).**
   - Mesh hashes are unchanged against the base `aacf37c` on both domains:
     tile `11741a81…`, quarter `ccebf96a…`. Quality is identical: 0.6296° /
     0.3955°, max degree 74 / 18, 0 Delaunay violations.
   - Refine time is -1.3 to +1.0 % over every cell of the sweep.
     Process time is -0.9 to +1.9 %.
   - At 1 / 20 threads: tile refine -0.3 / -0.8 %, process +0.4 / +0.9 %;
     quarter refine -0.3 / -0.4 %, process +0.7 / +1.3 %.
   - One balanced batch: base, 23c-2, 23c-2, base. That is 2 runs and 10
     samples a cell per side (`23c-2-smoke/tables.md`).
2. **The Velhas smoke, ANADEM from the cache, EPSG:31983, 10 m,
   `--pieces 64`.**
   - **Completion.** The run completed (exit 0). The 6 × 11 partition gives
     **51 piece files** and `index.json` (the other cells are empty). The
     cut run writes nothing at `--out`, as designed.
   - **Conformity holds over the piece files.**
     - 3,305 constraint segments appear in exactly two pieces, with
       identical end points and identical z (0 mismatches). No segment
       appears in more than two.
     - The 7,428 segments that only one piece has all lie on the domain's
       outline: their midpoints are within 0.66 mm of the outline projected
       into EPSG:31983 (`once_on_outline.txt`).
     - Points the pieces share have identical z (largest difference 0.0).
   - **The independent check finds 0 over tolerance.** It is 15c-2's
     `run_geo.independent_check`, run on the merged pieces. It covers
     13,791,972 interior source nodes (largest error 9.99998 m) and 27,806
     strip nodes (9.99616 m), with 0 outside the mesh. Its control, the
     mesh moved 30 m east, finds 844,330 and 1,873 over.
   - **What the check does not cover.** It checks the source (ANADEM)
     nodes. The resampled-grid nodes were not checked separately. The
     run's own figure for them is in each piece's report line: largest
     error 9.99991 m over the 51 pieces.
   - **0 constrained-Delaunay violations** with the seams and the outline
     as constraints: 10,733 constraint edges, 1,065,391 edges checked, 8,147
     of them by the exact test. Merged mesh: 361,185 points and 714,940
     triangles. Worst angle 0.0283°, max degree 17.
   - **Time and peak memory (recorded, not judged):** 6.17 s wall, peak
     footprint 1.64 GB (`/usr/bin/time -l`). The per-piece refine calls
     sum to 0.39 s (`velhas_refine_calls.jsonl`).
   - **Start-pass inserts per piece:** 19,866 over the 51 pieces; median
     439, range 0 to 727 (`velhas_pieces_stderr.txt`, one line per piece,
     and `velhas_refine_calls.jsonl`).

## Method

- **Bench.**
  - Builds: Release `_core` from `bench.build()` with this worktree's
    `.venv` (CPython 3.13): base `5815b450…` in a detached worktree, 23c-2
    `1299a76d…`.
  - Runs: `23c-2-smoke/pairs.sh` (log `pairs.log`), 5 repeats, threads 0
    and 1 to 20, both domains, `caffeinate -i`. Run directories:
    `23c2-base-r{1,2}/` and `23c2-r{1,2}/`.
- **The smoke run.**
  - Command: `smoke_driver.py` runs this branch's package with
    `MallocLargeCache=0` under `/usr/bin/time -l` and records each
    `cli.refine` call: `mesh --dem anadem-v1 --cache … --domain
    bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983
    --tolerance 10 --pieces 64 --ascii`. Its output is `velhas.{out,err,pmset}`.
  - Checks: `check_pieces.py` (→ `check_pieces.json`) does conformity,
    the merge, `bench.quality()` and the independent check. Its merged mesh
    is a measurement artefact, not 23d's stitcher. `once_on_outline.py`
    covers the unshared segments.
  - Interpreter: the checks need matplotlib, which this worktree's 3.13
    venv lacks. So they ran on 15f-3's 3.14 venv, with `QUALITY_PKG`
    pointing at master's 3.14 Release package for `tin_engine.stats.quality`
    only. The piece files themselves are this branch's.
- **Memory at the start.** The compressor held about 12-14 GB of RAM
  (Safari's WebKit processes), and swap stood at 10.37 of 11.26 GB. The
  smoke run's 1.64 GB fit.
- The piece files (ASCII), the merged mesh and the base worktree were in the
  scratchpad and are deleted. To regenerate them, rerun the command above.

## Log
- 2026-10-04T04:34:57Z: AC 100 %; compressor ~13.9 GB occupied, swap 10.45 of 11.26 GB.
