# 23c-4 acceptance: cuts on a cell size not exact in binary (@perf, 2026-10-19)

`worktree-23c-cellsize` at `e970bde4` (code review approved on `0adf1f46`;
the later commits are prose). Base: 23c-3, `worktree-23c-3` at `bb5b2fcc`
(no code differs from `09c8aad4`, accepted in `../23c-3/`), run from a
detached scratch worktree. The rule is `docs/increments/23-basin-scale.md`,
"23c-4", "`@perf` applies": bench hashes unchanged (K1), and the Velhas
piece's files at 10 m equal 23c-3's (K6). No timing claim is asked; the
brief adds short runs at 9.7 m and 9.9 m.

AC power, 100 %, charged, for every run (`pmset -g batt` in each `run.json`,
`pairs.log`, `logs/*.pmset`, `cellsize.pmset`). Swap 6,006.81 MB used of
7.0 GB at the start and the end. Apple M1 Max, 8 P + 2 E cores, 32 GiB,
macOS 27.0.1. Nothing of mine ran beside a timing run.

## Verdict

**ACCEPTED.** Every whole-metre output is byte-identical to 23c-3's, and
every 9.7 m and 9.9 m cut that 23c-3 refused now runs with no oracle finding.

## 1. Bench (`tools/bench.py`, the uncut path), 02:42-02:50 UTC

- **K1: mesh hashes unchanged** in all three runs (head-r1, base-r1,
  head-r2) and equal to 23c-3's accepted `head-r4`: tile `11741a81…`,
  quarter `a60fb597…`. Quality the same: worst angle 0.6296° / 0.3955°,
  max degree 74 / 18, within tolerance, 0 Delaunay violations of 692,056 /
  641,792 edges.
- Same `_core` (sha256 `52390c7f…`, Release, bounds checks `libc++ fast`,
  bench.py's `--hardening on`) on both trees and in 23c-3's run; bench.py
  blob `16ad2ff5…` on both. `git diff bb5b2fcc e970bde4 -- include src
  bindings CMakeLists.txt tools` is empty.
- bench.py's verdicts: head-r1 against 23c-3's `head-r4` (an earlier
  session) REGRESSION on one cell, quarter refine at 12 threads 0.1658 ->
  0.1744 s (+5.2 %); base-r1 against head-r1 ACCEPTED; head-r2 against
  base-r1, back to back, ACCEPTED (quarter t=12 0.1668 -> 0.1680 s). The
  branch changes no uncut-path code; the one cell is taken as noise, not
  measured further.
- `head-r4`'s tile is in 23c-3's slow process layout (`../23c-3/README.md`,
  part 1); this batch's tile runs are all in the fast one.

| run | tile t0 | tile t1 | tile t12 | tile ceiling | quarter t0 | quarter t1 | quarter t12 | quarter ceiling |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 23c-3 head-r4 (accepted) | 0.2478 | 0.5231 | 0.2369 | 2.11x | 0.1679 | 0.4220 | 0.1658 | 2.51x |
| 23c-4 head-r1 | 0.1934 | 0.4887 | 0.1939 | 2.53x | 0.1739 | 0.4191 | 0.1744 | 2.41x |
| 23c-3 base-r1 | 0.1928 | 0.4791 | 0.1904 | 2.48x | 0.1697 | 0.4270 | 0.1668 | 2.52x |
| 23c-4 head-r2 | 0.1913 | 0.4767 | 0.1914 | 2.49x | 0.1668 | 0.4217 | 0.1680 | 2.53x |

Refine medians in seconds; t0 is the core's default thread count; the
ceiling is t1 / t0. Every thread count is in each run's `README.md` under
`2026-10-19/`.

## 2. The Velhas piece: ANADEM from the cache, EPSG:31983, 10 m, 02:51 UTC

`cmp` of every file in `velhas.vtk.pieces/`, `index.json` included, 23c-3
(`bb5b2fcc`) against 23c-4, run back to back (`cmp.txt`):

| `--pieces` | files | differ | wall 23c-3 | wall 23c-4 | peak 23c-3 | peak 23c-4 |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 52 | 0 | 8.88 s | 8.85 s | 154.7 MB | 155.2 MB |
| 16 | 14 | 0 | 5.55 s | 5.51 s | 406.9 MB | 407.0 MB |
| 256 | 155 | 0 | 9.44 s | 9.47 s | 111.1 MB | 111.2 MB |

Each `index.json` is also byte-identical to the one 23c-3's acceptance
kept (`../23c-3/logs/head_t10_p*.index.json`, from `09c8aad4`), so the
piece digests equal 23c-3's on the same inputs.

## 3. Cuts at 9.7 m and 9.9 m, 02:52-02:57 UTC

`cellsize_runs.py`, one process per tree, each tree's Release package, the
oracles of 23c-4's `tests/python/pieces_fixtures.py` for both: conformity,
one z per (x, y), a valid union against the noded outline, no overlap,
exact constrained Delaunay, the tolerance; seams within 1e-6 m of a lattice
line at 9.7 and 9.9 m, on an exact lattice index at 10 m. `square` is the
design's scene (41 × 41 nodes); `big` the same hills at 201 × 201 nodes and
the square scaled by 5. `--tolerance 1`, default threads. Wall is the
in-process `rasputin mesh` call; the oracle time is separate.

| scene | h (m) | `--pieces` | 23c-3 exit | 23c-4 exit | 23c-4 wall s | piece files | triangles | oracle findings | oracle s | bytes vs 23c-3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| square | 9.7 | 4 | 2 | 0 | 0.068 | 4 | 1369 | 0 | 0.223 | 23c-3 wrote none |
| square | 9.7 | 16 | 2 | 0 | 0.051 | 16 | 1260 | 0 | 0.199 | 23c-3 wrote none |
| square | 9.9 | 4 | 2 | 0 | 0.04 | 4 | 1302 | 0 | 0.204 | 23c-3 wrote none |
| square | 9.9 | 16 | 0 | 0 | 0.05 | 16 | 1226 | 0 | 0.193 | identical |
| square | 10.0 | 4 | 0 | 0 | 0.038 | 4 | 1163 | 0 | 0.154 | identical |
| square | 10.0 | 16 | 0 | 0 | 0.048 | 16 | 1163 | 0 | 0.152 | identical |
| big | 9.7 | 4 | 2 | 0 | 0.104 | 4 | 34781 | 0 | 15.053 | 23c-3 wrote none |
| big | 9.7 | 16 | 2 | 0 | 0.134 | 16 | 34801 | 0 | 15.024 | 23c-3 wrote none |
| big | 9.7 | 64 | 2 | 0 | 0.178 | 64 | 34797 | 0 | 15.122 | 23c-3 wrote none |
| big | 9.9 | 4 | 2 | 0 | 0.127 | 4 | 33452 | 0 | 13.987 | 23c-3 wrote none |
| big | 9.9 | 16 | 2 | 0 | 0.136 | 16 | 33456 | 0 | 14.155 | 23c-3 wrote none |
| big | 9.9 | 64 | 2 | 0 | 0.171 | 64 | 33215 | 0 | 13.973 | 23c-3 wrote none |
| big | 10.0 | 4 | 0 | 0 | 0.068 | 4 | 30967 | 0 | 12.052 | identical |
| big | 10.0 | 16 | 0 | 0 | 0.092 | 16 | 30898 | 0 | 11.943 | identical |
| big | 10.0 | 64 | 0 | 0 | 0.138 | 64 | 30882 | 0 | 11.929 | identical |

23c-3's exit 2 is the strip refusal the design names ("a seam strip, rows
0-20, columns 17-20, is outside its piece" at square 9.7 m, 4 pieces).
"Bytes" is one sha256 over every file of `<out>.pieces`. Square 9.9 m with
16 pieces, which 23c-3 already cut, is unchanged byte for byte. The scale is
small; nothing here is a basin-scale timing.

## Method

- **Bench:** `pairs.sh` (W B W), 5 repeats, threads default and 1 to 20,
  domains `tile` and `quarter`, `--hardening on`, `caffeinate -i`; each tree
  built fresh by bench.py into its own `build-bench`. The head runs are
  marked dirty: the only untracked files were this directory's.
  `summarize.py` writes `tables.md`.
- **Velhas:** `velhas_all.sh` runs `velhas.sh` (23c-3's, the venv path
  changed) for each case through `velhas_driver.py` (23c-3's, unchanged):
  `mesh --dem anadem-v1 --cache ~/.cache/rasputin --domain
  bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983
  --tolerance 10 --pieces P --binary`, under `/usr/bin/time -l`.
- **Decimal cells:** `cellsize_runs.py PKG TESTS RECORD SCRATCH`, records in
  `cellsize_{head,base}.jsonl`, table in `cellsize_table.md`.
- **Sanitizers:** not applicable; no C++ differs from the base.
- **Not kept:** the meshes and piece files (session scratchpad, deleted) and
  the scratch worktree. Rerun `pairs.sh W B S`, `velhas_all.sh W B S` and
  `cellsize_runs.py` with a worktree of `bb5b2fcc` as B to regenerate them.
