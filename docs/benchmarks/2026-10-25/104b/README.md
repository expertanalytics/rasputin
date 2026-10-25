# 104b acceptance: the error summary, `--rms` (@perf, 2026-10-25)

Head `worktree-tri-budget-rms` at `0d751480`. Base: master `d886847d`, the
branch's fork point (`git diff d886847d origin/master -- include src_python
bindings CMakeLists.txt tools` is empty), run from a detached scratch
worktree. The rule is `docs/increments/104-triangle-budget.md`, section 11
(104b): the 1 m benchmark's mesh hashes unchanged; without `--rms` the whole
`rasputin mesh` time within the quick check's noise; with `--rms` the
"error summary" row within a factor 1.5 of section 4.3's measured pass and
the run total grown by that row and no more; on Kultsjön the budget run pays
one pass, under 0.08 s.

Apple M1 Max (8 P + 2 E cores), 32 GiB, macOS 27.0.1. `tools/bench.py`
blob `16ad2ff5`, `--hardening on` (bounds checks `libc++ fast`). Every run
on AC (`pmset -g batt`: "Now drawing from 'AC Power'", 89-97 %, charging),
read before and after each run and kept in each `run.json`, `runs.jsonl`
line and quick-check log. Nothing of mine ran beside a timing run; another
Claude session's process was running throughout.

**Sanitizer run first.** A scratch `_core` of the head built with
`-fsanitize=undefined -fno-sanitize-recover=all` and libc++ extensive
hardening (`RelWithDebInfo`) passed `tests/python/test_core_error_summary.py`
and `tests/python/test_cli_rms.py` (53 passed) before any timing.

## Verdict

**ACCEPTED** on the gated measures: every mesh is byte-identical to
master's, the run without `--rms` is within noise of master on every quick
check case, and the Kultsjön budget run pays one pass of 0.05 s (under
0.08 s), its total 0.046 s over master's.

**Not met: section 11's factor 1.5 on the `--rms` row.** The row is 1.6 to
3.8 times section 4.3's measured pass on five of the seven quick check
cases (table in section 3). Profiled cause: the pass builds the mesh's
adjacency (`LatticeMesh::build` in `summary_lattice`, `measure.hpp:85`)
serially before the parallel scan; section 4.3 timed the scan alone. The
`--rms` option is off by default, so its time is reported here, not gated
(`.claude/agents/perf.md`, section 4); the miss is handed back for a ruling.

## 1. Bench (`tools/bench.py`, the uncut path, 1 m tolerance)

| run | power | commit | `_core` sha256 | tile hash | quarter hash | tile t0 | tile t1 | tile t12 | tile t20 | tile ceiling | quarter t0 | quarter t1 | quarter t12 | quarter ceiling | bench.py verdict |
|---|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| head-ac-r1 | AC 92 % | 0d751480 | a353ac1e | 11741a81 | a60fb597 | 0.1868 | 0.4795 | 0.1908 | 0.1891 | 2.54x | 0.1660 | 0.4207 | 0.1666 | 2.39x | ACCEPTED |
| base-ac-r1 | AC 93 % | d886847d | 52390c7f | 11741a81 | a60fb597 | 0.1890 | 0.4769 | 0.1937 | 0.1924 | 2.48x | 0.1661 | 0.4220 | 0.1676 | 2.42x | ACCEPTED |
| head-ac-r2 | AC 94 % | 0d751480 | a353ac1e | 11741a81 | a60fb597 | 0.1942 | 0.4770 | 0.1880 | 0.1947 | 2.45x | 0.1656 | 0.4223 | 0.1660 | 2.42x | ACCEPTED |

Refine medians in seconds, 5 repeats; t0 is the core's default thread
count; the ceiling is t1 over t20 (best: tile 2.55x, quarter 2.53-2.54x).
`bench.py compare` of each head run against `base-ac-r1`, and against
104a's `head-ac-r2`, is ACCEPTED (5 % threshold). Head-ac-r2 is marked
"dirty" for this folder, which was untracked when it ran.

- **Mesh hashes unchanged**: tile `11741a81…`, quarter `a60fb597…`, as in
  104a and 23c-4. Quality identical: worst angle 0.6296° / 0.3955°, max
  degree 74 / 18, within tolerance, 0 Delaunay violations of 692,056 /
  641,792 edges.
- The `_core` differs from master's (`a353ac1e` against `52390c7f`): the
  accumulator template is new code; `refine`'s instantiation gives the same
  meshes in the same time.

## 2. The quick check without `--rms` (`tools/bench_quick.py`)

The stored AC baseline (2026-10-19) was refused, "NO BASELINE: machine
os_version: 27.0.1 vs 27.0", so master was measured back to back with
`--save-baseline` and the head twice against it. The cases file named an
option renamed since (`--features-map`, now `--features-legend`), which
made the numedalslagen and lagan cases exit 2 on master and head alike;
`docs/benchmarks/quick/cases.toml` is corrected in this commit, and the
stored baseline is now master `d886847d` on macOS 27.0.1
(`quick/baseline-master-d886847d.json` is a copy).

| case | master total median (min-max) s | head, without `--rms` (`rms_runs.py`, 5 runs) | quick check, head r1 | head r2 |
|---|---|---|---|---|
| tile | 0.4387 (0.4229-0.4765) | 0.420 | no change | no change |
| tile, 1 thread | 0.7612 (0.7160-0.7692) | 0.702 | FASTER, 0.76 → 0.71 s (−7.3 %, band 7.0 %) | FASTER, −7.7 % |
| quarter | 0.3641 (0.3557-0.3697) | 0.358 | no change | no change |
| numedalslagen | 8.8292 (8.8169-8.8973) | 8.873 | no change | no change |
| lagan | 14.5337 (14.5214-14.6235) | 14.622 | no change | no change |
| geilo-al-ramp | 0.3362 (0.3346-0.3374) | 0.338 | no change | no change |
| romsdal-slope | 3.5794 (3.5759-3.6510) | 3.636 | no change | no change |
| velhas | 0.1276 (0.1263-0.1280) | 0.128 | no change | no change |

Both head runs: exit 0, no `SLOWER`, no `MESH CHANGED` line, so every
case's mesh sha256 equals master's. The tile one-thread `FASTER` is thought
to be master's first-run spread (its range is 7 % wide), not profiled.

## 3. With `--rms`: the error summary row

`rms_runs.py`: each quick check case with and without `--rms`, 5 repeats
interleaved, head only; `--stats` rows are rounded to 1 ms.

| case (threads) | section 4.3's pass (s) | "error summary" row, median (min-max) s | row / pass | total without (s) | total with (s) | grown by | nodes measured | RMS (m) |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| tile (all) | 0.0265 | 0.057 (0.056-0.060) | 2.2 | 0.420 | 0.474 | 0.054 | 25,233,600 | 0.0675 |
| tile (1) | 0.1567 | 0.219 (0.218-0.219) | 1.4 | 0.702 | 0.918 | 0.216 | 25,233,600 | 0.0675 |
| quarter (all) | 0.0160 | 0.052 (0.052-0.052) | 3.3 | 0.358 | 0.408 | 0.050 | 6,857,099 | 0.1333 |
| velhas (all) | 0.0038 | 0.009 (0.009-0.010) | 2.4 | 0.128 | 0.136 | 0.008 | 54,411 | 0.2471 |
| romsdal-slope (all) | 0.1201 | 0.456 (0.452-0.463) | 3.8 | 3.636 | 4.036 | 0.400 | 23,963,316 | 2.7306 |
| numedalslagen (all) | 0.1144 | 0.185 (0.184-0.187) | 1.6 | 8.873 | 9.028 | 0.155 | 55,214,155 | 2.9553 |
| lagan (all) | 0.1643 | 0.195 (0.193-0.198) | 1.2 | 14.622 | 14.811 | 0.189 | 8,202,663 | 3.2862 |
| geilo-al-ramp (all) | | 0.008 (0.008-0.008) | | 0.338 | 0.344 | 0.006 | 2,637,728 | 4.5713 |

- **The total grows by the row and no more**: on every case the growth is
  at or under the row's median.
- **`--rms` changes no mesh**: each case's mesh sha256 is the same with and
  without it, and the largest error is unchanged.
- **The cause of the miss, profiled** (`profile/`, one `sample` run of
  romsdal-slope with `--rms` on an unstripped `-O3` build of the head, row
  0.539 s): of 379 samples in the binding, 258 are the serial lattice
  build (`summary_lattice` → `LatticeMesh::build`, its edge sort 50 of
  them) and 110 the parallel scan (`summarise`, 106 of them the main thread
  joining workers). The scan's share matches section 4.3's 0.120 s; the
  lattice build is what section 4.3 did not time.

### The pass against thread count (tile, `--rms --threads N`, 5 repeats)

| threads | 1 | 2 | 4 | 8 | 10 |
|---|---:|---:|---:|---:|---:|
| "error summary" row, median (s) | 0.218 | 0.129 | 0.087 | 0.058 | 0.056 |
| speed-up over 1 thread | 1.0 | 1.7 | 2.5 | 3.8 | 3.9 |
| run total, median (s) | 0.918 | 0.653 | 0.532 | 0.472 | 0.478 |

The flattening from 8 threads is consistent with the serial lattice build
above; the sums are bit-identical at every count (RMS 0.06749796106317188
on each run).

## 4. Kultsjön: `--tolerance 10`, `--rms` and `--triangles 65893` (AC)

Section 7's command (GLO-30, catchment `kultsjon.geojson`, CORINE
features, 4 threads), to `.gpkg` with `--stats` and `--record`;
`kultsjon_runs.py`, 5 repeats interleaved over five runs, then one `.vtk`
each. Wall is the whole process; total is `--stats`' total. All 30 runs on
AC (97 %, charging).

| run | wall median (min-max) s | `--stats` total median (min-max) s | error summary row (s) | tolerance (m) | largest error (m) | RMS (m) | nodes measured | `.vtk` sha256 |
|---|---|---|---|---:|---:|---:|---:|---|
| master `--tolerance 10` | 3.092 (3.070-3.305) | 2.809 (2.790-3.020) | | 10 | 9.99984 | | | `a24ae336…` |
| head `--tolerance 10` | 3.081 (3.068-3.094) | 2.800 (2.789-2.805) | | 10 | 9.99984 | | | `a24ae336…` |
| head `--tolerance 10 --rms` | 3.141 (3.121-3.150) | 2.857 (2.837-2.870) | 0.049-0.050 | 10 | 9.99984 | 3.0511 | 2,106,535 | `a24ae336…` |
| master `--triangles 65893` | 3.867 (3.840-4.102) | 3.583 (3.560-3.817) | | 9.7885 | 9.78846 | | | `6c77a035…` |
| head `--triangles 65893` | 3.910 (3.887-3.918) | 3.629 (3.606-3.636) | 0.049-0.050 | 9.7885 | 9.78846 | 3.0032 | 2,106,399 | `6c77a035…` |

- **Byte-identical**: every `.vtk` equals master's and 104a's. The
  GeoPackages of master and head at 10 m differ only in the `run` table's
  `command` row (the output path). With `--rms` or a budget the `run` table
  gains `rms_error_m` and `error_nodes_measured` and nothing else differs.
- **The budget run pays one pass**: 0.05 s, under 0.08 s; its total is
  0.046 s over master's budget run.
- **P1's RMS met**: 3.003 m, within 2.98-3.02 m. The 10 m mesh's 3.051 m
  equals section 7's probe (3.05 m).

## Files and regeneration

- `../104b-*/`: bench.py's evidence. Meshes were in
  `$TMPDIR/rasputin-bench/<label>`, not kept; `python tools/bench.py run
  --label <label> --tree <tree>` regenerates them.
- `quick/`: the three quick check logs and the saved master baseline.
  `RASPUTIN_DATA` must be set in a worktree (the default resolves beside
  the worktree, not the main checkout).
- `rms_runs.py` (`rms/runs.jsonl`, `rms/stats/`), `kultsjon_runs.py`
  (`kultsjon/`): the drivers; they need the tree's `build-bench/pkg` from a
  bench.py run. Meshes went to `/tmp/perf-104b-meshes/`, not kept.
- `profile/`: the call-tree excerpt and the `--stats` of the profiled run.
  The full `sample` output (1.9 MB) is not kept: build `_core` with
  `-DCMAKE_BUILD_TYPE=RelWithDebInfo "-DCMAKE_CXX_FLAGS_RELWITHDEBINFO=-O3
  -DNDEBUG -g"` and run `sample <pid> 8 1` on the romsdal-slope case with
  `--rms`.
