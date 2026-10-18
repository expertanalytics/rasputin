# Increment 67 P2 acceptance: C++ comment sweep and the `as_written` split (@perf, 2026-10-18)

Branch `worktree-audit-cuts-p2`, head `eba95909` (code review approved),
against its merge base, master `9c80b7de` (PR #284). The rules:
`docs/increments/67-audit-cuts.md`, section 5.3;
`docs/increments/README.md`, "Acceptance".

## Verdict: ACCEPTED

| gate | measured | limit | result |
|---|---|---|---|
| mesh sha256, tile (ASCII quality run) | base `11741a81…0e005c`, head `11741a81…0e005c` | unchanged | met |
| mesh sha256, quarter (ASCII quality run) | base `a60fb597…bb60a34`, head `a60fb597…bb60a34` | unchanged | met |
| mesh sha256, `--binary` record run | tile `f5261521ad9d9024…` both sides; quarter `de5bde173c216fe2…` both sides | unchanged | met |
| Delaunay violations (`bench.py`'s check) | tile 0 of 692,056 edges, quarter 0 of 641,792, both sides | no loss | met |
| `edges_not_delaunay_as_written` (run record) | tile 0 / 0, quarter 0 / 0 (base / head) | no loss | met |
| `triangles_flat_as_written` (run record) | tile 0 / 0, quarter 0 / 0 (base / head) | no loss | met |
| 1 m benchmark and thread sweep | 42 cells (2 domains x 21 thread counts): median time change -6.1 % to +2.4 %, median of the changes -0.4 % | 5 % slower (`bench.py`'s band) | `ACCEPTED` by `bench.py` |
| quality | worst angle 0.6296° / 0.3955°, max degree 74 / 18 (tile / quarter), within tolerance: the same on both sides | no loss | met |

Section 5.3 requires both hashes unchanged (a refactor); they are, and they
are the same meshes as increments 65 and 66 accepted
(`../../2026-10-17/66-coincident-written-vertices`).

## Method

`tools/bench.py run`, default flags (DEM `7908_3_10m_z33.tif`, tolerance 1,
tile and quarter domains, threads 0 and 1 to 20, 5 repeats interleaved,
hardening on: libc++ fast mode on both sides), Release builds into each
tree's `build-bench`. Base then head, back to back: `bench-base/` is a
scratch clone of `9c80b7de` measured with `--tree` (it printed
`NO BASELINE`, no stored run for that commit), and `bench-head/` is this
worktree, judged with `--baseline bench-base`. Both trees clean; the same
`bench.py` (blob `16ad2ff5`) ran both.

Apple M1 Max, 8 P + 2 E cores, 32 GiB, macOS 27.0.1. Power: AC, 100 %, read
before and after each run (raw `pmset` text in each `run.json`) and before
and after the record runs. Base started 02:24, head 02:27. `_core` sha256:
base `9cc0ea4d…` (the same as increment 66's head), head `52390c7f…`.

`records/` holds one extra run per side and domain at the CLI's default
thread count, for the run record's `as_written` counts, which `bench.py`
does not report:
`bench.py _child --pkg <tree>/build-bench/pkg --threads 0 -- mesh --dem <DEM>
--tolerance 1 --out <path> --binary [--domain docs/benchmarks/2026-09-26/quarter.geojson]
--record records/<side>-<domain>.json`.

`delta.md` is every cell (`scripts/delta.py`, copied unchanged from
increment 66). Ceiling, speed-up from 1 thread to 20 (best in brackets):
tile 2.65x (2.69x at 17) base, 2.52x (2.55x at 12) head; quarter 2.46x
(2.51x at 19) base, 2.42x (2.49x at 19) head.

| domain | threads | base median s | head median s | change |
|---|---:|---:|---:|---:|
| tile | 0 (default) | 0.1874 | 0.1908 | +1.8 % |
| tile | 1 | 0.5039 | 0.4732 | -6.1 % |
| tile | 20 | 0.1898 | 0.1879 | -1.0 % |
| quarter | 0 (default) | 0.1676 | 0.1716 | +2.4 % |
| quarter | 1 | 0.4198 | 0.4182 | -0.4 % |
| quarter | 20 | 0.1704 | 0.1725 | +1.2 % |

The base's tile 1-thread median (0.5039 s, samples 0.479 to 0.511 s) is
above increment 66's runs of the same code (0.476 to 0.478 s); that cell is
what lowers the base-to-head tile ceiling. No profile was taken, so no cause
is claimed for the scatter.

## Meshes

Not committed. They were in the session scratchpad (`meshes-base/`,
`meshes-head67/`, `rec67/`). Regenerate: rerun the commands above.
