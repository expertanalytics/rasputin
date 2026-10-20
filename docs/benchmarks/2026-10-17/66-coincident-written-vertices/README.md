# Increment 66 acceptance: coincident written vertices (@perf, 2026-10-17)

Branch `worktree-coincident-vertices`, head `40a4ddce` (code review round 2
approved), against its merge base, master `61755c04` (PR #282). The rules:
`docs/increments/66-coincident-written-vertices.md`, section 5;
`docs/increments/README.md`, "Acceptance".

## Verdict: ACCEPTED

| gate | measured | limit | result |
|---|---|---|---|
| mesh sha256, tile (ASCII quality run) | base `11741a81…0e005c`, head `11741a81…0e005c` | unchanged | met |
| mesh sha256, quarter (ASCII quality run) | base `a60fb597…bb60a34`, head `a60fb597…bb60a34` | unchanged | met |
| mesh sha256, `--binary` record run | tile `f5261521…8e6a7c` both sides; quarter `de5bde17…83702e` both sides | unchanged | met |
| Delaunay violations (`bench.py`'s check) | tile 0 of 692,056 edges, quarter 0 of 641,792, both sides | no loss | met |
| `edges_not_delaunay_as_written` (run record) | tile 0 / 0, quarter 0 / 0 (base / head) | no loss | met |
| `triangles_flat_as_written` (run record) | tile 0 / 0, quarter 0 / 0 (base / head) | no loss | met |
| 1 m benchmark and thread sweep | 42 cells (2 domains x 21 thread counts): median time change -2.4 % to +4.2 %, median of the changes +0.7 % | 5 % (`bench.py`'s band) | `ACCEPTED` by `bench.py` |
| quality | worst angle 0.6296° / 0.3955°, max degree 74 / 18 (tile / quarter), within tolerance: the same on both sides | no loss | met |

Section 5 expected both hashes unchanged; they are, and they are the same
meshes as increment 65's acceptance (`../65-delaunay-as-written`).

## Method

`tools/bench.py run`, default flags (DEM `7908_3_10m_z33.tif`, tolerance 1,
tile and quarter domains, threads 0 and 1 to 20, 5 repeats interleaved,
hardening on: libc++ fast mode on both sides), Release builds into each
tree's `build-bench`. Base then head, back to back: `bench-base/` is a
scratch clone of `61755c04` measured with `--tree` (it printed
`NO BASELINE`, no stored run for that commit), and `bench-head/` is this
worktree, judged with `--baseline bench-base`. Both trees clean; the same
`bench.py` (blob `16ad2ff5`) ran both.

Apple M1 Max, 8 P + 2 E cores, 32 GiB, macOS 27.0.1. Power: AC, 100 %, read
before and after each run (raw `pmset` text in each `run.json`) and before
and after the record runs. Base started 15:31, head 15:34. `_core` sha256:
base `461eec90…`, head `9cc0ea4d…`.

`bench.py` does not report the run record's `triangles_flat_as_written` or
`edges_not_delaunay_as_written`, so `records/` holds one extra run per side
and domain, at the CLI's default thread count, through the same child and
each side's `build-bench/pkg`:
`bench.py _child --pkg <tree>/build-bench/pkg --threads 0 -- mesh --dem <DEM>
--tolerance 1 --out <path> --binary [--domain docs/benchmarks/2026-09-26/quarter.geojson]
--record records/<side>-<domain>.json`.

`delta.md` is every cell (`scripts/delta.py`, increment 65's script, copied
unchanged). Ceiling, speed-up from 1 thread to 20 (best in brackets): tile
2.48x (2.58x at 12) base, 2.47x (2.57x at 10) head; quarter 2.46x (2.52x at
14) base, 2.43x (2.51x at 15) head.

| domain | threads | base median s | head median s | change |
|---|---:|---:|---:|---:|
| tile | 0 (default) | 0.1913 | 0.1913 | +0.0 % |
| tile | 1 | 0.4761 | 0.4783 | +0.5 % |
| tile | 20 | 0.1921 | 0.1940 | +1.0 % |
| quarter | 0 (default) | 0.1724 | 0.1696 | -1.6 % |
| quarter | 1 | 0.4162 | 0.4184 | +0.5 % |
| quarter | 20 | 0.1695 | 0.1719 | +1.4 % |

No profile was taken, so no cause is claimed for the scatter.

## Meshes

Not committed. They were in the session scratchpad (`meshes-base/`,
`meshes-head/`, `rec/`) and are gone with it. Regenerate: rerun the
commands above.
