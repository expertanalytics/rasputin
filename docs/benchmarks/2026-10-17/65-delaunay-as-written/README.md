# Increment 65 acceptance: the last Delaunay check uses the written points (@perf, 2026-10-17)

Branch `worktree-refine-world`, head `5d36b15d` (code review round 2
approved), against its merge base, master `db12c3e1` (PR #276). The rules:
`docs/increments/65-delaunay-as-written.md`, section 5;
`docs/increments/README.md`, "Acceptance".

## Verdict: ACCEPTED

| gate | measured | limit | result |
|---|---|---|---|
| mesh sha256, tile | base `11741a81…0e005c`, head `11741a81…0e005c` | unchanged (base Delaunay count 0) | met |
| mesh sha256, quarter | base `a60fb597…bb60a34`, head `a60fb597…bb60a34` | unchanged (base Delaunay count 0) | met |
| Delaunay violations | tile 0 of 692,056 edges, quarter 0 of 641,792, on both sides | no loss | met |
| 1 m benchmark and thread sweep | 42 cells (2 domains x 21 thread counts): median time change -4.3 % to +3.7 %, median of the changes -0.4 % | 5 % (`bench.py`'s band) | `ACCEPTED` by `bench.py` |
| quality | worst angle 0.6296° / 0.3955°, max degree 74 / 18, within tolerance: the same on both sides | no loss | met |

Section 5 expected both hashes unchanged, since the base count is 0 on both
domains; they are. The meshes are the same as increment 34's acceptance
(`2026-10-10/slope-tol`).

## Method

`tools/bench.py run`, default flags (DEM `7908_3_10m_z33.tif`, tolerance 1,
tile and quarter domains, threads 0 and 1 to 20, 5 repeats interleaved,
hardening on: libc++ fast mode on both sides), Release builds into each
tree's `build-bench`. Base then head, back to back: no stored run for the
previous merge commit exists (`bench-base` printed `NO BASELINE`), so
`bench-base/` is a scratch clone of `db12c3e1` measured with `--tree`, and
`bench-head/` is this worktree, judged with `--baseline bench-base`. Both
trees clean; the same `bench.py` (blob `16ad2ff5`) ran both.

Apple M1 Max, 8 P + 2 E cores, 32 GiB, macOS 27.0.1. Power: AC, 100 %, read
before and after each run (raw `pmset` text in each `run.json`). Base
started 05:32, head 05:35. `_core` sha256: base `d54e0d40…`, head
`461eec90…`.

`delta.md` is every cell (`scripts/delta.py`, the 2026-10-10 script, copied
unchanged). Ceiling, speed-up from 1 thread to 20 (best in brackets): tile
2.62x (2.62x at 20) base, 2.58x (2.69x at 13) head; quarter 2.55x (2.62x at
16) base, 2.66x (2.68x at 11) head.

| domain | threads | base median s | head median s | change |
|---|---:|---:|---:|---:|
| tile | 0 (default) | 0.1859 | 0.1880 | +1.2 % |
| tile | 1 | 0.4923 | 0.4975 | +1.1 % |
| tile | 20 | 0.1881 | 0.1929 | +2.6 % |
| quarter | 0 (default) | 0.1660 | 0.1715 | +3.4 % |
| quarter | 1 | 0.4343 | 0.4421 | +1.8 % |
| quarter | 20 | 0.1701 | 0.1661 | -2.4 % |

The changes are larger than section 5's 0.5 % estimate of the check's cost
but scatter both ways across thread counts (the 1-thread spread within one
side is 0.426 to 0.455 s on quarter); no profile was taken, so no cause is
claimed for them.

## Meshes

Not committed. They were in the session scratchpad (`meshes-base/`,
`meshes-head/`) and are gone with it. Regenerate: rerun the commands above,
or one child per domain with `--ascii --out PATH`, as `bench-head/README.md`
states.
