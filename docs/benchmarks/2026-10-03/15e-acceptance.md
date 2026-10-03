# 15e acceptance: four memory fixes on the reprojected path (@perf, 2026-10-03)

**Verdict: ACCEPTED.** The 1 m benchmark and thread sweep do not change, and
the meshes are identical. On sub-basin 761 at 10 m, the peak drops from
14.3-14.5 GB to 6.9 GB with `MallocLargeCache=0`, and from 15.4-15.9 GB to
10.1 GB as shipped. The store is charged 16.4 B per point, down from
34-36 B. One measure is slower: the store's `freeze()` takes +0.26 to +0.39 s
(+7 to +11 %) as shipped, inside a run that is 6.7-8.7 s faster overall.

The increment is `docs/increments/15e-memory-fixes.md`, and its "Acceptance by
@perf" section is the brief. Branch `worktree-agent-ab601ec06c6508cb9` at
`78df916` (approved by `@reviewer` round 2), measured against its base
`7810cf8`. Master `390b516` differs from `7810cf8` only in governance
documents (`git diff --stat 7810cf8 390b516 -- include src src_python bindings tools/bench.py`
is empty), so `7810cf8` stands for master. The base was checked out with
`git worktree add --detach` in the session scratchpad and measured back to
back with the branch.

Machine: Apple M1 Max (8 P + 2 E cores), 32 GiB, macOS 27.0, Python 3.14.7,
numpy 2.5.3. **AC power, battery 100 % and charged, for every run.** The
`pmset -g batt` output from before and after each batch is in
`15e-acceptance/pmset_*.txt`, and each run's own record is in its `run.json`
or `runs/<tag>.json`. Builds: Release (`-O3 -DNDEBUG`, AppleClang 21), built
by `tools/bench.py` into each tree's `build-bench/`. The `_core` sha256 is
`41a22a6168db…` for the base (the same `.so` the 2026-10-02 basin-phases runs
used) and `a6309dc2b621…` for the branch. The bench records say the
branch's tree is `dirty`. The only cause is this evidence, untracked at
the time: `git status --porcelain` listed nothing but
`docs/benchmarks/2026-10-03/`.

## 1. The 1 m benchmark and thread sweep (`tools/bench.py`)

There are two back-to-back pairs, base then branch each time, with the
defaults: DEM `7908_3_10m_z33.tif`, tolerance 1 m, 5 repeats, 1 to 20
threads, and the tile and quarter domains. The evidence is in
`2026-10-03/15e-base-7810cf8/`, `15e/`, `15e-base-7810cf8-r2/` and
`15e-r2/`. `bench.py` printed ACCEPTED for `15e` against
`15e-base-7810cf8`, and for `15e-r2` against `15e-base-7810cf8-r2`.

`refine` median, seconds:

| domain | threads | base #1 | 15e #1 | base #2 | 15e #2 | Δ #1 | Δ #2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| tile | 0 (default) | 0.1953 | 0.1850 | 0.1857 | 0.1851 | -5.3 % | -0.3 % |
| tile | 1 | 0.5297 | 0.4822 | 0.4732 | 0.4714 | -9.0 % | -0.4 % |
| tile | 2 | 0.3383 | 0.3168 | 0.3019 | 0.3098 | -6.4 % | +2.6 % |
| tile | 4 | 0.2383 | 0.2282 | 0.2287 | 0.2273 | -4.3 % | -0.6 % |
| tile | 8 | 0.1958 | 0.1865 | 0.1848 | 0.1891 | -4.7 % | +2.3 % |
| tile | 10 | 0.1977 | 0.1852 | 0.1871 | 0.1872 | -6.4 % | +0.1 % |
| tile | 20 | 0.1949 | 0.1897 | 0.1876 | 0.1877 | -2.7 % | +0.0 % |
| quarter | 0 (default) | 0.1735 | 0.1660 | 0.1674 | 0.1648 | -4.3 % | -1.5 % |
| quarter | 1 | 0.4778 | 0.4126 | 0.4219 | 0.4148 | -13.6 % | -1.7 % |
| quarter | 2 | 0.2990 | 0.2680 | 0.2707 | 0.2707 | -10.4 % | -0.0 % |
| quarter | 4 | 0.2101 | 0.1992 | 0.2005 | 0.1995 | -5.2 % | -0.5 % |
| quarter | 8 | 0.1706 | 0.1639 | 0.1631 | 0.1641 | -4.0 % | +0.7 % |
| quarter | 10 | 0.1737 | 0.1634 | 0.1665 | 0.1635 | -6.0 % | -1.8 % |
| quarter | 20 | 0.1746 | 0.1665 | 0.1669 | 0.1658 | -4.6 % | -0.6 % |

Speed-up from 1 thread to 20 (best speed-up over 1 thread in brackets):

| domain | base #1 | 15e #1 | base #2 | 15e #2 |
|---|---|---|---|---|
| tile | 2.72 (2.74) | 2.54 (2.60) | 2.52 (2.57) | 2.51 (2.56) |
| quarter | 2.74 (2.80) | 2.48 (2.53) | 2.53 (2.59) | 2.50 (2.54) |

The second pair agrees to within ±2.6 %. In the first pair the base's run
was slow everywhere, by up to 13.6 % at 1 thread. That run was the
session's first, and against the stored 2026-10-02 runs `bench.py` judged
it a REGRESSION (`15e-base-7810cf8/README.md`). The same base commit then
ran within 1.7 % of the branch at 1 thread in pair 2. So pair 1's gap
belongs to that one run, not to the code. This was not profiled. The 1 m path does not use the
store or `resample` (the design says so). In all four runs the quality
results are the same and the mesh sha256 is the same for each domain:

- tile: worst angle 0.6296°, max degree 74, error within tolerance, 0
  Delaunay violations;
- quarter: worst angle 0.3955°, max degree 18, error within tolerance, 0
  Delaunay violations.

## 2. Memory phase by phase, sub-basin 761 at 10 m

The method and command are those of `2026-10-02/basin-phases/README.md`,
with the same domain (`761_epsg4674.geojson`), `--out-crs`, cache and
tolerance. Each of four configurations ran twice, interleaved (base and
branch, each as shipped and with `MallocLargeCache=0`): `15e-acceptance/run_all.sh`.
Each tree ran from a scratch venv whose only `.pth` puts that tree's
`build-bench/pkg` first, followed by the branch venv's site-packages. The
editable finder is not loaded, and the driver prints the `tin_engine` and
`_core` paths at start-up (`runs/*.log`).

The scripts are copied from `2026-10-02/basin-phases/`. Two of them were
changed:

- `phase_driver.py`. 15e's `resample` ends in `DemTile._adopt`, not the
  public constructor. The old wrap replaced `target_grid.DemTile` with a
  function and would break that call. The driver now wraps `_adopt` under
  the old marker name, so `analyse.py` reads both trees. It also records
  each tree's rows per block (256 for the base, `rows_per_block(27,786)` =
  37 for the branch).
- `analyse.py`. It reads the rows per block from that marker; before, the
  value was fixed at 256.

`run_phases.py`, `counts.py` and `probe_store.py` are unchanged.

All eight meshes are byte-identical (one sha256, `90a273d87ee3e784…`), with
2,940,024 triangles and 83,744 phase-2 insertions. Worst angle is 0.00584°,
max degree 17, and the achieved error is 9.9999994 m. So fix 1's smaller
blocks did not change a value (J3). The meshes were written to
`../rasputin_scratch/15e-acceptance/` and are not in the repository;
`run_all.sh` regenerates them.

Peak footprint per phase, in GB, rounds 1 / 2. "s" marks a sampled peak
(every 0.25 s); the other peaks are exact. The 2026-10-02 column is the
basin-phases README's 10 m column. Full tables: `15e-acceptance/tables.md`.

| phase | 2026-10-02 as shipped | base as shipped | 15e as shipped | 2026-10-02 `MLC=0` | base `MLC=0` | 15e `MLC=0` | design (live) |
|---|---:|---:|---:|---:|---:|---:|---:|
| P1 decode + assemble | 7.41 | 6.84 / 7.41 | 7.41 / 7.40 | 6.81 | 6.81 / 6.81 | 6.81 / 6.81 | 6.81 unchanged |
| P2 resample (blocks) | **17.83** | **15.36 / 15.91** | 8.76 / 8.82 | **13.97** | **14.26 / 14.49** | 5.32 s / 5.61 s | ~5.8 |
| P2 end (the old copy) | 10.19 s | 10.15 s / 10.35 s | 5.19 s / 5.18 s | 6.63 s | 6.64 s / 6.63 s | 4.48 s / 4.48 s | 4.4 |
| P3 refine phase 1 | 10.13 s | 7.64 s / 10.28 s | 5.19 s / 5.19 s | 4.70 s | 4.70 s / 4.75 s | 4.76 s / 4.76 s | |
| P4 store fill + freeze | 11.51 s | 11.11 s / 11.72 s | **9.76 / 9.75** | under-counted | under-counted | **6.93 / 6.92** | |
| P5 phase 2 | 12.16 s | 11.75 s / 12.36 s | **10.13 / 10.12** | under-counted | under-counted | 5.20 s / 5.20 s | |
| P6 trim + write | 7.84 s | 7.89 s / 7.80 s | 5.42 s / 5.41 s | | 2.82 s / 2.82 s | 4.72 s / 4.72 s | |
| `time -l` peak | 17.83 | 15.36 / 15.91 | **10.13 / 10.12** | 13.97 | 14.26 / 14.49 | **6.93 / 6.92** | ~6.8 (P1) |
| wall s | 48.1 | 54.3 / 51.8 | 45.2 / 44.0 | 49.5 | 59.8 / 57.3 | 51.3 / 50.3 | |

Against the design's expectations (`MallocLargeCache=0`, live figures):

- **Resample: 5.3-5.6 GB sampled, against ~5.8 GB derived.** Measured
  against the base, which reached 14.3-14.5 GB. The rise during resample
  is 3.0-3.3 GB over the mosaic, against 11.9-12.2 GB for the base. These
  peaks are sampled, so a transient shorter than 0.25 s may be missing
  (the lifetime peak was already set in P1).
- **End of P2: 4.48 GB against 4.4 GB derived** (mosaic 2.24 + canvas
  2.15). Before, it was 6.63 GB: the copy is gone (fix 2).
- **The overall peak is 6.92-6.93 GB, not at the load (6.81) but in P4,
  0.1 GB above it.** That is the store filling (4.53 GB) while the source
  mosaic (2.24 GB) is still held by the check-point iterator, plus the
  iterator's transients. The design's "the peak moves to P1's load,
  ~6.8 GB" holds to within 0.12 GB, though the phase that sets the peak is
  P4.
- **The tile is dropped before phase 2 (fix 3).** P4 starts at 2.39 GB on
  the branch and at 4.55 GB on the base: −2.16 GB, the 2.15 GB canvas.
- **The store in the 761 run: 4.53 GB rise in P4 = 16.8 B per point**
  (270.4 M points), against ~4.5 GB / 16.5-17 B expected. As shipped, it is
  4.62 GB = 17.1 B. With `MallocLargeCache=0`, the footprint now counts the
  store, because the slabs are large allocations. So the 2026-10-02 caveat
  "under-counted with the variable" does not apply to the branch. It still
  applies to the base, whose P4 and P5 figures above are not store sizes.

As shipped (no variable), the branch's peak is 10.12-10.13 GB in phase 2,
against 15.36-15.91 GB in resample for the base. The freed P1 copies stay
charged, as before: P1 ends at 7.41 GB, and P4 starts at 5.13 GB where
`MallocLargeCache=0` gives 2.39. The 64 MiB slabs then appear not to reuse
that charged memory, since P4 rises 4.62 GB on top of it. That
explanation is "thought to be"; `vmmap` was not run during the 761 runs.

**Surprising: the base re-measured lower than the 2026-10-02 record as
shipped.** It peaked at 15.36 and 15.91 GB, where the record says 17.83,
with the same `_core` sha256. With the variable it matched the record
(14.26 and 14.49 against 13.97). Two as-shipped base runs differ by
0.55 GB from each other, so the as-shipped peak varies from run to run
with what malloc keeps. The cause was not measured. The comparison above
is the back-to-back one, not the 2026-10-02 record.

## 3. The check-point store (`probe_store.py`)

`probe_store.py` is unchanged: 100 M random points in a 10,000 × 10,000
store, then frozen. It ran twice per tree, as shipped and with
`MallocLargeCache=0` (`runs/probe_store_*.out`).

That probe generates its inputs inside the measured interval, so as
shipped the freed input batches stay in the footprint. A second probe,
`probe_store_isolated.py`, uses the same seed and points but generates
every batch before the baseline reading, so the rise is the store's alone
(`runs/probe_store_isolated_*.out`). Results:

| | base, as shipped | 15e, as shipped | base, `MLC=0` | 15e, `MLC=0` |
|---|---:|---:|---:|---:|
| `probe_store.py`: process footprint after freeze | 3.46 / 3.89 GB (34.6 / 38.9 B per point) | 2.03 / 2.04 GB (20.3 / 20.4 B) | 0.28 / 0.28 GB (under-counted) | 1.66 / 1.66 GB (16.6 B) |
| isolated: store after add | 32.7 / 32.8 B per point | 16.43 / 16.43 B | 9.2 / 9.2 B (under-counted) | 16.42 / 16.42 B |
| isolated: store after freeze | **34.0 / 36.3 B per point** | **16.44 / 16.43 B** | 2.6 / 2.6 B (under-counted) | **16.42 / 16.42 B** |

The design expected **~16.5-17 B per point; measured 16.4**, against
34-36 B for the base in the same probe. That is 30 B in the 2026-10-02
single run; the base varies between runs as shipped. Freezing no longer
changes the footprint. On the base it raises it by 1.3-3.5 B per point.
Live data is 16 B per point. The rest is part-filled last chunks (≤
10,000 × 16 KiB) and slab rounding (25 slabs of 64 MiB, `vmmap`'s
"Malloc Large", 1.6 G). The plain probe's 20.3 B as shipped includes the
probe's freed input batches. The isolated probe shows the store's own
16.4 B.

## 4. Time: resample and the check-point store

Seconds, rounds 1 / 2. `resample` is the driver's span around
`dem_input.resample`. "check points: store" is the CLI's own `--stats`
row: the time inside `store.add` plus `freeze()` (`final_check.run`). The
freeze alone is the driver's span around its `PhaseClock` block, and the
adds are the row less the freeze.

| measure | base as shipped | 15e as shipped | Δ | base `MLC=0` | 15e `MLC=0` | Δ |
|---|---:|---:|---:|---:|---:|---:|
| resample | 20.04 / 19.94 | 13.49 / 12.98 | **−33 %** | 26.29 / 24.67 | 18.14 / 18.39 | **−28 %** |
| check points: store (row) | 5.07 / 5.23 | 5.46 / 5.42 | +5.6 % | 5.92 / 5.98 | 5.49 / 5.43 | −8.3 % |
| … freeze | 3.47 / 3.57 | 3.86 / 3.83 | **+9.2 %** | 3.71 / 3.73 | 3.75 / 3.77 | +0.9 % |
| … adds | 1.60 / 1.66 | 1.60 / 1.59 | −2 % | 2.21 / 2.25 | 1.74 / 1.66 | −24 % |
| check points: project | 8.80 / 8.78 | 8.89 / 8.61 | | 9.78 / 9.75 | 10.11 / 9.67 | |
| `--stats` total | 52.99 / 50.58 | 44.25 / 42.76 | **−16 %** | 58.41 / 56.03 | 50.21 / 49.31 | **−13 %** |

The 2026-10-02 records give resample 17.66 s (as shipped) and 17.40 s (`MLC=0`),
and the store row 5.12 s at 10 m as shipped.

- **Resample is faster, not slower.** It needs 524 blocks and as many
  transformers instead of 76 (design), yet it runs 28-33 % faster. This is
  thought to come from the smaller per-block working set (fewer page
  faults and less memory traffic). It was not profiled.
- **The freeze is 0.26-0.39 s slower as shipped (+7 to +11 % per pair,
  +9.2 % on the means)** and unchanged
  with the variable. This is the gather and scatter the design names.
  The adds are as fast or faster, so the whole store row moves by +0.3 s
  as shipped and −0.5 s with the variable. The whole run is 6.7-8.7 s
  faster either way. I record this rather than call it a regression: the
  acceptance rule's time threshold is `bench.py`'s, the 1 m path does not
  reach this code, and the design named this cost.

## Files

- `2026-10-03/15e-base-7810cf8/`, `15e/`, `15e-base-7810cf8-r2/`, `15e-r2/`: `bench.py`'s
  evidence (`run.json`, `raw.tsv`, generated `README.md`).
- `15e-acceptance/bench_*.out`: `bench.py`'s console output.
  `pmset_*.txt`: power before and after each batch.
- `15e-acceptance/run_all.sh`, `run_phases.py`, `phase_driver.py`,
  `analyse.py`, `counts.py`, `basin_out_crs.wkt`: the 761 runs. Output is
  in `runs/761_t10_{base,new}[_nolargecache]_r{1,2}.*`, and the tables are
  in `tables.md`.
- `15e-acceptance/probe_store.py`, `probe_store_isolated.py`: the store probes.
  Output is in `runs/probe_store*_{base,new}*.out`.

## Reproduce

Check out the base: `git worktree add --detach <dir> 7810cf8`. Then build
both trees and run the bench pairs from the branch's venv:

```
python tools/bench.py run --label 15e-base-7810cf8 --tree <dir> --out-root docs/benchmarks
python tools/bench.py run --label 15e --tree . --out-root docs/benchmarks \
  --baseline docs/benchmarks/2026-10-03/15e-base-7810cf8
```

The pair-2 runs are the same with `-r2` labels. Make one venv per tree, with
no packages, and a `.pth` that lists `<tree>/build-bench/pkg` and then the
branch venv's site-packages. With `D=../rasputin_data`, run from
`15e-acceptance/`:

```
./run_all.sh <venv-base> <venv-new> $D <mesh dir>
<venv-new>/bin/python analyse.py runs 761_t10_base_r1 … > tables.md
for v in base new; do <venv-$v>/bin/python probe_store.py; MallocLargeCache=0 <venv-$v>/bin/python probe_store.py; done
(and the same for probe_store_isolated.py)
```
