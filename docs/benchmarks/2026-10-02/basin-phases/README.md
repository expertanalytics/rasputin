# Where the memory goes, phase by phase: one BHO level-3 sub-basin (@perf, 2026-10-02)

The measurement that `docs/research/basin-memory-options.md` (branch
`worktree-basin-memory`, §1) asks for: mesh one level-3 sub-basin of the São
Francisco from the cache, sample memory over time, cut it by phase, and check
the note's per-phase model before it is extrapolated to the whole basin.

Branch `worktree-basin-phases`, cut from `worktree-15c-2` at `0b42f5b`. No
production code changed. Apple M1 Max (8 P + 2 E cores, 32 GiB), macOS 27.0,
Python 3.14.7, **AC power** for every run (`pmset -g batt` before and after, in
each `runs/<tag>.json`). `_core` built Release (`-O3 -DNDEBUG`, AppleClang 21)
by `tools/bench.py`'s `build()` into `build-bench/` of this tree (sha256
`41a22a6168db…`, the same `.so` as `../basin-anadem/`), used from a venv of
its own with `build-bench/pkg` on its path. One run per configuration, no
repeats: these are memory figures, not timings, and the timings below are
given only to place the phases.

## The domain

BHO level 3 is not served as a layer (ANA's `Divisao_de_bacias` stops at
level 2), so `fetch_level3.py` builds each `76k` outline as the union of the
BHO 2017 **50k** drainage areas whose `COBACIA` starts with `76k`
(`SPR/BHO2017_50K_AREADRENAGEM`, EPSG:4674), and plans each one's target grid
in the basin's `--out-crs` (`runs/basin_out_crs.wkt`, the one
`../basin-anadem/` pasted). The nine areas sum to 633.6 k km², the level-2
basin's 633.6 k km². Table: `runs/level3_units.json`; outlines in
`../rasputin_data/sao_francisco_piece/bho2017_50k_level3/` (not in the
repository; rerun the script to regenerate). A measurement domain per B13 (c),
not a product.

**761** was chosen: 209,316 km², 36,983 vertices, target grid 19,377 × 27,786
= **0.538 G nodes** at 30 m, the only unit in the asked 0.3-0.6 G band (the
next largest, 769, is 0.196 G). `rasputin fetch anadem-v1 --domain … --out-crs
…` found all 2,223 blocks present (`runs/fetch_761.out`).

## Method

- `run_phases.py` runs `/usr/bin/time -l python phase_driver.py …` and
  samples the CLI process every 0.25 s with `proc_pid_rusage`: resident size,
  `phys_footprint` (what `time -l` calls peak memory footprint; what the OS
  charges, compressed and swapped pages included) and the lifetime peak
  footprint. Swap in use is watched (kill at +8 GB; never reached).
- `phase_driver.py` is the CLI, in-process and unchanged, with markers written
  at the start and end of `assemble`, `resample` (and the `DemTile`
  constructor inside it, the copy at its end), `refine` (phase 1),
  `final_check.run` (store fill and phase 2), `refine_points` (phase 2) and
  every `PhaseClock` row. Each marker records the footprint and the lifetime
  peak at that instant, so a phase in which the lifetime peak rose has an
  **exact** peak; elsewhere the peak is the largest sample (0.25 s).
- `counts.py` takes the counts the CLI printed (final triangles; phase 2's
  check points and insertions from the mesh header). Phase 1's triangle count
  is not printed; it is taken as final less 2 per phase-2 insertion.
- `analyse.py` prints the tables below from `runs/`.

```
rasputin mesh --dem anadem-v1 --cache ../rasputin_data/cache \
  --domain ../rasputin_data/sao_francisco_piece/bho2017_50k_level3/761_epsg4674.geojson \
  --out-crs "$(cat runs/basin_out_crs.wkt)" --tolerance {50,10,5} --binary --out <scratch> --stats runs/<tag>.stats.md
```

Meshes were written to `../rasputin_scratch/basin-phases/` (not kept in the
repository; the command above regenerates them).

## Finding first: macOS malloc keeps freed memory charged to the process

Two things the options note's model does not have, both measured here:

1. **Freed large blocks stay in the footprint.** `probe_large.py`: allocate
   and free 4 GB of NumPy arrays; the footprint stays at 4.02 GB, and
   `malloc_zone_pressure_relief` does not return it. With the environment
   variable `MallocLargeCache=0` it drops to 0.02 GB (`runs/probe_large.out`).
   `vmmap` names the region "Malloc Large (empty)": 4.7 GB of it after the
   761 source window is loaded and dropped (`runs/probe_tm.out`).
2. **The check-point store fragments the small-block zone.** `probe_store.py`:
   100 M points in a `CheckPoints` (10,000 rows), frozen. Live bytes 1.5 GB
   (`vmmap`'s allocated, 16 B per point as designed), but the zone holds
   3.0-3.9 GB dirty and the footprint is **3.0-4.1 GB, 30-41 B per point**
   (three runs; `runs/probe_store.out`). `freeze()`'s `shrink_to_fit` raised
   the footprint (2.80 → 3.04 GB), it did not lower it.

So each run was made twice: as shipped (the figure that decides whether a run
fits), and with `MallocLargeCache=0` (`*_nolargecache`). With that variable
the large-array phases (P1, P2) give live sizes, but `proc_pid_rusage`'s
footprint then **under-counts** the store: the probe shows 0.27 GB footprint
for 1.5 GB live (`vmmap` says 1.6 GB, `runs/probe_store_nolargecache.out`).
Its P4-P6 figures are therefore not used for the store. The mesh's rises
(P3, P5) are given from both kinds of run, and both are **lower bounds**: as
shipped, a phase can reuse freed memory that is already charged; with the
variable, the footprint may under-count as it does for the store.

## Per phase, sub-basin 761 (0.538 G target nodes, 0.561 G source nodes, 270.4 M check points)

Peak footprint per phase, GB. "s" marks a sampled peak (0.25 s), the rest are
exact. Full tables: `python analyse.py runs <tag>…`; raw data in `runs/`.

| phase | 50 m as shipped | 10 m as shipped | 10 m `MallocLargeCache=0` | 5 m `MallocLargeCache=0` |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 7.41 | 7.41 | 6.81 | 6.81 |
| P2 resample (blocks) | **18.27** | **17.83** | **13.97** | **14.10** |
| P2 end (`DemTile` copy) | 10.42 s | 10.19 s | 6.63 s | 6.63 s |
| P3 refine phase 1 | 10.32 s | 10.13 s | 4.70 s | 5.03 s |
| P4 store fill + freeze | 11.07 s | 11.51 s | (under-counted) | (under-counted) |
| P5 phase 2 | 11.10 s | 12.16 s | (under-counted) | (under-counted) |
| P6 trim + write | 6.89 s | 7.84 s | | |
| `time -l` peak footprint | 18.27 | 17.83 | 13.97 | 14.10 |
| wall s | 39.7 | 48.1 | 49.5 | 66.7 |
| final triangles | 335,681 | 2,940,024 | 2,940,024 | 7,791,814 |

The peak is set in P2 at every tolerance; tolerance moves it by < 0.5 GB.

## Per unit, against the options note

| item | note's estimate | measured on 761 | how |
|---|---|---|---|
| source mosaic, held | 4 B per source node | **4.0 B** (2.24 GB, float32 19,685 × 28,496) | array size |
| source mosaic, while loading | not in the model | **12.0 B** per source node (3 copies: `decode_window`'s buffer, `DemTile`'s copy, one more), 6.76 GB | P1 rise, `MallocLargeCache=0` |
| … charged after loading, as shipped | not in the model | **13.1 B** per source node (7.35 GB): the two freed copies stay | P1 end footprint |
| target canvas | 4 B per node, + 4 B copy at P2's end | **4 B + 4 B** (2.15 GB each; mosaic + 2 canvases = 6.63 GB) | P2 end |
| resample transients | 149 B per node × 256 rows × cols × threads | **133-135 B** (11.6-11.8 GB rise less the canvas; 10 threads, 27,786 cols) | P2 peak, `MallocLargeCache=0` |
| … as shipped | | +4.3 GB on top (the freed mosaic copies) | 18.27 vs 13.99 at 50 m |
| check-point store | 16 B per point frozen, up to 32 B while growing | **16 B live; 30-41 B charged** after freeze | `probe_store.py` |
| … in the 761 run | | P5 peak less mosaic and tile: 7.8 GB = **29 B** per point | 10 m as shipped |
| phase 1 mesh | 310 B per triangle | **≥ 73-79 B** per phase-1 triangle (refine's own rise, transients included; lower bound) | P3, `MallocLargeCache=0` |
| phase 2, two meshes | 2 × 310 = 620 B per triangle | **≥ 170-220 B** per final triangle (5 m `MallocLargeCache=0`: 1.35 GB / 7.8 M; 10 m as shipped: 0.65 GB / 2.9 M; lower bounds) | P5 rise |
| check points per box node `f` | 0.35-0.44 (basin) | 0.50 for 761 (fill 0.43); **1,292 points per km² of domain** | `store_size` |

Where the note was wrong:

- **P1 is not 4 B per source node but 12 B at the load, and 13 B stays
  charged** (the freed copies are kept by malloc). Not in the model at all.
- **Freed memory does not leave**: what the note costs as transient (P1's
  copies, P2's copy, the store's growth) stays in the footprint as shipped.
- **The store is charged 2-2.5× its design**, from small-zone fragmentation,
  not from vector growth alone; freezing does not give it back.
- **The mesh rose far less than costed**: ≥ ~75 B per phase-1 triangle and
  ≥ ~170-220 B per final triangle in phase 2, against 310 and 620. These are
  lower bounds (see above), so "the mesh is cheaper" is not settled; but on
  761 the mesh rise never exceeded 1.4 GB, and at basin scale even 620 B per
  triangle is not what decides.
- The resample transient was close (133-135 B against 149, −10 %), and the
  canvas and its copy exact.
- The model's P4 (mosaic + canvas + store) holds as shipped only with the
  store at ~29-41 B, not 16-32 B.

## The whole basin, extrapolated with these figures

Basin, from `../basin-anadem/`: 2.13 G source nodes, target grid 50,315 ×
40,943 = 2.06 G nodes, 633.6 k km² (so ~0.82 G check points at 1,292 per
km²). Triangles scaled by area from 761 (×3.03): ~8.9 M at 10 m, ~23.6 M at
5 m, costed at 200 B (measured, a lower bound) to 620 B (the note) each. *Est.* throughout: one sub-basin's per-unit figures, scaled. Whether
malloc's large cache keeps 8.5 GB blocks as it kept 2.2 GB ones on 761 was not
measured; "as shipped" assumes it does, "live" is what the code holds.

| phase | live (GB) | as shipped (GB) |
|---|---:|---:|
| P1 end: mosaic | 8.5 (25.6 while loading) | 27.9 |
| P2 peak: mosaic + canvas + 135 B × 256 × 40,943 × 10 | **31.0** | ~48 |
| P2 end: mosaic + 2 canvases | 25.0 | ~42 |
| P5, 10 m: mosaic + canvas + store + mesh (1.8-5.5) | 16.7 + 13.1 + 1.8-5.5 = **31.6-35.3** | 16.7 + 24-34 + 1.8-5.5 = **42-56** |
| P5, 5 m: same, mesh 4.7-14.6 | **34.5-44.4** | **45-65** |

So on 32 GiB the basin does not fit at any tolerance as shipped, and P1's load
alone (25.6 GB live while loading) leaves no room for P2. That agrees with the
basin run dying in P2, and with the note's conclusion, but for more reasons
than the note gives: the note's P1, the kept-freed memory and the store's
charge are each larger than costed; the mesh is smaller.

What this means for the note's option (a) is @architect's to weigh; three
measured facts bear on it. Releasing freed memory needs either fewer copies or
an allocator that gives it back: `MallocLargeCache=0` cut 761's peak from 17.8
to 14.0 GB with no code change. The store needs its rows reserved from a
counting pass or allocated as one block, or it is charged ~2.5× its 16 B.
With the mosaic, the copies and P2's block transients gone, P5 at 5 m is
≈ canvas 8.2 + store 13.1 + mesh 4.7-14.6 ≈ 26-36 GB live, *est.*

## Files

- `fetch_level3.py`: level-3 outlines and their grid sizes → `runs/level3_units.json`.
- `run_phases.py`, `phase_driver.py`: the sampled, marked run → `runs/<tag>.{json,mem.csv,markers.jsonl,log,stats.md}`.
- `counts.py`, `analyse.py`: counts and tables → `runs/<tag>.{counts,phases}.json`.
- `probe_alive.py` → `runs/761_probe_alive.{json,out}`; `probe_tm.py` → `runs/probe_tm.out`;
  `probe_large.py` → `runs/probe_large.out`; `probe_store.py` → `runs/probe_store*.out`.
- The 50 m as-shipped run was made before the driver recorded the mosaic's
  shape; its `assemble` end marker carries the shape from `761_t10` (same
  plan), marked `shape_added_from`.
