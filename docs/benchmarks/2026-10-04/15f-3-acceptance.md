# Increment 15f-3 acceptance: the edge strip (@perf, 2026-10-04)

Branch `worktree-15f-3` at `83c7fd2` (`@reviewer` round 2 APPROVED), against
its merge base `c193cb1`, back to back. AC power for every run (100 %,
charged; `pmset -g batt` before and after each). The rule is
`docs/increments/15f-edge-strip.md`, "Acceptance (`@perf`)", 15f-3: (1)
`bench.py`, (2) Bygdin with CORINE and the independent strip check, (3) a basin
piece. Swap was 10.54 to 10.57 GB used of 11.26 GB throughout. That is
history: it did not grow during any run, and `memory_pressure` read 51 % free
at the start.

## Verdict

**ACCEPTED, with one cost finding.**

- **Correctness and geometry are as designed.**
  - The tile's mesh is byte-identical to the base's. The quarter's differs
    by one vertex (the strip's), as the design expects.
  - The strip checks find 0 points over tolerance on 15f-3 everywhere, where
    the base has thousands. Refine time is within noise.
- **The cost finding: the strip costs more end to end than its phases show.**
  - On the 1 m benchmark, the whole process takes +36 to +56 % longer, and
    the time inside the CLI (`app_s`) +52 to +111 %.
  - This holds on the tile too, where the mesh is unchanged.
  - cProfile puts it in `_core.refine_strip`: 0.341 s on the tile, against
    refine's 0.198 s. The strip's timed phases cover 0.003 s of that. The
    rest lands in `--stats`' `other` row.
  - Bygdin at 1 m: `refine_strip` takes 1.51 s, of which 0.034 s is in the
    timed phases. The process takes +1.5 s (+34 %) and peak RSS rises by
    +0.29 GiB.
  - Where inside `refine_strip` the time goes was not measured: the Release
    `_core` has no symbols for `sample`. It is thought to be the set-up
    (rebuilding the mesh from arrays) or the DEM rescan (F2).
  - The reprojected path does not call `refine_strip`. There the cost is
    +2.6 / +6.1 / +7.5 % of process time on Velhas at 20 / 10 / 5 m.
  - The design gives no threshold for end-to-end time, so this goes to
    `@architect` rather than deciding the verdict.

## (1) `tools/bench.py`, 1 m benchmark and thread sweep

Two balanced batches, B N N B then N B B N (`15f-3-acceptance/pairs.sh`,
`pairs.log`), 5 repeats, threads 0 and 1 to 20, both domains, Release `_core`
from `bench.build()`: base `0b3362e5…` (the same binary as `24-on-r2`) and
15f-3 `0ee8cd40…`. The tables come from `summarize_bench.py`
(`tables-bench.md`). Pooled medians over 4 runs a side:

| domain | measure | 1 thread | 20 threads | default | 2-20 range |
|---|---|---:|---:|---:|---:|
| tile | refine_s | +0.4 % | +1.0 % | +0.6 % | |
| quarter | refine_s | -1.6 % | -0.2 % | -0.9 % | |
| both | refine_s | | | | -2.7..+1.6 % |
| tile | app_s | +51.9 % | +99.6 % | +100.9 % | |
| quarter | app_s | +53.8 % | +110.4 % | +111.0 % | |
| tile | proc_s | +36.6 % | +54.8 % | +55.8 % | |
| quarter | proc_s | +36.2 % | +55.1 % | +54.9 % | |

- **Geometry.**
  - Tile: identical, `11741a81…` in all eight runs, as E7 and the
    grid-line remark predict.
  - Quarter: base `ccebf96a…`, 15f-3 `a60fb597…`, with one more vertex
    (Delaunay checked 641,791 against 641,792). Refine's counters are
    identical (rounds 41, inserted 213,464).
  - Quality is the same on both domains: worst angle 0.6296° / 0.3955°, max
    degree 74 / 18, 0 Delaunay violations, within tolerance.
- **The strip's phases** in a 15f-3 `--stats` run
  (`profile/bench_tile_15f3.stats.md`): generate 0.001 s, scan 0.001 s,
  split 0.001 s, against `other` 0.341 s and refine 0.198 s.
  `profile/bench_{tile,quarter}_15f3.txt` hold the cProfile lines for
  `refine_strip`.

## (2) Bygdin with CORINE, 10 m and 1 m

16c's inputs and options, through each tree's Release package
(`bygdin.sh`). Order B N N B, 2 timed `--binary` runs per block and
tolerance (4 per side), then one `--ascii` run per side for `quality()` and
the check. Logs, `--stats` and `/usr/bin/time -l` output are in `bygdin/`.
The tables, from `summarize_bygdin.py`, are in `tables-bygdin.md`.

| | 10 m base | 10 m 15f-3 | 1 m base | 1 m 15f-3 |
|---|---:|---:|---:|---:|
| process wall, median of 4 (s) | 2.36 | 2.32 | 4.32 | 5.81 |
| peak RSS (GiB) | 0.35 | 0.35 | 0.82 | 1.11 |
| strip points (`line_points_checked`) | - | 175,246 | - | 180,436 |
| `line_points_on_nodata` | - | 0 | - | 0 |
| strip inserted / DEM nodes the rescan inserted | - | 199 / 40 | - | 9,012 / 2,001 |
| `line_points_refused` (largest error) | - | 0 (0) | - | 0 (0) |
| `line_max_error_m` | - | 9.98676 | - | 0.999843 |
| refine rounds / inserted | 17 / 17,879 | 17 / 17,879 | 24 / 549,527 | 24 / 549,527 |
| strip time: generate + scan + split (s) | - | 0.008 | - | 0.041 |
| `other` (s, median) | 0.005 | 0.048 | 0.005 | 1.420 |
| worst angle | 0.00408° | 0.00875° | 0.00408° | 0.02177° |
| share under 1° | 0.070 % | 0.053 % | 0.100 % | 0.009 % |
| max degree | 16 | 14 | 29 | 18 |
| Delaunay violations (bench.py) | 0 | 0 | 0 | 3, all exact ties (below) |

The design asks for the share under 10°, and `quality()` does not compute
it. From each ASCII run's `--stats` quality table: 1.70 % (base) and 1.68 %
(15f-3) at 10 m, and 1.36 % and 0.96 % at 1 m.

**The independent check** is `strip_check.py`. It does not import
`tin_engine`.

- **Its inputs.** It reads the domain with shapely, CORINE from the
  GeoPackage with sqlite3, the DTM10 tiles with tifffile and their `.tfw`
  files (nodes at pixel centres), and the ASCII mesh with its own parser.
- **What it checks.** The crossings of the input outline and the CORINE
  borders (clipped, noded) with the DEM's grid lines, and the midpoints
  between neighbouring crossings. Each gets z bilinear from the DEM and is
  compared with the mesh's z along the nearest constraint line, within 1 cm.
  The input lines are at most 0.69 mm from the mesh's lines (1 mm snap).
- **The alignment control.** 33,334 to 561,792 mesh vertices lie on DEM
  nodes, and each carries the node's value exactly (largest difference 0.0).
- **The plant control.** With the DEM shifted 10 m east, the check finds
  44,186 crossings over at 1 m.

| mesh | crossings over (of 83,184) | largest error | midpoints over (of 91,555) | largest |
|---|---:|---:|---:|---:|
| base 10 m | 575 | 62.61 m | 521 | 61.05 m |
| 15f-3 10 m | 0 | 9.98687 m | 0 | 9.86647 m |
| base 1 m | 18,567 | 72.93 m | 18,066 | 70.11 m |
| 15f-3 1 m | 1 (by 1.3e-5 m) | 1.000013 m | 1 (by 1.2e-4 m) | 1.000117 m |

- **The base's figures are the projected path's Surprise 3, measured for
  the first time:** 575 crossings over at 10 m and 18,567 at 1 m, up to
  72.9 m.
- **On 15f-3**, the check run on the mesh's own constraint lines (the same
  script with `--lines-from-mesh`) finds **0 over at both tolerances**. Its
  largest errors, 9.98676 m and 0.999843 m, equal `line_max_error_m` to the
  printed digits.
- **The one point over on the input lines** lies 0.29 mm from the mesh's
  line, where the surface slope makes 1.3e-5 m. That is the 1 mm snap
  between input and mesh, not a missed point. `line_points_refused` is 0.
  So by the rule's terms ("0 over, or exactly the refused points") the
  result is 0 over.
- The check's output is in `bygdin-check/*.json`.

**The 3 Delaunay "violations" at 1 m** (15f-3 only) are exactly cocircular
quadrilaterals. Their corners sit on 2.5 m multiples (midpoints and
crossings). Exact incircle on the positions rounded to the millimetre gives 0
for all three. The written world coordinates are off by about 2e-10 m
(`x0 + col dx` in double), which tips `bench.py`'s float test to "inside".
They are ties in the lattice where the mesher decides, not violations. This
is measured on three cases only.

**Q1: how many midpoints are inserted at all.** `q1_vertices.py`
(`q1_vertices.txt`) sorts each mesh's constraint-line vertices into those on
a grid line and those off every grid line, and compares 15f-3 with the base:

| run | strip inserted | at crossings | at midpoints |
|---|---:|---:|---:|
| Bygdin 10 m | 199 | 189 | 10 |
| Bygdin 1 m | 9,012 | 8,391 | 621 |
| Velhas 20 m | 6 | 5 | 1 |
| Velhas 10 m | 45 | 43 | 2 |
| Velhas 5 m | 275 | 253 | 22 |

The two columns sum to `line_points_inserted` in every row. At Bygdin they
also agree with `strip_check.py`'s matching of vertices to its own candidate
points (8,391 and 622 at 1 m). **Midpoints are 4 to 17 % of the strip's
insertions: 6.9 % at Bygdin 1 m and 8.0 % at Velhas 5 m, the finest
tolerances run.** The 17 % is 1 of 6.

## (3) The Velhas piece (BHO 76949) on ANADEM, EPSG:31983

`velhas.py run` (15c-2's `run_geo.py`, rewritten because its VTK-header
regex no longer matches either tree). Order B N N B, 2 timed runs per
tolerance per run (4 per side), then one `--ascii` run each. ANADEM came from
the cache with no network fetch. The tables, from `summarize_velhas.py`, are
in `tables-velhas.md`.

| | 20 m base | 20 m 15f-3 | 10 m base | 10 m 15f-3 | 5 m base | 5 m 15f-3 |
|---|---:|---:|---:|---:|---:|---:|
| process wall, median (s) | 3.04 | 3.12 | 4.07 | 4.32 | 6.79 | 7.30 |
| peak RSS (GiB) | 1.91 | 1.89 | 1.91 | 1.86 | 1.88 | 1.98 |
| strip points / inserted | - | 73,449 / 6 | - | 73,511 / 45 | - | 73,657 / 275 |
| `line_points_refused`, `line_points_on_nodata` | - | 0, 0 | - | 0, 0 | - | 0, 0 |
| `line_max_error_m` | - | 19.7943 | - | 9.9575 | - | 4.9985 |
| `max_error_m` | 19.99997 | 19.99997 | 9.999996 | 9.999996 | 4.9999997 | 4.9999997 |
| final check inserted (rounds) | 5,686 (8) | 5,691 (8) | 29,368 (11) | 29,406 (11) | 142,643 (12) | 142,874 (12) |
| strip generate + final check scan + split (s) | 0.218 | 0.237 | 0.411 | 0.442 | 0.909 | 0.964 |
| 15c's source-node check: interior / strip over | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| strip check: crossings over | 15 | 0 | 118 | 0 | 633 | 0 |
| strip check: largest error (m) | 23.83 | 19.79 | 23.83 | 9.958 | 23.83 | 4.998 |
| worst angle | 0.1135° | 0.1135° | 0.0418° | 0.0525° | 0.0303° | 0.0558° |
| Delaunay violations | 0 | 0 | 0 | 0 | 0 | 0 |

- **15c's source-node check stays at 0 over.** It covers 13,791,972
  interior nodes and 27,806 strip nodes (`velhas.py check`, which runs
  `run_geo.independent_check`). Its control, the mesh moved 30 m east,
  finds 111,789 / 845,911 / 3,013,718 over.
- **The independent strip check on the resampled grid finds 0 over on
  15f-3.**
  - How it works: `strip_check.py --resampled` rebuilds the target grid
    independently: node (R, K) at (30 K, -30 R), with z bilinear from the
    ANADEM window at the node moved by pyproj, rounded to float32. It then
    checks the mesh's constraint lines against that grid.
  - Alignment: 784,578 mesh vertices on target nodes agree with the rebuilt
    values to within 6.1e-5 m (one float32 ulp at these heights).
  - Its control, the grid's values taken from 30 m further west, finds
    4,374 crossings over at 5 m.
  - On the base it measures the reprojected path's Surprise 3: 15 / 118 /
    633 crossings over, up to 23.8 m.
  - Its largest errors on 15f-3 equal `line_max_error_m`.
- **The strip's time and memory beside the final check's.** At 5 m: strip
  generate 0.003 s; final check scan 0.786 s (base 0.742) and split 0.175 s
  (base 0.167). The final check now carries the strip in its own loop. Peak
  RSS is within the run-to-run spread (1.86 to 1.98 GiB).
- **Not run: BHO level-3 unit 766 at 20 and 5 m.** This run's brief named
  the Velhas piece only. With swap nearly full (0.7 GB free) a 2.8 GB-peak
  run was not added unasked.

## Where refine_strip's untimed time goes

Measured after the verdict above, still alone and on AC.

**It is a per-call fixed cost that scales with the size of refine's output
mesh, not with the strip's points.**
- 90 % of it is spent rebuilding the lattice mesh from refine's output
  arrays (`detail::to_lattice`, and in it `LatticeMesh::build`).
- Another 5 to 6 % is `legalise_all` over the whole rebuilt mesh.
- The binding's copies and the GIL take about 0.1 %.

**Method.**
- **Build.** `_core` from `588e879` (15f-3's code at `83c7fd2`), built as
  RelWithDebInfo with the Release level: `-O3 -DNDEBUG -g`, hardening on, in
  `build-prof/`. pybind11 keeps the symbols in that build type, and its
  package is assembled as `bench.build()` does.
- **Driver.** `15f-3-acceptance/prof_strip.py` runs `rasputin mesh` from
  that package and wraps `edge_strip.run`. After the real call it calls
  `_core.refine_strip` K more times on the same inputs, then K times with an
  **empty strip** (0 points, from `constraint_check_points` on no edges). It
  prints each call's wall time and its timed phases (`profile/refine_strip_calls_*.json`).
- **Sampling.** macOS `sample` at 1 ms over the repeated calls. The call
  tree under `terrain::refinement::refine_strip` is summed by
  `sample_tree.py` (`profile/refine_strip_sample_*.txt`).
- **Inputs.** The 1 m benchmark tile (K = 15, 472,374 triangles,
  23,664 strip points) and Bygdin 1 m with CORINE (K = 4,
  1,145,434 triangles, 180,436 strip points).

**Per call, with and without the strip's points** (medians while sampled;
an unsampled trial on the tile gave 0.30 to 0.33 s for both):

| case | strip points | wall per call | timed phases (scan + split) | same call, empty strip |
|---|---:|---:|---:|---:|
| tile 1 m | 23,664 | 0.391 s | 0.003 s | 0.384 s (timed 0.002 s) |
| Bygdin 1 m | 180,436 | 1.665 s | 0.041 s | 1.624 s (timed 0.005 s) |

Removing every strip point saves only the timed phases: 0.007 s on the tile
and 0.041 s on Bygdin. So the untimed time does not depend on the strip. It
grows with the mesh refine_strip receives: 0.81 µs per triangle on the tile,
1.42 µs on Bygdin. Why the cost per triangle rises with size was not measured.

**Where it goes**, as a share of the samples under `refine_strip`
(tile 5,838 samples, Bygdin 8,284):

| function (inclusive) | tile | Bygdin |
|---|---:|---:|
| `detail::point_loop` (all of refine_strip) | 100 % | 100 % |
| `detail::to_lattice`: rebuild refine's output as a `LatticeMesh` | 90.4 % | 90.9 % |
| `LatticeMesh::build` | 70.7 % | 74.7 % |
| `std::unordered_map<uint64, uint32>::emplace` in `build` (the directed-edge map) | 39.6 % | 49.8 % |
| `free` (`_xzm_free_tc`, `_free`) in `build`: the map's nodes released | 15.7 % | 14.5 % |
| `to_lattice`'s own code: the `lattice_position` loop and the `std::map` constraint lookups, inlined | about 18 % | about 15 % |
| `mesh::legalise_all` over the whole rebuilt mesh (`must_flip`, `lattice_incircle`) | 6.3 % | 4.7 % |
| sorting (`std::__introsort`) | 1.4 % | under 1 % |
| the binding (`start_mesh`'s copies, result conversion) | 0.1 % | 0.1 % |

- **The map.** The `unordered_map` in `LatticeMesh::build` gets no
  `reserve` (`include/terrain/mesh/lattice_mesh.hpp`, `build`). It holds 3
  entries per triangle, each a separate allocation, all freed when `build`
  returns.
- **The constraint map.** In `to_lattice`
  (`include/terrain/refinement/refine.hpp`), constraint edges are looked up
  in a `std::map`, three times per triangle.
- **Why refine does not pay this.** `refine` runs the same `to_lattice` on
  the start mesh, which has 15,614 triangles on Bygdin. `refine_strip` runs
  it on refine's output, which has 1,145,434.
- **The GIL is not a factor.** The binding releases it before the call, and
  nothing else was running.
- **About 553 / 1,805 samples sit at one unsymbolised address** of `_core`
  (`+0x7218`, a stub region). They are counted inside the `build` subtree
  above. `atos` does not resolve the address.

Nothing was changed. This is a measurement for `@architect`. The design
question is whether refine_strip should take refine's `LatticeMesh` instead
of rebuilding it from arrays.

## Files and clean-up

- Run directories in `bench.py`'s format: `15f-3-base-r{1..4}/` and
  `15f-3-r{1..4}/`.
- `15f-3-acceptance/`:
  - the drivers: `pairs.sh`, `bygdin.sh` and `velhas.sh` with their logs,
    and `velhas.py`;
  - the checks: `strip_check.py` and `q1_vertices.py`;
  - the summarizers and their `tables-*.md`;
  - `bygdin/`: logs, `--stats` and `time -l` output;
  - `velhas/<label>/`: `results.json`, logs and `--stats`;
  - `bygdin-check/` and `velhas-check/`: the check's JSON;
  - `profile/`: the cProfile extracts and two `--stats` files.
- The meshes (Bygdin ASCII up to 71 MB each, Velhas ASCII up to
  99 MB), the base worktree and its `build-bench/` were in this session's
  scratchpad and have been deleted. To regenerate them:
  1. `git worktree add --detach B c193cb1`.
  2. Build each tree with `bench.build(bench.make_runner(), Path(T), "on")`.
  3. Run `pairs.sh W B S`, `bygdin.sh W B S` and `velhas.sh W B S`.
  4. Then run `strip_check.py` as in the log, `velhas.py check
     velhas/<label> --shift 30`, and the summarizers.
- The sanitize-first rule does not apply: no C++ was patched for this run.
  The profiling build (`build-prof/`, RelWithDebInfo at `-O3 -g`) is the
  branch's code unchanged. It is gitignored and stays in the worktree.

## Log
- 2026-10-04T00:12:15Z: start. AC 100 %, charged. Swap 10.57 GB used of 11.26 GB (history; memory_pressure 51 % free). Load ~16 (dasd).
- 2026-10-04T00:36:16Z: (1) bench done, AC 100 % all 8 runs. Refine within noise (tile -0.0..+1.0 %, quarter -1.6..-0.2 % at the listed cells). Tile mesh identical; quarter mesh differs (+1 vertex), as the design expects. app_s +52..+111 %, proc_s +36..+56 % (tile too, where nothing changes): the strip's cost end to end. Next: Bygdin.
- 2026-10-04T00:41:16Z: (2) Bygdin done (AC 100 % every run, swap flat). Independent check on base: 575 crossings over at 10 m (max 62.6 m), 18,567 at 1 m (max 72.9 m). On 15f-3: 0 over on the mesh's own lines at both tolerances (max 9.98676 / 0.999843 m, equal to rasputin's line_max_error_m). On the input lines, 1 crossing over by 1.3e-5 m at 1 m, at a 0.29 mm offset between input line and mesh line. 15f-3 at 1 m: 3 bench.py Delaunay "violations", all exactly cocircular (exact incircle 0 at mm-rounded coordinates; written coordinates off by ~2e-10 m). Cost: refine_strip 1.51 s at 1 m, of which 0.034 s in its timed phases (cProfile); on the bench tile 0.341 s against refine's 0.198 s. Next: Velhas.
- 2026-10-04T00:46:14Z: (3) Velhas runs done (AC, swap flat at 10.54 GB). Strip check on the rebuilt resampled grid: base 15 / 118 / 633 crossings over at 20 / 10 / 5 m (max 23.8 m); 15f-3 0 over at all three; control (+30 m grid shift) 4,374 over. Running 15c's source-node check now.
- 2026-10-04T00:52:01Z: write-up done; verdict ACCEPTED with the refine_strip cost finding.
- 2026-10-04T00:56:57Z: refine_strip profiled (RelWithDebInfo -O3 -g, sample): 90 % to_lattice rebuild of refine's output, fixed per call, independent of strip points.
