# Fitting the São Francisco basin in memory: where it goes, and the options

Status: research note, 2026-10-02, `@architect`, on `worktree-15c-2` at
`c47bff6`. Analysis only: the choice between the options is Ola's. Figures
marked *est.* are worked out from the code, not measured.

## 1. Where the memory goes

The run being explained is @perf's 50 m attempt on BHO level-2 basin 76
(`docs/benchmarks/2026-10-02/basin-anadem/README.md`). Its peak footprint
was 33.2 GB, and it was stopped at 189 s without finishing. The target grid
is 50,315 rows × 40,943 columns (2.06 G nodes); the source box is 2.13 G
ANADEM nodes. The code path is `dem_input._open_reprojected`, then
`cli._dem_mesh`, then `final_check.run`.

| phase | what is alive (on top of the previous phase, unless freed) | size, basin |
|---|---|---:|
| P1 assemble | source mosaic, float32 canvas of the whole source box (`mosaic.assemble`) | 8.5 GB |
| P2 resample | mosaic + target canvas + **per-block temporaries** | 8.5 + 8.2 + up to ~15.6 GB |
| P2 end | `resample` returns `DemTile(meta=meta, array=canvas)` (`target_grid.py:205`). The public constructor copies (`io/models.py:108`), so for a moment there are two canvases | 8.5 + 8.2 + 8.2 GB |
| P3 refine | mosaic (the lazy `check_point_blocks` iterator holds `TileWindows(mosaic.tile)`, `dem_input.py:208,215`) + target tile + phase-1 mesh | 16.7 GB + mesh |
| P4 store | everything from P3 (D1 says "drop the target tile", but `opened.tile` stays referenced through `final_check.run`, `cli.py:1495`) + `CheckPoints` filling up | 16.7 + mesh + 12-29 GB *est.* |
| P5 phase 2 | store, frozen and shrunk to 16 B per point (`check_points.hpp:88`) + two mesh copies | ≥ 12-14 GB + mesh |

- **The per-block temporaries were measured here** (tracemalloc, one
  256-row block of the basin grid, one thread; the probe is
  `/tmp/blockprobe.py`, not kept). `resample` peaked at **149 B per node**,
  1.56 GB per block. `_pool` runs `os.cpu_count()` blocks at once, which is
  10 on Ola's M1 Max, so up to ~15.6 GB of transients *est.*. Together with
  the mosaic and the canvas that comes to ~32 GB. That matches the 33.2 GB
  measured at 189 s, so **the run most likely died in P2, before refine
  started**. This is an inference; which phase it was was not measured.
- Even with P2 fixed, P4 does not fit: 16.7 GB + ~0.73-0.9 G points × 16 B,
  with up to 2× more while the per-row vectors grow (D4 says so itself).
- Mesh, from the accuracy note's triangle counts (with the final check) and
  the 23 file's 310 B per triangle (a max-RSS fit, refine transients
  included): ~1 GB at 50 m *est.*, ~10 GB at 5 m (33 M triangles), ~30-34 GB
  at 2 m (~100-110 M triangles, extrapolated).

**@perf measurement to confirm** (one run, ~1 h, no code change): mesh one
level-3 sub-basin of 76, chosen so that `target_grid_for` gives a
200-400 M node grid. Run at 50 m and at 5 m, with a sampler that records
`phys_footprint` every 0.25 s, aligned to the `--stats` phase rows. The
prediction to confirm or refute, per phase: P2 ≈ mosaic + canvas + threads ×
149 B × 256 × cols; P2 end ≈ mosaic + 2 canvases; P4 ≈ mosaic + canvas +
16-32 B × check points. If the peaks land within ±25 %, the formulas
extrapolate to the basin.

## 2. The options

**(a) Stream within one process.** Cheap fixes first:
- size resample blocks by node count, not row count (~3 lines, −15 GB);
- adopt the canvas instead of copying it (~3 lines, −8 GB);
- actually drop the target tile before phase 2, as D1 says (~15 lines,
  −8 GB in P4);
- reserve the store's rows from a counting pass (~20 lines with C++, up to
  −12 GB).

Then the stream itself:
- a `SourceWindows` backed by the cache. A window is a windowed `assemble` of
  the sub-plan that meets it, which honours 23a-1's decision 1 and J7.
  `resample` and `check_point_blocks` read through it, so the mosaic never
  exists (~60-90 lines, −8.5 GB).

Total ~120-200 lines in one PR. What remains is the target grid (8.2 GB,
phase 1), then the store (12-14 GB, phase 2), plus the mesh. *Est.* peaks:
50 m ~14 GB, 5 m ~25-30 GB (marginal on 32 GiB), 2 m ~45-50 GB (does not
fit).

Streaming the target grid as well would mean a windowed raster inside the
C++ refine, whose rounds touch triangles anywhere. That is a much larger
change, and it is what 23's pieces avoid. Streaming phase 2's points per
round from the cache (no store) is ~150-250 lines of C++ and Python, and it
decodes the source ~10-12 times.

What (a) keeps of 23: all of it. D1 already names a cached `SourceWindows`
as "a change of provider, not of algorithm", and 23c's per-piece windows
need exactly this provider. **Risks:** decode work goes up (each source block
is read by resample bands and by check blocks; an LRU of decoded blocks
bounds it); 2 m stays out of reach on 32 GiB; it is still one serial phase.

**(b) Increment 23's decomposition as designed** (23b-23g; ~1,555 lines
estimated, ~2,160 at +39 %; plus 23e). Each piece stays within
`--memory-budget` (B14, default 16 GiB), so the basin at 2 m, and below,
fits. It also takes refine's serial phase off the critical path. It keeps
everything, and it is the design. **Risks:**
- the longest road before any basin mesh exists (six PRs, two of them C++
  with vertex removal);
- **`b(T)` was calibrated on the projected path.** The reprojected path adds
  per node the mosaic (~4 B), the copy (~4 B) and the store (16-32 B per
  source node), and the resample transient grows with a piece's columns ×
  threads, not with its area. At coarse tolerances a piece can therefore
  overrun the budget about 2-3× *est.*, unless (a)'s cheap fixes land first
  or `b(T)` is re-measured on the reprojected path.

**(c) Separate sub-basins as a stopgap** (BHO level-3 within 76, or 23e's
DEM-derived units later). About 0 production lines: a script loop over
`rasputin mesh --domain <unit> --out-crs <one basin CRS>`.
- What it gives: N independent files. BHO is an exact coverage (measured in
  23), so the outlines meet without gap or overlap. All units share one
  global lattice (J6) as long as `--out-crs` is the same and
  `default_spacing` rounds to 30 m at every centroid (likely across
  8-21° S; not checked).
- The seams are not conforming: each side splits outline edges on its own,
  in the edge strip and in phase 2 (J9), so there are hanging vertices, and
  z along a seam can differ by up to ~2T.
- What it keeps of 23: the cache and the fetch. It contradicts B1/B13 (BHO
  only for measurement domains), so it is a measurement device, not a
  product.
- **Risks:** the largest Pfafstetter unit sets the peak, and these units are
  unbalanced (23's level-6 table: 93× between smallest and largest). With
  today's ~35-45 B per grid node *est.* (mosaic + canvas + copy + store,
  before the mesh) plus the transients, a unit with a box over roughly
  400-500 M nodes will not fit at any tolerance. Level-3 box sizes were not
  measured; the BHO attributes file on disk covers 63k km², not the basin.

**(d) Other routes.**
- (d1) (a)'s four cheap fixes alone (~40 lines). Basin *est.* ~25 GB at
  50 m, so it likely fits at 20-50 m; with (c) the units get ~2× larger
  headroom.
- (d2) `np.memmap` the mosaic and the target canvas to scratch SSD (~10
  lines). The OS pages them out without growing swap, which is what tripped
  the guard. Speed collapses if refine's working set exceeds RAM, and it
  hides the cause rather than removing it.
- (d3) A ≥ 64-128 GB machine (0 lines). By B14 ("on limited resources, you
  should not expect to reproduce huge catchments") this is legitimate. At
  2 m, ~60 GB *est.* without fixes.
- (d4) A leaner mesh per triangle: 310 B is well above the intrinsic
  ~60-80 B. That is a profile item for 21d, not a quick win.

Prior art: out-of-core and streaming Delaunay construction exists (Isenburg,
Liu, Shewchuk, Snoeyink, SIGGRAPH 2006; Agarwal, Arge, Yi, ESA 2005;
TerraStream, Danner et al., ACM GIS 2007; cited from memory, **unverified**).
None of them is a tolerance-refined TIN checked against a reprojected
source, so (a) and (b) are engineering here, not novelty.

## 3. What a 2-5 m tolerance makes necessary at basin scale

With 2-5 m tolerance as the target, today's code is out at every tolerance,
because its tolerance-independent arrays alone exceed 32 GiB. **(a), or at
least (d1), is necessary at any tolerance**, and it is needed under (b) too:
otherwise `b(T)` underestimates the reprojected path.

At **5 m**, (a) reaches the basin on 32 GiB, but only marginally (~25-30 GB
*est.*; the mesh is ~10 GB). At **2 m** the mesh alone (~30 GB) plus the store
does not fit one process. That needs (b), or (d3), or a store-free streamed
phase 2 plus a leaner mesh. (c) reaches 5 m if the largest unit fits after
(d1), and 2 m only with small units, and in every case without matching
edges.

| option | prod. lines | basin at 5 m / 2 m on 32 GiB *est.* | edges match | keeps 23 | main risk |
|---|---:|---|---|---|---|
| (a) stream in one process | 120-200 | marginal / no | n/a (one mesh) | yes, a prerequisite of 23c | 2 m out of reach; more decoding |
| (b) 23b-23g | ~1,555-2,160 | yes / yes | yes (after 23g) | is 23 | time to first basin mesh; `b(T)` low on this path |
| (c) separate sub-basins | ~0 | largest unit decides / small units only | no (≤ ~2T, hanging vertices) | cache and fetch only | unbalanced units; contradicts B1/B13 |
| (d1) cheap fixes only | ~40 | no / no (fits 20-50 m) | n/a | yes | partial |
| (d3) bigger machine | 0 | yes / yes *est.* | n/a | yes | cost; does not help 32 GiB users |

**Recommendation (the decision is Ola's):**
1. First, @perf's one measurement (§1), to confirm the per-phase model.
2. Then (a), as one PR, with the cheap fixes first. It is needed on every
   route, it makes 23c's per-piece windows real, and it probably gives a
   basin mesh at 5 m. Use (c) only as an interim measurement on DEM-derived
   or BHO units if a basin-wide figure is wanted sooner.
3. Keep (b) as the route to 2 m and to parallel speed, and re-measure
   `b(T)` on the reprojected path after (a).
