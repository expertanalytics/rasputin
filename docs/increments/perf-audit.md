# Audit: where `rasputin mesh` is slow, and why — a whole-pipeline look

Status: **audit, not a design.** Nothing here is built. Written by
`@architect` on 2026-10-07, against master `91262c73`; every `file:line`
below is pinned to that commit. It answers Ola's question of 2026-10-07:
are there design patterns in the pipeline that are known to be slow, and is
Python doing work that C++ should do? Each finding that is worth building
becomes its own increment with its own design.

## The answer in short

1. **The two slow Swedish phases are one GEOS call, not decoding and not
   JSON.** On Lagan, `decode` (35.7 s) and `features read` (30.1 s) are
   almost all GEOS `buffer` (growing a polygon by a distance) applied to the
   catchment outline, three times. Reading and resampling the DEM takes
   0.7 s; reading the 62 MB GeoJSON takes 1.5 s. The outline is SMHI's
   raster outline: 56,261 vertices, every edge 10 m long and axis-parallel (a
   staircase). GEOS's buffer gets very slow on a staircase once the distance
   is several steps wide (finding F1).
2. **Two of the three calls can go with the mesh byte-identical**, measured on
   Lagan and Numedalslågen: one is the same buffer computed twice, and the
   other only needs the convex hull of the result, which can be grown from the
   hull instead. That takes Lagan from 73.2 s to 27.5 s. A third change, the
   outline grown in pieces of 1,000 edges and then united, gives the same
   region in a fifth of the time and the same `.vtk`; with all three, Lagan
   ran in 17.3 s.
3. **Python is not doing per-element work that C++ should do in any hot
   path.** Resampling, reprojection and the clip are already vectorised NumPy,
   pyproj or GEOS calls. What remains in Python that grows with the mesh is
   three sort-or-loop steps (F4, F5, F6) that cost under a second each today
   and grow to minutes at basin scale. None of the fixes needs the I/O
   boundary (`CLAUDE.md` §2) changed.

## Method

- **Code:** master `91262c73`'s `src_python/tin_engine`, copied to a scratch
  directory and run with the main checkout's `.venv` (Python 3.14.7, shapely
  2.1.2 with GEOS 3.13.1, pyproj 3.8.0) and its compiled `_core`. No C++ was
  built. No C++ changed between that `_core` (built from `ed125121`) and
  `91262c73`: `git diff --stat ed125121 91262c73 -- bindings include src
  CMakeLists.txt` prints nothing. The launcher drops the editable install's
  import finder so the copy loads (`perf-audit-probes/run.py`; it prints
  `tin_engine.__file__`, and it printed the scratch copy every run).
- **Inputs:** Lagan and Ljungan above Flåsjö as in Ola's Swedish runs (GLO-30
  from `rasputin_data/sweden_glo30_cache` resampled to EPSG:3006, CORINE
  GeoJSON in EPSG:3035, `--tolerance 10 --binary`); Numedalslågen as in the
  bottleneck profile (`docs/benchmarks/2026-10-06/bottlenecks/README.md`:
  DTM10, CORINE GeoPackage).
- **Power:** AC for every run here (`pmset -g batt` before and after each:
  "Now drawing from 'AC Power'", 100 %). Another agent's session was running
  at the same time, so every timing is **indicative only**: one run, no
  repeats, and not a baseline for `tools/bench.py`.
- **Identity:** the prototype fixes were checked by the SHA-256 of the
  written `.vtk` against master's on the same input
  (`b2738acf…` Lagan, `34f7117e…` Numedalslågen; equal for every variant
  below).
- **Profiles:** cProfile of the whole command (`run.py --profile`) to find
  the calls, then uninstrumented `--stats` runs for the figures. The probe
  scripts are in `docs/increments/perf-audit-probes/`; they need the scratch
  layout `run.py` describes (the package copy beside them, in `pkg/`).

## What was measured

### Lagan, phase table, master and the prototype fixes (AC, one run each)

| variant | decode | features read | features clip | total | `.vtk` |
|---|---|---|---|---|---|
| master `91262c73` | 35.70 s | 30.07 s | 3.13 s | 73.20 s | reference |
| (1) the same buffer computed once | 18.45 s | 29.98 s | 3.12 s | 55.68 s | identical |
| (1) + (2) features region grown from the hull | 18.58 s | 1.51 s | 3.15 s | 27.48 s | identical |
| (2) + (3) the outline grown in pieces | 8.73 s | 1.50 s | 3.14 s | 17.51 s | identical |

(1), (2) and (3) are the patches in `perf-audit-probes/run_patched.py`. In
the last row (1) was not active: (3) replaces the method (1) wraps, so both
grown regions were computed, in pieces. With (1) as well, decode would lose
one of its two 3.7 s piecewise buffers; about 14 s total, estimated and
not run.

Ola's own run of 2026-10-06 (before #200 and #204) was 88.6 s: decode
40.1 s, features read 30.3 s, clip 13.0 s. The clip has since fallen to
3.1 s (30b, #200).

### Where Lagan's decode goes (`perf-audit-probes/decode_probe.py`, AC)

| step | seconds |
|---|---|
| `target_grid_for`: the domain grown by √2·31 m, mitred (`target_grid.py:107`) | 19.54 |
| the same buffer again (`dem_input.py:213`) | 18.47 |
| `source_region`: that region moved to EPSG:4326 and grown again (`target_grid.py:144-145`) | 0.22 |
| plan, cache check | 0.004 |
| `assemble`: read and decode the six GLO-30 windows, seams | 0.16 |
| `resample`: 17.4 M target nodes, bilinear, 10 threads | 0.55 |

`features read` on Lagan: `feature_input.source_region` 28.9 s in one call
(cProfile), of which GEOS `buffer` is all but milliseconds; JSON parse and
`shape()` 1.45 s.

### GEOS buffer cost depends on the outline's shape

(`perf-audit-probes/buffer_probe.py`, `outline_probe.py`, AC)

| outline | vertices | median edge | mitred, √2·h | round, 100 m | hull, then 100 m |
|---|---|---|---|---|---|
| Lagan, SMHI SVAR (staircase) | 56,261 | 10.0 m, all axis-parallel | 18.60 s | 34.46 s | 0.004 s |
| Ljungan above Flåsjö, SVAR (staircase) | 23,260 | 10.0 m | 2.35 s | 2.38 s | 0.002 s |
| Numedalslågen, NVE | 14,093 | 48.1 m | 0.005 s | 0.018 s | 0.001 s |
| Rio das Velhas, BHO 5k (São Francisco piece) | 7,310 | 101.5 m, none axis-parallel | 0.003 s | | |
| São Francisco basin, BHO level 2 | 77,310 | 107.1 m, none axis-parallel | 0.032 s | | |

On one staircase (Ljungan, `buffer_scaling.py`), the time against the
distance d: 0.013 s at 5 m and 0.011 s at 10 m (up to one step), 0.54 s at
20 m, 2.6 s at 43.8 m, 12.2 s at 100 m. From Ljungan to Lagan (2.4 times the
vertices, same 10 m step, same d) the time grows 8 times. The cost is in how
far the grown outline folds over itself, which grows with the number of steps
inside the distance. Removing collinear vertices does not help
(`simplify(0)` drops 48 of 56,261), and nor does simplifying by 1 m.

### Numedalslågen (DTM10, smooth NVE outline, AC)

Master 7.09 s, with (1) + (2) 6.70 s, `.vtk` identical. Here the buffers cost
milliseconds; what is left is spread out (rows below). Ljungan, master: 9.97 s,
of which decode 5.61 s and features read 2.97 s — the same shape as Lagan.

## Ranked findings

"Cost now" is measured on AC unless marked *est.* (estimated). "São
Francisco" is how the cost grows at the target scale (~637,000 km², 30 m
DEM, increment 23's pieces). Size is production lines, rough.

| # | finding | where (`91262c73`) | cost now | São Francisco | fix pattern | size | needs |
|---|---|---|---|---|---|---|---|
| F1a | The same mitred buffer of the domain computed twice on the reprojected path | `target_grid.py:107`, `dem_input.py:213` | 18.5 s Lagan, ~2.3 s Ljungan, ~0 Norway | depends on the outline (F1c) | compute once, pass the grown polygon to `target_grid_for` (or return it) | ~5 | nothing; output byte-identical (measured) |
| F1b | `feature_input.source_region` grows the whole domain by 100 m with round joins, only to take a convex hull; called up to three times per source (GeoPackage box, region, geographic source's bound) | `feature_input.py:161`; calls at `:294`, `:311`, `:316` | 28.9 s Lagan, ~2.4 s Ljungan, 0.08 s Norway | as F1c | grow the hull, not the domain: conv(P ⊕ B) = conv(P) ⊕ B; compute once per source | ~5 | nothing; `.vtk` byte-identical on Lagan and Numedalslågen (measured) |
| F1c | GEOS buffer of a staircase outline, superlinear in vertices and in distance per step | the one buffer left after F1a: `target_grid.py:107`; also `_domain_plan`, `dem_input.py:250` on the projected path | 18.5 s Lagan after F1a; 3.7 s in pieces | none for BHO outlines (0.03 s for the whole basin); minutes or worse for any raster-traced outline with steps finer than the grid | grow the outline in short overlapping pieces and unite (measured: same region, same `.vtk`), or replace the grown polygon by distance tests on the nodes that need it | ~25 | a short design: what equality the region needs (section "F1c") |
| F2 | Features clip: one GEOS `intersection` per cut chain against the whole, unprepared domain, serially | `feature_input.py:354` | 3.1 s Lagan (56k-vertex domain), 1.8 s Numedalslågen | grows with cut chains × domain vertices; pieces bound the domain | vectorised `shapely.intersection` over all chains on a thread pool (GEOS releases the GIL); or clip in the core's noder (the domain ring is already noded with the lines) | ~15 / ~150 | thread pool: none; noder clip: a design |
| F3 | Refine's mesh crosses to Python as arrays, then `refine_strip` and `refine_points` rebuild it (copy, lattice mesh, adjacency, a full Lawson pass) outside their timers | `bindings/core.cpp:350-378`, `refine_points.hpp:248`, `:265` | `other` row: 0.60 s Numedalslågen, 0.65 s Lagan (and 1.42 s Skiensvassdraget in the 2026-10-06 profile) | linear in triangles, every piece | time it first (a row for it); then keep refine's mesh in its outcome and hand it on | ~60 C++ | design: an addition to the binding surface (`core.cpp`'s rule) |
| F4 | Constraint edges joined to their feature bits with Python sets and dicts, one triangle side and one chain edge at a time | `cli.py:514-557` | 0.34 s Numedalslågen, 0.43 s Lagan | linear in start triangles and constraint edges; *est.* tens of seconds per run with basin-scale features | NumPy: edges from the mask bits, `np.unique` on packed keys, `searchsorted` join; same output | ~30 | nothing; pin by identical `.vtk` |
| F5 | Land-cover components: an `argsort` of every triangle side (3T keys), then array union-find | `landcover.py:52` | `regions` 0.35 s Lagan, 0.54 s Numedalslågen under cProfile (`land cover` row, uninstrumented: 0.47 / 0.70 s) | T log T time and ~5 int64 arrays of 3T: *est.* a minute and several GB at 10⁸ triangles | flood fill on the core's triangle adjacency in C++ (pure compute, within the boundary) | ~80 C++ | design: an addition to the binding surface |
| F6 | ASCII `.vtk` (the default) formatted one number at a time in Python | `vtk_legacy.py:168` | 0.97 s on the Rio das Velhas piece (`write: encode`, `--ascii`; `rasputin_scratch/sao_francisco_piece_meshes/rio_das_velhas_76949_anadem_tol10m_stats.md`) | linear; *est.* minutes at 10⁸ triangles | binary as the default, or a vectorised text encoder with the same round-trip text | ~10 | Ola's call on the default (question 3) |
| F7 | GeoJSON features: whole file parsed by `json`, each geometry built by `shape()` | `feature_input.py:419-425` | 1.45 s for 62 MB, 1.74 M vertices (0.70 parse + 0.75 build); GEOS's own reader is no faster (0.65 s, and the properties still need parsing) | linear in the file, not the domain | read the GeoPackage the windows were cut from (supported: R-tree box query, `io/geopackage.py`) | 0 | Ola's call on the script (question 2) |
| F8 | Small Python loops and repeated preparation | `domain.py:69-70` (`check_extent`, per vertex), `target_grid.py:244` (domain WKB decoded and prepared once per 512² block), `outline.py:72-80` (`catchment`'s ring walk, per segment) | 0.01 s, part of `check points: project` 0.48 s Lagan, not measured | linear; `target_grid.py:244` grows with blocks × domain vertices | vectorise; prepare once per thread | ~10 each | nothing |
| F9 | Noder's edge-property merge in a `std::map` | `include/terrain/noding/node.hpp:374` | inside `start mesh: node` 0.65 s Lagan; its share not measured (no native profile, no build in this audit) | n log n with a node allocation per key | sorted vector of (key, bits), then one reduce | ~15 C++ | measure first |

Checked and not slow: resampling (vectorised, threaded: 17.4 M nodes in
0.55 s); reprojection (pyproj on arrays everywhere; 119 `Transformer`
constructions on Lagan cost 0.16 s in total); the DEM decode itself (six
windows, 0.16 s); every long C++ call releases the GIL (`node`,
`triangulate`, `sample`, `refine`, `refine_points`, `refine_strip`,
`refine_seam`, `constraint_check_points`, `upstream`, `accumulate`,
`reduce_ring`).

## F1 in more detail

**Why it is slow.** GEOS buffers a polygon by offsetting every edge, noding
the result and keeping the outer region (JTS's `BufferBuilder`, Martin
Davis; GEOS's port). On a staircase whose step (10 m) is smaller than the
distance (44 m or 100 m), each step's offset overlaps several neighbours',
so the offset curve has many self-intersections to node and many small
loops to discard. JTS's own input simplifier (`BufferInputLineSimplifier`)
removes only concavities shallower than 1 % of the distance, so 10 m steps
survive it. Recalled, not reread: the 1 % default.

**F1a and F1b are free.** F1a removes a repeat. F1b uses that the convex hull
of a polygon grown by a disc is its hull grown by the disc (Minkowski sum
with a convex set commutes with taking the convex hull). The two differ only
in where GEOS puts the vertices on the arcs: on Lagan the two regions differ
by 11,596 m² of 8,583 km², and both lie within 0.5 m of the true region
(chord error of 8 segments per quarter circle at 100 m). R5's region only has
to contain the domain with its margin, so the swap keeps R5's invariant. The
mesh was byte-identical on both catchments; an increment would add a test
that the region contains the domain grown by the margin less the chord error.

**F1c needs a decision.** The grown polygon from `target_grid.py:107` is used
for the grid's bounds (snapped outward to the 31 m lattice) and, moved to the
DEM's CRS, as the "needed" region: the source nodes that must be covered
(`mosaic._uncovered`), and the overlap nodes the seam report counts. Two
routes:

- *Pieces, then a union.* A polygon grown by d is the polygon united with its
  boundary grown by d. Growing the boundary in pieces of 1,000 edges (each
  overlapping the next by one edge, flat ends) and uniting them with the
  polygon took 3.67 s instead of 18.53 s on Lagan and 1.46 s instead of
  2.25 s on Ljungan, with a symmetric difference of 0 m² and equal bounds
  (`chunked_buffer.py`). With 250-edge pieces it was faster still, but differed
  by 0.019 m² on Lagan and its bounds were not equal — so the equality
  is measured, not proven, and depends on the piece size. A design must say
  whether "the same region to within rounding" is enough (it is for coverage
  and the seam count; the grid's bounds are snapped to the lattice, which
  absorbs rounding unless a bound lies within rounding of a lattice line).
- *No grown polygon.* What the code needs is "which nodes lie within d of the
  domain" and the bounds. Prepared GEOS `dwithin` on node arrays answers the
  first per node against an index of the domain's edges, and the bounds of a
  mitred offset can be taken from the offset vertices directly. This changes
  `mosaic`'s `needed` argument from a polygon to a predicate, a wider change.

Simplify-then-grow-wider (`simplify(ε)`, then `d + ε`) is fast (0.07 s on
Lagan) but did **not** contain the exact region with mitred joins
(`decode_probe.py`: "covers: False" for ε = 31 m and 10 m), so it is not a
safe shortcut as it stands.

**At São Francisco.** With BHO outlines (smooth, ~100 m edges) F1 costs
milliseconds even for the whole basin's 77,310 vertices. It bites on any
outline traced from a raster finer than the DEM, as SMHI's are, and on any
cut piece's outline that keeps such a staircase. Its growth (8 times for
2.4 times the vertices; about 4 times for each doubling of the distance once
it passes one step) means a basin-sized staircase would not finish in a
working day; that is extrapolation, not measured.

## What 30d should look at first

30d is `@perf`'s re-run and profile of Lagan on master. This audit's profile
says decode and features read are the three buffers above; 30d should first
confirm that without cProfile (it costs little: the `--stats` rows with the
patches of `run_patched.py` are in the table above), then take the next
increment's baseline from the remaining rows: `features clip` (F2), `start
mesh: node`, and `other` (F3). The decode fix itself (DEM read, 30c) is not
where Lagan's time goes.

## Prior art: legacy and literature

### Literature

Only F1b, F1c and F5 propose a method; none claims novelty.

- F1b: the convex hull of a Minkowski sum is the Minkowski sum of the convex
  hulls (de Berg, Cheong, van Kreveld and Overmars, *Computational Geometry:
  Algorithms and Applications*, 3rd ed., Springer 2008, ch. 13 on Minkowski
  sums). Recalled, not reread: the chapter.
- F1c, pieces: the buffer as a union of the polygon and the buffered
  boundary, unions computed as a cascaded union (JTS `CascadedPolygonUnion`,
  Martin Davis; `shapely.union_all` calls GEOS's unary union). The JTS buffer
  algorithm and its input simplifier are as named above. Recalled, not
  reread.
- F5: flood fill or union-find over triangle adjacency (Tarjan, *J. ACM*
  22(2), 1975), as 16c and 30a cite; only where it runs moves.

### Legacy

`git grep -n "\.buffer(" legacy-archive -- legacy/rasputin` returns
`legacy/rasputin/geometry.py:259` (a `GeoPolygon` grown by a value) and
`legacy/rasputin/gml_repository.py:204` (a domain grown by `buffer_size`):
plain GEOS buffers, with nothing about cost to carry across.

## Questions for Ola

1. **F1c's route.** Grow the domain in pieces (same region to within
   rounding, measured; about 25 lines), or drop the grown polygon for
   distance tests on the nodes (exact, a wider change to `mosaic`)?
   *Default: pieces, with an increment test that the region is the GEOS
   buffer's to within 1e-6 m on the staircase fixtures.*
2. **CORINE for Swedish runs.** Your script cuts GeoJSON windows from the
   CORINE GeoPackage and passes those; `mesh` can read the GeoPackage itself
   (`--features <gpkg> --features-layer <layer>`, R-tree query of the
   domain's box). Switch the script to the GeoPackage? *Default: yes; no code
   change. Saves about 1.4 s on Lagan; the window step goes.*
3. **ASCII as the default `.vtk`.** The default writes text, one number at a
   time in Python (F6); at basin scale that is minutes. Make `--binary` the
   default, or keep text and vectorise its encoder? *Default: keep text as
   the default and vectorise the encoder when a basin run shows it.*

No finding needs the I/O boundary changed: F3 and F5 move pure computation
into the core, and F7 changes only which file is read.
