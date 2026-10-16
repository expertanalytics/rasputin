# Where `rasputin mesh` spends its time on two real catchments

@perf, 2026-10-06, master at `8199f30`. This is a profile for @architect's
design of the fixes, not an acceptance run: production code is unchanged, and
every rewrite timed below exists only in the scratch scripts in `scripts/`.

## Method

- **Machine:** Apple M1 Max (8 performance + 2 efficiency cores), 32 GB,
  macOS 27.0. Python 3.13.15, shapely 2.1.2 (GEOS 3.13.1), numpy 2.5.3,
  tifffile 2026.9.20, imagecodecs 2026.8.16.
- **Power: battery for every run** (`pmset -g batt` before and after each
  `--stats` run: `raw/stats/*_power.txt`).
  The percentages and the profiles hold on battery. Absolute seconds compare only
  with other battery runs. The `--stats` figures from earlier this morning
  (`../rasputin_scratch/norway/*/*_stats.md`) have no recorded power state,
  so they are quoted, not compared.
- **Build:** `uv venv --python 3.13 && uv pip install ".[codecs]"` in the
  worktree: scikit-build-core Release, bounds checks on (libc++ fast mode),
  as `rasputin mesh` prints. The same build mode as the earlier runs.
- **Command:** `scripts/stats.sh <catchment> <repeat>`. It runs `rasputin mesh` with
  DEM `rasputin_data/DTM10_UTM33_20260925`, domain
  `rasputin_scratch/norway/<c>/<c>_outline_nve.geojson` (EPSG:25833),
  `--features rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer
  corine2018 --features-map corine --tolerance 10 --binary`, with 10 threads.
  There were 3 repeats per catchment, run alternately, and the tables give the median.
- **Profiles:** cProfile over the whole command
  (`raw/profiles/*_cprofile_cumulative.txt`), line_profiler over the hot
  functions (`scripts/lineprof.py`, `raw/profiles/*_lines.txt`), and then
  uninstrumented `perf_counter` timings of each step, in scratch copies of the
  phase code (`scripts/clip_bench.py`, `lc_bench.py`, `decode_bench.py`,
  medians of 3 or 5). The line_profiler figures inflate Python loops
  (`_runs` shows 3.0 s there against 0.6 s uninstrumented), so the savings
  below come from the uninstrumented scripts only.
- **Identity checks:** every rewrite is compared with the production function's
  output on the same input: the clipped lines as WKB bytes per feature, in
  order; the land-cover code array and its three stderr counts; the
  assembled DEM canvas (`array_equal`, NaN equal) and the seam report. These
  checks hold on these two inputs and are not a proof for every input.
- Meshes (`.vtk`, 58 and 134 MB) are not committed. They were written to the session scratchpad and
  deleted. Rerun `scripts/stats.sh` to regenerate them. The scripts write to
  a scratchpad path set at their top (`S=`), so change it before rerunning.

## 1. The phase table, reproduced (median of 3, battery)

| phase | Numedalslagen s | share | Skiensvassdraget s | share |
|---|---|---|---|---|
| decode | 2.048 | 14.3 % | 1.331 | 6.3 % |
| features read | 0.433 | 3.0 % | 0.279 | 1.3 % |
| **features clip** | **5.833** | **40.6 %** | **5.679** | **26.8 %** |
| start mesh (build + node + triangulate + constraint edges) | 0.980 | 6.8 % | 1.963 | 9.3 % |
| refine | 1.161 | 8.1 % | 3.153 | 14.9 % |
| (of it) refine: start quality | 0.537 | 3.7 % | 1.140 | 5.4 % |
| edge strip (generate + scan + split) | 0.149 | 1.0 % | 0.287 | 1.4 % |
| trim | 0.051 | 0.4 % | 0.119 | 0.6 % |
| **land cover** | **2.868** | **20.0 %** | **6.694** | **31.6 %** |
| write (encode + disk) | 0.031 | 0.2 % | 0.104 | 0.5 % |
| other | 0.612 | 4.3 % | 1.474 | 7.0 % |
| **total** | **14.371** | | **21.174** | |

Per-repeat values: `raw/stats/`. Earlier this morning (power state not
recorded): 14.7 s and 22.4 s; clip 5.9 and 5.8 s, land cover 3.0 and 7.0 s, decode 2.2
and 1.6 s. Output sizes match those runs: 1,287,334 and 3,029,295
triangles.

## 2. Features clip (`feature_input._Tally._take`)

### Input sizes

| | Numedalslagen | Skiensvassdraget |
|---|---|---|
| domain area / outline vertices | 5,548 km² / 14,093 | 10,807 km² / 16,008 |
| hull of the domain (`source_region`) area / its bounding box | 12,887 / 28,075 km² | 13,103 / 20,088 km² |
| CORINE features read (R-tree query of the hull's box) | 7,932 | 5,522 |
| of them meeting the hull | 3,297 | 3,113 |
| of them meeting the domain (kept) | 1,611 | 2,628 |
| rings read / vertices read | 29,183 / 5,417,093 | 26,156 / 4,997,951 |
| segments tested against the hull, one by one | 5,387,910 | 4,971,795 |
| rings with no segment kept | 24,722 | 21,880 |
| chains out of `pre_clip` | 4,672 (603,400 vertices) | 4,445 (581,865 vertices) |
| chains the domain already covers / chains it cuts | 1,760 / 2,912 | 3,268 / 1,177 |

The clip time follows the vertices read (5.4 M and 5.0 M), not the catchment
area. The vertices read follow the hull's bounding box, which is 5.1× the domain
for the long, thin Numedalslagen and 1.9× for Skiensvassdraget. This is why the
clip time stays about the same when the area doubles. The query box is
also widened by each row's own width plus height (`query_features`, R5 "Long
edges"), so it reads more again. That last share was not measured on its own.

### Where the time goes (uninstrumented, `raw/clip_*.json`)

| step | Numedalslagen s | Skiensvassdraget s |
|---|---|---|
| `pre_clip`, all features | 3.54 | 2.70 |
| (of it) `shapely.intersects(segments, hull)`, 29k / 26k calls | 1.73 | 1.07 |
| (of it) segment `LineString`s built per ring | 0.92 | 0.84 |
| (of it) `_runs`, a Python loop over every vertex of every ring | 0.61 | 0.57 |
| `shapely.intersection(chain, domain)`, on chains the domain covers | 1.14 | 2.32 |
| `shapely.intersection(chain, domain)`, on chains it cuts | 0.87 | 0.46 |
| the rest (classify, append) | 0.22 | 0.13 |
| total of the scratch copy (the phase measures 5.83 / 5.68) | 5.75 | 5.58 |

The intersection is GEOS overlay against the full, unprepared domain (14k and 16k
vertices) for each chain, so a chain that lies wholly inside still pays for
the whole outline.

### Rewrites tried (identical WKB, every feature, both catchments)

| variant | Numedalslagen s | Skiensvassdraget s |
|---|---|---|
| production | 5.75 | 5.58 |
| A: endpoint test first (`intersects_xy` with the prepared hull), exact segment test only where neither end is inside; `_runs` returns at once when no edge is kept | 5.28 | 5.12 |
| A + one prepared `hull.intersects(feature)` first, skipping `pre_clip` for a feature that misses the hull | 3.23 | 3.65 |
| A + no `intersection` when the prepared domain `covers` the chain (the chain is kept as is) | 4.14 | 2.83 |
| **A + both** | **2.11** | **1.48** |

## 3. Land cover (`landcover.label_triangles`)

### Input sizes (`raw/lc_*.json`)

| | Numedalslagen | Skiensvassdraget |
|---|---|---|
| triangles / triangle sides | 1,287,334 / 3,862,002 | 3,029,295 / 9,087,885 |
| constraint edges | 156,919 | 275,572 |
| components (areas between lines) = points looked up | 1,860 | 2,881 |
| coded polygons (whole, unclipped CORINE polygons) | 1,611 | 2,628 |
| polygon vertices: total / largest polygon / median | 1,056,147 / 370,070 / 70 | 1,372,887 / 370,070 / 69 |
| point–polygon pairs whose boxes meet | 9,574 | 14,285 |
| polygon vertices walked by those tests (sum over pairs) | 920 M | 1.38 G |

### Where the time goes, and rewrites (uninstrumented, identical codes and counts)

| step | production s (N / S) | rewrite | rewrite s (N / S) |
|---|---|---|---|
| `_lookup`: `STRtree.query(points, predicate="intersects")` | 1.42 / 2.16 | candidates by box, then `shapely.prepare` on the candidate polygons and a vectorised `intersects` (same pairs; a tree of points queried with the polygons is as fast) | 0.09 / 0.12 |
| best triangle per component: `np.lexsort` over all triangles | 0.31 / 1.06 | `np.maximum.at` for each component's largest r, then `np.minimum.at` for the lowest index at it | 0.017 / 0.047 |
| `regions`: `np.isin(keys, cut keys)` then stable `argsort` | 1.01 / 3.07 | one default `argsort` of all keys, then cut keys removed with `searchsorted` | 0.52 / 1.43 |
| **`label_triangles`, whole** | **2.88 / 6.60** | **all three** | **0.68 / 1.79** |

The lookup tests each point against an unprepared polygon, which walks all its
vertices: one 370,070-vertex polygon is a candidate for many points. In both
catchments, what remains after the rewrites is mostly the `argsort` in
`regions`. Inferred, not measured: the C++ core already has triangle
adjacency, and a flood fill there would avoid the sort.

## 4. Decode (`dem_input.open_dem`)

### Input sizes

| | Numedalslagen | Skiensvassdraget |
|---|---|---|
| tiles in the DEM directory (every header read) | 254 | 254 |
| tiles in the plan | 16 | 12 |
| canvas | 17,882 × 15,667 float32, 1,121 MB | 13,675 × 14,653 float32, 802 MB |
| overlap nodes tested against the domain for the seam report | 5,226,633 | 3,066,621 |
| seams reported | 6 pairs, 101,520 nodes | none |

### Where the time goes, and rewrites (median of 5, interleaved; `raw/decode_*.json`)

| step | production s (N / S) | rewrite | rewrite s (N / S) |
|---|---|---|---|
| read 254 headers (`footprints`) | 0.087 / 0.085 | (none tried) | |
| domain plan | 0.105 / 0.076 | (none tried) | |
| decode the windows (`decode_window`, tiles one after another, 4 threads per tile) | 0.875 / 0.651 | 10 threads per tile | 0.482 / 0.361 |
| `_covered`: every overlap node tested against the domain, only for the seam report | 0.411 / 0.246 | test only the nodes where both tiles are valid and the gap reaches the seam threshold, the only nodes `_seam` counts | ~0 |
| `_decide`, `_seam`, canvas fill and copies | ~0.2 / ~0.15 | | |
| **`open_dem`, whole** | **1.607 / 1.143** | **both** | **0.868 / 0.611** |

Every run gave an identical canvas (sha256 prefix `a186d215264b3b85` and
`784d8c7a700e3461`) and identical seams. Not explained: the CLI's
`decode` row (2.05 / 1.33 s) is higher than `open_dem` measured on its own
(first call 1.72 s for Numedalslagen, then 1.60 s; `raw/firstcall.txt`).
Importing tifffile and imagecodecs early does not change it.

## 5. `other` and `start quality` (cheap checks only)

- **`other`** (0.61 / 1.47 s): cProfile shows `_core.refine_strip` at
  0.59 / 1.49 s in total, while its own timers book only 0.025 / 0.063 s
  (`edge strip: scan` + `split`). The difference, 0.57 / 1.42 s, has no row and
  lands in `other`, which covers almost all of `other`. It is thought to be
  the C++ setup and output around the strip run (`start_mesh` in
  `bindings/core.cpp` copies the triangles into a new mesh). A native profile
  has not been taken.
- **`refine: start quality`** (0.54 / 1.14 s) is inside C++ `refine`. Not
  profiled here, because it needs a native profiler.

## 6. Ranked for @architect

Savings are measured on the scratch copies, all with identical output. The
whole-run figure is a projection by subtraction, not an end-to-end run.

| rank | phase | fix | saves (N / S) | must stay byte-identical |
|---|---|---|---|---|
| 1 | land cover | prepared polygons in `_lookup`; group max instead of `lexsort`; one `argsort` in `regions` | 2.2 s / 4.8 s | `land_cover_code` per triangle, and the stderr counts (outside, overlapped, thin) |
| 2 | features clip | whole-feature hull test; endpoint prefilter and early exit in `pre_clip`/`_runs`; no `intersection` for a chain the domain covers | 3.6 s / 4.1 s | each kept feature's lines (vertices, order), and so the noded graph and the mesh; the kept / cut / outside counts |
| 3 | other (`refine_strip` setup and output) | to be profiled first; no row for it today | up to 0.57 s / 1.42 s | the mesh |
| 4 | decode | more decode threads per tile (or tiles in parallel, which needs a memory check); seam mask on qualifying nodes only | 0.74 s / 0.53 s | the canvas (every node), the seam report |

Projected totals with ranks 1, 2 and 4: about 7.8 s for Numedalslagen (from 14.4 s) and
11.7 s for Skiensvassdraget (from 21.2 s). After that, refine and the start mesh
lead.

Points for the design to settle, not settled by these runs:

- Skipping `intersection` for a covered chain keeps the chain object itself.
  GEOS does not promise that `intersection` returns a covered line unchanged,
  so the identity above is an observation on these inputs.
- The endpoint prefilter assumes `intersects_xy` on a vertex and
  `intersects` on a segment agree for a vertex on the hull's boundary. They
  agreed here.
- 10 decode threads were tried on battery on a 10-core machine. The right count,
  and whether tiles can decode in parallel within the memory budget
  (`assemble`'s docstring: one load peaks at 2 to 3 tiles), is open.
