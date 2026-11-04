# Increment 137: each triangle's DEM height histogram in the mesh file

Status: **designed** (`@architect`, 2026-11-04, on public master `98fe80c3`;
branch of worktree `dem-stats`). Replaces the plan-area design of commit
`579d193d`: plan area is three vertices away for any reader, so no option
writes it. Full lane: a per-triangle number a downstream reader takes as
input, and a wrong ownership rule gives wrong numbers without a crash. New
C++ lives in `include/terrain/statistics/`; nothing under
`include/terrain/refinement/` or `include/terrain/mesh/` changes, and the
refine path is unchanged with the option off, so no `@perf` acceptance run
is required (`docs/increments/README.md`, *Acceptance*); section 8 asks for
timings, reported and not gated (`perf.md` section 4, an option off by
default).

## 0. The plain picture

- **Today** a mesh file says nothing about the DEM under a triangle except
  through its three corner heights.
- **After**, `rasputin mesh --dem-heights` writes, for each triangle, the
  share of its area at each height above or below the triangle's own
  plane: 23 classes, 1 m wide, centred on −10 … +10 m, with the two end
  classes open (Ola, 2026-11-04: "Cell stats: means simply
  histogramming"). Shares, not pixel counts, so the numbers do not depend
  on the DEM's pixel size once the grid resolves the triangle; a triangle
  holding fewer than 50 valid DEM nodes is flagged unresolved and its
  histogram is all zeros. Without the option every file is byte for byte
  today's.

```
 DEM file --(Python: decode, CRS, NoData)--> DemTile (mesh's CRS, metres)
 mesh as written (Trimmed: x, y, written z) -----------------+
                                                             v
 _core.dem_height_shares(view, vertices, triangles, threads)    [C++]
     per triangle: owned valid nodes -> residual z_node - plane -> class
     -> shares (T, 23) float32, nodes (T,) uint32
                                                             v
 dem_heights.height_shares(...)  ->  HeightShares               [Python]
     area_m2 = nodes * dx * dy;  resolved = nodes >= 50;  unresolved rows -> 0
                                                             v
 write_vtk / write_ply / write_triangles  (stable names, section 3.6)
```

## 1. Prior art: legacy and literature

*Literature.* A triangle's share of area per height class is a hypsometric
(area-altitude) distribution (Strahler 1952, *Hypsometric (area-altitude)
analysis of erosional topography*, GSA Bulletin 63(11): 1117-1142), taken
here per triangle and about the triangle's plane rather than about sea
level, so it holds the relief the mesh smooths away. Sub-grid elevation
classes per coarse cell, as fractions of the cell's area, are the
established way to hand sub-grid relief to a coarse grid (Leung and Ghan
1995, *A subgrid parameterization of orographic precipitation*, Theor.
Appl. Climatol. 52: 95-118, doi 10.1007/BF00865510); they use absolute
elevation bands, this uses fixed bands of height above the plane, because
the triangle's plane already carries its mean height and tilt. Node
ownership is an edge-function test with a tie rule (Pineda 1988, *A
parallel algorithm for polygon rasterization*, Computer Graphics 22(4):
17-20), made exact by a symbolic shift (Edelsbrunner and Mücke 1990,
*Simulation of Simplicity*, ACM TOG 9(1): 66-104). **Departure considered
and not taken:** the land-cover fractions weight each cell by its exact
overlap with the triangle (`include/terrain/land_cover/overlap.hpp`). For
heights that would share a node between triangles and still give a
one-class histogram to a triangle smaller than a pixel; whole-node
ownership is kept (the main session's brief) and the threshold of 3.3
handles the small triangle. No novelty is claimed; none searched.

*Legacy.* `git grep -n -i -E "hist|hypso|rough" legacy-archive -- legacy`
finds no per-triangle DEM statistic; `legacy/rasputin/writer.py` wrote named
per-triangle arrays through meshio. Carried across: a per-triangle quantity
is a named cell array of the mesh file.

*Ideas file.* Takes up the "carry each triangle's sub-grid distribution"
candidate of "A terrain mesh that follows the physics, not only the height"
(heights only; slope and aspect not). Overlaps "A real storage format for
meshes" (23 more numbers per cell for its reader to carry), "Mean, RMS and
percentiles of the vertical error in --stats" and "unbiased height
statistics under the sup-norm" (a resolved triangle's histogram gives its
mean signed misfit to within half a class); takes up neither.

## 2. Measured

Scratch probes, master `98fe80c3`, M1 Max on mains power, on
`tests/fixtures/dem_archive/7908_3_10m_z33.tif` (5,051 × 5,051 nodes,
10 m, EPSG:25833, 42,476 NoData nodes). Ownership as in 3.2; 20,000
triangles sampled per mesh. 95 % of the tile's mesh area is sea, every node
exactly 0 m, so residuals there are exactly 0; figures below are for the
**land triangles** (some residual not 0), except the counts of triangles.

| mesh | triangles | land area in triangles with < 50 nodes | land nodes with \|r\| > 5 m | \|r\| > 10 m | max \|r\| |
|---|---|---|---|---|---|
| `--tolerance 1` | 464,290 | about 82 % | 0 | 0 | 1.00 m |
| `--tolerance 5` | 116,389 | about 38 % | 0 | 0 | 5.00 m |
| `--tolerance 10` | 69,571 | about 18 % | 8.9 % | 0 | 10.0 m |
| `--tolerance 20` | 50,147 | about 7 % | 28 % | 8.1 % | 20.0 m |
| `--stride 30` | 56,594 | 0 | 51 % | 35 % | 195 m |

(The < 50 column interpolates the probe's 25 and 64 node rows.)

**Resolution.** The same triangles measured on the 10 m DEM and on its
2 × 2 block mean (20 m, nodes at the block centres): the total variation
distance between the two share vectors (half the sum of absolute
differences; 0 is equal, 1 is disjoint), 1 m classes (inner range ±5 m in
the probe), by the triangle's
node count on the 20 m grid, medians over tolerances 5, 10 and 20 m:

| 20 m nodes | 9-15 | 16-24 | 25-48 | 49-99 | 100-399 |
|---|---|---|---|---|---|
| median distance | 0.16-0.27 | 0.11-0.20 | 0.07-0.15 | 0.04-0.09 | 0.02-0.05 |
| 90th percentile | 0.29-0.42 | 0.21-0.32 | 0.13-0.23 | 0.09-0.15 | 0.05-0.09 |

The distance falls roughly as 1/√n (median times √n: 0.92 at 9-15 nodes,
0.59 at 100-399, `--tolerance 10`), as the sampling noise of a histogram
of n nodes would. With 0.5 m classes it is 1.5-2.0 times larger in every
row from 25 nodes up.

**Cost.** `--tolerance 10 --rms` on the fixture: the error summary, the same
walk over 25.5 M nodes, 0.037 s at 10 threads; the run 0.229 s. A
one-thread C++ scratch over 25.5 M residuals: the error summary's per-node
term 0.034 s, the same plus a class index and a count 0.051 s, so the new
per-node term is about half the existing walk's.

## 3. Design

### 3.1 The classes

Class k, for k = −11 … +11, holds the residual r (metres) with
(k − ½) ≤ r < (k + ½); class −11 is r < −10.5, class +11 is r ≥ 10.5.
Index k + 11 in every array (0 … 22). The 22 edges −10.5, −9.5, …, +10.5
are exact in double, so the class of r is exact in r. Why:

- **Centred on 0.** Flat ground under a flat triangle (water, the sea of
  section 2) gives residuals of exactly 0 or rounding noise about 0; a class
  edge at 0 would split that noise between two classes. Centred, it lies in
  class 0.
- **1 m wide.** The width bounds what the histogram can say: about the
  vertical error of a 30 m global DEM (several metres) and above a 10 m
  lidar DEM's (decimetres). 0.5 m classes cost 1.5-2.0 times the sampling
  noise of section 2 for detail most DEMs do not hold.
- **±10 m inner range.** A tolerance mesh bounds |r| by its tolerance, so up
  to `--tolerance 10` every node lands in an inner class (section 2); at 20 m
  8 % of land nodes fall into the open classes, at `--stride 30` 35 %.
  The open classes still count them, so shares always sum to 1.
- **Fixed, not tolerance-relative:** one meaning per name across meshes, as
  a reader by name needs. The edges are written into every file (3.6), so a
  later change of classes breaks no reader that reads them.

Constant scale: metres of height; checked up to a 20 m tolerance and a
195 m residual (section 2).

### 3.2 Which DEM nodes a triangle owns

A DEM node is the centre of its pixel. Node p belongs to the triangle that
holds p moved by (ε, ε²) in world axes (x east, y north), ε → 0⁺. For a
counter-clockwise triangle in world axes, p is in when, for each edge a → b,
orient(a, b, p) > 0, or = 0 and the edge is *inclusive*: it runs south
(b.y < a.y), or due east (b.y = a.y and b.x > a.x). In lattice indices
(column right, row down) that is b.row > a.row, or b.row = a.row and
b.col > a.col. A vertex that is a node goes through the same test (two of
its orients are 0). orient is the existing exact predicate; the tie compares
coordinates. Each node in the meshed region belongs to exactly one triangle,
across piece seams too, since the rule reads positions, not triangle numbers
(the scratch brute force recorded in commit `579d193d`: 300 random lattice
triangulations with off-node fans, 0 violations; flipping the column tie or
dropping the edge case gave 3,000 and 10,690). A node-aligned rectangle of
C × R nodes owns exactly (C − 1)(R − 1) nodes: the shift drops its east
column and north row.

A NoData node is owned but not counted. The residual is
r = z_node − plane(p), the plane through the triangle's three corners
**as written** (x, y and the file's z, so a piece's seam heights are the
written ones), evaluated in double.

### 3.3 Resolved or not

n = the triangle's owned valid nodes. **Resolved when n ≥ 50
(`MIN_NODES`).** Then share_k = n_k / n (float32), summing to 1 within
23 · 2⁻²⁴ ≈ 1.4e-6. Unresolved: every share 0.0 and `dem_height_resolved`
0. Why 50: by section 2 a triangle with about 50 nodes on the coarser of two
grids has a median distance between the two histograms of 0.04-0.09 and a
90th percentile of 0.09-0.15, that is, up to about a tenth of its area
moves class between a 10 m and a 20 m DEM; at 25 nodes it is about 1.7
times that. A smaller
triangle on a tolerance mesh is one the plane already fits to within the
tolerance, so a reader loses little by treating it as all class 0 (the flag
says so). The threshold is resolution-dependent by intent: it is the test
"the grid resolves the triangle". Constant scale: nodes per triangle;
checked at 10 m and 20 m grids, up to 1,000 nodes (section 2).

`dem_area_m2` = n · dx · dy (square metres in the plane of the mesh's CRS,
the DEM's pixel area times the valid nodes owned), written for every
triangle, resolved or not. It approaches the plan area as the pixels shrink,
less the NoData, so it does not depend on the pixel size either, and a
reader aggregating triangles weights each histogram by it.

### 3.4 The interfaces (`@tester` calls these)

C++, `include/terrain/statistics/height_shares.hpp`, namespace
`terrain::statistics`:

```cpp
struct HeightClasses { double step_m = 1.0; std::uint32_t half = 10; };  // 2*half + 3 classes
template <raster::RasterSource R>
void height_shares(const R& dem, std::span<const Point2> xy, std::span<const double> z,
                   std::span<const std::array<std::uint32_t, 3>> triangles, HeightClasses classes,
                   std::span<float> shares /* T * K */, std::span<std::uint32_t> nodes /* T */,
                   unsigned threads);
```

Per triangle: the lattice positions of its corners (as
`refinement/measure.hpp`'s `summary_lattice` places them: included, not
edited), `mesh::for_each_row_span` and `raster::for_each_row_segment` for
the closed node set, the tie test of 3.2 at each span's ends, on a flat
edge's row and at node vertices; counts in a local array of K; then its row
of shares (zeros when n = 0) and n. One output row per triangle and no
reduction, so the bytes do not depend on the thread count. Refusals:
`std::invalid_argument` for a vertex outside the node rectangle or not
finite, a triangle that is not counter-clockwise, or spans of the wrong
length. No policy: the threshold is Python's.

Binding (`bindings/core.cpp`, `_core.pyi`), GIL released:

```python
def dem_height_shares(dem: RasterView, vertices: NDArray[float64],  # (M, 3)
                      triangles: NDArray[uint32], *, threads: int = 0,  # (T, 3)
                      ) -> tuple[NDArray[float32], NDArray[uint32]]: ...  # (T, 23), (T,)
```

Python, `src_python/tin_engine/dem_heights.py` (no CLI, no I/O):

```python
STEP_M: Final = 1.0; HALF: Final = 10; MIN_NODES: Final = 50
def edges_m() -> NDArray[float64]                      # (22,): -10.5 ... 10.5
@dataclass(frozen=True, slots=True)
class HeightShares:
    shares: NDArray[float32]     # (T, 23); zero rows where not resolved
    area_m2: NDArray[float64]    # (T,)  nodes * dx * dy
    resolved: NDArray[bool_]     # (T,)
def height_shares(tile: DemTile, mesh: Trimmed, threads: int) -> HeightShares
```

`height_shares` refuses a `tile.meta.geographic` tile (`ValueError`): its
dx · dy would be square degrees.

### 3.5 Where it runs

The tile the mesh was built on, in the mesh's CRS (the DEM, or the target
grid a reprojected or geographic DEM is resampled onto), on the mesh as
written, after the trim:

- `--stride` path: after `trim`, the tile held.
- `_refined`: when `how.dem_heights`, the tile is kept past the trim and
  `height_shares` runs there; `Refined` gains `heights: HeightShares | None`.
  On the final-check path (a resampled grid) this holds the target tile
  through the final check: the tile's bytes, 4 per target node for
  float32, for that phase only. With the option off, nothing moves.
- `--triangles`: trials run with `dem_heights=False`; the kept mesh is
  measured once with the held tile, as its error summary is.
- Cut runs: each piece through `_refined`; `Piece` gains `heights`, and
  `_write_piece` writes them. A piece's window already holds every node its
  triangles cover (refinement scans them), and 3.2 gives a seam node to one
  piece.

`How` gains `dem_heights: bool = False`. With a resampled grid, the
histogram is of the grid the mesh was refined against, not of the source's
own nodes (question 3).

### 3.6 The files (stable names)

| name | `.vtk` | `.ply` | `.gpkg` |
|---|---|---|---|
| `dem_height_shares` | cell array, 23 components, `float` | face `list uchar float`, always 23 | table `dem_height_shares` (`triangle`, `class`, `share`), resolved and non-zero rows only |
| `dem_area_m2` | cell array, `double` | face `double` | `triangles` column `REAL` |
| `dem_height_resolved` | cell array, `unsigned_char` | face `uchar` | `triangles` column `INTEGER` |
| `dem_height_edges_m` | `FieldData` array, 22 `double` | comment `dem_height_edges_m -10.5 … 10.5` | table `dem_height_classes` (`class`, `lower_m`, `upper_m`, NULL at the open ends) |

`.vtk`: every cell array covers every cell; lines get zeros. Placed after
the existing cell arrays. `.ply`: after the existing face properties, binary
and ASCII. Record entries (so in the file's dataset strings):
`dem_height_min_nodes` 50 and `dem_height_resolved_area_share`, the share of
the written plan area in resolved triangles. The names are the contract:
renamed only with a deprecation by a later increment.

**Reader.** `MeshArrays` gains `dem_height_shares` `(T, 23)` float32,
`dem_height_edges_m` `(22,)` float64, `dem_area_m2` `(T,)` float64 and
`dem_height_resolved` `(T,)` bool, each the triangle entries only, or None.
A wrong type, entry count, or component count other than edges + 1 is a
`ValueError` naming the array.

### 3.7 The option, refusals and exclusions

`--dem-heights`, a flag, off by default. Help: "Write, per triangle, the
share of its area at each height about its own plane, in 1 m classes from
−10 to +10 m with open ends (dem_height_shares), from the DEM nodes it
holds; with dem_area_m2 and dem_height_resolved (50 nodes or more)." Works
with `--stride`, `--tolerance`, `--triangles`, `--domain` and cut runs, in
`.vtk`, `.ply` and `.gpkg`; mirrored in the MCP tool's arguments.

Not built: slope and aspect distributions, hollow and knoll shares (they
need a window past the triangle), moments, a tolerance-relative class
width, plan area (three vertices away).

## 4. Invariants

1. **Off is today:** without `--dem-heights`, `.vtk`, `.ply` and `.gpkg`
   bytes equal master `98fe80c3`'s for the same command.
2. **Shares sum to 1 or the triangle is flagged:** resolved rows sum to 1
   within 1e-5 (float32 rounding of 23 terms: at most 23 · 2⁻²⁴ ≈ 1.4e-6,
   a bound, not a measurement);
   unresolved rows are all exactly 0.0.
3. **Each node once:** over all triangles of a run (all pieces of a cut
   run), every valid DEM node whose shifted point lies in the written mesh
   is owned by exactly one triangle, and no other node by any.
4. **The areas add up:** Σ `dem_area_m2` = dx · dy · (the nodes of 3),
   exactly for integer counts times an exact pixel area; for a
   node-aligned rectangle without NoData it equals Σ plan area.
5. **Thread count changes nothing:** the arrays are byte-identical for any
   `--threads`.

## 5. Tests (`@tester`)

Figures are from the exact scratch oracle (`Fraction` ownership by 3.2,
residuals in double) unless said. "Shares" are listed as {class in metres:
count}; share = count / n.

1. **11 × 11 nodes, 10 m, CLI.** `micro_tiff`, z = 2 m on odd columns, 0
   elsewhere, `--stride 10 --dem-heights`, `.vtk` and `.ply`: 2 triangles.
   Whichever diagonal the mesh takes, one triangle has n = 55, shares
   {0: 30, 2: 25} (nw-se diagonal) or {0: 25, 2: 30} (ne-sw), resolved,
   `dem_area_m2` 5,500; the other n = 45, unresolved, all shares 0.0,
   `dem_area_m2` 4,500. Σ `dem_area_m2` = 10,000 = Σ plan area. The test
   reads the diagonal from the file.
2. **NoData.** Test 1 with node (row 5, col 3) NoData (world (30, 50)):
   the 55 triangle becomes n = 54, {0: 30, 2: 24} (nw-se) — or for ne-sw
   the 45 triangle becomes 44, {0: 25, 2: 19}; Σ `dem_area_m2` 9,900.
3. **The threshold, Python level** (`height_shares` with a hand-built
   `Trimmed`, nw-se triangles of test 1): NoData at col 0, rows 5-9 gives
   n = 50, resolved; rows 4-9 gives n = 49, unresolved.
4. **Class edges and open ends, binding level.** 6 × 6 nodes, spacing 1,
   upper-left node at world (0, 5); one triangle (0.5, 0.5), (4.5, 0.5),
   (0.5, 4.5), written z 0 at all corners; DEM 0 except world (1, 1) 3.0,
   (1, 2) 10.5, (1, 3) −10.6, (2, 1) −10.5, (2, 2) −0.5, (3, 1) 0.5, and
   (3, 4) 7.0 outside. n = 6, one count each at indices 14, 22, 0, 1, 11,
   12 (classes +3, open upper, open lower, −10, 0, +1). The four nodes on
   the long edge (x + y = 5, running north-west) are not its own.
5. **Ties and vertices, binding level.** The same grid, z 0 except world
   (1, 1) 3.0, (2, 2) −0.5, (3, 1) 0.5, (3, 4) 7.0; triangles
   (0, 0), (4, 0), (0, 4) and (4, 0), (4, 4), (0, 4), z 0: n = 10
   {0: 9, 3: 1} and n = 6 {0: 5, 1: 1}; 16 = 4 × 4. Vertex (0, 0) is the
   first triangle's.
6. **Plane as written.** Test 5's lower triangle with written z 1.0 at all
   corners: all 10 nodes in class −1, except (1, 1) in +2.
7. **Resolution, binding level.** z = 3 sin(2πx/400) cos(2πy/300) m on
   [0, 1000]² m, at 10 m (101²) and 5 m (201²) nodes; the 32 triangles of a
   250 m grid split nw-se, corners at z(x, y). n is 300-325 at 10 m and
   1,225-1,275 at 5 m; Σ n 10,000 and 40,000 (= 100² and 200²). Classes
   −5 … +5 occur. The distance of section 2 per triangle: measured max 0.029,
   median 0.017; the test asserts max ≤ 0.05.
8. **Each node once, property.** C++ (`tests/cpp/property`): 300 random
   lattice triangulations with off-node fans and the brute-force oracle of
   3.2: every node in the union owned once, none outside; Σ n equals the
   oracle's count.
9. **Cut run.** `pieces_fixtures` `hills()`, `SQUARE`, tolerance 1,
   `--pieces 4 --dem-heights`: 4 piece files, 1,163 triangles; Python oracle
   over all pieces: each valid node in the union owned once (invariant 3),
   Σ `dem_area_m2` = 100 m² × that count.
10. **Off is byte-identical.** Tests 1 and 9 and a `.gpkg` of test 1
    without the flag: sha256 equal to digests recorded at `98fe80c3`
    (constants in the test; that commit named in a comment).
11. **Threads.** Test 7's 5 m case through the binding, and test 9's run,
    at threads 1, 2 and 7: equal bytes.
12. **Formats agree.** Test 1 as `.vtk`, `.ply` and `.gpkg`: the same
    shares, areas and flags; the edges array, comment and table give the
    same 22 edges; `.vtk` lines get zeros; `read_mesh` returns the triangle
    part.
13. **Reader refusals.** A hand-made `.vtk` with `dem_height_shares` of
    type `double`, 22 components, or one cell short, and an edges array of
    21: `ValueError` naming the array; a file without them gives None.
14. **CLI.** `--dem-heights` in `rasputin mesh --help` and the MCP tool's arguments.

**Invariant-critical:** the ownership rule. Mutation targets in
`height_shares.hpp`: each tie branch (south, due east), the vertex case,
and the two class-edge comparisons.

**Mutants** the suite must kill: the column tie flipped (tests 1, 5, 8);
the row tie flipped (1, 5, 8); node vertices left out, as `for_each_row_span`
reports them (5, 8); `<=` for `<` at a class edge (4); classes not centred,
`floor(r)` (4); the open ends not clamped (4); the plane from the DEM, not
the written z (6); NoData counted in n (2, 3); `>` for `>=` at 50 (3); an
unresolved row left non-zero (1); a shared count array across threads (11);
the flag ignored in `_write_piece` (9); anything written when off (10);
lines not zero-padded (12).

## 6. Degeneracies

A triangle owning no node: n = 0, unresolved, area 0. A zero-area
triangle never reaches the file (the engine writes none). A vertex outside
the tile: `invalid_argument` (a defect upstream, never a user's input).
All-NoData triangle: as n = 0.

## 7. Cost and memory at 637,000 km², 30 m

7.1e8 DEM nodes, 28 times the fixture.
- **Time:** the walk is the error summary's plus half again (section 2):
  an estimate of 0.055 s on the fixture and 1.5 s for the whole basin at
  10 threads, spread over its pieces. Prediction for the fixture command at
  `--tolerance 10`: +0.06 s for the walk and +0.03 s to encode 7 MB, on a
  0.23 s run.
- **Memory:** 92 bytes of shares, 4 of node count, 8 of area and 1 of flag:
  105 bytes per triangle while held, about twice that while one file is
  encoded. Per DEM node that is 105 × (triangles per node): on the fixture's
  land at `--tolerance 10` about 0.02 triangles per node, 2 bytes against
  the piece budget's estimate of 37 per node (`decompose.bytes_per_node(10)`); at the
  every-node bound (2 triangles per node) 210 bytes against 637. Plus, with a
  resampled grid, 4 bytes per target node through the final check. The piece
  budget's estimate is not changed; `@perf` reports the peak with and without.
- **Files:** 101 bytes per cell. `--tolerance 10` on the fixture: 7.0 MB on
  a 2.6 MB `.vtk`.

## 8. Acceptance (timings reported, not gated)

On mains power, median of 3, against `98fe80c3`: the fixture command
(`--tolerance 10`, `.vtk`) with and without `--dem-heights`, the same with
`--out-crs EPSG:3035` (the resampled path), and cut by `--pieces 4`; time
and peak resident memory each.

## 9. Size

Estimated +440 counted lines: `height_shares.hpp` +120, binding and
`_core.pyi` +45, `dem_heights.py` +50, `cli.py` +60 (flag, `_refined`,
budget, stride, piece and single writes, record), `refine_job.py` and
`basin_run.py` +10, `io/vtk_legacy.py` +35, `io/ply.py` +30,
`io/geopackage_write.py` +40, `io/vtk_read.py` +45, MCP arguments +5. One PR.

## 10. Questions for Ola

1. Classes: 1 m wide, centred on whole metres from −10 to +10 m, two open
   end classes (23). Default: **this**; 0.5 m classes from −5 to +5 m was
   the other proposal, with more noise and a narrower range.
2. A triangle with fewer than 50 valid DEM nodes is flagged and its
   histogram left empty. Default: **50**.
3. A reprojected or geographic DEM (the large basin's case): measure
   against the grid the mesh is refined on, holding it a little longer.
   Default: **yes**; measuring against the source's own nodes waits for a
   reader who needs it.
