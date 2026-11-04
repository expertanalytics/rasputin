# Increment 137: each triangle's DEM height histogram in the mesh file

Status: **designed, revised** (`@architect`, 2026-11-04, on public master
`98fe80c3`; branch of worktree `dem-stats`). The design review of `6f0e20df`
returned CHANGES REQUESTED: its 50-node threshold left 72 % of Bygdin's
triangles without a histogram. This revision weights each DEM pixel by its
exact overlap with the triangle (the main session's ruling between the two
reviewers' fixes) and reads a reprojected DEM's own pixels. Full lane: a
per-triangle number a downstream reader takes as input, wrong without a
crash if the overlap or the plane is wrong. New C++ lives in
`include/terrain/statistics/`; `include/terrain/land_cover/overlap.hpp` is
refactored so both uses share one clip; nothing under
`include/terrain/refinement/` or `include/terrain/mesh/` changes and the
refine path is unchanged with the option off, so no `@perf` acceptance run
is required (`docs/increments/README.md`, *Acceptance*); section 8 asks for
timings, reported and not gated (`perf.md` section 4, an option off by
default).

## 0. The plain picture

- **Today** a mesh file says nothing about the DEM under a triangle except
  through its three corner heights.
- **After**, `rasputin mesh --dem-heights` writes, for each triangle, the
  share of its area at each height above or below the triangle's own
  plane: 23 classes, 1 m wide, centred on −10 … +10 m, the two end classes
  open (Ola, 2026-11-04: "Cell stats: means simply histogramming"), and the
  share of its area over NoData. Each DEM pixel counts with its exact area
  inside the triangle, so every triangle gets a histogram, however small,
  and the 24 numbers sum to 1. A reprojected or geographic DEM is read on
  its own pixels, not on the resampled grid. Without the option every file
  is byte for byte today's.

```
 mesh as written (Trimmed: x, y, written z, in the mesh's CRS)
 DEM windows (SourceWindows: the tile, or the source DEM on the reprojected path)
        |
        v  [Python, dem_heights.height_shares]
 corners -> the DEM's index space (reprojected first when the CRSs differ),
 per chunk: the window under its corners, the area of one pixel per row
        |
        v  [C++, _core.dem_height_shares]
 per triangle: clip to each pixel -> piece area, piece centroid
   -> residual z_pixel - plane(centroid) -> class; NoData -> its own share
   -> shares (T, 23) float32, nodata (T,) float32
        |
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
the plane already carries the triangle's mean height and tilt. The weights
are the exact overlap of each pixel with the triangle, the construction of
conservative remapping and of `exactextract`'s coverage fractions, by the
clip of Sutherland and Hodgman 1974 (*Reentrant polygon clipping*, CACM
17(1): 32-42), all cited and compared in `26-land-cover.md` section 1 for
`overlap.hpp`, which this reuses. **Departures:** (a) a pixel is read as
one height over its whole square (area-registered); on a slope s a pixel
cut by an edge then differs from the plane by up to s times half a pixel
at its far corner, which a point reading would not show; (b) the residual
of a cut pixel is taken at its piece's centroid, the mean of the plane over
the piece, not at the pixel centre (section 2 measures the difference);
(c) on the reprojected path the triangle is carried into the source's index
space by its corners, with straight edges, and its plane by barycentric
weights, which assumes the map is affine across one triangle (section 2
measures the error). No novelty is claimed; none searched.

*Legacy.* `git grep -n -i -E "hist|hypso|rough" legacy-archive -- legacy`
finds no per-triangle DEM statistic; `legacy/rasputin/writer.py` wrote named
per-triangle arrays through meshio. Carried across: a per-triangle quantity
is a named cell array of the mesh file.

*Ideas file.* Takes up the "carry each triangle's sub-grid distribution"
candidate of "A terrain mesh that follows the physics, not only the height"
(heights only; slope and aspect not). Overlaps "A real storage format for
meshes" (24 more numbers per cell for its reader to carry) and "Mean, RMS
and percentiles of the vertical error in --stats" (the histogram gives a
triangle's mean signed misfit to within half a class); takes up neither.

## 2. Measured

Scratch probes on master `98fe80c3`, M1 Max on mains power, a C++ probe
built on `overlap.hpp`'s clip with `-O3` and the build's hardening, and an
exact `Fraction` oracle for the test figures. Residual classes as in 3.1.

**Residuals at the nodes.** `tests/fixtures/dem_archive/7908_3_10m_z33.tif`
(5,051 × 5,051 nodes, 10 m, EPSG:25833, 95 % sea at exactly 0 m), land
triangles only:

| mesh | triangles | land nodes with \|r\| > 5 m | \|r\| > 10 m | max \|r\| |
|---|---|---|---|---|
| `--tolerance 1` | 464,290 | 0 | 0 | 1.00 m |
| `--tolerance 5` | 116,389 | 0 | 0 | 5.00 m |
| `--tolerance 10` | 69,571 | 8.9 % | 0 | 10.0 m |
| `--tolerance 20` | 50,147 | 28 % | 8.1 % | 20.0 m |
| `--stride 30` | 56,594 | 51 % | 35 % | 195 m |

**Bygdin.** `--dem DTM10_UTM33_20260925 --domain
docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson --tolerance 10
--land-cover corine2018_dtm10_utm33.gpkg --land-cover-classes corine
--edges all`: 54,429 triangles on a 2,195 × 3,314 node tile, no NoData.
71.8 % of the triangles have under 50 pixels of plan area, holding 18.5 %
of it. Σ of each triangle's pixel overlaps equals its plan area to
1.6e-12 (relative, largest); rows sum to 1 within 2.4e-7.

**Cost against the node walk.** The same Bygdin tile and mesh, median of 5:
the error summary's walk (`refinement::error_summary`, 3,025,569 nodes)
against the overlap histogram:

| threads | node walk | overlap walk | ratio |
|---|---|---|---|
| 10 | 13.7 ms | 16.4 ms | 1.2 |
| 1 | 75.0 ms | 117 ms | 1.56 |

Peak resident memory 40 MB and 39 MB: the inputs; the walk allocates 23
doubles per thread and writes 96 bytes per triangle. On the São Francisco
piece below (1,858,576 triangles, 30.9 M nodes, 10 threads): 0.42-0.47 s
against 0.35-0.46 s, and reprojecting its 933,232 vertices 0.05 s.

**Centroid or pixel centre** (departure b), Bygdin, the distance between
the two histograms (total variation: half the sum of absolute differences;
0 equal, 1 disjoint), median by plan area in pixels: under 1, 1.00; 1-4,
0.70; 4-16, 0.29; 16-50, 0.10; 50 and over, 0.03. Area-weighted over the
catchment, the share outside class 0 moves from 78.8 % to 78.1 %. A plane
evaluated at a centre outside a small triangle is an extrapolation; the
centroid's is the piece's own mean.

**Resampled grid or the source's pixels.** The reprojected path refines
against a bilinearly resampled grid. Histograms from it against histograms
from the source's own pixels (the nodes `final_check.stored` files), same
mesh, area-weighted over triangles:

| fixture | triangles | distance | spread, resampled / source | class-0 share, resampled / source | resampled narrower in |
|---|---|---|---|---|---|
| `tests/fixtures/velhas` ANADEM, `--out-crs EPSG:31983 --tolerance 1` | 7,386 | 0.25 | 0.674 / 0.705 m | 0.604 / 0.556 | 52 % of triangles |
| same, `--tolerance 5` | 969 | 0.15 | 1.219 / 1.289 m | 0.225 / 0.213 | 70 % |
| São Francisco piece `bho2017_5k_76949_anadem_window_epsg4674.tif`, `--tolerance 5` | 1,858,576 | 0.24 | 1.347 / 1.463 m | 0.216 / 0.194 | 68 % |

(Spread: the standard deviation of the class, in metres.) The sign is the
smoothing's: the resampled histograms are narrower, by 4-8 % in spread,
with more area in class 0; the mean class is unchanged (0.116 m both on
the São Francisco piece). **The faithful path is chosen** (3.5).
Departure (c) measured on Velhas: a triangle's centroid mapped to the
source's index space lies at most 1.2e-3 m (`--tolerance 1`, edges up to
474 m) and 3.7e-3 m (`--tolerance 5`, up to 849 m) from the mean of its
mapped corners.

**Resolution.** z = 3 sin(2πx/400) cos(2πy/300) m on [0, 1000]² m, 32
triangles of a 250 m grid: overlap histograms from a 10 m and a 5 m DEM
differ by 0.011 (median) and 0.017 (largest). Node-count histograms on the
fixture, 10 m against its 20 m block mean, differed by a median of 0.04-0.09
at 49-99 nodes and 1.5-2.0 times that with 0.5 m classes (commit
`6f0e20df`, section 2); this is the class-width argument of 3.1.

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
- **1 m wide.** About the vertical error of a 30 m global DEM (several
  metres) and above a 10 m lidar DEM's (decimetres). 0.5 m classes cost
  1.5-2.0 times the sampling noise of section 2 for detail most DEMs do not
  hold.
- **±10 m inner range.** A tolerance mesh bounds |r| at the nodes by its
  tolerance, so up to `--tolerance 10` every node lands in an inner class
  (section 2); at 20 m 8 % of land nodes fall into the open classes, at
  `--stride 30` 35 %. The open classes still count them.
- **Fixed, not tolerance-relative:** one meaning per name across meshes, as
  a reader by name needs. The edges are written into every file (3.6).

Constant scale: metres of height; checked up to a 20 m tolerance and a
195 m residual (section 2).

### 3.2 What a triangle covers

The DEM's **index space**: x is the column and y the row, pixel (r, c) is
the unit square [c, c + 1] × [r, r + 1], and node (r, c) is its centre, as
`overlap.hpp` and `land_cover_fractions.index_space` already have it. A
mesh vertex goes there by `col = (x − x_min) / dx + ½`, `row = (y_max − y)
/ dy + ½` when the mesh is in the DEM's CRS, and by
`land_cover_fractions.index_space` (one reprojection to the DEM's CRS)
when it is not.

Per triangle, `overlap.hpp`'s clip gives every pixel's piece inside it:
whole interior pixels without a clip, cut pixels as polygons. Each piece
has an area a (index units) times the row's pixel area A_row (m²:
dx · dy projected, the ellipsoidal cell area of the row for a geographic
DEM, `land_cover_fractions.row_areas`'s rule). Its residual is
r = z_pixel − plane(q), with q the piece's centroid (the pixel centre for
a whole pixel) and plane the affine function through the triangle's three
corners in index space with their **written** z (so a piece's seam heights
are the written ones). In the DEM's own CRS that is the triangle's plane;
on the reprojected path it is the plane by barycentric weights (departure
c). A pixel the view's NoData names, or a NaN, adds its weight to the
NoData share and to no class.

No pixel is owned: a pixel cut by an edge gives each triangle its own
piece, and the pieces of a pixel tile it to rounding, across piece seams
too, since each crossing is computed from the edge's endpoints in one
canonical order (`overlap.hpp`'s header).

### 3.3 The shares

W = Σ a · A_row over the triangle's pieces. share_k = (weight in class k) /
W and nodata = (weight on NoData) / W, each float32; the 23 shares and
nodata sum to 1 within 24 · 2⁻²⁴ ≈ 1.4e-6 (a bound). W > 0 for every
triangle of positive area inside the window, so every row is defined:
**no threshold, no empty row**. A triangle whose pixels are all NoData has
nodata = 1.0 and every share 0.0; the run record counts such triangles
(`dem_height_triangles_without_data`, 3.6), and `--stats` shows the count.
A reader wanting the heights of the valid part divides the shares by
1 − nodata; a reader aggregating triangles weights each by its plan area
times 1 − nodata.

A triangle of a pixel or less has a one- or two-class histogram: that is
what the DEM says at its resolution, not a defect.

### 3.4 The interfaces (`@tester` calls these)

`overlap.hpp` refactor: `detail::overlap_one` becomes
`detail::for_each_piece(tri, rows, cols, visit)`, calling
`visit(r, c, const Polygon* piece)` with `nullptr` for a whole pixel;
`class_areas` passes a visitor doing today's arithmetic in today's order,
so its bytes do not change. A helper `detail::area_centroid_local(const
Polygon&)` gives a piece's area and centroid in the local coordinates
`area_local` uses.

C++, `include/terrain/statistics/height_shares.hpp`, namespace
`terrain::statistics`:

```cpp
struct HeightClasses { double step_m = 1.0; std::uint32_t half = 10; };  // 2*half + 3 classes
struct HeightShares { std::vector<float> shares; /* T*K */ std::vector<float> nodata; /* T */ };
template <class T>  // float or double, the DEM's promoted dtype
[[nodiscard]] HeightShares height_shares(std::span<const Point2> corners,  // 3 per triangle, index space
                                         std::span<const double> z,        // 3 per triangle, written
                                         const raster::RasterView<T>& dem, std::span<const double> row_area,
                                         HeightClasses classes, unsigned threads);
```

Per triangle: the canonical rotation `class_areas` applies (corners and
their z together, so the vertex order changes nothing), the plane, the
pieces, a local array of K + 1 doubles, then its row. One output row per
triangle and no reduction, so the bytes do not depend on the thread count.
Refusals, as `class_areas`'s and before any work: corner and z spans not 3
per triangle, row_area not one per row, a non-finite corner or z, a
zero-area triangle (`std::invalid_argument`), a corner past the window
(`land_cover::OutsideWindow`).

Binding (`bindings/core.cpp`, `_core.pyi`), GIL released, the
`class_areas` binding's checks:

```python
def dem_height_shares(corners: NDArray[float64],   # (T, 3, 2) index space of `dem`
                      z: NDArray[float64],         # (T, 3)
                      dem: NDArray[float32 | float64], nodata: float | None,  # (rows, cols)
                      row_area: NDArray[float64],  # (rows,)
                      *, threads: int = 0) -> tuple[NDArray[float32], NDArray[float32]]: ...  # (T, 23), (T,)
```

Python, `src_python/tin_engine/dem_heights.py` (no CLI, no file I/O):

```python
STEP_M: Final = 1.0; HALF: Final = 10; CHUNK: Final = 2**18
def edges_m() -> NDArray[float64]                      # (22,): -10.5 ... 10.5
@dataclass(frozen=True, slots=True)
class HeightShares:
    shares: NDArray[float32]     # (T, 23)
    nodata: NDArray[float32]     # (T,)
    without_data: int            # rows with nodata == 1.0
def height_shares(windows: SourceWindows, vertices: NDArray[float64], triangles: NDArray[uint32],
                  *, mesh_crs: str | None, threads: int, chunk: int = CHUNK) -> HeightShares
```

`mesh_crs` None: the mesh is in `windows.meta`'s CRS, no reprojection;
otherwise the vertices are reprojected from it. Per chunk of `chunk`
triangles: the pixels under its corners (floor of the least to ceil of the
greatest, clipped to the raster; a window past a `HeldWindows` part is its
`ValueError`, a defect), `windows.window(...)`, the corners shifted to it,
row areas for it. `row_areas`'s body moves to a helper over (top edge,
dx, dy, rows, CRS) that it and `dem_heights` call.

### 3.5 Where it runs

On the mesh as written, after the trim, against:

- **the DEM's CRS** (no `--out-crs`, or the DEM's own): the tile the mesh
  was refined on, `TileWindows(tile)`, `mesh_crs=None`. `_refined` keeps
  `tile` past the trim when `how.dem_heights`; the `--stride` path has it.
- **the reprojected path** (`--out-crs` another CRS, every geographic DEM):
  the **source's** windows, never the resampled grid (section 2). `DemInput`
  and `PieceDem` gain `source: SourceWindows | None`, the `TileWindows` or
  `HeldWindows` they already build for the check points, set only when the
  option is on; `_refined` and the stride path take it, `mesh_crs` the
  target CRS. Today that source tile is freed once the check-point store is
  built (as the code reads: the exhausted block generator held its last
  reference); with the option it is held through the final check.
- `--triangles`: trials run with `dem_heights=False`; the kept mesh is
  measured once, as its error summary is.
- Cut runs: each piece through `_refined`; `Piece` gains `heights`, and
  `_write_piece` writes them. A piece's tile, or its `HeldWindows` part
  snapped outward to source nodes, holds every pixel its triangles touch.

`How` gains `dem_heights: bool = False`. With the option off nothing is
held and nothing moves.

### 3.6 The files (stable names)

| name | `.vtk` | `.ply` | `.gpkg` |
|---|---|---|---|
| `dem_height_shares` | cell array, 23 components, `float` | face `list uchar float`, always 23 | table `dem_height_shares` (`triangle`, `class`, `share`), non-zero rows only |
| `dem_nodata_share` | cell array, `float` | face `float` | `triangles` column `REAL` |
| `dem_height_edges_m` | `FieldData` array, 22 `double` | comment `dem_height_edges_m -10.5 … 10.5` | table `dem_height_classes` (`class`, `lower_m`, `upper_m`, NULL at the open ends) |

`.vtk`: every cell array covers every cell; lines get zeros. Placed after
the existing cell arrays. `.ply`: after the existing face properties, binary
and ASCII. Record entry (so in the file's dataset strings and `--stats`):
`dem_height_triangles_without_data`. No NaN is written, so ASCII readers
need no special value. The names are the contract: renamed only with a
deprecation by a later increment.

**Reader.** `MeshArrays` gains `dem_height_shares` `(T, 23)` float32,
`dem_height_edges_m` `(22,)` float64 and `dem_nodata_share` `(T,)` float32,
each the triangle entries only, or None. A wrong type, entry count, or
component count other than edges + 1 is a `ValueError` naming the array.

### 3.7 The option, refusals and exclusions

`--dem-heights`, a flag, off by default. Help: "Write, per triangle, the
share of its area at each height about its own plane, in 1 m classes from
−10 to +10 m with open ends (dem_height_shares), and the share over NoData
(dem_nodata_share), from the DEM's own pixels." Works with `--stride`,
`--tolerance`, `--triangles`, `--domain`, `--out-crs` and cut runs, in
`.vtk`, `.ply` and `.gpkg`; mirrored in the MCP tool's arguments.

Not built: slope and aspect distributions, hollow and knoll shares (they
need a window past the triangle), moments, a tolerance-relative class
width, plan area (three vertices away), a histogram of the resampled grid.

## 4. Invariants

1. **Off is today:** without `--dem-heights`, `.vtk`, `.ply` and `.gpkg`
   bytes equal master `98fe80c3`'s for the same command.
2. **Every row sums to 1:** Σ shares + nodata = 1 within 1e-5 (the bound
   of 3.3 is 1.4e-6), for every triangle.
3. **Each pixel's area once:** over all triangles of a run (all pieces of a
   cut run), the weights a pixel gives sum to its area inside the written
   mesh, to rounding (1e-9 relative in double; 1e-6 through the float32 shares).
4. **Source, not resample:** on the reprojected path no share depends on
   the resampled grid's values.
5. **Thread count changes nothing:** the arrays are byte-identical for any
   `--threads`.
6. **Land cover unchanged:** `class_areas` returns today's bytes after the
   refactor; the land-cover suite passes untouched.

## 5. Tests (`@tester`)

Figures are from the exact `Fraction` oracle (clip, area, centroid and
class in exact arithmetic) unless said. Shares are {class in metres:
share}; unlisted classes and nodata are 0.

1. **11 × 11 nodes, 10 m, CLI.** `micro_tiff`, z = 2 m on columns 1 and 2
   (x 10 and 20 m), 0 elsewhere, `--stride 10 --dem-heights`, `.vtk` and
   `.ply`: 2 triangles, plane 0. Whichever diagonal the mesh takes, the
   triangle holding both west corners has {0: 33/50, 2: 17/50} and the other
   {0: 47/50, 2: 3/50}. The test reads the diagonal from the file. (Cut
   pixels on every side: the half columns at x = 0 and the diagonal.)
2. **NoData.** Test 1 with node (row 5, col 3) NoData (world (30, 50)): the
   triangle holding both west corners becomes {0: 16/25, 2: 17/50}, nodata 1/50; the
   other is unchanged; `dem_height_triangles_without_data` 0.
3. **All NoData, Python level** (`height_shares` with `TileWindows` over a
   hand-built 3 × 3 tile, 10 m, centre node NoData, and a hand-built mesh):
   a triangle inside the centre pixel (index corners (1.1, 1.1), (1.4, 1.1),
   (1.1, 1.4)) gets nodata 1.0, shares all 0.0, and `without_data` 1.
4. **Class edges and open ends, binding level.** 6 × 6 pixels, row_area 1;
   one triangle (0, 0), (6, 0), (0, 6) in index space, z 0 at its corners;
   DEM 0 except (row, col) (0, 0) 3.0, (0, 1) 10.5, (0, 2) −10.6, (1, 0)
   −10.5, (1, 1) −0.5, (2, 0) 0.5, and (5, 5) 7.0 outside. Area 18; shares
   at index 0, 1, 12, 14, 22 each 1/18 (open lower, −10, +1, +3, open
   upper), index 11 (class 0) 13/18.
5. **Piece centroid.** 3 × 3 pixels of 0; triangle (1.1, 1.1), (1.4, 1.1),
   (1.1, 1.4) with z 0, 3, 0: plane 1.0 at the centroid, {−1: 1}. (At the
   pixel centre it would be 4.0, class −4.)
6. **Plane as written.** Test 4 with z 1.0 at all corners: index 10
   (class −1) 13/18, index 0 2/18, index 11, 13, 21 each 1/18.
7. **Row areas.** 2 × 2 pixels, row 0 at 0 m, row 1 at 2 m, row_area
   (1, 3); triangles (0, 0), (2, 0), (0, 2) and (2, 0), (2, 2), (0, 2):
   {0: 1/2, 2: 1/2} and {0: 1/10, 2: 9/10}.
8. **Resolution, binding level.** z = 3 sin(2πx/400) cos(2πy/300) m on
   [0, 1000]² m, 10 m (101²) and 5 m (201²) nodes, the 32 triangles of a
   250 m grid split nw-se, corners at z(x, y): classes −5 … +5 occur; the
   distance between the two histograms, measured median 0.011, largest
   0.017 (probe, double); the test asserts at most 0.03.
9. **Each pixel's area once, property.** C++ (`tests/cpp/property`): 300
   random triangulations of a rectangle in index space with off-node fans;
   for 20 random pixels each, the pixel set to 100 m and the rest 0 with
   corner z 0: Σ over triangles of (plan area × share of class +11) equals
   the pixel's area inside the rectangle within 1e-6 relative (float32 shares).
10. **Cut run.** `pieces_fixtures` `hills()`, `SQUARE`, tolerance 1,
    `--pieces 4 --dem-heights`: 4 piece files, 1,163 triangles; a Python
    oracle (shapely clip in index space, centroid residual) over every
    triangle of every piece: shares and nodata within 1e-5.
11. **Reprojected path, CLI.** `tests/fixtures/velhas`, `--out-crs
    EPSG:31983 --tolerance 5 --dem-heights`: 969 triangles; the same
    oracle in the source's index space (corners reprojected to EPSG:4674,
    `row_areas` per row): shares within 1e-4. Area-weighted spread 1.289 m
    and class-0 share 0.213 (probe), asserted within 0.01 and 0.005; the
    resampled grid gives 1.219 m and 0.225, outside both.
12. **Off is byte-identical.** Tests 1, 10 and 11, a `.gpkg` of test 1
    without the flag, and `class_areas`'s three arrays on `TestThreads`' input
    in `tests/python/test_land_cover_overlap.py`: sha256 equal to digests
    recorded at `98fe80c3` (constants in the test; that commit named in a
    comment).
13. **Threads.** Test 8's 5 m case through the binding, and test 10's run,
    at threads 1, 2 and 7: equal bytes.
14. **Formats agree.** Test 1 as `.vtk`, `.ply` and `.gpkg`: the same
    shares and nodata shares; the edges array, comment and table give the
    same 22 edges; `.vtk` lines get zeros; `read_mesh` returns the triangle
    part.
15. **Reader refusals.** A hand-made `.vtk` with `dem_height_shares` of
    type `double`, 22 components, or one cell short, `dem_nodata_share` of
    type `double`, and an edges array of 21: `ValueError` naming the
    array; a file without them gives None.
16. **CLI.** `--dem-heights` in `rasputin mesh --help` and the MCP tool's
    arguments.

**Invariant-critical:** the overlap weights and the residual. Mutation
targets: in `height_shares.hpp` the two class-edge comparisons, the open-end
clamp, the centroid evaluation, the NoData branch and the denominator; in
`overlap.hpp` the whole-pixel test (`full_begin`, `full_end`) and the piece
area.

**Mutants** the suite must kill: `<=` for `<` at a class edge (4); classes
not centred, `floor(r)` (4); open ends not clamped (4); residual at the
pixel centre (5); a cut pixel counted whole, or a whole pixel clipped
wrongly (1, 9); the plane from the DEM, not the written z (6); row_area
ignored (7); NoData counted in class 0 (2, 3); the denominator without
NoData (2: 32/49 for 16/25); the resampled grid on the reprojected path
(11); a shared accumulator across threads (13); the flag ignored in
`_write_piece` (10); anything written when off (12); lines not
zero-padded (14); `class_areas`'s visitor reordered (12).

## 6. Degeneracies

A triangle smaller than a pixel: one or a few pieces, a histogram as any
other. All its pixels NoData: nodata 1.0, shares 0.0, counted (3.3). A
zero-area triangle never reaches the file (the engine writes none) and is
refused by the binding. A corner past the DEM window: `OutsideWindow`, a
defect upstream, never a user's input. A reprojected vertex with no image
(`inf`): the chunk window's `ValueError` (the mesh reaches outside the DEM),
before the binding, which would refuse it as not finite.

## 7. Cost and memory at 637,000 km², 30 m

7.1e8 DEM nodes, 23 times the São Francisco piece of section 2.

- **Time (estimate):** the walk scales with pixels and cut pixels; 0.35-0.46
  s on the piece at 10 threads gives about 10 s for the basin, plus 1.2 s
  to reproject its vertices, spread over its pieces, against the
  piece's 12 s mesh run times 23. Prediction for the fixture command at
  `--tolerance 10`: +0.05 s for the walk (1.2 times the error summary's
  0.037 s) and +0.03 s to encode 6.7 MB, on a 0.23 s run.
- **Memory:** 96 bytes per triangle while held (92 of shares, 4 of nodata
  share), about twice that while one file is encoded, plus 72 bytes per
  triangle of a chunk (corners and z, 2^18 triangles, 19 MB) while
  measured. On the fixture's land at `--tolerance 10` about 0.02 triangles
  per node: 2 bytes per node against the piece budget's estimate of 37
  (`decompose.bytes_per_node(10)`). On the reprojected path, the source
  tile held through the final check: 4 bytes per source node (float32)
  beside the check-point store's 16 (`check_points.hpp`, `Packed`). The
  piece budget's estimate is not changed; `@perf` reports the peak with
  and without.
- **Files:** 96 bytes per cell in binary `.vtk`, 97 in `.ply`.

## 8. Acceptance (timings reported, not gated)

On mains power, median of 3, against `98fe80c3`: the fixture command
(`--tolerance 10`, `.vtk`) with and without `--dem-heights`, the same cut by
`--pieces 4`, and the Velhas fixture with `--out-crs EPSG:31983`; time and
peak resident memory each.

## 9. Size

Estimated +430 counted lines: `height_shares.hpp` +90, `overlap.hpp`
refactor +15, binding and `_core.pyi` +45, `dem_heights.py` +55,
`land_cover_fractions.py` row-area helper +8, `cli.py` +65 (flag,
`_refined`, stride, piece and single writes, record), `dem_input.py` +8,
`refine_job.py` and `basin_run.py` +10, `io/vtk_legacy.py` +30, `io/ply.py`
+25, `io/geopackage_write.py` +35, `io/vtk_read.py` +40, MCP arguments +5.
One PR.

## 10. Questions for Ola

1. Classes: 1 m wide, centred on whole metres from −10 to +10 m, two open
   end classes (23). Default: **this**; 0.5 m classes from −5 to +5 m was
   the other proposal, with more noise and a narrower range.
