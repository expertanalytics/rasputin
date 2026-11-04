# Increment 137: each triangle's DEM height histogram in the mesh file

Status: **designed, revised twice** (`@architect`, 2026-11-04, on public
master `98fe80c3`; branch of worktree `dem-stats`). The design review of
`6f0e20df` returned CHANGES REQUESTED: its 50-node threshold left 72 % of
Bygdin's triangles without a histogram. The first revision weights each DEM
pixel by its exact overlap with the triangle (the main session's ruling
between the two reviewers' fixes) and reads a reprojected DEM's own pixels.
The second follows Ola's ruling on the classes (2026-11-04: "The range must
depend on the DEM vtol", vtol being `--tolerance`): the classes are a tenth
of the tolerance wide (3.1), and so a pixel's residual is taken at its
centre, the point the tolerance bounds (3.2, section 2). Full lane: a
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
  plane (Ola, 2026-11-04: "Cell stats: means simply histogramming"): 23
  classes, each a tenth of the run's tolerance wide, centred on 0, ±s, …,
  ±10 s with s = tolerance / 10, so the inner classes span the tolerance
  and half a class past it; the two end classes are open. At `--tolerance
  10` that is 1 m classes from −10.5 to +10.5 m; at `--tolerance 1`,
  0.1 m classes. Plus the share of its area over NoData. Each DEM pixel
  counts with its exact area inside the triangle, so every triangle gets a
  histogram, however small, and the 24 numbers sum to 1. A reprojected or
  geographic DEM is read on its own pixels, not on the resampled grid. The
  file states its class edges. Without the option every file is byte for
  byte today's.

```
 mesh as written (Trimmed: x, y, written z, in the mesh's CRS)
 DEM windows (SourceWindows: the tile, or the source DEM on the reprojected path)
        |
        v  [Python, dem_heights.height_shares]
 corners -> the DEM's index space (reprojected first when the CRSs differ),
 per chunk: the window under its corners, the area of one pixel per row;
 class width s = tolerance / 10
        |
        v  [C++, _core.dem_height_shares]
 per triangle: clip to each pixel -> piece area
   -> residual z_pixel - plane(pixel centre) -> class of width s; NoData -> its own share
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
elevation bands, this uses bands of height above the plane, as wide as a
tenth of the mesh's tolerance, because the plane already carries the
triangle's mean height and tilt and the tolerance bounds what is left. The weights
are the exact overlap of each pixel with the triangle, the construction of
conservative remapping and of `exactextract`'s coverage fractions, by the
clip of Sutherland and Hodgman 1974 (*Reentrant polygon clipping*, CACM
17(1): 32-42), all cited and compared in `26-land-cover.md` section 1 for
`overlap.hpp`, which this reuses. **Departures:** (a) a pixel's height
counts with its exact piece area, but (b) its residual is taken at the
pixel centre, the DEM node, against the triangle's plane, for every piece
of the pixel: that is the point the tolerance bounds, so like is compared
with like (the pixel's value and the plane's mean over the pixel's square).
For a cut pixel whose centre lies outside the triangle the plane is
extrapolated by at most 0.71 pixel; section 2 measures how much area that
puts past the tolerance. The piece centroid, the first revision's choice,
adds a staircase term of up to slope × 0.71 pixel that the tolerance does
not scale, and fills the open classes at small tolerances (section 2);
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
exact `Fraction` oracle for the test figures. The first revision's figures
(Bygdin, cost, resampled grid) were taken with 1 m classes and the piece
centroid; the second revision's (classes against the tolerance, resolution,
Velhas at 3.1's classes) are marked so. Probes of the second revision:
`rasputin_scratch/arch137-range`.

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

**Classes against the tolerance** (second revision). The fixture of the
first table, land triangles, area-weighted; "open" is the area in the two
open classes, "used" the classes holding at least 1 % of the area:

| `--tolerance` | median triangle, pixels | s = tol/10, centre: open / used | s = tol/10, centroid: open / used | 1 m, centre: open / used |
|---|---|---|---|---|
| 1 | 1.0 | 0.72 % / 17 | 15.7 % / 21 | 0.00 % / 3 |
| 5 | 3.0 | 0.12 % / 15 | 0.83 % / 15 | 0.01 % / 9 |
| 10 | 6.5 | 0.04 % / 14 | 0.13 % / 14 | 0.04 % / 14 |
| 20 | 25.5 | 0.01 % / 14 | 0.01 % / 14 | 7.43 % / 23 |

Fixed 1 m classes put 83 % of the area in class 0 at `--tolerance 1` and
7.4 % in the open classes at 20; classes a tenth of the tolerance wide use
14-17 classes at every tolerance. With the residual at the piece centroid,
15.7 % of the area at `--tolerance 1` lies past ±1.05 m although every node
is within 1 m: in triangles of under one pixel the median triangle has 75 %
of its area there, the slope across the part of a pixel between its centre
and the piece. At the pixel centre 0.72 %. On Velhas at `--tolerance 5`
(s = 0.5 m) the open classes hold 0.19 % (centroid) and 0.07 % (centre).
**Centroid or pixel centre** in the first revision, Bygdin, 1 m classes,
the distance between the two histograms (total variation: half the sum of
absolute differences; 0 equal, 1 disjoint), median by plan area in
pixels: under 1, 1.00; 1-4, 0.70; 4-16, 0.29; 16-50, 0.10; 50 and over,
0.03; the share outside class 0 moved from 78.8 % to 78.1 %.

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
triangles of a 250 m grid, 1 m classes: overlap histograms from a 10 m and
a 5 m DEM differ by 0.012 (median) and 0.018 (largest) with the residual at
the pixel centre (second revision; 0.011 and 0.017 at the centroid).
Velhas at `--tolerance 5` with 3.1's classes (s = 0.5 m, pixel centre),
969 triangles, each triangle's own spread (the standard deviation of its
class centres k·s, open classes at ±11 s) weighted by its plan area:
1.212 m from the source's pixels, 1.141 m from the resampled grid; class-0
share 0.118 and 0.148 (second revision; the first table's 1.289 / 1.219 m
and 0.213 / 0.225 are 1 m classes at the centroid, reproduced by the same
probe). Node-count histograms on the
fixture, 10 m against its 20 m block mean, differed by a median of 0.04-0.09
at 49-99 nodes and 1.5-2.0 times that with 0.5 m classes (commit
`6f0e20df`, section 2); this is the class-width argument of 3.1.

## 3. Design

### 3.1 The classes

The run's tolerance τ (metres) sets the class width s = τ / 10 (one
quotient in double). Class k, for k = −11 … +11, holds the residual r
(metres) with e(k − 1) ≤ r < e(k), where e(j) = (j + ½) · s is one product
in double; class −11 is r < e(−11) = −10.5 s, class +11 is r ≥ e(10) =
10.5 s. Index k + 11 in every array (0 … 22). C++ and Python form the 22
edges e(−11) … e(10) by the same two operations, so the class of r is
exact against the edges the file states. At τ = 10, s = 1 and the edges
are −10.5 … +10.5 exactly. Why:

- **Centred on 0.** Flat ground under a flat triangle (water, the sea of
  section 2) gives residuals of exactly 0 or rounding noise about 0; a class
  edge at 0 would split that noise between two classes. Centred, it lies in
  class 0.
- **The range follows the tolerance** (Ola's ruling). Every DEM node inside
  a tolerance mesh is within τ of it, and the largest is τ to 4 digits
  (section 2's first table), so the residual at a pixel centre inside the
  triangle is within τ. The inner classes reach half a class past it
  (±1.05 τ): the classes k·s ± s/2 are centred for every k, as class 0 is,
  so the nodes the refinement leaves at τ lie mid-class ±10, not on an
  edge where the last bit of z would choose between an inner and an open
  class. No wider margin: past it lie only cut pixels whose centre is
  outside the triangle, 0.72 % of the land area at τ = 1 and 0.04 % at 10
  (section 2), which the open classes count.
- **A fixed count, so the width scales.** The refinement stops once the
  residual is within τ, so the residual's spread scales with τ: classes
  a tenth of τ wide use 14-17 classes at τ = 1 … 20, where fixed 1 m
  classes use 3 at τ = 1 (83 % in class 0) and spill 7.4 % into the open
  classes at 20 (section 2). The count adds no parameter: it is the
  reviewed design's 23 classes, and at τ = 10, where Bygdin was measured,
  the classes are the reviewed design's 1 m classes byte for byte. A fixed
  count also keeps every array (T, 23) at every tolerance.
- **Per run, stated in the file.** Meshes at different tolerances have
  different edges, and a class index means k · s only with the file's s.
  The edges are written into every file (3.6); a reader reads them, never
  assumes them, and a reader comparing two meshes re-bins by their edges.
  The open classes have no stated width; a reader that needs one takes s,
  which the edges give.

**τ for each kind of run.** `--tolerance T`: τ = T. `--triangles N`: the
tolerance the search chose and the run record states (`tolerance_m`; a
floor answer is stated at the floor's largest error, which bounds every
node). `--tolerance-near` and `--tolerance-slope` lower the tolerance
locally and never raise it (the CLI refuses an N above `--tolerance`), so
`--tolerance` bounds every node and τ is it: one set of classes per run,
narrower residuals near lines and on steep ground. A cut run refines every
piece and seam at one `--tolerance`, so every piece file states the same
edges. Refused, before any work: `--dem-heights` without `--tolerance` or
`--triangles` ("needs --tolerance or --triangles: the height classes are a
tenth of the tolerance"), since a `--stride`-only mesh bounds nothing
(195 m at `--stride 30`, section 2); and `--tolerance 0`. A `--triangles`
run whose chosen tolerance is 0 (every node on the mesh) exits non-zero
after the search, naming `--tolerance`. No run has a tolerance that varies
upward, so none needs a second set of classes.

Constant scale: the count (10 per side) is scale-free; the edges are
metres of height through τ, checked at τ = 1, 5, 10 and 20 on the 7908
fixture and 5 on Velhas (section 2).

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
r = z_pixel − plane(c + ½, r + ½), at the pixel centre for every piece,
whole or cut, with plane the affine function through the triangle's three
corners in index space with their **written** z (so a piece's seam heights
are the written ones). In the DEM's own CRS that is the triangle's plane;
on the reprojected path it is the plane by barycentric weights (departure
c), and the centre is a source node the final check held within τ. A
centre inside the triangle is a node the refinement measured against this
plane, so |r| ≤ τ there; a centre outside it (a cut pixel) extrapolates
the plane by at most 0.71 pixel (departure b). A pixel the view's NoData
names, or a NaN, adds its weight to the NoData share and to no class.

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
so its bytes do not change. A piece's area is `area_local`'s; no centroid
is needed.

C++, `include/terrain/statistics/height_shares.hpp`, namespace
`terrain::statistics`:

```cpp
inline constexpr int half = 10;  // classes -half-1 ... half+1: 2*half + 3 = 23
struct HeightShares { std::vector<float> shares; /* T*23 */ std::vector<float> nodata; /* T */ };
template <class T>  // float or double, the DEM's promoted dtype
[[nodiscard]] HeightShares height_shares(std::span<const Point2> corners,  // 3 per triangle, index space
                                         std::span<const double> z,        // 3 per triangle, written
                                         const raster::RasterView<T>& dem, std::span<const double> row_area,
                                         double step_m,                    // s = tolerance / 10, from Python
                                         unsigned threads);
```

The class of r: k = floor(r / s + ½), then lowered while r < (k − ½) · s
and raised while r ≥ (k + ½) · s, then clamped to ±11, so it agrees with
the edges e(j) = (j + ½) · s bit for bit (3.1).

Per triangle: the canonical rotation `class_areas` applies (corners and
their z together, so the vertex order changes nothing), the plane, the
pieces, a local array of K + 1 doubles, then its row. One output row per
triangle and no reduction, so the bytes do not depend on the thread count.
Refusals, as `class_areas`'s and before any work: corner and z spans not 3
per triangle, row_area not one per row, a non-finite corner or z, a
step_m not finite and above 0, a zero-area triangle
(`std::invalid_argument`), a corner past the window
(`land_cover::OutsideWindow`).

Binding (`bindings/core.cpp`, `_core.pyi`), GIL released, the
`class_areas` binding's checks:

```python
def dem_height_shares(corners: NDArray[float64],   # (T, 3, 2) index space of `dem`
                      z: NDArray[float64],         # (T, 3)
                      dem: NDArray[float32 | float64], nodata: float | None,  # (rows, cols)
                      row_area: NDArray[float64],  # (rows,)
                      *, step_m: float, threads: int = 0) -> tuple[NDArray[float32], NDArray[float32]]: ...  # (T, 23), (T,)
```

Python, `src_python/tin_engine/dem_heights.py` (no CLI, no file I/O):

```python
HALF: Final = 10; CHUNK: Final = 2**18
def step_m(tolerance_m: float) -> float                # tolerance_m / HALF; ValueError unless finite and > 0
def edges_m(tolerance_m: float) -> NDArray[float64]    # (22,): (arange(-HALF-1, HALF+1) + 0.5) * step_m
@dataclass(frozen=True, slots=True)
class HeightShares:
    shares: NDArray[float32]     # (T, 23)
    nodata: NDArray[float32]     # (T,)
    edges_m: NDArray[float64]    # (22,), what the writers state
    without_data: int            # rows with nodata == 1.0
def height_shares(windows: SourceWindows, vertices: NDArray[float64], triangles: NDArray[uint32],
                  *, tolerance_m: float, mesh_crs: str | None, threads: int,
                  chunk: int = CHUNK) -> HeightShares
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
  `tile` past the trim when `how.dem_heights`.
- **the reprojected path** (`--out-crs` another CRS, every geographic DEM):
  the **source's** windows, never the resampled grid (section 2). `DemInput`
  and `PieceDem` gain `source: SourceWindows | None`, the `TileWindows` or
  `HeldWindows` they already build for the check points, set only when the
  option is on; `_refined` takes it, `mesh_crs` the
  target CRS. Today that source tile is freed once the check-point store is
  built (as the code reads: the exhausted block generator held its last
  reference); with the option it is held through the final check.
- `--triangles`: trials run with `dem_heights=False`; the kept mesh is
  measured once, as its error summary is, with τ the chosen tolerance.
- The no-refinement path (`--stride` alone) never runs it: refused (3.1).
- Cut runs: each piece through `_refined`; `Piece` gains `heights`, and
  `_write_piece` writes them. A piece's tile, or its `HeldWindows` part
  snapped outward to source nodes, holds every pixel its triangles touch.

`How` gains `dem_heights: bool = False`; τ is `how.tolerance`, or the
chosen tolerance under `--triangles`. With the option off nothing is held
and nothing moves.

### 3.6 The files (stable names)

| name | `.vtk` | `.ply` | `.gpkg` |
|---|---|---|---|
| `dem_height_shares` | cell array, 23 components, `float` | face `list uchar float`, always 23 | table `dem_height_shares` (`triangle`, `class`, `share`), non-zero rows only |
| `dem_nodata_share` | cell array, `float` | face `float` | `triangles` column `REAL` |
| `dem_height_edges_m` | `FieldData` array, 22 `double` | comment `dem_height_edges_m` and the 22 edges, `repr` | table `dem_height_classes` (`class`, `lower_m`, `upper_m`, NULL at the open ends) |

The edges are the run's: e(−11) … e(10) of 3.1 for its τ, the same in
every piece file of a cut run. `.vtk`: every cell array covers every cell; lines get zeros. Placed after
the existing cell arrays. `.ply`: after the existing face properties, binary
and ASCII. Record entry (so in the file's dataset strings and `--stats`):
`dem_height_triangles_without_data`. No NaN is written, so ASCII readers
need no special value. The names are the contract: renamed only with a
deprecation by a later increment.

**Reader.** `MeshArrays` gains `dem_height_shares` `(T, 23)` float32,
`dem_height_edges_m` `(22,)` float64 and `dem_nodata_share` `(T,)` float32,
each the triangle entries only, or None. A wrong type, entry count, or
component count other than edges + 1, or edges not strictly increasing,
is a `ValueError` naming the array. The reader returns the file's edges; it
checks no value against a constant.

### 3.7 The option, refusals and exclusions

`--dem-heights`, a flag, off by default. Help: "Write, per triangle, the
share of its area at each height about its own plane, in 23 classes a
tenth of the tolerance wide, centred on 0, with open ends past 1.05 times
the tolerance (dem_height_shares; the edges in dem_height_edges_m), and the
share over NoData (dem_nodata_share), from the DEM's own pixels. Needs
--tolerance or --triangles." Works with `--tolerance` (and `--stride` as
its start grid), `--triangles`, `--domain`, `--out-crs` and cut runs, in
`.vtk`, `.ply` and `.gpkg`; mirrored in the MCP tool's arguments.

Not built: slope and aspect distributions, hollow and knoll shares (they
need a window past the triangle), moments, fixed classes in metres,
plan area (three vertices away), a histogram of the resampled grid.

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
7. **The classes follow the tolerance:** every file states the edges
   (j + ½) · τ / 10, j = −11 … 10, of its run's τ (3.1), the same in every
   piece of a cut run.
8. **Inside the tolerance, inside the range:** a piece whose pixel centre
   lies inside its triangle is in an inner class; only pieces of pixels
   centred outside the triangle reach the open classes.

## 5. Tests (`@tester`)

Figures are from the exact `Fraction` oracle (clip, area and class in
exact arithmetic, the residual at the pixel centre) unless said. Shares
are {class centre in metres: share}; unlisted classes and nodata are 0.
Binding-level tests pass `step_m=1.0` (τ = 10) unless said, so their
classes are 1 m. Second revision: tests 1, 2 and 14 gain `--tolerance 10`
(measured: still 2 triangles, the same two as `--stride 10` alone, largest
difference 2 m, with and without test 2's NoData node), so their figures
stand; tests 5, 8, 10, 11 and 12 are recomputed; 17 and 18 are new.

1. **11 × 11 nodes, 10 m, CLI.** `micro_tiff`, z = 2 m on columns 1 and 2
   (x 10 and 20 m), 0 elsewhere, `--stride 10 --tolerance 10 --dem-heights`,
   `.vtk` and `.ply`: 2 triangles, plane 0, τ = 10, s = 1 m. Whichever diagonal the mesh takes, the
   triangle holding both west corners has {0: 33/50, 2: 17/50} and the other
   {0: 47/50, 2: 3/50}. The test reads the diagonal from the file. (Cut
   pixels on every side: the half columns at x = 0 and the diagonal.)
2. **NoData.** Test 1 with node (row 5, col 3) NoData (world (30, 50)): the
   triangle holding both west corners becomes {0: 16/25, 2: 17/50}, nodata 1/50; the
   other is unchanged; `dem_height_triangles_without_data` 0.
3. **All NoData, Python level** (`height_shares` with `TileWindows` over a
   hand-built 3 × 3 tile, 10 m, centre node NoData, a hand-built mesh and
   `tolerance_m=10`):
   a triangle inside the centre pixel (index corners (1.1, 1.1), (1.4, 1.1),
   (1.1, 1.4)) gets nodata 1.0, shares all 0.0, and `without_data` 1.
4. **Class edges and open ends, binding level.** 6 × 6 pixels, row_area 1;
   one triangle (0, 0), (6, 0), (0, 6) in index space, z 0 at its corners;
   DEM 0 except (row, col) (0, 0) 3.0, (0, 1) 10.5, (0, 2) −10.6, (1, 0)
   −10.5, (1, 1) −0.5, (2, 0) 0.5, and (5, 5) 7.0 outside. Area 18; shares
   at index 0, 1, 12, 14, 22 each 1/18 (open lower, −10, +1, +3, open
   upper), index 11 (class 0) 13/18.
5. **Pixel centre, outside the triangle.** 3 × 3 pixels of 0; triangle
   (1.1, 1.1), (1.4, 1.1), (1.1, 1.4) with z 0, 3, 0: the plane is 4.0 at
   the pixel centre (1.5, 1.5), {−4: 1}. (At the piece centroid it would be
   1.0, class −1.)
6. **Plane as written.** Test 4 with z 1.0 at all corners: index 10
   (class −1) 13/18, index 0 2/18, index 11, 13, 21 each 1/18.
7. **Row areas.** 2 × 2 pixels, row 0 at 0 m, row 1 at 2 m, row_area
   (1, 3); triangles (0, 0), (2, 0), (0, 2) and (2, 0), (2, 2), (0, 2):
   {0: 1/2, 2: 1/2} and {0: 1/10, 2: 9/10}.
8. **Resolution, binding level.** z = 3 sin(2πx/400) cos(2πy/300) m on
   [0, 1000]² m, 10 m (101²) and 5 m (201²) nodes, the 32 triangles of a
   250 m grid split nw-se, corners at z(x, y): classes −5 … +5 occur; the
   distance between the two histograms, measured median 0.012, largest
   0.018 (probe, double, pixel centre); the test asserts at most 0.03.
9. **Each pixel's area once, property.** C++ (`tests/cpp/property`): 300
   random triangulations of a rectangle in index space with off-node fans;
   for 20 random pixels each, the pixel set to 100 m and the rest 0 with
   corner z 0: Σ over triangles of (plan area × share of class +11) equals
   the pixel's area inside the rectangle within 1e-6 relative (float32 shares).
10. **Cut run.** `pieces_fixtures` `hills()`, `SQUARE`, tolerance 1,
    `--pieces 4 --dem-heights`: 4 piece files, 1,163 triangles (the mesh
    is the first revision's), τ = 1, s = 0.1 m; every piece file states the
    edges (j + ½) · 0.1; a Python oracle (shapely clip in index space,
    residual at the pixel centre, classes of 3.1 at s = 0.1) over every
    triangle of every piece: shares and nodata within 1e-5; and invariant
    8: each triangle's open-class share is at most the oracle's share of
    its pieces whose pixel centre lies outside it.
11. **Reprojected path, CLI.** `tests/fixtures/velhas`, `--out-crs
    EPSG:31983 --tolerance 5 --dem-heights`: 969 triangles, s = 0.5 m; the
    same oracle in the source's index space (corners reprojected to
    EPSG:4674, `row_areas` per row): shares within 1e-4. Spread is each
    triangle's own (the standard deviation of its class centres k · s, the
    open classes at ±11 s), weighted by its plan area in EPSG:31983; with
    the class-0 share weighted the same way: 1.212 m and 0.118 (probe,
    section 2), asserted within 0.01 and 0.005; the resampled grid gives
    1.141 m and 0.148, outside both.
12. **Off is byte-identical.** Tests 1 (with its `--tolerance 10`), 10 and
    11, a `.gpkg` of test 1 without the flag, and `class_areas`'s three
    arrays on `TestThreads`' input in
    `tests/python/test_land_cover_overlap.py`: sha256 equal to digests
    recorded at `98fe80c3` (constants in the test; that commit named in a
    comment). The `.gpkg` digest skips the run-record rows the host fills
    (`command`, `threads`, `rasputin_version`). File digests are keyed by
    (`sys.platform`, `platform.machine()`) and skip on a platform without
    a recorded one, as `R11_SHA256` in
    `tests/python/test_cli_catchment_geographic.py` does; `class_areas`'s
    arrays are compared on every platform.
13. **Threads.** Test 8's 5 m case through the binding, and test 10's run,
    at threads 1, 2 and 7: equal bytes.
14. **Formats agree.** Test 1 as `.vtk`, `.ply` and `.gpkg`: the same
    shares and nodata shares; the edges array, comment and table give the
    same 22 edges, −10.5 … +10.5 (τ = 10; the run's largest difference is
    2 m, so edges from it would be ±2.1); `.vtk` lines get zeros;
    `read_mesh` returns the triangle part.
15. **Reader refusals.** A hand-made `.vtk` with `dem_height_shares` of
    type `double`, 22 components, or one cell short, `dem_nodata_share` of
    type `double`, an edges array of 21, and edges not increasing:
    `ValueError` naming the array; a file without them gives None.
16. **CLI.** `--dem-heights` in `rasputin mesh --help` and the MCP tool's
    arguments.
17. **Refusals.** Test 1's DEM with `--stride 10 --dem-heights` (no
    `--tolerance`) and with `--tolerance 0 --dem-heights`: non-zero exit
    before any file, the message naming `--tolerance`; the binding with
    `step_m` 0, −1 and NaN: `ValueError`.
18. **The range follows the tolerance, CLI.** Test 1 with `--tolerance 2`
    instead (measured: still 2 triangles, largest difference 2 m): s =
    0.2 m, edges (j + ½) · 0.2; the triangle holding both west corners
    {0: 33/50, 2.0: 17/50} at index 21 (class +10, the outermost inner
    class: the 2 m nodes at exactly τ lie mid-class, not open), the other
    {0: 47/50, 2.0: 3/50}.

**Invariant-critical:** the overlap weights and the residual. Mutation
targets: in `height_shares.hpp` the two class-edge comparisons, the
division by s, the open-end clamp, the pixel-centre evaluation, the NoData
branch and the denominator; in `dem_heights.py` `step_m` and `edges_m`; in
`overlap.hpp` the whole-pixel test (`full_begin`, `full_end`) and the piece
area.

**Mutants** the suite must kill: `<=` for `<` at a class edge (4); classes
not centred, `floor(r)` (4); open ends not clamped (4); residual at the
piece centroid (5, 10); s fixed at 1 m (10, 11, 18); no half-class margin,
the inner range ±τ (18: the 2 m nodes go open); the range from the largest
difference, not τ (14); edges not stated per run (10, 18); a cut pixel counted whole, or a whole pixel clipped
wrongly (1, 9); the plane from the DEM, not the written z (6); row_area
ignored (7); NoData counted in class 0 (2, 3); the denominator without
NoData (2: 32/49 for 16/25); the resampled grid on the reprojected path
(11); a shared accumulator across threads (13); the flag ignored in
`_write_piece` (10); anything written when off (12); lines not
zero-padded (14); `class_areas`'s visitor reordered (12).

## 6. Degeneracies

A triangle smaller than a pixel: one or a few pieces, a histogram as any
other; its pixel centres usually lie outside it, so its residuals are the
plane extrapolated (departure b) and may fall in the open classes. A
tolerance of 0 has classes of no width: refused (3.1). All its pixels NoData: nodata 1.0, shares 0.0, counted (3.3). A
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

Estimated +430 counted lines (second revision: the refusals and τ's
threading +10, the centroid and the stride path −10): `height_shares.hpp` +90, `overlap.hpp`
refactor +15, binding and `_core.pyi` +45, `dem_heights.py` +55,
`land_cover_fractions.py` row-area helper +8, `cli.py` +65 (flag,
`_refined`, refusals, piece and single writes, record), `dem_input.py` +8,
`refine_job.py` and `basin_run.py` +10, `io/vtk_legacy.py` +30, `io/ply.py`
+25, `io/geopackage_write.py` +35, `io/vtk_read.py` +40, MCP arguments +5.
One PR.

## 10. Questions for Ola

1. Classes a tenth of the tolerance wide, 23 at every tolerance (1 m at
   `--tolerance 10`, 0.1 m at 1), the inner ones reaching half a class past
   the tolerance. Default: **this**; fixed 1 m classes over ±tolerance was
   the other way, which leaves 3 classes at tolerance 1 and 41 at 20.
