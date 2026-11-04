# Increment 137: each triangle's DEM height histogram in the mesh file

Status: **designed, revised three times** (`@architect`, 2026-11-04, on
public master `98fe80c3`; branch of worktree `dem-stats`). The design
review of `6f0e20df` returned CHANGES REQUESTED: its 50-node threshold left
72 % of Bygdin's triangles without a histogram. The first revision weights
each DEM pixel by its exact overlap with the triangle and reads a
reprojected DEM's own pixels. The second tied the classes to the tolerance
(Ola, 2026-11-04: "The range must depend on the DEM vtol") by a fixed count
of 23, a tenth of the tolerance wide. Ola rejected that ("You can't be
serious that on a 1m vtol mesh, you propose 23 buckets": 0.1 m classes are
finer than the DEM's own accuracy). The third, this one, rests on four
choices, each with its owner:

- **Ola:** the range follows the tolerance, the count shrinks with the
  triangle, and it converges ("We need less statistical buckets the smaller
  the triangle is. So this should also converge").
- **Main session:** a width equal to the DEM's published vertical accuracy,
  the count following from the tolerance; Ola replied "Exactly." to it.
- **Main session:** relative, not absolute, accuracy (a height difference
  cancels systematic offsets).
- **This design:** GLO-30's gentle-slope figure and the 1 m default (3.1).

Each triangle carries only the classes its own DEM nodes reach, so the
count shrinks to one as the mesh refines (3.2, invariant 9). Full lane: a per-triangle number a downstream
reader takes as input, wrong without a crash if the overlap or the plane is
wrong. New C++ lives in `include/terrain/statistics/`;
`include/terrain/land_cover/overlap.hpp` is refactored so both uses share
one clip; nothing under `include/terrain/refinement/` or
`include/terrain/mesh/` changes and the refine path is unchanged with the
option off, so no `@perf` acceptance run is required
(`docs/increments/README.md`, *Acceptance*); section 8 asks for timings,
reported and not gated (`perf.md` section 4, an option off by default).

## 0. The plain picture

- **Today** a mesh file says nothing about the DEM under a triangle except
  through its three corner heights.
- **After**, `rasputin mesh --dem-heights` writes, for each triangle, the
  share of its area at each height above or below the triangle's own
  plane (Ola, 2026-11-04: "Cell stats: means simply histogramming"), in
  classes of one fixed height width w: the DEM's published relative
  (point-to-point) vertical accuracy (1.216 m for ANADEM and GLO-30), else
  1 m (DTM10's figure). Classes are
  centred on 0 (…, −w, 0, +w, …). A triangle carries only the run of
  classes from the lowest to the highest that its own DEM nodes fall in:
  one class when its ground is within ±w/2 of its plane. Finer meshes
  have smaller residuals and so fewer classes; once the tolerance is below
  w/2 every triangle has one. The tolerance caps the count: at most
  2n + 1 classes, n the least whole number with (n + ½) w > tolerance
  (3 at tolerance 1 m and w = 1 m, 21 at 10 m). Plus the share over
  NoData. Each pixel counts with its exact area inside the triangle; the
  shares and the NoData share sum to 1. The file states the class edges
  and, per triangle, where its classes start and how many there are.
  Without the option every file is byte for byte today's.

```
 mesh as written (Trimmed: x, y, written z, in the mesh's CRS)
 DEM windows (SourceWindows: the tile, or the source DEM on the reprojected path)
 class width w (the DEM source's relative_vertical_accuracy_m, else 1 m) and tolerance
        |
        v  [Python, dem_heights.height_shares]
 n = half_count(tolerance, w); corners -> the DEM's index space (reprojected
 first when the CRSs differ); per chunk: the window under its corners, row areas
        |
        v  [C++, _core.dem_height_shares(…, class_m=w, half=n)]
 per triangle: clip to each pixel -> piece area
   -> residual z_pixel - plane(pixel centre) -> class k, |k| <= n
   -> the span [lo, hi] of the classes of the nodes it holds; other pieces fold to the nearer end
   -> its shares (hi - lo + 1 floats), first class, NoData share
        |
        v
 sparse table: offsets (T+1), first (T), shares (S), nodata (T), edges (2n+2)
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
Appl. Climatol. 52: 95-118, doi 10.1007/BF00865510); they use a fixed set
of absolute elevation bands per cell. This uses bands of height above the
plane, one fixed width for the run, and per triangle only the bands its
relief reaches. The weights are the exact overlap of each pixel with the
triangle, the construction of conservative remapping and of
`exactextract`'s coverage fractions, by the clip of Sutherland and Hodgman
1974 (*Reentrant polygon clipping*, CACM 17(1): 32-42), all cited and
compared in `26-land-cover.md` section 1 for `overlap.hpp`, which this
reuses. The sparse table per triangle is the compressed-row layout the
land-cover fractions already use (`ClassFractions`, `io/models.py`).

*The width.* A class narrower than the DEM's error resolves noise, so the
width is the **relative** (point-to-point) error the DEM's distributor
publishes, as one standard deviation in metres (main session's ruling,
2026-11-04): relief inside a triangle is a difference of heights, and an
absolute figure carries systematic offsets that cancel in a difference.
A relative LE90 is already the 90th percentile of a difference of two
points, so its σ is LE90 / 1.645, with no √2:

- **Kartverket DTM10**: "± 1 meter standardavvik i høyde i områder
  etablert i NDH" (± 1 m standard deviation in areas established in the
  national detailed height model), "± 2 til 6 meter standardavvik i høyde
  i områder som ikke er etablert i NDH ennå" (data.norge.no, *DTM 10
  Terrengmodell UTM33 2021*, read 2026-11-04). That is an absolute height
  figure; Kartverket publishes no relative one for DTM10 (searched
  2026-11-04: data.norge.no, Geonorge, the FKB-Laser specifications). An
  absolute σ bounds the relative one from above, so **1 m** (the current
  DTM10 is made from hoydedata.no's laser model) is kept, as an upper
  bound, and the entry says so.
- **Copernicus GLO-30**: relative vertical accuracy "< 2 m LE90 for
  slopes ≤ 20 %, < 4 m for slopes > 20 %", point to point within 1° × 1°
  (Product Handbook v5.0, 29.11.2022, sections 2.1-2.2; quoted in
  `docs/research/anadem-accuracy.md` section 2), so 2 / 1.645 =
  **1.216 m** (the conversion is ours, for normal errors). The gentle-slope
  figure, not the worst: one number per product, and a width set by steep
  ground would hide real relief on the flats, which are most of the São
  Francisco basin. On ground steeper than 20 % GLO-30's relative σ is
  about 2.4 m (4 / 1.645), twice the class width, so there a triangle's
  outer classes partly resolve noise.
- **ANADEM**: publishes only absolute figures (against ICESat-2: bias
  1.50 m, RMSE 6.99 m, standard deviation 6.79 m; Laipelt et al. 2024,
  *Remote Sens.* 16(13): 2321, Table 3), so it takes the relative figure
  of the GLO-30 it is built from, **1.216 m**. Its 3 × 3 median filter
  and canopy correction may change that; nobody has measured it
  (`docs/research/anadem-accuracy.md` section 2).
- **Anything else**: 1 m, the DTM10 figure, the finest product the
  pipeline meshes.

**Departures:** (a) a pixel's height counts with its exact piece area, but
(b) its residual is taken at the pixel centre, the DEM node, against the
triangle's plane, for every piece of the pixel: that is the point the
tolerance bounds. (c) A triangle's classes span only the classes of the
nodes it **holds** (centres inside it or on its edges); a piece of a pixel
centred outside it, whose residual extrapolates the plane by up to 0.71
pixel, counts in the nearer end class of that span. Without (c) such pieces
add a class to 29-70 % of the fixture's land triangles and stop the count
from converging (section 2); (c) moves 0.1-0.7 % of the land area by one
class or more. (d) On the reprojected path the triangle is carried into the
source's index space by its corners, with straight edges, and its plane by
barycentric weights, which assumes the map is affine across one triangle
(section 2 measures the error). No novelty is claimed; none searched.

*Legacy.* `git grep -n -i -E "hist|hypso|rough" legacy-archive -- legacy`
finds no per-triangle DEM statistic; `legacy/rasputin/writer.py` wrote named
per-triangle arrays through meshio. Carried across: a per-triangle quantity
is a named cell array of the mesh file.

*Ideas file.* Takes up the "carry each triangle's sub-grid distribution"
candidate of "A terrain mesh that follows the physics, not only the height"
(heights only; slope and aspect not). Overlaps "A real storage format for
meshes" (a sparse per-cell table for its reader to carry) and "Mean, RMS
and percentiles of the vertical error in --stats" (the histogram gives a
triangle's mean signed misfit to within half a class); takes up neither.

## 2. Measured

Scratch probes on master `98fe80c3`, M1 Max on mains power, C++ probes
built on `overlap.hpp`'s clip with `-O3`, and an exact `Fraction` oracle
for the test figures. The first revision's figures (Bygdin, cost, resampled
grid) were taken with 1 m classes and the piece centroid, and stand for the
weights and the cost. The third revision's (marked so) are in
`rasputin_scratch/arch137-classes` (`span.cpp`, `span2.cpp`, `an2.py`,
`an3.py`, `t8.py`, `mk.py`); they weight each piece by its index area
(equal row areas), which on Velhas changes no figure by more than 5e-4
relative.

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

Peak resident memory 40 MB and 39 MB: the inputs. On the São Francisco
piece below (1,858,576 triangles, 30.9 M nodes, 10 threads): 0.42-0.47 s
against 0.35-0.46 s, and reprojecting its 933,232 vertices 0.05 s.

**Classes per triangle** (third revision). The 7908 fixture, w = 1 m (its
DEM is a path: no source entry), classes by the span of held nodes
(departure c). "One" is the share of triangles with one class; "mean" the
mean count, per triangle and weighted by plan area; "folded" the area moved
by departure (c):

| `--tolerance` | n, bound 2n+1 | triangles, land | land: one | land: mean (area-weighted) | largest | all: one | all: mean (area-weighted) | land folded |
|---|---|---|---|---|---|---|---|---|
| 20 | 20, 41 | 50,147, 10,904 | 2.7 % | 15.07 (15.79) | 40 | 78.5 % | 4.08 (1.71) | 0.11 % |
| 10 | 10, 21 | 69,571, 30,068 | 8.9 % | 7.41 (9.12) | 21 | 60.3 % | 3.79 (1.38) | 0.23 % |
| 5 | 5, 11 | 116,389, 76,263 | 19.5 % | 3.78 (5.14) | 11 | 46.9 % | 2.83 (1.19) | 0.40 % |
| 2 | 2, 5 | 261,221, 219,199 | 41.4 % | 1.88 (2.49) | 5 | 50.6 % | 1.74 (1.06) | 0.64 % |
| 1 | 1, 3 | 464,290, 419,922 | 69.6 % | 1.33 (1.63) | 3 | 72.4 % | 1.30 (1.03) | 0.70 % |

On land both measures move monotonically to 1, as do the means over all
triangles. The one-class share over all triangles does not (78.5 → 46.9 →
72.4 %): sea triangles always have one class and are 78 % of the triangles
at 20 and 10 % at 1, a change of mix, not of any triangle. With the span
over **all** pieces instead (no departure c), the land figures are 1.6,
5.9, 14.1, 30.5, 51.3 % one class and means 19.43, 10.47, 5.53, 2.64,
1.73: pieces centred outside raise the count of 69.8, 67.1, 57.9, 41.4 and
29.0 % of land triangles. Every triangle holds a node here (its corners lie
on nodes).

`pieces_fixtures.hills(n=201)` (10 m, uniform noise 0-3 m on smooth hills),
`--stride 200`, w = 1 m: at `--tolerance` 8, 4, 2, 1, 0.4: 552, 1,546,
11,279, 42,457, 65,379 triangles; one class 0.7, 0.6, 20.6, 73.1, 100.0 %;
mean 7.783, 4.875, 2.492, 1.283, 1.000; largest 16, 9, 5, 3, 1 (bounds 17,
9, 5, 3, 1); folded 0.17-1.23 %. The span over all pieces gives 96.0 % and
1.042 at 0.4: it does not converge. Each run 0.4 s.

Velhas (`tests/fixtures/velhas`, ANADEM, `--out-crs EPSG:31983 --tolerance
5`, 969 triangles, source pixels): at w = 1 m, 11.7 % one class, mean 3.99,
largest 10 (bound 11), folded 1.42 %; at w = 1.216 m (ANADEM's, 3.1), n = 4,
13.1 % one class, mean 3.47 (area-weighted 4.40), largest 8 (bound 9),
folded 1.26 % (for comparison, at ANADEM's absolute σ 6.79 m: 73.9 %, 1.26). 22 of its triangles
hold no source node (corners are not source nodes on this path).

**Centroid or pixel centre** (second revision, fixture, classes a tenth of
the tolerance): with the residual at the piece centroid, 15.7 % of the land
area at `--tolerance 1` lies past ±1.05 m although every node is within
1 m; at the pixel centre 0.72 %. Bygdin, 1 m classes, the two histograms'
distance (half the sum of absolute differences), median by plan area in
pixels: under 1, 1.00; 1-4, 0.70; 4-16, 0.29; 16-50, 0.10; 50 and over,
0.03.

**Resampled grid or the source's pixels.** The reprojected path refines
against a bilinearly resampled grid. Histograms from it against histograms
from the source's own pixels (the nodes `final_check.stored` files), same
mesh, area-weighted over triangles (first revision, 1 m classes, centroid):

| fixture | triangles | distance | spread, resampled / source | class-0 share, resampled / source | resampled narrower in |
|---|---|---|---|---|---|
| `tests/fixtures/velhas` ANADEM, `--out-crs EPSG:31983 --tolerance 1` | 7,386 | 0.25 | 0.674 / 0.705 m | 0.604 / 0.556 | 52 % of triangles |
| same, `--tolerance 5` | 969 | 0.15 | 1.219 / 1.289 m | 0.225 / 0.213 | 70 % |
| São Francisco piece `bho2017_5k_76949_anadem_window_epsg4674.tif`, `--tolerance 5` | 1,858,576 | 0.24 | 1.347 / 1.463 m | 0.216 / 0.194 | 68 % |

(Spread: the standard deviation of the class, in metres.) The resampled
histograms are narrower, with more area in class 0. **The faithful path is
chosen** (3.5). Third revision, Velhas at `--tolerance 5`, w = 1 m, pixel
centre, span rule, each triangle's own spread (the standard deviation of
its class centres k·w) weighted by plan area: 1.164 m from the source's
pixels, 1.121 m from the resampled grid; class-0 share 0.228 and 0.262.
Departure (d) measured on Velhas: a triangle's centroid mapped to the
source's index space lies at most 1.2e-3 m (`--tolerance 1`) and 3.7e-3 m
(`--tolerance 5`) from the mean of its mapped corners.

**Resolution.** z = 3 sin(2πx/400) cos(2πy/300) m on [0, 1000]² m, 32
triangles of a 250 m grid, w = 1 m, n = 10: overlap histograms from a 10 m
and a 5 m DEM differ by 0.012 (median) and 0.018 (largest), classes −5 …
+5 occur, and every triangle has the same count from both (5-8; third
revision, span rule).

## 3. Design

### 3.1 The classes

The width w (metres) is the DEM source's `relative_vertical_accuracy_m`, else
`DEFAULT_CLASS_M` = 1.0. Class k holds the residual r (metres) with
(k − ½) · w ≤ r < (k + ½) · w; each bound is one product in double. The
run's n is the least whole number ≥ 0 with (n + ½) · w > τ (τ the
tolerance, below), computed as `floor(τ / w + ½)` raised while
(n + ½) · w ≤ τ. Every class index is clamped to −n … n; index c = k + n
(0 … 2n) in every array. The run's 2n + 2 edges are (j + ½) · w,
j = −n − 1 … n, formed by the same product in C++ and Python, so the class
of r is exact against the edges the file states. At w = 1 m: n = 1 at
τ = 1 (edges ±0.5, ±1.5), n = 10 at τ = 10 (edges −10.5 … +10.5, the
reviewed design's), n = 0 below τ = 0.5. Why:

- **The width is the DEM's, not the tolerance's** (the main session's
  proposal, to which Ola replied "Exactly."; owners in the status). A
  class narrower than the DEM's own error splits noise; the width is the
  distributor's published relative figure (section 1). **No user override**: the
  width is a property of the DEM product, stated once where the product
  is registered, so two meshes of one DEM always share their classes and
  no run tunes it. A DEM given by path, or registered without the field,
  takes 1 m, and the run record says so; a user with a DEM of known
  accuracy registers it with `relative_vertical_accuracy_m` in their sources file
  (`sources.py` already reads local entries). Python resolves the source
  and passes the core two numbers, `class_m` and `half`: no key, path or
  CRS crosses.
- **Centred on 0.** Flat ground under a flat triangle (water, the sea of
  section 2) gives residuals of exactly 0 or rounding noise about 0; centred,
  it lies in one class.
- **The count per triangle** (Ola's third ruling). Each triangle carries
  the classes lo … hi of the non-NoData nodes it holds (departure c),
  where every residual is within τ, so |k| ≤ n by the margin (n + ½) w >
  τ: the bound 2n + 1 holds by construction, without open end classes.
  The strict margin keeps a node at exactly τ inside: at τ = 2.5, w = 1,
  n = 3, not 2. A triangle holding no node (a feature vertex or a
  reprojected corner off the nodes; 22 of Velhas's 969) has one class, the
  class of its area-weighted mean residual clamped to ±n: it carries no
  node the tolerance checked, so nothing supports a spread. Inside a span
  a class may hold 0.
- **Convergence** (invariant 9). Below τ = w/2, n = 0 and every triangle
  with data has one class. Above it the count follows the triangle's own
  relief, and section 2 measures it falling toward 1 as τ falls.
- **The degenerate run does not collapse further.** At n = 0 every
  triangle still writes its one share and its NoData share: readers meet
  one layout at every tolerance, and the NoData share is still
  information. 12 bytes per triangle more than nothing.
- **Per run, stated in the file.** The edges and w are written into every
  file (3.6); a reader reads them, never assumes them.

**τ for each kind of run.** `--tolerance T`: τ = T. `--triangles N`: the
tolerance the search chose and the run record states (`tolerance_m`; a
floor answer is stated at the floor's largest error, which bounds every
node). `--tolerance-near` and `--tolerance-slope` only lower the tolerance
locally, so `--tolerance` bounds every node and τ is it. A cut run refines
every piece at one `--tolerance` against one source, so every piece file
states the same edges. Refused, before any work: `--dem-heights` without
`--tolerance` or `--triangles` ("needs --tolerance or --triangles: the
tolerance bounds the height classes"), since a `--stride`-only mesh bounds
nothing (195 m at `--stride 30`); `--dem-heights` with `--tolerance 0`;
and n > `MAX_HALF` = 127 ("--tolerance 200 at class width 1 m needs 401
height classes; at most 255"), since a `.ply` list holds at most 255. Under
`--triangles` the last two are checked after the search, exiting non-zero
and naming `--tolerance`. Without the flag all run as today.

Constant scale: `DEFAULT_CLASS_M` is metres of height, DTM10's published
absolute σ (an upper bound on its relative one), checked at τ = 1 … 20 on the 7908 fixture; `MAX_HALF` is a count (the
`.ply` `uchar` list), reached at τ ≥ 127.5 w; the held-node margin 1e-9
(3.2) is index units, checked at 5,051 columns.

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
has an area a (index units) times the row's pixel area A_row (m²: dx · dy
projected, the ellipsoidal cell area of the row for a geographic DEM,
`land_cover_fractions.row_areas`'s rule). Its residual is ρ = z_pixel −
plane(c + ½, r + ½), at the pixel centre for every piece, with plane the
affine function through the triangle's three corners in index space with
their **written** z. A pixel the view's NoData names, or a NaN, adds its
weight to the NoData share and to no class.

A node is **held** by the triangle when each of its three edge functions,
divided by its edge's length, is ≥ −1e-9 at the node (index units): inside
or on an edge, with a margin for the rounding of a corner placed on a node
(≤ 5,051 · 2⁻⁵² ≈ 1e-12 on the fixture). A node on a shared edge is held
by both triangles; none is owned. In the DEM's own CRS a held node is one
the refinement measured against this plane, so |ρ| ≤ τ there.

No pixel is owned: a pixel cut by an edge gives each triangle its own
piece, and the pieces of a pixel tile it to rounding, across piece seams
too, since each crossing is computed from the edge's endpoints in one
canonical order (`overlap.hpp`'s header).

### 3.3 The shares

Per triangle, a local array of 2n + 1 doubles takes each piece's weight
a · A_row at its class clamped to ±n; the held nodes' classes give lo and
hi. Then the weights below lo are added to lo and those above hi to hi, in
ascending class order, and the span lo … hi is the row. W = Σ a · A_row
over all pieces, NoData included. share_i = (weight in class lo + i) / W
and nodata = (weight on NoData) / W, each float32; a row's count + 1
numbers sum to 1 within (count + 1) · 2⁻²⁴ (a bound; 2.5e-6 at 41). A
triangle whose pixels are all NoData has count 0 and nodata 1.0; the run
record counts such triangles (`dem_height_triangles_without_data`), and
`--stats` shows the count. A reader wanting the heights of the valid part
divides the shares by 1 − nodata.

### 3.4 The interfaces (`@tester` calls these)

`overlap.hpp` refactor: `detail::overlap_one` becomes
`detail::for_each_piece(tri, rows, cols, visit)`, calling
`visit(r, c, const Polygon* piece)` with `nullptr` for a whole pixel;
`class_areas` passes a visitor doing today's arithmetic in today's order,
so its bytes do not change.

C++, `include/terrain/statistics/height_shares.hpp`, namespace
`terrain::statistics`:

```cpp
struct HeightShares {
    std::vector<std::uint32_t> offsets;  // T + 1: triangle t's shares are [offsets[t], offsets[t + 1])
    std::vector<std::int32_t> first;     // T: class index (0 … 2·half) of t's first share; 0 when it has none
    std::vector<float> shares;           // offsets[T]
    std::vector<float> nodata;           // T
};
template <class T>  // float or double, the DEM's promoted dtype
[[nodiscard]] HeightShares height_shares(std::span<const Point2> corners,  // 3 per triangle, index space
                                         std::span<const double> z,        // 3 per triangle, written
                                         const raster::RasterView<T>& dem, std::span<const double> row_area,
                                         double class_m, int half, unsigned threads);
```

The class of ρ: k = floor(ρ / w + ½), lowered while ρ < (k − ½) · w and
raised while ρ ≥ (k + ½) · w, then clamped to ±half. Per triangle: the
canonical rotation `class_areas` applies (corners and their z together),
the plane, the pieces, the local array, then its row. Rows go to per-block
buffers concatenated in triangle order, so the bytes do not depend on the
thread count. Refusals, before any work: corner and z spans not 3 per
triangle, row_area not one per row, a non-finite corner or z, class_m not
finite and above 0, half outside 0 … 127, a zero-area triangle
(`std::invalid_argument`), a corner past the window
(`land_cover::OutsideWindow`).

Binding (`bindings/core.cpp`, `_core.pyi`), GIL released, the
`class_areas` binding's checks:

```python
def dem_height_shares(corners: NDArray[float64],   # (T, 3, 2) index space of `dem`
                      z: NDArray[float64],         # (T, 3)
                      dem: NDArray[float32 | float64], nodata: float | None,  # (rows, cols)
                      row_area: NDArray[float64],  # (rows,)
                      *, class_m: float, half: int, threads: int = 0,
                      ) -> tuple[NDArray[uint32], NDArray[int32], NDArray[float32], NDArray[float32]]: ...
                      # offsets (T+1,), first (T,), shares (S,), nodata (T,)
```

`io/models.py` gains `HeightTable`, frozen, the `ClassFractions` pattern:
`offsets` uint32 (T + 1), `first` int32 (T), `shares` float32 (S),
`nodata` float32 (T), `edges_m` float64 (2n + 2), `class_m` float. Its
validator refuses: another dtype or rank; offsets not starting at 0,
decreasing, or ending other than at S; `first` or `nodata` not T long;
edges fewer than 2, of odd length or not increasing; a row with
`first + count > len(edges_m) − 1` or `first < 0`. Property `count`
(`diff(offsets)`).

Python, `src_python/tin_engine/dem_heights.py` (no CLI, no file I/O):

```python
DEFAULT_CLASS_M: Final = 1.0; MAX_HALF: Final = 127; CHUNK: Final = 2**18
def half_count(tolerance_m: float, class_m: float) -> int      # 3.1's n; ValueError unless both finite and > 0, or n > MAX_HALF
def edges_m(half: int, class_m: float) -> NDArray[float64]     # (2·half + 2,): (arange(-half - 1, half + 1) + 0.5) * class_m
def height_shares(windows: SourceWindows, vertices: NDArray[float64], triangles: NDArray[uint32],
                  *, tolerance_m: float, class_m: float, mesh_crs: str | None, threads: int,
                  chunk: int = CHUNK) -> tuple[HeightTable, int]  # the table and its triangles without data
```

`mesh_crs` None: the mesh is in `windows.meta`'s CRS, no reprojection;
otherwise the vertices are reprojected from it. Per chunk of `chunk`
triangles: the pixels under its corners (floor of the least to ceil of the
greatest, clipped to the raster; a window past a `HeldWindows` part is its
`ValueError`, a defect), `windows.window(...)`, the corners shifted to it,
row areas for it; chunk tables are joined by shifting each chunk's offsets.
`row_areas`'s body moves to a helper over (top edge, dx, dy, rows, CRS)
that it and `dem_heights` call.

`sources.py`: `Source` gains `relative_vertical_accuracy_m: PositiveFloat |
None` ("one standard deviation of the difference of two nearby heights, as
the distributor states it or as its LE90 / 1.645, in metres"), refused unless `values == "heights"`, as `position_error_m` is
for non-features. `data/sources.toml`: `anadem-v1` and `glo30` 1.216,
each with the handbook's quotation and, for ANADEM, that the figure is
its parent's (section 1). `cli.py` takes w from the entry
`_dem_sources` returns, else `DEFAULT_CLASS_M`. The main session adds
`relative_vertical_accuracy_m = 1.0` to Ola's local `dtm10-norway` entry
(outside the repository; equal to the default), with a comment that it is
Kartverket's absolute σ, kept as an upper bound.

### 3.5 Where it runs

On the mesh as written, after the trim, against:

- **the DEM's CRS** (no `--out-crs`, or the DEM's own): the tile the mesh
  was refined on, `TileWindows(tile)`, `mesh_crs=None`. `_refined` keeps
  `tile` past the trim when `how.dem_heights`.
- **the reprojected path** (`--out-crs` another CRS, every geographic DEM):
  the **source's** windows, never the resampled grid (section 2). `DemInput`
  and `PieceDem` gain `source: SourceWindows | None`, the `TileWindows` or
  `HeldWindows` they already build for the check points, set only when the
  option is on; `_refined` takes it, `mesh_crs` the target CRS. With the
  option the source tile is held through the final check.
- `--triangles`: trials run with `dem_heights=False`; the kept mesh is
  measured once, with τ the chosen tolerance.
- The no-refinement path (`--stride` alone) never runs it: refused (3.1).
- Cut runs: each piece through `_refined`; `Piece` gains `heights`, and
  `_write_piece` writes them. A piece's tile, or its `HeldWindows` part
  snapped outward to source nodes, holds every pixel its triangles touch.

`How` gains `dem_heights: bool = False` and `dem_class_m: float = 1.0`;
τ is `how.tolerance`, or the chosen tolerance under `--triangles`. With
the option off nothing is held and nothing moves.

### 3.6 The files (stable names)

| name | `.vtk` | `.ply` | `.gpkg` |
|---|---|---|---|
| `dem_height_shares` | dataset array, S `float` | face `list uchar float` (its length is the count) | table `dem_height_shares` (`triangle`, `class`, `share`), every class of every span, zeros inside a span kept |
| `dem_height_offsets` | dataset array, T + 1 `unsigned_int` | (the list lengths) | (the rows) |
| `dem_height_first` | dataset array, T `int` | face `int` | (the `class` column) |
| `dem_nodata_share` | cell array, `float`, 0 on lines | face `float` | `triangles` column `REAL` |
| `dem_height_edges_m` | dataset array, 2n + 2 `double` | comment `dem_height_edges_m` and the edges, `repr` | table `dem_height_classes` (`class`, `lower_m`, `upper_m`) |

**How a reader finds a triangle's classes:** triangle t (the t-th
triangle, lines not counted, as `land_cover_offsets` counts them) has
count = offsets[t + 1] − offsets[t]; its i-th share is
shares[offsets[t] + i], the share of class index first[t] + i, between
edges[first[t] + i] and edges[first[t] + i + 1]. In `.ply` the face's list
is its shares and `dem_height_first` its first class; in `.gpkg` the rows
of `dem_height_shares` for the triangle, `class` being the class index.

Why this layout: the count varies by triangle, from 1 to 2n + 1, and the
mean is far below the bound (3.79 against 21 at τ = 10 on the fixture), so
a padded (T, 2n + 1) array would hold mostly zeros and its width would
still change from run to run. The sparse table costs 4 bytes per share and
12 per triangle, and is the layout `land_cover_offsets`/`_classes`/
`_fractions` already put in all three formats, so every writer and reader
extends a pattern it has. Contiguous spans store one first class instead
of a class per share. Record entries (so in the file's dataset strings and
`--stats`): `dem_height_triangles_without_data`, `dem_height_class_m`, and
`dem_height_class_from` (`relative_vertical_accuracy_m of source <key>`, or
`default: the DEM states no vertical accuracy`). No NaN is written. The
names are the contract: renamed only with a deprecation by a later
increment.

**Reader.** `MeshArrays` gains `dem_heights: HeightTable | None`, from the
arrays above; a wrong type or length is the model's `ValueError`, naming
the array. The reader checks no value against a constant.

### 3.7 The option, refusals and exclusions

`--dem-heights`, a flag, off by default. Help: "Write, per triangle, the
share of its area at each height about its own plane, in classes as wide
as the DEM's published relative vertical accuracy (1 m when its source states
none), centred on 0; each triangle carries the classes its own DEM nodes
reach (dem_height_shares, dem_height_offsets, dem_height_first; the edges
in dem_height_edges_m), and the share over NoData (dem_nodata_share),
from the DEM's own pixels. Needs --tolerance or --triangles." Works with
`--tolerance` (and `--stride` as its start grid), `--triangles`,
`--domain`, `--out-crs` and cut runs, in `.vtk`, `.ply` and `.gpkg`;
mirrored in the MCP tool's arguments.

Not built: a class-width option, slope and aspect distributions, hollow
and knoll shares, moments, plan area, a histogram of the resampled grid.

## 4. Invariants

1. **Off is today:** without `--dem-heights`, `.vtk`, `.ply` and `.gpkg`
   bytes equal master `98fe80c3`'s for the same command.
2. **Every row sums to 1:** Σ shares + nodata = 1 within (count + 1) · 2⁻²⁴
   for every triangle, the bound of 3.3 (2.5e-6 at 41 classes, 2⁻¹⁶ ≈
   1.53e-5 at 255).
3. **Each pixel's area once:** over all triangles of a run (all pieces of a
   cut run), the weights a pixel gives sum to its area inside the written
   mesh, to rounding (1e-9 relative in double; 1e-6 through float32).
4. **Source, not resample:** on the reprojected path no share depends on
   the resampled grid's values.
5. **Thread count changes nothing:** the arrays are byte-identical for any
   `--threads`.
6. **Land cover unchanged:** `class_areas` returns today's bytes after the
   refactor; the land-cover suite passes untouched.
7. **The classes are the DEM's and the tolerance's:** every file states the
   edges (j + ½) · w, j = −n − 1 … n, w its source's width and n 3.1's for
   its τ, the same in every piece of a cut run.
8. **The span is the held nodes':** each triangle's classes run from the
   lowest to the highest class of the non-NoData nodes it holds (one class,
   the mean residual's, when it holds none; none when it has no data), so
   count ≤ 2n + 1.
9. **Convergence:** at τ < w/2 every triangle with data has exactly one
   class; on a refined sequence of meshes of one DEM the mean count falls
   toward 1 (measured, section 2).

## 5. Tests (`@tester`)

Figures are from the exact `Fraction` oracle (clip, area and class in
exact arithmetic, the residual at the pixel centre, the span rule) and
were reproduced by the probe `span2.cpp` unless said. A triangle's shares
are written {class k in metres / w: share} with its first class index and
count; classes inside the span not listed hold 0, and nodata is 0 unless
listed. Binding-level tests pass `class_m=1.0, half=10` unless said. The
CLI fixtures are paths, so w = 1 m. Third revision: tests 1, 2, 4-9, 10,
11, 14, 15, 17, 18 change; 19-22 are new.

1. **11 × 11 nodes, 10 m, CLI.** `micro_tiff`, z = 2 m on columns 1 and 2
   (x 10 and 20 m), 0 elsewhere, `--stride 10 --tolerance 10
   --dem-heights`, `.vtk` and `.ply`: 2 triangles, plane 0, n = 10.
   Whichever diagonal the mesh takes, the triangle holding both west
   corners has first 10, count 3, {0: 33/50, 2: 17/50}, and the other
   first 10, count 3, {0: 47/50, 2: 3/50} (it holds the diagonal's
   column-1 and -2 nodes). The test reads the diagonal from the file.
2. **NoData.** Test 1 with node (row 5, col 3) NoData (world (30, 50)): the
   triangle holding both west corners becomes {0: 16/25, 2: 17/50}, nodata
   1/50, first 10, count 3; the other is unchanged;
   `dem_height_triangles_without_data` 0.
3. **All NoData, Python level** (`height_shares` with `TileWindows` over a
   hand-built 3 × 3 tile, 10 m, centre node NoData, `tolerance_m=10,
   class_m=1`): a triangle inside the centre pixel (index corners (1.1,
   1.1), (1.4, 1.1), (1.1, 1.4)) gets count 0, first 0, nodata 1.0, and
   1 triangle without data.
4. **Class edges and clamp, binding level.** 6 × 6 pixels, row_area 1;
   one triangle (0, 0), (6, 0), (0, 6) in index space, z 0 at its corners
   (it holds the 21 nodes with r + c ≤ 5); DEM 0 except (row, col) (0, 0)
   3.0, (0, 1) 10.5, (0, 2) −10.6, (1, 0) −10.5, (1, 1) −0.5, (2, 0) 0.5,
   and (5, 5) 7.0 outside. Area 18: first 0, count 21, {−10: 2/18 (−10.6
   clamped, −10.5 on its lower edge), 0: 13/18 (−0.5 on its lower edge),
   +1: 1/18 (0.5 on its lower edge), +3: 1/18, +10: 1/18 (10.5 clamped)}.
5. **No held node.** 3 × 3 pixels of 0; triangle (1.1, 1.1), (1.4, 1.1),
   (1.1, 1.4) with z 0, 3, 0: it holds no node; the plane is 4.0 at the
   pixel centre (1.5, 1.5): first 6, count 1, {−4: 1}. (At the piece
   centroid it would be 1.0, class −1.)
6. **Plane as written.** Test 4 with z 1.0 at all corners: first 0, count
   21, {−10: 2/18, −1: 13/18, 0: 1/18, +2: 1/18, +10: 1/18}.
7. **Row areas.** 2 × 2 pixels, row 0 at 0 m, row 1 at 2 m, row_area
   (1, 3); triangles (0, 0), (2, 0), (0, 2) and (2, 0), (2, 2), (0, 2),
   z 0: each first 10, count 3, {0: 1/2, 2: 1/2} and {0: 1/10, 2: 9/10}.
8. **Resolution, binding level.** z = 3 sin(2πx/400) cos(2πy/300) m on
   [0, 1000]² m, 10 m (101²) and 5 m (201²) nodes, the 32 triangles of a
   250 m grid split nw-se, corners at z(x, y): classes −5 … +5 occur; each
   triangle's count is the same from both DEMs; the distance between the
   two histograms on the shared classes, measured median 0.012, largest
   0.018 (probe); the test asserts at most 0.03.
9. **Each pixel's area once, property.** C++ (`tests/cpp/property`): 300
   random triangulations of a rectangle in index space with off-node fans;
   for 20 random pixels each, the pixel NoData and the rest 0 with corner
   z 0: Σ over triangles of (plan area × nodata share) equals the pixel's
   area inside the rectangle within 1e-6 relative. (NoData, not a high
   pixel: departure c folds a high pixel's outside pieces into class 0.)
10. **Cut run.** `pieces_fixtures` `hills()`, `SQUARE`, tolerance 1,
    `--pieces 4 --dem-heights`: 4 piece files, 1,163 triangles, n = 1;
    every piece file states the edges −1.5, −0.5, 0.5, 1.5; a Python
    oracle (shapely clip in index space, residual at the pixel centre,
    held nodes by 3.2's test, the span rule) over every triangle of every
    piece: the same first and count, shares and nodata within 1e-5.
11. **Reprojected path, CLI.** `tests/fixtures/velhas`, `--out-crs
    EPSG:31983 --tolerance 5 --dem-heights`: 969 triangles, w = 1 m,
    n = 5; the same oracle in the source's index space (corners reprojected
    to EPSG:4674, `row_areas` per row): shares within 1e-4, and first and
    count agree for every triangle (both map the corners the same way). Spread
    (each triangle's standard deviation of its class centres k · w) and
    class-0 share, weighted by plan area in EPSG:31983: 1.164 m and 0.228
    (probe), asserted within 0.01 and 0.005; the resampled grid gives
    1.121 m and 0.262, outside both.
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
    offsets, first classes, shares and nodata shares; the edges array,
    comment and table give the same 22 edges, −10.5 … +10.5; `.vtk` lines
    get nodata 0 and no rows; `read_mesh` returns the `HeightTable`;
    `dem_height_class_m` 1.0 and `dem_height_class_from` the default's
    words in the record.
15. **Reader refusals.** A hand-made `.vtk` with `dem_height_shares` of
    type `double`, offsets ending one short of the shares, offsets
    decreasing, `dem_height_first` one short, a row whose first + count
    passes the classes, `dem_nodata_share` of type `double`, edges of odd
    length, and edges not increasing: `ValueError` naming the array; a
    file without them gives None.
16. **CLI.** `--dem-heights` in `rasputin mesh --help` and the MCP tool's
    arguments.
17. **Refusals.** Test 1's DEM with `--stride 10 --dem-heights` (no
    `--tolerance`), with `--tolerance 0 --dem-heights`, and with
    `--tolerance 200 --dem-heights` (401 classes): non-zero exit before
    any file, the message naming `--tolerance`; the binding with
    `class_m` 0, −1 and NaN and with `half` −1 and 128: `ValueError`;
    `relative_vertical_accuracy_m` on a features source and on a class source:
    the source file's `ValueError`.
18. **The count follows the tolerance, CLI.** Test 1 with `--tolerance 2`
    (measured: still 2 triangles, largest difference 2 m): n = 2, edges
    −2.5 … 2.5 (6); the triangle holding both west corners first 2, count
    3, {0: 33/50, 2: 17/50}, class +2 the outermost (the 2 m nodes at
    exactly τ inside), the other first 2, count 3, {0: 47/50, 2: 3/50}.
19. **Pieces centred outside fold, binding level.** 5 × 5 pixels of 0
    but (0, 0) 2.0 and (1, 3) 5.0; triangle (0.2, 0.2), (4.0, 0.2),
    (0.2, 4.0), z 0, area 7.22: it holds node (0, 0) (class 2) and nodes
    of class 0; pixel (1, 3)'s piece of area 0.02 is centred outside
    (3.5, 1.5) and folds into class 2: first 10, count 3, {0: 328/361,
    2: 33/361}. Without departure c: count 6, a share 1/361 at +5.
20. **The width is the source's, CLI.** Test 1's DEM registered as a local
    source with `relative_vertical_accuracy_m = 2.5` (a sources file under
    `RASPUTIN_DATA` in `tmp_path`), `--dem <key> --stride 10 --tolerance
    10 --dem-heights`: n = 4, edges (j + ½) · 2.5, j = −5 … 4 (10, ±11.25);
    the 2 m nodes are class +1: the west triangle first 4, count 2,
    {0: 33/50, +1: 17/50}, the other first 4, count 2, {0: 47/50, +1:
    3/50}; record `dem_height_class_m` 2.5 and `dem_height_class_from`
    naming the key.
21. **Convergence, CLI.** `write_dem(hills(n=201))`, `--stride 200
    --dem-heights` at `--tolerance` 4, 2, 1 and 0.4 (`.vtk`): 1,546,
    11,279, 42,457, 65,379 triangles; one-class share 0.6, 20.6, 73.1,
    100.0 % and mean count 4.875, 2.492, 1.283, 1.000 (probe); asserted:
    the mean count strictly falls, the one-class share never falls, every
    count is at most 2n + 1 (9, 5, 3, 1), and at 0.4 every count is 1.
    The one-class share is monotone only from tolerance 4 down: at 8 and
    4 it counts 4 and 9 triangles of 50-350 m² (1-4 pixels; probe), so it
    moves by whole slivers (0.7 → 0.6 %), not by resolution.
22. **`half_count` and `edges_m`.** (τ, w) → n: (1, 1) 1, (0.4, 1) 0,
    (0.5, 1) 1, (2.5, 1) 3, (10, 1) 10, (10, 2.5) 4, (5, 1.216) 4,
    (127.4, 1) 127; (127.5, 1), (1, 0), (0, 1), (NaN, 1): `ValueError`.
    `edges_m(1, 1.0)` is (−1.5, −0.5, 0.5, 1.5) exactly.

**Invariant-critical:** the overlap weights, the residual and the span.
Mutation targets: in `height_shares.hpp` the two class-edge comparisons,
the division by w, the clamp to ±half, the pixel-centre evaluation, the
held-node test, the fold to the span's ends, the no-node class, the NoData
branch and the denominator; in `dem_heights.py` `half_count` and
`edges_m`; in `overlap.hpp` the whole-pixel test (`full_begin`,
`full_end`) and the piece area.

**Mutants** the suite must kill: `<=` for `<` at a class edge (4); classes
not centred, `floor(r)` (4); no clamp to ±half (4); residual at the piece
centroid (5, 10); the span over all pieces (19, 21); folding to the wrong
end (19); the no-node triangle given class 0, or its pieces' spread (5);
held-node test strict, no margin (1: the corner and diagonal nodes); the
width from the tolerance, τ / 10 (1, 18, 20); the source's width ignored
(20); n by a non-strict margin (22, (2.5, 1)); edges not stated per run
(10, 18, 20); offsets or first off by one (14, 15); a cut pixel counted
whole, or a whole pixel clipped wrongly (1, 9); the plane from the DEM,
not the written z (6); row_area ignored (7); NoData counted in class 0
(2, 3, 9); the denominator without NoData (2: 32/49 for 16/25); the
resampled grid on the reprojected path (11); a shared accumulator across
threads (13); the flag ignored in `_write_piece` (10); anything written
when off (12); `class_areas`'s visitor reordered (12).

## 6. Degeneracies

A triangle smaller than a pixel: one or a few pieces; it usually holds no
node, so it has one class (3.1). A tolerance below w/2: n = 0, one class
everywhere, written as at any tolerance (3.1). A tolerance of 0: refused
with `--dem-heights` (3.1); without the flag it runs as today. All its
pixels NoData: count 0, nodata 1.0, counted (3.3). A held node on NoData
joins no span; a triangle whose held nodes are all NoData but which has
valid pieces gets the no-node class. A zero-area triangle never reaches the
file and is refused by the binding. A corner past the DEM window:
`OutsideWindow`, a defect upstream. A reprojected vertex with no image
(`inf`): the chunk window's `ValueError`, before the binding.

## 7. Cost and memory at 637,000 km², 30 m

7.1e8 DEM nodes, 23 times the São Francisco piece of section 2.

- **Time (estimate):** the walk scales with pixels and cut pixels; the
  held-node test adds three products per cut piece and the fold 2n + 1
  adds per triangle (63 at n = 20 against tens of pixels), so the first
  revision's figures stand to within noise: 0.35-0.46 s on the piece at 10
  threads, about 10 s for the basin, plus 1.2 s to reproject its vertices,
  against the piece's 12 s mesh run times 23. Prediction for the fixture
  command at `--tolerance 10`: +0.05 s for the walk and +0.01 s to encode
  2 MB, on a 0.23 s run.
- **Memory:** 12 bytes per triangle (offset, first, nodata) plus 4 per
  share while held: 27 bytes per triangle on the fixture at τ = 10 (mean
  3.79), at most 96 at n = 10 and 176 at n = 20; twice that while the
  per-block buffers are joined, plus 72 bytes per triangle of a chunk
  (corners and z, 2^18 triangles, 19 MB) while measured. On the reprojected
  path, the source tile held through the final check: 4 bytes per source
  node beside the check-point store's 16 (`check_points.hpp`, `Packed`).
  The piece budget's estimate is not changed; `@perf` reports the peak with
  and without.
- **Files:** binary `.vtk` and `.ply`, 12 bytes per triangle plus 4 per
  share.

## 8. Acceptance (timings reported, not gated)

On mains power, median of 3, against `98fe80c3`: the fixture command
(`--tolerance 10`, `.vtk`) with and without `--dem-heights`, the same cut by
`--pieces 4`, and the Velhas fixture with `--out-crs EPSG:31983`; time and
peak resident memory each.

## 9. Size

Estimated +507 counted lines (third revision: the held-node test, the fold
and the sparse join +35, `HeightTable` +30, the source field and width
resolution +15, the record entries +5; the fixed-count arrays −15):
`height_shares.hpp` +120, `overlap.hpp` refactor +15, binding and
`_core.pyi` +50, `dem_heights.py` +60, `io/models.py` +30,
`land_cover_fractions.py` row-area helper +8, `sources.py` +6, `cli.py` +70
(flag, width, `_refined`, refusals, piece and single writes, record),
`dem_input.py` +8, `refine_job.py` and `basin_run.py` +10,
`io/vtk_legacy.py` +30, `io/ply.py` +25, `io/geopackage_write.py` +35,
`io/vtk_read.py` +35, MCP arguments +5. One PR.

## 10. Questions for Ola

None. The width's kind (relative, not absolute) was the main session's
ruling, 2026-11-04.
