# Increment 34: a tighter vertical tolerance where the ground is steep

**Status:** designed (`@architect`, 2026-10-09). Design review round 1 asked
for changes (five blockers, ten suggestions), round 2 for two one-line fixes
and five suggestions, round 3 for two holes around the NoData class; all are
answered in this file and its probes; round 4 found one figure, fixed;
round 5 approved. Built on branch `worktree-slope-tol`: red `3ce92a7e`,
green `952bd1ee` (512 counted lines against 376-533), red-step scaffolding
removed `f68a9440`, mutation round `7d661ae9` (M1-M9 all killed); the
as-built record is sections 9.1-9.3 and 10's *As built* table. Next: code
review, then `@perf`'s acceptance run, which adds the `romsdal-slope` case to
`docs/benchmarks/quick/cases.toml` and so turns the quick-check pin
(`test_the_shipped_cases_are_sections_3_table`) green. Questions for Ola in
section 12; the design is written on their defaults, taken while Ola was
out.

## 1. What Ola asked

Ola, 2026-10-08: "Another idea; if the slope of a surface is steeper than a
certain angle, we want a finer refinement." Asked whether finer means a tighter
vertical tolerance or smaller triangles, and how the slope should be measured:
"vtol means finer. Second question is very good, and needs to be researched.
Deriving slope from a coarse triangle might miss local steeps." 2026-10-09:
"perhaps you can work on the slope adaptive refinement in the meantime?"

"vtol" is the vertical tolerance, `--tolerance`: the largest allowed height
difference between the mesh and the DEM at a DEM node. Increment 33 made it
depend on the distance to named lines. This increment makes it depend on the
slope of the ground as well: on ground steeper than a set angle the mesh is
held closer to the DEM.

The answer to Ola's research question (section 3): the slope is measured at
every DEM node, from the DEM at its own spacing, and **every node is held to
the tolerance of its own slope**. No triangle's slope is used, so no triangle,
however coarse, can hide a steep node. The probe shows what the triangle's own
slope would miss: 0.45 to 1.0 % of the nodes end above the tolerance their
slope asks for (section 3, "The rules compared").

## 2. Prior art: legacy and literature

### Literature, read

- **Greedy insertion against a sup-norm bound**, the method `refine` is
  (increment 14): Garland and Heckbert, *Fast Polygonal Approximation of
  Terrains and Height Fields*, CMU-CS-95-181, 1995
  (`https://www.mgarland.org/files/papers/scape.pdf`, read as text through
  `pdftotext`). Its section 3 compares "importance measures": "We explored four
  categories of importance measures: local error, curvature, global error, and
  products of selected other measures." Its section 3.4, *Product Measures*:
  "The method we used was to combine one or more of the importance measures
  given above with some bias measures. [...] Using these product measures, we
  were able to achieve results which were only slightly poorer than those
  produced by the local error measure. However, product measures are more
  complex, and hence more expensive, than any of the measures discussed so far
  with the exception of global error." Its footnote 2: "We use products rather
  than sums because the units of measure of the constituent terms are
  generally unrelated."
  *What this increment takes from it:* the node inserted is the one with the
  largest error divided by its node's allowed error, which is a product
  measure in their sense (the error times a weight, 1 / allowed). *What
  differs:* they judged products by how well they approximate under one
  global threshold, and found them slightly worse; here the threshold itself
  varies by node, and the product is what makes the inserted node one that
  breaks its own bound. Their stopping rule, one threshold for all, becomes
  "no node above its own allowed error".
- **Curvature rather than slope.** The same report, section 3.2: "The
  piecewise-linear reconstruction effected by T approximates nearly planar
  functions well, but does more poorly on curved surfaces." A plane at any
  slope is fitted exactly by one triangle; the vertical error a TIN makes is
  driven by curvature, not slope. Slope here does not predict error; it
  selects where the user wants a smaller error. That is Ola's intent as
  stated, and the reason the driver is a tolerance, not an importance
  measure.
- **Mesh size from the slope of the bed, in coastal ocean meshing.** Roberts,
  Pringle and Westerink, *OceanMesh2D 1.0*, Geoscientific Model Development
  12:1847-1868, 2019 (`https://gmd.copernicus.org/articles/12/1847/2019/`,
  the PDF read as text through `pdftotext`), section 3.3.3, the topographic slope
  mesh size function `h_slp = (2 pi / alpha_slp) b / |grad b|`: "where
  2π/αslp is the number of elements that resolve the topographic slope, and
  ∇b is the gradient of the bathymetry evaluated on a structured grid". On
  noise: "Typically the gradient of the bathymetry often contains a high
  degree of noise, which results in high mesh refinement" and "Therefore, we
  propose low-pass filtering the bathymetry before calculating the gradient."
  The filter's cutoff there is a physical length of the ocean (the Rossby
  radius). *What differs:* they size elements horizontally from the gradient;
  here the vertical tolerance follows it, and the triangle size follows from
  greedy insertion. Smoothing is not taken over: the noise probe (section 3,
  "Noise") finds it does not matter on DTM10.
- **The slope formula.** Horn, *Hill shading and the reflectance map*, Proc.
  IEEE 69:14-47, 1981, the 3 by 3 weighted difference used by GIS packages
  (found, not read; the formula is written out in section 3). GDAL's `gdaldem`
  manual (`https://gdal.org/en/stable/programs/gdaldem.html`, read): "The
  literature suggests Zevenbergen & Thorne to be more suited to smooth
  landscapes, whereas Horn's formula to perform better on rougher terrain."
  Jones, *A comparison of algorithms used to compute hill slope as a property
  of the DEM*, Computers & Geosciences 24:315-323, 1998 (found, not read):
  the search engine's summary says it compares eight algorithms and ranks the
  second-order finite difference first. The probe measured the two on DTM10
  (section 3, "Which slope"): their shares of nodes at or above each angle tried
  differ by at most 0.8 percentage points.
- **Slope as a quality measure of terrain simplification.** Biniaz, *Slope
  Preserving Terrain Simplification — An Experimental Study*, CCCG 2009
  (`https://cccg.ca/proceedings/2009/cccg09_16.pdf`, read): "This paper
  introduces a new measure of quality for terrain surface simplification
  aiming at preserving the slope of the surface"; "a simplification is
  considered “good” by slope measure if for most points p of output terrain,
  the slope at p on the simplified terrain is not significantly different
  from the corresponding slope on the original terrain." It also names the
  exact vertical bound this project keeps: "specifying an exact numerical
  error tolerance ε such that the simplified terrain must lie within vertical
  distance ε of the original, at every point". *What differs:* that work
  preserves slope as the output quality; this increment keeps the vertical
  bound and only lets its value depend on slope.
- **Hydrological TINs.** Vivoni, Ivanov, Bras and Entekhabi, *Generation of
  Triangulated Irregular Networks Based on Hydrological Similarity*, Journal
  of Hydrologic Engineering 9(4):288, 2004 (the abstract, as OpenAlex
  rebuilds it from its word index, which can drop words; the full text was
  not reachable): three ways of building TINs, "(1) Traditional, (2)
  hydrographic features (3) hydrological similarity TINs", the last "using
  the concept of hydrological similarity provided by a topographic or wetness
  index". The wetness index of Beven and Kirkby (1979) contains the local
  slope (recalled, not re-read). This is resolution set by a terrain
  attribute that includes slope, by a different mechanism (a sampling filter
  on an index, not a vertical bound).
- **DEM error grows on steep ground.** Hutchinson, Xu and Stein, *Recent
  Progress in the ANUDEM Elevation Gridding Procedure*, Geomorphometry 2011
  (`http://www.geomorphometry.org/uploads/pdf/pdf2011/HutchinsonXu2011geomorphometry.pdf`,
  read): "all remotely sensed elevation data have significant random errors
  that depend on the inherent limitations of the observing instrument, as
  well as surface slope and roughness". ANUDEM itself grids elevation data
  ("by minimizing the sum of a user-specified roughness penalty and a
  weighted sum of squares of the residuals"); it builds no TIN and has no
  slope-dependent tolerance. Its relevance here: a tolerance tighter than the
  DEM's own error on steep ground buys fidelity to the DEM's noise, not to
  the ground. Section 12, Q5, says so to Ola.
- **Distance-driven tolerance**, increment 33
  (`docs/increments/33-feature-tolerance.md`, section 2): Gmsh's `Threshold`
  and `Min` fields, view-dependent level of detail. The ramp below is the
  same ramp with the slope in place of the distance, and "the tightest of the
  drivers" is Gmsh's `Min`.
- Found, not read: *3D Simplification Methods and Large Scale Terrain
  Tiling*, Remote Sensing 12(3):437, 2020 (MDPI refused the fetch); the
  search engine's summary says greedy insertion with a height error adds
  more vertices on very steep ground. The probe measures that directly
  (section 3: uniform 10 m leaves 17.7 % of Romsdalen's nodes above the
  slope tolerance, 2.4 % of Geilo-Ål's).

### Novelty

None is claimed. A tolerance or a mesh size that varies over the domain is in
all the work of increment 33's section 2; one that follows the slope is in the
coastal meshing work above; a product of the error and a weight is Garland and
Heckbert's. What is specific here is small: a sup-norm greedy insertion in
which every DEM node is held to the tolerance at its own DEM slope, exactly
(an exact comparison per node, not a ratio test; section 4.3), composed with a
per-triangle distance tolerance. Searches for an error-bounded TIN whose
vertical bound follows the slope ("TIN generation from DEM vertical error
threshold depends on slope adaptive tolerance steep terrain greedy insertion";
"terrain TIN maximum vertical error varies with slope tolerance per point
adaptive mesh steep slopes tighter tolerance") found nothing of the kind; that
is weak evidence and is not a claim. If Ola publishes, this is an application
of known ideas.

### Legacy

Nothing to carry across. The legacy tree computed slopes of the output's
triangle normals, to colour the mesh; no refinement read them:

```
$ git grep -l -i -E "slope|steep|gradient" legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/geo_tiff_reader.py
legacy-archive:legacy/rasputin/geometry.py
legacy-archive:legacy/rasputin/mesh_colouring.py
legacy-archive:legacy/rasputin/solar_position.h
legacy-archive:legacy/rasputin/triangulate_dem.h
```

The hits are `compute_slopes` from face normals
(`legacy/bindings.cpp@legacy-archive:372`), a `-slope` output field in
`legacy/rasputin/geo_tiff_reader.py@legacy-archive:26`, a cached `slopes`
property, and colour bands in `legacy/rasputin/mesh_colouring.py@legacy-archive:26-36`
(dark from 55 degrees, red from 30). Those bands are the only domain hint: the
legacy viewer marked 30 degrees, which is also this design's example
threshold.

## 3. The slope, and the rule each node is held to

**The slope of a node** is Horn's: with the node's eight neighbours `a b c /
d . f / g h i` (rows top to bottom), `gx = ((c + 2f + i) - (a + 2d + g)) /
(8 dx)`, `gy = ((g + 2h + i) - (a + 2b + c)) / (8 dy)`, slope `atan(hypot(gx,
gy))`, from the DEM refine runs on, at its own spacing (`dx`, `dy` may differ;
the CLI fixtures' grid is 10 by 5 m).

*A missing neighbour* (outside the grid, or NoData) is filled so that a plane
reads its own slope:

- an edge neighbour (`b`, `d`, `f`, `h`) by its reflection through the node,
  `2 z(node) - z(opposite)`, when the opposite neighbour is there; else by
  `z(node)`;
- a corner neighbour (`a`, `c`, `g`, `i`) by its reflection through the node
  when the opposite corner is there; else by `z(row neighbour) + z(column
  neighbour) - z(node)`, from the edge neighbours as filled above.

On any plane this is exact at every node that has a neighbour on at least one
side along each axis: the border, the four corners, and nodes beside NoData
(`34-probes/README.md`, item 1b: a plane at 37.3 degrees on the 10 by 5 m grid,
with a NoData node in it, reads 37.300 degrees at every valid node). *The one
departure:* a node with neither edge neighbour along an axis. Its two edge
neighbours on that axis take `z(node)`, so that axis's part of the slope comes
from the corner neighbours alone, at half Horn's weight, where they are valid
(a 40-degree plane facing north, NoData above and below the node, reads
22.8 degrees), and is lost where they are missing too, as in a one-row or
one-column grid (36.0 degrees for a 40-degree plane facing 30 degrees off the
row); both in item 1b. Either way it reads low. Round 1's rule (a missing neighbour
takes the node's own z) read a 40-degree plane facing 30 degrees as 30.3
degrees at a border node, 12.9 at a corner and 34.5 beside a NoData node
(item 1b). Round 1's probes computed neither rule (NumPy's edge padding
repeats the outermost row, not the node); on the Romsdalen window that put
27.0 % of the outer ring at 30 degrees or more against 41.0 % with this rule
and 38.2 % inside, and the two differ only on the outer ring (item 1b). Every
probe now computes this rule (`34-probes/slope_stats.py`, `horn`). A NoData
node gets the reserved class 255 (4.5); its error is never measured.

**Classes.** The slope is stored as a class, one byte per node: the slope in
half degrees, **rounded up**, 0 to 180. Rounding up means the stored slope is
never below Horn's, so the tolerance it gives is never looser than the ramp
at Horn's slope. *Constant:* 0.5 degrees, an angle, independent of scale;
checked on DTM10 at 2.55e7 nodes (tile `6901_3`, `34-probes/README.md`,
item 5). The classes are computed by comparing `gx^2 + gy^2` with a table of
`tan^2(c / 2 degrees)`, no `atan` per node; that is `@developer`'s choice to
make or not, the rule is "the smallest class whose angle is at least the
slope".

**The ramp.** With the steep tolerance `N`, the general one `F`
(`--tolerance`) and two angles `START <= END` in degrees, the allowed error at
slope `s` is

```
t(s) = N                                      for s >= END
t(s) = F                                      for s <= START (and s < END)
t(s) = F + (N - F) (s - START) / (END - START) for START < s < END
```

tested in that order, so a step (`START = END`) gives `N` from `END` up and
`START = END = 0` holds every node to `N`; evaluated at `s = class / 2`.
Requires `0 < N <= F` and `0 <= START <= END < 90`. `N = F` is today's mesh
(G2). It is increment 33's ramp turned round: tight at the far
end of the variable, not the near one.

**The rule: every DEM node is held to the tolerance at its own slope.** A
triangle converges when every node of its scan set (the DEM nodes in the
closed triangle, as today) has error at most `t(s(n))`, and every check point
or strip point (the resampled path and the edge strip, increments 15c and 15f)
has error at most `t` of the **cell it is filed in**: the largest class of the
cell's valid corners, the DEM nodes at its corners that are not NoData (0 when
none is valid; 4.5), in the cell `CheckPoints` files the point in
(`include/terrain/refinement/check_points.hpp@9ba38490:10-13`). When a
triangle does not converge, the node inserted is the one with the largest
`error / t(s(n))`; ties go to the larger error, then to the first in today's
order (section 4.3). With one tolerance everywhere that is today's worst node,
so `N = F` is today's mesh bit for bit (G2).

*Why per node.* Three other rules were simulated (`34-probes/greedy_sim.py`, a
Python model of greedy insertion; `34-probes/README.md`, item 2). On two
10 by 10 km windows of DTM10, `N = 2`, `F = 10`, a step at 30 degrees:

| rule | Romsdalen (37 % of nodes at 30° or more) | Geilo-Ål (4.8 %) | nodes above their own tolerance |
|---|---|---|---|
| uniform `F` (today, `--tolerance 10`) | 56 193 | 16 460 | 17.7 %, 2.4 % |
| **per node** (this design) | **279 337** | **42 016** | 0 |
| per triangle, the steepest node in it | 321 830 (+15 %) | 59 394 (+41 %) | 0 |
| per triangle, the steepest node of every cell it meets | 344 937 (+23 %) | 73 816 (+76 %) | 0 |
| per triangle, the slope of its own plane | 271 457 | 39 019 | 1.0 %, 0.45 % |
| uniform `N` (`--tolerance 2`) | 412 641 | 197 623 | 0 |

The simulation's own counts (it overcounts rasputin by 28 to 57 %; item 3);
the relative order is what is used. The per-triangle rules keep the
guarantee but hold a whole triangle to the tolerance of its steepest node,
which costs most where steep ground is patchy, as around Geilo: +41 % and
+76 %. The last of them is the rule needed to cover check points
on the resampled path, so it is the per-triangle rule one would actually
build. The triangle's own plane is cheapest and is what Ola warned against: a
coarse triangle across a steep slope and a flat shelf reports a middling slope,
and 1.0 % (Romsdalen) and 0.45 % (Geilo-Ål) of the nodes end above the
tolerance their own slope asks for.

This departs from increment 33's sketch of a slope driver
(`docs/increments/33-feature-tolerance.md@9ba38490:414-418`: a policy whose
per-triangle `at(m, t)` reads "its steepest cell, say"). That sketch is one of
the two per-triangle rows of the table; the per-node rule needs fewer
triangles for the same guarantee.

*Why a byte per node, not a per-triangle query.* The scan already walks every
node of the triangle in contiguous row segments
(`include/terrain/refinement/scan.hpp@9ba38490:204-250`); the class row is a
second contiguous byte row at the same index, so the weight rides on the same
pass. Recomputing Horn's slope per node per scan would read nine heights per
node in every scan; a precomputed class costs one byte per node once.

**Which slope** (`34-probes/README.md`, item 1). On Romsdalen and on 33's whole
Geilo-Ål section, the share of nodes at or above 15 to 45 degrees by Horn's
formula and by centred differences (Zevenbergen and Thorne's gradient) differ
by at most 0.8 percentage points (Geilo-Ål at 15 degrees: 18.8 against
19.6 %). The largest bilinear cell gradient (the bound `foot_epsilon` uses,
`include/terrain/refinement/refine.hpp@9ba38490:240`) is an upper bound, not a
slope: it puts 51 % of Romsdalen at 30 degrees or more against Horn's 37 %.
Horn it is: the two agree on DTM10, and GDAL's manual (section 2) leans to
Horn on rough ground.

**Noise and smoothing** (item 1). Gaussian noise added to DTM10's z: 0.1 m
moves Horn's slope by 0.16 to 0.19 degrees on average, 1.0 m by 1.57 to 1.94
degrees, and the share at 30 degrees or more by at most 0.2 percentage points.
Horn's differences are linear in z, so the slope noise scales as noise over
spacing: 1.0 m on the 10 m grid stands for 0.1 m on a 1 m grid (DTM1). A 3 by
3 mean filter before Horn changes the shares by at most 1.5 percentage points
(Geilo-Ål at 15 degrees). The per-node rule is also the one that noise hurts
least: a node's tolerance depends on its own slope only, where any
per-triangle rule takes the largest slope over many nodes, and the largest of
many noisy values is pulled up by the noise. No smoothing in this increment;
a smoothing radius for DTM1 is a later flag if a DTM1 run shows the need
(section 11).

**The resampled path** (`34-probes/README.md`, item 1c). With `--out-crs`
refine runs on the target grid, square at the source's north-south spacing
rounded to whole metres (`src_python/tin_engine/target_grid.py@1a754deb:91-102`,
`default_spacing`), and the slope is that grid's. For a geographic source the
east-west spacing is finer than that by the cosine of the latitude, so the
target grid smooths east-west. On Copernicus GLO-30 at N47 E011 (the
Wetterstein and Karwendel; source 20.96 m east-west by 30.88 m, target 31 m
square in EPSG:25832, bilinear from the source): at 30 degrees or more 41.0 %
of the target's nodes against 42.2 % of the source's, at 40 degrees 15.2
against 15.5 %, the median 26.7 against 27.2 degrees. The resampled path
therefore tightens slightly less ground than the source would; the target is
the grid the check points are filed in, so it is the one used (Q4).

**Step or ramp.** On Romsdalen the ramp from 25 to 35 degrees gives 253 061
simulated triangles against 279 337 for the step at 30, 9 % fewer
(item 2). The flag takes both angles, so either is one command (section 5).

## 4. Where it sits

### 4.1 Data flow

```
--tolerance F, --tolerance-slope N START END             (cli.py)
        |  ToleranceSlope (Pydantic, frozen): near_m, start_deg, end_deg
        v
_core.SlopeTolerance(view, N, F, START, END, threads)    (binding; the GIL released)
        |  terrain::raster::steepness(view, threads): one byte per node, parallel by rows
        |  two 256-entry tables: t(class / 2) for classes 0-180, F for 181-255, and the inverses
        v  immutable C++ object on the view's geometry, built once per run
refine(..., field=..., slope=...) / refine_strip(...) / refine_points(...)
   scan: per node, error vs t(class) exactly, and the largest error / t(class)
```

The C++ core receives no path, no CRS and no file (CLAUDE.md §2): the slope is
computed from the raster view refine already reads. On the resampled path
(`--out-crs`, increment 15c) it is the target grid's view, and the check points
of the final check use the target grid's classes (Q4; section 3, *The
resampled path*).

### 4.2 C++: `include/terrain/raster/steepness.hpp` (new)

```cpp
namespace terrain::raster {
inline constexpr std::size_t kSteepnessClasses = 181;  // half degrees, 0 to 90
inline constexpr std::uint8_t kNoDataClass = 255;     // a NoData node (4.5)
// Horn's slope of every node, as section 3's class (half degrees, rounded up),
// missing neighbours filled as section 3 says, kNoDataClass for a NoData node,
// row-major, rows() * cols() bytes. Parallel over blocks of rows; each byte is
// a function of the node's 3 by 3 neighbourhood only, so the result does not
// depend on the thread count.
template <RasterSource R>
[[nodiscard]] std::vector<std::uint8_t> steepness(const R& dem, unsigned threads);
}
```

It is a raster function, not a refinement one: it reads heights and the
geometry and nothing of the mesh, and it is the module a later slope output
(a per-triangle slope array, increment 10's "per-triangle datasets") would
reuse.

### 4.3 C++: `include/terrain/refinement/slope_tolerance.hpp` (new)

```cpp
namespace terrain::refinement {

struct SlopeRamp {  // metres and degrees; 0 < near <= far, 0 <= start <= end < 90, all finite
    double near, far, start_deg, end_deg;
    [[nodiscard]] double at(double slope_deg) const noexcept;  // section 3's t(s)
};

class SlopeTolerance {
public:
    // nullopt with one plain sentence in `why` for a ramp value that breaks
    // its bound or is not finite. Computes steepness(dem, threads).
    template <raster::RasterSource R>
    static std::optional<SlopeTolerance> make(const R& dem, SlopeRamp ramp, unsigned threads,
                                              std::string& why);
    const raster::RasterGeometry& geometry() const noexcept;
    std::span<const std::uint8_t> row(std::size_t r) const noexcept;    // the classes of one row
    std::uint8_t cell_class(mesh::MeshVertex p) const noexcept;         // largest valid corner class of p's cell, 0 if none
    double allowed(std::uint8_t c) const noexcept;  // ramp.at(c / 2.0) for c <= 180, else F; from the table
    double weight(std::uint8_t c) const noexcept;   // 1 / allowed(c), from the table
    std::array<std::size_t, 256> histogram() const;  // for the tests; entry 255 counts NoData
    double near() const noexcept;
    double far() const noexcept;
private:
    raster::RasterGeometry geometry_;
    SlopeRamp ramp_;
    std::vector<std::uint8_t> classes_;
    std::array<double, 256> allowed_, weight_;  // entry 255 (NoData): F
};

// A per-triangle policy (33's TolerancePolicy) with the per-node slope beside it.
// lowest(), highest() and at() are the triangle part's; refine reads `slope` in the scan.
template <TolerancePolicy P>
struct Sloped {
    P triangle;
    const SlopeTolerance* slope;
    [[nodiscard]] double lowest() const { return triangle.lowest(); }
    [[nodiscard]] double highest() const { return triangle.highest(); }
    [[nodiscard]] double at(const mesh::LatticeMesh& m, std::uint32_t t) const { return triangle.at(m, t); }
};
}
```

`Sloped<P>` models 33's `TolerancePolicy`, so the entry points keep their
signatures; refine tells it apart at compile time (`if constexpr`, a trait
`is_sloped<P>`), and `detail::varies` judges its triangle part, so
`Sloped<UniformTolerance>` keeps no `allowed` vector. Without a slope, every
instantiation is today's: G1.

`cell_class` is the cell `CheckPoints` would file `p` in, `(floor(row),
floor(col))` clamped to the last cell, so a check point and the class it is
judged by come from one rule. For a node it is still the cell's largest
valid corner, not the node's own class; nodes are judged by `row()` in the scan.

**The scan with a slope** (`scan.hpp`). A second result, filled only when the
policy is `Sloped`, kept in its own vector beside `results` so the default
path's `ScanResult` and its memory stay as they are:

```cpp
struct SlopeScan {
    bool over = false;                        // some node's error > allowed(class), compared exactly
    double ratio = 0.0;                       // the largest error * weight(class)
    std::optional<mesh::LatticeVertex> node;  // where; ties: larger error, then today's order
    NodeLocation where = NodeLocation::Inside;
};
```

- `over` is the exact comparison `error > allowed(class)`, never `ratio > 1`:
  `error * (1 / t)` can round to 1 for an error one unit in the last place
  above `t`. The guarantee rests on `over` alone; `ratio` only chooses the node
  (M1).
- The ranking key is the pair `(error * weight, error)`, compared
  lexicographically, strictly larger replacing, in today's row-major walk.
  Rounding is monotone, so with one weight for every node the product never
  reverses the order of two errors, and where it makes two different errors
  equal, the second key picks the larger: the node chosen is today's (G2, M2).
- The void branch (a NoData corner) is unchanged and fills no `SlopeScan`.
- The node walk keeps its frozen-edge unswitching
  (`include/terrain/refinement/scan.hpp@9ba38490:222-226`); the slope adds a
  byte row read beside each `RowSegment` and the ranking, under `if
  constexpr`, so the scan refine runs without a slope is compiled as today.

**The serial phase** (`refine.hpp`, and `point_loop` in `refine_points.hpp`):

- a triangle splits when `slope_scan.over`, or when today's test fails against
  the triangle part (`needs_split(r, limit(t))`,
  `include/terrain/refinement/refine.hpp@9ba38490:395`);
- the node inserted: `slope_scan.node` when `over`, else today's `r.node`;
- **laziness** (33's §4.4): when `over`, the scan stores `-infinity` in the
  `allowed` slot, where the triangle part keeps one (split, no query), as
  for an error above `highest()`; the triangle part is asked only for a
  triangle that no node's slope already splits. With `Sloped<UniformTolerance>` nothing is ever asked;
- **the foot's epsilon** (20b's `eps = clamp(tol / G, ...)`): `tol` is the
  smaller of the triangle part's allowed error and `allowed(class)` of the
  node being inserted. Where the slot holds `-infinity` (the scan skipped the
  query, because of `over` or an error above `highest()`), the foot asks the
  triangle part's `at(m, t)` there, as 33 does at
  `include/terrain/refinement/refine.hpp@9ba38490:417`; so an `over` slot
  costs one query when, and only when, its node has a constraint within the
  foot cap. Skipping that query would change the foot beside a line and break
  G8;
- `RefineOutcome` gains `max_slope_share`: the largest `ratio` over the final
  triangles (0 without a slope). By the guarantee it is at most 1, to
  rounding (the weight is rounded): the test bound is `1 + 1e-12`.

**Check points and strip points** (`refine_points.hpp`, `strip_scan.hpp`).
The same pair of results per triangle: a stored point's or a strip point's
class is `cell_class` of its position; a DEM node rescanned by `refine_strip`
is the scan's. Across the three sets the rule stays "first void wins", then
the larger key `(error * weight, error)`, ties to the earlier set
(`include/terrain/refinement/strip_scan.hpp@9ba38490:51-53`), and `over` is
the OR of the sets'.

### 4.4 Composition with increment 33's lines

A node is held to the smaller of the two: the lines' ramp at its triangle's
distance (33, per triangle) and the slope's ramp at its own slope (per node).
`--tolerance-near` and `--tolerance-slope` together build
`Sloped<LineTolerance>`. The two tests are independent and both must pass: a
node above its slope tolerance splits the triangle whatever the distance; a
triangle whose largest error is above the lines' allowed error splits whatever
the slopes. The tightest driver wins, as 33's §4.6 asked, without the slope
becoming a per-triangle number.

### 4.5 Python

`src_python/tin_engine/tolerance_field.py` gains the data, beside 33's
`ToleranceLines`:

```python
class ToleranceSlope(BaseModel, frozen=True):
    near_m: float     # N, > 0
    start_deg: float  # START, 0 <= START <= END
    end_deg: float    # END, < 90
```

`cli.py` checks the flag (section 5), builds the spec, and in `_dem_mesh`
builds `_core.SlopeTolerance(view, N, F, START, END, threads)` once from the
view refine reads (the target grid's on the resampled path), timed as a new
`--stats` phase, `slope`. It passes `slope=` to `refine`, `edge_strip.run`
and `final_check.run`, as 33 passes `field=`. Without the flag nothing is
built and every call is today's.

The binding's `with_policy` (`bindings/core.cpp@9ba38490:386-391`) dispatches
four ways (field or not, slope or not), and refuses a slope built on another
raster geometry: `ValueError("the slope tolerance was built on another raster
geometry")`, 33's pin B for the slope.

**The run record** (increment 25) gains, only with the flag:

| name | value | plain label |
|---|---|---|
| `tolerance_slope_m` | `N`, a number | Tolerance on steep ground |
| `tolerance_slope_deg` | `"25 to 35"`, text | Slope where the tolerance starts to tighten, and where it reaches the steep value |
| `slope_nodes_tightened` | `"<k> of <n> DEM nodes (<p> %)"`, text | DEM nodes held to less than the general tolerance because of their slope |
| `max_error_slope_share_of_tolerance` | a number, at most 1 | Largest error as a share of the tolerance its slope allows |

`slope_nodes_tightened` counts the DEM nodes inside the mesh (the domain, not
the canvas around it) whose class gives `t < F`, out of all valid nodes inside
the mesh. C++ counts it once, after the loop, over the final triangles: each
node in a triangle's closed scan set counted by one triangle only (a node on
an edge by the triangle of lower index of the two across it, or by the only
one), and each vertex that is a node once, from the vertex list. *Per node:*
a node strictly inside a triangle is visited once, a node on an edge twice
(by both triangles, counted by one), so at most two byte reads per node,
plus one exact on-edge test for the nodes at the ends of a row span. *In
parallel:* over fixed blocks of triangle slots, one pair of integer counts
per block, summed at the end; integer sums do not depend on the order, so
the result is the same at any thread count (G7). *NoData:* class 0 does not
tell NoData from flat ground, so the classes reserve the value 255
(`kNoDataClass`) for a NoData node: the count skips it, `cell_class` takes the
largest of a cell's valid corners (0 when none is), and the two tables have
256 entries. `make` fills entries 0-180 from the ramp and 181-255 with `F`
(and the weight `1 / F`): `steepness` writes only 0-180 and 255, so 181-254
are never read, and entry 255 is read only for a NoData node, which the scan
never ranks (its error is set to 0 before, as today). If `cell_class` took
the plain largest corner, 255 included (M9), a point beside NoData would be
held to `F`; hence the valid-corner rule and test 6 (b). It is
`RefineOutcome::slope_nodes` (two counts); on the resampled path the final
check computes it over the target grid's nodes, from the same classes.
`histogram()` stays, for the tests. On the DEM path
`max_error_slope_share_of_tolerance` is the larger of refine's and the edge
strip's `max_slope_share`, on the resampled path the final check's, as
`max_error_m` is today.

**The summary line** (`run_record.summary`) says today "Every DEM node inside
the mesh is within F m of it", which stays true. With the flag it gains one
sentence: "Nodes on steep ground are within the tolerance their slope allows
(largest error <s> of it)." with `<s>` the share as a percentage; and when the
share is above 1 a `Warning:` line, as `max_error_m` above the tolerance
gets today.

### 4.6 Locality, threads, pieces, memory

- The classes are computed once, read-only after, shared by every scan thread;
  the scan's extra work is per node and local. No global operation is added.
- **Pieces** (increment 23): each piece computes the classes of its own
  raster. They equal the global ones except on the raster's outer ring of
  nodes, where a missing neighbour is filled (section 3), which is exact
  only on a plane. A piece whose
  raster holds one node beyond its domain on every side gets the global
  classes everywhere it meshes. The seam pass's slope (on the seam strip's
  raster) goes in with whichever of 23c-2 and 34 lands second, as 33's seam
  distance does (33, §4.6); nothing on master calls `refine_seam` yet.
- **Memory**: one byte per node of the grid refine runs on, plus one
  `SlopeScan` per triangle slot with the flag (sizes below). Numedalslagen's
  canvas (2.80e8 nodes) needs 0.28 GB beside a run that peaks at 2.60 GB
  today (+11 %); the whole Romsdalen tile 25.5 MB; the São Francisco basin
  at 30 m about 0.71 GB in one piece (637 000 km² over 900 m²), less per
  piece. *Per triangle slot, as built* (`sizeof`, Apple clang 21, arm64,
  at `7d661ae9`; independent of scale; the slot vectors are as long as the
  mesh's triangle count):
  - `refine` with the flag: one `SlopeScan`, 40 bytes (`SlopeRank`'s flag,
    two `double`s, an optional node and its location); at 3.5 million
    triangles (section 6's high estimate) about 0.14 GB. `ScanResult`
    stays 32 bytes.
  - `refine_points` (the final check) and `refine_strip` (the edge strip)
    with the flag: a second `PointScan` per slot (`steep`), 120 bytes; at
    3.5 million triangles about 0.42 GB.
  - Every run, flag or not: `PointScan` grew from 104 to 120 bytes
    (`SlopeRank` is its base, 9.3 pin 4), so the final check and the edge
    strip use 16 bytes more per slot than before 34; at 3.5 million
    triangles about 0.056 GB.

### 4.7 Against the four criteria

1. *Data and execution:* `ToleranceSlope` and `SlopeRamp` are data; the
   classes are computed from them once and only read.
2. *State:* no global state; the classes and the tables are immutable after
   `make`; the scan stays pure.
3. *Dependencies:* none added. The probe's SciPy is not a rasputin dependency
   and nothing in the product calls it.
4. *Async:* the classes are computed in the binding with the GIL released, on
   the thread that then calls `refine`, as the field is today.

## 5. Interface

```
rasputin mesh --dem DTM10 --domain basin.geojson \
    --tolerance 10 --tolerance-slope 2 25 35 --out basin.vtk
```

- `--tolerance-slope N START END`: "Hold DEM nodes on steep ground to N metres:
  the general --tolerance up to START degrees of slope, N from END degrees,
  linear between (START = END: a step). The slope is the DEM's own, at each
  node." Three numbers (a Typer tuple, as `--tolerance-ramp`).
- Refusals, in the CLI's existing style (33's §9.3, pin 5), checked in this
  order: "applies only with --dem"; "--tolerance-slope needs --tolerance";
  "must be finite, got <values>"; "N must be above 0, got <n>"; "N <n> is above
  --tolerance <F>"; "START <S> is below 0"; "START <S> is above END <E>";
  "END <E> must be below 90". Each with `--tolerance-slope` as the flag named.
- With `--tolerance-near` too: both apply (section 4.4).
- `N = F`: the CLI builds no slope (every node would get `F`), so the mesh is
  `--tolerance F`'s by construction, and says so on stderr: "--tolerance-slope:
  N equals --tolerance, so slope changes nothing". C++ still has G2's test.

## 6. The demonstration case: Romsdalen

Tile `6901_3` of DTM10 (`../rasputin_data/DTM10_UTM33_20260925`, 2 550 km²,
5051 by 5051 nodes): Romsdalen, Isfjorden and the walls below Trollveggen; 22 %
of its nodes are at 30 degrees or more, 31 % at 25 or more (item 5). One
command, no domain file:

```
rasputin mesh --dem DTM10_UTM33_20260925/6901_3_10m_z33.tif \
    --tolerance 10 --tolerance-slope 2 25 35 --out romsdal-slope.vtk
```

**Expected size and time** (`34-probes/README.md`, items 3, 5 and 6; Mac on AC):

| case | uniform 10 m | uniform 2 m | with `--tolerance 10 --tolerance-slope 2`, estimate |
|---|---|---|---|
| Romsdalen window, 100 km² | 38 935 triangles | 321 985 | about 195 000 (ramp 25-35), 216 000 (step 30) |
| Geilo-Ål window, 100 km² | 10 502 | 143 783 | about 29 000 (step 30) |
| tile `6901_3`, 2 550 km² | 615 459, refine 0.38 s, 0.68 s in all, 0.47 GB | 5 916 344, refine 3.74 s, 6.32 s in all, 2.19 GB | 1.4 to 3.5 million |

The estimates place each rule where the simulation places it between its own
uniform counts, applied to rasputin's uniform counts (`estimate.py`). The
tile's range runs from Geilo-Ål's position (0.14) to Romsdalen's (0.55); its
steep share lies between the two windows'. At the high end, refine at the
uniform 2 m rate (0.63 µs per output triangle, wall) and section 10's
per-triangle slowdown give about 2.4 to 2.6 s of refine; the whole run about
4 to 5 s. These are estimates; `@perf` measures the real ones.

**With 33's lines** (item 6). On the Geilo-Ål window, with the quick check's
`geilo-al-ramp` flags (`--tolerance 20`, the Bergen Line at 1 m, ramp 0 to
3000 m) and `--tolerance-slope 2 25 35`: rasputin meshes the lines alone in
3 687 triangles; the simulation puts the lines with the slope at 24 288
against its 5 667 for the lines alone and 473 736 at uniform 1 m, which
places the combined mesh at about 18 000 triangles. Here the slope does most
of the work: 11.1 % of the window lies within 3 km of the line and 2.1 % of
its nodes are at 35 degrees or more, many far from the line, where the lines
alone allow 20 m. The simulation measures the lines' distance from the
triangle's nodes and corners, at least 33's exact distance, so it errs low.

The tests use synthetic DEMs (section 9), not this data.

## 7. Guarantees and the checks that enforce them

| # | guarantee | checked by (section 9) |
|---|---|---|
| G1 | Without `--tolerance-slope` every output byte is today's | test 1, the bench run |
| G2 | With `N = F` the mesh is the `--tolerance F` mesh, bit for bit (C++); the CLI builds no slope then | tests 4 (a), 12 |
| G3 | Every valid DEM node's error is at most `t(class(n) / 2)`, and with lines at most the smaller of that and 33's ramp at its triangle's distance | test 5 |
| G4 | Every check point's and strip point's error is at most `t` of its cell's class (the largest valid corner, 0 if none), except a point on a frozen edge (increment 23's seam: never named, counted in `on_frozen`, as today) | test 6 |
| G5 | `class(n)` is the smallest half degree at or above Horn's slope of n, missing neighbours filled as section 3 says; at a class boundary, within rounding, the higher class | test 2 |
| G6 | The slope reaches every comparison: `START = END = 0` (every node steep) gives the `--tolerance N` mesh, bit for bit, feet included | test 4 (b) |
| G7 | The output is the same for any thread count | test 7 |
| G8 | With lines, the slope adds only its own condition: slope with `N = F` and lines equals the lines alone, bit for bit | test 4 (c) |

## 8. Degeneracies

- A flat DEM (every class 0): `t = F` everywhere, unless `END = 0` (then
  `START = 0` too), which holds every node to `N` (G6).
- A NoData node: class 255, never measured, never counted; its neighbours
  fill it as section 3 says (by reflection, exact on a plane). A void
  triangle (a NoData corner) is carved as today; no `SlopeScan`.
- A one-row or one-column grid, or a one-node-wide strip of data between
  NoData: no neighbour along one axis on either side, so that axis's part of
  the slope is lost (`gy = 0` or `gx = 0`) and the slope reads low: section
  3's one named departure (item 1b: 36.0 degrees for a 40-degree plane).
- A check point exactly on a grid line: filed in one cell by `CheckPoints`'
  floor rule; judged by that cell's corners only, as stored.
- `END = 90` is refused: no finite DEM gives a slope of 90 degrees, and the
  ramp would never reach N.
- `N > 0` is required (33 allows `N = 0` for lines): the weight is `1 / N`,
  and a zero tolerance would make it infinite. `--tolerance 0` without the
  slope is unchanged.

## 9. Tests `@tester` writes red first

Fixtures (`34-probes/fixture_figures.py`; figures in `34-probes/README.md`, item 4):

- **V1**, 65 by 65 nodes at 10 m: a flat floor (x below 200 m), a 40-degree
  wall to x = 440 m, a plateau, plus 3 m bumps `3 sin(2 pi y / 80) sin(2 pi x
  / 130)`. 35.4 % of its nodes are at 30 degrees or more (classes); with
  `N = 2`, `F = 10` and a step at 30, 1 496 of 4 225 nodes are held to 2 m.
  rasputin's uniform meshes from the node rectangle's outline: 28 triangles
  at 10 m, 796 at 2 m. Simulated: per node 408 triangles, uniform 2 m 1 142.
- **V2**, the same valley on the CLI fixtures' `micro_tiff` grid, 33 rows by
  41 columns, 10 m by 5 m: 48.9 % at 30 degrees or more; 661 of 1 353 nodes
  held to 2 m by the step at 30; rasputin uniform 6 triangles at 10 m, 112 at
  2 m.

- **V1c**, V1 with its south-east corner cut by an off-node domain edge
  from (640, -300) to (150, -640) m off the north-west node, across the wall;
  a constraint, so the feet act on it. rasputin master at `--tolerance 2`
  (the mesh test 4 (b) must equal): 339 vertices, 10 feet; 3 DEM-node
  vertices on the wall lie between `eps(2)` and `eps(10)` of the edge
  (`34-probes/mutant_fixtures.py`, item 4b), so under M6, whose epsilon is
  `eps(10)` there, each would be offered to the foot search as it went in,
  which the correct run does not do, and a foot taken changes the mesh.
- **The check points of test 6**, 4 per cell of V1 at dyadic offsets, z the
  bilinear surface plus noise of up to 4 m, drawn as
  `scattered(65, 4, 7)` draws them (`tests/python/test_core_refine_points.py@4d4cd740:109-120`):
  129 cells have corners on both sides of 30 degrees; 253 points in them are
  nearest a corner below 30 degrees, which M3 would hold to 10 m instead of 2, and 136 of those
  carry noise above 2 m (item 4b). Whether one of them ends between 2 and
  10 m under M3 depends on the mesh; hence M3's escape in the mutant list.
- **V1n**, V1 with one NoData node on the wall, row 32, column 30 (x =
  300 m). Its four cells (rows 31-32, columns 29-30) each have that NoData
  corner, and their largest valid corner classes are 83, 86, 83 and 81 (40.5
  to 43 degrees); `scattered(65, 4, 7)` puts 4 check points in each, 8 of the
  16 with noise above 2 m. Over the node rectangle 4 224 nodes are valid and
  1 495 are held to 2 m by the step at 30 (item 4b).
- **V1's line** (tests 5 and 7, with lines): a segment from (100, -50) to
  (600, -600) m off V1's north-west node, `N = 1`, ramp 0 to 200 m, `F = 10`,
  the slope's step at 30 with `N = 2`: 2 871 of 4 225 nodes are within 200 m
  of it; the line holds 1 752 nodes tighter than the slope does, the slope
  1 341 tighter than the line (item 4b), so both drivers decide somewhere.

C++ (`tests/cpp/unit`, `tests/cpp/property`):

1. **Today's mesh, bit for bit**: increment 33's test 1 stands; the new
   headers change no default instantiation. On V1 with constraints and feet
   on, `refine(..., options)` and `refine(..., options, UniformTolerance{t})`
   equal, for all three entry points.
2. **Steepness** against the test's own Horn in `long double`: every class on
   V1 and V2 is the smallest half degree at or above the test's slope (equal,
   or one class above where the test's slope lies within 1e-9 degrees of a
   class boundary: never lower, and at a boundary within rounding the higher
   class, pinned); **P1**, a plane at 37.3 degrees facing 30 degrees off the
   row on V2's 10 by 5 m grid, 21 by 21 nodes, with one NoData node at its
   centre: every valid node is class 75, the border nodes, the four corners
   and the eight nodes around the NoData one included (item 1b measures
   37.300 degrees at every valid node; round 1's rule gave 12.9 degrees at a
   corner of a 40-degree plane facing 30 degrees); a NoData node (class 255); a one-row grid,
   the named departure (36.0 degrees for a 40-degree plane facing 30
   degrees, item 1b);
   `dx != dy` (V2: a mutant using `dx` for both axes, M4, changes classes on
   the wall); rounding up (M5); 1 and 8 threads equal.
3. **The ramp**: `t` at 0, START, the middle, END, 89.5; the step; never
   increasing with slope over classes 0-180; entries 181-255 equal `F`.
4. **Equalities, bit for bit, on all three entry points, on V1c (V1 cut by
   the diagonal edge, above; constraints and feet on)**: (a)
   `Sloped<UniformTolerance{10}>` with `N = F = 10` equals
   `UniformTolerance{10}`; (b) `START = END = 0` with `N = 2`, `F = 10`
   equals `UniformTolerance{2}` (the slope must reach the feet's epsilon
   and `refine_points`' comparison, M6, M7); (c) `Sloped<LineTolerance>` with
   `N = F` equals `LineTolerance` alone.
5. **The guarantee, oracle independent of the code**: after refine on V1 and
   on V1c, every valid DEM node's `|z - plane|` is at most the ramp
   at the test's own class (test 2's oracle); with a line too (V1's line,
   above), at most the smaller of that and 33's ramp at the brute-force
   distance.
6. **Check points and the strip**: on 33's resampled-path setup (a check-point
   store over V1's geometry, 4 points per cell at dyadic offsets, z the
   bilinear surface plus noise of up to 4 m, as
   `tests/python/test_core_refine_points.py`'s `scattered(65, 4, 7)`), with
   `N = 2`, `F = 10` and a step at 30 degrees, every point's
   error is at most `t` of the largest class of its cell's valid corners
   (0 when none is), computed in the test (M3); and the edge strip along the
   diagonal edge: every strip point within its cell's `t`, on V1c.
   (b) **V1n** (*Fixtures* above), the same store and flags: the 16 check
   points in the four cells around the NoData node, each cell with that
   NoData corner and valid corners at 37 degrees or more, end within
   `N = 2`, whatever triangle holds them (none has the NoData node as a
   corner: refine never inserts one); and a `cell_class` unit test on V1n's
   cell (31, 30) gives 86, not 255 (M9). (c) On V1n, `slope_nodes` reports
   4 224 valid nodes and 1 495 tightened: the NoData node is in neither.
7. **Threads**: 1 and 8 threads give the same mesh, with the slope alone and
   with lines (the steepness pass and the scan on the TSan job's list).
8. **Exactness of `over`** (M1): a scan of one triangle whose single inner
   node has an error one unit in the last place above its `allowed(class)`,
   chosen by the test so that `error * weight` rounds to exactly 1; the
   triangle must split.
9. `SlopeTolerance::make`'s refusals, one each ("near", "far", "start",
   "end"); the binding's other-geometry refusal.

Python (`tests/python`):

10. The binding: `SlopeTolerance(view, ...)`'s `histogram()` sums to the
    node count and matches the test's own classes on V1; `slope=None` is
    today's call.
11. CLI refusals of section 5, each in its words.
12. CLI on V2: the record's four fields; `max_error_slope_share_of_tolerance <=
    1 + 1e-12`; more vertices on the wall (x from 200 to 440 m) than the
    `--tolerance F` mesh has there; `N = F` gives the `--tolerance F` mesh
    and the stderr line; without the flag the record has none of the four
    fields.
13. The resampled path (`--out-crs`): the slope is built on the target grid
    and passed to the final check (its record value is the final check's).

**Invariant-critical suite: tests 2, 4, 5, 6 (with (b) and (c)), 8.** Mutation targets the kill
record must cover:

- M1 converge on `ratio <= 1` instead of the exact `over` (test 8);
- M2 rank by `error * weight` alone, no second key (test 4 (a), on a fixture
  where two different errors give one product; `@tester` finds or builds one,
  or reports M2 as equivalent on the fixtures with the evidence);
- M3 a check point's class from its nearest node instead of its cell's
  largest valid corner (test 6; layout and figures under *Fixtures* above; if no point ends between
  N and F under the mutant, `@tester` reports M3 as equivalent on the fixture
  with that evidence);
- M4 Horn with `dx` for both axes (test 2, V2);
- M5 the class rounded down (test 2);
- M6 the foot's epsilon with the triangle part's tolerance only (test 4 (b)
  on V1c; figures under *Fixtures* above; if the mutant's mesh is
  equal, `@tester` reports M6 as equivalent on the fixture with that
  evidence);
- M7 `point_loop` ignores the slope (test 4 (b) on `refine_points`, test 6);
- M9 `cell_class` takes the plain largest of the four corners, 255 included
  (test 6 (b): the unit case always; the 16 points, 8 of them with noise
  above 2 m, if one ends between N and F under the mutant);
- M8 the lazy test asks the triangle part even when `over` (no output change:
  an efficiency mutant, killed only by counting `at()` calls; `@tester` adds a
  counting test policy, or reports it as equivalent in output).

### 9.1 Pins after the red step

`@tester`'s red commit `3ce92a7e` settled what this design left open. Listed
from that commit's message and from the tests at `7d661ae9`; recorded, not
ruled.

1. *`RefineOutcome::slope_nodes`* has two members, `valid` and `tightened`
   (`tests/cpp/property/prop_refinement_slope_tolerance.cpp`, header,
   "CHOSEN HERE"); test 6 (c) also pins that only nodes inside the mesh are
   counted (V1c's pentagon).
2. *Python types:* `SlopeTolerance.histogram()` returns 256 ints, and
   `RefineOutcome.max_slope_share` is a float property
   (`tests/python/test_core_slope_tolerance.py`, header).
3. *`slope=`* is keyword-only and last, after 33's `field=`, on `refine`,
   `refine_strip` and `refine_points`: the stub's keyword pins in
   `tests/python/test_core_refine_points.py`, `test_core_edge_strip.py` and
   `test_cli_tolerance_near.py` gain it.
4. *A class boundary:* z = x on a 10 m grid (45 degrees, between classes 90
   and 91) is class 91 at every node, the border included
   (`tests/cpp/unit/test_raster_steepness.cpp`, the case "a slope exactly on a
   class boundary takes the higher class (pinned)").
5. *Section 3's two departures as test values:* the one-row grid gives class
   73 (36.005 degrees); NoData above and below a node gives 46 (22.760
   degrees).
6. *`cell_class`:* V1n's cell (31, 30) is 86, not 255; a cell with no valid
   corner is 0; a position past the last cell clamps to it
   (`tests/cpp/unit/test_refinement_slope_tolerance.cpp`).
7. *`make`'s refusals (test 9):* each sentence contains the name of the value
   broken ("near", "far", "start" or "end"); the rest of the words are not
   pinned. The edges of the bounds are accepted with `why` empty.
8. *Test 8's fixture:* an error one unit in the last place above
   `allowed(class)`, chosen so that `error * (1 / N)` is exactly 1.
9. *M8:* a counting policy, `CountingAsks`, counts the triangle part's `at()`
   calls for triangles whose nodes' slope already splits them.
10. *Test 6 in Python* runs on `scattered(65, 4, 7)`'s own draw, with a guard
    that the draw gives 129 mixed cells, 253 points and 136 with noise above
    2 m.
11. *Test 11's order:* for four triples that break two bounds, the first in
    section 5's order is named and the later one is not said.
12. *Test 12's record and output:* `tolerance_slope_m` is a JSON float;
    `tolerance_slope_deg` reads "30 to 30" for a step and "25 to 35" for a
    ramp; `slope_nodes_tightened` is matched as "661 of 1353 DEM nodes (<p>
    %)", with p read as a number within 0.05 of the count's share ("CHOSEN
    HERE"); the share is a float above 0 and at most `1 + 1e-12`; the summary
    contains "Nodes on steep ground are within the tolerance their slope
    allows (largest error"; N = F builds no slope and prints
    "--tolerance-slope: N equals --tolerance, so slope changes nothing";
    `--help` contains "--tolerance-slope" and "steep".
13. *Test 13:* the record's share equals the final check's last
    `max_slope_share`, and the `--stats` report has a `slope` phase.
14. *The quick-check pin:* `test_the_shipped_cases_are_sections_3_table` in
    `tests/python/test_bench_quick.py` gains `("romsdal-slope", [0], 1, 3)`,
    `--tolerance 10`, `--tolerance-slope 2 25 35`, the tile, and the input
    key `$RASPUTIN_DATA/DTM10_UTM33_20260925`. Departure from section 10,
    which put the case in `docs/benchmarks/quick/cases.toml` in the same
    commit: that file was outside `@tester`'s write limit, so the pin stays
    red until `@perf`'s acceptance commit adds the case. At `7d661ae9` it
    fails (`pytest tests/python/test_bench_quick.py -k shipped_cases`).

### 9.2 Kill record

Summary (not part of the record): `@tester`'s mutation round, commit
`7d661ae9`, planted M1-M9 in a scratch copy of `f68a9440` and killed all
nine; M2 is killed by a fixture new in that commit, two adjacent-double
errors with one product, where only the second key picks today's node. The
table below is the record, pasted as `@tester` gave it. The commit message of
`7d661ae9` lists test 3 among M5's kills; this table does not, and the table
is the record. The commit is left as it is (a default taken while Ola was
away).

| mutant | planted in | change | killed by |
|---|---|---|---|
| M1 | refinement/scan.hpp, SlopeRank::rank | `over` from `!(err * weight <= 1)` instead of `err > allowed` | test 8 |
| M2 | refinement/scan.hpp, SlopeRank::rank | rank by `q > ratio` alone, no second key | test 4 (a), the M2 fixture (new here) |
| M3 | refinement/refine_points.hpp, scan_points | a stored point's class from its nearest node's row() instead of cell_class | test 6 (a) check points; test 6 (b) V1n |
| M4 | raster/steepness.hpp | dy8 = 8 dx | test 2: V1 and V2 classes, P1, the M4 case |
| M5 | raster/steepness.hpp | class = count of bounds, no +1 (rounded down) | test 2 (8 cases); test 6 (b) cell_class unit cases (2) and the histogram case; tests 5, 6 (a), 6 (b), 6 (c) |
| M6 | refinement/refine.hpp, the foot's tol | min with the triangle part's highest() instead of allowed(node's class) | test 4 (b) |
| M7 | refinement/refine_points.hpp, point_loop | `over` forced false: points compared with the triangle part only | test 4 (b); test 6 (a) points and strip; test 6 (b) |
| M8 | refinement/refine.hpp, the scan's slot | allowed_at (asks the triangle part) evaluated even when `over`; output unchanged | the M8 counting policy |
| M9 | refinement/slope_tolerance.hpp, cell_class | the plain largest corner, 255 included | test 6 (b) cell_class unit cases; test 6 (b) V1n's 16 points |

### 9.3 Pins after the green step

`@developer`'s green commit `952bd1ee` settled the rest. Listed from that
commit's message and from the code at `7d661ae9`; recorded, not ruled.

1. *The classes, from a table with a relative lowering.*
   `raster::steepness_bounds()` holds tan² of 0.5 to 89.5 degrees, each
   multiplied by `1 - 1e-12` (about 1e-11 degrees); a node's class is 1 plus
   the number of bounds at or below gx² + gy², and 0 when that sum is 0; no
   atan per node. So a slope on a boundary, to a relative 1e-12, takes the
   higher class (pin 4 above). The constant is relative, so it holds at any
   grid spacing. A NaN height counts as missing, as the NoData value does
   (`include/terrain/raster/steepness.hpp`).
2. *`make`'s checks and words, in this order:* each value finite (near, far,
   start, end: "<name> must be finite"), then "near must be above 0", "near
   must be at most far", "start must be >= 0", "start must be at most end",
   "end must be below 90 degrees". The binding raises
   `ValueError("SlopeTolerance: " + that sentence)`.
3. *The tables:* 256 entries each; 0-180 from the ramp at `c / 2` degrees,
   181-255 equal to F, weight `1 / allowed`.
4. *`SlopeRank`, shared by two scans.* `{over, ratio, error}` and its
   `rank()` are the base of the scan's `SlopeScan` (`scan.hpp`) and of
   `PointScan` (`strip_scan.hpp`). The design's `SlopeScan` had no `error`;
   the second key needs it. `PointScan`'s own `error` moved into the base, so
   `PointScan` grows by `over` and `ratio`: 104 to 120 bytes, by `sizeof` at
   `da3d1d63` and `7d661ae9` (Apple clang 21, arm64). That is on the default
   path too: `refine_points` and `refine_strip` keep one `PointScan` per
   triangle slot with or without a slope. `ScanResult`, which section 4.3
   kept as it was, is unchanged.
5. *Compile-time switches:* `detail::varies<P>` judges
   `SlopeParts<P>::Triangle`; `slope_of(policy)` gives the slope or a null
   `const NoSlope*`; `NoSlope` is the default slope type of `scan`,
   `scan_points` and `scan_strip`, so without a slope each compiles as
   before.
6. *`refine`:* when `over`, the slot takes the slope's node and its location,
   and `-inf` in `allowed`; the foot's epsilon takes the smaller of the
   triangle part's tolerance and `allowed` of the inserted node's class.
7. *`refine_points` and `refine_strip`:* a second `PointScan` per slot,
   `steep`. Stored and strip points are ranked at `cell_class` of their
   position, rescanned DEM nodes at their own class; `offer_steep` ORs `over`
   and keeps the strictly larger `(ratio, error)`, ties to the earlier set. A
   void triangle ignores the slope's `over` (`over = sloped &&
   !results[t].is_void && s->over`); when `over`, `results[t].take(*s)`
   inserts the slope's point and keeps the triangle's own `max_error`,
   `uncovered`, `on_frozen` and `frozen_error`.
8. *The count of tightened nodes* lives in `slope_tolerance.hpp`
   (`count_slope_nodes`, blocks of 4096 slots), not in `refine.hpp`. It runs
   at the end of `refine` and of `point_loop`, so the edge strip and the
   final check each count again over their final mesh. The record takes
   `refine`'s count on the DEM path (`inside, tight = out.slope_nodes`) and
   the final check's on the resampled path.
9. *The binding:* `Sloped<const LineTolerance&>` (`P` may be a const
   reference, so the field is not copied) and `Sloped<UniformTolerance>`,
   dispatched four ways; `RefineOutcome.slope_nodes` is a Python tuple
   `(valid, tightened)`; `SlopeTolerance(view, near, far, start, end,
   threads=0)` releases the GIL while it computes the classes; the stub's
   `histogram()` returns `list[int]`.
10. *The CLI:* `_tolerance_slope` runs right after `_tolerance_lines`, before
    the other refusals of `mesh`; refusal values are printed with `:g`. The
    flag sits in its own `--help` panel, "Tolerance on steep ground"; the
    green commit's reason: in the main table its wider metavar truncated
    `--tolerance-near-crs` at 80 columns (not re-run for this record). The
    slope is built in `_dem_mesh` in the `--stats` phase `slope`, with F =
    `--tolerance` and the default `threads` (0, all cores).
11. *The run record:* `tolerance_slope_deg` is `f"{start:g} to {end:g}"`;
    `slope_nodes_tightened` is `"<k> of <n> DEM nodes (<p> %)"`, p to one
    decimal, `n` floored at 1 in the division. `run_record._entry` now records
    any float as a number (before: only names ending `_m` or `_deg`), so the
    share is a JSON number; the rule reaches every float-valued entry,
    whatever its name. The summary: with the share at most 1, the sentence
    "Nodes on steep ground are within the tolerance their slope allows
    (largest error <s> % of it)." with s to 4 significant digits; above 1,
    no such sentence, and a separate line "Warning: the largest error on
    steep ground is <s> % of what its slope allows."
12. *N = F:* `_tolerance_slope` returns nothing after the stderr line, so no
    slope is built, no `slope` phase is timed, and the record has none of
    the four slope fields.
13. *TSan:* `test_raster_steepness` and `prop_refinement_slope_tolerance`
    added to both lists of the TSan job (`.github/workflows/main.yaml`), as
    the red commit asked.

## 10. Size, split point, speed

| part | counted lines (estimate) | basis |
|---|---|---|
| `raster/steepness.hpp`: Horn, the fill of missing neighbours, classes, parallel rows | 55-75 | `line_tolerance.hpp`'s ramp and `make` (33): about 40 of its 143 |
| `refinement/slope_tolerance.hpp`: ramp, tables, `make`, `cell_class`, `histogram`, `Sloped`, trait | 70-95 | 33's `line_tolerance.hpp`, 143, less its search |
| `scan.hpp`: `SlopeScan`, the ranking, the byte row | 40-60 | the node and off-node branches, about 50 lines today |
| `refine.hpp`: the slope vector, the split test, the node choice, the foot, `max_slope_share`, the count of tightened nodes inside the mesh | 45-65 | 33 added 59 to it, with `max_error_near` (not repeated here) |
| `refine_points.hpp`, `strip_scan.hpp`: the point classes, the choice across sets | 40-60 | 33 added 33 to `refine_points.hpp` |
| binding and `_core.pyi`: the class, `slope=` on three functions, four-way dispatch | 50-70 | 33: 54 |
| `tolerance_field.py`: `ToleranceSlope` | 8-12 | |
| `cli.py`: one flag, the refusals, building and timing, the record | 50-70 | 33's three flags and the field: 105; one flag and no file here |
| `edge_strip.py`, `final_check.py`, `run_record.py` (with the summary sentence) | 18-26 | 33: 14 |
| **total** | **376-533** | |

**As built**, `python3 tools/count_loc.py da3d1d63 7d661ae9` (the design's
last commit to the mutation round; `9ba38490`, the branch's base, gives the
same 512, since the design commits are prose only). The green commit alone,
`python3 tools/count_loc.py 3ce92a7e 952bd1ee`, also gives 512: the red
step, the scaffolding removal and the mutation round touch no counted file.

| part | estimate | built | why it differs |
|---|---|---|---|
| `raster/steepness.hpp` | 55-75 | 78 | 3 over |
| `refinement/slope_tolerance.hpp` | 70-95 | 160 | holds the count of tightened nodes (`SlopeNodes` and `count_slope_nodes`, 45 lines), which the estimate put in `refine.hpp`'s row; without it 115, 20 over (`make` alone, with its checks and words, is 24) |
| `scan.hpp` | 40-60 | 30 | under |
| `refine.hpp` | 45-65 | 37 | under: the count moved to `slope_tolerance.hpp`; with it, 82 |
| `refine_points.hpp`, `strip_scan.hpp` | 40-60 | 38 + 18 = 56 | within |
| `bindings/core.cpp` and `_core.pyi` | 50-70 | 38 + 19 = 57 | within |
| `tolerance_field.py` | 8-12 | 4 | under |
| `cli.py` | 50-70 | 63 | within |
| `edge_strip.py`, `final_check.py`, `run_record.py` | 18-26 | 3 + 3 + 21 = 27 | 1 over: the summary sentence and its warning line, and the `_entry` change (9.3, pin 11) |
| **total** | **376-533** | **512** | within; under 700, one PR |

Under 700 in one PR. **Split point** if the red suite pushes the estimate past
600: PR A the C++ headers, the entry points, the binding and stub (nothing a
user sees changes); PR B the model, the CLI and the record.

**Cost of every added step** (the quick check's catchments are measured at
default flags in item 3: Numedalslagen 11.12 s, Lagan 21.56 s; they gain
nothing, G1):

| step | default flags | with `--tolerance-slope` | basis |
|---|---|---|---|
| classes (`steepness`) | none: not built | one pass over the grid: at most 24 ns per node on one thread, twice NumPy's float64 plain Horn on the tile (12.2 ns, item 5), which is the arithmetic of every interior node; the fill of section 3 runs only on the outer ring and beside NoData (NumPy applies it to every node: 40.4 ns, item 5); tile `6901_3` at most 0.61 s of CPU, Numedalslagen's canvas at most 6.7 s of CPU, Lagan's target grid at most 0.42 s; about a tenth of that in wall time on 10 threads | estimate, not measured (no C++ built) |
| the scan | none: compiled as today | per node one byte read, two table reads, a multiply, two comparisons; the scan is 24 % of refine at uniform 2 m on the tile (0.891 of 3.744 s); if the scan slows by 30 to 60 %, refine per output triangle slows by 7 to 15 % | estimate, not measured |
| the resampled path (`--out-crs`) | none | the classes of the target grid (Lagan's 1.74e7 nodes: the steepness row above), and per check point one `cell_class` (four byte reads beside the 16-byte point); the final check's scan was 0.34 s of Lagan's 21.56 s (item 3), so even doubled it adds about 1.6 % | estimate, not measured |
| the count of tightened nodes | none | one parallel pass over the final triangles' node sets, after the loop: at most two byte reads per node (4.5), against the scan's four-byte height and plane per node; so at most one scan's worth (0.89 s at uniform 2 m on the tile, item 3, on 10 threads) | estimate |
| the triangle part's query | none | none with the slope alone (`UniformTolerance`); with lines fewer than 33's, since `over` skips them | section 4.3 |
| memory | 16 bytes more per triangle slot in the final check and the edge strip (`PointScan` 104 to 120 bytes) | one byte per node; per triangle slot 40 bytes in `refine` and a second 120-byte `PointScan` in the final check and the edge strip (4.6) | `sizeof` at `7d661ae9` |

**Speed judgment.** This touches refine and scan code, so `@perf`'s
acceptance applies (`docs/increments/README.md`): `tools/bench.py`'s 1 m
benchmark and thread sweep at default flags show no change beyond the bench's
band, with the mesh's SHA-256 unchanged (G1). The quick check gains the case
`romsdal-slope`: `mesh --dem $RASPUTIN_DATA/DTM10_UTM33_20260925/6901_3_10m_z33.tif
--tolerance 10 --tolerance-slope 2 25 35` (threads `[0]`, runs 3, warmup 1),
no baseline at first. Its `inputs` entry is the directory,
`["$RASPUTIN_DATA/DTM10_UTM33_20260925"]`, the key `numedalslagen` and
`geilo-al-ramp` already use, not the tile file: `_mismatch` in
`tools/bench_quick.py@4d4cd740:106-117` refuses a baseline whose input keys
differ, so a new key would make every case's baseline stop matching. The
tile lies in that directory, so the directory's entry covers it. Adding it breaks the pin on the shipped cases,
`test_the_shipped_cases_are_sections_3_table`
(`tests/python/test_bench_quick.py@1a754deb:399-411`), which lists every
case's name, threads, warm-up, runs and `--tolerance`: `@tester` adds
`("romsdal-slope", [0], 1, 3)` and `"romsdal-slope": 10.0` to it in the red
step, and the case to `docs/benchmarks/quick/cases.toml` in the same commit,
so the pin and the file never disagree on a commit. `@perf` records the
case's first baseline. Accepted when the default runs are unchanged and the
case's refine time per output triangle is at most 1.25 times the uniform
`--tolerance 2` run's on the same tile (0.633 µs: 3.744 s for 5 916 344, one
run, item 3; `@perf` re-measures it beside the slope run), and its `slope`
phase is under 0.5 s of wall time (Q7). "Refine time" is `--stats`' `refine`
phase; the classes are built before it, in their own phase.

## 11. Not in scope

A slope-dependent largest triangle size (Ola ruled "vtol means finer");
curvature as a driver; smoothing the DEM before the slope (a later flag for
DTM1, section 3); the slope from the original DEM on the resampled path (Q4);
the seam pass's slope with increment 23's pieces (4.6); a per-node rule for
33's lines (a distance query per node is what 33 avoided); a slope output on
the mesh.

## 12. Questions for Ola

Each with the default this design is written on.

- **Q1. The flag.** `--tolerance-slope N START END`: the general tolerance up
  to START degrees, N from END degrees, linear between; a step when START =
  END. The other shape would be `--tolerance-steep N ANGLE`, a step only.
  *Default: the three-number flag; your "steeper than a certain angle" is
  `--tolerance-slope N A A`.*
- **Q2. Each node held to its own slope's tolerance**, not each triangle to
  its steepest node's. Measured in the simulation: the triangle rule needs
  15 to 23 % more triangles in Romsdalen and 41 to 76 % more around Geilo,
  for the same guarantee. *Default: per node.*
- **Q3. Which slope.** Horn's 3 by 3 slope on the DEM's own grid, no
  smoothing (on DTM10 the choice of formula and 1 m of added noise move the
  steep share by under 1 percentage point at 30 degrees). *Default: Horn, no
  smoothing; a smoothing option later if DTM1 shows noise.*
- **Q4. With `--out-crs`**, the slope comes from the grid rasputin meshes on
  (the resampled one), not from the original DEM. That grid is square, at the
  original's north-south spacing rounded to whole metres
  (`src_python/tin_engine/target_grid.py@1a754deb:91-102`); a geographic
  DEM's east-west spacing is finer (by the cosine of the latitude), so the
  resampled grid sees less of the steepest detail and its slope reads a
  little lower. Measured on GLO-30 in the Wetterstein and Karwendel (section
  3, *The resampled path*): 41.0 % of the nodes at 30 degrees or more against
  the source's 42.2 %, the median 26.7 against 27.2 degrees. *Default: the
  resampled grid, the one the check points are filed in; slopes read about
  half a degree low there.*
- **Q5. How tight is useful on steep ground?** The DEM's own error grows with
  slope (Hutchinson et al., section 2). A steep tolerance below the DEM's
  error there follows the DEM's noise. *Default: no floor is enforced; the
  example uses 2 m on DTM10.*
- **Q6. The demonstration and quick-check case**: the Romsdalen tile with
  `--tolerance 10 --tolerance-slope 2 25 35`. *Default: that case.*
- **Q7. How much slower per triangle may a slope run be?** Estimated 7 to 15
  % slower per output triangle than a uniform run. *Default, taken while you
  were out: accept up to 1.25 times, and the slope pass under 0.5 s on the
  tile; above that the scan's slope work is made cheaper before merging.*

## 13. ROADMAP

Row 34, added with this design; it carries the Status line's state.

Follow-up options from code review round 1's suggestions, not in this PR:

- *S2:* `@perf`'s acceptance times both point-scan paths (the final check and
  the edge strip) at default flags, since `PointScan` grew by 16 bytes on
  every run (4.6). If that shows a cost, the slope's winner (`steep`) gets a
  type of its own, so the default `PointScan` returns to 104 bytes.
- *S3:* the binding does not check that the slope's F equals the `tolerance`
  passed beside it; the CLI builds both from `--tolerance`. A check in
  `with_policy` would refuse a mismatch from Python callers.

## Review

Increment 34 design review round 1, 2026-10-09, @reviewer (`9ba38490..1a754deb`, prose and probes only, 0 counted lines): CHANGES REQUESTED. Five blockers, in one pass: (B1) Horn's slope is understated by about half at the grid border and next to NoData under the "missing neighbour takes the node's own z" rule (a 40° plane reads 22.8° at a border node, 18.4° at a corner, below the 30° step, so those nodes get F); the probes use `np.pad(mode="edge")`, not the stated rule; fix with a rule exact on planes (e.g. 2·z(node) − z(opposite) when that exists) or name the departure, make the probes compute the stated rule, and add a steep border node and a NoData-adjacent node to test 2. (B2) Q4's reason is false: the target grid is square at the source's north-south spacing rounded to whole metres (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/src_python/tin_engine/target_grid.py@1a754deb:91-95), so its slope can read lower than the source's; the default stands, the sentence and its effect need stating. (B3) M3 and M6 kills unbacked: M6 dies in test 4(b) only if a constraint crosses V1's wall with a foot search reaching it (foot epsilon `clamp(tol / G, dx/100, dx/2)`, /Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/include/terrain/refinement/refine.hpp@1a754deb:240-256); M3 only if check points fall in mixed-class cells ending between N and F; state the layouts or give the equivalence escape. (B4) "A sum of ridges" (tests 4, 5) unnamed and without figures. (B5) the quick-check case `romsdal-slope` breaks the shipped-case pin (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/tests/python/test_bench_quick.py@1a754deb:399-411); name the change and its owner. Sound: per-node Horn slope at the DEM's spacing, 1 byte per node; `Sloped<P>` and slope = --tolerance giving today's mesh bit for bit; the split-node rule; the default path byte-identical; Q1-Q3, Q5, Q6 defaults; quotations verbatim where read, and unread sources stated. Probes re-run: fixture_figures.py and estimate.py match README items 4 and 5. Suggestions S1-S10: imagecodecs in the probe environment; resampled-path cost; Q7 about 1.25 times; pin rounding to the higher class; `slope_nodes_tightened` scope; summary-line wording; section 6 header; foot path asks `at()` on −∞ slots; G4 excludes frozen-edge check points; an estimate for slope with lines on geilo-al-ramp.

Increment 34 design review round 2, 2026-10-09, @reviewer (`75711c70..4d4cd740`, prose and probes only, 0 counted lines): CHANGES REQUESTED. Round 1's five blockers are answered; border_check.py, resample_check.py, mutant_fixtures.py and lines_sim.py re-run word for word against the probes README (items 1b, 1c, 4b, 6). Two one-line fixes: (R1) the size table's total says 380-545 but its rows add to 376-533, and ROADMAP row 34 repeats the wrong range; (R2) test 6's check-point figures (129 mixed cells, 246 points, 136 with noise above 2 m) come from the probe's own random draw, not from `scattered(65, 4, seed 7)`'s draw order (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/tests/python/test_core_refine_points.py@4d4cd740:109-120), which gives 129, 253 and 136; compute from that order or say the draw differs. Met: B1 (fill exact on planes, the one-row departure stated at 36.005°), B2 (Wetterstein 41.0 vs 42.2 %, median 26.7 vs 27.2°), B3 (V1c's three wall vertices between eps(2) ≈ 2.4 m and the 5 m cap; M3's layout and fallback), B4 (2 871 / 1 752 / 1 341 reproduced), B5 (the `romsdal-slope` case and its pin in the same red commit); S3's 1.25 times consistent with the 7-15 % estimate; S5's ownership rule sound. Suggestions: bound the node count after the loop per node and say whether it runs in parallel (it also needs each node's NoData flag, which class 0 does not give); name `romsdal-slope`'s input as the directory key so the baseline still matches (`_mismatch` in /Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/tools/bench_quick.py@4d4cd740); narrow section 3's "loses that axis's part" sentence; spell out N = 2, F = 10 and the step at 30° in tests 4 (b) and 6; the case table in docs/increments/perf-quick-check.md section 3 already lacks `geilo-al-ramp`.

Increment 34 design review round 3, 2026-10-09, @reviewer (`e8d14626..838650f8`, prose and probes only, 0 counted lines): CHANGES REQUESTED. R1 and R2 fixed (376-533 in the table and ROADMAP row 34; mutant_fixtures.py re-run prints 129 / 253 / 136); suggestions taken correctly (border_check.py 22.760° = atan(tan 40° / 2); the parallel node count; the directory input key per `_mismatch` at /Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/tools/bench_quick.py@4d4cd740:106-117). Class 255 is sound against G1, G2 and G6 (the scan zeroes a NoData node's error, /Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/include/terrain/refinement/scan.hpp@838650f8:194-197) but leaves two holes: (H1) four places still say "the largest of the 4 corners" while 4.5 says "the largest valid corner" (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/docs/increments/34-slope-tolerance.md@838650f8:247-248, :394, :766, :800); with 255 a class value, the plain largest would hold points beside NoData to F and break G4; each needs "valid" (0 when no corner is valid); also :342's "181-entry table" and the `allowed(c)` comment without the 255 exception. (H2) The valid-corner rule has no test and no mutation target: add a test 6 case or a `cell_class` unit test (one NoData corner and one corner at 30° or more, the point inside an all-valid triangle, held to N = 2), named mutant M9 "`cell_class` takes the plain largest of the four corners, 255 included". Suggestions: say what `make` puts in entries 181-254 or that they are never read; pin in a test that a NoData node is absent from both `slope_nodes_tightened` numbers.

Increment 34 design review round 4, 2026-10-09, @reviewer (`12c30aa4..e0f31c8a`, prose and probes only, 0 counted lines): CHANGES REQUESTED, one figure: (F1) the V1n paragraph gave classes 83, 86, 83 and 81 as "41.5 to 43 degrees"; class 81 is 40.5°, so 40.5 to 43 (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/docs/increments/34-slope-tolerance.md@e0f31c8a:729-730); fixed by the main session in this recording commit. Otherwise round 3's holes are answered: every corner rule says "largest valid corner, 0 when none is valid" (lines 248-250, 396, 430, 550, 672, 778, 819); the data-flow figure, the `allowed(c)` comment and section 4.5 agree on 0-180 from the ramp and 181-255 as F; test 3 covers both; mutant_fixtures.py `m9` re-run matches README item 4b; test 6 (b)'s unit case kills M9 on every run; test 6 (c) excludes the NoData node from both counts; the 376-533 estimate holds. Suggestions for @tester's brief: test 6 (c)'s 4 224 / 1 495 are a count over DEM nodes, reported by `refine`'s `RefineOutcome`; M9 belongs after M8 in the list.

Increment 34 design review round 5, 2026-10-09, @reviewer (`e0f31c8a..cffb87bb`, prose only, 0 counted lines): APPROVED. F1 fixed: the V1n paragraph reads "83, 86, 83 and 81 (40.5 to 43 degrees)" (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/docs/increments/34-slope-tolerance.md@cffb87bb:729-730); the commit otherwise adds only the round 4 record. The 376-533 estimate stands. The stale Status line and ROADMAP row 34 are updated in this recording commit. Round 4's two suggestions go to @tester's brief.

Increment 34 code review round 1, 2026-10-09, @reviewer (`da3d1d63..61238e7f`, 512 counted lines; not pushed, no CI): CHANGES REQUESTED, four small blockers, no correctness defect. (B1) leftover red-step indirection: /Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/tests/python/test_core_slope_tolerance.py@61238e7f:49-56 reaches `_core.SlopeTolerance` through an unused `# type: ignore[attr-defined]` and `Factory`. (B2) memory figures stale: sections 4.6 and 10 say about 32 bytes per slot (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/docs/increments/34-slope-tolerance.md@61238e7f:590-597, :1055); as built `SlopeScan` is 40 bytes, with the flag `refine_points` and `refine_strip` keep a second 120-byte `PointScan` per slot, and `PointScan` grew 104 → 120 bytes on every run. (B3) /Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/tests/cpp/unit/test_refinement_slope_tolerance.cpp@61238e7f:19 names tests 3, 6 and 9 as invariant-critical; section 9 names 2, 4, 5, 6, 8. (B4) six line citations shifted by this branch: 15f-edge-strip.md:87 and :358, 23-basin-scale.md:3263, 27-node-sampling.md:148 and :172, docs/retrospectives/2026-10-07-night-20c-1.md:477; pin each to `@da3d1d63`. Holds: steepness, tables, `cell_class`, `Sloped<P>`, the split rule, the serial phase and the parallel count match the design (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/include/terrain/raster/steepness.hpp@61238e7f); default path byte-identical, re-run by the reviewer on tile 6901_3 at `da3d1d63` and `61238e7f`, DEM path SHA-256 `48ebbf92…` and resampled path `d10bff3a…` on both; pins match the code, `run_record._entry` changes no existing entry; TSan lists (/Users/skavhaug/projects/rasputin/.claude/worktrees/slope-tol/.github/workflows/main.yaml@61238e7f:161,177); kill record M1-M9 with each killing test present; 512 inside 376-533. Suggestions: (S1) a share of 1 + 1e-12 from rounding alone prints a false Warning; (S2) @perf covers both point-scan paths at default flags, optionally `steep` gets its own type; (S3) the binding does not check the slope's F equals `tolerance`.
