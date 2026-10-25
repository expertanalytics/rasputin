# Increment 105: the cell table — hydrological cells on a coarse mesh

**Status: designed, not built.** By `@architect`, 2026-10-25, on master
`074e6322`, worktree `shyft-cells`, from Ola's two sentences of 2026-10-25:
"For a coarse triangulation, we should provide the typical shyft-cell data.
Like lake fraction, forest fraction, etc." and "So we can do this already?
So, not use the land types as constraints, but simply partition the cell into
its types?" Number 105: 104 is the highest under `docs/increments/` over
every ref (`git for-each-ref` + `git ls-tree`, 2026-10-25).

**Size: bigger than the one-page threshold.** About 305 net production
lines in three parts (section 9), each its own PR under 150 lines: **105a**
the table (works today for a class raster, so for São Francisco with
MapBiomas), **105c** fractions from polygons (what Kultsjön and Norway need:
their land cover is vector), **105b** reservoirs. The probe (section 8) was
run before any of it is built.

## 1. The answer to Ola's question

Yes. `rasputin mesh --land-cover <class raster>` already partitions every
triangle exactly into its classes without making the classes constraints
(`src_python/tin_engine/land_cover_fractions.py`, increment 26). What is
missing: the five Shyft categories summed from the classes, the per-cell
geometry Shyft wants, and a GIS table holding them. Today the fractions reach
only `.vtk` (three packed arrays) and `.ply` (face lists); the GeoPackage
`triangles` layer has one class per triangle and no fractions.

## 2. Prior art: legacy and literature

**Shyft.** Its cell is `geo_personality` (`cpp/shyft/hydrology/cell/geo_personality.h`
on `gitlab.com/shyft-os/shyft`, master, read 2026-10-25; `geo_cell_data.h`
only aliases it: `using geo_cell_data = hydrology::cell::geo_personality;`).
`land_type_fractions` stores glacier, lake, reservoir and forest;
`unspecified()` is `return 1.0 - glacier_ - lake_ - reservoir_ - forest_;`,
`snow_storage()` is `1.0 - lake_ - reservoir_`; `set_fractions` allows a sum
up to `1.0 + 1.0e-3` ("Allow small round off errors") and otherwise throws
"LandTypeFractions:: must be >=0.0 and sum <= 1.0". The grid constructor is
`geo_personality(mid_point, area = default_area_m2, catchment_id = -1,
radiation_slope_factor = default_radiation_slope_factor, land_type_fractions,
routing_info)` with `double constexpr default_radiation_slope_factor = 0.9;`.
**The TIN constructor** takes three `geo_point` vertices, computes
`mid_point_ = calc_mid_point(); area_m2 = calc_tin_horizontal_projected_area();`
and sets `radiation_slope_factor_ = 1.0; // no factor required if we have TIN`;
`slope()` and `aspect()` come from `normal_up(v1, v2, v3)`. `cell/model.h`:
the factor is "used by the radiation interpolation to adjust the local
radiation with respect to idw sources". The Python tutorial
(shyft.readthedocs.io, *The Shyft api*) builds `LandTypeFractions(glacier,
lake, reservoir, forest, unspecified)` from a NetCDF with `forest-fraction`,
`reservoir-fraction`, `lake-fraction`, `glacier-fraction`, `x, y, z, area,
catchment_id`, with `rad_fx = 0.9`. **So the brief's premise is off:** Shyft
has no slope formula for the factor. It is a tuning multiplier on
interpolated radiation, 0.9 for a 1 km grid cell and **1.0 for a TIN cell**,
whose slope and aspect Shyft takes from the vertices. This design writes
slope and aspect, and no factor column (section 4.3). Shyft and rasputin
together: Silantyeva et al., *Shyft and Rasputin*, `README.md` *Publications*.

**Which classes are forest, lake, glacier** is not Shyft's: it takes the
fractions as given. The groupings below follow each source's own
nomenclature: CORINE level 1 to 3 (3.1 Forests = 311 to 313; 5 Water bodies
= 511 to 523; 335 Glaciers and perpetual snow), the AR50 `arealtype` codes
(30 Forest, 70 Snow and glacier, 81 Freshwater, 82 Ocean; labels in
`ar50.toml`), and MapBiomas Collection 11's legend
(`legend_code_mapbiomas_brazil_collection_11.csv`, downloaded 2026-10-25,
whose row order follows its level-1 groups: the first six rows 3, 6, 4, 7, 5,
49 are its *Forest* group; 9 Forest Plantation sits under *Farming*; 33 and
31 are *Water*).

**Legacy.** `git grep -n -i -E 'shyft|land_type|lake_fraction|forest_fraction|glacier_fraction|reservoir' legacy-archive -- legacy`
returns `legacy/rasputin/application.py:125-147` (`land_type_partition`),
`globcov_repository.py`, `gml_repository.py`, `wfs_repository.py`,
`web_visualize.py` and two tests. Intent, read in `application.py`: one class
per face, **the class at the face's centre point** (`lt_repo.land_cover(...,
geo_points=GeoPoints(xy=mesh.cell_centers[:, :2]))`), saved as the face field
`cover_type` for the Shyft coupling of the paper. A point sample, not a
partition; increment 26's exact fractions supersede it. Nothing carried over.

## 3. What a cell carries, and where each value comes from

| column | unit | from |
|---|---|---|
| `fid` | | the triangle's `fid` in the `triangles` layer (join key) |
| `geom` | POINT Z, run CRS | plan centroid of the triangle; z the mean corner height (= `z_mean`) |
| `catchment_id` | integer | 1 on every row: one run is one domain is one catchment (Q2) |
| `area_m2` | m² | plan area, as `triangles.area_m2` and Shyft's `calc_tin_horizontal_projected_area` |
| `z_m` | m | the midpoint's z, repeated as a column for tables without geometry |
| `slope_deg` | degrees | `atan2(hypot(n_x, n_y), n_z)` of the upward normal, as the `triangles` layer |
| `aspect_deg` | degrees clockwise from grid north, downslope | `atan2(n_x, n_y) mod 360`; NULL when `n_x = n_y = 0` |
| `glacier`, `lake`, `reservoir`, `forest` | plan-area share, 0 to 1 | section 4.1 and 4.2 |
| `unspecified` | share | `1 - glacier - lake - reservoir - forest` (Shyft's own definition) |

Grid north, not true north: the meridian convergence is not applied (Shyft's
`aspect()` does not either). Fractions are float64 sums of the file's float32
class fractions; the four never exceed 1 by more than `44 * 6e-8` (44 classes
at most in a shipped legend), so `unspecified` is clamped at 0 and the row
sums to 1 exactly by construction.

## 4. The rule

**4.1 Categories live in the legend files.** Each `[[class]]` in
`src_python/tin_engine/data/legends/*.toml` may carry
`category = "forest" | "lake" | "glacier"`; a class without one is
unspecified; a `[[split.class]]` inherits its parent's; class 0 (no data,
constraint lines) is unspecified. Legend check: a `water = true` class with
no category is refused naming the class, so water never silently falls into
unspecified. The shipped mapping:

| legend | forest | lake | glacier |
|---|---|---|---|
| `corine`, `corine-water`, `clc18_kode` | 311, 312, 313 | 511, 512, 521, 522, 523 | 335 |
| `ar50` | 30 (and its splits 3031, 3032, 3033, 3099) | 81, 82 | 70 |
| `mapbiomas-c11` | 3, 4, 5, 6, 7, 49 (its *Forest* group) and 9 Forest Plantation (Q1) | 33, 31 | none |

Left unspecified on purpose: CORINE 324 Transitional woodland-shrub and 244
Agro-forestry; AR50 60 Marsh. **A user's own mapping** is a legend file:
copy a shipped one, edit `category`, pass it with `--land-cover-legend
my.toml` (exists today). No new option. A user holding the EEA's own 100 m
CLC raster (grid codes 1 to 44, uint8) needs only such a legend file over
those codes: no code change.

**4.2 Reservoir (105b).** CORINE, AR50 and MapBiomas cannot tell a regulated
lake from a lake. `--reservoirs <catalogue vector source or file>` gives
reservoir polygons; per triangle `r_t` is the plan share inside them
(shapely intersection, polygons moved to the mesh CRS once), then
`reservoir_t = min(lake_t, r_t)` and `lake_t -= reservoir_t`. Exact when the
polygons are the lake outlines at the regulated level; otherwise bounded by
the lake share. Without `--reservoirs`, reservoir is 0 and the run record
says so (no warning: nothing was found wrong). Sources: **Norway**, NVE's
reservoir layer `VassdragsreguleringVannkraft/MapServer/9` (registered in
Ola's sources file as `nve-reservoirs`, cached 2026-10-17: 1,227 polygons
with `magasinNr`, `vatnLnr`, HRV and LRV; NLOD 1.0); NVE's lake database in
the catalogue (`nve-lakes`) collects no regulation field. **Sweden**: SMHI's
dam register is points, not polygons; the probe used the CORINE water polygon
holding the "Kultsjön" dam point, which is what a user would do by hand.
**São Francisco**: ANA's reservoir polygons (SNIRH); not fetched, licence
not checked; to be registered as a source before the basin run.

**Glacier** comes from the land-cover class (CORINE 335, AR50 70; MapBiomas
has none, and the basin has none). Outline overlays are left out (section 7):
NVE's Glacier Area Outline 2018-2019 (`../rasputin_data/nve_glaciers`, 6,736
outlines, shapefile, UTM 33N) is **CC BY-NC-SA 4.0, non-commercial**, so it
cannot ship in a product run without NVE's permission; RGI 7.0
(`../rasputin_data/rgi7`, regions 08 Scandinavia and 15 South Asia East,
shapefiles, WGS 84) is **CC BY 4.0** ("may be used freely under the
Creative Commons Attribution 4.0 International License", rgidata.org user
guide, read 2026-10-25). Both are shapefiles, which rasputin does not read
(`.gpkg`, `.geojson`, `.gml` only); a conversion step is the user's for now.

**4.3 Slope and aspect, no factor.** Section 2: Shyft's TIN cells use 1.0.
A consumer that wants Shyft's grid-style factor has slope and aspect to make
one; writing a constant column would be noise.

**4.4 The GeoPackage.** Two additions to the `.gpkg` `rasputin mesh` writes,
both whenever `--land-cover` is given and `--out` ends in `.gpkg`; no new
option:

- layer `cells`: POINT Z, one row per triangle in `fid` order, the columns of
  section 3, registered and R-tree indexed like `triangles`. A point, not a
  second polygon copy: QGIS joins it to `triangles` on `fid`, and the file
  does not double.
- attributes table `land_cover_fractions (fid, triangle INTEGER, code
  INTEGER, fraction REAL)`, the sparse fractions table in long form, one row
  per (triangle, class) entry, with an index on `triangle`. This is how the
  existing class fractions reach GeoPackage. Scale: the Kultsjön probe mesh
  has 2,254 triangles and the delivery's 10 m mesh 65,893; the CLI test scene
  has 1,825 entries over 1,590 triangles (1.15 per triangle). A São Francisco
  mesh at 10 m (about 25 million triangles by the Kultsjön ratio) would give
  about 30 million rows, about 1 GB of SQLite: acceptable for a one-off, and
  the coarse mesh this increment is for is three orders smaller.

Without a legend that has categories (a `--land-cover file.tif` with no
`--land-cover-legend`), the fractions table is written, the `cells` layer is
not, and the record says "cells: none, the land cover has no legend".
`.vtk` and `.ply` are unchanged.

**4.5 Fractions from polygons (105c).** `--land-cover` with a vector value
(a `.gpkg`, `.geojson` or `.gml` file, or a catalogue vector source with a
legend such as `clc2018` or Ola's `ar50-innlandet`): the class polygons are
read through `feature_input.read_source` with the legend, moved to the mesh
CRS, and each triangle is cut against the polygons an `STRtree` returns
(`shapely.intersection` areas), giving the same `ClassFractions` table as the
raster route. Uncovered area is class 0. This is the probe's route (section
8): 15.7 s for 2,254 triangles against 805 polygons, single-threaded. The
cost per triangle falls as triangles shrink (fewer polygons each), so the
65,893-triangle case is not 30 times longer; a scratch timing on it belongs to
105c's design update. A coded `--features` legend beside it stays refused, as
today.

## 5. Blueprint

```
legend .toml ──legends.load──▶ Legend(categories: Mapping[int, Category])        [+15 lines]
class raster ──land_cover_fractions.assign (26)──▶ ClassFractions (CSR)           [exists]
class polygons ──polygon_fractions.assign (105c)──▶ ClassFractions                [105c]
ClassFractions × categories ──cells.categorise──▶ (T, 4) float64                  [cells.py, pure]
vertices, triangles ──cells.geometry──▶ midpoint xyz, area, slope, aspect        [cells.py, pure]
reservoir polygons ──cells.overlay──▶ r_t, then min(lake, r_t)                   [105b]
CellTable (frozen dataclass of arrays) ──geopackage_write.write_cells──▶ `cells`
ClassFractions ──geopackage_write.write_fractions──▶ `land_cover_fractions`
cli.mesh: option checks, the two calls, three record lines                       [wiring only]
```

`cells.py` imports numpy and the two models only: no `_core`, no file, no
CRS object, so every test runs on hand arrays. `geopackage_write` keeps its
contract: an open connection in, no path, no `_core`. Pydantic holds the
legend (data); numpy holds the sums (execution). No new dependency; no C++;
blocking as `_write_geopackage` is today, inside `clock.phase("write:
geopackage")`.

## 6. The red suite (`@tester`), with measured fixtures

- **T1 legends.** Every shipped coded legend loads with the categories of
  section 4.1, checked as sets per category. A `.toml` with class 512
  `water = true` and no `category` is refused naming 512; `category = "sea"`
  is refused.
- **T2 categorise.** Offsets `[0, 2, 3, 5]`, classes `[311, 512, 335, 0,
  312]`, fractions `[0.25, 0.75, 1.0, 0.5, 0.5]`, `corine` categories →
  forest `[0.25, 0, 0.5]`, lake `[0.75, 0, 0]`, glacier `[0, 1, 0]`,
  reservoir 0, unspecified `[0, 0, 0.5]`; every row sums to 1 exactly.
- **T3 geometry.** F1 of `tests/python/test_io_geopackage_write.py`
  (vertices (0,0,0), (1,0,1), (1,1,2), (0,1,3); triangles (0,1,2), (0,2,3)):
  midpoints (2/3, 1/3, 1) and (1/3, 2/3, 5/3); areas 0.5, 0.5; slopes
  54.7356° and 72.4516°; aspects **225.0° and 161.5651°** (normals (−1,−1,1)
  and (1,−3,1)). Plus (0,0,0), (1,0,0), (0,1,1): slope 45°, aspect 180°
  (south); and the flat (0,0,0), (1,0,0), (0,1,0): slope 0, aspect NULL.
  Figures computed 2026-10-25 with numpy from the formulas of section 3.
- **T4 writer.** In a `:memory:` GeoPackage: `cells` registered as `POINT`
  with z, R-tree rows = triangles, the columns of section 3 with REAL and
  INTEGER types; `land_cover_fractions` rows = entries, index on `triangle`;
  `io.geopackage` reads the rows back in `fid` order.
- **T5 CLI.** The scene of `tests/python/test_cli_mesh_land_cover_fractions.py`
  (DEM 40 × 40 at 60 m in EPSG:31983, the synthetic class cache,
  `--tolerance 3`) with `--out x.gpkg --land-cover-legend mapbiomas-c11`,
  measured 2026-10-25 at `074e6322`: **1,590 triangles, 1,825 entries,
  classes {0, 9, 15, 21, 25, 30, 41, 75}, 1,367 triangles of one class**;
  `cells` has 1,590 rows, `land_cover_fractions` 1,825; forest equals the
  class-9 share per row (the only forest class present), lake, glacier and
  reservoir are 0 on every row, row sums 1 within 1e-6; `cells.area_m2` equals
  `triangles.area_m2` row by row, total 4,490,161 m² (relative 1e-9). Without
  `--land-cover` the content list stays `triangles, lines, run` (plus
  `land_cover_classes` when coded); `.vtk` output gains no array.
- **Mutation targets** (invariant-critical): `unspecified = 1 − sum` (a
  mutant summing the five), the aspect quadrant (`atan2(n_y, n_x)`), the
  category of class 0, and in 105b `min(lake, r)`.

## 7. Left out

Glacier outline overlays and a shapefile reader (4.2); `--catchment-id` and
multi-catchment tables (Q2); face area (Shyft's `rarea`, = `area_m2 /
cos(slope)`); true-north aspect; cells in `.vtk` or `.ply`; a Shyft NetCDF
writer (Shyft reads `cell_data` from NetCDF; netCDF4 would be a new
dependency, and the GeoPackage is one `sqlite3` read away for the user's own
script); any change to `--features` (constraints stay the exact route for
fine meshes, as the 104 showcase uses them).

## 8. The probe: Kultsjön's table before it is built

Commands and outputs: `docs/increments/105-probes/probe_cells.py`,
`probe_cells.out`. Run 2026-10-25 on this worktree's build, the delivery's
outline (`~/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon.geojson`,
1,702.9 km², 867 vertices, EPSG:3006), GLO-30 from the fetch cache.

1. **`--triangles 500` is refused**: "the fewest triangles this DEM gives is
   2254 triangles, at a largest error of 536 m; ask for at least 2254". The
   floor is the outline's own triangulation after the 25° start-quality pass;
   with `--start-min-angle 0` it is 865 (the outline's 867 vertices).
   About 500 cells needs a coarser outline; nothing in this increment.
2. **`--triangles 2254`** writes 2,254 triangles (mean 0.76 km²), tolerance
   536 m found in one trial, 753.7 kB GeoPackage.
3. **The land cover is CORINE polygons** (`~/.cache/rasputin/clc2018/features.gpkg`,
   EPSG:3035); the catalogue has no CORINE raster, so `--land-cover` cannot be
   run on Kultsjön today (4.5 is the fix). The probe cut the 2,254 triangles
   against the 805 polygons in the box with shapely: 15.7 s.

The first rows of the table (x, y in EPSG:3006; shares of plan area):

```
cell  x_m      y_m       z_m    area_km2 slope aspect glacier lake  resv  forest unspec
   1    456119   7230281  1158.3    0.019  13.4   91.8 0.000 0.000 0.000 0.000 1.000
   4    457161   7229225  1162.3    0.094  11.1    1.7 0.000 0.184 0.000 0.000 0.816
   7    456996   7227977  1284.8    0.166   7.2  100.4 0.000 0.000 0.000 0.000 1.000
```

Totals over 1,702.9 km², area-weighted: **glacier 0; lake 59.43 km²
(3.49 %); reservoir 55.85 km² (3.28 %, the Kultsjön water polygon, 4.2);
forest 480.77 km² (28.23 %); unspecified 1,106.81 km² (65.00 %, mostly 322
Moors and heathland, 803.8 km²)**. 1,596 of 2,254 triangles are one category
(share ≥ 0.999); every row sums to 1.000000. Check against the delivery's
constraint mesh (`meshes/kultsjon.gpkg`, `SELECT code, sum(area_m2) FROM
triangles GROUP BY code`): forest 311 + 312 + 313 = 378.16 + 79.56 + 23.05 =
480.77 km² and water 512 = 115.28 km², both equal to the probe's at the
printed 0.01 km², so the coarse partition loses no area against the exact
fine mesh.

## 9. Size

Estimates, each driver named: **105a** about 135 lines (legend model and
check +15, `cells.py` 50, `write_cells` and `write_fractions` 45, CLI wiring
and record lines 25; TOML is not counted). **105b** about 50 (polygon read
through the existing reader 15, overlay 25, option 10): short path. **105c**
about 120 (vector branch of `_class_source` 30, `polygon_fractions.py` 60,
record and refusals 30). Order: 105a, 105c, 105b, so São Francisco's table
exists first and Kultsjön's second (Q4). Each PR has its own red and green
commit; none touches refine or mesh code, so no `@perf` acceptance run.

## 10. Ideas-file entries taken up or overlapped

*A real storage format for meshes* (Ola, 2026-10-20): overlapped, the
fractions table moves one more array into the GeoPackage. The research note
in `docs/research/raster-to-vector.md` section 5 ("Burkhart et al. 2021:
cells carry land-type fractions ... a fraction per cell is a third option")
is what this builds. No ideas entry is taken up whole.

## 11. Questions for Ola

1. MapBiomas *Savanna Formation* (4, the cerrado) counts as forest, as
   MapBiomas's own level-1 group does. Default: yes.
2. `catchment_id` is 1 on every row; one run, one catchment. Default: yes.
3. No new option: the table is written whenever `--land-cover` meets a
   `.gpkg` output. Default: yes.
4. Build order 105a (raster, São Francisco), 105c (polygons, Kultsjön and
   Norway), 105b (reservoirs). Default: yes.
5. Glacier outlines, when wanted: RGI 7.0 (CC BY 4.0) rather than NVE's
   (CC BY-NC-SA 4.0). Default: RGI.

## Review

- **Design review, round 1, 2026-10-25.** `@reviewer`, head `bb7c8e893e8b4793f8515fc42b67449246531d73`, range `074e6322..bb7c8e89`. Verdict: CHANGES REQUESTED. LOC: 0 (design); about 305 estimated in three PRs. Holds: Shyft quotes match `geo_personality.h` (TIN slope factor 1.0, grid 0.9); category mapping matches the shipped legends; licences (RGI CC BY 4.0, NVE glaciers CC BY-NC-SA, reservoir layer NLOD 1.0, 1,227 polygons); no GDAL or shapefile reader; T2, T3 recompute; refusal at 2,254 (865 with `--start-min-angle 0`) re-run; probe sums agree. Blocking: (1) no timings: say no step runs without `--land-cover` and time 105a writes and the 105c cut on the delivery's Kultsjön mesh; (2) no test cases for 105b and 105c (the `min(lake, r)` mutant unkilled): add lake 0.5/r 0.8 -> reservoir 0.5, lake 0, and lake 0.5/r 0.2 -> reservoir 0.2, lake 0.3, and a polygon-cut case; (3) refusing an uncategorised water class at load breaks `tests/python/legend_fixtures.py:35,61` and `test_legends.py:761`: amend them in red or refuse only when the cells table is written; (4) wordings: the delivery mesh has 86,537 triangles, not 65,893 (4.4, 4.5), São Francisco about 32 million triangles and 37 million rows; "sums to 1 exactly" is up to 44 x 6e-8 over; the 15.7 s includes reading and reprojecting; Q3 "no new option" but 105b adds `--reservoirs`; "the 104 showcase" -> "the triangle-budget showcase". Not pushed; no CI.
