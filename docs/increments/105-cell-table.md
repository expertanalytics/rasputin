# Increment 105: the cell table — hydrological cells on a coarse mesh

**Status: designed, not built; design round 2.** By `@architect`, 2026-10-25,
on master `074e6322`, worktree `shyft-cells`; round 2 the same day after the
design review's round 1 (section *Review*) and the rulings the main session
relayed from Ola (section 11). From Ola's two sentences of 2026-10-25:
"For a coarse triangulation, we should provide the typical shyft-cell data.
Like lake fraction, forest fraction, etc." and "So we can do this already?
So, not use the land types as constraints, but simply partition the cell into
its types?" Number 105: 104 is the highest under `docs/increments/` over
every ref (`git for-each-ref` + `git ls-tree`, 2026-10-25).

**Size: bigger than the one-page threshold.** About 340 net production
lines in three parts (section 9), each its own PR under 150 lines, built in
this order: **105c** fractions from polygons (what Kultsjön and Norway need:
their land cover is vector), **105a** the table (for a class raster too, so
for São Francisco with MapBiomas), **105b** reservoirs. The letters are the
round-1 names, kept so the review record still reads. The probes (sections 8
and 8a) were run before any of it is built.

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
| `unspecified` | share | `max(0, 1 - glacier - lake - reservoir - forest)` (Shyft's own definition, clamped) |
| `savanna` | share | the savanna classes (4.1); **inside** `unspecified`, never beside it: `savanna ≤ unspecified` on every row |

Savanna is a sixth column because Ola ruled it a class of its own ("I'd take
Savanna as a separate class for sure", relayed 2026-10-25), and the cerrado
is most of the São Francisco basin. Shyft has no such fraction, so for a
Shyft reader the five Shyft columns still sum to 1 and savanna sits in
`unspecified`; a reader that wants it reads the sixth column. No other
category gets a column: the per-class shares are in `land_cover_fractions`.

Grid north, not true north: the meridian convergence is not applied (Shyft's
`aspect()` does not either). Fractions are float64 sums of the file's float32
class fractions. A row's class fractions each carry a float32 rounding of at
most `2^-24` (6e-8), so a row of `k` entries sums to 1 within `k * 6e-8`; `k`
is at most the legend's class count, 44 in the largest shipped legend
(`grep -c '^code = ' src_python/tin_engine/data/legends/*.toml`: corine,
corine-water and clc18_kode 44, mapbiomas-c11 33, ar50 20), so the four
categories exceed 1 by at most 2.7e-6. `unspecified` is clamped at 0, and the
five then sum to 1 within 2.7e-6: **not exactly**, and the tests say "within
3e-6". A user legend with more classes widens the bound by 6e-8 per class.

## 4. The rule

**4.1 Categories live in the legend files.** Each `[[class]]` in
`src_python/tin_engine/data/legends/*.toml` may carry
`category = "forest" | "lake" | "glacier" | "savanna" | "unspecified"`; a
class without one is unspecified; a `[[split.class]]` inherits its parent's;
class 0 (no data, constraint lines) is unspecified. Any other value is refused at load
(a Pydantic `Literal`). `Legend` gains two fields filled at load:
`categories: Mapping[int, Category]` (class code to its named category) and
`water_codes: frozenset[int]` (the classes with `water = true`).

**The legend check runs only when the cells table is written.** Loading never
refuses a class for lacking a category: the test legends
(`tests/python/legend_fixtures.py:35,61`, a `water = true` Lake and a
`water = true` 512 with no category) and every user legend written before
this increment keep loading, so no existing fixture changes. The pure
function `cells.categories_of(legend)` returns `legend.categories` or raises
`ValueError` naming every code in `legend.water_codes` that has no category
("class 512 is water but has no category; give it lake or unspecified"), so
water never silently falls into unspecified. `cli.mesh` calls it right after
`_class_source`, before any DEM is read, whenever `--land-cover` comes with a
legend and `--out` ends in `.gpkg`; the refusal is a usage error on
`--land-cover-legend`. The same legend with a `.vtk` output runs as today.
The shipped legends never trip it: their `water = true` classes are exactly
the lake and sea classes below (`grep -B3 'water = true'` over the five
files, 2026-10-25), and every one carries a category. The shipped mapping:

| legend | forest | lake | glacier | savanna | unspecified, named |
|---|---|---|---|---|---|
| `corine`, `corine-water`, `clc18_kode` | 311, 312, 313 | 511, 512, 521, 522 | 335 | none | 523 Sea and ocean |
| `ar50` | 30 (and its splits 3031, 3032, 3033, 3099) | 81 | 70 | none | 82 Ocean |
| `mapbiomas-c11` | 3 Forest Formation, 5 Mangrove, 6 Floodable Forest, 49 Wooded Sandbank Vegetation, 9 Forest Plantation | 33, 31 | none | 4 Savanna Formation, 7 Flooded Savanna (beta) | none |

Savanna is its own category (Ola, relayed 2026-10-25, section 11): MapBiomas
4 is the cerrado, scattered trees over grass and shrub, which MapBiomas files
under its level-1 *Forest* group but Shyft's forest (a canopy that intercepts
rain and snow) is not. 7 Flooded Savanna goes with it by its own name (a
default taken, section 11). Sea is unspecified, not lake (same ruling): Shyft's
lake share is water that stores and evaporates inside the catchment, and a
coastal triangle over the sea is outside it. The two sea classes carry
`category = "unspecified"` so the check passes them on purpose. Left
unspecified without a word: CORINE 324 Transitional woodland-shrub and 244
Agro-forestry; AR50 60 Marsh. **A user's own mapping** is a legend file:
copy a shipped one, edit `category`, pass it with `--land-cover-legend
my.toml` (exists today). No new option. A user holding the EEA's own 100 m
CLC raster (grid codes 1 to 44, uint8) needs only such a legend file over
those codes: no code change.

**4.2 Reservoir (105b).** CORINE, AR50 and MapBiomas cannot tell a regulated
lake from a lake. **One new option**, `--reservoirs <catalogue vector source
or file>` (Ola, relayed 2026-10-25: wanted), gives reservoir polygons, read
through `feature_input.read_source` and moved to the mesh CRS once. Per
triangle, `r_t = area(t ∩ R) / area(t)` where `R = shapely.union_all` of the
polygons, so two outlines of one reservoir (an upper and a lower regulated
level, say) count once; then `reservoir_t = min(lake_t, r_t)` and
`lake_t -= reservoir_t`. Exact when the polygons are the lake outlines at the
regulated level; otherwise bounded by the lake share, so a reservoir polygon
over land the land cover calls forest adds nothing. `--reservoirs` without
`--land-cover`, or with a `.vtk` output, is a usage error ("needs
--land-cover and a .gpkg output"). Without `--reservoirs`, reservoir is 0 and
the run record says so (no warning: nothing was found wrong). Sources: **Norway**, NVE's
reservoir layer `VassdragsreguleringVannkraft/MapServer/9` (registered in
Ola's sources file as `nve-reservoirs`, cached 2026-10-17: 1,227 polygons
with `magasinNr`, `vatnLnr`, HRV and LRV; NLOD 1.0); NVE's lake database in
the catalogue (`nve-lakes`) collects no regulation field. **Sweden**: SMHI's
dam register is points, not polygons; the probe used the CORINE water polygon
holding the "Kultsjön" dam point, which is what a user would do by hand.
**São Francisco**: ANA's reservoir polygons (SNIRH); not fetched, licence
not checked; to be registered as a source before the basin run.

**Glacier** comes from the land-cover class (CORINE 335, AR50 70; MapBiomas
has none, and the basin has none). Outline overlays are left out (section 7).
When they are built, Ola's ruling (relayed 2026-10-25, "5: default" to the
main session's proposal) is: **Norway from NVE's Glacier Area Outline
2018-2019** (`../rasputin_data/nve_glaciers`, 6,736 outlines, shapefile, UTM
33N; Sentinel-2, newer and finer than RGI's year-2000 outlines for Norway),
which is **CC BY-NC-SA 4.0, non-commercial**: a run for a paying customer
cannot use it without NVE's permission, and the run record must say so;
**elsewhere RGI 7.0** (`../rasputin_data/rgi7`, regions 08 Scandinavia and
15 South Asia East, shapefiles, WGS 84), **CC BY 4.0** ("may be used freely
under the Creative Commons Attribution 4.0 International License",
rgidata.org user guide, read 2026-10-25). Both are shapefiles, which rasputin
does not read (`.gpkg`, `.geojson`, `.gml` only); a conversion step is the
user's for now.

**4.3 Slope and aspect, no factor.** Section 2: Shyft's TIN cells use 1.0.
A consumer that wants Shyft's grid-style factor has slope and aspect to make
one; writing a constant column would be noise.

**4.4 The GeoPackage.** Two additions to the `.gpkg` `rasputin mesh` writes,
both whenever `--land-cover` is given and `--out` ends in `.gpkg` (Ola,
relayed 2026-10-25: yes to both tables, no option to ask for them):

- layer `cells`: POINT Z, one row per triangle in `fid` order, the columns of
  section 3, registered and R-tree indexed like `triangles`. A point, not a
  second polygon copy: QGIS joins it to `triangles` on `fid`, and the file
  does not double.
- attributes table `land_cover_fractions (fid, triangle INTEGER, code
  INTEGER, fraction REAL)`, the sparse fractions table in long form, one row
  per (triangle, class) entry, with an index on `triangle`. This is how the
  existing class fractions reach GeoPackage. Scale: the Kultsjön probe mesh
  has 2,254 triangles and the delivery's 10 m mesh 86,537
  (`sqlite3 ~/rasputin-output/opal-se2-2026-10-21/meshes/kultsjon.gpkg
  'SELECT count(*) FROM triangles'`, 2026-10-25); the CLI test scene has
  1,825 entries over 1,590 triangles (1.15 per triangle). A São Francisco
  mesh at 10 m (about 32 million triangles: 86,537 per 1,702.9 km² is 51 per
  km², times 637,000 km²) would give about 37 million rows, about 1 GB of
  SQLite: acceptable for a one-off, and the coarse mesh this increment is for
  is three orders smaller. Section 8a has the measured write times.

Without a legend that has categories (a `--land-cover file.tif` with no
`--land-cover-legend`), the fractions table is written, the `cells` layer is
not, and the record says "cells: none, the land cover has no legend".
105a leaves `.vtk` and `.ply` unchanged (105c changes the type of one
existing array in each, 4.5).

**4.5 Fractions from polygons (105c).** `--land-cover` with a vector value
(a `.gpkg`, `.geojson` or `.gml` file, or a catalogue vector source with a
legend such as `clc2018` or Ola's `ar50-innlandet`): the class polygons are
read through `feature_input.read_source` with the legend, moved to the mesh
CRS, **subdivided**, and each triangle is cut against the pieces an `STRtree`
returns (`shapely.intersection` areas without `grid_size`, so no snapping, summed per class), giving the same
`ClassFractions` table as the raster route. Uncovered area is class 0. A coded
`--features` legend beside it stays refused, as today.

*Subdivision.* Each class polygon is halved along the longer side of its bounding box with `shapely.clip_by_rect` until every piece has at most **64 coordinates** (`shapely.get_num_coordinates`: closing points and holes count; PostGIS's `ST_Subdivide` idea). `clip_by_rect` output is not guaranteed valid; pieces are not repaired, and T6's area sums guard them; a piece keeps its polygon's
class, and the sum over pieces equals the sum over the polygon because the
pieces tile it. Without it the cut is 7 ms per triangle at any mesh size,
paid in the vertices of the CORINE polygon each triangle meets (section 8a:
610.7 s for the delivery's 86,537 triangles); with it, 2.7 s plus 1.0 s to
subdivide. The constant 64 assumes polygons of CORINE's and AR50's kind
(hundreds to thousands of vertices each, the largest 803.8 km²) and was
checked at 86,537 triangles against 805 polygons (25,041 pieces); 256 gave
the same time within noise on the 2,254 mesh (section 8a), so the choice is
not sharp. `polygon_fractions.subdivide` is its own pure function, tested by
itself (T6).

**Class numbers widen to int32 (105c).** `ClassFractions.classes` is `uint8`
today (`src_python/tin_engine/io/models.py`), written as `unsigned_char` in
`.vtk` (`vtk_legacy.py`) and `list uchar uchar` in `.ply` (`ply.py`): a class
raster is 8-bit. CORINE's `Code_18` (111 to 523) and AR50's splits (3031 to
3099) do not fit, and the Code type allows up to `2^31 - 1`
(`legends.py:26`). 105c makes `classes` `int32`, the type `land_cover_code`
already has, so one class number has one type in every file: `.vtk` writes
`int`, `.ply` `list uchar int`, the raster route casts `class_areas`'s uint8
once in `fractions_of`, and the readers follow the declared type
(`vtk_read.py` knows `int`; `tests/python/vtkread.py` and `plyread.py` read
the declared type). About 10 lines; the `.vtk` and `.ply` arrays change
type, and a reader outside this repository (the private package's) must read
the declared type rather than assume `unsigned_char` (section 11, item 6).
Cost: 3 more bytes per entry, about 110 MB on a São Francisco `.vtk` of 37
million entries.

## 5. Blueprint

```
legend .toml ──legends.load──▶ Legend(categories, water_codes)                   [+15 lines]
Legend ──cells.categories_of──▶ Mapping[int, Category], or ValueError            [cells.py, pure]
class raster ──land_cover_fractions.assign (26)──▶ ClassFractions (CSR)           [exists]
class polygons ──polygon_fractions.subdivide, .assign (105c)──▶ ClassFractions    [105c]
ClassFractions × categories ──cells.categorise──▶ (T, 5) float64                  [cells.py, pure]
   (glacier, lake, forest, savanna, unspecified; reservoir 0 until 105b)
vertices, triangles ──cells.geometry──▶ midpoint xyz, area, slope, aspect        [cells.py, pure]
triangles, reservoir union ──cells.overlay──▶ r_t                                 [105b, pure]
(lake, r_t) ──cells.split_reservoir──▶ (lake − min(lake, r), min(lake, r))        [105b, pure]
CellTable (frozen dataclass of arrays) ──geopackage_write.write_cells──▶ `cells`
ClassFractions ──geopackage_write.write_fractions──▶ `land_cover_fractions`
cli.mesh: option checks, the two calls, three record lines                       [wiring only]
```

`cells.py` imports numpy, shapely and the two models only: no `_core`, no
file, no CRS object, so every test runs on hand arrays and hand polygons.
`categorise` sums with `np.isin` per category and `np.bincount` over the
row's triangle index, so it is indifferent to the class number's width.
`overlay` takes triangles as corner arrays and the reservoir union as one
shapely geometry already in the mesh CRS; the reading and the CRS move are
`feature_input`'s, as for `--features`. `geopackage_write` keeps its
contract: an open connection in, no path, no `_core`. Pydantic holds the
legend (data); numpy holds the sums (execution). No new dependency; no C++;
blocking as `_write_geopackage` is today, inside `clock.phase("write:
geopackage")`.

**Nothing runs without `--land-cover`.** Every new call sits under the
existing `if classes is not None` branch of `cli.mesh` or inside
`_write_geopackage` behind a `fractions is not None` test; a run without
`--land-cover` executes no new statement but those two `None` tests, and its
GeoPackage, `.vtk` and `.ply` are byte for byte what they are today. T5 pins
the GeoPackage content list and the `.vtk` arrays of such a run.

## 6. The red suite (`@tester`), with measured fixtures

- **T1 legends.** Every shipped coded legend loads with the categories of
  section 4.1, checked as sets per category (sea in `unspecified`,
  `mapbiomas-c11` savanna `{4, 7}`), and
  `water_codes` equal to its `water = true` classes. A `.toml` with class 512
  `water = true` and no `category` **loads** (the fixture of
  `legend_fixtures.py:61` is such a legend) and `cells.categories_of` refuses
  it naming 512; `category = "sea"` is refused at load; a split's sub-classes
  inherit the parent's category (`ar50`: 3031 is forest).
- **T2 categorise.** Offsets `[0, 2, 3, 5]`, classes `[311, 512, 335, 0,
  312]`, fractions `[0.25, 0.75, 1.0, 0.5, 0.5]`, `corine` categories →
  forest `[0.25, 0, 0.5]`, lake `[0.75, 0, 0]`, glacier `[0, 1, 0]`,
  reservoir 0, unspecified `[0, 0, 0.5]`; every row sums to 1 exactly (the
  fixture's fractions are binary fractions, exact in float32). A fourth row
  with classes `[523, 311]` and fractions `[0.75, 0.25]` → lake 0, forest
  0.25, unspecified 0.75: sea is unspecified. With `mapbiomas-c11`
  categories, a row with classes `[3, 4]` and fractions `[0.25, 0.75]` →
  forest 0.25, savanna 0.75, unspecified 0.75 (savanna inside it), and
  glacier + lake + reservoir + forest + unspecified = 1. A row whose float32 fractions
  `[0.4, 0.6]` (classes 311, 512) sum to `1 + 3.0e-8` in float64 (measured
  with numpy 2026-10-25) → unspecified 0, not `-3.0e-8`, and the five sum to
  1 within 3e-6.
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
  class-9 share per row (the only forest class present), lake, glacier,
  savanna and reservoir are 0 on every row (no class 4 or 7 in the scene),
  the five Shyft columns sum to 1 within 3e-6; `cells.area_m2` equals
  `triangles.area_m2` row by row, total 4,490,161 m² (relative 1e-9). Without
  `--land-cover` the content list stays `triangles, lines, run` (plus
  `land_cover_classes` when coded); `.vtk` output gains no array. **The legend
  check:** the same scene with `--land-cover-legend lake.toml`, a `coded_map`
  legend whose Lake class (code 2, `water = true`) has no category, and
  `--out x.gpkg` exits 2 with a usage error naming `--land-cover-legend` and
  class 2 before any DEM is read (no `x.gpkg.part` exists afterwards); the
  same with `--out x.vtk` writes the mesh with its three fraction arrays as
  today; `--reservoirs` (105b) with `--out x.vtk` or without `--land-cover`
  is a usage error.
- **T6 polygon fractions (105c).** `polygon_fractions.assign` on triangles
  `A = (0,0), (2,0), (0,2)` and `B = (2,0), (2,2), (0,2)` (plan area 2 each)
  against class polygons 311 = the strip `[0,1] × [0,2]` and 512 = the square
  `[1,2] × [0,1]`: A is 1.5 in the strip and 0.5 in the square (the corner
  triangle `(1,0), (2,0), (1,1)`); B is 0.5 in the strip (`(0,2), (1,1),
  (1,2)`), 0.5 in the square (`(1,1), (2,0), (2,1)`) and 1.0 uncovered. So
  offsets `[0, 2, 5]`, classes `[311, 512, 0, 311, 512]` (ascending within a
  row, as `ClassFractions` requires), fractions `[0.75, 0.25, 0.5, 0.25,
  0.25]`, dominant codes `[311, 0]`; `classes.dtype` is `int32`. Figures by
  hand, checked with shapely's `intersection(...).area` 2026-10-25, as were
  T7's. A polygon of code 9999 absent from the legend is refused
  naming 9999; two polygons of one class over one triangle give one entry.
  `subdivide`: a 1,000-gon (`shapely.Point(1, 1).buffer(1, quad_segs=250)`,
  1,001 coordinates, inscribed in the square the triangles A and B tile) with
  `most = 64` gives 32 pieces of at most 58 coordinates whose areas sum to
  the 1,000-gon's within 1e-9 relative, and the share of A and of B is
  0.785393 both ways, over the pieces and over the whole (half the 1,000-gon
  each; π/4 = 0.785398 for the true circle). Figures from `probe_cells.py`'s
  `subdivide` with shapely, 2026-10-25. A polygon already under `most` comes
  back as itself, one piece.
  The CLI case: the T5 scene with `--land-cover polys.geojson` (two
  rectangles in the scene's EPSG:31983, written with
  `feature_fixtures.write_geojson`) and a `coded_map` legend on their
  attribute: the written fractions equal shapely's intersection areas over
  the `triangles` layer read back (bound 1e-6), rows sum to 1, and
  `land_cover_code` is the largest share; the `.vtk` array
  `land_cover_classes` is declared `int`.
- **T7 reservoir split (105b).** `split_reservoir(lake, r)`, the reviewer's
  two cases and two more: lake 0.5, r 0.8 → reservoir 0.5, lake 0.0 (kills
  `reservoir = r`); lake 0.5, r 0.2 → reservoir 0.2, lake 0.3 (kills
  `reservoir = lake` and a dropped `lake -= reservoir`); lake 0.0, r 0.6 →
  both 0 (a reservoir polygon over land); lake 0.5, r 0.0 → unchanged.
  `overlay` on the triangles A and B of T6: the square `[0,1]²` lies in A
  (its far corner `(1,1)` is on A's hypotenuse) → `r = [0.5, 0]`; the square
  `[1,2]²` lies in B → `r = [0, 0.5]`; the square `[0.5,1.5]²`, cut in half by
  the line `x + y = 2` → `r = [0.25, 0.25]`; the union of `[0,1]²` given
  twice → `[0.5, 0]`, once. The CLI case: the T5 scene with `--reservoirs
  res.geojson` holding one rectangle over the water class of a `coded_map`
  legend: `reservoir + lake` equals the lake share without `--reservoirs`
  row by row, and `reservoir ≤ lake` without it on every row; the run record
  says `reservoirs: none` when the option is absent.
- **Mutation targets** (invariant-critical): `unspecified = 1 − sum` (a
  mutant summing the five, and one dropping the clamp), the aspect quadrant
  (`atan2(n_y, n_x)`), the category of class 0, the sea class's category
  (lake instead of unspecified), savanna subtracted from unspecified (the
  T2 row `[3, 4]` would give unspecified 0), the `water_codes − categories`
  set in `categories_of` (a mutant checking all classes, which refuses
  `corine`), in 105c the subdivision's piece count or an unclipped half
  (areas no longer sum to the ring's), and in 105b `min(lake, r)` and
  `lake -= reservoir`.

## 7. Left out

Glacier outline overlays and a shapefile reader (4.2); `--catchment-id` and
multi-catchment tables (Q2); face area (Shyft's `rarea`, = `area_m2 /
cos(slope)`); true-north aspect; cells in `.vtk` or `.ply`; a Shyft NetCDF
writer (Shyft reads `cell_data` from NetCDF; netCDF4 would be a new
dependency, and the GeoPackage is one `sqlite3` read away for the user's own
script); any change to `--features` (constraints stay the exact route for
fine meshes, as the triangle-budget showcase uses them).

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
   against the 805 polygons in the box with shapely: 15.7 s in round 1, of
   which (round 2's re-run, `probe_cells.out`, phases printed) reading the
   mesh is 0.0 s, reading and reprojecting the 805 polygons 0.3 s and the
   STRtree query plus cut 15.6 s; 4,034 (triangle, class) entries, 1.79 per
   triangle. Section 8a has the same probe on the 86,537-triangle mesh.

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

## 8a. Timings: what the increment adds to a run

All runs 2026-10-25 on this worktree's venv, **on battery** (`pmset -g batt`:
94 % discharging at the start, 90 % at the end; an AC run would be faster),
single-threaded Python, one run each unless said. Outputs under
`docs/increments/105-probes/`. Nothing below runs without `--land-cover`
(section 5).

**105a, the two writes** (`probe_write.py`, random geometry of the stated
size through the product's `write_triangles` and the design's two tables;
`probe_write.out`):

| triangles | fraction entries | `triangles` layer today | `cells` layer | `land_cover_fractions` |
|---|---|---|---|---|
| 2,254 | 4,034 | 0.01 s | 0.01 s | 0.00 s |
| 86,537 (run twice) | 117,333 | 0.46 s, 0.45 s | 0.45 s, 0.44 s | 0.09 s, 0.09 s |

The `cells` layer costs what the `triangles` layer costs (one row and one
R-tree entry per triangle); the fractions table a fifth of that. Against the
delivery's whole `rasputin mesh` run (minutes of refinement at 10 m) both
are noise; at São Francisco's 32 million triangles, about 170 s and 35 s by
proportion. The file size the probe prints is for random, overlapping
triangles and says nothing about a real mesh.

**105c, the polygon cut** (`probe_cells.py`, the CORINE polygons of the
Kultsjön box, 805 polygons; `probe_cells.out` and `probe_cells_86537.out`):

| triangles | read mesh | read and reproject polygons | STRtree query and cut | entries |
|---|---|---|---|---|
| 2,254 | 0.0 s | 0.3 s | 15.6 s | 4,034 (1.79 per triangle) |
| 86,537 | 0.3 s | 0.3 s | **610.7 s** | 117,333 (1.36 per triangle) |

Both runs give the same category totals to the printed 0.01 km² (lake 59.43,
reservoir 55.85, forest 480.77, unspecified 1,106.81 km²), so the cut is
right at both sizes; it is the cost that is wrong. Round 1's hope that the
cost per triangle falls with triangle size was false: it is 7 ms per triangle
at both sizes, because each `shapely.intersection` is paid in the vertices of
the CORINE polygon it meets (the 322 Moors polygon alone is 803.8 km², one
ring of thousands of vertices), not in the triangle. Ten minutes for the
delivery mesh is not acceptable: the exact raster route does the same job in
C++ in seconds, and the constraint route cut these very polygons into the
mesh for the delivery.

**The fix, measured** (`probe_cells.py MESH MOST_VERTICES`,
`probe_cells_subdivided.out`): subdivide every polygon into pieces of at most
`MOST_VERTICES` vertices first (section 4.5).

| triangles | at most | pieces | subdivide | STRtree query and cut | whole probe |
|---|---|---|---|---|---|
| 2,254 | 256 | 5,909 | 0.8 s | 0.3 s | 1.3 s |
| 2,254 | 64 | 25,041 | 1.0 s | 0.2 s | 1.5 s |
| 86,537 | 64 | 25,041 | 1.0 s | **2.7 s** | 4.3 s |

Every run gives the totals above to 0.01 km² and the same entry counts
(4,034 and 117,333): the pieces tile the polygons. 610.7 s becomes 2.7 s
(225 times); the subdivision itself is a fixed 1 s for the box's 805
polygons, independent of the mesh. By proportion, a São Francisco mesh of 32
million triangles would cut in about 17 minutes single-threaded with
MapBiomas as polygons, which nobody will do: its land cover is a raster, and
the raster route is C++. The vector route is for CORINE and AR50 domains of
up to a few hundred thousand triangles, where it is seconds.

## 9. Size and order

Estimates, each driver named: **105a** about 140 lines (legend fields +15,
`cells.py` with `categories_of` 55, `write_cells` and `write_fractions` 45,
CLI wiring, the legend check and record lines 25; TOML is not counted).
**105b** about 55 (polygon read through the existing reader 15, `overlay`
and `split_reservoir` 25, the `--reservoirs` option and its two refusals
15): over the short path's 50 by a hair, so the full path. **105c** about
145 (vector branch of `_class_source` 30, `polygon_fractions.py` 75 of which
`subdivide` 15, the int32 widening in `models.py`, `vtk_legacy.py`, `ply.py`
and `land_cover_fractions.py` 10, record and refusals 30). **Order: 105c
first** (Ola, relayed 2026-10-25: "Polygons first"), then 105a, then 105b:
after 105c a CORINE or AR50 polygon source gives every triangle its exact
class fractions in `.vtk` and `.ply`, as a class raster does today; 105a
adds the two GeoPackage tables for both routes; 105b the reservoir split.
Each PR has its own red and green commit; none touches refine or mesh code,
so no `@perf` acceptance run.

## 10. Ideas-file entries taken up or overlapped

*A real storage format for meshes* (Ola, 2026-10-20): overlapped, the
fractions table moves one more array into the GeoPackage. The research note
in `docs/research/raster-to-vector.md` section 5 ("Burkhart et al. 2021:
cells carry land-type fractions ... a fraction per cell is a third option")
is what this builds. No ideas entry is taken up whole.

## 11. Questions for Ola, and his rulings so far

**Ruled** (relayed by the main session, 2026-10-25, Ola's words quoted; the
design follows them): 1 savanna is a class of its own ("I'd take Savanna as a
separate class for sure": the `savanna` column, section 3); 2 `catchment_id`
is 1 on every row ("ok"); 3 the table is written whenever `--land-cover`
writes a `.gpkg` ("yes"); 4 polygons first ("Aha, I understand. Polygons
first.": 105c before 105a, section 9); 5 glacier outlines, when built, from
NVE for Norway and RGI 7.0 elsewhere ("default", to the main session's
proposal; section 4.2); 7 the per-class `land_cover_fractions` table is
wanted ("yes"); 8 a `--reservoirs` option is wanted ("yes"); 9 sea goes to
unspecified, not lake ("unspecified").

**Defaults taken without a ruling**, each flippable by a legend file:
MapBiomas 7 *Flooded Savanna (beta)* is savanna with 4, by its name; 9
*Forest Plantation* stays forest; CORINE 324 Transitional woodland-shrub and
244 Agro-forestry, and AR50 60 Marsh, stay unspecified.

**Default taken, not put to Ola:**

6. 105c widens class numbers in `.vtk` and `.ply` to `int` for every land-cover run, raster runs too (4.5). Default: yes. `git grep -n land_cover_classes private/master` returns nothing, and the private package's `tests/python/vtkread.py` reads the declared type (`"int": ">i4"`), so no private reader breaks.

## Review

- **Design review, round 1, 2026-10-25.** `@reviewer`, head `bb7c8e893e8b4793f8515fc42b67449246531d73`, range `074e6322..bb7c8e89`. Verdict: CHANGES REQUESTED. LOC: 0 (design); about 305 estimated in three PRs. Holds: Shyft quotes match `geo_personality.h` (TIN slope factor 1.0, grid 0.9); category mapping matches the shipped legends; licences (RGI CC BY 4.0, NVE glaciers CC BY-NC-SA, reservoir layer NLOD 1.0, 1,227 polygons); no GDAL or shapefile reader; T2, T3 recompute; refusal at 2,254 (865 with `--start-min-angle 0`) re-run; probe sums agree. Blocking: (1) no timings: say no step runs without `--land-cover` and time 105a writes and the 105c cut on the delivery's Kultsjön mesh; (2) no test cases for 105b and 105c (the `min(lake, r)` mutant unkilled): add lake 0.5/r 0.8 -> reservoir 0.5, lake 0, and lake 0.5/r 0.2 -> reservoir 0.2, lake 0.3, and a polygon-cut case; (3) refusing an uncategorised water class at load breaks `tests/python/legend_fixtures.py:35,61` and `test_legends.py:761`: amend them in red or refuse only when the cells table is written; (4) wordings: the delivery mesh has 86,537 triangles, not 65,893 (4.4, 4.5), São Francisco about 32 million triangles and 37 million rows; "sums to 1 exactly" is up to 44 x 6e-8 over; the 15.7 s includes reading and reprojecting; Q3 "no new option" but 105b adds `--reservoirs`; "the 104 showcase" -> "the triangle-budget showcase". Not pushed; no CI.
- **Design review, round 2, 2026-10-25.** `@reviewer`, head `123820cb0f6f08d2128c99733c84f40b9065ecd7`, range `4163bc45..123820cb`. Verdict: CHANGES REQUESTED, only for Ola's ruling on widening class numbers to `int` in `.vtk` and `.ply` for every land-cover run (now question 6, default yes; `git grep land_cover_classes private/master` finds no reader). LOC: 0 (design); estimate 140 + 55 + 145 = 340, each part under 150. Round 1 closed: section 8a timings match the probes (610.7 s -> 2.7 s, same totals); T6, T7 recomputed; legend check only when the cells table is written (call before the DEM read, `cli.py:907`); five wordings fixed. Four wordings applied in this commit: section 9 heading restored; "64 coordinates" with what counts; clip_by_rect validity and intersection without snapping; question 6 replaces the pin. Not pushed; no CI.
- **Design review, round 3, 2026-10-25.** `@reviewer`, head `f8621a02a3988deade98d479b38ba7433cc40487`, range `123820cb..f8621a02`. Verdict: APPROVED WITH WORDINGS. LOC: 0 (design); estimate 140 + 55 + 145 = 340 in three PRs. Nothing needs Ola to rule; the `int` widening is a default taken, not put to Ola. Applied in this commit: the no-snapping note moved to the per-triangle intersection; the subdivision sentence rewritten; "(section 11, item 6)"; section 11 heading "Default taken, not put to Ola"; the private reader check (`tests/python/vtkread.py` reads the declared type). Not pushed; no CI.
