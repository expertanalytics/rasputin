# Audit: does rasputin work anywhere on the globe, or only on Ola's machine and data?

Status: **audit, not a design.** Written by `@architect` against master
`c317c7e4` (2026-10-12), after the land-cover fractions (#235) merged. Nothing
here changes code; every fix named is a candidate for a ROADMAP row Ola picks.
Where a claim needed a run, the run is in the appendix and its output is quoted
in the text.

The question, Ola's words (2026-10-12): "see if we're making a system that is
good only on my computer or datasets, and not generally applicable everywhere
on the globe"; "at places, the implementation is extremely concrete, perhaps
too much"; and his first finding, "the number of cores should of course not be
hardcoded".

## Verdict

**The geometry core and the coordinate-system layer are general; the
concreteness sits in four places around them.**

1. **Nothing lets a user set the thread count.** No command has a `--threads`
   option (`grep -n -- '--threads' src_python/tin_engine/cli.py` finds none).
   The C++ takes a `threads` argument everywhere and defaults to the machine's
   core count, which is the right default; the Python pools default to
   `os.cpu_count()`; one of them caps itself at 8 whatever the machine. Every
   run uses every core, and the stats file prints the core count as if it were
   the thread count. Severity: degrades (a shared box, a laptop on battery, a
   container, an agent running several meshes at once).
2. **Land-cover schemes live in Python code, and all of them are CORINE.**
   Three class maps (`corine`, `corine-water`, `clc18_kode`), one palette, the
   water codes 511-523, the CORINE copyright notice. Another scheme (MapBiomas
   polygons, ESA WorldCover, NLCD (the United States' National Land Cover
   Database), Norway's AR5 (its land resource map), the Austrian cadastre) needs a
   pull request, and schemes whose codes are not integers are refused. The
   remote-source catalogue is the same shape: three DEM sources and one station
   list, in code. Severity: blocks use with any other scheme.
3. **The land-cover clean-up defaults are CORINE-on-a-10-m-DEM numbers
   presented as general.** Repair 1 m, simplify 50 m, snap 5 m, start angle 15
   degrees, read margin 100 m; the gauge and river-burn constants (30 m floor,
   100 m, 500 m, 2000 m window margin) are DTM10 figures in metres. On a 1 m DEM
   with 1:5000 land cover the 50 m simplification destroys the input; on a 30 m
   DEM the 30 m floor is one cell. Severity: degrades.
4. **`catchment` refuses a DEM in degrees** (GLO-30, ANADEM), because it has no
   `--out-crs` path: the core boundary refuses degrees, and only `mesh` resamples
   first. The São Francisco basin is on ANADEM and GLO-30, so this is the one
   finding that blocks the stated goal. A DEM whose tiles come in two CRSs is
   refused too. Severity: blocks.

What is general: hemisphere and zone are handled through PROJ with no sign
assumption anywhere (the probe below runs southern UTM zones, an
equator-crossing box, both poles and a feet-based CRS through the CRS layer);
`suggest_crs` picks a conformal family for any box; the lattice frame is
`y = -row`, independent of where north is; locale appears nowhere; the C++ core
holds no CRS, no unit and no path. The macOS-only code is the benchmark tooling,
not the product.

**50 rows: 2 block use elsewhere, 23 degrade it, 25 are cosmetic.** The ordered
fix list and a PR split are in section 6.

## 1. How this was checked

Every row below is pinned to `c317c7e4` as `path@c317c7e4:line`, which
`python3 tools/check_citations.py` resolves. The sweeps:

- Core count and threads: `grep -rn -E 'cpu_count|hardware_concurrency|WORKERS|max_workers|threads' src_python/tin_engine include/terrain bindings tools/bench.py tools/bench_quick.py`.
- CRS and hemisphere: `grep -rn -E 'EPSG|25833|3006|3035|4326|UTM|hemisphere|northing|zone|geographic|is_geographic' src_python/tin_engine ...`.
- Land cover, DEM sources, data root: `grep -rn -i -E 'corine|mapbiomas|CORINE_NATURAL|palette|class_map|anadem|glo.?30|dtm10|kartverket|nve|RASPUTIN_DATA' ...`.
- Module-level constants: `grep -n -E '^[A-Z_]+(, [A-Z_]+)* *(: *[^=]+)? *= ' src_python/tin_engine/*.py src_python/tin_engine/io/*.py src_python/tin_engine/fetch/*.py` (every hit judged in section 2).
- C++ numeric constants: `grep -rn -E 'constexpr .*= *[0-9]|1 << |<< 18' include/terrain bindings`.
- Platform: `grep -rn -E 'darwin|macOS|pmset|caffeinate|sysctl|sw_vers|platform\.' src_python tools`, `grep -rn -E 'ffp-contract|Accelerate|APPLE|march=' CMakeLists.txt tests/cpp/CMakeLists.txt`.
- Locale and dates: `grep -rn -E 'locale|setlocale|strftime|%d\.%m|decimal' src_python/tin_engine`.
- Norwegian: `grep -rn -E '[æøåÆØÅ]|navn|vatnlnr|norsk|Norway' src_python/tin_engine`.
- The probe (appendix): `crs.suggest_crs`, `crs.reprojector`, `feature_input._bound`'s
  applicability test, `target_grid.default_spacing`, `dem_input._reprojected`
  and `crs.srs_entry` run on this tree's code in a scratch interpreter with
  `PROJ_NETWORK=OFF`, bypassing the editable finder (the first run silently
  imported another worktree's package; the script now asserts the path).

Severity: **blocks** means a user elsewhere is refused or gets a wrong answer
and only a code change fixes it; **degrades** means there is a work-around or a
worse result; **cosmetic** means naming, tooling, or a documented limit that is
rare in practice. Line estimates are net production lines (`CLAUDE.md` §2).

## 2. Findings

### A. Coordinate reference systems and hemisphere

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| A1 | `src_python/tin_engine/raster.py@c317c7e4:26-27`; `src_python/tin_engine/cli.py@c317c7e4:2324-2327` | The core takes only projected metres, and `catchment` has no `--out-crs` or resample step, so a DEM in degrees is refused at the core boundary. | Every geographic DEM: GLO-30 (worldwide), ANADEM (South America), so the São Francisco basin. | blocks | `catchment` runs on the resampled target grid the way `mesh --out-crs` does (`dem_input._open_reprojected`), with the same CRS suggestion on refusal (~120). |
| A2 | `src_python/tin_engine/crs.py@c317c7e4:142-151` | All DEM tiles share one CRS (`single_crs`), even with `--out-crs`. | A basin whose national tiles span two UTM zones or two countries' CRSs (the Rhine, the Danube, any border catchment). | degrades | Group footprints by CRS and resample each group onto the target grid (~150). Work-around today: reproject outside rasputin. |
| A3 | `src_python/tin_engine/target_grid.py@c317c7e4:140-143` | The source region does not cross ±180° or a pole. | Fiji, Chukotka, Kiribati, the Antarctic interior on a geographic DEM. | degrades | Documented limit, refused with a message; a per-hemisphere source region would lift it (~60). Low priority. |
| A4 | `src_python/tin_engine/target_grid.py@c317c7e4:91-101` | The resampled grid's spacing is a whole number of metres, floored at 1. | A sub-metre DEM in degrees gets a 1 m grid (probe: 0.0000045° at 50° N gives 1 m). 1" and 1/9" DEMs are fine (31 m, 3 m). | degrades | Spacing as a dyadic fraction of a metre, or `--out-spacing` (~20). |
| A5 | `src_python/tin_engine/feature_input.py@c317c7e4:16-18`, `src_python/tin_engine/feature_input.py@c317c7e4:607-613` | The fast pre-clip bound holds for a Transverse Mercator DEM only; other projections move every feature first (exact, slower). | Lambert-93 (France), LAEA Europe (Lambert azimuthal equal-area, EPSG:3035), polar stereographic: correct, slower. Southern UTM zones take the fast path (probe). | cosmetic | None needed; a bound for conic projections if the clip ever shows up as a hotspot. |
| A6 | `src_python/tin_engine/io/geotiff.py@c317c7e4:59-73` | A projected CRS given by parameters (GeoKey 3072 = 32767) is read only for Transverse Mercator and Lambert Conic Conformal 2SP. | Polar stereographic (Antarctica, Greenland), LAEA, Mercator, oblique Mercator (Switzerland, Malaysia), Albers (Alaska, NLCD) headers without an EPSG code are refused. | degrades | One table row per EPSG method (9829, 9820, 9804, 9822) with its parameter order (~20 each). |
| A7 | `src_python/tin_engine/io/geotiff.py@c317c7e4:276-277`, `src_python/tin_engine/io/geotiff.py@c317c7e4:294` | North-up rasters only: a `ModelTransformationTag` or a non-positive pixel scale is refused. | Rotated or south-up GeoTIFFs; rare in DEM distribution. | cosmetic | Documented limit. |
| A8 | `src_python/tin_engine/io/geotiff.py@c317c7e4:339-344` | A projected DEM's horizontal unit is the metre. | US State Plane in feet (EPSG:2229 and kin), refused at header read; probe: a feet `--out-crs` is refused with a clear message too. | degrades | Read feet CRSs and require `--out-crs` in metres, as for degrees (~30). |
| A9 | `src_python/tin_engine/catchment.py@c317c7e4:225-228`; `src_python/tin_engine/cli.py@c317c7e4:2361-2365` | `--rivers` must be in the DEM's CRS. | Any river file in another CRS (OpenStreetMap in EPSG:4326 over a national DEM). `--lakes` and `--seed-crs` accept any CRS. | degrades | `--rivers-crs` and a reprojection, as the lakes do (~25). |

No finding on the hemisphere itself: the probe moves points into UTM 23S
(SIRGAS 2000 and WGS 84) and 60S, `suggest_crs` returns Transverse Mercator for
the São Francisco box and an equator-crossing Amazon box, polar stereographic
with `lat_0 = -90` for the Antarctic coast; the fine outline's world frame is
`y = -row` (`src_python/tin_engine/outline.py@c317c7e4:4-6`), and no sign test
on northing exists in `src_python/` (`grep -rn 'northing' src_python/tin_engine`
hits prose and GeoKey names only). GLO-30 tile names take N/S and E/W
(`src_python/tin_engine/fetch/run.py@c317c7e4:273-275`). The suggested CRS has
no EPSG code and is written to a GeoPackage as the custom row 100000
(`src_python/tin_engine/crs.py@c317c7e4:164-186`; probe).

### B. Land-cover schemes

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| B1 | `src_python/tin_engine/feature_input.py@c317c7e4:59-66`, `src_python/tin_engine/feature_input.py@c317c7e4:106-123` | The class maps are Python objects: three, all CORINE (`Code_18`, `clc18_kode`), with CORINE's water codes and copyright notice. | Any other scheme needs a pull request: MapBiomas polygons, ESA WorldCover, NLCD, AR5 `artype`, the Austrian cadastre. The `property` map (`src_python/tin_engine/feature_input.py@c317c7e4:107-111`) takes any file already carrying vocabulary names, so edges can be had by pre-processing; no class code, palette or fraction can, hence blocks. | blocks | Class maps as data: a TOML file given to `--features-map path.toml` (name, attribute, classes, otherwise, code system, notice, palette, vocabulary), the CORINE maps shipped as package data under `tin_engine/data/` (~120). |
| B2 | `src_python/tin_engine/palettes.py@c317c7e4:15-65`; `src_python/tin_engine/cli.py@c317c7e4:1462-1476` | One palette, CORINE level 3; `rasputin palette` knows only `corine`; unknown codes are magenta. | Another scheme's mesh writes, but every class is magenta and unnamed. | degrades | The palette as a table in the class-map file (B1) (~60). |
| B3 | `src_python/tin_engine/cli.py@c317c7e4:803-842`, `src_python/tin_engine/cli.py@c317c7e4:240`; `src_python/tin_engine/feature_input.py@c317c7e4:56`, `src_python/tin_engine/feature_input.py@c317c7e4:660` | Clean-up defaults: repair 1 m, simplify 50 m, outline snap 5 m (under a 100 m cap), start angle 15°. Tuned on CORINE (minimum mapping unit 25 ha, minimum width 100 m) over a 10 m DEM. | A 1 m DEM with 1:5000 land cover: 50 m simplification removes the detail that data has; a 1 km DEM: 1 m repair is noise. | degrades | Defaults derived from the DEM cell (repair 0.1 cell, simplify 5 cells, snap 0.5 cell reproduce today's numbers at 10 m), or a `[clean_up]` block per scheme in the class-map file; the metre figures stay as the documented CORINE case (~40). Ola's ruling needed (section 7, Q1). |
| B4 | `src_python/tin_engine/land_cover_fractions.py@c317c7e4:33`, `src_python/tin_engine/land_cover_fractions.py@c317c7e4:45-73` | The class raster is in degrees on WGS 84: cell areas by `Geod(ellps="WGS84")`, index space by longitude and latitude. | True of MapBiomas (EPSG:4326). A class raster on another datum or a projected one (CORINE raster in EPSG:3035, NLCD in Albers) gets wrong cell areas. Latent: only catalogue sources reach it (B5). | degrades | Ellipsoid from `parse_crs(meta.crs).get_geod()`; constant cell area for a projected raster, or a refusal (~15). |
| B5 | `src_python/tin_engine/cli.py@c317c7e4:843-852` | `--land-cover` takes a cached catalogue source only. | A class raster on disk (a national land-cover GeoTIFF) cannot be used. | degrades | Accept a path as `--dem` does (~40). |
| B6 | `src_python/tin_engine/feature_input.py@c317c7e4:579-587` | Class codes are integers in 1 .. 2^31 - 1. | Schemes with string codes (UK LCM letters, "Forest"); AR5's `artype` is numeric and fine. | degrades | A code table in the class-map file maps strings to integers (in B1). |
| B7 | `src_python/tin_engine/io/cog.py@c317c7e4:16-17`, `src_python/tin_engine/io/cog.py@c317c7e4:167` | An empty block of a class raster is 0, "MapBiomas's convention"; the catalogue's `nodata` is not consulted. | A class source whose no-data is not 0 (ESA WorldCover uses 0 too; a future one may not). | cosmetic | Fill with the source's `nodata` (~5). |
| B8 | `src_python/tin_engine/features.py@c317c7e4:170-190` | The edge vocabulary (river, road, railway, coastline, contour, wall, ditch, land_cover, water) is "a default, not a schema", but the CLI offers no way to supply another. | A user whose features are "levee", "canal", "glacier edge" must map them onto these names. | cosmetic | The vocabulary as part of the class-map file (in B1). |

### C. DEM sources, tile layouts and the station set

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| C1 | `src_python/tin_engine/sources.py@c317c7e4:98-153`, `src_python/tin_engine/sources.py@c317c7e4:172-187` | The remote catalogue is Python code: ANADEM, GLO-30, MapBiomas, and NVE's station list. | Adding Kartverket DTM10, Lantmäteriet, USGS 3DEP, EU-DEM, FABDEM or another station service needs a pull request. Local files by path work anywhere already. | degrades | Catalogue as TOML in package data plus `$RASPUTIN_DATA/sources.toml` user entries, read into the existing Pydantic models (~80). |
| C2 | `src_python/tin_engine/fetch/run.py@c317c7e4:42-45` | GLO-30's tile planning box assumes a spacing of 0.001°, "more than any tile's below 80°". | Above 80° (Svalbard's north, Greenland) the longitude spacing is 5"-10", so blocks at a tile's edge can be under-planned and refused later as "not in the cache". | degrades | 0.003° (1 line). |
| C3 | `src_python/tin_engine/fetch/plan.py@c317c7e4:33`, `src_python/tin_engine/fetch/plan.py@c317c7e4:122` | A degree is 111 320 m, for the fetch margin of a geographic source. | East-west a degree is shorter, so the margin is larger than asked: the safe side. | cosmetic | None. |
| C4 | `src_python/tin_engine/io/repository.py@c317c7e4:74` | DEM tiles are GeoTIFFs (`.tif`, `.tiff`). | ASCII grid, netCDF, HGT, GeoPackage rasters are not read. The fetch route gives COGs. | cosmetic | Documented limit (INSTALL names it). |
| C5 | `src_python/tin_engine/data/nve_hrd_2025.csv`; `src_python/tin_engine/fetch/nve.py@c317c7e4:45-57`; `src_python/tin_engine/reference.py@c317c7e4:63-64` | `fetch-stations` and `catchment-batch` compare against NVE's Norwegian reference set; the result fields are `nve_in_ours`, `ours_in_nve`. | Norway-only by design; the `StationSource` model is general but has one entry, and the field names carry it. | cosmetic | Rename to `reference_in_ours`/`ours_in_reference`; a second station source when there is one (~30). |

Nothing in `src_python/` names DTM10 as a default (`grep -rn DTM10
src_python/tin_engine` hits comments only: `src_python/tin_engine/mosaic.py@c317c7e4:32-33`,
`src_python/tin_engine/gauge.py@c317c7e4:31`). The default set lives in the
quick check's cases (E1) and INSTALL.

### D. The machine: cores, threads, memory

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| D1 | `src_python/tin_engine/cli.py` (no `--threads` option: `grep -n -- '--threads' src_python/tin_engine/cli.py` is empty); `src_python/tin_engine/_core.pyi@c317c7e4:764-775`; `include/terrain/parallel_util/chunks.hpp@c317c7e4:45-46`; `include/terrain/land_cover/overlap.hpp@c317c7e4:273` | Every run uses every core. The C++ takes `threads` and defaults to `hardware_concurrency` (right); nothing above it lets a user set it. The only way is `tools/bench.py@c317c7e4:788-790`'s monkeypatch of `refine`. | A shared server, CI, a laptop on battery, an agent running several meshes at once, a container with a CPU quota (`hardware_concurrency` reports the host). | degrades | `--threads N` (0 = all cores) on `mesh` and `catchment`, `RASPUTIN_THREADS` as the default, carried in one `Resources` model to every `_core` call and every Python pool (~60). |
| D2 | `src_python/tin_engine/crs.py@c317c7e4:35-39` | Reprojection threads: `min(8, os.cpu_count())`, chunks of 250 000 points ("about 60 ms of PROJ work each", one machine's figure), fixed at import. | A 32-core box reprojects on 8 threads; the cap is arbitrary. | degrades | A `threads` argument to `reprojector` from D1; drop the 8 (~10). |
| D3 | `src_python/tin_engine/dem_input.py@c317c7e4:218`; `src_python/tin_engine/io/cog.py@c317c7e4:124`, `src_python/tin_engine/io/cog.py@c317c7e4:203` | Resample, check points and block decoding use `os.cpu_count()`; the parameter exists but no caller sets it. | Same as D1. | degrades | In D1. |
| D4 | `include/terrain/parallel_util/chunks.hpp@c317c7e4:36`, `include/terrain/parallel_util/chunks.hpp@c317c7e4:39-41`; `include/terrain/refinement/check_points.hpp@c317c7e4:56-57`; `include/terrain/refinement/refine.hpp@c317c7e4:310`; `src_python/tin_engine/land_cover_fractions.py@c317c7e4:31` | Block sizes tuned on one machine: `inline_below` 256 (from the 1 m benchmark), blocks of n / (16 threads), 16 KiB chunks and 64 MiB slabs, 4096 triangles, 2^18 triangles per class window. | Another machine may want other sizes; results do not depend on them. | cosmetic | None; each states its scale. |
| D5 | `src_python/tin_engine/cli.py@c317c7e4:1560`; `src_python/tin_engine/stats.py@c317c7e4:241-242` | The stats file prints `os.cpu_count()` as "Threads: N (hardware concurrency)". | True only because nothing can ask for fewer; a benchmark child running on 1 thread reports every core. | cosmetic | Print the count used (in D1). |
| D6 | `src_python/tin_engine/catchment.py@c317c7e4:263`; `src_python/tin_engine/mosaic.py@c317c7e4:55-57` | The flood window may use half the physical memory, read by `os.sysconf`; not overridable. | A container with a cgroup limit reads the host's memory: the cap is too high and the kernel kills the run instead of a refusal. | degrades | `--memory-budget` overriding the cap; read `/sys/fs/cgroup/memory.max` when present (~25). |
| D7 | `src_python/tin_engine/decompose.py@c317c7e4:19-20` | A 16 GiB default budget and a bytes-per-node table from one machine's sweep. `partition` is called by nothing in `src_python/` or `tools/` (`grep -rn decompose src_python tools` finds the module and one docstring, `src_python/tin_engine/io/mesh_index.py@c317c7e4:37`; its only importer is `tests/python/test_decompose_partition.py`). | Dead module carrying a machine constant. | cosmetic | Delete, or wire into D6 (-80 or +30). |
| D8 | `docs/increments/36-geopackage-output.md@c317c7e4:343-344` | Acceptance: peak RSS under 1 GB and 3 GB on the two test meshes. | Acceptance figures for Ola's machine, in docs, not code. | cosmetic | None; they are labelled as acceptance. |

### E. Fixed catchments and fixtures that double as defaults

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| E1 | `docs/benchmarks/quick/cases.toml@c317c7e4:24-78`; `tools/bench_quick.py@c317c7e4:205`, `tools/bench_quick.py@c317c7e4:213-215`; `docs/benchmarks/quick/hotspots.toml@c317c7e4:7-26` | Four of six quick-check cases read `$RASPUTIN_DATA` (Numedalslågen, the Bergen Line and Romsdalen the DTM10 set; Lagan a GLO-30 cache, Swedish SVAR (SMHI's water archive) catchments and CORINE); a missing input is a refusal, not a skip; the hotspot rulings are per case. | The quick check cannot run on another machine. It is a benchmark fixture, not product. | cosmetic | Skip absent cases; one public-data case (GLO-30 over the Velhas fixture) so another machine has a check (~10). |
| E2 | `tools/bench.py@c317c7e4:45-46` | The 1 m benchmark's DEM is the committed DTM10 tile 7908_3 and its quarter. | A committed fixture; fine. | cosmetic | None. |
| E3 | `tests/fixtures/dem_archive/7908_3_10m_z33.tif` (EPSG:25833); `tests/fixtures/velhas/anadem_velhas.tif` (EPSG:4674, southern, used by `tests/python/test_cli_mesh_geographic.py`) | The suites cover one northern projected DEM and one southern geographic DEM. | No fixture for a southern projected DEM, a conic one, or a feet one; the probes in this audit stand in. | cosmetic | Fixtures and tests only (0 production lines). |
| E4 | `src_python/tin_engine/gauge.py@c317c7e4:29-38`; `src_python/tin_engine/burn.py@c317c7e4:25-31`; `src_python/tin_engine/catchment.py@c317c7e4:63`; `src_python/tin_engine/cli.py@c317c7e4:2419-2420` | Gauge and burn constants in metres: U floor 30 ("three DTM10 cells"), lake gap 30, confluence 100, down extra 100, end cap 500, window margin 2000, map radius 500, reach up 1000. | On GLO-30 (30 m) the floor is one cell; on a 1 m DEM it is thirty. Tuned on NVE rivers over DTM10. | degrades | Express in cells where a cell is meant (the gauge module is pure shapely, so the cell size is passed in), keep the metre figures as the 10 m case (~15). |

### F. Hand-tuned constants in the algorithms

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| F1 | `src_python/tin_engine/cli.py@c317c7e4:225-236` | Snap spacing 1 mm "in the input's units", a millimetre on the ground at UTM magnitudes; scale and overflow bound stated. | Correct for any metre CRS; a parameter. | cosmetic | None. |
| F2 | `src_python/tin_engine/cli.py@c317c7e4:238-243` | Start angles 25° and 15° (with land cover), gain 0: parameters with ruled defaults. | Scale-free. | cosmetic | None. |
| F3 | `src_python/tin_engine/grid_domain.py@c317c7e4:23`, `src_python/tin_engine/grid_domain.py@c317c7e4:35` | 256 and 129 nodes per side for the start grids. | Node counts, scale-free; `--stride` overrides. | cosmetic | None. |
| F4 | `src_python/tin_engine/mosaic.py@c317c7e4:32-38` | Seam threshold 1 mm "DEM units (metres for DTM10)"; alignment 1e-6 cells. | Vertical metres assumed (H1); alignment scale-free. | cosmetic | Follows H1. |
| F5 | `include/terrain/mesh/lawson.hpp@c317c7e4:104`; `include/terrain/raster/steepness.hpp@c317c7e4:37` | Exact lattice incircle within 2^14 cells, else the filtered kernel; 181 half-degree slope classes. | Derived quantities with a fallback; scale-free. | cosmetic | None. |
| F6 | `src_python/tin_engine/feature_input.py@c317c7e4:54-56` | Features are read from the domain plus 100 m, densified at 1 km. | On a 1 km DEM (GMTED) 100 m is a tenth of a cell: a feature just outside the domain is dropped before the DEM's own one-cell grow. | degrades | `max(100 m, 2 cells)` (~2). |
| F7 | `src_python/tin_engine/reference.py@c317c7e4:28-33` | Match 95 %, close 80 %, offset 3 cells. | Fractions and cells: scale-free. | cosmetic | None. |

### G. Paths and naming

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| G1 | `src_python/tin_engine/cli.py@c317c7e4:1727-1733`, `src_python/tin_engine/cli.py@c317c7e4:1376`, `src_python/tin_engine/cli.py@c317c7e4:1649`, `src_python/tin_engine/cli.py@c317c7e4:1752` | The tile cache is `--cache` or `$RASPUTIN_DATA/cache`; without either, a refusal. | A first run anywhere (and an agent through the MCP server) must set an environment variable before `fetch` works. | degrades | Default to `$XDG_CACHE_HOME/rasputin` (`~/.cache/rasputin`), `RASPUTIN_DATA` kept as the override (~10). |
| G2 | `tools/bench_quick.py@c317c7e4:205`; `INSTALL.md@c317c7e4:160`, `INSTALL.md@c317c7e4:234` | Data lives in `../rasputin_data` beside the checkout. | A convention of Ola's disk, stated in INSTALL; tooling only. | cosmetic | None; the ideas file's "Structure for rasputin_data" covers it. |
| G3 | `src_python/tin_engine/hydrography.py@c317c7e4:53`; `src_python/tin_engine/sources.py@c317c7e4:179-183` | `vatnlnr`, `navn`, "Kilde: NVE", "Totalnedbørfelt til målestasjon". | These are NVE's field names and the licence's own words, read from NVE's service: correct, not a defect. | cosmetic | None (C5 covers the result-field names). |

### H. Locale and units

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| H1 | `src_python/tin_engine/io/geotiff.py@c317c7e4:175-178`, `src_python/tin_engine/io/geotiff.py@c317c7e4:193`; `src_python/tin_engine/cli.py@c317c7e4:1120` | Heights are metres: a `VerticalUnitsGeoKey` other than metre is refused, an absent one is assumed metres and reported as assumed. | A DEM in feet without the key (US county lidar, some state products): `--tolerance 1` reads as 1 foot, silently 3.28 times off. | degrades | `--vertical-unit {m,ft,us-ft}` scaling at decode, and `VerticalCSTypeGeoKey` (4096: EPSG 5703 metres, 6360 feet) read when present (~30). |
| H2 | `src_python/tin_engine/dem_input.py@c317c7e4:164-172`; `src_python/tin_engine/cli.py@c317c7e4:682-689` | The mesh is computed in a projected CRS in metres; `--tolerance` is metres. | A rule, applied everywhere, with a refusal naming the unit (A8). | cosmetic | Documented limit. |

Locale: no finding. `grep -rn -E 'locale|setlocale|strftime|%d\.%m|decimal' src_python/tin_engine`
hits nothing that formats by locale; numbers go through Python's
locale-independent formatting, dates through `isoformat`
(`src_python/tin_engine/fetch/nve.py@c317c7e4:211`, `src_python/tin_engine/fetch/run.py@c317c7e4:157`).

### I. macOS and Apple silicon

| # | Where | What is assumed | Where it breaks | Severity | Fix shape (lines) |
|---|---|---|---|---|---|
| I1 | `tools/bench.py@c317c7e4:51`, `tools/bench.py@c317c7e4:587-597`; `tools/bench_quick.py@c317c7e4:217`, `tools/bench_quick.py@c317c7e4:267` | `pmset -g batt`, `sysctl machdep.cpu.brand_string`, `hw.perflevel*`, `sw_vers`: the benchmark and the quick check run on macOS only; the acceptance rule needs the power state. | No acceptance run on Linux; CI cannot run the benchmark. | degrades | A machine-and-power provider per platform (`/proc/cpuinfo`, `/sys/class/power_supply`) (~60, tooling). |
| I2 | `src_python/tin_engine/mosaic.py@c317c7e4:55-57`; `pyproject.toml@c317c7e4:18-19` | `os.sysconf` for physical memory: POSIX. Linux and macOS are the declared platforms. | Windows is not claimed. | cosmetic | Documented limit. |
| I3 | `docs/increments/27-node-sampling.md@c317c7e4:202-216`, `docs/increments/27-node-sampling.md@c317c7e4:388-394`; `CMakeLists.txt` (no `-ffp-contract` flag: `grep -n contract CMakeLists.txt` is empty) | Apple clang fuses `x_min + col * delta_x` into a multiply-add; NumPy does not; GCC on x86-64 without `-mfma` does not. On a non-integer, non-dyadic grid origin or spacing the C++ and Python node coordinates differ in the last bit, and so can the mesh between a Mac and a Linux build. | Cross-platform bit reproducibility on such grids; invisible on integer grids (DTM10, GLO-30's target grid). | degrades | `-ffp-contract=off` on the targets, timed by `@perf` (2 lines of CMake); increment 27 names this as its own increment. |
| I4 | `.github/workflows/main.yaml@c317c7e4:56` | CI builds and tests on `ubuntu-latest` and `macos-latest`. | Positive: the product's gates run on Linux; only the perf tooling is Mac-only. | cosmetic | None. |

## 3. What is already general

So the reader knows where the boundary is, the things the sweep found right:

- The C++ core holds no CRS, no unit and no path; "metres" appears in four
  comments (`grep -rn metres include/terrain`:
  `include/terrain/refinement/refine.hpp@c317c7e4:100`,
  `include/terrain/refinement/refine_points.hpp@c317c7e4:109`,
  `include/terrain/refinement/line_tolerance.hpp@c317c7e4:45`,
  `include/terrain/refinement/slope_tolerance.hpp@c317c7e4:43`) and in no
  computation. Every parallel loop takes `threads`, 0 meaning the machine
  (`include/terrain/parallel_util/chunks.hpp@c317c7e4:17`).
- One reprojection site, `always_xy` (`src_python/tin_engine/crs.py@c317c7e4:3-6`);
  `suggest_crs` chooses polar stereographic, Transverse Mercator, Mercator or
  Lambert conic by latitude and shape on the DEM's own datum
  (`src_python/tin_engine/crs.py@c317c7e4:202-248`); probe: São Francisco
  0.23 % scale error, Svalbard 0.44 %, Antarctic coast 0.69 %.
- A geographic DEM is refused with the CRS to paste
  (`src_python/tin_engine/dem_input.py@c317c7e4:174-189`), and the resampled
  grid's spacing comes from the geodesic at the domain
  (`src_python/tin_engine/target_grid.py@c317c7e4:91-101`): 31 m for 1" at
  -15°, 0° and 80° (probe).
- `catchment --seed-crs` and `--lakes` take any CRS; the outline tolerance
  defaults to two cells (`src_python/tin_engine/catchment.py@c317c7e4:245-247`).
- The southern-hemisphere fixture exists (`tests/fixtures/velhas/`, ANADEM,
  EPSG:4674) and the geographic path is tested on it.

## 4. Prior art: legacy and literature

**Legacy.** `git grep -l -i -E 'utm|hemisphere|\bsouth\b' legacy-archive -- legacy`
(tests and Markdown dropped) returns one file, `legacy/rasputin/reader.py`: it
decoded Kartverket's private 160zz/161zz UTM codes into a zone and a south
flag (`legacy/rasputin/reader.py:252-262`), hemisphere-aware but Norway-specific,
and `legacy/rasputin/application.py:40` defaulted the target to `EPSG:32633`,
so the legacy hardcoded UTM 33N too. Ours reads the CRS through GeoKeys and
PROJ, which is the general form; nothing carried.
`git grep -l -i -E 'cpu_count|hardware_concurrency|omp_get|num_threads' legacy-archive -- legacy`
returns nothing: the legacy was serial, so no thread setting to inherit.
`git grep -l -i -E 'land_cover|landcover|globcov|corine' legacy-archive -- legacy`
returns `legacy/rasputin/application.py`, `globcov_repository.py`,
`gml_repository.py`, `land_cover_repository.py`, `tin_repository.py`,
`web_visualize.py`, `wfs_repository.py`: the legacy had an abstract
`LandCoverRepository` with a per-scheme `LandCoverBaseType` enum and `color`,
`describe`, `material` methods (`legacy/rasputin/land_cover_repository.py:8-46`),
implemented for GlobCover 2009 (global, EPSG:4326) and for a Norwegian WFS. That
is the scheme-as-swappable-component shape B1 proposes, in data instead of
classes. Nothing carried as code.

**Literature.** No novelty is claimed anywhere in this audit, so there is
nothing to check; the general rules the fixes lean on are standard practice:
the conformal-family choice by latitude band and extent that `suggest_crs`
already follows (Snyder 1987, *Map Projections: A Working Manual*, USGS
Professional Paper 1395); CORINE's minimum mapping unit of 25 ha and minimum
width of 100 m (EEA, *CLC2018 Technical Guidelines*, 2017), which is where the
50 m default comes from and why a general rule derives the clean-up distances
from the source's minimum width or the DEM cell; and "hardware concurrency as
the default, a user override, container limits honoured" for threads and
memory, as `OMP_NUM_THREADS`, `make -j` and the JDK's cgroup awareness do.

## 5. The ideas file

Entries this audit takes up or overlaps: *Catchments on a geographic DEM
(GLO-30)* (A1); *Product review (b), generality audit* (this file); *Structure
for rasputin_data and rasputin_scratch* (G1, G2); *Release hygiene with
priority: DTM10 download* (C1); *Repair CORINE once per area of interest* (B3:
the cached extract carries the clean-up distances used); *MCP server* (D1 and
G1: an agent needs a thread setting and a cache that exists without an
environment variable). Ideas this audit leaves out are listed in the handback.

## 6. Fixes in order, and a PR split

What I would ship first, each under 700 net lines:

1. **Threads and memory as settings** (D1, D2, D3, D5, D6, D7): `--threads`,
   `RASPUTIN_THREADS`, one `Resources` model through every `_core` call and
   Python pool, the stats row fixed, `--memory-budget` on `catchment` with the
   cgroup limit read, `decompose.py` deleted. About 130 net lines (about 210
   added, 80 removed). Ola's first finding.
2. **Catchment on a geographic DEM** (A1, A9): the resampled target grid route
   with the CRS suggestion, `--rivers-crs`. About 150 lines; hydrology, not
   refine or mesh, so no `@perf` acceptance by the rule, but a timing.
3. **Land-cover schemes as data** (B1, B2, B6, B8, and B3 if Ola rules for
   derived defaults): class-map files with palette and vocabulary,
   `--features-map path.toml`, the CORINE maps as package data. About 250 lines.
4. **Catalogue as data, class rasters by path** (C1, B5, B4, B7, C2): a TOML
   catalogue plus user entries, `--land-cover path.tif`, the ellipsoid from the
   raster's CRS. About 160 lines.
5. **Units and header coverage** (H1, A8, A6, A4, F6): vertical unit option,
   feet CRSs, the four projection methods, dyadic spacing, the read margin in
   cells. About 120 lines.
6. **Cache default, names, tooling on Linux** (G1, C5, I1, E1): the XDG cache
   default, `reference_in_ours`, a Linux machine-and-power provider, a public
   quick-check case. About 120 lines.
7. **Tiles in two CRSs** (A2), its own PR, about 150 lines; and
   **`-ffp-contract=off`** (I3), its own two-line PR with `@perf`.

## 7. Questions for Ola

Each with the default this audit assumes.

- **Q1** Should the land-cover clean-up defaults (1 m, 50 m, 5 m, 15°) be
  derived from the DEM cell, or stay absolute with the CORINE case documented?
  Default: derive (0.1, 5 and 0.5 cells reproduce today's numbers at 10 m); the
  start angle stays 15°.
- **Q2** Class maps and palettes as TOML files given by path, with the three
  CORINE maps shipped as package data and still reachable by name? Default: yes.
- **Q3** `--threads N`, 0 meaning every core, with `RASPUTIN_THREADS` as the
  default source? Default: yes.
- **Q4** Delete `decompose.py` (dead, with its 16 GiB budget) or wire it into
  the catchment's memory cap? Default: delete; the cap gets `--memory-budget`.
- **Q5** For `catchment` on a geographic DEM, resample onto the suggested
  projected grid as `mesh --out-crs` does, or flood in degrees? Default:
  resample.
- **Q6** Order: default as in section 6, threads first.

## Appendix: the probe

`probe_crs.py`, run as
`<python with pyproj 3.8.0, PROJ 9.8.1, shapely, numpy, pydantic> probe_crs.py <tree>/src_python`
with `PROJ_NETWORK=OFF`. The `sys.meta_path` line matters: an interpreter with
an editable install of another checkout imports that checkout's package
otherwise, and the first run of this probe did.

```python
import os, sys
os.environ["PROJ_NETWORK"] = "OFF"
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, sys.argv[1])
import tin_engine
assert tin_engine.__file__.startswith(sys.argv[1]), tin_engine.__file__
import numpy as np
from pyproj import CRS
from tin_engine import crs as C
from tin_engine.io.models import RasterMeta
from tin_engine.target_grid import default_spacing
from tin_engine.dem_input import DemRequest, OutCrsRequiredError, _reprojected
from pathlib import Path

for name, box, geo in [
    ("Sao Francisco (ANADEM datum)", (-47.0, -36.0, -21.0, -7.0), "EPSG:4674"),
    ("equator-crossing Amazon", (-62.0, -58.0, -2.0, 2.0), "EPSG:4326"),
    ("Svalbard", (10.0, 30.0, 76.0, 81.0), "EPSG:4326"),
    ("Antarctic coast", (160.0, 170.0, -78.0, -72.0), "EPSG:4326"),
]:
    s = C.suggest_crs(box, geo)
    print(f"{name:32s} {s.family:26s} scale {100*s.max_scale_error:.3g} %  areal {100*s.max_areal_error:.3g} %  {s.proj4[:60]}")
for dst in ("EPSG:31983", "EPSG:32723", "EPSG:32760"):
    f = C.reprojector("EPSG:4326", dst)
    print(dst, CRS(dst).name, f(np.array([[-44.0, -15.0], [-44.0, -5.0], [200.0, 95.0]])).round(2).tolist())
for code in ("EPSG:31983", "EPSG:32723", "EPSG:25833", "EPSG:3006", "EPSG:2154", "EPSG:3035", "EPSG:3031"):
    op = CRS(code).coordinate_operation
    names = {p.name for p in op.params} if op else set()
    ok = op is not None and op.method_name == "Transverse Mercator" and {"Longitude of natural origin", "Scale factor at natural origin"} <= names
    print(f"{code:11s} method={op.method_name if op else None!s:36s} pre-clip bound applies: {ok}")

def meta(dy, epsg=4326):
    return RasterMeta(x_min=0.0, y_max=0.0, delta_x=dy, delta_y=dy, cols=2, rows=2, epsg=epsg,
                      geographic=True, nodata=None, nodata_source="absent", pixel_is_area=False,
                      vertical_unit_assumed=True)
for lat in (-15.0, 0.0, 60.0, 80.0):
    print("1 arcsec at lat", lat, "->", default_spacing(meta(1/3600), (-44.0, lat)), "m")
print("1/9 arcsec at lat 40 ->", default_spacing(meta(1/32400), (-100.0, 40.0)), "m")
print("0.0000045 deg at lat 50 ->", default_spacing(meta(0.0000045), (10.0, 50.0)), "m")
try:
    _reprojected(DemRequest(sources=(Path("x.tif"),), target_crs="EPSG:2229"), [meta(1/3600)])
except ValueError as e:
    print("feet --out-crs:", str(e)[:120])
m = RasterMeta(x_min=-44.0, y_max=-14.0, delta_x=1/3600, delta_y=1/3600, cols=3600, rows=3600, epsg=4674,
               geographic=True, nodata=None, nodata_source="absent", pixel_is_area=False, vertical_unit_assumed=True)
try:
    _reprojected(DemRequest(sources=(Path("x.tif"),)), [m])
except OutCrsRequiredError as e:
    print("geographic DEM, no --out-crs:", e.head[:160])
s = C.suggest_crs((-47.0, -36.0, -21.0, -7.0), "EPSG:4674")
print("label of suggested:", C.crs_label(s.proj)[:60], "| srs_id", C.srs_entry(s.proj).srs_id)
```

Output on `c317c7e4` (2026-10-12):

```
Sao Francisco (ANADEM datum)     transverse Mercator        scale 0.229 %  areal 0.458 %  +proj=tmerc +lat_0=0 +lon_0=-41.5 +k=0.997712 +x_0=0 +y_0=0
equator-crossing Amazon          transverse Mercator        scale 0.0307 %  areal 0.0614 %  +proj=tmerc +lat_0=0 +lon_0=-60.0 +k=0.999693 +x_0=0 +y_0=0
Svalbard                         polar stereographic        scale 0.439 %  areal 0.881 %  +proj=stere +lat_0=90.0 +lon_0=20.0 +k_0=0.989479 +x_0=0 +y_
Antarctic coast                  polar stereographic        scale 0.689 %  areal 1.38 %  +proj=stere +lat_0=-90.0 +lon_0=165.0 +k_0=0.982259 +x_0=0 +
EPSG:31983 SIRGAS 2000 / UTM zone 23S [[607512.23, 8341431.16], [610859.93, 9447251.38], [inf, inf]]
EPSG:32723 WGS 84 / UTM zone 23S [[607512.23, 8341431.16], [610859.93, 9447251.38], [inf, inf]]
EPSG:32760 WGS 84 / UTM zone 60S [[5270117.24, -7830071.12], [5486390.72, -9263196.22], [inf, inf]]
EPSG:31983  method=Transverse Mercator                  pre-clip bound applies: True
EPSG:32723  method=Transverse Mercator                  pre-clip bound applies: True
EPSG:25833  method=Transverse Mercator                  pre-clip bound applies: True
EPSG:3006   method=Transverse Mercator                  pre-clip bound applies: True
EPSG:2154   method=Lambert Conic Conformal (2SP)        pre-clip bound applies: False
EPSG:3035   method=Lambert Azimuthal Equal Area         pre-clip bound applies: False
EPSG:3031   method=Polar Stereographic (variant B)      pre-clip bound applies: False
1 arcsec at lat -15.0 -> 31 m
1 arcsec at lat 0.0 -> 31 m
1 arcsec at lat 60.0 -> 31 m
1 arcsec at lat 80.0 -> 31 m
1/9 arcsec at lat 40 -> 3 m
0.0000045 deg at lat 50 -> 1 m
feet --out-crs: --out-crs EPSG:2229 is a Projected CRS with axes in ['US survey foot']; the mesh is computed in it, so it must be a proj
geographic DEM, no --out-crs: the DEM is geographic (EPSG:4674), so --out-crs is required; suggested for this box (transverse Mercator, +proj=tmerc +lat_0=0 +lon_0=-43.5 +k=0.999982 +x_0=0 +
label of suggested: PROJCRS["undefined",BASEGEOGCRS["SIRGAS 2000",DATUM["Sistema | srs_id 100000
```

## Review

Round 1, 2026-10-12, @reviewer (`c317c7e4..6d468ba6`, 0 production lines): CHANGES REQUESTED, three wordings, no finding wrong in substance; applied in this recording commit as specified: /Users/skavhaug/projects/rasputin/.claude/worktrees/generality-audit/docs/increments/generality-audit.md@6d468ba6:154 (E1: Lagan reads a GLO-30 cache, SVAR and CORINE, not DTM10), @6d468ba6:147 (D7: the grep that finds `decompose`'s one docstring and its only importer, a test), @6d468ba6:113 (B1: why the `property` map's work-around does not lower "blocks" to "degrades"). Applied with them, the reviewer's human-reader idea: NLCD, LAEA, SVAR and AR5 expanded on first use. Verified by running: `check_citations.py` green and every one of the about 110 pinned ranges dumped from `c317c7e4` and re-read as a quotation; the six section-1 greps reproduce; the appendix probe re-run on this tree with `tin_engine.__file__` asserted gives the saved 23 lines verbatim; A1's refusal reproduces on `tests/fixtures/velhas/anadem_velhas.tif`, and `catchment` has no `--out-crs`; the severity tally 2/23/25 and the seven-PR split hold, nothing rated lower blocks the São Francisco goal; a spot-check of `io/rivers.py`, `io/station_set.py`, `io/geojson.py`, `grow.py` and `sensitivity.py` found nothing of consequence missed. Ideas left for `@architect` at the first fix: list `src_python/tin_engine/reference.py@c317c7e4:30` under C; a cosmetic row for `src_python/tin_engine/io/rivers.py@c317c7e4:24` (a 1 cm vertex merge in metres). Lesson: a stale `_core` in the main venv blocks CLI runs; the Python-side refusal is reachable through `raster.to_core` directly. Pushed on the reviewer's wording without a second round, as the recording-commit practice allows.
