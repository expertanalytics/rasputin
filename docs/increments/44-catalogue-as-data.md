# Increment 44: the source catalogue as data, class rasters by path

**Status: designed; design review round 1 applied, not yet re-reviewed.** By
`@architect`, 2026-10-13, on `b1353075` (branch `worktree-catalogue-data`),
while Ola was away; the defaults taken are section 9's questions. Python only,
one PR.

## 1. What Ola asked, and what this fixes

Ola, 2026-10-13, answered "Yes" to queuing the generality audit's fixes 4 and
then 6 after the docstring sweeps (as the main session reported it). Fix 4 is
"Catalogue as data, class rasters by path": audit rows C1, B5, B4, B7 and C2
(`docs/increments/generality-audit.md`, sections 2 and 6).

- **C1, the catalogue is Python code.** ANADEM, GLO-30, MapBiomas and NVE's
  station list are Python objects
  (`src_python/tin_engine/sources.py@b1353075:98-156`,
  `src_python/tin_engine/sources.py@b1353075:175-190`). Adding a source
  (Kartverket, USGS 3DEP, FABDEM, ESA WorldCover) takes a pull request.
- **B5, `--land-cover` takes only a cached catalogue source**
  (`src_python/tin_engine/cli.py@b1353075:1477-1512`). A land-cover GeoTIFF on
  disk cannot be used.
- **B4, the class raster is assumed to be in degrees on WGS 84**
  (`src_python/tin_engine/land_cover_fractions.py@b1353075:33`,
  `src_python/tin_engine/land_cover_fractions.py@b1353075:45-57`). Measured
  (appendix, P3): on a projected class raster (EPSG:3035, 100 m cells) every
  row area is NaN and `assign` dies with an `IndexError`. On another ellipsoid
  the areas are slightly wrong, but fractions are ratios and the effect is
  below float32 (P4: the row-to-row ratio over ten rows differs by 7e-11
  between WGS 84 and International 1924). The projected case is the defect;
  latent until B5 opens the path.
- **B7, an empty block of a class raster is 0** and the source's no-data is
  never consulted (`src_python/tin_engine/io/cog.py@b1353075:15-17`,
  `src_python/tin_engine/io/cog.py@b1353075:168`). Measured (P2): a raster
  tagged no-data 255 gives 255 as a class in the fractions, while its empty
  blocks give 0.
- **C2, GLO-30's tile planning spacing, 0.001°**
  (`src_python/tin_engine/fetch/run.py@b1353075:42-45`,
  `src_python/tin_engine/fetch/run.py@b1353075:258-262`). Measured (P6,
  against the Copernicus DEM Product Handbook's Table 3): the box is too small
  only in the 85° to 90° band, where the longitude spacing is 10"; the audit's
  "above 80°" overstates it, since the 80° to 85° band (5") is covered. A
  domain near 87° S whose east edge is 0.0063° short of a tile edge gets one
  tile where the tile's own spacing needs two. The tile naming is GLO-30's in
  code as well (`src_python/tin_engine/fetch/run.py@b1353075:273-275`).

What the user sees after this increment:

```bash
rasputin fetch glo30 --domain d.geojson                 # as today
rasputin fetch my-dem --domain d.geojson                # an entry of $RASPUTIN_DATA/sources.toml
rasputin mesh --dem dem.tif --land-cover mapbiomas-c11  # as today, same bytes
rasputin mesh --dem dem.tif --land-cover clc2018_100m.tif --land-cover-map ./clc-grid.toml
```

## 2. Prior art: legacy and literature

*Legacy.* `git grep -l -i -E 'catalog|url_template|tile_list|\.toml|\.yaml|\.ini\b|configparser' legacy-archive -- legacy`
returns nothing: the legacy had no source catalogue. `git grep -l -i -E 'globcov|land_cover_repository|nodata|no_data' legacy-archive -- legacy`
(tests and Markdown dropped) returns `legacy/rasputin/application.py`,
`globcov_repository.py`, `gml_repository.py`, `tin_repository.py`. The legacy
read a class raster by path: `GlobCovRepository(path=...)` opened
`GLOBCOVER_L4_200901_200912_V2.3.tif` from a directory, with the CRS
hard-coded as EPSG:4326, sampled one class per vertex (not area fractions), and
called class 230 "no data" in its enum, not in a tag. Carried: the intent, a
class raster given by path (B5). Not carried: the fixed CRS (this design reads
it from the file) and per-vertex sampling (26b's exact fractions stay).
GlobCover's 230 shows that a scheme's no-data can be a class value with no
tag; section 4.6 keeps that case working as a class, named by the class map.

*Literature.* No novelty is claimed. The catalogue is a small, fetch-specific
cousin of the SpatioTemporal Asset Catalog (STAC 1.0.0, <https://stacspec.org>,
the OGC community standard for describing geospatial assets) and of a file
registry such as `pooch`'s. It departs from STAC on purpose: the fields here
are the ones the fetcher executes (`kind`, URL template, tile list, tile name,
planning spacing) and the attribution the mesh carries (credit, licence,
citations), which STAC holds only in extensions; a STAC reader would add a
schema of hundreds of fields for six we use. Cell areas on the ellipsoid are
GeographicLib's geodesic polygon areas (Karney, C. F. F. (2013), *Algorithms
for geodesics*, J. Geodesy 87, 43-55), through `pyproj.Geod`, as today. The
GLO-30 spacings are the Copernicus DEM Product Handbook, issue 5.0
(GEO.2018-1988-2, 2022-11-29), Table 3, page 15: latitude 1" everywhere,
longitude 1", 1.5", 2", 3", 5" and 10" in the bands from 0°-50° to 85°-90°.

## 3. The ideas file

Entries this design takes up or overlaps: *Product review (b), generality
audit* (this is its fix 4); *Release hygiene with priority: ... DTM10 download
instructions* (a user entry is how a Kartverket source would be added, if
Geonorge serves COGs; this design does not check that); *Structure for
rasputin_data and rasputin_scratch* (the user file sits at the data root,
`$RASPUTIN_DATA/sources.toml`; a restructure moves it with the rest); *MCP
server over the rasputin commands* (an agent can add a source without a pull
request). Left out and handed back: a `rasputin sources` command listing the
merged catalogue; tile layouts other than GLO-30's 1° grid named by the
south-west corner (USGS 3DEP names tiles by the north-west corner, lower case).

## 4. The design

### 4.1 Data flow

```
tin_engine/data/sources.toml ─┐
$RASPUTIN_DATA/sources.toml ──┴─> sources.catalogue(user_file) -> Catalogue
                                     │ (pure: text in, frozen models out; refusals name the file)
cli: data_root(environ) ─────────────┤
                                     ├─> fetch / fetch-stations / mesh --dem KEY / mesh --land-cover KEY
mesh --land-cover VALUE ─> key in catalogue? ── yes: CacheRepository(root, id, nodata=entry.nodata, z_unit="metre")
                                             └─ no:  TiffDemRepository([path], z_unit="metre")  (a file; refused if absent)
                       both satisfy land_cover_fractions.ClassRaster
                       ─> assign(): index_space (unchanged), row_areas(window, crs) (B4),
                          repository.load_class_window -> decode_class_window (B7)
```

The C++ core is untouched: `_core.class_areas` takes the same corners, classes
and row areas. No path, CRS or TOML crosses into it.

### 4.2 The catalogue file (C1)

One TOML file in the package, `src_python/tin_engine/data/sources.toml`, holds
the four shipped entries, each a table keyed by its id; the Python data in
`sources.py` goes. The table's key is the id; an `id` key inside a table is
refused. The fields are `RemoteSource`'s and `StationSource`'s as today, plus
two for tiled sources (4.4). The comments on the GLO-30 licence articles and
MapBiomas's credit format move into the TOML file as comments.

```toml
[source.glo30]
kind = "cog-tiles"
url_template = "https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"
tile_list_url = "https://copernicus-dem-30m.s3.amazonaws.com/tileList.txt"
tile_name = "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM"
tile_spacing = [0.0027778, 0.00027778]   # largest cell (lon, lat) of any tile, degrees, rounded up
crs = "EPSG:4326"
z_unit = "metre"
credit = "..."
licence_note = "..."

[station.nve-hrd]
service_url = "https://kart.nve.no/enterprise/rest/services"
...
```

Model changes, all in `sources.py`:

- `RemoteSource` and `StationSource` get `extra="forbid"`: a misspelt key in a
  user file is refused, not ignored (Pydantic ignores unknown keys by default).
- `RemoteSource.nodata` gets the default `None`: TOML has no null, so an
  absent key is "no no-data stated" (GLO-30's case).
- Format fields are checked with `string.Formatter().parse`: `credit` may use
  `{date}` and `{url}` only; `url_template` `{year}` for `one-cog-per-year`,
  `{name}` for `cog-tiles`, none for `one-cog`; `tile_name` `{ns}`, `{ew}`,
  `{lat}`, `{lon}` only. A dotted or indexed field (`{date.year}`,
  `{url[0]}`) is refused. This keeps a user file from reaching object
  attributes through `str.format`, and catches a template typo at load rather
  than at fetch.
- A `z_unit` on a `values = "classes"` entry is refused: class numbers are
  not heights. `RemoteSource.read_unit` is the unit every reader of the
  source's files is given: `"metre"` for classes, `z_unit` for heights
  (4.7). `fetch/run.py`'s `parse_prefix` call and the cached land-cover route
  pass it in place of `source.z_unit`.

### 4.3 Loading, and the user's entries

`SOURCES` and `STATION_SOURCES` stay module constants, now read from the
package file at import (`importlib.resources`, as the NVE station list and the
class maps are). Every command that takes a catalogue key calls
`catalogue(user_file)` once, which returns the shipped entries plus the user
file's:

- `cli.data_root(environ)` becomes the one reader of `RASPUTIN_DATA`;
  `cache_root` and the user file (`<data root>/sources.toml`) both derive from
  it. No new option or variable.
- No `RASPUTIN_DATA`, or no file there: the shipped entries alone, no message.
  `--cache` alone, without `RASPUTIN_DATA`, loads no user entries: the cache
  directory says where tiles are, not where the user file is.
- The user file holds `[source.<id>]` tables only. A `[station.<id>]` table is
  refused: fetching a station list follows NVE's map-service layout, which is
  code (`fetch/nve.py`), so a second station service is a pull request anyway.
- A user id equal to a shipped one is refused, naming both files (section 9,
  Q2).
- A user entry's `class_map`, when it is a relative path, is taken relative
  to the user file's directory.
- `catalogue` reads `SOURCES` from the module at call time, not at import, so
  the test suites that put an entry into `SOURCES` with `monkeypatch`
  (`test_cli_fetch.py`, `test_cli_mesh_cache.py`, `test_cli_mesh_geographic.py`,
  `test_cli_units.py`, `test_cli_mapbiomas.py`) keep working unchanged.
- Each refusal is a `ValueError` naming the file and the key path
  (`sources.toml: source.my-dem.crs: Field required`); the CLI turns it into a
  usage error. The wording comes from 41's reader, shared (4.8).
- Help texts that list sources list the shipped ones and add "or an entry of
  $RASPUTIN_DATA/sources.toml"; they are built at import and cannot see the
  user file.

### 4.4 Tiled sources as data (C2)

Two fields, required for `kind = "cog-tiles"` and refused for the other kinds:

- `tile_name`: the object name of the 1° tile whose south-west corner is
  `(lat, lon)`; `{ns}` is `N` or `S`, `{ew}` `E` or `W`, `{lat}` and `{lon}`
  the corner's absolute degrees. `RemoteSource.tile_for(lat, lon)` renders it;
  `fetch/run.py`'s `_tile_name` goes. P5: GLO-30's template gives
  `_tile_name`'s string on six corners, the poles and the antimeridian among
  them.
- `tile_spacing`: the largest `(lon, lat)` cell size of any tile, in degrees.
  `_tiles` plans its box with these as the nominal spacing; `TILE_SPACING`
  goes. GLO-30's is `(0.0027778, 0.00027778)`: 10" and 1", each rounded up in
  the fifth significant digit so that the box grown from it is never smaller
  than the one a tile's own spacing gives (P6: no band exceeds it, where 0.001
  is exceeded in the 85°-90° band). The box shrinks in latitude (6 cells of
  1", 0.00167°, against 0.006°) and grows slightly in longitude (0.00667°
  against 0.006°): on P6's 60° N case the tiles chosen go from 9 to 3.

Tiles stay 1° squares named by the south-west corner (section 3 hands back
other layouts).

### 4.5 `--land-cover` by path (B5)

`--land-cover VALUE` follows `--dem`'s rule
(`src_python/tin_engine/cli.py@b1353075:1888-1893`): a value that is a key of
the merged catalogue is that cached source, as today; anything else is a path,
so `./mapbiomas-c11` is a file. For a path:

- The file must exist and its header must read as a single-band `uint8`
  GeoTIFF with a CRS; it is read at the option check, before any meshing, and
  each refusal names the file. A missing file: "neither a catalogue source
  (anadem-v1, glo30, mapbiomas-c11, ...) nor a file". `--land-cover-year` with
  a file is refused ("a file has no years"). No cache is needed.
- The repository is `TiffDemRepository([path], z_unit="metre")` (4.7), which gains
  `load_class_window(name, window)`: its `load_window` with
  `decode_class_window` in place of `decode_window`. The tile name is the file
  name.
- The record: `land_cover_source` is the file name, `land_cover_codes`
  "`<file name>` class number; 0 = no data, and every constraint line"; no
  `land_cover_credit`, `land_cover_licence_note` or `land_cover_cite` (a file
  carries none; section 9, Q6).
- New option `--land-cover-map NAME|FILE.toml`: the class map that names and
  colours the classes in the `.gpkg` (41's loader). Default: the catalogue
  entry's `class_map`, none for a file (the `.gpkg` then carries codes
  without names, as a catalogue source without a class map does today). A map
  with no `code_system` (`property`, any uncoded map) is refused. Without
  `--land-cover`, refused ("needs --land-cover").
- Unchanged: "needs --dem", the one-`land_cover_code`-system refusal beside a
  coded `--features-map`, and a mesh reaching past the raster refused.

`land_cover_fractions.assign` takes `repository: ClassRaster`, a protocol of
the two methods it calls (`footprints`, `load_class_window`), in place of
`CacheRepository`; both repositories satisfy it.

### 4.6 No-data in a class raster (B7)

`decode_class_window` consults `meta.nodata`, which is the catalogue entry's
`nodata` for a cached source and the file's `GDAL_NODATA` tag for a path
(`read_page` already refuses a tag a `uint8` cannot hold: P2, `-9999`):

- No-data absent or 0: as today, byte for byte. Empty blocks are 0.
- No-data `n` from 1 to 255: the window starts filled with `n`, the stored
  blocks are decoded into it, then every cell equal to `n` becomes 0. Class 0
  stays rasputin's "no data" in every file it writes, so the record's
  `land_cover_codes` sentence stays true.
- No-data `n` not 0 and a stored cell holding 0: refused, a `GeoTiffError`
  naming the file, the value 0 and `n` ("class 0 means no data in rasputin's
  output; this raster's no-data is 255 and it holds class 0"). Merging a real
  class 0 into no data silently would be a wrong answer without a word
  (section 9, Q5). The check sees only the window read for the chunk at
  hand, not the whole file: a raster whose stored 0 lies outside every
  window the mesh needs is read without a word, which is right, since that 0
  reaches no triangle.
- A raster with no tag whose scheme uses a class value for no data (GlobCover's
  230) keeps that value as a class; the class map names it.

### 4.7 Cell areas from the raster's CRS (B4)

`row_areas(window, crs)` takes the raster's CRS (`meta.crs`; for a raster in
feet, its metre twin, whose `delta_x` and `delta_y` are already in metres):

- Projected (`parse_crs(crs).is_geographic` false): every row
  `delta_x * delta_y` m². These are areas on the projection's grid plane,
  not on the ellipsoid (they differ by the projection's area scale, which
  for a conformal projection is the square of the point scale). The scale
  cancels in a triangle's fractions, up to its change across the triangle's
  few cells, and rasputin writes fractions, not the m². P3: on the
  EPSG:3035 fixture the one triangle gets class 10 at 0.75 and 20 at 0.25,
  areas 240,000 and 80,000 m².
- A class raster reads no heights. Increment 43's height-unit check refuses a
  file whose coordinates are in feet and which states no height unit (P9:
  "the DEM's coordinates are in US survey feet (EPSG:2229) and the file does
  not say what unit its heights are in; give --z-unit us-ft, ft or m"), and
  `--z-unit` does not apply to `--land-cover`. So both land-cover routes
  build their repository with `z_unit="metre"` (4.1), and the fetcher reads a
  classes entry with `read_unit` (4.2): a stated unit for values that are not
  heights, which changes nothing else the header gives (P9: crs EPSG:26945,
  cell 30.480060960121918 m, the same as a DEM in that CRS). For a file in
  degrees or metres that states no height unit, giving `"metre"` returns the
  same meta as giving nothing (P9b: the cached MapBiomas header, which states
  no vertical GeoKey), so MapBiomas keeps its bytes. Left as it is: a class
  raster whose GeoKeys state a height unit other than metres
  (VerticalUnitsGeoKey 4099 = 9003) is refused, naming the file ("the file
  says its heights are in US survey feet; --z-unit says metres", P9). Letting
  it through would need a "no heights" switch threaded through the reader's
  three entry points and both repositories, for a key a class raster has no
  reason to carry.
- Geographic: today's per-row polygon, with a `Geod` built from the CRS's
  ellipsoid as `Geod(a=semi_major_metre, rf=inverse_flattening)`, or `f=0`
  for a sphere (`inverse_flattening` 0). Not `CRS.get_geod()`: for EPSG:4326
  its flattening differs from `Geod(ellps="WGS84")`'s in the last bits, and
  8,862 of 20,929 MapBiomas rows sampled get different areas (P7); built from
  `a` and `rf` it is bit-equal, 0 rows differ, so every MapBiomas mesh keeps
  its bytes. For the same reason `geodesy.cell_areas_by_row` is not reused
  (P4b: 2,301 of 12,000 MapBiomas rows differ in the last bits; its rows are
  centred on nodes, not cut from the cell corner).
- `_GEOD` goes.

### 4.8 Shared TOML reading

41's reader (`src_python/tin_engine/class_maps.py@b1353075:191-200`) moves into
a new module, `tin_engine/toml_data.py`, `parse(text, model, label) -> model`:
`tomllib.loads`, `model.model_validate`, and the same error text
(`"<label>: <key path>: <message>"`, the TOML error with its line). Class maps
call it with `label = f"class map file {shown}"`, so every 41 refusal keeps
its words; the catalogue with `label = <file name>`. Imports: `tomllib`,
`pydantic` and nothing of `tin_engine`.

### 4.9 What does not change

The C++ core; `_core.class_areas`; the fetcher's request, cache and manifest
format; the class-map file format; the tile cache layout; every DEM mesh;
every MapBiomas mesh from the cache (4.6 and 4.7 keep the bytes; T12).
`rasputin fetch glo30` plans a different box (4.4): this is intended, and
the region hash in the manifest, which is of the region as given, does not
change.

## 5. Interfaces and types

```python
# tin_engine/toml_data.py (new; layer L0, imports nothing of tin_engine)
def parse[M: BaseModel](text: str, model: type[M], label: str) -> M: ...

# tin_engine/sources.py
class RemoteSource(BaseModel):            # extra="forbid", frozen
    ...                                   # every field as today
    nodata: float | None = None
    tile_name: str | None = None          # cog-tiles only
    tile_spacing: tuple[PositiveFloat, PositiveFloat] | None = None  # cog-tiles only, (lon, lat) degrees
    def tile_for(self, lat: int, lon: int) -> str: ...
    @property
    def read_unit(self) -> Literal["metre", "foot", "US survey foot"] | None: ...
                                          # "metre" for classes, else z_unit; the Literal
                                          # as written today, since io.models imports NumPy (T4)
class StationSource(BaseModel): ...       # extra="forbid", frozen; fields as today
@dataclass(frozen=True, slots=True)
class Catalogue:
    sources: Mapping[str, RemoteSource]   # MappingProxyType
    stations: Mapping[str, StationSource]
SOURCES: Mapping[str, RemoteSource]       # shipped, read at import
STATION_SOURCES: Mapping[str, StationSource]
def catalogue(user_file: Path | None) -> Catalogue: ...   # ValueError naming the file

# tin_engine/cli.py
def data_root(environ: Mapping[str, str]) -> Path | None: ...   # the one reader of RASPUTIN_DATA
def cache_root(option: Path | None, environ: Mapping[str, str]) -> Path | None: ...  # via data_root
# mesh: --land-cover-map NAME|FILE.toml (str | None)

# tin_engine/land_cover_fractions.py
class ClassRaster(Protocol):
    def footprints(self) -> tuple[TileFootprint, ...]: ...
    def load_class_window(self, name: str, window: IndexWindow) -> ClassWindow: ...
def row_areas(window: ClassWindow, crs: str) -> npt.NDArray[np.float64]: ...   # no default
def assign(vertices, triangles, *, mesh_crs: str, repository: ClassRaster, name: str,
           chunk: int = CHUNK, threads: int = 1) -> Assigned: ...

# tin_engine/io/repository.py
class TiffDemRepository:
    def load_class_window(self, name: str, window: IndexWindow) -> ClassWindow: ...

# tin_engine/io/cog.py: decode_class_window's signature unchanged; behaviour per 4.6
```

Layering rows (`tests/python/test_layering.py`, `@tester`'s red commit):
`toml_data` in L0 with no imports; `sources: "toml_data"`;
`class_maps` gains `toml_data`; `land_cover_fractions` keeps its row (it
already imports `crs`). `sources.py` stays free of `tin_engine.fetch`, NumPy,
pyproj and tifffile, so the mesh path still names a key without them.

## 6. Constants

- `tile_spacing = (0.0027778, 0.00027778)` for GLO-30: degrees, the
  handbook's 10" and 1" rounded up; assumes GLO-30's published bands (any
  latitude, up to the poles); checked at 87° S, 82° N and 60° N (P6).
- No other new constant. The no-data mapping uses the file's value; class 0
  as "no data" is 26b's convention, kept.

## 7. Red suite for `@tester`

Fixtures and figures are the appendix's probes on `b1353075`. Expected values
that come from master (the four entries) are snapshotted in the red commit,
before the Python entries go.

**T1 shipped catalogue unchanged.** `SOURCES` has exactly `anadem-v1`,
`glo30`, `mapbiomas-c11`; `STATION_SOURCES` exactly `nve-hrd`. Each entry's
`model_dump(mode="json")` equals master's, saved under `tests/python/data/` in
the red commit from `b1353075`, except that `glo30` adds `tile_name` (the
template of 4.4) and `tile_spacing` `[0.0027778, 0.00027778]`. `notice()`
for each of the four gives master's text byte for byte.
**T2 file refusals**, each a `ValueError` naming the file and the key path:
an unknown key; an `id` key inside a table; `crs` missing; `kind` not one of
the three; `one-cog-per-year` without `years` or without `{year}`;
`cog-tiles` without `tile_list_url`, `tile_name`, `tile_spacing` or `{name}`;
`tile_name` or `tile_spacing` on a `one-cog` entry; a `tile_spacing` of 0 or
negative; `credit` with `{when}`, with `{date.year}`, with `{url[0]}`;
`tile_name` with `{x}`; `z_unit` on a `values = "classes"` entry; a TOML
syntax error (the message carries the line). `read_unit` is `"metre"` for
`mapbiomas-c11` and `z_unit` (`"metre"`) for the two others.
**T3 user entries.** With `RASPUTIN_DATA` unset, or set to a directory with
no `sources.toml`, the catalogue is the shipped one. A file with
`[source.my-dem]` (`one-cog`, the URL of `fetch_fixtures`' served projected
file, `crs = "EPSG:25833"`, a credit and licence) makes `rasputin fetch
my-dem` copy it as `entry()` does in `test_fetch_run.py`, and `mesh --dem
my-dem` read it from the cache. `[source.glo30]` in the user file is refused
naming both files; a `[station.x]` table is refused; a `class_map =
"maps/x.toml"` resolves against the user file's directory. The existing
`monkeypatch` suites (4.3's list) pass unchanged.
**T4 imports.** `test_sources_imports_pydantic_only`
(`tests/python/test_cli_mesh_cache.py@b1353075:420-428`) becomes: `sources`
imports no `tin_engine.fetch`, NumPy, pyproj or tifffile; standard library,
Pydantic and `tin_engine.toml_data` allowed.
**T5 41's refusals.** 41's class-map refusal suite passes with no edit (the
shared reader keeps every message).
**T6 tile names.** `SOURCES["glo30"].tile_for(lat, lon)` for `(60, 10)`,
`(-15, -45)`, `(0, 0)`, `(-1, -1)`, `(89, 179)`, `(-90, -180)` is
`Copernicus_DSM_COG_10_N60_00_E010_00_DEM`, `..._S15_00_W045_00_DEM`,
`..._N00_00_E000_00_DEM`, `..._S01_00_W001_00_DEM`, `..._N89_00_E179_00_DEM`,
`..._S90_00_W180_00_DEM` (P5).
**T7 tile planning.** With GLO-30's `tile_spacing`, the boxes of P6 choose:
87° S, `x_max` 165.9937 (box 165.2..165.9937 × -87.8..-87.2): 2 tiles, where
0.001 chose 1; 60° N, every edge 0.003° inside the 10°-11° E, 60°-61° N tile:
3 tiles, where 0.001 chose 9. `TestF8Glo30`'s entry
(`tests/python/test_fetch_run.py@b1353075:443-451`) gains the two fields in
the red commit; its expectations (2 listed, 2 unlisted, no request for N60)
stay.
**T8 projected class raster.** P3's fixture: 8 × 8 cells of 100 m in
EPSG:3035, upper-left corner (4,000,000, 3,000,000), `PixelIsArea`, columns
0-3 class 10 and 4-7 class 20. `row_areas` gives 10,000 for all 8 rows;
`assign` on the triangle (4,000,000, 3,000,000), (4,000,000, 2,999,200),
(4,000,800, 3,000,000) in EPSG:3035, repository `TiffDemRepository([file])`,
gives classes `[10, 20]`, fractions `[0.75, 0.25]` exactly, code 10. Today:
NaN rows and an `IndexError`.
*In feet* (P9's fixture, 656 bytes): 8 × 8 `uint8` cells of 100 US survey
feet in EPSG:2229, `PixelIsArea`, tie (0, 0) at (6,500,000, 2,100,000) ftUS,
columns 0-3 class 10 and 4-7 class 20, no vertical GeoKey. Today
`TiffDemRepository([file]).footprints()` refuses it with the "give --z-unit"
message; with `z_unit="metre"` its meta has `crs` EPSG:26945, `file_crs`
EPSG:2229 and `delta_x` 30.480060960121918. `row_areas(window, meta.crs)`
gives 929.0341161327483 m² for all 8 rows; `assign` on the triangle whose
corners are the index-space points (col, row) (2, 0), (6, 0), (2, 4) in EPSG:26945 metres
(x from `meta.x_min - delta_x / 2`, y from `meta.y_max + delta_y / 2`) gives
classes `[10, 20]`, float32 fractions `[0.75, 0.25]` (the float64 ratio is
0.7499999999995 after the round trip through the CRS), code 10. The same
file through the CLI, `mesh --land-cover <feet>.tif`, is not refused at the
option check. With VerticalUnitsGeoKey 4099 = 9003 added, the land-cover
route refuses it with a `GeoTiffError` naming the file (4.7). At the
fetcher, `parse_prefix(<feet bytes>, nodata=None, z_unit=entry.read_unit)`
for a `values = "classes"` entry reads it.
**T9 ellipsoid.** On `class_fixtures`' EPSG:4326 raster, `row_areas` equals
`Geod(ellps="WGS84")` row by row, bit for bit, for all 40 rows; on a window
with the MapBiomas header's spacing (0.00026949458523585647°, top
5.422904791408522°), every seventh row of 146,501, bit for bit (P7: 0 of
20,929 differ). On the same 1/128° grid moved to 60° N in EPSG:4230, the top
row's cell is 379,523.462526 m² (International 1924), not 379,488.217823
(WGS 84), relative tolerance 1e-9 (P4).
**T10 no-data.** P2's fixture: 32 × 32, tiles 16 × 16, classes 7 (columns
0-15) and 9 (16-31), the top-left 4 × 4 cells 255, `GDAL_NODATA` "255", tile 3
empty. The whole window holds classes {0, 7, 9}, 272 zeros (16 from 255, 256
from the empty tile), no 255; today {0, 7, 9, 255} with 256 zeros. With cell
(20, 2) set to 0 as well (P2b): a `GeoTiffError` naming the file, 0 and 255.
`class_fixtures.class_tiff(empty=(0,))` (no tag): equal to its oracle, 356
zeros, as today (P1).
**T11 the CLI by path.** The scene of `test_cli_mesh_land_cover_fractions.py`
(its DEM, `class_fixtures`' raster): `--land-cover <file>.tif` holding the
2025 raster's bytes writes `land_cover_offsets`, `land_cover_classes`,
`land_cover_fractions` and `land_cover_code` equal, array for array, to
`--land-cover test-classes` from the cache, and passes that suite's oracle;
the record's `land_cover_source` is the file name and has no
`land_cover_credit`. Refused, naming the option, before any phase is timed:
`--land-cover-year 2025` with a file; a missing path (the message lists the
catalogue keys); a float32 GeoTIFF (classes must be `uint8`);
`--land-cover-map property`; `--land-cover-map` without `--land-cover`. With
`--land-cover-map mapbiomas-c11` the `.gpkg`'s triangles carry MapBiomas
names and colours, as the cached route's do.
**T12 byte identity.** The existing land-cover suites
(`test_land_cover_overlap.py`, `test_io_class_window.py`,
`test_cli_mesh_land_cover_fractions.py`, `test_cli_mapbiomas.py`) pass with
no edit except T7's entry, T4's import test, and the four `row_areas(window)`
calls of `test_land_cover_overlap.py` (on its lines 230, 242 and 247, two on
the last), which gain `"EPSG:4326"` (no default: a default of WGS 84 is what
hid B4). The signature quoted in `tests/python/fraction_fixtures.py:28-30`
follows: `row_areas(window: ClassWindow, crs: str)`, the area on the
raster's own ellipsoid, or on its grid plane for a projected raster.

The MCP schema test that quotes `mesh`'s help (`tests/python/test_mcp_schema.py`)
follows the new `--land-cover` help and the new option. No mutation round: no
topology or invariant-critical code; the risks are a typo in moved data (T1
against master's snapshot) and the no-data branch (T10).

## 8. Size, acceptance, order

| Where | Added | Removed |
|---|---|---|
| `sources.py`: loader, `catalogue`, two fields, format checks, `tile_for`, `read_unit` and its refusal | 60 | 80 (the entries) |
| `toml_data.py` new; `class_maps.py` calls it | 15 | 12 |
| `fetch/run.py`: `tile_for`, `tile_spacing`, `read_unit` | 5 | 10 |
| `cli.py`: `data_root`, the catalogue threaded through every site that reads `SOURCES`, `--land-cover` by path, `--land-cover-map`, the record, `z_unit="metre"` on both land-cover routes | 86 | 31 |
| `land_cover_fractions.py`: `ClassRaster`, `row_areas` | 15 | 6 |
| `io/repository.py`, `io/cog.py` | 15 | 2 |
| **Total** | **about 195** | **about 145** |

About **50 net lines**; the data moves to a TOML file, which does not count
(`CLAUDE.md` §2). With 41's worst case, +39 % on what is added
(195 × 1.39 = 271), about **125 net**. Design review round 1 added about
6 lines and 2 removed (the height-unit route of 4.7); the estimate does not
move. The audit's 160 counted the data as
staying. One PR, well under 700. Not the short path: over 50 lines in the
worst case.

**Refine or mesh code: none**, so no `@perf` acceptance run.
`land_cover_fractions.py` runs after the mesh is made.

**Speed (P8, measured on an Apple M1 Max on AC power, Python 3.14.7,
median of 15 fresh interpreters).** Every run of any command that imports
`tin_engine.cli` gains the catalogue load: reading the shipped file through
`importlib.resources` (0.71 ms for a 4.8 kB file of the same package) and
parsing and validating the four entries (0.31 ms for a 3,036-byte draft in
4.2's shape), about **1 ms against 189 ms** for importing `tin_engine.cli`
today. Part of the 0.31 ms is validation the Python entries already pay at
import. `tomllib` adds nothing: `tin_engine.cli` already imports it, through
`class_maps`. A user file adds one more read and parse of its size. A mesh
with `--land-cover` gains one CRS read and one `Geod` per chunk of 2^18
triangles in `row_areas`: 0.70 ms (EPSG:4326) or 1.0 ms (EPSG:3035) on the
first call in a fresh process, 0.032 ms on later ones; today's per-row loop
on a 1,000-row MapBiomas window takes 2.2 ms. Neither is worth caching.

Order: red (`@tester`, T1-T12 and the layering rows), green (`@developer`),
code review. The master snapshot for T1 is taken in the red commit, from
`b1353075` or the branch's merge base if master has moved.

## 9. Defaults taken while Ola was away

Each is a question for Ola, with the default the design takes.

- **Q1** Where do your own sources go? Default: a file `sources.toml` at the
  data root (`$RASPUTIN_DATA/sources.toml`), read whenever `RASPUTIN_DATA`
  is set; no new option.
- **Q2** May your own file redefine a shipped source, say `glo30` pointing at
  a mirror? Default: no, it is refused; give the mirror its own name. This
  keeps a shipped name's credit and licence fixed.
- **Q3** May your own file add gauging-station lists? Default: no; station
  fetching follows NVE's map-service layout, which is code.
- **Q4** How does a land-cover file get class names and colours in the
  GeoPackage? Default: a new option, `--land-cover-map`, taking a shipped
  class-map name or a `.toml` file; without it a file's classes are numbers
  only.
- **Q5** A land-cover raster whose stated no-data is not 0 (say 255) but which
  also holds the value 0: refuse it, or count its 0 as no data? Default:
  refuse, with a message naming both values, since rasputin's output uses 0
  for no data.
- **Q6** A land-cover file carries no credit or licence: should the mesh
  record say so in words? Default: no; the credit fields are simply left out,
  and the source is the file name.

## Appendix: the probes

Scripts and saved output in `docs/increments/44-probes/`, run on `b1353075`
with the worktree's `.venv` (Python 3.14.7, the editable install of this
tree, asserted by each script):

- `probe_catalogue.py <tree>` (P1-P5, P2b), output `probe_catalogue.out`;
- `probe_rows.py <tree> /Users/skavhaug/projects/rasputin_data` (P4b: today's
  `row_areas` against `geodesy.cell_areas_by_row` on the cached MapBiomas
  header; P6a: each GLO-30 band's own growth against the nominal one), output
  `probe_rows.out`;
- `probe_tiles.py <tree>` (P6: tiles chosen on three boxes), output
  `probe_tiles.out`;
- `probe_geod.py` (P7: a `Geod` from the CRS against `ellps="WGS84"`),
  output `probe_geod.out`.
- `probe_speed.py <tree>` (P8: the catalogue load and the CRS read in
  `row_areas`, timed), output `probe_speed.out`; run on `10170928`, which
  holds the same code as `b1353075`;
- `probe_feet.py <tree> <scratch dir> /Users/skavhaug/projects/rasputin_data`
  (P9: a class raster in a feet CRS through increment 43's height-unit check,
  as merged; P9b: the cached MapBiomas header), output `probe_feet.out`; run
  on `10170928`.

The figures quoted above, abridged (`...` and the box growths in brackets are
in the saved files):

```
P2  meta nodata 255.0 tag | cells 255: 16 cells 0 (empty tile today): 256 classes: [0, 7, 9, 255]
P2  nodata -9999 on uint8: GeoTiffError GDAL_NODATA (42113) = '-9999' cannot be held by a uint8 cell
P2b cells 0: 1 at [(20, 2)]
P3  row_areas today on it: [nan, nan, nan]
P3  assign today: IndexError index 2 is out of bounds for axis 0 with size 2
P3  with constant 1e4 m2 rows: [10, 20] [0.75, 0.25] areas m2 [240000.0, 80000.0]
P4  EPSG:4326 WGS 84 ... cell 379488.217823 m2, 10 rows down 380377.918340, ratio 1.002344474676
P4  EPSG:4230 International 1924 ... cell 379523.462526 m2, 10 rows down 380413.219975, ratio 1.002344406966
P4b brazil_coverage-col11_2025: rows checked 12000, rows differing in a bit (worst window) 2301, max rel diff 1.2e-10
P5  template equals _tile_name on 6 corners: True | TILE_SPACING today 0.001
P6a today 0.001, 0.001: grows x 0.006000, y 0.006000; bands whose own growth exceeds it: ['85-90']
P6a design 0.0027778, 0.00027778: grows x 0.006667, y 0.001667; bands whose own growth exceeds it: none
P6  S87, east edge 0.0063 deg west of 166 E (10" band): today 1 tiles; own spacing 2 tiles; design 2 tiles
P6  N60, every edge 0.003 deg inside 10..11 E, 60..61 N: today 9 tiles; own spacing 1 tiles; design 3 tiles
P7  CRS(EPSG:4326).get_geod(): MapBiomas rows sampled 20929, differing 8862
P7  Geod(a=, rf=) from CRS(EPSG:4326).ellipsoid: MapBiomas rows sampled 20929, differing 0
P8  parse and validate the catalogue 0.314 ms | read a 4.8 kB shipped file 0.712 ms | import tin_engine.cli 188.6 ms | tomllib already imported: True
P8  CRS read and Geod: EPSG:4326 first 0.702 ms, EPSG:3035 first 1.005 ms, later 0.032 ms | row_areas, 1,000 rows today 2.150 ms
P9  read_page, no z_unit: GeoTiffError the DEM's coordinates are in US survey feet (EPSG:2229) and ... give --z-unit us-ft, ft or m
P9  z_unit='metre': crs EPSG:26945 file_crs EPSG:2229 delta 30.480060960121918
P9  4099 = 9003, z_unit='metre': GeoTiffError the file says its heights are in US survey feet; --z-unit says metres
P9  rows of delta_x * delta_y m2: 929.0341161327483 | classes [10, 20] fractions [0.7499999999995035, 0.25000000000049644]
P9b vertical GeoKeys: none | meta with z_unit='metre' equals meta with none: True
```

## Review

**Design review, round 1, 2026-10-13.** Range `b1353075..ae71ea85`, 0 production lines, @reviewer: CHANGES REQUESTED. Blocking: (1) a land-cover file in a feet CRS is refused by increment 43's height-unit check (`read_page`: a uint8 class raster in EPSG:2229 fails with "the DEM's coordinates are in US survey feet … give --z-unit"; `--z-unit` does not apply to `--land-cover`), contrary to section 4.7: say how class rasters skip the height-unit check (for example the land-cover path passes `z_unit="metre"`) and add a test to T8; (2) no speed figure for what every mesh run gains (the TOML catalogue parsed at import, the CRS parse in `row_areas`); the reviewer measured about 4 ms to import `tomllib` and 0.3 ms to parse a 3 kB file. Probes re-run and reproduced (NaN and IndexError; `Geod(a, rf)` 0 of 20,929 rows changed against 8,862; the 85°-90° band). Suggestions: T12 has four `row_areas` calls on three lines, not three; `tests/python/fraction_fixtures.py:28` quotes the old one-argument signature; the `class_maps.py:191-199` citation leaves out line 200; 4.6 should say the "stored 0" refusal sees only the window read; 4.3 should say `--cache` alone (without `RASPUTIN_DATA`) loads no user entries; 4.7 should say the m² from projected rows are grid-plane areas (fractions unaffected).
