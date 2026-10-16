# Increment 55: an agent finds the local elevation data and passes it by name

**Status: designed** (branch `worktree-mcp-datasets`), design review round 1 APPROVED WITH WORDINGS, applied by the main session.

## 1. What Ola asked, and why it failed

Ola ran Codex against rasputin's MCP server (the Model Context Protocol
server, `rasputin mcp`). Codex fetched Copernicus GLO-30 heights but could not
pass them to catchment delineation, and never found out that a 10 m terrain
model of all Norway already sits on the machine. Ola asks that an agent with
no access to the source or the docs can discover the local elevation datasets
and their properties, choose one by coverage, resolution and purpose, prefer
suitable local data over downloading worse data, pass dataset names between
tools, and run catchment then mesh; with the metadata exposed for the agent
to reason on, no choice logic for Norway, existing code reused, and
integration tests for discovery, choice and use. Ola, verbatim, is in the
main session's brief; the words that bind this design are "Reuse existing
functionality and avoid unnecessary abstractions or duplicated data
management logic" and "Ensure that the MCP exposes enough metadata for
independent agent reasoning, rather than embedding dataset-selection
decisions specific to Norway."

The causes, checked against the code at `9a7c312d`:

1. **`catchment` and `station-catchments` take paths only.** `mesh --dem KEY`
   turns a catalogue key into a cached source through `cli._dem_sources`;
   `catchment --dem` and `station-catchments --dem` are `list[Path]` and go
   straight to `repository_for(tuple(dem))`
   (`src_python/tin_engine/cli.py`, `def catchment`, `def station_catchments`),
   so `--dem glo30` is read as the file `./glo30` (Codex's call failed with
   "cannot read /Users/skavhaug/glo30", main session's reading of the
   session log). Even with the key resolved, delineation loads whole tiles
   (`catchment._grow`: `assemble(plan, repository.load)`), so a cache that
   holds only the fetched box is refused. Measured on the fixture of T6
   below: a cache lacking the 13 blocks of its northern block row, which the
   catchment never reaches, is refused today with "13 of the 325 blocks
   needed are not in the cache".
2. **No local dataset is known anywhere.** The catalogue
   (`src_python/tin_engine/sources.py`, `data/sources.toml`, and the user's
   `$RASPUTIN_DATA/sources.toml`) holds remote sources only, and no command
   lists sources or what they cover. The DTM10 set is a bare directory of
   254 GeoTIFFs (`ls ../rasputin_data/DTM10_UTM33_20260925/*.tif | wc -l`).
3. **Without `RASPUTIN_DATA` there is no user catalogue at all.** Codex starts
   the server with the client configuration a user writes, `command` and
   `args`, no `env`; the MCP Python SDK's own client passes only `HOME`,
   `LOGNAME`, `PATH`, `SHELL`, `TERM` and `USER` on POSIX
   (`mcp/client/stdio.py`, `DEFAULT_INHERITED_ENV_VARS`). `cli.data_root`
   reads `RASPUTIN_DATA` only, so the user's sources file is never read,
   and the server's instructions ("A catalogue DEM source must be fetched
   first: run fetch over the domain, then mesh") point the agent at a
   download.

## 2. Prior art: legacy and literature

Written before the design.

*Literature and practice.*

- **STAC**, the SpatioTemporal Asset Catalog specification
  (github.com/radiantearth/stac-spec). Its Collection's spatial extent is a
  bounding box whose "coordinate reference system ... is WGS 84
  longitude/latitude", south-west corner first (`collection-spec.md`,
  *Spatial Extent Object*); its common metadata has `gsd`, "the nominal
  Ground Sample Distance for the data, as measured in meters on the ground"
  (`commons/common-metadata.md`); its projection extension v2 names the CRS
  `proj:code` (`"EPSG:32659"`), `proj:epsg` having been removed in v2.0.0
  (github.com/stac-extensions/projection). The listing here reports the same
  three things in the same conventions: a lon/lat box west, south, east,
  north; a cell size in metres; the CRS as `EPSG:n` where there is one.
  Departure: no STAC documents, no items per tile and no STAC library. A
  catalogue of a few sources on one machine does not need the format, and
  the agent reads JSON either way. Nothing of STAC's guarantees is dropped,
  since it gives none beyond its schema.
- **The XDG Base Directory Specification** (specifications.freedesktop.org,
  *basedir-spec*): user configuration goes in `$XDG_CONFIG_HOME`, else
  `$HOME/.config`; a relative value is ignored. `cli.cache_root` already
  follows its cache half; this increment follows its configuration half for
  the user's sources file.
- **MCP** (modelcontextprotocol.io, specification 2025-06-18): a server's
  `instructions` are shown to the model at initialise; a tool result may carry
  `structuredContent`. Both are used as the server already uses them.
- **GDAL's `gdalinfo -json`** reports a `wgs84Extent` per raster: the same idea
  for one file. GDAL is prohibited (`CLAUDE.md` §2); nothing is taken but the
  idea, which STAC states anyway.

No novelty is claimed.

*Legacy.* The legacy tree read a data directory from `RASPUTIN_DATA_DIR`,
falling back to the working directory with a warning, and took fixed
subfolders under it (`dem_archive`, `tin_archive`, `globcov`, `corine`):

```
$ git grep -l -i -E 'RASPUTIN_DATA|dem_path|list_datasets|available.*(dem|dataset)' legacy-archive -- legacy
legacy-archive:legacy/rasputin/__init__.py
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/geometry.py
legacy-archive:legacy/rasputin/web_visualize.py
legacy-archive:legacy/tests/test_gml_repository.py
legacy-archive:legacy/tests/test_land_cover_repository.py
legacy-archive:legacy/tests/test_raster_repository.py
```

Read: `legacy/rasputin/__init__.py` lines 4-6 and `application.py` lines
55-61 and 179-181 hold the variable and the fixed folder names; the other
files use the resulting paths. Nothing is carried across: fixed folder names
are the machine-specific layout this increment avoids, and there was no
listing.

## 3. Ideas file

Taken up: **"A `rasputin sources` command, and tile layouts beyond GLO-30's"**
(the command; other tile layouts left out) and **"A default place for the user
source list; Intel hybrid core split"** (the default place; the core split
left out). Overlapping, left out: "Structure for rasputin_data and
rasputin_scratch" (this increment registers datasets in a file instead of
scanning folders, so no folder structure is needed for discovery), "Intent
layer over rasputin", "Two small fixes from 37's code review" (the heartbeat
await), "Data-dependent tests never run in a worktree".

## 4. The design

### 4.1 Data flow

```
user's sources.toml  ──┐                     (located by sources.user_source_file:
package sources.toml ──┴─> sources.catalogue   RASPUTIN_DATA, else XDG_CONFIG_HOME, else HOME)
                              │ Mapping[id, RemoteSource]  (kind "local" or a remote kind)
          ┌───────────────────┼─────────────────────────────┐
          v                   v                             v
  cli._dem_sources      cli `sources`                 cli `fetch`
  (mesh, catchment,     -> source_report.report        refuses kind "local"
   station-catchments)     reads headers only:
          │                local: TiffDemRepository
          v                remote: CacheRepository
  local:  (path,)  ──> repository_for -> TiffDemRepository  (unchanged)
  remote: CachedSource ─> repository_for -> CacheRepository (unchanged)
          │
          v
  mesh: open_dem (unchanged)   catchment: delineate, now windowed + checked
```

No second registry: a local dataset is one more entry in the existing
catalogue, read by the existing readers. The only new module is the report,
because the listing needs pyproj and shapely and `sources.py` must not import
them (its docstring: the mesh path names a key without them).

### 4.2 A local entry in the catalogue (`sources.py`)

`RemoteSource` gains `kind = "local"`, a `path` and a `description`. The class
keeps its name: renaming it touches every user of it for no behaviour, and
increment 54 rewrites the same lines of `cli.py` (section 8). Its docstring
says it holds local entries too.

```toml
# In the user's sources.toml
[source.dtm10-norway]
kind = "local"
path = "/data/DTM10_UTM33_20260925"   # a directory of GeoTIFFs, or one GeoTIFF
z_unit = "metre"                       # for files that do not state it
credit = "..."
licence_note = "..."
description = "Bare-ground terrain model, 10 m, ..."
```

Fields: `path: Path | None = None`; `description: str = ""` (every kind; one
or two plain sentences: what the cells represent, e.g. bare ground or the
surface with vegetation and buildings, and the nominal resolution and
region, as the distributor states them); `crs` becomes `str | None = None`.
Rules in `_consistent`, each a `ValueError` naming the field:

- `kind = "local"`: `path` needed; `url`, `url_template`, `tile_list_url`,
  `tile_name`, `tile_spacing`, `years` and `crs` refused (the CRS is the files'
  own, read from their headers; one source of truth); `values = "classes"`
  refused ("a local class raster is given to --land-cover as a file").
- any other kind: `path` refused, `crs` needed (as today, so `fetch/run.py`
  narrows it at each of its three reads (`_header`, `_tiles`, `_region`)).
- `catalogue()` resolves a relative `path` from the sources file's directory,
  as it already does for `class_map`. It does not check that the path
  exists: the listing reports a missing path as that source's problem, and a
  run refuses it as any unreadable `--dem`.

The shipped `data/sources.toml` gains a `description` for each source, from
the distributors' own pages, cited in a TOML comment: GLO-30 is a surface
model (DSM) including buildings and vegetation, 1 arc-second, global
(Copernicus DEM Product Handbook); ANADEM a bare-ground model for South
America derived from GLO-30 (Laipelt et al. 2024); MapBiomas Collection 11 a
30 m land-cover map of Brazil. `@developer` quotes each from its source; the
review checks the quotations.

### 4.3 Where the user's sources file is found, with no environment variable

**Decision: yes, the server finds it with no variable set**, in the XDG
configuration directory, and the order is:

1. `$RASPUTIN_DATA/sources.toml`, when `RASPUTIN_DATA` is set and the file exists;
2. else `$XDG_CONFIG_HOME/rasputin/sources.toml` when that variable is an
   absolute path, otherwise `$HOME/.config/rasputin/sources.toml` when `HOME`
   is absolute, if the file exists;
3. else none.

The first file that exists is the one read; one file only, never merged,
so an id is defined in one place. Reasons:

- It is the platform's own rule for user configuration (XDG), the same rule
  `cache_root` already applies to the cache, so it adds no convention.
- `HOME` is among the six variables the MCP SDK's client passes by default,
  and is set in every desktop client's child; Codex's server process had it
  (it fetched into `~/.cache/rasputin`). `RASPUTIN_DATA` is not passed unless
  the user writes it.
- It names no machine and no country: the path is the user's own config
  directory on any POSIX system.
- "Exists" rather than "is set" for step 1 means a user who keeps the file in
  `~/.config/rasputin` sees the same catalogue in a shell where
  `RASPUTIN_DATA` is set for the cache and in an agent's server where it is
  not.

Rejected: scanning `RASPUTIN_DATA` or the home directory for GeoTIFF folders
(an unbounded walk, no credit or licence, and a guessed answer to "what is
this data"); a new `data_root` key or pointer file (one more indirection,
and the cache already has a default without one); any default path under
the checkout (`../rasputin_data` is Ola's layout, not a product rule).
The cache keeps its rule: `--cache`, else `$RASPUTIN_DATA/cache`, else the XDG
cache. A cache kept under `RASPUTIN_DATA` is therefore not seen by a server
started without the variable; the listing says which cache it read, and
every command that reads one takes `--cache` (section 4.5).

Interface, in `sources.py` (pure: an environment mapping in, no NumPy, no
pyproj):

```python
def user_source_candidates(environ: Mapping[str, str]) -> tuple[Path, ...]:
    """The places looked at, in order: at most two."""

def user_source_file(environ: Mapping[str, str]) -> Path | None:
    """The first candidate that is a file."""
```

`cli.cache_root` is unchanged. `cli._known_sources` reads `user_source_file(environ)`.
`cli.data_root` stays, now used by `cache_root` and
`user_source_candidates` only.

### 4.4 `--dem KEY` everywhere (`cli.py`)

`_dem_sources(dem, root, known)` returns a third value, the catalogue entry
(None for paths). For a local key it returns `((entry.path,), None, entry)`:
the same directory-or-file path `--dem` already takes, so `repository_for`,
`open_dem` and `TiffDemRepository` are unchanged. For a remote key it
returns `((), CachedSource(...), entry)` as today. A key with other `--dem`
values, and a class source, are refused as today.

- **`mesh`**: the run record's `dem_credit`, `licence_note` and `cite` come
  from the entry whenever there is one (today only for a cached source), so
  a mesh from a local source carries its credit too. The height unit given
  to the reader is `--z-unit`, else the entry's `z_unit`, as today.
- **`catchment` and `station-catchments`**: `--dem` becomes `list[str]`, with
  the help of `mesh --dem` ("a GeoTIFF, several, one directory, or a source
  from `rasputin sources`"); each gains `--cache` (the same option object
  `mesh` uses); each calls `_dem_sources` and passes `paths` and `cached` to
  `repository_for`, with the entry's `z_unit` as for `mesh`.
- **`fetch`** refuses a local key: "`<id>` is a local source, read in place;
  it needs no fetch".
- **Catchment reads only what it needs.** `catchment._grow` calls
  `repository.check(plan)` before `assemble`, and passes
  `load_window=repository.load_window`, as `open_dem` does. The window's
  nodes are the same, so a catchment from files is unchanged node for node
  (the catchment suites stay green, which is the check); a cache need hold
  only the blocks the windows meet, and decoding reads only those blocks
  (also on files). A window whose blocks are missing raises `NotCached`
  carrying the window's node box in the DEM's CRS: `NotCached` gains an
  optional keyword `box: tuple[float, float, float, float] | None`, and its
  message then ends "for the window x0 y0 x1 y1 in the source's CRS".
  `catchment` turns it into a usage error ending
  "run: rasputin fetch <key> --bbox x0 y0 x1 y1" (the fetch box is in the
  source's CRS when `--out-crs` is not given, which is the DEM's CRS here),
  so an agent can fetch the missing area and rerun; `station-catchments`
  stops with the message, as it stops on any `ValueError` today.
- The repository protocol is unchanged: `load_window` and `check` are
  already in `DemRepository`.

### 4.5 The listing: `rasputin sources`

```
rasputin sources [--bbox WEST SOUTH EAST NORTH] [--cache PATH]
```

Prints one JSON document on stdout on one line, so the tool result's `stdout` holds it as one string. It is a tool through the
existing derivation, named `sources`, read-only in its MCP annotations.
`--bbox` is in degrees, WGS 84 longitude and latitude; a point is a box with
west equal to east and south equal to north. Exit 0 whenever the document is
written, even when a source has a problem or the user's file is refused;
the problem is in the document.

New module `src_python/tin_engine/source_report.py` (layer: Python I/O through
`io/repository.py`, pyproj through `crs.py`, shapely; no `_core`):

```python
class SourceReport(BaseModel):            # frozen
    id: str
    kind: str                             # "local", "one-cog", "cog-tiles", "one-cog-per-year"
    values: Literal["heights", "classes"]
    path: str | None                      # the files (local) or the cache directory (remote)
    present: bool                         # local: the path exists; remote: a cache manifest exists
    tiles: int                            # files read (local) or cached objects (remote)
    crs: str | None                       # crs.crs_label of the headers' CRS
    crs_name: str | None                  # pyproj's name, e.g. "ETRS89 / UTM zone 33N"
    cell_m: tuple[float, float] | None    # smallest and largest cell, longer side, metres
    z_unit: ZUnit | None
    z_unit_from: Literal["file", "catalogue", "assumed"] | None
    nodata: float | None
    coverage: tuple[float, float, float, float] | None   # lon/lat box W S E N
    box_covered: float | None             # with --bbox: share of the box covered, 0..1
    description: str
    credit: str
    licence_note: str
    cite: tuple[str, ...]
    problem: str | None                   # plain words when the headers could not be read

class CatalogueReport(BaseModel):
    user_file: str | None                 # the sources file read
    looked_at: tuple[str, ...]            # every candidate path, in order
    register_at: str | None               # where to add a local source: the file read, else the first candidate
    user_file_problem: str | None         # a refused file: its error; the package's sources are still listed
    cache: str | None
    note: str | None                      # see below
    sources: tuple[SourceReport, ...]     # package then user, in file order

def report(entries: Mapping[str, RemoteSource], cache: Path | None,
           box: Bounds | None) -> tuple[SourceReport, ...]:
```

What each figure is, all from the headers through the existing readers:

- **Local**: `repository_for((entry.path,), z_unit=entry.read_unit)` then
  `footprints()`; **remote**: `CacheRepository(cache, id,
  z_unit=entry.read_unit).footprints()` when `<cache>/<id>/manifest.json`
  exists, else nothing read (`present` false, figures None; the description
  carries what the source covers). For a remote source the figures describe
  the cached objects, whose blocks may be only the fetched part; the
  document's `note` says so in one sentence.
- `cell_m`: per tile, `LocalFrame.at(crs, tile centre).cell_m(meta)` (the
  longer side in metres; identity on a projected CRS); min and max over the
  tiles, rounded to 0.1 m.
- `coverage`: each tile's node box, densified to 8 segments per edge
  (`shapely.segmentize` at an eighth of the longer edge), reprojected to
  EPSG:4326 with `crs.reprojector`, unioned (`shapely.union_all`); the
  union's bounds, rounded to 1e-6 degrees. Constant: 8 segments, a chosen
  value at which edge curvature is negligible against a cell; checked at 254
  DTM10 tiles of 50 km (0.29 s to read the headers, 0.02 s for the union,
  measured with the main checkout at `9a7c312d`). A source crossing the
  antimeridian reports a box from -180 to 180: a stated limitation, not
  handled.
- `box_covered`: `area(union ∩ box) / area(box)` in square degrees (a chosen
  measure, not an area on the ellipsoid: it ranks sources for one box, which
  is all it is for), rounded to 1e-4; for a point or line box, 1.0 when it is
  covered by the union, else 0.0. None without `--bbox` or when nothing was read.
- `z_unit_from`: "file" when the header states it, "catalogue" when the entry's
  `z_unit` was used, "assumed" when neither (metres); None for classes.
- `problem`: any `OSError` or `ValueError` from reading that source, in its
  own words (`cli._words` where it applies); the other sources are still read.
- `note`: when no user file exists, "No user sources file: looked at
  <paths>. To use a local DEM, add a [source.<id>] table with kind = "local",
  path, credit, licence_note and description to <register_at>."; when a
  remote source is cached, the partial-blocks sentence above. Plain words,
  no labels.

The CLI command is thin: it builds the catalogue (catching the user file's
`ValueError` into `user_file_problem` and listing the package's sources),
calls `report`, fills the places, and prints `model_dump_json()`.

### 4.6 The MCP side

- `mcp/schema.py`: `READ_ONLY` gains `"sources"`; `RunResult` and the output schemas are unchanged.
- `INSTRUCTIONS` in `mcp/server.py`: the line "A catalogue DEM source must be
  fetched first: run fetch over the domain, then mesh." is deleted, and
  these lines lead:

  ```
  Call sources first: it lists the elevation and land-cover sources this machine knows, local ones read in place and remote ones read from a cache, each with its coverage as a longitude/latitude box, cell size in metres, CRS, height unit, tile count, description and credit. Give it bbox (west south east north, degrees) to see what share of your area each source covers.
  Choose a source by its coverage of your area, its cell size and what its description says the heights represent.
  Pass the chosen source's id as dem to mesh, catchment and station_catchments. A local source needs no fetch. A remote source is read from the cache: fetch the area first; a run that lacks data says what to fetch.
  sources also says which sources file and cache it read, and where a local dataset would be registered.
  ```

  General guidance only: no source, country or resolution is named.

### 4.7 Docs

- `README.md`, *Connect an agent*: the Codex block shows an optional
  `env = { RASPUTIN_DATA = "/path/to/data" }` line, marked as optional (for a
  cache kept there); a short paragraph says local datasets are registered in
  `~/.config/rasputin/sources.toml` (or `$RASPUTIN_DATA/sources.toml`) with
  the TOML example of section 4.2, and that `rasputin sources` lists what is
  registered. No test reads the README's `env` line.
- `project_structure.md`: an entry for `source_report.py`; the `sources.py`
  entry says it holds local entries and finds the user's file.
- `data/sources.toml`'s header comment: where the user's file is looked for.
- Help strings that say "$RASPUTIN_DATA/sources.toml" (`mesh`, `fetch`) say
  "your sources file (see rasputin sources)".

## 5. Deleted

- The `INSTRUCTIONS` line that sends every agent to `fetch`.
- `catchment`'s and `station-catchments`' private path-only handling of
  `--dem` (`repository_for(tuple(dem))` with no key), replaced by the shared
  `_dem_sources`.
- `catchment._grow`'s whole-tile `assemble(plan, repository.load)`.
- The `cached is not None` condition on the record's credit lines, replaced
  by the entry.

## 6. Left out

- Local land-cover entries (`values = "classes"` with `kind = "local"`):
  refused with a pointer to `--land-cover PATH`.
- A command that writes the user's sources file (`sources add`). Registering
  is a one-time act by the person who owns the data.
- Listing the remote tiles that exist but are not cached (needs the network;
  the description says what a remote source covers).
- Coverage as a polygon, and per-tile listings. The box and `box_covered`
  answer "does it cover my area"; a footprint polygon can follow if an agent
  needs it.
- Merging two user files, and a configurable cache location in the file.
- Antimeridian-crossing coverage.

## 7. Tests

Every test is `@tester`'s to write; the fixtures below are measured with the
main checkout at `9a7c312d` (same tree as this branch), by scratch probes run
through the CLI and the library.

**Fixture L, a local dataset**: `catchment_fixtures.bowl()` (400 rows x 200
columns, 100 m, EPSG:25833, corner node (500000, 6600000)) cut by
`mosaic_fixtures.quadrants(row_cut=200, col_cut=100, overlap=1)` and written
with `test_cli_mesh_mosaic.write_tiles`: four files, `ne.tif` 201 x 100,
`nw.tif` 201 x 101, `se.tif` 200 x 100, `sw.tif` 200 x 101 nodes, 80-82 kB
each. Measured: coverage lon 15.0 to 15.351879, lat 59.179575 to 59.538349;
`box_covered` 1.0 for (15.15, 59.33, 15.20, 59.38), 0.3748 for
(15.2, 59.3, 15.6, 59.4), 0.0 for (16, 60, 16.1, 60.1); heights not stated by
the files (the reader says "the file does not say"). `catchment` from the seed
`lat(100, 200)` = (510000, 6580000) in EPSG:25833, no lake: 623 nodes,
6.225000 km2 outline; with `lake_box()`: 766 nodes, 7.655000 km2, three
windows; `mesh --domain` of that outline, `--tolerance 5`: exit 0, 1064
triangles, largest difference 4.527 m. Each run under 0.1 s.

**Fixture C, a cached remote-style source**: the same `bowl()` as one tiled
GeoTIFF (16 x 16 blocks, Deflate, 25 x 13 = 325 blocks, 294,648 bytes,
`geotiff_fixtures.micro_tiff` with scale 100 m) written by
`cog_fixtures.write_cache` under a test catalogue key with `crs =
"EPSG:25833"`, registered in the test's user sources file. Measured with
`delineate` on `CacheRepository`: full cache, 623 nodes, 6,225,000 m2;
with the first block row (13 blocks, rows 0-15) skipped, refused today
with "13 of the 325 blocks needed are not in the cache".

Isolation: `conftest.py`'s autouse fixture gives every test its own empty
`XDG_CONFIG_HOME`, as it already does for `XDG_CACHE_HOME`, so no test reads
the real `~/.config/rasputin`.

- **T1 catalogue rules** (`sources.py`): a local entry with `path` loads, a
  relative path resolved from the file's directory; each refused field of
  section 4.2 is a `ValueError` naming it; a remote entry without `crs`, or
  with `path`, is refused; `description` defaults to "".
- **T2 places** (pure, environment mappings): the candidate order of section
  4.3, the existence rule ("set" but missing falls through), relative
  `XDG_CONFIG_HOME` and `HOME` ignored, `cache_root` unchanged for the existing table of cases.
- **T3 discovery through MCP** (in process, `mcp_driver.call` on
  `build_server(cli.app)` with the production child): `HOME` set to a temp
  directory holding `.config/rasputin/sources.toml` that registers fixture L,
  `RASPUTIN_DATA` and `XDG_CONFIG_HOME` deleted. `sources` returns
  `json.loads(structured_content["stdout"][0])` with: the local entry, `present` true, `tiles` 4, `crs` "EPSG:25833", `crs_name` naming UTM zone
  33N, `cell_m` (100.0, 100.0), `z_unit` "metre" with `z_unit_from`
  "catalogue", the coverage box of fixture L to 1e-6, its description and
  credit; `user_file` under the temp home; the shipped sources listed with
  `present` false and their descriptions; the tool's annotations read-only.
- **T4 choice from the listing**: with `bbox` (15.15 59.33 15.20 59.38), the
  local entry has `box_covered` 1.0 and every uncached remote source None;
  with (15.2 59.3 15.6 59.4), 0.3748; with (16 60 16.1 60.1), 0.0. A second
  local entry, fixture L written as one 200 m file (every second node), lists
  `cell_m` (200.0, 200.0), so the two are told apart by cell size alone.
- **T5 selection then use**: `sources` with bbox (15.15 59.33 15.20 59.38) over T4's two local entries; the test takes the id with `box_covered` 1.0 and the smallest `cell_m`, asserts it is the 100 m entry, then calls `catchment` with `dem=[that id]`, the seed
  and the lake of fixture L: exit 0, the GeoJSON's `nodes` 766; then `mesh`
  with `dem=["<key>"]`, `domain` that file and `tolerance` 5: exit 0, and the
  run record's `dem_credit` is the entry's credit. `fetch` with the local key
  is refused with "needs no fetch", before any network use.
- **T6 cached key in catchment** (no network; CLI is enough): fixture C,
  `catchment --dem <key> --cache <root>` from `lat(100, 200)`: 623 nodes,
  also with the first block row skipped (the red test of the windowed read);
  with the blocks of block row 12 (the seed's) skipped: exit 2, the message
  names the missing blocks and ends with
  `run: rasputin fetch <key> --bbox` and four numbers that contain the seed.
  `station-catchments --dem <key>` resolves the key (a run on fixture C with
  one gauge, or a refusal that is not "cannot read").
- **T7 the listing's own refusals**: a user file with an unknown field gives
  exit 0, `user_file_problem` naming the file, and the shipped sources; a
  local entry whose path is missing gives that entry `present` false and a
  `problem`, the others unaffected; no user file gives `note` naming
  `register_at` and `looked_at`.
- **T8 instructions**: the initialise result's instructions contain
  "Call sources first" and "A local source needs no fetch", and no longer
  contain "must be fetched first".
- **T9 no environment but `HOME`, over stdio**: the production server started
  through `StdioServerParameters` with `env={"HOME": <temp home>, "PATH":
  os.environ["PATH"]}` (what the SDK passes by default, less the user's own
  home), `sources` lists fixture L's entry. This is the Codex condition.
- Updated, not new: `test_mcp_schema.py`'s snapshots (the new tool), `test_catalogue_data.py`'s check of the help text
  naming `$RASPUTIN_DATA/sources.toml`, and `catchment_fixtures.MemoryRepository`,
  which gains `load_window` (a slice, recorded in `loads` as `load` is) and a
  `check` that passes.

No mutation round: no topology or invariant-critical code changes.

## 8. Overlap with increment 54

Increment 54 (open-core cleanup, branch `open-core-cleanup`, designed, not
built) rewrites every command docstring and `help=` string in `cli.py`, moves
command registration to one block at the end of the file, and regenerates the
MCP snapshots. This increment touches `_dem_sources`, the `catchment`,
`station-catchments`, `fetch` and `mesh` signatures and bodies, and adds a
command. Whichever builds second rebases on the other; the conflicts are in
text, not in design. If 54 lands first, `sources` is registered in its block,
under the panel *Data* beside `fetch`. The help texts written here already
follow 54's rule: no history labels.

## 9. Size

Net production lines (`CLAUDE.md` §2; TOML and docs not counted):

| Part | Lines |
| --- | --- |
| `sources.py`: local kind, `path`, `description`, rules, places | 45 |
| `cli.py`: `_dem_sources` entry, mesh credit, catchment and station-catchments keys and `--cache`, NotCached message, fetch refusal, `sources` command | 85 |
| `catchment.py` check and window; `io/cog.py` `NotCached` box; `fetch/run.py` assertion | 15 |
| `source_report.py` | 140 |
| `mcp/schema.py`, `mcp/server.py` | 4 |
| **Total** | **about 300**; worst case 400 |

One PR, the full path (over 50 lines). It touches no refine or mesh code, so
no `@perf` acceptance run.

## 10. Files

`src_python/tin_engine/sources.py`, `source_report.py` (new), `cli.py`,
`catchment.py`, `io/cog.py`, `fetch/run.py`, `mcp/schema.py`, `mcp/server.py`,
`data/sources.toml`; `README.md`, `project_structure.md`, `ROADMAP.md`;
tests under `tests/python/` (T1 to T9, `conftest.py`, the snapshot and helper
updates).

## Review

- **Design review, round 1, 2026-10-16.** Commit `8a29e876`, 0 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: `local` (repeats `kind`), `RunResult.data`/`JSON_STDOUT` and the `Origin`/`Place`/`cache_place` origin tracking dropped; T5 selects from the listing before use; the CRS narrowing at three reads. Held: no second registry; the windowed cache read is needed (`CacheRepository.load` reads whole tiles); the SDK passes HOME by default.
