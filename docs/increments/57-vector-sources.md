# Increment 57: vector datasets in the catalogue, the options each fits, output folders made

**Status: built, code review approved** (branch `worktree-vector-sources`, from master `a0a1ad45`): code green at `08cb70c2` except the output folders, whose design (section 4.6) is revised to Ola's ruling on `--out-parent`: the guard moves into the MCP server.

## 1. What Ola asked, and what failed

Ola's Codex run after increment 55 (`docs/increments/55-dataset-discovery.md`,
merged as #273) found the 10 m terrain model through `sources`, but found the
lake polygons (`rasputin_data/nve_hrd/lakes.geojson`) and the CORINE land
cover (`rasputin_data/corine2018_dtm10_utm33.gpkg`, layer `corine2018`) only
with its own shell (`rg`, `ogrinfo`). With only the MCP tools it could not
have. It also had to `mkdir` an output folder that `catchment` refused as
missing. Ola approved the follow-up ("They all sound good to me"): vector
datasets in the same catalogue, a per-dataset "fits this option" field, tools
that create missing output folders, a Codex test with its shell turned off;
and plugins registering their own datasets ("Yes, do it"). His standing words
from 55 bind this design too: "Ensure that the MCP exposes enough metadata for
independent agent reasoning, rather than embedding dataset-selection decisions
specific to Norway" and "Reuse existing functionality and avoid unnecessary
abstractions or duplicated data management logic."

The causes, checked against the code at `a0a1ad45`:

1. **The catalogue holds rasters only.** `RemoteSource.values` is
   `Literal["heights", "classes"]` and a local entry is "a local GeoTIFF or
   directory of them" (`src_python/tin_engine/sources.py`, `class RemoteSource`);
   `source_report._one` reads every local entry with `repository_for`, a
   GeoTIFF reader.
2. **Nothing tells an agent which option takes a file.** `mesh --features`
   is `list[Path]`, `catchment --lakes` and `--rivers` are `Path`, and
   `station-catchments --stations`, `--rivers`, `--lakes`, `--reference` are
   `Path` (`cli.py`, `def mesh`, `def catchment`, `def station_catchments`).
   A path is enough (MCP passes it as a string); what is missing is the
   listing of the file and of the options it fits. The options stay as they
   are.
3. **A missing output folder is refused.** `cli_kit.destination` ends with
   `if not target.parent.is_dir(): raise typer.BadParameter(f"{target.parent}
   is not an existing directory")`; every writer of `mesh`, `catchment` and
   `station-catchments` goes through it (seven calls in `cli.py`), and
   `palette --out` writes with no check at all (`out.write_text(text)`).

## 2. Prior art: legacy and literature

Written before the design.

*Literature and practice.*

- **STAC's table extension** (github.com/stac-extensions/table): a tabular or
  vector asset is described by `table:columns` ("A list of Column Objects
  describing each column"), `table:primary_geometry` ("The primary geometry
  column name") and `table:row_count` ("The number of rows in the dataset").
  The listing reports the same three: attribute names, the GeoPackage layer
  (whose geometry column the reader finds), and a feature count.
- **STAC's classification extension** (github.com/stac-extensions/classification):
  `classification:classes` lists the coded values present in an asset, each a
  Class Object with `value` and `name`. The listing reports a land-cover
  file's coded attribute and its distinct values; the names stay in the class
  map, as today.
- **STAC's Collection spatial extent** (as cited in 55): coverage as a WGS 84
  longitude/latitude box, west, south, east, north. A vector source's box is
  its layer's extent, reprojected the way 55 reprojects a tile.
- **The GeoPackage specification** (OGC 12-128r18), as the reader already
  follows it (`io/geopackage.py`): `gpkg_contents` carries an optional
  bounding box per table, `gpkg_geometry_columns` a declared geometry type
  name (which may be the generic `GEOMETRY`), and the R-tree extension a box
  per row. The extent is read from the R-tree, else from `gpkg_contents`.
- **OGC API Features**' collection description (extent, CRS, item count) is the
  same idea for a web service; nothing is taken from it beyond what STAC says.
- **XDG Base Directory Specification** (as in 55) for `register_at`.

Departures: no STAC documents and no STAC library, as in 55. The "fits" field
has no counterpart in STAC or OGC API Features: those describe data, not which
program options accept it. It is a statement about rasputin's own readers,
derived from their rules (section 4.4); no claim is made beyond that.

No novelty is claimed.

*Legacy.*

```
$ git grep -l -i -E 'list_layers|layer_names|feature_count|geometry_type|gpkg_contents|makedirs|mkdir' legacy-archive -- legacy
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/geometry.py
legacy-archive:legacy/tests/test_tin_repository.py
```

Read: the three hits are `mkdir` calls: `legacy/rasputin/application.py`
line 76, `tin_archive.mkdir(parents=True)`; `legacy/rasputin/geometry.py`
line 181, a `textures` folder; and the test's two. The legacy tree made its own fixed archive folders; it never
described a vector file. Nothing is carried across: the folder creation here
is for a path the caller names, under the existing boundary rules.

## 3. Ideas file

Read: `../rasputin_scratch/ideas/ideas.md`. Overlapping, left out:
"Structure for rasputin_data and rasputin_scratch" (registration, not a
folder layout, is how data is found); "Atomic mesh outputs" (this touches the
output path but not how a file is written); "Data-dependent tests never run in
a worktree"; "Intent layer over rasputin"; "`--stats` swallows the next option
as its path" (`--stats` becomes a writable path that allows `-`; the
ambiguity stays). Taken up from outside the file: the
three ideas from 55's code review (section 4.8), which the main session's
brief names but the file does not yet hold.

## 4. The design

### 4.1 Data flow

```
package sources.toml ─┐
user's sources.toml  ─┼─> sources.catalogue(user_file, plugin_files)
plugin sources.toml  ─┘     Catalogue(sources, origins, refused)
 (entry point group               │
  rasputin.sources)               │ Mapping[id, RemoteSource]
        ┌─────────────────────────┤
        v                         v
 cli: --dem KEY             cli `sources`
 (_dem_sources, as 55)      source_report.listing
                              raster: repository_for
                              vector: io/geopackage, io/geojson headers
                              facts ──> source_report.fits (pure)
                                 │
                                 v
                 each source's path, layer and fits (with class_map)
                                 │  the agent passes the path
                                 v
      --features, --lakes, --rivers, --stations, --reference PATH:
      the existing options and readers, unchanged
```

No second registry: a vector dataset is one more `kind = "local"` entry of the
same catalogue, read by the existing readers. No new module: the pure rules of
section 4.4 are a function in `source_report.py`, tested without a file.

### 4.2 A vector entry in the catalogue (`sources.py`)

```toml
[source.my-lakes]
kind = "local"
path = "/data/hydrography/lakes.geojson"    # .gpkg, .geojson or .json
credit = "..."
licence_note = "..."
description = "Lake polygons, 1:50 000, ..."

[source.my-land-cover]
kind = "local"
path = "~/data/landcover.gpkg"
layer = "landcover"            # a GeoPackage with more than one features table
class_map = "corine"           # optional: the map the file is coded in; its fits name only it
credit = "..."
licence_note = "..."
description = "Land-cover polygons, minimum mapping unit 25 ha, ..."
```

Changes to `RemoteSource`:

- `values: Literal["heights", "classes", "features"]`.
- `layer: str | None = None`.
- A `mode="before"` validator: a local entry whose `path` ends in `.gpkg`,
  `.geojson` or `.json` (case-insensitive; `VECTOR_SUFFIXES`, the suffixes the
  lake and feature readers share) and that states no `values` gets
  `values = "features"`. What the entry is comes from its file, so the user
  writes no redundant field.
- `_consistent`, each a `ValueError` naming the field:
  - `values = "features"`: `kind = "local"` needed; the path must have a
    vector suffix ("values: features need a .gpkg, .geojson or .json path");
    `z_unit` refused ("features hold no heights").
  - a vector suffix with `values` "heights" or "classes": refused ("values: a
    .gpkg path holds features").
  - `layer`: only with a `.gpkg` path.
  - `.gml` is not a vector suffix here: a GML file has no header to list, and
    `--features` still takes one as a path.
- Two `Literal`s widen to `["heights", "classes", "features"]`:
  `RemoteSource.values` (above) and `SourceReport.values`
  (`source_report.py`). `check_values`' body and its `wanted` parameter
  (heights or classes) stay as they are: with `RemoteSource.values` widened,
  `_dem_sources`' existing call refuses a features key given to `--dem` in
  its existing words ("--dem takes a source of heights; my-lakes holds
  features"), and `_class_source`'s refuses one given to `--land-cover`.
  No caller asks for features: the vector options take paths.

The catalogue's relative paths are resolved by one helper,
`_beside(file, value) = file.parent / Path(value).expanduser()`, used for
`path` and a path-like `class_map` in the user's file and in a plugin's file:
`~` is expanded (the third idea of section 4.8), and an absolute path after
expansion wins the join, as `Path` does.

### 4.3 What the listing reads from a vector file

`SourceReport` (in `source_report.py`) gains these fields; the raster fields
are unchanged and stay None for a vector source, and the vector fields stay at
their defaults for a raster:

```python
class ClassField(BaseModel):        # frozen
    attribute: str                  # the coded class map's attribute, as the file spells it
    values: tuple[str, ...]         # its distinct values, sorted as text, at most CLASS_VALUES
    distinct: int                   # how many distinct values there are

class Fit(BaseModel):               # frozen
    command: str                    # "mesh", "catchment", "station-catchments"
    option: str                     # "--dem", "--features", "--lakes", ...
    tool: str                       # the MCP tool: "station_catchments"
    argument: str                   # the MCP argument: "features", "lakes"
    class_map: str | None = None    # with --features: the --features-map to pass

# new SourceReport fields
origin: str                         # "package", "user", or "plugin <distribution>"
layer: str | None = None            # the GeoPackage layer read
geometry: tuple[str, ...] = ()      # GeoJSON type names present; "Geometry" for a
                                    #   GeoPackage column declared generic
features: int | None = None         # feature count
attributes: tuple[str, ...] = ()    # attribute names, file order, less key and geometry
classes: tuple[ClassField, ...] = ()
fits: tuple[Fit, ...] = ()          # every source, rasters included (section 4.4)
```

Read through the existing I/O modules (`io/repository.py` opens every file):

- **GeoPackage**, `io/geopackage.py` gains
  `layer_summary(conn, layer: GpkgLayer) -> LayerSummary` (declared geometry
  type name, feature count by `SELECT count(*)`, columns less the key and
  geometry column, extent: the R-tree's `min(minx), min(miny), max(maxx),
  max(maxy)` when the layer has one, else `gpkg_contents`' four bounds when
  all are set, else None) and `distinct_values(conn, layer, column, limit) ->
  tuple[tuple[str, ...], int]`. The layer is `layer_info(conn, entry.layer)`,
  so a file with two features tables and no `layer` is that entry's
  `problem` in the reader's own words ("name one layer of the features tables
  [...]").
- **GeoJSON**, `io/geojson.py` gains a pure `collection_summary(doc)`: the
  geometry type names, count, attribute names in first-seen order, the total
  bounds (`shapely.total_bounds` of the shapes), and whether a `crs` member is
  present; the CRS text comes from `read_collection(doc,
  default_crs=RFC7946_CRS)`, as the feature reader takes it.
- **Coverage and `box_covered`**: the extent box in the file's CRS goes through
  the same densify-and-reproject step as a raster tile: `_figures`' loop body
  becomes `_lonlat(crs, box) -> Polygon`, called per tile and once per vector
  source. For a vector source the box is the layer's extent, not its
  polygons; the listing's `note` says so in one sentence.
- **`classes`**: one `ClassField` per distinct attribute of the coded class
  maps (`ClassMap.codes` not empty) that fit the source (section 4.4).
  GeoPackage: `SELECT DISTINCT`; GeoJSON: from the parsed features.

Constants:

- `CLASS_VALUES = 100`, the most distinct values listed: a chosen value, above
  the 44 classes of CORINE level 3 and the class count of the shipped
  MapBiomas map. Checked at 34 distinct `code_18` values in 124,054 features.

Measured on Ola's files with `a0a1ad45` (scratch probes, warm cache):
`corine2018_dtm10_utm33.gpkg` (440 MB, two features layers): headers 0.001 s,
R-tree extent 0.027 s, `count(*)` 0.56 s, `SELECT DISTINCT code_18` 0.25 s
(34 values); `nve_hrd/rivers.geojson` (18.7 MB, 11,648 lines) parsed and
bounded in 0.24 s, `lakes.geojson` (7.9 MB, 81 polygons) in 0.11 s. A GeoJSON
is read whole, so a far larger one makes `sources` slower in proportion: a
stated limitation; no cap.

### 4.4 Which options a source fits (`source_report.fits`, pure)

```python
@dataclass(frozen=True, slots=True)
class Facts:
    values: Literal["heights", "classes", "features"]
    projected: bool | None          # heights: the CRS is projected; None when unknown
    suffix: str                     # ".gpkg", ".geojson", ".json"; "" for a raster
    geometry: frozenset[str]
    attributes: tuple[str, ...]
    crs_stated: bool                # GeoJSON: a crs member; GeoPackage: always
    class_map: str                  # the entry's, "" when none

def fits(facts: Facts, maps: Mapping[str, ClassMap]) -> tuple[Fit, ...]: ...
```

`maps` is every shipped map plus `property` (`class_maps.load` over
`shipped_names()`), passed in by `listing`; `fits` reads no file. Each
rule restates the acceptance condition of the reader the option calls, cited
here, so the derivation is from what the entry is and names no dataset:

| Fit | When | The reader it mirrors |
| --- | --- | --- |
| `mesh --dem`, `catchment --dem` | heights | `_dem_sources` |
| `station-catchments --dem` | heights, `projected` not False | `station_catchments`' refusal of a DEM in degrees |
| `mesh --land-cover` | classes | `_class_source` |
| `mesh --features`, one per map, `class_map` set | features, and the map's `attribute` is among the attributes (case-insensitive for `.gpkg`, exact for GeoJSON); only the entry's `class_map` when it names one | `query_features` lower-cases the column; `read_source` reads `p.get(attribute)` |
| `catchment --lakes` | features, geometry includes `Polygon`, `MultiPolygon` or `Geometry` | `read_lake_polygons` keeps polygons, skips the rest |
| `catchment --rivers`, `station-catchments --rivers` | GeoJSON, `crs_stated`, geometry exactly `{LineString}`, attribute `objectid` | `read_segments` |
| `station-catchments --lakes` | GeoJSON, `crs_stated`, geometry within `{Polygon, MultiPolygon}` | `read_nve_lakes` |
| `station-catchments --stations` | GeoJSON, `crs_stated`, geometry exactly `{Point}`, attribute `station` | `read_stations` |
| `station-catchments --reference` | GeoJSON, `crs_stated`, geometry within `{Polygon, MultiPolygon}`, attribute `station` | `read_references` |

A fit is a necessary condition read from the header, not a promise: the
reader may still refuse rows the header does not show (a missing or repeated
`objectid`, a repeated `station`, a polygon with no area). The `Fit`
docstring says so in those words. What the features mean (lakes or land
cover) is in the description, and the agent judges it. So a land-cover
GeoPackage also fits `catchment --lakes`: stated in the `Fit` docstring and the
instructions, not filtered by a guess. The `tool` and `argument` of a fit are
the MCP spellings, derived as `mcp/schema.py` derives them (`-` to `_`, no
leading dashes). An entry whose file was not read (missing, or a `problem`)
has no vector fits; a remote heights entry not yet cached still fits `--dem`
(the run says what to fetch).

### 4.5 The vector options take paths (unchanged)

Dropped in design review round 1: catalogue keys in the seven vector
options. The listing already gives each source's `path`, its `layer` and,
per fit, the `class_map`; an agent passes the path as a string, with the
layer as `--features-layer` or `--lakes-layer` and the fit's `class_map` as
`--features-map`. `cli.py`'s vector options, their help and the MCP tool
schemas are unchanged.

### 4.6 Output folders: made by the CLI, fenced by the MCP server

Ola's ruling (2026-10-16), on removing `--out-parent`: "We could remove it,
but we should have some kind of guard against writing everywhere"; on the
guard in the server with a fixed default root: "Yes, go with (a)". So the
fence moves from three CLI options to the one process an agent talks to, and
the CLI keeps only the checks every caller needs.

```
rasputin mcp [--out-root DIR]            DIR, else ~/rasputin-output
   │  build_server(app, out_root=...)    root = DIR.expanduser().resolve(), once
   v
call_tool: validate the arguments (as today)
   │  tool.writes(arguments)             schema.py, pure: (argument, path) for each
   │                                     value set on a writable path parameter
   │  _outside_root(writes, root)        server.py: await AsyncPath(p).resolve();
   │                                     not is_relative_to(root) -> refusal
   v  (refused: an error result, no child, nothing made, no job slot taken)
child: python -m tin_engine <argv>, env + RASPUTIN_MCP_OUT_ROOT=<root>
   │  cli_kit.destination(out, hint)     refuses a symlinked last part, resolves,
   v                                     makes the missing folders (the root too)
writer
```

**Which arguments are output paths.** The marker is Click's own: a path
parameter declared `writable=True` (`typer.Option(..., writable=True)`;
Click's `Path.writable`, which Typer passes through). `DerivedCommand` gains

```python
def writes(self, arguments: BaseModel) -> list[tuple[str, str]]:
    """(argument, path) for each value, given or default, of a parameter whose
    type is a writable path, options and positional arguments alike, every item
    of a repeatable one; a default is checked because the command writes there
    all the same. "-" is skipped where the parameter allows a dash (stdout, not
    a file)."""
```

derived from `command.params`, as `argv` is: no list per tool. Typer ignores
`writable` on a `str` annotation (probe: `StringParamType` has no
`writable`), so `--stats` and `--record` become `Path | None` with
`writable=True, allow_dash=True`; their `== "-"` tests compare to
`Path("-")`. The JSON schema is the same (`path` and `str` are both
`"string"`; probed on a toy app). Marked in the core:

| Tool | Arguments marked `writable=True` |
| --- | --- |
| `mesh` | `out`, `out_edges`, `stats`, `record_path` |
| `palette` | `out` |
| `catchment` | `out` |
| `station_catchments` | `out_dir` |
| `fetch_stations` | `out_dir` |
| `fetch` | `cache` (it writes tiles there; the read-only `--cache` of `mesh`, `catchment`, `station-catchments` and `sources` is not marked) |

A plugin's command is a tool of the same server, so its output options are
recognised the same way: the plugin declares them `writable=True`
(`docs/plugins.md`, *Output paths*). An unmarked option is not fenced; that
is stated there, and is the cost of deriving rather than listing.

**The check** (`server.py`, in `call_tool`, after validation and before the
job slot, so a refusal never waits behind a running job). Each path is
resolved with `AsyncPath.resolve()` (non-strict: a missing tail is kept),
which removes `..` and follows every symlink, the last part included, and is
judged against the root resolved once at build. A relative path is resolved
against the server's working directory, which the child inherits, so both
read it alike; `~` is not expanded, as the CLI does not expand it. The first
path outside is an error result, exit code None (as a validation error), in
these words:

```
out: /tmp/x/../escape/m.vtk resolves to /private/tmp/escape/m.vtk, outside the output root /Users/o/rasputin-output; write inside it (rasputin mcp --out-root sets the root)
```

The root itself may be a symlink (`~/rasputin-output` linked to a data disk):
it is resolved once, and paths through either spelling pass. Not checked,
stated limitations: a hard link inside the root to a file outside (rasputin
writes `.vtk`, `.ply` and JSON in place; an agent with only the MCP tools
cannot make a link); a symlink swapped in between the check and the child's
write (no tool makes one). The guard is about what an agent can name; a
plugin's own code can write anywhere, as any installed code can.

**`build_server(app, *, jobs=1, child=CHILD, out_root: Path | None = None)`**:
None means `Path.home() / DEFAULT_OUT_ROOT`, `DEFAULT_OUT_ROOT =
"rasputin-output"` (a chosen name). `rasputin mcp` gains `--out-root DIR`
(`Path | None`, help: "Tool calls may write only inside this folder, through
`..` and links too; missing folders inside it are made. Default:
~/rasputin-output.") and passes it to `serve`. Neither the given nor the
default root is made or checked at start: the first write makes it, through
`destination`, so an explicit root and the default behave alike, and a file
in the way is that first call's refusal naming the folder.

**Missing folders are made by the CLI**, in `destination`, not by the
server: the server then does no file writes beyond its sidecars, a shell
user gets the same folders, and a plugin writer through `destination` gets
them too. The containment check has already run, so only folders inside the
root are made from an MCP call. `fetch-stations --out-dir` already makes its
own (`mkdir(parents=True)`); `station-catchments --out-dir` makes its last
folder, `destination` its parents.

**`cli_kit.destination` becomes `destination(out: Path, hint: str = "--out")
-> Path`**: the final-part symlink refusal (kept: a link inside the root that
points elsewhere inside it passes containment and still overwrites a file the
caller never named), `resolve()`, then
`target.parent.mkdir(parents=True, exist_ok=True)`, an `OSError` a
`typer.BadParameter` naming `target.parent` and `hint`. `out_parent` goes
with the option, and `name` goes because the body never read it. It stays in
`plugin_api`; the change is a promised name changing, so the release that
carries it says so (`docs/plugins.md` table). No plugin caller is known:
`rasputin-radiation` is a copy of the core, not a plugin: it imports nothing from `plugin_api` (`grep -rn plugin_api`); the private package
is not on this machine, so it is assumed to have none.

**`--out-parent` is removed** from `mesh`, `catchment` and
`station-catchments`, with its parameter through `_report_target`,
`_record_target` and the seven `_destination` calls. From a shell there is
no fence: the shell user can write anywhere already.

**`sources` reports the root.** The server puts `RASPUTIN_MCP_OUT_ROOT=<root>`
in every child's environment (beside `CHILD_ENV`); `OUT_ROOT_VARIABLE` is
defined in `sources.py` (standard library only, so the server imports it
without pyproj). `listing(environ, ...)` already takes the environment and
sets the new `CatalogueReport.out_root: str | None` from it: the root under
the server, null from a shell. The variable only informs the listing; the
CLI fences nothing with it, and the help of `sources` does not offer it.

Accepted cost (as before): the folders are made when the path is checked,
before the run, so a run that fails later leaves an empty folder.

### 4.7 Plugins register sources

A plugin package declares an entry point in the new group `rasputin.sources`
whose value names a package (no `:attr`); the package holds `sources.toml`,
in the user file's format:

```toml
[project.entry-points."rasputin.sources"]
mydata = "my_plugin.data"
```

A sibling group rather than `rasputin.commands`: a plugin may ship data and no
command, and a sources file is found with no Typer app loaded. In `sources.py`
(standard library only, so the module's import rule holds):

```python
SOURCE_GROUP = "rasputin.sources"

@dataclass(frozen=True, slots=True)
class PluginFile:
    dist: str
    path: Path | None       # importlib.resources.files(value) / "sources.toml"
    error: str | None       # why there is none: an import error, an attr in the value, not on disk

def plugin_source_files() -> tuple[PluginFile, ...]: ...   # distribution-name order; never raises

def catalogue(user_file: Path | None, plugins: Sequence[PluginFile] = ()) -> Catalogue: ...
```

`Catalogue` gains `origins: Mapping[str, str]` ("package", "user", "plugin
<dist>") and `refused: tuple[str, ...]` (one line per plugin file skipped).
Order and clashes: package, then user, then plugins. The user's file is
refused whole on a clash with the package, as today. A plugin's file is
skipped whole, as one `refused` line, when it cannot be read or validated, or
when one of its ids is already taken (package, user, or an earlier plugin): a
plugin never replaces an entry and never stops the core, as `plugins.py` rules
for commands. `_known_sources` and `listing` pass `plugin_source_files()`;
`_known_sources` prints each `refused` line on stderr as `rasputin: <line>`,
the form `plugins._refuse` uses; the listing carries them in a new
`CatalogueReport.plugin_problems`.

### 4.8 The three ideas from 55's review

1. **`register_at`** is the XDG candidate (`$XDG_CONFIG_HOME/rasputin/sources.toml`,
   else `~/.config/rasputin/sources.toml`), the file a server started with no
   environment but `HOME` reads, whenever there is one; else the first
   candidate. When the file read is another (`$RASPUTIN_DATA/sources.toml`),
   the `note` says in one sentence that a server started without
   `RASPUTIN_DATA` reads `register_at` instead.
2. **Quoted fetch hints**: in `mesh`'s two hints (`--domain`, `--cache`, and
   the key) and `catchment`'s (`--cache` and the key), each value goes
   through `shlex.quote`, as `--out-crs` already does.
3. **`~` in a sources file**: `_beside` (section 4.2).

### 4.9 The MCP side

`mcp/server.py`'s `INSTRUCTIONS` gains, after the `sources` paragraph,
general guidance only (no dataset, country or resolution):

```
sources also lists vector datasets (polygons, lines, points) with their path, layer, geometry, feature count, attributes and, for land cover, the class values. Each source's fits names the tools and arguments it can be given to, with the class map to pass as features_map; what the features are is in its description. Pass its path to those arguments, its layer as features_layer or lakes_layer, and the fit's class map as features_map.
Output folders are created as needed, inside the output root that sources reports as out_root; give output paths in full: a path outside it, through .. or a link too, is refused.
```

`READ_ONLY`, `RunResult` and the schema derivation are unchanged; the tool
schemas lose `out_parent` (from `mesh`, `catchment`, `station_catchments`);
`test_mcp_schema.py`'s `MESH_SNAPSHOTS` hold none of the changed properties,
so they stand. The root is not in the instructions: they are fixed text, and
the listing carries it.

### 4.10 Docs

Written with this revision by `@architect`, except the last:

- `README.md`, *Connect an agent*: `--out-root` in the args example, the
  output root paragraph, a vector example table beside the raster one
  (generic names and paths), one sentence on `fits`.
- `docs/plugins.md`: *Ship datasets* (the entry point group, the package
  holding `sources.toml`, the clash rule); *Output paths* (`writable=True`);
  the `destination` row.
- `project_structure.md`: `sources.py`, `source_report.py`, `cli_kit.py`, and
  an entry for `mcp/`, which it lacked.
- `src_python/tin_engine/data/sources.toml`'s header comment is outside
  `@architect`'s write limit; `@developer` replaces its third sentence with:
  "Your own sources, local datasets among them (kind = "local", path: GeoTIFFs,
  or vector files .gpkg, .geojson or .json), go in your sources file in the
  same [source.<id>] form, under ids of their own: [the three places, as now].
  A plugin may ship a sources.toml of its own (entry point group
  rasputin.sources)."

### 4.11 Repeating Ola's Codex test without a shell

Not verified. Codex has no documented setting that removes its shell while
keeping MCP tools: openai/codex issue #6049 ("restrict to MCP tools") is open,
labelled enhancement. A commenter there reports that `shell_tool = false`
under `[features]` in `~/.codex/config.toml` removes the `exec_command` and
`write_stdin` tools; no maintainer confirms it. The same thread says
`--sandbox read-only` "gates tool effects but doesn't remove tools", so it is
not the switch. Proposed protocol, for Ola: set `[features] shell_tool =
false`, start a fresh session, ask Codex first to list its tools (the run is
void if a shell tool is listed), then repeat the original task; pass when it
finds the DEM, lakes and land cover through `sources` and runs catchment then
mesh with an output folder that does not exist yet. The automated
counterpart is T7.

## 5. Deleted

- `destination`'s refusal "is not an existing directory", replaced by the
  creation of the folders (section 4.6).
- `--out-parent` on `mesh`, `catchment` and `station-catchments`, and
  `destination`'s `out_parent` and unused `name` parameters, replaced by the
  server's output root (section 4.6).
- `palette`'s unchecked `out.write_text`, replaced by `destination`.
- `listing`'s rule `register = user or looked[0]`, replaced by section 4.8.1.

## 6. Left out

- Local class rasters by key (`kind = "local"`, `values = "classes"`): still
  refused with the pointer to `--land-cover PATH`.
- GML entries; layers of other formats (Shapefile, FlatGeobuf): no reader.
- A source's credit in the catchment file or the mesh record when a vector file the catalogue lists is passed: the class map's notice already reaches the mesh record, as today; carrying vector credits is its own change.
- Catalogue keys in the vector options (section 4.5): the path is passed.
- Footprints finer than the extent box; feature counts inside `--bbox`.
- Telling lakes from land cover in `fits`: the description says it.
- `station-catchments --only`, `--reference` beyond the key: unchanged.
- A fence for shell users: `--out-parent` is gone, and `RASPUTIN_MCP_OUT_ROOT`
  fences nothing in the CLI (section 4.6).
- Relative output paths taken from the root (running the child in it): a
  relative input path would then change meaning, and the server's sidecars,
  `_read` and `_outputs` resolve against its own working directory
  (`mcp/server.py`:176-211).
- Hard links and a symlink swapped in after the check (section 4.6).

## 7. Tests

Every test is `@tester`'s to write. Fixtures were measured with the main
checkout at `a0a1ad45` (the same tree as this branch) by the scratch probes
named, through the CLI and the library.

**Fixture L** is 55's (`dataset_fixtures.write_local`: the bowl, four tiles,
EPSG:25833, coverage lon 15.0 to 15.351879, lat 59.179575 to 59.538349).
**Fixture V, vector sources over it**:

- `lakes.geojson`: `dataset_fixtures.lake_file`, one Polygon
  (`catchment_fixtures.lake_box()`, x 509650 to 510350, y 6579650 to 6580350),
  a `crs` member EPSG:25833, attribute `name`. Measured: coverage (15.169718,
  59.355479, 15.182063, 59.361781); `box_covered` 0.031 for (15.15, 59.33,
  15.20, 59.38), 0.0 for (15.2, 59.3, 15.6, 59.4).
- `landcover.gpkg`: `gpkg_fixtures.write_gpkg` with two features layers in
  EPSG:25833, `landcover` (rows: 1, the box x 500000 to 520000, y 6560000 to
  6600000 less the lake, `Code_18` "311"; 2, the lake box, "512") and
  `tiles` (one row), so the entry needs `layer = "landcover"`. The fixture
  writer declares `GEOMETRY` and leaves `gpkg_contents`' bounds unset, so the
  extent comes from the R-tree: (500000, 6560000, 520000, 6600000). Measured:
  coverage (15.0, 59.178672, 15.353647, 59.538349); `box_covered` 1.0,
  0.3792 and 0.0 for 55's three boxes; columns `OBJECTID`, `Shape`, `Code_18`,
  so `attributes` ("Code_18",); distinct values "311", "512".
- Measured end to end (probe `p57/probe.py`, paths, not keys): `catchment
  --dem <L> --seed <SEED> --seed-crs EPSG:25833 --lakes lakes.geojson`: exit 0,
  766 nodes, 7.655000 km²; `mesh --dem <L> --domain <that file> --tolerance 5
  --features landcover.gpkg --features-layer landcover --features-map corine
  --out m.gpkg`: exit 0, 1083 triangles, largest difference 4.6722 m, record
  `features` "landcover.gpkg layer landcover, class map corine: 2 features, 2
  lines, 12 vertices"; triangles by `code`: 311 1004 (7,172,482 m²), 512 79
  (482,518 m²). Each run under 0.11 s.

Isolation as in 55 (`conftest.py` gives each test an empty `XDG_CONFIG_HOME`;
no test reads `~/.config/rasputin`). Plugin tests use `plugin_fixtures`'
installed-looking distribution with a `[rasputin.sources]` entry point and a
package directory holding `sources.toml`.

- **T1 catalogue rules** (`sources.py`): a `.gpkg`, `.geojson` or `.json` path
  with no `values` loads as features; each refusal of section 4.2 names its
  field; `layer` refused off a `.gpkg`; `~/x.gpkg` in `path`, and
  `~/m.toml` in `class_map`, expand under a monkeypatched `HOME`; a relative
  path is still taken from the file's directory; `mesh --dem <features key>`
  is refused with `check_values`' words of section 4.2.
- **T2 plugin files**: `plugin_source_files` finds the fake distribution's
  file; its entries list with `origin` "plugin fakeplug"; a plugin id equal to
  a package id, a user id or an earlier plugin's skips that plugin's whole file
  with one `refused` line naming the distribution and the id, the rest of the
  catalogue unchanged; an entry point value with `:attr`, or one whose package
  fails to import, is a `refused` line, never an exception; `mesh` with such a
  plugin prints the line on stderr and runs.
- **T3 vector headers** (`io/geopackage.py`, `io/geojson.py`): on fixture V,
  `layer_summary` gives `GEOMETRY`, 2, the R-tree extent and the columns less
  `OBJECTID` and `Shape`; with the R-tree dropped and `gpkg_contents`' bounds
  set, those bounds; with neither, None; `distinct_values` gives ("311",
  "512"), 2, and with `limit=1`, ("311",), 2; `collection_summary` of the lake
  file gives {"Polygon"}, 1, ("name",), the lake box, crs member present.
- **T4 fits, pure** (`source_report.fits`): a table test over `Facts`, one row per
  rule of section 4.4 and one row per exclusion (a `.gpkg` for `--rivers`; a
  GeoJSON with no `crs` member for `--stations`; `MultiLineString` for
  `--rivers`; a geographic DEM for `station-catchments --dem`; an entry
  `class_map` whose attribute is absent: no `--features` fit). The MCP
  spelling of each fit equals `mcp.schema.derive(cli.app)`'s tool and
  argument names for that command and option.
- **T5 fits agree with the readers**: for each fixture file (fixture V's two,
  plus `nve`-style stations, rivers and reference GeoJSON written by the
  existing station fixtures, and one GeoJSON without a `crs` member): where
  the listing names no fit, the option's reader refuses the file; where it
  names one, the reader accepts that fixture. A fit is a header-level
  necessary condition, so the test does not claim the converse for every
  file: the readers also refuse a missing or repeated `objectid`, a repeated
  `station`, a polygon with no area, which no header shows. This is the
  guard against drift between section 4.4 and the readers.
- **T6 the listing** (`rasputin sources --bbox 15.15 59.33 15.20 59.38` over
  fixture V and L registered in the test's user file): the figures of fixture
  V above; the landcover entry's `classes` one `ClassField` ("Code_18",
  ("311", "512"), 2) and `fits` holding `mesh --features` with `class_map`
  "corine" and with "corine-water", and `catchment --lakes`; the lakes
  entry's `fits` holding `catchment --lakes` and `station-catchments --lakes`;
  fixture L's entry fitting `--dem` on all three commands; every entry's
  `origin` "user", the shipped ones "package"; the `note` sentence on vector
  extents present. The landcover entry without `layer`: `problem` names both
  tables, the other entries unaffected.
- **T7 acceptance over MCP only** (`build_server(cli.app)` with the production
  child, in process; `HOME` a temp directory whose
  `.config/rasputin/sources.toml` registers fixture L, the lake file and the
  land-cover GeoPackage under ids that name nothing real; `RASPUTIN_DATA` and
  `XDG_CONFIG_HOME` unset). After that setup the test acts only through tool
  calls, on a server built with **no `out_root`**, so the default root
  `<HOME>/rasputin-output` applies and is absent before the first call:
  `sources` with `bbox` (15.15 59.33 15.20 59.38); it reads `out_root` off
  the listing (equal to the resolved `<HOME>/rasputin-output`), takes the DEM
  id from the entry fitting `catchment` `dem` with `box_covered` 1.0, and
  the lake source's `path` and the land-cover source's `path`, `layer` and
  fit `class_map` from the entries whose `fits` name `lakes` and `features`;
  `catchment` with `dem=[dem id]`, `seed`, `seed_crs`, `lakes=<lake path>`,
  `out` = `<out_root>/new/deeper/c.geojson`: exit 0, the output listed;
  `mesh` with `dem=[dem id]`, `domain` the catchment's output path taken from
  the tool result, `tolerance` 5, `features=[land-cover path]`,
  `features_layer=[its layer]`, `features_map=[the fit's class_map]`, `out` =
  `<out_root>/meshes/run1/m.gpkg`: exit 0, the run record (structured data)
  naming layer `landcover` and map `corine`. The oracle may then read the
  files: 766 nodes, 1083 triangles, codes 311:1004 and 512:79 (the figures
  measured above; the folders do not change them).
- **T8 output folders and the fence.** Replaces the eight red tests at
  `08cb70c2` that named `out_parent` or expected it (TestOutputFolders'
  three, the two palette tests, T7, T10's two): the palette tests, T7 and
  T10 stay with the changes stated here; the rest are below. Every refusal
  asserts that nothing was made (`sorted(tmp_path.rglob("*"))` before and
  after) and, through the server, that the child never ran (`mcp_driver.fake_child` with an `argv_file`,
  which stays absent).
  - **T8a `destination`** (`test_cli_kit.py`): a missing chain is made and
    the file is not; an existing folder is left as it is; `tmp/a/../b/x.ply`
    makes `tmp/b`, not `tmp/a`; a file in the way is refused naming the folder,
    hint `--out`; a symlinked `--out` whose target's folder is absent is
    refused before any folder is made. `palette corine --out <absent>/p.json`
    writes; `palette` refuses a symlinked `--out` (both as at `08cb70c2`).
    `TestOutputFolders`' three `out_parent` tests are deleted with the
    parameter.
  - **T8b recognition, pure** (`test_mcp_schema.py`): on a toy app, `writes`
    returns a writable option's value, each item of a repeatable one and a
    writable positional argument's, and the default of one left unset (a plugin's default outside the root is refused like a given path); skips one whose value is None, `-` on a parameter
    that allows a dash, and a path parameter not marked writable. On
    `cli.app`, the set of (tool, argument) with a writable path equals the
    table of section 4.6 (a drift guard: the list is the test's, not the
    code's), and no tool has an `out_parent` argument.
  - **T8c the fence** (`test_mcp_vector_sources.py`, `build_server(app,
    out_root=<tmp>/root, child=fake_child(...))`, the root absent): `mesh` with `out`
    `<root>/a/m.vtk` runs the child; refused, each naming the argument and the
    root: `out` = `<root>/../escape/m.vtk`; `out` through a symlinked folder
    `<root>/link` (made by the test) pointing at `<tmp>/outside`; `out` =
    `<root>/l.vtk`, a symlink to `<tmp>/outside/real.vtk`; `record_path`,
    `out_edges` and `stats` each outside with `out` inside. Passes: `stats`
    `-`; a relative `out` with the working directory (`monkeypatch.chdir`)
    inside the root; a root given as a symlink to `<tmp>/real` with `out`
    spelled through either name. `fetch` with `cache` outside is refused.
  - **T8d a plugin's output** (`test_plugins.py`, `plugin_fixtures`' fake
    distribution gaining a command with a `writable=True` `--out` that writes
    a line): over `build_server(plugins.application(...), out_root=<root>)`,
    inside writes, outside is refused and the child never runs.
  - **T8e the root reported**: `sources` over a server with
    `out_root=<tmp>/root` lists `out_root` as its resolved path; `rasputin
    sources` from the CLI with no `RASPUTIN_MCP_OUT_ROOT` lists `out_root`
    null; `rasputin mcp --out-root <dir>` (through `mcp_stdio_launcher.py`, the
    real stdio entry) reports it the same way; `mesh --out-parent x` is a
    usage error (no such option).
  - Existing tests that name `out_parent` go with the option:
    `test_mcp_server.py` (`test_out_parent_passes_through...`, the
    `with_out_parent` parametrisation, `test_out_outside_out_parent_...`),
    `test_cli_mesh_stats.py` (`test_out_parent_applies_to_the_report`,
    `test_stats_outside_out_parent`), `test_cli_mesh_plain_output.py`
    (`test_outside_out_parent`). Every other `build_server` call that writes
    outside the default root passes `out_root=tmp_path` (`mcp_driver.py`,
    `test_mcp_server.py`'s builder, `mcp_stdio_launcher.py`,
    `test_plugins.py`, `test_mcp_datasets.py` where it writes).
- **T9 the review ideas**: `register_at` is the XDG file with
  `RASPUTIN_DATA` set and its file present, and the `note` says so; the
  fetch hints of `mesh` (missing cache with `--domain` and `--cache` holding a
  space) and `catchment` (fixture C of 55, its seed block row skipped, a
  cache path with a space) parse back with `shlex.split` to the same values.
- **T10 instructions**: the initialise result contains "fits", "Pass its
  path", "Output folders are created as needed" and "out_root", and no
  source id.
- Unchanged: `test_mcp_schema.py`'s snapshots. 55's `register_at`
  assertions (`test_cli_sources.py`, `test_mcp_datasets.py`) all run with
  the XDG file as the only candidate, so they hold unchanged.

No mutation round: no topology or invariant-critical code changes.

## 8. Size

Net production lines (`CLAUDE.md` §2; TOML and docs not counted):

| Part | Lines |
| --- | --- |
| `sources.py`: features value, `layer`, before-validator, rules, `_beside`, plugin files, `catalogue` with origins and refusals | 70 |
| `io/geopackage.py` `layer_summary`, `distinct_values`; `io/geojson.py` `collection_summary` | 45 |
| `source_report.py`: vector branch, `_lonlat` split, `ClassField`, new fields, `register_at`, notes, plugin problems; `Facts`, `Fit`, `fits` | 125 |
| `cli.py`: plugin files to `_known_sources`, refusals on stderr, quoted hints, palette through `destination` | 15 |
| `cli_kit.destination`: containment out, folders made | 0 |
| `cli.py`: `--out-parent` out of three commands and two helpers; `writable=True` and `Path` on the marked options; `mcp --out-root` | -10 |
| `mcp/schema.py` `writes`; `mcp/server.py` root, fence, child variable, instructions | 35 |
| `sources.py` `OUT_ROOT_VARIABLE`; `source_report.py` `out_root` | 4 |
| **Total** | **about 315**: 286 measured at `08cb70c2` before the folder part (`python3 tools/count_loc.py da133bce 08cb70c2`), plus about 30 for it; worst case 400 |

Round 1 of the design review cut about 40 lines: `_vector_source` and the
seven options' retyping (section 4.5), and a module's overhead.

One PR, the full path (over 50 lines). It touches no refine or mesh code
(`include/terrain/`), so no `@perf` acceptance run.

## 9. Files

`src_python/tin_engine/sources.py`, `source_report.py`, `io/geopackage.py`, `io/geojson.py`, `cli.py`, `cli_kit.py`,
`mcp/server.py`, `mcp/schema.py`, `data/sources.toml` (comment); `README.md`, `docs/plugins.md`,
`project_structure.md`, `ROADMAP.md`; tests under `tests/python/` (T1 to T10,
the fixture helpers).

## Review

- **Design review, round 1, 2026-10-16.** Commit `46edf5b6`, 0 counted lines, @reviewer: CHANGES REQUESTED, applied by @architect: catalogue keys in the seven vector options dropped (section 4.5; the listing gives `path`, `layer` and each fit's `class_map`; old T7 cut); `source_fits.py` folded into `source_report.py`; section 4.2 names the two widened `Literal`s; T5 and the `Fit` docstring call a fit a header-level necessary condition. Held: one catalogue, readers reused, fits per option, the `rasputin.sources` group, the CORINE figures, both folder defaults.
- **Design review, round 2, 2026-10-16.** Commit `8f440523`, 0 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: section 6 credit bullet ("a vector key" becomes "a vector file the catalogue lists is passed"). Held: only `RemoteSource.values` and `SourceReport.values` widen; `check_values` refuses a features key for `--dem` and `--land-cover` before the path is used.
- **Design review, round 3, 2026-10-16.** Commit `28547032`, 0 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: `rasputin-radiation` is a copy, not a plugin; instructions say to give output paths in full; relative paths stay against the server's working directory (sidecars and output reading depend on it). Held: Typer keeps `writable` and `allow_dash` on a `Path` option with the same schema; only `fetch` writes the cache. Open to Ola: holding `fetch --cache` inside the root.
- **Code review, round 1, 2026-10-16.** Range `da133bce..9f634f94`, 342 counted lines, @reviewer: CHANGES REQUESTED: `DerivedCommand.writes` skipped a writable option left unset, so a plugin's default outside the root was written unchecked (fixed: red `dd11957f`); README states the default tile cache as the guard's one exception (`095ffa1b`). Held: `..`, symlinks, `-`, `fetch --cache`, refusal before any child; fits per option; `--out-parent` gone. `DemInput.label` goes in a follow-up (two plugin-visible signatures).
- **Code review, round 2, 2026-10-16.** Range `9f634f94..8648c039`, 340 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: the design sketch's `writes` docstring (defaults checked). Held: the hole is closed for repeatable and positional parameters; all nine core writable parameters default to None, so no core call is newly refused. Idea for later: a plugin default computed at run time could differ from the one checked.
