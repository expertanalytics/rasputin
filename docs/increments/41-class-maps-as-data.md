# Increment 41: land-cover class maps and palettes as data files

**Status: designed** by `@architect`, 2026-10-12, on `f09b9863` (branch
`worktree-class-maps`), while Ola was away; the defaults taken are listed in
section 9. Next: `@reviewer` on the design, then `@tester`'s red suite
(section 7). Python only, one PR.

## 1. What Ola asked, and what this fixes

The generality audit (`docs/increments/generality-audit.md`, rows B1, B2, B3,
B6, B8) found that rasputin knows one land-cover scheme, CORINE, and knows it
as Python code:

- **B1.** The class maps are three Python objects, all CORINE, with CORINE's
  water codes and copyright notice
  (`src_python/tin_engine/feature_input.py@f09b9863:59-66`,
  `src_python/tin_engine/feature_input.py@f09b9863:106-123`). Any other scheme
  (MapBiomas, ESA WorldCover, the US National Land Cover Database, Norway's
  AR5) needs a pull request.
- **B2.** One palette, CORINE's, as a Python table
  (`src_python/tin_engine/palettes.py@f09b9863:15-65`); `rasputin palette`
  knows only `corine` (`src_python/tin_engine/cli.py@f09b9863:1491-1512`).
- **B3.** The clean-up defaults, repair 1 m, simplify 50 m, outline snap 5 m
  (`src_python/tin_engine/cli.py@f09b9863:833-872`), were tuned on CORINE over
  a 10 m DEM and are wrong for a 1 m DEM with fine land cover or a 1 km DEM.
- **B6.** A class code must be an integer
  (`src_python/tin_engine/feature_input.py@f09b9863:579-588`), so a scheme with
  text codes ("Forest") cannot label triangles.
- **B8.** The edge vocabulary (river, road, ..., land_cover, water) is fixed
  (`src_python/tin_engine/features.py@f09b9863:176-190`); a user with "levee"
  or "canal" lines must rename them to fit.

Ola's rulings, 2026-10-12. To "Should land-cover class maps and palettes
become TOML data files a user points at with an option, with the CORINE ones
shipped inside the package, instead of Python code?" he answered "4: yes". On
the clean-up defaults: derive them from the land-cover source's stated
positional accuracy, with the DEM cell as a floor, and keep the command-line
options as overrides ("On 3: Probably, but we should have the option of
override").

What the user sees after this increment:

```bash
rasputin mesh ... --features clc.gpkg --features-map corine          # as today, same bytes
rasputin mesh ... --features lcm.gpkg --features-map ./uk-lcm.toml   # any scheme, by path
rasputin palette corine            # as today, same bytes
rasputin palette ./uk-lcm.toml     # the preset of any class-map file
```

## 2. Prior art: legacy and literature

*Legacy.* `git grep -l -i -E 'land_cover|landcover|globcov|corine|palette|colou?r' legacy-archive -- legacy`
(tests and Markdown dropped) returned `legacy/rasputin/application.py`,
`geometry.py`, `globcov_repository.py`, `gml_repository.py`,
`land_cover_repository.py`, `material_specification.py`, `mesh_colouring.py`,
`mesh_utils.py`, `tin_repository.py`, `web_visualize.py`, `wfs_repository.py`.
The relevant two: `land_cover_repository.py` declares an abstract
`LandCoverMetaInfoBase` with `color`, `describe` and `material` per class of a
per-scheme `Enum`, and `globcov_repository.py` implements it for GlobCover 2009
with the colours as a Python dict inside a method. That is the shape this
increment wants, a scheme as a swappable component, but held as code: one
subclass per scheme. Carried across: the intent (a scheme names, colours and
classifies its classes); not the code, because a data file needs no subclass
and no pull request. `material` (rendering properties) is not carried; nothing
in rasputin renders.

*Literature.* No novelty is claimed. Three pieces of practice the design leans
on:

1. **Legends as data.** Land-cover products ship their legend as a table:
   MapBiomas publishes `legend_code_mapbiomas_brazil_collection_11.csv`
   (code, label, colour; `brasil.mapbiomas.org/en/codigos-de-legenda/`, read
   2026-10-12), QGIS stores categorised styles as data (`.qml`, and the
   GeoPackage `layer_styles` table increment 36 already writes). A class-map
   file is the same: a legend, plus which edges a class's borders carry.
2. **String codes through a crosswalk.** Translating one scheme's class labels
   to another's codes is the "legend translation" of the FAO Land Cover
   Classification System (Di Gregorio and Jansen, *LCCS*, FAO 2000). The
   `[codes]` table (section 4) is the simplest crosswalk: text to integer.
3. **Detail no finer than the data's positional accuracy.** Tobler's rule
   ("Measuring spatial resolution", Proc. Land Resources Information Systems
   Conference, Beijing, 1987): detail below a source's resolution is not
   information. CORINE states a geometric accuracy better than 100 m, a
   minimum mapping unit of 25 ha and a minimum width of 100 m (EEA, *CLC2018
   Technical Guidelines*, 2017; Copernicus Land Monitoring Service product
   page for CLC2018, read 2026-10-12). Today's 50 m simplification is half
   that accuracy.

**Departure, named.** The three ratios of section 5 are not from the
literature. Only "simplify by half the stated accuracy" has the backing of
point 3; the repair ratio (1/100) and the snap ratio (1/20) are 20c-3's
tuning on CORINE (`docs/increments/20c-soft-quality.md`, M8, for the repair;
`docs/increments/32-landcover-simplify.md` for the simplification), restated as
fractions of CORINE's accuracy so that a CORINE run is unchanged. They are
a scaling rule, not a derived optimum.

## 3. The ideas file

Taken up: *Product review (b), generality audit* (rows B1, B2, B3, B6, B8 and
Ola's two rulings). Overlaps: *Repair CORINE once per area of interest*: that
cache's key must carry the clean-up distances used, which after this
increment are derived, not fixed; the cache design reads them from the run
record (they are recorded, section 5.4). Left out, handed back: none new.

## 4. The class-map file

One TOML file is one map. Read with the standard library's `tomllib` (Python
3.11+; `pyproject.toml` requires 3.12), so no new dependency. Validated by a
Pydantic model with `extra="forbid"`, so a misspelt key is refused naming the
file and the key.

```toml
name = "corine"                                # what the mesh file's record says
attribute = "Code_18"                          # the feature attribute read
code_system = "CORINE Land Cover level-3 code" # present: a coded map, labels triangles
positional_accuracy_m = 100                    # required with code_system
notice = "Contains modified CORINE Land Cover 2018 data ..."   # verbatim today's
licence = ""                                   # optional; written only when not empty
preset = "rasputin CORINE natural"             # ParaView preset title
unlisted = "land"                              # a value not listed: "land", "drop" or "refuse"

[edges]                       # the vocabulary names a class's borders carry
water = ["land_cover", "water"]
land = ["land_cover"]         # or "drop": no constraint

[no_class]                    # code 0: no polygon, and every constraint line
label = "No polygon; constraint lines"
colour = "#3c3c3c"

[[class]]
code = 111
label = "Continuous urban fabric"
colour = "#8c5f5a"
water = false                 # optional, default false; picks [edges].water or .land
# edges = ["land_cover"]      # optional per-class override of the two above

[codes]                       # optional (B6): attribute text -> integer code
# "Forest" = 1

[vocabulary]                  # optional (B8): names beyond the default vocabulary
# extra = ["levee", "canal"]
```

**Two kinds of map.** A map with `code_system` is *coded*: every `[[class]]`
has `code` (1 to 2^31 - 1, unique), `label` and `colour` (`#rrggbb`, lower
case); its polygons label triangles and are cleaned up. A map without
`code_system` is *uncoded*, the shape of today's `property` map: every
`[[class]]` has `value` (the attribute text) and `edges`; no label, colour or
accuracy. The loader refuses a field of the other kind.

**How it becomes today's `ClassMap`.** The loader (`tin_engine/class_maps.py`,
new) builds the existing frozen `ClassMap`
(`src_python/tin_engine/feature_input.py@f09b9863:75-90`), so
`open_features` and everything after it keep their shape:

| file | `ClassMap` field |
|---|---|
| each class, keyed `str(code)` (or each `[codes]` string whose code is a listed class) | `classes[key]` = the class's `edges`, else `[edges].water` or `.land`; `"drop"` becomes `None` |
| `unlisted = "land"` / `"drop"` / `"refuse"` | `otherwise` = `[edges].land` (or `"drop"` when that is `"drop"`) / `"drop"` / `"refuse"` |
| `code_system`, `notice`, `name`, `attribute` | `codes`, `notice`, `name`, `attribute` |
| `[no_class]` and the classes, in file order | `palette: Mapping[int, tuple[str, str]] \| None` (the table itself, replacing today's palette *name*) |
| `preset`, `licence`, `positional_accuracy_m`, `[codes]`, `[vocabulary].extra` | new fields `preset`, `licence`, `accuracy_m`, `string_codes`, `extra_names` |

Two small changes in `feature_input.py`: `_mask` returns None for a key mapped
to `None` (a listed class that is dropped, which `corine-water` needs for its
land classes); `_code` looks the value up in `string_codes` when the map has
them, and otherwise parses it as today. An unlisted text with no integer code
is refused naming the file, the feature and the value, unless `unlisted =
"drop"`.

**Why CORINE stays byte-identical in the mapping.** Today `corine` keys only
the five water codes, as strings, and gives everything else `("land_cover",)`
(`src_python/tin_engine/feature_input.py@f09b9863:93-104`). The file keys all
44 classes: water classes to `land_cover, water`, the rest to `land_cover`,
and `unlisted = "land"` gives an unlisted value (CORINE's `999`, a stray
`"0511"`) `land_cover` as today. The key is still `str(value)`, so `"0511"` is
unlisted in both. Same masks, same codes, for every value.

**Shipped maps,** under `src_python/tin_engine/data/class_maps/`, read with
`importlib.resources` as `fetch/nve.py` reads its packaged list
(`src_python/tin_engine/fetch/nve.py@f09b9863:62`):

| name | attribute | edges | note |
|---|---|---|---|
| `corine` | `Code_18` | water: land_cover + water; land: land_cover | today's `corine` |
| `corine-water` | `Code_18` | water: water; land: drop; unlisted drop | today's `corine-water` |
| `clc18_kode` | `clc18_kode` | as `corine` | today's `clc18_kode` |
| `mapbiomas-c11` | `class` | water (classes 26, 33, 31): land_cover + water; land: land_cover | new; the legend CSV's codes, labels and colours; accuracy 30 m |

The three CORINE files repeat the 44-class table (TOML data is not counted
production code, `CLAUDE.md` §2); a test pins that the three tables are equal,
so a colour changed in one is caught. No `extends` key: one file stays
readable on its own, which is what a user copying `corine.toml` to start a
scheme of their own needs.

`property` stays in code: it is the identity map over the vocabulary's names,
so it is built from the run's vocabulary (with any extras, section 6), not
shipped.

**MapBiomas.** Increment 26b's `--land-cover mapbiomas-c11` writes MapBiomas
class numbers as `land_cover_code` but has no names or colours
(`docs/increments/26-land-cover.md`, "Not in scope"). The catalogue entry
gains `class_map: str = ""` (`src_python/tin_engine/sources.py@f09b9863:133-155`;
still Pydantic-only), set to `"mapbiomas-c11"`; the `--land-cover` route then
passes that map's palette to the GeoPackage writer, as the CORINE route does
(`src_python/tin_engine/cli.py@f09b9863:1210`). `rasputin palette
mapbiomas-c11` gives the ParaView preset. Colours: MapBiomas's own legend
colours, not natural colours (default D5). The `.vtk` and `.ply` of a
`--land-cover` run do not carry a palette, so they are unchanged.

**The option.** `--features-map` keeps its name (the audit's, and the one
users have). Its value is a shipped name or a path: a value ending in `.toml`
or containing a path separator is a path; anything else is a shipped name or
`property`. An unknown name is refused listing the shipped names; a missing
or unreadable path is refused naming it; a TOML syntax error is refused with
`tomllib`'s line and column. `rasputin palette` takes the same: a shipped
name or a path, and refuses an uncoded map ("has no classes to colour").
The help texts list the shipped names from the package directory, not from a
literal.

## 5. Clean-up distances derived from the source's accuracy

### 5.1 The formulas

With `a` the coded source's stated positional accuracy (`positional_accuracy_m`)
and `h` the DEM cell (section 5.2), the length scale is

    s = max(a, h)

and the three defaults are

    repair         = s / 100
    simplification = s / 2
    outline snap   = min(s / 20, 50)

Division, not multiplication by 0.01 or 0.05, so CORINE's values are exact in
floating point.

**CORINE** (`a = 100`), on any DEM with a cell up to 100 m: `s = 100`, so
repair 1 m, simplification 50 m, snap 5 m, today's three numbers exactly
(checked: `100.0/100 == 1.0`, `100.0/20 == 5.0`, `100.0/2 == 50.0` in
CPython 3.12). Every CORINE case in use has a cell of at most 100 m:
Numedalslagen and the test fixtures at 10 m (`DELTA_X = 10.0`,
`tests/python/geotiff_fixtures.py@f09b9863:99`), Lagan on GLO-30 resampled to
EPSG:3006 at 31 m (`default_spacing` rounds the 1-arc-second north-south step,
`src_python/tin_engine/target_grid.py@f09b9863:91-101`), the Swedish SE2
catchments of increment 39 at 10 m.

Other cases, for scale: MapBiomas over ANADEM at 30 m (`a = 30`, `h = 30`):
0.3, 15, 1.5 m. A 1:5000 cadastre (`a = 2.5`) over a 1 m DEM: 0.025, 1.25,
0.125 m. CORINE over a 1 km DEM: `s = 1000`: 10, 500, 50 m.

**The snap cap.** The outline snap must stay under the read region's 100 m
margin (`MARGIN`, `src_python/tin_engine/feature_input.py@f09b9863:56`;
refused at or above it by `snap_to_outline`). `s / 20` reaches 100 at
`s = 2000`; the derived value is capped at 50 m, half the margin, which binds
from `s = 1000` (a 1 km DEM). Scale assumed: metres of a projected CRS;
checked against the margin only, at cells up to 2 km. A user's
`--features-outline-snap` keeps today's check (under 100).

**The floor.** "The DEM cell as a floor" is read as: the length scale is
never finer than one cell. A source more accurate than the DEM cannot place a
border more precisely than the mesh resolves it. The floor never binds for
CORINE on any DEM in use (above).

### 5.2 What `h` is

The larger of the two node spacings of the tile the features meet,
`max(meta.delta_x, |meta.delta_y|)` of `opened.tile`
(`src_python/tin_engine/cli.py@f09b9863:1113`), in the mesh CRS. That tile is
always projected: a geographic DEM needs `--out-crs`
(`src_python/tin_engine/cli.py@f09b9863:770-778`), and then the tile is the
resampled target grid. Units are the CRS's; the clean-up distances are in the
same units already (a CRS in feet is the audit's row H1, not this increment).

### 5.3 Several sources, none coded, and overrides

- **Several coded sources**: `a` is the largest of their accuracies (the
  coarsest source sets the scale). One mesh has one code system already
  (`src_python/tin_engine/cli.py@f09b9863:1856-1866`), so in practice they are
  one scheme.
- **No coded source**: the clean-up never runs (it acts on coded polygons
  only, `src_python/tin_engine/feature_input.py@f09b9863:363`), but the record
  writes the three values for every features run
  (`src_python/tin_engine/cli.py@f09b9863:1197`). `a = 100` there, so those
  runs record today's values.
- **Overrides**: `--features-repair`, `--features-tolerance` and
  `--features-outline-snap` default to `None` (given: used as is, with
  today's checks; not given: derived). Each replaces only its own distance.
  The help texts say "Default: from the class map's positional accuracy and
  the DEM cell (CORINE: 1)", and so on.
- The 15° start angle for simplified land cover
  (`src_python/tin_engine/cli.py@f09b9863:1120-1124`) reads the resolved
  simplification, not the option.

### 5.4 Where it is computed

A pure function in `class_maps.py`:

```python
def derived_cleanup(accuracy_m: float, cell: float) -> CleanUpDistances: ...
```

`CleanUpDistances` is a frozen model of the three numbers. The CLI calls it
after the DEM is open and before `_open_features`
(`src_python/tin_engine/cli.py@f09b9863:1116-1119`), merges the given
overrides, and passes the result to `FeatureRequest` as today. The early
range checks on given values stay where they are
(`src_python/tin_engine/cli.py@f09b9863:1000-1007`). The record's
`features_repair_m`, `features_tolerance_m` and `features_outline_snap_m`
carry the values used, so a run is reproducible from its record.

## 6. The vocabulary (B8)

`[vocabulary].extra` names edge kinds beyond the default nine. The run's
vocabulary is the default one plus every map's extras, in order of first
appearance, on bits 10 upward; bit 9 stays free because the catchment cut's
`seam` uses it (`src_python/tin_engine/features.py@f09b9863:195-197`).
Refused: an extra that repeats a default name, or more extras than bits 10 to
31 hold (22). With no extras the run's vocabulary *is*
`DEFAULT_VOCABULARY`, so its fingerprint and every file's vocabulary table are
unchanged.

The CLI builds the run's vocabulary once from the maps and passes it where it
passes `DEFAULT_VOCABULARY` today: `FeatureRequest`, `write_vtk`, the edges
`write_ply`, the GeoPackage lines writer and `start_chains`
(`src_python/tin_engine/cli.py@f09b9863:1283`, `@f09b9863:1307`,
`@f09b9863:1343`, `@f09b9863:2084`). `start_chains` already validates masks
against whatever vocabulary it is given. The gallery's roles
(`src_python/tin_engine/cli.py@f09b9863:206-219`) are untouched: the gallery
has no class maps.

## 7. Red suite for `@tester`

Fixtures: the CLI tests' `rough_dem(16)`, 17 by 21 nodes, cells 10 m by 5 m,
so `h = 10` and CORINE's distances are 1, 50, 5 m. Expected values that come
from master (the CORINE palette, the `rasputin palette corine` bytes, masks
per value) are snapshotted from `f09b9863` into `tests/python/data/` in the
red commit, before `palettes.CORINE_NATURAL` and `CLASS_MAPS` go.

**T1 shipped maps.** The shipped names are exactly `corine`, `corine-water`,
`clc18_kode`, `mapbiomas-c11`, each loaded by name; each file's `name` equals
its stem; the three CORINE class tables are equal.
**T2 CORINE behaviour unchanged.** For `corine`, `corine-water` and
`clc18_kode`, `_mask` and `_code` over the values `"111"`, `"511"`, `"523"`,
`"999"`, `"0511"`, `511` (an integer), `None`, `["311", "512"]` equal the
snapshot (water codes: bits 7 and 8, mask 384; other codes: bit 7, 128;
`corine-water`'s land codes: dropped).
**T3 palette unchanged.** `load("corine").palette` equals the snapshot,
including key order and code 0; `rasputin palette corine` writes the
snapshot's bytes; the preset title is `rasputin CORINE natural`.
**T4 file refusals**, each naming the file: missing `attribute`; a colour
not `#rrggbb`; a duplicate code; code 0 or 2^31 in `[[class]]`; an unknown
key; `unlisted` not one of the three; an edges name not in the vocabulary or
its extras; `code_system` without `positional_accuracy_m`; an accuracy of 0,
negative or not finite; a `[codes]` value outside 1 .. 2^31 - 1 or naming no
listed class; a coded class with `value`, or an uncoded class with `colour`;
a TOML syntax error (the message carries the line).
**T5 name or path.** `corine` resolves to the shipped file; `./x.toml` and
`dir/x` to paths; an unknown name is refused listing the shipped names; a
missing path is refused naming it. `rasputin palette` takes both and refuses
an uncoded map.
**T6 derived distances.** `derived_cleanup(a, h)` for `(100, 10)` and
`(100, 31)` gives `(1, 50, 5)` exactly; `(100, 200)` gives `(2, 100, 10)`;
`(30, 30)` gives `(0.3, 15, 1.5)`; `(100, 1500)` gives a snap of 50.
**T7 overrides and record.** On `rough_dem` with a coded map: no clean-up
option gives record values 1, 50, 5; `--features-tolerance 20` alone gives
1, 20, 5; a user map with `positional_accuracy_m = 40` gives 0.4, 20, 2; two
coded sources at 40 and 100 give 1, 50, 5; an uncoded-only run records 1, 50,
5.
**T8 text codes.** A GeoJSON with attribute `cover` in `"Forest"`, `"Lake"`
and a map with `[codes]` `Forest = 1`, `Lake = 2`, Lake water: triangles carry
codes 1 and 2, the lake border carries `water`; an unlisted `"Bog"` is
refused with `unlisted = "refuse"` or `"land"`, and dropped with `"drop"`.
**T9 vocabulary extras.** A map with `extra = ["levee"]` and a class whose
`edges = ["levee"]`: the `.vtk` vocabulary table has `levee` on bit 10 and
the levee lines carry it; without extras the vocabulary fingerprint equals
`DEFAULT_VOCABULARY`'s; an extra named `river` is refused; 23 extras are
refused.
**T10 MapBiomas.** The shipped map's anchor entries match the Collection 11
legend CSV (3 Forest Formation `#1f8d49`, 4 Savanna Formation `#7dc975`,
5 Mangrove `#04381d`; `@tester` reads the rest of the anchors from the CSV);
the catalogue's `mapbiomas-c11` names it; a `--land-cover` run's `.gpkg`
triangles carry class names and colours, and its `.vtk` carries no new field.
**T11 licence.** A map with `licence` writes it to the record; the CORINE
maps (empty) write no such field.
**T12 package data.** The four files are found through
`importlib.resources.files("tin_engine")`, as the NVE list is.

The existing CORINE suites (`test_cli_mesh_landcover.py`,
`test_cli_mesh_gpkg.py`, `test_cli_features_cleanup.py`,
`test_cli_mesh_multi_features.py`, `test_feature_input.py`, `test_palettes.py`,
`test_io_geopackage_write.py`) change only their imports of
`CORINE_NATURAL`, `PALETTES` and `CLASS_MAPS`, in the red commit; their
assertions stay. The MCP schema test that quotes `--features` help
(`tests/python/test_mcp_schema.py`) follows the new help text.

No mutation round: no topology or invariant-critical code; the risk is a
typo in a table, which T2 and T3 catch against master's snapshot.

## 8. Byte-identity on real data, speed, size

**Real-data check** (after green, before review; the main session or `@perf`,
about 5 minutes): on master `f09b9863` and on the branch, the quick check's
`numedalslagen` and `lagan` cases (`docs/benchmarks/quick/cases.toml`), each
written once as binary `.vtk` and compared with `cmp`; `numedalslagen` once
more as `.gpkg`, its `triangles`, `land_cover` and `layer_styles` tables
dumped with `sqlite3 .dump` and compared; and `rasputin palette corine`
compared with `cmp`. Any difference is a defect. These two cases cover both
CORINE paths in use: a pre-cut GeoPackage in the DEM's CRS at 10 m, and the
European GeoPackage in EPSG:3035 moved onto a 31 m resampled grid.

**Speed: no effect.** The clean-up runs with the same three numbers; the
added work is reading one TOML file of about 50 classes per `--features`
source, well under a millisecond. By the quick check's rule
(`docs/increments/perf-quick-check.md`, section 3) no timing run; the
byte-identity check above uses its cases, not its clock.

**Size.** By `CLAUDE.md` §2, TOML files are not counted.

| file | added | removed |
|---|---|---|
| `class_maps.py` (new: file model, loader, name or path, `derived_cleanup`, property map, vocabulary) | 125 | 0 |
| `feature_input.py` (CORINE constants and `CLASS_MAPS` out; `ClassMap` fields; `_mask`, `_code`) | 15 | 45 |
| `palettes.py` (the table and `PALETTES` out) | 0 | 52 |
| `cli.py` (map option, clean-up resolution, palette command, `--land-cover` palette, vocabulary threading, help) | 60 | 15 |
| `sources.py` (`class_map`) | 2 | 0 |
| **total** | **202** | **112** |

About **90 net lines**; with the +39 % worst case on what is added (202 x
1.39 = 281 added), about **170 net**. One PR, well under 700.

**Refine or mesh code: none.** No `@perf` acceptance run.

## 9. Defaults taken while Ola was away

- **D1** `--features-map` keeps its name and takes a name or a path.
- **D2** `s = max(a, h)`; repair `s/100`, simplification `s/2`, snap
  `min(s/20, 50)`.
- **D3** Several coded sources: the largest accuracy sets the scale.
- **D4** Runs with no coded source record the CORINE values (`a = 100`).
- **D5** MapBiomas shipped with its official legend colours (not natural
  colours, as Ola chose for CORINE in 16c), positional accuracy 30 m (its
  Landsat cell; MapBiomas states no separate positional accuracy), attribute
  `class` for polygon files.
- **D6** New edge names on bits 10 upward; bit 9 kept for the cut's `seam`.
- **D7** The three CORINE files each hold the full table; a test keeps them
  equal.
- **D8** `licence` written to the record only when not empty; CORINE's is
  empty (its notice already carries the terms), so CORINE files are unchanged.
