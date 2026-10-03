# Increment 25: plain output — what a run writes, in named fields and plain words

Status: **designed by `@architect`, 2026-10-03; revised the same day after
Ola's review of the first field table** ("This was a surprisingly long list of
things to put into a vtk file!"); **ruled by Ola the same day** ("Ruled by
Ola", at the end): both field tables approved as they stand, `--record PATH`
yes, units in the field names, `features_notice` kept in the file. Ready for
`@tester` once `@reviewer`'s round-1 findings are fixed (done, "Changes
after the design review, round 1", before "Review"; two defaults taken while
Ola was away are listed there). The mesh file carries only what a user of the mesh needs;
everything else goes to `--stats` and `--record`. Ola's ruling of
2026-10-03: before 15f-3's code step, everything `rasputin mesh` writes for a
reader is made plain. His words on the `elevation_source` sentence: "This is
tribal language". He also did not recognise "DEM holes": a raster covers its
rectangle, and what is meant is cells marked NoData. The text must say that.

Python only. No C++ changes, no change to any mesh's geometry, so no `@perf`
acceptance run is needed (`docs/increments/README.md`, "Acceptance": the diff
touches neither `include/terrain/refinement/` nor `include/terrain/mesh/`).

**Why "25".** Increments are numbered by subject, and this is not a part of
15's DEM work. 23 is the basin work and 24 is release hardening
(`docs/increments/24-release-hardening.md`, on its own branch), so 25 is the
first free number. It is placed before
15f-3 in the order of work, which is a scheduling fact, not a number.

## Prior art: legacy and literature

### Literature

There is no algorithm here; the prior art is metadata conventions for
scientific data files.

- **CF Conventions 1.11, section 2.6.2** ("Description of file contents"):
  global attributes `title`, `institution`, `source`, `history`,
  `references`, `comment`, each one named string. CF says this information is
  "mainly for the benefit of human readers", and that metadata should be
  readable by humans and parsable by programs. Units are a separate `units`
  attribute on each variable.
- **ACDD 1.3** (Attribute Convention for Data Discovery): adds `license`,
  `acknowledgement`, `processing_level` ("a textual description of the
  processing (or quality control) level of the data") and
  `geospatial_vertical_units`.

What this design takes: one named field per fact, snake_case, a short plain
value; the human-readable account (`summary`) kept separate from the facts a
script reads. What differs, and why:

- **Units go in the field name** (`tolerance_m`, `start_min_angle_deg`), not in
  a separate attribute. A VTK FieldData string array and a PLY comment have no
  slot for an attribute, and the catchment command's GeoJSON already names its
  properties this way (`outline_tolerance_m`, `fine_area_m2`, `cli.py:1730-1739`).
- **`licence_note` and `cite` keep their shipped names** (23a-2, Ola's ruling that a downloaded source's licence notes travel in the mesh file)
  rather than ACDD's `license` and `references`. Renaming them gains nothing a
  reader needs.

No novelty is claimed, so no search beyond these two conventions was made.

### Legacy

```
$ git grep -l -i -E "attrs\[|FieldData|comment |tolerance|max_error" legacy-archive -- legacy
legacy-archive:legacy/rasputin/tin_repository.py
```

`tin_repository.py` writes HDF5 attributes `timestamp`, `projection` and
`land_cover_type`: one named attribute per fact, no error or tolerance record,
no sentence. Nothing is carried across beyond that shape, which this design
also uses.

## What a run writes today (inventory, from master `586fbc1`)

Read from the code, and checked by running `rasputin mesh` on three committed
fixtures (projected: `tests/fixtures/dem_archive/7908_3_10m_z33.tif` at 5 m;
reprojected: `tests/fixtures/velhas/anadem_velhas.tif` with its catchment and
`--out-crs EPSG:31983` at 5 m; mosaic: `tests/fixtures/dtm10/seam` at 5 m;
features: the CORINE extract over the quarter-circle domain at 10 m).
Only `cli.py` writes to stderr (`grep -rn "err=True\|sys.stderr"
src_python/tin_engine` finds nothing outside it).

### The mesh file (`.vtk` FieldData; `.ply` header comments)

| field | produced at | what it is | when |
|---|---|---|---|
| `feature_bits`, `feature_names`, `feature_vocabulary` | `io/vtk_legacy.py:123-126` | the edge vocabulary: bit numbers, names, a fingerprint | always; structural, unchanged by this design |
| `crs` | `cli.py:876` | the mesh's CRS | always |
| `elevation_source` | `cli.py:876`, built at `cli.py:1496`, `:1518-1553`, `:1568-1570`, prefixed at `:870-875` | one sentence; its clauses are listed below | always (`none (z=0, --flat)` for a fixture, `:468`, `:960`) |
| `source_crs`, `source_transform`, `computation_grid` | `cli.py:880-885` | reprojected path: the DEM's own CRS, PROJ's name for the transform, the resampled grid (`square 30 m grid in EPSG:31983, node (R, K) at (30 K, -30 R), resampled bilinear from EPSG:4674`) | reprojected only |
| `licence_note`, `cite` | `cli.py:886-892` | the remote source's licence and works to cite | DEM from the cache |
| `dem_tiles`, `dem_seams` | `cli.py:893-897`, `mosaic.py:111-114` | tile names; overlaps that disagree (`a.tif \| b.tif: nodes 1, max 4, median 4`) or `none` | several files, a directory, or the cache |
| `domain`, `domain_crs`, `domain_transform` | `cli.py:898-903`, `:1461` | `catchment.geojson, 1 ring 0 holes, 27 vertices`; its CRS; transform or `none` | `--domain` |
| `features`, `features_crs`, `features_transform`, `features_notice` | `cli.py:904-934` | `clc2018_7908_3.gpkg:U2018_CLC2018_V2020_20u1, map corine, 60 features, 87 chains, 10266 vertices`; CRS; transform; the CORINE notice | `--features` |
| `land_cover_codes` | `cli.py:1045-1046`, `io/vtk_legacy.py:118-119` | what `land_cover_code` holds | a coded `--features-map` |

The `.ply` carries only some of these as comments (`cli.py:877`, `:892`,
`:900`, `:929`, `:933`, `:1009`): `crs`, `elevation` (the same sentence under
another name), `licence_note`, `cite`, `domain`, `features`,
`features_notice`, `land_cover_codes`. It lacks `source_crs`,
`source_transform`, `computation_grid`, `dem_tiles`, `dem_seams`,
`domain_crs`, `domain_transform`, `features_crs` and `features_transform`.
That is a defect this increment fixes (D2: the `.ply` carries the same fields as the `.vtk`).

### The clauses of `elevation_source`

| clause (as printed) | variable | what it counts or means | unit | non-trivial when |
|---|---|---|---|---|
| `mosaic of 4 tiles, 9 x 13 nodes;` | `cli.py:870` | the stitched grid's size, **rows x columns** (the text does not say which) | nodes | several tiles |
| `anadem-v1, <credit>;` | `cli.py:875` | the cache key and the source's credit text | | DEM from the cache |
| `refined from DEM nodes, constrained Delaunay` | `:1549` | the method: DEM nodes inserted until within tolerance, Delaunay except across lines | | `--tolerance` |
| `bilinear from DEM, stride 8` | `:1496` | no tolerance: every 8th DEM node triangulated, z interpolated | | no `--tolerance` |
| `tolerance 5 m` | `:1549` | the tolerance asked for | m | always on the refined path |
| `achieved max error 4.9997 m` | `:1550`, `out.max_error` | the largest \|DEM − mesh\| over DEM nodes in triangles whose three vertices have data; exact today (`refine.hpp:390-397`). **On the reprojected path this is phase 1's figure against the resampled grid, not against the source DEM** | m | always |
| `against the resampled grid; checked against 21615 source nodes: 37 inserted in 5 rounds, max error 4.9926 m at source nodes, 0 coincident with a start vertex (max 0 m)` | `:1518-1524`, `final_check.py:22-48` | reprojected path: the source DEM's nodes in the grown domain's rectangle (`store.size`, after duplicates); vertices the final check inserted and its rounds; the largest error at source nodes (`final.max_error`); source nodes equal to a start vertex, never checked, and their largest \|z − vertex z\| | count; m | reprojected only |
| `start stride 40` / `start domain boundary, boundary z bilinear` / `start domain boundary and features, vertex z bilinear` | `:1462`, `:1465`, `:1475` | what the starting mesh was: every 40th node, or the domain outline (and feature lines), with z interpolated at those vertices | | always |
| `start min angle 25 deg` / `start quality off` | `:1545-1547` | the start-quality setting (increment 20) | deg | always |
| `constraint feet on` / `off` | `:1551` | increment 20b: a worst DEM node very close to a line is replaced by the nearest point on the line | | always |
| `0 valid DEM nodes not covered` | `:1552`, `final.uncovered` | DEM nodes with data that lie in triangles with a NoData corner after refinement. By the stopping rule there are none (`refine.hpp:390-391`), so this is a **self-check, always 0** | count | never, by construction |
| `199 vertices without data dropped` | `:1568`, `elevation.py:25`, `:64` | mesh vertices where the DEM gives no height, removed with every triangle that uses them. With `--tolerance` a vertex at a DEM node has none when that node is NoData, and a vertex between nodes when one of the four nodes around it is. Without `--tolerance` every vertex is sampled bilinearly, and the sampler refuses a cell with any NoData corner even at zero weight (`12-dem-to-mesh.md`, R2), so a vertex on a valid node next to a NoData node is removed too: one cell of trim around NoData ("NoData on the no-tolerance path", below) | count | the DEM has NoData cells inside the area |
| `vertical unit assumed metres` | `:1569-1570` | the GeoTIFF has no `VerticalUnitsGeoKey`; any other unit than metres is refused (`io/geotiff.py:139-143`); always set for cache blocks (`fetch/run.py:249`) | | the key is missing |

### stderr

| line | at | content |
|---|---|---|
| the refine report | `cli.py:1554-1563`, `:1571` | `503 start quality nodes inserted, 252 start quality skips, 0 constraint feet, 0 feet refused, 53 rounds, 45683 points inserted, 92875 flips, 116389 triangles, achieved max error 4.9997 m, 0 valid DEM nodes not covered, 32258 start triangles, 0 start vertices off-node, 199 vertices without data dropped`. On the reprojected path its max error is phase 1's (same defect as above) |
| stride path | `:1571` | `199 vertices without data dropped` |
| mosaic | `:872` | `mosaic of 2 tiles, 256 x 563 nodes` |
| features read | `:1364-1368` | `60 features kept, 9 dropped outside, 19 clipped, 0 empty skipped` |
| no index | `:1370` | `<table>: no R-tree index, table scanned` |
| lines noded | `:1468` | `10802 input vertices, 5627 noded vertices` |
| land cover | `:1040-1044`, `landcover.py:100-108` | `land cover: 68 regions, 0 outside every polygon, 0 in more than one, 0 thinner than the snap` ("regions": groups of triangles not separated by a line; "thinner than the snap": a group whose widest triangle is narrower than twice the snap spacing, so its label may be on the wrong side) |
| `--out-crs` suggestion | `:1299` | a line to paste; plain already |
| fetch progress | `:1222` | `N of M bytes`; plain already |
| catchment command | `:1690-1729` | `window k: ... flood 0.12 s, contained`; `seed: ...`; `catchment: N nodes, X km2 of node area`; `fine outline: ... rings dropped ... holes filled ...`; `reduced outline: ...` |

### `--stats` (`stats.py`)

| section | rows | note |
|---|---|---|
| Sizes (`stats.py:186-206`) | DEM nodes (rows × columns, spacing), domain vertices, start vertices, start triangles, output vertices, output triangles, constraint edges, vertices without data dropped, file sizes | on the reprojected path "DEM nodes" is the resampled grid, not the DEM |
| DEM seams (`:286-288`) | one row per disagreeing tile pair | |
| Quality | min angle, vertex degree | unchanged by this design |
| Refinement (`:227-245`) | one wide row: tolerance, achieved max error, rounds, inserted, carved, flips, uncovered, quality inserted, quality skipped, feet | **on the reprojected path its "achieved max error" is phase 1's, and none of the final check's figures appear**; `feet refused` and `start vertices off-node` reach stderr only |
| Timings | phase rows | unchanged by this design |

The geographic run above shows it: the file's sentence says `max error
4.992618429102549 m at source nodes`, while stderr and `--stats` say
`achieved max error 4.992000034877265 m`. The two are different numbers about
different grids, and only one is about the DEM the user gave.

## The design

### D1. One record, printed in three places (four with `--record`)

A new pure module, `src_python/tin_engine/run_record.py` (no numpy, no
`_core`, no typer), holds the whole vocabulary and the record of one run:

```python
@dataclass(frozen=True, slots=True)
class Entry:
    name: str             # plain snake_case words, the unit as a suffix (_m, _deg)
    wording: str          # the label --stats prints, as in "The fields, for Ola"
    value: str            # ASCII, already formatted (numbers through stats._exact)
    number: float | int | None  # the bare number, for --record's JSON; None for text
    in_file: bool         # True only for the mesh-file fields of D2
    in_inputs: bool = False  # --stats "Inputs" (else "Result"); see "Settled after the red step"

@dataclass(frozen=True, slots=True)
class RunRecord:
    entries: tuple[Entry, ...]   # in table order; an omitted entry is absent

def file_fields(record) -> list[tuple[str, str]]      # .vtk FieldData and .ply comments
def stats_rows(record) -> list[tuple[str, str, str]]   # (wording, value, name)
def summary(record) -> str                             # one stderr line for a person
def as_json(record, version, command) -> str           # --record PATH (D5)
```

Builders (`refined_record`, `stride_record`, `flat_record`) take plain numbers
and strings that `cli.py` copies out of `RefineOutcome`, `PointRefineOutcome`,
`Trimmed`, `RasterMeta` and the input descriptions, as `stats.Refinement` is
copied today. They apply the omission rules and return a `RunRecord`.

```
_dem_mesh / mesh: outcomes, trimmed mesh, raster metadata, input descriptions
   -> numbers and strings copied out (cli.py)
   -> run_record.*_record(...)        pure; the only place wording lives
   -> RunRecord
        -> file_fields  -> write_vtk(fields=...) and write_ply(comments=...)   short subset
        -> stats_rows   -> stats.render: "Result" and "Inputs" sections        everything
        -> summary      -> stderr                                              one line
        -> as_json      -> --record PATH                                       everything, as JSON
```

Why one record: the three defects in the inventory (the reprojected run
reporting the resampled grid's error, the `.ply` missing fields, the
"DEM nodes" label on a resampled grid) are all one output drifting from
another. With one record they cannot drift, and the wording is tested without
the extension. It also shrinks `_dem_mesh`, which 15f-3 is about to grow.

### D2. What the mesh file carries

Ola's rule: only what a user of the mesh needs. The `.vtk` (FieldData) and
the `.ply` (header comments, `name value`) carry the same list.

| field | when | why it is in the file |
|---|---|---|
| `crs` | always | x and y mean nothing without it |
| `tolerance_m` | `--tolerance` | the accuracy that was asked for |
| `max_error_m` | `--tolerance` | the accuracy that was reached |
| `dem_source` | a DEM was used | where the heights came from |
| `dem_credit`, `licence_note`, `cite` | downloaded data only | the data's owner requires them to travel with the data (Ola's ruling on licence notes in the mesh file, increment 23) |
| `nodata_vertices_removed` | only when > 0 | the mesh has gaps where the DEM has none |
| `heights` | `--flat` only | the file must say that z is not real |

**`dem_source`**, one field for both cases: the DEM file names used,
separated by `; ` (one name for one file; the selected tiles for a
directory), or, for downloaded data, the dataset's name (`anadem-v1`). The
cache block names go to `--stats`.

**`max_error_m`** is defined as a bound from the start: no DEM node inside
the mesh is further than this from the mesh, in height. Today it is also the
exact largest error; after 15f-3, on the projected path, it is an upper bound
under the same name. On a reprojected run it is the error at the original
DEM's nodes (the final check's figure), which fixes the inventory's defect.

*Default taken while Ola was away:* **the nodes that sit on a vertex are
included.** Two kinds of DEM node are not compared by refinement: on the
reprojected path, an original-DEM node equal to a start vertex, whose height
can differ from the vertex's by any amount (`refine_points.hpp:28-31`, and
15c's note on these points); after 15f-3, on both paths, a node within
rounding of a vertex it is not (15f's L14 and L16). The run counts both and
gives their largest difference (`coincident`, `coincident_max_error`). So
that the definition above stays true, the record sets

`max_error_m = max(the refinement's or final check's max error, the largest difference at those nodes)`

and keeps the count and the difference itself in `--stats`
(`dem_nodes_at_vertices`, `dem_nodes_at_vertices_max_error_m`). When that
difference exceeds `tolerance_m`, stderr gets a warning naming the count and
the difference, and `max_error_m` then reads above `tolerance_m`: the file
does not claim a tolerance the mesh misses. This **answers 15f's open
question** (under its L14: whether the mesh file states the rounding
exception): the file states it through `max_error_m`, with no extra field.
At `--tolerance 0` after 15f-3, a node within rounding of a strip vertex can
make `max_error_m` a few units in the last place above 0; the tests compare
`max_error_m` with `max(tolerance_m, dem_nodes_at_vertices_max_error_m)`
(Tests, below). The safer alternative, reporting only the refinement's
figure and the exception in `--stats`, was not taken because it makes the
file's one accuracy number untrue for those nodes.

**`--flat`**: a gallery fixture has no DEM, so `heights` = `none: every z is
0 (--flat)` stays. Without it a flat mesh would look like real terrain.

**Not written any more**: `elevation_source` (and the `.ply`'s `elevation`
comment), `dem_tiles`, `dem_seams`, `source_crs`, `source_transform`,
`computation_grid`, `domain`, `domain_crs`, `domain_transform`, `features`,
`features_crs`, `features_transform`. All of them move to `--stats` (D4).
Nothing reads them from a file ("Compatibility").

#### Fields that stay because something needs them

These are not run results; they are written by the mesh writers and describe
arrays in the same file. Dropping them would make the file's own data
unreadable or break a legal requirement, so they stay, unchanged.

| field | written by | why it must stay |
|---|---|---|
| `feature_bits`, `feature_names` | `io/vtk_legacy.py:124-125`; `.ply` `feature_bit <bit> <name>` comments, `io/ply.py:111` | the key to the `feature_mask` cell array and the per-feature 0/1 arrays: without it a mask of 5 does not say "river and railway" |
| `feature_vocabulary` | `io/vtk_legacy.py:126`, `io/ply.py:112` | a digest of that key, so two files can be checked to use the same bits (increment 13) |
| `land_cover_codes` | `io/vtk_legacy.py:118-119`, `cli.py:1009` | says which code system the `land_cover_code` cell array holds (CORINE level 3); the ParaView preset from `rasputin palette corine` colours those codes and assumes that system |
| `features_notice` | `cli.py:930-933`, text at `feature_input.py:61-65` | the CORINE attribution ("Contains modified CORINE Land Cover 2018 data ... (c) European Union ..."), which the Copernicus data policy asks for on data derived from CORINE, as `licence_note` is for a downloaded DEM. Kept by Ola's ruling (2026-10-03) under the same rule as `dem_credit` |

### D3. The rules

1. **One name, one value, everywhere.** A file field appears in `--stats`
   (third column) and in `--record` with the same name and value; stderr's
   summary uses the same numbers.
2. **Units in the name**, numbers bare: `tolerance_m 5` (ruled by Ola,
   2026-10-03). Numbers are formatted by `stats._exact`, so a printed
   maximum never reads above the tolerance it met.
3. **A count of zero is omitted from the file**, never from `--stats`.
   This applies to counts only (today `nodata_vertices_removed`), never to a
   measured value: `--tolerance 0` is accepted (`cli.py:795`) and then
   `tolerance_m 0` and `max_error_m 0` are written
   (`test_cli_mesh_refine.py:133` asserts that maximum is 0.0).
4. **Self-checks** live in `--stats`. If one is ever non-zero it is also
   printed to stderr as a warning, so a broken invariant is never silent.
5. **No sentence in the file.** The one sentence, the summary, goes to stderr.
6. **Banned words** in values and on stderr: `uncovered`, `feet`, `stride`,
   `chains`, `coincident`, `carved`, `noded`, `R-tree`, `snap` (except in the
   field name `snap_to_lines`), `DEM holes`, `void`. A test enforces the list.
7. **NoData is called NoData.** Where the DEM has no height, the text says
   "NoData" or "where the DEM has no height", never "holes" or "void".

### D4. `--stats` holds the rest

`stats.render` gains two sections built from `stats_rows`:

- **Inputs**: the DEM's grid, tiles, seams, its own CRS and the resampling,
  the domain and the features (today's file fields, reworded).
- **Result**: every file field, then the counts and checks that are not in
  the file (the second table under "The fields, for Ola").

The wide one-row "Refinement" table goes; the Quality section and the Timings
stay as they are. Sizes keeps its rows (DEM nodes, domain vertices, start
vertices, start triangles, output vertices, output triangles, constraint
edges, file sizes) except "vertices without data dropped", which becomes
`nodata_vertices_removed` in Result; `start_vertices` and `start_triangles`
are printed in Sizes and not repeated in Result. On a reprojected run, "DEM
nodes" in Sizes becomes "resampled grid nodes". A table cell with `|` (a tile
name) is escaped `\|`, as today (`stats.py`, `_table`); tests compare values
after unescaping.

Renames, in `--stats` and `--record` only:

| old | new | why |
|---|---|---|
| `constraint_feet` (the clause "constraint feet on") | `snap_to_lines` | plain words. The CLI flag `--no-constraint-feet` keeps its name in this increment; renaming a flag is a separate change for Ola |
| `source_crs`, `source_transform` | `dem_crs`, `dem_transform` | pairs with `domain_crs` and `features_crs` |
| `computation_grid` | `resampled_grid` | says what it is |
| `uncovered` | `dem_nodes_outside_mesh` | says what it counts |

### D5. `--record PATH` (ruled yes by Ola, 2026-10-03)

`rasputin mesh ... --record PATH` writes the whole record of the run as a
small JSON file, for reproducibility. It is optional; without it nothing is
written. It changes nothing else about the run.

**Content.** One JSON object, keys in this order:

1. `"rasputin_version"`: `installed_version()`, as `rasputin version` prints it.
2. `"command"`: the command line, built as `--stats` builds its `command`
   (`cli.py`, `_write_report`: the program's base name, then the arguments
   exactly as given, `shlex.join`ed).
3. Every entry of the record that `--stats` prints in its Inputs and Result
   sections, in the same order, under the same names (the second table of
   "The fields, for Ola"), and so every mesh-file field of the first table
   that comes from the record, `features_notice` and `land_cover_codes`
   included. The writers' feature key (`feature_bits`, `feature_names`,
   `feature_vocabulary`) is the fixed vocabulary, not a run result, and is not
   in the record.

An entry omitted from the record (for example the 15f-3 counts on a run
without a tolerance) is absent from the JSON, not `null`. An entry that is in
`--stats` with 0 is in the JSON with 0.

**Values.** A numeric entry is a JSON number: an integer count as an integer,
a measured value (anything ending `_m` or `_deg`) as the float `Entry.number`
holds, which `json.dumps` writes as the shortest text that reads back to the
same float, so a whole number keeps its `.0` (`"tolerance_m": 5.0`). That float equals the one the `--stats` text parses to, so the two
agree to the bit. A text entry (`crs`, `start_mesh`, `snap_to_lines` as
`"on"`/`"off"`, `dem_source`, …) is a JSON string equal to its `--stats`
value.

**Determinism.** The same command in the same directory, on the same
installed version, gives the same bytes:

- no timings, no thread count, no date or time, no host or user name, no
  absolute path that the command line did not contain;
- key order fixed by the record, never by a dict or a set;
- `json.dumps(obj, indent=1, ensure_ascii=True)` plus one trailing newline;
  ASCII only, so no platform encoding enters;
- every value it holds is already deterministic: refinement's output is
  bit-identical for any thread count (increments 14 and 21), and the file
  lists, CRS labels and transform names come from the inputs.

**Example** (the projected run of the inventory, shortened):

```json
{
 "rasputin_version": "0.2.0.dev0",
 "command": "rasputin mesh --dem 7908_3_10m_z33.tif --tolerance 5 --out a.vtk --record a.json",
 "dem_grid": "5051 columns x 5051 rows, 10 m apart",
 "dem_tiles": "7908_3_10m_z33.tif",
 "dem_vertical_unit": "metres",
 "start_mesh": "every 40th DEM node",
 "start_min_angle_deg": 25.0,
 "snap_to_lines": "on",
 "crs": "EPSG:25833",
 "tolerance_m": 5.0,
 "max_error_m": 4.999692612800061,
 "dem_source": "7908_3_10m_z33.tif",
 "nodata_vertices_removed": 199,
 "refinement_rounds": 53,
 "...": "..."
}
```

**Where and when.**

- `PATH` is resolved like `--stats`'s (`_destination` with `--out-parent`).
  It is refused, before any file is written, if it resolves to the mesh file,
  the `--out-edges` file or the `--stats` file, with the same kind of usage
  error `_report_target` gives (`cli.py:1073-1085`). `-` is refused: standard
  output is `--stats -`'s.
- It is written after the mesh files and the `--stats` report, and its path is
  echoed on stdout like theirs. A refused or failed run writes no record.
- It works on every path the record covers: tolerance, no tolerance, and the
  gallery's `--flat`.

**Cost**: about 20 production lines (the option, its path checks, `as_json`).
The command line holds local paths; that is why it is in the record and not in
the mesh file.

### D6. 15f-3's figures

15f-3's counts go to `--stats` (and `--record`), not the file. The file is
affected only through `max_error_m`.

| from 15f-3 | field | where |
|---|---|---|
| strip points | `line_points_checked` | `--stats` |
| the store's NoData count | `line_points_on_nodata` | `--stats` |
| strip max error | `line_max_error_m` | `--stats` |
| refused points and their largest error | `line_points_refused`, `line_points_refused_max_error_m` | `--stats`; a stderr warning when > 0 |
| strip points inserted | `line_points_inserted` | `--stats` |
| DEM nodes inserted by the line check | `line_check_dem_nodes_inserted` | `--stats` |
| duplicate strip points | `line_points_duplicate` | `--stats` |
| points within rounding of a vertex, not checked, and their largest difference | `dem_nodes_at_vertices`, `dem_nodes_at_vertices_max_error_m` | `--stats`; **absent from the record on a path that does not produce them** (default taken while Ola was away): today only the reprojected path counts them, so before 15f-3 a projected run has neither entry, in `--stats` or `--record`, rather than a 0 nobody measured. Where a path produces them, they are always in `--stats` and `--record`, 0 included (rule 3 omits zero counts from the file only) |
| the "at most" bound | `max_error_m` | file |

**What changes in 15f-3's red suite** (branch `worktree-15f-3`, `4157dab`):

- The fourteen swaps of `achieved max error` for `max error at DEM nodes at
  most` (in `test_cli_constraint_feet.py`, `test_cli_mesh_domain.py`,
  `test_cli_mesh_domain_crs.py`, `test_cli_mesh_features.py`,
  `test_cli_mesh_refine.py`, `test_cli_start_quality.py`) are dropped: 25's red
  step rewrites those lines as `float(file_field(vtk, "max_error_m"))`, and 15f-3
  takes 25's version at the rebase. In that amendment,
  `test_cli_mesh_refine.py:133` (`--tolerance 0`, which today asserts the
  maximum is exactly 0.0) compares `max_error_m` with
  `max(tolerance_m, dem_nodes_at_vertices_max_error_m)`, read from the file
  and from `--stats`, because after 15f-3 a node within rounding of a strip
  vertex can lift it a few units in the last place above 0 (D2).
- `test_cli_mesh_edge_strip.py`:
  - the sentence-clause regex (`:106-110`) and its `clause()` reader become
    reads of the `--stats` Result rows above (or of the `--record` JSON);
    `:250-254`, `:297`, `:315` and `:396-399` read the same figures that way;
  - `:255` becomes `float(file_field(vtk, "max_error_m")) <= TOLERANCE`;
  - `:256` and `:403` (which phrase appears on which path) are dropped:
    there is no sentence; instead, `max_error_m` is in the file on both
    paths, and on the reprojected path it equals the final check's figure
    raised by the at-vertex figure (D2);
  - `:400-401` ("against the resampled grid", "checked against N source
    nodes") become `--stats` reads of `resampled_grid` and
    `dem_nodes_checked`; the order check `:402` is dropped;
  - the stderr check (`:258-261`) becomes a `--stats` check of
    `line_points_inserted` and `line_check_dem_nodes_inserted`;
  - `:90` imports `NUMBER`, `field` and `sentence` from
    `test_cli_mesh_refine.py`. **25 removes `sentence`** (there is no
    sentence) and replaces `field(text, pattern)` with a reader of one named
    field from the file and one of a `--stats` row, kept in
    `test_cli_mesh_refine.py` under the names `file_field(vtk, name)` and
    `stats_row(report, name)` so 15f-3's import changes by one line;
  - the module docstring's paragraph at `:28-37` (what the sentence pins,
    and the open question on the rounding exception) is rewritten: it pins
    fields, and the open question is answered by D2;
  - the `--stats` phase rows (`edge strip: generate`, `edge strip: scan
    (parallel)`, `edge strip: split + flip (serial)`) are unchanged.
- `test_cli_mesh_refine.py`: the docstring paragraph that 4157dab adds
  (lines 11-14 there, "Since increment 15f-3 ... the projected path says `max
  error at DEM nodes at most <e> m`") is dropped at the rebase; 25 rewrites
  that docstring to name the fields.
- The tests that call `_core` directly, and `test_edge_strip.py`, do not read
  output text and are unchanged.
- One `@tester` amendment commit on 15f-3 after its rebase onto 25, with the
  reason in the message (README, "Steps 2 and 3 are not strictly once each").

### D7. stderr

While the run goes, a person sees the input lines (reworded, table below),
then one summary line, then the output path. The long refine report goes; its
counts are in `--stats`, and `--stats -` prints them in the terminal.

Summary examples (rounded to 5 significant figures; the fields keep every
digit):

- `116389 triangles. Every DEM node inside the mesh is within 5 m of it
  (largest difference 4.9997 m). 199 vertices on NoData cells were removed
  with their triangles.`
- Reprojected: `983 triangles. Every node of the original DEM inside the mesh
  is within 5 m of it (largest difference 4.9926 m).`
- After 15f-3, if any line point could not be added: `Warning: 1 point along
  the lines could not be added; its difference is 5.4 m.`

### D8. Out of scope

- Refusal and error messages (64 `typer.BadParameter(` calls in `cli.py`):
  a separate, later pass.
- Renaming the `--no-constraint-feet` flag.
- The catchment command's GeoJSON properties: already plain, with units in
  the names.

## The fields, for Ola

Example values are from the runs in the inventory.

### In the mesh file (`.vtk` and `.ply` alike)

| field | plain wording | what it means | example | written when |
|---|---|---|---|---|
| `crs` | Coordinate system | The coordinate system of x and y | `EPSG:25833` | always |
| `tolerance_m` | Tolerance | The largest height difference allowed between mesh and DEM | `5` | with a tolerance |
| `max_error_m` | Largest height error | No DEM node inside the mesh is further than this from it (nodes sitting on a vertex included) | `4.999692612800061` | with a tolerance |
| `dem_source` | DEM | The DEM files used, or the downloaded dataset's name | `7908_3_10m_z33.tif` / `anadem-v1` | with a DEM |
| `dem_credit` | Credit | The credit the data's owner asks for | `Agencia Nacional de Aguas ... doi.org/10.5069/G9736P4G` | downloaded data |
| `licence_note` | Licence | The data's licence | `CC BY 4.0 (Creative Commons Attribution 4.0 International)` | downloaded data |
| `cite` | Please cite | Works the data's owner asks to be cited | `Laipelt, L., et al. (2024). ANADEM ...` | downloaded data that asks for it |
| `nodata_vertices_removed` | Vertices removed on NoData | Vertices where the DEM has no height (NoData cells), removed with their triangles: the mesh has gaps there | `199` | only when > 0 |
| `heights` | Heights | Says the heights are not real | `none: every z is 0 (--flat)` | `--flat` only |

Kept because the file's own arrays or a licence need them (D2): the feature
key (`feature_bits`, `feature_names`, `feature_vocabulary`),
`land_cover_codes`, and `features_notice` (the CORINE attribution).

### In `--stats` and `--record`

Every file field above, plus:

| field | plain wording | example |
|---|---|---|
| `start_mesh` | What refinement started from | `every 40th DEM node` / `the domain outline and the feature lines` |
| `start_min_angle_deg` | Starting mesh improved to this smallest angle (0 = off) | `25` |
| `snap_to_lines` | DEM nodes very close to a line were moved onto it | `on` |
| `dem_grid` | DEM grid size and spacing | `563 columns x 256 rows, 10 m apart` |
| `dem_tiles` | The tiles or downloaded blocks used | `6400_1_10m_z33.tif; 6400_4_10m_z33.tif` |
| `dem_seams` | Overlapping tiles that disagree, and by how much | `none: the tiles agree where they overlap` |
| `dem_vertical_unit` | Unit of z | `metres (assumed: the DEM file does not say)` |
| `dem_crs`, `dem_transform` | A reprojected DEM's own coordinate system, and the conversion | `EPSG:4674`, `UTM zone 23S` |
| `resampled_grid` | The grid a reprojected DEM was interpolated onto | `30 m square grid in EPSG:31983, 158 columns x 130 rows` |
| `resampled_grid_max_error_m` | Largest error against that grid, before the check against the original DEM | `4.992000034877265` |
| `dem_nodes_checked` | Nodes of the original DEM compared with the mesh | `21615` |
| `dem_check_points_inserted`, `dem_check_rounds` | Points that check added, and its passes | `37`, `5` |
| `dem_nodes_at_vertices`, `..._max_error_m` | DEM nodes on a vertex (to rounding), not compared, and their largest difference | `0`, `0` |
| `dem_nodes_outside_mesh` | Self-check: DEM nodes with data left outside the mesh (always 0) | `0` |
| `line_points_checked`, `line_max_error_m` | Points along lines (grid-line crossings and halfway between) compared with the mesh, and their largest difference (15f-3) | `1834`, `4.21` |
| `line_points_on_nodata`, `line_points_refused`, `line_points_refused_max_error_m`, `line_points_inserted`, `line_check_dem_nodes_inserted`, `line_points_duplicate` | The line check's other counts (15f-3) | `0`, ... |
| `refinement_rounds`, `points_inserted`, `points_inserted_on_nodata`, `edge_flips` | Refinement's work: passes, points added, of them into triangles with a NoData corner, edge swaps | `53`, `45683`, `7887`, `92875` |
| `start_quality_points_inserted`, `start_quality_points_skipped` | Points added to improve the starting mesh; tries skipped | `503`, `252` |
| `points_snapped_to_lines`, `snaps_refused` | Points moved onto lines; moves refused | `0`, `0` |
| `start_vertices_between_dem_nodes` | Starting-mesh vertices not on a DEM node (the starting mesh's size stays in the Sizes section; "Settled after the red step", item 6) | `0` |
| `domain`, `domain_crs`, `domain_transform` | The domain file and shape; its coordinate system and conversion | `catchment.geojson: 1 outline, 0 holes, 27 vertices` |
| `features`, `features_crs`, `features_transform` | Each features file: layer, class map, features, lines, vertices | `clc2018_7908_3.gpkg layer U2018_CLC2018_V2020_20u1, class map corine: 60 features, 87 lines, 10266 vertices` |

### stderr, reworded

| today | new |
|---|---|
| the refine report (one long line) | the summary line (D7) |
| `199 vertices without data dropped` (no tolerance) | the summary line |
| `mosaic of 2 tiles, 256 x 563 nodes` | `DEM: 2 files, 563 columns x 256 rows` |
| `60 features kept, 9 dropped outside, 19 clipped, 0 empty skipped` | `features: 60 kept (19 cut at the domain outline), 9 outside the domain, 0 empty` |
| `<table>: no R-tree index, table scanned` | `<table>: this layer has no spatial index, so every row was read` |
| `10802 input vertices, 5627 noded vertices` | `lines: 10802 vertices read, 5627 after joining shared edges and adding crossings` |
| `land cover: 68 regions, 0 outside every polygon, 0 in more than one, 0 thinner than the snap` | `land cover: 68 areas between lines; 0 in no polygon, 0 in more than one (the smallest wins), 0 too narrow to label with certainty` |
| catchment: `window 1: ..., flood 0.12 s, contained` / `grown on north, east` | `window 1: ..., searched in 0.12 s, catchment inside the window` / `window widened to the north, east` |
| catchment: `seed: the pour node at (x, y) (a pour point must lie on the flow line; it is not snapped)` | `start: the outlet node at (x, y) (an outlet must lie on the flow line; it is not moved there)` |
| catchment: `seed: lake of X km2, N nodes` | `start: lake of X km2, N DEM nodes` |
| catchment: `catchment: N nodes, X km2 of node area` | `catchment: N DEM nodes, X km2` |
| catchment: `fine outline: ..., 2 rings dropped (14 nodes), 1 holes filled (...)` | `outline along DEM cells: ..., 2 separate patches left out (14 nodes), 1 enclosed gap filled (...)` |

## Compatibility

Who reads today's text (`git grep -n elevation_source`, and the clauses'
words, over the tree):

- **Tests** (rewritten by `@tester` in 25's red step). The file-field reads
  become reads of `max_error_m`, `tolerance_m`, `dem_source` and the rest;
  reads of fields that move (`dem_tiles`, `dem_seams`, `source_crs`,
  `computation_grid`, `domain*`, `features*` but `features_notice`) become
  reads of the `--stats` Inputs section:
  `test_cli_mesh_dem.py` (`:119-272`), `test_cli_mesh_refine.py` (`:59` and
  its `field`/`sentence` helpers), `test_cli_mesh_mosaic.py` (`:140-518`, the
  `mosaic of` prefix, and `:266-271` `dem_seams`), `test_cli_mesh_geographic.py`
  (`:232-674`, `checked against`, `source_crs`, `computation_grid`),
  `test_cli_mesh_cache.py` (`:211-220`), `test_cli_fetch.py` (`:297-308`),
  `test_cli_mesh_vtk.py` (`:193`), `test_cli_mesh_domain.py`,
  `test_cli_mesh_domain_crs.py`, `test_cli_mesh_features.py`,
  `test_cli_start_quality.py`, `test_cli_constraint_feet.py`,
  `test_stats.py` (the Refinement table, `:235`, `:324-330`, `:399`),
  `test_cli_mesh.py:302` (the `.ply` comment `elevation none (z=0, --flat)`,
  which becomes `heights none: every z is 0 (--flat)`),
  `test_cli_mesh_landcover.py:13`, `:62-66` (the land-cover and features
  stderr lines, reworded), `test_cli_mesh_multi_features.py` (`features` at
  `:163-270`, `:480-516`, and `features_crs` at `:215`, `:528`: both move to
  `--stats`), `test_cli_mesh_stats.py:262`, `:285-304` (the Refinement
  section and the "vertices without data dropped" Sizes row).
  `test_io_ply.py:155` uses `elevation none (z=0, --flat)` only as a sample
  comment for the writer; no change.
  `test_io_vtk_legacy.py` and `test_io_vtk_readback.py` use `elevation_source`
  only as a sample field name for the writer; they need no change.
- **`tools/bench.py`**: does not read the text. Its numbers come from its own
  `BENCH` line, written from `refine`'s outcome (`tools/bench.py:747-752`), and
  its mesh hash covers the bytes from `POINTS` on (`:333`, `:350`), so header
  fields do not change it. Stored baselines stay comparable.
- **Dated evidence scripts** that parse today's text, each for its own
  day's logs, which they keep reading; none is changed:
  `docs/benchmarks/2026-09-26/bench1m/summarize.py:18` (`achieved max
  error`), `docs/benchmarks/2026-09-29/bygdin/analyse.py:125` and
  `docs/benchmarks/2026-09-29/bygdin-landcover/analyse.py:193` (the same),
  `docs/benchmarks/2026-10-01/basin-piece/run_sweep.py:77` and
  `docs/benchmarks/2026-10-02/basin-piece-anadem/run_sweep.py:77` (the
  `--stats` Sizes row "vertices without data dropped"). A new run of any of
  them on 25's output would need its field names.
- **`refine` stays a name in `cli.py`'s namespace**, imported from `_core`
  as today: `tools/bench.py:731-742` replaces `vars(cli)["refine"]` to count
  and time it, and 25's determinism test patches it the same way.
- **Palettes** (`palettes.py`): read nothing from the file.
- **Design records** that quote the sentence (12, 13, 14, 14b, 15, 15c, 16,
  16b, 20, 20b, 23) describe what shipped then and stay as they are. The one
  design not yet built that prescribes a clause, `15f-edge-strip.md` D7,
  points to D6 here.
- **Old files** stay readable: nothing in rasputin reads a mesh file back, and
  ParaView and QGIS show unknown FieldData and comments as text.

## Tests for `@tester`

Not invariant-critical, so no mutation round (README, "Cost constraints").

- **The record alone** (`tests/python/test_run_record.py`, no `_core`): every
  name is plain snake_case and unique; values are ASCII without control
  characters; `file_fields` holds exactly D2's fields for each path (projected,
  reprojected, downloaded, stride, `--flat`) and no other; zero
  `nodata_vertices_removed` is absent from the file and present in
  `stats_rows`; a non-zero self-check produces a stderr warning;
  `max_error_m` printed through `_exact` never reads above `tolerance_m`
  when it is not above it; every number in `summary` equals a field rounded to
  5 significant figures; no banned word (D3, rule 6) in any value or the summary.
- **One vocabulary** (CLI, a projected fixture with NoData): every file field
  is in `--stats` with the same value; the `.ply` comments equal the `.vtk`
  fields name for name, apart from the writer's own (`feature_bit` lines and
  the like); no `elevation_source`, no `elevation` comment; the kept fields of
  D2 are still written where they were.
- **The values are right**, against independent counts:
  `nodata_vertices_removed` equals the vertices whose bilinear sample is NoData,
  counted in NumPy from the DEM array, on the no-tolerance `--stride` path
  (there every start vertex is a DEM node at a known index, so the count is
  exact without re-running refinement); `max_error_m` is at least the largest
  height difference at DEM nodes inside the written triangles, computed in the
  test, and at most `max(tolerance_m, dem_nodes_at_vertices_max_error_m)`;
  `--tolerance 0` writes `tolerance_m 0` and `max_error_m 0` (rule 3 does not
  omit measured values).
- **The at-vertex rule** (D2), with the record alone: a record whose
  at-vertex difference exceeds the refinement's maximum reports that
  difference as `max_error_m`; above `tolerance_m` it also produces the stderr
  warning.
- **Reprojected** (the Velhas fixture): `max_error_m` is the error at the
  original DEM's nodes (checked independently, as 15c-2's tests do), not the
  resampled grid's; `--stats` has `dem_crs`, `resampled_grid`,
  `dem_nodes_checked` and "resampled grid nodes".
- **The other paths**: no tolerance (no tolerance fields), `--flat`
  (`crs` when given, `heights`), tiles (`dem_source` lists the files;
  `dem_seams` in `--stats` both ways), downloaded (`dem_source`, `dem_credit`,
  `licence_note`, `cite`), features (`features_notice` and `land_cover_codes`
  in the file; `features` in `--stats`).
- **stderr**: each reworded line by `re.search` on its parts, and the banned
  words over the whole of stderr of successful runs (a refusal's text quotes
  the user's own input and is out of scope, D8).
- **`command`**: tests that read it (in `--stats` or `--record`) set
  `sys.argv` with `monkeypatch`, because `CliRunner` does not, and the
  `command` is built from `sys.argv` (`cli.py:1118`).
- **`--stats` values with `|`** are compared after unescaping `\|`.
- **`--record`** (D5):
  - the file parses as JSON, is ASCII, and ends in one newline;
  - its keys are `rasputin_version`, `command`, then exactly the `--stats`
    Inputs and Result names in the `--stats` order;
  - counts are JSON integers; every `_m` and `_deg` value is a number equal to
    `float` of its `--stats` text; text values equal their `--stats` text;
  - every mesh-file field from the record is in it with the same value;
  - determinism: one command run twice gives the same bytes, and so does a
    run with `tin_engine.cli.refine` monkeypatched to a wrapper that passes
    `threads=1` (the binding's keyword, `_core.pyi:437`; the CLI has no thread
    option); no key or value contains a time, a date, the host name, or the
    word `seconds`;
  - an omitted entry is absent, not `null`;
  - refusals, each before any file exists: `PATH` equal to the mesh, to
    `--out-edges`, to the `--stats` file, and `-`;
  - a refused run (for example a bad `--tolerance`) leaves no record file;
  - `--flat` and the no-tolerance path write a record too.

## LOC estimate

`CLAUDE.md` §2's unit; tests excluded.

| file | what | net est. |
|---|---|---:|
| `run_record.py` (new) | `Entry`, `RunRecord`; about 45 entries with their wording; three builders and the omission and at-vertex rules; `file_fields`, `stats_rows`, `summary`, `as_json` | 170 |
| `cli.py` | the sentence, the refine report and most of the field assembly replaced by copying numbers into the builder and one short field list for both formats; reworded stderr lines and warnings; catchment wording; `--record` (the option, its path checks, the write) | 0 |
| `stats.py` | the Inputs and Result sections in place of the Refinement table; the Sizes rows | 30 |
| **total** | | **200** (+39 %: 278; +60 %: 320) |

Revised after `@reviewer`'s round 1, which expected 200 to 280: the first
estimate (135) counted the vocabulary at a line per entry and missed the
builders' plumbing. `cli.py` is put at zero net: it loses about 90 lines (the
sentence, the report, the two parallel field and comment lists) and gains
about as many (copying out the outcome numbers, the warnings, `--record`).
One PR. 15f-3 then adds its counts as builder arguments, about 10
lines fewer than its own estimate for the sentence clause.

## Ruled by Ola

Ola, 2026-10-03, on the revised design:

- **Both field tables are approved as they stand**: the mesh-file list,
  including the four fields kept with their reasons (the feature key,
  `land_cover_codes`, `features_notice`), and the `--stats` list.
- **`--record PATH`: yes.** Specified in D5.
- **Units in the field names: yes** (`tolerance_m 5`, a bare number).
- **`features_notice` stays in the mesh file**: the CORINE credit travels with
  the mesh, as `dem_credit` does for a downloaded DEM.

No question is open.

## Changes after the design review, round 1

`@reviewer`'s five blocking findings, each fixed above:

1. Rule 3 (D3) now omits zero **counts** only; a measured value of 0 is
   written (`--tolerance 0`).
2. D5's example shows measured values as floats (`5.0`, `25.0`), as
   `json.dumps` writes them.
3. Compatibility lists the test files the first grep missed, the dated
   evidence scripts, and the `refine` name that `bench.py` patches.
4. D6 lists every 15f-3 red-test line that reads the sentence, says that 25
   removes the `sentence` helper and what replaces it, and covers the two
   docstrings.
5. `max_error_m` includes the nodes that sit on a vertex (D2), which answers
   15f's open question on stating the rounding exception in the file.

The suggestions taken: the wrong section reference; internal labels (a
ruling's code, a test group's code) replaced by what they name; the NoData
count test on the no-tolerance path; `|` escaping in `--stats`; the banned
words checked on successful runs only; the Sizes rows that stay; the
catchment's lake line; `sys.argv` for `command`; the LOC estimate.

**Defaults taken while Ola was away** (for Ola to confirm or overrule):

- **`max_error_m` includes the at-vertex nodes** (D2). The alternative was
  to keep it as the refinement's figure and state the exception only in
  `--stats`; it was not taken because the file's one accuracy number would
  then be untrue for those nodes.
- **`--record -` is refused** (D5), so standard output stays `--stats -`'s.
- **`dem_nodes_at_vertices` and its largest difference are absent from the
  record on a path that does not measure them** (D6): before 15f-3, the
  projected path. A 0 there would claim a check that did not run.

## Settled after the red step (73863c3)

`@tester`'s red step pinned choices this design had not made. Each is
confirmed or corrected here before `@developer` starts. Two are corrected
(items 7 and 8), and `@tester` changes the tests named there; the rest stand
as pinned.

1. **Builders: confirmed.** `refined_record`, `stride_record` and
   `flat_record` take keyword arguments only, named after the record's
   entries, plus `triangles` (the output triangle count, for the summary). An
   argument not passed means the path does not produce that entry, and the
   entry is absent. `refined_record`'s `max_error_m` argument is the measured
   figure (refinement's, or the final check's on a reprojected run); the
   builder applies the on-vertex rule of D2 and may raise it.
2. **Warnings: confirmed.** `summary(record)` returns the summary sentence,
   then one line per warning, each starting `Warning:`. `cli.py` prints it to
   stderr as it is.
3. **`--stats` layout: confirmed.** A section `## Inputs` and a section
   `## Result`, each one table with rows `| wording | value | name |`; the
   name may be in backticks; `|` in a cell is escaped `\|`. A section with no
   rows is not printed (the gallery run has no Inputs). The order of sections,
   which the tests leave open, is: Sizes, Inputs, DEM seams, Quality (plan
   view, x/y), Result, Timings.
   Which entries are Inputs: `dem_grid`, `dem_tiles`, `dem_seams`,
   `dem_vertical_unit`, `dem_crs`, `dem_transform`, `resampled_grid`,
   `domain`, `domain_crs`, `domain_transform`, `features`, `features_crs`,
   `features_transform`, and the settings `start_mesh`,
   `start_min_angle_deg`, `snap_to_lines`. Everything else is Result, which
   starts with the file fields in D2's order. `Entry` gains `in_inputs: bool =
   False` (defaulted, so the five-argument `Entry(...)` in the tests holds);
   `stats_rows(record)` stays one list in record order, Inputs entries first,
   and the builders emit them first.
4. **`stats.Report`: confirmed.** Its `refinement` argument and the
   `Refinement` dataclass go, and so does `Sizes.dropped` with its row. The
   new sections reach it as two fields, `inputs` and `result`, each a
   sequence of `(wording, value, name)`, both defaulting to empty; `cli.py`
   splits `stats_rows` by `in_inputs` to fill them.
5. **`--record` key order: the test plan wins.** Keys are
   `rasputin_version`, `command`, then the record in `--stats` order, Inputs
   then Result. D5's example was wrong (it put `crs` first and `dem_grid`
   after `snap_to_lines`) and is corrected above.
6. **`start_vertices` and `start_triangles`: confirmed**, in Sizes only, so
   absent from Result and from `--record`. The `--stats` table under "The
   fields, for Ola" listed them; `start_vertices_between_dem_nodes` is the
   only one of the three in Result.
7. **`start_mesh`: corrected wording, same meaning.** It is present on the
   no-tolerance path too (that mesh is a start mesh with no refinement), and
   for a domain it reads `the domain outline`, with features `the domain
   outline and the feature lines`. For a stride `n` it reads `every DEM node`
   when `n` is 1 and otherwise `every <n><ordinal> DEM node` with the English
   ordinal: `every 2nd`, `every 3rd`, `every 21st`, `every 40th`, `every
   112th`. "every 2th" would not be plain. `recordread.py`'s parser
   (`every (?:(\d+)(?:st|nd|rd|th) )?DEM node`) already accepts this; any
   test that builds the expected text as `f"every {n}th DEM node"` must use
   the ordinal instead (`@tester`).
8. **Wording.** Confirmed: `dem_seams` with no disagreement reads `none: the
   tiles agree where they overlap`; `domain` reads `<file>: 1 outline, <h>
   hole(s), <n> vertices`; `features` reads `<file>[ layer <l>], class map
   <m>: <n> features, <c> lines, <v> vertices`, sources joined by `; `;
   `dem_grid` starts `<c> columns x <r> rows, `; `dem_vertical_unit`
   contains `metres`, and `assumed` only when the GeoTIFF has no vertical
   unit. **Corrected:**
   - **A disagreeing seam** is not today's `Seam.entry` text (`a.tif | b.tif:
     nodes 1, max 4, median 4`), which is the terse form Ola ruled out. It
     reads `<a> and <b> disagree at <n> node(s), by up to <max> m (median
     <median> m)`, for example `ne.tif and nw.tif disagree at 1 node, by up to
     4 m (median 4 m)`, pairs joined by `; ` in today's order. Numbers are
     formatted as `Seam.cells` formats them (`:g`). `@tester` changes the
     tests that pin the old text (`test_cli_mesh_mosaic.py:278`, `:302`,
     `:331-336`). The separate "## DEM seams" table (item 10) keeps its
     columns.
   - **Plurals are correct English** everywhere: `1 hole`, `2 holes`,
     `1 feature`, `1 line`, `1 vertex`; in the catchment, `1 separate patch`,
     `2 separate patches`, `1 enclosed gap`, `2 enclosed gaps`. The catchment
     and domain regexes already accept both forms; the `features` regex
     (`(\d+) features, (\d+) lines, (\d+) vertices`) must accept `feature`,
     `line` and `vertex` for a count of 1 (`@tester`).
9. **`--record` refusals: confirmed.** Each is a usage error (exit 2) naming
   `--record`, raised before any file is written; the record's path is the
   last line on stdout.
10. **"## DEM seams": confirmed**, kept as its own section.
11. **`tests/python/recordread.py`: confirmed**, holding `file_field` and
    `stats_row`; `test_cli_mesh_refine.py` re-exports them, so 15f-3's
    import from it keeps working (D6).
12. **The stride golden hashes the file from `POINTS` on: confirmed.** The
    header now holds fields whose wording this increment changes, so a
    whole-file digest would pin wording instead of the mesh. The bytes from
    `POINTS` on are the mesh and its arrays, which this increment must not
    change; the fields are pinned by name in the other tests. This is the
    same span `tools/bench.py` hashes (`:333`, `:350`).

## NoData on the no-tolerance path (found at the green step, a2c635d)

**What was found.** `test_cli_mesh_plain_output.py::TestTheValuesAreRight::
test_nodata_vertices_removed_counts_the_strided_nodes_without_data[1]`
expects 25 and gets 36. On the path without `--tolerance` the vertices are
sampled by `_core.sample`, which refuses a cell with any NoData corner, even a
corner of zero weight. So a vertex on a valid DEM node next to a NoData node
gets no height and is removed. `@developer` probed it: invalid rows 5-10 x
columns 7-12 against a NoData block at rows 6-10 x columns 8-12. This is
increment 12's one-cell trim, recorded there as "accepted for now ... changing
the sampler is not this increment's job" (`12-dem-to-mesh.md`, R2). Increment
14 notes that it does not apply with `--tolerance`, where z is read at nodes.
The inventory row above described only the tolerance path; it is corrected.

**Recommended: treat today's behaviour as a defect, to be fixed later in its
own small C++ increment, and pin 25's test to today's behaviour meanwhile.**
The sampler drops valid data: a DEM node with a value gets no height because a
neighbour has none. Increment 12 accepted this as a stopgap and did not argue
it was right. The fix (skip a corner whose weight is exactly zero, so a vertex
on a node reads that node) is C++ in `raster/sample.hpp` with its own tests.
Its acceptance has to show that the tolerance path is unchanged, and that
does not belong in a Python-only output increment. Describing it as the
intended rule instead would put a known loss of data into the field's
definition.

Until that fix:

- **The field keeps its name and file wording.** Its meaning on the
  no-tolerance path is "vertices on or next to a NoData cell". The inventory
  row above says so.
- **The summary on the no-tolerance path says it plainly.** `stride_record`'s
  summary reads `... vertices on or next to NoData cells were removed with
  their triangles.` instead of `on NoData cells`. With `--tolerance` the
  wording is unchanged. That is a one-line change for `@developer`.
- **What `@tester` changes** (one amendment commit, the reason in the
  message). In `test_nodata_vertices_removed_counts_the_strided_nodes_without_data`,
  the expected count is no longer "picked nodes that are NoData". It is the
  picked nodes `(r, c)` whose bilinear cell has a NoData corner, computed in
  NumPy with the sampler's rule. The cell's top-left node is
  `(min(r, rows - 2), min(c, cols - 2))`; the corners are that node and the
  nodes one row down, one column right, and one of each. The vertex count
  assertion uses the same expected count. The docstring says the count pins
  increment 12's one-cell trim, and that the sampler follow-up will change it
  back to the NoData nodes alone. `step` 1, 2 and 3 stay.
- **A follow-up row in `ROADMAP.md`**: the sampler reads a vertex on a node
  from that node. It is proposed, not designed.

## As built (green, a2c635d)

Where `@developer`'s code departs from, or fills in, the text above. Each is
accepted as described unless marked otherwise.

- **No `dem_grid` on a reprojected run.** There `resampled_grid` carries the
  grid's size and spacing, and the original DEM has no single grid in the
  target CRS. Accepted: one grid, named for what it is.
- **`land_cover_codes` is a record entry with `in_file=False`.** It appears
  in `--stats` and `--record`, and the writers keep putting it into the file
  as before (`write_vtk(land_cover_codes=...)`, the `.ply` comment). The file
  is unchanged; D2's kept-fields table holds.
- **`RunRecord.triangles`**, the output triangle count that the summary
  names, is a field of the record, not an entry. It is not in `--stats` or
  `--record` as an entry; Sizes already prints "output triangles".
- **`Sizes.resampled`** (a bool) switches the Sizes label to "resampled grid
  nodes" (D4).
- **Summary above the tolerance:** `The largest difference at a <DEM node |
  node of the original DEM> inside the mesh is <e> m, above the tolerance of
  <t> m.`, in place of the "within" sentence. Accepted.
- **The two warnings:** `Warning: <n> DEM nodes on a vertex differ from it by
  up to <d> m, more than the tolerance of <t> m.` and `Warning: <n> DEM
  nodes with data lie outside the mesh; there should be none.` Accepted, with
  one correction: both must use the singular for a count of 1 (`1 DEM node
  ... differs`, `1 DEM node with data lies`), as "Settled after the red step",
  item 8 rules. They do not today. This is a production fix for
  `@developer`; a test for it is `@tester`'s.
- **`--flat` summary:** `<n> triangles. The heights are not real: every z
  is 0 (--flat).` Accepted.
- **Singular counts** go through `run_record.plural(n, one, many)`. Accepted;
  the two warnings above should use it too.
- **`mosaic.py`'s `Seam.entry`** returns the settled wording (`<a> and <b>
  disagree at <n> node(s), by up to <max> m (median <median> m)`), so the
  wording lives with the seam rather than in `run_record.py`. Accepted:
  `Seam.cells` sits beside it and formats the same numbers.

## Review

**Design review, round 1, 2026-10-03.** Range `origin/master` (6cdc8cc) `..97d1e9c` (8de9714, 2462028, 16b3a5d, 97d1e9c). Verdict: CHANGES REQUESTED. LOC: 0 (design only); the reviewer thinks 135 is optimistic and expects 200 to 280, still far under 700. Inventory, the three bugs, bench.py, the 15f-3 red commit's fourteen swaps and the --record spec's buildability all check out. Blocking: (1) rule 3 (zeros omitted from the file) contradicts D2 at `--tolerance 0`; restrict it to counts; (2) D5's example JSON writes measured values as integers (`5`), against its own float rule; (3) Compatibility misses four test files (`test_cli_mesh.py`, `test_cli_mesh_landcover.py`, `test_cli_mesh_multi_features.py`, `test_cli_mesh_stats.py`); (4) D6 is incomplete on how 15f-3's red tests change (more lines in `test_cli_mesh_edge_strip.py`, its `sentence` import, docstrings, the extra paragraph in `test_cli_mesh_refine.py`); (5) `max_error_m`'s definition is not true for DEM nodes on or within rounding of a vertex; take the larger figure or state the exception, warn above tolerance, and say what happens to 15f's open question on stating it. Not pushed; no CI.

**Design review, round 2, 2026-10-03.** Range `333eb7d..ef52e8e` (231eefa merge of origin/master, ef52e8e fixes). Verdict: APPROVED. LOC: 0 (design only); estimate now about 200, plausible. All five round-1 blockers closed; the two defaults taken while Ola was away (`max_error_m` includes DEM nodes on a vertex, with a stderr warning above tolerance; `--record -` refused) are coherent and testable. Suggestions: say how 15f-3 amends the `--tolerance 0` test at `test_cli_mesh_refine.py:133`; use `file_field` consistently in D6; say whether the on-vertex entries are omitted or 0 on the projected path before 15f-3. Not pushed; no CI.
