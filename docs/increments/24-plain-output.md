# Increment 24: plain output — what a run writes, in named fields and plain words

Status: **designed by `@architect`, 2026-10-03, for Ola's approval of the
field table ("The fields, for Ola") before `@tester` starts.** Ola's ruling of
2026-10-03: before 15f-3's code step, everything `rasputin mesh` writes for a
reader is made plain. His words on the `elevation_source` sentence: "This is
tribal language". He also did not recognise "DEM holes": a raster covers its
rectangle, and what is meant is cells marked NoData. The text must say that.

Python only. No C++ changes, no change to any mesh's geometry, so no `@perf`
acceptance run is needed (`docs/increments/README.md`, "Acceptance": the diff
touches neither `include/terrain/refinement/` nor `include/terrain/mesh/`).

**Why "24".** Increments are numbered by subject; 23 is the last number in use
(`ROADMAP.md`), and this is not a part of 15's DEM work. It is placed before
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
- **`licence_note` and `cite` keep their shipped names** (23a-2, Ola's B16)
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
That is a defect this increment fixes (D3).

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
| `199 vertices without data dropped` | `:1568`, `elevation.py:25`, `:64` | mesh vertices where the DEM has no height (the node is NoData, or one of the four nodes around an off-node vertex is), removed with every triangle that uses them | count | the DEM has NoData cells inside the area |
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

### D1. One record, three renderings

A new pure module, `src_python/tin_engine/run_record.py` (numpy-free, no
`_core`, no typer), holds the vocabulary and the record:

```python
@dataclass(frozen=True, slots=True)
class Entry:
    name: str        # snake_case, ^[a-z][a-z0-9_]*$, the unit as a suffix (_m, _deg)
    wording: str     # the plain label, as in "The fields, for Ola"
    value: str       # ASCII, already formatted (numbers through stats._exact)
    in_file: bool    # False: --stats only (self-checks, process counters)

@dataclass(frozen=True, slots=True)
class RunRecord:
    entries: tuple[Entry, ...]   # in table order; omitted entries are absent

def file_fields(record) -> list[tuple[str, str]]   # .vtk FieldData and .ply comments
def summary(record) -> str                          # the `summary` field and stderr's line
def stats_rows(record) -> list[tuple[str, str, str]]  # (wording, value, name) for --stats
```

A builder, `refined_record(...)` / `stride_record(...)` / `flat_record(...)`,
takes plain numbers and strings (copied out of `RefineOutcome` and
`PointRefineOutcome` by `cli.py`, as `stats.Refinement` is today), applies the
omission rules of D2, and returns a `RunRecord`. `cli.py` then:

- writes `file_fields(record)` plus the input-description fields (`domain`,
  `features`, …) into the `.vtk` and, as `name value` comments, the `.ply`;
- prints `summary(record)` to stderr in place of the refine report;
- passes the record to `stats.render`, which prints `stats_rows` in a
  "Result" section (item | value | field) in place of the wide Refinement row.

Data flow:

```
_dem_mesh: RefineOutcome / PointRefineOutcome / Trimmed / RasterMeta
   -> numbers copied out (cli.py)
   -> run_record.*_record(...)          pure, no _core: the only place wording lives
   -> RunRecord
        -> file_fields -> write_vtk(fields=...) / write_ply(comments=...)
        -> summary     -> stderr, and the `summary` field
        -> stats_rows  -> stats.render (Result section)
```

Why a module and not more f-strings in `_dem_mesh`: the vocabulary is then in
one place, the three outputs cannot drift apart (D5's test reads all three
from one run), and the wording is unit-tested without the extension. It also
shrinks `_dem_mesh`, which 15f-3 is about to grow.

### D2. The rules

1. **One name, one value, everywhere.** A field written to the file has the
   same name and value in `--stats` (the third column) and, where stderr
   mentions it, the same number in `summary`.
2. **Units in the name**, numbers bare: `tolerance_m 5`. Numbers are
   formatted by `stats._exact`, so a printed maximum never reads above the
   tolerance it met.
3. **Omit what is trivially zero** from the file: `nodata_vertices_removed`,
   `dem_nodes_at_vertices*`, `line_points_refused*`,
   `line_points_on_nodata`. `--stats` always shows them, with 0.
4. **Self-checks and process counters go to `--stats` only**: DEM nodes left
   outside the mesh (today's `uncovered`), rounds, flips, insertions, start
   quality, feet counts, start sizes, the final check's insertions, phase 1's
   error against the resampled grid, the line check's insertions. If a
   self-check is ever non-zero, it is written to the file too and a line is
   printed to stderr, so a broken invariant is never silent.
5. **No sentence a script must parse.** The one sentence left, `summary`, is
   built from the fields and adds no fact of its own.
6. **Words that are banned from values and from stderr**: `uncovered`,
   `feet` (the field name `constraint_feet` matches the CLI flag; its value
   is `on` or `off`), `stride`, `chains`, `coincident`, `carved`, `noded`,
   `R-tree`, `snap`, `DEM holes`, `void`. D5 tests this list.
7. **NoData is called NoData.** Where the DEM has no height, the text says
   "NoData" and "where the DEM has no height", never "holes" or "void".
8. **`elevation_source` is not written any more**, and the `.ply`'s
   `elevation` comment goes with it. Nothing reads them (see
   "Compatibility").

### D3. The `.ply` carries every field

`.ply` comments become `name value` for every field the `.vtk` carries, in
the same order (one list, `comments = [f"{k} {v}" for k, v in fields]`). This
fixes the nine missing fields of the inventory and removes the parallel
`comments` list in `cli.py`.

### D4. Renames, and why each

| old | new | why |
|---|---|---|
| `elevation_source` | the fields `heights`, `tolerance_m`, `max_error_m`, … and `summary` | Ola's ruling |
| `source_crs`, `source_transform` | `dem_crs`, `dem_transform` | "source" is our word; this pairs with `domain_crs`/`domain_transform` and `features_crs`/`features_transform` |
| `computation_grid` | `resampled_grid` | says what it is; value rewritten without `node (R, K) at (30 K, -30 R)` |
| `--stats` "DEM nodes" on the reprojected path | "resampled grid nodes" | it is not the DEM |

`dem_tiles`, `dem_seams`, `licence_note`, `cite`, `domain*`, `features*`,
`land_cover_codes` keep their names; some values are reworded (table below).
`dem_tiles` is now written for a single file too: a mesh travels without its
command line.

### D5. The 15f-3 fields, from the start

15f-3 adds the edge strip. Its figures land as fields, not as D7's clause in
`15f-edge-strip.md` (that section is marked superseded by this one):

| `PointRefineOutcome` / store | field | where |
|---|---|---|
| `strip_points` (`refine_points.hpp:90`) | `line_points_checked` | file |
| `ConstraintCheckPoints::no_data` (`constraint_points.hpp:70`) | `line_points_on_nodata` | file when > 0; `--stats` |
| `strip_max_error` | `line_max_error_m` | file |
| `strip_refused`, `strip_refused_max_error` | `line_points_refused`, `line_points_refused_max_error_m` | file when > 0; `--stats` |
| `strip_inserted` | "line points inserted" | `--stats` only |
| `nodes_inserted` (projected path) | "DEM nodes inserted by the line check" | `--stats` only |
| `ConstraintCheckPoints::duplicates` | "line points that coincide" | `--stats` only |
| `coincident`, `coincident_max_error` (L14, L16: points within rounding, `r(g)`, of a vertex they are not; both paths after 15f-3, the reprojected path today) | `dem_nodes_at_vertices`, `dem_nodes_at_vertices_max_error_m` | file when > 0; `--stats` |

**"Max error at DEM nodes at most …" (15f D7).** The field is `max_error_m`,
and its meaning is defined from the start as a bound: "no DEM node inside the
mesh is further than this from it". Today it is also the exact maximum; after
15f-3, on the projected path, it is `max(out.max_error, final.max_error)`, an
upper bound, with the same name and the same wording. On the reprojected
path it is the final check's figure at the source DEM's nodes, today and
after (this fixes the inventory's defect).

**What changes in 15f-3's red suite** (branch `worktree-15f-3`, `4157dab`):

- The fourteen swaps of `achieved max error` for `max error at DEM nodes
  at most` (`test_cli_constraint_feet.py`, `test_cli_mesh_domain.py`,
  `test_cli_mesh_domain_crs.py`, `test_cli_mesh_features.py`,
  `test_cli_mesh_refine.py`, `test_cli_start_quality.py`) are dropped: 24's
  red step already rewrites those lines as `float(field(vtk, "max_error_m"))`,
  so 15f-3 takes 24's version of each at the rebase.
- `test_cli_mesh_edge_strip.py`: the `CLAUSE` regex (`:107`) becomes field
  reads of D5's names; the ordering assert (`:402`, the final check before the
  strip clause) is dropped, since fields have no sentence order; `:256` and
  `:403` (which phrase is where) become: `max_error_m` is present on both
  paths, and on the reprojected path equals the final check's error; the stderr
  `REPORT` regex (`:258-260`) becomes a `--stats` check on "line points
  inserted" and "DEM nodes inserted by the line check"; PY5's phase rows are
  unchanged.
- `test_core_*` and `test_edge_strip.py` do not read the output text and are
  unchanged.
- This lands as one `@tester` amendment commit on 15f-3 after it is rebased
  onto 24, with the reason in the message (README, "Steps 2 and 3 are not
  strictly once each").

### D6. stderr

After the run, stderr has the `summary` sentence and the input lines, reworded
(table below). The long refine report goes; its counters are in `--stats`
(the "Result" section). `rasputin mesh ... --stats -` prints them to stdout
for anyone who wants them in a terminal.

### D7. Out of scope

- Refusal and error messages (64 `typer.BadParameter(` calls in `cli.py`). Many are plain
  already; a pass over them is a separate, later item.
- `feature_bits`, `feature_names`, `feature_vocabulary`: structural, read by
  ParaView users through the arrays they name.
- The catchment command's GeoJSON properties: already snake_case with units
  in the names.

## The fields, for Ola

Example values are from the runs listed in the inventory. "File" means the
`.vtk` FieldData and the `.ply` header; every file field is also in
`--stats`.

### The mesh file

| field | plain wording | what it means | example | written when |
|---|---|---|---|---|
| `summary` | Summary | One sentence for a person, made from the fields below | see "The summary" | always |
| `crs` | Coordinate system | The coordinate system of x and y | `EPSG:25833` | always (unchanged) |
| `heights` | Heights | How each vertex got its height | `DEM value at DEM nodes; interpolated (bilinear) from the four nearest DEM nodes at other vertices` | always |
| `tolerance_m` | Tolerance | The largest height difference allowed between the mesh and the DEM | `5` | `--tolerance` |
| `max_error_m` | Largest height error at DEM nodes | No DEM node inside the mesh is further than this from the mesh, in height | `4.999692612800061` | `--tolerance` |
| `nodata_vertices_removed` | Vertices on NoData, removed | Mesh vertices where the DEM has no height (a NoData cell), removed together with their triangles | `199` | > 0 |
| `start_mesh` | Starting mesh | What refinement started from | `every 40th DEM node` / `the domain outline` / `the domain outline and the feature lines` | `--tolerance` |
| `start_min_angle_deg` | Smallest angle of the starting mesh | The starting mesh was improved until no triangle had a smaller angle; 0 means not done | `25` | `--tolerance` |
| `constraint_feet` | Points near lines moved onto them | `on`: a DEM node very close to a line was replaced by the nearest point on the line, to avoid thin triangles | `on` | `--tolerance` |
| `dem_tiles` | DEM files | The DEM files (or cache blocks) the heights came from | `7908_3_10m_z33.tif` | always with a DEM (now also for one file) |
| `dem_grid` | DEM grid | Size and spacing of the DEM grid used | `563 columns x 256 rows, 10 m apart` | projected DEM |
| `dem_vertical_unit` | Height unit | The unit of z | `metres` / `metres (assumed: the DEM file does not say)` | always with a DEM |
| `dem_seams` | Where DEM tiles disagree | Pairs of overlapping tiles that give different heights, and by how much; the overlap is split down the middle | `ne.tif and nw.tif: 1 node, up to 4 m apart (median 4 m)` / `none: the tiles agree where they overlap` | several tiles |
| `dem_source` | DEM source | The catalogue name of a downloaded DEM | `anadem-v1` | DEM from the cache |
| `dem_credit` | DEM credit | The credit the data's owner asks for | `Agência Nacional de Águas ... https://doi.org/10.5069/G9736P4G.` (ASCII-escaped) | DEM from the cache |
| `licence_note` | Licence | The data's licence | `CC BY 4.0 (Creative Commons Attribution 4.0 International)` | DEM from the cache (unchanged) |
| `cite` | Please cite | Works the data's owner asks to be cited | `Laipelt, L., et al. (2024). ANADEM ...` | DEM from the cache (unchanged) |
| `dem_crs` | DEM's own coordinate system | The DEM's coordinate system before it was reprojected | `EPSG:4674` | reprojected DEM |
| `dem_transform` | How the DEM was reprojected | PROJ's name for the conversion | `UTM zone 23S` | reprojected DEM |
| `resampled_grid` | Grid the DEM was resampled to | The square grid the reprojected DEM was interpolated onto | `30 m square grid in EPSG:31983, 158 columns x 130 rows, lines at whole multiples of 30 m; heights interpolated (bilinear) from the DEM` | reprojected DEM |
| `dem_nodes_checked` | DEM nodes checked | Nodes of the original DEM, in and around the area, compared with the mesh | `21615` | reprojected DEM |
| `dem_nodes_at_vertices` | DEM nodes on a vertex, not checked | DEM nodes that sit on a mesh vertex (to rounding) and so are not compared | `2` | > 0 |
| `dem_nodes_at_vertices_max_error_m` | Their largest height difference | The largest difference between such a node and the vertex it sits on | `1e-09` | > 0 |
| `line_points_checked` | Points checked along lines | Points on the domain outline and feature lines where they cross a DEM grid line, and halfway between crossings, each compared with the mesh (15f-3) | `1834` | `--tolerance` (after 15f-3) |
| `line_max_error_m` | Largest height error along lines | The largest difference at those points | `4.21` | as above |
| `line_points_on_nodata` | Line points on NoData | Points along lines where the DEM has no height, so not checked | `3` | > 0 |
| `line_points_refused` | Line points that could not be added | Points over the tolerance that could not be inserted without folding a triangle | `1` | > 0 |
| `line_points_refused_max_error_m` | Their largest height error | The largest difference at those points | `5.4` | > 0 |
| `domain` | Domain | The file and shape of the area meshed | `catchment.geojson: 1 outline, 0 holes, 27 vertices` | `--domain` |
| `domain_crs`, `domain_transform` | Domain's coordinate system; how it was converted | as today | `EPSG:31983`; `none` | `--domain` (unchanged) |
| `features` | Features | Each features file: layer, class map, and how many features, lines and vertices went in | `clc2018_7908_3.gpkg layer U2018_CLC2018_V2020_20u1, class map corine: 60 features, 87 lines, 10266 vertices` | `--features` |
| `features_crs`, `features_transform`, `features_notice` | as today | as today | | `--features` (unchanged) |
| `land_cover_codes` | Land cover codes | What the per-triangle `land_cover_code` numbers mean | `CORINE Land Cover level-3 code, attribute Code_18, class map corine; 0 = no polygon, and on every line` | coded class map |

On the gallery path (`--flat`), the file has `heights` = `none: every z is 0
(--flat)` and `summary`. On the stride path (no `--tolerance`), `heights` is
`interpolated (bilinear) from the DEM at every vertex; vertices are every 8th
DEM node`, with no tolerance fields.

### The summary

Made from the fields above and nothing else. Examples:

- Projected, today: `116389 triangles from a 5051 x 5051 node DEM (10 m).
  Every DEM node inside the mesh is within 5 m of it (largest difference
  4.9997 m). 199 vertices on NoData cells were removed with their triangles.`
- Reprojected: `983 triangles. Every node of the original DEM inside the mesh
  (21615 checked) is within 5 m of it (largest difference 4.9926 m).`
- After 15f-3 a sentence is added: `So is every checked point along the lines
  (1834 points, largest difference 4.21 m).` An exception adds a sentence
  naming its count and largest difference, for example `1 line point could not
  be added; its difference is 5.4 m.`

Numbers in the summary are rounded to 5 significant figures; the fields keep
every digit.

### `--stats` only (the "Result" section)

| plain wording | today's name | means |
|---|---|---|
| DEM nodes with data left outside the mesh (self-check, always 0) | `uncovered`, "valid DEM nodes not covered" | DEM nodes with a height that ended up in triangles removed for NoData; the refinement's stopping rule makes it 0 |
| Refinement rounds | `rounds` | passes over the mesh |
| Points inserted | `inserted` | vertices refinement added |
| … of them into triangles with a NoData corner | `carved` | |
| Edge flips | `flips` | |
| Points added to improve the starting mesh; attempts skipped | `quality inserted`, `quality skipped` | |
| Points placed on lines instead of beside them; refused | `feet`, `feet refused` | |
| Starting mesh vertices, triangles; vertices not on a DEM node | start vertices, start triangles, `start vertices off-node` | the last is on stderr only today |
| Largest error against the resampled grid | phase 1's `max_error` | reprojected path only |
| Points the DEM check inserted; its rounds | `37 inserted in 5 rounds` | reprojected path only |
| Line points inserted; DEM nodes inserted by the line check; line points that coincide | `strip_inserted`, `nodes_inserted`, `duplicates` | after 15f-3 |

### stderr, reworded

| today | new |
|---|---|
| the refine report (one long line) | the `summary` |
| `199 vertices without data dropped` (stride path) | the `summary` |
| `mosaic of 2 tiles, 256 x 563 nodes` | in the `summary` (`a 563 x 256 node DEM from 2 files`) |
| `60 features kept, 9 dropped outside, 19 clipped, 0 empty skipped` | `features: 60 kept (19 cut at the domain outline), 9 outside the domain, 0 empty` |
| `<table>: no R-tree index, table scanned` | `<table>: this layer has no spatial index, so every row was read` |
| `10802 input vertices, 5627 noded vertices` | `lines: 10802 vertices read, 5627 after joining shared edges and adding crossings` |
| `land cover: 68 regions, 0 outside every polygon, 0 in more than one, 0 thinner than the snap` | `land cover: 68 areas between lines; 0 in no polygon, 0 in more than one (the smallest wins), 0 too narrow to label with certainty` |
| catchment: `window 1: x ..., flood 0.12 s, contained` | `window 1: x ..., 130 x 158 nodes, searched in 0.12 s, catchment inside the window` (`grown on north, east` → `window widened to the north, east`) |
| catchment: `seed: the pour node at (x, y) (a pour point must lie on the flow line; it is not snapped)` | `start: the outlet node at (x, y) (an outlet must lie on the flow line; it is not moved there)` |
| catchment: `catchment: N nodes, X km2 of node area` | `catchment: N DEM nodes, X km2` |
| catchment: `fine outline: ..., 2 rings dropped (14 nodes), 1 holes filled (...)` | `outline along DEM cells: ..., 2 separate patches left out (14 nodes), 1 enclosed gap filled (...)` |

## Compatibility

Who reads today's text (`git grep -n elevation_source`, and the clauses'
words, over the tree):

- **Tests** (rewritten by `@tester` in 24's red step):
  `test_cli_mesh_dem.py` (`:119-272`), `test_cli_mesh_refine.py` (`:59` and
  its `field`/`sentence` helpers), `test_cli_mesh_mosaic.py` (`:140-518`, the
  `mosaic of` prefix, and `:266-271` `dem_seams`), `test_cli_mesh_geographic.py`
  (`:232-674`, `checked against`, `source_crs`, `computation_grid`),
  `test_cli_mesh_cache.py` (`:211-220`), `test_cli_fetch.py` (`:297-308`),
  `test_cli_mesh_vtk.py` (`:193`), `test_cli_mesh_domain.py`,
  `test_cli_mesh_domain_crs.py`, `test_cli_mesh_features.py`,
  `test_cli_start_quality.py`, `test_cli_constraint_feet.py`,
  `test_stats.py` (the Refinement table, `:235`, `:324-330`, `:399`). `test_io_vtk_legacy.py` and
  `test_io_vtk_readback.py` use `elevation_source` only as a sample field name
  for the writer; they need no change.
- **`tools/bench.py`**: does not read the text. Its numbers come from its own
  `BENCH` line, written from `refine`'s outcome (`tools/bench.py:747-752`), and
  its mesh hash covers the bytes from `POINTS` on (`:333`, `:350`), so header
  fields do not change it. Stored baselines stay comparable.
- **`docs/benchmarks/2026-09-26/bench1m/summarize.py:18`** parses `achieved
  max error` from that day's logs. It is dated evidence for those logs and
  keeps working on them; it is not changed.
- **Palettes** (`palettes.py`): read nothing from the file.
- **Design records** that quote the sentence (12, 13, 14, 14b, 15, 15c, 16,
  16b, 20, 20b, 23) describe what shipped then and stay as they are. The one
  design not yet built that prescribes a clause, `15f-edge-strip.md` D7, gets a
  pointer to D5 here in this branch.
- **Old files** stay readable: nothing in rasputin reads a mesh file back, and
  ParaView and QGIS show unknown FieldData and comments as text. New files do
  not carry `elevation_source`, so a reader that wants both looks for
  `max_error_m` first and falls back to the sentence.

## Tests for `@tester`

Not invariant-critical, so no mutation round (README, "Cost constraints").

- **R1, the record alone** (`tests/python/test_run_record.py`, no `_core`):
  every name matches `^[a-z][a-z0-9_]*$` and is unique; every value is ASCII
  with no control character; D2's omission rules (each zero-omitted field
  absent from `file_fields` at 0, present at 1, always in `stats_rows`); a
  self-check at 1 is in `file_fields`; `max_error_m` printed through
  `_exact` never reads above `tolerance_m` when it is not above it; every
  number in `summary` equals a field's value rounded to 5 significant
  figures; no banned word (D2 rule 6) in any value or in `summary`.
- **R2, one vocabulary** (CLI, projected fixture with NoData): every file field
  is in the `--stats` Result section with the same value; stderr's summary
  equals the `summary` field; the `.ply` comments equal the `.vtk` fields,
  name for name (D3); no `elevation_source`, no `elevation` comment.
- **R3, the values are right**, each against an independent count, not the
  producer's record: `nodata_vertices_removed` equals the number of input
  vertices whose bilinear sample is NoData, counted in NumPy from the DEM
  array; `max_error_m` ≥ the largest \|DEM − mesh\| at DEM nodes inside the
  written triangles, computed in the test; `dem_grid` matches the DEM's shape.
- **R4, reprojected** (velhas): `max_error_m` is the error at the source DEM's
  nodes (an independent check over the fixture's nodes, as 15c-2's tests
  build it), not phase 1's; `dem_crs`, `dem_transform`, `resampled_grid`,
  `dem_nodes_checked` present; `--stats` "resampled grid nodes".
- **R5, the other paths**: stride (no tolerance fields), `--flat` (`heights`
  and `summary` only), mosaic (`dem_tiles`, `dem_seams` both ways), cache
  (`dem_source`, `dem_credit`, `licence_note`, `cite`), features
  (`features`, `land_cover_codes` reworded).
- **R6, stderr wording**: each line in "stderr, reworded", by `re.search` on
  its parts, and the banned-word list over the whole of stderr.
- Existing tests in "Compatibility" are rewritten to read fields.

## LOC estimate

`CLAUDE.md` §2's unit; tests excluded.

| file | what | net est. |
|---|---|---:|
| `run_record.py` (new) | `Entry`, `RunRecord`, three builders, `file_fields`, `summary`, `stats_rows` | 130 |
| `cli.py` | `_dem_mesh`'s sentence and report replaced by the builder call; the field assembly in `mesh` (one list for `.vtk` and `.ply`); renames; reworded stderr lines; catchment wording | −20 |
| `stats.py` | the Result section from `stats_rows` in place of `_refinement`; the Sizes label | 10 |
| **total** | | **120** (+39 %: 167; +60 %: 192) |

One PR. 15f-3 then grows by its five `refined_record` arguments rather than by
a clause, about 10 lines fewer than its D7 estimate.

## Questions for Ola

None blocks the design; the table above is what Ola approves.

- **Units in the name** (`tolerance_m 5`), as the catchment file already does,
  or in the value (`tolerance 5 m`)? Recommended: in the name, because a script
  then reads a bare number.
- **The refine counters on stderr** (rounds, flips, insertions) move to
  `--stats` only. Recommended. If Ola wants them on stderr still, they print
  after the summary as one line per counter, in the Result section's words.
