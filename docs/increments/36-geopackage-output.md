# Increment 36: GeoPackage output

**Status: **built** on branch `worktree-gpkg-out` (design `01118660`, fixes
`0233d1e3`, design review rounds 1 and 2 recorded; red `87f41c3d` (`@tester`,
42 tests), green `c59246e5` + `0943a7a5` (`@developer`), layering rows
`4e7901ee`, mutation round `dddce8e8` (ten mutants, all killed, record in its
commit message), fix-up `7085dddf` (the `.part` housekeeping moved to the
CLI); 377 net production lines by `python3 tools/count_loc.py 81a551b3 HEAD` (the merge base; `3e6f87eb` gives the same count)
against the design's about 305; section 11 is the as-built record; `@perf`'s
acceptance `57c25074` ACCEPTED). Code review round 1 CHANGES REQUESTED (one blocker, Review section), fixed: red test `55cc3c3e` (`@tester`), guard `18442b27` (`@developer`), 378 net lines; code review round 2 APPROVED. Next: the push on Ola's yes, then a merge of master (#229, increment 35) before the queue.

## 1. What Ola asked

`rasputin mesh --out X.gpkg` writes a GeoPackage, chosen by the suffix of
`--out` as `.vtk` and `.ply` are, so QGIS opens a mesh as any vector layer
(Ola, 2026-10-05, ROADMAP row 36). The file is also for agents that cannot
read VTK (Ola, 2026-10-11, "intent-driven triangulation"): a script with
`sqlite3` and shapely, or an agent, checks a mesh from the file alone —
triangle count, per-triangle heights, land-cover code, plan area, and the
run record (tolerance, DEM source). Specification by example: Ola's stopgap
converter `../rasputin_data/lands/scripts/vtk_to_gpkg.py` (outside git),
which this increment replaces.

## 2. Prior art: legacy and literature

*Literature.* The file format is OGC GeoPackage 1.3 (OGC 12-128r18): the
geometry blob header (section 2.1.3, which `io/geopackage.py` already
decodes), the R-tree spatial index extension (Annex F.3), `gpkg_contents`,
`gpkg_geometry_columns`, `gpkg_spatial_ref_sys` and attribute tables
(`data_type = 'attributes'`, section 2.3). The spatial index is SQLite's
R*Tree module (sqlite.org/rtree.html), fed row by row; the GeoPackage
extension's triggers are left out, since the file is written once and the
extension's `scope` column says `write-only`, as the stopgap does. The stored
QGIS style is QGIS's own `layer_styles` table (QGIS Desktop User Guide,
"Save style to the database"), as the stopgap writes it and QGIS 3.28+ reads
it. Bulk-loading an SQLite R-tree is known and faster: Even Rouault's
`sqlite_rtree_bulk_load` builds the tree in memory and writes the `_node`,
`_parent` and `_rowid` shadow tables directly (SQLite forum,
sqlite.org/forum/forumpost/af937a36b3, 2023-10-22: 3.2 million rows, 26 s
to 6 s; shipped in GDAL 3.8), packing after Leutenegger, Lopez and Edgington,
"STR: a simple and efficient algorithm for R-tree packing", ICDE 1997. Not
built here (section 7: 36b), because the row-by-row index costs 2.5 s on the
quick check's mesh, and the shadow-table format is SQLite's private layout.
No novelty is claimed: this is an encoder for a published format.

*Legacy.* Nothing. The archived tree has no GeoPackage, SQLite or QGIS code:

    $ git grep -l -i -E "gpkg|geopackage|sqlite|qgis" legacy-archive -- legacy
    (no files; exit 1)

*In this tree.* `src_python/tin_engine/io/geopackage.py` (16b R2, R3) reads
the blob header, `gpkg_contents`, the SRS table and the R-tree; the writer is
its mirror and the reader reads the writer's files back in the tests.
`tests/python/gpkg_fixtures.py` (`write_gpkg`) already writes small test
GeoPackages with `sqlite3` and is kept as a test helper, not promoted.

## 3. Measured on this worktree (`docs/increments/36-probes/gpkg_write_probe.py`)

AC power, Python 3.14, SQLite 3.53.1, shapely 2.2.0; single runs. The probe
reads a binary `.vtk` and writes the `triangles` and `lines` layers with
section 4's attributes; `--stopgap` runs Ola's converter on the same mesh.

Every "candidate" figure is one command, from the worktree root:

    /usr/bin/time -l .venv/bin/python docs/increments/36-probes/gpkg_write_probe.py \
        <mesh>.vtk out.gpkg --rtree executemany --cache-mb 256 --stream --per-batch --vtk-time

`--rtree executemany` feeds the index rows straight into the virtual table,
`--cache-mb 256` sets the page cache for that build, `--stream --per-batch`
builds rows and blobs per batch of 100,000 (`--batch`'s default), `--vtk-time`
re-times the binary `.vtk` write beside it, and `/usr/bin/time -l`'s
`maximum resident set size` is the peak RSS. The stopgap column adds
`--stopgap ../rasputin_data/lands/scripts/vtk_to_gpkg.py`; the "same bytes"
claim below is `--twice`. **The probe's own defaults (`--rtree select`, no
`--cache-mb`, no `--stream`) are the rejected variants**, so a bare run does
not reproduce the table. The meshes, all binary `.vtk` from `rasputin mesh`:
row 1 `rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif
--tolerance 10 --out <mesh>.vtk --binary`; row 2 the `numedalslagen` case of
`docs/benchmarks/quick/cases.toml` with `--out <mesh>.vtk --binary` appended;
row 3 `../rasputin_scratch/norway/romsdalen/romsdalen_slope_2m_25_35.vtk`
(the `romsdal-slope` case, kept outside git). Row 1 re-run with that line
when the fix was written, 2026-10-11: total 0.40 s, R-tree 0.24 s, 20.0 MB,
160 MB peak RSS.

| mesh | triangles | `.vtk` | `.vtk` write (`--stats`) | GeoPackage, candidate | of which R-tree | size | peak RSS | stopgap |
|---|---|---|---|---|---|---|---|---|
| test tile `7908_3_10m_z33.tif`, `--tolerance 10` | 69,571 | 2.6 MB | 0.005 s | 0.40 s | 0.23 s | 20.0 MB | — | — |
| quick check's Numedalslagen (`--tolerance 10`, CORINE) | 653,039 | 29.2 MB | 0.28 s | 4.0 s | 2.5 s | 203 MB | 0.56 GB | 5.1 s, 163 MB |
| Romsdalen slope mesh (`romsdalen_slope_2m_25_35.vtk`) | 3,096,202 | 111.5 MB | about 1 s | 21.5 s | 14.2 s | 922 MB | 2.1 GB | 28.7 s, 779 MB, 4.1 GB |

Per triangle: a row costs about 2.2 µs and 242 bytes, the R-tree entry about
4.6 µs and 56 bytes. The stopgap's file is smaller only because it carries
fewer attributes. The file is 7 to 8 times the binary `.vtk`, as the ROADMAP
row foresaw (Graz: 678 MB against 111 MB): a row stores the triangle's whole
ring (four corners) and a bounding box, where `.vtk` stores each vertex once.

What the writer must do, each measured as the difference it made:

- **One transaction, `executemany` in batches of 100,000 rows** (`BEGIN` once,
  `COMMIT` once; `journal_mode=OFF`, `synchronous=OFF`, as the stopgap). The
  batch size is not sensitive: `--batch 10000` gave the same time on
  Numedalslagen.
- **The R-tree after the rows, into the virtual table directly** (14.2 s on
  Romsdalen) — not through a temporary table (`--rtree select`, 16.6 s).
  Feeding it in Morton order (`--rtree-order morton`) was slower (4.5 s
  against 3.8 s on Numedalslagen): the kernel's triangle order is already
  local. `PRAGMA cache_size` of 256 MB during the index build: 3.8 s to 2.9 s
  on Numedalslagen (without `--cache-mb`, the default cache is 2 MB).
- **Blobs and rows built per batch from array slices**, no whole-mesh list of
  tuples and no whole-mesh blob array: peak RSS on Romsdalen 5.8 GB with a
  list of all rows (neither flag), 3.2 GB streaming rows only (`--stream`),
  2.1 GB streaming blobs too (`--stream --per-batch`; 0.3 GB of that is the
  probe's own VTK read). Checked at 3.1 million
  triangles; at the São Francisco scale (about 10⁸ triangles) the rows would
  take about 4 minutes and the index 8, for a 30 GB file (section 7).
- **Same run, same bytes**: the second write of each mesh was SHA-256
  identical to the first. Nothing calls `now()`: `gpkg_contents.last_change`
  and `layer_styles.update_time` are written as the fixed text
  `1970-01-01T00:00:00.000Z` (question 3).

The blob encoding is one numpy structured array per batch (`GP`, version 0,
flags `0b0011` = little-endian with a 32-byte XY envelope, `srs_id`, then ISO
WKB: byte order 1, type 1003 PolygonZ, one ring of four points), as the
stopgap's `tri_blobs`; `io/geopackage.py`'s `decode_geometry` reads it.

## 4. The file

All layers in the run's CRS. Row `fid` = index + 1, so `fid - 1` is the
triangle's or edge's index in the `.vtk` of the same run; rows are inserted
in that order and never reordered.

**`triangles`** (PolygonZ, `gpkg_geometry_columns.z = 1`), one row per
triangle, in triangle order. Columns: `code INTEGER` (the `.vtk`'s
`land_cover_code`: 0 in no polygon; NULL on every row when the run had no
coded land cover), `land_cover TEXT` (the code's name from the palette, NULL
when the palette lacks the code), `colour TEXT` (`#rrggbb` from the palette;
`#ff00ff` for a code the palette lacks, the preset's NaN colour), `z_min`,
`z_max`, `z_mean REAL` (the three corner heights), `area_m2 REAL` (plan
area), `slope_deg REAL` (the triangle plane's slope, Horn's quantity on a
flat triangle: `atan2(|n_xy|, |n_z|)`). An R-tree index.

**`lines`** (LineStringZ), one row per constraint edge, in edge order, with
`length_m REAL` (plan length) and one `INTEGER` 0/1 column per vocabulary
property set on at least one edge, named as `write_vtk` names its cell
arrays (ruling 6 of increment 13). Written only when the mesh has edges.
No R-tree (67,240 rows on Numedalslagen draw in one scan).

**`land_cover`** (MultiPolygon, 2D), one row per coded polygon of the
features after the clip and the clean-up (`FeatureSet.features` with a
`polygon` and a `code`, the polygons `_land_cover` labels from): `code`,
`land_cover`, `colour`, `area_m2`. Written only when a coded source ran and
at least one polygon lies in the domain; otherwise absent (as built, review
round 1).
22 rows on Numedalslagen; cheap, so kept (the ROADMAP's constraint-polygons
layer). Constraint lines from non-coded sources are already in `lines`.

**`land_cover_classes`** (an attributes table), one row per code present in
`triangles`: `fid INTEGER PRIMARY KEY` (OGC 12-128 Requirement 119: every
attributes table has an integer primary key, the rowid alias), `code`,
`land_cover`, `colour`, `triangles INTEGER`, `area_m2 REAL`. The legend, and the class shares an agent reads with one
query. Written only with coded land cover.

**`run`** (an attributes table), the run record (increment 25) as rows:
`fid INTEGER PRIMARY KEY` (Requirement 119 again), `name TEXT UNIQUE`,
`wording TEXT`, `value TEXT`, `number REAL` (NULL for text entries), in
record order, straight from `RunRecord.entries`
(`name`, `wording`, `value`, `number`); first `rasputin_version` and
`command` as `--record` writes them, then every entry `as_json` writes. So the tolerance
(`tolerance_m`), the largest error (`max_error_m`), the DEM source
(`dem_source`) and the CRS are in the file. `--record` and `--stats` are
unchanged and still work beside a `.gpkg`.

**`layer_styles`** (QGIS), one row each for `triangles` and, when that layer
exists, `land_cover`,
`useAsDefault = 1`: a categorised style on `code` with the palette's colours
and labels, the stopgap's QML with its three `<prop>`/`<Option>` pairs, as
a raw-string template. Written only when the class map names a palette.
QGIS opens both coloured; any other reader colours by `colour`.

**`gpkg_spatial_ref_sys`**: the two reserved rows (-1, 0), EPSG:4326 (the
spec requires it) and the run's CRS: for an EPSG code, `organization =
'EPSG'`, `srs_id` = the code, `definition` = pyproj's WKT1_GDAL (what QGIS
and GDAL read); for a CRS pyproj parses but cannot name, `srs_id = 100000`,
`organization = 'rasputin'`, `definition` = its WKT2; for no CRS or free
text pyproj refuses (`--flat` gallery fixtures with `--crs`), `srs_id = -1`,
the reserved undefined cartesian row, and the text goes in `run` only. The
reader's `layer_info` turns each back into its `crs` text.

## 5. Modules and boundaries

- **`src_python/tin_engine/io/geopackage_write.py`** (new, about 250
  production lines). Pure in the reader's sense: takes an open
  `sqlite3.Connection`, numpy arrays and plain data; opens no file, knows no
  path, no `_core`, no pyproj. Functions: `polygon_blobs(xyz, srs_id)` and
  `line_blobs(xyz, srs_id)` (structured-array encoders, the reader's
  mirror); `create_geopackage(conn, srs: SrsEntry)` (the four `gpkg_*`
  tables and the SRS rows); `write_triangles(conn, vertices, triangles,
  codes, palette, srs_id, batch=100_000)`; `write_lines(conn, vertices,
  edges, edge_masks, vocabulary, srs_id)`; `write_land_cover(conn,
  polygons, palette, srs_id)` (shapely WKB via `shapely.to_wkb(...,
  flavor="iso")`; its header envelope is `minx, maxx, miny, maxy`, the order
  flag `0b0011` declares in 12-128 section 2.1.3 — not the stopgap's
  `poly_layer`, which reorders the bounds twice (its lines 229 and 118) and
  stores `minx, miny, maxx, maxy`); `write_classes`, `write_run(conn, rows)`,
  `write_styles(conn, palette, layers)`; one `register_layer(conn, name,
  kind, srs_id, bounds, z)` for the `gpkg_contents`/`gpkg_geometry_columns`
  rows, and `build_rtree(conn, table, boxes)` with the 256 MB cache set and
  reset. `SrsEntry` is a frozen Pydantic model (`srs_id`, `organization`, `code`,
  `name`, `definition`).
- **`src_python/tin_engine/crs.py`**: `srs_entry(crs: str | None) -> SrsEntry`
  (about 15 lines), the only place pyproj is asked for WKT.
- **`src_python/tin_engine/feature_input.py`**: `ClassMap.palette: str | None
  = None`, set to `"corine"` on the three CORINE maps (4 lines), so the
  writer looks the palette up in `palettes.PALETTES` by name and never by
  parsing the code-system text.
- **`src_python/tin_engine/cli.py`** (about 35 lines): `MESH_SUFFIXES` gains
  `.gpkg`; with a `.gpkg`, `--out-edges` is refused ("edges are the `lines`
  layer of the GeoPackage") and `--ascii` is refused ("applies to .vtk and
  .ply"); the dispatch branch opens `<out>.part` with `sqlite3.connect`,
  calls the writers in one `clock.phase("write: geopackage")`, closes, and
  `os.replace`s onto `out` (so an existing file is replaced, never appended
  to, and a failed run leaves no half file under the asked name). The
  composition root hands over what it already has: `vertices`,
  `surface_mesh.triangles`, `edges`, `edge_masks`, `codes`, `found.features`,
  `record.entries` with `installed_version()` and `_command()`, `srs_entry(crs)`.
  The `land_cover` polygons' shapely objects stay in Python; the C++ core
  never sees the file, the path or the CRS (I/O boundary, CLAUDE.md §2).
- **Dependencies:** none added; `sqlite3` is the standard library. No GDAL,
  OGR or Fiona (`tools/check_prohibited_deps.py` stays green).
- **Async-readiness:** the writer is one synchronous function over a
  connection, so a worker runs it under `asyncio.to_thread`; SQLite
  connections are thread-bound, so the connection is opened and closed in
  the same thread, as the CLI branch does.

## 6. Byte-identical `.vtk` and `.ply`

The increment adds a writer and a suffix; it changes no refinement, mesh or
encoder code, and the `.vtk`/`.ply` branches of the dispatch are untouched.
Pinned by `tests/python/test_refine_golden.py`
`test_the_kartverket_stride_vtk_is_unchanged_by_node_sampling` (SHA-256 of a
whole `.vtk` against `GOLDEN_STRIDE_VTK`) and by the quick check's
`mesh_sha256` per case (`tools/bench_quick.py`, "MESH CHANGED"), run as
section 9 says.

## 7. Left out, with reasons

- **The source-polygon ID per triangle (cell ID).** Not in the mesh today
  (ROADMAP row 36: its own small step, which also gives `.vtk` and `.ply` the
  attribute). Not folded in here (question 4). Meanwhile a reader joins
  `triangles` to `land_cover` by containment.
- **36b, bulk R-tree loading** (section 2's method): the index is 60 % of
  the write time at every scale measured; at basin scale (10⁸ triangles) it
  would be about 8 minutes and the file about 30 GB, so a basin is written
  per piece or not as a GeoPackage at all. Built only if Ola wants basin
  meshes in QGIS; measured expectation from Rouault's figure: the index 4
  times faster.
- **No `--no-index`, no 2D variant, no attribute selection flags.** The
  measured savings are 17 % of the file (index) and 9 % (2D), not worth a
  flag each; a reader that wants less drops tables with SQLite.
- **A second file (`--out-edges`-like).** Refused; edges are a layer.
- **Reading a mesh GeoPackage back as `--features`.** Its `land_cover` layer
  already reads (16b's reader, `--features-layer land_cover --features-map
  corine` needs the attribute `Code_18`, which it is not named); an idea,
  not this increment.

## 8. Tests (`@tester`), red first

Fixtures with measured figures: **F1**, a hand-made mesh of 4 vertices and
2 triangles on the unit square, vertices (0,0,0), (1,0,1), (1,1,2), (0,1,3),
triangles (0,1,2) and (0,2,3), one constraint edge (the diagonal 0-2) with one
property bit, codes 312 and 0, written with `srs_id` 25833 (F2's EPSG:25833;
not 0 or -1, which read the same in either byte order): plan areas 0.5 and
0.5, z_mean 1.0 and 5/3, slopes 54.7356° (`atan(√2)`) and 72.4516°
(`atan(√10)`), computed with the probe's formulas on this worktree; **F2**, the DEM tile
`tests/fixtures/dem_archive/7908_3_10m_z33.tif` with `--tolerance 10`:
69,571 triangles, 383 constraint edges, 38,921 points, 0.48 s to mesh, the
GeoPackage 20.0 MB and 0.40 s (section 3); **F3**, the land-cover case
`test_cli_mesh_landcover.py` already runs (its `bumpy` DEM with hand-made
CORINE polygons), for the coded path.

Invariant-critical suite (mutation round required), with its targets:

- **T1 round trip** (F1, `:memory:`): `layer_info` + `query_features` from
  `io/geopackage.py` read back each triangle as a PolygonZ equal to the
  corner coordinates, `fid` in triangle order, `srs_id` 25833; and, because
  `decode_geometry` skips the envelope (it reads from byte 8 + 32), T1 also
  unpacks each blob's header itself: bytes 0-4 `GP`, 0, `0b0011`; bytes 4-8
  `<i` = 25833; bytes 8-40 `<4d` = the triangle's own `minx, maxx, miny,
  maxy`. Oracle per mutant: byte-order flag 0 — the reader (`>i` on the
  little-endian bytes gives -379322368, not 25833); envelope order swapped
  (`minx, miny, maxx, maxy`) — the header unpack, since the reader never
  sees it; type 3 for 1003 — the reader (Polygon, `has_z` false, where
  PolygonZ is asserted); ring not closed, and point count 3 (the first three
  corners do not close) — the reader raises `GeoPackageError` (GEOS: "Points
  of LinearRing do not form a closed linestring").
- **T2 index matches rows** (F2, written with `batch=30_000`, so the 69,571
  triangles cross two batch boundaries and leave a partial batch of 9,571;
  under the default 100,000 no boundary is crossed and the batch mutant
  is never reached): for every `fid`, the R-tree box equals the triangle's
  bounds to float32 rounding, and the table and the index have the same row
  count. Mutants and their oracles: `fid` off by one in the index — the
  per-`fid` box comparison; the last partial batch dropped — the row count
  (60,000 against 69,571); `miny`/`maxx` swapped — the box comparison.
- **T3 attributes** (F1 and F3): `code`, `land_cover`, `colour`, `z_*`,
  `area_m2`, `slope_deg` equal the figures computed from the arrays; a code
  the palette lacks gives NULL name and `#ff00ff`; a run with no coded land
  cover gives NULL `code` on every row and no `land_cover`,
  `land_cover_classes` or `layer_styles`. Mutant: areas without the ½.
- **T4 determinism** (F2): two runs give SHA-256-identical files;
  `last_change` and `update_time` are the fixed text; no row of
  `gpkg_contents` or `layer_styles` differs.
- **T5 `run` table**: every entry of the record is a row, numbers in
  `number`, `rasputin_version` and `command` present; `land_cover_classes`
  sums to the triangle count and to the sum of `area_m2`.
- **T6 `lines`**: row count = edge count; each property column matches
  `edge_masks` through the vocabulary. Mutant: a column for a property set
  on no edge.
- **T7 SRS**: EPSG → organization `EPSG`, `definition` non-empty WKT and
  `layer_info(...).crs == "EPSG:25833"`; a WKT CRS without a code → 100000;
  `--flat --crs "free text"` → -1 and the text in `run`.
- **T8 CLI** (F2): `--out x.gpkg` writes the file and echoes its path;
  `--out-edges` and `--ascii` with `.gpkg` are refused with the sentences of
  section 5; an existing `x.gpkg` is replaced, not appended to (write twice,
  one `triangles` table of F2's count); `--record` and `--stats` beside it
  work; the `.vtk` of the same run is byte-identical to a run without the
  GeoPackage (section 6's pins stay).
- **T9 style**: the QML in `layer_styles` parses as XML with one `<category>`
  per palette code and `attr="code"`.

Not invariant-critical, no mutation round: T7 to T9.

## 9. Estimate, speed check, acceptance

- **Production lines:** about 250 (`io/geopackage_write.py`) + 15 (`crs.py`)
  + 4 (`feature_input.py`) + 35 (`cli.py`) ≈ **305 net**, under 700. The
  QML template is a raw literal and does not count.
- **Speed:** could change speed? No for `.vtk` and `.ply` runs (one
  `if` on the suffix); the new path is its own figure. The quick check runs
  before the push and its `mesh_sha256` must not change.
- **Acceptance:** (a) `@perf`'s quick check unchanged, every `mesh_sha256`
  equal; (b) the GeoPackage of the quick check's Numedalslagen case written in
  **at most 5 s** on AC power (candidate 4.0 s, stopgap 5.1 s) and of the
  Romsdalen mesh in at most 25 s (candidate 21.5 s), peak RSS under 1 GB and
  3 GB respectively, measured by `@perf` with `/usr/bin/time -l` on the built
  command; (c) QGIS opens the Numedalslagen file with `triangles` coloured by
  land cover and `land_cover` beneath it — Ola checks this by hand; (d) an
  agent's check from the file alone: `SELECT count(*) FROM triangles`,
  `SELECT value FROM run WHERE name = 'tolerance_m'`, and
  `SELECT code, sum(area_m2) FROM triangles GROUP BY code` agree with
  `--stats` (`@tester` scripts it in T5).

## 10. Questions for Ola (defaults taken while away)

1. **R-tree index always on?** It is 60 % of the write time and 17 % of the
   file; without it QGIS redraws 3 million triangles by full scans. *Default:
   on, no flag.*
2. **The attribute set** `code, land_cover, colour, z_min, z_max, z_mean,
   area_m2, slope_deg` on every triangle (the text columns cost 25 MB of 203
   on Numedalslagen). *Default: yes, all eight.*
3. **Fixed timestamps** (`1970-01-01T00:00:00.000Z`) so the same run gives
   the same bytes, instead of the write time. *Default: fixed.*
4. **The source-polygon ID** stays its own small step before or after this,
   not folded in. *Default: its own step.*
5. **Stored QGIS style** kept in 36 (about 30 logic lines plus the template).
   *Default: kept.*

## 11. As built (recorded by the main session from the handbacks)

**Red** `87f41c3d` (`@tester`): 42 tests, 0 passing, in
`tests/python/test_io_geopackage_write.py` (24: T1 round trip and header
bytes 0-40 and WKB bytes 40-53; T2 at batch=30_000 on F2's 69,571 triangles;
T3 by hand; T5 `write_run`; T6 `write_lines`; T7 `srs_entry`) and
`tests/python/test_cli_mesh_gpkg.py` (18: T3 figures against the `.vtk`
arrays, `land_cover` layer and envelope order; T4 SHA-256 and 1970 stamps; T5
`run` against `--record`/`--stats`; T6 lines via `feature_mask`; T7 EPSG and
`--flat --crs "free text"` giving -1; T8 refusals, replace, fid-1 equals the
`.vtk` index; T9 QML). Ten pins beyond the design, all in the tests:
`write_triangles` registers and indexes its layer, `write_lines` registers;
`palette` is the code table, not its name; `geometry_type_name` uppercase;
`write_run(conn, rows)` takes `(name, wording, value, number)` tuples, `run`
and `land_cover_classes` registered as `attributes`; `srs_entry` never raises
and -1 adds no SRS row; refusal wording "layer of the GeoPackage" + "lines",
"applies to .vtk and .ply"; an uncoded run's `gpkg_contents` is {triangles,
lines, run}; the QML category and symbol shape; `run`'s CRS row named `crs`;
WKB order 1, types 1003/1002.

**Green** `c59246e5` + `0943a7a5` (`@developer`), fix-up `7085dddf`: new
module `src_python/tin_engine/io/geopackage_write.py`; `.gpkg` dispatch and
`_write_geopackage` in `cli.py` (writes `<target>.part`, `os.replace` to the
target, the stale `.part` unlinked before `create_sqlite`); `create_sqlite` in
`io/repository.py` is `sqlite3.connect(path)` alone, since
`test_io_geopackage.py::TestBoundaries` allows only the repository to connect
and `test_io_repository.py::TestQ3OneModuleOpensFiles` forbids it to delete
(the design's section 5 had the CLI connecting; the tests won). `SrsEntry`,
`UNDEFINED_SRS`, `CUSTOM_SRS_ID` live in `io/models.py`, since the layering
table forbids `crs` importing an io codec and the red suite forbids the writer
importing `crs`. `write_classes(conn)` derives the legend from `triangles` by
GROUP BY. 377 net counted lines (the SQL and QML template bodies count:
`count_loc.py`'s raw-literal rule is C++'s `R"..."`). Whole Python suite 7162
passed, 17 skipped at `7085dddf`; mypy strict, ruff, the governance gates
clean; `test_refine_golden.py` passes, so `.vtk` and `.ply` are untouched.
One real write, Numedalslagen on AC: `write: geopackage` 3.98 s, 203.3 MB,
653,039 triangles, 67,240 lines, 22 `land_cover` rows (design 4.0 s /
203 MB). Peak RSS of the whole command 2.79 GB, meshing included; acceptance
(b)'s 1 GB is the writer's own share, which `@perf` judges.

**Layering** `4e7901ee` (`@tester`): `io.geopackage_write: features
io.geopackage io.models` in L2; `crs` gains `io.models`; `cli` gains
`io.geopackage_write`; each edge checked against the import lines.

**Mutation round** `dddce8e8` (`@tester`), scratch copy of `0943a7a5`, green
before, between and after: ten design-named mutants, all killed, each by the
design's named oracle. The full record, one row per mutant with the planted
line, the killing test and assertion, is the commit message of `dddce8e8`
(the lean rule: the kill record may live in the commit the handback names).
Short form: M1 flag 0 (reader: srs_id -379322368), M2 envelope swapped (header
unpack), M3 type 3, M4 ring unclosed, M5 points 3 (reader `GeoPackageError`),
M6 fid off by one, M7 partial batch dropped (69571 against 60000 at
batch=30_000), M8 miny/maxx swapped, M9 area without the half, M10 a column
for an unused property. Deviation: M3 dies by GEOS "ring not closed" and the
WKB type assertion, not the design's predicted flat Polygon. Red scaffolding
removed in the same commit: the `importlib.import_module` fixtures and the
red notes in both test files.

**Review-round-1 fix** `55cc3c3e` (`@tester`) + `18442b27` (`@developer`):
`TestLandCover::test_a_coded_run_with_no_polygon_in_the_domain_writes_without_a_land_cover_layer`
(a corine-mapped GeoJSON holding one LineString across the domain, so the clip
leaves no polygon): the `.gpkg` run exits 0, `land_cover` is absent as a table
and from `gpkg_contents`, `triangles.code` is 0 everywhere and equals the
`.vtk`'s, `land_cover_classes` has exactly the code-0 row. Guard in `cli.py`
only: `write_land_cover` runs `if cover:`, and `write_styles` gets
`"land_cover"` only then, so no `layer_styles` row names an absent layer (the
default taken; the test pins neither way). Whole suite 7163 passed, 17
skipped; gates clean; 378 net counted lines (+1).

**Decided on defaults while Ola was away** (section 10): R-tree always on;
all eight triangle attributes; fixed 1970 timestamps; the source-polygon ID a
later step; the stored QGIS style kept.

**Lessons for the retrospective:** a design naming where `sqlite3.connect`
lives should check the existing boundary guards first; housekeeping (`unlink`,
rename) belongs in the caller from the start, as the `.vtk`/`.ply` writers do.


**Acceptance** `57c25074` (`@perf`, ACCEPTED against section 9;
`docs/benchmarks/2026-10-11/gpkg-out/README.md`), AC power, head `e2470ae7`
against base `81a551b3` back to back (the worktree's stored baseline was on
numpy 2.5.3, so the base was measured with `--tree`), medians of 3: (a) all
7 quick-check meshes' SHA-256 equal to base and to master's refreshed
baseline; quick check `NO CHANGE`, -4.3 % to +1.6 %; (b) Numedalslagen
`write: geopackage` 3.932 s (limit 5 s), the writer's memory share +32 MB on a
2.77 GB command peak set by meshing (`.vtk` run 2.73 GB); Romsdalen write
19.368 s (limit 25 s), 810 MB file, writer's share +188 MB on 2.91 GB (`.vtk`
2.78 GB). Against master's refreshed baseline Lagan reads +37.7 %: the base has
the same features clip (11.58 s), so that is increment 35 missing from this
branch, not a regression. Pinned: the writer's share is RSS growth inside the
`write: geopackage` phase from a 0.1 s `ps` trace; the production Romsdalen
file is 810 MB where the design probe wrote 922 MB, not investigated. Lesson: a
branch based before a baseline refresh needs the back-to-back base run anyway.

## Review

Design review round 1, 2026-10-11, @reviewer (01118660..759504ed, prose only, 0 counted lines): CHANGES REQUESTED. Four design fixes, each one or two lines: (1) T1's kills do not hold at /Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/docs/increments/36-geopackage-output.md@759504ed:245-249, since the reader skips the envelope bytes (/Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/src_python/tin_engine/io/geopackage.py@759504ed:67-73) and the byte-order mutant survives an `srs_id` of 0 or -1; fix: T1 unpacks bytes 8-40 against the bounds, F1 uses `srs_id` 25833. (2) T2's batch-boundary mutant at `:252` is never crossed by F2's 69,571 triangles under a 100,000 batch; fix: T2 passes a small `batch=` (the signature at `:164` allows it). (3) `run` with `name TEXT PRIMARY KEY` (`:132`) and `land_cover_classes` without a key (`:126`) contradict OGC 12-128 Requirement 119 (attributes tables need an INTEGER primary key as rowid alias); fix: `fid INTEGER PRIMARY KEY`, `name TEXT UNIQUE`. (4) Section 3's figures (`:53-63`) name no command line, and the probe's defaults are the rejected variants; fix: the flags per table row. Verified sound: the figures match 01118660 and the derived rates; F1's areas, means and slopes by hand; the GP header flags 0b0011 agree with the reader and the stopgap; the byte-identity pins exist (tests/python/test_refine_golden.py@759504ed:190,202; tools/bench_quick.py@759504ed:132); the GDAL 3.8 R-tree figure confirmed; dependency gate green. Lesson: a mutant claimed against a reader that skips the mutated bytes is not a kill; name the oracle per mutant. Idea: the stopgap packs polygon-layer envelopes as minx,miny,maxx,maxy (wrong order); `write_land_cover` should state its order rather than copy it.

Design review round 2, 2026-10-11, @reviewer (38adbb87..0233d1e3, prose only, 0 counted lines): APPROVED. All four round-1 fixes in and correct: T1 unpacks the header bytes itself with an oracle per mutant matching the reader (/Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/src_python/tin_engine/io/geopackage.py@0233d1e3:63-82), F1 at srs_id 25833 (/Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/docs/increments/36-geopackage-output.md@0233d1e3:266-267, :279-290); T2 with batch=30_000 crosses two boundaries and a partial batch (:291-298); `land_cover_classes` and `run` have `fid INTEGER PRIMARY KEY`, `run.name` TEXT UNIQUE (:154, :160); section 3 gives the probe command and flags per row with the rejected variants' flags beside their figures (:61-80, :96-111); the row-1 re-run is stated as such (:79-81); `write_land_cover` states its envelope order (:196-199). check_citations exit 0. Idea for @tester: word T1's type-code oracle as "raises or `has_z` false", since a type-code-only mutant dies by GeoPackageError.

Code review round 1, 2026-10-11, @reviewer (81a551b3..e50ea446, 377 net counted lines; not pushed, no CI yet): CHANGES REQUESTED. One blocker: a coded run whose source leaves no polygon inside the domain crashes the GeoPackage write; `label_triangles` accepts an empty polygon list (/Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/src_python/tin_engine/landcover.py@e50ea446:143), so `cover` is `[]` (/Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/src_python/tin_engine/cli.py@e50ea446:1146) and `write_land_cover(conn, [], ...)` reaches `_bounds` on a zero-row array (/Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/src_python/tin_engine/io/geopackage_write.py@e50ea446:172, :360), `ValueError: zero-size array to reduction operation minimum`, probed in the venv; the same input writes `.vtk` fine. Fix: guard with `if cover:` at cli.py:1273 or return early in `write_land_cover`, plus one test (T3: coded source, no polygon in the domain; `land_cover` layer absent or empty, `land_cover_classes` still written); section 4's "written only when a coded source ran" gets the same word. Verified sound: both gates green; check_citations exit 0 with the at-risk quotations re-read; 83 tests pass; no red scaffolding; the kill record in dddce8e8 covers all ten targets and every named test exists; acceptance figures equal the raw records; layering rows match the imports; `.vtk`/`.ply` writers untouched. Question decided on default (yes): acceptance (b)'s memory limits read as the writer's own share (+32 MB, +188 MB; the command peaks at 2.77/2.91 GB from meshing alone). Two prose nits applied in this recording commit: the Status line's base is the merge base `81a551b3`; `SrsEntry` is a frozen Pydantic model, not a dataclass. Idea for the retrospective: `docs/increments/30d-outline-buffer-speed.md:103` cites `feature_input.py:158-167` unpinned, now off by three. Lesson: a layer written conditionally needs the empty case in the red suite.

Code review round 2, 2026-10-11, @reviewer (e50ea446..3839914f on 81a551b3, 378 net counted lines: cli.py 70, geopackage_write.py 271, crs.py 21, models.py 11, feature_input.py 3, repository.py 2; not pushed, no CI yet): APPROVED. The guard removes the crash path and nothing else (`cover` built at /Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/src_python/tin_engine/cli.py@3839914f:908 and :1146, guarded at :1275 and :1279; a non-empty cover takes the former calls byte for byte, so the acceptance figures of 57c25074 stand without a re-run); the test pins the reworded section 4 (/Users/skavhaug/projects/rasputin/.claude/worktrees/gpkg-out/tests/python/test_cli_mesh_gpkg.py@3839914f:466-490); section 4's sentences match the code; no red-step scaffolding; check_citations exit 0; 43 tests in the two files pass. Ideas, non-blocking: an assert that `layer_styles` has no `land_cover` row in the empty case; the round-1 record's unpinned `cli.py:1273` could carry `@e50ea446`.
