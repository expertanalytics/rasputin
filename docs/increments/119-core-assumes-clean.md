# Increment 119: the core takes valid polygons only

Status: **green** (`@architect`, 2026-11-01, on public master
`ae33748e`; branch `worktree-core-assumes-clean`; amended the same day:
115's switch done, and Ola's ruling that the refusal names the catalogue's
clean sibling, section 3; amended after the red suite: `clean_sibling` takes the CRS label, so `sources` keeps its layer, and `@tester`'s pins adopted, sections 3 and 7). Full lane (it deletes
the land-cover repair on the path to the mesh), Python only: no C++, no
`@perf` acceptance run (section 8). Pipe item `mesh/core-assumes-clean`,
ruled by Ola 2026-10-29 (the core states its contract, refuses invalid land
cover by name, and its repair code goes), after `mesh/clean-land-cover`
(private increment 115, merged 2026-11-01 as private PR #91). Number 119:
118 is the highest under `docs/increments/` over every local and
remote-tracking ref and every worktree.

## 0. The plain picture

- **Today** every `rasputin mesh` run repairs broken polygons of
  `--land-cover`, `--prefer` and `--reservoirs` again, near its domain,
  after moving them into the DEM's CRS, and records how many it repaired.
- **After** the core checks those polygons and repairs none. A polygon that
  reaches the read region and is not valid in the DEM's CRS stops the run
  before the mesh is built: exit code 2, the file named, and every bad
  polygon listed by feature id with GEOS's reason. Valid input gives the
  same mesh, byte for byte.
- **Ola's runs** read clean products made once by the private clean-up
  (increment 115): every polygon in them is valid in the CRS they were made
  for. A run (or an MCP agent) that names a raw source by catalogue id,
  where the catalogue holds its clean copy for the DEM's CRS, is told to
  use that copy ("use ar50-innlandet-clean-25833"). A public user with a
  broken file gets the list and the fix (repair the file in the DEM's CRS),
  never a silent change of their data.

## 1. Prior art: legacy and literature

*Literature.* Validity is the OGC Simple Features definition (ISO
19125-1:2004, section 6.1.11: rings simple and closed, holes inside the
shell, interiors connected), as GEOS's `IsValidOp` checks it and shapely
2.2 exposes it (`is_valid`, `is_valid_reason`; PostGIS's `ST_IsValid` and
`ST_IsValidReason` are the same code). Refusing at ingest with the reason,
instead of repairing, is the choice this item makes; the repair moves
upstream to the private product, which uses the even-odd rule of Ledoux,
Arroyo Ohori and Meijers, *Computers & Geosciences* 66 (2014), 121-131
(increment 115, section 1). No departure, no novelty claimed, none searched.

*Legacy.* `git grep -n -i -E "make_valid|buffer\(0|is_valid|isValid|explain_validity|repair" legacy-archive -- legacy`
returns nothing. Nothing to carry.

*Ideas file.* Completes "Repair CORINE once per area of interest" (taken up
by 115; this item makes the product the only way an invalid source reaches
the mesh). Overlaps "20c-3: close gaps between land-cover polygons": not
taken up; the gap closing (`--edges-repair`) stays (section 2).

*The rulings read.* The brief relays Ola's "1-2: yes" (2026-11-01). This
design reads it as his answer to increment 115's questions 1 and 2: the
border clean-up stays in the core, and lake and reservoir files are cleaned
by the product and refused by the core like land cover. If it answered
something else, section 2's `--prefer` and `--reservoirs` rows are the
ones to revisit.

## 2. What goes and what stays

Every place on the polygon input path that repairs or tolerates an invalid
polygon (`git grep -n -E "make_valid|buffer\(0|orient\(|is_valid|clip_by_rect|repair" -- src_python`,
at `ae33748e`):

| Where | What it does | Verdict |
|---|---|---|
| `src_python/tin_engine/feature_input.py@ae33748e:547` `_repaired` | even-odd repair of an invalid polygon inside a box; refuses a spike GEOS cannot rebuild | **goes** |
| `src_python/tin_engine/feature_input.py@ae33748e:716` `read_polygons` (`--prefer`, `--reservoirs`) | repairs each polygon read; returns the count; `tell` callback | **repair goes**, refusal in its place; returns the list only; `tell` goes |
| `src_python/tin_engine/land_cover_map.py@ae33748e:89` `_fixed` | checks each polygon reaching the box, repairs the invalid ones, counts them | **check and repair go** (the check moves to the reader); the cut to the box stays |
| `LandCoverMap.repaired`, `build`'s count | the count of repaired polygons | **goes**; `build`'s keyword `fix` becomes `cut` (it now only cuts) |
| `src_python/tin_engine/cli.py@ae33748e:2483` `_say_repaired` and its three calls | the stderr line "N of M polygons crossed itself and was repaired" | **goes** |
| run record keys `land_cover_repaired`, `reservoirs_repaired` (`cli.py`, `src_python/tin_engine/run_record.py@ae33748e:100`) | the counts | **go** |
| `shapely.coverage_clean` in `feature_input.py` (`--edges-repair`) | closes gaps and overlaps between valid polygons, at a distance set from the tolerance | **stays** (a border clean-up, not validity; 115 question 1) |
| `sub_catchments.repair` (`--catchments` units) | `coverage_clean` of the units | **stays** (not land cover; section 10) |
| `snap_to_outline`'s and `_loops`' `make_valid` (`src_python/tin_engine/feature_input.py@ae33748e:851`, `:876`) | fixes polygons the outline rule itself rebuilt | **stays** (the core's own output, not input) |
| `orient` in `domain.py:65`, `io/domain_file.py:90`, `cli.py:1836` | sets ring orientation of valid polygons | **stays** (orientation is not validity) |
| `io/domain_file.py:88` | already refuses an invalid domain with `explain_validity` | **stays**; the precedent this refusal follows |
| `polygon_fractions.subdivide`'s `clip_by_rect` | halves valid polygons, can give invalid halves | **stays**, left to `fractions-intersection` (section 9) |

## 3. The contract and the refusal

**The contract**, stated in `read_class_polygons`' and `read_polygons`'
docstrings and in the README's land-cover bullet: *every polygon of
`--land-cover`, `--prefer` and `--reservoirs` whose bounding box meets the
read box is valid (OGC Simple Features, as `shapely.is_valid` checks it)
after its move into the DEM's CRS.* The read box is the one `_fixed` uses
today: the bounds of `source_region(area, dem_crs, dem_crs, cell=cell)`
widened by `MARGIN` (100 m) on each side; for `read_polygons`, the bounds
of its points `xy` widened by `MARGIN`, the box `_repaired` used.

**Why "reaching the read box", not every polygon read.** A polygon whose
box misses the read box is never overlaid: land cover keeps it in
`outside`, read only by its rings (111, section 3, step 6), and
`read_polygons` drops an invalid one there as `_repaired` did (its repair
inside the box was empty). Checking it would refuse a Sula run for a
broken lake 100 km away, and check a whole GeoJSON file that has no index.
So the set checked is exactly the set checked today, at the same cost.

**Where.** In the readers, the I/O boundary, after the move into the DEM's
CRS (a polygon valid in its own CRS can break in the move: 61 of
Innlandet's AR50 polygons do from EPSG:25832 to EPSG:25833, 115 section
3), before `land_cover_map.build`, before `simplify_borders` (the first
`_core` call that sees a border) and before the mesh. `--reservoirs` is
read today after the mesh, with the mesh's vertices; it moves up to be read
and checked right after the land cover, with the read area's exterior
points (`area.polygon.exterior`, as `--prefer` uses), and the overlay with
the triangles stays where it is. The domain's box holds the mesh's box, so
every polygon the overlay meets is still read; one more polygon read can
only meet no triangle.

**One helper**, in `feature_input.py`:

```python
def _refuse_invalid(
    name: str,
    fids: Sequence[Any],
    moved: Sequence[BaseGeometry],      # in mesh_crs
    as_read: Sequence[BaseGeometry],    # in own, for the broken-by-the-move note
    box: tuple[float, float, float, float],  # in mesh_crs
    own: str,
    mesh_crs: str,
    threads: int,
) -> npt.NDArray[np.bool_]:
    """Which of ``moved`` reach ``box``; a FeatureError listing every one of
    those that is not valid."""
```

It takes the vertex-count grouping and the thread pool from today's
`_fixed` (`GROUPS`, so one 100k-vertex polygon does not hold a thread
while the rest wait), runs `shapely.is_valid` on the polygons reaching
`box`, and returns that mask, which `read_polygons` uses to drop the far
ones that are invalid. `read_class_polygons` gains a keyword
`near: tuple[float, float, float, float] | None = None`; None checks every
polygon read (the strict reading; the CLI always passes the box).

**What it says.** One `FeatureError`, raised once every polygon has been
checked, so a run lists all of them, sorted by fid:

    {name}: {n} of the {m} polygons reaching the read region are not valid in
    {mesh_crs}, the DEM's CRS; the mesh takes valid polygons only. Repair them
    in {mesh_crs} (shapely.make_valid, or QGIS's Fix geometries) and run again:
      feature {fid}: {shapely.is_valid_reason}
      feature {fid}: {reason} (valid in the file's CRS, {own}: broken by the move)

The second form is for a polygon valid as read and invalid only after the
move; it tells a user who made a clean product in the wrong CRS what
happened. `{name}` is the file's name (`path.name`), as the readers' other
refusals have it; `{own}` is the file's CRS text as the reader got it (a
GeoJSON `crs` member `EPSG:25832` gives `EPSG:25832`). `{m}` counts the
non-empty polygons reaching the read region: empty rows are skipped before
the check (section 6) and are not among them. Where the head's lines break
is free; each bad polygon is a line of its own that, stripped, reads
`feature <fid>: <reason>` with an optional note. A thin `--prefer` spike,
refused by `_repaired` today as one GEOS cannot rebuild, is now one more
invalid polygon: "lake.geojson: 1 of the 1 polygons …". `fid` is the reader's: the GeoPackage key, or the 0-based
feature position in a GeoJSON file (as every other core message names it).
Every polygon is listed, none held back: the raw Telemark AR50 file lists
460 on Numedalslågen (section 5). The words must not contain "given CRS":
the CLI sends a message holding that phrase to `--land-cover-crs`
(`src_python/tin_engine/cli.py@ae33748e:2275`).

The helper raises `InvalidPolygons`, a subclass of `FeatureError` in
`feature_input.py` (two lines), so every existing `except FeatureError`
still holds and the CLI can tell this refusal from the reader's others.

**The clean sibling** (Ola's ruling, 2026-11-01, relayed by the main
session: when the refusal names an invalid source by catalogue id and the
catalogue holds `<id>-clean-<EPSG code>` for the DEM's CRS, the refusal
names that sibling as the fix). One pure function in `sources.py`:

```python
def clean_sibling(key: str | None, label: str, known: Mapping[str, RemoteSource]) -> str | None:
    """`<key>-clean-<n>` when `label` is `EPSG:n` and `known` holds that id;
    None for no key, a label that is not `EPSG:<digits>`, or no such entry."""
```

`label` is the DEM CRS as `crs.crs_label` gives it; the caller, `cli._invalid`,
computes it. `sources` is layer 0 of the import map
(`tests/python/test_layering.py`: its row is `toml_data` alone, and its
docstring promises no pyproj), so it reads the label as text and never
imports `tin_engine.crs`. A PROJ string or WKT label does not start
`EPSG:` and gives None, so the design's test 11 holds unchanged. Not moved
to `cli.py`: the `<id>-clean-<n>` convention is the catalogue's, and test 11
stays a pure call that loads neither the CLI nor pyproj.

The lookup is by id only: the key is the catalogue id the user gave
(`_Vector.source.id` for `--land-cover` and `--prefer`; the
`--reservoirs` value when it is in `known`), never a path, file name or
directory searched or guessed. A path source has key None and keeps
today's text. The sibling's own entry is not inspected (its kind, classes
or path): the id is the whole convention. An id that already ends in
`-clean-<n>` gets no stripping: `ar50-innlandet-clean-25832` on an
EPSG:25833 DEM looks up `ar50-innlandet-clean-25832-clean-25833`, finds
none, and keeps today's text with its "broken by the move" lines.

In `cli.py`, one helper serves the three `except InvalidPolygons` sites
(placed before each site's `except FeatureError`):

```python
def _invalid(exc: InvalidPolygons, key: str | None, dem_crs: str,
             known: Mapping[str, RemoteSource], flag: str) -> typer.BadParameter:
```

It calls `clean_sibling(key, crs_label(dem_crs), known)` (`cli` already
imports `crs`). With a sibling, the message is the refusal with one
sentence in front:

    use {sibling} in place of {key} (its clean copy for {EPSG:n} in the catalogue). {refusal}

so "Invalid value for --land-cover: use ar50-innlandet-clean-25833 in
place of ar50-innlandet (its clean copy for EPSG:25833 in the catalogue).
ar50-innlandet: 7 of the 1454 polygons …". The fix goes first because the
list can run to hundreds of lines (460 on Telemark) and an agent reads the
head; the MCP server passes the child's stderr whole
(`src_python/tin_engine/mcp/server.py@ae33748e:310`), so the same words
reach an MCP agent. Without a sibling the message is the refusal alone.
`_land_cover_map` gains the parameter `known` (the run's catalogue, already
in hand at `src_python/tin_engine/cli.py@ae33748e:1026`) for the lookup.

**The exit.** The CLI's existing `FeatureError` handling turns it into
`typer.BadParameter` with `param_hint` `--land-cover`, `--prefer` or
`--reservoirs`: "Invalid value for --land-cover: …", exit code 2, no
output file written. The MCP tools get the same words as their tool error.

## 4. Blueprint: the data flow after

    read_source (own CRS) -> _moved (DEM CRS)
      -> _refuse_invalid(box)       # FeatureError: every invalid one reaching the box
      -> ClassPolygons.areas         # land cover: all rows, valid where they reach the box
      -> build(..., cut=True)        # _fixed: cut to the box, the rest to `outside`; no check
      -> map.areas                   # valid, overlaid unchecked as today
    read_polygons -> list[BaseGeometry]   # valid near and far, file order; invalid far dropped

The interfaces the tests call: `read_class_polygons(..., near=box)`,
`read_polygons(path, layer, mesh_crs, xy, named=..., lines=...) -> list[BaseGeometry]`,
`land_cover_map.build(base, preferred, classes, region, *, cut=True, outside=(), threads=1)`,
`LandCoverMap` without `repaired`, `sources.clean_sibling(key, label, known)`,
`feature_input.InvalidPolygons`, and `rasputin mesh` through `CliRunner`.

## 5. Ola's runs, and a public user

Measured with a prototype of the check (scratch `refusal_probe.py`: the
core's `read_source` with the CLI's read box, `_moved`, the `_fixed` box
test, `is_valid` on one thread; shapely 2.2.0, GEOS 3.14.1, M1 Max on
battery, tree `ae33748e`). Invalid polygons reaching the read box, raw
source against the product the private `clean-land-cover` writes (made in
the scratchpad for this probe, nothing written under `$RASPUTIN_DATA`):

| Case (DEM CRS EPSG:25833 unless said) | Raw source | Product in EPSG:25833 |
|---|---|---|
| Numedalslågen quick check, `corine2018_dtm10_utm33.gpkg` | 0 of 7,239 (is_valid 1,080 ms) | not needed |
| Lagan quick check, EEA CORINE (EPSG:3035 to EPSG:3006) | 0 of 7,150 (1,681 ms) | not needed |
| Sula, `clc2018` cache | 4 of 31 | 0 of 31 |
| Sula, `nve-lakes` | 0 of 116 | not needed |
| Bygdin, `clc2018` | 6 of 191 | 0 of 192 |
| Bygdin, `ar50-innlandet` | 7 of 1,454 | 0 of 1,599 |
| Bygdin, `nve-reservoirs` | 0 of 6 | not needed |
| Numedalslågen, `ar50-buskerud` / `ar50-telemark` | 296 of 59,518 / 460 of 85,639 | 0 of 60,958 / 0 of 89,944 |

So both quick-check land-cover cases pass unchanged, and every product
passes. Whole-file `nve-lakes` and `nve-reservoirs` caches hold no invalid
polygon at all; `smhi-lakes` (the file `smhi-reservoirs` filters) holds 2
of 42,178. The products hold more polygons near the box than the raw files
(Innlandet 1,599 against 1,454): the raw file is queried with the read
region's bounds in EPSG:25832, which at the corners of a box rotated
between UTM zones misses polygons the EPSG:25833 box test counts. None of
those 146 meets the read region (scratch `near_sets.py`: 0 of 146), so no
mesh change follows from it; 115's proof run (below) checks the mesh.

**Before this merges**: increment 115's switch, which this blocked on, is
done (main session, 2026-11-01). Checked the same day: 16 products under
`$RASPUTIN_DATA/clean/<id>_25833_20261101/<id>.gpkg` (`ls` gives the 15
AR50 counties and `clc2018`), 16 `[source.<id>-clean-25833]` entries in
`~/.config/rasputin/sources.toml` (`grep -c`), and the backup
`sources.toml.bak-20261101`. 115's proof run, Numedalslågen at tolerance
10 with `--land-cover ar50-buskerud ar50-telemark` against the two
products: the two `.vtk` files are byte-identical (`cmp`), 2,046,854
triangles (the file's `POLYGONS` line). What remains before merge is this
branch's own check on Sula and Bygdin with the `-clean-` ids, in the
acceptance step (section 8). His commands now name `clc2018-clean-25833`,
`ar50-innlandet-clean-25833` and so on; one that names the raw id is
refused with the sibling as its fix (section 3); the MCP server reads the
same catalogue. A product made in EPSG:25832 and
meshed on DTM10 (EPSG:25833) passed on Bygdin's box (0 of 1,454), but in
general it can break in the move and is then refused with the second line
form above: make one product per DEM CRS, as 115 says.

**A public user** with a broken file, the shipped `clc2018` source in
Norway among them (534 invalid polygons in the whole cache), gets the
refusal: the file, the count, each fid and reason, and the fix, to repair
the file in the DEM's CRS. Nothing public names the private clean-up.

## 6. Degeneracy policy

- Empty and missing geometries: skipped and counted as today, not checked.
- A polygon with a vertex that has no image in the DEM's CRS: refused by
  `_moved` as today, before the validity check.
- A MultiPolygon is checked whole (`is_valid` on the multipolygon: parts
  overlapping each other are invalid), as the product writes them.
- Near and far: a polygon whose bounding box only touches the read box
  (closed intervals, as `_fixed`'s `<=` and `>=`) is near and checked.

## 7. Tests `@tester` writes red

Fixtures, measured in the probes above (shapely 2.2.0, GEOS 3.14.1):

- the bow-tie `Polygon([(0,0),(4,4),(4,0),(0,4)])`: reason
  `Self-intersection[2 2]`;
- the hole across its shell, shell the square (0,0)-(10,10), hole
  (5,5),(15,5),(15,6),(5,6): reason `Self-intersection[10 5]`;
- AR50 Innlandet fid 191902, 12 vertices, valid in EPSG:25832, invalid in
  EPSG:25833 (reason `Self-intersection[359920.000048706 6825130.00009623]`
  on PROJ as installed; match the prefix `Self-intersection[` only). Its
  ring, in EPSG:25832, is in scratch `move_broken.py`'s output; the test
  carries it inline with NIBIO's NLOD credit added to
  `tests/fixtures/ar50/NOTICE`.

1. `read_class_polygons` with `near`: a file (EPSG:25833) holding a valid
   square and the bow-tie raises `FeatureError` whose text holds the file
   name, "1 of the 2 polygons", "EPSG:25833", the bow-tie's fid and
   `Self-intersection[2 2]`, and not the square's fid.
2. Every one named: the bow-tie and the hole across its shell in one file:
   both lines, in fid order, "2 of the 3"; no first-only stop.
3. Far is not checked: the bow-tie moved 10,000 m east, its box missing
   `near`: no error; `build` puts it in `outside` as the same object.
4. Broken by the move: fid 191902 in a file in EPSG:25832, mesh CRS
   EPSG:25833: refused, the line ends "(valid in the file's CRS,
   EPSG:25832: broken by the move)"; with mesh CRS EPSG:25832: accepted.
5. `read_polygons`: the bow-tie near is refused naming file and fid; the
   bow-tie far is dropped without error; a valid far polygon is returned,
   in file order with the near ones; the return is a list, and `tell` is
   no parameter. The thin spike of `test_cli_mesh_prefer.py`'s test 14
   (invalid, `Self-intersection[506380 8314890.90000001]`, measured) is
   refused as "lake.geojson: 1 of the 1 polygons" with its reason.
6. CLI `--land-cover` with the bow-tie file: exit 2, "Invalid value for
   --land-cover", the fid and reason; the `--out` file does not exist.
7. CLI `--prefer lake=` the bow-tie file: exit 2, "Invalid value for
   --prefer".
8. CLI `--reservoirs` the bow-tie file: exit 2, "Invalid value for
   --reservoirs", and refused before the mesh: the DEM-mesh function the
   CLI calls (`cli._dem_mesh`), patched to fail if reached, is not reached.
9. Gone: `feature_input._repaired`, `cli._say_repaired`,
   `LandCoverMap.repaired`, `build(fix=)` (now `cut=`, keyword-only,
   default True); a valid
   `--land-cover` run's record has no `land_cover_repaired`, a
   `--reservoirs` run's no `reservoirs_repaired`, and `run_record` describes
   neither.
10. Unchanged on valid input: the existing tests that read valid polygons
    stay as they are; the repair tests (`test_read_repair.py`, the repair
    classes of `test_feature_input.py`, `test_land_cover_map.py` and
    `test_cli_mesh_edges.py`) become refusal tests or are deleted, each
    deletion named in the red commit.
11. `clean_sibling`, pure: with `known` holding `bad` and
    `bad-clean-31983`, `clean_sibling("bad", "EPSG:31983", known)` is
    `"bad-clean-31983"`; None for key None, for a `known` holding only
    `bad-clean-25833`, and for a label that is not `EPSG:n` (a PROJ string
    `+proj=tmerc +lon_0=7.3 +ellps=GRS80`, which `crs_label` gives as text);
    `x-clean-25832` on `EPSG:25833` finds none. The tests pass labels, so
    `sources` needs no `crs` import (layering, section 3).
12. CLI, the sibling named (red now: today the run repairs and exits 0).
    The scene of `test_cli_land_cover_polygons.py` (DEM in EPSG:31983,
    domain bounds 505489.1 8315090.9 507609.3 8317210.7) and a user
    `sources.toml` (the pattern of
    `test_a_catalogue_vector_source_uses_its_legend`) holding `bad`, a
    local `classes = "corine"` entry whose GeoPackage holds the bow-tie
    moved by (506500, 8316100) (bounds 506500 8316100 506504 8316104, reason
    `Self-intersection[506502 8316102]`, both measured), and
    `bad-clean-31983`, any valid cover file. `--land-cover bad`: exit 2,
    the output holds "use bad-clean-31983 in place of bad" before
    `Self-intersection[506502 8316102]`. The same run with no
    `bad-clean-31983` entry, and with the GeoPackage given by path: exit 2,
    no "-clean-" in the output. `--reservoirs bad` with the sibling (a
    valid `--land-cover` file beside it): "use bad-clean-31983". The
    sentence is matched with whitespace removed (Rich wraps).
13. Local only (`test_cli_mesh_edges.py`, skipped without Ola's data): Sula
    with `--land-cover clc2018-clean-25833` gives Sol's mesh of 28 October,
    3,271 triangles, 30.4411 km², lake 8.935 % (the figures hold on the
    clean copy, checked by `@tester` on a scratch copy with the check); Sula
    with the raw `clc2018` exits 2, 4 of its 31 polygons refused (section
    5), and the message names `clc2018-clean-25833` first.

No suite is named invariant-critical: no mutation round.

## 8. Acceptance

No C++ and nothing under `include/terrain/refinement/` or
`include/terrain/mesh/` changes, so no `@perf` acceptance run. Before
review, the green branch runs the quick check's two land-cover cases
(Numedalslågen and Lagan, `tools/bench_quick.py`): each `.vtk` byte-identical
to master's, as section 5's zero counts predict. The same step, on Ola's
catalogue: Sula with `--land-cover clc2018-clean-25833` (and the raw
`clc2018`: section 7's test 13 runs both, so `pytest` on Ola's machine
covers Sula) and Bygdin with
`ar50-innlandet-clean-25833` (and `clc2018-clean-25833`) exit 0 (section
5's zero counts for the products); Bygdin with the raw `ar50-innlandet`
exits 2 and its message begins "use ar50-innlandet-clean-25833", the
ruling seen on a real case. This belongs in acceptance, not review: it
runs the green branch on Ola's data, and the reviewer reads its output.
No timing gate: the
check is today's `is_valid` over today's set of polygons, moved from
`_fixed` to the reader with the same grouping and threads; for `--prefer`
and `--reservoirs` it checks only the polygons reaching the box, where
`_repaired` checked every polygon read.

Run by `@reviewer` at 8b458b8e: Numedalslågen and Lagan `.vtk` byte-identical to ae33748e's; Bygdin with raw `ar50-innlandet` and raw `clc2018` exits 2, naming the `-clean-25833` copy first; both clean ids exit 0.

## 9. `fractions-intersection`

`polygon_fractions.subdivide` halves valid polygons with `clip_by_rect`,
which can give invalid halves and lose area (pipe item 1.2). This
increment does not touch it: its input, `map.areas`, is valid before and
after, and the fix (halve with `shapely.intersection`) is independent of
where validity is enforced. It stays its own fast-lane item, next.

## 10. Left out

- `--catchments` units and `catchment --lakes`: polygon inputs neither
  repaired nor checked today; not land cover, so not in Ola's ruling. The
  probe found none invalid in `smhi-svar2022-subcatchments` (26,340) or
  `nve-gauge-catchments` (2,132). An idea, not this item.
- The GeoJSON fid mismatch: the core names a GeoJSON feature by its 0-based
  position, the product's `clean_log` by its 1-based one (115, section
  4.1). Only a raw GeoJSON source sees both; noted for the handback.

## 11. Size

Estimate, counted lines (`tools/count_loc.py`): about 55 removed
(`_repaired` 21, `read_polygons`' repair and `tell` 9, `_fixed`'s check
and repair and the count 10, `_say_repaired` with its calls and the record
keys 15) and about 60 added (`_refuse_invalid` 22, the `near` keyword and
calls 5, the reservoir read moved up 8; the amendment's sibling:
`InvalidPolygons` 2, `clean_sibling` 7, `_invalid` 8, its three call sites
and `known` passed to `_land_cover_map` 8): **net about +5** (an
estimate, before the amendment −20). One PR, far under 700.

## 12. Questions for Ola

None. 115's switch is done (section 5); the sibling rule is his ruling of
2026-11-01 (section 3).

## Review

- **Code review, round 1, 2026-11-01.** `@reviewer`, head `8b458b8e49f21cce66eea234cbff7ba09883ca45`, range `ae33748e..8b458b8e`. Verdict: APPROVED WITH WORDINGS. LOC: 31 net (117 added, 86 removed; estimate about +5; `_refuse_invalid` about 40 lines with the fid sort, the message and the moved-note; the land-cover map's repair removal saved less). The checking box is `_fixed`'s cutting box; only the CLI calls `build` and `read_polygons`; reservoirs read before the mesh change no run that worked; `clean_sibling` stays in layer 0. 557 tests in the touched suites pass, both local Sula tests included. Three wordings applied in this commit: the README's land-cover contract, `project_structure.md`'s `read_class_polygons` entry, and section 8's run record with the Status line. Follow-ups: one helper for the vertex-count grouping `_refuse_invalid` and `_fixed` share; with DEM cells over 50 m an invalid far polygon in the strip between the read region and the `--prefer` box is dropped silently (as on master; outside the domain).
