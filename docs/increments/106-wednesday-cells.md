# Increment 106: four fixes for Wednesday's Swedish cells

**Status: designed** by `@architect`, 2026-10-26, on master `f17b3e11`,
worktree `demo-fixes`. Four small PRs, built in the order 106a, 106b, 106c,
106d. Number 106: 105 is the highest under `docs/increments/` over every ref
(`git for-each-ref` + `git ls-tree`, 2026-10-26).

Ola, 2026-10-26: "I want to use the MCP to build the catchments our Swedish
client wants Wednesday morning." He chose these four fixes the same day. On
Wednesday 2026-10-28 he builds Shyft cells for the client's catchments
through the MCP server, from the Opal outlines
(`~/rasputin-output/opal-se2-2026-10-21/`, read only). The MCP server builds
each tool from the CLI's own options, so every option below is an MCP
parameter with no extra work.

Evidence: the rehearsal note `docs/research/mcp-swedish-cells-rehearsal.md`
(brought into this branch from `worktree-shyft-cells-se`, commit `9d7af892`)
and the gap analysis before it (commit `e4b393b7` on the same branch, not
pushed; its figures that matter are repeated here). The probes in
`docs/increments/106-probes/` were run for this design, each with its command
in its docstring; their outputs went to `/tmp/p106/` (not kept).

## 0. The plain picture

| fix | what goes wrong today | after the fix | net lines | path | `@perf` |
|---|---|---|---|---|---|
| 106a | a self-intersecting land-cover or reservoir polygon crashes the cell table (AR50 Innlandet, Norway) | the polygon is repaired when read, and the run says how many were | about 20 | short | no |
| 106b | `network-catchment --network` takes only a file path, not the catalogue name `smhi-svar2022-subcatchments` | `--network`, `--lakes` and `--outlets` take a catalogue name or a file | about 15 | short | no |
| 106c | every cell has `catchment_id` 1; tagging cells afterwards puts 18.9 % of Kultsjön's area in the wrong one of its 75 SVAR units | `mesh --catchments` repairs SMHI's units, simplifies them with rasputin's own border simplification, makes their borders constraints, so no cell straddles two, and writes each cell's unit id | about 90 | full loop | **yes** |
| 106d | no Swedish reservoir polygons: the reservoir share is 0 | a catalogue entry can keep only the lakes a dam names; `smhi-reservoirs` gives Kultsjön 82.73 km² of reservoir | about 35 | short | no |

About 160 net production lines in all. No C++ in any of them.

`@perf` applies to 106c only: it adds constraint lines to what drives the
mesh (*Acceptance* in `docs/increments/README.md`). 106a changes only the
polygons the land-cover table and the reservoir share read after the mesh is
made; 106b only how a file is found; 106d only which reservoir polygons are
read. None of those three touches the mesh.

Under default flags 106a and 106c add 0 s: the quick check's cases
(`docs/benchmarks/quick/cases.toml`) use no `--land-cover`, no
`--reservoirs` and no `--catchments`, and 106a's two readers have no other
caller (`src_python/tin_engine/cli.py`, `read_polygons` and
`read_class_polygons`).

**Round 2** (2026-10-26). The choice of fixes is Ola's answer of
2026-10-26 to the main session's numbered question: the four fixes in
order, default taken. After round 1's review the main session cut 106c to
what Wednesday needs (each cell's catchment id, the borders as constraints,
the repair of SMHI's units); choosing fine units or main sub-basins per
call, by the file given, is the main session's default for Ola's open
question. So 106c drops the `catchments` table, `--catchments-simplify`,
`--catchments-layer` and grouping by an attribute, and simplifies with
rasputin's own border simplification; 106d's `named_by` keeps only what the
`smhi-reservoirs` entry uses.

## 1. Prior art: legacy and literature

*Literature.* Polygon repair: the OGC Simple Features validity rules, and
GEOS's `MakeValid` ("linework" method, shapely's `make_valid` default), which
rebuilds a ring that crosses itself into the polygons its line work encloses.
Coverage repair: GEOS's `CoverageCleaner` (shapely 2.2, `coverage_clean`),
used as published. Simplification: rasputin's own border simplification
(`src_python/tin_engine/border_simplify.py`, design
`docs/increments/32-landcover-simplify.md`, which carries its literature),
each border moved at most a band, each polygon keeping its area, junctions
fixed; round 1's GEOS `coverage_simplify` is no longer used. The cell
format: Shyft's `geo_personality`
(`cpp/shyft/hydrology/cell/geo_personality.h`, gitlab.com/shyft-os/shyft,
read for increment 105), whose `catchment_id` groups cells into
sub-catchments. Nothing here is new and no novelty is claimed.

*Legacy.* `git grep -n -i -E "catchment_id|make_valid|buffer\(0\)|reservoir"
legacy-archive -- legacy` returns nothing; `git grep -n -i -E
"catchment|valid" legacy-archive -- legacy` returns two unrelated lines
(`legacy/rasputin/mesh.py:125`, `legacy/rasputin/web_visualize.py:64`).
Nothing to carry across.

*Ideas file.* 106c takes up "Neighbouring catchments simplified together"
(shared borders stay shared). It serves "Swedish SE2 river-system
catchments" (the client's sub-catchments). 106a overlaps "Repair CORINE once
per area of interest" only in name: that entry caches the land-cover clean-up
for speed; 106a repairs invalid polygons on the cell-table path, which that
clean-up never sees.

## 2. 106a: repair invalid polygons when read

**The fault.** `polygon_fractions.assign` and `cells.overlay` intersect each
triangle with the polygons through GEOS, which throws a `TopologyException`
on a polygon whose ring crosses itself. Probe `ar50_invalid.py`: of the 2,843
AR50 Innlandet polygons read for the Nedre Sjodalsvatn catchment (NVE gauge
2.13.0, 480.61 km²), **32 are invalid**, all "Self-intersection"; cutting one
triangle against each throws on one of them (class 60, 19,503 m²).
`make_valid` turns each into 1 to 4 polygons with the same area to 0.1 m²,
and drops no polygonal part. The gap analysis repaired them by hand and the
run then went through (5,940 cells).

**Behaviour.** `read_class_polygons` and `read_polygons`
(`src_python/tin_engine/feature_input.py`) repair each polygon after moving
it into the mesh's CRS (a move can make a ring cross itself too): a valid
polygon is kept as it is; an invalid one is replaced by the polygonal parts
of `shapely.make_valid(g)` (default method), and dropped when it has none.
Both functions also return how many they repaired or dropped. Their only
callers are `mesh`'s land-cover and reservoir steps
(`src_python/tin_engine/cli.py`).

**Output.** When the count is above 0, one stderr line, because the data was
found wrong (Ola's rule: warnings only from findings), e.g.
`land cover: 1 of 2843 polygons crossed itself and was repaired`, and the same
with `reservoirs:`. The run record gets `land_cover_repaired` and
`reservoirs_repaired` (counts) whenever polygons are read. No refusal.

**Tests** (fixture: the bowtie `Polygon([(0,0),(4,4),(4,0),(0,4)])`, class 1,
in a GeoJSON file in a projected CRS; the triangle `(0,0),(4,0),(0,4)`, plan
area 8):
1. Today `assign` on the bowtie throws `GEOSException` ("TopologyException:
   side location conflict at 2 2"), and so does `overlay` (measured): the red
   step.
2. `read_class_polygons` returns the bowtie as two triangles, valid, area 8.0,
   and the count 1; `assign` on its output gives class 1 a fraction of 0.5
   and class 0 (uncovered) 0.5.
3. `read_polygons` and `overlay` give the triangle a reservoir share of 0.5.
4. A valid polygon comes back unchanged (`equals_exact`, tolerance 0) and
   counts 0.
5. `Polygon([(0,0),(1,1),(2,2),(0,0)])` (no area; `make_valid` gives a
   multilinestring) is dropped and counts 1.
6. A `mesh` run with the bowtie as `--land-cover` completes, prints the
   stderr line and records `land_cover_repaired: 1`.

**Size.** About 20 net lines: the repair helper (8), the two readers (4), the
two callers' lines and record entries (8). Short path.

## 3. 106b: catalogue names for `network-catchment`

**Behaviour.** `--network`, `--lakes` and `--outlets` of `network-catchment`
take a catalogue source of features (a key from the sources file, such as
`smhi-svar2022-subcatchments`, `smhi-lakes`, `smhi-outlet-points`) or a
`.gpkg`, `.geojson` or `.json` file, exactly as `mesh --reservoirs` does today.
A catalogue source's layer is its `layer`, unless `--network-layer` (or
`--lakes-layer`, `--outlets-layer`) is given. The resolution is one helper,
taken out of `_reservoir_source` (`src_python/tin_engine/cli.py`), which
`--reservoirs`, these three options and 106c's `--catchments` call. The
options become `str` in place of `Path`.

**Refusal.** A value that is neither: `X is neither a catalogue source of
features nor a .gpkg, .geojson, .json file` (the `--reservoirs` wording),
headed by the option. A source that is not cached names the fetch to run, as
`--reservoirs` does.

**Tests.** With a temporary sources file (`RASPUTIN_DATA`) registering the
network fixture of `tests/python/test_cli_network_catchment.py`, under a
key:
1. `--network <key>` writes the same output bytes as `--network <path>`;
   the same for `--lakes` and `--outlets`.
2. The layer: a key whose entry has `layer` reads that layer; an explicit
   `--network-layer` wins.
3. `--network no-such-thing` is refused with the wording above, headed
   `--network`.

**Size.** About 15 net lines. Short path.

## 4. 106c: sub-catchments as constraints, and each cell's id

### What it does

`mesh --domain D --catchments UNITS [--features-tolerance METRES]`. Each
polygon of `UNITS` is one sub-catchment. Fine units or main sub-basins are
chosen per call by the file given (section 6); there is no grouping option.

1. **Select.** The units are the polygons of `UNITS` (a catalogue source of
   features or a file, through 106b's helper; a catalogue entry's layer is
   its `layer`, a file's the reader's default) whose inner point
   (`point_on_surface`) lies inside the domain `D`, moved into the units'
   CRS vertex by vertex (`DomainPolygon.to_crs`). The units are read with the
   domain's box, so a GeoPackage's R-tree is used.
2. **Repair, in the units' CRS.** `coverage_clean(snapping_distance=1,
   gap_width=1)` (the land-cover repair's 1 m): SMHI's units as published
   are not an exact coverage (*Measured*).
3. **Simplify with rasputin's own border simplification.** The units and a
   frame (their union's bounding box, grown by 1 km, less the union) go to
   `border_simplify.simplify_borders` as one coverage, band `METRES`,
   clearance 1 m (the repair distance, as for land cover). The frame makes
   the union's outline a border like the others: it is simplified with
   them, every unit keeps its area, and no gap can open. `METRES` is the
   existing `--features-tolerance` (its help gains "and --catchments
   borders"); with `--catchments` and the option not given it is 200 m;
   0 repairs only.
4. **The domain is their union** (`coverage_union_all`), which replaces `D`
   for the rest of the run, before the DEM is opened, so the DEM window is
   chosen from it. `D` only selects.
5. **The borders are constraints.** All units' rings, moved into the DEM's
   CRS by the same vertex-by-vertex transform as the domain (so the
   outline's vertices are the domain's, exactly), become one
   `TerrainFeature` with a new edge property `catchment`, added to the run's
   feature set (with or without `--features`). The property joins the run's
   vocabulary only when `--catchments` is given, so a run without it writes
   the same bytes as today.
6. **Each cell's id.** Each triangle takes the unit holding its centroid
   (an STRtree `within` query); a centroid in no unit (not expected: no
   triangle crosses a border) takes the nearest unit, and a count above 0 is
   a stderr line, `catchments: N cells had their centroid in no unit and
   took the nearest`. `catchment_id` is the unit's fid in `UNITS` as the
   reader gives it (`feature_input.read_source`): a GeoPackage's `fid`
   (SVAR 2022: 1 to 26,340), a GeoJSON feature's `id`, or, with no `id`,
   its position from 0; so a GeoJSON made for `--catchments` gives each
   feature an `id`.
7. **Written.** The `cells` layer's `catchment_id` column carries the ids
   (today the constant 1). Without `--catchments` the cells layer is as
   today. No other table.

Stderr, one line: `catchments: 75 units of smhi-svar2022-subcatchments,
borders simplified at 200 m: 1558 vertices, 1702.67 km2`. The run record:
`catchments` (the source's name), `catchments_units` (the count) and
`catchments_band_m` (the band used, given or 200). The start mesh's minimum
angle stays 25 degrees: the 15-degree rule (`cli.py`, "15 degrees with land
cover simplified") reads only land cover, and the probes below ran at 25.

This also answers the rehearsal's item 1 for cells: the outline is the
simplified union of the units, so `network-catchment`'s unsimplified outline
(19,556 vertices for Kultsjön) is no longer the floor.

### Refusals, each headed by its option

1. `--catchments` without `--domain`: "needs --domain: the units whose inner
   point lies in the domain are the sub-catchments".
2. Without a cells layer: "needs --land-cover with a legend and a .gpkg
   output" (the `--reservoirs` wording; the ids are written only there).
3. Units in a geographic CRS: "UNITS is in degrees (EPSG:4326); the border
   band is in metres, so give the units in a projected CRS".
4. No unit's inner point in the domain: "no unit of UNITS has its inner
   point inside the domain".
5. A selected unit whose fid is not an integer from 1 to 2,147,483,647:
   "UNITS: feature N has no id Shyft's catchment_id can hold".
6. The units' union in more than one part: "the units join into N parts that
   do not touch: a, b km2" (`network.join`'s wording).

`--features-tolerance`'s existing check (finite, at least 0) covers the band.

### Measured (probes, master `f17b3e11`'s package)

Kultsjön and Ljungan, Opal outline as `D`, SVAR 2022 units, `--dem glo30
--out-crs EPSG:3006`; the units simplified as 106c will
(`units_borders.py`) then passed as `--features`, their union as
`--domain`, which is what 106c builds in memory; `straddle.py` counts.

- Selection: 428 units near the Opal outline, **75** with the inner point
  inside it, the same set `network-catchment` joins. Ljungan: 2,536 near,
  **809** inside, the same set as its join.
- SMHI's units as published are **not an exact coverage**: 74 of Kultsjön's
  75 units have edges that do not match their neighbour's (10,851 m of
  them; Ljungan 758 of 809, 63,590 m). `coverage_clean` at 1 m makes them
  one with no area change (0.0 m²).
- Round 2's first try, rejected: keep `D` as the domain and clip the units
  to it (`opal_gaps.py`). The Opal outline and the units' union differ by
  1.74 km² inside `D` in 2,750 slivers (Kultsjön; Ljungan 16.24 km², 9,082),
  and the land-cover outline rule, under 100 m, leaves 12 of Ljungan's.
  Hence the domain is the units' union, as in round 1.

| units | band | vertices | union's outline | unit area change | fewest triangles | simplify |
|---|---|---|---|---|---|---|
| Kultsjön, 75 | 0 m | 142,628 | 19,556 | 0 | not run | 0 |
| | 100 m | 2,724 | 390 | under 3e-12 | 5,694 | 0.5 s |
| | **200 m** | 1,558 | 196 | under 3e-12 | **3,498** | 0.6 s |
| | 300 m | 1,249 | 143 | under 3e-12 | 3,121 | 0.7 s |
| Ljungan, 809 | 0 m | 1,411,364 | 111,249 | 0 | not run | 0 |
| | **200 m** | 16,855 | 1,110 | under 4e-12 | **36,595** | 6.6 s |
| | 300 m | 13,253 | 758 | under 4e-12 | 31,851 | 7.9 s |

Every row: coverage valid, the union one polygon with no hole, its area
unchanged. Round 1's GEOS `coverage_simplify` at 300 m changed one unit by
36 % (Kultsjön) and 75 % (Ljungan) and left 2,002 vertices and 3,717 fewest
triangles for Kultsjön; rasputin's simplification keeps every area and
needs fewer. Reading, selecting and repairing take 0.3 s (Kultsjön) and
3.0 s (Ljungan); the frame's box comes from the union's bounds (a 1 km
buffer of Ljungan's 111,249-vertex outline took minutes).

- **No cell straddles.** Kultsjön at 200 m, `--triangles 4000`: 3,972
  triangles, 1.3 s, none with more than 4.2e-5 of its area outside the unit
  holding most of it, every centroid in a unit, 5 to 466 cells per unit. At
  300 m: 3,964 triangles, at most 2.0e-4 outside. Ljungan at 200 m,
  `--triangles 40000`: 40,175 triangles, 7.0 s, at most 3.0e-4 outside, 2
  to 446 cells per unit. The residue is the mesh's edges being straight in
  the DEM's degrees and the units' straight in metres. Control, the same
  check able to fail: Kultsjön's 200 m union meshed **without** the
  borders, `--triangles 4000`: 3,935 triangles, **1,383** with more than
  1e-4 outside, some with no unit holding half of them.
- Main sub-basins (`baro_groups.py`, round 1): by `BARONR`, Kultsjön's 75
  units are 1 sub-basin (1,702.7 km²); Ljungan's 809 are 4 (1,412.7 to
  4,315.2 km²), each one polygon. With all 809 units Ljungan cannot go
  under 36,595 cells at 200 m; the rehearsal's Ljungan run without
  sub-catchments had 9,985.

**Constants.** The 1 m repair distance assumes a projected CRS in metres and
SVAR's mismatches; checked on Ljungan's 809 units (the largest input), where
it changes no area. The 1 km frame margin only has to exceed the band
(checked at 300 m). The 200 m default band assumes units of 0.004 to
1,000 km² and a few thousand cells per 1,000 km²: it is the smallest band
probed that lets Kultsjön take the rehearsal's 4,000; checked on Kultsjön
and Ljungan (above).

### Tests (key properties)

Figures measured with shapely 2.2.0 (GEOS 3.14.1).

1. **Select** (pure): three units, the squares x in [0, 1000], [1000, 2000],
   [2000, 3000] m, y in [0, 1000]; domain `box(0, 0, 2050, 1000)`: the first
   two are selected (the third's inner point (2500, 500) is outside). And a
   U-shaped unit, `box(0, 0, 3000, 1000)` less `box(1000, 300, 2000, 1000)`
   (2,300,000 m²), whose centroid (1500, 454.35) lies in its notch, outside
   it, and whose inner point is (500, 650): against `box(900, 0, 2100,
   1000)` (holds the centroid, not the inner point, and intersects the unit)
   it is not selected; against `box(0, 0, 1000, 1000)` (holds the inner
   point, not the centroid) it is. A centroid or an intersects rule fails one
   of the two.
2. **Repair** (pure): `A = Polygon([(0,0),(1000,0),(1000,500),(1000,1000),(0,1000)])`
   and `B = Polygon([(1000,0),(2000,0),(2000,1000),(1000,1000),(1000.1,500)])`
   (a 50 m² slit): as given their coverage union is in two parts with
   1,999,950 m²; after the step it is one polygon of 2,000,000 m², B
   1,000,000 m² (measured, round 1).
3. **Ids** (pure): units with fids 7 and 9 get `catchment_id` 7 and 9; a
   selected unit with fid 0, 2^31 or `"x"` is refused (refusal 5).
4. **Simplify keeps the coverage and the outline moves** (pure): three units
   tiling about [0, 3000] x [0, 1000] m, the two shared borders zigzags at
   x = 1000 and 2000 (20 m out every other 50 m step) and the top edge a
   zigzag (20 m up every other 50 m step); areas 1,020,000, 1,010,000 and
   1,000,000 m², 148 vertices, the union's outline 66. At band 50 m: each
   unit's area unchanged (relative 0.0), coverage valid, the union one
   polygon with no hole, 110 vertices, the union's outline **28**, and that
   union is the domain written to the mesh. Without the frame the outline
   stays 66 (the "frame left out" mutant).
5. **End to end** (synthetic DEM, `.gpkg` with a class-raster land cover,
   projected CRS, domain `box(0, 0, 3000, 3000)`), three units with two
   borders on no grid line or grid diagonal:
   `A = Polygon([(0,0),(3000,0),(3000,700),(0,1900)])` (3,900,000 m²),
   `B = Polygon([(0,1900),(3000,700),(3000,3000),(1700,3000)])`
   (4,165,000 m²), `C = Polygon([(0,1900),(1700,3000),(0,3000)])`
   (935,000 m²); a valid coverage whose union is the square; inner points
   (750, 1300), (2250, 1300), (425, 2450). With `--features-tolerance 0` so
   the borders stay as given: no triangle has a share above 1e-6 in two
   units; each cell's `catchment_id` is the fid of the unit holding its
   centroid, and every fid appears; stderr has **no** "centroid in no unit"
   line; the run record entries.
6. **Unchanged without it**: the same run without `--catchments` writes
   `catchment_id` 1 on every row and the same `triangles` and `cells` rows
   as master.
7. **Refusals** 1, 2, 3, 4 and 6 above.

**Invariant-critical suite** (mutation round required): tests 1, 4, 5 and 6.
Mutation targets: the inner-point rule (centroid or intersects in its place:
test 1's U), the frame left out (test 4), the borders not added as
constraints (test 5's slanted borders), the ids taken from the first
triangle vertex in place of the centroid (test 5's absent stderr line: a
vertex on a border is within no unit), the fid range check (test 3).

**`@perf`**: the 1 m benchmark (no `--catchments`) unchanged, meshes
byte-identical to master's; Kultsjön at `--triangles 4000` and Ljungan at
`--triangles 40000` with `--catchments smhi-svar2022-subcatchments`, the
predicted counts (3,972 and 40,175, within the budget's 5 %), at most 3e-4
of any triangle outside its unit, and the time against the probes: mesh
1.3 s and 7.0 s, plus reading and repairing 0.3 s and 3.0 s, plus the
simplification 0.6 s and 6.6 s, plus 10 % noise.

**Size.** About 90 net lines: the option and its checks (10), a pure module
`sub_catchments.py` (select, repair, the frame and `simplify_borders`, the
union, the borders as a feature, the ids: 40), the wiring in `mesh` (20),
the writer's `catchment_id` column (8), stderr and record (12). Full loop.

## 5. 106d: Swedish reservoirs as a catalogue entry

**The data** (probe `se_reservoirs.py`): SMHI's dam register (10,702 dams)
names the water body a dam holds in `vy_eucd` on 1,632 dams (1,356 distinct);
`vy_eucd` is `SE` followed by the lake's `SJOID` in SMHI's lakes (42,178
polygons). 1,317 lakes are named so (by `VYID` only 20). No match without
removing the prefix. Inside Kultsjön: Satsjön 1.07, Ransarn 28.26 and
Kultsjön 53.40 km², **82.73 km²**. `regl_typ` is no filter (Kultsjön's dam
has 0).

**Behaviour.** A local features entry in a sources file may say which of its
features to keep:

```toml
[source.smhi-reservoirs]
kind = "local"
path = ".../vattenytor_2016_3006.geojson"
named_by = { source = "smhi-dams", field = "vy_eucd", key = "SJOID", prefix = "SE" }
credit = "Källa: SMHI, SVAR 2016 (vattenytor) and dammregistret, hämtade 2026-10-21."
licence_note = "CC BY 4.0 (SMHI's terms of use)."
description = "SMHI's lakes that a dam in SMHI's dam register names: reservoirs for mesh --reservoirs."
```

A feature is kept when `prefix + str(its key)` equals the `field` of some
feature of the source `source`. `named_by` is only as wide as this entry:
its four fields are all required text, with no defaults and no other match
rule, and both the entry and `source` are local `.geojson` features entries
(no layer, no GeoPackage, no remote source). `--reservoirs
smhi-reservoirs` reads the lakes (as today, with the run's box) and keeps the
named ones. The join is data, so the MCP call is `reservoirs:
"smhi-reservoirs"`.

**Refusals.** `named_by` on an entry that is not a local `.geojson` source
of features, or with a field missing (a sources-file error, as the other field checks; the `RemoteSource` model in
`src_python/tin_engine/sources.py`); `named_by.source` not in the catalogue
or not a local `.geojson` of features: `--reservoirs` refused naming it; any other option given an
entry with `named_by` (106b's helper): "X keeps only the features another
source names; only --reservoirs reads it".

**Tests.**
1. The model accepts `named_by` on a local `.geojson` entry and refuses it on
   a heights source, an `arcgis-features` one, a local `.gpkg` one, and with
   `prefix` left out.
2. The filter (pure): keys `1-1`, `2-2`, `3-3`; names `SE1-1`, `SE3-3`,
   none, `NO9`; prefix `SE`: keeps `1-1` and `3-3`.
3. `mesh --reservoirs <entry>` with a temporary sources file over two lake
   squares in the domain, one named by a dam point: the reservoir area is the
   named lake's (relative 1e-9); the other lake adds none.
4. `network-catchment --network <entry with named_by>` is refused.

**Registration.** Ola's sources file is outside git; after the merge the
main session adds the entry above to it (the path as for `smhi-lakes`), and
runs `--reservoirs smhi-reservoirs` on Kultsjön once.

**Size.** About 35 net lines: the model and its check (10), the filter and
its read (17), the helper's refusal and the record (8). Short path.

## 6. What Wednesday's call looks like

Per Opal outline, through the MCP `mesh` tool (CLI spelling), SVAR's fine
units:

```
rasputin mesh --out kultsjon.gpkg --dem glo30 --out-crs EPSG:3006
  --domain ~/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon.geojson
  --catchments smhi-svar2022-subcatchments
  --land-cover clc2018 --reservoirs smhi-reservoirs --triangles 4000
```

Each cell's `catchment_id` is then its unit's `fid` in SMHI's GeoPackage,
whose row also carries `ARO_UUID`, `BARONR` and the name, so the client can
join them. For SMHI's main sub-basins Ola gives a different file to
`--catchments`: a GeoJSON of SVAR's units joined by `BARONR`, each feature's
`id` its `BARONR` (1 to 11,700,118, within Shyft's range). The main session
writes it once with a short script after 106c merges, a stopgap outside the
product; units with no `BARONR` are left out of it.

Left out, as in the gap analysis: NMD land cover, glacier outlines, routing,
a Shyft read-back example, `--outline-tolerance` on `network-catchment`
(106c makes it unnecessary for cells).

## 7. Questions for Ola (defaults taken unless he says otherwise)

1. The border band's default, 200 m: every unit keeps its area exactly, but
   a border may move up to 200 m from SMHI's line; Kultsjön can then go down
   to 3,498 cells, Ljungan with all 809 units to 36,595. Default: 200 m,
   overridden per call with `--features-tolerance`.

## Review record

- **Design review, round 1, 2026-10-26.** `@reviewer`, head `89e22445b28122ba6b1915535a632cc1c48a3a73`. Verdict: CHANGES REQUESTED. Blocking: (1) 106c test 1's square units let a centroid rule pass: add a non-convex unit; (2) 106c test 5 lets "ids from the first vertex" and "borders not added" pass: assert the nearest-unit stderr line absent and slant a border; (3) scope: cut 106c to Wednesday's needs, `named_by` to the entry's. Not pushed; no CI.
