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
| 106c | every cell has `catchment_id` 1; tagging cells afterwards puts 18.9 % of Kultsjön's area in the wrong one of its 75 SVAR units | `mesh --catchments` makes the sub-catchment borders constraints, so no cell straddles two, and writes each cell's sub-catchment id | about 130 | full loop | **yes** |
| 106d | no Swedish reservoir polygons: the reservoir share is 0 | a catalogue entry can keep only the lakes a dam names; `smhi-reservoirs` gives Kultsjön 82.73 km² of reservoir | about 45 | short | no |

About 210 net production lines in all. No C++ in any of them.

`@perf` applies to 106c only: it adds constraint lines to what drives the
mesh (*Acceptance* in `docs/increments/README.md`). 106a changes only the
polygons the land-cover table and the reservoir share read after the mesh is
made; 106b only how a file is found; 106d only which reservoir polygons are
read. None of those three touches the mesh.

## 1. Prior art: legacy and literature

*Literature.* Polygon repair: the OGC Simple Features validity rules, and
GEOS's `MakeValid` ("linework" method, shapely's `make_valid` default), which
rebuilds a ring that crosses itself into the polygons its line work encloses.
Coverage repair and simplification: GEOS's `CoverageCleaner` (shapely 2.2,
`coverage_clean`) and `CoverageSimplifier` (shapely 2.1+,
`coverage_simplify`), the latter Visvalingam and Whyatt's area-based line
simplification (Visvalingam, M. and Whyatt, J. D. (1993), *Line
generalisation by repeated elimination of points*, The Cartographic Journal
30(1), 46-51) applied to each shared edge once, so neighbours stay gap-free;
both are used as published. The cell format: Shyft's `geo_personality`
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

`mesh --domain D --catchments UNITS --catchment-attribute NAME
[--catchments-layer L] [--catchments-simplify METRES]`:

1. **Select.** The units are the polygons of `UNITS` (a catalogue source of
   features or a file, through 106b's helper) whose inner point
   (`point_on_surface`) lies inside the domain `D`, moved into the units'
   CRS vertex by vertex (`DomainPolygon.to_crs`). The units are read with the
   domain's box, so a GeoPackage's R-tree is used.
2. **Group.** Units with the same value of `NAME` are one sub-catchment.
   With `ARO_UUID` every SVAR unit is one; with `BARONR` (SMHI's main
   sub-basin number) units are joined into sub-basins.
3. **Clean and simplify as one coverage, in the units' CRS.**
   `coverage_clean(snapping_distance=1, gap_width=1)` (the land-cover repair's
   1 m), the groups joined (`coverage_union_all` per group), then
   `coverage_simplify(tolerance)`, `METRES` defaulting to 300; 0 cleans only.
4. **The domain is their union** (`coverage_union_all`), which replaces `D`
   for the rest of the run, before the DEM is opened, so the DEM window is
   chosen from it. `D` only selects.
5. **The borders are constraints.** Each group's rings, moved into the DEM's
   CRS by the same vertex-by-vertex transform as the domain (so the outline's
   vertices are the domain's, exactly), become one `TerrainFeature` with a new
   edge property `catchment`, added to the run's feature set (with or without
   `--features`). The property joins the run's vocabulary only when
   `--catchments` is given, so a run without it writes the same bytes as
   today.
6. **Each cell's id.** Each triangle takes the group holding its centroid
   (an STRtree `within` query); a centroid in no group (not expected: no
   triangle crosses a border) takes the nearest group, and a count above 0 is
   a stderr line. `catchment_id` is the group's value when every value is an
   integer in 1 .. 2^31 - 1 (`BARONR`, `HARO`), else 1 .. K in the values'
   text order (`ARO_UUID`).
7. **Written.** The `cells` layer's `catchment_id` column carries the ids
   (today the constant 1), and a new attributes table `catchments` holds one
   row per group: `catchment_id`, `value` (the attribute's text), `cells`,
   `area_m2` (the sum of its cells' plan areas) and `source_area_m2` (the
   group's area as read, cleaned and unsimplified, in the units' CRS). Without
   `--catchments` the cells layer is as today and there is no `catchments`
   table.

Stderr, one line: `catchments: 75 units of smhi-svar2022-subcatchments in 75
groups by ARO_UUID, simplified at 300 m: 2002 vertices, 1702.13 km2 (1702.67
as published)`. The run record: `catchments` (the source's name),
`catchment_attribute`, `catchments_simplify_m`, `catchments_groups`.

This also answers the rehearsal's item 1 for cells: the outline is the
simplified union of the units, so `network-catchment`'s unsimplified outline
(19,556 vertices for Kultsjön) is no longer the floor.

### Refusals, each headed by its option

1. `--catchments` without `--domain`: "needs --domain: the units whose inner
   point lies in the domain are the sub-catchments".
2. Without a cells layer: "needs --land-cover with a legend and a .gpkg
   output" (the `--reservoirs` wording; the ids are written only there).
3. `--catchments` without `--catchment-attribute`, or either of
   `--catchment-attribute`, `--catchments-layer`, `--catchments-simplify`
   without `--catchments`: "needs --catchments" / "needs
   --catchment-attribute".
4. Units in a geographic CRS: "UNITS is in degrees (EPSG:4326);
   --catchments-simplify is in metres, so give the units in a projected CRS".
5. No unit's inner point in the domain: "no unit of UNITS has its inner
   point inside the domain".
6. A unit with no value for `NAME`: "UNITS: feature N has no NAME".
7. The units' union in more than one part: "the units join into N parts that
   do not touch: a, b km2" (`network.join`'s wording).
8. `--catchments-simplify` negative or not finite: "must be finite and at
   least 0".

### Measured (probes, master `f17b3e11`'s package)

Kultsjön, Opal outline as `--domain`, SVAR 2022 units, `--dem glo30
--out-crs EPSG:3006 --land-cover clc2018`; the units passed as `--features`
and their union as `--domain`, which is what 106c builds in memory
(`units_simplify.py`, `straddle.py`):

- Selection: 428 units near the Opal outline, **75** with the inner point
  inside it, the same set `network-catchment` joins. Ljungan: 2,536 near,
  **809** inside, the same set as its join.
- SMHI's units as published are **not an exact coverage**: 74 of Kultsjön's
  75 units have edges that do not match their neighbour's (10,851 m of
  them; Ljungan 758 of 809, 63,590 m), so `coverage_simplify` on the raw
  units gives a union in pieces. `coverage_clean` makes it valid with no
  area change (0.0 m²); then every tolerance below gives one polygon, no hole.

| units | simplify | vertices | largest unit area change | units over 5 % | total area change | fewest triangles |
|---|---|---|---|---|---|---|
| Kultsjön, 75 (smallest 0.011 km²) | 0 m | 142,650 | 0 | 0 | 0 | not run |
| | 100 m | 5,851 | 35.85 % | 1 | −0.049 km² | 12,559 |
| | 200 m | 2,939 | 35.85 % | 3 | −0.274 km² | 5,714 |
| | 300 m | 2,002 | 35.85 % | 4 | −0.546 km² (0.03 %) | 3,717 |
| Ljungan, 809 (smallest 0.004 km²) | 300 m | 19,568 | 75.31 % | 64 | −1.632 km² (0.01 %) | 32,802 |

Median unit area change at 300 m: 0.71 % (Kultsjön), 0.90 % (Ljungan). The
largest relative changes are on the smallest units.

- **No cell straddles.** Kultsjön at 300 m, `--triangles 4000`: 3,916
  triangles, **0** straddling two units (a triangle's share in a second unit
  above 1e-6), every centroid in a unit, 5 to 480 cells per unit, 11.6 s.
  At 200 m, `--triangles 6000`: 5,714 triangles, 0 straddling. Ljungan at
  300 m, `--triangles 40000`: 40,171 triangles, 0 straddling, 1 to 641 cells
  per unit, 22.4 s. Control, the same check able to fail: Kultsjön's
  300 m union meshed **without** the borders, `--triangles 4000`: 3,888
  triangles, **1,346** straddling.
- The price of the borders at a given count: Kultsjön at about 3,900
  triangles, tolerance found 120 m with borders against 58.5 m without
  (root-mean-square 27.67 m against 16.82 m).
- Grouping (`baro_groups.py`): by `BARONR`, Kultsjön's 75 units are 1 group
  (1,702.7 km²); Ljungan's 809 are 4 (1,412.7 to 4,315.2 km²), each one
  polygon. With all 809 units Ljungan cannot go under 32,802 cells; the
  rehearsal's Ljungan run without sub-catchments had 9,985.

**Constants.** The 1 m clean distance assumes a projected CRS in metres and
SVAR's mismatches; checked on Ljungan's 809 units (the largest input), where
it changes no area. The 300 m default assumes 10 to 1,000 km² units and a few
thousand cells per 1,000 km²; checked on Kultsjön and Ljungan (above).

### Tests (key properties)

1. **Select** (pure): three units, the squares x in [0, 1000], [1000, 2000],
   [2000, 3000] m, y in [0, 1000]; domain `box(0, 0, 2050, 1000)`: the first
   two are selected (the third's inner point (2500, 500) is outside).
2. **Clean** (pure): `A = Polygon([(0,0),(1000,0),(1000,500),(1000,1000),(0,1000)])`
   and `B = Polygon([(1000,0),(2000,0),(2000,1000),(1000,1000),(1000.1,500)])`
   (a 50 m² slit): as given their coverage union is in two parts with
   1,999,950 m²; after the step it is one polygon of 2,000,000 m², B
   1,000,000 m² (measured).
3. **Group and ids** (pure): values 7, 7, 9 give 2 groups with ids 7 and 9;
   values "y", "x" give ids 2 and 1 (text order); a value 0 or 2^31 makes
   the text rule apply.
4. **Simplify keeps the coverage** (pure): on a fixture of three units
   with wiggly shared borders, after `--catchments-simplify` the groups'
   union equals the domain written to the mesh, and the groups still tile it
   (coverage valid, no gap, areas summing to the union's within 1e-9).
5. **End to end** (synthetic DEM, three units in a projected CRS, `.gpkg`
   with a class-raster land cover): no triangle has a share above 1e-6 in two
   units; each cell's `catchment_id` is the unit holding its centroid; the
   `catchments` table has 3 rows whose `cells` sum to the triangle count and
   whose `area_m2` sum to the cells' total; the run record entries.
6. **Unchanged without it**: the same run without `--catchments` writes
   `catchment_id` 1 on every row, no `catchments` table, and the same
   `triangles` and `cells` rows as master.
7. **Refusals** 1, 2, 4, 5 and 7 above.

**Invariant-critical suite** (mutation round required): tests 1, 5 and 6.
Mutation targets: the inner-point rule (centroid or intersects in its place),
the borders not added as constraints, the ids taken from the first triangle
vertex in place of the centroid, the integer-value rule.

**`@perf`**: the 1 m benchmark (no `--catchments`) unchanged, meshes
byte-identical to master's; Kultsjön and Ljungan with `--catchments
smhi-svar2022-subcatchments --catchment-attribute ARO_UUID` at the budgets
above, the predicted counts, 0 straddling, and the time against the probe's
11.6 s and 22.4 s plus the selection and simplification (the probe script
took 10.3 s on Ljungan, reading 1.4 M unit vertices from the GeoPackage and
simplifying at five tolerances).

**Size.** About 130 net lines: three options and their checks (25), a pure
module `sub_catchments.py` (select, clean-group-simplify, union, the borders
as features, the ids: 50), the wiring in `mesh` (20), the writer's
`catchment_id` column and `catchments` table (20), stderr and record (15).
Full loop.

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
feature of the source `source` (also a local features entry). `--reservoirs
smhi-reservoirs` reads the lakes (as today, with the run's box) and keeps the
named ones. The join is data, so the MCP call is `reservoirs:
"smhi-reservoirs"`.

**Refusals.** `named_by` on an entry that is not a local source of features
(a sources-file error, as the other field checks; the `RemoteSource` model in
`src_python/tin_engine/sources.py`); `named_by.source` not in the catalogue
or not features: `--reservoirs` refused naming it; any other option given an
entry with `named_by` (106b's helper): "X keeps only the features another
source names; only --reservoirs reads it".

**Tests.**
1. The model accepts `named_by` on a local `.geojson` entry and refuses it on
   a heights source and an `arcgis-features` one.
2. The filter (pure): keys `1-1`, `2-2`, `3-3`; names `SE1-1`, `SE3-3`,
   none, `NO9`; prefix `SE`: keeps `1-1` and `3-3`.
3. `mesh --reservoirs <entry>` with a temporary sources file over two lake
   squares in the domain, one named by a dam point: the reservoir area is the
   named lake's (relative 1e-9); the other lake adds none.
4. `network-catchment --network <entry with named_by>` is refused.

**Registration.** Ola's sources file is outside git; after the merge the
main session adds the entry above to it (the path as for `smhi-lakes`), and
runs `--reservoirs smhi-reservoirs` on Kultsjön once.

**Size.** About 45 net lines: the model and its check (15), the filter and
its read (20), the helper's refusal and the record (10). Short path.

## 6. What Wednesday's call looks like

Per Opal outline, through the MCP `mesh` tool (CLI spelling):

```
rasputin mesh --out kultsjon.gpkg --dem glo30 --out-crs EPSG:3006
  --domain ~/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon.geojson
  --catchments smhi-svar2022-subcatchments --catchment-attribute ARO_UUID
  --land-cover clc2018 --reservoirs smhi-reservoirs --triangles 4000
```

Left out, as in the gap analysis: NMD land cover, glacier outlines, routing,
a Shyft read-back example, `--outline-tolerance` on `network-catchment`
(106c makes it unnecessary for cells).

## 7. Questions for Ola (defaults taken unless he says otherwise)

1. Which units are the sub-catchments on Wednesday: every SVAR unit
   (`ARO_UUID`: Ljungan 809, at least 32,802 cells) or SMHI's main
   sub-basins (`BARONR`: Ljungan 4)? 106c serves both; default: Ola picks
   per call, nothing to build.
2. The simplification default, 300 m: Kultsjön's total area changes by
   0.03 %, but its smallest unit (0.011 km²) by 36 %. Default: 300 m, with
   each group's published area in the `catchments` table.
