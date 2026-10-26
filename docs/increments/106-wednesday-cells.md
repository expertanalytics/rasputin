# Increment 106: four fixes for Wednesday's Swedish cells

**Status: designed** by `@architect`, 2026-10-26, on master `f17b3e11`,
worktree `demo-fixes`. Five small PRs, built in the order 106a, 106b, 106c,
106d, 106e (106e designed on 106c's head `a9366418`, worktree `demo-fixes-e`). Number 106: 105 is the highest under `docs/increments/` over every ref
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
| 106c | every cell has `catchment_id` 1; tagging cells afterwards puts 18.9 % of Kultsjön's area in the wrong one of its 75 SVAR units | `mesh --catchments` repairs SMHI's units, simplifies them with rasputin's own border simplification at a band set by the triangle budget, makes their borders constraints, so no cell straddles two, and writes each cell's unit id | about 95 | full loop | **yes** |
| 106d | no Swedish reservoir polygons: the reservoir share is 0 | a catalogue entry can keep only the lakes a dam names; `smhi-reservoirs` gives Kultsjön 82.73 km² of reservoir | about 35 | short | no |
| 106e | an official outline cannot be meshed at a small budget (Ljungan's forces at least 190,387 triangles), and NVE's one-part MultiPolygon outlines are refused: both took scripts outside rasputin | `mesh --simplify-domain` simplifies the outline with 106c's band from the budget; a one-part MultiPolygon is read as its polygon; the refusals, `sources` and the MCP instructions name the route, and a missing-blocks refusal names a fetch that covers the simplified outline | about 60 | full loop | **yes** |

About 225 net production lines in all. No C++ in any of them.

`@perf` applies to 106c and 106e only: 106c adds constraint lines to what
drives the mesh, 106e changes the domain when its flag is given (*Acceptance* in `docs/increments/README.md`). 106a changes only the
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

**Round 3** (2026-10-26), 106c only, on 106a's head `ee1dda88`. Ola's
ruling: the border band must not be a fixed distance (round 2's 200 m). In
his words: "let's say I want to mesh the northern hemisphere. Then I guess
this catchment would get maybe 8 vertices and edges. Or 5, or 10". The
triangle budget sets the resolution, so the band follows from it, and a
border gets about as many vertices as the triangle edges along it, at any
scale. `--features-tolerance` stays an explicit override. Section 7's
question is withdrawn.

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
   them, every unit keeps its area, and no gap can open. The frame's margin
   is 2 x band + 1 km, so it always exceeds the band. `METRES` is the
   existing `--features-tolerance` when given (its help gains "and
   --catchments borders"), and then the land cover's band too; 0 repairs
   only. When not given, the band comes from the triangle budget:

   `h = sqrt(4 A / (sqrt(3) N))`, `band = h / 8`

   with `A` the units' union area after the repair (m², in the units' CRS)
   and `N` the `--triangles` value. `h` is the edge of an equilateral
   triangle of area `A / N`, the mean triangle's edge at that budget. On
   SVAR's borders the simplification leaves segments about 8 bands long
   (*Measured*, round 3), so `h / 8` makes a border's segments about one
   triangle edge long: a border of length `L` gets about `L / h` vertices,
   at any scale. Kultsjön (1,702.67 km²) at 4,000 triangles: h 991 m, band
   124 m; at 400: 3,135 m and 392 m. The land cover then keeps its band
   derived from the legend (`legends.derived_cleanup`), so `--catchments`
   changes neither the land cover's band nor the 15-degree start-angle rule.
   Without `--triangles` (the `--tolerance` path) there is no budget and so
   no edge length to follow: `--catchments` then needs
   `--features-tolerance` (refusal 7). A pure function
   `band_m(area_m2, triangles)` in `sub_catchments.py` holds the formula.
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
borders simplified at 101 m (a triangle's edge at 6000 triangles, 809 m,
over 8): 2701 vertices, 1702.67 km2`; with `--features-tolerance` given, the
bracket reads `(--features-tolerance)`. The run record: `catchments` (the
source's name), `catchments_units` (the count) and `catchments_band_m` (the
band used, given or derived). The start mesh's minimum
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
7. Neither `--triangles` nor `--features-tolerance`: "needs --triangles,
   whose budget sets the border band, or --features-tolerance".

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

**Round 3: the band from the budget** (106a's head `ee1dda88`, the
extension rebuilt; probes `budget_band.py` and `budget_borders_plot.py`).
SMHI's main sub-basins (SVAR's units joined by `BARONR`, after the 1 m
repair): Kultsjön is one (1,702.67 km², its borders the outline, 19,555
vertices, 359.6 km), Ljungan four (12,838.64 km², 136,911 vertices,
2,346.3 km). Segment length against band, Kultsjön: band h/4 gives segments
1.95 to 2.23 h, h/2 gives 3.81 to 4.32 h, h gives 7.01 to 8.90 h (budgets
400, 4,000, 40,000), so segments run about 8 bands; hence h/8. At h/8:

| sub-basins | `--triangles` | h | band | border vertices | border length / h | triangles | mesh vertices on the outline |
|---|---|---|---|---|---|---|---|
| Kultsjön, 1 | 400 | 3,135 m | 392 m | 90 | 96 | 400 | 98 |
| | 4,000 | 991 m | 124 m | 290 | 312 | 3,963 | 329 |
| | 40,000 | 314 m | 39 m | 862 | 1,014 | not run | |
| Ljungan, 4 | 400 | 8,610 m | 1,076 m | 216 | 198 | 408 | |
| | 4,000 | 2,723 m | 340 m | 732 | 666 | 4,050 | |
| | 40,000 | 861 m | 108 m | 2,328 | 2,210 | not run | |

Border vertices against length / h: 0.85 to 1.10 over two areas eight times
apart and a hundredfold in budget. Every row: each sub-basin's area
unchanged (relative under 1e-14), coverage valid, the union one polygon with
no hole and its area unchanged, so no gap. The whole probe, reading and three
budgets, takes 1.0 s (Kultsjön) and 10.1 s (Ljungan); meshing 1.2 to
4.1 s. No triangle has more than 1.5e-5 of its area in a
second sub-basin (`straddle.py`), every centroid in one. The picture,
`docs/increments/106-probes/kultsjon-budget-borders.png`: SMHI's outline,
the simplified border and the mesh at 400 and 4,000.

SVAR's **fine units** at h/8 need a larger budget, because each unit's
junctions and short borders are fixed: Kultsjön's 75 at `--triangles 4000`
(band 124 m) is refused by the existing check, "the fewest triangles this
DEM gives is 4675"; at 6,000 (band 101 m, 2,701 vertices) it gives 6,057
triangles in 1.2 s. Ljungan's 809 at 40,000 (band 108 m) is refused at
54,645; at 80,000 (band 76 m, 36,970 vertices) it gives 78,046 in 30.9 s.
The refusal already names the budget to ask for.

**Constants.** The 1 m repair distance assumes a projected CRS in metres and
SVAR's mismatches; checked on Ljungan's 809 units (the largest input), where
it changes no area. The frame margin, 2 x band + 1 km, only has to exceed
the band; checked at bands 39 m to 3,135 m. The factor 8 assumes natural
catchment borders, rough at every scale probed (SVAR's); it is scale-free, so a hemisphere's budget gives a hemisphere's band
by the same rule; checked from 1,702.67 to 12,838.64 km² and 400 to 40,000
triangles. A border smoother than SVAR's (a straight line) gets fewer
vertices, never more.

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
7. **Refusals** 1, 2, 3, 4, 6 and 7 above.
8. **The band from the budget.** Pure: `band_m(4_000_000, 400)` is
   18.996 m (h 151.967 m), relative 1e-9; `band_m(4_000_000, 40_000)` is a
   tenth of it. End to end: test 5's run without `--features-tolerance`,
   with `--triangles 200` (A 9,000,000 m²): the record's
   `catchments_band_m` is 40.296 m (h 322.371 m), relative 1e-9, and the
   stderr line says "a triangle's edge at 200 triangles"; the same run
   with `--features-tolerance 0` records 0.

**Invariant-critical suite** (mutation round required): tests 1, 4, 5, 6
and 8.
Mutation targets: the inner-point rule (centroid or intersects in its place:
test 1's U), the frame left out (test 4), the borders not added as
constraints (test 5's slanted borders), the ids taken from the first
triangle vertex in place of the centroid (test 5's absent stderr line: a
vertex on a border is within no unit), the fid range check (test 3), the
band's factor (h or h/4 in place of h/8, and A/N in place of the
equilateral triangle's area: test 8).

**`@perf`**: the 1 m benchmark (no `--catchments`) unchanged, meshes
byte-identical to master's; with `--catchments smhi-svar2022-subcatchments`
and the band from the budget, Kultsjön at `--triangles 6000` and Ljungan at
`--triangles 80000`, the predicted counts (6,057 and 78,046, within the
budget's 5 %), at most 8.1e-5 and 3.0e-4 of any triangle outside the unit
holding most of it (unmeasured at these budgets; `@perf` sets them), and the time against the probes: mesh 1.2 s and
30.9 s, plus reading and repairing 0.3 s and 3.0 s, plus the
simplification 0.6 s and 6.1 s (unmeasured; `@perf` sets them), plus 10 % noise; and Kultsjön at
`--triangles 4000 --features-tolerance 200` gives round 2's 3,972.

**Size.** About 95 net lines: the option and its checks (12), a pure module
`sub_catchments.py` (select, repair, the frame and `simplify_borders`, the
band from the budget, the union, the borders as a feature, the ids: 43),
the wiring in `mesh` (20), the writer's `catchment_id` column (8), stderr
and record (12). Full loop.

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
  --land-cover clc2018 --reservoirs smhi-reservoirs --triangles 6000
```

With the band from the budget, Kultsjön's 75 fine units need about 6,000
triangles (4,000 is refused, naming 4,675); its one main sub-basin meshes
at 400. Each cell's `catchment_id` is then its unit's `fid` in SMHI's GeoPackage,
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

None open. Round 2's question (a 200 m default band) is answered by Ola's
ruling of round 3: the band follows the triangle budget.

## 8. 106e: the canonical route, outline to a 500-triangle mesh, with no step outside rasputin

**Round 4** (2026-10-26), on 106c's head `a9366418`. Ola: "We either
autocatchment or use an official outline. Then, we get 500 cell meshes. I
want this to be useful for Sol, only using the MCP. So we better make sure
there are no hidden secrets." On simplifying the outline to the budget: "No,
mesh --triangles N should not _always_ do that. But it _could_ be an option
to do that." And: "They do not need to be in the same geopackage": each
basin is its own run and file, so 106e has no assembly step and no
per-basin budget in one run.

**Round 2** (2026-10-26), on `f2a0abfe`: after round 1's review, item 4
(the missing-blocks refusal names a `--bbox` fetch, for the flag and for
`--catchments`), its probe and test 9; two wordings.

### The two routes, and where each stops today

- **Terrain:** `catchment --dem glo30 --seed ...` writes one `Polygon`, no
  holes, in the DEM's CRS (EPSG:4326 for glo30), its outline reduced at twice
  the DEM cell. Kultsjön lake-seeded (Opal's
  `catchments/kultsjon-lake-seeded-check.geojson`): 617 vertices; `mesh
  --domain` it `--dem glo30 --out-crs EPSG:3006 --triangles 500` is
  **refused**: "the fewest triangles this DEM gives is 1555".
- **Official:** `network-catchment --network smhi-svar2022-subcatchments`
  writes one `Polygon`, holes filled, in the network's CRS (EPSG:3006),
  unsimplified. Ljungan (Opal's `catchments-original/ljungan.geojson`):
  111,248 vertices; the same `mesh` call is **refused** at 190,387. Opal's
  delivery simplified its outlines with its own script ("topology-preserving
  simplified by 50 m", its `METHODS.txt`).
- **NVE's gauge outlines** (HydrologiskeData3, total catchment to a gauge)
  come as a `MultiPolygon` of one part; `--domain` refuses it ("the domain
  must be one Polygon, got MultiPolygon"), and Glomma's delivery unwrapped
  it by hand (its `METHODS.txt`: "The original one-part MultiPolygon was
  unwrapped to Polygon").

### What changes

1. **`mesh --simplify-domain`** (a flag, off by default). The `--domain`
   outline, holes included, is simplified as 106c simplifies the units' union:
   `sub_catchments.simplify([polygon], band)`, the polygon and its frame as
   one coverage, so the area is kept and every ring moves at most the band.
   The band is `--features-tolerance` when given (0: unchanged), else
   `sub_catchments.band_m(A, N)` (h/8, *What it does* in section 4) with `A`
   the domain's area in the CRS of the simplification and `N` the
   `--triangles` value. The CRS: the domain's own when projected; when in
   degrees, the Lambert azimuthal equal-area projection centred on the
   domain's centroid, on the domain CRS's ellipsoid, there and back vertex by
   vertex with `DomainPolygon.to_crs` (an equal-area CRS, so the area kept is
   the true area). It runs right after the domain is read, before the DEM is
   opened, so the window, the land cover's clip and the outline rule all see
   the simplified outline. Stderr, one line: `domain: outline simplified at
   963 m (a triangle's edge at 500 triangles, 7701 m, over 8): 111248 to 195
   vertices, 12838.64 km2`, with `(--features-tolerance)` in the bracket when
   given. Record: `domain_band_m`. Help: "Simplify the --domain outline, each
   ring moved at most a band and the area kept: --features-tolerance, or else
   an eighth of a triangle's edge at --triangles. --catchments always
   simplifies its outline."
   **With `--catchments`:** the domain there only selects units, and the
   outline is the units' union, already simplified at the same band; the
   flag is accepted and changes nothing (simplifying `D` first could change
   which inner points it holds).
   Refusals, headed `--simplify-domain`: without `--domain`, "needs
   --domain"; with neither `--triangles` nor `--features-tolerance`, "needs
   --triangles, whose budget sets the outline band, or --features-tolerance".
   Today's band from the land-cover legend is untouched: the flag reads
   `--features-tolerance` only as 106c's `--catchments` does.
2. **A one-part `MultiPolygon` domain** is read as its one polygon
   (`io/domain_file.py`); more parts are refused: "the domain must be one
   Polygon, got a MultiPolygon of 2 parts".
3. **The route is named where Sol looks** (the hidden steps found):
   - the fewest-triangles refusal, on a `--domain` run without the flag or `--catchments`, adds
     "; or give --simplify-domain to simplify the outline to the budget"
     (`cli.py`, `UnderFloor`);
   - `sources`' `fits` for polygon sources lists `mesh --reservoirs`, and
     `mesh --catchments` when the source's CRS is projected; today it names
     only `catchment --lakes` and `network-catchment --network`, while the
     MCP instructions tell the model to read `fits`
     (`source_report.fits`; vector `Facts.projected` is None today, set from
     the CRS the report already reads);
   - one sentence in the MCP instructions (`mcp/server.py`, `INSTRUCTIONS`):
     "catchment and network_catchment write a polygon that mesh reads as
     domain; with simplify_domain its outline follows the triangle budget."

   Checked and not hidden, so left out: the fetches without the flag or
   `--catchments` (each refusal names its `fetch`; with either, item 4),
   `network-catchment`'s seed (a unit's `ARO_UUID`, or a point with
   `--outlets smhi-outlet-points`; the `sources` descriptions name both
   attributes), the cells layer needing `--land-cover` and `.gpkg` (its
   refusal says so), and the CRS of `catchment`'s output (any CRS is read).
4. **The missing-blocks refusal names a fetch that covers the outline
   meshed** (round 2). `fetch --domain D` copies the DEM blocks of `D`'s
   bounds grown by 6 cells (`fetch/plan.py`, `source_box`: a margin of 4
   plus 2; about 186 m for glo30). With `--simplify-domain` the outline can
   pass `D`'s bounds by up to the band, and with `--catchments` the units'
   union passes them wherever `D` cuts a unit whose inner point it holds
   (*Measured*: up to 477 m and 868 m). Today mesh's refusal (`cli.py`, the
   `except NotCached` after `assemble_dem`) names `fetch` with the original
   `--domain`, which adds no block: an agent with only the MCP loops. With
   `--simplify-domain` or `--catchments` given and the exception's box set,
   the refusal names `rasputin fetch KEY --bbox X0 Y0 X1 Y1 [--cache C]`,
   the box being `NotCached.box` (the window's node box in the source's
   CRS) and no `--out-crs` (a `fetch --bbox` is in the source's CRS
   without it), as `catchment` does today (`cli.py`, its `except
   NotCached`). That hint moves into one helper both commands call, so
   `catchment`'s wording stays byte for byte. Without either option the
   refusal is as today. The feature sources (`--catchments`' units,
   vector land cover, reservoirs) need nothing: a catalogue fetch keeps
   every polygon meeting `D`'s box whole, and the units selected meet it.

### Measured (probes, 106c's head `a9366418`, extension built 2026-10-26)

`outline_budget.py` (each Opal outline as SMHI's network joins it) and
`outline_degrees.py` (an outline in degrees, via the centred equal-area
projection), both in `docs/increments/106-probes/`:

| outline | km² | vertices | N | band | after | area change | simplify |
|---|---|---|---|---|---|---|---|
| Kultsjön, official | 1,702.67 | 19,555 | 500 | 350.5 m | 106 | 2.4e-7 m² | 0.21 s |
| Ljungan, official | 12,838.64 | 111,248 | 500 | 962.6 m | 195 | 3.2e-5 m² | 2.43 s |
| | | | 5,000 | 304.4 m | 632 | 1.3e-5 m² | 1.24 s |
| Umeälven, official | 26,771.60 | 133,783 | 500 | 1,390.0 m | 138 | 3.1e-5 m² | 3.90 s |
| Kultsjön, terrain (EPSG:4326) | 1,097.30 | 617 | 500 | 281.4 m | 122 | relative 2.4e-8 (geodesic) | under 0.01 s |

All seven Opal outlines: valid, one polygon, no hole, at 500 and 5,000.
Meshed (`--dem glo30 --out-crs EPSG:3006 --triangles 500`): Kultsjön
official 490 triangles in 2.6 s; Ljungan 503 in 4.0 s (both runs, the
simplification apart); Kultsjön terrain 502; and with `--land-cover clc2018`
to `.gpkg` Kultsjön official writes 490 cells (46.5 s, nearly all 106a's
land-cover repair). Umeälven's glo30 blocks are not cached (not run).

`fetch_reach.py` (round 2; 3 min 18 s for all eight): how far each
outline's bounds pass `D`'s, the box `fetch --domain D` grows by about
186 m for glo30; N 500, EPSG:3006, the SVAR 2022 units from
`rasputin_data/sweden_smhi_svar/wfs-20261021`:

| `D` | `--simplify-domain` | `--catchments` (units) |
|---|---|---|
| Åselealven | 428.0 m | 428.0 m (726) |
| Faxälven | 334.9 m | 567.3 m (471) |
| Fjällsjöälven | 329.8 m | 114.1 m (464) |
| Indalsälven | 158.9 m | 868.1 m (1,311) |
| Kultsjön | 97.3 m | 97.3 m (75) |
| Ljungan | 242.5 m | 221.0 m (809) |
| Umeälven | 476.8 m | 302.0 m (1,436) |
| Kultsjön, terrain | 101.6 m | 107.8 m (51) |

Past 186 m: 6 of the 7 official outlines with the flag, 5 of 7 with
`--catchments`; so 106c's `--catchments` has the same gap, and item 4
covers both.

**Constants.** None new: the band is 106c's h/8, the frame and the 1 m
clearance are 106c's. The equal-area projection keeps area exactly and
stretches distance by under 1 % within about 1,800 km (16 degrees of arc)
of its centre, so a
band in degrees-CRS domains is right to that scale; checked on a 1,097 km²
outline, and the projected path up to 26,771.60 km².

### Tests (key properties)

Fixture `ZZ`: in `test_cli_mesh_catchments.py`'s scene (56 x 56 nodes at
60 m, EPSG:31983, `local` frame), the square [0, 3000]² with each side a
zigzag of 60 steps of 50 m, every other vertex pushed 20 m outward: 240
vertices, valid, 9,120,000 m². `ZZH`: `ZZ` with a hole, the square [1200,
1800]² zigzagged the same way, 12 steps a side, every other vertex 20 m
away from the hole's centre: 288 vertices,
8,736,000 m².

1. **Read** (pure): a one-part `MultiPolygon` of `ZZ` reads as `ZZ`,
   oriented (`equals_exact`, tolerance 0); a two-part one is refused with
   the wording above; a `Polygon` reads as today.
2. **Simplify, projected** (pure): `ZZ` at N 100, band 57.366336 m
   (relative 1e-9): 4 vertices, area relative change under 1e-12, valid.
   `ZZH` at N 100 (band 56.145636 m): outer ring 4, hole 4, area kept
   (the "frame left out" mutant leaves the hole's 48 and the outer 240).
3. **Simplify, in degrees** (pure): `ZZ` moved to EPSG:4326 (centroid
   -44.93, -15.24): the result is in EPSG:4326, 4 vertices, geodesic area
   relative change under 1e-6 (measured 3e-8), band 57.4 m (relative 1e-3).
   The band from an area in square degrees (about 5e-4) fails this.
4. **End to end**, `ZZ` as `--domain`, `--triangles 100`: without the flag
   refused, "the fewest triangles this DEM gives is 470", and the refusal
   names `--simplify-domain`; with it, 101 triangles (within the budget's
   5 %), the stderr line "240 to 4 vertices", `domain_band_m` 57.366336
   (relative 1e-9); with `--features-tolerance 30`, `domain_band_m` 30.
5. **Unchanged without it:** `ZZ` simplified by hand (the probe's output at
   N 400) meshed at `--triangles 400` without the flag gives 396 triangles
   and the same rows as 106c's head; no `domain_band_m` in the record.
6. **With `--catchments`**: test 5's units run of section 4 with the flag
   writes the same `triangles` and `cells` rows as without it.
7. **Refusals:** the two above.
8. **Named where Sol looks:** `fits` of a projected polygon `.geojson` lists
   `mesh --reservoirs` and `mesh --catchments`, of one in EPSG:4326 only
   `--reservoirs`, of a points file neither; `INSTRUCTIONS` holds
   "simplify_domain".
9. **The fetch named covers the outline** (round 2): `test_cli_mesh_cache.py`'s
   `TestW5TheFetchCommandInTheMessage.test_with_domain` scene (the cache of
   `test-utm33` without block 6, the triangle domain) with
   `--simplify-domain --features-tolerance 5` added: refused, "1 of the 6
   blocks needed are not in the cache for the window 500150 6599880 500400
   6599960", and the output holds `run: rasputin fetch test-utm33 --bbox
   500150 6599880 500400 6599960 --cache <cache>` and no `--domain` (the
   window measured without the flag on 106c's head; a triangle keeps its
   three vertices at any band). The same run with `--catchments U` in place of `--simplify-domain`, `U` a GeoJSON of one unit (the triangle, `id` 1), names the same `--bbox` run. The existing test without the flag keeps
   its `--domain` run unchanged, and `catchment`'s hint tests pass
   unchanged.

**Invariant-critical suite** (mutation round required): tests 1, 2, 3, 4.
Mutation targets: the one-part check (any `MultiPolygon` taken as its first
part: test 1), the frame left out (test 2's hole), the simplification in
degrees (test 3's band), the band's factor and its area (test 4's
`domain_band_m`), the flag ignored (test 4's count), the refusal's
`--domain` kept with the flag (test 9), or with `--catchments` (test 9's second run).

**`@perf`** (the domain drives the mesh): the 1 m benchmark unchanged and
byte-identical to 106c's head (no flag); Ljungan official at `--triangles
500` with the flag: 503 triangles, the time against the probes, 4.0 s mesh
plus 2.4 s simplification plus reading the outline (unmeasured;
`@perf` sets it), plus 10 % noise; Kultsjön terrain-traced: 502.

**Size.** About 60 net lines: the flag and its two checks (14), the
simplification with its CRS choice, stderr and record (18), the one-part
read (3), the refusal's hint (4), `fits` and vector `Facts.projected` (8),
the instruction sentence and the record wording (8), the `--bbox` fetch
helper shared with `catchment` and its use in mesh (about 5). Full loop: over 50
lines, and the domain drives the mesh.

*Prior art.* As section 1: rasputin's border simplification
(`docs/increments/32-landcover-simplify.md`) through 106c's frame; the
Lambert azimuthal equal-area projection (Snyder, *Map Projections: A Working
Manual*, USGS Professional Paper 1395, 1987, ch. 24), used as published.
*Legacy:* `git grep -n -i -E "simplif|MultiPolygon" legacy-archive -- legacy`
returns CGAL's surface-mesh edge collapse (`legacy/rasputin/mesh.py:102`,
legacy-archive's `legacy/rasputin/triangulate_dem.h`, line 243) and a raster-polygon intersection type
(the same file, lines 73 and 410); no outline simplification and no domain
reader. Nothing to carry across. *Ideas file:* overlaps "MCP tool descriptions: how to mesh
a domain too big for one run" (the instruction sentence only; the basin path
is not touched).

## Review record

- **Design review, round 1, 2026-10-26.** `@reviewer`, head `89e22445b28122ba6b1915535a632cc1c48a3a73`. Verdict: CHANGES REQUESTED. Blocking: (1) 106c test 1's square units let a centroid rule pass: add a non-convex unit; (2) 106c test 5 lets "ids from the first vertex" and "borders not added" pass: assert the nearest-unit stderr line absent and slant a border; (3) scope: cut 106c to Wednesday's needs, `named_by` to the entry's. Not pushed; no CI.
- **Design review, round 2, 2026-10-26.** `@reviewer`, head `ff1ab7fd7ab227c24d0ddc51e8c1901a9aec533b`, range `89e22445..ff1ab7fd`. Verdict: APPROVED WITH WORDINGS. Round 1's three findings resolved; one wording applied (which band the land cover gets with `--catchments`).
- **Design review, round 3, 2026-10-26.** `@reviewer`, head `706a3f6b29d8f6637f3ba398f644ac2389c93563`, range `ee1dda88..706a3f6b`. Verdict: APPROVED WITH WORDINGS. Band from the budget (h/8) checked and reproduced; three unrecorded acceptance figures marked unmeasured for `@perf`.
- **Design review of 106e, round 1, 2026-10-26.** `@reviewer`, head `f2a0abfe29780c9cb32f075a52fd251d727219b6`. Verdict: CHANGES REQUESTED. Blocking: a fetch loop, the simplified outline passing the fetched box (up to 477 m, Umeälven at 500) while the missing-blocks refusal names the original `--domain`; fix: name `--bbox` from the missing-blocks box, as `catchment` does, one test. Wordings: the table's "about 55", point 3's "or `--catchments`".
- **Design review of 106e, round 2, 2026-10-26.** `@reviewer`, head `a5399b83dad9c56f7ea5006dcbfcbb1edb7986e5`, range `f2a0abfe..a5399b83`. Verdict: APPROVED WITH WORDINGS. Fetch loop closed (refusal names the planned window's `--bbox`); probe reproduced; two wordings applied (test 9's `--catchments` run and its mutation target).
