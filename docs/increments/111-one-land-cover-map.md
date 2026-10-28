# Increment 111: one land-cover map, and the edges chosen from it

**Status: designed, approved in design review round 6; section 8 (the Python interface and the wordings tests read) added for the test migration** by `@architect`, 2026-10-28, on master `75ab2b4f`,
worktree `land-cover-edges`. One PR, full loop, Python only. Number 111: 110
is the highest under `docs/increments/` over every ref and worktree.

Ola, 2026-10-28, on the split below: "Yes, this is logical. Now we need to
make it consitent :)"

## 0. The plain picture

Three questions, three sets of options:

1. **What is on the ground?** Areas: the land-cover map, `--land-cover`
   (CORINE's polygons, AR50, a class raster), with `--prefer KEY=SOURCE`
   taking one class from another source (NVE's lakes and glaciers).
   `--land-cover` takes several class-polygon files under one class file
   (AR50 comes one file per county), read as one map. Lines:
   `--lines KEY=SOURCE` (rivers, roads), with or without a land-cover map.
2. **Which of those must the mesh follow?** The edges: `--edges`,
   chosen per class or line kind (`--edges lake --edges river`; `--edges
   all`). Edges come only from the areas and lines above. `--edges-*` shape
   them (today's `--features-*`).
3. **What does the mesh file say about the ground?** Each triangle's
   exact share of each class, and in a `.gpkg` the cells Shyft reads, as
   `--land-cover` gives them today. A run that follows class borders
   (`--edges` naming a class) writes a class code per triangle instead, as
   a `--features` run with CORINE does today, at today's cost; `--shares`
   asks for the shares as well. So every recorded run that does not use
   `--features` keeps its command, and none gets slower.

The catchment outline (`--domain`) and sub-catchment borders (`--catchments`)
stay their own concern. `--features` goes. What it carried moves to the
areas and lines; the one route with no home, a mesh's own `lines` layer fed
back, is answered by the direct command (section 4, two-stage probe).
Sula's direct command becomes

    rasputin mesh ... --land-cover clc2018 --prefer lake=nve-lakes --edges all
        --shares --edges-repair 0.1 --edges-tolerance 20 --edges-outline-snap 0

with no rules file to upload. A mesh that follows rivers only:

    rasputin mesh ... --lines river=rivers.gpkg --edges river

Today's quick-check case on Numedalslågen becomes `--land-cover <file>
--land-cover-layer corine2018 --land-cover-classes corine --edges all` and
costs what it costs today (section 6).

## 1. Prior art: legacy and literature

*Legacy.* `git grep -l -i -E "land.?cover|corine|constraint" legacy-archive -- legacy`
returns `application.py`, `land_cover_repository.py`, `globcov_repository.py`,
`gml_repository.py`, `wfs_repository.py`, `tin_repository.py`,
`web_visualize.py`, `triangulate_dem.h` and two tests: one land-cover source
looked up per triangle after meshing, never a constraint. Nothing to carry.

*Literature.* Composing one map from several sources, theme by theme, is map
conflation (Saalfeld, "Conflation: automated map compilation", *Int. J. GIS*
2(3), 1988, 217–228); the rule is 107's and 110's, unchanged. Keeping the
borders of chosen classes and dropping the rest is *selection* among the
generalisation operators (McMaster and Shea, *Generalization in Digital
Cartography*, AAG, 1992). No novelty claimed; none searched.

*Ideas file.* Overlaps "Line layers in the ArcGIS fetch": rivers from the
catalogue need it; this increment reads line files and cached line layers,
and fetches none.

## 2. Inventory: what assumes the old way, and what it becomes

| Today | Becomes |
|---|---|
| `mesh --features` (repeatable) | **deleted**: class polygons go to `--land-cover` (repeatable, section 3), lakes and glaciers to `--prefer`, rivers and roads to `--lines` |
| `--features-crs`, `--features-layer` (paired, repeatable) | `--land-cover-crs`, `--land-cover-layer`: one value each, for every `--land-cover` file. No layer or CRS option for `--prefer` or `--lines` files: none of the 732 recorded runs gives one (section 3) |
| `--features-rules`, built-in `property` (also `--features`' default, `cli.py:991`; `source_report.py:384`) | **deleted**: what classes mean is `--land-cover-classes`, what the mesh follows `--edges`. No recorded run uses `property` or `--features` without rules (732 commands) |
| a mesh's own `lines` layer fed back (`--features <mesh>.gpkg --features-layer lines` with rules on its `land_cover` and `water` columns; Ola's Sula runs of 27 October at 1, 2, 5, 10 and 20 m) | **no route**: the direct command gives the same borders and the same shares (section 4, two-stage probe) |
| `--features-repair`, `--[no-]features-merge-same-class`, `--features-tolerance`, `--features-outline-snap` | `--edges-repair`, `--[no-]edges-merge-same-class`, `--edges-tolerance`, `--edges-outline-snap`; same meanings; defaults from the land-cover class file's accuracy, as for a coded source today. `--edges-tolerance` also works without `--land-cover`: it sets the `--catchments` and `--simplify-domain` band, as `--features-tolerance` does today |
| `--land-cover`, `--land-cover-year` | kept; `--land-cover`, now repeatable, says what is on the ground, and gives the shares unless the run follows class borders without `--shares` (section 3) |
| `--land-cover-rules` | `--land-cover-classes` |
| `--prefer` | kept, area keys only (lake, glacier), still needs `--land-cover`; its sentence on `--features` lines deleted |
| (none) | `--lines KEY=SOURCE`, repeatable: rivers, roads and the other line kinds, with or without `--land-cover` (section 3) |
| (none) | `--shares`: the shares as well, in a run that follows class borders (section 3) |
| `--reservoirs`, `--catchments`: "needs --land-cover with meshing rules and a .gpkg output" (`cli.py:1689`) | unchanged without class borders in `--edges`; with them, also "needs --shares" |
| (none) | `--edges SELECTOR`, repeatable (section 3) |
| `RENAMED`, `_RenamedOptions` (`--features-legend`, `--land-cover-legend`) | **deleted**: their targets go |
| meshing-rules files (`data/meshing_rules/`, `meshing_rules.py`, "meshing rules" in help and docs) | **land-cover class files**: `data/land_cover_classes/`, `land_cover_classes.py` (beside `land_cover_fractions.py`), "land-cover class file" |
| land-cover class file keys `[edges]`, `class.edges`, `class.value`, `[vocabulary] extra` | **deleted** (edge names follow from the class, section 3) |
| `unlisted = "land"` | `unlisted = "keep"`; `drop`, `refuse` kept; `code_system` required |
| `water`, `[cells]`, `[[split]]`, `codes`, `accept_no_class`, `positional_accuracy_m`, `notice`, `licence`, `legend`, `name`, `attribute` | kept |
| `corine-water.toml` | **deleted**: `--edges water` |
| catalogue key `rules =` (shipped `sources.toml`: 2 entries; Ola's `~/.config/rasputin/sources.toml`: 15, all `"ar50"`) | `classes =` (Ola, question 1: rename); `rules` becomes an unknown key, refused (`extra="forbid"`, `sources.py:106`). The shipped file changes in the PR. Ola's own file is outside the tree: the main session edits it (`rules =` to `classes =`, 15 lines, a `.bak` copy first) with Ola's yes, right after the merge and the main checkout's fast-forward, before any command or MCP server reads it again; the PR description says so |
| `sources.py:152`: a `legend` key refused, naming `rules` | **deleted**: an unknown key, refused like any other |
| record `features`, `features_joined` | `edges`, `edges_joined` |
| record `features_crs`, `_transform`, `_accuracy_m`, `_notice`, `_licence` | `land_cover_crs`, `_transform`, `_accuracy_m`, `_notice`, `_licence` |
| record `features_repair_m`, `_merge_same_class`, `_tolerance_m`, `_outline_snap_m` | `edges_repair_m`, `edges_merge_same_class`, `edges_tolerance_m`, `edges_outline_snap_m` |
| record `land_cover_prefer` | kept; a new `lines` field, one row per `--lines` key ("river from X: N lines, L km"); a new `land_cover_overlap` field with more than one `--land-cover` file (section 3) |
| MCP `mesh` arguments `features*`, `land_cover_rules` | follow the CLI (the schema is derived): `edges`, `edges_*`, `lines`, `shares`, `land_cover_classes`, `land_cover_crs`, `land_cover_layer` |
| MCP `INSTRUCTIONS`: "features_rules", "features_layer" | class polygons to `land_cover` with the fit's classes; lakes, glaciers to `prefer`; rivers, roads to `lines`; `edges`; `shares` for class shares and cells |
| `sources` fits: `mesh --features` with rules; `mesh --land-cover` (class raster) | `mesh --land-cover` with classes (class polygons), `mesh --prefer` (other polygons), `mesh --lines` (lines); the class raster's `mesh --land-cover` unchanged (it gives shares, as today); `Fit.rules` → `Fit.classes` |
| `FeatureSource`, `open_features`, `_Tally.source`, `_take`, `pre_clip`, `_runs`, `_bound`, `_linework`, `R_MIN`, `R_MAX` | **deleted**: reading is `read_class_polygons` (reading as `_take`'s cleaned-polygon branch does: one query, one move, no validity check; its `_repaired` moves to the shares step, section 3), `read_polygons`, new `read_linestrings` |
| `FeatureRequest` | `EdgeCleanup(repair_m, merge_same_class, tolerance_m, outline_snap_m)` |
| `MeshingRules` | `LandCoverClasses`, less `classes` (value to edge names), `otherwise`'s names, `extra_names`, `cleaned` |
| `meshing_rules.run_vocabulary` | the default vocabulary, plus `catchment` with `--catchments` |
| #387: `cli._cover_lines` and its two refusals, `Prefer`, `open_features(prefer=)`, `_Tally._preferred` | **deleted**; the join moves into the map (section 3) |
| `preferred.overlay` | **deleted**; `preferred.replace` kept |
| `landcover.label_triangles`, `cli._land_cover` (codes from coded `--features`) | kept: without `--shares`, the code per triangle; with it, the codes come from the shares, as `--land-cover`'s do today |
| `gml-to-gpkg` help "for --features" | "for --land-cover" |
| `docs/benchmarks/quick/cases.toml` numedalslagen, lagan | `--land-cover <file> --land-cover-layer/-crs ... --land-cover-classes corine --edges all`, no `--shares`; names and inputs unchanged, so the stored baseline still judges them |
| `cells.py` (`MeshingRules`, "the meshing rules' `[cells]` table"), `legends.py` (docstring: legends named by meshing rules) | `LandCoverClasses`, "the land-cover class file's" |
| `io/domain_file.py`, `io/geojson.py` docstrings ("`--features`") | "`--land-cover`, `--prefer`, `--lines`" |
| `docs/GLOSSARY.md` *Meshing rules* (`--features-rules`, sources `rules` key, `property`) | *Land-cover class file*; a new *Edges* entry (`--edges`, against `--out-edges`, which writes them) |
| README, INSTALL, `project_structure.md`, ROADMAP | reworded |

## 3. Design

```
--land-cover (+ -classes, -layer, -crs)  --prefer KEY=SOURCE                     --lines KEY=SOURCE
        │ read near the read region, in the DEM's CRS (CLI, I/O)                     │ read_linestrings
        ▼                                                                            ▼
land_cover_map.build(base, preferred, classes, region) -> LandCoverMap | None   lines: ((key, line), ...)
        │ areas: (code, polygon)...                                     │            │
--edges, --edges-* → EdgeChoice, EdgeCleanup                            │            │
        ▼                                                               │            │
feature_input.edges_of(map, lines, choice, cleanup, domain, vocabulary) -> FeatureSet   [pure]
        │ + sub-catchment borders (unchanged) → mesh (C++ unchanged)    │
        ▼                                                               ▼
after the mesh: shares → assign_polygon_fractions(repaired map.areas) → shares, codes, cells
                class borders without --shares → label_triangles(chosen classes' areas) → codes
```

*Why the shares follow from the run, not from `--land-cover` alone.* Round
2 merged `--features`' class polygons (edges, a code per triangle) and
`--land-cover`'s (shares) into one `--land-cover`, and kept `--land-cover`'s
meaning: "the exact share of each class in every triangle", which made
every border-following run pay for shares (Ola, question 3: no). The edges
are drawn from the polygons before the mesh, the shares after it, and
neither needs the other. Round 3 made the shares opt-in, which changed the
12 recorded runs that pass `--land-cover` without `--features`, among them
the seven 26 October Opal runs (Åseleälven, Faxälven, Fjällsjöälven,
Indalsälven, Kultsjön, Ljungan, Umeälven, each with `--reservoirs`).
This round keeps both recorded habits: `--land-cover` alone gives the
shares, as those runs expect; a run that follows class borders gives
codes, as the 696 recorded `--features` runs without `--land-cover` (681
of them the Hafslund, Skagerak, Opal SE2 and Glomma batches) expect;
`--shares` asks for both, as Sula's 17 runs with `--features` and
`--land-cover` do. The one step the edges and the shares have in common, reading the polygons, does carry a cost trap
(section 4, read probe): `read_class_polygons` checks each whole polygon's
validity (`_repaired`), which the share overlay needs and the edges do not
(the clean-up clips each polygon to the read region and repairs it,
`_Tally._clean`). So the map is read without that check, and the check runs
in the shares step.

- **The map** (`tin_engine/land_cover_map.py`, frozen dataclass
  `LandCoverMap(areas, classes, values)`; None without `--land-cover`). The
  base is the class polygons over the read region (`source_region`), read
  as today's cleaned `--features` polygons are (`_Tally._take`: every
  polygon moved in one `shapely.transform` call, no validity check), or no
  areas for a class raster, read after the mesh as today. Each `--prefer`
  is 107's `replace` near the read region, first given wins, and each
  refilled piece is unioned with the same-class polygons it touches (110's
  join, moved from `_Tally._preferred`, on `(code, polygon)` pairs), on
  unchecked polygons as `_Tally._preferred` does today. The record's
  replaced, given and refilled areas are measured inside the domain.
- **Several `--land-cover` files** (repeatable): class-polygon files under
  the one land-cover class file (`--land-cover-classes`, or the catalogue
  entries' `classes =`, which must agree), as Ola's AR50 runs pass one
  file per county (Bygdin, Innlandet and ar50_46; Øyangen, ar50_32 and
  ar50_33). Their polygons, in the order given, are the base: one map, so
  the clean-up and the shares see one coverage, and `--prefer` cuts it as
  it cuts one file. Where two files overlap, the first given keeps the
  area (as `--prefer`'s first given wins): a later file's polygon is cut
  by the earlier files' polygons it shares area with (one
  `STRtree.query(predicate="intersects")` over each later file, keeping
  the pairs whose intersection has area above zero; `"overlaps"` would
  miss a polygon wholly inside another, and two equal polygons). The record's `land_cover_overlap`
  gives that area inside the domain and the file kept; `land_cover_source`,
  `land_cover_crs`, `_notice` and `_licence` list the files joined by `; `,
  as `features*` does today. Today `--features` cleans each file on its
  own, so the county line is an edge wherever it runs; here it is one
  coverage, so with `--edges-merge-same-class` (on for AR50 by its
  accuracy) the county line is an edge only where the class changes
  across it. On Bygdin the two files tile the domain exactly (0 m²
  overlap: of 345 intersecting pairs none shares area; 0 m² gap) and the county line runs 7,614 m inside it, all of it
  with one class on both sides: today's mesh follows it, the new one does
  not (question 4, section 7). Refused, naming `--land-cover`: a class
  raster with any other `--land-cover` source; two catalogue entries or
  files whose class files differ.
- **Lines** (`--lines KEY=SOURCE`, repeatable). `read_linestrings`
  (LineString parts, moved vertex by vertex) gives `((key, LineString), ...)`,
  separate from the map, so lines need no `--land-cover`. A key is an
  edge-vocabulary name other than `land_cover` and `water` (river, road,
  railway, coastline, contour, wall, ditch). A file works today; a
  catalogue key works once cached (the fetch of line layers is the ideas
  entry, left out).
- **Keys and refusals.** A `--prefer` key is one of the land-cover class
  file's `[cells.prefer]` (lake, glacier), as today. Refused, naming the
  option: `--prefer` without `--land-cover` or with a class raster (kept,
  `cli.py:2017`, `:2024`); `--prefer` with a line key (the message names
  `--lines`); `--lines` with an area key (naming `--prefer`); a polygon
  source under `--lines`; a line source under `--prefer`; a second source
  for one `--lines` key.
- **Layer and CRS of a file.** A catalogue key carries its entry's layer;
  the `--land-cover` file has `--land-cover-layer` and `--land-cover-crs`
  (the quick check's Numedalslågen and Lagan cases use them), passed to
  `read_source`'s `layer` and `crs` as today (`feature_input.py:656`). A
  `--prefer` or `--lines` file gets neither: no recorded run gives a layer
  or CRS for anything but the class polygons (`--features-layer`: `lines`,
  the fed-back route of section 4; `triangles` and `land_cover`, a mesh's
  classes, now `--land-cover-layer`). Such a file needs one features layer
  and a CRS of its own, or `read_source` refuses it as today, naming the
  layers. Refused, naming the option: either `--land-cover-*` option with a
  catalogue key or without `--land-cover`, a layer for a file that is not
  a GeoPackage, a CRS that is not the file's own.
- **`--edges`.** Selectors: `all` (every area class given and every
  `--lines` key), `water` (the land-cover class file's water classes), a
  category of its `[cells.categories]` (lake, forest, glacier, ...), a class
  code, a line key. Resolved at the CLI into `EdgeChoice(codes, lines)`. A
  border between two areas is an edge when either side's class is chosen;
  it carries `land_cover`, and `water` when that side is a water class (the
  masks every shipped land-cover class file but `corine-water` gives today).
  So `--edges water` sets both bits where `corine-water` set `water` only: a
  reader of the `land_cover` bit sees these lines too. A chosen line key's
  lines carry its name. Default: no edges. Refused, naming `--edges`: an
  unknown selector (the message lists the accepted ones), an area selector
  without class polygons from `--land-cover`, a line key without its
  `--lines`. A land-cover class file whose category is named `all`, `water`
  or a line key is refused when it loads, naming the category, so a
  selector means one thing. Help text: `--edges` says "which borders and
  lines the mesh must follow"; `--out-edges` says "write the mesh's edges to
  this PLY".
- **`edges_of`.** The chosen classes' areas go into the clean-up
  (`_Tally._clean`) as one coverage, and give the lines (`_add`); the other
  classes' areas stay out of it, as a class file that drops them
  (`corine-water`'s `land = "drop"`) keeps them out today. So `--edges all`
  cleans what `corine` cleans and `--edges water` what `corine-water`
  cleans. Lines are clipped to the domain, not cleaned, as `--features`
  lines are today; the `--edges-*` clean-up options shape area borders only.
  With no map, `edges_of` gives the lines alone.
- **The shares** are computed when `--land-cover` is given and either no
  `--edges` selector names a class or `--shares` (flag) asks. The map's
  areas are repaired near the mesh (`_repaired`, today's
  `read_class_polygons` step, now here), then `assign_polygon_fractions`
  gives each triangle's shares; or the class raster's route, as today. The
  codes per triangle come from the shares, and a `.gpkg` with the class
  file's `[cells]` gets the cells layer, both as `--land-cover` gives them
  today. `--shares` without `--land-cover` is refused, naming it.
- **Class borders without `--shares`**: a code per triangle from
  `label_triangles` over the chosen classes' cleaned polygons (0 where
  none), as a coded `--features` run does today. No shares, no cells; the
  record says "class shares: not computed (add --shares)". `--reservoirs`
  and `--catchments`' cell ids, which need the cells, are refused there,
  naming `--shares`; elsewhere they need `--land-cover` as today.
- **No class map** (`--lines` alone): no codes, no shares.
- **A line nothing uses** is refused, naming it and `--edges`: a
  `--lines KEY` that no `--edges` selector takes, as an unused input is
  refused elsewhere. `--land-cover` is never unused (it gives the shares or
  the borders).
- **A cut run** (`--pieces`) refuses the shares (pieces carry no class
  shares yet, as today: `--land-cover` without class borders is refused
  there, naming `--land-cover`, and so is `--shares`, naming `--shares`;
  both at the point where a cut run refuses `--stats` and `--out-edges`
  today, `cli.py`'s `part is not None` check: after the DEM's headers are
  read and the cut planned, before any DEM block is decoded, the map read
  or a piece meshed. Not moved earlier: only the plan says whether a run
  is cut, as `--piece-budget` or a small DEM may give one piece);
  `--land-cover` with class borders in `--edges`, and `--lines`, work as
  in an uncut run, and each piece gets its codes as `_write_piece` gives
  them today (`cli.py:2129`, codes by `label_triangles`). `--prefer` is
  allowed there: the map is built once before the cut, so the pieces' codes
  come from the preferred map, at no added code. Today a cut run refuses
  `--land-cover`, so `--prefer` never reached one; none of the 732
  recorded runs is cut (`--pieces`, `--piece-budget`: 0 runs).

## 4. Probe on Sula (NVE 73.27.0), 20 m

`docs/increments/111-probes/one_map_probe.py`, on AC power: step 1 runs
today's direct command (Sol's record,
`~/rasputin-output/mcp/ola/work/sula-mcp-direct-retest-20261028b/sula-20m-record.json`)
and keeps the map its lines are drawn from; step 2 meshes that map as the
only land cover, no `--prefer`; step 3 builds the map as above (repaired
CORINE, cut, refill, join) and meshes it.

| Run | Triangles | Not shared with today | Line vertices moved | Fractions, largest difference |
|---|---|---|---|---|
| today | 3,271 | | | |
| step 2 | 3,271 | 0 | 0 | 3.7e-9 |
| step 3, round 2's order (repair, then cut) | 3,271 | 14 | 2 of 1,201, at most 1.4 mm | 1.1e-11 |

Step 3 differs at two vertices only: its cut runs on the repaired polygons
(`_repaired` clips each self-crossing ring at the read box), today's line
copy on the unrepaired ones. Repairing today's map after the cut
gives today's mesh exactly, so the repair's order, not its effect on the map,
moves them; both maps give the same class areas inside the domain to 1 m².
That run is this round's order (cut on the unchecked polygons, repair in
the shares step, section 3): today's mesh exactly.
Every run's cells total 30.4411 km² with 8.935 % lake.
Planted check: the same command with CORINE's own lakes gives 2,025
triangles, so the probe tells maps apart.

*Two-stage probe* (scratch scripts, not kept; trees `132a89cf` and
`369bac65`, on AC power). On 27 October Ola's Sula runs at 1, 2, 5, 10 and 20 m took their
borders from a prepared mesh's `lines` layer
(`sula-73-27-0-features20m-fixed-20261027/prepared-constraints-20m.gpkg`,
cleaned at 20 m from a 1 m mesh's classes), a route the new options do not
take. The direct command (Sol's, above; `--tolerance 5`) against Ola's 5 m
run:

Round 3's table counted the catchment outline (34.36 km, the same in
both runs) among the edge lines, which raised every share of length
within a distance; the figures below are the class borders alone
(`lines` rows with `land_cover` set), sampled every metre, re-measured on
tree `369bac65` with the direct command rerun at `--tolerance 5` (10,283
triangles, as round 3).

| Lines, 5 m | Two-stage (27 Oct) | Direct | Within 1 / 5 / 20 m of the other run's, two-stage's; direct's | Beyond 5 m, two-stage's; direct's | Farthest |
|---|---|---|---|---|---|
| lake (`water` set) | 45.67 km | 45.66 km | 86.2 / 96.2 / 100 %; 86.1 / 96.3 / 100 % | 1,743 m; 1,700 m | 25.3 m; 25.5 m |
| other class borders | 25.29 km | 25.38 km | 97.3 / 99.3 / 100 %; 97.0 / 99.1 / 100 % | 171 m; 237 m | 14.5 m; 19.7 m |

A lake line is measured against the other run's lake lines, any other
border against the other run's nearest class border (against other
borders only, one stretch of the direct run's lies 49.8 m off: its
neighbour in the two-stage run is a lake line). Triangles 11,479 and
10,283; lake share 8.935 % in both. The length
beyond 5 m is mostly on lake lines: 1.7 km of 45.7 km (3.8 %) of lake
border lies 5 to 25.5 m from the other run's, against 0.2 km of 25.3 km
of the other class borders. Each route cleans the borders on its own at
20 m (`--features-tolerance 20`): the two-stage run cleaned a 1 m mesh's
classes, the direct run cleans CORINE with NVE's lakes. Each keeps its
lines within 20 m of what it cleaned, so the two can lie up to 40 m apart;
the farthest found is 25.5 m. So the lake lines differ by more than a few
metres along 1.7 km. Against the source both routes clean towards, NVE's
lake border inside the domain (49.63 km):

| Lake lines, 5 m | Within 1 / 5 / 20 m of NVE's lake border | Farthest from it | NVE's border within 20 m of the run's lake lines | Farthest |
|---|---|---|---|---|
| two-stage | 22.0 / 71.9 / 100 % | 18.8 m | 100 % | 19.9 m |
| direct | 21.9 / 71.8 / 100 % | 18.8 m | 100 % | 19.7 m |

Both routes keep the lake borders within 20 m of NVE's, and equally
close: neither is the more faithful. That 28 % of each run's lake border
lies 5 to 19 m from NVE's is the 20 m clean-up Ola's Sula commands ask
for, in both routes alike; a smaller `--edges-tolerance` keeps it closer.
The edge lines are inputs, the same at every tolerance: the 20 m pair
gives the same figures (4,943 and 3,271 triangles, round 3). The
triangle counts differ by the start mesh's smallest angle (25 degrees in
the two-stage run, 15 in the direct), not by the borders; the direct run
takes 7.3 s and no first stage. Planted check: CORINE's own lakes (Ola's
27 October `diagnostic-corine-only.gpkg`, 13.48 km of lake line) against
the direct run: of the direct run's lake lines, 11.2 % lie within 5 m of
CORINE's, the farthest 2,234 m. So the fed-back route gets no home.

*Read probe* (scratch script, not kept: `read_class_polygons` with the
`corine` rules on the quick check's two files, the read box from the domain
outline; tree `f4c7be7a`). Numedalslågen, 7,932 polygons: 3.8 s per read,
of which 90 % is `is_valid` in `_repaired` (cProfile: 17.2 of 19.1 s over
three reads), none repaired. Today's `features read` phase on that case is
under 0.2 s (stored quick-check baseline). Lagan, 9,173 polygons: 18 to
30 s per read, timed while the machine's load average stood at 70, so an
order of size only. Reading the map through `read_class_polygons` as it
stands would have cost most of a border-following run's time again; hence
the read without the check (section 3).

## 5. Tests

Fixture: `test_cli_mesh_prefer.py`'s scene (`cover_polygons`: 311, 512, 312
with a hole; `preferred_lake`), `--edges-repair 0 --edges-tolerance 0
--edges-outline-snap 0`. Measured with `preferred.overlay` and shapely: the
map holds 311 1,409,620 m², 312 509,495 m², 512 420,090 m². Every clean-up
switch is 0, so cleaning only the chosen classes (section 3) leaves the
line figures below as round 2 measured them.

1. `build`: those areas, to 1 m²; 110's join test (unit square) moved to it.
2. Line length inside the domain, off the outline, ±0.01 m: `--edges all`
   9,082.7 m; `--edges lake` 2,600.4 m; `--edges forest` 7,482.8 m;
   `--edges 311` 5,600.0 m; none without `--edges` (the run gives shares).
   With `all`, the `water` lines total 2,600.4 m.
3. A river (`(200.3, −2500.7)`, `(1350.9, −1100.2)`, `(2900.6, −300.1)`
   from `dem_origin`): `--land-cover ... --lines river=... --edges river
   --shares`, the `river` lines total 2,448.4 m, and the shares equal the
   run without `--lines`; `--lines river=... --edges lake` is refused,
   naming `--lines`.
4. Case 7's shares (existing test) unchanged under `--edges all --shares`.
5. Each refusal of section 3, exit 2, naming its option (the land-cover class file's
   category named `all` or `river`: naming the category); `--features` is
   an unknown option; a sources file with a `rules` key is refused.
6. MCP: the `mesh` tool lists `edges`, `lines`, `shares` and
   `land_cover_classes`, and no `features` argument.
7. Sula (data-dependent): the command of section 0 gives 3,271 triangles,
   area 30.4411 km² and lake 8.935 %, as Sol's record.
8. No land-cover map: test 3's river, `--lines river=... --edges river` and
   `--edges all`, no `--land-cover`: the `river` lines total 2,448.4 m, no
   other edge, no code and no share per triangle.
9. When the shares are computed: `--land-cover ...` alone writes the
   share table, as today (case 7's command, unchanged, gives case 7's
   shares). `--land-cover ... --edges all` writes no share table, its
   `--stats` has no `land cover fractions` phase, and its code per triangle
   equals the codes today's tree gives for `--features` with the scene's
   class file (pinned from today's tree before the change). Adding
   `--shares` writes the share table, and the codes are its dominant
   classes. `--reservoirs` with `--edges all` and no `--shares` is refused,
   naming `--shares`; `--shares` without `--land-cover` is refused.
   `--edges river` alone (test 8's run) gives no share table.
10. Several `--land-cover` files (the scene's polygons split in two files;
    `P` a planted box of class 211, arable land, which `--prefer lake=`
    does not replace (it replaces 512 only, the class file's
    `[cells.prefer]`), inside forest 311, from `dem_origin`
    `(420.3, −1319.8)` to `(620.4, −1119.1)`: 40,160.1 m², perimeter
    801.6 m, 379.8 m west of the preferred lake). The three runs below
    have `--prefer lake=...`. Measured on tree `12616b91` (scratch script, not
    kept): the files joined first given wins with the overlap rule of
    section 3, then `preferred.replace` (near the read region; no piece
    refilled), the areas by shapely; the line lengths from today's CLI,
    the built map as one `--features` file with `corine` and the clean-up
    switches above, inside the domain, and equal to the map's borders
    inside the domain by shapely. Planted check: `P` against forest 311
    gives no pair under `"overlaps"` and one, of 40,160.1 m², under the
    rule of section 3.
    - `a.geojson` (311, 312) then `b.geojson` (512, `P`), `--edges all
      --shares`: `P` is cut away whole, so the map is test 1's (311
      1,409,620 m², 312 509,495 m², 512 420,090 m²; ±1 m²); lines
      9,082.7 m, `water` lines 2,600.4 m (test 2's), shares equal the
      one-file run's (case 7) to 1e-6; `land_cover_overlap` gives
      40,160 m² (±1) and names `a.geojson`; `land_cover_source` names both
      files.
    - The same files, `b` first: `build` gives 311 1,369,460 m², 211
      40,160 m², 512 420,090 m², 312 509,495 m² (±1 m²); `--edges all`
      lines 9,884.3 m (`P`'s border added), `water` 2,600.4 m (unchanged:
      211 is no water class).
    - 311 cut at `x = x0 + 750.5` into a west half (with 312, in one file)
      and an east half (with 512, in the other), no `P`: `--edges all
      --no-edges-merge-same-class` lines 10,782.1 m (the 1,699.4 m cut
      added); `--edges-merge-same-class` 9,082.7 m: the county-line case.
    - Refused, exit 2, naming `--land-cover`: a class raster (the scene's
      catalogue key) with `cover.geojson`; two catalogue entries whose
      `classes` differ.
11. A cut run (`--pieces 4`, `.vtk` out, the scene and clean-up switches
    above; the `mesh` helper's `--tolerance 3`). Pinned from today's tree
    (`75ab2b4f`'s code) with `--features cover.geojson --features-rules
    corine`, the same switches and `--pieces 4`: pieces `0-0-0`, `0-1-0`,
    `1-0-0`, `1-1-0` with 447, 484, 437 and 438 triangles; codes per piece
    {0: 171, 311: 276}, {0: 198, 311: 15, 312: 87, 512: 184}, {0: 159,
    311: 278}, {0: 258, 311: 14, 512: 166}; triangle area by code over all
    pieces 311 1,529,799.8 m², 512 910,920.4 m², 312 144,160.4 m² (the
    polygons' areas inside the domain: 1,529,799.9, 910,920.2, 144,160.6).
    - `--land-cover cover.geojson --land-cover-classes corine --edges all
      --pieces 4` gives those four pieces' triangles and codes exactly.
    - With `--prefer lake=lake.geojson` added: area by code over all
      pieces 311 1,409,620 m², 512 420,090 m², 312 144,161 m², each ±1 m²
      (the map of test 1, inside the domain).
    - Refused, exit 2, at the cut-run check of section 3 (after the DEM's
      headers are read and the cut planned, before any block is decoded
      or piece written; no `.vtk` file is left): `--land-cover cover.geojson
      --pieces 4` (no class borders, so shares), naming `--land-cover`;
      `--edges all --shares --pieces 4`, naming `--shares`.

## 6. Size, timing, left out

Net production lines: negative (about 300 deleted beyond what is added:
round 2's 400 less `--lines`, `--shares` and the kept `label_triangles`
route; estimated, not counted). Tests change widely (about 270 `--features`
invocations), all mechanical.

*Timing, border-following runs.* The quick check's stored baseline
(`docs/benchmarks/quick/baseline-apple-m1-max-ac.json`, tree `d886847d`,
medians in seconds) and round 2's run of today's code (tree `af82b7d3`):

| Case, today's command | Total, stored baseline | Total, round 2 | Read, round 2 | Clean-up, round 2 | Codes, round 2 |
|---|---|---|---|---|---|
| numedalslagen | 8.83 | 9.45 | under 0.2 (stored baseline) | 5.94 | 0.30 |
| lagan | 14.53 | 15.30 | 0.93 | 5.39 | 0.21 |
| numedalslagen, `corine-water` | (no case) | 3.40 | | 0.74 | 0.27 |

The new commands (`--land-cover ... --edges all`, so no shares) run the
same steps on the same polygons: the read as `_Tally._take` does it (one
`read_source` query, one move), the clean-up of every class (`all` chooses
every class, as `corine` does), and the codes by `label_triangles`. `build`
with no `--prefer` passes the base through. No new term is estimated above
noise, so both cases are expected at today's time: about 9 s on
Numedalslågen and 15 s on Lagan, within the quick check's band. `--edges
water` costs what `corine-water` costs today (3.4 s). The quick check's
cases stay at about 2 minutes.

*Gate.* `tools/bench_quick.py run` on the branch, on AC power: no `SLOWER`
verdict for numedalslagen or lagan. The cases keep their names and inputs,
so the stored baseline judges them (`_mismatch` compares inputs, not
arguments); phases renamed from `features ...` are not judged by name, the
total is. No `MESH CHANGED` verdict is expected either: the lines are
today's (section 4, step 3), so the meshes are. `@developer` runs it and
names both verdicts in the handback.

*Timing, several files.* One new term, only with a second `--land-cover`
file: the overlap query with its intersection areas, 0.23 s on Bygdin
(1,697 and 247 polygons, `read_source` over each file's `source_region`
box; 345 intersecting pairs, none sharing area; scratch probe on tree
`12616b91`, AC power). Cleaning the
two files as one coverage costs what cleaning each costs (`coverage_clean`
over the read region's 591 and 12 polygons: 0.249 s per file, 0.245 s
joined, best of three). One file adds nothing.

*Timing, the shares.* What `--land-cover` costs today, in the runs that
compute them (section 3): round 2's run measured 22.8 s of shares on Numedalslågen and 52 s on Lagan
(653,039 triangles on Numedalslågen), the polygon read with its validity
check included.

No `@perf` acceptance run (no refine or mesh code). Left out: fetching line
layers from the catalogue (the ideas entry above); shares in cut runs;
faster class shares.

## 7. Questions for Ola

Answered, 2026-10-28: (1) the sources-file key `rules` becomes `classes`
(section 2). (2) Rivers, and roads the same way, without a land-cover map:
`--lines` (section 3). (3) No slowdown for a run that follows borders: it
computes no shares unless `--shares` asks (sections 3 and 6).

Open:

4. *A county line inside a catchment.* Bygdin's two AR50 files meet on a
   7.6 km county line inside the catchment, with the same class on both
   sides all along it. Today's mesh follows that line (each file is
   cleaned on its own); with the files read as one map it does not, since
   the line is no border on the ground. Default: it is not followed. Say
   if county lines should stay edges.

Settled without asking, by Ola's 732 recorded `rasputin mesh` commands
(every distinct command in a JSON file under `~/rasputin-output`),
counting the sources given to each option the design changes:

| Option | Runs | Given twice | What the second source is |
|---|---|---|---|
| `--features` | 713 | 17 | 9: NVE's lakes beside CORINE (Sula), now `--prefer`; 5: the same mesh's `lines` layer under two rule files (Sula, fed back), no route; 3: a second AR50 county file (Bygdin twice, Øyangen), now a second `--land-cover` |
| `--features-rules`, `--features-legend` | 13, 700 | 8, 9 | the rule or legend of those second files |
| `--features-layer` | 8 | 5 | `lines` twice, the fed-back runs |
| `--land-cover`, `--prefer` | 29, 18 | 0, 0 | |
| `--features-crs`, `--pieces`, `--piece-budget`, `--catchments` | 0 | | |

No run gives any option three times. The `--features-*` clean-up options
(706 runs each), `--land-cover-year` (2), `--reservoirs` (7) and
`--simplify-domain` (9) are given once where given.

- What `--features` took and the new options do not: a line file whose
  `property` column names each row's kind, and polygons kept only for
  their outline (`property = wall`). No recorded run uses `property` or
  `--features` without rules; they go. The one route that runs did use, a
  mesh's own `lines` layer fed back (Sula, 27 October), gives the borders
  the direct command gives (section 4, two-stage probe).
- `--land-cover` alone keeps meaning the shares, so the 12 recorded runs
  that pass it without `--features` keep their commands (section 3), and
  nothing is refused as unused but a `--lines` key no `--edges` takes.

## 8. Python interface, and the wordings tests read

Added after design review round 6, for `@tester`'s migration of the
existing tests (`Migration (111)` commits, up to `d54e7cb8`). Interface
only: section 3's behaviour stands. Where this section and section 3's
diagram differ (`build` never returns None; the map's `values`), this
section is the pin.

*Names.* In `tin_engine.feature_input` (which keeps `read_source`,
`SourceRows`, `source_region`, `read_lake_polygons`, `read_polygons`,
`split_key`, `Split`, `FeatureError`, `FeatureSet`, `TerrainFeature`'s
import, `snap_to_outline`, `OutlineSnap`, `IN_LINE`, `REST_MIN_WIDTH_M`,
`_untouched`, `_outline_lines`, `_chains`, and the module name
`simplify_borders` that tests patch):

```python
class LandCoverClasses(BaseModel):  # frozen, extra="forbid": MeshingRules less what goes
    name: str
    attribute: str
    codes: str                         # the code system; required (no uncoded file)
    legend: Legend                     # required with a code system, as `load` checks today
    classes: frozenset[str]            # the listed values: MeshingRules.classes' keys
    otherwise: Literal["refuse", "drop", "keep"] = "refuse"   # "keep": the value is the code
    notice: str = ""; licence: str = ""; accuracy_m: float | None = None
    string_codes: Mapping[str, int] = {}; splits: tuple[Split, ...] = ()
    categories: Mapping[int, Category] = {}; water_codes: frozenset[int] = frozenset()
    prefer: Mapping[Category, int] = {}; accept_no_class: bool = False

class EdgeChoice(BaseModel):         # frozen, extra="forbid"; the CLI resolves --edges into it
    every_class: bool = False          # `all`: every class in the map
    codes: frozenset[int] = frozenset()   # classes whose borders are edges
    lines: frozenset[str] = frozenset()   # --lines keys whose lines are edges

class EdgeCleanup(BaseModel):        # frozen, extra="forbid"; all off by default, as FeatureRequest
    repair_m: float = 0.0
    merge_same_class: bool = False
    tolerance_m: float = 0.0
    outline_snap_m: float = 0.0

@dataclass(frozen=True, slots=True)
class ClassPolygons:                 # what read_class_polygons read
    areas: tuple[tuple[int, BaseGeometry], ...]  # (code, polygon), file order, in mesh_crs, unchecked
    crs: str; layer: str | None; scanned: bool; empty: int   # empty: rows skipped as empty

@dataclass(frozen=True, slots=True)
class SourceLines:                   # what read_linestrings read
    lines: tuple[LineString, ...]      # in mesh_crs, file order, a MultiLineString as its parts
    crs: str; layer: str | None; scanned: bool; empty: int

def read_class_polygons(path: Path, layer: str | None, classes: LandCoverClasses, mesh_crs: str,
                        box_for: Callable[[str], tuple[float, float, float, float]],
                        *, crs: str | None = None, threads: int = 1) -> ClassPolygons
def read_linestrings(path: Path, layer: str | None, mesh_crs: str,
                     box_for: Callable[[str], tuple[float, float, float, float]],
                     *, threads: int = 1) -> SourceLines
def repaired_areas(areas: Sequence[tuple[int, BaseGeometry]], xy: npt.NDArray[np.float64],
                   tell: Callable[[int, int], None] | None = None,
                   ) -> tuple[list[tuple[int, BaseGeometry]], int]
def edges_of(land_cover: LandCoverMap | None, lines: Sequence[tuple[str, LineString]],
             choice: EdgeChoice, cleanup: EdgeCleanup, domain: DomainPolygon,
             vocabulary: EdgeVocabulary, *, region: Polygon) -> FeatureSet
```

In `tin_engine.land_cover_map`, as `test_land_cover_map.py` pins it:

```python
@dataclass(frozen=True)
class LandCoverMap:
    areas: tuple[tuple[int, BaseGeometry], ...]   # unchecked polygons, DEM's CRS
    classes: LandCoverClasses
    replaced: preferred.Replaced | None           # the --prefer cut; None without --prefer
    overlap: tuple[tuple[int, BaseGeometry], ...] # (index of the file kept, area cut from a later one)

def build(base: Sequence[Sequence[tuple[int, BaseGeometry]]],      # one per --land-cover file
          preferred: Sequence[tuple[int, Sequence[BaseGeometry]]],  # as preferred.replace's chosen
          classes: LandCoverClasses, region: Polygon) -> LandCoverMap
```

`replaced` and `overlap` are section 3's `values`: the CLI measures them
inside the domain for `land_cover_prefer` and `land_cover_overlap`, since
`build` sees the read region only. A class raster gives
`LandCoverMap((), classes, None, ())`; no `--land-cover` gives no map (the
CLI passes None).

- `box_for` is `read_source`'s: the CLI passes `lambda own:
  source_region(domain, dem_crs, own, cell=cell).bounds`, as `_Tally.region`
  does today. `region` in `edges_of` is that region in the DEM's CRS.
- `read_class_polygons` reads as `_Tally._take`'s cleaned branch: one query,
  every polygon moved in one `shapely.transform` call, no validity check.
  Refused (`FeatureError` naming file and feature): a value no class takes
  under `otherwise`; a vertex with no image in the DEM ("has a vertex with
  no image in the DEM", as today); a point, a collection, and a line (lines
  go to `--lines`). The CLI's hint is `--land-cover`.
- `read_linestrings` refuses a polygon row (the CLI's hint `--lines`; the
  message names `--prefer`), a point, a collection, and a vertex with no
  image. It takes no `crs`: no `--lines-crs` (section 3).
- `repaired_areas` is today's `read_class_polygons` repair half
  (`_repaired` near `xy`, emptied polygons dropped, the count returned and
  told), run in the shares step on `map.areas`; its stderr line stays
  `_say_repaired("land cover", ...)`'s.
- `edges_of` is pure. A chosen area's mask is `land_cover`, plus `water`
  when its code is in `water_codes`; its `TerrainFeature` keeps `code` and
  its cleaned polygon (for `label_triangles`), `fid` its index in
  `land_cover.areas` (a merged group, its first). A line's mask is its key's
  bit, `fid` `f"{key} {i}"` with `i` its index in `lines`, `code` None.
  Features: the areas first, then the lines in order. The clean-up runs, as
  `open_features` decides today, when a switch is on or the map was cut
  (`replaced` not None), with snapping `repair_m or -1.0` after a cut; with
  every switch off and no cut it does not run and `cover_vertices` is None,
  so test 11's pinned meshes hold. `FeatureSet` drops `crs`, `layers`,
  `scanned` and `clip_seconds` (the readers and the CLI's clock hold them);
  `counts` has one entry per source edges were taken from: the map, when
  the choice names a class, then each chosen `--lines` key in the order
  its first line appears.
- The CLI imports `edges_of` by name, so a spy patches `cli.edges_of`, as
  tests patch `cli.open_features` today.
- `land_cover_classes` (module renamed from `meshing_rules`) keeps `load`,
  `shipped_names`, `known_codes`, `derived_cleanup` and `CleanUpDistances`
  (its fields in `EdgeCleanup`'s names). `load("property")` is an unknown
  name. A class file without `code_system` is refused naming
  `code_system`.

*What replaces each removed name, for the unstarted files.*

| Removed | Becomes |
|---|---|
| `open_features` on polygons | `read_class_polygons` → `build` → `edges_of` |
| `open_features` on lines | `read_linestrings` → `edges_of(None, ...)` |
| `FeatureSource` | the reader's arguments (path, layer, classes, crs) |
| `FeatureRequest` | `EdgeCleanup`, plus `EdgeChoice`; `threads` goes to the readers |
| `MeshingRules` | `LandCoverClasses` |
| `Prefer`, `open_features(prefer=)` | `build`'s `preferred` |
| `_Tally`, `pre_clip`, `_runs`, `_linework`, `_bound`, `R_MIN`, `R_MAX`, `_mask` | go, with their tests |
| `meshing_rules.run_vocabulary`, `FIRST_EXTRA_BIT`, `MAX_EXTRAS`, `NO_EDGES` | go: `DEFAULT_VOCABULARY`, plus `catchment` with `--catchments` |
| `preferred.overlay`, `Overlaid`, `cli._overlay_preferred` | go; `replace` stays |

Per file:

- `test_feature_input`: `TestSources`, `TestOneGeojsonRule` and the
  `read_source` cases stay. `property` tests, `the_request_vocabulary`,
  name-per-value maps, `corine_water`'s, and the per-source layer/CRS
  pairing go. CORINE and `clc18_kode` masks: `edges_of` with
  `EdgeChoice(every_class=True)` (water classes `land_cover|water`, the
  rest `land_cover`). `TestGeometry`, `TestCrs`, `TestAnEmptyGeometryIsNull`:
  the two readers (an empty row counts in `empty`; a self-crossing ring is
  read, not refused). `TestClip`: `edges_of`.
- `test_legends`: `fi._mask`, `TestVocabularyExtras`, the uncoded base,
  `class.edges` overrides, `property` go; `t2` masks through `edges_of`;
  the rest load through `land_cover_classes`.
- `test_cli_legends`: `--features-rules` → `--land-cover-classes`;
  "palette refuses a rules file" → "a land-cover class file".
- `test_cli_mesh_features`: the gallery's lines → `--lines` per kind
  (`road`, `river`, `wall`), its polygons → `--land-cover` with
  `--edges all`. `--edges` without `--domain`: "needs --dem, --domain and
  --tolerance", hint `--edges`. The `--features-*` option tests go or
  rename (`--land-cover-layer`, `--land-cover-crs`, each refused without
  `--land-cover`); `FEATURES_FIELD` follows the `edges` row below.
- `test_cli_mesh_multi_features`: positional pairing, two maps over one
  file, and codeless-beside-coded go; two coded files are two
  `--land-cover` files under one class file (test 10); the shared seam
  is a land-cover border and a `--lines` line; two CORINE files give one
  `land_cover_notice` (de-duplicated, as today).
- `test_cleanup_every_source`, `test_feature_repair`,
  `test_cli_features_cleanup`: `request(...)` → `EdgeCleanup`, the read →
  `read_class_polygons`/`build`, the spy → `cli.edges_of`.
  `TestProperty` and `rp5` (an uncoded map) go: lines are never cleaned
  (section 3; test 8's river length checks it at default switches).
- `test_source_fits`: `read_class_polygons` with `land_cover_classes.load(fit.classes)`.
- `test_preferred`: the `overlay` tests go; `replace`'s stay.

*Wordings.*

- `--stats` row `edges`: one entry per `counts` source, joined by `; `:
  the map `<file> + <file>[ layer <L>], classes <C>: <n> features, <c>
  lines, <v> vertices`; a line key `<file>, lines <key>: <n> features,
  <c> lines, <v> vertices`. Today's form, with `rules <R>` replaced by
  `classes <C>` or `lines <key>`; the counts unchanged.
- stderr, when `edges_of` runs: `edges: <n> kept (<c> cut at the domain
  outline), <o> outside the domain, <e> empty`, `<e>` the readers' `empty`
  plus `FeatureSet.empty`. The no-index line is unchanged.
- Phases: `land cover read` (the `--land-cover` files and `--prefer`
  sources read, and `build`), `lines read` (the `--lines` files), `edges`
  (`edges_of`), `edges: clean-up` (its clean-up). `features read`,
  `features clip` and `features clip: clean-up` go. The two ruled
  entries of `docs/benchmarks/quick/hotspots.toml` (lagan, numedalslagen)
  change `phase` to `edges` in the PR, their rulings verbatim: it is the
  same clean-up, less the vertex move, now under `land cover read`.

*Pins checked* (from the migration commits): `--lines` as two files for
the crossing's road and river (one source per key); line lengths from the
`lines` layer to 0.01 m (test 2's bound); "meshing rules" becoming
"land-cover class file"; `--stats` row `features` becoming `edges`; a class
file without `code_system` refused naming it, hint `--land-cover-classes`;
`classes =` for `rules =`, on both shipped entries. All six confirmed.

## Review

- Design round 1: CHANGES REQUESTED, 75ab2b4f..a2f4dcaa, 0 (no speed estimate for the map and per-class edges; line files lose layer and CRS; Ola's sources file breaks on the rules → classes rename; six wordings and inventory gaps).
- Design round 2: APPROVED WITH WORDINGS, af82b7d3..fd566170, 0; question 3 reworded and question 4 added in this commit.
- Design round 3: CHANGES REQUESTED, fd566170..0da191d5, 0 (question 4's premise false: Ola's 27 Oct Sula 2 m and 5 m runs fed the 20 m mesh's `lines` layer back, which `--lines` refuses; question 5 narrowed to the `--land-cover` half, worded from the 26 Oct Swedish runs; `--lines-crs` and `--prefer-crs` have no user; three wordings).
- Design round 4: CHANGES REQUESTED, 132a89cf..498e940a, 0 (several land-cover files in one run dropped unsaid: bygdin-ar50-vtol5 and oyangen-ringkollen pass two AR50 county files; Sula two-stage probe's 25.7 m exceeds the stated 20 m band and the 2 % beyond 5 m is not located; no cut-run (`--pieces`) test, `--prefer` in cut runs unstated).
- Design round 5: CHANGES REQUESTED, 369bac65..bf9b88de, 0 (Sula and cut-run findings fixed; §3's overlap query `STRtree.query(predicate="overlaps")` returns nothing for a polygon wholly inside another, so test 10's planted box `P` inside forest 311 gives no overlap: query with "intersects" and keep pairs of shared area above zero, and re-measure Bygdin's 0 m² if it used "overlaps"; test 10 bullet 2 ignores that `--prefer lake=` replaces `P` (class 512) as 311, measured 311 1,409,620 m², 512 420,089.9 m², so the stated 1,369,460 m², 460,250 m² and 3,402.0 m water lines are wrong; wording: test 11's "before any phase" vs the refusal after DEM open and piece split, `cli.py:1186-1203`).
- Design round 6: APPROVED, 12616b91..3ea5ceb8, 0 (all three round-5 fixes hold: Bygdin re-probed, 345 intersecting pairs, none sharing area, 0.24 s; test 10 figures checked by hand, P 40,160.07 m², 311 1,369,459.9 m², lines 9,884.3 m; test 11 refusal after `_open_dem` and `partition`, before `assemble_dem` and `_open_features`).
