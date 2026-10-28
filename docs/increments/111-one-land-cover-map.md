# Increment 111: one land-cover map, and the edges chosen from it

**Status: designed, revised after Ola's answers to questions 1 to 3 (lines without a map, shares only when asked); design review round 3 due** by `@architect`, 2026-10-28, on master `75ab2b4f`,
worktree `land-cover-edges`. One PR, full loop, Python only. Number 111: 110
is the highest under `docs/increments/` over every ref and worktree.

Ola, 2026-10-28, on the split below: "Yes, this is logical. Now we need to
make it consitent :)"

## 0. The plain picture

Three questions, three sets of options:

1. **What is on the ground?** Areas: the land-cover map, `--land-cover`
   (CORINE's polygons, AR50, a class raster), with `--prefer KEY=SOURCE`
   taking one class from another source (NVE's lakes and glaciers). Lines:
   `--lines KEY=SOURCE` (rivers, roads), with or without a land-cover map.
2. **Which of those must the mesh follow?** The edges: `--edges`,
   chosen per class or line kind (`--edges lake --edges river`; `--edges
   all`). Edges come only from the areas and lines above. `--edges-*` shape
   them (today's `--features-*`).
3. **What does the mesh file say about the ground?** `--shares`: each
   triangle's exact share of each class, and in a `.gpkg` the cells Shyft
   reads. Without it, a run that follows class borders writes a class code
   per triangle, as a `--features` run with CORINE does today, at today's
   cost. The shares are computed only when `--shares` asks.

The catchment outline (`--domain`) and sub-catchment borders (`--catchments`)
stay their own concern. `--features` goes. What it carried moves to the
areas and lines, except two inputs that stop working (question 4).
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
| `mesh --features` (repeatable) | **deleted**: class polygons go to `--land-cover`, lakes and glaciers to `--prefer`, rivers and roads to `--lines` |
| `--features-crs`, `--features-layer` (paired, repeatable) | `--land-cover-crs`, `--land-cover-layer`: one value, the `--land-cover` file's; `--prefer-crs KEY=CRS`, `--prefer-layer KEY=LAYER`, `--lines-crs KEY=CRS`, `--lines-layer KEY=LAYER`: per key, for a file (section 3) |
| `--features-rules`, built-in `property` (also `--features`' default, `cli.py:991`; `source_report.py:384`) | **deleted**: what classes mean is `--land-cover-classes`, what the mesh follows `--edges`. Lost with `property`: the two inputs of question 4 |
| `--features-repair`, `--[no-]features-merge-same-class`, `--features-tolerance`, `--features-outline-snap` | `--edges-repair`, `--[no-]edges-merge-same-class`, `--edges-tolerance`, `--edges-outline-snap`; same meanings; defaults from the land-cover class file's accuracy, as for a coded source today. `--edges-tolerance` also works without `--land-cover`: it sets the `--catchments` and `--simplify-domain` band, as `--features-tolerance` does today |
| `--land-cover`, `--land-cover-year` | kept; `--land-cover` now says what is on the ground, and computes no shares by itself (`--shares` does) |
| `--land-cover-rules` | `--land-cover-classes` |
| `--prefer` | kept, area keys only (lake, glacier), still needs `--land-cover`; its sentence on `--features` lines deleted |
| (none) | `--lines KEY=SOURCE`, repeatable: rivers, roads and the other line kinds, with or without `--land-cover` (section 3) |
| (none) | `--shares`: each triangle's class shares, the codes from them, and the `.gpkg` cells layer (section 3) |
| `--reservoirs`, `--catchments`: "needs --land-cover with meshing rules and a .gpkg output" (`cli.py:1689`) | "needs --shares, a land-cover class file with `[cells]` and a .gpkg output" |
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
| record `land_cover_prefer` | kept; a new `lines` field, one row per `--lines` key ("river from X: N lines, L km") |
| MCP `mesh` arguments `features*`, `land_cover_rules` | follow the CLI (the schema is derived): `edges`, `edges_*`, `lines`, `lines_*`, `shares`, `land_cover_classes`, `land_cover_crs`, `land_cover_layer` |
| MCP `INSTRUCTIONS`: "features_rules", "features_layer" | class polygons to `land_cover` with the fit's classes; lakes, glaciers to `prefer`; rivers, roads to `lines`; `edges`; `shares` for class shares and cells |
| `sources` fits: `mesh --features` with rules | `mesh --land-cover` with classes (class polygons), `mesh --prefer` (other polygons), `mesh --lines` (lines); `Fit.rules` → `Fit.classes` |
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
--land-cover (+ -classes, -layer, -crs)  --prefer KEY=SOURCE (+ -layer, -crs)   --lines KEY=SOURCE (+ -layer, -crs)
        │ read near the read region, in the DEM's CRS (CLI, I/O)                     │ read_linestrings
        ▼                                                                            ▼
land_cover_map.build(base, preferred, classes, region) -> LandCoverMap | None   lines: ((key, line), ...)
        │ areas: (code, polygon)...                                     │            │
--edges, --edges-* → EdgeChoice, EdgeCleanup                            │            │
        ▼                                                               │            │
feature_input.edges_of(map, lines, choice, cleanup, domain, vocabulary) -> FeatureSet   [pure]
        │ + sub-catchment borders (unchanged) → mesh (C++ unchanged)    │
        ▼                                                               ▼
after the mesh: --shares → assign_polygon_fractions(repaired map.areas) → shares, codes, cells
                no --shares → label_triangles(chosen classes' areas) → codes, as today
```

*Why round 2 tied the shares to the edges, and why nothing needs to.* Round
2 merged `--features`' class polygons (edges, a code per triangle) and
`--land-cover`'s (shares) into one `--land-cover`, and kept `--land-cover`'s
meaning: "the exact share of each class in every triangle". That was a
naming choice, not a constraint of the data: the edges are drawn from the
polygons before the mesh, the shares after it, and neither needs the other.
The one step they share, reading the polygons, does carry a cost trap
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
  a file has `--land-cover-layer`/`--land-cover-crs` (the `--land-cover`
  file), `--prefer-layer KEY=LAYER`/`--prefer-crs KEY=CRS` and
  `--lines-layer KEY=LAYER`/`--lines-crs KEY=CRS` (each key's file), so a
  river layer in a multi-layer GeoPackage and a GML file naming no `srsName`
  keep the route `--features-layer` and `--features-crs` gave them; both go
  to `read_source`'s `layer` and `crs` as today (`feature_input.py:656`).
  Refused, naming the option: either with a catalogue key, without its
  `--land-cover`, `--prefer KEY` or `--lines KEY`, a layer for a file that
  is not a GeoPackage, a CRS that is not the file's own (`read_source`
  refuses it today).
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
- **`--shares`** (flag). Needs `--land-cover`. The map's areas are repaired
  near the mesh (`_repaired`, today's `read_class_polygons` step, now
  here), then `assign_polygon_fractions` gives each triangle's shares; or
  the class raster's route, as today. The codes per triangle come from the
  shares, and a `.gpkg` with the class file's `[cells]` gets the cells
  layer, both as `--land-cover` gives them today. `--reservoirs` and
  `--catchments`' cell ids need `--shares` where they needed `--land-cover`.
- **Without `--shares`**, a run with area edges writes a code per triangle
  from `label_triangles` over the chosen classes' cleaned polygons (0 where
  none), as a coded `--features` run does today; with no area edges, no
  codes. No shares, no cells.
- **An input nothing uses** is refused, naming it and what would use it:
  `--land-cover` with neither an area selector in `--edges` nor `--shares`;
  a `--lines KEY` that no `--edges` selector takes (question 5).
- **A cut run** (`--pieces`) refuses `--shares` (pieces carry no class
  shares yet, as today); `--land-cover` with `--edges` and `--lines` work as
  in an uncut run, and each piece gets its codes as today.

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
   `--edges 311` 5,600.0 m; none without `--edges` (with `--shares`, so
   `--land-cover` is used). With `all`, the `water` lines total 2,600.4 m.
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
8. Test 3's river in a GeoPackage with a second layer, `--lines river=FILE
   --lines-layer river=rivers`, and in a GML 2 file (`gml:LineString`,
   `gml:coordinates`, as OGR writes; the 3.2 reader refuses curves) naming
   no `srsName`, `--lines-crs river=<the scene's CRS>`: the `river` lines
   total 2,448.4 m each (same geometry as test 3); without the option,
   refused.
9. No land-cover map: test 3's river, `--lines river=... --edges river` and
   `--edges all`, no `--land-cover`: the `river` lines total 2,448.4 m, no
   other edge, no code and no share per triangle.
10. Shares only when asked: `--land-cover ... --edges all` without
    `--shares` writes no share table, its `--stats` has no `land cover
    fractions` phase, and its code per triangle equals the codes today's
    tree gives for `--features` with the scene's class file (pinned from
    today's tree before the change). With `--shares`, the share table is
    written and the codes are its dominant classes. `--reservoirs` without
    `--shares` is refused.

## 6. Size, timing, left out

Net production lines: negative (about 300 deleted beyond what is added:
round 2's 400 less `--lines`, `--shares` and the kept `label_triangles`
route; estimated, not counted). Tests change widely (about 270 `--features`
invocations), all mechanical.

*Timing, border-following runs.* The quick check's stored baseline
(`docs/benchmarks/quick/baseline-apple-m1-max-ac.json`, tree `d886847d`,
medians in seconds) and round 2's run of today's code (tree `af82b7d3`):

| Case, today's command | Baseline total | Round 2 run | Read | Clean-up | Codes |
|---|---|---|---|---|---|
| numedalslagen | 8.83 | 9.45 | under 0.2 | 5.94 | 0.30 |
| lagan | 14.53 | 15.30 | 0.93 | 5.39 | 0.21 |
| numedalslagen, `corine-water` | | 3.40 | | 0.74 | 0.27 |

The new commands (`--land-cover ... --edges all`, no `--shares`) run the
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
total is. `@developer` runs it and names the verdict in the handback.

*Timing, `--shares`.* What `--land-cover` costs today, and only when asked:
round 2's run measured 22.8 s of shares on Numedalslågen and 52 s on Lagan
(653,039 triangles on Numedalslågen), the polygon read with its validity
check included.

No `@perf` acceptance run (no refine or mesh code). Left out: fetching line
layers from the catalogue (the ideas entry above); shares in cut runs;
faster class shares.

## 7. Questions for Ola

Answered, 2026-10-28: (1) the sources-file key `rules` becomes `classes`
(section 2). (2) Rivers, and roads the same way, without a land-cover map:
`--lines` (section 3). (3) No slowdown for a run that follows borders: the
shares are computed only with `--shares` (sections 3 and 6).

4. Two kinds of input that `--features` takes today stop working.
   (a) A line file whose rows each name their kind in a column called
   `property` (one row `river`, the next `road`): after the change, a line
   file gives one kind, `--lines river=FILE` makes every line in it a river,
   so such a file must first be split into one file per kind. Your answer
   to question 2 does not change this: rivers and roads each get their own
   `--lines`. (b) Polygons that are not land-cover classes but whose outline
   the mesh should follow (a dam or a building, tagged `property = wall`):
   no way to give them after the change. Neither is used by your
   `~/.config/rasputin/sources.toml` (its 15 `rules` entries are all
   `ar50`, class polygons; its one line source, `nve-main-rivers`, has no
   `property` column, so `--features` refuses it today, checked on a file
   with its columns; after the change it enters as `--lines
   river=nve-main-rivers` once cached), nor by the shipped catalogue, nor by
   the quick check (CORINE with the `corine` class file). Drop both?
   Default: yes.
5. An input that nothing uses: `--land-cover` with no `--shares` and no
   class border in `--edges`, or `--lines river=...` with no `--edges` that
   takes rivers. Today `--land-cover` alone means "compute the shares", so a
   command that passes `--land-cover` alone for the shares would now be
   refused with the message "add --shares for each triangle's class shares, or
   --edges for its borders". Refuse such inputs, naming the fix? Default:
   yes (the other choice is that `--land-cover` alone keeps meaning
   `--shares`, and unused lines are read and ignored).

## Review

- Design round 1: CHANGES REQUESTED, 75ab2b4f..a2f4dcaa, 0 (no speed estimate for the map and per-class edges; line files lose layer and CRS; Ola's sources file breaks on the rules → classes rename; six wordings and inventory gaps).
- Design round 2: APPROVED WITH WORDINGS, af82b7d3..fd566170, 0; question 3 reworded and question 4 added in this commit.
