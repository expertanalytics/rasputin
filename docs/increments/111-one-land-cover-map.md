# Increment 111: one land-cover map, and the edges chosen from it

**Status: designed** by `@architect`, 2026-10-28, on master `75ab2b4f`,
worktree `land-cover-edges`. One PR, full loop, Python only. Number 111: 110
is the highest under `docs/increments/` over every ref and worktree.

Ola, 2026-10-28, on the split below: "Yes, this is logical. Now we need to
make it consitent :)"

## 0. The plain picture

Two questions, two sets of options:

1. **What is on the ground?** The land-cover map: `--land-cover` (CORINE's
   polygons, AR50, a class raster) and `--prefer KEY=SOURCE` (NVE's lakes
   and glaciers as areas, rivers and roads as lines). `--prefer` decides which
   source gives each class. The result is one map, built once, before the mesh.
2. **Which of those must the mesh follow?** The meshing rules: `--edges`,
   chosen per class (`--edges lake --edges river`; `--edges all`). Edges come
   only from the map. `--edges-*` shape them (today's `--features-*`).

The catchment outline (`--domain`) and sub-catchment borders (`--catchments`)
stay their own concern. `--features` goes: everything it carried is land cover.
Sula's direct command becomes

    rasputin mesh ... --land-cover clc2018 --prefer lake=nve-lakes --edges all
        --edges-repair 0.1 --edges-tolerance 20 --edges-outline-snap 0

with no rules file to upload.

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
| `mesh --features` (repeatable) | **deleted**: class polygons go to `--land-cover`, lakes, glaciers, rivers, roads to `--prefer` |
| `--features-crs`, `--features-layer` (paired, repeatable) | `--land-cover-crs`, `--land-cover-layer`: one value, the `--land-cover` file's |
| `--features-rules`, built-in `property` | **deleted**: what classes mean is `--land-cover-classes`, what the mesh follows `--edges` |
| `--features-repair`, `--[no-]features-merge-same-class`, `--features-tolerance`, `--features-outline-snap` | `--edges-repair`, `--[no-]edges-merge-same-class`, `--edges-tolerance`, `--edges-outline-snap`; same meanings; defaults from the class file's accuracy, as for a coded source today; the tolerance still sets the `--catchments` and `--simplify-domain` band |
| `--land-cover`, `--land-cover-year` | kept |
| `--land-cover-rules` | `--land-cover-classes` |
| `--prefer` | kept; also takes line sources (section 3); its sentence on `--features` lines deleted |
| (none) | `--edges SELECTOR`, repeatable (section 3) |
| `RENAMED`, `_RenamedOptions` (`--features-legend`, `--land-cover-legend`) | **deleted**: their targets go |
| meshing-rules files (`data/meshing_rules/`, `meshing_rules.py`, "meshing rules" in help and docs) | **class files**: `data/classes/`, `classes.py`, "class file" |
| class file keys `[edges]`, `class.edges`, `class.value`, `[vocabulary] extra` | **deleted** (edge names follow from the class, section 3) |
| `unlisted = "land"` | `unlisted = "keep"`; `drop`, `refuse` kept; `code_system` required |
| `water`, `[cells]`, `[[split]]`, `codes`, `accept_no_class`, `positional_accuracy_m`, `notice`, `licence`, `legend`, `name`, `attribute` | kept |
| `corine-water.toml` | **deleted**: `--edges water` |
| catalogue key `rules =` (shipped `sources.toml`, Ola's own file: 15 entries) | `classes =`; the old key is an unknown key |
| record `features`, `features_joined` | `edges`, `edges_joined` |
| record `features_crs`, `_transform`, `_accuracy_m`, `_notice`, `_licence` | `land_cover_crs`, `_transform`, `_accuracy_m`, `_notice`, `_licence` |
| record `features_repair_m`, `_merge_same_class`, `_tolerance_m`, `_outline_snap_m` | `edges_repair_m`, `edges_merge_same_class`, `edges_tolerance_m`, `edges_outline_snap_m` |
| record `land_cover_prefer` | kept, plus a row per line key ("river from X: N lines, L km") |
| MCP `mesh` arguments `features*`, `land_cover_rules` | follow the CLI (the schema is derived): `edges`, `edges_*`, `land_cover_classes`, `land_cover_crs`, `land_cover_layer` |
| MCP `INSTRUCTIONS`: "features_rules", "features_layer" | class polygons to `land_cover` with the fit's classes; lakes, glaciers, rivers, roads to `prefer`; `edges` |
| `sources` fits: `mesh --features` with rules | `mesh --land-cover` with classes (class polygons), `mesh --prefer` (polygons or lines); `Fit.rules` → `Fit.classes` |
| `FeatureSource`, `open_features`, `_Tally.source`, `_take`, `pre_clip`, `_runs`, `_bound`, `_linework`, `R_MIN`, `R_MAX` | **deleted**: reading is `read_class_polygons`, `read_polygons`, new `read_lines` |
| `FeatureRequest` | `Cleanup(repair_m, merge_same_class, tolerance_m, outline_snap_m)` |
| `MeshingRules` | `ClassRules`, less `classes` (value to edge names), `otherwise`'s names, `extra_names`, `cleaned` |
| `meshing_rules.run_vocabulary` | the default vocabulary, plus `catchment` with `--catchments` |
| #387: `cli._cover_lines` and its two refusals, `Prefer`, `open_features(prefer=)`, `_Tally._preferred` | **deleted**; the join moves into the map (section 3) |
| `preferred.overlay` | **deleted**; `preferred.replace` kept |
| `landcover.label_triangles`, `cli._land_cover` (codes from coded `--features`) | **deleted**: codes come from the fractions; `landcover.regions` kept |
| `gml-to-gpkg` help "for --features" | "for --land-cover" |
| `docs/benchmarks/quick/cases.toml` numedalslagen, lagan | `--land-cover <file> --land-cover-layer/-crs ... --land-cover-classes corine --edges all` |
| README, INSTALL, `project_structure.md`, ROADMAP | reworded |

## 3. Design

```
--land-cover (+ -classes, -layer, -crs)   --prefer KEY=SOURCE ...
        │ read near the read region, in the DEM's CRS (CLI, I/O)
        ▼
land_cover_map.build(base, preferred, classes, region) -> LandCoverMap        [pure]
        │ areas: (code, polygon)…   lines: (key, line)…            │
--edges, --edges-* → EdgeChoice, Cleanup                           │
        ▼                                                          │
feature_input.edges_of(map, choice, cleanup, domain, vocabulary) -> FeatureSet [pure]
        │ + sub-catchment borders (unchanged) → mesh (C++ unchanged)│
        ▼                                                          ▼
                     after the mesh: assign_polygon_fractions(map.areas) → fractions, codes, cells
```

- **The map** (`tin_engine/land_cover_map.py`, frozen dataclass
  `LandCoverMap(areas, lines, classes, values)`). The base is
  `read_class_polygons` over the read region (`source_region`), or no areas
  for a class raster, read after the mesh as today. Each area `--prefer` is
  107's `replace` near the read region, first given wins, and each refilled
  piece is unioned with the same-class polygons it touches (110's join, moved
  from `_Tally._preferred`, on `(code, polygon)` pairs). The record's
  replaced, given and refilled areas are measured inside the domain. A line
  `--prefer` is `read_lines` (LineString parts, moved vertex by vertex) and
  becomes the map's lines of that key.
- **Keys.** An area key is one of the class file's `[cells.prefer]` (lake,
  glacier), as today. A line key is an edge-vocabulary name other than
  `land_cover` and `water` (river, road, railway, coastline, contour, wall,
  ditch). Refused, naming `--prefer`: a line source for an area key, a
  polygon source for a line key, a second source for one line key. A line
  key works with a class raster too.
- **`--edges`.** Selectors: `all`, `water` (the class file's water classes),
  a category of its `[cells.categories]` (lake, forest, glacier, ...), a class
  code, a line key. Resolved at the CLI into `EdgeChoice(codes, lines)`. A
  border between two areas is an edge when either side's class is chosen;
  it carries `land_cover`, and `water` when that side is a water class (the
  masks every shipped class file gives today). A chosen line class's lines
  carry its name. Default: no edges. Refused, naming `--edges`: an unknown
  selector (the message lists the accepted ones), an area selector with a
  class raster, any selector without `--land-cover`.
- **`edges_of`.** Every area goes into the clean-up (`_Tally._clean`) as one
  coverage, so the repair and the simplification see each border's real
  neighbour; only chosen classes' polygons give lines (`_add`). Lines are
  clipped to the domain, not cleaned, as `--features` lines are today.
- **Fractions** come from `map.areas`: the map the edges were drawn from.
- **A cut run** (`--pieces`) takes `--land-cover` only with `--edges`, for
  the edges; `--pieces`' help says pieces carry no class fractions yet.

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
| step 3, as designed | 3,271 | 14 | 2 of 1,201, at most 1.4 mm | 1.1e-11 |

Step 3 differs at two vertices only: its cut runs on the repaired polygons
(`_repaired` clips each self-crossing ring at the read box), today's line
copy on the unrepaired ones. Repairing today's map after the cut
gives today's mesh exactly, so the repair's order, not its effect on the map,
moves them; both maps give the same class areas inside the domain to 1 m².
Every run's cells total 30.4411 km² with 8.935 % lake.
Planted check: the same command with CORINE's own lakes gives 2,025
triangles, so the probe tells maps apart.

## 5. Tests

Fixture: `test_cli_mesh_prefer.py`'s scene (`cover_polygons`: 311, 512, 312
with a hole; `preferred_lake`), `--edges-repair 0 --edges-tolerance 0
--edges-outline-snap 0`. Measured with `preferred.overlay` and shapely: the
map holds 311 1,409,620 m², 312 509,495 m², 512 420,090 m².

1. `build`: those areas, to 1 m²; 110's join test (unit square) moved to it.
2. Line length inside the domain, off the outline, ±0.01 m: `--edges all`
   9,082.7 m; `--edges lake` 2,600.4 m; `--edges forest` 7,482.8 m;
   `--edges 311` 5,600.0 m; none without `--edges`. With `all`, the `water`
   lines total 2,600.4 m.
3. A river (`(200.3, −2500.7)`, `(1350.9, −1100.2)`, `(2900.6, −300.1)`
   from `dem_origin`): `--prefer river=... --edges river`, the `river` lines
   total 2,448.4 m; with `--edges lake` instead, no `river` line, and the
   fractions equal the run without the river.
4. Case 7's shares (existing test) unchanged under `--edges all`.
5. Each refusal of section 3, exit 2, naming its option; `--features` is an
   unknown option.
6. MCP: the `mesh` tool lists `edges` and `land_cover_classes` and no
   `features` argument.
7. Sula (data-dependent): the command of section 0 gives 3,271 triangles,
   area 30.4411 km² and lake 8.935 %, as Sol's record.

## 6. Size, timing, left out

Net production lines: negative (about 400 deleted beyond what is added;
estimated, not counted). Tests change widely (about 270 `--features`
invocations), all mechanical. Timing: the land cover is read and cut once,
not twice. No `@perf` acceptance run (no refine or mesh code); the quick
check's numedalslagen and lagan cases change command, so their stored
timings stop being a baseline and the next quick check starts one. Left out:
fetching line layers from the catalogue (the ideas entry above); fractions
in cut runs.
