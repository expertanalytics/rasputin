# Increment 111: one land-cover map, and the edges chosen from it

**Status: built** (red `bce16583`, green `a8246816`, code review round 1's fixes `6cd8bea3` and `77b3def6`; at `77b3def6`, 19 net production lines, 1,057 added and 1,038 removed, by `tools/count_loc.py 75ab2b4f 77b3def6`; full suite 11,232 passed, 18 skipped, 1 xfailed; README, INSTALL, GLOSSARY and `project_structure.md` reworded; next `@reviewer`'s code review round 2). Design rounds 8 and 9 (section 3, *Deciding on the input polygon*; test 13) amend round 7's polygon rule after code review round 5, before more code; round 9 (after design review round 8) checks `is_valid` on the polygon's parts near the read box before any cut, and sends a valid polygon with an edge through a corner of the box to `shapely.intersection`. Designed, approved in design review round 6, section 8 (the Python interface and the wordings tests read) added for the test migration, by `@architect`, 2026-10-28, on master `75ab2b4f`,
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
and fetches none. Design round 8 overlaps "Repair CORINE once per area of
interest": a cached, repaired extract would need no ring check. Not taken
up.

## 2. Inventory: what assumes the old way, and what it becomes

| Today | Becomes |
|---|---|
| `mesh --features` (repeatable) | **deleted**: class polygons go to `--land-cover` (repeatable, section 3), lakes and glaciers to `--prefer`, rivers and roads to `--lines` |
| `--features-crs`, `--features-layer` (paired, repeatable) | `--land-cover-crs`, `--land-cover-layer`: one value each, for every `--land-cover` file. No layer or CRS option for `--prefer` or `--lines` files: none of the 732 recorded runs gives one (section 3) |
| `--features-rules`, built-in `property` (also `--features`' default, `src_python/tin_engine/cli.py@75ab2b4f:991`; `src_python/tin_engine/source_report.py@75ab2b4f:384`) | **deleted**: what classes mean is `--land-cover-classes`, what the mesh follows `--edges`. No recorded run uses `property` or `--features` without rules (732 commands) |
| a mesh's own `lines` layer fed back (`--features <mesh>.gpkg --features-layer lines` with rules on its `land_cover` and `water` columns; Ola's Sula runs of 27 October at 1, 2, 5, 10 and 20 m) | **no route**: the direct command gives the same borders and the same shares (section 4, two-stage probe) |
| `--features-repair`, `--[no-]features-merge-same-class`, `--features-tolerance`, `--features-outline-snap` | `--edges-repair`, `--[no-]edges-merge-same-class`, `--edges-tolerance`, `--edges-outline-snap`; same meanings; defaults from the land-cover class file's accuracy, as for a coded source today. `--edges-tolerance` also works without `--land-cover`: it sets the `--catchments` and `--simplify-domain` band, as `--features-tolerance` does today |
| `--land-cover`, `--land-cover-year` | kept; `--land-cover`, now repeatable, says what is on the ground, and gives the shares unless the run follows class borders without `--shares` (section 3) |
| `--land-cover-rules` | `--land-cover-classes` |
| `--prefer` | kept, area keys only (lake, glacier), still needs `--land-cover`; its sentence on `--features` lines deleted |
| (none) | `--lines KEY=SOURCE`, repeatable: rivers, roads and the other line kinds, with or without `--land-cover` (section 3) |
| (none) | `--shares`: the shares as well, in a run that follows class borders (section 3) |
| `--reservoirs`, `--catchments`: "needs --land-cover with meshing rules and a .gpkg output" (`src_python/tin_engine/cli.py@75ab2b4f:1689`) | unchanged without class borders in `--edges`; with them, also "needs --shares" |
| (none) | `--edges SELECTOR`, repeatable (section 3) |
| `RENAMED`, `_RenamedOptions` (`--features-legend`, `--land-cover-legend`) | **deleted**: their targets go |
| meshing-rules files (`data/meshing_rules/`, `meshing_rules.py`, "meshing rules" in help and docs) | **land-cover class files**: `data/land_cover_classes/`, `land_cover_classes.py` (beside `land_cover_fractions.py`), "land-cover class file" |
| land-cover class file keys `[edges]`, `class.edges`, `class.value`, `[vocabulary] extra` | **deleted** (edge names follow from the class, section 3) |
| `unlisted = "land"` | `unlisted = "keep"`; `drop`, `refuse` kept; `code_system` required |
| `water`, `[cells]`, `[[split]]`, `codes`, `accept_no_class`, `positional_accuracy_m`, `notice`, `licence`, `legend`, `name`, `attribute` | kept |
| `corine-water.toml` | **deleted**: `--edges water` |
| catalogue key `rules =` (shipped `sources.toml`: 2 entries; Ola's `~/.config/rasputin/sources.toml`: 15, all `"ar50"`) | `classes =` (Ola, question 1: rename); `rules` becomes an unknown key, refused (`extra="forbid"`, `src_python/tin_engine/sources.py@75ab2b4f:106`). The shipped file changes in the PR. Ola's own file is outside the tree: the main session edits it (`rules =` to `classes =`, 15 lines, a `.bak` copy first) with Ola's yes, right after the merge and the main checkout's fast-forward, before any command or MCP server reads it again; the PR description says so |
| `src_python/tin_engine/sources.py@75ab2b4f:152`: a `legend` key refused, naming `rules` | **deleted**: an unknown key, refused like any other |
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
- **What every polygon in the map is** (code review round 4; this rule
  replaces the "unchecked" of the bullet above). The *read box* is the read
  region's bounding box widened by `MARGIN` (100 m). After `build` the map
  holds two sets, as two fields of `LandCoverMap`:
  - `areas`: the polygons that reach the read box, each a valid, non-empty
    Polygon or MultiPolygon, inside it except a preferred polygon kept whole
    (`Replaced.given`). Which of the two forms a polygon takes is decided on
    the input polygon (design round 8, *Deciding on the input polygon*,
    below): its cut to the box (`shapely.clip_by_rect`), or its rebuild
    (`_repaired`'s construction: its rings cut to the box, faces kept by
    parity, valid by construction). Checked and rebuilt in one place only,
    `_fixed` in `build` with `fix` on, in the thread pool, as now.
  - `outside`: the polygons whose bounding box misses the read box, and
    those whose cut is empty: as read, whole, never checked, never repaired,
    never overlaid. Only the edges step reads them, and only by their rings
    (below).

  `build` with `fix` off (the `--prefer` build in the CLI) takes its base as
  `areas` already so and is given the base map's `outside` (a keyword,
  default empty), which it passes through. The preferred polygons are valid
  as `read_polygons` gives them (`_repaired` each). GEOS's overlays of valid
  polygons give valid polygons (the OverlayNG documentation: results are
  always valid; not probed here), so every polygon `_joined`, `replace` and
  the refill join make is valid too. So any step after `build` may overlay
  `areas` without a check, and no step overlays `outside`. A test asserts
  `shapely.is_valid` on every polygon of `areas` (section 5, test 12).

  *Outside polygons, counted and kept.* Without clean-up (the pre-clip
  path of `edges_of`), each chosen class's outside polygon goes through
  today's ring test (`region.intersects(g.boundary)`, then `pre_clip`), so a
  hole lying in the region whose shell misses it (`test_feature_input.py`'s
  p5 cases) still gives its ring, exactly as given. Its feature carries no
  polygon: it labels no triangle, since by its cut it has no area in the
  box. One with no ring kept is counted outside. With clean-up on, or with
  the map cut by `--prefer`, each chosen class's outside polygon is counted
  outside, with no overlay: the clean-up works on area inside the read
  region, and an outside polygon has none there. The outside polygons come
  after the `areas` in the edges step's feature order. A valid outside
  polygon gives no line, so the meshes of valid inputs do not change. The
  shares and the codes never see `outside`.

  *Each step that met unchecked polygons, and why it is now safe:*
  1. `_joined` (the overlap of several files): sees `areas` only. The `far`
     index set, and the bookkeeping that excluded it, are deleted.
  2. `preferred.replace` and its `_cut` (the round-4 crash,
     `src_python/tin_engine/preferred.py@435d407b:70`):
     the base is `areas`, and the cutters are preferred polygons, both
     valid. A lake reaching an outside bowtie cuts nothing outside.
  3. `preferred._neighbours`' `subdivide`: the kept polygons come from
     `_cut` of valid polygons.
  4. `_refilled_joined`'s `union_all`: kept polygons and refilled pieces,
     all from valid overlays.
  5. `areas_inside` (from `_prefer_record`): `land.base.areas`, valid.
  6. The edges step's ring test: reads `outside` by rings and boundaries
     only. A predicate on a ring is safe on any ring.
  7. `_overlap_record`: the overlap pieces are intersections of valid pairs.
  8. The edges step's clean-up, `_Tally.clean`'s `shapely.intersection(g,
     region)` (`src_python/tin_engine/feature_input.py@435d407b:362`), not on `@developer`'s list: today
     p5's L-shaped polygon with `CLEANUP` raises GEOSException
     ("TopologyException: unable to assign free hole to a shell").
     Probed on tree `435d407b`. It now gets `areas` only.
  9. The shares, `assign_polygon_fractions(map.areas)` in `cli.py`: today it
     receives the outside polygons unchecked. Now it gets `areas` only.
  10. `label_triangles` over the pre-clip path's polygons: today p5's
      invalid polygon is a labelling polygon. Now an outside feature has
      none.

  *Why the outside polygons are neither checked whole nor repaired*
  (`@developer`'s proposal: check each outside polygon whole, repair the
  invalid ones with `make_valid`). Measured on tree `435d407b`, AC power,
  one thread, best of three: each case read as the CLI reads it, outside
  being the polygons whose cut to the read box is empty.

  | Input | Polygons read | Outside | Their vertices | Outside invalid | Check outside whole | Check the cuts |
  |---|---|---|---|---|---|---|
  | Bygdin, Innlandet AR50 | 1,697 | 244 | 781,220 | 16 | 145 ms | 48 ms |
  | Bygdin, ar50_46 | 247 | 113 | 484,609 | 5 | 92 ms | 2 ms |
  | Sula, clc2018 | 184 | 156 | 5,508,059 | 52 | 1,351 ms | 0.5 ms |
  | Numedalslågen, CORINE | 7,932 | 695 | 3,361,282 | 0 | 2,503 ms | 238 ms |

  `make_valid` of Sula's 52 invalid outside polygons takes 32.8 s. Checking
  them whole would add 2.5 s of work to the quick check's Numedalslågen
  case (about 9 s today): at least 0.58 s on any number of threads, since
  one polygon cannot be split. On Sula it would add 34 s of work (1.35 s
  of checks, 32.8 s of repair). Repair also
  breaks p5. `make_valid` gives the hole back as a shell, clockwise and
  from another vertex, with method `linework` and with `structure`. With
  the proposal patched into `_fixed` in a scratch plugin, p5's first case
  fails under both methods.

  The anomaly, settled: both timings hold, for different polygons. The
  cost of `is_valid` is set by holes, and by whether the polygon is valid.
  A valid polygon is checked through: CORINE's largest near Numedalslågen,
  409,242 vertices and 3,195 holes, takes 578 ms. An invalid one stops at
  its first fault: the same polygon after Sula's move from EPSG:3035,
  409,239 vertices, 3,180 holes, takes 63 ms. A ring without holes is fast:
  a 50,001-vertex circle takes 0.23 ms, the kind of polygon behind
  `@developer`'s 0.5 ms. So `_fixed`'s docstring holds on real inputs. Its
  sentence about the edges step is reworded to this rule.

  *Deciding on the input polygon* (design round 8, after code review
  round 5). Round 7 cut first and checked the cut. But `clip_by_rect` is
  only right on a valid polygon. GEOS cuts each ring at the box and closes
  the pieces along the box's edge, in a direction taken from the whole
  ring's orientation (`Orientation::isCCW` on the shell and on each hole,
  `src/operation/intersection/RectangleIntersection.cpp` on GEOS `main`).
  A ring that crosses or runs back over itself breaks that in two ways:
  - Its pieces in the box collapse to a 3-point ring, and the call raises
    (code review round 5; `@tester`'s 10 cases in
    `TestACutToAThreePointRing`, `8c5b4e2a`).
  - A stretch of it is wound the other way from the whole ring (a
    clockwise spike on a counter-clockwise square), so the pieces close the
    wrong way round. The cut comes out valid but inside out.

  The inside-out case: the square `EAST` with a spike to `(150, 50)` and
  back along `y = 50 + w`, which ends by running back over the square's west
  side. With `w` from 1e-9 to 10 m, its cut covers 488,639 to 490,000 m²,
  nearly the whole read box. `_repaired`, and `make_valid` then the cut,
  both give 136.11·`w` m² for `w` from 1e-6 to 10 m (1,361.1 m² at
  `w = 10`); at 1e-9 m both give under 2e-7 m². Master repaired every
  polygon on read, so it gave the right area. The cut-first order came in
  at `a4bcab1e`, as a speed fix.

  Design round 8 checked each ring near the box with
  `shapely.is_simple` and then cut the input. Design review round 8 broke
  it: two simple rings can still make the cut wrong (a hole poking out
  through its shell across the box's edge gives a valid cut of 220,000 m²
  where the rebuild keeps nothing), an empty cut can hide a rebuild that is
  not empty, and the cut raised on simple rings in 50 of its random cases.
  So validity is now decided on the polygon, before any cut, but only on
  its rings near the box.

  The rule (design round 9), per polygon reaching the read box, in
  `_fixed`:
  1. *Reach.* A polygon whose bounding box misses the read box (closed
     comparison) goes to `outside`, unchecked: no ring of it can enter the
     box or enclose it.
  2. *The near subset.* Of the input keep each part whose shell's bounding
     box meets the read box, each with only its holes whose bounding box
     meets it. These parts make one MultiPolygon, the *subset* (built with
     `shapely.polygons` and `shapely.multipolygons` over the group's rings,
     so it runs in the pool). `shapely.is_valid` is called on the subset,
     never on the whole input.
  3. *Subset invalid* (a ring near the box crosses or touches itself or
     another, a hole outside its shell, two parts overlapping): the polygon
     is rebuilt (`_faces`, step 6), and `clip_by_rect` never sees it. A
     rebuild that keeps something goes to `areas`, counted. One that keeps
     nothing goes to `outside`, as given and not counted, when every ring of
     the subset is simple (`shapely.is_simple`, on these polygons only):
     p5's L-shaped shell, whose hole lies outside it, stays as round 7 has
     it. Otherwise it is dropped and counted, as master's repair on read
     drops it (the 8 cases of `NO_AREA_IN_THE_BOX`).
  4. *Subset valid, and a corner of the read box on its boundary*
     (`shapely.intersects` of the subset's boundary with the four corners):
     the area is `shapely.intersection(subset, read box)`, its polygonal
     part. GEOS's cut is wrong there (below).
  5. *Subset valid otherwise:* the subset, not the input, is cut. An empty
     cut goes to `outside`, the input as given. A valid cut goes to
     `areas`. An invalid cut is replaced by `shapely.intersection(subset,
     read box)`, as in step 4, and not counted: the input is valid near the
     box. This is a valid polygon whose cut GEOS leaves with a
     self-touching ring; Bygdin has two, where the polygon leaves the box
     and comes back to touch its edge. In steps 4 and 5 an intersection
     with no area goes to `outside`.
  6. *Rebuilt*: `_faces(g, xy)` is `_repaired`'s construction without its
     `is_valid` shortcut. `_repaired(g, xy)` becomes `g if g.is_valid else
     _faces(g, xy)`, unchanged for its other callers. The map's `repaired`
     count is the number of polygons counted in step 3. The record's
     wording `land_cover_repaired` becomes "Land-cover polygons rebuilt
     near the mesh (invalid there)".

  `clip_by_rect` raising on a valid subset is not caught. None was found:
  0 of the 28,203 valid subsets of the searches below, and none on
  Numedalslågen, Sula or Bygdin. Catching it would hide a new kind of fault
  the way round 7's order hid this one. Steps 2 to 5 run in `check(ix)`, in
  the pool's 32 groups, as round 7's cut check does.

  *Why checking the subset is enough.* Inside the read box, the input's
  parity (what `_faces` keeps: a point in an odd number of parts, a part
  holding a point in its shell and in none of its holes) equals the
  subset's, whatever the left-out rings do:
  - A left-out ring's bounding box misses the read box, so no point of the
    box is inside it or on it.
  - A left-out part has such a shell, so it holds no point of the box. Its
    holes cannot add one, even a hole near the box lying outside its shell
    (possible only in an invalid input).
  - In a kept part, a point of the box is in the part when it is in the
    shell and in no hole. The left-out holes hold no point of the box, so
    the part and its subset part agree there. A hole far off that crosses
    its shell far off is such a hole.
  A valid subset's parity is the subset itself, so its part in the box is
  the right answer, and GEOS cuts it right except at the corners (below).
  The subset is cut, not the input, because GEOS reads the left-out rings
  too: shell `square(−50, 50, 300, 400)` with hole `square(−300, −300,
  −250, 500)`, outside its shell and left of the box. The subset is the
  shell, cut to 122,500 m², the right area. The input's cut is invalid,
  82,500 m², the hole subtracted. In the search below, the input's cut
  differed from the subset's (in area, or by being invalid) in 90 cases.

  *The search* (scratch `r9_search.py`, five seeds of 40,000 draws, read box
  (−200, −200) to (500, 500), coordinates on a 50 m grid from −500 to
  1,000): each draw is a shell with up to 3 holes, 30 % with a second part,
  each ring a rectangle, a random 3- to 8-gon, or a square with a spike of
  width 0, 1e-9, 1e-3 or 10 m. Of 179,456 polygons reaching the box,
  154,076 have an invalid subset and are rebuilt. Of the 25,380 with a
  valid subset (4,528 touching a corner, step 4; 79 with an invalid cut,
  step 5), none raises in the cut, and none gives an area (cut or
  intersection) whose symmetric difference with `_faces` of the whole
  input exceeds 1e-6 m² (in 45 cases `_faces` itself raises, and the cut
  was compared with `shapely.intersection(subset, box)` instead, or not at
  all for the 8 of them that touch a corner). 1,255 cuts are empty, and no
  rebuild of those keeps anything. The
  reviewer's own generator (`probe.py`, 20,000 draws of rectangles; 16,509
  reach the box, 2,823 valid subsets) gives the same: no raise, no area
  differing, no empty cut hiding an area.

  *GEOS's cut of a valid polygon with an edge through a corner of the box.*
  The search found one valid polygon cut wrong: an edge passing exactly
  through a corner of the read box, not at a vertex. The cut comes out
  inside out. The polygon `(−450, 650), (800, 900), (300, 400), (800, 650),
  (800, −500), (−450, −500)`, whose edge from (300, 400) to (800, 650)
  passes through (500, 500): cut 5,000 m², valid; its intersection with the
  box 485,000 m². In a search of valid polygons only (scratch
  `r9_corner.py`, 300,000 draws on 50, 100 and 250 m grids, 21,510 valid
  and reaching the box), 3 of the 3,136 with a corner on the boundary were
  cut wrong, and 0 of the 18,374 without. Moving that vertex 1e-9 m off makes
  the cut right, and `shapely.intersects` finds the contact exactly. Master
  overlaid whole polygons, so it never met this. On real data, no subset
  touches a corner: a corner is at the DEM points' minimum or maximum
  ±`MARGIN`, so an edge needs exact collinearity with it.

  *Two faults noted, not fixed here:*
  - `_repaired` and `make_valid` differ on a hole that encloses its shell:
    shell `square(0, 0, 100, 100)`, hole `square(−300, −300, 600, 600)`.
    `_repaired` keeps nothing (no point is in the shell and in no hole);
    `make_valid` then the cut gives 480,000 m² (the hole less the shell).
    The same on master and here: the same construction. Under the rule it
    goes to `outside` (step 3: rings simple, nothing kept).
  - The rebuild itself raises on a spike of near-zero width: `_repaired`
    of `POLYGON ((−150 −300, 0 −299.999999999, −100 0, 0 −300, 0 −500, −150
    −500, −150 −300))`, a spike 1e-9 m wide, raises GEOSException
    "CoverageUnion cannot process overlapping inputs"; at 1e-7 m it gives
    2.2e-6 m². In the search, 492 of the 154,076 rebuilds raise, 430 of
    them with a 1e-9 m spike. Master calls the same construction on every
    invalid polygon on read, so this is not new.

  *Why the many-holed polygons stay affordable.* `is_valid` pays for
  holes: GEOS tests each hole against the shell and against the other
  holes. The subset holds only the holes near the box. Numedalslågen's
  largest polygon reaching the box has 370,070 vertices and 905 holes. Its
  whole `is_valid` takes 465 ms; its subset's takes 29 ms (design review
  round 8, scratch `tim.py`). One polygon cannot be split across threads,
  so 465 ms would be the floor of a whole check on any number of threads;
  the subset check's floor is 29 ms.

  *Cost and effect on real data.* Scratch probes on battery power (compare
  only with battery runs), best of three. "Today" is `_fixed` at
  `8c5b4e2a`; "round 8" the ring check above, as a scratch version, timed
  on `8c5b4e2a`; "round 9" the rule above (scratch `r9_rule.py`), timed on
  `6ea1e369` (the same code). Each file is read as `read_class_polygons`
  reads it, with the read box from `source_region` (cell 10 m, DEM CRS
  EPSG:25833).

  | Input | Polygons read | Reaching the box | Today, 1 / 4 threads | Round 8, 1 / 4 threads | Round 9, 1 / 4 threads | Whole `is_valid` of those with a cut, 1 thread |
  |---|---|---|---|---|---|---|
  | Numedalslågen, CORINE (`corine2018`) | 7,932 | 7,239 | 319 / 133 ms | 533 / 195 ms | 748 / 331 ms | 1,084 ms |
  | Sula, clc2018 | 184 | 31 | 10 / 4 ms | 26 / 22 ms | 43 / 38 ms | 2 ms |
  | Bygdin, AR50 Innlandet + ar50_46 | 1,697 + 247 | 1,588 | 84 + 4 / 50 + 2 ms | 140 + 11 / 72 + 5 ms | 156 + 13 / 78 + 6 ms | 103 ms |

  Round 9's terms on Numedalslågen, one thread: building the subsets
  25 ms, their `is_valid` 388 ms, the corner test 4 ms (only the 15 subsets
  whose bounding box holds a corner; 10 ms over all), the cut of the
  subsets 71 ms (the inputs' cut took 73). The sum of the new terms,
  417 ms, is the measured 429 ms more than today, to noise.

  What changes in the map, against today (round 9 gives every area equal
  by WKB to round 8's on all four files):
  - Numedalslågen: nothing. All 7,237 areas are equal by WKB, 695
    outside, 0 rebuilt, no subset invalid, none touching a corner. So RP8's
    re-pin and the quick check's meshes hold.
  - Sula: 27 of 28 areas are equal by WKB. The 28th has a self-touching
    ring near the box: it is now rebuilt, with the same area. One polygon
    that was `outside` has a non-simple ring near the box and nothing in
    it: it is now dropped. `repaired` goes from 0 to 2.
  - Bygdin, Innlandet: 1,449 of 1,453 areas are equal. Two polygons with
    a self-crossing ring near the box were kept as their cut today, and are
    now rebuilt: the same area in the box, to 0 m². Two valid polygons
    whose cut GEOS left invalid were kept whole today (3.08 and 0.43 km²).
    Now each is its intersection with the box (1.50 and 0.26 km²), the same
    there to 0 m². `repaired` goes from 5 to 7 (round 8 counted those two
    as well: 9). ar50_46: unchanged.
  - Every kept cut equals `shapely.intersection(subset, box)` to 0 m² on
    all four files. On these inputs, no cut today differs in area from the
    rebuild: the flaw is real but did not fire here.
  - With the rule patched into a scratch pytest plugin (`r9_plugin.py`),
    all 668 tests in `test_feature_input.py`, `test_land_cover_map.py`,
    `test_read_repair.py`, `test_cli_mesh_edges.py` and `test_preferred.py`
    pass, test 7 among them (Sula, 3,271 triangles), and today's 30
    failing ones (`TestACutToAThreePointRing`). Planted check: the plugin
    dropping every area fails 123, test 7 among them (2,466 triangles).

  Rejected: `is_valid` on each whole polygon that reaches the box. It costs
  1,084 ms on Numedalslågen, with a 465 ms floor on any number of threads,
  against the subset check's 429 ms more on one thread and 198 ms more on
  four, with a 29 ms floor. The quick check runs Numedalslågen with every
  core (`threads = [0]`).
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
  `src_python/tin_engine/cli.py@75ab2b4f:2017`, `src_python/tin_engine/cli.py@75ab2b4f:2024`); `--prefer` with a line key (the message names
  `--lines`); `--lines` with an area key (naming `--prefer`); a polygon
  source under `--lines`; a line source under `--prefer`; a second source
  for one `--lines` key.
- **Layer and CRS of a file.** A catalogue key carries its entry's layer;
  the `--land-cover` file has `--land-cover-layer` and `--land-cover-crs`
  (the quick check's Numedalslågen and Lagan cases use them), passed to
  `read_source`'s `layer` and `crs` as today (`src_python/tin_engine/feature_input.py@75ab2b4f:656`). A
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
  them today (`src_python/tin_engine/cli.py@75ab2b4f:2129`, codes by `label_triangles`). `--prefer` is
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
12. What every polygon in the map is (section 3's rule; unit tests,
    `test_feature_input.py`'s `BOX`, read region `BOX` grown by 100 m, read
    box that grown by 100 m more). Measured on tree `435d407b`:
    - `build` on one file of `INNER` (311, 10,000 m², inside the box) and
      a bowtie 50 km off (`at(50000, 50000)`, `at(50100, 50100)`,
      `at(50100, 50000)`, `at(50000, 50100)`, invalid, cut empty): `areas`
      is `INNER` alone, 10,000.0 m², and `outside` is the bowtie, the same
      object as given; `repaired == 0`; `shapely.is_valid` holds on every
      polygon of `areas`. Today `areas` holds both, the bowtie invalid
      (0.0 m² by GEOS's reading).
    - p5's two polygons with `CLEANUP` (`one(..., switches=CLEANUP)`): no
      feature, `outside == 1`, no exception, for both. Today the first
      gives that, and the second (shell `ELL`) raises GEOSException in
      `_Tally.clean`.
    - p5 without clean-up: its existing assertions unchanged, and the one
      feature's `polygon` is None.
    - The round-4 red test (`test_cli_mesh_edges.py`, the far bowtie under
      the long lake) and round 3's (the overlapping out-of-box pair), as
      they stand.
    - `shapely.is_valid` also holds on every polygon of `areas` in both
      runs, and on the p1 bowtie after repair.
13. Deciding on the input polygon (section 3, design rounds 8 and 9; unit
    tests on `build`, test 12's `BOX`, read box (−200, −200) to (500, 500)
    from `at`'s origin; areas to 1e-6 m²). Measured on tree `8c5b4e2a`
    (round 9's cases on `6ea1e369`, the same code, with the rule in scratch
    `r9_rule.py`) with shapely 2.2.0, GEOS 3.14.1.
    - `TestACutToAThreePointRing` (`8c5b4e2a`) stands as written. Each
      case has a ring near the box that is not simple, so its subset is
      invalid, it is rebuilt, and `clip_by_rect` never sees it:

      | Case | Rule gives |
      |---|---|
      | the 8 of `NO_AREA_IN_THE_BOX` (one spike into, across, or through the corner of the box; two spikes; a tip one ulp off; a self-crossing square with a spike; a flat ring across the box's edge; two points and back) | rebuilt to nothing, dropped: `areas` holds the neighbour `SMALL` alone, `outside` is empty, `repaired == 1` |
      | `HOLE_SPIKE` (the spike as a hole of a 2 km square) | rebuilt: one area, the whole read box, 490,000 m² |
      | `PART_SPIKE` (`SMALL` with the spike as a second part) | rebuilt: one area, `SMALL`, 1,600 m² |

      Each test's premise (the cut raises) stays: it is what makes the
      order matter.
    - The inside-out case: `ring(*EAST, (600, 50), (150, 50), (600, 50 +
      w))`, invalid, for `w` = 10 and 1e-3. Premise: its cut is valid and
      covers 488,638.9 m² at `w = 10`. `build` gives one area of 136.1111·`w`
      m² (1,361.1111 and 0.1361 m²), and its symmetric difference with
      `_repaired`'s result is under 1e-6 m². `repaired == 1`, `outside`
      empty. Red today: one area of 488,638.9 m².
    - Its valid twin, the spike drawn the other way (`ring(*EAST, (600, 50 +
      w), (150, 50), (600, 50))`, valid): kept as its cut, the same
      136.1111·`w` m², `repaired == 0`. So the input decides, not the area.
    - A spy on `shapely.clip_by_rect` sees no polygonal input from any
      invalid case above (lines from `_faces` are allowed); the valid twin
      is cut, as its subset. Planted check: round 7's
      `_fixed` makes the spy fire on every case.
    - Faults between simple rings: a hole across its shell in the box
      (shell `square(0, 0, 300, 300)`, hole `square(250, 100, 350, 200)`).
      Its subset is invalid, so it is rebuilt: 85,000 m², `repaired == 1`
      (its cut is invalid, 80,000 m²). The same fault beyond the box (shell
      `square(−1000, −1000, 1000, 1000)`, hole `square(900, 100, 1100,
      200)`): the hole is left out of the subset, which is valid; kept as
      its cut, the whole read box, 490,000 m², `repaired == 0`.
    - Design review round 8's cases, each with simple rings and an invalid
      subset, so rebuilt; the spy sees none of them:
      - A hole poking out through its shell across the box's edge: shell
        `square(−50, −100, 200, 900)`, hole `square(−100, −200, 950, 550)`.
        Premise: its cut is valid, 220,000 m² (`make_valid` then the cut
        270,000). The rebuild keeps nothing and its rings are simple, so
        `areas` is empty, `outside` is the polygon as given, `repaired ==
        0`. Red today: one area of 220,000 m².
      - An empty cut whose rebuild is not empty: shell `square(−250, −500,
        950, 550)`, holes `square(−500, −200, 850, 600)`, `square(−100,
        −300, 850, 450)` and `square(500, 200, 700, 700)`. Premise: its cut
        is empty. `build` gives one area of 390,000 m², `repaired == 1`,
        `outside` empty. Red today: no area, `outside` holds it.
      - A cut that raises on simple rings: shell `square(−300, −300, 900,
        100)`, hole `square(−350, −300, 850, 100)` (sharing two sides).
        Premise: the cut raises GEOSException (a ring of 3 points). The
        rebuild keeps nothing: `outside` holds it, `repaired == 0`. Red
        today: GEOSException.
    - The subset is cut, not the input: shell `square(−50, 50, 300, 400)`,
      hole `square(−300, −300, −250, 500)` (outside its shell, left of the
      box). Premise: the input's cut is invalid, 82,500 m². `build` gives
      one area of 122,500 m², `repaired == 0`, and the spy sees the shell
      alone. Red today: the same area, but `repaired == 1`.
    - A valid polygon with an edge through a corner of the box:
      `ring((−450, 650), (800, 900), (300, 400), (800, 650), (800, −500),
      (−450, −500))`, valid, its edge from (300, 400) to (800, 650) through
      the corner (500, 500). Premise: its cut is valid and covers 5,000 m².
      `build` gives one area of 485,000 m², whose symmetric difference with
      `shapely.intersection` of the polygon and the read box is under 1e-6
      m², `repaired == 0`; the spy sees no polygonal input. Red today: one
      area of 5,000 m².
    - p5's two polygons and the far bowtie: as test 12, `outside`,
      `repaired == 0`. The far bowtie misses the box; p5's first shell is
      left out of its subset, which is empty and valid; p5's L-shaped shell
      makes an invalid subset (its hole lies outside it) whose rebuild keeps
      nothing, with simple rings.
    - The record's `land_cover_repaired` wording is section 3's.

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

*Timing, the polygons' check.* Three new terms, design round 9's subset
check (section 3, *Deciding on the input polygon*): building the subsets
(25 ms), their `is_valid` (388 ms) and the corner test (4 ms), on
Numedalslågen, one thread, battery. Measured whole: 429 ms more than today
on one thread and 198 ms more on four, against a total of about 9 s, so
under 5 % on one thread and under 2.5 % on four, inside the quick check's
5 % band (`tools/bench_quick.py`, the band is at least 5 %). The quick
check runs on every core. Sula and Bygdin add under 90 ms. The outside
polygons are only set aside. The rejected whole checks are timed in
section 3.

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

## As built

- RP8's re-pin to the cut kept at a4bcab1e (Ola 2026-10-29), checked on the
  quick check's Numedalslågen case against 2f70ce84 (whole polygons):
  `rasputin mesh --dem $RASPUTIN_DATA/DTM10_UTM33_20260925 --domain
  $RASPUTIN_DATA/numedalslagen_outline_nve.geojson --land-cover
  $RASPUTIN_DATA/corine2018_dtm10_utm33.gpkg --land-cover-layer corine2018
  --land-cover-classes corine --edges all --tolerance 10`, with an empty
  `XDG_CONFIG_HOME` (Ola's sources.toml still says `rules =`). Both meshes:
  653,039 triangles, 334,190 vertices; 28 vertices moved, 22 by under
  0.015 m, the farthest 8.745 m. Cause: the clean-up's repair step, not the
  lines. With every clean-up switch at 0 the meshes are identical
  (1,140,546 triangles), as with only simplification or only the outline
  snap on. The cut changes the vertex order (part order or ring start) of
  22 of the 3,297 polygons, and 2 more by 1.5e-11 m; `shapely.coverage_clean`
  (1 m snap) depends on that order, so 23 cleaned polygons differ, their
  borders by at most 0.86 m.

## Review

- Design round 1: CHANGES REQUESTED, 75ab2b4f..a2f4dcaa, 0 (no speed estimate for the map and per-class edges; line files lose layer and CRS; Ola's sources file breaks on the rules → classes rename; six wordings and inventory gaps).
- Design round 2: APPROVED WITH WORDINGS, af82b7d3..fd566170, 0; question 3 reworded and question 4 added in this commit.
- Design round 3: CHANGES REQUESTED, fd566170..0da191d5, 0 (question 4's premise false: Ola's 27 Oct Sula 2 m and 5 m runs fed the 20 m mesh's `lines` layer back, which `--lines` refuses; question 5 narrowed to the `--land-cover` half, worded from the 26 Oct Swedish runs; `--lines-crs` and `--prefer-crs` have no user; three wordings).
- Design round 4: CHANGES REQUESTED, 132a89cf..498e940a, 0 (several land-cover files in one run dropped unsaid: bygdin-ar50-vtol5 and oyangen-ringkollen pass two AR50 county files; Sula two-stage probe's 25.7 m exceeds the stated 20 m band and the 2 % beyond 5 m is not located; no cut-run (`--pieces`) test, `--prefer` in cut runs unstated).
- Design round 5: CHANGES REQUESTED, 369bac65..bf9b88de, 0 (Sula and cut-run findings fixed; §3's overlap query `STRtree.query(predicate="overlaps")` returns nothing for a polygon wholly inside another, so test 10's planted box `P` inside forest 311 gives no overlap: query with "intersects" and keep pairs of shared area above zero, and re-measure Bygdin's 0 m² if it used "overlaps"; test 10 bullet 2 ignores that `--prefer lake=` replaces `P` (class 512) as 311, measured 311 1,409,620 m², 512 420,089.9 m², so the stated 1,369,460 m², 460,250 m² and 3,402.0 m water lines are wrong; wording: test 11's "before any phase" vs the refusal after DEM open and piece split, `cli.py:1186-1203`).
- Design round 6: APPROVED, 12616b91..3ea5ceb8, 0 (all three round-5 fixes hold: Bygdin re-probed, 345 intersecting pairs, none sharing area, 0.24 s; test 10 figures checked by hand, P 40,160.07 m², 311 1,369,459.9 m², lines 9,884.3 m; test 11 refusal after `_open_dem` and `partition`, before `assemble_dem` and `_open_features`).
- Code review round 1: CHANGES REQUESTED, 75ab2b4f..dc9969ea, +14 (1,052 added, 1,038 removed; all seven gates pass; full suite 11,229 passed, 18 skipped, 1 xfail; Sula test 7 passed; Bygdin two-file AR50 run 199,510 triangles, overlap none; the -300 estimate is what was wrong (round 2, before several files, the overlap record and section 8), nothing beyond the design; the six test amendments each supported by the design. Blocking: the inventory's docs rewording never done: README.md:101 and INSTALL.md:269 give refused `--features ... --features-rules corine` commands, "meshing rules" and `--features-*` remain in README, INSTALL, docs/GLOSSARY.md:256, project_structure.md:220,306; a self-crossing polygon in an overlapping pair of two `--land-cover` files crashes with GEOSException in `land_cover_map._joined` (bowtie reproduced): repair the sharing pair or catch and repair, with a test. Follow-ups: `run_record.WORDING["land_cover_accuracy_m"]` never written; `LandCoverClasses.accuracy_m` optional though the loader requires it. The main session adds: the `rules =` refusal names the file and key but not the fix; it should say to rename it to `classes`).
- Code review round 2: CHANGES REQUESTED, dbfa0bdd..57c9a0a2, +19 (+5 this round; round-1 findings closed: docs match the code, no user-facing `--features*` or "meshing rules" left, `rules` refusal names the rename; green commit touches no test file; seven gates pass; full suite 11,232 passed, 18 skipped, 1 xfail; merges cleanly with 9c317c16; both polygons of each pair repaired, three-file probe passed. Blocking, a regression against master: a polygon in no overlapping pair (a single file's, or the earlier file's) stays unrepaired, and `preferred.replace` crashes on it with GEOSException at `preferred.py:70` `_cut` (one `--land-cover` file with a self-crossing polygon plus `--prefer lake=`; master passes, its reader repaired first): repair every polygon near the region once in `build`/`_joined`, keep the repair count in `repaired_areas` and `land_cover_repaired`, test single file plus `--prefer`, fix the unit docstring "an earlier file's bowtie is not checked"; the same fix removes a crash master has, `--edges` with clean-up on such a polygon (`src_python/tin_engine/feature_input.py@435d407b:362`). Follow-ups from round 1 still open).
- Code review round 3: CHANGES REQUESTED, 57c9a0a2..b605b1c8, +43 (1,081 added, 1,038 removed; round-2 findings closed; touched suites 298 passed, mypy and ruff clean; blocking: `_fixed` keeps a polygon outside the region's box whole and unchecked, so two overlapping out-of-box polygons, one a bowtie, raise GEOSException through the CLI with two GeoJSON `--land-cover` files; `THREADS = 8` ignores `--threads`, making README.md:87-88 false; RP8's re-pin evidence (master comparison, Numedalslågen mesh at a4bcab1e, 28 vertices moved, at most 8.7 m, and why) is not recorded in the branch).
- Code review round 4: CHANGES REQUESTED, b605b1c8..31fb6608, +52 (+9 this round; round-3 findings closed: `--threads` reaches the build, the overlap join skips out-of-box polygons, RP8 evidence recorded; new tests fail on b605b1c8; 170 tests pass; blocking: an out-of-box polygon is still kept whole and unchecked, so a bowtie 50 km off plus a `--prefer lake=` lake reaching it crashes `preferred.py:70` `_cut` with GEOSException, where master passes).
- Design round 7: APPROVED WITH WORDINGS, 435d407b..03b0247d, 0 (the map's polygon rule: `areas` valid and checked once in `_fixed`, `outside` whole and read only by rings; all 10 consumers confirmed plus `cli.py:1342`'s `cover`; p5 and the round-3/4 counts hold; the rejection of whole-polygon repair stands on the timings and on `make_valid` breaking p5; three wordings applied in this commit).
- Code review round 5: CHANGES REQUESTED, dfcf2e94..2b678e9f, +59 (+7 this round; all 10 consumers keep round 7's rule; 319 touched tests pass; RP8 re-pin re-checked, 60 rows equal by WKB, 59 fids moved; blocking: `_fixed`'s `clip_by_rect` raises GEOSException on a self-crossing polygon with a zero-width spike into the read box (a 3-point ring after the cut), for every `--land-cover` run, since a4bcab1e; master gives one feature; fix: catch per polygon and send it to `_repaired`, plus a test).
- Design round 8: CHANGES REQUESTED, 8c5b4e2a..8a523604, 0 (validate before the cut: two simple rings still break the cut-then-check step; a hole poking out through its shell across the box edge gives a valid cut of 220,000 m² against 0 rebuilt and 270,000 by make_valid, and a random search found empty cuts whose rebuild is not empty; the cut also raises on simple rings (50 cases found); fix: `is_valid` before the cut on the parts whose shell is near the box, each with only its holes near the box, in place of `is_simple`, 29 ms on Numedalslågen's largest polygon and 380 + 177 ms over its 7,239 polygons; the four pins hold except the uncaught raise).
