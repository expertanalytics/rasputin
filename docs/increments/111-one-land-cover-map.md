# Increment 111: one land-cover map, and the edges chosen from it

**Status: built** (red `bce16583`, green `a8246816`, code review round 1's fixes `6cd8bea3` and `77b3def6`; at `77b3def6`, 19 net production lines, 1,057 added and 1,038 removed, by `tools/count_loc.py 75ab2b4f 77b3def6`; full suite 11,232 passed, 18 skipped, 1 xfailed; README, INSTALL, GLOSSARY and `project_structure.md` reworded; since then code review rounds 2 to 5 and design rounds 7 to 11; the branch stands at +75 by `tools/count_loc.py 75ab2b4f 486b0e40`). Design round 10 (section 3, *What every polygon in the map is*; tests 12 and 13) replaces rounds 7 to 9's polygon rule: each polygon reaching the read box is repaired whole, as master repairs on read, then cut with `shapely.intersection`. Design round 11 closes round 10's
review: a rebuild that raises is refused inside `_repaired` on all three
paths (section 3, step 6). The shares' area loss in `subdivide` stays out
(Ola 2026-10-29, "split"; pipe item `fractions-intersection`). Under h39 (the harness change by which a design fix after a review goes to code review, not to another design review; approved by Ola 2026-10-29, merged in PR #402) it has no design review round: code review checks it.
Designed, approved in design review round 6, section 8 (the Python interface and the wordings tests read) added for the test migration, by `@architect`, 2026-10-28, on master `75ab2b4f`,
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
and fetches none. Design round 10 overlaps "Repair CORINE once per area of
interest": a cached, repaired extract would need no check per run; Ola's
private cleaning plugin is where it goes. Not taken up.

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
| `FeatureSource`, `open_features`, `_Tally.source`, `_take`, `_bound`, `_linework`, `R_MIN`, `R_MAX` | **deleted**: reading is `read_class_polygons` (reading as `_take`'s cleaned-polygon branch does: one query, one move; its `_repaired` moves to `land_cover_map.build`, section 3), `read_polygons`, new `read_linestrings`. `pre_clip` and `_runs` stay (`feature_input.py`): `edges_of`'s path without clean-up draws an outside polygon's rings with them (section 3, *Outside polygons*) |
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
after the mesh: shares → assign_polygon_fractions(map.areas) → shares, codes, cells
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
`--land-cover` do. The one step the edges and the shares have in common is
reading the polygons. Both get them repaired once, in `build` (*What every
polygon in the map is*, below), and the check skips the polygons whose
bounding box misses the read box, where master spends most of its check
(section 4, read probe).

- **The map** (`tin_engine/land_cover_map.py`, frozen dataclass
  `LandCoverMap(areas, outside, classes, values)`; None without
  `--land-cover`). The base is the class polygons over the read region
  (`source_region`): one `read_source` query, every polygon moved in one
  `shapely.transform` call, then repaired and cut in `build` (next bullet);
  or no areas for a class raster, read after the mesh as today. Each
  `--prefer` is 107's `replace` near the read region, first given wins, and
  each refilled piece is unioned with the same-class polygons it touches
  (110's join, moved from `_Tally._preferred`, on `(code, polygon)` pairs),
  on the map's valid `areas`. The record's replaced, given and refilled
  areas are measured inside the domain.
- **What every polygon in the map is** (design round 10). Rounds 7 to 9
  cut first and checked less than the whole polygon, to save time; each
  check left a polygon the cut got wrong. Ola ruled on 2026-10-29 that
  repair leaves the open core for a private cleaning plugin, and that the
  core will then refuse invalid input. Until then the map repairs as master
  does, then cuts with the overlay that is right on any valid polygon. The
  *read box* is the read region's bounding box widened by `MARGIN` (100 m,
  the margin `_repaired` and the GeoPackage query already use; checked on
  catchments up to Lagan's, 7,866 polygons read). In `_fixed`, in `build`
  with `fix` on, in the thread pool's groups as now, per polygon `g`:
  1. *Misses the box.* `g`'s bounding box misses the read box (closed
     comparison): to `outside`, as read, unchecked. No ring of it can enter
     the box or enclose it.
  2. *Repaired whole.* Otherwise `_repaired(g, xy)`, master's repair on
     read: `g` itself when `shapely.is_valid(g)`, else its rebuild in the
     read box (its rings cut to the box, faces kept by parity). Each invalid
     `g` is counted in `repaired`, whatever its rebuild keeps (master's
     count, and master's record wording, "Land-cover polygons that crossed
     themselves, repaired").
  3. *Cut.* A valid `g` whose bounding box lies inside the read box is kept
     as read: its cut would be itself, and keeping it keeps master's map
     there (6,893 of Numedalslågen's 7,237 areas). A valid `g` that crosses
     the box's edge is cut with `shapely.intersection(g, read box)`. A
     rebuild already lies in the box and is not cut. Never
     `shapely.clip_by_rect`: GEOS's rectangle clip is right only on a valid
     polygon with no edge through a corner of the box (design review round
     9: 3,920 of 5,505 valid polygons with an edge exactly through a corner
     cut wrong, 5,697 m² against 484,303 m²).
  4. *Its polygonal part.* The result itself when it is a Polygon or
     MultiPolygon of area above zero. From a GeometryCollection (a cut that
     touches the box along a line or at a point beside its area), its
     Polygon parts, with each MultiPolygon's parts, as one MultiPolygon, or
     the Polygon when there is one. GEOS's overlay gives valid polygons, and
     its polygonal parts do not overlap, so the part is valid (test 12).
  5. *Kept or set aside.* A part with area goes to `areas`. Otherwise (a
     rebuild that keeps nothing, or a cut that is empty, only lines or only
     points) `g` goes to `outside`, as read.
  6. *A rebuild that raises* (design round 11). GEOS raises
     GEOSException in the rebuild on some spikes up to about 1e-7 m wide
     (below). `_repaired` itself catches it and raises `FeatureError`:
     "<what> near (x, y) in the DEM's CRS cannot be repaired: it has a spike
     too thin to rebuild; clean the file", (x, y) the polygon's first vertex
     to 0.1 m, `<what>` its new keyword `what`. So every caller refuses
     it, and no traceback is left:
     - `_fixed` passes "a land-cover polygon of class N"; the CLI wraps
       `build` (`cli.py`'s `base = build(...)`, today outside any `try`) in
       the `except FeatureError` of the reads above it, hint `--land-cover`;
     - `read_polygons` passes "<file>: feature <fid>" (it keeps each row's
       fid beside its polygon for this). Its two CLI callers already turn a
       `FeatureError` into exit 2: `--prefer` (hint `--prefer`, `cli.py`'s
       `_land_cover_map`) and `--reservoirs` (hint `--reservoirs`, after
       the mesh). Round 10 found both ending in a traceback.
     It neither drops the polygon nor guesses its area, since inflow between
     river points needs the classes where they are, and it is the refusal
     the core will give every invalid polygon once repair moves to the
     plugin.

  `build` with `fix` off (the `--prefer` build in the CLI) takes its base as
  `areas` already so and is given the base map's `outside` (a keyword,
  default empty), which it passes through. The preferred polygons are valid
  as `read_polygons` gives them (`_repaired` each). GEOS's overlays of valid
  polygons give valid polygons (the OverlayNG documentation: results are
  always valid; not probed here), so every polygon `_joined`, `replace` and
  the refill join make is valid too. So any step after `build` may overlay
  `areas` without a check, and no step overlays `outside`. A test asserts
  `shapely.is_valid` on every polygon of `areas` (section 5, test 12).
  `feature_input.repaired_areas`, which no production code calls, is
  deleted; the shares step takes `map.areas` as they are.

  *Outside polygons, counted and kept.* Without clean-up (the pre-clip
  path of `edges_of`), each chosen class's outside polygon goes through
  today's ring test (`region.intersects(g.boundary)`, then `pre_clip`), so a
  hole lying in the region whose shell misses it (`test_feature_input.py`'s
  p5 cases) still gives its ring, exactly as given. Its feature carries no
  polygon: it labels no triangle, since it has no area in the box. One with
  no ring kept is counted outside. With clean-up on, or with the map cut by
  `--prefer`, each chosen class's outside polygon is counted outside, with
  no overlay: the clean-up works on area inside the read region, and an
  outside polygon has none there. The outside polygons come after the
  `areas` in the edges step's feature order. The shares and the codes never
  see `outside`.

  So an invalid polygon with no area in the box (p5's L-shaped shell, the
  spikes of `NO_AREA_IN_THE_BOX`) is drawn as master's `--features` route
  drew it, which never repaired, and left out of the shares as master's
  `--land-cover` route left it out, which dropped it on read. Probed on
  master `43f1b2e9` with master's `open_features` under `corine`, beside
  `SMALL` (class 312): without clean-up each spike's ring is a line of
  class 311 (150 m for the spike into the box; 300 m for the spike across
  it; 70.711 m through the corner; two lines, 300 m, for two spikes) and
  nothing is counted outside; with clean-up only `SMALL` is a feature and 1
  is counted outside. One difference: design review round 8's hole poking
  out through its shell, with clean-up, raises GEOSException on master
  ("side location conflict"); here it is counted outside.

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
     region)` (`src_python/tin_engine/feature_input.py@435d407b:362`):
     today p5's L-shaped polygon with `CLEANUP` raises GEOSException
     ("TopologyException: unable to assign free hole to a shell").
     Probed on tree `435d407b`. It now gets `areas` only.
  9. The shares, `assign_polygon_fractions(map.areas)` in `cli.py`: it gets
     `areas` only.
  10. `label_triangles` over the pre-clip path's polygons: an outside
      feature has none.

  Out of scope (Ola 2026-10-29, "split", section 7): `subdivide`'s
  `clip_by_rect` can still lose area on valid polygons, so the shares
  (step 9) and `areas_inside` (step 5) can come out short (design review
  round 10: 10 of 2,817 synthetic polygons, up to 57,515 m²). Pipe item
  `fractions-intersection` fixes it, after `core-assumes-clean`.

  *Known faults of the repair, left to the cleaning plugin:*
  - `_repaired` and `make_valid` differ on a hole that encloses its shell:
    shell `square(0, 0, 100, 100)`, hole `square(−300, −300, 600, 600)`.
    `_repaired` keeps nothing (no point is in the shell and in no hole);
    `make_valid` then the cut gives 480,000 m² (the hole less the shell).
    The same on master: the same construction. It goes to `outside`.
  - The rebuild raises on a thin spike: GEOSException "CoverageUnion cannot
    process overlapping inputs". Scratch `r10_spike.py`: round 9's spike
    (from (0, −300) to a random tip in the box and back, its sides `w` apart
    at the base), 400 tips per width, at the origin and at
    (500,000.3, 6,600,000.7). Of the invalid ones, at the origin and at
    those coordinates: `w` = 1e-9 m, 100 of 103 and 56 of 103 raise;
    1e-8 m, 76 of 101 and 62 of 95; 1e-7 m, 2 of 104 and 1 of 108; 1e-6 m,
    none. Master raises the same in `read_class_polygons` and
    `read_polygons`, so a master `--land-cover`, `--prefer` or
    `--reservoirs` run on such a polygon ends in a traceback; master's
    `--features` runs never repaired, so they did not meet it. Step 6 turns
    each traceback into a refusal. None of the four real cases holds one:
    their 11 invalid polygons reaching the box (Sula 4, Bygdin 7) all
    rebuild.

  *Cost and effect on real data* (scratch `r10_time.py` and
  `r10_terms.py`, tree `492f5d8c`, battery power: compare only with battery
  runs; best of three). Each file is read as `read_class_polygons` reads
  it, with the read box from `source_region` (cell 10 m in EPSG:25833;
  Lagan 30 m in EPSG:3006). "Master" is master's repair on read,
  `_repaired` of every polygon read, one thread (master runs it so, and
  only in `--land-cover` runs; its `--features` runs check nothing).
  "Today" is the branch's `_fixed` at `492f5d8c`.

  | Input | Read | Bounding box meets the read box | Master, 1 thread | Today, 1 / 4 / 10 threads | Round 10, 1 / 4 / 10 threads |
  |---|---|---|---|---|---|
  | Numedalslågen, CORINE (`corine2018`) | 7,932 | 7,239 | 3,478 ms | 323 / 133 / 107 ms | 1,359 / 655 / 539 ms |
  | Lagan, CORINE | 7,866 | 7,150 | 11,783 ms | 336 / 222 / 201 ms | 1,896 / 572 / 508 ms |
  | Sula, clc2018 | 184 | 31 | 1,391 ms | 10 / 4 ms | 76 / 66 ms |
  | Bygdin, AR50 Innlandet + ar50_46 | 1,697 + 247 | 1,454 + 134 | 236 + 106 ms | 87 + 4 / 50 + 2 ms | 188 + 27 / 62 + 11 ms |

  Round 10's terms on Numedalslågen, one thread: the whole `is_valid` of
  the 7,239 polygons 1,026 ms, the intersection of the 346 that cross the
  box's edge 313 ms; their sum less today's cut and check (309 ms) is
  1,030 ms, the measured 1,036 ms more than today, to noise. On every core
  the floor is one polygon's `is_valid`: the largest, 465 ms, which no
  thread count splits. Master spends most of its 3,478 ms on the 693
  polygons whose bounding box misses the read box, which round 10 does not
  check.

  What changes in the map, against today: on all four files the same
  number of areas, the same codes in the same order, and each area equal to
  today's inside the read box to 0 m² (symmetric difference); `outside`
  holds the same number; every area is valid. `repaired` goes from 0 to 4
  on Sula and from 5 to 7 on Bygdin's Innlandet file (each invalid polygon
  reaching the box now counts), and stays 0 elsewhere. Today's cut gives
  each polygon new WKB even inside the box; round 10 keeps 6,893 of
  Numedalslågen's 7,237 areas as read (Lagan 6,916 of 7,147, Bygdin 1,450
  of 1,587, Sula 9 of 28). With round 10 patched into a scratch pytest
  plugin (`r10_plugin.py`), 657 of the 681 tests in
  `test_feature_input.py`, `test_land_cover_map.py`, `test_read_repair.py`,
  `test_cli_mesh_edges.py`, `test_preferred.py` and `test_cli_mesh_prefer.py`
  pass, test 7 among them (Sula, 3,271 triangles); the 24 that fail are the
  `TestACutToAThreePointRing` cases test 13 changes. Planted check: the
  plugin dropping every area fails 125, test 7 among them.
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
  `--edges` selector names a class or `--shares` (flag) asks.
  `assign_polygon_fractions` gives each triangle's shares from the map's
  `areas`, repaired in `build`; or the class raster's route, as today. The
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
order of size only. Most of master's check is spent on polygons whose
bounding box misses the read box (Numedalslågen: 693 of 7,932, over 2.4 s
of 3.5 s, design round 10), so `build` checks only the polygons reaching
it (section 3).

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
    Stands as written under round 10 (`r10_plugin.py`: every test above
    passes). Added: `INNER` comes back as the very object given (`is`),
    being inside the read box; and a valid U-shaped polygon crossing the
    box's edge, `ring((450, 0), (700, 0), (700, 300), (500, 300), (500,
    200), (600, 200), (600, 100), (450, 100))`, whose intersection with the
    read box is a GeometryCollection of a 5,000 m² Polygon and a 100 m
    LineString along x = 500 (premise), gives one area, that Polygon,
    5,000 m², valid, `repaired == 0`.
13. Repair whole, then cut (section 3, design round 10; unit tests on
    `build`, test 12's `BOX`, read box (−200, −200) to (500, 500) from
    `at`'s origin; areas to 1e-6 m²). Measured on tree `492f5d8c` with the
    rule in scratch `r10_time.py` (`r10_cases.py`), shapely 2.2.0, GEOS
    3.14.1. Design rounds 8 and 9's subset, corner and `is_simple` cases are
    deleted with the rule they tested, except where kept below.
    - `TestACutToAThreePointRing` (code review round 5) changes. Its
      premises stand (each cut raises; each of the 8 `NO_AREA_IN_THE_BOX`
      rebuilds to nothing), and `HOLE_SPIKE` (one area, the whole read box,
      490,000 m²) and `PART_SPIKE` (one area, `SMALL`, 1,600 m²) stand,
      `repaired == 1`. The 8 `NO_AREA` cases now go to `outside`, as given,
      `repaired == 1`, `areas` `SMALL` alone. Their runs (`open_one`,
      beside `SMALL` of class 312) give master's `--features` figures
      (section 3, probed on `43f1b2e9`; the same under `r10_plugin.py`):
      without clean-up, features 312 then 311 (outside polygons come after
      the areas), 311's lines totalling 150 m for the spike into the box, the
      one-ulp tip, the self-crossing square, the flat ring and the two
      points and back; 300 m for the spike across the box; 70.711 m
      through the corner; two lines, 300 m, for two spikes; `outside == 0`.
      With `CLEANUP`: feature 312 alone, `outside == 1`.
    - The inside-out case: `ring(*EAST, (600, 50), (150, 50), (600, 50 +
      w))`, invalid, for `w` = 10 and 1e-3. Premise: its `clip_by_rect` is
      valid and covers 488,638.9 m² at `w = 10`. `build` gives one area of
      136.1111·`w` m² (1,361.1111 and 0.1361 m²), whose symmetric difference
      with `_repaired`'s result is under 1e-6 m²; `repaired == 1`,
      `outside` empty. Its valid twin (`ring(*EAST, (600, 50 + w), (150,
      50), (600, 50))`): the same area, `repaired == 0`.
    - A spy on `shapely.clip_by_rect` sees no polygonal input from any case
      here (`_repaired`'s cut of the rings is lines). Planted check:
      today's `_fixed` makes it fire on every case.
    - Faults between simple rings, all invalid, so rebuilt and counted
      (`repaired == 1`): a hole across its shell in the box (shell
      `square(0, 0, 300, 300)`, hole `square(250, 100, 350, 200)`):
      85,000 m² (its `clip_by_rect` is invalid, 80,000 m²); the same fault
      beyond the box (shell `square(−1000, −1000, 1000, 1000)`, hole
      `square(900, 100, 1100, 200)`): the whole read box, 490,000 m²; a
      hole outside its shell, left of the box (shell `square(−50, 50, 300,
      400)`, hole `square(−300, −300, −250, 500)`): 122,500 m² (its
      `clip_by_rect` is invalid, 82,500 m²); an empty `clip_by_rect` whose
      rebuild is not (shell `square(−250, −500, 950, 550)`, holes
      `square(−500, −200, 850, 600)`, `square(−100, −300, 850, 450)`,
      `square(500, 200, 700, 700)`): 390,000 m².
    - Rebuilt to nothing, so `outside` as given, `repaired == 1`: a hole
      poking out through its shell across the box's edge (shell
      `square(−50, −100, 200, 900)`, hole `square(−100, −200, 950, 550)`;
      its `clip_by_rect` is valid, 220,000 m²); a hole sharing two sides
      with its shell (shell `square(−300, −300, 900, 100)`, hole
      `square(−350, −300, 850, 100)`; its `clip_by_rect` raises); a hole
      enclosing its shell (shell `square(0, 0, 100, 100)`, hole
      `square(−300, −300, 600, 600)`).
    - A valid polygon with an edge through a corner of the box:
      `ring((−450, 650), (800, 900), (300, 400), (800, 650), (800, −500),
      (−450, −500))`, its edge from (300, 400) to (800, 650) through the
      corner (500, 500). Premise: its `clip_by_rect` is valid and covers
      5,000 m². `build` gives one area of 485,000 m², `repaired == 0`. Red
      today: 5,000 m².
    - A valid polygon touching the box along an edge only (`ring((500, 0),
      (600, 0), (600, 100), (500, 100))`): `outside`, as given,
      `repaired == 0`.
    - A rebuild that raises: `Polygon([at(−150, −300), at(0,
      −299.99999999), at(−8.7, 41.9), at(0, −300), at(0, −500), at(−150,
      −500)])`, a spike 1e-8 m wide at its base. Premise: invalid.
      `_repaired(g, xy, what="x")` raises `FeatureError` whose message
      starts "x near (", holds "cannot be repaired" and the first vertex to
      0.1 m (red at `c5cd80c6`: GEOSException, "CoverageUnion cannot
      process overlapping inputs"). `build` raises `FeatureError` naming
      class 311 and "cannot be repaired"; the same polygon in a
      `--land-cover` GeoJSON through `rasputin mesh ... --edges all` exits 2
      with that message and the hint `--land-cover`, and no traceback.
    - p5's two polygons and the far bowtie: as test 12, `outside`. The far
      bowtie and p5's first polygon miss the box (`repaired == 0`); p5's
      L-shaped shell reaches it, is invalid (its hole lies outside it) and
      rebuilds to nothing (`repaired == 1`).
    - `test_land_cover_map.py`'s `counted` helper and `test_read_repair.py`'s
      `classes` helper call `repaired_areas`, which is deleted: `counted`
      becomes the literal count (1 in its one use), and `classes` reads
      through `land_cover_map.build([read.areas], [], classes, region)`,
      giving `(areas, repaired)`, with `region` a polygon whose read box
      holds `TRIANGLE`. Their figures are to stand; `@tester` confirms them
      on the change (not probed here: the plugin run kept the helpers).
    - The record's `land_cover_repaired` wording stays master's.
14. A spike under `--prefer` and in `read_polygons` (section 3, step 6;
    design round 11). Measured on tree `c5cd80c6` (scratch
    `a11_prefer_spike2.py`). In `test_cli_mesh_prefer.py`'s scene
    (`dem_origin()` (505,380, 8,317,320) in `MESH_CRS`; `read_polygons`'
    box from `domain_ring()` spans (9.1, −2,329.1) to (2,329.3, −9.3) from
    it), with `b = y0 − 2429.1`, the lake `Polygon([(x0 + 850, b), (x0 +
    1000, b + 1e-8), (x0 + 991.3, y0 − 1987.2), (x0 + 1000, b), (x0 + 1000,
    b − 200), (x0 + 850, b − 200)])`: invalid, 30,000 m², its spike from
    100 m below the box to 341.9 m into it. Written with `write_lake`:
    - `read_polygons(path, None, MESH_CRS, np.array(domain_ring()),
      lines="lines go to --lines")` raises `FeatureError` naming
      `lake.geojson`, the feature and "cannot be repaired" (red at
      `c5cd80c6`: GEOSException).
    - `rasputin mesh` with `corine(cover, "--prefer", f"lake={lake}")`
      (the scene's `mesh` helper) exits 2, the message holds "cannot be
      repaired" and the hint `--prefer`, and the output holds no
      "Traceback". `--reservoirs` reads through the same `read_polygons`,
      so the unit case covers it; no CLI case (its read box is the mesh's
      vertices, which the test would have to pin first).

*Functions the tests call* (tests 12 to 14; section 8 has the
signatures): `land_cover_map.build`; `feature_input._repaired(geometry,
xy, *, what=...)`, `feature_input.read_polygons`,
`feature_input.edges_of`, `feature_input.read_class_polygons`,
`feature_input.FeatureError`; `preferred.replace` (test 12's round-4 run
reaches it through the CLI); through the CLI, `rasputin mesh` by `test_cli_mesh_prefer.py`'s `mesh` helper and
`test_cli_mesh_edges.py`'s runs; and `test_feature_input.py`'s `one` and
`open_one`.

## 6. Size, timing, left out

Net production lines: +75 at `486b0e40` (`tools/count_loc.py 75ab2b4f 486b0e40`: 1,127 added, 1,052 removed); design rounds 10 and 11 added 16 of them. Tests change widely (about 270
`--features` invocations), all mechanical.

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
total is. A `MESH CHANGED` verdict is possible: the clean-up depends on
the polygons' vertex order (*As built*, RP8), and round 10 keeps the
polygons inside the read box as read but cuts those crossing its edge
another way than the branch does today. `@developer` runs it, names both
verdicts in the handback, and, if the mesh differs from master's, gives the
vertices moved and the farthest, as RP8's record does.

*Timing, several files.* One new term, only with a second `--land-cover`
file: the overlap query with its intersection areas, 0.23 s on Bygdin
(1,697 and 247 polygons, `read_source` over each file's `source_region`
box; 345 intersecting pairs, none sharing area; scratch probe on tree
`12616b91`, AC power). Cleaning the
two files as one coverage costs what cleaning each costs (`coverage_clean`
over the read region's 591 and 12 polygons: 0.249 s per file, 0.245 s
joined, best of three). One file adds nothing.

*Timing, the polygons' check.* Two new terms against master's
`--features` route, which checks nothing (section 3, *Cost and effect on
real data*): the whole `is_valid` of the polygons reaching the read box
and the intersection of those crossing its edge, 1,026 + 313 ms on
Numedalslågen on one thread, battery. On every core (the quick check's
`threads = [0]`, 10 on this machine) round 10 takes 539 ms against today's
107: 432 ms more, set by one polygon's `is_valid` (465 ms). Against the
stored baseline's 8.83 s that is 4.9 %, at the edge of the quick check's
5 % band (`tools/bench_quick.py`: at least 5 %, and more than 0.05 s), so
a `SLOWER` verdict on numedalslagen is likely (question 5). Lagan: 307 ms
more on every core, 2.1 % of 14.53 s. Sula and Bygdin add under 90 ms.
Against master's `--land-cover` runs, which check every polygon read on
one thread, round 10 is faster: 539 against 3,478 ms on Numedalslågen.

*Timing, the shares.* What `--land-cover` costs today, in the runs that
compute them (section 3): round 2's run measured 22.8 s of shares on Numedalslågen and 52 s on Lagan
(653,039 triangles on Numedalslågen), the polygon read with its validity
check included.

No `@perf` acceptance run (no refine or mesh code). Left out: fetching line
layers from the catalogue (the ideas entry above); shares in cut runs;
faster class shares; `subdivide`'s area loss (pipe item
`fractions-intersection`, section 3).

## 7. Questions for Ola

Answered, 2026-10-28: (1) the sources-file key `rules` becomes `classes`
(section 2). (2) Rivers, and roads the same way, without a land-cover map:
`--lines` (section 3). (3) No slowdown for a run that follows borders: it
computes no shares unless `--shares` asks (sections 3 and 6).

Answered, 2026-10-29: whether the shares' `subdivide` fix (pipe item
`fractions-intersection`, taken in by design round 11) stays in 111 or
leaves it. Ola, verbatim: "split". It leaves: 111 keeps the refusal of a
rebuild that raises; repair leaves the core in items 1.6 (the private
cleaning plugin) and 1.7 (`core-assumes-clean`, the core assumes clean
data).

Open:

4. *A county line inside a catchment.* Bygdin's two AR50 files meet on a
   7.6 km county line inside the catchment, with the same class on both
   sides all along it. Today's mesh follows that line (each file is
   cleaned on its own); with the files read as one map it does not, since
   the line is no border on the ground. Default: it is not followed. Say
   if county lines should stay edges.
5. *Half a second more on a Numedalslågen-sized run that follows class
   borders.* Repairing as master does means checking each polygon near the
   mesh whole. Your answer 3 asked that no border-following run gets slower.
   This check makes the quick check's Numedalslågen run about 0.43 s slower
   on every core (about 5 %, at the edge of its band); Lagan about 0.3 s
   (2 %). One CORINE polygon's check takes 0.47 s by itself, so more cores
   do not help. The refusing core that comes after the cleaning plugin
   needs the same check, unless the plugin's output is marked as checked.
   Default: accept it, and record the quick check's verdict on
   numedalslagen with this reason, without re-baselining.

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
def edges_of(land_cover: LandCoverMap | None, lines: Sequence[tuple[str, LineString]],
             choice: EdgeChoice, cleanup: EdgeCleanup, domain: DomainPolygon,
             vocabulary: EdgeVocabulary, *, region: Polygon) -> FeatureSet
```

In `tin_engine.land_cover_map`, as `test_land_cover_map.py` pins it:

```python
@dataclass(frozen=True)
class LandCoverMap:
    areas: tuple[tuple[int, BaseGeometry], ...]   # valid, repaired and cut in build; DEM's CRS
    outside: tuple[tuple[int, BaseGeometry], ...] # as read: no area in the read box
    classes: LandCoverClasses
    replaced: preferred.Replaced | None           # the --prefer cut; None without --prefer
    overlap: tuple[tuple[int, BaseGeometry], ...] # (index of the file kept, area cut from a later one)
    repaired: int                                 # invalid polygons reaching the read box

def build(base: Sequence[Sequence[tuple[int, BaseGeometry]]],      # one per --land-cover file
          preferred: Sequence[tuple[int, Sequence[BaseGeometry]]],  # as preferred.replace's chosen
          classes: LandCoverClasses, region: Polygon, *, fix: bool = True,
          outside: Sequence[tuple[int, BaseGeometry]] = (), threads: int = 1) -> LandCoverMap
```

`replaced` and `overlap` are section 3's `values`: the CLI measures them
inside the domain for `land_cover_prefer` and `land_cover_overlap`, since
`build` sees the read region only. A class raster gives
`LandCoverMap((), (), classes, None, (), 0)`; no `--land-cover` gives no map (the
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
- The repair is `build`'s (section 3), counted in `LandCoverMap.repaired`;
  its stderr line stays `_say_repaired("land cover", ...)`'s.
- Design round 11 (section 3, step 6):

  ```python
  # feature_input
  def _repaired(geometry: BaseGeometry, xy: npt.NDArray[np.float64], *,
                what: str = "a polygon") -> BaseGeometry | None
      # FeatureError "<what> near (x, y) in the DEM's CRS cannot be repaired: ..."
  ```

  `read_polygons`' signature is unchanged.
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
| `_Tally.source`, `_take`, `_linework`, `_bound`, `R_MIN`, `R_MAX`, `_mask` | go, with their tests; `_Tally`, `pre_clip` and `_runs` stay, for `edges_of` |
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
- Design round 9: CHANGES REQUESTED, 6ea1e369..07c11c28, 0 (`clip_by_rect` cuts valid polygons wrong when a box corner lies on an edge within float rounding: 3,920 of 5,505 probed, e.g. 5,697 m² against 484,303 m², and the corner test caught none; fix: cut every valid subset with `shapely.intersection` against the box, one step, no corner test; name how the polygon part of a mixed result is taken; time Lagan. Timings, 1/4 threads, round 9 against intersection: Numedalslågen 749/343 against 1,010/398 ms, Sula 44/39 against 45/38, Lagan 779/448 against 872/411, Bygdin 156/79 against 242/104; same areas on all four. Ola 2026-10-29 then moved repair out of the open core into a private cleaning plugin, so the per-run repair design of rounds 7 to 9 is withdrawn).
- Design round 10: CHANGES REQUESTED, 492f5d8c..692f005b, 0 (branch +59; round 9's corner finding closed: intersection right on 20,000 of 20,000 near-corner polygons at two origins, `clip_by_rect` control wrong on 2,033 and 8,454; repair leaves 0 invalid in 23,800 rebuilds, so nothing unrepaired reaches the overlap join, `--prefer` or the edges clean-up; blocking: an invalid polygon with a zero-width spike makes the rebuild raise a GEOS error, refused in `_fixed` but a traceback through `read_polygons` for `--prefer` (cli.py:2197) and `--reservoirs` (cli.py:1362), fix: raise `FeatureError` inside `_repaired` plus a `--prefer` test; `clip_by_rect` in `subdivide` still cuts valid polygons wrong for the shares and `areas_inside` (10 of 2,817 synthetic, up to 57,515 m²), steps 5 and 9 must say so or take pipe item `fractions-intersection` in; line 101 says `pre_clip` and `_runs` are deleted, still at feature_input.py:204 and :222).
- Code review round 6: APPROVED WITH WORDINGS, 75ab2b4f..486b0e40, +75 (+16 from design rounds 10 and 11; `_fixed`, `_cut`, `_repaired` and the CLI's `try` around `build` follow section 3, steps 1 to 6; probes: empty input, a polygon touching the box along a line goes to `outside`, a crossing multipolygon keeps its parts in order, 3,000 random near-corner polygons match `shapely.intersection` with none invalid, a crossing bowtie rebuilt to 20,000 m², a far bowtie to `outside` and not counted repaired, the thin spike refused with the planned message; 751 touched tests pass, full suite 11,366 passed; gates pass; three wordings applied in this commit; follow-up: `_fixed`'s docstring says "master's repair on read", to read "the repair on read").
