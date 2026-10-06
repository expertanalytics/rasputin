# Increment 16c — a land-cover class per triangle, in natural colours

Status: **built (green `0487ed0`), under review.** Designed by `@architect`
2026-09-29, night. Written before
`@tester`, per `docs/increments/README.md` step 1, on branch
`increment16c-landcover-labels` off master `b4847d7` (16b-1/2 merged, 22 not).
Ola was asleep while this was written; every question a design would
otherwise put to Ola is decided below as a **Default (2026-09-29), for Ola to
confirm** (collected under "Defaults for Ola to confirm").

**Closes.** `ROADMAP.md`'s 16c row, 16b's Q2 (a) as ruled on 2026-09-28
("a label per triangle by flood fill across unconstrained edges ..., written
as a cell field `land_cover` with the CLC code"), and Ola's request of
2026-09-29: *"For the example, I would like the land cover polygons colored in
a natural way for the vtk-vizualisation."* The example is increment 22's
Bygdin catchment meshed with `--features <Norway CORINE extract>
--features-map corine`. After this increment,

```sh
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --domain bygdin_reduced_t20.geojson --tolerance 10 \
    --features ../rasputin_data/corine2018_dtm10_utm33.gpkg \
    --features-layer corine2018 --features-map corine --out bygdin.vtk
rasputin palette corine --out corine_natural.json
```

writes a `.vtk` whose triangles carry their CORINE code in the cell array
`land_cover_code`, and a ParaView colour preset that paints those codes in
natural colours (forest green, rock grey, heath light green, bog brown, water
blue, glacier white-blue).

**Not closed.** "Water on exactly one side" as an edge property (16b's
Q2 note): the labels make it computable, but nothing here writes it.
Hydro-flattening lakes. A raster land-cover source (MapBiomas). Labels from a
`--features-map` without class codes (`property`; see R2). A colour table for
any code system other than CORINE.

## Ola's direction this design is built on

1. **16b's Q2 (a)**, ruled 2026-09-28: a separate increment, next after 16b;
   a label per triangle by flood fill across unconstrained edges, the method
   of Shewchuk's Triangle (`-A`), a cell field with the CLC code, about 80
   lines of pure Python over the trimmed mesh and the features.
2. **Natural colours, not the official CLC legend** (2026-09-29, via the
   main session's brief): the official legend colours peat bogs blue
   (`077-077-255`) and the sea almost white (`230-242-255`)
   (`Legend/CLC_legend.csv` in Ola's download, read for this design).
3. **The I/O boundary** (`CLAUDE.md` §2): file decoding and CRS stay in
   Python; nothing here goes near `_core`.
4. **Increment 13, ruling 5**: every cell array covers every cell; per-triangle
   data gives lines a fill value.
5. **Lean** (2026-09-29): red, green and review in one night; no throwaway
   implementations.

## Prior art: legacy and literature

### Literature

- **Regional attributes on a constrained triangulation.** Shewchuk,
  "Triangle: Engineering a 2D quality mesh generator and Delaunay
  triangulator", *Applied Computational Geometry*, LNCS 1148, 203-222, 1996.
  Triangle's `-A` takes a list of (seed point, attribute), finds the triangle
  containing each seed, and spreads the attribute across every edge that is
  not a segment; a triangle no seed reaches gets 0. Recalled (16b cited it the
  same way), not reread tonight.
  **What differs, and why.** R1 spreads first and looks up second: it finds
  the connected components of triangles across non-constraint edges with no
  seeds at all, then asks which input polygon each *component* lies in, using
  one point per component. Triangle's order needs one seed per polygon piece,
  chosen from the input, and that fails in three cases this project meets:
  a polygon cut into two pieces by the domain's boundary (the quarter circle
  does that), a polygon cut in two by a polyline running across it (a road),
  and a seed that lands in the wrong triangle because the polygon is thinner
  there than the noder's snap. Components-first has no seeds to place, so it
  has none of the three. The spread itself (constant across non-constraint
  edges, blocked by constraints, 0 where nothing applies) is Triangle's.
- **Connected components, vectorised.** Union-find with path compression:
  Tarjan, "Efficiency of a good but not linear set union algorithm", *J. ACM*
  22(2):215-225, 1975. The array form R1 uses, hooking the larger root onto
  the smaller and then pointer-jumping until every label is a root, is the
  hook-and-shortcut scheme of Shiloach and Vishkin, "An O(log n) parallel
  connectivity algorithm", *J. Algorithms* 3(1):57-67, 1982. Both recalled.
  **What differs:** Shiloach-Vishkin hooks conditionally and bounds its
  rounds by O(log n); R1 hooks every root to the smallest root it is joined
  to and loops to a fixed point, which is correct (labels only decrease, and
  hooking larger onto smaller cannot make a cycle) but claims no round bound.
  The acceptance records the rounds.
- **Point in polygon** is GEOS's (shapely 2 `STRtree.query(...,
  predicate="intersects")`). Nothing is claimed about it
  beyond what shapely documents. (This line said "over prepared polygons";
  as built, the polygons were the tree and were not prepared. Increment 30a,
  `docs/increments/30a-landcover-speed.md`, makes each polygon the query, which
  GEOS prepares.)
- **Where the labels came from before.** The legacy labelled cell centres,
  one Python point test per cell (see Legacy). The per-centre test is kept,
  as the **test oracle**, not as production (R1).
- **The CORINE nomenclature.** Kosztra, Büttner, Hazeu and Arnold, *Updated
  CLC illustrated nomenclature guidelines*, European Environment Agency, 2019:
  the 44 level-3 classes and their codes. Recalled. The class names used in
  R4 are the `LABEL3` column of `Legend/CLC_legend.csv` in Ola's download
  (`../rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/Legend/`),
  read for this design.
- **Natural colour for land cover.** Patterson and Kelso, "Hal Shelton
  revisited: designing and producing natural-color maps with satellite land
  cover data", *Cartographic Perspectives* 47:28-55, 2004: land cover drawn in
  the colours the ground has when seen from above (dark greens for conifers,
  tans and greys for bare ground, blues for water, white for ice), instead of
  a categorical legend. Recalled. R4's table follows that convention by hand;
  it is taste, and is not claimed to reproduce their palette.
- **The ParaView preset file.** ParaView's colour-map presets are a JSON list
  of objects; a categorical preset carries `"Name"`, `"Annotations"` (value,
  label, value, label, ...) and `"IndexedColors"` (r, g, b in 0-1, one triple
  per annotated value, in order), plus an optional `"NanColor"`. Recalled
  from ParaView's own presets file, **not checked against a running ParaView
  tonight** (none is installed). It is the acceptance's manual item, and the
  fallback if it fails is in R4.

**Novelty: none claimed.** Labelling a constrained triangulation's faces by
region is standard (Triangle, 1996). No search was needed because nothing new
is claimed; if a later increment claims a guarantee for region labels on
snapped input, it searches first (16b listed the queries).

### Legacy

```sh
$ grep -rlE 'cover_type|cover_color|def color\(' legacy/ | grep -v pycache | sort
legacy/rasputin/application.py
legacy/rasputin/globcov_repository.py
legacy/rasputin/gml_repository.py
legacy/rasputin/land_cover_repository.py
legacy/rasputin/tin_repository.py
legacy/rasputin/wfs_repository.py
legacy/tests/test_gml_repository.py
```

What matters, read directly:

- `legacy/rasputin/gml_repository.py:184-226`, `land_cover`: every cell
  centre tested against every polygon, one `shapely.Point` at a time (`#
  TODO: Move to C++ for speed!`), raising on a centre in no polygon.
  **Carried as the test oracle only** (R1, "the oracle"), with two changes:
  the test is vectorised, and a centre in no polygon is 0, not an error,
  because a general input need not cover the domain.
- `legacy/rasputin/gml_repository.py:69-117`, `LandCoverMetaInfo.color`: the
  official CLC RGB triples. **Not carried**: Ola asked for natural colours
  (direction 2).
- `legacy/rasputin/application.py:126-148`: writes the face fields
  `cover_type` (the code) and `cover_color` (RGB). **Carried:** a per-face
  code. **Not carried:** a per-face colour array (R4 says why) and the names
  (`cover_type` collides with nothing, but R3's `land_cover_code` says what it
  holds and sits next to 16b's `land_cover` bit).

`@migration-expert` is not needed: nothing numeric is ported.

## The blueprint: data flow and boundaries

```
cli.mesh (composition root)
  │
  ├─ feature_input.open_features(...)                     16b, extended [R5]
  │     a feature from a map with class codes also keeps
  │       code:    int      (its attribute value, e.g. Code_18 "322" -> 322)
  │       polygon: Polygon | MultiPolygon | None, in the DEM's CRS
  │     a coded polygon with no boundary inside the domain but covering it
  │       is kept, with no lines
  │
  ├─ _dem_mesh(...) -> Trimmed                            unchanged
  │     vertices (N, 3), triangles (T, 3), edges (E, 2), edge_masks
  │
  ├─ landcover.label_triangles(                           NEW, pure [R1, R2]
  │       vertices, triangles, edges,
  │       polygons=[(f.polygon, f.code) ...], margin=2 * snap spacing)
  │     -> CoverLabels(codes (T,) int32, regions, outside, overlapped, thin)
  │       regions(triangles, edges) -> (T,) component ids   numpy only
  │       one point per component -> STRtree query -> smallest-area rule
  │     timed as the phase "land cover"; one stderr line
  │
  ├─ io.vtk_legacy.write_vtk(..., triangle_codes=codes)   [R3]
  │     cell array land_cover_code: 0 on every LINES cell, the code on each
  │     triangle; FieldData string land_cover_codes names the code system
  ├─ io.ply.write_ply(..., face_codes=codes)              [R3]
  │
  └─ rasputin palette corine [--out FILE]                 NEW command [R4]
        palettes.paraview_preset(CORINE_NATURAL) -> JSON
```

Boundaries:

- `landcover.py` imports numpy and shapely, nothing first-party, never
  `_core`, never a path. It takes arrays and `(geometry, code)` pairs and
  returns arrays and counts. It is unit-testable on a hand-made mesh of three
  triangles.
- `palettes.py` is data plus one function; it imports nothing first-party and
  knows nothing about meshes. The CLI command and the acceptance's render
  script both read it, so there is one colour table.
- The writers stay pure (`bytes` out) and gain one optional array each.
- No C++ changes, no bindings, no new dependency (numpy, shapely and, for
  the picture only, the existing `viewer` extra's `vtk`).

## Rulings

### R1. The label: components first, one point per component

**Production.** For a trimmed mesh with triangles `T`, constraint edges `E`
and coded polygons `P`:

1. **Components** (`regions`). Each triangle's three edges as undirected keys
   `min(a, b) * N + max(a, b)` (int64; N is the vertex count). A key in `E`
   is blocked (`np.isin` on the same keys). The unblocked keys are sorted;
   two equal neighbours in the sorted order are one interior edge and join
   its two triangles. (A key appears at most twice: the trimmed mesh is a
   manifold.) Components by the array union-find of the literature section:
   `parent = arange(T)`; repeat { `ru, rv = parent[u], parent[v]`;
   stop if `ru == rv` everywhere; `np.minimum.at(parent, max(ru, rv),
   min(ru, rv))`; pointer-jump `parent = parent[parent]` until it is fixed }.
   The component id of a triangle is its component's smallest triangle
   index, so ids do not depend on the order edges are visited.
2. **One point per component.** Each triangle's incentre and inradius `r`
   (incentre `(a·A + b·B + c·C)/(a + b + c)` with `a = |BC|` and so on, `r =
   2·area/(a + b + c)`, in x and y). Per component, the triangle with the
   largest `r`, ties to the lowest triangle index; its incentre is the
   component's point.
3. **Which polygon.** `STRtree([p for p, _ in P]).query(points,
   predicate="intersects")`. A point in no polygon gives its component 0. A
   point in several gives the code of the one with the **smallest area**, ties
   to the smaller code (Default D2).
4. **Spread.** Every triangle gets its component's code.

**Why one point is enough, and when it is not.** A component is bounded by
constraint edges, and the input polygon boundaries lie within `δ` of the
constraint edges made from them. `05-noder.md` (fix 8, the one-way Hausdorff
form) bounds the other direction: every output edge derived from an input
segment lies within `h/√2` of it per noding round, `k·h/√2` after `k` rounds,
for snap spacing `h`. The output chain runs, inside that tube, from near one
end of the input segment to near the other, so it crosses the perpendicular
through any point of the segment within the same distance: `δ <= k·h/√2`
(0.7 mm per round at the default `h = 1 mm`). `margin = 2h` covers two
rounds. Refinement's points on constraint segments (20b's feet) are rounded
to doubles, which adds nanometres. This bound says when the lookup is exact;
the tests rest on the oracle, not on it.
The incircle of a triangle lies inside the triangle, and no constraint edge
enters a triangle, so the incentre is at least `r` from every constraint edge
and at least `r - δ` from every input boundary. With `margin = 2h > δ`, a
component whose largest `r` exceeds `margin` has a point that is strictly on
one side of every input boundary, so the lookup is exact up to GEOS's own
arithmetic. A component whose largest `r` is at most `margin` is counted as
`thin`: it is still labelled by its point, but the stderr line says how many
there were. (If the noder took more than two rounds, `δ` can exceed
`margin`; the count is then optimistic, and the oracle is the check.) On CORINE such a component is a sliver the snap made between two
boundaries that nearly coincide.

**Is the polygon the same object the constraints came from?** Yes, by R5: the
polygon tested is the feature's geometry moved into the DEM's CRS vertex by
vertex with the same transform that moved its chains, so its edges are the
straight lines the engine draws. The domain clip is linework on the chains
and is not repeated on the polygon; the triangles only exist inside the
domain anyway.

**The oracle (tests only).** Each triangle's **centroid** is tested against
`P` with the same predicate and the same overlap rule, for every triangle with
`r > 1.5·margin`. The centroid is at least `2r/3` from the triangle's sides
(its distance to side `a` is `2·area/(3a)`, and `a <= perimeter/2`), so for
those triangles it is further than `margin > δ` from any input boundary and
the oracle's answer is exact. Production and oracle must agree on every such
triangle. This is the legacy's per-centre test; it borrows the producer's
*predicate* (a point against the input polygons) and none of its records (no
components, no chosen points), per the `computational-geometry` skill's rule.
It is O(T) point queries, which is why it is the oracle and not production.

**Degeneracy and edge cases.**

| case | outcome |
|---|---|
| overlapping polygons (not CORINE; any GeoJSON) | smallest area wins, ties to the smaller code; the component is counted in `overlapped` (D2) |
| two polygons sharing a boundary that the noder merged (every CORINE neighbour pair) | one constraint edge between them; its two sides are different components, each labelled by its own point |
| boundaries closer than the snap, merged or not | at most a sliver component, labelled by its point, counted `thin` |
| a polygon with holes | shapely's test excludes the hole; the hole's components get whatever polygon fills it, else 0 |
| a lake inside a forest, the forest with a hole (CORINE) | lake components 512, forest components 31x |
| a lake inside a forest, the forest without a hole (hand-made input) | lake components the lake's code, by the smallest-area rule |
| a polygon clipped by the domain into several pieces | each piece is its own component and is labelled |
| a polyline across a polygon (a road through a forest) | two components, both the forest's code |
| a polygon covering the whole domain, no boundary inside | kept with no lines (R5); the single component gets its code |
| a triangle outside the domain | none exist: the engine's mesh is the domain's |
| a triangle dropped by `trim` (no DEM data) | not in the file, not labelled; components are computed after trimming |
| a component in no polygon | 0, counted in `outside` |
| no coded polygons at all (`property` map, or no `--features`) | no labelling, no array, no field (R2) |
| a constraint edge with mask 0 (the domain boundary, an unclassified edge) | blocks the spread like any constraint (it is in `E`) |

**Determinism.** The codes are a function of the triangle array, the
constraint edges and the polygons as a set: component ids are
smallest-index, the chosen point is largest-`r` then smallest-index, and the
overlap rule is order-free. 16b's I3 (rows in reverse order give a
bit-identical file) therefore extends to the new array with no new code; its
existing test covers it.

### R2. Where it runs: Python, over the trimmed mesh's arrays

In `src_python/tin_engine/landcover.py`, on numpy arrays, after `trim` and
before the writers. Not in the C++ core, because:

- the class codes and polygons are attribute data read in Python, in a CRS
  handled in Python; the core would need neither the polygons nor the codes,
  only the component pass, and that is about 25 lines of numpy;
- a C++ pass would need bindings, a rebuild, a C++ suite and, touching
  `include/terrain/mesh/`, `@perf`'s bench acceptance, which is not what a
  night is for.

**Cost at scale.** One sort of `3T` int64 keys, one `np.isin`, the union-find
rounds (each O(T)), and one point query per component, not per triangle. At
1 m on Bygdin (1.13 M triangles without features) that is a few million keys;
**not measured** (no prototype, per the retrospective agenda of 2026-09-29).
The acceptance records the phase time and the round count at 10 m and at 1 m.
If the phase costs more than the refine it follows at 1 m, moving `regions`
into the core is the next step, and the Python version stays as its oracle.

**Which maps label.** Only a class map that declares a code system does
(`ClassMap.codes`, R5): `corine` and `clc18_kode` (and `corine-water`, whose
kept features are the water polygons only, so land is 0). The `property` map's
values are vocabulary names, not codes; it writes no labels (Default D3).

### R3. What the files carry

**`.vtk`** (`io/vtk_legacy.py`, `write_vtk(..., triangle_codes=None)`):

- A cell array **`land_cover_code`**, `int`, one value per cell in the file's
  cell order: **0 on every `LINES` cell**, then each triangle's code (13's
  rulings 4 and 5: lines first, every array covers every cell, lines get the
  fill value). It is written in the `FIELD` block beside the per-feature
  arrays, not as the `SCALARS` block, which stays `feature_mask` (13, ruling
  5 is not reopened). The `FIELD` count includes it; with no feature arrays
  the block is written for it alone.
- **Not named `land_cover`**: 16b's vocabulary bit 7 is `land_cover`, and
  ruling 6 already writes a 0/1 cell array under that name. The writer refuses
  `triangle_codes` with a vocabulary that names a property `land_cover_code`.
- A dataset string **`land_cover_codes`** in `FieldData`: the code system,
  the attribute, the map, and what 0 means, e.g.
  `CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no
  polygon, and every constraint line`. `land_cover_codes` joins `RESERVED`.
- `triangle_codes` must have one entry per triangle and fit `int32`; else
  `ValueError`.

**`.ply`** (`io/ply.py`, `write_ply(..., face_codes=None)`): on the face file
only, a face property `int land_cover_code` after `vertex_indices`, and the
same text as a header comment `land_cover_codes ...`. Refused with `edges`
(the edge file has no faces). MDAL is expected to offer it as a face dataset
in QGIS; recalled, not measured, and not an acceptance item.

**stderr**, one line after the `trim` line:
`land cover: <R> regions, <O> outside every polygon, <V> in more than one,
<S> thinner than the snap`. The phase `land cover` appears in `--stats`'s
phase table through the existing clock; no other report row (16b's ruling
that a table in a run report is output, not data).

### R4. Natural colours: a ParaView preset, generated from one table

**Where the table lives.** `src_python/tin_engine/palettes.py`:
`CORINE_NATURAL: dict[int, tuple[str, str]]`, code → (CLC `LABEL3`, `#rrggbb`),
all 44 CLC level-3 classes plus 0, and `paraview_preset(table, name) ->
list[dict[str, object]]`. It ships in the package, is type-checked, and the
tests pin it. Not a JSON data file in the package: `wheel.packages` would
carry it, but a Python table is one object that the CLI and the render
script both import, with no file-location lookup.

**How a user gets it.** A new command, `rasputin palette corine [--out
FILE]`: the ParaView preset as JSON, to `FILE` or stdout. An unknown name is a
usage error listing the known ones. `mesh` does not write it beside the mesh
(Default D4): one more file per run, and a second name to keep from
colliding with `--out-edges` and `--stats`, for a table that never changes
between runs.

**How it is applied in ParaView** (the manual acceptance item):

1. Open `bygdin.vtk`; *Color By* `land_cover_code` (cell data).
2. Colour Map Editor → *Choose Preset* → *Import*, pick
   `corine_natural.json`, select "rasputin CORINE natural", *Apply*. Imported
   presets persist in ParaView's settings, so this is once per machine.
3. If the categorical mode is not switched on by the preset, tick
   *Interpret Values As Categories*.
4. Optionally *Save current colour map as default for arrays named
   `land_cover_code`* so every later file opens coloured.

Constraint lines carry 0 and draw in the colour for 0, a dark grey, which
reads as a class boundary. A code the table lacks draws in `NanColor`,
magenta, so it is seen rather than blended in.

**Rejected: a per-cell RGB array in the file** (the legacy's `cover_color`).
It needs no import, but it bakes one taste into every mesh file, adds three
bytes per cell, and a second viewer would want its own. It stays the fallback
if the preset fails the manual check.

**The table.** Natural colours, chosen by hand after Patterson and Kelso; the
official CLC colour is not used for any class. Hex is the value `@tester` pins
by family (Tests), not by exact triple.

| code | class (`LABEL3`) | colour | in the Norway extract |
|---|---|---|---|
| 0 | no polygon; constraint lines | `#3c3c3c` | (lines) |
| 111 | Continuous urban fabric | `#8c5f5a` | yes |
| 112 | Discontinuous urban fabric | `#a88d86` | yes |
| 121 | Industrial or commercial units | `#8e8a96` | yes |
| 122 | Road and rail networks and associated land | `#6e6e6e` | yes |
| 123 | Port areas | `#7d8796` | yes |
| 124 | Airports | `#a8a8a8` | yes |
| 131 | Mineral extraction sites | `#b49b78` | yes |
| 132 | Dump sites | `#8a7d64` | yes |
| 133 | Construction sites | `#bcae98` | yes |
| 141 | Green urban areas | `#86b86e` | yes |
| 142 | Sport and leisure facilities | `#a3cf7e` | yes |
| 211 | Non-irrigated arable land | `#e6d58c` | yes |
| 212 | Permanently irrigated land | `#d9cc6e` | |
| 213 | Rice fields | `#cdd89a` | |
| 221 | Vineyards | `#9c7a44` | |
| 222 | Fruit trees and berry plantations | `#a9ad5e` | yes |
| 223 | Olive groves | `#8f9a52` | |
| 231 | Pastures | `#b4d47a` | yes |
| 241 | Annual crops associated with permanent crops | `#dccf94` | |
| 242 | Complex cultivation patterns | `#d2c47c` | yes |
| 243 | Land principally occupied by agriculture, with significant areas of natural vegetation | `#bcc47e` | yes |
| 244 | Agro-forestry areas | `#a9b574` | |
| 311 | Broad-leaved forest | `#4f8f3f` | yes |
| 312 | Coniferous forest | `#1f5a2e` | yes |
| 313 | Mixed forest | `#357438` | yes |
| 321 | Natural grasslands | `#c2d68a` | yes |
| 322 | Moors and heathland | `#a7c47f` | yes |
| 323 | Sclerophyllous vegetation | `#8a9658` | |
| 324 | Transitional woodland-shrub | `#7ea65a` | yes |
| 331 | Beaches, dunes, sands | `#e9ddb2` | yes |
| 332 | Bare rocks | `#8f8f8f` | yes |
| 333 | Sparsely vegetated areas | `#b9b8a0` | yes |
| 334 | Burnt areas | `#4b3f3a` | yes |
| 335 | Glaciers and perpetual snow | `#eef6fb` | yes |
| 411 | Inland marshes | `#6e9470` | yes |
| 412 | Peat bogs | `#8a6642` | yes |
| 421 | Salt marshes | `#7f9f8c` | |
| 422 | Salines | `#d8d6cc` | |
| 423 | Intertidal flats | `#b3bdb3` | yes |
| 511 | Water courses | `#4c8ec4` | yes |
| 512 | Water bodies | `#3e7bb6` | yes |
| 521 | Coastal lagoons | `#5b93b3` | |
| 522 | Estuaries | `#5188b4` | yes |
| 523 | Sea and ocean | `#2b5d8e` | yes |

The "in the Norway extract" column is from
`sqlite3 -readonly ../rasputin_data/corine2018_dtm10_utm33.gpkg "select
code_18, count(*) from corine2018 group by code_18"`, run for this design: 34
codes, all present in the table. The unclassified codes 990, 995 and 999 do
not occur there and are left out; they draw in `NanColor`.

### R5. What `feature_input` keeps for labelling

- `ClassMap` gains `codes: str = ""`, the name of the code system. Non-empty
  means the map's attribute values are integer class codes and labels are
  written. `corine`, `corine-water` and `clc18_kode` set it to `CORINE Land
  Cover level-3 code`; `property` leaves it empty.
- `TerrainFeature` gains `code: int | None = None` and `polygon: Polygon |
  MultiPolygon | None = None`, both set only for a coded map: `code` is
  `int(str(value))` of the feature's attribute, refused (`FeatureError`,
  naming the feature and value) unless it is an integer in 1 .. 2³¹-1 — a
  list value is refused too; `polygon` is the feature's polygonal parts,
  moved into the DEM's CRS with the **same** transform that moved its chains
  (on the moved-first path it is `moved` itself; on 16b's geographic
  pre-clip path the whole polygon is moved once more). Lines have no polygon.
- **A coded polygon with no boundary in the domain** is today counted
  `outside` and dropped. For a coded map it is kept with `lines=()` if it
  covers the domain, tested by one point, `domain.polygon.point_on_surface()`:
  no boundary crosses the domain's interior, so the interior is wholly inside
  the polygon or wholly outside it. Without this, a catchment lying entirely
  inside one large CORINE polygon would be labelled 0 everywhere. It counts as
  kept, not clipped.
- Uncoded maps behave exactly as in 16b: no polygon is moved or kept, so no
  cost is added to them.
- Codes are checked only on the features the map keeps: under
  `corine-water`, an unlisted value is dropped before its code is read, as
  16b drops any unlisted value.

### R6. The CLI

- `mesh`: labelling runs whenever `--features` is given with a coded map and
  the output is `.vtk` or `.ply`. No new flag (Default D1).
- `palette`: R4.
- The `mesh` help text gains one sentence naming `land_cover_code` and
  pointing at `rasputin palette`.

## Invariants

- **I1 (spread).** For every interior edge that is not a constraint edge, its
  two triangles have the same code. Checked on the arrays alone, no geometry.
- **I2 (oracle).** For every triangle with `r > 1.5·margin`, the code equals
  the centroid oracle's (R1).
- **I3 (file shape).** `land_cover_code` has one value per cell, 0 on every
  `LINES` cell, and the triangle part equals `CoverLabels.codes`.
- **I4 (determinism).** 16b's reversed-rows run gives a bit-identical file.
- **I5 (area).** On a fixture whose polygons are known, the triangle area
  (x, y) per code equals the area of `polygon ∩ domain` per code within
  `margin × (boundary length inside the domain)` plus rounding.

## Defaults for Ola to confirm

Each marked **Default (2026-09-29), for Ola to confirm**.

- **D1.** Labels are always written for a coded map; no `--no-labels` flag.
- **D2.** Overlapping polygons: the smallest area wins, ties to the smaller
  code. The alternative, the last feature in source order wins (painter's
  rule), makes the result depend on row order, which 16b's I3 forbids.
- **D3.** The `property` map writes no labels. Giving its names integer ids
  would need an id table in the file; nobody has asked for it.
- **D4.** `rasputin palette corine` writes the preset on request; `mesh`
  does not write it beside the mesh.
- **D5.** Names: the cell array `land_cover_code`, the field
  `land_cover_codes`, the command `palette`, the preset "rasputin CORINE
  natural". 16b's Q2 said "a cell field `land_cover`"; that name is taken by
  bit 7's 0/1 array.
- **D6.** The colours in R4's table. Taste, and Ola's to change; the tests
  pin only the families (green forest, grey rock, light-green heath, brown
  bog, blue water, white-blue glacier).
- **D7.** No suite is named invariant-critical, so no mutation round: every
  labelling fixture is checked against the independent oracle (I2) and I1,
  which is what a mutation round would test (Ola's lean-brief rule).

## Tests for `@tester` (the red suite)

Lean: no throwaway implementation, no mutation round (D7). Two new files and
additions to three.

**`tests/python/test_landcover.py`** — pure, no `_core`, hand-made meshes:

- `regions`: two triangles sharing an unconstrained edge are one component;
  the same with the edge in `E` are two; a strip of 1 000 triangles with no
  constraints is one component (the loop reaches a fixed point); ids are the
  smallest triangle index; permuting `E`'s rows changes nothing.
- `label_triangles`: a square cut by a constrained diagonal with a polygon
  on each side gives each side its code; no polygons gives all 0 and
  `outside` = regions; two nested polygons give the inner one's code
  (smallest area) and count `overlapped`; a sliver component (largest `r` <=
  margin) is counted `thin`; `codes` is int32 of length T.
- The oracle as a test helper (`landcover_oracle(vertices, triangles,
  polygons, margin)`, per R1), used by the CLI tests below.

**`tests/python/test_cli_mesh_landcover.py`** — through `rasputin mesh` on the
synthetic `bumpy` DEM of `test_cli_mesh_features.py`, GeoJSON features with a
`Code_18` attribute and `--features-map corine`. For each fixture: I1, I2,
I3, and the expected codes named below.

1. **Two squares side by side** (311 and 512) sharing an edge: the triangles
   on either side of the shared edge carry 311 and 512; nothing is 0 inside
   the squares, the rest of the domain is 0.
2. **A hole**: a forest (312) with a square hole and nothing in it: the hole
   is 0.
3. **A lake inside a forest**: (a) the forest holed and the lake filling the
   hole: 512 and 312; (b) the forest not holed: still 512 in the lake
   (D2), and the stderr line says `in more than one` > 0.
4. **A polygon clipped by the domain into two pieces** (a band crossing a
   concave notch of the domain): both pieces carry its code.
5. **A road across a polygon** (a `LineString` feature with a code, under
   `corine`): both sides carry the polygon's code.
6. **A polygon covering the domain** with no boundary inside: every triangle
   carries its code; the feature is counted kept.
7. **Boundaries 0.1 mm apart** (two squares whose shared side is offset
   by less than the snap): both squares labelled correctly away from the
   seam (I2), and the run succeeds.
8. **Refusals**: a `Code_18` of `"forest"` under `corine` is a usage error
   naming the feature and value; `--features-map property` writes no
   `land_cover_code` and no `land_cover_codes`.
9. **`.ply`**: the face file has `property int land_cover_code` and the
   comment; the edge file has neither.
10. **VTK readback** (`importorskip("vtk")`, the `viewer` extra's CI step):
    `vtkPolyDataReader` sees `land_cover_code` with `GetNumberOfCells()`
    values, 0 on the line cells.

**Additions to `test_cli_mesh_features.py::TestCommittedExtract`** (the
quarter circle on the committed tile and extract, at 10 m): I1, I2, I3 on the
real mesh; I5 per code against shapely's `intersection` of each extract
polygon (moved to EPSG:25833) with the quarter circle, within 0.01 % of the
domain's area; the reversed-rows test (I4) passes unchanged.

**Additions to `test_io_vtk_legacy.py` and `test_io_ply.py`**: the array's
place and fill in ASCII and binary, the length refusal, the refusal of a
vocabulary naming `land_cover_code`, and the reserved field name.

**`tests/python/test_palettes.py`**: every one of the 34 Norway codes (R4's
list) and 0 is in `CORINE_NATURAL`; the 44 CLC codes are all there; families:
311-313 have g > r and g > b; 322 has g > r and g > b and is lighter than
312; 332 is grey (r = g = b); 412 has r > g > b (brown); every 5xx has b > r
and b > g; 335 has every channel >= 0xe0 and b >= r; `paraview_preset`
gives one object with `Name`, `Annotations` of length 2 × (entries),
`IndexedColors` of length 3 × (entries) in [0, 1], in the same order;
`rasputin palette corine` writes that JSON to stdout and to `--out`; an
unknown name is a usage error.

## Acceptance (`@perf`, AC power recorded)

On this branch, with the Bygdin reduced catchment taken from 22's branch
(no merge needed; 16c does not depend on 22's code):

```sh
git show increment22-autocatchment:docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson \
    > $SCRATCH/bygdin_reduced_t20.geojson
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --domain $SCRATCH/bygdin_reduced_t20.geojson --tolerance 10 \
    --features ../rasputin_data/corine2018_dtm10_utm33.gpkg \
    --features-layer corine2018 --features-map corine \
    --out $SCRATCH/bygdin_lc.vtk --stats bygdin_lc.stats.md
rasputin palette corine --out corine_natural.json
```

Evidence under `docs/benchmarks/2026-09-29/16c-bygdin/`:

1. **It opens in VTK with the array**: `vtkPolyDataReader` reads the file;
   `land_cover_code` is present with one value per cell and 0 on every line.
2. **Class shares against CORINE clipped to the catchment.** Triangle area
   (x, y) per code, as shares of the mesh's area, against 22's table
   (`bygdin/README.md` on 22's branch: 333 42.34 %, 332 21.32 %, 322
   17.00 %, 512 16.46 %, 335 2.44 %, 412 0.33 %, 142 0.10 %). **Pass:** the
   same seven codes, no 0 triangles, and every share within 0.01 percentage
   points. The expected gap is `δ` times the boundary length, well under a
   hectare; anything larger is a leak between regions.
3. **I1 and I2** on the Bygdin mesh (the oracle over all ~83 000 triangles).
4. **Time**: the `land cover` phase and the union-find rounds at 10 m, and
   one run at `--tolerance 1` with the same features for the scale. Recorded,
   not gated (R2).
5. **A picture**, `bygdin_landcover.png`: `render.py` in the evidence
   directory, using the `viewer` extra's `vtk` offscreen (`vtkRenderWindow`
   with `SetOffScreenRendering(1)`, `vtkWindowToImageFilter`, `vtkPNGWriter`;
   this design checked that the chain writes a PNG with vtk 9.7.0 in the
   repository's venv on this Mac), a map view from above with a
   `vtkLookupTable` in indexed mode built from `palettes.CORINE_NATURAL`, and
   a legend of the codes present. No new dependency.
6. **Manual, for Ola**: import the preset in ParaView and record the version
   and whether it switched to categorical by itself (R4, steps 2-3).

`@perf`'s bench and thread sweep are not required: nothing under
`include/terrain/refinement/` or `include/terrain/mesh/` or what drives them
changes; labelling runs after the mesh is finished.

## LOC

Production lines as `CLAUDE.md` §2 counts them (tests excluded):

| file | estimate | as built (green `0487ed0`) |
|---|---|---|
| `landcover.py` (`regions`, incentres, `label_triangles`, `CoverLabels`) | 75 | 83 |
| `palettes.py` (45-entry table, `paraview_preset`) | 60 | 64 |
| `feature_input.py` (`codes`, `code`, `polygon`, coded refusals, covering polygon) | 30 | 22 |
| `io/vtk_legacy.py` (array, fill, field, refusals) | 15 | 19 |
| `io/ply.py` (face property, comment, refusal) | 12 | 16 |
| `cli.py` (label call and phase, stderr line, fields and comments, `palette` command, help) | 40 | 45 |
| **total** | **~230** | **net 249 (262 added, 13 removed)** |

The estimate's worst case was about 320. Well under the 700-line ceiling; one
PR.

**As-built writer API.** `write_vtk(..., triangle_codes=None,
land_cover_codes="")`: the codes array and the `land_cover_codes` field text
are separate arguments, the text written only with the array.
`write_ply(..., face_codes=None)`.

## Not in scope

- An edge property "water on exactly one side" (computable from the labels).
- Labels for the `property` map, or for any map without integer codes (D3).
- Colour tables other than CORINE's; a colour array in the mesh file (R4).
- Moving `regions` into the C++ core (R2 says when it would be).
- Changing which array is the `.vtk` file's active `SCALARS` (13, ruling 5).
