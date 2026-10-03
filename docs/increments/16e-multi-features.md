# Increment 16e — several `--features` sources in one mesh

Status: **designed** (`@architect`, 2026-09-29). Written before `@tester`, per
`docs/increments/README.md` step 1, on branch `increment16e-multi-features` off
`origin/master` (`7d9882d`). This is a CLI-and-plumbing increment: the internal
model already carries several sources; only the command line takes one today.

**Closes.** `ROADMAP.md`'s (new) 16e row, and the tail of 16b's Q4 ruling —
"one `--features` source per run ... a later CLI change" (16b, "Ruled by Ola",
2026-09-28, Q4). After this increment,

```sh
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --domain holleia.geojson --tolerance 10 \
    --features ../rasputin_data/corine2018_dtm10_utm33.gpkg \
        --features-layer corine2018 --features-map corine \
    --features holleia_parcel.geojson --features-map property \
    --out holleia.vtk
```

meshes one domain with **both** the CORINE land-cover boundaries and the
hunting-parcel boundary as constraints, each source read through its own class
map, CRS and layer. The motivating case is real: this session identified the
Holleia hunting terrain and meshed it, but `rasputin mesh --features` took one
source, so CORINE and the parcel could not both constrain one mesh. Ola wants
both.

**Not closed.** A per-edge record of which source an edge came from (that is
16d, proposed by Ola 2026-09-28, still to be designed). A class map read from a
file (16b's Q4, still open). Any change to how edges merge, are noded, are
clipped or are labelled — 16e adds sources; it changes no geometry.

## Ola's direction this design is built on

1. **Both constraints in one mesh** (2026-09-29): CORINE land cover *and* the
   hunting-parcel boundary, from two separate files, in the same run. The
   engine already merges every chain into one PSLG; the only thing missing is a
   command line that names more than one source.
2. **16b's Q4** (Ola, 2026-09-28): "one `--features` source per run plus the
   built-in maps", with the note that several sources in one run "are the API's
   (`FeatureRequest.sources` is a tuple) and a later CLI change (Q4)". 16e is
   that later CLI change; nothing about it reopens Q4's other half (the built-in
   maps).
3. **Each source keeps its own class map, CRS and layer** (16b R1, R2): a
   `FeatureSource` already carries `path`, `class_map`, `layer`, `crs`. CORINE
   needs `--features-map corine --features-layer corine2018`; the parcel needs
   the default `property` map and no layer. The two must not share a map.

## Prior art: legacy and literature

### Literature

**N/A, and none claimed.** This increment adds no algorithm. It lets the
command line name several feature sources; the geometry — noding, shared-edge
merge, clipping, labelling — is unchanged from 16b and 16c and carries their
citations (Chew 1989 for the CDT, Hobby 1999 and Goodrich et al. 1997 for the
snap-rounding noder, Shewchuk 1996 for the region labels). The pairing scheme
(R1) is a command-line convention, not a published method. No novelty is
claimed, so no search was made; if a later increment claims something new about
multi-coverage constraint sets it searches first, as 16b's list requires.

### Legacy

```sh
$ grep -rliE 'multi.?source|--features|append.*source|list.*of.*(source|repositor)' legacy/ | grep -v pycache | sort
(no output: nothing matched)
```

Nothing. The legacy had no `--features` flag and no notion of several vector
sources in one mesh: `legacy/rasputin/gml_repository.py` read one GML land-cover
delivery and made no constraint of it (`constraints()` returned `[]`, cited in
16b's Legacy section). There is nothing to carry across, and
`@migration-expert` is not needed.

## What is already there (verified on this base, `7d9882d`)

The plumbing below `cli.py` already models several sources. Read directly:

- `src_python/tin_engine/feature_input.py:133` — `FeatureRequest.sources:
  tuple[FeatureSource, ...]`.
- `feature_input.py:121-128` — each `FeatureSource` carries its own `path`,
  `class_map`, `layer` and `crs`.
- `feature_input.py:235-248` — `open_features` **loops over
  `request.sources`**, accumulates every source's features into one `_Tally`,
  and returns per-source tuples `crs=(...)` and `layers=(...)` alongside the
  flat `features` tuple.
- `feature_input.py:152-169` — `FeatureSet.crs` and `.layers` are **already
  per-source tuples**; only `cli.py` reads them as `[0]` today
  (`src_python/tin_engine/cli.py@7d9882d:836` and `src_python/tin_engine/cli.py@7d9882d:842`; on master 16e's per-source loop in `mesh`,
  `cli.py:905-911`, replaces both).
- `src_python/tin_engine/cli.py:1433` (`_dem_mesh`; line 1205 at `7d9882d`) —
  `_dem_mesh` hands `features.features`
  (a flat `tuple[TerrainFeature, ...]`, source-agnostic) to `start_chains`, and
  the noder (`build_pslg → node → triangulate`) merges every chain into one
  PSLG. It never sees which source a chain came from.

So the change is confined to `cli.py`: turn the four single options into
repeatable ones, pair them into a list of `FeatureSource`, build one
`FeatureRequest(sources=(...))`, and stop indexing the per-source tuples at `0`.

## The blueprint: data flow and boundaries

```
cli.mesh (composition root)
  │  --features PATH ...            repeatable  ->  list[Path]           [R1]
  │  --features-crs / -layer / -map repeatable  ->  list[str]            [R1]
  │
  ├─ _feature_sources(paths, crss, layers, maps) -> tuple[FeatureSource, ...]  [R1,R2]
  │     index-paired; degeneracy -> typer.BadParameter                  [R3]
  │     (today's _feature_source, generalised to N and reused per index)
  │
  ├─ open_features(FeatureRequest(sources=...), domain, dem_crs)        16b, UNCHANGED
  │     loops sources, one FeatureSet: flat features, per-source crs/layers
  │
  ├─ start_chains(domain, found.features, vocabulary)                   16b, UNCHANGED
  │     every source's chains in one PSLG; the noder merges shared and
  │     cross-source edges alike (R4)
  │
  ├─ _dem_mesh(...) -> Trimmed                                          UNCHANGED
  ├─ _land_cover(...) once, over found.features from every coded source 16c, R5 [R5]
  └─ file fields: one `features` / `features_crs` / `features_transform`
        line PER SOURCE, joined; one notice per distinct notice         [R6]
```

Boundaries, unchanged and each still checkable by the existing import tests:

- **Nothing new below `cli.py`.** `feature_input.py` already takes a
  `FeatureRequest` of many sources; no new module, no new import, no `_core`,
  no path handling moves. The I/O boundary (`CLAUDE.md` §2) is untouched.
- **Declarative still holds.** The CLI composes a frozen `FeatureRequest` from
  the flags and hands it to a pure `open_features`. 16e makes the tuple longer;
  it adds no state and no side effect.

## Rulings

### R1. The flags become repeatable, and pair by position

The four options become repeatable `list`s, and the Nth `--features` takes the
Nth of each paired option:

- **`--features PATH`** — repeatable; `list[Path]`, in the order given. Each is
  one source, read by suffix exactly as 16b R1 (`.geojson`/`.json`, `.gpkg`,
  `.gml`).
- **`--features-map NAME`** — repeatable; `list[str]`. The Nth applies to the
  Nth `--features`. **A missing entry** (the list is shorter than
  `--features`, including the empty list) defaults to `property`, exactly as
  the single option defaulted to `property` when omitted (16b R1;
  formerly `_feature_source`, `src_python/tin_engine/cli.py@7d9882d:1109`; on master `_feature_sources`,
  `cli.py:1317`).
- **`--features-crs TEXT`** — repeatable; `list[str]`. The Nth applies to the
  Nth `--features`. A missing entry is `None` (the source reads its own CRS:
  16b R1's file rule per suffix).
- **`--features-layer NAME`** — repeatable; `list[str]`. The Nth applies to the
  Nth `--features`. A missing entry is `None` (a single-`features`-table
  GeoPackage takes the only one; GeoJSON and GML refuse a layer: 16b R1,
  `src_python/tin_engine/cli.py@7d9882d:1114`; on master `_feature_sources`, `cli.py:1325`).

**Why positional pairing, not a compound flag.** Typer maps a repeated
`Option` to a `list`, so `--features a.gpkg --features b.geojson` is
`["a.gpkg", "b.geojson"]` with no new parsing (`python-development` skill §3,
Typer). The alternative — one flag carrying `path,map,crs,layer` as a
delimited string — would need its own parser, its own quoting rules for a WKT
CRS with commas, and its own error messages, and it hides the per-field
validation Typer gives for free. Positional pairing reads left to right the way
the flags are typed, and each `--features` with the `--features-map` written
next to it (as in the example above) pairs correctly because both lists grow in
type order.

**No shared-map convenience.** One `--features-map` is **not** broadcast to
several `--features`. It is a footgun: a user who writes `--features corine.gpkg
--features parcel.geojson --features-map corine` almost certainly wants `corine`
for the first and `property` for the second, and a silent broadcast would apply
`corine` to the parcel and refuse it (a parcel attribute is not a CLC code) or,
worse, mislabel it. The rule is one map per source or the default; a shared map
must be written out per source. (This is the only place a convenience was
weighed and rejected; R3 makes the mismatch that would tempt it an explicit
error.)

**Backward compatible.** A single `--features x --features-map corine`
still means one source with the `corine` map, because a one-element list pairs
with a one-element list. Every 16b/16c invocation keeps its meaning; the 16b
and 16c CLI suites pass unchanged (verified by @tester as a regression, R7).

### R2. Building the sources

`cli.py` gains `_feature_sources(paths, crss, layers, maps) ->
tuple[FeatureSource, ...]`, which replaces today's single `_feature_source`
(`src_python/tin_engine/cli.py@7d9882d:1104-1116`; on master `_feature_sources` is
`cli.py:1308-1342`). It:

1. Runs R3's degeneracy checks and raises `typer.BadParameter` on any failure.
2. For each index `i`, pads the shorter paired lists: `map = maps[i]` or
   `"property"`; `crs = crss[i]` or `None`; `layer = layers[i]` or `None`.
3. Reuses today's per-source validation unchanged: an unknown map name is a
   usage error naming it (`src_python/tin_engine/cli.py@7d9882d:1110`; on master `cli.py:1318`); a
   `--features-layer` on a non-`.gpkg` source is a usage error
   (`src_python/tin_engine/cli.py@7d9882d:1114`; on master `cli.py:1325`). The message names **which**
   `--features` it is about (its 1-based index and the path), because with
   several sources "applies to a .gpkg only" alone no longer identifies the
   source.
4. Returns `tuple(FeatureSource(path=p, class_map=CLASS_MAPS[m], layer=l,
   crs=c) ...)`.

`FeatureRequest(sources=_feature_sources(...))` then goes to `open_features`
as one call (`src_python/tin_engine/cli.py@7d9882d:1126` passed a one-element tuple; on master
`_open_features`, `cli.py:1353`), which already loops.

### R3. Degeneracy policy — every case a clear `typer.BadParameter`

Consistent with the existing guard "applies only with --features"
(`src_python/tin_engine/cli.py@7d9882d:726-728`), which is kept and generalised (on master the R3
length rule in `mesh`, `cli.py:750-763`). Each row raises
`typer.BadParameter` with the named `param_hint`:

| case | outcome |
|---|---|
| **a paired option is longer than `--features`** (e.g. two `--features-map`, one `--features`) | refused: `"N --features-map given for M --features; each pairs with one source by position"`, `param_hint` the offending flag. This subsumes the old "applies only with --features" when `--features` is absent (M = 0). |
| **a paired option with no `--features` at all** (`--features-map corine`, no `--features`) | the same rule with M = 0: refused naming the flag. Preserves `src_python/tin_engine/cli.py@7d9882d:726-728` (on master `cli.py:753-763`). |
| **a paired option shorter than `--features`** | **allowed**: the missing entries take their defaults (R1). This is the common case (`--features-crs` for one source only). |
| **`--features` given with no `--domain`** | refused, unchanged: `"needs --dem, --domain and --tolerance"`, `param_hint` `--features` (`src_python/tin_engine/cli.py@7d9882d:729-730`; on master `mesh`, `cli.py:764-765`). Applies whenever at least one `--features` is present. |
| **empty `--features` list** (no `--features` at all) | not a feature run: `found` stays `None`, the mesh has no constraints, exactly as today when `--features` is omitted. **But** any paired option present with no `--features` is the row above (refused). |
| **the same path given twice** (`--features x.gpkg --features x.gpkg`) | **allowed**: two sources reading the same file, e.g. once as `corine` and once as `corine-water`, or two different `--features-layer` of one GeoPackage. The noder merges the duplicated edges (R4); it is the user's stated intent, and refusing it would block a real use (two layers of one file). |
| **an unknown `--features-map` at any index** | refused naming the index and the name, listing the known maps (as `src_python/tin_engine/cli.py@7d9882d:1110`; on master `cli.py:1318-1323`). |
| **a `--features-layer` on a non-`.gpkg` source at any index** | refused naming the index and the path (as `src_python/tin_engine/cli.py@7d9882d:1114`; on master `cli.py:1325-1329`). |

The one length rule — *a paired option may not be longer than `--features`* —
covers the empty list, the map-without-features and the count-mismatch cases in
one check per paired option. A shorter paired list is always the defaults.

### R4. The noder across sources — nothing new needed

**Cross-source shared and crossing edges are already handled**, by the same
mechanism 16b measured within one source. Verified by reading the flow:

- `_dem_mesh` flattens every source's features into one chain list before the
  engine sees them (`_dem_mesh`, `cli.py:1433`, line 1205 at `7d9882d`:
  `start_chains(domain, features.features, ...)`), so the noder never distinguishes a within-source shared edge from a
  cross-source one. Both are two chains with the same endpoints.
- 16b R7 and its measurement M3 established that **the noder merges every
  doubled edge exactly by the node-id edge key and unions the masks** (5b's
  amended guarantee 14(a)). That merge is keyed on node ids, which are a
  function of the node set, not of the source. So a CORINE lake edge that
  coincides with the hunting-parcel boundary merges into one edge carrying the
  union of both masks — the same outcome 16b R4 describes for a CORINE edge
  between a lake and a forest.
- Two edges from different sources that **cross** (not coincide) are snapped
  and split at the crossing like any other crossing (05b step 4, 16b R6's
  "Crossing points ... No special case"). Nothing about crossing depends on the
  source either.

**One note for `@tester`, not a code change:** 16b's shared-edge measurement
(M2, M3) was on CORINE, whose neighbours share edges bit-for-bit after
reprojection. Two *different* sources reprojected from *different* CRSs need not
land on the same coordinates: a parcel boundary in EPSG:4326 and a CORINE edge
in EPSG:3035, both moved into UTM33 vertex by vertex, may differ by the
reprojection's own rounding even where they are "the same" line on the ground.
That is not a defect — they are different data — and the noder snaps both to its
1 mm grid, merging them iff they round to the same nodes, exactly as 16b's
crossing rule intends. **The design claims no exact merge across sources of
different CRSs**, only that the noder treats cross-source geometry by the same
snap-and-node rule as within-source geometry. R7's suite pins that two sources
whose edges coincide *in the DEM's CRS* merge; it does not assert coincidence
between reprojected sources.

### R5. Land-cover labelling over several sources

`_land_cover` (16c, `cli.py:1025-1040`; 949-964 at `7d9882d`) already gathers coded polygons across the
whole `FeatureSet`: `polygons = [(f.polygon, f.code) for f in found.features
if f.polygon is not None and f.code]` (`cli.py:1030`; 954 at `7d9882d`). Because `open_features`
merges every source's features into that one tuple, a coded source's polygons
are labelled whether or not other, uncoded sources are also present. **No
change.** Two subtleties, decided here:

- **The class map used for labelling** was read at `src_python/tin_engine/cli.py@7d9882d:855` as
  `CLASS_MAPS[features_map or "property"]` — the single map name (on master
  `mesh` takes the first coded source instead, `cli.py:931`). With several
  sources this is wrong: labelling must run when **any** source's map has codes
  (`ClassMap.codes` non-empty), and the `land_cover_codes` field text (16c R3)
  must name that coded map. **Ruling:** labelling runs iff at least one source's
  `class_map.codes` is non-empty; the field text names the first coded source's
  map. If two sources carry *different* code systems (`corine` and some future
  non-CORINE coded map), that is refused as a usage error — one mesh has one
  `land_cover_code` cell array and one code system (Default D1, for Ola). In
  practice every coded map today is CORINE (`corine`, `corine-water`,
  `clc18_kode` all set `codes = "CORINE Land Cover level-3 code"`), so mixing
  them is fine (same system) and the refusal never fires on today's maps.
- **`f.code` is `None` for an uncoded source's features**, so the parcel
  (map `property`) contributes no polygon to labelling — its boundary is a
  constraint that blocks the flood fill (16c R1) but names no class. That is
  correct: the parcel says *where*, CORINE says *what*.

### R6. What the file records — one line per source

Today `src_python/tin_engine/cli.py@7d9882d:834-854` writes one `features`, one `features_crs`, one
`features_transform` and at most one `features_notice`, reading `found.crs[0]`,
`found.layers[0]` and the single map name (on master the per-source loop in
`mesh`, `cli.py:898-927`). With several sources:

- **`features`** — one entry per source, joined (e.g. `; `-separated), each as
  today's text (`<name>[:layer], map <map>, <k> features, <c> chains, <v>
  vertices`) but per source. The per-source feature and chain counts are the
  new work: today `dem_run.feature_counts` is one `(chains, vertices)` pair for
  the whole run. **Ruling:** keep the whole-run totals for `feature_counts`
  (they drive the stderr "input vertices / noded vertices" line, which is about
  the merged PSLG), and derive the per-source `k features` from `found.features`
  grouped by source. Grouping needs a source index per feature, which
  `TerrainFeature` does not carry today. **Ruling (D2, for Ola):** rather than
  add a field to the frozen `TerrainFeature`, `open_features` returns a new
  per-source count tuple `FeatureSet.counts: tuple[int, ...]` (features kept per
  source), filled in the existing per-source loop (`feature_input.py:245-248`).
  This is ~3 lines in `feature_input.py` and keeps the file record honest per
  source. Chains and noded-vertex counts stay whole-run (they are properties of
  the merged PSLG, not of a source).
- **`features_crs`** — one per source, joined, from `found.crs` (already a
  tuple; stop indexing `[0]`).
- **`features_transform`** — one per source, joined, computed per source's own
  CRS as today (`transform_description(own, dem_crs)`), `none` where equal.
- **`features_notice`** — the **distinct** notices across the sources' maps,
  joined; a notice appears once even if two sources carry it (both CORINE maps
  carry the same Copernicus notice). This keeps the attribution correct without
  repeating it.
- **stderr** — `_open_features`'s "N features kept, ... dropped, ... clipped,
  ... empty" line already reports whole-run totals from the merged `FeatureSet`
  (`_open_features`, `cli.py:1358-1362`; 1131-1134 at `7d9882d`); it stays
  whole-run. The R-tree-scan warning (`cli.py:1363-1364`; 1136-1137 at
  `7d9882d`) already iterates `found.scanned`, which spans sources.

### R7. Tests, and the regression floor

No suite is named **invariant-critical** for 16e: it adds no algorithm, so
there is nothing whose topology decision a mutation round would probe (per
`docs/increments/README.md`, "Mutation testing is required only for the
invariant-critical suite"). What the red suite pins:

**`tests/python/test_cli_mesh_features.py`** (the 16b/16c CLI suite), new cases:

1. **Two GeoJSON sources, one mesh** over the synthetic `bumpy` DEM: a
   forest-and-lake CORINE-coded source and an uncoded `property` boundary
   source; the mesh carries both sets of constraint edges, the CORINE source's
   `land_cover_code`, and the boundary source's edges block the flood fill.
2. **Positional pairing**: `--features a --features b --features-map corine`
   (one map, two sources) gives source `a` the `corine` map and source `b` the
   default `property` map (the shorter list defaults; R1). `--features-crs`
   given for the second source only pairs with the second.
3. **The same GeoPackage twice** as two sources with two layers (or two maps),
   accepted; the edges merge (R4).
4. **Two sources whose edges coincide in the DEM's CRS** (both authored in the
   DEM's CRS so no reprojection rounding): the shared edge merges to one edge
   whose mask is the union of the two maps' bits (R4). This is the cross-source
   analogue of 16b's within-source M3 merge, pinned in CI on a hand-made
   fixture.
5. **Degeneracy**, one assertion per R3 row: two `--features-map` with one
   `--features` refused naming the flag and the counts; `--features-crs` with no
   `--features` refused (the generalised old guard); an unknown map at index 2
   refused naming the index; a `--features-layer` on a `.geojson` at index 2
   refused naming the index and path.
6. **Regression**: every existing single-`--features` case in this file and in
   `test_cli_mesh_landcover.py` passes unchanged (R1's backward compatibility).
7. **File record** (R6): with two sources the `.vtk` `features` field has two
   entries, `features_crs` two, and a single Copernicus `features_notice` even
   when both sources use CORINE maps; per-source `k features` counts are
   correct.

**`tests/python/test_feature_input.py`** (if `FeatureSet.counts` is added per
R6/D2): `open_features` on a two-source request returns `counts` of length 2
matching the features kept per source, and `crs`/`layers` of length 2 in source
order (this part already holds; the test pins it against regression).

The 16b/16c geometry suites (`test_cli_mesh_features.py::TestCommittedExtract`,
`test_landcover.py`) are unchanged: 16e touches no geometry.

## Invariants

- **I1 (pairing).** For sources `0..n-1`, source `i` is built from
  `features[i]`, `maps[i]` (or `property`), `crss[i]` (or `None`) and
  `layers[i]` (or `None`); a paired list longer than `features` is refused.
  Checked on `_feature_sources`'s output alone, no geometry.
- **I2 (merge across sources).** Two constraint chains from different sources
  with endpoints that round to the same nodes in the DEM's CRS become one noded
  edge whose mask is the union of the two chains' masks — 16b's I within-source
  merge, extended to cross-source. Checked on a hand-made two-source fixture in
  the DEM's CRS (R7 case 4). No exactness is claimed across differing source
  CRSs (R4).
- **I3 (labelling spans sources).** `land_cover_code` is computed from every
  coded source's polygons together; an uncoded source's boundaries block the
  flood fill but add no code. Checked by R7 case 1.
- **I4 (backward compatible).** A single `--features`/`--features-map`/etc.
  invocation produces a byte-identical mesh to the pre-16e code. Checked by the
  regression cases (R7 case 6).
- **I5 (file record per source).** The `.vtk`/`.ply` `features` and
  `features_crs` fields have one entry per source in source order;
  `features_notice` lists each distinct notice once. Checked by R7 case 7.

## Defaults for Ola to confirm

Each marked **Default (2026-09-29), for Ola to confirm**.

- **D1.** Two sources with **different code systems** in one mesh is refused
  (one mesh, one `land_cover_code` system). Two sources with the *same* system
  (both CORINE maps) is fine. Alternative: keep only the first code system's
  labels — rejected as silent data loss.
- **D2.** `FeatureSet` gains `counts: tuple[int, ...]` (features kept per
  source), so the file's per-source `features` line is honest, rather than
  adding a source index to the frozen `TerrainFeature`. ~3 lines in
  `feature_input.py`.
- **D3.** No shared-`--features-map` broadcast (R1): a map applies to one source
  or the source takes `property`. A shared map must be written per source.
  Rejected the broadcast as a footgun (R1).
- **D4.** The same file may be given twice as two sources (R3): two layers or
  two maps of one GeoPackage is a real use; the edges merge.
- **D5.** Whole-run counts (chains, noded vertices, "kept/dropped/clipped"
  stderr) stay whole-run; only per-source `features`, `features_crs`,
  `features_transform` and `counts` become per-source (R6). The merged PSLG has
  no per-source chain count.

## Not `@perf`'s

No file under `include/terrain/refinement/`, `include/terrain/mesh/`, or what
drives them changes — no C++ at all. The mesh core, the noder and the refiner
are byte-for-byte unchanged; 16e only lengthens the tuple of sources the CLI
composes before the engine runs. **`@perf`'s bench and thread sweep are not
required** (`docs/increments/README.md`, "Acceptance"). The one performance
observation worth recording is that N sources cost N reads plus one noder run
over the union — the noder's cost is already governed by the total edge count,
which 16b-0 made near-linear; adding a second small source (a parcel boundary)
to a large one (CORINE) adds its own edges to that total and nothing else.

## LOC

Production lines as `CLAUDE.md` §2 counts them (blank, comment, docstring and
raw-literal-body lines excluded; tests excluded):

| file | change | estimate |
|---|---|---|
| `cli.py` | four `Path\|None`/`str\|None` options → repeatable `list` options; `_feature_source` → `_feature_sources` (R2, R3 checks); `found.crs[0]`/`layers[0]`/single-map reads → per-source loops for the file fields (R6); coded-map detection over sources (R5) | ~40 |
| `feature_input.py` | `FeatureSet.counts` and its fill in the per-source loop (R6/D2) | ~5 |
| `open_features` / `_dem_mesh` / `start_chains` / the noder | none | 0 |
| **total** | | **~45** |

Worst case ~70 if the per-source file-field formatting (R6) is wordier than
expected. Well under the 700-line ceiling; one PR, one merge commit (never a
squash: `docs/increments/README.md`).

## Not in scope

- A per-edge record of the source an edge came from (16d, proposed by Ola).
- A class map read from a file (16b's Q4, still open).
- Any change to noding, the shared-edge merge, clipping or labelling geometry.
- Broadcasting one `--features-map` to several sources (D3, rejected).
- A compound single-flag syntax carrying path+map+crs+layer (R1, rejected).

## Review

**Round (citation-only fix), 2026-10-03.** Range `390b516..b213dc0` (6394c45, b213dc0). Verdict: CHANGES REQUESTED. LOC: 0 production lines (docs + test docstrings only; test ASTs unchanged modulo docstrings). Blocking: the branch shifts `## LOC` from 418-434 to 428-444, breaking `docs/retrospectives/2026-09-29-orchestrator-and-hooks-audit.md:229`; pin it as `16e-multi-features.md@390b516:418-434`. All 30 re-pointed or pinned citations read as quotations on master, `7d9882d` and `cb6f78b`. Not pushed; no CI.
