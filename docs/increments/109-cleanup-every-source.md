# Increment 109: the land-cover clean-up reaches every source, and says when it does not

**Status: designed** (`@architect`, 2026-10-27, on `5de1144d`, written against
increment 108's names). One PR, cut from master once 108 has merged, built
the same day. About 70 net production lines, no C++ file, no `@perf` run
(section 5).

## 1. What it fixes

Increment 108, *3. The finding that motivates increment 109*: on Sula
(73.27.0) Sol's own files ignored `--features-tolerance` (Sol's runs gave 4,167 triangles
with Sol's CORINE file and 29,395 with Sol's NVE lake file at 1, 20 and 50 m,
against 2,926 / 1,988 / 1,656 with the shipped `corine`). Three causes, each
probed there: a file without class codes skips the clean-up; a source that
does not cover the domain keeps every unshared border, so lakes are never
simplified; the run record writes the merged line and vertex totals into
each source's entry. Ola's requirement, through the brief: a missing field is
refused or given a stated default, never silently ignored.

## 2. Prior art: legacy and literature

*Literature.* The simplifier is increment 32's coverage simplification
(`docs/increments/32-landcover-simplify.md`, *2. Prior art: legacy and literature*): it fixes a
coverage's outer boundary by design. Closing a partial coverage with its
complement so that every border is shared is the device
`sub_catchments.simplify` already uses (a frame polygon around the union,
`src_python/tin_engine/sub_catchments.py@5de1144d:56-65`); here the
complement is the domain's uncovered rest, so the domain outline stays the
one fixed boundary. No departure from the method, no novelty claimed.

*Legacy.* `git grep -l -i -E 'simplif|coverage' legacy-archive -- legacy`
returns six files; the hits are mesh decimation by edge collapse
(`legacy/rasputin/mesh.py:102-108`, Lindstrom-Turk), not polygon borders.
Nothing carried.

*Ideas file.* "Repair CORINE once per area of interest" overlaps (it would
cache this stage's output); not taken up.

## 3. The design

1. **Every source's polygons are cleaned.** In `_rows`
   (`src_python/tin_engine/feature_input.py@5de1144d:396-400`), every
   Polygon or MultiPolygon of a non-`property` source that its rules do not
   drop goes to the clean-up, coded or not. The same-class merge groups an
   uncoded source by its class value text (the cover items carry that key
   beside the code). An uncoded polygon still labels no triangle: its
   `TerrainFeature` keeps `code` and `polygon` None.
2. **The uncovered rest is simplified with them.** In `_clean`, with the
   band on, `rest = domain − coverage_union_all(inside)`; its parts that an
   inward buffer of 1 cm empties are left out, the others join the polygons
   given to `simplify_borders` and are dropped after. The constant
   `REST_MIN_WIDTH_M = 0.01` is chosen, not derived: it assumes a metric CRS
   and was checked on Sula's 30 km² domain, where CORINE's rest is 18 parts
   of zero area, all under 1 cm wide (`/tmp/legend-split-probe/corine_rest.py`).
   A source that tiles the domain therefore adds nothing.
3. **`positional_accuracy_m` is required in every rules file**, refused when
   missing: "positional_accuracy_m missing: the source's positional accuracy
   in metres, which sets the default clean-up distances". Sol's two files
   are refused until they state one. `property` has none and is never
   cleaned (its help text says so).
4. **Default distances from the sources given.** The clean-up scale is the
   largest accuracy over the non-`property` sources, not over coded ones
   with CORINE's 100 m as fallback
   (`src_python/tin_engine/cli.py@5de1144d:1207`); with none, no clean-up
   runs.
5. **Unused options are said.** When the stage ran on no polygon, the record
   leaves the `features_repair_m`, `features_merge_same_class`,
   `features_tolerance_m` and `features_outline_snap_m` keys out; an
   explicitly given clean-up option prints one stderr line, e.g.
   "features: --features-tolerance not used: no polygon in any --features
   source" (a finding about the data).
6. **Per-source record counts.** Each source's `features` entry counts its
   own lines and vertices after clip and clean-up (`FeatureSet.features` are
   in source order; slice by `counts`); the merged totals move to one new key
   `features_joined`, "N lines, M vertices after joining shared edges and
   adding crossings", described in `run_record.py`.

7. **The model takes the file's names** (audit
   `docs/reviews/2026-10-27-concept-audit.md` on `ed0f6660`, finding 2.11):
   `MeshingRules.codes` becomes `code_system`, `string_codes` becomes
   `codes`, `otherwise` becomes `unlisted`, `accuracy_m` becomes
   `positional_accuracy_m` (no 0.0 fallback: it is required, item 3).
8. **One clean-up default** (audit finding 1.5): `FeatureRequest`'s four
   clean-up fields lose their all-off defaults; the derived distances
   (`derived_cleanup`) and merge on become a constructor,
   `FeatureRequest.derived(sources, cell)`, that the command calls, so a
   caller building a request itself gets the command's clean-up or states
   every value.

## 4. Prediction and tests

Prediction: shipped sources that tile the domain give byte-identical meshes
(gate T6); `corine-water` and every lake source now simplify. On Sula's 44
lakes (`/tmp/legend-split-probe/frame_probe.py`, `simplify_borders` on
master): 6,972 vertices alone at bands 1, 20 and 50 m; with the rest piece
3,947 / 502 / 362, area unchanged.

1. **Uncoded source cleaned.** Two lakes, circles of radius 100 m
   (`quad_segs=64`) at (300, 500) and (700, 500) in a 1,000 m square domain,
   uncoded rules: 514 vertices in. Measured with `simplify_borders` and the
   rest piece, repair 0.01 m: band 1 m → 54, 5 m → 28, 20 m → 15; each
   lake's area 31,412.773 m² kept. Without the rest piece: 514 at every band.
2. **Coded lakes:** the same fixture with codes, the same counts (cause 2).
3. **Tiling source:** two halves of the square, rest empty; the existing
   CORINE clean-up tests pass unchanged.
4. **Refusal:** a rules file without `positional_accuracy_m`, message as
   above.
5. **Record:** two sources give their own line and vertex counts and one
   `features_joined`; a `property`-only run with `--features-tolerance 5`
   prints the not-used line and no clean-up keys.
6. **Gate, no change for CORINE:** the quick check's `numedalslagen` and
   `lagan` meshes equal the sha256 values recorded in increment 108's gate 6
   (both cover their domain), one run each, and their wall times stay within
   10 % of 108's figures.
7. **Sula, reported, not gated:** Sol's two files with an accuracy added
   (100 m CORINE, 10 m lakes) at `--features-tolerance` 1, 20, 50: triangle
   counts reported to Ola beside the shipped `corine`'s.

## 5. Size

About 70 net production lines: the clean-up gate and merge key +10, the
rest piece +10, required accuracy and the default +5, the not-used line and
omitted keys +15, the per-source counts +20, the derived request +10, field renames +0. No C++; CORINE meshes are
unchanged by gate 6, so no `@perf` acceptance run.
