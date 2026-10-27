# Increment 109: the land-cover clean-up reaches every source

**Status: designed, round 2** (`@architect`, 2026-10-27, on `1100fe2c`).
One PR from master, no C++ file, no `@perf` run (section 6). About 60 net
production lines.

## 0. What Ola asked

Ola, verbatim (2026-10-27): "Is there any way of relaxing the inner
constrains similarly as the catchment polygon?"; on the legend mix-up:
"Design, write, build, merge."; and "yes, push the Sula fix when it's
reviewed".

Scope for tonight, a default the main session took while Ola is away (his demo is tomorrow; not his ruling): the
simplification reaches every source (section 3, items 1, 2 and 4), the
per-source counts in the run record (item 5), and a stated default for a
missing positional accuracy (item 3), not a refusal, so Sol's demo files
keep running. Deferred to a later increment: section 4.

## 1. What it fixes

Ola's Sula case (NVE station 73.27.0, `--tolerance 20`): border
simplification (`--features-tolerance`) changes nothing for Sol's own files.
Re-run on master `1100fe2c` (`/tmp/cleanup-every-source-probe/run.sh`, one
source per run), triangles at `--features-tolerance` 1 / 20 / 50:

| source, rules | triangles | land-cover vertices after clean-up |
|---|---|---|
| Sol's CORINE file (no `code_system`) | 4,167 / 4,167 / 4,167 | no clean-up ran |
| Sol's NVE lake file (no `code_system`) | 29,395 / 29,395 / 29,395 | no clean-up ran |
| the lake file with codes added | 25,763 / 25,763 / 25,763 | 6,987 at every band |
| shipped `corine-water` | 1,930 / 1,930 / 1,930 | 220 at every band |
| shipped `corine` | 2,926 / 1,988 / 1,656 | 1,203 / 677 / 369 |

Three causes (increment 108, *3. The finding that motivates increment 109*,
re-checked at `1100fe2c`):

1. **A rules file without `code_system` is never cleaned.**
   `src_python/tin_engine/feature_input.py@1100fe2c:396-399`: a feature gets a
   class code only if the rules have a code system, and only coded polygons
   go to the clean-up. Repair, same-class merge, simplification and outline
   snap all skip it; the record still prints distances derived from CORINE's
   100 m (`src_python/tin_engine/cli.py@1100fe2c:1225-1227`).
2. **A border only one polygon has is never moved.** The simplifier fixes
   every border not shared by two polygons (increment 32, *7. Guarantees and
   the checks that enforce them*, item 4), so a lake shore, which no
   polygon of the same source shares, stays as read. This hits coded lake
   sources and the shipped `corine-water` too.
3. **The record repeats run totals per source**
   (`src_python/tin_engine/cli.py@1100fe2c:1324-1336`).

## 2. Prior art: legacy and literature

*Literature.* The simplifier is increment 32's coverage simplification
(`docs/increments/32-landcover-simplify.md`, *2. Prior art: legacy and
literature*); it fixes a coverage's outer boundary by design. Closing a
partial coverage with its complement so that every border is shared is the
device `sub_catchments.simplify` already uses (a frame around the union,
`src_python/tin_engine/sub_catchments.py@1100fe2c:56-65`); here the
complement is the domain's uncovered rest, so the domain outline stays the
one fixed boundary. No departure from the method, no novelty claimed.

*Legacy.* `git grep -l -i -E 'simplif|coverage' legacy-archive -- legacy`
returns six files; the hits are mesh decimation by edge collapse
(`legacy/rasputin/mesh.py`, Lindstrom-Turk), not polygon borders. Nothing
carried.

*Ideas file.* "Repair CORINE once per area of interest" overlaps (it would
cache this stage's output); not taken up.

## 3. The design (tonight)

In `feature_input.py`, `meshing_rules.py`, `cli.py` and `run_record.py`.

1. **Every source's polygons are cleaned, except `property`'s.** In
   `_Tally._take`, every Polygon or MultiPolygon of a source whose rules have
   `cleaned = True` (a new `MeshingRules` field, `bool`, default True; False
   only for `property`) and that its rules do not drop goes to the clean-up,
   coded or not. The same-class merge groups an uncoded polygon by its class
   value text. An uncoded polygon still labels no triangle: its
   `TerrainFeature` keeps `code` and `polygon` None. `--features-rules`' help
   gains the sentence "property: each feature's own edge names, never
   cleaned or simplified." (today it says nothing of it,
   `src_python/tin_engine/cli.py@1100fe2c:779-786`).
2. **The uncovered rest is simplified with them.** In `_clean`, with the band
   on, `rest = domain − coverage_union_all(inside)`; its parts narrower than
   `REST_MIN_WIDTH_M = 0.01` (a part that `buffer(-REST_MIN_WIDTH_M / 2)`
   empties) are left out, the others join the polygons given to
   `simplify_borders` and are dropped after. The constant is chosen, not
   derived: it assumes a metric CRS. Checked on `1100fe2c`: Sula with shipped
   `corine` (30 km²), 18 rest parts, all of zero area; numedalslagen and
   lagan (the quick check's land-cover cases, the largest checked), one part
   of zero area each; Sula's coded lakes, 16 parts, 27.7 km², all kept. A
   source that tiles the domain adds nothing, so its simplifier input, and
   its mesh, are unchanged.

   *Failure policy.* If `simplify_borders` refuses the input with the rest
   piece added (it raises `ValueError` saying why; the rest comes from two
   overlay calls, and the outline can cut a lake), the run is **refused**,
   as a refusal of the coded clean-up is today
   (`src_python/tin_engine/feature_input.py@1100fe2c:361-363`): "<file>:
   border simplification refused: <why>; --features-tolerance 0 meshes this
   source with its borders as read". Why not a fallback: a silent or
   warned fallback to unsimplified borders is the defect this increment
   fixes, and the triangle count would again not follow the option. Sula's
   lakes simplify without a refusal (section 5, prediction).
3. **A missing positional accuracy gets a stated default: the DEM cell.**
   `MeshingRules.accuracy_m` becomes `float | None` (None: the file states
   none; `property`'s is None and never read, item 1). In the command, a
   cleaned source with None counts as the DEM cell size, which leaves the
   clean-up scale `s = max(accuracy, cell)` at the cell, so it neither
   raises nor lowers the scale. Said once per such source on stderr and in
   the record key `features_accuracy_m` (one entry per source, `; `-joined):
   "sula-nve-lake-constraints.toml states no positional_accuracy_m; the DEM
   cell, 10 m, is used". A finding about the data, so it is said (Ola's rule
   on CLI warnings).
4. **Default distances from the sources given.** The clean-up scale is the
   largest accuracy over the cleaned sources (item 3's default included),
   not over coded ones with CORINE's 100 m as fallback.
5. **Per-source record counts.** Each source's `features` entry counts its own
   lines and vertices after clip and clean-up (`FeatureSet.features` are in
   source order; slice by `counts`); the merged totals move to one new key,
   `features_joined`, "N lines, M vertices after joining shared edges and
   adding crossings", described in `run_record.py`.

## 4. Later, not tonight

1. **The model takes the file's names.** Finding: the file says
   `code_system`, `codes`, `unlisted` and `positional_accuracy_m`
   (`src_python/tin_engine/meshing_rules.py@1100fe2c:92-100`), the model it
   becomes says `codes`, `string_codes`, `otherwise` and `accuracy_m`
   (`src_python/tin_engine/feature_input.py@1100fe2c:111-119`), so `codes`
   means two things across one conversion
   (`src_python/tin_engine/meshing_rules.py@1100fe2c:211-221`).
2. **One clean-up default.** Finding: `FeatureRequest`'s four clean-up fields
   default to all off (`src_python/tin_engine/feature_input.py@1100fe2c:142-148`),
   while the command derives 1, on, 50 and 5
   (`src_python/tin_engine/cli.py@1100fe2c:1225-1230`). Proposed: the fields
   lose their defaults and a constructor `FeatureRequest.derived(sources,
   cell)` gives the command's.
3. **Unused options said.** When the stage ran on no polygon, the record
   leaves the clean-up keys out and an explicitly given clean-up option
   prints "features: --features-tolerance not used: no polygon in any
   --features source".
4. **Required accuracy.** Whether a missing `positional_accuracy_m` is later
   refused instead of defaulted (item 3) is Ola's call.

## 5. Prediction and tests

Prediction: sources that tile the domain give byte-identical meshes (gate 7);
`corine-water` and every lake source now simplify. On Sula's 44 lakes
(`/tmp/cleanup-every-source-probe/frame_probe.py`, `simplify_borders` on
`1100fe2c`): 6,972 vertices alone at bands 1, 20 and 50 m; with the rest piece
3,947 / 502 / 362, area unchanged (0.0 m²).

1. **Uncoded source cleaned.** Two lakes, circles of radius 100 m
   (`quad_segs=64`) at (300, 500) and (700, 500) in a 1,000 m square domain,
   uncoded rules: 514 vertices in. Measured on `1100fe2c` with
   `simplify_borders` and the rest piece, repair 0.01 m: band 1 m → 54,
   5 m → 28, 20 m → 15; each lake's area 31,412.773 m² kept. Without the rest
   piece: 514 at every band.
2. **Coded lakes:** the same fixture with codes, the same counts (cause 2).
3. **Tiling source:** two halves of the square; rest empty; the existing
   CORINE clean-up tests pass unchanged.
4. **The 1 cm filter.** The 1,000 m square; a source of two polygons, the
   lake at (300, 500) as in test 1 and `box(0, 0, 1000, 1000 − w)` less the
   lake: 519 vertices. The rest is one strip, 5 vertices. `w = 0.008`: area
   8.000 m², left out, the simplifier gets 2 polygons; `w = 0.012`: area
   12.000 m², kept, it gets 3. Band 5 m, repair 0.01 m: the source's 519
   vertices become 33 in both, lake area 31,412.773 m² kept
   (`/tmp/cleanup-every-source-probe/sliver.py`).
5. **Refusal message:** `simplify_borders` patched to raise
   `ValueError("x")` with an uncoded source: the run is refused with the
   message of item 2, naming the file.
6. **Record:** two sources give their own line and vertex counts and one
   `features_joined`; an uncoded source with no `positional_accuracy_m` on a
   10 m DEM gives item 3's sentence on stderr and in `features_accuracy_m`,
   and the derived tolerance 5 m (s/2 with s = 10).
7. **Gate, CORINE meshes byte-identical and no slower.** The quick check's
   `numedalslagen` and `lagan` cases (`docs/benchmarks/quick/cases.toml`, plus
   `--out X.vtk --binary`), sha256 of the `.vtk` must stay
   `1b077a1a5ce8170f7cfebe9a01cdcfe551ecc50d22048f2e73f95a83b32b4705` and
   `4c9f885a5948444e1c2c9a695faedbf3e163d59a3490d6a4f12f793123d83965`
   (increment 108's gate 6; re-measured on `1100fe2c`, same hashes).
   Speed: wall time, one run each after, within 11 % of the same power
   state's figures (10 % noise plus the new term). Before, on `1100fe2c`,
   battery, 2026-10-27: numedalslagen 10.58 / 10.42 s, lagan 17.20 / 17.11 s;
   on AC, increment 108's 10.41 / 9.98 s and 16.26 s. The new term, timed in
   the clean-up (`/tmp/cleanup-every-source-probe/rest_probe.py`): 0.075 s
   (numedalslagen) and 0.186 s (lagan), about 1 %; the rest is one part of
   zero area in both, so it is left out and the simplifier's input is as
   today.
8. **Sula, reported, not gated:** Sol's two files as they are (no accuracy:
   item 3's default) and `corine-water`, at `--features-tolerance` 1, 20,
   50: triangle counts reported to Ola beside section 1's table.

## 6. Size

About 60 net production lines: the clean-up gate, `cleaned` flag and merge
key +10, the rest piece and its filter +12, the refusal message +3, the
accuracy default, its stderr line and record key +12, the scale over cleaned
sources +3, the per-source counts and `features_joined` +20. No C++; CORINE
meshes unchanged by gate 7, so no `@perf` acceptance run. Test helpers that
build `MeshingRules` with `accuracy_m=0.0` for "none" may need `@tester`'s
amendment in the red step.

## Review

- Design round 1: CHANGES REQUESTED, 1100fe2c..b243f7e5, 0 net production lines (Ola's words not quoted, no 1 cm filter test, no failure policy, a false help-text claim).
- Design round 2: APPROVED WITH WORDINGS, ..606c90b8, 0; the one wording applied in this commit.
