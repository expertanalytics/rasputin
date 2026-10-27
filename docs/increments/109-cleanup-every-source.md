# Increment 109: the land-cover clean-up reaches every source, and says when it does not

**Status: designed** (`@architect`, 2026-10-27, on `1100fe2c`). One PR from
master, no C++ file, no `@perf` run (section 5). About 70 net production
lines.

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

## 3. The design

All in `feature_input.py`, `meshing_rules.py`, `cli.py` and `run_record.py`.

1. **Every source's polygons are cleaned.** In `_Tally._take`, every Polygon
   or MultiPolygon of a non-`property` source that its rules do not drop goes
   to the clean-up, coded or not. The same-class merge groups an uncoded
   polygon by its class value text. An uncoded polygon still labels no
   triangle: its `TerrainFeature` keeps `code` and `polygon` None.
2. **The uncovered rest is simplified with them.** In `_clean`, with the band
   on, `rest = domain − coverage_union_all(inside)`; its parts that a 1 cm
   inward buffer empties are left out, the others join the polygons given to
   `simplify_borders` and are dropped after. `REST_MIN_WIDTH_M = 0.01` is
   chosen, not derived: it assumes a metric CRS. Checked on `1100fe2c`: Sula
   with shipped `corine` (30 km²), 18 rest parts, all of zero area and under
   1 cm wide; numedalslagen and lagan (the quick check's land-cover cases,
   the largest checked), one part of zero area each; Sula's coded lakes, 16
   parts, 27.7 km², all kept.
   A source that tiles the domain adds nothing, so its input to the
   simplifier, and its mesh, are unchanged.
3. **`positional_accuracy_m` is required in every rules file**, refused when
   missing: "positional_accuracy_m missing: the source's positional accuracy
   in metres, which sets the default clean-up distances". Today it is required
   only with `code_system` (`src_python/tin_engine/meshing_rules.py@1100fe2c:121-122`).
   Every shipped file has one; Sol's two files are refused until they state
   one. `property` has none and is never cleaned (its help text says so).
4. **Default distances from the sources given.** The clean-up scale is the
   largest accuracy over the non-`property` sources, not over coded ones with
   CORINE's 100 m as fallback; with none, no clean-up runs.
5. **Unused options are said.** When the stage ran on no polygon, the record
   leaves the `features_repair_m`, `features_merge_same_class`,
   `features_tolerance_m` and `features_outline_snap_m` keys out; an
   explicitly given clean-up option prints one stderr line, e.g.
   "features: --features-tolerance not used: no polygon in any --features
   source" (a finding about the data).
6. **Per-source record counts.** Each source's `features` entry counts its own
   lines and vertices after clip and clean-up (`FeatureSet.features` are in
   source order; slice by `counts`); the merged totals move to one new key,
   `features_joined`, "N lines, M vertices after joining shared edges and
   adding crossings", described in `run_record.py`.
7. **The model takes the file's names.** Finding: the file says
   `code_system`, `codes`, `unlisted` and `positional_accuracy_m`
   (`src_python/tin_engine/meshing_rules.py@1100fe2c:92-100`), the model it
   becomes says `codes`, `string_codes`, `otherwise` and `accuracy_m`
   (`src_python/tin_engine/feature_input.py@1100fe2c:111-119`), so `codes`
   means two different things across one conversion
   (`src_python/tin_engine/meshing_rules.py@1100fe2c:211-221`).
   `MeshingRules` takes the file's four names; `positional_accuracy_m` has no
   0.0 fallback (item 3). No aliases.
8. **One clean-up default.** Finding: `FeatureRequest`'s four clean-up fields
   default to all off (`src_python/tin_engine/feature_input.py@1100fe2c:142-148`),
   while the command derives 1, on, 50 and 5
   (`src_python/tin_engine/cli.py@1100fe2c:1225-1230`), so a caller building a
   request itself gets no clean-up without a word. The fields lose their
   defaults, and `derived_cleanup` plus merge-on becomes a constructor,
   `FeatureRequest.derived(sources, cell)`, that the command calls.

## 4. Prediction and tests

Prediction: sources that tile the domain give byte-identical meshes (gate 6);
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
4. **Refusal:** a rules file without `positional_accuracy_m`, message as above.
5. **Record:** two sources give their own line and vertex counts and one
   `features_joined`; a `property`-only run with `--features-tolerance 5`
   prints the not-used line and no clean-up keys.
6. **Gate, CORINE meshes byte-identical and no slower.** The quick check's
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
7. **Sula, reported, not gated:** Sol's two files with an accuracy added
   (100 m CORINE, 10 m lakes) and `corine-water`, at `--features-tolerance`
   1, 20, 50: triangle counts reported to Ola beside section 1's table.

## 5. Size

About 70 net production lines: the clean-up gate and merge key +10, the rest
piece +10, required accuracy and the default +5, the not-used line and omitted
keys +15, the per-source counts +20, the derived request +10, field renames
(13 sites in `src_python`) +0. No C++; CORINE meshes unchanged by gate 6, so
no `@perf` acceptance run. Test helpers that write uncoded rules files
(`tests/python/legend_fixtures.py`) and that build `FeatureRequest` without
clean-up fields need `@tester`'s amendment in the red step.
