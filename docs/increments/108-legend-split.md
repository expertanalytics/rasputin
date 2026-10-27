# Increment 108: a legend is the colour scheme; meshing rules are a file of their own

**Status: designed** (`@architect`, 2026-10-27, on `5de1144d`). One PR, branch
`worktree-legend-split`, about 110 net production lines, no C++ file, no
`@perf` run (section 6). Ola ordered it built and merged today, 2026-10-27.

## 1. What Ola asked

Ola, 2026-10-27: "WHAT?! The legend is the color scheme!!" and "Legend and
meshing rules are so different it hurts you have mixed them." Today one
legend file (`src_python/tin_engine/data/legends/*.toml`, read by
`src_python/tin_engine/legends.py`) holds both the colours and labels of a
source's classes and the rules that turn the source into constraint lines.
This splits them into two files, with no alias for any old name.

## 2. Prior art: legacy and literature

*Literature.* Cartography and GIS keep a data set's styling apart from its
data and schema: OGC *Symbology Encoding Implementation Specification* 1.1.0
(OGC 05-077r4, 2006) styles features in a document separate from the feature
type, as QGIS keeps a layer style in a `.qml` file or the `layer_styles`
table. The FAO LCCS (Di Gregorio and Jansen, 2000, chapter "From
classification to legend") is the source of the word "legend" in increment 58;
there the legend is the classes of one data set as shown, which fits the
colour-scheme file below. The meshing rules (which borders become which
constraint edges, and at what accuracy) have no counterpart in either; they
are ours. No novelty is claimed.

*Legacy.* `git grep -l -i -E 'legend|palette|colou?r' legacy-archive -- legacy`
returns ten files; the land-cover one, `legacy/rasputin/land_cover_repository.py`,
gives colours from a separate `color(cls, *, land_cover_type)` method (line
16), apart from the classification: the split this increment makes. Nothing
is carried across; the colours come from today's shipped files.

*History.* Before increment 41 the two were apart: "class maps" (rules) and
`palettes.py` (colours). Increment 41 merged them into one file; increment 58
renamed it "legend" (Ola: "Class-map is a very bad name for this") and
already noted the mismatch ("our legend also says which edge names ... a
cartographic legend has no such part"). The name "class map" is ruled out,
so the rules take Ola's own words: **meshing rules**.

*Ideas file.* No entry takes this up. "Product review: readable code docs, and
general applicability" overlaps (its rulings made these files data, increment
41).

## 3. The two causes of the Sula result, and the record bug

Sol's run on Sula (73.27.0, `--tolerance 20`) used two user files without
`code_system` (`sula-corine-constraints.toml`, value-and-edges classes;
`sula-nve-lake-constraints.toml`, no classes, `unlisted = "land"`).
Re-run from master's venv into `/tmp/legend-split-probe` (one source per run):

1. **An uncoded file is never cleaned.** `src_python/tin_engine/feature_input.py@5de1144d:396`
   gives a feature a class code only `if legend.codes`; line 397 makes it a
   polygon only with a code; line 399 sends only those to the clean-up. So
   repair, same-class merge, simplification and outline snap all skip it, and
   the run record still prints the three distances, derived from CORINE's
   100 m (`src_python/tin_engine/cli.py@5de1144d:1207`, `default=100.0`). Probe: Sol's CORINE file
   gives 53 lines, 1,250 vertices at `--features-tolerance` 1, 20 and 50;
   the built-in `corine`, 1,203 / 677 / 369 land-cover vertices after
   clean-up.
2. **A source that does not tile the domain is never simplified, even coded.**
   The simplifier keeps every border not shared by two polygons
   (`docs/increments/32-landcover-simplify.md`, *7. Guarantees and the checks that enforce them*, item 4),
   so a lake shore stays as read. Probe: the lake file rewritten with codes
   (`lakes-coded.toml`, accuracy 10 m): 9,818 vertices after the clip,
   6,987 after clean-up at 1, 20 and 50 m alike. The shipped `corine-water`
   has the same defect.
3. **The record repeats run totals per source.** `src_python/tin_engine/cli.py@5de1144d:1307-1316`
   writes the merged PSLG's lines and vertices into every source's entry
   (Sol's record: "112 lines, 8180 vertices" twice). Same code path; fixed
   here.

## 4. The design

**Legend file** (`data/legends/<name>.toml`): the colour scheme only.
`name`, `preset` (the ParaView preset name), optional `[no_class]` (label,
colour of code 0), and `[[class]]` with `code`, `label`, `colour`. Every key
required, nothing else allowed; a repeated code or code 0 in `[[class]]` is
refused. A legend may list codes no rules file uses.

**Meshing-rules file** (`data/meshing_rules/<name>.toml`): `name`,
`attribute`, **`positional_accuracy_m` required in every file** (refused when
missing: "positional_accuracy_m missing: the source's positional accuracy in
metres, which sets the default clean-up distances"), `unlisted`, optional
`code_system`, `notice`, `licence`, `[edges]`, `[[class]]` (`code` or `value`,
`water`, `category`, `edges`), `[[split]]` (sub-classes `value`, `code` only),
`codes`, `vocabulary.extra`, `[prefer]`, and **`legend`**: a shipped legend
name or a path, relative to the rules file. With `code_system`, `legend` is
required and every class and sub-class code must be in it (refused naming the
missing codes); without, `legend` is refused. `label` and `colour` are refused
in a rules file.

**Python.** `legends.py` keeps only the legend: a frozen `Legend` (`name`,
`preset`, `colours: Mapping[int, tuple[str, str]]`), `load_legend(value)`.
New `meshing_rules.py` takes today's model, loader, `derived_cleanup` and
`run_vocabulary`; the frozen model in `feature_input.py` is renamed
`MeshingRules`, its `palette` and `preset` fields replaced by
`legend: Legend | None` (composed, not copied). `FeatureSource.legend`
becomes `.rules`. Readers of colours (`palette`, the GeoPackage writer's
classes and style, `known_codes` at `src_python/tin_engine/cli.py@5de1144d:1916`) read
`rules.legend.colours`, or for `known_codes` the rules' own class codes.
`property` stays built in, has no accuracy and no legend, and is never
cleaned (help text says so).

**The clean-up, fixed.** (a) Every polygon of a non-`property` source that is
not dropped goes to the clean-up, coded or not; the same-class merge groups
an uncoded source by its class value. An uncoded polygon still labels no
triangle. (b) Before simplifying, `_clean` adds the domain's uncovered rest
(`domain − union`) as one more polygon and drops it afterwards, so every
shore is a shared border and only the domain outline stays fixed. Rest parts
an inward buffer of 1 cm empties are left out: on Sula, CORINE's rest is 18
parts of zero area, all under 1 cm wide (`corine_rest.py`), so a source that
tiles the domain adds nothing and its mesh is byte-identical. The constant
assumes a metric CRS and was checked on Sula's 30 km² domain. (c) The default
accuracy is the largest over the non-`property` sources; with none, no
clean-up runs. (d) When the stage ran on no polygon, the record leaves the
four `features_*` clean-up keys out, and an explicitly given clean-up option
prints one stderr line: "features: --features-tolerance not used: no polygon
in any --features source" (a finding about the data, so allowed).

**Record.** Each source's entry counts its own lines and vertices after
clip and clean-up (from `FeatureSet.features`, which are in source order,
sliced by `counts`); the merged totals move to one new key `features_joined`
("N lines, M vertices after joining shared edges and adding crossings").

**Names.** No aliases; old ones refused:

| Today | After |
|---|---|
| `mesh --features-legend` | `--features-rules` |
| `mesh --land-cover-legend` | `--land-cover-rules` |
| MCP `features_legend`, `land_cover_legend` | `features_rules`, `land_cover_rules` (derived from the options; `mcp/server.py`'s instruction text by hand) |
| sources-file key `legend` | `rules` |
| `sources` fit field `legend` | `rules` |
| `rasputin palette NAME` (any legend) | a legend name or legend file only |

**Shipped files.** Legends: `corine` (one file: the three CORINE files'
colours are equal, checked with `load`), `mapbiomas-c11`, `ar50` (sub-class
colours as plain classes). Rules: `corine`, `corine-water`, `clc18_kode`
(each `legend = "corine"`), `mapbiomas-c11`, `ar50`, accuracies unchanged
(100, 100, 100, 30, 2 m).

No `--legend` override option in this PR: a user picks colours by naming a
legend file in a rules file.

## 5. Prediction and tests

Prediction: every shipped coded source that tiles the domain gives a
byte-identical mesh; Sol's two files are refused until they state an accuracy;
with one stated, both simplify. Probe of (b) on Sula's 44 lakes (`frame_probe.py`):
lake vertices 6,972 alone at bands 1, 20, 50; with the rest piece 3,947 /
502 / 362, area unchanged.

Tests (red first):

1. **Uncoded source is cleaned.** Two lakes, circles of radius 100 m
   (`quad_segs=64`) at (300, 500) and (700, 500) in a 1,000 m square domain,
   uncoded rules: 514 vertices in. Measured with `simplify_borders`: band
   1 m → 54, 5 m → 28, 20 m → 15, each lake's area 31,412.773 m² kept.
2. **Coded lakes, same fixture with codes:** same counts (cause 2).
3. **Tiling source unchanged:** two halves of the square, rest empty; the
   existing CORINE clean-up and mesh snapshots pass unchanged.
4. **Rules refusals:** missing `positional_accuracy_m`; `label` or `colour`
   in rules; coded without `legend`; a class or sub-class code missing from
   the legend (message names it); `legend` without `code_system`.
5. **Legend refusals:** extra key, repeated code, code 0 in `[[class]]`, bad
   colour; `rasputin palette corine` bytes unchanged; `palette corine-water`
   refused.
6. **Old names refused:** `--features-legend`, `--land-cover-legend`,
   sources key `legend`, MCP `features_legend` (beside the class-map checks
   in `tests/python/test_legend_name.py`).
7. **Record:** two sources give their own line and vertex counts and one
   `features_joined`; a `property`-only run with `--features-tolerance 5`
   prints the not-used line and no clean-up keys.

## 6. Size and checks

About 110 net production lines (`python3 tools/count_loc.py`): rules model
moved (+0), legend model and loader +35, the two clean-up fixes +25, record
+15, CLI names and palette +15, defaults and the not-used line +20. TOML is
not counted; the rename touches about 150 production and 300 test lines
modified. No C++, so no `@perf` acceptance run; the byte-identity of the
shipped coded sources is pinned by test 3 and the existing snapshots.
