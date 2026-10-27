# Increment 108: a legend is the colour scheme; meshing rules are a file of their own

**Status: designed, round 3** (`@architect`, 2026-10-27, on `5de1144d`). One
PR, branch `worktree-legend-split`, about 70 net production lines, no C++
file, no `@perf` run (section 6). The split only: it changes no mesh. The
Sula defects it surfaced (section 3) are increment 109, designed on its own
branch from master once this merges (one increment file per PR). Ola ordered both built and merged today, 2026-10-27.

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
returns twelve files; the land-cover one, `legacy/rasputin/land_cover_repository.py`,
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

## 3. The finding that motivates increment 109

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
   writes the merged lines and vertices into every source's entry
   (Sol's record: "112 lines, 8180 vertices" twice).

All three are fixed in increment 109, not here.

## 4. The design

**Legend file** (`data/legends/<name>.toml`): the colour scheme only.
`name`, `preset` (the ParaView preset name), optional `[no_class]` (label,
colour of code 0), and `[[class]]` with `code`, `label`, `colour`. Every key
required, nothing else allowed; a repeated code or code 0 in `[[class]]` is
refused. A legend may list codes no rules file uses.

**Meshing-rules file** (`data/meshing_rules/<name>.toml`): today's fields less
the colours: `name`, `attribute`, `unlisted`, optional `code_system`,
`positional_accuracy_m` (required with `code_system`, as today; increment 109
makes it required in every file), `notice`, `licence`, `[edges]`,
`[[class]]` (`code` or `value`, `water`, `edges`), `[[split]]`
(sub-classes `value`, `code` only), `codes`, `vocabulary.extra`, `[cells]`
(below), and **`legend`**: a shipped legend name or a path, relative to the
rules file. With `code_system`, `legend` is required and every class and
sub-class code must be in it (refused naming the missing codes); without,
`legend` is refused. `label`, `colour`, `preset` and `[no_class]` are
refused in a rules file, each naming the legend as their place.

**`unlisted = "land"` needs `[edges]`.** Today an uncoded file may state
`unlisted = "land"` with no `[edges]` table, and every unlisted feature is
then dropped without a word: `src_python/tin_engine/legends.py@5de1144d:168-173`
turns `unlisted = "land"` into `"drop"` when `[edges]` is absent, so the file
says land and the run drops. It is refused: "unlisted = \"land\" needs an [edges] table: its land list names
the edges an unlisted feature gets". No shipped file and neither of Sol's
files hits it (all have `[edges]` or are coded).

**`[cells]`: the cell mapping, apart from the edges.** What the cell table
and `--prefer` read is not a meshing rule: a class's `category` and the
`[prefer]` table are read only by the cell table
(`src_python/tin_engine/cells.py@5de1144d:30-37`) and by `--prefer`
(`src_python/tin_engine/cli.py@5de1144d:1985`), never by the building of
constraint lines. So it moves out of `[[class]]` into one table, allowed only with `code_system`:

```toml
[cells.categories]          # each code a cell category counts
lake = [512]
glacier = [335]
forest = [311, 312, 313]
[cells.prefer]              # the one class a --prefer source's polygons take
lake = 512
glacier = 335
```

Keys of `categories` are the `Category` names (`forest`, `lake`, `glacier`,
`savanna`, `unspecified`). Refused: a code in two categories, a code not
listed (class or sub-class), `category` inside `[[class]]` or a top-level
`[prefer]` (each naming `[cells]`). Split sub-classes still take their
parent's category. A class's `water` flag stays in `[[class]]`: it picks the
water edge list; the cells' check that every water class has a category
(`src_python/tin_engine/cells.py@5de1144d:30-37`) is unchanged. The
resulting `categories` and `prefer` mappings equal today's for every shipped
file (test 3).

**Not in 108, in 109:** the model's field names aligned with the file's
(`src_python/tin_engine/feature_input.py@5de1144d:111-118`: model `codes` is
the file's `code_system`, `string_codes` its `codes`, `otherwise` its
`unlisted`, `accuracy_m` its `positional_accuracy_m`), and the request
model's clean-up defaults, all off
(`src_python/tin_engine/feature_input.py@5de1144d:139-145`) where the
command derives 1, on, 50 and 5.

**Python.** `legends.py` keeps only the legend: a frozen `Legend` (`name`,
`preset`, `colours: Mapping[int, tuple[str, str]]`), `load_legend(value)`.
New `meshing_rules.py` takes today's model, loader, `derived_cleanup` and
`run_vocabulary`; the frozen model in `feature_input.py` is renamed
`MeshingRules`, its `palette` and `preset` fields replaced by
`legend: Legend | None` (composed, not copied). `FeatureSource.legend`
becomes `.rules`. Readers of colours (`palette`, the GeoPackage writer's
classes and style) read `rules.legend.colours`. **`known_codes`**
(`src_python/tin_engine/cli.py@5de1144d:1916`, today the palette's keys)
becomes the rules' own codes: every listed class code, every split
sub-class code, and 0 when the legend has `[no_class]`; this equals today's
palette keys for every shipped file, so AR50's sub-classes (3031 ... 5099)
stay known. `property` stays built in, with no legend.

**Names.** No aliases. Each old name is refused with a message naming the
new one:

| Today | After | Refusal |
|---|---|---|
| `mesh --features-legend` | `--features-rules` | "--features-legend is now --features-rules; a legend is the colour scheme only" (one table `RENAMED` in `cli.py`, appended to click's no-such-option error) |
| `mesh --land-cover-legend` | `--land-cover-rules` | as above |
| MCP `features_legend`, `land_cover_legend` | `features_rules`, `land_cover_rules` (derived from the options) | the schema's `additionalProperties: false`; `mcp/server.py`'s instruction text names the new arguments |
| sources-file key `legend` | `rules` | "legend is now rules (a source names its meshing rules)", before the strict model's error |
| `sources` fit field `legend` | `rules` | output only |
| run record `features`: "legend X"; `src_python/tin_engine/run_record.py@5de1144d:35` "layer, legend, ..." | "rules X"; "layer, meshing rules, ..." | output only |
| mesh field `land_cover_codes`: "legend X" | unchanged: X is the rules file's name | output only, kept so gate 6 holds |
| `rasputin palette NAME` (any legend) | a legend name or legend file only | `palette corine-water`: "unknown legend corine-water; use ar50, corine, mapbiomas-c11, or a .toml file" |

The `land_cover_codes` text is written into the mesh
(`src_python/tin_engine/cli.py@5de1144d:2153`, today `legend.name`, the one
file's name), so renaming it would change gate 6's sha256. It takes the rules
file's name: `corine` for both gate cases, as today; `corine-water` and
`clc18_kode` keep their own names although both name the legend `corine`.

Prose renamed in the same PR: `docs/benchmarks/quick/cases.toml` (lines 32
and 46), `INSTALL.md:270`, `README.md` (lines 31, 101, 132, 213, 232, 247 and
the wording at 44, 47, 128, 136), `docs/GLOSSARY.md:246-254` (the legend
entry cut to the colour scheme; a new "meshing rules" entry), and
`project_structure.md:300-317`. The guard test
`tests/python/test_legend_name.py` gains the old option names, key and MCP
arguments.

**Shipped files.** Legends: `corine` (one file: the three CORINE files'
colours are equal, checked with `load` on `5de1144d`), `mapbiomas-c11`,
`ar50` (sub-class colours as plain classes). Rules: `corine`, `corine-water`,
`clc18_kode` (each `legend = "corine"`), `mapbiomas-c11`, `ar50`; every other
field unchanged.

**Ola's sources file at merge.** `~/.config/rasputin/sources.toml` has 15
`legend = "ar50"` keys, all a shipped name that is also a shipped rules
name. The moment 108 merges, the main session copies the file to
`sources.toml.bak-108` and rewrites each `legend =` key at the start of a
line to `rules =`, nothing else; it then checks the file with `rasputin
sources` and tells Ola, who restarts his MCP server. The main session
restarts nothing.

No `--legend` override option in this PR: a user picks colours by naming a
legend file in a rules file.

## 5. Tests and gates

Tests (red first):

1. **Rules refusals:** `label`, `colour`, `preset` or `[no_class]` in rules
   (each naming the legend); coded without `legend`; a class or sub-class
   code missing from the legend (message names it); `legend` without
   `code_system`; an uncoded file with `unlisted = "land"` and no `[edges]`
   (message as above); `category` in `[[class]]`, a top-level `[prefer]`, a
   code in two categories, a `[cells]` code not listed, `[cells]` without
   `code_system`.
2. **Legend refusals:** extra key, repeated code, code 0 in `[[class]]`, bad
   colour; `rasputin palette corine` bytes unchanged; `palette corine-water`
   refused as in the table.
3. **Same colours:** for each of the five shipped names, the composed
   legend's colours and preset, and the rules' `categories`, `prefer` and
   `water_codes`, equal today's `load(name)` fields, written by `@tester` from `5de1144d` into a snapshot
   `tests/python/data/legends_5de1144d.json` in the red commit (the older
   `legends_f09b9863.json` predates `ar50`).
4. **AR50 known codes:** `known_codes` for `ar50` holds 3031, 3032, 3033,
   3099, 5051, 5052, 5054, 5055, 5099 and equals today's palette keys; an AR50
   class-polygon run with a sub-class polygon is not refused.
5. **Old names refused with the new name:** `--features-legend`,
   `--land-cover-legend`, sources key `legend`; MCP `features_legend` absent
   from the schema.

Gates, run by `@developer` on the green head and recorded in its handback:

6. **No mesh changes.** The quick check's `numedalslagen` and `lagan` cases
   (`docs/benchmarks/quick/cases.toml`, args as there with `--features-rules`,
   plus `--out X.vtk --binary`), one run each, sha256 of the `.vtk`. Before,
   measured on master's venv (source identical to `5de1144d`), AC power,
   2026-10-27: numedalslagen
   `1b077a1a5ce8170f7cfebe9a01cdcfe551ecc50d22048f2e73f95a83b32b4705` (two
   runs, same hash), lagan
   `4c9f885a5948444e1c2c9a695faedbf3e163d59a3490d6a4f12f793123d83965`. After
   must equal both.
7. **No time cost.** The same runs' wall time: before 10.41 s and 9.98 s
   (numedalslagen), 16.26 s (lagan). The split adds one small TOML read per
   run, so the expected cost is zero; one run each after is judged within
   noise (10 %), with `@developer` still running beside it.

## 6. Size and checks

About 70 net production lines (`python3 tools/count_loc.py`): rules model
moved (+0), legend model and loader +35, the refusals naming new names +10,
`known_codes` and the colour readers +10, `[cells]` and its refusals +12,
the `unlisted` refusal +3. TOML is not counted; the rename
modifies about 150 production and 300 test lines. Over the short path's 50,
so the full path. No C++ and no mesh change (gate 6), so no `@perf`
acceptance run.

## Review

- Design round 1: CHANGES REQUESTED, origin/master..8196a44f, 0 net production lines (scope beyond the ruling, no speed or byte-identity pin, missed renames, AR50 sub-codes).
- Design round 2: CHANGES REQUESTED, ..87ff7c00, 0 (mesh text land_cover_codes unpinned; 109 on the same branch; audit citations not on master).
- Design round 3: APPROVED WITH WORDINGS, ..24fedc3d, 0; the two status wordings applied in this commit.
- Code round 1: CHANGES REQUESTED, 5de1144d..62e321f9, 106 net production lines (three messages still say legend where meshing rules decide: polygon_fractions.py:90, cli.py:1381, cli.py:1675).
