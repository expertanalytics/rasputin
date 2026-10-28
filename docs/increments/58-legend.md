# Increment 58: "class map" becomes "legend"

**Status: built, awaiting code review** (`@architect` design, 2026-10-16;
design review round 1 approved with wordings; red `a002d924`, green
`f9587095`, prose in the docs commit). Short path
(`docs/increments/README.md@aabac548:76`): about 0 net production lines,
no C++ file, nothing under *Acceptance*. One public PR, branch
`worktree-legend`, cut from increment 57's head `ed0d313d`; it merges after 57.

## 1. What Ola asked

Ola, 2026-10-16: "Class-map is a very bad name for this. It could mean a lot
of things, depending on context." "Keep the guard tight, and let's use
legend." "We don't want to pile up legacy so hard no to aliases."

So: one name, `legend`, everywhere a user, an agent or a reader of the code
meets the thing; the old spellings removed, not kept beside the new ones. No
alias, no deprecation warning, no reader that accepts the old key.

## 2. Prior art: legacy and literature

Written before the design.

*Literature.* "Legend" is the cartographic and land-cover term for this object.
The FAO Land Cover Classification System separates a *classification* (the
abstract scheme of classes) from a *legend* (that classification applied to
one data set at one scale, with its map units): Di Gregorio and Jansen,
*Land Cover Classification System (LCCS): Classification Concepts and User
Manual*, FAO, 2000, chapter "From classification to legend"
(<https://www.fao.org/4/X0596E/X0596e00.htm>); ISO 19144-1:2009 builds on the
same system. Our file fits the legend side: it is tied to one source's
attribute and codes, carries labels and colours for display, and states the
source's positional accuracy. The distributors use the word too: the EEA
calls CORINE's code table the CLC legend, and MapBiomas ships
`legend_code_mapbiomas_brazil_collection_11.csv`, which our
`mapbiomas-c11.toml` cites. One departure: our legend also says which edge
names (water, land, extra names) each class's borders carry in the mesh; a
cartographic legend has no such part. That is ours, and is named as such in
the glossary entry (section 6). No novelty is claimed.

*Legacy.* Nothing: the legacy tree had no class tables or legends. The grep,
at the tag, returned no file (exit status 1):

```
$ git grep -l -i -E 'class[_ -]?maps?|legend' legacy-archive -- legacy
$ echo $?
1
```

*Ideas file.* No entry asks for this rename. Overlapping entries, not taken
up: "Product review: readable code docs, and general applicability" (its
generality audit made the class maps data files, increment 41; this renames
what 41 built); "A check that project_structure.md lists every package module"
(this PR renames one module and one data folder in that file by hand).

## 3. The rename

Counted at `ed0d313d` with
`git grep -n -i -P 'class[_ -]?maps?|features[-_]map\b|land[-_]cover[-_]map\b|\bcmap\b'`:
390 lines in 48 files outside the history records (section 5): 125 in
`src_python/`, 250 in `tests/`, 6 in `README.md`, 1 in `INSTALL.md`, 5 in
`project_structure.md`, 3 in `docs/research/`. The history records hold 1 477
more and do not change. The brief's narrower pattern
(`class[_ -]?maps?`) finds 997 lines in the whole tree.

### 3.1 Names users and agents see

| Old | New |
|---|---|
| `rasputin mesh --features-map` | `--features-legend` |
| `rasputin mesh --land-cover-map` | `--land-cover-legend` |
| MCP `mesh` arguments `features_map`, `land_cover_map` | `features_legend`, `land_cover_legend` (derived from the options, `mcp/schema.py`; no hand edit) |
| sources file key `class_map = "..."` | `legend = "..."` |
| `sources` listing (command and MCP tool): each entry's `class_map`, each fit's `class_map` | `legend` |
| `rasputin palette` argument help "A class map: ..." | "A legend: ..." |
| run record, `features` entry wording "Each features file: layer, class map, features, lines, vertices" and value "`x.gpkg`, class map corine: ..." | "layer, legend, ..." and "legend corine" |
| errors "class map file X cannot be read", "unknown class map X", "class map X: ...", "... is not in the map X" | "legend file X cannot be read", "unknown legend X", "legend X: ...", "... is not in the legend X" |
| `--features-map` help "The Nth --features' map from attributes to edge bits: ..., or a class-map file (.toml)" | "The Nth --features' legend: ..., or a legend file (.toml)" |
| shipped legends folder `tin_engine/data/class_maps/` | `tin_engine/data/legends/` (file names unchanged: `corine.toml`, `corine-water.toml`, `clc18_kode.toml`, `mapbiomas-c11.toml`) |

The legend file's own keys (`name`, `attribute`, `code_system`, `edges`, ...)
contain no "map" and do not change, so a user's own legend `.toml` passed by
path keeps working. The names of the legends (`corine`, `property`, ...) do not
change.

### 3.2 Internal names (no half-renamed code)

- Module `tin_engine/class_maps.py` → `tin_engine/legends.py` (`git mv`, so
  `count_loc.py`'s `-M` sees a rename). `load` keeps its name; `cli.py`'s
  import alias `load_class_map` → `load_legend`.
- `feature_input.ClassMap` → `Legend`; `FeatureSource.class_map` → `legend`.
- `sources.RemoteSource.class_map` → `legend`, and `_read`'s relative-path
  handling of that key.
- `source_report.Fit.class_map`, `Facts.class_map` → `legend`; `fits(facts,
  maps)` and every `maps: Mapping[str, ClassMap]` parameter → `legends`.
- `cli.py`: parameters `features_map`, `land_cover_map` → `features_legend`,
  `land_cover_legend`; locals `feature_maps`, `maps`, `map_value`, `coded_map`,
  `cmap` → `feature_legends`, `legends`, `legend_value`, `coded`, `legend`;
  `_ClassSource.class_map` → `legend` (`_ClassSource` keeps its name: it is a
  source of classes, not a map); the `--features` help sentence that lists
  "--features-crs, --features-layer and --features-map".
- `feature_input.py`, `source_report.py`: local `cmap` → `legend`.
- `mcp/server.py` `INSTRUCTIONS`: "the class map to pass as features_map" →
  "the legend to pass as features_legend" (twice).
- Docstrings and comments in `palettes.py`, `legends.py`, `feature_input.py`,
  `sources.py`, `source_report.py`, `cli.py` say "legend", including a bare
  "map" that means the legend (`cli.py` lines 252, 784, 902, 985, 1760;
  `feature_input.py` lines 71, 564).
- The three CORINE `.toml` headers point at `docs/increments/41-...`; the
  pointer stays (it is a history record's file name, section 5).
- `data/sources.toml`, MapBiomas entry: `class_map` → `legend`.
- `plugin_api.py` exports none of these names: no change. `docs/plugins.md`
  names no class map: no change beyond what section 4 says.

### 3.3 What happens to old inputs (no fallback reader)

Each case was checked against the code at `ed0d313d`:

- **Old option on the command line** (`--features-map corine`): Typer refuses
  it, "No such option", exit 2. Click's suggestion comes from
  `difflib.get_close_matches` over the option names; computed against the
  renamed options it offers `--features`, `--features-layer`,
  `--features-crs` for `--features-map` (not `--features-legend`), and
  `--land-cover`, `--land-cover-year`, `--land-cover-legend` for
  `--land-cover-map`. Accepted: no hint code is added (it would be the
  legacy Ola ruled out); README and the PR description carry the new names.
- **Old key in a user's sources file** (`class_map = "corine"`):
  `RemoteSource` is `extra="forbid"` (`sources.py`, line 80), so the whole
  file is refused, as today for any unknown key, naming the file and the key.
- **Old key in a plugin's `sources.toml`**: that plugin's file is skipped as
  one `refused` line (`catalogue`, unchanged); the core runs. No plugin in
  Ola's repositories uses the key (`git grep` in `rasputin-radiation`, the
  one plugin checked out beside this repository: no hit).
- **Old MCP argument** (`features_map`): not in the tool's schema, whose
  argument model is `extra="forbid"` (`mcp/schema.py`, line 153); refused.
- **Mesh files written before**: the run record's field *names* do not change
  (the entry is `features`); only the value text does ("class map corine" →
  "legend corine"). Nothing in `tin_engine` reads that field back:
  `io/vtk_read.read_mesh` reads points and cells, not FieldData. An old
  `.vtk`, `.ply` or `.gpkg` reads exactly as before; its text keeps the old
  words.
- **A `--record` JSON written before**: same keys, old value text.

## 4. Prose that changes

- `README.md` (6 lines: the CORINE example command, the land-cover paragraph,
  the `sources` paragraph, the sources-file example's `class_map` key and its
  sentence).
- `INSTALL.md` line 270 (`--features-map corine`).
- `project_structure.md`: the `class_maps.py` and `data/class_maps/` entries,
  `palettes.py`'s pointer, `feature_input.py`'s "by a ClassMap".
- `docs/GLOSSARY.md`: a new section "Land cover" with one entry, **Legend**:
  a TOML file that names one land-cover source's attribute, its classes'
  codes, labels and colours, its positional accuracy, its notice and licence,
  and, which a cartographic legend lacks, the edge names each class's borders
  carry; shipped ones in `tin_engine/data/legends/`, chosen with
  `--features-legend` and `--land-cover-legend` or a sources-file `legend`
  key.
- `docs/research/raster-to-vector.md` line 597 ("a class map may coarsen it"
  → "a legend may coarsen it"); `docs/research/nmd-landcover.md` heading
  "Class mapping" → "Matching NMD's classes to CORINE's" (it means a
  correspondence between two class systems, not our file).
- `docs/research/lagan-slivers.md` line 71 (`--features-map corine` →
  `--features-legend corine`).

`@architect` writes these after green, in a docs commit, as in increment 56.

## 5. What does not change: the history records

`docs/increments/`, `docs/benchmarks/`, `docs/retrospectives/` and the
ROADMAP rows of earlier increments record what was built, with the names of
the time; rewriting them would make them false. They keep "class map" and
`--features-map`. Increment 41's file keeps its name,
`41-class-maps-as-data.md`, and the pointers to it stay. Benchmark stats files
are not compared as text by `tools/bench.py` (no hit for the pattern in
`tools/`), so their old wording breaks nothing.

## 6. Tests

### 6.1 The guard (new): `tests/python/test_legend_name.py`

One test. It lists the tracked files (`git ls-files`, run from the repository
root the test file sits under) and fails naming file and line for every line
matching, case-insensitively,
`class[ _-]?maps?|features[-_]map\b|land[-_]cover[-_]map\b|\bcmap\b`.
Tight: exactly three exceptions, each a named rule, no allow-list of lines:

1. paths under `docs/increments/`, `docs/benchmarks/`,
   `docs/retrospectives/`, and `ROADMAP.md` (history records, section 5);
2. the token `41-class-maps-as-data` is removed from a line before matching
   (a history record's file name);
3. the guard's own file.

It fails, not skips, when `git ls-files` lists no file. Fixture figures at
`ed0d313d`: 390 matching lines in 48 files outside rule 1; the token of rule 2
on 9 lines (3 shipped `.toml` headers, `class_maps.py`, 5 test docstrings).
Its red state is those 390 lines; after green it must report 0. Made able to
fail: `@tester` plants one `class_map` in a scratch line of a counted file,
sees the guard name it, removes it.

### 6.2 Existing tests that change (red step)

`@tester` renames, in 30 files under `tests/python/` (250 matching lines),
every option, key, field, import and message the tests use to the new names,
so they fail against `ed0d313d` and pass after green:

- files renamed: `test_class_maps.py` → `test_legends.py`,
  `test_cli_class_maps.py` → `test_cli_legends.py`, `class_map_fixtures.py` →
  `legend_fixtures.py`, `class_map_snapshot.py` → `legend_snapshot.py`,
  `data/class_maps_f09b9863.json` → `data/legends_f09b9863.json` (its
  provenance string names `feature_input`'s old table by a description, not by
  its old identifier);
- `data/catalogue_b1353075.json`: the three `"class_map"` keys spelled
  `"legend"`, as the model dumps them now; `test_catalogue_data.py`'s
  docstring says the snapshot's keys are spelled as the model spells them;
- the rest (options and fields only): `feature_fixtures.py`,
  `fraction_fixtures.py`, `vector_fixtures.py`, `test_catalogue_data.py`,
  `test_cli_features_cleanup.py`, `test_cli_land_cover_file.py`,
  `test_cli_mesh_features.py`, `test_cli_mesh_gpkg.py`,
  `test_cli_mesh_land_cover_fractions.py`, `test_cli_mesh_landcover.py`,
  `test_cli_mesh_multi_features.py`, `test_cli_mesh_plain_output.py`,
  `test_cli_sources_vector.py`, `test_feature_chains_bits.py`,
  `test_feature_input.py`, `test_grow.py`, `test_io_geopackage_write.py`,
  `test_io_gml.py`, `test_layering.py`, `test_mcp_schema.py`,
  `test_mcp_vector_sources.py`, `test_palettes.py`, `test_source_fits.py`,
  `test_sources_vector.py`.

New pins, one each, added to the existing files:

- `--features-map corine` and `--land-cover-map mapbiomas-c11` are refused
  with exit 2 and "No such option" (no alias);
- a user sources file with `class_map = "corine"` is refused naming the file
  and `class_map` (no fallback reader);
- the MCP `mesh` tool schema has `features_legend` and `land_cover_legend` and
  neither old name.

No mutation round: a rename's failure mode is a missed spelling, which the
guard catches.

## 7. Size

A rename is modified lines: each counts once added and once removed, so the
net is near 0. Measured at `ed0d313d` with `count_loc.py`'s own line
classifier (`counted_lines`) over `src_python/` and `tools/`: 113 counted
lines carry an old name (`cli.py` 50, `source_report.py` 25,
`feature_input.py` 19, `class_maps.py` 11, `mcp/server.py` 4, `sources.py` 3,
`run_record.py` 1), plus 28 docstring and comment lines that do not count.
Expected `count_loc.py` result: about +113 / -113, net 0 to +10 (a longer
name may push a line past ruff's width and split it). The `.toml` data files
and the tests do not count.

## 8. What a user must change (for the PR description)

- Command lines and scripts: `--features-map` → `--features-legend`,
  `--land-cover-map` → `--land-cover-legend`. The old options are refused.
- Your sources file (`~/.config/rasputin/sources.toml` or
  `$RASPUTIN_DATA/sources.toml`): the key `class_map` → `legend`. A file with
  the old key is refused whole until it is changed.
- MCP clients: the `mesh` arguments `features_map` / `land_cover_map` →
  `features_legend` / `land_cover_legend`; the `sources` listing's
  `class_map` fields are now `legend`.
- Your own legend `.toml` files: no change. Meshes written before: no change.

## 9. Defaults taken (Ola away)

1. "The guard" is read as a test that fails on any old spelling outside the
   history records, with three named exceptions (section 6.1).
2. History records (`docs/increments/`, `docs/benchmarks/`,
   `docs/retrospectives/`, past ROADMAP rows) keep the old words.
3. No hint for a mistyped old option; Click's suggestion stands even though it
   does not name the new option.
4. `_ClassSource`, `load`, the legend names and the legend file keys keep
   their names.
5. Short path: one `@reviewer` round over design, red and green.

## Review

- **Design review, round 1, 2026-10-16.** Commit `60700288`, 0 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: `lagan-slivers.md` line 71; the `--features-map` help row; bare "map" meaning the legend in `cli.py` and `feature_input.py`. Held: 390 lines in 48 files; strict models refuse old keys and arguments; old mesh files read; no alias or hint survives. Not blocking: about 15 more counted lines (locals `maps`, `feature_maps`, `map_value`, `coded_map`); the guard could also match tracked file paths.
