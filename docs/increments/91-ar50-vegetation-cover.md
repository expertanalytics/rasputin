# Increment 91: a legend splits classes by a second attribute (AR50's open firm ground by vegetation cover, forest by tree type)

**Status: designed** by `@architect`, 2026-10-22, on `2600a740` (master),
branch `worktree-ar50-vegetation`; amended after `@reviewer`'s design review
round 1 (rows stay 3-tuples, one rule for a missing second attribute,
whole-number matching, the quick check's cost, no listing change) and after
Ola's ruling to split forest too (several splits per legend). Full path:
about 70 counted lines (section 7), over the short path's 50; no C++ file, no
refine or mesh code, so no `@perf` acceptance run. Public.

Approved by Ola today: open firm ground (AR50 `arealtype` 50) split by
vegetation cover (`vegetasjonsdekke`), and, by his later ruling, forest
(`arealtype` 30) split by tree type (`treslag`) as well. His reason for the
forest split: AR50's coniferous areas can later bound where the spruce of
NIBIO's forest resource map (SR16) is taken, and the two approaches (AR50's
tree type alone; SR16 within AR50's coniferous forest) are to be compared.
That comparison is a later step of the forest plan, not this increment. Codex
found that the shipped `ar50` legend
(`src_python/tin_engine/data/legends/ar50.toml`) reads only `arealtype`, so a
mesh loses AR50's vegetation cover and tree type.

## 1. Prior art: legacy and literature

*Literature and standards.* No method is claimed and no claim of novelty is
made: a class refined by a second attribute is a table lookup. The published
precedent is NIBIO's own: its presentation rules for AR50 (*Presentasjonsregler:
Arealressurskart - AR50*, 12 June 2026,
`https://dokument.geonorge.no/tegneregler/ar50/20260615/presentasjonsregler-ar50.pdf`,
downloaded 2026-10-22 and read with `pdftotext`) draw the map layer
"Bonitet, 18 klasser" with queries of exactly this form, `"ARTYPE" = 50 AND
"ARVEGET" = 55` and `"ARTYPE" = 30 AND "ARTRESLAG" IN (31,32,33,99)`. Its
table "Snaumark (ARVEGET)" colours the vegetation-cover codes and its table
"Treslag (ARTRESLAG)" the tree-type codes 31, 32 and 33; each states that its
attribute covers the same ground as its arealtype (50, 30). The codes and
their meaning: the schema `ar50.xsd` (SOSI AR50 2022-05-23,
`ArealressursVegetasjonsdekkeEnumerationType`,
`ArealressursTreslagEnumerationType`) and NIBIO's product specification,
section 5.1.15 `ArealressursVegetasjonsdekke` (both cited in
`docs/increments/73-gml-vector-source.md`, section 1). Departure from NIBIO:
its 18-class layer lumps 54 with 99 and 52 with 53, and does not split forest
by tree type; this legend keeps each code apart (section 3), and a model that
wants NIBIO's lumping maps the codes itself.

*Legacy.* Nothing. The greps, and the file lists they returned (none; exit
status 1; the same command with "triangle" lists `legacy/bindings.cpp`, so
the ref resolves and the probe can match):

```
$ git grep -n -i -l "vegetasjon\|arealtype\|ar50\|snaumark\|treslag\|barskog" legacy-archive -- legacy
$ git grep -n -i -l "sub.class\|subclass\|second attribute" legacy-archive -- legacy
```

*Earlier legend increments.* 41 (`docs/increments/41-class-maps-as-data.md`,
section 4) set the file format and its loader; 58 renamed it "legend"; 73
(section 4.5) shipped `ar50` over `arealtype` alone and listed the other
coded columns as left out (section 7 there). This increment extends 41's
format by one optional list of tables; nothing in 41's mapping changes.

*Ideas.* Overlaps, not taken up: the private package's "Evaporation per
triangle from land-cover classes" (a consumer of these codes) and "spruce
apart from pine" (SR16; this increment gives it AR50's coniferous areas to
bound it, and the comparison Ola asked for belongs there). Nothing in the core
ideas file is taken up.

## 2. What the data holds (measured)

`sqlite3` over Innlandet's GeoPackage
(`$RASPUTIN_DATA/AR50_UTM32_Innlandet_20261019/ar50.gpkg`, 287,730 features;
the three columns `INTEGER`, every value `typeof` integer):

| `arealtype` | `vegetasjonsdekke` (count) | `treslag` (count) |
|---|---|---|
| 30 | 98 only | 31 (90,706), 32 (12,067), 33 (17,440), 99 (14,241) |
| 50 | 51 (3,421), 52 (13,267), 54 (17,519), 55 (8,370), 99 (2,095) | other values |
| every other | 98 only | other values |

(`treslag` on land that is not forest holds 31, 32, 33, 39 or 98: tree cover
of marsh and open ground. The forest split reads it only for class 30.)

The same queries over the 14 other county files Ola downloaded today
(`AR50_Norway_20261022/ar50_{11,15,18,39,40,42,46,50,55,56}.gpkg`,
`AR50_UTM32_OsloViken_20261022/ar50_{03,31,32,33}.gpkg`): arealtype 50 always
has one of 51, 52, 54, 55, 99, and no other arealtype has anything but 98;
arealtype 30 always has one of 31, 32, 33, 99. **53 (lichen cover) occurs in
none of the 15 files**, nor does forest with `treslag` 39 or 98, though the
schema lists them.

## 3. The legend file: an optional list of `[[split]]` tables

A coded legend may carry any number of `[[split]]` tables, each with the
class it splits, its second attribute, and one `[[split.class]]` per value of
that attribute.

```toml
[[split]]
code = 30                        # the class split; a listed [[class]] code
attribute = "treslag"            # its second attribute, read only for that class

[[split.class]]
value = 31                       # the second attribute's value, matched by its key
code = 3031
label = "Forest, coniferous: at least 50 % of the forested area covered by conifers"
colour = "#895944"

[[split]]
code = 50
attribute = "vegetasjonsdekke"

[[split.class]]
value = 51
code = 5051
label = "Open firm ground, no vegetation"
colour = "#c0c0c0"
```

Rules, each refused by the loader with a `ValueError` naming the file (as
every legend refusal is):

1. `[[split]]` only in a coded legend (one with `code_system`).
2. Each split's `code` is a listed `[[class]]` code, and no two splits name
   the same class.
3. Each split's `attribute` differs from `attribute` (case-insensitively, as
   the GeoPackage reader matches columns). Two splits may share an attribute.
4. Every sub-class has `value`, `code`, `label`, `colour`, nothing else
   (`extra="forbid"`): no `water`, no `edges`.
5. Within one split the values are distinct by key (below); one value may
   recur in another split (99 in both AR50 splits). The sub-class codes of
   all splits are distinct and none is a listed `[[class]]` code. This joins
   the existing "a class code is listed twice" check over the union of all
   lists.

The loader keeps the splits sorted by class code, whatever their order in
the file.

**The key of a value.** The same function keys a sub-class's `value` in the
loader and a feature's second value at run time, so 52, 52.0, `"52"` and
`"52.0"` all match a sub-class `value = 52`:

- an integer (not a bool): its decimal text, `"52"`;
- a finite float with no fractional part: the integer's text, so a GeoJSON
  52.0 gives `"52"`;
- a text that `float()` reads as a finite whole number: the same, after
  stripping spaces (`" 52 "`, `"52.0"`, a GML element's text);
- any other text: the text itself (a code system with letters, `"5a"`);
- None, a bool, a non-whole number (52.5), NaN, a list: no key, so no
  sub-class matches and the parent's code stays.

So a loader that sees `value = 52` and `value = "52.0"` in one split refuses
them as one value listed twice (rule 5).

**What a sub-class inherits.** Its borders carry the parent's edge names
(`_mask` reads the first attribute only, so it does not change); it is never
water unless the parent is. Its code, label and colour are its own, and all
three enter the palette beside the parent's.

**How a feature is coded.** The first attribute gives the code as today
(`_code`). If a split names that code, the key of that split's attribute's
value is looked up among its sub-classes: found, the sub-class code; not
found (a key the split does not list, a value with no key, a null, a missing
GeoJSON property or GML element, a GeoPackage without the column), **the
parent's code stays**. A feature reads only its own class's split attribute:
forest ignores `vegetasjonsdekke`, open firm ground ignores `treslag`, and a
class no split names ignores both. A source with no second attribute at all
is section 6's case. A first-attribute value whose code equals a sub-class
code (an unlisted value under `unlisted = "land"`, say a stray 5051) is
refused naming the file, the feature and the value, so a sub-class code
never means two things in one mesh.

**Numbering, so codes never collide.** The loader does not impose a scheme:
the file states each code, and rule 5 refuses a collision with any listed
code, and the run-time refusal above one with any unlisted code. The AR50
file numbers each sub-class `100 × arealtype + value`: 3031 to 3099 for
forest, 5051 to 5099 for open firm ground. AR50's arealtype codes are two
digits (10 to 99, `ar50.xsd`), so a four-digit code cannot be one; and the
number reads back as its two parts.

**The 99 cases** ("not registered": 2,095 open-firm-ground and 14,241 forest
polygons in Innlandet): each its own sub-class, 5099 and 3099, not folded
into the parent. NIBIO's "Snaumark (ARVEGET)" table gives vegetation cover 99
a class and colour of its own, and its 18-class layer lists forest 99 beside
31, 32 and 33, so this copies the source rather than asserting a lumping; a
model may still map 5099 as 50 and 3099 as 30. The parent codes 30 and 50
stay listed for polygons whose second value the split does not list (98, a
value outside the schema's codes, or none): in all 15 county files, no such
polygon exists.

**Left out of the lists**: 53 (lichen cover) and forest's 39 and 98, in no
polygon of any of the 15 files (section 2). If a later release has one, it
keeps the parent's code by the rule above; adding it is one table in the data
file.

### The shipped `ar50.toml`

Header and the nine `[[class]]` tables unchanged (so `code_system` stays
"AR50 arealtype code"), then the forest split:

```toml
[[split]]
code = 30
attribute = "treslag"
```

| `value` | `code` | `label` | `colour` | NIBIO's RGB ("Treslag (ARTRESLAG)") |
|---|---|---|---|---|
| 31 | 3031 | Forest, coniferous: at least 50 % of the forested area covered by conifers | `#895944` | 137-89-68 |
| 32 | 3032 | Forest, deciduous: less than 20 % of the forested area covered by conifers | `#8b786f` | 139-120-111 |
| 33 | 3033 | Forest, mixed: 20 to 50 % of the forested area covered by conifers | `#8b786e` | 139-120-110 |
| 99 | 3099 | Forest, tree type not registered | `#9ecc73` | none; forest's own colour |

then the open-firm-ground split:

```toml
[[split]]
code = 50
attribute = "vegetasjonsdekke"
```

| `value` | `code` | `label` | `colour` | NIBIO's RGB ("Snaumark (ARVEGET)") |
|---|---|---|---|---|
| 51 | 5051 | Open firm ground, no vegetation | `#c0c0c0` | 192-192-192 |
| 52 | 5052 | Open firm ground, barren to non-continuous vegetation cover | `#ffe8c2` | 255-232-194 |
| 54 | 5054 | Open firm ground, barren to good vegetation cover | `#e1c790` | 225-199-144 |
| 55 | 5055 | Open firm ground, good to lush non-continuous vegetation cover | `#c7a971` | 199-169-113 |
| 99 | 5099 | Open firm ground, vegetation cover not registered | `#f8f7ba` | 248-247-186 |

- **Labels**: the parent's label ("Forest", "Open firm ground") and the
  English definition of each code in `ar50.xsd` (and, for vegetation cover,
  specification section 5.1.15), ASCII, so a label read alone in the
  GeoPackage says what it is.
- **Colours**: NIBIO's, converted from RGB. NIBIO gives 32 and 33 colours one
  step apart in blue (111, 110), so the two read as one in a map; they are
  kept as NIBIO gives them. Forest 99 has no NIBIO colour and takes forest's
  own, `#9ecc73`. NIBIO's vegetation table labels 54 "Ikke registrert
  vegetasjonsdekke", the same words as 99; its 18-class table calls the 54
  colour "middels frisk vegetasjon", so the row's label is NIBIO's slip and
  its colour is taken as 54's.
- The file's header comment gains one sentence naming the sources of the
  splits' codes and colours, and this file, section 3.

## 4. Blueprint: data flow and boundaries

```
ar50.toml --tomllib--> legends._File (+ _Split, _SubClass)   validation, rules 1-5
                              | .legend()
                              v
feature_input.Legend (+ splits: tuple[Split, ...])  palette += sub-classes
                              |
read_source(path, layer, attribute, box_for, crs, extra=<the splits' attributes>)
   .gpkg  query_features(..., extra)  one more column per name in the SELECT,
                                      NULL for a name the layer has no column for
   .gml   read_gml(stream, attribute, extra)   None for an absent element
   .json  p.get(name) per name                 None for an absent property
   rows:  (fid, geometry, value)    unchanged, 3-tuples
   extra: ((v1, v2, ...), ...)      one tuple per row, in `extra`'s order;
                                    () when no extra attribute is asked
                              |
_Tally.source: pairs each row with its extra tuple when the legend splits
_Tally._take --> _mask(value)            unchanged: edges from the first attribute
             --> _code(legend, value, what, extra by name)   parent -> sub-class code
                              |
clean-up (merge by code: sub-classes stay apart) -> label_triangles -> codes
                              |
writers: .vtk/.ply land_cover_code; .gpkg triangles, land_cover,
land_cover_classes, QML styles, all from the palette (unchanged code)
```

Types (Python, Pydantic V2, frozen):

```python
class Split(BaseModel):            # feature_input.py, beside Legend
    model_config = ConfigDict(frozen=True, extra="forbid")
    code: int                      # the class split
    attribute: str                 # its second attribute
    codes: Mapping[str, int]       # split_key(value) -> sub-class code

class Legend(BaseModel):
    ...
    splits: tuple[Split, ...] = () # sorted by code; no two with one code
```

- **`legends.py`**: `_SubClass(_Strict)` (`value: int | str`, `code: Code`,
  `label`, `colour` with the existing pattern), `_Split(_Strict)` (`code: Code`,
  `attribute: str`, `classes: list[_SubClass] = Field(alias="class")`),
  `_File.split: list[_Split] = []`; rules 1 to 5 in `_consistent`;
  `legend()` adds each sub-class to the palette and builds the sorted
  `splits`.
- **`feature_input.py`**: `Split`, `Legend.splits`, and `split_key(value) ->
  str | None` (section 3's key, used by the loader too); `read_source` takes
  `extra: tuple[str, ...] = ()`; `SourceRows` gains `extra: tuple[tuple[Any,
  ...], ...] = ()`, one tuple per row in row order, filled only when `extra`
  is given. **`rows` stays a list of 3-tuples**, so no caller that does not
  ask for extra attributes changes. `_Tally.source` passes the splits'
  attributes, each once (a name repeated case-insensitively is read once),
  and hands `_take` each row's extra values by name (none when the legend
  has no split); `_code` takes them and applies section 3's rule and refusal.
- **Every caller of `read_source`** (`git grep -n "read_source(" --
  src_python`), and what it does after this change:
  - `feature_input.py` `_Tally.source`: asks for the splits' attributes when
    the legend has splits; the only reader of `extra`.
  - `feature_input.py` `read_lake_polygons`: unchanged (no extra; unpacks
    3-tuples).
  - `tolerance_field.py` `shapes_in` (lines 67-68): unchanged, same reason.
  - `io/network_file.py` `units_at` (lines 120-121): unchanged, same reason.

  Their existing tests stay green unchanged, which is the check that the rows
  kept their shape: `tests/python/test_sources_arcgis.py`,
  `tests/python/test_cli_tolerance_near.py`,
  `tests/python/test_cli_network_catchment.py`, and the direct
  `read_source` tests in `test_feature_input.py`, `test_gml_import.py` and
  `test_border_simplify.py`.
- **`io/geopackage.py`**: `query_features(..., extra: tuple[str, ...] = ())`;
  `RawFeature.extra: tuple[Any, ...] = ()`. A name the layer has no column
  for is selected as `NULL` (section 6), matched case-insensitively as the
  first attribute is. The first attribute's missing column is still refused,
  as today.
- **`io/gml.py`**: `read_gml(stream, attribute, extra=())`;
  `GmlFeature.extra: tuple[str | None, ...] = ()`. The GeoJSON branch reads
  `tuple(p.get(name) for name in extra)`.
- **`cli.py` `_land_cover`**: the `land_cover_codes` text gains, per split in
  code order, `, class 30 split by treslag into 3031, 3032, 3033, 3099` (the
  split's code, attribute and sorted sub-codes), between the attribute and
  the legend name. The whole text stays ASCII (`checked_ascii` already guards
  it).
- **`source_report.py`**: unchanged. `rasputin sources` (and the MCP
  `sources` tool) does not list the split attributes' values in this
  increment (Ola's default, taken at design review).

Boundaries: no C++; nothing crosses into the core but the polygons and their
integer codes, as today; the legend stays data (TOML) checked once, the reader
stays generic (it reads named columns, not "AR50"). Each part is testable
alone: the loader from a TOML string, `split_key` and `_code` from values and
a `Legend`, the readers from a small file.

## 5. What a mesh carries

- **`.vtk` / `.ply`**: `land_cover_code` per triangle, now 3031 to 3099 on
  forest and 5051 to 5099 on open firm ground where the second value is
  registered, 30 and 50 where it is not given; the dataset string
  `land_cover_codes` with one split clause per split (section 4).
- **`.gpkg`**: the `triangles` layer's `code`, `land_cover` (label) and
  `colour`; the `land_cover` polygon layer per code; `land_cover_classes` with
  triangles and area per code; the QML style's categories. All come from the
  palette, so no writer changes.
- **`rasputin palette ar50`**: 19 entries (no-class, nine classes, four
  forest and five open-firm-ground sub-classes).
- **The run record**: `land_cover_codes` as above; notice and licence
  unchanged.

**Cost (measured, scratch probe, AC power, this branch's base).** Merging by
class keeps sub-classes apart, so their shared borders become constraint
lines. Each probe meshes the current `ar50` against copies of the GeoPackages
with arealtype 50 recoded `5000 + vegetasjonsdekke` ("vegetation"), and also
arealtype 30 recoded `3000 + treslag` ("both"), under legends listing those
codes: what this design computes. `rasputin mesh` on DTM10
(`DTM10_UTM33_20260925`); wall time is the median of three runs after one
warm-up.

| case | | `ar50` today | vegetation | both |
|---|---|---|---|---|
| Jotunheimen, 15 km square (EPSG:25832 x 457,500-472,500, y 6,812,500-6,827,500; 62 % open firm ground, 3 of 248 polygons forest), `--tolerance 5` | land-cover vertices after clean-up | 10,638 | 25,550 | 25,550 |
| | output triangles | 143,549 | 156,671 (+9 %) | 156,671 (+9 %) |
| | wall time | 1.31 s | 1.62 s | 1.59 s |
| Elverum, 15 km square (x 630,000-645,000, y 6,745,000-6,760,000; 920 of 1,467 polygons forest), `--tolerance 5` | land-cover vertices | 31,807 | 31,845 | 34,906 (+10 %) |
| | output triangles | 60,678 | 60,728 | 64,235 (+6 %) |
| | wall time | 2.01 s | 2.01 s | 2.08 s |
| Numedalslågen, the quick check's case with its flags (`--tolerance 10`), AR50 files 32, 33, 39, 40 in place of CORINE | land-cover vertices | 785,120 | 1,182,054 | 1,316,609 |
| | output triangles | 1,587,559 | 2,041,638 (+29 %) | 2,178,870 (+37 %) |
| | wall time | 80.1 s | 91.7 s (+15 %) | 91.8 s (+15 %) |

The Numedalslågen probe needed repaired copies: `ar50` over the four county
files as downloaded stops today with a GEOS `TopologyException`, both in the
files' own CRS and after the move to the DEM's (1,229 polygons are invalid in
EPSG:25832, and 54 more become invalid when moved to EPSG:25833). The
copies were made valid with `shapely.make_valid` and stored in EPSG:25833, so
no run reprojects; the three columns compare like with like.

**The quick check's own cases** (`docs/benchmarks/quick/cases.toml`) use the
`corine` legend or none; no shipped legend but `ar50` has a split. With no
split, `read_source` gets `extra=()`: the GeoPackage `SELECT` is the same
statement, GML and GeoJSON read nothing more, and `_Tally` pairs each row
with an empty tuple. The expected cost is none; a `tools/bench_quick.py` run
on the code is the check.

## 6. A source without the second attribute, and a mesh made with the old legend

**One rule for all three readers: a missing second attribute is no value.**
A GeoPackage layer without the column (matched case-insensitively), a GML
feature without the element, a GeoJSON feature without the property: each
gives None as that feature's value for it, so the feature keeps the parent's
code, 30 or 50 (section 3). Nothing is refused for it. This keeps every AR50
export that meshes under `ar50` today meshing, with the codes it gets today
(an export carrying only `arealtype`, say). The first attribute's rules do
not change: its missing GeoPackage column is still refused.

What tells Ola whether a split took effect is the mesh, not an error: its
codes (3031 to 3099, 5051 to 5099 present or not) and the
`land_cover_classes` table. A misspelt `attribute` in a hand-written legend's
`[[split]]` is the cost of this rule: it splits nothing, silently. The
shipped `ar50` names the columns NIBIO's schema names, so this applies only
to a hand-written legend.

**A mesh made with the old legend.** Unchanged on disk and still read by
`read_mesh`; nothing migrates and no alias or version key is added. Its
`land_cover_codes` text has no split clause, which is how a reader tells the
two apart. Its codes stay valid under the new legend: 30 and 50 are still
listed with the same labels and colours, so `rasputin palette ar50` still
colours it. A mesh with the splits needs a re-run. The private package has no
code that maps AR50 codes (`git grep -n -i "ar50\|arealtype" private/master`
finds only the ideas file, research notes and one research script reading
the GeoPackage directly).

## 7. Size, and where the prose changes

| File | Counted lines |
|---|---|
| `legends.py`: `_SubClass`, `_Split`, the list field, rules 1 to 5, palette and `splits` in `legend()` | about 30 |
| `feature_input.py`: `Split`, `Legend.splits`, `split_key`, `read_source`'s `extra` and `SourceRows.extra`, `_Tally.source`'s pairing, `_code`'s lookup and refusal | about 28 |
| `io/geopackage.py`: the extra columns, `NULL` when absent | about 6 |
| `io/gml.py`: the extra elements | about 3 |
| `cli.py` `_land_cover`: the split clauses | about 3 |
| `tolerance_field.py`, `io/network_file.py` | 0 (rows stay 3-tuples) |
| `data/legends/ar50.toml` | 0 (TOML is not counted) |
| Total | **about 70**, under the 700 of `CLAUDE.md` section 2 |

Prose, in the same PR (not counted):

- **README.md**, the "Features as constraints" item: the `ar50` legend reads
  its land type (`arealtype`) and splits forest by tree type (`treslag`, codes
  3031 to 3099) and open firm ground by vegetation cover (`vegetasjonsdekke`,
  codes 5051 to 5099).
- **docs/GLOSSARY.md**, *Legend*: one sentence: a legend may split classes by
  a second attribute into sub-classes with codes of their own.
- **project_structure.md**: the `data/legends/` row says `ar50` splits 30
  and 50; the `legends.py` row names the splits.
- **MCP tool descriptions**: no change. They are derived from the CLI's help
  (`mcp/schema.py`), which names legends only by name; the `sources` tool
  is unchanged (section 4). Ola's own `sources.toml` entry `ar50-innlandet`
  already says "vegetation cover" and "tree type"; it is not in git and needs
  no edit.
- **ROADMAP.md**: the row for 91.

## 8. The red suite

No invariant-critical suite: no mutation round. Fixture figures as measured
on this branch's base.

1. **Loader, happy path** (`tests/python/test_legends.py`, new class): a TOML
   string with three classes (10, 30, 50) and two splits, listed 50 first:
   50 by `veg` (51 → 5051, 99 → 5099) and 30 by `tre` (31 → 3031,
   99 → 3099). `load` gives `splits` in code order, 30 then 50, with
   `codes == {"31": 3031, "99": 3099}` and `{"51": 5051, "99": 5099}` (99 in
   both, allowed); the palette has 0, 10, 30, 50, 3031, 3099, 5051, 5099 with
   their labels and colours; `classes` keys are only `"10"`, `"30"`, `"50"`.
2. **Loader refusals**, one test per case, each checking the `ValueError`
   names the file: a split in an uncoded legend; a split's `code` not listed;
   two splits of class 30; a split's `attribute` equal to `attribute` (case
   differs); a sub-class with `water = true`; two sub-classes of one split
   with value 51 and `"51.0"` (one key); a sub-class code 10 (a listed code);
   two sub-classes of one split with one code; one sub-class code (3031) in
   two splits.
3. **Shipped `ar50`** (amend `TestAr50`): header unchanged; palette is the
   nine-row table plus section 3's four forest and five open-firm-ground rows
   exactly, each colour its RGB in hex (3099: `#9ecc73`); `splits` is
   (30, `treslag`, `{"31": 3031, "32": 3032, "33": 3033, "99": 3099}`),
   (50, `vegetasjonsdekke`,
   `{"51": 5051, "52": 5052, "54": 5054, "55": 5055, "99": 5099}`);
   `classes` unchanged (nine keys).
4. **Coding a feature** (`tests/python/test_feature_input.py`), a GeoJSON
   source under `ar50` and the suite's small domain, one feature per line:
   arealtype 50 with vegetasjonsdekke 52 → 5052; 50 with 52.0 → 5052; 50 with
   `"52"` → 5052; 50 with 52.5 → 50; 50 with `true` → 50; 50 with 98 → 50; 50
   with no `vegetasjonsdekke` property → 50; 50 with treslag 31 and no
   vegetasjonsdekke → 50; 30 with treslag 31 → 3031; 30 with treslag 39 → 30;
   30 with vegetasjonsdekke 52 and no treslag → 30; 60 with treslag 31 and
   vegetasjonsdekke 52 → 60; arealtype 5051 → refused naming the file, the
   fid and 5051. Masks unchanged (`land_cover`). `split_key` alone, as a
   table: 52, 52.0, `"52"`, `" 52 "`, `"52.0"` → `"52"`; `"5a"` → `"5a"`;
   None, `True`, 52.5, NaN, `"nan"`, `[52]` → None.
5. **GeoPackage reader** (`tests/python/test_feature_input.py` or the
   GeoPackage reader's test file): fixture A
   (`tests/fixtures/ar50/ar50_sample.gml`, 6 features in fid order: arealtype
   60, 60, 60, 60, 60, 30; vegetasjonsdekke 98 on all six; treslag 31 on the
   forest feature and on four of the marsh features, 32 on one marsh
   feature), converted with `gml_to_geopackage`:
   `read_source(..., extra=("vegetasjonsdekke", "treslag"))` returns `extra`
   of six tuples, each `(98, t)` with `t` the feature's treslag, and `rows`
   equal to a read without `extra`. Under `ar50` the codes are 60 five times
   and **3031** for the forest feature. `extra=("VEGETASJONSDEKKE",)` reads
   the same 98s (case-insensitive); `extra=("vegetasjonsdekke",
   "no_such_column")` gives `(98, None)` six times, not a refusal
   (section 6); a read without `extra` gives `extra == ()`.
6. **GML 2 reader** (`tests/python/test_io_gml.py`): a two-feature document,
   `read_gml(stream, "a", extra=("b", "c"))`, `b` on the first feature only
   and `c` on neither: `(b's text, None)` and `(None, None)`. The GeoJSON
   branch of `read_source` likewise: a two-feature collection, one with
   property `b`, gives `extra == ((that value,), (None,))` for `extra=("b",)`.
7. **Through the CLI** (amend `TestAr50ThroughTheCli` in
   `tests/python/test_cli_legends.py`): `palette ar50` has `2 * 19`
   annotations, `"5052 Open firm ground, barren to non-continuous vegetation
   cover"` and `"3031 Forest, coniferous: at least 50 % of the forested area
   covered by conifers"` among them; a mesh of WEST arealtype 30 / treslag 31
   and EAST 50 / vegetasjonsdekke 52: triangle codes include 3031 and 5052,
   the border between them carries `land_cover` only, the record's
   `land_cover_codes` contains `"class 30 split by treslag into 3031, 3032,
   3033, 3099, class 50 split by vegetasjonsdekke into 5051, 5052, 5054,
   5055, 5099"`; the `.gpkg` triangles carry the palette's label and colour
   for each code. The same mesh with `treslag` and `vegetasjonsdekke`
   removed from the GeoJSON meshes, with codes 30 and 50 (section 6). The
   existing test there (WEST arealtype 30, EAST 81, no second attribute) is
   unchanged and still gets 30 and 81.

The callers that do not ask for extra attributes are covered by their
existing tests, unchanged (section 4). `TestAr50Listing` in
`tests/python/test_cli_sources_vector.py` is unchanged: the listing does not
change. No existing test meshes fixture A under `ar50` and checks its codes
(`git grep -n '"ar50"' tests`: the CLI legend tests above and the listing),
so no existing figure changes.

## 9. Left out

- **The comparison of the two forest approaches** (AR50's tree type alone;
  SR16's spruce within AR50's coniferous forest): Ola's later step of the
  forest plan, with SR16.
- **The other columns** (`skogbonitet`, `jordbruk`): not asked for. The
  format takes them as further `[[split]]` tables, one per class; a class
  split by two attributes at once is not designed here.
- **53 (lichen cover), forest's 39 and 98**: in no county file (section 2).
- **A split for land-cover rasters** (`--land-cover`): a raster holds one
  code per cell; a split legend passed there only lends its palette, as any
  coded legend does today.
- **Invalid AR50 polygons**: found while measuring, not changed here; handed
  back as an idea (1,716 of Innlandet's 287,730 polygons are invalid,
  self-intersecting; `rasputin mesh` with `ar50` stops with a raw GEOS
  `TopologyException` over a 20 km box at Bygdin and over Numedalslågen,
  where 54 polygons valid in EPSG:25832 also turn invalid when moved to
  EPSG:25833).
