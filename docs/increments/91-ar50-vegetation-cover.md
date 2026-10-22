# Increment 91: a legend splits one class by a second attribute (AR50's open firm ground by vegetation cover)

**Status: designed** by `@architect`, 2026-10-22, on `2600a740` (master),
branch `worktree-ar50-vegetation`. Full path: about 60 counted lines
(section 7), over the short path's 50; no C++ file, no refine or mesh code,
so no `@perf` acceptance run. Public.

Approved by Ola today, with the defaults: open firm ground (AR50 `arealtype`
50) split by vegetation cover (`vegetasjonsdekke`); forest not split now,
since spruce will come from NIBIO's forest resource map (SR16) later. Codex
found that the shipped `ar50` legend
(`src_python/tin_engine/data/legends/ar50.toml`) reads only `arealtype`, so a
mesh loses AR50's vegetation cover.

## 1. Prior art: legacy and literature

*Literature and standards.* No method is claimed and no claim of novelty is
made: a class refined by a second attribute is a table lookup. The published
precedent is NIBIO's own: its presentation rules for AR50 (*Presentasjonsregler:
Arealressurskart - AR50*, 12 June 2026,
`https://dokument.geonorge.no/tegneregler/ar50/20260615/presentasjonsregler-ar50.pdf`,
downloaded 2026-10-22 and read with `pdftotext`) draw the map layer
"Bonitet, 18 klasser" with queries of exactly this form, `"ARTYPE" = 50 AND
"ARVEGET" = 55`. Its table "Snaumark (ARVEGET)" colours the vegetation-cover
codes, and states that the attribute covers the same ground as arealtype 50.
The codes and their meaning: the schema `ar50.xsd` (SOSI AR50 2022-05-23,
`ArealressursVegetasjonsdekkeEnumerationType`) and NIBIO's product
specification, section 5.1.15 `ArealressursVegetasjonsdekke` (both cited in
`docs/increments/73-gml-vector-source.md`, section 1). Departure from NIBIO:
its 18-class layer lumps 54 with 99 and 52 with 53; this legend keeps each
code apart (section 3), and a model that wants NIBIO's lumping maps the codes
itself.

*Legacy.* Nothing. The greps, and the file lists they returned (none; exit
status 1):

```
$ git grep -n -i -l "vegetasjon\|arealtype\|ar50\|snaumark" legacy-archive -- legacy
$ git grep -n -i -l "sub.class\|subclass\|second attribute" legacy-archive -- legacy
```

*Earlier legend increments.* 41 (`docs/increments/41-class-maps-as-data.md`,
section 4) set the file format and its loader; 58 renamed it "legend"; 73
(section 4.5) shipped `ar50` over `arealtype` alone and listed the other
coded columns as left out (section 7 there). This increment extends 41's
format by one optional table; nothing in 41's mapping changes.

*Ideas.* Overlaps, not taken up: the private package's "Evaporation per
triangle from land-cover classes" (a consumer of these codes) and "spruce
apart from pine" (SR16; Ola: forest is not split now). Nothing in the core
ideas file is taken up.

## 2. What the data holds (measured)

`sqlite3` over Innlandet's GeoPackage
(`$RASPUTIN_DATA/AR50_UTM32_Innlandet_20261019/ar50.gpkg`, 287,730 features;
both columns `INTEGER`, every value `typeof` integer):

| `arealtype` | `vegetasjonsdekke` (count) |
|---|---|
| 50 | 51 (3,421), 52 (13,267), 54 (17,519), 55 (8,370), 99 (2,095) |
| every other (10, 20, 30, 60, 70, 81, 99) | 98 only |

The same query over the 14 other county files Ola downloaded today
(`AR50_Norway_20261022/ar50_{11,15,18,39,40,42,46,50,55,56}.gpkg`,
`AR50_UTM32_OsloViken_20261022/ar50_{03,31,32,33}.gpkg`): arealtype 50 always
has one of 51, 52, 54, 55, 99, and no other arealtype has anything but 98.
**53 (lichen cover) occurs in none of the 15 files**, though the schema
lists it.

## 3. The legend file: one optional `[split]` table

A coded legend may carry one `[split]`: the class it splits, the second
attribute, and one `[[split.class]]` per value of that attribute.

```toml
[split]
code = 50                        # the class split; a listed [[class]] code
attribute = "vegetasjonsdekke"   # the second attribute, read only for that class

[[split.class]]
value = 51                       # the second attribute's value, matched as text
code = 5051
label = "Open firm ground, no vegetation"
colour = "#c0c0c0"
```

Rules, each refused by the loader with a `ValueError` naming the file (as
every legend refusal is):

1. `[split]` only in a coded legend (one with `code_system`).
2. `split.code` is a listed `[[class]]` code.
3. `split.attribute` differs from `attribute` (case-insensitively, as the
   GeoPackage reader matches columns).
4. Every sub-class has `value`, `code`, `label`, `colour`, nothing else
   (`extra="forbid"`): no `water`, no `edges`.
5. The values are distinct as text (`str(value)`); the sub-class codes are
   distinct and none is a listed `[[class]]` code. This joins the existing
   "a class code is listed twice" check over the union of both lists.

**What a sub-class inherits.** Its borders carry the parent's edge names
(`_mask` reads the first attribute only, so it does not change); it is never
water unless the parent is. Its code, label and colour are its own, and all
three enter the palette beside the parent's.

**How a feature is coded.** The first attribute gives the code as today
(`_code`). If that code is `split.code`, the second attribute's value, as
text, is looked up among the sub-classes: found, the sub-class code; not found
(a value the split does not list, a null, a missing GeoJSON property), **the
parent's code stays**. A feature of any other class ignores the second
attribute. A first-attribute value whose code equals a sub-class code (an
unlisted value under `unlisted = "land"`, say a stray 5051) is refused naming
the file, the feature and the value, so a sub-class code never means two
things in one mesh.

**Numbering, so codes never collide.** The loader does not impose a scheme:
the file states each code, and rule 5 refuses a collision with any listed
code, and the run-time refusal above one with any unlisted code. The AR50
file numbers each sub-class `100 × 50 + vegetasjonsdekke`: 5051 to 5099. AR50's
arealtype codes are two digits (10 to 99, `ar50.xsd`), so a four-digit code
cannot be one; and the number reads back as its two parts.

**The 99 case** ("vegetation cover not registered", 2,095 polygons in
Innlandet): its own sub-class, 5099, not folded into the parent. NIBIO's
"Snaumark (ARVEGET)" table gives 99 a class and colour of its own, so this
copies the source rather than asserting a lumping; a model may still map 5099
as it maps 50. The parent code 50 stays listed for open firm ground whose
vegetation cover the split does not list (98, a value outside the schema's
codes, or none): in all 15 county files, no such polygon exists.

**53 is left out**: no polygon in any of the 15 files (section 2). If a later
release has one, it keeps the parent's code, 50, by the rule above; adding it
is one table in the data file.

### The shipped `ar50.toml`

Header and the nine `[[class]]` tables unchanged (so `code_system` stays
"AR50 arealtype code"), then:

```toml
[split]
code = 50
attribute = "vegetasjonsdekke"
```

and one `[[split.class]]` per row:

| `value` | `code` | `label` | `colour` | NIBIO's RGB ("Snaumark (ARVEGET)") |
|---|---|---|---|---|
| 51 | 5051 | Open firm ground, no vegetation | `#c0c0c0` | 192-192-192 |
| 52 | 5052 | Open firm ground, barren to non-continuous vegetation cover | `#ffe8c2` | 255-232-194 |
| 54 | 5054 | Open firm ground, barren to good vegetation cover | `#e1c790` | 225-199-144 |
| 55 | 5055 | Open firm ground, good to lush non-continuous vegetation cover | `#c7a971` | 199-169-113 |
| 99 | 5099 | Open firm ground, vegetation cover not registered | `#f8f7ba` | 248-247-186 |

- **Labels**: "Open firm ground" (the parent's label) and the English
  definition of each code in `ar50.xsd` and specification section 5.1.15,
  ASCII, so a label read alone in the GeoPackage says what it is.
- **Colours**: NIBIO's, converted from RGB. NIBIO's table labels 54 "Ikke
  registrert vegetasjonsdekke", the same words as 99; its 18-class table calls
  the 54 colour "middels frisk vegetasjon", so the row's label is NIBIO's slip
  and its colour is taken as 54's.
- The file's header comment gains one sentence naming the sources of the
  split's codes and colours, and this file, section 3.

## 4. Blueprint: data flow and boundaries

```
ar50.toml --tomllib--> legends._File (+ _Split, _SubClass)   validation, rules 1-5
                              | .legend()
                              v
feature_input.Legend (+ split: Split | None)  palette += sub-classes
                              |
read_source(path, layer, attribute, box_for, crs, second=legend.split.attribute)
   .gpkg  query_features(..., second)  one more column in the SELECT
   .gml   read_gml(stream, attribute, second)
   .json  p.get(second)
   rows: (fid, geometry, value, second value)   second None when not asked
                              |
_Tally._take --> _mask(value)            unchanged: edges from the first attribute
             --> _code(legend, value, what, second)   parent -> sub-class code
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
    attribute: str                 # the second attribute
    codes: Mapping[str, int]       # str(value) -> sub-class code

class Legend(BaseModel):
    ...
    split: Split | None = None
```

- **`legends.py`**: `_SubClass(_Strict)` (`value: int | str`, `code: Code`,
  `label`, `colour` with the existing pattern), `_Split(_Strict)` (`code: Code`,
  `attribute: str`, `classes: list[_SubClass] = Field(alias="class")`),
  `_File.split: _Split | None = None`; rules 1 to 5 in `_consistent`;
  `legend()` adds each sub-class to the palette and builds `Split`.
- **`feature_input.py`**: `Split`, `Legend.split`; `read_source` takes
  `second: str | None = None` and its rows become 4-tuples (`SourceRows.rows`,
  `read_lake_polygons`' unpacking follows); `_Tally.source` passes
  `legend.split.attribute` when there is a split; `_code` takes the second
  value and applies section 3's rule and refusal.
- **`io/geopackage.py`**: `query_features(..., second: str | None = None)`;
  `RawFeature.second: Any = None`; a missing second column is refused with the
  existing words ("layer ... has no column ...").
- **`io/gml.py`**: `read_gml(stream, attribute, second=None)`;
  `GmlFeature.second: str | None = None`. The GeoJSON branch reads
  `p.get(second)`.
- **`cli.py` `_land_cover`**: the `land_cover_codes` text gains
  `, class 50 split by vegetasjonsdekke into 5051, 5052, 5054, 5055, 5099`
  (generic: the split's code, attribute and sorted sub-codes) between the
  attribute and the legend name. The whole text stays ASCII
  (`checked_ascii` already guards it).
- **`source_report.py` `_coded`**: a fitting legend's split attribute is
  listed too, so `rasputin sources` (and the MCP `sources` tool) shows
  `vegetasjonsdekke`'s values beside `arealtype`'s.

Boundaries: no C++; nothing crosses into the core but the polygons and their
integer codes, as today; the legend stays data (TOML) checked once, the reader
stays generic (it reads a named column, not "AR50"). Each part is testable
alone: the loader from a TOML string, `_code` from a `Legend`, the readers
from a small file.

## 5. What a mesh carries

- **`.vtk` / `.ply`**: `land_cover_code` per triangle, now 5051 to 5099 on
  open firm ground with its cover registered, 50 where it is not given; the
  dataset string `land_cover_codes` with the split clause (section 4).
- **`.gpkg`**: the `triangles` layer's `code`, `land_cover` (label) and
  `colour`; the `land_cover` polygon layer per code; `land_cover_classes` with
  triangles and area per code; the QML style's categories. All come from the
  palette, so no writer changes.
- **`rasputin palette ar50`**: 15 entries (no-class, nine classes, five
  sub-classes).
- **The run record**: `land_cover_codes` as above; notice and licence
  unchanged.

**Size of the effect (measured, scratch probe, AC power).** Merging by class
keeps sub-classes apart, so their shared borders become constraint lines.
`rasputin mesh` on DTM10 (`DTM10_UTM33_20260925`), a 15 km square in
Jotunheimen (EPSG:25832, x 457,500 to 472,500, y 6,812,500 to 6,827,500; 62 %
of it arealtype 50), `--tolerance 5`, the current `ar50` against a copy of the
GeoPackage with arealtype 50 recoded `5000 + vegetasjonsdekke` and a legend
listing those codes (what this design computes):

| | `ar50` today | split |
|---|---|---|
| land-cover vertices after clean-up | 10,638 | 25,550 |
| output triangles | 143,549 | 156,671 (+9 %) |
| wall time | 2.35 s | 2.69 s |

A shapely-only probe of the land-cover vertices after merge and a 5 m border
simplification gave the same order: ×2.28 on a 20 km box at Bygdin (76 % open
firm ground), ×2.01 on a 50 km box there, ×1.00 on a 20 km lowland box near
Hamar (1 %). In mountains the land-cover lines about double; the mesh grows
by about a tenth at a 5 m tolerance.

## 6. A mesh made with the old legend

Unchanged on disk and still read by `read_mesh`; nothing migrates and no
alias or version key is added. Its `land_cover_codes` text has no split
clause, which is how a reader tells the two apart. Its codes stay valid under
the new legend: 50 is still listed with the same label and colour, so
`rasputin palette ar50` still colours it. A mesh with the split needs a
re-run. The private package has no code that maps AR50 codes
(`git grep -n -i "ar50\|arealtype" private/master` finds only the ideas file,
research notes and one research script reading the GeoPackage directly).

## 7. Size, and where the prose changes

| File | Counted lines |
|---|---|
| `legends.py`: `_SubClass`, `_Split`, the field, rules 1 to 5, palette and `Split` in `legend()` | about 25 |
| `feature_input.py`: `Split`, `Legend.split`, `read_source`'s `second`, `_take`'s unpacking, `_code`'s lookup and refusal, `read_lake_polygons` | about 20 |
| `io/geopackage.py`: the second column | about 6 |
| `io/gml.py`: the second element | about 3 |
| `cli.py` `_land_cover`: the split clause | about 3 |
| `source_report.py` `_coded`: the split attribute | about 2 |
| `data/legends/ar50.toml` | 0 (TOML is not counted) |
| Total | **about 60**, under the 700 of `CLAUDE.md` section 2 |

Prose, in the same PR (not counted):

- **README.md**, the "Features as constraints" item: the `ar50` legend reads
  its land type (`arealtype`) and splits open firm ground by vegetation cover
  (`vegetasjonsdekke`) into codes 5051 to 5099.
- **docs/GLOSSARY.md**, *Legend*: one sentence: a legend may split one class
  by a second attribute into sub-classes with codes of their own.
- **project_structure.md**: the `data/legends/` row says `ar50` splits 50;
  the `legends.py` row names the split.
- **MCP tool descriptions**: no change. They are derived from the CLI's help
  (`mcp/schema.py`), which names legends only by name, and the `sources` tool
  lists the split attribute's values through `_coded`. Ola's own
  `sources.toml` entry `ar50-innlandet` already says "vegetation cover"; it is
  not in git and needs no edit.
- **ROADMAP.md**: the row for 91.

## 8. The red suite

No invariant-critical suite: no mutation round. Fixture figures as measured
on this branch's base.

1. **Loader, happy path** (`tests/python/test_legends.py`, new class): a TOML
   string with two classes (10, 50), `[split] code = 50`, attribute
   `veg`, sub-classes 51 → 5051 and 99 → 5099. `load` gives `split.codes ==
   {"51": 5051, "99": 5099}`; the palette has 0, 10, 50, 5051, 5099 with their
   labels and colours; `classes` keys are only `"10"` and `"50"`.
2. **Loader refusals**, one test per rule 1 to 5, each checking the
   `ValueError` names the file: split in an uncoded legend; `split.code` not
   listed; `attribute` equal to the split's (case differs); a sub-class with
   `water = true`; two sub-classes with value 51 (one as integer, one as text
   `"51"`); a sub-class code 10 (a listed code); two sub-classes with one code.
3. **Shipped `ar50`** (amend `TestAr50`): header unchanged; palette is the
   nine-row table plus section 3's five rows exactly, each colour its RGB in
   hex; `split` is code 50, attribute `vegetasjonsdekke`, codes
   `{"51": 5051, "52": 5052, "54": 5054, "55": 5055, "99": 5099}`;
   `classes` unchanged (nine keys).
4. **Coding a feature** (`tests/python/test_feature_input.py`), a GeoJSON
   source under `ar50` and the suite's small domain: arealtype 50 with
   vegetasjonsdekke 52 → 5052; 50 with 98 → 50; 50 with no
   `vegetasjonsdekke` property → 50; 30 with 52 → 30; arealtype 5051 →
   refused naming the file, the fid and 5051. Masks unchanged (`land_cover`).
5. **GeoPackage reader** (`tests/python/test_feature_input.py` or the
   GeoPackage reader's test file): fixture A
   (`tests/fixtures/ar50/ar50_sample.gml`, 6 features: arealtype 30 once and
   60 five times, vegetasjonsdekke 98 on all six), converted with
   `gml_to_geopackage`: `read_source(..., second="vegetasjonsdekke")` returns
   98 as each row's fourth item; the codes under `ar50` are still 30 and 60.
   A `second` naming no column is refused naming the column.
6. **GML 2 reader** (`tests/python/test_io_gml.py`): a two-feature document,
   `read_gml(stream, "a", second="b")` gives each feature's `b` text, None
   where absent.
7. **Through the CLI** (amend `TestAr50ThroughTheCli` in
   `tests/python/test_cli_legends.py`): `palette ar50` has `2 * 15`
   annotations, `"5052 Open firm ground, barren to non-continuous vegetation
   cover"` among them; a mesh of WEST arealtype 50 / vegetasjonsdekke 52 and
   EAST 50 / 55: triangle codes include 5052 and 5055, the border between
   them carries `land_cover` only, the record's `land_cover_codes` contains
   `"class 50 split by vegetasjonsdekke into 5051, 5052, 5054, 5055, 5099"`;
   the `.gpkg` triangles carry the palette's label and colour for each code.
8. **Listing** (amend `TestAr50Listing.test_14` in
   `tests/python/test_cli_sources_vector.py`): fixture A's `classes` gains
   `{"attribute": "vegetasjonsdekke", "values": ["98"], "distinct": 1}` after
   `arealtype`'s entry, in the file's column order.

## 9. Left out

- **Forest by tree type (`treslag`)**, and the other columns: Ola, forest is
  not split now (SR16 later). The format allows one split per legend; a
  second would be a list of `[split]` tables, not designed here.
- **53 (lichen cover)**: in no county file (section 2).
- **A split for land-cover rasters** (`--land-cover`): a raster holds one
  code per cell; a split legend passed there only lends its palette, as any
  coded legend does today.
- **Invalid AR50 polygons**: found while measuring, not changed here; handed
  back as an idea (1,716 of Innlandet's 287,730 polygons are invalid,
  self-intersecting, and `rasputin mesh --features ar50.gpkg` over a 20 km box
  at Bygdin stops with a raw GEOS `TopologyException`).
