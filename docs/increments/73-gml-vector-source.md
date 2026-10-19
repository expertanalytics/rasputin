# Increment 73: GML 3.2 land-resource maps as a vector source

**Status: designed** by `@architect`, 2026-10-19, on `9ba6c5ec` (master,
increment 72 merged as #294), branch `worktree-gml-vector`. Full path: about
190 counted lines (section 8), no C++ file, no refine or mesh code, so no
`@perf` acceptance run.

Ola, 2026-10-19: "Add DOM10, start on the GML reader asap." The file is
NIBIO's AR50 land-resource map for Innlandet,
`$RASPUTIN_DATA/AR50_UTM32_Innlandet_20261019/34_25832_ar50_gml.gml`
(1,766,264,204 bytes, NLOD licence). The aim: register it in Ola's
`sources.toml` and use it wherever a vector source is used today
(`mesh --features`, `rasputin sources`, the MCP tools).

## 1. Prior art: legacy and literature

*Literature and standards.* The format is OGC GML 3.2.1 (ISO 19136:2007);
the file follows the SOSI application schema AR50 2022-05-23
(`https://skjema.geonorge.no/SOSI/produktspesifikasjon/AR50/20220523/ar50.xsd`).
The target is OGC GeoPackage 1.3 (OGC 12-128r18), which
`io/geopackage_write.py` already writes. Streaming a large XML document with a
pull parser that drops each finished element (`xml.etree.ElementTree.iterparse`
with the element cleared) is standard practice and the `geospatial-data-formats`
skill's rule. No method is claimed, and no claim of novelty is made.

*Open readers checked.* OGR's GML driver (GDAL) reads this file, but GDAL and
OGR are prohibited (`CLAUDE.md` section 2). `pygml` 0.3.2 (MIT, geopython;
PyPI metadata: `requires_dist` is `lxml`) parses a single GML 3.1/3.2
geometry element to `__geo_interface__`; its documentation shows geometry
parsing only (no features, attributes or streaming of a document), and it
would add `lxml`. The ring parsing it would save is about 15 lines, so nothing
is adopted. Departures from OGR's behaviour, which
users of QGIS will compare against: the feature's `gml:id` is kept as a column
`gml_id` (as OGR does); nested properties are not flattened (OGR flattens
them), see section 7.

*Legacy.* The legacy tree read CORINE's GML2 with `lxml`, whole-document.
The grep, and the files it returned:

```
$ git grep -n -i -l "gml\|posList" legacy-archive -- legacy
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/tests/test_gml_repository.py
```

`gml_repository.py` reads `gml:outerBoundaryIs` / `gml:coordinates` (GML2)
only; increment 16b already carried that across as `io/gml.py`'s `read_gml`.
Nothing of GML 3.2 is in the legacy tree, so nothing more is carried across.

## 2. What the file holds (measured)

Counted over the whole file with streaming `grep -o` (5.9 s) and a streaming
parse (`docs/increments/73-probes/parse_probe.py`, 21.6 s):

| Figure | Value |
|---|---|
| Features (`gml:featureMember`, one `app:ArealressursFlate` each) | 287,730 |
| Geometry | one `gml:Polygon` per feature (no `gml:Surface`, no `gml:MultiSurface`), `gml:exterior`/`gml:interior` rings of `gml:LinearRing` with `gml:posList`; no `srsDimension` (2-D) |
| Interior rings (holes) | 44,605 |
| Positions | 40,865,758; largest feature 159,686 |
| `srsName` | `http://www.opengis.net/def/crs/EPSG/0/25832` on every polygon (UTM 32N, ETRS89; `parse_crs` reads it as EPSG:25832, axes east, north) |
| `gml:boundedBy` | absent: the extent is known only after a full read |
| Longest exterior edge | 4,395 m |

Attributes, with the schema's code lists (`ar50.xsd`) and the counts in the
file:

| Element | Meaning | Codes in the file (count) |
|---|---|---|
| `arealtype` | land type | 10 developed (5,007), 20 agriculture (16,426), 30 forest (134,454), 50 open firm ground (44,672), 60 marsh (78,442), 70 snow and glacier (382), 81 fresh water (7,981), 99 not surveyed (366); 82 ocean is in the schema, not in this file |
| `treslag` | tree type (crown cover) | 31 coniferous (126,402), 32 deciduous (13,243), 33 mixed (19,556), 39 no tree cover (84,126), 98 not relevant (30,162), 99 not registered (14,241) |
| `skogbonitet` | forest site quality | 11 impediment (112,757), 12 low (38,114), 13 medium (29,722), 18 high (18,286), 98 (74,834), 99 (14,017) |
| `vegetasjonsdekke` | vegetation cover above the tree line | 51 (3,421), 52 (13,267), 54 (17,519), 55 (8,370), 98 (243,058), 99 (2,095) |
| `jordbruk` | agriculture subclass | 24 (11,209), 25 (5,217), 98 (271,304) |
| `kartstandard`, `informasjon`, `oppdateringsdato` | map standard (`AR50` everywhere), a note, the update time | text |

AR50 has **no soil field**: the schema's `ArealressursFlate` has `område`,
`arealtype`, `skogbonitet`, `treslag`, `jordbruk`, `vegetasjonsdekke`,
`kartstandard` and the common properties, nothing on soil. NIBIO's product
sheet says AR50 suits visualisation, not geometric analysis or area statistics
below the tree line, and areas under 15 decares are merged into their
neighbours; the catalogue description in section 6 says so.

## 3. The decision: a one-time conversion to GeoPackage

Two layouts were weighed on the real file:

| | Read the `.gml` on every run | Convert once to `.gpkg` |
|---|---|---|
| `mesh --features`, 10 km box | stream all 1.77 GB, about 22 s, every run | R-tree query: 1,054 rows in 0.10 s (`read_probe.py`) |
| `rasputin sources` listing | a full read for the extent, count and classes (no `gml:boundedBy`) | `gpkg_contents` and `distinct_values`, as today |
| Code that changes | `feature_input.read_source`, `source_report._vector`, `sources.py`'s suffixes | none of them |
| New code | a streaming reader plus a box filter in the read path | a streaming reader plus a converter and a command |

**Converting once wins**: every consumer of vector sources already reads
GeoPackage through an R-tree, and the arcgis-features fetch (increment 63)
already stores what it fetches as a GeoPackage through the same writer
functions (`create_features_table`, `insert_features`). The converted file is
an ordinary local `.gpkg` source: no new catalogue kind, no new reader path.

The full conversion was run as a probe through the existing writer
(`docs/increments/73-probes/convert_probe.py`, on AC power): 287,730 features
in **27.0 s, peak resident memory 126 MiB**, an 832,618,496-byte GeoPackage
that `layer_info`, `layer_summary`, `distinct_values` and
`feature_input.read_source` read unchanged (layer `ArealressursFlate`,
EPSG:25832, R-tree present; arealtype's 8, treslag's 6 and skogbonitet's 6
distinct values as in section 2).

## 4. Blueprint: data flow and boundaries

```
rasputin gml-to-gpkg SOURCE --out OUT [--crs TEXT]          (cli.py, L5)
  opens SOURCE "rb", OUT.part via io.repository.create_sqlite
        |
        v
gml_import.gml_to_geopackage(stream, conn, crs=None)        (new, L4)
  checks the CRS (crs.parse_crs, same_crs), types the columns,
  writes in chunks of 1000 through io.geopackage_write
        ^
        | Iterator[Gml32Feature], one at a time
io.gml.iter_gml32(stream)                                    (io/gml.py, L2)
  iterparse, element cleared after each gml:featureMember;
  imports no tin_engine module (NumPy, shapely, stdlib)
        |
OUT.part renamed onto OUT  -> an ordinary local .gpkg source:
  mesh --features OUT --features-legend ...   (feature_input, unchanged)
  rasputin sources                            (source_report, unchanged)
```

Data and execution stay apart: `iter_gml32` is a pure decoder of a byte
stream into frozen values; `gml_to_geopackage` is a pure function of a stream
and a connection (opens nothing, knows no path); the command owns the files.
All three are blocking (expat, GEOS, SQLite); an async caller runs
`gml_to_geopackage` in `asyncio.to_thread`. No state is kept between calls.
The MCP server derives the new command as a tool automatically (every shown
command is one); `--out` is marked `writable=True`, so the output root fences
it, and `SOURCE` is judged by the HTTP read limit like every other input.

### 4.1 `io/gml.py`: the GML 3.2 decoder

Added beside `read_gml`, sharing `GmlError`:

```python
GML32 = "{http://www.opengis.net/gml/3.2}"

@dataclass(frozen=True, slots=True)
class Gml32Feature:
    gml_id: str
    type_name: str                  # the feature element's local name
    geometry: Polygon | MultiPolygon | None
    srs_name: str | None            # the geometry's srsName, inherited by a MultiSurface's members
    values: Mapping[str, str | None]  # leaf properties by local name, document order

def iter_gml32(stream: IO[bytes]) -> Iterator[Gml32Feature]: ...
```

Rules, each a `GmlError` naming the `gml:id` where there is one:

- **Members.** A feature is the one element inside a `gml:featureMember`
  (GML 3.2 namespace). After it is yielded the root is cleared, so memory
  holds one feature at a time. A document with no `gml:featureMember` is
  refused ("no GML 3.2 gml:featureMember"): a WFS 2.0 `wfs:member` document,
  or a GML2 one, is not read silently as empty. Malformed XML is refused
  naming the parser's line, as `read_gml` does.
- **Geometry.** The geometry is the child of a property element of the
  feature that is a `gml:Polygon`, `gml:Surface` or `gml:MultiSurface`.
  `gml:Polygon` is a shapely `Polygon` of its `gml:exterior` ring and every
  `gml:interior` ring (holes kept, in order). `gml:Surface` is its
  `gml:patches/gml:PolygonPatch` children read the same way: one patch a
  `Polygon`, more a `MultiPolygon`. `gml:MultiSurface` is a `MultiPolygon` of
  every `gml:surfaceMember` and `gml:surfaceMembers` child, a member Surface's
  patches flattened in. Rings are not repaired, reoriented or validated:
  stored as served, as the arcgis fetch stores CORINE. A feature with two
  surface properties is refused; a curve or point geometry (`gml:LineString`,
  `gml:Curve`, `gml:Point`, `gml:MultiCurve`, ...) is refused ("surfaces
  only; gml:LineString is not read"). A feature with no geometry is yielded
  with `geometry=None`.
- **Rings.** A ring is a `gml:LinearRing` with one `gml:posList`. Its
  dimension is the `srsDimension` of the `posList`, else of the geometry
  element, else 2; the third coordinate onwards is dropped. A count not a
  multiple of the dimension, fewer than four positions, or a coordinate that
  is not finite is refused. `gml:pos` and `gml:coordinates` inside a
  `LinearRing`, and `gml:Ring` (curve members, arcs), are refused naming the
  element.
- **Values.** Every child of the feature with no element children is a
  value: its local name, its text stripped, None when empty or `xsi:nil`.
  Nested properties (`identifikasjon`, `kopidata`, `kvalitet`) are skipped.

`read_gml` (GML2) reads the namespace declarations before the tree
(`iterparse` events `start-ns` then `start`) and refuses a document that
declares the GML 3.2 namespace: "a GML 3.2 file: convert it with `rasputin
gml-to-gpkg` and give the .gpkg". Today such a file fails only after
`ET.parse` has built the whole tree, which for AR50 is several GB, or, with
`--features-crs`, reads as zero features. The rest of `read_gml` keeps its
behaviour (the tree is still built whole: OGR's GML2 files are small).

### 4.2 `gml_import.py`: the converter

```python
@dataclass(frozen=True, slots=True)
class GmlImport:
    layer: str                      # the feature type, also the table name
    crs: str                        # as stored: the srsName text, or the given crs
    features: int                   # rows written
    without_surface: int            # features skipped for having no geometry
    columns: tuple[str, ...]        # gml_id first, then the value columns

#: Features per transaction, a chosen value: 1000 kept the probe's peak at
#: 96 MiB on 16,140 features and 126 MiB on AR50's 287,730.
CHUNK = 1000

def gml_to_geopackage(stream: IO[bytes], conn: sqlite3.Connection, *, crs: str | None = None) -> GmlImport: ...
```

- **One table.** The first feature names the table (`type_name`) and its
  columns: `fid INTEGER PRIMARY KEY` (1, 2, ... in document order), `geom
  MULTIPOLYGON` (the writer wraps a `Polygon` as one part), `gml_id TEXT`,
  then one column per value of the first feature, in its order: `INTEGER`
  when the text is a decimal integer (`[+-]?[0-9]+`), else `TEXT`. Values are
  inserted as text; SQLite's integer affinity stores `"60"` as 60 and keeps
  any later non-integer text as text, so no value is lost or refused for its
  type. A later feature of another type, or with a value the columns lack, is
  refused naming its `gml_id` and the name; a missing value is NULL. A value
  named `fid`, `geom` or `gml_id` is refused. The column types are passed as
  the esri names `create_features_table` takes
  (`esriFieldTypeInteger`, `esriFieldTypeString`).
- **The CRS.** Every geometry's `srsName` must be the same text; a second one
  is refused naming the feature and both names. With none anywhere, `crs` is
  required ("no geometry names its srsName; give --crs"); with both, they must
  be `same_crs`. The CRS must parse (`parse_crs`) and its first axis must
  point east: a CRS listed latitude first (EPSG:4326 as an OGC URI, where GML
  3.2 writes latitude, longitude) is refused, not swapped, since no fixture
  of such a file exists. The GeoPackage stores `crs.srs_entry` of it
  (EPSG:25832 for AR50). Nothing is reprojected here.
- **Writes.** `create_geopackage` and `create_features_table` on the first
  feature, then `insert_features` every `CHUNK` features and at the end
  (each a transaction; the R-tree and the layer bounds grow with them). A
  document with no feature with a surface is refused, so no empty layer is
  written.

### 4.3 The command

`rasputin gml-to-gpkg SOURCE --out OUT [--crs TEXT]`, help "Convert a GML
3.2 file of surfaces to a GeoPackage, once, for --features and the
catalogue". `OUT` must end in `.gpkg`; it is written as `OUT.part` (a stale
one removed first) and renamed onto `OUT` only when the conversion finished,
as `_write_geopackage` does; any refusal removes the `.part` and is a usage
error (exit 2) naming SOURCE and the feature. On success one line on stderr,
in plain words:

```
ar50.gpkg: 287730 features of ArealressursFlate in EPSG:25832, 0 without a surface; columns gml_id, informasjon, oppdateringsdato, arealtype, skogbonitet, treslag, jordbruk, vegetasjonsdekke, kartstandard
```

and `OUT` on stdout, as `palette --out` prints its file.

### 4.4 The catalogue

`sources.py` keeps its suffixes. A local entry whose path ends in `.gml` is
refused at load with "a .gml file is read after `rasputin gml-to-gpkg`
converts it; register the .gpkg". The shipped `sources.toml` header comment
and the module docstring say the same in one sentence each.

## 5. The CRS: EPSG:25832 into a UTM 33 run

Nothing new: `feature_input` already moves a projected source's vertices into
the DEM's CRS one by one and pre-clips there with no widening, which it states
is exact for any pair but geographic-over-Transverse-Mercator. The only error
is that a straight edge in UTM 32 is slightly curved in UTM 33. Measured with
pyproj on a 4.4 km edge (AR50's longest) at UTM 32 eastings 500,000 and
700,000 m, northing 6,800,000 m (about 61° N): the moved edge strays at most
**12.2 mm and 13.0 mm** from its chord (`read_probe.py`). AR50 claims 2 m or
worse; the stray is negligible and needs no densifying.

## 6. What the catalogue lists

Ola's own `$RASPUTIN_DATA/sources.toml` (not in git; written by the main
session once the PR merges, since it is Ola's data file) gains:

```toml
[source.ar50-innlandet]
kind = "local"
path = "AR50_UTM32_Innlandet_20261019/ar50.gpkg"
credit = "Inneholder data under Norsk lisens for offentlige data (NLOD) tilgjengeliggjort av NIBIO. Kilde: NIBIO, AR50 Arealressurskart 1:50 000, Innlandet, lastet ned 2026-10-19; konvertert til GeoPackage av rasputin."
licence_note = "Norsk lisens for offentlige data (NLOD) 1.0, https://data.norge.no/nlod/no/1.0."
description = "NIBIO's land-resource map at 1:50 000 for Innlandet: land type (arealtype), tree type (treslag), forest site quality (skogbonitet), vegetation cover and agriculture subclasses, in EPSG:25832. Areas under 15 decares are merged into neighbours; NIBIO: for visualisation, not geometric analysis or area statistics below the tree line."
```

The licence is NLOD 1.0, as data.norge.no states for NIBIO's AR50 land-type
dataset (its distributions' licence field: `NLOD_1_0`). `rasputin sources` then
lists, from the GeoPackage's headers: CRS EPSG:25832, geometry MultiPolygon,
287,730 features, the nine columns of section 4.2 as attributes, and the
coverage box. `classes` stays empty until a coded legend names `arealtype`
(section 7, question 1).

## 7. Left out

- **An AR50 legend file** (codes to labels, colours and edge names, which
  would make `classes` list and let `--features-legend ar50` code triangles):
  a coded legend needs `positional_accuracy_m`, which drives the clean-up
  distances, and NIBIO's product sheet says only "2 and worse". A question
  for Ola (handback), and its own short-path increment (TOML only, 0 counted
  lines).
- **Reading from the `.zip`** Ola downloaded (536 MB): `zipfile` streams a
  member, so `SOURCE` could be a zip; not asked for.
- **Nested properties** (`lokalId` inside `identifikasjon`): NIBIO's object
  number. Flattening rules are a design of their own; `gml_id` identifies the
  feature.
- **Curves and points**, `gml:pos` rings, arcs, latitude-first CRSs, several
  feature types per file, WFS 2.0 `wfs:member`: refused in words, not read.
- **Reprojection at conversion** (`--to-crs`): the run moves vertices already.
- **A direct `.gml` read on every run** (section 3).

## 8. Size

| File | Counted lines |
|---|---|
| `io/gml.py`: `Gml32Feature`, `iter_gml32`, ring, polygon, surface and multi-surface helpers, the namespace refusal in `read_gml` | about 80 |
| `gml_import.py`: `GmlImport`, `gml_to_geopackage`, column typing, the CRS checks | about 70 |
| `cli.py`: the command | about 35 |
| `sources.py`: the `.gml` refusal | about 3 |
| Total | **about 190**, under the 700 of `CLAUDE.md` section 2 |

Also, not counted: `project_structure.md` rows for `gml_import.py` and the
new part of `io/gml.py`; the `sources.toml` comment; `test_layering.py` gains
`"gml_import": "crs io.geopackage_write io.gml io.repository"` at L4 (or
without `io.repository` if the developer passes no helper from it) and
`gml_import` in `cli`'s list; `io.gml` stays `""`.

## 9. The red suite

No invariant-critical suite: no mutation round. New files
`tests/python/test_io_gml32.py`, `tests/python/test_gml_import.py`; amended
`test_sources_vector.py`, `test_io_gml.py`, `test_mcp_schema.py`,
`test_layering.py`, and the CLI test file the tester picks.

**Fixture A, cut from the real file:** `tests/fixtures/ar50/ar50_sample.gml`,
the real file's first 720 bytes (declaration and the root element with its
namespaces, verbatim), an XML comment after the declaration carrying the NLOD
credit of section 6, then six whole `gml:featureMember` elements copied
byte for byte, then `</gml:FeatureCollection>`. About 24 KB. The cut is
reproducible with `docs/increments/73-probes/slice.py` and `pick.py`. Its
figures, measured (`pick.py` on the first 100 MB):

| `gml:id` | holes | exterior positions | interior positions | area m² | arealtype | treslag | skogbonitet |
|---|---|---|---|---|---|---|---|
| `id551e0bc0-b63d-45b7-a7ba-0ecdff9d7c95` | 0 | 37 | 0 | 15,315.0 | 60 | 32 | 11 |
| `idd58d2d03-6841-4a8e-8043-10e00987a82a` | 0 | 10 | 0 | 815.1 | 60 | 31 | 13 |
| `idfc6f8e2f-9c4b-445a-b840-1693a8483aec` | 0 | 15 | 0 | 4,899.6 | 60 | 31 | 11 |
| `id9dd0239a-c844-4fc2-a42a-1a6e70666a08` | 0 | 50 | 0 | 15,966.2 | 60 | 31 | 11 |
| `id40a10b70-c9f7-4ac4-b74b-30bd33ed551d` | 0 | 55 | 0 | 27,579.9 | 60 | 31 | 11 |
| `id0ed7f868-8d50-4dc3-bb29-53b1d9ecff7d` (feature 12,843 of the file) | 2 | 192 | 60 | 353,418.5 | 30 | 31 | 18 |

All six are valid shapely polygons.

**Fixture B, hand-made strings in the tests**, each a few hundred bytes:
a `gml:Surface` with two `gml:PolygonPatch`; a `gml:MultiSurface` of two
`gml:surfaceMember` polygons, one with a hole; a `posList` with
`srsDimension="3"`; and one string per refusal below.

Tests:

1. *Fixture A decodes.* `iter_gml32` yields six features in order with the
   `gml_id`s, hole counts, position counts (exterior including the closing
   position) and areas of the table (areas to 0.1 m²), `type_name`
   `ArealressursFlate`, `srs_name` the EPSG 25832 URI, and `values` exactly
   `informasjon, oppdateringsdato, arealtype, skogbonitet, treslag, jordbruk,
   vegetasjonsdekke, kartstandard` in that order (no `identifikasjon`,
   `kopidata`, `kvalitet`).
2. *Surfaces map.* Fixture B's Surface is a `MultiPolygon` of 2 parts; the
   MultiSurface a `MultiPolygon` of 2 parts with 1 hole; a one-patch Surface
   a `Polygon`; the 3-D `posList` gives 2-D coordinates with the same x, y.
3. *Refusals, each a `GmlError` naming the `gml:id`:* a curve geometry, two
   surface properties, a `gml:pos` ring, a `gml:Ring`, a count not a multiple
   of the dimension, a three-position ring, a `NaN` coordinate, malformed XML
   (naming the line), no `gml:featureMember` at all (a GML2 document and a
   `wfs:member` one).
4. *Memory is bounded.* 20,000 synthetic features (a 5,229,019-byte
   document, the generator of `mem_probe.py`) consumed through `iter_gml32`
   under `tracemalloc`: peak under **1 MiB**. Measured with the element
   cleared: 167 KiB (236 KiB at 2,000 features); without the clear, 30,154 KiB
   at 20,000 and 3,149 KiB at 2,000, so the test fails on a reader that keeps
   the tree. 1 MiB is a chosen bound, about 6 times the measured peak.
5. *Round trip through the existing readers.* `gml_to_geopackage` of fixture
   A into a fresh connection, closed, then read back: `layer_info` gives table
   `ArealressursFlate`, CRS EPSG:25832, an R-tree; `layer_summary` 6 features,
   `MULTIPOLYGON`, attributes `gml_id` then the eight of test 1;
   `feature_input.read_source(path, None, "arealtype", box)` with a box
   holding all six returns 6 rows whose areas match the table and values 60,
   60, 60, 60, 60, 30; `distinct_values(..., "arealtype", 100)` gives
   `("30", "60")` and 2. The returned `GmlImport` says 6 features, 0 without
   a surface.
6. *Converter rules.* Column types: `arealtype` INTEGER, `informasjon` TEXT
   (from `PRAGMA table_info`); a later value the first feature lacked, a
   second feature type, a value named `gml_id`, two `srsName`s, no `srsName`
   and no `crs`, a `crs` that differs from the `srsName`, an EPSG:4326 URI
   (latitude first): each refused naming what it names in section 4.2; no
   `srsName` with `crs="EPSG:25833"` stores srs_id 25833; a feature with no
   geometry is skipped and counted; an `xsi:nil` value is NULL; a document
   whose features all lack a surface is refused.
7. *The GML2 reader refuses GML 3.2 early.* `read_gml` on fixture A raises a
   `GmlError` whose text names `rasputin gml-to-gpkg`; the 200-feature CORINE
   GML2 fixture still reads as before (the existing suite).
8. *The catalogue.* A local entry with path `/data/ar50.gml` is refused with
   the text of section 4.4, whatever its `values`; the existing
   `test_features_need_a_vector_path` case `lines.gml` is amended to that
   text.
9. *The command.* `rasputin gml-to-gpkg A --out x.gpkg` exits 0, prints
   `x.gpkg` on stdout and the section 4.3 line (with `6 features of
   ArealressursFlate in EPSG:25832`) on stderr; an existing `x.gpkg` is
   replaced; `--out x.sqlite` is a usage error; a refused document exits 2
   leaving neither `x.gpkg` nor `x.gpkg.part`.
10. *MCP pins.* `TOOL_NAMES` gains `gml_to_gpkg` in `rasputin --help` order,
    `REAL_POSITIONALS` gains `("gml_to_gpkg", "source", {"out": "x.gpkg"})`,
    `WRITABLE` gains `("gml_to_gpkg", "out")`.
