# Increment 63: `rasputin fetch` downloads lake polygons and land cover for a box

**Status: designed** (branch `worktree-fetch-vectors`, from master `3cb6d416`); not yet reviewed.

## 1. What Ola asked

Ola wants rasputin to fetch and prepare its own input data, so that a fresh
machine gets what Codex found by hand in Ola's `rasputin_data` folder.
Increment 57 (`docs/increments/57-vector-sources.md`, merged as #275) put
vector files into the source catalogue as local entries, each listed with the
options it fits. Codex's run used two such files:

- `rasputin_data/nve_hrd/lakes.geojson`: 81 lake polygons near the HRD
  gauging stations (HRD is NVE's Hydrological Reference Dataset), written by
  `rasputin fetch-stations nve-hrd`; not all of Norway's lakes.
- `rasputin_data/corine2018_dtm10_utm33.gpkg`: CORINE Land Cover 2018
  polygons cut from the Europe-wide file (version v2020_20u1) by a one-off
  script (`corine2018_dtm10_utm33.README.txt`: "Later this becomes a tested
  rasputin command (Ola, 2026-09-28)").

Ola's yes (2026-10-17) to "fetchable vector data: make NVE's lake polygons
and CORINE land cover sources that rasputin fetch downloads and registers",
with a licence question for him (section 8).

What is missing today, checked against `3cb6d416`:

1. **The remote kinds are rasters only.** `RemoteSource.kind` is
   `"local" | "one-cog" | "cog-tiles" | "one-cog-per-year"`, and every remote
   kind is a Cloud-Optimised GeoTIFF read by byte ranges
   (`src_python/tin_engine/fetch/run.py`, `_header` parses a TIFF prefix).
   `_consistent` refuses `values = "features"` for any kind but `local`
   ("values: features need a local source").
2. **The cache holds raster blocks only.** `CacheManifest.objects` are
   `CachedObject`s (content length, header hash, block shape), and
   `CacheRepository` reads TIFF blocks.
3. **The one network reader of NVE's lakes asks per station.**
   `fetch/nve.py` queries `Innsjodatabase2/MapServer/5` in a 100 m envelope
   round each HRD station; there is no request for a box.

## 2. Prior art: legacy and literature

Written before the design.

*Literature and practice.*

- **The ArcGIS REST API, Query (Map Service/Layer)**
  (developers.arcgis.com/rest/services-reference/enterprise/query-map-service-layer/,
  read 2026-10-17). The design uses four of its rules, quoted: `returnIdsOnly`
  ("there is no limit on the number of object IDs returned in the ID array
  response"); `objectIds` ("The object IDs of this layer or table to be
  queried"); `outSR` ("If outSR is not specified, the geometry is returned in
  the spatial reference of the map"), so every request names it; and the
  caution to "Avoid using objectIds when a layer/table reports true for
  uniqueIdInfo.OIDFieldContainsHashValue" (both layers report no
  `uniqueIdInfo`; probed). The page's own paging route, `resultOffset` with
  `exceededTransferLimit`, is not used: section 4.3 says why.
- **OGC API Features, Part 1** is the open standard for the same job (items
  of a collection by `bbox`, paged by `next` links). Neither provider offers
  it for these layers (NVE's Geonorge record lists no WFS or OGC API; the
  CLMS product page names a WMS and the ArcGIS REST service only). Nothing is
  taken from it.
- **The GeoPackage specification** (OGC 12-128r18), as `io/geopackage.py`
  and `io/geopackage_write.py` already follow it: a features table with an
  integer primary key, the binary header of section 2.1.3 before ISO WKB,
  and the R-tree of Annex F.3.
- **STAC's Collection extent**, as 55 and 57 use it: the listing's
  coverage box.

Departures: none from the API's documented rules. No novelty is claimed.

*Legacy.*

```
$ git grep -l -i -E 'arcgis|objectIds|resultOffset|innsjo|corine|discomap|returnIdsOnly' legacy-archive -- legacy
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/tests/test_gml_repository.py
```

Read: every hit is the legacy tree reading CORINE from a local GML archive
(`application.py` line 59, `Path(os.environ["RASPUTIN_DATA_DIR"]) /
"corine"`; line 127, `gml_repository.GMLRepository(path=corine_archive)`;
the test's line 18). The legacy tree never fetched land cover or lakes.
Nothing is carried across.

## 3. Ideas file

Read: `../rasputin_scratch/ideas/ideas.md`. Overlapping, left out:

- "Repair CORINE once per area of interest": this increment caches the
  polygons as served, unrepaired; a repaired extract could later sit beside
  them in the same cache folder.
- "Structure for rasputin_data and rasputin_scratch": fetched data goes to
  the existing tile cache, not to a new folder layout.
- "Release hygiene ... DTM10 download instructions": a fetcher for
  Kartverket's DTM10 is a raster kind, not this one.

None is taken up.

## 4. The services, verified

Every figure below was measured on 2026-10-17 from this machine, anonymous,
with `curl`; the box is the bounding box of Numedalslagen
(`rasputin_data/numedalslagen_outline_nve.geojson`: 72,423 to 229,070 m east,
6,554,563 to 6,733,362 m north, EPSG:25833; about 157 km by 179 km, reaching
the coast at Larvik). It is the largest box measured.

### 4.1 NVE's lake database (Innsjødatabasen)

| Fact | Value | Where verified |
| --- | --- | --- |
| Service | `https://kart.nve.no/enterprise/rest/services/Innsjodatabase2/MapServer/5`, layer `Innsjodatabase`, `esriGeometryPolygon`, ArcGIS 12 | `?f=json` |
| Account or key | none | anonymous queries answered |
| Formats | `JSON, geoJSON, PBF` | `supportedQueryFormats` |
| CRS | EPSG:25833 (the layer's `extent.spatialReference.latestWkid`) | `?f=json` |
| Id field | `objectid` (`esriFieldTypeOID`) | `?f=json` |
| Page limit | `maxRecordCount` 2000 | `?f=json` |
| Box query | yes: an envelope with `inSR`; `returnIdsOnly` gave 26,538 ids for the box, 267,194 for the layer, no truncation flag | probes |
| Size for the box | 26,538 polygons; 127,234,169 bytes of GeoJSON in 14 requests of 2000, 44 s, asking `objectid,vatnlnr,navn,areal_km2,hoyde` | probe |
| `objectIds` by GET | 1000 ids (4,906-character URL) answered; 2000 ids (10,906 characters) dropped the connection | probe |
| Change date | none: the layer states no `editingInfo` | `?f=json` |
| Download of all Norway | `http://nedlasting.nve.no/gis`, Shapefile or SOSI, a form | Geonorge record below |
| Licence | "Norsk lisens for offentlige data (NLOD)", linked to `http://data.norge.no/nlod/no/1.0`; "Åpne data"; NVE "accepts no liability for damage or injuries due to errors in data and incorrect use of data"; scale 1:20,000; updated "Annually" | Geonorge record `823b8639-9a49-41bf-8571-3608435eb149` ("Hydrology - lakes"), WebFetch |
| Credit NVE asks | NVE names no credit line for this layer. For `api.nve.no` it asks: "Når du bruker data skal du så langt det lar seg gjøre lenke til den aktuelle tjenesten." | nve.no/om-nve/apne-data-og-api-fra-nve/, WebFetch |
| Credit NLOD asks when the licensor names none | «Inneholder data under Norsk lisens for offentlige data (NLOD) tilgjengeliggjort av [navnet på lisensgiver]»; "Er informasjonen endret skal lisenstaker tydelig angi at endringer er blitt gjort av lisenstaker" | data.norge.no/nlod/no/1.0, section 5 |

The service describes itself: "Innsjødatabasen til NVE inneholder
informasjon om vann og innsjøer. Til sammen inneholder basen ca 243 000
innsjøer.Alle innsjøene større enn 2500 m²/2,5 dekar har fått tildelt et
unikt nasjonalt løpenummer." (`serviceDescription`; the layer holds 267,194
features, more than the description's count.)

### 4.2 CORINE Land Cover 2018, vector

| Fact | Value | Where verified |
| --- | --- | --- |
| Product | "CORINE Land Cover 2018 (vector/raster 100 m), Europe, 6-yearly", release V2020_20u1, EPSG:3035; vector DOI `https://doi.org/10.2909/71c95a07-e296-44fc-b22b-415f42acfdf0` | land.copernicus.eu/en/products/corine-land-cover/clc2018, WebFetch |
| The product page's own web service | "REST API: image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CLC2018_WM/MapServer" (Web Mercator); its sibling `CLC2018_LAEA` serves the same layer in EPSG:3035 | product page; `Corine?f=json` lists both |
| Service used | `https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CLC2018_LAEA/MapServer/0`, "Corine Land Cover 2018 vector", `esriGeometryPolygon`, ArcGIS 11.5, `copyrightText` "European Environment Agency (EEA)" | `?f=json` |
| Account or key | none for the map service. The CLMS download API needs one: "Authentication is handled using EU Login" (land.copernicus.eu how-to guide "How to create API tokens to communicate with the CLMS APIs") | anonymous queries answered; guide fetched |
| Fields | `OBJECTID` (OID), `Code_18` (text, length 3), `ID`, `Remark`, `Area_Ha`, `Shape_Length`, `Shape_Area` | `?f=json` |
| Page limit | `maxRecordCount` 1000 | `?f=json` |
| Box query | yes: 8,297 ids for the box; 2,375,056 for the layer, no truncation flag | probes |
| Size for the box | 8,297 polygons; about 92 MB of GeoJSON in EPSG:3035: eight pages of 1000 gave 52,740,107 bytes, the last 297 about 39 MB, two of them sea polygons of 17 and 19 MB each | probe |
| Answer-size limit | a page holding both sea polygons (offset 8100, 200 records) is an HTTP 500 every time; 25-record pages holding one each are answered (17,092,659 and 19,082,916 bytes) | probe |
| `objectIds` by GET | 1000 ids (7,999-character URL) answered | probe |
| Change date | none: no `editingInfo` | `?f=json` |
| Data policy | "full, open and free access" under Commission Delegated Regulation (EU) No 1159/2013; users "shall inform the public of the source of these products and services"; where adapted or modified, "the user shall clearly state this"; and "shall make sure not to convey the impression to the public that the user's activities are officially endorsed by the European Union" | land.copernicus.eu/en/data-policy, fetched |
| Credit CLMS asks | for a dataset: "Dataset title. European Union's Copernicus Land Monitoring Service information , <provide the name and URL of data access point> (Accessed on DD.MM.YYYY) . <DOI: dataset DOI, if applicable>" | data policy page |
| EEA's terms for its own material | CC BY 4.0, "on condition that the EEA is acknowledged as the source"; nothing on map services or automated access | eea.europa.eu/en/legal-notice, WebFetch |

**To verify: the version the map service holds.** The same `OBJECTID`
differs between the service and Ola's v2020_20u1 file: 1880545 is `Code_18`
`112`, `ID` `EU-1880545` in the service and `324`, `EU_1880545` in
`U2018_CLC2018_V2020_20u1.gpkg`; the counts differ by 350 (2,375,056 and
2,375,406). The service states no version. So the catalogue entry names no
version and cites no DOI (the DOI is the V2020_20u1 release's), and the
shipped class maps' `notice`, which says "version 2020_20u1", is not reused
for this source (section 5.6).

## 5. The design

### 5.1 Data flow

```
rasputin fetch nve-lakes --bbox X0 Y0 X1 Y1 [--out-crs C] [--dry-run] [--refresh]
   │ FetchRequest (unchanged: box or domain, out_crs, dry_run, refresh)
   v
cli.fetch: kind == "arcgis-features" ──> fetch/features.py: fetch_features(request, source,
   │                                       get_text, cache)       (else fetch/run.py, unchanged)
   │  1. GET <url>?f=json                 layer facts: geometry type, id field, fields,
   │                                      maxRecordCount, CRS; checked against the entry
   │  2. box -> source CRS                transform_bounds (crs.py), no growth
   │  3. GET <url>/query?returnIdsOnly    the ids meeting the box
   │  4. ids not yet cached, in chunks    GET <url>/query?objectIds=...&outFields=<fields>
   │     (halved on refusal)              &outSR=<epsg>&f=geojson, one at a time
   v
FeatureCache (io/repository.py, under CacheWriter's lock)
   <cache>/<source>/features.gpkg   one table `features`, pk = the layer's id field,
                                    R-tree; each chunk one sqlite transaction
   <cache>/<source>/manifest.json   FeatureManifest: layer facts, requests
   <cache>/<source>/NOTICE.txt      notice(source, credits)
   │
   v
rasputin sources: a remote features entry is read from features.gpkg as a local
GeoPackage is (source_report._vector); path, fits, classes, coverage listed
   │  the agent passes the listed path
   v
--features PATH --features-map corine, catchment --lakes PATH (unchanged readers)
```

No second registry: the entry is a shipped `[source.<id>]` of the existing
catalogue; "registered" means the cache holds its file, which the listing
reads, exactly as a fetched raster source is listed through its cache today.
No user file is written.

### 5.2 The catalogue entry (`sources.py`)

`RemoteSource` gains a fifth kind and one field:

- `kind: Literal["local", "one-cog", "cog-tiles", "one-cog-per-year", "arcgis-features"]`.
- `fields: tuple[str, ...] = ()`: the attributes asked for and kept, an
  allow-list sent as `outFields` (never `*`, as `fetch/nve.py` rules for
  NVE's data); the id field is always kept, as the primary key, and is not
  listed.
- The `mode="before"` validator also gives an `arcgis-features` entry that
  states no `values` the value `"features"`: the kind says what it holds.
- `_consistent` for `arcgis-features`, each a `ValueError` naming the field:
  `url` needed; `crs` needed (as for every remote kind); `fields` needed and
  without repeats; `url_template`, the tile fields, `years`, `path`,
  `layer`, `nodata` and `z_unit` refused; `values` must be `"features"`.
  `fields` is refused for every other kind. `values = "features"` is
  allowed for `local` (as now) and `arcgis-features`.
- `_URL_FIELDS["arcgis-features"] = set()`; `url_for()` gives `url`.
- The module docstring and the `RemoteSource` docstring say the fifth kind.

The two shipped entries, in `data/sources.toml`:

```toml
[source.nve-lakes]
kind = "arcgis-features"
url = "https://kart.nve.no/enterprise/rest/services/Innsjodatabase2/MapServer/5"
crs = "EPSG:25833"
# The lake layer's fields asked; municipality, reservoir and globalid fields are never collected.
fields = ["vatnlnr", "navn", "areal_km2"]
credit = "Inneholder data under Norsk lisens for offentlige data (NLOD) tilgjengeliggjort av Norges vassdrags- og energidirektorat (NVE). Kilde: NVE, Innsjødatabasen, {url}, hentet {date}; lagret som GeoPackage av rasputin."
licence_note = "Norsk lisens for offentlige data (NLOD) 1.0, https://data.norge.no/nlod/no/1.0. NVE accepts no liability for damage or injuries due to errors in data and incorrect use of data."
description = "NVE's lake database: \"Innsjødatabasen til NVE inneholder informasjon om vann og innsjøer.\" Lakes larger than 2500 m² have a national number (vatnlnr); scale 1:20 000 (Geonorge)."

[source.clc2018]
kind = "arcgis-features"
url = "https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CLC2018_LAEA/MapServer/0"
crs = "EPSG:3035"
fields = ["ID", "Code_18"]
class_map = "corine"
credit = "CORINE Land Cover 2018 (vector/raster 100 m), Europe, 6-yearly. European Union's Copernicus Land Monitoring Service information, the European Environment Agency's map service {url} (accessed on {date}). Modified: the polygons meeting a box, stored as a GeoPackage by rasputin. Not endorsed by the European Union."
licence_note = "Copernicus Land Monitoring Service: full, open and free access (Commission Delegated Regulation (EU) No 1159/2013). Name the source, state any modification, and do not imply endorsement by the European Union."
description = "Land cover of Europe in 2018 in the 44 classes of the 3-level CORINE nomenclature (Code_18); \"minimum mapping unit (MMU) for status layers is 25 hectares\", \"minimum width of linear elements is 100 metres\" (the service's description)."
```

Defaults taken (Ola away; each follows existing behaviour): the NVE allow-list
is `fetch/nve.py`'s for layer 5 less `objectid`, which becomes the key; the
credit's date is ISO (`credit_for`), not CLMS's DD.MM.YYYY; the NVE credit is
NLOD's own default sentence in Norwegian, as the record links the Norwegian
1.0 text. The ids `nve-lakes` and `clc2018` are chosen names.

### 5.3 The fetch (`fetch/features.py`, new)

Network logic only, as `fetch/nve.py` is: `get_text` is injected
(`RangeClient().get_text`, which already retries timeouts and 5xx after 1, 2
and 4 s and names rasputin in `User-Agent`), the requests go out one at a
time, and URLs are built with `http.query_url`. No `urllib` import.

```python
class LayerFacts(BaseModel):          # frozen; from <url>?f=json
    geometry: Literal["esriGeometryPolygon"]
    id_field: str                     # the field of type esriFieldTypeOID
    fields: Mapping[str, str]         # name -> esri field type
    max_records: int                  # maxRecordCount
    epsg: int                         # extent.spatialReference.latestWkid, else wkid

def layer_facts(source: RemoteSource, get_text) -> LayerFacts: ...
# `fetch/plan.source_box` is split: its frame, domain bounds, `transform_bounds` and ±180°/pole refusal become
# `box_in_crs(request, crs) -> Bounds`. `source_box` adds the raster margin and two-cell growth to it, and the
# features fetch calls `box_in_crs` without growth.
def chunks(ids: Sequence[int], size: int) -> list[list[int]]: ...               # pure
def fetch_features(request, source, get_text, cache: FeatureCache) -> FeatureReport: ...
```

Rules, each a `FetchError` naming the source and URL:

- **The layer is checked against the entry** before any query: a geometry
  type other than `esriGeometryPolygon` is refused ("only polygon layers are
  fetched"; both shipped layers are polygons and nothing reads other kinds
  from a service yet); `epsg` must be `source.crs` (`same_crs`), as the COG
  header check does; every name in `source.fields` must be a layer field
  (named in the refusal); a layer with
  `uniqueIdInfo.OIDFieldContainsHashValue` true is refused (the API's
  caution above).
- **The box.** `request.box` in `out_crs` (else the source CRS), or the
  domain's bounds, moved into the source CRS with `crs.transform_bounds`
  (densified, outward). No growth: a feature meeting the box comes whole, so
  a mesh's reads (which take the features meeting its own box) find every
  feature of a domain inside the requested box. A box crossing ±180° or a
  pole, or with no image in the source CRS, is refused as
  `fetch/plan.source_box` refuses it.
- **Ids.** `query?where=1=1&geometry=<x0,y0,x1,y1>&geometryType=esriGeometryEnvelope&inSR=<epsg>&spatialRel=esriSpatialRelIntersects&returnIdsOnly=true&f=json`.
  An answer with `exceededTransferLimit` true is refused.
- **Chunks.** The ids not in the cache, ascending, in chunks of
  `min(max_records, MAX_IDS)`. `MAX_IDS = 500`, a chosen value: NVE dropped
  a 2000-id GET (10,906 characters) and answered 1000 (4,906); 500 keeps
  NVE's five-digit ids near 2,400 characters and CORINE's seven-digit ids
  near 4,000, under the 8,000-character URL CORINE answered. Checked at
  NVE's and CORINE's id lengths only; a layer with longer ids is not checked.
  Each chunk: `query?objectIds=<ids>&outFields=<fields>&outSR=<epsg>&returnGeometry=true&f=geojson`.
- **Halving.** A chunk whose request fails (a `FetchError` from `get_text`,
  which has already retried) is split in two and each half asked in turn; a
  chunk of one id that fails is the fetch's `FetchError`, naming the id. The
  measured case: CORINE's two sea polygons in one 200-record answer are a
  500 every time, each alone is answered. Cost: a failing chunk spends the
  client's 7 s of retries before it is split; an outage costs about nine
  levels of that before the fetch stops (500 halved to 1).
- **Every answer is checked.** The ids returned must be exactly the ids
  asked (a missing or extra one: "the layer changed during the fetch; run
  again"); the answer's `crs` member must name `outSR`; each geometry must be
  `Polygon` or `MultiPolygon` with finite coordinates. A feature with no
  geometry is refused, naming its id.
- **Geometry stored.** Each polygon as `MultiPolygon` (`shapely`), 2D, as
  served: no repair, no simplification, no clipping to the box. The
  GeoPackage column is declared `MULTIPOLYGON`.
- **Dry run.** Steps 1 to 3 only: the report gives the ids in the box, those
  already cached, and the requests the real run would send. Writes nothing.
- **`--refresh`** discards the source's cache folder first, as for rasters
  (`CacheWriter.discard`).
- **`--year`** is refused by the existing `url_for(year)` check ("has no
  years"); `--connections` does not apply (help text says "byte-range
  sources"), since requests go one at a time.

`FeatureReport` (frozen): `in_box`, `present`, `fetched`, `requests`,
`bytes` (answer bytes), `seconds`. The CLI prints one line in the raster's
form: `nve-lakes: 26,538 features in the box, 0 present, fetched 26,538 in
56 requests, <bytes> bytes, <seconds> s into <cache>/nve-lakes` (54 chunks
of at most 500 ids, plus the layer and id requests).

### 5.4 The cache (`io/repository.py`, `io/geopackage_write.py`)

- **`FeatureManifest`** (Pydantic, frozen), beside `CacheManifest`:
  `source`, `kind = "arcgis-features"`, `url`, `crs`, `id_field`, `fields`,
  `rasputin_version`, `first_fetched: date`, `requests:
  tuple[CachedRequest, ...]` (the existing model: region hash and date,
  hashed as `fetch/run._region` hashes it; that helper moves where both
  fetches reach it).
  `CacheManifest` is unchanged; `manifest.json`'s `kind` tells the two
  apart, and `CacheRepository` refuses a features manifest in plain words
  ("<source> holds features; pass its features.gpkg to a vector option").
- **`CacheWriter` gains `features(srs: SrsEntry, pk: str, columns:
  Sequence[tuple[str, str]]) -> sqlite3.Connection`**, inside `with writer:` only (the lock
  holds), opening `<dir>/features.gpkg`, creating it on first use with
  `create_geopackage` and a features table, then setting
  `journal_mode=DELETE` and `synchronous=NORMAL`: `create_geopackage`'s
  `journal_mode=OFF` suits a file written once, not a cache added to across
  runs, where a crash mid-chunk must leave the last committed chunk. Each
  chunk's features, R-tree rows and the `gpkg_contents` bounds update are one
  transaction, so a rerun resumes from the ids already present
  (`SELECT <id_field> FROM features`). `put_manifest` takes either manifest.
  `NOTICE.txt` through `put_notice(notice(source, (credit,)))`, the credit
  with the first fetch's date.
- The `CacheWriter` docstring's 'Every file is written as `<name>.part` and renamed' gains the exception: `features.gpkg` is written in place by SQLite transactions, whose journal guards a crash.
- **`io/geopackage_write.py` gains** `create_features_table(conn, name,
  srs_id, kind, pk, columns: Sequence[tuple[str, str]])` (columns
  `TEXT`, `REAL` or `INTEGER` from the esri field types; others refused),
  the blob code of `write_land_cover` (`shapely.bounds`, the `<2sBBi4d` header, ISO WKB) moves into `multipolygon_blobs(geometries, srs_id) -> (blobs, boxes)`, called by `write_land_cover` and `insert_features`. `build_rtree` is split into `create_rtree(conn, table)` (the virtual table and its `gpkg_extensions` row) and a fill that takes row ids. `create_features_table` calls `register_layer`;
  and `insert_features(conn, name, rows)` writing rows and their R-tree
  entries (id = the primary key). The R-tree table is created empty with
  the table (the `build_rtree` statement, without its bulk fill) and
  registered in `gpkg_extensions`. This module still opens nothing.

### 5.5 The listing (`source_report.py`)

`_one` for an `arcgis-features` entry: `where = <cache>/<id>/features.gpkg`,
`present` when `<cache>/<id>/manifest.json` is a file (the existing rule);
when present, `_vector` reads it as a local GeoPackage, so the entry gets
`path`, `geometry`, `features`, `attributes`, `classes` and `fits` exactly
as 57's local entries do. `_vector` takes the path as a parameter instead of
reading `entry.path`. The existing note "A remote source's figures describe
its cached objects, which may hold only the part fetched" covers it; its
coverage is the cached layer's extent box, as for any vector source. A
remote features entry not yet fetched lists no fits (nothing read), and
`present` false.

Fits that follow from 57's rules: `clc2018` fits `mesh --features` with
`class_map` `corine` (the entry's), whose attribute `Code_18` the table has;
both fit `catchment --lakes` (polygons). Neither fits the `station-catchments`
GeoJSON options: those read GeoJSON only, and their lakes come from
`fetch-stations`.

### 5.6 What a mesh made from it says

`--features PATH --features-map corine` carries the class map's `notice`,
which names "version 2020_20u1", the release of Ola's local file. For the
service's polygons that version is unverified (section 4.2). The design
makes no change to the class map: the cache's `NOTICE.txt` carries the
source's own credit, and the version question is for Ola (section 8).

### 5.7 MCP

`fetch` is already a tool derived from the CLI (`mcp/schema.py`,
`OPEN_WORLD`); nothing changes in the server. The instructions' sentence
"A remote source is read from the cache: fetch the area first" holds for the
new kind too.

## 6. Tests (for `@tester`)

No test reaches the network: `get_text` is a fake serving recorded answers
(small hand-made JSON shaped as the probes' answers), and the cache is
`tmp_path`. Fixtures are hand-made: three to six polygons, ids such as 3, 7,
11, 12; their figures are stated in each test.

1. **Catalogue.** Both shipped entries load; `arcgis-features` gets
   `values = "features"`; each refusal of section 5.2 by its field name
   (`fields` missing, `fields` on a `one-cog`, `layer`, `z_unit`, `path`,
   `years` on the new kind).
2. **Pure helpers.** `box_in_crs` with the same CRS (unchanged box) and
   with EPSG:4326 to EPSG:25833 (outward); `chunks` (7 ids, size 3 →
   3, 3, 1).
3. **Layer check.** Refusals: a line layer; a CRS mismatch; a field the
   layer lacks (named); `OIDFieldContainsHashValue` true.
4. **Fetch, happy path.** Four ids in the box, chunk size 2: two chunk
   requests whose `objectIds`, `outFields` and `outSR` are asserted from the
   URLs; the GeoPackage has four rows, keyed by id, `MULTIPOLYGON`, the
   allow-listed columns only, an R-tree whose extent equals the polygons';
   `manifest.json` is a features manifest; `NOTICE.txt` holds the credit
   with the date.
5. **Resume and overlap.** A second fetch of an overlapping box asks only
   the ids not cached (asserted on the URLs); a fetch interrupted by a
   `FetchError` on the second chunk leaves the first chunk's rows, and the
   rerun asks only the second.
6. **Halving.** A fake refusing any chunk holding ids 11 and 12 together:
   the fetch splits, both are stored; a fake refusing id 12 alone fails
   naming 12.
7. **Answer checks.** An id missing from an answer; an answer whose `crs`
   is not `outSR`; a feature with null geometry; `exceededTransferLimit` on
   the ids answer.
8. **Dry run** sends only the three metadata and id requests and writes
   nothing (the cache folder does not exist after).
9. **Listing.** After a fetch, `sources` lists the entry `present`, `path`
   the cached `features.gpkg`, `features` 4, `fits` with
   `mesh --features` (`corine`) for a `Code_18` fixture and `catchment
   --lakes`; before a fetch, `present` false and no fits.
10. **Through the readers.** The cached file passed to `catchment --lakes`'s
    reader (`read_lake_polygons`) and to `--features` with `corine` reads
    the fixture's polygons and codes: the cache is a GeoPackage the
    existing readers take.
11. **CLI.** `fetch clc2018 --bbox ... --dry-run` prints the plan line;
    `fetch nve-lakes --year 2020` is a usage error; a raster source's
    fetch is unchanged (the existing suite).

Not in CI: a live fetch. `@tester` may add one marked like the existing
data-dependent tests, off by default.

Invariant-critical suite: none (no topology or predicate); no mutation round.

## 7. Size and siting

| File | Net lines (estimate) |
| --- | --- |
| `fetch/features.py` (new) | 140 |
| `io/repository.py` (manifest, `features()`, refusal) | 45 |
| `io/geopackage_write.py` (table, blob, insert) | 35 |
| `sources.py` (kind, `fields`, rules) | 25 |
| `source_report.py` (`_one`, `_vector` path) | 15 |
| `cli.py` (branch, report line, help) | 30 |
| **Total** | **about 295**, one PR |

`data/sources.toml` is not counted (`CLAUDE.md` §2). No new dependency:
`sqlite3` and `json` are the standard library, `shapely` and `pyproj` are in
the core stack; `tools/check_prohibited_deps.py` stays green. No C++ is
touched, and no refine or mesh code, so no `@perf` acceptance run.

Assessment: data (the catalogue entry, `FetchRequest`) stays apart from
execution (`fetch_features`); the fetch holds no state beyond the cache it
is given; its network call is injected, so every rule is tested offline;
it is a blocking function like `fetch_station_set`, which an API worker runs
with `asyncio.to_thread`.

Stated limitations: no change detection (neither layer states an edit
date), so boxes fetched in different years may mix NVE's annual vintages;
the manifest's request dates say when each box was fetched, and
`--refresh` starts anew. A whole-country fetch is possible (267,194 lake
ids, 2,375,056 CORINE ids, Europe) but was not measured; by the box's
figures all Norway's lakes are about ten boxes, of the order of 1 GB of
answers (to verify).

## 8. Questions for Ola

1. **Licence and credit.** rasputin would download, for a box the user
   names, NVE's lake polygons (Norwegian open-data licence NLOD, which asks
   us to name NVE and the licence) and CORINE land cover from the European
   Environment Agency's public map service (free and open under EU rules,
   which ask us to name the source, say that we changed the data and not
   suggest EU endorsement). Neither needs an account. The data is never
   shipped with rasputin; each user's own fetch downloads it and writes the
   credit into a NOTICE.txt beside it. May rasputin ship these two download
   entries with the credit texts in section 5.2? Default: yes.
2. **CORINE version.** The EEA map service seems to hold a different
   release of CORINE 2018 than the v2020_20u1 file on your disk (one polygon
   checked has a different class, and the counts differ by 350 of 2.4
   million); the service does not say which. Should the fetched CORINE be
   used as it is, with no version named in its credit? Default: yes, and the
   CORINE class map's notice (which says "version 2020_20u1") stays as it
   is until the version is known.

## Review

- **Design review, round 1, 2026-10-17.** Commit `fe0259c4`, 0 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: `fetch/plan.source_box` split into `box_in_crs`; the `write_land_cover` blob code and `build_rtree` reused; the `CacheWriter` docstring exception; about 295 lines. Held: no second registry, no user file written; service facts (record limits, CRSs, fields); licence duties in NOTICE.txt. Idea: drop `ID` from the clc2018 fields. Open with Ola: shipping the two entries; the CORINE version.
