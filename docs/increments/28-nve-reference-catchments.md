# Increment 28 — NVE reference catchments: our catchments against NVE's, station by station

Status: **designed, awaiting review and Ola's answers** (`@architect`,
2026-10-04), branch `worktree-nve-catchments` off master `d20126b`. Every
choice Ola would normally make is marked "Default (@architect, 2026-10-04)"
and repeated, with its alternative, under "Questions for Ola".

**Closes.** Catchments for real Norwegian gauging stations, computed from the
DEM by `rasputin`, one polygon per station, each usable as `--domain`, and a
measured comparison with NVE's own catchment polygon for every station. After
this increment:

```sh
rasputin fetch-stations nve-hrd --out-dir ../rasputin_data/nve_hrd
rasputin catchments --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --stations ../rasputin_data/nve_hrd/stations.geojson \
    --reference ../rasputin_data/nve_hrd/reference.geojson \
    --out-dir ../rasputin_scratch/hrd_catchments
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --domain ../rasputin_scratch/hrd_catchments/2.11.0.geojson --tolerance 1 --out narsjo.vtk
```

It generalises increment 22's Bygdin acceptance (one catchment, area within
2 %, node overlap both ways) to 140 stations, and it adds what 22 left out:
flow accumulation, and snapping a gauge's coordinates onto the DEM's flow line.

**Not closed.** Discharge (only the keys to join it later are stored; see
"Discharge"). Lake gauges whose outlet lies farther than the snap radius from
the station point (a known limit, measured by the acceptance). Holes, as in 22.
Parallel stations (one station at a time; see "The batch"). Stream networks,
Strahler order. Stations outside DTM10's coverage, and any other DEM.

## Ola's ask, quoted

Ola, 2026-10-04: "I think we should also start making some actual
hydrological catchments for Norway soon. NVE should have a list of
unregulated catchments, where they also record the discharge." Approved the
same day as a new increment.

## What the data says (measured 2026-10-04)

Measured before the design, by throwaway scripts in
`../rasputin_scratch/nve28/` (not committed; the scratch folder is not a
channel, and these numbers are re-measured by the acceptance run).

### NVE's list: the Hydrological Reference Dataset (HRD)

- **The list exists and is the one Ola means.** NVE's streamflow part of the
  HRD, "Norwegian streamflow reference dataset for climate change studies"
  (the PDF at
  `https://www.nve.no/media/19932/norwegian-streamflow-reference-dataset-for-climate-change-studies-2021-aip.pdf`,
  fetched 2026-10-04, sha256 `62dc7741...7907`, "last updated in March 2026"),
  says: "All active streamflow time series stored in NVE's Hydra II database
  classified as active and unregulated constitute the basis", selected by six
  criteria (under 10 % of the area affected, no significant regulation, at
  least 20 years of record, active, good data, adequate metadata), and "140
  active and unregulated streamflow stations in Norway have been included in
  the 2025-version of the HRD". Its Table 1 gives, per station, the REGINE
  area number, the main number, the version of the discharge series, the name,
  and per-analysis usability flags. Parsing Table 1 (pdftotext, one regular
  expression) gives exactly 140 rows.
- **No coordinates or polygons in it.** The station number is
  `<regine>.<main>.0`: the HRD row "2 11 0 Narsjø" is NVE station `2.11.0`,
  and its third column is the discharge series version (parameter 1001,
  version 0), not part of the station number. Checked on Knappom (row
  "2 142 1"), which is station `2.142.0`, and on Narsjø, whose layer 14
  record (below) has discharge series versions 0, 1 and 2.
- **The wider set.** NVE's map service lists 862 stations with an active
  discharge series; 500 of them have regulation degree 0 for both area and
  reservoir (`reguleringsgradareal`, `reguleringsgradmagasin`), and 127 of the
  140 HRD stations are among those 500. The HRD adds the quality review; the
  500 do not have it. Default (@architect, 2026-10-04): **the 140 HRD
  stations**; the 500 are a later option (Question 1).

### NVE's polygons: "Totalnedbørfelt til målestasjon"

- **The service.** NVE's ArcGIS map service `HydrologiskeData3`,
  `https://kart.nve.no/enterprise/rest/services/HydrologiskeData3/MapServer`:
  layer 0 `Malestasjoner` (points: `stasjonnr`, `stasjonnavn`,
  `totalt_feltareal_km2`, `stasjonstatus`, ...), layer 14 `Vannforing_aktiv`
  (active discharge series: `stasjonnr`, `versjon`, regulation degrees), and
  layer 38 `Malest_totalnedb` (polygons: `stasjonnr`, `nedborfeltaareal_km2`,
  `oppdateringsdato`). All in EPSG:25833, the CRS of DTM10. The query
  endpoint returns GeoJSON in the asked CRS with
  `.../38/query?where=stasjonnr in ('2.11.0',...)&outFields=*&outSR=25833&f=geojson`,
  `maxRecordCount` 2000. The service description (Norwegian) calls the
  polygon "the catchment upstream of the gauging station".
- **All 140 found.** Layer 0 has every one of the 140 station numbers (all
  `stasjonstatus` 1, active), layer 38 has a polygon for every one: 130
  stations have one polygon, 10 have three, which differ by their update date
  and slightly in area (Narsjø: 119.74, 119.43, 119.74 km²). Rule: **the
  newest `oppdateringsdato` wins**, ties to the larger `objectid`; the chosen
  date is recorded.
- **Shape.** All 140 valid by shapely; one multipart (`19.79.0` Gravå, two
  parts), none with holes; median 983 vertices. Polygon area by shapely
  against the layer's `nedborfeltaareal_km2`: at most 0.005 km² apart. Against
  the station layer's `totalt_feltareal_km2`: from −2.8 % to +6.9 %
  (`11.4.0` the largest), so the two NVE areas are not the same number, and the
  polygon, not the attribute, is the reference here.
- **What the polygons are.** Not documented on the layer. NVE's catchment
  tool NEVINA (user guide, August 2025,
  `https://publikasjoner.nve.no/diverse/2025/Brukerveiledning.i.Nevina2025.pdf`)
  generates catchments "either fully automatically from the base data or a
  combination of automatically generated and catchments upstream of the point
  taken from REGINE", on a 20 m DEM since version 4; it warns that REGINE
  boundaries "can also be wrong". So the reference is NVE's best boundary,
  not ground truth, built on a coarser DEM than ours. Disagreements of a cell
  or two along the divide are expected; the acceptance is built for that
  (below).
- **Station point against polygon.** 98 stations lie inside their polygon,
  42 outside it, every one of the 42 within 233 m of it. Inside, the distance
  to the polygon's boundary has quartiles 20 m, 78 m, 351 m (max 4.4 km:
  lake gauges and stations in the middle of a wide valley floor). So the
  station point is near the outlet, but rarely on the DEM's flow line.

### Coverage by DTM10_UTM33_20260925

- **The DEM.** 254 tiles of 5051 × 5051 nodes at 10 m, EPSG:25833, spaced
  50 km apart so neighbours overlap by 51 nodes; together they span x
  −100 km to 1150 km, y 6400 km to 7950 km. NoData is −32767.
- **Which stations it covers.** Sampling each tile every 10th node (100 m)
  and each polygon at those samples: **138 of 140 polygons meet no NoData
  sample**. Two do: `244.2.0` Neiden (2,945 km², 2.6 % of its samples NoData)
  and `234.18.0` Polmak nye (14,171 km², 0.02 %). Both cross into Finland.
  The flood refuses a catchment that reaches NoData (increment 22, "The
  window", step a), so these two are expected to be refused, and are
  reported as such, not as failures. A 100 m sample can miss a NoData strip
  narrower than 100 m; the run itself is the check.
- **Tiles per catchment**, for the 138: 70 within one tile's footprint, 39
  in two, 8 in three, 18 in four (across a tile corner), 1 in five
  (`139.35.0` Trangen), 2 in nine (`212.10.0` Masi, `311.6.0` Nybergsund).
  The overlap strips count as both tiles, so "two" includes catchments that
  only reach into a 510 m overlap. 15a's mosaic handles all of these; the
  acceptance reports the result by tile count, so a seam artefact would show
  as misses gathered at multi-tile stations.
- **Sizes.** From 0.44 km² (`20.11.0` Tveitdalen, about 4,400 nodes) to
  14,171 km²; median 135 km²; 12 under 10 km², 44 from 10 to 100, 75 from 100
  to 1000, 9 over 1000. Sum 61,016 km², about 610 M catchment nodes at 10 m.

### Licence

- **NLOD.** Geonorge's metadata for "Totalnedbørfelt til målestasjon"
  (uuid `ac1c71db-9850-4e89-8162-2baba8b980e7`, fetched 2026-10-04) says
  "Åpne data" under "Norsk lisens for offentlige data (NLOD)" and adds that
  NVE disclaims liability for errors in the data and their use. NVE's HydAPI
  documentation says its data is under NLOD, "compatible with CC Navngivelse
  3.0 Norge (CC BY 3.0)". Increment 22's acceptance recorded the same terms
  for the reservoir layer, and NVE's request to credit it and link its
  services (`docs/benchmarks/2026-09-29/bygdin/README.md`, "Data and
  licence").
- **What is committed.** Default (@architect, 2026-10-04):
  - the list of 140 HRD rows (station number, series version, name, the HRD
    start year of daily data), as a data file in the package, credited to NVE
    in `NOTICE.md` and in the file's header: small, factual, and the one
    piece that cannot be fetched from a service (it lives in a PDF);
  - the per-station result table of the acceptance, which holds NVE's areas
    as numbers;
  - **not** the polygons (8.2 MB for all 160 versions) or the station points:
    `rasputin fetch-stations` fetches them, and records the fetch date, the
    URLs and each file's sha256 in a manifest, as 22's acceptance recorded its
    one polygon. Question 2.

### Discharge

NVE's HydAPI (`https://hydapi.nve.no/api/v1/`) serves the discharge series;
it needs a free API key (`/Stations` answers 401 without one) and is under
NLOD. **Out of scope now** (Default, @architect, 2026-10-04): no discharge
is fetched or stored. Each station keeps the two keys that join it later: the
station number (`2.11.0`, HydAPI's `StationId`) and the HRD's discharge series
(parameter 1001 and its version, `1001.0`). Question 3.

## Prior art: legacy and literature

**Legacy.** Nothing to carry over.

```
$ git grep -liE "nve|hydra|gauge|gauging|station|discharge|vannf|snap|accumulat" legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/avalanche.py
legacy-archive:legacy/rasputin/reader.py
legacy-archive:legacy/rasputin/solar_position.h
legacy-archive:legacy/rasputin/writer.py
```

Every hit but one is "convert" matching `nve`. `avalanche.py` calls NVE's
avalanche-forecast API, which has nothing to do with catchments or gauges.
Increment 22's grep for catchment terms still returns no files.

**Literature.** Each reference was checked against Crossref on 2026-10-04
(title, authors, venue, DOI), except where marked; abstracts or the pages
named were read, no paper in full.

- **Jenson 1991**, "Applications of hydrologic information automatically
  extracted from digital elevation models", *Hydrological Processes*
  5(1):31-44, doi:10.1002/hyp.3360050104. The "snap pour point" rule: move
  the gauge to the cell of highest flow accumulation within a radius.
  **The snap below builds on it.** What differs: the accumulation is the
  Priority-Flood's own drainage (each node drains to the node that flooded
  it, increment 22), not D8, so the snapped node's upstream set is exactly
  what `upstream()` labels from it; and the window is grown until every node
  in the radius has its whole upstream area inside it (below), which a
  fixed-extent accumulation grid does not guarantee at a window's edge.
- **Lindsay, Rothwell and Davies 2008**, "Mapping outlet points used for
  watershed delineation onto DEM-derived stream networks", *Water Resources
  Research* 44(8):W08442, doi:10.1029/2007WR006507. States the problem this
  increment meets ("outlet point positions taken from hydrometric stations
  commonly do not coincide with stream locations extracted from DEMs") and
  proposes AORA, which also uses the reported catchment area and stream
  names. **Departure, and why:** no reported area is used to place the seed
  by default, because the acceptance compares against NVE's area, and
  placing by it would make that comparison partly circular. The cost is the
  failure Jenson's rule is known for: in a radius that reaches a larger
  river below a confluence, the snap jumps to it. The acceptance measures how
  often (Question 4 offers the area-guided rule instead).
- **Lehner 2012**, GRDC Report 41, "Derivation of watershed boundaries for
  GRDC gauging stations based on the HydroSHEDS drainage network" (BfG,
  not on Crossref; the GRDC page
  `https://grdc.bafg.de/products/basin_layers/watershed_boundaries/` read):
  candidate outlets within 5 km, chosen by "the reported catchment area of
  the GRDC station and the distance between the station coordinates and the
  derived outlet", failures inspected by hand, at 15 arc-seconds (about
  500 m). **Färber et al. 2025**, "GRDC-Caravan: extending Caravan with data
  from the Global Runoff Data Centre", *ESSD* 17(9):4613-4625,
  doi:10.5194/essd-17-4613-2025, is where the GRDC page points for the
  current method (its section 2.3; not read here). The area-guided
  alternative of Question 4 is this rule.
- **Seppä, Gonzales Inca, Uusikivi and Alho 2026**, "CAMELS-FI:
  hydrometeorological time series and landscape properties for 320
  catchments in Finland", *ESSD* 18(7):4745-4769,
  doi:10.5194/essd-18-4745-2026. The closest prior work found: catchments for
  320 gauges delineated from Finland's 10 m DEM with WhiteboxTools, streams
  from an accumulation threshold, gauges snapped to them ("some gauges had to
  be moved slightly"), official river and basin data as guiding features,
  then every boundary inspected by eye and corrected with "virtual walls"
  until satisfactory. No agreement figure against official polygons is
  reported (as read from the article page). What differs here: no guiding
  data, no manual correction, and the agreement with the official polygon is
  measured and published per station.
- **Valseth, Valnes, Lappegard, Silantyeva and Mardal 2025**, "Development
  of CAMELS-Nordic, a large-scale hydrometeorological and catchment
  properties dataset for Norway and Sweden", EGU General Assembly 2025,
  EGU25-10411 (abstract read; not on Crossref as an article). It includes
  "catchment files"; the abstract does not say how they were made.
- **Liu, Zhang, Que, Tian, Yang and Hou 2025**, "Automating nested watershed
  delineation over the world considering endorheic basins and islands",
  *International Journal of Digital Earth* 18(1):2513044,
  doi:10.1080/17538947.2025.2513044 (search summary only): evaluated against
  4,301 GRDC stations over 1000 km², by relative area error alone (93 % below
  5 %). Area error alone is weaker than the overlap measured here: two
  catchments of equal area can share little.
- **The UK's NRFA** (catchment boundary page,
  `https://nrfa.ceh.ac.uk/data/about-data/catchment-information/catchment-boundary-and-areas`):
  boundaries from a 50 m terrain model, stations snapped automatically ("can
  lead to errors"), manual digitising where it fails, and the warning that
  under 100 km² "a small absolute difference in catchment area can cause
  large percentage differences". That warning is why the classes below have
  a size-independent second test.
- **Barnes, Lehman and Mulla 2014a** (doi:10.1016/j.cageo.2013.04.024) and
  **O'Callaghan and Mark 1984** (doi:10.1016/S0734-189X(84)80011-0), as cited
  in increment 22: the flood the accumulation reuses, and D8, which it
  departs from for the reason 22 gives.

**Novelty.** None is claimed. Searched: the six works above, Crossref and
web searches for automatic delineation of gauge catchments compared with
official polygons (overlap, Jaccard, area ratio), NVE's NEVINA guide, and
the HRD report. Not found in that search: a published, per-station,
two-way overlap comparison of a fully automatic 10 m delineation against
NVE's polygons for the HRD stations. That is a gap in what was searched, not
a claim; a write-up would search again (Scandinavian journals, NVE reports,
theses) before saying more.

## The design

### Data flow

```
rasputin fetch-stations nve-hrd --out-dir D                       [network]
  fetch/nve.py   reads the packaged HRD list (140 rows), queries NVE's
                 layers 0 and 38 through fetch/http.py, newest polygon per
                 station, writes D/stations.geojson, D/reference.geojson,
                 D/NOTICE.txt, D/manifest.json (URLs, date, sha256)

rasputin catchments --dem ... --stations D/stations.geojson
                    [--reference D/reference.geojson] --out-dir O   [offline]
  cli.py         paths stop here: io/station_set.py reads both files into
                 Station models and reference polygons; repository_for(dem)
  catchment_batch.run_batch(request, repository, stations, references, sink)
     for each station, in file order, one at a time:
       catchment.delineate(CatchmentRequest(seed=station, seed_crs=crs,
                                            snap_radius=R), repository)
         window loop seeded by the disc of radius R round the station  (22's loop)
         _core.accumulate(window)                                      [C++, new]
         snap: the disc node of most upstream nodes
         _core.upstream(window, that node) -> outline -> reduce         (22)
       reference.agreement(result, reference polygon)  [pure, shapely]
       reference.classify(agreement)                    [pure]
       sink.catchment(station, result)  -> O/<station>.geojson  (io/geojson.py)
       sink.row(StationResult)          -> O/results.csv
     reference.summarise(rows) -> O/summary.json, and stderr
```

No path, file, URL or CRS crosses into C++: `accumulate` sees the same
`RasterView` `upstream` does and returns an array (CLAUDE.md §2, I/O
boundary). No new dependency: the network is `urllib` inside
`fetch/http.py` (its existing `get_text`), polygons are shapely, the HRD list
is read with `csv` and `importlib.resources`. `tools/check_prohibited_deps.py`
covers the rest.

### Accumulation, from the same flood (C++)

`include/terrain/hydrology/accumulate.hpp`, header-only, next to
`upstream.hpp`, a template on `RasterSource`:

```cpp
namespace terrain::hydrology {
struct AccumulateOutcome {
    std::vector<std::uint32_t> count; // row-major; 0 on NoData, else the number of
                                      // nodes that drain through the node, itself included
};
template <raster::RasterSource R>
[[nodiscard]] AccumulateOutcome accumulate(const R& z);
}
```

- **The same drainage as `upstream`, by construction.** 22's flood labels a
  node *in* when it is a seed or was flooded from an *in* node; so each node
  drains to the node that flooded it (its "flooder"), and outlets drain out
  of the window. `accumulate` runs the same flood (same outlets, same key
  (level, push counter), same neighbour order) and records, per node, which
  of its 8 neighbours flooded it (one byte), and the order in which nodes
  were reached (the push order: a node is always reached after its flooder,
  so the reverse of that order visits every node before its flooder). One
  pass over the reverse order adds each node's count into its flooder's.
- **To keep the two floods from drifting**, the outlet set-up and the
  queue loop move out of `upstream.hpp` into `hydrology/flood.hpp`,
  `detail::flood(z, state, on_reach)`, where `on_reach(i, j)` is called when
  popped node `i` reaches node `j`; `upstream` labels in it, `accumulate`
  records in it. `upstream`'s behaviour and suite are unchanged.
- **The oracle is exact**: for every node `c` with data,
  `accumulate(z).count[c] == upstream(z, {c}).nodes_in`. Both are "the nodes
  whose chain of flooders passes through `c`". A test can check it node for
  node on small random DEMs with pits, flats and NoData, which is the
  invariant-critical suite of this increment (below). Also: the outlets'
  counts sum to the number of nodes with data.
- **Limits.** Refused with `std::length_error` when the raster has 2^32 nodes
  or more (the count and the order are 32-bit); the binding maps it to
  `ValueError`. Memory per node: 1 byte for the flooder, 4 for the order, 4
  for the count, beside the elevation array and the queue. Serial,
  deterministic; the binding releases the GIL and returns a numpy `uint32`
  array of the raster's shape that the result owns.

### The snap (Python, `catchment.py`)

`CatchmentRequest` gains `snap_radius: float | None = None`, metres in the
DEM's CRS. `None` keeps increment 22's behaviour exactly (no snap). A radius
with `lakes` is refused (`ValueError`: a lake is already the seed); a negative
or non-finite radius is refused.

With a radius `R`:

1. The station point moves into the DEM's CRS, as 22's seed does.
2. **The disc** is every DEM node at distance at most `R` from the point
   (closed; computed from node coordinates, not from a buffered polygon).
3. **The window loop is 22's, seeded by the whole disc**: the first window is
   the disc's bounds plus `WINDOW_MARGIN_M`, and it grows until the disc's
   catchment is clear of the window's edge, or refuses as 22 does (NoData,
   the data's edge, the memory cap). This is what makes the snap sound: the
   disc's catchment is the union of every disc node's upstream area, so once
   it is inside the window, every disc node's count in that window is its
   whole count, up to 22's known limit (a closed depression across the
   window's edge). A fixed window could cut the main river off at its edge
   and pick a side stream instead.
4. `_core.accumulate` over that window. **The snapped node is the disc node
   with data of the largest count**; ties go to the one nearest the station,
   then to the smaller (row, column). No disc node with data:
   `CatchmentError` ("no DEM node with data within R m of the station").
5. `_core.upstream` over the same window, seeded by the snapped node alone.
   Its catchment lies inside the disc's, so no further growth is needed. From
   here on it is 22's pour-point path: the outline is the ring around the
   snapped node, then the reduction.

The memory cap of 22's step 5 counts 9 more bytes per node when snapping.

`Catchment` gains `snap: Snap | None`:

```python
@dataclass(frozen=True, slots=True)
class Snap:
    station: tuple[float, float]  # the station point, in the DEM's CRS
    node: tuple[float, float]     # the snapped node
    distance: float               # metres between them
    radius: float                 # R
    upstream_nodes: int           # the snapped node's count
    nearest_upstream_nodes: int   # the count of the node nearest the station (no snap)
    disc_nodes: int               # disc nodes with data
```

**How a snap is reported.** `rasputin catchment --snap-radius METRES` (no
default: without it, 22's behaviour) prints one line on stderr: "moved the
start 143 m, from the station (x, y) to the DEM node (x, y) with the most
water within 250 m: 1,190,021 nodes (119.0 km²) drain through it; through
the node nearest the station, 312 (0.031 km²)". The GeoJSON's properties get
the same numbers. The batch writes them to every row.

**Default radius: 250 m.** Default (@architect, 2026-10-04). All 42
stations outside their NVE polygon lie within 233 m of it; a radius much
larger reaches bigger rivers more often. The acceptance re-runs every miss
at 100 m and 500 m, so the radius's part in each miss is measured rather than
guessed (Question 4).

**Known limits, measured not fixed.** (a) A disc that reaches a larger river
below a confluence snaps to it: the result is the larger river's catchment,
a miss with "ours in NVE's" low. (b) A lake gauge whose lake outlet is
farther than `R` from the station point snaps to a lake node near the
station, which drains only part of the lake's catchment: a miss with "NVE's
in ours" low. (c) A disc that reaches a large river makes the window loop
flood that river's catchment, to be thrown away after the snap: slow, and
refused if over the memory cap. Each shows in the acceptance's table by its
signature.

### The station set (`fetch/nve.py`, `io/station_set.py`)

**The packaged list.** `src_python/tin_engine/data/nve_hrd_2025.csv`: a
header naming the source PDF, its date and NVE's credit as comment lines,
then 140 rows `station,series_version,name,hrd_start_daily`
(`2.11.0,0,Narsjø,1931`). Data, not code; shipped in the wheel (the green step
checks `importlib.resources.files("tin_engine") / "data"` from an installed
wheel, since `wheel.packages` is the package directory). Written once by
hand-checked extraction from the PDF; `@tester` pins its row count (140), the
uniqueness of the station numbers, and the three spot rows quoted above.

**`rasputin fetch-stations nve-hrd --out-dir DIR [--refresh]`.** A catalogue
entry in `sources.py`, beside the DEM sources: `StationSource(id="nve-hrd",
service_url=..., list_file="nve_hrd_2025.csv", crs="EPSG:25833", credit=...,
licence_note=...)`. `fetch/nve.py`:

- queries layer 0 for the listed stations' points and attributes, and layer
  38 for their polygons, 40 station numbers per `where ... in (...)` query
  (measured to work; under the service's 2000-record cap and URL limits),
  with `outSR=25833&f=geojson`, through `RangeClient.get_text` (22a-2's
  retries and refusals; `fetch/http.py` stays the only module importing
  `urllib`, rule F11);
- refuses, naming the station, when a listed station has no point or no
  polygon; keeps the newest polygon per station (above) and records its
  update date and how many versions there were;
- writes `stations.geojson` (one Point feature per station, properties
  `station`, `name`, `series` = `["1001.0"]`, `nve_area_km2` from layer 0,
  `hrd_start_daily`), `reference.geojson` (one feature per station, the
  polygon, `station`, `reference_area_km2`, `reference_updated`,
  `versions`), both with a `crs` member, `NOTICE.txt` (the catalogue's
  credit, as 23a-2's `notice` does), and `manifest.json` (the query URLs,
  the fetch time in UTC, each file's sha256). Deterministic order: the list
  file's.
- Pure parts (the query URLs, choosing the newest version, building the
  files' content) are functions with no network, tested on canned answers
  as `fetch_fixtures.py` does; the network call is injected.

**`io/station_set.py`** reads both files back: `read_stations(path) ->
(tuple[Station, ...], crs)`, `read_references(path) -> (Mapping[str,
Polygon | MultiPolygon], crs)`. It refuses a file without a `crs` member
(the skill's rule; NVE's files always have one), duplicate station
numbers, and non-point or non-polygon geometry. `Station` is a frozen
Pydantic model: `station: str` (pattern `^\d+\.\d+\.\d+$`), `name`, `x`,
`y`, `series: tuple[str, ...]`, `nve_area_km2: float | None`. Any GeoJSON
of points with a `station` property works, so a user can bring their own
list.

### Agreement and classes (`reference.py`, pure)

Per station, against the reference polygon, on the DEM's node lattice (the
final window's origin and spacing, extended as far as either polygon
reaches; no DEM is read):

- `ours`: lattice nodes strictly inside our **fine** outline (the reduction
  is checked by 22's own guarantees, not here); `ref`: nodes strictly inside
  NVE's polygon (all parts); `both`. Counted with `shapely.contains_xy` in
  row bands, so memory stays bounded on Polmak-size polygons. This is the
  Bygdin comparison of 22's acceptance, unchanged.
- **NVE's in ours** = both / ref, **ours in NVE's** = both / ours.
- **Area ratio** = our fine area / NVE's polygon area (shapely).
- **Mean divide offset** = (ref + ours − 2·both) × cell area / NVE polygon's
  perimeter, in metres: how far, on average, our divide sits from NVE's.
  Size-independent, which is what NRFA's warning about small catchments asks
  for.

`classify` (Default, @architect, 2026-10-04; Question 5):

| class | rule | counts as |
|---|---|---|
| `match` | both overlaps ≥ 95 %, **or** mean divide offset ≤ 3 cells (30 m on DTM10) | pass |
| `close` | both overlaps ≥ 80 % | finding |
| `miss` | anything else | finding |
| `refused` | `delineate` refused (NoData, the data's edge, memory cap, no data in the disc), with its message | reported, not a failure |

Bygdin, the one case measured so far, would be `match` (99.12 % and 99.33 %).
The offset test is there for the 12 catchments under 10 km², where a one-cell
disagreement along the whole divide is several per cent of the nodes.

`summarise(rows) -> Summary`: counts per class; for the area ratio and both
overlaps, the minimum, 10th, 25th, 50th, 75th, 90th percentile and maximum,
over all accepted stations and per size band (under 10, 10-100, 100-1000,
over 1000 km²) and per tile count (1, 2, 3-4, 5+). Deterministic JSON.

### The batch (`catchment_batch.py`)

```python
class BatchRequest(BaseModel):      # frozen
    snap_radius: float = 250.0
    outline_tolerance: float | None = None
    only: tuple[str, ...] = ()      # station numbers; empty is all

class BatchSink(Protocol):
    def catchment(self, station: Station, result: Catchment) -> None: ...
    def row(self, row: StationResult) -> None: ...

async def run_batch(request: BatchRequest, repository: DemRepository,
                    stations: Sequence[Station], stations_crs: str,
                    references: Mapping[str, BaseGeometry] | None,
                    sink: BatchSink) -> Summary
```

- **One station at a time**, each `delineate` in `asyncio.to_thread`: a
  flood of a large catchment takes gigabytes, and two at once would race for
  the memory cap. Sequential order is the file's, so the output is
  deterministic. A GUI or API worker awaits it and can cancel between
  stations. Parallel stations are a later option.
- A refusal (`CatchmentError`) becomes a `refused` row and the batch goes
  on; any other exception stops it (a bug is not a data refusal).
- `StationResult` (frozen): station, name, class, refusal message, snap
  numbers, nodes, fine and reduced area, NVE's polygon area and the station
  layer's area, the agreement numbers, tile count, windows, seconds.
- Without `--reference`, no agreement and no class: the batch just makes the
  catchments (the product for any list of stations).
- No paths: the sink is how files get written; `cli.py` passes a directory
  sink, tests pass a list.

**`rasputin catchments`**:

```
rasputin catchments --dem PATH [--dem PATH ...] --stations FILE
                    [--reference FILE] [--snap-radius METRES] [--only ID ...]
                    [--outline-tolerance METRES] --out-dir DIR [--out-parent DIR]
```

writes `DIR/<station>.geojson` (22's catchment file, plus the station's
number, name, series and the snap; `mesh --domain` reads it), `DIR/results.csv`,
`DIR/summary.json`, and one stderr line per station (number, name, class,
the three numbers) and the summary at the end.

**The GeoJSON writer moves** from `cli.py` to `io/geojson.py`,
`catchment_geojson(polygon, crs, properties) -> bytes`, no path, as
`project_structure.md` already recommends ("The catchment GeoJSON writer is
in `cli.py`"); both commands call it, and that paragraph is replaced by the
module's entry.

### New and changed files

| File | What | Estimate |
|---|---|---|
| `include/terrain/hydrology/flood.hpp` | the shared flood, moved out of `upstream.hpp` | 45 (moved, ~10 net) |
| `include/terrain/hydrology/upstream.hpp` | uses it | −35 |
| `include/terrain/hydrology/accumulate.hpp` | `accumulate`, `AccumulateOutcome` | 45 |
| `bindings/core.cpp`, `_core.pyi` | `accumulate` | 30 |
| `catchment.py` | `snap_radius`, disc, `Snap`, pick, re-flood, cap | 75 |
| `cli.py` | `catchment --snap-radius` and its line | 20 |
| **PR 1** | | **about 190** |
| `data/nve_hrd_2025.csv` | the list (data, not counted) | 0 |
| `sources.py` | `StationSource`, the `nve-hrd` entry | 30 |
| `fetch/nve.py` | queries, newest version, files, manifest | 120 |
| `io/station_set.py` | readers, `Station` | 50 |
| `io/geojson.py` | the moved writer | 25 (cli.py −25) |
| `reference.py` | agreement, classes, summary | 100 |
| `catchment_batch.py` | `BatchRequest`, `BatchSink`, `run_batch`, `StationResult` | 80 |
| `cli.py` | `fetch-stations`, `catchments`, the directory sink | 100 |
| `NOTICE.md`, `project_structure.md` | NVE's credit; the new modules | docs |
| **PR 2** | | **about 480** |

Both under the 700 ceiling of CLAUDE.md §2. PR 1 touches C++ (one build
round), PR 2 is Python only, so PR 2's red step can be written while PR 1
is in review.

### The PR split

- **PR 1, the snap.** Accumulation and the snap, behind
  `rasputin catchment --snap-radius`. Answers on its own: "the catchment of
  this gauge". Red, green, review.
- **PR 2, the stations and the batch.** Fetch, read, compare, run. Red,
  green, review, then the acceptance run below.

Neither touches refine or mesh code, so `tools/bench.py`'s benchmark and
scaling sweep do not apply (`docs/increments/README.md`, "Acceptance").

## The red suites

Lean, as 22's were: no throwaway implementation. **The invariant-critical
suite** (mutation testing required, `docs/increments/README.md`, "Cost
constraints") is `accumulate`'s exact oracle, because every snap and
therefore every station result rests on it.

**PR 1, C++ (`tests/cpp/unit/test_hydrology_accumulate.cpp`):**

- *The oracle*: on a few hundred small random DEMs (with pits, flats, equal
  heights, NoData holes and NoData borders), for every node with data,
  `count == upstream(z, {node}).nodes_in`; NoData nodes count 0.
- *Conservation*: the outlets' counts sum to the nodes with data; every
  count is at least 1 on data.
- *By hand*: a V-valley (the outlet counts the whole valley); a single
  cell; a flat plateau draining over one rim node; a 1 × n strip.
- *Determinism*: twice gives equal arrays. *Refusal*: the size check, by a
  geometry stub reporting 2^32 nodes (no allocation).
- `upstream`'s existing suite, unchanged, after the move to `flood.hpp`.

**PR 1, Python:**

- `test_core_accumulate.py`: shape, dtype `uint32`, ownership (the array
  outlives the view), the oracle on two DEMs through the binding.
- `test_catchment.py` gains: a station 3 cells beside a synthetic river
  snaps to the river node of largest count, and the catchment equals one
  flood from that node; a tie goes to the nearer node, then the smaller
  (row, column); a disc that reaches past the window's first extent grows the
  window (the river's upstream reaches far beyond the disc), and the result
  equals a whole-raster run; a disc all NoData is refused; `snap_radius`
  with lakes, negative, NaN are refused; radius 0 picks the nearest node and
  reports distance as its offset; no radius gives 22's result bit for bit.
- `test_cli_catchment.py` gains: `--snap-radius` prints the snap line and
  writes the snap properties.

**PR 2, Python** (no network anywhere):

- `test_fetch_nve.py`: the packaged list (140 rows, unique numbers, the
  three spot rows); the query URLs (chunks of 40, `outSR=25833`,
  `f=geojson`); newest version wins, ties to the larger `objectid`; a
  missing station or polygon is refused by name; the written files, from
  canned service answers, read back through `io/station_set.py`; the
  manifest's sha256 matches the files.
- `test_station_set.py`: missing `crs`, duplicate numbers, wrong geometry
  types, a bad station number are refused; a user's own points file reads.
- `test_reference.py`: agreement on hand-built polygons on a 10 m lattice
  (identical: 100 %, offset 0; shifted by one cell: the offset is one cell;
  disjoint: 0 %); the class boundaries at exactly 95 %, 80 % and 30 m; the
  summary's percentiles on a known list; band and tile grouping.
- `test_catchment_batch.py`, on the synthetic tiled DEM of
  `test_cli_catchment.py` with three stations (one matching a reference
  drawn from its own flood, one with a reference shifted to make it a miss,
  one on NoData, refused): the rows, the classes, the order, the summary;
  a bug-type exception stops the batch.
- `test_cli_catchments.py`: the files in `--out-dir`, `--only`, the stderr
  lines; `--reference` absent gives catchments and no classes; a station
  file without `crs` is refused.

## Acceptance: every covered HRD station

Run after PR 2 is green, by `@perf` (it measures, and owns the evidence
layout), under `docs/benchmarks/<date>/nve-hrd/` with a `run.sh`, the
commit, `pmset -g batt`, the manifest of the fetched station set (its sha256
values, since the service can change) and the outputs that are not NVE's
data:

1. `rasputin fetch-stations nve-hrd` into `../rasputin_data/nve_hrd`.
2. `rasputin catchments` over all 140 at the default radius (250 m) and
   tolerance; wall time and peak memory per station and in total.
   Expected (estimate, not a measurement): the catchments total about 610 M
   nodes, windows about three times that, so tens of minutes on the M1 Max.
3. **The table**: `results.csv` committed, and the summary by class, size
   band and tile count in the README.
4. **Every finding explained**: each `miss` gets a line (snap jumped to
   another river, lake gauge, DEM artefact, NVE polygon disagrees with the
   DEM, or "not explained"), with a re-run at 100 m and 500 m; `close` rows
   are summarised by cause. The two Finnish-border stations are expected
   `refused` (NoData); a refusal for any other reason is a finding.
5. **What passes the increment**: the batch completes for every station with
   a row each; 22's guarantees hold on every accepted reduced outline (area
   kept, simple, start inside); every finding has its line. No share of
   `match` is required this time: this run is the baseline the next
   increments improve (Question 5).
6. Bygdin is not in the HRD (it is regulated); 22's run stays its record.

## What each persona reads

`@tester` and `@developer`: this file, then `docs/increments/22-auto-catchment.md`
("The seed", "Flow and membership", "The window"), and for PR 2
`docs/increments/23-basin-scale.md` ("The fetch step and the tile cache")
for `fetch/`'s rules. `@developer` also reads `include/terrain/hydrology/upstream.hpp`
and `src_python/tin_engine/fetch/http.py`.

## Questions for Ola

(pending)
