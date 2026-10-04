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

(pending)

## Questions for Ola

(pending)
