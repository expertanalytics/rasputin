# Increment 29 — NVE reference catchments: our catchments against NVE's, station by station

Status: **design approved by `@reviewer` (round 5, 2026-10-04); PR 1 merged as #173; PR 3 merged as #178; PR 2 (the gauge on the river): green `193079d`, 569 production lines; code review round 3 (2026-10-05) approved; `@perf`'s placement figures next, then the push; questions 5 and 6 open for Ola, written to their defaults**
(`@architect`, 2026-10-04), branch `worktree-nve-catchments` off master
`d20126b`. Ola's rulings of 2026-10-04 are in the section below. Round 2 closed the burn's drainage claim
(checked node by node, not assumed), the ELVIS data cases, the PR order, and
added "Data use". Round 3 makes the burnt chain taut (no two of its nodes
that are not consecutive may be neighbours), corrects the lowering bound and
the tile sizes, and records Ola's ruling on the nine shifted-tile stations.
Round 4 states the end extension's stopping rule, qualifies what the burn
leaves unchanged, words the mixed-grid refusal for every cause the mosaic
has, re-measures `7707_1`'s overlap against all four neighbours and
`7707_3`'s against all seven, has the figure script's survey only propose
cases that the full run must confirm, and records Ola's ruling on the figure
script. Round 5 approved the design, and Ola ruled the last three questions
the same day. The gauge
placement is redesigned to his direction (the mapped river, then the DEM's
flow path along it, never an area objective), with a per-station sensitivity
check. Five PRs (merge order 1, 3, 2, 4, 5). Questions 1 to 4 under "Questions for Ola" are
ruled (see "Ola's rulings"); questions 5 and 6, from PR 2's code review,
are open, and the design is written to their defaults.

**Closes.** Catchments for real Norwegian gauging stations, computed from the
DEM by `rasputin`, one polygon per station, each usable as `--domain`, and a
measured comparison with NVE's own catchment polygon for every station. After
this increment:

```sh
rasputin fetch-stations nve-hrd --out-dir ../rasputin_data/nve_hrd
rasputin station-catchments --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --stations ../rasputin_data/nve_hrd/stations.geojson \
    --rivers ../rasputin_data/nve_hrd/rivers.geojson \
    --reference ../rasputin_data/nve_hrd/reference.geojson \
    --out-dir ../rasputin_scratch/hrd_catchments
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --domain ../rasputin_scratch/hrd_catchments/2.11.0.geojson --tolerance 1 --out narsjo.vtk
```

It generalises increment 22's Bygdin acceptance (one catchment, area within
2 %, node overlap both ways) to 140 stations. It also adds what 22 left out:
flow accumulation; placing a gauge where it physically is (on NVE's mapped
river, then on the DEM's flow path along it); and a per-station check of
whether the catchment area is well defined at the gauge at all.

**Not closed.** Discharge (only the keys to join it later are stored; see
"Discharge"). Residual inflow between gauges on one river (the interfaces are
sketched so that this increment does not block it; see "Residual inflow,
later"). Burning the whole mapped network into the DEM (only the gauge's own
reach is burnt). The nine stations whose catchments straddle DTM10's
half-cell-shifted tiles (expected refusals, reported apart as known
refusals, not failures; a later increment fixes the tiles) and
Femundsenden (no mapped river within 500 m; refused until the fallback).
Holes, as in 22. Parallel stations (one station at a time; see "The batch").
Stations outside DTM10's coverage, and any other DEM.

**Merge note.** `ROADMAP.md:54` holds this increment's row (29) on this
branch and increment 28's row (no fused multiply-add contraction) on branch
`worktree-fpc`, on the same line of the same base. Whichever merges second
gets a textual conflict there; the resolution keeps both rows, 28 above 29.

## Ola's ask, quoted

Ola, 2026-10-04: "I think we should also start making some actual
hydrological catchments for Norway soon. NVE should have a list of
unregulated catchments, where they also record the discharge." Approved the
same day as a new increment.

Ola, the same day, on the gauge placement: "At some point, it would be great
to see how the gauge placement algorithm works. Preferably a few examples
with maps "before" and "after", as pngs." Taken as a deliverable of PR 2
("Placement figures", below).

## Ola's rulings (2026-10-04)

- **Questions 1, 2, 3, 5 and 6** (numbered as in the first version of this
  file; "Questions for Ola" below is numbered afresh): the defaults stand ("defaults on all six",
  Ola, for these five). Question 1's text is corrected below (364 + 136, not
  "500 unregulated"); the default of the 140 is unchanged. Question 5 gains
  the class `uncertain`, which Ola's direction on the gauge asks for.
- **Question 4, the gauge placement, reopened.** Ola: "On the gauge, we
  actually need to be very careful." and "Some of the points of interest
  will be saddle point problems, right? One of my plans is to also compute
  the _residual inflow_ to rivers, and then moving a point 100m downstream
  will always give you more inflow. So that rule is broken." He agreed ("This
  matches well") to the direction this revision follows. The gauge is placed
  where it physically is, never by an area or flow-count objective. A
  sensitivity check per station marks ill-posed stations as uncertain.
  Residual inflow is later work, and nothing here may block it.

**Ruled later the same day (2026-10-04).**

- **The nine stations on the half-cell-shifted tiles** (Question 2 under
  "Questions for Ola", now closed). Ola: "Yes, keep them as failures - this
  is a success story!" They stay **expected refusals** in this increment, and
  a later increment fixes the eight tiles (a resample; see "The
  half-cell-shifted tiles").
- **How they are counted.** Ola: "Exactly, the nine are not failues on their
  own, they are not available for us due to other limitations." They are
  **reported apart, as known refusals, not as failures**: each of the nine
  rows names the mixed grids as the cause (`refusal_cause = "mixed_grid"`
  and one tile from each grid, "The batch"), and the summary has one line for them;
  the acceptance README's line adds that a later increment fixes the tiles.
- **The placement-figure script** ("Placement figures"). Ola, 2026-10-04
  (10:51 UTC), answering two questions from the main session: "map questions:
  1: No, 2: Yes." Question 1: does the 130-line script count toward PR 2's
  700-line limit? No. Question 2: may it plot with matplotlib, in that
  evidence script only and in no dependency list? Yes.
- **Questions 1, 3 and 4 under "Questions for Ola"**, now closed. Ola,
  2026-10-04 (11:52 UTC): "Yes to all three.", answering the main session's
  restatement of the three with their defaults:
  - **Question 1, when a station is too uncertain to score:** the area is
    tested up and down the river as far as the gauge's coordinates are from
    the river, at least 30 m, against a 5 % bar (`SWING_MAX = 0.05`).
  - **Question 3, how much river to burn:** only the gauge's own stretch now;
    burning the whole network is a later increment.
  - **Question 4, the nearest-stream fallback:** included, as the last PR
    (PR 5).

**Ruled for PR 3's green step (2026-10-04, 18:19 UTC).** Ola: "defaults on
both", answering two questions `@tester` raised with PR 3's red step
(`213a2ef`); the third point below is `@architect`'s, not a question to Ola.

- **The station's name comes from NVE's station layer** (layer 0,
  `stasjonnavn`), copied as served, not from the packaged list. The list's
  `name` column is extracted from the PDF's Table 1, where long names wrap
  over two lines; the list stays the source of which stations, their series
  version and `hrd_start_daily`, and its names only serve to check the
  extraction. Measured 2026-10-04 against layer 0 for all 140, with the
  first-pass extraction of the design rounds: all 140 have a non-blank name in
  layer 0 (no feature in the whole layer has a null or blank one), and 12
  differ from the extracted names, each because the PDF wraps or splits it
  (`311.4.0`: "(Femunden)" extracted, "Femundsenden (Femunden)" in layer 0;
  `27.15.0`: "t)" against "Austrumdal (Austrumdalsvatnet)"; `234.13.0`:
  "Iesjokka" against "Veahkkava, Iesjokka"). **Red amendment:** the fake
  service gives one station a layer 0 name that differs from its list row,
  and `stations.geojson` must carry layer 0's.
- **A river segment served as a MultiLineString is refused**, naming its
  `objectid`, as "The station set" already says (one LineString per
  segment; "the wrong geometry type" is refused). Not split, not merged: a
  segment is one link of one river, and the chain of "Placing the gauge"
  joins segments end to start. The refusal is `read_segments`'s; the fetch
  writes the feature as served, so a user's file and a fetched one are
  refused alike. **Red amendment, before green:** `test_rivers.py` adds a
  MultiLineString segment among good ones and expects a `ValueError` naming
  its `objectid` and "MultiLineString".
- **Where the count of dropped copies goes** (`@tester`'s pin 15, changed).
  `io.rivers.drop_copies(segments) -> (kept, dropped)` stays as pinned: pure,
  public, applied by `read_segments`. But a two-part `read_segments` loses
  the count at the reader, and nothing downstream can recover it, since
  `drop_copies` on the reader's output finds nothing left. So
  **`read_segments(path) -> (segments, crs, dropped)`**, a three-part tuple
  (not a field on `RiverSegment`: the count describes the file, not a
  segment). It reaches Ola as one stderr line from each command that reads a
  river file, after reading it: `rasputin catchment --rivers` (PR 2) and
  `rasputin station-catchments` (PR 4), worded, for example, "rivers: 4,812 segments read,
  25 exact copies dropped (same river, same vertices to 1 cm)" (segments read
counts the file's features, copies included: kept plus dropped); PR 4 also
  writes it to `summary.json` as `river_copies_dropped`. `fetch-stations`
  writes copies as served and reports nothing about them. **Red amendment:**
  the two places that unpack `read_segments` (`test_rivers.py`'s `read` and
  `test_fetch_nve.py`'s read-back) take three parts, and
  `test_the_count_dropped` also asserts the reader's count is 25.

**PR 3, after code review round 1 (2026-10-04): the lake number's 0.** Not
a ruling: a correction of the data, found by `@reviewer` and re-measured by
`@architect` ("The data are not clean"). NVE's river service sends
`vatnlnr` = 0 for "no lake", so in `kind`'s mapping (see "Placing the
gauge") "set" means **not null and not 0**. The design said "set", and the
green step's `io.rivers.kind_of` reads it as "not null", so a blank-type
river with 0 (both blank-type features of the layer) comes out `lake`.

- **`@tester`, red first:** the fake service sends 0 as NVE does: the
  blank-type segment of `tests/python/nve_fixtures.py` (`BLANK_TYPE`, today
  `vatnlnr=None`) gets `vatnlnr=0`; `test_rivers.py` adds that `vatnlnr` = 0
  with a blank type and with a null type each give `river`, and keeps a
  positive number with a null type giving `lake`.
- **`@developer`, green, in `io/rivers.py`'s `kind_of`:** null or blank type
  is `lake` only when `vatnlnr` is not null and not 0. With it, the
  reviewer's two suggestions on error wording: a malformed user file (station
  or river) is refused with a `ValueError` naming what is wrong, not a
  `TypeError` or `KeyError` escaping from the reader; and `fetch-stations`
  reports a missing field in plain words rather than a raw `KeyError` text
  (`cli.py` catches `KeyError` with `FetchError` and `ValueError` and prints
  it as is).

**PR 3, after PR #178's CI (2026-10-04): the suite does not read `NOTICE.md`.**
Ruled by Ola: "yes, remove the NOTICE.md check." h11's prose-read hook failed
#178's three Python legs (3.12, 3.13, 3.14) because `test_fetch_nve.py`'s
`test_notice_md_credits_nve_under_nlod` read `NOTICE.md`. The main session put
it to Ola that a unit test checking the repository's prose is the wrong place
for it: the suites check what rasputin produces (`NOTICE.txt` of every fetch,
the packaged CSV's header), and `NOTICE.md` stays prose, credited by hand and
seen by `@reviewer` on any PR that touches it. `@tester` removes that one test;
`NOT_PROSE` in `tools/ci_changes.py` is unchanged. NVE's credit in `NOTICE.md`
(its section "Data in the package"; this file's §"Licence" and §"Data use")
stands.

**PR 2's red step: three design questions (2026-10-05).** Not Ola's
rulings: `@architect`'s answers to `@tester`'s questions on red `3a06ffc`,
written into the design where they apply. None needs Ola; the first follows
his ruling on the gauge (Question 4 above).

1. **The valley floor biased the gauge downstream.** "The least elevation
   within 30 m" picks, on a falling floor, the node about 30 m downstream of
   each point, so the placed node too: a downstream bias, which Ola's ruling
   forbids. **Changed**: the least elevation is taken across the line, in a
   cross-section one step wide ("Following the river", step 2, with the
   measurement). The tests that pinned the old rule's nodes change:
   `test_burn.py`'s off-floor, tie and placed-node tests (chains on rows 5
   to 40, not 6 to 41 or 7 to 42; placed node `(15, 20)`, `node_offset_m`
   28.0), and `test_cli_catchment.py`'s placed node (row 150, not 152),
   offset, swing rows and stderr sentence ("moved ... onto the DEM's valley
   floor", "Following the river" and `rasputin catchment` above). The
   embankment tests, and every test at a 5 m corridor with its points on
   nodes, are unchanged.
2. **`U` between nodes.** Read to `U` means every sample within `U` trusted
   and a chain node at or past `U` ("Sensitivity", step 2; `D` in "The
   window loop"), not `checked_down_m >= U`. `@tester`'s reading, adopted.
3. **The corridor's boundary is inside** (`<=`, with 1e-6 m of slack for
   rounding; step 2).

Also pinned there: the prefix tier ("Placing the gauge", step 2); the count
of segments read (the file's features, copies included); a station in
another CRS is transformed into the river file's CRS, which must be the
DEM's; `GaugePath.placed` is one `int` in 29 ("Residual inflow, later").

**PR 2's green step: `@developer`'s choices beyond the design
(2026-10-05).** `@architect`'s reading of green `91857df`, not Ola's
rulings. Adopted and written where they apply: the extension's tie (a
straight step first, "Following the river", step 5); NoData never a
candidate, in the cross-section or the extension (steps 2 and 5); the
resample's end and `dropped_m` (step 1); a tenth of the chain (step 3);
equal or falling counts are `chain_not_draining` ("Sensitivity", step 3);
the reach's CRS, stage B's repeated window, `GaugeResult` without
`reach_fork`, the `direction` cause joined by the command, and the memory
figure ("The window loop"). **One is changed, provisionally**: an upstream
flag within `U` shortened `checked_up_m` and nothing else; the design's
"counts as unread" now means `drains` false and the cause
`chain_not_draining` ("Sensitivity", step 2, with the proof). It is
`@architect`'s reading, not a ruling: question 6 for Ola, default yes. **Red amendment** (`@tester`):
`test_sensitivity.py`'s `test_a_flag_upstream_stops_the_read_there` and
`test_the_read_stops_at_the_first_flag_not_at_a_position` also assert
`drains` false, `"chain_not_draining"` in `causes` and `well_posed` false;
`test_a_flag_past_u_does_not_shorten_the_read` (a flag at -40 m with `U` =
30 m, well posed) stays as it is and must stay green. **Green** (`@developer`, `sensitivity.py`):
`drains` is also false when the upstream read stopped at an untrusted
sample whose arc length is within `-U` (a few lines).

**PR 2's code review, round 1 (2026-10-05): a NoData node on the chain.**
`@reviewer` found that the chain can pass through NoData ("Following the
river", step 2, "No NoData on the chain"); the design now refuses such a
station, provisionally (question 5 for Ola, default refuse). **Red**
(`@tester`, lean): in `test_burn.py`, one test parametrized over two float32
sentinels, −32767 and 3.4e38: a floor that runs down column 20 above row 20
and down column 22 from row 20 on, `z[20, 21]` NoData, the mapped line down
column 20, corridor 30 m. The test first asserts, from the fixture alone,
that `(20, 21)` is the straight join's middle node between the floor nodes
`(19, 20)` and `(20, 22)`, so it reaches the join-through-NoData path; then
`burn_reach` raises `ValueError` and the message names a gap in the DEM.
The same DEM with data at `(20, 21)` burns without error (the control: the
check is on the chain, not on NoData anywhere in the window). In
`test_catchment.py`, one test: the same kind of gap **downstream** of the
placed node (so that stage A's own NoData and seed refusals cannot fire
first), and `delineate` with that reach raises `CatchmentError` with the
gap message. Each must fail on `0588bfa` by not raising, not by another
refusal; say so in the handback. **Green** (`@developer`): in `burn.py`,
after `_taut` and before the direction check, a taut-chain node without
data raises `ValueError` with the plain message of step 2 and the node's
position in the DEM's CRS; in `catchment.py`, the burn's `ValueError`
reaches the caller as a `CatchmentError` with the same words (as `_plan`
turns a `MosaicError` into one); in `sensitivity.py`, the comment at lines
71-72 says what "Sensitivity", step 2, now says about `D` (past `U` its
flags do not matter; a node exactly at `U` is a sample and must be
trusted). A few lines; no other rule changes.

**PR 2's code review, round 2 (2026-10-05): NaN cells, and the refusal's
own type.** `@reviewer` found that `burn_reach` builds its data mask as
`raw != m.nodata`, all true when `nodata` is None
(`src_python/tin_engine/burn.py@96b3881:134`), while the core counts a NaN
cell as NoData whatever the sentinel (`include/terrain/raster/raster.hpp`,
`is_nodata`) and `RasterMeta` never holds a NaN sentinel
(`src_python/tin_engine/io/models.py`, `nodata` with `allow_inf_nan=False`;
`io/geotiff.py` reads a NaN tag as no sentinel). So a float DEM gapped with
NaN reaches `burn_reach` with `nodata` None, every NaN cell is "data", and
the reviewer's probe gets no refusal, `lowered_max_m` NaN, the sensitivity
well posed and NaN in the catchment file. The design now says a NaN cell is
NoData ("Following the river", step 2, "No NoData on the chain"). Separately,
`_burnt_flood` turns any `ValueError` from `burn_reach` into a
`CatchmentError`, and stage B's `except CatchmentError: pass`
(`src_python/tin_engine/catchment.py@96b3881:290`) would then swallow an
ordinary bug that happens to raise `ValueError`. **Red** (`@tester`, lean):
in `test_burn.py`, the parametrized NoData-on-the-chain test gains a third
case, float32 with `nodata` None and `z[20, 21]` NaN, refused with the gap
message; the same NaN case at the `delineate` level in `test_catchment.py`
(the gap downstream of the placed node, `CatchmentError` with the gap
message); and the refusals of `burn_reach` (the gap, the placed position
outside the window, no data node near the reach) are pinned as
`burn.BurnRefusal`, a `ValueError` subclass. The NaN cases must fail on
`96b3881` by not raising; the subclass pins fail on the missing name. Say
so in the handback. **Green** (`@developer`): in `burn.py`, the mask
excludes NaN, for example `~np.isnan(raw)`, combined with `raw != m.nodata`
when `nodata` is set; it is the one mask the cross-section, the chain
check and the extension already share, so nothing else changes there.
`class BurnRefusal(ValueError)` in `burn.py`, raised by those three
refusals; `_burnt_flood` catches `BurnRefusal` only and turns it into a
`CatchmentError` with the same words. A few lines; no other rule changes.

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
  the 2025-version of the HRD". Parsing its Table 1 (pdftotext, one regular
  expression) gives exactly 140 rows.
- **What Table 1's columns mean** (from the PDF's "Explanations to Table 1"):
  "Regine area" is the river basin number and "Main no" the main number.
  "Version" is the version of the discharge parameter ("1 = 1001.1"). Then
  come the station name, the first year of daily data in the series ("Record
  start daily data"), NVE's recommended first year for analysis ("HRD start
  daily data"), and the same two for fine-resolution data. Last come flags
  marking a series, or some of its years, as not recommended for one kind of
  analysis: spring floods, floods, summer and winter low flow, monthly flow
  or annual flow.
- **No coordinates or polygons in it.** The station number is
  `<regine>.<main>.0`: the HRD row "2 11 0 Narsjø" is NVE station `2.11.0`.
  Its third column is the discharge series version (parameter 1001, version
  0), not part of the station number. Checked on Knappom (row "2 142 1"),
  which is station `2.142.0`, and on Narsjø, whose layer 14 record (below)
  has discharge series versions 0, 1 and 2.
- **The wider set.** NVE's map service lists 862 stations with an active
  discharge series (layer 14). 364 of them have regulation degree 0 for both
  area and reservoir (`reguleringsgradareal`, `reguleringsgradmagasin`). 136
  have no regulation degree recorded at all, and none has one recorded and
  the other missing. All 140 HRD stations are in layer 14: 127 of them among
  the 364, none among the 136. The HRD adds NVE's quality review; the others
  do not have it. Ruled (Question 1 of the first version): **the 140 HRD stations**; the 364 are a
  later option.

### NVE's polygons: "Totalnedbørfelt til målestasjon"

- **The service.** NVE's ArcGIS map service `HydrologiskeData3`,
  `https://kart.nve.no/enterprise/rest/services/HydrologiskeData3/MapServer`:
  layer 0 `Malestasjoner` (points: `stasjonnr`, `stasjonnavn`,
  `totalt_feltareal_km2`, `stasjonstatus`, `vassdragsnr`,
  `elvenavnhierarki`, ...), layer 14 `Vannforing_aktiv` (active discharge
  series: `stasjonnr`, `versjon`, regulation degrees), and layer 38
  `Malest_totalnedb` (polygons: `stasjonnr`, `nedborfeltaareal_km2`,
  `oppdateringsdato`). All in EPSG:25833, the CRS of DTM10. The query
  endpoint returns GeoJSON in the asked CRS with
  `.../38/query?where=stasjonnr in ('2.11.0',...)&outFields=stasjonnr,nedborfeltaareal_km2,oppdateringsdato&outSR=25833&f=geojson`,
  `maxRecordCount` 2000. The service description (Norwegian) calls the
  polygon "the catchment upstream of the gauging station".
- **All 140 found.** Layer 0 has every one of the 140 station numbers (all
  `stasjonstatus` 1, active). Layer 38 has a polygon for every one: 130
  stations have one polygon, 10 have three, which differ by their update date
  and slightly in area (Narsjø: 119.74, 119.43, 119.74 km²). Rule: **the
  newest `oppdateringsdato` wins**, ties to the larger `objectid`; the chosen
  date is recorded.
- **Shape.** All 140 valid by shapely; one multipart (`19.79.0` Gravå, two
  parts), none with holes; median 969.5 vertices (newest polygon per
  station). Polygon area by shapely against the layer's
  `nedborfeltaareal_km2`: at most 0.005 km² apart. Against the station
  layer's `totalt_feltareal_km2`: from −2.8 % to +6.9 % (`11.4.0` the
  largest). So the two NVE areas are not the same number, and the polygon,
  not the attribute, is the reference here.
- **What the polygons are.** Not documented on the layer. NVE's catchment
  tool NEVINA (user guide, August 2025,
  `https://publikasjoner.nve.no/diverse/2025/Brukerveiledning.i.Nevina2025.pdf`)
  generates catchments "either fully automatically from the base data or a
  combination of automatically generated and catchments upstream of the point
  taken from REGINE", on a 20 m DEM since version 4. It warns that REGINE
  boundaries "can also be wrong". So the reference is NVE's best boundary,
  not ground truth, built on a coarser DEM than ours. Disagreements of a cell
  or two along the divide are expected, and the acceptance is built for that
  (below).
- **Station point against polygon.** 98 stations lie inside their polygon,
  42 outside it, every one of the 42 within 233 m of it. Inside, the distance
  to the polygon's boundary has quartiles 20 m, 75 m and 355 m (Python's
  `statistics.quantiles`, default method; max 4.4 km: lake gauges and
  stations in the middle of a wide valley floor). So the station point is
  near the outlet, but rarely on the DEM's flow line.

### NVE's river network: ELVIS ("Elvenett")

- **The service.** NVE's map service `Elvenett1`,
  `https://kart.nve.no/enterprise/rest/services/Elvenett1/MapServer`, same
  host and query interface as layers 0, 14 and 38, EPSG:25833,
  `maxRecordCount` 2000. Layer 2 `elvenett` is the complete network, as
  polylines: 1,954,539 segments (`returnCountOnly`, 2026-10-04). Its fields
  include `objekttype` (single-line river `ElvBekk`, river centreline
  `ElvBekkMidtlinje`, lake centreline `InnsjøMidtlinje`; see "The
  data are not clean" below), `strekninglnr` (the segment's national serial number),
  `elvid` (one id per branch of the network), `vassdragsnr` (the watercourse
  number of the REGINE unit), `elvenavn`, `elvenavnhierarki`,
  `elveordenstrahler`, `vatnlnr` (the lake number, for lake centrelines) and
  `til_utlop`. Layer 1 `hovedelv` holds only the main rivers.
- **What it is.** NVE's product sheet (`Produktark: Elvenett - ELVIS`, NVE,
  12.05.2017, from Geonorge's register, sha256 `58c3b956...1aeb`) says that
  ELVIS is "derived from the water theme of N50 map data", at scale
  1:50,000. Lakes and two-line rivers get "a mathematical centreline". The
  network "contains information about the direction of flow", and it is
  meant for upstream and downstream analyses and "as an aid when generating
  hydrologically correct elevation models" (translated). Version 2 is being
  rebuilt from newer N50 data (about 65 % done in May 2017).
- **Licence.** The same sheet: distributed under NLOD, with the source text
  "Kilde: NVE" required on publication. The same terms as layers 0 and 38.
- **How it relates to the stations.** The station layer carries no segment
  or branch id, only `vassdragsnr` and `elvenavnhierarki`. Measured over the
  140 stations, querying layer 2 within a 1 km square round each station
  point (28 s for all 140): 139 have a line within 500 m. The distance from
  the station point to the nearest line has median 20 m, 75th percentile
  53 m, 90th 158 m, 95th 252 m and max 461 m; 31 are within 10 m, 76 within
  25 m, 101 within 50 m, 116 within 100 m and 132 within 250 m. The one
  without, `311.4.0` Femundsenden, is a lake gauge 612 m from the nearest
  line (a 6 km square). For 113 of the 139 the nearest line carries the
  station's own `vassdragsnr`. The other 26 pass a looser test (a shared
  watercourse-number prefix, or the river name), which was not checked one
  by one. For 46 of the 139 the nearest line is a lake centreline (lake
  gauges). For 16 stations, a second branch (another `elvid`) lies within
  100 m of the station point, and for 5 within 50 m: these are the
  confluence cases. Measured from the mapped position `P` instead (the foot
  of the station on the nearest line, which is what the design's flag
  uses), 19 stations have a second branch within 100 m.
- **The data are not clean** (measured 2026-10-04; `outStatistics` grouped
  by `objekttype` over the whole layer, and the 1 km samples of the 140
  stations). (a) `objekttype` has **25 distinct values**, including null
  (94 features), a blank (2) and strays (`SK`, `20.08.2014`): eight lake
  spellings (`InnsjøMidtlinje`, `InnsjoMidtlinje`, `InnsjøMidtlinjeReg`,
  `InnsjøMitlinje`, `InnsjøMidtloinje`, `Innsjømidtlinje`, `InnsjoMidtlin`,
  `InnsjøRegulert`), several river and fictive-link spellings
  (`ElvBekkRegulert`, `FiktivElv`, `ElvelinjeFiktiv`, `BreMidtlinje` (a
  glacier), ...). A nearest line of the 1 km samples has a null type
  (`152.4.0`: five features, all with `vatnlnr`, the lake number, set). The
  service sends `vatnlnr` = 0 for "no lake" (956,447 features of the layer
  have 0, 444,793 a positive number, 551,648 none), so **set means not null
  and not 0**. In the samples a positive lake number is on 419 of the 442
  lake lines (9 have 0, 14 none) and on 1 of the 949 river lines (489 have
  0, 458 none); over the whole layer, on 58 of the 94 null-type features
  (none has 0) and on 1,779 river-typed ones. Both blank-type features
  (`objectid` 11166506 "Cap'pirjåkka" and 11513676) have 0 and are rivers.
  So the lake number decides only where the type is null or blank.
  (Re-measured 2026-10-04 by `@architect`: `returnCountOnly` queries on
  `Elvenett1/MapServer/2` with `where` `vatnlnr = 0`, `vatnlnr > 0`,
  `vatnlnr IS NULL`, and each crossed with `objekttype IS NULL` and
  `objekttype = ' '`; the samples re-fetched as squares of half-side 500 m
  round the 140 layer-0 points, 1,396 distinct features. The design rounds
  counted 0 as set, which gave "428 of the 442" and "496 river features".) (b) **Exact copies**: of the `strekninglnr` values
  shared by more than one `objectid` in the samples (14), 13 are exact
  copies of one geometry (9 pairs at `82.4.0`, 4 groups of five at
  `139.35.0`; 25 extra features in all), the copies sharing one `elvid`; the
  14th, at `79.3.0`, is two different geometries (a lake centreline, 334 m
  apart by Hausdorff distance), the same `strekninglnr`, one `elvid`. (c)
  **Forks**: in the sampled lines, where the chain of the nearest line's
  `elvid` is followed end to start within 1 m, one station (`19.80.0`) meets
  a point with two different continuations (a sample, not the reach
  envelopes). The design takes each case as a rule below.
- **Direction of digitising.** For 138 of the nearest lines (100 m or
  longer, both ends on one tile), the DEM is lower at the line's last vertex
  than at its first for 97, higher for 2, and within 0.5 m for 39 (mostly
  flat lake surfaces). So the lines are digitised downstream, as the sheet
  says. The design still checks each reach against the DEM (below).
- **The DEM does not follow the mapped line.** Along the 90 of those lines
  that are rivers rather than lakes and drop from start to end, the DEM was
  sampled every 5 m. The largest climb along the direction of flow has
  median 1.6 m, 75th percentile 3.7 m and 90th 5.0 m, with a maximum of
  51.6 m; 42 of the 90 climb more than 2 m somewhere. A line at 1:50,000 sits
  off the 10 m DEM's valley floor by tens of metres in places, and a road
  embankment in the DEM can dam a mapped river. Both decide the design of
  the burn below: the line is first moved to the valley floor, and the DEM
  is lowered only where it still climbs.
- **Fetched like the polygons.** Per station, one envelope query on layer 2
  with `outSR=25833&f=geojson` and a short `outFields` list. The whole
  network is not fetched (1.95 M segments).

### Coverage by DTM10_UTM33_20260925

- **The DEM.** 254 tiles at 10 m, EPSG:25833, spaced 50 km apart; 243 of
  them are 5051 × 5051 nodes (columns × rows), so neighbours overlap by 51
  nodes. Eleven differ: seven of the eight shifted tiles (below) are
  5052 × 5053, the eighth, `7507_4`, is 5052 × 4103, and three normal tiles
  are smaller: `7305_3` 5051 × 2881, `7405_1` 3521 × 5051 and `7405_2`
  3511 × 5051 (image sizes read with `tifffile`, 2026-10-04). Together they
  span x −100 km to 1150 km, y 6400 km to 7950 km. NoData is −32767.
- **Which stations it covers.** Sampling each tile every 10th node (100 m)
  and each polygon at those samples: **138 of 140 polygons meet no NoData
  sample**. Two do: `244.2.0` Neiden (2,945 km², 2.6 % of its samples NoData)
  and `234.18.0` Polmak nye (14,171 km², 0.02 %). Both cross into Finland.
  The flood refuses a catchment that reaches NoData (increment 22, "The
  window", step a), so these two are expected to be refused, and are
  reported as such, not as failures. A 100 m sample can miss a NoData strip
  narrower than 100 m; the run itself is the check.
- **The half-cell-shifted tiles.** Eight of the 254 tiles (`7304_1`,
  `7507_4`, `7606_2`, `7707_1`, `7707_3`, `7807_2`, `7807_3`, `7808_3`) have
  their nodes 5 m off the others' in x, and in x only. In the world files
  (`.tfw`, the centre of the first cell) the normal tiles' x is a multiple of
  10 m (`6400_1`: 49750) and the shifted tiles' is 5 m off one (`7304_1`:
  449745); the y values of both are multiples of 10 m. **The cause is a
  resample, not a label**: all 254 tiles declare PixelIsArea and their
  GeoTIFF tags agree with their `.tfw` files, and where a shifted tile
  overlaps a normal neighbour, the median elevation difference is smallest
  at the declared position. **The method**: the shifted tile's node centres
  are moved by −5, 0 and +5 m in x; the neighbour is interpolated linearly in
  x at those positions (the rows share the y lattice); the median absolute
  difference is taken over up to 400 evenly spaced rows of the overlap,
  NoData excluded. **`7707_1`**, against all four of its normal neighbours
  (those whose bounding box overlaps it by at least 200 m both ways:
  `7707_2`, `7707_4`, `7708_3`, `7708_4`): 0.02 to 0.11 m at the declared
  position, against 0.19 to 0.45 m moved either way. **`7707_3`**, against its seven
  normal neighbours by the same rule (re-run by @architect): three of its four
  long overlaps favour the declared position clearly (`7607_4` 0.02 m,
  `7706_2` 0.07 m, `7707_4` 0.05 m, against 0.50 to 1.02 m moved), but the
  fourth, `7707_2`, is a **near tie**: 0.25 m declared, 0.24 m moved −5 m,
  0.73 m moved +5 m; of the three 500 m corner overlaps, `7607_1` favours the
  declared position (0.07 against 0.13 and 0.27) and two are near ties
  (`7706_1` 0.23 against 0.70 and 0.24; `7606_1` 1.27 against 3.91 and
  1.32). The other
  tiles were compared against one neighbour each, by the main session
  (2026-10-04, not re-run here except `7507_4`, which @reviewer re-measured
  in round 4): `7807_2` 0.11 m against 0.38 and 0.38; `7507_4` 0.00 against
  0.06 and 0.06. Three tiles touch neighbours only on flat ground, where the
  test cannot tell, and `7304_1` differs by about 2 m at any shift. The
  near ties do not change the conclusion: all four of `7707_1`'s overlaps
  and three of `7707_3`'s four long ones favour the declared position
  clearly, and none favours a moved one clearly. So the data were resampled onto
  a shifted grid, and the later fix resamples them onto the common one;
  moving the label would misplace them by 5 m. Increment 15a's mosaic refuses a
  selection that mixes the two lattices. Nine HRD polygons meet both kinds of tile:
  `156.15.0`, `196.11.0`, `206.3.0`, `208.2.0`, `208.3.0`, `209.4.0`,
  `212.49.0`, `213.2.0`, `223.2.0` (tile footprints intersected with NVE's
  polygons, newest version each). They are **expected refusals** in this
  increment (ruled by Ola, 2026-10-04; see "Ola's rulings"), and a later
  increment resamples the eight tiles onto the common lattice. Two more
  (`156.24.0`, `213.4.0`) lie on shifted tiles only, one lattice, and run.
  A window can reach a shifted tile where the catchment does not; such a
  refusal is a finding to explain, not an expected one.
- **Tiles per catchment**, for the 138: 70 within one tile's footprint, 39
  in two, 8 in three, 18 in four (across a tile corner), 1 in five
  (`139.35.0` Trangen), 2 in nine (`212.10.0` Masi, `311.6.0` Nybergsund).
  The overlap strips count as both tiles, so "two" includes catchments that
  only reach into a 510 m overlap. 15a's mosaic handles all of these except
  the nine above. The acceptance reports the result by tile count, so a seam
  artefact would show as misses gathered at multi-tile stations.
- **Sizes.** From 0.44 km² (`20.11.0` Tveitdalen, about 4,400 nodes) to
  14,171 km²; median 131 km² (130.98, newest polygons); 12 under 10 km², 44
  from 10 to 100, 75 from 100 to 1000, 9 over 1000. Sum 61,016 km², about
  610 M catchment nodes at 10 m.

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
- **What is committed.** Ruled (Question 2 of the first version):
  - the list of 140 HRD rows, four columns (station number, discharge series
    version, name, HRD start year of daily data), as a data file in the
    package, credited to NVE
    in `NOTICE.md` and in the file's header: small, factual, and the one
    piece that cannot be fetched from a service (it lives in a PDF);
  - the per-station result table of the acceptance, which holds NVE's areas
    as numbers;
  - **not** the polygons (8.2 MB for all 160 versions), the station points
    or the river lines: `rasputin fetch-stations` fetches them, and records
    the fetch date, the URLs and each file's sha256 in a manifest, as 22's
    acceptance recorded its one polygon. The river lines fall under the same
    ruling as the polygons (fetched, not committed); their credit is "Kilde:
    NVE", as the product sheet asks.

### Discharge

NVE's HydAPI (`https://hydapi.nve.no/api/v1/`) serves the discharge series;
it needs a free API key (`/Stations` answers 401 without one) and is under
NLOD. **Out of scope now** (ruled, Question 3 of the first version): no discharge
is fetched or stored. Each station keeps the two keys that join it later: the
station number (`2.11.0`, HydAPI's `StationId`) and the HRD's discharge series
(parameter 1001 and its version, `1001.0`).

## Data use

Ola asked that nothing dodgy is done with NVE's data. The rules, each held by
a test of `fetch/nve.py` (PR 3), except the `NOTICE.md` credit below:

- **Sources.** The HRD report (the PDF named above); `HydrologiskeData3`
  layers 0 (`Malestasjoner`) and 38 (`Malest_totalnedb`); `Elvenett1` layer 2
  (`elvenett`). HydAPI is not used and no API key is stored anywhere.
- **Licence.** NLOD; "Kilde: NVE" is in `NOTICE.txt` of every fetch, in the
  packaged CSV's header (each tested) and in `NOTICE.md` (checked by hand by
  `@reviewer`, block "PR 3, after PR #178's CI" above). NVE disclaims liability for errors
  in the data and their use; the README of the acceptance says so too.
- **Committed**: the 140-row HRD list (four columns, from a PDF) and the
  acceptance's `results.csv` (NVE's areas as numbers). **Fetched** into
  `../rasputin_data/nve_hrd/` and never committed: polygons, station points,
  river lines, manifest.
- **Only the fields the design needs, by explicit allow-list per layer, in
  every request (never `outFields=*`).** Layer 0: `stasjonnr`, `stasjonnavn`,
  `totalt_feltareal_km2`, `stasjonstatus`, `vassdragsnr`, `elvenavnhierarki`.
  Layer 38: `stasjonnr`, `nedborfeltaareal_km2`, `oppdateringsdato`,
  `objectid`. Layer 2: `objectid`, `objekttype`, `strekninglnr`, `elvid`,
  `vassdragsnr`, `elvenavn`, `vatnlnr`. **Never collected**: `stasjoneier` (the
  owner), ELVIS's `oppdatertav` (editor ids, some look like personal
  initials), `globalid`, and layer 38's discharge normals.
- **Query volume, kept small.** Only `fetch-stations` uses the network: about
  8 batched queries for layers 0 and 38 (40 stations each) and 140 ELVIS
  envelope queries, sent one at a time, with 23a-2's retries and back-off. A
  fetch is reused, not repeated: if the files exist, nothing is requested
  unless `--refresh`. `station-catchments` and `catchment` never touch the
  network. The client identifies itself (`User-Agent: rasputin/<version>`)
  instead of the library's default.

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
Increment 22's grep for catchment terms still returns no files. For the
river network and the burn:

```
$ git grep -liE "burn|elvenett|elvis|hydrography|river" legacy-archive -- legacy
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/gml_repository.py
```

Both hits are land-cover class names for burnt ground (`"No data (burnt
areas, clouds,…)"`, `burnt = 334`), not stream burning or rivers.

**Literature.** Each reference was checked against Crossref on 2026-10-04
(title, authors, venue, DOI), except where marked; abstracts or the pages
named were read, no paper in full.

- **Jenson 1991**, "Applications of hydrologic information automatically
  extracted from digital elevation models", *Hydrological Processes*
  5(1):31-44, doi:10.1002/hyp.3360050104. The source of the rule that moves
  an outlet to the **nearest** stream cell within a distance. Round 1 of
  this design misattributed to it the rule of the **largest** accumulation
  within a radius (the "snap pour point" of common GIS tools), and chose
  that rule. Both are point-to-raster rules: neither looks at where the
  river is mapped. **Used here only as the fallback** for a station with no
  mapped river line within the map radius (one of the 140, Femundsenden, at the
  default 500 m map radius; any station list a user brings without a river file). The
  acceptance also runs it as a **comparison**, on the stations that end up
  `miss` or `uncertain`, so that its effect is measured.
- **Lindsay, Rothwell and Davies 2008**, "Mapping outlet points used for
  watershed delineation onto DEM-derived stream networks", *Water Resources
  Research* 44(8):W08442, doi:10.1029/2007WR006507 (abstract read, from
  OpenAlex). It states this increment's problem: "Outlet point positions
  taken from hydrometric stations commonly do not coincide with stream
  locations extracted from digital elevation models". It proposes AORA,
  which "uses water body names to identify locations for outlet
  repositioning" and had "the fewest repositioning errors" against "two
  existing automated techniques" over 993 stations. The abstract names
  neither of the two. The first author's own documentation of the tools
  (WhiteboxTools, `JensonSnapPourPoints`, doc comment in its source) says
  the nearest-stream rule "should be preferred" over the largest-stream
  rule, which near a confluence "may re-position outlets on the main-trunk
  stream". **What is taken:** AORA's principle, that a gauge goes where the
  mapped water body says it is, not where a raster quantity peaks. **What
  differs:** the mapped river's geometry places the gauge, not only its
  name; the station's watercourse number filters the candidate lines (the
  name is kept as a second filter); and the DEM is conditioned along the
  mapped reach (below). AORA also checks the reported area; that is
  rejected here (next item).
- **Lehner 2012**, GRDC Report 41, "Derivation of watershed boundaries for
  GRDC gauging stations based on the HydroSHEDS drainage network" (BfG, not
  on Crossref; the GRDC page
  `https://grdc.bafg.de/products/basin_layers/watershed_boundaries/` read):
  candidate outlets within 5 km, chosen by "the reported catchment area of
  the GRDC station and the distance between the station coordinates and the
  derived outlet", failures inspected by hand, at 15 arc-seconds (about
  500 m). **Färber et al. 2025**, "GRDC-Caravan: extending Caravan with data
  from the Global Runoff Data Centre", *ESSD* 17(9):4613-4625,
  doi:10.5194/essd-17-4613-2025, is where the GRDC page points for the
  current method (its section 2.3; not read here). **Rejected, and why**
  (Ola, 2026-10-04): the flow accumulation grows monotonically downstream,
  so any rule that maximises area, or that matches a reported area among
  candidates near a confluence, biases the gauge downstream. That bias is
  fatal for residual inflow (the inflow between two points on one river).
  Matching a reported area also borrows the reference answer that the
  acceptance compares against.
- **Hellweger 1997**, "AGREE — DEM surface reconditioning system", Center
  for Research in Water Resources, University of Texas at Austin (a web
  report; not on Crossref, and not found online in this round, so
  **unverified**: cited as the usual name for the method). Stream burning:
  lower the DEM along the vector stream lines, with a buffer that slopes
  towards them, so that derived flow follows the mapped network.
- **Lindsay 2016**, "The practice of DEM stream burning revisited", *Earth
  Surface Processes and Landforms* 41(5):658-668, doi:10.1002/esp.3888
  (abstract read, from OpenAlex). Names the artefact of common burning:
  "topological errors resulting from the mismatched scales of the
  hydrography and DEM data sets", in particular "erroneous stream piracy
  caused by the rasterization of multiple stream links to the same DEM grid
  cell". Its TopologicalBreachBurn prunes the network to the DEM's
  resolution, "restricts flow within individual stream reaches", and gives
  the larger stream priority where two share a cell. A plain burn
  (FillBurn) lost accuracy at coarse resolution (kappa 0.953 down to 0.490,
  against 0.952 to 0.921). **What is taken, and what differs:** only one
  chain of the network is burnt, the gauge's own reach, so no two links are
  rasterised together and no piracy between links can arise. The reach is
  first moved onto the DEM's valley floor, inside a narrow corridor round
  the mapped line, because ELVIS is at 1:50,000 against a 10 m DEM (the
  scale mismatch Lindsay names). The DEM is then lowered by breaching only
  where it still climbs along the reach. This is a local, minimal form of
  burning, not AGREE over the network. **The departure drops** the effect a
  network-wide burn has on divides (where the DEM puts a divide the mapped
  network crosses). Burning the whole network is a later increment (Question
  3, ruled by Ola, 2026-10-04).
- **Soille, Vogt and Colombo 2003**, "Carving and adaptive drainage
  enforcement of grid digital elevation models", *Water Resources Research*
  39(12):1366, doi:10.1029/2002WR001879 (Crossref and abstract checked
  2026-10-04). Carving: instead of filling a pit, lower the terrain along a
  descending path from it; and "adaptive drainage enforcement", where known
  river networks are imposed on the DEM "only in places where the automatic
  river network extraction deviates substantially from the known networks".
  **What is taken:** the principle that the burn lowers the DEM only where
  it disagrees with the map (the descent lowers a node only where the DEM does
  not already fall). **What differs:** one chain, not the network; the chain
  is first moved onto the valley floor; and, unlike carving's paths found by
  a flood from the outlets, the descent here is along a given chain, so
  whether the flood then drains along it is checked (`flow_to`), not assumed.
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

**Novelty.** None is claimed. Searched: the works above, Crossref,
OpenAlex and web searches for automatic delineation of gauge catchments
compared with official polygons (overlap, Jaccard, area ratio), outlet
placement on mapped hydrography, and stream burning; NVE's NEVINA guide, the
ELVIS product sheet, and the HRD report. The per-station sensitivity of the
catchment area along the river (below) was not searched as a method of its
own. It is presented as a diagnostic, not a contribution, and a write-up
would search for it first. Not found in that search: a published, per-station,
two-way overlap comparison of a fully automatic 10 m delineation against
NVE's polygons for the HRD stations. That is a gap in what was searched, not
a claim; a write-up would search again (Scandinavian journals, NVE reports,
theses) before saying more.

## The design

### Data flow

```
rasputin fetch-stations nve-hrd --out-dir D                       [network]
  fetch/nve.py   reads the packaged HRD list (140 rows), queries NVE's
                 layers 0 and 38 and ELVIS layer 2 through fetch/http.py,
                 newest polygon per station, writes D/stations.geojson,
                 D/reference.geojson, D/rivers.geojson, D/NOTICE.txt,
                 D/manifest.json (URLs, date, sha256)

rasputin station-catchments --dem ... --stations D/stations.geojson
        --rivers D/rivers.geojson [--reference D/reference.geojson] --out-dir O
                                                                    [offline]
  cli.py         paths stop here: io/station_set.py and io/rivers.py read the
                 files into Station and RiverSegment models and reference
                 polygons; repository_for(dem)
  catchment_batch.run_batch(request, repository, stations, segments, references, sink)
     for each station, in file order, one at a time:
       gauge.place(Gauge(station), segments) -> Placement | None   [pure, shapely, no DEM]
       catchment.delineate(CatchmentRequest(seed=station, seed_crs=crs,
                                            reach=placement.reach), repository)
         window loop (22's), seeded by the placed node:
           burn.burn_reach(window, reach) -> burnt window, GaugePath  [numpy, no shapely]
           _core.upstream(burnt window, placed node) -> outline -> reduce   (22)
         _core.accumulate(burnt window)                                 [C++, new]
         sensitivity.assess(counts, flags, path, uncertainty)         [pure, numpy]
       reference.agreement(result, reference polygon)   [pure, shapely]
       reference.classify(agreement, sensitivity)       [pure]
       sink.catchment(station, result)  -> O/<station>.geojson  (io/geojson.py)
       sink.row(StationResult)          -> O/results.csv
     reference.summarise(rows) -> O/summary.json, and stderr
```

No path, file, URL or CRS crosses into C++: `accumulate` sees the same
`RasterView` `upstream` does and returns arrays (CLAUDE.md §2, I/O boundary).
The burn edits an elevation array in Python and hands the result to the same
binding; it takes plain arrays, so it imports no shapely, and `gauge.py`
imports shapely but sees no DEM. No new dependency: the network is `urllib`
inside `fetch/http.py` (its existing `get_text`), lines and polygons are
shapely, the HRD list is read with `csv` and `importlib.resources`.
`tools/check_prohibited_deps.py` covers the rest.

### Accumulation, from the same flood (C++)

`include/terrain/hydrology/accumulate.hpp`, header-only, next to
`upstream.hpp`, a template on `RasterSource`:

```cpp
namespace terrain::hydrology {
struct AccumulateOutcome {
    std::vector<std::uint32_t> count; // row-major; 0 on NoData, else the number of
                                      // nodes that drain through the node, itself included
    std::vector<std::uint8_t> reach;  // bit 0: that catchment touches the window's edge,
                                      // bit 1: it touches NoData (upstream's two flags)
    std::vector<std::uint8_t> flow_to; // the neighbour the node drains to (its flooder), as
                                       // 3*(dr+1)+(dc+1) with dr, dc in -1..1 (4 is never
                                       // used); 255 for an outlet and on NoData
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
  pass over the reverse order adds each node's count into its flooder's. The
  flooder byte is returned as `flow_to`: it is the drainage tree itself, and
  it lets a caller check, node by node, that a chain of nodes drains along the
  chain (the burn, below) and, later, that one gauge's chain of flooders
  passes through another's (residual inflow). It is the byte the pass already
  holds, so returning it costs no memory beyond handing the array to the caller.
- **The flags ride the same pass.** A node starts with bit 0 when it lies
  within one node of the window's edge, and with bit 1 when it is an outlet
  beside NoData or is an 8-neighbour of such an outlet; the pass ORs each node's bits
  into its flooder's. These are exactly the conditions under which
  `upstream`'s `touches_edge` and `touches_nodata` fire for a seed at that
  node, so a count whose bits are clear is the node's whole catchment, and a
  count with a bit set is a lower bound. Sensitivity (below) needs this to
  know which counts it can trust.
- **To keep the two floods from drifting**, the outlet set-up and the
  queue loop move out of `upstream.hpp` into `hydrology/flood.hpp`,
  `detail::flood(z, state, on_reach)`. The outlets are every node with data
  on the window's edge or beside NoData (an 8-neighbour of a NoData node);
  `on_reach(j, j)` announces each outlet `j` as it is queued, and
  `on_reach(i, j)` is called when popped node `i` reaches node `j` first. A
  reached node is set to *out* before the call, which may relabel it. The
  flood returns the outlets beside NoData, ascending, since both callers
  need them for `touches_nodata`. `upstream` labels in it, `accumulate`
  records in it. `upstream`'s behaviour and suite are unchanged.
- **The oracle is exact**: for every node `c` with data,
  `accumulate(z).count[c] == upstream(z, {c}).nodes_in`, and the two bits
  of `reach[c]` equal `touches_edge` and `touches_nodata` of the same call.
  Both are "the nodes whose chain of flooders passes through `c`". `flow_to`
  is held to it too: `count[c] == 1 + sum(count[d])` over the nodes `d` whose
  `flow_to` points at `c`, and `upstream(z, {n}).mask[d] == 1` for the
  neighbour `n = flow_to[d]` of every node `d` that has one. A test can
  check it node for node on small random DEMs with pits, flats and NoData,
  which is the invariant-critical suite of this increment (below). Also: the
  outlets' counts sum to the number of nodes with data.
- **Limits.** Refused with `std::length_error` when the raster has 2^32 nodes
  or more (the count and the order are 32-bit); the binding maps it to
  `ValueError`. Memory per node: 1 byte for the flooder (`flow_to`), 4 for the order, 4
  for the count, 1 for the bits, beside the elevation array and the queue.
  Serial, deterministic; the binding releases the GIL and returns numpy
  arrays of the raster's shape that the result owns.

### Placing the gauge (`gauge.py`, pure)

**The rule (Ola, 2026-10-04).** A gauge is placed by where it physically is:
on NVE's mapped river line, then on the DEM's flow path along that line,
which is burnt where the DEM and the line disagree. **No rule here looks at
an area, a flow count or a reported catchment size.** Flow accumulation grows
downstream, so any rule that prefers more area (the largest count within a
radius, the area closest to a reported one) is biased downstream, and a
catchment pushed downstream carries inflow the gauge never saw: fatal for
residual inflow (below). The accumulation count is used only to *read* the
area at a place already chosen (the sensitivity), never to choose the place.

`gauge.place(gauge, segments, *, map_radius=500.0, reach_up=1000.0)`
returns a `Placement` or `None`; `Gauge` is a frozen (x, y, watercourse, river)
so that `gauge.py` does not import the station reader. `RiverSegment` is the frozen model of one
ELVIS line: `objectid`, `elvid`, `vassdragsnr`, `name`, `objekttype` (the raw
value, kept as served, `None` allowed), `kind`, and the line's vertices in the
file's CRS. **`kind` is a total mapping** (in `io/rivers.py`, where the model
is built): `lake` when the casefolded `objekttype` starts with `innsj` (all
eight lake spellings of "The data are not clean"), or when `objekttype` is
null or blank and `vatnlnr` is set, which means not null and not 0, since
NVE sends 0 for "no lake" ("The data are not clean"; `152.4.0`'s five lines
have 495; both blank-type features have 0 and are rivers); `river` for
every other value, the fictive links, the glacier lines and the strays
included. A stray never raises; the raw value stays in the model and in the
station's row (`objekttype` of the chosen line), so a surprising class is
visible, and `@tester` pins the mapping on all 25 values.

1. **Candidates**: segments within `map_radius` of the station point.
   Default 500 m: 139 of the 140 stations have a line within 461 m (the one
   without, `311.4.0` Femundsenden, is a lake gauge 612 m from the nearest).
2. **Which river**: the station's own watercourse number says which river it
   is on. Tiers: the segment's `vassdragsnr` equals the station's (113 of the
   139 nearest lines), then shares its prefix or its river name (26 more,
   not checked one by one), then any. **Prefix** (pinned 2026-10-05 for
   PR 2's green step): the segment's number, not equal to the station's, is
   a leading part of it as a string (`station.startswith(segment)`), and the
   two share the main number before the first dot (`"002.A"` is a prefix of
   `"002.AB"`; `"002"` is not, and `"002.AB"` is not one of `"002.A"`: the
   station's area drains into the segment's river, not the other way). The
   prefix tier comes before the name tier. **The nearest segment of the best
   non-empty tier wins**, and the tier is reported (`placed_on`: `number`,
   `prefix`, `name`, `any`), so a station placed on a line that is not its
   own is visible in the table, not hidden in a choice.
3. **The mapped position `P`** is the nearest point of the chosen line to the
   station point (the foot of the perpendicular, or an end vertex). The
   station-to-line distance `d` is reported. **Position uncertainty
   `U = min(max(d, 30 m), map_radius)`**: the gauge's coordinates cannot
   place it along the river better than their own distance from the river,
   and 30 m (three cells) is the floor the DEM's resolution gives. The
   measured `d` has median 20 m, 75th percentile 53 m and 90th 158 m, so
   half the stations get the 30 m floor and a tenth get 160 m or more.
4. **The reach** is the chain of segments with `P`'s `elvid`, joined where an
   end of one lies within 1 m of the start of the next (a segment is a link
   of one river, digitised downstream). **Exact copies are dropped first**:
   within one `elvid`, segments whose vertex lists are equal to 1 cm are one
   segment, the smallest `objectid` kept (`io/rivers.py`'s `read_segments`
   does it, so a user's file is cleaned as a fetched one is, and returns the
   count dropped, which the commands print; see "Ola's rulings", last block). **A fork stops the chain**: where, after that, two
   different segments of the `elvid` continue the chain (or two lead into its
   first one), the chain stops there, flags `reach_fork` and reports the
   metres it has; it does not pick a branch. Downstream, a fork before
   `U` is a shortfall (`reach_down_m < U`), so the station is `uncertain`
   by "Sensitivity". The chain runs from `reach_up` metres upstream of
   `P` (default 1000 m: the nearest line's median length is 751 m, so one
   segment is often not enough, and a bridge embankment a few hundred metres upstream
   dams the DEM river as much as one at the gauge) to `U + 100 m` downstream
   of it. The chain stops where no segment of that `elvid` continues it; the
   metres actually available are reported (`reach_up_m`, `reach_down_m`). The
   result is `Reach(line, at, uncertainty)` (below), in the file's CRS.
5. **Flags carried to the result, not decided here**: `lake` (the chosen
   segment is a lake centreline: 46 of 139 nearest lines), `confluence_near`
   (a segment of another `elvid` within 100 m of `P`, measured from `P`, the
   mapped position, not from the station point: **19 of 139**; from the
   station point it is 16, and 5 within 50 m) and the two distances. `None` (no segment within `map_radius`) is the case for the
   fallback.

Measured on 139 stations, 2026-10-04 (1 km square envelopes; `strekninglnr`
is shared by 14 groups of features in them, 13 of which are exact copies and
one (at `79.3.0`) two different geometries, so `objectid` is the key and exact
copies are dropped by geometry, above): the nearest line's length has quartiles 392 m, 751 m and 1196 m, maximum
4673 m. The envelopes fetched for the reach are larger (below).

```python
class Reach(BaseModel):  # frozen; in the request's seed_crs
    line: tuple[tuple[float, float], ...]  # downstream order, at least two vertices
    at: float           # metres along `line` from its first vertex: the mapped position
    uncertainty: float  # metres: the half-width of the sensitivity window
    corridor: float = 30.0  # metres: the burn's half-width round the line
```

### Following the river in the DEM (`burn.py`, numpy)

The mapped line is 1:50,000 against a 10 m DEM (the scale mismatch Lindsay
2016 names), and the DEM does not follow it: along the 90 measured river
lines the DEM climbs along the flow by 1.6 m at the median and 5.0 m at the
90th percentile, up to 51.6 m. A gauge placed on the line's nearest node
would often sit on a hillside or behind an embankment. So the reach is first
moved to the valley floor, then lowered where it still climbs:

1. **Resample** the reach at the DEM's node spacing, by arc length, from its
   first vertex: points at 0, `step`, 2 `step`, ... up to the line's length,
   with no extra point at the last vertex when the length is not a multiple
   of `step` (the chain can then end less than one step short of the line's
   end). Points outside the window are dropped, all of them, and
   `GaugePath.dropped_m` is their count times `step` (pinned at PR 2's green
   step; the window holds the whole reach plus `WINDOW_MARGIN_M` unless the
   data ends, so this is 0 in the usual case, and it is not carried into
   `GaugeResult`).
2. **Valley floor, across the line.** For each resampled point `p`, the node
   of least elevation in its **cross-section**: the nodes within `corridor`
   of `p` (distance `<= corridor + 1e-6 m`: a node exactly at the corridor's
   distance is inside) whose offset from `p` along the line is at most half
   the resample step (`|(n - p) · t| <= step/2 + 1e-6 m`, `t` the unit
   direction of the line's segment that holds `p`; a point on a vertex takes
   the segment that starts there, the last point the last segment). Ties: the
   nearer to `p`, then the smaller (row, column). A cross-section with no node
   with data (when `corridor` is below half a cell's diagonal, 7.07 m on
   DTM10, or when every node in it is NoData) takes the data node nearest
   `p`, same ties. NoData nodes are never chosen: they
   are left out of the cross-section and of the nearest-node fallback, and a
   point with no data node within its search box is a `ValueError` (pinned at
   PR 2's green step). **Why across and not the whole
   disc** (PR 2's red step, `@tester`): the least elevation anywhere within
   `corridor` lies, on a falling floor, at the downstream edge of the disc, so
   every point, and the placed node with them, moves downstream by about the
   corridor. Measured 2026-10-05 with a throwaway script (not committed) on a
   straight floor falling 0.5 m per 10 m, at 0°, 30° and 45° to the lattice:
   the disc rule moved the chosen node 6 to 8 m downstream on average at a
   15 m corridor, 22 to 26 m at 30 m and 52 to 56 m at 60 m. The
   cross-section bounds the move along the line to half a step (5 m on
   DTM10) by construction, at every corridor; on a floor that is flat
   across the line and falls along it, the least node of each cross-section
   sits towards its downstream side, a downstream lean of 3 to 4.5 m on
   average at 30° to the lattice (`@reviewer`'s probe, PR 2 code review,
   round 1). A
   downstream bias is what Ola's ruling forbids ("moving a point 100m
   downstream will always give you more inflow"), and it would grow with the
   15 m / 60 m corridor re-runs of the acceptance; across the line, the
   corridor sets how far off the line the floor may be, not how far along.
   This is the usual definition of a thalweg, the line joining the lowest
   points of successive cross-sections. The 1e-6 m is about 500 times
   float64's spacing at 1e7 m, the largest coordinate a UTM northing reaches;
   checked at Norway's northings (up to 7.9e6 m) only by that arithmetic.
   The chosen nodes are joined into one 8-connected chain
   (a straight node-to-node line between consecutive chosen nodes).
   **The chain is then made taut**, before the descent: no two of its nodes
   that are not consecutive may be equal or 8-neighbours. From the first node
   on, for each kept node `chain[k]`, let `j` be the largest index after `k`
   whose node is `chain[k]` itself or one of its 8 neighbours. If `j > k+1`,
   the nodes `chain[k+1..j-1]` are dropped, and `chain[j]` too when it is
   `chain[k]` itself (a loop: the chain went away and came back); the next
   kept node is then the one after the cut. Dropping nodes creates no new
   neighbours, so one pass leaves the chain taut, and it stays 8-connected:
   `chain[k]` and the next kept node are neighbours. **Why taut**: the flood
   names a node's flooder when it first pushes it, and a lower chain node
   pops first. On a chain with a square corner, `(2,4), (3,4), (3,5)`, the
   lower `(3,5)` pushes `(2,4)` before `(3,4)` can, so `(2,4)` drains past
   its successor and `drains` fails though the burn is correct (@reviewer,
   round 3, with `upstream.hpp`'s flood). The cut (here, `(3,4)`) removes
   every such shortcut: a chain node's only chain neighbours are then its
   predecessor, higher, and its successor. Checked 2026-10-04 with a Python
   copy of that flood (a throwaway script, not committed) on 500 random
   8-connected walks burnt into high ground: 190 of the raw walks failed
   `drains`, none of the taut ones. A straight join between two chosen nodes
   has no such corner; corners come from the joins between joins, from loops
   and from the extension (step 5), so the rule covers all three. Arc lengths
   are computed after the cut, along the chain as it stands (10 m for a
   straight step, 14.1 m for a diagonal one on DTM10), and a resample point
   whose node was dropped maps to `chain[k]`, the kept node the cut starts
   from (for a loop, the first visit).
   **No NoData on the chain** (PR 2's code review, round 1; provisional,
   question 5 for Ola, default refuse). The chosen nodes have data, but the
   straight join between two of them does not look at the nodes it passes,
   and the taut cut can map the placed node onto one of those. A NoData node
   on the chain is burnt with its sentinel as its elevation: with −32767,
   every node below it is lowered to about −32767 m; with a positive
   sentinel (3.4e38) the node itself is lowered to the floor's height and
   `lowered_max_m` reports about 3.4e38 m, which reaches the catchment file;
   the direction check's means read the sentinel too. So **after the taut
   cut, before the direction check, every node of the taut chain must have
   data**; if one has none, the station is refused, with the plain message
   "the river line crosses a gap (NoData) in the DEM at (x, y): the station
   is refused" (`BurnRefusal`, a `ValueError` subclass, from `burn_reach`, a
   `CatchmentError` from `delineate`, like the window's NoData refusal;
   only `BurnRefusal` is turned into a `CatchmentError`, so stage B's
   `except CatchmentError` cannot swallow an ordinary bug). Nodes the cut dropped are
   never burnt or read and are not checked. **A NaN cell is NoData**, as
   in the core (`include/terrain/raster/raster.hpp`, `is_nodata`: NaN
   whatever the sentinel), for the cross-section's choice, this check and
   the extension alike: one mask, a node with data being one that is
   neither NaN nor equal to the sentinel when there is one. `RasterMeta`
   never holds a NaN sentinel (a NaN tag is read as no sentinel), so a
   NaN-gapped float DEM arrives with `nodata` None and its gaps only as NaN
   cells (PR 2's code review, round 2). One check covers every reader:
   the placed node is a taut-chain node, the direction check and the
   descent read only taut-chain nodes, and the extension never takes a
   NoData node (step 5). Stage B's windows give the same taut chain as stage A's
   (the first window already holds the reach and its corridor; stage B's
   extension can run further, but the extension never takes NoData), so the
   refusal comes from stage A. Refused, not rerouted: going round the gap
   (which side, how far) would be a new placement rule. Probe
   (`@architect`, 2026-10-05, a throwaway script, not committed), float32,
   the floor down column 20 above row 20 and down column 22 from row 20,
   `z[20, 21]` NoData, the line down column 20, `U` and corridor 30 m, at
   200 m: on `0588bfa` the chain passes `(20, 21)`, which is the placed
   node; `lowered_max_m` is 33,057 m with −32767 (the chain's last node
   burnt to −32767.04 m) and 3.4e38 with 3.4e38; with data at `(20, 21)`,
   9.95 m. How many of the 140 stations this refuses is not measured; the
   acceptance reports it.
   Corridor, 30 m: a guess from the
   "tens of metres" above; the acceptance re-runs every finding at 15 m and
   60 m.
3. **Direction check.** If the mean elevation over the chain's last tenth is
   more than 2 m above that of its first tenth (reach at least 100 m), the
   line runs against the DEM's slope (2 of 138 measured lines did): the
   result is flagged `direction_disagrees` and the station becomes
   `uncertain` (a refusal would hide it). Lakes and flat reaches (within
   0.5 m: 39 of 138) pass. Pinned at PR 2's green step: a tenth is
   `ceil(n / 10)` nodes of the taut chain of `n` nodes (at least one), raw
   elevations, before the extension of step 5; "at least 100 m" is the taut
   chain's own arc length, not the mapped line's.
4. **Descent.** Along the chain from its first node, `z'[k] = min(z[k],
   z'[k-1] - 0.001 m)`, computed in the array's own dtype: the elevation is
   lowered only where it does not already fall, and the chain falls strictly
   (0.001 m is more than float32's spacing below 16384 m, where it is 0.00098 m,
   and Norway's highest ground is 2469 m, where it is 0.00024 m, so
   `z'[k] < z'[k-1]` holds in float32 as in float64). **The intent** is
   that the flood's flooder of each chain node is the next one down, so water
   that reaches the chain stays on it. **The descent does not by itself
   guarantee that**: the flood (`upstream.hpp`) pushes a node at the higher of
   its own height and the level it was reached at, equal levels first in,
   first out, so a neighbour off the chain that is lower than the next chain
   node, or lies on a lower route to the outlet, can flood the next chain node
   before the chain does. The guarantee is therefore **checked, not assumed**:
   `accumulate`'s `flow_to` says, for each chain node, which neighbour it
   drains to, and the chain holds where `flow_to[k]` is the direction of
   `chain[k+1]` for every `k` in the read stretch (`drains`, in "Sensitivity").
   Where it does not hold, the station is `uncertain` with that cause. A flat
   lake reach is lowered by 0.001 m per node, so by `0.001 m × chain_nodes`
   at most, extension included. For the longest chain the defaults give
   (`reach_up` 1000 m, `U + 100 m` = 600 m downstream at the 500 m map radius,
   and the 500 m end cap: about 2100 m, about 210 nodes at one node per 10 m)
   that is about 0.21 m. A chain that wanders inside its corridor has more
   nodes than its length in tens of metres, so 0.21 m is the usual bound, not
   a hard one; `chain_nodes` and the largest lowering are reported
   (`lowered_nodes`, `lowered_max_m`), and an embankment shows as a large one.
5. **The end of the chain.** A chain whose last node is still lowered (its
   `z'` is below its raw `z`) ends in what the flood sees as a pit: it is
   filled to its spill level and ordered breadth-first, so the counts near the
   end say nothing about the river. So the chain is **extended downstream
   along the raw valley floor** past such an end: from the last node, step to
   the least-elevation 8-neighbour that is not on the chain and is not an
   8-neighbour of any chain node but the last (raw `z`; ties: a straight step
   before a diagonal one, then the smaller (row, column)), apply the descent
   rule to it, and repeat until a node needs no lowering (`z[k] <= z'[k-1] -
   0.001`: the raw ground falls on its own there), or the cap `end_cap` is
   reached, or no neighbour qualifies. A NoData node or one outside the window
   is never a candidate (NoData has no elevation to compare), so the
   extension stops at them only when no other neighbour qualifies. **Why
   straight first** (pinned at PR 2's green step, from `@developer`'s code;
   `test_an_embankment_in_the_last_50_m_extends_the_chain_until_the_ground_falls`
   meets a flat crest of three nodes at one height below the end): (row,
   column) alone turns every tie towards the north-west, whatever the
   river's direction; a straight step first keeps a chain running along a
   lattice axis on that axis, is symmetric left and right of it, and spends
   10 m of `end_cap`, not 14.1 m. It is still a lattice rule: on a flat
   reached by a diagonal step the extension turns 45° onto an axis, and a tie
   between two straight steps still goes to the smaller (row, column). **It
   cannot move the gauge**: the placed node is fixed by step 6 on the chain
   before the extension, and the extension only appends nodes downstream of
   the reach's last one. What a tie decides is which flat nodes the
   extension crosses, and so `end_extended_m`, `end_closed`, and `D` when `U`
   reaches past the reach's last node; on a flat the counts are a filled
   pit's, ordered breadth-first (this step's first sentence), and say nothing
   about the river whichever way the tie goes.
   **The extension keeps the chain
   taut by that choice**, rather than by a second pass after it: a node cut
   after the descent would stay lowered while off the chain. `end_cap` is
   500 m of arc length along the extension, and the metres govern, not a
   node count. **The stopping rule**: a step is taken only if the extension's
   arc length after it is at most `end_cap` (a step of exactly 500 m total is
   taken); the step that would pass 500 m is not taken, and if the last node
   is still lowered there, the cap is hit. On DTM10 that allows at most 35
   diagonal steps (495.0 m; a 36th would make 509.1 m) or 50 straight ones
   (exactly 500 m, summed from 10 m steps, so no rounding decides it). The
   extension then ends at most 500 m of arc length, and so at most 500 m in a
   straight line, from the reach's last node, inside the window, whose margin
   round the reach is `WINDOW_MARGIN_M`, 2000 m. The metres added are reported (`end_extended_m`, 0 when the end
   was not lowered). **A cap hit, or no neighbour with data, inside the
   window, that keeps the chain taut, before the raw ground falls, marks the chain
   `monotone = False`** (the chain's end is not closed), which makes the
   station `uncertain`. A flat floor downstream of
   the gauge, such as a lake, ends here, as it should.
6. **The placed node** is the chain node that the resample point nearest `at`
   chose (after the taut cut, see step 2); its distance from `P` is reported
   (`node_offset_m`, at most about the corridor, and, by step 2's
   cross-section, mostly across the river: along it, at most half a step
   before a taut cut moves it).

The burn is a pure function of the window array, its georeference and the
`Reach`, with no state; it returns a copy (the DEM repository's array is
never written) and the `GaugePath`: the chain's (row, column) array
(extension included), the index of the placed node, the arc length of each node
from it, and `end_extended_m` and `end_closed`. Only one chain of
the network is burnt, the gauge's own reach, so two links never share a
cell (the piracy Lindsay names cannot occur). **The burn can still move
drainage away from the chain**: lowering a chain node that was a
depression's spill point changes where that depression drains, and with it
the drainage of the depression's whole contributing area, however far that
reaches. Only elevations on the chain change; their effect on the flood does
not stop at the chain.

**Known limits, measured not fixed.** (a) The least-elevation node in a
30 m corridor can lie in a parallel gully; the chain then follows that, and
the acceptance's lowered-node counts and the miss analysis show it. (b) The
burn does nothing for a mapped line that is wrong (NEVINA warns REGINE can
be, and ELVIS is derived from N50 at a coarser scale). (c) A lake reach gets
an artificial channel 0.001 m deep per node (about 0.21 m over the longest
usual chain, step 4): the catchment is unchanged except where a chain node
was another depression's spill point (the drainage move just above), but the
count along
a lake's flat floor jumps where the flat's nodes join the chain, which is
what the sensitivity reports as `uncertain`.

### The window loop and the catchment (`catchment.py`)

`CatchmentRequest` gains `reach: Reach | None = None`. `None` keeps
increment 22's behaviour exactly. With a reach, `seed` is the station point
(kept for the report), the pour point is the placed node, and a reach with
`lakes` is refused (`ValueError`: a lake is already the seed). `Catchment`
gains `gauge: GaugeResult | None` (below).

Stage A decides the catchment and every refusal. It is 22's loop, with these
changes: the first window is the reach's bounds (both ends and the corridor)
plus `WINDOW_MARGIN_M`; each window is assembled, burnt, and floods from the
placed node alone; it grows until the placed catchment is clear of the edge,
or refuses as 22 does (NoData, the data's edge, the memory cap). **A refusal
is a refusal of the placed node's catchment, nothing else.** The burn runs
again in each larger window; it depends only on elevations round the reach,
which every window holds.

Stage B reads the area where the sensitivity needs it. The downstream end of
the sensitivity window, `D`, is the first chain node at or past `U`
downstream of the placed node (arc length `>= U`; `U` need not be a multiple
of the node step, and a diagonal step is 14.1 m). Stage B continues from stage A's window: it floods from `D` and grows as
stage A does, so that `D`'s own catchment is clear of the edge. Whether a count
along the chain can be trusted is **not** inferred from `D`'s catchment
containing the placed one (it does so only where the chain drains, which is
what the sensitivity checks): it is decided node by node, by the flag bits of
the count itself, on both sides. If stage B succeeds, `_core.accumulate` runs
in its window and both sides are read. If it refuses (NoData, the edge, the
cap), `accumulate` runs in stage A's window and each side is read as far as
its counts have no flag bit set; `downstream_checked` says `whole` (read to
`U`), `partly` or `none`, and a downstream side not read to `U` makes the
station `uncertain` (next section). Stage B never turns a
catchment into a refusal. In the final window the floods are three: `D`'s
(to size the window), `accumulate`, and the placed node's (the mask the
outline comes from); the earlier windows have one each. The memory cap's
per-node figure grows by the elevation item size (the burnt copy) and 10
bytes (accumulate's four arrays): `item size × 2 + 2 + 10` bytes per node
with a reach, against 22's `item size + 2`.

Pinned at PR 2's green step (`@developer`'s choices, adopted):

- **A reach is in the DEM's CRS.** `delineate` refuses a request with a
  reach whose `seed_crs` is not the DEM's (`CatchmentError`), as the CLI
  refuses a river file in another CRS: the burn reads the reach in the
  window's coordinates, and a library caller gets the same refusal as the
  command.
- **Stage B starts from stage A's last window's bounds**, so its first
  window is that window again; `Catchment.windows` lists it once, not twice
  (the code drops stage B's first entry only when its bounds equal stage A's
  last).
- **`GaugeResult` has no `reach_fork`.** The fork is the placement's
  (`Placement.reach_fork`), and `delineate` never sees the placement; the
  command that has both writes it into the catchment file from the
  placement (`catchment --rivers` now, `station-catchments` in PR 4).
- **The `direction` cause is joined outside `assess`.** `Sensitivity.causes`
  holds the four causes `assess` can see; `direction` comes from
  `GaugeResult.direction_ok`, and `catchment --rivers` appends it to the
  causes it writes and to its verdict. `Sensitivity.well_posed` therefore
  does not reflect the direction check; the station's verdict is "no cause
  in the joined list". **PR 4, as the second reader, moves the join into
  `catchment.py`** (one `causes` on `GaugeResult`, read by both commands),
  as the GeoJSON writer moves when it gets its second caller.

```python
@dataclass(frozen=True, slots=True)
class GaugeResult:
    node: tuple[float, float]  # the placed node, in the DEM's CRS
    chain: tuple[tuple[float, float], ...]  # the burnt chain's nodes, DEM CRS, downstream
    node_offset_m: float       # from the mapped position P to the node
    chain_nodes: int           # nodes of the burnt chain
    lowered_nodes: int
    lowered_max_m: float
    direction_ok: bool
    end_extended_m: float      # the chain's end run along the raw valley floor
    end_closed: bool           # False: cap, NoData or window edge before the ground fell
    downstream_checked: Literal["whole", "partly", "none"]
    sensitivity: Sensitivity   # next section
```

### Sensitivity: is the area well defined at this gauge? (`sensitivity.py`, pure)

A gauge's coordinates are uncertain by roughly their distance from the river,
and the river's catchment area changes along the river. If moving the gauge
within its own uncertainty changes the area by more than the disagreement the
comparison tolerates, then a `match` or a `miss` would be decided by where
along the river we put the node, not by how well we delineate. Ola's
"saddle point problems" are this: a confluence just below the gauge, or a
flat valley floor or lake where the drainage line is arbitrary. The check
reads the area along the burnt chain at the placed node, never at any other
choice of node, and never uses NVE's polygon or area.

`sensitivity.assess(count, reach_bits, flow_to, path, cell_area_km2,
uncertainty, reach_down_m)` (arrays in, a frozen `Sensitivity` out; no DEM, no
shapely):

1. **Samples**: every chain node whose arc length from the placed node lies in
   `[-U, +U]`, one per node (10 m apart on DTM10; the chain advances a node
   at a time). The area at a sample is `count × cell area`, read from the
   accumulation that stage B ran: one gather, no further flood.
2. **Which counts to trust, by their flags on both sides.** A sample is
   trusted when its `reach_bits` is zero: its count is then the node's whole
   catchment. Containment is not used (the chain upstream of the placed node
   lies inside the placed catchment only where it drains, which is checked
   below, not assumed). Each side is read outward from the placed node up to
   the first untrusted sample; the distances read are `checked_up_m` and
   `checked_down_m` (the arc length of the farthest trusted sample of that
   unbroken run, so at most `U`). **Downstream, the read must reach `U`**:
   every sample in `(0, U]` is trusted **and the chain has a node at arc
   length `>= U`** (that node is `D`, stage B's; its own count is not a
   sample and its flags do not matter when it lies past `U`; a node exactly
   at `U`, within the 1e-6 m slack, is both `D` and the last sample in
   `(0, U]`, and must be trusted like the others: a flag there gives
   `downstream_unread`). This is not `checked_down_m >= U`:
   where `U` is not a multiple of the step, the farthest sample lies short of
   `U` (with `U` = 35 m on straight steps, the samples are at 10, 20 and
   30 m), and the strict test would mark every such station unread, every
   diagonal river at the 30 m floor among them. The mapped river must also
   reach it (`reach_down_m >= U`, from "Placing
   the gauge"; the extension of "Following the river" is not mapped river). A
   side that stops short is a gauge whose neighbourhood was not looked at, and
   the station is **`uncertain`** (cause `downstream_unread`), not scored.
   Upstream, a river that starts within `U` ends the read without penalty (the
   area above a source is the hillside's), but a flagged count upstream stops
   it and counts as unread. **What "counts as unread" means** (provisional:
   `@architect`'s reading after PR 2's green step, which gave it no
   consequence; question 6 for Ola, default yes): a
   flagged sample within `U` upstream makes **`drains` false**, so the
   station is `uncertain` with the cause `chain_not_draining`. The reason is
   a proof, not a guess: `accumulate` ORs each node's bits into its
   flooder's, so a node whose flooder chain passes through the placed node
   has bits only if the placed node has them; the placed node's are clear
   (stage A refuses otherwise), so the flagged node does not drain through
   the placed node, nor (the same argument) into the read stretch next to
   it. The river within `U` upstream then lies outside the catchment, which
   is exactly the position sensitivity this check exists to report.
   `checked_up_m` still stops before the flagged sample, and the swing is
   still taken over the trusted run.
3. **Does the chain drain along itself?** `drains` is true when, for each
   consecutive pair of read samples, `flow_to` of the upper one is the
   direction of the lower one: each chain node drains into the next, node by
   node, in the flood the counts come from. The descent of the burn intends
   this and does not guarantee it ("Following the river", step 4). `monotone`
   is true when the counts increase strictly (`count[k+1] > count[k]`, which
   `drains` implies: the lower node's count includes the upper one's and
   itself) and the chain's end is closed (`end_closed`, step 5 there); it is
   the cheap cross-check, and the two failures are reported as separate causes
   (`chain_not_draining`, `chain_end_open`). Counts that are equal or fall
   along the read stretch give `chain_not_draining`, whether or not the end
   is closed (pinned at PR 2's green step: by the implication just stated
   they cannot occur where `drains` holds, so they are a drainage failure;
   an open end adds `chain_end_open` beside it).
4. **The measures**: `A0`, the area at the placed node; `A_up` at the farthest
   trusted sample upstream and `A_down` at the farthest downstream; the
   **swing** `max(A0 - A_up, A_down - A0) / A0`, **one-sided** on purpose: the
   question is how far the area at the true position can be from the area at
   the placed node, which is the larger of the two ends. Adding both ends would
   count a move up and a move down at once, which one position cannot make, and
   would hold the station to a bar stricter than the `match` bar it is
   measured against; and
   the **largest step**, the biggest area increase between two neighbouring
   samples, with its position in metres from the placed node (negative
   upstream), which names the confluence when there is one.
5. **The rule.** The station is **well posed** when the swing is at most
   `SWING_MAX = 0.05`, `drains` and `monotone` hold, and the downstream side
   was read to `U`. `SWING_MAX` is the `match` class's own bar (95 % overlap,
   below): a position error inside the gauge's uncertainty is allowed to cost
   no more than the comparison allows. A station that is not well posed is
   **`uncertain`**: reported with its agreement numbers, not scored `match` or
   `miss`, and counted apart, by its causes (`swing`, `downstream_unread`,
   `chain_not_draining`, `chain_end_open`, and, from the burn, `direction`).
   Ruled by Ola, 2026-10-04 (Question 1, the default): `SWING_MAX = 0.05`,
   one-sided, `U` as in "Placing the gauge".

A confluence 40 m below the gauge with a tributary of 20 % of the placed
area gives a swing of at least 0.2 and a step at about +40 m (arithmetic, not a measurement). A flat lake
floor gives a step wherever the flat's nodes join the chain. A river in a
well-defined valley gains area smoothly: over 2 × 30 m, expected (an estimate
not yet measured) a fraction of a per cent of a catchment of 100 km². Small
catchments should more often be uncertain (a 0.44 km² catchment may gain
several per cent over 60 m), and that would be true of them, not an artefact: the acceptance reports the share
uncertain per size band.

**Cost and novelty.** One gather from arrays that already exist. It is a
diagnostic, not a method claimed (see Novelty).

### Residual inflow, later: what 29 keeps open

**Not built here.** Ola plans to compute the residual inflow to rivers: for
two gauges A above B on one river, the catchment of B without that of A, the
area the river gathers between them. That needs the gauges to *nest* (A's
catchment lies inside B's) and the difference to be only what lies between
their physical positions. The check that they nest needs **no NVE polygon**:
it is a property of our own two delineations, and a failure is a finding in
itself.

Why 29's placement suits it: each gauge sits where it physically is, so the
difference is the area between two real positions. A rule that moves each
gauge to where the area is largest moves each by an unrelated distance
downstream, and the difference then contains the moves.

What 29 must not do: delineate each station on its own burn and then
subtract. Two stations' burns differ upstream of both (each chain's descent
starts from its own first node), so their catchments need not nest exactly. The
later increment delineates **all gauges of one river from one burn and one
flood**: the counts then come from one tree of flooders, and nesting holds
where A's chain of flooders passes through B, which is **checkable node by
node in the shared flood** (`flow_to` along the joint chain, as the sensitivity
checks `drains`; it is not guaranteed by the burn's descent). Where it holds,
the residual is exactly `count_B - count_A` nodes; where it does not, the pair
is reported as not nested. The sketch, not built:

```python
# Reach.at becomes a tuple: one chain, several gauges (29 passes one)
def delineate_river(request: RiverRequest, repository: DemRepository) -> tuple[Catchment, ...]
def residual(upstream: Catchment, downstream: Catchment) -> Residual
    # polygon difference, its area, and the nesting check: the part of the
    # upstream outline outside the downstream one, in nodes (0 when nested)
```

What 29 does so that this is not blocked:

- The placed position is an index into the chain, not a coordinate, and
  the chain and its flood are the unit, not the station. In 29,
  `GaugePath.placed` is one `int` (PR 2's red step pins it so); the residual
  increment widens it to an index array, a change local to `burn_reach` and
  `sensitivity.assess` (amended 2026-10-05 from "take an index array" here,
  which no 29 caller needs).
- `Placement` and `StationResult` record the river (`elvid`), the segment
  (`objectid`) and the position along it, so gauges can be grouped and
  ordered along a river without geometry.
- `place` builds a chain from `P` alone; a joint chain is the same
  construction from the most upstream gauge's `reach_up` to the most
  downstream gauge's end, with no change to the burn.
- Every station's fine catchment polygon is written (22's file), so
  `residual` can run from the output files.
- The sensitivity window of each gauge is recorded, so a residual between two
  gauges closer than the sum of their windows can be flagged as uncertain.

### The station set (`fetch/nve.py`, `io/station_set.py`, `io/rivers.py`)

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
  38 for their polygons, each with its explicit `outFields` allow-list (see
  "Data use"), 40 station numbers per `where ... in (...)` query
  (measured to work; under the service's 2000-record cap and URL limits),
  with `outSR=25833&f=geojson`, through `RangeClient.get_text` (23a-2's
  retries and refusals; `fetch/http.py` stays the only module importing
  `urllib`, rule F11);
- queries ELVIS layer 2 once per station, by envelope: the station point
  plus `map_radius + reach_up + 500 m` (2 km) each way, `outFields=objectid,
  objekttype,strekninglnr,elvid,vassdragsnr,elvenavn,vatnlnr`, `outSR=25833&f=geojson`
  (an explicit allow-list, see "Data use"; `vatnlnr` is there for the lake
  rule of "Placing the gauge", and `objekttype` is written as served).
  The service returns whole features that meet the envelope, so a chain
  longer than the envelope is complete where it meets it. Segments seen from
  several stations are kept once, by `objectid`. **Unverified:** only 1 km
  squares were queried (28 s for 140); at 4 km the answer is larger, and a
  reply flagged `exceededTransferLimit` is refused, naming the station, so a
  truncated river never reaches `place`;
- refuses, naming the station, when a listed station has no point or no
  polygon; keeps the newest polygon per station (above) and records its
  update date and how many versions there were;
- writes `stations.geojson` (one Point feature per station, properties
  `station`, `name` = layer 0's `stasjonnavn` (not the list's),
  `series` = `["1001.0"]`, `nve_area_km2` from layer 0,
  `hrd_start_daily`, `watercourse` = layer 0's `vassdragsnr`, `river` = the
  first name of `elvenavnhierarki`; no other layer 0 field is copied),
  `reference.geojson` (one feature per station, the polygon, `station`,
  `reference_area_km2`, `reference_updated`, `versions`), `rivers.geojson`
  (one LineString per segment, the seven fields above; copies are kept as
  served and dropped by `read_segments`), all with a `crs`
  member, `NOTICE.txt` (the catalogue's credit, "Kilde: NVE", as 23a-2's
  `notice` does), and `manifest.json` (the query URLs, the fetch time in UTC,
  each file's sha256). Deterministic order: the list file's, then `objectid`.
- Pure parts (the query URLs, choosing the newest version, de-duplicating
  segments, building the files' content) are functions with no network,
  tested on canned answers as `fetch_fixtures.py` does; the network call is
  injected.

**`io/station_set.py`** reads the stations and references back:
`read_stations(path) -> (tuple[Station, ...], crs)`, `read_references(path)
-> (Mapping[str, Polygon | MultiPolygon], crs)`; **`io/rivers.py`**
`read_segments(path) -> (tuple[RiverSegment, ...], crs, dropped)`, where
`dropped` is the count of exact copies removed by the pure
`drop_copies(segments) -> (kept, dropped)`. All refuse a file
without a `crs` member (the skill's rule; NVE's files always have one),
duplicate station numbers or segment ids, and the wrong geometry type (a
segment that is a MultiLineString is refused, naming its `objectid`).
`Station` is a frozen Pydantic model: `station: str` (pattern
`^\d+\.\d+\.\d+$`), `name`, `x`, `y`, `series: tuple[str, ...]`,
`nve_area_km2: float | None`, `watercourse: str | None`, `river: str | None`.
Any GeoJSON of points with a `station` property works, so a user can bring
their own list; without `watercourse` and `river` every candidate line is in
the lowest tier of "Placing the gauge".

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
  perimeter, in metres: the area between the two outlines divided by the
  length of NVE's. Exact when one outline lies a constant distance outside
  the other (a square grown by one cell on every side: one cell, up to the
  corner term), and size-independent, which is what NRFA's warning about
  small catchments asks for. **A shift along one axis scores half the shift**,
  because two of a square's four sides do not move (a square shifted one
  cell along x: half a cell); the measure is the average over the whole
  divide, not the largest displacement.

`classify` takes the agreement numbers and the sensitivity and returns a
class and, for a match, the test that passed. Precedence: `refused`, then
`uncertain`, then the rest. (Default, @architect, 2026-10-04.)

| class | rule | counts as |
|---|---|---|
| `refused` | `delineate` refused (NoData, the data's edge, memory cap, tiles on two grids, no river line within the map radius, a NoData node on the river's chain), with its message and `refusal_cause` | reported, not a failure; `mixed_grid` counted apart as known refusals |
| `uncertain` | not well posed: swing over 5 %, the downstream side not read to `U` (a flag stopped it, or the mapped river ends sooner), the chain not draining along itself or its end open; or the line runs against the DEM's slope | reported with its agreement numbers, not scored |
| `match` | both overlaps ≥ 95 % (`match_by = "overlap"`), **or** mean divide offset ≤ 3 cells, 30 m on DTM10 (`match_by = "offset"`; the overlap test is tried first) | pass |
| `close` | both overlaps ≥ 80 % | finding |
| `miss` | anything else | finding |

Bygdin, the one case measured so far, would be `match` by overlap (99.12 %
and 99.33 %). The offset test is there for the 12 catchments under 10 km²,
where a one-cell disagreement along the whole divide is several per cent of
the nodes; `match_by` is in every row, so the share that passes only by the
offset test is visible and is not mistaken for the overlap bar.

`summarise(rows) -> Summary`: counts per class and per `match_by`; for the area
ratio and both overlaps, the minimum, 10th, 25th, 50th, 75th, 90th percentile
and maximum, over the scored stations (`match`, `close`, `miss`) and per size
band (under 10, 10-100, 100-1000, over 1000 km²) and per tile count (1, 2,
3-4, 5+); the share `uncertain` per size band, and the causes of `uncertain`, each
counted (a station can have several): `swing` (split by its largest step's
position: a confluence step, or a flat floor or lake), `downstream_unread`,
`chain_not_draining`, `chain_end_open`, `direction`; the refusals counted by
`refusal_cause`, with the `mixed_grid` ones apart as **known refusals**: one
summary line, "N stations refused because their windows select tiles on
two different grids, which rasputin does not combine, and neither grid
covers the window alone (known refusals, not failures)", and their station numbers. The line is general because the
refusal is (below): on DTM10 its only cause is the half-cell shift, and the
acceptance README, not the summary, says that a later increment resamples
those tiles. Deterministic JSON.

### The batch (`catchment_batch.py`)

```python
class BatchRequest(BaseModel):      # frozen
    map_radius: float = 500.0
    reach_up: float = 1000.0
    outline_tolerance: float | None = None
    only: tuple[str, ...] = ()      # station numbers; empty is all

class BatchSink(Protocol):
    def catchment(self, station: Station, result: Catchment) -> None: ...
    def row(self, row: StationResult) -> None: ...

async def run_batch(request: BatchRequest, repository: DemRepository,
                    stations: Sequence[Station], stations_crs: str,
                    segments: Sequence[RiverSegment], segments_crs: str,
                    references: Mapping[str, BaseGeometry] | None,
                    sink: BatchSink) -> Summary
```

- **One station at a time**, each `delineate` in `asyncio.to_thread`: a
  flood of a large catchment takes gigabytes, and two at once would race for
  the memory cap. Sequential order is the file's, so the output is
  deterministic. A GUI or API worker awaits it and can cancel between
  stations. Parallel stations are a later option.
- A refusal (`CatchmentError`) becomes a `refused` row and the batch goes
  on; any other exception stops it (a bug is not a data refusal). The row's
  `refusal_cause` is `mixed_grid` when the refusal is the mosaic's
  mixed-lattice one, `no_river` for `place` returning `None`, and `other`
  for the rest (NoData, the data's edge, the memory cap: their message says
  which). **`mixed_grid` is told by type, not by message**: `mosaic.py`
  raises a new `MixedGridError(MosaicError)`, carrying two tiles' names,
  where `_covering` today raises the plain `MosaicError` of `_mixed(...)`
  (same message, so increment 15a's suite is unchanged), and
  `catchment._plan` maps it to a new `MixedGridRefusal(CatchmentError)` with
  the same `tiles`. **What it covers**: every refusal `_mixed` words, so tiles
  that differ in CRS, spacing, registration or NoData value, as well as a
  lattice offset; `mixed_grid` therefore means "tiles on two grids", not
  "half-cell-shifted tiles". On DTM10 only the offset occurs (all 254 tiles
  have 10 m spacing, NoData −32767, EPSG:25833 and PixelIsArea, read with
  `tifffile`, 2026-10-04). **Which two tiles**: the ones `_covering` passes to
  `_mixed`, `chosen[0].selected[0]` and `chosen[1].selected[0]`, the first
  selected tile of each of the first two grids: one tile from each grid,
  named so the cause can be found, and not necessarily a tile the catchment
  crosses (the window's margin can select a tile the catchment never
  reaches). The row names those two tiles as such. A station
  `place` returns `None` for is `refused` ("no mapped river line within 500 m
  of the station") until the fallback of "The fallback" lands.
- `StationResult` (frozen): station, name, class, `match_by`, refusal message,
  `refusal_cause` and, for `mixed_grid`, one tile from each grid;
  the placement (`placed_on`, station-to-line distance, `U`, `lake`,
  `confluence_near`, `elvid`, segment `objectid`, `reach_up_m`,
  `reach_down_m`); the gauge numbers (`node_offset_m`, `lowered_nodes`,
  `lowered_max_m`, `direction_ok`, `downstream_checked`); the sensitivity
  (`area_up`, `A0`, `area_down`, swing, largest step and where, `checked_up_m`,
  `checked_down_m`, `drains`, `monotone`, the causes);
  the burn's `end_extended_m` and `end_closed`, and `reach_fork`;
  nodes, fine and reduced area, NVE's polygon area and the station layer's
  area, the agreement numbers, tile count, windows, seconds.
- Without `--reference`, no agreement and no class beyond `refused` and
  `uncertain`: the batch just makes the catchments (the product for any
  list of stations).
- No paths: the sink is how files get written; `cli.py` passes a directory
  sink, tests pass a list.

**`rasputin station-catchments`** (named for what it makes, a catchment per
station; the earlier name `catchments` differs from `catchment` by one letter):

```
rasputin station-catchments --dem PATH [--dem PATH ...] --stations FILE
                    --rivers FILE [--reference FILE] [--map-radius METRES]
                    [--reach-up METRES] [--only ID ...]
                    [--outline-tolerance METRES] --out-dir DIR [--out-parent DIR]
```

writes `DIR/<station>.geojson` (22's catchment file, plus the station's
number, name, series, the placement and the sensitivity; `mesh --domain`
reads it), `DIR/results.csv`, `DIR/summary.json`, and one stderr line per
station (number, name, class, the three numbers) and the summary at the end.

**`rasputin catchment`** gains `--rivers FILE` (with `--map-radius` and
`--reach-up`): the seed is placed on the nearest line as above, with no
watercourse number, and one stderr line says where: "placed on the river line
31 m from the station (river Nea, line 8841), moved 6 m onto the DEM's valley
floor; 119.0 km² drain through it, and the area changes by 0.3 % within
30 m up and down the river: well defined". The 6 m is `node_offset_m`,
which step 2's cross-section makes mostly a move across the river, so the
sentence no longer says "farther on" (amended 2026-10-05). A station given
in another CRS (`--seed-crs`, WGS84 by default as in 22) is transformed into
the river file's CRS before `place`, and the request's seed is passed in
that CRS; the river file's CRS must be the DEM's, or it is refused (the
burn takes the reach in the window's coordinates).

**The GeoJSON writer moves** in PR 4, where `station-catchments` becomes
the second command to write a catchment file (PR 3 writes none), from
`cli.py` to `io/geojson.py`,
`catchment_geojson(polygon, crs, properties) -> bytes`, no path, as
`project_structure.md` already recommends ("The catchment GeoJSON writer is
in `cli.py`"); both commands call it, and that paragraph is replaced by the
module's entry.

### New and changed files

Estimates are production lines under CLAUDE.md §2's counting. 22's
`catchment.py` ran 44 % over its estimate, so each PR's second figure adds
that margin. Round 2 changed them: PR 3 takes `io/rivers.py` and
`RiverSegment` from PR 2 and gains the `User-Agent` header and the field
allow-lists (PR 3 from about 300 to 380); PR 2 gains the chain end, loop cuts,
forks, the `drains` check and the causes, and loses the river reader (375 to
415); PR 1 gains `flow_to` (150 to 155). Round 3 adds the taut pass to
`burn.py` (inside its 105: it replaces the loop cut, which it subsumes) and
the known-refusal cause to PR 4 (295 to 310). After PR 3's code review,
round 1, the moved GeoJSON writer goes from PR 3 to PR 4, which is the
first to need it (PR 3 380 to 355, PR 4 310 to 335). PR 2's green step
(`91857df`) came to 556 net lines, 34 % over its 415 and inside the
44 % margin (598): `burn.py` took 62 more than estimated (the
cross-section of step 2 and the extension's candidate rules) and `cli.py` 55
more (the river file's reading and CRS check, the options' checks and help,
the placement report). With the second red amendment's green
(`0588bfa`) it is 558 net (578 added, 20 removed; `burn.py` 167,
`gauge.py` 117, `sensitivity.py` 85, `catchment.py` 94, `cli.py` 95);
with the NoData fix of code review round 1 (green `96b3881`) it is 568 net
(588 added, 20 removed; `burn.py` 174, `catchment.py` 97, the others as
before); with the NaN fix of code review round 2 (green `193079d`) it is
569 net (589 added, 20 removed; `burn.py` 175, the others as before).

| File | What | Estimate |
|---|---|---|
| `include/terrain/hydrology/flood.hpp` | the shared flood, moved out of `upstream.hpp` | 45 |
| `include/terrain/hydrology/upstream.hpp` | uses it | 10 (35 removed) |
| `include/terrain/hydrology/accumulate.hpp` | `accumulate`, `AccumulateOutcome`, the bits | 60 |
| `bindings/core.cpp`, `_core.pyi` | `accumulate`, with `flow_to` | 40 |
| **PR 1, accumulation** | | **about 155** |
| `data/nve_hrd_2025.csv` | the list (data, not counted) | 0 |
| `sources.py` | `StationSource`, the `nve-hrd` entry | 30 |
| `fetch/http.py` | the `User-Agent: rasputin/<version>` header | 5 |
| `fetch/nve.py` | queries with the field allow-lists, newest version, segments, files, manifest | 160 |
| `io/station_set.py` | readers, `Station` | 50 |
| `io/rivers.py` | `RiverSegment`, `read_segments`, the `kind` mapping, exact copies | 65 |
| `cli.py` | `fetch-stations` | 45 |
| `NOTICE.md`, `project_structure.md` | NVE's credit; the new modules | docs |
| **PR 3, the stations and rivers** | | **about 355 (511 with the margin)** |
| `gauge.py` | `Gauge`, `Placement`, `place`, forks, `Reach` | 100 (117 at green) |
| `burn.py` | valley floor, taut pass, descent, chain end, direction, `GaugePath` | 105 (175 at green) |
| `sensitivity.py` | `assess`, `Sensitivity`, `drains`, the causes | 75 (85 at green) |
| `catchment.py` | request field, stages A and B, burn per window, `GaugeResult` | 95 (97 at green) |
| `cli.py` | `catchment --rivers` and the placement line | 40 (95 at green) |
| `docs/benchmarks/<date>/nve-placement/render.py` | the placement figures (evidence script, not counted, Ola's ruling; "Placement figures") | (130) |
| **PR 2, the gauge on the river** (needs PRs 1 and 3) | | **about 415 (598); 569 net (589 added) at green `193079d`** |
| `reference.py` | agreement, classes, `match_by`, summary | 115 |
| `catchment_batch.py` | `BatchRequest`, `BatchSink`, `run_batch`, `StationResult`, `refusal_cause` | 105 |
| `mosaic.py`, `catchment.py` | `MixedGridError`, `MixedGridRefusal` (round 3, Ola's ruling on counting) | 10 |
| `cli.py` | `station-catchments`, the directory sink | 80 |
| `io/geojson.py` | the moved writer (moved from PR 3 after its code review, round 1) | 25 (cli.py −25) |
| **PR 4, the batch and the comparison** | | **about 335 (482)** |
| `catchment.py`, `cli.py`, `catchment_batch.py` | the fallback (below) | 80 |
| **PR 5, the fallback** | | **about 80 (115)** |

Every PR is under the 700 ceiling of CLAUDE.md §2, even with the margin. No
module needs a new dependency.

### The PR split

Each PR answers a question on its own, and the seams are where a module's
inputs are files or arrays, not another PR's types.

- **PR 1, accumulation.** C++ and its binding: "how many nodes drain through
  each node". One build round; the oracle against `upstream` is its whole
  suite.
- **PR 3, the stations and rivers.** `fetch-stations`, the readers
  (`io/station_set.py`, and `io/rivers.py` with `RiverSegment`, the `kind`
  mapping and the dropping of exact copies), the packaged list. Answers "give
  me NVE's 140 stations, rivers and polygons, with a manifest". The seam with
  PRs 2 and 4 is the three files it writes and the models it reads them into
  (so their tests use hand-written files). Python only; independent of PRs 1
  and 2, so its red step can be written while PR 1 is in review.
- **PR 2, the gauge on the river.** Place, burn, sensitivity and the window
  stages, behind `rasputin catchment --rivers FILE` on a river file the user
  has. Answers "the catchment of this gauge, and is it well defined". Python
  only; **needs PR 1 (`accumulate`) and PR 3 (`RiverSegment`,
  `read_segments`)**, so the merge order is 1, 3, 2, 4, 5. `gauge.place`
  takes its own small `Gauge` (point, watercourse number, river name) rather
  than PR 3's `Station`, so the two meet only at `RiverSegment`. Ends with
  the placement figures (four stations, before and after, PNG; "Placement
  figures").
- **PR 4, the batch and the comparison.** `reference.py`,
  `catchment_batch.py`, `station-catchments`; then the acceptance run.
  Needs PRs 1 to 3.
- **PR 5, the fallback.** Nearest stream for a station with no river line
  near (below). Small, and last because it is the only part Ola may rule out.

Neither touches refine or mesh code, so `tools/bench.py`'s benchmark and
scaling sweep do not apply (`docs/increments/README.md`, "Acceptance").

### Placement figures, before and after (PR 2)

Ola asked to see how the placement works, on a few real stations, as PNG
maps before and after. **A deliverable of PR 2**, the first PR in which a
gauge is placed, burnt and delineated (`catchment --rivers`); it runs on
PR 3's fetched files, which merge before PR 2 (order 1, 3, 2).

**The script.** `docs/benchmarks/<date>/nve-placement/render.py`, beside the
figures and their README, as Bygdin's `render.py` sits beside its pictures
(`docs/benchmarks/2026-09-29/bygdin-landcover/`). It composes the package's
public functions and adds no logic of its own to the placement:
`read_stations`, `read_segments`, `read_references`, `gauge.place`,
`catchment.delineate` (whose `GaugeResult.chain` gives the burnt chain), and,
for the flow paths, one window read round the reach, `burn.burn_reach` on it
and `accumulate` on the raw and the burnt arrays. **Plotting with matplotlib**,
in this evidence script only, and not added to the package's dependencies or
to any extra: ruled by Ola, 2026-10-04 ("map questions: 1: No, 2: Yes.",
question 2; see "Ola's rulings"). Increment 6's ruling keeps it out of the
package, and other evidence scripts already use it outside the package
(`docs/benchmarks/2026-10-02/15c-2-acceptance/run_geo.py`'s independent
check). The repository's own picture route, `tin_engine.viz.svg`, draws a
triangulation as SVG and writes no PNG; VTK (the `viewer` extra, which
Bygdin's PNGs used) draws a 3-D scene, not a map with lines and a legend.
The README records the matplotlib version, and `run.sh` the commands.

**Per station, two PNGs** (`<station>_before.png`, `<station>_after.png`),
same map extent: the reach's bounds plus 300 m.

- *Before*: a hillshade of the raw DEM; the DEM's flow paths, the nodes whose
  raw `accumulate` count is at least 0.1 km² (1000 nodes on DTM10), drawn
  darker as the count grows; NVE's mapped river line (the chosen segment and
  the rest of its reach); NVE's station point; the mapped position `P`.
- *After*: the same hillshade; the flow paths of the **burnt** array; the
  burnt chain (lowered nodes marked), the placed node, the sensitivity
  window (`U` up and down); and an inset at catchment scale with our fine
  outline over NVE's polygon, with the class numbers of the station
  (overlaps, area ratio, swing, causes if `uncertain`).

**Which stations: four cases**, chosen by rule, not by eye, and the README
gives the numbers that chose each:

1. **Placed straight onto the DEM's flow path**: well posed, no node lowered
   (`lowered_nodes == 0`), and of those the one with the smallest NVE polygon
   area (the reference layer's area, since the survey computes no catchment;
   quick to run).
2. **The line burnt in**: well posed, with the largest `lowered_max_m`.
3. **Marked uncertain**: cause `swing`, with a confluence step; from the 19
   stations with `confluence_near`, the first in list order that is.
4. **A lake gauge** (`lake`), the first in list order that is well posed.

Cases 1 and 2 are proposed by a survey mode of the script (`--survey`: place
and burn every station on its first window, no catchment, so seconds per
station). The survey's "well posed" is only the first window's (no
`drains` failure, a closed end, the swing within its bar there), so a
proposed case is **confirmed by the full run** (`catchment --rivers`): if
the full run makes it uncertain or refuses it, the next candidate in the
survey's order is run, until one is confirmed. Cases 3 and 4 are chosen by
running `catchment --rivers` on the candidates in order until one
qualifies. A case no station qualifies for is reported as such, not filled
by a near miss; the README lists every candidate run and why each one that
was passed over failed.

**Lines.** About 130 (the survey, the two maps, the inset). **Not counted
toward PR 2's 700-line ceiling, by Ola's ruling** of 2026-10-04 ("map
questions: 1: No, 2: Yes.", question 1; see "Ola's rulings"): it is an
evidence script, outside the package. PR 2 stays at about 415 (598 with the
44 % margin). **Run by `@perf`** after PR 2 is green, as it owns evidence.

### The fallback: nearest stream (Jenson), PR 5

**When.** Only where the river line cannot place the gauge: no segment within
`map_radius` (one of the 140, Femundsenden, a lake gauge), or no river file
(a station list a user brings without one). **And as a comparison** in the
acceptance: the same rule run on every station that ends `miss` or
`uncertain`, so that what the mapped river buys is measured. It is never used
to repair a station the river rule placed.

**The rule** is Jenson 1991's (above): the nearest DEM node within `R` (default
250 m: all 42 stations outside their NVE polygon lie within 233 m of it) whose
accumulation count is at least `stream_min_km2` (default 0.1 km²; a guess,
and the comparison run measures it). `accumulate` runs over a fixed window,
the disc's bounds plus 1 km: a count with a flag bit set is a lower bound, so
reaching the threshold is conclusive, and a flagged node under it is treated
as not a stream (a stated limit, only at the window's edge). No node
qualifies: `CatchmentError` ("no DEM stream within R m"). The found node is a
plain pour point of 22's `delineate`, then. There is no chain, so **no
sensitivity is computed**: the row says `placed_by = "nearest stream"`,
the class is decided by the agreement numbers alone and the summary counts
these apart, never among the stations the sensitivity cleared.

**Not a rule that looks at area.** The threshold is a floor on the count, not
a search for the largest: the nearest qualifying node wins, as Lindsay et al.
2008's reading of Jenson says it should.

## The red suites

Lean, as 22's were: no throwaway implementation. **The invariant-critical
suite** (mutation testing required, `docs/increments/README.md`, "Cost
constraints") is `accumulate`'s exact oracle, because every sensitivity and
every station result rests on it. **The mutation target is `accumulate`'s own
code**: its `on_reach` callback (the flooder, which is `flow_to`, and order records) and its
reverse pass (the count addition and the flag OR). `flood.hpp` is the code
moved out of `upstream.hpp` and is covered by `upstream`'s existing suite,
unchanged and run after the move; it gets no mutation round of its own. The
Python suites below are not invariant-critical.

**PR 1, C++ (`tests/cpp/unit/test_hydrology_accumulate.cpp`):**

- *The oracle*: on a few hundred small random DEMs (with pits, flats, equal
  heights, NoData holes and NoData borders), for every node with data,
  `count == upstream(z, {node}).nodes_in`, and the two bits of `reach` equal
  that call's `touches_edge` and `touches_nodata`; NoData nodes count 0 and
  carry no bits.
- *Conservation*: the outlets' counts sum to the nodes with data; every
  count is at least 1 on data.
- *By hand*: a V-valley (the outlet counts the whole valley); a single
  cell; a flat plateau draining over one rim node; a 1 × n strip; a node
  two cells from the edge (bit 0 clear) beside one at one cell (set).
- *`flow_to`*: for every node with data, 255 exactly when it is an outlet,
  else a neighbour with data; `count[c] == 1 + sum(count[d])` over the `d`
  with `flow_to[d]` pointing at `c`; the pointed-at neighbour `n` has
  `upstream(z, {n}).mask[d] == 1`; following `flow_to` from any node reaches
  an outlet (no cycle); a two-node ridge with equal heights (the first-in
  rule decides, the same way twice).
- *Determinism*: twice gives equal arrays. *Refusal*: the size check, by a
  geometry stub reporting 2^32 nodes (no allocation).
- `upstream`'s existing suite, unchanged, after the move to `flood.hpp`.

**PR 1, Python:** `test_core_accumulate.py`: shapes, dtypes (`uint32`,
`uint8`; `flow_to` is `uint8`), ownership (the arrays outlive the view), the oracle on two DEMs
through the binding.

**PR 2, Python** (hand-built DEMs and lines; no network; `RiverSegment`s are
built directly or read through PR 3's `read_segments`):

- `test_gauge.py`: tiers (own watercourse number beats a nearer line of
  another river; prefix and name tiers; `any`); the nearest of the best tier;
  exact copies do not make a fork, a real fork (two different continuations,
  as at `19.80.0`) stops the chain, flags `reach_fork` and reports the
  metres, upstream as downstream; `kind` of the chosen line reaches the
  `lake` flag;
  `P` as a perpendicular foot and as an end vertex; `U` at `d` = 10 m (30),
  100 m (100) and the cap; the chain joins ends within 1 m and stops at a gap
  of 2 m; `reach_up` cuts at the metre; a river that ends gives a shorter
  reach and the metres are reported; `lake` and `confluence_near` at the
  100 m boundary, measured from `P` (a second branch 90 m from `P` and 110 m from
  the station point counts); no line within the radius gives `None`; a segment file in
  another CRS is refused.
- `test_burn.py`: an embankment across a valley: `upstream` on the raw array
  from a node below it counts the strip below the embankment only, and on
  the burnt array counts the whole valley; the burnt array is `<=` the input
  everywhere, equals it off the chain, and on the chain equals the input
  wherever the input is already below the carried level minus 0.001 m (computed
  in the array's dtype); the chain falls strictly, `z'[k+1] < z'[k]`, in
  float32 and in float64 (a flat chain at 3000 m, where float32's spacing is
  0.00024 m); an embankment in the **last 50 m** of the chain: the end is
  lowered, the chain is extended along the raw valley floor until the raw
  ground falls (`end_extended_m` > 0, the added nodes are 8-connected and
  fall), and the burnt array's `accumulate` has `flow_to[k]` pointing at
  `chain[k+1]` along the whole chain; **a cornered chain drains**: a reach
  whose valley-floor nodes join with a square corner (`(2,4), (3,4), (3,5)`,
  the round 3 case) gives a taut chain (`(3,4)` dropped), the burnt
  array's `flow_to` follows it node by node, and on any chain the burn
  returns no two non-consecutive nodes are equal or 8-neighbours, the
  extension's nodes included (an extension whose least-elevation neighbour
  would make a corner takes the next one, and with none left the end is
  open); a flat floor below the end hits the cap
  (500 m) and gives `end_closed = False`, with exactly 50 extension nodes
  when the floor runs straight (`end_extended_m` 500) and exactly 35 when it
  runs diagonally (`end_extended_m` about 495.0; the stopping rule of step 5),
  and NoData or the window's edge
  before the ground falls does the same; a lower node beside the chain, off
  it, that takes the flow: `drains` is reported false at the first such node
  (the burn is not wrong in that case, the chain is just not the drainage);
  a line three cells off the valley floor gives a chain on the
  floor, and with a corridor of two cells it stays within two; ties; a chain
  that revisits a node (a line that doubles back along its own valley) drops
  the loop and stays 8-connected and taut, with arc lengths recomputed, a
  dropped resample point mapped to the node the cut starts from, and the
  placed node still on the chain; the direction check at 2 m and
  at a reach of 100 m, a flat lake reach passing it; the input array is not
  written; the placed node is the chain node the resample point nearest
  `at` chose (step 6) and its offset is reported; **the valley floor is taken
  across the line** (step 2): on a straight floor falling along the line the
  chosen node is the floor node abeam each point, not one downstream; a node
  exactly `corridor` away is inside; a cross-section with no node takes the
  nearest node; a reach end outside the window is dropped and
  reported, a placed position outside it is refused; **a NoData node on
  the straight join** between two chosen floor nodes refuses the station,
  with either sentinel or with a NaN cell and no sentinel, and with data
  there the reach burns (PR 2's code review, rounds 1 and 2); its refusals
  are `BurnRefusal`, a `ValueError` subclass.
- `test_sensitivity.py`: a smooth gain; a confluence step of known size at a
  known position; a flat floor; a swing of exactly 0.05 is well posed and
  just above it is not; the swing is one-sided (areas 0.04 below and 0.04
  above the placed one give 0.04, well posed); flagged counts trim the
  window on both sides by their bits, not by position, and `checked_up_m` and
  `checked_down_m` say so; **a downstream side not read to `U` is
  `uncertain` (cause `downstream_unread`)**, whether a flag stopped it or the
  mapped river ends before `U` (`reach_down_m < U`), while a river that
  starts within `U` upstream is not penalised; with `U` between nodes (35 m
  on 10 m steps), trusted samples to 30 m and a chain node at or past 35 m
  read to `U`, a chain whose last node is at 30 m does not; `monotone` is strict (two equal
  neighbouring counts give `False`) and a fall or an open end gives it
  `False`; `drains` is false when `flow_to` of a read sample points off the
  chain; each cause is counted apart.
- `test_catchment.py` gains: a synthetic valley with an embankment and a
  gauge beside the river, from a reach: the catchment equals one flood from
  the hand-burnt placed node, and the whole valley is in it; **a refusal
  belongs to the placed node**: NoData reached only by the catchment of a
  node below it (a tributary from the NoData) gives a catchment,
  `downstream_checked` of `partly` or `none` and the cause `downstream_unread`,
  while NoData in the placed
  catchment itself is refused; a NoData node on the chain downstream of the
  placed node is refused with the gap message, and so is a NaN node there
  on a DEM with no sentinel; stage B grows the window past stage A's, and
  the result equals a whole-raster run; a reach with `lakes` is refused;
  `reach=None` gives 22's result bit for bit.
- `test_cli_catchment.py` gains: `--rivers` prints the placement line and
  writes the placement and sensitivity properties; without it, unchanged.

**PR 3, Python** (no network anywhere; the tests of `fetch/` also pin "Data
use" below):

- `test_fetch_nve.py`: every request's `outFields` is exactly the layer's
  allow-list ("Data use") and never `*`, and none names `stasjoneier`,
  `oppdatertav`, `globalid` or a layer 38 discharge-normal field; the
  requests go out one at a time, with the `User-Agent` of `fetch/http.py`
  (tested there on a stub); nothing is requested when the output files exist
  and `--refresh` is absent; the packaged list (140 rows, unique numbers, the
  three spot rows); each station's `name` is layer 0's, not the list's;
  the query URLs (chunks of 40, the station envelopes,
  `outSR=25833`, `f=geojson`); newest version wins, ties to the larger
  `objectid`; a missing station or polygon is refused by name; a reply
  flagged `exceededTransferLimit` is refused by station; a segment seen from
  two stations is kept once by `objectid`; exact geometry copies under two
  `objectid`s are both written (the dropping is the reader's); `objekttype` is
  written as served, null and blank included, with `vatnlnr`; the written files, from canned service answers,
  read back through `io/station_set.py` and `io/rivers.py`; no owner or
  contact field is copied; the manifest's sha256 matches the files.
- `test_station_set.py`, `test_rivers.py`: missing `crs`, duplicate numbers
  or ids, wrong geometry types, a bad station number are refused; a user's
  own points file reads. `test_rivers.py` also pins `kind` on all 25
  `objekttype` values of "The data are not clean" (each lake spelling is
  `lake`; null and blank are `lake` with `vatnlnr` set (not null and not 0)
  and `river` with it null or 0; the
  strays, `FiktivElv` and `BreMidtlinje` are `river`; case and `ø`/`o`
  variants), the raw value kept; exact copies within an `elvid` give the
  segment with the smallest `objectid` and a count of those dropped
  (`read_segments`'s third part and `drop_copies`'s second), a
  MultiLineString segment is refused by `objectid`, while
  equal geometry in two different `elvid`s, and two different geometries
  under one `strekninglnr` (`79.3.0`), are both kept.

**PR 4, Python**:

- `test_reference.py`: agreement on hand-built polygons on a 10 m lattice
  whose edges lie halfway between nodes, so node counts equal areas
  (identical: 100 %, offset 0; a 100-cell square shifted by one cell along x:
  offset **0.5 cell** and overlaps 99 %; grown by one cell on every side:
  404 / 400 = 1.01 cells; disjoint: 0 %). The classes are tested on
  `classify` with the numbers given, not through geometry: exactly 95 %,
  80 %, 30 m and a swing of 0.05; a downstream side not read to `U` gives
  `uncertain` whatever the overlaps are; `refused` beats `uncertain` beats the rest;
  the overlap test is tried first and `match_by` says which passed; the
  summary's percentiles on a known list; band and tile grouping; the
  causes of `uncertain` are counted each (a station with two causes counts in
  both); `mixed_grid` refusals are counted apart as known refusals, with
  their summary line and station numbers, and never among the scored.
- `test_catchment_batch.py`, on the synthetic tiled DEM of
  `test_cli_catchment.py` with a river file and five stations (one matching a
  reference drawn from its own flood, one with a reference shifted to make it
  a miss, one with a confluence just below it, `uncertain`, one whose downstream river ends
  within `U` of it, `uncertain` although its overlap is 100 %, and one with no
  river line near, refused): the rows, the classes, the order, the summary; a
  bug-type exception stops the batch; a sixth station on a DEM of two tiles
  half a cell apart, whose catchment crosses both, is `refused` with
  `refusal_cause = "mixed_grid"` and one tile name from each grid, and the batch goes on.
  `test_mosaic.py` gains: the mixed-lattice refusal is a `MixedGridError`
  with the two names, its message unchanged, and two tiles that differ only
  in spacing raise it too (the cause is general, "The batch").
- `test_cli_station_catchments.py`: the files in `--out-dir`, `--only`, the
  stderr lines; `--reference` absent gives catchments and no scored classes;
  a station or river file without `crs` is refused.

**PR 5, Python**: `test_nearest_stream.py`: with two streams, one nearer and
smaller and one farther and larger, the nearer wins (the rule does not look
at area); the threshold at its boundary; no qualifying node is refused; a
flagged lower bound that already reaches the threshold counts; the row has
no sensitivity and `placed_by = "nearest stream"`; the batch falls back for
a station with no line and not for one with a line.

## Acceptance: every covered HRD station

Run after PR 4 is green (the fallback comparison after PR 5), by `@perf` (it
measures, and owns the evidence layout), under
`docs/benchmarks/<date>/nve-hrd/` with a `run.sh`, the commit, `pmset -g
batt`, the manifest of the fetched station set (its sha256 values, since the
service can change) and the outputs that are not NVE's data:

1. `rasputin fetch-stations nve-hrd` into `../rasputin_data/nve_hrd`.
2. `rasputin station-catchments` over all 140 at the defaults (map radius
   500 m, reach 1000 m, corridor 30 m, outline tolerance); wall time and peak
   memory per station and in total. Expected (estimate, not a measurement):
   the catchments total about 610 M nodes, the final window holds three
   floods, so tens of minutes to a few hours on the M1 Max.
3. **The table**: `results.csv` committed, and the summary by class, size
   band and tile count in the README, with the share `uncertain` per band.
4. **Every finding explained**: each `miss` gets a line (placed on another
   river, `placed_on` not `number`; lake gauge; burn artefact, many lowered
   nodes; NVE polygon disagrees with the DEM; or "not explained"), with a
   re-run at corridor 15 m and 60 m and map radius 250 m and 1000 m;
   each `uncertain` gets its causes from the sensitivity (confluence step and
   where, flat floor or lake, downstream side not read to `U`, chain not
   draining, chain end open, line against the slope);
   `close` rows are summarised by cause. **Expected refusals**: the two
   Finnish-border stations (NoData), the nine shifted-tile stations
   (`refusal_cause = "mixed_grid"`, each row naming one tile from each grid; reported
   apart as known refusals, not failures, Ola's ruling; the README's line
   says a later increment resamples the eight tiles onto the common grid),
   and Femundsenden, which has no river line (until PR 5). A refusal for any
   other reason is a finding, including a window that reaches a shifted tile
   where the catchment does not. **Likely to need that explanation**:
   `156.24.0` and `213.4.0`, whose catchments lie on shifted tiles only, but
   whose 2 km window margin may reach a normal tile; if either is refused
   `mixed_grid`, its line says so and it is counted with the known refusals,
   by cause, not as a failure.
5. **The comparison with nearest stream** (after PR 5): the fallback at 250 m
   on every `miss` and `uncertain` station, the class under each rule side by
   side, and the fallback on Femundsenden. It measures what the mapped river
   buys; it replaces nothing.
6. **Checks that need no NVE data**, on every station that is not refused: the
   placed node's count equals the catchment's node count before reduction
   (the exact oracle, through the real path); the counts along the chain do
   not fall downstream (`monotone`, strict) and each chain node drains into the
   next (`drains`, from `flow_to`); the burn lowered no node off the chain.
   A failure of the first is a defect; a failure of the second or third is a
   station whose burn did not hold, listed (it is `uncertain` by then, and the
   list says how many of the `uncertain` it explains).
7. **What passes the increment**: the batch completes for every station with
   a row each; 22's guarantees hold on every accepted reduced outline (area
   kept, simple, start inside); every finding has its line; the checks of
   step 6 hold. No share of `match` is required this time: this run is the
   baseline the next increments improve (Ola, 2026-10-04).
8. Bygdin is not in the HRD (it is regulated); 22's run stays its record.

## What each persona reads

`@tester` and `@developer`: this file, then `docs/increments/22-auto-catchment.md`
("The seed", "Flow and membership", "The window"), and for PR 3
`docs/increments/23-basin-scale.md` ("The fetch step and the tile cache")
for `fetch/`'s rules. `@developer` also reads `include/terrain/hydrology/upstream.hpp`
and `src_python/tin_engine/fetch/http.py`. `@perf`, for PR 2's figures:
"Placement figures" and `docs/benchmarks/2026-09-29/bygdin-landcover/`
(the layout of a picture's evidence).

## Questions for Ola

Each has a default; the design above is written to the defaults, so the
round can start on them and a different answer changes only the part named.
The earlier questions (which stations, what NVE data to commit, discharge,
what counts as a good catchment, gauges on lakes) are ruled, under "Ola's
rulings"; the way a gauge is placed follows your direction. Questions 1 to 4
below are ruled and closed (2026-10-04): question 2 as its entry says, questions 1,
3 and 4 with their defaults ("Yes to all three."). Questions 5 and 6 came
from PR 2 (2026-10-05) and are open; the design and the code follow their
defaults until Ola rules.

1. **When is a station too uncertain to score?** A gauge's coordinates are
   usually a few tens of metres off the river, and further for some. For each
   station we read the catchment area at points along the river, up and down
   from where the gauge is placed, as far as the gauge's coordinates are off
   the river (at least 30 m). If the area at either end of that stretch
   differs from the area at the placed point by more than 5 % (the same 5 % as
   the "match" bar), or the stretch downstream could not be read that far, or
   the river's flow path does not hold together, the station is reported as
   "uncertain" and is not marked good or bad: a confluence just below it, or
   a flat valley floor, makes the answer depend on where we put the point.
   *Ruled 2026-10-04, closed: the default* (see "Ola's rulings").
   *Default: 5 %, and the distance the coordinates are off the river, at
   least 30 m.* Alternative: a fixed 30 m for every station (fewer
   uncertain, but a gauge 200 m from its river would be treated as exact),
   or a different percentage, or adding the two ends together (stricter). The run shows how many stations each choice
   makes uncertain.
2. **The nine stations on the shifted tiles.** *Ruled 2026-10-04, closed*
   (see "Ola's rulings"): refused in this increment and reported apart as
   known refusals, not failures; a later increment resamples the eight tiles
   onto the common grid. The number is kept so that the others keep theirs.
3. **Burn the whole river network, or only the gauge's own stretch?** Only
   the stretch of river round the gauge (about a kilometre upstream and a
   little downstream) is lowered into the elevation model, so water follows
   the mapped river there. Burning every river would also move divides
   across the map, but it needs care where rivers are mapped less precisely
   than the elevation model. *Ruled 2026-10-04, closed: the default* (see
   "Ola's rulings"). *Default: only the gauge's stretch now; the whole
   network is a later increment, and the residual-inflow work may want it.*
   Alternative: lower every mapped river inside each catchment's window now.
   That also corrects where divides run, but needs the river lines for the
   whole window (far more to fetch than one reach), a pruning step so two
   lines never share a cell, and a way to cope with lines that are wrong; it
   would add an increment before the batch.
4. **The nearest-stream fallback.** For a station with no mapped river within
   500 m (Femundsenden, a lake gauge, is the one of the 140), and as a
   comparison on every station that ends uncertain or a miss, the gauge can
   be moved to the nearest elevation-model stream within 250 m (Jenson's
   rule). It never uses area. It is the last, small piece (about 80 lines)
   and decides nothing for stations the river placed. *Ruled 2026-10-04,
   closed: the default* (see "Ola's rulings"). *Default: include it,
   as the last piece.* Alternative: leave it out; Femundsenden is then
   refused and the comparison is not made.
5. **The river's path in the DEM crosses a gap.** The path we burn joins
   the lowest DEM points beside the mapped river with straight steps, and a
   step can pass through a cell with no elevation (NoData). Left alone, the
   no-data value is burnt in as if it were a height, and the station's
   numbers come out absurd. *Open (2026-10-05). Default: refuse the station,
   with a message saying the river line crosses a gap in the DEM.*
   Alternative: route the path round the gap; that is a new placement rule
   (which side, how far) and would be designed first. How many stations it
   touches is not known yet; the acceptance run counts them.
6. **A river point just upstream of the gauge that does not drain through
   it.** Within the gauge's position uncertainty upstream, a river point
   whose own catchment reaches the edge of the area read, or a gap in the
   DEM, cannot drain through the gauge (the gauge's catchment has neither,
   or the station is refused). So the river there lies outside the gauge's
   catchment, which is the kind of position sensitivity the check exists to
   catch. *Open (2026-10-05). Default: yes, mark the station uncertain*
   (cause "chain not draining"). Alternative: only stop reading the area
   there, with no other effect, as the first code did; such stations are
   then scored.

## Review

(Renumbered 28 to 29 on 2026-10-04: increment 28 is `28-no-fma-contraction.md` on branch `worktree-fpc`. The rounds below were recorded as "28" and are left so.)

**28, design review, round 1, 2026-10-04.** Range `104c883..682c36c` (design and ROADMAP row only). Verdict: CHANGES REQUESTED. LOC: 0 production; estimates PR 1 about 190, PR 2 about 480, both under the ceiling, but 22's `catchment.py` ran 44 % over its estimate and PR 2 names no split seam. Re-measured and true: the HRD PDF (sha256, quotations, 140 unique rows); layers 0/14/38 (fields, EPSG:25833, all 140 stations, 130×1 + 10×3 polygons, 42 points outside, farthest 233 m); areas, size bands and tiles per catchment; no NoData within 250 m of any station; NLOD; HydAPI 401; all seven DOIs; the accumulation oracle (the flood's visit order does not depend on the seed). Blocking: (1) the "500 unregulated" are 364 with regulation 0 plus 136 with none recorded (`:68-73`, Question 1); (2) Jenson 1991 is the nearest-stream-cell snap, and Lindsay et al. 2008 prefer it over the max-accumulation rule chosen here (`:193-214`, Question 4); the departure is unnamed and Ola ruled (a) without that option; (3) nine stations (`156.15.0`, `196.11.0`, `206.3.0`, `208.2.0`, `208.3.0`, `209.4.0`, `212.49.0`, `213.2.0`, `223.2.0`) straddle the eight half-cell-shifted tiles and will be refused by the mosaic's mixed-grid rule, contrary to `:133-135`, acceptance step 4 and the ROADMAP row; (4) the divide offset of a square shifted one cell along an axis is half a cell, not one (`:640`); (5) median vertices 969.5 (not 983), inside-distance quartiles 20/75/355 m (not 20/78/351), median area 131 km² (not 135). Not pushed; no CI.

**28, design review, round 2, 2026-10-04.** Range `26588b1..fe01d37` (design and ROADMAP row only). Verdict: CHANGES REQUESTED. LOC: 0 production; estimates PR 1 about 150, PR 2 about 375 (540 with margin), PR 3 about 300, PR 4 about 295, PR 5 about 80, all under the ceiling. Round 1 blockers 1, 2, 4 and 5 are closed; blocker 3 is closed in this file but not in the ROADMAP row. Re-measured and true: ELVIS layer 2 (fields, EPSG:25833, 1,954,539 features), nearest-line lengths 392/751/1196 m (maximum 4673), distances (median 20 m, maximum 461), 113 own watercourse number, 46 lake lines, 16 and 5 second branches from the station point; the three new DOIs and the Lindsay and WhiteboxTools quotations; the accumulate flags equal `upstream`'s `touches_edge` and `touches_nodata`. Blocking: (1) the strict descent does not make each chain node drain into the next (a pit at a lowered chain end is filled at its spill level and ordered breadth-first; a lower neighbour off the chain can take the flow), so `:626-630`, `:672-673`, `:718-720` and "nesting by construction" `:770-773` overclaim: extend the chain past a lowered end, trust counts by their flags, make `monotone` strict; (2) the loop rule `:615-616` breaks 8-connectivity: cut the loop; (3) a station whose downstream side was not read to U can be scored: make it `uncertain`; (4) ELVIS has 25 `objekttype` spellings including null (station `152.4.0`) and serves exact duplicate segments under different `objectid`s (2 copies at `82.4.0`, 5 at `139.35.0`): a total lake/river mapping, de-duplication by geometry, a rule for forks; (5) `confluence_near` within 100 m of P is 19 stations, not 16 (`:584-585`); (6) PR 3 needs PR 2's `io/rivers.py` (`:1035-1036`); (7) `ROADMAP.md:54` omits the nine shifted-tile refusals and Femundsenden. Suggested: a short "Data use" section with explicit field allow-lists (never `stasjoneier` or `oppdatertav`); a one-sided swing measure, or say why it adds both sides; Soille et al. 2003 carving as prior art. Not pushed; no CI.

**29, design review, round 3, 2026-10-04.** Range `fe01d37..a14e7e6` (design and ROADMAP row only, including the renumbering). Verdict: CHANGES REQUESTED. LOC: 0 production. Estimates add up: PR 1 155; PR 3 380 (547 with the 44 % margin); PR 2 415 (598); PR 4 295 (425); PR 5 80 (115). All are under the 700-line ceiling. Round 2 blockers 1 to 7 and its three suggestions are closed as written. No "28" referring to this increment is left (`git grep` over the tree). Re-measured and true: ELVIS `objekttype` has 25 values summing to 1,954,539, including null (94) and a blank (2, a single space), with eight lake spellings, all starting `innsj` when casefolded. Exact copies: 9 pairs at `82.4.0` and 4 groups of five at `139.35.0`. At `79.3.0`, one `strekninglnr` has two geometries under one `elvid`. `152.4.0` has five null-type lines. On the nine straddling stations, neither lattice's tiles cover the NVE polygon. `156.24.0` and `213.4.0` are covered by the shifted lattice alone. The shift is in x only (world files). Soille et al. 2003 checked against Crossref and the abstract quotation. The float32 spacings. `WINDOW_MARGIN_M` = 2000. The three at-risk citations from `check_citations.py` were re-read and hold. Blocking: (1) A chain with a square corner, or any two non-consecutive chain nodes that are 8-neighbours, breaks `drains` and `monotone` even after a correct burn. The flood names a node's flooder when it is first pushed. So the lower node `chain[k+2]` pushes `chain[k]` before `chain[k+1]` can. Simulated with the flood of `upstream.hpp`: the path (2,4),(3,4),(3,5) gives `flow_to` (2,4)→(3,5) and (3,4)→(4,5); the same path with a diagonal corner drains node by node. Straight lines between consecutive valley-floor nodes, loop cuts and the end extension all make such corners. Make the chain taut before the descent and after the extension: where `chain[k]` is an 8-neighbour of `chain[j]` with `j > k+1`, drop `chain[k+1..j-1]`, as the loop cut does. Add a test that a cornered chain drains. (2) `:1474-1476` sends PR 2 to 23's fetch rules, but `fetch/` is PR 3 since round 2. (3) `:749-750` and `:786`: "at most 0.14 m" no longer holds once the end extension exists. The chain can reach 1000 + 600 + 500 m, which is about 210 nodes, so about 0.21 m. (4) `:239` "254 tiles of 5051 × 5051" is false for 11 tiles. Seven shifted tiles are 5052 × 5053 and `7507_4` is 5052 × 4103. `7305_3` is 5051 × 2881, `7405_1` 3521 × 5051 and `7405_2` 3511 × 5051. This claim predates the range. While fixing it, state the cause of the shift, measured by the main session on 2026-10-04: the data were resampled onto a shifted grid, not mislabelled, so the later fix is a resample. Not pushed; no CI.

**29, design review, round 4, 2026-10-04.** Range `a14e7e6..e7b5923` (design and ROADMAP row only). Verdict: CHANGES REQUESTED. LOC: 0 production. Estimates: PR 1 155; PR 3 380 (545); PR 2 415 (598); PR 4 310 (446); PR 5 80 (115), all under the 700-line ceiling. The 130-line figure script is not counted, by Ola's ruling of 2026-10-04 ("map questions: 1: No, 2: Yes."). Round 3's blockers 1, 2 and 4 and its four suggestions are closed; blocker 3 is closed except for the node count in blocker 1 below. The taut rule was checked with an independent copy of `upstream.hpp`'s flood on 2,000 random chains with repeats, loops and turn-backs: 1,738 raw chains failed `drains`, no taut one did, and every output was taut and 8-connected. On 400 runs of the end extension across an embankment, the step rule kept every chain taut, and no node was flooded from a non-consecutive chain node; the remaining failures all start at a lower node beside the chain, which `drains` reports. Re-measured and true: the eleven tile sizes; the x shift on exactly the eight tiles, with y on the lattice; PixelIsArea and tag/`.tfw` agreement on all 254; the overlap differences of `7707_3` and `7507_4`; the 0.21 m bound; 16c's LOC table without `render.py`; matplotlib in no dependency list; `mosaic.py:385` and `catchment._plan` as described; 15a's tests catch by `MosaicError`, so a subclass leaves them unchanged. Blocking: (1) `:832-838`: "36 nodes if every step is diagonal" and "ends at most 500 m from the reach's last node" disagree, since 36 diagonal steps are 509 m. State the stopping rule (stop before a step that would pass 500 m: 35 diagonal, 50 straight), as the red suite's cap test depends on it. (2) `:864-866`: "the catchment is unchanged" was edited in this range and is contradicted by the new `:853-858` (lowering a depression's spill point moves its drainage). Qualify it: the catchment is unchanged except where a chain node was another depression's spill point. Not pushed; no CI.

**29, design review, round 5, 2026-10-04.** Range `e7b5923..786584e` (design file only). Verdict: APPROVED. LOC: 0 production. Estimates: PR 1 155; PR 3 380 (547); PR 2 415 (598); PR 4 310 (446); PR 5 80 (115), all under the 700-line ceiling. The 130-line figure script is not counted (Ola's ruling of 2026-10-04). Round 4's two blockers and five suggestions are closed. 786584e is complete: status line, round 4 entry, no half-edited section. Re-measured and true: 35 diagonal steps 494.97 m, 36 steps 509.12 m, 50 straight steps exactly 500.0 m; on all 254 DTM10 tiles, spacing 10 m, NoData −32767, EPSG:25833 and PixelIsArea (`tifffile`); `_mixed`'s five causes and the two tiles `_covering` passes it (`mosaic.py:385`, `:571-593`); the neighbour sets of `7707_1` (4) and `7707_3` (7) under the stated rule; `7707_3`'s eleven values, to 0.01 m, including both near ties; `7707_1` declared 0.02 to 0.11 m, moved 0.19 to 0.45 m (the file says 0.44). Required with this record: `ROADMAP.md:54` still says "awaiting round 4". Suggested: 0.45 m without the parenthetical; "10:51 UTC"; the status sentence names `7707_3` and the survey change; the mixed-grid summary line adds "and neither grid covers the window alone". Not pushed; no CI.

**29 PR 1, code review, round 1, 2026-10-04.** Range `51275f5..edef966` (red 0624bbf and 8d87f77, green e43020b, scaffolding removal edef966). Verdict: CHANGES REQUESTED. LOC: 182 added and 53 removed, 129 net. Estimate "about 155" added: 17 % over, inside the 44 % margin, far under 700. Own Release build ctest 932/932; ASan/UBSan pass; full pytest 4069 passed, 17 skipped; gates green. Mutation round: callback 8/8 killed; reverse pass 7/8 (survivor equivalent: order[0] is always the first outlet); starting bits and size check 7/8 (survivor: limit 0xFFFFFFFE, needs a 40 GB raster); flood.hpp's new code 4/4. Blocking: (1) project_structure.md:384 says no accumulation pass and :75-78 lacks flood.hpp and accumulate.hpp; fix in PR 1. (2) Status line (29-nve-reference-catchments.md:3) and ROADMAP.md:54 say "ready for the red step of PR 1". Before push: merge master (120 behind, ROADMAP conflict), re-run check_citations. Suggestions: bindings/core.cpp:1268 say "the node's flooder"; tests/python/test_core_accumulate.py:22 "not on CI runners"; tests/cpp/unit/test_hydrology_accumulate.cpp:14-28 drop restated interface; design :682-685 add on_reach(j, j) and NoData-adjacent outlets. Not pushed; no CI.

**29 PR 1, code review, round 2, 2026-10-04.** Range `edef966..9736bfd` (81392f5 round 1 recorded; 0b86517 project_structure, status line and ROADMAP; 725c8cd design states on_reach(j, j) and the NoData-adjacent outlets; c59ed6b merge of origin/master 99093af; 5dbfec0 keeps project_structure.md:166 in place; f5e8e76 accumulate docstring; 9736bfd two test comments). Whole PR `99093af..9736bfd`. Verdict: CHANGES REQUESTED, on citations alone. LOC: 181 added and 53 removed, 128 net (round 1: 182/129; the only production edit since then, f5e8e76, sits inside a raw docstring, and the one-line difference comes from how the count treats the diff's alignment at the upstream docstring's closing line). Estimate "about 155" added: 17 % over, inside the 44 % margin, far under 700. Fresh Release build ctest 994/994 (test_hydrology_accumulate: 11 cases, 163,676 assertions); extension rebuilt into the worktree venv, full pytest 4420 passed, 17 skipped; mypy, ruff check, ruff format, prohibited-deps, detria boundary and check_citations all exit 0. All round-1 blockers and suggestions closed, and each new claim was checked against flood.hpp, accumulate.hpp and the binding. No red-step scaffolding. The merge left hydrology/ and raster/ untouched, so round 1's mutation record stands. Blocking: the branch adds 8 lines to bindings/core.cpp above line 927 (one include at line 16, the BoundAccumulate struct at about line 316), so three live citations now point 8 lines too high: `15f-edge-strip.md:565` and `:1537` cite `bindings/core.cpp:1174` (the gil_scoped_release is now at 1182), and `27-node-sampling.md:133` cites `bindings/core.cpp:927-932` (the sample docstring is now at 935-940). Re-cite them and re-run check_citations. Not pushed; no CI.

**29 PR 1, code review, round 3, 2026-10-04.** Range `9736bfd..36c2370` (one docs-only commit: the round-2 record, `15f-edge-strip.md:565` and `:1537`, `27-node-sampling.md:133`, the status line and ROADMAP row 29). Whole PR `99093af..36c2370`. Verdict: APPROVED. LOC: 181 added and 53 removed, 128 net, unchanged from round 2 (no production file in this range). Estimate "about 155" added: 17 % over, inside the 44 % margin, far under 700. `git diff --stat 9736bfd..36c2370` lists only ROADMAP.md, 15f-edge-strip.md, 27-node-sampling.md and 29-nve-reference-catchments.md (7 added, 5 removed lines), and nothing under bindings/, include/, src_python/ or tests/. So round 2's build, test, gate and mutation results stand. Round 2's three blockers are closed. Both 15f citations now give `bindings/core.cpp:1182`, which is the `const py::gil_scoped_release unlocked;` in the `refine_points` binding (lines 1173-1183). `27-node-sampling.md:133` now gives `bindings/core.cpp:935-940`, which is the `sample` docstring, "Bilinear z at each of the (N, 2) points" through "never NaN.". `check_citations.py --base origin/master` exits 0. All 34 lines on its at-risk list were re-read as quotations. The live ones hold: `05b-noder-driver.md:1749` (`tests/cpp/CMakeLists.txt:189`), `test_features.py:583` (`project_structure.md:166`), `29-nve-reference-catchments.md:55` (`ROADMAP.md:54`), and the three just re-cited. The rest sit in dated review records or a retrospective, which are history and left as written. The status line and ROADMAP row 29 match the tree. The branch contains origin/master `99093af`. No red-step scaffolding. No @perf run is needed: PR 1 touches neither refine nor mesh code (design, line 1395). Suggestion: `15f-edge-strip.md:1537` says the number was recomputed when "increment 29 PR 1 took in master". In fact the 8-line shift comes from PR 1's own additions to `bindings/core.cpp` (the include and `BoundAccumulate`), not from a master merge. Not pushed; no CI.

**29 PR 3, code review, round 1, 2026-10-04.** Range `a1445c6..aae91bb` (rulings `b519bd9`, red amendment `05f348a` and `c18c8d5`, green `aae91bb`; the red step `213a2ef` before it). Verdict: CHANGES REQUESTED. LOC: 348 added and 7 removed, 341 net, under the estimate of about 380 (355 once the GeoJSON writer moves to PR 4) and far under 700. Blocking: (1) NVE's river service sends `vatnlnr` = 0 for "no lake" (live query: 956,447 features with 0; both blank-type features, `objectid` 11166506 "Cap'pirjåkka" and 11513676, have 0 and are rivers), but `io/rivers.py`'s `kind_of` treats any non-null value as set, so a blank-type river with 0 is a lake; the design's "set" (`:772`, `:1678` at `aae91bb`) must say not null and not 0, and "set on 496 river features" (`:300-301`) counted the zeros (the reviewer's sample has 2 river features with a positive number); needs a red test and a fix ("PR 3, after code review round 1", under "Ola's rulings"). (2) Status line and `ROADMAP.md:54` still say a red amendment is to come. (3) `15f-edge-strip.md:591` cites `cli.py:1561` and `:1580`, which PR 3 moved to `:1594` and `:1613`. (4) The PR 3 file table lists the GeoJSON writer's move to `io/geojson.py`, which PR 3 neither does nor needs. (5) `project_structure.md` lacks `io/station_set.py`, `io/rivers.py`, `repository.py`'s `read_json`, the `data/` folder, and the `fetch/` and `sources.py` entries missing since 23a-2. Suggestions: a malformed user file refused with `ValueError`, not `TypeError` or `KeyError`; `fetch-stations` reports a missing field in plain words, not a raw `KeyError`; `@tester`'s present-tense "HOW THIS FILE GOES RED" paragraphs (left to `@orchestrator`, not this round). Not pushed; no CI.

**29 PR 3, code review, round 2, 2026-10-04.** Range `aae91bb..a4e4726`; whole PR `d926644..a4e4726`. Verdict: CHANGES REQUESTED. LOC for the whole PR: 387 added, 7 removed, 380 net (tokenize count; green `c4e50ff` alone 36 added, 15 removed, 21 net), against the estimate of about 355 (511 with the 44 % margin): inside the margin, far under 700. Full pytest at HEAD on a rebuilt `_core`: 4499 passed, 118 skipped, 8 failed, all eight `test_settings_wiring`'s executable-hook check, which fails only in the git-archive scratch copy (27/27 in the real worktree). mypy, ruff check, ruff format, prohibited deps, detria boundary and check_citations all clean. Red `b59c648`: 21 tests fail for the intended reasons; all pass at HEAD; green touches no test file. All round-1 blockers and suggestions closed. No mutation testing owed: only `accumulate`'s oracle is invariant-critical, and its kill record is with PR 1. No `@perf` run: PR 3 touches no refine or mesh code. Blocking: (A, `@tester`) `tests/python/test_features.py:583` cites `project_structure.md:166`, now `:190`; (B, `@architect`) the status line and `ROADMAP.md:54` still say the red test and the fix "come next"; (C, `@tester`) the present-tense "HOW THIS FILE GOES RED" paragraphs in `test_fetch_nve.py`, `test_rivers.py` and `test_station_set.py` describe the tree before green and are false now. The PR 3 round-1 record above stays as written. Suggestions, for PR 4 or later: a feature that is a JSON list gives `AttributeError`, and a malformed reference polygon gives shapely's `TypeError`, not a `ValueError` naming the fault; a one-vertex LineString is accepted as a river segment (PR 2's chain may want it refused); `fetch-stations` catches `KeyError` over more than the service answer (narrow the `try`); stale `cli.py` citations in 15e and 25 predate this branch (for `@orchestrator`'s re-citation sweep). Not pushed; no CI.

Fixes for round 2: A and C by `@tester` in `f91cdaa` (re-cited to `project_structure.md:190`; the three paragraphs deleted); B by `@architect` in the commit that records this round (the status line and `ROADMAP.md:54`). The PR 1 records above (round 2, round 3) cite `project_structure.md:166` as the file was then and are left as written.

**29 PR 3, code review, round 3, 2026-10-04.** Range `a4e4726..c7d427a` (`f91cdaa`, `c7d427a`); whole PR `d926644..c7d427a`. Verdict: CHANGES REQUESTED. Delta: 7 lines added, 16 removed, in `ROADMAP.md`, this file and three test files; no production code, so the count stays 380 net and round 2's build, test and gate results stand. Round 2's blockers A, B and C are closed. `check_citations.py --base d926644` exits 0; the at-risk lines this delta could move re-read as quotations and hold. `ruff check` and `ruff format --check` clean. Blocking: the status line and `ROADMAP.md:54` said round 2 asked for "two citations"; it asked for one (`test_features.py:583`). Fixed by the main session in the commit that records this round.

**29 PR 3, code review, round 4, 2026-10-04.** Range `c7d427a..5bdcc93` (`c1eeec4`, `5bdcc93`; main session, docs only); whole PR `d926644..5bdcc93`. Verdict: APPROVED. Delta: 5 lines added, 2 removed, in `ROADMAP.md` and this file; production count stays 380 net (estimate about 355, 511 with the margin). Round 3's blocker closed: the status line and `ROADMAP.md:54` say one citation. The round-3 record's range, commits and counts match `git log` and `git diff --shortstat a4e4726..c7d427a`. `check_citations.py --base d926644` exits 0; no at-risk line moved by this delta. No mutation record or `@perf` run owed. CI is owed after the push: `gh pr checks` green before merge.

**29 PR 3, code review, round 5, 2026-10-04.** Range `3dff7bc..ea489e5` (master merge `af2bc38`; `4e877ae` and its revert `09a44c6`; Ola's ruling `9607b98`; `ea489e5` deletes the `NOTICE.md` test). Verdict: CHANGES REQUESTED. Production lines unchanged since round 4. The merge is clean (`git merge-tree --write-tree 3dff7bc 791abb6` gives the tree of `af2bc38`). Full suite on a rebuilt `_core`: 4572 passed, 118 skipped, exit 0, so the prose-read hook stayed silent. `ruff check` and `ruff format --check` clean. The remaining NVE credit tests exist (CSV header, catalogue, every fetch's `NOTICE.txt`); `NOTICE.md` still credits NVE. `check_citations.py` exits 0. Blocking: the ruling said only the 3.14 leg failed, but all three Python legs did; "Data use" said every rule is held by a test, including the `NOTICE.md` credit. Suggestions: name whose §"Licence" and §"Data use" are meant; update the status line. All four fixed by the main session in the commit that records this round.

**29 PR 3, code review, round 6, 2026-10-04.** Range `ea489e5..741504f` (one docs commit, main session; `ROADMAP.md` +1/-1, this file +9/-5). Verdict: APPROVED. Round 5's two blockers and two suggestions closed: the ruling names all three Python legs (`gh pr checks 178` agrees); "Data use" excepts the `NOTICE.md` credit, checked by hand, and the tested credits are tested (`test_fetch_nve.py`: CSV header, catalogue, `NOTICE.txt`); `NOTICE.md`'s "Data in the package" exists; status line and ROADMAP row 29 match the history. `check_citations.py` exits 0; the lines this commit moves are cited only by dated records, except this file's `:55` on `ROADMAP.md:54`, which holds. Master's `accd52a` (h13) touches no file on this branch; `git merge-tree --write-tree HEAD origin/master` is clean, so no merge before the push. CI is owed on the new head; red CI there voids this approval.

**29 PR 2, code review, round 1, 2026-10-05.** Range `529613a..0588bfa` (red `3a06ffc`, amendment `5357785`, green `91857df`, `@architect`'s reading `6f55af3`, second amendment `3c845da`, green `0588bfa`). Verdict: CHANGES REQUESTED. LOC by tokenize over the five `src_python` files: 578 added, 20 removed, 558 net (`burn.py` 167, `gauge.py` 117, `sensitivity.py` 85, `catchment.py` 94, `cli.py` 95 net). This record uses the net figure, as PR 3's did; against the estimate of about 415 (598 with the 44 % margin) it is inside the margin on either basis, and far under 700. Full suite on a rebuilt `_core`: exit 0, 4799 passed, 17 skipped, the prose-read hook silent; mypy, ruff check, ruff format, prohibited deps, detria boundary and check_citations green. Red before green for both pairs: `5357785` 21 failed and 103 errors, then `91857df` 171 passed; `3c845da` on `91857df` 3 failed, then `0588bfa` 172 passed. The suite can fail on both new rules: the whole-disc floor rule fails 8 tests, and dropping the straight-before-diagonal tie fails the embankment test. No mutation testing owed (only `accumulate`'s oracle is invariant-critical, killed with PR 1); no `bench.py` run (no refine or mesh code). `@perf`'s placement figures are still owed before PR 2 is complete. Blocking: (1) a NoData node on the chain: `_floor_node` never picks NoData, but `_join` (`src_python/tin_engine/burn.py@0588bfa:86-97`) links chosen nodes by straight 8-connected runs without looking at the nodes between, and the taut cut can map the placed node onto one. Probe: float32, NoData −32767, floor in column 20 above row 20 and column 22 below, `z[20, 21]` NoData, line down column 20, corridor 30 m: the chain runs through `(20, 21)` and everything below is burnt to about −32767.002 m; with a positive sentinel (3.4e38) the hole is burnt to 290.499 m and `lowered_max_m` of about 3.4e38 reaches the catchment file; the direction check's means read the sentinel too. Write the rule into step 2, covering the taut cut and the direction check, with the refusal as the default for Ola. (2) "Sensitivity", step 2: "its own count is not a sample and its flags do not matter" holds only where `D` lies past `U`; a node exactly at `U` is a sample in `(0, U]` and must be trusted (a flag at +30 m with `U` = 30 gives `downstream_unread`, which the code does correctly); fix the design, and the comment at `src_python/tin_engine/sensitivity.py@0588bfa:71-72` with the green step. (3) The status line and `ROADMAP.md:54` behind the tree; the estimate block's 556; `25-plain-output.md:53`'s citation of the catchment properties. (4) The upstream-flag rule presented as ruled by `@architect`: mark it provisional, an open question for Ola. Suggestions: state the cross-section's bound (half a step; a downstream lean of 3 to 4.5 m on average on a floor flat across, at 30°) rather than "at most 4.6 m either way"; "a cross-section with no node" also happens when every node in it is NoData; `25-plain-output.md`'s field table (lines 89-98, 110-111, 128) cites `cli.py` lines that no longer quote it: pin them to the revision described. Not pushed; no CI.

Fixes for round 1, by `@architect` in the commit that records it: (1) "No NoData on the chain" in step 2 (refuse after the taut cut, before the direction check; question 5, default refuse), the red and green asks in the block "PR 2's code review, round 1", the refusal in the class table and the red suites; (2) the wording of `D` in "Sensitivity", step 2; (3) the status line, `ROADMAP.md:54`, the estimate block (558 net, 578 added; `sensitivity.py` 85). `25-plain-output.md:53` cites `cli.py:1925-1928`, and at `0588bfa` `"fine_area_m2"` is at line 1925 and `"outline_tolerance_m"` at 1928 (`grep -n`), so the citation holds and is left as it is; the round's 1929 did not reproduce. (4) The upstream-flag rule marked provisional, question 6, default yes. Suggestions taken: the cross-section's bound, the all-NoData cross-section, and the field table's `cli.py` citations pinned to `586fbc1`, the revision the inventory describes. The production fixes (the NoData refusal and the `sensitivity.py` comment) come with the red amendment and green step next.

**29 PR 2, code review, round 2, 2026-10-05.** Range `0588bfa..96b3881` (`84b5ede` round 1 recorded, red `6c336f4`, green `96b3881`); the whole PR is `529613a..96b3881`. Verdict: CHANGES REQUESTED. LOC over the five `src_python` files of `529613a..96b3881`: 588 added, 20 removed, 568 net (`burn.py` 174, `gauge.py` 117, `sensitivity.py` 85, `catchment.py` 97, `cli.py` 95 net), counted by round 1's method, written here so it can be rerun: `tokenize` decides which lines are code (no blank lines, no comment-only lines, no docstrings; a multi-line string counts its first line only); the added lines are the `+` ranges of `git diff -U0 529613a 96b3881` hunks, judged as code in the file at `96b3881`, and the removed lines are the `-` ranges, judged in the file at `529613a`. Inside the 598 margin and under 700. Red before green: `6c336f4`, 3 tests fail with DID NOT RAISE; `96b3881`, 132 pass. Full suite on a freshly built `_core`: 4795 passed (7 `test_settings_wiring` failures come from the `git archive` copy the run used; 27 of 27 pass in the worktree); the prose-read hook silent; mypy, ruff check, ruff format, prohibited deps, detria boundary and check_citations green. The gap check covers every reader for finite sentinels. No mutation testing owed, no `bench.py` run (no refine or mesh code); `@perf`'s placement figures still owed. Blocking: (1) NaN cells: `burn_reach`'s mask is `raw != m.nodata`, all true when `nodata` is None (`src_python/tin_engine/burn.py@96b3881:134`); the core counts NaN as NoData whatever the sentinel and `RasterMeta` never holds a NaN sentinel, so a NaN-gapped float DEM is burnt with no refusal, `lowered_max_m` NaN, the sensitivity well posed and NaN in the catchment file. (2) The status line, `ROADMAP.md:54` and the estimate block behind the tree. (3) "Stage B's windows give the same chain as stage A's" overstates: stage B's extension can run further; the same taut chain is what holds. Suggestion: `burn_reach` raises its own `ValueError` subclass and only that becomes a `CatchmentError`, so stage B's `except CatchmentError: pass` (`src_python/tin_engine/catchment.py@96b3881:290`) cannot swallow an ordinary bug. Not pushed; no CI.

Fixes for round 2, by `@architect` in the commit that records it: (1) "No NoData on the chain" in step 2 says a NaN cell is NoData, as in the core, with one mask for the cross-section, the chain check and the extension; the red and green asks in the block "PR 2's code review, round 2"; the red suites name the NaN cases; (2) the status line, `ROADMAP.md:54` and the estimate block (568 net, 588 added; `burn.py` 174, `catchment.py` 97); (3) "the same taut chain", with the reason. The suggestion is taken: `BurnRefusal` in step 2 and in the green ask, pinned by `@tester` in the red step. The LOC count above was rerun by `@architect` with a script written from the rule (not committed) and gave the same figures, file by file.

**29 PR 2, code review, round 3, 2026-10-05.** Range `96b3881..193079d` (`cc52b8f` round 2 recorded, red `9625478`, green `193079d`); the whole PR is `529613a..193079d`. Verdict: APPROVED. LOC by round 2's rule over the five `src_python` files: 589 added, 20 removed, 569 net (`burn.py` 175, `gauge.py` 117, `sensitivity.py` 85, `catchment.py` 97, `cli.py` 95 net), one net line more than round 2; inside the 598 margin and under 700. Round 2's blockers are closed: (1) the mask at `src_python/tin_engine/burn.py@193079d:140` excludes NaN whatever the sentinel, matching the core's `is_nodata`; (2) the status line, `ROADMAP.md:54` and the estimate block matched the tree as round 2 left it; (3) the design says "the same taut chain". The suggestion is taken: `BurnRefusal` is raised by all three of `burn_reach`'s refusals, and `_burnt_flood` catches only it. Red `9625478` fails 8 tests for the intended reasons; all pass at `193079d`. Full suite on a freshly built `_core`: exit 0, 4783 passed, 18 skipped (`test_settings_wiring` 27 of 27 run in the worktree); the prose-read hook silent; mypy, ruff check, ruff format, prohibited deps, detria boundary and check_citations green. No other code in the PR builds its own DEM mask: `gauge.py` reads no array, `sensitivity.py` reads `accumulate`'s output, `catchment.py` uses only the array's `.shape`. A probe with NaN cells and a finite sentinel is refused at `193079d`. No mutation testing owed, no `bench.py` run (no refine or mesh code); `@perf`'s placement figures are still owed before PR 2 is complete. The status line lagging the review does not block; it is set in the commit that records this round. Suggestions, not blocking: a fourth `GAPS` case, `(math.nan, -32767.0)`, in `test_burn.py`, to pin NaN beside a sentinel; `burn_reach`'s docstring (`src_python/tin_engine/burn.py@193079d:135-137`) names two of its three refusals (not "no DEM node with data near the reach"). Not pushed; no CI.

Recorded by `@architect`, with the status line, `ROADMAP.md:54` and the estimate block (569 net, 589 added; `burn.py` 175). The LOC count was rerun from round 2's rule and gave the same figures, file by file. Both suggestions are left for a later step, to avoid another review loop for a test case and a docstring: whoever next touches `burn.py` or `test_burn.py` (PR 4 at the latest) adds the `(math.nan, -32767.0)` case and names the third refusal in the docstring.
