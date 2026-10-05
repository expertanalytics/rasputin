# Increment 29 — NVE reference catchments: our catchments against NVE's, station by station

Status: **PR 1, 2, 3 and 4 merged (#173, #179, #178, #182); the acceptance run over all 140 stations is in, at `docs/benchmarks/2026-10-05/nve-hrd/` (as fixed in `cd8b5be`): 74 match, 5 close, 6 miss, 39 uncertain, 16 refused (14 of them expected), 33 minutes and at most 7.7 GB; evidence review round 2 approved, the push waits for Ola; PR 5 (the nearest-stream fallback) and acceptance step 5 (the comparison with nearest stream) not started; questions 6 and 10 go back to Ola with the run's figures**
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
ruled (see "Ola's rulings"). Questions 5 and 6 came from PR 2's code
review, 7 from its placement figures, and 8 from PR 4's green step. On
2026-10-05 Ola ruled that lake gauges get increment 22's Bygdin method, in
PR 4 ("Lake gauges"); questions 9 and 10 come from that design. The same day
Ola ruled questions 5 to 10 as built (their defaults), to be reassessed after
the full 140-station run ("Questions 5 to 10: Ola's ruling").

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
    --lakes ../rasputin_data/nve_hrd/lakes.geojson \
    --out-dir ../rasputin_scratch/hrd_catchments
rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 \
    --domain ../rasputin_scratch/hrd_catchments/2.32.0.geojson --tolerance 1 --out atnasjo.vtk
```

The example station is `2.32.0` Atnasjø (question 7, default), not `2.11.0`
Narsjø: PR 2's placement figures
(`docs/benchmarks/2026-10-05/nve-placement/README.md`) give Atnasjø a well
defined 459.8 km² catchment that agrees with NVE's to 98.4 % / 99.2 %, and
Narsjø an uncertain 0.0012 km² one.

It generalises increment 22's Bygdin acceptance (one catchment, area within
2 %, node overlap both ways) to 140 stations. It also adds what 22 left out:
flow accumulation; placing a gauge where it physically is (on NVE's mapped
river, then on the DEM's flow path along it); and a per-station check of
whether the catchment area is well defined at the gauge at all. A gauge on a
lake is seeded with the whole lake, as 22 seeded Bygdin ("Lake gauges").

**Not closed.** Discharge (only the keys to join it later are stored; see
"Discharge"). Residual inflow between gauges on one river (the interfaces are
sketched so that this increment does not block it; see "Residual inflow,
later"). Burning the whole mapped network into the DEM (only the gauge's own
reach is burnt). The nine stations whose catchments straddle DTM10's
half-cell-shifted tiles (expected refusals, reported apart as known
refusals, not failures; a later increment fixes the tiles) and
Femundsenden (no mapped river within 500 m; refused until the fallback).
Seeding a lake together with the river from its outlet down to a gauge
below it (question 10; such a gauge stays on the river, "Lake gauges").
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

**PR 4's red step: `@tester`'s choices, read by `@architect` (2026-10-05).**
Not Ola's rulings. Red `2b9b39f` (`test_reference.py`,
`test_catchment_batch.py`, `test_cli_station_catchments.py`,
`test_io_geojson.py` and `batch_fixtures.py`, new; additions to
`test_mosaic.py`, `test_catchment.py` and `test_burn.py`). **PR 4 is built on
PR 2's branch, which is not pushed yet**, so PR 4's base moves when PR 2 is
merged; and questions 5, 6 and 7 for Ola are still open. Their defaults are
in the code PR 4 builds on; a different answer changes PR 2's code, and PR
4's suite only where a station meets the rule. None of the batch fixture's
stations does: none has a gap on its chain (question 5), the causes the
suite expects were measured on PR 2's code, and with question 6's
alternative causes can only disappear, which leaves the expected causes as
they are (`1.140.0` and `1.135.0` have none, `1.130.0` has only
`downstream_unread`, `1.150.0` keeps `swing`); question 7 changes only the
example under "Closes".

*The one departure, adopted.* The batch runs on `gauge_fixtures.two_basins()`
(the stage B terrain of `test_catchment.py`), not on the valley of
`test_cli_catchment.py` that "The red suites" names: that valley has no
tributary, so the station with a confluence just below it cannot be built
on it. The suite's intent (a synthetic tiled DEM, a river file, the five
stations and the sixth on two grids) is unchanged; "The red suites" now
says so.

*Adopted as the design* (each is now part of the design, as the suite
states it):

1. `reference.agreement(fine, reference, meta) -> Agreement`: our fine
   outline, NVE's polygon (Polygon or MultiPolygon, every part counted), and
   the final window's `RasterMeta`, of which only the origin and the spacing
   are read. `Agreement` is a frozen dataclass: `ours`, `ref`, `both`,
   `nve_in_ours`, `ours_in_nve`, `area_ratio`, `divide_offset_m`, `cell_m`.
2. `reference.classify(agreement, gauge, refused=False) -> (class,
   match_by)`, with `gauge` a `catchment.GaugeResult` or None. The class is
   `refused`, `uncertain`, `match`, `close`, `miss`, or None for a station
   with no cause and no reference; `match_by` is `overlap` or `offset` for a
   match and None otherwise. With no gauge (the fallback of PR 5) the
   agreement alone decides. Classify reads the gauge's joined `causes`
   (next item), not `Sensitivity.causes`, so `direction` counts.
3. `GaugeResult.causes`: the sensitivity's causes in their order, then
   `direction` when `direction_ok` is false. This is the move "The window
   loop" already announced; `catchment --rivers` reads it instead of joining
   the list itself, and writes the same output.
4. A swing of exactly 0.05 raises no cause: `swing <= SWING_MAX` is well
   posed, as "Sensitivity", step 5, says.
5. `summarise(rows) -> Summary`, a frozen Pydantic model, read by
   attribute from the rows, and its JSON keys: `stations`; `classes` (all
   five); `match_by` (both); `scored`, holding `area_ratio`, `nve_in_ours`
   and `ours_in_nve`, each with `min`, `p10`, `p25`, `p50`, `p75`, `p90`,
   `max`, numpy's default percentile (linear between the closest ranks);
   `by_size`, bands `under 10`, `10-100`, `100-1000`, `over 1000` (km²; a
   band's lower bound is in it), decided by NVE's polygon area, or by our
   fine area without a reference, each band with `stations`, `classes`,
   `uncertain_share` and the three measures; `by_tiles`, groups `1`, `2`,
   `3-4`, `5+`, with the same contents; `uncertain_causes`;
   `refusal_causes`; `known_refusals` with `count`, `stations` (in row
   order) and `line`.
6. The known-refusal line, word for word as "Agreement and classes" gives
   it, with N the count.
7. `StationResult` field names: `station_class` (a field cannot be called
   `class`; the CSV column is headed `class`), `refusal_message`,
   `grid_tiles` (one tile from each grid, `mixed_grid` only, else None),
   `distance_m`, `uncertainty_m`, areas in km² (`fine_area_km2`,
   `reference_area_km2` for NVE's polygon, `nve_area_km2` for the station
   layer's), `tiles` (how many of the repository's tiles have a node box,
   from `footprints()`, that meets our fine outline), `windows` (how many
   windows), `causes` a tuple; a refused row has None in every field it
   could not have (a refusal by `delineate` keeps its placement; one by
   `place` has none).
8. `run_batch`: references are in the river file's CRS; `sink.catchment`
   is called before the station's row, never for a refused station; a
   number in `only` that is not in the list raises `ValueError` naming it
   before any station runs; `only` keeps the file's order.
9. Only `CatchmentError` is a refusal. Any other exception, a plain
   `ValueError` included, stops the batch, and the rows already given to
   the sink stay given.
10. `MixedGridError.tiles` is the two names in the message's order;
    `catchment._plan` raises `MixedGridRefusal` from it (`raise ... from`),
    with the same `tiles` and words. Other refusals of the plan are not
    `MixedGridRefusal`.
11. `io.geojson.catchment_geojson(polygon, crs, properties) -> bytes`:
    UTF-8 JSON of a `FeatureCollection` with the `crs` member `{"type":
    "name", "properties": {"name": crs}}` and one `Feature`, the polygon's
    exterior ring; `json.dumps` as 22 called it, so `rasputin catchment`'s
    files keep their bytes. `cli.py` no longer builds the document, and
    `io/geojson.py` opens no file.
12. `results.csv`: a header, one row per station in file order, the columns
    `StationResult`'s fields in their order; a tuple is joined by `;`, None
    is an empty cell, any other value is written as `csv` writes it.
13. `summary.json` is `Summary`'s JSON plus `river_copies_dropped`, which
    the command adds (the batch never sees the river file).
14. stderr: the copies line of "Ola's rulings" (PR 3's green step), then
    one line per station with its number, its name and its class (the word
    `refused` for a refusal; for a scored station also the two overlaps and
    the area ratio, as "The batch" says), then the summary. A stations or
    rivers file without a `crs` member is refused naming its option
    (`--stations`, `--rivers`) and the missing CRS, and nothing is written.

*Added where the suite is silent* (new design; red amendment below):

- **Fixed lists are written in full.** `classes`, `match_by`,
  `uncertain_causes` (`swing`, `downstream_unread`, `chain_not_draining`,
  `chain_end_open`, `direction`) and `refusal_causes` (`mixed_grid`,
  `no_river`, `other`) always hold every key, in that order, 0 where none.
  Causes are counted over `uncertain` rows only.
- **Empty groups.** A measure over no scored station is `null` in place of
  its seven-value object (in `scored`, a band or a tile group).
  `uncertain_share` is the band's `uncertain` rows over its rows that are
  not `refused` (a refused station was never assessed, so it must not
  dilute the share), and `null` when there are none. A band's `stations`
  counts every row in it, refused ones included. A row with no area (no
  reference and no catchment) is in no band, and a row with `tiles` None in
  no tile group, so the bands' counts may sum to less than `stations`.
- **The known-refusal line for one station**: "1 station refused because
  its window selects tiles on two different grids, which rasputin does not
  combine, and neither grid covers the window alone (a known refusal, not a
  failure)". With none, `line` is `null`.
- **The reference file's CRS must be the river file's** (equal as
  `pyproj.CRS`, through `crs.parse_crs`), since `run_batch` takes the
  references in that CRS; otherwise `station-catchments` refuses, naming
  `--reference` and both CRSs, and writes nothing. No reprojection of
  polygons.
- **`cell_m` on unequal spacing** is `sqrt(dx × dy)`, the side of a square
  of one cell's area, which is what the offset (an area over a length) is
  measured against; on DTM10 it is 10 m. Not tested: every DEM of this
  increment has square cells.
- **`StationResult` is a frozen dataclass** with slots, as `GaugeResult`
  is; the CSV's columns come from `dataclasses.fields` in order. Its
  `station_class` is typed `Literal["refused", "uncertain", "match",
  "close", "miss"] | None`.

*The swing is not split.* "Agreement and classes" said the summary splits
the `swing` cause by its largest step's position, a confluence or a flat
floor or lake, but gave no rule, and the counts alone have none: a flat
floor also gives a step where the flat's nodes join the chain ("Sensitivity",
after step 5), so a step does not tell a confluence from a flat. A rule would
need a new threshold with nothing measured behind it. **`swing` stays one
count in `summary.json`**; each row carries `largest_step`,
`largest_step_at_m` and `lake`, and the acceptance README tells the causes
apart station by station (acceptance step 4 already asks for that).

**Red amendment** (`@tester`, lean, one commit before green):
in `test_reference.py`, (a) a summary over rows that lack some causes and
refusal causes still has all five and all three keys, in order, the absent
ones 0; (b) a band whose rows are one `uncertain` and one `refused` has
`stations` 2 and `uncertain_share` 1.0, a band with only refused rows has
`uncertain_share` null and the three measures null; (c) a row with neither
area is in no band, and one with `tiles` None in no tile group; (d) one
`mixed_grid` refusal gives the singular line above, and none gives `line`
null. In `test_cli_station_catchments.py`, a `--reference` file in another
CRS (EPSG:32633) is refused, naming `--reference`, and nothing is written.
The amendment is `7cd56a4`; like the rest of the red suite, it failed on the
missing modules until green.

**PR 4's green step: `@developer`'s choices, read by `@architect`
(2026-10-05).** Green `1c39ef7`, on red `2b9b39f` and `7cd56a4`. Not Ola's
rulings. Each choice the design left open is adopted as the design, or
changed with a red and a green ask below.

*The area, as built (question 8, default).* `fine_area_km2` is `ours` (the
lattice nodes strictly inside our fine outline) × cell area, with or
without a reference, so the size band of a station with no reference and
the area ratio both read the same count as the overlaps;
`reduced_area_km2` is the reduced outline's shapely area;
`reference_area_km2` is NVE's polygon's shapely area; `nve_area_km2` is the
station layer's own figure, copied. `area_ratio` = `ours` × cell area /
`reference_area_km2` ("Agreement and classes" says why it is half a cell
above the traced outline's area). The catchment file keeps 22's
`fine_area_m2` (the traced outline's shapely area) and `reduced_area_m2`,
so `results.csv`'s `fine_area_km2` and the file's `fine_area_m2` differ by
half a cell; the README of the acceptance run says so.

*Adopted:*

1. **`StationResult`'s columns**, in the order of
   `catchment_batch.py`: the station and its class (`station`, `name`,
   `class`, `match_by`, `refusal_message`, `refusal_cause`, `grid_tiles`);
   the placement; the gauge and the burn (the burn's `end_extended_m` and
   `end_closed` with the gauge, not after the sensitivity as "The batch"
   lists them); the sensitivity and the joined `causes`; the catchment and
   the agreement; `tiles`, `windows`, `seconds`. Three columns beyond the
   design's list: `chain_nodes` (`GaugeResult`'s own field, which the row
   takes whole), `reduced_area_km2` and `divide_offset_m` (both in "The
   batch"'s words, "fine and reduced area" and "the agreement numbers",
   now with names).
2. **A refused row keeps `reference_area_km2`, `nve_area_km2` and
   `seconds`.** They are known before the refusal, and the size band of a
   refused row needs the reference area (the red amendment's band with an
   `uncertain` and a `refused` row counts both).
3. **The catchment file is written when its row arrives.** The batch still
   calls `sink.catchment` before `sink.row`; the directory sink holds the
   catchment until the row, because the file's properties come from the
   row. The properties: `station`, `name`, `series`, then 22's (`nodes`,
   `fine_vertices`, `fine_area_m2`, `reduced_vertices`, `reduced_area_m2`,
   `outline_tolerance_m`, `windows`), then the row's fields from
   `placed_on` to `causes` (the placement, the gauge, the burn, the
   sensitivity). No `seed` or `seed_crs` (the station's position is in the
   stations file) and no agreement numbers or class (they are in
   `results.csv`). A batch stopped by a bug writes no file for the station
   in flight, which never got its row.
4. **The stderr lines.** Per station, `<station> <name>: <word>`, the word
   being the class, or "well defined, no reference" for a station with no
   cause and no reference polygon; a refusal adds `: <message>`, an
   `uncertain` row its causes in brackets, and a row with both overlaps
   adds "; NVE's in ours 98.4 %, ours in NVE's 99.2 %, area ratio 0.984"
   (one decimal, three for the ratio). At the end, `summary: N stations:
   refused a, uncertain b, match c, close d, miss e`, the known-refusal
   line when there is one, and the output directory on stdout.
5. **The reference CRS refusal** names `--reference` and says "the
   reference file's CRS, X, is not the river file's, Y; polygons are not
   reprojected".
6. **`nve_in_ours` and `ours_in_nve` are 0.0 when their divisor is 0.**
   `ours` is never 0 (the outline holds the placed node); `ref` is 0 when
   NVE's polygon holds no node. 0.0 fails both overlap bars, and the offset
   test still decides, which is the right test for a polygon smaller than
   the lattice resolves. The third divisor is not covered: see change (b).
7. **`MixedGridError.tiles` defaults to `()`.** `dem_input.py:248`
   re-raises a mosaic error as `type(exc)(message)` with more words, which
   calls the constructor with the message alone; the mesh path reads no
   tiles, and the catchment path, the only reader, gets them from
   `_covering` directly.

*Changed* (`@tester` red, then `@developer` green; one round, lean):

- **(a) An unknown `--only` is refused by the command, naming `--only`,
  before anything is written.** Today `run_batch`'s `ValueError` ("not in
  the stations file: 9.9.9") comes out as a plain `Error:` line, without the
  option's name, after `--out-dir` has been created. *Red*
  (`test_cli_station_catchments.py`): `--only 9.9.9` exits non-zero, the
  output names `--only` and `9.9.9`, and `--out-dir` does not exist.
  *Green* (`cli.py`): check the numbers against the stations read, before
  `target.mkdir`, as `typer.BadParameter(..., param_hint="--only")`.
  `run_batch` keeps its own check (red-step item 8) for callers other than
  the command.
- **(b) A reference polygon with no area is refused when the file is
  read.** `agreement` divides by NVE's polygon's area and perimeter: a
  zero-area ring raises `ZeroDivisionError`, and an empty polygon a
  `ValueError` from its NaN bounds ("cannot convert float NaN to
  integer"), both checked on `1c39ef7`. Either stops a batch, possibly
  hours in, at that station. *Red* (`test_station_set.py`): a reference
  feature whose polygon has a zero-area ring, and one with empty
  coordinates, each make `read_references` raise `ValueError` naming the
  station's number. *Green* (`io/station_set.py`): `read_references`
  refuses a geometry that is empty or has area 0, naming the station.
  The command already refuses a `ValueError` there naming `--reference`.
- **(c) `results.csv` is written row by row.** Today the table is written
  only when the batch ends, so a bug at the 139th station of a run of
  hours loses the 138 rows (stderr keeps them, in less detail). The
  directory sink opens `results.csv` when it is made, writes the header,
  and writes and flushes each row as it arrives; `summary.json` stays at
  the end. *Red* (`test_cli_station_catchments.py`): with
  `catchment_batch.delineate` replaced by a stub that delineates the first
  station and raises `RuntimeError` at the second, the command fails, and
  `results.csv` holds the header and the first station's row, and
  `summary.json` does not exist. *Green* (`cli.py`): the sink as above;
  the bytes of a completed run's `results.csv` do not change.

Changes (a) to (c) are red `6ebff1b` and green `33b1f2d`.

**PR 4's code review, round 1: the suggestions, ruled by `@architect`
(2026-10-05).** Not Ola's rulings. The round's one blocker,
`project_structure.md`, waits on Ola (see the round's record under
"Review"). Of its four suggestions, three become changes (d) to (f),
`@tester` red then `@developer` green, one round, lean; the fourth is a
test-only commit. Together about 20 production lines, so PR 4 stays under
700 (about 560 net) and is not split.

- **(d) A river file not in the DEM's CRS is refused before anything runs,
  naming `--rivers`, and before `--out-dir` is created.** Before (d),
  `station-catchments` checked only the reference file against the river
  file's CRS. The DEM's CRS is checked per station, inside `delineate`
  ("the river reach must be in the DEM's CRS"), so every placed station
  becomes a `refused` row with cause `other`, the coverage test compares
  station points and tile boxes in two CRSs, and the run takes its whole
  length to say one thing. `catchment --rivers` refuses the same file at
  once (`_placed` in `cli.py`). The references need no check of their own:
  they are refused unless in the river file's CRS, so once the river file is
  in the DEM's CRS, so are they. *Red* (`test_cli_station_catchments.py`):
  a river file and a reference file, both in EPSG:32633, over a DEM in
  EPSG:25833: non-zero exit, the output names `--rivers` and both CRSs, and
  `--out-dir` does not exist. (`test_catchment_batch.py`): `run_batch` with
  `segments_crs` EPSG:32633 over that DEM raises `ValueError` naming both
  CRSs, and the sink has received nothing. *Green*: one rule in one place,
  `catchment.py`'s `check_reach_crs(crs, repository) -> None`, raising
  `ValueError("the river file's CRS, X, is not the DEM's, Y")` when
  `parse_crs(crs)` is not the first tile's CRS (the comparison `_placed`
  makes today). `_placed` calls it, mapping the error to
  `typer.BadParameter(..., param_hint="--rivers")` as now;
  `station_catchments` calls it the same way right after the repository is
  opened, before the request is built and before `mkdir`; `run_batch`
  calls it before its first station, beside its `--only` check, for callers
  other than the command. A DEM whose tiles are in several CRSs keeps
  `delineate`'s own refusal ("the tiles are in N CRSs").
- **(e) A failed write names the output, not `--dem`.** The
  `except OSError` around `run_batch` in `station_catchments` exists for
  the DEM's reads, but the directory sink writes inside `run_batch`, so a
  catchment file that cannot be written reads "Invalid value for --dem:
  cannot read .../out/1.140.0.geojson: [Errno 21] Is a directory". The
  `results.csv` open and the `summary.json` write sit outside every
  `except`, so their failure is a traceback. *Red*
  (`test_cli_station_catchments.py`, three cases, each with a directory made
  beforehand where the command writes a file: `<out-dir>/<station>.geojson`,
  `<out-dir>/results.csv`, `<out-dir>/summary.json`): non-zero exit; the
  output names `--out-dir` and the path; it contains neither `--dem` nor
  "cannot read"; no traceback. *Green* (`cli.py`): every write the command
  makes (making `--out-dir`, the `results.csv` open, its header, each row's
  write and flush, its close, each catchment file, `summary.json`) turns an `OSError` into
  `typer.BadParameter(f"cannot write {path}: {exc.strerror or exc}",
  param_hint="--out-dir")`, through one small helper. `BadParameter` is
  neither an `OSError` nor a `ValueError`, and `run_batch` stops on any
  exception but `CatchmentError`, so it reaches the user unchanged; the
  `except OSError` around `run_batch` is then about reads alone.
- **(f) A station with no name prints no stray space.** The stderr line is
  `f"{row.station} {row.name or ''}: {word}"`, so a station with no name
  prints `1.2.0 : match`. *Red* (`test_cli_station_catchments.py`): a
  station with no `name` gives a stderr line starting `<station>: `. *Green*
  (`cli.py`): the name and its space only when the name is set and not
  empty.
- **Test-only: plain imports.** The `importlib` fixtures `cb`, `gj` and
  `ref` in `test_catchment_batch.py`, `test_io_geojson.py` and
  `test_reference.py` let the red suite be collected before the modules
  existed; now they exist, they become plain imports. `@tester`, in its own
  commit ahead of the red one, so the suite passes unchanged at `33b1f2d`
  in between. `test_mosaic.py`'s `importlib` fixtures are older than PR 4
  and are left.

The test-only commit is `e8e93d0`; changes (d) to (f) are red `bc78be3`
and green `cabff74`. They came to 31 production lines against "about 20"
(`cli.py` 25, of which 7 are ruff splitting the `tin_engine.catchment`
import over several lines once `check_reach_crs` joined it;
`catchment.py` 4; `catchment_batch.py` 2), so PR 4 is 572 net, not about
560, and still not split.

Code review round 2 found that (e) missed one write: making `--out-dir`
(`target.mkdir(exist_ok=True)`), so an `--out-dir` that is an existing
file, or under a parent that cannot be written, ended in a
`FileExistsError` or `PermissionError` traceback. Two more red cases
(`44dc954`: an `--out-dir` that is a file; a read-only parent, skipped when
run as root) and green `f952f2d` put the `mkdir` under the same helper;
the green also moves the `results.csv` header write and its close under
it, which no test exercises. PR 4 is then 573 net.

**Lake gauges: Ola's ruling (2026-10-05).** PR 2's chain fails on lake
gauges such as Narsjø and Tingvatn (catchments under 1 km² against NVE's
119 and 272 km², "Placement figures"). Ola, 05:07 UTC: "What we did in
Bygdin, which was a lake, worked really well. Can't we do the same now?",
and 05:11 UTC: "Fold into PR4". The design is "Lake gauges", below; the evidence is "NVE's
lakes: Innsjødatabasen"; questions 9 and 10 for Ola come from it.

**Lake gauges' green step: `@developer`'s choices, read by `@architect`
(2026-10-05).** Red `a0831c9` (`@tester`), green `fff6ac9` (`@developer`,
no test file touched); `@developer` reports the full suite on a rebuilt
`_core` at exit 0, 4989 passed, 17 skipped. Counted by PR 2's round-2 rule,
the whole of PR 4 over `9e666f4..fff6ac9` is 691 net (744 added, 53
removed), 8 lines of room under 700: `cli.py` 201, `catchment_batch.py` 195,
`reference.py` 187, `io/station_set.py` 32, `gauge.py` 27, `io/geojson.py`
21, `catchment.py` 13, `fetch/nve.py` 10, `mosaic.py` 5, `burn.py` 0. The
lake work alone, over `48d1315..fff6ac9`, is 118 net (150 added, 32
removed) against "about 100": `cli.py` 28, `io/station_set.py` 28,
`gauge.py` 27, `catchment_batch.py` 22, `fetch/nve.py` 10, `reference.py`
3. So the PR 4b seam (under "New and changed files") does not fire, but
**any change that takes PR 4 to 700 net or more forces it** (8 net lines
of room are left, since the rule is "under 700"): PR 4 is then
published as it stood at `48d1315` and the lake work becomes PR 4b. A
fix that a later review asks for in the lake work counts toward the 700
like any other line. `@developer`'s choices beyond the design, each
checked against the code and adopted:

- **Several lakes contain the point**: `lake_seed` keeps them all
  (`LakeSeed.lakes`, all seeded together) and, on the `lake_line` rule,
  measures the gap to the nearest of them (`min` of the distances).
- **The row's `lake_number` and `lake_name` come from the first lake** in
  `LakeSeed.lakes`, which is the lakes file's order.
- **A lake row that is refused keeps `seeded_by` = `lake`**, with its
  `lake_rule`, `lake_number`, `lake_name` and `lake_distance_m`: the seed
  was chosen before `delineate` refused, and the row says which path ran.
- **`fetch/nve.py`'s `LAKE_ENVELOPE_HALF` (100 m) is not in its
  `__all__`**, while `ENVELOPE_HALF` is; only the fetch reads it, and
  the tests keep their own copy in `nve_fixtures.py`.

**Questions 5 to 10: Ola's ruling (2026-10-05).** Ola, quoted by the main
session: "We keep increment 29 as is, and reassess after the full run." So
every question under "Questions for Ola" is ruled: 1 to 4 as recorded above,
and 5 to 10 kept as built, which is each one's default. A gap in the DEM on
the river's burnt path refuses the station (5); a river point just upstream
of the gauge that does not drain through it makes the station uncertain (6);
the example meshes Atnasjø, `2.32.0` (7); the table's area is the DEM point
count times a cell's area (8); Femundsenden, with no river line near it,
stays refused until the nearest-stream fallback, PR 5 (9); a gauge a little
below a lake's outlet keeps the river method (10). All ten are looked at
again once the acceptance run over all 140 stations is in.

**Deferred from PR 4's code review, round 4.** PR 4 has 8 production lines
of room left under the 700-line ceiling, so the two non-blocking
suggestions wait for a later PR: `read_lakes` refusing a malformed lake
(too few points in a ring, a `vatnlnr` that is not a number) with the
underlying library's message, which does not name the lake, and reading a
string `"0"` as lake number 0 rather than none; and `Lake` moving out of the
reader module `io/station_set.py`, next to `Reach`, so that `gauge.py` does
not depend on a reader.

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
  gauges; PR 2's code on the 2026-10-05 data finds 48, with the nearest
  placement and with the tiered one alike:
  `docs/benchmarks/2026-10-05/nve-placement/README.md`, "Counts that
  differ"). For 16 stations, a second branch (another `elvid`) lies within
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

### NVE's lakes: Innsjødatabasen (measured 2026-10-05, for "Lake gauges")

Measured by `@architect` with throwaway scripts in the session's scratchpad
(not committed), on the fetch of 2026-10-05 in `../rasputin_data/nve_hrd`
and PR 2's code at this branch's head.

- **What 22 seeded Bygdin with.** CORINE Land Cover 2018, class 512 (water
  bodies), feature `fid` 54101 of `corine2018_dtm10_utm33.gpkg`, a local
  extract (22's "What the data says about Bygdin"). CORINE's minimum
  mapping unit is 25 ha with a 100 m minimum width
  (`docs/research/raster-to-vector.md`, "CORINE's minimum mapping unit"), so
  it leaves out small lakes: 7 of the 48 lakes below are under 0.25 km²
  (the smallest 0.023 km²). It is also not fetched by rasputin.
- **NVE's lake database is a layer like ELVIS.** `Innsjodatabase2`,
  `https://kart.nve.no/enterprise/rest/services/Innsjodatabase2/MapServer`,
  layer 5 `Innsjodatabase` (polygons), EPSG:25833, `maxRecordCount` 2000, on
  the same host and query interface as layers 0, 38 and ELVIS (read with
  `?f=json` on the service and on layer 5, 2026-10-05). 267,194 features
  (`returnCountOnly`). Its fields include `objectid`, `vatnlnr` (the
  national lake number), `navn`, `areal_km2`, `hoyde`, `magasinnr`,
  `vassdragsnr`, `kommune` and `globalid`. The service description says
  every lake larger than 2500 m² has a unique national serial number.
  Geonorge's record "Innsjødatabase" (uuid
  `823b8639-9a49-41bf-8571-3608435eb149`, read through
  `kartkatalog.geonorge.no/api/getdata/` on 2026-10-05): "Åpne data" under
  NLOD, scale 1:20,000, NVE as the organisation.
- **ELVIS's lake number is the lake database's.** For each of the 140
  stations, layer 5 was queried with the station point ± 1000 m as the
  envelope (`outFields=objectid,vatnlnr,navn,areal_km2,vassdragsnr,hoyde,magasinnr`,
  `outSR=25833&f=geojson`; 140 queries, one at a time, 34 s, 8.3 MB, no
  answer truncated). PR 2's tiered placement (`gauge.place`, the defaults)
  puts 48 stations on a lake line; 47 of those lines carry a lake number
  (`vatnlnr`, not null and not 0), and for all 47 a polygon with that number
  was in the answer, and the mapped position `P` lies inside it. The 48th,
  `97.1.0` Fetvatn, is on a lake line with no number; its `P` lies inside
  Fitjavatnet.
- **Where the stations are.** 18 station points lie inside a lake polygon
  (16 of them placed on a lake line, and `62.18.0` Svartavatn and `191.2.0`
  Øvrevatn placed on a river line), 55 within 10 m of one, 64 within 30 m.
  No station is within 30 m of two lakes (one is within 50 m of two, three
  within 100 m). Of the 48 placed on a lake line, the station's distance to
  the lake its `P` lies in has quartiles 0.0, 1.3 and 6.0 m; two are
  farther than 30 m: `12.197.0` (52 m) and `22.16.0` Myglevatn ndf. (460 m:
  "ndf." is "below"; its line is a pond's, Tveitevatnet, while Myglevatnet is
  87 m away). `127.11.0` Veravatn is placed on a pond's line 118 m away but
  lies inside Veresvatnet. Femundsenden (`311.4.0`, no river line within
  500 m) is 14 m from Femunden. 24 stations placed on a river line have a
  lake on their reach upstream of `P` (within `reach_up`, 1000 m): gauges
  below a lake outlet.
- **The lake seed works on them.** 22's path (`catchment.delineate` with
  `lakes`, the lake's polygon, and a seed point inside it), run on the 48
  stations the rule of "Lake gauges" picks, compared with NVE's polygons by
  `reference.agreement` and classed by "Agreement and classes": **41
  `match`** (all by overlap), 3 `close`, 1 `miss`, 3 refused on tiles of
  two grids (`203.2.0`, `213.2.0`, `191.2.0`); 299 s in all, the longest
  56 s (`62.5.0`, 1091 km²). Over the 45 scored, NVE's in ours has 10th
  percentile 96.7 % and median 98.8 %; ours in NVE's 96.2 % and 98.8 %;
  the area ratio 0.989 and 1.000. Narsjø: 98.2 % / 97.5 %, ratio 1.007
  (PR 2's path: 0.0012 km²); Tingvatn: 99.4 % / 99.5 % (0.48 km²); Atnasjø:
  98.5 % / 99.1 % (98.4 % / 99.2 % by PR 2's path). Of the 8 lake gauges
  PR 2's tiered placement ran in full (`full_runs_tiers.json`), 5 got
  catchments under 1 km²; the lake seed gives 4 of them `match`. The fifth
  is the `miss`: `83.2.0` Viksvatn (Hestadfjorden), whose point lies inside
  Hestadfjorden (watercourse number `083.C2`; the station's is `083.C1`),
  19.3 km² against NVE's 508 km², 99.3 % of ours inside NVE's: a lake
  whose own catchment is a small part of the gauge's, so less area, not
  more (why the station point lies in it is not looked at). The `close`
  rows: `16.66.0` (93.6 % / 98.7 %, 6.5 km²), `35.9.0` (98.7 % / 93.3 %)
  and `26.29.0` Refsvatn (99.3 % / 84.0 %, ratio 1.18: 18 % more area than
  NVE's; where the extra lies is not looked at, and the acceptance must say).

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
  (`elvenett`); `Innsjodatabase2` layer 5 (`Innsjodatabase`, the lakes, PR 4,
  "Lake gauges"). HydAPI is not used and no API key is stored anywhere.
- **Licence.** NLOD; "Kilde: NVE" is in `NOTICE.txt` of every fetch, in the
  packaged CSV's header (each tested) and in `NOTICE.md` (checked by hand by
  `@reviewer`, block "PR 3, after PR #178's CI" above). NVE disclaims liability for errors
  in the data and their use; the README of the acceptance says so too.
- **Committed**: the 140-row HRD list (four columns, from a PDF) and the
  acceptance's `results.csv` (NVE's areas as numbers). **Fetched** into
  `../rasputin_data/nve_hrd/` and never committed: polygons, station points,
  river lines, lake polygons, manifest.
- **Only the fields the design needs, by explicit allow-list per layer, in
  every request (never `outFields=*`).** Layer 0: `stasjonnr`, `stasjonnavn`,
  `totalt_feltareal_km2`, `stasjonstatus`, `vassdragsnr`, `elvenavnhierarki`.
  Layer 38: `stasjonnr`, `nedborfeltaareal_km2`, `oppdateringsdato`,
  `objectid`. Layer 2: `objectid`, `objekttype`, `strekninglnr`, `elvid`,
  `vassdragsnr`, `elvenavn`, `vatnlnr`. Lake layer 5: `objectid`, `vatnlnr`,
  `navn`, `areal_km2`. **Never collected**: `stasjoneier` (the
  owner), ELVIS's `oppdatertav` (editor ids, some look like personal
  initials), `globalid`, layer 38's discharge normals, and the lake layer's
  municipality fields (`kommnr`, `kommune`), reservoir fields and `globalid`.
- **Query volume, kept small.** Only `fetch-stations` uses the network: about
  8 batched queries for layers 0 and 38 (40 stations each), 140 ELVIS
  envelope queries and 140 lake envelope queries (PR 4), sent one at a time, with 23a-2's retries and back-off. A
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

**Lake gauges (2026-10-05).** *Legacy*: nothing to carry over.

```
$ git grep -liE "innsj|vatnlnr|lake" legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/rasputin/material_specification.py
legacy-archive:legacy/rasputin/triangulate_dem.h
legacy-archive:legacy/rasputin/web_visualize.py
```

Every hit is about drawing lakes: `extract_lakes` splits a mesh's faces into
lake and terrain for the web viewer, and the land-cover readers give lake
classes a material. None seeds a catchment. *Literature*: the method is
increment 22's lake seed (its "The seed" and "Prior art"), every DEM node
inside the lake polygon a seed of one Priority-Flood labelling, which 22
measured on Bygdin (99.12 % of NVE's nodes in ours). What differs here: the
polygon is NVE's lake database rather than CORINE, and the lake is chosen by
the station's position and the mapped river, not by a point the user gives.
**HydroLAKES** (Messager, Lehner, Grill, Nedeva and Schmitt 2016,
"Estimating the volume and age of water stored in global lakes using a
geo-statistical approach", *Nature Communications* 7:13603,
doi:10.1038/ncomms13603, checked on Crossref 2026-10-05; the data page
`https://www.hydrosheds.org/hydrolakes` known only from a search summary)
ties each lake to the HydroSHEDS river network by one pour point and reads
the lake's upstream area there. Seeding the whole lake needs no pour point,
so no rule has to choose one; how HydroLAKES chooses its pour points was not
read.

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
                 layers 0 and 38, ELVIS layer 2 and the lake layer 5
                 through fetch/http.py, newest polygon per station, writes
                 D/stations.geojson, D/reference.geojson, D/rivers.geojson,
                 D/lakes.geojson (PR 4), D/NOTICE.txt,
                 D/manifest.json (URLs, date, sha256)

rasputin station-catchments --dem ... --stations D/stations.geojson
        --rivers D/rivers.geojson [--reference D/reference.geojson]
        [--lakes D/lakes.geojson] --out-dir O
                                                                    [offline]
  cli.py         paths stop here: io/station_set.py and io/rivers.py read the
                 files into Station and RiverSegment models and reference
                 polygons; repository_for(dem)
  catchment_batch.run_batch(request, repository, stations, segments, references, sink)
     for each station, in file order, one at a time:
       gauge.place(Gauge(station), segments) -> Placement | None   [pure, shapely, no DEM]
       gauge.lake_seed(gauge, placement, lakes) -> LakeSeed | None  [pure, PR 4]
       a lake seed: catchment.delineate(CatchmentRequest(seed=lake_seed.point,
                        lakes=lake_seed.lakes), repository)          (22's path)
       otherwise:   catchment.delineate(CatchmentRequest(seed=station, seed_crs=crs,
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
   segment is a lake centreline: 46 of 139 nearest lines; 48 by PR 2's code
   on the 2026-10-05 data, under either placement, per the placement
   figures' README, "Counts that differ"), `confluence_near`
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
    causes: tuple[str, ...]    # PR 4: the sensitivity's causes, then "direction" if not direction_ok
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

### Lake gauges: the lake is the seed (PR 4)

Ola's ruling of 2026-10-05 ("Ola's rulings", last block). A gauge on a lake
is seeded as increment 22 seeded Bygdin: **every DEM node inside the lake's
polygon is a seed**, through 22's own path (`CatchmentRequest` with `lakes`
and `lakes_crs`, `reach` None), so no chain is burnt across the lake and no
sensitivity is read along it. `catchment.py` does not change. The evidence is
"NVE's lakes: Innsjødatabasen" (41 `match` of the 48 stations the rule
below picks; of the 8 lake gauges PR 2's path ran in full, 2 were `match`).

**Which stations: `gauge.lake_seed`** (pure, shapely, no DEM; beside
`place`, whose `Placement` it reads):

```python
LAKE_GAP_M = 30.0  # metres in the river file's CRS; see below

@dataclass(frozen=True, slots=True)
class Lake:                  # io/station_set.py; one polygon part
    number: int | None       # NVE's vatnlnr; None for null or 0
    name: str | None         # navn
    polygon: Polygon         # in the lake file's CRS

class LakeSeed(BaseModel):   # frozen; arbitrary types allowed
    rule: Literal["inside", "lake_line"]
    point: tuple[float, float]  # the station (inside) or P (lake_line), river file's CRS
    lakes: tuple[Lake, ...]     # the lakes whose polygon contains `point`
    distance_m: float           # from the station to the lake; 0 inside

def lake_seed(gauge: Gauge, placement: Placement | None, lakes: Sequence[Lake],
              *, gap: float = LAKE_GAP_M) -> LakeSeed | None
```

1. **`inside`**: the station point lies inside a lake polygon
   (`Polygon.contains`, strict, as 22's `_lake` tests it). The point is the
   station. This holds with or without a placement.
2. **`lake_line`**: otherwise, the placement exists, its chosen line is a
   lake line (`Placement.lake`), the mapped position `P` lies inside a lake
   polygon, and the station's distance to that polygon is at most `gap`
   (`<= gap + 1e-6`, the corridor's slack). The point is `P`. The lake is
   found by geometry, not by the lake number, so a lake line without one
   (`97.1.0`) still finds its lake; the number is only reported.
3. Otherwise `None`: PR 2's river path, unchanged, including the `no_river`
   refusal when there is no placement.

`lakes` in `LakeSeed` are all the lakes containing `point`: one with NVE's
data (none of the 140 stations is within 30 m of two). Two, from a user's
file with overlapping polygons, reach 22's `_lake`, which refuses ("the
seed point ... is in 2 lakes; give one", a `LakeError`, so a `CatchmentError`):
a `refused` row with cause `other`. A MultiPolygon lake is split into its
parts by the reader, each a `Lake` with the same number and name (22: "a
multipolygon contributes the part containing the point").

**`LAKE_GAP_M` = 30 m** is `U`'s floor ("Placing the gauge", step 3), three
cells of DTM10, set before the lake seeds were run, not tuned to their
agreement. It assumes a metric CRS of metres (EPSG:25833) and was checked on
the 140 HRD stations only: at 30 m the rule picks 48; at 52 m it would add
`12.197.0`, and at 460 m `22.16.0`, whose line is another lake's ("NVE's
lakes"). The gap applies to `lake_line` only; `inside` needs none.

**Why it stays position-faithful** (Ola's residual-inflow goal: never
area-maximising). No rule reads an area, a count or the reference: the lake
is chosen by containment and one distance. A gauge on a lake measures the
lake's outflow, and every lake node drains to the outlet, so the gauge's
catchment is everything that drains into the lake: the set the seed labels,
the same from any position on the lake. The seed can add area the gauge
does not see in two ways, both measurable: (a) a polygon that reaches past
the real outlet (22's Bygdin: CORINE's polygon ended at x 168460, about
370 m past the dam at x 168087), which adds what drains into that stretch;
(b) a station on a river flowing into the lake whose nearest line in its tier
is the lake's, which would add the lake's other inflows. The tiers make (b)
unlikely (44 of the 48 lake lines carry the station's own watercourse number,
the other 4 its river's name), and the acceptance lists every lake row whose
ours-in-NVE's is under 95 % with the cause (`26.29.0` is one; "NVE's
lakes"). The seed misses area in two ways, which are deficits, never extra
area: a station outside the polygon (at most 30 m) on the river below the
outlet loses what drains into those metres; and a station inside a lake that
is not the one it gauges gets that lake's smaller catchment (`83.2.0`). **A
gauge farther below the outlet** than the gap, or inside no lake, stays on
the river path: its catchment is the placed node's, which holds the lake
through the DEM's own drainage and nothing below the gauge. 24 stations have
a lake on their reach within `reach_up` above `P`; the four of them in PR 2's
full runs (`2.633.0`, `55.4.0`, `83.12.0`, `101.1.0`) agree with NVE's
polygons to at least 97.1 % both ways. Seeding the lake together with the
river from its outlet down to such a gauge is not built (question 10).

**What a lake row has, and has not.** No chain, no burn and no sensitivity:
`Catchment.gauge` is None, so the row's gauge, burn and sensitivity columns
are empty, `causes` is empty, and `classify(agreement, None)` decides from
the agreement alone (the rule PR 4 already has for a station with no gauge),
so a lake row is never `uncertain`. The sensitivity measures how the area
changes along the river within `U`; on a lake every position gives the same
seed, so the swing is 0 by construction, and the counts PR 2 read along a flat
lake (Narsjø's −33 %) measured the burn's 1 mm-per-node channel, not the
gauge. What remains uncertain, whether the station is on this lake, is
reported (`lake_rule`, `lake_distance_m`), not scored. The window loop and
the memory cap are 22's (the lake's bounds plus the margin; item size + 2
bytes per node).

**The row and the summary.** `StationResult` gains, after `reach_fork` (so
the catchment file's properties "from `placed_on` to `causes`" carry them):
`seeded_by` (`"river"` or `"lake"`; None when `place` refused), `lake_rule`,
`lake_number`, `lake_name`, `lake_distance_m` (None on a river row).
`Summary` gains `by_seed`, the groups `river` and `lake` with a band's
contents (`stations`, `classes`, `uncertain_share`, the three measures). The
stderr line of a lake row adds `(seeded by the lake Narsjøen)` after the
class word, or the lake's number when it has no name, or "its lake" when it
has neither. The lakes file must be in the river file's CRS, as the
references must ("The batch"); otherwise `station-catchments` refuses,
naming `--lakes` and both CRSs, and writes nothing. Without `--lakes` there
is no lake path and the batch is PR 4's as it stands.

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
- **queries the lake layer once per station** (PR 4, "Lake gauges"), by
  envelope: the station point ± 100 m (`LAKE_ENVELOPE_HALF`; any lake within
  `LAKE_GAP_M` of the station meets it, and the service returns whole
  features), `outFields=objectid,vatnlnr,navn,areal_km2`,
  `outSR=25833&f=geojson`, through the same `_url` and `_features` (so a
  truncated answer is refused, naming the station); the ELVIS envelope query
  and this one share one envelope-URL function. Lakes seen from several
  stations are kept once, by `objectid`;
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
  served and dropped by `read_segments`), `lakes.geojson` (PR 4: one feature
  per lake, the geometry as served, Polygon or MultiPolygon, the four fields
  above, in `objectid` order), all with a `crs`
  member, `NOTICE.txt` (the catalogue's credit, "Kilde: NVE", as 23a-2's
  `notice` does), and `manifest.json` (the query URLs, the fetch time in UTC,
  each file's sha256). Deterministic order: the list file's, then `objectid`.
- Pure parts (the query URLs, choosing the newest version, de-duplicating
  segments, building the files' content) are functions with no network,
  tested on canned answers as `fetch_fixtures.py` does; the network call is
  injected.

**`io/station_set.py`** reads the stations and references back:
`read_stations(path) -> (tuple[Station, ...], crs)`, `read_references(path)
-> (Mapping[str, Polygon | MultiPolygon], crs)`, and (PR 4) `read_lakes(path)
-> (tuple[Lake, ...], crs)`, which splits a MultiPolygon into its parts, reads
`vatnlnr` null or 0 as no number, and refuses, naming the feature's
`objectid` (or its index when it has none), a geometry that is not a Polygon
or MultiPolygon, or is empty or has area 0; **`io/rivers.py`**
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
- **Area ratio** = `ours` × cell area / NVE's polygon area (shapely): our
  area **counted on the lattice**, the same `ours` as the overlaps, over
  NVE's exact polygon area (as built in PR 4's green step, `1c39ef7`;
  question 8). It is not the traced fine outline's shapely area
  (`Catchment.fine_area`, the catchment file's `fine_area_m2`): the traced
  ring runs through the midpoints between in- and out-nodes and cuts each
  corner, so `ours` × cell area exceeds it by half a cell: 50 m² on DTM10,
  under 0.01 % of a 1 km² catchment (six random-walk node sets traced with
  `outline.trace`, outer ring against the nodes strictly inside it,
  2026-10-05: 0.5 cell each time). `ours` counts the nodes of filled holes
  and leaves out those of dropped rings, so it can differ from `nodes`.
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
counted (a station can have several): `swing` (one count, not split by
cause; "PR 4's red step"), `downstream_unread`,
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
                    sink: BatchSink,
                    lakes: Sequence[Lake] | None = None) -> Summary  # PR 4, "Lake gauges"
```

- **The path, per station** (PR 4, "Lake gauges"): after `place`, with
  `lakes` given, `gauge.lake_seed(gauge, placement, lakes)`; a `LakeSeed`
  gives `CatchmentRequest(seed=lake_seed.point, seed_crs=segments_crs,
  lakes=tuple(l.polygon for l in lake_seed.lakes), lakes_crs=segments_crs,
  outline_tolerance=...)`, and `None` gives the reach request as before. A
  station with no placement goes to the `no_river` refusal only when
  `lake_seed` is `None`. `lakes` are in the river file's CRS, as the
  references are; the command checks it.

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
                    --rivers FILE [--reference FILE] [--lakes FILE] [--map-radius METRES]
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
`project_structure.md` recommended (its paragraph "The catchment GeoJSON
writer is in `cli.py`"); both commands call it, and PR 4 replaced that
paragraph with one sentence naming the module, beside its `io/` entry.

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
PR 4's green step (`1c39ef7`) came to 530 net (567 added, 37 removed),
counted by PR 2's round-2 rule over `9e666f4..1c39ef7`: `reference.py` 184,
`catchment_batch.py` 171, `cli.py` 140 net, `io/geojson.py` 21,
`catchment.py` 9, `mosaic.py` 5, `burn.py` 0 (a docstring). That is 58 %
over its 335 and 48 lines past the margin (482), and under 700. With
changes (a) to (c) (red `6ebff1b`, green `33b1f2d`) it is 541 net (580
added, 39 removed) over `9e666f4..33b1f2d`: `reference.py` 184,
`catchment_batch.py` 171, `cli.py` 147, `io/geojson.py` 21,
`catchment.py` 9, `mosaic.py` 5, `io/station_set.py` 4, `burn.py` 0;
59 past the margin and under 700. With changes (d) to (f) (red `bc78be3`,
green `cabff74`) it is 572 net (614 added, 42 removed) over
`9e666f4..cabff74`: `reference.py` 184, `catchment_batch.py` 173,
`cli.py` 172, `io/geojson.py` 21, `catchment.py` 13, `mosaic.py` 5,
`io/station_set.py` 4, `burn.py` 0; 90 past the margin and 128 under 700.
With code review round 2's fix (red `44dc954`, green `f952f2d`) it is 573
net (615 added, 42 removed) over `9e666f4..f952f2d`, `cli.py` 173 and the
others as at `cabff74`; 91 past the margin and 127 under 700.
`@developer`'s account of the excess: `StationResult`'s field list, about
50 lines (one line per column, 48 columns, which "The batch" lists
in words rather than counts); the summary's models and its fixed key
lists, about 40 (`Measure`, `Measures`, `Group`, `KnownRefusals`,
`Summary`, and the five lists written in full by the red step's
additions); and `station-catchments`' ten option declarations with their
help, about 45. None of it is a rule the design did not ask for; the
estimate counted the rules and not the declarations. The changes asked
under "PR 4's green step" add a few lines each. **PR 4 is not split**:
the file names no split seam for PR 4, and the ceiling, not the margin,
is what would make one fire. If a later round took it past 700, the seam
is `reference.py` (pure; its own suite, `test_reference.py`) with
`MixedGridError` and `MixedGridRefusal`, about 200 lines, as one PR, and
the batch, the command and the moved writer as the next.

**Lake gauges (Ola, "Fold into PR4") add about 100 lines**: the six rows
marked "Lake gauges" in the table. PR 4 then comes to about 673 net, 27
under 700; with the 44 % margin on the new lines, about 717. PR 4's own
green ran 58 % over its estimate, so passing 700 is likely, not certain. The
ceiling is the rule (CLAUDE.md §2), so **the seam is named now**: if the lake
green step's count over the whole PR passes 700, PR 4 is published as it
stood at `48d1315` (573 net, its review done), and the lake work, whose red
and green commits all come after that commit, becomes **PR 4b, lake gauges**,
on top of it (about 100 to 145 lines, its own review). Ola's "fold into PR
4" is then kept in substance, one branch and one review round of the lake
work before the acceptance run, but not as one pull request; the main
session tells Ola so in the round's recap. Nothing is split before the count
says so. The lake green step (`fff6ac9`) came to 118 net (150 added, 32
removed), 18 % over its 100 and inside the margin, and PR 4 to 691 net
(744 added, 53 removed), 8 lines of room under 700: not split, but a change that
takes it to 700 net or more splits it ("Lake gauges' green step", under
"Ola's rulings").

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
| `reference.py` | agreement, classes, `match_by`, summary | 115 (184 at green) |
| `catchment_batch.py` | `BatchRequest`, `BatchSink`, `run_batch`, `StationResult`, `refusal_cause` | 105 (171 at green; 173 at `cabff74`) |
| `mosaic.py`, `catchment.py` | `MixedGridError`, `MixedGridRefusal` (round 3, Ola's ruling on counting); `check_reach_crs` (change (d)) | 10 (14 at green; 18 at `cabff74`) |
| `cli.py` | `station-catchments`, the directory sink | 80 (140 net at green, the writer's lines removed; 147 at `33b1f2d`; 172 at `cabff74`; 173 at `f952f2d`) |
| `io/station_set.py` | a reference polygon with no area refused (change (b)) | (4 at `33b1f2d`) |
| `io/geojson.py` | the moved writer (moved from PR 3 after its code review, round 1) | 25 (cli.py −25; 21 at green) |
| `fetch/nve.py` | lake layer 5: its allow-list and layer, one envelope-URL function for ELVIS and lakes, the per-station lake query, `lakes.geojson` ("Lake gauges") | 17 (10 at green) |
| `io/station_set.py` | `Lake`, `read_lakes` | 20 (28 at green) |
| `gauge.py` | `LAKE_GAP_M`, `LakeSeed`, `lake_seed` | 18 (27 at green) |
| `catchment_batch.py` | the `lakes` argument, the path per station, five columns | 18 (22 at green) |
| `reference.py` | `by_seed` | 6 (3 at green) |
| `cli.py` | `--lakes`: option, reading, CRS check, the stderr words | 20 (28 at green) |
| **PR 4, the batch and the comparison** | | **about 335 (482); 530 net (567 added) at green `1c39ef7`; 541 net (580 added) at `33b1f2d`; 572 net (614 added) at `cabff74`; 573 net (615 added) at `f952f2d`; with lake gauges about 673 (about 100 more; 717 with the margin on them); 691 net (744 added) at lake green `fff6ac9`** |
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
  Needs PRs 1 to 3. Two findings of PR 2's placement figures
  (`docs/benchmarks/2026-10-05/nve-placement/README.md`) to look at before
  its batch run: (a) **Narsjø (`2.11.0`), a known uncertain case for the
  tier design.** Its nearest line, 83.9 m off, is a lake centreline of an
  unnamed river (`elvid` 002-34-2695) that carries the station's own
  watercourse number, `002.Q1B`, as Nøra's line does at 127.9 m; the tiered
  rule cannot tell two lines with the same number apart, so it picks the
  same line, and the catchment is 0.0012 km² against NVE's 119 km² (marked
  uncertain). Whether the river's name as a tie-break within a tier would
  help is untested. **Resolved for Narsjø by the lake seed** (its point lies
  inside Narsjøen; 98.2 % / 97.5 %, "NVE's lakes"). PR 4 also takes the lake
  gauges ("Lake gauges"), which changes PR 3's merged `fetch/nve.py` (the
  lake layer, `lakes.geojson`) and `io/station_set.py` (`read_lakes`): PR 4
  touches PR 3's code, not its rules, and "Data use" gains the lake layer. (b) **Burns that lower a node by tens of metres**:
  `88.4.0` Lovatn by 41.4 m, `2.284.0` Sælatunga by 22.7 m, `62.10.0`
  Myrkdalsvatn by 20.0 m (`survey.csv`, column `lowered_max_m`); why is not
  looked at.
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
- `test_catchment_batch.py`, on a synthetic tiled DEM (the two-basin
  terrain of `gauge_fixtures`, not `test_cli_catchment.py`'s valley, which
  has no tributary; "PR 4's red step") with a river file and five stations (one matching a
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

**PR 4, lake gauges** (`@tester`, lean, one red commit after `48d1315`, so
the split of "New and changed files" stays mechanical; hand-built polygons
and the existing fixtures, no network):

- `test_gauge.py`, `lake_seed`: a station inside a lake gives `inside`, the
  station as the point and that lake, with a placement on a river line and
  with no placement at all; a station outside every lake, placed on a lake
  line whose `P` lies in a lake 20 m from the station, gives `lake_line`
  with `P` as the point and `distance_m` 20; the same at exactly 30 m gives
  it, at 30.5 m gives `None`; placed on a river line 10 m below a lake
  (station outside it) gives `None`; placed on a lake line whose `P` is in
  no polygon gives `None`; no placement and a lake 14 m away gives `None`
  (Femundsenden, question 9's default); a station inside two overlapping
  polygons gives a `LakeSeed` holding both; no rule reads an area (two
  lakes of very different size, the station inside the smaller, gives the
  smaller).
- `test_station_set.py`, `read_lakes`: a MultiPolygon feature gives one
  `Lake` per part, same number and name; `vatnlnr` 0 and null give number
  None; no `crs` member, a LineString, an empty geometry and a zero-area
  ring are each refused naming the feature's `objectid`.
- `test_fetch_nve.py`: the fake service answers layer 5 (in
  `nve_fixtures.py`); each lake request names exactly
  `objectid,vatnlnr,navn,areal_km2` (never `*`, `globalid` or `kommune`) and
  the station point ± 100 m as its envelope, one request per station; a lake
  seen from two stations is written once; a truncated lake answer is refused
  naming the station; `lakes.geojson` is in the files and the manifest, and
  reads back through `read_lakes`; an output directory holding the four
  older files but no `lakes.geojson` is fetched again without `--refresh`.
- `test_catchment_batch.py`: on `gauge_fixtures.two_basins()`, a lake
  polygon drawn over one basin's lowest part and a station inside it gives
  a row with `seeded_by` `lake`, `lake_rule` `inside`, the gauge, burn and
  sensitivity columns None, `causes` empty, and a catchment equal to
  `delineate` with that lake (22's path) run directly; with a reference
  drawn from that catchment it is `match`, never `uncertain`; the same
  station without `lakes` keeps its river row; `summary.by_seed` counts one
  of each; a station inside two overlapping lakes is `refused` with cause
  `other` and the batch goes on.
- `test_cli_station_catchments.py`: `--lakes` writes the five columns and
  the stderr words "(seeded by the lake <name>)"; a lakes file in another
  CRS (EPSG:32633) is refused naming `--lakes` and both CRSs, before
  `--out-dir` is created; without `--lakes`, every row is a river row and
  the five new columns are empty.

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
2. `rasputin station-catchments` over all 140, with `--lakes`, at the defaults (map radius
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
   `close` rows are summarised by cause. **Lake rows** (`seeded_by` `lake`)
   are reported apart, by class; each with ours-in-NVE's under 95 % gets a
   line saying where the extra area lies (a polygon past the outlet, a
   station on an inflow, a divide; `26.29.0` is expected), and each with
   NVE's-in-ours under 95 % one saying why (`83.2.0` is expected). The river
   rows with a lake on their reach above `P` (24 in "NVE's lakes") are listed
   apart with their classes: that list answers question 10. **Expected refusals**: the two
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
   by cause, not as a failure. So are `191.2.0` and `203.2.0`, which the lake
   seed's probe saw refused on two grids ("NVE's lakes") though neither is
   among the nine.
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
from PR 2 (2026-10-05), 7 from its placement figures (the same day), and 8
from PR 4's green step (the same day), 9 and 10 from the lake-gauge design
(the same day). Ola ruled all six on 2026-10-05: kept as built, which is
each one's default, to be reassessed after the full 140-station run ("Questions
5 to 10: Ola's ruling", under "Ola's rulings").

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
   numbers come out absurd. *Ruled 2026-10-05, closed: the default, to be reassessed after the full run.*
   *Default: refuse the station,
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
   catch. *Ruled 2026-10-05, closed: the default, to be reassessed after the full run.*
   *Default: yes, mark the station uncertain*
   (cause "chain not draining"). Alternative: only stop reading the area
   there, with no other effect, as the first code did; such stations are
   then scored.
7. **Which station the example meshes.** The example under "Closes" meshed
   Narsjø (`2.11.0`), the first station in NVE's list. PR 2's placement
   figures show that Narsjø gets a 0.0012 km² catchment (NVE's is 119 km²):
   its nearest mapped line is a small unnamed river's line across the lake,
   not Nøra's, and the station is marked uncertain. *Ruled 2026-10-05, closed: the default, to be reassessed after the full run.*
   *Default: mesh Atnasjø (`2.32.0`) instead*, a lake gauge whose catchment
   (459.8 km²) agrees with NVE's to 98.4 % one way and 99.2 % the other.
   Alternative: keep Narsjø and fix its placement first (PR 4's tier
   design; see "The PR split"), or pick another station. The lake seed
   ("Lake gauges") now gives Narsjø a catchment that agrees with NVE's to
   98.2 % and 97.5 %, so keeping Narsjø is possible once that lands.
8. **Which area is "our area" in the comparison table.** Two areas of our
   catchment are close to each other but not equal: the number of DEM
   points inside our outline times the area of one cell, and the area of
   the outline itself. The outline cuts the corners between the DEM
   points, so it is half a cell smaller, 50 m² on the 10 m DEM, whatever
   the catchment's size. The table's `fine_area_km2` and the area ratio
   against NVE's polygon use the point count, the same count the two
   overlap percentages use; the catchment file's `fine_area_m2` keeps the
   outline's area, as increment 22 wrote it. *Ruled 2026-10-05, closed: the default, to be reassessed after the full run.*
   *Default:
   keep the point count*, so the area ratio and the overlaps are
   measured on the same thing. Alternative: use the outline's area in the
   table, so that it equals the file's figure; the area ratio then drops
   by half a cell over NVE's area (0.005 % on a 1 km² catchment), and a
   new test pins it.
9. **A lake gauge with no river line near it.** Femundsenden (`311.4.0`)
   has no mapped river line within 500 m, but its point is 14 m from the
   lake Femunden. A station is seeded with a lake when its point lies in
   the lake, or when the river map puts it on the lake's line and it is
   within 30 m of the lake. Femundsenden meets neither, because there is no
   river line to say which way the water runs past it: a station on a river
   flowing into a lake would wrongly get the whole lake's catchment.
   *Ruled 2026-10-05, closed: the default, to be reassessed after the full run.*
   *Default: no; Femundsenden stays refused until the nearest-stream
   fallback (PR 5).* Alternative: seed with a lake within 30 m when there
   is no river line (a few lines; Femundsenden would then get
   Femunden's catchment).
10. **A gauge a little below a lake's outlet.** 24 stations sit on the river
   within 1 km below a lake. They keep the river method: the river is
   followed from the gauge up onto the lake, and the catchment is
   everything that drains through the gauge, the lake included. The four of
   them already run in full agree with NVE's catchments to 97 % or better.
   *Ruled 2026-10-05, closed: the default, to be reassessed after the full run.*
   *Default: keep the river method for them, and let the full run's list of
   these 24 show whether any fails.* Alternative: seed them with the lake
   plus the river from the outlet down to the gauge. That needs the
   catchment code to combine a lake and a river stretch, which it refuses
   today, and a new rule for the position check, about 60 more lines, past
   PR 4's limit, so a later PR.

## Review

(Renumbered 28 to 29 on 2026-10-04: increment 28 is `28-no-fma-contraction.md` on branch `worktree-fpc`. The rounds below were recorded as "28" and are left so.)

**28, design review, round 1, 2026-10-04.** Range `104c883..682c36c` (design and ROADMAP row only). Verdict: CHANGES REQUESTED. LOC: 0 production; estimates PR 1 about 190, PR 2 about 480, both under the ceiling, but 22's `catchment.py` ran 44 % over its estimate and PR 2 names no split seam. Re-measured and true: the HRD PDF (sha256, quotations, 140 unique rows); layers 0/14/38 (fields, EPSG:25833, all 140 stations, 130×1 + 10×3 polygons, 42 points outside, farthest 233 m); areas, size bands and tiles per catchment; no NoData within 250 m of any station; NLOD; HydAPI 401; all seven DOIs; the accumulation oracle (the flood's visit order does not depend on the seed). Blocking: (1) the "500 unregulated" are 364 with regulation 0 plus 136 with none recorded (`:68-73`, Question 1); (2) Jenson 1991 is the nearest-stream-cell snap, and Lindsay et al. 2008 prefer it over the max-accumulation rule chosen here (`:193-214`, Question 4); the departure is unnamed and Ola ruled (a) without that option; (3) nine stations (`156.15.0`, `196.11.0`, `206.3.0`, `208.2.0`, `208.3.0`, `209.4.0`, `212.49.0`, `213.2.0`, `223.2.0`) straddle the eight half-cell-shifted tiles and will be refused by the mosaic's mixed-grid rule, contrary to `:133-135`, acceptance step 4 and the ROADMAP row; (4) the divide offset of a square shifted one cell along an axis is half a cell, not one (`:640`); (5) median vertices 969.5 (not 983), inside-distance quartiles 20/75/355 m (not 20/78/351), median area 131 km² (not 135). Not pushed; no CI.

**28, design review, round 2, 2026-10-04.** Range `26588b1..fe01d37` (design and ROADMAP row only). Verdict: CHANGES REQUESTED. LOC: 0 production; estimates PR 1 about 150, PR 2 about 375 (540 with margin), PR 3 about 300, PR 4 about 295, PR 5 about 80, all under the ceiling. Round 1 blockers 1, 2, 4 and 5 are closed; blocker 3 is closed in this file but not in the ROADMAP row. Re-measured and true: ELVIS layer 2 (fields, EPSG:25833, 1,954,539 features), nearest-line lengths 392/751/1196 m (maximum 4673), distances (median 20 m, maximum 461), 113 own watercourse number, 46 lake lines, 16 and 5 second branches from the station point; the three new DOIs and the Lindsay and WhiteboxTools quotations; the accumulate flags equal `upstream`'s `touches_edge` and `touches_nodata`. Blocking: (1) the strict descent does not make each chain node drain into the next (a pit at a lowered chain end is filled at its spill level and ordered breadth-first; a lower neighbour off the chain can take the flow), so `:626-630`, `:672-673`, `:718-720` and "nesting by construction" `:770-773` overclaim: extend the chain past a lowered end, trust counts by their flags, make `monotone` strict; (2) the loop rule `:615-616` breaks 8-connectivity: cut the loop; (3) a station whose downstream side was not read to U can be scored: make it `uncertain`; (4) ELVIS has 25 `objekttype` spellings including null (station `152.4.0`) and serves exact duplicate segments under different `objectid`s (2 copies at `82.4.0`, 5 at `139.35.0`): a total lake/river mapping, de-duplication by geometry, a rule for forks; (5) `confluence_near` within 100 m of P is 19 stations, not 16 (`:584-585`); (6) PR 3 needs PR 2's `io/rivers.py` (`:1035-1036`); (7) `ROADMAP.md:54` omits the nine shifted-tile refusals and Femundsenden. Suggested: a short "Data use" section with explicit field allow-lists (never `stasjoneier` or `oppdatertav`); a one-sided swing measure, or say why it adds both sides; Soille et al. 2003 carving as prior art. Not pushed; no CI.

**29, design review, round 3, 2026-10-04.** Range `fe01d37..a14e7e6` (design and ROADMAP row only, including the renumbering). Verdict: CHANGES REQUESTED. LOC: 0 production. Estimates add up: PR 1 155; PR 3 380 (547 with the 44 % margin); PR 2 415 (598); PR 4 295 (425); PR 5 80 (115). All are under the 700-line ceiling. Round 2 blockers 1 to 7 and its three suggestions are closed as written. No "28" referring to this increment is left (`git grep` over the tree). Re-measured and true: ELVIS `objekttype` has 25 values summing to 1,954,539, including null (94) and a blank (2, a single space), with eight lake spellings, all starting `innsj` when casefolded. Exact copies: 9 pairs at `82.4.0` and 4 groups of five at `139.35.0`. At `79.3.0`, one `strekninglnr` has two geometries under one `elvid`. `152.4.0` has five null-type lines. On the nine straddling stations, neither lattice's tiles cover the NVE polygon. `156.24.0` and `213.4.0` are covered by the shifted lattice alone. The shift is in x only (world files). Soille et al. 2003 checked against Crossref and the abstract quotation. The float32 spacings. `WINDOW_MARGIN_M` = 2000. The three at-risk citations from `check_citations.py` were re-read and hold. Blocking: (1) A chain with a square corner, or any two non-consecutive chain nodes that are 8-neighbours, breaks `drains` and `monotone` even after a correct burn. The flood names a node's flooder when it is first pushed. So the lower node `chain[k+2]` pushes `chain[k]` before `chain[k+1]` can. Simulated with the flood of `upstream.hpp`: the path (2,4),(3,4),(3,5) gives `flow_to` (2,4)→(3,5) and (3,4)→(4,5); the same path with a diagonal corner drains node by node. Straight lines between consecutive valley-floor nodes, loop cuts and the end extension all make such corners. Make the chain taut before the descent and after the extension: where `chain[k]` is an 8-neighbour of `chain[j]` with `j > k+1`, drop `chain[k+1..j-1]`, as the loop cut does. Add a test that a cornered chain drains. (2) `:1474-1476` sends PR 2 to 23's fetch rules, but `fetch/` is PR 3 since round 2. (3) `:749-750` and `:786`: "at most 0.14 m" no longer holds once the end extension exists. The chain can reach 1000 + 600 + 500 m, which is about 210 nodes, so about 0.21 m. (4) `:239` "254 tiles of 5051 × 5051" is false for 11 tiles. Seven shifted tiles are 5052 × 5053 and `7507_4` is 5052 × 4103. `7305_3` is 5051 × 2881, `7405_1` 3521 × 5051 and `7405_2` 3511 × 5051. This claim predates the range. While fixing it, state the cause of the shift, measured by the main session on 2026-10-04: the data were resampled onto a shifted grid, not mislabelled, so the later fix is a resample. Not pushed; no CI.

**29, design review, round 4, 2026-10-04.** Range `a14e7e6..e7b5923` (design and ROADMAP row only). Verdict: CHANGES REQUESTED. LOC: 0 production. Estimates: PR 1 155; PR 3 380 (545); PR 2 415 (598); PR 4 310 (446); PR 5 80 (115), all under the 700-line ceiling. The 130-line figure script is not counted, by Ola's ruling of 2026-10-04 ("map questions: 1: No, 2: Yes."). Round 3's blockers 1, 2 and 4 and its four suggestions are closed; blocker 3 is closed except for the node count in blocker 1 below. The taut rule was checked with an independent copy of `upstream.hpp`'s flood on 2,000 random chains with repeats, loops and turn-backs: 1,738 raw chains failed `drains`, no taut one did, and every output was taut and 8-connected. On 400 runs of the end extension across an embankment, the step rule kept every chain taut, and no node was flooded from a non-consecutive chain node; the remaining failures all start at a lower node beside the chain, which `drains` reports. Re-measured and true: the eleven tile sizes; the x shift on exactly the eight tiles, with y on the lattice; PixelIsArea and tag/`.tfw` agreement on all 254; the overlap differences of `7707_3` and `7507_4`; the 0.21 m bound; 16c's LOC table without `render.py`; matplotlib in no dependency list; `mosaic.py:385` and `catchment._plan` as described; 15a's tests catch by `MosaicError`, so a subclass leaves them unchanged. Blocking: (1) `:832-838`: "36 nodes if every step is diagonal" and "ends at most 500 m from the reach's last node" disagree, since 36 diagonal steps are 509 m. State the stopping rule (stop before a step that would pass 500 m: 35 diagonal, 50 straight), as the red suite's cap test depends on it. (2) `:864-866`: "the catchment is unchanged" was edited in this range and is contradicted by the new `:853-858` (lowering a depression's spill point moves its drainage). Qualify it: the catchment is unchanged except where a chain node was another depression's spill point. Not pushed; no CI.

**29, design review, round 5, 2026-10-04.** Range `e7b5923..786584e` (design file only). Verdict: APPROVED. LOC: 0 production. Estimates: PR 1 155; PR 3 380 (547); PR 2 415 (598); PR 4 310 (446); PR 5 80 (115), all under the 700-line ceiling. The 130-line figure script is not counted (Ola's ruling of 2026-10-04). Round 4's two blockers and five suggestions are closed. 786584e is complete: status line, round 4 entry, no half-edited section. Re-measured and true: 35 diagonal steps 494.97 m, 36 steps 509.12 m, 50 straight steps exactly 500.0 m; on all 254 DTM10 tiles, spacing 10 m, NoData −32767, EPSG:25833 and PixelIsArea (`tifffile`); `_mixed`'s five causes and the two tiles `_covering` passes it (`mosaic.py:385`, `:571-593`); the neighbour sets of `7707_1` (4) and `7707_3` (7) under the stated rule; `7707_3`'s eleven values, to 0.01 m, including both near ties; `7707_1` declared 0.02 to 0.11 m, moved 0.19 to 0.45 m (the file says 0.44). Required with this record: `ROADMAP.md:54` still says "awaiting round 4". Suggested: 0.45 m without the parenthetical; "10:51 UTC"; the status sentence names `7707_3` and the survey change; the mixed-grid summary line adds "and neither grid covers the window alone". Not pushed; no CI.

**29 PR 1, code review, round 1, 2026-10-04.** Range `51275f5..edef966` (red 0624bbf and 8d87f77, green e43020b, scaffolding removal edef966). Verdict: CHANGES REQUESTED. LOC: 182 added and 53 removed, 129 net. Estimate "about 155" added: 17 % over, inside the 44 % margin, far under 700. Own Release build ctest 932/932; ASan/UBSan pass; full pytest 4069 passed, 17 skipped; gates green. Mutation round: callback 8/8 killed; reverse pass 7/8 (survivor equivalent: order[0] is always the first outlet); starting bits and size check 7/8 (survivor: limit 0xFFFFFFFE, needs a 40 GB raster); flood.hpp's new code 4/4. Blocking: (1) project_structure.md:384 says no accumulation pass and :75-78 lacks flood.hpp and accumulate.hpp; fix in PR 1. (2) Status line (29-nve-reference-catchments.md:3) and ROADMAP.md:54 say "ready for the red step of PR 1". Before push: merge master (120 behind, ROADMAP conflict), re-run check_citations. Suggestions: bindings/core.cpp:1268 say "the node's flooder"; tests/python/test_core_accumulate.py:22 "not on CI runners"; tests/cpp/unit/test_hydrology_accumulate.cpp:14-28 drop restated interface; design :682-685 add on_reach(j, j) and NoData-adjacent outlets. Not pushed; no CI.

**29 PR 1, code review, round 2, 2026-10-04.** Range `edef966..9736bfd` (81392f5 round 1 recorded; 0b86517 project_structure, status line and ROADMAP; 725c8cd design states on_reach(j, j) and the NoData-adjacent outlets; c59ed6b merge of origin/master 99093af; 5dbfec0 keeps project_structure.md:166 in place; f5e8e76 accumulate docstring; 9736bfd two test comments). Whole PR `99093af..9736bfd`. Verdict: CHANGES REQUESTED, on citations alone. LOC: 181 added and 53 removed, 128 net (round 1: 182/129; the only production edit since then, f5e8e76, sits inside a raw docstring, and the one-line difference comes from how the count treats the diff's alignment at the upstream docstring's closing line). Estimate "about 155" added: 17 % over, inside the 44 % margin, far under 700. Fresh Release build ctest 994/994 (test_hydrology_accumulate: 11 cases, 163,676 assertions); extension rebuilt into the worktree venv, full pytest 4420 passed, 17 skipped; mypy, ruff check, ruff format, prohibited-deps, detria boundary and check_citations all exit 0. All round-1 blockers and suggestions closed, and each new claim was checked against flood.hpp, accumulate.hpp and the binding. No red-step scaffolding. The merge left hydrology/ and raster/ untouched, so round 1's mutation record stands. Blocking: the branch adds 8 lines to bindings/core.cpp above line 927 (one include at line 16, the BoundAccumulate struct at about line 316), so three live citations now point 8 lines too high: `15f-edge-strip.md:565` and `:1537` cite `bindings/core.cpp:1174` (the gil_scoped_release is now at 1182), and `27-node-sampling.md:133` cites `bindings/core.cpp:927-932` (the sample docstring is now at 935-940). Re-cite them and re-run check_citations. Not pushed; no CI.

**29 PR 1, code review, round 3, 2026-10-04.** Range `9736bfd..36c2370` (one docs-only commit: the round-2 record, `15f-edge-strip.md:565` and `:1537`, `27-node-sampling.md:133`, the status line and ROADMAP row 29). Whole PR `99093af..36c2370`. Verdict: APPROVED. LOC: 181 added and 53 removed, 128 net, unchanged from round 2 (no production file in this range). Estimate "about 155" added: 17 % over, inside the 44 % margin, far under 700. `git diff --stat 9736bfd..36c2370` lists only ROADMAP.md, 15f-edge-strip.md, 27-node-sampling.md and 29-nve-reference-catchments.md (7 added, 5 removed lines), and nothing under bindings/, include/, src_python/ or tests/. So round 2's build, test, gate and mutation results stand. Round 2's three blockers are closed. Both 15f citations now give `bindings/core.cpp:1182`, which is the `const py::gil_scoped_release unlocked;` in the `refine_points` binding (lines 1173-1183). `27-node-sampling.md:133` now gives `bindings/core.cpp:935-940`, which is the `sample` docstring, "Bilinear z at each of the (N, 2) points" through "never NaN.". `check_citations.py --base origin/master` exits 0. All 34 lines on its at-risk list were re-read as quotations. The live ones hold: `05b-noder-driver.md:1749` (`tests/cpp/CMakeLists.txt:189`), `tests/python/test_features.py@36c2370:583` (`project_structure.md:166`), `29-nve-reference-catchments.md:55` (`ROADMAP.md:54`), and the three just re-cited. The rest sit in dated review records or a retrospective, which are history and left as written. The status line and ROADMAP row 29 match the tree. The branch contains origin/master `99093af`. No red-step scaffolding. No @perf run is needed: PR 1 touches neither refine nor mesh code (design, line 1395). Suggestion: `15f-edge-strip.md:1537` says the number was recomputed when "increment 29 PR 1 took in master". In fact the 8-line shift comes from PR 1's own additions to `bindings/core.cpp` (the include and `BoundAccumulate`), not from a master merge. Not pushed; no CI.

**29 PR 3, code review, round 1, 2026-10-04.** Range `a1445c6..aae91bb` (rulings `b519bd9`, red amendment `05f348a` and `c18c8d5`, green `aae91bb`; the red step `213a2ef` before it). Verdict: CHANGES REQUESTED. LOC: 348 added and 7 removed, 341 net, under the estimate of about 380 (355 once the GeoJSON writer moves to PR 4) and far under 700. Blocking: (1) NVE's river service sends `vatnlnr` = 0 for "no lake" (live query: 956,447 features with 0; both blank-type features, `objectid` 11166506 "Cap'pirjåkka" and 11513676, have 0 and are rivers), but `io/rivers.py`'s `kind_of` treats any non-null value as set, so a blank-type river with 0 is a lake; the design's "set" (`:772`, `:1678` at `aae91bb`) must say not null and not 0, and "set on 496 river features" (`:300-301`) counted the zeros (the reviewer's sample has 2 river features with a positive number); needs a red test and a fix ("PR 3, after code review round 1", under "Ola's rulings"). (2) Status line and `ROADMAP.md:54` still say a red amendment is to come. (3) `15f-edge-strip.md:591` cites `cli.py:1561` and `:1580`, which PR 3 moved to `:1594` and `:1613`. (4) The PR 3 file table lists the GeoJSON writer's move to `io/geojson.py`, which PR 3 neither does nor needs. (5) `project_structure.md` lacks `io/station_set.py`, `io/rivers.py`, `repository.py`'s `read_json`, the `data/` folder, and the `fetch/` and `sources.py` entries missing since 23a-2. Suggestions: a malformed user file refused with `ValueError`, not `TypeError` or `KeyError`; `fetch-stations` reports a missing field in plain words, not a raw `KeyError`; `@tester`'s present-tense "HOW THIS FILE GOES RED" paragraphs (left to `@orchestrator`, not this round). Not pushed; no CI.

**29 PR 3, code review, round 2, 2026-10-04.** Range `aae91bb..a4e4726`; whole PR `d926644..a4e4726`. Verdict: CHANGES REQUESTED. LOC for the whole PR: 387 added, 7 removed, 380 net (tokenize count; green `c4e50ff` alone 36 added, 15 removed, 21 net), against the estimate of about 355 (511 with the 44 % margin): inside the margin, far under 700. Full pytest at HEAD on a rebuilt `_core`: 4499 passed, 118 skipped, 8 failed, all eight `test_settings_wiring`'s executable-hook check, which fails only in the git-archive scratch copy (27/27 in the real worktree). mypy, ruff check, ruff format, prohibited deps, detria boundary and check_citations all clean. Red `b59c648`: 21 tests fail for the intended reasons; all pass at HEAD; green touches no test file. All round-1 blockers and suggestions closed. No mutation testing owed: only `accumulate`'s oracle is invariant-critical, and its kill record is with PR 1. No `@perf` run: PR 3 touches no refine or mesh code. Blocking: (A, `@tester`) `tests/python/test_features.py@a4e4726:583` cites `project_structure.md:166`, now `:190`; (B, `@architect`) the status line and `ROADMAP.md:54` still say the red test and the fix "come next"; (C, `@tester`) the present-tense "HOW THIS FILE GOES RED" paragraphs in `test_fetch_nve.py`, `test_rivers.py` and `test_station_set.py` describe the tree before green and are false now. The PR 3 round-1 record above stays as written. Suggestions, for PR 4 or later: a feature that is a JSON list gives `AttributeError`, and a malformed reference polygon gives shapely's `TypeError`, not a `ValueError` naming the fault; a one-vertex LineString is accepted as a river segment (PR 2's chain may want it refused); `fetch-stations` catches `KeyError` over more than the service answer (narrow the `try`); stale `cli.py` citations in 15e and 25 predate this branch (for `@orchestrator`'s re-citation sweep). Not pushed; no CI.

Fixes for round 2: A and C by `@tester` in `f91cdaa` (re-cited to `project_structure.md:190`; the three paragraphs deleted); B by `@architect` in the commit that records this round (the status line and `ROADMAP.md:54`). The PR 1 records above (round 2, round 3) cite `project_structure.md:166` as the file was then and are left as written.

**29 PR 3, code review, round 3, 2026-10-04.** Range `a4e4726..c7d427a` (`f91cdaa`, `c7d427a`); whole PR `d926644..c7d427a`. Verdict: CHANGES REQUESTED. Delta: 7 lines added, 16 removed, in `ROADMAP.md`, this file and three test files; no production code, so the count stays 380 net and round 2's build, test and gate results stand. Round 2's blockers A, B and C are closed. `check_citations.py --base d926644` exits 0; the at-risk lines this delta could move re-read as quotations and hold. `ruff check` and `ruff format --check` clean. Blocking: the status line and `ROADMAP.md:54` said round 2 asked for "two citations"; it asked for one (`tests/python/test_features.py@c7d427a:583`). Fixed by the main session in the commit that records this round.

**29 PR 3, code review, round 4, 2026-10-04.** Range `c7d427a..5bdcc93` (`c1eeec4`, `5bdcc93`; main session, docs only); whole PR `d926644..5bdcc93`. Verdict: APPROVED. Delta: 5 lines added, 2 removed, in `ROADMAP.md` and this file; production count stays 380 net (estimate about 355, 511 with the margin). Round 3's blocker closed: the status line and `ROADMAP.md:54` say one citation. The round-3 record's range, commits and counts match `git log` and `git diff --shortstat a4e4726..c7d427a`. `check_citations.py --base d926644` exits 0; no at-risk line moved by this delta. No mutation record or `@perf` run owed. CI is owed after the push: `gh pr checks` green before merge.

**29 PR 3, code review, round 5, 2026-10-04.** Range `3dff7bc..ea489e5` (master merge `af2bc38`; `4e877ae` and its revert `09a44c6`; Ola's ruling `9607b98`; `ea489e5` deletes the `NOTICE.md` test). Verdict: CHANGES REQUESTED. Production lines unchanged since round 4. The merge is clean (`git merge-tree --write-tree 3dff7bc 791abb6` gives the tree of `af2bc38`). Full suite on a rebuilt `_core`: 4572 passed, 118 skipped, exit 0, so the prose-read hook stayed silent. `ruff check` and `ruff format --check` clean. The remaining NVE credit tests exist (CSV header, catalogue, every fetch's `NOTICE.txt`); `NOTICE.md` still credits NVE. `check_citations.py` exits 0. Blocking: the ruling said only the 3.14 leg failed, but all three Python legs did; "Data use" said every rule is held by a test, including the `NOTICE.md` credit. Suggestions: name whose §"Licence" and §"Data use" are meant; update the status line. All four fixed by the main session in the commit that records this round.

**29 PR 3, code review, round 6, 2026-10-04.** Range `ea489e5..741504f` (one docs commit, main session; `ROADMAP.md` +1/-1, this file +9/-5). Verdict: APPROVED. Round 5's two blockers and two suggestions closed: the ruling names all three Python legs (`gh pr checks 178` agrees); "Data use" excepts the `NOTICE.md` credit, checked by hand, and the tested credits are tested (`test_fetch_nve.py`: CSV header, catalogue, `NOTICE.txt`); `NOTICE.md`'s "Data in the package" exists; status line and ROADMAP row 29 match the history. `check_citations.py` exits 0; the lines this commit moves are cited only by dated records, except this file's `:55` on `ROADMAP.md:54`, which holds. Master's `accd52a` (h13) touches no file on this branch; `git merge-tree --write-tree HEAD origin/master` is clean, so no merge before the push. CI is owed on the new head; red CI there voids this approval.

**29 PR 2, code review, round 1, 2026-10-05.** Range `529613a..0588bfa` (red `3a06ffc`, amendment `5357785`, green `91857df`, `@architect`'s reading `6f55af3`, second amendment `3c845da`, green `0588bfa`). Verdict: CHANGES REQUESTED. LOC by tokenize over the five `src_python` files: 578 added, 20 removed, 558 net (`burn.py` 167, `gauge.py` 117, `sensitivity.py` 85, `catchment.py` 94, `cli.py` 95 net). This record uses the net figure, as PR 3's did; against the estimate of about 415 (598 with the 44 % margin) it is inside the margin on either basis, and far under 700. Full suite on a rebuilt `_core`: exit 0, 4799 passed, 17 skipped, the prose-read hook silent; mypy, ruff check, ruff format, prohibited deps, detria boundary and check_citations green. Red before green for both pairs: `5357785` 21 failed and 103 errors, then `91857df` 171 passed; `3c845da` on `91857df` 3 failed, then `0588bfa` 172 passed. The suite can fail on both new rules: the whole-disc floor rule fails 8 tests, and dropping the straight-before-diagonal tie fails the embankment test. No mutation testing owed (only `accumulate`'s oracle is invariant-critical, killed with PR 1); no `bench.py` run (no refine or mesh code). `@perf`'s placement figures are still owed before PR 2 is complete. Blocking: (1) a NoData node on the chain: `_floor_node` never picks NoData, but `_join` (`src_python/tin_engine/burn.py@0588bfa:86-97`) links chosen nodes by straight 8-connected runs without looking at the nodes between, and the taut cut can map the placed node onto one. Probe: float32, NoData −32767, floor in column 20 above row 20 and column 22 below, `z[20, 21]` NoData, line down column 20, corridor 30 m: the chain runs through `(20, 21)` and everything below is burnt to about −32767.002 m; with a positive sentinel (3.4e38) the hole is burnt to 290.499 m and `lowered_max_m` of about 3.4e38 reaches the catchment file; the direction check's means read the sentinel too. Write the rule into step 2, covering the taut cut and the direction check, with the refusal as the default for Ola. (2) "Sensitivity", step 2: "its own count is not a sample and its flags do not matter" holds only where `D` lies past `U`; a node exactly at `U` is a sample in `(0, U]` and must be trusted (a flag at +30 m with `U` = 30 gives `downstream_unread`, which the code does correctly); fix the design, and the comment at `src_python/tin_engine/sensitivity.py@0588bfa:71-72` with the green step. (3) The status line and `ROADMAP.md:54` behind the tree; the estimate block's 556; `25-plain-output.md:53`'s citation of the catchment properties. (4) The upstream-flag rule presented as ruled by `@architect`: mark it provisional, an open question for Ola. Suggestions: state the cross-section's bound (half a step; a downstream lean of 3 to 4.5 m on average on a floor flat across, at 30°) rather than "at most 4.6 m either way"; "a cross-section with no node" also happens when every node in it is NoData; `25-plain-output.md`'s field table (lines 89-98, 110-111, 128) cites `cli.py` lines that no longer quote it: pin them to the revision described. Not pushed; no CI.

Fixes for round 1, by `@architect` in the commit that records it: (1) "No NoData on the chain" in step 2 (refuse after the taut cut, before the direction check; question 5, default refuse), the red and green asks in the block "PR 2's code review, round 1", the refusal in the class table and the red suites; (2) the wording of `D` in "Sensitivity", step 2; (3) the status line, `ROADMAP.md:54`, the estimate block (558 net, 578 added; `sensitivity.py` 85). `25-plain-output.md:53` cites `cli.py:1925-1928`, and at `0588bfa` `"fine_area_m2"` is at line 1925 and `"outline_tolerance_m"` at 1928 (`grep -n`), so the citation holds and is left as it is; the round's 1929 did not reproduce. (4) The upstream-flag rule marked provisional, question 6, default yes. Suggestions taken: the cross-section's bound, the all-NoData cross-section, and the field table's `cli.py` citations pinned to `586fbc1`, the revision the inventory describes. The production fixes (the NoData refusal and the `sensitivity.py` comment) come with the red amendment and green step next.

**29 PR 2, code review, round 2, 2026-10-05.** Range `0588bfa..96b3881` (`84b5ede` round 1 recorded, red `6c336f4`, green `96b3881`); the whole PR is `529613a..96b3881`. Verdict: CHANGES REQUESTED. LOC over the five `src_python` files of `529613a..96b3881`: 588 added, 20 removed, 568 net (`burn.py` 174, `gauge.py` 117, `sensitivity.py` 85, `catchment.py` 97, `cli.py` 95 net), counted by round 1's method, written here so it can be rerun: `tokenize` decides which lines are code (no blank lines, no comment-only lines, no docstrings; a multi-line string counts its first line only); the added lines are the `+` ranges of `git diff -U0 529613a 96b3881` hunks, judged as code in the file at `96b3881`, and the removed lines are the `-` ranges, judged in the file at `529613a`. Inside the 598 margin and under 700. Red before green: `6c336f4`, 3 tests fail with DID NOT RAISE; `96b3881`, 132 pass. Full suite on a freshly built `_core`: 4795 passed (7 `test_settings_wiring` failures come from the `git archive` copy the run used; 27 of 27 pass in the worktree); the prose-read hook silent; mypy, ruff check, ruff format, prohibited deps, detria boundary and check_citations green. The gap check covers every reader for finite sentinels. No mutation testing owed, no `bench.py` run (no refine or mesh code); `@perf`'s placement figures still owed. Blocking: (1) NaN cells: `burn_reach`'s mask is `raw != m.nodata`, all true when `nodata` is None (`src_python/tin_engine/burn.py@96b3881:134`); the core counts NaN as NoData whatever the sentinel and `RasterMeta` never holds a NaN sentinel, so a NaN-gapped float DEM is burnt with no refusal, `lowered_max_m` NaN, the sensitivity well posed and NaN in the catchment file. (2) The status line, `ROADMAP.md:54` and the estimate block behind the tree. (3) "Stage B's windows give the same chain as stage A's" overstates: stage B's extension can run further; the same taut chain is what holds. Suggestion: `burn_reach` raises its own `ValueError` subclass and only that becomes a `CatchmentError`, so stage B's `except CatchmentError: pass` (`src_python/tin_engine/catchment.py@96b3881:290`) cannot swallow an ordinary bug. Not pushed; no CI.

Fixes for round 2, by `@architect` in the commit that records it: (1) "No NoData on the chain" in step 2 says a NaN cell is NoData, as in the core, with one mask for the cross-section, the chain check and the extension; the red and green asks in the block "PR 2's code review, round 2"; the red suites name the NaN cases; (2) the status line, `ROADMAP.md:54` and the estimate block (568 net, 588 added; `burn.py` 174, `catchment.py` 97); (3) "the same taut chain", with the reason. The suggestion is taken: `BurnRefusal` in step 2 and in the green ask, pinned by `@tester` in the red step. The LOC count above was rerun by `@architect` with a script written from the rule (not committed) and gave the same figures, file by file.

**29 PR 2, code review, round 3, 2026-10-05.** Range `96b3881..193079d` (`cc52b8f` round 2 recorded, red `9625478`, green `193079d`); the whole PR is `529613a..193079d`. Verdict: APPROVED. LOC by round 2's rule over the five `src_python` files: 589 added, 20 removed, 569 net (`burn.py` 175, `gauge.py` 117, `sensitivity.py` 85, `catchment.py` 97, `cli.py` 95 net), one net line more than round 2; inside the 598 margin and under 700. Round 2's blockers are closed: (1) the mask at `src_python/tin_engine/burn.py@193079d:140` excludes NaN whatever the sentinel, matching the core's `is_nodata`; (2) the status line, `ROADMAP.md:54` and the estimate block matched the tree as round 2 left it; (3) the design says "the same taut chain". The suggestion is taken: `BurnRefusal` is raised by all three of `burn_reach`'s refusals, and `_burnt_flood` catches only it. Red `9625478` fails 8 tests for the intended reasons; all pass at `193079d`. Full suite on a freshly built `_core`: exit 0, 4783 passed, 18 skipped (`test_settings_wiring` 27 of 27 run in the worktree); the prose-read hook silent; mypy, ruff check, ruff format, prohibited deps, detria boundary and check_citations green. No other code in the PR builds its own DEM mask: `gauge.py` reads no array, `sensitivity.py` reads `accumulate`'s output, `catchment.py` uses only the array's `.shape`. A probe with NaN cells and a finite sentinel is refused at `193079d`. No mutation testing owed, no `bench.py` run (no refine or mesh code); `@perf`'s placement figures are still owed before PR 2 is complete. The status line lagging the review does not block; it is set in the commit that records this round. Suggestions, not blocking: a fourth `GAPS` case, `(math.nan, -32767.0)`, in `test_burn.py`, to pin NaN beside a sentinel; `burn_reach`'s docstring (`src_python/tin_engine/burn.py@193079d:135-137`) names two of its three refusals (not "no DEM node with data near the reach"). Not pushed; no CI.

Recorded by `@architect`, with the status line, `ROADMAP.md:54` and the estimate block (569 net, 589 added; `burn.py` 175). The LOC count was rerun from round 2's rule and gave the same figures, file by file. Both suggestions are left for a later step, to avoid another review loop for a test case and a docstring: whoever next touches `burn.py` or `test_burn.py` (PR 4 at the latest) adds the `(math.nan, -32767.0)` case and names the third refusal in the docstring.

**29 PR 2, evidence review, round 1, 2026-10-05.** Commit `cdc0805` (`@perf`'s placement figures; range `e0ccdd9..cdc0805`): 25 files, all under `docs/benchmarks/2026-10-05/nve-placement/`, no production lines, nothing outside that directory, no data committed. Verdict: CHANGES REQUESTED. Re-measured and true: the checksums in `provenance.txt` recomputed and matching; the survey recounted (140 stations, 136 burnt, 4 refused; 5 well posed on the first window; lowered nodes fewest 3, median 81, most 187, none zero; 19 stations with another river line within 100 m under the nearest placement, 26 under the tiered one; 48 lake lines); cases 1 and 3 correctly empty, case 2 `101.1.0`, case 4 `2.32.0`; the full runs, 9 well defined (at least 87.9 % both ways) and 9 uncertain (5 under 1 km², 4 at about 98.5 %); Narsjø re-run end to end, matching `cli_confirm.log`. Findings: (1) the README's "Counts that differ" says six stations were refused on two grids where "the design expects nine", but two of the six (`156.24.0`, `213.4.0`) are not among the nine: this file predicted them as findings to explain ("Coverage by DTM10", and acceptance step 4 names both), so the prediction came true and up to 11 can be refused on two grids; Ola's ruling covers "the nine"; (2) the example command under "Closes" meshes Narsjø, which this evidence gives a 12-node, 0.0012 km² polygon (its nearest line is an unnamed lake centreline, `elvid` 002-34-2695, sharing watercourse number `002.Q1B` with Nøra's): switch to Atnasjø `2.32.0` as the default of a new question for Ola, and keep Narsjø as a known uncertain case for PR 4's tier design (the tiered rule cannot tell lines sharing a number apart; whether the river name as a tie-break would help is untested); (3) "46 lake lines" (in "What the data says" and in "Placing the gauge", "46 of 139") needs a note that PR 2's code on the 2026-10-05 data finds 48 under either placement, pointing to the README; (4) the status line and `ROADMAP.md:54` do not say the figures are in. Suggestion: record `88.4.0` Lovatn (a burn lowering a node by 41.4 m; also 22.7 m and 20.0 m elsewhere) as an item to look at before PR 4's batch run. Not pushed; no CI.

Fixes for evidence round 1, by `@architect` in the commit that records it: (1) the README's "Counts that differ" says the prediction for `156.24.0` and `213.4.0` came true, up to 11 on two grids, and that Ola's ruling covers the nine; (2) the example meshes `2.32.0` (question 7, default Atnasjø), and Narsjø is item (a) under PR 4 in "The PR split"; (3) both "46" places carry the note; (4) the status line and `ROADMAP.md:54`. The suggestion is taken as item (b) under PR 4 in "The PR split".

**29 PR 2, evidence review, round 2, 2026-10-05 (copied from `@reviewer`'s handback).** Range `cdc0805..1d58e21` (one docs commit, 3 files, +57/-11, 0 production lines; PR 2 stays at 569). Verdict: APPROVED. Round 1's four findings closed: four of the six two-grid refusals are among the nine, `156.24.0` and `213.4.0` were predicted (9 + 2 = 11); the example meshes `2.32.0` (459.8478 km², well posed, no causes; Narsjø 0.0012 km² against 119.43 km², nearest line `objectid` 12252329, `elvid` 002-34-2695, `002.Q1B`, Nøra's `elvid` 002-34-1299 at 127.9 m with the same number); 48 lake stations in both surveys; status line and `ROADMAP.md:54` current. The PR 4 burn depths match `survey.csv`. `check_citations.py` exits 0. Suggestions: this file's line about `156.24.0` and `213.4.0` "on shifted tiles only … and run" could point to the README; with the tiered placement `22.22.0` Søgne (21.95 m) replaces Myrkdalsvatn among the three deepest burns.

**29 PR 4, code review, round 1, 2026-10-05 (from the main session's summary of `@reviewer`'s handback).** Range `9e666f4..33b1f2d` (red `2b9b39f`, `@architect` `6475806`, red amendment `7cd56a4`, green `1c39ef7`, `@architect` `33c2dc4`, red `6ebff1b`, green `33b1f2d`). Verdict: CHANGES REQUESTED. LOC by PR 2's round-2 rule: 580 added, 39 removed, 541 net (`reference.py` 184, `catchment_batch.py` 171, `cli.py` 147, `io/geojson.py` 21, `catchment.py` 9, `mosaic.py` 5, `io/station_set.py` 4, `burn.py` 0); 59 past the 482 margin, under 700, so no split (the fallback seam applies only past 700). The reviewer's count reproduces 530 at `1c39ef7` and 569 for PR 2. Gates clean. Full suite on a freshly built `_core`: 4900 passed; 8 failures are artefacts of the scratch copy and pass in the worktree. Red before green for both pairs (`2b9b39f`/`7cd56a4` then `1c39ef7`; `6ebff1b` then `33b1f2d`). No mutation testing owed and no `bench.py` run (the 140-station acceptance run comes after PR 4). Prose claims checked: the half cell between `ours` × cell area and the traced outline's area, exact on 8 traced sets; classify and summarise as designed; a refusal against a bug (only `CatchmentError` becomes a row); `MixedGridError` and `MixedGridRefusal`, and `dem_input.py:248`'s re-raise; catchment files byte-identical across `9e666f4..33b1f2d` for three runs, and `results.csv` identical but for `seconds`. `33c2dc4`'s re-citations hold. PR 2's two leftovers are closed. PR 4 does not depend on questions 5 to 7. Blocking: `project_structure.md` is stale: line 203 says "the catchment GeoJSON writer is in cli.py (22)", the paragraph at line 515 records that exception as still open, and there are no entries for `reference.py`, `catchment_batch.py` and `io/geojson.py`. That file is at the repository root, outside `@architect`'s write limit; it waits on Ola's ruling on who edits root files (default: yes, `@architect` commits its drafted text). Suggestions: (1) `station-catchments` should refuse a river file (and with it the references) in a CRS other than the DEM's up front, naming `--rivers`, before creating `--out-dir`; today every placed station becomes an `other` refusal ("the river reach must be in the DEM's CRS, EPSG:25833"), while `catchment --rivers` refuses at once; (2) the `except OSError` around `run_batch` also wraps the sink's writes, so a write failure reads "Invalid value for --dem: cannot read .../out/1.140.0.geojson: [Errno 21] Is a directory"; a write failure should name the output, in plain words; (3) the new suites' `importlib` fixtures (`cb`, `gj`, `ref`) can become plain imports (`@tester`, test-only); (4) a station with no name prints a stray space in its stderr line. Not pushed; no CI.

Rulings on round 1, by `@architect` in the commit that records it: the blocker stays open until Ola rules on root-file ownership; suggestions (1), (2) and (4) are adopted as changes (d), (e) and (f), and (3) as a test-only commit, under "PR 4's code review, round 1"; the status line, `ROADMAP.md:54`, the estimate block and PR 4's table rows are set to `33b1f2d`.

**29 PR 4, code review, round 2, 2026-10-05 (from the main session's summary of `@reviewer`'s handback).** Range `33b1f2d..72b62ce` (test-only `e8e93d0`, red `bc78be3`, green `cabff74`, `@architect` `72b62ce`). Verdict: CHANGES REQUESTED. LOC reproduced: 614 added, 42 removed, 572 net. Red before green: `bc78be3`'s 8 tests fail at `33b1f2d` for the right reasons and pass at `cabff74`. `e8e93d0` changes no assertion (test counts 92, 87 and 14 equal before and after). Full suite on a freshly built `_core`: 4909 passed; 7 failures are artefacts of the scratch copy. Gates clean. Every re-citation in `72b62ce` holds. The round-1 record above, written from the main session's summary, matches what the reviewer would have written. Blocking: (1) "today" clauses in three test docstrings, written at the red step, were made false by the green step; (2) change (e) did not cover `target.mkdir(exist_ok=True)`: an `--out-dir` that is an existing file ends in a `FileExistsError` traceback, and one under a read-only parent in a `PermissionError` traceback, which makes false the design's "every write the command makes" and `ROADMAP.md:54`'s "a failed write reported as a write to the output directory"; (3) round 1's blocker, `project_structure.md`, is still open, waiting on Ola. Not blocking, a follow-up item: the `catchment` command's own write of its catchment file (`cli.py:1956`, `target.write_bytes(...)`) is not under any `except` either; it is older than PR 4. Suggestion, adopted: the ruling for (d) said "Today `station-catchments` checks only the reference file…", now "Before (d), …". Not pushed; no CI.

Fixes for round 2: (1) `@tester`'s `a98c805` rewords the three docstrings to "before change (d)" / "before change (e)", docstrings only; (2) red `44dc954` (`@tester`) adds the cases `out_dir_is_a_file` and `out_dir_parent_read_only` (the latter skipped when run as root) to change (e)'s parametrised test, and green `f952f2d` (`@developer`) makes `--out-dir` under the same `_writing` helper, and also moves the `results.csv` header write and its close under it, which no test exercises; full suite 4919 passed. Recorded by `@architect` in the commit that records this round, with `mkdir` added to (e)'s list of writes under "PR 4's code review, round 1", the reword of (d), the status line, `ROADMAP.md:54`, the estimate block and PR 4's table rows (573 net, 615 added, 42 removed; `cli.py` 173; 91 past the margin, 127 under 700). Blocker (3) stays open.

**29 PR 4, code review, round 3, 2026-10-05 (copied from `@reviewer`'s handback).** Range `72b62ce..07793ab` (`a98c805`, red `44dc954`, green `f952f2d`, `07793ab`); whole PR `9e666f4..07793ab`. Verdict: APPROVED, on one condition: round 1's open blocker, `project_structure.md`, which waits on Ola's ruling on who edits root files. Nothing else blocks. LOC by PR 2's round-2 rule: 615 added, 42 removed, 573 net (reference 184, catchment_batch 173, cli 173, io/geojson 21, catchment 13, mosaic 5, io/station_set 4, burn 0); 91 past the 482 margin, 127 under 700. Round 2's blockers closed: `a98c805` changes docstrings only; making `--out-dir` is under `_writing` (`44dc954` fails two cases, `f952f2d` passes all five; the read-only-parent case ran on macOS). Every filesystem call in `station_catchments` ends in a plain refusal; the untested header write and close were probed with a failing `Path.open` wrapper (exit 2, refusal naming `--out-dir`, no traceback) and a no-failure control. Full suite on a rebuilt `_core` exit 0, 4919 passed, 17 skipped; prose hook silent; gates clean. The round-2 record matches. No mutation testing or `bench.py` owed. Suggestion: a failing-close test case when `test_cli_station_catchments.py` is next touched.

**29 PR 4, code review, round 4, 2026-10-05 (copied from `@reviewer`'s handback).** Range `07793ab..c9a40c7`; whole PR `9e666f4..c9a40c7`. The handback follows word for word, by Ola's ruling that review verdicts are recorded so; only its level-2 headings are lowered to level 4, so that they do not end this section.

#### Result

**Verdict: CHANGES REQUESTED.** Both blocking issues are prose fixes in Markdown, and both are small. The code, the tests, the size count and the gates all pass. No production line has to change, so the count stays at 691 net.

**Size Metrics**
- This round: `07793ab..c9a40c7`. Whole PR: `9e666f4..c9a40c7`. I counted with PR 2's round-2 method: `tokenize` decides which lines are code, added lines are the `+` ranges of the `git diff -U0` hunks judged at the new revision, and removed lines are the `-` ranges judged at the old one. The counting script is in my scratchpad and is not committed.
- **Whole PR: 744 added, 53 removed, 691 net.** This matches the record exactly. Per file (net): `cli.py` 201, `catchment_batch.py` 195, `reference.py` 187, `io/station_set.py` 32, `gauge.py` 27, `io/geojson.py` 21, `catchment.py` 13, `fetch/nve.py` 10, `mosaic.py` 5, `burn.py` 0.
- **Lake work (`07793ab..c9a40c7`): 150 added, 32 removed, 118 net.** Per file: `cli.py` 28, `io/station_set.py` 28, `gauge.py` 27, `catchment_batch.py` 22, `fetch/nve.py` 10, `reference.py` 3.
- **Against the estimate:** the six lake rows add up to 99, so the lake work came in 18 % over, inside the 44 % margin. The split point the design names ("PR 4b, lake gauges") does not apply, because 691 is under 700. Only 8 net lines of room remain.
- **Note:** 744 *added* is over 700. The project has used the net figure since PR 2 and PR 3, and I followed that. Ola may want to know how close this is.
- No new code packed under `# fmt: skip` or `# fmt: off`. No C++ in the PR.
- Focus of the round: the lake path (`gauge.lake_seed`, `read_lakes`, the lake fetch, `run_batch(lakes=)`, `by_seed`, `--lakes`) and `project_structure.md`.

**CI:** the branch `worktree-29-pr4` is not pushed, so there is no PR and no CI yet. This is the review before the first push. CI still has to go green after the push. `git merge-tree` against master `e3203f1` merges cleanly.

**Blocking Issues**
1. **Duplicated text in `docs/increments/29-nve-reference-catchments.md`, introduced by `efb3894`.** Four passages appear twice in a row:
   - lines 22–26: "On 2026-10-05 Ola ruled that lake gauges get increment 22's Bygdin method …"
   - **lines 39–40: `--lakes ../rasputin_data/nve_hrd/lakes.geojson \` appears twice in the example command under "Closes"**, the command the acceptance run copies.
   - lines 56–58: "A gauge on a lake is seeded with the whole lake …"
   - lines 68–71: "Seeding a lake together with the river from its outlet …"
   
   Delete the second copy of each.
2. **`ROADMAP.md:54` (row 29) says something that is now false.** `ef403bd` edited this row, and it still says "four questions open for Ola, with defaults (…)" and lists only questions 5 to 8. The increment file's status line says questions 5 to 10 are open. Questions 9 and 10, from the lake design, need adding: Femundsenden stays refused, and gauges a little below a lake outlet keep the river method.

**What holds (checked against the code)**
- **The lake path matches the design:**
  - `LAKE_GAP_M = 30.0`.
  - `LakeSeed` is frozen and holds `rule`, `point`, `lakes` and `distance_m`.
  - The `inside` rule uses strict `contains`, with or without a placement.
  - The `lake_line` rule needs the placement on a lake line, finds the lake containing the mapped position by geometry rather than by lake number, and allows `<= gap + 1e-6`.
  - A non-finite station gives no seed.
  - `run_batch`: a lake seed gives 22's request (`lakes`, `lakes_crs`, no reach). The `no_river` refusal applies only when there is neither a placement nor a seed. A lake row gets no gauge columns, and `classify` receives `gauge=None`.
  - The five columns sit after `reach_fork`.
  - `by_seed` has the groups `river` and `lake`.
  - `--lakes` uses a shared CRS-check helper (`_read_beside`) that also replaces the `--reference` block.
- **The lake fetch:** layer 5 asks for exactly the four allowed fields, uses a ±100 m envelope through one function shared with the river query, keeps each lake once by `objectid`, refuses a truncated answer, and writes `lakes.geojson` into `FILES`. Because the file is in `FILES`, an old fetch directory without it is fetched again.
- **`project_structure.md`:**
  - The new entries for `gauge.py`, `burn.py`, `sensitivity.py`, `reference.py`, `catchment_batch.py`, `io/geojson.py`, `read_lakes`, the `catchment.py` pour point and `nve.py` (lakes) each match the code's imports and behaviour.
  - The new `io/` exception for `results.csv` and `summary.json` matches the `csv`/`json` writes in `cli.py`.
  - The replacement paragraph is true: both commands call `catchment_geojson`.
- **Leftover red-step comments:** none. The new tests use the agreed "Before the change" wording. The "PR 4's red step" mentions are section citations.
- **Test coverage:** every item in the design's lake test list has a test, including the 30 m / 30.5 m boundary, two overlapping lakes refused with cause `other`, a lake catchment equal to 22's path run directly, the other-CRS refusal before `--out-dir` is made, and the refetch when `lakes.geojson` is missing.
- **Gates run locally, all green:** mypy, ruff check, ruff format --check, prohibited deps, detria boundary, and `check_citations` (it resolves every citation and lists 94 at risk).
- **At-risk citations re-read as quotations:**
  - The live `tests/python/test_features.py@c9a40c7:583` → `project_structure.md:208` still quotes "never imports _core".
  - `29…md:75` → `ROADMAP.md:54` is still row 29.
  - Neither the docs nor the tests have unpinned line citations into this round's changed modules.
  - The rest of the at-risk list are historical review records.
- **Python suite:** full run gave 4989 passed and 17 skipped. I used the worktree's existing venv and `_core` (built Oct 5 08:14; no C++ changed since `c59ed6b`), with bytecode, cache and coverage output sent outside the tree.
- **Mutation testing:** none owed. Only `accumulate`'s oracle is named invariant-critical, and its kill record is with PR 1.
- **`@perf` run:** none owed, because no refine or mesh code is touched.

**Suggestions (non-blocking)**
- `read_lakes`: a ring with too few points, or a `vatnlnr` like `"x"`, is refused with shapely's or `int()`'s bare message, which does not name the lake. The message does reach the user as a `--lakes` error. A string `"0"` gives lake number 0 rather than None. NVE serves integers, so this only matters for files users make themselves.
- `gauge.py` is described as "pure shapely" but now imports `io/station_set.py` for `Lake`. The design chose this. If a later PR touches it, `Lake` could move next to `Reach` so the pure module does not depend on a reader module.

#### Pinned or assumed beyond the design
- The size ceiling is judged on net lines, as in PR 2's and PR 3's reviews.
- The suite ran on the existing `_core` without a rebuild, because I may not build C++ in this run. PR 4 changes no C++ and the `.so` is newer than the last C++ commit.

#### Questions for Ola
- Should the 700-line ceiling count lines added rather than net? PR 4 is 691 net but 744 added. Default: keep counting net, as every review so far has.

#### Lessons
- A design commit (`efb3894`) doubled four passages, probably an edit applied twice, and three later rounds of red, green and record commits did not notice. A check for identical adjacent lines in `docs/increments/*.md` would catch this mechanically. I am passing it on for `@orchestrator`.

#### ASK OLA lines
none

#### GUARD FALSE POSITIVE lines
none

Files: `/Users/skavhaug/projects/rasputin/.claude/worktrees/29-pr4/docs/increments/29-nve-reference-catchments.md` (lines 22–26, 39–40, 56–58, 68–71) and `/Users/skavhaug/projects/rasputin/.claude/worktrees/29-pr4/ROADMAP.md` (line 54). My note file is `/Users/skavhaug/projects/rasputin/.claude/current-task/reviewer-083758.md`; deleting it is the spawner's job.

Fixes for round 4, by `@architect` in the commit that records it: (1) the second copy of each of the four doubled passages deleted (lines 22–26, 39–40, 56–58 and 68–71 as they were), and the whole file scanned for any other repeated adjacent passage, none found; (2) `ROADMAP.md:54`, this file's status line and "Questions for Ola" record Ola's ruling of 2026-10-05 on questions 1 to 10: all kept as built, to be reassessed after the full 140-station run ("Questions 5 to 10: Ola's ruling", under "Ola's rulings"). Both suggestions are deferred to a later PR, since PR 4 has 8 production lines of room left ("Deferred from PR 4's code review, round 4"). The reviewer's question on counting added or net lines is the main session's to put to Ola and is not answered here. No production line changed.

**29 PR 4, code review, round 5, 2026-10-05 (copied from `@reviewer`'s handback).** Range `c9a40c7..bde3b8e`; whole PR `9e666f4..bde3b8e`. The handback follows word for word, by Ola's ruling that review verdicts are recorded so; only its level-2 headings are lowered to level 4, so that they do not end this section. Lines already at level 4 inside its code block are as the handback has them.

#### Result

**Verdict: CHANGES REQUESTED.** There is one blocking issue, and it is in the prose. The record of round 4 says it copies my handback word for word, with only the headings changed. That is not true: the last seven lines of the handback are missing. Everything else in this round passes, and no production line changed.

**Size Metrics**
- **Range reviewed:** `c9a40c7..bde3b8e`, one `@architect` commit.
- **Files changed:** two, `ROADMAP.md` and `docs/increments/29-nve-reference-catchments.md`, both prose. `git diff --name-only` lists nothing under `src_python/`, `tools/`, `tests/` or C++.
- **PR size:** unchanged at 744 added, 53 removed, **691 net**. I judged the 700-line ceiling on the net figure, as rounds 1 to 4 and PR 2's and PR 3's reviews did. Whether the ceiling should count added lines instead is still a question for Ola. Under the net count there are 8 lines of room left.
- **Focus of the round:** the round-4 record, the doubled passages, how Ola's rulings are recorded, and citations.

**CI:** the branch `worktree-29-pr4` has not been pushed, so there is no PR and no CI yet. This is the review before the first push, and CI still has to go green after the push.

**Blocking Issues**
1. **The round-4 record is not word for word, but says it is.** The sentence is at `docs/increments/29-nve-reference-catchments.md:3140`: "The handback follows word for word … only its level-2 headings are lowered to level 4". I pulled my round-4 handback from the run's transcript and compared it with the recorded block, which runs from "#### Result" to just before "Fixes for round 4". After lowering `##` to `####`, the two match up to "#### Lessons". The handback's last seven lines are missing:
   ```
   #### ASK OLA lines
   none

   #### GUARD FALSE POSITIVE lines
   none

   Files: `/Users/.../29-pr4/docs/increments/29-nve-reference-catchments.md` (lines 22–26, 39–40, 56–58, 68–71) and `/Users/.../29-pr4/ROADMAP.md` (line 54). My note file is `.../reviewer-083758.md`; deleting it is the spawner's job.
   ```
   The fix is one of two:
   - **(a)** add those lines after the Lessons block, keeping the `####` level (`tools/brief.py` splits sections only at `## `, so this is safe); or
   - **(b)** change the sentence at line 3140 to say the closing "none" sections and the file list were left out.

   Ola's ruling asks for a verbatim record, so (a) is the better choice.

**What holds (checked)**
- **Production code did not change.** The diff touches only the two Markdown files.
- **The doubled passages are gone.** All four were deleted, including the doubled `--lakes` line in the example command under "Closes". I scanned the whole file three ways: identical neighbouring lines, repeated pairs of long lines, and repeated sentences longer than 60 characters. Nothing doubled is left:
  - The only repeated pair (lines 1140 and 1281) is the same file name in the output of two different `git grep` commands.
  - The two sentences that appear three times are on purpose: the ruling marker "*Ruled 2026-10-04, closed: the default*", and PR 1's estimate, repeated in each of its three review records.
- **Ola's ruling is recorded correctly in all three places.** The status line, `ROADMAP.md:54` (row 29) and "Questions for Ola" all say questions 1 to 10 are ruled, kept as built, and looked at again after the 140-station run.
  - The new paragraph "Questions 5 to 10: Ola's ruling" sits under "## Ola's rulings (2026-10-04)", where the cross-references say it is.
  - Each of questions 5 to 10 now carries "*Ruled 2026-10-05, closed: the default, to be reassessed after the full run.*".
  - The opening paragraph of "Questions for Ola" no longer calls them open.
- **Questions 9 and 10 mean what the lake design says:**
  - 9: Femundsenden (no river line within 500 m, 14 m from Femunden) stays refused until the nearest-stream fallback, PR 5. This matches the lake-seed rule (inside the lake, or on the lake's line within 30 m) and line 2839.
  - 10: a gauge a little below a lake's outlet keeps the river method, and seeding a lake together with the river below it is not built. This matches lines 2027–2033.
  - Both ROADMAP and the ruling paragraph state them that way.
- **The deferred-suggestions paragraph** gives my two round-4 suggestions accurately: `read_lakes` messages that do not name the lake, `"0"` read as lake 0, and `Lake` moving out of `io/station_set.py`. "8 production lines of room" is right: 691 net, and the PR must stay under 700.
- **Citations:** `python3 tools/check_citations.py` exits 0 and lists 99 citations to re-read; I re-read the live ones.
  - `ROADMAP.md:54` is still row 29, so the citation at line 70 holds.
  - `tests/python/test_features.py@bde3b8e:583` → `project_structure.md:208` is unchanged in this range.
  - Nothing cites line numbers in the increment file from outside it, so the 5 lines removed near its top move no citation.
  - The rest of the list is dated review records, which stay as written.
- **Leftover red-step comments:** none possible, since no code changed.
- **Mutation record:** none owed in this round; the one invariant-critical suite's record is with PR 1.
- **`@perf` run:** none owed, because no refine or mesh code is touched.

**Suggestions (non-blocking)**
- The ROADMAP row says "Ola ruled all ten of the increment's questions on 2026-10-05". Questions 1 to 4 were ruled on 2026-10-04, and the 2026-10-05 ruling confirmed them. The increment file has it right ("1 to 4 as recorded above, and 5 to 10 kept as built"). The row could say the same.

#### Pinned or assumed beyond the design
- The size ceiling is judged on net lines, as in every earlier review.
- I took my round-4 handback text from this session's transcript of the round-4 reviewer run. It is the `SubagentHandback` call in subagent `aa531c608e07d23ea`.

#### Questions for Ola
- Should the 700-line ceiling count lines added rather than net (added minus removed)? PR 4 is 691 net but 744 added. Default: keep counting net, as every review so far has. This was carried over from round 4 and is still unanswered.

#### Lessons
- When a handback is copied "word for word", the copy can be checked mechanically. Before writing that claim, run `diff` between the handback text (taken from the transcript) and the recorded block. A hand copy dropped the closing lines here.
- `tools/brief.py` counts recorded review rounds as lines that contain "APPROVED" or "CHANGES REQUESTED". A verbatim handback that uses either word on more than one line makes the count too high. Today the count is still right only because my round-4 text used the verdict word on one line. Passing both lessons on for `@orchestrator`.

#### ASK OLA lines
none

#### GUARD FALSE POSITIVE lines
none

Files: `/Users/skavhaug/projects/rasputin/.claude/worktrees/29-pr4/docs/increments/29-nve-reference-catchments.md` (line 3140, and the record block that ends before "Fixes for round 4") and `/Users/skavhaug/projects/rasputin/.claude/worktrees/29-pr4/ROADMAP.md` (line 54). My note file is `/Users/skavhaug/projects/rasputin/.claude/current-task/reviewer-085322.md`; deleting it is the spawner's job.

Fixes for round 5, by `@architect` in the commit that records it: (1) option (a): the round-4 record now ends with the handback's closing lines (the two "none" sections and the file list), at level 4; both records were checked by `diff` against the handback text, with their headings raised back to level 2, and both diffs are empty; (2) the suggestion taken: `ROADMAP.md:54` says questions 1 to 4 were ruled on 2026-10-04 and confirmed on 2026-10-05, and 5 to 10 were ruled on 2026-10-05. The question on counting added or net lines, asked in rounds 4 and 5, is closed: Ola ruled net (2026-10-05: "net, and write it into CLAUDE.md"), and `CLAUDE.md` §2 now states the ceiling as net production lines, added minus removed, counted by the method these reviews use. PR 4's 691 net therefore leaves 8 lines of room under the ceiling. No production line changed.

**29 PR 4, code review, round 6, 2026-10-05 (copied from `@reviewer`'s handback).** Range `bde3b8e..9bb1723`; whole PR `9e666f4..9bb1723`. The handback follows word for word, by Ola's ruling that review verdicts are recorded so; only its level-2 headings are lowered to level 4, so that they do not end this section.

#### Result

**Verdict: APPROVED.** Nothing blocks. Both commits change prose only. The two review records match my handbacks word for word, and the new wording in `CLAUDE.md` §2 states the counting method these reviews have used. Measured that way, PR 4 is still 691 net lines.

**Size Metrics**
- **Range reviewed:** `bde3b8e..9bb1723`, two `@architect` commits (`3fab6c6`, `9bb1723`). The whole PR is `9e666f4..9bb1723`.
- **Files changed in this round:** `CLAUDE.md`, `ROADMAP.md` and `docs/increments/29-nve-reference-catchments.md`. Nothing under `src_python/`, `tests/`, `tools/`, C++ or `.github/` changed since round 4 (`git diff --stat c9a40c7..9bb1723` over those paths prints nothing).
- **PR size:** 744 lines added, 53 removed, **691 net**. I measured this again with my own counter, in the scratchpad and since removed. It takes the `+` and `-` ranges of `git diff -U0 9e666f4 9bb1723 -- src_python`. It drops blank lines, comment-only lines, docstrings and the bodies of raw strings. Added lines are judged at the new revision and removed lines at the old one. By file: cli.py +201, catchment_batch.py +195, reference.py +187, station_set.py +32, gauge.py +27, geojson.py +21, catchment.py +13, nve.py +10, mosaic.py +5, burn.py 0. The total matches rounds 4 and 5. Since the PR has to stay under 700, there are 8 lines of room.
- **Focus of the round:** the verbatim records, the new ceiling wording, ROADMAP row 29 and the status line, and citations.

**CI:** the branch `worktree-29-pr4` is not pushed (`git ls-remote` and `gh pr list --head` both print nothing). So there is no CI yet. This is the review before the first push, and CI still has to go green after it.

**Checks**
1. **The records are word for word.** I took both handbacks from the `SubagentHandback` calls in the transcripts: round 4 from subagent `aa531c608e07d23ea`, round 5 from `aadbad41b7cbf3dc0`. For each round I took the recorded block from "#### Result" up to the "Fixes for round N" paragraph and raised unindented `#### ` headings back to `## `. Then I ran `diff` against the handback. Both diffs are empty (round 4: 74 lines on each side; round 5: 73 lines on each side). Each handback has 6 level-2 headings, so the heading change really was tested. The level-4 lines inside round 5's indented code block stayed as the handback has them, as the record's header says. Ola's quotation in the round-5 fix paragraph, "net, and write it into CLAUDE.md", matches Ola's message in this session's transcript exactly.
2. **`CLAUDE.md` §2.**
   - *Does it state the rule the reviews used?* Yes. It defines net lines as added minus removed, with the same exclusions on both sides. It takes the lines from the PR's `git diff -U0` hunks against its base, judging added lines at the new revision and removed lines at the old one. That is the method above, and it gives the same 691.
   - *Is it still the only statement of the rule?* Yes. The paragraph still ends "This is the only statement of the rule; everywhere else points here". I searched every rule file for the ceiling. `.claude/agents/reviewer.md`, `developer.md` and `tester.md`, the three skills that mention it, and `docs/PRINCIPLES.md` E4 all point to §2 and do not repeat the number or the unit. Nothing conflicts with it. `docs/increments/README.md` does not state it. No `tools/` script counts lines, so no code can disagree with it.
   - *The gap @architect left open (a multi-line string that is not raw and not a docstring):* it does not matter for this PR. As worded, such lines are not blank lines, comments, docstrings or raw-string bodies, so they count. The only open question is a blank line inside such a string. In the PR's ten production files, all 103 multi-line strings are docstrings: a tokenizer pass found none elsewhere and no triple-quoted f-string. Leaving it open is fine until a PR has such a string.
3. **ROADMAP row 29, the status line, citations.**
   - Row 29 is still `ROADMAP.md:54`.
   - The row now says questions 1 to 4 were ruled on 2026-10-04 and confirmed on 2026-10-05, and 5 to 10 were ruled on 2026-10-05. That matches the markers in the increment file (lines 2974–3002 carry "Ruled 2026-10-04", lines 3010–3070 carry "Ruled 2026-10-05") and the paragraph at line 702.
   - The status line and the row both record round 5, its fix, and the net ruling, and both name round 6 as next. That matches the tree.
   - `python3 tools/check_citations.py` exits 0 and lists 105 citations to re-read. The live ones hold: `ROADMAP.md:54` is row 29, the `:3140` quoted inside the round-5 record is still the round-4 header line, and nothing cites line numbers in `CLAUDE.md` except a test fixture, so the three lines added to §2 move no citation.
4. **Production code did not change since round 4.** It is 691 net, as above.

**Other checks:** red-step comments: none possible, since no code changed. Mutation record: none owed this round; the one suite the increment names invariant-critical has its record with PR 1. `@perf` acceptance run: none owed, since no refine or mesh code is touched.

**Suggestions (non-blocking)**
- The round-5 fix paragraph (line 3295) says "691 net is therefore 8 under the ceiling". ROADMAP row 29 still says "691 production lines, 9 under the 700 ceiling". Both can be defended: 691 is 9 below 700, and there are 8 lines of room because the PR must stay *under* 700. A reader may still see a contradiction. "8 lines of room" in both places would remove it.
- `.claude/REQUIRED-READING.md` says to restart before spawning a persona after a change to `CLAUDE.md` is merged. This PR carries one. Once it merges, the main session should restart before its next spawn, so that every persona reads the net wording.

#### Pinned or assumed beyond the design
- I counted only `src_python/`. The whole PR touches no C++, and tests are excluded by §2.
- I identified docstrings with `ast` (the first statement of a module, class or function, if it is a string) and comment lines with `tokenize`. That is my reading of §2's exclusions, and it gives the figure earlier rounds recorded.
- I took the round-4 and round-5 handback texts from this session's transcripts of those runs, as the task directed.

#### Questions for Ola
none

#### Lessons
- A "word for word" claim can be checked with `diff` in a minute once the handback is pulled from the transcript. The diff can only fail if the heading change is applied to the record alone and the handback has headings for it to change. Both did here, so the empty diffs mean something.

#### ASK OLA lines
none

#### GUARD FALSE POSITIVE lines
none

Files: `/Users/skavhaug/projects/rasputin/.claude/worktrees/29-pr4/CLAUDE.md` (section 2), `/Users/skavhaug/projects/rasputin/.claude/worktrees/29-pr4/ROADMAP.md` (line 54), `/Users/skavhaug/projects/rasputin/.claude/worktrees/29-pr4/docs/increments/29-nve-reference-catchments.md` (line 3, records at lines 3140–3217 and 3219–3295). My note file is `/Users/skavhaug/projects/rasputin/.claude/current-task/reviewer-090320.md`; deleting it is the spawner's job.

Fixes for round 6, by `@architect` in the commit that records it: suggestion 1 taken: the round-5 fix paragraph above, `ROADMAP.md:54`, and the two estimate passages of this file that said "9 under 700" (under "Lake gauges' green step" and under "New and changed files") now all say PR 4's 691 net leaves 8 lines of room under the 700 ceiling. Suggestion 2, restarting before the next spawn once this PR's change to `CLAUDE.md` merges, is the main session's to act on, as `.claude/REQUIRED-READING.md` already requires; no rule changes. The status line and `ROADMAP.md:54` record PR 4 as approved in code review round 6, with the push waiting for Ola. No production line changed.

### Acceptance evidence, round 1: `@reviewer`, `bc01cd8..7e7e4ce`

`@reviewer`'s record, word for word:

**29 acceptance, evidence review, round 1, 2026-10-05 (copied from `@reviewer`'s handback).** Range `bc01cd8..7e7e4ce` (3 commits by `@perf`, 170 files, +3462/-0, all under `docs/benchmarks/2026-10-05/nve-hrd/`). Verdict: CHANGES REQUESTED.

**Size and CI.** Production lines: 0. The 874 added lines of `.py` and `.sh` are evidence scripts under `docs/benchmarks/`, which `pyproject.toml:117` excludes from the gates. The branch is not pushed: no PR, no CI. `tools/ci_changes.py bc01cd8 7e7e4ce` says `code=true`, so the push will run every CI job, not only the governance gates. `check_citations.py` and `check_prohibited_deps.py` pass.

**Re-measured and true:**
- **Classes.** 140 rows: match 74, close 5, miss 6, uncertain 39, refused 16. 85 stations are scored, and 74 of them (87 %) match.
- **Size bands.** By NVE's area: under 10 km²: 12; 10-100: 44; 100-1000: 75; over 1000: 9. The class counts in each band match the README. The share uncertain is 33/17/38/60 %.
- **Refusals.** By cause: `mixed_grid` 13, `other` 2, `no_river` 1. Each `mixed_grid` row names one tile from each grid.
- **Seeding.** Lake-seeded: 48 in all, 45 not refused, giving 41 match / 3 close / 1 miss / 3 refused, the same as the probe in "NVE's lakes". Every `uncertain` row is river-seeded.
- **Step 6 checks.**
  - The exact oracle (placed count equals catchment nodes) holds on 79 of 79.
  - `drains` is false on 35 rows and `monotone` on 32; the 32 are a subset of the 35, so the union is 35. All 35 are `uncertain`, and 22 of them have NVE's in ours under 10 %.
  - The burn lowered no node off the chain: 0 of 79.
  - Outlines simple and seed inside: 124 of 124. The largest area difference is 2.6e-6 m².
- **Uncertain rows.**
  - The 9 `drains`-false rows and the 11 uncertain rows with both overlaps at or above 95 % are the ones listed.
  - Causes: 35/25/9/8/0.
  - Swing-only rows are `2.284.0`, `38.1.0`, `124.2.0` and `237.1.0`.
  - 11 rows have `confluence_near`, and one is on a lake line (`22.16.0`).
  - Bypass: 17 of the 22 have a node 10 times the placed count, and the 5 exceptions are the ones named.
- **Reruns.** All 24 re-runs of the six misses kept their class. Node counts are within 0.0098 %. `83.2.0` is placed by `name` at 250 m.
- **Window check.** Only `105.1.0` grows in the ±12 km window.
- **Generated files.** `analyse.py .` regenerates `analysis.md` byte for byte. The 124 committed outlines, `results.csv` and `summary.json` are byte-identical to the scratch batch output.
- **Totals.** 194.5 M nodes; median 4.9 s per station; the longest is `234.13.0` at 145 s.
- **`105.1.0`.** `windows_105.1.0.txt` shows the last window at 9.3 × 9.3 km with `edge False`. `trace_exit` shows the path leaving at column 941 of 942 over heights 11.69-11.89 m, 3.27 km east. `chain_counts_by_window` gives the placed node 0.75, 30.3, 56.5, 66.6 and 119.7 km² at ±2.6/6/8/10/12 km. So the window does decide the flow.
- **Sagafoss.** The fourth window is clear of its edge. Its in-nodes reach x 849.77 km, plus 2 km is 851.8 km, past the window's 850.9 km. The margin then doubles to 32 km, and the refusal is raised at `catchment.py:259`. NVE's polygon lies on `7708_4` alone and is 16 km from the nearest shifted tile (`7808_3`, by the world files).
- **Knappom.** 368.87 km plus a 32 km margin is 400.87 km, which matches the refusal message. The data end at 400.26 km, 26 km past NVE's x of 374 km.

**Blocking:**
1. **The lake-above-P list is incomplete by construction, and the README's reason is false.** `lake_above.py` reads `lakes.geojson`. The fetch takes lakes only within ±100 m of each station (`fetch/nve.py:55`, `LAKE_ENVELOPE_HALF = 100.0`), so the file holds 81 lakes. A lake 100-1000 m up the reach is missing unless its polygon comes within that envelope. Two of the four stations the design names (`:2030-2032`), `2.633.0` and `55.4.0`, are missing: the nearest lake in the file is 100 km and 51 km away. The README says "the probe's rule is not recorded", but `:952-954` records it: a ±1000 m query of layer 5, and a lake on the reach within `reach_up` above `P`. Ola's ruling on question 10 (`:3070-3072`) asks for "the full run's list of these 24". Required:
   - re-run the list with lakes queried along each reach (the design's ±1000 m envelope), giving the classes of all 24, or of as many as the rule finds;
   - replace "not looked at" with the cause;
   - update the uncertain section's "five have a lake on the reach above P".
2. **Sagafoss (`212.48.0`) is listed under "Known and expected refusals … (not failures)".** Acceptance step 4 (`:2912-2913`) says: "A refusal for any other reason is a finding, including a window that reaches a shifted tile where the catchment does not." That is this case, and it is not one of the four stations the design allowed into the known count. Move it to "Not as expected" beside Knappom and call it a finding. The known two-grid refusals are then 12, plus one finding.
3. **Knappom: "the catchment's nodes stayed inside x 329-369 km" is not what the probe shows.** In the fourth window `edge True`: the in-nodes end at x 368.87 km, at the window's own east edge (368.88). So the catchment was cut by the window, not settled. Say so: a fifth window was needed (NVE reaches 374 km), and its 32 km margin went past the data's edge at 400.26 km where about 6 km would have stayed inside. The conclusion ("the window rule refuses, not the catchment") stands.
4. **`run.sh` contradicts the commit.** Its header says the catchment files go to `$WORK/batch` "(not committed; regenerate with this script)", but `catchments/` is committed, and no script copies it.
   - Say the outlines are copied from `$WORK/batch` and by which command, or add the copy to `run.sh`.
   - Its last comment, "Steps 3 and 6", runs only `analyse.py`, before `checks.py` exists. Name the order: `run.sh`, then `checks.py`, then `analyse.py` again.
5. **The mechanism in "The burnt path does not hold" names the wrong rule.** It says the 1 mm drop "does not stop the DEM's steepest descent from choosing a lower node beside the path". The flow here is not steepest descent: `flow_to` is each node's flooder in the Priority-Flood (`accumulate.hpp:23`, `upstream.hpp`). A node drains to the neighbour that first pushes it, the lowest one already flooded, and that neighbour can be off the path. Reword to match the code.
6. **Acceptance step 3 asks for the summary by tile count in the README.** It is only in `analysis.md`. Add the table, or a line pointing to it.
7. **The status line (`docs/increments/29-nve-reference-catchments.md:3`) and `ROADMAP.md:54` must change before the push.** Both still say the pushes of PR 2 and PR 4 wait for Ola, though they merged as #179 and #182, and neither says the acceptance run is in. This is the same finding as evidence round 1 of PR 2. The run is `@perf`'s evidence; these lines are `@architect`'s.

**Suggestions (non-blocking):**
- **`196.11.0`.** It is one of the nine expected two-grid refusals, but it ran (22 nodes). It belongs under "Not as expected" beside Polmak, with the reason: the catchment never grew to the shifted tiles.
- **The window check.** It is cut by its own ±12 km edge on 33 of the 79 rows (`edge_cut`), so it could see a cut-short catchment on 46 at most. Say "1 in 46 testable" rather than "one in 79 is the measured rate".
- **`105.1.0` at ±12 km.** The count is still cut there (`edge_cut` 1), so 119.7 km² is a lower bound. Say so.
- **`2.279.0`.** The path also carries nothing again at +42 m (`chain_counts_2.279.0.txt:29`), so "rejoins it below the gauge" is only part of the picture.

**On the 2.4 MB of outlines.** It is in order. They are our reduced outlines (ours, derived from DTM10, with no NVE data), byte-identical to the batch output. Earlier evidence committed reduced outlines (`2026-09-29/bygdin/`), and `2026-10-04` is 6.9 MB. They let step 6 be re-checked and PR 5's comparison be run without the batch. Only blocker 4's wording is needed.

Not pushed; no CI.

Fixes: items 1-6 and the four suggestions by `@perf` in `cd8b5be` (the window check worded as conclusive on 46, finding nothing there, and as finding `105.1.0` among the 33 edge-cut rows); item 7 by `@architect`.

### Acceptance evidence, round 2: `@reviewer`, `7e7e4ce..d9bbfe1`

`@reviewer`'s record, word for word:

**29 acceptance, evidence review, round 2, 2026-10-05 (copied from `@reviewer`'s handback).** Range `7e7e4ce..d9bbfe1` (3 commits: `cd8b5be` by `@perf`, `0130889` by the main session, `d9bbfe1` by `@architect`; 14 files, +457/-107). Verdict: APPROVED.

**Size and CI.**
- Production lines: 0. The scripts that changed are evidence scripts under `docs/benchmarks/`, and the only other files are `ROADMAP.md` and this increment file.
- The branch is not pushed, so there is no PR and no CI. `tools/ci_changes.py bc01cd8 d9bbfe1` says `code=true`, so the push will run every CI job.
- `check_citations.py` exits 0. I re-read its at-risk list as quotations:
  - `ROADMAP.md:54` is still row 29.
  - `:3` is still the status line.
  - `chain_counts_2.279.0.txt:29` is still the +42.4 m row that carries 0.000 km².
- `check_prohibited_deps.py` passes.
- No mutation record is owed, and no `@perf` acceptance run of `bench.py`: this round is evidence, not refine or mesh code.

**Round 1's items, each checked:**
1. **Lakes along each reach (fixed).**
   - I re-ran `lake_query.py` against NVE today. It made 93 queries and returned 187 lakes. Its output file has SHA-256 `480b7130…`, exactly the value in `provenance.txt`, and its printout is byte-identical to `step4/lake_query.txt`.
   - Only `109.9.0`'s reach leaves its ±1000 m box (by 111 m), and that reach meets no lake.
   - `lake_above.py` on that output reproduces `step4/lake_above.txt` byte for byte: 26 stations. The probe can fail: run on the fetch's own `lakes.geojson` (81 lakes) it gives 16.
   - From `results.csv`, both populations are 91: the stations placed on a river line, and the stations seeded on the river. Their union is 93.
   - On the river-line population the rule finds 24, the design's count, with `2.633.0` (Skjølja) and `55.4.0` (Røykenesvatnet) among them. All four of PR 2's full runs are in it, at 97.1 % or better both ways.
   - On the river-seeded population it finds 25: 11 match, 4 miss, 7 uncertain, 3 refused. That is the README's table name for name. The difference is `12.197.0` and `22.16.0` (placed on a lake line, seeded on the river) against `62.18.0` (on a river line, seeded by its lake).
   - The uncertain section's "seven have a lake on the reach above P" matches the seven names.
2. **Sagafoss (fixed).** `212.48.0` is now under "Not as expected" as a finding, so the known two-grid refusals are 12 and the expected refusals 14.
   - `windows_212.48.0.txt`: fourth window `edge False`, in-nodes x 819.41-849.77 km.
   - NVE's polygon is x 819.41-849.73, y 7766.77-7783.68 km, as stated.
3. **Knappom (fixed).** `windows_2.142.0.txt` shows `edge True` in all four windows.
   - The fourth window runs x 321.40-368.88 km, and its in-nodes reach x 368.87 km.
   - The margin doubles each window (4, 8, 16 km), so the fifth window's 32 km margin reaches x 400.87 km. That is the refusal's x range, against data that end at 400.26 km.
4. **`run.sh` (fixed).** It now copies the outlines into `catchments/`, then runs `analyse.py`, `lake_query.py`, `checks.py` and `analyse.py` again, in that order. The header and the README's "Commands" row say the copy, `checks.py` and the second `analyse.py` were run by hand on the day, with the same commands.
   - `checks.py` re-run here with the new lakes (on `results.csv` and the 124 committed outlines) gives a `checks.csv` byte-identical to the committed one.
   - Against `7e7e4ce`, only the `lake_above_p` column changed, on 9 stations (`2.633.0`, `12.215.0`, `18.10.0`, `19.104.0`, `55.4.0`, `79.3.0`, `105.1.0`, `153.1.0`, `237.1.0`): 13 + 9 = 22 river rows not refused.
   - `analyse.py` and `findings.py` regenerate `analysis.md` and `step4/uncertain_table.md` byte for byte.
5. **The flow rule (fixed).** The new wording matches `flood.hpp`:
   - keys are (level, push counter), so the lowest level pops first and equal levels pop first in, first out;
   - `on_reach(i, j)` names the popped node `i` that reached `j` first;
   - `accumulate.hpp` stores that node as `flow_to`, "its flooder".
   - "Rejoins it, more than once" holds in `chain_counts_2.279.0.txt`: 45.9 km² again at +14 m and at +250 m, after nothing at -24 to 0 m and at +42 m.
6. **The tile-count table (fixed).** The README's table is `analysis.md`'s, and `results.csv` gives the same counts. All 22 uncertain rows under a tenth of NVE's area meet one tile.
7. **Status line and ROADMAP (fixed).** The status line and `ROADMAP.md:54` match `summary.json`: 74/5/6/39/16 and 14 expected refusals. `batch.log` gives 1994.36 s (33 minutes) and a peak of 7,745,503,232 bytes (7.7 GB). PRs #173, #179, #178 and #182 are named as merged.

**Round 1's suggestions:**
- `196.11.0` has its own line: 22 nodes, 0.002 km² of 637, one tile, two windows.
- `105.1.0`'s 119.7 km² is called a lower bound.
- `2.279.0`'s second loss of flow is in.
- **The window check.** `@perf` was right not to use "1 in 46". In `window_check_12km.csv`, 33 rows are edge-cut and 46 are not. All 46 have a ratio of exactly 1.0, and the only row that grows (`105.1.0`, ratio 160) is one of the 33. So "1 in 46" would have been false. "Conclusive on 46, finding none; one of the 33 cut; 1 in 79, with 32 not cleared" is exact.

**Suggestions (non-blocking):**
- **Knappom.** "A margin of about 6 km would have covered NVE's 374 km" is true in x only, and the wording came from round 1. NVE's polygon reaches y 6790.47 km. That is 34.7 km north of the fourth window's in-nodes (6755.75 km, at that window's north edge), and past even the fifth window's 6787.75 km. Say that the catchment still had to grow north. The margin grows on all four sides, so the east side ran off the data while the catchment needed more room only to the north. The conclusion, that the window rule refuses and not the catchment, stands.
- **The flow rule.** "The flood rises from the window's edge": `flood.hpp` also makes every node beside NoData an outlet. Add "or beside NoData".
- **The status line and `ROADMAP.md:54`** say "its evidence review round 2 is next". Update both when this round is recorded.
- **`ROADMAP.md:54`** no longer gives the plan's outline: the five PRs and their estimates, including PR 5's 80 lines, and the expected refusals. They are now only in the increment file. Restore one clause if the row should still show PR 5's size.

Not pushed; no CI.

Taken in the recording commit: the status line and `ROADMAP.md:54` now say round 2 approved and the push waits for Ola. The Knappom, flow-rule and ROADMAP-outline suggestions are left open.
