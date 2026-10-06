# Increment 29 acceptance: all 140 NVE reference stations

`@perf`, 2026-10-05. Steps 1-4, 6 and 7 of "Acceptance: every covered HRD
station" in `docs/increments/29-nve-reference-catchments.md`. Step 5 (the
nearest-stream fallback) waits for PR 5. HRD is NVE's Hydrological
Reference Dataset, the 140 gauging stations listed in
`src_python/tin_engine/data/nve_hrd_2025.csv`.

## At a glance

One run of `rasputin station-catchments` over all 140 stations, with
`--lakes`, at the defaults. **Every station got a row; the run took 33
minutes (1994 s) and at most 7.7 GB of memory.**

| class | stations | what it means |
|---|---:|---|
| match | **74** | both overlaps with NVE's polygon at least 95 % (all 74 by overlap; none needed the 30 m divide-offset test) |
| close | 5 | both overlaps at least 80 % |
| miss | 6 | anything else |
| uncertain | 39 | the area changes too much near the gauge, or the burnt river path does not hold; not scored |
| refused | 16 | no catchment: 13 on tiles of two grids (12 known, 1 a finding), 2 off the DEM's edge (1 known, 1 a finding), 1 with no river line |

Scored stations (match, close, miss) are 85; 74 of them (87 %) match.

**Known and expected refusals, apart (not failures): 14.**

- 12 stations refused because their windows select tiles on two different
  grids, which rasputin does not combine, and neither grid covers the
  window alone: `156.15.0`, `156.24.0`, `191.2.0`, `203.2.0`, `206.3.0`,
  `208.2.0`, `208.3.0`, `209.4.0`, `212.49.0`, `213.2.0`, `213.4.0`,
  `223.2.0`. Eight tiles of the DEM sit half a cell (5 m) east of the
  others; a later increment resamples them onto the common grid. Of the
  nine the design expected, eight are here; the ninth, `196.11.0`, ran
  (below, "Not as expected"). The four the design named as likely
  (`156.24.0`, `213.4.0`, `191.2.0`, `203.2.0`) are all here. A thirteenth
  refused on two grids, `212.48.0` Sagafoss, is not a known refusal: the
  design's step 4 makes a window that reaches a shifted tile where the
  catchment does not a finding (below).
- `244.2.0` Neiden: refused, its catchment needs DEM nodes where there is no
  tile (Finland; 131 km² of NVE's polygon lies outside every tile). Expected.
- `311.4.0` Femundsenden: no river line within 500 m. Expected until PR 5.

**Not as expected: two refusals that are findings, two runs expected refused.**

- `212.48.0` Sagafoss, refused on two grids, is a finding: its NVE polygon
  lies on one normal tile, 16 km from any shifted one. Measured window by
  window (`probe_windows.py`, `step4/windows_212.48.0.txt`): the fourth
  window's catchment is clear of its edge (`edge False`) and spans x
  819.4-849.8 km, y 7766.8-7783.3 km (NVE's polygon: x 819.4-849.7, y
  7766.8-7783.7). Its in-nodes plus 2 km still reach past that window (x
  851.8 against 850.9 km), so the loop grows it, by the margin doubled to
  32 km, and that window's tiles include the shifted `7707_1` (the call
  stack: `_grow`, `src_python/tin_engine/catchment.py@44fa7f5:259`, then `_plan`). So the catchment was
  found, on one grid, and the window rule refused it.
- `2.142.0` Knappom was refused off the DEM's edge although its NVE polygon
  lies 26 km inside the tiles. Measured window by window
  (`probe_windows.py`, `step4/windows_2.142.0.txt`): every one of the four
  windows touches its edge (`edge True`); in the fourth, x 321.4-368.88 km,
  the catchment's nodes reach x 368.87 km, one cell from the window's east
  edge, so a fifth window was needed (NVE's polygon reaches x 374 km). That
  window's margin, 32 km (doubled each window), took it to x 400.87 km,
  past the data's east edge at 400.26 km; a margin of about 6 km would have
  covered NVE's 374 km and stayed 25 km inside the data. So the window rule
  refuses, not the catchment. A finding.
- `234.18.0` Polmak nye was expected to be refused (its polygon crosses into
  Finland) but ran: the gauge's burnt river path does not carry the
  river's flow, so the flood stayed at 6.3 km² of NVE's 14,171 km² and
  never reached the border. It is `uncertain` (chain not draining).
- `196.11.0` Lille Rostavatn was expected to be refused on two grids but
  ran, in the same way: its catchment is 22 nodes (0.002 km² of NVE's 637),
  on one tile after two windows, so it never grew to the shifted tiles. It
  is `uncertain` (swing, chain not draining).

**By size band** (NVE's polygon area):

| band, km² | stations | match | close | miss | uncertain | refused | share uncertain (of not refused) |
|---|---:|---:|---:|---:|---:|---:|---:|
| under 10 | 12 | 5 | 3 | 0 | 4 | 0 | 33 % |
| 10-100 | 44 | 29 | 2 | 3 | 7 | 3 | 17 % |
| 100-1000 | 75 | 38 | 0 | 3 | 25 | 9 | 38 % |
| over 1000 | 9 | 2 | 0 | 0 | 3 | 4 | 60 % |

**By tile count** (the DEM tiles our fine outline meets; refused stations
have none, so the 124 not refused; from `analysis.md`):

| tiles | stations | match | close | miss | uncertain | share uncertain |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 81 | 40 | 4 | 4 | 33 | 41 % |
| 2 | 26 | 22 | 0 | 1 | 3 | 12 % |
| 3-4 | 17 | 12 | 1 | 1 | 3 | 18 % |
| 5 or more | 0 | 0 | 0 | 0 | 0 | |

All 22 `uncertain` catchments under a tenth of NVE's meet one tile: they
are 22 of the one-tile row's 33 `uncertain`.

**Lake gauges against river gauges.** Seeded with their lake (45 not
refused): 41 match, 3 close, 1 miss (`83.2.0`, expected). Seeded on the
river (79 not refused): 33 match, 2 close, 5 miss, **39 uncertain**. Every
`uncertain` is a river gauge.

**What to look at first.**

1. **The 35 river gauges whose burnt path does not drain** (`drains` false
   in `results.csv`): they are 35 of the 39 `uncertain`, and 22 of them got
   a catchment under a tenth of NVE's (`NVE's in ours` under 10 %), several
   of a few nodes (`196.11.0`: 22 nodes; `311.6.0` Nybergsund: 0.00 km² of
   4,418). This is the largest single loss: the burn in the DEM does not
   make water follow the mapped river at these gauges. Nine of them have
   both overlaps at or above 95 % all the same (`12.70.0`, `22.16.0`,
   `22.22.0`, `73.21.0`, `83.12.0`, `112.8.0`, `122.11.0`, `122.14.0`,
   `205.6.0`): the path fails within the position uncertainty, but the
   catchment from the placed node is right. With `2.284.0` and `38.1.0`
   (swing only), eleven `uncertain` stations have both overlaps at or
   above 95 %.
2. **The window can decide the flow** (`105.1.0` Osenelv, 0.75 km² of
   NVE's 138): across flat ground, the river's water leaves through the
   window's edge, so the catchment is clear of the edge and the window
   stops growing, though a ±12 km window gives the same node at least
   119.7 km² (a lower bound: that window still cuts the catchment). One
   case in the 79 checked, and the check is conclusive on only 46 of them
   (33 are cut by the ±12 km window's edge, `105.1.0` among them); it is
   increment 22's window loop, not this increment's
   code. And `153.1.0` Storvatn (2.7 km² of 49): the DEM drains
   the lake Storvatnet 1 km away from the mapped river. Both passed the
   sensitivity check, so they are `miss`, not `uncertain`.
3. **Three catchments larger than NVE's** with NVE's polygon almost wholly
   inside ours (98.5-98.8 %): `15.49.0` Halledalsvatn (ratio 1.73),
   `12.215.0` Storeskar (1.32), `307.7.0` Landbru (1.28).
4. `2.142.0` Knappom and `212.48.0` Sagafoss: refused by the window rule
   (the growth margin off the DEM's edge; a shifted tile in the window),
   not by their catchments, which were inside the data and on one grid.
   Both are findings.

## Findings, one line each (step 4)

Measured with the scripts named; "not looked at" where nothing was measured.
"Ours" is our catchment, "NVE's" NVE's polygon.

### The six `miss` rows

Each was re-run at corridor 15 m and 60 m and at map radius 250 m and
1000 m (`rerun.sh`; `reruns/*.csv`). **No class changed in any of the 24
re-runs**, and the catchments are the same to within 0.01 % of their nodes:
the corridor changes the burnt path (other chain and lowered-node counts in
every river row), not the catchment. At map radius 250 m, `83.2.0` is placed
on another line (by name, not number), and is seeded by the same lake.

| station | name | NVE's in ours / ours in NVE's | explanation |
|---|---|---|---|
| `12.215.0` | Storeskar | 98.7 % / 74.7 % | **NVE's polygon disagrees with the DEM.** Ours holds one extra piece of 39.5 km², 12 km from the station, far beyond the 1 km of river the burn touches; along the 13.5 km where it meets NVE's divide the ground is not flat (1327-1697 m, no sample flat; `divide.py`). Why NVE's polygon leaves it out is not looked at |
| `15.49.0` | Halledalsvatn | 98.8 % / 57.1 % | **NVE's polygon disagrees with the DEM.** One extra piece of 43.4 km², 2.2 km from the station; 7 % of the 7.4 km shared edge flat. The station is 47 m from the lake Halldalsvatnet but placed on the river line, so it is seeded on the river (a lake on the reach above P). Why NVE's polygon leaves the piece out is not looked at |
| `83.2.0` | Viksvatn (Hestadfjorden) | 3.8 % / 99.3 % | **Lake gauge, expected** ("NVE's lakes"): the station point lies inside Hestadfjorden, whose own catchment (19.3 km²) is a small part of the gauge's 508 km². Why the point is in that lake is not looked at |
| `105.1.0` | Osenelv v/Øren | 0.4 % / 77.6 % | **Not one of the listed causes: the window decides the flow.** The batch's last window (9.3 × 9.3 km round the reach; `probe_windows.py`) gives the placed node 0.75 km², clear of the window's edge, so the window stops growing. In a ±12 km window the same node carries 119.7 km² (NVE: 138.1; `chain_counts.py`), a lower bound, as that window still cuts the catchment (`edge_cut` 1 in `window_check_12km.csv`). Traced (`trace_exit.py`) in a ±4.7 km window, which gives the same 0.75 km²: the river's water 1.5 km above the gauge (118.7 km² in the large window) runs 3.3 km east across ground flat to 0.2 m (11.7-11.9 m) and leaves through the window's edge, which the flow routing treats as an outlet. So "clear of the edge" does not prove the catchment complete. See "The window and flat ground" below |
| `153.1.0` | Storvatn | 5.4 % / 98.0 % | **The DEM drains the lake another way than the mapped river.** The station is 5 m from Storvatnet (3.8 km²) but outside it and on a river line, so it is seeded on the river. Within 300 m of P the DEM's largest count is 3.0 km²; the DEM's river with 48.2 km² (NVE: 49.1) runs 971 m from the mapped line (`explain.py`, ±12 km window). The same 2.7 km² at window sizes 4.7 to 12 km, so not the window |
| `307.7.0` | Landbru | 98.5 % / 77.3 % | **NVE's polygon disagrees with the DEM.** One extra piece of 17.0 km², 4.7 km from the station; the 9.6 km shared edge is not flat (no sample). Why is not looked at |

None is "placed on another river": all six were placed on a line with the
station's own watercourse number (`placed_on` = `number`), at 1 to 20 m
(`83.2.0`: 322 m, then seeded by its lake). None is a burn artefact of
the kind "many lowered nodes": the extra or missing areas lie 2 to 12 km
from the station, past the 1 km the burn changes, or (`105.1.0`) the burn
changes nothing at the placed node (raw and burnt counts equal there).

### The window and flat ground (from `105.1.0`)

`window_check.py` re-reads the placed node's count in a ±12 km window for
all 79 river-seeded stations not refused (`window_check_12km.csv`). Only
`105.1.0` gets a larger count there than in the batch (whose catchments are
all clear of their windows' edges by construction). The check sees a
catchment cut short only when the ±12 km window holds more of the true
one, so on catchments larger than that window it can miss one: 33 of the
79 are cut by the ±12 km window's edge (`edge_cut` 1). On the 46 not cut
the check is conclusive and finds none; among the 33 cut it finds one,
`105.1.0` (its 119.7 km² a lower bound), and clears none of the other 32.
So the measured rate is 1 in 79, with 32 not cleared; the cause (the window's edge as an
outlet across a flat) is a property of the window loop of increment 22,
not of this increment's code.

### The five `close` rows, by cause

| cause | stations |
|---|---|
| a divide differs on ground that is not flat (largest differing piece 0.6 to 9.5 km², 0 to 8 % of the shared edge flat; `divide.py`) | `6.10.0` Gryta (ours larger by 0.64 km², 2.5 km from the station), `19.79.0` Gravå (NVE's larger by 0.95 km², 0.3 km from the station), and the lake rows `26.29.0`, `35.9.0` and `16.66.0` below |

### The lake rows (seeded by their lake), apart

Of 48 lake-seeded: 41 `match`, 3 `close`, 1 `miss`, 3 refused on tiles of two
grids (`191.2.0`, `203.2.0`, `213.2.0`). The same counts as the lake seed's
probe in "NVE's lakes".

| station | name | class | NVE's in ours / ours in NVE's | where the difference lies |
|---|---|---|---|---|
| `26.29.0` | Refsvatn | close | 99.3 % / **84.0 %** | **a divide**: one extra piece of 9.5 km², 5.0 km from the station, its 6.8 km edge with NVE's polygon on ground that is not flat (no sample flat). It is far from the outlet, so neither a polygon past the outlet nor a station on an inflow |
| `35.9.0` | Osali (Botnavatnet) | close | 98.7 % / **93.3 %** | **a divide**: one extra piece of 1.16 km², 1.7 km from the station, not flat |
| `16.66.0` | Grosettjern | close | **93.6 %** / 98.7 % | NVE's larger by 0.42 km² in all, the largest piece 0.15 km² at the station (0.0 km): NVE's polygon reaches further round the outlet than ours; the station is 25 m inside NVE's edge. A 6.5 km² catchment, so 0.42 km² is 6 % |
| `83.2.0` | Viksvatn (Hestadfjorden) | miss | **3.8 %** / 99.3 % | expected (above): the lake is not the gauge's whole catchment |

### River rows with a lake on the reach above P (question 10)

The rule, `@architect`'s probe's ("NVE's lakes" in the increment file):
the mapped reach from its upstream end (up to `reach_up`, 1000 m above P)
to P meets one of NVE's lake polygons (layer 5 of `Innsjodatabase2`).
`lake_above.py`, output `step4/lake_above.txt`, which names each lake.

**The input.** The fetch's `lakes.geojson` holds only the lakes within
±100 m of each station point (`fetch/nve.py`, `LAKE_ENVELOPE_HALF`, sized
for the lake seed's 30 m gap): 81 lakes. A lake higher up the reach is not
in it, so the first version of this list, read from that file, found 16,
not the design's 24; that was the cause. `lake_query.py` now queries layer
5 once per station seeded on the river or placed on a river line (93),
with the envelope that holds the station point ± 1000 m (the probe's) and
the whole reach above P (`step4/lake_query.txt`: every reach lies inside
its ±1000 m box but one, `109.9.0`, 111 m past it, with no lake on it).
187 lakes, no answer truncated, 2026-10-05 12:00 CEST; the file is NVE's
data, not committed (`$WORK/lakes_reach.geojson`, SHA-256 `480b7130…`, in
`provenance.txt`). `step4/lake_above.txt` marks each lake the ±100 m file
lacks.

**The count.** On the probe's population, the 91 stations placed on a
river line, the rule finds **24**, the design's count, `2.633.0` and
`55.4.0` among them. The question 10 asks is about the stations that keep
the river method, the 91 seeded on the river; that population differs by
three (`12.197.0` and `22.16.0` are placed on a lake line but seeded on the
river, as their station lies 52 and 460 m from the lake their P lies in,
past the seed's 30 m; `62.18.0` Svartavatn
is placed on a river line but seeded by the lake it lies in), and on it the
rule finds **25**:

| class | stations |
|---|---|
| match (11) | `2.633.0` Stortorp, `12.197.0` Grunke, `16.127.0` Viertjern, `18.10.0` Gjerstad, `18.11.0` Tjellingtjernbekk, `19.104.0` Songedalsåi, `55.4.0` Røykenes, `101.1.0` Engsetvatn, `148.2.0` Mevatnet, `168.3.0` Lakså bru, `189.3.0` Tennevikvatn |
| miss (4) | `12.215.0` Storeskar, `15.49.0` Halledalsvatn, `105.1.0` Osenelv, `153.1.0` Storvatn (all four above) |
| uncertain (7) | `12.178.0` Eggedal, `12.188.0` Langtjernbekk, `19.96.0` Storgama ovf., `22.16.0` Myglevatn ndf., `79.3.0` Nessedalselv, `83.12.0` Haukedalsvatn ndf., `237.1.0` Båtsfjord |
| refused, two grids (3) | `156.24.0` Bogvatn, `212.49.0` Halsnes, `223.2.0` Lombola |

Of the 22 not refused, 11 match, 4 miss, 7 are uncertain. The four misses'
causes as measured above: `12.215.0`'s and `15.49.0`'s extra pieces lie 12
and 2.2 km from the station (ours holds 98.7 and 98.8 % of NVE's),
`105.1.0`'s water leaves through the window's edge, and `153.1.0`'s DEM
river (48.2 km²) runs 971 m from the mapped line. Whether seeding the lake with the river (question 10's alternative) would
change any of the 25 is not measured. In the table under "Every
`uncertain` station" and in `analysis.md`, `lake above P` now reads these
lakes (`checks.py`'s fifth argument); re-run so, `checks.csv` changed in
that column only, on nine stations.

### The `uncertain` rows

All 39 are river-seeded. Their causes, counted (a station can have
several): chain not draining 35, swing 25, downstream side not read 9,
chain end open 8, line against the slope 0.

- **The burnt path does not hold (35).** Exactly the 35 with `drains`
  false; 32 of them also fail `monotone`. 22 of them got a
  catchment under a tenth of NVE's. For 17 of those 22, in the burnt DEM a
  node at most 3 nodes from the placed one carries at least ten times its
  count (`bypass.py`, ±6 km window, `bypass.csv`): the river's water passes
  beside the placed node. Measured along one path (`2.279.0` Kråkfoss;
  `chain_counts.py`, ±6 km window): in the burnt DEM the path node 34 m
  above the placed node carries 45.9 km², the three path nodes from 24 m
  above down to the placed node carry nothing, the nodes 14 and 28 m
  below carry 45.9 km² again, and the node 42 m below carries nothing
  again: the water leaves the path and rejoins it, more than once. The
  flow is the Priority-Flood's (`flood.hpp`, `accumulate.hpp`,
  `upstream.hpp`): the flood rises from the window's edge, popping the
  lowest level first, and each node drains to its flooder, the neighbour
  whose pop reached it first. The burn lowers each path node only 1 mm
  (`DROP_M`) below the one above it, which does not make the next path node
  that flooder: a neighbour beside the path that pops earlier (a lower
  level, or an equal one pushed earlier) floods the path node first, and
  the path node's water goes there. That this is what happens at
  `2.279.0`'s three empty nodes is read from the counts, not traced node by
  node. In 5 of the 22 (`41.8.0`, `87.10.0`, `88.11.0`, `212.10.0`,
  `234.18.0`) no node that close carries ten times the placed node's
  count; where the DEM's river runs there is not looked at.
- **Swing alone (4)**: `2.284.0`, `124.2.0`, `237.1.0`, and `38.1.0` (with
  the downstream side not read). In each, the whole catchment joins the
  burnt path between 20 m above the placed node and the node itself (the
  largest step equals the catchment, at -20 to 0 m), and the path's upper
  end at `U` carries almost nothing (0.0002 to 0.004 km²): the mapped line
  above the gauge is not where the DEM's river runs. Two of the four have
  both overlaps at or above 95 % (`2.284.0`, `38.1.0`).
- **Lake or flat**: one of the 39 is on a lake line (`22.16.0`), seven have
  a lake on the reach above P (`12.178.0`, `12.188.0`, `19.96.0`,
  `22.16.0`, `79.3.0`, `83.12.0`, `237.1.0`; list above). A confluence (another river line
  within 100 m of P, `confluence_near`) is near 11 of them.

The full table, one row per station, with the largest step and where, `U`,
the lake flags and the bypass count, is "Every `uncertain` station" below.

## The checks that need no NVE data (step 6)

On the 124 stations not refused (79 seeded on the river, 45 with a lake).
The burn check re-runs `burn.burn_reach` on stage A's first window
(`checks.py`); the rest read the run's own rows and catchment files.

| check | failures |
|---|---:|
| the placed node's count equals the catchment's node count (the exact oracle; a failure would be a defect) | **0** of 79 |
| the counts along the chain rise strictly downstream and the end is closed (`monotone`) | 32 of 79 |
| each chain node drains into the next (`drains`) | 35 of 79 |
| the burn lowered no node off the chain | **0** of 79 |
| reduced outline simple (increment 22) | **0** of 124 |
| seed strictly inside the reduced outline (increment 22) | **0** of 124 |
| area kept, reduced against fine within 1e-9 (increment 22) | **0** of 124; the largest difference 2.6e-6 m² |

The stations whose burn did not hold (`monotone` or `drains` false) are
35, the same 35 as `drains` false; they are 35 of the 39 `uncertain`. The
other four are `uncertain` by swing alone (`2.284.0`, `124.2.0`,
`237.1.0`) or swing and the downstream side not read (`38.1.0`).

## Step 7: does the increment pass?

| criterion | result |
|---|---|
| the batch completes for every station, a row each | yes: 140 rows |
| 22's guarantees on every accepted reduced outline | yes: 124 of 124 |
| every finding has its line | yes: the six `miss`, five `close`, 39 `uncertain`, the lake rows and every refusal ("Findings") |
| step 6's checks hold | the oracle and the burn's no-off-chain check hold everywhere; 35 stations whose burn did not hold are listed, as step 6 asks |
| a share of `match` | none required (Ola, 2026-10-04); this run is the baseline |

## Step 4 status

Done: every `miss` explained and re-run at the four settings; every
`uncertain` with its causes; the `close` rows by cause; the lake rows apart;
the river rows with a lake above P; every refusal. Not done: none of
step 4 is left. Step 5 waits for PR 5.

### Every `uncertain` station, with its causes

Causes from the sensitivity (`causes`): `swing` the area changes over 5 % within the position uncertainty `U` up or down the river; `downstream_unread` the river was not read to `U` below the gauge; `chain_not_draining` a node of the burnt path does not drain into the next, or the counts do not rise; `chain_end_open` the path's lowered end found no way down within the cap. Largest step: the largest rise in area between two path nodes within `U`, and where (m, + downstream). Bypass: in the burnt DEM, the largest count within 30 m of the placed node, against the placed node's own (`bypass.py`).

| station | name | NVE km² | ours km² | NVE's in ours % | causes | largest step km² at m | U m | lake line / lake above P | bypass km² |
|---|---|---:|---:|---:|---|---|---:|---|---:|
| 2.279.0 | Kråkfoss | 434.5 | 0.000 | 0.0 | swing, chain_not_draining, chain_end_open | 4.89 at 14 | 84 |  | 45.96 |
| 2.284.0 | Sælatunga | 457.5 | 460.994 | 99.1 | swing | 460.99 at 0 | 30 |  | 41.06 |
| 2.303.0 | Dombås | 493.9 | 0.003 | 0.0 | swing, downstream_unread, chain_not_draining | 0.00 at -14 | 30 |  | 26.06 |
| 3.22.0 | Høgfoss | 299.8 | 330.409 | 99.2 | chain_not_draining, chain_end_open | 330.41 at 38 | 40 |  | 24.56 |
| 12.70.0 | Etna | 568.5 | 565.452 | 98.8 | swing, downstream_unread, chain_not_draining | 565.45 at 0 | 130 |  | 54.37 |
| 12.178.0 | Eggedal | 310.6 | 0.000 | 0.0 | swing, chain_not_draining | 309.97 at 28 | 32 | lake above P | 34.07 |
| 12.188.0 | Langtjernbekk | 4.7 | 0.001 | 0.0 | swing, chain_not_draining | 0.00 at 24 | 31 | lake above P | 4.71 |
| 19.80.0 | Stigvassåi | 14.5 | 15.197 | 98.9 | chain_not_draining | 0.00 at 0 | 30 |  | 15.20 |
| 19.96.0 | Storgama ovf. | 0.6 | 0.000 | 0.1 | swing, chain_not_draining | 0.00 at 28 | 30 | lake above P | 0.60 |
| 20.2.0 | Austenå | 277.2 | 290.883 | 99.6 | chain_not_draining | 290.88 at 0 | 30 |  | 36.19 |
| 20.11.0 | Tveitdalen | 0.4 | 0.003 | 0.7 | swing, chain_not_draining | 0.42 at 10 | 33 |  | 0.45 |
| 22.16.0 | Myglevatn ndf. | 182.2 | 182.405 | 99.4 | swing, chain_not_draining | 182.42 at 153 | 462 | lake line, lake above P | 33.91 |
| 22.22.0 | Søgne | 203.0 | 204.644 | 96.5 | chain_not_draining, chain_end_open | 204.71 at 310 | 454 |  | 16.36 |
| 38.1.0 | Holmen | 116.7 | 117.488 | 99.2 | swing, downstream_unread | 117.49 at 0 | 30 |  | 28.71 |
| 41.8.0 | Hellaugvatn | 27.5 | 0.124 | 0.4 | swing, chain_not_draining | 0.13 at 91 | 110 |  | 0.13 |
| 48.5.0 | Reinsnosvatn | 120.3 | 0.001 | 0.0 | swing, downstream_unread, chain_not_draining, chain_end_open | 0.01 at -10 | 129 |  | 0.35 |
| 62.15.0 | Kinne | 510.9 | 0.195 | 0.0 | swing, chain_not_draining | 0.19 at 0 | 30 |  | 56.66 |
| 73.21.0 | Frostdalen | 25.8 | 25.521 | 97.2 | chain_not_draining | 0.01 at -14 | 30 |  | 22.11 |
| 73.27.0 | Sula | 30.4 | 0.003 | 0.0 | swing, chain_not_draining | 0.00 at 14 | 30 |  | 19.82 |
| 79.3.0 | Nessedalselv | 30.2 | 0.006 | 0.0 | swing, chain_not_draining | 30.11 at 28 | 30 | lake above P | 27.94 |
| 83.6.0 | Byttevatn | 104.5 | 0.000 | 0.0 | swing, chain_not_draining, chain_end_open | 0.06 at 24 | 31 |  | 5.71 |
| 83.12.0 | Haukedalsvatn ndf. | 205.3 | 206.980 | 99.5 | chain_not_draining | 0.03 at 0 | 30 | lake above P | 36.43 |
| 87.10.0 | Gloppenelv v/Bergheim | 218.5 | 0.092 | 0.0 | swing, downstream_unread, chain_not_draining | 0.08 at 0 | 93 |  | 0.09 |
| 88.11.0 | Strynsvatn | 485.1 | 0.677 | 0.1 | chain_not_draining, chain_end_open | 0.73 at 177 | 222 |  | 0.68 |
| 109.9.0 | Driva v/Risefoss | 744.4 | 0.000 | 0.0 | downstream_unread, chain_not_draining | 0.00 at 0 | 291 |  | 41.00 |
| 112.8.0 | Rinna | 87.9 | 88.228 | 99.2 | chain_not_draining | 88.23 at 0 | 30 |  | 31.15 |
| 122.11.0 | Eggafoss | 655.2 | 655.632 | 99.4 | chain_not_draining | 0.00 at -10 | 30 |  | 63.68 |
| 122.14.0 | Lillebudal bru | 168.1 | 168.951 | 99.3 | chain_not_draining | 0.00 at -14 | 30 |  | 49.85 |
| 122.17.0 | Hugdal bru | 545.9 | 0.013 | 0.0 | swing, chain_not_draining | 546.06 at 24 | 33 |  | 63.88 |
| 124.2.0 | Høggås bru | 494.5 | 643.919 | 99.4 | swing | 643.92 at -10 | 30 |  | 38.42 |
| 139.35.0 | Trangen | 852.3 | 0.001 | 0.0 | swing, downstream_unread, chain_not_draining | 0.00 at 14 | 87 |  | 19.26 |
| 150.1.0 | Sørra | 6.6 | 5.859 | 85.6 | chain_not_draining | 0.01 at 0 | 30 |  | 5.87 |
| 196.11.0 | Lille Rostavatn | 637.3 | 0.002 | 0.0 | swing, chain_not_draining | 0.31 at 119 | 186 |  | 0.06 |
| 200.4.0 | Skogsfjordvatn | 136.0 | 0.000 | 0.0 | swing, chain_not_draining, chain_end_open | 0.00 at 28 | 30 |  | 5.18 |
| 205.6.0 | Didnojokka | 111.0 | 111.171 | 98.4 | swing, downstream_unread, chain_not_draining | 111.17 at -14 | 186 |  | 38.23 |
| 212.10.0 | Masi | 5618.3 | 3.515 | 0.0 | swing, chain_not_draining | 3.72 at 318 | 337 |  | 3.52 |
| 234.18.0 | Polmak nye | 14171.0 | 6.279 | 0.0 | downstream_unread, chain_not_draining | 6.28 at 62 | 158 |  | 6.28 |
| 237.1.0 | Båtsfjord | 23.1 | 22.512 | 94.6 | swing | 22.51 at -20 | 30 | lake above P | 22.52 |
| 311.6.0 | Nybergsund | 4418.1 | 0.001 | 0.0 | chain_not_draining, chain_end_open | 0.32 at -10 | 73 |  | 34.25 |

## What was run

| | |
|---|---|
| Code | master `bc01cd8` (PR #182 merged), this worktree |
| Extension | built fresh from this tree into the worktree's own venv (`uv pip install --reinstall-package rasputin --no-cache ".[codecs]"`, Release, bounds checks on: `rasputin version` says `on (libc++ fast)`); its SHA-256 in `provenance.txt` |
| Machine, power | Apple M1 Max, 32 GB; on AC power, battery 100 % charged, at the start and at the end (`provenance.txt`). Nothing else ran beside the batch |
| Software | Python 3.13.15, numpy 2.5.3, shapely 2.1.2, pyproj 3.8.0, tifffile 2026.9.20 |
| Station set | `rasputin fetch-stations nve-hrd`, 2026-10-05 10:27-10:28 CEST, into `../rasputin_data/nve_hrd`: stations, NVE's polygons, ELVIS river lines and (new since the placement figures) NVE's lakes. The fetch's own `manifest.json` is copied here; the SHA-256 of every file is in `provenance.txt`. The earlier fetch (02:44, no lakes file) had identical stations, polygons and rivers (same SHA-256) |
| DEM | `../rasputin_data/DTM10_UTM33_20260925`, 254 tiles, 10 m |
| Settings | the defaults: map radius 500 m, reach up 1000 m, corridor 30 m, outline tolerance twice the cell (20 m) |
| Commands | in this order: `run.sh` (steps 1-3: the fetch, the batch, the copy of results and catchments, `analyse.py`), then `checks.py` (step 6), then `analyse.py` again (the tables in `analysis.md`, now with the checks). On 2026-10-05 the catchment copy, `checks.py` and the second `analyse.py` were run by hand after `run.sh`; they are now `run.sh`'s last lines, the same commands. `lake_query.py` (NVE's lakes along each reach, network) was added in step 4's review and feeds `checks.py`'s fifth argument, the `lake above P` column only. Then `step4.sh` (step 4: `rerun.sh` with `rerun.py` for the re-runs, then `explain.py`, `divide.py`, `chain_counts.py`, `trace_exit.py`, `probe_windows.py`, `window_check.py`, `bypass.py`, `lake_above.py` and `findings.py`; each script's docstring says what it measures). Run twice, `step4.sh`'s re-runs gave the same rows to the bit except `seconds` |

## Time and memory

Wall time 1994 s for the 140 stations in one process (`/usr/bin/time -l`,
in `batch.log`'s last lines), peak resident set 7.75 GB. Per station
(`seconds` in `results.csv`): median 4.9 s, the longest 145 s
(`234.13.0` Veahkkava, 2,077 km²). The design's estimate was tens of
minutes to a few hours for about 610 M catchment nodes; the run's
catchments hold 194.5 M (the sum of `nodes`): the refused and the wrongly
small `uncertain` catchments, among them the three largest of the 140,
made almost none.

Per-station peak memory is sampled (every 0.2 s, `rss.tsv.gz`) and
attributed to the station whose line came next (`analyse.py` says how);
it is in the "peak GB" column of `analysis.md`. Median 1.82 GB, largest
7.21 GB. The sampler can miss a peak shorter than 0.2 s, and memory the
allocator keeps from an earlier station counts again; why two small
lake stations near the Swedish border (`307.5.0`, `308.1.0`) show 7.2 GB
is not looked at.

## Files

| file | what |
|---|---|
| `results.csv` | one row per station, the batch's own table |
| `summary.json` | the batch's summary (classes, bands, tile counts, percentiles, causes) |
| `analysis.md` | every table, generated from the files here by `analyse.py`, with every station's row |
| `checks.csv` | step 6's checks per station (`checks.py`) |
| `batch.log` | the batch's stderr, a time stamp per line, and `/usr/bin/time -l` at the end |
| `rss.tsv.gz` | the memory samples |
| `catchments/` | each accepted catchment's reduced outline (GeoJSON, ours, not NVE's data), copied from the batch's `$WORK/batch/<station>.geojson` (`run.sh`; checked identical with `diff -rq`, 124 files) |
| `manifest.json`, `provenance.txt` | the station set's manifest; commit, versions, SHA-256, power |
| `reruns/` | step 4's re-runs, one CSV and log per setting |
| `step4/` | step 4's probe outputs, one file per probe (`step4.sh`) |
| `window_check_12km.csv`, `bypass.csv` | the window check and the bypass count, per station |

NVE's own data (stations, polygons, rivers, lakes, and `lake_query.py`'s
lakes along each reach) is not committed; `run.sh` fetches it, and
`provenance.txt` has its SHA-256 as fetched.
