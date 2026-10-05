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
| refused | 16 | no catchment: 13 on tiles of two grids, 2 off the DEM's edge, 1 with no river line |

Scored stations (match, close, miss) are 85; 74 of them (87 %) match.

**Known and expected refusals, apart (not failures).**

- 13 stations refused because their windows select tiles on two different
  grids, which rasputin does not combine, and neither grid covers the
  window alone: `156.15.0`, `156.24.0`, `191.2.0`, `203.2.0`, `206.3.0`,
  `208.2.0`, `208.3.0`, `209.4.0`, `212.48.0`, `212.49.0`, `213.2.0`,
  `213.4.0`, `223.2.0`. Eight tiles of the DEM sit half a cell (5 m) east
  of the others; a later increment resamples them onto the common grid.
  Of the nine the design expected, eight are here; the ninth, `196.11.0`,
  ran (its catchment is a wrong 22-node one, below). The four the design
  named as likely (`156.24.0`, `213.4.0`, `191.2.0`, `203.2.0`) are all
  here. One was not foreseen: `212.48.0` Sagafoss, whose NVE polygon lies
  on one normal tile 16 km from any shifted one; the window's growth
  margin reached it (see "Findings").
- `244.2.0` Neiden: refused, its catchment needs DEM nodes where there is no
  tile (Finland; 131 km² of NVE's polygon lies outside every tile). Expected.
- `311.4.0` Femundsenden: no river line within 500 m. Expected until PR 5.

**Not as expected.**

- `234.18.0` Polmak nye was expected to be refused (its polygon crosses into
  Finland) but ran: the gauge's burnt river path does not carry the
  river's flow, so the flood stayed at 6.3 km² of NVE's 14,171 km² and
  never reached the border. It is `uncertain` (chain not draining).
- `2.142.0` Knappom was refused off the DEM's edge although its NVE polygon
  lies 26 km inside the tiles. Measured window by window
  (`probe_windows.py`): the catchment's nodes stayed inside x 329-369 km,
  NVE's polygon is x 324-374 km, and the fifth window's growth margin (32
  km, doubled each window) reached x 400.9 km, past the last tile. So the
  window rule refuses, not the catchment. A finding.

**By size band** (NVE's polygon area):

| band, km² | stations | match | close | miss | uncertain | refused | share uncertain (of not refused) |
|---|---:|---:|---:|---:|---:|---:|---:|
| under 10 | 12 | 5 | 3 | 0 | 4 | 0 | 33 % |
| 10-100 | 44 | 29 | 2 | 3 | 7 | 3 | 17 % |
| 100-1000 | 75 | 38 | 0 | 3 | 25 | 9 | 38 % |
| over 1000 | 9 | 2 | 0 | 0 | 3 | 4 | 60 % |

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
2. **Two tiny catchments marked well defined**: `105.1.0` Osenelv
   (0.75 km² of NVE's 138) and `153.1.0` Storvatn (2.7 km² of 49). Their
   sensitivity check passed, so they are scored `miss`, not `uncertain`.
3. **Three catchments larger than NVE's** with NVE's polygon almost wholly
   inside ours (98.5-98.8 %): `15.49.0` Halledalsvatn (ratio 1.73),
   `12.215.0` Storeskar (1.32), `307.7.0` Landbru (1.28).
4. `2.142.0` Knappom and `212.48.0` Sagafoss: refused by the window's growth
   margin, not by their catchments.

## Findings, one line each (step 4)

Filled in below as the re-runs (corridor 15 m and 60 m, map radius 250 m
and 1000 m) finish; see "Step 4 status" for what is done.

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
| every finding has its line | see "Step 4 status" |
| step 6's checks hold | the oracle and the burn's no-off-chain check hold everywhere; 35 stations whose burn did not hold are listed, as step 6 asks |
| a share of `match` | none required (Ola, 2026-10-04); this run is the baseline |

## Step 4 status

In progress at the time of this commit.

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
| Commands | `run.sh` (steps 1-3), `checks.py` (step 6), `analyse.py` (the tables in `analysis.md`), `rerun.sh` with `rerun.py` (step 4's re-runs), `probe_windows.py` (window extents) |

## Time and memory

Wall time 1994 s for the 140 stations in one process (`/usr/bin/time -l`,
in `batch.log`'s last lines), peak resident set 7.75 GB. Per station
(`seconds` in `results.csv`): median 4.9 s, the longest 145 s
(`234.13.0` Veahkkava, 2,077 km²). The design's estimate was tens of
minutes to a few hours for about 610 M catchment nodes; the run made fewer
(the refused and the wrongly small `uncertain` catchments, among them
the three largest of the 140, made almost none).

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
| `catchments/` | each accepted catchment's reduced outline (GeoJSON, ours, not NVE's data) |
| `manifest.json`, `provenance.txt` | the station set's manifest; commit, versions, SHA-256, power |
| `reruns/` | step 4's re-runs, one CSV and log per setting |

NVE's own data (stations, polygons, rivers, lakes) is not committed; `run.sh`
fetches it, and `provenance.txt` has its SHA-256 as fetched.
