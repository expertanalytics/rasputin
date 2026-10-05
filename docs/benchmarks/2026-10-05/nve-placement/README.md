# Placement figures, before and after (increment 29, PR 2)

`@perf`, 2026-10-05. Ola asked to see how the gauge placement works, "a few
examples with maps "before" and "after", as pngs". The design is
`docs/increments/29-nve-reference-catchments.md`, "Placement figures, before
and after (PR 2)". No refine or mesh code is touched, so no timing is
reported; the seconds in `full_runs*.json` are for orientation only.

## What was run

| | |
|---|---|
| Code | worktree `29-pr2` at `e0ccdd9`, whose code is green `193079d` (the one commit between them changes only `ROADMAP.md` and the increment file) |
| Extension | rebuilt (`cmake --build build-pyext --target _core`) and copied into the venv before the runs; both copies have the SHA-256 in `provenance.txt` |
| Machine, power | Apple M1 Max, 32 GB; on AC power, battery 100 % charged, at the start and at the end (`provenance.txt`) |
| Software | Python 3.14.7, numpy 2.5.3, shapely 2.1.2, pyproj 3.8.0, matplotlib 3.11.2 (installed into a scratch directory on `PYTHONPATH`, in no dependency list: Ola's ruling, 2026-10-04) |
| Stations, NVE polygons, river lines | `rasputin fetch-stations nve-hrd`, fetched 2026-10-05 00:44 UTC into `../rasputin_data/nve_hrd`; SHA-256 in `provenance.txt`. A second fetch 13 minutes later gave byte-identical files |
| DEM | `../rasputin_data/DTM10_UTM33_20260925` |
| Commands | `run.sh`, from the repository root; the script is `render.py`. Its last two steps (the check through the command, the `full_runs*.json` summaries) were added after the main run and run on their own, on the same inputs and work directory |

`render.py` composes the package's functions and adds nothing to the
placement: `read_stations`, `read_segments`, `read_references`,
`gauge.place`, `catchment.delineate`; for the flow paths, one window round
the reach (stage A's first window: the reach's bounds plus 30 m and 2 km),
`burn.burn_reach` on it and `_core.accumulate` on the raw and the burnt
arrays. The agreement numbers (NVE's in ours, ours in NVE's, area ratio) are
counted as the design's "Per station" section says, in `render.py`; PR 4's
`reference.py` will be the official count.

**Which placement.** The design confirms each case with `catchment
--rivers`, which places the gauge on the **nearest** river line (no
watercourse number). That is the default here and gives the committed
figures. The batch of PR 4 will use the **tiered** placement (the station's
own watercourse number first, then a prefix of it, then the river's name,
then any line); `--tiers` runs the same survey and choice that way, for
comparison (`survey_tiers.csv`, `cases_tiers.json`, `full_runs_tiers.json`,
and one figure pair, `12.70.0_tiers_*`).

**Checked against the command.** The five stations drawn with the nearest
placement were also run through `rasputin catchment --rivers`
(`cli_confirm.log`): the same segment, the same node counts and the same
verdicts as `render.py`'s runs. For every drawn station the burnt chain of
the first window equals `delineate`'s chain node for node (printed by
`draw`).

## The four cases, chosen by rule

The survey (`survey.csv`) placed and burnt all 140 stations on their first
window: 136 burnt, 3 refused on tiles of two grids (`208.3.0`, `212.49.0`,
`213.2.0`) and one with no river line within 500 m (`311.4.0`,
Femundsenden), as the design expects. Only 5 of the 136 are well posed on
the first window: the first window holds only small catchments, and a
catchment that reaches the window's edge flags its counts, so the
downstream side is "not read".

| Case | Rule | Station | Why |
|---|---|---|---|
| 1. Placed straight onto the flow path | well posed on the first window, no node lowered; smallest NVE polygon | **none qualifies** | No station lowers zero nodes: the fewest is 3 (`19.79.0`), the median 81, the most 187 (the full runs' chains have 75 to 167 nodes). The burn lowers every node that does not already fall by 1 mm, so flat stretches and lake surfaces are always lowered |
| 2. The line burnt in | well posed on the first window; largest `lowered_max_m` | **`101.1.0` Engsetvatn** | 2.52 m, the largest of the five; the full run confirms it well defined |
| 3. Marked uncertain: swing with a confluence step | among the 19 stations with another river line within 100 m of `P`, the first in list order that qualifies | **none qualifies** | All 19 were run in full (`cases.json`): 7 well defined, 9 uncertain, 3 refused (two grids: `156.24.0`, `206.3.0`, `213.4.0`). The three with `swing` (`73.27.0`, `83.6.0`, `87.10.0`) all have `chain_not_draining` too: the area jumps because the burnt chain does not carry the river's flow, not because a tributary joins. With the tiered placement none qualifies either: its nine with `swing` all have `chain_not_draining` or `downstream_unread` as well (`cases_tiers.json`) |
| 4. A lake gauge | the first lake gauge in list order that is well posed | **`2.32.0` Atnasjø** | `2.11.0` Narsjø, first in the list, is uncertain (below) |

**"A confluence step", pinned by `@perf`**: `swing` is the station's only
cause (the chain drains, its end is closed, the downstream side was read, so
the gain is real) and one step between neighbouring chain nodes is over 5 %
of the area by itself. The design does not define the term. With the looser
reading (`swing` among the causes, the step over 5 %) the first qualifying
station is `73.27.0` Sula with the nearest placement and `12.70.0` Etna with
the tiered one; both are drawn below as extras, and neither shows a
tributary.

## The figures

Each pair has the same extent, the reach's bounds plus 300 m. *Before*: a
hillshade of the raw DEM, its flow paths (nodes draining at least 0.1 km²,
darker as more drain), the chosen segment and the rest of the reach in green,
the other NVE river lines dashed, the station and the mapped position `P`.
*After*: the burnt DEM's flow paths, the burnt chain (lowered nodes marked),
the placed node, the sensitivity window (`U` up and down), and our outline
over NVE's polygon with the numbers.

| Station | Figures | Segment, distance | Placed node from `P` | Lowered | Area at the node | Swing | Verdict | NVE's in ours / ours in NVE's | Area ratio |
|---|---|---|---|---|---|---|---|---|---|
| `101.1.0` Engsetvatn (case 2) | `101.1.0_*.png` | 12155230, 3.4 m | 4.3 m | 106 nodes, max 2.52 m | 40.64 km² | 0.0 % | well defined | 98.9 % / 97.1 % | 1.018 |
| `2.32.0` Atnasjø (case 4) | `2.32.0_*.png` | 10794796 (lake line), 93.8 m | 11.6 m | 133, max 0.22 m | 459.8 km² | 0.0 % | well defined | 98.4 % / 99.2 % | 0.993 |
| `2.11.0` Narsjø (extra) | `2.11.0_*.png` | 12252329 (lake line), 83.9 m | 2.0 m | 74, max 1.08 m | 0.0012 km² | -33.3 % | uncertain: chain not draining, chain end open | 0.0 % / 100 % | 0.00001 |
| `73.27.0` Sula (extra) | `73.27.0_*.png` | 11789696, 12.0 m | 12.5 m | 13, max 0.87 m | 0.0033 km² | 63.6 % | uncertain: swing, chain not draining | 0.0 % / 0.0 % | 0.0001 |
| `24.9.0` Tingvatn (extra) | `24.9.0_*.png` | 10746896 (lake line), 164.8 m | 5.1 m | 164, max 0.22 m | 0.48 km² | -3.0 % | uncertain: chain not draining, chain end open | 0.2 % / 99.5 % | 0.0018 |
| `12.70.0` Etna (extra, tiered placement) | `12.70.0_tiers_*.png` | 12392197, 129.8 m | 28.0 m | 3, max 0.59 m | 565.5 km² | 100 % | uncertain: swing, downstream not read, chain not draining | 98.8 % / 99.3 % | 0.995 |

## What the figures show

- **Where it works** (Engsetvatn, Atnasjø): the chain follows the mapped
  line on the valley floor or the lake, the placed node sits on the DEM's
  main flow path, the area does not change within `U`, and the catchment
  agrees with NVE's to 97 to 99 %.
- **Narsjø (`2.11.0`), the increment's own example station, gets a 12-node
  catchment.** The nearest line (83.9 m) is a lake centreline of an unnamed
  river (`elvid` 002-34-2695) whose reach comes in from the north-east, not Nøra, whose lake line is 127.9 m
  away; both carry the station's watercourse number, so the tiered placement
  picks the same line. The chain then runs onto the flat lake, the end
  extension runs 496 m (the cap is 500 m) without finding lower ground, and
  the flow does not reach the placed node. The sensitivity marks it
  uncertain, so it would not be scored, but the polygon written is 0.001 km²
  against NVE's 119 km². Tingvatn (`24.9.0`, a lake gauge 165 m from its
  line) fails the same way: chain on a flat lake, end open, 0.48 km² against
  272 km².
- **Sula (`73.27.0`)**: the mapped line and the DEM's channel run visibly
  apart upstream (tens of metres, read by eye from the figure); near the gauge the chain runs one node beside
  the channel, and the placed node drains 0.0033 km². Uncertain, correctly.
- **Etna (`12.70.0`), tiered placement only**: the station's watercourse
  number (`012.EF61`) is carried by an unnamed side stream (129.8 m) and by
  another segment of Randselva's line (130.9 m); the tier rule takes the
  nearer of the two, by 1.1 m, over Randselva's own line at 21.5 m (numbered
  `012.EF51`). The mapped reach is then the side stream, ending on the main
  river; the placed node lands on the main river, so the catchment still
  agrees with NVE's (98.8 %), but the chain above it drains 0.0003 km² and
  the station is uncertain.
- **The sensitivity check caught every wrong catchment seen here.** Of the
  18 full runs with the nearest placement (the candidates of cases 3 and 4,
  not a random sample): 9 well defined, all agreeing with NVE's polygon to
  at least 87.9 % both ways; 9 uncertain, of which 5 are catchments under
  1 km² against NVE's 30 to 272 km² (`2.11.0`, `24.9.0`, `73.27.0`,
  `83.6.0`, `87.10.0`) and 4 agree to at least 98.5 % (`38.1.0`,
  `83.12.0`, `122.11.0`, `122.14.0`: uncertain by `chain_not_draining` or
  `chain_end_open`, cause not looked at).
- **Swing can be negative.** When the area at the placed node is smaller
  than upstream or larger than downstream, the one-sided swing is below 0
  (`2.11.0`: -33 %, so `catchment` prints "the area changes by -33.3 %");
  such a station is uncertain through `chain_not_draining` instead.

Not looked at: why `88.4.0` Lovatn lowers a node by 41.4 m, `2.284.0`
Sælatunga by 22.7 m and `62.10.0` Myrkdalsvatn by 20.0 m (`survey.csv`).

## Counts that differ from the design's

With the nearest placement, 19 stations have another river line within
100 m of `P`, as the design says; with the tiered placement, 26. Lake lines:
48 with either placement (the design says 46 of 139). Six stations were
refused on tiles of two grids in these runs (`208.3.0`, `212.49.0`,
`213.2.0` on the first window; `156.24.0`, `206.3.0`, `213.4.0` when their
window grew); the design expects nine, which only the batch over all 140 will
show.

## Files

| File | What |
|---|---|
| `render.py`, `run.sh` | the script and the commands |
| `provenance.txt` | commit, extension and data checksums, power state at start and end |
| `survey.csv`, `survey_tiers.csv` | every station on its first window: placement, burn, first-window sensitivity, refusal |
| `cases.json`, `cases_tiers.json`, `choose.log`, `choose_tiers.log` | the candidates run for each case, and why each was passed over |
| `full_runs.json`, `full_runs_tiers.json` | every full run: placement, burn, sensitivity, area, agreement numbers (the outline and the chain left out) |
| `cli_confirm.log` | the drawn stations through `rasputin catchment --rivers` |
| `*_before.png`, `*_after.png` | the figures |

The full-run outlines and chains (`$WORK/*/runs/*.json`, about 7 MB in all) and the
command's GeoJSON outputs stay out of the repository; `run.sh` regenerates
them (not timed; the full runs' own seconds are in `full_runs*.json`).
