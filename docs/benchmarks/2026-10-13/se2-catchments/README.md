# SE2 catchments: six Swedish river systems meshed sub-catchment by sub-catchment

Increment 39 (`docs/increments/39-se2-catchments.md`). Every SMHI SVAR 2022
sub-catchment of the Umeälven, Ångermanälven (with Faxälven, Fjällsjöälven and
Åseleälven), Indalsälven and Ljungan, meshed with `rasputin mesh` on GLO-30
with CORINE 2018 land cover at a 10 m vertical tolerance, the polygon itself as
the outline. Two levels: the 31 BARO areas first, then the 5,338 finest units;
the four whole systems after, if the night allows.

Källa: SMHI. The SVAR 2022 polygons (CC BY 4.0 SE) are used as mesh outlines;
the BARO areas are dissolved from the units (a modification). The meshes carry
GLO-30's and CORINE's credits in their `run` layer (`dem_credit`,
`features_notice`).

## Files here

| file | what |
|---|---|
| `fetch_svar.py` | the SVAR download (committed with the design; already run) |
| `se2_batch.py` | the runner: `prepare`, `run --level baro|fine|huvud`, `report` |
| `summary_baro.csv`, `summary_fine.csv`, `summary_huvud.csv` | the ledgers, one row per polygon (copied from scratch by `report --copy-to`) |
| `overview_<system>_<level>.svg` | one map per system and level, filled by status |
| `report.md` | the runner's own tables, as written by `report` |

The meshes themselves (`.gpkg` per unit, `.gpkg` and `.vtk` per BARO area)
are under `../rasputin_scratch/sweden/se2/out/<system>/<level>/`, never
committed. Regenerate: `prepare`, then `run --level baro` and
`run --level fine --jobs 3` from a worktree with a `.venv` (the runner finds
`<worktree>/.venv/bin/rasputin`; `--rasputin` names another).

## Method

- **Code**: branch `worktree-se2-catchments` at `0355e898` (code-identical to master `f09b9863`
  plus the design commit; `src_python`, `bindings`, `include` and
  `CMakeLists.txt` unchanged from master). `rasputin version` 0.2.0.dev0;
  `_core` built from this tree by `tools/new_worktree.py --existing` on
  2026-10-12 (uv wheel build, Release, bounds checks on, libc++ fast).
- **Machine**: Apple M1 Max, 10 cores, 32 GB, macOS 27.0; AC power
  throughout (`pmset -g batt`: "Now drawing from 'AC Power'"). Nothing else of
  ours ran beside the batch; the one-minute load before the BARO run was 4.8
  from the desktop's own processes (iTerm, loginwindow, Spotlight), no
  `rasputin` process.
- **Python**: 3.14.7; shapely 2.2.0, pyproj 3.8.0, numpy 2.5.4. No plotting
  library: the maps are SVG from shapely's `svg()`.
- **Inputs**: `$RASPUTIN_DATA/sweden_smhi_svar/se2/` (5,338 units, 31 BARO
  numbers, 4 main catchments; `fetch_svar.py`), the GLO-30 cache
  `$RASPUTIN_DATA/sweden_glo30_cache` (642 MB, complete for the four systems;
  no network during the batch), CORINE 2018 from the European GeoPackage the
  Lagan quick-check case reads.
- **The command per polygon**, as `docs/benchmarks/quick/cases.toml`'s Lagan
  case: `rasputin mesh --dem glo30 --cache ... --out-crs EPSG:3006 --domain
  <polygon>.geojson --features <CORINE.gpkg> --features-layer
  U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine
  --tolerance 10 --out <id>.gpkg --stats <id>_stats.md --record
  <id>_record.json`; BARO areas and whole systems get a second run to `.vtk`.
- **Timing**: wall clock of the subprocess from the runner (`seconds`;
  interpreter start-up and imports included, unlike the `--stats` total), one
  run each, no repeats; `peak_gb` is the child's own `ru_maxrss` from
  `wait4`. BARO areas ran one at a time; the finest level three at a time, so
  its per-unit seconds are under contention and are not comparable with the
  BARO rows or with the design's probe.
- **Quality**: `max_error_m` from the run record (the largest difference
  against GLO-30's own nodes); a row over 10 m is counted in the table's last
  column. The `--stats` file beside each mesh has the angle and degree
  tables.

## Preparation (2026-10-12 18:22 CEST)

`prepare`: 38.9 s. 5,338 unit files, 31 BARO files, 4 whole-system files; every
BARO dissolve is one `Polygon` (no `note` in `index_baro.csv`). Check:
`wc -l index_*.csv` gives 5,339, 32 and 5 lines.

## BARO level (31 areas, one process at a time)

`run --level baro --jobs 1`, 2026-10-12 16:24:55 to 16:36:30 UTC, 694 s in
all: **31 of 31 meshed, 98,095 km², 6,149,260 triangles, 366.4 s for the
`.gpkg` runs and 327.9 s for the `.vtk` runs, no refusal, no row over 10 m**
(largest `max_error_m` 9.999995 m), peak 3.26 GB (the largest, 38009). Three
outputs opened after the run: `sqlite3 <id>.gpkg "select count(*) from
triangles"` gives the row's triangle count for 42002, 28008 and 40004, their
`run` layers carry `max_error_m` and the tiles, and the `.vtk` headers'
`POLYGONS` counts match. The ledger is `summary_baro.csv` (copied here by
`report --copy-to`; the original is `../rasputin_scratch/sweden/se2/summary_baro.csv`).

The `_stats.md` and `_record.json` beside each BARO mesh are the `.vtk` run's
(the second run wrote over the first's; the ledger's figures were read from the
first). The runner tags the second run's files `_vtk_stats.md` from the fix in
the next commit on, so the whole-system runs keep both.

| system | BARO | units | km² | triangles | vertices | max error m | gpkg s | vtk s | peak GB |
|---|---|---|---|---|---|---|---|---|---|
| angermanalven | 38001 | 206 | 3,543 | 212,431 | 131,115 | 10.000 | 13.54 | 11.90 | 2.85 |
| angermanalven | 38002 | 173 | 2,942 | 184,903 | 113,447 | 10.000 | 11.34 | 10.27 | 2.90 |
| angermanalven | 38003 | 190 | 3,546 | 261,755 | 164,364 | 10.000 | 17.07 | 15.48 | 3.24 |
| angermanalven | 38004 | 50 | 1,077 | 74,834 | 49,120 | 10.000 | 6.07 | 5.62 | 2.52 |
| angermanalven | 38005 | 276 | 5,332 | 287,821 | 174,227 | 10.000 | 17.54 | 15.91 | 3.24 |
| angermanalven | 38006 | 51 | 869 | 68,505 | 43,538 | 10.000 | 4.99 | 4.64 | 2.46 |
| angermanalven | 38007 | 365 | 5,508 | 344,406 | 207,836 | 10.000 | 20.41 | 17.43 | 3.07 |
| angermanalven | 38008 | 153 | 2,477 | 176,907 | 107,469 | 10.000 | 10.56 | 9.56 | 2.87 |
| angermanalven | 38009 | 318 | 6,549 | 400,281 | 229,447 | 10.000 | 19.05 | 16.70 | 3.26 |
| indalsalven | 40001 | 119 | 2,972 | 178,131 | 105,447 | 10.000 | 9.03 | 7.95 | 2.87 |
| indalsalven | 40002 | 57 | 922 | 77,922 | 47,942 | 10.000 | 5.20 | 4.74 | 2.35 |
| indalsalven | 40003 | 160 | 3,083 | 195,111 | 119,272 | 10.000 | 11.51 | 10.23 | 2.67 |
| indalsalven | 40004 | 217 | 4,329 | 327,921 | 201,986 | 10.000 | 19.45 | 17.44 | 3.09 |
| indalsalven | 40005 | 70 | 1,324 | 97,852 | 62,624 | 10.000 | 6.72 | 6.15 | 2.60 |
| indalsalven | 40006 | 129 | 2,297 | 146,254 | 88,584 | 10.000 | 8.51 | 7.65 | 2.87 |
| indalsalven | 40007 | 223 | 4,127 | 248,873 | 155,807 | 10.000 | 14.17 | 12.60 | 3.04 |
| indalsalven | 40008 | 194 | 3,998 | 212,926 | 124,524 | 10.000 | 11.45 | 10.19 | 2.89 |
| indalsalven | 40009 | 58 | 1,037 | 68,478 | 42,837 | 10.000 | 5.21 | 4.82 | 2.58 |
| indalsalven | 40010 | 84 | 2,553 | 124,420 | 77,942 | 10.000 | 8.05 | 7.28 | 2.78 |
| ljungan | 42001 | 213 | 3,384 | 243,936 | 145,436 | 10.000 | 12.18 | 10.73 | 2.81 |
| ljungan | 42002 | 289 | 4,315 | 307,820 | 180,287 | 10.000 | 13.99 | 12.11 | 2.88 |
| ljungan | 42003 | 204 | 3,727 | 185,704 | 113,212 | 10.000 | 10.96 | 9.69 | 2.91 |
| ljungan | 42004 | 103 | 1,413 | 119,987 | 72,290 | 10.000 | 6.86 | 6.15 | 2.55 |
| umealven | 28001 | 169 | 2,678 | 192,650 | 125,940 | 10.000 | 14.09 | 12.90 | 2.92 |
| umealven | 28002 | 69 | 839 | 56,711 | 37,986 | 10.000 | 5.27 | 4.93 | 2.38 |
| umealven | 28003 | 180 | 3,576 | 214,514 | 137,460 | 10.000 | 14.92 | 13.61 | 2.99 |
| umealven | 28004 | 107 | 2,973 | 154,652 | 97,341 | 10.000 | 9.62 | 8.91 | 2.86 |
| umealven | 28005 | 284 | 4,423 | 256,856 | 159,041 | 10.000 | 16.03 | 14.40 | 2.97 |
| umealven | 28006 | 215 | 4,391 | 256,822 | 162,503 | 10.000 | 16.81 | 15.26 | 3.13 |
| umealven | 28007 | 58 | 1,083 | 74,001 | 46,978 | 10.000 | 5.06 | 4.64 | 2.43 |
| umealven | 28008 | 354 | 6,809 | 395,876 | 235,271 | 10.000 | 20.78 | 18.01 | 3.11 |

"10.000" is the three-decimal rounding of values between 9.9997 and 9.99999 m;
the exact values are in the CSV. Against the design's probe (42004: 7.42 s,
119,987 triangles, 2.5 GB with `/usr/bin/time`), this run gives 42004 the same
triangle count, 6.86 s and 2.55 GB. Seconds grow with area roughly as
2.3 s + 0.0027 s/km² at this level (6,809 km² in 20.8 s), under the design's
0.0045 s/km² projection; the 20 minutes projected became 11.6.

## Finest level: the launch record

Launched detached at 2026-10-12 16:37:16 UTC (18:37 CEST) by
`../rasputin_scratch/sweden/se2/launch_fine.sh` (`nohup`, pid in `fine.pid`:
79982, the `sh` that runs the two steps), log `../rasputin_scratch/sweden/se2/fine.log`:

1. `run --level fine --jobs 3`: the 5,338 units, largest first, three
   subprocesses at once;
2. then `run --level huvud --jobs 1`: the four whole systems, `.gpkg` and
   `.vtk`, one at a time (design step 5; a refusal on size is a row, not a
   stop).

At 16:37:57 UTC, 41 s in: 25 units landed, all `ok`, 25 `.gpkg` on disk,
the largest units (190 to 500 km²) at 3.0 to 7.6 s each with three at once (23 of the first 25 within 3.7 to 5.4 s)
(the design's probe: 4.2 s for a 250 km² unit alone), 1.6 s per unit of wall
clock at the top of the size range. **Projection** (not a measurement): the
median unit is 12 km² and 3.1 s alone in the probe, so about 1.1 s per unit
wall with three at once, giving roughly 1.6 to 2 h for the level, a finish
around 18:30 UTC (20:30 CEST), and the four whole systems (26,642 to 31,843
km²; Lagan's 6,441 km² took 112 s in the quick check) perhaps 40 min more.
Resumable: if the machine slept or the process died,
`se2_batch.py run --level fine --jobs 3` from this worktree skips the `ok`
rows and continues; `ps -p $(cat ../rasputin_scratch/sweden/se2/fine.pid)`
says whether it is still running.

The per-unit `seconds` at this level are under contention (three processes,
each with 10 refine threads) and are not comparable with the BARO rows.

## The morning

From this worktree (the main session):

```sh
tail -3 ../rasputin_scratch/sweden/se2/fine.log          # "huvud exit 0" when both steps are done
.venv/bin/python docs/benchmarks/2026-10-13/se2-catchments/se2_batch.py report \
    --copy-to docs/benchmarks/2026-10-13/se2-catchments
```

`report` writes `report.md` (the per-system table by status, triangles, mesh
seconds, the largest and the slowest polygon, the rows over 10 m, and every
refusal with its message), one `overview_<system>_<level>.svg` per system and
level, and `summary_<system>_<level>.geojson` (every polygon with its row's
fields, for QGIS) under the scratch root, and copies the SVGs, the ledgers and
`report.md` here. The three checks the design names for `@reviewer`: the
counts by status add up to 1,436 / 1,782 / 1,311 / 809 (`report`'s table);
every `ok` row's `.gpkg` exists (`ls out/*/fine/*.gpkg | wc -l` equals
`grep -c ,ok, summary_fine.csv`); the ledger has no duplicate id
(`cut -d, -f3 summary_fine.csv | sort | uniq -d` is empty).

## Result (the batch finished 2026-10-12 18:05:23 UTC)

| level | polygons | ok | refused | triangles | wall time | largest max error |
|---|---:|---:|---:|---:|---:|---:|
| finest SVAR units | 5,338 | 5,338 | 0 | 17,889,867 | 4,865 s, three at once | under 10 m |
| BARO areas | 31 | 31 | 0 | 6,149,260 | 694 s, one at a time | 9.999995 m |
| whole systems | 4 | 4 | 0 | 4,914,027 | 421 s | under 10 m |

Every finest unit's `.gpkg` is in `../rasputin_scratch/sweden/se2/` (5,338
files, the largest error 9.999996 m); the BARO areas and the whole systems
have a `.vtk` beside the `.gpkg`. The per-system tables are in `report.md`
(its "mesh s" column adds the per-unit seconds, which ran three at a time,
so it is process time, not wall time), the ledgers in the three
`summary_*.csv`, and one overview map per system and level in the twelve
`overview_<system>_<level>.svg`, each unit filled by its status.
