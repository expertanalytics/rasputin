# Increment 39 — Six Swedish river systems: every SMHI sub-catchment meshed in one night

Status: **designed** (`@architect`, 2026-10-12), branch `worktree-se2-catchments`
off master `5d2e706d`. No production code: a batch of `rasputin mesh` runs over
SMHI's sub-catchment polygons, driven by a script under `docs/benchmarks/`, with
a table and a map per river system by the morning of 2026-10-13.

## Ola's ask

Ola, 2026-10-12: "Yes, do the swedish catchments, high-pri." The six river
systems are the Ljungan, Indalsälven, Faxälven, Fjällsjöälven, Åseleälven and
Umeälven, in Sweden's SE2 electricity bidding area. The first wording of the
task was a delineation to about 70 forecast points we do not have the
coordinates of; the main session relayed Ola's re-scope the same afternoon:
"I'll not get the list. My goal is to have each (sub)catchment TINed by
tomorrow morning." and "Should be easy peacy", then his rulings: "Yes to all,
but thousands of sub-catchments is quite extensive." So:

- every sub-catchment polygon of SMHI's SVAR 2022 for the six systems is
  meshed with `rasputin mesh`, the polygon itself as the outline, on GLO-30
  with CORINE land cover, at the vertical tolerance of the Lagan quick-check
  case (10 m), the way that case already does (`docs/benchmarks/quick/cases.toml`);
- no delineation (`rasputin catchment` is not used), so the open question of a
  catchment on a DEM in degrees (generality audit, row A1) does not arise;
- the granularity is decided by measurement (below): SMHI's finest level in
  principle, a coarser level first if the finest does not fit the night.

Nothing from the document that prompted the task is in this file beyond the
river names and the counts.

## What the data says (measured 2026-10-12, M1 Max, 10 cores, 32 GB, on AC power)

### SMHI's SVAR 2022: what the units are

SVAR is SMHI's national water archive (Svenskt Vattenarkiv). Its 2022 edition
(delivery description `SVAR 2022:1.2`, 2025-04-09,
`https://opendata-download.smhi.se/svar/SVAR_Basprodukter_2022_1_2.pdf`) has
three nested polygon levels that matter here:

- **Huvudavrinningsområden** (main catchments, numbered `HARO`): the whole
  river system to the sea. The six rivers are four of them, because
  Faxälven, Fjällsjöälven and Åseleälven are the Ångermanälven's branches:
  Umeälven `28000`, Ångermanälven `38000`, Indalsälven `40000`, Ljungan `42000`.
- **BARO** areas (`BARONR`): "polygons that divide Sweden into catchments of a
  manageable size, with a comprehensible extent and a name that describes the
  area well; a larger main catchment is divided into several". 4 to 10 per
  system here.
- **Delavrinningsområden** (sub-catchments, `ARO_UUID`): the finest level,
  "polygons representing catchments on land and at sea", generated for every
  water body, every gauge of SMHI's basic network and the stations and dams
  SMHI's hydrological model S-HYPE needs, from Lantmäteriet's 1 m height model
  averaged to a 10 m grid (downloaded October 2020) snapped to the national land
  cover grid; Norwegian parts come from NVE's REGINE units and Kartverket's
  DTM10, linked by `MAINDOWN` (each unit names its downstream neighbour;
  `-9999` = none). Fields: `AREA` (m²), `AREA_UPSTREAM` (m²), `NAME`, `HARO`,
  `BARONR`, `VAROID` (the water-body catchment the unit belongs to), `COUNTRY`.
  A unit's outlet is the point where it drains into `MAINDOWN`; SMHI publishes
  the outlets as `Utloppspunkter_2022`, not needed here.

There is a fourth level, the water-body catchments (`VAROID`), but it is nearly
the finest one (5,179 against 5,338 units), so it is not a coarser option.

**Counts, from the download** (the finest level; `BARONR` distinct values from
the same file):

| system | `HARO` | area km² (SVAR) | sub-catchments | BARO areas | sub-catchment area km²: min / median / mean / max | under 1 km² | over 100 km² |
|---|---|---|---|---|---|---|---|
| Umeälven | 28000 | 26,772 | 1,436 | 8 | 0.00 / 12.6 / 18.6 / 485 | 91 | 22 |
| Ångermanälven (with Faxälven, Fjällsjöälven, Åseleälven) | 38000 | 31,843 | 1,782 | 9 | 0.00 / 12.7 / 17.9 / 383 | 102 | 20 |
| Indalsälven | 40000 | 26,642 | 1,311 | 10 | 0.02 / 13.7 / 20.3 / 495 | 68 | 26 |
| Ljungan | 42000 | 12,839 | 809 | 4 | 0.00 / 12.1 / 15.9 / 250 | 47 | 7 |
| **all** | | **98,095** | **5,338** | **31** | | 308 | 75 |

Every unit is a single `Polygon` (two have one hole); the largest outline has
28,644 vertices. Ljungan's four BARO areas are 1,413 to 4,315 km², dissolved
from 103 to 289 units each, and each dissolves to one `Polygon`. Combined
extent in SWEREF99 TM (EPSG:3006): easting 345,540 to 763,710, northing
6,882,190 to 7,383,010; in degrees 11.93 to 20.36 E, 62.07 to 66.57 N.

Check: `python3 -c` over the four files below, or the count query
`.../Delavrinningsomraden_2022_SE2_INT/wfs?...&resultType=hits&cql_filter=HARO=42000`
(returns `numberMatched="809"`).

### Where the SVAR data is, and how it was fetched

Fetched 2026-10-12 into `$RASPUTIN_DATA/sweden_smhi_svar/se2/` with a
generic User-Agent (`Mozilla/5.0 (generic)`), paged 100 units a request:

| file | units | bytes |
|---|---|---|
| `umealven_28000_delaro_svar2022_3006.geojson` | 1,436 | 52.7 MB |
| `angermanalven_38000_delaro_svar2022_3006.geojson` | 1,782 | 61.7 MB |
| `indalsalven_40000_delaro_svar2022_3006.geojson` | 1,311 | 47.5 MB |
| `ljungan_42000_delaro_svar2022_3006.geojson` | 809 | 27.2 MB |
| `huvud_<system>_<haro>_svar2022_3006.geojson` (4 files) | 1 each | 2.6 to 3.1 MB |
| `NOTICE.txt` | | credit and licence |

169 MB in 243 s over 54 requests, by
`docs/benchmarks/2026-10-13/se2-catchments/fetch_svar.py` (also copied to
`$RASPUTIN_DATA/lands/sweden/scripts/svar_haro_fetch.py`, the country recipe's
folder). The WFS is `https://opendata-view.smhi.se/SMHI_vatten/Delavrinningsomraden_2022_<district>/wfs`,
one layer per water district: Umeälven is in `SE1_INT` (the Bothnian Bay water
district, not the electricity area), the other three in `SE2_INT`, and the
Norwegian headwater units of the Ångermanälven (45) and Indalsälven (49) and one
of the Umeälven in `NO_FI`; a per-system file joins the districts. Two quirks:
the service's chunked responses break above about 5 MB (`IncompleteRead`
after 5 to 9 MB, three times in a row at 500 units a page), so the page is
100 units; and CQL geometry literals are northing first (the Sweden recipe,
section 3).

**Licence** (checked 2026-10-12 with the same generic User-Agent,
`https://www.smhi.se/data/om-smhis-data/villkor-for-anvandning`): "Med våra
öppna data följer licensvillkoren Creative commons Erkännande 4.0 SE ... Detta
är tillåtet även för kommersiella ändamål. Du ska ange SMHI som källa och även
ange om du har ändrat i licensmaterialet." So every output states "Källa:
SMHI" and that the polygons were used as mesh outlines (dissolving BARO areas
from units is a modification, and is said so). The mesh files already carry
GLO-30's and CORINE's credits (`dem_credit`, `features_notice`).

### GLO-30 and CORINE

The four main catchments meet **36 one-degree GLO-30 tiles** (Umeälven 17,
Ångermanälven 15, Indalsälven 16, Ljungan 11; overlapping); 4 were cached
before today (`N62`–`N63`, `E012`–`E013`). `rasputin fetch glo30 --out-crs
EPSG:3006 --cache $RASPUTIN_DATA/sweden_glo30_cache` with the Ljungan box
(`--bbox 368340 6880190 628520 7005290`) fetched 12 tiles, 54 blocks,
101.9 MB in 13.6 s; with `--domain` on the main-catchment polygon,
Indalsälven fetched 54 blocks, 100.2 MB in 12.5 s, Ångermanälven 51 blocks,
99.5 MB in 13.5 s, and Umeälven 105 blocks, 203.9 MB in 23.8 s. In all, 505 MB
in about 64 s; the cache (`sweden_glo30_cache/`, 642 MB) is complete for the
four systems after those four calls, so the night's batch touches no network.

CORINE 2018 for Sweden is read straight from the European GeoPackage the
Lagan case uses (`$RASPUTIN_DATA/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg`,
EPSG:3035); the read is windowed by the domain, so a small unit costs little
(`features read` 1.8 s for a 1,489 km² polygon against 62.7 s for Lagan's
6,441 km² in the quick-check baseline). `$RASPUTIN_DATA/sweden_corine/` holds
only three old clips and is not used.

### One unit timed: today's command works on these polygons

The per-unit command is the Lagan case's, with the unit as `--domain` and a
`.gpkg` output (increment 36; QGIS opens it as a vector layer with the land
cover per triangle):

```sh
rasputin mesh --dem glo30 --cache $RASPUTIN_DATA/sweden_glo30_cache --out-crs EPSG:3006 \
  --domain <unit>.geojson \
  --features $RASPUTIN_DATA/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg \
  --features-layer U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine \
  --tolerance 10 --out <unit>.gpkg --stats <unit>_stats.md --record <unit>_record.json
```

Measured on master's code (`_core` built from master, same C++ as `5b8bb3e0`),
`/usr/bin/time -l`, one run each, outputs in
`../rasputin_scratch/sweden/se2/probe/`:

| polygon | km² | triangles | wall s | peak GB | output |
|---|---|---|---|---|---|
| Ljungan's smallest unit, `72478BAD…` (Mynnar i Hållstasjön) | 0.004 | 15 | 2.35 | 2.1 | `.gpkg` |
| Ljungan's median unit, `77F7CD61…` (Mynnar i Holmsjön) | 12.1 | 3,247 | 3.08 | 2.15 | `.gpkg` 1 run |
| the same unit to `.vtk` | 12.1 | 3,247 | 2.34 | 2.15 | `.vtk` |
| Ljungan's largest unit, `A7647994…` (Utloppet av Revsundssjön) | 250 | 31,435 | 4.22 | 2.3 | `.gpkg` |
| Ljungan BARO area 42004 (103 units dissolved) | 1,413 | 119,987 | 7.42 | 2.5 | `.gpkg` |
| Flåsjön's gauge catchment (89 units dissolved, on disk since 2026-10-06) | 1,489 | 82,217 | 7.77 | 2.6 | `.gpkg` 29 MB |

Every run met the 10 m tolerance against GLO-30's own nodes (largest
difference 9.9994 m or less; 0.62 m for the 15-triangle unit). The cost per
unit is a fixed 2.3 s (interpreter, imports, cache headers, the CORINE read)
plus about 0.004 s per km²; the 2.1 GB peak is the same fixed cost and does
not grow with the unit.

**Projection** (sum over the 5,338 units of 3.1 s plus 0.0045 s per km²
beyond 12 km²): **about 4.6 h for the finest level, one process at a time**;
**under 20 min for the 31 BARO areas** (7.4 s at 1,413 km², scaled by area to
the largest at 4,315 km²); both together under 5 h, against a window of nearly
15 h from the design (16:15 UTC) to 07:00 UTC. So the finest level fits the
night, and the BARO level is run first anyway (the main stems whole, within
20 minutes of the start). A pool of three processes (3 × 2.5 GB of 32 GB) is
expected to cut the finest level to about 2 h, because the fixed cost is
serial Python; the runner measures it rather than assumes it (first ten units
serial, then the pool).

## Prior art: legacy and literature

**Literature.** No new method: the increment orchestrates `rasputin mesh`
(increments 14, 16, 20c, 36) over published polygons. The polygons' own
method is SMHI's, stated above from the delivery description; the comparison
with reference catchments that increment 29 did for Norway
(`docs/increments/29-nve-reference-catchments.md`, "Agreement and classes")
does not apply, because nothing is delineated: the SVAR polygon *is* the
outline, so there is no second polygon to compare it with. What this run does
show is `mesh` on 5,338 real polygons at once, which is the largest batch the
command has had; the table of refusals is its finding.

**Legacy.** Nothing on Swedish data or SMHI:

```
$ git grep -n -i -l "svar\|smhi\|sweden\|sverige\|delavrinning" legacy-archive -- legacy
(no output)
```

The nearest in the current tree is `catchment_batch.py` (increment 29, PR 4:
one station at a time, each a row for a sink, refusals as rows). This
increment borrows that shape for its script but adds nothing to the package.

## The design

### What is built

No production line. Three files under `docs/benchmarks/2026-10-13/se2-catchments/`
(the `docs/` exclusion of `CLAUDE.md` §2; not counted), written by `@developer`:

1. `fetch_svar.py` — the SVAR download above (already run; committed with the
   design so the data's origin is in git).
2. `se2_batch.py` — the runner (about 150 lines; stdlib, `shapely`, `pyproj`
   only; no prohibited dependency).
3. `README.md` — the morning's record: the table, the maps, the refusals, the
   timings, the credits; written by the runner's `report` step and then edited
   by hand for the human reader.

Data flow, all under `../rasputin_scratch/sweden/se2/` (never committed):

```
sweden_smhi_svar/se2/<system>_<haro>_delaro_svar2022_3006.geojson   (SMHI, on disk)
        │ split                                   │ dissolve by BARONR
        ▼                                         ▼
units/<system>/<ARO_UUID>.geojson        baro/<system>/<BARONR>.geojson
        │                                         │
        ▼  rasputin mesh (the command above), one subprocess per polygon
out/<system>/fine/<ARO_UUID>.gpkg (+ _stats.md, _record.json)
out/<system>/baro/<BARONR>.gpkg and .vtk (+ _stats.md, _record.json)
        │
        ▼  ledger, appended one row per finished polygon
summary_fine.csv, summary_baro.csv
        │
        ▼  report
overview_<system>.svg, summary_<system>.geojson, README.md
```

### The runner, `se2_batch.py`

Subcommands, each a pure function of the files on disk:

- `prepare`: writes one GeoJSON per polygon (EPSG:3006, a `crs` member, the
  SVAR properties kept), the finest level under `units/`, the BARO level
  under `baro/` (dissolved with `unary_union(...).buffer(0)`; a dissolve that
  is not one `Polygon` keeps its largest part and says so in the ledger).
  Nothing is fetched: the SVAR files and the GLO-30 cache are inputs, and a
  missing one is a refusal naming the path.
- `run --level baro|fine [--jobs N] [--only SYSTEM]`: runs the command above
  for every polygon of the level **not already in the ledger with status
  `ok`**, largest `AREA` first across all four systems (Ola's order: the main
  stems in the morning table), `--jobs` subprocesses at once (default 1; the
  night uses 3 after the first ten units, see the projection). One row per
  polygon: `system, level, id, name, area_km2 (SVAR), triangles, vertices,
  max_error_m, seconds, peak_gb, status, message`, from the `_stats.md`'s
  `output triangles` and `output vertices` rows, the record's `max_error_m`,
  the subprocess's wall clock and `ru_maxrss`, and the exit code: `ok`, or
  `refused` with the first line of stderr (a non-zero exit from `mesh` is a
  refusal in its own words), or `error` with the exception (a crash of the
  runner itself stops the batch; a refusal does not). Resumable: a killed
  batch restarts with the same command and skips what is done; `.gpkg` is
  written through `.part` by increment 36, so a half-written file never
  carries the asked name, and a row is appended only after the subprocess
  exits. `--retry` re-runs the `refused` and `error` rows.
- `report`: the per-system summary GeoJSON (every polygon with its row's
  fields, for QGIS), an SVG map per system (outlines simplified to 100 m,
  filled by status: meshed, refused, not run; drawn with `shapely`'s own
  `svg()`, no plotting library — `matplotlib` is not installed and is not
  added), and the README tables: per system the count by status, the sum of
  triangles, the total time, the largest and the slowest unit, and the
  refusals listed with their messages.

Settings in the runner are the command's arguments above, the tolerance
(10 m), the data paths from `$RASPUTIN_DATA`, and the output root; nothing is
hard-coded to this machine beyond those two environment-derived paths.

### Which `rasputin`

The runner calls `<repo>/.venv/bin/rasputin` or `python -m tin_engine` from a
worktree whose `_core` is built from master's C++ (`tools/new_worktree.py`
does both). The measurements above used the `mapbiomas-26b` worktree's venv
with this branch's `src_python` on the path, because that venv's extension is
master's C++ (`git diff --stat master 15b77a67 -- bindings include
CMakeLists.txt` is empty); the batch does the same or builds its own — either
way the README records `rasputin version` and the commit.

### The night's order, and the morning

1. `prepare` (minutes).
2. `run --level baro` (about 20 min): 31 meshes, `.gpkg` and `.vtk`.
3. `run --level fine --jobs 3` (about 2 to 5 h): 5,338 meshes, `.gpkg` only.
4. `report`, then the README by hand; the SVGs and the two CSVs are
   committed with it (about 1 MB; the `.gpkg` files, about 5 GB in all, stay
   in scratch).
5. If time remains: one mesh per whole system from the `huvud_*` polygons,
   `.vtk` and `.gpkg`, on Lagan's scaling (6,441 km² in 112 s) about 10 min
   each; a refusal on size is accepted and recorded.

### Why not a package command

A `mesh-batch` command would be the right home once two such batches exist
(the São Francisco pieces are the other candidate); tonight's script is the
evidence for its design, not a substitute for it. Listed under Ideas.

### Tests

None in the package: no production code changes. The runner's own checks are
in its code and in the README (the counts by status add up to the counts in
the table above; every `ok` row's `.gpkg` exists; the ledger has no duplicate
id). `@reviewer` checks those three by running them.

### Speed statement

No refine or mesh code is touched; no `@perf` acceptance run. The runner's
timings go in the README as a record of `mesh` at small sizes, which the
quick check does not cover (its smallest case is a 10 m tile).

### Estimate

0 counted lines. About 150 lines of script and 40 of fetch script under
`docs/`, one README, four SVGs, two CSVs.

## Questions for Ola, with defaults

1. **Vertical tolerance.** Default: 10 m, the Lagan quick-check case's value
   (ruled "yes to all" through the main session; recorded here so the number
   is in one place). A 5 m run over the finest level would roughly double the
   triangles and add perhaps an hour; not started unless asked.
2. **Land cover.** Default: CORINE 2018 on every mesh, as the Lagan case
   (ruled the same way).
3. **Which output per unit.** Default: `.gpkg` only for the 5,338
   sub-catchments (QGIS; 3.1 s each), both `.gpkg` and `.vtk` for the 31 BARO
   areas and the whole systems (ParaView-sized). Writing `.vtk` too for every
   unit is a second `mesh` run each (2.3 s more per unit, about 3.5 h serial).
4. **Granularity.** Default: both levels, BARO first then the finest, since
   both fit the night; the whole systems after, if time remains.
5. **The GeoPackages.** Default: they stay in `../rasputin_scratch/sweden/se2/`
   (about 5 GB); say if a copy should go somewhere you can hand over.

## What each persona reads

`@developer`: this file ("The runner", "The night's order"), the Sweden
recipe (`$RASPUTIN_DATA/lands/sweden/RECIPE.md`, sections 3 and 4) and
`docs/benchmarks/quick/cases.toml` (the Lagan case). `@reviewer`: "Tests" and
the README's claims, run. No `@tester` step: there is no red suite to write.
