# Increment 85: a catchment joined from a published sub-catchment network

**Status: designed** by `@architect`, 2026-10-21, on master `5d8904ee`,
branch `worktree-network-catchment`. Full path (about 450 counted lines, no
C++ file, no refine or mesh code, so no `@perf` run). Number 85 because 84 is
the highest number found in `docs/increments/` over every local branch and
remote-tracking ref (`git ls-tree --name-only <ref> docs/increments/`, the
merge-queue refs left out) and every worktree under `.claude/worktrees/`,
checked 2026-10-21; 82 to 84 are the private package's.

## 1. What is asked, and why

Ola said yes on 2026-10-21 to building this now, ahead of a planned MCP
tool for querying features; his words were relayed by the main session and
are not quoted here.

An outside assistant (Codex, through rasputin's MCP tools) must build the
natural catchment of each forecast point on regulated Swedish rivers. The
terrain route (`rasputin catchment`, which floods the DEM) misses whole
branches where a dam embankment in a surface model blocks the flow. The pilot
at Kultsjön (Ångermanälven, SMHI main catchment HARO 38000), run by the main
session on 2026-10-21, output in
`/Users/skavhaug/projects/rasputin_scratch/kultsjon-pilot-20261021/`:

- `rasputin catchment --dem glo30 --lakes <SMHI lakes> --seed 508443.75
  7211434.12 --seed-crs EPSG:3006`: 1,096.43 km² in 6 s.
- SMHI's own network above Kultsjön's outlet (unit
  `11B6DEC8-0A79-40FE-AF91-7413AE9B4703`): 75 units, 1,702.67 km².
- Shared 1,092.63 km²; only in SMHI's 610.04 km²; only in the terrain
  catchment 3.79 km². The 24 units of SMHI's with less than half their area
  inside the terrain catchment hold 605.9 km²: the Ransarån branch, the
  Ransaren reservoir and above. A dam embankment in GLO-30 (a surface model)
  is the likely cause; unconfirmed.
- Without `--lakes`, the point seed gave 0.0016 km² (the nearest DEM node in a
  flat lake).

All but the first figure were measured for this design (section 7.3 gives the
script's steps), with the pilot polygon moved from EPSG:4326 into EPSG:3006.

Sweden publishes the answer to "what drains here" as data: SMHI's SVAR 2022
sub-catchments, each naming the unit it drains to. Where such a network
exists, joining its units is exact to the publisher's own map, takes under a
second, and needs no DEM. This increment adds that route, and a cross-check
that runs the terrain route from the same seed and reports where the two
disagree, which would have caught Kultsjön.

## 2. Prior art: legacy and literature

*Literature.* Building a catchment by collecting the unit catchments upstream
of an outlet along a network, then dissolving them, is the standard vector
method of watershed delineation. Heberger, "Fast, accurate watershed
delineation with a hybrid of raster and vector methods" (preprint submitted to
Environmental Modelling & Software, dated 2025-12-14,
https://mghydro.com/pages/Heberger_delineation_2025.pdf, read 2026-10-21),
section 1.2, describes it: find the "terminal" or "home" unit catchment that
holds the outlet by an overlay, find all units upstream by network analysis,
merge and dissolve them. It names the method's bias: "watershed boundaries
created with vector unit catchments will usually be downstream of the
requested outlet, and the watershed will be larger than it should be". His
hybrid (section 1.3, and his `delineator.py` on MERIT-Basins) then splits the
home unit at the outlet with a raster flow-direction trace. He traces the
hybrid to Djokic and Ye (1999, "Fast Watershed Delineation") and to the
USGS's StreamStats and its Network-Linked Data Index (NLDI), whose
split-catchment process (https://api.water.usgs.gov/docs/nldi-pygeoapi/splitcatchment,
read 2026-10-21) cuts the local NHDPlus catchment at a point on a flowline,
snapping points within 200 m to the flowline. His optional `MATCH_AREAS`
moves the outlet to the reach whose upstream area best matches an estimate.

The network form this increment reads, one row per unit with its own ID and
its downstream unit's ID, is the HYPE model's GeoData form (`SUBID`,
`MAINDOWN`, as documented by the HYPEtools R package's
`AllDownstreamSubids` and import vignette, read 2026-10-21; bifurcations are
a separate file there, `BranchData.txt`). SMHI's SVAR 2022 uses the same
column name, `MAINDOWN`. HydroBASINS (Lehner and Grill 2013, "Global river
hydrography and network routing", Hydrological Processes 27(15), 2171-2186)
has the same form (`HYBAS_ID`, `NEXT_DOWN`, with `UP_AREA` in km²), as
MERIT-Basins does.

What this increment does and where it departs:

- **Whole units, the home unit only when the seed is at its outlet.** The
  vector method's downstream bias is exactly what Ola's residual-inflow work
  cannot have (a seed must be placed where it is, never moved to maximise
  area). So a point seed is accepted only when it, or the lake under it, lies
  within a radius of the home unit's published outlet point; otherwise it is
  refused with the two faithful choices in words (section 3.3). This departs
  from Heberger's whole-unit fallback and keeps the method's guarantee (the
  catchment is the publisher's own) without its bias.
- **No split of the home unit here.** The hybrid's split is the faithful
  answer for a point inside a unit, and it is the terrain route clipped to
  the home unit; it is left out of this increment (section 9) and the
  refusal says so.
- **No area matching.** `MATCH_AREAS`-style placement is ruled out by Ola's
  standing direction (position, never area).
- **The cross-check** (the terrain route from the same seed, compared unit by
  unit) is not a step of the published methods read; it is a report, not a
  new method, and no novelty is claimed for it or for anything here.

*Norway.* NVE's REGINE units do not fit this form. NVE's property list for
the catchment data (`ObjekttyperEgenskaper_Nedborfelt.pdf`, 2015-04-27,
http://gis3.nve.no/metadata/objekttyperEgenskaper/ObjekttyperEgenskaper_Nedborfelt.pdf,
read 2026-10-21) gives a REGINE unit its `vassdragsnummer`, the parent
catchment's number (`nedborfeltVassdragNrOverordnet`), the upstream area
(`nedborfeltOppstromAreal_km2`) and the names of its lowest and highest
points, but no downstream unit's ID: up- and downstream follow from the
hierarchical numbering. A REGINE network needs a reader that derives the
downstream link from the number; left out (section 9).

*Legacy.* Nothing. The command, run on 2026-10-21:

```
$ git grep -l -i -E 'maindown|next_down|downstream|subcatchment|sub-catchment|regine|svar' legacy-archive -- legacy
$ echo $?
1
```

It returned no file. The legacy tree never read a published catchment
network.

## 3. The design

### 3.1 The command and the MCP tool

`rasputin network-catchment`, served by the MCP server as the tool
`network_catchment` with no server change: the server derives one tool per
Typer command (`src_python/tin_engine/mcp/schema.py@5d8904ee:202`).

| Option | Meaning |
|---|---|
| `--out PATH` | `.geojson` or `.json`, as `catchment --out` |
| `--network PATH` | the unit polygons: `.gpkg`, `.geojson` or `.json`, read in place |
| `--network-layer TEXT` | the GeoPackage's features table |
| `--id-attribute NAME` | the unit's own ID |
| `--down-attribute NAME` | the downstream unit's ID |
| `--upstream-area-attribute NAME` | optional: the area upstream of the unit's outlet, the unit's own included |
| `--area-unit m2\|km2` | the unit of that attribute; needed with it, refused without it |
| `--name-attribute NAME` | optional: a name for the units in the report |
| `--unit ID` | seed: this unit and every unit upstream of it |
| `--seed X Y`, `--seed-crs` | seed: a point, as `catchment --seed` (default CRS EPSG:4326) |
| `--lakes PATH`, `--lakes-layer` | with `--seed`: the lake under the point is the seed, as in `catchment` |
| `--outlets PATH`, `--outlets-layer` | with `--seed`: one point per unit, keyed by the same ID attribute |
| `--outlet-radius METRES` | with `--seed`: how near the seed must be to its unit's outlet; default 200 |
| `--check-dem SOURCE` | with `--seed`: also run the terrain route and report the difference |
| `--cache`, `--threads`, `--memory-budget`, `--z-unit` | with `--check-dem`, as for `catchment` |

Exactly one of `--unit` and `--seed`. Generic: no attribute name, CRS or ID
format of SMHI's is written in the code. For SMHI the call is

```
rasputin network-catchment --out kultsjon.geojson \
  --network .../SVAR2022_delavrinningsomraden.gpkg \
  --id-attribute ARO_UUID --down-attribute MAINDOWN \
  --upstream-area-attribute AREA_UPSTREAM --area-unit m2 --name-attribute NAME \
  --seed 508443.75 7211434.12 --seed-crs EPSG:3006 \
  --lakes .../vattenytor_2016_3006.geojson \
  --outlets .../utloppspunkter_2022_3006.geojson --check-dem glo30
```

and `--unit 11B6DEC8-0A79-40FE-AF91-7413AE9B4703` in place of the four seed
options when the ID is known.

**IDs are compared as text**: each attribute value is turned into text, an
integral float written as an integer (`2120008490.0` → `2120008490`), so a
HydroBASINS integer ID, a GeoJSON number and a `--unit` argument compare
equal. A null or empty own ID is refused.

### 3.2 What is computed

1. **Read the links.** Every row's own ID, down ID and, when named, upstream
   area and name; no geometry. A GeoPackage is read with one `SELECT` of
   those columns; a GeoJSON file is parsed whole. The network's CRS is the
   file's own, read as `--lakes` reads one
   (`src_python/tin_engine/feature_input.py@5d8904ee:535-540`: a GeoJSON
   file without a `crs` member is EPSG:4326, RFC 7946).
2. **Build the network** (`Network`, pure): refuse a duplicate own ID. A down
   value that names no unit is a network end (the sea, a border, the edge of
   the file); the ends are counted by value and reported ("71 units drain
   out of the file, all to -9999" on SVAR).
3. **Find the seed unit.** `--unit`: the ID must be a unit. `--seed`: move
   the point into the network's CRS; the units whose polygon covers it
   (GeoPackage: R-tree query at the point) must be exactly one. With
   `--lakes`, the lake under the point is read as `catchment` reads it
   (`read_lake_polygons`, `src_python/tin_engine/feature_input.py@5d8904ee:543`)
   and is the seed geometry; the home unit is still the one under the point.
   Then the outlet test of section 3.3.
4. **Collect upstream.** From the seed unit, follow every unit whose down ID
   is a collected unit, breadth first. Reaching a unit twice is a cycle:
   refused, naming the units on it.
5. **Join.** Read the collected units' geometries only (GeoPackage: by
   primary key), union them (`shapely.union_all`). The result must be one
   polygon; a multipart union means units that do not touch, refused with
   the parts' areas. Holes are counted and their area reported, then left
   out of the file, as the terrain route fills them
   (`catchment_geojson` writes the exterior ring only,
   `src_python/tin_engine/io/geojson.py@5d8904ee:109-111`). No reduction: the
   outline is the publisher's, vertex for vertex.
6. **Check the area.** Areas through `area_m2`
   (`src_python/tin_engine/geodesy.py@5d8904ee:131`: planar in a projected CRS,
   geodesic on the ellipsoid in degrees). With an upstream-area attribute, the
   joined polygon's area (holes included, as the publisher's area includes
   them) must be within 1 % of the seed unit's stated upstream area;
   otherwise refused, naming both and the difference ("the network file lacks
   units, or its links are wrong"). This is what catches a clipped extract
   or a dangling link, which no structural check can tell from a network
   end.
7. **Cross-check** (`--check-dem`), section 3.5.
8. **Write** the polygon with `catchment_geojson` in the network's CRS, so
   `mesh --domain` reads it as it reads `catchment`'s output, and print the
   path. Properties: the seed (unit, or point and CRS), the home unit and its
   outlet distance, `units` (the joined IDs, sorted), `unit_count`,
   `area_km2`, `stated_upstream_area_km2` (or null), holes and their area,
   the network ends, and the check's block when run.

stderr carries one plain line per step, as `catchment` does: the network
read, the seed unit (name, own area, distance to its outlet), the join
(units, area, holes), the area check, the check's figures.

### 3.3 A point inside a unit: the whole unit only at its outlet

A point seed is accepted when the seed geometry (the point, or the lake
under it) lies within `--outlet-radius` metres of the home unit's outlet
point. Then the whole home unit drains through the seed and the catchment is
the home unit and everything upstream: position-faithful, nothing added
below the seed but at most the strip within the radius.

Otherwise refused, in words: the home unit's ID and name, its own area, the
distance from the seed to its outlet, and the two choices that stay
faithful, with their figures: `--unit <home>` (includes up to the home
unit's own area below the seed) and the terrain route (`catchment`), which
cuts at the point. Without `--outlets`, or when the outlet file has no point
for the home unit, a point seed is refused the same way ("give --outlets, or
--unit"). Two outlet points for one ID are refused.

Why not "whole unit always" (Heberger's vector fallback): it moves the seed
to the unit's outlet, which is area-maximising by construction; on Kultsjön
the bare point lies 13,988 m from its unit's outlet, inside a unit of
195.70 km². Why not the border with the downstream unit as the outlet: it
is no point; over all 26,258 linked SVAR units that share a border with
their downstream unit (11 do not) that border has a median length of
5,560 m and a median 23 % of the unit's perimeter, so "near the border" is
not "at the outlet". SMHI's outlet points lie on that border: median 5 m
from it, 99 % within 15 m (56 of 25,407 beyond 50 m, 13 beyond 500 m). 933
SVAR units have no outlet point; a point seed in one of them is refused, and
`--unit` still works.

Measured on Kultsjön: the lake under the pilot's seed is 53.40 km² and lies
85 m from the outlet point of `11B6DEC8-...`, 9 m from the border with its
downstream unit; accepted at 200 m. The bare point (no `--lakes`) lies
13,988 m away; refused.

### 3.4 Refusals

Each is a `typer.BadParameter` naming the option, exit code 2, in plain
words, as `catchment`'s are:

| Case | Option named |
|---|---|
| both or neither of `--unit`, `--seed` | `--unit` |
| `--lakes`, `--outlets`, `--outlet-radius` or `--check-dem` without `--seed`; `--area-unit` without the area attribute or the reverse; a `--*-layer` without its file | the option |
| a missing attribute in the network or outlet file | `--network` / `--outlets` |
| a null or empty own ID, a duplicate own ID | `--network` |
| `--unit` names no unit | `--unit` |
| the seed point in no unit, or on the border of two | `--seed` |
| the seed not within the radius of its unit's outlet, or no outlet point for it | `--seed` |
| a cycle upstream of the seed | `--network` |
| a multipart join | `--network` |
| the joined area more than 1 % off the stated upstream area | `--network` |

### 3.5 The cross-check against the terrain route

`--check-dem SOURCE` (with `--seed`, so the terrain route has a position)
runs `delineate` (`src_python/tin_engine/catchment.py`) with the same
point, seed CRS and lakes, exactly as `catchment` builds its request
(`src_python/tin_engine/cli.py@5d8904ee:2722-2728`), moves its reduced outline
into the network's CRS, and reports, never refuses on, the difference:

- `terrain_area_km2`, `shared_km2`, `network_only_km2`, `terrain_only_km2`;
- `units_mostly_outside`: each joined unit with less than half its area
  inside the terrain catchment (ID, name when given, area), and their count
  and total.

A failure of the terrain route itself (not cached, cut by NoData) is
reported as the check's failure in plain words, and the network catchment is
still written: the network result does not depend on the DEM. Expected on
Kultsjön (section 7.3): 1,096.43, 1,092.63, 610.04, 3.79 km²; 24 units,
605.9 km².

### 3.6 Discovery

- `sources`' `fits` (`src_python/tin_engine/source_report.py@5d8904ee:165`)
  lists `network-catchment --network` for any polygon source, as it lists
  `catchment --lakes` (the reader's acceptance condition is "polygons"), and
  `network-catchment --outlets` for any point source.
- The MCP server's instructions
  (`src_python/tin_engine/mcp/server.py@5d8904ee:89`) gain two sentences:
  where a polygon source's description says each unit names its downstream
  unit, `network_catchment` joins the units upstream of a unit or a seed at
  an outlet; `check_dem` compares it with the terrain route.
- The command's help says the same in one paragraph.

## 4. Ideas in the private file

Read `/Users/skavhaug/projects/rasputin_scratch/ideas/ideas.md` on
2026-10-21.

- **Taken up in part:** "Catchment from a lake seed behind a dike" (detect
  a lake seed cut off by a dike): the cross-check detects it wherever a
  network exists; without one, nothing changes.
- **Overlaps, left out:** "The catalogue cannot mark a surface model" (GLO-30
  is one, the likely Kultsjön cause); "Neighbouring catchments simplified
  together" (joined units share borders exactly, but nothing here
  simplifies); "Line layers in the ArcGIS fetch" (river networks as lines,
  not polygons).

## 5. Module siting and interfaces

No C++, no `_core` call, no new dependency (`shapely`, `pyproj`, `sqlite3`
only). Separation: the network logic is pure and tested without files; the
reading is I/O; the command is the thin shell.

`src_python/tin_engine/network.py` (new, pure, no I/O):

```python
class NetworkError(ValueError): ...           # plain words for the person asking

@dataclass(frozen=True, slots=True)
class Link:
    id: str; down: str; upstream_area_m2: float | None; name: str | None

@dataclass(frozen=True, slots=True)
class Network:
    links: Mapping[str, Link]                  # by own ID
    ends: Mapping[str, int]                    # down values naming no unit -> count
    @classmethod
    def build(cls, links: Sequence[Link]) -> Network: ...   # refuses empty or duplicate IDs
    def upstream(self, unit: str) -> tuple[str, ...]: ...   # unit first, breadth first; refuses a cycle

def unit_id(value: object) -> str: ...         # the text rule of section 3.1
def home_unit(point: Point, covering: Sequence[tuple[str, BaseGeometry]]) -> str: ...
def outlet_distance_m(seed: BaseGeometry, outlet: Point, crs: str) -> float: ...  # LocalFrame in degrees
def join(units: Sequence[BaseGeometry], crs: str) -> Joined: ...   # refuses multipart
def check_area(joined_m2: float, stated_m2: float) -> None: ...    # refuses beyond AREA_TOLERANCE
def compare(network: Polygon, terrain: Polygon, units: Mapping[str, BaseGeometry],
            crs: str) -> TerrainCheck: ...
```

`Joined` and `TerrainCheck` are frozen dataclasses of the figures in
sections 3.2 and 3.5.

`src_python/tin_engine/io/network_file.py` (new; reads, refuses with
`FeatureError` naming the file):

```python
def read_links(path, layer, id_attr, down_attr, area_attr, area_unit, name_attr) -> tuple[str, list[Link]]
def units_at(path, layer, id_attr, point) -> list[tuple[str, BaseGeometry]]   # point in the file's CRS
def unit_geometries(path, layer, id_attr, ids) -> dict[str, BaseGeometry]
def outlet_point(path, layer, id_attr, unit, crs) -> Point | None              # in `crs`
```

The command `network_catchment` in `src_python/tin_engine/cli.py`, beside
`catchment`, reusing `cli_kit.destination`, `cli_kit.cli_resources`,
`_dem_repository` and `catchment_geojson`. `project_structure.md` gains the
two modules.

Async-readiness: every step is blocking: under 1 s on Kultsjön, 54 s at SVAR's largest set (section 7.3); the MCP
server runs the command as a child process, as every tool.

## 6. Constants

| Constant | Value | Scale assumed | Checked at |
|---|---|---|---|
| `--outlet-radius` default | 200 m | the publisher's outlet points lie on the unit's downstream border; a chosen value, the NLDI's flowline snap distance | SVAR: 99 % of 25,407 outlet points within 15 m of that border; Kultsjön lake at 85 m |
| `AREA_TOLERANCE` | 1 % of the stated upstream area | stated areas rounded per unit | SVAR, all 26,340 units: the sum of `AREA` upstream is at most 7.1e-7 off `AREA_UPSTREAM`; Kultsjön's union area equals it to 1e-6 km². HydroBASINS (km², one decimal) not checked |
| `MOSTLY_OUTSIDE` | under 1/2 of a unit's area inside the terrain catchment | a chosen threshold for the report only | Kultsjön: 24 units, 605.9 km² |

The radius and the tolerance are in metres and in a fraction; on a network
in degrees the radius is measured in a `LocalFrame` at the outlet and areas
are geodesic, so no constant assumes a projected CRS. Largest network
checked: SVAR 2022, 26,340 units.

## 7. Tests

### 7.1 Invariant-critical suite (mutation testing required)

`tests/python/test_network.py`, pure, on hand-built `Link` lists and
shapely squares of 1 km side in EPSG:3006 (each unit 1.000000 km²) unless
stated. Mutation targets the kill record must cover: the upstream walk
(direction, the visited test, the seed included), the cycle refusal, the
duplicate-ID refusal, the network-end count, the outlet radius comparison
(`<=` against `<`), the area tolerance comparison, the multipart refusal,
`unit_id`'s integral-float rule.

1. Chain C→B→A (A at the end, down `-9999`): `upstream("A") == ("A","B","C")`,
   `upstream("B") == ("B","C")`; ends `{"-9999": 1}`.
2. Confluence: B and C both drain to A, D to C: `upstream("A")` holds all
   four, A first.
3. A cycle B→C→B upstream of A: refused naming B and C; `upstream` of a unit
   outside the cycle still works.
4. Duplicate own ID: refused naming it. Empty and null IDs: refused.
5. `unit_id`: `7`, `7.0`, `"7"` give `"7"`; `7.5` gives `"7.5"`.
6. Join of a 2 × 1 row of squares: one polygon, 2.000000 km², no holes.
   Join of a ring of eight squares round a missing ninth: one polygon,
   8.000000 km², one hole of 1.000000 km². Two squares 1 m apart: refused,
   naming the two parts' areas.
7. Area check: stated 2.000000 km², joined 2.000000: passes; joined 1.979
   (1.05 % short): refused; joined 1.981 (0.95 %): passes.
8. Home unit: a point inside one square; a point in no square (refused); a
   point on the shared edge of two (refused, both named).
9. Outlet distance: 199.9 m and 200.0 m accepted at radius 200, 200.1 m
   refused; a lake polygon whose edge is 85 m from the outlet accepted.
   Degrees: the same squares moved into EPSG:4326 give the same distance
   within 0.1 % and the same areas within 0.1 %.
10. Compare: network = two squares, terrain = one of them plus a 0.1 km²
   strip outside: shared 1.0, network-only 1.0, terrain-only 0.1, one unit
   mostly outside.

### 7.2 Command tests

`tests/python/test_cli_network_catchment.py`, on a GeoPackage and a GeoJSON
network written by the test (the five-unit network of test 2 plus an end),
an outlet GeoJSON and a lake GeoJSON. The happy path writes a GeoJSON that
`mesh --domain`'s domain reader accepts in the network's CRS, with `units`,
`unit_count`, `area_km2` and `stated_upstream_area_km2` as computed; the same
network as `.gpkg` and `.geojson` gives the same file. One test per row of
the section 3.4 table. `--check-dem` on the existing small DEM fixture of the
`catchment` tests, with a network of two squares drawn round its known
catchment, asserts the block's keys and that the network is written when the
terrain route is refused (DEM not cached). The `fits` test gains the two new
fits; the MCP schema snapshot gains the tool.

### 7.3 The Kultsjön case, recorded

Not a pytest: the test session hides the data folder
(`tests/python/conftest.py`), and the SVAR file is 830 MB. After green, the
main session runs the section 3.1 command on the registered SMHI files and
GLO-30, and records in section 11 that it gives exactly the figures below,
measured for this design on 2026-10-21 with a script reading the same files
(links read 0.23 s, upstream walk 0.01 s, union 0.33 s):

| Figure | Value |
|---|---|
| units joined | 75 |
| joined area, stated upstream area | 1,702.67275 km², 1,702.67275 km² |
| parts, holes, vertices | 1, 0, 19,556 |
| home unit, own area | `11B6DEC8-0A79-40FE-AF91-7413AE9B4703` "Utloppet av Kultsjön,Gåaltoe", 195.70 km² |
| lake under the seed, distance to the outlet | 53.40 km², 85 m: accepted |
| bare point (no `--lakes`), distance | 13,988 m: refused |
| check: terrain, shared, network only, terrain only | 1,096.43, 1,092.63, 610.04, 3.79 km² |
| units mostly outside the terrain catchment | 24, 605.9 km² |
| SVAR network ends | 71 units, all down value `-9999` |

The largest upstream set in SVAR, at the unit with the largest stated
upstream area ("N Kvarkens S kustvatten"): 8,453 units, 168,699.1 km² joined
and stated, one polygon of 198,843 vertices; links read 0.30 s, union
53.7 s with `shapely.union_all` (shapely 2.2.0, GEOS 3.14.1). The faster
`shapely.coverage_union_all` is not usable: SVAR's neighbours do not share
vertices exactly, and it raised "CoverageUnion cannot process incorrectly
noded inputs" on that set and returned 122 parts for Kultsjön's 75 units.
The join uses `union_all`; a faster join at basin scale is left out
(section 9).

## 8. Size

About 450 counted lines, worst 580, one PR under the 700-line limit:
`network.py` about 150, `io/network_file.py` about 70, the command about 210
(sixteen options at about seven lines each, and the body), `fits` and the MCP
instructions about 10, `project_structure.md` not counted. If the developer's
count passes 650, the cross-check (section 3.5, about 90 lines with its four
options) moves to a second PR on the same increment.

## 9. Left out

- **Splitting the home unit at a point** (the hybrid method's step): the
  terrain route clipped to the home unit. The faithful answer for a point
  inside a unit; a later increment, since the refusal of section 3.3 is safe
  meanwhile.
- **REGINE**: a reader deriving the downstream link from NVE's hierarchical
  numbers.
- **Bifurcations** (a unit draining to two units, HYPE's `BranchData.txt`):
  the network has one down ID per unit.
- **A network in the catalogue** with its attribute names, so `--network`
  takes a source name: vector options take paths today.
- **Reducing the joined outline**: the outline is written as published
  (19,556 vertices on Kultsjön); `mesh --domain` meshed SVAR polygons
  directly in increment 39.
- **A faster join at basin scale**: 54 s for SVAR's largest set of 8,453
  units (section 7.3); a pairwise or cascaded union, if a run of that size
  is wanted interactively.

## 10. Decided on defaults (no question for Ola)

- The command is a new one, not a mode of `catchment`: its inputs (a network
  and attribute names) and its result (no DEM, no reduction) differ, and a
  mode would double `catchment`'s refusal table.
- The output CRS is the network's, as asked.
- A point seed away from its unit's outlet is refused, not cut and not
  widened (section 3.3).
- The cross-check reports and never refuses.
