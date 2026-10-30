# Increment 114: every --bbox names its coordinate system

**Status: ready for `@developer`** (red suite in; "Rulings on the red suite" settles its pins). Designed by `@architect`, 2026-10-29, on
master `f0dfe61b`, worktree `bbox-crs`; design review round 1 applied on
master `fbeaf649` (no second design round: its fixes are checked at code
review). One PR, full loop, Python only. Number 114: 113 is the highest
under `docs/increments/` over every worktree (`lsm-thin-soil`), 112 the
highest on the private repo's refs. Every `file:line` below is at
`fbeaf649`. The branch `worktree-prefer-holes` may shift these line numbers.

## 0. The plain picture

Today a box is four bare numbers, and each command guesses their CRS:
`sources` assumes degrees, `mesh` and `fetch` assume `--out-crs`, or else the
DEM's or the source's CRS. After this change a box carries its CRS. `--bbox
XMIN YMIN XMAX YMAX` needs `--bbox-crs CRS`, and a run without it is refused.
A box in a CRS other than the one the command works in is moved into that
CRS first.

    rasputin sources --bbox 10 59 11 60 --bbox-crs EPSG:4326
    rasputin mesh --dem dtm10 --bbox 10 59 11 60 --bbox-crs EPSG:4326 --out a.vtk

## 1. Prior art: legacy and literature

*Legacy.* `git grep -n -i bbox legacy-archive -- legacy` returns only
lines 202-203 of the archive's gml_repository.py, a domain polygon's bounds used for a buffer
size. That box is internal, and its CRS is the domain's. Nothing to carry.

*Literature.* Two published forms exist, and both fall back to WGS 84
longitude and latitude when no CRS is given:
- OGC API - Features (Part 2, *CRS by reference*) has a companion `bbox-crs`
  parameter. Without it the values are read as CRS84. See the OGC API Maps
  discussion quoting the rule (github.com/opengeospatial/ogcapi-maps/issues/82)
  and pygeoapi's docs, which move the box into the storage CRS
  (docs.pygeoapi.io/en/latest/crs.html).
- WFS 1.1 and 2.0 key-value requests put a CRS URI as an optional fifth value
  of `BBOX`. Without it the values are WGS 84 degrees (OGC 04-094, quoted at
  discourse.osgeo.org, "WFS 1.1.0 and BBOX with KVP").

We take OGC API's form and name. There are two departures, each with its
reason:
1. **No default.** A missing CRS is refused. Ola's rule is that every
   coordinate has a CRS and none is implied. This drops the convenience of a
   default and nothing else.
2. **Axis order is always x then y.** x is easting or longitude whatever the
   CRS's own axis order. This is the `always_xy` convention of
   `src_python/tin_engine/crs.py:3-6`. Under WFS 1.1, `urn:...:EPSG::4326`
   means latitude first. Here `EPSG:4326` means longitude first, and the help
   text says so.

The fifth-value form (`--bbox X0 Y0 X1 Y1 CRS`) was probed and rejected.
`114-probes/click_probe.txt`: `mesh --bbox 1 2 3 4 --out a.vtk` swallows
`--out` as the CRS and fails with "Missing option '--out'", which points at
the wrong thing. The companion option gets refusal words of our own. It also matches the
repository's `--domain-crs`, `--seed-crs` and `--tolerance-near-crs`. No
novelty is claimed.

## 2. Every box today (master `fbeaf649`)

| Where | Assumed CRS today | Stated at |
|---|---|---|
| `mesh --bbox` (MCP `mesh.bbox`) | `--out-crs` if given, else the DEM's | `cli.py:596-604`, `cli.py:725` (`--out-crs` help); box built at `cli.py:2930-2936`; applied in `dem_input.py:214` (projected path) and `dem_input.py:344-346` (target CRS) |
| `fetch --bbox` (MCP `fetch.bbox`) | `--out-crs` if given, else the source's | `cli.py:2634-2640`, box built at `cli.py:2696`; `fetch/plan.py:38-40`, `_frame_box` `fetch/plan.py:114-122`, `box_in_crs` `fetch/plan.py:134-141`, `region_sha256` `fetch/plan.py:163-172` |
| `sources --bbox` (MCP `sources.bbox`) | WGS 84 degrees | `cli.py:2782-2790`, check `cli.py:2802-2803`; `source_report.py:373-377`, `_cover` `source_report.py:226-234`; MCP `INSTRUCTIONS`, `mcp/server.py:93`; `README.md:213` |
| hint `_fetch_run` (mesh refusal: "run: rasputin fetch ...") | mesh's own rule, passed on | `cli.py:1635-1649`, used at `cli.py:998-1038`, `1245` |
| hint `_bbox_fetch` (window box of a `NotCached`) | the source's, unstated | `cli.py:1652-1659`; box set at `io/repository.py:381`; `NotCached` `io/cog.py:56-75`; used at `cli.py:1244`, `3738` (catchment) |

No other command or MCP argument takes a box. `git grep -n -- '--bbox'
src_python` lists only these, plus the wording in `mosaic.py:221,614`. The MCP
tools are derived from the click parameters (`mcp/schema.py:76-93`), so
`mesh`, `fetch` and `sources` are the three MCP tools with a `bbox`.

## 3. Design

**Type.** `CrsBox` goes in `domain.py`, beside `DomainPolygon`, the other
geometry with a CRS. It is a frozen Pydantic model:

    class CrsBox(BaseModel):
        crs: str                                    # validated by parse_crs
        box: tuple[float, float, float, float]      # x_min, y_min, x_max, y_max
        def outline(self, dst: str | CRS) -> Xy     # (4 * (EDGE_STEPS + 1), 2)
        def to_crs(self, dst: str | CRS) -> tuple[float, float, float, float]

- Invariants: every value is finite; `x_min <= x_max` and `y_min <= y_max`
  (words: "XMIN must not exceed XMAX, nor YMIN YMAX"). Equal values are
  allowed because `sources` takes a point. In a geographic CRS,
  `-90 <= y_min` and `y_max <= 90`, else "latitude must lie within ±90° in
  <crs>". `mesh` and `fetch` still refuse an empty box through `Bounds`,
  after the move.
- `outline(dst)`: each side of the box sampled at `EDGE_STEPS + 1` equal
  steps, mapped point by point by `crs.reprojector(crs, dst)`. Built with
  NumPy, not `shapely.segmentize`: GEOS refuses a point box
  (`114-probes/box_probe.py` hit "point array must contain 0 or >1
  elements").
- `to_crs(dst)` returns the numbers untouched (bit for bit) when
  `same_crs(crs, dst)`. Otherwise it returns the outline's bounds. Refused,
  as `ValueError` naming both CRSs:
  - a non-finite point: "the box (...) in A has no image in B";
  - with `dst` geographic, two neighbouring outline points more than 180°
    of longitude apart: "the box (...) in A crosses ±180° or holds a pole
    in B; split it". On the probe, a box across 180° (EPSG:32660) and one
    round the pole (EPSG:3413) step 360.0°, the quad 1e-6° and the São
    Francisco box 0.0064° (`114-probes/box_probe.txt`, `edge_probe.txt`).
- **What the move gives.** The outline's points lie on the true image of
  the box's edges, so its bounds lie inside or on the true image's bounds,
  never outside. How far inside, against 200 001 points an edge
  (`114-probes/edge_probe.txt`), on a 12° x 14° box (São Francisco basin,
  `(-48, -21, -36, -7)`, to and from EPSG:32723):

  | Move | degrees to UTM | UTM to degrees |
  |---|---|---|
  | pyproj `transform_bounds`, 21 points an edge (today's) | 8.7 m | 26.4 m |
  | `outline`, `EDGE_STEPS = 2000` | 0.0 m | 1.6 mm |

  The 0.0 m is luck: the edge's extreme, on the zone's central meridian,
  falls on a sample. The error falls with the square of the step. The mesh then snaps the box
  outward to whole nodes (30 m on São Francisco's DEM), so a node line is
  lost only when the true edge lies within 1.6 mm past it.
- **Not used:** pyproj's `transform_bounds` at more points. Its cost grows
  with the square of the count (0.06, 1.3, 118 ms at 21, 101, 1001 points),
  where `outline` takes 1.1 ms at 2000.

**CLI.** `mesh`, `fetch` and `sources` each gain the option
`--bbox-crs TEXT`, help: "The --bbox's CRS, anything pyproj reads; EPSG:4326
is longitude and latitude in degrees. Never assumed. Required with --bbox."
One helper in `cli.py` turns `(bbox, bbox_crs)` into a `CrsBox | None` for
all three, replacing the `Bounds` built at `cli.py:2696` and
`cli.py:2930-2936` and the check at `cli.py:2802-2803`. It refuses:

| Case | Exit | Words (param `--bbox-crs`) |
|---|---|---|
| `--bbox` without `--bbox-crs` | 2 | "--bbox needs --bbox-crs, the CRS its four numbers are in (EPSG:4326 for longitude and latitude in degrees); a box's CRS is never assumed" |
| `--bbox-crs` without `--bbox` | 2 | "applies only with --bbox" (as `--domain-crs`, `cli.py:1078-1079`) |
| CRS pyproj cannot read | 2 | `parse_crs`'s own words |
| a broken invariant (param `--bbox`) | 2 | `CrsBox`'s words |
| no image, or across ±180° | 2 | `to_crs`'s words |

The help of each `--bbox` drops its assumption and says "XMIN YMIN XMAX YMAX
in --bbox-crs (required with it); x is easting or longitude whatever the
CRS's axis order". `mesh` adds "in another CRS than the mesh's, the box
holding the image of its edges there is meshed, snapped outward to the
DEM's nodes, or with --out-crs to its grid's". The phrase "snapped outward
to the DEM's nodes" is kept because `test_mcp_schema.py:433` reads it. On
the `--out-crs` path the box in the target CRS becomes the domain of
`target_grid_for` (`dem_input.py:344-352`), whose bounds, grown by a cell's
diagonal, snap outward to the target grid. The `--out-crs` help loses
"--bbox is then in this CRS." (`cli.py:725`). `fetch --out-crs` becomes
"The mesh's CRS, which --domain or --bbox is moved into, as mesh moves it."

**The CRS each command works in (the frame).** The box is moved into it,
then used as today:
- `mesh`: `DemRequest` keeps `bounds: Bounds | None`, in the DEM's CRS (in
  `target_crs` on the reprojected path), since it is plugin API
  (`test_plugin_api.py:47-50`; a plugin calls `DemRequest(bounds=grown)`). It gains `box: CrsBox | None`, last, so
  existing calls are unchanged. At most one of `bounds`, `box` and `domain`
  is given. `plan_dem`, once it has the footprints, turns `box` into
  `bounds` by `to_crs` into the frame: `target_crs` on the reprojected path,
  else `single_crs(footprints)`. With no footprints the box is not moved
  and the plan refuses as today. The `SIGNATURES` entry for `DemRequest`
  (`test_plugin_api.py:47-50`) gains `box=None` in its own commit, named in
  the PR body, as that file asks.
- `fetch`: `FetchRequest.box: CrsBox | None`. `_frame_box` moves it into the
  frame (`out_crs`, else the source) by `to_crs`. `region_sha256` hashes the
  box's CRS as well as the frame. `box_in_crs`'s docstring says what its
  move gives: "moved with `transform_bounds` (21 points an edge, so the
  image may fall inside the true one by metres on a box of a few hundred
  kilometres)", dropping "outward".
- `sources`: `listing(..., box: CrsBox | None)`. `box_covered` is measured
  on the area: the share of the box's lon/lat outline,
  `shapely.Polygon(box.outline("EPSG:4326"))`, that the coverage holds. A
  point or line box (zero area) keeps today's covered-or-not test, on
  `box.to_crs("EPSG:4326")`. `to_crs` runs first in either case, so an
  image across ±180° is refused with its words. `_lonlat`
  (`source_report.py:218-223`) stays for the tiles' coverage, where its 8
  segments a side are enough (`SEGMENTS`, `source_report.py:47`).

**Printed boxes name their CRS.** Each refusal or hint that prints a box
ends it with "in <CRS>" (`crs_label`): `fetch/plan.py:130` and `:192`,
`io/repository.py:365`, `mosaic.py:179` (the tiles' CRS), and `NotCached`.

**Hints.** `_fetch_run` writes `--bbox ... --bbox-crs <the box's CRS>`.
`NotCached` gains `crs: str | None`, set at `io/repository.py:381` from
`plan.meta.crs`; its words say "in EPSG:25833" (the CRS) instead of "in the
source's CRS" (`io/cog.py:75`). `_bbox_fetch` writes `--bbox-crs` from it.
Each printed command must run as printed (test 8).

**MCP.** No change in `schema.py`. Each of the three tools gains a
`bbox_crs` property, `anyOf [string, null]`, whose description is the help
above. `bbox`'s description says "required with it". JSON Schema has
`dependentRequired`, but click has no notion of it, and deriving it from
names would wrongly bind `domain_crs`. So an agent reading only the schema
learns the pairing from the two descriptions. A call that omits `bbox_crs`
gets the refusal's words above, naming `--bbox-crs`, as every refusal names
its flag today. `INSTRUCTIONS` (`server.py:93`) becomes "Give it bbox (west
south east north) and bbox_crs (EPSG:4326 for degrees) ...".

**Docs.** `README.md:24,26` and `:213` ("with `bbox` and `bbox_crs`"),
`INSTALL.md:239`, `project_structure.md:143` (plus the `domain.py` entry
naming `CrsBox`) and the script `tests/fixtures/velhas/extract.py`.
Benchmark records under `docs/benchmarks/` are history and stay as they are.

**Exclusions.** pyproj's area of use is not checked: a box far outside a
projection's valid area moves to finite nonsense, as a domain does today. A
box crossing ±180° in its own CRS stays inexpressible (`x_min <= x_max`),
and one whose image crosses it, or holds a pole, is refused, not split.

**Constants.** `EDGE_STEPS = 2000` (`domain.py`): steps along each side of
the box. It assumes a box up to São Francisco's size (12° x 14°, about
1300 x 1550 km), and was checked at that box, where the edge falls short by
at most 1.6 mm. The error grows with the square of the box's side, so a box
four times larger falls short by about 26 mm (an estimate, not probed).
**Timing.** One `to_crs` per run: 1.1 ms for the outline at São Francisco
scale, 8.9 ms with the transformer's setup on the quad
(`114-probes/box_probe.txt`). No timing gate. A box given in the frame's
CRS takes today's numbers bit for bit, so meshes are unchanged. `@perf`'s
acceptance run is not needed; that is the main session's call.

**Ideas file.** This overlaps the entry "`--bbox` help: 'its nodes' is
ambiguous with `--out-crs`". The help keeps "the DEM's nodes", already
reworded in `2f297c6f`, and adds "or with --out-crs to its grid's", so that
entry can close.

## 4. Tests (`@tester`)

Fixtures are measured by `docs/increments/114-probes/box_probe.py`, which
applies the design's move (`EDGE_STEPS = 2000`); `edge_probe.py` measures
the São Francisco error.
- Quad: `quadrants(whole(9, 13))`, 13 x 9 nodes, x 500000..500120 (dx 10),
  y 6599960..6600000 (dy 5), EPSG:25833.
- Sources: `dataset_fixtures.py` `COVERAGE`, `INSIDE`, `OUTSIDE`.
- Fetch: `test_cli_fetch.py` `BBOX`, EPSG:25833.

1. **Unit, `CrsBox`.** The same CRS returns the tuple bit for bit
   (`EPSG:4326` to `OGC:CRS84` too: `same_crs` is True). EPSG:25833
   `(500022, 6599969, 500087, 6599991)` to EPSG:4326 gives
   `(15.000389013, 59.538070744, 15.001538378, 59.53826831)` to 1e-9. The
   point box `(500000, 6600000)` in EPSG:25833 gives the point
   `(15.0, 59.53834913)` in EPSG:4326, to 1e-8. Refused, each naming both
   CRSs where two are involved: a non-finite value; `x_min > x_max`; an
   unreadable CRS (`ValueError`); `y_max = 90.5` in EPSG:4326; the box
   `(700000, 5000000, 900000, 5100000)` in EPSG:32660 to EPSG:4326 (across
   180°; pyproj gives 179.54 to -177.84); the box `(-500000, -500000,
   500000, 500000)` in EPSG:3413 to EPSG:4326 (holds the pole).
2. **Mesh, same CRS.** Quad, `--bbox 500022 6599969 500087 6599991
   --bbox-crs EPSG:25833`, stride 1. Expect 56 points, x 500020..500090,
   y 6599965..6599995 (cols 2-9, rows 1-7), equal to master's run without the
   option.
3. **Mesh, box in degrees.** The same run with the lon/lat box of test 1 and
   `--bbox-crs EPSG:4326`. Expect the same 56 points, bit for bit. The
   margin: the round trip moves the box's edges by at most 0.9 mm, and its
   edges lie 2-9 m from a node line.
4. **Refusals.** For each of `mesh`, `fetch`, `sources`, parametrized:
   `--bbox` alone exits 2 with the words of section 3, and nothing is
   written or requested. `--bbox-crs` alone exits 2, "applies only with
   --bbox". `--bbox-crs EPSG:0` exits 2, param `--bbox-crs`. The EPSG:32660
   box of test 1 to `sources` exits 2 with "crosses ±180°".
5. **Sources.** `INSIDE` in EPSG:25833, `(508522.7, 6576808.3, 511380.3,
   6582383.5)`: `box_covered` 1.0. `OUTSIDE` in EPSG:25833 `(555607.8,
   6651832.7, 561353.6, 6663057.1)`: 0.0. A box over `COVERAGE`'s east
   edge, `(517060.4, 6573496.4, 522780.6, 6579094.1)` in EPSG:25833: 0.4964
   (to 1e-4; the share by bounds would be 0.4961, so the test tells area
   from bounds). Both are measured against the union of the four tiles'
   `_lonlat` outlines, which is what the listing holds; see "Rulings on the
   red suite", 1. The UTM point `(508522.7, 6576808.3)` gives 1.0.
6. **Fetch.** A dry run with `BBOX` and `--bbox-crs EPSG:25833` prints the
   same plan lines as master. The same box in EPSG:4326, `(15.001791152,
   59.536663469, 15.005904106, 59.537800446)`, prints identical plan lines.
   The round trip moves the box's edges by at most 13.5 mm; the grown fetch
   box's edges lie at least 1.3 m from a node line.
7. **MCP schema.** `mesh`, `fetch` and `sources` each have `bbox_crs`:
   `anyOf [string, null]`, default null, and a description containing
   "Never assumed" and "Required with --bbox.". `bbox`'s description
   contains "--bbox-crs". `INSTRUCTIONS` contains "bbox_crs".
8. **Hints run.** A mesh refused for a missing cache prints a `rasputin
   fetch ... --bbox ... --bbox-crs EPSG:25833` command, and the `NotCached`
   words name "in EPSG:25833". Run as printed, the command exits 0, and the
   mesh then runs (extend `test_cli_sources.py:136,166`).
9. **Migration.** Every existing `--bbox` call in `tests/` (64 lines in 18
   files, `git grep -c -- '--bbox\|"bbox"\|bbox=' -- tests`) gains its
   fixture's CRS. `FetchRequest(box=Bounds(...))` (14 to 17 calls,
   `git grep -n -A4 'FetchRequest(' -- tests src_python | grep box=`)
   becomes `box=CrsBox(...)`. `DemRequest(bounds=...)` is unchanged.

## 5. Size

These are estimates of counted lines.

| Module | Net |
|---|---|
| `domain.py` (`CrsBox`, `outline`, refusals, `EDGE_STEPS`) | 32 |
| `cli.py` (3 options, helper, 3 call sites, hints; less two `Bounds` and a check) | 28 |
| `dem_input.py` (`box` field, one-of check, move in `plan_dem`) | 10 |
| `fetch/plan.py` | 6 |
| `source_report.py` (`listing`, `_cover` on the outline) | 5 |
| `io/cog.py`, `io/repository.py`, `mosaic.py` (CRS in words and `NotCached`) | 6 |
| `mcp/server.py` | 1 |
| **Total** | **about 90; worst case 120** |

## Rulings on the red suite

`@architect`, 2026-10-30, on `@tester`'s four pins in `b89da23d`,
`abccbc7f` and `986e7bab`.

1. **Test 5's figure: 0.4964, confirmed.** The design's 0.5183 was measured
   against `COVERAGE` as a plain degree rectangle. The listing measures
   against the union of the tiles' own `_lonlat` outlines
   (`source_report.py:249`, `:260`), so the test must too. `114-probes/share_probe.py`
   prints all four figures for the edge box: tiles by area 0.496418, tiles by
   bounds 0.496067, rectangle by area 0.518333, rectangle by bounds 0.517871
   (`share_probe.txt`). `INSIDE` and `OUTSIDE` give 1.0 and 0.0 every way.
   0.4964 and 0.4961 are 3.5e-4 apart, so `abs=1e-4` still tells area from
   bounds. Section 4, test 5 now says 0.4964.
2. **`docs/benchmarks/quick/cases.toml`: `@developer` edits it**, in the
   commit that adds `--bbox-crs` to `mesh`, so that no commit has a case the
   CLI refuses. It is a live input of the quick check, not a benchmark record,
   so section 3's "records stay as they are" does not cover it. The edit is
   the one `test_bench_quick.py:437-441` pins: `--bbox-crs EPSG:31983` after
   the velhas box. Precedent: `7974b15f` (`@developer` added the case) and
   `62e321f9` (`@architect` renamed options in it). `@perf` need not act. The
   quick check compares input fingerprints, not arguments
   (`perf-quick-check.md`, section 3), and a box in the frame's CRS meshes bit
   for bit, so the saved baselines stay valid. A quick check run with
   `--tree` at a commit before this PR now fails on velhas, since that tree
   has no `--bbox-crs`.
3. **Wordings: all three confirmed.** The non-finite refusal contains
   "finite" (section 3 lists the invariant but gave no words). The latitude
   check is a `CrsBox` invariant, so it names `--bbox` (section 3's table,
   "a broken invariant"). The ±180° and no-image refusals come from `to_crs`,
   so they name `--bbox-crs` (the table's header). Both hints print `--bbox
   X0 Y0 X1 Y1 --bbox-crs C` as one group, then `--out-crs` where
   `_fetch_run` has one, then `--cache`. `fetch` moves that box into its
   `--out-crs` as `mesh` does, so the printed command fetches the same window.
4. **Import edges: all three allowed.** `test_layering.py`'s rule is that a
   module imports from its own layer or below. `mosaic` to `crs` is within L1
   (values and pure algorithms). `io.repository` (L2) to `crs` (L1) points
   down, and L2 already imports `crs` (`io.geotiff`, `io.geojson`).
   `source_report` (L4) to `domain` (L1) points down. No cycle: `crs` imports
   only `io.models` and `run_record`, and `domain` only `crs` and `io.models`.
   Nothing new reaches `_core` or C++.

## Review

- Design round 1: CHANGES REQUESTED, f0dfe61b..0d64a829, 0 (sound at its core; fixes go to code review, no second design round; estimate low: about 90 lines, migration 64 CLI lines in 18 files plus about 17 `FetchRequest(box=)` and about 10 `DemRequest` sites; holds: same-CRS bit for bit, pole boxes, a projected point box, the MCP schema; blocking: `DemRequest` is plugin API (`test_plugin_api.py:47-50`, the lsm radiation search calls `DemRequest(bounds=grown)`), so keep `bounds` and add `box: CrsBox | None`; a box whose geographic image crosses ±180° comes back west > east (EPSG:32660 to 4326: 179.54, -177.75), `to_crs` must refuse naming both CRSs; edges are not held, 21 points per edge falls 11.2 m short on a São Francisco-sized box into EPSG:32723, so "outward" is false; with `--out-crs` the box snaps to `target_grid_for`, not the DEM's nodes; `box_covered` must use the area, reusing `_lonlat` (`source_report.py:218`); refuse |lat| > 90 in a geographic CRS; map gap `README.md:213`; cli.py lines after about 900 moved in fbeaf649, +21 to +33).
- Code round 1: CHANGES REQUESTED, 1242528e..252ad6ac, 107 (blocking: a `--bbox` with no area in the mesh's CRS is refused naming `--dem`, not `--bbox`, on mesh, mesh --out-crs and MCP mesh; refuse it under `--bbox`, with a test. Wording: project_structure.md:142 `#   imports io.models, crs and resources first-party`. Follow-ups: line box off the central meridian measured by its bounding box in `sources`; `NotCached` prints "in None"; `_plan_reprojected`'s "--out-crs needs a domain or a box in the target CRS" no longer fits). `Review: CHANGES REQUESTED 252ad6accabed50d018e01f6a67c730c47085ea4`
- Code round 2: APPROVED, 1242528e..600e0c4e, 120 (round 1's blocker fixed on mesh, mesh --out-crs and MCP mesh by the reviewer's own probe; @perf not needed, a box in the mesh's CRS meshes bit for bit). `Review: APPROVED 600e0c4e3c798ae959eda4291e418bd61fc3d4f3`

## As built

- Green at `19e41026`, master merged at `c3309f5e`, code round 1 fixes red `eee26bb1`, green `600e0c4e`; 120 counted lines against master `1242528e` (estimate 90, worst case 120).
- Pinned beyond the design by `@tester`, judged sound in code round 2: a box with no area is refused under `--bbox` naming the CRS it ended up in (`NoAreaError`); `--out-crs` with neither domain nor box is refused under `--out-crs` (`TargetRegionError`); `NotCached` refuses a box without a CRS (`TypeError`).
- Open follow-up: in `sources`, a line box off the central meridian is measured by its bounding box.
