# Increment 114: every --bbox names its coordinate system

**Status: designed** by `@architect`, 2026-10-29, on master `f0dfe61b`,
worktree `bbox-crs`. One PR, full loop, Python only. Number 114: 113 is
the highest under `docs/increments/` over every worktree (`lsm-thin-soil`),
112 the highest on the private repo's refs. Every `file:line` below is at
`f0dfe61b`. PR #404 (`fbeaf649`, merged since) changes `cli.py`, but none of
its hunks touch a line cited here. The branch `worktree-prefer-holes` may
shift these line numbers.

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

## 2. Every box today (master `f0dfe61b`)

| Where | Assumed CRS today | Stated at |
|---|---|---|
| `mesh --bbox` (MCP `mesh.bbox`) | `--out-crs` if given, else the DEM's | `cli.py:596-603`, `cli.py:725` (`--out-crs` help); applied in `dem_input.py:214` (projected path) and `dem_input.py:344-347` (target CRS) |
| `fetch --bbox` (MCP `fetch.bbox`) | `--out-crs` if given, else the source's | `cli.py:2599-2607`; `fetch/plan.py:38-40`, `_frame_box` `fetch/plan.py:114-122`, `region_sha256` `fetch/plan.py:163-172` |
| `sources --bbox` (MCP `sources.bbox`) | WGS 84 degrees | `cli.py:2748-2756`; `source_report.py:373-377`; MCP `INSTRUCTIONS`, `mcp/server.py:93` |
| hint `_fetch_run` (mesh refusal: "run: rasputin fetch ...") | mesh's own rule, passed on | `cli.py:1602-1616`, used at `cli.py:977-1017`, `1212` |
| hint `_bbox_fetch` (window box of a `NotCached`) | the source's, unstated | `cli.py:1619-1626`; box set at `io/repository.py:381`; `NotCached` `io/cog.py:56-69`; used at `cli.py:1211`, `3705` (catchment) |

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
        def to_crs(self, dst: str | CRS) -> tuple[float, float, float, float]

- Invariants: every value is finite; `x_min <= x_max` and `y_min <= y_max`.
  Equal values are allowed because `sources` takes a point. `mesh` and
  `fetch` still refuse an empty box through `Bounds`, after the move.
- `to_crs` returns the numbers untouched (bit for bit) when
  `same_crs(crs, dst)`. Otherwise it returns `transform_bounds(crs, dst, box)`
  (densified with 21 points, outward, `crs.py:97-105`).
- A non-finite result is refused with a `ValueError` naming both CRSs: "the
  box (...) in EPSG:x has no image in EPSG:y".

**CLI.** `mesh`, `fetch` and `sources` each gain the option
`--bbox-crs TEXT`, help: "The --bbox's CRS, anything pyproj reads; EPSG:4326
is longitude and latitude in degrees. Never assumed." One helper in `cli.py`
turns `(bbox, bbox_crs)` into a `CrsBox | None` for all three. It refuses:

| Case | Exit | Words (param `--bbox-crs`) |
|---|---|---|
| `--bbox` without `--bbox-crs` | 2 | "--bbox needs --bbox-crs, the CRS its four numbers are in (EPSG:4326 for longitude and latitude in degrees); a box's CRS is never assumed" |
| `--bbox-crs` without `--bbox` | 2 | "applies only with --bbox" (as `--domain-crs`, `cli.py:1047-1048`) |
| CRS pyproj cannot read | 2 | `parse_crs`'s own words |
| no image in the frame | 2 | `to_crs`'s words |

The help of each `--bbox` drops its assumption and says "XMIN YMIN XMAX YMAX
in --bbox-crs (required with it); x is easting or longitude whatever the
CRS's axis order". `mesh` adds "in another CRS than the mesh's, the box
holding its image there is meshed, snapped outward to the DEM's nodes". The
phrase "snapped outward to the DEM's nodes" is kept because
`test_mcp_schema.py:433` reads it. The `--out-crs` help loses "--bbox is
then in this CRS." (`cli.py:725`). `fetch --out-crs` becomes "The mesh's CRS,
which --domain or --bbox is moved into, as mesh moves it."

**The CRS each command works in (the frame).** The box is moved into it,
then used as today:
- `mesh`: `DemRequest.bounds: Bounds` becomes `box: CrsBox | None`.
  `plan_dem` moves the box into `single_crs(footprints)` on the projected path
  and into `target_crs` on the reprojected path. With no footprints the box
  is not moved and the plan refuses as today.
- `fetch`: `FetchRequest.box: CrsBox | None`. `_frame_box` moves it into the
  frame (`out_crs`, else the source). `region_sha256` hashes the box's CRS as
  well as the frame.
- `sources`: `listing(..., box.to_crs("EPSG:4326"))`. Its west-east check
  (`cli.py:2769`) runs after the move, so an image that crosses ±180° is
  refused, with "in EPSG:4326" added to the words.

**Hints.** `_fetch_run` writes `--bbox ... --bbox-crs <the box's CRS>`.
`NotCached` gains `crs: str | None`, set at `io/repository.py:381` from
`plan.meta.crs`. `_bbox_fetch` writes `--bbox-crs` from it. Each printed
command must run as printed (test 8).

**MCP.** No change in `schema.py`. Each of the three tools gains a
`bbox_crs` property, `anyOf [string, null]`, whose description is the help
above. `bbox`'s description says "required with it". JSON Schema has
`dependentRequired`, but click has no notion of it, and deriving it from
names would wrongly bind `domain_crs`. So an agent reading only the schema
learns the pairing from the two descriptions. A call that omits `bbox_crs`
gets the refusal's words above, naming `--bbox-crs`, as every refusal names
its flag today. `INSTRUCTIONS` (`server.py:93`) becomes "Give it bbox (west
south east north) and bbox_crs (EPSG:4326 for degrees) ...".

**Docs.** `README.md:24,26`, `INSTALL.md:239`, `project_structure.md:143`
(plus the `domain.py` entry naming `CrsBox`) and the script
`tests/fixtures/velhas/extract.py`. Benchmark records under `docs/benchmarks/`
are history and stay as they are.

**Exclusions.** pyproj's area of use is not checked: a box far outside a
projection's valid area moves to finite nonsense, as a domain does today. A
box crossing ±180° stays inexpressible (`x_min <= x_max`).

**Constants.** None new. `transform_bounds`'s 21 points are existing.
**Timing.** One `transform_bounds` per run, 1.9 ms
(`114-probes/box_probe.txt`). No timing gate. A box given in the frame's CRS
takes today's numbers bit for bit, so meshes are unchanged. `@perf`'s
acceptance run is not needed; that is the main session's call.

**Ideas file.** This overlaps the entry "`--bbox` help: 'its nodes' is
ambiguous with `--out-crs`". The help keeps "the DEM's nodes", already
reworded in `2f297c6f`, so that entry can close.

## 4. Tests (`@tester`)

Fixtures are measured by `docs/increments/114-probes/box_probe.py`.
- Quad: `quadrants(whole(9, 13))`, 13 x 9 nodes, x 500000..500120 (dx 10),
  y 6599960..6600000 (dy 5), EPSG:25833.
- Sources: `dataset_fixtures.py` `COVERAGE`, `INSIDE`, `OUTSIDE`.
- Fetch: `test_cli_fetch.py` `BBOX`, EPSG:25833.

1. **Unit, `CrsBox.to_crs`.** The same CRS returns the tuple bit for bit
   (`EPSG:4326` to `OGC:CRS84` too: `same_crs` is True). EPSG:25833
   `(500022, 6599969, 500087, 6599991)` to EPSG:4326 gives
   `(15.000389013, 59.538070744, 15.001538378, 59.53826831)` to 1e-9. A point
   box stays a point. A non-finite value and `x_min > x_max` are refused. An
   unreadable CRS is a `ValueError`.
2. **Mesh, same CRS.** Quad, `--bbox 500022 6599969 500087 6599991
   --bbox-crs EPSG:25833`, stride 1. Expect 56 points, x 500020..500090,
   y 6599965..6599995 (cols 2-9, rows 1-7), equal to master's run without the
   option.
3. **Mesh, box in degrees.** The same run with the lon/lat box of test 1 and
   `--bbox-crs EPSG:4326`. Expect the same 56 points, bit for bit. The
   margin: the round trip grows the box by at most 0.9 mm, and its edges lie
   2-9 m from a node line.
4. **Refusals.** For each of `mesh`, `fetch`, `sources`, parametrized:
   `--bbox` alone exits 2 with the words of section 3, and nothing is
   written or requested. `--bbox-crs` alone exits 2, "applies only with
   --bbox". `--bbox-crs EPSG:0` exits 2, param `--bbox-crs`.
5. **Sources.** `INSIDE` in EPSG:25833, `(508522.7, 6576808.3, 511380.3,
   6582383.5)`, moves back to `(15.14978, 59.329933, 15.200294, 59.380067)`,
   inside `COVERAGE`. Expect `box_covered` 1.0. `OUTSIDE` in EPSG:25833
   `(555607.8, 6651832.7, 561353.6, 6663057.1)` gives 0.0. A UTM point gives
   1.0 or 0.0 the same way.
6. **Fetch.** A dry run with `BBOX` and `--bbox-crs EPSG:25833` prints the
   same plan lines as master. The same box in EPSG:4326, `(15.001791152,
   59.536663469, 15.005904106, 59.537800446)`, prints identical plan lines.
   The round trip grows the box by at most 13.5 mm; the grown fetch box's
   edges lie at least 1.3 m from a node line.
7. **MCP schema.** `mesh`, `fetch` and `sources` each have `bbox_crs`:
   `anyOf [string, null]`, default null, and a description containing
   "Never assumed". `bbox`'s description contains "--bbox-crs".
   `INSTRUCTIONS` contains "bbox_crs".
8. **Hints run.** A mesh refused for a missing cache prints a `rasputin
   fetch ... --bbox ... --bbox-crs EPSG:25833` command. Run as printed, it
   exits 0, and the mesh then runs (extend `test_cli_sources.py:136,166`).
9. **Migration.** Every existing `--bbox` call in `tests/` (about 60 in 17
   files, `git grep -c -- '--bbox\|"bbox"\|bbox=' -- tests`) gains its
   fixture's CRS. `DemRequest(bounds=...)` (4 calls) becomes `box=CrsBox(...)`.

## 5. Size

These are estimates of counted lines.

| Module | Net |
|---|---|
| `domain.py` (`CrsBox`) | 22 |
| `cli.py` (3 options, helper, 3 call sites, hints) | 30 |
| `dem_input.py` | 8 |
| `fetch/plan.py` | 6 |
| `io/cog.py`, `io/repository.py` | 4 |
| `mcp/server.py` | 1 |
| **Total** | **about 70; worst case 110** |

## Review

(none yet)

- Design round 1: CHANGES REQUESTED, f0dfe61b..0d64a829, 0 (sound at its core; fixes go to code review, no second design round; estimate low: about 90 lines, migration 64 CLI lines in 18 files plus about 17 `FetchRequest(box=)` and about 10 `DemRequest` sites; holds: same-CRS bit for bit, pole boxes, a projected point box, the MCP schema; blocking: `DemRequest` is plugin API (`test_plugin_api.py:47-50`, the lsm radiation search calls `DemRequest(bounds=grown)`), so keep `bounds` and add `box: CrsBox | None`; a box whose geographic image crosses ±180° comes back west > east (EPSG:32660 to 4326: 179.54, -177.75), `to_crs` must refuse naming both CRSs; edges are not held, 21 points per edge falls 11.2 m short on a São Francisco-sized box into EPSG:32723, so "outward" is false; with `--out-crs` the box snaps to `target_grid_for`, not the DEM's nodes; `box_covered` must use the area, reusing `_lonlat` (`source_report.py:218`); refuse |lat| > 90 in a geographic CRS; map gap `README.md:213`; cli.py lines after about 900 moved in fbeaf649, +21 to +33).
