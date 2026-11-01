# Increment 122: NVE's lakes overrule the land cover, and no triangle is left without a class

Status: **design** (`@architect`, 2026-11-01, on public master `ae33748e`;
branch `worktree-lakes-overrule`; amended the same day for Ola's ruling
on coverage gaps: a part of the domain no input file covers stops the run,
and only `--land-cover-partial` fills it, sections 3.6 and 7). Full lane (it changes the class each
triangle carries: a wrong class gives wrong numbers without a crash),
Python only: no C++, nothing under `include/terrain/refinement/` or
`include/terrain/mesh/`, so no `@perf` acceptance run (section 8). Pipe
item `mesh/lakes-overrule`. Number 122: 121 is the highest under
`docs/increments/` over every local and remote-tracking ref and every
worktree (`git ls-tree` over each; 120 and 121 are private increments on
`worktree-lsm-*` branches). Built **after** increment 119 merges (public
PR 419; section 5).

## 0. The plain picture

- **Today**, a run with a land-cover map and NVE's lakes,
  `--land-cover clc2018-clean-25833 --prefer lake=<NVE lakes>`, already
  lets each NVE lake decide inside it (increment 107). What the land cover
  calls water outside every NVE lake takes the class of the neighbour it
  shares the longest border with, which may itself be water (a river);
  an island inside an NVE lake where the land cover had water takes the
  class of the nearest polygon, water or not; and a piece with neither is
  left with no class, code 0. Any part of the domain no land-cover
  polygon covers is code 0 too (`land_cover_codes`: "0 = in no polygon").
- **After**, those pieces never become water, and **no triangle of a
  finished mesh carries code 0**: each connected group of triangles that
  no polygon covers takes the class of the nearest triangle whose class is
  not water (Ola, 2026-11-01). Each such triangle is marked
  (`land_cover_filled` = 1 in the file), and the run record counts them and
  their area. Code 0 is left only on the constraint lines a `.vtk` writes
  as cells.
- **Missing data stops the run** (Ola, 2026-11-01: "If you have missing
  data, that should be an error. If you by some strange reason want to
  continue making geometries with only partial land cover, you need a
  manual override."). A part of the domain that no `--land-cover` or
  `--prefer` polygon covers is a coverage gap: the run exits 2 before the
  mesh is built, naming each gap's area and place in the DEM's CRS and the
  catalogue file that covers it (Bygdin with AR50 Innlandet alone: "add
  ar50-vestland-clean-25833"). `--land-cover-partial` lets the run go on;
  the gap is then filled and marked like an island, and recorded. An
  island (land cover had water, NVE has none) is not a gap: some file
  covers it, so it fills without the flag.
- **Ola's Bygdin run without `--prefer`** (CORINE alone): 394 of the 1,180
  triangles inside NVE's lakes carry a CORINE land class, 5.61 km². With
  `--prefer lake=` they are 143 triangles, 0.41 km², every one within
  28.5 m of NVE's shore, inside the border clean-up's 50 m shift
  (section 2). The overrule needs NVE's lakes named in the run; this
  increment does not add them by itself.

## 1. Prior art: legacy and literature

*Literature.* Filling cells with no class from the nearest classed cell is
the categorical nearest-neighbour fill of raster GIS: Esri's Nibble
"replaces cells of a raster corresponding to a mask with the values of the
nearest neighbors" ([ArcGIS Pro, Nibble](https://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/nibble.htm)),
and its `zones` option restricts the donors as we restrict them to
non-water classes. The tie-break by the longest shared border is the rule
increment 107 took from map generalisation (AR50: "Areas under 15 decares
are merged into neighbours"; 107, section 1). Departure: the unit filled
is a connected group of triangles, not a cell, so one island takes one
class; and the distance is between triangles as closed sets, not cell
centres. No novelty is claimed; none searched.

*Legacy.* `git grep -n -i -E "nearest|unclassified|no_class|fill" legacy-archive -- legacy`
finds only `ExposedHeightFill` in `legacy/rasputin/web_visualize.py`, a
display parameter. Nothing to carry.

*Ideas file.* Overlaps "20c-3: close gaps between land-cover polygons"
(the slits it closes would be filled here if any were left, and one wider
than 2 cm is now refused, section 3.6); not taken up. Takes up no entry.

*The gap check.* GEOS's coverage validation (shapely 2.2's
`coverage_invalid_edges`, "the `gap_width` parameter allows to specify
the maximum width of gaps to detect") judges gaps by width, as ours does,
but the other way round: it flags gaps narrower than the width as slivers
and allows wider holes; ours measures against the domain, refuses every
part wider than `GAP_NOISE` and drops only the narrower noise. It also
needs exactly matching edges, which two files need not have. Tiling an
overlay into threads is ordinary; no novelty claimed, none searched.

## 2. Bygdin, measured

Scratch probes `lake_probe.py`, `map_probe.py`, `island_count.py` and
`fill_proto3.py` (scratchpad of this session), on master `ae33748e`, M1
Max: `rasputin mesh --dem DTM10_UTM33_20260925 --domain
docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson --tolerance
10 --edges all --threads 8`, NVE's lakes the 201 polygons of
`rasputin_scratch/temp-check-20261021/nve-lakes-25833.geojson`, all in
EPSG:25833. "In an NVE lake" is a triangle's centroid inside their union.

| `--land-cover`, `--prefer` | Triangles | In an NVE lake | of them land | Code 0 |
|---|---|---|---|---|
| `clc2018-clean-25833`, none | 54,429 | 1,180 | 394 (5.61 km²) | 0 |
| `clc2018-clean-25833`, `lake=` NVE | 59,623 | 3,684 | 143 (0.41 km², all within 28.5 m of the shore) | 0 |
| `ar50-innlandet-clean-25833`, `lake=` NVE | 125,811 | 11,046 | 63 (0.003 km²) | 2,245 (7.86 km², 5 groups) |
| the same plus `ar50-vestland-clean-25833` | 125,973 | 11,046 | 63 | 0 |

- **The map is whole.** With CORINE and NVE, `land_cover_map.build`'s
  areas cover the domain (304.91 km²) with no gap over 1 m², no land
  inside NVE's lakes and no water outside them. The 143 and 227 border
  triangles (land in a lake, water outside one) straddle the shore: half
  of each lies on either side (median share inside NVE 0.54 and 0.24); the
  border clean-up (`edges_max_shift_m` 50) placed the constraint line, not
  the class rule.
- **Islands.** `preferred._islands` fills 61 pieces today (0.108 km², 143
  triangles), all with land classes; on the 29 October mesh, before it
  existed, the same islands were code 0 (145 triangles, 0.108 km²). Left
  out of the map instead (prototype: `_islands` returning nothing), the
  mesh has the same vertices and the same set of triangles, and exactly
  145 triangles turn 0; the prototype fill of section 3.2 gives 142 of
  them the class `_islands` gives today. No refilled piece took a water
  class on Bygdin (0 of 1,011).
- **Coverage gaps.** AR50 Innlandet alone leaves 7.86 km² west of the
  county line uncovered, and two islands lie there: code 0 today. Adding
  Vestland's AR50 covers them.
- **The gap, as geometry** (scratch `gap_probe.py`: the CLI's reads and
  `build`, then the domain less the union of the joined base areas and
  NVE's lakes; EPSG:25833 throughout). CORINE: 0 parts, 0 m². AR50
  Innlandet: 5 parts, 7,863,437 m²: 5.820 km² (largest inscribed circle
  radius 590 m, centre 136216 6830714), 1.993 km² (453 m), 43,928 m²
  (73 m), 6,161 m² (25 m), 570 m² (10.6 m). Innlandet plus Vestland:
  **0 parts, 0 m² exactly**, so no sliver along the county line. The
  union and difference took 1.57 to 1.77 s with AR50, 0.22 s with
  CORINE (read and `build` 0.75 to 1.25 s). Of the 16 clean land-cover
  products in Ola's catalogue, two cover the largest gap's centre:
  `ar50-vestland-clean-25833` and `clc2018-clean-25833`; reading each at
  a 2 m box took 2.58 s for all 16.
- **Noise, and the cost at scale** (scratch `gap_numedal.py` and
  `gap_tiles.py`, no `--prefer`). Numedalslågen, 5,548 km², with the
  clean AR50 of Buskerud, Telemark, Vestfold, Vestland and Akershus
  (168,659 polygons after `build`; read and `build` 17.4 s): one union
  and difference took 111.5 s and left 6 parts, 0.00009 m² in all, every
  inscribed radius under 0.00005 m. By 5 km tiles in 8 threads (section
  3.6): 3.5 s, 48 pieces (the tile cuts split the noise), 0.000001 m² in
  all, every radius under 0.000001 m. On Bygdin the tiles give the one
  union's parts to the mm² (Innlandet alone, no lakes: 6,505,443.209,
  2,193,468.263 and 43,927.675 m², in 0.1 s against 1.57 s), and with
  Vestland one noise piece, dropped. The first tiled run shared one
  `shapely.prepare`d domain among the threads and died with a
  segmentation fault (exit 139); without the `prepare` it ran clean. One
  run each, so not proven the cause; the design prepares nothing shared.

## 3. Design

### 3.1 The refill never makes water; islands leave the map (`preferred.py`)

`replace(base, chosen, near, *, water)` gains `water: frozenset[int]`,
the class file's water codes (`LandCoverClasses.water_codes`, passed by
`land_cover_map.build`). A leftover piece of a replaced class takes the
class of the **non-water** kept polygons sharing the longest border with
it (107's rule, the same `BAND` and `DECIMALS`); a piece with no such
neighbour is left out of the map. `_islands` goes: an island is now a piece
left out, and the mesh fill gives it its class, marked. One rule for every
`--prefer` key: a refill never invents water, for a glacier key too.

### 3.2 The fill on the finished mesh (`landcover.py`)

```python
@dataclass(frozen=True)
class Filled:
    codes: npt.NDArray[np.int32]     # (T,), no 0
    filled: npt.NDArray[np.bool_]    # (T,), True where the code came from the fill
    groups: int
    largest: tuple[float, float, float]  # area, x, y of the largest group (0, nan, nan: none)

def fill_no_class(vertices: npt.ArrayLike, triangles: npt.ArrayLike,
                  codes: npt.ArrayLike, *, water: frozenset[int],
                  asking: npt.ArrayLike | None = None) -> Filled:
```

Pure (numpy and shapely, no file, no CRS object). `asking` (default
`codes == 0`) marks the triangles to fill; every other triangle keeps its
code. The rule, exactly:

1. **Groups.** The asking triangles split into groups connected through
   shared edges (constraint edges do not cut a group; `regions` on the
   asking triangles with no edges).
2. **Donors.** The triangles not asking whose code is neither 0 nor in
   `water`. If there are none, every triangle not asking whose code is not
   0 (a mesh of water and 0 only). If there are none either, `ValueError`:
   "the land cover gives no triangle of the mesh a class, so code 0
   cannot be filled; it does not cover the domain".
3. **Nearest.** A group's distance to a donor triangle is the smallest
   2-D Euclidean distance between their closed triangles (`shapely.distance`
   semantics; touching gives 0), in the mesh's own x and y: the CRS the
   file's `crs` string names, always projected (a DEM in degrees needs
   `--out-crs`, `OutCrsRequiredError`). The group takes the code of a
   donor at the smallest distance, exact float equality as GEOS returns it.
4. **Ties** among donors at that distance: the code whose donors share the
   most edge length with the group (edges common to a group triangle and a
   donor, summed per code, compared rounded to `preferred.DECIMALS`), then
   the lower code. So a group the shore touches takes 107's
   longest-border class, and an island across water takes the nearest
   shore's.

The result depends on no triangle or vertex order. `largest` is the group
of largest total area, ties to the group holding the lowest triangle
index, with the centroid of its largest triangle (same tie) as x, y.
Cost: only when some code is 0 (one `any` otherwise). The prototype took
0.17 s on 59,623 triangles (all 3,768 water triangles asking, 201 groups)
and 0.44 s on 125,973 (10,994 asking); querying the nearest per group's
rim triangle, not per group polygon, took it from 8.1 s to 0.44 s, so the
tree is queried with each group's rim triangles (those with an edge on
the group's border or on the mesh's) and the minimum taken per group.
The donors' tree may hold only donors with a vertex shared with a
non-donor or on the mesh border: the nearest point of the donor set lies
on its rim. An estimate, not measured: 10⁷ triangles cost tens of
seconds when some code is 0, in the edge table and the tree; nothing
when none is.

### 3.3 The shares

In the shares path (`--shares`, `_class_fractions`), a triangle's class-0
share (`polygon_fractions.assign`'s uncovered rest, or a raster's no-data
class) joins the share of one class, merged into its row
(`land_cover_fractions.merge_zero_share(table, into) -> ClassFractions`,
pure; `into[t]` the code triangle `t`'s 0 share joins): the triangle's
code after the fill when that is not water; otherwise the code
`fill_no_class(..., asking=<water-coded triangles with a 0 share>)` gives
it. The dominant code is filled first by section 3.2. A triangle whose 0
share moved is marked like one whose code was filled. So no
`land_cover_fractions` row and no cell carries class 0, and the Shyft
`cells` layer no longer counts an uncovered island as unspecified.

The gap rule reaches the shares path at two points. With class polygons,
the map is the one section 3.6 checked before the mesh, so a 0 share
left is an island's, a lone piece's or noise under `GAP_NOISE`, and
merges without the flag. With a raster, every 0 share inside the domain
is no data: without `--land-cover-partial` the run is refused (section
3.6) before `merge_zero_share`; with it, the shares merge as above.

### 3.4 The mark, the record, the words

- `.vtk`: an `int` cell array `land_cover_filled`, 1 where filled, 0
  elsewhere and on every line cell, written whenever `land_cover_code` is.
  `.ply`: a face property `uchar land_cover_filled` after
  `land_cover_code`. `.gpkg`: a `triangles` column `land_cover_filled`
  (`INTEGER` 0 or 1; NULL where `code` is NULL). `io/vtk_read.read_mesh`
  returns it (`MeshArrays.land_cover_filled`, None when absent).
- Record keys `land_cover_filled` (triangles marked) and
  `land_cover_filled_m2` (their area, plus the area of the 0 shares moved,
  in the mesh CRS's plane), both 0 on a run with nothing to fill.
- One stderr line when any is filled:
  `land cover: 145 triangles (0.1078 km2) in no polygon took the nearest non-water class: 61 groups, the largest 0.0162 km2 at 147268 6814797 in EPSG:25833`
  (Bygdin with CORINE and NVE, by the prototype; the words are free, the
  figures and the CRS are not). The largest group points a user at a forgotten
  file (Innlandet alone: 5.82 km² past the county line).
- `land_cover_codes` ends "0 = every constraint line; no triangle"
  instead of "0 = in no polygon, and every constraint line" (both paths,
  `src_python/tin_engine/cli.py@ae33748e:2148` and `:2565`).
- `accept_no_class` in the class files is unchanged: it still lets a legend
  colour code 0 for the lines.

### 3.5 Where it runs

`cli.py`: after `label_triangles` in `_land_cover` (whole mesh and each
`--pieces` piece through `_write_piece`), and after the codes of
`_class_fractions` (polygons and raster). A `ValueError` from the fill
becomes `typer.BadParameter(param_hint="--land-cover")`, exit 2, no file
written. A piece is filled on its own triangles: a group cut by a seam
may take a different class on each side (the pieces' promise covers the
vertices along a seam, not codes).

### 3.6 A coverage gap stops the run (`land_cover_map.py`, `cli.py`)

**What a gap is.** The part of the domain polygon (the `area` that
`_land_cover_map` reads around: the `--domain`, else the DEM window's
box, in the DEM's CRS) that lies in no polygon of any `--land-cover`
file and no polygon of any `--prefer` file, as read. It is computed from
the inputs, not from the map after `preferred.replace`: an island or a
lone water piece that `replace` leaves out lies inside a land-cover
polygon of a water class, so it is never part of a gap, whatever its
size. No area threshold decides; only geometry.

```python
@dataclass(frozen=True)
class Gap:
    parts: tuple[tuple[float, float, float], ...]  # (area m², x, y), area descending
    area: float                                    # sum of the parts' areas

def coverage_gap(domain: Polygon, covers: Sequence[BaseGeometry], *,
                 threads: int = 1) -> Gap:
```

Pure (shapely, no file, no CRS; x, y and areas in the domain's own
plane, which the caller names). `covers` is the joined base map's
polygons (`base.areas`, cut near the region by `build`) and the
`--prefer` polygons. The rule:

1. **Tiles.** The domain's box is cut into square tiles of `GAP_TILE`
   (5,000 m, chosen; section 2's probes); per tile, in `threads` threads,
   the tile's part of the domain less the union of the `covers` the tile
   meets (one `STRtree` query per tile), each clipped to the tile
   (`clip_by_rect`); no geometry shared among the threads is
   `shapely.prepare`d (section 2). Measured: 3.5 s on Numedalslågen
   against 111.5 s for one union; 0.1 s on Bygdin. `GAP_TILE` assumes
   projected metres; checked at 5,548 km² (1,152 tiles). An estimate,
   not measured: São Francisco's 637,000 km² with polygon land cover,
   about 25,500 tiles, some 80 s at Numedalslågen's polygon density (its
   tile count times 3.5 s / 1,152), against the hours its mesh takes.
2. **Parts.** The pieces of every tile, merged across tile borders
   (`union_all` of the pieces only), split into connected parts.
3. **Noise.** A part whose largest inscribed circle
   (`shapely.maximum_inscribed_circle(part, GAP_NOISE / 10)`: the
   default tolerance, a thousandth of the part's box, is coarser than the
   noise on a long strip) has radius at most `GAP_NOISE` (0.01 m,
   chosen) is dropped. Why a tolerance at all: Numedalslågen's
   five AR50 counties, which share borders cut from one national map,
   leave 6 parts of 0.00009 m² in all, the widest of radius under
   0.00005 m: overlay rounding, not missing data; refusing them would
   refuse a run with every file. Why this value: it assumes projected
   metres with coordinates up to 10⁷ (UTM northings); it is 200 times the
   widest noise measured (5,548 km², the largest input checked) and 1,000
   times below the narrowest real gap part measured (Bygdin, radius
   10.6 m). It is a width, not an area: a long thin strip with no data
   wider than 2 cm is still a gap.
4. **Place.** Each part's x, y is the centre of that circle: the point of
   the part farthest from data, where a user looks.

**Where and when.** In `_land_cover_map`, right after `built`, before the
lines are read, before the border clean-up and before the mesh, on the
whole domain once (so `--pieces` refuses before any piece is written).
This is the first moment the uncovered area is known, as 119 refuses
invalid polygons at the read. A raster land cover (`--land-cover` a class
GeoTIFF or a cached year) has no map before the mesh: there the gap is
the class-0 share (the raster's no-data) of the triangles inside the
domain, known only after `assign_fractions`; it is refused there, before
any file is written, with the same words and parts made of the connected
groups of triangles holding a 0 share (area: the 0 shares' area; place:
the centroid of the group's largest such triangle).

**The refusal.** `CoverageGap(FeatureError)` in `land_cover_map.py` (two
lines), turned by `cli.py` into `typer.BadParameter(param_hint=
"--land-cover")`: "Invalid value for --land-cover: …", exit 2, no file
written; the MCP tools get the words as their tool error. The words
(line breaks free; each part a line of its own; every part listed,
none held back):

    the land cover leaves 7.8634 km2 of the domain in no polygon of any
    --land-cover or --prefer file (5 parts); missing land cover is refused.
    Add the file that covers it, or pass --land-cover-partial to fill each
    gap from the nearest non-water class and mark it:
      5.8199 km2 around 136216 6830714 in EPSG:25833; covered by ar50-vestland-clean-25833
      1.9929 km2 around 141916 6830806 in EPSG:25833; covered by ar50-vestland-clean-25833
      ...

The CRS is `crs.crs_label(dem_crs)` on every line; areas to 4 decimals
of km², places to the metre. "covered by" lists the hint for that part,
and is absent when there is none.

**The hint.** It reads files, so it lives in `cli.py`, not in the
file-free `sources.py`: `_covering(gap, used, classes_name, known,
dem_crs) -> list[list[str]]` (per part, the ids covering its centre),
called only on refusal. Candidates: the catalogue's
local `features` entries whose `classes` is the run's class file by its
shipped name, not among the run's `--land-cover` ids; of a raw `X` and
its `clean_sibling` for the DEM CRS (119), only the sibling. Each
candidate is read with `read_class_polygons` at the box of all the
parts' centres widened by 1 m; a candidate is named for each centre one
of its polygons covers. A candidate that fails to read (missing file,
invalid polygons, `FeatureError`) is skipped: the hint is help, not a
check. No candidate is fetched (an `arcgis-features` entry is not
local). A class file given by path, or a raster land cover, gets no
hint; a `--land-cover` given by path still does (`used` is then empty). Bygdin: 15 candidates, about 2.5 s (section 2), all on the
refusal path.

**The override `--land-cover-partial`** (a boolean flag, default off;
help: "Mesh even where no land-cover file covers the domain: each gap
takes the nearest non-water class and is marked"). With it, the gap is
reported, not refused: one stderr line
`land cover: 7.8634 km2 of the domain in no file (5 parts, the largest 5.8199 km2 around 136216 6830714 in EPSG:25833), filled by --land-cover-partial`;
then the run goes on, and section 3.2's fill gives the gap's triangles
the nearest non-water class and marks them, exactly as an island's.
Nothing else differs: the mark already tells a reader which triangles
carry no data of their own. Record keys, written whenever a land cover is
given: `land_cover_partial` ("on" or "off"), `land_cover_gap_m2` (the
gap's area in the DEM CRS's plane, 0 when none) and `land_cover_gaps`
(its parts). The MCP mesh tool gains no field: a manual override is a
person's, and the agent gets the refusal naming the file to add
(question 1).

## 4. Blueprint

    read (119: valid polygons) -> land_cover_map.build -> preferred.replace(..., water=)
      -> coverage_gap(domain, base.areas + --prefer polygons, threads=) -> Gap
         (parts and no --land-cover-partial: _covering hint, CoverageGap, exit 2)
      -> map.areas (islands, lone water pieces and, with the flag, gaps uncovered)
      -> edges_of -> mesh -> label_triangles | assign/dominant (codes, 0 where uncovered)
      -> fill_no_class(vertices, triangles, codes, water=) -> Filled
      -> [shares: raster 0 shares and no flag: CoverageGap; merge_zero_share(table, into)]
      -> write_vtk/write_ply/write_triangles(..., filled=Filled.filled), record, stderr

Interfaces the tests call: `land_cover_map.coverage_gap`,
`land_cover_map.Gap`, `land_cover_map.CoverageGap`,
`landcover.fill_no_class`, `landcover.Filled`,
`land_cover_fractions.merge_zero_share`,
`preferred.replace(base, chosen, near, *, water)`,
`land_cover_map.build` (unchanged signature after 119), `read_mesh(...)
.land_cover_filled`, and `rasputin mesh` through `CliRunner`.

## 5. Depends on increment 119

119 (`core-assumes-clean`, PR 419, approved) makes every polygon reaching
`build` valid. This design relies on it: the refill and the islands left
out are overlays (`difference`, the shared-border bands) of valid input,
with no repair that can leave slivers to fill; and it is written against
119's code: `build(..., cut=...)` not `fix=`, `read_polygons` returning a
list, `_class_fractions` without `land_cover_repaired`. So the branch is
rebased on master after PR 419 merges, before `@tester` starts; begun on
`ae33748e`, the red suite would conflict in `cli.py`'s `_land_cover_map`
and `_class_fractions` and test 119's removed `repaired`.

## 6. Degeneracy policy

- No code 0 on the mesh: the fill does nothing and marks nothing.
- No donor at all: refused (section 3.2, step 2). Without
  `--land-cover-partial` such a mesh is a gap over the whole domain and is
  refused by section 3.6 first; step 2 is reached only with the flag.
- A domain no polygon reaches: one part, the whole domain; with no
  `--land-cover` file at all, no check (no land cover asked for).
- A gap part touching the domain border: still a gap (no data there).
- A gap with no candidate covering it: the words without "covered by".
- A mesh with no class file (a raster without one): `water` is empty, so
  every coded triangle is a donor.
- A group touching a donor only at a vertex: distance 0, shared length 0.
- Water-only donors: the group becomes water, marked like any other.

## 7. Tests `@tester` writes red

Unit grids: unit squares, `c` columns by 2 rows, each square two
triangles, codes per column (CORINE: 322 heath, 333 sparse, 511 river,
512 lake). Figures measured with the prototype (scratch `fixtures.py`).

1. `fill_no_class`, `[333, 512, 0, 512, 512, 322]` (24 triangles): one
   group of 4, all 333 (distance 1 against 2); `filled` true on those 4
   only; `groups` 1.
2. A tie at distance 1, `[333, 512, 0, 512, 322]` (20): 322, the lower.
3. The border decides, `[333, 0, (322 below, 512 above)]` (12): 333
   (shared 2 against 1), over the lower 322.
4. Water only, `[512, 0]` (8): 512, marked.
5. Nothing to fill from, `[0, 0]`: `ValueError` with the step-2 words.
6. Two groups, `[0, 333, 512, 0, 322]` (20): 333 on the west group, 322
   on the east (8 filled); and every case of 1 to 6 gives the same code per
   triangle with its triangles and vertices permuted (seed 0).
7. `replace`, a water neighbour is passed over: heath on [0,100]², lake
   512 [20,80]², river 511 [80,90] x [10,90], NVE lake [20,78] x [20,80]:
   the piece [78,80] x [20,80] (120 m²) shares 60 m with the river and 4 m
   with heath; today it takes 511, after it takes 322.
8. `replace`, an island is left out: 107's case 5 (ring [30,70]² less
   [45,55]²): heath 8,400 m², lake 1,500, the island's 100 m² uncovered
   (today heath 8,500). `test_preferred.py`'s two island tests and
   `test_land_cover_map.py::test_a_map_cut_by_prefer_covers_the_region_once`
   are amended to this (their region covered less the island in the
   base's water, 100 m²).
9. End to end, `test_cli_mesh_prefer.py`'s scene (its left-out water is
   code 0 today), with `--land-cover-partial`, since its
   `cover_polygons()` leave 1,905,280.4 m² of its domain uncovered (test
   17): `.gpkg`, `.vtk` and `.ply` hold no triangle with code 0,
   `land_cover_filled` is 1 exactly where today's code is 0, the record's
   `land_cover_filled` equals that count, and no `land_cover_fractions`
   row has code 0; the stderr line is printed once.
10. A refusal: a land cover whose only polygon lies outside the domain
    gives exit 2 naming `--land-cover` and section 3.6's words, one part,
    the domain's area (the red premise, checked first: today it exits 0
    with every triangle code 0). 10b: the same with
    `--land-cover-partial` gives exit 2 and the step-2 words.
11. On Ola's data (skipped without it), Bygdin with `clc2018-clean-25833`
    and `--prefer lake=` NVE: no triangle code 0; 145 marked, 0.1078 km²
    (±0.0001); every triangle whose centroid lies in an NVE lake farther
    than `edges_max_shift_m` (50 m) from its shore carries 512.
12. `merge_zero_share`: rows `{0: 0.25, 512: 0.75}` and `{0: 0.4, 322:
    0.6}`, both into 322, give `{322: 0.25, 512: 0.75}` and `{322: 1.0}`
    (float32, offsets rebuilt); a table with no 0 share comes back equal.

Coverage gaps (section 3.6). Figures of 13 to 16 measured with shapely
2.2.0 (scratch, this session); of 17 and 18 with the test suite's own
fixtures (`fraction_fixtures.domain_ring`, `cover_polygons()`).

13. `coverage_gap`, domain [0,100]², covers [0,100] x [0,60] and [0,60]
    x [60,100]: one part, 1,600 m² exactly, at (80, 80) ±0.01 (inscribed
    radius 20); `area` 1,600.
14. Not gaps: 107's case 5 (heath, ring [0,100]² less [30,70]²; lake 512
    [30,70]²; `--prefer` ring [30,70]² less [45,55]²): no part, `area` 0,
    though the island [45,55]² is left out of the map (test 8). Test 7's
    scene: no part either.
15. Noise: covers [0,50] x [0,100] and [50.001,100] x [0,100] (a slit of
    0.1 m², radius 0.0005 m): no part. With [50.1,100] (10 m², radius
    0.05 m): one part, 10 m² (±1e-9).
16. Tiles: domain [0,12000] x [0,1000], covers [0,4000] and [6000,12000]
    across: one part (not two, though `GAP_TILE` cuts it at 5,000),
    2,000,000 m², radius 500; `threads=1` and `threads=4` give equal
    `Gap`s, and so does the covers' order reversed.
17. CLI refusal, `test_cli_mesh_prefer.py`'s scene without the flag
    (`--land-cover` `cover_polygons()`, domain `domain_ring()`, EPSG:31983):
    exit 2, "Invalid value for --land-cover", 2 parts listed, 1.9053 km²
    in all, the larger 1.8751 km², each line ending "in EPSG:31983", no
    "covered by", and no output file (the red premise: today exit 0).
18. The override: the same with `--land-cover-partial`: exit 0, no
    triangle code 0, the stderr gap line once, record
    `land_cover_partial` "on", `land_cover_gaps` 2, `land_cover_gap_m2`
    1,905,280.4 (±0.5); without land cover the three keys are absent;
    with full cover `land_cover_gaps` 0 and "off".
19. An island fills without the flag: domain `domain_ring()`; heath 322
    over the DEM's box less L, lake 512 on L =
    `box(x0+800.3, y0-1600.7, x0+1600.9, y0-800.1)`, `--prefer lake=` L
    less H = `box(x0+1100.2, y0-1300.4, x0+1300.6, y0-1100.8)`
    (39,999.84 m², `(x0, y0)` = `dem_origin()`): exit 0, record
    `land_cover_gap_m2` 0, every triangle whose centroid lies in H has
    code 322 and `land_cover_filled` 1, every other has 0 marked.
20. The hint: a `$RASPUTIN_DATA/sources.toml` in `tmp_path` with local
    corine-class entries `east` and `east-clean-31983` (each a GeoJSON
    of the gap of test 17 as one 211 polygon) and `other` (the same
    polygon, classes `ar50`): test 17's run names
    `east-clean-31983` on both lines, never `east` or `other`.
21. Raster: test_cli_mesh_land_cover_fractions.py's cached class raster,
    whose class-0 cells the mesh reaches (its fixture says so): without
    the flag exit 2 naming `--land-cover`; with it exit 0 and no
    `land_cover_fractions` row with code 0.

**7a. Existing tests that change.** A prototype (scratch `fillplug.py`:
`_islands` returning nothing, the section 3.2 fill on the CLI's codes, no
shares merge) fails 27 of master's 11,750 tests, all of which pass on
master: `test_preferred.py` (2) and `test_land_cover_map.py` (1), test 8
above; and 24 CLI tests whose scenes leave land uncovered, in
`test_cli_mesh_landcover.py` (11), `test_cli_mesh_edges.py` (4),
`test_cli_mesh_gpkg.py` (2), `test_cli_mesh_land_cover_fractions.py` (2),
`test_cli_mesh_multi_features.py` (2), `test_cli_land_cover_polygons.py`,
`test_cli_mesh_catchments.py` and `test_cli_mesh_simplify_domain.py` (1
each). `@tester` amends each in the red commit to the filled codes, says
per test what changed, and keeps each test's own subject (a test about
holes still checks the hole's border). The shares merge may move more;
those are amended the same way and listed.

The refusal reaches further: every CLI test whose land cover leaves a
gap now exits 2, whatever its subject. Not counted (no test run, disk);
a prediction: every test on `cover_polygons()` (a 1.9 km² gap) and on the
fractions fixture's raster, so dozens across the 199 tests of the 28
files passing `--land-cover`. `@tester` runs the suite on the red commit
first, and adds `--land-cover-partial` to each refused test whose
subject is not coverage (its codes are then the fill's, as above); a
test about coverage gets its own case instead. The list goes in the red
commit's message.

## 8. Acceptance

No `@perf` run (the Status line). Before review, on the green branch:
Numedalslågen and Lagan (`tools/bench_quick.py`) give master's vertices
and triangles, and master's codes wherever master's code is not 0 (not
measured whether they hold any 0; the record's `land_cover_filled` says);
Bygdin with `clc2018-clean-25833 --prefer lake=` NVE gives master's
vertices and triangle set (corners sorted), codes equal but on the 145
island triangles, of which the prototype puts 3 on another class than
`_islands`, and record `land_cover_gaps` 0; Bygdin with
`ar50-innlandet-clean-25833 --prefer lake=` NVE exits 2 before the mesh
phase starts, listing 5 parts, 7.8634 km² in all, the largest 5.8199 km²
around 136216 6830714 in EPSG:25833 (±1 m), each of the five naming
`ar50-vestland-clean-25833` and no other id (prediction from section 2's
probe); the same with `--land-cover-partial` prints the gap line and the
fill line with 2,623 triangles, 8.006 km² and 92 groups (the
prototype's figures: the county gap plus 87 islands), its largest
5.82 km² (±0.01), and records `land_cover_gap_m2` 7,863,437 (±1); with
`ar50-vestland-clean-25833` added it exits 0 with `land_cover_gaps` 0.
Numedalslågen's quick check gives `land_cover_gaps` 0, and its
`land cover read` phase grows by under 1 s (prediction; CORINE file). `@reviewer` reads the outputs; this runs Ola's data, so
it belongs here, not in review.

## 9. Left out

- Adding NVE's lakes when a run does not name them, and `--reservoirs`
  deciding classes: the overrule stays `--prefer lake=`.
- The shore band (143 and 227 triangles on Bygdin): placement by the
  border clean-up, set by `--edges-max-shift`.
- Closing slits inside one product (ideas "20c-3"): a slit wider than
  `GAP_NOISE` is now refused, and closing it stays the clean-up's job.
- `--land-cover-partial` on the MCP mesh tool (question 1).

## 10. Size

Estimate, counted lines (`tools/count_loc.py`): `fill_no_class` and
`Filled` 45; `preferred.py` −17 (`_islands`) +4; the shares merge 15;
`cli.py` (three call sites, stderr, record, refusal, words) 30; writers
and reader (`vtk_legacy`, `ply`, `geopackage_write`, `vtk_read`) 30; run
record keys 3. The gap (section 3.6): `coverage_gap`, `Gap`,
`CoverageGap` and the two constants 45; `_covering` 25; the flag, the
refusal words, the gap line and three record keys in `cli.py` 30; the
raster refusal and its groups 15; record keys 3: about +118. **Net
about +230**. One PR, far under 700. Tests are not counted.

## 11. Questions for Ola

Answered 2026-11-01: a coverage gap is an error; a manual override
continues (section 3.6).

1. Should the MCP mesh tool also take the override? Default: no; the
   override stays a person's flag on the command line, and an MCP agent
   gets the refusal naming the file to add.
