# Increment 107: lake and glacier from a preferred source

**Status: designed** by `@architect`, 2026-10-26, on master `752b7601`,
worktree `prefer-classes`. One PR, full loop (about 140 net production lines),
no C++, no refine or mesh code, so no `@perf`. Number 107: 106 is the highest
under `docs/increments/` over every ref (`git for-each-ref` + `git ls-tree`).

Ola, 2026-10-26: "So, can we specify different polygons from different data
sources yet? Like "use Corine, but not for lakes and glaciers. There you should
prefer xyz"?", then "Yes, do that. I think that if you --perfer, that
disambigous this issue?". On case (b) below: "I think we should have warnings if
we get less water from this. It shouldn't be that hard to build a heuristic to
choose land cover types in this case. But it would mostly signal a fault."

## 0. The plain picture

Today `mesh --land-cover clc2018` takes every class from CORINE, whose 25 ha
minimum mapping unit drops small lakes: Sula (NVE station 73.27.0, 30.44 km²)
gets 5.82 % lake against NVE's lake database's 8.94 % (NVE's station record says
9.5 %). After this increment:

    rasputin mesh ... --land-cover clc2018 --prefer lake=nve-lakes --prefer glacier=nve-glaciers

takes lakes and glaciers from NVE and everything else from CORINE. The MCP
server builds its `mesh` tool from the CLI's options, so `prefer` is an MCP
parameter (a list) with no extra work.

## 1. Prior art: legacy and literature

*Legacy.* `git grep -n -i -E "glacier|lake" legacy-archive -- legacy`: 6
files; `legacy/rasputin/gml_repository.py` and `globcov_repository.py` map one
source's classes to a lake or terrain material, `extract_lakes`
(`legacy/rasputin/triangulate_dem.h`) splits faces by that one source. One
source per mesh; nothing to carry across.

*Literature.* Replacing a coarse map's water class with a dedicated water
product is established practice: ESA's CCI land cover ships a separate water
bodies map (150 m, from Envisat ASAR and other water products) "consistent
with the Water Bodies class" of its annual maps
([fedeo record](https://fedeo.ceos.org/series/eo:platform/Envisat/7e139108035142a9a1ddd96abcdfff36)).
The fill rule of (b) is the generalisation rule CORINE and AR50 apply
themselves: an area too small to keep is merged into a neighbour (AR50's own
words, in the `ar50-*` entries of Ola's `~/.config/rasputin/sources.toml`: "Areas under
15 decares are merged into neighbours"). No novelty is claimed; none searched.

*Ideas file.* Overlaps the parked NMD row in `ROADMAP.md` (a finer land cover
for Sweden): `--prefer` lets a fine source supply chosen classes without
replacing CORINE. Takes up no ideas-file entry.

## 2. What "prefer" means where the sources disagree

The key is a cell category (`lake`, `glacier`, as in the legends'
`category`). Each legend names, in a new `[prefer]` table, the one class a
preferred polygon takes and the only class it replaces: CORINE `lake = 512`
(water bodies), `glacier = 335`; AR50 `lake = 81`, `glacier = 70`. So CORINE's
rivers (511), lagoons (521), estuaries (522) and sea (523) stay CORINE's, as
NVE's lake database holds none of them. A legend without the key is refused,
naming the legend.

- **(a) preferred lake where CORINE has heath: lake.** Preferred polygons are
  cut out of every CORINE polygon, whatever its class. At Sula 3.49 % of the
  catchment moves this way, mostly from sparse vegetation (2.44 %) and heath
  (1.02 %).
- **(b) CORINE water where the preferred source has none.** CORINE's 512 is
  dropped wholly ("not for lakes"). A leftover piece is one connected part
  (`shapely.get_parts`) of what is left of a replaced polygon after the cut.
  Each takes the class whose kept CORINE polygons share the longest border
  with it, lengths summed per class (a tie: the lower code). A piece with no
  kept neighbour (it borders only preferred polygons) is left out: its area is
  uncovered, which `assign` already gives class 0 and the cell table counts as
  unspecified. Why: it is what CORINE's own generalisation does with an
  area it cannot keep, and the shoreline CORINE drew too wide borders the
  land it swallowed. Then a **warning** (a finding from the data, so allowed)
  when, inside the mesh's triangles, the preferred source gives less of a category than the
  replaced class did, by more than 1e-9 of it (rounding):
  `lake: nve-lakes gives <A> km² where clc2018 gave <B> km²; <C> km² of clc2018 water had no nve-lakes lake and took a neighbouring class`.
  Every run with `--prefer` records the three areas per key. Sula: at most 0.37 %
  refilled (the probe counts all CORINE water, not 512 alone), no warning. Nigardsbrevatn (76.5.0, 65.23 km²): glacier 76.19 %
  CORINE against 70.17 % NVE, so the warning fires there (the glacier has
  shrunk; CORINE is generalised).
- **(c) two `--prefer` sources overlap: the first on the command line wins**,
  the later one cut out of it. Order is visible and deterministic; the
  overlap is small where measured (NVE lake inside NVE glacier: 0.004 % of
  Nigardsbrevatn). The same key twice unions the two sources (Norway's and
  Sweden's lakes in one border basin).
- **(d) `--reservoirs`** is unchanged: it splits whatever lake share the cells
  have, now the preferred one. NVE's reservoirs name NVE lakes (`vatnLnr`), so
  the two now match better.

`--prefer` with a raster `--land-cover` (MapBiomas, a GeoTIFF) is refused, as
is `--prefer` without `--land-cover`; rasters are later work if asked.

## 3. Sources for Norway

| key | source | state |
|---|---|---|
| lake | `nve-lakes` (NVE Innsjødatabase, `Innsjodatabase2/MapServer/5`, lakes over 2,500 m², NLOD) | shipped catalogue entry, `rasputin fetch nve-lakes` |
| glacier | `nve-glaciers`: NVE `Bre2/MapServer/9` (`Breflate`, glacier units, 6,743 polygons, EPSG:25833) | **needs an entry** in `src_python/tin_engine/data/sources.toml`; licence not stated by the service (question 1) |
| glacier | Ola's file copy `rasputin_data/nve_glaciers/GlacierAreaOutline_NO_2018_2019_N` | shapefile, which rasputin does not read; CC BY-NC-SA 4.0 |
| glacier | RGI 7 Scandinavia (`rasputin_data/rgi7/`) | zipped shapefile; not readable as is |
| lake | `nve_hrd/lakes.geojson` | 81 gauge lakes only; not a lake map |

## 4. Blueprint

```
read_class_polygons(land cover)  ──┐
read_polygons(each --prefer src) ──┼─> preferred.overlay(base, chosen, within) -> Overlaid
legend.prefer {category: code}  ───┘        │ polygons [(code, geom)], areas per key
                                            ├─> polygon_fractions.assign (unchanged)
                                            └─> warning + run-record entries (CLI)
```

- `tin_engine/preferred.py`, pure (shapely in, shapely out, no file, no CRS):
  `overlay(base: Sequence[tuple[int, BaseGeometry]], chosen: Sequence[tuple[int, Sequence[BaseGeometry]]], corners: Float64) -> Overlaid`,
  where `chosen` is (replaced code, polygons) in command-line order,
  `corners` the mesh's triangles (T, 3, 2) and `Overlaid` a frozen
  dataclass: the polygons, and per code the replaced, given and refilled
  areas inside the triangles. No outline is needed: the polygon route has
  only the mesh's points, and the triangles are the mesh. About 75 lines.
  Steps, as the timing probe runs them (it tests against the domain outline
  where the code tests against the triangles):
  1. keep the preferred polygons that meet a triangle (STRtree of the
     triangles, `predicate="intersects"`); unclipped, as the base is;
  2. cut each base polygon by the union of the preferred polygons it meets;
  3. leftover pieces of the replaced class that meet a triangle;
  4. **the shared border:** the kept polygons near the pieces subdivided by
     `polygon_fractions.subdivide` (64 coordinates, as `assign` does), each
     piece's boundary buffered by **0.01 m**, and the length of that band's
     intersection with each subdivided piece's boundary summed per class.
     Why a band: halving inserts points on the cut lines, so a border both
     sides share is only nearly equal after subdivision; exact GEOS
     intersection then finds points, length 0. Measured (`prefer_timing.py`):
     exact finds a neighbour for 9,110 of 10,344 pieces in Numedalslågen and
     3,965 of 4,126 in Lagan; 1 mm, 1 cm and 1 m all find 9,263 and 3,987.
     1 cm sits in the middle of that plateau, far above coordinate rounding
     in projected metres (about 1e-9 m at 7e6 m) and far below CORINE's
     100 m minimum width, so a polygon that only touches a corner adds at
     most about 2 cm. Checked on those two catchments (projected metres,
     6,441 km² the larger);
  5. areas: replaced and refilled from `cells.overlay(corners, polygons)`
     times the triangle areas, given from the class fractions.
- Legend: `[prefer]` table and `Legend.prefer: Mapping[Category, int]`, about
  10 lines; TOML data uncounted.
- CLI: repeatable `--prefer KEY=SOURCE` (catalogue name or file, fetched and
  read as `--reservoirs` is), checks, warning, record entries: about 60 lines.
- `--prefer` needs no `.gpkg` output: it changes the class fractions, which
  every mesh format carries; only the cell table needs `.gpkg`.
- Each preferred catalogue source's credit, licence note and citation go into
  the mesh metadata beside `land_cover_credit`, as
  `land_cover_prefer_credit` and the like, one `key: text` part per key;
  a file gives its name only, as `--land-cover` does.

## 5. Tests (synthetic squares, metres, one CRS)

1. Base 100×100 m class 322 with a 512 square [40,60]²; preferred lake
   [30,50]×[40,60]: lake 400 m² as 512, the leftover 200 m² becomes 322, 322
   totals 9,600 m²; no warning (400 is not less than 400).
2. Base left half 322, right half 312, a 512 square [44,60]×[40,60] across
   the line; no preferred lake: the lake borders 322 for 32 m and 312 for
   40 m, so all 320 m² become 312, and the warning fires (0 < 320).
3. Preferred lake [0,20]² and glacier [10,30]²: lake first gives lake 400,
   glacier 300 m²; glacier first, glacier 400, lake 300.
4. Base 100×100 m class 322 with a 511 river strip [0,100]×[45,55]
   (1,000 m²; 322 is 9,000 m²); preferred lake [40,60]²: lake 400 m² as
   512, 511 800 m², 322 8,800 m²; no 512 in the base, so no leftover and no
   warning (400 is not less than 0).
5. A piece bordering only a preferred lake: base 322 square with a 512
   square [40,60]², preferred lake the ring [30,70]² less [45,55]²: the
   100 m² piece [45,55]² has no kept neighbour and is left out (class 0,
   100 m²).
6. Refusals: unknown key, legend without `[prefer]`, raster `--land-cover`,
   `--prefer` alone, `KEY=` with no source.
7. CLI end to end: the cell table's lake share comes from the preferred
   polygons, `--reservoirs` splits it, class fractions still sum to 1.

## 6. Size and timing

About 145 net production lines (`tools/count_loc.py`), far under 700. One
PR, buildable Tuesday 2026-10-27 for Wednesday: lakes need nothing new; the
glacier entry is one catalogue entry beside it. Probes:
`docs/increments/107-probes/prefer_probe.py` (shares),
`docs/increments/107-probes/prefer_timing.py` (timing; Apple M1 Max on
battery, one run each, no mesh build).

Expected added time of `--prefer lake=...` on the quick check's two
land-cover catchments (the quick check's own cases pass neither option, so
its timings do not change):

| step | Numedalslågen, NVE lakes | Lagan, SMHI lakes |
|---|---|---|
| read the preferred source | 2.9 s (26,538 lakes in the box) | 3.1 s (42,176) |
| cut (lakes meeting the domain) | 3.7 s (5,550 lakes, 900 polygons) | 1.8 s (530, 724) |
| neighbour search | 2.7 s (10,344 pieces) | 3.0 s (4,126) |
| area accounting, two `cells.overlay` calls | about 1 s (105b: 0.45 s each at 86,537 triangles) | about 1 s |
| `assign`'s input | 201,118 pieces become 249,788 (+24 %) | 636,492 become 635,241 (0 %) |
| **total** | **about 10 s, plus up to 24 % of the land-cover fractions phase** | **about 9 s** |

The fractions phase's own time on these catchments is not measured (it needs
a mesh); `@developer` reports it, with and without `--prefer`, in the green
handback. The first prototype, which cut by the union of every lake in the
box and searched neighbours polygon by polygon, took 68 s for the search
alone: the steps above are pinned because they are the measured ones.

## 7. Question for Ola

1. NVE's glacier map service states no licence, and the 2018–2019 archive copy
   is CC BY-NC-SA (non-commercial). Default: add `nve-glaciers` with that
   stated in its licence note, and the Wednesday demo uses `--prefer` for
   lakes only unless you confirm the licence.
