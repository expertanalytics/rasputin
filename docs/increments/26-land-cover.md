# Increment 26: land cover for the São Francisco basin — water bodies as constraints, class fractions per triangle

Status: **designed by `@architect`, 2026-10-03, at Ola's request; not
ruled.** Written before `@tester` per `docs/increments/README.md` step 1.
Nothing is built. The questions for Ola are at the end, each with a
recommended default. Six pull requests, listed under "The PR split".

**Why "26".** 25 is plain output (`docs/increments/25-plain-output.md`, on
its own branch); 26 is the first free number on every branch
(`git log --all --name-only --format='' -- 'docs/increments/2[5-9]*'` lists
only `25-plain-output.md`). In the basin order of work this is the
"basin's own inputs" item (`docs/increments/23-basin-scale.md`, "Order of
work", item 12), the land-cover half of it; sub-catchments and rivers from
the DEM are the other half and get their own design.

**Names used here.** Earlier increments, by number, and what they did:
10, mesh output (its LOC overrun is the +39 % worst case used for
estimates); 11, decoding GeoTIFFs; 16b, polygons and lines from
GeoPackage or GeoJSON as constraint lines, with class maps naming what each
line is; 16c, a land-cover code per triangle (regions between constraint
lines, one point-in-polygon test per region); 16e, several feature sources
in one mesh, one code system per mesh; 21, parallel refinement and Ola's
determinism ruling; 22, a catchment from the DEM, its traced outline reduced
by area-preserving segment collapse; 23, basin scale: 23a-1 reads DEM blocks
from the local tile cache, 23a-2 is `rasputin fetch`, which fills the cache,
23b to 23g cut a large domain into pieces, mesh them in parallel and stitch
them with the seams cleaned up; 25, plain output (what a run writes, in
named fields, `--stats` and `--record`); 15f-3, the last part of the edge
strip. This increment's own pull requests are 26a to 26e. **BHO** is Brazil's
official coding of drainage basins; a "level-3 unit" is one of the nine
sub-basins of the São Francisco at its third level, used here only as
measurement domains (Ola's ruling on BHO).

## In short

MapBiomas maps Brazil's land cover as a 30 m raster of class numbers, one
file per year. This increment puts it into the mesh in two ways, as Ola ruled
on 2026-10-01:

1. **Water bodies become constraint lines.** Every lake or reservoir larger
   than a minimum area is traced along the cell edges, moved into the mesh's
   coordinates, and simplified without changing its area (the method of
   increment 22's catchment outline, extended to many rings, islands, and
   lines it must not cross). The polygons go into `rasputin mesh` as
   features, like CORINE's do today.
2. **Every other class becomes a list of fractions per triangle**, for
   example "62 % soybean, 31 % savanna, 7 % pasture". The fractions are the
   exact areas of the 30 m cells inside each triangle. Fractions too small to
   matter are dropped, and the area they held is not lost: it is carried to
   the next triangles along a space-filling curve and placed there, which
   is Floyd-Steinberg error diffusion done on class areas (Ola's "local
   out-of-balance ledger"). Class totals over the whole mesh then stay exact
   to about a hectare, and those of any compact region within a bound the
   ledger itself reports.

Measured on real basin meshes for this design (the probes are in
`docs/increments/26-probes/`, results under "What was measured"):

- dropping small fractions the way land-surface models do (renormalising
  inside each triangle) loses 8 to 45 % of the open water and 1.4 to 11 %
  of the coffee in the Rio Corrente sub-basin, at a 5 % cutoff; the ledger
  keeps every class's total within 0.05 %;
- the cutoff must be **relative and absolute together**: an entry is dropped
  only if it is under 5 % of its triangle **and** under 1 ha. A purely
  relative cutoff drops whole fields from the large flat triangles (the
  largest is 67 km² at 20 m), the ledger then carries square kilometres
  across the map, and local errors are as bad as renormalising;
- MapBiomas Collection 11 (August 2026, years 1985 to 2025) is on a public
  bucket as one tiled GeoTIFF per year; the tiles meeting the basin's
  outline are 72 MB per year (177 MB for its whole box), against 1.71 GiB of
  ANADEM.

## Ola's rulings this builds on

All from `docs/research/raster-to-vector.md` ("Open questions for Ola" and "A
further ruling, on the fractions themselves"), 2026-10-01, unless marked:

- **The hybrid.** "Water bodies and rivers are constraints; every other class
  is carried as a fraction per triangle."
- **A minimum area, no merging.** "A water body is a constraint only above a
  minimum area tied to the tolerance (a few times tolerance²); a smaller one
  stays a fraction of the triangles it falls in, and nothing is merged into a
  neighbour."
- **Rivers from the DEM.** BHO is dropped as a geometry source ("I'm not
  interested in archaic maps", `23-basin-scale.md`, ruling on BHO); river
  lines will come from DEM-derived drainage, designed with sub-catchments.
  Until then a river line is whatever `--features` gives.
- **Exact area.** "Each lake is traced on cell edges, reprojected to metres,
  and reduced by increment 22's area-preserving segment collapse, with a
  no-crossing check against the other lakes, the river lines and the domain
  outline, plus the swept-region test against their vertices." River-lake
  crossings are kept as nodes on both lines; only new crossings are errors.
- **8-connected water, free pinch points.** Water cells touching at a corner
  are one water body; the pinch may open in the reduction, the two sides may
  never cross.
- **The legend.** "We need to keep high resolution on vegetation types and
  crop farming types in Brazil": MapBiomas's full legend, crop types
  included; "year as an option, default latest". The collection was not
  ruled (question 1 below).
- **The cutoff and the ledger.** "We could even have a cutoff on the
  fractions. 0.1% soybean does not carry so much information." "95% corn,
  5% soybean _could_ become 100% corn." "A local out-of-balance ledger,
  trying to compensate for missing covers", kept by area; "the general idea,
  Floyd-Steinberg dithering, applied to class areas instead of pixel
  intensities, and related publications should be used to resolve this."
- **Determinism** (increment 21's ruling): same input, same output, for any
  thread count. And from the basin design: the mesh does not depend on the
  machine.

## What was measured for this design

By `@architect`, 2026-10-03, on AC power, Apple M1 Max, with the probes in
`docs/increments/26-probes/` (throwaway, not production code; each file says
how to run it). MapBiomas Collection 11, year 2025, read straight from its
public bucket by range requests.

### The MapBiomas file

`fetch_size.py` and a `HEAD` request, against
`https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection11/lulc/coverage/brazil_coverage/brazil_coverage-col11_2025.tif`:

- 763,286,575 bytes, `Accept-Ranges: bytes`, `Last-Modified: Tue, 18 Aug
  2026 21:38:54 GMT`, and a real `ETag` (an MD5-style digest, unlike
  ANADEM's placeholder);
- one page, no overviews: 146,501 rows × 154,470 columns of `uint8`, LZW
  with the horizontal predictor (so it needs the `codecs` extra,
  `imagecodecs`, which `pyproject.toml` already offers), tiles of 256 × 256;
- 346,092 tiles, 194,841 of them empty (byte count 0: outside Brazil); the
  header and tile offsets end at byte 4,153,648, so 23a-2's header rule
  (1 MiB, doubled until complete) reads 8 MiB;
- EPSG:4326, `PixelIsArea` (a value is the class of the whole cell), cell
  0.00026949458523585647° (30.0 m at the equator), upper-left corner at
  (-74.02099974839176°, 5.42303953870114°); no `GDAL_NODATA` tag: 0 means
  no data by MapBiomas's convention;
- three tiles over the basin (Três Marias, Sobradinho, western Bahia) were
  range-read and decoded alone by tifffile with imagecodecs, as 23a-1 decodes
  ANADEM blocks. Classes seen: 3, 4, 9, 11, 12, 15, 21, 24, 25, 30, 33, 39,
  41, 48, 62, 75;
- **cell area** (WGS 84 ellipsoid, `pyproj.Geod`): 887.5 m² at 7° S, 864.3 m²
  at 15° S, 836.0 m² at 21° S: the basin's range. "A 30 m cell" below means
  about 864 m².

What a fetch would cost (`fetch_size.py`, tiles meeting the box of each
outline, BHO outlines as measurement domains as 23's ruling on them allows):

| domain | tiles in its box | bytes | empty tiles | tiles meeting the outline | bytes |
|---|---:|---:|---:|---:|---:|
| the basin (BHO level 2) | 32,835 | 177 MB | 5,640 | 11,714 | 72 MB |
| Rio das Velhas piece | 540 | 5.1 MB | 0 | 265 | 2.5 MB |
| unit 761 (lower basin, Sobradinho, Itaparica) | 8,658 | 42 MB | 438 | 3,893 | 20 MB |
| unit 769 (upper basin, Três Marias) | 3,243 | 30 MB | 0 | 2,078 | 19 MB |

Per year. ANADEM over the basin's outline is 1.71 GiB (23's measurement).

### Fractions on real basin meshes

`fractions_probe.py` on the BHO level-3 unit 764, the Rio Corrente in
western Bahia (34,243 km², soybean and cotton on the plateau, savanna
elsewhere), meshed from ANADEM at 20, 10, 5 and 2 m by `@perf` on
2026-10-02 (`docs/benchmarks/2026-10-02/basin-level3/README.md`; the files
are in `../rasputin_data/sao_francisco_piece/meshes/level3/`). For each
triangle, the exact area of every 30 m cell inside it, by class; then the
cutoff and the ledger as designed below, with variants. The exact areas are
the reference; every error below is against them.

The meshes:

| tolerance | triangles | mean triangle | in 30 m cells | median / 90th pct / 99th pct / largest | classes per triangle, exact: mean (max) |
|---:|---:|---:|---:|---|---:|
| 20 m | 115,578 | 29.7 ha | 342 | 4.7 ha / 63 ha / 3.9 km² / 67 km² | 1.96 (11) |
| 10 m | 275,421 | 12.5 ha | 143 | 2.5 ha / 26 ha / 1.4 km² / 67 km² | 1.79 (10) |
| 5 m | 682,531 | 5.0 ha | 58 | 1.1 ha / 10 ha / 0.56 km² / 24 km² | 1.60 (10) |
| 2 m | 2,329,553 | 1.5 ha | 17 | 0.32 ha / 2.8 ha / 16 ha / 9.2 km² | 1.40 (10) |

Triangle area spans four to five orders of magnitude inside one mesh. Where a
crop is the largest class, the median triangle is 2.4 ha at 20 m and 2.8 ha
at 2 m, against 5.0 ha and 0.30 ha elsewhere: at 5 m and finer the flat
cropland gets much larger triangles than the rest, as the research note
guessed; at 20 m it does not.

Three ways of dropping small entries, all at a 5 % cutoff. "Renormalise" is
what CTSM, WRF's Noah mosaic and SWAT do: drop the small entries and scale the
rest up inside the triangle. "Ledger" is this design (the rule "present",
below). "Misplaced" is the share of a square's area whose class is wrong,
half the sum over classes of |written − exact| divided by the square's area,
over squares fully inside the domain; the 95th percentile over the squares.

| tolerance | method | water total | coffee total | worst class total | misplaced, 5 km squares | misplaced, 25 km squares | largest ledger | classes per triangle |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 20 m | renormalise | −45.4 % | −10.7 % | −45.4 % (water) | 5.0 % | 2.6 % | | 1.60 |
| 20 m | ledger, 5 % only | +0.05 % | +0.04 % | +0.05 % | 5.3 % | 0.78 % | 588 ha | 1.48 |
| 20 m | **ledger, 5 % and 1 ha** | 0.00 % | 0.00 % | 0.00 % | **0.15 %** | **0.013 %** | **4.4 ha** | 1.65 |
| 10 m | renormalise | −28.6 % | −8.3 % | −28.6 % (water) | 3.6 % | 1.7 % | | 1.52 |
| 10 m | ledger, 5 % only | +0.02 % | +0.01 % | −0.09 % | 3.2 % | 0.42 % | 535 ha | 1.47 |
| 10 m | **ledger, 5 % and 1 ha** | 0.00 % | +0.01 % | +0.01 % | **0.17 %** | **0.016 %** | **4.5 ha** | 1.50 |
| 5 m | renormalise | −18.4 % | −3.2 % | −18.4 % (water) | 2.4 % | 1.2 % | | 1.42 |
| 5 m | ledger, 5 % only | 0.00 % | +0.01 % | −0.04 % | 1.8 % | 0.20 % | 182 ha | 1.27 |
| 5 m | **ledger, 5 % and 1 ha** | 0.00 % | 0.00 % | −0.04 % | **0.16 %** | **0.011 %** | **4.6 ha** | 1.36 |
| 2 m | renormalise | −7.9 % | −1.4 % | −8.5 % (other non-vegetated) | 1.3 % | 0.68 % | | 1.29 |
| 2 m | ledger, 5 % only | 0.00 % | 0.00 % | +0.03 % | 0.83 % | 0.09 % | 68 ha | 1.16 |
| 2 m | **ledger, 5 % and 1 ha** | +0.01 % | +0.01 % | +0.02 % | **0.14 %** | **0.013 %** | **4.7 ha** | 1.14 |

What this shows:

1. **Renormalising loses scattered classes.** Open water in the Corrente is
   river channels and small ponds, below 5 % of almost every triangle it
   touches, so it vanishes. Coffee loses 1.4 to 10.7 %; soybean, which
   dominates where it occurs, gains 1.1 to 3.4 %. This is the systematic bias
   Saura 2002 and Moody and Woodcock 1995 describe, measured. Adding the 1 ha
   floor to renormalisation cuts the loss by half or more but does not remove
   it (water −4 to −12 %).
2. **The ledger keeps every class's total**: the worst class over all four
   tolerances is within 0.05 % with the floor and 0.09 % without it, and the
   remainder left in the ledger at the end of the curve is at most 0.95 ha.
3. **A purely relative cutoff makes the ledger travel.** 5 % of a 67 km²
   triangle is 3.3 km², a real field. The ledger then holds up to 588 ha and
   places it where the curve next finds the class, so 5 km squares are as
   wrong as with renormalisation. With the 1 ha floor no single drop exceeds
   1 ha, the ledger stays at 4.4 to 4.7 ha **at every tolerance**, and the
   5 km squares are 0.14 to 0.17 % misplaced: a ninth to a thirtieth of
   plain renormalisation's, and a half to a third of renormalisation with the
   same floor (0.32 to 0.45 %, in `summarise.py`'s full tables). The ledger's
   size is set by the floor, not by the mesh.
4. **The cutoff saves little storage**: 1.96 → 1.65 entries per triangle at
   20 m, 1.40 → 1.14 at 2 m. It is worth having for what Ola asked it for
   (no 0.1 % soybean), not for size.
5. **Capping the number of classes per triangle** (as WRF's mosaic keeps the
   three largest) was tried at 4 and 3 with the same rules: it brings the
   ledger back to 20 mean triangles and 5 km errors of 5 % at 20 m, because
   the large triangles must drop real fields again. Not used.
6. **"Anywhere" versus "present"** (whether the ledger may put a class into a
   triangle where the raster has none of it): "anywhere" is slightly better
   locally (5 km: 0.10 % against 0.14 % at 2 m) and leaves a smaller end
   remainder, but writes 8 to 18 % more entries (20 m to 2 m) and can put
   soybean on a triangle
   that has none; "present" never does. Ola's first proposal ("kept":
   only into triangles that already hold the class above the cutoff) is
   worse than both: with the floor, a ledger of 11.5 ha against 4.4 ha at
   20 m and 0.30 % misplaced on 5 km squares against 0.15 %. The full
   tables (cutoffs 1, 5 and 10 %; floors none, 1 ha and 10 ha; four rules,
   "kept" run at 20 and 2 m) are printed by `summarise.py`.

Time, single-threaded probe code, first runs of 2026-10-03: the exact areas
took 6.4, 7.0, 8.2 and 11.3 s at 20, 10, 5 and 2 m (the JSON files kept are
from later re-runs, 6.6 s at 20 m and 11.9 s at 2 m: run-to-run noise with
other work on the machine) (about 0.15 µs per cell covered plus 2.2 µs per
triangle); the ledger took 2.4 s for the 2.33 M triangles at 2 m (1.0 µs per
triangle, scanning all 256 class slots each time).

### Water bodies in three sub-basins

`water_probe.py`: MapBiomas class 33 (river, lake and ocean) over the box of
each BHO level-3 unit, cut into 8-connected bodies (cells touching at a
corner are one body). A body counts for the unit when any of its cells'
centres is inside the outline. (The first version of this probe used the
body's first cell in raster order, which lost Três Marias from unit 769:
that body's first cell lies north of the outline. `@reviewer` found it in
design review round 1; the figures below are from the corrected probe.)
Every body of at least 1.44 ha (four times (60 m)², the minimum area of
question 7) was traced on its cell edges, moved to the basin's Transverse
Mercator (the CRS of the level-3 meshes), its pinch points opened by 1 nm,
and reduced by increment 22's `reduce_ring` as it stands (outer rings
only; islands not reduced). Areas use one cell area for the whole unit,
868 m² (the Corrente's mean), so they are good to a few percent across the
basin's latitudes.

| unit | area | bodies | ≥ 1.44 ha | their share of the water | ≥ 10 ha | ≥ 10 km² | traced vertices | reduced, 30 m | 60 m | 120 m |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 761, lower basin (Sobradinho, Itaparica) | 209,316 km² | 36,567 | 3,457 | 98.4 % | 518 | 6 | 383,224 | 73,689 | 38,243 | 25,340 |
| 769, upper basin (Três Marias) | 106,394 km² | 30,181 | 2,836 | 94.8 % | 345 | 4 | 335,672 | 57,465 | 32,335 | 24,375 |
| 764, Rio Corrente | 34,243 km² | 1,474 | 96 | 88.8 % | 16 | 1 | 9,604 | 1,472 | 876 | 664 |

The water inside the outlines, counted by cell centre: 4,947 km² in 761,
1,543 km² in 769, 18 km² in 764. "Their share" is of the water in the
bodies that count for the unit, whole, including the parts of a body
outside the outline. "Traced vertices" counts corners where the outline
turns; "reduced" is after `reduce_ring` at the horizontal tolerance shown.

**The largest body in each unit is the river network, not a lake.** In 761
it is the São Francisco itself, with Sobradinho and Itaparica joined by the
channel: 4,425 km² of water cells, a box of roughly 450 × 670 km, 3,159
holes (islands), and an outer ring of 157,138 traced vertices (4,799 km²
with the islands inside it), reduced to 12,655 vertices at 60 m in 0.84 s.
In 769 it is Três Marias joined to the upper channel: 1,290 km², about
455 × 270 km, 620 holes, 182,794 traced vertices, 12,925 at 60 m in 1.0 s.
Even in the Corrente the largest body (19.7 km², 15 holes) is the river.
All three reach the edge of their unit's window, so each was cut there:
on the whole basin the channel very likely joins the reservoirs and the
river down to the sea into **one body with one very large ring** (the sum
of these pieces and more: of the order of half a million to a million
traced vertices and several thousand holes). Section 2.4 and 26d's tests
are sized for that.

Three more things the probe settles:

- **The minimum area is what keeps the count down.** 761 has 36,567 bodies,
  but 3,457 of them hold 98 % of its water; in 769, the 2,836 bodies of at
  least 1.44 ha hold 95 %. Scaled by area from these three units (55 % of
  the basin, and a body crossing two units is counted in both, so a rough
  figure), the basin has about 11,600 water bodies of at least 1.44 ha and
  about 130,000 constraint vertices at 60 m.
- **Area is kept.** The worst relative change over all 6,389 rings and three
  tolerances is 3.3e-11 (a 5.3 ha pond, 1.8e-6 m²). Measuring it needed
  care: the textbook shoelace as two dot products (Σ x·y' − Σ x'·y) cancels
  catastrophically on the 761 channel ring and reported a change of 7 m²
  (1.5e-9); one cross product per edge summed with `math.fsum` reports
  none. The tests must measure area the second way.
- **Pinches are common.** 1.3 % of the traced vertices in 769 and 761 are
  pinch points (one corner visited twice: 4,518 and 5,162 of them). With a
  pinch closed, increment 22's crossing check refuses any collapse whose new
  edge touches the other visit of the pinch; on the research note's
  diagonal strip that blocked every collapse. Opened by 1 nm, these rings
  reduce. The design below makes the check accept a closed pinch instead,
  so that no vertex is moved by hand and the area stays exact.

## Prior art: legacy and literature

### Legacy

```
$ git grep -il -E "mapbiomas|land.?cover|fraction|dither|globcov|water.?bod|lake" legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/rasputin/land_cover_repository.py
legacy-archive:legacy/rasputin/material_specification.py
legacy-archive:legacy/rasputin/tin_repository.py
legacy-archive:legacy/rasputin/triangulate_dem.h
legacy-archive:legacy/rasputin/web_visualize.py
legacy-archive:legacy/rasputin/wfs_repository.py
legacy-archive:legacy/tests/test_gml_repository.py
legacy-archive:legacy/tests/test_land_cover_repository.py
```

The research note read all but three of these (its "Legacy" section): the
legacy answer to a land-cover raster (GlobCover) was the class at each
triangle's centre, no constraints and no fractions; CORINE as vector GML is
16b's. The three new hits, read from the tag:

- `legacy/rasputin/triangulate_dem.h:846`, `extract_lakes`, bound in
  `legacy/bindings.cpp:371`:
  the triangles whose slope is under 0.01 are "lakes", returned as a
  separate face list. A display heuristic (flat is not water: a flat field
  is flat, a river in a gorge is not). **Not carried.**
- `material_specification.py`: a three.js material for those faces. **Not
  carried.**

`git grep -n -i -E "fraction|mapbiomas|dither" legacy-archive -- legacy`
finds no line. Nothing
on fractions, cutoffs or water tracing is carried over; `@migration-expert`
is not needed.

### Literature

Most of the prior art is collected in `docs/research/raster-to-vector.md`
(sections 1 to 5 and "Dropped small fractions"), with how each citation was
checked; it is not repeated here. What this design takes from it, what it
adds, and where it departs:

**Exact cell areas in a triangle.** Conservative remapping between grids
(Jones 1999, SCRIP; Ullrich and Taylor 2015) computes each target cell's
value from the exact overlap areas with the source cells, so integrals are
kept by construction; Taylor 2024 shows the conservation fails when cell
shapes are misrepresented. The clipping of a convex polygon by axis-parallel
lines is Sutherland and Hodgman 1974 ("Reentrant polygon clipping",
*CACM* 17(1):32-42, Crossref checked 2026-10-03). `exactextract` (Baston,
ISciences, Apache 2.0; documentation read 2026-10-03,
`isciences.github.io/exactextract/background.html`) computes the covered
fraction of every raster cell under a polygon by traversing each ring once
"making note of when it enters or exits a raster cell", and knows every
untouched cell is wholly in or out. **What differs:** our polygons are
triangles, so the clip is cheap and needs no ring traversal; rows are clipped
as strips and each strip's interior cells, wholly covered, take their full
area without a clip (exactextract's observation); and every edge-gridline
crossing is computed from the edge's two endpoints in one canonical order, so
the two triangles on an edge compute the same crossing bit for bit and the
pieces of a cell tile it to rounding. No library is used: exactextract's
library requires GEOS and its command line GDAL (its `CMakeLists.txt` at
`3b4926b`), and GDAL is prohibited (`CLAUDE.md` §2).

**The cutoff and the ledger.** Floyd and Steinberg 1976 (from a secondary
source; see the research note) is the mechanism: quantise, carry the signed
error forward to unvisited neighbours, total kept. **Error diffusion along a
space-filling curve** replaces the raster scan by a curve: Witten and Neal
1982 ("Using Peano curves for bilevel display of continuous-tone images",
*IEEE CG&A* 2(3):47-52) and Velho and Gomes 1991 ("Digital halftoning with
space filling curves", *SIGGRAPH '91*:81-90), both Crossref checked
2026-10-03; Asano 1996 partitions the curve into squares. **Vector error
diffusion** (Damera-Venkata and Evans 2001, Crossref checked) carries the
error as one vector, here over classes, summing to zero. **Stability**: Fan
1993 and Eschbach and Pedersen 2017 show multi-level error diffusion can
build large local errors. **What models do**: CTSM, WRF's Noah mosaic and
SWAT drop small tiles and renormalise inside the cell (the research note
read their source and documentation); Johnson and Clarke 2021 keep region
totals by quota. **What this design takes:** the one-dimensional carry along
a Hilbert curve, as Witten and Neal and Velho and Gomes do, because it gives
an identity no two-dimensional filter gives (the error over any stretch of
the curve is the difference of two ledger states, "What it guarantees");
and the error as a zero-sum vector of areas. **Where it departs, and why:**

- *Triangles of unequal area*, not pixels: the quantum is the triangle's
  area, so a triangle can absorb at most its own area. Nothing in the
  halftoning literature found covers irregular cells (the research note's
  search).
- *No two-dimensional weights* (Floyd-Steinberg's 7/16, 3/16, 5/16, 1/16 to
  the forward neighbours): on a mesh each triangle has one to three forward
  edge-neighbours of any size, the weights have no published rule, the
  multi-neighbour form can go unstable (Fan 1993), and it loses the
  identity above. The 1-D curve's next triangle is almost always an edge- or
  vertex-neighbour.
- *A threshold that is relative and absolute together*, which no source
  found uses: CTSM's `toosmall_*` thresholds are relative (percent of the
  grid cell), WRF's is a count. Measured above: a relative threshold alone,
  on a mesh whose triangles span five orders of magnitude, makes the ledger
  carry hundreds of hectares.
- *No cap on the number of classes* per triangle (CTSM's `n_dom_pfts`, WRF's
  `mosaic_cat`): measured above, it brings the large errors back.
- *A class is placed only where the raster has some of it* in that
  triangle. Halftoning places a dot anywhere the accumulated error says;
  that would put soybean on a triangle with none, so it is a departure, and
  measured to cost little (point 6 above).
- *The ledger's size is not bounded by a theorem.* The chairman-assignment
  bound (Tijdeman 1980) is for equal quanta and one choice per step; here it
  is measured (4.4 to 4.7 ha on the Corrente at every tolerance) and
  reported on every run, so a run where it grows is seen.

**Water bodies.** Labelling 8-connected components by runs of cells, two
scans with union-find: He, Chao and Suzuki 2008 ("A run-based two-scan
labeling algorithm", *IEEE TIP* 17(5):749-756, Crossref checked
2026-10-03); runs, not cells, are stored, which is what makes the basin's
2.15 G-cell box affordable (water is a few percent of it). Crack following
(Kovalevsky 1989; the research note's section 1) gives the outline on cell
edges. The reduction is increment 22's area-preserving segment collapse
(Kronenfeld, Stanislawski, Buttenfield and Brockmeyer 2020), extended as the
research note's section 2 argued: several rings in one check, and the
swept-region test of Saalfeld 1999 and de Berg, van Kreveld and Schirra 1998
against every vertex that must keep its side. **The pinch rule** treats a
ring that touches itself at a vertex as a *weakly simple* polygon: Chang,
Erickson and Xu, "Detecting weakly simple polygons", *SODA 2015*
(proceedings issued December 2014), pp. 1655-1670, and Akitaya, Aloupis,
Erickson and Tóth, "Recognizing weakly simple polygons", *Discrete Comput.
Geom.* 58:785-821, 2017 (both Crossref checked 2026-10-03; neither read
beyond the record). They decide whether a polygon with touching or
overlapping edges can be perturbed into a simple one, in general; our case
is the simplest one, two visits of one vertex with no shared edge, where the
condition is that the two visits' edge pairs do not interleave around the
vertex. **What differs from increment 22:** the check accepts that touch
(today it refuses it, which blocks the reduction next to every pinch), lines
that must not be crossed are added to the check, and crossings that exist in
the input (a river entering a reservoir, a lake on the domain's outline)
become fixed vertices of the ring instead of errors.

**MapBiomas.** Souza et al. 2020, *Remote Sensing* 12(17):2735 (Crossref
checked; the research note read its section 2.3.5, the spatial filter). The
Collection 11 algorithm theoretical basis document, the method's own
description (`ATBD-General-Collection-11-versao-1.pdf`
on `brasil.mapbiomas.org`, read 2026-10-03): the legend of 33 mapped
classes with their numbers (its Table 3), a spatial filter that "removes
isolated pixels", and a minimum mapping unit of "6 pixels (approximately
0.5 ha)" on transitions. The Collection 11 factsheet (version 12/08/2026,
same site, read 2026-10-03): years 1985 to 2025, four new classes, the
licence and how to cite (under "Getting MapBiomas").

### Novelty

The research note's two searches (2026-10-01) found every piece and no
error-diffusion scheme for land-cover fractions, no error diffusion of class
areas over an irregular mesh, and no land-surface or hydrological model that
compensates trimmed fractions in neighbouring cells. One more web search,
2026-10-03 ("Hilbert curve error diffusion sub-grid land cover fractions
threshold conserve area triangulated mesh hydrological model"), found nothing
new. This design adds the relative-and-absolute threshold, the "present"
rule and the prefix identity as its own choices; they are small, and **no
novelty is claimed**. Before any claim, the research note's list of checks
stands (grey literature on ORCHIDEE, SURFEX, ISBA, VIC and mHM tile
handling; LUH2; full reads of Asano et al. 2003, Tokuyama 2007 and Brunton
et al. 2015), plus a read of Velho and Gomes 1991 for any treatment of
unequal quanta.

## 1. Class fractions per triangle

### 1.1 What is computed, and in which coordinates

The input is a finished mesh (vertices in the mesh's CRS, in metres;
triangles) and the MapBiomas window covering it, for one year. The output is,
per triangle, a short list of (class, fraction) pairs whose fractions sum to
one.

**The coordinates are the raster's own.** Every mesh vertex is moved once to
longitude and latitude (pyproj through `crs.reprojector`, the one place a
transformer is built, `always_xy`), then to the raster's index space,
`col = (lon − lon₀) / s`, `row = (lat₀ − lat) / s` with `s` the cell size and
(`lon₀`, `lat₀`) the window's upper-left corner. There every cell is the unit
square `[c, c+1] × [r, r+1]`. The alternative, the mesh's CRS, makes every
cell a different quadrilateral, so every clip is a general polygon clip and
every cell corner must be reprojected (about 735 M cells under the basin
against 119 M vertices at 1 m). Rejected.

**A triangle's edges are straight in index space**, which is not quite the
image of a straight edge in the mesh's CRS. Measured (`crs_probe.py`, in UTM
23S and in the basin's Transverse Mercator, at four places in the basin): a
straight 1 km edge lies within 8.7 mm of its straight image, 5 km within
0.22 m, 20 km within 3.4 m. The two triangles on an edge use the same
straight segment, so the triangles still tile the domain's image without gap
or overlap; what the bend moves is a sliver of a cell from one triangle to
its neighbour (at most 0.22 m × 5 km, under 0.01 % of a triangle with 5 km
sides). No area is lost.

**Areas are true areas.** Each overlap piece's index-space area is
multiplied by the ellipsoidal area of a cell in its row (WGS 84, computed
once per row by `pyproj.Geod`: 887.5 m² at 7° S, 836.0 m² at 21° S). So
`a[t, k]`, the area of class `k` in triangle `t`, is in square metres on the
ellipsoid, and the exact fractions are `f[t, k] = a[t, k] / Σₖ a[t, k]`.

**The ledger works in the mesh's own measure**, the planar area `A[t]` of the
triangle in the mesh's CRS, which is what a consumer multiplies a fraction
by: the exact class areas it starts from are `f[t, k] · A[t]`. The ratio of
planar to ellipsoidal area (the projection's scale squared) changes by
under 1e-3 across 20 km in either CRS (`crs_probe.py`: at most 7.3e-4 at the
four places), so the choice changes
no fraction measurably; it decides only which total is kept exact, and the
one the user of the mesh sees is the planar one.

**Class 0** (no data in MapBiomas) is a class like any other: it gets a
fraction where it occurs. Inside units 761, 769 and 764 there is none (761's
window holds 29 M class-0 cells, every one outside its outline, at sea;
checked by cell centre against the BHO outlines).

### 1.2 The exact overlap

C++, a new header `include/terrain/land_cover/overlap.hpp`, pure, no I/O:

```
struct ClassAreas {                      // a sparse table, one row per triangle
    std::vector<std::uint32_t> offsets;  // T + 1
    std::vector<std::uint8_t>  classes;  // ascending within a row
    std::vector<double>        areas;    // m², > 0
};
ClassAreas class_areas(std::span<const Point2> corners,     // 3 per triangle, index space
                       const raster::RasterView<std::uint8_t>& classes,
                       std::span<const double> row_area);   // m² per cell, one per row
```

Per triangle: for each row it spans, clip it to the strip `r ≤ y ≤ r + 1`
(two clips by horizontal lines, Sutherland and Hodgman); in the strip, every
column wholly between the strip polygon's left and right chains is a full
cell (area 1, no clip), and only the columns its edges cross are clipped by
vertical lines. Pieces of zero area are skipped. Every crossing of a
triangle edge with a grid line is computed from the edge's endpoints taken in
one canonical order (the lexicographically smaller first), so the two
triangles on an edge get the same crossing bit for bit, and the pieces of a
cell sum to its area up to the rounding of the shoelace sums (in coordinates
local to the cell's corner, which are below 2 in magnitude). Triangles are
independent: chunks run in parallel and their rows are concatenated in chunk
order, so the result does not depend on the thread count. A triangle outside
the window is a refusal (`OutsideWindow`), never a clamp: the window is
built to hold every vertex plus one cell. `raster::RasterView<T>` is
already a template (`include/terrain/raster/view.hpp`); what is new is a
binding that hands it a `uint8` NumPy buffer without a copy.

The probe's version of this (one clip per cell, a heap allocation per clip)
ran at about 0.15 µs per covered cell plus 2.2 µs per triangle on one
thread; the full-cell spans remove most of the first term.

### 1.3 The visiting order: a Hilbert curve over the mesh

Each triangle's centroid, in the mesh's CRS, is quantised to 2³² steps on
each axis of the mesh's bounding square and given its 64-bit Hilbert index.
Triangles are visited by (index, centroid x, centroid y, then the
triangle's three vertices sorted by (x, y), compared coordinate by
coordinate). The centroids of two triangles of a valid mesh are never equal
in exact arithmetic (each lies inside its own triangle, and the triangles do
not overlap), but computed in doubles two can round to the same values; the
vertices, which are distinct triangles' input data, settle that tie. So the
order is total and a function of the mesh alone: not of the thread count,
the order of rows in the file, or how the mesh was cut into pieces. For the
same reason each triangle's corners are handed to the clipper in a fixed
rotation (starting at its smallest vertex by (x, y), counter-clockwise), so
the triangle's areas do not depend on how its row lists the vertices. Every
stretch of a Hilbert curve covers a compact region, and the triangles whose
centroids fall in an aligned square of the curve's grid are one stretch:
that is the property the guarantee below uses.

The exact overlap runs over the curve in chunks of consecutive triangles
(2¹⁸ by default), each with its own raster window: the chunk's index-space
box plus one cell, decoded from the cache. A chunk is compact, so its window
is small, and memory is bounded by the chunk, not the basin.

### 1.4 The cutoff and the ledger

**Parameters.** A relative cutoff `c` (`--land-cover-cutoff`, percent,
default 5) and an absolute floor `m` (`--land-cover-floor`, hectares,
default 1). An entry of area `x` in a triangle of area `A` is **small** when
`x < c·A` **and** `x < m`. The rule for where a class may be placed is fixed
at "present" (question 4). `--land-cover-cutoff 0` writes the exact
fractions, with no ledger.

**State.** The ledger `L`: one signed area in m² per class, starting at
zero. Positive means "this much of class k is owed", negative "this much too
much was written". `Σₖ L[k] = 0` always, to rounding.

**Per triangle `t`, in curve order**, with exact class areas `a[k]` summing
to `A`:

1. `v[k] = a[k] + L[k]` for every class in either list.
2. If `t` is forced to water (section 2.6), the output is `{33: A}`; go to 6.
3. The candidates `S` are the classes with `a[k] > 0`, `v[k] > 0` and `v[k]`
   not small. If `S` is empty, `S` is the one class with `a[k] > 0` and the
   largest `v[k]` (ties to the smaller class number).
4. `o[k] = A · v[k] / Σ_{j∈S} v[j]` for `k` in `S`.
5. If some `o[k]` is small and `S` has more than one class, remove every
   small one but the largest and repeat step 4. This ends: `S` shrinks.
6. `L[k] ← v[k] − o[k]` for every class seen (`o[k] = 0` outside `S`).
7. Write `o[k] / A` as the fractions, in ascending class order; the largest
   (ties to the smaller class) is also the triangle's `land_cover_code`.

What it does in words: a triangle receives what the raster put in it plus
whatever the triangles before it on the curve owe or overdrew; it keeps the
classes that are big enough and present in it, scaled to fill its area
exactly; what it could not take, or took too much of, goes on to the next
triangle. Ola's "95 % corn, 5 % soybean could become 100 % corn" is step 3
with the soybean small; the 5 % then waits in the ledger until the curve
reaches a triangle that has soybean in it, and is written there.

**Implementation.** C++, `include/terrain/land_cover/ledger.hpp`: a
`Ledger` object with `push(chunk of class areas, planar areas, forced mask)
-> chunk of fractions`, called by the driver for the chunks in curve order;
its state is the dense vector `L` (256 doubles) and the list of classes
touched, so a step costs the classes of the triangle and of the ledger, not
256. Serial by nature (each step reads the state the previous one left), and
cheap: the probe, scanning all 256 slots, took 1.0 µs per triangle.
Determinism follows from the order (1.3) and from serial evaluation.

### 1.5 What it guarantees, and what it does not

By construction, and each tested (under "Invariants"):

- **Per triangle**: the fractions sum to one (to the float32 rounding of the
  file); no written entry is small, unless it is the triangle's only class;
  a class appears only where the raster has some of it in that triangle,
  except water in a forced triangle.
- **Over any stretch of the curve**, for every class, the area written minus
  the exact area equals the ledger before the stretch minus the ledger after
  it (sum steps 1 and 6 over the stretch). So the error of a class over a
  stretch is at most twice the largest ledger entry seen, whatever the
  stretch's length: over the whole mesh it is the ledger left at the end;
  over the triangles whose centroids fall in any square of the curve's grid
  it is at most twice the largest entry.
- **Over any region**, the error is at most twice the largest entry times
  the number of stretches the region cuts the curve into. For a catchment
  that number grows with its perimeter, so for irregular regions the bound
  is weak and the measured figure (the 5 km and 25 km squares above, which
  are not aligned to the curve) is the useful one.

Not guaranteed, measured and reported on every run: **the size of the
ledger.** No published bound covers unequal quanta with a threshold (see
the literature). On the Corrente it stayed between 4.4 and 4.7 ha at 20, 10,
5 and 2 m with the 1 ha floor, and the end remainder was at most 0.95 ha. A
run reports its largest entry and the end remainder, per class, in `--stats`
and `--record`, so a run where the ledger grows is seen, not hidden.

**That figure was measured without water forcing** (2.6): the Corrente
meshes have no lake constraints. With lakes forced to 100 % water, every
land cell inside a lake polygon (the shore displacement of the reduction,
within `τ` of the outline, and every filled island) enters the ledger as
owed land and every water cell outside it as owed water. Along a shore the
curve passes between the two sides often and pays them back; but inside a
large lake, where every triangle is forced, nothing can be paid back until
the curve leaves the lake, so the ledger may grow with a lake's filled
islands and shoreline, not with the floor. How much is not known: it is a
measurement in `@perf`'s acceptance (section 11, item 2b), and the cost is
named in question 8.

### 1.6 What the mesh file carries

Fractions are needed by the mesh's user, so they go in the mesh file
(increment 25's rule: the file carries what a user of the mesh needs). Not
every triangle has the same number of classes (1 to 11 on the Corrente),
and a fixed width would pad every triangle to the largest: six to ten times
the data, at 1.1 to 1.7 entries per triangle. So the file stores the sparse
table once, compressed-row style:

**`.vtk`** (legacy, as today):

- cell array `land_cover_code` (`int`): the largest class after the cutoff,
  0 on every line cell. The name and the fill are 16c's, so the ParaView
  route of `rasputin palette` works once a MapBiomas palette exists (not in
  this increment).
- FieldData arrays (dataset-level, any length; ParaView reads them, a
  script uses them): `land_cover_offsets` (`unsigned_int`, one more than the
  number of triangles: triangle `i`'s entries are `offsets[i]` up to
  `offsets[i+1]`, triangles counted in file order after the line cells),
  `land_cover_classes` (`unsigned_char`) and `land_cover_fractions`
  (`float`, a type the writer gains). A mesh with more than 2³² − 1 entries
  is refused (the basin at 1 m has about 0.3 G).
- FieldData strings: `land_cover_codes` (16c's field; here "MapBiomas
  Collection 11 class number; 0 = no data, and every constraint line"),
  `land_cover_source` ("MapBiomas Collection 11, 2025"), and the licence
  fields the data's owner requires, `land_cover_credit`,
  `land_cover_licence_note` and `land_cover_cite`, as 23a-2 wrote them for a
  downloaded DEM.

**`.ply`**: the face element gains `property list uchar uchar
land_cover_classes`, `property list uchar float land_cover_fractions` and
`property int land_cover_code` (PLY lists are sparse natively), and the same
strings as header comments.

**`--stats` and `--record`** (increment 25's record): the cutoff
(`land_cover_cutoff_pct`), the floor (`land_cover_floor_ha`), triangles
forced to water, and per class the exact area and the written area in km²,
their difference, the largest ledger entry in m² and the end remainder in m².
A largest entry above 10 times the floor prints a warning to stderr (25's
rule for self-checks).

## 2. Water bodies as constraints

### 2.1 Which cells, and how small is too small

- **Water** is MapBiomas class 33, "river, lake and ocean". Class 31,
  aquaculture (fish ponds), stays a fraction unless Ola rules otherwise
  (question 6).
- **A water body** is a set of water cells connected through edges or
  corners (8-connected, as ruled). Two bodies are therefore always at least
  one land cell apart, and never share a point.
- **The horizontal tolerance** `τ` (`--water-tolerance`, metres, default
  60 m, twice the cell, as increment 22's `--outline-tolerance` defaults to
  twice the DEM's cell) is how far the reduced outline may stand from the
  traced one, at vertices (22's guarantee).
- **The minimum area** is `4τ²` (`--water-min-area-factor`, default 4):
  1.44 ha at the default, about 17 cells. A body below it is not traced; its
  cells stay class 33 in the fractions of the triangles they fall in, as
  ruled ("a smaller one stays a fraction ... nothing is merged"). On the
  three measured units, the bodies above 1.44 ha hold 89 to 98 % of the
  water and are 6.5 to 9.5 % of the bodies.
- **Islands.** A hole in a water body (land, connected through edges only,
  since land is the complement of 8-connected water) of at least the minimum
  area is kept as a hole; a smaller one is filled, so its land becomes part
  of the lake polygon. Its area then reaches the fractions through the
  forced-water rule (2.6) and the ledger, like every other land cell inside
  the polygon.

### 2.2 Labelling and tracing, on cell edges

C++, `include/terrain/raster_vector/water.hpp`, pure: rows of class values in,
rings out.

- **Labelling by runs** (He, Chao and Suzuki 2008): each row's water cells
  as runs; a run joins every run of the row above that it touches through an
  edge or a corner (8-connectivity: the column ranges overlap after
  widening by one); union-find over runs; area per body accumulated per run
  as cell areas by row. The input arrives a strip of rows at a time (a row of
  256-row tiles decoded from the cache). Labelling needs only the previous
  row's runs; tracing needs the runs of the bodies it traces; so memory is
  the water's runs, not the window: the basin's box is 2.15 G cells, its
  water a few percent of that (2.7 % of unit 761, 0.3 % of 769).
- **Tracing.** For every body at or above the minimum area, its boundary on
  the cell edges ("cracks"): an edge between a water cell and a land cell,
  oriented with the water on its left. Followed corner to corner; at a
  corner with two ways out (the diagonal case, water at two opposite
  corners) the turn that keeps the water connected is taken (the right turn,
  with the water on the left), so the ring passes that corner twice: a
  **pinch**. Corners where the direction does not change are not output.
  Output: closed rings of integer cell corners `(col, row)`, outer rings
  counter-clockwise and holes clockwise in a frame with y up, each with its
  body number and area.
- **Determinism.** Bodies are numbered by their first cell in raster order;
  each ring starts at its lowest corner in (row, col) order. A function of
  the raster alone.
- **The window** is the domain's box grown by `τ` and two cells. A body that
  reaches the window's edge is closed along the edge. That piece of ring
  lies outside the domain's box, so it never meets the mesh, and the area
  inside the domain is exact (2.4).

`water_probe.py` is this, in Python and with labelling by cells, run on the
three units above: 6,383 bodies traced, every ring closed, the pinches
visited twice.

### 2.3 Into metres

Every ring corner is moved once: the raster's affine to longitude and
latitude, then pyproj to the mesh's CRS (the same `crs.reprojector`). A
corner shared by two rings (a hole touching its outer ring at a pinch)
moves identically for both. Then all rings are shifted by one fixed origin,
the window's lower-left corner in metres, as 22 does, so coordinates stay
below about 2 × 10⁶ m and every ring is in the same frame (the check below
compares rings with each other).

### 2.4 The reduction: increment 22's collapse, for many rings and fixed lines

`reduce_rings` in `include/terrain/vector_simplify/area_collapse.hpp`,
generalising `reduce_ring` (which stays, as the one-ring case):

```
struct FixedLine { std::span<const Point2> points; };   // a polyline that may not be crossed
ReduceManyOutcome reduce_rings(std::span<const Point2> vertices,
                               std::span<const std::uint32_t> ring_offsets,  // R + 1
                               std::span<const FixedLine> fixed,
                               double tolerance);
// -> rings (same layout), status per ring, counts (collinear, collapses,
//    rejected for crossing, for a fixed vertex inside, for tolerance, anchors)
```

**Fixed lines** are every line the mesh will also hold and the lakes must not
newly cross: the domain's outline and holes, and the lines of every
`--features` source given to `rasputin water` (today the river lines; the
DEM-derived ones when they exist).

The steps, against increment 22's ("The reduction" in
`22-auto-catchment.md`):

1. **Anchors (new).** Every crossing or touch of a ring with a fixed line is
   inserted into the ring as an **anchor**: a vertex at the crossing point
   (computed in doubles from the two segments; detected with the exact
   predicates) that is never collapsed. It may be the `A` or `D` of a
   collapse, never `B` or `C`. A stretch where a ring runs along a fixed
   line (collinear overlap) gets anchors at both ends, and every vertex in
   between is frozen too. The fixed line itself is not changed; the noder
   nodes the crossing later, and the anchor lies on the fixed segment to
   rounding, well inside its 1 mm snap.
2. **Collinear pass**, as 22, never dropping an anchor or a pinch.
3. **Candidates and deviation**, as 22, per ring, against that ring's traced
   vertices: the same area rule for `E`, the same deviation, the same
   tolerance test at vertices.
4. **One order over all rings**: a heap on (deviation, ring, vertex id),
   least first. Serial and deterministic, as 22.
5. **Checks before a collapse**, with the exact kernel:
   - *No crossing*, now against every current edge of every ring and every
     fixed segment, through one edge grid. Allowed touches: at `A` and `D`
     with their own neighbouring edges (22); at a **pinch** (new): a new edge
     may end at a point another visit of the ring (or another ring of the
     same body) also passes, provided the two visits do not interleave
     there: around the point, the two edges of one visit lie on one side of
     the other visit's two edges (four exact orientation tests). This is the
     "weakly simple" condition for the simplest case, two visits of one
     vertex. A new edge that passes through such a point, rather than ending
     at it, is refused. At an **anchor** (new): the new edge leaving it must
     lie on the same side of the fixed line as the edge it replaces (one
     orientation test), so the ring still crosses there, once.
   - *Swept region* (generalising 22's keep-point test): no vertex of
     another ring, of another visit of a pinch, or of a fixed line lies
     strictly inside the loop `A-B-C-D-E-A` (winding number) or on a new
     edge. Candidates come from the same grid.
6. **Apply**, re-evaluate the four neighbours, repeat until nothing is
   admissible or each ring is at its floor (four vertices, or its anchors
   and pinches).

**Guarantees**, each tested:

- each ring's area equals the traced ring's to rounding, so each body's
  area (outer minus holes) does too;
- the area of each body **inside the domain** is kept too: a collapse's
  loop never contains a fixed vertex and its new edges cross no fixed line,
  so the area it moves stays on one side of the outline;
- every ring is weakly simple (simple, except pinches that stay closed),
  no two rings cross, and no ring crosses a fixed line except at its
  anchors, which are the traced ring's own crossings, one for one;
- every traced vertex is within `τ` of the reduced ring and every reduced
  vertex within `τ` of the traced ring (22's guarantee; the Hausdorff
  distance is measured, not guaranteed, as in 22);
- deterministic, for any thread count (serial).

**Size.** The largest input is not a lake but the river network: in unit
761 one ring of 157,138 traced vertices with 3,159 islands, in 769 one of
182,794 with 620, each cut at its unit's window; on the whole basin very
likely one body from the upper reservoirs to the sea, of the order of a
million traced vertices and several thousand rings. The reduction is
serial (one heap: 22's determinism) and `n log n`: today's `reduce_ring`
took 0.8 and 1.0 s on those two rings at 60 m, so about 10 s for a
million is the expectation, and the edge grid's bucket side (the tolerance,
at least the mean traced edge) keeps each query to a few buckets however
many rings share the grid. Memory is a few hundred bytes per traced vertex,
a few hundred MB at worst. If the serial heap proves too slow, the
published way out is 22's own note: cut the rings into stretches fixed at
their ends and reduce stretches in parallel; not designed here.

**The rejected alternative**: open every pinch by a nanometre before the
reduction and run 22's check unchanged. The probe did that and it works,
but it moves vertices by hand (so the area is no longer the traced area
bit for bit), the size of a safe nudge depends on the coordinates' magnitude
(1 nm is about 70 units in the last place at 10⁵ m), and the output then
holds vertex pairs a few nanometres apart that the noder must merge.

### 2.5 Into the mesh

A new command, like increment 22's `rasputin catchment`:

```
rasputin water --source mapbiomas-c11 [--year 2025] --domain D [--domain-crs ...] \
               --out-crs C [--features RIVERS ...] [--water-tolerance 60] \
               [--water-min-area-factor 4] --out water.geojson
rasputin mesh --dem anadem-v1 --domain D --out-crs C --tolerance T \
              --features water.geojson --features-map mapbiomas-water \
              --land-cover mapbiomas-c11 [--land-cover-year 2025] --out basin.vtk
```

- `water` writes a GeoJSON `FeatureCollection` in the mesh's CRS (with its
  `crs` member, as `catchment` writes), one feature per water body: a
  `MultiPolygon` (a body whose reduced ring still touches itself at a pinch
  is split there into parts that touch at a point, which is valid), with
  properties `code` 33, `year`, the source and its credit, the traced and
  reduced vertex counts, and the area in m². Coordinates at `repr`
  precision so the area survives the round trip (22's rule).
- `--features-map mapbiomas-water` is a new class map (16b's `ClassMap`):
  every feature gets the existing `water` edge property (bit 8 of the
  vocabulary) and, being a coded map, its polygon and code for 16c's labels
  (code system "MapBiomas Collection 11 class number"). Nothing else in 16b
  changes: the rings are linework, clipped to the domain, noded with the
  other features.
- Why a separate command and not a flag of `mesh`: the polygons can be
  looked at and reused (several tolerances of mesh on one set of lakes), the
  step is tested without meshing, and `mesh` keeps reading features from
  files only. The cost is one more command in the recipe.

### 2.6 Triangles inside a lake are water

16c's labelling gives every triangle inside a water polygon the code 33
(components across unconstrained edges, one point-in-polygon test per
component). In the fractions pass such a triangle is **forced**: written as
100 % water, its exact areas still entering the ledger (step 2 of 1.4). The
land cells the reduced outline put inside the lake (within `τ` of the shore)
and the water cells it left outside are equal in area, because the polygon
keeps the traced area (filled islands apart, which are land inside by
definition). Either way the ledger carries what forcing displaces, so class
totals stay exact.

**The cost, not yet measured.** The ledger can only pay owed land back in a
triangle that is not forced. Along a shore the curve alternates between
lake and land and pays it back within a few triangles; but a large lake's
interior is a long run of forced triangles, and the land owed from its
filled islands and its shoreline waits until the curve leaves the lake. The
4.4 to 4.7 ha of 1.5 were measured on meshes without lakes; with the
basin's river network (3,159 islands in unit 761's river body alone; how
many are below the minimum area, and so filled, was not counted) the
largest ledger entry may be much
larger, and the land it carries is written on the next land triangles the
curve meets, which can be some way from the island it came from. `@perf`
measures this (section 11, item 2b); if it is large, the remedy to weigh is
to keep more islands as holes (a smaller minimum area for islands than for
lakes), not to weaken the ledger.

Forcing applies only when the land-cover year equals the water polygons'
year (the `year` property). With another year (a reservoir that did not yet
exist in 1985, say), lake triangles get their exact fractions like any
other, the lake outline stays a constraint, and stderr says so once.

## 3. Getting MapBiomas

### 3.1 The catalogue entry

`rasputin fetch` (23a-2) copies what a mesh will read into the cache, and
the mesh never touches the network. MapBiomas fits it with four additions
to the catalogue's data model (`sources.py`, `RemoteSource`), none to the
download logic:

| field | ANADEM, GLO-30 today | MapBiomas |
|---|---|---|
| `kind` | `one-cog`, `cog-tiles` | `one-cog-per-year` (new): `url_template` with `{year}` |
| `years` (new) | none | 1985 to 2025; the default is the last |
| `values` (new) | `heights` | `classes`: `--dem` refuses it, `--land-cover` and `water --source` refuse a `heights` source |
| empty tiles (new) | refused when read (23a-1) | read as 0, no data |
| `crs`, `nodata` | as now | EPSG:4326, 0 |

```
mapbiomas-c11
  url_template  https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection11/
                lulc/coverage/brazil_coverage/brazil_coverage-col11_{year}.tif
  object id     brazil_coverage-col11_{year}           (one per year in the cache)
  credit        MapBiomas Project - Collection 11 of the Annual Series of Land Use and Land
                Cover Maps of Brazil, accessed on {date} through the link: {url}
  licence_note  question 2 (CC BY or CC BY-SA 4.0)
  cite          Souza et al. (2020), Reconstructing Three Decades of Land Use and Land Cover
                Changes in Brazilian Biomes with Landsat Archive and Earth Engine.
                Remote Sensing 12(17), 2735. https://doi.org/10.3390/rs12172735
```

The credit format is MapBiomas's own ("MapBiomas Project- Collection
[version] of the Annual Series of Land Use and Land Cover Maps of Brazil,
accessed on [year] through the link: [LINK]", terms of use page,
`sites.mapbiomas.org/conheca-o-mapbiomas/termos-de-uso/`, read 2026-10-03;
the factsheet asks for the date, so the date of the fetch, from the cache's
manifest, is filled in). The date `{date}` is the day the object's header
was first fetched; the manifest gains that date per object.

`rasputin fetch mapbiomas-c11 [--year 2025] --domain D --out-crs C` plans
the tiles meeting the domain's box, as for ANADEM. Sizes per year, measured:
177 MB for the basin's box, 5.1 MB for the Velhas piece's. The header is
complete within 8 MiB (23a-2's doubling rule finds it). Decoding needs
`imagecodecs` (LZW with the horizontal predictor), the existing `codecs`
extra; without it, 11's refusal names the extra.

**Reading.** A class window is a `uint8` array with its corner, cell size
and shape (`PixelIsArea`: the tie point is a cell's corner, not a node as
for the DEMs), read through 23a-1's block cache (`decode_window`), with
empty tiles as zeros. A new frozen model, `ClassWindow`, beside the DEM's
`RasterMeta`; the DEM path is untouched.

### 3.2 Which collection

| collection | released | years | classes | notes |
|---|---|---|---|---|
| 11 | 12 Aug 2026 (factsheet); files dated 18-19 Aug 2026 | 1985-2025 | 33 mapped (four new: flooded savanna, salt marsh, herbaceous and shrub formation, wind farm, all marked beta; cotton, a crop the Corrente grows, is beta too) | the latest; on the bucket as `collection11/` |
| 10.1 | 9 Feb 2026 | 1985-2024 | as 10 | fixed river, lake and ocean in the Amazon for 1985-2003, and no-data pixels in the Amazon, Pampa and the coast |
| 10 | before 10.1 (date not checked) | 1985-2024 | | on the bucket as `collection_10/` (`brazil_coverage_2024.tif` answered 200) |

**Recommended: Collection 11, year 2025.** It is the latest, it has the
latest year (Ola's default), and its new classes are small in area (0.84
Mha together over Brazil, against 66.5 Mha of agriculture; wind farms and
the herbaceous and shrub formation may occur in the basin, which was not
checked). An 11.1 may follow as 10.1 did; the
catalogue key carries the collection (`mapbiomas-c11`), so a later key is an
added entry and an old mesh's record still names what it used.

## 4. Pieces and seams

The basin is meshed in pieces when it does not fit the memory budget
(increment 23): the domain is cut along lattice lines, each piece is meshed
on its own with its seams frozen, and the stitched file has the seams
removed. How land cover fits:

- **Water bodies are input geometry**, computed once for the whole domain
  before it is cut (the "global vector step" of 23, which never reads the
  DEM). A lake crossing a cut is two pieces of linework like any feature;
  where a cut meets a lake outline the noder makes a vertex, which 23 keeps
  ("a corner where a seam meets ... a feature is an input-constraint vertex
  and is kept"). Nothing new.
- **Fractions and the ledger are a function of the finished mesh**: one
  Hilbert curve over all its triangles, one ledger. In a cut run they are
  computed on the stitched mesh, after the seam cleanup, never per piece.
  So the result does not depend on how many pieces there were or where the
  cuts ran, and there is no ledger to hand from piece to piece. Piece files
  (`--no-stitch`) carry no land cover (question 9). A piece-by-piece ledger
  was considered and rejected: its triangles near seams change in the
  cleanup, so it would have to be redone there, and its result would depend
  on the cut.
- **The cost** is one pass over the stitched mesh, which must be read back
  (the stitcher writes it streaming): about 11 GB at 1 m for the whole
  basin (section 5), under the 16 GB budget, and serial in its ledger part.
  This waits for the stitcher (23d, 23g); until then every mesh is one
  piece, and the land-cover pass runs inside `mesh` on the mesh in memory.
- **Forced water** uses 16c's labels on the same mesh, so it too is
  computed once, on the stitched mesh.

## 5. Cost, at 1 to 50 m

Triangle counts for the whole basin: at 2, 5, 10 and 20 m the sums of
`@perf`'s nine ANADEM level-3 meshes (633,624 km²,
`docs/benchmarks/2026-10-02/basin-level3/README.md`); at 1 and 50 m, which
were not meshed on ANADEM, the GLO-30 estimate from 200 random boxes
(635,194.5 km², `docs/benchmarks/2026-10-01/basin-piece/README.md`; a surface
model, so high at fine tolerances). A 30 m cell is about 864 m².

| tolerance | basin triangles | mean triangle | in 30 m cells | classes per triangle written | land cover in the file, `.vtk` / `.ply` | memory of the pass, one piece | ledger, serial |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 237 M (GLO-30) | 0.27 ha | 3.1 | ~1.1 | 3.2 / 2.7 GB | ~6 GB | ~50 s |
| 2 m | 120.2 M | 0.53 ha | 6.1 | 1.14 | 1.6 / 1.4 GB | ~3 GB | ~25 s |
| 5 m | 34.6 M | 1.8 ha | 21 | 1.36 | 0.51 / 0.44 GB | 0.9 GB | ~7 s |
| 10 m | 13.1 M | 4.8 ha | 56 | 1.50 | 0.20 / 0.18 GB | 0.3 GB | ~3 s |
| 20 m | 5.0 M | 12.7 ha | 147 | 1.65 | 81 / 71 MB | 0.12 GB | ~1 s |
| 50 m | 1.0 M (GLO-30) | 64 ha | 735 | ~1.8 | 17 / 15 MB | 25 MB | < 1 s |

How each column is made:

- **Mean triangle**: the area over the triangles. At 1 m it is 3.1 cells
  over the basin; the research note's 1.2 cells is the Velhas piece's, which
  is steeper than the basin (2.4 times the basin's density at 1 m, per the
  basin-piece README). Either way, at 1 m a triangle holds a few cells, and
  the cutoff there mostly removes the slivers of cells cut by its edges.
- **Classes per triangle written**: the Corrente's, with the 5 % cutoff and
  the 1 ha floor (section "What was measured"); 1 and 50 m extrapolated from
  the trend (2 m: 1.14; 20 m: 1.65). One sub-basin only: steeper and more
  mixed land will have more.
- **File**: binary. `.vtk`: 4 bytes of offset, 5 per entry (class and
  float fraction) and 4 for `land_cover_code` per triangle; `.ply`: 2 bytes
  of list counts, 5 per entry and 4 for the code. Today's binary mesh is
  about 36 bytes per triangle (unit 761 at 2 m: 977 MB for 27.1 M
  triangles), so land cover adds about a third. Text files are larger
  (not measured).
- **Memory of the pass**, over what meshing already holds: the curve order
  (a 64-bit key and a 32-bit index per triangle, 12 bytes) and the output
  (as in the file), about 25 bytes per triangle, against the 310 bytes per
  triangle meshing takes (the basin-piece fit), so under 10 %. The exact
  areas live only per chunk. **For a cut run**, where the pass reads the
  stitched mesh back, add the mesh itself (vertex x and y, 16 bytes per
  vertex, about 8 per triangle; triangles 12): about 45 bytes per triangle,
  **about 11 GB at 1 m for the whole basin**, under the 16 GB budget.
- **Ledger**: serial, at the 0.2 µs per triangle a sparse implementation
  should reach (the probe's dense one did 1.0 µs). An estimate for `@perf`
  to replace.
- **The exact overlap** is not in the table: it is parallel, and the probe's
  cost (0.15 µs per covered cell plus 2.2 µs per triangle, one thread) gives,
  for the basin's roughly 735 M covered cells and 237 M triangles at 1 m,
  about 10 minutes on one thread in probe code, so a minute or two on 8
  cores before the full-cell spans make it cheaper. Moving 119 M vertices
  through pyproj is the other cost not measured here.
- **Fetch**: per year, 72 MB for the tiles meeting the basin's outline and
  177 MB for its box, against ANADEM's 1.71 GiB for the blocks meeting the
  outline.

**Water bodies add constraint vertices.** Scaled from the three measured
units, the basin has about 11,600 bodies of at least 1.44 ha and about
130,000 reduced vertices at 60 m (the river network's ring in 761 alone
12,655, in 769 12,925). Each vertex of a
constraint forces a few triangles where the mesh would otherwise be coarse,
so the added triangles are a few times 130,000: perhaps 5 to 15 % at 20 m
(5.0 M triangles), well under 1 % at 1 m (237 M). 16b's 3.6 times on CORINE
came from every class border; water alone is a small part of that. An
estimate; the acceptance measures it.

## 6. The blueprint

```
rasputin fetch mapbiomas-c11 [--year Y] --domain D --out-crs C           (26a)
   catalogue entry -> 23a-2's planner and downloader, unchanged     -> cache

rasputin water --source mapbiomas-c11 [--year Y] --domain D --out-crs C [--features F] ...  (26d)
   water.py (Python, pure apart from the cache read)
     ClassWindow strips from the cache (io/repository, 23a-1)          ── uint8 rows
     tracer = _core.WaterTracer(water classes, cell_area_by_row, min_area)
     tracer.push(strip) for each strip; tracer.rings() ──────────────── C++, raster_vector/water.hpp
        -> rings of integer corners, body ids, areas                       (labels by runs, cracks)
     corners -> lon/lat -> mesh CRS (crs.reprojector), minus one origin
     fixed lines = domain outline and holes + lines of F, same CRS
     _core.reduce_rings(rings, fixed, tolerance) ─────────────────────── C++, vector_simplify/area_collapse.hpp
     MultiPolygons (split at closed pinches) -> water.geojson, with code, year, credit

rasputin mesh ... --features water.geojson --features-map mapbiomas-water
                  --land-cover mapbiomas-c11 [--land-cover-year Y] [--land-cover-cutoff 5]
                  [--land-cover-floor 1]                                   (26b, 26c)
   the mesh as today (water rings are 16b linework; 16c labels lake triangles 33)
   land_cover_fractions.py (Python driver)
     vertices -> lon/lat -> index space of the class grid (crs.reprojector)
     _core.hilbert_order(centroids) ─────────────────────────────────── C++, land_cover/ledger.hpp
     for each chunk of 2^18 triangles along the curve:
        ClassWindow of the chunk's box, from the cache
        _core.class_areas(corners, window, row_area) ───────────────── C++, land_cover/overlap.hpp (parallel)
        ledger.push(areas, planar areas, forced mask) ───────────────── C++, serial
     -> offsets, classes, fractions, land_cover_code, the ledger report
   writers: .vtk FieldData + cell array, .ply face lists; run record entries
================================ _core boundary ================================
C++ sees: integer class rows and a cell area per row; triangle corners in the
class grid's index space; planar areas; rings and fixed lines in metres in a
local frame. No CRS, no path, no year, no class names.
```

**Boundaries.** `land_cover_fractions.py` and `water.py` import NumPy, shapely (water
output only) and first-party I/O through `io/repository.py`, never a path
below the CLI; the C++ headers are pure and know no Python. The class
numbers' meaning (33 is water) lives in the catalogue entry and the class
map, in Python; C++'s forced mask and water class are numbers it is handed.
Every step is a function of its inputs, testable alone: the overlap on a
three-triangle mesh and a 4 × 4 class array; the ledger on hand-made class
areas; the tracing on a 6 × 6 mask; the reduction on hand-made rings.
The module is `land_cover_fractions.py`, not `land_cover.py`, so that it
cannot be confused with 16c's `landcover.py` (labels from polygons), which
it calls for the lake labels and does not replace.
**Async**: `land_cover_fractions.assign` and `water.extract` are blocking with the C++
calls releasing the GIL; an API worker runs them in `asyncio.to_thread`, as
it runs `catchment.delineate`.

## 7. Invariants

Fractions (26b, 26c):

- **F1 (one per triangle).** Every triangle's fractions sum to 1 within
  1e-6 (float32 in the file), with classes ascending and none repeated.
- **F2 (exact areas).** Before the cutoff: for each class, the areas summed
  over the triangles equal the class's area inside the mesh's footprint to
  1e-9 relative; checked against an independent clip (shapely's
  intersection of each triangle with each cell box, on small fixtures), and
  per cell: the pieces of a cell over all triangles sum to the cell's area
  covered, to 1e-12 of a cell.
- **F3 (the threshold).** No written entry is small (below `c` of its
  triangle and below the floor) unless it is the triangle's only class; no
  class is written where the raster has none of it in that triangle, except
  water in a forced triangle; a forced triangle is exactly `{33: 1}`.
- **F4 (the stretch identity).** For every prefix of the curve and every
  class, written minus exact equals minus the ledger after the prefix, to
  1e-9 of the prefix's area. The test pushes chunks of one triangle and
  reads the ledger after each.
- **F5 (order).** Permuting the triangle rows, or the vertex numbering,
  gives the same fractions per triangle; so does any thread count.
- **F6 (no cutoff).** `--land-cover-cutoff 0` writes the exact fractions and
  leaves the ledger at zero.

Water (26d):

- **W1 (bodies).** Every 8-connected body of at least the minimum area is
  traced, none smaller; two cells touching at a corner are one body; a
  hole below the minimum is filled.
- **W2 (area).** Each reduced ring's area equals its traced ring's to 1e-9
  relative, measured with one cross product per edge summed exactly
  (`math.fsum`; the two-dot-product shoelace is not precise enough, see
  "What was measured"); each body's area inside the domain likewise.
- **W3 (topology).** Each ring weakly simple (shapely `is_valid` after
  splitting at pinches, and a brute-force pairwise test with exact
  orientation on small cases); no two rings cross; the number of crossings
  of each ring with each fixed line is the traced ring's.
- **W4 (sides).** Every vertex of every fixed line and of every other ring
  has the same winding number with respect to each ring before and after.
- **W5 (tolerance at vertices)** as 22.
- **W6 (determinism)** as 22: the same bits for the same input.

## 8. Degeneracy policy

| case | outcome |
|---|---|
| a triangle vertex exactly on a grid line or corner | zero-area pieces skipped; the canonical crossing makes both neighbours agree |
| a triangle wholly inside one cell | one entry, that class, fraction 1 |
| a triangle partly outside the class window | refusal `OutsideWindow` (the window is built to contain every vertex, so this is a defect) |
| class 0 (no data) inside the mesh | a class like any other, reported in `--stats` |
| a triangle of zero planar area | none in a valid mesh; refused if met |
| every class of a triangle small (more than 1/c classes, all under the floor) | the largest kept alone (step 5 keeps one) |
| the ledger's largest entry above 10 times the floor | written, with a warning on stderr |
| two water cells touching at a corner | one body; the ring visits the corner twice (a pinch) |
| a hole touching its outer ring at a corner | both rings visit it; the pinch rule covers rings of one body |
| a body reaching the window's edge | closed along the edge, outside the domain; area inside the domain exact |
| a ring passing through a fixed line's vertex | an anchor there, as for a crossing |
| a ring running along a fixed line | anchors at both ends, the vertices between frozen |
| a fixed line touching a ring without crossing | an anchor, so the touch is kept and not turned into a crossing |
| two fixed lines crossing on a ring | one anchor, at the shared point |
| a reduced ring that would fall below four vertices | stops at its floor (22) |
| the land-cover year differs from the water's | no forcing; the outlines stay; one stderr line |
| no `--land-cover` | no arrays, no fields, as today |
| `--land-cover` with a CORINE map in the same mesh | refused (16e's rule: one code system per mesh) |

## 9. The PR split

Counted in `CLAUDE.md` §2's unit. Each estimate with the worst case at
+39 % (increment 10's overrun), as 23 does; every PR stays under 700 there.

| PR | what | estimate | +39 % |
|---|---|---:|---:|
| **26a** | **MapBiomas in the catalogue, and class windows** | | |
| | `sources.py`: `one-cog-per-year`, `years`, `values`, empty tiles as no data, the entry | 35 | |
| | `fetch/plan.py`, `fetch/run.py`: the year, the object id per year, the fetch date per object | 30 | |
| | `io/models.py`, `io/cog.py`: `ClassWindow`, `uint8` blocks, empty tiles as zeros | 45 | |
| | `io/repository.py`: a class window over a box, in strips | 40 | |
| | `cli.py`: `fetch --year`, refusals by `values` | 20 | |
| | **26a total** | **170** | **236** |
| **26b** | **Exact fractions per triangle** (no cutoff) | | |
| | `land_cover/overlap.hpp`: strips, full-cell spans, canonical crossings, chunks in parallel | 150 | |
| | `bindings/core.cpp`, `_core.pyi`: `class_areas`, the `uint8` view | 55 | |
| | `land_cover_fractions.py`: index-space corners, row areas, chunks, dominant class | 100 | |
| | `io/vtk_legacy.py`, `io/ply.py`: FieldData arrays, `float`, face lists | 60 | |
| | `cli.py`: `--land-cover`, `--land-cover-year`, record entries | 55 | |
| | **26b total** | **420** | **584** |
| **26c** | **The cutoff and the ledger** | | |
| | `land_cover/ledger.hpp`: Hilbert keys and order, the quantiser, `Ledger::push`, the report | 170 | |
| | `bindings/core.cpp`, `_core.pyi` | 45 | |
| | `land_cover_fractions.py`: curve order, forced water from 16c's labels, the year check | 60 | |
| | `cli.py`, `run_record.py`: cutoff and floor options, the ledger report, the warning | 45 | |
| | **26c total** | **320** | **445** |
| **26d-1** | **Water bodies, traced** | | |
| | `raster_vector/water.hpp`: labels by runs, areas, cracks, pinches, holes filled | 230 | |
| | `bindings/core.cpp`, `_core.pyi` | 45 | |
| | `water.py`: strips from the cache, corners to metres, polygons | 80 | |
| | `cli.py`: `rasputin water`, writing the traced polygons (no reduction yet) | 70 | |
| | **26d-1 total** | **425** | **591** |
| **26d-2** | **Water bodies, reduced** | | |
| | `area_collapse.hpp`: many rings, one heap and grid, the pinch rule, anchors, fixed lines, swept region | 230 | |
| | `bindings/core.cpp`, `_core.pyi`: `reduce_rings` | 40 | |
| | `water.py`: fixed lines from the domain and `--features`, split at pinches | 50 | |
| | `feature_input.py`: the `mapbiomas-water` class map | 15 | |
| | `cli.py`: `--water-tolerance`, `--water-min-area-factor` | 20 | |
| | **26d-2 total** | **355** | **493** |
| **26e** | **Cut runs, and land cover for another year** (after 23g) | | |
| | read back our own binary `.vtk` and `.ply` (vertices, triangles, labels) | 120 | |
| | the pass on the stitched mesh; `rasputin land-cover MESH --year Y` rewriting it | 90 | |
| | **26e total** | **210** | **292** |

Order: 26a, then 26b and 26c (fractions on today's one-piece meshes: the
level-3 units and the Velhas piece already exist and can be given land cover
without remeshing once 26e's reader exists, or remeshed), then 26d-1 and
26d-2, then 26e when the stitcher exists. 26b and 26d-1 can be built in
parallel once 26a is in. Modules to watch: `area_collapse.hpp` grows from
341 lines to about 570, past which the ring bookkeeping should move to its
own header; `water.hpp` at 230.

**Acceptance class** (`docs/increments/README.md`, "Acceptance"): no PR
touches `include/terrain/refinement/` or `include/terrain/mesh/` or what
drives them, so the 1 m benchmark and thread sweep are not required by the
rule. The water polygons change what is meshed, so 26d-2's acceptance
measures the mesh with and without them anyway (below).

**Mutation testing**: the invariant-critical suites are the ledger's (F3,
F4: a wrong sign or a missed class in the ledger update conserves nothing
and still produces plausible fractions) and the reduction's pinch and
anchor rules (W2-W4). Those two suites get a mutation round; the rest not.

## 10. Tests `@tester` can write red

Lean (Ola's rule): fixtures by hand, no throwaway implementations, the
mutation round only where named above.

**26a** (Python): the catalogue entry's fields and URL for 2025 and 1985; a
year outside the range refused; `--dem mapbiomas-c11` refused and a heights
source refused as land cover; planning the tiles for a box against a
synthetic tiled `uint8` GeoTIFF on the local test server (23a-2's fixture
pattern), with some tiles empty; empty tiles read as zeros; `PixelIsArea`
corner arithmetic (a cell's corner, not its centre); the credit with the
fetch date; LZW with predictor through `imagecodecs`, skipped without it.

**26b** (C++ Catch2 and Python): the overlap on hand cases (a triangle in one
cell; a unit right triangle on a 2 × 2 grid, areas 0.5 per cell; a triangle
with a vertex on a grid corner; a long sliver along a row); F2 against
shapely on random triangulations of a 20 × 20 class grid (the oracle uses
GEOS's clipping, not ours: the producer's predicate, not its records);
per-cell sums; parallel chunks equal one chunk bit for bit; the `.vtk`
FieldData arrays and the `.ply` lists read back (VTK reader, as 16c's test);
`land_cover_code` 0 on line cells.

**26c**: F3, F4 and F6 on hand-made class areas (a triangle of 95 % corn and
5 % soybean followed by one with soybean present: the second receives the
5 %; a class below the cutoff everywhere and present everywhere comes out
whole; a large triangle with a 3 km² minority above the floor keeps it); the
Hilbert order on a known 4 × 4 grid of centroids; F5 by permutation; forced
triangles; the year mismatch line.

**26d-1**: masks by hand: a single cell (four corners); two cells touching at
a corner (one body, one ring, the corner twice); a ring with a hole; a hole
below the minimum filled; a body cut by the window edge; runs across a
strip boundary; areas by row; W1.

**26d-2**: the research note's diagonal strip, pinches closed, reduces (22's
check refuses it today: red by construction); a lake crossed by a river
line keeps both crossings and gains none (W3), the anchors fixed; a lake on
the domain outline keeps its inside area (W2); a small island of another
ring inside a collapse's loop refuses it (W4); a 2,000-vertex traced ring
is a fixture. **Size**: the whole basin's river network is very likely one
body (see "What was measured"), so one ring of the order of a million
traced vertices with thousands of holes. A synthetic channel (a meandering
band 3 to 20 cells wide with a few thousand islands, 10⁶ traced vertices,
generated by the test, not stored) is reduced and traced under a time bound
set from the probe's figures (157,138 and 182,794 vertices reduced in about
1 s each by today's `reduce_ring`), marked slow and run in CI's Release leg
only. It pins that the one heap and the edge grid stay `n log n` with
thousands of rings in the grid.

## 11. `@perf`'s acceptance: a basin piece

On AC power, recorded, against the existing meshes where possible.

1. **Fractions, without remeshing** (26b, 26c): the Rio Corrente (unit 764)
   at 20, 10, 5 and 2 m, the Velhas piece at 10 and 5 m, and unit 769 at
   10 m. Report per mesh: class areas exact against MapBiomas's own count
   (the areas of the cells whose centres fall inside the outline, an
   independent and coarser reference, which may differ by up to the area
   of the cells the outline crosses); written against exact per
   class, with crops listed; the misplaced share on 1, 5 and 25 km squares;
   the largest ledger and the end remainder; classes per triangle; time of
   the overlap (by threads 1, 2, 4, 8) and of the ledger; peak memory.
   **Pass:** every class within 0.1 % or 1 ha of exact; the largest ledger
   within 10 times the floor; the 5 km 95th percentile within 0.5 %;
   numbers close to this design's probe on the same mesh and year (the
   probe kept ellipsoidal, not planar, areas in its ledger), or the
   difference explained.
2. **Water** (26d): `rasputin water` on units 769 and 761: bodies, traced and
   reduced vertices at 30, 60 and 120 m against "What was measured", the
   area check per body, validity, time and memory; the same on the whole
   basin's outline for the counts (the fetch is 177 MB), including the
   size of the largest body and the time to trace and reduce it.
   2b. **Fractions with water forced** (26c with 26d): unit 761 (or 769) at
   10 m meshed with its water polygons, fractions with forcing on and off:
   the largest ledger entry, the end remainder, and the 5 km and 25 km
   misplaced shares, against the no-lake figures of 1.5. **Pass:** reported;
   a largest entry above 10 times the floor goes to Ola with the island
   remedy of 2.6 before 26e is built.
3. **The mesh with water** (26d-2): unit 769 at 20 and 10 m, and the Velhas
   piece at 10 m, meshed with and without `water.geojson`: triangles, refine
   time, peak memory, worst angle. This is the number section 5 only
   estimates.
4. Evidence under `docs/benchmarks/<date>/26-land-cover/`, as `bench.py`'s
   runs are kept.

## 12. Questions for Ola

Each with the recommendation, which is the default if Ola does not rule.

1. **Which MapBiomas collection?** Collection 11 (August 2026, 1985-2025),
   10.1 (February 2026, 1985-2024) or 10. *Recommended: 11, year 2025 by
   default*, as the latest with the latest year; the key names the
   collection, so a later one is an added entry.
2. **The licence of what rasputin makes from MapBiomas.** The Collection 11
   factsheet
   (https://brasil.mapbiomas.org/wp-content/uploads/sites/3/2026/08/Factsheet-Colecao-11-12082026-1.pdf,
   public, read 2026-10-03) says "licença Creative Commons CC-BY";
   MapBiomas's terms of use page
   (https://sites.mapbiomas.org/conheca-o-mapbiomas/termos-de-uso/) says
   "CC-BY-SA" (Attribution-ShareAlike 4.0). Under ShareAlike, three outputs
   are adaptations and must be shared under the same licence: **mesh files
   carrying fractions**; **`water.geojson`**, the traced and reduced water
   polygons; and **meshes constrained by those polygons**, even without
   fractions, since their lake edges are MapBiomas's outlines. rasputin's
   code (MIT) is not affected. *Recommended: write CC BY-SA 4.0 in the
   licence note of all three, the stricter reading, until MapBiomas says
   otherwise*; asking them is Ola's call, as licensing is.
3. **The cutoff.** *Recommended: an entry is dropped only when it is under
   5 % of its triangle and under 1 ha*, both options. Measured: a relative
   cutoff alone makes the ledger carry hundreds of hectares across the map.
   The floor could be 0.5 ha, MapBiomas's own minimum mapping unit; 1 ha
   measured better than 10 ha and is the proposal.
4. **Where may the ledger put a class?** Three rules were measured: only
   into triangles already holding it above the cutoff (Ola's first
   proposal), into triangles holding any of it ("present"), or anywhere
   (halftoning's rule). *Recommended: "present"*: soybean is only ever
   written where the raster has some soybean. With the 1 ha floor, Ola's
   first proposal left a larger ledger (11.5 ha against 4.4 ha at 20 m, 6.7
   against 4.7 ha at 2 m) and more misplaced area on 5 km squares (0.30 %
   against 0.15 % at 20 m, 0.19 against 0.14 % at 2 m), and it lets a crop
   thinly scattered below the cutoff everywhere vanish; "anywhere" is a
   little better locally but writes 8 to 18 % more entries and puts classes
   where none was mapped.
5. **Which regions must be right?** *Recommended: the guarantee as designed
   (any stretch of the Hilbert curve, so any square of its grid, within
   twice the largest ledger entry), and the acceptance measured on 1, 5 and
   25 km squares*, until DEM-derived sub-catchments exist; then on those.
6. **Is aquaculture (class 31) water for the constraints?** *Recommended:
   no*: only class 33; fish ponds stay fractions.
7. **The water tolerance and minimum area.** *Recommended: 60 m (twice the
   cell, as the catchment outline) and four times its square, 1.44 ha*. On
   the measured units that keeps 6.5 to 9.5 % of the bodies, holding 89 to
   98 % of the water; the basin would have about 11,600 bodies and 130,000
   constraint vertices, a large share of them on the river network, which
   is very likely one body from the reservoirs to the sea.
8. **Are triangles inside a lake 100 % water?** *Recommended: yes, when the
   land-cover year is the water's year*; the land inside the outline goes to
   the shore through the ledger. The cost: inside a large lake nothing can
   be paid back, so filled islands and shoreline land may make the ledger
   much larger than the 4.4 to 4.7 ha measured without lakes, and move that
   land some way from where it was mapped. Not measured yet; `@perf`
   measures it before the cut-run work, and keeping more islands as holes is
   the remedy to weigh if it is large.
9. **Piece files in a cut run.** *Recommended: no land cover in piece files,
   only in the stitched file*, so the result does not depend on the cut.
10. **`rasputin water` as its own command**, writing a GeoJSON that `mesh`
    reads as features. *Recommended: yes*, like `catchment`; the polygons can
    be looked at and reused.
11. **When.** The basin order puts "the basin's own inputs" after the basin
    run. 26a-26c touch no refinement code and serve meshes that exist now
    (the nine level-3 units). *Recommended: 26a, 26b and 26c after 15f-3 and
    25, before the remaining basin PRs (23b on); 26d after them; 26e with
    23g.*

## Not in scope

- A natural-colour palette for MapBiomas classes (`rasputin palette
  mapbiomas`); `land_cover_code` is ready for one.
- Rivers from the DEM, and sub-catchments (the other half of "the basin's
  own inputs").
- Class polygons for classes other than water (the research note's full
  coverage simplification), which the hybrid ruling does not need.
- A class map coarsening the legend (Ola: "a class map may coarsen it");
  the fractions keep the full legend, and coarsening them is a sum a reader
  can do.
- Several years in one mesh file.
- A cap on classes per triangle (measured harmful).


## Review

**Design review, round 1, 2026-10-03.** Range `origin/master...` 77f35b4, 9df3822. Verdict: CHANGES REQUESTED. LOC: 0 production lines (design and throwaway probes). The rulings of 2026-10-01 are kept; the fractions-and-ledger figures, the fetch sizes, the CRS figures, the sources and the licence conflict (terms of use CC BY-SA 4.0, factsheet CC BY) all re-checked and hold. Blocking, doc fixes only: (1) unit 769's water row omits its largest body (Três Marias joins the channel, and the probe assigns a body to a unit by its first cell, which lies outside the outline), so the basin-wide lake counts are low; (2) the "Sobradinho" ring is the main São Francisco river across unit 761 (670 × 442 km, 3,159 holes), and "346 km" is half the extent, not a length; on the whole basin the channel likely joins the reservoirs into one body; (3) the ledger measurement excludes lakes forced to 100 % water, which @perf's acceptance must cover; (4) the licence question must also name the water polygons and meshes constrained by them. Not pushed; no CI.
