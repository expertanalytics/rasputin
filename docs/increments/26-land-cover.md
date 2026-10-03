# Increment 26: land cover for the São Francisco basin — water bodies as constraints, class fractions per triangle

Status: **designed by `@architect`, 2026-10-03, at Ola's request; not
ruled.** Written before `@tester` per `docs/increments/README.md` step 1.
Nothing is built. The questions for Ola are at the end, each with a
recommended default. Five pull requests, listed under "The PR split".

**Why "26".** 25 is plain output (`docs/increments/25-plain-output.md`, on
its own branch); 26 is the first free number on every branch
(`git log --all --name-only --format='' -- 'docs/increments/2[5-9]*'` lists
only `25-plain-output.md`). In the basin order of work this is the
"basin's own inputs" item (`docs/increments/23-basin-scale.md`, "Order of
work", item 12), the land-cover half of it; sub-catchments and rivers from
the DEM are the other half and get their own design.

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
   out-of-balance ledger"). Class totals over the basin then stay exact to a
   few hundred square metres, and so does any compact region, within a
   bound the ledger itself reports.

Measured on real basin meshes for this design (the probes are in
`docs/increments/26-probes/`, results under "What was measured"):

- dropping small fractions the way land-surface models do (renormalising
  inside each triangle) loses 5 to 45 % of the open water and up to 13 % of
  the coffee in the Rio Corrente sub-basin; the ledger keeps every class's
  total within 0.05 %;
- the cutoff must be **relative and absolute together**: an entry is dropped
  only if it is under 5 % of its triangle **and** under 1 ha. A purely
  relative cutoff drops whole fields from the large flat triangles (the
  largest is 67 km² at 20 m), the ledger then carries square kilometres
  across the map, and local errors are as bad as renormalising;
- MapBiomas Collection 11 (August 2026, years 1985 to 2025) is on a public
  bucket as one tiled GeoTIFF per year; the whole basin's box is 177 MB per
  year, a tenth of ANADEM's.

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
   remainder, but writes 18 % more entries and can put soybean on a triangle
   that has none; "present" never does. The full tables (cutoffs 1, 5 and
   10 %; floors none, 1 ha and 10 ha; three rules) are printed by
   `summarise.py`.

Time, single-threaded probe code: the exact areas took 6.4, 7.0, 8.2 and
11.3 s at 20, 10, 5 and 2 m (about 0.15 µs per cell covered plus 2.2 µs per
triangle); the ledger took 2.4 s for the 2.33 M triangles at 2 m (1.0 µs per
triangle, scanning all 256 class slots each time).

### Water bodies in three sub-basins

`water_probe.py`: MapBiomas class 33 (river, lake and ocean) over the box of
each BHO level-3 unit, cut into 8-connected bodies (cells touching at a
corner are one body), each counted inside the unit when its first cell in
raster order is inside the outline (a crude test: a body crossing the
outline counts in one unit only, which is why Três Marias, whose northern
end is outside 769, is not among 769's largest). Every body of at least
1.44 ha (four times (60 m)², the minimum area of question 7) was traced on
its cell edges, moved to the basin's Transverse Mercator (the CRS of the
level-3 meshes), its pinch points opened by 1 nm, and reduced by increment
22's `reduce_ring` as it stands (outer rings only; islands not reduced).

| unit | area | water | bodies | ≥ 1.44 ha | their share of the water | ≥ 10 ha | ≥ 10 km² | traced vertices | reduced, 30 m | 60 m | 120 m |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 761, lower basin (Sobradinho, Itaparica) | 209,316 km² | 5,582 km² | 36,561 | 3,454 | 98.4 % | 517 | 5 | 381,202 | 73,254 | 38,034 | 25,231 |
| 769, upper basin | 106,394 km² | 353 km² | 30,174 | 2,834 | 75.9 % | 344 | 3 | 152,862 | 30,440 | 19,406 | 16,144 |
| 764, Rio Corrente | 34,243 km² | 8.6 km² | 1,471 | 95 | 63.2 % | 15 | 0 | 3,156 | 740 | 496 | 417 |

"Traced vertices" counts corners where the outline turns; "reduced" is after
`reduce_ring` at the horizontal tolerance shown. The largest ring,
Sobradinho in 761, has 157,138 traced vertices over 346 km and 4,799 km²;
it reduces to 12,655 vertices at 60 m. Three things the probe settles:

- **The minimum area is what keeps the count down.** 761 has 36,561 bodies,
  but 3,454 of them hold 98 % of its water; in 769, the 2,834 bodies of at
  least 1.44 ha hold 76 %. Scaled by area from these three units (55 % of
  the basin, so a rough figure), the basin has about 11,600 water bodies of
  at least 1.44 ha and about 105,000 constraint vertices at 60 m.
- **Area is kept.** The worst relative change over all 6,383 rings and three
  tolerances is 3.3e-11 (a 5.3 ha pond, 1.8e-6 m²). Measuring it needed
  care: the textbook shoelace as two dot products (Σ x·y' − Σ x'·y) cancels
  catastrophically on Sobradinho's ring and reported a change of 7 m²
  (1.5e-9); one cross product per edge summed with `math.fsum` reports
  none. The tests must measure area the second way.
- **Pinches are common.** 2.3 % of the traced vertices in 769 and 1.3 % in
  761 are pinch points (one corner visited twice: 3,545 and 5,126 of them).
  With a pinch closed, increment 22's crossing check refuses any collapse
  whose new edge touches the other visit of the pinch; on the research
  note's diagonal strip that blocked every collapse. Opened by 1 nm, these
  rings reduce. The design below makes the check accept a closed pinch
  instead, so that no vertex is moved by hand and the area stays exact.

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

- `triangulate_dem.h:846`, `extract_lakes`, bound in `bindings.cpp:371`:
  the triangles whose slope is under 0.01 are "lakes", returned as a
  separate face list. A display heuristic (flat is not water: a flat field
  is flat, a river in a gorge is not). **Not carried.**
- `material_specification.py`: a three.js material for those faces. **Not
  carried.**

`grep -n -i "fraction\|mapbiomas\|dither"` over the tag finds no line. Nothing
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
Collection 11 algorithm description (`ATBD-General-Collection-11-versao-1.pdf`
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
planar to ellipsoidal area (the projection's scale squared) changes by at
most 7e-4 across 20 km in either CRS (`crs_probe.py`), so the choice changes
no fraction measurably; it decides only which total is kept exact, and the
one the user of the mesh sees is the planar one.

**Class 0** (no data in MapBiomas) is a class like any other: it gets a
fraction where it occurs. Inside the basin it occurs only outside Brazil's
outline, which the basin does not reach.

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
Triangles are visited by (index, centroid x, centroid y). The centroids of
two triangles of a valid mesh are never equal (each lies inside its own
triangle, and the triangles do not overlap), so the order is total and a
function of the mesh alone: not of the thread count, the order of rows in the
file, or how the mesh was cut into pieces. Every stretch of a Hilbert curve
covers a compact region, and an aligned square of the curve's grid is one
stretch: that is the property the guarantee below uses.

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
  over any square of the curve's grid it is at most twice the largest entry.
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

