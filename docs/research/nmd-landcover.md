# Proposal: Sweden's national land cover (NMD) in `rasputin mesh`

Status: **proposal, parked**. Written by `@architect`, 2026-10-06, at Ola's
request: "yes, write up the NMD proposal, but we save it for a rainy day (or a
backlock for you to work on while I sleep)". Not a design: it says what NMD
would give, what it costs, and what Ola would rule before a design. Nothing is
built. Its ROADMAP row is unnumbered and marked parked.

## Why

Ola meshed two Swedish catchments on 2026-10-06, Lagan (6,441 km²) and
Ljungan above Flåsjö, with CORINE 2018 as `--features`, and found the land
cover in Lagan "really coarse". It is coarse by design. Lantmäteriet's
product description of the Swedish CORINE says CLC 2018 is the revised 2012
map merged with the changes found by comparing Sentinel-2 images from 2017
with the 2012 images, and that "The surfaces are generalized to a minimum size
of 25 hectares" (source 1, section 1.1 and "Lineage"). In Ola's Lagan extract
the median polygon edge is 75 m (measured 2026-10-06 over the outer rings of
all 7,319 features of `lagan_clc2018_3035.geojson`).

NMD (Nationella marktäckedata, Naturvårdsverket) maps the same country at
10 m. This note is about using it instead of CORINE, or beside it.

## What NMD is

Checked against Naturvårdsverket's own documents and download server on
2026-10-06 (sources 2 to 5):

- **NMD 2018, version 1.1**: a raster of 10 m cells, "a minimum mapping unit
  of 0,01 hectare" (one cell), 25 classes in three levels (16 forest classes,
  open wetland, arable land, two kinds of other open land, three artificial,
  inland water, marine water), code 0 outside Sweden. SWEREF 99 TM
  (EPSG:3006), 8-bit GeoTIFF. Extra layers: object height and crown cover,
  forest productivity, land use (pasture, power lines, landscaped areas), low
  mountain forest. Licence CC0; "NMD, Naturvårdsverket can be given as the
  source" (source 2, section 1.3).
- **NMD 2023** also exists, newer than the version Ola named. Version 2.1
  (March 2026) has **54 classes in four levels**, 10 m, 16-bit, EPSG:3006,
  CC0, credit "NMD2023 v2.1, Naturvårdsverket"; it is "produced continuously,
  starting in southern Sweden", so it does not yet cover the whole country.
  Version 0.x has 16 classes and covers the areas with new forest laser data.
  Unlike NMD 2018, NMD 2023 is classified per pixel, not per image segment
  (source 3, sections 1.2 and 2).
- A generalised NMD 2018 existed as version 1.0 and is still on the server,
  but the current Swedish product description lists only the ungeneralised
  map (source 4, footnote 1).

One file checked by hand: Kronoberg county's NMD 2018 file is a single
9,792 × 15,471 `uint8` GeoTIFF, LZW-compressed, in 128-pixel tiles, cell
corners on whole 10 m coordinates, EPSG:3006 in its GeoKeys, NoData tag "0"
(tifffile, 2026-10-06). The English description says PackBits; this county
file is LZW.

## What a user would get, in numbers

One 10 km square near Ljungby in Lagan (x 430,000 to 440,000, y 6,295,000 to
6,305,000 in EPSG:3006; 86 % of it inside Lagan), NMD 2018 against Ola's
CORINE extract. Probe: `docs/research/nmd-landcover-probes/window_probe.py`,
run 2026-10-06 in the repository's venv (shapely 2.1.2, GEOS 3.13.1). The
square includes the town of Ljungby, so it is likely more fragmented than an
average square of Lagan; it is one sample.

| | CORINE 2018 | NMD 2018, raw | NMD, sieved to 1 ha | NMD, sieved to 25 ha |
|---|---:|---:|---:|---:|
| regions in the square | 92 after the clip | 47,344 | 1,926 | 87 |
| smallest region | 25.1 ha (whole polygons) | 0.01 ha (19,758 regions of one cell) | 1 ha | 25.2 ha |
| class boundary inside the square | 298 km | 4,602 km | 1,798 km | 479 km |
| classes present | 16 | 24 | | |

"Sieved" means each region under the area is merged into the neighbour class
it shares most boundary with, smallest first (GRASS `r.reclass.area`'s rule,
roughly; the probe's docstring says how). NMD sees what CORINE generalises
away: inland water is 0.9 % of the square in NMD against 0.2 % in CORINE, and
roads and railways are 4.8 % as their own
class. Forest totals agree roughly (NMD's forest classes and clear-fells
59 %, CORINE's forests and transitional woodland 61 %).

What it costs as constraint lines, on the square's north-west 3 km × 3 km
corner: vertices of the polygons after `shapely.coverage_simplify` (shared
boundaries simplified once, the square's frame kept) at four tolerances.
Vertices are counted per polygon ring, as CORINE's are, so a shared vertex
counts once per neighbour in both columns.

| on 9 km² | polygons | vertices, unsimplified | 5 m | 10 m | 20 m | 40 m |
|---|---:|---:|---:|---:|---:|---:|
| CORINE (as delivered) | 9 | 303 | | | | |
| NMD raw | 3,845 | 71,483 | 59,773 | 39,861 | 33,429 | 31,679 |
| NMD sieved to 0.25 ha | 625 | 38,704 | 29,170 | 16,178 | 9,994 | 7,376 |
| NMD sieved to 1 ha | 221 | 26,228 | 19,410 | 10,326 | 5,890 | 3,744 |
| NMD sieved to 25 ha | 8 | 5,334 | 4,244 | 2,336 | 1,412 | 972 |

`coverage_simplify` stands in here for the area-preserving reduction a
design would use (option A below); it does not keep areas, and the reduction's
own counts were not measured.

Scaled to Lagan, if this corner were typical (it is a single 9 km² sample):
Ola's Lagan run took CORINE's 274,426 vertices and made 863,897 triangles at
`--tolerance 10` (his `lagan_clc_glo30_tol10m_stats.md`). NMD sieved to 1 ha
and simplified at 10 m is 34 times CORINE's vertex count in the corner,
about 9 million ring vertices over Lagan, perhaps half of them distinct; at
about two triangles per vertex that is still on the order of ten million
triangles from the land cover alone. Sieved to 25 ha and simplified at 40 m it is about 3 times CORINE.

## The raster-to-polygon problem, and the options

`--features` takes vector polygons (16b), and 16c labels each triangle with
the one class it lies in. NMD is a raster, so using it that way means
turning it into polygons first. Raw cell outlines are staircases with a
vertex at every step (the tables above), and every vertex is a constraint the
mesh must keep. The project's input model holds: features are discrete
polygons, and tolerances set the resolution. The options:

- **A. Every class as constraint polygons.** Sieve on the raster, trace the
  cell boundaries once per shared arc, simplify each arc with its area kept,
  hand the arcs to `--features`. This is the method of
  `docs/research/raster-to-vector.md` ("What combines into a method for our
  case"): increment 22's area-preserving collapse, generalised to a coverage.
  Fidelity: the class areas are exact up to the sieve; the boundaries are
  within the simplification tolerance. Cost: the mesh grows with the land
  cover, not the terrain (above: from 3 to 34 times CORINE's constraints for
  Lagan, depending on sieve and tolerance). NMD's roads and railways are a
  special case: a road is one or two cells wide and connected over the whole
  catchment, so a sieve by area keeps it, and its two staircase sides become
  long thin polygons. It needs a width rule (Monmonier's, in the research
  note) or merging into its neighbours before tracing, and roads as lines
  belong to a road source, not to land cover.
- **B. A sieve and tolerance tied to `--tolerance`.** Option A, with the
  sieve area and the simplification tolerance derived from the mesh's own
  resolution instead of chosen by hand, so the land cover cannot ask for
  smaller triangles than the terrain does. The rule (for example a minimum
  area of a few times the tolerance squared) is a design question; in the
  tables the sieve cuts more vertices than the simplification does.
- **C. No constraints: class fractions per triangle.** Each triangle gets the
  exact areas of the 10 m cells inside it ("62 % spruce, 31 % arable, 7 %
  road"), and small fractions are carried to neighbouring triangles instead
  of dropped. This is **increment 26's design** for MapBiomas
  (`docs/increments/26-land-cover.md`, sections 1.1 to 1.6), ruled by Ola
  2026-10-04 and not yet built. The mesh is the one meshing without
  `--features` gives, so it should be smaller than today's 863,897
  triangles on Lagan, since CORINE's constraints go (not measured).
  Fidelity: no class boundary is an edge, but every class area is exact per
  triangle, which is what a hydrological model reads: HYPE, SMHI's model,
  divides each sub-basin into soil and land-use classes with an area each
  (source 6). A colour picture shows the dominant class per triangle, so it
  looks blockier than the raster.
- **D. Majority labels, refined.** Sample the class at each triangle and
  refine until each triangle is mostly one class, as `mesher` does for the
  Canadian Hydrological Model (research note, section 5). Fidelity sits
  between A and C; the class boundaries pull small triangles in anyway, and
  rasputin's refinement is driven by height error today, so it means a new
  refinement criterion in the C++ core.
- **C plus a few constraints**, which is increment 26 as ruled: fractions
  for every class, and only water traced, reduced and kept as constraint
  polygons (26d). For NMD that would be classes 61 (lakes and watercourses)
  and 62 (sea), perhaps also 2 (open wetland) if Ola wants bogs as edges.

**What this proposal recommends**: C plus water, that is, NMD as a second
source for increment 26's machinery, not a new pipeline. A alone is not
recommended for catchments the size of Lagan.

## Class mapping

NMD's classes do not map one-to-one onto CORINE's 44 (CORINE has
agricultural mosaics, NMD splits forest by species and by wetland; NMD's 42,
vegetated other open land, covers much of what CORINE splits into grassland,
heath and pasture).
A crosswalk was searched for (2026-10-06) and none was found from
Naturvårdsverket. So the proposal is **a code system of its own** (`nmd2018`,
and later `nmd2023`, which has different classes), with NMD's own codes as
`land_cover_code` and a natural-colour palette beside `rasputin palette
corine`. NMD's files carry a QGIS legend with a colour per class, which could
be the starting point. Today one mesh has one code system
(`src_python/tin_engine/cli.py@ed12512:1436`), so NMD and CORINE would not
go into one mesh, which is no loss if NMD replaces CORINE in Sweden. A
CORINE crosswalk could come later as an option for anyone comparing with
other countries.

## Fetching

Naturvårdsverket serves plain files over HTTPS
(`https://geodata.naturvardsverket.se/nedladdning/marktacke/`), with byte
ranges honoured (checked by listing the zip archives below with ranged reads):

- NMD 2018, the whole country: one 1.4 GB zip, whose map is one 1.57 GB
  GeoTIFF (deflate-compressed inside the zip, so it must be downloaded and
  unpacked whole, not read in ranges).
- NMD 2018, per county: 21 zips of 17 to 237 MiB in `NMD2018/bas_lan_ogen/`
  (Kronoberg 38 MiB, Jönköping 45 MiB, Halland 26 MiB, Skåne 36 MiB, Jämtland
  137 MiB, Västernorrland 73 MiB). Lagan lies in Jönköping, Kronoberg and
  Halland counties, perhaps also Skåne (not checked against county outlines).
- NMD 2023 version 2.1: one 2.7 GB zip, whose map is a 10.9 GB GeoTIFF.

None of these is a cloud-optimised GeoTIFF on a server, so the fetch cannot
read blocks remotely as `rasputin fetch` does for GLO-30 and ANADEM. A fetch
would be: `rasputin fetch nmd2018 --domain X` picks the counties the domain
meets (from county bounding boxes recorded once in the catalogue), downloads
their zips, unpacks the one GeoTIFF from each into the cache, checks its CRS
and cell size, and records the fetch date. Reading a window from the cached
file is the existing local block reader (`LocalTiffBlocks` and
`decode_window`, `src_python/tin_engine/io/cog.py@ed12512:73` and
`src_python/tin_engine/io/cog.py@ed12512:104`), which should handle a tiled
LZW file like Kronoberg's; not yet tried on NMD. The catalogue's source kinds
are COGs only today (`src_python/tin_engine/sources.py@ed12512:30`), so
"zipped files per region" is a new kind.

## Size

Counted in `CLAUDE.md` §2's unit, rough, assuming increment 26's PRs 26a to
26d have landed (they build the class windows, the fractions, the ledger and
the water tracing):

| PR | what | lines |
|---|---|---:|
| N1 | the zipped-files source kind, the NMD 2018 entry, county boxes, `rasputin fetch nmd2018` | about 200 |
| N2 | the `nmd2018` code system, class names and palette, water classes for 26d, credit and licence in the mesh record | about 150 |

If NMD came before increment 26, it would have to carry 26a to 26c itself
(about 890 lines estimated there), so it would be increment 26 with a Swedish
first source. Option A instead of C would add 26d's tracing and reduction for
every class, a width rule for roads and the tolerance rule (about 900 lines
over two PRs).

## Questions for Ola before a design

1. **What is the land cover for?** A picture that looks right, or a model
   input (class areas per triangle)? Default: model input, option C plus
   water.
2. **NMD 2018 or NMD 2023?** Default: NMD 2018 version 1.1 (whole country,
   25 classes, 8-bit), with NMD 2023 added once it covers the catchments Ola
   meshes.
3. **Which classes become edges?** Default: inland water and sea only (61,
   62), as increment 26 does for MapBiomas.
4. **Own codes or CORINE's?** Default: NMD's own code system and palette.
5. **Per county or the whole country?** Default: per county, chosen from the
   domain.
6. **When?** Default: after increment 26's 26a to 26c, as night work, since
   it reuses them.

## Prior art and novelty

The literature on turning a categorical raster into constraint polygons, and
on land cover in hydrological meshes, is in `docs/research/raster-to-vector.md`
(sections 1 to 6, and "Novelty"); the fractions and the ledger are increment
26's (its "Prior art: legacy and literature"). Nothing here is new beyond
those: NMD is another input to methods already chosen. No novelty is claimed.

## Sources

Each read on 2026-10-06.

1. Lantmäteriet, *Product description: CORINE Land Cover* (Swedish CLC),
   https://www.lantmateriet.se/contentassets/703ae721445f4398a2fd3c890472bbf8/e_clcshmi.pdf
   (text extracted with `pdftotext`; sections 1.1 and 2.1.1).
2. Naturvårdsverket, *National Land Cover Database (NMD)*, product description
   for NMD 2018,
   https://geodata.naturvardsverket.se/nedladdning/marktacke/NMD2018/NMD2018_ProductDescription_ENG.pdf
   (sections 1.1, 1.3, 2.3, 2.4, tables 1 and 2).
3. Naturvårdsverket, *Nationella marktäckedata 2023, Basskikt NMD2023 version
   2.1, Produktbeskrivning*,
   https://geodata.naturvardsverket.se/nedladdning/marktacke/NMD2023/Basskikt_v2_x/NMD2023_Produktbeskrivning_Basskikt_NMD2023_v2_1.pdf
   (sections 1.2, 2.1 to 2.4, and the terms of use).
4. Naturvårdsverket, *NMD 2018 basskikt, Produktbeskrivning* (Swedish),
   https://geodata.naturvardsverket.se/nedladdning/marktacke/NMD2018/NMD_Produktbeskrivning_NMD2018Basskikt.pdf
   (section 1.1 and its footnote 1).
5. The download server's directory listings,
   https://geodata.naturvardsverket.se/nedladdning/marktacke/ (NMD2018/,
   NMD2018/bas_lan_ogen/, NMD2023/), and the zip directories of the three
   whole-country archives, read with HTTP range requests.
6. SMHI, *HYPE model description*,
   https://hypeweb.smhi.se/wp-content/uploads/sites/11/2023/03/hype_model_description.pdf
   ("Land routines", "Basic assumptions": "A subbasin in HYPE is divided into
   classes depending on land use, soil type etc."; the class area comes from
   `slc_nn` and `area` in GeoData.txt).
