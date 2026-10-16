# Land-cover simplification probe (step 1), 2026-10-08

Ola: "we need to implement is a CORINE simplification, similar to what we do on
the outline of the auto-catchments. We get too many triangles (150k) and still
only 777m vtol." This probe measures what the switch already on master does
(`--features-tolerance`, PR #218), before step 2 (research and design) starts.
No production code changes here.

## Setup

- Code: master `b39426c0` (#218 merged), worktree `clc-simplify`, its own
  `.venv` and `_core` (`tin_engine.__file__` checked to point into it).
  shapely 2.2.0 on GEOS 3.14.1. Mac on AC power (`pmset -g batt`, in the log).
- Case: the German fused catchment (Isar above Krün + Walchensee + Loisach
  above Kochelsee), 1 195.3 km², 515 outline vertices, GLO-30 resampled to
  31 m in EPSG:25832, CORINE 2018 window (1 881 polygons, 18 classes inside).
  Same inputs as `rasputin_scratch/germany/isar-loisach/run_meshes_clc.sh`.
- Land-cover clean-up at its defaults: repair 1 m, same-class merge, outline
  rule 5 m. Varied: `--features-tolerance` (FT below) 0, 10, 30, 100 m;
  vertical tolerance 50, 20, 10 m with `--start-min-angle` 25 (the default) and
  0 (the lowest the CLI accepts); and the "minimal" mesh (`--tolerance 1e6
  --start-min-angle 0`, no height refinement). 28 runs, each 2.6 to 4.1 s.
- Files here: `probe.sh` (the 28 runs), `summarise.py` (writes
  `summary.md`), `cleanup_split.py` (writes `cleanup_split.txt`), and per run
  `ft<FT>_<run>_stats.md` (`--stats`) and `_errors.json` (Ola's
  `lands/scripts/mesh_error_stats.py`: every GLO-30 node inside the outline).
  Meshes and records: `rasputin_scratch/germany/isar-loisach/clc_simplify/`.

## Results

Full table: `summary.md`. Triangles, with RMS height error in brackets
(max error is always the tolerance, to 0.1 m; minimal: max error in metres):

| FT | minimal | 50 m, 25° | 50 m, 0° | 20 m, 25° | 20 m, 0° | 10 m, 25° | 10 m, 0° | land-cover vertices after clean-up |
|---|---|---|---|---|---|---|---|---|
| 0 | 150 697 (59.4; max 777) | 303 470 (9.2) | 161 272 (12.7) | 336 943 (5.5) | 216 973 (6.3) | 462 120 (3.4) | 374 775 (3.6) | 223 408 |
| 10 | 124 554 (59.4; max 775) | 257 713 (9.4) | 135 121 (12.8) | 292 730 (5.6) | 191 011 (6.3) | 422 235 (3.5) | 348 712 (3.6) | 184 903 |
| 30 | 73 205 (59.5; max 775) | 134 077 (10.7) | 84 027 (12.9) | 177 944 (6.1) | 141 178 (6.4) | 327 945 (3.6) | 304 346 (3.7) | 108 146 |
| 100 | 26 315 (62.5; max 812) | 59 081 (13.1) | 39 580 (14.0) | 120 198 (6.5) | 107 385 (6.6) | 292 543 (3.6) | 285 950 (3.6) | 41 066 |
| no land cover (NOTES.md, master `64481aae`) | 513 (545.0; max 1 381) | 23 912 (15.0) | | 96 189 (6.5) | | 267 748 (3.6) | | |

1. **The 777 m is not the land cover's doing.** The minimal mesh has no
   height refinement, so its error is whatever the triangles between border
   and outline vertices miss (mean signed error -10 m at FT 0 and 100;
   where the 777 m sits was not checked). FT 100 cuts its triangles 5.7 times
   (150 697 to 26 315) and its error stays (max 812 m, RMS 62.5 m against
   777 m, 59.4 m). Border vertices on the terrain do help against no land
   cover at all (513 triangles, max 1 381 m, RMS 545 m), but the first
   26 000 triangles of borders do that; the next 124 000 buy almost nothing.
   A height target needs `--tolerance`, not more borders.
2. **The 25° start pass is the largest cost with land cover at 50 m**: it
   doubles the mesh at FT 0 (303 470 against 161 272) and adds half at FT 100
   (59 081 against 39 580), by adding points along the borders. At 0° there
   are more slivers (table: 71 against 27 under 1° at FT 100, 50 m).
3. **At 10 m the terrain decides**: FT 100 saves only 24 % (285 950 against
   374 775 at 0°) and lands 7 % above the mesh without land cover.
4. **FT 100 with 0° at 50 m**: 39 580 triangles, 66 % above the plain mesh
   (23 912), against 12.7 times it (303 470) today at FT 0 and 25°.
5. **Slivers**: under 1° falls from 1 134 (FT 0) to 155 (FT 10) on the
   minimal mesh, then 53 and 89. The worst angle at 10 m is 0.0121° at FT 0,
   10 and 30 alike: one triangle with a 1 cm edge and a 65 m apex, class 332
   (bare rock), 988 m inside the outline, at 685 969 E, 5 251 307 N; the
   simplification does not touch it. Its cause was not checked.
6. **Land cover moved**: the stats' moved area (7 945 to 8 747 m²) is only
   the outline rule's, not the simplification's; the record has no figure for
   the latter. Measured here instead (below).

### Do classes keep their shape at FT 100?

Area per class, minimal mesh against the source clipped to the outline
(`summary.md`, exact polygon intersection): class totals stay within 2 %
(+0.5 % for 332 bare rock; -1.8 to -1.9 % for 121 industrial, 324 shrub and
335 glacier: the small classes lose). But **58.1 km², 4.9 % of the catchment, is
labelled with a different class than the source's** (FT 30: 8.5 km², 0.7 %;
FT 10: 0.6 km², 0.05 %). The totals hold because errors cancel, not because
the shapes stay.

The tolerance is not a distance bound (`cleanup_split.txt`, Hausdorff
distance per merged class before and after `coverage_simplify`, outline rule
left out): FT 10 moves a class border up to 33 m, FT 30 up to 147 m, FT 100
up to 325 m (median 227 m). GEOS documents the tolerance as the square root
of the Visvalingam-Whyatt area threshold, "roughly" a distance. CORINE's own
positional accuracy is 100 m or better (minimum mapping unit 25 ha, minimum
width 100 m), so FT 100 moves borders up to three times past the source's accuracy.

### Phases at 40 % or more

`features clip` is 40.6 to 50.8 % of the total in 23 of the 28 runs, and its
sub-row `features clip: clean-up` 40.3 to 42.6 % in 7 (FT 30 and 100, where
the mesh is small). The clean-up takes 1.1 to 1.4 s whatever FT is. A cProfile
of one run (FT 100, minimal; profiler overhead included): clean-up 1.21 s, of
which the outline rule (`snap_to_outline`) 0.53 s, `coverage_clean` 0.35 s,
clipping 0.18 s, `coverage_simplify` 0.12 s, the class merge 0.03 s. The
simplification itself is cheap; the outline rule is the largest part.

## What the existing switch achieves, and where it falls short

Achieves: a valid coverage (shared borders simplified once, no gaps, no
overlaps), the outer boundary kept, triangle counts down 5.7 times on the
minimal mesh and 4 to 5 times at 50 m (FT 100 against FT 0), slivers down, for 0.1 s.

Short of the outline-like reduction Ola asked for (the auto-catchment's
`reduce_ring`, `docs/increments/22-auto-catchment.md`, "The reduction"):

- **No area guarantee per class.** Visvalingam-Whyatt drops vertices by
  triangle area; class areas drift up to 1.9 % at FT 100 and 4.9 % of the
  ground changes class. `reduce_ring` keeps area exactly (APSC).
- **No distance bound.** The tolerance is an area threshold's square root;
  borders move 3 times the stated value. `reduce_ring` holds a band against
  the fine ring.
- **Terrain-blind and class-blind.** One planar tolerance for every border;
  nothing keeps vertices where a border crosses relief, or where a class is
  narrow or small (glacier, water, urban). The height refinement adds points
  back along borders where the terrain needs them, which is why 10 m meshes
  barely shrink.
- **Small pieces by planar area only.** The JTS documentation of the same
  class: rings smaller than the area tolerance (FT squared: 1 ha at FT 100)
  "are removed where possible", and the largest part of each input is kept.
  So a piece is dropped by size alone, at a scale tied to FT rather than to
  CORINE's 25 ha minimum mapping unit, and to whichever neighbour the
  coverage gives it. Not counted here.

## For step 2 (research), not a design

Questions to answer: can `reduce_ring`'s area-preserving collapse run per
shared border (between two junctions, ends fixed), where one collapse keeps
both neighbours' areas, with a band against the source border and a
topology check across the whole coverage? Should the band be measured in 3-D
on the border draped on the DEM, so vertices stay where the border crosses
relief? Should small pieces under a stated area be merged into a neighbour
first? And what should the default start angle be with land cover?

Literature to read (none of the papers read yet; venue and DOI found by
search today, or taken from `docs/increments/22-auto-catchment.md` for
Kronenfeld, Buchin, Saalfeld and Visvalingam; Estkowski and Mitchell from
memory, not checked):

- Kronenfeld, Stanislawski, Buttenfield, Brockmeyer 2020, "Simplification of
  polylines by segment collapse: minimizing areal displacement while
  preserving area", *Int. J. Cartography* 6(1):22-46,
  doi:10.1080/23729333.2019.1631535. APSC; its abstract says the areas of
  adjoining polygons are kept exactly, so it fits shared borders.
- Buchin, Meulemans, van Renssen, Speckmann 2016, "Area-preserving
  simplification and schematization of polygonal subdivisions", *ACM TSAS*
  2(1), doi:10.1145/2818373. The subdivision case directly.
- de Berg, van Kreveld, Schirra 1998, "Topologically correct subdivision
  simplification using the bandwidth criterion", *CaGIS* 25(4):243-257,
  doi:10.1559/152304098782383007. A distance band plus topology on a
  subdivision.
- Saalfeld 1999, "Topologically consistent line simplification with the
  Douglas-Peucker algorithm", *CaGIS* 26:7-18, doi:10.1559/152304099782424901.
- Estkowski and Mitchell 2001, "Simplifying a polygonal subdivision while
  keeping it simple", *SoCG 2001*; the hardness of the fewest-vertex version.
- Haunert and Wolff 2010, "Area aggregation in map generalisation by
  mixed-integer programming", *IJGIS* 24(12):1871-1897,
  doi:10.1080/13658810903401008. Merging small pieces into neighbours by class.
- Davis 2023, JTS/GEOS `CoverageSimplifier` (what `--features-tolerance`
  calls), and Visvalingam and Whyatt 1993, *Cartographic J.* 30:46-51,
  doi:10.1179/000870493786962263: the baseline.
- Garland and Heckbert 1995, "Fast polygonal approximation of terrains and
  height fields", CMU-CS-95-181: greedy insertion, for the terrain-aware side.
- CORINE Land Cover 2018 product page (Copernicus Land Monitoring Service):
  minimum mapping unit 25 ha, width 100 m, positional accuracy 100 m or
  better; sets a natural scale for the band.
