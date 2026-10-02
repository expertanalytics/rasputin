# What vertical tolerance does ANADEM's accuracy support?

Status: research note, 2026-10-02, `@architect`. Input to Ola's open decision
B11, the basin tolerance (`docs/increments/23-basin-scale.md`, B11). Research
only: no design, and no recommendation in place of Ola's choice. Every source
below was read on 2026-10-02 unless marked **unverified**; figures marked
*computed* were worked out here and are not measurements.

## 1. ANADEM's own validation (Laipelt et al. 2024)

Read in full from the authors' institutional copy (UFRGS LUME, the published
MDPI version, 17 pages; MDPI itself refused automated download).

- **Method.** Copernicus GLO-30 minus a vegetation bias predicted per pixel by
  gradient tree boosting (200 trees, Huber loss) trained on GEDI ground
  returns and Landsat-8/Sentinel-2 indices, in 93 tiles of 5° × 5°. Afterwards
  a **3 × 3 focal median filter** "to reduce noise"; water bodies keep
  Copernicus's values. No urban correction (§2.8, §4).
- **Reference data.** ICESat-2 ground heights only: 20,000 random points per
  tile, 600,000 in the model comparison (§2.9, §3.4). **No airborne lidar, no
  GNSS, no local reference** was used; the conclusions name lidar validation
  as future work. GEDI trained the model and is not a validation set.
- **Overall (Table 3, against ICESat-2):** ANADEM bias 1.50 m, RMSE 6.99 m,
  STD 6.79 m, median 0.75 m. Copernicus: bias 9.56 m, RMSE 12.40 m.
  FABDEM: 1.76 m, 6.81 m.
- **By land cover (Table 4, MapBiomas classes), ANADEM RMSE / STD:**
  croplands 2.26 / 2.23 m, grasslands 3.39 / 3.36 m, urban 3.78 / 3.31 m,
  pastures 4.70 / 4.48 m, savanna 4.75 / 4.12 m, forests 7.09 / 7.08 m.
  (The text calls these "percentage RMSE"; the figures it quotes are the RMSE
  column.)
- **By tree cover (Table 2):** RMSE 5.82 m at 0-20 % cover, 5.27 m at 21-40 %,
  rising to 8.88 m at 61-80 %.
- **By slope (§3.3, Figure 7):** the error "showed a proportional increase
  with increasing slope". Under 3° it fell from 3.7 m (Copernicus) to 0.57 m
  (ANADEM); above 21°, from 9.5 m to 5.3 m. Figure 7's quantity is not named;
  its aspect panel gives 1.30-1.63 m for ANADEM, which matches the overall
  bias rather than the RMSE, so **0.57 m is most likely a mean (bias), not an
  RMSE** (inferred from the figures, not stated by the paper).

What it does not report: any relative (point-to-point, slope) accuracy, which
is what flow routing depends on.

## 2. The Copernicus GLO-30 accuracy it builds on

Copernicus DEM Product Handbook, v5.0 (29.11.2022), §1 table and §2.1-2.2:
absolute vertical accuracy **< 4 m LE90**, relative vertical accuracy
**< 2 m LE90 for slopes ≤ 20 %, < 4 m for slopes > 20 %** (point-to-point
within 1° × 1°), validated by DLR against ICESat GLAS; a global mean, "local
deviations can occur". Table 12: mean LE90ABS 1.92 m outside Greenland and
Antarctica. These are radar-surface figures: on open ground they hold; under
canopy the surface is the canopy (Copernicus's 14.3 m forest bias in
Laipelt's Table 4).

*Computed:* for normal errors, LE90 = 1.645 σ, so the 2 m relative figure is
σ ≈ 1.2 m per point, or ≈ 0.9 m if the LE90 is of the difference of two
independent points. Whether ANADEM keeps, improves or worsens Copernicus's
relative accuracy (the ML correction varies pixel to pixel; the median filter
smooths) is **not measured by anyone I found**.

One independent Brazilian test exists: Salgado et al. 2026 compare six DEMs
including ANADEM against lidar for two Brazilian dam-break studies. Only the
abstract was readable: Copernicus and FABDEM had the best "altimetric
accuracy"; ANADEM and FABDEM over-estimated flooded area. **Its numbers are
unverified** (paywall).

## 3. Our own measurements

- **The resampling floor** (`docs/benchmarks/2026-10-02/basin-piece-anadem/README.md`,
  pinned at `aa5dcb1`): projecting ANADEM to the 30 m EPSG:31983 grid and
  reading it back at the source nodes is off by median 0.28 m, p99 2.51 m;
  12.9 % of nodes over 1 m, 2.3 % over 2 m. Our own pipeline step is already
  at the size of a 1 m tolerance.
- **Triangles** (same README, basin by the 200-box estimate): 175 M at 1 m,
  84 M at 2 m, 28 M at 5 m, 11.0 M at 10 m: 2-3 times more at each step
  from 10 m down to 1 m.
- **With the guaranteed tolerance** (15c-2's final check, branch
  `worktree-15c-2` at `999559a`, `docs/benchmarks/2026-10-02/15c-2-acceptance/README.md`,
  section (b); not merged, and its verdict is a regression on an unrelated
  crash): the Velhas piece needs 1.595 × the triangles at 1 m, 1.18 × at
  5 m, 1.09 × at 10 m (2 m not run). *Computed* by applying those ratios to
  the basin estimate: about 280 M at 1 m, 33 M at 5 m, 12 M at 10 m. At 1 m
  the piece's mesh holds about half as many vertices as there are source
  nodes; most of what a 1 m tolerance adds is spent following sub-metre
  structure, which section 1 says is within the DEM's error.
- ANADEM is about 20 % smoother than Copernicus on the piece (same README,
  Surprise 4), consistent with the median filter and canopy removal.

## 4. Literature on tolerance relative to DEM error

- **Vivoni et al. 2004** (J. Hydrol. Eng. 9(4):288-302) mesh with a
  maximum-deviation tolerance (the drop heuristic, ArcInfo Latticetin, as
  ours is in spirit). At 4 m tolerance the TIN-to-DEM RMSE stayed between 1.2
  and 2 m over basins of very different relief. With the same 4 m tolerance a
  noisier DEM (SRTM) gave a **denser TIN in a flat canyon bottom** than a
  smoother one (USGS): a tolerance near the noise level spends triangles
  tracing noise.
- **Vivoni et al. 2005** (Hydrol. Process. 19(11):2101-2122), on a 64 km²
  catchment with tRIBS: TIN RMSE is well below the tolerance (1 m → RMSE
  0.53 m; 6.8 m → 2.3 m). The model response was **weakly sensitive to
  refinement finer than the resolution it was calibrated at and strongly
  sensitive to coarsening below it**; what mattered was resolving the
  near-stream variable source areas, not the global tolerance. They place the
  break near 10 % of the DEM nodes kept (6.8 m tolerance in that basin).
- I found **no published rule** that sets a TIN tolerance as a multiple of
  DEM error. "Do not mesh below the DEM's noise" is a reasonable inference
  from the above, not a cited result. Lee (1991, IJGIS), the drop-heuristic
  source, is known here only through Vivoni's citation (**unverified**).

*Computed*, error budget: from Vivoni's figures a tolerance T adds an RMSE of
roughly 0.3-0.53 T. Combined in quadrature with ANADEM's absolute error on
open land (pasture/savanna, σ ≈ 4.5-4.7 m), the total grows by at most 1 % at
1 m, 1-3 % at 2 m, 5-16 % at 5 m and 19-55 % at 10 m. Against a relative
(slope-relevant) error of σ ≈ 0.9-1.2 m (section 2, assumed to carry over to
ANADEM), it grows by 3-16 %, 12-55 %, 60-211 % and 2.7-6 times.

## 5. What the evidence says about 1, 2, 5 and 10 m

For hydrological use in the São Francisco basin, with the decision Ola's:

- **1 m** sits inside ANADEM's error everywhere (smallest RMSE 2.26 m, in
  croplands), inside the ~1 m relative noise, and at the size of our own
  resampling error. It costs ~280 M triangles and mostly buys noise. The
  evidence does not support it as a terrain-accuracy requirement.
- **2 m** is about the relative noise level and below every land-cover RMSE.
  It adds almost nothing to the absolute error and somewhat to the relative.
  It is the finest value the data can arguably justify, at ~3 × the triangles
  of 5 m.
- **5 m** is about one absolute RMSE on open land, so it barely moves the
  absolute error (5-16 %). It is coarser than the relative noise, so local
  slopes on gentle terrain are smoothed. That matters on the flats, where
  ANADEM's mean error is ~0.6 m and channels and floodplain relief are a few
  metres.
- **10 m** clearly degrades the surface relative to the data (absolute error
  up 19-55 %), and on flat terrain it can merge valleys a few metres deep.

Uncertainty: every accuracy figure is continental (ICESat-2 over South
America), not per basin. Nothing measures ANADEM's relative accuracy. The
basin's land-cover mix is not measured here (MapBiomas is a 23e input). The
error-budget factor comes from two US basins on other DEMs. On the evidence,
the defensible band is **2 to 5 m**. Where the basin is flat, the hydrology
depends more on getting the channels into the mesh (rivers as constraints,
23e) than on the global tolerance, so the choice within that band is likely
to matter less than whether channels are constrained. That last point is an
inference from Vivoni 2005, not measured here.

## Sources

- Laipelt, L.; Comini de Andrade, B.; Collischonn, W.; de Amorim Teixeira, A.;
  Paiva, R.C.D.; Ruhoff, A. "ANADEM: A Digital Terrain Model for South
  America." *Remote Sens.* 2024, 16(13), 2321.
  [doi:10.3390/rs16132321](https://doi.org/10.3390/rs16132321); full text read
  at <https://lume.ufrgs.br/bitstream/10183/280581/1/001206574.pdf>.
- OpenTopography, ANADEM dataset page (vertical datum EGM2008, 30 m; artefact
  disclaimer): <https://portal.opentopography.org/datasetMetadata?otCollectionID=OT.082025.4674.1>.
- Copernicus DEM Product Handbook, GEO.2018-1988-2, v5.0, 29.11.2022:
  <https://dataspace.copernicus.eu/sites/default/files/media/files/2024-06/geo1988-copernicusdem-spe-002_producthandbook_i5.0.pdf>.
- Vivoni, E.R.; Ivanov, V.Y.; Bras, R.L.; Entekhabi, D. "Generation of
  Triangulated Irregular Networks Based on Hydrological Similarity."
  *J. Hydrol. Eng.* 2004, 9(4), 288-302,
  doi:10.1061/(ASCE)1084-0699(2004)9:4(288);
  <http://vivoni.asu.edu/pdf/VivoniJHE2004.pdf>.
- Vivoni, E.R.; Ivanov, V.Y.; Bras, R.L.; Entekhabi, D. "On the effects of
  triangulated terrain resolution on distributed hydrologic model response."
  *Hydrol. Process.* 2005, 19(11), 2101-2122,
  [doi:10.1002/hyp.5671](https://doi.org/10.1002/hyp.5671);
  <http://vivoni.asu.edu/pdf/VivoniHP2005.pdf>.
- Salgado, S.R.T.; Carvalho, E.M.S.; Viseu, M.T.; de Oliveira, O.F.
  "Comparative Analysis of Global DEMs for Dam-Break Flood Modeling and
  Inundation Mapping." *Water Resour. Manag.* 2026, 40(3),
  [doi:10.1007/s11269-025-04427-9](https://doi.org/10.1007/s11269-025-04427-9).
  Abstract only (via OpenAlex); figures **unverified**.
- Lee, J. (1991), comparison of TIN-building methods, IJGIS: cited only
  through Vivoni et al. 2004; **unverified**.
