# The São Francisco basin piece on ANADEM: triangle counts against tolerance (@perf, 2026-10-02)

The GLO-30 measurement of `../../2026-10-01/basin-piece/` re-run on **ANADEM**,
the terrain model it was meant to use. Everything else is the same: the piece
(BHO ottobasin 76949, upper and middle Rio das Velhas, 11,667.6 km², used as a
measurement domain only), the outline, the projection and resampling
(EPSG:31983, 30 m square grid, bilinear in the source's index space), the
tolerances (1, 2, 5, 10, 20, 50 m), the 200-box basin sample (same seed and
centres, so the boxes pair one to one) and the 30-box check inside the piece.
The Release `_core` is byte-identical to the GLO-30 run's: sha256
`eb5f06e88a7c…` in both `build_and_machine.json` files, built by
`tools/bench.py`'s `build()` from `5a57793`. Measured on AC power (`pmset -g
batt` before and after every block, in the JSON), Apple M1 Max (8P + 2E,
32 GiB), macOS 27.0, Python 3.14.7. No production code was changed.

## What changed from the GLO-30 scripts

- **Source.** ANADEM's own host still answers 403, so `prep_dem.py` reads
  OpenTopography's single continental COG
  (`https://opentopography.s3.sdsc.edu/raster/ANADEM/ANADEM_be/anadem_v1_compressed_COG.tif`,
  70.9 GB, BigTIFF, 512 × 512 Deflate tiles, NoData -9999). One 8 MiB range
  read gets every IFD (as `docs/increments/23-probes/anadem_cog.py` showed),
  parsed once per process. The script stops unless the full-resolution page's
  offsets number its blocks, and unless the GeoKeys say area-registered
  EPSG:4674. After that it range-reads only the blocks under each window.
- **CRS.** ANADEM is SIRGAS 2000 geographic (EPSG:4674), not WGS 84. The source
  window is written in its own CRS (`*_anadem_window_epsg4674.tif`).
  `prep_dem.py` takes the source CRS from its `SOURCE_EPSG` table (for
  ANADEM, `fetch` stops unless the COG's GeoKeys say EPSG:4674), and
  `resample` and `check` use that table too. Only `run_sweep.py`'s final check
  reads the CRS from the window file. Neither assumes EPSG:4326.
- `align_check.py` (new): checks that ANADEM's registration matches GLO-30's
  on the piece. `compare.py` (new): the side-by-side tables below.
- `run_sweep.py`, `sample_boxes.py`, `analyse.py`, `scaling.py`,
  `plant_check.py`: copied, with only the paths in their docstrings changed
  (plus the CRS read in `run_sweep.py`).

**The fetch was cheap.** The piece took 160 blocks (103 MiB) in 34 s. The 200
basin boxes took 573 blocks (324 MiB) in 325 s in total, about 1.6 s a box.
The whole sample ran in 14.5 min. No window had a NoData node.

## Results: ANADEM against GLO-30

From `compare.py` (`comparison.md`). The single-DEM tables, in the GLO-30
README's form, are in `analysis.md`.

### Triangles per tolerance

| tolerance | piece, GLO-30 | piece, ANADEM | ratio | basin, GLO-30 | basin, ANADEM (95 %) | ratio | paired box ratio, median (p10-p90) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 10,431,955 | 8,699,303 | 0.83 | 237 M | **175 M** (155-196 M) | 0.74 | 0.69 (0.54-0.83) |
| 2 m | 5,736,406 | 4,584,784 | 0.80 | 114 M | **84 M** (74-95 M) | 0.74 | 0.72 (0.58-0.83) |
| 5 m | 1,993,312 | 1,576,681 | 0.79 | 36 M | **28 M** (24-32 M) | 0.77 | 0.78 (0.65-0.86) |
| 10 m | 807,337 | 649,845 | 0.80 | 13.8 M | **11.0 M** (9.4-12.7 M) | 0.80 | 0.81 (0.72-0.91) |
| 20 m | 314,507 | 264,833 | 0.84 | 4.9 M | **4.2 M** (3.5-4.9 M) | 0.85 | 0.88 (0.77-1.00) |
| 50 m | 82,517 | 74,841 | 0.91 | 1.0 M | **0.9 M** (0.8-1.1 M) | 0.90 | 0.94 (0.79-1.04) |

The basin figure is the 200 boxes' mean density times the basin's area of
635,194.5 km² (BHO level 2). The interval covers the sampling error of the
mean only, as in the GLO-30 README's method, step 7.

### The final check's excess (interior source nodes off the mesh by more than the tolerance)

Each DEM is checked at its own source nodes. ANADEM's lattice
(0.000269495°) is about 3 % finer than GLO-30's (1/3600°), so the piece has
13.79 M interior ANADEM nodes against 13.0 M GLO-30 ones.

| tolerance | piece, GLO-30 | piece, ANADEM | basin nodes, GLO-30 | basin nodes, ANADEM | share of basin vertices, GLO-30 | ANADEM |
|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 22.30 % | 16.64 % | 52.4 M | 33.0 M | 44 % | 38 % |
| 2 m | 8.44 % | 5.52 % | 17.0 M | 9.9 M | 30 % | 24 % |
| 5 m | 1.42 % | 0.85 % | 2.62 M | 1.45 M | 15 % | 10 % |
| 10 m | 0.31 % | 0.18 % | 0.53 M | 0.30 M | 8 % | 5 % |
| 20 m | 0.06 % | 0.03 % | 0.12 M | 0.07 M | 5 % | 3 % |

The resampled grid itself, read back at the piece's source nodes (B6-style,
`data_meta/*_check.json`), is off by:

| | median | p99 | max | share over 1 m | over 2 m |
|---|---:|---:|---:|---:|---:|
| GLO-30 | 0.36 m | 3.14 m | 27.8 m | 18.8 % | 4.6 % |
| ANADEM | 0.28 m | 2.51 m | 31.2 m | 12.9 % | 2.3 % |

The boundary strip (within 30 m of the piece's outline) is as unchecked as on
GLO-30. Its worst source node is 14.1-36.0 m off on ANADEM over the six
tolerances, against 17.9-38.5 m on GLO-30 (`comparison.md`).

### The piece on ANADEM: time, memory, quality

| tolerance | triangles | vertices / grid nodes | refine | process wall | max RSS |
|---:|---:|---:|---:|---:|---:|
| 1 m | 8,699,303 | 33.6 % | 5.24 s | 8.27 s | 3.11 GiB |
| 2 m | 4,584,784 | 17.7 % | 2.68 s | 4.65 s | 1.69 GiB |
| 5 m | 1,576,681 | 6.1 % | 0.90 s | 1.98 s | 0.94 GiB |
| 10 m | 649,845 | 2.5 % | 0.42 s | 1.24 s | 0.78 GiB |
| 20 m | 264,833 | 1.1 % | 0.22 s | 0.96 s | 0.69 GiB |
| 50 m | 74,841 | 0.3 % | 0.09 s | 0.77 s | 0.62 GiB |

Medians of 5 runs, 10 threads. At every tolerance: worst angle 0.113°, max
degree 14, **0 constrained-Delaunay violations**, within tolerance, and the
control finds 0 of 12,957,257 grid nodes over tolerance. Max RSS is 0.56 GiB
+ 304 B per triangle (GLO-30: 0.55 GiB + 310 B). With Q9's 8.6 GiB canvas,
the basin's one-process memory floor comes to **58 GiB at 1 m, 32 GiB at
2 m, 16 GiB at 5 m** (GLO-30: 77, 41, 19 GiB). The basin's refine time,
extrapolated linearly, is 1.8 min at 1 m.

Refine thread scaling at 1 m (medians of 3, interleaved): 10.95 s at 1
thread, 5.22 s at 10, a **speed-up of 2.10**, against 2.05 on GLO-30 (12.95 s
to 6.32 s). That is the same input piece on a different DEM, so it is a
like-for-like comparison of the same binary on AC, not a baseline for an
increment.

The box method still holds: at every tolerance the piece's whole-catchment
density lies inside the 95 % interval of its 30 boxes (1 m: 745.6 against
786.8, 682.5-892.2).

## Surprises

1. **ANADEM's saving shrinks as the tolerance grows.** Over the basin it cuts
   triangles by 26 % at 1 m but only by 10 % at 50 m. The 2026-10-01 README
   expected fewer triangles "at small tolerances", which is now measured.
   The piece saves less than the basin at 1 m (17 % against 26 %). Its 2-10 m
   saving (20 %) is close to the basin's. Why the piece saves less at 1 m is
   not measured. Vegetation share (ANADEM removes the canopy bias) and Belo
   Horizonte's built-up area are candidate causes, untested.
2. **Even so, the basin at 1 m does not fit in one 32 GiB process.** The
   memory floor falls from 77 to 58 GiB at 1 m, and 2 m lands exactly at
   32 GiB. The GLO-30 conclusion stands: memory, not time, limits the basin.
3. **The final check is still not a touch-up at 1-2 m.** At 1 m it would find
   33 M basin source nodes, 38 % of the mesh's vertex count. Resampling ANADEM
   loses less than resampling GLO-30 (median 0.28 m against 0.36 m), but still
   about the size of a 1 m tolerance.
4. **Registration is right, and ANADEM sits lower.** On the piece's 30 m
   grid, ANADEM minus GLO-30 is smallest at zero shift (rms 3.53 m, against
   at least 5.6 m for a one-cell shift; `runs/align_check.json`), so the
   area-registered reading (node at the cell centre) is the correct one.
   ANADEM is 1.42 m lower on average and about 20 % smoother (mean
   |discrete Laplacian| 2.47 m against 3.07 m), as a canopy-removed surface
   should be. Median box relief barely moves (191 m against 192 m).
5. **The box sample does not need a tile mosaic.** OpenTopography's single
   COG removes the "tile selection by position" that the GLO-30 README listed
   under "Not done". The whole basin would be about 3,061 blocks (1.71 GiB,
   the probe's count). The two rates measured here put that at about 11 min
   (the piece's 160 blocks in 34 s, eight requests in flight) to 29 min (the
   boxes' 0.57 s a block, per-box overhead included). Neither is a measured
   basin fetch.

## Data and credits

Derived from **ANADEM**: "Agência Nacional de Águas e Saneamento Básico.
(2025). ANADEM: A Digital Terrain Model for South America. Distributed by
OpenTopography. https://doi.org/10.5069/G9736P4G." Licensed CC BY 4.0. As
OpenTopography's dataset acknowledgement asks, also cite
Laipelt, L.; Comini de Andrade, B.; Collischonn, W.; de Amorim Teixeira, A.;
Paiva, R.C.D.; Ruhoff, A. "ANADEM: A Digital Terrain Model for South America."
Remote Sens. 2024, 16, 2321. ANADEM is projected and resampled here, so
modified. ANADEM is itself derived from
Copernicus GLO-30 ("produced using Copernicus WorldDEM-30 © DLR e.V.
2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under
COPERNICUS by the European Union and ESA; all rights reserved"). The outlines
are the **BHO 2017** catchments of the Agência Nacional de Águas e Saneamento
Básico (ANA). None of these datasets is in the repository (`NOTICE.md`).

## Reproduce

From the repository root, the venv's python as `PY`,
`D=../rasputin_data/sao_francisco_piece`,
`S=../rasputin_scratch/basin-piece-anadem`,
`E=docs/benchmarks/2026-10-02/basin-piece-anadem`,
`E0=docs/benchmarks/2026-10-01/basin-piece`. The outlines come from
`$E0/fetch_bho.py` and the `curl` in `$E0/README.md` "Reproduce". The data
folder's `README.md` lists sources and licences, and `data_meta/sha256.txt`
has the inputs' checksums.

```bash
$PY -c "import sys; sys.path.insert(0,'tools'); import bench; from pathlib import Path; \
print(bench.build(bench.make_runner(), Path('.').resolve()))"          # Release build-bench/pkg
O=$D/bho2017_5k_76949_outline_epsg4674.geojson
$PY $E/prep_dem.py fetch anadem $O $D && $PY $E/prep_dem.py resample anadem $O $D \
  && $PY $E/prep_dem.py check anadem $O $D                              # ~35 s
G=$D/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif; W=$D/bho2017_5k_76949_anadem_window_epsg4674.tif
$PY $E/align_check.py $G $D/derived/bho2017_5k_76949_glo30_epsg31983_30m.tif > $E/runs/align_check.json
$PY $E/run_sweep.py velhas76949 $G $O $W $E/runs/velhas76949_anadem $S --repeats 5   # ~12 min
$PY $E/scaling.py $G $O $E/runs/velhas76949_anadem $S                   # ~2.5 min
$PY $E/plant_check.py $G $O $W $S 20 15 > $E/runs/checks.txt            # the control must fail when shifted
$PY $E/sample_boxes.py anadem $D/bho2017_level2_76_raw.geojson $D $E/runs/basin_boxes_anadem $S --n 200  # ~15 min
$PY $E/sample_boxes.py anadem $O $D/piece_boxes $E/runs/piece_boxes_anadem $S --n 30
$PY $E/analyse.py $E/runs/velhas76949_anadem/results.json $O $D/bho2017_level2_76_raw.geojson \
  $E/runs/basin_boxes_anadem/boxes.json $E/runs/piece_boxes_anadem/boxes.json \
  $E/runs/velhas76949_anadem/scaling.json > $E/analysis.md
$PY $E/compare.py $E0/runs $E/runs $D/bho2017_level2_76_raw.geojson > $E/comparison.md
```

## Checks behind the claims

- **The control can fail** (`runs/checks.txt`). On the piece at 20 m, as
  meshed, 0 grid nodes are over tolerance (max 19.99976 m). With the same mesh
  moved 15 m east, 20,101 are over.
- **Registration**: `align_check.py` exits non-zero unless the zero shift fits
  best. It did fit best (Surprises 4).
- **Header**: `prep_dem.py` stops unless the 8 MiB prefix gives the
  full-resolution page an offset for every block. The probe recorded that
  this check fails on a 1 MiB prefix.
- **Paired boxes**: `compare.py` asserts that the two samples have the same
  box centres.
- Not re-run here: determinism (the GLO-30 run's `runs/checks.txt` covers the
  same binary) and the brute-force point-location cross-check.

## Not done

- The whole basin's 3,061 blocks were not fetched; the box sample stands in
  for them.
- Large outputs: the binary piece meshes are in
  `../rasputin_scratch/basin-piece-anadem/`; the ASCII meshes were deleted
  after reading. The ANADEM window, grids and boxes are in
  `../rasputin_data/sao_francisco_piece/` (names contain `anadem`).

## Files

`prep_dem.py`, `run_sweep.py`, `sample_boxes.py`, `analyse.py`, `scaling.py`,
`plant_check.py` (adapted from the 2026-10-01 measurement), and the new
`align_check.py` and `compare.py`. These are one-off measurement scripts.
`analysis.md` and `comparison.md`: every table. `runs/`: raw JSON, the
piece's per-run logs and `--stats` files, and the box samples' logs in
`logs.tar.gz`. `data_meta/`: the window's and grid's metadata, the B6-style
check and checksums. `build_and_machine.json`: commit, build and machine.
