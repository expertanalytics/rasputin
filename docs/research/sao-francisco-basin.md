# The São Francisco basin as rasputin's target

Status: research note, 2026-09-27, main session. Ola's aim (2026-09-27):
"sub second is very good for smaller catchments, and of course faster is
better. But the important thing is huge catchments. I want to construct the
San Fransisco Basin in Brazil eventually, which is about the size of Norway.
And I think we can get 50m rasters for the DEM."

Facts below are from the sources listed at the end, checked by web search on
2026-09-27. Figures marked *computed* were computed here with pyproj; figures
marked *inferred* are extrapolations, not measurements.

## The basin

- **Area: 636,920 km²** (OAS). Norway is about 385,000 km² with Svalbard and
  Jan Mayen and about 324,000 km² without (from memory, not checked here), so
  the basin is roughly 1.65 to 2 times Norway.
- It drains parts of Minas Gerais, Goiás, Bahia, Pernambuco, Alagoas, Sergipe
  and the Federal District.
- Roughly latitude 7°S to 21°S and longitude 36°W to 48°W (approximate
  bounding box, used for the computations below). *Computed:* about
  1,280 km × 1,620 km in a projected CRS.

## Elevation data

| DEM | spacing | what it is | licence |
|---|---|---|---|
| **ANADEM** (ANA / UFRGS, 2024) | 30 m (1″) | Copernicus GLO-30 with vegetation bias removed for South America (Landsat-8, Sentinel-2, GEDI lidar). Mean bias 9.6 m (COPDEM) → 1.5 m; in forest 14.3 m → 0.4 m | free and open source per its authors; exact licence to be read from the repository |
| Copernicus GLO-30 | 30 m (1″) | a surface model (tree canopy and buildings included) | free (Copernicus licence) |
| Copernicus GLO-90 | 90 m (3″) | the same at 3″ | free |
| FABDEM | 30 m | Copernicus with forests and buildings removed (Bristol / Fathom) | CC BY-NC-SA 4.0: **non-commercial only** |
| IBGE MDE | 1:25,000 and 1:50,000 scale | national mapping agency, aerial restitution | coverage is partial |

- **No Brazil-wide 50 m DEM turned up.** A 50 m grid would most likely be
  ANADEM or Copernicus resampled from 30 m, or IBGE's 1:50,000 material where
  it exists. Where Ola's 50 m source is from is a question for Ola.
- **For hydrology, ANADEM is the obvious first choice.** It is a terrain model
  rather than a surface model, it was built for South America by the national
  water agency's partners, and it is free. It comes in MGRS tiles (GitHub) and
  on Google Earth Engine.
- All of these are in **geographic coordinates** (EPSG:4326, arc-seconds), not
  a projected grid. *Computed* (pyproj, WGS84): at 15°S, 1″ is 30.7 m
  north–south and 29.9 m east–west, so the cells are not square in metres.

## Other inputs for the basin

- **Basin outline:** ANA's *Base Hidrográfica Ottocodificada* (BHO) has
  Otto-Pfafstetter basins at several levels as GeoPackage, including the São
  Francisco. That gives the confining polygon directly, without auto-catchment.
  BHO's drainage lines (`GEOFT_BHO_TRECHO_DRENAGEM`) are candidate river
  polylines for 16b.
- **Land cover:** CORINE does not cover Brazil. **MapBiomas** (annual, 30 m,
  1985 onward, CC BY 4.0) is the Brazilian equivalent. It is a **raster**, not
  vector polygons like CORINE's GeoPackage, so 16b would need a
  raster-to-polygon step, or a separate path, for it.

## Size, against rasputin today

| grid | nodes inside the basin | float32 DEM |
|---|---:|---:|
| 30 m | 708 M | 2.8 GB |
| 50 m | 255 M | 1.0 GB |
| 90 m | 79 M | 0.3 GB |

*Computed* from the area. The bounding box holds more (about 830 M nodes at
50 m). Today's 1 m benchmark tile is 25.5 M nodes, so the basin at 50 m is
about **10×** the tile in nodes and at 30 m about **28×**.

*Inferred, not measured:* the scan visits every domain node about 5.6 times
over a refine (serial profile) at roughly 9 ns per visit on one thread. At
50 m that is on the order of 10-15 s of scan on one thread and a few seconds on
8. The triangle count, and so the serial phase and output size, depends on
relief and on the tolerance, not on the node count. It cannot be extrapolated
from the tile, which is a coastal fjord tile at 10 m with 1 m tolerance.
Measuring it on a real piece of the basin is the first thing to do.

The DEM fits in memory on a 32 GB machine at 30 m and at 50 m. Memory does not
force streaming at this size; parallel scaling and the serial phase set the
time.

## What this means for the roadmap

In roughly the order the basin needs them:

1. **A DEM in several tiles (ROADMAP gap 6).** ANADEM comes as MGRS tiles and
   Copernicus as 1°×1° tiles. The basin spans about 14° × 12°, so on the order
   of 150 1° tiles. Gap 6 today assumes aligned tiles in one CRS.
2. **Geographic DEMs and the computation CRS** (ROADMAP's "inputs in their own
   CRS"). rasputin's core expects a projected grid. Two routes:
   - resample the DEM onto a projected grid in Python before meshing (pyproj
     plus NumPy, chunked; no GDAL), or
   - mesh in a projected CRS while sampling the geographic grid.

   **One CRS for the whole basin is a real choice.** *Computed:* UTM zone 23S
   reaches a scale error of 1.07 % at the basin's east edge, and the Brazil
   Polyconic (EPSG:5880) 4.6 %. A Lambert conformal conic fitted to the basin
   (standard parallels 10°S and 18.5°S, central meridian 42°W) stays within
   **0.46 %** over the whole bounding box.
3. **21b's integer incircle does not apply to non-square cells.** QW2 needs
   `dx == dy`. On a geographic grid, or a projected grid resampled with
   non-square cells, refine falls back to the filtered kernel. That is correct
   but loses 21b's gain. Resampling to a square projected grid keeps it.
4. **Parallel refine (21d) matters more here than on the benchmark.** At
   10-28× the nodes, the scan and the serial phase both grow. Ola's
   determinism ruling of 2026-09-27 ("very large areas … we should not rely on
   bit-identical outputs") was made with exactly this case in mind.
5. **Domain decomposition** (21's option D, deferred to the mosaic work) is
   the natural fit for a basin of this size: split the basin into parts, mesh
   them in parallel, and keep shared boundaries. The tile-parallel greedy
   insertion prior art (Remote Sensing 12(3):437, 2020) is the reference.
6. **16b with the basin's own data:** BHO rivers as polylines, and MapBiomas
   land cover, which needs the raster step above.

## Open questions for Ola

1. Where does the 50 m DEM come from? If it is a resampling of ANADEM or
   Copernicus, meshing the 30 m source directly may be as cheap: refinement
   visits nodes, and the output size is set by the tolerance.
2. What tolerance is intended for the basin? It sets the triangle count more
   than anything else.
3. Commercial use? FABDEM is non-commercial only; ANADEM and Copernicus are
   the safe choices if that matters.

## Sources

- OAS, *São Francisco River Basin* brochure:
  http://www.oas.org/en/sedi/dsd/iwrm/Past_Projects/Documents/Sao_Francisco_Brochure.pdf
- "ANADEM: A Digital Terrain Model for South America", Remote
  Sensing 16(13):2321, 2024: https://www.mdpi.com/2072-4292/16/13/2321 ;
  data: https://github.com/HGE-IPH/anadem , https://hge-iph.github.io/anadem/
- Copernicus DEM: https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM ;
  https://registry.opendata.aws/copernicus-dem/
- FABDEM licence: https://www.fathom.global/insight/fabdem-download/
- IBGE digital elevation model:
  https://www.ibge.gov.br/en/geosciences/digital-surface-models/digital-surface-models/19081-digital-elevation-model.html
- ANA BHO: https://metadados.snirh.gov.br/geonetwork/srv/api/records/0c698205-6b59-48dc-8b5e-a58a5dfcc989
- MapBiomas Collection 9 (FAO catalogue):
  https://data.apps.fao.org/catalog/iso/416a2581-b5e2-43a5-bb93-a2fc26cd1d68
- Tile-parallel terrain simplification: https://www.mdpi.com/2072-4292/12/3/437
