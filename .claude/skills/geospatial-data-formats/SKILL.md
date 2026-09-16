# Agent Skill: Geospatial Data Formats & CRS Management

Enforces zero-GDAL, high-performance geospatial data ingestion and processing using Shapely, PyProj, and clean metadata models.

## 1. Dependency Boundaries (No GDAL)
* **Strict Prohibition:** Never introduce `GDAL`, `OGR`, or `Fiona` as dependencies. 
* **Allowed Libraries:** Use python-native and lightweight C-bound libraries:
  * **Vector:** `shapely` (for geometry and predicates) and `geojson` or `ujson` (for fast JSON parsing).
  * **Raster/TIN:** `rasterio` (for GeoTIFF I/O) or direct binary parsing where applicable.
  * **Projections:** `pyproj` (for CRS definition and transformations).

## 2. Ingestion & Performance (Large Datasets)
* **GeoJSON:** Process large GeoJSON files using streaming or iterative parsing (e.g., `ijson`) to keep memory footprints low. Convert geometries directly into `shapely` objects.
* **TIN in TIFF:** Handle large GeoTIFF-based TIN surfaces efficiently. Use windowed reading (`rasterio.windows`) to process data chunks asynchronously without loading the entire raster into memory.
* **Custom XML Parsers:** Build light, iterative XML parsers using `xml.etree.ElementTree.iterparse` to ingest custom XML terrain formats line-by-line, enforcing strict schema validation.

## 3. Explicit CRS Enforcement
* **Mandatory Definition:** Every spatial dataset, file input, and API payload *must* explicitly define its Coordinate Reference System (CRS) at all levels.
* **Validation:** Reject any geospatial input that lacks a verifiable CRS (EPSG code, WKT, or PROJ string).
* **PyProj Integration:** Use `pyproj.CRS` to validate incoming coordinate systems and `pyproj.Transformer` (with `always_xy=True`) for all reprojections before passing data to the C++ core.

## 4. Output Architecture (Surfaces & Metadata)
* **2D Surface Generation:** Outputs must represent valid 2D planar subdivisions or meshed surfaces mapped to `shapely` geometry collections.
* **Granular Metadata:** Every output must encapsulate metadata at three explicit tiers via Pydantic models:
  * **Global Level:** CRS, overall bounding box (`bbox`), total element counts, and transformation history.
  * **Edge Level:** Structural topology tags (e.g., hard breakline, soft constraint, outer boundary).
  * **Cell/Triangle Level:** Reference DTM error metrics, slope, aspect, and source file identifiers.

## 5. Testing & Change Limits
* **Mocking I/O:** Test geospatial ingestion using small, hand-crafted GeoJSON snippets, micro-TIFFs, and minimal XML strings.
* **PR Constraint:** Maximum **700 LOC** per pull request (excluding test data).

