---
name: geospatial-data-formats
description: Geospatial ingestion and CRS handling in Python. Use when reading or writing rasters and vectors, parsing GeoTIFF with tifffile, working with shapely geometry, pyproj transformations, CRS metadata models, or deciding what crosses the Python/C++ I/O boundary.
---

# Agent Skill: Geospatial Data Formats & CRS Management

Enforces high-performance geospatial data ingestion and processing using Shapely, PyProj, and clean metadata models.

## 1. Dependency Boundaries
* **Allowed Libraries:** Use python-native and lightweight C-bound libraries:
  * **Vector:** `shapely`.
  * **Raster/TIN:** direct binary parsing, or a pure-Python reader such as `tifffile`. The project's own reader will live under `src_python/tin_engine/io/`.
  * **Projections:** `pyproj` (for CRS definition and transformations).

## 2. Ingestion & Performance (Large Datasets)
* **TIN in TIFF:** Handle large GeoTIFF-based TIN surfaces efficiently. Read by window/tile so chunks are processed asynchronously, without loading the entire raster into memory.
* **Custom XML Parsers:** Build light, iterative XML parsers using `xml.etree.ElementTree.iterparse` to ingest custom XML terrain formats line-by-line, enforcing strict schema validation.

## 3. Explicit CRS Enforcement
* **Mandatory Definition:** Every spatial dataset, file input, and API payload *must* explicitly define its Coordinate Reference System (CRS) at all levels.
* **Validation:** Reject any geospatial input that lacks a verifiable CRS (EPSG code, WKT, or PROJ string).
* **PyProj Integration:** Use `pyproj.CRS` to validate incoming coordinate systems and `pyproj.Transformer` (with `always_xy=True`) for all reprojections before passing data to the C++ core.

## 5. Testing
* **Mocking I/O:** Test geospatial ingestion using small, hand-crafted GeoJSON snippets, micro-TIFFs, and minimal XML strings.

