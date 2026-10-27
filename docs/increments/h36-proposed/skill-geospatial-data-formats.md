---
name: geospatial-data-formats
description: Geospatial ingestion and CRS handling in Python. Use when reading or writing rasters and vectors, parsing GeoTIFF with tifffile, working with shapely geometry, pyproj transformations, CRS metadata models, or deciding what crosses the Python/C++ I/O boundary.
---

# Geospatial Data Formats & CRS

## 1. Libraries
`shapely` for vectors; `tifffile` or direct binary parsing for rasters, the
project's reader under `src_python/tin_engine/io/`; `pyproj` for CRS.

## 2. Ingestion
Read large rasters by window or tile, never whole. Parse XML iteratively
(`xml.etree.ElementTree.iterparse`) with strict schema validation.

## 3. CRS
Every dataset, file and payload states its CRS; reject input without a
verifiable one. Validate with `pyproj.CRS`; reproject with
`pyproj.Transformer(always_xy=True)` before the data reaches C++.

## 5. Testing
Small hand-made GeoJSON, micro-TIFFs and minimal XML.
