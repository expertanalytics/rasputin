---
name: geospatial-data-formats
description: Geospatial input in Python: rasters, vectors, GeoTIFF with tifffile, shapely, pyproj and CRS, and what crosses the I/O boundary.
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
