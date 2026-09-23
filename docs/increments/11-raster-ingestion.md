# Increment 11: raster ingestion, the Python side

**Status.** Design. Written by `@architect` before `@tester` is spawned, per
`docs/increments/README.md` step 1. No production code and no tests here.

**Input.** `docs/increments/11-raster-ingestion-prior-art.md`, the
`@migration-expert` report, already reviewed and corrected. That file holds the
evidence; this file holds the rulings. Where they disagree, this file wins, and
where this file overturns a recommendation it says so by section number.

**Closes.** Half of `ROADMAP.md` MVP gap 1: nothing decodes a GeoTIFF. After
this increment a GeoTIFF becomes a validated Python object. It still does not
reach `_core`; ruling 2 says why, and names the increment that finishes the job.

---

## 1. The rulings, in one list

1. **Decode only. No C++, no bindings, no `_core` import.** The adapter and the
   zero-copy view are increment 12.
2. **`tin_engine/io/geotiff.py` takes a binary stream, not a path.**
3. **The reader's job is refusal.** Eighteen named refusals, section 5.
4. **Area-registered files are converted, not refused** — the project's only DEM
   fixture is one. An *absent* registration key is refused.
5. **Use `tifffile`'s own GeoKey decoding.** Do not port `GeoKeysInterpreter`.
6. **CRS comes from `ProjectedCSTypeGeoKey` alone**, resolved through
   `pyproj.CRS.from_epsg`. No proj4 reassembly, no free-text ellipsoid regex.
7. **Projected-and-metre is tested on the constructed CRS**, never on which
   GeoKeys are present.
8. **`always_xy=True` is a hard requirement**, enforced by a grep test that
   lands in this increment, before the first transformer exists.
9. **NoData: absent tag means no sentinel.** No value is ever guessed. An
   explicit caller override exists and contradicting the tag is refused.
10. **The mosaic walk is a separate increment.** So is `GeoPolygon`.
11. **`crop_image_to_polygon` is not ported.** Its index arithmetic is, in
    increment 12, where the window has a caller.
12. **`project_structure.md`'s CRS paragraph is wrong and is rewritten here.**
13. **`pyproject.toml`'s codec comment is false as installed** and is corrected
    in this increment's PR.

---

## 2. Ruling 1 and 2: what this increment is, and what it is not

`project_structure.md` names `tin_engine/raster.py` as the only adapter from
decoded data into `_core`. This increment does not write it.

There are no raster bindings today. Check it:

```
$ grep -in raster bindings/core.cpp src_python/tin_engine/_core.pyi
(no output)
```

Writing the adapter therefore means writing `include/terrain/raster/view.hpp`,
and `project_structure.md` is explicit that `RasterView` must bring contiguous
row access into the `RasterSource` concept **in the same change** — which edits
`raster.hpp`, `sample.hpp` and every model and test double behind them. That
change has no caller until refinement exists, and a concept requirement added
for no consumer is a requirement nobody can check.

The alternative, an adapter that copies into the existing owning `Raster<T>`,
contradicts the boundary contract's "numpy owns the buffer" and would be thrown
away. So:

- **Increment 11**: `tin_engine/io/geotiff.py` and `tin_engine/io/models.py`.
  Pure Python. Imports `tifffile`, `pyproj`, `numpy`, `pydantic`, nothing
  first-party, and **never `_core`**. Testable with no compiled extension in
  the process — the same shape as `io/ply.py`, which is why increment 10 was
  cheap.
- **Increment 12**: `view.hpp` plus the row-access concept change, the pybind11
  buffer binding, `_core.pyi`, and `tin_engine/raster.py`. Also the index
  window (`window_for`, and the surviving arithmetic of
  `crop_image_to_polygon`), which by then has a caller.

Say this plainly in the roadmap row: gap 1 is half closed, and a decoded tile
that reaches nothing is still progress, because everything downstream needs the
tile and nothing downstream exists.

### `geotiff.py` takes a stream

`decode_dem(source: BinaryIO) -> DemTile`. `tifffile.TiffFile` accepts a file
object, so nothing is lost, and three things are gained: the module opens no
file and resolves no path, `@tester` builds fixtures as `io.BytesIO` micro-TIFFs
with `tifffile.imwrite`, and path resolution stays where `cli.py` already keeps
it. `io/ply.py` is the precedent on the write side; this is its mirror.

---

## 3. Ruling 12: the CRS paragraph

`project_structure.md`'s `raster` section justified Python-side decoding partly
with *"CRS interpretation means PROJ — GDAL's own dependency."*

**That is false.** PROJ is a standalone library. GDAL depends on PROJ, not the
other way round. `CLAUDE.md` §2 lists PyProj in the core stack, so this project
already depends on PROJ and always did — the legacy used `pyproj` for every
transformation and never linked GDAL (prior art §2). The sentence argued that
we must avoid the thing we are built on.

The conclusion it was defending is still right, for a different reason, and the
paragraph tangled two rules that need separating:

- **Where a CRS string may live.** The C++ core has no use for one. Nothing in
  `geometry.hpp`, `raster.hpp` or `sample.hpp` reads a CRS, nothing can, and a
  field nothing reads is a field that rots — the legacy's proved it by carrying
  a proj4 string into the mesh for no consumer. This is an argument about
  dependency surface and unused state. It is not an argument about PROJ.
- **What the numbers must mean.** This is the load-bearing one, and it is the
  one the evidence supports. Nothing on either side of the boundary is
  unit-aware. Measured, at this tree's `.venv` (pyproj 3.8.0, PROJ 9.8.1), for
  a 0.0002777° cell at 60°N reprojected to EPSG:32633:

  ```
  $ .venv/bin/python -c "
  import pyproj, math
  t = pyproj.Transformer.from_crs('EPSG:4326', 'EPSG:32633', always_xy=True)
  d = 0.0002777
  x0, y0 = t.transform(8.5, 60.0)
  x1, y1 = t.transform(8.5 + d, 60.0)
  x2, y2 = t.transform(8.5, 60.0 + d)
  ew, ns = math.hypot(x1-x0, y1-y0), math.hypot(x2-x0, y2-y0)
  print('cell metres: dx=%.3f dy=%.3f ratio=%.3f' % (ew, ns, ns/ew))"
  cell metres: dx=15.514 dy=30.977 ratio=1.997
  ```

  `dx` is east-west, `dy` north-south.

  A 2:1 anisotropy that no code anywhere can see. The refinement error metric
  would then compare a horizontal distance in degrees against a vertical one in
  metres. `sample.hpp`'s bilinear weights `tx` and `ty` would be computed in
  degrees and applied to metres.

So the rule is not "no CRS in C++, because PROJ". It is two rules: **no CRS
string crosses the boundary, because nothing there reads one**, and **every
number that crosses is metres in a projected CRS, because nothing there can
tell**. The second is what the reader enforces.

The paragraph is rewritten in `project_structure.md` on this branch, in the
`raster` section. Both halves are stated, each with its own justification, and
the PROJ sentence is gone.

### The geographic refusal has never been implemented anywhere

Prior art §7 established this and it is worth restating as a ruling rather than
a finding: the legacy never rejected a geographic CRS. A geographic GeoTIFF
raised `CRSError` during proj4 assembly and died before reaching the mesher.
The defence was a typo (`_GeoGraphicTypeGeo` for `_GeographicTypeGeoKey`).

So `refuses_geographic_crs` is not a port. It is the first implementation of a
rule this project has asserted since before it had a raster module. `@tester`
should treat it as greenfield and `@reviewer` should not look for prior art.

---

## 4. Ruling 4: area-registered files are converted

Prior art §4.2 and `project_structure.md` both record that the legacy defined
`GTRasterTypeGeoKey` (1025) and read it nowhere, so an area-registered file
lands half a cell off. The brief proposed refusing such files.

**Overturned.** The project's only DEM fixture is area-registered. Measured:

```
$ .venv/bin/python -c "
import tifffile
m = tifffile.TiffFile('tests/fixtures/dem_archive/7908_3_10m_z33.tif').geotiff_metadata
print(repr(m['GTRasterTypeGeoKey']))"
<RasterPixel.IsArea: 1>
```

`repr`, not `print`: tifffile returns an `IntEnum`, and printing one gives a
bare `1` — which is the value `RasterPixelIsArea`, but a reader checking this
would see a number and have to look it up.

Refusing area registration would refuse the only DEM in the tree, and with it
Kartverket's 10 m product, which is what this engine is for. A rule whose first
act is to reject the reference input is the wrong rule.

The defect was never the registration. It was the silent half cell. And the
conversion is exact, not a heuristic: an area-registered file's tie point names
the upper-left **corner** of the upper-left pixel, so the node grid is that
corner shifted inward by half a cell in each axis.

```
x_min = tie_x + delta_x / 2
y_max = tie_y - delta_y / 2      # area-registered only
```

`cols`, `rows`, `delta_x` and `delta_y` are unchanged. The node extent is then
inset by half a cell on each side and spans `(n - 1)` spacings, which is what
`RasterGeometry` means.

The fixture corroborates the arithmetic. Its tie point is
`(0, 0, 0, 799745.0, 7950255.0, 0.0)` at a 10 m scale, and the shift gives
`x_min = 799750.0`, `y_max = 7950250.0` — both exact multiples of the cell
size, which a corner offset by half a cell is and a node grid mis-shifted by
half a cell is not.

**What is refused is not knowing.** If `GTRasterTypeGeoKey` is absent, or is
32767 (user-defined), the reader raises and names the tag. Half a cell is the
entire quantity in dispute and there is no way to recover it from the data.
Guessing a default here is exactly the legacy's failure mode wearing a
specification citation. `RasterMeta` records which convention the file declared
and the shift applied, so the decision is auditable after the fact.

---

## 5. Ruling 3: the refusals

Prior art §8.1: seven of eleven legacy defects are "returned something instead
of raising". Each entry below is one named test. The names are the test names.

**Container and georeferencing**

1. `refuses_missing_georeferencing` — `ModelTiepointTag` (33922) or
   `ModelPixelScaleTag` (33550) absent. The legacy defaulted them to zeros and
   `(1.0, 1.0)` and produced a unit-spaced raster at the origin (prior art
   §4.5). This is the one place the prior art calls refusal unambiguously right.
2. `refuses_model_transformation` — `ModelTransformationTag` (34264) present. A
   rotated or sheared transform cannot be represented by `RasterGeometry`, and
   the legacy never looked. Refuse even if a tie point and scale are also
   present: a file carrying both is contradicting itself.
3. `refuses_multiple_tiepoints` — `ModelTiepointTag` with more than six values.
   That is a ground-control-point list, not an affine placement.
4. `refuses_nonzero_pixel_scale_z` — `ModelPixelScaleTag` ScaleZ is not zero.
   The file is asserting a vertical scaling this reader drops. Dropping it
   silently is the legacy's whole class of defect.
5. `refuses_nonpositive_pixel_scale` — ScaleX or ScaleY not finite and
   positive. `RasterGeometry`'s constructor demands positive spacings, and a
   negative ScaleY in the wild means a file that has already flipped north-up
   somewhere upstream.
6. `refuses_unknown_raster_type` — `GTRasterTypeGeoKey` (1025) absent or 32767.
   Section 4.
7. `refuses_degenerate_shape` — fewer than two rows or two columns. `bilinear`
   returns `nullopt` for a raster too small to hold a 2×2 neighbourhood, so
   such a tile can be decoded and then never sampled. Refuse at decode.
8. `refuses_ambiguous_pages` — more than one page where the extra pages are not
   flagged reduced-resolution. Page 0 is used; being silently one of several
   full-resolution images is not acceptable.
9. `refuses_multi_sample` — `SamplesPerPixel` is not 1. A DEM has one band.
10. `refuses_unsupported_dtype` — anything outside the promotion table in
    section 7.
11. `refuses_missing_codec` — the page's compression needs `imagecodecs`. The
    message names the compression scheme and the `codecs` extra. A bare
    `KeyError` from inside `tifffile` is not a diagnostic; see section 8.

**CRS**

12. `refuses_missing_crs` — no `ProjectedCSTypeGeoKey` (3072), or its value is
    32767, or it is not a resolvable EPSG code. Ruling 6: there is no fallback
    path. The legacy had two, one of them a regex over human prose.
13. `refuses_geographic_crs` — the constructed CRS has `is_projected` false.
    Section 3.
14. `refuses_non_metre_linear_unit` — any horizontal axis has
    `unit_conversion_factor != 1.0`; or `ProjLinearUnitsGeoKey` (3076) is
    present and is not 9001; or the two disagree. The legacy's version of this
    was `{9001: "m"}[value]`, a `KeyError` with no message (prior art §4.8).
    Also `refuses_non_metre_vertical_unit` — `VerticalUnitsGeoKey` (4099)
    present and not 9001. Elevation feet with horizontal metres is a real USGS
    product and would produce a mesh wrong by a factor of 3.28 in one axis
    only. Absent means metres, and `RasterMeta` records that it was assumed.

**NoData**

15. `refuses_unparseable_nodata` and `refuses_nodata_not_representable` —
    section 6.
16. `refuses_contradictory_nodata_override` — section 6.

### What is explicitly *not* refused

- An unknown or unhandled GeoKey. It is ignored. The legacy raised `ValueError`
  on any key id outside its own incomplete enum, so a valid GeoTIFF carrying
  `GeogLinearUnitsGeoKey` could not be opened (prior art §5.10).
- `GeographicTypeGeoKey` being present. Ruling 7, and this is not hypothetical:
  the fixture carries `GeographicTypeGeoKey` 4258 (EUREF89, geographic) **and**
  `ProjectedCSTypeGeoKey` 25833 (ETRS89 / UTM 33N, projected). A refusal keyed
  on the presence of a geographic key refuses the reference DEM. Measured:
  `pyproj.CRS.from_epsg(25833).is_projected` is `True`;
  `pyproj.CRS.from_epsg(4258).is_projected` is `False`. Test the constructed
  CRS, never the key list.
- `GDAL_METADATA` (42112), citation keys, datum and ellipsoid keys. Ignored.
  The EPSG code carries all of it and `pyproj` resolves it.
- No refusal is a bare `assert`. Prior art §5.9: the legacy validated file
  content with five of them, all stripped by `python -O`. Every refusal here
  raises a real exception.

### The exception type

One exception, `GeoTiffError(ValueError)`, in `io/models.py`. Every refusal
raises it with a message that names the tag or GeoKey by number and by name,
and the file's value. Subclasses are not worth it: no caller will branch on
which refusal fired, and the one that might — the missing codec — is
distinguishable by its message naming the extra to install.

---

## 6. Ruling 9: NoData

This is the largest genuinely new decision in the increment. There is no prior
art: `GDAL_NODATA` (42113) appears nowhere in the legacy and a void in a DEM
became an elevation (prior art §1.3). `sample.hpp` already implements the
semantics — any of the four bilinear corners NoData returns `nullopt` — and
`raster.hpp`'s `is_nodata` also catches NaN unconditionally via `v != v`.
Discovery is all that is missing.

**Tag present.** Tag 42113 is ASCII. The fixture's is the string `'-32767'`.
Parse with `float()`. Then:

- text that parses to NaN (`'nan'`) yields `nodata = None`. `is_nodata` already
  catches NaN, and a NaN sentinel compared with `==` would match nothing, so
  passing it would be worse than passing nothing.
- otherwise convert to the **promoted array dtype** and require the round trip
  to be exact. `sample.hpp`'s comparison is `v == *nodata_`, so a sentinel that
  does not land on a representable value matches no cell and silently disables
  NoData entirely. `refuses_nodata_not_representable`.
- text that parses as neither: `refuses_unparseable_nodata`.

**Tag absent, which for most DEM products it is. Ruling: `nodata = None`. No
value is guessed, and the array is not scanned.**

The tempting alternative is to sniff for -9999 or -32767, or to take the array
minimum. Reject it. A wrong sentinel deletes real terrain, and it does so
invisibly and in the opposite direction from the failure it was meant to
prevent: a void that becomes an elevation makes a visible spike, while a valid
elevation that becomes a void makes a hole the mesh quietly interpolates
across. Between two silent failures, prefer the loud one. -9999 is also a
legitimate value in a bathymetric product.

Float products that mark voids with NaN are covered for free by `is_nodata`,
with no tag and no configuration. That covers a large share of the cases the
sniffing was meant to catch.

**The escape hatch.** `decode_dem(source, nodata=...)` takes an optional
explicit sentinel. It is the caller asserting a fact about their data, which is
different in kind from the reader guessing one. If the tag is also present and
the two differ, `refuses_contradictory_nodata_override`: one of the two is
wrong and the reader cannot tell which.

`RasterMeta` records whether the sentinel came from the tag, from the caller,
or is absent. An absent sentinel is a fact about the tile, and a later stage
that cares can ask.

---

## 7. Types

All Pydantic V2, all frozen. `io/models.py`.

```
RasterMeta
    x_min, y_max, delta_x, delta_y : float        # node grid, metres
    cols, rows                     : int          # from array.shape only
    epsg                           : int          # projected, metre
    nodata                         : float | None
    nodata_source                  : "tag" | "caller" | "absent"
    pixel_is_area                  : bool         # what the file declared
    vertical_unit_assumed          : bool         # section 5, refusal 14

DemTile
    meta  : RasterMeta
    array : numpy 2-D, C-contiguous, float32 or float64
```

`RasterMeta` is what crosses into increment 12's adapter, and it is
deliberately a *superset* of the boundary contract: `epsg`, `nodata_source`,
`pixel_is_area` and `vertical_unit_assumed` exist for diagnostics and provenance
and **must not** be forwarded to `_core`. The boundary contract still stands: one
array, four keyword-named affine scalars, one optional sentinel, nothing else.
`epsg` living on `RasterMeta` is not a CRS crossing the boundary; `raster.py`
not passing it is what makes that true, and increment 12 owns that test.

**Dimensions come from `array.shape`.** `cols` and `rows` are derived in the
model validator, never accepted as input. `project_structure.md` already rules
this and prior art §5.7 is the incident: the legacy built a shape out of
`numpy.float64` values and its guarding assertion compared `(12.0, 16.0)` to
`(12, 16)` and passed.

**Promotion table**, from `project_structure.md`:

| file dtype | array dtype |
|---|---|
| int8, uint8, int16, uint16 | float32 |
| int32, uint32 | float64 |
| float32 | float32 |
| float64 | float64 |
| anything else | `refuses_unsupported_dtype` |

`int32` into `float32` would silently quantise elevations; `int64` and `uint64`
do not fit float64's mantissa either, so they are refused rather than promoted.
Complex and boolean are refused.

The array is made C-contiguous and set read-only at construction. Read-only
now, in Python, rather than in increment 12's adapter, because
`project_structure.md`'s concurrency rule — every refinement thread sampling
the same raster — is a property of the buffer, and a buffer that was ever
writeable is one somebody can hold a writeable handle to.

---

## 8. Ruling 5 and 13: the decoder, and a false comment in `pyproject.toml`

**Use `tifffile.TiffFile.geotiff_metadata`.** It resolves the GeoKey directory,
including the 34736 double and 34737 ASCII indirections that prior art §2 warns
a `tag_v2`-derived reader would miss. Measured on the fixture, it returns
`GeogSemiMajorAxisGeoKey` (from 34736) and `GTCitationGeoKey` (from 34737)
alongside the inline shorts. `GeoKeysInterpreter` is 145 lines of reflection,
regex and proj4 assembly replaced by an attribute lookup, and prior art §6.4
already rules it must not be carried across.

The five tags prior art enumerates are still the five tags that matter: 33922
and 33550 read directly from `page.tags`, and 34735/34736/34737 read through
`geotiff_metadata`. 34264 is read only to refuse.

**`pyproject.toml`'s codec comment is false as installed.** It claims
"tifffile covers uncompressed, Deflate, PackBits and LZW with no extra
dependency". Measured, tifffile 2026.9.15, this tree's `.venv`, no
`imagecodecs`:

```
$ .venv/bin/python -c "
import tifffile, numpy as np, io
a = np.arange(16, dtype=np.float32).reshape(4, 4)
for comp in [None, 'deflate', 'packbits', 'lzw']:
    b = io.BytesIO()
    try:
        tifffile.imwrite(b, a, compression=comp)
        b.seek(0)
        print(comp, 'ok', tifffile.imread(b).dtype)
    except Exception as e:
        print(comp, 'FAIL', type(e).__name__, str(e)[:80])"
None ok float32
deflate ok float32
packbits FAIL KeyError "<COMPRESSION.PACKBITS: 32773> requires the 'imagecodecs' package"
lzw FAIL KeyError "<COMPRESSION.LZW: 5> requires the 'imagecodecs' package"
```

And the fixture `tests/fixtures/dem_archive/7908_3_10m_z33.tif` is LZW, so the
repository's only DEM cannot be read without the `codecs` extra. Two
consequences, both this increment's PR to land, per principle C3:

- Correct the comment: without `imagecodecs`, tifffile covers uncompressed and
  Deflate. PackBits and LZW need the extra, along with JPEG, ZSTD, LERC and
  WebP.
- `@tester` builds micro-TIFF fixtures with `compression=None` or `'deflate'`
  so the suite runs with no optional dependency, and marks any test that reads
  the real Kartverket fixture `skipif` on `imagecodecs` being importable.
  `refuses_missing_codec` is then testable in both directions: a PackBits
  micro-TIFF is refused with a useful message when the extra is absent, and
  read when it is present.

**Whole-page read.** `page.asarray()`. Windowed decode needs `page.aszarr()`
and therefore `zarr`, which is not a dependency and is not being added here.
The 5051² float32 fixture is 102 MB, which is acceptable for one tile, and
`project_structure.md`'s two meanings of "window" stay separate: the I/O window
is deferred with its `zarr` question, the index window is increment 12's.

---

## 9. Ruling 8: `always_xy=True`, enforced rather than asserted

Prior art §5.1 measured it: the same transform without `always_xy` puts a point
6 000 km away, silently, and the legacy was correct only because every CRS it
built was spelled `+init=`, which forces lon/lat order and is deprecated in
PROJ 6+. Re-measured here:

```
$ .venv/bin/python -c "
import pyproj
a = pyproj.Transformer.from_crs('EPSG:4326', 'EPSG:32633')
b = pyproj.Transformer.from_crs('EPSG:4326', 'EPSG:32633', always_xy=True)
print('no always_xy', a.transform(8.5, 60.9))
print('always_xy   ', b.transform(8.5, 60.9))"
no always_xy (6165244.414816713, 1344368.4624163064)
always_xy    (147736.9380934023, 6769135.917888515)
```

This increment performs no transform at all — one tile, no polygon, no
reprojection. That is precisely why the rule lands **now**, before the first
`Transformer` exists, rather than after somebody writes one.

**The rule.** Every `pyproj.Transformer.from_crs` in `src_python/` passes
`always_xy=True`. No CRS anywhere in `src_python/` is spelled `+init=`.

**The enforcement**, and it is a test rather than prose: a test greps
`src_python/` for `from_crs(` and for `+init=`, and asserts that every
`from_crs` call site passes `always_xy=True` and that `+init=` does not occur.
It passes vacuously today, so per principle A3 `@tester` must show it failing
against a planted violation, in the same commit, and record that it did.
`legacy/` is out of scope; it is frozen and every one of its CRS strings is
`+init=`.

A test, not a gate in `tools/`, because `tools/*.py` has no suite behind it and
this rule has exactly one consumer directory.

When increment 12 or the mosaic increment needs transforms, they go through one
helper so the grep has one site to find.

---

## 10. Ruling 10: the mosaic, and `GeoPolygon`

**The mosaic walk is a separate increment, and it is not this one.** Both
mosaic defects live there — a mosaic takes its CRS from the first tile in
`Path.glob` order and never checks the rest (prior art §5.8), and the stop
threshold `1e-10` is an absolute area in the squared units of whatever CRS the
domain happens to carry (§4.12). Neither is a decode defect. `ROADMAP.md` gap 1
needs one tile, and the refusal list above is already the size of an increment.

That increment's first refusal is `refuses_mixed_crs_mosaic`, and its
threshold is relative to the domain area rather than absolute.

**`GeoPolygon` is worth porting**, as prior art §8.4 recommends: polygon plus
CRS, with `transform`, `intersects`, `intersection`, `difference`, `buffer`,
minus `to_cpp` and minus the CGAL `SimplePolygon`/`Polygon` conversion. Keep
the `not touches` rule — two abutting DEM tiles share an edge, and without it
every neighbour is pulled in to contribute a zero-area sliver. That is real
knowledge about tiled archives, not an implementation detail.

**But it lands with the mosaic**, in a new `tin_engine/geo.py`, not here.
It is geometry plus CRS, not I/O, so it does not belong in `io/`; and it has no
caller in a single-tile decode. A module written one increment before anything
uses it is a module whose interface is guessed.

**`crop_image_to_polygon` is not ported**, per prior art §8.6. Its arithmetic
is sound and its quirks are all Pillow's: the `+1` exists because
`Image.crop` excludes the last index, and `tifffile` does not read by box. The
surviving content is the index arithmetic — rows run downward from `y_max`,
columns rightward from `x_min` — and `geometry.hpp`'s `clamped_cell_of` already
holds it. Its one genuine defect goes to increment 12 as a named refusal: a
polygon entirely west of the raster clipped to a one-column strip of the west
edge and claimed it was the requested region (prior art §5.6).

---

## 11. Prior art: what is carried across

Required by `docs/increments/README.md` step 1. The full analysis is
`docs/increments/11-raster-ingestion-prior-art.md`; this is the disposition.

```
$ git grep -ni 'gdal\|ogr\|fiona\|rasterio\|osgeo' -- legacy
legacy/bindings.cpp
legacy/rasputin/gml_repository.py
legacy/rasputin/reader.py
legacy/rasputin/solar_position.h
legacy/rasputin/triangulate_dem.h
```

Every hit is a substring false positive — `geographic`, `geographical`,
`geographic_latitude`, `bg::cs::geographic`, and four `ogr` XML-namespace
prefixes in GML produced by `ogr2ogr`. **No prohibited dependency is used**, so
nothing here has to be replaced for `CLAUDE.md` §2 compliance. The prior-art
report ran this command and lists the same five files.

**Carried across as knowledge, not as code:** north-up encoded as a positive
`delta_y` subtracted with the row index; grid registration spanning `n - 1`
spacings; the tie point unpack order `(col, row, _, x, y, _)`, whose legacy
variable names are the opposite of its own convention everywhere else; Pillow's
`size` being `(width, height)` against numpy's `(rows, cols)`; the `not touches`
rule for abutting tiles; the `160zz`/`161zz` UTM zone encoding, documented and
not implemented because the EPSG path covers it.

**Not carried across:** `GeoKeysInterpreter` and its proj4 assembly; the
`[20000, 32760]` EPSG window; the free-text ellipsoid regex; `+init=`; Pillow;
`crop_image_to_polygon`; clamping lookups; the `raster_data_float` /
`raster_data_double` pair; `make_mesh`'s proj4 parameter.

---

## 12. Testing notes for `@tester`

- **No suite here is invariant-critical, so no mutation round.** The failure
  mode of this module is a wrong refusal or a missing one, and each is one
  named test against one crafted micro-TIFF. `docs/increments/README.md`'s cost
  argument puts mutation rounds where topology decisions are; there are none
  here.
- **Fixtures are `io.BytesIO` micro-TIFFs built with `tifffile.imwrite`**, one
  per refusal, 4×4 or smaller, uncompressed or Deflate. Section 8.
- **Two tests against the real fixture**, both `skipif` on `imagecodecs`: a
  happy-path decode, and an assertion of the half-cell shift —
  `x_min == 799750.0` and `y_max == 7950250.0` from a tie point of
  `(799745.0, 7950255.0)` at a 10 m scale. Those two numbers are the only
  check that ruling 4 is implemented and not merely documented.
- **The `always_xy` test must be shown failing** against a planted violation.
  Section 9, principle A3.
- **The promotion table is a parametrised test**, one case per row, asserting
  the resulting `array.dtype` — not that it is "a float".

---

## 13. LOC estimate

**An estimate, not a measurement** (principle B4). Non-comment production
lines, per `CLAUDE.md` §2.

| file | estimate |
|---|---|
| `src_python/tin_engine/io/models.py` | 80 |
| `src_python/tin_engine/io/geotiff.py` | 230 |
| **total** | **310** |

`geotiff.py` is mostly the refusal list: eighteen named refusals at roughly
three lines each, plus tag extraction, the CRS resolution, the NoData path and
the promotion table. Comfortably inside the ceiling, which is deliberate —
increment 12 carries the C++ concept change and needs the room.

Tests are excluded from the ceiling and will be larger than the production
code; one micro-TIFF per refusal is the point of the design.
