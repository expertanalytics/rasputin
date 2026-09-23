# Increment 11, prior art: legacy raster ingestion

**Status.** Analysis, not design. Written by `@migration-expert` before
`@architect` writes `11-*.md` proper and before `@tester` is spawned, per
`docs/increments/README.md` step 1. No design, no code, no tests here.

Read at master `fd945cb`. Every line number below is a line in `legacy/`, which
is frozen and does not move.

Code blocks are verbatim from the cited line, with two exceptions marked where
they occur: `...` is an elision, and one over-long line in
`legacy/rasputin/triangulate_dem.h` is wrapped to fit. The trailing
`# legacy/...:NNN` comment in each block is the citation, not source text.
Verified by matching each quoted line against its cited file.

---

## 1. Findings first

1. **The legacy never used GDAL, OGR, Fiona or Rasterio.** It read GeoTIFF with
   Pillow's `tag_v2` and decoded the GeoKey directory by hand. The only `ogr`
   string in the tree is an XML namespace prefix inside CORINE GML files.
   Nothing in the raster path has to be replaced for §2 compliance.
2. **The legacy read exactly two kinds of georeferencing, and both give a
   projected CRS.** A file carrying `ProjectedCSTypeGeoKey` works only if that
   code lies in `[20000, 32760]`; a file carrying `ProjectionGeoKey` instead
   works through a different path that has no such window. Measured by
   executing the legacy `GeoKeysInterpreter`:

   ```
   ProjectionGeoKey 16033, no ProjectedCSTypeGeoKey
     -> '+proj=utm +zone=33 +units=m +ellps=WGS84 +no_defs'   pyproj: accepted
   ProjectedCSTypeGeoKey 32633
     -> '+init=EPSG:32633 +no_defs'                            pyproj: accepted
   ProjectedCSTypeGeoKey 3857  (below 20000)
     -> '+units=m +no_defs'                                    pyproj: CRSError
   ```

   An earlier revision of this section said the legacy could read *only*
   `[20000, 32760]`, which is the window on one path presented as the whole
   story. Increment 11 designs its refusals against this list, so the error
   would have become a test. Two bugs narrow the `ProjectedCSTypeGeoKey` path,
   both reproduced below. The practical consequence is that the "CRS never crosses
   into C++" question has an evidence answer: the legacy only ever handled
   projected CRS, by accident rather than by design.
3. **There is no NoData handling anywhere in the legacy DEM path.** Not the tag,
   not a sentinel, not a mask. `GDAL_NODATA` (42113) appears nowhere. A void in
   a DEM became an elevation.
4. **The C++ `raster/` module already absorbs most of the legacy's geometric
   intent**, including three defects it names explicitly. What is left to port is
   the Python decode, and it is mostly a list of things the legacy did *not* do.
5. **The single most dangerous thing to carry across is the axis order.** The
   legacy is correct only because every CRS it constructs is spelled
   `+init=EPSG:...`, and it never passes `always_xy=True`. A rewrite that spells
   the same CRS `EPSG:4326` and keeps the rest of the code unchanged puts every
   coordinate 6 000 km away. Measured in §5.1.
6. **A multi-tile mosaic takes its CRS from the first tile and never checks the
   others.** `legacy/rasputin/mesh.py:34`.

---

## 2. The greps, and what they returned

Run at `fd945cb`.

```
$ git grep -ni 'nodata\|no_data\|42113' -- legacy
legacy/rasputin/globcov_repository.py:38:    no_data = 230
legacy/rasputin/globcov_repository.py:68:            LandCoverType.no_data: "No data (burnt areas, clouds,…)"}
legacy/rasputin/globcov_repository.py:103:            LandCoverType.no_data: (0, 0, 0)
```

All three are the GlobCover land-cover class value 230. Nothing in the DEM path.

```
$ git grep -n 'always_xy' -- legacy
```

No output. Zero occurrences in the whole legacy tree.

```
$ git grep -n '+init=' -- legacy
legacy/rasputin/application.py:90:                                                    crs=pyproj.CRS.from_string("+init=EPSG:4326"))
legacy/rasputin/application.py:96:                                  crs=pyproj.CRS.from_string("+init=EPSG:4326"))
legacy/rasputin/application.py:101:    target_crs = pyproj.CRS.from_string(f"+init={res.target_coordinate_system}")
legacy/rasputin/globcov_repository.py:120:            self.data_crs = CRS.from_string("+init=EPSG:4326")
legacy/rasputin/gml_repository.py:159:        self.data_crs = CRS.from_string(f"+init={pstr}")
legacy/rasputin/reader.py:210:            proj4_str = f"+init=EPSG:{epsg}"
legacy/rasputin/web_visualize.py:54:    input_crs = pyproj.CRS.from_string("+init=EPSG:4326")
legacy/tests/test_gml_repository.py:19:    input_crs = pyproj.CRS.from_string("+init=EPSG:4326")
legacy/tests/test_land_cover_repository.py:29:    crs = pyproj.CRS.from_string("+init=EPSG:4326")
legacy/tests/test_land_cover_repository.py:44:    input_crs = pyproj.CRS.from_string("+init=EPSG:4326")
legacy/tests/test_mesh.py:38:                      coordinate_system=f"+init=epsg:{epsgid}",
legacy/tests/test_mesh.py:56:                      coordinate_system="+init=epsg:32633",
legacy/tests/test_mesh.py:76:                       coordinate_system="+init=epsg:32633",
legacy/tests/test_mesh.py:84:                       coordinate_system="+init=epsg:32633",
legacy/tests/test_raster_repository.py:17:    input_crs = pyproj.CRS.from_string("+init=EPSG:4326")
```

Every CRS in the legacy, production and test, is spelled with the deprecated
`+init=` prefix. That is not decoration; see §5.1.

```
$ git grep -ni 'gdal\|ogr\|fiona\|rasterio\|osgeo' -- legacy
```

Returns hits in five files — `legacy/bindings.cpp`,
`legacy/rasputin/gml_repository.py`, `legacy/rasputin/reader.py`,
`legacy/rasputin/solar_position.h`, `legacy/rasputin/triangulate_dem.h`. Every
one is a false positive on the substring:

- `geographic`/`geographical` identifiers in `solar_position.h`, `reader.py`,
  `triangulate_dem.h`, and `geographic_latitude`/`geographic_longitude`
  parameter names in `bindings.cpp` (380, 383, 387, 388);
- `bg::cs::geographic` in `legacy/rasputin/triangulate_dem.h:650`;
- four `ogr` XML-namespace lookups in `gml_repository.py` (145, 168, 170, 173),
  element names in GML produced by `ogr2ogr`, not an import.

**No prohibited dependency is used.** An earlier revision of this paragraph
listed four files where the command returns five, omitting `bindings.cpp`.
`docs/increments/README.md` asks for the command *with the file list it
returned* precisely so a reader can tell a run from a retelling, and a list that
is not what the command prints defeats that.

```
$ git grep -n 'Image.open' -- legacy
legacy/rasputin/geometry.py:157:        img = Image.open(rasputin_data_dir / self.material_spec["texture_file_name"])
legacy/rasputin/geometry.py:215:        with Image.open(filepath) as image:
legacy/rasputin/globcov_repository.py:117:        with Image.open(self.path) as image:
legacy/rasputin/globcov_repository.py:146:        with Image.open(self.path) as image:
legacy/rasputin/reader.py:409:    with Image.open(filepath) as image:
legacy/tests/test_read_raster_file.py:34:        tmp = Image.open(fp)
```

```
$ git grep -n 'tag_v2' -- legacy
legacy/rasputin/globcov_repository.py:121:            self.model_tie_point = image.tag_v2.get(GeoTiffTags.ModelTiePointTag.value)
legacy/rasputin/globcov_repository.py:122:            self.model_pixel_scale = image.tag_v2.get(GeoTiffTags.ModelPixelScaleTag.value)
legacy/rasputin/reader.py:96:    # Using .tag_v2 instead of legacy .tag
legacy/rasputin/reader.py:97:    image_tags = image.tag_v2
legacy/rasputin/reader.py:384:    j_tag, i_tag, _, x_tag, y_tag, _ = image.tag_v2.get(tiepoint_idx, (0, 0, 0, 0, 0, 0))
legacy/rasputin/reader.py:387:    delta_x, delta_y, _ = image.tag_v2.get(scale_idx, (1.0, 1.0, 0.0))
```

Five TIFF tags are read. Three through `tag_v2` — 33922, 33550, 34735 — and
two more through `image_tags[...]`, which is why the grep above does not show
them: 34736 for double-valued GeoKeys and 34737 for ASCII (§3.2). A reader
built from the `tag_v2` grep alone cannot decode either.

```
$ git grep -n GTRasterTypeGeoKey -- legacy
legacy/rasputin/reader.py:36:    GTRasterTypeGeoKey = 1025

$ git grep -n '34264\|ModelTransformation' -- legacy
```

The second returns nothing. Both confirm the gaps `project_structure.md` already
claims under "Obligations this puts on the reader".

---

## 3. What the legacy actually did, step by step

The whole DEM path is `legacy/rasputin/reader.py`, 475 lines, driven from
`legacy/rasputin/application.py:104-108`.

### 3.1 Open

`legacy/rasputin/reader.py:409` — `with Image.open(filepath) as image:`. Pillow. Lazy: the
array is only materialised at `legacy/rasputin/reader.py:420`, `np.array(image)`, after any
crop. So the legacy already read a window rather than the whole file, though
by accident of Pillow's laziness rather than by design.

### 3.2 Read the georeferencing

Three tags, and nothing else (`legacy/rasputin/reader.py:382-400`):

- 33922 `ModelTiePointTag`, unpacked as `(j, i, _, x, y, _)` — a raster point
  `(j, i)` and the model point `(x, y)` it maps to.
- 33550 `ModelPixelScaleTag`, unpacked as `(delta_x, delta_y, _)`.
- 34735 `GeoKeyDirectoryTag`, decoded by `extract_geo_keys` (`legacy/rasputin/reader.py:93-144`)
  into a name→value dict, with short values inline, doubles indirected through
  tag 34736 and ASCII through 34737.

The origin is then derived (`legacy/rasputin/reader.py:392-393`):

```python
    x_min = x_tag - delta_x * j_tag
    y_max = y_tag + delta_y * i_tag
```

and the far corner (`legacy/rasputin/reader.py:395-396`):

```python
    x_max = x_tag + delta_x * (n - 1 - j_tag)
    y_min = y_tag - delta_y * (m - 1 - i_tag)
```

### 3.3 Turn the GeoKeys into a CRS

`GeoKeysInterpreter` (`legacy/rasputin/reader.py:147-291`) reflects on its own method names:
for each GeoKey it looks up `self._<KeyName>`, calls it, and merges the returned
proj4 fragments, raising on conflicts (`legacy/rasputin/reader.py:186-191`). `to_proj4`
(`legacy/rasputin/reader.py:203-222`) short-circuits to `+init=EPSG:<code>` if any handler
produced an `EPSG` entry, otherwise concatenates the fragments, and always
appends `+no_defs`.

### 3.4 Choose the region

Two levels.

**Tile selection**, `RasterRepository.get_intersections`
(`legacy/rasputin/reader.py:434-456`): glob `*.tif` in a directory; for each, build a
`GeoPolygon` from its extent and CRS (`legacy/rasputin/geometry.py:212-220`); keep it if the
target polygon intersects; transform the target polygon *into that tile's CRS*
and crop; then subtract the tile from the remaining target and stop when the
remainder is small. Tiles are consumed in `Path.glob` order, which is
filesystem order — not sorted.

**Windowing inside a tile**, `crop_image_to_polygon` (`legacy/rasputin/reader.py:347-380`):
convert the polygon's bounding box to row/column indices, clamp to the image,
and call Pillow's `Image.crop`. The new extent is recomputed from the clamped
indices (`legacy/rasputin/reader.py:375-378`) rather than re-read.

Only the polygon's **bounding box** is used. The polygon itself never masks the
raster; it is carried separately into C++ and used there to clip the mesh.

### 3.5 Arrays

`legacy/rasputin/reader.py:420` — `image_array = np.array(image)`. Whatever dtype Pillow gives.
No promotion, no dtype check, no contiguity check, no NoData mask. The comment
above it ("Store the array to ensure the pointer to the data stays alive") is
the only nod to lifetime; the actual anchor is `keep_alive<1,2>` in
`legacy/bindings.cpp:156`.

### 3.6 NoData

Nothing. See §2.

### 3.7 Coordinate systems and reprojection

**The raster is never reprojected.** `legacy/rasputin/application.py:101-112` does the
opposite: it picks the CRS of the first DEM tile that intersects the domain
(`legacy/rasputin/reader.py:458-465`), transforms the *user's polygon* into that CRS, meshes in
that CRS, and then transforms the *finished mesh's vertices* into the target CRS
(`legacy/rasputin/application.py:116-123`). The proj4 string of the meshing CRS is pushed into
C++ and stored on the mesh (`legacy/bindings.cpp:188-193`).

This matters for §7.

---

## 4. Domain knowledge encoded in constants and conditionals

Each of these is a decision a rewrite gets silently wrong.

### 4.1 North-up, encoded as a sign

```python
    y_max = y_tag + delta_y * i_tag          # legacy/rasputin/reader.py:393
    y_min = y_tag - delta_y * (m - 1 - i_tag) # legacy/rasputin/reader.py:396
```

`delta_y` is stored **positive** and subtracted as the row index grows. Row 0 is
the north edge. The GeoTIFF `ModelPixelScaleTag` ScaleY is positive by
specification and the north-up-ness is implicit in it; the legacy bakes that
implication in and never checks it. Matched by
`include/terrain/raster/geometry.hpp`'s `y_min()` and `node()`.

### 4.2 Grid-registered: the extent spans `n - 1` spacings, not `n`

```python
    @property
    def x_max(self):
        return self.x_min + self.delta_x*(self.shape[1] - 1)   # legacy/rasputin/reader.py:308-309
```

Samples are points on a grid, not areal cells. This is the half-pixel decision.
`GTRasterTypeGeoKey` (1025), which is the tag that says whether the file is
pixel-is-point or pixel-is-area, is **defined at `legacy/rasputin/reader.py:36` and read
nowhere**. An area-registered file therefore lands half a cell north-west and one
cell too small in each axis, with no warning.

### 4.3 Tie point unpack order is `(j, i, _, x, y, _)` — column first

```python
    j_tag, i_tag, _, x_tag, y_tag, _ = image.tag_v2.get(tiepoint_idx, (0, 0, 0, 0, 0, 0))
                                                        # legacy/rasputin/reader.py:384
```

The GeoTIFF tie point is `(I, J, K, X, Y, Z)` where `I` is the **column** (raster
x) and `J` the row. The legacy names them `j_tag`, `i_tag` — the opposite of its
own `(i, j) = (row, col)` convention used everywhere else in the same file. The
values are used consistently with the spec (`j_tag` multiplies `delta_x`), so the
code is right and the names are wrong. A rewrite that trusts the names transposes
the origin.

### 4.4 Pillow's `size` is `(width, height)`

```python
    n, m = image.size          # legacy/rasputin/reader.py:389
    ...
    return ImageExtents(shape=(m, n), ...)   # legacy/rasputin/reader.py:398
```

`n` is columns, `m` is rows, and `shape` is `(rows, cols)` to match numpy. Any
replacement reader must reproduce this swap at its own boundary.

### 4.5 Defaults on missing georeferencing

```python
    j_tag, i_tag, _, x_tag, y_tag, _ = image.tag_v2.get(tiepoint_idx, (0, 0, 0, 0, 0, 0))
    ...
    delta_x, delta_y, _ = image.tag_v2.get(scale_idx, (1.0, 1.0, 0.0))
                                                        # legacy/rasputin/reader.py:384, 387
```

A TIFF with no georeferencing at all becomes a unit-spaced raster at the origin
and proceeds. This is the one place where "refuse" is unambiguously right and the
legacy chose "default".

### 4.6 The projected-EPSG window `[20000, 32760]`

```python
    @staticmethod
    def _ProjectedCSTypeGeoKey(value):
        if _isinteger(value) and 20000 <= value <= 32760:
            return {"EPSG": value}          # legacy/rasputin/reader.py:224-227
```

The upper bound excludes 32767, the GeoTIFF "user-defined" sentinel. The lower
bound has no such justification: it silently drops every projected EPSG code
below 20000. See §5.2.

### 4.7 The UTM zone encoding `160zz` / `161zz`

```python
            zone = value % 16000
            south = bool(zone // 100)
            if south:
                zone %= 100
            return dict(proj="utm", zone=zone, south=south)   # legacy/rasputin/reader.py:257-262
```

`ProjectionGeoKey` values 16001-16060 are northern UTM zones and 16101-16160
southern. Worth keeping as documented knowledge; not worth keeping as code, since
the EPSG path covers it whenever the file carries a `ProjectedCSTypeGeoKey`.

### 4.8 Linear units: metres only

```python
        unit_name = {9001: "m"}[value]       # legacy/rasputin/reader.py:267
```

EPSG unit code 9001 is metre. Anything else — 9002 international foot, 9003 US
survey foot — raises `KeyError` from a dict lookup, which is a crash rather than a
wrong answer. That is the better failure, but the message is useless.

### 4.9 Ellipsoid parsed out of a free-text citation with a regex

```python
        gcs_re = {"GRS":    ("GRS[ ,_]{0,1}(19|)\\d\\d", "\\d\\d$"),
                  "WGS":    ("WGS[ ,_]{0,1}(19|)\\d\\d", "\\d\\d$"),
                  "sphere": ("sphere", None)}    # legacy/rasputin/reader.py:275-277
```

`flags=2` at `legacy/rasputin/reader.py:280` is `re.IGNORECASE` spelled as a magic number. This
reads `GTCitationGeoKey`-style human prose such as `"GCS_WGS_1984"` and produces
`+ellps=WGS84`. It is a guess dressed as a parse, and it is the fallback the
broken paths in §5.2 fall into.

### 4.10 The crop window's `+1`

```python
    j_max = np.clip(np.ceil((x_max_p - x_min)/delta_x) + 1, 1, n)
    i_max = np.clip(np.ceil((y_max - y_min_p)/delta_y) + 1, 1, m)   # legacy/rasputin/reader.py:364-365
```

The comment three lines below says why: Pillow's `Image.crop` excludes the last
index (`legacy/rasputin/reader.py:367-370`). The `+1` is a Pillow artefact, not geometry. A
`tifffile`-based reader must not copy it blindly — but it must keep *some*
margin, because bilinear sampling at the window edge needs the next node outside
the requested box.

### 4.11 Row index runs from `y_max` downward

```python
    i_min = np.clip(np.floor((y_max - y_max_p)/delta_y), 0, m - 1)   # legacy/rasputin/reader.py:363
```

`y_max - y` for rows, `x - x_min` for columns. The asymmetry is the whole
north-up convention in one line. Reproduced in
`include/terrain/raster/geometry.hpp`'s `clamped_cell_of`.

### 4.12 The mosaic stop threshold, `1e-10`

```python
                if target_polygon.polygon.area < 1e-10:
                    break                     # legacy/rasputin/reader.py:454
```

An absolute area, in the squared units of whatever CRS the domain currently
carries. In UTM metres that is 1e-10 m², which is fine. In degrees it is about
1.2 km² at the equator. Another constant that only works because the legacy
always ended up in a projected CRS.

### 4.13 Touching tiles do not count as intersecting

```python
    def intersects(self, other: "GeoPolygon") -> bool:
        op = other.transform(target_crs=self.crs).polygon
        sp = self.polygon
        return sp.intersects(op) and not sp.touches(op)   # legacy/rasputin/geometry.py:248-251
```

Two DEM tiles that abut share an edge. Without the `not touches`, every
neighbouring tile would be pulled in and contribute a zero-area sliver. Keep the
rule; it is real domain knowledge about tiled DEM archives.

### 4.14 The boundary epsilon `sqrt(dx² + dy²) * 1e-10`

```cpp
        double eps = pow(pow(delta_x, 2) + pow(delta_y, 2), 0.5) * 1e-10;
        // legacy/rasputin/triangulate_dem.h:421
```

Already carried across, and already corrected, in
`include/terrain/raster/geometry.hpp::boundary_epsilon`, which records that the
legacy value is about 1.5 ulp at Norwegian UTM33 northings and therefore
indistinguishable from rounding noise.

### 4.15 Pillow's decompression-bomb limit is disabled

```python
Image.MAX_IMAGE_PIXELS = None     # legacy/rasputin/globcov_repository.py:12
```

Needed because the GlobCover global raster is larger than Pillow's default cap.
Not a DEM concern, but it records a real fact: land-cover rasters in this domain
are of an order where a naive full read is not an option.

---

## 5. Defects worth knowing about

### 5.1 No `always_xy`, and correctness resting on `+init=`

`always_xy` appears **zero** times in `legacy/` (§2). Every
`Transformer.from_crs` therefore honours the authority's declared axis order.
The code survives only because every CRS is constructed from a `+init=` string,
and `+init=` forces longitude-latitude order.

Measured, at `fd945cb`, with this tree's `.venv` (pyproj 3.8.0, PROJ 9.8.1):

```
$ .venv/bin/python -c "..."
pyproj 3.8.0 PROJ 9.8.1
init= axis: ['lon', 'lat']
plain axis: ['Lat', 'Lon']
no always_xy, EPSG:4326->32633, transform(8.5,60.9) = (6165244.414816713, 1344368.4624163064)
always_xy                        transform(8.5,60.9) = (147736.9380934023, 6769135.917888515)
```

The correct UTM33N easting/northing for that point is the second line. The first
is what `legacy/rasputin/geometry.py:237` would produce if `+init=EPSG:4326` at
`legacy/rasputin/application.py:90` were modernised to `EPSG:4326` and nothing else changed —
6 million metres of error, and no exception.

`+init=` is deprecated in PROJ 6+. The rewrite cannot keep it. So the rewrite
**must** pass `always_xy=True`, as `geospatial-data-formats/SKILL.md` §3 already
requires. This is the one place where modernising the CRS spelling and *not*
changing the transformer call is catastrophic, and where a test written against
re-derived intent would happily pin the wrong answer.

### 5.2 A geographic GeoTIFF could not be read at all — two bugs

**Bug one: the handler name does not match the GeoKey name.** The lookup is
`getattr(self, f"_{geokey_name}")` (`legacy/rasputin/reader.py:175`). The GeoKey is named
`GeographicTypeGeoKey` (`legacy/rasputin/reader.py:40`). The method is named:

```python
    @staticmethod
    def _GeoGraphicTypeGeo(value):        # legacy/rasputin/reader.py:229-230
```

`GeoGraphic` for `Geographic`, and no `Key` suffix. The handler is never called;
the key lands in `ignored_keys` and is logged at debug level.

**Bug two: `[20000, 32760]` drops valid projected codes** (§4.6).

Reproduced by executing the legacy `GeoKeysInterpreter` directly (the class was
extracted by `exec`-ing `legacy/rasputin/reader.py` up to `identify_projection`, so that the
Pillow and CGAL imports do not have to be satisfied):

```
handler present for GeographicTypeGeoKey: False
method actually defined                 : ['_GeoGraphicTypeGeo']
geographic tiff (EPSG:4326)  -> '+ellps=WGS84 +no_defs'
UTM33N tiff     (EPSG:32633) -> '+init=EPSG:32633 +no_defs'
EPSG:3857 tiff               -> '+ellps=WGS84 +no_defs'
pyproj.CRS.from_proj4('+ellps=WGS84 +no_defs')
  -> CRSError: Invalid projection: +ellps=WGS84 +no_defs +type=crs
```

The string has no `+proj=` at all, because the only handler that would have
supplied one was skipped. `GeoPolygon.from_raster_file` calls
`pyproj.CRS.from_proj4` on it (`legacy/rasputin/geometry.py:220`), so
`RasterRepository.get_intersections` dies on the first such file in the archive
directory — including files it was not asked about, because it builds a
`GeoPolygon` for *every* `*.tif` before testing intersection
(`legacy/rasputin/reader.py:442-445`).

This is the answer to §7: the legacy handled only projected CRS because
everything else crashed.

### 5.3 The transposed read, with its root cause

`project_structure.md` and `include/terrain/raster/geometry.hpp` both cite a
transposed row/column read. Here is where it lives.

```cpp
        indices.emplace_back(std::array<unsigned int, 2>{(unsigned int)((pt[0]-x0)/dx),
                                                         (unsigned int)((y1 - pt[1])/dy)});
        // legacy/rasputin/triangulate_dem.h:908
```

`coordinates_to_indices` emits `(col, row)`. Its consumer:

```cpp
        result.emplace_back(ptr[idx[0]*N + idx[1]]);  // TODO: Implement range check?
        // legacy/rasputin/triangulate_dem.h:921
```

`extract_buffer_values` treats `idx[0]` as the **row**. The two disagree. On a
non-square raster this reads the wrong cell or runs off the buffer; on a square
one it silently returns the transpose.

The legacy noticed the symptom and never found the cause:

```python
        # TODO: This does not work, perhaps figure out why?
        #all_land_types = np.asarray(td.extract_uint8_buffer_values(indices, image))
        all_land_types = np.asarray([image.getpixel(tuple(idx)) for idx in indices])
        # legacy/rasputin/globcov_repository.py:164-166
```

The fallback is correct, because Pillow's `getpixel` takes `(x, y)` = `(col,
row)`, which is what `coordinates_to_indices` produces. So the working path
agreed with the producer and the fast path did not. The named `CellIndex` in
`include/terrain/raster/geometry.hpp` already makes this unwriteable — this
section is here so the record of *why* exists somewhere.

Two further defects in the same two functions: `M` and `N` are parameters of
`coordinates_to_indices` and are never used, so there is no clamping; and the
cast to `unsigned int` turns any point west or north of the origin into a huge
positive index. The `TODO` on line 921 is the author saying so.

### 5.4 Out-of-bounds read in bilinear interpolation

```cpp
        int i = std::min<int>(std::max<int>(static_cast<int>((y_max-y) /delta_y), 0), num_points_y - 1);
        int j = std::min<int>(std::max<int>(static_cast<int>((x-x_min) /delta_x), 0), num_points_x - 1);
        // legacy/rasputin/triangulate_dem.h:369-370
```

clamps to the last index, and then `get_interpolated_value_at_point` reads
`data[(i + 1)*num_points_x + j + 1]` (`legacy/rasputin/triangulate_dem.h:392`) unconditionally.
Already recorded and fixed by `bilinear_cell_of` in the shipped
`include/terrain/raster/geometry.hpp`. Noted here only to confirm the shipped
comment against the source.

### 5.5 Clamping means a query outside the raster returns a plausible number

Same lines. A point anywhere on Earth gets an elevation. Already fixed by
`cell_of` returning `std::nullopt`.

### 5.6 A window is returned for a polygon that misses the raster entirely

```python
    j_min = np.clip(np.floor((x_min_p - x_min)/delta_x), 0, n - 1)
    ...
    j_max = np.clip(np.ceil((x_max_p - x_min)/delta_x) + 1, 1, n)
    # legacy/rasputin/reader.py:362-365
```

If the polygon lies entirely west of the raster, both quantities are negative
before clipping; afterwards `j_min = 0` and `j_max = 1`. The function returns a
one-column strip of the raster's west edge with an extent claiming to be the
requested region. It does not refuse. `get_intersections` normally shields this
by testing intersection first, but `crop_image_to_polygon` is public and
`legacy/tests/test_read_raster_file.py:83` calls it directly.

### 5.7 The window is computed from floats and passed to Pillow as floats

`j_min`, `i_min`, `j_max`, `i_max` are `numpy.float64` (`np.clip` of `np.floor`),
handed straight to `image.crop(box=...)` at `legacy/rasputin/reader.py:372`, and then used to
build `ImageExtents(shape=(i_max - i_min, j_max - j_min), ...)` at
`legacy/rasputin/reader.py:375`. The shape of a raster is a pair of floats. The assertion that
is supposed to catch a shape/array mismatch,

```python
        shapes = self.shape, self.array.shape
        assert shapes[0] == shapes[1], ...    # legacy/rasputin/reader.py:331-332
```

compares `(12.0, 16.0)` to `(12, 16)` and passes. `project_structure.md`'s rule
that dimensions are derived from `array.shape` and never passed alongside it
removes this whole class.

### 5.8 A mosaic takes its CRS from the first tile and never checks the rest

```python
            proj4_str = data[0].coordinate_system
            for raster in data:
                rasterdata_cpp.add_raster(raster.to_cpp())
            # legacy/rasputin/mesh.py:34-36
```

Tiles arrive from `Path.glob` in filesystem order (`legacy/rasputin/reader.py:440`). Two tiles in
different UTM zones are merged into one coordinate space with no complaint.

### 5.9 Bare `assert` for input validation

`legacy/rasputin/reader.py:110` (`assert KeyDirectoryVersion == 1`), `legacy/rasputin/reader.py:142`
(`assert len(geo_keys) == NumberOfKeys`), `legacy/rasputin/reader.py:320-321`, `legacy/rasputin/reader.py:332`,
`legacy/rasputin/reader.py:407`. All are stripped by `python -O`. All are validating *file
content*, which is exactly what must not be validated by `assert`.

### 5.10 The GeoKeys enum is incomplete, and an unknown key is fatal

```python
        key_name = GeoKeys(key_id).name     # legacy/rasputin/reader.py:117
```

`GeoKeys(...)` raises `ValueError` for any key id not in the enum. The enum at
`legacy/rasputin/reader.py:33-86` omits 2052 (`GeogLinearUnitsGeoKey`), among others. A
perfectly valid GeoTIFF carrying that key cannot be opened. The correct
behaviour is to keep the unknown key's numeric id and carry on.

### 5.11 `extract_geo_keys` can leave `key_value` unbound

`legacy/rasputin/reader.py:119-141` is three independent `if` statements over `KeyValueTags`, not
an `if/elif/else`, and there is no final `else`. The binding survives from the
previous loop iteration if none matches. Unreachable today because
`KeyValueTags(location)` raises first, but it is one enum extension away from
writing a stale value under a fresh key name.

---

## 6. What must not be carried across

1. **`make_mesh(raster_data, polygon, proj4_str)`** — `legacy/bindings.cpp:187-194`,
   reached from `legacy/rasputin/mesh.py:43`. A CRS string as a C++ parameter.
   Forbidden by `CLAUDE.md` §2 and `docs/PRINCIPLES.md` E3. Nothing in the new
   `raster/` module has a CRS field, and nothing should acquire one.
2. **Pillow as the GeoTIFF decoder.** Not prohibited, but `CLAUDE.md` §2's Core
   Stack names `tifffile`, and Pillow forced three of the defects above: the
   `+1` in §4.10, `MAX_IMAGE_PIXELS` in §4.15, and "cropping discards tiff tags"
   (`legacy/rasputin/reader.py:370`) which is why `read_raster_file` has to extract tags before
   cropping (`legacy/rasputin/reader.py:410-412`).
3. **`+init=EPSG:...`.** Deprecated, and removing it without adding
   `always_xy=True` is §5.1.
4. **`GeoKeysInterpreter`'s proj4 assembly.** Reflection on method names, regex
   over free text, `+no_defs` appended unconditionally. `pyproj.CRS` accepts an
   EPSG code directly; the GeoKey directory only has to yield that code.
5. **Clamping lookups.** Already replaced by `std::optional` in the shipped
   header; the Python side must not reintroduce a clamping window.
6. **CGAL `SimplePolygon`/`Polygon` conversion** in `legacy/rasputin/geometry.py:261-277`,
   including its note that "CGAL polygons have no vertex repeated and
   orientation matters". The new PSLG has its own chain rules.
7. **The `triangulate_dem` module name and its `raster_data_float` /
   `raster_data_double` pair.** Dtype promotion is now decided in Python at
   decode time (`project_structure.md`), so the binding does not need a
   per-dtype class exposed to the caller.

---

## 7. What the legacy did about projection — the evidence, not the argument

`project_structure.md` (the `raster` section, "Decided (@architect): GeoTIFF
decoding lives in Python") justifies the rule partly with *"CRS interpretation
means PROJ — GDAL's own dependency."* That justification is backwards: PROJ is a
standalone library, GDAL depends on PROJ, and `CLAUDE.md` §2 lists PyProj as
approved. This report does not rewrite that paragraph. It supplies the evidence
whoever does will need.

**What the legacy actually did:**

- It never reprojected a raster. Not once. There is no resampling code in the
  legacy tree.
- It meshed in **the DEM's own CRS**, chosen by whichever tile
  `RasterRepository.coordinate_system` found first (`legacy/rasputin/reader.py:458-465`), and
  reprojected the finished mesh's vertices afterwards
  (`legacy/rasputin/application.py:116-123`).
- It moved vector data — the domain polygon, the land-cover query points —
  across CRS boundaries, never grids. `legacy/rasputin/geometry.py:235-241` and
  `legacy/rasputin/globcov_repository.py:141-143`.
- It never rejected a geographic CRS, and there is no check anywhere that the
  meshing CRS is projected. The protection was accidental: §5.2 shows that a
  geographic GeoTIFF raised `CRSError` during CRS construction and never reached
  the mesher. The defence was a typo.

**What this supports and what it does not.**

It supports the *hazard* half of the paragraph. A degrees-based raster
interpolates happily: nothing in `legacy/rasputin/triangulate_dem.h:376-395`, and nothing in the
shipped `sample.hpp`, is aware of units. `delta_x = delta_y = 0.0002777°` yields
cells about 31 m by 15 m at 60°N, a 2:1 anisotropy that no code on either side of
the boundary can see. The refinement error metric, which compares an interpolated
height against a mesh height, would be measuring a quantity whose horizontal
units are degrees and whose vertical units are metres. That is a real and
unguarded silent failure, and §4.12's `1e-10` area threshold is a second
instance of the same hazard.

It does **not** support "CRS interpretation means PROJ, and PROJ is GDAL's
dependency". The legacy used `pyproj` directly for every transformation and
never linked GDAL (§2).

It also does not support the rule as *stated*. "CRS never crosses into C++" and
"reject a geographic CRS" are two different rules. The first is about where a
string lives; the second is about what the numbers mean. The evidence above
argues strongly for the second. The first is defensible on its own terms — the
C++ core has no use for a CRS, and a field nothing reads is a field that rots —
but that is an argument about dependency surface and unused state, not about
PROJ. Whoever rewrites the paragraph should separate the two.

---

## 8. What this changes about the next increment

1. **The Python reader's job is mostly refusal.** Of the eleven defects in §5,
   seven are "returned something instead of raising". The increment's design
   should enumerate what the reader refuses — area-registered files, rotated
   transforms (`ModelTransformationTag`, absent from the legacy entirely),
   missing georeferencing, a geographic CRS, a non-metre linear unit, a mosaic
   whose tiles disagree on CRS — and treat each as a named test.
2. **NoData has no prior art.** §1.3. `sample.hpp` already implements NoData
   *semantics*; discovery is greenfield. The design must say what the reader
   does when `GDAL_NODATA` (42113) is absent, which for most DEM products it is.
   There is no legacy answer to consult and nothing to preserve. This is the
   largest genuinely new decision in the increment.
3. **The CRS rule needs `always_xy=True` stated as a hard requirement in the
   increment file**, not left to the skill. §5.1 is the incident.
4. **`GeoPolygon` is worth porting, roughly as-is, minus `to_cpp`.**
   `legacy/rasputin/geometry.py:205-260` — polygon plus CRS, with `transform`, `intersects`
   (keeping the `not touches` rule of §4.13), `intersection`, `difference`,
   `buffer`. It is small, it is the right abstraction, and the mosaic walk in
   `legacy/rasputin/reader.py:434-456` is built entirely out of it. Pydantic V2 plus
   `shapely` covers it.
5. **The mosaic walk is a separate concern from the single-file decode**, and
   probably a separate increment. §5.8 and §4.12 are both mosaic bugs, and
   nothing in `ROADMAP.md`'s gap 1 requires more than one tile.
6. **`crop_image_to_polygon` should not be ported.** Its arithmetic is sound
   (§4.10, §4.11) but every one of its quirks is a Pillow quirk, and `tifffile`
   reads by page and tile rather than by box. Port the *index arithmetic*, drop
   the function.
7. **Two unrelated meanings of "window"** already flagged in
   `project_structure.md`. §3.4 confirms the legacy had only the I/O meaning,
   and confirms it used the polygon's bounding box and never the polygon.

---

## 9. Files read

- `legacy/rasputin/reader.py` (all 475 lines)
- `legacy/rasputin/geo_tiff_reader.py`, `legacy/rasputin/land_cover_repository.py`
- `legacy/rasputin/globcov_repository.py`, `legacy/rasputin/geometry.py`
- `legacy/rasputin/application.py`, `legacy/rasputin/mesh.py`
- `legacy/rasputin/triangulate_dem.h` (`RasterData`, `interpolate_boundary_points`,
  `coordinates_to_indices`, `extract_buffer_values`)
- `legacy/bindings.cpp` (raster and mesh bindings)
- `legacy/tests/test_read_raster_file.py`, `legacy/tests/test_raster_repository.py`
- `include/terrain/raster/{geometry,raster,sample}.hpp`
