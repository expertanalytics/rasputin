# GeoTIFF: a projected CRS given by parameters, read as the EPSG code it is

Status: **code review round 2 approved at `35c44c8`; +159 net. Next: push
on Ola's yes, after audit PR B (branch `worktree-audit-crs`) merges.** Design
approved at round 2 (`900701c`), red steps `0f6c271` and `5270a64`, green
steps `c8984e8` and `35c44c8`.
Written by `@architect` on branch
`worktree-geotiff-param-crs`, from audit PR B's green head `65cd528`
(`same_crs`, `docs/increments/python-audit.md` section 9). It lands after PR B.
Not refine or mesh code, so no `@perf` run.

Ola's ruling, verbatim: "D9: a" (2026-10-05). D9 was the main session's
question, from `docs/increments/python-audit.md` section 9, question 1: should
the GeoTIFF reader accept a projected CRS given by parameters (GeoKeys, no EPSG
code) when it matches an EPSG code, so that Ola's Austrian openDEM file reads?
Option a was "yes, as a later small PR after PR B". It overturns increment 11's
ruling 6 ("only `ProjectedCSTypeGeoKey` can yield an accepted CRS, resolved
through `pyproj.CRS.from_epsg`").

## 1. What changes, in one paragraph

Today a file whose `ProjectedCSTypeGeoKey` (3072) is 32767 ("user-defined":
the CRS is spelt out by its parameters rather than named by a code) is refused
at `src_python/tin_engine/io/geotiff.py@44fa7f5:285-288`. After this PR the
reader builds that CRS from the GeoKeys, on the datum the file names by EPSG
code in `GeographicTypeGeoKey` (2048), looks for EPSG projected CRSs on that
same datum that PROJ calls equivalent to it, and, when it finds one (or
several that are one CRS under several numbers), reads the file **as that
EPSG code**: `RasterMeta.epsg` is the code and `meta.crs` is `EPSG:n`, exactly
as if the file had carried the code. Nothing downstream changes. Everything
else stays refused, with a sentence that says why.

The file that motivates it: `../rasputin_data/austria_dgm10/dhm_at_lamb_10m_2018.tif`
(openDEM Austria, 10 m, float32, 58 061 x 31 793 cells, 256 x 256 tiles,
Deflate). Its GeoKeys, read with tifffile 2026.9.20 (probe in section 3):

| GeoKey | value |
|---|---|
| GTModelTypeGeoKey (1024) | 1, projected |
| GTRasterTypeGeoKey (1025) | 1, pixel is area |
| GTCitationGeoKey (1026) | `MGI_Austria_Lambert` |
| GeographicTypeGeoKey (2048) | 4312, MGI |
| GeogCitationGeoKey (2049) | `MGI` |
| GeogAngularUnitsGeoKey (2054) | 9102, degree |
| GeogSemiMajorAxisGeoKey (2057) | 6377397.155 |
| GeogInvFlatteningGeoKey (2059) | 299.1528128000033 |
| ProjectedCSTypeGeoKey (3072) | 32767, user-defined |
| ProjectionGeoKey (3074) | 32767, user-defined |
| ProjCoordTransGeoKey (3075) | 8, Lambert conic conformal, two standard parallels |
| ProjLinearUnitsGeoKey (3076) | 9001, metre |
| ProjStdParallel1GeoKey (3078) | 46.0 |
| ProjStdParallel2GeoKey (3079) | 49.0 |
| ProjFalseOriginLongGeoKey (3084) | 13.33333333300013 |
| ProjFalseOriginLatGeoKey (3085) | 47.5 |
| ProjFalseOriginEastingGeoKey (3086) | 400000.0 |
| ProjFalseOriginNorthingGeoKey (3087) | 400000.0 |

Tie point (0, 0) at (108 875, 586 555), pixel scale 10 x 10, GDAL_NODATA
`-3.4028234663852886e+38`. EPSG:31287 (MGI / Austria Lambert) has the same
parameters, with the longitude of origin 13.3333333333333: the file's differs
by 3.3e-10 degrees, about 0.03 mm on the ground.

## 2. Prior art: legacy and literature

*Literature.* The method is the one the GeoTIFF standard and its two reference
readers already use; nothing new is claimed.

- **OGC GeoTIFF 1.1** (OGC 19-008r4, clause 7,
  <https://docs.ogc.org/is/19-008r4/19-008r4.html>): 3072 = 32767 means a
  user-defined projected CRS, which then needs `GeodeticCRSGeoKey` (2048, called
  `GeographicTypeGeoKey` in 1.0) and `ProjectionGeoKey` (3074) populated; the
  ellipsoid keys (2056-2059) are for a *user-defined* ellipsoid; angular map
  projection parameters are in the unit of `GeogAngularUnitsGeoKey` (2054); and
  a user-defined CRS always has the axis order east, north. GeoTIFF 1.0 tied
  no parameter list to a method, and 1.1 names no EPSG method or parameter
  either ("a future version of GeoTIFF might do this").
- **libgeotiff**, `GTIFGetDefn` in `libgeotiff/geo_normalize.c`
  (<https://github.com/OSGeo/libgeotiff>): reads the method from 3075 and each
  parameter from its GeoKey, **falling back** between alternatives (for Lambert
  2SP, the longitude of origin from `ProjNatOriginLongGeoKey`, then
  `ProjFalseOriginLongGeoKey`, then `ProjCenterLongGeoKey`) and setting an
  absent parameter to 0. It also lets `GeogSemiMajorAxisGeoKey` and
  `GeogInvFlatteningGeoKey` override the ellipsoid of the coded geographic CRS.
- **GDAL**, `GTIFGetOGISDefnAsOSR` in `frmts/gtiff/gt_wkt_srs.cpp`
  (<https://github.com/OSGeo/gdal>): builds an `OGRSpatialReference` from
  libgeotiff's normalised definition. It does **not** look for an EPSG code
  for a user-defined projected CRS when it reads the file; its only
  `FindBestMatch(100)` there is for a compound CRS whose two parts both carry
  codes. Identification is left to the caller. It **keeps
  `GeogTOWGS84GeoKey` (2062)** for 3072 = 32767: its `OSR_STRIP_TOWGS84`
  strip applies only when the CRS was imported from the 3072 code
  (`bGotFromEPSG`, which needs 3072 not 32767), and otherwise it calls
  `SetTOWGS84`, so the result is a bound CRS
  (`gt_wkt_srs.cpp` at GDAL commit `d6cd774`, lines 1268-1289 and 1327-1371). When writing,
  it sets 2062 when the CRS it is given carries a TOWGS84 other than the one
  GDAL guesses for the EPSG code, for a user-defined geographic CRS or in
  GeoTIFF 1.0 mode, its default for a 2D CRS (same commit, lines 2030-2034
  and 3412-3460). So a GDAL-written file with 3072 = 32767 can carry 2062.
- **The callers' identification**: PROJ's `ProjectedCRS::identify`
  (`src/iso19111/crs.cpp`, <https://github.com/OSGeo/PROJ>), exposed as
  pyproj's `CRS.list_authority` and `CRS.to_epsg(min_confidence=70)`, and used by
  rasterio's `CRS.to_epsg` with the same default threshold of 70. PROJ's
  documented confidences: 100 name and definition match; 90 equivalent, names
  differ; 70 equivalent base CRS, conversion and coordinate system; 50
  "equivalent base ellipsoid and conversion, but the coordinate system do not
  match (e.g. different axis ordering)"; 25 similar names only. Its candidate
  search for a projected CRS (`AuthorityFactory::createProjectedCRSFromExisting`,
  `src/iso19111/factory.cpp`) is an SQL query that compares the *i*-th
  parameter of the conversion with the *i*-th parameter in the database (the
  two Lambert 2SP standard parallels excepted, which it takes in either order),
  so **the order of the parameters decides whether a candidate is found at
  all**. PROJ says itself that identify "might miss legitimate matches".
- **PROJ's equivalence tolerance**: `DEFAULT_MAX_REL_ERROR = 1e-10`
  (`include/proj/common.hpp`), the relative error up to which two parameter
  values are "equivalent".

What differs from that prior art, and why:

- **A missing `ProjectionGeoKey` (3074) is accepted**, which GeoTIFF 1.1
  does not allow for 3072 = 32767 (section 4, pinned with the checks).
- **`GeogTOWGS84GeoKey` (2062) is refused, not kept** (GDAL keeps it and
  builds a bound CRS). rasputin's result is an EPSG code, which has no room
  for a file's own datum shift; question 4 asks Ola whether to ignore it
  instead.
- **No fallbacks and no defaults** between parameter GeoKeys (libgeotiff's).
  Increment 11's ruling 3 is that the reader's job is refusal; a missing
  parameter is refused, naming the key. GDAL writes Lambert 2SP with the
  false-origin keys and transverse Mercator with the natural-origin keys, which
  is what this design reads.
- **The ellipsoid keys never override the coded datum** (libgeotiff's
  behaviour). They are checked against it, and a disagreement is refused.
- **Identification at read time** (GDAL leaves it to the caller), because the
  rest of rasputin speaks `EPSG:n`; and **not at PROJ's confidence threshold**:
  the Austrian file, built as below, gets confidence 50 for EPSG:31287 (the
  code's axes are north, east), so pyproj's and rasterio's default of 70
  returns no code for it (probe, section 3). Every candidate PROJ offers is
  judged by the rule in section 4 instead, whatever its confidence.

Not checked: what GDAL plus rasterio return for this file. GDAL is a
prohibited dependency (`CLAUDE.md` section 2), so it was not run.

*Legacy.* Nothing is carried across.
`git grep -nE "32767|[Uu]ser.?[Dd]efined|ProjCoordTrans|StdParallel|FalseOrigin|LambertConf|to_epsg|from_proj4|identify" legacy-archive -- legacy`
returned:

```
legacy/rasputin/geometry.py:213, 216, 220
legacy/rasputin/mesh.py:43, 46
legacy/rasputin/reader.py:60, 63, 64, 69, 70, 71, 72, 127, 128, 156, 294, 338, 412
legacy/rasputin/tin_repository.py:71, 100
legacy/rasputin/triangulate_dem.h:367
legacy/tests/test_gml_repository.py:44
```

`legacy/rasputin/reader.py@legacy-archive:60-72` only names the Lambert GeoKeys
in an enum; `GeoKeysInterpreter` has no handler for `ProjCoordTransGeoKey` or
any projection parameter (its handlers, `git grep -n "def _" legacy-archive --
legacy/rasputin/reader.py`, are lines 225-271: the EPSG windows, the
ellipsoid size, a UTM-zone decoder and a regex over the citation text). It
builds a PROJ string with no datum, which increment 11 ruled out and PR B's
rule calls never the same as an EPSG code. The rest of the hits are
`from_proj4` call sites of that string.

## 3. Finding: PR B's `same_crs` does not accept this file

PR B's question 1 (`docs/increments/python-audit.md`, section 9, "Questions
for Ola") says the CRS built from these GeoKeys "is then the same as
EPSG:31287 by B's rule (probe: EPSG:31287 as WKT2 with the file's
`lon_0=13.33333333300013` is equivalent and a `noop`)". That probe passed
only because the WKT2 kept EPSG's name, "MGI / Austria Lambert", and PROJ
built the operation from the database entry of that name. Renamed "unknown",
as any CRS built from GeoKeys is, the same WKT2 is **not** the same as
EPSG:31287: PROJ calls the two equivalent (leg 1 of `same_crs` holds), but
the operation between them is an inverse Lambert followed by a forward
Lambert, not `proj=noop` (leg 2 fails), because the two longitudes print
differently at PROJ's 15 significant digits (`13.3333333330001` and
`13.3333333333333`). Over the file's whole footprint that pipeline moves no
point by more than 2.6e-5 m. Probe (pyproj 3.8.0, PROJ 9.8.1), with
`same_crs` at `src_python/tin_engine/crs.py@65cd528:73-85`:
`docs/increments/geotiff-crs-by-parameters-probes/match_probe.py` prints
`B's same_crs(built Austrian, EPSG:31287): False`.

So the reader cannot decide "matches EPSG:n" by `same_crs` alone, or the file
this PR exists for stays refused. `same_crs` still does all the work
downstream: once the file reads as `EPSG:31287`, every comparison with a
river file, a domain or `--out-crs` is `same_crs` on that text, unchanged.
(Whether PR B's question-1 sentence is corrected on PR B's branch is the main
session's call; this file does not edit `python-audit.md`.)

## 4. The design

### Data flow

```
GeoKeys (tifffile's geotiff_metadata, a dict)          io/geotiff.py, layer L2
  │  3072 == 32767
  ▼
_parametric_epsg(geokeys) ── checks, each a GeoTiffError naming key, number, value
  │  base   = pyproj.CRS.from_epsg(2048)                (the datum, by code)
  │  method = METHODS[3075]                              (a table: data, not code)
  │  values = [geokeys[k] for k in method's keys, in EPSG's parameter order]
  │  built  = _projected_crs(base, method, values)       (PROJJSON -> pyproj.CRS)
  ▼
crs.epsg_matches(built) -> tuple[int, ...]             crs.py, layer L1 (cached by key)
  │  PROJ's EPSG candidates, every confidence
  │  kept: on built's own base geographic CRS, and equivalent once x-then-y
  ▼
none          -> refused: no EPSG code has these parameters
several, not
mutually same -> refused: ambiguous (same_crs between them)
otherwise     -> the lowest code n:  RasterMeta(epsg=n), meta.crs == "EPSG:n"
```

No path, no file and no CRS crosses into C++ (`CLAUDE.md` section 2, I/O
boundary): all of it is in `io/geotiff.py` and `crs.py`, and the core sees
only the decoded array, as today.

### `crs.py`: one new function, and `same_crs` shares its first leg

```python
def epsg_matches(crs: CRS) -> tuple[int, ...]:
    """The EPSG codes PROJ offers for `crs`, at every confidence, that sit on
    `crs`'s own base geographic CRS and that PROJ calls equivalent to `crs`
    once both are in x-then-y order (the first leg of `same_crs`), lowest
    first. Empty when there are none."""
```

- Candidates: `crs.list_authority(auth_name="EPSG", min_confidence=0)`.
- **On the same base**: `CRS.from_epsg(n).geodetic_crs.to_epsg(min_confidence=100)`
  equals `crs.geodetic_crs.to_epsg(min_confidence=100)`, the 2048 code. Without
  this filter PROJ's equivalence lets a CRS on another realisation through:
  UTM 35N by parameters on ETRS89 (4258) matched both EPSG:25835 (on 4258) and
  EPSG:3067 (on EUREF-FIN, 10690), and those two are not `same_crs` (the
  operation between them is not a `noop`). The file named 4258; it reads as a
  code on 4258.
- **Equivalent once x-then-y**: the first leg of `same_crs`, on the
  transformer's copies of the two CRSs. It is factored out of `same_crs` as a
  private helper so the two cannot drift; `same_crs` is unchanged in what it
  returns. The transformer stays in `crs.py`, the one `Transformer.from_crs`
  site (increment 11 ruling 8 and `tests/python/test_always_xy.py` hold).
- Why the first leg only, here and nowhere else: leg 2 (the operation is
  `proj=noop`) exists to catch a parameter PROJ keeps only in a PROJ string's
  remark (`python-audit.md` section 9). A CRS built from GeoKeys is a
  structured PROJJSON object with no remark, so there is nothing for leg 2 to
  catch, and section 3 shows what leg 2 does to it instead: it refuses a
  rounding at the 15th digit.

### `io/geotiff.py`: the method table

The table is data. Each method lists its parameters **in the EPSG method's own
parameter order**, because PROJ's candidate search compares parameters by
position (section 2); a parameter order of its own finds no candidate at all
(probe: pyproj's `LambertConformalConic2SPConversion`, which puts the two
standard parallels first, gives EPSG:31287 no candidate; the same values in
EPSG's order give it at confidence 50).

| 3075 | method (EPSG) | parameters in EPSG order: EPSG parameter, from GeoKey, unit |
|---|---|---|
| 1, transverse Mercator | Transverse Mercator (9807) | latitude of natural origin 8801 ← `ProjNatOriginLatGeoKey` (3081), degree; longitude of natural origin 8802 ← `ProjNatOriginLongGeoKey` (3080), degree; scale factor at natural origin 8805 ← `ProjScaleAtNatOriginGeoKey` (3092), unity; false easting 8806 ← `ProjFalseEastingGeoKey` (3082), metre; false northing 8807 ← `ProjFalseNorthingGeoKey` (3083), metre |
| 8, Lambert conic conformal with two standard parallels | Lambert Conic Conformal (2SP) (9802) | latitude of false origin 8821 ← `ProjFalseOriginLatGeoKey` (3085), degree; longitude of false origin 8822 ← `ProjFalseOriginLongGeoKey` (3084), degree; latitude of 1st standard parallel 8823 ← `ProjStdParallel1GeoKey` (3078), degree; latitude of 2nd standard parallel 8824 ← `ProjStdParallel2GeoKey` (3079), degree; easting at false origin 8826 ← `ProjFalseOriginEastingGeoKey` (3086), metre; northing at false origin 8827 ← `ProjFalseOriginNorthingGeoKey` (3087), metre |

Two methods only (question 1 for Ola): the Austrian file's, and transverse
Mercator, which covers UTM and the Gauss-Krüger grids (MGI's own Austrian
zones among them). Another method is refused, naming it; adding one is a row
and a red test.

`_projected_crs(base, method, values)` builds PROJJSON: `"type":
"ProjectedCRS"`, `"name": "unknown"`, `"base_crs"` the base's own
`to_json_dict()` (so its datum and ID are EPSG's), the conversion with the
method's name and EPSG ID and each parameter with its name, EPSG ID, value
and unit, and a Cartesian coordinate system of easting then northing in
metres (GeoTIFF 1.1's fixed order for a user-defined CRS). Then
`pyproj.CRS.from_json_dict`.

### `io/geotiff.py`: the checks, and the order they run in

`_header` dispatches on 3072: absent is the geographic path, unchanged;
32767 is `_parametric_epsg(geokeys)`; any other value is `_projected_epsg`,
unchanged. A 2048 of 32767 on the geographic path stays refused as today: a
geographic CRS by parameters is out of scope.

The checks run in the order 1, 2, 3, 4, 5, **8**, 6, 7, 9, 10, 11, 12, 13:
the angular unit (row 8) is checked before the ellipsoid and prime meridian
(row 6), because `GeogPrimeMeridianLongGeoKey` (2061) is in the unit 2054
names, so a file in grads is refused for its unit rather than for a meridian
that "disagrees". The row numbers stay as they are, since the red tests and
the code name them.

Every message, all thirteen rows, begins with the same words `P`, so the two
existing cases of
`test_refuses_missing_crs` that put 32767 in 3072
(`tests/python/test_io_geotiff.py@44fa7f5:884-888`, which look for `3072`,
`ProjectedCSTypeGeoKey` and `32767`) keep passing whatever check fires:

`P` = `ProjectedCSTypeGeoKey (3072) = 32767 (a CRS given by parameters)`

| # | condition | message |
|---|---|---|
| 1 | 3074 present and not 32767 | `{P}, but ProjectionGeoKey (3074) = {v} names its projection by EPSG code, which is not read; only a projection given by its parameters is` |
| 2 | 3075 absent | `{P} has no ProjCoordTransGeoKey (3075), so its projection method is unknown` |
| 3 | 3075 not in the table | `{P} uses ProjCoordTransGeoKey (3075) = {v} ({tifffile's name}); only transverse Mercator (1) and Lambert conic conformal with two standard parallels (8) are read` |
| 4 | 2048 absent, 32767, not a resolvable code, or not a geographic 2D CRS | `{P}: GeographicTypeGeoKey (2048) {is absent / = v, which is not an EPSG code / = v, a {type_name}}; its datum must be given there as an EPSG geographic CRS code` |
| 5 | `GeogGeodeticDatumGeoKey` (2050), `GeogPrimeMeridianGeoKey` (2051) or `GeogEllipsoidGeoKey` (2056) present | `{P}: {name} ({id}) is present; the datum is read only from GeographicTypeGeoKey (2048)` |
| 6 | `GeogSemiMajorAxisGeoKey` (2057), `GeogSemiMinorAxisGeoKey` (2058), `GeogInvFlatteningGeoKey` (2059) or `GeogPrimeMeridianLongGeoKey` (2061) present and not one finite number (`N` below); or present and not within 1e-10 relative of the coded CRS's own value | `{P}: {name} ({id}) = {v!r}, which is not one finite number`; or `{P}: {name} ({id}) = {v} disagrees with EPSG:{g}'s {ellipsoid name} ({expected})` (for 2061: its prime meridian) |
| 7 | `GeogTOWGS84GeoKey` (2062) present | `{P}: GeogTOWGS84GeoKey (2062) = {v} gives the file's own datum shift to WGS 84, which is not read; the datum is read only from GeographicTypeGeoKey (2048)` |
| 8 | `GeogAngularUnitsGeoKey` (2054) present and not 9102 | `{P}: GeogAngularUnitsGeoKey (2054) = {v}; only degrees (9102)`: the geographic path's words, after `P` |
| 9 | `ProjLinearUnitsGeoKey` (3076) absent | `{P}: ProjLinearUnitsGeoKey (3076) is absent, so the unit of its false easting and northing is unknown` |
| 10 | 3076 not 9001 | `{P}: ProjLinearUnitsGeoKey (3076) = {v}; only metres (9001)`: the coded path's words, after `P` |
| 11 | a parameter GeoKey of the method absent, or present and not one finite number (`N` below) | `{P}: {method} needs {name} ({id}), which is absent`; or `{P}: {name} ({id}) = {v!r}, which is not one finite number` |
| 12 | no EPSG code matches | `{P}: it is {method} on EPSG:{g} ({base name}), and no EPSG projected CRS on that datum has these parameters, so it cannot be named. Reproject the file to a CRS with an EPSG code first` |
| 13 | several match and are not all `same_crs` with the lowest | `{P}: its parameters match {EPSG:a, EPSG:b}, which are not the same CRS, so it cannot be named` |

Pinned with the checks:

- **`N`, one finite number** (rows 6 and 11): the value tifffile gives is an
  `int` or a `float`, not a `bool`, and `math.isfinite` holds. tifffile gives
  an inline SHORT as an `int`, one double from `GeoDoubleParamsTag` (34736) as
  a `float`, several as a `tuple`, and a key stored in `GeoAsciiParamsTag`
  (34737) as a `str`; the last two fail `N`. Without it a NaN, an infinity or
  a tuple reaches `float()` or PROJ and escapes as `TypeError` or pyproj's
  `CRSError`, which a caller catching `GeoTiffError` (a mosaic) does not
  catch. `{v!r}` prints `nan`, `inf` or the tuple. The `N` check of each key
  runs before that key's comparison (row 6) or its use (row 11).
- **The keys read as codes** (3072, 3074, 3075, 2048, 2054, 3076) go through
  `int()`, as 3072 does on the coded path (`_projected_epsg`). A value in one
  of them that is not an `int` (a tuple, a NaN or infinite double, a string)
  escapes as `ValueError`, `OverflowError` or `TypeError`, not
  `GeoTiffError`, on both paths alike. That is out of this PR; it is put to
  Ola as a small follow-up (section 10, question 5).
- **2054 absent reads as degrees**, as libgeotiff does (its angular unit factor
  starts at 1 degree in `GTIFGetDefn`). `GeogAzimuthUnitsGeoKey` (2060) is not
  read: neither method has an azimuth.
- **`GeogTOWGS84GeoKey` (2062) present is refused** (row 7; question 4 for
  Ola, default refuse). A TOWGS84 is the file's own datum shift to WGS 84: with
  it, the CRS is an ISO 19111 bound CRS, which is what GDAL builds for 3072 =
  32767 (section 2). Reading the file as the bare EPSG code would drop a shift
  the file states and let PROJ choose another when a reprojection runs; that
  is a silent change of datum handling, so it is refused rather than
  ignored. The Austrian file has no 2062. Matching the shift against EPSG's
  transformations instead is not in this PR.
- **Citations (1026, 2049, 3073) are never read**, as increment 11's ruling 6
  says: no free text decides a CRS.
- **Lowest code first** among several that are one CRS (question 2 for Ola):
  EPSG holds 239 such groups among its transverse Mercator codes (section 6);
  the choice changes only the label, since `same_crs` holds between them.
- **Rows 5, 6, 7 and 9 have no counterpart on the other two paths.** A file
  with an EPSG code in 3072 is read as on master; its GeoKeys beyond 3072, 3076
  and the vertical unit are not consulted, and it may still omit 3076.
- **A missing 3074 is accepted, which goes beyond OGC GeoTIFF 1.1** (clause 7
  wants 3074 populated for 3072 = 32767). Only 3074 present with a code is
  refused (row 1); the method is read from 3075 either way, so an absent 3074
  hides nothing the reader uses.

### The result, and the cost of reading it

`_parametric_epsg` returns the code, and `_header` builds `RasterMeta(epsg=n,
...)` as for a coded file. No new field on `RasterMeta`. The decoded pixels
are not moved: the file is read as EPSG:n to within PROJ's equivalence
(section 5).

A mosaic reads every tile's header, so the match is cached:
`functools.cache` on a private function, `_epsg_matches_cached`, keyed by the hashable tuple (2048
code, 3075 code, the parameter values in table order) returning
`epsg_matches`' tuple; the refusals stay outside it. One distinct CRS costs
about 1 to 30 ms in the probe (identify plus one transformer per candidate); every
further tile with the same keys costs a dictionary lookup. Pure and
thread-safe in effect (`functools.cache` may compute a key twice under a race,
with the same result), so it runs unchanged inside `asyncio.to_thread`.

### Boundaries

- `io.geotiff` imports `crs` (layer L2 importing L1: downward, allowed). Its
  row in `tests/python/test_layering.py@65cd528:57` becomes
  `"io.geotiff": "crs io.models"`, a `@tester` edit in the red commit. No
  `UPWARD` entry.
- The module docstring of `io/geotiff.py` (`@44fa7f5:8-12`) changes: it imports
  `io.models` and `crs`; the CRS is resolved through `from_epsg`, or for 3072 =
  32767 through this file's rule.
- `crs.py`'s public surface gains `epsg_matches`.

### Increment 11's rulings

- **Ruling 6 is amended** (`docs/increments/11-raster-ingestion.md`, ruling
  6 and refusal 12): a 3072 of 32767 may yield an accepted CRS, by this
  file's rule, and the result is always an EPSG code. For 3072 = 32767 the
  amendment overrides "only `ProjectedCSTypeGeoKey` can yield an accepted
  CRS, resolved through `pyproj.CRS.from_epsg`" (the CRS is built from the
  parameters, then named by a code); "`GeographicTypeGeoKey` is read only
  when 3072 is absent" (2048 is read, as the datum); refusal 12's "or its
  value is 32767"; and round 2's reading (d), "When 3072 is present, whatever
  its value, 2048 is not consulted" and that a 32767 message "must not claim
  2048 was consulted" (rows 4, 6 and 12 name 2048's code). "No proj4
  reassembly, no free-text ellipsoid regex" still holds.
- **Ruling 7 is amended for this path.** Projected is by construction
  (`built` is a `ProjectedCRS`), and metres are decided by 3076: present and
  9001 (rows 9 and 10), after which `built`'s axes are metres because the
  reader writes them so. Testing `built`'s axis units would test the reader's
  own literal. The matched code's axes are then metres too, because PROJ's
  equivalence compares axis units (probe: the Austrian CRS rebuilt with foot
  axes matches no code). The coded path keeps ruling 7 as it is.
- **Ruling 8 holds**: no `Transformer` in `io/geotiff.py`.
- Ruling 1 holds (no `_core`, no C++).

## 5. Constants and what "matches" tolerates

There is one constant, and it is PROJ's: **1e-10 relative**
(`DEFAULT_MAX_REL_ERROR`), used twice.

- In `epsg_matches`, through PROJ's equivalence. At the magnitudes of a
  national grid it is at most about 1 mm of false northing (10 000 000 m, a
  southern-hemisphere UTM zone), 0.04 mm of false easting at 400 000 m, and
  1.3e-9 degrees (about 0.1 mm) of a longitude of origin at 13.3, 1.8e-8
  degrees (2 mm at the equator) at 180. Probe at the
  Austrian file: a false easting off by 1e-5 m still matches EPSG:31287, off
  by 1 mm it matches nothing. Checked at the Austrian file's parameters and
  over EPSG's own 3 941 codes (section 6), not beyond.
- In check 6, with `math.isclose(..., rel_tol=1e-10, abs_tol=1e-10)`
  against the coded CRS's `ellipsoid.semi_major_metre`, `semi_minor_metre`,
  `inverse_flattening` and `prime_meridian.longitude`. The absolute 1e-10
  (metres, unity or degrees, the key's own unit) applies to all four keys and
  only matters for a value near 0, which in practice is a prime meridian of 0
  (a relative tolerance at 0 accepts only 0); at an axis of 6.4e6 m or an
  inverse flattening of 299 the relative term is the larger by far. At the
  Earth's semi-major axis that is 0.6 mm; the Austrian file's
  inverse flattening differs from Bessel 1841's by 1.1e-14 relative and passes.

No constant of rasputin's own.

## 6. The probe and the sweep

`docs/increments/geotiff-crs-by-parameters-probes/match_probe.py`, run from
the repository root with the project venv (pyproj 3.8.0, PROJ 9.8.1),
implements the rule above independently of the production code and prints:

| case | matches |
|---|---|
| the Austrian file (MGI, 4312), parallels 46 then 49, or 49 then 46 | EPSG:31287 |
| the same, false easting + 1e-5 m / + 1 mm | EPSG:31287 / none |
| the same on ETRS89 (4258) / on WGS 84 (4326) | EPSG:3416 / none |
| transverse Mercator at 15 E, 0.9996, 500 000: on ETRS89 / on WGS 84 | EPSG:25833 / EPSG:32633 |
| the same with scale factor 1 on ETRS89 | none (PROJ offers six candidates; none is equivalent) |
| MGI Gauss-Krüger M31 (13.333 E, 1, 450 000, -5 000 000) | EPSG:31258 |
| UTM 35N on ETRS89 | EPSG:25835 only (EPSG:3067, on EUREF-FIN, is filtered out) |
| Xian 1980 Gauss-Krüger CM 75E (4610) | EPSG:2338 and 2370, `same_crs` to each other: reads as 2338 |

With `--sweep`, it rebuilds every non-deprecated EPSG projected CRS whose
method is 9807 or 9802, on a geographic 2D base with an EPSG code, Greenwich,
parameters in degrees, metres or unity (3 941 codes; 839 skipped for other
units, 20 for another prime meridian), from its own parameters as a GeoTIFF
would carry them, and matches it back (27 s):

- 3 648 match only themselves;
- 54 match one other code that `same_crs` calls the same (the code's axis
  order differs: EPSG:3045, ETRS89 / UTM 33N with north first, reads as
  25833, east first, which is what GeoTIFF's fixed east-north order means);
- 239 match two or three codes, and all 245 pairs within those groups are
  `same_crs` (the lowest is chosen);
- **none matches nothing, and no group is ambiguous.**

So check 13 is never reached by any EPSG definition at PROJ 9.8.1; it stays
as the guard for a database where that changes, and its red test reaches it
by substituting `crs.epsg_matches` (section 8).

## 7. The NoData text: not in this PR

The Austrian file's `GDAL_NODATA` is `-3.4028234663852886e+38`, which is
float32's lowest value exactly. rasputin reads the tag itself
(`src_python/tin_engine/io/geotiff.py@44fa7f5:345-368`) and accepts it: probe,
`_nodata` on the file's page returns `(-3.4028234663852886e+38, 'tag')`.
tifffile parses the same tag for its own `page.nodata`, decides it is "not
castable to float32" (its `numpy.min_scalar_type` check), sets it to 0 and
logs a WARNING on the `tifffile` logger; with no logging configured, Python
prints that line on stderr on every open of the file. rasputin never reads
`page.nodata`; tifffile uses it only to fill sparse blocks, which `io/cog.py`
refuses and this file does not have (0 of its 28 375 blocks are empty).

It does not belong here: it is not about the CRS, it affects every float32
GDAL file with that sentinel whatever its CRS, and its fix (a filter on that
one tifffile message) is a separate concern with its own test. Question 3
asks Ola whether to do it as a separate small PR.

## 8. Red tests (`@tester`, one commit, before any code)

The fixtures extend `tests/python/geotiff_fixtures.py`, whose `micro_tiff`
writes GeoKeys only as inline SHORTs today: it needs `GeoDoubleParamsTag`
(34736) for the double-valued keys and `GeoAsciiParamsTag` (34737) for the
citations. **`AUSTRIA_KEYS`** is the table in section 1, all eighteen keys, on
the usual 3 x 4 baseline; each refusal fixture is it with exactly one change,
checked by `test_geotiff_fixtures.py` with tifffile alone, as the suite does
for every refusal. Not the 2 GB file.

1. **The Austrian GeoKeys read as EPSG:31287**: `read_meta`, `read_header`,
   `read_page` and `decode_dem` each give `epsg == 31287`, `crs ==
   "EPSG:31287"`, `geographic` False, and `decode_dem`'s array equal to the
   baseline's. Red at `0f6c271`: refused at 3072.
2. **Either order of the standard parallels** (46, 49 and 49, 46): 31287.
3. **Transverse Mercator by parameters**: 15 E, 0.9996, 500 000, 0 on ETRS89
   (4258) reads as 25833; on WGS 84 (4326) as 32633; MGI Gauss-Krüger M31 on
   4312 as 31258.
4. **One CRS under two codes reads as the lower**: Xian 1980 Gauss-Krüger CM
   75E's parameters (from `CRS.from_epsg(2338)`) on 4610 read as 2338.
5. **A CRS by parameters that matches nothing stays refused**, naming `3072`,
   `32767`, `no EPSG projected CRS` and the base code: the Austrian
   parameters on 4326; transverse Mercator at 15 E with scale factor 1 on
   4258 (PROJ offers candidates, so this pins that candidates are judged, not
   trusted); the Austrian false easting + 1 mm. And its pair: + 1e-5 m reads
   as 31287 (PROJ's tolerance, both sides).
6. **The datum filter**: UTM 35N on 4258 reads as 25835 (unfiltered, EPSG:3067
   would make it ambiguous).
7. **Ambiguous stays refused**: with `crs.epsg_matches` substituted
   (monkeypatch) to return `(3067, 25835)`, the message names both codes and
   `not the same CRS`. The one test that substitutes production code: no
   EPSG definition reaches this branch (section 6).
8. **Each check of section 4 fires on its one defect**, and its message holds
   the key's number, name and value: 3074 = 16033; 3075 absent; 3075 = 11
   (Albers); 2048 absent, 32767, 9999 (unresolvable), 25833 (projected);
   2050 present; 2057 = 6378137.0; 2059 = 299.0; 2062 = (577.326, 90.129,
   463.919, 5.137, 1.474, 5.297, 2.4232), seven doubles (row 7); 2054 = 9105
   (grad); 3076 absent; 3076 = 9002; 3084 absent. Every message starts with
   `P` (`startswith`, rows 8 and 10 included).
9. **The two existing user-defined cases** of `test_refuses_missing_crs` are
   unchanged and still pass (3072 = 32767 with nothing else, and beside 2048 =
   4326; check 2 fires).
10. **Cached**: the test first calls `cache_clear()` on the cached matcher
    (`tin_engine.io.geotiff._epsg_matches_cached`, the one private name the
    suite touches, so an earlier test's read cannot make it pass), then two
    reads of the Austrian fixture call PROJ's identify once (count calls to
    `pyproj.CRS.list_authority` by monkeypatch).
11. **The layering row**: `"io.geotiff": "crs io.models"`.
12. **The real file, when present** (skipped otherwise, as
    `tests/python/test_cli_catchment.py`'s Bygdin test is): `read_meta` on
    `../rasputin_data/austria_dgm10/dhm_at_lamb_10m_2018.tif` gives
    `EPSG:31287`, 58 061 x 31 793, `x_min` 108 880 (the tie point 108 875 plus
    half a cell: pixel is area). Header
    only, so it reads in well under a second.

Added at code review round 1 (`c8984e8`), red before `@developer` changes the
code:

13. **A parameter that is not one finite number is refused**, as a
    `GeoTiffError` (`pytest.raises(GeoTiffError)`, so a `TypeError` or
    `CRSError` fails the test) whose message starts with `P` and holds the
    key's name, number and value and `not one finite number`: row 11 with
    `ProjFalseOriginEastingGeoKey` (3086) = NaN, = inf, and = (400000.0,
    400000.0), two doubles; row 6 with `GeogSemiMajorAxisGeoKey` (2057) =
    NaN, = inf, and = (6377397.155, 6377397.155). Each fixture is
    `AUSTRIA_KEYS` with that one change, its witness checked with tifffile
    alone (the value survives the write as NaN, inf or a 2-tuple).
14. **A file in grads is refused for its unit**: `AUSTRIA_KEYS` with 2054 =
    9105 and 2061 = 2.5969213 (Paris, in grads) gives row 8's message
    (`GeogAngularUnitsGeoKey`, `2054`, `9105`) and not `disagrees`. Red at
    `5270a64`: row 6 ran first and refused the meridian.

Not invariant-critical: no mutation round (`docs/increments/README.md`, cost
constraints).

## 9. Net production lines

Estimate: **about +135, between 110 and 150**
(`python3 tools/count_loc.py <base> <head>`); measured at the green step,
`python3 tools/count_loc.py 65cd528 c8984e8`: +154. Code review round 1 adds
about 8 (the `N` check, shared by rows 6 and 11); measured at `35c44c8`: +159. Well under the 700 of
`CLAUDE.md` section 2. Counted as `ruff format` lays it out: a refusal whose
f-string passes 100 columns wraps to three to five lines, as
`_projected_epsg`'s do today.

- `crs.py`: `epsg_matches` about 10, the first leg factored out of `same_crs`
  about +3.
- `io/geotiff.py`: the method table about 15 (one line per parameter), the
  datum-key tuples about 6, `_parametric_epsg` with its thirteen checks about
  65 (about 4 lines per refusal, plus the lookups), `_projected_crs` about 25
  (the PROJJSON literal: type, name, base, conversion with method and
  parameters, and the two-axis coordinate system), the cached matcher about
  6, the dispatch in `_header` about +3.

The method table may be packed one parameter per line under `# fmt: off`
(`CLAUDE.md` section 2); the review says so if it is.

## 10. Questions for Ola (defaults hold until he answers)

1. **Which projection methods?** Only transverse Mercator (UTM, Gauss-Krüger)
   and Lambert conic conformal with two standard parallels (your Austrian
   file's)? Default: yes, those two; any other method is refused with a
   sentence naming it, and adding one later is a table row and a test.
2. **When the parameters match one CRS that EPSG lists under two numbers**
   (for example EPSG:2338 and 2370, two names for the same Chinese grid),
   read the file as the lower number? Default: yes. The alternative is to
   refuse such a file.
3. **The stray warning line.** Opening your Austrian file prints a tifffile
   warning on stderr about its NoData value, although rasputin reads that
   value correctly. Silence that one message in a separate small PR after
   this one? Default: yes.
4. **A datum shift written in the file.** Some files given by parameters
   also carry `GeogTOWGS84GeoKey` (2062), the file's own shift from its
   datum to WGS 84 (GDAL writes it in some cases; your Austrian file has
   none). Refuse such a file, with a sentence naming the key? Default: yes,
   refuse. The alternative is to ignore the key and read the file as the
   EPSG code, which lets PROJ pick its own shift when reprojecting, so a
   reprojected position need not be where the file meant it.
5. **A malformed code key.** A GeoKey that should hold a code (for example
   the projected CRS code, 3072, or the linear unit, 3076) but holds
   something else, such as two numbers, a non-finite number or text, makes
   the reader fail with a Python error instead of rasputin's sentence about
   the file. This is so for files with and without an EPSG code, and is
   older than this PR. Fix it in a separate small PR, so every such file is
   refused with a sentence naming the key? Default: yes, after this one.

## Citations this PR pins

`src_python/tin_engine/io/geotiff.py@44fa7f5` (byte-identical at `65cd528`:
`git diff --stat 44fa7f5 65cd528 -- src_python/tin_engine/io/geotiff.py` is
empty), `src_python/tin_engine/crs.py@65cd528`,
`tests/python/test_layering.py@65cd528`, `tests/python/test_io_geotiff.py@44fa7f5`,
`legacy/rasputin/reader.py@legacy-archive`.

## Review

- Design round 1 (`@reviewer`, at `fbbfb56`): three fixes (rows 7 and 9 did not start with `P`; the TOWGS84 prior-art claim was wrong for 3072 = 32767; the increment 11 amendments did not name what they override) and three suggestions (LOC recount, `cache_clear()` in red test 10, 3074 beyond GeoTIFF 1.1), all taken in the commit after `fbbfb56`.
- Design round 2 (`@reviewer`, at `900701c`): approved.
- Code round 1 (`@reviewer`, at `c8984e8`): changes requested: a NaN, infinite or multi-valued parameter escaped as `TypeError` or `CRSError` (rows 6 and 11 now refuse it; red test 13), five citations in increments 12, 15 and 25 unpinned (pinned to `65cd528`), a docstring tense (`@tester`); suggestions (2054 before row 6, red test 14; section 5's absolute tolerance) taken in the commit after `c8984e8`.
- Code round 2 (`@reviewer`, at `35c44c8`): approved; round-1 findings closed; +159 net against about +162.
