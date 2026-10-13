# Increment 43: height units, feet coordinates and more projection methods

**Status: built; code review round 1 asked for changes** (2026-10-13, on
`3fc5b60f`, branch `worktree-units`): the history labels in the new docstrings
and comments go, and `tools/scratch_copy.py`'s citation follows its moved
line, by `@developer`; this section and section 4.8 by the main session. Then
code review round 2, push and enqueue on Ola's yes. Python only, one PR.
Not refine or mesh code, so no `@perf` acceptance run (section 7).

## 1. What Ola asked, and what this fixes

Ola, 2026-10-13: "On the audit, go ahead with 1 and 5." In the main session's
list, 1 was the generality audit's row H1 and 5 its rows A6 and A8
(`docs/increments/generality-audit.md`, section 2). The audit's fix 5 also
names A4 and F6; they fit the same PR and are in it (section 7).

- **H1, heights in feet.** A `VerticalUnitsGeoKey` (4099) other than metre is
  refused, and an absent one is taken as metres
  (`src_python/tin_engine/io/geotiff.py@97c10c4d:175-179`). A DEM in feet whose
  file does not say so (US county and state lidar) is read as metres, so its
  heights and `--tolerance` are 3.28 times off. The record says "metres
  (assumed: the DEM file does not say)"
  (`src_python/tin_engine/cli.py@97c10c4d:1173-1175`), but nothing is printed.
  It is the one audit finding that gives a wrong answer without a word.
- **A8, coordinates in feet.** A projected CRS whose axes are in feet (US State
  Plane, EPSG:2229 and the other 956 EPSG feet CRSs) is refused
  (`src_python/tin_engine/io/geotiff.py@97c10c4d:338-345`).
- **A6, projection methods.** A projected CRS given by parameters (3072 =
  32767) is read for transverse Mercator and Lambert conic conformal 2SP only
  (`src_python/tin_engine/io/geotiff.py@97c10c4d:55-74`): polar stereographic,
  Lambert azimuthal equal area, Mercator, Albers and oblique Mercator files
  without an EPSG code are refused.
- **A4, sub-metre resampling.** The resampled grid's spacing is a whole number
  of metres, at least 1 (`src_python/tin_engine/target_grid.py@97c10c4d:91-101`):
  a 0.5 m DEM in degrees is resampled to 1 m.
- **F6, the feature read margin.** Features are read from the domain plus
  100 m (`src_python/tin_engine/feature_input.py@97c10c4d:54-56`); on a 1 km
  DEM that is a tenth of a cell, so a feature just outside the domain is
  dropped before the DEM's own one-cell grow.

## 2. Prior art: legacy and literature

*Literature.* Nothing new is claimed; every rule below is the GeoTIFF
standard's or its reference readers'.

- **OGC GeoTIFF 1.1** (OGC 19-008r4, <https://docs.ogc.org/is/19-008r4/19-008r4.html>):
  `VerticalGeoKey` (4096, `VerticalCSTypeGeoKey` in 1.0) names the vertical
  CRS by EPSG code, `VerticalUnitsGeoKey` (4099) its unit;
  `ProjLinearUnitsGeoKey` (3076) is the unit of the projected axes and of the
  false easting and northing; units are EPSG unit codes (9001 metre, 9002
  foot, 9003 US survey foot).
- **The EPSG Geodetic Parameter Dataset** and IOGP Guidance Note 7-2 (IOGP
  publication 373-7-2): each method's parameters, codes and order, which the
  table in section 4.5 copies (polar stereographic variants A 9810 and B 9829,
  Lambert azimuthal equal area 9820, Mercator variants A 9804 and B 9805,
  Albers 9822, Hotine oblique Mercator variants A 9812 and B 9815). The US
  survey foot is 1200/3937 m and the foot 0.3048 m (EPSG units 9003, 9002;
  pyproj gives `0.30480060960121924` for the first).
- **GDAL's writer**, `GTIFSetFromOGISDefnEx` in `frmts/gtiff/gt_wkt_srs.cpp`
  (<https://github.com/OSGeo/gdal>, commit `b0844ad1`, lines 2309-2343 Albers,
  2447-2483 Mercator, 2547-2576 polar stereographic, 2578-2654 Hotine,
  2892-2918 Lambert azimuthal): which GeoKey each parameter is written to.
  That is the convention this design reads, as the earlier design did for
  transverse Mercator and Lambert 2SP
  (`docs/increments/geotiff-crs-by-parameters.md`, section 2).
- **libgeotiff's reader**, `GTIFGetDefn` in `libgeotiff/geo_normalize.c`
  (<https://github.com/OSGeo/libgeotiff>, commit `75cfca53`, lines 1717-1772
  Mercator, 2064-2104 polar stereographic): falls back between alternative
  keys and defaults an absent one (scale 1, origin 0); GDAL then calls
  `SetPS(lat, long, scale)`, which is variant A at a latitude of plus or minus
  90 and variant B otherwise.

What differs, and why:

- **No fallbacks or defaults between keys** (libgeotiff has them). The
  earlier design ruled that a missing parameter is refused, naming its key;
  that holds for the new methods.
- **Feet coordinates are rescaled, not resampled.** The audit's fix shape for
  A8 was to require `--out-crs` in metres, as for a DEM in degrees. A
  projected CRS in feet is the same projection as its twin in metres with
  every coordinate times the unit factor (probe in section 4.4: the feet CRS
  and the twin built from it agree to 1.2e-9 m on EPSG:2229, 2263, 6539 and
  2240). So the grid is relabelled with no interpolation and no option
  needed. The figure is for the built twin; naming it by an EPSG code can
  move the grid by 0.1 mm (section 4.4).

*Legacy.* Nothing is carried across.
`git grep -n -i -E 'feet|foot|ftUS|VerticalUnits|LinearUnits|9002|9003|PolarStereo|AlbersEqual|LambertAzim|ObliqueMercator|CT_Mercator' legacy-archive -- legacy`
returned:

```
legacy/rasputin/reader.py:44, 61, 62, 86, 265
```

Lines 44-86 are an enum of GeoKey names; `_ProjLinearUnitsGeoKey`
(`legacy/rasputin/reader.py@legacy-archive:264-268`) maps 9001 to "m" and
raises `KeyError` for any other unit. The legacy read metres only.

## 3. The ideas file

Taken up: *Product review: readable code docs, and general applicability
beyond Ola's machine and data* (rows H1, A6, A8, A4, F6 of its audit).
Overlaps: *Catchments on a geographic DEM (GLO-30)*, which shares the
resampled grid that A4 changes; nothing in it is built here.

## 4. The design

### 4.1 Data flow

```
GeoTIFF header (GeoKeys)                                 io/geotiff.py
  │ horizontal: 3072 code, or 3072 = 32767 + METHODS row (4.5)
  │   axes in feet? ──> metre twin (crs.metre_twin, 4.4), factor f_xy
  │ vertical: 4096 / 4099 / caller's z_unit / assumption (4.3), factor f_z
  ▼
RasterMeta: x_min, y_max, delta_* already times f_xy; crs = the twin;
            file_crs = the file's own CRS when they differ;
            z_unit = the file's height unit; vertical_unit_assumed
  │
  ▼ pixels: decode_dem, cog.decode_window
array: valid cells times f_z, in place (NoData untouched)  ─> metres from here on
  │
  ▼ mosaic, resample, check points, _core: unchanged, all in metres
```

The conversion happens once, where the file is decoded. Nothing below `io/`
learns that a file was in feet; nothing crosses into C++ but metres
(`CLAUDE.md` section 2, I/O boundary). The caller's choice reaches `io/` as
a value next to `nodata`, the same way the NoData sentinel does today.

### 4.2 Units

One table in `io/geotiff.py`, data rather than code:

```python
#: EPSG unit code to (name pyproj uses, metres per unit, option spelling).
UNITS = {9001: ("metre", 1.0, "m"), 9002: ("foot", 0.3048, "ft"),
         9003: ("US survey foot", 1200 / 3937, "us-ft")}
ZUnit = Literal["metre", "foot", "US survey foot"]       # in io/models.py
```

A unit code outside the table is refused, naming the key and listing the
three. A unit read from a CRS's axes (pyproj's `unit_conversion_factor`) is
matched to a row by its factor, compared with `math.isclose(rel_tol=1e-12)`;
a factor outside the table (Clarke's foot, Indian yard, German legal metre:
32 non-deprecated EPSG projected CRSs) is refused.

### 4.3 Heights (H1)

**The file's unit**, read in `_header`. Key 4096 is read for the first time:
today the reader ignores it (it is named only in the refusal of a compound
CRS in 3072, `src_python/tin_engine/io/geotiff.py@3406a453:333-337`), and
4099 is read only to refuse a unit other than metre.

| 4096 (vertical CRS) | 4099 (unit) | Result |
|---|---|---|
| absent or 32767 | absent | silent, see below |
| absent or 32767 | in `UNITS` | that unit |
| absent or 32767 | other | refused, naming 4099 and the three accepted codes |
| an EPSG vertical CRS, axis up | absent | the axis unit, if in `UNITS`; else refused |
| an EPSG vertical CRS, axis up | in `UNITS` | refused unless the two agree |
| not an EPSG code, not a vertical CRS, or axis down (a depth) | any | refused, naming 4096 and what it holds |

Probe (pyproj 3.8.0, PROJ 9.8.1): 5703 NAVD88 height is metre, 6360
NAVD88 height (ftUS) is US survey foot, 8228 NAVD88 height (ft) is foot,
5941 NN2000 height and 3855 EGM2008 height are metre, 6358 NAVD88 depth
(ftUS) points down.

**The caller's unit**, `z_unit: ZUnit | None`, from `--z-unit {m,ft,us-ft}`
on `mesh`, `catchment` and `station-catchments` (every command with `--dem`),
or from the catalogue entry of a cached source (`RemoteSource.z_unit`, new,
`"metre"` for `anadem-v1` and `glo30`, `None` for class sources). It travels
as `nodata` does: `DemRequest.z_unit`, `repository_for`, `TiffDemRepository`,
`CacheRepository`, `read_header`, `read_page`, `decode_dem`, `read_meta`, and
`fetch/plan.py`'s `parse_prefix`, which calls `read_page` on a remote file's
first bytes. It is one more keyword beside `nodata` (and beside `threads`,
which increment 38 added to `TiffDemRepository`, `CacheRepository` and
`cog.decode_window`); the two do not interact.

| file says | caller says | horizontal unit | Result | `vertical_unit_assumed` |
|---|---|---|---|---|
| a unit | nothing or the same | any | the file's | False |
| a unit | another | any | refused: "the file says its heights are in X; --z-unit says Y" | |
| silent | a unit | any | the caller's | True |
| silent | nothing | metre or degree | metre | True |
| silent | nothing | foot or US survey foot | refused (below) | |

The refusal: "the DEM's coordinates are in US survey feet (EPSG:2229) and the
file does not say what unit its heights are in; give --z-unit us-ft, ft or m".
A feet DEM whose height unit is unstated is the case H1 names, so it is never
guessed. A metre or degree DEM with silent heights stays metres, as today,
because every DEM on Ola's disk is one (all 254 DTM10 tiles, the Austrian
openDEM and BEV files, ANADEM and GLO-30 carry neither 4096 nor 4099; probe
with tifffile 2026.9.20 over `../rasputin_data`).

**Scaling.** After the pixels are decoded and promoted, in `decode_dem` and
in `io.cog.decode_window` (the two places a DEM array is made from file
bytes), one helper. It multiplies in float64 and casts to the array's dtype
once (a float32 DEM is promoted to float64 for the product, then stored back),
so the figures of section 6 are float64 products rounded once to float32:

```python
def to_metres(array: npt.NDArray[Any], meta: RasterMeta) -> None:
    """Valid cells of `array` times the height unit's metres, in place; NoData
    cells keep the sentinel. Nothing when the unit is the metre. A converted
    height that equals the sentinel is refused (GeoTiffError)."""
```

The metre case returns before touching the array, so every metre DEM decodes
to the same bytes as today. `RasterMeta.nodata` keeps the file's sentinel.

**Mosaic.** `plan_mosaic` refuses selected tiles whose `z_unit` differ
("tiles give heights in different units: A in feet, B in metres"); their
arrays would be right, but the record names one unit. `_grid` and `resample`
copy `z_unit` and `file_crs` into the meta they build, so the single-tile path
of `assemble` still matches (`plan.tiles[0].meta == plan.meta`).

**Speed.** For a metre DEM the new work is two GeoKey lookups and one axis
unit check per file header, and `to_metres` returns at its first line. For a
feet DEM it is one pass over the array: 0.07 s for a 5000 x 5000 float32 tile
(the masked float64 multiply and cast, numpy 2 on this Mac, 2026-10-13), once
per decoded window.

**What the user sees.** The record's `dem_vertical_unit` becomes one of
"metres", "metres (assumed: the DEM file does not say)", "metres, converted
from US survey feet" or "metres, converted from feet (given with --z-unit)".
When the unit is assumed, `mesh` prints one line to stderr:
"DEM heights read as metres; the file does not say. If they are feet, give
--z-unit ft or us-ft." Ola ruled to print it (section 8).

### 4.4 Coordinates in feet (A8)

A new function in `crs.py`, the one module that builds CRSs:

```python
def metre_twin(crs: CRS) -> tuple[CRS, float]:
    """`crs` with its axes in metres and its linear conversion parameters
    (false easting and northing, easting and northing at an origin) converted
    to metres, named "<name>, in metres"; and the metres per unit of its axes.
    `crs` itself and 1.0 when its axes are already in metres."""
```

`_projected_epsg` (an EPSG code in 3072) and `_parametric_epsg` (3072 =
32767) both end the same way: when the axes are in a `UNITS` foot, the meta
gets the twin, named by `crs.epsg_matches(twin)` when that gives one CRS
(`EPSG:n`), else carried as its own text (`crs_label(twin)`, PROJJSON) with
`epsg` None.

**A named twin is PROJ's equivalent, not an exact copy.** `epsg_matches`
takes a code when PROJ calls the two CRSs equivalent, and PROJ compares
parameters to 1e-10 relative. US State Plane feet CRSs round their false
origin to a thousandth of a foot: EPSG:2229's 6561666.667 US ft is
2000000.0001016 m, where EPSG:26945 has 2000000 m. So labelling the grid
EPSG:26945 moves every node by 0.10 mm in easting and 0.10 mm in northing
(the twin's own PROJJSON would not). Probe (pyproj 3.8.0, PROJ 9.8.1, false
origin in feet times the factor, minus the metre code's): 2229 to 26945 and
6420 to 6419, +1.016e-4 m in both; 2240 to 26967, -1.016e-4 m in easting,
0 in northing; 2263 to 32118, 5.8e-11 m. Accepted, and stated here: 0.1 mm is
under a two-thousandth of the finest cell read (0.5 m), and the EPSG name
lets `--out-crs`, the GeoPackage row and every tool downstream recognise the
grid. Where PROJ's test fails (2881 below, 5e-10 relative), the twin keeps its
own text and moves nothing. `x_min`, `y_max`, `delta_x` and `delta_y` are each multiplied by
the factor once, after `_placement`. `file_crs` keeps the file's own label
(`EPSG:2229`).

On the parametric path, 3076 = 9002 or 9003 is now read: the linear
parameters are converted to metres before the CRS is built, and the axes are
metres, so `epsg_matches` sees a metre CRS. 3076 present on the EPSG-code path
must agree with the code's axes, else refused.

Probe (`allfeet.py`, run from the scratchpad over every non-deprecated EPSG
projected CRS): 957 CRSs have axes in foot or US survey foot; for 458 the
twin is named by `epsg_matches` (2229 to 26945, 2263 to 32118, 2240 to
26967, 6420 to 6419), and 499 have no EPSG twin (2881, NAD83(HARN) / Florida
East (ftUS), is one: its false easting, 656166.667 ftUS, is 200000.0001 m,
5e-10 relative from EPSG:2777's 200000 m, past PROJ's 1e-10). The unnamed
twin is carried as PROJJSON: `same_crs` of it with itself is True, a point
moved through it lands on the same longitude and latitude as through the feet
CRS to 0.0 degrees, and `srs_entry` gives the GeoPackage row 100000 named
"NAD83(HARN) / Florida East (ftUS), in metres". Downstream it is the path a
CRS without an EPSG code already takes (the suggested `--out-crs`).

The record gains `dem_file_crs`, written only when the file's CRS differs
from the CRS the grid is read as (a feet file); for every other DEM the entry
is absent, as today: "EPSG:2229, in US survey feet; read as
EPSG:26945, the same projection in metres". `--out-crs` and `--tolerance`
stay metres; feet `--out-crs` is still refused (H2).

### 4.5 Projection methods given by parameters (A6)

`METHODS` gains rows; the key becomes the method's own name, and a small
function picks the row from 3075 and, for two codes, the keys present:

| 3075 | Chosen when | EPSG method | Parameters, in EPSG order, from GeoKey |
|---|---|---|---|
| 1 | always | Transverse Mercator 9807 | (unchanged) |
| 8 | always | Lambert Conic Conformal (2SP) 9802 | (unchanged) |
| 3 | always | Hotine Oblique Mercator (variant A) 9812 | 8811 from 3089 `ProjCenterLat`, 8812 from 3088 `ProjCenterLong`, 8813 from 3094 `ProjAzimuthAngle`, 8814 from 3096 `ProjRectifiedGridAngle`, 8815 from 3093 `ProjScaleAtCenter`, 8806 from 3082, 8807 from 3083 |
| 9815 | always | Hotine Oblique Mercator (variant B) 9815 | as variant A, but 8816 from 3082 and 8817 from 3083 (GDAL writes the projection-centre easting and northing to the false easting keys) |
| 7 | 3092 `ProjScaleAtNatOrigin` present, 3078 absent | Mercator (variant A) 9804 | 8801 from 3081 `ProjNatOriginLat`, 8802 from 3080 `ProjNatOriginLong`, 8805 from 3092, 8806, 8807 |
| 7 | 3078 `ProjStdParallel1` present, 3092 absent | Mercator (variant B) 9805 | 8823 from 3078, 8802 from 3080, 8806, 8807; 3081, if present, must be 0 |
| 10 | always | Lambert Azimuthal Equal Area 9820 | 8801 from 3089, 8802 from 3088, 8806, 8807 |
| 11 | always | Albers Equal Area 9822 | 8821 from 3081, 8822 from 3080, 8823 from 3078, 8824 from 3079, 8826 from 3082, 8827 from 3083 |
| 15 | 3081 is plus or minus 90 | Polar Stereographic (variant A) 9810 | 8801 from 3081, 8802 from 3095 `ProjStraightVertPoleLong`, 8805 from 3092, 8806, 8807 |
| 15 | 3081 is not plus or minus 90 | Polar Stereographic (variant B) 9829 | 8832 from 3081, 8833 from 3095, 8806, 8807; 3092, if present, must be 1 |

Refused, each naming the keys: 3075 = 7 with both or neither of 3078 and
3092; 3075 = 7 variant B with 3081 other than 0 (variant B has no latitude of
origin parameter, so a value other than 0 would otherwise be dropped silently);
3075 = 15 not at a pole with 3092 other than 1; any listed key absent;
any other 3075 (the message lists the ten readable methods).

**Matching widened for polar axes.** EPSG's polar CRSs label their axes by
meridian ("north along 90°E"), so PROJ never calls a CRS built with east and
north axes equivalent to them, though the operation between the two is
`proj=noop`. `crs.epsg_matches` therefore accepts a candidate on the same
base CRS **whose axes carry a meridian** when the operation between the two
is `noop`, as well as on today's equivalence. Candidates without meridian
axes are judged exactly as today, so the Austrian file and every existing
match are unchanged. The widened rule is tested twice (section 6): test 10's
polar files reach their code only through it, and test 10b calls
`epsg_matches` directly on a polar CRS and on the same CRS moved 1 m.

Probe (`rule.py` in the scratchpad: parameters taken from the EPSG entry,
written to the GeoKeys above, built as the reader will, matched): every one
of 24 CRSs reads as its own code: 3031, 3413, 3995, 3996 (polar B); 32661
and 5041 both read as 5041, 32761 and 5042 as 5042 (polar A: each pair is one
CRS with two axis orders, `same_crs` True, the lowest code taken as the
earlier design rules); 3035, 6931, 6932, 3571, 3574 (Lambert azimuthal); 3395,
3832 (Mercator A); 3994 (Mercator B); 3338, 5070, 3577 (Albers); 3375, 3079
(Hotine A); 2056, 21781, 29873 (Hotine B). With the longitude parameter moved
by 1 degree, all 24 are refused (no match).

### 4.6 Resampled spacing below a metre (A4)

`default_spacing` keeps `max(1, round(s))` for a source spacing `s` of at
least 1 m and returns `2 ** round(log2(s))` below it, at least 2^-3 m
(0.125 m; was 2^-10, section 9.5). `TargetGrid.spacing` becomes a float that
must be a whole number from 1 to 1000 or 2^-k for k from 1 to 3. Powers of two
keep the lattice exact: node `(R, K)` at `(K h, -R h)` is a multiple of
2^-3 m, exact in float64, and for any coordinate below 2^42 m the start mesh's
1 mm snap leaves it in place (measured, section 9.5; being a whole number of
millimetres alone is not enough, section 9.3). The record prints the spacing with `:g`, so a whole
spacing prints as today ("10 m").

| source spacing | today | after |
|---|---|---|
| 0.0000045 degrees at 50° N (0.500531 m) | 1 m | 0.5 m |
| 1/9 arc-second at 50° N (3.432996 m) | 3 m | 3 m |
| 1 arc-second at 50° N (30.896962 m) | 31 m | 31 m |

### 4.7 Feature read margin (F6)

`MARGIN` becomes `max(100 m, 2 cells)`, with the cell the larger of the
DEM's `delta_x` and `delta_y` in its CRS (the target grid's on the
reprojected path). `open_features` and `feature_input.source_region` (in
`src_python/tin_engine/feature_input.py`; `target_grid.source_region` is
another function and does not change) take `cell: float`.
At 10 m (DTM10) and 30 m (GLO-30 resampled at 31 m) the margin stays 100 m,
so those meshes do not change; at 1 km it is 2 km. `DENSIFY` (1 km) stays: it
bounds the chord of a curved edge after the move, not the margin.

### 4.8 What does not change

`_core`, its bindings and every C++ file. The mosaic's seam threshold (1 mm)
is now in metres for every DEM (audit F4 follows). For geographic and metre
DEMs, the decoded arrays, byte for byte. The mesh files too, with two
exceptions, both intended: a resampled source (geographic or projected) whose
spacing falls below about 0.71 m now gets a power-of-two spacing (section
4.6), and 0.125 m below 2^-3 m (section 9.5); and with features given, a DEM
whose cells are over 50 m gets a feature read margin of two cells instead of
100 m (section 4.7). The quick check's six meshes are unchanged.

The run record does change in two places. A cached source (`glo30`,
`anadem-v1`) now states metres in its catalogue entry, so its
`dem_vertical_unit` reads "metres" where it read "metres (assumed: the DEM
file does not say)". A local DEM's record is unchanged: `dem_file_crs` is
written only for a feet file.

## 5. Constants

| Constant | Value | Scale assumed | Checked at |
|---|---|---|---|
| foot, US survey foot | 0.3048, 1200/3937 m | exact definitions | n/a |
| unit factor match | `rel_tol=1e-12` | factors of order 0.3 to 1 | pyproj's 957 feet CRSs |
| pole test for 3075 = 15 | abs(lat) == 90 exactly | GDAL writes the EPSG value | 32661, 5041, 32761, 5042 |
| sub-metre spacing floor | 2^-3 m (section 9.5) | the 1 mm snap; coordinates below 2^42 m | 0.5 m (the A4 probe); 2^-1 to 2^-3 at origins to 7.7e6 m; multiples of 2^-3 m below 2^42 m (section 9.5) |
| feature margin | max(100 m, 2 cells) | DEM cells 1 m to 1 km | 10 m and 31 m unchanged |
| EPSG name for a feet twin | PROJ's equivalence, 1e-10 relative | false origins up to 1e7 m, so 1 mm at most | 2229, 2240, 6420 (0.10 mm), 2263 (6e-11 m) |

## 6. Red suite for `@tester`

Fixtures are micro-GeoTIFFs written with tifffile (no GDAL, no rasterio),
each 3 x 3 nodes unless said; the figures below were computed with numpy
2 and pyproj 3.8.0 / PROJ 9.8.1. Heights are float64 products (cell times
the factor) cast once to float32, as `to_metres` does (section 4.3); a
float32 multiply gives different values in four cells of tests 1 and 2.

**The invariant-critical suite** (mutation testing required) is
`tests/python/test_io_units.py`, tests 1 to 8. Mutation targets: each
factor in `UNITS`; the precedence of file over caller; the contradiction
refusal; the NoData mask in `to_metres`; the factor applied to `delta_*` and
not only to the corner; the metre case returning without a write.

1. **Heights in US survey feet.** int16 `[[100,200,300],[400,500,600],[700,800,-9999]]`,
   GDAL_NODATA `-9999`, EPSG:25833, 4099 = 9003: `decode_dem` gives float32
   `[[30.480060577392578, 60.960121154785156, 91.440185546875], [121.92024230957031, 152.40029907226562, 182.88037109375], [213.3604278564453, 243.84048461914062, -9999.0]]`,
   `z_unit == "US survey foot"`, `vertical_unit_assumed` False. The same
   through `read_page` plus `cog.decode_window`.
2. **Heights in feet**, 4099 = 9002: `[[30.479999542236328, 60.959999084472656, 91.44000244140625], [121.91999816894531, 152.39999389648438, 182.8800048828125], [213.36000061035156, 243.83999633789062, -9999.0]]`.
   **A converted height equal to NoData**: int16 `[[-10000, 1], [2, 3]]`
   (2 x 2), GDAL_NODATA `-3048`, 4099 = 9002 is refused (GeoTiffError):
   -10000 ft is -3048.0 m, in float64 and in float32.
3. **4096 alone**: 6360 gives test 1's array; 8228 test 2's; 5703 the
   input unchanged. **Refused**: 4096 = 6360 with 4099 = 9001; 4096 = 6358
   (a depth); 4096 = 4326 (not vertical); 4096 = 1 (no such code); 4099 =
   9005; each message names the key by number.
4. **Caller's unit**: silent file plus `z_unit="US survey foot"` gives test
   1's array with `vertical_unit_assumed` True; 4099 = 9001 plus
   `z_unit="foot"` is refused naming both.
5. **Silent metre file**: bytes identical to today's decode, `z_unit ==
   "metre"`, `vertical_unit_assumed` True; the array is not written to (pass a
   read-only array to `to_metres`).
6. **Feet coordinates, named twin.** EPSG:2229 in 3072, 4099 = 9003, tie
   (0, 0) at (6 500 000, 2 100 000) ftUS, scale 3 x 3 ftUS, pixel is area:
   `crs == "EPSG:26945"`, `file_crs == "EPSG:2229"`, `x_min ==
   1981204.419608839`, `y_max == 640080.8229616459`, `delta_x == delta_y ==
   0.9144018288036575` (each `(file value) * 1200/3937`, corner shifted half a
   cell first).
7. **Feet coordinates, unnamed twin.** EPSG:2881: `epsg is None`, the crs
   text parses to a CRS named "NAD83(HARN) / Florida East (ftUS), in metres",
   and the point (700 000, 600 000) ftUS and its metre image map to the same
   longitude and latitude (difference 0.0).
8. **Feet coordinates, silent heights**: test 6's file without 4099 is
   refused, the message naming `--z-unit`; with `z_unit="US survey foot"` it
   reads.
9. **Parametric in feet.** 2229's parameters by GeoKeys (3075 = 8, 2048 =
   4269, 3076 = 9003, false origin 33.5, -118, parallels 35.4666666666667 and
   34.0333333333333, 6561666.667 and 1640416.667 ftUS): reads as
   `EPSG:26945`.
10. **Methods** (`tests/python/test_io_geotiff_methods.py`), one file per
    row, parameters copied from the EPSG entry into the keys of section 4.5:
    3031, 5041, 3035, 6931, 3395, 3994, 3338, 3375, 2056 each read as their
    own code; 32661's keys read as 5041. Each with its longitude parameter
    moved 1 degree is refused as "cannot be named". Refused by shape: 3075 =
    7 with both 3078 and 3092, and with neither; 3075 = 15 at -71 with 3092 =
    0.97; 3075 = 7 with 3078 and 3081 = 10 (Mercator B with a latitude of
    origin); 3075 = 12 (azimuthal equidistant), its message listing the ten.
10b. **The noop rule of `epsg_matches`** (in `tests/python/test_crs.py`):
    `CRS.from_proj4("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0
    +y_0=0 +datum=WGS84 +units=m +no_defs")` (east and north axes; the
    operation to EPSG:3031 is `proj=noop`, probe 2026-10-13) gives `(3031,)`;
    the same with `+x_0=1` (operation a two-step pipeline) gives `()`.
11. **`epsg_matches` unchanged off the poles**: the Austrian keys of
    `docs/increments/geotiff-crs-by-parameters.md` still read as 31287, and
    UTM 35N by parameters on 4258 still reads as 25835 alone.
12. **Mosaic**: two aligned tiles, one with 4099 = 9001 and one 9002, are
    refused naming both files. A single feet tile through `assemble` returns
    the loaded tile (no canvas), as a metre one does.
13. **A4**: `default_spacing` of a geographic meta with `delta_y =
    0.0000045` at (10, 50) is 0.5; with 1/9 arc-second 3, with 1 arc-second
    31; a projected 10 m meta 10. `TargetGrid(spacing=0.75)` and `spacing=0`
    are refused; `0.5` and `2**-3` accepted, `2**-4` refused (section 9.5).
14. **F6**: `source_region` with `cell=1000` reaches 2 km beyond the domain
    (a point 1.5 km out is inside the region), with `cell=10` 100 m (the same
    point is outside).
15. **CLI** (`tests/python/test_cli_units.py`, on a 4 x 4 fixture): `mesh`
    on test 8's file exits 2 with the `--z-unit` sentence; with `--z-unit
    us-ft` it writes a mesh whose z values are the metre heights and whose
    record says "metres, converted from US survey feet (given with --z-unit)";
    on a silent metre file stderr carries the one assumption line once;
    `--z-unit ft` on a 4099 = 9001 file exits 2. `catchment --help` and
    `station-catchments --help` list `--z-unit`.
16. **Catalogue**: `SOURCES["glo30"].z_unit == "metre"`; a cached-source
    run's record says "metres" with no "assumed".

## 7. Size, acceptance, order

About 250 net production lines (`python3 tools/count_loc.py`), well under
700: `io/geotiff.py` about 135 (units, heights, twin wiring, six method rows
packed one parameter per line under `# fmt: off` as the existing table is,
the selector), `crs.py` 20, `cli.py` 45 (the option on three commands, the
stderr line, two record entries), `io/repository.py` 10, `mosaic.py` 8,
`target_grid.py` 14, `feature_input.py` 8, `sources.py` 4, `io/models.py` 5,
`dem_input.py` 4, `io/cog.py` 2. Worst case 330.

No `@perf` acceptance: nothing under `include/terrain/refinement/` or
`include/terrain/mesh/` changes, and a metre DEM decodes to the same bytes.
The review checks that the cases of `docs/benchmarks/quick/cases.toml`
(`tools/bench_quick.py`) give unchanged meshes, as increment 41's did.

**Speed.** No new step runs on the quick check's catchments: all are metre
DEMs (DTM10 tiles and folder) or `glo30` through the cache, whose catalogue
states metres, so `to_metres` returns before touching the array, and their
cells (10 m, 31 m) keep the 100 m feature margin. Per file header the reader
adds two GeoKey lookups and one axis unit check. A feet DEM pays one masked
float64 multiply per decoded window (0.07 s per 25 million cells, section
4.3).

Order: increment 41 is in the base. Increment 38 is merged into this branch
(`3406a453`); it added `threads` to the repositories and to
`cog.decode_window`, the call chain `z_unit` travels, so the red suite is
written against that signature. Increment 40 (PR #245, in the merge queue)
edits the catchment commands in `cli.py` where `--z-unit` is added: a clash
in text, not in meaning. If 40 lands first, merge master before the red
step; otherwise 40's branch resolves it.

## 8. Pinned beyond the audit, and questions for Ola

Pinned: the option is `--z-unit {m,ft,us-ft}` (the audit wrote
`--vertical-unit`); a feet DEM with silent heights is refused rather than
guessed; feet coordinates are rescaled to the metre twin rather than resampled
behind a required `--out-crs`; Hotine oblique Mercator is in (the audit listed
it as broken but left it out of the fix); a mosaic of tiles in two height
units is refused.

**Ruled: the stderr line is printed.** Every DEM on Ola's disk (DTM10, the
Austrian files, ANADEM) leaves its height unit unstated, so `rasputin mesh`
prints one extra stderr line on each of those runs: "DEM heights read as
metres; the file does not say. If they are feet, give --z-unit ft or us-ft."
Ola, 2026-10-13: "great idea regarding std-err. Use it." GLO-30 and ANADEM
through the cache do not print it, because their catalogue entries state
metres.

Pinned after design review round 1: a feet CRS's twin is named by the EPSG
code PROJ calls equivalent, accepting the 0.1 mm false-origin offset of
section 4.4.

## 9. Addendum after the green step

`@developer`'s green commit `0f67a567` (300 net lines, Python only) left
tests red. Full suite in the worktree's venv at `0f67a567`
(`.venv/bin/pytest -q --no-cov`, 2026-10-13): 7 failed, 7655 passed, 18
skipped. The commit message reports 8 failed and 7654 passed, but names only
the 7 below; the eighth is in 9.4. Each is settled here; the changes are
`@tester`'s (tests) unless marked `@developer`.

### 9.1 The word "feet" on stderr: the ban narrowed

Five cases of `test_cli_mesh_plain_output.py::TestStderr::test_no_banned_word_on_a_successful_run`
fail on the line of section 4.3 that Ola ruled is printed ("If they are feet,
give --z-unit ft or us-ft."). Increment 25's rule 6 bans `feet` because it
meant the refinement's jargon ("0 constraint feet, 0 feet refused").

**Ruled** (Ola, 2026-10-13, recorded by the main session; his words were
"yes, narrow the ban"): the ban covers that jargon only, and "feet" as a unit
is allowed. In both copies of the pattern,
`tests/python/test_cli_mesh_plain_output.py` (`BANNED`) and
`tests/python/test_run_record.py` (`BANNED`, which a feet DEM's record value
"metres, converted from US survey feet" would trip), replace `\bfeet\b` by
`\bconstraint feet\b|\bconstraint-feet\b|\bfeet refused\b` (the hyphenated
spelling too, so the jargon cannot come back as a compound):

```python
BANNED = re.compile(
    r"(?i)\buncovered\b|\bconstraint feet\b|\bconstraint-feet\b|\bfeet refused\b|"
    r"\bstride\b|\bchains\b|"
    r"\bcoincident\b|\bcarved\b|\bnoded\b|R-tree|\bsnap|\bDEM holes\b|\bvoid\b"
)
```

The comment above each says the narrowing is increment 43's, by Ola's ruling.
No production change.

### 9.2 The catchment fingerprint: updated

`test_catchment.py::test_no_reach_gives_22s_result_bit_for_bit` hashes
`repr` of the result's fields, `result.meta` among them, and `RasterMeta`
gained `z_unit='metre'` and `file_crs=None`. **Update the two constants**;
do not exclude the fields. The fingerprint is meant to trip on any change to
the result, the meta included; an exclusion list would hide a later wrong
`z_unit`. The evidence that nothing else moved, re-run here (a scratch copy of
`digest` that removes the text `, z_unit='metre', file_crs=None` from the repr,
and asserts it was there): the stripped digests equal the old constants, so
the mask, outline, counts and every other field are unchanged.

| constant | old | new |
|---|---|---|
| `BOWL_LAKE_DIGEST` | `0efdc506ba2c0d63409f13da1bd761388b34762d2cc10bd26e189755cd902c5f` | `10971eaf66ab83436437537731c04863980e617341f1373761f18b7a70ed9e29` |
| `BOWL_POINT_DIGEST` | `02ca6c38e46f445b5bf01996531616c130a0942766478a26d69b24346ac1740c` | `ea4f31d3aa5e85fa7f7d0c0090b8415f5de743ef5ae3072330a81a3fbe00a407` |

A comment beside them names the two fields and this check.

### 9.3 A feet DEM's mesh: the 1 mm snap, its own increment

`test_cli_units.py::TestTheFeetDem::test_given_us_ft_the_mesh_is_in_metres`
cannot pass. The noder snaps every start-mesh vertex to a 1 mm lattice
anchored at zero (`DEFAULT_SNAP_SPACING`,
`src_python/tin_engine/cli.py@0f67a567:280`; the anchor is fixed by
`include/terrain/core/snap_grid.hpp@0f67a567:6-11`, which forbids an offset).
The fixture's nodes (3 ftUS cells, 0.9144018 m) are off that lattice, so
boundary vertices move outward by up to 0.5 mm, leave the raster and are
dropped. Increment 27 recorded the effect ("A limit of the exact rule"),
then for made-up grids only.

Measured (NumPy emulation of `SnapGrid::snap` and `world`, `llround` half
away from zero; nodes as `x_min + c * delta_x`):

| grid | anchor at zero: outside / nodes | anchored at the DEM's origin: outside / nodes |
|---|---|---|
| the test's 4 x 4, 3 ftUS cells | 7 / 16 | 0 / 16 |
| metre 4 x 4 at x = 500000.0004 | 4 / 16 | 0 / 16 |
| metre, 10 m, 1 m and 0.5 m cells, whole-metre origin, 300 x 300 | 0 / 90000, none moved | 0 / 90000, none moved |
| 300 x 300, 1 ftUS cells, origin (300000, 200000) ftUS | 300 / 90000 | 0 / 90000 |
| 300 x 300, 3 ft cells, origin (300000, 200000) | 599 / 90000 | 599 / 90000 |

So **anchoring at the DEM's origin is not the fix**: it saves the fixture and
the first grid by luck, but a far edge whose coordinate is off the lattice
still moves out (last row). The counts that leave depend on the origin as
well as the cell.

**The follow-up's scope is any DEM whose nodes are not fixed points of the
snap**, metre DEMs included: a node `x` moves when
`llround(x / 0.001) * 0.001`, in float64, is not `x`, and it moves by one
float step even when `x` is a whole number of millimetres, because 0.001 is
not exact in float64. `@reviewer` measured it through `rasputin mesh` at
`0f67a567` (design review round 3): 4 x 4 grids of 10 m cells keep 12 of 16
points with the first node at x = 5440353.066 and 9 of 16 at x = 3275842.3.
The emulation agrees: of 200 000 random whole-millimetre values below
10^7 m, 13.5 % move; of 200 000 random whole metres below 2^30 m, none; of
200 000 random multiples of 2^-3 m in each binade from 2^20 to 2^42 m, none,
and from 2^42 m on, some do (97 113 in the binade from 2^42 m).

Which inputs that reaches in practice, measured over every GeoTIFF under
`../rasputin_data` (each file's node rows and columns as `read_meta` gives
them, the same emulation): the projected DEMs (DTM10 UTM33, both sets of 254
tiles; the Austrian 1 m and 10 m files; the projected GLO30 and São Francisco
30 m files, 462 of them) have whole-metre origins and spacings, and no node
moves. The geographic files' degree nodes do move, but a geographic DEM is
always resampled (`--out-crs` is required), and the resampled canvas of
increment 15c has its origin at `col0 * h` and `-row0 * h` with `h` a whole
number of metres or 2^-1 to 2^-3 m (once this addendum's floor of 2^-3 m is built; the code at `0f67a567` still accepts 2^-10 m), so its nodes are multiples of 2^-3 m and
do not move below 2^42 m. So today the reach is: feet grids, and projected
DEMs with an origin or spacing that is not a multiple of 2^-3 m (decimetres,
centimetres, millimetres), none of which is in the data on disk. The real fix decides how the noder's lattice and a DEM's
lattice meet (keeping DEM nodes unsnapped, a spacing that divides the cell, or
clamping a sampled vertex onto the raster's edge); it touches the noder or
the sampler, so it is not Python-only and is outside this design. **It goes to
its own increment**, with the main session to propose it to Ola.

In 43, `@tester` splits the test:

- **`test_given_us_ft_the_mesh_is_in_metres`, on a fixture whose nodes are on
  the lattice**: tie point (6 499 987, 2 098 421) ftUS, scale 3937 x 3937
  ftUS (3937 ftUS is 1200 m exactly), pixel is area, EPSG:2229, no 4099, the
  same `HEIGHTS`. Measured through `rasputin mesh --z-unit us-ft` at
  `0f67a567`: exit 0, 16 points, x from 1981800.0 to 1985400.0, y from
  635400.0 to 639000.0 (all exact), the z values equal the float64 products
  cast to float32, `dem_file_crs` as in the test, no assumption line. Assert
  `points[:, 0].min() == 1981800.0` and `points[:, 1].max() == 639000.0` in
  place of the half-cell expressions, and keep every other assertion.
- **`test_a_feet_grid_off_the_millimetre_lattice_keeps_its_edge`**, the old
  3 ftUS fixture and assertions, marked
  `@pytest.mark.xfail(strict=True, reason=...)`, the reason naming the 1 mm
  snap anchored at zero and saying a later increment fixes it. Strict, so the
  fix turns it into a failure that makes someone remove the mark.

Size: no production line. No `@perf` run: nothing under
`include/terrain/refinement/` or `include/terrain/mesh/` changes.

### 9.4 The eighth failure

Not reproduced. Two full runs at `0f67a567` each gave the same 7 failures
and 7655 passes; the commit's message names 7 and leaves the eighth
unnamed, so it was intermittent; its cause was not found. `@tester` names any test outside 9.1 to 9.3 that fails
in the red-to-green run, with its output, before `@developer` is spawned.

### 9.5 Found here: the A4 floor below 2^-3 m is off the same lattice

The same emulation over A4's resampled grids (300 x 300 nodes, origin a
multiple of the spacing, at x = 500000, 4000000 and -1234567 plus a few
cells): spacings 2^-1, 2^-2 and 2^-3 m move no node, and multiples of
2^-3 m stay fixed below 2^42 m (section 9.3's sampling), so section 4.6 states
the bound as 2^42 m. From 2^-4 m (62.5 mm, not a whole millimetre) down,
67500 to 89996 of 90000 nodes move and, at these origins, 300 to 599 leave the
raster (how many leave depends on the origin); at 2^-10 m (0.98 mm, under the snap spacing) 4151 nodes
also merge (85849 distinct of 90000). Section 4.6 shipped 2^-10 as the floor,
so a source finer than about 9 cm would mesh wrongly. **`@developer`:** the
floor becomes 2^-3 m: in `src_python/tin_engine/target_grid.py`,
`default_spacing`'s `2.0**-10` becomes `2.0**-3`, and `_lattice_step` accepts
2^-k for k from 1 to 3 (its message says so). No net line. **`@tester`:** in
`tests/python/test_spacing_and_margin.py`, the accepted case `2.0**-10`
becomes `2.0**-3`, and `2.0**-4` is refused; `default_spacing` of a
geographic meta whose `delta_y` is 1e-7 degrees (about 1 cm) is 0.125.
Sections 4.6, 5 and test 13 are updated to match.

### 9.6 Ruling: push and enqueue

Ola, 2026-10-13, recorded by the main session: yes to pushing increment 43
and enqueuing it once it is review-approved. The main session holds this yes
under `.claude/REQUIRED-READING.md`'s approval rules (one push and one
enqueue; a yes recorded on disk is asked again after a context loss).

## Review

**Design review, round 1, 2026-10-13.** Range `97c10c4d..192d7245`, 0 counted lines, @reviewer: CHANGES REQUESTED. A named feet-CRS twin is 0.1 mm off, not exact (EPSG:2229's false easting is 2000000.0001016 m against EPSG:26945's 2000000 m; state the offset and that PROJ's equivalence decides the name); the height test figures don't say float32 or float64 (default float64, cast once); the order paragraph is stale (38, not 41, shares the call chain; merge master before the red step); no speed line; the noop test of the widened `epsg_matches` is unnamed; section 4.8 overclaims unchanged bytes (the run record changes for cached sources). Suggestions: refuse Mercator variant B's key 3081 unless 0; add `fetch/plan.py`'s `read_page` to the `z_unit` list; a test where a converted height equals NoData; say key 4096 is now read; name the module of `source_region`.

**Design review, round 2, 2026-10-13.** Range `3406a453..93423d2d`, 0 counted lines, @reviewer: APPROVED. All six round-1 blockers and all five suggestions resolved; figures and probes re-run (pyproj 3.8.0, PROJ 9.8.1). Suggestions: pin whether the outline snap limit, which shares `MARGIN` (`src_python/tin_engine/feature_input.py@3406a453:634`), stays at 100 m when `MARGIN` follows the cell size; the 2263 to 32118 offset is -1.7e-10 m in easting by the probe, not 5.8e-11 m (section 4.4) or 6e-11 m (constants table).

**Design review, round 3, 2026-10-13.** Range `0f67a567..92e0c3c5`, 0 counted lines, @reviewer: CHANGES REQUESTED. Section 9.3's "no node moves" is false for metre DEMs with whole-millimetre origins (first node at x = 5440353.066 keeps 12 of 16 points, at x = 3275842.3 keeps 9 of 16, 10 m cells, through `rasputin mesh`), so the follow-up covers any DEM whose nodes are not exactly whole 0.001 m steps in float64; sections 4.6 and 9.5's 2^50 m bound becomes 2^42 m. Fingerprints, the 1200 m grid, 599 and 4151, and the narrowed pattern confirmed. Suggestions: give the 1 ftUS row's origin ((300000, 200000) ftUS reproduces 300); "300 to 599 leave" depends on the origin; add `\bconstraint-feet\b` to the narrowed pattern.

**Design review, round 4, 2026-10-13.** Range `7f524956..d1028d1a`, 0 counted lines, @reviewer: CHANGES REQUESTED. Round 3's blockers fixed; the new "no multiple of 2^-3 m moves below 2^43 m" is false (4552287362515.375 m moves; 97,742 of 200,000 from 2^42 to 2^43 m move), so 2^43 becomes 2^42 at the constants row and section 9's three statements. 13.6 %, whole metres, the 977 projected files and the pattern confirmed. Suggestion: the resampled-canvas spacings of 2^-1 to 2^-3 m hold once the floor change is built. Both applied as wording fixes by the main session.

**Code review, round 1, 2026-10-13.** Range `origin/master..3fc5b60f`, 300 counted lines, @reviewer: CHANGES REQUESTED. (1) 26 new history labels such as "(43, 4.3)" and "increment 43, 4.5" in docstrings and comments across 12 production files break the docstring rule merged with increment 42 (#246); the checker misses them because its scope is empty. (2) Section 4.8's "geographic and metre DEMs, mesh files byte for byte" is false for resampled sources finer than about 0.71 m (projected ones too) and for cells over 50 m when features are given. (3) The Status line is stale. (4) `tools/scratch_copy.py:19` cites `tests/python/test_io_geotiff.py:1271`, which this branch moved to 1275. Verified: the six quick-check meshes are unchanged, the kill record covers every target and the equivalent mutant holds, no leftover red-step code, check_citations exits 0. Suggestion: `CODES` one entry per line.
