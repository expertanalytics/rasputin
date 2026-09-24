# Increment 11: raster ingestion, the Python side

**Status.** Implemented, in review, not merged. Written by `@architect`
before `@tester` was spawned, per `docs/increments/README.md` step 1. The red
suite and the green implementation are on branch `increment11-raster`
(`git log --oneline -- src_python/tin_engine/io/ tests/python/test_io_geotiff.py`).
`@reviewer`'s first pass led to round 3 below. Its rulings are not in the
code or the suite yet. This file holds no code and no tests.

*Amended (round 3), the status line.* It said "Design" after the code was
green.

**Amended after the red suite.** `@tester`'s red commit
(`git log --oneline -1 -- tests/python/test_io_geotiff.py`) turned up seven
design problems. Each one is ruled in place below, marked *Amended (problem N)*.
Two of them went to the user, who chose A1 and B1; §14 records both.

**Amended again (round 2).** Updating the suite to match turned up two more
gaps and four readings `@tester` had to pin. They are ruled in place, marked
*Amended (round 2)*, and §12 lists the tests they change.

**Amended a third time (round 3).** One user ruling (refuse a compound CRS in
3072) and `@reviewer`'s design-side findings on the green code. They are ruled
in place, marked *Amended (round 3)*. §12 lists the tests they change. One
question went to the user: errors raised from inside tifffile or a codec.
The user chose option C on 2026-09-24; §14 records it.

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
3. **The reader's job is refusal.** Nineteen named refusals, section 5.
   *Amended (round 3):* eighteen before `refuses_compound_crs` (13a).
4. **Area-registered files are converted, not refused** — the project's only DEM
   fixture is one. An *absent* registration key is refused.
5. **Use `tifffile`'s own GeoKey decoding.** Do not port `GeoKeysInterpreter`.
6. **Only `ProjectedCSTypeGeoKey` can yield an accepted CRS**, resolved through
   `pyproj.CRS.from_epsg`. No proj4 reassembly, no free-text ellipsoid regex.
   `GeographicTypeGeoKey` is read only when 3072 is absent, and only so that
   the refusal names the real reason (§5, problem 3). *Amended (round 2):* a
   code reached through 2048 is always refused, whatever it resolves to (§5,
   refusal 13). *Amended (round 3), user ruling:* the 3072 CRS must also be
   two-dimensional. A compound CRS in 3072 is refused (§5, refusal 13a).
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
  first-party except `io/models.py`, and **never `_core`** — the same shape as
  `io/ply.py`, which is why increment 10 was cheap.

  *Amended (problem 2).* This section used to say the modules are "testable
  with no compiled extension in the process". That is false as the tree
  stands. `tin_engine/__init__.py` re-exports `Point2`, `Point3`, `cross` and
  `dot` from `_core`, so importing any submodule loads the extension first.
  Check it:

  ```
  $ .venv/bin/python -c "import sys, tin_engine.io.ply; print('tin_engine._core' in sys.modules)"
  True
  ```

  The claim that holds is about the source: neither module imports `_core`,
  and `test_module_never_imports_core` checks that by walking the AST. Whether
  to make the stronger, process-level claim true went to the user, who chose
  A1 (§14): the source-level claim is the one this increment makes.
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

### A tie point that is not at pixel (0, 0)

*Amended (problem 7).* The formulas above assume the tie point's raster
coordinates `(I, J)` are `(0, 0)`. That is true of every file GDAL writes, but
nothing required it. **Ruling: convert, do not refuse.** A tie point says that
raster point `(I, J)` maps to model point `(X, Y)`. With a north-up affine
placement, that fixes pixel `(0, 0)` exactly, so there is nothing to guess:

```
x_min = X - I * delta_x  [+ delta_x / 2 if area-registered]
y_max = Y + J * delta_y  [- delta_y / 2 if area-registered]
```

`I` is the **column** and `J` the row. This is a port, not a new rule: the
legacy did the same offset at `legacy/rasputin/reader.py:392-393` (prior art
§3.2), under variable names that had the axes the wrong way round (prior art
§4.3). `I` and `J` must be finite, like every other tie-point value; see
refusal 1.

The tie point's `K` and `Z` are ignored. Refusal 4 already requires ScaleZ to
be zero, and with ScaleZ zero the file declares no vertical mapping for `Z` to
offset, so dropping them discards nothing.

*Amended (round 2), reading (a) confirmed.* "Ignored" covers non-finite
values too. A NaN or infinite `K` or `Z` is not refused, because the reader
never uses them, and refusing a value that is not used would be refusing on
something the file does not assert. Refusal 1's finiteness rule covers exactly
`I`, `J`, `X` and `Y`, the four values the placement uses. Measured, tifffile
2026.9.20: a tie point written as `(0, 0, nan, X, Y, inf)` reads back with
the `nan` and `inf` intact, so a micro-TIFF can carry this case.

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
   `ModelPixelScaleTag` (33550) absent, or any of the tie point's `I`, `J`,
   `X`, `Y` not finite (*amended, problem 7*: a NaN origin is as unusable as a
   missing one, and nothing else caught it). The legacy defaulted them to zeros and
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

   *Amended (round 3), `@reviewer` blocking finding 1: a malformed tie point
   escaped refusals 1 to 3.* `decode_dem` calls `tif.geotiff_metadata` before
   it reads the georeferencing tags. `geotiff_metadata` reshapes 33922 into
   rows of 6 and 34264 into 4×4, so a wrong length raises tifffile's bare
   `ValueError` before any refusal can run. Measured, tifffile 2026.9.20, on
   micro-TIFFs through `decode_dem`:

   - 33922 with 1, 3, 5, 7 or 9 values:
     `ValueError: cannot reshape array of size N into shape (6)`. With 0 or 12
     values the refusal fires as designed.
   - 34264 with 1, 3, 12 or 17 values:
     `ValueError: cannot reshape array of size N into shape (4,4)`. So
     refusal 2 was delivered only for exactly 16 values.
   - 33550 with 1 value: `TypeError: 'float' object is not iterable`. This
     one is in the reader, not tifffile: a DOUBLE tag of count 1 reads back as
     a bare `float`, not a tuple (33922 does the same).

   **Ruling: the three tags are checked from `page.tags` before
   `geotiff_metadata` is read.** The order in `decode_dem` becomes:

   1. `_single_page`, then `_check_page` (refusals 7 to 11). Neither reads
      GeoKeys.
   2. The tag checks, reading `page.tags` only, in this order:
      - 34264 present, with any length or value: refusal 2.
      - 33922 or 33550 absent: refusal 1.
      - 33922 with more than 6 values: refusal 3.
      - 33922 with fewer than 6 values, or 33550 with any count but 3:
        refusal 1.
      - Then the existing finiteness, sign and ScaleZ checks.

      A scalar tag value is a sequence of length 1.
   3. `tif.geotiff_metadata`, and only then the steps that need GeoKeys:
      registration (refusal 6), then CRS and units (12 to 14).

   So `_placement` splits in two: the tag part moves before
   `geotiff_metadata`, and the registration part stays after it. No new
   refusal name is needed. The messages change, because the counts are now
   the diagnosis:

   - refusal 3: "`ModelTiepointTag (33922)` has N values; exactly 6 (one tie
     point) are read". The old text said "a GCP list", which is false for 7
     or 9 values.
   - refusal 1, short tie point: "`ModelTiepointTag (33922)` has N values;
     need 6". Short scale: "`ModelPixelScaleTag (33550)` has N values; need
     3".

   Moving a 2-value 33550 from refusal 5 to refusal 1 changes no existing
   test: every `test_refuses_nonpositive_pixel_scale` case has 3 values.

   Measured, the same run: a corrupt `GeoKeyDirectoryTag` (34735), or one
   that points into a 34736 that is absent, does **not** raise.
   `geotiff_metadata` logs an error and returns the keys it could read, and
   the reader then refuses on the missing key (1025 or 3072) with a true
   message. So 34735 needs no pre-check.
6. `refuses_unknown_raster_type` — `GTRasterTypeGeoKey` (1025) absent or 32767.
   Section 4.
7. `refuses_degenerate_shape` — fewer than two rows or two columns. `bilinear`
   returns `nullopt` for a raster too small to hold a 2×2 neighbourhood, so
   such a tile can be decoded and then never sampled. Refuse at decode.
8. `refuses_ambiguous_pages` — more than one page where the extra pages are not
   flagged reduced-resolution. Page 0 is used; being silently one of several
   full-resolution images is not acceptable. *Amended (round 2), reading (c)
   confirmed:* page indices in the message are 0-based, to match "page 0"
   here and tifffile's `pages[i]`. The page after the one that is read is
   page 1.

   *Amended (round 3): an extra page that is a `TiffFrame` is refused.* The
   green code tests `isinstance(extra, TiffPage) and not extra.is_reduced`,
   so a `TiffFrame` passes without a check. A `TiffFrame` is a page tifffile
   chose not to parse in full. It has no `subfiletype` and no `is_reduced`
   (measured), so the reader cannot know whether it is reduced-resolution.
   Skipping it is the silent case this refusal exists to stop. **Ruling:
   refuse any extra page that is not a `TiffPage`.** The message names the
   page index and the page count. It says that `NewSubfileType` (254) could
   not be read, because tifffile returned a frame. It does not print a value.

   When does this happen? `TiffFile(stream)` with no private arguments makes
   frames only for LSM, NDPI and ScanImage files (`tifffile/tifffile.py`,
   `TiffFile.__init__`, the `_lsm_load_pages`, `_ndpi_load_pages` and
   `_load_virtual_frames` branches). All three are microscopy formats, never
   a DEM. Measured on a two-page micro-TIFF: both pages are `TiffPage`, and
   still are after `geotiff_metadata`. With the private `_useframes=True`,
   page 1 is a `TiffFrame`. So refusing costs no real input. The alternative
   was to force a full parse (`tif.pages.useframes = False` before the loop).
   It is rejected because it reaches into tifffile's page cache for a case
   no DEM produces.
9. `refuses_multi_sample` — `SamplesPerPixel` is not 1. A DEM has one band.
10. `refuses_unsupported_dtype` — anything outside the promotion table in
    section 7.
11. `refuses_missing_codec` — the page's compression, **or its predictor**,
    needs `imagecodecs`. The message names the tag (`Compression` (259) or
    `Predictor` (317)), the scheme, and the `codecs` extra. A bare exception
    from inside `tifffile` is not a diagnostic; see section 8. *Amended
    (problem 1):* the check is a capability probe, run before any decode, and
    not a list of scheme names. §8 gives the probe, and why the predictor is
    included. *Amended (round 2), reading (b) confirmed:* the scheme name is
    tifffile's enum member name, `COMPRESSION(v).name` or
    `PREDICTOR(v).name`. For predictor 3 that is `FLOATINGPOINT` (measured,
    tifffile 2026.9.20), and a case-insensitive match on "floating" finds
    it. A value that tifffile's enum does not know has no name; the message
    then gives the number alone, which it names anyway.

    *Amended (round 3): the message wording.* `@reviewer` found three faults
    in the green message (`geotiff.py`, the `value not in known` branch).
    The ruling for each:

    - **Private API.** It looks names up with `enum._value2member_map_`.
      **Ruling:** use `value in enum`. That is public on Python 3.12 and
      later, which is the project's floor. Measured, Python 3.14.7:
      `5 in tifffile.COMPRESSION` is `True`, and `60000 in ...` is `False`.
      `COMPRESSION(60000)` raises `ValueError`, so it must not be the test.
    - **"(unknown)".** For a number the enum does not know, it prints
      `= 60000 (unknown)`. **Ruling:** print `Compression (259) = 60000`,
      with no brackets, as this refusal already said.
    - **Install advice that may be false.** It tells the user to install
      `codecs` whatever the scheme is. Measured, tifffile 2026.9.20 with
      imagecodecs 2026.8.16 installed in a scratch target: 28 `COMPRESSION`
      members are still not in `DECOMPRESSORS`, e.g. `THUNDERSCAN` (32809),
      `PIXARLOG`, `SGILOG` and `JBIG`. Every `PREDICTOR` member is. So with
      the extra present, the green code advises installing what is already
      installed. **Ruling:** the advice depends on whether `imagecodecs` can
      be imported (`importlib.util.find_spec("imagecodecs") is not None`).
      No import is needed.
      - **Absent:** "…cannot be decoded as installed; the `codecs` extra
        (imagecodecs) may decode it: pip install 'rasputin[codecs]'". The
        word is "may", because the reader cannot know without the extra.
        LZW and the floating-point predictor are the known cases, and §8's
        probe checks any other.
      - **Present:** "…cannot be decoded by tifffile, even with imagecodecs
        installed". No install advice.

    The tag, its number and value, and the scheme name when there is one
    stay as they are. `test_refuses_missing_codec`'s assertions still hold,
    because it runs only without the extra.

**CRS**

12. `refuses_missing_crs` — no `ProjectedCSTypeGeoKey` (3072), or its value is
    32767, or it is not a resolvable EPSG code. Ruling 6: there is no fallback
    path. The legacy had two, one of them a regex over human prose.
13. `refuses_geographic_crs` — the constructed CRS has `is_projected` false.
    Section 3.

    *Amended (problem 3).* As first written, CRS came from 3072 alone. A
    normally encoded geographic file has `GTModelTypeGeoKey` (1024) = 2,
    `GeographicTypeGeoKey` (2048) = 4326, and no 3072. Such a file therefore
    hit refusal 12 ("no CRS") and never refusal 13. The message was false: the
    file does declare a CRS, just a geographic one. So refusal 13 could only
    fire for a file that puts a geographic code in the projected key, which is
    the rare case.

    **Ruling.** CRS resolution tries 3072 first. If 3072 is absent, it tries
    2048. Whichever code it uses is resolved with `from_epsg`, and the result
    goes through the same `is_projected` test. So the geographic file is
    refused as geographic, and the message names 2048, `GeographicTypeGeoKey`
    and its value. Ruling 7 still holds: the refusal comes from the
    constructed CRS, not from which keys are present. ~~Ruling 6 still holds
    too. 2048 is always a geographic CRS, so it can never produce an accepted
    tile.~~ (Struck in round 2: that was a claim about the file, not a rule
    the reader enforced. See below.) It is not a fallback path. Refusal 12
    now means that neither key yields a resolvable code. That covers 3072
    absent with no 2048, and 3072 set to 32767 or unresolvable. The message
    names 3072, plus 2048 if it was consulted. `GTModelTypeGeoKey` is still
    not read.

    *Amended (round 2).* The problem-3 rule put a 2048 code through
    `is_projected` and accepted it if that passed. So a file with no 3072 and
    2048 = 25833 was accepted, and ruling 6 ("only 3072 can yield an accepted
    CRS") was false. Nothing in the reader stopped it. Measured, tifffile
    2026.9.20: such a micro-TIFF reads back with `GeographicTypeGeoKey` =
    `25833`, a plain int because 25833 is not in tifffile's `GCS` enum, and
    `pyproj.CRS.from_epsg(25833).is_projected` is `True`.

    **Ruling.** A code reached through 2048 is never accepted. Which refusal
    fires depends on the CRS it resolves to, so the message stays true:

    - it resolves and `is_geographic` is true: `refuses_geographic_crs`,
      naming 2048, `GeographicTypeGeoKey` and the code. This is the case
      problem 3 was about.
    - it resolves to anything else (projected, geocentric, vertical,
      compound): `refuses_missing_crs`. The message names 3072 as absent,
      and 2048 with its code and the CRS's `type_name`, e.g. "Projected
      CRS". It says that a projected CRS belongs in 3072.
    - it does not resolve, or is 32767: `refuses_missing_crs`, as before.

    Measured, pyproj 3.8.0: `is_geographic` is `True` for 4326, 4258 and
    4979, and `False` for 25833 (projected), 4978 (geocentric), 5972
    (compound) and 3855 (vertical).

    The alternatives, and why not:

    - *Refuse every 2048 code as geographic, whatever it is* (the main
      session's suggestion). It enforces ruling 6 just as well, and costs
      one line less. But for 2048 = 25833 the message would say "geographic"
      about a CRS that pyproj says is projected. That is the kind of false
      message problem 3 was raised to remove, and it breaks ruling 7: the
      refusal would come from which key was present, not from the CRS.
    - *Accept it.* A lenient reader could say the file clearly means
      25833. But the GeoTIFF spec puts a projected code in 3072, GDAL never
      writes one in 2048, and a file that does has broken its own key
      directory. Ruling 3 is that the reader's job is refusal. Guessing what
      a malformed file meant is the legacy's failure mode.

    These do not set two sound principles against each other. The chosen
    rule meets both ruling 6 and ruling 7, so it is ruled here and not sent
    to the user.

    *Also round 2, found while checking the above.* Refusal 13 fires when
    `is_projected` is false. Through 3072 that includes a geocentric or
    vertical code, and a message that says "geographic" would then be false
    for the same reason. So the refusal 13 message names the constructed
    CRS's `type_name` (e.g. "Geographic 2D CRS", "Geocentric CRS") instead of
    asserting "geographic". The test name stays as it is. No test pins the
    word, so no test changes.

    *Round 2, reading (d) confirmed.* When 3072 is present, whatever its
    value, 2048 is not consulted. For 3072 = 32767 the message must name
    "3072" and "32767". It need not mention 2048, and must not claim 2048
    was consulted. It may mention that 2048 is present; the suite does not
    pin that, and it should not.

13a. `refuses_compound_crs`: *added in round 3, by user ruling.* The 3072
    CRS passes `is_projected` but is not a 2-D horizontal CRS.

    **Why.** 3072 = 5972 (ETRS89-NOR / UTM 32N + NN2000 height) was
    accepted, with `epsg=5972`. pyproj reports a compound CRS as projected
    when its horizontal part is projected, so refusal 13 does not fire. All
    three axes are metres, so refusal 14 does not fire either. GeoTIFF puts
    the vertical part in `VerticalGeoKey` (4096), not in 3072. The user
    ruled that it must be refused.

    **Detection.** Refuse when `crs.is_compound` **or**
    `len(crs.axis_info) != 2`. Check it right after the `is_projected` test
    and before the unit checks. Measured, pyproj 3.8.0 / PROJ 9.8.1, every
    EPSG code that `from_epsg` resolves and that has `is_projected` true:

    | `type_name` | axes | count |
    |---|---|---|
    | Projected CRS | 2 | 5362 |
    | Compound CRS | 3 | 320 |
    | Projected CRS | 3 | 1 |

    The one 3-axis projected CRS is 9895, "LUREF / Luxembourg TM (3D)",
    with an ellipsoidal-height axis. It is not compound, so `is_compound`
    alone misses it, and `type_name` alone ("Projected CRS") misses it too.
    The axis count catches both shapes. `is_compound` is kept anyway: it is
    what the user ruled on, and it is what gives the message a sub-CRS to
    name. The same scan found nothing else that passes `is_projected`.

    It is one refusal, not two. The 3-D projected CRS has the same defect: a
    vertical axis in the key meant for a horizontal CRS. It is a compound
    CRS in all but its PROJ encoding.

    **Order against refusal 13.** `is_projected` is tested first, so nothing
    that refusal 13 already covers moves. 4978 (geocentric, 3 axes) and 9707
    (WGS 84 + EGM96 height, a compound with a geographic horizontal part,
    `is_projected` false) both stay refusal 13, and its message already
    names "Compound CRS" for 9707. If this test ran first, it would take
    over the geocentric case that `test_refuses_geographic_crs_names_the_crs_type`
    pins.

    **Message.** It names:
    - `ProjectedCSTypeGeoKey (3072)` and the code;
    - the constructed CRS's `type_name`;
    - the axis count;
    - the horizontal sub-CRS's code, if there is one: `sub_crs_list[0].to_epsg()`
      when `sub_crs_list` is not empty and `to_epsg()` is not `None`;
    - that the horizontal code belongs in 3072 and the vertical one in
      `VerticalGeoKey (4096)`.

    Example: "ProjectedCSTypeGeoKey (3072) = 5972 is a Compound CRS with 3
    axes; horizontal part EPSG 11022. Put the horizontal CRS in 3072 and the
    vertical CRS in VerticalGeoKey (4096)".

    Measured: the horizontal code for 5972 is **11022**, "ETRS89-NOR
    [EUREF89] / UTM zone 32N". It is not 25832. So the message gives what
    pyproj gives and does not promise a particular code. `to_epsg()` can
    return `None`, which is why "if there is one" is part of the rule. For
    9895, `sub_crs_list` is empty and the message names "Projected CRS" and
    "3 axes".

    Through 2048, a compound code is already refused as missing CRS, naming
    "Compound CRS" (round 2). Nothing changes there.
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

*Amended (problem 4).* The report said that four refusals have no tag to
name: degenerate shape, ambiguous pages, unsupported dtype and missing codec.
Only half of that is true. Each of the four comes from a baseline TIFF tag,
and the rule stands for all of them. **Ruling: every refusal names the tag
that its decision reads, by number and by name, plus the file's value.** Where
a tag's raw value means little on its own, the message adds the value derived
from it:

| refusal | names | and adds |
|---|---|---|
| `refuses_degenerate_shape` | `ImageLength` (257) and/or `ImageWidth` (256), with its value | "need at least 2" |
| `refuses_ambiguous_pages` | `NewSubfileType` (254) of the first offending page, with its value (0 when absent) | that page's index and the page count |
| `refuses_unsupported_dtype` | `SampleFormat` (339) and `BitsPerSample` (258), with their values | the numpy dtype, e.g. `int64` |
| `refuses_missing_codec` | `Compression` (259) or `Predictor` (317), with its value | the scheme name and the `codecs` extra |

Measured, tifffile 2026.9.20: an `int64` micro-TIFF reads back as
`bitspersample 64, sampleformat 2`, and a `bool` one as `1, 1`. The numpy
dtype is kept in the message because `SampleFormat 2, BitsPerSample 64` is
not how anyone would search for the problem. Both columns are required.

---

## 6. Ruling 9: NoData

This is the largest genuinely new decision in the increment. There is no prior
art: `GDAL_NODATA` (42113) appears nowhere in the legacy and a void in a DEM
became an elevation (prior art §1.3). `sample.hpp` already implements the
semantics — any of the four bilinear corners NoData returns `nullopt` — and
`raster.hpp`'s `is_nodata` also catches NaN unconditionally via `v != v`.
Discovery is all that is missing.

**Tag present.** Tag 42113 is ASCII. The fixture's is the string `'-32767'`.
**Parse the tag's text, never `page.nodata`.** Measured, tifffile 2026.9.20:
for the text `'void'` or `'12abc'`, `page.nodata` logs a warning and returns
`0`. That turns a malformed tag into a sentinel which deletes every
sea-level cell. Read `page.tags[42113].value` and parse it with `float()`,
refusing text that contains `_` (Python's `float()` accepts `'1_000'`; no TIFF
writer means that). Then:

- text that `float()` rejects: `refuses_unparseable_nodata`.
- everything else goes through **one representability check** (below). If it
  passes and the value is NaN, the result is `nodata = None`. `is_nodata`
  already catches NaN, and a NaN sentinel compared with `==` would match
  nothing.

*Amended (problem 5).* The single check that follows replaces the earlier
"round trip through the promoted dtype". It answers all four of the
questions `@tester` raised.

**The representability check, against the *file* dtype.** A sentinel is
accepted only if a cell of the file's own dtype can hold exactly that value:

- **integer file**: the value is finite, integral, and inside that dtype's
  `[min, max]`;
- **float file**: the value is NaN, or finite and survives an exact round
  trip through that dtype.

Anything else is `refuses_nodata_not_representable`, with the message naming
the file dtype.

Why the file dtype and not the promoted one (5c). The check exists because
`sample.hpp` compares `v == *nodata_`, so a sentinel that no cell can equal
switches NoData off without a word. Cells can only hold values of the file's
dtype. The promoted dtype is wider, so it accepts sentinels that can never
match: `0.5` or `-9999` on a `uint8` file pass a float32 round trip and still
match no cell. Every row of the promotion table in §7 is exact by
construction, so a value that is representable in the file dtype stays equal
after promotion. A NaN on an integer file is refused for the same reason: no
integer cell is NaN.

Why no infinities (5d). `float('1e400')` is `inf`. The text names a finite
number that no IEEE type can hold, but the float is indistinguishable from a
file that meant infinity. The earlier rule accepted it, because inf survives a
round trip. The ruling is that the check requires *finite or NaN*, so
`'1e400'`, `'inf'` and `'-inf'` are all refused. A file whose voids really are
±inf cannot be read. That trade is chosen on purpose: the refusal is loud,
no such DEM product is known, and the rule can be loosened when one appears.

**The caller's sentinel passes the same check (5b).** `v == *nodata_` does
not care where the sentinel came from, so the check cannot care either. A
caller sentinel that fails it is refused under the same name,
`refuses_nodata_not_representable`. The message names the `nodata=`
argument, not tag 42113.

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

*Amended (round 3): a `bool` is not a sentinel.* Measured on the green code:
`decode_dem(s, nodata=True)` gives `nodata == 1.0`, and `nodata=False`
gives `0.0`. `np.True_` behaves the same way. The type checker does not
catch it, because `bool` is an `int` and an `int` is accepted where a
`float` is expected. `False` would delete every sea-level cell, which is
the §6 failure written as a typo. **Ruling:** `decode_dem` raises
`TypeError` when `nodata` is a `bool` or a `numpy.bool_`. The message
names the `nodata=` argument and the type. It is `TypeError` and not
`GeoTiffError`, for the reason in §7 (B1): it is a programming error in the
caller, not a fact about a file. Any other real number is still converted
with `float()`, as before.

`RasterMeta` records whether the sentinel came from the tag, from the caller,
or is absent. An absent sentinel is a fact about the tile, and a later stage
that cares can ask.

**`nodata_source` records who made the NoData declaration, not whether a
number crosses the boundary (5a).**

| tag | caller | `nodata` | `nodata_source` |
|---|---|---|---|
| absent | absent | `None` | `"absent"` |
| absent | `v` | `v`, or `None` if NaN | `"caller"` |
| `t` | absent | `t`, or `None` if NaN | `"tag"` |
| `t` | equal to `t` | `t`, or `None` if NaN | `"tag"` |
| `t` | differs from `t` | refused | `refuses_contradictory_nodata_override` |

A `'nan'` tag is a declaration: the file says its voids are NaN. Recording it
as `"absent"` would tell a later stage that the file said nothing, which is
false. When the caller agrees with the tag, the answer is `"tag"`, because
the file would have given the same result without the caller. "Equal" means
the two parsed values are equal, with NaN treated as equal to NaN. The
contradiction check runs after both values have passed the representability
check, so a message never compares a valid value with an invalid one.

The model validator enforces one implication: `nodata_source == "absent"`
means `nodata is None`. Its converse does not hold, because of the NaN rows.

*Amended (round 2): a NaN sentinel becomes `None` whoever declared it.* The
rule above was stated for the tag only, and the table said `v` for a caller
with no tag. The reason for the rule never depended on the source: a NaN
sentinel matches no cell under `==`, `is_nodata` already catches NaN, and so
the two tiles behave the same. So:

- **Caller `nodata=nan`, no tag, float file:** `nodata = None`,
  `nodata_source = "caller"`. The caller made a declaration, and the file
  did not. This is the table's second row, now with the NaN clause.
- **Caller `nan`, integer file:** refused as not representable, as before.
  The NaN check runs after the representability check, as for the tag.
- **Caller `nan`, tag `nan`:** equal, so `nodata = None` and `"tag"`, the
  fourth row. This agrees with the row above: both give `None`, and only the
  source differs, because only the source should.
- **Caller `nan` with a finite tag, or a finite caller value with a `nan`
  tag:** not equal, so `refuses_contradictory_nodata_override`. NaN is equal
  to NaN here and to nothing else.

**Enforced at the type.** `RasterMeta.nodata` is declared
`float | None = Field(allow_inf_nan=False)`, so a `RasterMeta` holding a NaN
or an infinite sentinel cannot be built. The reader cannot forget the rule on
one path and keep it on another. It also keeps `RasterMeta` equality
working, which a NaN field breaks (`nan != nan`). Measured, pydantic 2.13.5:
`None` and `-9999.0` are accepted, and `nan` and `inf` both raise
`ValidationError` with `type=finite_number`. As in §7 (B1), that is a
programming error in whoever built the model, so it is Pydantic's error and
not `GeoTiffError`.

---

## 7. Types

All Pydantic V2, all frozen. `io/models.py`.

```
RasterMeta
    x_min, y_max, delta_x, delta_y : float        # node grid, metres
    cols, rows                     : StrictInt    # DemTile checks == array.shape (B1, §14)
    epsg                           : int          # projected, metre
    nodata                         : float | None # finite; allow_inf_nan=False (§6, round 2)
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

**Dimensions come from `array.shape`.** `project_structure.md` already rules
this and prior art §5.7 is the incident: the legacy built a shape out of
`numpy.float64` values and its guarding assertion compared `(12.0, 16.0)` to
`(12, 16)` and passed.

*Amended (problem 6).* As first written, `cols` and `rows` were to be
"derived in the model validator, never accepted as input". They cannot be
derived by `RasterMeta`'s own validator, because the array lives on `DemTile`.
There are two sound readings, and they conflict. The user chose B1 (§14,
choice B):

- **Chosen (B1): `RasterMeta` keeps `rows` and `cols` as `StrictInt`.
  `DemTile`'s `model_validator(mode="after")` refuses a tile where
  `(meta.rows, meta.cols) != array.shape`.** `StrictInt` rejects `12.0`, which
  kills the legacy incident at the type. The `DemTile` check makes a mismatch
  impossible to construct. So "comes from the shape" is enforced, even though
  it is no longer literally "derived". The reason to prefer this:
  `RasterMeta` then describes the whole node grid on its own, including the
  far corner, which needs `cols` and `rows`. The mosaic increment's footprint
  walk (§10) has to read many tiles' extents from their headers without
  decoding 100 MB of pixels each. With B1 it can build a `RasterMeta` from
  `ImageWidth` and `ImageLength` and has nothing to reshape.
- **Alternative (B2): remove `rows` and `cols` from `RasterMeta`, and make
  them read-only properties on `DemTile` that return `array.shape`.** One
  source of truth, and no validator at all. The cost: `RasterMeta` stops
  being a complete grid description, and the mosaic increment will either add
  the fields back or invent a second metadata type.

The error for a B1 mismatch is Pydantic's `ValidationError`, not
`GeoTiffError`. It is a programming error in whoever built the model. It does
not describe a file.

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

*Amended (round 3): `DemTile` copies the array.* The green validator returns
`np.ascontiguousarray(value).view()` with the view's flag cleared. For an
input that is already C-contiguous, that is a view of the **caller's**
buffer. Measured, numpy 2.5.3: `np.shares_memory(a, tile.array)` is `True`,
and after `a[0, 0] = 42` the tile reads `42.0`. The paragraph above says
this must not happen. There is a second hole. The view's base is writeable,
so `tile.array.flags.writeable = True` succeeds, and the tile hands out a
writeable handle after all.

**Ruling: copy.** The validator makes an owned C-contiguous copy
(`np.array(value, order="C", copy=True)`), clears the copy's writeable
flag, and stores a view of the copy. Measured: setting `writeable = True`
on such a view raises
`ValueError: cannot set WRITEABLE flag to True of this array`. The other
choice was to keep the view and write down that the tile shares memory. That
is rejected, because the concurrency rule above is the reason the array is
read-only at all, and a shared buffer breaks that rule where no one can see
it.

**Cost.** `decode_dem` holds two copies of the pixels for a moment: the
array from `page.asarray()`, and the tile's copy. For the 102 MB fixture the
peak is about 204 MB. That is accepted for one tile. `decode_dem` must not
avoid the cost with `model_construct`, because that skips every validator,
including B1's shape check. Python cannot stop someone who reaches for
`tile.array.base` and sets its flag. The rule is only that the tile never
*hands out* a writeable handle.

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
dependency". Measured, tifffile 2026.9.20 (*amended, round 3*: this said
2026.9.15, and 2026.9.20 is what is installed), this tree's `.venv`, no
`imagecodecs`, **measuring writes** (see the amendment below):

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
repository's only DEM cannot be read without the `codecs` extra.

*Amended (problem 1).* **The measurement above is of `imwrite`, and the reader
only reads.** Reading and writing do not match. Measured, tifffile
2026.9.20, no `imagecodecs`: a real PackBits strip, spliced into a micro-TIFF
by hand, **reads** correctly through tifffile's built-in
`tifffile/_imagecodecs.py:packbits_decode`. LZW still fails on read, with
`ValueError: <COMPRESSION.LZW: 5> requires the 'imagecodecs' package`. The
read failure is a `ValueError`, not the `KeyError` shown above. So "PackBits
needs the extra" was wrong for reading, and a PackBits codec refusal could
never fire.

Which schemes need the extra depends on the installed tifffile, and on the
Python under it. `50000` (ZSTD) resolves without the extra here only because
tifffile's built-in module imports `compression.zstd`, which is new in Python
3.14 (`tifffile/_imagecodecs.py`, the `from compression import zstd` lines).
So the rule is written as a probe, not as a list:

```
$ .venv/bin/python -c "
import tifffile
d, p = tifffile.TIFF.DECOMPRESSORS, tifffile.TIFF.PREDICTORS
print({c: c in d for c in (1, 5, 8, 32773, 50000)}, {c: c in p for c in (1, 2, 3)})"
{1: True, 5: False, 8: True, 32773: True, 50000: True} {1: True, 2: True, 3: False}
```

**Ruling.** Before `asarray`, the reader checks
`page.compression in tifffile.TIFF.DECOMPRESSORS` and
`page.predictor in tifffile.TIFF.PREDICTORS`. If either is false, it raises
`refuses_missing_codec`. Membership tries to resolve the codec and returns
`False` instead of raising (`CompressionCodec.__contains__` and
`PredictorCodec.__contains__` in `tifffile/tifffile.py`), so this is
tifffile's own answer and not a copy of it. The reader never matches
tifffile's exception text.

**The predictor is new here.** `Predictor` (317) = 3, the floating-point
predictor, needs `imagecodecs` exactly as LZW does. Measured by patching a
Deflate micro-TIFF's predictor to 3:
`ValueError: <PREDICTOR.FLOATINGPOINT: 3> requires the 'imagecodecs' package`.
GDAL writes `PREDICTOR=3` for float DEMs, so this is not an edge case.

~~Not verified: that `TIFF.DECOMPRESSORS.__contains__` behaves the same across
the whole range that `pyproject.toml` pins (`tifffile>=2024.1`).~~

*Amended (round 3): the floor is checked.* `@reviewer` ran the suite at
tifffile 2024.1.30, and it was re-run for this amendment the same way. The
old tifffile goes into a scratch target, put ahead of the venv's copy on the
path:

```
$ uv pip install --python .venv/bin/python --target $SCRATCH/tf 'tifffile==2024.1.30' --no-deps
$ PYTHONPATH=$SCRATCH/tf .venv/bin/python -m pytest -q tests/python/test_io_geotiff.py tests/python/test_geotiff_fixtures.py
183 passed, 3 skipped, 143 warnings
```

That is the same count as at 2026.9.20 (183 passed, 3 skipped). Check that
the run used the old copy:
`PYTHONPATH=$SCRATCH/tf .venv/bin/python -c "import tifffile; print(tifffile.__version__)"`
prints `2024.1.30`. This is the lowest release the pin allows. The pin stays.
It is one data point at the low end, not a check of every release in the
range.

Consequences, both this increment's PR to land, per principle C3:

- Correct the `pyproject.toml` comment. Name no scheme list that rots. Say
  that without `imagecodecs`, tifffile reads uncompressed and Deflate plus
  whatever its built-ins cover; that LZW (the Kartverket fixture) and the
  floating-point predictor need the extra; and give the one-line probe above
  as the way to check any other scheme.
- `@tester` builds micro-TIFF fixtures with `compression=None` or `'deflate'`,
  so the suite runs with no optional dependency. Any test that reads the real
  Kartverket fixture is marked `skipif` on `imagecodecs` being importable.
  `refuses_missing_codec` is testable in both directions with **LZW**, not
  PackBits: refused with a useful message when the extra is absent, and read
  when it is present. The suite already does this.

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
  `(799745.0, 7950255.0)` at a 10 m scale. ~~Those two numbers are the only
  check that ruling 4 is implemented and not merely documented.~~

  *Amended (round 3).* That sentence was false twice:
  - `TestPlacement.test_area_registered_nodes_are_shifted_inward_half_a_cell`
    checks the half-cell shift on a micro-TIFF. `delta_x` is not equal to
    `delta_y` there, so a swapped shift cannot pass. It runs everywhere.
  - The two Kartverket tests never run in CI. `main.yaml` installs `.[dev]`
    only (`grep -n 'pip install' .github/workflows/main.yaml`), so
    `imagecodecs` is absent and `needs_codecs` skips them.

  The Kartverket tests check that the micro-TIFF result also holds on the one
  real product. They are a local check, not the check of ruling 4.
  `test_kartverket_fixture_is_shifted_half_a_cell`'s docstring repeats the
  false sentence, so `@tester` corrects it (§12, round 3). Whether CI should
  install `codecs` is not ruled here. Leaving it out is what proves the
  suite runs without optional dependencies.
- **The `always_xy` test must be shown failing** against a planted violation.
  Section 9, principle A3.
- **The promotion table is a parametrised test**, one case per row, asserting
  the resulting `array.dtype` — not that it is "a float".

### Tests that round 2 changes (`tests/python/test_io_geotiff.py`)

Added or extended. No existing assertion is reversed.

- `test_refuses_missing_crs`: add a case `absent_geographic_projected`,
  `{PROJECTED_CS_TYPE: None, GEOGRAPHIC_TYPE: EPSG_UTM33}`, naming
  `(*_PROJECTED, *_GEOGRAPHIC, "25833")`. This is the case that was
  accepted before. Optionally add `absent_geographic_geocentric` with 4978,
  which takes the same branch. `REFUSALS` in `geotiff_fixtures.py` needs no
  new entry, because the refusal is not new.
- `TestNoData.test_caller_nan_without_tag_yields_no_sentinel` (new): a
  float32 micro-TIFF with no tag, `nodata=math.nan`, gives `nodata is None`
  and `nodata_source == "caller"`.
- `test_refuses_contradictory_nodata_override`: add two parametrised cases,
  a NaN caller with tag `"-32767"`, and `nodata=-32767.0` with tag `"nan"`.
  Together with `test_caller_nan_agreeing_with_nan_tag_is_accepted` they pin
  "NaN equals NaN and nothing else" from both sides.
- `TestShapeAgreement` or a sibling class,
  `test_meta_refuses_non_finite_nodata` (new), parametrised over `nan` and
  `inf`: `RasterMeta.model_validate` with that `nodata` raises Pydantic's
  `ValidationError`.
- `TestPlacement.test_tie_point_k_and_z_are_ignored`: parametrise, adding a
  case with `K = nan` and `Z = inf`. Reading (a).

Unchanged, readings confirmed: `test_refuses_missing_codec` (b),
`test_refuses_ambiguous_pages` (c), and the `user_defined` and
`user_defined_beside_geographic` cases of `test_refuses_missing_crs` (d).

### Tests that round 3 changes

Files: `tests/python/test_io_geotiff.py` and `geotiff_fixtures.py`. No
existing assertion is reversed. Each item below fails against the green code
as it stands, unless it says otherwise.

**Refusal 13a, compound CRS (user ruling).**
- `test_refuses_compound_crs` (new), parametrised:
  - `compound_projected`: 3072 = 5972. Names `3072`,
    `ProjectedCSTypeGeoKey`, `5972`, `Compound CRS`, `11022` and `4096`.
  - `projected_3d`: 3072 = 9895. Names `3072`, `9895`, `Projected CRS`
    and `3`.
- `geotiff_fixtures.py`: add `EPSG_COMPOUND = 5972`, and a `REFUSALS`
  entry `refuses_compound_crs` whose defect check reads 3072 = 5972 back
  with tifffile alone.
- `test_refuses_geographic_crs_names_the_crs_type`: add a case
  `compound_geographic`, 9707 → `Compound CRS`. This pins the order: it
  must stay refusal 13. It passes today, and it guards the ordering.

**Refusals 1 to 3, checks moved before `geotiff_metadata`.**
- `test_refuses_missing_georeferencing`: add cases `tiepoint_1_value` (the
  scalar path), `tiepoint_3_values`, `tiepoint_5_values`,
  `scale_1_value` and `scale_2_values`. Each names the tag, its number
  and the count.
- `test_refuses_multiple_tiepoints`: parametrise over 7, 9 and 12 values.
  Each names `33922`, `ModelTiepointTag` and the count. 12 is today's
  catalogue case.
- `test_refuses_model_transformation`: parametrise over 16 values (today's
  case), 3 and 12. The 3 and 12 cases also leave out the tie point and
  scale, so the test shows that 34264 is refused before either is looked at.

**Refusal 8, `TiffFrame`.**
- `test_refuses_ambiguous_pages`: add a case `unparsed_frame`. A two-page
  micro-TIFF whose page 1 is **full-resolution** (subfiletype 0), opened with
  `monkeypatch.setattr(tifffile, "TiffFile", functools.partial(<orig>, _useframes=True))`
  so that page 1 comes back as a `TiffFrame`. It is refused, naming `254`,
  `NewSubfileType`, page `1` and the count `2`. It uses a private
  tifffile argument. That is accepted in a test, and the test should say
  why: no public route yields a frame without an LSM, NDPI or ScanImage
  file. Measured on the green code: this file is **accepted**. The second
  full-resolution image is dropped without a word, which is exactly what
  refusal 8 is for. A reduced-resolution frame is refused too under the
  ruling. That is the accepted cost, and the test does not need to pin it.

**Refusal 11, codec message.**
- `test_refuses_missing_codec`: add a case `unknown_compression`,
  Compression = 60000. It names `259`, `Compression` and `60000`, and
  asserts that `(unknown)` does not appear.
- `test_undecodable_scheme_with_extra_present_gives_no_install_advice`
  (new, `needs_codecs`): Compression = 32809 (`THUNDERSCAN`). It names
  `259`, `32809` and `THUNDERSCAN`, and asserts that `pip install` does
  not appear. This test does not run in CI (see above). It runs wherever the
  extra is installed.

**§6, `bool` sentinel.**
- `TestNoData.test_caller_bool_is_refused` (new), parametrised over `True`,
  `False` and `np.True_`: raises `TypeError` naming `nodata` and `bool`.

**§7, the copy.**
- `TestArray.test_tile_does_not_share_the_callers_buffer` (new): build a
  `DemTile` from a writeable C-contiguous float32 array. Assert that
  `np.shares_memory` is false, and that writing to the caller's array
  leaves `tile.array` unchanged.
- `TestArray.test_tile_array_cannot_be_made_writeable` (new): on a decoded
  tile, `tile.array.flags.writeable = True` raises `ValueError`.

**§12 text.**
- `test_kartverket_fixture_is_shifted_half_a_cell`: its docstring only.
  Remove "The only check that ruling 4 is implemented", and point at
  `TestPlacement.test_area_registered_nodes_are_shifted_inward_half_a_cell`.

The decoder-error ruling (§14, choice C) adds the tests listed there under
"Tests, if C".

---

## 13. LOC estimate

**An estimate, not a measurement** (principle B4). Non-comment production
lines, per `CLAUDE.md` §2.

| file | estimate |
|---|---|
| `src_python/tin_engine/io/models.py` | 90 |
| `src_python/tin_engine/io/geotiff.py` | 250 |
| **total** | **340** |

*Amended.* The earlier figure was 310. The amendments add about 30 lines:
the predictor probe (+3), the 2048 path in CRS resolution (+5), NoData's
file-dtype check, the finiteness rule and the provenance table (+10), the tie
point's `(I, J)` offset and finiteness (+4), the table of tag names for
tagless refusals (+3), and `DemTile`'s shape validator under B1 (+5). Choice A2
in §14, not taken, would have added about 10 lines to
`src_python/tin_engine/__init__.py`.

*Amended (round 2).* About +3, to **343**: the `is_geographic` branch for a
2048 code (+2), and `allow_inf_nan=False` on `nodata` (+1, and it replaces no
line). Treating a caller NaN like a tag NaN adds nothing if both go through
one normalising helper, which is how `@developer` should write it. The table
above is left at its round-1 figures.

*Amended (round 3).* About +15, to roughly **358**. This is still an
estimate, and it is not the reconciliation against the measured count, which
is left to the main session. The parts:
- refusal 13a: +4;
- the tag checks moved before `geotiff_metadata`, the scalar-to-sequence
  step and the scale-count branch: +4;
- the `TiffFrame` branch of refusal 8: +1;
- the codec message (the `find_spec` import and branch; the public `in` in
  place of `_value2member_map_` costs nothing): +3;
- the `bool` guard on `nodata`: +2;
- the copy in `DemTile`: +1.

The decoder-error ruling (§14, choice C) adds about 8 more, so about 366.

`geotiff.py` is mostly the refusal list: eighteen named refusals at roughly
three lines each, plus tag extraction, the CRS resolution, the NoData path and
the promotion table. Comfortably inside the ceiling, which is deliberate —
increment 12 carries the C++ concept change and needs the room.

Tests are excluded from the ceiling and will be larger than the production
code; one micro-TIFF per refusal is the point of the design.

---

## 14. Ruled by the user

Two amendments set two sound principles against each other, so they went to
the user rather than being settled quietly. **The user chose A1 and B1 on
2026-09-24.** The options are kept below so the reasons stay on record.

**Choice A: the "no compiled extension in the process" claim (§2, problem 2).**

- **A1, chosen: make the claim about the source.** `geotiff.py` and
  `models.py` never import `_core`, and `test_module_never_imports_core` checks
  this on the AST. No code changes, and nothing outside the increment is
  touched. `viz/` and `features.py` already make exactly this claim
  (`project_structure.md`, "never imports _core"), so A1 matches the rest of
  the tree. The cost: nobody can run `io/` without first building the
  extension. CI always builds it, so today that costs nothing.
- **A2: make the process-level claim true.** Change `tin_engine/__init__.py`
  to re-export lazily through a module `__getattr__` (PEP 562), with a
  `TYPE_CHECKING` import so that `mypy --strict` still sees `Point2` and the
  others. Then a subprocess test can assert that
  `'tin_engine._core' not in sys.modules` after `import tin_engine.io.geotiff`.
  The gain: isolation that is enforced, not a convention, and it covers
  `viz/` too. The cost: about 10 lines in a file this increment does not own,
  a public-import behaviour change that `tests/python/test_core.py` relies on
  (it should keep passing, but it has to be re-run), and a new failure mode.
  An `ImportError` for a missing extension would move from `import tin_engine`
  to first attribute access.

**Choice B: where `rows` and `cols` live (§7, problem 6).**

- **B1, chosen:** keep them on `RasterMeta` as `StrictInt`, with
  `DemTile` refusing any disagreement with `array.shape`. `RasterMeta` stays a
  complete grid description, which is what the mosaic increment's header-only
  footprint walk needs.
- **B2:** remove them from `RasterMeta`, and make them properties on `DemTile`
  that return `array.shape`. One source of truth, no validator. The mosaic
  increment then re-adds them or builds a second metadata type.

The other five problems were ruled in place by `@architect`. Problem 3 amends
ruling 6, and says so in §1. The user may still want to look at that
amendment: it reads a GeoKey that ruling 6 used to exclude, though only to name
the correct refusal.

### Choice C: errors raised from inside tifffile or a codec

*Added in round 3.* This widens `@reviewer`'s open question 2. **The user
chose C on 2026-09-24.** The options are kept below so the reasons stay on record.

**What escapes.** Every refusal in §5 raises `GeoTiffError`. A file that
tifffile or a codec cannot parse raises whatever that library raises.
Measured through `decode_dem` on the green code, tifffile 2026.9.20, with
and without imagecodecs 2026.8.16:

| input | without imagecodecs | with imagecodecs | a `ValueError`? |
|---|---|---|---|
| empty, or not a TIFF | `tifffile.TiffFileError` | same | yes |
| cut off inside the IFDs | `tifffile.TiffFileError` | same | yes |
| cut off inside the 8-byte header | `struct.error` | same | **no** |
| cut off inside a strip | `ValueError` ("failed to read") | same | yes |
| corrupt Deflate payload | `zlib.error` | `imagecodecs.DeflateError` (a `RuntimeError`) | **no** |
| corrupt LZW payload | refusal 11 (no codec) | `imagecodecs.ImcdError` (a `RuntimeError`) | **no** |
| 33922 or 34264 of the wrong length | bare `ValueError` | same | yes; ruled by round 3 (§5) |

The last row is ruled by round 3 (§5, after refusal 5). It is listed
because it is the same kind of escape, and the same kind of fix.

A corrupt 34735 does not escape. tifffile logs it and returns fewer keys,
and a refusal then fires (§5).

**Why it matters.** "The reader's job is refusal" (ruling 3) is a promise
about what a caller must catch. Today a CLI that catches `GeoTiffError`, or
even `ValueError`, still crashes on a corrupt Deflate strip. And the type it
would have to catch *changes when `imagecodecs` is installed*.

**Options.**

- **A. Wrap everything at the boundary.** One `try` around the whole body
  of `decode_dem`; any `Exception` that is not already `GeoTiffError`
  becomes `GeoTiffError("could not decode: <type>: <message>")`, chained
  `from` the original. About 4 lines. The caller catches one type. The
  cost: it also wraps the reader's **own** bugs. The measured
  `TypeError` for a 1-value 33550, before round 3, would have looked like
  a file refusal, and nobody would have filed it as a bug. It breaks §5's
  "every refusal names the tag", because no tag is known here.
- **B. Wrap a named list** (`TiffFileError`, `struct.error`,
  `zlib.error`, the imagecodecs error classes). About 6 lines. It never
  wraps a reader bug. The cost: the list depends on which codec backend is
  installed, as the Deflate row shows. imagecodecs has one error class per
  codec, and they are `RuntimeError`s. The list is exactly the kind of
  enumeration that §8 refused to keep for codec names, because it rots
  without anyone seeing.
- **C. Wrap by call site, not by type** *(chosen)*. Wrap only the
  three places that hand control to tifffile: `TiffFile(source)` together
  with the page walk, `geotiff_metadata`, and `asarray`. Catch
  `Exception` there, re-raise `MemoryError` untouched, and raise
  `GeoTiffError` chained `from` the original. The message names the stage
  ("TIFF structure", "GeoKey directory", "pixel data"), the original type
  and the original text. About 8 lines. The caller catches one type, and
  the list of exception types never has to be known. The reader's own code
  runs outside the wrapped regions, so its bugs still surface as bugs.
  The cost: the most lines of the three. The "names the tag" rule gets an
  explicit exception for this one class of failure, which names the stage
  instead. And a tifffile bug would be reported as a bad file.
- **D. Document it.** `decode_dem`'s docstring and §5 say that refusals are
  `GeoTiffError`, and that a file tifffile cannot parse raises what tifffile
  raises. 0 lines. The cost: it gives up on ruling 3 for exactly the files
  most likely to arrive broken (truncated downloads). The caller must catch
  `Exception`, which is option A done at every call site instead of once.

**Why C.** A and C both keep the promise to the caller. C keeps it without
turning the reader's own defects into refusals, and the round-3 `TypeError`
shows that those defects exist. B keeps a list that rots. D breaks the
promise. What C costs is lines, and the budget has room (§13).

Where C and A differ is the principle "a bug should surface as a bug". A
holds the boundary and ignores that principle. C holds both, at the price
of a few lines and a softer "names the tag" rule. It is sent to the user
because D is also a coherent choice: a reader that refuses what it
understands and lets through what it does not.

**Tests (C was chosen).** `test_undecodable_input_is_a_geotiff_error`
(new), parametrised over `empty`, `not_tiff`, `truncated_header`,
`truncated_ifd`, `truncated_strip` and `corrupt_deflate`. Each asserts
`GeoTiffError`, and that `__cause__` is set. Add a `needs_codecs` case,
`corrupt_lzw`. Under C only: `test_reader_bug_is_not_wrapped`, which
monkeypatches `_placement` (or its round-3 successor) to raise `TypeError`
and asserts `TypeError` is what comes out. Under D: no tests, and a
docstring change.
