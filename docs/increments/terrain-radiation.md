# Terrain radiation: a radiation geometry per mesh, and time series on demand

Status: **designed; TR1 built by increment 46's PR 2**, the rest not started (2026-09-28, `@architect`, branch
`design-terrain-radiation`; brought to master by increment 46's PR 2). TR1 has code
and tests in that PR; the rest is design only. The work is split
into seven increments, TR1-TR7 ("Increments and LOC"). ROADMAP places it after
auto-catchment. TR1 (the sun's position) depends on nothing and can move earlier
if Ola wants. Written before `@tester`, per `docs/increments/README.md` step 1.

## Ruled by Ola

- **2026-09-28, the purpose**, verbatim: the sun position should be used "to
  compute shades and reflections in a catchment, important for hydrology and
  especially snow melt, but also evapotranspiration. Later, we will also use
  the DOM10 for this." "Reflections" means radiation reflected off the
  surrounding terrain.
- **2026-09-28, the consumers**, verbatim: "the most imminent would be either
  shyft or some ML model for hydrological forcasting. I would also like to use
  it for solar park analysis, or even crops."
- **2026-09-28, time resolution**, verbatim: "Monthly averages is not really
  useful, as they blur out the max/mins who are important for the energy
  balance calculations. Sun on a cold winter day leads to sublimation."
- **2026-09-28, the shape** (the main session proposed it and Ola said yes):
  - a fixed **radiation geometry** per mesh: slope, aspect, horizon angles per
    direction, sky-view factor, terrain-reflection factor, and later a second
    set from DOM10;
  - a fast **evaluator** in the C++ core that computes per-triangle factors
    for any timestamps on demand, without storing them;
  - exact **lit intervals** per triangle and day, derived from the rest as a
    convenience;
  - thin **consumer adapters**: Shyft, ML tensors, and site polygons for solar
    parks and crops;
  - **forcing** (clear-sky or measured radiation) comes from the consumer.
- **2026-09-28, the hard requirements** Ola raised. Each has a ruling below;
  "Ola's requirements, one by one" maps each one to it.
  1. Azimuth matters at high latitudes. Horizon angles per direction are the
     primary data, the sun's path is followed through time, and lit intervals
     are derived from those two. Also: azimuthal resolution tied to latitude,
     atmospheric refraction, polar night and midnight sun, several lit
     intervals in one day, and diffuse and reflected radiation dominating at
     low sun.
  2. The globe is spherical and our CRSs are Cartesian. This needs a local
     frame per triangle (grid convergence, the local vertical), a horizon
     search corrected for Earth curvature and refraction, a search that
     reaches tens of km beyond the catchment, and handling for UTM scale error
     and for the basin's geographic DEM (15c).
  3. Size. 10 M triangles × hourly × a year is about 350 GB, so evaluation is
     on demand. The design states the evaluator's cost and the geometry's
     size.
- **2026-09-28, on Q1-Q10: "go with c, and yes to the rest."** Q6 was
  re-posed after Ola pointed at Shyft's Penman-Monteith evapotranspiration,
  which takes net radiation from Shyft's own slope-and-aspect radiation routine
  (thought to follow Allen, Trezza and Tasumi 2006; to verify in Shyft's
  source). Option **(c)**: Shyft keeps its radiation routine and
  Penman-Monteith; rasputin supplies per cell a **shading series** (lit, or
  the lit fraction per time step, from the horizons) and the static
  **sky-view** and **terrain-reflection** factors, which Shyft's routine
  applies to its beam, diffuse and reflected terms. Q1-Q5, Q7, Q9, Q10 as
  recommended. Q8 as recommended: a pyranometer station in a steep Norwegian
  valley first (still to be found), Rjukan second.
- **Ola, same day, on extending Shyft:** "Shyft can be extended, if we can
  come up with a good argument this is important to get the hydrological
  balance improved." So the design owes **evidence**, not just a mechanism:
  an acceptance probe on a real Norwegian catchment that compares Shyft-style
  slope-and-aspect radiation with and without rasputin's horizon shading,
  sky view and terrain reflection, and reports the difference in the
  radiation and in the snowmelt-season energy balance per cell (the winter
  low-sun months above all, where the design expects the largest effect). The
  size of that difference is the argument.

## Scope

**In.**
- The sun's position from a UTC timestamp. This is the legacy
  `solar_position.h`, ported and verified.
- A **receiver set**: points with a surface normal, a geodetic position and a
  local frame. The triangles of a mesh are one receiver set. Solar-panel sites
  and crop fields are others.
- Per receiver: the horizon in N directions, found by searching the DEM,
  including tiles outside the domain, with Earth curvature and refraction
  corrections. Also the sky-view factor and the terrain configuration factor.
- The evaluator: per receiver and timestamp, whether the receiver is lit,
  cos(incidence) and the direct-beam factor, either at an instant or as a mean
  over an interval.
- Lit intervals per receiver over a UTC window.
- Adapters for Python consumers: NumPy tensors for ML, aggregation over site
  polygons, plane-of-array receivers for solar panels, and a Shyft adapter
  once Shyft's input API has been read (Q6).

**Not in.** Computing clear-sky irradiance, atmospheric transmittance, or
splitting global radiation into beam and diffuse. That is forcing, and forcing
is the consumer's (Ola). Also out: multiple reflections between terrain
patches (radiosity), canopy transmissivity, longwave radiation, and DOM10 (a
hook only; see R12). Each is under "Not in scope" with the reason it waits.

## Prior art: legacy and literature

Written before the design. **No web search tool was available to this round.**
Every citation below is from memory unless it says otherwise. Each one is
listed under "Facts to verify" with what needs checking. The numbers this
design depends on were computed here instead of taken from the literature,
and "What was measured" says how.

### Literature

**The sun's position.**
- **Reda, I., Andreas, A. (2004).** *Solar position algorithm for solar
  radiation applications.* Solar Energy 76(5):577-589; also NREL/TP-560-34302
  (revised January 2008). This is what we use, unchanged. Its stated
  uncertainty is ±0.0003° for years −2000 to 6000, given UT1 and ΔT. The
  legacy header follows its equation numbering (4, 24-47). Its test case
  (2003-10-17 12:30:30 local time, UTC−7, Golden CO) is the TR1 oracle.
  **What differs:** we evaluate the geocentric part once per timestep and the
  topocentric part per receiver (R5). This is a regrouping of the same
  equations, not a new method. The 8.8″ parallax term is kept only in the
  point-by-point reference path (R2).
- **Meeus, J. (1998).** *Astronomical Algorithms*, 2nd ed., Willmann-Bell.
  SPA's source, for the equations of time, sidereal time and nutation.
- **Espenak, F., Meeus, J. (2006).** *Five Millennium Canon of Solar
  Eclipses*, NASA/TP-2006-214141. The ΔT polynomials. For 2005-2050:
  ΔT = 62.92 + 0.32217·t + 0.005589·t² s, with t = year − 2000. **What
  differs:** nothing, except that the consumer may pass ΔT. The polynomial
  gives 75.1 s for 2026. IERS measured about 69 s (to verify). The 6 s
  difference moves the sun by 7·10⁻⁵°, which is negligible (R2).
- **Refraction.** SPA's eq. 42 is Bennett's formula (Bennett, G.G. 1982, *The
  calculation of astronomical refraction in marine navigation*, J. Navigation
  35:255-259), scaled by pressure and temperature. SPA applies it only when
  the true elevation is at least −0.8333° (the solar radius plus 0.5667°).
  **What differs:** we need it slightly lower, for receivers on summits whose
  horizon dips below 0°. R2 extends it to −1.8° and holds it constant below.
  The reason is measured under "What was measured".

**Horizons and shading from a DEM.**
- **Dozier, J., Bruno, J., Downey, P. (1981).** *A faster solution to the
  horizon problem.* Computers & Geosciences 7(2):145-151. Horizons for every
  node of a grid profile in O(n) per profile, by walking the profile's upper
  convex hull.
- **Dozier, J., Frew, J. (1990).** *Rapid calculation of terrain parameters
  for radiation modeling from digital elevation data.* IEEE Trans. Geoscience
  and Remote Sensing 28(5):963-969. The 1981 sweep along grid lines in N
  directions, and **the definitions of the sky-view factor V_d and the terrain
  configuration factor C_t used here** (R7).
- **Stewart, A.J. (1998).** *Fast horizon computation at all points of a
  terrain with visibility and shading applications.* IEEE TVCG 4(1):82-93.
  Horizons on a grid or a TIN, all points at once, in O(n log² n) per set of
  directions. This is the reference for "horizon on the TIN".
- **Corripio, J.G. (2003).** *Vectorial algebra algorithms for calculating
  terrain parameters from DEMs and solar radiation modelling in mountainous
  terrain.* Int. J. GIS 17(1):1-23. The sun as a unit vector and shading by
  scanning along the sun's direction. We use the vector form for incidence
  (R5).
- **Hofierka, J., Šúri, M. (2002).** *The solar radiation model for Open
  Source GIS: implementation and applications* (GRASS `r.sun`, `r.horizon`).
  Horizons per direction, stored as rasters. Also **PVGIS**, which serves a
  DEM horizon per point for solar-park work (Šúri, Huld, Dunlop and co-authors,
  from 2005). PVGIS is an independent comparator for TR3 (Q7).
- **Marsh, C.B., Pomeroy, J.W., Spiteri, R.J. (2012).** *Implications of
  mountain shading on calculating energy for snowmelt using unstructured
  triangular meshes.* Hydrological Processes 26(12):1767-1778. Also **Marsh et
  al. (2020)**, *The Canadian Hydrological Model (CHM) v1.0*, Geosci. Model
  Dev. 13:225-247, and its mesher. **This is the closest prior art to this
  whole design:** a variable-resolution TIN for snow hydrology, with terrain
  shading per triangle. What method they use for shading (a horizon table or a
  per-timestep projection) was not checked (to verify). **What differs, as far
  as is known:** they do not store a horizon table with a reusable evaluator
  or compute lit intervals, and do not search beyond the mesh (to verify).
- **Montero, G. et al. (2009).** *Solar radiation and shadow modelling with
  adaptive triangular meshes.* Solar Energy 83(7):998-1012 (authors and pages
  to verify). Shadows on an adaptive TIN by projecting triangles along the sun
  direction, per timestep. This is the legacy's approach and cost model (per
  timestep, over the whole mesh), which R4 replaces.
- **Line-of-sight curvature and refraction.** The drop of a sight line over
  distance d is d²(1−k)/(2R), with a refraction coefficient k ≈ 0.13 (Gauss's
  value, standard in geodesy and in GIS viewshed tools). k varies near the
  ground, and strongly under inversions (for example Hirt, C. et al. 2010,
  *Monitoring of the refraction coefficient in the lower atmosphere using a
  controlled setup of simultaneous reciprocal vertical angle measurements*,
  J. Geophys. Res. 115:D21102). We take k as a parameter (R4, Q3).

**Sky view, diffuse and reflected radiation.**
- **Dozier and Frew (1990)**, above: V_d, the cosine-weighted fraction of the
  isotropic sky a tilted surface sees, and C_t ≈ (1 + cos S)/2 − V_d.
- **Zakšek, K., Oštir, K., Kokalj, Ž. (2011).** *Sky-view factor as a relief
  visualization technique.* Remote Sensing 3(2):398-415. **Their SVF is a
  different quantity**: the unweighted solid-angle fraction, 1 − mean(sin h).
  We use V_d, because V_d is what multiplies diffuse irradiance, and we name
  it so that the two are not confused (R7).
- **Hay, J.E., Davies, J.A. (1980)**, the circumsolar/isotropic split, and
  **Perez, R. et al. (1990)**, *Modeling daylight availability and irradiance
  components from direct and global irradiance*, Solar Energy 44(5):271-289.
  These are diffuse models the combiner may offer (R9). The anisotropy inputs
  are forcing, so they come from the consumer.
- **Helbig, N., Löwe, H., Lehning, M. (2009).** *Radiosity approach for the
  shortwave surface radiation balance in complex terrain.* J. Atmos. Sci.
  66(9):2900-2912. Multiple reflections between patches, a matrix method. It
  is not in scope, but it is the reference for what C_t approximates.
- **Allen, R.G., Trezza, R., Tasumi, M. (2006).** *Analytical integrated
  functions for daily solar radiation on slopes.* Agric. For. Meteorol.
  139:55-73. Integrates cos(incidence) over hour angle in closed form. We
  believe Shyft's radiation routine is built on it (to verify). The interval
  mean (R5) uses the same closed-form integral between the lit endpoints.
- **Iqbal, M. (1983).** *An Introduction to Solar Radiation*, Academic Press;
  and **Duffie, J.A., Beckman, W.A.**, *Solar Engineering of Thermal
  Processes*: textbook background for incidence on tilted planes, for the
  solar-park adapter.

**Novelty: none claimed.** Two things could look new:
- **Earth curvature in the convex-hull horizon sweep.** Along a profile, the
  curvature-corrected tangent to point p from observer o is
  tan h = (w_p − w_o)/d + 2c·s_o, where w = z − c·s², c = (1−k)/(2R) and s is
  arc length along the profile. The derivation is two lines: substitute
  z = w + c·s² into (z_p − z_o − c·d²)/d with d = s_p − s_o. The added term
  does not depend on p, so the maximising point is the upper convex hull in
  (s, w), and Dozier's 1981 sweep stays exact with curvature. It is recorded
  under R4's rejected alternative. It is not claimed. Search before anyone
  claims it: "horizon convex hull earth curvature DEM", "viewshed curvature
  correction line sweep".
- **Certified lit-interval search** (R6): a sampling step bounded by the
  sun's and the horizon's rates of change, so that no lit interval is missed.
  This is standard interval root isolation, and it is not claimed. Search:
  "sunrise sunset horizon profile root finding guaranteed".

Before any publication, search Google Scholar and Web of Science for
"horizon angle TIN snowmelt", "terrain shading unstructured mesh hydrology",
"sky view factor triangulated irregular network", and CHM's shading papers.

### Legacy

```sh
$ grep -rliE 'solar|shade|shadow|sun|horizon|radiation|irradiance|sky.?view|albedo|azimuth|zenith' legacy/ | sort
legacy/bindings.cpp
legacy/rasputin/application.py
legacy/rasputin/calculate.py
legacy/rasputin/geo_tiff_reader.py
legacy/rasputin/mesh.py
legacy/rasputin/reader.py
legacy/rasputin/solar_position.h
legacy/rasputin/tin_repository.py
legacy/rasputin/triangulate_dem.h
legacy/rasputin/writer.py
legacy/tests/test_mesh.py
legacy/tests/test_tin_repository.py
```

`reader.py` matches only on GeoTIFF key names (`GeogAzimuthUnitsGeoKey`,
`ProjAzimuthAngleGeoKey`) and an empty `read_sun_posisions`. The rest is
two pieces of logic and the plumbing around them.

**1. `legacy/rasputin/solar_position.h` (Ola, 2019-08-26, 752 lines): NREL SPA,
in SPA's own notation.** Carried, after a fix to ΔT and a guard on refraction.

- The algorithm is SPA, equation by equation (eqs. 4, 24-47 are labelled). The
  periodic-term tables are inlined as named constants (`L00` ... `L063`,
  `B00`, `R00`, and 63 nutation terms). Time comes in as a
  `std::chrono::system_clock::time_point` and becomes a Julian day through
  `sys_days{January/1/1970}`, which is C++20 `<chrono>` with no external date
  library.
- **Checked against NREL's published test values** (the probe under "What was
  measured"): JD 2452930.312847, L 24.0182616917°, B −0.0001011219°,
  R 0.9965422974 AU, topocentric zenith 50.11162°, azimuth 194.34024°,
  incidence 25.18700° (slope 30°, azimuth rotation −10°). All agree to every
  printed digit. One date does not check every periodic term, so TR1 checks
  the tables term by term as well (R2).
- **Defect: ΔT is wrong by about 56 minutes after 2009.** From 1955 to 2009,
  `coarse_*_calc` look ΔT up in a 5-year table. Outside those years they fall
  back to `parabolic_dt_calc`, which computes `-20 + 32*pow((year - 1820)/20.0, 2)`.
  The Morrison-Stephenson parabola divides by 100, not 20. The code gives
  3375 s for 2026, where the parabola gives 116 s and the true value is about
  69 s. Measured effect at Oslo, 2026-01-15 12:00 UTC: elevation off by
  0.010° and azimuth by 0.038°. That is small, but it is a defect, and it
  grows quadratically with the year.
- **Defect: refraction applied without SPA's guard.**
  `topocentric_zenith_angle` applies eq. 42 at every elevation. Below about
  −1.9° the formula turns non-physical: 0.74° at −2°, 0.18° at −4°, and a
  singularity near −5.1°. SPA itself applies it only above −0.8333°.
- **Refraction was never used for shading.** The collector
  `azimuth_and_elevation()` returns the *uncorrected* elevation e0, and
  `shade` uses that value. So the legacy shading has no refraction: 0.48° at
  the horizon, which at 65°N in December is about 14 minutes at each end of a
  2.9-hour day.
- `jd_from_clock` truncates to whole seconds, which is harmless. Its
  `system_clock` is Unix time, UTC without leap seconds, which is right for
  SPA's UT to within |UT1 − UTC| < 0.9 s (R2).

**2. Per-face shading: `legacy/rasputin/triangulate_dem.h` lines 599-727
(`compute_shadow`, `is_shaded`, `shade`, `compute_shadows`), via CGAL's AABB
tree.** The idea is carried. The mechanism is replaced.

- What it does. For each face, per timestamp: the face is shaded if its
  normal faces away from the sun (self-shadow). Otherwise it is shaded if a
  ray from the face's centroid toward the sun hits any other face of the
  mesh (`CGAL::AABB_tree::first_intersection`, which skips the face itself).
  `shade` gets each face's lat/lon from `boost::geometry::srs` using the
  mesh's proj4 string, then runs SPA per face and per timestamp. Output: one
  bool per face per timestamp. `legacy/rasputin/application.py:159-198` loops over
  timestamps and writes every array to HDF5/XDMF (`ShadeRepository`).
- **Carried:** the two tests and their order, self-shadow (the normal against
  the sun) and then occlusion; a boolean lit status as one output; SPA at each
  face's own lat/lon; UTC `time_point`s.
- **Replaced, and why:**
  - *CGAL AABB ray casting per timestep* is replaced by **a horizon table per
    receiver, built once** (R3, R4). This is required by the No-CGAL mandate
    anyway. The cost model is the bigger reason: ray casting costs
    O(faces × timesteps × log faces), and storing the result is the 350 GB
    problem.
  - *Occluders limited to the mesh.* A ray leaving the catchment hits nothing,
    so every mountain outside the domain was ignored. R4 searches the DEM tens
    of km past the domain.
  - *No grid convergence.* `is_shaded` builds the sun vector as
    `(sin az, cos az, tan el)` directly in the projected frame, treating grid
    north as true north. In UTM33 that is an azimuth error of −8.9° on the
    Vestlandet coast and +15.2° at Vardø (measured below). R3 corrects it.
  - *No Earth curvature, no refraction* (above).
  - *`compute_shadow(mesh, direction)` rebuilds the AABB tree inside its
    per-face loop*, O(n² log n) per call. `compute_shadows` returns an empty
    vector (`// TODO: Rewrite totally!`).
  - *Per-face proj4 transform* through Boost.Geometry. That is CRS work in
    C++, which the I/O boundary now forbids. R10 moves it to Python, once per
    receiver.
- `compute_slope` and `compute_aspect` (lines 757-780) take slope as
  atan2(|n_xy|, n_z) and aspect as atan2(n_x, n_y), in grid axes. Carried as
  definitions, now in true ENU axes (R3).

**`@migration-expert`:** no domain constant is re-derived here. The SPA
coefficients are copied from the published tables, and TR1's oracle is NREL's
published values (plus an independent implementation, Q2), not the legacy.
Recommended: one short `@migration-expert` pass before TR1's red step. It
diffs the legacy's inlined terms against the paper's tables A4.2 and A4.3
term by term, because a transcription error in a small term can survive a
single-date check.

## What was measured

Run on 2026-09-28 on Ola's Mac (Apple M1 Max, on battery), with Apple clang
21.0.0 and the project venv (pyproj 3.8.0, PROJ 9.8.1). The probes were
throwaway scripts in a scratch directory. They are not committed, since this
round is design only. Each recipe is below so that a reader can re-run it.

**P1. Legacy SPA against NREL's test case.** Compile a file that includes
`legacy/rasputin/solar_position.h` with `-std=c++20` and call
`calendar_solar_position(2003, 10, 17 + 19.508333/24, 39.742476, -105.1786,
1830.14, …)` with ΔT = 67 s, then `corrected_solar_elevation(e0, 820, 11)`
and `indicence_angle(z, 30, -10, Γ)`. Output: `JD 2452930.312847`,
`L 24.0182616917 B -0.0001011219 R 0.9965422974`,
`zenith 50.11162 azimuth 194.34024 incidence 25.18700`. These equal the
published values quoted in the Legacy section. Those were recalled from the
paper and not reread (to verify). Two independent sources agree on all ten
significant digits, so the recollection is very unlikely to be wrong.
It compiled unchanged under C++20. `<chrono>`'s calendar types
(`sys_days`, `January/1/1970`) are in Apple's libc++.

**P2. Legacy ΔT.** `parabolic_dt_calc()(2026)` = 3374.9 s; with /100 it
would be 115.8 s. SPA at Oslo (59.91°N 10.75°E), 2026-01-15 12:00 UTC:
ΔT = 69 s gives e0 = 8.72791°, az 187.92892°; ΔT = 3374.9 s gives
e0 = 8.73784°, az 187.89079°.

**P3. SPA's refraction formula, eq. 42 at P = 1010 hPa, T = 10 °C**, as a
function of true elevation e0: −4° → 0.184°, −3° → 0.518°, −2° → 0.742°,
−1° → 0.647°, −0.5667° → 0.573°, 0° → 0.483°, 1° → 0.362°, 3° → 0.228°,
10° → 0.090°. Refraction grows as e0 falls until about −1.9°, then falls,
which is not physical. **It is monotone on [−1.8°, 90°].** The dip of the
horizon seen from Galdhøpiggen (2469 m) is 1.49° with k = 0.13 (1.60°
geometric), so a Norwegian receiver never needs the formula much below
about −2.2°. R2 holds the refraction constant below −1.8°. The error this
costs is under 0.1°, less than what the atmosphere's day-to-day variation
does at the horizon anyway.

**P4. Grid convergence in EPSG:25833**, from
`pyproj.Proj(CRS.from_epsg(25833)).get_factors(lon, lat)`, cross-checked by
transforming a 100 m step along the meridian:

| site | lon, lat | `meridian_convergence` | grid azimuth of true north | point scale m |
|---|---|---:|---:|---:|
| Vestlandet coast (Solund) | 4.85, 61.08 | −8.906° | +8.906° | 1.00326 |
| Rjukan | 8.59, 59.88 | −5.550° | +5.550° | 1.00118 |
| Oslo | 10.75, 59.91 | −3.679° | +3.679° | 1.00029 |
| central meridian | 15.00, 65.00 | 0.000° | 0.000° | 0.99960 |
| Tromsø | 18.96, 69.65 | +3.714° | −3.714° | 0.99989 |
| Vardø | 31.11, 70.37 | +15.218° | −15.219° | 1.00397 |

This agrees with Ola's measurements (−8.9°, −3.7°, 0°) and with
(λ − 15°)·sin φ to within 0.05°. Vardø comes out at +15.2° here against
Ola's +15.1°, most likely from a slightly different point. **The sign
convention, pinned:** pyproj's `meridian_convergence` γ is the true azimuth
of grid north, so **A_grid = A_true − γ**. A test pins this at Solund (R3).
The scale error is 0.3-0.4 % at the zone's edges, as Ola said.

**P5. Sight-line drop, d²(1−k)/(2R)** with R = 6371 km and k = 0.13:

| d | drop | angle | with k = 0 | with k = 0.25 |
|---:|---:|---:|---:|---:|
| 1 km | 0.1 m | 0.004° | 0.004° | 0.003° |
| 5 km | 1.7 m | 0.020° | 0.022° | 0.017° |
| 20 km | 27.3 m | 0.078° | 0.090° | 0.067° |
| 50 km | 170.7 m | 0.196° | 0.225° | 0.169° |
| 100 km | 682.8 m | 0.391° | 0.450° | 0.337° |

These equal Ola's figures (27 m / 0.08°, 171 m / 0.20°, 683 m / 0.39°). The
angle is the drop divided by d, which is d(1−k)/(2R).

**P6. How fast the sun moves at the horizon** (spherical formulas,
finite-differenced and compared with the closed forms
de/dt = 15°/h · cos φ · sin A and dA/dt = 15°/h · sin φ, which hold at
e = 0; they agree to 0.02°/h):

| lat | declination | sunrise azimuth | de/dt | dA/dt | day length | noon elevation |
|---:|---:|---:|---:|---:|---:|---:|
| 60° | −23.44° | 142.7° | 4.5°/h | 13.0°/h | 5.5 h | +6.6° |
| 65° | −23.44° | 160.3° | 2.1°/h | 13.6°/h | 2.9 h | +1.6° |
| 70° | −23.44° | polar night | | | | −3.4° |
| 70° | −10° | 120.5° | 4.4°/h | 14.1°/h | 8.1 h | +10.0° |
| 78° | −10° | 146.6° | 1.7°/h | 14.7°/h | 4.5 h | +2.0° |

(geometric, no refraction or solar disc.) At high latitude, low sun moves
three to eight times faster in azimuth than in elevation. **An error δh in a
horizon angle becomes an error of about δh / (de/dt) in when the sun
appears.** At 65°N in December: 0.1° → 2.8 min; refraction ignored (0.48°) →
14 min per end; curvature ignored at 20 km (0.078°) → 2.2 min; the solar
semi-diameter (0.267°) → 7.5 min. Near local noon on a day when the sun barely
clears the horizon, de/dt → 0 and the sensitivity has no bound. That is
exactly the "sun on a cold winter day" case, and why R6 certifies its root
search instead of sampling.

**P7. The evaluator's inner loop, prototyped.** A 60-line C++ loop over
1 M receivers × 24 timesteps (`-O2`, one thread, float32, N = 64 horizon
directions stored as int16 in units of 0.01°). Per receiver-timestep it does
three dot products against a precomputed east/north/up basis, `asin` for
the elevation, the refraction formula, `atan2` for the azimuth, a linear
interpolation in the horizon table, the lit test, and cos(incidence) as one
more dot product. **28.7 ns and 28.9 ns per receiver-timestep** on two runs,
on battery. This is a prototype of R5's kernel, not the kernel. `@perf`
measures the real one in TR5.

**P8. Dozier and Frew's V_d against brute force.** The formula in R7,
integrated numerically over 10⁵ azimuths, against a direct cosine-weighted
integration over the sky hemisphere (1500 × 3000 cells):

| case | V_d (formula) | brute force | (1 + cos S)/2 |
|---|---:|---:|---:|
| flat, horizon 0° | 1.00000 | 1.00000 | 1.00000 |
| flat, horizon 10° all round | 0.96985 | 0.96973 | 1.00000 |
| 30° slope, no terrain | 0.93301 | 0.93301 | 0.93301 |
| 35° slope facing S, 15° ridge to the N | 0.90736 | 0.90736 | 0.90958 |
| 35° slope, 20° horizon all round | 0.73403 | 0.73425 | 0.90958 |

**The formula is right only if each direction's horizon is first raised to
the slope's own plane and clamped at ≥ 0°.** Without the clamp, the 30° slope
with no terrain gives 1.0 instead of 0.933. R7 states the clamp. The
remaining differences, up to 2·10⁻⁴, come from the brute force's grid.

## Rulings

### R1. The blueprint: components, data flow, and the boundary

Four stages. Each is a pure function of its inputs, and each can be tested
without the stages after it.

```
 PYTHON (CRS, files, time zones, forcing)            C++ CORE (numbers only; no path, no CRS)
 ─────────────────────────────────────────           ─────────────────────────────────────────
 RadiationRequest (frozen Pydantic)
   mesh, DEM source, n_directions, k, ...
        │
 [A] receivers_from_mesh(mesh, frame) ──────────▶ ReceiverSet (SoA spans)            [TR2]
     centroid → lat, lon      (crs.reprojector)     x, y, z_obs in the frame; lat, lon;
     γ, m per receiver        (crs.frame_factors)   J (frame → local east/north metres);
     J = (1/m)·R(γ); normal in ENU                  normal (E, N, U); exclusion radius
        │
 [B] plan_search(receivers, D_near, D_far) ──────▶ horizon_search(receivers,           [TR3]
     near mosaic: 15a repository, full res          near RasterView, far RasterView,
     far grid: block means of the same tiles        c, N, ...) → int16[n × N], reach
        │                                          sky_factors(normal, horizon)        [TR4]
        │                                            → V_d, C_t
        ▼
 RadiationGeometry  (persisted, .npz; R11)
        │
 [C] evaluate(geometry, times, mode) ───────────▶ SunEphemeris(times, ΔT)              [TR5]
     times: tz-aware → UTC int64 ms                 evaluate_instant / evaluate_mean
     asyncio.to_thread, GIL released                → lit, cos_i, sin_e per receiver·time
     lit_intervals(geometry, window) ───────────▶ lit_intervals → CSR (offsets, starts, [TR6]
                                                    ends)
        │
 [D] adapters: combine with forcing, tensors, ◀── (factors are per unit inclined area)  [TR7]
     site polygons, panels, Shyft
```

- **The solar position** (`radiation/solar_position.hpp`, TR1) is used by
  [C] and by the reference path in tests. It does not depend on any mesh.
- **The receiver set is the abstraction, not the mesh.** Triangles, panel
  sites and crop points are all receivers. Stages [B] and [C] never see a
  triangle index. The evaluator can be unit-tested on three hand-made
  receivers, and the horizon search on a synthetic raster with one receiver.
- **Config and execution are separated.** `RadiationRequest` and
  `EvaluationRequest` are data (Pydantic, frozen). The functions that execute
  them live in their own modules. The CLI (`rasputin radiation ...`, TR3b and
  TR5) only parses flags into requests.
- **Homes.** C++: `include/terrain/radiation/` (`time.hpp`,
  `solar_position.hpp`, `receivers.hpp`, `horizon.hpp`, `sky.hpp`,
  `evaluate.hpp`, `lit_intervals.hpp`), header-only like the rest of the
  core. Bindings: a new `bindings/radiation.cpp`, so that `bindings/core.cpp`
  does not grow (`tools/check_detria_boundary.py` is unaffected; it only
  confines `detria.hpp`). Python: `src_python/tin_engine/radiation/`
  (`request.py`, `frame.py`, `search.py`, `geometry.py`, `evaluate.py`,
  `adapters/`). Also `crs.py` gets `frame_factors`.

### R2. Time and the sun's position [TR1]

- **UTC throughout, as Unix time.** In C++,
  `std::chrono::sys_time<std::chrono::milliseconds>`. In NumPy,
  `datetime64[ms]`. At the Python boundary a timezone-aware `datetime` is
  converted to UTC. **A naive `datetime` or a string without an offset is
  refused**, because guessing a zone is how a one-hour shift in melt enters.
  Unix time has no leap seconds. SPA wants UT1, and |UT1 − UTC| < 0.9 s
  moves the Earth's rotation angle by at most 0.004°. That is below
  anything else here, so UTC is used as UT1 and this is written down.
- **Only `<chrono>`.** `sys_days`, `year_month_day` and duration arithmetic,
  as the legacy already did (P1 compiled under Apple libc++). No time-zone
  database is needed, because nothing leaves UTC in C++. No external date
  library (CLAUDE.md §2).
- **ΔT** is a parameter with a default. The default is Espenak and Meeus's
  polynomial for the timestamp's year (their 1986-2005 and 2005-2050 pieces,
  and the long-term parabola with **/100** outside them). A consumer who has
  IERS values passes them. The legacy defect (P2) is fixed and pinned by a
  test at 2026.
- **The port.** SPA's periodic terms become `constexpr` tables, one row per
  term (A, B, C), instead of the legacy's named constants. The equations stay
  in SPA's order and notation, with SPA's equation numbers in comments, as
  Ola wrote them. The layout splits in two:
  - `geocentric_sun(t, ΔT) -> GeocentricSun {α, δ, ν, R}` depends only on
    time. It is the expensive part: about 200 trigonometric terms.
  - `topocentric_sun(GeocentricSun, lat, lon, h) -> {e0, A}` is the cheap,
    per-location part, including the parallax correction (eqs. 33-40).
  - `apparent_elevation(e0, P, T)`: eq. 42, applied for e0 ≥ −1.8° and held
    at its −1.8° value below (P3). SPA cuts it off at −0.8333°. We extend it
    because a summit's horizon dips below 0°, down to −1.5° at Galdhøpiggen
    (P3).
  - `solar_position(t, lat, lon, h, ΔT, P, T)`: the full SPA, point by point.
    This is the **reference path** that tests compare everything else against.
- **Pressure and temperature** default to a standard atmosphere at the
  receiver's height (P from the barometric formula, T = 10 °C). The consumer
  may pass P and T per timestep. At the horizon, −20 °C instead of 10 °C adds
  12 % to refraction (0.06°), which at 65°N in December is about 2 minutes.

### R3. The receiver frame: convergence, scale, the local vertical [TR2]

The mesh lives in a Cartesian computation frame: a projected CRS for Norway,
and 15c's lattice frame for a geographic DEM. The sun lives in true
directions on an ellipsoid. **Each receiver gets a local frame, computed once
in Python and passed as numbers:**

- **Position:** the receiver's (x, y) in the computation frame, and its
  geodetic latitude and longitude from `crs.reprojector(frame CRS,
  EPSG:4326)`. ETRS89 against WGS84 is under 1 m, which is 10⁻⁵° and
  negligible for the sun. The core treats lat/lon as numbers on an ellipsoid.
  It never learns which CRS they came from, so the I/O boundary holds
  (CLAUDE.md §2, R10).
- **J, a 2×2 matrix per receiver**, maps a displacement in the frame to
  local true (east, north) metres.
  - Projected CRS: J = (1/m)·Rot(γ), where Rot(γ) turns azimuths clockwise
    by γ, γ is pyproj's `meridian_convergence`, and m is the point scale
    (written m because k is the refraction coefficient in R4). Both come
    from a new `crs.frame_factors(crs, lon, lat)` that wraps
    `Proj(crs).get_factors`, next to the one `Transformer.from_crs` site.
    The sign is pinned by P4: A_grid = A_true − γ, so a grid vector with grid
    azimuth A_g has true azimuth A_g + γ. A test pins the Solund value
    (−8.906°) and the Vardø value (+15.218°) through J, not only through γ.
  - The lattice frame (15c): one column step is du degrees of longitude and
    hx frame units, one row step dv degrees of latitude and hy frame units,
    so J = diag(N(φ)·cos φ·(π/180)·du/hx, M(φ)·(π/180)·dv/hy), which is
    about diag(N cos φ / a, M / a). N and M are the ellipsoid's radii of
    curvature. It is anisotropic, as 15 R8 says. Only the receiver frame
    needs to know that.
  - **One representation for both.** The evaluator never sees γ or m. The
    geometry build uses J for two things: the normal, and the direction of
    each horizon ray in the frame.
- **The surface normal in true ENU.** Transform the triangle's edge vectors'
  horizontal parts by J, keep z, and take the cross product, oriented upward.
  Slope S = atan2(|n_EN|, n_U) and aspect = atan2(n_E, n_N), true azimuth
  (the legacy definitions, now in true axes). m enters slope as a factor of
  at most 1.004 on tan S, which is small, but applying it costs nothing.
- **The local vertical is the ellipsoid normal at (lat, lon).** DEM heights
  are orthometric (NN2000). The deflection of the vertical in Norway is some
  arcseconds, 10⁻³° at most, and is ignored and written down.
- **Earth radius for the curvature term:** Python passes one number, R =
  √(M·N) (the Gaussian mean radius) at the receivers' mean latitude. In the
  direction of the ray the radius varies by ±0.3 % between M and N. The drop
  at 50 km then changes by 0.5 m, which is negligible.
- **Exclusion radius r_i:** for a triangle, the largest ground distance from
  its centroid to a vertex. For a point receiver, a parameter (default half a
  DEM cell). R4 says why.
- **Observer height:** z_obs = the TIN's z at the centroid (the triangle's
  plane) plus a per-receiver offset. The offset is 0 for terrain, the panel
  height for panels, and the canopy height for crops if wanted.

### R4. The horizon: searched on the DEM, per receiver, with curvature and refraction [TR3]

**On the DEM grid, not on the TIN.** There are three reasons:
1. **The occluders are mostly outside the mesh.** The mesh stops at the
   catchment. The horizon of a valley floor is set by ridges that are often
   outside it, and at low sun by mountains tens of km away (P5, P6). No TIN
   exists there, and meshing the surroundings only to shade the inside costs
   more than searching the grid.
2. **The DEM is the reference.** The TIN is the DEM within `--tolerance`.
   Searching the reference means the horizon does not depend on the
   tolerance the mesh was built with.
3. **Regular grids are fast to search and easy to bound.** Row-major memory,
   bilinear sampling (`raster/sample.hpp` exists already), and a max-pyramid
   for pruning. Stewart's TIN method (1998) is cleverer and asymptotically
   better per point set, but it needs the whole terrain as a TIN.

**Per receiver, not per DEM node.** Dozier and Frew (1990) sweep each grid
line and get every node's horizon in O(n) per direction. That is the right
tool for a horizon raster. It is not the right tool here, because:
- **our receivers are not DEM nodes.** They are triangle centroids, panel
  sites and crop points.
- **each receiver has its own near-field exclusion** (below), and a sweep
  has no per-observer start distance.
- interpolating node horizons to centroids brings in the DEM's micro-relief
  inside the triangle's footprint. That is the TIN-versus-DEM mismatch the
  exclusion radius exists to avoid.

*Rejected alternative, recorded:* the Dozier sweep with curvature. It stays
exact under Earth curvature: in heights w = z − c·s² the maximising point is
still the upper convex hull ("Novelty", above). If someone later wants a
horizon raster (every DEM node, for a gridded product), that is the method.

**The reference definition. Tests pin this; the fast path must equal it bit
for bit.** Receiver i, direction k:
- True azimuth A_k = k·360°/N, k = 0..N−1, so directions are **equally
  spaced in true azimuth, not grid azimuth.** The ray's direction in the
  frame is u_ik = J_i⁻¹·(sin A_k, cos A_k), so one unit of the ray parameter
  is one ground metre.
- Samples at ground distance d_j = r_i + j·Δ_near for d_j ≤ D_near, read
  bilinearly from the **near raster** (the DEM at full resolution), with
  Δ_near = half the near cell. From D_near to D_far, Δ_far = half the far
  cell, on the **far raster** (R4, "Search radius").
- tan h_ik = max over j of (z(d_j) − z_obs,i − c·d_j²) / d_j, with
  c = (1 − k)/(2R) (R3, Q3), and h_ik = atan of that. With no valid sample
  at all, h_ik = −90°, and `reach` says so.
- Stored as **int16 in units of 0.01°** (±327°). The quantisation error is
  ≤ 0.005°, which is 8 s at 65°N in December (P6).
- **Refraction, both kinds, consistently.** The terrain ray bends by
  terrestrial refraction k, which is folded into c. The sun's apparent
  elevation gets astronomical refraction (R2). Lit means apparent sun above
  apparent horizon. Using both is right: the two light paths are different
  paths through the same atmosphere.

**The fast path: a max-pyramid, with pruning that must not change the
result.** Each raster gets a pyramid of maxima over 2^L × 2^L node blocks.
The march along a ray skips a stretch [d_lo, d_hi] whose covering pyramid
cell has maximum Z when this upper bound cannot beat the current best:
a/d_lo − c·d_lo if a = Z − z_obs ≥ 0, and a/d_hi − c·d_lo if a < 0.
A bilinear value never exceeds its four nodes, and the bound holds for every
d in the stretch, so **every sample the exhaustive loop would take either is
taken, with the same arithmetic, or provably cannot raise the max.** The
result equals the exhaustive one bit for bit, whatever the thread count.
This is TR3's invariant-critical suite. The same bound at the top level ends
the ray early, once even the raster's global maximum cannot beat the best.

**Near field: the triangle is its own plane.** Samples start at the
exclusion radius r_i, not at 0. Inside its own footprint the triangle *is*
the surface: self-shadowing is its normal against the sun (R5), exactly as in
the legacy. Starting at 0 would make DEM nodes that sit up to `--tolerance`
above the triangle's plane cast horizons onto it: 1 m at 10 m is 5.7°. Beyond
r_i the DEM and the neighbouring triangles agree within the tolerance.

**Search radius and the far field.**
- The furthest anything can matter is where even the highest terrain drops
  below the observer's horizontal: D_max = √((z_max − z_obs)/c). For 2469 m
  seen from sea level that is 190 km. The pyramid's early stop reaches that
  bound on its own, so D_far is a cap, not a guess.
- **Two rasters.** The **near raster** is the DEM at full resolution over the
  domain's box grown by D_near (default 20 km). It is assembled by 15a's
  `plan_mosaic`/`assemble` with that box as `bounds`, reading tiles outside the
  domain. A 30 km catchment at 10 m gives a 70 km box: 49 M nodes, 196 MB
  float32, plus a third of that for the pyramid. The **far raster** covers the
  box grown by D_far (default 200 km) at a coarse spacing (default 10 cells,
  100 m). It is built by **block means**, streaming one tile at a time
  through the same repository, so no whole-area 10 m canvas ever exists
  (for the 30 km catchment, a 430 km box at 100 m: 18.5 M nodes, 74 MB). Block means, not maxima: a mean is
  the surface sampled at 100 m, unbiased. A peak at 50 km is lowered by tens
  of metres at most, which is about 0.03°.
- **Resampling is allowed here, and only here.** 15 R8 forbids resampling
  because it moves the surface that the **mesh's sup-norm tolerance** is
  measured against. The far raster feeds only an angle, and has no tolerance
  to lose. The near raster is not resampled.
- **Coverage.** DTM10 ends at the border and the coast. Beyond it lie
  Swedish, Finnish and Russian mountains and the sea. A sample outside the
  mosaic or on NoData is **no obstacle**, and the receiver records, per
  direction, the distance searched before the first such sample: **`reach`,
  uint8 in km** (255 = the full D_far). A receiver near the Swedish border
  thus says "east horizon searched to 12 km" instead of pretending. What to
  use beyond Kartverket's coverage is Q4.
- **Why not ray-march the whole thing at 10 m?** It is correct, and the
  pyramid makes it affordable, but the far field at 10 m needs a 430 km
  canvas of 1.85 G nodes (7.4 GB). The two-raster split exists to fit memory,
  not to save work.

**Cost, estimated and to be measured by `@perf` in TR3.** Exhaustive, per
ray: 20 km / 5 m + 180 km / 50 m ≈ 7600 samples. At about 3 ns per bilinear
sample that is 23 µs per ray and 1.5 ms per receiver at N = 64. For 10 M
receivers that is 4 h on one core, about 30 min on 8. The pyramid should
remove most of it, since in the usual case the far field is far below the
best near-field tangent, but the factor is unmeasured. It is a one-off per
mesh and DEM. At Shyft's scale (10⁵-10⁶ receivers) it is seconds to minutes
even without pruning.

**Deterministic and parallel.** Receivers are independent: a parallel-for
over receivers (`parallel_util/chunks.hpp`), no shared mutable state, the
rasters read-only (14 R7's concurrency rule). The output does not depend on
the thread count.

### R5. The evaluator: factors for any timestamps, on demand [TR5]

**Per timestep, once for all receivers:** the geocentric sun (R2) as a unit
vector in Earth-fixed coordinates, s_t = (cos δ cos G, −cos δ sin G, sin δ),
with G = ν − α the Greenwich hour angle. The cost is a few µs per timestep and
does not grow with the mesh.

**Per receiver, once at construction:** the east, north and up unit vectors
E_i, N_i, U_i in the same Earth-fixed frame, from (lat, lon), and the normal
n_i expressed in them.

**Per receiver and timestep** (the inner loop prototyped in P7):
- sin e = s_t·U_i, A = atan2(s_t·E_i, s_t·N_i) (true azimuth), and
  e_app = e + refraction(e; P, T).
- h = the horizon at A, linear between the two neighbouring directions.
- **lit** = e_app + limb > h **and** the apparent sun is in front of the
  triangle's plane (n_i·s_app > 0). `limb` is 0 (the sun's centre, default)
  or +0.2667° (the upper limb, for comparing with almanac sunrise times).
- **cos_i** = n_i·s_app, unclipped, for the consumer's own use. The beam
  factor is **b = lit · max(0, cos_i)**, the multiplier on direct normal
  irradiance (DNI).
- **sin_e_app**, the multiplier that turns DNI into beam on the horizontal,
  so that a consumer holding only horizontal beam can form the ratio.
- The parallax term (≤ 8.8″ = 0.0024°) is dropped in this fast path. The
  evaluator is tested against the reference path (R2) within 0.003° in
  elevation and azimuth.

**Two modes.**
- `instant`: the factors at each timestamp. Cost: about 29 ns per
  receiver-timestep on one core (P7). 10 M receivers × 8760 hours is
  8.8·10¹⁰ evaluations: 42 min on one core, about 6 min on 8 if it scales.
  One timestep for 10 M receivers takes 0.3 s on one core.
- `mean` over [t_j, t_j + Δ): the time-mean of b, of lit and of sin_e_app.
  **This is what an hourly energy balance needs.** An instantaneous value at
  the top of the hour misstates every hour that holds a sunrise or a
  shadow edge, and at 65°N in December that is most of the day's hours.
  Method: take the lit sub-intervals from R6 restricted to the step, then
  integrate cos_i over each with 4-point Gauss-Legendre. The pieces are
  smooth, because the lit endpoints are exactly where the kinks are.
  Allen et al.'s closed form (2006) is exact only for fixed declination and
  no refraction, and the quadrature costs about the same. Cost: R6's search
  once per receiver-day (estimated 1.5-6 µs) plus about 4 × 30 ns per
  receiver-step. For 10 M receivers × a year of hours that is roughly 7 h
  on one core, 1-2 h on 8. For Shyft's 10⁵ cells it is about a minute on
  one core.

**Loop order for long series:** chunks of receivers on the outside (parallel)
and time on the inside. One receiver's horizon row (128 B at N = 64) then
stays in L1 across the whole series.

**What is not stored.** Nothing per timestep. The consumer asks for the
times it needs, and gets arrays of shape (T, n) as float32 (lit as uint8),
or chunks of them from an async iterator (R13).

### R6. Lit intervals: derived, certified, for any window [TR6]

Lit intervals say *when*, not *how much* (Ola). They are derived from the
horizon and the sun's path. They are the building block of R5's `mean`
mode, and a convenience in their own right, for example "hours of direct sun
on this slope on 15 January".

- **Windows, not "days".** `lit_intervals(geometry, t0, t1)` returns, for
  each receiver, the sorted disjoint intervals of [t0, t1) in which it is lit
  (R5's lit test), as CSR: `offsets[n+1]`, `start[m]`, `end[m]`, all UTC
  ms. A "day" is a window the adapter builds (a UTC day, or local solar
  midnight to midnight). **Polar night** gives zero intervals. **Midnight
  sun** gives one interval equal to the window, or several if the horizon
  cuts the sun's circle. **Several intervals per day** fall out of the same
  code; there is no special case.
- **The function searched.**
  F(t) = min( e_app(t) + limb − h(A(t)),  asin(n·s_app(t)) ),
  in degrees. The first term is the sun against the horizon, the second the
  sun against the triangle's own plane. Lit ⇔ F > 0. Both terms are
  continuous: h is piecewise linear and continuous in A, and A is continuous
  across north (mod 360°).
- **Certified stepping, so no interval is missed.** Both terms have bounded
  rates:
  - |de_app/dt| ≤ 1.2·ω, with ω = 15°/h. |de/dt| ≤ ω, and on
    [−1.8°, 90°] refraction's slope makes 1 + dR/de lie between about 0.8
    and 1 (P3), so 1.2 is a safe margin.
  - |dh/dt| ≤ L_i·|dA/dt|, where L_i is the receiver's largest horizon slope
    between neighbouring directions, and |dA/dt| ≤ ω·(sin φ + cos φ·tan e)
    (P6 at the horizon; the tan e term bounded over the band below).
  - the angle to the plane changes at most as fast as the sun direction.
  With K_i the sum of these bounds, from any t with F(t) ≠ 0 the next
  evaluation can safely be at t + |F(t)|/K_i: no root can lie closer.
  Outside the band [h_min,i − 1°, h_max,i + 1°] the sign of the first term
  is known without the horizon, and the steps are long (a polar-night day
  costs a handful of evaluations). Near a root the steps shrink. Stepping
  stops at **1 s**, and the root is reported to within 1 s.
- **Grazing** (the sun skims a peak at local noon, F touches 0): the 1-s
  floor bounds the work. An interval is reported only where F > 0 was
  actually observed. A lit stretch shorter than the floor may be missed;
  I9 allows that and nothing longer than 2 s.
- **Cost:** estimated 50-200 evaluations of about 30 ns per receiver-day,
  so 1.5-6 µs. `@perf` measures it in TR6.
- This is standard interval root isolation with a Lipschitz bound. It is not
  claimed as new ("Novelty").

### R7. Sky-view factor and terrain configuration factor [TR4]

Static per receiver, computed once from the horizon table and the normal.

- **V_d (Dozier and Frew 1990): the fraction of isotropic diffuse sky
  irradiance that reaches the tilted surface**, relative to an unobstructed
  horizontal surface. Slope S, aspect Ā (true), horizon h_k at A_k:

  V_d = (1/N) Σ_k [ cos S · sin² H_k + sin S · cos(A_k − Ā) · (H_k − sin H_k cos H_k) ]

  with **H_k = 90° − max(h_k, h_plane,k, 0°)**, where
  h_plane,k = atan(−tan S · cos(A_k − Ā)) is the slope's own plane seen in
  direction A_k. **The clamp at 0° is required.** Without it a 30° slope with
  no terrain gives 1.0 instead of (1 + cos 30°)/2 = 0.933 (P8). Negative
  horizons (a summit) are still stored and used by the sun test in R5. The
  sky is only the upper hemisphere.
- **C_t = (1 + cos S)/2 − V_d**: the terrain configuration factor, the part
  of the surface's view that is terrain rather than sky, in the same
  cosine-weighted sense. On flat ground with horizon h all round,
  V_d = cos² h and C_t = sin² h.
- **Named in the files as `sky_view_factor_vd`**, with the definition in the
  metadata, because Zakšek et al.'s SVF (unweighted, 1 − mean sin h) is
  widely used under the same name and differs a lot: flat ground with a 20°
  horizon all round has V_d = cos² 20° = 0.883 and Zakšek's SVF =
  1 − sin 20° = 0.658.
- **Checked against brute force** as in P8, for every synthetic terrain in
  TR4's suite, within 10⁻³.
- **Why these matter most in winter** (Ola: diffuse and reflected dominate at
  low sun): on a north-facing slope steeper than the noon sun (6.6° at
  60°N at the winter solstice, P6) the beam factor b is 0 all day. What
  arrives is DHI·V_d plus the reflected term, with snow albedo about 0.8.
  C_t is not small on valley sides: 0.18 for P8's 35° slope with a 20°
  horizon.

### R8. How many horizon directions: tied to latitude, set by measurement

N is a parameter, a power of two. **The default is 64** (5.625°) until TR3's
acceptance experiment sets it. What drives the choice:
- **Rays spread with distance:** neighbouring rays are d·2π/N apart: at
  N = 64, 98 m at 1 km, 1 km at 10 km, 4.9 km at 50 km. A narrow peak between
  two rays at 10 km is missed. At low sun at high latitude, peaks at that
  distance are exactly what decides the lit intervals.
- **The error budget comes from latitude.** A horizon error δh becomes a
  timing error δh/(de/dt), and de/dt = 15°/h · cos φ · sin A falls toward the
  pole and toward winter's southern sunrises (P6). For a target timing error
  τ: δh ≤ τ · 15°/h · cos φ · sin A. At 70°N with sunrise near A = 150°
  (sin A = 0.5) and τ = 5 min, δh ≤ 0.21°. At 60°N with the same sin A:
  δh ≤ 0.31°.
- **Acceptance experiment (TR3, `@perf` with `@architect`):** on three real
  DTM10 sites (a fjord valley at 60°N, an inland valley at 65°N, a Finnmark
  site at 70°N), compute the horizon at N = 1024 as the reference, then at
  32, 64, 128 and 256. Measure the 99th-percentile horizon error after
  linear interpolation, and the lit-interval timing error over a year. The
  smallest N meeting the budget at the site's latitude becomes the default
  for that band. That is how "tied to latitude" becomes a rule rather than
  an adjective. If the answer is 128 north of 66°N, the default is chosen
  from the receivers' mean latitude.

### R9. Forcing is the consumer's; the combiner is a small, separate step [TR7]

The core returns geometry factors. **It does not know irradiance.** A
combiner, in Python and NumPy, elementwise over (T, n), turns factors and the
consumer's forcing into irradiance on each receiver:

- isotropic: G = DNI·b + DHI·V_d + ρ·G_surr·C_t
- Hay-Davies: the diffuse term becomes DHI·[(1 − A_i)·V_d + A_i·b/sin e],
  with the anisotropy index A_i = DNI/E₀ (extraterrestrial normal irradiance)
  from the consumer's forcing. The circumsolar part is thus shaded exactly
  like the beam.
- ρ (albedo of the surrounding terrain) and G_surr (irradiance on it) are
  forcing. Dozier and Frew's convention, G_surr = the receiver's unshaded
  GHI, is the default. The error in that approximation is what a radiosity
  method fixes (not in scope).
- A consumer with only GHI must split it into DNI and DHI itself (an Erbs-
  or Perez-type decomposition). That is forcing, not geometry.
- **Per unit of what area.** Every factor above is per unit of **inclined**
  surface. A model that keeps its fluxes per unit of **plan** area, as a
  gridded or cell-based snow model does, multiplies by the triangle's
  3D-to-plan area ratio 1/cos S. The adapters say which they return. The
  factor is 1.10 at 25° and 1.41 at 45°, so mixing the two up is a 10-40 %
  error on steep slopes. That is why it is written down.

### R10. Where the Python/C++ boundary sits

| crosses into C++ | as | computed in Python by |
|---|---|---|
| receiver positions (x, y) in the frame, z_obs | float64 arrays | the mesh; `raster` bilinear for point receivers |
| lat, lon | float64 degrees | `crs.reprojector(frame, EPSG:4326)` |
| J per receiver | float64 (n, 2, 2) | `crs.frame_factors` (γ, m) or the lattice frame |
| normal in ENU, exclusion radius | float64 | `radiation/frame.py` |
| near and far rasters | `RasterView` (existing) | 15a repository, `raster.to_core` |
| c = (1 − k)/(2R), N, D_near, D_far, limb | numbers | `RadiationRequest` |
| timestamps | int64 ms since the Unix epoch, UTC | `radiation/time.py` (tz-aware only) |
| ΔT, P, T | float64 | request defaults or the consumer |

**What never crosses:** a CRS, a path, a time zone, forcing. Convergence is
computed in Python and enters only inside J. That answers the brief's
question directly: yes, convergence (with scale) is computed in Python and
passed per receiver, folded into one 2×2 matrix. CLAUDE.md §2's I/O boundary
holds as written. 15c's frame rule ("every number that crosses is in the
computation frame", 15 R8) is extended by one clause: **geodetic lat/lon and
J cross as numbers describing each receiver, not as a frame the core
reasons in.**

### R11. What is stored, and how big it is

**The radiation geometry, per receiver:**

| field | type | bytes |
|---|---|---:|
| lat, lon | float64 × 2 | 16 |
| normal (E, N, U) | float32 × 3 | 12 |
| V_d, C_t | float32 × 2 | 8 |
| horizon | int16 × N | 2N |
| reach | uint8 × N | N |

That is 36 + 3N bytes: **228 B at N = 64, so 2.3 GB for 10 M triangles**
(1.3 GB at N = 32, 4.2 GB at N = 128). A time series of one float32 factor
at 10 M × 8760 h is 350 GB, and three factors are 1 TB. So the geometry is
about 150-800× smaller than one year of what it generates, and it serves any
year and any time step. At Shyft's scale (10⁵ cells, N = 64) it is 23 MB.

Built on load, not stored: the E/N/U basis (36 B per receiver as float32),
which is cheaper to recompute than to read.

- **Format: `.npz`** (NumPy's zip of arrays, no new dependency), with a JSON
  metadata member: mesh identity (a hash of vertices and triangles), DEM
  provenance (tile names, far-field source, D_near, D_far, block size), N,
  k, R, the SPA and ΔT model versions, the frame (CRS label or lattice), and
  the V_d definition. A geometry refuses to evaluate against a mesh whose
  hash differs. Alternative: HDF5 (Q5).
- **Also written into the `.vtk`** (increment 13), when asked: slope,
  aspect, V_d, C_t and the mean horizon as cell fields, for looking at in
  ParaView. Not the horizon table.

### R12. DOM10: a hook, not a design

Ola: "Later, we will also use the DOM10 for this." The geometry build takes
the surface to search as a **named input**: `surface = "dtm" | "dom"`, one
repository each. A second geometry, searched over DOM10, sits beside the
first, and the evaluator does not care which it has. What waits for its own
design: whether a DTM receiver under canopy should see the DOM at all (at
10 m, a crown next to the receiver becomes a 45° horizon, and what is needed
is a transmissivity, not a horizon), and receivers on the canopy or roof
surface. Whether DOM10 shares DTM10's lattice is to verify. Hydrology keeps
meshing the DTM (ROADMAP).

### R13. Adapters [TR7]

Thin, pure Python over R5 and R6. None of them adds C++.
- **ML tensors:** `factor_tensor(geometry, times, fields, groups=None) ->
  (array (T, G, F) float32, metadata)`, where `groups` optionally aggregates
  receivers (sub-catchments, elevation bands) with plan-area weights. Also
  an async variant, `aiter_factor_tensor(...)`, that yields time chunks
  through `asyncio.to_thread`, with the GIL released in C++, so a GUI or an
  API worker is not blocked.
- **Site polygons (solar parks, crops):** `SiteSet` from shapely polygons in
  their own CRS, reprojected by `crs.reprojector` (as 15b does for domains).
  Each site gets its triangles with plan-area intersection weights, and the
  aggregation is area-weighted.
- **Panels:** `panel_receivers(points, tilt, azimuth_true, height)`
  builds point receivers whose normal is the panel's (plane of array),
  with z_obs = DEM + height. Horizon and evaluator are unchanged. A solar
  park is then a receiver set with no mesh at all.
- **Shyft:** waits on Q6. Shyft needs either per-cell factor series or the
  static geometry for its own radiation routine, and which one depends on
  its API.

## Ola's requirements, one by one

| requirement | where | in one line |
|---|---|---|
| 1. horizon per direction is the primary data | R4, R11 | int16 × N per receiver, in true azimuth |
| 1. the sun's path followed through time | R5, R6 | azimuth and elevation per timestamp against the horizon at that azimuth |
| 1. lit intervals derived, not primary | R6 | from the horizon and the path; they give when, R5 gives how much |
| 1. azimuthal resolution tied to latitude | R8 | the δh budget scales with cos φ; N set by a measured acceptance run |
| 1. atmospheric refraction | R2, R4 | astronomical for the sun (eq. 42, extended to −1.8°), terrestrial k for the terrain ray |
| 1. polar night, midnight sun, several intervals a day | R6 | windows, not days; zero, one or many intervals from the same code |
| 1. diffuse and reflected dominate at low sun | R7, R9 | V_d and C_t, definitions checked against brute force (P8) |
| 2. per-triangle local frame, grid convergence | R3 | J = (1/m)·Rot(γ), γ from pyproj, sign pinned (P4) |
| 2. local vertical | R3 | ellipsoid normal at the receiver's lat/lon; deflection ignored, stated |
| 2. curvature and refraction in the horizon search | R4 | c·d² with c = (1−k)/(2R), in the reference definition itself |
| 2. search tens of km beyond the catchment | R4 | near raster to 20 km at full resolution, far raster to 200 km at 100 m, through 15a's repository |
| 2. UTM scale error | R3 | inside J; ≤ 0.4 % on tan S |
| 2. the geographic-DEM basin (15c) | R3 | the lattice frame's J; nothing else changes |
| 3. size | R5, R11 | nothing per timestep is stored; 228 B per receiver at N = 64; about 29 ns per receiver-timestep |
| legacy `solar_position.h` | Legacy, R2 | carried, verified against NREL (P1), ΔT and refraction guard fixed |
| legacy ray-cast shading | Legacy, R4 | idea carried, CGAL AABB replaced by horizon tables |
| horizon on DEM or TIN | R4 | DEM, per receiver, with the triangle as its own near field |
| std::chrono only, UTC | R2 | `sys_time<milliseconds>`, naive times refused |
| the Python/C++ boundary | R10 | convergence in Python, folded into J |

## Invariants

- **I1.** The reference path (R2) reproduces NREL's published test values to
  the published digits, and its periodic-term tables equal the paper's, term
  by term.
- **I2.** Every timestamp is UTC. A naive `datetime` is refused.
  Conversions between `datetime64[ms]`, `sys_time<milliseconds>` and Julian
  day are exact to the millisecond.
- **I3.** A_grid = A_true − γ with pyproj's sign. Pinned through J at Solund
  and at Vardø (P4).
- **I4. The pruned horizon search equals the exhaustive one bit for bit**,
  for any thread count (R4). *Invariant-critical.*
- **I5.** The horizon never decreases when D_near or D_far grows, and never
  increases when c grows. Adding a constant to both the terrain and z_obs
  leaves it unchanged.
- **I6.** 0 ≤ V_d ≤ 1, 0 ≤ C_t ≤ 1, V_d + C_t = (1 + cos S)/2, and V_d
  matches brute-force hemisphere integration within 10⁻³.
- **I7.** The fast evaluator agrees with the reference path within 0.003° in
  elevation and azimuth. Its lit status agrees except where the sun is within
  0.003° of the horizon.
- **I8.** All outputs are independent of the thread count.
- **I9. Lit intervals are complete and tight.** At every 1-s sample of the
  window, the lit test agrees with interval membership, except within 1 s of
  an endpoint. No lit stretch longer than 2 s is missed. *Invariant-critical.*
- **I10.** In `mean` mode, the mean of b over a step equals a 1-s brute-force
  average within 10⁻⁴.
- **I11.** A geometry refuses to evaluate against a mesh whose hash differs
  from its own.
- **I12.** No CRS, path or time zone enters C++. This is enforced by the
  existing structure and reviewed, not by a new gate.

## Degeneracy policy

- **The sun exactly on the horizon** (F = 0): unlit. Lit is strict, F > 0.
- **Azimuth exactly on a stored direction:** that value, with no
  interpolation.
- **A ray that leaves coverage or meets NoData:** no obstacle past that
  point. `reach` records the distance (R4). A receiver whose every sample is
  missing has h = −90° and reach 0 in that direction.
- **A zero-area triangle** has no normal. Refused at build, naming the
  triangle. The CDT does not produce one; the guard is cheap.
- **A vertical triangle** (n_U = 0): allowed, slope 90°.
- **|lat| > 89.9°:** refused, since the ENU basis is undefined at the pole.
  Svalbard (≤ 81°N) is inside.
- **A window with t1 ≤ t0:** refused. **No timestamps:** empty arrays.
- **A negative observer offset:** refused. **A point receiver outside the
  near raster:** refused.
- **Timestamps far from the present:** allowed. SPA holds for years −2000 to
  6000. The ΔT default extrapolates by the /100 parabola, and the metadata
  records which ΔT was used.

## Verification plan

1. **The sun (TR1).** NREL's test case (P1's values). SPA's published
   sunrise, transit and sunset for the same case (to verify: 06:12:43,
   11:46:04 and 17:20:19 local time, from memory). The term tables diffed
   against the paper. If Q2 is yes, 10⁴ random (time, lat, lon, height)
   against pvlib's independent SPA within 10⁻⁶°. ΔT at 2026 against the
   legacy defect. The refraction guard at −1.8° and below.
2. **The frame (TR2).** J against a finite-difference transform (P4's method)
   at Solund, Oslo, the central meridian and Vardø. A synthetic plane of known
   true slope and aspect, placed off the central meridian, gives the right
   true aspect (this catches the legacy's missing convergence).
3. **Horizons on synthetic analytic terrains (TR3a).** Each terrain is
   piecewise linear on the grid, so bilinear sampling is exact and the closed
   form is the exact reference:
   - **V-valley**: axis on a grid line, side slopes β, ridges at half-width
     W. For a receiver on the axis, tan h(A) = |sin A|·tan β − c·r_i. With
     c = 0 that is |sin A|·tan β exactly.
   - **Square crater with planar ramps**: floor to W, ramp to height H at
     W + L. Away from the corner diagonals, tan h(A) = (H − c·d²)/d with
     d = (W + L)/max(|sin A|, |cos A|). This holds while c < H·W/(L·d²),
     that is, while the ramp's rising tangent outruns the curvature term.
     Test dimensions (H = 500 m, W = L = 1 km) satisfy it by three orders
     of magnitude.
   - **Isolated cone**, apex on a node and at a sample distance r₀ from the
     receiver: tan h toward the apex = (H − c·r₀²)/r₀. This checks the
     curvature and refraction term at 20, 50 and 100 km against P5's table.
   - **Flat plane**: horizon = −atan(c·r_i). The sun test then reproduces
     astronomical sunrise and sunset (with `limb` = upper and the standard
     atmosphere) against the reference path's own sunrise.
   - The same terrains rotated by a J with γ = ±15° (the Vardø case), so
     that a missing convergence shows up as a 15° shift in h(A).
4. **Sky factors (TR4).** V_d and C_t against brute force (P8's method) on
   each synthetic terrain.
5. **Lit intervals (TR6).** Against 1-s brute sampling on adversarial
   horizon profiles: combs, a horizon that touches the sun's noon elevation
   (grazing), polar night (70°N, 21 Dec), midnight sun (70°N, 21 Jun, a
   north notch giving two intervals), and a day with three intervals. Closed
   form for a constant horizon h₀ with declination held fixed:
   cos H₀ = (sin h₀ − sin φ sin δ)/(cos φ cos δ).
6. **A Norwegian site with known winter sun (TR3b acceptance). The source
   must be found; these are candidates to verify:**
   - **Rjukan** (59.88°N, 8.59°E): the town centre is in the mountains'
     shadow from about late September to about mid-March, which is why the
     sun mirrors (Solspeilet, 2013) were built. The documented dates and the
     exact spot the mirrors light are needed.
   - **Tromsø, "Soldagen" 21 January:** the day the sun is first seen from
     the town over the mountains, against an astronomical return around 15
     January. The difference is terrain. The documented viewpoint is
     needed.
   - **Measured radiation**: a pyranometer station in a steep valley with
     10-min or hourly global radiation (NIBIO's agrometeorological network
     or MET Norway's Frost API; to find). On clear days, first light at the
     station against predicted lit intervals is an objective, dated test.
     This is the strongest candidate, and the one to use for any
     publication.
   - Longyearbyen's "Solfestuka" (the sun reaches the old hospital steps
     around 8 March against an astronomical return around 16 February) is
     the classic, but Svalbard is outside DTM10 and needs Norsk
     Polarinstitutt's DEM.
7. **External comparator (optional, Q7).** PVGIS's horizon profile for a
   handful of points, compared within its own DEM's resolution.

## Increments and LOC

Estimates against the 700-line ceiling (CLAUDE.md §2: production code only;
comments, docstrings and raw-literal bodies do not count; tests excluded).
**SPA's coefficient tables are counted.** They are `constexpr` array rows,
not raw literals.

| # | what | depends on | C++ | bind | Python | total |
|---|---|---|---:|---:|---:|---:|
| TR1 | time and the sun: `time.hpp`, `solar_position.hpp` (tables ≈ 260 rows: L 129, B 7, R 59, nutation 63), ΔT, refraction guard; `radiation/time.py` (UTC only) | nothing | 440 | 40 | 60 | **≈ 540** |
| TR2 | the receiver frame: `receivers.hpp` (SoA views, ENU basis); `radiation/frame.py` (centroids, lat/lon, J, normals, exclusion radius), `crs.frame_factors`, `RadiationRequest` | 15b (shipped) | 80 | 40 | 230 | **≈ 350** |
| TR3a | the horizon kernel: `horizon.hpp` (max-pyramid, exhaustive reference, pruned search, curvature, two rasters, reach) | TR2 | 300 | 50 | 0 | **≈ 350** |
| TR3b | the horizon's inputs: `radiation/search.py` (near mosaic through 15a, far grid by streamed block means), CLI `rasputin radiation geometry`; N acceptance run (R8) | TR3a, 15a (shipped) | 0 | 0 | 320 | **≈ 320** |
| TR4 | sky factors and the geometry file: `sky.hpp`; `radiation/geometry.py` (.npz, metadata, mesh hash); `.vtk` cell fields | TR3a | 60 | 20 | 200 | **≈ 280** |
| TR5 | the evaluator: `evaluate.hpp` (ephemeris, `instant`), bindings with the GIL released, `radiation/evaluate.py` (sync and async), CLI | TR1, TR4 | 220 | 60 | 160 | **≈ 440** |
| TR6 | lit intervals and `mean` mode: `lit_intervals.hpp` (certified stepping, CSR), Gauss-Legendre means | TR5 | 280 | 40 | 60 | **≈ 380** |
| TR7a | combiner (isotropic, Hay-Davies), ML tensors, groups, async iterator | TR5, TR6 | 0 | 0 | 280 | **≈ 280** |
| TR7b | site polygons and panel receivers | TR7a | 0 | 0 | 250 | **≈ 250** |
| TR7c | Shyft adapter | Q6 | | | | to size after Q6 |


**Built by increment 46** (`docs/increments/46-hbv-snow.md`, its section 10), to this design and not as a second version:

- TR1: all of it, in increment 46's PR 2 (built).
- TR2: all of it, in increment 46's PR 4 (built).
- TR3a: all of it, in increment 46's PR 5 (built), with the far raster optional (absent means
  D_far = D_near) and the pins of `46-hbv-snow.md` section 10.8 (reach found
  from coverage before the march, floored to whole km; the early stop over
  both rasters; far samples from D_near + Δ_far; int16 rounded half away
  from zero).
- TR3b: only the near raster's call in `radiation/search.py`, in PR 6; no
  far raster, no `rasputin radiation geometry` CLI, no N-by-latitude run.
- TR5: the C++ `instant` evaluator, its binding and a synchronous Python
  call, in PR 6, reading the receiver set and horizon table as arrays
  before TR4 exists; no CLI, no async call, no `mean` mode.
  `RadiationGeometry` gets only TR2's and TR3a's fields; TR4 adds the sky
  factors, the `.npz` file, TR5's CLI and its async call to the same class
  and module.

PRs 4, 5 and 6 each mark their line built when they land.

About 3200 lines in all. The largest is TR1 at about 540. If the review
counts differently and TR1 goes over, the tables move into their own PR
(TR1a: tables with a transcription test; TR1b: the equations).

- **Order:** TR1 → TR2 → TR3a → TR3b → TR4 → TR5 → TR6 → TR7. TR1 depends
  on nothing and could be built now (Q9). ROADMAP places the rest after
  auto-catchment.
- **The lattice frame's J** (R3) arrives with 15c or after it. TR2 does the
  projected case only, and its J interface is the same.
- **Invariant-critical suites, mutation-tested:** **TR3a** (I4, pruned
  equals exhaustive, plus the closed forms) and **TR6** (I9, complete and
  tight intervals). TR1's tables are data. A second, independent
  implementation (Q2) and the term-by-term diff test them better than
  mutants would.
- **Template cross products:** none. The rasters are `RasterView<float>`,
  the only instantiation used.
- **`@perf`:** this touches no refine or mesh code, so rule 2's acceptance
  does not apply. But R4, R5 and R6 quote cost estimates. `@perf` measures
  TR3a, TR5 and TR6 on a real mesh, and the measurements replace the
  estimates in this file before anyone quotes them elsewhere.
- **`@migration-expert`:** one short pass before TR1's red step (Legacy).

## Not in scope

- **Irradiance itself**: clear-sky models, transmittance, splitting GHI.
  That is forcing, and Ola ruled forcing to be the consumer's.
- **Multiple reflections (radiosity, Helbig et al. 2009).** C_t with one
  bounce is the standard hydrological approximation. A view-factor matrix is
  O(n²) before pruning. It waits until a consumer shows that one bounce is
  not enough.
- **Longwave radiation** (terrain emission, sky longwave × V_d). V_d is the
  factor it would need, so the door is open. Nobody has asked yet.
- **Canopy transmissivity and DOM10's own design** (R12).
- **A horizon raster for every DEM node** (Dozier-Frew sweep, R4's
  rejected alternative). This is for gridded products; nobody has asked.
- **Partial shading of a large triangle** (a fraction lit, not a bool). One
  horizon per receiver gives one answer per triangle. A 1 km lake triangle
  that is half in shadow reads as lit or unlit at its centroid. Sub-sampling
  receivers per triangle, several points with area weights, would fix it
  with no new machinery, at k× the cost. Not designed until a consumer needs
  it.

## Facts to verify

No web access in this round. Each item says what to check.

1. **Reda and Andreas 2004** (Solar Energy 76(5):577-589; NREL/TP-560-34302,
   2008 revision): the citation; the test values quoted in Legacy (JD, L, B,
   R, zenith 50.11162°, azimuth 194.34024°, incidence 25.18700°); the
   sunrise, transit and sunset for the test case (06:12:43, 11:46:04,
   17:20:19 local, from memory); the ±0.0003° uncertainty and its year
   range; and that eq. 42 is applied only for e0 ≥ −(0.26667° + 0.5667°).
2. **Table numbers** A4.2 (L, B, R terms) and A4.3 (nutation) in the paper.
3. **Espenak and Meeus 2006** ΔT polynomials (NASA/TP-2006-214141): the
   2005-2050 coefficients 62.92, 0.32217, 0.005589; the long-term parabola
   −20 + 32·((y − 1820)/100)². Also the IERS-measured ΔT for 2026 (≈ 69 s, from
   memory).
4. **Bennett 1982** as the source of SPA's eq. 42 (J. Navigation 35:255-259).
5. **Dozier, Bruno and Downey 1981** (Computers & Geosciences 7(2):145-151)
   and **Dozier and Frew 1990** (IEEE TGRS 28(5):963-969): citations; **the
   V_d formula in R7 exactly as printed** (P8 confirms the recalled formula
   against brute force, with the clamp at 0°, but not that the paper prints
   it this way); C_t ≈ (1 + cos S)/2 − V_d.
6. **Stewart 1998** (IEEE TVCG 4(1):82-93): citation, and its complexity
   claim.
7. **Corripio 2003** (Int. J. GIS 17(1):1-23).
8. **Marsh, Pomeroy and Spiteri 2012** (Hydrol. Process. 26(12):1767-1778)
   and **Marsh et al. 2020, CHM** (GMD 13:225-247): citations; **which
   shading method they use** (horizon table or per-timestep); whether they
   search beyond the mesh. Closest prior art, so it matters for any
   publication.
9. **Montero et al. 2009** (Solar Energy 83(7):998-1012): authors, pages and
   method.
10. **Hofierka and Šúri 2002** (GRASS `r.sun`); **PVGIS's horizon**: which
    DEM it uses north of 60°N, and its azimuth step.
11. **Zakšek, Oštir and Kokalj 2011** (Remote Sensing 3(2):398-415): the SVF
    definition quoted in R7.
12. **Hay and Davies 1980; Perez et al. 1990** (Solar Energy 44(5):271-289);
    **Helbig, Löwe and Lehning 2009** (J. Atmos. Sci. 66(9):2900-2912);
    **Allen, Trezza and Tasumi 2006** (Agric. For. Meteorol. 139:55-73).
13. **Hirt et al. 2010** (J. Geophys. Res. 115:D21102) on the variability of
    k near the ground; k = 0.13 as the conventional value.
14. **Shyft**: that its radiation routine follows Allen et al. 2006; what its
    cell data carries (slope, aspect, a radiation slope factor?); whether its
    fluxes are per plan area; its input API for radiation (Q6).
15. **Rjukan**: the documented dates when the town centre loses and regains
    direct sun, and the spot. **Tromsø Soldagen** (21 January) and the
    viewpoint it refers to. **Longyearbyen**: 8 March at the old hospital
    steps, 16 February astronomical.
16. **A pyranometer station** in a steep Norwegian valley with sub-hourly
    global radiation (NIBIO LMT, MET Norway Frost API).
17. **DOM10** shares DTM10's lattice and tiling (R12).
18. **DTM10 at sea and at the border**: NoData, 0, or absent tiles (R4's
    `reach` depends on it; this can be read from Ola's archive with 15a's
    probe, without the web).

## Questions for Ola (ruled 2026-09-28, see "Ruled by Ola")

Each has a recommendation. None blocks TR1.

- **Q1. Horizon directions.** Default N = 64 until TR3b's acceptance run
  (R8) measures what each latitude band needs, then a default by band.
  *Recommend: yes.* At N = 64, 10 M triangles take 2.3 GB of geometry.
- **Q2. pvlib as a test-only dependency**, as an independent second SPA
  implementation over 10⁴ random cases. It is BSD-licensed, pure Python plus
  NumPy, and is not on the prohibited list. NREL's own C code has a
  restrictive licence and should not be vendored. *Recommend: yes, tests
  only.*
- **Q3. The refraction coefficient k.** A constant 0.13 (the geodetic
  convention) with a request override. Winter inversions over snow can push
  k well above that, and a per-season or per-timestep k would make the
  horizon depend on time. *Recommend: a constant, overridable per geometry;
  revisit only if the station comparison (verification 6) shows a bias at
  low sun.*
- **Q4. Beyond Kartverket's coverage** (Sweden, Finland, Russia, the sea).
  (a) record `reach` and treat the rest as no obstacle, or (b) also fill the
  far raster from Copernicus GLO-90 (EPSG:4326, resampled into the far grid,
  which R4 allows for the far field only). *Recommend: (a) in TR3b, with the
  far raster's source a named input so (b) is additive later.* Border
  catchments in Troms, Finnmark and Hedmark will see a truncated eastern
  horizon until (b).
- **Q5. Geometry file format.** `.npz` (no new dependency) or HDF5 (h5py; the
  legacy used HDF5/XDMF for shades). *Recommend: .npz*, with the `.vtk`
  cell fields for viewing.
- **Q6. Shyft.** Which does Shyft want: (a) per-cell factor time series from
  our evaluator, fed into its energy balance, or (b) the static geometry
  (slope, aspect, V_d, horizon) for its own radiation routine to use? Which
  Shyft version and API? *Recommend: (a) first. It keeps the horizon logic
  in one place, and Shyft's own routine does not shade by terrain, as far
  as is known (to verify).*
- **Q7. PVGIS as an external comparator** in an acceptance probe (network,
  not CI). *Recommend: yes, for the solar-park use.*
- **Q8. The Norwegian verification site.** *Recommend: a pyranometer station
  in a steep valley first (objective, dated), with Rjukan as the famous
  second.* Do you know a station, or should `@architect` search MET's and
  NIBIO's station lists?
- **Q9. Build TR1 early?** It depends on nothing, fixes a legacy defect, and
  is useful on its own (sun position per timestamp for Shyft and ML).
  *Recommend: yes, as a small PR whenever there is a gap in the Norwegian
  work; the rest waits for auto-catchment as ROADMAP says.*
- **Q10. The solar disc at the horizon.** The lit test uses the sun's centre
  by default, with the upper limb as an option (R5). A fractional disc,
  where lit goes from 0 to 1 over 0.53°, is more physical and costs one
  more interpolation, but it makes "lit" a fraction instead of a bool.
  *Recommend: centre, bool, for now.* The disc matters for minutes a day,
  and the station data in verification 6 will say whether it is worth it.
