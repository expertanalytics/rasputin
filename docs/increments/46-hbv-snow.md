# Increment 46: HBV snow on the triangles, melt corrected for the sun, and a melt animation

**Status: designed; design review round 1 answered** (2026-10-13,
`@architect`, branch `worktree-snow-melt`). Design only: no code, no tests.
Three PRs (section 8). Ola was away when it was written; every default taken
without him is listed in section 9 as one of seven questions with its
default. Every measured figure comes from the saved probe,
`docs/increments/46-probes/probe_spring.py`, whose output is
`46-probes/probe_spring.out`, except where a section names another source.

## 1. What Ola asked

Ola, 2026-10-13, picking the idea: "schedule the implementation of the
HBV-snow routine to live on our triangles. My goal is to get an animation of
snow melt, first on synthetic data (eg by assuming uniform snow cover and slope
dependent melting based on the sun). See how well this will look for a quite
fine mesh."

His rulings on the idea the same day, from the private ideas file:

- the sun's correction is available **both** through the threshold temperature
  and through the melt factor, so the two can be compared ("I think both
  approaches could be interesting");
- a prototype candidate for Shyft, written to rasputin's standards, not
  Shyft's;
- wind-distributed snowfall is a later increment, not this one.

What this increment delivers, in one line: `rasputin snow --mesh mesh.vtk
--out DIR` runs an HBV snow routine on every triangle of a rasputin mesh
through a synthetic spring, with the melt corrected for the sun by each
triangle's slope and aspect, and writes an animation ParaView plays.

## 2. Prior art: legacy and literature

Written before the design.

### 2.1 Literature

**The HBV snow routine.**
- Bergström, S. (1976). *Development and application of a conceptual runoff
  model for Scandinavian catchments.* SMHI Reports RHO 7. The origin of HBV.
- Lindström, G., Johansson, B., Persson, M., Gardelin, M., Bergström, S.
  (1997). *Development and test of the distributed HBV-96 hydrological model.*
  J. Hydrology 201:272-288. The routine used here: precipitation is snow below
  a threshold temperature; degree-day melt above a threshold; refreezing of
  the meltwater the pack holds at a fraction of the degree-day factor; the pack
  holds liquid water up to a fraction of its solid part and releases the rest.
- Shyft's HBV snow, read from its source on 2026-10-13
  (`gitlab.com/shyft-os/shyft`, `cpp/shyft/hydrology/methods/hbv/snow/`,
  `calculator.h`, `parameter.h`, `state.h`). Its parameters are `tx` (rain or
  snow), `ts` (melt onset), `cx` (degree-day factor), `lw` (liquid water
  fraction) and `cfr` (refreeze: `cfr·cx·(ts − t)`), plus a snow distribution:
  redistribution factors `s` at quantile points `intervals` (default five
  points, all factors 1, so uniform). The per-point update is HBV-96's.
  **What differs here:** (a) the same five parameters under the same names,
  so a Shyft port is a rename-free move; (b) the distribution is equal-area
  bins with a constant factor each, not factors at quantile points with
  piecewise-linear interpolation, because Shyft's interpolated snow-covered
  area (SCA) needs a fix-up after every snowfall and its own comment says it
  can break the mass balance; bins balance mass exactly (section 4.2, I1);
  (c) the parameter defaults differ (section 5): Shyft's `cx = 1` and
  `cfr = 0.5` are placeholders, not values for Norway.
- Sælthun, N. R. (1996). *The "Nordic" HBV model.* NVE Publication 7/1996.
  Snowfall log-normally distributed within a unit, with a fixed coefficient of
  variation (CV), which gives a snow-covered area that falls continuously as
  the pack melts. The CV form of the distribution used here is this one. (Its
  section on the distribution was not read; the description is from papers
  citing it, found 2026-10-13. To verify before a publication.)
- Skaugen, T. (2004, Annals of Glaciology 38) argues that a distribution fixed
  in time misstates melt floods. Noted, not adopted: a fixed distribution is
  what HBV and Shyft's HBV use.
- Liston, G. E. (2004). *Representing subgrid snow cover heterogeneities in
  regional and global models.* J. Climate 17:1381-1397. A log-normal
  distribution with a closed-form depletion curve (SCA from mean, CV and
  accumulated melt). **Considered and rejected** for this increment: the
  closed form holds for one accumulation followed by melt, and a spring with
  snowfall during melt breaks it. The CV values it tabulates by terrain type
  were not read (to verify); the default CV here is a chosen value.

**Melt corrected for the sun.**
- Hock, R. (1999). *A distributed temperature-index ice- and snowmelt model
  including potential direct solar radiation.* J. Glaciology 45(149):101-111.
  Melt = (MF + a·I)·T for T above a threshold, where I is the potential
  clear-sky direct solar radiation on the slope, computed from the sun's
  position, slope, aspect and topographic shading, with I = I0·(Rm/R)²·
  Ψa^(P/(P0·cos Z))·cos θ, at 15-minute steps aggregated. (Formula as given by
  the hydrobricks documentation's account of Hock 1999, read 2026-10-13; the
  paper itself was not read. To verify before a publication.) **This is the
  melt-factor form, built on as is**, with three departures, each named:
  1. **The anomaly, not the absolute value.** The correction uses
     ΔI = I(slope) − I(horizontal), not I. A horizontal triangle then melts
     exactly as plain HBV, in both forms, so `cx` keeps its HBV meaning and
     the comparison Ola asked for is between slopes only. Algebraically this
     is Hock's model with a base factor MF(t) = cx − a·I_h(t) that varies with
     the season; Hock's MF is constant. Named, not hidden.
  2. **No cast shadows.** Self-shading only (a slope facing away from the sun
     gets nothing); shadows cast by other terrain need horizons, which the
     terrain radiation design computes (its TR3). The synthetic mountain is a
     single massif with slopes up to 51°, so cast shadows matter only at low
     sun; named as the first thing the radiation work adds.
  3. **One pressure for the domain.** Ψa's exponent uses P/P0 = 1 at the
     reference height, not each triangle's pressure. At the probe's summit
     (1924 m) the direct beam would be about 12 % stronger at 30° sun
     (Ψa^−0.4 with Ψa = 0.75). Named; the per-triangle pressure is one
     array more when wanted.
- **The threshold form** (Ola's original idea: the sun lowers the melt
  threshold on exposed slopes). Melt = cx·(T − (ts − d·ΔI)). Rearranged,
  melt = cx·(T − ts) + cx·d·ΔI: a degree-day term plus a radiation term that
  does not scale with temperature. That is the shape of Pellicciotti et al.'s
  enhanced temperature-index model (Pellicciotti, F. et al. 2005, J.
  Glaciology 51(175):573-587, melt = TF·T + SRF·(1 − α)·G, cited from memory,
  to verify), with the difference that here the radiation term only acts while
  T is above the shifted threshold. No published HBV variant that shifts the
  melt threshold by radiation was found in the searches made (below); it is
  not claimed as new either, because the search was short.
- Cazorzi, F., Dalla Fontana, G. (1996). *Snowmelt modelling by combining air
  temperature and a distributed radiation index.* J. Hydrology 181:169-187
  (from memory, to verify). The other classic radiation-index melt model; it
  scales the melt factor by a radiation index relative to a reference, which
  is close to departure 1 above.
- Hock, R. (2003). *Temperature index melt modelling in mountain areas.* J.
  Hydrology 282:104-115. The review that places all of the above.

**Snow models on triangulated terrain.**
- Marsh, C. B. et al. (2020). *The Canadian Hydrological Model (CHM) v1.0.*
  Geosci. Model Dev. 13:225-247. Snow processes per triangle of a
  variable-resolution TIN, with shading. The closest prior art to running snow
  on rasputin's triangles.
- tRIBS (Rinehart, Vivoni and Brooks 2008, Ecohydrology 1:253-270, from
  memory, to verify) runs a snow model on a TIN with radiation sheltering.

**Novelty: none claimed.** A snow model on a TIN exists (CHM, tRIBS);
radiation-corrected degree-day melt exists (Hock 1999 and after). Searches
made 2026-10-13: "Hock 1999 distributed temperature-index ... potential
direct solar radiation"; "shyft hydrology hbv_snow parameter tx cx ts lw cfr";
"Sælthun 1996 Nordic HBV model snow distribution lognormal"; "Liston 2004
subgrid snow CV". Before any publication: search "HBV threshold temperature
radiation aspect", "degree-day threshold shift solar radiation", and CHM's
snow papers.

### 2.2 Legacy

```sh
$ git grep -l -i -E "snow|hbv|melt|degree.day|solar|sun|radiation" legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/calculate.py
legacy-archive:legacy/rasputin/geo_tiff_reader.py
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/rasputin/reader.py
legacy-archive:legacy/rasputin/solar_position.h
legacy-archive:legacy/rasputin/triangulate_dem.h
```

- **No snow model.** The `snow` hits in `globcov_repository.py` and
  `gml_repository.py` are land-cover class names ("Permanent snow and ice",
  `glacier_and_snow`). Nothing is carried.
- **The sun's position**, `solar_position.h` (752 lines, Ola's SPA port,
  bound in `legacy/bindings.cpp:379-392`). Carried by the terrain radiation design's
  TR1, which this increment reuses (section 4.4); not carried twice.
- **Shading** (`legacy/rasputin/triangulate_dem.h:600-672`, `legacy/rasputin/calculate.py:9-11`,
  `legacy/rasputin/geo_tiff_reader.py:22-71`, `legacy/rasputin/reader.py:471`): a face is shaded when its
  normal faces away from the sun, or when a CGAL ray from its centroid hits
  another face. The first test is this increment's self-shading,
  `max(0, n·s)`. The ray cast is CGAL and is not carried; cast shadows come
  with the radiation design's horizons.
- **Slope** (`legacy/rasputin/triangulate_dem.h:757-758`): `atan2(|n_xy|, n_z)`, the same
  definition as here.

## 3. The ideas file

Read `../rasputin_scratch/ideas/ideas.md` before the design.

- **Taken up:** "Snowpack on the TIN", parts 1 (HBV snow per triangle,
  continuous SCA), 3 (melt corrected for the sun) and 4 (a synthetic spring
  animated), as Ola's pick says.
- **Overlaps, reused not duplicated:** the terrain radiation design (branch
  `design-terrain-radiation`, `docs/increments/terrain-radiation.md`, ruled
  2026-09-28). Section 4.4 says which part.
- **Left out, handed back:** part 2 (wind-distributed snowfall, the next snow
  increment); "Dual-use applications" and "Intent layer" (no overlap).

## 4. The design

### 4.1 The blueprint

```
 rasputin synthetic-terrain --out synth.tif          (new; PR 3)
     snow/terrain.py  mountain(spec) -> float32 grid    pure numpy
     io/geotiff_write.py  encode_geotiff(grid, origin, step, epsg) -> bytes
 rasputin mesh --dem synth.tif --tolerance 1 --out mesh.vtk   (exists, unchanged)
 rasputin snow --mesh mesh.vtk --out DIR [--sun-correction none|melt-factor|threshold] ...
     cli.snow (composition root; the only place a path exists)
       │  flags -> SnowJob (frozen Pydantic: HbvSnowParameters, SnowDistribution,
       │          SunCorrection, SyntheticSpring, FrameSpec, Resources)
       ├─ io/vtk_read.py  read_mesh(bytes) -> MeshArrays(vertices, triangles, crs)
       ├─ snow/geometry.py  triangle_geometry(mesh, convergence) -> normals (ENU),
       │                    centroid z, horizontal area                 pure numpy
       ├─ crs.py  centre (lat, lon) and grid convergence at the centre  pyproj
       ├─ snow/forcing.py  SyntheticSpring -> per step: temperature per triangle,
       │                   precipitation per triangle, S sun directions (ENU),
       │                   S direct-beam values           numpy + _core sun (TR1)
       ├─ snow/run.py  async run(job, geometry) -> AsyncIterator[Frame]
       │      per step:  _core.radiation_anomaly(normals, sun, beam) -> ΔI (n)
       │                 _core.SnowState.step(params, forcing, threads)     GIL released
       │      per frame: copy swe, sca out of the state
       └─ io/ensight.py  pure encoders: geometry bytes, variable bytes, case text
 C++ (include/terrain/snow/, header-only; no file, no path, no CRS, no time)
     hbv_snow.hpp         parameters, distribution, state, step
     radiation_index.hpp  ΔI per triangle from normals, sun directions, beam
 bindings/snow.cpp        bind_snow(py::module_&), called from bindings/core.cpp
```

- **Data against execution.** Every parameter is a frozen Pydantic model
  (`snow/params.py`); the C++ receives plain numbers and arrays. The CLI only
  turns flags into a `SnowJob`.
- **State.** The only mutable thing is the `SnowState` object Python holds;
  the C++ step writes only its own triangle's slots, so it runs in parallel
  over triangles through `parallel_util::for_each_block` with the run's
  `--threads` (increment 38's `Resources`). No global state, no time inside
  C++: the step sees a length in seconds, never a date.
- **Swappable pieces.** The HBV kernel never sees the sun: it takes ΔI as an
  array, empty when the correction is off. The radiation kernel never sees
  time or the date: it takes sun directions. The sun's position is TR1's.
  When the radiation design's horizons arrive (its TR3-TR5), ΔI gains a lit
  factor and the HBV kernel does not change.
- **The I/O boundary holds** (`CLAUDE.md` §2): the C++ gets arrays from
  Python. Files are read and written only in `cli.py`.
- **Async-ready.** `run` is an async generator; each frame's batch of steps
  runs in `asyncio.to_thread`, and the binding releases the GIL, so a GUI or
  an API worker can drive it.
- **Testable alone.** The HBV kernel is tested with no mesh (arrays of
  triangles that need not exist); the radiation kernel on hand-made normals;
  the writers on small arrays; the run on a 2-triangle mesh.

### 4.2 The HBV kernel (`include/terrain/snow/hbv_snow.hpp`)

```cpp
namespace terrain::snow {
struct HbvSnowParameters {  // Shyft's names
    double tx;   // °C: precipitation falls as snow below it
    double ts;   // °C: melt threshold
    double cx;   // mm per °C per day: degree-day factor, >= 0
    double lw;   // liquid water the pack holds, fraction of its solid part, >= 0
    double cfr;  // refreeze coefficient: refreeze = cfr·cx·dt·(ts − T), >= 0
};
enum class SunCorrection : std::uint8_t { none, melt_factor, threshold };
struct SunCorrectionParameters {
    SunCorrection mode;
    double r;    // melt_factor: mm per °C per day per W/m²
    double d;    // threshold: °C per W/m²
};
// N equal-area bins, 1 <= N <= 64; factors finite, f_k >= 0, ascending,
// mean within 1e-9 of 1. The constructor throws std::invalid_argument naming
// the first rule broken; otherwise it divides the factors by their mean, so
// they average to 1 to rounding (N·eps), not to 1e-9.
class SnowDistribution {
  public:
    explicit SnowDistribution(std::span<const double> factors);
    std::span<const double> factors() const;
};
class SnowState {  // owns n·N solid (sp) and liquid (sw) values, a copy of
                   // the factors, and per triangle swe, sca, outflow of the
                   // last step
  public:
    // throws std::invalid_argument on a negative or non-finite initial swe
    static SnowState uniform(std::size_t triangles, const SnowDistribution&,
                             double initial_swe_mm);  // sp_k = f_k·swe, sw_k = 0
    std::span<const double> swe() const, sca() const, outflow() const;
};
struct StepForcing {
    std::span<const double> temperature;    // °C, n
    std::span<const double> precipitation;  // mm over the step, n, >= 0
    std::span<const double> radiation_anomaly;  // W/m², n; empty iff mode none
    double dt_seconds;                       // > 0
};
enum class StepStatus : std::uint8_t {
    ok, invalid_parameters, size_mismatch, non_finite_forcing,
    negative_precipitation, bad_step_length };
StepStatus step(const HbvSnowParameters&, const SunCorrectionParameters&,
                SnowState&, const StepForcing&, unsigned threads);
}
```

**Where each refusal lives.** Two homes, by when the input is fixed:

- **Construction throws.** `SnowDistribution`'s constructor (bad N, a
  negative, non-finite or descending factor, a mean off by more than 1e-9) and
  `SnowState::uniform` (a negative or non-finite initial SWE) throw
  `std::invalid_argument`; pybind11 turns it into `ValueError` with the same
  message. No object exists to be left half-made.
- **`step` returns a status.** `invalid_parameters`: `cx`, `lw` or `cfr`
  negative, or any of `tx`, `ts`, `cx`, `lw`, `cfr`, `r`, `d` non-finite, or a
  mode outside the enum. `size_mismatch`: temperature or precipitation not of
  length n, or the radiation anomaly not of length n under a correction or
  not empty under `none`. `non_finite_forcing`, `negative_precipitation`,
  `bad_step_length` (dt not finite and > 0) as named.
  Any status other than `ok` leaves the state untouched (validation runs over
  the whole input before any slot is written); the binding raises
  `ValueError` naming the status and, for a non-finite value, the first
  triangle.

`params.py` refuses the same values earlier, in Pydantic; the C++ checks stay,
because the binding can be called without `params.py`.

**One step, per triangle i and bin k**, with Δt = dt_seconds / 86400 days:

1. snow = P if T < tx, else 0; rain = P − snow.
2. sp_k += f_k·snow.
3. The melt factor m and threshold θ: `none`: m = cx, θ = ts.
   `melt_factor`: m = max(0, cx + r·ΔI_i), θ = ts. `threshold`: m = cx,
   θ = ts − d·ΔI_i.
4. If the bin has snow (sp_k > 0): sw_k += rain. Otherwise rain leaves at
   once (out_k += rain).
5. If T > θ: melt_k = min(m·Δt·(T − θ), sp_k); sp_k −= melt_k;
   sw_k += melt_k. Else if T < ts: refreeze_k = min(cfr·cx·Δt·(ts − T), sw_k);
   sp_k += refreeze_k; sw_k −= refreeze_k. (Between θ and ts, possible in the
   threshold form on a shaded slope, neither happens.)
6. out_k += max(0, sw_k − lw·sp_k); sw_k = min(sw_k, lw·sp_k).
7. swe_i = mean_k(sp_k + sw_k); sca_i = (number of k with sp_k > 0) / N;
   outflow_i = mean_k(out_k).

The sun corrects melt only. Refreezing and the rain/snow split keep `cx`, `ts`
and `tx`: the radiation term is a daytime energy gain, and night refreezing on
a sunny slope is not slower for it. Named as a choice (section 9, Q4).

**Invariants** (the suite tests each):
- I1, mass: swe_after = swe_before + P − outflow, per triangle per step, to
  1e-9 × max(1 mm, swe_before + P). Each bin's update conserves sp + sw + out
  exactly up to rounding; the snowfall spread over the bins sums to
  P·mean(f), and the constructor's renormalisation makes mean(f) = 1 to N·eps
  (7e-15 at N = 64), so the error is rounding of order 1e-14 relative to
  swe + P, not the 1e-9 the factors are first checked to.
- I2, bounds: sp_k ≥ 0, 0 ≤ sw_k ≤ lw·sp_k (+1e-12), 0 ≤ sca ≤ 1,
  outflow ≥ 0.
- I3, flat is plain: with ΔI_i = 0, all three modes give bit-identical
  states.
- I4, N = 1 is plain HBV-96 for one point.
- I5, the step is independent of `threads` (bit-identical states for 1 and 8).

### 4.3 The sun correction (`include/terrain/snow/radiation_index.hpp`)

```cpp
// out_i = (1/S) Σ_s beam_s · (max(0, n_i·u_s) − max(0, u_s.up))
StepStatus radiation_anomaly(std::span<const double> normals,  // 3n, unit, ENU, up >= 0
                             std::span<const double> sun,      // 3S, unit, ENU
                             std::span<const double> beam,     // S, W/m², >= 0
                             std::span<double> out, unsigned threads);
```

- u_s is the sun's direction at the midpoint of sub-interval s of the step
  (S equal sub-intervals; S = 4 for an hourly step, Hock's 15 minutes), and
  beam_s the clear-sky direct normal irradiance then, 0 with the sun below the
  horizon.
- The beam is forcing and is computed in Python (`snow/forcing.py`), as the
  radiation design's R9 rules (forcing is the consumer's, never the core's).
  The formula is this design's, Hock's form with P/P0 = 1 (R9 does not state
  one): beam = I0·(1/R)²·Ψa^(1/sin e) with
  I0 = 1368 W/m², Ψa = 0.75, R the Earth-Sun distance in AU and e the sun's
  elevation, both from TR1, and P/P0 = 1 (departure 3).
- The normals are per unit inclined area, as Hock's I is: ΔI is irradiance on
  the slope, and the melt it drives is depth on the slope. SWE in the output
  is per unit inclined area too (mm on the slope). Named, so that a reader
  comparing with a raster model knows the convention.
- **Normals in true east-north-up.** `snow/geometry.py` takes the triangle
  normal in the mesh's grid axes and turns it by the grid convergence γ at the
  mesh's centre (pyproj's `meridian_convergence`): one γ per mesh. On the
  synthetic domain (on UTM 33's central meridian) γ is 0.0000° at the centre,
  −0.0984° at the west edge and +0.0984° at the east edge
  (`46-probes/probe_sizes.out`). For a 100 km catchment
  away from a central meridian the spread is about 1°; the radiation design's
  per-receiver frame (TR2) replaces the single γ when it lands.

### 4.4 The sun's position: the radiation design's TR1, not a second copy

The radiation design's TR1 (SPA, Reda and Andreas 2004; the legacy port,
UTC as `std::chrono::sys_time<milliseconds>`, `geocentric_sun` /
`topocentric_sun`, about 540 lines) is exactly what this increment needs
from it, and nothing more: no horizons (TR3), no sky view (TR4), no
evaluator (TR5). Ola ruled in Q9 of that design that TR1 may be built early.
**PR 2 of this increment is TR1, as that design states it**, with one
addition to its binding: a vectorised call for one place and many times,

```python
def sun_positions(times_utc_ms: npt.NDArray[np.int64], lat: float, lon: float,
                  height_m: float) -> SunPositions  # elevation, azimuth (deg), distance (AU)
```

returning the topocentric elevation without refraction (refraction matters
within a degree of the horizon, where the beam is near zero) and the true
azimuth from north, clockwise. The radiation design file is on the unmerged
branch `design-terrain-radiation`; PR 2 brings it to `master` (section 9,
Q2). A shorter, less precise sun formula (NOAA's, about 40 lines, 0.01°)
would let PR 3 start without PR 2, at the price of two sun algorithms in the
tree; rejected for that reason, and offered to Ola as Q1.

### 4.5 The synthetic case

**Terrain** (`snow/terrain.py`, `rasputin synthetic-terrain`): a 12 km
square, 10 m grid (1201 × 1201 nodes), EPSG:25833, node (0, 0) at
(494000, 6806000), so the centre (500000, 6800000) is on UTM 33's central
meridian at 61.334° N, 15.000° E (probe: pyproj). Height, with x, y from the
centre and r, φ polar:

z = 300 + 1600·exp(−(r/2500)²)·(1 + 0.15·(1 − exp(−(r/1500)²))·cos 5φ)
  + 250·exp(−((x − 3500)² + (y + 3000)²)/900²) + noise

- a massif to 1900 m with five ridges (so every aspect occurs at every
  height) and a second hill to the south-east;
- noise: a random-phase surface with power spectrum k^−3, no wavelengths under
  40 m, scaled to 25 m RMS, from `numpy.random.default_rng(46)` (the exact
  construction is `46-probes/probe_spring.py`'s, which PR 3's `terrain.py`
  reproduces bit for bit); it gives the mesh the triangle density of real
  terrain (without it, 1 m tolerance gives 31 838 triangles; with it,
  646 870);
- measured (`46-probes/probe_spring.out`): z from 234.0 to 1924.4 m; slope
  from the grid's central differences mean 14.0°, max 51.2°.

**Mesh** (`46-probes/probe_sizes.out`; wall times include process start and
the write): `rasputin mesh --dem synth.tif --tolerance 1`, "quite fine":
646 870 triangles, 324 265 vertices, mean plan area 223 m² (an equilateral of
23 m sides), 0.73 s. At 2 m: 294 770 triangles, 147 936 vertices, 0.49 s; at
0.5 m: 1 279 127 triangles, 1.17 s.

**Forcing** (`SyntheticSpring`, all defaults chosen, not derived):
- 2027-04-01T00:00Z to 2027-07-01T00:00Z (91 days), hourly steps, S = 4;
- temperature at 300 m: −3 °C rising linearly to +14 °C, plus a daily cycle of
  ±4 °C (a cosine) peaking at 14:00 mean solar time at 15° E (13:00 UTC);
  lapse −0.6 °C per 100 m from the triangle's centroid height; evaluated at
  each step's midpoint;
- no precipitation (the kernel handles it; the synthetic spring is dry);
- initial snow: 600 mm SWE on every triangle, spread over N = 10 bins with
  CV = 0.4;
- HBV parameters: tx = 0, ts = 0 °C, cx = 3 mm/°C/day, lw = 0.1, cfr = 0.05;
- corrections: r = 0.015 mm/°C/day per W/m², d = 0.02 °C per W/m². Chosen so
  that (a) mid-April at 61° N, a 30° south slope's daily-mean ΔI of +68.9 W/m²
  (section 6) raises the melt factor by about a third (0.015 · 68.9 = 1.03), and (b) the two forms give the same
  extra melt at 4 °C above the threshold: r·ΔI·4 = cx·d·ΔI ⇔ d = 4r/cx.
- Frames: every 24 hours (at 24:00 UTC); 91 frames plus the start.

### 4.6 The animation format: EnSight Gold, geometry once

Three candidates, measured on the 1 m mesh (646 870 triangles) for 92 frames
of two per-triangle values (SWE and SCA):

| format | how a frame is stored | size, 92 frames | needs |
|---|---|---:|---|
| legacy `.vtk` series (numbered files, or a `.vtk.series` list) | the whole mesh again in every frame | ≈ 2.6 GB (23.3 MB measured mesh file + 2 × 2.6 MB arrays, × 92) | the existing writer plus float cell arrays |
| **EnSight Gold binary** (`.case` + `.geo` + one file per value per frame) | geometry once (11.7 MB), then 2 587 724 bytes per value per frame | **488 MB (465 MiB) measured** | a new writer, ≈ 60 lines, no dependency |
| VTKHDF with a static mesh | geometry once | about EnSight's | `h5py`, a new dependency |

**EnSight Gold.** The probe writer (`46-probes/probe_sizes.py`, its output
`probe_sizes.out`) wrote the 92 frames in 0.33 s: geometry 11 654 352 bytes,
2 587 724 bytes per value per frame, 487 796 213 bytes in all. ParaView
6.1.1's `pvpython` (`OpenDataFile`) read it back with 92 time steps (0, 24,
..., 2184 hours), 646 870 cells and the planted value at step 50. ParaView
opens `snow.case` and plays it; the time axis is in hours since the start.
**The case file lists the time values one per line**: VTK's EnSight reader
reads a case-file line of at most 256 characters, and with the 92 values on
one line the probe's first run failed with "'time values' not found". EnSight Gold coordinates are 32-bit floats, so the writer
stores x − x0 and y − y0 (x0, y0 the mesh's minimum corner, written in the
geometry file's description line and in `snow.json`): at 12 km a float's step
is 1 mm; at UTM northings it would have been 0.5 m.

Also written: `slope_deg`, `aspect_deg` and `radiation_anomaly_mean`
(ΔI averaged over the run) as variables without a time set, and
`melt_out_day` (the first frame whose SCA is below 0.5, in days; −1 if
never). A `snow.json` beside them records the job (the Pydantic model dumped)
and the origin. Run the three corrections into three directories and open
them side by side to compare.

`.gpkg` for QGIS is not an animation format at this size (60 million
features); a melt-out map as a GeoPackage layer is offered as an idea.

### 4.7 Is a continuous SCA worth it on small triangles?

Yes, as a distribution within the triangle, with a default CV of 0.4, for two
reasons:

- **Snow varies below any mesh's resolution.** Drifts and micro-relief vary
  snow depth over metres, so a 223 m² triangle is not one snow depth. The
  Nordic HBV's fixed-CV distribution says nothing about unit size, and here
  the same CV holds for every triangle. Whether the CV should grow with
  triangle area (variance grows with the support up to a correlation length of
  some hundred metres) is left out and handed back as an idea; it matters when
  triangles get large or are pooled into response units.
- **The melt looks right.** With CV = 0 (one bin) a triangle goes from full
  to bare in one step, so the snow line is a sharp edge along triangle edges.
  With CV = 0.4 it is a band. Probe (section 6, 2 m mesh): on 30° south
  slopes at 900-1100 m under the melt-factor correction, the median triangle
  takes 30 days from its first bare bin to bare, against 0 days (one step)
  with one bin. Over the domain under plain HBV, the snow-covered fraction at
  day 60 is 0.630 with ten bins and 0.888 with one: with one bin the cover
  only falls where a whole triangle has melted out, so the domain stays white
  until late and then goes bare in a rush (0.888 at day 60, 0.080 at day
  90).

**The price.** SCA moves in steps of 1/N (0.1 at N = 10), not continuously.
Shyft's interpolation between quantile points makes it continuous, at the
mass-balance cost named in 2.1. State memory is 2N + 3 doubles per triangle:
119 MB at N = 10 on the 1 m mesh. N is a parameter (1 to 64); N = 1 is
plain single-point HBV.

### 4.8 The Python surface

- `snow/params.py`: `HbvSnowParameters`, `SnowDistribution.lognormal(cv,
  bins)` (bin factors = the log-normal's quantiles at the bin midpoints
  (k + ½)/N, from `statistics.NormalDist`, divided by their mean; CV = 0
  gives one bin), `SunCorrection` (mode, r, d), `SyntheticSpring`,
  `FrameSpec`, `SnowJob`. Frozen; bounds as section 5.
- `snow/geometry.py`: `TriangleGeometry` (normals ENU, centroid z, area,
  slope, aspect) from vertices, triangles, γ.
- `snow/forcing.py`: `SyntheticSpring` → per step k: temperature (n),
  precipitation (n), sun (S × 3), beam (S). Times are timezone-aware UTC
  `datetime`s; a naive one is refused, as the radiation design rules.
- `snow/run.py`: `async def run(job, geometry) -> AsyncIterator[Frame]`;
  `Frame(time_hours, swe, sca)` with float32 copies.
- `io/vtk_read.py`: `read_mesh(data: bytes) -> MeshArrays`. Reads rasputin's
  own binary legacy VTK (header line `rasputin mesh`): the `crs` string from
  FieldData, POINTS, POLYGONS; skips LINES. Refuses ASCII ("write the mesh
  without --ascii"), another program's file, and a file with no `crs`.
- `io/geotiff_write.py`: `encode_geotiff(grid, x0, y0, step, epsg) ->
  bytes`, point-registered, float32, through tifffile, read back by
  `io/geotiff.py`.
- `io/ensight.py`: `geometry_bytes`, `variable_bytes`, `case_text`. Pure.
- CLI: `rasputin synthetic-terrain --out PATH` and `rasputin snow --mesh PATH
  --out DIR` with `--sun-correction`, `--start`, `--days`, `--step-minutes`,
  `--frame-hours`, `--initial-swe`, `--snow-cv`, `--bins`, `--threads`, and
  the HBV parameters as options with section 5's defaults. `--out` must be a
  new or empty directory.

### 4.9 Degeneracy and refusal policy

- A mesh whose CRS is geographic (degrees) is refused: normals need metres.
  ("snow needs a mesh in a projected CRS in metres; this one is in degrees").
- A degenerate triangle (zero area) gets the vertical normal and ΔI = 0; it
  is counted and reported on stderr, never a NaN.
- A triangle whose normal points down (none in a rasputin mesh, which is a
  height field) is turned up.
- Steep slopes: no limit. Snow does not slide off in this increment (the
  research note's 50-60° sloughing is part 2's).
- Non-finite forcing, negative precipitation, a step of zero or negative
  length, mismatched array lengths: refused, state untouched.
- Polar night and midnight sun: nothing special; the beam is 0 below the
  horizon and the sum handles 24 h of sun.

### 4.10 Not in this increment

Wind-distributed snowfall (part 2); cast shadows, sky view and reflected
radiation (the radiation design's TR3-TR7); diffuse radiation; canopy;
glaciers; per-triangle pressure in the beam; a CV that grows with triangle
area; real forcing (station or gridded weather); outflow routing; Shyft
adapters; a GeoPackage melt-out layer.

## 5. Constants

Each states the scale it assumes and the largest input it was checked at.

| constant | value | scale assumed | checked at |
|---|---|---|---|
| tx, ts | 0 °C | HBV-96 practice; chosen | synthetic spring, 61° N |
| cx | 3 mm/°C/day | Norwegian HBV range 2-5; chosen | probe, 294 770 triangles × 2184 steps |
| lw | 0.1 | Shyft's default | same |
| cfr | 0.05 | HBV-light's default (to verify); Shyft uses 0.5 | same |
| CV, N | 0.4, 10 | chosen (2.1, Liston's table not read) | same |
| r | 0.015 mm/°C/day per W/m² | 61° N, April; chosen (4.5) | probe |
| d | 0.02 °C per W/m² | = 4r/cx; chosen | probe |
| I0, Ψa | 1368 W/m², 0.75 | Hock 1999's form; Ψa chosen | 0-90° sun |
| S | 4 per hour | Hock's 15 minutes | probe |
| distribution mean tolerance | 1e-9 on input, then renormalised | doubles, N ≤ 64 | N = 64 |
| mass-balance tolerance | 1e-9 × max(1 mm, swe_before + P) | swe ≤ 10⁴ mm, P ≤ 100 mm per step, 10⁴ steps; factors renormalised | suite (P1); the probe's dry spring at 294 770 triangles, section 6 |
| float32 offset | x0, y0 = mesh minimum | meshes up to 10³ km (step 0.06 m) | 12 km |

## 6. What was measured

The probes are committed beside this file, in `docs/increments/46-probes/`,
each with its output: `probe_spring.py` (terrain, the 2 m mesh, the spring in
NumPy; `probe_spring.out`), `probe_sizes.py` (the other meshes, γ, the
EnSight set; `probe_sizes.out`), `probe_sun.py` (`probe_sun.out`) and
`noaa_sun.py`, the sun the other two share. Run 2026-10-13 on an M1 Max on AC
power, with the worktree's built `_core` for `rasputin mesh` and NumPy 2.5 for
the rest. The sun is NOAA's solar calculator (the spreadsheet's Meeus-based
formulas, about 0.01°), standing in for TR1; its highest elevations at
61.334° N, 15° E: 28.51° (2027-03-20), 38.44° (04-15) and 52.10° (06-21),
azimuth there within 0.15° of 180°. (A first run used NOAA's shorter
fractional-year formula, whose declination is about 0.3° low in April:
38.13° on 04-15. Its outputs were discarded.)

- **Mesh sizes and times:** section 4.5.
- **Daily-mean direct beam on 30° planes, 2027-04-15, 61.334° N 15° E,
  96 samples of 15 min, Ψa = 0.75, I0 = 1368:** horizontal 166.5 W/m²; south
  235.4; north 53.0; east 153.9; west 154.4. So ΔI = +68.9 (S), −113.5 (N),
  −12.6 (E), −12.1 (W).
- **Synthetic spring, 2 m mesh (294 770 triangles), hourly over 91 days
  (2184 steps), dry, section 4.5's defaults.** "Domain" figures are weighted
  by each triangle's inclined area; "class S" and "class N" are the triangles
  with centroid at 900-1100 m, slope 25-35° and aspect within 22.5° of south
  (782 triangles) or north (845); a triangle's melt-out day is the first
  daily frame whose SCA is below 0.5, −1 if none.

  | run | domain SCA, day 30 / 60 / 90 | domain SWE mm, day 30 / 60 / 90 | median melt-out day, class S / N | class S, days first bare bin to bare (median) |
  |---|---|---|---|---|
  | none, N = 10 | 1.0000 / 0.6298 / 0.0829 | 565.29 / 203.58 / 23.56 | 83.0 / 82.0 | not reached by day 91 |
  | melt_factor, N = 10 | 1.0000 / 0.6559 / 0.1152 | 565.14 / 226.07 / 36.05 | 75.5 / −1 | 30 |
  | threshold, N = 10 | 1.0000 / 0.6337 / 0.0932 | 558.13 / 206.93 / 28.22 | 75.0 / −1 | 35 |
  | none, N = 1 | 1.0000 / 0.8884 / 0.0804 | 566.15 / 149.13 / 18.10 | 83.0 / 83.0 | 0 |
  | melt_factor, N = 1 | 1.0000 / 0.7460 / 0.1102 | 566.36 / 186.88 / 29.31 | 76.0 / −1 | 0 |

  What this says. **The sun moves melt-out on 30° south slopes at 1000 m
  about a week earlier** (75 against 83 days), in both forms. **North
  slopes at that height do not melt out by 1 July** under either correction
  (their run-mean ΔI is −104.8 W/m², the south class's +43.7). **The domain
  keeps more snow with the correction than without** (226 against 204 mm at
  day 60 in the melt-factor form): the anomaly is relative to a horizontal
  plane, and over a domain of slopes facing every way the shaded loss
  outweighs the sunlit gain. Whether `cx` should then be recalibrated
  against the corrected run is a calibration question, outside this
  increment. **The threshold form melts earlier in the cold weeks** (558
  against 565 mm at day 30): it lets sunny slopes melt below 0 °C, where the
  melt-factor form, which keeps θ = ts, does not. **Mass:** the largest
  per-triangle residual of 600 − SWE − outflow over a run was 6.9e-12 mm.
- **Cost:** the five runs together, with ΔI computed once per step, took
  130.2 s in NumPy, one thread, 58.7 ms median per step. The C++ is expected
  to be several times faster per thread and to scale with threads; the 1 m
  mesh has 2.2 times the triangles. PR 3's acceptance run times the 1 m
  spring for real (section 8); no C++ figure is quoted before it.
- **Animation size:** section 4.6.

## 7. Red suites for `@tester`

Figures come from section 6's probes or are computed by hand here.

**PR 1, the HBV kernel** (Catch2 for the header; pytest for the binding and
`params.py`). Invariant-critical, so mutation-tested. Each target with the
fixture that kills it:

| mutant | killed by |
|---|---|
| rain/snow `T < tx` → `T <= tx` | H10 |
| the min dropped in melt | H6 |
| the min dropped in refreeze | H11 |
| the retention cap dropped or `lw·sp` → `lw·(sp + sw)` | H1, H4 |
| the bare-bin rain path (rain into an empty bin's sw) | H5b (in H5 the mutant is equivalent: the cap releases the rain anyway) |
| the SCA count `sp_k > 0` → `sp_k >= 0` | H5, H6, H7 |
| the mode switch: `none` applying ΔI | H8 (`none` melts 12, not 18) |
| the mode switch: melt_factor and threshold swapped | H12 (H8 cannot: the two forms agree there by design) |
| the clamp `m = max(0, …)` dropped | H8b |
| the refreeze condition `T < ts` → `T < θ` | H8b |
| the refreeze `else` dropped (refreeze also after melt) | H12 |

`T < ts` → `T <= ts` in the refreeze condition is an equivalent mutant
(refreeze is cfr·cx·Δt·0 = 0 there) and is not a target.

- H1 melt: one bin, sp = 100, T = 5, cx = 3, ts = 0, Δt = 1 day, P = 0 →
  sp 85, sw 8.5, outflow 6.5, swe 93.5, sca 1.
- H2 refreeze, from H1's state: T = −2, cfr = 0.05 → refreeze 0.3: sp 85.3,
  sw 8.2, outflow 0, swe 93.5.
- H3 snowfall: T = −1, P = 10, tx = 0, N = 4 factors (0.5, 0.8, 1.2, 1.5):
  sp rises by (5, 8, 12, 15), swe by 10.
- H4 rain on snow: sp = 100, T = 2, P = 10 → melt 6: sp 94, sw 9.4,
  outflow 6.6, swe 103.4.
- H5 rain on bare ground: sp = 0, T = 2, P = 10 → outflow 10, swe 0, sca 0.
- H5b rain on a bare bin below ts: one bin, sp = 0, sw = 0, tx = −1, ts = 0,
  T = −0.5, P = 10, cx = 3, cfr = 0.05, Δt = 1 day → outflow 10, sp 0,
  sca 0. Rain kept in the empty bin would refreeze 0.075 and give sca 1.
- H6 melt-out within a step: sp = 10, T = 10 → sp 0, outflow 10, sca 0.
- H7 bins: N = 4, factors (0.5, 0.8, 1.2, 1.5), initial 100 mm, T = 20,
  cx = 3, Δt = 1 day → sp (0, 20, 60, 90), sw (0, 2, 6, 9), sca 0.75,
  swe 46.75, outflow 53.25.
- H8 the corrections at their break-even: T = 4, ts = 0, cx = 3, ΔI = +100:
  melt_factor with r = 0.015 → melt 18 mm/day; threshold with d = 0.02 →
  θ = −2, melt 18 mm/day; `none` → melt 12 mm/day (one bin, sp = 100,
  Δt = 1 day each).
- H8b a shaded slope with water in the pack: one bin, sp = 100, sw = 5,
  T = 4, ts = 0, cx = 3, cfr = 0.05, Δt = 1 day, ΔI = −300. melt_factor
  (r = 0.015): m = max(0, 3 − 4.5) = 0, no melt; threshold (d = 0.02):
  θ = +6 > T, no melt, and no refreeze because T > ts. Both: sp 100, sw 5,
  outflow 0 (cap 10), swe 105. Without the clamp, m = −1.5 gives melt −6 and
  sp 106; refreezing on T < θ gives sp ≠ 100 (any sign of (θ − T) moves sp,
  since sw = 5 > 0).
- H10 rain at the threshold: one bin, sp = 100, sw = 0, T = tx = ts = 0,
  P = 10, cx = 3, Δt = 1 day: rain (T < tx is false), no melt, no refreeze:
  sp 100, sw 10, outflow 0 (cap 10), swe 110. `<=` would give sp 110, sw 0.
- H11 refreeze capped by the water: one bin, sp = 100, sw = 1, T = −10,
  cfr = 0.05, cx = 3, Δt = 1 day: the potential refreeze 1.5 > sw, so
  refreeze 1: sp 101, sw 0, outflow 0, swe 101. Without the min: sp 101.5,
  sw −0.5.
- H12 melt below ts in the threshold form: one bin, sp = 100, sw = 0,
  T = −1, ts = 0, cx = 3, cfr = 0.05, d = 0.02, ΔI = +100, Δt = 1 day:
  θ = −2 < T, melt 3: sp 97, sw 3, outflow 0 (cap 9.7). With the refreeze
  `else` dropped, a refreeze of 0.15 follows: sp 97.15, sw 2.85. The same
  input in melt_factor form (r = 0.015) melts nothing (T < ts) and refreezes
  min(0.15, 0) = 0: sp 100.
- H9 sub-day step: H1 with dt = 3600 s → melt 0.625 mm.
- P1-P5, properties over random forcing (Hypothesis in Python, Catch2
  generators in C++; 10⁴ steps, swe up to 10⁴ mm, T in [−30, 30],
  P in [0, 100], ΔI in [−600, 600]): I1-I5.
- V1 refusals, each in the home section 4.2 gives it:
  - `SnowDistribution` throws `std::invalid_argument` (Python: `ValueError`):
    a factor negative or NaN, factors descending, mean off by 1e-6, N = 0 or
    65. `SnowState::uniform` throws on initial SWE −1 or NaN.
  - `step` returns `invalid_parameters` for cx, lw or cfr = −0.1 and for NaN
    in any of the seven parameters; `size_mismatch` for temperature or
    precipitation of length n ± 1, ΔI of length n with mode `none`, ΔI empty
    with a correction; `non_finite_forcing` for NaN or ±inf in temperature,
    precipitation or ΔI; `negative_precipitation` for P = −1e-12;
    `bad_step_length` for dt = 0, −3600 and NaN. Each: the state's sp, sw,
    swe, sca and outflow compared before and after, bit for bit unchanged;
    from Python, `ValueError` naming the status.
- V2 renormalisation: factors (0.5, 0.8, 1.2, 1.5 + 4e-10) (mean off by
  1e-10, accepted) come back with mean 1 to 1e-15; P1's mass check then holds
  at the 1e-9 × max(1 mm, swe_before + P) tolerance with P up to 100 mm.
- D1 `SnowDistribution.lognormal(0.4, 10)` factors, to 1e-4: 0.4974, 0.6287,
  0.7228, 0.8080, 0.8930, 0.9838, 1.0873, 1.2154, 1.3973, 1.7664 (probe);
  mean 1 within 1e-12; CV 0 → one factor 1.0.

**PR 2, TR1:** the radiation design's own TR1 suite (its "Verification plan"),
plus `sun_positions` against the NOAA figures above (`46-probes/probe_sun.out`):
the day's highest elevation at 61.334° N 15° E within 0.02° of 38.44°
(2027-04-15) and 52.10° (2027-06-21), azimuth 180° within 0.2° at the
elevation maximum.

**PR 3, correction, synthetic case, animation.**
- R1 hand planes, one sample, sun at 30° elevation due south, beam 800:
  horizontal 400; ΔI = +292.8 (south 30°: the sun is 30° from the slope's
  normal, 800·0.866 − 400), −400 (north 30°: the sun grazes the slope),
  −53.6 (east 30°: 800·(0.5·0 + 0.866·0.5) − 400); sun below horizon
  (beam 0) → 0. Probe line `R1` in `46-probes/probe_spring.out`.
- R2 day means (section 6): through `sun_positions` and `snow/forcing.py`,
  2027-04-15 at 61.334° N 15° E: ΔI = +68.9, −113.5, −12.6, −12.1 W/m² for
  30° S, N, E, W within 0.5 W/m² (TR1 and NOAA's calculator differ by about
  0.01°; probe line `R2`).
- R3 γ: normals turned by +0.0984° at (506000, 6800000) EPSG:25833 (pyproj).
- R4 `read_mesh(write_vtk(...))` returns the same vertices, triangles and crs
  bit for bit; ASCII, foreign and crs-less files refused in plain words.
- R5 `encode_geotiff` → `io/geotiff.py` round trip: shape 1201 × 1201, step
  10, EPSG:25833, z range 234.0 to 1924.4 within 0.1 m for the default
  terrain.
- R6 EnSight: byte sizes exact (per value per frame: 244 + 4·T bytes; for
  T = 646 870, 2 587 724); with the `vtk` test extra, read back through
  `vtkGenericEnSightReader` (time values, a planted value at step 50, the
  cell count); skipped without it.
- R7 the run on a 2-triangle mesh (one 30° south facet, one 30° north)
  through 10 April days: the south facet's SWE falls faster under both
  corrections and equally under `none`; I1 over the run.
- E1 end to end (marked slow): `synthetic-terrain`, `mesh --tolerance 2`
  (294 770 triangles), `snow` for each of the three corrections with
  N = 10; against section 6's table, rows 1-3, with its definitions (inclined
  area weights, classes S and N, melt-out on daily frames): domain SCA at
  days 30, 60 and 90 within 0.01, domain SWE within 2 %, the median melt-out
  day of classes S and N within 1 day (class N: −1 under both corrections).
  The fixture's figures are measured, not derived; a miss is first checked
  against `46-probes/probe_spring.py` run on the same mesh.

## 8. Size, acceptance and order

| PR | what | C++ | bind | Python | net |
|---|---|---:|---:|---:|---:|
| 1 | `hbv_snow.hpp`, `bindings/snow.cpp` (`SnowState`, `step`, GIL released), `_core.pyi`, `snow/params.py` | 120 | 70 | 80 | ≈ 270 |
| 2 | TR1 as the radiation design states it, plus `sun_positions` | — | — | — | ≈ 560 (its estimate + 20) |
| 3 | `radiation_index.hpp` + binding; `snow/terrain.py`, `geometry.py`, `forcing.py`, `run.py`; `io/vtk_read.py`, `io/geotiff_write.py`, `io/ensight.py`; two CLI commands | 40 | 25 | 480 | ≈ 545 |

Each under 700 (`python3 tools/count_loc.py`). PR 1 and PR 2 touch disjoint
files and can be built at the same time by two writers; PR 3 needs both.
`project_structure.md` gains a `snow` module (header-only, depends on
`parallel_util` only; nothing depends on it) in PR 1.

- **`@perf`:** no refine or mesh code is touched, so the acceptance rule does
  not apply. PR 3's acceptance instead runs the default spring on the 1 m mesh
  (646 870 triangles) at 1 and 8 threads, records wall time, peak memory and
  output size in this file, and checks the animation in ParaView (one
  screenshot under `docs/benchmarks/<date>/` for Ola, the only image).
- **Mutation testing:** PR 1's kernel only (section 7).
- **Template cross products:** none.

## 9. Questions for Ola (defaults taken while he was away)

1. **The sun's position: the full SPA from the radiation design (PR 2, about
   560 lines, a PR on the path to the animation), or a short formula (about
   40 lines, 0.01°, accurate enough for melt) so the animation comes one PR
   sooner, at the price of two sun formulas in the tree until the SPA lands?**
   Default: the full SPA, built at the same time as PR 1.
2. **The radiation design is on an unmerged branch. May PR 2 bring that design
   file to `master`?** Default: yes, unchanged except TR1's status.
3. **Snow-covered area moves in steps of a tenth (ten bins). Good enough, or do
   you want Shyft's interpolation, which is smooth but can lose mass?**
   Default: steps of 1/N, N = 10.
4. **Should the sun also slow refreezing on shaded slopes and speed it on
   sunny ones?** Default: no, the sun changes melt only.
5. **Melt and SWE per square metre of slope (as Hock), or per square metre of
   map?** On a 30° slope they differ by 15 %. Default: per square metre of
   slope.
6. **The synthetic spring: dry, −3 to +14 °C at 300 m over April-June,
   600 mm of snow everywhere. Any snowfall or rain events in it?** Default:
   dry.
7. **Within-triangle snow variation (CV 0.4) the same on every triangle,
   whatever its size?** Default: yes; a size-dependent CV is an idea for
   later.

## Review

**Design review, round 1, 2026-10-13.** Range `b1353075..0e3247f7`, 0 production lines, @reviewer: CHANGES REQUESTED. Blocking: (1) placeholders left in the file: section 6 says `[RESULTS]` and section 4.7 `[single-bin figures in section 6]`; the melt-out and snow-cover figures the architect's handback gave (about a week earlier on 30° south slopes at 1000 m; cover 0.62 against 0.88 at day 60) are not in the file, and test E1 checks against them; (2) test R1's south case: with the sun at 30° and the slope at 30°, the angle from the slope's normal is 30°, so 800·0.866 − 400 = +292.8, not +400 (north −400 and east −53.6 are right); (3) test V1 asks the C++ to refuse negative cx, lw or cfr and bad bin factors, but `StepStatus` has no status for bad parameters, and the `SnowDistribution` constructor's refusal (exception or factory) is not said; (4) listed mutants not all killed: `<=` for `<` in the rain/snow test needs a fixture at T = tx with P > 0; dropping the min on refreeze needs a refreeze larger than sw; refreezing when T < θ instead of T < ts needs H8 to start with sw > 0; add the m ≥ 0 clamp and the refreeze condition as targets; (5) the mass-balance tolerance 1e-9·max(1, swe) is tighter than the error bin factors with a mean checked only to 1e-9 can give (P·1e-9): renormalise in the constructor, or scale the tolerance with P; (6) section 9 has seven questions, the ROADMAP row says eight. Wordings: `parallel_util::parallel_for` does not exist, the function is `for_each_block`; R9 of the radiation design does not state the beam formula; 465 MB is MiB (488 MB decimal, the unit of the 2.6 GB). Verified: the D1 factors, tests H1-H9 except R1, section 6's ΔI values (NOAA formula within 0.5 W/m²), latitude and γ, the EnSight byte sizes, TR1's estimate and split, the literature citations, the I/O boundary, std::chrono only; the full sun-position PR over the 40-line formula justified as the default. Idea: pass the step length as `std::chrono::duration<double>`.
