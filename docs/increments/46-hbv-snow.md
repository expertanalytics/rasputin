# Increment 46: HBV snow on the triangles, melt corrected for the sun, and a melt animation

**Status: designed; design review round 2 approved; addendum for cast
shadows (section 10) written 2026-10-13, not yet reviewed** (`@architect`,
branch `worktree-snow-melt`). PR 1's kernel is built. Six PRs (section 10.7,
which supersedes section 8's table). Ola was away when it was written; every default taken
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
  2. **No cast shadows** (withdrawn from PR 6 on, section 10.2). Self-shading only (a slope facing away from the sun
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
                          (dropped: the combiner is NumPy, section 10.5)
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

*Amended by section 10.5: the formula below is computed in NumPy, not in a C++
kernel; from PR 6 the factor max(0, n_i·u_s) comes from the radiation
evaluator, with cast shadows.*

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

Wind-distributed snowfall (part 2); sky view and reflected
radiation (the radiation design's TR4-TR7; cast shadows are in, section 10); diffuse radiation; canopy;
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
| the mode switch: `none` applying ΔI | V1 (a ΔI given under `none` is `size_mismatch`; H8's `none` row is given no ΔI) |
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
  *Retired in PR 6 with the single-γ turn it tests (section 10.5); TR2's
  per-triangle frame and PR 4's J tests replace it.*
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

*Superseded for the PR split and sizes by section 10.7; the acceptance and
mutation rules below still hold for PRs 1-3.*

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

## 10. Addendum: cast shadows in the melt season (2026-10-13)

Written by `@architect` on the snow branch at `85aec66f`, after PR 1's kernel
was built. Sections 1-9 stand except where this section names a change; each
change is marked in place with a pointer here. Every figure in this section
comes from `46-probes/probe_shadow.py`, whose output is
`46-probes/probe_shadow.out` (run 2026-10-13 on an M1 Max on AC power, NumPy
and the main checkout's `tin_engine` for `open_dem`), except the lines named
`W`, `F`, `M` and `S`, which come from `46-probes/probe_review3.py` and its
`.out` (design review round 3, same day and machine, no DEM read).

### 10.1 What Ola ruled

Ola, 2026-10-13, verbatim: "I need shades for the melting season." "I never
asked for a synthetic mountain. I asked for synthetic weather." "Let's do
Bygdin. On q2, I don't want to do a lot of extra work." "So - real alg, synth
data is ok." "Go for it."

Read with the main session's questions he answered:

- **Cast shadows are in this increment.** A triangle in the shadow of other
  terrain gets no direct sun, as the radiation design defines it: its horizon
  search (that design's TR3) and its lit test (its TR5), not a shortcut that
  a later radiation increment would throw away.
- **Synthetic data is fine where it saves work.** The weather is synthetic.
  The synthetic terrain (section 4.5) stays for the tests and the first
  animation.
- **Bygdin is the real terrain**: the Bygdin catchment on DTM10, meshed
  before (ROADMAP rows 22 and 16c), is the run that shows shadows on real
  relief. No data handling beyond what `rasputin mesh --dem` already does.

**A correction to this design's record.** The earlier version of this design
invented a synthetic mountain (section 4.5) where Ola had asked only for
synthetic weather. The mountain stays, because it is cheap, already probed,
and the tests need a terrain of known shape; Bygdin is the demonstration.

### 10.2 Prior art

- **Literature.** The horizon search, the receiver frame and the lit test
  are the radiation design's (branch `design-terrain-radiation`,
  `docs/increments/terrain-radiation.md`, its *Prior art* section: Dozier and
  Frew 1990, Stewart 1998, Corripio 2003, Marsh et al. 2012 and 2020); this
  addendum adds nothing to it. Hock 1999's I (section 2.1) already includes
  topographic shading, so **departure 2 of section 2.1 is withdrawn from PR 6
  on**: with shadows on, the correction is Hock's I with cast shadows, and
  departures 1 (the anomaly against an open horizontal plane) and 3 (one
  pressure) remain. No novelty is claimed.
- **Legacy.** Section 2.2's grep covers it: the legacy shaded a face by a CGAL
  ray from its centroid (`legacy/rasputin/triangulate_dem.h:600-672`). The
  idea is carried as the radiation design carries it (horizon tables, no
  CGAL); no legacy code is carried.

### 10.3 The ideas file

Read again for this addendum. Taken up: "Snowpack on the TIN" part 3 (melt
corrected for the sun), now with shadows; the terrain radiation design, which
it overlaps, is reused and not duplicated (10.4). Left out, as before:
"Snow follow-ups" (size-dependent variation, melt-out-day map, outflow per
frame).

### 10.4 Which radiation steps, and how much of each

The smallest set that gives a correct lit or shadowed answer per triangle at
each sun sample, each built as the radiation design states it:

| radiation step | what the snow needs from it | built here |
|---|---|---|
| TR1, the sun | the sun at each sample time | all of it: this increment's PR 2 (unchanged) |
| TR2, the receiver frame | per triangle: latitude and longitude, the normal in true east-north-up, the map from grid to true directions (its J), the exclusion radius, the observer height | all of it: PR 4 |
| TR3a, the horizon kernel | the horizon in N = 64 true directions per triangle, with Earth curvature and refraction; the exhaustive reference and the pruned search that must equal it bit for bit | all of it: PR 5 (built) |
| TR3b, the horizon's inputs | the **near raster**: the DEM within 20 km of the mesh, at full resolution | **the near raster only**, in PR 6 |
| TR4, sky factors and the geometry file | nothing: the lit test needs no sky view, and the horizons are computed in the run, not stored | none |
| TR5, the evaluator | the `instant` mode: per triangle and sample time, lit, cos(incidence), the beam factor b = lit·max(0, cos i) and the sine of the apparent sun elevation | the C++ evaluator, its binding and a synchronous Python call, in PR 6; not its CLI |
| TR6, TR7 | nothing (exact lit intervals and the `mean` mode; adapters) | none |

**What is left out of TR3b, and why that is safe in the melt season.** TR3b
also builds a **far raster** (the DEM to 200 km, as 100 m block means), a CLI
`rasputin radiation geometry`, and the run that sets N by latitude. The snow
needs none of them:

- The far field changes almost nothing between April and August. Probe
  `B3`, on 2000 random DTM10 nodes inside the Bygdin catchment, horizons
  searched to 20 km on the 10 m DEM and from 20 to 60 km on 100 m block
  means: terrain beyond 20 km stands above the near horizon in 4.36 % of the
  (node, direction) pairs, and nowhere higher than 2.87° above the
  horizontal. On 2027-04-15, 05-15 and 06-21 it changes the lit state of 11,
  18 and 10 (node, 15-minute) samples, all with the sun at or below 2.97°,
  where the direct beam before the cosine is at most 5.1 W/m² on 21 June
  (5.3 W/m² at 1 AU; 0.00 and 0.23 W/m² at the highest flips on 04-15 and
  05-15, 1.02° and 1.90°: probe `W`, `46-probes/probe_review3.out`); the beam
  energy it removes rounds to 0.0000 % of the lit beam on each day. It
  matters in winter, at low sun, which is the radiation work's, not the melt
  season's.
- **Beyond 60 km**, where probe `B3` did not search, a bound replaces a
  measurement: terrain at ground distance d gives at most
  tan h = (z_max − z_obs)/d − c·d, which falls with d, so d = 60 km bounds all of
  it. With z_max = 2469 m (Galdhøpiggen, Norway's highest summit; B0's
  near raster tops out at 2468 m) and the lowest receiver of probe `B1`
  (1054 m), the horizon from beyond 60 km is at most 1.116°; at 1500 m,
  0.691°; at 2240 m, below the horizontal (probe `F`). The direct beam at
  1.116° is 0.001 W/m² (probe `W`). The bound assumes DTM10's coverage, as
  `reach` does (the radiation design's Q4); terrain outside it is no
  obstacle.
- So `RadiationRequest` keeps the radiation design's fields, and PR 6 runs
  with the search reach `D_far` equal to `D_near` (20 km), with no far raster.
  TR3a's kernel takes the far raster as optional (`None` from Python) for
  exactly this; `reach` records where each direction's search stopped, as its
  R4 rules.
- N = 64, the radiation design's default until its acceptance run (its Q1).

**The near raster comes from the existing mesh path.** `radiation/search.py`
calls `open_dem(DemRequest(sources=..., bounds=<mesh box grown by D_near>))`,
the call `rasputin mesh --dem` makes, and `raster.to_core` for the view.
Probe `B0`: for Bygdin that is 7314 × 6195 nodes (181 MB as float32),
assembled in 0.9 s from `DTM10_UTM33_20260925`, no NaN, heights 15 to
2468 m. No new reader, no new data. A synthetic DEM works the same way; its
rays leave the raster at its edge, and `reach` says so.

**What TR5 is built without.** The radiation design has TR5 read a stored
geometry (TR4's file). Here its C++ evaluator takes the receiver set (TR2)
and the horizon table (TR3a) as arrays, and PR 6 adds
`radiation/geometry.py`'s `RadiationGeometry` with the fields TR2 and TR3
produce (latitude, longitude, normal, horizon, reach, N, the request). TR4
later adds the sky-view and terrain factors and the `.npz` save and load to
the same class, and TR5's CLI lands with it. The evaluator is not written
twice.

### 10.5 How the lit factor enters the sun correction

The HBV kernel (PR 1) does not change: it takes ΔI per triangle and per step,
as before.

**ΔI from the evaluator's beam factor.** With S sample times in a step (S = 4
for an hourly step), beam_s the clear-sky direct normal irradiance at sample
s (section 4.3's forcing, unchanged), and b_is and sin e_is the evaluator's
beam factor and the sine of the apparent sun elevation for triangle i at
sample s:

  ΔI_i = (1/S) Σ_s beam_s · (b_is − max(0, sin e_is))

- The reference stays an **open horizontal plane** (section 2.1, departure
  1), so a flat triangle in the open melts as plain HBV, and a flat triangle
  on a shaded valley floor gets ΔI = −beam·sin e while it is shaded: it melts
  less than plain HBV, which is the effect asked for.
- **This is the radiation design's combiner** (its R9, the beam term) minus
  the horizontal reference; forcing stays in Python, as R9 rules.
- The combiner is NumPy, elementwise over (S, n). Section 4.3's C++ kernel
  `radiation_index.hpp` and its binding are therefore **dropped from PR 3**
  (amendment): PR 3 computes the same formula in NumPy with the
  self-shading-only factor b_is = max(0, n_i·u_s) and sin e_s from the
  domain's one sun (section 6's 58.7 ms per step in NumPy included this).
  PR 6 replaces that source of b by the evaluator, and drops PR 3's single-γ
  normal turn (section 4.3) for TR2's per-triangle frame. After PR 6, b always
  comes from the evaluator: with `--shadows off` the evaluator gets a horizon
  of −90° in every direction (no terrain), which leaves only self-shading.
  One source of b, one source of normals. E1 (section 7) then runs with
  `--shadows off`; if it misses its tolerances, the cause to check first is
  the evaluator's refraction and per-triangle sun, and the fixture is
  re-measured in an amendment commit (`docs/increments/README.md`).
- `--shadows on|off`: default `on` when `--dem` is given, `off` without it;
  `--shadows on` without `--dem` is refused ("shadows need the DEM the mesh
  was made from: give --dem"). `SnowJob` gains `shadows: RadiationRequest |
  None`, the radiation design's frozen request, not a copy of its fields.
- **The sign of γ, pinned for PR 3.** Section 4.3's single-γ turn, while it
  lives: true azimuth = grid azimuth + γ, since pyproj's
  `meridian_convergence` is the true azimuth of grid north (the radiation
  design's P4); a grid-north vector (0, 1) becomes (sin γ, cos γ). Test R3
  means this. (`probe_spring.py` turns the other way; at its γ = 0 that
  changes nothing.) **Test R3 is retired in PR 6**, which deletes the turn;
  PR 4's J tests (10.8) pin the same sign through TR2.

**Async and threads.** The evaluator runs in the same `asyncio.to_thread`
batch as the HBV step; its binding releases the GIL (TR5 as designed) and
takes the run's `--threads`. The horizon search runs once per run, before the
first step, in its own `to_thread` call.

**Cost, estimated, for `@perf` to measure in PR 6's acceptance.** From the
radiation design's unmeasured estimates (its R4 and R5): about 3 ns per
sample and 4000 samples per 20 km ray at 5 m, so about 0.8 ms per triangle at
N = 64 before pruning; Bygdin's 10 m mesh (52 198 triangles) is then about
40 s on one thread before the pyramid prunes it. The evaluator: about 29 ns
per triangle and sample; the synthetic 2 m mesh (294 770 triangles) over 91
days at S = 4 is 2.6·10⁹ evaluations, about 75 s on one thread. For scale
only: the probe's NumPy search of 2000 Bygdin nodes to 20 km took 26.6 s.

### 10.6 What the shadows change, measured

Probes `S1` and `B2`: 3000 random nodes of the synthetic DEM and 2000 random
DTM10 nodes inside the Bygdin catchment, each node a receiver with the DEM's
normal (central differences, turned to true north by Bygdin's γ = −5.724°),
exclusion radius 10 m, N = 64, NOAA's sun, the beam of section 4.3, 96
samples of 15 minutes a day, no refraction. "Removed" is the share of the
day's direct beam on the slopes (self-shading already applied) that cast
shadows take away.

| site | 2027-04-15 | 2027-05-15 | 2027-06-21 |
|---|---|---|---|
| synthetic massif: beam removed / mean ΔI change / 5th percentile | 0.7 % / −1.1 / −5.6 W/m² | 0.5 % / −1.2 / −5.8 W/m² | 0.5 % / −1.3 / −6.1 W/m² |
| Bygdin: beam removed / mean ΔI change / 5th percentile | 1.6 % / −2.7 / −14.3 W/m² | 1.3 % / −3.0 / −14.9 W/m² | 1.2 % / −3.3 / −16.5 W/m² |
| Bygdin: nodes losing more than 20 W/m² of daily-mean ΔI | 3.1 % | 2.9 % | 3.5 % |

What this says. **In the melt season the shadows are a local effect, not a
domain-wide one**: 1-2 % of the beam over the catchment, but about 3 % of
Bygdin's points lose more than 20 W/m² of daily-mean beam, about 30 % of the
gain a 30° south slope has in mid-April (section 4.5's 68.9 W/m²). The
synthetic massif, a single dome, shades itself little; Bygdin's valleys and
cirques shade more (median of the mean horizon 9.45° against the massif's
6.87°). The figures are on DEM nodes, not triangles; PR 6's acceptance
measures them on the mesh. The probe places node (0, 0) at the mosaic's
(x_min, y_max), which `io/geotiff.py`'s `_placement` already moves to the
first cell's centre for an area-registered tile such as DTM10's.

### 10.7 The PRs

| PR | what | depends on | C++ | bind | Python | net |
|---|---|---|---:|---:|---:|---:|
| 1 | the HBV kernel (built on this branch) | — | 120 | 70 | 80 | ≈ 270 |
| 2 | TR1 and `sun_positions` (unchanged) | — | — | — | — | ≈ 560 |
| 3 | the correction in NumPy (10.5), the synthetic case, the animation; no `radiation_index.hpp` | 1, 2 | 0 | 0 | 495 | ≈ 495 |
| 4 | TR2, as the radiation design states it | 2 | 80 | 40 | 230 | ≈ 350 |
| 5 | TR3a, as the radiation design states it, the far raster optional | 4 | 300 | 55 | 0 | ≈ 355 |
| 6 | TR5's evaluator (`instant`, no CLI), TR3b's near raster (`radiation/search.py`), `RadiationGeometry`'s fields, the snow wiring (`--dem`, `--shadows`, b from the evaluator, normals from TR2) | 3, 5 | 220 | 50 | 220 | ≈ 490 |

PR 3's 495 is section 8's 545 less the radiation kernel and its binding (65)
plus the NumPy combiner (15). If PR 3's count passes 650 at its green step, the EnSight writer
and its test R6 split off as PR 3b. PR 6's Python: the evaluator's synchronous call
50, `search.py` 60, `RadiationGeometry` 30, the snow wiring 80 (job and CLI
fields, the run's evaluator call, the geometry from TR2). PR 5's binding is
55 lines against the radiation design's 50 for TR3a: the 5 are the far
raster taken as optional (an `std::optional` view, `None` from Python, and
the check that sets D_far to D_near when it is absent, 10.8's pin 4). Each
under 700 (`python3 tools/count_loc.py`).

**The record in the radiation design's file.** PR 2 brings
`docs/increments/terrain-radiation.md` to `master` (question 2's default) and,
in the same commit as TR1's status, adds under its *Increments and LOC* table
a paragraph "Built by increment 46", which says exactly this:

- TR2: all of it, in increment 46's PR 4.
- TR3a: all of it, in PR 5, with the far raster optional (absent means
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

PRs 4, 5 and 6 each mark their line done in that file when they land, so the
file never calls a built step unbuilt and no second version is written. If
Ola says no to question 2, PR 2 makes the same edit as a commit on
`design-terrain-radiation` instead.

**Order.** PR 1 and PR 2 at the same time (as section 8). Then PR 3 and PR 4
at the same time: PR 3 no longer touches C++ or `_core.pyi`, so they share
only `project_structure.md`. PR 5 after PR 4, PR 6 last. The first animation
(PR 3) has self-shading only; PR 6 adds the shadows and the Bygdin run.

**Not built here, and still the radiation design's**: TR3b's far raster, its
CLI and the N-by-latitude run; TR4; TR5's CLI; TR6; TR7. They follow in that
design's order.

**Mutation testing.** PR 1's kernel (section 7) and **PR 5**, which the
radiation design names invariant-critical (pruned equals exhaustive, its I4,
plus the closed forms); its targets and the fixture that kills each are in
10.8. PR 6: none.

**`@perf`.** No refine or mesh code is touched. PR 3's acceptance as section
8 says. PR 6's acceptance is the Bygdin run (10.9) and the synthetic 2 m run
with shadows: horizon-search time, evaluator time per step, peak memory, at 1
and 8 threads, recorded in this file, and one screenshot of Bygdin for Ola.

### 10.8 Red suites for `@tester`

**PR 4, TR2.** The radiation design's verification 2 (J against a
finite-difference transform at Solund, Oslo, the central meridian and Vardø;
a synthetic plane off the central meridian gives its true aspect; Solund
−8.906° and Vardø +15.218° through J), plus Bygdin: γ = −5.724° at
61.3694° N, 8.4852° E in EPSG:25833, through J, within 0.001° (probe `B0`;
−5.994° at the catchment's west edge, −5.454° at its east edge).

**PR 5, TR3a.** The radiation design's verification 3 and its I4 suite,
except the flat plane's sunrise, which needs the evaluator and moves to PR 6
(test S); the flat plane's horizon, −atan(c·r_i), stays here (−0.0000391° at
r_i = 10 m, probe `S`). I4 is run at 1 and 8 threads. The
V-valley, measured on the reference definition (probe `V1`): side slopes 30°,
ridges at half-width 1000 m, 10 m grid with the axis on a grid line, receiver
on the axis floor, exclusion radius 5 m, N = 64. With c = 0, tan h equals
|sin A|·tan 30° to 4.7·10⁻¹⁴ in every direction (h(45°) = 22.2077°,
h(90°) = 30.0000°, h(0°) = h(180°) = 0); with c = 6.828·10⁻⁸ m⁻¹
(k = 0.13, R = 6371 km), tan h equals |sin A|·tan 30° − 5c to 5.1·10⁻¹⁴
(h(0°) = −0.000020°). Also: a far raster of `None` gives the same table, bit
for bit, as a far raster that is never reached (D_far = D_near).

*Pinned for PR 5, where the radiation design leaves a choice open* (each
pin is what the mutation fixtures below assume):

1. `reach` is found from coverage alone before the height march: the first
   sample off the raster or reading a NaN node; its distance floored to
   whole km; 255 when no sample is missing. The march ends there. So the
   pruned and the exhaustive search agree on `reach` bit for bit too, and I4
   covers both arrays.
2. The ray's early stop bounds the rest of the ray with the largest maximum
   of the rasters still ahead on it: both while on the near raster, the far
   raster's alone after D_near.
3. Far samples are at D_near + j·Δ_far, j ≥ 1.
4. An absent far raster sets D_far to D_near.
5. The int16 is 100·h rounded half away from zero.

*Mutation targets, each with the fixture that kills it.* Checked on a
one-dimensional model of the pruned search as R4 describes it
(`46-probes/probe_review3.py`, lines `M` of its `.out`): a ray due east
along a grid row, a max-pyramid over 2^L-node blocks, the bounds and pins
above; the exhaustive loop is the reference, and both agree on every
fixture. In the C++ test each fixture is a raster whose heights depend on x
only, the receiver on a node, J the identity, and the check is direction
k = 16 (east) plus I4 in all 64. "h" is the int16 table entry.

| fixture | terrain and receiver | reference h, reach |
|---|---|---|
| F1, above all terrain | near 10 m, flat at 0 with a 2000 m ridge node at 10 km, to 21 km; far 100 m, flat at 1000 m, to 200 km; z_obs 3000 m, r_i 10 m, D_near 20 km, D_far 200 km, c = 6.828·10⁻⁸ | −134 (−1.3388°), 255 |
| F2, far peak at 25 km | near 10 m, flat at 0; far 100 m, flat at 0 with a 100 m node at 25 km; z_obs 0, r_i 10 m, D_near 20 km, D_far 200 km, same c | 13 (0.1314°), 255 |
| F4, spike at r_i | near 10 m, flat at 0 with a 100 m node at 10 m; no far raster; z_obs 20 m, r_i 7.5 m, D_near 2 km, c = 0 | 8223 (82.2348°), 255 |
| F5, near and far disagree | near 10 m, flat with a 1000 m node at 3 km (past D_near); far 100 m, flat with a 1000 m node at 1.5 km (inside D_near) and a 400 m node at 4 km; z_obs 0, r_i 10 m, D_near 2 km, D_far 5 km, c = 0 | 571 (5.7106°), 255 |
| F6, coverage gap | near 10 m, flat at 0 with a NaN node at 3.7 km and a 1000 m node at 5 km; no far raster; z_obs 0, r_i 10 m, D_near 20 km, c = 0 | 0, 3 |
| F7, far raster `None` | near 10 m, flat at 0; far `None`, D_far asked 60 km; z_obs 100 m, r_i 10 m, D_near 20 km, c = 6.828·10⁻⁸ | −36 (−0.3647°), 255 |

| target (the mutant) | killed by | mutant's h, reach |
|---|---|---|
| the pruning bound's a < 0 branch: a/d_lo for a/d_hi | F1 | −575, 255 |
| the a ≥ 0 branch: a/d_hi for a/d_lo | F2 | 0, 255 |
| the curvature term at d_lo: c·d_hi for c·d_lo | F1, F2 | −136; 0 |
| the early stop with the current raster's maximum only (pin 2) | F2 | 0, 255 |
| the search starting at r_i + Δ_near instead of r_i | F4 | 7720, 255 |
| the search starting at 0 (d_j = j·Δ_near) | F4 | 8287, 255 |
| the switch: the near raster read past D_near | F5 (also F1, F2) | 1843, 4 |
| the switch: the far raster read from r_i | F5 (also F2, F4, F6, F7) | 3080, 255 |
| recording `reach`: the march not ended at the first missing sample | F6 | 1131, 3 |
| recording `reach`: rounded, not floored (pin 1) | F6 | 0, 4 |
| recording `reach`: not recorded (255) | F6 | 0, 255 |
| the empty far raster: `None` read as zeros to the asked D_far (pin 4) | F7 | −30, 255 |
| the int16: truncated, or floored (pin 5) | V-valley h(45°); F7 | 2220 for 2221; −37 for −36 |

**PR 6, the evaluator, the near raster and the snow wiring.**

- L1, the lit test (the evaluator's per-sample function, on a hand-made
  horizon row): the V-valley's row stored as the radiation design stores it,
  int16 in 0.01°, round(100·atan(|sin A_k|·tan 30°)), so k = 16 (90°) is
  3000, k = 21 (118.125°) 2698 and k = 22 (123.75°) 2564; floor normal up.
  - Sun at azimuth 90°: apparent elevation 29.99° → unlit; 30.01° → lit;
    exactly 30.00° → unlit (lit is strict, the radiation design's degeneracy
    policy).
  - Sun at azimuth 120°: the horizon interpolated between k = 21 and 22 is
    26.5333° (the closed form is 26.5651°; probe `V2` gives 26.5372° before
    quantisation): 26.52° → unlit, 26.55° → lit.
  - Self-shading: a 30° slope facing north, horizon −90° all round, sun at
    azimuth 180° and 20° → unlit, b = 0.
- L2, the combiner (probe `V3`): the V-valley floor, sun at azimuth 90°,
  beam 800 W/m², one sample: at 20° elevation b = 0, reference 0.3420,
  ΔI = −273.6 W/m²; at 40°, b = reference = 0.6428 and ΔI = 0 within
  10⁻⁹ W/m². With a horizon of −90° all round (`--shadows off`), a flat
  triangle's ΔI is 0 within 10⁻⁹ W/m² at every sample.
- L3, shadows only take away: on the synthetic 2 m mesh, 2027-04-15, for
  every triangle and sample, b with the searched horizon ≤ b with −90°, and
  no triangle's daily-mean ΔI is higher with shadows than without.
- L4, the near raster: `search.py` on `synth.tif` asks `open_dem` for the
  mesh box grown by 20 km and gets the whole 1201 × 1201 grid (coverage ends
  there); a receiver at (497200, 6800000) has the raster's west edge 3.2 km
  away, and its `reach` toward 270° is 3 km.
- L5, refusals: `--shadows on` without `--dem`; a DEM in another CRS than the
  mesh; a DEM that does not cover the mesh. Each in plain words, before any
  step runs.
- I7, the evaluator against the reference sun (the radiation design's I7):
  receivers at Bygdin's centre (61.3694° N, 8.4852° E, 1500 m), Solund and
  Vardø, every 15 minutes of 2027-04-15, 05-15 and 06-21 with the sun above
  −1°: the evaluator's apparent elevation and true azimuth within 0.003° of
  TR1's `solar_position` (the parallax it drops is at most 0.0024°, that
  design's R5). On L1's V-valley horizon row, the evaluator's lit agrees
  with lit formed from the reference sun and the same row, except where the
  apparent sun is within 0.003° of the interpolated horizon.
- I8, thread count (the radiation design's I8): the synthetic 2 m mesh
  (294 770 triangles) on 2027-04-15, 96 samples of 15 minutes, horizons
  searched and evaluated at 1 and at 8 threads: horizon, `reach`, lit,
  cos i, b and sin e equal bit for bit, and so the day's ΔI.
- S, the flat plane's sunrise (moved from PR 5; the radiation design's
  verification 3): a flat DEM, a vertical receiver facing east at Bygdin's
  centre (so self-shading does not decide), r_i = 10 m; its horizon is
  −atan(c·r_i) = −0.0000391°, stored as 0 in every direction. With `limb`
  the upper limb, the test bisects TR1's reference path for the time the
  apparent upper limb meets the stored horizon, and the evaluator is unlit
  4 s before it and lit 4 s after it. NOAA's formula, standing in for TR1
  (probe `S`): 04:06:00 UTC at azimuth 67.84° on 2027-04-15, the sun rising
  0.111°/min, so 0.003° is 1.6 s; 01:44:30 UTC at 31.04° on 2027-06-21,
  0.062°/min, 2.9 s. Both days are tested.
- E2 (marked slow): `synthetic-terrain`, `mesh --tolerance 2`, `snow --dem
  synth.tif --shadows on`: the domain-mean change in daily-mean ΔI on
  2027-04-15 between −0.5 and −2.5 W/m² (probe `S1`, on 3000 DEM nodes:
  −1.1); the figure on the triangles is recorded at the PR's acceptance, and
  the band then narrowed to it in an amendment commit.

### 10.9 The Bygdin run (PR 6's acceptance, not CI)

```sh
rasputin mesh --dem $DTM10 --domain docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson \
    --tolerance 10 --out bygdin.vtk                      # 52 198 triangles (ROADMAP row 22)
rasputin snow --mesh bygdin.vtk --dem $DTM10 --shadows on --days 122 --out bygdin_snow/
```

`$DTM10` is `../rasputin_data/DTM10_UTM33_20260925`. The weather is section
4.5's synthetic spring run to 1 August (122 days), because Bygdin lies from
1054 to 2240 m (probe `B1`'s nodes) and at the spring's lapse rate much of it
would still be white on 1 July. `@perf` records the times and memory (10.7),
the beam removed by shadows on the three dates against 10.6's node figures,
and one screenshot. A second run with `--shadows off` beside it shows the
difference in ParaView.

### 10.10 The radiation design's open questions

Its Q1-Q10 were all ruled by Ola on 2026-09-28 ("go with c, and yes to the
rest"). None blocks this work. Those that touch it, with the ruling it takes:
Q1, N = 64; Q3, refraction coefficient k = 0.13, constant; Q4, terrain
outside the DEM's coverage is no obstacle and `reach` records it (Bygdin's
near raster has no NaN, probe `B0`); Q10, the sun's centre and a lit test
that is true or false. Q2 (pvlib in TR1's tests) is PR 2's, as before. Q5-Q9
do not touch it. Its *Facts to verify* hold nothing this work depends on.

### 10.11 Questions for Ola (defaults taken)

8. **Mountains more than 20 km away are left out of the shadow search for
   now.** On Bygdin they only matter with the sun below 3°, where the direct
   sun is too weak to melt anything; they matter in winter. Default: leave
   them for the radiation work.
9. **The first animation (on the synthetic mountain) comes without cast
   shadows, and the shadows come three PRs later with the Bygdin run.** The
   alternative holds the first animation back until the shadows are built.
   Default: the first animation first.
10. **The Bygdin run uses the same synthetic weather, run to 1 August instead
    of 1 July,** because Bygdin lies above 1050 m. Default: yes.

## Review

**Design review, round 1, 2026-10-13.** Range `b1353075..0e3247f7`, 0 production lines, @reviewer: CHANGES REQUESTED. Blocking: (1) placeholders left in the file: section 6 says `[RESULTS]` and section 4.7 `[single-bin figures in section 6]`; the melt-out and snow-cover figures the architect's handback gave (about a week earlier on 30° south slopes at 1000 m; cover 0.62 against 0.88 at day 60) are not in the file, and test E1 checks against them; (2) test R1's south case: with the sun at 30° and the slope at 30°, the angle from the slope's normal is 30°, so 800·0.866 − 400 = +292.8, not +400 (north −400 and east −53.6 are right); (3) test V1 asks the C++ to refuse negative cx, lw or cfr and bad bin factors, but `StepStatus` has no status for bad parameters, and the `SnowDistribution` constructor's refusal (exception or factory) is not said; (4) listed mutants not all killed: `<=` for `<` in the rain/snow test needs a fixture at T = tx with P > 0; dropping the min on refreeze needs a refreeze larger than sw; refreezing when T < θ instead of T < ts needs H8 to start with sw > 0; add the m ≥ 0 clamp and the refreeze condition as targets; (5) the mass-balance tolerance 1e-9·max(1, swe) is tighter than the error bin factors with a mean checked only to 1e-9 can give (P·1e-9): renormalise in the constructor, or scale the tolerance with P; (6) section 9 has seven questions, the ROADMAP row says eight. Wordings: `parallel_util::parallel_for` does not exist, the function is `for_each_block`; R9 of the radiation design does not state the beam formula; 465 MB is MiB (488 MB decimal, the unit of the 2.6 GB). Verified: the D1 factors, tests H1-H9 except R1, section 6's ΔI values (NOAA formula within 0.5 W/m²), latitude and γ, the EnSight byte sizes, TR1's estimate and split, the literature citations, the I/O boundary, std::chrono only; the full sun-position PR over the 40-line formula justified as the default. Idea: pass the step length as `std::chrono::duration<double>`.

**Design review, round 2, 2026-10-13.** Range `df70bc3b..27fab48d`, 0 production lines, @reviewer: APPROVED. All six blockers and three wordings fixed; the four probes re-run, every figure reproduced except wall times (the five-run table, melt-out medians 83, 75.5 and 75.0, north −1, 226 against 204 mm, the 6.935e-12 residual, section 4.7's 0.630/0.888, R1 +292.8/−400/−53.6, the EnSight bytes, pvpython's 92 steps); five new fixtures (H5b, H8b, H10, H11, H12) hand-checked as killing; the one undetectable mutant (`T < ts` to `T <= ts`) equivalent, since the refreeze is 0 there; renormalised bin factors make 1e-9·max(1, swe+P) sound. New assumptions sound: NOAA's longer sun calculator (the short one is 0.31° low on 04-15), the temperature peak at 14:00 solar time, case-file times one per line. Suggestion: H11 should state ts = 0.

**Design review, round 3 (the addendum, section 10), 2026-10-13.** Range `85aec66f..91d65f3c`, 0 production lines, @reviewer: CHANGES REQUESTED. Blocking: (1) PR 5, the horizon kernel, is invariant-critical but names no mutation targets or killing fixtures; name at least the pruning bound's a < 0 branch (a receiver above all terrain), the curvature term at d_lo, the early stop, the search starting at r_i, the switch from near to far raster, recording `reach`, the empty far raster; (2) PR 6 builds the radiation design's evaluator but leaves out that design's checks for it, I7 (within 0.003° of the reference sun, lit agreeing) and I8 (same result at any thread count); verification 3's flat-plane sunrise needs the evaluator, so it moves from PR 5 to PR 6; (3) say which PR records in `terrain-radiation.md` what is built here and how (TR2; TR3a with the far raster optional; TR3b's near raster only; TR5 without its CLI), since question 2's default would leave that file calling them unbuilt and risk a second version. Wordings: the beam at 2.97° is 5.1-5.3 W/m², not 5.6; the far-field probe covers 20-60 km only, so add the bound beyond 60 km (at most 1.1° from the lowest receiver); say test R3 (the single-γ turn) is retired in PR 6; PR 5's binding is 55 lines against the radiation design's 50 with no reason given. Held: `probe_shadow.py` re-run on DTM10, every figure reproduced (timings within 3 %); γ's direction as pinned; the lit-test figures, V-valley closed forms and combiner values; PR sizes and PR 3's 495; Bygdin's 52,198 triangles; questions 8-10.

**Design review, round 4 (the addendum), 2026-10-13.** Range `52a3bddf..9793bc83`, 0 production lines, @reviewer: APPROVED WITH WORDINGS (pin 2 says the maxima of the rasters still ahead; a PR 3b split point at 650), applied by the main session. All three blockers and four wordings of round 3 fixed; `probe_review3.py` re-run, output identical; section 10.8's 13 targets on six fixtures, pruned equal to exhaustive on every fixture, every mutant killed in the one-dimensional model; I7, I8 and test S in PR 6, the vertical east-facing receiver right for test S; pins 1, 3, 4, 5 hold, and coastal DTM10 tiles (`67m1_1`, `67m1_2`, `6700_1`, `66m1_1`) hold no NoData, the sea stored as 0, so pin 1's stop at the first missing sample bites only at the mosaic's edge. PR 3 at 495 against a 37 % overrun is about 678, under 700, so the split point is a margin, not a requirement.

**Code review, round 1 (PR 1, the HBV kernel), 2026-10-13.** Range `origin/master..ced40c13`, 369 net production lines (about 270 estimated: C++ 168 against 120, mostly the full validation; binding 125 against 70, mostly five repeated property lambdas; no split seam named, far under 700), @reviewer: APPROVED WITH WORDINGS, applied by the main session: section 7's kill-table row for "`none` applying ΔI" names V1, not H8; `project_structure.md`'s `snow.cpp` entry says the ValueError names the status, and the triangle for a non-finite value. The kernel follows section 4.2 step by step and checks every refusal before it writes a slot; `step()` checks parameters, sizes, step length, non-finite values, negative precipitation in that order, the size check before any array is indexed; the binding passes arrays only and releases the GIL; "triangle i" is the 0-based index; the kill record in `ced40c13` covers all 11 targets; no red-step markers. Not pushed, no CI yet. Suggestions: a refusal case with NaN at triangles 1 and 2 expecting "triangle 1"; one helper for the five array views (about 15 lines); `inline` on `step`'s forward declaration. Idea: PR 3 should name a split point before it starts, since the same 37 % overrun would take its 545 to about 745.

**Code review, round 1 (PR 2, TR1 sun position), 2026-10-13.** Range `cde70666..e08034b7`, 568 net production lines (about 540-560 estimated; tables 301), @reviewer: CHANGES REQUESTED. (1) The sunset test's comment in `tests/cpp/unit/test_radiation_sun.cpp` is wrong: SPA's published sunset, 17:20:19 local, is the UT day's (16 October's evening, 00:20:19 UTC on the 17th), and the reference path crosses h0' at 17:20:18.40 (pvlib's SPA, bisected), so the published sunset is checkable and must be checked: `crossing(e0, kNrelDay - hour, kNrelDay + 2 * hour)` within 2000 ms of `kNrelDay + 1219000`. (2) The branch conflicts with master in `pyproject.toml`: merged by the main session, both paths kept. Wording: the radiation headers, binding and package in `project_structure.md`: added by the main session. Held: equations 4-46 match SPA; ΔT's decimal-year rule (64.508 s for October 2003 by hand); the 195 L/B/R rows and 63 nutation rows identical to pvlib 0.16.1's arrays; std::chrono only; the I/O boundary; the GIL released; layer 1; the terrain-radiation.md paragraph. Ideas: refuse NaT in `to_utc_ms`; stop `sun_positions` truncating float input silently; name pvlib (BSD-3) as the tables' source; the default pressure and temperature belong in TR5.

**Code review, round 2 (PR 2, TR1 sun position), 2026-10-14.** Range `e08034b7..6b795249`, 568 net production lines (unchanged), @reviewer: APPROVED. Round 1's sunset blocker is fixed: the crossing is checked within 2000 ms of 17:20:19, the published sunset of the UT day. The merge keeps both sides of `pyproject.toml` (`git merge-tree`: only the conflict markers differ). The project_structure.md entries for radiation are true. CI is still to be checked after the push.

**Code review, round 1 (PR 4, TR2 receiver frame), 2026-10-14.** Range `e08034b7..59910025`, 216 net production lines (about 350 estimated: C++ 59 against 80, binding and stub 40 against 40, Python 117 against 230, since TR2 covers the projected case only), @reviewer: APPROVED WITH WORDINGS, applied by the main session after merging PR 2's later commits (`791532f6`): `project_structure.md` gains `receivers.hpp`, the receiver check on `radiation.cpp`, `frame.py`, `request.py` and `frame_factors` on `crs.py`; its rule that every number crossing is metres in a projected CRS gains the increment 46 exception (latitude and longitude in degrees and the per-receiver J cross as numbers describing each place); `terrain-radiation.md` marks TR2 built. J is (1/m)·Rot(γ) with true = grid + γ; the exclusion radius a ground distance through J; the zero-area test |n| > 1e-9 times the edge-length product, edges from the first vertex; "receiver i" errors; defaults 64 directions, k 0.13, 20 km and 200 km, centre limb; no CRS crosses; no red-step scaffolding; no acceptance run set for PR 4. Not pushed, CI pending.

**Code review, round 1 (PR 5, TR3a horizon kernel), 2026-10-14.** Range `b1971263..8ecc5596`, 354 net production lines (about 355 estimated: C++ 310 against 300, binding and stub 44 against 55, Python 0), @reviewer: APPROVED WITH WORDINGS: the red step's "Interface assumed" text in both test files becomes "Interface under test" (@tester); `horizon.hpp`'s header says the sample spacing is half the smaller cell side (@developer); TR3a marked built here and in `terrain-radiation.md` (main session); applied in `b8194219`, `e859ec70` and `2eaadb27`. Held: the pruning bound, the early stop over both rasters, far samples from D_near + Δ_far, D_far = D_near with no far raster, rounding half away from zero, pyramid blocks covering their boundary nodes, the pruned search keeping samples through the exhaustive one's functions; @developer's pins (spacing, reach at most 254, the bound's slack, levels 1 and up with a NoData byte per block, d_far checked only with a far raster); the kill record `8ecc5596` covers all 13 targets of section 10.8. Not pushed, CI pending.
