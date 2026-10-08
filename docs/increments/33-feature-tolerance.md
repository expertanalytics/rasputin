# Increment 33: a vertical tolerance that varies with distance to named lines

**Status:** designed (`@architect`, 2026-10-08), round 1's findings answered,
not built. Next: design review round 2. Questions for Ola in section 12; the design is written on their
defaults.

## 1. What Ola asked

Ola, 2026-10-08: "It's feature adaptive refinement, like 1m vtol resolution
close to the railway from Oslo to Bergen, with 20m vtol resolution from 3 km
radially away. It's all about stating my intent, and you using rasputin to
answering the challange." The same day: "I'd like this one to be worked on
tonight".

"vtol" is the vertical tolerance, today's `--tolerance`: the largest allowed
height difference between the mesh and the DEM at any DEM node. Today it is
one number for the whole domain. This increment lets it depend on how far a
place is from a set of lines the user names: tight near them, loose far away.

## 2. Prior art: legacy and literature

### Literature, read

- **Greedy insertion against a sup-norm bound**, the method `refine` already
  is (increment 14): Garland and Heckbert, *Fast Polygonal Approximation of
  Terrains and Height Fields*, CMU-CS-95-181, 1995
  (`https://www.mgarland.org/files/papers/scape.pdf`, read as text). Its
  stopping rule is one global error threshold. Section 3.4, "Product
  Measures", tried weighting a point's importance by "bias measures" such as
  absolute height, and found the results "only slightly poorer than those
  produced by the local error measure" and slower. That is a weight on the
  *choice* of point, not a bound that varies by place; nothing in the report
  bounds the error differently in different places. What this increment
  keeps from it: the worst node of a triangle is still the one inserted; only
  the threshold it is compared with changes.
- **Terrain at a resolution that varies over the domain.** Cignoni, Puppo and
  Scopigno, *Representation and visualization of terrain surfaces at variable
  resolution*, The Visual Computer 13(5):199-217, 1997: the abstract (read at
  `https://science.ulysseus.eu/records/h2rpk-xjm84`) describes a model
  "encoding a history of either refinement or simplification of a
  triangulation"; the search engine's summary of the paper says it extracts a
  surface "at a resolution variable over the domain according to an
  application-dependent threshold function". The full text was behind
  Springer's login and was **not read**, so how the threshold is tested
  against a triangle there is not known to this design. De Floriani, Magillo
  and Puppo, *Building and traversing a surface at variable resolution*, IEEE
  Visualization 1997, 103-110 (abstract read through OpenAlex): "an
  interruptible algorithm for extracting a representation at a resolution
  variable over the surface". Their system paper, *VARIANT*, GeoInformatica
  4(3):287-315, 2000, was found but not read. These are multiresolution
  models queried for a variable threshold; this increment builds one mesh to
  one variable threshold directly, with no hierarchy.
- **View-dependent level of detail**, where the allowed error grows with
  distance from the viewer: Lindstrom et al., *Real-time, continuous level of
  detail rendering of height fields*, SIGGRAPH 1996, 109-118 (abstract read):
  "employs a variable screen-space threshold to bound the maximum error of
  the projected image"; Duchaineau et al., *ROAMing terrain*, IEEE
  Visualization 1997, 81-88 (abstract read): "optimizes flexible
  view-dependent error metrics, produces guaranteed error bounds". The
  structure is the same as here (a bound that is a function of distance to
  something), with a viewpoint where we have lines.
- **Mesh sizing by distance to geometry.** Gmsh's `Threshold` field (manual,
  section "Gmsh mesh size fields", read at
  `https://gmsh.info/doc/texinfo/gmsh.html`): "Return F = SizeMin if
  Field[InField] <= DistMin, F = SizeMax if Field[InField] >= DistMax, and the
  interpolation between SizeMin and SizeMax if DistMin < Field[InField] <
  DistMax", linear unless `Sigmoid` is set, with the input "usually a
  distance"; its `Min` field takes the smallest of several fields. Persson,
  *Mesh size functions for implicit geometries and PDE-based gradient
  limiting*, Engineering with Computers 22(2):95-109, 2006 (found, not read).
  Shewchuk's Triangle (LNCS 1148, 1996) takes a user function deciding
  whether a triangle is too large (`-u`). This increment's ramp and its
  `--tolerance-ramp START END` are Gmsh's `Threshold` applied to the vertical
  tolerance instead of the element size, and its "tightest of the drivers"
  (section 4.6) is Gmsh's `Min`.

### Novelty

None is claimed. A bound that varies over the domain, set by distance to
geometry, is in all four lines of work above. What is specific here is small:
a sup-norm greedy insertion whose per-triangle threshold is the ramp at the
triangle's exact distance to the lines, so the bound holds at every DEM node
(section 3). A search for an error-bounded TIN whose vertical bound follows
distance to a road or railway ("adaptive TIN vertical error tolerance varies
with distance to road railway corridor DEM simplification greedy insertion")
found nothing of the kind; that is weak evidence and is not a claim. If Ola
publishes, this is an application of known ideas, not a contribution.

### Legacy

Nothing to carry across. The legacy tree coarsened with CGAL's
Lindstrom-Turk edge collapse to an edge-count ratio or a maximum size; it had
no tolerance and no distance field:

```
$ git grep -l -i -E "tolerance|distance|max_size" legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/geo_tiff_reader.py
legacy-archive:legacy/rasputin/mesh.py
legacy-archive:legacy/rasputin/triangulate_dem.h
legacy-archive:legacy/rasputin/web_visualize.py
```

The hits are `max_size` (an edge count for the collapse), a "signed
distance" comment in `legacy/rasputin/triangulate_dem.h@legacy-archive:456`, and a viewer's `-dx` "Distance
in meters"; none is an error bound.

## 3. The tolerance field and the rule each triangle is held to

**The ramp.** With near tolerance `N`, far tolerance `F` (today's
`--tolerance`), and a ramp from `S` to `E` metres (`0 <= S <= E`), the allowed
error at distance `d` from the nearest line is

```
t(d) = N                              for d <= S
t(d) = N + (F - N) (d - S) / (E - S)  for S < d < E
t(d) = F                              for d >= E
```

Linear, as Gmsh's default. A step is `S = E`. Ola's case is `N = 1`,
`F = 20`, `S = 0`, `E = 3000`. Requires `0 <= N <= F`: the ramp only ever
tightens near the lines, so `t` never decreases with distance. `N = F` is
today's mesh (section 9, test 3).

*Why linear, not a step or a smooth curve.* A step puts every triangle that
reaches inside `E` at 1 m: the probe (`33-probes/README.md`, item 3) puts the
Geilo-Ål section at 630 967 triangles for the step against 53 154 for the
linear ramp, for no stated need. A sigmoid adds a shape nobody asked for; the
tolerance is already continuous with the linear ramp. A geometric ramp
(`N (F/N)^s`) would hold the tolerance tighter for longer; it is a one-line
change if Ola wants it, and is not the default.

**"Near"** is the distance to the line itself (`S = 0` by default): 1 m on
the line, 1.63 m at 100 m, 4.2 m at 500 m. Holding 1 m out to a band first is
`S > 0`; its cost is in section 12, Q1.

**The rule per triangle: the ramp at the triangle's distance.** Triangle `T`
is allowed `t(d(T))`, where `d(T)` is the smallest distance from any point of
the closed triangle to any line (0 when a line crosses or touches it). The
scan's worst node and its error are computed as today; `T` converges when
that error is at most `t(d(T))`.

This keeps the guarantee honest **at every DEM node**, in the strong form:
for a node `n` in `T`, `d(n) >= d(T)` and `t` does not decrease, so
`t(d(T)) <= t(d(n))`. When refinement stops, every valid node's error is at
most the ramp at that node's own distance. The two alternatives:

- *The tolerance at the worst node* (`t(d(n*))`): the scan would have to
  compare each node's error with that node's own tolerance to find the node
  that breaks its bound, which means a distance query per DEM node per scan,
  in the hot loop that increment 18 made a contiguous row walk. Testing only
  the worst node against its own tolerance is not a guarantee: a node with a
  smaller error closer to the line can break a tighter bound.
- *The tolerance at the nearest point* is the per-triangle rule above; "at
  the centroid" would not be a guarantee (a corner can be nearer the line
  than the centroid).

The per-triangle rule costs one distance query per triangle scanned, not per
node, and leaves the scan's inner loop and its tie-break as they are, so the
point inserted is the one today's code would insert. Its price is some extra
refinement where a large triangle reaches towards the line: such a triangle
is held to the bound at its nearest point. Near the line triangles are small
(at 1 m, 4 944 per km² within 100 m of the Geilo-Ål line: about 200 m² each);
a triangle about 20 m across spans 0.13 m of tolerance there; far out, a
triangle about 200 m across near the 3 km mark spans 1.3 m (19 m over
3 000 m, times 200 m), and is held to the low end of that. The probe's
counts do not include this (section 6).

**Lines are simplified, and the distance is corrected so the bound stays
true.** Python simplifies each line by `margin` metres (Douglas-Peucker) and
passes `margin` with the segments. Every point of the original line is within
`margin` of the simplified one, so `d_original >= d_simplified - margin`;
C++ uses `max(0, d_simplified - margin)`, which is never more than the true
distance, so the allowed error is never more than the true ramp's. Default
`margin` 1 m. *Constant:* metres of a projected CRS; checked on Bane NOR's
network at 69 479 vertices (mean segment 14 m), which it reduces to 13 628
(`33-probes/README.md`, item 1); the tolerance it costs is
`(F - N) margin / (E - S)`, 0.0063 m in Ola's case.

## 4. Where it sits

### 4.1 Data flow

```
--tolerance-near FILE N, --tolerance-ramp S E, --tolerance F      (cli.py)
        |  ToleranceLines (Pydantic, frozen): path, crs, near, start, end
        v
tolerance_field.line_segments(spec, window, mesh_crs)             (Python, pure)
   read_source (feature_input.py, unchanged)  -> line geometries in the file's CRS
   pyproj Transformer(file CRS -> mesh CRS, always_xy=True)
   shapely.simplify(margin, preserve_topology=False)
   keep whole segments whose box meets the window grown by E + margin
        |  float64 array (k, 4): x0 y0 x1 y1, in the mesh CRS; margin
        v
_core.LineTolerance(geometry, segments, near, far, start, end, margin)   (binding)
        |  immutable C++ object, built once per run (once per piece, 23c)
        v
refine(..., field) / refine_strip(..., field) / refine_points(..., field)
   per scanned triangle: allowed = field.at(mesh, t), only when needed (4.4)
```

The C++ core never sees the file, its path or its CRS (CLAUDE.md §2's I/O
boundary): it gets an array of coordinates in the frame it already meshes in.
The mesh CRS is the DEM's, or `--out-crs`'s target grid on the resampled path
(increment 15c), the CRS the domain is already transformed into.

### 4.2 C++: `include/terrain/refinement/line_tolerance.hpp` (new)

```cpp
namespace terrain::refinement {

struct ToleranceRamp {      // metres; 0 <= near <= far, 0 <= start <= end, all finite
    double near, far, start, end;
    double margin = 0.0;    // subtracted from every distance, >= 0
    [[nodiscard]] double at(double distance) const noexcept;  // section 3's t(max(0, d - margin))
};

// A tolerance policy: what each triangle is allowed. lowest() and highest()
// bound at() over every triangle; refine uses them to skip a query (4.4).
template <class P>
concept TolerancePolicy = requires(const P& p, const mesh::LatticeMesh& m, std::uint32_t t) {
    { p.lowest() } -> std::convertible_to<double>;
    { p.highest() } -> std::convertible_to<double>;
    { p.at(m, t) } -> std::convertible_to<double>;
};

struct UniformTolerance {   // today's behaviour
    double value;
    [[nodiscard]] double lowest() const noexcept { return value; }
    [[nodiscard]] double highest() const noexcept { return value; }
    [[nodiscard]] double at(const mesh::LatticeMesh&, std::uint32_t) const noexcept { return value; }
};

class LineTolerance {       // the ramp at a triangle's distance to the segments
public:
    // nullopt with a plain reason for a ramp that breaks its bounds or a
    // non-finite coordinate; zero segments is allowed (every triangle gets far).
    static std::optional<LineTolerance> make(const raster::RasterGeometry& g,
                                             std::span<const std::array<double, 4>> segments,
                                             ToleranceRamp ramp, std::string& why);
    double lowest() const noexcept;   // ramp.near, or ramp.far with no segments
    double highest() const noexcept;  // ramp.far
    double at(const mesh::LatticeMesh& m, std::uint32_t t) const;
    double distance(const mesh::LatticeMesh& m, std::uint32_t t) const;  // d(T), capped at end + margin
private:
    ToleranceRamp ramp_;
    std::vector<Segment2> segments_;  // in the lattice-metre frame (col dx, -row dy)
    noding::BroadPhase index_;        // reused: kernel-free bucket index (noding/broad_phase.hpp)
};
}
```

- **The frame.** Segments are converted once to the lattice-metre frame
  `(x - x_min, y - y_max)`, which is `(col dx, -row dy)`, the frame
  `lattice_frame` already uses, so a triangle's corners convert by two
  multiplications and no large UTM offsets enter the distance arithmetic.
- **Triangle to segment distance** in doubles: 0 if an end of the segment is
  inside the closed triangle or the segment crosses an edge; otherwise the
  smallest of the four point-to-segment distances (the two ends to the
  nearest edge, the three corners to the segment). It is a threshold, not a
  topology decision, so no exact predicate is needed; a rounding error of
  1e-9 m moves the allowed error by 1e-11 m. It is a pure function of its
  inputs, evaluated the same way on every thread, so the output stays
  bit-identical for any thread count.
- **The search.** `noding::BroadPhase` (k = ceil(sqrt(n)) buckets a side over
  the segments' box; its contract: every segment whose closed box meets the
  query box is visited) is queried with the triangle's box grown by `g`,
  starting at `g = E / 16` and doubling, the last value being `E + margin`
  itself (`E/16, E/8, E/4, E/2, E, E + margin`: at most six index queries).
  Any segment within `g` of the triangle has its box inside the grown box, so
  once the best distance found is at most `g` (the distance this query's box
  grew by, not a multiple of it) it is the true minimum and the search stops;
  after the query at `E + margin` finds nothing that near, the distance is
  "far" (capped at `E + margin`). A segment between `g` and `2g` away can lie
  outside the box grown by `g` while a farther one lies inside it (a box
  reaches `g` straight out but about `1.41 g` at its corners), so stopping at
  "best at most `2g`" would miss it (M4). *Constant:* `E / 16`
  assumes `E` in metres of a projected CRS; checked at `E = 3000` on the
  corridor's 12 129 segments, where it is 187 m against buckets of about
  2.6 by 1.2 km (111 a side over the segments' box). Reuse of
  `BroadPhase` is the default; if its header pulls in more than boxes and
  segments, `@developer` writes a 40-line bucket index in
  this header instead and says so.
- **`make` refuses**: `near < 0`, `near > far`, `start < 0`, `start > end`,
  `margin < 0`, any non-finite value, any non-finite coordinate. Each with
  one plain sentence.

### 4.3 The refinement entry points

Each gains an overload taking a policy; the existing signature forwards to it
with `UniformTolerance{options.tolerance}` and keeps its refusals:

```cpp
template <raster::RasterSource R, TolerancePolicy P>
RefineOutcome refine(const R&, const IndexedMesh2&, edges, masks, const RefineOptions&, const P&);
// and the same for refine_points and refine_strip (refine_points.hpp)
```

`options.tolerance` is not read by the policy overload; the binding sets it
to `far` so a message or a record that prints it still says the far value.
`RefineOutcome` gains one figure, `max_error_near` (section 4.5).

### 4.4 Inside the loop: when the allowed error is computed

- **In the parallel scan**, after the scan of `t`: if `max_error <= lowest()`
  the triangle converges and if `max_error > highest()` it splits, whatever
  its distance; only in between is `allowed = at(m, t)` computed, and stored
  in a vector beside `results` (one `double` per slot). The serial phase
  compares `max_error` with it. With `UniformTolerance`, `lowest() ==
  highest()`, so no query ever runs and the comparisons are today's; and the
  `allowed` vector exists only for a policy whose bounds can differ, chosen
  at compile time (`if constexpr` on the policy type, or a policy trait), so
  the default path allocates nothing it did not allocate before.
- **Constraint feet** (increment 20b's `eps(n) = clamp(tol / G, ...)`): `tol`
  becomes the triangle's allowed error. When the scan skipped the query
  (error above `highest()`), the serial phase computes it there, only when a
  constraint lies within the foot cap (the existing early test); the triangle
  is unchanged at that point (a touched triangle is skipped before).
- **The edge strip and the final check** (`refine_points`, `refine_strip`,
  increments 15f and 15c-2): the same rule per triangle holding a check
  point, the same laziness.
- **The quality start** (increment 20, 20c) reads no tolerance and is
  unchanged.

### 4.5 Python

`src_python/tin_engine/tolerance_field.py` (new):

```python
class ToleranceLines(BaseModel, frozen=True):
    path: Path
    crs: str | None          # --tolerance-near-crs; None: the file's own
    near_m: float            # N
    start_m: float           # S
    end_m: float             # E
    margin_m: float = 1.0

def line_segments(
    spec: ToleranceLines, window: Box, mesh_crs: str
) -> tuple[npt.NDArray[np.float64], float]:
    """(k, 4) segments in mesh_crs, and the margin they were simplified by:
    every LineString and MultiLineString of the file, transformed, simplified
    by margin_m, kept whole when its box meets window grown by end_m +
    margin_m. A line that simplification empties (a closed loop smaller than
    the margin) is kept as a zero-length segment at its first vertex, never
    dropped. Polygons and points refused."""
```

It reads through `feature_input.read_source`, so GeoJSON (with or without a
`crs` member), GeoPackage and GML 2 work as `--features` files do, with the
same CRS checks. The margin returned is the one the binding passes to `LineTolerance`, so the
correction C++ applies is the simplification Python did. `cli.py` builds the spec, `_dem_mesh` builds the field once
from the tile's geometry and passes it to `refine`, `edge_strip.run` and
`final_check.run`. The record (increment 25) gains `tolerance_near_m`,
`tolerance_ramp_m` ("0 to 3000"), `tolerance_lines` (the file's name and its
segment count after simplification) and `max_error_near_lines_m`: the largest
error over the final triangles whose allowed error is `N` (with `S = 0`, the
triangles a line crosses or touches). `tolerance_m` stays the far value.
Without `--tolerance-near`, none of these appear and nothing is built.

### 4.6 Locality, parallel refinement, and pieces

- The field is read-only after `make`, shared by every scan thread; a query
  touches only the buckets near one triangle. No global operation is added to
  refinement.
- **Per tile or piece** (increment 23's pieces, when 23c-2 lands): each piece
  builds its own field from the segments whose boxes meet the piece's window
  grown by `E + margin`. Beyond that distance a segment cannot change any
  triangle's allowed error, so the field per piece equals the global one on
  that piece. Segments are kept whole, never cut at the window, so two pieces
  beside a seam hold the same segments near it.
- **Seam identity comes from the seam pass, not from the piece fields.** Each
  piece converts the segments into its own lattice-metre frame, measured from
  the corner of its own raster (increment 23's N17,
  `docs/increments/23-basin-scale.md`, section N17), so two pieces' fields
  give the same distance to rounding, not bit for bit. That is enough inside
  a piece, which only its own field ever judges; a 1e-9 m difference in
  distance moves the allowed error by about 1e-11 m. Anchoring every piece's
  segments to one shared origin was not chosen: the triangle corners would
  then need the piece's offset added, which rounds in the same way unless the
  offset and the cell size are exact binary fractions. At a seam, both pieces
  call `refine_seam` on the same strip raster (23's N17: the strip is a
  function of the edge alone), so the seam pass builds its field in the
  strip's frame from the segments selected by the strip's window, the same
  inputs for both pieces, and gets the same bits. That needs the allowed
  error of a piece of seam, `t(d(segment piece))`, a segment-to-segment
  distance: a `LineTolerance::at(Point2 a, Point2 b)` of about 15 lines. It
  is **not in this increment**, because nothing on master calls
  `refine_seam` yet; it goes in with whichever of 23c-2 and 33 lands second.
  Its test checks the seam pass, not the piece fields: two pieces beside a
  seam, each running the seam pass on it, get the same allowed error for
  every piece of seam, bit for bit; and the two pieces' own fields agree on
  a triangle near the seam to within 1e-9 m of distance.
- 23c's memory estimate (`decompose.bytes_per_node`, by tolerance) would
  take, per piece, the node count in each distance band times the band's
  `b(t)`; also with 23c-2.
- **A second driver, slope** (an idea of Ola's, not built here): a policy
  whose `at(m, t)` reads the DEM under `t` (its steepest cell, say) is local
  and per tile in the same way. The policy refine takes would then be the
  smallest of its drivers (`lowest` and `highest` the smallest of theirs), as
  Gmsh's `Min` field; nothing in this design changes for it.

### 4.7 Against the four criteria

1. *Data and execution:* `ToleranceLines` and `ToleranceRamp` are data; the
   field is built from them and only queried.
2. *State:* no global state; the field is immutable and thread-safe; the
   scan stays pure.
3. *Dependencies:* none added (shapely, pyproj and NumPy are in the stack;
   the bucket index exists).
4. *Async:* the field is built before `refine`, in the same thread that
   calls it, and the binding releases the GIL as it does today.

## 5. Interface

```
rasputin mesh --dem DTM10 --domain corridor.geojson \
    --tolerance 20 --tolerance-near bergen_line.geojson 1 --tolerance-ramp 0 3000
```

- `--tolerance F`: unchanged; with `--tolerance-near` it is the far value.
- `--tolerance-near FILE N` (a file and a number): "Hold the mesh to N metres
  on these lines (.geojson, .gpkg or .gml). Needs --tolerance-ramp."
- `--tolerance-near-crs CRS`: the file's CRS when it does not say, as
  `--features-crs`.
- `--tolerance-ramp START END`: "From START metres from the lines, where N
  stops, to END metres, where --tolerance starts." Required with
  `--tolerance-near`; there is no natural default, so none is invented.
- Refusals, each a usage error in plain words: `--tolerance-near` without
  `--tolerance` or without `--dem`; `--tolerance-ramp` without
  `--tolerance-near`; N above F; START above END; a negative or non-finite
  number; a file with no lines ("bergen_line.geojson has no lines; polygons
  and points are not used here"). Lines wholly beyond END of the window are
  not an error: stderr says "no tolerance lines within 3000 m of the domain;
  every triangle is held to 20 m".
- **The lines do not become constraint lines** (default, Q3). A railway's
  centre line is not a break line of the terrain at 10 m (the track bed is an
  embankment or a cut about that wide), forcing mesh edges along it adds
  vertices at bilinear heights, and a tunnel's line would be pulled to the
  surface. A user who wants both gives the same file to `--features` as
  well. One flag, one effect.
- Single feature file in this increment; several files with different `N`
  are the `Min` of several fields (4.6) and a later, small step.

## 6. The demonstration case: the Bergen Line

**The line.** Bane NOR's *Jernbane - Banenettverk* (Geonorge metadata
`c3da3591-cded-4584-a4b1-bc61b7d1f4f2`): centre lines of every railway link,
with name (`banenavn`), status (`banestatus`, `I` in operation) and medium
(`medium`, `U` tunnel). Licence: *Norsk lisens for offentlige data (NLOD)*
1.0, `http://data.norge.no/nlod/no/1.0`, whose section 5 asks for the line
"Contains data under the Norwegian licence for Open Government data (NLOD)
distributed by [name of licensor]", a link to the licence, and a note that
the data were changed.
National GML 3.2 in EPSG:25833, 12 MB zipped. **Preferred to
OpenStreetMap** (default, Q6): the owner's data, attributes for status and
tunnels, and NLOD has no share-alike clause, where the ODbL's section 4.4 says
"Any Derivative Database that You Publicly Use must be only under the terms
of" the ODbL or a compatible licence. rasputin's GML reader reads GML 2, not
this file's GML 3.2 (`posList`, `gml:id`), so a preparation script converts
the selected links to GeoJSON in EPSG:25833 with a `crs` member.

**The route** as the trains run (default, Q5): Oslo S to Drammen
(Drammenbanen, Askerbanen), Drammen to Hokksund (Sørlandsbanen), Hokksund to
Hønefoss (Randsfjordbanen), Hønefoss to Bergen (Bergensbanen, 371.7 km by the
probe), links in operation, tunnels included (default, Q2). Bane NOR's own
"Bergensbanen" starts at Hønefoss.

**The DEM.** DTM10 at `../rasputin_data/DTM10_UTM33_20260925`: the probe's
5 km corridor of Drammenbanen, Randsfjordbanen and Bergensbanen touches 14
tiles and is fully covered (Drammen to Hokksund was not in it). DTM1 is not
held locally; the whole corridor at 1 m would be about 4.7e9 nodes inside
the domain and a canvas of the domain's box far beyond memory, so it waits
for increment 23's pieces (Q4).

**The domain.** The line buffered by 5 km (the ramp's 3 km plus 2 km of
20 m ground to show the change), simplified by 20 m, one polygon (Q7). The
probe's corridor without Drammen-Hokksund was 4 196 km² with 1 154 vertices
and three holes; the three lines' corridor was 4 689 km² before
Drammen-Hokksund is added.

**Expected size and time** (probe, `33-probes/README.md`; Mac on AC):

| case | uniform 20 m | uniform 1 m | the ramp (estimate) |
|---|---|---|---|
| Geilo-Ål section, 267 km² | 5 967 triangles, 0.8 s | 970 037, 1.5 s | about 53 000 |
| Hokksund-Bergen corridor, 4 196 km² | 229 611, 2.5 s, 2.4 GB | 21 669 768, 37 s, 11.8 GB | about 1.37 million |

At 1.37 million triangles, refine's measured rate at 1 m on the corridor
(21.7 million in 17.9 s) gives about 1.1 s; the field's queries add 0.9 to
3.1 s of wall time on 10 threads (section 10, from `33-probes/query_cost.py`);
with the 1.6 s decode the whole corridor should take about 4 to 7 s. The
triangle count is an estimate, likely low: section 3's extra refinement where
a large triangle reaches towards the line is not in it, and was not
measured.

**The first, small case: Geilo to Ål**, 23 km of x (127 000 to 150 000 in
EPSG:25833) of Bergensbanen through Hallingdal, the 5 km buffer cut to that
range: 267 km², DTM10 tiles `6701_2` and `6701_3`, about one second. It is
the quick check's new case (section 10) and the first run Ola looks at. The
tests use synthetic DEMs and the repository's fixtures, not this data.

**Where things go.** The Bane NOR file and the prepared GeoJSON under
`../rasputin_data/banenor_banenettverk/` (not committed); the preparation
script and its README (source URL, licence, the selection) under
`docs/benchmarks/<date>/bergen-line/`, with the run's `--stats` and record.
A `rasputin fetch` source for Bane NOR is a later step (Ola's wish that
rasputin fetch its own data); the script is the stopgap.

## 7. Guarantees and the checks that enforce them

| # | guarantee | checked by (section 9) |
|---|---|---|
| G1 | Without `--tolerance-near` every output byte is today's | test 1, the bench run |
| G2 | With `N = F`, or with no segment within reach, the mesh is the `--tolerance F` mesh, bit for bit | tests 2, 3 |
| G3 | Every valid DEM node's error is at most `t(d(n))`, `d` to the original lines: after the transform to the mesh CRS, before simplification | test 5 |
| G4 | Every edge-strip and final-check point's error is at most `t` at its triangle's distance | test 8 |
| G5 | The output is the same for any thread count | test 9 |
| G6 | The lazy query changes nothing | test 7 |
| G8 | The field reaches every comparison: a ramp that holds every triangle at `N` gives the `--tolerance N` mesh, bit for bit, feet included | test 3 (b) |
| G7 | Simplification never makes the allowed error larger | test 11 |

## 8. Degeneracies

- A zero-length segment (two equal points): a point; its distance is the
  point's.
- A segment along a triangle edge, or through a vertex: distance 0.
- A line wholly outside the window grown by `E + margin`: dropped in Python;
  zero segments is a valid field (every triangle gets `F`).
- `S = E`: a step; the ramp formula is not evaluated between (no division by
  zero).
- `N = 0`: allowed, as `--tolerance 0` is; refinement terminates for the same
  reason (a finite lattice).
- Lines in a geographic CRS: transformed to the mesh CRS first; the mesh CRS
  is always projected (a geographic DEM needs `--out-crs`).
- A line with a NaN or infinite coordinate: refused by `make`.

## 9. Tests `@tester` writes red first

C++ (`tests/cpp/unit`, `tests/cpp/property`):

1. **Today's mesh, bit for bit.** On two DEM fixtures with constraints and
   feet on, `refine(..., options)` and `refine(..., options,
   UniformTolerance{t})` give identical vertices, z, triangles, edges, masks
   and counters; the same for `refine_points` and `refine_strip`.
2. `LineTolerance` with zero segments equals `UniformTolerance{F}`, bit for
   bit, for all three entry points.
3. Two equalities with segments, bit for bit, on all three entry points, on
   test 1's fixtures with constraints and feet on:
   (a) `LineTolerance` with `N = F` equals `UniformTolerance{F}`;
   (b) `LineTolerance` with `N` well below `F` and `S` beyond the largest
   distance any triangle can have to the lines (a line crossing the domain,
   `S` above the domain's diagonal), so the ramp holds every triangle at `N`,
   equals `UniformTolerance{N}`; the `LineTolerance` run has
   `options.tolerance = F`, as the binding sets it. (b) is the case where the
   field, not `options.tolerance`, must reach the foot's epsilon and
   `refine_points`' comparison (M6, M7).
4. **The ramp:** `t` at 0, `S`, the middle, `E` and beyond; the step `S = E`;
   the margin; never decreasing (a sweep of 10 000 distances).
5. **The guarantee, oracle independent of the field.** A synthetic DEM (a
   sum of ridges) and three lines (one crossing the domain, one along a
   constraint, one ending inside a triangle); after refine, every valid DEM
   node's `|z - plane|` is at most `t` of its brute-force distance to the
   original segments, computed in the test.
6. **Distance:** triangle and segment crossing with both ends outside (0), an
   end inside (0), touching an edge (0), along an edge (0), parallel at a
   known offset, nearest at a corner, nearest at an end, a zero-length
   segment; the search's stop: a small triangle, one segment about `1.3 g`
   away diagonally (inside the box grown by `g = E / 16`) and one about
   `1.1 g` away straight out (outside it), where the answer is the second
   (M4); then 10 000 random triangles against random segment sets: the
   indexed distance equals the brute-force minimum exactly, with and without
   the cap.
7. **Laziness changes nothing:** a test-only policy that reports the full
   range `[0, inf)` as its bounds (so every triangle is queried) gives the
   same output as `LineTolerance`, bit for bit.
8. **The strip:** a constraint edge running along a line, `N` well below
   `F`; every strip point within `S` of the line ends within `N`.
9. **Threads:** 1 and 8 threads give the same mesh with a field.
10. `make`'s refusals, one each.

Python (`tests/python`):

11. `line_segments`: transformation from EPSG:4326 to a UTM CRS; MultiLineString
    split into segments; polygons and points refused in words; a line 2 999 m
    outside the window kept and one at 3 002 m dropped (`E = 3000`, margin 1);
    on 1 000 random points, the distance to the returned segments less the
    margin `line_segments` returns is never above the distance to the
    original line (transformed, not simplified); the same for a closed loop
    smaller than the margin, which comes back as at least one segment, not
    as nothing.
12. CLI: each refusal of section 5, in its words; the stderr line for lines
    out of reach; the record's four new fields; `max_error_near_lines_m <= N`.
13. The resampled path (`--out-crs`): the lines are transformed to the
    target grid's CRS, not the DEM's (a line given in the DEM's CRS ends up
    where the mesh's own coordinates put it).
14. `rasputin mesh` with the flags on a fixture: more vertices within 100 m
    of the line than the `--tolerance F` mesh has there (greedy insertion is
    not monotone in the tolerance, so no count against the `--tolerance N`
    mesh is asserted).

**Invariant-critical suite: tests 1-8, 11** (3 includes its new part (b);
8 is the strip). Mutation targets the kill record
must cover:

- M1 the ramp interpolates from `E` instead of `S` (or drops the margin);
- M2 `d(T)` measured from the centroid instead of the closed triangle;
- M3 a segment crossing a triangle with both ends outside gets a positive
  distance;
- M4 the search stops when the best is at most `2g` instead of `g`, the
  distance the box grew by (killed by test 6's stop case);
- M5 the laziness swapped: converge when the error is at most `highest()`;
- M6 the foot's epsilon uses `F` instead of the triangle's allowed error
  (killed by test 3 (b));
- M7 `refine_points` compares with `options.tolerance` instead of the field
  (killed by test 3 (b));
- M8 Python drops the margin correction;
- M9 Python selects with the window not grown by `E`.

## 10. Size, split point, speed

| part | counted lines (estimate) | basis |
|---|---|---|
| `line_tolerance.hpp`: ramp, policy concept, uniform policy, distance, search, `make` | 110-140 | `seam.hpp` counts 102, `constraint_points.hpp` 133 |
| `refine.hpp`: the policy overload, lazy allowed, foot epsilon, `max_error_near` | 30-40 | the parts listed; no close comparable |
| `refine_points.hpp`: the policy overloads, lazy allowed | 25-35 | 23b's green commit added 22 to it for one entry point |
| binding (`LineTolerance` class, optional field argument on three functions) and `_core.pyi` | 60-85 | 23b's green commit: 46 in `bindings/core.cpp` and 31 in `_core.pyi`, 77, for one entry point |
| `tolerance_field.py` | 45-65 | `outline.py` counts 50, `decompose.py` 55 |
| `cli.py` (three flags, refusals, passing the field, record) | 40-60 | the parts listed; no close comparable |
| **total** | **310-425** | |

The counts are `tools/count_loc.py`'s rule applied to the files at
`700f57a6` (`counted_lines` on each file whole), and to 23b's green commit
(`python3 tools/count_loc.py 3c464ec3~1 3c464ec3`).

Under 700 in one PR. **Split point** if the red suite pushes the estimate past
600: PR A the C++ header, the entry-point overloads, the binding and stub
(nothing a user sees changes); PR B `tolerance_field.py`, the CLI and the
record.

**Cost of every added step.**

| step | quick check's cases, default flags | Hokksund-Bergen corridor with the ramp | basis |
|---|---|---|---|
| read the line file (`read_source`, GeoJSON) | not run | under 0.5 s (69 479 vertices) | estimate; `json` plus `shapely.shape` |
| transform (`pyproj.Transformer.from_crs(..., always_xy=True).transform`, PROJ's default operation, no tolerance) | not run | none if prepared in EPSG:25833; about 0.02 s for 69 479 points otherwise | estimate |
| simplify (`shapely.simplify(geom, 1.0, preserve_topology=False)`, GEOS Douglas-Peucker, tolerance 1 m) | not run | 0.0 s measured (to 13 628 vertices) | probe item 1 |
| select segments (NumPy box test) | not run | milliseconds | O(n) |
| build the field (`BroadPhase` over about 12 000 segments) | not run | milliseconds | O(n) |
| per-triangle query in the scan | none: `UniformTolerance` never queries | 0.7 to 2.7 µs per index query (one to four occupied buckets of about 45 segments, about 15 ns per distance); a triangle beyond `E` runs all six queries of the doubling search, one within `E/16` runs one; per scanned triangle 1.8 to 6.4 µs, weighted by triangles per distance band; 3.5 scans per final triangle, 1.37 million final triangles: 8.8 to 30.8 s of CPU, 0.9 to 3.1 s of wall time on 10 threads | estimate, not measured (no C++ exists): `33-probes/query_cost.py`, its assumptions stated there; errs high, as it does not count the queries the laziness skips |
| `max_error_near_lines_m` at the end | not run | one query per final triangle, under 0.3 s wall | estimate |

**Library calls**, with method and tolerance: `pyproj.Transformer.from_crs(
src, dst, always_xy=True)` (PROJ's best available operation; no tolerance);
`shapely.simplify(..., tolerance=margin_m, preserve_topology=False)`
(Douglas-Peucker in GEOS); `shapely.get_coordinates` (exact copy). No GEOS
distance or buffer runs in the product; the distances are the C++ header's.

**Speed judgment.** This touches refine code, so `@perf`'s acceptance
applies (`docs/increments/README.md`): `tools/bench.py`'s 1 m benchmark and
thread-scaling sweep at default flags must show no change beyond the bench's
band, with the mesh's SHA-256 unchanged (G1). The quick check
(`tools/bench_quick.py`) gains the Geilo-Ål case with the ramp
(`docs/benchmarks/quick/cases.toml`), which has no baseline at first; `@perf`
records one. Accepted when the default-flag runs are unchanged and the
Geilo-Ål ramp run's refine time per output triangle is at most 4 times the
uniform 1 m run's on the same section. The basis: the uniform 1 m corridor
refines at 0.83 µs of wall time per output triangle (17.89 s for 21.7
million), and the queries add 0.64 to 2.26 µs per output triangle (3.5 scans
at 1.8 to 6.4 µs, on 10 threads), so the estimate is 1.8 to 3.7 times; the
first draft's 1.5 was below the estimate's own low end. Above 4 times, the
query is what `@developer` cuts first (a bucket index sized to `E / 16`, or
asking only whether a segment lies within the distance the error needs). The
corridor run is timed once for the record, not judged.

## 11. Not in scope

The seam pass's allowed error and 23c's memory estimate (4.6); several
`--tolerance-near` files; polygons as near features (a town, a lake); slope as
a driver; a `rasputin fetch` source for Bane NOR; DTM1.

## 12. Questions for Ola

Each with the default this design is written on.

- **Q1. Hold 1 m for a band before the ramp starts?** Estimated from the
  measured uniform meshes on the corridor (`33-probes/README.md`, item 4): linear from the line out to 3 km, about 1.37 million triangles; 1 m
  held to 100 m first, about 1.75 million; to 500 m, about 3.1 million.
  *Default: no, the ramp starts at the line (`--tolerance-ramp 0 3000`); the
  flag lets you choose per run.*
- **Q2. Tunnels:** keep the tunnel sections in the line, so the ground above
  a tunnel is also held to 1 m? About 168 of Bergensbanen's 743 km of drawn
  links are tunnels. *Default: keep them; the line is the line as Bane NOR
  draws it.*
- **Q3. Should the line also become a mesh edge?** *Default: no; give the same
  file to `--features` too if you want that.*
- **Q4. DEM:** DTM10 for the whole corridor now, and the Geilo-Ål section on
  DTM1 later, once the 1 m tiles are fetched; the whole corridor on DTM1 waits
  for the basin-scale pieces (increment 23). *Default: DTM10 now.*
- **Q5. Route:** Oslo S to Bergen as the trains run (via Drammen, Hokksund and
  Hønefoss), or only Bane NOR's "Bergensbanen" (Hønefoss to Bergen)?
  *Default: as the trains run.*
- **Q6. Source:** Bane NOR's network (open licence, the owner's data) or
  OpenStreetMap (share-alike licence)? *Default: Bane NOR.*
- **Q7. Corridor width:** 5 km each side, so 2 km of 20 m ground shows beyond
  the ramp. *Default: 5 km.*

## 13. ROADMAP

Row 33, added with this design: "designed, not built".

## Review

Increment 33 design review round 1, 2026-10-08, @reviewer (`4cd7e050..c22fb74c`, prose and probes only, 0 counted lines): CHANGES REQUESTED. Five fixes before the red step, all in /Users/skavhaug/projects/rasputin/.claude/worktrees/vtol-field/docs/increments/33-feature-tolerance.md: (1) §4.6 says pieces compute the same distance bit for bit, false while segments are converted into each piece's own frame (/Users/skavhaug/projects/rasputin/.claude/worktrees/vtol-field/docs/increments/23-basin-scale.md@c22fb74c:2476-2487, N17); anchor the segment frame to a shared origin, or say the fields agree to rounding and seam identity comes from the seam pass; (2) M4 as worded is an equivalent mutant (a stricter stop test); reword it as stopping when the best is at most 2g; (3) no listed test kills M6 or M7: add a test where a `LineTolerance` ramp puts every triangle at N and equals `UniformTolerance{N}` bit for bit on all three entry points with feet on, and add it and test 8 to the invariant-critical list; (4) §10's per-query cost contradicts its basis (0.7-2.7 µs, not "at most about 1.5 µs"; 8.2 s, not "under 8 s") and leaves out the doubling search (about 6 index queries per triangle beyond E); restate, and carry into §6's "4 to 6 s"; (5) §6's "a lower bound on triangles" is not shown; say "an estimate, likely low", and drop "Measured estimate" in Q1. Checked and sound: the per-triangle rule keeps the guarantee; `max(0, d_simplified - margin)` never exceeds the true distance; segment selection across pieces; default runs byte-identical (tests 1, 2, G1); literature quotations checked word for word; the seven questions have defaults; library calls named; code claims exist. Suggestions: give §10's line counts a basis; say the default path allocates no `allowed` vector (`if constexpr`); test 11 uses `line_segments`' margin and adds a closed loop smaller than the margin; G3's "original lines" means after the CRS transform.
