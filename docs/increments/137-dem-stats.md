# Increment 137: each triangle's plan area in the mesh file

Status: **designed** (`@architect`, 2026-11-04, on public master `98fe80c3`;
branch of worktree `dem-stats`). Full lane (a per-triangle number a
downstream reader takes as input: a wrong value gives wrong numbers without
a crash). Python only: no C++, nothing under `include/terrain/refinement/`
or `include/terrain/mesh/`, so no `@perf` acceptance run is required
(`docs/increments/README.md`, *Acceptance*); section 7 asks for timings,
reported and not gated, as `perf.md` section 4 says for an option off by
default.

## 0. The plain picture

- **Today** a `.gpkg` mesh already carries each triangle's plan area
  (`area_m2` in the `triangles` and `cells` layers). A `.vtk` or `.ply`
  mesh does not; a reader recomputes it from the corners.
- **After**, `rasputin mesh --plan-area` writes the same number into a
  `.vtk` (cell array `area_m2`) and a `.ply` (face property `area_m2`), for
  one mesh and for every piece of a cut run, and `io.vtk_read.read_mesh`
  returns it. Without the option every file is byte for byte today's.
- **Not built:** statistics of the DEM pixels under each triangle (the
  spread of the DEM about the triangle's face, its higher moments, a
  quadratic-fit residual, hollow and knoll shares, the 3D surface area).
  Section 4 says why, and holds a design for the pixel walk should a user
  need it (question 1).

## 1. Prior art: legacy and literature

*Literature.* A triangle's plan area is half the cross product of two edge
vectors (the shoelace formula); nothing to cite. It is the area in the
plane of the mesh's projected CRS, not on the ellipsoid: a conformal
projection scales areas by k², and in a UTM zone |k² − 1| stays under about
0.2 % (k from 0.9996 on the central meridian to about 1.001 at the zone's
edge; Snyder 1987, *Map Projections: A Working Manual*, USGS PP 1395,
section 8). The `.gpkg` writer has written this plane area since increment
36; this increment does not change its meaning.

The statistics not built (section 4) have their literature: Jenness 2004
(*Calculating landscape surface area from digital elevation models*,
Wildlife Society Bulletin 32(3): 829-839) for surface area over plan area
as roughness; Grohmann, Smith and Riccomini 2011 (*Multiscale analysis of
topographic surface roughness in the Midland Valley, Scotland*, IEEE TGRS
49(4): 1200-1213) for the spread of residual heights after a fitted surface;
Pineda 1988 (*A parallel algorithm for polygon rasterization*, SIGGRAPH,
Computer Graphics 22(4): 17-20) for edge functions with a tie rule for a
pixel on an edge, which section 4.2 makes exact by a symbolic shift
(Edelsbrunner and Mücke 1990, *Simulation of Simplicity*, ACM TOG 9(1):
66-104). No novelty is claimed; none searched.

*Legacy.* `git grep -n -i -E "area|rugos|rough" legacy-archive -- legacy`
finds no per-triangle area; `legacy/rasputin/writer.py` wrote named
per-triangle arrays (`add_cell_data(data, name)`) through meshio. Carried
across: a per-triangle quantity is a named cell array of the mesh file.

*Ideas file.* Overlaps "A real storage format for meshes" (one more VTK
array for its reader to carry), "A terrain mesh that follows the physics,
not only the height" (its "carry each triangle's sub-grid distribution"),
"Mean, RMS and percentiles of the vertical error in --stats" and "unbiased
height statistics under the sup-norm" (its candidate 3, the mean signed
error per triangle, is section 4.2's `dem_misfit_mean_m`). Takes up none.

## 2. Measured

Scratch probes in this session, master `98fe80c3`, M1 Max on mains power,
`rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif
--tolerance 1` (EPSG:25833, 10 threads):

- **The mesh:** 464,290 triangles, 872 constraint lines, `.vtk` 16,869,273
  bytes (36.3 per cell). Plan areas from 50 m² to 96,000 m², median 100 m²;
  their sum 2,546.01 km² (the 50.5 km square less the NoData trim).
- **Cost:** `landcover.plan_areas` over the 464,290 triangles took 0.032 s.
  The array adds 8 bytes per cell (lines included), 3.72 MB, +22 % on this
  file.
- **One number everywhere:** `landcover.plan_areas`, `cells.geometry`'s
  area and `write_triangles`' inline formula agree bit for bit on all
  464,290 triangles (`np.array_equal`).
- **Pieces** (`tests/python/pieces_fixtures.py`: `hills()`, 41 × 41 nodes
  at 10 m, domain `SQUARE`, 115,382.4 m², tolerance 1). Uncut: 1,211
  triangles, the areas sum to the domain's area to −1.1e-12 relative.
  `--pieces 4`: 4 files, 1,163 triangles, sum 115,382.386 m², −1.2e-7
  relative. The pieces overlap by 1.5e-11 m² and their union is their sum;
  the 0.014 m² is the outline as the cut run nodes it (0.074 m² of the
  domain outside the union, 0.060 m² of the union outside the domain), not
  the areas.
- **Repeatable:** two runs of the fixture command give identical `.vtk`
  files; the file holds no version or time, so "byte for byte today's" is
  testable by digest.

## 3. Design

**The option.** `--plan-area`, a flag, off by default. Help: "Write each
triangle's plan area, area_m2 (square metres in the plane of the mesh's
CRS), into a .vtk or .ply; a .gpkg always has it." With `.gpkg` it changes
nothing. It is no refinement option: it works with `--stride`,
`--tolerance`, `--triangles` and cut runs alike, and is not on the cut
run's list of refused options.

**The number.** `landcover.plan_areas(vertices[:, :2], triangles)` of the
mesh as written (after the trim, in written triangle order), float64.
The mesh's CRS is always projected in metres: `dem_input.py` refuses an
`--out-crs` that is not, and a DEM in feet is read as its metre twin
(`crs.metre_twin`). So the unit is m² on every path, and the file's `crs`
string names the plane.

**The files.**
- `.vtk`: `write_vtk(..., triangle_areas=None)`; given, the cell array
  `area_m2` of type `double` in the `features` FIELD block, one entry per
  cell, **0.0 on every constraint line** (every cell array covers every
  cell, as `land_cover_code` does). Placed after the existing cell arrays.
- `.ply`: `write_ply(..., face_areas=None)`; given, the face property
  `double area_m2` after the existing face properties, binary and ASCII.
- `.gpkg`: unchanged.
- No record entry: the array's presence says it was asked for, and a
  `.gpkg` run with the flag stays the same file as one without.

**The reader.** `MeshArrays` gains `area_m2: NDArray[float64] | None`,
`(T,)`, the triangle entries only. `read_mesh` reads a `double` array named
`area_m2`; another type, or a count other than lines + triangles, is a
`ValueError` naming the array, as for `land_cover_code`.

**Cut runs.** `_write_piece` gets the flag and passes each piece's
`plan_areas` to its writer. Each piece's array covers that piece's
triangles.

**Stable names.** `area_m2` is the downstream contract: the same name the
`.gpkg` layers use for the same number. `--plan-area` and `area_m2` are
renamed only with a deprecation by a later increment.

**Exclusions.** No 3D surface area (it is plan area over the cosine of the
triangle's slope, from the same corners); no DEM pixel statistics
(section 4); no change to `.gpkg`, `--stats` or `--record`.

## 4. The statistics not built, and the walk they would need

### 4.1 Which earn a place

Judged by what reads them and what they say beyond plan area:

| statistic | needs | verdict |
|---|---|---|
| plan area | corners | **built**; the one a downstream reader asks for |
| 3D surface area, rugosity | corners; DEM | out: from the corners it is the slope again; from the DEM, a private probe found the micro-relief adds under 1 % of area (median ratio 1.00-1.01 in every 100 m height band) |
| misfit count, mean, SD | DEM walk | held (4.2): no reader yet. A private probe on a 305 km² mesh (10 m DEM, tolerance 10 m) found the SD's median 1.93-2.08 m in all twelve 100 m height bands (p10-p90 about 0.9-3.7 m): on a tolerance mesh it is set by the tolerance |
| skewness, kurtosis, max \|misfit\| | DEM walk | out: the same probe found them level across bands and no better than plan area at sorting the outcomes it tested |
| residual after a quadratic fit | DEM walk, 12+ pixels | out: correlates 0.71 with the SD, adds nothing in that probe |
| hollow and knoll shares | a window past the triangle | out: needs DEM rows outside the triangle, so a halo per piece |

### 4.2 The held design (built only if question 1 says so)

- **Fields:** `dem_nodes` (`unsigned_int`, valid DEM nodes owned),
  `dem_misfit_mean_m` and `dem_misfit_sd_m` (`double`; the misfit is DEM
  height minus the triangle's face, the face the plane through its written
  corner heights; SD the population SD; mean and SD NaN when `dem_nodes` is
  0). A count and a mean let any union of triangles get its exact mean
  misfit.
- **Which nodes belong to a triangle:** node p (pixel centre) belongs to
  the triangle holding p shifted by (ε, ε²) in world axes, ε → 0⁺. Per
  edge a → b of a counter-clockwise triangle, in lattice (col, row): p is
  in when orient(a, b, p) > 0, or = 0 and the edge is *inclusive*,
  `b.row > a.row`, or `b.row == a.row` and `b.col > a.col`. Exact: orient
  is the existing exact predicate, and the edge test compares coordinates.
  Each node in the meshed region belongs to exactly one triangle, vertices
  included, across piece seams too, since the rule reads geometry, not
  triangle numbers (scratch brute force: 300 random lattice triangulations
  with off-node fans, 0 violations; flipping the column tie or dropping
  the edge case gives 3,000 and 10,690).
- **Where:** C++, a new header beside `refinement/measure.hpp` reusing
  `mesh::for_each_row_span` and `raster::for_each_row_segment`, plus the
  three vertices, with Welford's update in row-major order, one slot per
  triangle: bit-identical for any thread count. Called on the trimmed mesh
  while the tile is held (in `_refined`, trim moves before `del tile`).
  Refused at first with `--out-crs` (the tile is gone before the final
  check) and in cut runs until a piece's window is shown to hold every
  node its triangles own.
- **Cost:** the error summary (`--rms`), the same walk, measured 0.079 s
  for 25,233,600 nodes at 10 threads on the fixture above; the Welford
  update is estimated at no more than doubling it. A 637,000 km² basin on a
  30 m DEM has 7.1e8 nodes, 28 times the fixture: about 2-5 s at 10
  threads, and 16 bytes per triangle.
- **Invariants and tests:** `dem_nodes` sums to the valid DEM nodes whose
  shifted point lies in the written mesh; hand-worked 3 × 3 and 4 × 4 DEMs
  with both diagonals, an off-node vertex and a NoData node; thread counts
  1, 2 and 7 give equal bytes; the brute force above as the oracle.
  Estimated +250 counted lines.

## 5. Invariants

1. **Off is today:** without `--plan-area`, `.vtk`, `.ply` and `.gpkg`
   bytes equal those of master `98fe80c3` for the same command.
2. **One number:** with it, `area_m2` per triangle equals
   `landcover.plan_areas` of the file's own vertices and triangles, bit for
   bit, and equals the `.gpkg` `triangles.area_m2` of the same command.
3. **Covers the mesh:** the triangle entries sum to the area of the union
   of the written triangles (to 1e-9 relative; the rounding of a sum of
   10⁶ positive doubles is under 10⁶ · 2⁻⁵³ ≈ 1.1e-10, checked at 464,290).
4. **Every line is 0**, and the array has lines + triangles entries.

## 6. Tests (`@tester`)

Fixtures from `geotiff_fixtures.micro_tiff` and `pieces_fixtures`; figures
are hand-worked or from section 2.

1. **Hand-worked areas.** A 3 × 3-node DEM, x spacing 10 m, y spacing
   20 m, z = 3x (sloped, so a 3D area differs), `--stride 1 --plan-area`:
   8 triangles, each `area_m2` exactly 100.0 (½ · 10 · 20); in `.vtk` and
   `.ply`, binary and ASCII.
2. **Off-node domain.** `pieces_fixtures` `hills()`, `SQUARE`, tolerance 1,
   `--plan-area`: 1,211 triangles; sum within 1e-9 relative of 115,382.4 m²;
   each entry equals `plan_areas` of the file bit for bit.
3. **Off is byte-identical.** The commands of tests 1 and 2 without the flag:
   sha256 of the `.vtk` and `.ply` equal digests the tester records at
   `98fe80c3` (constants in the test, with that commit named in a comment
   of the test file); the `.gpkg` `triangles` table rows equal with and
   without the flag.
4. **Lines get 0.** Test 2's domain with `--lines` (or a hole), so the
   `.vtk` has lines: the raw array has lines + triangles entries, the first
   `lines` exactly 0.0, and `read_mesh` returns the triangle part only.
5. **Across formats.** One command written as `.vtk`, `.ply` and `.gpkg`:
   the three `area_m2` columns are equal bit for bit.
6. **Cut runs.** Test 2 with `--pieces 4`: 4 piece files, 1,163 triangles,
   each piece's `area_m2` equals `plan_areas` of that piece bit for bit; the
   sum over pieces equals the area of the union of their triangles to 1e-9
   relative (not the domain's: section 2).
7. **Reader refusals.** `read_mesh` on a hand-made file whose `area_m2` is
   `int`, or has one entry too few: `ValueError` naming `area_m2`; a file
   without it gives `None`.
8. **CLI.** `--plan-area` appears in `rasputin mesh --help` and in the MCP
   tool's argument model, which mirrors the click parameters.

**Mutants** the suite must kill: the 3D cross product's length (test 1,
z = 3x); a missing ½ (test 1); `dx²` for `dx · dy` (test 1, 10 by 20);
float32 storage (tests 2, 5); lines not zero-padded or the offset by
`lines` dropped (test 4); the flag ignored in `_write_piece` (test 6); the
array written when off (test 3). Not invariant-critical: no mutation run.

## 7. Acceptance (timings reported, not gated)

On mains power, median of 3, against `98fe80c3`: the fixture command of
section 2 with and without `--plan-area` (`.vtk` and `.ply`), and the same
tile cut by `--pieces 4`. Prediction: +0.03 s for the areas plus the
encode and write of 3.7 MB, within run-to-run noise of the 1.0-1.2 s run.
At the 637,000 km² basin's bound (every one of 7.1e8 DEM nodes a vertex,
1.4e9 triangles) the array is 11 GB over all pieces and about 100 s of
`plan_areas`, spread over pieces each inside its memory budget; a
tolerance mesh has far fewer triangles.

## 8. Size

Estimated +70 counted lines: `cli.py` +20 (flag, two call sites),
`io/vtk_legacy.py` +12, `io/ply.py` +15, `io/vtk_read.py` +18,
`write_triangles` reusing `plan_areas` −5. One PR.

## 9. Questions for Ola

1. Build the DEM pixel walk of section 4.2 (count, mean and spread of the
   DEM about each triangle) now, or only plan area? Default: **plan area
   only**; 4.2 waits for a reader that uses it.
