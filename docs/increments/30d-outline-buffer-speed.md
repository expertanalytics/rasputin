# Increment 30d — growing the catchment outline, made faster

Status: **designed** by `@architect`, 2026-10-07, on branch
`worktree-buffer-speed-2` at `d9f9c025` (`@perf`'s baseline on master
`f81b20b7`). Every `file:line` below is pinned to `f81b20b7`. Ola ruled
earlier on 2026-10-07: grow the outline in pieces. One question is still open
(section 9); the design takes its default, and a "no" changes section 3.3
only.

**What this is.** `rasputin mesh` grows the catchment outline ("buffers" it:
the polygon widened by a distance) in four places. On a raster-traced outline
(SMHI's: every edge 10 m long and axis-parallel, a staircase) GEOS's buffer
is very slow once the distance spans two or more steps. On Lagan it is 98 %
of `decode` and 95–97 % of `features read`
(`docs/benchmarks/2026-10-07/30d-buffer/README.md`). The audit that found it
is `docs/increments/perf-audit.md` on branch `worktree-perf-audit`
(`1da4a144`, PR #208, not merged), findings F1a–F1c.

**What changes for a user.** Nothing in any output: the mesh stays
byte-identical. Lagan with the European CORINE GeoPackage goes from 103.4 s
to about 17 s; Norwegian runs (smooth NVE outlines) do not change.

## 1. Prior art: legacy and literature

### Literature

No novelty is claimed; each change uses a textbook identity.

- **The features region from the hull (3.2).** The convex hull of a Minkowski
  sum is the Minkowski sum of the convex hulls, so conv(P ⊕ B) = conv(P) ⊕ B
  for a disc B (de Berg, Cheong, van Kreveld and Overmars, *Computational
  Geometry: Algorithms and Applications*, 3rd ed., Springer 2008, ch. 13).
  Recalled, not reread. The two computed regions differ only in where GEOS
  places vertices on the arcs (8 segments per quarter circle).
- **Growing in pieces (3.3).** For a closed region P, P ⊕ B = P ∪ (∂P ⊕ B):
  every point within d of P is in P or within d of its boundary. The boundary
  is grown as overlapping open pieces and the pieces united with P by GEOS's
  unary union (cascaded union, JTS `CascadedPolygonUnion`, Martin Davis).
  GEOS's own buffer (JTS `BufferBuilder`) offsets every edge and nodes the
  result; on a staircase whose step is shorter than the distance, the offsets
  of neighbouring steps overlap and the noding grows faster than linearly.
  Short pieces keep each noding small. Recalled, not reread.
- **Departure.** The pieced region equals GEOS's only to rounding on the
  outlines measured, not by construction (section 4). That is why the route
  is gated (3.3).

### Legacy

`git grep -n "\.buffer(" legacy-archive -- legacy` returns
`legacy/rasputin/geometry.py:259`, `legacy/rasputin/gml_repository.py:204`,
`legacy/tests/test_mesh.py:129` and `legacy/tests/test_mesh.py:167`: plain
GEOS buffers, nothing about cost. Nothing is carried across.

## 2. What is slow (measured, `@perf`, AC, median of 3)

| call (`f81b20b7`) | what | Lagan | Ljungan | Numedalslågen |
|---|---|---|---|---|
| `target_grid.py:107` | domain grown by √2·h, mitred, for the grid's bounds | 17.3 s | 2.1 s | not on this path |
| `dem_input.py:213` | **the same buffer again**, for the needed region | 16.7 s | 2.1 s | not on this path |
| `dem_input.py:250` | projected path: domain grown by the cell diagonal, mitred | not on this path | not on this path | 0.005 s |
| `feature_input.py:161` | domain grown by 100 m, round, **only to take its convex hull**; called at `:294` (GeoPackage box), `:311` (region), `:316` (geographic source) | 28 s per call, 2 calls | 2.2 s per call | 0.017 s |

Candidates in full runs, Lagan with the GeoPackage: master 103.4 s; (1) one
buffer instead of two, 61.6 s; (1)+(2) the features region from the hull,
31.1 s; (1)+(2)+(3) the outline grown in pieces, 16.4–18.1 s. Every `.vtk`
byte-identical on Lagan and Ljungan, both feature inputs.

## 3. The blueprint

All four changes are Python, inside the I/O boundary (`CLAUDE.md` §2): the
grown polygons serve the target grid, the mosaic plan and the features
pre-clip, all Python. No C++ (see 3.5).

```
dem_input._reprojected                      dem_input._domain_plan
  grid, grown = target_grid_for(...)  (3.1)    grown = grow_mitred(P, diag)  (3.3)
     └─ grow_mitred(P, √2·h)          (3.3)
  source_region(grid, meta, grown)  (unchanged; its own buffer at :145 stays GEOS)

feature_input.source_region(domain, dem_crs, source_crs)          (3.2)
  conv(P).buffer(MARGIN) → segmentize → reproject → convex hull

grow.py (new):  uses_pieces(polygon, distance) -> bool            (3.3, the gate)
                grow_mitred(polygon, distance) -> Polygon
```

### 3.1 One buffer, not two (F1a)

`target_grid_for(domain, target, spacing)` returns
`tuple[TargetGrid, Polygon]`: the grid and the grown domain it was snapped
around. `dem_input.py:212-213` unpack it; line 213's buffer goes. The rule
"grown by √2·h with mitred corners" stays in one function. `target_grid_for`
has one production caller (`dem_input.py:212`) and two test calls
(`tests/python/test_target_grid.py:148`, `:160`), which `@tester` amends.

### 3.2 The features region from the hull (F1b)

`feature_input.source_region` (`feature_input.py:158-167`) keeps its
signature and R5's rule (16b, R5: the region holds the domain's image with
close to 100 m to spare; anything the margin lets in is clipped exactly in
R6). Only line 161 changes:
`shapely.convex_hull(domain.polygon).buffer(MARGIN)` instead of
`domain.polygon.buffer(MARGIN)`. The hull of the result is then the same
region to within the arc chord error, 100·(1 − cos(π/32)) = 0.48 m
(scale: MARGIN = 100 m; measured 0.215–0.252 m boundary distance on Lagan,
Ljungan, Numedalslågen). It is 4 ms instead of 28 s, so the three calls per
source need no caching.

### 3.3 The outline grown in pieces, behind a gate (F1c)

New module `src_python/tin_engine/grow.py`, pure functions, no state:

- `grow_mitred(polygon: Polygon, distance: float) -> Polygon`: the polygon
  grown by `distance` with mitred corners (shapely's default mitre limit,
  5.0). When `uses_pieces` is false it is exactly
  `polygon.buffer(distance, join_style="mitre")`. When true: the exterior
  ring's vertices cut into pieces of `PIECE_EDGES` edges, each piece
  overlapping the next by one edge, the last wrapping past the start vertex;
  every piece grown with `shapely.buffer` on the array (mitre joins, flat
  caps), then `shapely.union_all` of the polygon and the pieces. As in the
  audit's prototype (`perf-audit-probes/run_patched.py`, patch 3).
- `uses_pieces(polygon: Polygon, distance: float) -> bool`, the gate. True
  only when all three hold:
  1. the polygon has no holes (a hole's ring would need its own pieces; no
     measured outline has one);
  2. the exterior has at least `PIECES_FROM` vertices (closing vertex not
     counted);
  3. the median edge length is at most `distance / 2` — the distance spans
     at least two steps, where GEOS's time jumps (Lagan: 0.03 s at one step,
     10 m; 1.8 s at two, 20 m).
- Constants: `PIECE_EDGES = 1000` and `PIECES_FROM = 5000`. Scale: 10 m
  steps, distances 14–44 m. Largest input checked: Lagan, 56,261 vertices,
  d = 43.84 m (at 1,000 edges the region equals GEOS's, 0 m² symmetric
  difference, equal bounds; at 250 and 500 it differed by up to 6.1e-7 m).
  Below 5,000 vertices there are at most five pieces and GEOS is fast.
- Callers: `target_grid_for` (`target_grid.py:107`) and `_domain_plan`
  (`dem_input.py:250`). Not `target_grid.py:145` (the grown polygon moved to
  the DEM's CRS and grown by a source cell: 0.2 s on Lagan, already no longer
  a fine staircase) and not `target_grid.py:243` (hull of 512² blocks).

**Why the gate, and why not by vertex count alone.** On a smooth outline the
pieced region is not GEOS's: on Numedalslågen's NVE outline it reaches 10.8 m
beyond GEOS's at one place. A probe for this design
(scratch, not committed; shapely 2.1.2, GEOS 3.13.1) found the difference is
**not caused by the pieces**: GEOS's buffer of the whole ring as one line,
united with the polygon, gives the same 101.27 m² and 10.83 m. It sits at a
sharp notch (vertices 7537–7538 turn by −162° and +149°), and it changes with
the mitre limit (101 m² at 5, 164 m² at 10). So a polygon's buffer and its
ring's buffer differ at sharp turns. Both regions still cover the round
buffer to 1e-5 m (checked on Numedalslågen and Ljungan). A count-only gate
would catch Numedalslågen (14,093 vertices) and the São Francisco basin's BHO
outline (77,310 vertices, edges ~107 m, GEOS 0.03 s), where pieces gain
nothing and lose byte-identity. Rule 3 keeps both on GEOS: Numedalslågen's
median edge is 48.1 m against d = 14.14 m; BHO's ~107 m against ~42 m.
A smooth outline with short edges would pass the gate; its region would then
be a superset of GEOS's by up to the mitre excess, still covering every
point within d (section 4), so correct, though not byte-identical with
master.

### 3.4 What does not change

Refine and mesh code: **untouched** (`include/terrain/refinement/`,
`include/terrain/mesh/` and the binding are not in the diff). So `@perf`'s
`tools/bench.py` acceptance does not apply; the short check of section 7
does. `mosaic`, `TargetGrid`, R5's pre-clip and every model stay as they are.

### 3.5 Against the assessment framework

1. Data and execution: no configuration added; the constants are
   module-level, not user options.
2. State: `grow_mitred` and `uses_pieces` are pure, thread-safe.
3. Dependencies: shapely and NumPy only, both present. Nothing from §2's list.
4. Async: unchanged; decode stays a synchronous call run off the loop.

**No C++.** GEOS is already compiled code; the cost is its algorithm on a
staircase, not Python overhead (about 56 Python iterations build the pieces
on Lagan). A C++ offset would duplicate GEOS for regions only the Python I/O
layer reads. Ola's broader question (are there slow patterns, is Python
doing C++'s work) is answered by the audit on PR #208: no per-element Python
work sits in a hot path today; its other findings (F2–F9) are separate
increments.

## 4. What makes a buffer "the same"

- **Containment, always (both routes).** The region covers the polygon grown
  by `distance − 1e-5 m` with round joins (`quad_segs=64`). This is what
  the grid bounds and the needed region need (every node bilinear z reads
  lies within the cell diagonal). Scale: coordinates up to 7e6 m (UTM
  northings), d 14–100 m; checked on Ljungan (23,259 vertices) and
  Numedalslågen (14,092). GEOS's own buffer falls short of `d − 1e-6` by
  8e-7 m on Numedalslågen, so 1e-6 is too tight for either route.
- **Equality with GEOS on the gated route.** Boundary Hausdorff distance to
  GEOS's mitred buffer at most 1e-6 m, and every bound within 1e-6 m.
  Measured 0 at 1,000-edge pieces on Lagan and Ljungan. Bounds snap outward
  to the lattice, so a bound within 1e-6 m of a lattice line could still
  move the grid by one line; that is a residual risk (section 8), not an
  error.
- **The features region (3.2).** Covers the domain grown by
  `MARGIN − 0.5 m` (round); its bounds within 0.5 m of master's region.

## 5. Tests for `@tester` (the red suite)

Small synthetic fixtures; no new data files. Suggested home:
`tests/python/test_grow.py`, plus the amended `test_target_grid.py`.

1. **Gate table** (`uses_pieces`): a staircase of ≥ 5,000 vertices with
   10 m steps at d = 43.84 → True; the same at d = 14.14 (under two steps)
   → False; a staircase of 4,999 vertices → False; a 6,000-gon circle with
   edges longer than d/2 → False; a staircase with a hole → False.
2. **Equality** on the gated staircase (radius about 6.5 km, 10 m steps,
   ≥ 5,000 vertices): `grow_mitred` against
   `buffer(d, join_style="mitre")`, as in section 4.
3. **Containment** for both routes: the staircase, and a smooth outline
   with a sharp notch (a −162° turn followed by +149°, like Numedalslågen's
   vertices 7537–7538), section 4's first bullet.
4. **Off the gate it is GEOS's:** on the smooth notched outline,
   `grow_mitred`'s WKB equals `buffer(d, join_style="mitre")`'s.
5. **One buffer (3.1):** `target_grid_for` returns the grid and the grown
   polygon, and the grid's bounds are the grown polygon's snapped outward.
   A spy on `grow.grow_mitred` over one reprojected decode (the `velhas`
   fixture) counts one call.
6. **Features region (3.2):** on the staircase, section 4's third bullet,
   against master's formula kept in the test.
7. **Byte-identical meshes**, platform-independent (no stored hash): the
   same `mesh` run twice in one test, once as master behaves (the gate forced
   false by monkeypatching `grow.PIECES_FROM` above any count, and master's
   `source_region` formula patched in), once as built; the `.vtk` bytes
   equal. Cases: (a) the `velhas` DEM (reprojected path) with a staircase
   domain in EPSG:31983 above 5,000 vertices (5 m steps, radius about 3.5 km,
   inside the fixture's ~9.6 km); (b) the `dtm10` fixture with the CORINE
   GeoPackage fixture (`tests/fixtures/corine/clc2018_7908_3.gpkg`) for the
   features region. Lagan and Ljungan do not fit in the suites (100 s and
   15 s per run, data not in CI); their identity is section 7's.

Not invariant-critical: no mutation round.

## 6. Size

| file | counted lines, net |
|---|---|
| `grow.py` (new) | about +25 |
| `target_grid.py` | about +1 |
| `dem_input.py` | about −1 |
| `feature_input.py` | 0 |
| total | **about +25**, far under 700 |

Check: `python3 tools/count_loc.py <base> <head>`.

## 7. The speed check (`@perf`, at most 15 minutes)

One case: **Lagan with the European GeoPackage** (as in the baseline's
`scripts/run_mesh.sh`), on the branch head, AC power, one warm-up and three
runs. Report the median of `total`, `decode` and `features read`, against
the baseline's master medians (103.38, 35.12, 57.23 s; same machine, same
inputs, so master is not re-run). Expected about 17 s, 4–5 s and 1–2 s.
The `.vtk`'s SHA-256 must equal master's in
`docs/benchmarks/2026-10-07/30d-buffer/raw/runs/vtk_sha256.txt`. Plus one
run of Numedalslågen (UTM33 GeoPackage, 7 s) for its SHA-256 only: the gate
keeps it on GEOS. No sweep, no profile unless the check misses.

## 8. Risks

- A gated outline whose grown bound lies within 1e-6 m of a lattice line
  could gain or lose one grid line against master. The mesh would differ,
  still correct. Not seen on Lagan or Ljungan.
- A staircase reprojected into a rotated CRS keeps short edges, so rule 3
  still gates it; nothing depends on axis-parallel edges.
- Why a polygon's buffer and its ring's buffer differ at sharp turns was
  located, not explained; the gate makes it not matter for byte-identity.

## 9. Questions for Ola

1. **Pieces only on long staircase outlines?** Use pieces only when the
   outline has at least 5,000 vertices *and* its median edge is at most half
   the distance (a staircase where GEOS is slow), and GEOS's single buffer
   otherwise. *Default: yes* (this design). A "no" (pieces on every
   outline of 5,000 vertices or more) changes section 3.3 only, and test 4
   becomes a containment test: Numedalslågen's region would then be up to
   10.8 m larger than master's and its mesh would no longer be
   byte-identical.

## 10. ROADMAP

Row 30d in `ROADMAP.md`.
