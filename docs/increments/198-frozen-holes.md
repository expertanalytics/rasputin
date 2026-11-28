# Increment 198: `--frozen-holes`, a domain's holes kept exactly as given

Status: **designed** (`@architect`, 2026-11-28, on master `949bc7d66`; amended after design
review: a hole vertex with no DEM under it is refused, section 3.4). Next: design review
(`@reviewer`), then red (`@tester`), green (`@developer`), test run, `@perf`'s acceptance run,
code review.

Lane: **full**. The change is in `rasputin mesh`'s refinement driver (`cli.py`'s `_refined` and
`_dem_mesh`), which `docs/increments/README.md` (*Acceptance*) counts as "what drives"
`include/terrain/refinement/`, so `@perf`'s acceptance run is part of the lane too (section 8).

## 0. The plain picture

Some users mesh a band around a mesh they already have, so that the two can be joined along
the hole boundary into one mesh. The band's domain is a polygon with a hole, and the hole is
the outline of the existing mesh. For the join to work, the band's mesh must keep the hole
exactly as given: every vertex on the hole ring stays, at its place, and no new vertex goes on
a hole edge. Today refinement splits a hole edge wherever the DEM bends under it, just as it
does every other constraint edge, so the band's hole ring no longer matches the outline of the
existing mesh.

`--frozen-holes` (off by default) asks for this: no vertex is added on a hole ring of
`--domain` and no edge of one is split. The meshing core already supports this through its
frozen edges, built for cut runs' seams (increment 23). This increment hands the hole edges to
that support. It changes no C++.

What it gives up: the DEM's error along a hole edge is not refined away. It is measured and
recorded (`frozen_edges_max_error_m`), not judged against the tolerance. Every DEM node off the
hole edges is still within the tolerance.

Two cases are refused, not repaired. One is a constraint line (a `--lines` line or a land-cover border)
that meets a hole edge between its vertices. The noder would put a vertex there, which splits
the hole edge. The run stops and names the point. The user either gives the line a vertex
there or meshes without the line. The other is a hole vertex with no DEM under it (nodata or
off the DEM): trimming would drop it and its edges, so the hole ring could not be joined.

## 1. Prior art: legacy and literature

### 1.1 Literature

- **Triangle's `-Y` switch** (Shewchuk, *Triangle: Engineering a 2D Quality Mesh Generator and
  Delaunay Triangulator*, 1996; switch pages at <https://www.cs.cmu.edu/~quake/triangle.switch.html>
  and <https://www.cs.cmu.edu/afs/cs/project/quake/public/www/triangle.Y.html>). `-Y`
  "prohibits the insertion of Steiner points on the mesh boundary", `-YY` "on any segment". The
  stated use is ours: keeping a boundary "so that it conforms to some adjacent mesh". Its stated
  cost is ours too: "you will probably sacrifice some of the quality of the mesh", with
  triangles of poor aspect ratio possible next to the kept boundary. TetGen's `-Y` is the same
  rule in 3D (<https://wias-berlin.de/software/tetgen/switches.Y.html>).
- **What this departs from:** Triangle's refinement is quality-driven; ours is error-driven
  against a DEM. Freezing a hole edge drops the guarantee that DEM nodes *on* that edge are
  within the tolerance. That guarantee is replaced by a recorded figure (section 3.6). Nodes
  *off* the frozen edges keep the guarantee, because refinement can still insert them (probed
  in section 2).
- **Novelty: none claimed.** Freezing boundary segments to join meshes is standard practice
  (Triangle `-Y`). The core's frozen edges already exist (increment 23, the seam protocol).

### 1.2 Legacy

`git grep -n -i -E "frozen|fixed.?edge|no.?steiner|stitch|merge.?mesh" legacy-archive -- legacy`
returns nothing. `git grep -n -i hole legacy-archive -- legacy` finds the CGAL polygon-with-holes
plumbing (`legacy/rasputin/triangulate_dem.h:83-139`, `legacy/rasputin/geometry.py:268-273`).
Holes there are subtracted from the domain and refined like any other boundary. No frozen
boundary exists, so nothing is carried across.

### 1.3 Ideas

None of the ideas-file entries is taken up or overlaps this one.

## 2. Measured (probes on master `949bc7d66`)

The probes are scratch scripts (`probe.py`, `probe2.py`, `probe3.py`, `probe3b.py` in the
scratch folder `arch-frozen-holes-20261128/`; not committed). All of them use the fixture of
section 5.

- **Today the hole edges are split.** With no flag, the output has **28 vertices on the hole
  edges**: the 10 hole vertices and 18 added by refinement. The output has 5,438 triangles.
- **Noding.** The noded hole chain (seen by a spy on `cli.node`) has **10 nodes** with no line,
  and 10 with the line x = 500.0004 m, which meets the hole only at its vertices. With the line
  x = 512.3 m it has **12 nodes**. The first added node in the chain's walk is **(512.3, 600.0)**;
  the other is (512.3, 400.0). The largest move of a hole vertex by the 1 mm snap is
  **0.0004 m**.
- **The flag, stood in for by monkeypatching** (`probe3.py`): bit 11 on the hole chains, the same
  bit as `frozen` to `_refined`, and the bit removed before `trim`. The output has **10 vertices
  on the hole edges** and 5,462 triangles. The record's `max_error_m` is 0.0999999829. There are
  **9,760** DEM nodes off the hole edges, and every one is within 0.1 m of the written mesh
  (largest 0.0999999829). On the **80** nodes on the hole edges the largest error is
  **0.6000000238**. The same 10 hole vertices result with `--start-min-angle 25`, with
  `--no-constraint-feet`, and with `--triangles 3000` (the budget path: 2,232 triangles, largest
  error off the frozen edges 0.2000122).
- **A constraint line inside a hole.** Lines are clipped to
  the domain polygon, holes included (`feature_input.py:349-358`,
  `shapely.intersection(line, dom)`). A line wholly inside the hole is dropped: the run prints
  "edges: 0 kept …, 1 outside the domain", and its mesh has 5,438 triangles, the same as with
  no line. A line crossing the hole is cut at the hole ring. Its pieces end on the hole edge,
  and no vertex lies inside the hole (0 in the x = 512.3 m run). So the mesher never meshes a
  hole's inside. Where a cut piece ends between two hole vertices, the noder adds that end as
  a node on the hole edge. That node is what section 3.3 refuses.
- **No DEM under a hole vertex** (`probe4.py`, `probe5.py`, with the stand-in of `probe3.py`).
  The fixture's DEM with nodata (-9999) at the nodes (400, 400) and (400, 600), two hole
  corners: the run exits 0, and of the **10** hole edges (bit 11 in `final.masks`) **7** survive
  `trim`; the 3 edges at the two corners are dropped silently. Both corners have no z
  (`final.valid` false). The first such end among the bit-11 edges in `final.edges` order is
  (400, 400) on the `--tolerance 0.1` run and (400, 600) on the `--triangles 3000` run's kept
  trial, which also keeps 7 of 10.
- **No sampled path with a domain.** `rasputin mesh --domain` without `--tolerance` or
  `--triangles` exits 2 with "--domain needs --tolerance or --triangles". So the flag, which
  needs `--domain`, always runs through `_refined`.
- **The core's frozen error does not cover DEM nodes.** On the runs above the edge strip's
  outcome has `on_frozen` 0 and `on_frozen_max_error` 0.0: the core counts only stored points
  (strip and check points) on a frozen edge, and its DEM node scan skips them uncounted
  (`scan.hpp:36-44`).

## 3. Design (`src_python/tin_engine/cli.py`, `run_record.py`; Python only)

### 3.1 The flag

`rasputin mesh --frozen-holes`, a boolean, off by default. Help text: "Keep each hole ring of
--domain as given: no vertex added on it and no edge of it split, so the mesh can be joined to
another along the hole. The DEM's error along the hole edges is recorded, not refined."

Refused, as usage errors and before any DEM is read where the inputs allow:
- without `--domain`: "--frozen-holes needs --domain";
- with `--simplify-domain`: "--frozen-holes keeps the holes as given; --simplify-domain would
  move them";
- on a cut run (`--pieces`, or a `--piece-budget` the window exceeds): "--frozen-holes does not
  run on pieces yet". The reason: the partition lines are clipped to the domain, so one that
  crosses a hole would put a node on the hole edge (section 3.3), and the user cannot give an
  automatic line a vertex. Joining frozen holes and seams is left out (section 9).

### 3.2 The bit

`FROZEN_HOLE = 1 << 11`, a module constant in `cli.py` (bit 9 is the cut's seam,
`src_python/tin_engine/features.py:180-182`; bit 10 is taken by `cli.py:1083`). With the flag, the hole
chains from `start_chains` (mask 0 today, `chains.py:51-52`) get this bit in `_dem_mesh`. The
noder unions masks (`_chain_masks`), so a hole edge that is also a land-cover border keeps
both. `_Tuning` gets a field `frozen: int = 0`, and `_refined` passes `frozen | tuning.frozen`
as its frozen mask. Every refinement call takes it this way: plain, budget trials, and the
edge strip and final check (they already take `frozen`, `cli.py:3717-3758`). The core then
"An edge whose mask meets frozen_mask gets no vertex: the scan skips the nodes on it, no foot
is taken on it, and the quality pass skips a node on it" (`refine.hpp:45-47`).

The bit is internal. `_refined` passes the masks to `trim` with the bit, counts it (section
3.4), and then removes it from `trim`'s output masks, so the written vocabulary, its fingerprint and
every file's edge columns are as today. A mask is cast as `np.uint32`, never as a negative
Python int. Without the flag, none of this code runs.

### 3.3 The check after noding

Right after `_engine`, with the flag, for each hole ring k (noded chain 1 + k; the noder keeps
chain order): let W be the noded chain's node walk (`indices_of`), and I the nodes of the
ring's input vertices (`node_of_input_vertex`). The run is accepted when `len(W)` equals the
ring's vertex count and every node of W is in I. Otherwise it is refused, writing no file:
- with a node of W not in I (the first in walk order), in the DEM's CRS: "a line meets the edge
  of hole {k + 1} of --domain between its vertices at ({x:.3f}, {y:.3f}); give the line a vertex
  there or mesh without it";
- with fewer nodes than vertices: "hole {k + 1} of --domain has {n} vertices but {m} after
  snapping to the {spacing} m grid; vertices closer than that merged".

### 3.4 The check after trim

In `_refined`, with the flag: let n be the number of edges in `final.edges` whose mask carries
the bit (the hole edges, as the core splits none) and m the number in `trim`'s output masks.
The run is accepted when m equals n. Otherwise it is refused, writing no file: "{n - m} of the
{n} hole edges of --domain have an end with no DEM under it (nodata or off the DEM), so the
mesh could not be joined along the hole; give holes that lie on DEM data; first at ({x:.3f},
{y:.3f})". The point is the first end with `final.valid` false among the bit's edges in
`final.edges` order; where none is invalid (a valid vertex left in no valid triangle), it is
the first dropped edge's first end. Every refinement call passes here, so a budget run refuses
on its first trial.

### 3.5 The contract

With `--frozen-holes`, every hole edge of the input appears in the output as one edge, with no
vertex on it but its two ends. Each end lies within the noder's snap of the input vertex: at
most half a lattice diagonal, `spacing * sqrt(2) / 2`, which is 0.71 mm at the default
`--snap-spacing` of 1e-3 m. Where a joining tool needs bit-equal coordinates, it matches within
that distance and takes the other mesh's (x, y).

### 3.6 Recorded

| entry | wording | where | value |
|---|---|---|---|
| `frozen_holes` | "Holes of --domain kept as given (--frozen-holes)" | mesh file, after `max_error_m` | the number of holes; None (not written) without the flag or with no hole |
| `frozen_edges_max_error_m` | "Largest height error along the kept hole edges, not refined" | mesh file, after `frozen_holes` | see below; None (not written) without the flag, or with no valid point |
| `hole_vertices_snapped_max_m` | "Farthest a hole vertex was moved by the snap to the grid" | counts section | max over hole vertices of the larger of \|x − x_in\| and \|y − y_in\|; None (not written) without the flag |

`frozen_edges_max_error_m` is computed in `_refined` right after `refine`, while the tile is
held. Points are placed along each frozen edge from end a to end b at steps of the tile's cell
size (the smaller of dx and dy), both ends included. Each point's DEM z comes from `sample`
(bilinear, as the mesh's vertex z). Each point's mesh z is the chord `z_a + t (z_b − z_a)`,
using refine's z of the two ends. The vertices on a frozen edge never move and none is added,
so the chord is the written mesh's z there. The value is the largest |DEM − chord| over valid
points. On the reprojected path the tile is the resampled grid, so the figure is measured
against that grid. These two entries go in the file because the file's `max_error_m` no longer
covers the hole edges; a reader of the file needs both. The core's `on_frozen_max_error` was
weighed as a simpler source and does not serve: it counts stored points only, and is 0 on the
tile path with 80 DEM nodes on the hole edges (section 2).

`max_error_m` and the stderr summary stay as they are. The scan skips nodes on a frozen edge
("never the argmax … not counted", `scan.hpp:36-44`), so `max_error_m` is the largest error
off the frozen edges.

### 3.7 Interfaces the tests call

Only the CLI: `rasputin mesh … --frozen-holes`, its exit code, its `.vtk` (`feature_names`,
points, triangles, lines and FieldData) and its `--record` JSON. No new public function.

## 4. Invariants

- **I1 (kept):** with the flag, the output's vertices on hole edges are exactly the hole
  rings' input vertices, one each, each within `spacing * sqrt(2) / 2` of its input. Every
  input hole edge is one output edge. Where that cannot hold (a crossing line, merged vertices,
  a hole vertex with no DEM), the run is refused.
- **I2 (off is identity):** without the flag no new code runs. Output bytes equal master's.
- **I3 (tolerance off the frozen edges):** every DEM node in the mesh that lies farther than the
  coincidence radius from a hole edge is within the tolerance. The error on the hole edges is
  recorded, not judged.
- **I4 (refusal):** a refusal writes no file. The option and crossing refusals (sections 3.1,
  3.3) come before refinement, the DEM refusal (section 3.4) after it.

## 5. Tests (`@tester`, red first)

**The fixture.** A DEM on a 10 m grid, 101 × 101 nodes over 0 to 1,000 m in x and y (local
coordinates; the files carry them offset by the lower-left corner, for example `micro_tiff`'s
tie point with y counted from 1,000 m below it). Its heights are z = (x − 500)² / 1,000 m,
exact at the nodes up to float32. The domain is the square 0 to 1,000 m with one hole, the
square 400 to 600 m. The hole's bottom and top edges carry the vertices x = 450.0004, 500.0004
and 550.0004 m, which are 0.4 mm off the 1 mm lattice. So the hole has 10 vertices: 4 corners
and 3 on each of two edges. Tolerance 0.1 m, no lines unless a test says so, defaults otherwise.

**Figures by hand.** The hole's bottom and top edges lie on grid rows y = 400 and 600, so DEM
nodes lie on them every 10 m. The side edges x = 400 and 600 have constant z = 10, so their
error is 0. On a 50 m piece of the bottom or top edge, for example from 450 to 500, the chord
runs from z = 2.5 to 0, and its error against the parabola is (x − 450)(500 − x) / 1,000. At
the nodes this is 0.4, 0.6, 0.6, 0.4 for x = 460, 470, 480, 490. The DEM is bilinear between
nodes, so the error is linear between nodes and largest at a node: **0.600**. The other three
pieces of each edge mirror this. Node counts: 101² = 10,201 nodes in all, 19² = 361 strictly
inside the hole, and 80 on the hole ring (the perimeter of a 21 × 21 block), which leaves
**9,760** off the frozen edges. All of these were probed equal (section 2).

- **F1 frozen.** With `--frozen-holes`:
  - Exactly 10 output vertices lie on the hole ring. Each is within 0.71 mm of its input vertex,
    and the three per edge at x = …0004 have moved 0.4 mm in x.
  - Each of the 10 input hole edges is one output line cell.
  - `hole_vertices_snapped_max_m` is 0.0004 to 1e-9, `frozen_holes` is 1, and
    `frozen_edges_max_error_m` is 0.600 to 1e-6, at x = 420, 430, 470, 480, 520, 530, 570
    and 580 on both edges. Measured: 0.6000000238.
  - `feature_names` is as without the flag (no bit 11 in any written mask).
  - With a `--lines wall=` line x = 500.0004 m from y = 0 to 1,000 m (with `--edges wall`),
    which meets the hole only at its vertices: accepted, and the same 10 hole vertices.
  - The same 10 hole vertices with `--triangles 3000` in place of `--tolerance`.
- **F2 within tolerance off the frozen edges, largest error on them recorded.** In the F1 run,
  each of the 9,760 DEM nodes off the hole ring is within 0.1 m (plus 1e-9) of the plane of the
  output triangle holding it, recomputed from the written vertices (measured: largest
  0.0999999829). The 80 nodes on the ring are not judged. Their largest error is F1's recorded
  figure.
- **F3 refusal.**
  - A `--lines wall=` line x = 512.3 m from y = 0 to 1,000 m meets the bottom and top hole edges
    between their vertices. With the flag: exit code 2, no output file, and a message naming
    "the edge of hole 1" and the point local (512.3, 600) or (512.3, 400), printed in the files' coordinates (the offset added) to 3 decimals.
    The walk order is not a contract; measured today, 600 comes first. Without the flag the
    same run is accepted.
  - `--frozen-holes` without `--domain` is refused with "needs --domain".
  - `--frozen-holes` with `--simplify-domain` is refused with "would move them" (exit 2 alone
    proves nothing: master exits 2 on the unknown option).
  - The fixture's DEM with nodata (-9999, `micro_tiff(..., nodata="-9999")`) at the nodes
    (400, 400) and (400, 600), two hole corners, with `--tolerance 0.1` and again with
    `--triangles 3000`: exit code 2, no output file, and a message holding "3 of the 10 hole
    edges", "no DEM" and the point local (400, 400) or (400, 600) in the files' coordinates to 3
    decimals (measured: (400, 400) first on the tolerance run, (400, 600) on the budget run; the
    order is not a contract). Without the flag both runs are accepted.
  - `--frozen-holes` with `--pieces 4` is refused with "does not run on pieces yet".
  - Each refusal writes no file.
- **F4 off is identity.**
  - Without the flag, the fixture's hole edges are split as today: more than 10 output vertices
    on the hole ring (measured on master: 28), and the record has none of the three entries.
  - A domain without a hole (the 0 to 1,000 m square) gives byte-equal `.vtk` files and equal
    `--record` JSON except its `command` entry, with and without the flag.
  - Byte equality with master's own output is `@perf`'s check (section 8). A suite cannot run
    master, and pinning a hash would tie the test to one platform's floating point.

No suite here is invariant-critical, so no mutation testing.

## 6. Degeneracy policy

- **A line along a hole edge,** ending at hole vertices: accepted, because the hole chain keeps
  its count and the edge's mask is the union. A line ending between hole vertices is refused
  (section 3.3).
- **A land-cover border meeting a hole edge** is refused like a line, because borders are
  clipped to the domain the same way. Land-cover borders across a frozen hole are out of scope:
  the use meshes the band without land cover.
- **Two hole vertices within the snap** merge. Refused (section 3.3, second message).
- **A domain in another CRS:** "as given" means as moved into the DEM's CRS. The snap figure is
  measured against those coordinates. A user joining meshes gives the domain in the DEM's CRS.
- **A hole vertex with no DEM under it:** refused (section 3.4). **Nodata under the middle of
  a hole edge** whose ends have z: the edge is kept, and the points there are left out of
  `frozen_edges_max_error_m`.

## 7. Speed

Off by default, so no run changes unless it asks. `bench_quick` has 4 cases with a `--domain`
(`docs/benchmarks/quick/cases.toml`), and all 4 are polygons with **0 holes** (checked with shapely). So the gated cost there is 0.

With the flag, the new work is one membership pass over each hole chain after noding, and one
sampling pass along the hole edges: (ring length / cell) points. That is, for example, 13,000
points for 130 km of ring at 10 m, one bilinear sample each, against seconds to minutes of
refinement. Not gated.

## 8. Acceptance (`@perf`)

On AC power, the merge head against `949bc7d66` (or the then base). `bench_quick` and the 1 m
acceptance benchmark, flag off: every output byte-identical to the base, and times within the
benchmark's noise. This is I2 on real cases. No timed run with the flag.

## 9. Left out

- Frozen holes on cut runs (refused, section 3.1).
- Moving or snapping a crossing line onto a hole vertex. It is refused instead; a user who
  needs borders across a frozen hole puts their vertices on it.
- A frozen outer ring. Only holes are frozen; the use joins along a hole.
- The joining tool itself. This increment only makes the band's mesh joinable.

## 10. Size

Python only, no C++, no binding change. Estimate, in counted lines:

| where | what | lines |
|---|---|---|
| `cli.py` `mesh` | the option and its help | 9 |
| `cli.py` | the three refusals | 9 |
| `cli.py` | the flag to `_dem_mesh`, `_Tuning.frozen`, the partial | 6 |
| `cli.py` `_dem_mesh` | the bit on hole chains | 3 |
| `cli.py` | the check after noding, its two messages, the snap figure | 16 |
| `cli.py` `_refined` | the frozen mask from `_Tuning`; the chord sampling; the check after `trim`, its message, the bit removed | 26 |
| `run_record.py` | three entries | 5 |
| **total** | | **about 73** |

One pull request, far under 700.

## 11. Questions for Ola

None. Defaults taken, each open to a no:
- The flag is refused on cut runs, not supported.
- The frozen bit stays internal, so written files carry no new edge column.
- `frozen_holes` and `frozen_edges_max_error_m` go in the mesh file, because the file's
  `max_error_m` no longer covers the hole edges.
