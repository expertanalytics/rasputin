# Increment 27: a vertex on a DEM node reads that node

Status: **implemented, in code review**, not pushed. Designed and approved
2026-10-03 (Ola ruled: nodes only, the rule in the C++ sampler); red `6605dfe`,
pins on the red step's choices `80990e0`, green `6fd4076` and `a0817ba`; 13 net
production lines; no `@perf` run needed (see "`@perf`"). The branch has master,
with 25 (#157), merged in. Decisions for Ola are in "For Ola".

**The defect.** Without `--tolerance`, z comes from `_core.sample`, which is
`terrain::raster::bilinear` over every mesh vertex
(`include/terrain/raster/sample.hpp:26-59`). `bilinear` refuses a point when
any of the four corners of its cell is NoData, also a corner whose weight is
zero (`sample.hpp:42-44`). A vertex exactly on a DEM node that has a value
then gets no height when a neighbour has none, and the trim
(`src_python/tin_engine/elevation.py:52-55`) removes it with its triangles.
So a ring of valid data, one cell wide, is lost along every NoData area. On
increment 25's fixture at stride 1 that is 36 vertices removed where 25 are
NoData (`25-plain-output.md`, "NoData on the no-tolerance path").

Increment 12 accepted this as a stopgap ("accepted for now ... changing the
sampler is not this increment's job", `12-dem-to-mesh.md`, R2). Increment 25
pinned today's behaviour in a test, and said a later C++ increment would
change it back. This is that increment.

## Prior art: legacy and literature

### Literature

The method is ordinary bilinear interpolation on a regular grid (Press,
Teukolsky, Vetterling and Flannery, *Numerical Recipes*, 3rd ed., 2007, §3.6,
"Interpolation on a grid in multidimensions"). At a grid node the bilinear
weights are 1 for that node and 0 for the other three corners, so the
interpolant equals the node value there and the other corners do not enter.
The method says nothing about missing values. That is a convention each tool
picks.

GDAL is a prohibited dependency (`CLAUDE.md` §2), but its convention is prior
art: its warper leaves NoData source pixels out of its resampling kernels, and
ticket #3658 records that its cubic-spline and Lanczos kernels renormalise
the weights when NoData pixels lower the total weight below 1
(<https://trac.osgeo.org/gdal/ticket/3658>). I did not read GDAL's bilinear
kernel source and make no claim about it.

**What differs here, and why.** This increment does not renormalise. It skips
only corners whose weight is exactly zero, and that happens only when the
point is a node. The value returned is then the bilinear value in exact
arithmetic, so no new surface is invented over missing data. Renormalising
would give a height wherever *any* corner has data. That is a larger change
in what a NoData area means for the mesh, and nobody has asked for it.

**Novelty:** none is claimed. Searched: "GDAL warp kernel bilinear nodata
source pixels excluded weights renormalized" and "scipy RegularGridInterpolator
nan neighbour zero weight". Found: the GDAL ticket above. Nothing found changes
the design.

### Legacy

The legacy sampler had no NoData handling at all (`raster.hpp:58-59` says so).
Grep against the archive tag:

```
$ git grep -l -i -E 'nodata|no_data|bilinear|interpolat' legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/triangulate_dem.h
```

At the tag, lines 375-388 of `legacy/rasputin/triangulate_dem.h` are
`get_interpolated_value_at_point`, bilinear with no sentinel test, and line 182
of `legacy/bindings.cpp` binds it; `globcov_repository.py` uses
`no_data` only as a land-cover class, with class code 230. Nothing is carried across:
there is no legacy NoData rule to keep.

## The rule

**A point that is a DEM node, bit for bit, reads that node and nothing else.**

In floating point, "zero weight" is defined by node identity, not by the
computed weights:

- Let `fc = clamp((p.x - x_min) / delta_x, 0, cols - 1)` and
  `fr = clamp((y_max - p.y) / delta_y, 0, rows - 1)`, and `n` the node
  `(round(fr), round(fc))`. The point is **on a node** when
  `RasterGeometry::node(n) == p`, both coordinates bit for bit.
- On a node, the other three corners of the cell have weight zero and are not
  read. z is `value_at(n)`, refused only when `n` itself is NoData or NaN.
- Every other point reads all four corners of its cell, exactly as today:
  refused if any of the four is NoData, and the same expression, so the same
  bits, as today.

That "every other point" includes:

- **A point one ulp from a node.** Its far-corner weights are about 1e-16,
  not zero. Strictly, those corners contribute, so a NoData there refuses the
  point. Exactness is the rule; there is no tolerance band.
- **A point on a cell side** (one coordinate on a node line, the other not).
  Two corners have weight zero, but this increment does not skip them; see
  "Why cell sides are not in this increment".

Why this test, and not `tx == 0.0` on the computed weights: `cell_of` can put
a node into the cell before it. When `(p.x - x_min) / delta_x` rounds just
below the integer, `tx` comes out as 1 - 1e-16 instead of 0, and the zero
weight is never seen. Measured with a test program that includes today's
`sample.hpp` (`c++ -std=c++20 -O2`, Apple clang 21) and calls `bilinear` at
`g.node(c)` for every node of a 200 x 300 grid of values drawn uniformly from [0, 1000) (`std::mt19937`, seed 1, so
the counts depend on the draws and the standard library):
with origin 500000.3 / 7900000.7 and spacing 0.7 / 0.3, 38 580 of 60 000
samples are not bit-equal to the node value (largest difference about
6.6e-7 m); with origin 0.1 / 100.1 and spacing 0.1, about 31 000. `@reviewer`'s
rerun with other draws got 31 116 and 6.72e-7. With origin 0 and spacing 1,
with a dyadic grid (500000.5 / 7900000.25, spacing 0.5 / 0.25), and with the
Kartverket fixture's geometry, there are none. Node identity is also the relation the producers use:
`subsample` builds stride vertices with `node`'s expression
(`grid_domain.py:68-69`), and refine decides that a start vertex is a node by
the same `node(round) == p` test (`lattice_position`, `refine.hpp:132-141`).

### Where it lives

- **`RasterGeometry::node_at(const Point2&) -> std::optional<CellIndex>`**, new,
  in `include/terrain/raster/geometry.hpp`: the node `p` is, bit for bit, or
  nullopt. nullopt also for a non-finite `p` and for a point outside the node
  rectangle (`cell_of`'s test). The clamp-and-round is written the way
  `lattice_position` writes it, so `node_at` says "node" exactly where
  refine's whole test at `refine.hpp:417-426` does; a test pins that (S5).
  `node_at` answers on any raster, a 1 x N one included; it is `bilinear`'s
  guard, not `node_at`, that keeps such a raster at nullopt.
- **`bilinear`** (`sample.hpp`): after the existing `bilinear_cell_of` guard,
  `if (const auto n = g.node_at(p)) return is_nodata(*n) ? nullopt :
  value_at(*n);`, then the existing four-corner code, unchanged. The guard
  stays first, so a raster with fewer than 2 rows or columns still answers
  nullopt everywhere, nodes included, as today.
- `bilinear_batch`, the `sample` binding and the Python side change in no
  code. The binding's docstring (`bindings/core.cpp:935-940`), the stub's
  (`src_python/tin_engine/_core.pyi`, `sample`) and the comment above
  `bilinear` say the new rule.
- `lattice_position` is **not** changed to call `node_at` here. That would
  touch `refine.hpp` for a refactor with no behaviour in it. It is a one-line
  follow-up once this lands, if wanted.

### Why cell sides are not in this increment

1. **No vertex on the path being fixed lies on a cell side.** Without
   `--tolerance` there is no `--domain` and no `--features` (`src_python/tin_engine/cli.py@44fa7f5:845`
   refuses `--domain` without `--tolerance`). The vertices are the stride
   nodes and the ring through them. The noder makes no crossings there, and
   the triangulation adds no points. Every vertex is a node.
2. **On the tolerance path the two samplers must agree.** refine measures and
   carves with `vertex_z` (`include/terrain/refinement/scan.hpp:75-99`) and
   writes the output z of an off-node start vertex with `raster::bilinear`
   (`refine.hpp:425-426`). If only `bilinear` skipped cell-side corners, a domain
   vertex on a lattice line next to NoData would be void to the scan but valid
   in the output. A cell-side rule must change both. That changes refine:
   carving, the golden digests, and `@perf` acceptance. Increment 23 cuts
   pieces along lattice lines, so at basin scale such vertices are common.
3. **The frames disagree on what "on a side" means.** In the world frame it is
   `p.x == node x`; in the fractional frame it is "col is an integer". Increment
   16 records that the two disagree within an ulp of a cell line
   (`16-domain-polygon.md`, under T-real). A cell-side rule has to choose a
   frame, which is design work of its own.

If Ola wants it, it is a separate increment touching `scan.hpp` and
`sample.hpp` together, with `@perf`. Recommended: propose it only if a
measured run shows ragged seams along NoData.

## Which paths change

| path | changes? | why |
|---|---|---|
| no `--tolerance`, stride grid (12's R6) | **yes** | every vertex is a node; a valid node next to NoData keeps its height and its triangles |
| no `--tolerance`, mosaic tile (15) | **yes**, the same way | same `sample` call on the assembled tile |
| `--tolerance`, start-boundary / domain output z (`refine.hpp:425-426`) | **no**, by construction | `bilinear` is called there only for a start vertex that is **not** a node by refine's whole test: `lattice_position`, then `node && g.node(c) == p` (`!given \|\| (node && g.node(c) == p) ? vertex_z : bilinear`). `node_at` answers the same as that whole test (S5), so the new branch is never taken from refine |
| `--tolerance`, `vertex_z` (scan, carving, feet, edge strip) | **no** | not touched; it already reads a node with `value_at` (`scan.hpp:86-87`) |
| `--tolerance`, reprojected (15c: `resample`, `refine_points`) | **no** | `resample` is Python with its own four-corner rule (`target_grid.py:163-196`); `refine_points` calls `vertex_z`, not `bilinear` |
| a reprojected tile meshed without `--tolerance`, if a run does so | the stride sampling of the resampled tile follows the new rule; `resample` does not change | the target grid's nodes are `col0 * h` with integer `h` (`target_grid.py:48`), so they are exact |

**What changes in the stride output, exactly.** (a) Valid nodes next to NoData
keep their z and their triangles; `nodata_vertices_removed` becomes the number
of NoData vertices. (b) A node vertex whose old sample was a few ulps off
(the measurement above) now gets the node value exactly. Where neither applies,
the output is byte-identical. On the Kartverket fixture
(`tests/fixtures/dem_archive/7908_3_10m_z33.tif`, 5051 x 5051, origin
799750 / 7950250, 10 m) both are empty. Sampling all 25 512 601 nodes with
today's extension: 42 476 refused, all of them NoData nodes, so no valid node
is refused; and every valid sample is bit-equal to the node value. Its stride
output should therefore be byte-identical before and after (test S8). The
measurement used the extension installed in the main checkout's venv, built
2026-09-29. `sample.hpp` has not changed since `3e2e9ec`, and `geometry.hpp`
has gained only `operator==` (`67df94c`), so that build samples as master
does.

### A limit of the exact rule

The noder writes vertices at `SnapGrid::world(snap(p))`, multiples of 1 mm
(`include/terrain/noding/node.hpp:355-358`, `DEFAULT_SNAP_SPACING = 1e-3`).
Where a node's coordinate is not already the double nearest a multiple of 1 mm,
the stride vertex moves off the node by up to 0.5 mm, is not a node any more,
and keeps today's four-corner rule. A NumPy emulation of the snap on the
Kartverket fixture's stride-1 nodes: 0 of 25 512 601 move. On made-up grids:
origin 500000.3 / spacing 0.7, 36 000 of 90 000 move; origin 1234.5678 /
spacing 30, all of them; origin 0 / spacing 30.87, 9 911.

A second source is in the C++ itself. Apple clang contracts
`x_min_ + col * delta_x_` in `RasterGeometry::node` into a fused multiply-add,
and NumPy (`subsample`) does not. With origin 0.1 and spacing 0.1, 95 of 300
`node` x coordinates differ from NumPy's (a test program built with
`c++ -std=c++20 -O2`, compared with Python's `math.fma` and with NumPy). On
such grids the stride vertex is not `node(n)` bit for bit, even before the
noder.

Real inputs are not like this. Kartverket and 1 m grids have integer origins
and spacings, and the reprojected target grid is `col0 * h` with integer `h`.
Both effects need a non-integer, non-dyadic origin or spacing. The
`subsample` docstring says that using `node`'s expression keeps the bits equal.
Under contraction that is not true, so this increment corrects the docstring
(a documentation defect found here, fixed in this PR). Test fixtures for this
increment use integer or dyadic geometry; see S1.

The same FMA effect means that `lattice_position` classifies the stride start
vertices of such a grid as off-node on the tolerance path. That was true before
this increment. It is reported to Ola under "For Ola" and not changed here.

## Tests for `@tester`

Red first, on this branch (master, with 25, merged in). No mutation round:
the change is one early return, and the golden digests are the backstop for
the one invariant that matters (S5, S7). The C++ cases go in
`tests/cpp/unit/test_raster_view.cpp` (batch and view) and
`tests/cpp/unit/test_raster.cpp` (geometry and `bilinear`), except S5, which
has its own target, `tests/cpp/unit/test_raster_node_at.cpp`, so that its
compile failure before green does not hide S1-S4 and S6.

New:

- **S1. On a node next to NoData, valid and exact** (integer geometry:
  origin 500000 / 7900000, spacing 10 / 5; S4 covers a non-dyadic grid). For a NoData node at each
  of its eight neighbours in turn, sentinel and NaN both: the node is valid
  and `z == value_at` bit for bit. Include an interior node, a node on the
  last row, one on the last column, and all four grid corners, because the
  clamped last cell makes the NoData corner sit on the other side there. In
  C++ the points are `g.node(c)`, computed by the same `node` that `node_at`
  compares against, so any geometry works.
  A Python test that builds coordinates in NumPy uses integer or dyadic
  geometry (see "A limit of the exact rule").
- **S2. On a NoData node, refused.** Sentinel and NaN.
- **S3. Off-node points keep today's rule.** Refused when any corner is NoData:
  a cell-interior point; a point on a cell side with NoData across the side
  (this pins the deferral); `std::nextafter` of a node toward the NoData
  neighbour, in x and in y (this pins exactness). On a NoData-free raster, a
  seeded set of off-node points gives `z` bit-equal to today's formula. The
  test keeps a copy of today's expression as its reference, so a reordered
  expression in `sample.hpp` fails.
- **S4. Exact at every node, on a non-dyadic grid.** All nodes `g.node(c)` of
  a 200 x 300 grid with origin 500000.3 / 7900000.7 and spacing 0.7 / 0.3:
  `z == value_at` everywhere, far edges included. Today about 38 500 of
  60 000 fail (38 580 in the measurement above; 38 520 with `@tester`'s
  values). The assertion does not depend on the count. The count may differ with another compiler or
  standard library, but the green assertion holds everywhere: the test's
  points and `node_at`'s comparison both come from `RasterGeometry::node`,
  whatever that compiles to. A dyadic grid does not work here: it is
  exact today too (measured: 0 failures).
- **S5. `node_at` agrees with refine's node test** (invariant-critical for
  "the tolerance path is unchanged"). The reference is refine's whole test at
  `refine.hpp:417-426`: `v = lattice_position(g, p)`, then `v.is_node() &&
  g.node(c) == p` with `c` the node of `v`. `lattice_position(g, p).is_node()`
  alone is **not** the reference. It is true for an off-node point whose
  fractional coordinates round to integers, so it disagrees with a correct
  `node_at` (`@tester` found 1300 such probes on a 0.1 grid, 350 of them from
  `nextafter`). Over finite points in the node rectangle (nodes, `nextafter`
  neighbours of nodes, cell-side points, random points, and NumPy-style
  unfused node coordinates), on an integer grid, a dyadic grid and the
  non-dyadic "tenths" grid (origin 0.1 / 100.1, spacing 0.1, 120 x 120):
  `g.node_at(p).has_value()` equals the reference, and the indices agree when
  both say node. Also `node_at` returns nullopt for NaN, infinities and points
  outside the rectangle.
- **S6. Degenerate rasters unchanged.** A 1 x N and an N x 1 raster still give
  nullopt at their nodes. Points outside, and non-finite points, still give
  nullopt.
- **S7. Tolerance-path digests, default flags.** In `test_refine_golden.py`,
  add digests for `rasputin mesh --dem KARTVERKET --tolerance 1` with default
  flags (start quality 25, constraint feet on), tile and quarter circle,
  through `_cli_outcome`'s spy. They are **recorded in the red commit, from the
  pre-change production code**, and the green commit must not change them.
  The existing increment-17 digests stay as they are and must stay green.
- **S8. Stride-path digest on Kartverket.** SHA-256 of the `.vtk` that
  `rasputin mesh --dem KARTVERKET` writes at its default stride, with the
  version string kept out of the bytes. Again recorded in the red commit from
  pre-change code. It must be unchanged after green (see "Which paths change").
  The red step also pins `nodata_vertices_removed` = 397 at that stride, and
  it stays 397 after green: every removed vertex there is a NoData node.
  `@needs_codecs`.

Changed (one amendment commit on top of 25's tests, the reason in the message):

- `test_raster_view.cpp`, "bilinear_batch: a zero-weight NoData corner still
  makes the point invalid" and "bilinear_batch: NaN corners are NoData even
  with no sentinel": flip to valid, with the node value, and rename them. S1
  covers the cases; keep or fold them as `@tester` judges.
- `test_raster_view.cpp`, "z is the node value, exactly inside, to 1e-9 on the
  far edges": exact everywhere (`==`), renamed.
- `test_core_raster.py`, `TestSample::test_sentinel_and_nan_corners_are_not_valid`
  and `TestToCore::test_forwards_the_sentinel`: each samples a node next to a
  sentinel or NaN and expects it to be refused. Move the sentinel or NaN onto
  the sampled node, or sample an off-node point in that cell, so both tests
  still prove what their names say.
- Increment 25's `test_cli_mesh_plain_output.py`,
  `test_nodata_vertices_removed_counts_the_strided_nodes_without_data`: the
  expected count goes back to "picked nodes that are NoData" (25 at stride 1
  on the `holed` fixture), and the docstring says so. Steps 1, 2 and 3 stay,
  with literal counts 25, 9 and 4. Only stride 1 is red before green; the
  trim reaches no picked node at strides 2 and 3. Drop
  `trimmed_by_the_sampler` and its assertion that at stride 1 the trim
  removes more than the NoData nodes (both done in the red step, `6605dfe`).
- `test_the_no_tolerance_summary_says_on_or_next_to`, in both files that had
  it: renamed `test_the_no_tolerance_summary_says_on_nodata_cells`
  (`test_cli_mesh_plain_output.py:279`, `test_run_record.py:443`). Both expect
  `vertices on NoData cells were removed`, and also check that "next to" is
  absent from stderr.
- `test_cli_mesh_dem.py`, `test_a_nodata_edge_row_is_dropped_and_counted`: the
  wording goes back to `6 vertices on NoData cells`. The count stays 6: that
  fixture's NoData row is row 0, and under today's rule a node's cell reaches
  down and right, not up.
- `test_run_record.py` (25's): the stride summary says `on NoData cells`, and
  the parametrisation over `projected` / `stride` collapses to one wording.

`@tester` greps for any other test that pins a zero-weight refusal
(`grep -rn -i 'zero.weight\|next to NoData\|on or next to' tests`). The scan
test "a NoData corner with zero weight still voids an off-node vertex"
(`test_refinement_scan_offnode.cpp:302`) is about `vertex_z` on a cell side.
It **does not change**.

## Production changes for `@developer`

- `geometry.hpp`: `node_at`, as specified above.
- `sample.hpp`: the early return in `bilinear`; the header comment states the
  rule.
- `bindings/core.cpp` and `_core.pyi`: the `sample` docstrings.
- 25's `run_record.py`: the stride summary goes back to `on NoData cells`,
  one wording for both paths (the `tolerance is not None` choice and its
  comment go).
- `grid_domain.py`: the `subsample` docstring correction above.

Docs in the code PR, each a line saying increment 27 changed it:
`12-dem-to-mesh.md` R2, the "NoData" and "On the last row or column" bullets;
`12-dem-to-mesh.md:345` (test 4, "next to a NoData corner even at zero
weight", and the 1e-9 on the last row and column) and `:443` (the exclusion
"Changing `bilinear`'s NoData rule"). In `25-plain-output.md`: the inventory
row for `vertices without data dropped` (around line 121), the
"`nodata_vertices_removed` on the path without `--tolerance`" bullet (around
751), and the "Until that fix" block in "NoData on the no-tolerance path"
(around 858-894). The 16 text on input vertices on a node (`16-domain-polygon.md`,
around "An input vertex that sits exactly on a node") is still true and is not
touched. `ROADMAP.md`: the proposed sampler row that 25 adds becomes row 27,
updated at merge as the README requires.

## `@perf`

**Not required.** The README's rule covers diffs that touch
`include/terrain/refinement/`, `include/terrain/mesh/`, or what drives them.
This diff touches `include/terrain/raster/`. refine calls `bilinear` only in its
output loop, once per off-node start vertex, and S5/S7 show the result is
unchanged. The added cost per call is one division, one round and one compare
per axis. The stride path's `sample` phase is O(vertices) and was never a
bottleneck. `@reviewer` may quote the `sample` phase from `--stats` on the
Kartverket stride run before and after, as a sanity figure, not as
acceptance.

## LOC estimate

About 15 production lines: `node_at` about 10, the early return 3, and
`run_record.py` a net -1. Docstrings, raw literals and comments are excluded
under `CLAUDE.md` §2. Far under the ceiling.

## For Ola

1. **Scope: nodes only, cell sides later if ever (recommended).** It fixes the
   reported defect completely, since every vertex without `--tolerance` is a
   node. The tolerance path stays byte-identical, and no `@perf` round is
   needed. The alternative, cell sides too, means changing refine's `vertex_z`
   with it, new golden digests, and `@perf` acceptance. Increment 23's
   lattice-line cuts are where that would matter.
2. **Where the rule lives: in the C++ sampler (recommended)** rather than in
   the stride path. The alternative is that the stride path reads `value_at`
   by rounded index in Python. That works on any grid geometry, including the
   ones in "A limit of the exact rule", but it leaves `bilinear`'s trap in
   place for the next caller and adds a second "is this a node" relation.
   Real DEMs (integer origin and spacing) are not affected by the limit.
3. **For information, not a decision now.** Apple clang fuses
   `RasterGeometry::node`'s expression and NumPy does not, so on a non-dyadic
   grid geometry C++ and Python disagree on node coordinates in the last bit.
   That bites refine's start-vertex classification too, today. It is invisible
   on real DEMs. If it ever matters, the fix is project-wide
   (`-ffp-contract=off`, or `node` written so it cannot be fused). That
   touches every geometric computation, so it would be its own increment with
   `@perf`.

Recommended defaults: 1 and 2 as recommended; 3 noted.

## Review

**Design review, round 1, 2026-10-03.** Commit c277ce1. Verdict: CHANGES REQUESTED. LOC: 0 (design only); ~15 plausible (a `node_at` prototype is 8 lines). The rule, the unchanged `--tolerance` path, the 38,580 / 60,000 and 95 / 300 measurements, Kartverket's 0 off-node samples and the test plan all hold; every 25 test named exists on master. Blocking: (1) merge master; `cli.py:811-812` is now :830; (2) cite `grid_domain.py` 66-67 for the stride expression; (3) the status line and line 219 assume 25 not yet merged; (4) `test_the_no_tolerance_summary_says_on_or_next_to` exists in two files: name both. Not pushed; no CI.

**Design review, round 2, 2026-10-03.** Range `1477da0..8d8e394` (merge 9fb2218, prose 8d8e394). Verdict: CHANGES REQUESTED, one line: the stride expression is at `grid_domain.py:68-69`, not 66-67 (the round-1 number was the reviewer's slip). Everything else checked: the merge is clean, the new citations read as quoted, the round-1 items done. APPROVED once 68-69 is in, with no further round. Not pushed; no CI.

**Design review, round 2 condition met, 2026-10-03.** b50c8d0 cites `grid_domain.py:68-69` as the reviewer required; per round 2 the design is APPROVED with no further round. Not pushed; no CI.

## Ruled by Ola, 2026-10-03

Ola: "yes to both", to the two decisions as recommended. (1) Nodes only now: a point exactly on a DEM node reads that node; cell sides are deferred until a measured run shows ragged edges along NoData. (2) The rule lives in the C++ sampler (`bilinear` with `RasterGeometry::node_at`), not only in the Python stride path.

**27, code review, round 1, 2026-10-04.** Range `7eb0b79..a0817ba` (red 6605dfe, pins 80990e0, green 6fd4076 and a0817ba). Verdict: CHANGES REQUESTED. LOC: 13 net (`node_at` in `geometry.hpp` 11, the early return in `sample.hpp` 3, `run_record.py` -1) against about 15. Code sound; the notes in `12-dem-to-mesh.md` and `25-plain-output.md`, the binding docstring and the stub match the code; citations into this increment's files re-read and hold. Blocking: (1) red-step text in `tests/cpp/unit/test_raster_node_at.cpp@9a47922:16-18` and `tests/cpp/CMakeLists.txt@9a47922:389-390` (say why the target is separate now, without the red-step story); (2) this file's status line says "designed … Not started"; (3) `ROADMAP.md:50`'s sampler row still says "not designed", and this design assigns it to become row 27 in this PR. Not pushed; no CI.

**27, code review, round 2, 2026-10-04.** Range `9a47922..11759aa` (6c1369a test comments, 11759aa status line and ROADMAP row 27); whole increment `7eb0b79..11759aa`. Verdict: APPROVED. LOC: 13 net, unchanged. All three round-1 blockers closed: both comments now give the reason the target is separate (it includes `refine.hpp`, so it needs the backend target and Threads), the status line and `ROADMAP.md:50` match the branch; citations unchanged in line count and hold. Not pushed; no CI. At merge, ROADMAP row 27 says shipped.
