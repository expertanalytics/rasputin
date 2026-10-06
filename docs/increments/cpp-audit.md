# Audit: the C++ core — layering, duplication, bindings, dead code

Status: **audit, not a design.** Written by `@architect` against master at
`44fa7f5`. Read-only on code. Every `file:line` citation is pinned to
`44fa7f5`. It is the counterpart of `python-audit.md` (branch
`worktree-python-audit`, `8241e14`) and follows its shape. It proposes a
sequence of refactor PRs; each still runs the normal loop
(`docs/increments/README.md`): design note, red tests where behaviour changes,
code, review, and `@perf`'s byte-identical run where the diff touches refine or
mesh code (marked **[refine/mesh]** below).

Ola's ask: audit the code, keeping in mind that tests are the lesser problem;
the concern is losing control of the code base at about 10,000 lines each of
C++ and Python.

**Prior art.** None applies: this audits shipped code and claims nothing new.
A PR below that designs a new module carries its own prior-art section.

## Verdict in five lines

1. The C++ core is **6,052 code lines**, not 10,000 (8,827 non-blank), and its
   include graph is already close to layered: against the six-layer table of
   section 5 it has **one** violation. Python is the larger control risk.
2. The control risk in C++ is concentrated in `refinement/` (1,452 code lines,
   29 % of the headers) and in `bindings/core.cpp` (985 code lines in one
   file): the refinement round loop exists twice, and lattice and segment
   arithmetic is spelt out by hand in four to six places.
3. About **230 production lines** are reached only by tests (point-in-ring,
   signed area, `Point3`, `for_each_chunk`, a legacy border test, ...), and the
   Python package exposes `Point2`, `Point3`, `dot`, `cross` that no Python
   module uses.
4. The bindings convert arrays seven ways with different rules, and one has
   become a silent fault: `refine` drops a constraint edge with a negative
   index, where `indexed_mesh` refuses it (probe in C3).
5. Bilinear interpolation is written three times, once in Python; the Python
   copy differs from the C++ one in the last bits at 53 % of random points and
   in its NoData rule beside a node (probe in C6).

## 1. Size, re-measured

Counted with `tools/count_loc.py`'s own `counted_lines` (the `CLAUDE.md` §2
rule: no blank, comment, docstring or raw-literal line), on `44fa7f5`,
`lib/detria/detria.hpp` excluded:

| Part | Code lines | Non-blank |
|---|---|---|
| `include/terrain/` (47 headers) | 4,948 | 7,276 |
| `src/` (2 files, the two detria translation units) | 119 | 222 |
| `bindings/core.cpp` | 985 | 1,329 |
| **C++ production** | **6,052** | **8,827** |
| `tests/cpp/` (84 files) | 22,292 | 28,962 |
| `src_python/` (for comparison, same counter) | 7,687 | 9,924 |

By component (code lines): `refinement/` 1,452, `core/` 923, `noding/` 802,
`mesh/` 634, `vector_simplify/` 286, `raster/` 248, `hydrology/` 186,
`predicates/` 171, `cdt/` 134, `parallel_util/` 93. Largest files:
`bindings/core.cpp` 985, `refinement/refine_points.hpp` 402,
`refinement/refine.hpp` 331, `noding/noded_pslg_builder.hpp` 326,
`vector_simplify/area_collapse.hpp` 286. A third of the non-blank header
lines are comments (1,637 comment lines); 92 of them cite design-document
labels such as `L12`, `N7`, `D4`, readable only beside the increment file.

Boundaries that hold today, checked: no header includes `pybind11` (two
mention it in comments); no file in `include/`, `src/` or `bindings/` includes
`<fstream>`, `<filesystem>`, `<iostream>`, `<sstream>` or `<cstdio>`, or names
a CRS; `detria.hpp` is included only by the two `src/` files
(`tools/check_detria_boundary.py`).

## 2. Findings, ranked by lines saved times risk

Each finding: evidence, proposed shape, lines that would go (an estimate from
the cited ranges, net production lines under §2), risk, and the PR in section 6.

### C1. Production code that only tests reach — about 230 lines, low risk (PR B)

The only production consumer of the headers is `bindings/core.cpp`. These are
reached from no production file (search of `include/`, `src/`, `bindings/`
for each name outside its own definition):

- `core/ring.hpp`: `PointRing` (`include/terrain/core/ring.hpp@44fa7f5:83-102`), `edge(R, i)`,
  `all_finite`, `bounding_box(R)`, `signed_area` (`include/terrain/core/ring.hpp@44fa7f5:142-192`),
  `PointInRing` and `point_in_ring` (`include/terrain/core/ring.hpp@44fa7f5:250-306`): 71 code lines.
  Production uses only `IndexedRing` and `orientation<K>`.
- `Point3`, its operators, `dot`, `cross`, hash and formatter
  (`include/terrain/core/point.hpp@44fa7f5:37-68, 89-100, 121-139`): 51 code lines. No header uses
  `Point3`.
- The Python surface `Point2`, `Point3`, `dot`, `cross`
  (`bindings/core.cpp@44fa7f5:398-480`, 52 code lines), re-exported by
  `src_python/tin_engine/__init__.py@44fa7f5:5-7` and used by no module in
  `src_python/` or `tools/`. A side effect: because `__init__.py` imports
  `_core`, no `tin_engine` module can be imported without the extension.
  Probe: with `sys.modules['tin_engine._core'] = None`, `import
  tin_engine.stats` raises `ImportError`; the control without the block
  imports. `stats.py`'s docstring says it is "testable without the
  extension" (`src_python/tin_engine/stats.py@44fa7f5:3-4`); it is not.
- `parallel_util::for_each_chunk` (`include/terrain/parallel_util/chunks.hpp@44fa7f5:41-75`, 32 lines):
  both callers use `for_each_block`.
- `RasterGeometry::contains_strict` and `boundary_epsilon`
  (`include/terrain/raster/geometry.hpp@44fa7f5:135-162`, 13 lines): the legacy border draping they served
  is gone.
- `Box2::center`, `Box2::intersects` (`include/terrain/core/bbox.hpp@44fa7f5:88-93, 140-146`),
  `is_left_turn`, `is_collinear` (`include/terrain/predicates/orientation.hpp@44fa7f5:79-88`),
  `reversed(Orientation)`, `reversed(Segment2)`, `is_degenerate(Segment2)`:
  about 15 lines.

Not dead, kept: `NodedPslg::edge_base` and `CheckPoints::reserved_points`
are test seams that observe real invariants; `FastKernel` is the measured
baseline the filtered kernel is tested against; `indexed_mesh`,
`refine_seam`, `SeamOutcome` and `frozen_mask` have no Python production
caller on master but are used by `src_python/tin_engine/basin_run.py` on the
in-flight `worktree-23c` branch.

Shape: delete, with their tests. `PointRing` goes with `point_in_ring`'s
tests, or the tests build an `IndexedRing` over `0..n-1`. Risk: low; nothing
in production changes. The Python surface is a public-API question for Ola.

### C2. The refinement round loop exists twice — about 60 lines, high risk [refine/mesh] (PR F)

`refine` (`include/terrain/refinement/refine.hpp@44fa7f5:297-433`) and `detail::point_loop`
(`include/terrain/refinement/refine_points.hpp@44fa7f5:248-458`) run the same skeleton: build the
lattice mesh and set the frozen mask, build the frame, `legalise_all`, seed
`active`, then per round a `for_each_block` scan into one result per slot, a
serial split pass (`touched`, `skipped`, `seeds`, the neighbour-touched skip,
`split_inside` or `split_edge`, `legalise_around`, counters), `rebuild_active`,
and after the loop the same output assembly (vertices to world, `z`, `valid`,
`triangles`, `constraint_edges`). `include/terrain/refinement/refine_points.hpp@44fa7f5:3-5` says why: "refine.hpp's loop
with a different scanner; refine.hpp itself is not edited (J1)", a rule of
increment 15c. The tolerance check is written three times
(`include/terrain/refinement/refine.hpp@44fa7f5:294-296`, `include/terrain/refinement/refine_points.hpp@44fa7f5:476-478, 496-498`), and `refine_seam`
refuses the same input by throwing instead (`include/terrain/refinement/seam.hpp@44fa7f5:66-67`).

`point_loop` is also two functions in one: `refine_points` and `refine_strip`
select their halves through null pointers and `if constexpr` on a `NoSet`
sentinel (`include/terrain/refinement/refine_points.hpp@44fa7f5:207, 225`).

Shape: one `round_loop(mesh, frame, options, scan_one, choose_split,
on_insert)` in `refinement/round.hpp`, used by both; `refine` passes `scan`,
`point_loop` its combined scanner. Risk: high. It is the hot path, so `@perf`
must show byte-identical meshes and no timing loss on the 1 m set and the
thread sweep. It needs Ola's ruling to lift J1 (question 2).

### C3. The bindings convert arrays seven ways, and one way drops a constraint silently — about 60 lines, medium risk (PR C)

- An (N, 2) float64 array is read by `as_points` (`bindings/core.cpp@44fa7f5:186-216`, refuses a
  path with `TypeError`), `xy_points` (`bindings/core.cpp@44fa7f5:324-332`), `start_mesh`
  (`bindings/core.cpp@44fa7f5:358-360`), `sample` (`bindings/core.cpp@44fa7f5:914-917`), `constraint_check_points`
  (`bindings/core.cpp@44fa7f5:1151-1157`) and `CheckPoints.add` (`bindings/core.cpp@44fa7f5:1099-1101`, exact dtype, no
  cast), with six different messages and two dtype rules.
- Index arrays are read as `uint32` with `forcecast` in `refine`
  (`bindings/core.cpp@44fa7f5:1043-1048`), `start_mesh` and `constraint_check_points`, but as checked
  `int64` in `indexed_mesh` (`bindings/core.cpp@44fa7f5:697-713`), whose comment gives the reason:
  "a negative index ... is refused rather than wrapped". The other entry points
  wrap. Probe (installed extension, 5 x 5 DEM, two triangles):
  `refine(..., edges=[[0, 1]])` returns `Ok` with constraint edge `[[0, 1]]`;
  `edges=[[0, -1]]` returns `Ok` with **no** constraint edge, and no error;
  `indexed_mesh` with a `-1` raises `ValueError`.
- `U32` and `U8` are declared at file scope (`bindings/core.cpp@44fa7f5:334-336`) and again
  inside `refine` and `upstream` (`bindings/core.cpp@44fa7f5:1043, 1309`).
- The pattern "release the GIL, `std::visit` over the float/double view" is
  written out at nine `std::visit` sites; the (T, 3) triangle view twice
  (`bindings/core.cpp@44fa7f5:657-669, 984-994`).
- Two comments are now false: "The GIL is released around exactly one call:
  triangulate" (`bindings/core.cpp@44fa7f5:385-391`), where the file has 12
  `gil_scoped_release`; and the binding surface of `06-cdt-viewer.md` being
  "exhaustive" (`bindings/core.cpp@44fa7f5:482-488`), when the file has since grown the raster,
  refine, check-point, seam, hydrology and ring-reduction entry points.

Shape: a private `bindings/marshal.hpp` with `points(obj, name)`,
`indices(obj, name, shape, bound)` (int64 read, negative and out-of-range
refused, one message form), `with_view(raster, f)` (GIL and visit), and
`readonly_rows<T, K>(owner, span)`. Optionally split `core.cpp` by component
(`bind_cdt.cpp`, `bind_refine.cpp`, `bind_hydrology.cpp`, `module.cpp`) so no
binding file passes about 400 lines. Risk: medium. The negative-index refusal
is a behaviour change, so it needs a red test first (`@tester`), and refusal
wordings pinned by suites either stay or change in a `@tester` amendment.

### C4. Lattice and segment arithmetic copied through refinement — about 40 lines, medium risk [refine/mesh] (PR D)

The C++ counterpart of the Python audit's F2:

- World point to lattice position: `detail::lattice_position`
  (`include/terrain/refinement/refine.hpp@44fa7f5:132-142`) and `RasterGeometry::node_at`
  (`include/terrain/raster/geometry.hpp@44fa7f5:98-108`) are the same clamp-and-round; the second's
  comment says it "is written as refine's lattice_position writes it".
  `CheckPoints::add` divides a third time (`include/terrain/refinement/check_points.hpp@44fa7f5:82`).
- Lattice position to world, `x_min + col dx`, `y_max - row dy`, for an
  off-node vertex: `include/terrain/refinement/refine.hpp@44fa7f5:424`, `include/terrain/refinement/refine_points.hpp@44fa7f5:452`, `include/terrain/refinement/seam.hpp@44fa7f5:131`.
- Projection of a point onto an edge in (col, row): `include/terrain/refinement/scan.hpp@44fa7f5:112-115`,
  `include/terrain/refinement/strip_scan.hpp@44fa7f5:149-151`, `include/terrain/refinement/strip_scan.hpp@44fa7f5:205-207`, `include/terrain/refinement/refine_points.hpp@44fa7f5:142-143`, and in world units
  `include/terrain/refinement/refine.hpp@44fa7f5:254-256`. The copies already cost code: `near_constraint` excludes
  frozen edges a second time because the scan's copy is "a copy of this
  expression that the compiler may contract differently"
  (`include/terrain/refinement/strip_scan.hpp@44fa7f5:196-198`). The build sets no `-ffp-contract`
  (no match for `contract` in `CMakeLists.txt`).
- Linear z along an edge: `along` (`include/terrain/refinement/strip_scan.hpp@44fa7f5:76-80`), the seam pass
  (`include/terrain/refinement/seam.hpp@44fa7f5:108`, "as the strip's along"), the frozen error (`include/terrain/refinement/refine_points.hpp@44fa7f5:145`).
- The cell of a fractional coordinate, `min(floor(f), n - 2)`:
  `include/terrain/refinement/check_points.hpp@44fa7f5:52` (`last_cell`), `include/terrain/refinement/scan.hpp@44fa7f5:90-91` (`vertex_z`),
  `include/terrain/raster/geometry.hpp@44fa7f5:130-131`.
- Layering cost: `constraint_points.hpp` and `seam.hpp` include `refine.hpp`
  for `detail::lattice_position` alone, and `scan.hpp` for `vertex_z` alone
  (`include/terrain/refinement/constraint_points.hpp@44fa7f5:102, 149`; `include/terrain/refinement/seam.hpp@44fa7f5:80`), so the edge-strip generator and the
  seam pass compile against the whole refine loop, `quality.hpp`, `lawson.hpp`
  and the thread pool.

Shape: `refinement/frame.hpp`, below `scan.hpp`: `lattice_position(g, p)`
(defined once; `RasterGeometry::node_at` keeps its own spelling or calls a
shared `raster` helper — the two are in different layers), `world(g, v)`,
`vertex_z(dem, v)`, and in `mesh/`: `project(a, b, p) -> {sigma, distance}`
and `lerp_along`. Risk: medium. Moving an expression into an inline function
can change floating-point contraction, so `@perf`'s run must show every mesh
byte-identical on the 1 m set and the reprojected golden suite.

### C5. Mesh traversal written several times — about 20 lines, medium risk [refine/mesh] (in PR D)

- "Which slot of neighbour u points back at t": `include/terrain/mesh/lattice_mesh.hpp@44fa7f5:248-250`,
  `include/terrain/mesh/lattice_mesh.hpp@44fa7f5:269-271`, `include/terrain/mesh/lawson.hpp@44fa7f5:134-136`, `include/terrain/refinement/strip_scan.hpp@44fa7f5:181-183`, `include/terrain/refinement/refine.hpp@44fa7f5:280-282`.
  Shape: `LatticeMesh::back_slot(u, t)` and `apex_across(t, e)`.
- The undirected edge key four ways: `std::uint64_t` in `include/terrain/refinement/strip_scan.hpp@44fa7f5:70-72`,
  `std::minmax` pairs in `include/terrain/refinement/refine.hpp@44fa7f5:164-175`, a pair over a run in
  `include/terrain/noding/node.hpp@44fa7f5:141-146`. One `edge_key` in `core/`.
- The flip test: `detail::must_flip`'s frame branch (`include/terrain/mesh/lawson.hpp@44fa7f5:156-159`) and
  `strip_fits`'s `flips` lambda (`include/terrain/refinement/strip_scan.hpp@44fa7f5:168-174`). The copy omits the
  integer fast path `lattice_incircle`; the answer is the same by
  `include/terrain/mesh/lawson.hpp@44fa7f5:83-91`'s argument, but the two can now drift. Shape:
  `mesh::detail::frame_flips(a, b, c, d, frame)` used by both.

### C6. Bilinear interpolation three times, one in Python, with measured drift — 0 C++ lines, about 25 Python lines; Ola's ruling (PR G)

- `raster::bilinear` (`include/terrain/raster/sample.hpp@44fa7f5:27-58`), on a world point, with
  increment 27's rule that a point exactly on a node reads that node alone.
- `refinement::vertex_z` (`include/terrain/refinement/scan.hpp@44fa7f5:79-98`), on a lattice position:
  the same formula and order. Two C++ copies are defensible (one takes world
  points, one lattice positions), but should share the four-corner formula.
- `target_grid.resample` (`src_python/tin_engine/target_grid.py@44fa7f5:176-196`), in numpy, with
  the factored form `(1 - fy)((1 - fx) a + fx b) + fy(...)`, the `isfinite`
  NoData rule, and no node exception. Probe (installed extension, 50 x 50
  float32 DEM of 100-2000 m, 100,000 random points): the numpy expression and
  `_core.sample` differ bitwise at 53,032 points, by at most 6.4e-12 m. And on
  a 3 x 3 grid with NaN at row 0, column 1, `_core.sample` at the node
  (row 0, column 0) is valid; `resample`'s stencil test refuses the same node.

Shape: `resample` calls `_core.sample` on each block's points over a
`raster_view` of the window (it already releases the GIL), or the Python audit's
F9 `valid_mask` ruling adopts the C++ rule and a cross-language test pins the
two together. Either way reprojected meshes change in the last bits, so it is
Ola's call (question 3) and needs `@perf`.

### C7. `Pslg` and `NodedPslg` repeat their accessors — about 25 lines, low risk (PR E)

`include/terrain/core/pslg.hpp@44fa7f5:102-156` and `include/terrain/core/noded_pslg.hpp@44fa7f5:77-115` (26 code lines
each) hold the same seven accessors over the same three vectors (`vertices`, `chains`, `chain_indices`,
`indices_of`, `ring`, `edge_count`, `edge`), line for line; the bindings
already bind them with one template for that reason (`bindings/core.cpp@44fa7f5:130-181`).
Shape: a `ChainSet` value type in `core/` holding the three vectors and the
accessors, a private member of both; both keep their distinct types, so
"holding a `NodedPslg` proves noding ran" survives. Risk: low; no arithmetic
changes. Not refine or mesh code.

### C8. Concepts and templates: where they earn their keep and where not — about 15 lines (in PRs B and D)

- **Earn it:** `RasterSource` (two production instantiations, `RasterView<float>`
  and `RasterView<double>`, plus `Raster<T>` in tests); `CdtBackend` behind
  `triangulate` (two test doubles, `FakeCdtBackend` and `FailingCdtBackend`);
  `ExactPredicates` (one model, kept so the exact backend can be swapped in
  one line, as the `modern-cxx` skill asks).
- **Nominal:** `GeometryKernel K` on `legalise_all`, `legalise_around` and
  `improve` has one production instantiation, `DefaultKernel`, while the same
  layer hard-wires `DefaultKernel` in `orient_sign`
  (`include/terrain/mesh/lattice_mesh.hpp@44fa7f5:107-109`) and `strip_fits` (`include/terrain/refinement/strip_scan.hpp@44fa7f5:164`), and
  `must_flip`'s integer shortcut is correct only for an exact `K`
  (`include/terrain/mesh/lawson.hpp@44fa7f5:139-143`). The parameter's one other use is a call-counting kernel in
  one test. Shape: leave the parameter (it costs nothing) but state on
  `legalise_*` that `K` must be exact; do not add a second kernel to mesh code.
- **Ring:** a concept with two models, of which production uses one
  (`IndexedRing`); goes with C1.
- **Raster:** `Raster<T>` (`include/terrain/raster/raster.hpp@44fa7f5:37-77`) and `RasterView<T>`
  (`include/terrain/raster/view.hpp@44fa7f5:14-46`) repeat `value_at`, `nodata`, `is_nodata`, `row`;
  the NoData rule is written a third time in the scan
  (`include/terrain/refinement/scan.hpp@44fa7f5:183`). Shape: `Raster<T>` owns a vector and hands out a
  `RasterView<T>`; one `is_missing(v, nodata)` (about 15 lines).
- **The opposite case:** C2's two loops are where an abstraction is missing.

### C9. Layering: one sibling edge, and `refine.hpp` as a hub (PR A)

Against the table of section 5, checked by a script over every
`#include <terrain/...>` line: one violation,
`vector_simplify/area_collapse.hpp` includes `noding/intersect.hpp`
(`include/terrain/vector_simplify/area_collapse.hpp@44fa7f5:18`), a sibling component. `intersect.hpp`
depends only on `core/` and `predicates/`, so it moves down to `core/`
(`core/segment_intersect.hpp`). Inside `refinement/`, `refine.hpp` is
included by four of its six siblings, two of them only for one helper (C4).
`NodedPslg` names `noding::NodedPslgBuilder` as a friend by forward
declaration (`include/terrain/core/noded_pslg.hpp@44fa7f5:65-67`); not an include, so the gate
allows it, but it is the one place a lower layer knows a higher one by name.

## 3. Tests (light; `@tester` audits the suites)

- **The stub has no drift gate.** `_core.pyi` is 687 lines written by hand,
  161 of them docstrings that repeat the binding docstrings. `mypy.stubtest`
  against the installed extension reports 88 errors, all pybind11 artefacts
  (metaclass, enum member types, `__init__(*args)`, `__index__`,
  `__members__`); none is a missing name or a wrong parameter. A stubtest run
  with an allowlist for those categories would catch a name or signature that
  drifts, at the cost of a few seconds in the Python job.
- **Cross-language contracts with no test.** Python and C++ must agree bit for
  bit on the node coordinate `x_min + col dx` (refine decides "node" by
  `g.node(c) == p`, `include/terrain/refinement/refine.hpp@44fa7f5:137-142`), and should agree on NoData and bilinear
  (C6). A property test calling both sides on random grids pins each.
- C1 deletes the tests of what it deletes.

## 4. The include-layering gate

Three ways to enforce section 5's table:

1. **A governance script, `tools/check_cpp_layers.py`** (recommended). It parses
   `#include <terrain/...>` directives — not prose, as
   `check_detria_boundary.py` and `check_prohibited_deps.py` already do —
   maps each header to its layer and component by a table in the script,
   and fails on an include upward or across siblings, on a header missing from
   the table, and on a `pybind11` or file-I/O include under `include/`. It runs
   in the governance job, which runs on every PR in seconds and needs no
   build. A pytest beside it plants a violation in a temporary tree and checks
   the script names it. Cost: about 80 lines of tooling. `tools/check_*` is a
   governed path, so Ola approves the file.
2. **A pytest only, `tests/python/test_cpp_layering.py`**, the twin of the
   Python audit's T2. Not governed, but it runs only in the Python job, after
   the extension build.
3. **CMake**, one interface target per layer with its own include root. It
   would make a violation fail to compile, but means moving every header and
   rewriting every include; a tool such as include-what-you-use adds a
   dependency. Not proposed.

A probe of the table (a 40-line script run against `44fa7f5`) found every
header listed and the one violation of C9, so the gate lands green once
`intersect.hpp` moves.

## 5. Target architecture (one page)

Six layers by file, not by directory: `core/` today spans three of them.
A header includes only from its own component or a lower layer; components on
one layer do not include each other.

```
L5  bindings        bindings/*.cpp: the only pybind11; marshal.hpp converts
                    arrays once; one file per component, module.cpp registers
L4  refinement      frame.hpp (new: lattice_position, world, vertex_z)
                    < scan.hpp < round.hpp (new: the one round loop)
                    < refine.hpp, check_points.hpp, constraint_points.hpp
                    < strip_scan.hpp < refine_points.hpp; seam.hpp
L3  algorithms      noding/ | cdt/ | mesh/ | hydrology/ | vector_simplify/
                    (siblings; each may use L0-L2 only)
L2  geometry        core/segment, core/segment_intersect (moved from noding/),
                    core/ring (IndexedRing, orientation), core/chain_set (new),
                    core/pslg, core/noded_pslg, core/pslg_builder,
                    core/indexed_mesh, raster/raster, raster/view,
                    raster/sample, raster/row_segments
L1  predicates      predicates/exact, kernel, detria_exact, default_kernel
                    (detria.hpp only in src/predicates and src/cdt)
L0  values          core/point, core/bbox, core/edge_properties,
                    core/snap_grid, predicates/orientation, raster/geometry,
                    parallel_util/chunks, build_info
```

What changes against today: `intersect.hpp` moves down; the strip generator and
the seam pass stop including `refine.hpp`; `refine` and the final check share
one loop; lattice and segment arithmetic is written once per layer; the
bindings convert arrays in one place and are split by component; about 230
lines that only tests reach are gone. `refinement/` should end near 1,300 code
lines and no binding file above about 400.

What keeps control as the code grows, beyond this round: the gate of section 4
on every PR, and one rule for new refinement work — a new pass adds a scanner
or a split policy to `round.hpp`'s loop, never a third copy of the loop.

## 6. PR order

Each PR is a refactor: behaviour-preserving unless it says otherwise, and
well under 700 net lines. None waits for the in-flight 23c-2 branch: its
commits touching C++ are merges only (`git log master..worktree-23c -- include
src bindings`), though PR C and 23c-2 both touch `_core.pyi`.

| # | PR (branch name) | Takes | Net production lines | Waits for | Gates beyond review |
|---|---|---|---|---|---|
| A | `cpp-layering-gate` | section 4, C9 | about +80 tooling (not counted), 0 production | Ola's approval of the governed file | planted-violation test |
| B | `cpp-dead-code` | C1, C8 (Ring) | about -230 | Ola's answer on the Python surface | none |
| C | `cpp-bindings-marshal` | C3 | about -60 | nothing | red test for the negative-index refusal |
| E | `cpp-chain-set` | C7 | about -25 | A (new header in the table) | none |
| D | `cpp-lattice-frame` | C4, C5, C8 (Raster) | about -55 | A | `@perf`: byte-identical meshes |
| F | `cpp-refine-round` | C2 | about -60 | D, Ola's ruling on J1 | `@perf`: byte-identical meshes, no timing loss |
| G | `bilinear-once` | C6, section 3 cross-language tests | about -25 Python | Ola's ruling; the Python audit's PR A | `@perf` on the reprojected path |

Total: about -450 production lines (C++ about -430), the bindings' silent
drop fixed, and the drift points (lattice position, edge projection, flip
test, NoData rule, bilinear) each written once.

## Questions for Ola

1. **The Python `Point2`, `Point3`, `dot`, `cross`.** No module uses them;
   they make every `tin_engine` import need the compiled extension. Default:
   remove them from the package and the bindings (PR B).
2. **Increment 15c's rule J1, "refine.hpp is not edited".** It is why the
   round loop exists twice. Default: lift it for PR F only, with `@perf`'s
   byte-identical run as the condition.
3. **One bilinear.** Making Python's resample use the C++ rule changes
   reprojected meshes in the last bits and at nodes beside NoData. Default:
   adopt the C++ rule (node exception, NaN-or-sentinel NoData) in PR G.
4. **The gate's home.** Default: `tools/check_cpp_layers.py` in the governance
   job (governed, so your approval per edit); the alternative is an ungoverned
   pytest.

## Rulings

2026-10-05, Ola, verbatim: "defaults on all, CI speed first". So, question
by question:

1. **Remove the Python `Point2`, `Point3`, `dot` and `cross`** from the
   package and the bindings (PR B).
2. **Increment 15c's rule J1** ("refine.hpp is not edited") is lifted for
   PR F only, and only if `@perf` shows byte-identical meshes on the 1 m set
   and the thread sweep, with no timing loss.
3. **One bilinear:** Python's resample adopts the C++ rule (a point on a
   node reads that node alone; NaN-or-sentinel NoData), in PR G, with
   `@perf` on the reprojected path.
4. **The include-layering gate** is `tools/check_cpp_layers.py`, in the
   governance job. It is a governed file, so Ola approves its edits.

CI speed comes first (`docs/increments/h17-ci-test-time.md`); these PRs
follow it.
