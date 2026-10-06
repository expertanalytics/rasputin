# Audit: the C++ core — layering, duplication, bindings, dead code

Status: **audit, not a design.** Written by `@architect` against master at
`44fa7f5`. Read-only on code. Every `file:line` citation is pinned to
`44fa7f5`. It is the counterpart of `python-audit.md` (branch
`worktree-python-audit`, `8241e14`) and follows its shape. It proposes a
sequence of refactor PRs; each still runs the normal loop
(`docs/increments/README.md`): design note, red tests where behaviour changes,
code, review, and `@perf`'s byte-identical run where the diff touches refine or
mesh code (marked **[refine/mesh]** below). Section 7 is the design of
PR B, written against `1035690`. PR B status: code and docs in; next `@reviewer`.

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
| B | `cpp-dead-code` | C1, C8 (Ring) | -291 (section 7) | nothing (ruling 1) | none (section 7) |
| C | `cpp-bindings-marshal` | C3 | about -60 | nothing | red test for the negative-index refusal |
| E | `cpp-chain-set` | C7 | about -25 | A (new header in the table) | none |
| D | `cpp-lattice-frame` | C4, C5, C8 (Raster) | about -55 | A | `@perf`: byte-identical meshes |
| F | `cpp-refine-round` | C2 | about -60 | D, Ola's ruling on J1 | `@perf`: byte-identical meshes, no timing loss |
| G | `bilinear-once` | C6, section 3 cross-language tests | about -25 Python | Ola's ruling; the Python audit's PR A | `@perf` on the reprojected path |

Total: about -450 production lines (C++ about -430), the bindings' silent
drop fixed, and the drift points (lattice position, edge projection, flip
test, NoData rule, bilinear) each written once.

## 7. PR B design: `cpp-dead-code` (finding C1, ruling 1)

Written by `@architect` against `1035690` (h17 PR 1's approved head; the
branch `worktree-cpp-dead` starts there and lands after it). Every
production file and every test file cited below at `44fa7f5` is
byte-identical at `1035690`: `git diff --quiet 44fa7f5 1035690 -- include
src bindings src_python tests/cpp/unit tests/cpp/support/ring_cases.hpp
tests/cpp/support/mesh_queries.hpp tests/cpp/property/prop_ring_invariants.cpp
tests/cpp/property/prop_cdt_invariants.cpp tests/python/test_core.py
tests/python/test_core_cdt.py tests/python/test_features.py
tests/python/test_io_vtk_legacy.py project_structure.md` exits 0. So the
citations stay pinned to `44fa7f5`, as in the rest of this file; the two
files h17 PR 1 changed, `tests/cpp/CMakeLists.txt` and the workflow, are
cited at `1035690`.

**Prior art.** None: a deletion claims nothing new and builds on no method.
Legacy: nothing is carried across. `git grep -l -E
'point_in_ring|contains_strict|for_each_chunk|Point3' legacy-archive --
legacy` returns one file, `legacy/rasputin/triangulate_dem.h`, and only for
CGAL's `Point_3` alias and its trait specialisations; the legacy border
draping `contains_strict` was written for is not in this tree (C1).

### 7.1 What goes, and why each has no production caller

"Production caller" means a use in `include/`, `src/`, `bindings/`,
`src_python/` or `tools/`, outside the item's own definition. Each was
searched by name at `1035690`, on every local branch with a worktree
(including `worktree-23c`, whose `basin_run.py` imports only
`IndexedMesh2, RefineOutcome, SeamOutcome, indexed_mesh, refine_seam` from
`_core`), and found only in comments. The one production user of the
Python names is `src_python/tin_engine/__init__.py@44fa7f5:5`, which this PR
rewrites. Counted lines are `tools/count_loc.py`'s `counted_lines` over the
range at `1035690`.

| # | Item | Range | Counted |
|---|---|---|---|
| R1 | the `Ring` concept and its comment | `include/terrain/core/ring.hpp@44fa7f5:48-56` | 5 |
| R2 | `PointRing` | `include/terrain/core/ring.hpp@44fa7f5:83-102` | 12 |
| R3 | `edge(R, i)`, `all_finite`, `bounding_box(const R&)`, `signed_area` | `include/terrain/core/ring.hpp@44fa7f5:142-192` | 31 |
| R4 | `template <Ring R>` on `detail::extreme_vertex` | `include/terrain/core/ring.hpp@44fa7f5:207` | 1 |
| R5 | `PointInRing`, `point_in_ring` | `include/terrain/core/ring.hpp@44fa7f5:250-304` | 27 |
| P1 | `Point3`, its operators, `dot(Point3)`, `cross(Point3)` | `include/terrain/core/point.hpp@44fa7f5:37-68` | 26 |
| P2 | `std::hash<Point2>` and its comment | `include/terrain/core/point.hpp@44fa7f5:74-87` | 8 |
| P3 | `std::hash<Point3>` | `include/terrain/core/point.hpp@44fa7f5:89-100` | 12 |
| P4 | `std::formatter<Point3>` | `include/terrain/core/point.hpp@44fa7f5:121-138` | 13 |
| B1 | `Box2::center` and its comment | `include/terrain/core/bbox.hpp@44fa7f5:83-92` | 3 |
| S1 | `is_degenerate(Segment2)` | `include/terrain/core/segment.hpp@44fa7f5:43-46` | 3 |
| O1 | `is_left_turn`, `is_collinear` | `include/terrain/predicates/orientation.hpp@44fa7f5:79-88` | 6 |
| C1 | `parallel_util::for_each_chunk` | `include/terrain/parallel_util/chunks.hpp@44fa7f5:41-75` | 32 |
| G1 | `RasterGeometry::contains_strict`, `boundary_epsilon` | `include/terrain/raster/geometry.hpp@44fa7f5:135-162` | 13 |
| K1 | `#include <pybind11/operators.h>` (only the point classes use `py::self`) | `bindings/core.cpp@44fa7f5:2` | 1 |
| K2 | `using terrain::Point3;` | `bindings/core.cpp@44fa7f5:54` | 1 |
| K3 | the `Point2`, `Point3`, `dot`, `cross` bindings | `bindings/core.cpp@44fa7f5:398-480` | 52 |
| Y1 | the `Point2`, `Point3` stubs | `src_python/tin_engine/_core.pyi@44fa7f5:15-56` | 36 |
| Y2 | the `dot`, `cross` stubs | `src_python/tin_engine/_core.pyi@44fa7f5:61-72` | 8 |
| I1 | the re-export | `src_python/tin_engine/__init__.py@44fa7f5:5` (line 7 becomes `__all__ = ["installed_version"]`) | 1 |

Two items are beyond C1's list, found while re-verifying it:

- **P2, `std::hash<Point2>`.** Its one production user is the Python
  `__hash__` (K3). With K3 gone only `tests/cpp/unit/test_point.cpp` reaches
  it, and two comments already forbid it as a dedup key
  (`include/terrain/core/snap_grid.hpp@44fa7f5:39-42`,
  `include/terrain/core/pslg_builder.hpp@44fa7f5:37-43`). `formatter<Point2>`
  stays: `pslg_builder.hpp` formats a vertex into a diagnostic.
- **R1, R4, the `Ring` concept** (C8). With `PointRing` gone it has one model.
  `detail::extreme_vertex` and `orientation<K>` take `const IndexedRing&`;
  `orientation` keeps its kernel parameter `K`: `PslgBuilder::build<K>` passes
  it through, and the suites run it under `FastKernel` and `DefaultKernel`.

Three items on C1's list **stay**, because a test of live code uses each as
its instrument, and moving them into `tests/cpp/support/` would split a
value type's vocabulary to save about ten lines:

- `reversed(Orientation)` (`include/terrain/predicates/orientation.hpp@44fa7f5:67-69`):
  the symmetry checks of the three kernels and of `prop_predicates_symmetry`.
- `reversed(Segment2)` (`include/terrain/core/segment.hpp@44fa7f5:39-41`): the
  noder's `classify` and `segment_meets_cell` symmetry tests.
- `Box2::intersects` (`include/terrain/core/bbox.hpp@44fa7f5:140-147`): the
  all-pairs oracle of the invariant-critical broad-phase suite
  (`prop_noding_broad_phase.cpp`).

Found and left for PR E: `Pslg::ring` and `NodedPslg::ring` have no
production caller either (production builds `IndexedRing` directly,
`include/terrain/core/pslg_builder.hpp@44fa7f5:320`,
`include/terrain/noding/node.hpp@44fa7f5:207`), but they are two of the
seven accessors PR E moves into `ChainSet`, and `prop_cdt_invariants.cpp`
uses one. Not checked at this grain: which of `Point2`'s own operators
(`operator/`, unary minus, `double * Point2`) production still uses; the
build would say, and the lines are few.

**Includes stay.** No `#include` line is removed except K1. `ring.hpp`'s
`bbox.hpp`, `segment.hpp`, `<cmath>` and `<concepts>`, `point.hpp`'s
`<functional>` and `geometry.hpp`'s `<cmath>` may become unneeded, but
suites compile through them transitively, and `@developer` may not edit a
test to add the include a suite was borrowing. Header hygiene belongs with
the layering gate (PR A).

**Comments that name a removed item** are rewritten in the same commit, each
in place and keeping its line count where a citation points below it:
`ring.hpp`'s header (`include/terrain/core/ring.hpp@44fa7f5:3-30`, on
`PointRing`, the two-model `bounding_box` ambiguity and the
`point_in_ring<DefaultKernel>` example), `include/terrain/core/segment.hpp@44fa7f5:5-9`
(`point_in_ring`'s parity rule), `include/terrain/core/pslg_builder.hpp@44fa7f5:16` (`PointRing`)
and `:305` (`signed_area`), `include/terrain/core/snap_grid.hpp@44fa7f5:39-42` (`std::hash<Point2>`;
**same line count**: `snap_grid.hpp:43`, `:63` and `:94` are cited unpinned),
`include/terrain/parallel_util/chunks.hpp@44fa7f5:3-30` (the header describes `for_each_chunk` first),
`src_python/tin_engine/_core.pyi@44fa7f5:1-6` (the docstring's reason is `cross`'s overloads), and
`src_python/tin_engine/stats.py@44fa7f5:3-4`, which becomes true and needs no edit.

`project_structure.md` lines 14, 17, 52, 446, 501 and 505 name `Point3`,
`PointRing`, `point_in_ring`, `for_each_chunk` and the `__all__` list;
`@developer` edits each **in place, one line for one line** (line 505 may go),
because `project_structure.md:166`, `:190`, `:208`, `:359` and `:384` are
cited unpinned from other files.

### 7.2 Refine or mesh code: no `@perf` run

Nothing under `include/terrain/refinement/` or `include/terrain/mesh/`
changes. Two edited headers are on refine's path: `chunks.hpp` (refine's
scan calls `for_each_block`) and `raster/geometry.hpp` (`RasterGeometry`).
In both the deletion is of something no translation unit instantiates or
calls (a function template, and two non-virtual inline members that change
no layout), so the compiled refine is the same code. `@reviewer` checks it
by `git diff 1035690 HEAD -- include/terrain/parallel_util/chunks.hpp
include/terrain/raster/geometry.hpp`: no `-` or `+` line outside the
comment header, `for_each_chunk` and the two members. If that diff shows
anything else, `@perf`'s byte-identical run is owed.

### 7.3 Lines

Net production lines by `CLAUDE.md` §2: **-291** (C++ -246, Python -45),
the table's sum. Rewritten comments, the modified signatures on R4 and
`orientation` (one line out, one in) and line 7 of `src_python/tin_engine/__init__.py` count zero.
Check: `python3 tools/count_loc.py 1035690 HEAD` on the branch prints -291;
a different number means an include or an operator went that this design
did not list, and the difference is explained in the review.

### 7.4 `@tester`: the suites

Rule: **a test that tests only removed code is deleted, not kept or
rewritten.** A test of live code that used a removed item as its instrument
is kept and given another instrument. Every change below lands in
`@tester`'s commit, before `@developer`'s, and must compile and pass against
`1035690`'s production code (it only stops using things).

Deleted outright:

- `tests/python/test_core.py` (all of it tests the Python `Point2`,
  `Point3`, `dot`, `cross`).
- `tests/cpp/unit/test_refinement_chunks.cpp` and its registration
  (`tests/cpp/CMakeLists.txt@1035690:224` and `:227`, with the comment at
  `:216-223` rewritten), and its name in the TSan job's build and run lists
  (`.github/workflows/main.yaml@1035690:155` and `:169`, each edited **in
  place** so `.github/workflows/main.yaml:306`, cited unpinned, does not move).
  `tests/python/test_ci_changes.py`'s `test_h3_tsan_builds_exactly_the_suites_it_runs`
  fails if the workflow names a suite CMake no longer registers, so the
  CMake and workflow edits go in one commit. `for_each_block` keeps its
  TSan coverage through `test_refinement_chunks_dynamic`.
- In `tests/cpp/unit/test_ring.cpp`: the concept cases ("the shipped models
  satisfy Ring", "Ring requires exactly size() and vertex(i)...", "the
  algorithms work through the concept...") with their test models
  `NoVertex`, `SizeIsInt`, `VertexByValue`, `ArrayTriangle`; every `edge`,
  `all_finite`, `bounding_box`, `signed_area` and `point_in_ring` case
  (`@44fa7f5:288-450` and `:564-750`, `:776-803`); the `PointRing`
  constructor cases that `IndexedRing` already pins ("rejects fewer than
  three", "rejects a stored closure", and `PointRing`'s half of "cannot be
  built from a temporary").
- In `tests/cpp/property/prop_ring_invariants.cpp`: every `point_in_ring`
  and `signed_area` case (`@44fa7f5:101-208`, `:278-357`) and both raster
  boundary cases (`:359-408`, `contains_strict` against `point_in_ring`).
- In `tests/cpp/unit/test_point.cpp`: every `Point3` case and both hash
  cases ("Point2 hash deduplicates...", "Point hashing is consistent with
  equality for signed zero"); the `Point3` lines of the traits and
  format-spec cases.
- In `test_bbox.cpp` the three `center` cases and the `center` line of the
  degenerate-box case (`@44fa7f5:131`); in `test_segment.cpp` "is_degenerate
  is exact" and the `is_degenerate` line at `:74`; in
  `test_predicates_orientation.cpp` the `is_left_turn` and `is_collinear`
  cases; in `test_raster.cpp` "contains_strict excludes the boundary",
  "boundary tolerance stays above rounding noise", and the
  `contains_strict` line at `:238`.

Kept with a new instrument:

- **`prop_cdt_invariants.cpp`'s domain test** ("no triangle lies in a hole
  or outside the domain", `@44fa7f5:369-391`) checks the live triangulator
  with `point_in_ring` as its oracle. `point_in_ring` and `PointInRing` move
  **verbatim** into a new `tests/cpp/support/point_in_ring.hpp`, namespace
  `terrain::test`, taking `const IndexedRing&`, and the test calls that.
  One case of the deleted unit suite moves with it, retargeted: "point_in_ring
  classifies a square" (`tests/cpp/unit/test_ring.cpp@44fa7f5:570-583`), so the oracle has a
  direct check of its own. `tests/cpp/support/mesh_queries.hpp@44fa7f5:21`'s comment then
  points at the support header.
- **`orientation<K>`'s cases** run over `IndexedRingCase` only:
  `ring_cases.hpp` drops `PointRingCase`, `RingCases` becomes the two
  `IndexedRing` cases, `ExactRingCases` the one. `describe(PointInRing)`
  moves to the support header or goes.
- **"orientation agrees with sign(signed_area)"**
  (`tests/cpp/property/prop_ring_invariants.cpp@44fa7f5:250-269`) keeps its orientation
  assertions and drops the area lines: the star generator builds
  counterclockwise rings by construction, which is the oracle.
- **`PointRing` policy cases with no `IndexedRing` twin** ("rejects a ring
  collapsed to a single point", "accepts the degeneracies that real data
  contains", "does not check finiteness", "a legal ring's cyclic rotation may
  be rejected as a closure") are retargeted to `IndexedRing` through
  `IndexedRingCase::Holder`. The closure and size checks they pin are
  `IndexedRing`'s too (`detail::check_ring_size`, `check_ring_closure`).

Prose in suites that becomes false: `tests/python/test_core_cdt.py@44fa7f5:967-968`
(points at `test_core.py`), `tests/python/test_features.py@44fa7f5:39-41` and
`:584-586` and `tests/python/test_io_vtk_legacy.py@44fa7f5:19-21` ("`__init__.py`
imports `_core`"), `tests/cpp/unit/test_refinement_refine.cpp@44fa7f5:8` and
`tests/cpp/unit/test_refinement_chunks_dynamic.cpp@44fa7f5:3-4` (`for_each_chunk`).
Master's T2 (`97eea35`, the Python audit's layering test) since pinned the
four increment files' citations of `tests/python/test_features.py@44fa7f5:583`, and moved the
`test_features.py` and `test_io_vtk_legacy.py` module checks into
`tests/python/test_layering.py`. The merge takes master's prose there, less
its reason "`__init__.py` imports `_core`", and moves `test_layering.py`'s
package-root row to layer L0 with no imports; the `fetch.http` ->
`tin_engine` exception goes with it (`python-audit.md`, F12).

**The red step: one behaviour test.** Deleting dead code changes no
behaviour a suite can observe, except one: today no `tin_engine` module
can be imported without the compiled extension, because `__init__.py`
imports `_core`. After this PR the package imports without it, and
`stats.py`'s docstring ("testable without the extension") becomes true.
`@tester` adds to `tests/python/test_stats.py` one test that runs, in a
subprocess with `sys.executable`, `import sys; sys.modules["tin_engine._core"]
= None; import tin_engine.stats`, and asserts exit status 0. Red at
`1035690` (`ModuleNotFoundError: import of tin_engine._core halted`), green
after I1. Probe, run for this design against a copy of the package with the
editable finder bypassed: with today's `__init__.py` the import fails with
that error; with `__init__.py` reduced to its `importlib.metadata` import it
succeeds. No test is written for "a removed name stays removed": in C++ an
absent name is a hard compile error, not something a `requires` clause can
probe, and a Python `hasattr` test would pin an absence nobody reads. The
red commit therefore holds the deletions and retargets above plus this one
failing test.

### 7.5 `@developer`: the production edits

Delete the table's ranges; turn `extreme_vertex` into an `inline` function
and `orientation` into `template <pred::GeometryKernel K>`, both on
`const IndexedRing&`; rewrite the comments of 7.1 in place; reduce
`src_python/tin_engine/__init__.py` to `installed_version`; edit `project_structure.md` in place.
No test file is touched. Gates: the C++ build under the compiler gate,
`ctest`, `pytest`, `mypy`, `ruff`, `tools/count_loc.py` (-291),
`tools/check_citations.py` (no `broken`; re-read the at-risk list as
quotations). The merge adds a `ROADMAP.md` row for PR B.

### 7.6 Citations pinned by this design

Unpinned citations into lines this PR deletes or moves were pinned to
`44fa7f5` in the commit that added this section (T2's lesson): `segment.hpp:68`
(`05-noder.md`, `kernel-sufficiency-audit.md`), `include/terrain/core/ring.hpp@44fa7f5:233`
(`05b-noder-driver.md`, twice), `_core.pyi:3` and `:437` (`05c-noder-wiring.md`,
`25-plain-output.md`), ten `bindings/core.cpp` citations
(`15f-edge-strip.md`, `23-basin-scale.md`, `27-node-sampling.md`,
`29-nve-reference-catchments.md`), `chunks.hpp:10` and
`tests/cpp/unit/test_refinement_chunks.cpp@44fa7f5:79` (`21-parallel-refine.md`),
`test_pslg_builder.cpp:961` and `prop_cdt_invariants.cpp:398`
(`h14-parallel-sanitizer-tests.md`). Left unpinned, because the edits above
keep their lines in place: `snap_grid.hpp:43`, `:63`, `:94`;
`project_structure.md:166`, `:190`, `:208`, `:359`, `:384`;
`.github/workflows/main.yaml:59`, `:306`;
`tests/cpp/CMakeLists.txt:3`, `:189`, `:199`. `@reviewer` re-reads each of
these as a quotation.

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

## Review

One line per step, PR B (branch `cpp-dead`):

- Red: `50bef89` (`@tester`) adds the `stats`-imports-without-`_core` test, deletes the tests of removed items, moves `point_in_ring` to test support; `90ba8be` makes `test_pslg`'s edge-vs-ring case compare against the ring's vertices, not `terrain::edge`.
- Green: `6530d30` (`@developer`) deletes the 20 ranges of 7.1; `tools/count_loc.py 4208932 6530d30` counts -291, as 7.3 planned.
- `ab79bfe` (`@tester`) deletes `test_refinement_chunks.cpp` and its CMake registration; ctest 843/843, pytest 5089 passed at this head.
- Conflict between 7.4 (the suite goes in the red commit) and 7.5 ("no test file is touched"): resolved by splitting it, `@developer` takes the suite off the TSan lists in place, `@tester` deletes the file next.
- Departure 1: `include/terrain/predicates/orientation.hpp@6530d30:9`'s comment now names `reversed`, not the removed `is_left_turn`.
- Departure 2: `_core.pyi`'s comment above the two `describe` overloads no longer points at the removed `cross` stubs.
- Departure 3: `bindings/core.cpp`'s module docstring no longer opens with "geometry primitives"; it names the triangulator and mesh surface.
- Docs (`@architect`): `project_structure.md` lines 14, 17, 52, 446, 501, 505, 506 rewritten one for one; `05-noder.md`'s two `segment.hpp` citations and two 7.6 names pinned to `44fa7f5`; `check_citations` exits 0.
