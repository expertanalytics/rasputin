# Increment 6 — the CDT viewer

Status: design settled; **6a and 6b-i have landed, and 6b-ii is on this
branch.** The contingency split under "Files and LOC" fired for 6a (the binding
surface, merged as `311459b`) and the 6b seam under "The firing condition" fired
too (6b-i, merged as `ab4681c`). This branch is `6b-ii`: `viz/style.py`,
`viz/svg.py`, `viz/fixtures.py`, `viz/__init__.py`'s re-exports and `cli.py`'s
`draw` command, against the two suites committed red in `5c9b13c`. With it, the
increment's stated requirement is met -- `rasputin draw catchment --out x.svg`
puts a picture in front of a person. Re-runnable:
`git log --oneline master..HEAD`, `ls src_python/tin_engine/viz/` and
`rasputin draw --list`.

**The second LOC gate fired as well, after the fact rather than in time to act
on it** -- see "A second gate" below for the measurement and for why the seam it
names could not be cut once the red suite existed.

**The number is authoring order, not ship order.** This increment is designed
sixth and **ships before `05b`** (`docs/increments/05b-noder-driver.md`, not yet
written). The ruling and its reasons are the first section below, because the
ordering is the substantive decision here and the rest follows from it.

Ships `bindings/core.cpp` additions, `src_python/tin_engine/_core.pyi`
additions, a new `src_python/tin_engine/viz/` package, one new Typer command,
and a rewrite of `ROADMAP.md`. Depends on increment 4 (`Pslg`, `triangulate`,
`IndexedMesh2`, `CdtOutcome`) and on nothing else. Depends on **no** increment-5
artefact, which is what makes the ordering free to choose.

Makes possible: **a person looking at a picture of a constrained Delaunay
triangulation and forming an opinion about it.** That is the whole requirement,
in the user's words — "so that I can use my intuition on the quality of our
work". Not coverage, not an assertion. Everything below is judged against
whether it gets a picture in front of a human that is worth having an opinion
about.

**Invariant-critical suite:** `tests/python/test_viz_scene.py`, and that suite
alone. Reasons, including the three suites that are deliberately *not*
invariant-critical and the golden-file test that is rejected outright, under
"What is worth testing".

## The ordering ruling: the viewer ships before 5b

Without the noder, input must already be noded — no crossing constraints, no
T-junctions, no coincident vertices — because `04-cdt.md` risk 1 makes
`DuplicatePointsFound` a hard failure and observes that real data is full of
coincident vertices. So a viewer landing now renders **hand-built synthetic
input** and not a real catchment. That is the cost. It is worth paying, for four
reasons, in descending weight:

1. **The thing the user wants intuition about already exists.** The quality
   questions a picture answers — how a sliver sits against a breakline, whether
   a hole is carved cleanly, what a near-collinear fan triangulates into, where
   the Delaunay property visibly stops at a constraint — are questions about the
   *triangulator*, which shipped in increment 4. The noder does not change what
   the triangulator produces; it widens the set of inputs the triangulator will
   accept. Waiting for 5b buys realism of *input*, not quality of *output*.
2. **5b is the increment that most needs an instrument, and building the
   instrument afterwards is building it too late.** `05-noder.md` records three
   oracles from increment 5 that agreed with the error they were written to
   catch: exact incidence where the producer used hot-pixel proximity, a
   "second" candidate generator that was the same broad phase, and the dedup's
   own records standing in for the input (guarantees 14 and 15; the pattern is
   generalised in `.claude/REQUIRED-READING.md` as the object-identity
   question). Every one of those is a failure a person spots in a second in a
   drawing — an unsplit T-junction is a line ending in the middle of another
   line. A renderer that exists before 5b makes 5b checkable by eye. A renderer
   that arrives after it does not.
3. **The adversarial cases most worth seeing are the synthetic ones.** A
   corner-touching hole, a fan of near-collinear constraints, a sliver at a
   4000:1 aspect ratio — these are shapes a person can hold in their head and
   compare the picture against. A 40 000-triangle river network is a picture
   whose correctness nobody can assess by looking; it is a picture you look at
   to see whether anything is *obviously* wrong, which is a weaker and later use.
4. **Cost and sequencing.** The viewer depends on nothing 5b produces, and 5b is
   estimated at ~480 non-comment lines of topology work (`05-noder.md`, "5b —
   the next PR"). Deferring a self-contained increment behind the largest
   remaining one delays every picture by a full round and de-risks nothing.

**The objection, stated rather than dismissed.** Intuition about a synthetic
sliver may not transfer to intuition about a river network. This is real and
partly unfixable. It is mitigated, not eliminated, by one fixture rule stated
under "The gallery": **fixtures are authored at realistic coordinate
magnitudes** — a UTM 33N-shaped domain with easting near 4.3e5 and northing near
6.9e6, not the unit square — because the robustness behaviour a person is being
asked to judge is magnitude-sensitive and a picture drawn in [0,1]² silently
exercises the easy case. Shape realism is cheap too: a catchment-like outline
with two lake holes and a braided breakline is thirty coordinates typed by hand.

The honest close-out is scheduled rather than assumed: **after 5b renders a real
catchment through this same renderer, re-ask the user whether the synthetic
gallery was predictive.** If it was not, that is a finding about the fixtures,
not about the renderer, and it is fixed by adding fixtures.

## The dependency ruling: hand-written SVG, and no new dependency

`CLAUDE.md` §2's stack is Pydantic V2, Typer, Shapely, PyProj, NumPy, tifffile.
**matplotlib is not in it, and this increment does not put it there.**

Ruled on the merits rather than inherited:

- What a plotting library would buy here is anti-aliased raster output, axes
  furniture, and interactive pan/zoom/pick. The first two we do not need — SVG
  is resolution-independent and this drawing has no axes — and the third is
  **explicitly out of scope** (see "Not in scope"). So the purchase is
  approximately nothing.
- What it would cost is a §2 amendment, a large transitive closure
  (`contourpy`, `kiwisolver`, `fonttools`, `pillow`, a backend selection
  problem, a font-availability problem), and a rendering pipeline whose output
  changes between versions — which is a poor foundation for a picture whose
  entire purpose is that a human trusts what they see.
- SVG costs one dependency-free module. It opens in any browser, scales
  losslessly, diffs as text in git, carries per-edge stroke classes and a legend
  natively, and can be pasted into a PR comment. It keeps the I/O boundary
  intact in the plainest possible way: C++ hands out a mesh as data, Python
  writes bytes.
- And if interactivity is ever genuinely wanted, the answer is still not
  matplotlib: it is the file we already emit, plus a little JavaScript, in the
  viewer the user already has open.

**Ruling: no new dependency.** NumPy is already in §2 and is used for the array
boundary; the renderer itself needs nothing beyond the standard library and
Pydantic V2 for its style model. `tools/check_prohibited_deps.py` is unaffected
— nothing here is on the prohibited list and nothing new is declared.

## The binding surface: what crosses, and what must not

This increment is the first time anything but a value type crosses the pybind11
firewall, and that is the part of it with permanent consequences. The rule from
`project_structure.md`'s `raster` section is absolute and is not relaxed here:
**the C++ core never opens a file, never sees a path, never links a codec, and
CRS never crosses into it.** A mesh crossing *out* is data, not a file. The
picture is written by Python, from data, to a path Python resolved.

### What crosses out

| Exposed as | From | Shape |
|---|---|---|
| `ChainRole` | `core/pslg.hpp` | `py::enum_`, values `Outer`, `Hole`, `Breakline` |
| `CdtStatus` | `cdt/result.hpp` | `py::enum_`, seven values, with `describe()` bound as `__doc__`-adjacent `describe(status) -> str` |
| `Pslg` | `core/pslg.hpp` | opaque, **not constructible from Python**, read-only accessors only |
| `IndexedMesh2` | `core/indexed_mesh.hpp` | opaque, read-only, three zero-copy arrays |
| `CdtOutcome` | `cdt/result.hpp` | read-only `status`, `message`, `mesh`, `ok()` |
| `build_pslg(...)` | `core/pslg_builder.hpp` | free function, declarative, returns a result object |
| `triangulate(pslg, delaunay=True)` | `cdt/triangulate.hpp` | free function, **GIL released** |
| `PslgError` | `core/pslg_builder.hpp` | `py::enum_`, eight values |

`PslgError` was missing from this table in the design as first written, and
`PslgDiagnostic.error` has to be *something*: 6a's suite requires only that four
independently broken chains give four distinct values, which an enum, an int and
a string all satisfy. **Ruled in 6a: an enum**, because `ChainRole` and
`CdtStatus` already establish that a C++ enumeration crosses this boundary as an
enumeration, a string would make every Python consumer compare magic text, and
an int would lose the names the diagnostic is read for. `PslgDiagnostic` and
`Chain` cross as frozen record classes for the same reason, and are named here
so the surface stays exhaustive.

`Pslg` accessors: `vertices` as a read-only `(N, 2)` `float64` array, and
`chains` as a list of frozen records `(begin, count, role, is_river)` plus
`chain_indices` as a read-only `(M,)` `uint32` array. `indices_of(c)` is bound
because the role map needs it and recomputing `begin:begin+count` in Python is
an off-by-one nobody should be given the opportunity to write twice.

`IndexedMesh2` accessors: `vertices` `(N, 2)` `float64`, `triangles` `(T, 3)`
`uint32`, `constrained_edges` `(T,)` `uint8`, plus `triangle_count` and
`empty`. All three arrays are **zero-copy views with `writeable = False`**, with
the owning mesh kept alive by the array's base object. That mirrors the raster
boundary contract's "the adapter sets the array read-only first" and it is the
one genuinely dangerous line in this increment — see risk 4.

### `build_pslg` is a free function, not a bound builder

`PslgBuilder` is **not** exposed. Three reasons, and the first is decisive:

1. `build()` is rvalue-ref-qualified (`PslgBuildResult build() &&`,
   `pslg_builder.hpp`). Binding it means either binding a consuming method on a
   Python object that remains reachable afterwards — a moved-from builder any
   Python name can still call — or wrapping it in a lambda that copies. Both are
   worse than not having the type.
2. A mutable stateful builder object on the Python side is exactly what §2 of
   the architect's brief calls state bleeding across a boundary. The Python
   layer should describe *what* PSLG it wants and receive one or a list of
   reasons why not.
3. It keeps the surface small enough to defend against accretion (risk 3).

So the signature is one declarative call:

```
build_pslg(vertices, chains) -> PslgBuildResult
```

where `vertices` is any `(N, 2)` float64-convertible array-like and `chains` is
a sequence of `(indices, role, is_river)`. The C++ side runs
`PslgBuilder::add_chain(span<const uint32_t>, role, is_river)` per entry and
`build<DefaultKernel>()`. The returned `PslgBuildResult` exposes `ok`, `pslg`
(`None` unless `ok`), and `diagnostics` as a list of frozen records
`(error, chain, vertex, message)` — the **whole** vector, because
`pslg_builder.hpp` is explicit that bulk-wrong input reported one failure at a
time turns one fix into N round trips, and a binding that returns only the first
diagnostic throws that design away at the boundary.

`PslgDiagnostic.vertex` is bound as-is, including its `kNoVertex` sentinel, and
the stub documents that it is **never** an out-of-range value — the header's
comment anticipates precisely a Python binding dereferencing it, so the stub
must carry that warning across rather than leave it in C++.

### What must not cross, ever

No path, no file object, no CRS, no `str` that means a filename. No raw pointer
and no `std::span` reaching Python unanchored. No mutable view of any core
buffer. No `detria` type, under any name. No `Segment2`, `IndexedRing`,
`PointRing`, `Chain` as a live C++ object, or `SnapGrid` — Python has Shapely
for geometry and nothing in this increment needs a C++ segment. No adjacency,
because `indexed_mesh.hpp` is explicit that it guarantees none and `mesh` owns
topology. No `PslgBuilder`. No kernel template parameter.

### GIL

`triangulate` wraps the call in `py::gil_scoped_release`. `bindings/core.cpp`'s
existing header comment already says this is required "the moment the
triangulation kernel is exposed"; this is that moment, and the comment is
updated to point at the call rather than at the future. Nothing else released:
the accessors are O(1) and the array views copy nothing.

## Where the code lives

```
bindings/core.cpp                       # +enums, +Pslg, +IndexedMesh2,
                                        #  +CdtOutcome, +build_pslg, +triangulate
src_python/tin_engine/
  _core.pyi                             # stubs for all of the above
  viz/
    __init__.py                         # re-exports render_svg, Scene, SvgStyle
    protocols.py                        # MeshLike, PslgLike -- typing.Protocol
    scene.py                            # build_scene(pslg, mesh) -> Scene  [the real logic]
    style.py                            # SvgStyle: Pydantic V2, frozen
    svg.py                              # (Scene, SvgStyle) -> str
    fixtures.py                         # the synthetic gallery, declarative
  cli.py                                # + the `draw` command
```

**`viz/` never imports `_core`.** It consumes `MeshLike` and `PslgLike`, two
`typing.Protocol`s in `protocols.py` describing exactly the read-only array
attributes listed above. This is not decoration: it is what lets the entire
renderer be unit-tested against a twelve-line fake mesh without a compiled
extension in the process, and it is what lets the two halves of the contingency
split be built in either order. The only module that touches `_core` is
`cli.py`, which is the composition root — the same shape as
`project_structure.md`'s rule that exactly one Python module constructs a core
raster.

`scene.py` and `svg.py` are **pure functions over immutable data**. No file is
written below `cli.py`; `render_svg` returns a `str`. That is what makes the
renderer trivially wrappable in `asyncio.to_thread` by a future GUI or API
worker without this increment carrying an async layer that has no I/O in it —
see "Async: ruled out, deliberately".

### The command

One line to a picture:

```
rasputin draw sliver-fan --out /tmp/sliver.svg
```

- `rasputin draw --list` prints the gallery: name, one-line description, what
  it is there to show.
- `--out PATH` is optional; with no `--out` the command writes to a temp file
  and prints the path, because the common case is "show me this now".
- `--labels` draws vertex indices; refused above `--label-limit` (default 500
  vertices) with a message naming the count, rather than emitting an
  unreadable file.
- `--no-delaunay` forwards `CdtOptions::delaunay = false`, so the user can see
  *both* triangulations of the same input and form an opinion about what the
  Delaunay property is buying. This is one bound bool and it is the cheapest
  intuition in the increment.
- `--title TEXT` is a free string drawn in the header band. A CRS name, if the
  user knows one, goes here — as text, on the Python side, having never been
  anywhere near C++.
- Path handling per the Python skill: `pathlib.Path.resolve()`, reject a path
  that resolves outside an explicitly provided parent, no symlink following.

There is no `--geojson`. Real data ingestion is out of scope; see "Not in
scope".

## What the picture shows

Drawn in this order, because later strokes must win where they overlap:

1. **Page background**, a light neutral that is *not* white. A hole in the mesh
   is simply the absence of triangles, so it reads as a hole only if the page
   reads as a page. This is the one stylistic decision that is load-bearing for
   comprehension rather than taste.
2. **Triangle fill**, a single pale tone. No data-driven colour: there is no
   datum yet — elevation is out of scope — and a gradient invented to look busy
   is a gradient someone will later read as meaning something.
3. **Unconstrained edges**, thin, light. Deduplicated: each interior edge is
   shared by two triangles and must be emitted once, or every interior edge is
   drawn at double weight and the picture lies about density.
4. **Constrained edges**, thicker and saturated, drawn on top, **coloured by the
   role of the input chain they came from** — outer boundary, hole boundary,
   breakline, river breakline. Four visually distinct strokes.
5. **Vertices**, small dots, off by default; **indices**, off by default.
6. **Legend** (the five stroke classes), **scale bar** in world units, and a
   **header band** carrying status, triangle count, vertex count,
   constrained-edge count and the bbox.

### Role colouring requires a second, independent derivation — on purpose

`IndexedMesh2`'s mask says *constrained*, not *which role*. Role lives in the
input `Pslg`'s chains. The scene therefore builds a map from a sorted vertex
index pair to a `ChainRole`, walking `pslg.chains()` and `pslg.chain_indices()`,
and joins it against the mesh's mask. This is sound because `triangulate.hpp`'s
semantic obligation 1 and `indexed_mesh.hpp` guarantee 4 both pin that the mesh
vertex array **begins with** the PSLG's, element-wise and in order, so index `k`
means the same point on both sides; and because `04-cdt.md`'s central finding is
that detria splits nothing, so each input constraint edge is exactly one output
edge and the join is 1:1.

#### Ring closure is the caller's policy, not the scene's — `closed_roles`

**Defect found and closed during 6b-i.** The paragraph above asks the scene to
colour ring edges *including their closure* — a `Pslg` never stores a ring's
closing edge — while `viz/protocols.py` types `ChainLike.role` as `object`
precisely so that `viz/` cannot name `_core.ChainRole`. Those two requirements
are jointly unsatisfiable inside `scene.py`: it cannot tell a ring from a
breakline, so it either closes every chain (inventing a constraint on every
breakline) or closes none (making every ring's closing edge a false
`MASKED_EDGE_WITHOUT_CHAIN`).

Resolved by moving the policy out rather than the type in:
`build_scene(..., closed_roles=())` takes the roles whose chains are rings, and
`cli.py` — the single composition root, which already imports `_core` and knows
the enum — supplies `(ChainRole.Outer, ChainRole.Hole)`. Membership is tested
with `==` only, which `object` provides, so the opacity survives intact. A chain
of two or fewer vertices is never closed, because fewer than three vertices is
not a ring — and, for the one-vertex case, because the closure would be
`_key(a, a)`, a self-loop violating `SceneEdge`'s `a < b` ordering. Not because
of double emission: an earlier revision of this line said "closing one would
emit its single edge twice", which the implementation contradicts — the join
accumulates into a dict keyed on `_key(a, b)`, so a repeated key only ORs
`is_river`.
`tests/python/test_viz_scene.py::TestRoleJoin::test_a_single_vertex_ring_is_never_closed_into_a_self_loop`
is what defends the guard.

`tests/python/test_viz_scene.py::TestRoleJoin::test_closure_policy_is_the_callers`
pins the defaulted behaviour: with no `closed_roles`, the closing edge is
reported as a disagreement rather than guessed at.

Both halves of that join can therefore disagree, and **the renderer draws the
disagreement instead of reconciling it**:

- a mesh edge whose mask bit is set but which matches no input chain pair is
  drawn in an alarm colour and counted in the header band;
- an input chain edge that appears in no triangle's mask is drawn as a dashed
  alarm stroke along its own endpoints, likewise counted.

Either count being non-zero is a finding about the backend's mask, visible at a
glance. This is the object-identity discipline of `.claude/REQUIRED-READING.md`
applied to a picture: the two sources are genuinely different objects — the
backend's mask and the input's chains — so their agreement means something.

### The mask convention is a rotation away from being silently wrong

`indexed_mesh.hpp` defines bit `e` as the edge `(v[e], v[(e+1)%3])` and says in
so many words that this is **not** CGAL's "edge `e` is opposite vertex `e`"
convention, which is a rotation of it. A renderer that rotates it draws a
perfectly plausible picture with the constraint strokes on the wrong edges — a
wrong picture that does not *look* wrong. The scene builder must read edge
endpoints through `IndexedMesh2::edge(t, e)`-equivalent arithmetic and nothing
else, and this is the single most important thing the invariant-critical suite
pins.

### Degeneracy and emptiness are the loudest thing on the page

`04-cdt.md` risk 5: `Ok` with an empty mesh is a *successful* backend call, and
in a release build with the assert compiled out, one test is the only guard. A
viewer is the natural second guard, and it only works if a failure is
conspicuous rather than blank.

So the renderer **never** emits a blank page:

- **Non-`Ok` outcome** — no mesh is drawn. The input PSLG is drawn alone, in its
  role colours, over a tinted background, with the status name and
  `describe(status)` set large in the header band and the backend's `message`
  beneath it. The user sees exactly the geometry they submitted and exactly what
  the engine said about it — and if a non-`Ok` run arrives carrying a full mesh,
  which `testing.md`'s `cdt` catalog says no outcome ever does, the status still
  wins over a mesh that cannot exist.
- **`Ok` with zero triangles** — same treatment, plus the header band reading
  `Ok BUT EMPTY` in the alarm colour. This is the specific silent failure risk 5
  names, and after this increment it is the most visually obvious state the tool
  can produce.
- **Zero-area bounding box** (all vertices coincident or collinear) — the
  viewport transform would divide by zero. The bbox is padded to a minimum
  extent and the header band says so, rather than the picture being an
  accidental point.
- **Non-finite coordinates** cannot reach here — the PSLG validator's stage 3
  rejects them — but finiteness is checked anyway, because the cost is one line
  and the failure mode without it is an SVG that renders as nothing in silence.

  **The check moved from the viewport transform (6b-ii) to `build_scene`
  (6b-i), during 6b-i.** The reason is that the bbox is computed before the
  transform ever runs, and a single NaN vertex makes `min`/`max` NaN, which
  makes the whole box NaN: by the time the transform could assert, the scene
  it would assert about is already meaningless, and every downstream comparison
  against a NaN bound is silently false. `build_scene` raises `ValueError` on
  any non-finite coordinate of the array it is about to bound — the mesh's
  vertices when a mesh is drawn, the PSLG's otherwise. On *any*, not on the
  first: `np.isfinite(vertices).all()` inspects the whole array and the message
  names no coordinate.
  `TestBoundingBox::test_a_non_finite_coordinate_is_refused` pins it over
  `nan`, `inf` and `-inf`. 6b-ii's transform therefore inherits a finite bbox as
  a precondition and does not re-check it.

**No drawable mesh means no findings.** Whenever the scene draws the input PSLG
alone — non-`Ok`, or `Ok` with zero triangles, or no mesh at all — there is no
constrained-edge mask for the chains to disagree with, so `findings` is empty
and every chain edge keeps its role. The alternative would be an alarm per chain
edge on a picture whose real alarm is already the status; it is also what makes
"Either count being non-zero is a finding about the backend's mask" true, since
a run with no mask contributes no counts.
`TestClassification::test_a_scene_without_a_mesh_reports_no_findings` pins it.

### Y-axis flip

World y increases upward; SVG y increases downward. The viewport transform flips
it. If it does not, every picture is mirrored, every winding reads backwards,
and the user's intuition is being trained on a reflection. It is one sign and it
is worth a named test.

## Async: ruled out, deliberately

`CLAUDE.md` §2 and the python skill both push toward async orchestration, and
this increment declines it. `render_svg` is CPU-bound pure computation over
in-memory arrays; the only I/O is one `Path.write_text` at the CLI edge. An
`async def` wrapper here would contain no `await` that does anything, and
ceremony in the one place a reader looks for a worked example is worse than
absence.

What the increment *does* owe async-readiness is discharged where it actually
matters: `triangulate` releases the GIL, and `render_svg` is a pure function of
its arguments. Those two facts are exactly what a future GUI backend needs to
write `await asyncio.to_thread(render_svg, mesh, pslg, style)` without any
change here. Stated so that nobody later "fixes" the absence.

## The gallery

Declarative data in `fixtures.py` — a name, a description, a vertex array, a
list of `(indices, role, is_river)`. No procedural generation, because a fixture
whose coordinates are computed is a fixture nobody can check by reading. All
authored at realistic coordinate magnitudes (see the ordering ruling).

Each exists to answer one question a person can ask of the picture:

| Name | Shows |
|---|---|
| `catchment` | outer ring with two lake holes and a braided breakline — the shape the project is actually for |
| `sliver-fan` | a fan of near-collinear constraints; what the triangulator does with extreme aspect ratios |
| `corner-hole` | a hole touching the outer ring at exactly one vertex — `InvalidTopology`'s neighbour, and a topology a person should look at |
| `hole-in-hole` | a hole nested directly inside a second hole. **Drawn as a third failure presentation, not as a mesh** -- the backend answers `InvalidTopology`, "A hole was directly inside another hole", which is precisely "what in-domain means" made visible. `test_viz_svg.py::TestGallery::test_hole_in_hole_really_nests_twice` requires two chains of role `hole` with one strictly inside the other, so every fixture satisfying the suite is a fixture the backend refuses; the row as first written expected a mesh |
| `breakline-chain` | an open breakline crossing the interior; where the Delaunay property visibly stops |
| `river` | the same, with `is_river` set, so the `is_river` stroke is exercised before 5b depends on it |
| `not-noded` | two crossing constraints — a **deliberate failure fixture**, rendering the non-`Ok` presentation, so that presentation is seen rather than assumed |
| `degenerate` | an all-collinear point set. **Refused before `triangulate` is ever reached, so this is *not* the `DegenerateGeometry` presentation** -- `build_pslg` rejects the ring with `PslgError.DegenerateRing`, "chain 0 is declared Outer but every vertex is collinear", and returns no `Pslg` at all. The design's intent survives: the input is drawn alone with the engine's own words in the header band. The words are a `PslgDiagnostic`'s. Found in 6b-ii's red step; see "The `PslgLike` is the fixture itself" below, which is the decision this forced |

The last two matter as much as the first six: a failure presentation nobody has
looked at is a failure presentation that is wrong.

### The `PslgLike` is the fixture itself -- and why the roles are strings

**Defect found and closed during 6b-ii**, in the same shape as 6b-i's
`closed_roles` finding: the design named a presentation without checking which
object would still exist when it fired.

`build_scene` takes a `PslgLike`. For a fixture the validator *rejects* there is
no `Pslg` -- `build_pslg` returns `ok == False` and `pslg is None` -- so on the
`degenerate` fixture there is nothing to hand it, and "the input PSLG is drawn
alone" has no input PSLG. Two ways out: give `build_scene` a nullable `Pslg` and
a separate copy of the vertices and chains, or make the fixture itself satisfy
`PslgLike`. **Ruled: the fixture is the `PslgLike`**, because the protocol exists
precisely so that anything with those four members can be drawn, and because the
alternative re-derives the chain structure at the composition root, which is one
more place for the role join to go wrong.

That forces the roles. `viz/` never imports `_core`, so `fixtures.py` cannot name
`ChainRole`; it authors roles as the strings `"outer"`, `"hole"` and
`"breakline"`, and `cli.py` -- which already knows the enum -- maps them for
`build_pslg` and passes `closed_roles=("outer", "hole")`. `ChainLike.role` is
typed `object` and the scene compares roles with `==` only, so both spellings
work unchanged; the renderer names a stroke class from either with
`getattr(role, "name", role)`. `test_viz_svg.py::TestStrokeClasses::test_a_string_role_names_its_own_stroke_class`
is the half of that only the string exercises, and
`test_cli_draw.py::TestFailurePresentation` is what makes the whole arrangement
necessary rather than merely possible.

## Files and LOC

The unit is **non-comment production lines**, per `CLAUDE.md` §2. Note that
pybind11 docstrings are string arguments and therefore *count*, which is why the
binding estimate is larger than its apparent complexity; the existing
`bindings/core.cpp` is the calibration.

| File | Contents | Est. LOC |
|---|---|---|
| `bindings/core.cpp` (additions) | two enums, `Pslg`, `IndexedMesh2`, `CdtOutcome`, `PslgBuildResult`, `build_pslg`, `triangulate`, array views + keep-alive | ~150 |
| `src_python/tin_engine/_core.pyi` (additions) | stubs for all of the above | ~60 |
| `src_python/tin_engine/viz/protocols.py` | `MeshLike`, `PslgLike`, `ChainLike` | ~25 |
| `src_python/tin_engine/viz/style.py` | `SvgStyle`, frozen Pydantic V2 | ~45 (**realised 33**) |
| `src_python/tin_engine/viz/scene.py` | edge dedup, mask decode, role join, disagreement findings, bbox, `Scene` | ~90 |
| `src_python/tin_engine/viz/svg.py` | viewport transform, element emission, legend, scale bar, header band | ~110 (**realised 304**) |
| `src_python/tin_engine/viz/fixtures.py` | the eight fixtures | ~80 (**realised 175**) |
| `src_python/tin_engine/viz/__init__.py` | re-exports | ~5 |
| `src_python/tin_engine/cli.py` (additions) | `draw`, `--list`, path validation | ~45 (**realised 139**) |

**~610 production LOC.** Under `CLAUDE.md` §2's 700, but the margin is thin and
the binding half is the half that overruns. Non-production changes on top:
`CMakeLists.txt` gains `terrain_cdt` **and `terrain_predicates`** on the `_core`
link line (the bindings target currently links headers only) -- the second
because `build_pslg` validates under `DefaultKernel`, which is the same pair,
for the same two reasons, that `tests/cpp/CMakeLists.txt`'s
`add_terrain_cdt_test` links. Omitting it is not caught by the C++ build: a
pybind11 MODULE links on macOS with undefined symbols allowed and fails at
`import`; `testing.md` gains a `viz` section marked
`[planned]` until merge, then `[live]`; `project_structure.md`'s directory
listing, `bindings/core.cpp` section and Python API surface section gain `viz/`
and the new binding surface.

**Pre-declared contingency split, dependency-ordered** — declared here rather
than discovered mid-PR, because the seam is already natural:

- **6a — the binding surface.** `bindings/core.cpp`, `_core.pyi`,
  `CMakeLists.txt`, plus `viz/protocols.py`. ~235 LOC. Independently valuable:
  after 6a a person can triangulate from a Python REPL and look at the arrays.
- **6b — the renderer.** `viz/{scene,svg,style,fixtures,__init__}.py` and the
  CLI command. ~375 LOC. Buildable and fully testable against the protocols
  without 6a, which is the point of `protocols.py` being in 6a's half.

**6a as shipped: 467 non-comment production lines against the ~235 estimated**
-- `bindings/core.cpp` +283 (est. ~150), `_core.pyi` +125 (est. ~60),
`viz/protocols.py` 52 (est. ~25), `viz/__init__.py` 7 (not estimated). Measured
on the tree this paragraph ships in, by counting added, non-comment, non-blank
lines:

```sh
for f in bindings/core.cpp src_python/tin_engine/_core.pyi \
         src_python/tin_engine/viz/protocols.py \
         src_python/tin_engine/viz/__init__.py; do
  git diff master...HEAD -- "$f" | grep '^+' | grep -v '^+++' | sed 's/^+//' \
    | grep -vcE '^\s*(//|#|\*|/\*|\*/|$)'
done
```

Record the figure against the tree it ships in, not against an earlier commit.
Three revisions failed that and went stale within one commit each: 441/264, then
437/260 at `bc53100`, then 457/278 at `1424300` once the `PslgError` value
strings and the shape formatter landed -- and this revision's own +10 is the
`chains`/`diagnostics` copy documentation. The command is the definition, so a
figure citing a commit older than the paragraph is a claim the reader's own run
will contradict. Under `CLAUDE.md` §2's 700 at all four.

The overrun is entirely docstring and stub body: risk 6 named pybind docstrings
as the pressure and named the right one, but the estimate was made against a
calibration file with far shorter ones. Two consequences, and they are
different in kind:

- **The halves cannot be recombined.** 437 + the ~375 estimated for the
  renderer is 812. 6a and 6b are two PRs, not one, and that is now arithmetic
  rather than caution.
- **6b's own estimate is suspect, but not by the same factor.** See the seam
  ruling below.

### 6b's seam: pre-declared, armed, with a firing condition

**Ruled: yes, 6b carries a pre-declared seam, and it cuts at `Scene`.**

- **6b-i — the geometry mapping.** `viz/scene.py` and
  `tests/python/test_viz_scene.py`, which is this increment's sole
  invariant-critical suite and therefore carries the mutation round. Est. ~90
  production lines. It consumes `MeshLike` and `PslgLike`, which 6a already
  shipped, and imports nothing else; it is reviewable end to end against a
  hand-built fake with no compiled extension in the process.
- **6b-ii — the renderer and the command.** `viz/svg.py`, `viz/style.py`,
  `viz/fixtures.py`, `viz/__init__.py`'s re-exports, `cli.py`'s `draw`, and the
  `test_viz_svg.py` / `test_cli_draw.py` suites. Est. ~285.

**The cut is dependency-ordered and cannot be reversed.** `svg.py` is
`(Scene, SvgStyle) -> str` and `Scene` is defined in `scene.py`, so 6b-i ships
first. `style.py` goes with 6b-ii, not 6b-i: `SvgStyle` is the renderer's input
and the scene builder never sees it, so pairing them keeps 6b-i purely the
mesh-and-PSLG-to-geometry mapping. This is the one place the split is *not*
symmetric with 6a's, where `protocols.py` was deliberately pushed into the
earlier half so either half could be built first.

The honest cost: **6b-i produces no picture.** It ends with a `Scene` a person
can inspect in a REPL and a mutation round, and the increment's stated
requirement -- a human looking at a drawing -- is met only at 6b-ii. That is
acceptable for a contingency that may not fire; it would not be acceptable as a
default, which is why the seam is armed rather than fired.

**Why not simply apply 6a's realised factor.** 437/235 is 1.86, and 1.86 x 375
is 697 -- at the ceiling, which is the argument for a seam. But that factor is
not transferable verbatim, and the counter-evidence is in this repo: increment
3's `include/terrain/core/pslg_builder.hpp` was estimated at 255 and came in at
**249 non-comment** (`grep -vcE '^\s*(//|\*|/\*|\*/|$)'
include/terrain/core/pslg_builder.hpp`; 457 raw), a factor of 0.98 on a file
that is nothing but algorithm. The 6a blow-up is concentrated in
declaration-surface files -- pybind docstrings, `.pyi` stub bodies,
`Protocol` members, each of which is one line of meaning per several lines of
text. `scene.py` and `svg.py` are algorithm and should behave like increment 3;
`fixtures.py` is declarative data and is the 6b file most likely to behave like
6a. So 6b's plausible range is ~400-600, not a confident ~697.

**The firing condition, because increment 3 named a seam and it went unused.**
Increment 3 said "no split" *and* named where the split would cut, and the seam
was never used because nobody re-measured -- a note to a future reader is not a
contingency. So the measurement here has a moment, an owner, a command and a
threshold:

> After `viz/scene.py` is green and before `viz/svg.py` is started,
> `@developer` runs the count command above over `src_python/tin_engine/viz/`.
> **If `scene.py` alone exceeds 150 non-comment lines** (est. 90; factor 1.67),
> 6b-ii is split into its own PR without further argument.

**Measured on the tree this paragraph ships in: `scene.py` is 193 non-comment
lines by the command, 194 in truth. The gate fires; 6b-ii is its own PR.**
Re-runnable:
`grep -vcE '^\s*(//|#|\*|/\*|\*/|$)' src_python/tin_engine/viz/scene.py`.
The one-line gap is the command's `\*` alternative — meant for a C block-comment
continuation — swallowing `build_scene`'s keyword-only marker, a bare `    *,`.
Noted rather than corrected. The command is not the definition of the unit —
`CLAUDE.md` §2 is, exclusively, and the unit is non-comment lines. The grep is
the **instrument** every figure in this document was taken with, and it
approximates that unit imperfectly, as this one line demonstrates. It is kept
because a Python-only variant would make the numbers incomparable to 6a's, and
because a one-line undercount at this size cannot flip a gate on its own. The
blind spot is exactly one line in `scene.py` and zero lines in all three 6a
Python files, so no published figure is affected. Re-runnable:
`grep -nE '^\s*\*' src_python/tin_engine/viz/*.py src_python/tin_engine/_core.pyi`
prints exactly `scene.py:203:    *,`.
Factor 2.16 on the ~90 estimate. The composition is the interesting part, and it
contradicts the "algorithm behaves like increment 3" prediction only in the
accounting: **68 of the 194 true lines are docstrings** and 126 are
executable — against `@tester`'s undocumented throwaway prototype at 130. The
algorithm came in *at* the prototype's size; what the count command sees on top
of it is documentation, since `grep -v '#'` excludes `#` comments and counts
docstring lines. So the gate fired on declaration-and-documentation surface
again, exactly as 6a did, and the split it forces is still the right call —
`svg.py`, `style.py` and `fixtures.py` would be documented to the same standard
and 500 is not a generous remaining budget. Re-runnable for the 68/126 split:
```sh
python3 - <<'PY'
import ast, pathlib, re
src = pathlib.Path("src_python/tin_engine/viz/scene.py").read_text()
doc = {n for e in ast.walk(ast.parse(src))
       if isinstance(e, ast.Expr) and isinstance(e.value, ast.Constant)
       and isinstance(e.value.value, str)
       for n in range(e.value.lineno, e.value.end_lineno + 1)}
counted = [i for i, l in enumerate(src.splitlines(), 1) if not re.match(r"^\s*(#|$)", l)]
print(len(counted), len([i for i in counted if i in doc]))
PY
```
(prints `194 68`; the 194 is the true count, per the note above, and
126 = 194 - 68 is the executable remainder.)
This is the measured evidence for the open question of whether `CLAUDE.md` §2's
unit should exclude docstrings; it is not settled here, and the gate is applied
as written.

150 rather than a number near 700 because the measurement has to happen while a
split is still cheap: `scene.py` is the only 6b module whose realised density
can be known before the largest module is written, and a split decided at 690
is a split decided too late. If `scene.py` lands at or under 150, the seam
stands down and 6b ships as one PR.

**A second gate, because the first samples the wrong file.** The trigger above
measures `scene.py`, which this document argues is *algorithm* and therefore
the 6b module least likely to inflate -- while naming `fixtures.py` as the one
most likely to behave like 6a, and `fixtures.py` sits in 6b-ii, the half the
first gate does not cover. A gate that samples the low-risk file and stands
down on its behalf can pass and still overrun for exactly the reason predicted.
So: **after `fixtures.py` is green, re-run the count over `viz/`. If the 6b
total exceeds 500 non-comment lines, 6b-ii splits again** -- at 500 rather than
700 because `cli.py` and the two ordinary suites still follow, and the point of
a pre-declared seam is that it fires before the ceiling, not at it.

**Measured on the tree this paragraph ships in: the 6b total is 709 by the
command, and the gate fires by 209 lines.** Per file, by
`grep -vcE '^\s*(//|#|\*|/\*|\*/|$)'`: `scene.py` 193 (6b-i's, unchanged),
`svg.py` 304, `fixtures.py` 175, `style.py` 33, plus 4 added lines in
`__init__.py`. The instrument's `\*` blind spot swallows four more real lines
(three in `svg.py`, one in `cli.py`), so the true figure is 713; it is recorded
the same way 6b-i's was, with the command's number first, because the command is
what a reader will re-run.

**The seam it names could not be cut, and that is the finding rather than an
excuse.** The gate's own wording places the measurement "after `fixtures.py` is
green" -- which is after `@tester` has committed `test_viz_svg.py` red, and that
single file covers `style.py`, `svg.py` **and** `fixtures.py` together. Cutting
6b-ii in two therefore means splitting a committed test file, and
`docs/increments/README.md` step 3 forbids the green step from touching one. So
the gate as written can only fire once the split it authorises has become
unavailable. **The fix is to the protocol, not to this increment: a LOC gate has
to sit before the red step, not after it**, because the red step is what fixes
the shape of the PR. Any future pre-declared seam should be measured against the
*design's* file list at the moment `@tester` is briefed.

What the ceiling in `CLAUDE.md` §2 actually bounds -- one pull request -- is not
breached: this PR is **655** production lines by the command (659 true), 45
under the 700. The 709 is the cumulative 6b figure across two PRs, which is what
the gate asked for and which no rule bounds. Composition, for the open question
about whether §2's unit should exclude docstrings: of `svg.py`'s 307 true lines,
107 are docstring and 200 executable; `fixtures.py` is 32 and 143, and its 143
is almost entirely typed coordinates; `style.py` is 23 and 10. Declaration and
documentation surface is where the lines went for the third time in this
increment, exactly as 6a and 6b-i recorded.

## What is worth testing

A renderer's failure mode is a wrong-looking picture, and a person catches that
instantly. Testing effort goes where the failure mode is a *right*-looking wrong
picture.

**Invariant-critical, mutation round applies — `tests/python/test_viz_scene.py`
only.** The mesh-and-PSLG-to-geometry mapping:

- the mask bit -> edge endpoint convention, including a mutant that rotates it
  to CGAL's (this is the mutant the suite exists for);
- interior edge dedup, including that a shared edge appears once and that a
  constrained edge shared by two triangles is not double-classified;
- the role join: an edge's role comes from the chain that contains it, an edge
  in two chains resolves by a stated precedence, and both disagreement findings
  fire when constructed;
- the degenerate and empty classifications, including `Ok`-with-zero-triangles;
- bbox including the zero-extent padding.

**Not invariant-critical, ordinary tests:**

- `tests/python/test_viz_svg.py` — that the output parses as XML, that element
  counts match the scene's primitive counts, that the viewport transform flips y
  and preserves aspect ratio, and that a non-`Ok` scene emits the status text.
  **No assertion on any style string.** The stylesheet gets no test at all.
- `tests/python/test_core_cdt.py` — array shapes, dtypes, `writeable is False`,
  that an array outlives the mesh object being dropped (risk 4), enum round
  trips, that a bad PSLG returns the full diagnostics list rather than raising,
  and that no `Ok` outcome carries an empty mesh.
- `tests/python/test_cli_draw.py` — `--list` names every fixture, `--out` writes
  a file, `--labels` above the limit exits non-zero with the count in the
  message, a path outside the permitted parent is refused.

**Rejected outright: a golden-file SVG comparison.** It would pin the
stylesheet, fail on every cosmetic improvement, and detect only the class of
defect a human sees instantly. Written down because it is the first thing
someone will propose.

## Not in scope

Elevation, z coordinates, hillshade, any 3D. Refinement and its visualisation.
Real data ingestion — no GeoJSON reader, no GeoTIFF reader, no `--geojson`
flag. CRS, in any form other than a `--title` string Python never interprets.
Interactivity: no pan, no zoom, no picking, no JavaScript. Raster output: no
PNG, no `cairosvg`. Animation. Mesh editing. A `NodedPslg` before/after view —
that belongs to 5b's own PR, which will consume this renderer rather than change
it.

## Risks

1. **Synthetic intuition may not transfer to real hydrology.** The ordering
   ruling's own cost. Mitigated by realistic coordinate magnitudes and
   catchment-shaped fixtures, not eliminated. Discharged by a scheduled
   re-ask of the user after 5b renders a real catchment.
2. **The role join is a second derivation of "which edges are constrained".**
   If it and the backend's mask were derived from the same source, their
   agreement would confirm a bug rather than test it. They are not — one is the
   input's chains, the other the backend's output mask — and the renderer
   *shows* disagreement rather than reconciling it. Stated so nobody later
   "simplifies" the join by reading the role off the mask.
3. **Binding-surface gravity.** Once `Pslg` and `IndexedMesh2` are reachable
   from Python, every later increment is tempted to add one more accessor, and
   the firewall erodes by accretion rather than by decision. The surface above
   is exhaustive; an addition is a design change and needs a reason in an
   increment file, not a line in a PR.
4. **Zero-copy lifetime.** A `py::array_t` over an `IndexedMesh2`'s vectors
   without a correct base object is a use-after-free whose symptom is a
   plausible-looking wrong picture. It is the one genuinely dangerous line here.
   **What guards it, as shipped in 6a:** every array-returning accessor goes
   through one helper, `readonly_view` in `bindings/core.cpp`, which takes the
   owning `py::object` as the array's base -- six surfaces in total
   (`Pslg.vertices`, `Pslg.chain_indices`, `Pslg.indices_of`,
   `IndexedMesh2.vertices`, `.triangles`, `.constrained_edges`). `@reviewer`
   read the mechanism on all six, and all six have now been falsified by
   substituting a decoy owner and confirming the suite goes red: the three
   `IndexedMesh2` arrays in `94f94e2`, the three `Pslg` arrays in `16826a4`.
   The second round also measured the gap rather than asserting it -- with the
   `Pslg` keep-alive removed, the new test failed and **all 77 others passed**,
   so the first round could not have caught it.
   The assertions are `TestZeroCopyLifetime` in `tests/python/test_core_cdt.py`:
   one test drops every reference to the mesh and the outcome, churns the
   allocator, and re-reads; a second asserts the base object and its refcount
   directly, so a failure says "no keep-alive" rather than "the numbers
   changed". That pair runs in the `python` job of
   `.github/workflows/main.yaml` on all three interpreter legs.

   **What does not guard it: asan.** No CI job compiles this increment's array
   code under a sanitizer. The `sanitizers` job
   (`.github/workflows/main.yaml:48-69`) configures `-DCMAKE_BUILD_TYPE=Debug`
   with `-fsanitize=address,undefined` and sets no other option, and
   `RASPUTIN_BUILD_PYTHON` defaults `OFF` (`CMakeLists.txt:73`), so
   `bindings/core.cpp` is never in that build at all. The `python` job builds
   the extension -- via `pip install -e`, which configures
   `RASPUTIN_BUILD_PYTHON=ON` (`pyproject.toml:77`) -- but with no sanitizer
   flags. Re-runnable in two greps:
   `grep -n RASPUTIN_BUILD_PYTHON .github/workflows/main.yaml` returns nothing,
   and `grep -n 'option(RASPUTIN_BUILD_PYTHON' CMakeLists.txt` shows the `OFF`.
   An earlier revision of this line claimed "and by the asan Debug job"; it was
   never true, and it is exactly the "X is verified by Y" where Y evaluates a
   different object that `.claude/REQUIRED-READING.md` rules on -- named here
   rather than deleted, because this increment is the one that cites that rule.

   Closing the gap means enabling `RASPUTIN_BUILD_PYTHON` on the sanitizer
   runner, which needs pybind11 and a Python dev environment installed there.
   That is its own change with its own justification and is **not** part of this
   increment; the `RASPUTIN_SANITIZER` option already listed as *planned* in
   `project_structure.md`'s build section is the natural place for it.
5. **`--labels` on a large mesh** produces an enormous unreadable file.
   Guarded by the refusal above `--label-limit`, which is a hard error and not a
   warning.
6. **The LOC margin is thin** (~610 of 700) and pybind docstrings count against
   it. Hence the pre-declared split rather than a discovered one. **Realised:**
   this risk fired. The split was used, 6a overran its estimate by 86%, and the
   cause was the one named here. See "6b's seam" for what follows from it.
7. **Scope creep toward a plotting library.** The first feature request that
   cannot be served by static SVG — hover, picking, a colour ramp over elevation
   — will be argued as a reason to add matplotlib. The dependency ruling above
   is the answer, and the escape hatch is the browser, not a plotting stack.

## Documentation this PR fixes

Per `docs/increments/README.md`, these are fixed in this PR or not recorded.
There is no ledger.

**This section is split by owing PR, because two different obligations were
being conflated.** `docs/increments/README.md`'s rule closes over *defects found
during* an increment -- pre-existing rot the increment happened to walk past.
Deferring one of those to the half of the increment that creates the thing it
documents is consistent with the rule. But a falsehood an increment *creates*
is not a found defect; nothing licenses shipping a document this branch made
untrue, and a deferred entry of that kind reads as work 6a skipped. So:

- **Owed by 6a**, and fixed in this PR: the two `project_structure.md`
  statements that `94f94e2` falsified.
- **Owed by 6b**: the `ROADMAP.md` rewrite and the `testing.md` `viz` section.
  Both name 6b's modules and `tests/python/test_viz_scene.py`, which did not
  exist when this was written; writing them at 6a would have meant a `[live]`
  catalog entry for an absent suite, which is the same class of defect as the
  risk-4 line above. Genuinely 6b's, and deferring them is what the README
  licenses.

  **The `viz` section is written in 6b-i, not deferred again**, because the
  reason for deferring it was that the suite it names did not exist and 6b-i is
  the branch that creates it — the deferral would otherwise be false about its
  own tree. The seam split is precisely why this obligation moved earlier. It
  lands with `test_viz_scene.py` marked `[live]` and `test_viz_svg.py` /
  `test_cli_draw.py` marked `[planned]`; the "What we do not test" entry for the
  stylesheet and the rejected golden-file comparison stays with 6b-ii, which is
  the half that writes the stylesheet.

### Owed by 6a

**`project_structure.md`** -- two statements `94f94e2` made false.

1. The `bindings/core.cpp` section said the module "Currently binds the
   `Point2`/`Point3` value types". False as of `94f94e2`. Rewritten for the
   surface that exists, including the statement that `PslgBuilder` is
   deliberately not exposed and why.
2. The directory listing showed `src_python/tin_engine/` without `viz/`, which
   now exists. Added, with the two modules that are actually there and the rest
   marked as 6b's.

The Python API surface section also gains a line for `tin_engine.viz`; the
`draw` command is 6b's and is named there as pending rather than described.

### Owed by 6b

**`ROADMAP.md` — rewritten in 6b-ii. Done.** The defect was larger than a missing line.

The file's last commit is `f4efb6b`, *Use meshio for writing to file*, dated
**9 November 2018** (`git log -1 --format='%ad %s' f4efb6b`). Its three bullets
describe an OFF reader, a sun-ray file format and a shading script — the legacy
CGAL project, not the post-CGAL rebuild. It mentions no increment, no PSLG, no
CDT and no noder. `grep -rn ROADMAP . --exclude-dir=legacy --exclude-dir=.git`
returns **no reference from any other document** — run before this file existed;
run after, this file is the only hit, which is the point.

So the observation that increments 1 through 5a landed without touching it is
not negligence. It is evidence that `ROADMAP.md` is not the file that states
sequence — `docs/increments/` is, and `README.md` there says as much. The defect
is that an unreferenced, seven-year-stale roadmap sits at the top level where a
new reader will find it first and be misled by it.

Two honest repairs: delete it, or rewrite it as an index. **Ruled: rewrite as an
index, and done on this branch** -- one table, one row per increment, carrying
the number, one line, the status with its merge commit, and the path to the
record; because the one thing `docs/increments/` genuinely lacks is a single
screen showing what shipped and what is next. The rewrite carries, per
increment, only: the number, one line, the status, and the path to the record —
and **nothing that duplicates an increment file**, because duplicated detail is
exactly how the current file rotted. It gains this increment and `05b` as the
two open entries, with the ordering ruling above stated in one line so the
number-versus-ship-order gap is not a puzzle for the next reader.

**`project_structure.md`** -- the remaining entry, which is 6b's because it
describes 6b's files.

1. The directory listing's `viz/` entry gains `scene.py`, `style.py`, `svg.py`
   and `fixtures.py` as they land, and loses the "(6b)" marker. **Done for
   `scene.py` in 6b-i**, which also corrected that line's signature: shipped is
   `build_scene(pslg, mesh=None, ...)`, not `(mesh, pslg)`. **The remaining
   three, and the `__init__.py` line that still promised re-exports as future
   work, are done in 6b-ii.**
2. The Python API surface section's `tin_engine.viz` line gains the `draw`
   command once `cli.py` carries it. **Done in 6b-ii**, and it states the two
   things about the command that are not guessable: that the role mapping onto
   `ChainRole` is the composition root's, and that the fixture itself is the
   `PslgLike`.

**`testing.md`**

3. Add a `viz` section to the invariant catalog naming `test_viz_scene.py` as
   the invariant-critical suite. **Done in 6b-i**, with `test_viz_scene.py`
   `[live]` and the two ordinary suites `[planned]`; the split-by-owing-PR note
   at the top of this section says why it moved earlier.
4. State explicitly in "What we do not test" that the stylesheet is not tested
   and that a golden-file SVG comparison is rejected — otherwise that section
   leaves the gap open for someone to fill. **Done in 6b-ii**, which is the half
   that writes the stylesheet; asserting at 6a or 6b-i that an unwritten
   stylesheet is untested would have been the same defect in the other
   direction. The `viz` section also moves from `[partly live]` to `[live]`,
   with `test_viz_svg.py` and `test_cli_draw.py` marked `[live]` and what they
   pin stated in one sentence.
