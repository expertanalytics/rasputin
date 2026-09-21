# Increment 5c — wiring the noder through

**Status: design.** Written before `@tester` is briefed, per
`docs/increments/README.md`. No production code and no tests are part of this
document.

**Branch point:** master at `93fc772` — increment 5b merged as PR #76. Written
against the **merged tree and the shipped headers**, not against 5b's
description of them, which is why several of 5b's forward statements are
corrected below rather than carried.

## What this increment is for

The user's destination, in their words: *"the case where a road crosses a river,
or enters into a forest."*

5b delivers that at the C++ level and **nothing calls it.** `node<K>` is
instantiated by tests and by no production caller;
`terrain::cdt::triangulate` still takes a `Pslg`, so the only route from Python
into the triangulator is the route that refuses a crossing. Today:

```
$ rasputin draw not-noded --out /tmp/x.svg
CdtStatus.NotNoded, 0 triangles, the input PSLG drawn alone
```

5c is the increment that makes that command draw a mesh. Its product is not a
type and not an algorithm; it is a **path** — `build_pslg` → `node` →
`triangulate` → `build_scene` — and the thing that makes it worth a round is
that at the end of it a person can open an SVG and judge a noded crossing by
looking at it. `test_noding_node.cpp:171` already proves the crossing splits;
no human has seen it.

Keep that visible while reading the rest: every ruling below is in service of
one command producing one picture.

## Prior art in `legacy/`

**Nothing to carry across, and `@migration-expert` is not spawned** — but the
negative answer is not the useful part, and it was not reached by assuming it.
5b's section answered the question for the *noder*; 5c's subject is different
(a binding surface and a composition root), so the greps are different.

**1. The noding capability itself.** 5b's greps, re-run on this tree:

```sh
grep -rln "noding\|snap\|intersect\|broad phase\|split" legacy/
```

returns six files — `legacy/rasputin/triangulate_dem.h`, `legacy/bindings.cpp`,
`legacy/rasputin/reader.py`, `legacy/rasputin/geometry.py`,
`legacy/rasputin/gml_repository.py`, `legacy/tests/test_polygons.py` — and

```sh
grep -rni "snap" legacy/ | wc -l
grep -rni "noding\|broad.phase\|hot.pixel" legacy/ | wc -l
```

both print **0**. Unchanged from 5b: the legacy tree's entire noding capability
was one CGAL template argument, `Exact_predicates_tag` at
`legacy/rasputin/triangulate_dem.h:49`.

**2. The binding surface, which is 5c's largest row.**

```sh
grep -rln "pybind11\|PYBIND11_MODULE" legacy/
```

returns **two** files: `legacy/bindings.cpp` (406 lines) and
`legacy/rasputin/triangulate_dem.h`. What they show is a **counter-example, and
it is worth more than a blank answer.** `legacy/bindings.cpp:187-194`:

```cpp
m.def("make_mesh",
    [] (const R& raster_data, const P polygon, const std::string proj4_str) {
        return rasputin::mesh_from_raster(raster_data, polygon, proj4_str);
    }, py::return_value_policy::take_ownership)
```

A **CRS string crossing into C++**, which `CLAUDE.md` §2 now forbids by name;
and `:252-316` bind `CGAL::SimplePolygon`, `CGAL::Polygon` and
`CGAL::MultiPolygon` directly, so the library's vocabulary *was* the Python API.
There is no firewall to port. The rule 5c works under — the binding layer maps
core types to Python and nothing else crosses — exists because of this file.

**3. The composition root and its knobs**, which is where `--snap-spacing` has
to be judged.

```sh
grep -rln "argparse\|click\|typer\|__main__\|entry_points\|console_scripts" legacy/
```

returns `legacy/rasputin/application.py`, `legacy/rasputin/geo_tiff_reader.py`,
`legacy/rasputin/web_visualize.py`. The whole legacy CLI surface for meshing is
`legacy/rasputin/application.py:37-50`: `-x`, `-y`, `-polyfile`,
`-target-coordinate-system` (default `EPSG:32633`), `-ratio` (default `0.4`),
`-override`, `-land-type-partition`, `uid`, `-silent`.

**There is no tolerance, epsilon, spacing or precision option anywhere in it.**

```sh
grep -rni "tolerance\|epsilon\|\beps\b\|precision" legacy/
```

returns only `legacy/rasputin/tin_repository.py`'s XML `Precision="8"`
attributes, `legacy/rasputin/solar_position.h`'s astronomical obliquity, and one
line that is the real find: `legacy/rasputin/triangulate_dem.h:421`

```cpp
double eps = pow(pow(delta_x, 2) + pow(delta_y, 2), 0.5) * 1e-10;
```

a hard-coded robustness constant inside a header, invisible to every caller —
with its counterpart in the tests, `legacy/tests/test_mesh.py:129`'s
`polygon.buffer(1e-10)  # Account for fixed float precision`, a test working
around a constant it could not name.

**That is the prior art, and it is an argument rather than an absence.** The
legacy product's geometric robustness policy lived in two places a user could
not reach: a template argument and a magic epsilon. 5c's `--snap-spacing` is the
same policy, made nameable, defaulted at the composition root, and pointed at by
three of `NodeStatus`'s nine rows. Nothing is carried across; the shape of what
was there is the reason for the shape of what replaces it.

## The blueprint: data flow and the four boundaries

`@architect`'s obligation before `@developer` writes anything. The arrows are
data; the horizontal rules are the boundaries, and each one is a rule stated
elsewhere in the project rather than invented here.

```
  viz/fixtures.py            a Fixture: float64 (N,2), string roles, int masks
  (or, later, a reader)                 |
============================ cli.py is the ONLY module below this line that
                                        |  imports _core AND viz.style
  cli.py  ── ROLES ──────────► build_pslg(vertices, chains)  -> PslgBuildResult
                                        |                       (Pslg | diagnostics)
          ── --snap-spacing ──► node(pslg, spacing)          -> NodeOutcome
                                        |                       (NodedPslg | status)
          ── --delaunay ──────► triangulate(noded, delaunay) -> CdtOutcome
                                        |                       (IndexedMesh2 | status)
          ── Attempt(source, closed_roles, mesh, ok, ...) ───►
============================ bindings/core.cpp is the ONLY pybind11 TU
  viz/scene.py   build_scene(PslgLike, MeshLike) -> Scene   [imports no _core]
  viz/svg.py     render_svg(Scene, SvgStyle)     -> str     [holds the stylesheet]
```

Four boundaries, none of which 5c may bend:

1. **The I/O boundary** (`CLAUDE.md` §2). No path, no file, no codec and **no
   CRS** crosses into C++. `node()` takes a `Pslg` and two numbers. The legacy
   counter-example is `make_mesh(raster, polygon, proj4_str)`, quoted above.
2. **The pybind11 firewall.** `include/terrain/**` gains no pybind11 include and
   no pybind11 type. Every line 5c adds to the binding layer is in
   `bindings/core.cpp`.
3. **`viz/` imports no `_core` and holds no vocabulary.** `viz/protocols.py`
   types `ChainLike.role` as `object` and `properties` as `int` for exactly this
   reason. 5c feeds `build_scene` a `NodedPslg` — a **core type crossing into
   `viz/` as a structural `PslgLike`**, which is legal precisely because `viz/`
   never names it.
4. **`cli.py` is the sole composition root.** Every policy choice 5c makes —
   which spacing, which kernel, which graph the picture is drawn from, what a
   failure looks like — is made here or in `bindings/core.cpp`, and nowhere
   deeper.

**Testability at each seam, which is what licenses the shape.** `node()` is
callable from pytest with a `Pslg` and a float and no CLI; `build_scene` is
callable with a hand-built `PslgLike` and no extension; `_triangulated`'s
product is a value, not a side effect, so the CLI's failure presentation is
testable without writing an SVG. Nothing in 5c requires the meshing kernel to
be instantiated in order to test the layer above it.

## The C++ signature change, and the file 5b's list does not contain

5b priced this as "~6 lines" across four files: `cdt/triangulate.hpp`,
`cdt/detria_backend.hpp`, `src/cdt/detria_backend.cpp`, plus the include swaps
(`05b-noder-driver.md:1631-1640`). **There is a fifth**, and it is the one with
a design decision in it:

```sh
grep -rn "const Pslg&" include/terrain/cdt/
```

returns `triangulate.hpp:52`, `triangulate.hpp:59` — and
`constrained_edges.hpp:48`, `ConstraintEdgeSet::ConstraintEdgeSet(const Pslg&)`,
which `src/cdt/detria_backend.cpp:128` constructs from its argument. Retyping
the backend without retyping this does not compile.

**Ruling: `ConstraintEdgeSet`'s constructor becomes a template over a small
concept; `triangulate` does not.** The two halves are deliberate and they point
in opposite directions.

```cpp
// include/terrain/cdt/constrained_edges.hpp
template <typename G>
concept ChainedGraph = requires(const G& g, std::size_t c) {
    { g.chains().size() } -> std::convertible_to<std::size_t>;
    { g.edge_count(c) } -> std::convertible_to<std::size_t>;
    { g.indices_of(c) } -> std::convertible_to<std::span<const std::uint32_t>>;
};
```

* **Why the helper is structural.** `ConstraintEdgeSet` reads exactly three
  members and `Pslg` and `NodedPslg` both have all three. Retyping it to
  `const NodedPslg&` instead would force
  `tests/cpp/unit/test_cdt_constrained_edges.cpp` — 249 lines, twenty-odd
  `Pslg` fixtures, several of them deliberately degenerate — to route every
  fixture through `node<K>()`, which is the architecture's own prohibition:
  *a component that cannot be unit-tested without instantiating the meshing
  kernel is rejected*. Worse, half those fixtures **cannot survive noding** (a
  ring whose vertices repeat, a chain sharing an edge with another), so the
  suite would not be adapted, it would be deleted. A concept keeps the mask
  computation testable on hand-built input forever.
* **Why the entry point is not.** `terrain::cdt::triangulate` takes
  `const NodedPslg&` and nothing else. The whole architectural product of 5b+5c
  is that **un-noded input becomes unrepresentable at the entry point rather
  than diagnosed inside it** (`core/noded_pslg.hpp:36-40`). Templating
  `triangulate` on `ChainedGraph` would restore `Pslg` as a legal argument and
  undo the increment while leaving every test green.
* **Why the concept lives in `constrained_edges.hpp` and not in a new header.**
  It has exactly one consumer. A `terrain/core/chained_graph.hpp` invites a
  second, and the second will be `triangulate`.

`ChainedGraph` is also the first *written* statement of the structural equality
that `NodedPslg` was given on purpose — `core/noded_pslg.hpp:41-45` gives two
reasons for that shape (the CDT wrapper, and `viz/protocols.py`'s `PslgLike`)
and both were prose. This is the compiler checking it.

**`src/cdt/detria_backend.cpp` needs no body change.** Verified line by line: it
touches `pslg.vertices()` (`:91`, `:154`), `pslg.chains()` (`:93`, `:95`),
`pslg.indices_of(c)` (`:94`), `pslg.edge_count(c)` (`:110`) and
`ConstraintEdgeSet{pslg}` (`:128`). `NodedPslg` has every one with the same
signature and the same meaning. The change is the parameter type, the include,
and nothing else — 5b's "the wrapper body is untouched" is correct.

## The binding surface

`06-cdt-viewer.md`'s rule stands: the binding surface is exhaustive and an
addition needs a reason in an increment file. This is that file for the noder's
half. Everything below goes in `bindings/core.cpp` and nowhere else.

### What crosses, and what does not

| Crosses | As | Why |
|---|---|---|
| `NodeStatus` | `py::enum_`, 9 rows | the precedent `ChainRole`, `PslgError` and `CdtStatus` set |
| `describe(NodeStatus)` | an **overload** of the existing `describe` | below |
| `NodedPslg` | opaque class, `PslgLike`-shaped, plus 3 accessors | `viz/` joins on it |
| `NodeOutcome` | `pslg` / `status` / `message` / `ok()` | mirrors `PslgBuildResult` and `CdtOutcome` |
| `node(pslg, spacing, max_rounds=4)` | free function, GIL released | the driver |

| Does **not** cross | Why |
|---|---|
| `SnapGrid` | a bound class would carry `snap`, `world`, `cell_min`, `cell_max` — a grid vocabulary in Python next to no consumer. `NodedPslg.grid_spacing -> float` is the one value a caller needs, and it is the value it was handed. |
| `GridPoint` | integer lattice coordinates are an implementation detail of the noder; Python sees `vertices`, which are already `world(g)`. |
| `NodeOptions` | two scalars. A bound options struct is a second spelling of two keyword arguments. `CdtOptions` is not bound either. |
| `NodedPslgBuilder` | it is the verifier, and Python constructing a candidate is Python constructing an unverified `NodedPslg` by another name. The private-constructor proof (`core/noded_pslg.hpp:145-148`) must not have a Python back door. |
| a composite `mesh(vertices, chains, spacing)` | composition is `cli.py`'s, per boundary 4. A convenience entry point here would be a second composition root, and the first thing it would swallow is which failure the caller got. |

### `describe` becomes an overload set, and that is a ruling

Two enumerations now have prose. The alternatives are one overloaded `describe`
or two names (`describe` / `describe_node_status`).

**Overload.** A caller holding `outcome.status` does not know or care which
enumeration it came from at the point of printing — `cli.py` prints a status
band the same way for all three refusal layers, and two names force it to branch
on type to choose a function, which is precisely the branch the band does not
otherwise need.

The safety question is real and it is answered by pybind11's two-pass
resolution: the first pass runs every overload with conversions **disabled**, and
`py::enum_` registers a distinct Python type per enumeration, so
`describe(CdtStatus.Ok)` and `describe(NodeStatus.Ok)` resolve exactly even
though both are integer `0` underneath. `@tester` pins that with a call of each
and a `describe("Ok")` that must still be a `TypeError` — and the existing
docstring's claim, "Raises TypeError for anything that is not a CdtStatus", must
be widened when the second overload lands or it becomes false on the day it is
written.

### `NodedPslg`'s binding, and the duplication it would otherwise create

`py::class_<Pslg>` (`bindings/core.cpp:394-436`) binds `vertices`, `chains`,
`chain_indices` and `indices_of`, the last with an `IndexError` guard replacing
a debug assert. `NodedPslg` needs the same four with the same semantics, because
`viz/protocols.py`'s `PslgLike` is exactly those four and the renderer must be
able to consume either.

**Ruling: one function template binds the shared four for both types.**

```cpp
template <typename T>
[[nodiscard]] py::class_<T> bind_pslg_like(py::module_& m, const char* name, const char* doc);
```

`Pslg`'s existing block is replaced by a call to it; `NodedPslg`'s extra
accessors are chained onto the returned `py::class_`. The reason is not the
~25 lines saved. It is that **`PslgLike` is a structural contract with three
implementations** — `Pslg`, `NodedPslg` and `viz.fixtures.Fixture` — and a
member added to one hand-written block and not the other is a drift that
`tests/python/test_viz_protocols.py` catches only for whichever type that suite
happens to name. One template makes the two bindings the same object.

The three accessors `NodedPslg` adds:

* `grid_spacing -> float`, from `grid().spacing()`. Not a `SnapGrid`, per the
  table above.
* `edge_properties -> NDArray[np.uint32]`, a read-only `(E,)` zero-copy view,
  built the way `point_view` is built: a `static_assert` that
  `sizeof(EdgeProperties) == sizeof(std::uint32_t)` and
  `std::is_standard_layout_v<EdgeProperties>`, then a reinterpret. **A view of
  bare masks, not of a bound `EdgeProperties`** — increment 7 ruled that
  `EdgeProperties` is deliberately not bound because Python already has an
  integer with `|`, `&` and `bit_count()`, and `Chain.properties`
  (`bindings/core.cpp:349-350`) already crosses as `bits()`. Consistency with that
  is the whole reason, and the `static_assert` is what stops the reinterpret
  being a silent lie if the type ever grows a member.
* `node_of_input_vertex -> NDArray[np.uint32]`, a read-only `(N,)` view. This is
  guarantee 9's replacement and the **only** thing that maps an input vertex id
  onto a node id. Without it a caller who built the `Pslg` cannot say where
  vertex 4 went, and every such question becomes a coordinate search.

`edge_base` is **not** bound. It is an accessor whose value a Python caller can
compute as a prefix sum, and the C++ header's own reason for existing
(`core/noded_pslg.hpp:133-137`: "no caller holding (c, k) recomputes a prefix
sum") applies to C++ callers walking spans, not to a Python caller who has the
dense array and `chains`. If a Python consumer turns up that indexes by `(c, k)`
it is one line to add; adding it now is a surface with no consumer, which is
what the "what does not cross" table exists to refuse. **Recorded as a decision
so it is not read as an oversight.**

### `node()`: the GIL, the kernel, and the two arguments

```python
def node(pslg: Pslg, spacing: float, max_rounds: int = 4) -> NodeOutcome
```

* **GIL released**, exactly as `triangulate` is (`bindings/core.cpp:535-543`),
  and 5b's header states the property that licenses it:
  `node<K>` "IS A PURE FUNCTION of (pslg, options). No statics, no caches, no
  global state" and two calls from two threads produce bit-identical output
  (`include/terrain/noding/node.hpp:5-11`). The header says that reasoning was
  done "so that it is safe to release the GIL around at 5c"; this is 5c taking
  it. Nothing Python-owned is touched between release and reacquire: the input
  `Pslg` is a C++ object the Python wrapper owns, and the outcome is converted
  after the lock returns.
* **`terrain::pred::DefaultKernel`**, the same kernel `build_pslg` picks
  (`bindings/core.cpp:517`). The kernel is a compile-time parameter and the
  binding layer is where this project chooses it; no kernel parameter crosses,
  by `06-cdt-viewer.md`'s ruling.
* **`spacing` has no default here**, mirroring `NodeOptions::spacing`, which has
  none "EVER" and says why (`include/terrain/noding/node.hpp:83-89`). The
  binding is not the composition root. `cli.py` is, and it defaults.
* **`max_rounds` defaults to 4**, mirroring the C++ default. That is not a
  policy choice being made twice: it is the binding reproducing the C++
  signature, which is what every other `py::arg(...) = ...` in this file does.
* **Failure is data, not an exception**, like `build_pslg` and unlike a
  marshalling error. A refused noding is a fact about the terrain; only a
  mis-shaped argument raises.

### `_core.pyi`

Stubs for all of the above, or `mypy --strict` fails at `cli.py`'s first call.
Two things in it are not mechanical:

* `describe` becomes **two `@overload` stubs plus an implementation stub**, the
  pattern `cross` already uses in this file (`_core.pyi:62-67`) and the reason
  the file is hand-written rather than generated (`_core.pyi:3-5`).
* `NodeOutcome.pslg` is `NodedPslg | None`, matching `PslgBuildResult.pslg`'s
  stub — the `Optional` is the type system carrying "engaged iff status == Ok",
  and it is what makes `cli.py`'s `if outcome.pslg is None` a narrowing rather
  than a defensive check mypy cannot see through.

## `cli.py`: the composition root makes four choices

### 1. `--snap-spacing`, and its default

`node.hpp:83-89` is explicit that C++ has no default spacing and that "a
default in `cli.py` at 5c is the composition root making a policy choice and is
a different thing". This is that choice, and it must be a named module constant
with its reasoning attached, not a literal inside a `typer.Option`:

```python
#: Metres. The gallery's coordinates are UTM-shaped (viz.fixtures.ORIGIN), so
#: this is a millimetre on the ground.
DEFAULT_SNAP_SPACING = 1e-3
```

Why a millimetre, stated as a range rather than as a magic number:

* **Upper bound — ring collapse and false merging.** Coarser spacing merges
  features that are genuinely distinct. `test_noding_node.cpp:378` is the
  demonstration: at a decimetre grid, a road four millimetres from a river
  becomes *the same edge* with both properties. That is correct behaviour and
  the wrong default for a gallery whose job is to show what the input said.
* **Lower bound — the grid's range.** `kMaxGridIndex` is `2^51`
  (`core/snap_grid.hpp:63`), so at spacing `s` the representable coordinate is
  `s * 2.25e15`. At `1e-3` that is `2.25e12` metres, four orders of magnitude
  above any projected coordinate; `CoordinateOutOfRange` is unreachable from
  terrain data at this default. It is *not* unreachable at `1e-12`, which is why
  the option exists rather than the constant being hidden.
* **It must be measured, not asserted.** `@tester` runs the whole gallery at
  the default and requires `NodeStatus.Ok` for every fixture the validator
  accepts. A default that refuses a shipped fixture is a broken default, and
  that is a test rather than a paragraph.

**A non-finite or non-positive `--snap-spacing` is a `typer.BadParameter`, not a
status.** The engine also refuses it — `InvalidSnapSpacing`, step 0 of the
driver — and that is not duplication, because the two answer different
questions to different callers. `typer` rejects a **usage error** and exits 2
without drawing anything, the same line `_destination` already draws between a
bad `--out` and a refused fixture. The C++ check is a **library precondition**
with callers that are not this CLI, and it stays exercised by the Python binding
suite calling `node()` directly with `0.0` and `float("nan")`. Neither can be
deleted in favour of the other: delete the CLI guard and `rasputin draw
--snap-spacing -1` silently produces a picture of a failure that is not about
the terrain; delete the C++ guard and the library has an unchecked precondition.

### 2. `--max-rounds` is **not** exposed, and `describe` loses a lever

The CLI exposes the spacing and not the cap, deliberately.

The cap is a bound on a loop, not a parameter of the answer — at any value the
outcome is "the same mesh or `NotConverged`, never a different mesh"
(`node.hpp:91-95`). So raising it can only ever convert a refusal into a
success, and `node.hpp:34-46` records the one case where it provably cannot: a
constraint edge through a lattice corner with nodes on both flanking cells has
**period two** under the closed-cell predicate, so `Ok` is unreachable *at any
cap*. Handing a user a knob whose visible effect on the one failure they are
most likely to meet is "the same refusal, slower" is worse than not having it.

**And the header's sentence goes with it. `describe(NodeStatus::NotConverged)`
today names two levers** — "use a finer spacing **or raise the cap**"
(`noded_pslg_builder.hpp:106-108`). **Ruling: it names one.** `@developer` makes
this edit in the green commit:

```cpp
            return "the constraint set did not settle within the round cap; "
                   "use a finer spacing";
```

**This is not a compromise struck between the CLI's knobs and the header's
prose, and reading it that way would get the next one wrong.** The sentence was
wrong on its own terms, for every audience, before the CLI existed. `Ok` is
unreachable at any cap on the corner graze — that is a property of the
predicate, not of who is calling — so a caller who *can* turn the cap is told to
turn it into the same refusal. The design already knew: `05b-noder-driver.md:1478`
reads "raise the cap — **except on risk 16's corner graze, where neither
helps**". The header states the unqualified form. The design was right and the
header drifted from it; 5c is where a human first reads the header's words, so
5c is where the drift surfaced. The CLI did not create the problem, it exposed
it.

The shortened text is still true for a library caller, because it states the
lever that holds for all of them and no more. What is dropped is advice that was
never universally true, not advice specific to the CLI.

**One implementation constraint, and it is a citation rather than a style
preference: the arm stays at exactly two source lines.**
`05b-noder-driver.md:1560` cites `noded_pslg_builder.hpp:217` as the site of a
recorded mutation, and that line number is downstream of this string. Collapsing
the two lines into one moves it and falsifies a mutation record this branch is
not re-running. The replacement above is two lines for that reason.

`@tester` and `@developer` should read the whole of this as: the presentation
text for `NotConverged` (below) **must not promise the cap**, because the CLI
does not turn it and, for the deferred case, turning it would not help anyone
who did.

### 3. What a user sees when the noder refuses — including the corner graze

5c is the first increment where a human meets `NodeStatus`. The presentation is
the one the CLI already has, extended by one layer rather than reinvented:

* **the picture**: the input drawn alone, in role colours, with findings
  suppressed — `SceneKind.FAILED`, exactly what `not-noded` renders today
  (`viz/scene.py:235-249`). A drawn failure is still a drawn picture and the
  exit code stays 0; a caller who wants the verdict reads the band.
* **the band**: `status` is the enum's `.name`, `message` is
  `describe(status)` followed by the engine's own `outcome.message` — the same
  two-part join `_triangulated` already builds for a `CdtStatus`
  (`cli.py:132`).

For the corner graze specifically, the band reads:

```
NotConverged   the constraint set did not settle within the round cap; use a
               finer spacing   the constraint set did not settle in 4 rounds
               at spacing 0.001: <the last round's refusal>
```

That is the shortened arm ruled above, printed. The band names one lever
because the header now names one.

**The default does not dodge 5d, and pretending otherwise would be the wrong
ruling.** `05d-corner-graze.md:30-45` shows the orbit arising from *ordinary GIS
input* — a 45° road and two rivers, every coordinate a round decimetre, at
`SnapGrid{0.1}`. The trigger is input already rounded to the grid's own
resolution, and **every** spacing has such inputs; choosing `1e-3` moves which
data set is unlucky, not whether one exists. So the presentation above is not a
stopgap for an exotic case, it is the answer a real user will read, and that is
why it is designed rather than defaulted to a traceback.

**The predicate does not change** — `segment_meets_cell`'s closed cells are 5a's
shipped behaviour and `docs/increments/05d-corner-graze.md` owns the fix. What
5c decides is only what the user is told, and the answer is: *the noder's own
words, presented as a diagnosis about the data, with the spacing as the lever.*
It is not presented as a bug, because for the ordinary snap-induced case it is
not one; and it is not presented as fatal, because a finer spacing genuinely
cures that case. The user who has hit the period-two orbit will refine, get the
same answer, and have found 5d — which is the correct outcome of an increment
that deliberately does not fix it.

### 4. The scene's source — and 5b's "one-line change" is wrong

`05b-noder-driver.md:1656-1661` says feeding `build_scene` the noded graph "is a
one-line change in `cli.py`". **It is a two-field change, and doing the
documented half alone turns the whole gallery red.**

`cli.py:212` today passes the **fixture** as the `PslgLike`, together with
`CLOSED_ROLES = ("outer", "hole")` — *strings*, because a fixture's roles are
strings. `viz/scene.py:167` closes a ring with

```python
if len(walk) > 2 and any(chain.role == closed for closed in closed_roles):
```

A `NodedPslg`'s chains carry `ChainRole.Outer`, an enum, which equals no string.
Pass the noded graph while leaving `CLOSED_ROLES` alone and **no ring closes**:
every ring's closing edge is a masked mesh edge with no chain behind it, one
`MASKED_EDGE_WITHOUT_CHAIN` finding per ring, drawn in the alarm colour over a
correct mesh. That is the same catastrophe 5b predicted, reached by a mechanism
5b did not name — its analysis was about the *index spaces* diverging, and this
one is about the *role vocabulary* diverging.

**Ruling: the source and its closed-role vocabulary are one value, and the
composition root may not hold them apart.**

```python
@dataclass(frozen=True, slots=True)
class Attempt:
    source: PslgLike          # the most-processed graph that exists
    closed_roles: tuple[object, ...]   # the vocabulary `source`'s roles speak
    mesh: IndexedMesh2 | None
    ok: bool
    status: str
    message: str
```

`_triangulated` returns one of these instead of its current four-tuple, and
constructs it at each of the three exit points with both fields set together. A
source without its roles is then unrepresentable, which is the only form of this
fix that a later edit cannot undo by accident.

The rule for `source` is one sentence with no branch on the CDT's status:
**the most-processed graph that exists.** No `Pslg` — the fixture, with string
roles. No `NodedPslg` — still the fixture (a `Pslg`'s indices *are* the
fixture's, so nothing is gained by switching, and the fixture is the only source
available on the validator's own refusal path). A `NodedPslg` — the
`NodedPslg`, with `(ChainRole.Outer, ChainRole.Hole)`, whether or not the
triangulation then succeeded. The last clause is deliberate: when `node()`
succeeds and `triangulate` fails, the noded graph is the *truer* picture, since
it shows the splits the CDT was actually handed.

Passing `ChainRole` values into `build_scene` does not breach boundary 3.
`closed_roles` is typed `Sequence[object]` and `_chain_edges` uses only `==` and
`__hash__` — `viz/protocols.py:33-38` says so in as many words, and it is the
composition root supplying the enum that keeps `viz/` from naming it.

## The gallery: `not-noded` keeps its name and gains two colours

Three small changes outside `cli.py`, and they are what make the increment's
product visible rather than merely true.

**The fixture keeps its name.** After 5c, `rasputin draw not-noded` draws a
noded mesh, and the name looks wrong. It is not: the name describes **the
input**, which is still not noded, and `05b-noder-driver.md:180-184` is right
that the fixture's value is being *a regression fixture with two states* — the
picture changing between two commits with the same name is the increment's
product. Renaming it destroys that.

**Its module docstring becomes wrong and is fixed in this PR**, per
`docs/increments/README.md`'s rule that a documentation defect found during an
increment is fixed in that increment's PR or not recorded.
`viz/fixtures.py:21-26` says "**Three** of the eight are deliberate failures...
``not-noded`` and ``hole-in-hole`` are backend refusals". After 5c there are
**two**, and `not-noded` is the showcase. The inline comment at `:208-209` —
"the noder (5b) is what fixes it" — becomes a statement about a fix that has
landed.

**Its two breaklines get property bits.** Today both carry mask `0`
(`viz/fixtures.py:207,210`), while the C++ suite's fixture of the *same
geometry* gives them road and river (`test_noding_node.cpp:162-165`). Two
changed literals — `1` for the river bit and `2` for the road bit under
`DEFAULT_VOCABULARY`'s numbering — and the picture becomes the user's sentence:
a road crossing a river, in two colours, meeting at a constructed node.

`fixtures.py` writes the masks as literals with a comment, because `viz/` may
not import a vocabulary (`fixtures.py:193-195` already does exactly this for
`BREAKLINE`). Note that `test_noding_node.cpp:76-77` numbers them the other way
round — `bit(0)` is its road. That is not a defect to reconcile: **the C++ holds
no vocabulary at all**, by increment 7's ruling, and its test constants are
local names. It is written down here only so that nobody "fixes" the two into
agreement and thereby creates the second vocabulary the project has twice
refused.

**One stylesheet rule and one precedence entry.** `viz/svg.py:60` has
`line.river` and nothing for `road`; `cli.py`'s `_PRECEDENCE` is `("river",)`.
Without both, a road edge draws identically to the row above it and the legend
gains a row a reader cannot tell apart — `cli.py:73-79` records that this was
*measured* by drawing the gallery. So: one CSS line in `svg.py`, and
`_PRECEDENCE = ("river", "road")`, water before infrastructure, keeping
`cli.py:64-67`'s stated ordering rule. Two lines, and they are the difference
between a picture that shows the case and a picture that shows a grey crossing.

**What is *not* added: a merge fixture.** The road-along-a-river case — one
edge, property mask 3 — needs the two chains inside half a cell of each other,
which at a millimetre default means four-micron separation, or a per-fixture
spacing. `Fixture` has no spacing field and adding one makes the gallery hold a
policy the composition root is supposed to own. Deferred, named, and reachable
today from Python by calling `node()` directly with a coarse spacing — which is
what the binding suite does.

## Which statuses are diagnoses and which are self-checks

5c is where both enumerations reach a human, so the classification has to be
written down once. It is **not** re-derived in Python.

**Diagnoses — the caller acts on these, and the lever is named.**

| Status | What the caller does |
|---|---|
| `NodeStatus::InvalidSnapSpacing` | fix the argument (the CLI refuses it first, as a usage error) |
| `NodeStatus::CoordinateOutOfRange` | coarsen the spacing, or re-project |
| `NodeStatus::RingCollapsed` | use a finer spacing |
| `NodeStatus::NonSimpleRing` | repair the polygon upstream — the input was already wrong |
| `NodeStatus::NotConverged` | use a finer spacing; if it persists, this is `05d-corner-graze.md` |
| `CdtStatus::InvalidTopology` | fix the ring nesting |
| `CdtStatus::BackendFailure` | unclassified; read the message |

**Self-checks — "this is our bug", and no input is known to reach them.**

| Status | Why it is a self-check |
|---|---|
| `NodeStatus::MalformedOutput` | the driver establishes 11, 12, 13 and 16 by construction |
| `NodeStatus::RingDegenerateAfterSnap` | a winding flip requires a graze, so `NonSimpleRing` fires first; established by a bounded exhaustive search over 4-gons and 5-gons (`node.hpp:159-167`) |
| `NodeStatus::NotRun` | a default-constructed outcome; it never crosses the binding |
| `CdtStatus::NotNoded` | guarantee 14 makes it unreachable from a `NodedPslg` — **as of this increment** |
| `CdtStatus::DegenerateGeometry` | guarantees 12 and 13, plus `world` injectivity |
| `CdtStatus::MalformedInput` | the `Pslg` validator already excludes every condition behind it |
| `CdtStatus::NotRun` | as above |

**`RingDegenerateAfterSnap` stays a self-check.** 5b demoted it from diagnosis to
self-check and 5c does not reopen that: the check remains, the status remains
reachable in the enumeration, and no fixture is owed for it — "deleting a cheap
check whose job is to be unreachable is how a later change makes it reachable in
silence" (`node.hpp:165-167`).

**Ruling: 5c adds no classification code.** `describe()` already carries the
distinction *as the text the band prints* — "this is our bug" for
`MalformedOutput`, "self-check, no input is known to reach it" for
`RingDegenerateAfterSnap`, "repair the polygon upstream" for `NonSimpleRing` —
so a user who hits one reads the right thing with zero lines added, and a
Python-side table would be a second authority for a C++ fact, drifting the first
time a row is reworded.

**Where the predicate goes when something needs to branch**, recorded so the
next increment does not invent a place: an `is_self_check(NodeStatus)` beside
`describe` in `noding/noded_pslg_builder.hpp`, bound like `describe`. The caller
that will want it is an API worker deciding whether to retry or to page someone
— not a CLI printing a band. It is **not** added in 5c, because a predicate with
no caller is a surface with no consumer.

## What becomes unreachable in `CdtStatus`

Unchanged from `05-noder.md` and `05b-noder-driver.md:1598-1616`, and this is
the increment where it becomes true. `04-cdt.md` gains **one paragraph appended
beneath its mapping table**, leaving the table intact as the record of what was
true at increment 4:

- `PointOnConstrainedEdge` and `ConstrainedEdgeIntersection` — guarantee 14
  makes both unreachable, so `CdtStatus::NotNoded` joins `MalformedInput` as a
  self-check.
- `DuplicatePointsFound` — guarantee 12 plus `world` injectivity.
  `PolylineDuplicateConsecutivePoints` — guarantee 13. `DegenerateGeometry`
  loses both its reachable rows; `AllPointsAreCollinear` was already
  unreachable.
- **Reachable from a `NodedPslg`: `Ok`, `InvalidTopology`, `BackendFailure`.**

`describe(CdtStatus::NotNoded)`'s text (`cdt/result.hpp:45-48`) keeps naming the
noder as the fix, because at that point it is naming the thing that should have
been called. No enumerator is deleted.

**And the remaining statuses stay testable without the noder** — which is the
half of this that was not obvious and had to be checked. `NodedPslgBuilder` is
a **public** entry point (5b made it public precisely so a hand-built candidate
could be refused), and it "is not the Pslg validator": it does not ask for an
outer chain, does not check ring winding and does not compute nesting
(`noded_pslg_builder.hpp:45-49`). So a suite can hand-build a `NodedPslg` with a
hole outside its outline, get an `Ok` from the builder, and hand it to the
backend to get `InvalidTopology` — no `node<K>()` call, no snapping, no
kernel-wide instantiation. Testability of the CDT's failure surface survives the
signature change; it survives *because* of a 5b decision made for another
reason, and that is worth recording rather than rediscovering.

What does **not** survive: a test that reaches `NotNoded` or
`DegenerateGeometry` through the backend. Those inputs cannot be expressed as a
`NodedPslg` at all — a crossing is refused by guarantee 14(a), and a collinear
3-gon by 14(b), since the middle node's cell meets the closing edge. Those tests
are **deleted, not adapted**, and the deletion is authorised here so that it is
not read at review time as a suite being weakened.

## Gate C, discharged

`05b-noder-driver.md:472-483` obliges 5c's designer to re-estimate
`bindings/core.cpp` and `_core.pyi` **against what `NodeStatus`, `NodedPslg` and
`NodeOutcome` actually turned out to be**, before `@tester` is briefed, and to
move `cli.py` and `--snap-spacing` out if the total passes **550**.

**Gate C's premise held, and that deserves saying before its defects.** The
surface 5c must bind is exactly what 5b assumed: `NodeStatus` shipped as the
nine-row enumeration with a total `describe`, `NodedPslg` as the ten-accessor
type, `NodeOutcome` as the three-field struct with `ok()`. Not one row moved.
What moved is the conversion factor, not the API — so the gate fires on the
arithmetic and not on a redesign.

### Defect 1: the name it reserves is taken

Gate C says "move to a 5d". `docs/increments/05d-corner-graze.md` now exists and
is a different increment — the closed-cell period-two orbit, whose fix changes
5a's shipped predicate. **The overflow increment is `5e`**, and it is named here
so that the reservation is unambiguous whether or not it is ever allocated.

### Defect 2: the C++ factor is three-for-four, and the re-derivation

The sentence declining to double the C++ rows (`05b-noder-driver.md:1734-1736`)
rests on "increments 3, 4 and 5a all came in at or under their C++ estimates".
With 5b in the sample:

| Increment | Estimated | Measured | Ratio |
|---|---|---|---|
| 3 (`pslg_builder.hpp`) | 255 | 249 | 0.98 |
| 4 | 342 | 267 | 0.78 |
| 5a | 270 | 179 | 0.66 |
| 5b | 445 | 675 | **1.52** |

* **pooled** (sum measured / sum estimated): 1370 / 1312 = **1.04**
* **arithmetic mean of ratios**: 0.98
* **geometric mean of ratios**: 0.94
* **range**: 0.66 to 1.52, a spread of 2.3×

5b flagged this itself, in its reconciliation section, and stopped one step
short: it recorded "three-for-four, not three-for-three", said the sentence "is
weaker than when it was written", and did not re-derive the factor. Re-deriving
it is Gate C's assignment and is what follows.

**What that changes is not the centre but the spread.** "C++ comes in at or
under" was a claim about the *centre*, and by the pooled and mean figures it is
still roughly right — the centre is 1.0, not 0.8. The defensible reading is that
**the C++ estimator is unbiased and imprecise**, so a point estimate is worth
quoting and a *margin* computed from it is not. Adopted for this document:
**central factor 1.0, adverse factor 1.52** — the observed maximum, which is 5b's
own — applied to the C++ rows. The corresponding Python/binding figures, from
increment 6a: `bindings/core.cpp` 283/150 = 1.89, `_core.pyi` 125/60 = 2.08,
`scene.py` 194/90 = 2.16, **pooled 602 / 300 = 2.01**.

One figure in 5b's argument for not doubling the C++ rows is independently
wrong, and it is worth separating from the factor question: "5c's C++ row is ~6
mechanical lines" (`:1833-1834`) omits `constrained_edges.hpp`. The row is **18**
across five files. That is still small enough for the conclusion to hold — 18 at
the adverse 1.52 is 27 — but it triples the row the conclusion rests on, and it
was found by grepping the tree rather than by re-reading the design.

### The re-estimate, by measured analogue rather than by sketch-and-double

The doubling factor was always a proxy for a known defect in the *instrument*: a
line-by-line sketch of a binding block omits structural elements, so it
systematically under-counts. **5c does not need the proxy**, because every block
it adds has a structural twin already shipped in the same file, in the same
house style, for a type of the same shape. Each row below is priced against a
named analogue, measured with this project's C++ instrument
(`grep -vcE '^\s*(//|$)'`) and Python instrument (`grep -vcE '^\s*(#|$)'`):

These five commands were run before the rows below were written, and their
outputs are the figures quoted. **The line ranges expire, and this increment is
what expires them** — 5c edits `bindings/core.cpp`, so after the green commit
they address different text. They are a design-time measurement, reproducible at
master `93fc772` and nowhere later; the rule they encode, "price a binding row
against the shipped block of the same shape", is what a later reader should
carry, not these numbers.

```sh
# reproduce any row at 93fc772: the analogue's line range, this file's instrument
sed -n '394,436p' bindings/core.cpp | grep -vcE '^\s*(//|$)'   # py::class_<Pslg>  -> 39
sed -n '483,497p' bindings/core.cpp | grep -vcE '^\s*(//|$)'   # py::class_<CdtOutcome> -> 13
sed -n '322,334p' bindings/core.cpp | grep -vcE '^\s*(//|$)'   # py::enum_<CdtStatus>   -> 12
sed -n '534,551p' bindings/core.cpp | grep -vcE '^\s*(//|$)'   # m.def("triangulate")   -> 13
sed -n '165,186p' src_python/tin_engine/_core.pyi | grep -vcE '^\s*(#|$)'  # class Pslg -> 17
```

**`bindings/core.cpp`:**

| Row | Analogue (measured) | Est. |
|---|---|---|
| `py::enum_<NodeStatus>`, 9 rows, no per-value prose | `CdtStatus`, 7 rows = 12 | 14 |
| `describe(NodeStatus)` overload | `describe(CdtStatus)` = 4 | 5 |
| `bind_pslg_like<T>`, **net** over the `Pslg` block it replaces | `py::class_<Pslg>` = 39 | +18 |
| `py::class_<NodedPslg>`: class, doc, 3 extra accessors | `py::class_<Chain>` (`:341-361`) = 19 | 42 |
| `properties_view` helper | `point_view` = 9 | 11 |
| `py::class_<NodeOutcome>` | `CdtOutcome` = 13 | 16 |
| `m.def("node", ...)` | `m.def("triangulate")` = 13 | 20 |
| retyped `triangulate`, includes, usings | — | 9 |
| **Total** | | **135** |

**`_core.pyi`:**

| Row | Analogue | Est. |
|---|---|---|
| `NodeStatus` enum, 9 rows | `CdtStatus` = 9 | 12 |
| `class NodedPslg` | `class Pslg` = 17 | 24 |
| `class NodeOutcome` | `CdtOutcome` = 9 | 11 |
| `describe`: two `@overload` + impl stub | `describe` = 2 | 7 |
| `node()` | `triangulate()` = 1 | 5 |
| `triangulate` retype | — | 1 |
| **Total** | | **60** |

**135 and 60, against 5b's ~215 and ~90.** The direction is down, and the reason
is the instrument rather than optimism: 5b priced a *sketch* and doubled it,
this prices *shipped blocks* and does not. Two rows 5b's sketch did not contain
appear elsewhere in the table below and more than absorb the difference.

### Verdict

**~256 against Gate C's 550. The gate does not fire. `cli.py` and
`--snap-spacing` stay in 5c, and `5e` is not allocated.**

The adverse case is stated too, because 5b's lesson was that a row named as at
risk must be *re-priced*, not merely watched: applying the full measured 2.01 to
every Python and binding row and 1.52 to the C++ row gives **~505**. Still under
550, and under `CLAUDE.md` §2's 700 with room. Both the point estimate and the
adverse case clear the gate, which is what makes the verdict robust rather than
a near-miss argued in one direction.

The at-risk row, named here and re-priced rather than watched: **`bind_pslg_like`
plus `py::class_<NodedPslg>`, 60 of the 135.** It is the only row that *modifies*
shipped code rather than adding beside it, and modification is where an estimate
built from an analogue is weakest — the analogue is the thing being changed. It
is priced 50 % above the block it replaces for that reason.

### Gate D — 5c's own, and it is a budget with an obeyable seam

The one gate in this project's record that worked fired **before its suite
existed** (`05b-noder-driver.md:430-432`); the one that failed did so because the
committed suite spanned the seam it would have had to cut. So Gate D is armed
with that property removed by construction:

> **`bindings/core.cpp`'s addition is budgeted at 135 non-comment lines.**
> `@developer` measures the binding block **after writing it and before starting
> `cli.py`**. If it exceeds **200**, `cli.py`, `--snap-spacing`, the `Attempt`
> record and the gallery changes move to **5e**, and 5c ships the type crossing
> with its Python binding suite alone.

That measurement happens after the red commit, which 6b-ii says is too late —
**for a suite split.** This seam does not split a suite: the CLI surface and the
binding surface are separate test files by Gate B's one-file-per-surface rule
(`tests/python/test_cli_draw.py` already exists beside `test_core_cdt.py`), so
the 5e cut removes whole files from the PR rather than dividing one. That is the
property that makes this gate cheap enough to be obeyed, and it is the only
reason it is armed at all.

## Files and LOC

Instruments as declared in `05b-noder-driver.md`: C++ `grep -vcE '^\s*(//|$)'`,
Python `grep -vcE '^\s*(#|$)'`. Docstrings and `R"doc(...)"` bodies count as
code under both, which is the definition `CLAUDE.md` §2's *non-comment lines*
approximates.

| File | Contents | Est. |
|---|---|---|
| `include/terrain/cdt/triangulate.hpp` | `CdtBackend` and `triangulate` take `const NodedPslg&` | 3 |
| `include/terrain/cdt/detria_backend.hpp` | the declaration | 2 |
| `include/terrain/cdt/constrained_edges.hpp` | `ChainedGraph` concept, templated constructor | 10 |
| `src/cdt/detria_backend.cpp` | signature and include; body untouched | 3 |
| `include/terrain/noding/noded_pslg_builder.hpp` | `describe(NotConverged)` drops the cap lever; two lines in, two lines out | 0 |
| `bindings/core.cpp` | `NodeStatus`, `describe` overload, `bind_pslg_like`, `NodedPslg`, `NodeOutcome`, `node()`, retyped `triangulate` | 135 |
| `src_python/tin_engine/_core.pyi` | stubs for all of it | 60 |
| `src_python/tin_engine/cli.py` | `node()` call, `--snap-spacing`, `Attempt`, the noder's failure presentation, `_PRECEDENCE` | 38 |
| `src_python/tin_engine/viz/fixtures.py` | two property literals, the module docstring's failure count | 4 |
| `src_python/tin_engine/viz/svg.py` | one `line.road` stylesheet rule | 1 |
| **Total** | | **~256** |

Against `CLAUDE.md` §2's 700, and against Gate C's 550: comfortable at the point
estimate and clear at the adverse ~505.

**Two rows 5b's list does not contain**, both found by reading the shipped tree
rather than the design: `constrained_edges.hpp` (the fifth file with
`const Pslg&` in it — without it the branch does not compile) and the
`fixtures.py`/`svg.py` pair (without which the increment's product is a grey
crossing nobody can read). Together they are 15 lines and they are the reason
5c's total did not fall as far as the binding re-estimate alone would suggest.

**Not production code, and the round is not free.** The test churn is the
largest uncosted part of 5c.

**What this list enumerates, so that its count is checkable rather than a
judgement: every test file that already exists and that 5c changes** — not the
subset someone judged expensive, and not the new suites, which are in "What is
worth testing" instead. `tests/python/test_core_noding.py` is therefore absent
by construction: it is new. The membership rule resolves to a command, which is
how it should be re-derived rather than trusted:

```sh
git diff --name-only <branch point>..HEAD -- tests/ | grep -v test_core_noding
```

**Seven files, four of them C++** — three suites and the `cdt_cases.hpp` support
header the C++ suites share. An eighth entry is listed below precisely because
it is *not* one of the seven.

The earlier count of five was not a scope line, it was a list with no membership
rule, which is how `test_cdt_backend_seam.cpp` stayed off it while being
unbuildable without a change. A rule is what makes the next omission visible.

* `tests/cpp/unit/test_cdt_detria_backend.cpp` (~500 lines) and
  `tests/cpp/property/prop_cdt_invariants.cpp` (~430) build `Pslg` fixtures and
  call `triangulate<DetriaBackend>`. Every call site must route through
  `node<DefaultKernel>` or through a hand-built `NodedPslgBuilder` candidate, and
  the cases pinning `NotNoded` and `DegenerateGeometry` are **deleted**, per the
  ruling above.
* `tests/cpp/support/cdt_cases.hpp` is where most of that churn actually lands,
  and it is a **support header rather than a suite** — counted here because the
  membership rule is "test file", and named because a reader looking for the two
  suites' diff will find half of it in a third file. It gains a
  `node_fixture<K>()` beside `build_fixture<K>()`, a fixture-level spacing
  constant, and a `node_of(noded, input)` helper for guarantee 9's array. That
  last one is the load-bearing addition: **node ids are not input indices** —
  the node set is sorted by `GridPoint` — so every assertion in the two suites
  written against a literal vertex index has to route through it.
* `tests/cpp/unit/test_cdt_constrained_edges.cpp` (~250) is the **eighth file,
  and it is not one of the seven**: it is **unchanged**. That is what the
  `ChainedGraph` concept buys, and it is the concrete return on the ten lines it
  costs. It is listed under a rule that excludes it because an unchanged file is
  only evidence if someone says in advance that it will be.
* `tests/cpp/unit/test_cdt_backend_seam.cpp` appeared in no list until this
  revision. It has to change, and not by choice: `CdtBackend` is spelled
  in terms of the entry point's parameter, so retyping `triangulate` retypes the
  concept, and every fake backend in the file declares that parameter. Leave
  them at `const Pslg&` and the fakes the suite asserts *are* backends stop
  satisfying `CdtBackend`, which `STATIC_REQUIRE` reports at compile time. That
  is the finding, not a compile accident — the suite exists to be the place
  where a change to the seam's signature announces itself.

  Two substantive changes `@tester` made here, both **endorsed**, recorded so
  `@reviewer` audits them against a written reason:

  - **`UnNodedCdtBackend` is new**: the pre-5c signature kept alive as a type
    that must now *fail* the concept. Without it, the retype is revertible with
    every other assertion in the file staying green — nothing else in the suite
    can tell a concept that rejects `Pslg` from one that accepts it. It is also
    the first place a compiler checks 5c's architectural product: **un-noded
    input unrepresentable at the entry point rather than diagnosed inside it.**
  - **`FailingCdtBackend` returns `BackendFailure` instead of
    `CdtStatus::NotNoded`**, required by 5c's own ruling above that `NotNoded`
    becomes a self-check no `NodedPslg` can reach. A fake handing one back would
    read as a claim that it still can.
* `tests/python/test_core_cdt.py`'s
  `test_reports_a_non_noded_input_as_a_failure_status` asserts that
  `crossing_pslg` fails. After 5c that input cannot reach `triangulate` at all,
  so the case splits in two: the noding suite asserts the crossing *succeeds*,
  and what stays in the CDT suite is that handing `triangulate` a `Pslg` is a
  `TypeError` rather than a status. No line number is cited, because this
  increment is what moves it — locate it by name, or by
  `grep -rn non_noded tests/python/`.
* `tests/python/test_cli_draw.py` pins the `not-noded` fixture's presentation,
  which is the picture this increment changes.
* `tests/python/test_viz_protocols.py` is the seventh, and the one the earlier
  count missed for a different reason than the seam suite did: its *new*
  assertion — `NodedPslg` satisfying `PslgLike`, the third implementation — is
  listed under "What is worth testing" as added coverage, and that made it look
  like a file with no churn. It has churn as well. Its shared `mesh` fixture
  builds through `triangulate(result.pslg)` today, so 5c routes it through
  `node()` and every test unpacking that fixture changes arity with it. Added
  coverage and forced churn in one file is exactly the case a list without a
  membership rule drops.

All of it is `@tester`'s, in its own commits with the reason in the message, per
`docs/increments/README.md`: no test change hides inside an implementation
commit.

## What is worth testing

**Gate B's mapping applies: one test file per production surface.** Four
surfaces, four files, and the C++ half of this increment adds no new suite —
the existing CDT suites are *amended* rather than replaced.

| Surface | Suite | Invariant-critical? |
|---|---|---|
| `ChainedGraph` + retyped entry point | `tests/cpp/unit/test_cdt_backend_seam.cpp` (amended) | no |
| the `node`/`NodedPslg`/`NodeOutcome` binding | `tests/python/test_core_noding.py` (new) | **yes** |
| `cli.py`'s pipeline and failure presentation | `tests/python/test_cli_draw.py` (amended) | no |
| the structural join | `tests/python/test_viz_protocols.py` (amended) | no |

**That row used to name `test_cdt_constrained_edges.cpp`, and it contradicted
"Files and LOC", which says that file is unchanged.** The contradiction is
resolved in favour of *unchanged*, because that is the load-bearing claim: the
stated return on the ten lines `ChainedGraph` costs is that not one of that
file's twenty-odd hand-built `Pslg` fixtures needs adapting, and the only way to
demonstrate it is to not touch the file. A case added there would spend the
claim to buy a case that fits elsewhere.

Elsewhere is the seam suite. `ChainedGraph` is a seam between two graph types
and one consumer, and the retyped entry point is the same suite's subject
already; `test_cdt_backend_seam.cpp` is also registered against the plain test
helper, which is what lets it assert on `triangulate` without linking a backend.
`@tester` placed the case there and its header comment at
`test_cdt_backend_seam.cpp:44-52` states this argument in the file itself.
Endorsed; the design is what was wrong.

**The invariant-critical suite is the binding one**, and the reason is where the
silent failures live. Everything in `bindings/core.cpp` that 5c adds is a
lifetime or a reinterpret:

* **`edge_properties`' reinterpret.** A `(E,)` uint32 view over a
  `std::vector<EdgeProperties>`. If the `static_assert` on size and layout is
  wrong, or the stride is, the array reads plausible garbage — masks that are
  *almost* right are exactly the failure a picture cannot show. Mutants owed:
  the stride off by a factor, the length taken from `chains().size()` instead of
  the array's own, the `static_assert` deleted.
* **Every view's base object.** `06-cdt-viewer.md` risk 4: a view without its
  owner is a use-after-free the moment the outcome is dropped. The suite must
  hold `noded.edge_properties` and `noded.node_of_input_vertex`, drop every other
  reference to the `NodeOutcome`, force a collection and read — the pattern the
  existing mesh-view tests use.
* **Writeability.** Both new arrays must have `WRITEABLE` cleared; a mutant that
  forgets it lets Python write through into a type whose entire contract is
  immutability.
* **The GIL release.** A `node()` call from two threads on the same `Pslg`, with
  bit-identical output required — the property `node.hpp:5-11` claims and the
  only thing that makes the release safe.
* **`describe`'s overload resolution.** `describe(CdtStatus.Ok)` and
  `describe(NodeStatus.Ok)` must return their own enumeration's sentence even
  though both are `0`, and `describe("Ok")` must still be a `TypeError`.
* **`node()`'s refusals as data.** `spacing=0.0` and `float("nan")` return
  `InvalidSnapSpacing` rather than raising; a path-like argument in `pslg`'s
  place is a `TypeError`, matching `build_pslg`'s boundary.

**The end-to-end acceptance criterion**, which is `05b-noder-driver.md:143-146`'s
criterion with the call it was waiting for:

```sh
.venv/bin/python -c "
import numpy as np
from tin_engine import _core as c
V=np.array([[0,0],[700,0],[700,700],[0,700],[100,100],[600,600],[100,600],[600,100]],
           dtype=float)+np.array([430000.,6900000.])
ch=[(list(range(4)),c.ChainRole.Outer,0),([4,5],c.ChainRole.Breakline,2),
    ([6,7],c.ChainRole.Breakline,1)]
r=c.build_pslg(V,ch)
n=c.node(r.pslg, 1e-3)
o=c.triangulate(n.pslg)
print(n.status, o.status, o.mesh.triangle_count)
"
```

must print `NodeStatus.Ok CdtStatus.Ok` and a triangle count above zero, and the
mesh's vertex array must contain the snapped image of `(350, 350) + ORIGIN`
**once**. Note the rebuild obligation before any `pytest` that measures C++:
`.claude/REQUIRED-READING.md` — `pytest` does not rebuild the extension, and a
green run from a stale `.so` is byte-identical to a real one.

**And the criterion a person checks**, which is the increment's actual product:

```sh
rasputin draw not-noded --out /tmp/noded.svg
```

draws a **mesh** with **zero findings**, a road and a river in two colours, and
`CdtStatus.Ok` in the band. Zero findings is the load-bearing half: a single
`MASKED_EDGE_WITHOUT_CHAIN` means the scene was built from the wrong graph or
with the wrong `closed_roles`, which is the defect this design spent a section
on and the one a passing status will not reveal.

## Risks

1. **The scene is fed a `NodedPslg` with string `closed_roles`.** Ruled out
   structurally by the `Attempt` record above rather than by care; the residual
   risk is `Attempt` being introduced and then bypassed. Mitigation:
   `_triangulated` returns `Attempt` and `draw` unpacks nothing else, so there
   is no second path to pass a source.
2. **The binding row overruns.** Gate D, with a seam that costs the deletion of
   whole test files rather than the division of one.
3. **A stale `.so` makes the round's central claim unfalsifiable.** This is the
   one hazard with a recorded incident behind it — a pre-fix `bindings/core.cpp`
   was restored on disk and `pytest` reported 94 passed. `@tester` and
   `@developer` rebuild and reinstall explicitly, per `.claude/REQUIRED-READING.md`.
4. **The default spacing refuses a gallery fixture.** Caught by the
   whole-gallery test rather than by argument; if it fires, the answer is a
   different default, not a per-fixture knob.
5. **The C++ suites' amendment quietly weakens coverage.** Two suites lose the
   cases that pinned `NotNoded` and `DegenerateGeometry`. Authorised above, with
   the mechanism that makes each unreachable stated, so `@reviewer` audits the
   deletion against a written reason rather than against a diff.

## What this PR corrects elsewhere

Per `docs/increments/README.md` — found during this increment, fixed in this
increment's PR, or not recorded:

* **`docs/increments/04-cdt.md`** gains the unreachability paragraph beneath its
  mapping table. The table stays.
* **`docs/increments/05b-noder-driver.md`** is corrected in three places: the
  C++ signature change is **five** files, not four (`constrained_edges.hpp` is
  missing from `:1631-1640`); feeding the scene a `NodedPslg` is **not** a
  one-line change (`:1660`) and the second field is `closed_roles`; and Gate
  C's overflow increment is **5e**, not 5d (`:479`), because `05d` is now the
  corner graze. A fourth: its Python-suites paragraph cited a line number in
  `tests/python/test_core_cdt.py` for a test this increment splits, so the
  number is dropped in favour of the test's name and a `grep`, per
  `.claude/REQUIRED-READING.md` on resolved values. The paragraph's claim is
  still true; only the citation expired.
* **`src_python/tin_engine/viz/fixtures.py:21-26`** — three deliberate failures
  become two, and `not-noded` changes role.

**Not corrected, because 5b corrected it itself:** the C++ estimate paragraph at
`:1734-1736` rests on a three-for-three record, and 5b's own reconciliation
section already says so (`:1829-1836`, "three-for-four, not three-for-three"). It
did not, however, *re-derive the factor* — it said the sentence "is weaker than
when it was written" and left it standing. That is the half Gate C hands to this
document, and the table above is it. Editing the shipped reconciliation section
would be editing a record of what that round concluded; the successor estimate is
where the corrected factor belongs.
