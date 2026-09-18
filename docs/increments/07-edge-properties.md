# Increment 7 — edge properties

Status: design settled. **This increment lands before 5b, by the user's
ruling**, and after `increment6b-ii-renderer-red` merges — see "Sequencing".

It exists because of one correction from the user, quoted verbatim because
every ruling below is downstream of it:

> we should classify edges and not nodes. An edge might have several
> properties. At a coarse resolution, the same line segment can be both a road
> and a river.

`docs/increments/05b-noder-driver.md` ruled on the *representation* while
designing the noder, and ruled that widening it was not 5b's to do. This
document owns those rulings now; they are restated below with their reasons
rather than referenced, because a ruling that lives only in the document of a
different increment is a ruling that gets re-derived. What 5b did **not**
settle — the blast radius and its sequencing, the shape of `Chain`, the
boundary model, the renderer — is settled here.

## Prior art in `legacy/`

**Nothing on edges, and the legacy tree's silence is the positive finding.**

```sh
grep -rlni "is_river\|river\|road\|railway\|coastline\|landcover\|land_cover" legacy/
```

returns **nine** files:

```
legacy/rasputin/application.py
legacy/rasputin/globcov_repository.py
legacy/rasputin/gml_repository.py
legacy/rasputin/land_cover_repository.py
legacy/rasputin/tin_repository.py
legacy/rasputin/web_visualize.py
legacy/rasputin/wfs_repository.py
legacy/tests/test_gml_repository.py
legacy/tests/test_land_cover_repository.py
```

and

```sh
grep -rlni "bitset\|bitmask\|flags\|attribute" legacy/
```

returns **four**: `legacy/rasputin/mesh.py`, `legacy/rasputin/reader.py`,
`legacy/rasputin/tin_repository.py`, `legacy/rasputin/wfs_repository.py`.

Every hit in both lists is one of two things, and neither is a per-edge
property:

- **Land cover, as an enumeration of area classes, attached per face.**
  `LandCoverBaseType` (`legacy/rasputin/land_cover_repository.py:8`) is an
  `Enum`; `legacy/rasputin/globcov_repository.py:15` declares **23** enumerators and
  `legacy/rasputin/gml_repository.py:13` (CORINE) declares **44**
  (`grep -c "^    [a-z_0-9]* = [0-9]" ...` prints `23` and `44`). The
  repository interface is
  `constraints(self, *, domain: GeoPolygon) -> List[GeoPolygon]`
  (`legacy/rasputin/land_cover_repository.py:42`, implemented at `legacy/rasputin/gml_repository.py:181` and
  `legacy/rasputin/globcov_repository.py:132`) — **every legacy constraint is a polygon**.
- **Colour per face, and geometry partitioned by it.**
  `Geometry.split_by_colors` (`legacy/rasputin/geometry.py:46`) partitions
  faces by a colour attribute. That is the face-based channel
  `parallel_refinement.md` describes, and it is sound for area features.

**There is no linear-feature channel anywhere in `legacy/`.** Not a per-edge
attribute, not a river polyline, not a road polyline — `grep -rni "river"
legacy/` hits only land-cover *descriptions* (`legacy/rasputin/globcov_repository.py:50-62`,
strings like "regularly flooded"), and the only road in the tree is
`legacy/rasputin/gml_repository.py:18`'s `road_and_rail = 122`, a CORINE **area** class.

**Nothing is carried across, and the reason is the useful one.** The legacy
product could classify a forest and could not classify a road, because its only
semantic channel was face-based and a road has no interior to be inside. That
is precisely the gap `parallel_refinement.md`'s corrected "Edge metadata"
section names, arrived at independently, and the legacy tree is the evidence
that it is a real historical loss rather than a hypothetical one.
`@migration-expert` is **not** spawned for this increment: there is no legacy
source to read intent from, because the intent was never expressed.

The one number carried across is a **bound**: 23 and 44 land-cover classes,
against under ten linear features. That asymmetry is the whole argument for
`kMaxProperties = 32` below, and it is measured here rather than assumed.

## What this increment ships, in one sentence

`Chain::is_river`, one `bool`, becomes `Chain::properties`, an opaque 32-bit
set that `terrain::` cannot name; the names live in one Pydantic model in
Python; and the renderer draws at most one property stroke, chosen by a
precedence the composition root supplies.

## The rulings inherited from 5b, with their reasons

These were argued in `05b-noder-driver.md`, "Edge properties". They are not
re-derived; they are restated because this is the increment that owns them.

### `EdgeProperties` is a value type over `std::uint32_t` with no named members

```cpp
// include/terrain/core/edge_properties.hpp, namespace terrain
class EdgeProperties {
public:
    static constexpr unsigned kMaxProperties = 32;

    constexpr EdgeProperties() = default;                       // the empty set
    [[nodiscard]] static constexpr EdgeProperties bit(unsigned i) noexcept;
    [[nodiscard]] constexpr bool empty() const noexcept;
    [[nodiscard]] constexpr bool contains(EdgeProperties) const noexcept;  // superset
    [[nodiscard]] constexpr std::uint32_t bits() const noexcept;

    friend constexpr EdgeProperties operator|(EdgeProperties, EdgeProperties) noexcept;
    friend constexpr bool operator==(const EdgeProperties&, const EdgeProperties&) = default;
};
```

**No enumeration, no `River`, no `Road`, no named member of any kind.**
`terrain::` never spells a feature name. The mapping from bit position to
feature name is a Python concern and lives in a Pydantic model at the boundary,
alongside CRS metadata and everything else `project_structure.md:185` keeps out
of C++. This is **stricter than what ships today**, where `Chain::is_river`
spells "river" inside the core, and it is the one place this correction makes
the architecture cleaner rather than merely wider.

The consequence a future refinement or hydrology policy must live with: it is
handed a **mask** by its caller and is never compiled against a vocabulary. A
policy that wants "refine harder near rivers" takes an `EdgeProperties` mask as
a parameter. That is the cost of the firewall and it is the right cost.

### Union, not priority

`|` is commutative, associative and idempotent, and the merge is a reduce over
an **unordered** set of contributing chains — the noder's broad phase visits
buckets in whatever order it visits them, and the answer may not depend on
that. A priority scheme ("highest-ranked contributor wins") has none of the
three and needs a tie-break whose only honest source is input order, which is
exactly what node ids are sorted to avoid.

Precedence is a *drawing* question, not a *data* question, and it is answered
in the renderer section below — at the one layer where a human is looking.

### 32 is a ceiling, named rather than discovered

One word per edge, trivially parallel-reducible, no allocation. The vocabulary
it must hold is **linear** features — river, road, railway, coastline, contour,
wall, ditch — which is under ten. Land cover stays face-based (the sound half
of the claim `parallel_refinement.md` corrects) and its legacy vocabulary alone
is 23 classes, measured above. If a 33rd linear feature ever arrives, widening
to 64 is a one-line change to a type nobody pattern-matches on, **because
nobody can**: it has no named members.

### Rejected, recorded so it is not re-proposed

A per-edge index into a caller-supplied attribute table. It is genuinely open
in cardinality, and it makes the merge **allocate** — the union of two table
rows is a new row — which turns a parallel reduce into a synchronised one and
puts an indirection on the hottest array in the noder's output. An open
vocabulary is worth having; an open cardinality is not, at that cost, on
evidence of ten.

## `Chain::properties`, and the third parameter of `add_chain`

### The field

```cpp
struct Chain {
    std::uint32_t   begin{};
    std::uint32_t   count{};
    ChainRole       role{ChainRole::Breakline};
    EdgeProperties  properties{};      // was: bool is_river{false}
};
```

**Named `properties`, not `edge_properties`.** A `Chain` has no edges of its
own that outlive it — the noder's output edges do. The relation is *inheritance
downward*: every output edge a chain contributes geometry to receives that
chain's set, and `NodedPslg::edge_properties()` is the union of the sets so
inherited. Spelling the chain's field `edge_properties` would suggest the two
arrays are the same object; they are not, and the noder's guarantee 15 is
exactly the statement of how one becomes the other.

The existing comment at `include/terrain/core/pslg.hpp:77-79` — "one bit per
chain, permitted on every role … and never validated, because it is data, not
structure" — survives its widening intact and is the one sentence in this
increment that does not change meaning. Its "permitted on every role" clause
matters more now than it did: a wide river or a lake is legitimately an `Outer`
or `Hole` ring, and so is a walled enclosure.

### The parameter

`PslgBuilder::add_chain`'s third parameter is a **bound positional** on the
merged tree: `bindings/core.cpp:442` passes it positionally,
`tests/cpp/support/pslg_cases.hpp:252` and
`tests/cpp/property/prop_pslg_invariants.cpp:456` pass it positionally, and
both overloads default it. **Every line number in this section and in "How the
red step handles two mutation-round suites" below is read against `c3210c1`**,
the merged tree this increment measured, and this PR is what moves them;
`tools/check_citations.py` flags them as at-risk for exactly that reason, and
they resolve — `git show c3210c1:tests/cpp/support/pslg_cases.hpp | sed -n
'216p;252p'` prints `bool is_river{false};` and the `add_chain` forwarding
line.

**Ruling: it stays third, stays positional, keeps its default, and changes
type.**

```cpp
PslgBuilder& add_chain(std::span<const std::uint32_t>, ChainRole,
                       EdgeProperties = {});
PslgBuilder& add_chain(std::span<const Point2>, ChainRole,
                       EdgeProperties = {});
```

The important property is what this does to the 14 call sites, and it is the
whole reason the blast radius is safe:

- **`add_chain(idx, role, true)` stops compiling.** `EdgeProperties` has a
  defaulted default constructor and a named static factory, and **no
  constructor from `bool`, `int` or `std::uint32_t`**. There is no implicit
  conversion and there must never be one. Every one of the 14 files fails
  loudly at the compiler or at `mypy`; none of them silently keeps compiling
  with a changed meaning.
- **`add_chain(idx, role)` keeps compiling and keeps its meaning.** The empty
  set is what `is_river = false` meant. The large majority of call sites —
  every generated hole and every generated outer ring in
  `pslg_cases.hpp:valid_chain_specs` — are in this class and change by a
  rename only.

A converting constructor would be a convenience worth exactly one line of typing
and would cost the property that makes the rename mechanical-and-safe rather
than mechanical-and-silent. It is rejected here so it is not added later as a
tidy-up.

### Rejected: a compatibility shim

Keeping `bool is_river() const { return properties_.contains(EdgeProperties::bit(0)); }`
as a deprecated accessor would let the 14 files migrate one at a time. It is
rejected on two grounds, either sufficient: it puts the name "river" in
`terrain::` **permanently**, which is the single thing this increment exists to
remove; and it makes every call site keep compiling, which destroys the loud
failure the previous section relies on.

## The boundary model

New module, `src_python/tin_engine/features.py`. It **imports nothing
first-party and never imports `_core`**, so it is constructible and testable
with no compiled extension in the process. **Corrected in this PR:** this
paragraph originally read "which is what lets `viz/` depend on it". That is
unsupportable and unnecessary. `test_viz_svg.py::TestModuleIsolation` pins that
`style.py` and `fixtures.py` import **no** first-party module at all and that
`svg.py`'s first-party imports are confined to `{.scene, .style, .protocols}`,
so no `viz/` module may import `features` — and under the renderer ruling below
none needs to: `svg.py` takes its precedence from `SvgStyle` and `cli.py`, the
composition root, is the one module that names both.

```python
class EdgeProperty(BaseModel):          # frozen, extra="forbid"
    name: Annotated[str, Field(pattern=r"^[a-z][a-z0-9_]*$")]
    bit:  Annotated[int, Field(ge=0, lt=32)]

class EdgeVocabulary(BaseModel):        # frozen, extra="forbid"
    properties: tuple[EdgeProperty, ...]

    # model_validator(mode="after"): names unique, bits unique.

    def mask(self, *names: str) -> int: ...
    def names(self, mask: int) -> tuple[str, ...]: ...
    def fingerprint(self) -> str: ...

DEFAULT_VOCABULARY = EdgeVocabulary(properties=(
    EdgeProperty(name="river",     bit=0),
    EdgeProperty(name="road",      bit=1),
    EdgeProperty(name="railway",   bit=2),
    EdgeProperty(name="coastline", bit=3),
    EdgeProperty(name="contour",   bit=4),
    EdgeProperty(name="wall",      bit=5),
    EdgeProperty(name="ditch",     bit=6),
))
```

### What it validates, and why each check is there

- **`bit` in `[0, 32)`.** The C++ ceiling, restated at the only layer that sees
  untrusted input.
- **Bits unique, names unique.** Two names on one bit is the defect risk 20 is
  about, in its cheapest and most detectable form. It is refused at
  construction, which is the only moment it is detectable for free.
- **`name` matches `^[a-z][a-z0-9_]*$`.** This is **not** style. On branch
  `increment6b-ii-renderer-red`, `svg.py`'s `_edge_classes` interpolates its tokens
  into a `class="…"` attribute **unescaped** — line 154,
  `f"role-{_role_name(edge.role)}"`; only `_text` escapes, at line 162. Run
  `git show increment6b-ii-renderer-red:src_python/tin_engine/viz/svg.py`.
  Today every token comes from an enum name and is safe by construction. A
  vocabulary read from a configuration file is not, and a name containing a
  quote ends the attribute. The pattern is that hole's only closure, and it is
  cheaper than escaping because a CSS class token has no legitimate use for any
  character the pattern excludes.

### How a producer and a consumer come to agree which bit means "river"

This is risk 20 (`05b-noder-driver.md`), inherited here with the obligation to
mitigate it. Nothing in C++ can check it: the core merges opaque bits, so two
producers that disagree about bit meanings produce a wrong mesh no C++ suite
can see. Three mechanisms, in descending strength:

1. **A bare mask is unreachable through the Python API.** `mask()` and
   `names()` are *methods on a vocabulary*, not module functions. You cannot
   obtain a mask without holding the vocabulary that produced it, and you
   cannot interpret one without holding a vocabulary you chose. This makes the
   dangerous state — a number with no units — costly to construct rather than
   the default.
2. **`fingerprint()`, a stable digest over the sorted `(bit, name)` pairs**, is
   a field of whatever serialized artifact carries a mesh or a noded PSLG,
   **exactly as CRS is**. An artifact read back with a fingerprint that does not
   match the vocabulary in hand is refused. The vocabulary is therefore *data
   travelling with the mesh*, not a constant in the source — which is the same
   ruling `project_structure.md:185` already makes about CRS, for the same
   reason, and the analogy is the argument.
3. **`names(mask)` raises on a bit no property names**, rather than dropping it.
   Dropping is the silent loss this entire increment exists to prevent, and
   this is the only moment the system can see it. Reaching this state means the
   mask came from a different vocabulary than the one in hand — risk 20's
   failure, caught.

**Residual, stated because it is real:** nothing stops a caller passing a raw
integer through pybind by hand. The mitigation makes the wrong thing awkward;
it does not make it impossible, and no C++ suite will ever see it. That is the
price of keeping feature names out of `terrain::`, and it is the right price.

### The one asymmetry, called out so it does not read as an inconsistency

The **data** path may not lose a property: `names()` raises. The **drawing**
path may decline to draw one: a property with no stroke in the style simply
contributes no token, and the edge is still drawn as the constrained breakline
it is. The failure modes are different — a lost property is a wrong answer, an
undrawn property is a picture with less in it — and `06-cdt-viewer.md`'s ruling
that "testing effort goes where the failure mode is a *right*-looking wrong
picture" is exactly this distinction.

## The C++ / Python boundary: `EdgeProperties` is not bound

**Ruling: `bindings/core.cpp` does not expose a `py::class_<EdgeProperties>`.**
It exposes the raw word.

- `Chain.properties` is an `int` (from `EdgeProperties::bits()`).
- `build_pslg`'s chain tuple third element is an `int` mask, cast to
  `EdgeProperties` on the C++ side.

Reason: Python already has an integer with `|`, `&` and `bit_count()`, and the
*semantics* live in `EdgeVocabulary`. A bound class would add a second,
competing vocabulary object next to the Pydantic one and would carry no meaning
the `int` does not. The pybind layer is a firewall, not a place to grow an API.

The cost is that an `int` at the boundary is unvalidated, so:

**`build_pslg` raises `ValueError` on a mask that is negative or has a bit at or
above 32**, naming the chain index. Not a `PslgDiagnostic`: `PslgError`
enumerates structural defects of a constraint set, every one of which a C++
caller can also commit, and this one no C++ caller can — `EdgeProperties::bit`
is the only route to a set bit and `i >= kMaxProperties` is a precondition
violation, not a datum. Putting a marshalling error into that enum would put a
Python type error inside a C++ vocabulary. It joins the mis-shaped vertex array
(`bindings/core.cpp:452`) as a marshalling `ValueError`, which is the existing
precedent and the right one.

**`EdgeProperties::bit(i)` has a narrow contract**: precondition
`i < kMaxProperties`, debug-asserted, unchecked in release. This is exactly
`SnapGrid`'s shape — `can_snap` is "the driver's one-pass admission check … and
the only place a coordinate is compared against `kMaxGridIndex`. After that
pass `snap()` is branch-free" (`include/terrain/core/snap_grid.hpp:94-96`), and
`snap` asserts (`:114`). Admission happens once, where untrusted data arrives;
everything after it is branch-free. The precedent is cited because it means this
increment is not inventing a policy.

## The renderer: the ruling 5b named and deferred

5b deferred this, proposing "a declared precedence in the style model with the
set still carried on the edge" and leaving it to "whoever owns `style.py`".
**Ruled here, in this increment's PR.** Two of 5b's premises about the renderer
are wrong and the correction changes the answer's shape; see "What I now find
false".

### The finding stands: two strokes on one line is unreadable

An edge that is both road and river cannot carry both stroke classes. At
gallery scale, two overlaid strokes on one polyline read as a rendering defect.

### The ruling

- **`SceneEdge.is_river: bool` becomes `SceneEdge.properties: int`**, carrying
  the **full** set. `scene.py` never ranks and never drops — it unions, exactly
  as `_chain_edges` already ORs, and for the same reason it types `role` as
  `object`: the scene builder is forbidden from knowing a vocabulary.
- **`svg.py`'s `_edge_classes` emits at most one property token**, chosen by the
  first match in an ordered list.
- **That list is a field of `SvgStyle`**, `style.py`, which is already "frozen
  Pydantic V2" and already the renderer's declarative input:

  ```python
  class PropertyStroke(BaseModel):    # frozen, extra="forbid"
      bit:   Annotated[int, Field(ge=0, lt=32)]
      token: Annotated[str, Field(pattern=r"^[a-z][a-z0-9_-]*$")]

  # on SvgStyle:
  property_strokes: tuple[PropertyStroke, ...] = ()
  ```

  In **descending draw priority**, first match wins. Priority is explicit and
  independent of bit numbering — a vocabulary is free to number by feature
  family while the stylesheet draws water over infrastructure.
- **The default is `()`, and `cli.py` supplies the gallery's list.** The
  precedent is exact and already in the tree: `scene.py`'s `closed_roles` is
  supplied by `cli.py` for precisely this reason — "``closed_roles`` is
  therefore supplied by ``cli.py``, the composition root that knows the enum"
  (`src_python/tin_engine/viz/scene.py:141-143`). `style.py`'s own docstring
  draws the same line: it models "the arithmetic the viewport divides by", and
  colours are "taste" that lives in `svg.py`'s CSS. A default naming `river`
  would be policy in the module that declares it holds none.

### One stale half of a comment, recorded rather than edited

`tests/python/test_viz_svg.py:132-136`, the comment above `GALLERY_STROKES`
(`:137`), claims both "Water over infrastructure" and "deliberately NOT in bit
order". With
`RIVER_BIT = 0` those cannot both hold: water-first *is* bit order here. **The
list itself is authoritative and correct** — it is road-first, which is what
keeps `test_the_declared_order_is_preserved` able to catch a model that sorts by
bit — so only the phrase is stale. It is recorded here rather than fixed because
commit 4 is `@developer`'s and touches no test file; the edit is one line and is
`@tester`'s to make.

### Why this is 7's and not 6's

The set on the edge and the single token in the document are the same decision
seen from two sides, and splitting them would ship a `SceneEdge` carrying a set
that the renderer reads as a bool for one PR. It is ~12 lines of the estimate
below. It lands in **this PR**.

## Sequencing

**6b-ii merges → 7 → 5b → 5c.**

- **7 before 5b** is the user's ruling. 5b's own recommendation agreed, for the
  reason that the noder is the first consumer that genuinely needs union
  semantics; `05b-noder-driver.md` is written so either order works, and this
  one saves 5b from shipping a one-member set it would immediately widen.
- **7 after `increment6b-ii-renderer-red`** is this document's ruling, on
  measured cost. That branch carries `is_river` in `svg.py`, `fixtures.py`,
  `cli.py`, `test_viz_svg.py` (24 occurrences, the largest single churn file in
  the whole change) and `test_cli_draw.py`. Landing 7 first means paying that
  conflict twice — once in 7 against an unmerged branch, and again when 6b-ii
  rebases through a type change. It also means ruling on `_edge_classes` and
  `SvgStyle.property_strokes` against files that do not exist on master.

  **Risk if violated:** 7's renderer section becomes unimplementable and must be
  deferred to a follow-up, at which point `SceneEdge.properties` ships as a set
  the renderer reads as a bool — the exact split the previous section rejects.

## The blast radius, measured

```sh
grep -rln is_river --include='*.hpp' --include='*.cpp' --include='*.py' \
     --include='*.pyi' . | grep -v '^./legacy'
```

returns **14** files on `increment5b-design` (identical to master for these
paths). With per-file occurrence counts, because the file list alone is what
made this change look uniformly expensive:

| File | hits | Production? | Invariant-critical? |
|---|---|---|---|
| `include/terrain/core/pslg.hpp` | 2 | yes | — |
| `include/terrain/core/pslg_builder.hpp` | 4 | yes | — |
| `bindings/core.cpp` | 4 | yes | — |
| `src_python/tin_engine/_core.pyi` | 1 | yes | — |
| `src_python/tin_engine/viz/protocols.py` | 1 | yes | — |
| `src_python/tin_engine/viz/scene.py` | 6 | yes | — |
| `tests/cpp/support/pslg_cases.hpp` | 2 | no | support header |
| `tests/cpp/support/cdt_cases.hpp` | 1 | no | support header |
| `tests/cpp/property/prop_pslg_invariants.cpp` | **1** | no | **yes** |
| `tests/cpp/unit/test_pslg.cpp` | 4 | no | no |
| `tests/cpp/unit/test_pslg_builder.cpp` | 6 | no | no |
| `tests/python/test_core_cdt.py` | 2 | no | no |
| `tests/python/test_viz_protocols.py` | 1 | no | no |
| `tests/python/test_viz_scene.py` | **10** | no | **yes** |

Plus, on `increment6b-ii-renderer-red` and therefore on master once it merges
(`git grep -c -i river increment6b-ii-renderer-red`):
`src_python/tin_engine/viz/svg.py` (6), `src_python/tin_engine/viz/fixtures.py`
(8), `src_python/tin_engine/cli.py` (1), `tests/python/test_viz_svg.py` (24),
`tests/python/test_cli_draw.py` (2).

## How the red step handles two mutation-round suites

This is the question the increment exists to answer before anyone writes a
line, because it is where the round's cost actually is.

### First: one of the two is not affected

`tests/cpp/property/prop_pslg_invariants.cpp` contains **exactly one**
`is_river` token, at line 456 of the pre-rename file:

```cpp
b.add_chain(std::span<const std::uint32_t>{idx}, s.role, s.is_river);
```

**Corrected in this PR: one token is not one site, and the cost was
understated.** `grep -c is_river` prints 1 because the other three sites pass
the field *positionally*, without naming it —
`git show c3210c1:tests/cpp/property/prop_pslg_invariants.cpp | grep -n
'ChainRole::[A-Za-z]*, false}'` prints lines 500, 506 and 511, three
`ChainSpec{…, ChainRole::X, false}` aggregate initialisers. Four sites, not
one. A token count read as a site count, which is the same error this document
records at "What I now find false" item 2 for `fixtures.py`.

The **ruling** the figure supports is unaffected and stands: it is a
**pass-through** of `ChainSpec::is_river` (`pslg_cases.hpp:216`) into the
builder, the three aggregate initialisers are the empty set spelled positionally,
**no assertion in that file mentions it**, and no mutant in increment 3's round
concerns it — that round is about vertices, windings and diagnostic sets.
`@tester` reconfirmed that independently in the red step.

**Ruling: `prop_pslg_invariants.cpp` gets no mutation round in this
increment.** A mutation round exists to prove that a suite's *assertions* kill
mutants; this change alters no assertion in that file. What it alters is a
field name in a shared support header, one forwarding expression and three
positional aggregate initialisers. Recording
this as a ruling rather than as an omission is the point: `05b-noder-driver.md`
counted two invariant-critical suites by file list, which is true and, as a cost
estimate, misleading by roughly the whole of the C++ half.

### Second: the affected one gets *cheaper*, and that is measurable

`tests/python/test_viz_scene.py`'s affected surface is bounded and enumerable:
the `FakeChain` field (`:101`), `CHAIN_MEMBERS` (`:216`), the fake builder
(`:255-256`), and five assertions (`:481`, `:484`, `:491`, `:528`, and the
comment at `:598`).

The load-bearing one is
`test_the_river_bit_is_the_or_over_every_chain_on_the_edge` (`:508-528`), whose
own comment records what one bit cost it:

> Both orders, because a single order is passed by an implementation that
> simply takes the last chain's bit — **measured: with only the river-last
> case, the overwrite mutant survived.**

Under one bit, `{river} | {river} == {river}`, so "union" and "take the last
contributor" are indistinguishable except by choosing fixture orderings that
make them differ, and the suite needed a `parametrize` over both. **Over sets
they come apart**: two chains carrying `bit(0)` and `bit(1)` produce
`{road, river}`, which equals neither operand, so **one** fixture in **one**
order kills the overwrite mutant. The `parametrize` collapses.

This is the same improvement `05b-noder-driver.md` records for the noder's
mutant 7 — the widening buys an assertion where the one-bit form needed a
fixture chosen to make the difference visible — arriving here first because
this is the increment that does the widening.

### Third: the commit shape, which is the actual answer

A rename does not "fail" the way a missing feature fails; it fails to compile.
So the round is **four commits in two red/green pairs**, and the pairing is what
keeps the invariant-critical files to a single, reviewable touch.

1. **Red (`@tester`) — the new thing only.**
   `tests/cpp/unit/test_edge_properties.cpp` and
   `tests/python/test_features.py`. Both are complete suites for units that do
   not exist: the C++ one fails to compile, the Python one fails at import.
   **This commit touches no existing file, test or production.**
2. **Green (`@developer`).** `include/terrain/core/edge_properties.hpp` and
   `src_python/tin_engine/features.py`. Both new files. Nothing else changes;
   the build and every existing suite stay green.
3. **Red (`@tester`) — the migration.** The rename across every test and support
   file, *plus the assertions the new width makes possible*: a chain carrying
   two disjoint properties in `test_pslg_builder.cpp`, a two-disjoint-property
   edge in `test_viz_scene.py`, the out-of-range mask `ValueError` in
   `test_core_cdt.py`, the single-token ruling in `test_viz_svg.py`. Red for a
   mechanical reason — it names `Chain::properties`, which does not exist — and
   visibly so.
4. **Green (`@developer`) — the migration.** `pslg.hpp`, `pslg_builder.hpp`,
   `bindings/core.cpp`, `_core.pyi`, `protocols.py`, `scene.py`, `style.py`,
   `svg.py`, `fixtures.py`, `cli.py`. **No test file.**

**The mechanism, and it is checkable rather than aspirational:**
`git show --stat` on commit 3 must name **only** paths under `tests/`, and on
commit 4 **no** path under `tests/`. That is `docs/increments/README.md`'s rule
("`@developer` does not edit tests, and no test change hides inside an
implementation commit") applied to a change whose whole character is a rename,
and it is the property that makes a 14-file diff auditable in two `git show`s.

Pairs 1–2 and 3–4 are **not** independent and must not be parallelised: commit
3 does not compile without commit 2's header.

### Fourth: what makes the migration safe rather than merely mechanical

The no-implicit-conversion ruling above. There is no call site in the 14 files
that can silently keep compiling with a changed meaning: either it passed a
`bool` literal and now fails, or it passed nothing and is unchanged in meaning.
`mypy --strict` gives the Python half the same property. **The compiler is the
migration's test**, which is why this round does not need a mutation budget
proportional to its file count.

## Gates

Per `.claude/REQUIRED-READING.md` and `05b-noder-driver.md`'s corrected
reasoning: **a gate's moment is before the red commit, because the red step
fixes the PR's shape, and a gate whose remedy costs more than the split it buys
is argued away rather than obeyed.** Three seams have been armed in this project
and one was obeyed; that record is an argument against arming many. **Two are
armed here, both cheap at the moment they fire.**

**Gate A — the shape gate, `@tester`'s, before *each* red commit.** Not a line
count. Commit 1 names only new test files; commit 3 names only paths under
`tests/`. `@tester` checks `git diff --cached --name-only` before committing and
raises a deviation rather than committing. The remedy is to move a file between
the two commits, which costs nothing, which is the whole reason this gate can be
obeyed at the moment it fires.

**Gate B — the sequencing gate, before commit 1.** If
`increment6b-ii-renderer-red` has not merged, the renderer section is
unimplementable. The remedy is to stop and ask, not to defer the section: see
"Sequencing", where deferring it is rejected by name.

**Not armed: a LOC gate.** The estimate is ~185 against 700 and the largest
single production row is 45. A threshold with 500 lines of headroom is a
trip-wire nobody can trip, and arming it would dilute the two above. Recorded
so the omission reads as a decision.

## Files and LOC

Instruments, per `05b-noder-driver.md` (the repo has two and they disagree by
15% on C++): C++ `grep -vcE '^\s*(//|$)'`, Python `grep -vcE '^\s*(#|$)'`.

| File | Contents | Est. |
|---|---|---|
| `include/terrain/core/edge_properties.hpp` | `EdgeProperties`, its six members, `kMaxProperties` | ~45 |
| `include/terrain/core/pslg.hpp` | field rename and retype; comment widened | ~3 |
| `include/terrain/core/pslg_builder.hpp` | two `add_chain` signatures, two aggregate inits | ~4 |
| `bindings/core.cpp` | `properties` readonly as `int`, tuple cast, mask range check + `ValueError`, docstrings | ~30 |
| `src_python/tin_engine/_core.pyi` | `Chain.properties: int`, `build_pslg` chain tuple | ~3 |
| `src_python/tin_engine/features.py` | `EdgeProperty`, `EdgeVocabulary`, `mask`, `names`, `fingerprint`, `DEFAULT_VOCABULARY` | ~55 |
| `src_python/tin_engine/viz/protocols.py` | `ChainLike.properties: int` | ~2 |
| `src_python/tin_engine/viz/scene.py` | `SceneEdge.properties: int`; `_chain_edges` union | ~10 |
| `src_python/tin_engine/viz/style.py` † | `PropertyStroke`, `SvgStyle.property_strokes` | ~16 |
| `src_python/tin_engine/viz/svg.py` † | `_edge_classes` single-token selection | ~8 |
| `src_python/tin_engine/viz/fixtures.py` † | the `RIVER` row's mask; `FixtureChain` field | ~4 |
| `src_python/tin_engine/cli.py` † | pass `property_strokes`; chain tuple | ~5 |

† arrives with `increment6b-ii-renderer-red`; see "Sequencing".

**~185 non-comment production lines.** No `src/` file: `edge_properties.hpp` is
`constexpr` throughout and header-only.

**This is higher than `05b-noder-driver.md`'s ~130, and the difference is
identified rather than absorbed.** That estimate omits `features.py` — the
Pydantic model it names as risk 20's only mitigation — and omits `style.py`,
`cli.py` and the `bindings/core.cpp` mask range check. Its `bindings/core.cpp`
row is also ~25 against a measured project factor of two on every binding
estimate so far (6a: +283 against ~150, +125 against ~60); ~30 here is the
sketch, and the honest statement is that this is the one row that has doubled
twice before. If it doubles again the total is ~215, still a third of the
ceiling, which is why no LOC gate is armed.

**Test churn, counted separately because tests are excluded from the ceiling and
are the expensive part of this round**: 8 existing test files on master, 5 more
once 6b-ii merges, 2 new suites. The largest single file is `test_viz_svg.py` at
24 occurrences, and it is **not** invariant-critical (`06-cdt-viewer.md:672-676`).

## What is worth testing

| Unit | Suite | Invariant-critical? |
|---|---|---|
| `core/edge_properties.hpp` | `tests/cpp/unit/test_edge_properties.cpp` | **no** |
| `features.py` | `tests/python/test_features.py` | **yes**, mutation round |
| the widened scene join | `tests/python/test_viz_scene.py` (existing) | **yes**, round re-run |
| the widened chain field | `tests/cpp/property/prop_pslg_invariants.cpp` | **no** — see above |
| the single-token stroke | `tests/python/test_viz_svg.py` (existing) | no |

**`test_edge_properties.cpp` — no mutation round, and the reason is the
README's rule rather than an omission.** It is a value type with no kernel
parameter and no exactness claim: `bit`, `empty`, `contains`, `bits`, `|`, `==`,
plus `static_assert`s that the whole surface is usable in a constant expression.
Its failure mode is a typo and it is loud. `05-noder.md` excluded `node_set.hpp`
in writing for the same reason and that is the precedent. The one property worth
asserting rather than enumerating: `|` is commutative, associative and
idempotent over generated masks, because those three are the entire load-bearing
claim and they are one `REQUIRE` each.

**`test_features.py` — invariant-critical, mutation round.** This is the only
place in the system where a wrong answer is **silent**: the C++ merges opaque
bits and a mis-mapped vocabulary produces a wrong mesh no C++ suite can see
(risk 20). Everything else in this increment fails at a compiler.

### Mutants the round must kill

1. **`mask()` shifting by the tuple index instead of by `EdgeProperty.bit`.**
   Killed by a vocabulary whose bits are not `0..n-1` — the one-property
   vocabulary `(name="road", bit=3)` is the cheapest.
2. **The bit-uniqueness validator dropped.** Killed by a two-names-one-bit
   vocabulary that must be refused at construction.
3. **The name-uniqueness validator dropped.** Killed by two properties sharing
   a name on different bits, where `mask("river")` becomes order-dependent.
4. **`names(mask)` silently dropping a bit no property names.** Killed by
   asserting the raise, naming the bit. This is the mutant the ruling exists
   for and the one whose survival is indistinguishable from correct behaviour
   in any picture.
5. **The name pattern relaxed to any `str`.** Killed by two refusals: a name
   containing `"` and a name containing a space — both of which
   `src_python/tin_engine/viz/svg.py:_edge_classes` would interpolate into a `class` attribute
   unescaped.
6. **`fingerprint()` hashing names only, or bits only.** A digest that cannot
   see which bit a name occupies is a digest that blesses exactly risk 20's
   failure, so this is the mutation round's reason for existing. Two cases are
   needed, and **the swap this entry originally named is not one of them**:

   - *Bits only* is killed by a rename on unchanged bits — `{river: 0}` against
     `{creek: 0}`.
   - *Names only* is **not** killed by "two vocabularies over the same names
     with two bits swapped", which is what this entry claimed when the design
     was written. Measured in the red step: a names-only digest emitted in
     sorted-`(bit, name)` order prints `river;road` for `{river: 0, road: 1}`
     and `road;river` for `{river: 1, road: 0}`, so the swap already
     fingerprints differently and the mutant survived all 72 cases of the
     suite's first draft. The case that separates them is two vocabularies
     agreeing on **every name and its order** while disagreeing about which bit
     one of them occupies — `{river: 0, road: 1}` against `{river: 0, road: 2}`,
     whose names-only digest is `river;road` both times. That is risk 20
     exactly: an artifact written by one is read by the other with a property
     shifted.

   Both cases are in `tests/python/test_features.py`
   (`test_fingerprint_sees_a_permutation_of_the_bits`, kept because it is a real
   property, and `test_fingerprint_sees_a_bit_moved_without_reordering_the_names`,
   which is the one that kills the mutant).

   Related, and not in the original entry: a `hash()`-based digest is killed
   only by a cross-process check. `PYTHONHASHSEED` salting makes an in-process
   digest look perfectly stable while refusing every artifact ever written, so
   the round needs `test_fingerprint_is_stable_across_processes`, which compares
   two children run under different seeds.
7. **`mask()` accepting an unknown name and returning 0.** Killed by asserting
   the raise. Silent zero is a constraint that quietly loses every property.

And in the re-run scene round:

8. **The union replaced by "first contributor wins" or "last contributor
   wins".** Killed by one two-disjoint-property edge, in one order — see "the
   affected one gets cheaper".
9. **`SceneEdge.properties` taking the first chain's set rather than the union
   across chains**, which is the role join's rule applied to the wrong field.
   Killed by the same fixture; the two fields deliberately have *different*
   merge rules (role: first wins; properties: union) and that difference is now
   assertable where under one bit it was not.

### Template spend

**Zero `TEMPLATE_TEST_CASE`, and there is nothing to arm.** `EdgeProperties` has
no kernel parameter, no template parameter of any kind, and no exactness claim.
`Chain` and `PslgBuilder` are already covered by increment 3's cross product and
this increment adds no kernel-dependent path to either. The README's rule is
that a cross product must prove something a later increment depends on; none
would.

## Degeneracy and failure policy

| Condition | Where | Response |
|---|---|---|
| Empty property set | everywhere | **Legal, and the default.** An unclassified constraint. Not an error, not a diagnostic, not a warning. |
| `EdgeProperties::bit(i)`, `i >= 32` | C++ | Precondition, debug `assert`, unchecked in release. `SnapGrid::snap`'s shape. |
| Mask with a bit `>= 32`, or negative, from Python | `build_pslg` | `ValueError`, naming the chain index. Marshalling, not a `PslgDiagnostic`. |
| Two vocabulary properties on one bit | `EdgeVocabulary` | Pydantic `ValidationError` at construction. |
| Two vocabulary properties with one name | `EdgeVocabulary` | Pydantic `ValidationError` at construction. |
| Name outside `^[a-z][a-z0-9_]*$` | `EdgeProperty` | Pydantic `ValidationError`. |
| `mask()` given an unknown name | `EdgeVocabulary` | Raises, naming the name. Never a silent 0. |
| `names()` given a bit no property names | `EdgeVocabulary` | Raises, naming the bit. Never a silent drop. |
| Artifact fingerprint ≠ vocabulary in hand | the reader | Refuse. Same policy as a CRS mismatch. |
| Edge with a property that has no stroke | `svg.py` | **Draw no extra token.** Deliberately *not* an error — see "the one asymmetry". |
| Edge with several properties that all have strokes | `svg.py` | Highest-priority token only, by `SvgStyle.property_strokes` order. |
| A chain's set on a role with no edges (1-vertex chain) | validator | Untouched. `is_river` was "never validated, because it is data, not structure" (`pslg.hpp:77-79`) and that survives. |

**Never silently.** Every row above either succeeds, or refuses with the value
that did not fit. `04-cdt.md`'s ruling ("never silently produces a mesh that
ignores a crossing") is this project's general form and this table is its
instance for semantics rather than geometry.

## Risks

1. **Risk 20 is inherited and cannot be closed, only narrowed.** Nothing in C++
   can check that two producers agree about bit meanings. Mitigated by the three
   mechanisms above, of which only the third (`names()` raising) actually
   *detects* a disagreement, and only when the two vocabularies differ in
   coverage rather than in assignment. **Two vocabularies that both name bit 3,
   differently, are undetectable except by fingerprint** — which is why the
   fingerprint is the mechanism that must travel with the data, and why mutant 6
   is in the round.
2. **`DEFAULT_VOCABULARY` becomes a de-facto standard nobody declared.** Bit 0
   is `river` so that today's one-bit data and the `river` gallery fixture keep
   meaning what they meant; that is a migration convenience and it will be read
   as a specification. Mitigated by pinning it with a test that asserts
   `DEFAULT_VOCABULARY.mask("river") == 1` *and* a comment saying the constant
   is a default and not a schema. Residual: real.
3. **The renderer's precedence is invisible in the picture.** An edge that is
   both road and river draws as a river and nothing says so. A reader concludes
   there is no road. Not mitigated in this increment — the honest fixes are a
   legend or a dashed overlay, both of which are 6's territory and neither of
   which is free. Named so it is not discovered. The counterweight is that the
   *data* is intact and queryable; only the drawing is lossy, which is the
   asymmetry ruled on above.
4. **The migration crosses a branch boundary.** Five of the affected files exist
   only on `increment6b-ii-renderer-red`. Mitigated by Gate B. Residual: if
   6b-ii's own review round changes `_edge_classes` or `fixtures.py`, this
   design's line estimates for those rows go stale and must be re-read, not
   re-used.
5. **`EdgeProperties` ships with no C++ consumer that merges anything.** Until
   5b, nothing in `terrain::` calls `operator|` outside its own suite; `Chain`
   merely carries the field. That is the third increment in a row shipping
   ahead of its caller (`05-noder.md` risk 6, `05b-noder-driver.md` risk 18),
   and it is the direct consequence of the ceiling being per PR. Discharged at
   5b.
6. **`int` at the pybind boundary is a type with no units.** `Chain.properties`
   is an `int` and so is a mis-typed vertex count. Mitigated only by the mask
   range check and by `mask()`/`names()` being the sole sanctioned converters.
   A `py::class_` would not have fixed this — it would have moved the problem
   into a second vocabulary object — which is the argument for the ruling, not
   a defence of its cost.

## Documentation this PR fixes

Per `docs/increments/README.md`: fixed here, or not recorded. There is no
ledger, and this list is closed.

- **`parallel_refinement.md`** — the corrected "Edge metadata" section is sound
  and is not rewritten. One line changes: `:190-192` points at
  `05b-noder-driver.md` for "the ruling and … which increment owns the type";
  the owner is this file.
- **`auto_catchments.md`** — `:12`, `:106`, `:107`, `:154` tag river polylines
  `is_river = true` and describe lakes as having "no `is_river` flag". Restated
  over property sets, with the Strahler-cutoff paragraph (`:154`) kept, since
  choosing when a creek becomes a river is a real caller decision that survives
  the widening unchanged.
- **`project_structure.md`** — `:266`'s `mesh` section says "the per-edge
  feature property set (increment 7; one bit, `is_river`, until then)"; the
  parenthesis goes. The directory listing at `:64-78` gains
  `src_python/tin_engine/features.py`, and the `core` module description gains
  `edge_properties.hpp`.
- **`docs/increments/02-core-geometry.md:101-102`** — "the noder's `is_river`
  merge cares" about `operator==` being ordered. The claim about ordering is
  true and stays; the spelling is stale.
- **`docs/increments/03-pslg.md`** — `:47`, `:65`, `:172-173` and `:593` carry
  the `bool is_river` declaration, its paragraph and the two `add_chain`
  signatures. Widened in place, with a note that increment 3 shipped the one-bit
  form and 7 replaced it — `03-pslg.md` remains the record of what was true at
  increment 3 and says so, following `05b-noder-driver.md`'s treatment of
  `04-cdt.md`'s mapping table.
- **`docs/increments/06-cdt-viewer.md`** — `:153`, `:187-188`, `:336`, `:444`
  and `:457` describe the chain record and the `river` fixture row over one bit.
  Widened, and `:457`'s row gains the two-property fixture the renderer ruling
  requires.
- **`docs/increments/05b-noder-driver.md`** — three corrections, all measured
  in "What I now find false" below: the `style.py` attribution, the
  `fixtures.py` row count, and the increment-7 LOC row. Its "Edge properties"
  section otherwise stands and is *not* rewritten; a pointer to this file
  replaces its "for increment 7 to rule on" framing, since the ruling now
  exists.
- **`docs/increments/05-noder.md`** — `:781`, `:840-853`, `:890-895`, `:926`,
  `:972-1003` are already marked by `05b-noder-driver.md` as superseded in
  width. **Not re-touched here**; that marking is 5b's PR's and duplicating it
  would be the ledger this project refuses to keep.
- **`testing.md`** — the `noding` invariant list's OR rule, if 5b's PR has not
  already replaced it by the time this lands. Checked at review; if 5b's text
  is in place it is already stated over sets and there is nothing to do.

## What I now find false in the documents this increment depends on

Each with the command that refutes it, run before being written down.

1. **`05b-noder-driver.md`: "whoever owns `style.py`" owns the river stroke.**
   `git grep -n -i river increment6b-ii-renderer-red -- src_python/tin_engine/viz/style.py`
   returns **nothing**. The stylesheet is `svg.py`'s `STYLESHEET`
   (`line.river { stroke: #0b8fb0; … }`, line 56 of
   `git show increment6b-ii-renderer-red:src_python/tin_engine/viz/svg.py`) and the
   class token is chosen in `_edge_classes`, lines 155-156 of the same file. `style.py` is `SvgStyle`, a frozen
   Pydantic model of canvas arithmetic whose own docstring says colours live in
   `svg.py` "because a class token on an element is structure … while the colour
   that token resolves to is taste". The ruling above *does* put the precedence
   list in `style.py` — but as structure, and for the opposite reason from the
   one implied.
2. **"All eight rows of `fixtures.py`."**
   `git grep -c -i river increment6b-ii-renderer-red -- src_python/tin_engine/viz/fixtures.py`
   prints `8`, but those are eight *occurrences*, of which exactly **one** is a
   gallery row: `RIVER`, at line 178 of
   `git show increment6b-ii-renderer-red:src_python/tin_engine/viz/fixtures.py`. The
   rest are the `FixtureChain` field (line 45), the builder loop (82-83), a
   comment (95), the fixture's description string (180) and the gallery tuple
   entry (218).
   `06-cdt-viewer.md:457` confirms one river row. A line count read as a row
   count.
3. **The 14-file list is complete for master and incomplete for 6b-ii.**
   `git grep -ln is_river increment6b-ii-renderer-red` also returns
   `src_python/tin_engine/cli.py` (`:95`) and `tests/python/test_cli_draw.py`,
   and `git grep -c -i river` on that branch returns **24** for
   `tests/python/test_viz_svg.py` — the largest single test-churn file in the
   change, absent from every list of it so far. It is not invariant-critical
   (`06-cdt-viewer.md:672-676`).
4. **"Two of which carry mutation rounds" is true by file list and misleading
   as a cost.** `grep -c is_river tests/cpp/property/prop_pslg_invariants.cpp`
   prints `1`, and that one line is a pass-through with no assertion attached.
   Ruled on above.
5. **`05b-noder-driver.md`'s ~130-line estimate for increment 7 omits its own
   mitigation.** Its row list is `edge_properties.hpp`, `pslg.hpp` +
   `pslg_builder.hpp`, `bindings/core.cpp`, `_core.pyi`, `protocols.py`,
   `scene.py`, and the three 6b-ii viz files. `features.py` — the Pydantic model
   its risk 20 names as the only mitigation — is not in it, nor are `style.py`,
   `cli.py` or the mask range check. ~185 is this document's figure.
6. **`docs/increments/02-core-geometry.md:101-102` still calls it "the noder's
   `is_river` merge".** Stale spelling; the substantive claim (`operator==` on
   `Segment2` is ordered and the merge depends on it) is true and stays.

**Nothing in `parallel_refinement.md`'s corrected "Edge metadata" section is
false.** It was re-read line by line against this design: the union semantics,
the commutativity/associativity/idempotence argument, the "geometrically
forgotten but not semantically forgotten" distinction, the area/linear split,
the 23-class land-cover figure (verified above at 23) and the "the C++ core does
not name the properties" ruling are all exactly what this increment implements.
Its only defect is the forward pointer in item 1 of the fix list.
