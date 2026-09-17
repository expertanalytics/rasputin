# Increment 4 — the constrained Delaunay triangulation

Status: design settled. One item the brief listed as already settled is
**overturned below on evidence from the vendored header** — the scratch index
buffer is deleted rather than protected — and one shipped function in
`core/pslg.hpp` loses its only reason to exist as a consequence.

Ships `include/terrain/core/indexed_mesh.hpp`,
`include/terrain/cdt/{result,constrained_edges,triangulate,detria_backend}.hpp`
and `src/cdt/detria_backend.cpp`. Namespaces `terrain` (the mesh) and
`terrain::cdt` (everything else). Depends on increment 3 (`Pslg`, `Chain`,
`ChainRole`), increment 2 (`Point2`, `Segment2`) and increment 1 only for the
test oracle, not for the production path.

Makes possible: **a real CDT over a polygon with holes** — the goal line of the
plan increment 1 was written against. Everything after this point (`mesh`,
`refinement`, `flip`) consumes triangles; nothing in the tree has produced one
yet.

**Invariant-critical suites:** `tests/cpp/unit/test_cdt_detria_backend.cpp`,
`tests/cpp/unit/test_cdt_constrained_edges.cpp` and
`tests/cpp/property/prop_cdt_invariants.cpp`. Mutation testing applies to those
three. `test_cdt_backend_seam.cpp` and `test_cdt_indexed_mesh.cpp` are not
invariant-critical; the first is a link-level test whose failure mode is a link
error and the second is a value type. Reasons under "The invariant-critical
suite".

## The finding that changes the design: detria closes its own polylines

The brief, increment 3's "seam to increment 4", and `pslg.hpp`'s comment on
`closed_index_buffer_size` all rest on one premise: that `addOutline` and
`addHole` need a **closed** index span, so the wrapper must materialise a
closing index somewhere, so it must own a scratch buffer, so that buffer's
lifetime is the hazard this increment owns.

The premise is false at the pinned SHA. `detria.hpp:3543-3560`,
`createConstrainedEdges`:

```cpp
// We allow polylines to have the same start and end vertex, e.g. [0, 1, 2, 3, 0]
if (polyline.front() == polyline.back()) { ... startingIndex = 1; }
else { prevVertexIdx = polyline.back(); startingIndex = 0; }
```

The `else` arm is the open case, and it seeds `prevVertexIdx` with the **last**
index so the first edge emitted is `(back, front)`. An open span is closed by
detria itself; an explicitly closed span is accepted as an alternative
spelling. Measured against the vendored header: the unit square as
`{0,1,2,3}` and as `{0,1,2,3,0}` both triangulate to the same 2 interior
triangles, and `{0,1,2}` — a triangle, open — gives 1.

So `Pslg::indices_of(c)` goes **straight into `addOutline`/`addHole`**. It is
already a contiguous span of `std::uint32_t` into a buffer the `Pslg` owns, and
the `Pslg` outlives the `Triangulation` by construction (below). Consequences,
all of them in the direction of less code:

- **No scratch buffer, no `ClosedChainSpans` type, no closing index anywhere in
  the project.** The rule increment 3 wanted ("the only place a closing index is
  materialised") becomes "there is no such place", which is the stronger form of
  the same rule and needs no enforcement mechanism.
- **The dangling-span hazard mostly evaporates.** It was a hazard about a buffer
  that no longer exists. What remains is the far simpler question of whether the
  `Pslg` outlives the `Triangulation`, asked once, at one call site.
- **`closed_index_buffer_size(const Pslg&)` has no caller and never will.** It is
  five lines of shipped production code whose header comment describes a caller
  this increment has now decided not to write.

**Ruling: delete `closed_index_buffer_size` and its property assertion in this
increment's PR.** `docs/increments/README.md` requires a defect found during an
increment to be fixed in that increment's PR or not recorded, and an exported
function whose documentation asserts a caller that does not exist is worse than
a missing one — the next reader writes the caller. The deletion touches
`include/terrain/core/pslg.hpp` and `tests/cpp/property/prop_pslg_invariants.cpp`;
per the protocol the test deletion lands as its own commit naming this finding,
never folded into the green commit.

The surrounding comment in `pslg.hpp` **stays, rewritten**. Its real content is
"there is deliberately no `closed_indices_of(c)`", and that rule survives intact
and is now free: nobody needs a closed span, so nobody needs an accessor that
returns one.

**What this costs: a dependency on library behaviour rather than on library
API.** Auto-closing is a property of the vendored source at
`8aa25f3e0ded…`, not a documented API contract. It is taken deliberately, and
paid for with a **characterisation test** in the backend suite — the same square
given open and closed must produce identical output — carrying a comment that
says: if a re-pin fails this test, the wrapper must materialise closing indices,
and the design for that is one scratch buffer sized by a restored
`closed_index_buffer_size`. That is a five-line change, and knowing where it
goes is the whole value of writing this down.

## Where `IndexedMesh2` lives

`include/terrain/core/indexed_mesh.hpp`, namespace `terrain`. **Not** in `cdt/`,
**not** in `mesh/`.

The CDT *returns* it and `mesh` later builds its ternary tree *from* it. If it
lives in `cdt/`, then `mesh` — and eventually `refinement`, `flip` and the
bindings — depend on `cdt`, which is the module they have the least to do with.
If it lives in `mesh/`, then `cdt` depends upward on the module that consumes
it. Sited in `core/`, both depend downward and neither depends on the other,
which is what lets a future backend and a future mesh be replaced independently.
This is the same argument increment 3 made for `Pslg` and it has the same shape:
a type owned by one of its two consumers acquires that consumer's vocabulary.

`project_structure.md`'s dependency diagram draws `cdt → mesh`. That arrow is
wrong and this PR redraws it; see "Documentation this PR fixes".

## The mesh

```cpp
namespace terrain {

using TriangleIndices = std::array<std::uint32_t, 3>;

// Bit e of a triangle's mask is set iff the edge (v[e], v[(e+1) % 3]) is a
// constraint edge of the input. Bit 0 is v0->v1, bit 1 is v1->v2, bit 2 is
// v2->v0. See "The constrained-edge mask" for what is NOT the convention.
inline constexpr std::uint8_t kEdgeMask01 = 1u << 0;
inline constexpr std::uint8_t kEdgeMask12 = 1u << 1;
inline constexpr std::uint8_t kEdgeMask20 = 1u << 2;

class IndexedMesh2 {
public:
    IndexedMesh2() = default;
    IndexedMesh2(std::vector<Point2> vertices,
                 std::vector<TriangleIndices> triangles,
                 std::vector<std::uint8_t> constrained_edges);

    [[nodiscard]] std::span<const Point2>          vertices() const noexcept;
    [[nodiscard]] std::span<const TriangleIndices> triangles() const noexcept;
    [[nodiscard]] std::span<const std::uint8_t>    constrained_edges() const noexcept;

    [[nodiscard]] std::size_t triangle_count() const noexcept;
    [[nodiscard]] bool        empty() const noexcept;   // triangle_count() == 0

    // Precondition: t < triangle_count(), e < 3.
    [[nodiscard]] Segment2 edge(std::size_t t, std::size_t e) const noexcept;
    [[nodiscard]] bool     is_constrained(std::size_t t, std::size_t e) const noexcept;

private:
    std::vector<Point2>          vertices_;
    std::vector<TriangleIndices> triangles_;
    std::vector<std::uint8_t>    constrained_edges_;
};

}  // namespace terrain
```

Structure-of-arrays, as settled: `triangles_` and `constrained_edges_` are
parallel and index-aligned, which is what lets a refinement pass stream the mask
without touching coordinates.

**The constructor is public, and that is a deliberate break from `Pslg`.**
`Pslg`'s constructor is private because it has exactly one legitimate producer
and its existence is a proof that validation ran. `IndexedMesh2` has plural
producers on purpose — `DetriaBackend` today, `FakeCdtBackend` in the seam test,
`mesh`'s flatten-to-final step later, and eventually an ingress that hands
Python an existing TIN. A private constructor with a growing friend list is a
worse design than an honest value type, and — decisively — **a fake backend that
cannot construct a mesh cannot prove the seam is real**, which is the one thing
this increment has to demonstrate.

What it guarantees:

1. `triangles().size() == constrained_edges().size()`.
2. Every index in `triangles()` is `< vertices().size()`.
3. Every mask is `< 8`.
4. `vertices()` **begins with** the producing PSLG's vertex buffer, element-wise
   and in order. Index `k` means the same point coming out as going in, which is
   what lets Python hold a parallel attribute array — the same guarantee, and
   the same reason, as `Pslg` guarantee 9. `DetriaBackend` appends nothing, so
   for it the two arrays are equal; the seam permits a future backend that
   inserts Steiner points to append them.
5. Immutable after construction; `const IndexedMesh2` is safe for concurrent
   read from any number of threads, which `refinement` requires.

1–3 are **debug asserts in the constructor, not release checks and not a status
channel**. The mesh is produced in bulk by code inside this repo, so a violation
is a programmer error, and the asan+ubsan Debug job runs every assert on every
fixture. The O(3T) scan is not worth paying in release on the one path where
throughput matters.

What it explicitly does **not** guarantee:

- **That every vertex is referenced by a triangle.** A `Pslg` may carry
  unreferenced vertices, and a vertex outside the outer ring or inside a hole is
  referenced only by triangles the wrapper discards. Verified: a square outline
  plus one point at `(9,9)` triangulates to 2 interior triangles over a
  5-element vertex array. The Euler check below counts *referenced* vertices for
  exactly this reason.
- **That the mesh is Delaunay, conforming, or non-empty.** Those are claims
  about a particular producer, asserted in that producer's suite.
- **Any adjacency.** No neighbour array, no half-edge, no vertex-to-triangle
  map. `mesh` owns topology; this is the flat handoff. Adding adjacency here
  would put the ternary tree's vocabulary in `core/`.

## The result channel: a status, a message, no exceptions

```cpp
namespace terrain::cdt {

enum class CdtStatus : int {
    Ok = 0,
    NotRun,               // default-constructed outcome; no backend was called
    NotNoded,             // constraints cross, or a vertex lies on a constraint
    DegenerateGeometry,   // coincident points, repeated indices, all-collinear input
    InvalidTopology,      // hole outside every outline, stacked or overlapping rings
    MalformedInput,       // a precondition the Pslg validator already excludes
    BackendFailure,       // the library failed in a way we do not classify
};

[[nodiscard]] constexpr std::string_view describe(CdtStatus) noexcept;

struct CdtOptions {
    bool delaunay{true};
};

struct CdtOutcome {
    CdtStatus     status{CdtStatus::NotRun};
    std::string   message;
    IndexedMesh2  mesh;

    [[nodiscard]] bool ok() const noexcept { return status == CdtStatus::Ok; }
};

}  // namespace terrain::cdt
```

**Status plus message, and no exception crosses the backend seam.** The reasons
are the ones the brief records and they are worth keeping in the header:
`std::expected` is C++23 and this project is C++20; a future backend may be a C
library, and a status is the only channel a C library can speak; and the
pybind11 layer must translate to a Python exception regardless, so an exception
here would be caught and rebuilt one frame later.

**`status`'s default member initializer is `NotRun`, a failure value.** Same
argument as `Chain::role` defaulting to `Breakline`: a default-constructed
`CdtOutcome` that reaches a reader must not be able to claim success. It also
gives the enumerator a real producer, so it is testable rather than decorative.

**This is a different channel from increment 3's diagnostics vector, and the
difference is the same one increment 3 named — keep it explicit in both
headers.** `PslgBuilder` validates *input authored outside the process*, where
several chains are commonly wrong at once for one underlying reason (a layer
digitised clockwise), so reporting one failure turns one fix into N round trips.
This is a *backend seam*: one call to a vendored library, which stops at its
first problem and has exactly one thing to say about it. Collecting a vector
here would be collecting a vector of length one. Malformed-by-construction input
is a precondition violation caught by the PSLG validator, not a status code —
which is precisely what `MalformedInput` exists to detect if it ever is not.

**`message` carries detria's own diagnosis.** `getErrorMessage()` produces
`"Point 4 is exactly on a constrained edge (between points 0 and 1)"` and
`"Multiple input points had the same position (1, 1), at index 5 and at index
4"`. Discarding that in favour of a bare enum throws away the most useful
information at the exact point of failure — the offending vertex indices, which
our enum has no field for.

**`getErrorMessage()` builds through `<sstream>`, so it is called only on the
failure path.** On success `message` is left **empty** — not detria's
`"The triangulation was successful"`. `message.empty()` on `Ok` is a pinned
invariant so the mutant that formats a message on the happy path is killed by
an assertion rather than by a profiler.

**`CdtOptions` has one member and that is not an oversight.** `delaunay` is
forwarded to `triangulate(bool)` and earns its place by being the knob the
property suite turns to prove the Delaunay property is *earned*: with
`delaunay=false` detria produces a valid constrained triangulation that fails
the visibility-Delaunay property, so the mutant "pass `false` unconditionally"
dies. The struct itself exists so that adding a second option later does not
change the `CdtBackend` signature and therefore does not touch every backend.
Nothing else goes in it — not a kernel, not a tolerance (this project has no
tolerances), not a "compute the mask" flag, since an optional output doubles
the test matrix to save an O(T) loop.

### Status mapping — exhaustive, by switch, never by cast

Every `detria::TriangulationError` maps to exactly one `CdtStatus`. The mapping
is grouped by **what the caller should do**, not by where in detria the failure
arose; the enumerator name and detria's message go into `message`, so nothing is
lost by the grouping.

| `detria::TriangulationError` | `CdtStatus` | Reachable from a valid `Pslg`? |
|---|---|---|
| `NoError` | `Ok` | yes |
| `TriangulationNotStarted` | `MalformedInput` | no — we always call `triangulate` |
| `LessThanThreePoints` | `MalformedInput` | no — an `Outer` chain forces `count >= 3` |
| `NonFinitePositionFound` | `MalformedInput` | no — PSLG stage 3, including unreferenced vertices |
| `PolylineIndexOutOfBounds` | `MalformedInput` | no — PSLG stage 2 |
| `PolylineTooShort` | `MalformedInput` | no — PSLG stage 1 |
| `DuplicatePointsFound` | `DegenerateGeometry` | **yes** — `Pslg` permits coincident vertices |
| `PolylineDuplicateConsecutivePoints` | `DegenerateGeometry` | **yes** — `Pslg` permits repeated indices |
| `AllPointsAreCollinear` | `DegenerateGeometry` | no — a collinear `Outer` ring is `DegenerateRing` |
| `PointOnConstrainedEdge` | `NotNoded` | **yes** — a T-junction; increment 5's job |
| `ConstrainedEdgeIntersection` | `NotNoded` | **yes** — crossing constraints; increment 5's job |
| `HoleNotInsideOutline` | `InvalidTopology` | **yes** — the nesting check increment 3 deferred |
| `StackedPolylines` | `InvalidTopology` | **yes** — same |
| `EdgeWithDifferentConstrainedTypes` | `InvalidTopology` | **yes** — a hole sharing an edge with its outline |
| `AssertionFailed` | `BackendFailure` | a bug in detria |

The `MalformedInput` rows are the interesting ones: **each is something the PSLG
validator already excludes, so `MalformedInput` firing means either the `Pslg`
did not come from the builder or our wrapper mis-built a span.** It is a
self-check with a name, not a user-facing status, and `describe()` says so.

Spelled as a `switch` naming every enumerator, **never a cast** — increment 1
established that rule for detria's orientation enums and mutation-tested the
failures a cast produces. Here the switch has a second job: a `return
BackendFailure` **after** the switch rather than in a `default:` arm, so
`-Wswitch` under `-Werror` turns a re-pin that adds an enumerator into a
compile error, while the runtime path stays total.

## The backend seam

```cpp
// include/terrain/cdt/triangulate.hpp
namespace terrain::cdt {

template <typename B>
concept CdtBackend = requires(const Pslg& pslg, const CdtOptions& options) {
    { B::triangulate(pslg, options) } -> std::same_as<CdtOutcome>;
};

template <CdtBackend B>
[[nodiscard]] CdtOutcome triangulate(const Pslg& pslg, const CdtOptions& options = {});

}  // namespace terrain::cdt
```

Spelled with a **qualified static call**, matching `GeometryKernel`: the concept
itself requires models to be usable without an instance, so statelessness is
part of the signature rather than a separate purity test that a future model
would not be covered by. A backend that needs scratch storage allocates it
inside the call; `parallel_refinement.md` has every refinement thread sharing one
constraint set, and a stateful backend object is where that stops being true.

**The generic entry point is not a forwarding wrapper.** It checks the seam's
postconditions, in debug, for every backend including ones not written yet:

- `out.ok() == !out.mesh.empty()` — no `Ok` with an empty mesh, no mesh
  alongside an error. The first half is not hypothetical: a `Triangulation` with
  points but **no** `addOutline` call triangulates successfully and yields
  **zero** interior triangles, so "forgot to add the outline" is a silent
  success without this assert.
- `out.ok() == out.message.empty()` — success says nothing, failure says
  something. Note the operator: an earlier draft wrote `!=`, which is false in
  *both* cases (on success both sides are true, on failure both are false), so
  the assert fired on every call. @tester transcribed it literally and three
  suites aborted before running.

Checking a contract at the seam once is the difference between a concept that
names a signature and a concept that means something. A backend's own suite can
still assert more.

**Semantic obligations of the concept**, stated in the header because a concept
cannot express them and a model that violates one is a wrong mesh, not a compile
error:

1. The returned vertex array begins with `pslg.vertices()`, element-wise, in
   order. Additional vertices may only be appended.
2. Only in-domain triangles are returned: nothing inside a `Hole`, nothing
   outside every `Outer`.
3. Triangles are counterclockwise under an exact orientation predicate.
4. The mask convention below is the backend's to honour.
5. The backend decides orientation and incircle at least as exactly as the
   kernel that validated the `Pslg`. Increment 3 pinned that a `Pslg` is valid
   *with respect to a kernel*; a backend using naive doubles could disagree with
   the validator about a sliver's winding. `DetriaBackend` satisfies this
   because `UseRobustOrientationTests` and `UseRobustIncircleTests` default to
   `true` and are the same Shewchuk arithmetic `DetriaExact` already exposes.
   **No kernel template parameter is threaded through the CDT.** detria does not
   accept one, threading a dead parameter to document an obligation is how a
   parameter later acquires a meaning nobody checked, and the obligation is
   stated where it can be read.

```cpp
// include/terrain/cdt/detria_backend.hpp -- declaration only, no detria include
namespace terrain::cdt {

struct DetriaBackend {
    [[nodiscard]] static CdtOutcome triangulate(const Pslg&, const CdtOptions&);
};

static_assert(CdtBackend<DetriaBackend>);

}  // namespace terrain::cdt
```

The `static_assert` sits in the header so the shipped backend is re-checked
against the concept in every TU that includes it — the cheapest possible
version of "the seam is not fake".

### Proving swappability instead of asserting it

`tests/cpp/unit/test_cdt_backend_seam.cpp` defines a `FakeCdtBackend` whose
`triangulate` ignores its arguments and returns a canned two-triangle mesh, and
runs `terrain::cdt::triangulate<FakeCdtBackend>` over a real `Pslg`.

What makes it a proof rather than a restatement is **the link line**. The target
is registered with the plain `add_terrain_test` helper: it links
`terrain_headers` and Catch2 and **nothing else** — not `terrain_cdt`, not
`terrain_predicates`, not `detria`. Therefore:

- If `triangulate.hpp` ever includes `detria_backend.hpp`, or the generic path
  ever names `DetriaBackend`, this target fails to **link**. A seam that
  survives only because both halves are always present is not a seam.
- The `Pslg` fixture is built with `std::move(b).build<pred::FastKernel>()`.
  `FastKernel` is header-only, so the suite needs no compiled target at all. On
  a well-separated unit square `FastKernel` and `DefaultKernel` cannot disagree,
  which is exactly the condition under which the cheaper kernel is legitimate —
  and it is also the fact that lets this file stay unlinked.
- The file carries the `#ifdef DETRIA_HPP_INCLUDED → #error` guard increment 1
  established, so a build change that puts `lib/detria` on a broad include path
  is caught here too.

The static properties are pinned alongside: `STATIC_REQUIRE(CdtBackend<FakeCdtBackend>)`,
and a `NotACdtBackend` with a wrong return type that must **fail** `CdtBackend`,
because a concept nothing fails is a concept that constrains nothing.

## The constrained-edge mask, and the convention that gets transposed

detria exposes no public per-edge query. `getTopology()` exists but is
documented "mostly used for testing" and returns an internal half-edge type, so
using it would drag detria's topology vocabulary across the seam for one bit per
edge. The mask is **ours to compute**, from the `Pslg` and the output triangles,
with no backend involvement — which is the good outcome, because it makes the
mask testable without a triangulation.

```cpp
// include/terrain/cdt/constrained_edges.hpp
namespace terrain::cdt {

// Every constraint edge of a Pslg, as an unordered index pair, sorted once.
class ConstraintEdgeSet {
public:
    explicit ConstraintEdgeSet(const Pslg&);

    [[nodiscard]] bool contains(std::uint32_t a, std::uint32_t b) const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;   // after dedup

    // (min(a,b) << 32) | max(a,b). Exposed so a test can pin the symmetry
    // rather than infer it: key(a,b) == key(b,a), for every pair.
    [[nodiscard]] static constexpr std::uint64_t key(std::uint32_t a, std::uint32_t b) noexcept;

private:
    std::vector<std::uint64_t> keys_;
};

[[nodiscard]] std::uint8_t constrained_mask(const ConstraintEdgeSet&,
                                            const TriangleIndices&) noexcept;

}  // namespace terrain::cdt
```

**A sorted `std::vector<std::uint64_t>` with `std::ranges::binary_search`, not a
hash set.** The key is exact and small, so hashing buys nothing; one allocation
of known size beats a rehashing container; the structure is read-only after
construction and therefore shareable by every refinement thread without
synchronisation, per `parallel_refinement.md`; and it is contiguous, which the
same document asks for by name. Sorting also dedups for free, which matters
because a breakline may legitimately repeat an edge.

**This is not "dedup" in the sense increment 3 prohibited.** That prohibition is
about *coordinate* keys, which need a snap grid that does not exist until
increment 5. These keys are vertex **indices**: exact integers with no tolerance
question. `constrained_edges.hpp` includes no `<unordered_map>` and no
`<unordered_set>`, and says why in one line.

**This works only because detria never splits a constraint edge.** It inserts no
Steiner points; a vertex lying exactly on a constraint is a hard
`PointOnConstrainedEdge` failure (`detria.hpp:3858`, `detria.hpp:3932`), not a
split. So every constraint edge appears in the output as **exactly one**
triangulation edge and a pair lookup is exact. That assumption belongs to
`DetriaBackend`, not to the seam — a backend that splits edges (poly2tri does
insert points) would need the mask computed by walking the split chain. Stated
in the header next to the class, because it is the assumption that silently
becomes false when the backend changes.

### The bit convention

**Bit `e` of a triangle's mask is set iff the edge from `v[e]` to
`v[(e + 1) % 3]` is a constraint edge.** Bit 0 is `v0→v1`, bit 1 is `v1→v2`,
bit 2 is `v2→v0`. `IndexedMesh2::edge(t, e)` returns that same segment, so the
mask and the accessor cannot drift.

**Explicitly rejected: the "edge `e` is opposite vertex `e`" convention**, which
CGAL uses and which a reader with that background will assume. It maps edge `e`
to `(v[e+1], v[e+2])` — our bit `e+1`. The two conventions are a **rotation of
each other**, so they agree on any triangle with zero or three constrained edges
and disagree on every other one. That is why the test is specific: a triangle
with **exactly one** constrained edge, in a known vertex order, asserting the
bit index by number. Anything weaker passes under both conventions.

Triangle winding is load-bearing here. `forEachTriangle`'s `cwTriangles`
parameter **defaults to `true`** (`detria.hpp:2917`), and the wrapper must pass
`false`. Under the default, `v1` and `v2` arrive swapped, which reflects the mask
— bits 0 and 2 exchange and bit 1 stays — so a forgotten argument is a silently
transposed mask *and* clockwise triangles. One line, two invariants, both
property-tested.

## Lifetimes: exactly what is handed to detria

`addOutline`, `addHole` and `setPoints` store a `detria::detail::ReadonlySpan`
and **do not copy** (`detria.hpp:1437-1452`, `2810`, `2826`, `2840`). The
caller's buffers must outlive `triangulate()`. With the scratch buffer gone, the
rule is short enough to state completely:

> The backend hands detria exactly two spans — `pslg.vertices()` and
> `pslg.indices_of(c)` — plus individual `std::uint32_t` values to
> `setConstrainedEdge`, which are copied into detria's own vector. **No buffer
> constructed inside the backend is ever passed to a detria call.**

Three things keep that true:

1. **The `Triangulation` is a local in the same function as the `const Pslg&`
   parameter** and is never returned, stored or moved out. It therefore cannot
   outlive the argument. It is held **by value on the stack**: `Triangulation`
   is move-only with allocator template parameters (`detria.hpp:2765-2767`;
   copy deleted, move-construct defaulted, no move-assign), so a
   `std::unique_ptr` buys nothing and a member would give the wrapper state the
   `CdtBackend` concept forbids.
2. **`Pslg` exposes no accessor that returns a closed index span**, so there is
   no `closed_indices_of(c)` to be tempted by and, after this increment, no
   reason to want one.
3. **`ConstraintEdgeSet` is constructed after `triangulate()` returns**, not
   before. It is the only buffer the backend owns, and building it on the far
   side of the call means it cannot be confused for something detria is holding.
   It also means it is not built at all on the failure path.

**This still cannot be a compile error**, and the honest statement is that it is
now a *smaller* thing that cannot be a compile error: a future edit that adds a
per-ring `std::vector<std::uint32_t>` and hands it to `addOutline` compiles,
dangles at `triangulate()`, and is caught only by the sanitizer job. So the
**multi-ring, multi-hole ASan test is still owed and is still the mechanism** —
with a single ring a dangling span usually still points at live memory and the
bug passes. The fixture is an outer ring, two disjoint holes and two breaklines,
run in the asan+ubsan Debug job, with a comment saying it passes by not
crashing.

## The wrapper, in order

`src/cdt/detria_backend.cpp` — the second and last TU in the project that
includes `detria.hpp`.

```
CdtOutcome DetriaBackend::triangulate(const Pslg& pslg, const CdtOptions& options):
    detria::Triangulation<terrain::Point2, std::uint32_t> tri;   // by value, on the stack
    tri.setPoints(pslg.vertices());                               // span, not copied

    for each chain c:
        Outer     -> tri.addOutline(pslg.indices_of(c))           // span, not copied, auto-closed
        Hole      -> tri.addHole(pslg.indices_of(c))
        Breakline -> for k in [0, pslg.edge_count(c)):
                         s = pslg.indices_of(c);
                         tri.setConstrainedEdge(s[k], s[k + 1])

    if (!tri.triangulate(options.delaunay)):
        return { map(tri.getError()), tri.getErrorMessage(), {} };   // sstream only here

    ConstraintEdgeSet constraints(pslg);
    vector<TriangleIndices> triangles;  vector<uint8_t> masks;
    triangles.reserve(tri.getMaxNumTriangles());                  // upper bound, one allocation
    tri.forEachTriangle([&](detria::Triangle<uint32_t> t) {
        TriangleIndices v{t.x, t.y, t.z};
        triangles.push_back(v);
        masks.push_back(constrained_mask(constraints, v));
    }, /*cwTriangles=*/false);                                    // NOT the default

    return { CdtStatus::Ok, {},
             IndexedMesh2{ {pslg.vertices().begin(), pslg.vertices().end()},
                           std::move(triangles), std::move(masks) } };
```

Notes that are rulings rather than commentary:

- **`terrain::Point2` is handed to detria directly.** detria's point access needs
  only `.x` and `.y` of a common scalar type; increment 1 established this for
  the predicates. No adapter, no coordinate copy, no `PointGetter`.
- **`forEachTriangle`, not `forEachTriangleOfEveryLocation`, not
  `forEachTriangleOfLocation`.** It filters to `TriangleLocation::Interior`
  (`detria.hpp:2917-2924`), so hole triangles and convex-hull triangles never
  reach the mesh. That is what makes the Euler check below meaningful; with the
  other variants it is unwritable.
- **Breaklines reach the backend edge by edge**, never as a polyline, and
  **`addPolylineAutoDetectType` is prohibited.** We know each chain's role from
  its `ChainRole`; letting the library re-derive it geometrically is a
  silent-wrong-answer channel. The concrete failure: a closed-loop breakline (a
  contour, a ring road) auto-detected as a hole would carve a void out of the
  domain. Verified in both directions — fed as individual constrained edges, a
  closed loop inside a square leaves the full 10-triangle interior intact.
- The loop indexes `s[k]` and `s[k+1]` over `k < edge_count(c)` directly rather
  than going through `Pslg::edge(c, k)`, because `edge` returns a `Segment2` of
  *points* and `setConstrainedEdge` wants *indices*. `edge_count(c)` is still
  the authority on how many edges an open chain has — which is the accessor's
  stated reason for existing.
- `getMaxNumTriangles()` counts interior, hole and convex-hull triangles, so it
  is an upper bound and the reserve never under-allocates.

## Degeneracy, verified against the vendored header

Each row below was measured against `lib/detria/detria.hpp` at the pinned SHA,
not inferred. Each is a named test in the backend suite, which is the mechanism a
sceptical reader re-runs.

| Input | detria | `CdtStatus` | Ruling |
|---|---|---|---|
| Two vertices with identical coordinates, **even unreferenced** | `DuplicatePointsFound` — the duplicate scan is over the whole point array (`detria.hpp:3178-3194`) | `DegenerateGeometry` | Not pre-checked by us. A coordinate key needs the snap grid; increment 5 owns it. This is the most likely failure on real data. |
| Repeated consecutive index in a chain (`{0,0,1,2,3}`) | `PolylineDuplicateConsecutivePoints` | `DegenerateGeometry` | `Pslg` accepts it (increment 2 §degeneracy 4), detria does not. Passed through as a status. |
| Hole touching its outer ring at **one shared index** | accepted; a square with a corner-touching triangular hole gives 5 interior triangles | `Ok` | Accepted, and the fixture is **excluded from the Euler check** — see below. |
| Hole touching its outer ring at **two coincident vertices** | `DuplicatePointsFound` | `DegenerateGeometry` | How you spell the touch decides whether it works. The header says: share the index. |
| Hole sharing a whole **edge** with its outer ring | `EdgeWithDifferentConstrainedTypes` | `InvalidTopology` | Rejected. A hole that shares an edge with its boundary is a notch in the boundary, not a hole. |
| A vertex lying exactly **on** a constraint edge it is not part of | `PointOnConstrainedEdge` | `NotNoded` | The T-junction `parallel_refinement.md` assigns to the noder. |
| **Collinear consecutive vertices within a ring** (`{0,0},{2,0},{4,0},{4,4},{0,4}`) | accepted, 3 interior triangles | `Ok` | Distinct from the row above and easily confused with it. Redundant ring vertices are fine; a *foreign* vertex on a constraint is not. |
| Crossing constraint edges | `ConstrainedEdgeIntersection` | `NotNoded` | The un-noded-input case; see the ruling below. |
| Hole outside every outline | `StackedPolylines` (not `HoleNotInsideOutline`, which fires on a different arm) | `InvalidTopology` | The nesting check increment 3 deferred, landing here as a status instead of a diagnostic. Both enumerators map to the same status, which is why the grouping is by action. |
| Outline inside a hole inside an outline (an island in a lake) | accepted, 9 interior triangles | `Ok` | Supported. detria computes nesting itself, which is why increment 3 stores no parent map. |
| Vertex outside the outer ring, or inside a hole | accepted; referenced by no interior triangle | `Ok` | Kept in the vertex array for index identity. |
| Empty PSLG | unreachable — stage 6 requires an `Outer` chain, which forces `count >= 3` | — | `LessThanThreePoints` maps to `MalformedInput` as a self-check. |
| Single-chain PSLG (one `Outer` triangle) | accepted, 1 interior triangle | `Ok` | The minimal valid input. A fixture. |
| Sliver outer ring (`{0,0},{4,0},{4,1e-13}`) | accepted, 1 triangle | `Ok` | No area threshold anywhere. "No near-zero-area triangle" is not an invariant this project can state — see `testing.md` fixes. |

**One `Pslg` guarantee we lean on and should name.** Stage 4 rejects a stored
closure, so `idx[begin] != idx[begin + count - 1]` as points. Combined with
auto-closing, that means the edge detria synthesises to close a ring can never be
a duplicate-consecutive pair. A concern discharged by an upstream check rather
than by a check here, which is what increment 3 exists for.

## The Euler check, stated so it can be written

For a triangulation of a polygon with holes:

```
triangles == 2n - b - 2 + 2h
```

with `n` the number of **referenced** vertices, `b` the number of vertices lying
on any boundary ring (outer or hole), and `h` the number of hole rings.
Unit square: `2·4 - 4 - 2 = 2`. Square with an interior point: `2·5 - 4 - 2 = 4`.
Square with a disjoint square hole: `2·8 - 8 - 2 + 2 = 8`. All three verified.

Two conditions on using it, both of which a naive test gets wrong:

1. **`n` counts referenced vertices, not `vertices().size()`.** A vertex outside
   the domain or inside a hole is in the array and in no triangle.
2. **The rings must be disjoint.** For the corner-touching hole the formula
   predicts 6 and the true answer is 5: the shared vertex is one vertex on two
   rings, and the domain is no longer a disk with a hole. The property helper
   therefore takes `n`, `b` and `h` from the fixture's declared structure, and
   the touching fixture is asserted with an explicit expected count and a comment
   naming this, rather than being fed to the general helper.

`testing.md` currently says "triangle count matches Euler's formula given vertex
count and boundary", which is not writable as a test. The replacement text is
below.

## Noded input: the ruling

**Increment 3's recommendation still holds. The signature reads `const Pslg&`
today and changes to `const NodedPslg&` mechanically at increment 5.**

Introducing a `NodedPslg` now — as an alias, a tag or a stub — would be a type
whose only constructor is a lie: there is no noder, so nothing could establish
its promise, and a name that asserts a property nobody checks is worse than the
honest weaker type. Increment 3 reserved the name; reserving it is the whole
mitigation.

**What the CDT does when handed a self-intersecting PSLG: it fails, loudly, with
`CdtStatus::NotNoded`, and returns no mesh.** It does not attempt repair, does
not insert intersection points, and — this is the load-bearing part — **never
silently produces a mesh that ignores a crossing**. detria detects the crossing
itself (`ConstrainedEdgeIntersection`) and so does the T-junction case
(`PointOnConstrainedEdge`); the wrapper's job is to translate, not to guess.
`describe(CdtStatus::NotNoded)` names the noder as the fix.

That makes the sequencing risk visible rather than latent: between increments 4
and 5 the wrapper is *correct* on all input and *useful* only on input that
happens not to cross. Its fixtures are hand-made non-crossing sets. When
increment 5 lands, `NotNoded` becomes unreachable by type and joins
`MalformedInput` as a self-check.

## Files and LOC

The unit is **non-comment production lines** — `CLAUDE.md` §2's unit, stated
explicitly because three increments of estimates were ambiguous about it. These
headers carry roughly as many comment lines as code lines, by house style, so
the raw `wc -l` will be about double.

| File | Contents | Est. LOC |
|---|---|---|
| `include/terrain/core/indexed_mesh.hpp` | `TriangleIndices`, the mask constants, `IndexedMesh2` | ~65 |
| `include/terrain/cdt/result.hpp` | `CdtStatus`, `describe`, `CdtOptions`, `CdtOutcome` | ~55 |
| `include/terrain/cdt/constrained_edges.hpp` | `ConstraintEdgeSet`, `constrained_mask` | ~55 |
| `include/terrain/cdt/triangulate.hpp` | the `CdtBackend` concept, the generic entry point | ~25 |
| `include/terrain/cdt/detria_backend.hpp` | `DetriaBackend` declaration, `static_assert` | ~12 |
| `src/cdt/detria_backend.cpp` | the wrapper and the status mapping | ~130 |

**~342 production LOC. No split needed.** Non-C++ changes on top, listed
separately because they are not production code but are part of the PR:

- `CMakeLists.txt`: a `terrain_cdt` STATIC target over the one TU, `terrain_headers`
  PUBLIC, `detria` **PRIVATE** — the privacy is what keeps `lib/detria` off the
  INTERFACE include path and makes the one-TU rule structural, exactly as
  `terrain_predicates` does.
- `tests/cpp/CMakeLists.txt`: an `add_terrain_cdt_test` helper (links
  `terrain_cdt` and `terrain_predicates`) plus five suite registrations. The seam
  suite deliberately uses the **plain** `add_terrain_test` helper; that is the
  test, not an oversight, and the CMake comment must say so or someone will
  "fix" it.
- `tools/check_detria_boundary.py`: `PERMITTED` gains
  `"src/cdt/detria_backend.cpp"`. **Extended, never bypassed** — the script's
  floor check (a permitted TU that stops including detria is also a finding)
  keeps working on both entries.
- `include/terrain/core/pslg.hpp`: delete `closed_index_buffer_size`, rewrite the
  comment above it. `prop_pslg_invariants.cpp` loses one assertion, in its own
  commit.

**Contingency split, dependency-ordered**, if implementation overruns — the
likely cause being the status mapping's fifteen `case` arms plus `describe`:

- **4a** — `core/indexed_mesh.hpp`, `cdt/result.hpp`, `cdt/constrained_edges.hpp`,
  `cdt/triangulate.hpp`, and the seam and mask suites. Everything with no detria
  dependency; buildable and testable against `terrain_headers` alone.
- **4b** — `cdt/detria_backend.hpp`, `src/cdt/detria_backend.cpp`, the CMake
  target, the boundary allowlist, and the backend and property suites.

That is the same seam increment 1 split on, for the same reason: 4a's code does
not depend on a library, and 4b's code is nothing but the dependency. I do not
expect to need it.

## The invariant-critical suite

**`tests/cpp/unit/test_cdt_detria_backend.cpp` — invariant-critical, mutation
round.** One named test per row of the degeneracy table above, asserting
`status` and — for the `Ok` rows — the exact interior triangle count. Plus:

- The **characterisation test** for auto-closing: the same square as an open and
  as an explicitly closed span produce identical meshes. Comment says what to do
  if a re-pin breaks it.
- `Ok` implies `message.empty()`; every failure implies a non-empty message.
  Never assert message *text* — the same rule increment 3 set, for the same
  reason, and here the text is a third party's.
- A default-constructed `CdtOutcome` is `NotRun` and not `ok()`.
- `describe` renders every enumerator, including `MalformedInput` and
  `BackendFailure`, which no reachable input produces.
- The polygon-with-holes happy path: outer ring, two disjoint holes, two
  breaklines. The goal-line fixture.
- The **multi-ring multi-hole ASan test**, run in the Debug sanitizer job,
  carrying a `READ BEFORE SIMPLIFYING` comment: it passes by not crashing.

**`tests/cpp/unit/test_cdt_constrained_edges.cpp` — invariant-critical,
mutation round.** No triangulation at all: a hand-built `Pslg`, a hand-built
`TriangleIndices`, and the bit convention asserted by index. `key(a,b) ==
key(b,a)`. `contains` on a ring's closing edge and on a breakline's last edge.
Cheap, and it is where a transposition dies.

**`tests/cpp/property/prop_cdt_invariants.cpp` — invariant-critical.**
Generators produce nested rectangle-with-holes families, ring resolutions and
breakline counts over a seeded range. Properties, each named:

- **Constraint preservation**: every ring edge and every breakline edge of the
  input appears as an edge of exactly one or two output triangles, and carries
  its mask bit on every triangle that has it.
- **Visibility-Delaunay**: for every interior non-constrained edge `(a,b)` with
  apexes `c` and `d`, either `incircle(a,b,c,d) != Inside` **or** the open
  segment `cd` crosses at least one constraint edge.
- **Euler**: `triangles == 2n - b - 2 + 2h`, with the two conditions above.
- **Orientation**: every output triangle is `CounterClockwise`.
- **In-domain**: no output triangle's centroid classifies `Inside` any hole ring,
  and every centroid classifies `Inside` some outer ring.
- **Index identity**: `mesh.vertices()` is element-wise equal to
  `pslg.vertices()` for this backend.
- **No status/mesh disagreement**: no input yields `Ok` with an empty mesh, and
  none yields a non-`Ok` status with a non-empty mesh.

**`tests/cpp/unit/test_cdt_backend_seam.cpp` — not invariant-critical.** Its
proof is the link line and the `STATIC_REQUIRE`s; a mutant of it is a build
failure, not a wrong answer. Mutation spend here buys nothing.

**`tests/cpp/unit/test_cdt_indexed_mesh.cpp` — not invariant-critical.**
Accessors, `edge`, `is_constrained`, copy and move. Its failure mode is a typo.

### Mutants the round must kill

1. `forEachTriangle` replaced by `forEachTriangleOfEveryLocation`, by
   `forEachHoleTriangle`, or by `forEachTriangleOfLocation(..., All)`.
2. `cwTriangles` left at its default `true`.
3. The mask read as "edge `e` is opposite vertex `e`", and the mask bits reversed
   `(2,1,0)`. Both die only on the one-constrained-edge fixture.
4. `triangulate(false)` instead of `triangulate(options.delaunay)`.
5. The status mapping replaced by a `static_cast`; any two of `NotNoded`,
   `DegenerateGeometry` and `InvalidTopology` collapsed into one.
6. `message` populated on the success path; `message` dropped on the failure
   path.
7. A `CdtOutcome` carrying a mesh alongside a non-`Ok` status, or `Ok` with an
   empty mesh — the "forgot `addOutline`" mutant, which is a *successful*
   detria call.
8. `addOutline` and `addHole` swapped for a role.
9. Breaklines fed via `addPolylineAutoDetectType`; killed by the closed-loop
   breakline fixture, which must not carve a hole.
10. The breakline loop running to `count - 1` instead of `edge_count(c)`, or a
    ring's closing edge omitted from `ConstraintEdgeSet`.
11. `ConstraintEdgeSet` built before `triangulate()` rather than after —
    behaviourally invisible, so it is **not** a mutant and is not pinned. Listed
    here so nobody spends the round trying.

### Template spend

**No cross product, and almost no templates.** The production path has exactly
one function template, `triangulate<B>`, instantiated twice in the whole suite —
once with `DetriaBackend`, once with `FakeCdtBackend` — in two different TUs.
The CDT takes no kernel parameter, so there is no kernel axis at all.

**The property suite runs under `DefaultKernel` only, and `FastKernel` is
forbidden there.** This is a ruling, not a default. The visibility-Delaunay
property evaluates `incircle` as an *oracle* against a mesh detria built with
its own robust predicates; an oracle less exact than the thing it judges reports
false failures on precisely the cocircular inputs the property is about. The
suite carries that sentence at the top, because `FastKernel` is header-only and
therefore the tempting choice for a suite that is slow to link.

## Corrections from the red suite

Four things the suite found, all measured against the pinned header rather than
argued. They are folded in above; recorded here because the reasoning is the
part that does not survive a diff.

1. **The entry point's second postcondition had `!=` where it needs `==`.** See
   above. It aborts on every call, success or failure.
2. **The auto-close characterisation test cannot be written as originally
   specified.** "The same square, given open and as an explicitly closed span,
   produces identical meshes" has no route: a `Pslg` cannot carry a closed ring
   (guarantee 5 — the validator rejects a stored closure), and no test target may
   include `detria.hpp`. What is written instead asserts the observable
   consequence on an L-shaped ring, chosen because its closing edge is a wall of
   the domain rather than a hull edge: status `Ok`, exactly 4 interior triangles
   (the closed hexagon's Euler value; every-location gives 5), and the closing
   edge present as a boundary edge used by exactly one triangle and flagged
   constrained. A detria that stopped auto-closing cannot produce that by
   accident.
3. **"Every constraint edge appears in the output" is false for legal input.** A
   breakline lying outside the outer ring, or inside a hole, is legal in a
   `Pslg` — which promises no disjointness and no nesting — triangulates `Ok`,
   and appears in **no** interior triangle. Measured. The property is therefore
   stated over in-domain constraints only, with the exception pinned by its own
   fixture. Unscoped it would be a coin toss on real data.
4. **Increment 3's generator cannot be reused.** `valid_chain_specs` scatters
   breaklines and holes over one disc independently, so generated breaklines
   routinely cross generated hole rings — `NotNoded`, no mesh. The CDT suite
   carries its own grid-cell generator with margins.

Two rows the degeneracy table was missing, both reachable from a valid `Pslg`:

- **`HoleNotInsideOutline` has two arms and the table had one.** A hole outside
  every outline gives `StackedPolylines`; `HoleNotInsideOutline` fires when a
  hole *contains* the outer ring, so the outermost polyline is a hole
  (`detria.hpp:4114`). Both map to `InvalidTopology`, so the mapping was right
  and only the table was short.
- **A closed-loop breakline must close on a shared index.** Spelled with a
  coincident closing vertex — which is what `pslg_cases.hpp`'s
  `closed_polyline()` produces — it is `DuplicatePointsFound`, not `Ok`. This is
  the corner-touching-hole lesson in a second place, and a fixture that gets it
  wrong silently stops testing what it was written for.

One correction to this document's own rationale for the multi-ring ASan test:
under ASan's quarantine a *single* ring also trips a per-ring temporary buffer,
and unsanitized the corrupted indices fail ordinary assertions too. The
multi-ring fixture is still right — a wrapper that reuses one buffer rather than
freeing per ring needs several rings to be caught — but "with a single ring the
bug passes" is not what was measured.

## Risks

1. **`DuplicatePointsFound` is a hard failure on input the `Pslg` accepts, and
   real data is full of coincident vertices** — two features digitised to the
   same corner, a lake boundary sharing a node with a river. Until the noder
   lands, every such input returns `DegenerateGeometry` and no mesh. This is the
   practical blocker on end-to-end use, it is larger than the crossing-constraint
   risk because it needs no pathology to trigger, and it is not fixable here: a
   coordinate-keyed dedup needs the snap grid. Named so increment 5 is scheduled
   with it in view.
2. **The auto-close dependency is on library behaviour, not library API.** A
   re-pin could change it. Mitigated by the characterisation test and by the
   five-line fallback being written down above.
3. **The mask assumes constraint edges are never split**, which is true of
   detria and false of at least one plausible replacement. It is a
   `DetriaBackend` assumption stated in `constrained_edges.hpp`; a splitting
   backend must compute the mask differently, and the seam does not promise
   otherwise.
4. **Un-noded input between increments 4 and 5**, inherited from increment 3 risk
   1 and discharged as far as it can be by `NotNoded` being a tested,
   named status rather than a wrong mesh.
5. **The `Ok`-with-empty-mesh mode is a *successful* backend call.** Only the
   entry point's debug assert and one named test stand between it and a silent
   empty pipeline; in a release build with the assert compiled out, the test is
   the only guard.
6. **`CdtOptions` will attract members.** Tolerances, a "repair" flag, a kernel.
   Each one multiplies the fixture matrix and the first two contradict rulings
   this project has already made. `delaunay` is there because a test turns it.

## Documentation this PR fixes

Per `docs/increments/README.md`, these are fixed in this PR or not recorded.
There is no ledger.

**`project_structure.md`**

1. The dependency diagram draws `cdt → mesh`. Wrong: `IndexedMesh2` is a `core`
   type, `cdt` produces one and `mesh` consumes one, so both arrows run down to
   `core` and neither module depends on the other. Redraw.
2. The directory listing under `include/terrain/` shows only `core/point.hpp`.
   It has been stale since increment 2. Add `bbox.hpp`, `segment.hpp`,
   `ring.hpp`, `pslg.hpp`, `pslg_builder.hpp`, `indexed_mesh.hpp` and the new
   `cdt/` headers, and un-mark `src/` as *(planned)* — `src/predicates/` exists
   and `src/cdt/` now does too.
3. The `cdt` section says "Detria preferred, poly2tri fallback" as though the two
   were interchangeable. poly2tri inserts Steiner points, which breaks the mask
   computation described here. State the real swappability claim — the
   `CdtBackend` concept, proved by `FakeCdtBackend` against an unlinked detria —
   and mark poly2tri as an unassessed option.

**`testing.md`**

4. The `cdt` section carries **no `[live]`/`[planned]` marker** while every
   sibling does. Mark it `[planned]` until this increment merges, then `[live]`.
5. Replace all four `cdt` bullets. The current text:
   - "Every constraint segment ... appears as a **connected chain** of
     triangulation edges" — understates and thereby hides a status: detria
     splits nothing, so each constraint edge is exactly **one** output edge, and
     anything that would need splitting is a `NotNoded` failure.
   - "The unconstrained sub-triangulation is Delaunay (incircle test passes for
     every interior edge)" — **false for every CDT with a constraint**, and
     increment 1 recorded the correction. Replace with the visibility form
     stated above.
   - "No degenerate triangle (zero or near-zero area below epsilon)" — there is
     no epsilon in this project (increment 2, "Tolerances: there are none") and
     a legitimate CDT over near-collinear terrain constraints produces slivers.
     Replace with: every output triangle is counterclockwise under an exact
     orientation predicate, which is the strongest true statement.
   - "Triangle count matches Euler's formula given vertex count and boundary" —
     not writable as a test. Replace with the exact formula, the definition of
     `n`, `b` and `h`, and the disjoint-rings condition.
   Add two bullets the current list lacks: the mask convention, and that no
   outcome carries a mesh with a non-`Ok` status or an empty mesh with `Ok`.
6. "What we do not test" says third-party libraries are trusted to their own
   suites. Add one sentence carving out **characterisation tests that pin
   behaviour we depend on at a pinned SHA** — the auto-close test is exactly
   that, and under the rule as written someone deletes it as out of scope.

**`parallel_refinement.md`**

7. "Library choices" says "The CDT library is used exactly once, at step 3".
   Step 3 is constraint noding; the CDT is step 4. One word. (The reference to
   "constraint edges from step 3" in the flip section is correct and stays.)

**`CLAUDE.md`**

8. §4 says `ctest` "runs test_point, test_raster, test_solar_position". That has
   been stale since increment 1 and this increment adds five more suites. Drop
   the enumeration rather than extending it every increment.
