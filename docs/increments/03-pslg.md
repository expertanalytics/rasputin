# Increment 3 — the PSLG

Status: design settled. Suite not yet written.

Ships `include/terrain/core/{pslg,pslg_builder}.hpp`. Header-only, namespace
`terrain`. Nothing goes in `src/`. Depends on increment 2 (`Ring`,
`IndexedRing`, `orientation<K>`, `Box2`, `Segment2`) and increment 1
(`pred::GeometryKernel`, `pred::Orientation`).

Makes possible: **a constraint set that has been checked once and never again.**
Every module above this one — the noder, the CDT wrapper, refinement — currently
has to either re-validate its input or assume it. This increment is the place
where the assumption is discharged, so that `Pslg` can be the type whose
existence is the proof.

**Invariant-critical suite:** `test_pslg_builder.cpp` and
`prop_pslg_invariants.cpp`. Mutation testing applies there, and the one
opt-in kernel instantiation described below. `test_pslg.cpp` covers accessors
and the span algebra; single instantiation, no mutation round.

## Why the PSLG type lives in `core/` and not in `noding/`

`project_structure.md` puts "PSLG construction" under `noding`, and that is
where the *noder* belongs. The PSLG *type* is not the noder's output alone: it
is also the noder's input, and the CDT wrapper's input. A type owned by one of
its two consumers acquires that consumer's vocabulary — a snap grid in the
validator, a river-override set in the accessors — and then the other consumer
depends upward.

So: `terrain::Pslg` is a core value type. The noder is a transformation
`Pslg -> NodedPslg`, and the dependency arrow runs `noding -> core` and
`cdt -> core`, never `core -> noding`.

## The data

```cpp
enum class ChainRole : std::uint8_t { Outer, Hole, Breakline };

[[nodiscard]] constexpr bool is_closed(ChainRole r) noexcept;   // Outer || Hole

struct Chain {
    std::uint32_t begin{};
    std::uint32_t count{};                       // distinct vertices; closure NOT counted
    ChainRole role{ChainRole::Breakline};
    bool is_river{false};
    friend constexpr bool operator==(const Chain&, const Chain&) = default;
};
```

`role`'s default member initializer is `Breakline` and that is deliberate.
`Breakline` is the role that promises the least — no winding contract, no
implied closure, no interior. A `Chain` that reaches a reader without having
been through the validator must not be able to claim it is a domain boundary,
and aggregate initialization of a missing enum member would otherwise give
`ChainRole(0)`. If the enumerator order is ever changed, this NSDMI is what
keeps the property; it is not decoration.

`count` counts distinct vertices. For `Outer` and `Hole` closure is implied and
the repeated last index is **not** stored, exactly as increment 2 specified for
`Ring`. For `Breakline` the chain is open and `count` is the vertex count of an
open polyline.

`is_river` is one bit per chain. It is permitted on every role — a wide river or
a lake is legitimately an area feature — and it is never validated, because it
is data, not structure. Only the noder (increment 5) can produce an edge whose
bit disagrees with its source chain, and it contributes a sparse override set
then, on `NodedPslg`, not here.

```cpp
class Pslg {
public:
    Pslg(const Pslg&) = default;   Pslg& operator=(const Pslg&) = default;
    Pslg(Pslg&&) = default;        Pslg& operator=(Pslg&&) = default;

    [[nodiscard]] std::span<const Point2>        vertices()      const noexcept;
    [[nodiscard]] std::span<const Chain>         chains()        const noexcept;
    [[nodiscard]] std::span<const std::uint32_t> chain_indices() const noexcept;

    // Sub-span of chain_indices(). Precondition: c < chains().size().
    [[nodiscard]] std::span<const std::uint32_t> indices_of(std::size_t c) const noexcept;

    // Zero-copy view of a closed chain. Precondition: is_closed(chains()[c].role).
    [[nodiscard]] IndexedRing ring(std::size_t c) const;

    [[nodiscard]] std::size_t edge_count(std::size_t c) const noexcept;  // count | count-1
    [[nodiscard]] Segment2    edge(std::size_t c, std::size_t k) const noexcept;

private:
    Pslg(std::vector<Point2>, std::vector<std::uint32_t>, std::vector<Chain>);
    friend class PslgBuilder;

    std::vector<Point2> vertices_;
    std::vector<std::uint32_t> chain_indices_;
    std::vector<Chain> chains_;
};
```

**`Pslg` owns all three buffers and is immutable after construction.** There is
no mutator, no non-const accessor, no lazy cache, no mutable member. A
`const Pslg&` is therefore freely shareable across threads with no
synchronisation, which is a requirement of `parallel_refinement.md` and not an
aspiration — every refinement thread reads the same constraint set.

**`Pslg` stores no span and no iterator into itself.** Every view is computed on
demand from the vectors. This is what makes the defaulted copy and move
constructors correct; a cached `IndexedRing` member would make a copied `Pslg`
point at the original's buffers, and that bug survives every test that never
copies. Do not add one.

The constructor is private with `PslgBuilder` a friend, so the only way to hold
a `Pslg` is to have passed validation. `build()` constructs the `Pslg` as a
local and moves it into the returned `std::optional` — `std::optional::emplace`
would construct from optional's context, which friendship does not reach.

`ring(c)` builds `IndexedRing{vertices_, indices_of(c)}`. Two spans, no
allocation, no gather; this is the whole reason `Ring` is a concept with a
non-owning indexed model rather than a type. It is **not** marked `noexcept`
even though it cannot throw: `IndexedRing`'s constructor rechecks size, closure
and the two boundary indices, and every one of those was established by the
validator, but the compiler cannot know that. A comment says so. That
constructor's checks being provably unreachable here is the cleanest evidence
that the validation order below is right.

`edge(c, k)` wraps for closed chains and does not for open ones, so
`edge(c, edge_count(c) - 1)` is the closing edge of a ring and the last segment
of a breakline. This is the accessor the noder's broad phase and the CDT
wrapper's `setConstrainedEdge` loop both consume, and it exists so neither of
them writes `(k + 1) % count` and gets the open case wrong.

```cpp
[[nodiscard]] std::size_t closed_index_buffer_size(const Pslg&) noexcept;  // sum(count + 1)
```

Discussed under "The seam to increment 4".

## The builder

```cpp
class PslgBuilder {
public:
    PslgBuilder() = default;
    explicit PslgBuilder(std::vector<Point2> vertices);          // takes ownership

    std::uint32_t append_vertices(std::span<const Point2>);      // returns first new index
    PslgBuilder& add_chain(std::span<const std::uint32_t>, ChainRole, bool is_river = false);
    PslgBuilder& add_chain(std::span<const Point2>,        ChainRole, bool is_river = false);

    template <pred::GeometryKernel K>
    [[nodiscard]] PslgBuildResult build() &&;
};
```

**A separate builder type, not a validating factory function.** A factory
`make_pslg(vertices, indices, chains)` requires the caller to have already built
the flat index buffer and computed every `begin`/`count` pair, which is the
arithmetic this increment exists to do once and correctly. The builder owns that
arithmetic: `add_chain` appends to the flat buffer and records
`begin = buffer.size() - indices.size()`, so `begin`/`count` are never written by
a caller and an off-by-one in them is unrepresentable rather than validated.

The index-taking `add_chain` is the primitive. The point-taking overload appends
its points to the vertex buffer **verbatim, with no dedup**, and emits the
consecutive index run; it exists so tests and simple Python ingress do not each
hand-roll the same loop. Callers whose geometry shares vertices between chains —
a hole touching its outer ring at a node, two breaklines meeting at a junction —
use `append_vertices` once and then the index overload, which is the encoding
this representation exists for.

`build()` is **rvalue-ref-qualified**: `std::move(b).build<DefaultKernel>()`.
The builder is consumed and its three buffers are moved into the `Pslg`, so a
successful build costs zero copies of the vertex data. It also makes "build
twice and get two `Pslg`s sharing nothing" a compile error rather than a subtle
question about what the second build sees.

`K` is the sole template parameter and is named at the call site. No default, for
the same reason `ring.hpp` gives none: defaulting it would make
`pslg_builder.hpp` include `default_kernel.hpp` and drag the compiled
`terrain_predicates` target into every consumer of a pure-header core type.

**`add_chain` never fails.** It is total: an empty span, a one-element span, an
out-of-range index are all recorded and rejected at `build()`. A builder that
can fail mid-accumulation needs a second failure channel, and two channels is
how half the failures end up unreported.

## The failure channel: diagnostics, not a status, not an exception

```cpp
enum class PslgError : int {
    NoOuterChain, ChainTooShort, IndexOutOfRange, NonFiniteVertex,
    StoredClosure, WrongWinding, DegenerateRing, VertexCountOverflow,
};

inline constexpr std::uint32_t kNoChain  = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::uint32_t kNoVertex = std::numeric_limits<std::uint32_t>::max();

struct PslgDiagnostic {
    PslgError error;
    std::uint32_t chain{kNoChain};    // index into chains(), or kNoChain
    std::uint32_t vertex{kNoVertex};  // index into vertices(), or kNoVertex
    std::string message;              // human-readable; NOT machine-parsed
};

struct PslgBuildResult {
    std::optional<Pslg> pslg;                      // engaged iff diagnostics.empty()
    std::vector<PslgDiagnostic> diagnostics;
    [[nodiscard]] bool ok() const noexcept { return pslg.has_value(); }
};

[[nodiscard]] std::string describe(const PslgDiagnostic&);
```

**This differs from the CDT-wrapper convention of a status enum plus one
message, and the difference is deliberate.** That convention exists to keep
exceptions from crossing a *backend seam* — the boundary to vendored third-party
code, where one failure is one call's failure. This is not a backend seam; it is
an *input validation* boundary, and its inputs are hand-authored or
Python-generated constraint sets in which several chains are wrong at once for
the same underlying reason (a whole layer digitised clockwise, a whole file with
a shifted index base). A channel that reports the first failure turns one fix
into N round trips through a pipeline whose cheapest stage is a DEM decode.

**So the validator is exhaustive: it never early-returns.** Every check runs on
every chain and every diagnostic is collected. "N independently broken chains
produce at least N diagnostics" is a named property test, and an early-return
mutant is one of the mutants the suite must kill.

**No exceptions of our own cross this boundary.** `std::vector` and
`std::string` may still throw `bad_alloc`; nothing else does.

This is a different channel from increment 2's, and the split is by error
*class*, not by taste. `PointRing`'s and `IndexedRing`'s constructors throw
`std::invalid_argument` because a malformed view is a **programmer** error on a
type built per query, where a returned status would be ignored at ten call
sites. A malformed constraint set is a **data** error arriving from outside the
process, where a status must be handled and a diagnostic must be legible.
`pslg_builder.hpp` carries that sentence as a comment, because the obvious
"simplification" is to make one of them match the other.

**Tests assert `error`, `chain` and `vertex`. Tests never assert `message`
text.** The message is `std::format`ted at diagnosis time and carries the
offending *values* — the out-of-range index and its position within the chain,
the observed winding — precisely because those do not fit the two index fields.
Pinning its wording in a test freezes prose and produces a suite that fails on
improvements to its own error messages.

## Validation: exactly what, in exactly what order

The order is not stylistic. Each stage establishes the precondition of the next,
and the last two stages evaluate geometry that is meaningless if an earlier one
failed. A chain that fails a stage is **excluded from the later stages** but does
not stop them running on other chains.

**0. Overflow.** `vertices.size()` and the flat index buffer's size must each fit
in `std::uint32_t`. `Chain::begin`/`count` are 32-bit; a 2^32-vertex DEM border
is not a use case, but a silent truncation is not an acceptable way to say so.
`VertexCountOverflow`, `chain = kNoChain`.

**1. Structure, per chain, O(1).** `count >= 3` for `Outer` and `Hole`;
`count >= 2` for `Breakline` — a one-vertex breakline is a point constraint, and
this pipeline has no such thing. `ChainTooShort`, naming the chain. (`begin` and
`count` need no bounds check against the flat buffer: the builder wrote them, and
the invariant that they partition the buffer is asserted by the property suite,
not re-derived at runtime.)

**2. Index range, O(total indices), once.** Every index in the flat buffer is
`< vertices.size()`. `IndexOutOfRange`, naming the chain, with the offending
value and its position in the message. **This is the check that makes
`IndexedRing::vertex(i)`'s unchecked read safe**, and increment 2 explicitly
deferred it here. It is done once over the whole buffer, never per query.

**3. Finiteness, O(vertices), once over the vertex buffer.** Not `all_finite(ring)`
per ring: a vertex shared by k chains would be scanned k times, and — decisively —
**unreferenced vertices must be checked too**, because the CDT wrapper hands
detria's `setPoints` the *entire* point array, so an unreferenced NaN reaches the
backend. `NonFiniteVertex`, `chain = kNoChain`, `vertex` = the buffer index.

Finiteness must precede every predicate call. `orient2d` on a NaN coordinate
returns `Collinear` by `orientation.hpp`'s documented total behaviour, so a
non-finite ring validated in the other order reports `DegenerateRing` — a
diagnostic that sends the reader to look for collinear geometry that is not
there.

**4. Stored closure, per closed chain, O(1).**
`vertices[idx[begin]] == vertices[idx[begin + count - 1]]` is rejected:
`StoredClosure`, naming the chain and the repeated vertex. Same rule and same
reason as increment 2 — closure is implied everywhere in this project, and
accepting a stored one gives every ring two encodings with every off-by-one bug
living in the gap. The comparison is on *points*, not indices: two distinct
indices onto coincident vertices close a ring just as surely as one index used
twice.

**Breaklines are exempt from this check.** A breakline is open, so coincident
first and last points are not an implied closure spelled twice — they are a
closed polyline, which is legitimate geometry (a contour, a ring road that is
not a domain boundary). Applying the check to `Breakline` is a mutant the suite
must kill in both directions.

**5. Winding, per closed chain, one `orientation<K>` call.** `Outer` must be
`CounterClockwise`, `Hole` must be `Clockwise`, `Collinear` is
`DegenerateRing` for both. A wrong winding is `WrongWinding` naming the chain,
its declared role and the observed orientation, and it is **never silently
reversed**: if the caller's geometry disagrees with the caller's declared role,
one of the two is a bug and guessing which hides it forever. Python normalises
with shapely, where reversal is cheap and testable; C++ asserts what it was
promised.

This stage runs on `Pslg::ring`-shaped `IndexedRing` views built from chains that
have already passed stages 1, 2 and 4 — which is exactly why `IndexedRing`'s
throwing constructor is unreachable inside the validator.

`signed_area` is used nowhere in this increment and must not be. Its sign
cancels to noise on a sliver at UTM33 magnitudes, and the winding of a hole is
precisely the topology decision increment 2 forbade it from driving.

**6. Domain, O(chains).** At least one chain has role `Outer`. `NoOuterChain`,
`chain = kNoChain`. A constraint set with no outer boundary has no bounded
domain to mesh, and `addOutline` is not optional for us.

## What `Pslg` guarantees, and what it refuses to

A downstream module may rely on all of the following **without re-checking**.
That list is the entire point of the type.

1. Every coordinate in `vertices()` is finite. Including unreferenced ones.
2. Every value in `chain_indices()` is `< vertices().size()`.
3. `chains()` is non-empty and at least one chain has role `Outer`.
4. Every closed chain has `count >= 3`; every breakline has `count >= 2`.
5. No closed chain stores its closure: its first and last vertices are distinct
   *points*.
6. Every `Outer` ring is counterclockwise and every `Hole` ring is clockwise,
   under the kernel that validated it. Neither is collinear.
7. The sub-spans `indices_of(0..n)` partition `chain_indices()` in chain order,
   contiguously, with no gap and no overlap: `begin_0 == 0` and
   `begin_{i+1} == begin_i + count_i`. A consumer may walk the flat buffer once.
8. `ring(c)` on a closed chain never throws.
9. The vertex buffer is element-wise equal to what the builder accumulated. **No
   dedup, no reordering, no reversal, no insertion.** A caller's index `k` means
   the same point going out as it did going in, which is what lets Python hold a
   parallel attribute array.
10. `const Pslg` is safe for concurrent read from any number of threads.

And it explicitly does **not** promise, and downstream must not assume:

- **That any chain is simple.** Self-intersection is not detected. Increment 2
  excluded segment intersection on purpose: intersection *construction* rounds,
  rounding needs the snap grid, and the snap grid is increment 5's vocabulary.
- **That two chains are disjoint.** No crossing, touching or overlap test of any
  kind runs here.
- **That a hole lies inside an outer ring**, or that breaklines lie inside the
  domain, or that the chains form a planar subdivision.
- **That vertices are distinct.** Duplicate points, coincident points, repeated
  consecutive indices and zero-length edges are all accepted.
- **That every vertex is referenced by a chain**, or that a vertex is referenced
  by at most one chain.

**Therefore a valid `Pslg` is not a valid CDT input.** It is a valid *noder*
input. The type asserts everything that can be decided without constructing a
point, and nothing that cannot. Increment 5 introduces `NodedPslg` — a distinct
type whose only producer is the noder, carrying the additional promise that no
two edges cross in their interiors and that coordinates lie on the snap grid.
Increment 4's wrapper signature should be written against `const Pslg&` today
and changed to `const NodedPslg&` then; that change is mechanical, and knowing it
is coming is cheaper than discovering it. It is listed under Risks.

## Grouping: neither declared nor computed, and not in this increment

The brief asks whether a polygon-with-holes is an explicit grouping or implied by
ordering, and whether a hole is associated with its outer ring by declaration or
by a `point_in_ring` nesting computation.

**The answer is neither, and no parent map is stored.** Roles are per-chain and
flat. There is no `Polygon` type, no group index on `Chain`, and no
`outer_of(c)` accessor. I am making this call against the obvious design, so
here is the reasoning in full:

- **detria computes nesting itself.** Increment 4 does not need a parent map, so
  it would be derived data with no consumer at the increment that produces it.
- **On un-noded input the computation is not merely expensive, it is wrong.**
  Nesting by representative vertex is well-defined only for non-crossing rings,
  and nothing here promises non-crossing rings. A lake polygon that pokes a
  metre over the DEM border — which the noder will clip without comment — would
  be rejected as "hole outside any outer ring". **A validator that rejects
  legitimate input is worse than one that misses an illegitimate one**, because
  the workaround for the former is to stop validating.
- **It would have to be recomputed after noding anyway**, since the noder splits
  and merges chains. Derived data with a short lifetime and two producers.
- It is the only check in the increment that costs O(C² · n) and the only one
  that would need the kernel for anything but winding.

The diagnostic value of catching a stray hole is real. It is paid for at the
right increment: a free function `nesting_forest<K>(const NodedPslg&)` returning
a parent forest with role alternation enforced (a `Hole`'s parent is an `Outer`;
an `Outer`'s parent is a `Hole` or nothing, which is how an island in a lake is
expressed), over input whose non-crossing property makes the answer meaningful.
Not in scope now, and named here so it is not re-derived.

Consequence to state plainly: **a hole inside no outer ring, and a hole inside
two, are both accepted by this increment without comment.** detria will resolve
them by its own nesting rules and the result may be a mesh the caller did not
intend. That is the cost of the deferral, taken knowingly.

## Duplicate and coincident vertices: the builder does not dedup

Three independent reasons, any one sufficient:

1. **There is no legal key.** `std::hash<Point2>` exists in `point.hpp` and
   `point.hpp` itself says it is for finding a *known* point in a hashed
   container, not for coordinate dedup. Hashing raw doubles makes two points a
   nanometre apart distinct and gives no control over what "the same point"
   means. Anything dedup-shaped in this project keys on **snapped integer
   coordinates**, and the snap grid does not exist until increment 5.
2. **Dedup is a mutation of caller-declared topology**, and this increment's
   entire posture is check, do not fix — the same posture that refuses to
   reverse a wrongly wound ring. Collapsing two coincident vertices silently
   changes which chains share a node, which changes the mesh.
3. **Exact duplicates are harmless to everything `Pslg` promises.** A repeated
   consecutive vertex gives a zero-length edge; `on_segment` reduces to point
   equality, `point_in_ring`'s half-open parity test is unaffected, and
   `orientation`'s collinear walk steps over it. Increment 2 specified all three
   and the noder is what makes them go away.

`pslg_builder.hpp` includes no `<unordered_map>`, `<unordered_set>` or `<map>`,
and specifies no key. That is checkable by grep and worth a comment.

**The rotation hazard, inherited from increment 2.** `{A, A, B, C}` is a legal
ring and its rotation `{A, B, C, A}` is an illegal stored closure. Increment 3 is
safe from this only because **it never rotates a chain**: the stored index order
is the order `add_chain` received, byte for byte, and guarantee 9 above says so.
Any future normaliser — one that rotates a chain to start at its lowest vertex,
say, to make two chains comparable — **must collapse an adjacent duplicate
first**. That sentence goes in the header above `add_chain`, not only here.

## The seam to increment 4

What the CDT wrapper reads, and nothing else:

| Wrapper needs | Reads | detria call |
|---|---|---|
| the point array | `vertices()` | `setPoints` |
| outer boundaries | `ring(c)` / `indices_of(c)` for `role == Outer` | `addOutline` |
| holes | `ring(c)` / `indices_of(c)` for `role == Hole` | `addHole` |
| breaklines | `edge(c, k)` for `k < edge_count(c)` | `setConstrainedEdge(a, b)` |
| provenance for diagnostics | `chains()[c]` | — |

Breaklines reach detria **edge by edge**, never as a chain. The chain structure
exists for our benefit — provenance, the `is_river` bit, validation — not the
backend's.

**The one-scratch-buffer rule.** `addOutline`, `addHole` and `setPoints` store a
`ReadonlySpan` and do not copy; the caller's buffers must outlive
`triangulate()`. `vertices()` and `indices_of(c)` are views into buffers the
`Pslg` owns, so they satisfy that for free — as long as the `Pslg` outlives the
triangulation, which is one lifetime question asked once instead of once per
ring. But detria wants a *closed* contiguous span, and `indices_of(c)` is not
closed. The wrapper therefore builds **one** scratch index buffer per
triangulation — each closed chain followed by its repeated first index — and
hands detria sub-spans of it. That is the only place in the project where a
closing index is materialised.

Two things make that hard to violate by accident, and one thing does not:

- **`Pslg` exposes no accessor that returns a closed index span.** There is
  deliberately no `closed_indices_of(c)`. The only way to get one is to allocate,
  and the only correct place to allocate is once, up front.
- **`closed_index_buffer_size(const Pslg&)` returns the exact element count**,
  `Σ (count + 1)` over closed chains. It exists so the wrapper's shape is
  "reserve exactly this, then fill", under which a per-ring `std::vector` inside
  the loop is visibly wrong rather than merely wrong. Five lines, and it is the
  only concession `core/` makes to a detria-shaped need.
- What does *not* enforce it: the rule cannot be a compile error the way
  increment 2's deleted rvalue constructors are. A per-ring local
  `std::vector<std::uint32_t>` handed to `addOutline` compiles, dangles at
  `triangulate()`, and is caught only by the asan job. Increment 4 owes a test
  that exercises a multi-ring PSLG under asan for exactly this reason.

The closed-span type itself (`ClosedChainSpans`, or whatever increment 4 calls
it) belongs to `terrain::cdt`, **not** to `core/`. Putting a type whose purpose
is to store a closure in the header whose purpose is to say closure is implied
guarantees the next reader uses it for something else.

## Files and LOC

| File | Contents | Est. LOC |
|---|---|---|
| `include/terrain/core/pslg.hpp` | `ChainRole`, `is_closed`, `Chain`, `Pslg`, `closed_index_buffer_size` | ~190 |
| `include/terrain/core/pslg_builder.hpp` | `PslgError`, `PslgDiagnostic`, `PslgBuildResult`, `describe`, `PslgBuilder`, `detail::validate` | ~250 |

**~440 production LOC. No split.** Under the 700-LOC gate with room, and there
is no dependency-ordered seam worth cutting: every check in the validator is
cheap and they only make sense as one ordered pass.

The two-header split is not arbitrary. `cdt` and `refinement` consume `Pslg` and
have no business instantiating the validator's templates; `pslg.hpp` therefore
includes `ring.hpp` (for `IndexedRing`) and nothing from `predicates/` beyond
what `ring.hpp` already pulls, while `pslg_builder.hpp` is where `K` appears.

If implementation overruns — the likely cause being `describe` growing one
`std::format` call per error enumerator — the seam is
**3a = `pslg.hpp` plus the kernel-free validation (stages 0–4)**,
**3b = winding (stage 5) and the domain check**, in that dependency order,
mirroring how increment 1 split at the point where a dependency arrived. I do
not expect to need it.

Nothing goes in `src/`. No existing header changes. `CMakeLists.txt` at
`tests/cpp/` gains three `add_terrain_backend_test` entries (all three name a
kernel, so none uses the plain helper).

## The invariant-critical suite

**`tests/cpp/unit/test_pslg_builder.cpp` — invariant-critical, mutation round.**
One named test per `PslgError` enumerator, each asserting the error code, the
chain index, and the vertex index where the error sets one. Plus:

- A build with three independently broken chains yields at least three
  diagnostics — the exhaustiveness property.
- A `Breakline` with coincident first and last points **builds successfully**,
  and the same index sequence declared `Outer` is rejected `StoredClosure`.
- A non-finite *unreferenced* vertex is rejected.
- A ring whose winding is correct but which contains a NaN is reported
  `NonFiniteVertex`, not `DegenerateRing` — the stage-ordering test.
- `ok()` and `pslg` agree in both directions; a failed build yields no `Pslg`.

Mutants this suite must kill, named so the round is not improvised:

1. `Outer`/`Hole` winding expectations swapped.
2. `Collinear` accepted rather than rejected for either role.
3. The closure check extended to `Breakline`; the closure check skipped for
   `Hole`; the closure check comparing indices instead of points.
4. The index range check using `<=` instead of `<`.
5. Finiteness restricted to referenced vertices.
6. `count >= 3` weakened to `>= 2` for closed chains, or `>= 2` to `>= 1` for
   breaklines.
7. Stage 5 moved before stage 3 (winding before finiteness).
8. The validator early-returning on the first diagnostic.
9. `NoOuterChain` dropped, or satisfied by a `Hole`.

**`tests/cpp/property/prop_pslg_invariants.cpp` — invariant-critical.**
Generators produce random chain sets over random vertex buffers, both valid and
deliberately corrupted. Properties: guarantees 1–9 above, each named; that
`indices_of` sub-spans partition the flat buffer exactly; that `ring(c)` never
throws for any `c` with a closed role; that
`closed_index_buffer_size(p) == Σ (count + 1)`; that a copied `Pslg` compares
equal element-wise and its spans point into *its own* buffers, not the
original's; that `edge(c, edge_count(c) - 1)` closes a ring and does not close a
breakline.

**`tests/cpp/unit/test_pslg.cpp` — not invariant-critical.** Accessors, span
algebra, `edge`/`edge_count`, copy and move. Single instantiation, no mutation
round; its failure mode is a typo.

**Template instantiations — one opt-in case, not a cross product.**
`Pslg` is kernel-free data and the only kernel-parameterised code is stage 5, so
the builder suite runs under `DefaultKernel` throughout. Exactly one
`TEMPLATE_TEST_CASE` over `{FastKernel, DefaultKernel}` exists, and it does one
job: a sliver `Hole` at UTM33 magnitudes whose winding `FastKernel` gets wrong
and `DefaultKernel` gets right, so the two kernels produce different build
outcomes from the same input. What it proves is that the validator's winding
decision actually flows through `K` — that it has not been accidentally written
kernel-independently, for instance by reaching for `signed_area` — and it pins
the fact increment 4 depends on: **a `Pslg` is valid with respect to a kernel,
and the CDT must be built with the same one.**

No `Ring`-model cross product. The validator builds only `IndexedRing`;
`PointRing` never appears in this increment, and instantiating over it would
test a code path production does not have.

## Risks

1. **Increment 4 consumes an un-noded `Pslg`.** The CDT wrapper lands before the
   noder, so between increments 4 and 5 the wrapper is safe only on inputs that
   happen not to cross. Mitigation: increment 4's fixtures are hand-made
   non-crossing sets, its header says in the first paragraph that it requires a
   noded PSLG, and the `NodedPslg` signature change is booked now rather than
   discovered later. This is the largest risk in the increment and it is a
   sequencing risk, not a design one.
2. **The deferred nesting check.** A stray hole reaches detria and produces a
   mesh nobody intended, with no diagnostic. Accepted; see the argument above.
   If it bites in practice before increment 5, the cheap partial mitigation is a
   *warning-level* diagnostic (a `PslgBuildResult` that is `ok()` and still
   carries entries) rather than a rejection — but that means introducing a
   severity field, and a severity field is how a validator starts having
   opinions. Do not add it speculatively.
3. **`Chain` is an aggregate and will be constructed by hand in a test**, with
   `begin`/`count` that do not match the flat buffer, producing a `Pslg`-shaped
   thing the validator never saw. Mitigated by the private constructor: there is
   no public path from a hand-built `std::vector<Chain>` to a `Pslg`. Tests that
   want a corrupt `Pslg` must go through the builder, which is the point.
4. **`build() &&` will be called on an lvalue** and produce a compile error that
   reads badly. Accepted; the fix is one `std::move` and the alternative is a
   silent copy of the whole vertex buffer.
5. **`describe` is the natural place for scope creep** — coordinates, ring
   previews, suggested fixes. It formats the error, the chain, the vertex and
   the offending value. Nothing else. Unbounded rendering of a ring was already
   excluded in increment 2 for the same reason.

## Documentation debt

**All three `project_structure.md` errors listed in increment 2 are still
present and still uncleared** (predicates under `src/geometry_predicates/`; the
"Shewchuk … header-only" prose, false in all three clauses; the directory
listing missing `include/terrain/predicates/` and `core/{segment,bbox,ring}.hpp`
and marking `lib/` as *(planned)*). They are not restated here — fix them there.

Increment 3 adds to that list:

- The directory listing needs `core/{pslg,pslg_builder}.hpp`, and the dependency
  diagram still needs the `core` node below `predicates` that increment 2 asked
  for. `core` is now the node both `noding` and `cdt` point at.
- **The `noding` section is wrong about the `is_river` bit.** It says the module
  "outputs a clean PSLG with one bit per edge (`is_river`)". The bit is **one per
  chain**, carried on `Chain`, and the noder contributes a *sparse override set*
  for the few edges whose bit disagrees with their source chain. Per-edge storage
  for a bit that is uniform over a chain is a per-edge allocation to represent
  nothing.
- The `noding` section should say the module *transforms* a `Pslg` into a
  `NodedPslg`, rather than reading as though it constructs the representation
  from nothing. `cdt`'s description — "translates between rasputin's PSLG
  representation and the library's API" — is now concrete and should name the
  type and the one-scratch-buffer rule.

**`testing.md`'s invariant catalog still starts at `noding` and still has no
entry for anything that exists.** It owes the retroactive `predicates` and
`core geometry` sections increment 2 identified, and it should gain a
`core geometry — PSLG` section carrying these, which are the guarantees above
stated as testable properties:

- A `Pslg` exists only if its build produced an empty diagnostics list; there is
  no partially valid `Pslg`.
- Every index in `chain_indices()` is `< vertices().size()`; every coordinate in
  `vertices()` is finite, including unreferenced ones.
- Every `Outer` ring is counterclockwise and every `Hole` ring is clockwise under
  the validating kernel; neither is collinear.
- No closed chain stores its closure; breaklines are exempt and may be closed
  polylines.
- `indices_of(i)` sub-spans partition `chain_indices()` contiguously in chain
  order, with no gap and no overlap.
- The vertex buffer is element-wise equal to the builder's input: no dedup, no
  reordering, no reversal.
- `ring(c)` never throws for a closed chain.
- Validation is exhaustive: N independently broken chains produce at least N
  diagnostics.
- A `const Pslg` is safe for concurrent read; no accessor mutates or caches.
- Negative, and equally load-bearing: a valid `Pslg` promises **no** simplicity,
  **no** pairwise disjointness and **no** nesting. Any test asserting one of
  those is testing the noder and belongs in its catalog.
