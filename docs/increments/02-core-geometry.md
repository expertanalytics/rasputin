# Increment 2 — 2D geometry value types

Status: design settled. Suite red. Implementation pending.

Ships `include/terrain/core/{bbox,segment,ring}.hpp`. Header-only, namespace
`terrain` directly — `terrain::raster` and `terrain::pred` are submodule
namespaces, `core` is not one. Nothing goes in `src/`.

Makes possible: an input domain exists, with a fixed winding contract and an
exact point-in-ring classification, before any triangulation code does.

**Invariant-critical suite:** `test_ring.cpp` and `prop_ring_invariants.cpp`.
Mutation testing and the kernel × model cross product apply there. `test_bbox`
and `test_segment` are value types; single instantiation, no mutation round.

## Two deviations from the original sketch

**No `Polygon` type and no `Breakline` type.** A polygon-with-holes is exactly
one `Outer` chain plus N `Hole` chains, and a breakline is one `Breakline`
chain — that is increment 3's `Chain`/`ChainRole` encoding. A second
representation now guarantees a conversion layer and a second ownership answer
later. Roles, grouping and `point_in_polygon` land with the PSLG.

**`Ring` is a concept with two models, not a type.** Increment 3's chain is
`vertices[chain_indices[begin..begin+count)]`. Those points are *not
contiguous*, so `std::span<const Point2>` cannot represent a chain without a
gather.

## bbox.hpp — `class Box2` (~120 LOC)

```cpp
constexpr Box2() noexcept;                    // EMPTY: lo = +inf, hi = -inf
Box2(const Point2& lo, const Point2& hi);     // throws: non-finite, or lo > hi
static constexpr Box2 empty() noexcept;

constexpr bool is_empty() const noexcept;
constexpr const Point2& lo() const noexcept;
constexpr const Point2& hi() const noexcept;
constexpr double width()  const noexcept;     // 0 when empty
constexpr double height() const noexcept;
constexpr Point2 center() const noexcept;     // precondition: !is_empty(); lo/2 + hi/2

constexpr void expand(const Point2& p) noexcept;   // precondition: p finite
constexpr void expand(const Box2& b)   noexcept;

constexpr bool contains(const Point2&) const noexcept;  // closed, exact
constexpr bool contains(const Box2&)   const noexcept;  // empty is contained
constexpr bool intersects(const Box2&) const noexcept;  // empty never intersects

friend constexpr bool operator==(const Box2&, const Box2&) = default;

Box2 bounding_box(std::span<const Point2>);   // validates; throws on non-finite
```

Invariant: either `is_empty()`, or `lo.x <= hi.x && lo.y <= hi.y` with all four
coordinates finite.

`center` is computed as `lo/2 + hi/2`, **not** `(lo + hi) / 2`: a box with
finite, legal corners can overflow the sum, and a box symmetric about the origin
hides it because `lo + hi` cancels to exactly zero.

Default construction is the **empty** box, not a degenerate box at the origin.
This is why `Box2` is a class and not an aggregate: an aggregate
`{Point2 min; Point2 max;}` value-initialises to a box at `(0,0)`, so
`Box2 b; for (p : pts) b.expand(p);` silently contains the origin. As
specified, that fold is correct and an empty input yields an empty box.

`expand(Point2)` carries a finiteness *precondition* and is `noexcept` — it is
the per-vertex inner loop and must not validate. **`expand(Box2)`'s precondition
is "finite *or empty*"**: the empty box's corners are ±inf, and expanding by an
empty sub-box must be a no-op, which is exactly what makes a fold over partial
boxes correct (a per-chunk or parallel bounding-box computation). An earlier
draft said "finite" for both overloads, which contradicted the empty sentinel. `bounding_box` is an entry point
and does validate. A NaN reaching `expand` corrupts the box silently; accepted,
and the reason finiteness is nailed down below.

`Box2`'s rendering is otherwise unspecified, with one requirement: the empty
box must render distinctly and **without `inf`**. Printing
`Box2(Point2(inf, inf), Point2(-inf, -inf))` reads as a bug report rather than
as the identity element.

Excluded: any tolerance, `intersection()` (constructive, no caller, needs an
empty-result convention), `expand(margin)`, `operator|` / `operator&`, any
`RasterGeometry` interop. A `Box2 box_of(const RasterGeometry&)` adapter belongs
in `raster/`, never `core/` — the dependency arrow runs raster → core.

## segment.hpp — `struct Segment2` (~100 LOC)

```cpp
struct Segment2 { Point2 a{}, b{};
                  friend constexpr bool operator==(const Segment2&, const Segment2&) = default; };

constexpr Segment2 reversed(const Segment2&) noexcept;
constexpr bool is_degenerate(const Segment2&) noexcept;   // a == b, exact
Box2 bbox(const Segment2&);                               // throws on non-finite

template <pred::GeometryKernel K>
bool on_segment(const Segment2&, const Point2&);
```

`operator==` is **ordered**: `{a,b} != {b,a}`. Ring edges are directed and the
noder's `is_river` merge cares. An unordered comparison, if ever needed, gets
its own name.

`on_segment` is the one predicate: point on closed segment, exact — collinearity
from `K::orient2d`, betweenness from closed comparison of input coordinates. No
division, no construction. On a degenerate segment it reduces to `p == s.a`,
which is correct and needs no special case.

## ring.hpp — the `Ring` concept and its algorithms (~230 LOC)

```cpp
template <typename R>
concept Ring = requires(const R& r, std::size_t i) {
    { r.size()    } -> std::same_as<std::size_t>;    // distinct vertices; closure NOT counted
    { r.vertex(i) } -> std::same_as<const Point2&>;  // precondition: i < size()
};
```

Two models, both **non-owning views**, both with rvalue-`std::vector`
constructors `= delete`d:

- `PointRing(std::span<const Point2>)` — contiguous. Tests, Python ingress, a
  draped DEM border.
- `IndexedRing(std::span<const Point2> vertices, std::span<const std::uint32_t> chain)`
  — `vertex(i)` is `vertices[chain[i]]`. The chain slice arrives already
  sub-spanned. What increment 3 hands down, zero-copy.

A concept rather than one type with an "empty index span means identity" mode:
two encodings behind one branch is the implicit-mode design this project
rejects elsewhere (detria's `addPolylineAutoDetectType` is the same shape of
mistake). Compile-time polymorphism is the house style — `RasterSource`,
`ExactPredicates`, `GeometryKernel` are all concepts.

Non-owning because increment 3 owns the vertex buffer and one flat index
buffer. An owning `Ring` forces a gather-and-copy per ring, which is both an
allocation per ring and the exact setup for the detria use-after-free:
`addOutline` stores a `ReadonlySpan` and does not copy, so a temporary ring
buffer built inside a loop dangles at `triangulate()`. Non-owning means the
buffer handed to detria is the same buffer the PSLG owns, and the lifetime
question is asked once at PSLG level rather than once per ring. The deleted
rvalue constructors turn the most likely form of that mistake into a compile
error.

For increment 4: detria's `addOutline` wants a *closed*, contiguous span, which
`IndexedRing` is not. The wrapper builds **one** scratch index buffer for the
whole PSLG — each chain followed by its repeated first index — and hands detria
sub-spans of that single long-lived buffer. One allocation per triangulation,
not per ring, and the only place a repeated closing index is ever materialised.

```cpp
template <Ring R> Segment2 edge(const R&, std::size_t i);   // {vertex(i), vertex((i+1) % size())}
template <Ring R> Box2 bounding_box(const R&);
template <Ring R> bool all_finite(const R&) noexcept;       // O(n), for increment 3's validator
template <Ring R> double signed_area(const R&) noexcept;    // magnitude only, kernel-free

template <pred::GeometryKernel K, Ring R> pred::Orientation orientation(const R&);

enum class PointInRing : int { Outside = -1, Boundary = 0, Inside = 1 };
template <pred::GeometryKernel K, Ring R> PointInRing point_in_ring(const R&, const Point2&);
```

The free `bounding_box(std::span<const Point2>)` and the
`template <Ring R> bounding_box(const R&)` coexist unambiguously **only because
no ring model is convertible to `std::span<const Point2>`**. `PointRing`'s span
constructor must stay `explicit` and must never gain an implicit conversion
back. Not testable except as a compile-failure test; it needs a header comment.

`edge` is a free function rather than a concept requirement so the concept stays
two lines — `size()` and `vertex(i)` are the whole obligation for a future third
model.

`K` is the **first** template parameter so call sites read
`point_in_ring<DefaultKernel>(ring, p)` with `R` deduced. No default argument
for `K`: defaulting it to `DefaultKernel` would make `ring.hpp` include
`default_kernel.hpp` and drag the compiled `terrain_predicates` target into
every consumer of a pure-header core type. Callers name their kernel. Only
`orient2d` is ever called; no narrower `Orient2dKernel` concept is introduced,
since a second name to keep in sync buys nothing.

Excluded: any owning ring, any `Polygon`, any simplicity or self-intersection
test, `is_convex`, centroid, `reversed(ring)` (a view cannot reverse without
owning), `std::formatter` for the models (unbounded length; diagnostics use
per-vertex `Point2` formatting), hashing of any kind.

## The winding contract

**Outer rings wind counterclockwise (positive area). Holes wind clockwise. The
interior is on the left of every directed edge, for both roles, uniformly.**

This polarity agrees with `Orientation::CounterClockwise = +1`, so
`orientation(outer) == CounterClockwise` needs no mental negation; it agrees
with OGC simple features and RFC 7946; and "interior on the left" is one
sentence covering both roles, which is what a hole flood-fill and the noder's
directed-edge merge consume. Holes-also-CCW makes the edge-direction statement
role-dependent and buys nothing.

**Checked, not enforced — and the check is not in increment 2.** Rings are
non-owning views: they cannot normalise a buffer they do not own, and reversing
requires a copy, the precise allocation this representation avoids. So
increment 2 ships the instrument, `orientation<K>(ring)`, and increment 3 runs
the check once per chain where the role is assigned: `Outer` requires
`CounterClockwise`, `Hole` requires `Clockwise`, `Collinear` is a hard error for
both.

**A ring wound the wrong way is rejected with a diagnostic naming the chain
index and role. It is never silently reversed.** If the caller's geometry
disagrees with the caller's declared role, one of the two is a bug, and guessing
which hides it forever. Python (shapely, where reversal is cheap and testable)
normalises before anything crosses the boundary; C++ asserts what it was
promised.

`point_in_ring` does **not** depend on winding. Winding is purely a
role-consistency contract; no geometric query here breaks if it is violated.

## signed_area, orientation, point_in_ring

**`signed_area` is approximate and kernel-free** — shoelace translated to
`vertex(0)` to cut cancellation, `0.5 * Σ cross(v_i - v_0, v_{i+1} - v_0)`. Its
**sign must never drive a topology decision**: over a sliver ring at UTM33
magnitudes the sum cancels to noise, and taking its sign is how a hole gets
classified as an outline.

**`orientation<K>` is exact and sums nothing.** It finds an extreme vertex
(min y, ties by min x, ties by lowest index — **non-normative**: any
convex-hull vertex is convex, so max-y and max-x tie-breaks are equally correct
and the choice is unobservable) and returns
`K::orient2d(prev, v, next)` — one predicate call, no accumulation, no
cancellation. If that triple is collinear it walks `prev` back and `next`
forward to the first non-collinear pair — **`prev` and `next` advance
independently, never in lockstep**. A lockstep walk returns `Collinear` for the
proper triangle `{A,A,B,C}`, where every symmetric pair about the extreme vertex
is collinear and the walk exhausts. An all-collinear ring returns `Collinear`. Exact only for a *simple* ring, which is accepted: a non-simple
ring has no well-defined winding to report.

The extreme-vertex search needs a lexicographic comparison of two `Point2`. It
is `detail::lexicographically_before` local to `ring.hpp`, and **`point.hpp`
gets no ordering and no `operator<=>`**. The comparison is only meaningful under
the finiteness precondition the ring carries and `Point2` does not; putting it
on the type makes it available exactly where it is unsafe, which is why it was
left off. This is the one place in the increment someone will be tempted to
reopen that decision.

**`point_in_ring` is three-valued and exact, even-odd:**

- p exactly on an edge, endpoints included → `Boundary`. Via `on_segment<K>` per
  edge, tested before the parity update, early return.
- p exactly on a vertex → `Boundary`, falling out of the above.
- otherwise even-odd crossing parity, half-open in y: count edge `(u,v)` when
  `(u.y <= p.y) != (v.y <= p.y)`, side decided by `K::orient2d(u, v, p)`,
  **never** by a computed x-intersection. Odd → `Inside`, even → `Outside`.
- zero-length edges contribute nothing — both endpoints fall the same side of
  the half-open test. Correct by construction, no special case.
- self-intersecting ring → the even-odd classification, which is total and
  deterministic. A *specified* answer, not unspecified. Even-odd is chosen over
  nonzero-winding precisely because it needs no consistent orientation, and a
  self-intersecting ring has none.

Three-valued rather than `bool` because CDT hole removal, DEM-border draping and
the noder all need `Boundary` distinguished; folding it into either bucket is
how boundary vertices get deleted along with their hole.

**`point_in_ring` performs no arithmetic of its own.** Every numeric decision is
a `K::orient2d` call or a comparison between two input coordinates. No division,
no intersection construction, no accumulation. Under `DefaultKernel` the result
is exact for every finite input. It is invariant under cyclic rotation and under
reversal of the ring; both are property tests.

## Tolerances: there are none

No epsilon parameter appears in any of the three headers, and none may be added
without an architectural decision. `Box2` containment, `on_segment` and
`point_in_ring` are all closed and exact.

The two tolerances in the neighbourhood answer different questions and are never
substituted for each other:

- `RasterGeometry::boundary_epsilon()` — relative 1e-12, ~8e-6 m at UTM33
  northings. Absorbs a few ulp of clipping noise against the DEM rectangle.
  `contains_strict(p)` is the authority for *"may I bilinear-sample here"*.
- The future snap grid — sized for planimetric fidelity and noding robustness,
  not tied to cell size; decimetres or centimetres for typical terrain work.
  See the snap rounding section of `parallel_refinement.md`. Orders of magnitude
  coarser than `boundary_epsilon` either way.

`point_in_ring(border_ring, p)` is the authority for *"is p in the meshing
domain"*, and has no epsilon. On a ring built from `RasterGeometry` corner nodes
its vertices are bit-identical to what `contains_strict` compares against, so
the border classifies as exactly `Boundary`.

Rule downstream must obey: a point classified `Boundary` on the DEM border ring
is not in the strict interior and must not be bilinear-sampled without a clamped
fallback. Where the snap grid and `boundary_epsilon` disagree — a vertex the
snap grid moved more than 8e-6 m across the border — **the snap grid wins**:
post-snap geometry is the geometry, and `boundary_epsilon` is noise absorption
for pre-snap clipping with no authority over a deliberate move.

Increment 2 never sees either number.

## Degeneracy policy

Checked at ring construction, throwing `std::invalid_argument`:

1. **Fewer than three vertices → rejected.** O(1). `point_in_ring` and
   `orientation` therefore never see a 2-vertex ring.
2. **Stored closure (`vertex(0) == vertex(size()-1)`) → rejected.** The
   highest-value check in the increment. Closure is implied everywhere in this
   project — increment 3 deliberately does not store the repeated index — and
   accepting a stored one means every ring has two encodings and every
   off-by-one ring bug lives in the gap. The message says closure is implied.
3. `IndexedRing` checks size ≥ 3 and closure, **plus an O(1) range check of
   exactly indices `0` and `size()-1`**. The closure check *is*
   `vertices[chain[0]] == vertices[chain[n-1]]`, so without that guard the one
   constructor documented as not validating indices performs two unchecked
   out-of-range reads. Range-checking those two costs nothing and does not touch
   the policy. The full per-index range check against `vertices.size()` remains
   O(n) and belongs to the PSLG's one-time validation, not to every view
   construction — a view may be built per query.

Accepted and fully specified:

4. **Repeated consecutive vertices.** Total everywhere: `on_segment` reduces to
   point equality, parity is unaffected, `orientation`'s collinear walk steps
   over them. Undesirable but not this layer's problem — the noder dedups.

   **Consequence: a legal ring's cyclic rotation can be illegal.** `{A,A,B,C}`
   is legal and its rotation `{A,B,C,A}` is a stored closure, which rule 2
   rejects. Both rules are individually right and the check cannot distinguish
   the cases. Any consumer that rotates a chain — normalising it to start at its
   lowest vertex, say — must collapse an adjacent duplicate first. Increment 3
   should know this before it writes one.
5. **Zero-area / all-collinear rings.** `orientation` → `Collinear`,
   `signed_area` → ~0, `point_in_ring` → `Boundary` on the spine and `Outside`
   elsewhere. Rejected at the PSLG when a role is attached, since neither role's
   winding can be satisfied.
6. **Self-intersecting rings → accepted, not detected.** No simplicity check and
   no segment-segment intersection in this increment: intersection
   *construction* rounds, rounding needs the snap grid, and the snap grid is
   increment 5's vocabulary. The noder is what makes self-intersection go away.

Not checked anywhere in increment 2:

7. **Finiteness.** Scanning n coordinates on every view construction is the
   wrong cost for a type built per query. It is a precondition of the vertex
   buffer, validated once at the PSLG and Python boundaries, exactly as
   `orientation.hpp` already states for the predicate kernel. `all_finite(ring)`
   exists so increment 3's validator and the tests do not each write their own.
   On a non-finite ring `point_in_ring` is total but meaningless and
   `bounding_box` is silently corrupt — documented, not guarded.

**What increment 2 may assume: nothing.** Every function is total on every ring
that survives construction, including non-simple ones. A layer that assumes
simplicity two increments before the noder exists is a layer tested on inputs it
will never see in production.

## Segment intersection: deliberately absent

`Segment2` gets exactly one predicate, `on_segment<K>`, and nothing else.

`on_segment` is in because `point_in_ring`'s `Boundary` answer is defined by it,
it is exact, and it constructs and rounds nothing. It is a classification, not
an intersection.

Out until the noder needs them: `segments_intersect` (a predicate with no caller
is a predicate with no meaningful test, and the noder's broad phase wants a
batched form), `intersection_point` (constructive — divides, rounds, and the
rounding target is the snap grid), distance-to-point, projection, length.
Segment intersection is the noder's defining operation and should be designed
with the snap grid in the same head.

## Risks

1. **Dangling views** — the failure mode this representation invites, mirroring
   detria's non-copying spans one layer down. Mitigated by ownership in the type
   names and `= delete`d rvalue-`vector` constructors, so the loop-local
   temporary is a compile error. Residual: a view over a `std::vector` that is
   *reallocated* while the view lives. Unguardable in C++; gets a header comment
   and an asan test.
2. **Four instantiations** of every algorithm per TU using both kernels and both
   models. Accepted, consistent with `RasterSource`, mitigated by the algorithms
   being short.
3. **`signed_area`'s sign will be misused** despite the comment. If review
   catches it once, rename to `approximate_area_magnitude`. Not pre-empting;
   the name matches every other library.

No production header needs to change. `point.hpp`, `raster/*` and `predicates/*`
are untouched. `std::hash<Point2>` is used nowhere here — nothing in this
increment is dedup-shaped.

## Documentation debt

Cleared in the governance-cleanup pass, not carried forward. `project_structure.md`,
`README.md`, `testing.md`, `CLAUDE.md` and the workflow were corrected there, and
`docs/increments/README.md` now requires a doc defect found during an increment to be
fixed in that increment's PR or not recorded — this section is the evidence for why.
