# Increment 5 — the constraint noder

Status: design settled, **and the increment is split**. The full noder — snap
rounding, broad phase, exact pairwise intersection, segment splitting, dedup,
`is_river` merge and the `NodedPslg` that carries the result — does not fit one
PR under `CLAUDE.md` §2. It ships as **5a** (this document, in full detail) and
**5b** (specified here to the level that fixes the seam, designed in full in its
own record).

Three premises inherited from the existing documents are **overturned below on
evidence**, not inherited:

1. `parallel_refinement.md:148` has the noder output "a list of `(p0, p1,
   is_river)` segments plus the deduplicated vertex array". That output cannot
   drive the CDT: `detria` needs `addOutline`/`addHole` per ring
   (`04-cdt.md:547-556`), and a flat edge list has thrown the `Outer`/`Hole`
   distinction away. The noder preserves chains and roles.
2. `03-pslg.md:65-70` and `project_structure.md`'s `noding` section both say the
   noder contributes `is_river` as "a sparse per-edge **override** set". There is
   nothing left to override: after noding an edge may descend from two chains at
   once, so it has no single source bit to differ from. Dense, one byte per
   edge.
3. `testing.md`'s `noding` invariant list has **five** bullets (lines 106-110);
   **three of them are false**, one of them false in a way that snapping
   guarantees rather than merely permits. The fifth, the `is_river` OR rule at
   line 110, is true and survives as guarantee 15. Replacements under
   "Documentation this PR fixes", covering all five.

**All three quotations above are of text that no longer stands, and that is the
point rather than a defect** — but they stopped standing in two separate waves,
which matters because only the first was increment 7's.

Increment 7 fixed items 1 and 3 at their sources in `e51ad75`:
`parallel_refinement.md`'s step 6 outputs property sets rather than "a list of
`(p0, p1, is_river)` segments", `03-pslg.md` describes a property set and the
C++ member is `Chain::properties`, and `testing.md` states the OR rule over sets
and union. Read items 1-3 as "what these documents said when this design was
written, and why", with `git show 4834568:` in front of each path to see it.

**Item 2 outlived that wave and was settled by 5b itself**, in `6f7f7c8`, when
the design was rebased onto the merged increment 7 and the conflict had to be
resolved one way or the other. `03-pslg.md` and `project_structure.md` both said
*sparse override set*; both now say **dense**, with item 2's own reason — after
noding an edge may descend from two chains at once, so it has no single source
value to be an exception to. So item 2 is no longer a premise this design
overturns against the documents; it is a ruling the documents now carry. It is
kept here because the argument is the reason, and deleting it would leave the
conclusion standing in two files with the reason in none. Only the line numbers were refreshed here — **nothing else
in this file was touched by increment 7**, and in particular the one-bit
spelling in this document's own prose is superseded in width by
`05b-noder-driver.md`, which is that branch's marking to make and not this
one's.

Nothing here re-derives an exactness claim that commit `e412a43` already
refuted. Snapped coordinates are **not** exactly representable, and the
fall-through path through `DetriaExact` is the common one on snapped breakline
data; both facts are load-bearing below rather than worked around.

### What this revision changed, and why

Revised against `docs/increments/kernel-sufficiency-audit.md` (commit
`be54efc`), which measured what the provisional draft of this document had
inherited. The audit's verdict is that **the predicates kernel is sufficient and
needs no change** — nothing in `include/terrain/predicates/` moves for the
noder. Three things follow, all of them local to this increment:

1. **The hot-pixel gap, and it is the substantive one.** `on_segment<K>`
   (`include/terrain/core/segment.hpp:68-74`) gates on
   `orient2d(...) == Collinear`, which is **exact incidence**. Snap rounding
   needs **hot-pixel proximity** — whether a segment passes *through* a cell —
   and at a decimal spacing the two differ on the great majority of snapped
   input. Every branch of the draft's `classify` that could find a T-junction
   went through `on_segment`, so T-junctions classified `Disjoint`, the host was
   never split, and guarantee 14's verification pass agreed with the error. 5a
   now ships `SnapGrid::cell_min`/`cell_max` and
   `noding::segment_meets_cell<K>` — four `orient2d` calls against a cell's
   corners plus a bbox overlap, exact, division-free, tolerance-free — and 5b's
   split pass is driven by that predicate rather than by `classify`. See
   "The hot-pixel predicate". `SegmentRelation` is unchanged.
2. **Ruling 4 rested on a false premise.** It quoted `01-predicates.md`'s
   incircle cost paragraph verbatim: "collinear-but-unresolvable triples are the common case in snapped
   breakline data". Measured, that is wrong in its most important half —
   `world(g) = g·spacing` is not affine for a non-dyadic spacing, so snapping
   does not preserve collinearity and only 207 of 2976 grid-collinear
   general-direction triples report `Collinear` at 0.1 m. What is common is
   fall-through with a **definite sign**. The ruling's *conclusion* survives
   unchanged and is better supported than before. The sentence in
   `01-predicates.md` is a documentation defect and this PR fixes it where it
   lives.
3. **Dyadic spacing becomes guidance.** `parallel_refinement.md`'s dyadic-lever
   paragraph bundled two properties and dismissed both. Only one of them needs local coordinates.
   Exact affine scaling — which is what the noder's topology depends on — needs
   **dyadic spacing alone**, measured at CRS-origin anchoring and Web Mercator
   magnitudes. New subsection "Choosing the spacing". It is guidance, not a
   default: risk 2 stands and there is still no default spacing in C++.

Revised a second time against `@reviewer`'s pass on this document. Four things
moved, all of them in the direction of not letting a check share a defect with
its subject: the split pass's candidate generator now carries a stated
conservativeness postcondition with a brute-force oracle, and guarantee 14 gains
an all-pairs check on small inputs ("the split pass"); the iteration's
termination condition is the guarantee-14 verification and explicitly not "no new
node"; the totality property's oracle no longer shares `on_segment<K>` with its
subject; and arc order is no longer described as construction-free, because it is
unfiltered `double` arithmetic and is now risk 10. The claim that 5a does not
depend on `on_segment<K>` was overstated and is corrected wherever it appeared.

Revised a third time against the same reviewer's second pass, in the same
direction. Guarantee 15's check had reproduced the self-confirming shape the
second revision removed from 14 — re-deriving `edge_is_river()` from the
provenance map the dedup built in the same traversal — and now has an oracle
computed from the input, in `prop_noding_no_crossings.cpp` beside 14's; risk 12
follows it. A citation of `parallel_refinement.md:41-51` still quoted the
unresolvable-triple claim that this document's ruling 4 corrected, and now
quotes the measured fall-through rate. Three smaller things: the endpoint-arm
analysis names where the zero-length-`s` case is handled, the endpoint arm's
unkilled `Touching` branch is listed as mutant 18, and `for_each_candidate` is
fixed as a declaration rather than described in prose.

What the audit found *survives* is not churned here: totality of `classify`,
`crossing_point` as the only constructive function in the project, the mandatory
clamp, and risk 4 as bounded-and-accepted are unchanged. The audit measured zero
non-tie cell mismatches in ~45,000 crossings and the clamp never firing, which
is stronger support than the draft had for any of them.

**5a ships** `include/terrain/core/snap_grid.hpp`,
`include/terrain/noding/intersect.hpp` and
`include/terrain/noding/node_set.hpp`. Namespaces `terrain` (the grid) and
`terrain::noding` (the rest). Header-only; nothing goes in `src/`. Depends on
increment 2 (`Point2`, `Segment2`) and increment 1 (`pred::GeometryKernel`,
`pred::Orientation`). It **does** include `core/segment.hpp` and `classify`
calls `on_segment<K>` in its endpoint arm — what 5a does not do is *detect
T-junctions or drive splitting* through that predicate, which is the narrower
and true form of the claim the provisional draft overstated; see "the hot-pixel
predicate". It depends on increment 3 for **nothing**, which is the seam — see
below.

**Invariant-critical suites (5a):**
`tests/cpp/unit/test_noding_snap_rounding.cpp`,
`tests/cpp/unit/test_noding_pairwise_intersection.cpp` and
`tests/cpp/property/prop_noding_snap_invariants.cpp`. Mutation testing applies
to those three. `tests/cpp/unit/test_noding_node_set.cpp` is not
invariant-critical. Reasons under "The invariant-critical suite".

Makes possible, at 5b: **the end of `CdtStatus::NotNoded` and of
`DuplicatePointsFound`** — increment 4's risk 1, which it names as "the
practical blocker on end-to-end use". 5a on its own makes nothing possible; it
is the arithmetic 5b is unwritable without, in the same relationship increment 1
had to increment 2.

## The split, and why the seam is here

The ceiling is per PR. Counting non-comment production lines, the whole noder is
roughly 270 (5a) + 480 (5b) ≈ 750, which is **over the ceiling outright**. The
provisional draft put it at 690 — nominally one PR, and not one in practice
because the estimate had no margin. The hot-pixel predicate added ~35 lines to
5a and ~20 to 5b and settled the question: the split is now forced by the
ceiling rather than argued past it. That is the useful reading of the audit's
one finding — the seam did **not** move, and the addition landed on the side of
the seam it belongs on, because a segment-versus-cell predicate is a question
about numbers.

More decisively, the two halves fail differently:

> **5a is everything whose correctness is a question about numbers. 5b is
> everything whose correctness is a question about topology.**

Snapping, intersection classification, intersection construction and dedup keys
are decided by signs, rounding and exact comparison; a defect in them is silent
and is killed by mutation testing on tiny fixtures. Broad phase, arc-order
splitting, chain reassembly, winding re-validation and the no-crossing
verification are decided by bookkeeping over the whole graph; a defect in them
is found by a property over generated constraint sets. Those are different test
spends, and the README's cost constraints exist to let them be paid separately.

**The seam has a mechanism, not just an argument: no file in 5a includes
`core/pslg.hpp`, and no 5a test target links anything but `terrain_headers`,
`terrain_predicates` and Catch2.** That is grep-checkable and link-checkable, in
the same spirit as increment 4's deliberately unlinked `test_cdt_backend_seam`
target. If a 5a header acquires a `Chain`, the split has been violated and the
build says so.

**What 5a does not do, and this is the honest cost:** it lands three headers with
no production caller. Only tests instantiate them until 5b merges. Increment 1
did exactly this and the alternative — dragging the driver in to give the
arithmetic a caller — is the thing the ceiling exists to prevent.

## 5a: the types

### The snap grid

```cpp
// include/terrain/core/snap_grid.hpp
namespace terrain {

struct GridPoint {
    std::int64_t ix{};
    std::int64_t iy{};

    friend constexpr bool operator==(const GridPoint&, const GridPoint&) = default;
    friend constexpr auto operator<=>(const GridPoint&, const GridPoint&) = default;
};

// |index| must not exceed this, on either axis. Derivation under "Injectivity".
inline constexpr std::int64_t kMaxGridIndex = std::int64_t{1} << 51;

// Finite and strictly positive. Pulled out as a free constexpr predicate for
// the same reason increment 3 pulled out detail::sizes_fit_u32: a check that
// can only be exercised through the object that asserts on it is a check that
// rots. 5b's driver calls this and returns a status; SnapGrid asserts it.
[[nodiscard]] constexpr bool is_valid_spacing(double spacing) noexcept;

class SnapGrid {
public:
    // Precondition: is_valid_spacing(spacing). Debug-asserted.
    explicit constexpr SnapGrid(double spacing) noexcept;

    [[nodiscard]] constexpr double spacing() const noexcept;

    // True iff snap(p) is representable and world() is injective there.
    // Precondition: p is finite -- Pslg guarantee 1, not re-checked.
    [[nodiscard]] bool can_snap(const Point2& p) const noexcept;

    // Precondition: can_snap(p). Debug-asserted, unchecked in release.
    [[nodiscard]] GridPoint snap(const Point2& p) const noexcept;

    // Total, and a pure function of g and spacing(). No precondition.
    [[nodiscard]] Point2 world(const GridPoint& g) const noexcept;

    // world(snap(p)). Same precondition as snap.
    [[nodiscard]] Point2 snapped(const Point2& p) const noexcept;

    // The closed cell (hot pixel) of g: [(ix-1/2)s, (ix+1/2)s] x likewise in y.
    // Low and high corners. Pure functions of g and spacing(), like world().
    // Precondition: |ix|, |iy| <= kMaxGridIndex. Debug-asserted.
    [[nodiscard]] Point2 cell_min(const GridPoint& g) const noexcept;
    [[nodiscard]] Point2 cell_max(const GridPoint& g) const noexcept;

private:
    double spacing_{};
};

}  // namespace terrain
```

**`SnapGrid` lives in `core/`, not in `noding/`,** and that is forced rather than
chosen. `NodedPslg` is a `core` type by the argument increment 3 made for `Pslg`
and increment 4 repeated for `IndexedMesh2` — it is the noder's output *and* the
CDT's input, so a type owned by either consumer acquires that consumer's
vocabulary. `NodedPslg` carries a `SnapGrid` by value, because
`parallel_refinement.md:134-137` observes that "coordinates lie on the snap
grid" is uncheckable unless the type carries the grid. A `core` header may not
include a `noding` header — the cardinal rule in `project_structure.md` —
therefore `SnapGrid` is `core`. One scalar, anchored at zero, exactly as that
document specifies; there is no offset member and there must never be one, since
the anchor is the entire reason `snap` is a pure function of `p`.

**`GridPoint` has `operator<=>` as well as `operator==`.** The ordering is not
geometric and means nothing spatially; it exists so `NodeSet` can sort. It is
defaulted, so it is lexicographic in `(ix, iy)` and total — `std::int64_t` has
no NaN, which is the whole reason the dedup key is integral. `point.hpp`'s
`std::hash<Point2>` is **prohibited as a dedup key** here, as
`parallel_refinement.md:126-128` rules: it cannot find a NaN and it hashes
near-coincident points apart.

### The snapping contract

```
ix = std::llround(p.x / spacing_)
iy = std::llround(p.y / spacing_)
world(g) = Point2{ g.ix * spacing_, g.iy * spacing_ }

half     = spacing_ * 0.5                       // exact: scaling by 2 is exact
cell_min(g) = Point2{ (2*g.ix - 1) * half, (2*g.iy - 1) * half }
cell_max(g) = Point2{ (2*g.ix + 1) * half, (2*g.iy + 1) * half }
```

Five rulings, each of which a reasonable implementation gets wrong:

**1. `std::llround`, never `std::rint` or `std::nearbyint`.** The latter two
honour the *dynamic* floating-point rounding mode, so the same input can snap to
different cells in two threads or two builds. `llround` is specified to round
half away from zero unconditionally. They differ at exactly the ties — `0.5`
snaps to index 1 under `llround` and to index 0 under `rint`'s round-half-even —
which is what makes the mutant killable by a single fixture rather than by luck.

**2. Divide by the spacing; do not multiply by a precomputed reciprocal.**
`p.x * (1.0 / spacing_)` is two roundings where `p.x / spacing_` is one, and the
two disagree on a set of inputs that is easy to find by search and impossible to
predict by reading. The reciprocal is the obvious optimisation and it buys a
division in a loop that is already dominated by predicate calls.

**3. `snap` is a pure function of `p` and `spacing_`, and nothing else.** No
bounding box, no raster origin, no accumulated state. That is the anchoring
decision's payoff, and `parallel_refinement.md:105-118` gives the reason: an
input-dependent anchor makes adding one vertex re-snap every other vertex, and
forces a global reduction across the front of a module designed as four
parallel-fors. The signature is the enforcement — there is no argument through
which an anchor could arrive.

**4. Exactness: `world(g)` is not exact, and that is settled, not open.** Commit
`e412a43` and `parallel_refinement.md:34-39` establish it: snapped world
coordinates are `index * spacing`, neither factor is generally a dyadic
rational, and `orient2d` then forms differences *and* products, so exactness
would need coordinate magnitudes under roughly 2^26 rather than 2^53.

The provisional draft continued from there by quoting `01-predicates.md`'s
incircle cost paragraph — "collinear-but-unresolvable triples are the common
case in snapped breakline data" — as settled evidence. **That sentence is wrong and this PR fixes it where
it lives.** The audit measured the underlying map: `world(g) = g·spacing` is not
an affine map when `spacing` is not a dyadic rational, because `fl(ix·s)`
perturbs each coordinate independently by up to half an ulp. Snapping therefore
does **not** preserve collinearity, and three points exactly collinear *on the
grid* are usually not collinear in world coordinates:

```
spacing  CRS       grid-collinear triples reported Collinear    filter fell through
                        general direction   axis-aligned
0.1 m    UTM33            207 / 2976        1024 / 1024              1517 / 4000
0.1 m    WebMerc          284 / 3003         997 /  997              1311 / 4000
0.01 m   UTM33             39 / 3010         990 /  990              1062 / 4000
2^-4 m   UTM33          3007 / 3007          993 /  993              4000 / 4000
```

(`kernel-sufficiency-audit.md` §5.1; the probe is in its appendix.) So what
snapped data actually produces is **fall-through with a definite sign** — 38 %
of triples at 0.1 m — not fall-through with a collinear answer. Axis-aligned
runs are the exception and are exactly collinear always, because one coordinate
is literally equal across the three points. An *exactly* collinear triple is in
fact the cheapest input measured, because the filter's degenerate answer
short-circuits.

**The conclusion is unchanged and is better supported than before: the design
absorbs this because `DefaultKernel` is total and decisive, not because the
fallback is rare.** It is now known to be *not* rare on a measurement rather
than on an inherited sentence. Do not write a test that asserts snapped input is
cheap, and do not reach for `FastKernel` on snapped geometry — 38 % of its
orientation queries on this data are exactly where it is permitted to be wrong.

Three places in the design depend on the corrected version rather than on the
old one, and each is handled where it lives: `classify`'s collinear arm (see
"the collinear arm"), what the suite may fixture as `Overlapping` (same place),
and the choice of spacing (see "Choosing the spacing").

**5. Cell corners are one rounding, computed in index space.** `cell_min`
and `cell_max` must be `(2*ix ± 1) * (spacing_ * 0.5)`, never
`world(g).x ± spacing_ * 0.5`. The second spelling is two roundings on top of an
already-rounded coordinate; the first is a single `fl` of an exact product of an
exact integer and an exactly-halved spacing, which makes a corner a pure
function of `g` and `spacing()` in precisely the sense `world` is. The integer
`2*ix ± 1` cannot overflow: `kMaxGridIndex` is `2^51` and `std::int64_t` holds
`2^52` with ten bits to spare. At a dyadic spacing the product is exact
outright — see "Choosing the spacing", which is one of the reasons that
subsection exists.

What snapping *does* buy is the only property the noder actually needs:

> **Two points snap to the same `GridPoint` if and only if their snapped world
> coordinates are bit-identical doubles**, because `world` is a pure function of
> `g`. Coincidence therefore has an exact key, and dedup is legal.

That is a statement about equality, not about accuracy, and it is the one to
repeat in the header.

### Injectivity, and where `kMaxGridIndex` comes from

The "if and only if" above needs `world` to be **injective** on the permitted
index range: two *distinct* grid points must not produce identical doubles, or
dedup would leave two nodes at one coordinate and hand `detria` the
`DuplicatePointsFound` this increment exists to eliminate.

For a normalized double `x`, `ulp(x) <= x * 2^-52`. The exact gap between
neighbouring lattice coordinates is `s`, and each computed `fl(k*s)` sits within
`k*s*2^-53` of its exact value, so two neighbours stay distinct when
`s > 2 * (k+1) * s * 2^-53`, i.e. `k + 1 < 2^52`.

`kMaxGridIndex = 2^51` therefore holds a full factor of two of margin over a
bound that is itself a bound. It is not tight and should not be: the practical
range is 1e8–1e9 indices at decimetre spacing over a large projected CRS
(`parallel_refinement.md:119-125`), and even a 1 mm grid over Web Mercator's
±2.0037e7 m reaches only 2e10 ≈ 2^34.2. `2^51` is 2.25e15.

This is also why the index type is `std::int64_t` and not `std::int32_t`, which
that same passage rules for a different reason (intermediate arithmetic
overflows it immediately). Both reasons are live; keep both.

**`can_snap` is the check, and it is a separate pass, not a branch inside
`snap`.** The driver runs it once over the vertex buffer and emits one
diagnostic per offending vertex, exactly as increment 3's stage 3 runs finiteness
once over the buffer rather than per ring. After that pass, `snap` is
branch-free in the parallel-for. `can_snap` is also the only place
`|p.x| / spacing` is compared against `kMaxGridIndex`; `snap` assumes it.

### Choosing the spacing: dyadic is worth more than the documents gave it

**Guidance, not a default. There is still no default spacing anywhere in C++**
(risk 2), and spacing remains a declared parameter entering at the Python
boundary. What follows is a recommendation with a measured reason attached, so
that whoever picks the number picks it knowing what it buys.

`parallel_refinement.md` offered dyadic spacing "together with local
coordinates" as the way to get "the exactness property" back; the unbundled
replacement is at lines 53-74 and this PR wrote it. Those are two
different properties and only one of them needs the local coordinates:

1. **`world` is an exact affine scaling, so grid-collinearity survives into
   world coordinates.** Needs **dyadic spacing alone**. `ix · 2⁻⁴` is exact for
   `|ix| < 2⁵³`, and Web Mercator at 6.25 cm needs only `|ix| ≈ 3.2e8 ≈ 2²⁸`.
   Measured at CRS-origin anchoring with no local coordinates anywhere:
   3007/3007 and 3016/3016 grid-collinear triples report `Collinear`, against
   207/2976 at 0.1 m. The `2^-4 m` row of ruling 4's table is this property.
2. **The determinants are exact in plain double, so the filter passes.** Needs
   coordinate magnitudes under ~2²⁶, hence local coordinates. This one stays
   dismissed, exactly as `parallel_refinement.md` dismisses it, and its absence
   costs nothing in correctness: at `2^-4 m` the filter falls through 100 % of
   the time and answers `Collinear` 100 % of the time, at 3.2 ns.

Property 1 is the one the noder's *topology* rests on, in four named places:
`classify`'s endpoint arm (which calls `on_segment<K>`, and is exact-incidence
for that reason), `classify`'s collinear arm, `segment_meets_cell`'s cell corners
(ruling 5's product becomes exact), and `Overlapping` being reachable at all on
diagonal input. Property 2 is only a performance property. **The first is free
at CRS-origin anchoring and at Web Mercator magnitudes**, which is the finding
`parallel_refinement.md` lost by bundling them.

So: prefer a negative power of two inside whatever decimetre-to-centimetre band
the input precision and boundary tolerance already imply — `2^-3 m` = 12.5 cm,
`2^-4 m` = 6.25 cm, `2^-5 m` = 3.125 cm. It costs nothing against the decimal
neighbour it replaces and it makes a class of topology answers exact instead of
perturbed.

This is not a licence to assume it. `SnapGrid` accepts any spacing
`is_valid_spacing` admits, no code branches on dyadicness, and nothing in the
suite may assert an exactness that only holds for the dyadic case — with the one
exception in "the collinear arm", which is a rule about which *fixtures* may
assert `Overlapping`.

### Pairwise intersection: classification, then construction, never together

```cpp
// include/terrain/noding/intersect.hpp
namespace terrain::noding {

enum class SegmentRelation : std::uint8_t {
    Disjoint,     // no common point
    Touching,     // exactly one common point, and it is an endpoint of at least one
    Crossing,     // exactly one common point, interior to BOTH
    Overlapping,  // collinear, sharing a sub-segment of positive length
};

// Exact and constructive-free: signs from K, betweenness from closed
// comparisons between coordinates the caller supplied. Total on every pair,
// including degenerate segments.
template <pred::GeometryKernel K>
[[nodiscard]] SegmentRelation classify(const Segment2& s, const Segment2& t);

// THE ONLY CONSTRUCTIVE FUNCTION IN THIS PROJECT.
// Precondition: classify<K>(s, t) == SegmentRelation::Crossing. Debug-asserted.
template <pred::GeometryKernel K>
[[nodiscard]] GridPoint crossing_point(const SnapGrid&, const Segment2& s, const Segment2& t);

// Does the closed cell of g meet segment s? This is the HOT-PIXEL question and
// it is NOT on_segment: on_segment asks exact incidence, this asks proximity to
// within half a cell. Exact, total, division-free, tolerance-free.
// No precondition beyond g being in range for cell_min/cell_max.
template <pred::GeometryKernel K>
[[nodiscard]] bool segment_meets_cell(const SnapGrid& grid, const Segment2& s, const GridPoint& g);

}  // namespace terrain::noding
```

**The classification/construction split is the design.** Increment 2 deferred
both `segments_intersect` and `intersection_point` with the note that
"intersection *construction* rounds, rounding needs the snap grid" and that
"segment intersection is the noder's defining operation and should be designed
with the snap grid in the same head"
(`docs/increments/02-core-geometry.md:349-350, 380-381`). Designed with it in the same
head, the answer is that they are two functions and only one of them rounds:

- `classify` is a **predicate**. It divides nothing, constructs nothing, and has
  no tolerance. Under `DefaultKernel` its answer is exact for every finite
  input, which is what makes it usable to decide topology.
- `crossing_point` is a **construction**. It divides, it rounds, and its output
  is a `GridPoint` — never a `Point2`. Returning a grid point rather than a
  world point is the type-level statement that the rounding target is the grid
  and that no un-snapped constructed coordinate may escape into the pipeline.

The three relations `Disjoint`, `Touching` and `Overlapping` need **no
construction at all**, and that is a consequence of ordering the pipeline
correctly: input vertices are snapped *before* pairwise testing, so every
endpoint is already a grid point, and every split point for a touch or an
overlap is one of the four endpoints. Construction happens once per crossing
pair and nowhere else.

#### `classify`, stated so it can be implemented

Let `o1 = K::orient2d(s.a, s.b, t.a)`, `o2 = K::orient2d(s.a, s.b, t.b)`,
`o3 = K::orient2d(t.a, t.b, s.a)`, `o4 = K::orient2d(t.a, t.b, s.b)`.

- If all four are `Collinear`, the pair is collinear: take the intersection of
  the two x-intervals and of the two y-intervals, using `std::min`/`std::max`
  and `<=` on the caller's coordinates and no arithmetic. If either interval is
  empty → `Disjoint`. If both intersections are single points → `Touching`. If
  either has positive length → `Overlapping`.
- Otherwise, if `o1 != o2` and `o3 != o4`, the segments meet at one point. It is
  `Crossing` iff none of `o1..o4` is `Collinear`; if any is, the meeting point
  is an endpoint and the relation is `Touching`.
- Otherwise, if any endpoint lies on the other segment (`on_segment<K>`, which
  is already shipped and already exact) → `Touching`.
- Otherwise `Disjoint`.

**The endpoint arm is a totality fallback, not a discriminator, and that must be
written down or it gets deleted.** Under an *exact* `K` the arm never returns
`Touching`: reaching it needs `o1 == o2` or `o3 == o4`, and both cases are
already decided. If `o1 == o2 != Collinear`, `t` lies strictly to one side of the
line through `s`, so the two closed segments share no point at all and
`Disjoint` is right; symmetrically for `o3 == o4 != Collinear`. If `o1 == o2 ==
Collinear` then both endpoints of `t` lie on the line through `s` — provided
`s` is non-degenerate; the zero-length `s` is not covered by this sentence and
is taken separately below, to the same verdict — which forces all four to be
`Collinear` — the collinear arm — unless `t` is zero-length, in which case all
four are `Collinear` again. So under `DefaultKernel` every `Touching` is decided
by the collinear arm or the four-orientation arm, and this arm returns
`Disjoint`.

It stays, for two reasons that are not style. First, `classify` is a template on
`pred::GeometryKernel` and must be total and defensible for **any** conforming
`K`, including `FastKernel`, whose wrong orientation is exactly what drops a
genuine meeting into this arm — measured, at Web Mercator magnitudes, as a
wrongly-`Collinear` answer dropping a genuine crossing there, see "Template
spend"; without it the answer there is `Disjoint`, which is
the silent wrong answer rather than the conservative one. Second, deleting it
makes the exhaustiveness of the classification rest on the case analysis in the
paragraph above rather than on the code, and that analysis is precisely the kind
of thing this project has already had to re-derive once.

Consequence for the suite, and it is a prohibition: **no fixture may assert that
a `Touching` answer came from this arm under `DefaultKernel`**, because none
does; a suite that tries will end up asserting a coin flip, exactly as a diagonal
`Overlapping` fixture would. What *does* reach the arm under `DefaultKernel` is
its `Disjoint` answer, and the zero-length-segment-off-the-other's-line fixture
in the mandatory degenerate block is the one that gets there — that fixture is
this arm's coverage, and the suite comment says so.

**Degenerate segments get no special case, and that is a ruling copied from
`segment.hpp` with its reason intact.** A `Pslg` permits a repeated consecutive
index, so a zero-length `Segment2` reaches here (`03-pslg.md` guarantee list,
"that vertices are distinct" — explicitly not promised). For `s = (p, p)` every
`orient2d(p, p, ·)` is `Collinear`, so `o1 == o2 == Collinear` always. If `p`
lies on the line through `t` the other two are `Collinear` too, the collinear arm
is taken, and the interval intersection degenerates to a point when `p` lies on
`t` and is empty when it does not — `Touching` and `Disjoint` respectively. If
`p` is off that line, `o3 == o4 != Collinear` and the endpoint arm returns
`Disjoint`, which is also right. All three are right without a branch for any of
them. `segment.hpp`'s comment says it best and it applies verbatim: "an
explicit special case here is how the degenerate answer gets subtly wrong". The
consequence for the suite is the opposite of relaxing: because there is no
branch to see in the code, the zero-length fixtures are **mandatory**, not
optional.

**The collinear arm must measure the shared extent rather than return
`Overlapping` on sight.** Two collinear segments sharing exactly one endpoint —
consecutive edges of a ring, which is the most common pair in the whole input —
are `Touching`, and splitting them at each other's endpoints is a no-op. An
implementation that reports `Overlapping` there would merge every ring edge with
its neighbour. It has its own fixture.

**But the collinear arm is near-dead code on diagonal snapped input, and the
provisional draft had this exactly backwards.** It named the defect above as
"the single most likely defect in the function". By ruling 4's measurement the
likely defect is on the other side: at 0.1 m, "all four orientations
`Collinear`" fires for axis-aligned pairs always and for **at most** ~7 % of
general-direction pairs — 207/2976 is a per-*triple* rate, and the arm needs all
four of a pair's orientation triples to report `Collinear`, so the pair rate is
lower still and 7 % is an upper bound rather than the measurement. Nothing here
extrapolates a pair figure that was never measured; the conclusion only needs the
bound. A genuinely overlapping grid-collinear diagonal pair therefore falls through to
the four-orientation arm and classifies as `Crossing` or `Disjoint` depending on
which way an ulp-level perturbation falls. **Both are reachable.** Three
consequences, and they are rulings rather than observations:

- **The suite may assert `Overlapping` only on fixtures where the arm provably
  fires**: axis-aligned pairs (one coordinate literally equal across the
  points), or any-direction pairs at a dyadic spacing. A fixture that builds a
  diagonal overlap at 0.1 m and asserts `Overlapping` is asserting a coin flip,
  and it will pass on the author's machine. Both kinds appear in the suite —
  the dyadic one to pin the arm, and a **0.1 m diagonal overlap asserted to be
  `Crossing`-or-`Disjoint`**, with a comment saying why the weaker assertion is
  the honest one.
- **`Overlapping` is no longer load-bearing for the `is_river` merge.**
  `parallel_refinement.md:176`'s road-along-a-river rule cannot depend on a
  relation that essentially never fires for a diagonal river. It does not have
  to: after 5b's hot-pixel split pass, two coincident edges are coincident as
  **node-id pairs**, and the merge is driven by exact integer edge-key dedup on
  `(min(u,v), max(u,v))`. See "the split pass". That is a stronger mechanism
  than the one the draft assumed, and it needs no collinear answer from anyone.
- **`Overlapping` stays in the enum anyway.** It is the exact truth when it is
  true, `classify` must stay total and exhaustive over the exact classification
  of two closed segments, and the chain-against-itself pass that produces
  `NonSimpleRing` reads it. Demoting a value from load-bearing to merely correct
  is not a reason to delete it.

**Both axes, not a dominant axis.** Choosing "the axis along which the segments
vary" is the natural optimisation and it needs a branch that is wrong for a
vertical pair, for a horizontal pair and for a zero-length segment. Intersecting
both intervals costs four comparisons and has no degenerate case.

#### `crossing_point`, and the clamp

```
d1 = s.b - s.a;  d2 = t.b - t.a
den = cross(d1, d2)                       // nonzero: Crossing implies non-parallel
u   = cross(t.a - s.a, d2) / den
p   = Point2{ s.a.x + u * d1.x, s.a.y + u * d1.y }
p   = clamp p into the intersection of the two segments' coordinate ranges
return grid.snap(p)
```

**The clamp is mandatory and it is not a tolerance.** `den` is near zero for a
near-parallel crossing, and the computed `u` can then be off by many orders of
magnitude, putting `p` arbitrarily far from either segment. The true
intersection provably lies in the intersection of the two segments' coordinate
ranges, so clamping each coordinate into `[max(minₛ, minₜ), min(maxₛ, maxₜ)]`
bounds the damage to a region that provably contains the answer, using only
`std::clamp` on coordinates the caller supplied. No new `Box2` member is
introduced for this — `bbox.hpp` has `intersects` but no `intersection`, and
adding one to a `core` header to serve one call site in `noding` is the wrong
direction of dependency.

**The clamp is also what establishes `snap`'s precondition**, and that is the
second reason it is mandatory rather than defensive. `crossing_point` ends in
`grid.snap(p)`, whose precondition is `can_snap(p)`; an unclamped `p` from a
near-parallel pair can be arbitrarily far out and can fail it, or — worse in
release, where the assert is gone — index past `kMaxGridIndex` and break `world`
injectivity. After the clamp, `p` lies inside the intersection of two segments'
coordinate ranges, every endpoint of which already passed `can_snap` in the
driver's one-pass check, so `can_snap(p)` holds by monotonicity of
`|p| / spacing` and no second check is needed. Remove the clamp and this function
acquires an unstated precondition on its inputs' magnitudes.

**What the clamp does not buy**, stated because the opposite is easy to assume:
it bounds the error, it does not make it small. For a near-parallel pair the
clamped point can be anywhere in a long thin overlap region, and snapping then
commits to whatever cell that is. The combinatorial answer in that case is
"whatever snap rounding says", which is the standing contract of snap rounding
and not a defect of this function. Risk 4.

**`crossing_point` takes `K` only to assert its precondition.** Its arithmetic is
kernel-free — a kernel decides signs, and this function decides a value. The
template parameter buys a debug assert that the caller classified first, which
is worth one instantiation axis in a function whose misuse is otherwise silent.

#### The hot-pixel predicate, and why `on_segment` cannot do this job

> Does the closed cell of `GridPoint g` — the square `[(ix±½)s] × [(iy±½)s]` —
> meet segment `s`?

**This is a different question from the one `on_segment` answers, and the
difference is the whole finding of the kernel audit.** `on_segment<K>`
(`segment.hpp:68-74`) returns `false` unless
`K::orient2d(s.a, s.b, p) == Collinear`: **exact incidence**. Snap rounding is
defined on **proximity**: Goodrich-Guibas-Hershberger-Tanenbaum route each
segment through every hot pixel it *passes through*, not through every pixel
whose centre it exactly contains. By ruling 4 those differ on the great majority
of snapped input, so a design that detects T-junctions through `on_segment`
reports `Disjoint` for most real ones, never splits the host, and — because
guarantee 14's verification pass asks the same question — agrees with its own
error. That is a mesh in which a breakline junction silently fails to be a
junction, which is the failure class this project spends its mutation budget on.

**The implementation, and it needs nothing the kernel does not already have:**

```
lo = grid.cell_min(g);  hi = grid.cell_max(g)
// 1. bbox reject, closed comparisons, no arithmetic
if max(s.a.x, s.b.x) < lo.x || min(s.a.x, s.b.x) > hi.x  -> false
if max(s.a.y, s.b.y) < lo.y || min(s.a.y, s.b.y) > hi.y  -> false
// 2. does the LINE through s separate the cell?
o_i = K::orient2d(s.a, s.b, corner_i)   for the four corners
if all four o_i are the same non-Collinear value          -> false
return true
```

Four `orient2d` calls and four comparisons. **Exact, total, division-free,
tolerance-free, and expressible in the kernel exactly as it stands** — no exact
construction, no new backend entry point, no change to `ExactPredicates` or
`GeometryKernel`. `include/terrain/predicates/` does not move for the noder.

Four rulings on it:

- **Both halves are necessary and neither is sufficient.** This is separating-axis
  between a segment and an axis-aligned box: the candidate axes are the box's two
  face normals (step 1) and the segment's normal (step 2). Dropping step 1 admits
  any cell the segment's *line* crosses, however far along it — the false
  positive a reader will not see in a hand-picked fixture, because a hand-picked
  fixture puts the cell near the segment. Dropping step 2 admits any cell in the
  segment's bounding box, which for a long diagonal is most of the domain.
- **"All the same non-`Collinear` value", not "all strictly positive or all
  strictly negative".** A corner exactly on the line makes one `o_i` `Collinear`,
  the cell is grazed, and the answer is `true`. Written as a strict sign test it
  becomes `false` and the predicate loses exactly the incidences it exists to
  find. This is the mutant with the most innocent-looking diff in 5a.
- **It is `noding`, not `core`.** `segment.hpp` may not acquire it: the question
  is about a snap grid, `core/segment.hpp` does not know there is one, and a
  `core` header may not include a `noding` header nor grow a `SnapGrid` parameter
  for one caller. The cell *geometry* is `SnapGrid`'s (`cell_min`/`cell_max`,
  ruling 5); the cell *predicate* is `noding`'s. That split also keeps
  `snap_grid.hpp` kernel-free.
- **`SegmentRelation` does not change, and must not.** The hot-pixel question is
  segment-versus-cell; `SegmentRelation` classifies segment-versus-segment. Risk
  8 already rules that anything further "is a different question and needs a
  different function name", and this is the first such question to actually
  arrive. A fifth enumerator, or a tolerance parameter on `classify`, is how this
  design fails.

**What this does to the pipeline is 5b's business and is fixed at the seam.**
`classify` remains the source of *new* nodes — a `Crossing` is where
`crossing_point` manufactures one — and `segment_meets_cell` becomes the source
of *incidences*: for each segment, the nodes whose cells it meets, split in arc
order. T-junctions, three-way meets and collinear overlaps all fall out of that
one pass, and none of them needs an exact-incidence answer. See "the split
pass".

### `NodeSet`: dedup by sorted grid keys

```cpp
// include/terrain/noding/node_set.hpp
namespace terrain::noding {

// Every distinct grid point in the noded input, sorted once. Node ids are
// indices into that sorted sequence.
class NodeSet {
public:
    // Sorts and uniques. Takes by value and consumes.
    explicit NodeSet(std::vector<GridPoint> points);

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::span<const GridPoint> points() const noexcept;   // sorted, unique

    // Precondition: g was among the constructor's arguments. Debug-asserted.
    [[nodiscard]] std::uint32_t id_of(const GridPoint& g) const noexcept;

    [[nodiscard]] const GridPoint& operator[](std::uint32_t id) const noexcept;

private:
    std::vector<GridPoint> points_;
};

}  // namespace terrain::noding
```

**A sorted `std::vector` with `std::ranges::lower_bound`, not a hash map**, for
the reasons increment 4 gave for `ConstraintEdgeSet` and which apply unchanged:
the key is exact and small so hashing buys nothing; one allocation of known size
beats a rehashing container; the structure is read-only after construction and
therefore shareable by every thread without synchronisation, which
`parallel_refinement.md` requires; and sorting dedups for free.
`node_set.hpp` includes no `<unordered_map>` and no `<unordered_set>`, and says
why in one line.

**Node ids are in lexicographic grid order, not first-appearance order.** This
is a ruling with a consequence. Lexicographic order is a function of the *set* of
points and of nothing else — not of input order, not of which thread found a
crossing first, not of the order the broad phase visited buckets in. So the node
numbering, and therefore every downstream index in the mesh, is reproducible
across runs and across thread counts, which is the bit-identity
`testing.md`'s parallelism section asks for. First-appearance order would make
the numbering depend on scheduling. It also matches
`parallel_refinement.md:207`'s "parallel sort + unique" exactly.

**The consequence, and it is not small: index identity with the input `Pslg` is
gone.** `Pslg` guarantee 9 promises a caller's index `k` means the same point
going out as going in, "which is what lets Python hold a parallel attribute
array". Noding must renumber — it deletes coincident vertices and inserts
crossing points — so that guarantee cannot survive, and pretending otherwise
would be the lying-type mistake increment 4 refused to make. **5b therefore
carries the forward map on `NodedPslg`:**

```cpp
// sized pslg.vertices().size(); node_of_input_vertex()[k] is the node input
// vertex k became. Many-to-one: coincident inputs share a node.
[[nodiscard]] std::span<const std::uint32_t> node_of_input_vertex() const noexcept;
```

That preserves the *intent* of guarantee 9 under an operation that cannot
preserve its letter. Named here, in 5a's record, because it is the thing a
reader of guarantee 9 will otherwise assume still holds.

## 5b: what the seam promises, in enough detail to be built against

Designed in full in `docs/increments/05b-noder-driver.md`. Fixed here:

```cpp
// include/terrain/core/noded_pslg.hpp, namespace terrain
class NodedPslg {
public:
    [[nodiscard]] std::span<const Point2>        vertices() const noexcept;
    [[nodiscard]] std::span<const Chain>         chains() const noexcept;
    [[nodiscard]] std::span<const std::uint32_t> chain_indices() const noexcept;
    [[nodiscard]] std::span<const std::uint32_t> indices_of(std::size_t c) const noexcept;
    [[nodiscard]] IndexedRing                    ring(std::size_t c) const;
    [[nodiscard]] std::size_t                    edge_count(std::size_t c) const noexcept;
    [[nodiscard]] Segment2                       edge(std::size_t c, std::size_t k) const noexcept;

    [[nodiscard]] SnapGrid                       grid() const noexcept;
    [[nodiscard]] std::span<const std::uint8_t>  edge_is_river() const noexcept;
    [[nodiscard]] std::span<const std::uint32_t> node_of_input_vertex() const noexcept;
private: /* private constructor, the noder its only friend */ };
```

**Same shape as `Pslg`, deliberately** — same accessors, same chain/role
structure, same closure-is-implied convention — so increment 4's wrapper changes
`const Pslg&` to `const NodedPslg&` and nothing else, which is the "mechanical"
change `03-pslg.md:503-506` and `04-cdt.md:659-660` both booked. It is **not** a
subclass of `Pslg` and there is no conversion between them: the whole point is
that the two types promise different things, and an implicit conversion would
let un-noded input reach the CDT through a signature that says it cannot.

`NodedPslg` adds these guarantees on top of `Pslg`'s 1–8 and 10 (9 is replaced by
`node_of_input_vertex`, above):

11. Every coordinate is `grid().world(g)` for some `GridPoint` with
    `|ix|, |iy| <= kMaxGridIndex`. Checkable: `grid().snapped(v) == v`
    element-wise, bitwise.
12. No two vertices are equal — dedup, and injectivity above.
13. No edge is zero-length; no chain has a repeated consecutive index.
14. **Two clauses, and they need two different checks.** (a) No two edges cross
    in their interiors — for every pair of distinct output edges `classify<K>`
    returns `Disjoint` or `Touching`, **or** returns `Overlapping` and the two
    edges have equal node-id pairs `(min(u,v), max(u,v))`. *(Amended by
    `docs/increments/05b-noder-driver.md`. The clause as first written —
    `Disjoint` or `Touching`, full stop — is false on the one input guarantee 15
    exists for: a road noded along a river leaves two chains carrying the same
    edge, and `classify` on two identical segments returns `Overlapping`. A
    partial overlap is still a violation, and after the split pass there is no
    third case, since each of two overlapping collinear edges is split at the
    other's nodes.)* (b) No node's **cell** meets an edge it is not an
    endpoint of — `segment_meets_cell<K>` is false for every (node, edge)
    candidate pair with that node not an endpoint. This is the promise the
    type's name makes, and 5b establishes it by **verification** — a second
    broad-phase pass over the output — not by construction.

    Clause (b) is the hot-pixel form and it replaces the draft's "no vertex lies
    in the interior of an edge", which was stated in terms of exact incidence
    and was therefore satisfiable by a mesh full of surviving T-junctions: the
    check would have asked the same question the splitting asked and agreed with
    it. In the hot-pixel form the check is strictly stronger than exact
    incidence, so a T-junction that the split pass missed is caught here rather
    than blessed. It is also exactly the fixpoint condition the iteration below
    converges to, which is not a coincidence — it is the condition classical
    snap rounding is defined by.

    **The verification is only worth what its candidate generation is worth, and
    that is the remaining self-confirming hole.** Splitting and verifying both
    draw candidates from `broad_phase.hpp`. A broad-phase defect that drops a
    (segment, node) pair — bucketing by endpoints instead of by every cell a
    segment touches, or querying a node by its centre instead of by its cell's
    `Box2`, both of which the broad-phase section below has to forbid explicitly
    — makes the split miss that pair **and** makes the verification miss it, and
    the mesh then verifies clean with a live T-junction in it. That is the
    defect this increment fixed at the predicate layer, one layer further down.
    It is closed at the seam, in two parts, below: a conservativeness contract on
    the broad phase that can be checked without the noder, and a brute-force
    verification of 14 itself on small inputs.
15. `edge_is_river()` is index-aligned with the flat edge enumeration
    `(c, k), k < edge_count(c)`, and each bit is the OR over every input chain
    that contributed geometry to that edge. **Verified, like 14, against an
    oracle built from the input rather than from the bookkeeping it checks.**

    The obvious check is the wrong one and must not be built: having 5b carry
    step 5's provenance map into the verification pass and re-derive
    `edge_is_river()` from it is the self-confirming shape of clause 14, one
    layer further up and with no oracle at all. The dedup produces the bit *and*
    the map in a single traversal, so re-deriving one from the other shows only
    that the OR was applied consistently to the contributors the dedup
    **recorded**. A dropped contributing chain — the defect this guarantee
    exists to exclude, and the road-along-a-river case at
    `parallel_refinement.md:176` that `is_river` turns on — is absent from both,
    and a mesh that has silently lost a river bit verifies clean. Naming a check
    that cannot fail is worse than naming none, because it converts "unchecked"
    into "apparently checked".

    The check is therefore `prop_noding_no_crossings.cpp`'s, same file, same
    round and same pattern as 14's all-pairs check — see "the seam" below for
    the assertion.

**Guarantee 14 is checked, not assumed, because snap rounding can create
crossings that were not in the input.** Two segments that pass within less than
a cell of each other can snap onto each other; that is the cascading case
Halperin & Packer's iterated snap rounding addresses
(`parallel_refinement.md:139`). 5b's policy: iterate the node-split-dedup round
**until the guarantee-14 verification passes**, with a cap, and report
`NodeStatus::NotConverged` at the cap rather than emit crossing edges under a
type that promises none. One pass is the common case; the cap exists so the
failure is a status instead of a hang.

**"No new node" is not the same fixpoint and must not be used as the test.** A
round can produce no node and still change the graph: if node `g` already exists
and its cell meets an edge `e` that `g` is not an endpoint of, the round splits
`e` at `g` — a new *edge*, no new node — so a loop spelled "until a round
produces no new node" can report converged in a state that violates 14(b). That
is the surviving-T-junction failure arriving through the loop condition instead
of through the predicate. The verification is the condition; node counts are a
progress metric and nothing more.

**Superseded in its representation, not in its reasoning:** an edge carries a
**property set** merged by **union**, not one bit merged by OR — at a coarse
resolution the same segment can be both a road and a river, and the face-based
fallback that excused dropping non-river semantics works only for area features.
Everything below about *density*, about the merge being keyed on node ids and
about there being no single source bit to override survives verbatim and is
strengthened; only the width changes. `docs/increments/05b-noder-driver.md`,
"Edge properties", carries the ruling and the scope.

**The `is_river` representation is dense, one `std::uint8_t` per edge**, which
overturns `03-pslg.md:65-70`'s "sparse override set" on the terms that document set
for it. An override set is sparse only if most edges agree with a single source
chain. After noding an edge can descend from two chains — a road snapped onto a
river is exactly `parallel_refinement.md:176`'s example, and the merged edge
keeps `is_river = true` while the road is geometrically forgotten — so there is
no single base bit for an override to be an exception to. A dense byte array is
one allocation, index-aligned with the edge enumeration the CDT wrapper already
walks, and testable by reading it.

**The merge is keyed on node ids, not on `Overlapping`.** Two edges are the same
edge iff their sorted node-id pairs `(min(u,v), max(u,v))` are equal — exact
integer comparison over `NodeSet` ids, with no geometry in it at all. The draft
implicitly read the merge off `classify` returning `Overlapping`, and ruling 4
shows that relation essentially never fires for a diagonal river at a decimal
spacing. It does not need to: the road's vertices snap into the river's hot
pixels, the split pass gives both chains the same nodes, and the duplicate edges
then collapse on an integer key. The geometry decides *where the nodes are*; the
integers decide *which edges are the same edge*. That is the general shape this
module should keep.

### The split pass, and what drives it

5b's step 4 — "split each input segment at all intersections lying on it, in
arc-order" — is **driven by `segment_meets_cell`, not by `classify`**. Fixed
here because it is the seam's most load-bearing consequence of the audit:

1. Broad phase → candidate segment pairs.
2. `classify<K>` per pair. A `Crossing` calls `crossing_point` and contributes a
   **new** node. `Touching` and `Overlapping` contribute nothing new — every
   split point for them is already one of the four endpoints, which are already
   nodes. This is the only step that *manufactures* coordinates.
3. `NodeSet` over every snapped input vertex and every constructed crossing.
4. Broad phase again, now segments against **nodes**. For each segment, the
   nodes `g` with `segment_meets_cell<K>(grid, seg, g)`, sorted in arc order
   along the segment, are its split points. Split there.
5. Dedup edges on `(min(u,v), max(u,v))`; OR the `is_river` bits.
6. Verify guarantee 14, both clauses, over candidates from the same index.
   Verification passes → done; otherwise iterate from 2, up to the cap.

**The seam fixes how step 6 is allowed to be trusted, because steps 1/4 and step
6 share a candidate generator.** Three requirements, and 5b may not weaken any
of them:

- **`broad_phase.hpp` carries a conservativeness postcondition, and it is stated
  on the query rather than on the noder.** The index exposes exactly one
  primitive, fixed here beside `NodedPslg` and for the same reason — guarantee
  14 rests on it, so it is the same category of binding commitment as the type
  it produces:

  ```cpp
  // include/terrain/noding/broad_phase.hpp, namespace terrain::noding
  // Visits every indexed segment whose closed bbox intersects q, and may visit
  // others. F is invoked as f(std::uint32_t segment_index).
  template <class F> void for_each_candidate(Box2 q, F&& f) const;
  ```

  Its contract is: *the visited set contains every segment whose closed bounding
  box intersects `q`*. False positives are free; a false negative is a bug in a
  function whose entire job is not to have them. Both passes go through that one
  call, with the caller computing `q` — an edge's bbox at step 2, a node's
  `cell_min`/`cell_max` box at step 4 and at 14(b). The contract mentions no
  segment pair, no node and no snap grid, so it has a **brute-force oracle that
  does not use the index**: for
  generated segment sets and generated query boxes, compare the visited set
  against all-pairs bbox intersection. `prop_noding_broad_phase.cpp`, 5b,
  invariant-critical. This is the part that is closed by construction — there is
  one function to get right, and its correctness is checkable without noding
  anything.
- **Guarantee 14 is additionally verified all-pairs, on small inputs, in the
  property suite.** `prop_noding_no_crossings.cpp` (the name `testing.md:301`
  already reserves) runs the full noder on generated constraint sets of a few
  dozen segments and then checks both clauses of 14 by **brute force over every
  (edge, edge) and (node, edge) pair**, with no broad phase anywhere in the
  check. A broad-phase defect and a matching verification blind spot cannot
  cancel against an oracle that does not share the generator. It is O(n²) and
  that is why it is small-input and a property rather than the production path:
  production keeps the indexed verification, and this property is what says the
  indexed verification agrees with the truth.
- **Guarantee 15 is verified in the same property, from the input, and never
  from the noder's provenance map — in the relation the producer used.** For
  every output edge `e` with node-id pair `(u, v)`: if some segment `(a, b)` of
  some input chain with `is_river` satisfies
  `segment_meets_cell<K>(grid, {snap(a), snap(b)}, u)` **and** the same against
  `v`, **and** `u` and `v` both lie between that segment's own snapped endpoints
  in the arc order of step 4 — the dot-product ordering, excluding a chain that
  passes near both nodes without spanning them — then `edge_is_river()[e]` must
  be true. No edge keys and no bookkeeping from step 5 anywhere in the check; the
  oracle borrows a *predicate*, the same standing as 14(b)'s all-pairs check,
  which already brute-forces `segment_meets_cell`. Spanning-plus-proximity *is*
  the contribution relation, and exactly because clause 14(b) holds: an output
  edge that a contributing segment was not split at would put a node strictly
  inside that segment. The same property checks 14 first and by brute force, so
  the implication the oracle rests on is checked rather than assumed.
  Invariant-critical, 5b.

  **Stated as collinearity plus interval containment this oracle would be red on
  correct output**, and the dual of the self-confirming invariant is worth naming
  once: a check that cannot *pass*. The split pass does not place nodes on
  segments, it places them where `segment_meets_cell` holds — hot-pixel
  proximity, up to the `h/√2` risk 9 accepts. Exact incidence would demand
  `Collinear` at **both** endpoints of a split edge and, by ruling 4's
  measurement, get it in ~7 % of cases; unsplit edges would verify trivially
  (`u = snap(a)`, `v = snap(b)`) while split edges failed — and split edges are
  the merged road-along-a-river case that is the only work guarantee 15 does.

  **Only that direction is asserted, not `iff`.** The hot-pixel relation is a
  proximity test, so a river running near `u` and near `v` and spanning them in
  arc order satisfies it without having contributed to `e`; betweenness narrows
  that window and does not close it, because the river need only run *near* the
  edge, not along it. Asserting the converse would fire on correct output, which
  is the same gap **hot-pixel dominance** records below and for the same reason:
  the converse does not hold and the property must not assert it. A spuriously
  set `is_river` bit is therefore not caught here, and 5b owes a mutant for the
  dedup's OR rather than a stronger oracle: strengthening this one back to exact
  incidence is the failure mode above.

Between them the answer to "who checks the checker" is a function with an
independent oracle plus end-to-end checks whose oracles are built from the input
rather than from the noder's own records. What is *not* claimed: that the
production verification is correct by construction. It is
correct conditional on a contract that is separately falsifiable, which is the
strongest available statement short of making the production path quadratic.

Step 4 is where T-junctions, three-way meets, collinear overlaps and
snap-induced incidences are all handled, by one predicate, with no exact
incidence anywhere. **`on_segment<K>` does not find a single incidence and does
not split a single segment**: no step of the split pass calls it, and no
node-versus-segment question in this module is asked through it. It is still
called, once, by `classify`'s endpoint arm — see "the endpoint arm is a
totality fallback", which rules on exactly what that call is worth and why it
stays.

The review grep is therefore not "does `on_segment` appear" — on a correct
implementation it does — but **`segment_meets_cell` must not be implemented in
terms of `on_segment`**, which is mutant 16 stated as a diff a reviewer can look
for.

Arc order is **not** construction-free, and calling it so was wrong. The nodes
met by a segment are ordered by the dot product of `(world(g) - s.a)` with
`(s.b - s.a)`: that is two subtractions and two products per coordinate in
`double`, on the caller's coordinates, outside `K` — unfiltered floating-point
arithmetic, at odds with `01-predicates.md:294`'s "signs only" which risk 4 cites
approvingly. It produces no new *point*, which is all the draft's phrasing was
entitled to claim, and the distinction matters because the error is real:

> **Named risk (risk 10): two nodes whose true projections onto `s` differ by
> less than the rounding error of the computed dot product can be ordered
> wrongly, and a non-monotone arc order yields a split chain that doubles back on
> itself — which guarantee 14 does not catch**, because a doubling-back chain can
> be free of interior crossings and of node-cell incidences alike.

**Accepted, with the bound stated.** The computed dot product carries a relative
error of a few ulps of the largest intermediate, so two nodes are at risk only
when their projections are within roughly `2^-50` of the segment's squared length
apart. The nodes on a segment are distinct grid points, hence at least one
`spacing` apart in some coordinate; their projections coincide to that precision
only when the segment is near-perpendicular to the line joining them, which is
the tie case below and is broken deterministically anyway. The residual is a pair
that is *near*-tied rather than tied, and for such a pair either order puts both
nodes on the chain between the same neighbours at sub-cell separation. That is
inside the `h/√2` deformation risk 9 already accepts, so it buys no new failure
mode — but it is a floating-point risk and it is now named as one rather than
denied.

Routing it through a sign is possible and is deliberately *not* done: it would
need an exact dot-product sign, i.e. a new kernel entry point, and the audit's
verdict is that `include/terrain/predicates/` does not move for the noder. 5b
revisits it only if a fixture shows a non-monotone order in practice.

**Exact ties are possible** — two distinct nodes within half a cell of each
other, on a line near-perpendicular to the segment, can project equally — and
they are broken by `GridPoint`'s defaulted `operator<=>`. That tie-break is
arbitrary geometrically and that is fine; what it must be is *total and
independent of input order*, which is the same reproducibility requirement node
ids carry, for the same reason.

## Degeneracy and failure policy

What the noder **fixes**, what it **rejects**, and what it **passes through**.
The first column is the condition; the fourth names what increment 4's
`CdtStatus` did with it before the noder existed.

| Input condition | Noder | Mechanism | Was, at increment 4 |
|---|---|---|---|
| Two vertices with identical coordinates | **fixes** | dedup: one cell, one node | `DegenerateGeometry` (`DuplicatePointsFound`) — increment 4 risk 1 |
| Two vertices within one cell but not identical | **fixes** | same | `Ok`, and a sliver in the mesh |
| Repeated consecutive index in a chain | **fixes** | the zero-length edge collapses at dedup and is dropped | `DegenerateGeometry` (`PolylineDuplicateConsecutivePoints`) |
| T-junction: a vertex near another constraint's interior | **fixes** | `segment_meets_cell` + split the host at that node — **not** `Touching`, which is exact incidence and mostly false here | `NotNoded` (`PointOnConstrainedEdge`) |
| Crossing constraints | **fixes** | `Crossing` + `crossing_point` + split both | `NotNoded` (`ConstrainedEdgeIntersection`) |
| Collinear overlap (road along a river) | **fixes** | both split at each other's nodes by `segment_meets_cell`; the duplicate edges dedup on their node-id pairs and the bits OR. Does **not** depend on `classify` reporting `Overlapping` | `Ok`, silently double-constrained |
| Three or more constraints meeting at a point | **fixes** | implicitly: all snap to one cell (`parallel_refinement.md:199`), and every segment through that cell is split at it | `DegenerateGeometry` |
| Hole touching its outer ring at one shared node | **passes through** | already one node after dedup; increment 4 measured 5 interior triangles and pinned it | `Ok` |
| Collinear redundant vertices within a ring | **passes through** | not a defect; increment 4 measured 3 interior triangles | `Ok` |
| A breakline that is a closed polyline | **passes through** | open chain, coincident ends; legal per `Pslg` stage 4's breakline exemption | `Ok` |
| Vertex outside the outer ring, or inside a hole | **passes through** | it still snaps, still dedups, still gets a node | `Ok` |
| Self-intersecting **breakline** | **fixes** | split at its own crossings like any other pair | `NotNoded` |
| Self-intersecting **closed ring** | **rejects**, `NonSimpleRing` | see below | `NotNoded` |
| A ring left with fewer than 3 distinct nodes after dedup | **rejects**, `RingCollapsed` | the spacing is too coarse for that ring | `Ok`, or `DegenerateGeometry` |
| A ring whose winding is no longer its declared role after snapping | **rejects**, `RingDegenerateAfterSnap` | see below | `Ok`, with a hole meshed as an island |
| A coordinate whose index exceeds `kMaxGridIndex` | **rejects**, `CoordinateOutOfRange` | `can_snap`, one pass over the vertex buffer | `Ok`, or nonsense |
| A spacing that is zero, negative or non-finite | **rejects**, `InvalidSnapSpacing` | `is_valid_spacing`, once | n/a |
| Snap-induced crossing that the iteration cap does not clear | **rejects**, `NotConverged` | guarantee 14's verification pass | n/a |
| Hole outside every outline; outline inside a hole | **passes through** | see "Nesting" | `InvalidTopology` |
| Non-finite coordinate | unreachable | `Pslg` guarantee 1 | `MalformedInput` self-check |

**Two rejections are new failure modes that the noder *creates*, and they are the
honest price of snapping.** They must be in the table and in the status enum
rather than discovered:

- **`RingCollapsed`.** A hole smaller than a cell becomes fewer than three
  distinct nodes. `Pslg` guarantee 4 said `count >= 3`; after dedup that can be
  false. The operational fix is a finer spacing, which is why the status names
  the ring and the message carries the spacing.
- **`RingDegenerateAfterSnap`.** `Pslg` guarantee 6 says every `Outer` ring is
  counterclockwise *under the kernel that validated it* — and that was true of
  the **pre-snap** coordinates. Snapping moves every vertex by up to half a cell
  diagonal, which can flatten or reverse a sliver ring. So **the winding must be
  re-evaluated after snapping**, and — following increment 3's stage 5 rule
  exactly — it is **never silently reversed**. This consequence appears in none
  of the existing documents and is the reason `NodedPslg` re-states guarantee 6
  rather than inheriting it.

**Self-intersecting closed rings are rejected, not repaired.** Splitting the
edges of a figure-eight ring produces a ring that still visits a node twice; its
"interior" is not well defined, so neither is its winding or its role.
Resolving it into simple components is polygon repair, which belongs upstream —
`parallel_refinement.md:198` already says "flag for resimplification or
auto-fix" and the flag is this status. Detect (the same machinery: a chain
paired against itself), reject, name the chain. **No severity field, no
warning-level diagnostic**: increment 3 risk 2 warns that a severity field is how
a validator starts having opinions, and that ruling stands.

**Nesting is still not computed.** `03-pslg.md:533-539` reserved
`nesting_forest<K>(const NodedPslg&)` for the increment where non-crossing input
makes the answer meaningful — which is this one. It is still deferred, and now
for a better reason than "the input might cross": `detria` computes nesting
itself and diagnoses the failures as `InvalidTopology` (increment 4's table,
three rows), so a second implementation here would be derived data with no
unique consumer. Named so it is not re-derived a third time.

### What replaces `NotNoded` and `DuplicatePointsFound`

After 5b, the CDT's entry point takes a `NodedPslg`. Then:

- `PointOnConstrainedEdge` and `ConstrainedEdgeIntersection` — guarantee 14 makes
  both unreachable, so `CdtStatus::NotNoded` **joins `MalformedInput` as a
  self-check**: it firing means the `NodedPslg` did not come from the noder or
  the wrapper mis-built a span. `04-cdt.md:676-680` predicted exactly this and it
  lands unchanged.
- `DuplicatePointsFound` — guarantee 12 plus `world` injectivity make it
  unreachable. `PolylineDuplicateConsecutivePoints` — guarantee 13. So
  `CdtStatus::DegenerateGeometry` loses both of its reachable rows and becomes a
  self-check too; `AllPointsAreCollinear` was already unreachable.
- **Reachable from a `NodedPslg`: `Ok`, `InvalidTopology`, `BackendFailure`.**
  Nothing else. 5b appends one paragraph saying so beneath `04-cdt.md`'s mapping
  table rather than rewriting the table, so the record of what was true at
  increment 4 stays legible.

Note the shape of this: the failures do not disappear, they **move**. Four CDT
statuses become three noder statuses (`RingCollapsed`,
`RingDegenerateAfterSnap`, `NonSimpleRing`) plus three that are about the grid
rather than the geometry. A claim that the noder makes the input "clean" would be
false; what it does is move every remaining failure to the stage that can say
something actionable about it.

## The broad phase: 5b, and it is not the snap grid

**In scope for the increment, not for 5a.** It is graph-level machinery by the
seam above, it has no meaning without segments to bucket, and it is the one part
whose cost model matters at scale.

**It is a spatial index and nothing more, and its spacing has no relationship to
the snap grid's.** `parallel_refinement.md:143` already says so; here is the
number that makes it non-negotiable. At 5 cm — the lower end of the recommended
decimetre-to-centimetre band, and the spacing the figure below is computed at — a
100 km domain has a snap lattice of 2e6 × 2e6 = **4e12 cells**. At a decimetre it
is 1e12, which does not change the conclusion by anything that matters.
A broad phase bucketed at that spacing would allocate one bucket per cell to hold
a few thousand segments. The two grids answer different questions: the snap grid
quantises *coordinates* and is sized by input precision and boundary tolerance;
the broad phase filters *candidate pairs* and is sized by segment density, which
is what makes the raster's own grid a convenient reuse rather than a required
one. Conflating them couples pair-testing cost to coordinate precision, which is
exactly backwards — refining the snap grid to reduce collapse would then explode
the index.

Consequences for 5b's design, stated now so they are not litigated later: the
broad phase takes its cell size from segment count and domain extent, not from
`SnapGrid::spacing()`; `broad_phase.hpp` does not include `snap_grid.hpp`; and a
pair that shares no bucket is never classified, so bucketing must be by every
cell a segment *touches*, not by its endpoints.

**That last one is not a style preference, it is the conservativeness
postcondition** the split pass's step 6 fixes at the seam: `for_each_candidate`
must visit every segment whose closed bbox meets the query box. Bucketing by
endpoints violates it for any segment longer than a bucket, and the violation is
invisible end to end, because the same omission silences guarantee 14's
verification. The contract is stated on the query so that it has an oracle — see
"the split pass".

**It answers two queries, not one.** Step 2 of the split pass asks it for
segment-versus-segment candidates; step 4 asks it for segment-versus-**node**
candidates, which is a point query against the same buckets and needs no second
structure. Note the one place the two grids do touch, and it is a containment
not a coupling: a node's hot pixel has extent `spacing`, so a node query must
return every broad-phase bucket its *cell* overlaps, not only the bucket its
centre falls in. Since the broad-phase cell is orders of magnitude larger than
the snap cell, that is almost always one bucket and never more than four — but
"almost always" is not "always", and the boundary case is exactly a node on a
bucket edge, which is where constraints congregate. `broad_phase.hpp` therefore
takes the query extent as a `Box2`, not a `Point2`, and stays ignorant of where
that box came from. A `Point2` overload must not exist: it is the second live way
to break conservativeness, and an overload that compiles is an overload someone
calls at the one site where the node's cell straddles a bucket edge.

## Files and LOC

The unit is **non-comment production lines**, `CLAUDE.md` §2's unit. These
headers carry roughly as many comment lines as code, by house style, so raw
`wc -l` will be about double.

**5a — this PR:**

| File | Contents | Est. LOC |
|---|---|---|
| `include/terrain/core/snap_grid.hpp` | `GridPoint`, `kMaxGridIndex`, `is_valid_spacing`, `SnapGrid` incl. `cell_min`/`cell_max` | ~80 |
| `include/terrain/noding/intersect.hpp` | `SegmentRelation`, `classify<K>`, `crossing_point<K>`, `segment_meets_cell<K>` | ~135 |
| `include/terrain/noding/node_set.hpp` | `NodeSet` | ~55 |

**~270 production LOC**, up from the provisional draft's ~235: `cell_min`/
`cell_max` are ~10 and `segment_meets_cell` is ~25, both of them the audit's
one finding. Nothing in `src/`. No existing header changes — 5a is purely
additive, which is part of why it is a safe first half. Note what did *not*
happen: `include/terrain/predicates/` is untouched, which is the audit's verdict
expressed as a diff.

Non-C++ changes in the same PR, listed separately because they are not
production code:

- `tests/cpp/CMakeLists.txt`: four suite registrations. `test_noding_node_set`
  names no kernel and uses the plain `add_terrain_test`; the other three use
  `add_terrain_backend_test`, which is the existing helper that links
  `terrain_predicates`.
- `tools/check_detria_boundary.py`: **no change**. 5a includes no detria and
  must not; the script's `include/` rule (zero occurrences, ever) already covers
  the new headers.
- The documentation fixes listed at the end.

**5b — the next PR**, estimated so the seam is honest rather than optimistic:

| File | Contents | Est. LOC |
|---|---|---|
| `include/terrain/core/noded_pslg.hpp` | `NodedPslg` and its accessors | ~120 |
| `include/terrain/noding/broad_phase.hpp` | uniform bucket index over segments, one `for_each_candidate(Box2, F&&)` query | ~100 |
| `include/terrain/noding/node.hpp` | `NodeStatus`, `describe`, `NodeOptions`, `NodeOutcome`, `node<K>`, the hot-pixel split pass and edge-key dedup | ~250 |
| `include/terrain/cdt/*`, `src/cdt/detria_backend.cpp` | `const Pslg&` → `const NodedPslg&` | ~10 |

**That last row is wrong and the whole 5b estimate with it.** `bindings/core.cpp:466`
calls `terrain::cdt::triangulate<DetriaBackend>(pslg, options)` on a
`const Pslg&`, so retyping the entry point breaks the Python extension, and the
only repair is to bind a producer of `NodedPslg` — which drags
`src_python/tin_engine/_core.pyi`, `src_python/tin_engine/cli.py`, the Python
suites and the renderer's scene join with it. No Python file appears in the
table above. `docs/increments/05b-noder-driver.md` re-estimates the whole of 5b
at ~786 non-comment production lines and splits it; the C++ figures below stand,
what was missing is everything on the other side of the binding.

**~480 production LOC**, header-only (the driver is a template on `K`, as
increment 3's validator is, so nothing goes in `src/noding/`). Under the
ceiling, without the margin 5a has. The ~20 over the draft is step 4 of the
split pass and guarantee 14's second clause; the further ~10 is the broad
phase's single-query surface and its stated postcondition.

**The seam does not move if the arithmetic moves.** 270 + 480 ≈ 750 is over the
ceiling, but the split is not being justified by 5.7 % — if a re-measure came in
at 690 the split would still stand, because the two halves fail differently, are
tested differently and are reviewable independently, which is the argument the
section above makes and the only one that does not depend on an estimate. Note
also that 5b's estimate contains a header nobody has designed in full yet
(`noded_pslg.hpp` at ~120), so the honest direction of travel for that number is
up. Nobody should relitigate the seam on a line count in either direction.

**Contingency split of 5b, dependency-ordered**, if it overruns. **Superseded:
`docs/increments/05b-noder-driver.md` splits 5b at a different seam — all of the
C++ in 5b, the signature change and the Python crossing in 5c — on the measured
grounds that C++ estimates here have come in at or under and binding estimates at
twice. The two below are recorded as what was proposed, not as a live
alternative.** The overrun cause anticipated was `describe` growing one
`std::format` call per status enumerator, which is exactly what nearly split
increment 3 and increment 4:

- **5b** — `core/noded_pslg.hpp` and `noding/broad_phase.hpp`, plus their suites.
  `NodedPslg` would then need a producer to exist at all, so the split point is a
  `NodedPslgBuilder` with the invariant checks and no noding — which is a real
  type with a real job (it is what verifies guarantee 14) and not a stub.
- **5c** — `noding/node.hpp`, the driver, and the CDT signature change.

That is the same seam increment 3 identified and did not need. It is more likely
to be needed here.

## The invariant-critical suite (5a)

**`tests/cpp/unit/test_noding_snap_rounding.cpp` — invariant-critical, mutation
round.** The arithmetic is three lines and every defect in it is silent.

- The tie: a coordinate exactly at a cell midpoint snaps away from zero on both
  sides of the origin. This is the fixture that separates `llround` from `rint`,
  and it is the only one that does.
- Purity: `snap(p)` for a point computed alone equals `snap(p)` computed with
  any other points in hand — asserted by calling it twice through different
  paths. Weak by construction, because the signature already forbids the
  alternative; kept because the *mutant* lives in 5b's driver, and a reader
  needs to find the statement of the property somewhere.
- Injectivity at the boundary: `world` of two adjacent grid points at
  `kMaxGridIndex` gives distinct doubles; at a deliberately over-large index —
  constructed directly, not via `snap` — they collide. The second half is what
  makes `kMaxGridIndex` a number with a reason rather than a number.
- `can_snap` at, just under and just over the bound, on each axis independently,
  all four corners, in the style increment 3 used for `sizes_fit_u32`. Cheap
  here, because unlike that check this one needs no 64 GiB.
- `is_valid_spacing` on `0.0`, `-0.0`, a negative, an infinity, a NaN, a
  subnormal and a normal value. `-0.0` is the one that gets written wrong:
  `spacing > 0.0` is false for it, which is the right answer, but
  `spacing != 0.0` is also false, and an implementation reaching for the second
  spelling passes every other case.
- **Cell corners, ruling 5.** `cell_min`/`cell_max` bracket `world(g)` and
  abut: `cell_max(g).x == cell_min(g + (1,0)).x` bitwise, which is what makes
  the cells a partition rather than a set of squares with gaps. Then the mutant
  fixture: a spacing and index where `(2*ix - 1) * (spacing * 0.5)` and
  `world(g).x - spacing * 0.5` differ in the last bit. Found by seeded search
  over the property generator at 0.1 m, not by hand, and the test comment says
  so — the same provenance rule as the reciprocal mutant below.
- **Not asserted: that a snapped coordinate is exactly `ix * spacing` as a real
  number.** `e412a43` settled that it is not. A test asserting it would fail, and
  a test asserting it "approximately" would be a tolerance, which this project
  does not have. The same goes for cell corners: at a non-dyadic spacing they are
  `fl` of the exact corner and nothing may assert otherwise. A **dyadic** fixture
  may assert corner exactness, and should, because that is the payoff "Choosing
  the spacing" claims and an unasserted claim is a claim nobody checked.

**`tests/cpp/unit/test_noding_pairwise_intersection.cpp` — invariant-critical,
mutation round.** One named test per `SegmentRelation`, plus the pairs that
distinguish the arms:

- Two ring edges sharing exactly one endpoint, collinear → `Touching`, **not**
  `Overlapping`. The single most likely defect.
- A crossing at a shared endpoint → `Touching`, not `Crossing`. The `Collinear`
  check inside the four-orientation arm is what decides this.
- A T-junction: one endpoint interior to the other segment → `Touching`.
- Collinear, overlapping over a positive length, in all four
  containment arrangements (partial each way, one inside the other, identical).
- Collinear, disjoint, sharing neither endpoint.
- Parallel, non-collinear, with overlapping coordinate ranges — the pair that
  passes an interval test and must fail the orientation test.
- **Vertical and horizontal collinear pairs**, which is what the both-axes
  ruling exists for and what a dominant-axis implementation gets wrong.
- **Zero-length segments**: against a segment that contains the point
  (`Touching`), against one that does not (`Disjoint`), and two zero-length
  segments equal and unequal. Mandatory, because there is no branch in the code
  to point at.
- `classify` is symmetric: `classify(s, t) == classify(t, s)` for every fixture,
  and `classify(s, t) == classify(reversed(s), t)`. A relation that depends on
  which argument came first would make the broad phase's pair ordering
  load-bearing.
- `crossing_point` on a crossing engineered to land exactly on a lattice point:
  the returned `GridPoint` is that point exactly.
- `crossing_point` on a near-parallel pair: the result lies within the
  intersection of the two segments' coordinate ranges. This is the unconditional
  invariant, and it is the one the clamp is for.
- **Collinear overlap, both sides of ruling 4.** A dyadic-spacing diagonal
  overlap and an axis-aligned decimal-spacing overlap assert `Overlapping`. A
  0.1 m diagonal overlap at UTM33 magnitudes asserts only
  `Crossing || Disjoint`, with the comment saying why anything stronger is a
  coin flip. The second is the fixture the provisional draft would not have
  written, and it is the one that stops a later reader "fixing" the arm.

**`segment_meets_cell` — the hot-pixel block, and it carries the audit's one
finding.** In the same file, because it shares the fixtures.

- **The finding itself, as one test.** A grid-collinear diagonal T-junction at
  0.1 m at UTM33 magnitudes: a node on a host segment by construction *in index
  space*. `on_segment<DefaultKernel>` returns **false** — asserted, so the
  difference between the two predicates is pinned rather than described — and
  `segment_meets_cell` returns **true**. If one test from this increment is read
  in five years, it should be this one.
- A segment crossing a cell's interior through no lattice point at all.
- A segment grazing exactly one cell corner: `true`. The strict-sign mutant.
- A segment whose *line* crosses the cell but whose extent stops short of it:
  `false`. The dropped-bbox mutant, and it needs a long segment far from the
  cell, which is the fixture nobody writes by accident.
- A segment inside the cell's bounding box but strictly to one side of its
  line — a long diagonal against a cell near a corner of its bbox: `false`.
- Axis-aligned segments along a cell edge, and along the line one half-cell out.
- Zero-length segments: at the cell centre, just inside a corner, just outside.
- The cell of each of a segment's own endpoints: always `true`, trivially, and
  asserted because step 4 of the split pass relies on it to not lose endpoints.
- Symmetry under reversing the segment, for the same reason `classify` has it.

**`tests/cpp/property/prop_noding_snap_invariants.cpp` — invariant-critical.**
Catch2 `GENERATE` over a seeded range, per `testing.md`'s framework section.

- **Displacement**: `|snapped(p) - p| <= spacing/√2` per point, asserted per axis
  as `<= spacing/2` plus one ulp. This is the numeric content of `testing.md`'s
  `√2·h` bullet, which is arithmetically right and attached to a false claim; see
  the fixes below.
- **Idempotence**: `snap(world(g)) == g` for every `g` within the bound, and
  `snapped(snapped(p)) == snapped(p)` bitwise.
- **Key equivalence**: `snap(p) == snap(q)` iff `snapped(p) == snapped(q)`
  bitwise. Both directions, over generated pairs including pairs deliberately
  placed on either side of a cell boundary. This is the property dedup rests on.
- **Order independence**: shuffling the input to `NodeSet` does not change
  `points()` nor any `id_of`. The reproducibility claim, and the reason ids are
  in grid order.
- **Relation totality**: over generated segment pairs at UTM33 magnitudes,
  `classify<DefaultKernel>` agrees with an independent oracle built from the four
  orientations written out longhand — the same trick increment 4 used to keep the
  mask honest. **The oracle may not call `on_segment<K>`**: `classify`'s endpoint
  arm calls it, so an oracle that shares it is tautological on exactly the
  exact-incidence cases the property is about. The oracle spells betweenness out
  instead, in the four comparisons `segment.hpp:68-74` uses —
  `min(a.x,b.x) <= p.x <= max(a.x,b.x)` and likewise in y, with the collinearity
  it is conditioned on coming from the oracle's own `orient2d` call. Four lines,
  no shared code, and the sharing is then only of `K` itself, which is the
  subject of the test rather than a component of it.
- `Crossing` implies `crossing_point` lies within both coordinate ranges, on
  every generated crossing.
- **Hot-pixel dominance**: over generated (segment, point) pairs,
  `on_segment<DefaultKernel>(s, p)` implies
  `segment_meets_cell<DefaultKernel>(grid, s, snap(p))` — the hot-pixel
  predicate is strictly weaker than exact incidence and never misses one. The
  converse does not hold and the property must not assert it; that gap *is* the
  finding.
- **Hot-pixel containment**: `segment_meets_cell(grid, s, g)` implies
  `world(g)` is within `spacing/√2` of `s`, per axis as `spacing/2` plus one
  ulp. This is the bound the split pass's deformation inherits and it is the
  same arithmetic as the displacement property above.

**`tests/cpp/unit/test_noding_node_set.cpp` — not invariant-critical.** Sort,
unique, `lower_bound` over a totally ordered integral key with a defaulted
`operator<=>`. Its failure mode is a typo or an off-by-one at a boundary, both of
which an exhaustive small-input suite catches; there is no topology decision
inside it, because the decision "equal keys collapse" is made by the key, and the
key is 5a's snapping. Mutation spend here buys nothing, on the same argument
increment 4 made for `test_cdt_indexed_mesh`. It still owes one adversarial
test — two points straddling a cell boundary by one ulp, which must remain two
nodes — because that is the case where a reader expects collapse.

### Mutants the round must kill

1. `std::rint` or `std::nearbyint` in place of `std::llround` — the tie fixture.
2. `snap` truncating (`static_cast<std::int64_t>(x / s)`) instead of rounding.
3. Multiplying by a cached reciprocal instead of dividing. Killed by a seeded
   search over the property generator, not by a hand-picked fixture; the test
   comment must say that is where the input came from.
4. `kMaxGridIndex` raised to `2^62`, or `can_snap` comparing `<` where it needs
   `<=` or vice versa — the four-corner boundary block.
5. `is_valid_spacing` written as `spacing != 0.0` — killed by exactly three
   fixtures, `-1.0`, `-0.1` and `-inf`, which it wrongly accepts. **Not** by
   `-0.0`, which compares equal to `0.0` in IEEE, and **not** by `NaN` or
   `+inf`, which the `<= max()` conjunct rejects under either spelling.
   Measured against the built mutant, not reasoned from the clause.
6. `Touching` and `Crossing` swapped when the meeting point is an endpoint; any
   one of the four orientation signs flipped.
7. The collinear arm returning `Overlapping` for a single shared point.
8. The collinear arm testing one axis instead of both — the vertical pair.
9. The clamp removed from `crossing_point` — killed by **a found near-parallel
   fixture**, and by that alone. This entry previously ruled the opposite: that
   the universal coordinate-range assertion, which holds for *every* crossing,
   killed the mutant without a specially-found pair, and that such a pair "may
   not exist at reasonable search effort". **Both halves of that are false, and
   the suite measured it.** A well-conditioned crossing lands inside the
   coordinate ranges with or without the clamp, so the universal invariant alone
   leaves the mutant alive through the whole suite; and the fixture does exist —
   a seeded 4e5-sample scan over near-parallel crossings at Web Mercator
   magnitudes found one in under a minute (`u = 1.5`, i.e. half a segment beyond
   `s.b`, putting the unclamped point 3.3e7 m outside the range box on a pair
   whose true crossing is interior to both). It is the last case in
   `test_noding_pairwise_intersection.cpp`'s "crossing_point lands in the
   intersection of the coordinate ranges", with its provenance at the fixture.
   The universal invariant stays, because it is the statement of the contract;
   it is simply not what kills this mutant.
10. `NodeSet` ids assigned in first-appearance order — killed by the shuffle
    property, and by nothing else.
11. `NodeSet` keyed on `world()` coordinates rather than on `GridPoint`.
    **Behaviourally invisible under injectivity, so it is not a mutant and is
    not pinned.** Listed so nobody spends the round on it. The reason to key on
    the grid point anyway is that it needs no injectivity assumption and no
    floating-point comparison at all.
12. An input-dependent anchor — subtracting a bounding-box origin before
    snapping. **Not expressible in 5a**: `snap`'s signature has no argument it
    could arrive through. It is a live mutant in 5b's driver and belongs to that
    round.
13. `segment_meets_cell` written with strict signs (`all > 0 || all < 0`)
    instead of "all the same non-`Collinear` value" — the grazed-corner fixture.
    The most innocent-looking diff in 5a.
14. The bbox half of `segment_meets_cell` dropped — the far-away-line fixture.
15. The orientation half dropped, leaving a bbox test — the diagonal-against-a-
    corner fixture.
16. `segment_meets_cell` implemented as `on_segment` — i.e. the provisional
    design. Killed by the T-junction fixture, and by nothing else in the file.
17. Cell corners computed as `world(g) ± spacing/2` — two roundings. Killed by
    the seeded last-bit fixture; note this mutant is *behaviourally invisible*
    at a dyadic spacing, so it must be hunted at 0.1 m.
18. The `on_segment<K>` call in `classify`'s endpoint arm replaced by
    `return Disjoint`. **Killed by nothing in the suite, and that is accepted,
    not an oversight.** Under `DefaultKernel` the arm never returns `Touching`
    (see "the endpoint arm is a totality fallback"), so no fixture can kill it
    there without asserting a coin flip. The single `TEMPLATE_TEST_CASE` *does*
    reach the arm under `FastKernel` — that correction is in "Template spend" —
    but its assertion is "not `Crossing`", which both `Touching` and the
    mutant's `Disjoint` satisfy, so it does not kill it either. The arm is
    defended by the case analysis, not by a test.
    Listed, like 11, so nobody spends the round hunting the fixture.

### Template spend

**One opt-in cross product, one instantiation otherwise.** `classify<K>` and
`segment_meets_cell<K>` are the two kernel-parameterised functions in 5a that
decide anything.

Exactly one `TEMPLATE_TEST_CASE` over `{FastKernel, DefaultKernel}`, doing one
job: a near-parallel crossing at **Web Mercator full-extent** magnitudes where
`FastKernel` gets the orientation wrong and answers something other than
`Crossing`, while `DefaultKernel` reports `Crossing`. What it proves is that the classification
actually flows through `K` — that nobody has written the four determinants
inline in `double` — which is the same job, and the same justification, as
increment 3's single sliver-winding cross product.

**This paragraph previously specified the fixture at UTM33 magnitudes and
predicted a mis-signed orientation reported as `Disjoint`. Both were wrong, and
the suite measured both.** At UTM33 magnitudes the pair does not exist:
`FastKernel`'s error in `orient2d` goes as `eps·|dx·dy|`, while the smallest
nonzero determinant the coordinates can express goes as `|dx|·ulp(y)` ≈
`eps·|dx|·|y|`, so the first exceeds the second only when a point's offset from
the segment is comparable to the absolute coordinate magnitude — at northing
7.9e6 over even a 100 km domain that ratio is 0.013. Measured: over 1.8e6
near-line triples about a 100 km UTM33 segment the two kernels disagreed on
**zero**. At Web Mercator full extent (±2.0037e7 with a 4e7 span) the ratio
reaches 2 and `FastKernel` does go wrong — but it answers `Collinear` rather
than a flipped sign, which drops a genuine crossing into `classify`'s **endpoint
arm**, exactly the case that arm is documented to exist for. The arm's answer
there is `Touching` without floating-point contraction and `Disjoint` with it,
so the assertion is the robust one they share: under `FastKernel` the answer is
not `Crossing`. That is a better outcome than the prediction, because the
endpoint arm now has a fixture behind it rather than only a case analysis — note
that this does *not* make mutant 18 killable, since the mutant's `return
Disjoint` is one of the two answers this fixture already tolerates.

**Everything else runs under `DefaultKernel` only, and `FastKernel` is forbidden
in the property suite.** Same ruling and same reason as increment 4: the property
suite's oracle is built from the kernel, and an oracle less exact than the
subject reports false failures on precisely the degenerate inputs the property is
about. `parallel_refinement.md:41-51` adds a second reason specific to this
module — snapped data is *where* the filter fall-throughs live, 38 % of
grid-collinear triples at 0.1 m, so `FastKernel` is at its worst on exactly this
increment's inputs. What snapping manufactures is fall-through with a definite
sign, not an unresolvable triple; ruling 4 of this document is what corrected
that paragraph, and this citation follows it.

`crossing_point<K>` is instantiated once, under `DefaultKernel`. Its arithmetic
does not depend on `K`.

`segment_meets_cell<K>` is instantiated once, under `DefaultKernel`, and does
**not** get a second cross product. The one it would buy — `FastKernel`
mis-signing a corner — proves nothing `classify`'s cross product has not already
proved about flow-through-`K`, and it would cost a second compile of the most
fixture-heavy file in 5a. The README's rule is that a cross product must prove
something a later increment depends on; this one would not.

## Risks

1. **Snap rounding can introduce crossings that were not in the input.** Two
   segments passing within a cell of each other snap onto one another. This is
   inherent to snap rounding, not a defect here, and it is why guarantee 14 is
   *verified* rather than assumed. Residual after mitigation: input that needs
   more rounds than the cap is rejected with `NotConverged` rather than meshed.
   The operational lever is a finer spacing, which trades this risk against
   risk 2.
2. **Snapping degenerates rings, and the two levers are opposed.** A coarse
   spacing collapses small holes (`RingCollapsed`) and flattens slivers
   (`RingDegenerateAfterSnap`); a fine spacing collapses fewer near-coincidences
   and so leaves more crossings for risk 1 and more work for the CDT. There is no
   value that is safe for all input, which is why there is **no default spacing
   anywhere in C++** — the same ruling, for the same reason, as `K` having no
   default in `ring.hpp` and `pslg_builder.hpp`. It is a declared parameter of
   input precision and boundary tolerance (`parallel_refinement.md:130-137`).
3. **The winding re-check is a consequence nothing in the tree anticipated.**
   `Pslg` guarantee 6 is a statement about pre-snap coordinates and this document
   is the first to say so. The risk is real: a reversed sliver hole that reaches
   `detria` meshes as an island with no diagnostic anywhere — a silent wrong
   mesh, which is the failure class this project spends its mutation budget
   avoiding.

   **The mitigation is not the winding re-check, and this entry said it was.**
   Measured during increment 5b's red round: four candidate slivers built to flip
   a hole's winding were all refused as `NonSimpleRing` before the winding check
   ran, because a winding flip puts a vertex inside the hot pixel of the opposite
   edge, the split pass is then *required* by guarantee 14(b) to split that edge
   at it, and the ring repeats a node id. Drop the winding re-check entirely and
   the ring is still refused. **The mitigation is guarantee 14(b) plus the
   node-repeat `NonSimpleRing` check**; the winding re-check is defence in depth
   behind it, and `NodeStatus::RingDegenerateAfterSnap` is a self-check with a
   name rather than a diagnosis. **The relocation is measured, not inferred**: a
   bounded exhaustive search over 4- and 5-gons finds that *every* ring whose
   winding flips under snapping is caught by 14(b) — in every block, under two
   independently written programs — so deleting the winding re-check refuses all
   of them anyway. The ruling, the mechanism, the
   search with its program and its honest bound, and what would reopen it are in
   `docs/increments/05b-noder-driver.md`, section "`RingDegenerateAfterSnap` is a
   self-check". This correction is increment 5b's, made here because this is
   where the claim lives.
4. **Construction error is unbounded for near-parallel pairs**, and the clamp
   bounds it only to the overlap region. For such a pair the combinatorial answer
   is whatever snap rounding says. Not fixable without exact constructions, which
   `01-predicates.md:290-295` prohibits by name ("signs only"). Named, bounded,
   accepted.
5. **Index identity with the input `Pslg` does not survive noding.** Guarantee 9
   is what lets Python hold a parallel attribute array, and 5b must deliver
   `node_of_input_vertex` or that capability is lost silently. Called out in 5a's
   record because it is a promise made two increments earlier and broken here.
6. **The sequencing risk persists one more increment.** The CDT still takes a
   `const Pslg&` until 5b merges, so increment 3 risk 1 and increment 4 risk 4
   are inherited unchanged by 5a and are not discharged by it. 5a improves
   nothing end to end; it is deliberately a foundation PR.
7. **The narrow phase is quadratic in bucket occupancy.** A dense breakline
   network — a river delta, an urban road layer — can put thousands of segments
   in one bucket. 5b owns the bucket sizing, and the mitigation is a bucket size
   driven by segment density; the residual is that adversarially dense input is
   slow, not wrong.
8. **`SegmentRelation` will attract a fifth enumerator.** "Nearly touching",
   "collinear within tolerance", "shares a bucket". Each one is a tolerance, and
   this project has none (increment 2, "Tolerances: there are none"). The four
   are exhaustive over the exact classification of two closed segments; anything
   further is a different question and needs a different function name. The
   first such question has now arrived and been answered that way:
   `segment_meets_cell` is a separate function, not a relation, not a parameter
   on `classify`. That is the precedent; the risk is that the next one is
   resolved differently because it looks smaller.
9. **Hot-pixel incidence is proximity, so splitting deforms.** A segment split at
   a node whose cell it merely *passes through* acquires a vertex that was not on
   it, moving it by up to `h/√2`. That is not a defect — it is the definition of
   snap rounding, it is bounded, and it is the reason `testing.md`'s length bullet
   has to become a one-way Hausdorff claim (fix 8). The consequence that bites is
   that a deformed segment can enter a cell it previously missed, which is what
   feeds risk 1's iteration. Residual: input where deformation keeps finding new
   cells past the cap is rejected `NotConverged`, not meshed. The measured
   counterweight is that the cascade is rarer than the theory allows — the audit
   found the near-parallel generator produced no surviving crossings at all at
   0.1 m, because two segments closer than a cell snap into coincidence rather
   than into a crossing.
10. **Arc order is unfiltered floating-point, not a sign.** The dot product that
    orders a segment's split points is `double` arithmetic outside `K`, so two
    nodes whose projections differ by less than its rounding error can be ordered
    wrongly, and a non-monotone order gives a chain that doubles back — which
    guarantee 14 does not catch. Bounded and accepted for the reason given under
    "the split pass": the affected pairs are sub-cell apart along the segment,
    which is inside the deformation risk 9 already accepts. The mitigation if it
    ever bites is an exact dot-product sign, which is a kernel change, which is
    why it is a named risk here instead of a design change.
11. **A dense cell attracts many splits.** Every segment through one hot pixel is
    split at it, so a busy junction produces a high-degree node and a burst of
    short edges. This is correct output and it is what the CDT wants; the cost is
    that step 4 of the split pass is quadratic in (segments × nodes) per bucket,
    which is risk 7 with a second multiplicand. 5b owns the bucket sizing for
    both.
12. **The broad phase is the design's one remaining shared-oracle surface.**
    Split and verify both draw candidates from it, so one dropped candidate
    hides itself. It is the one remaining because guarantee 15 nearly added a
    second — re-deriving `edge_is_river()` from the dedup's own provenance
    map — and that is ruled out above for an input-side oracle in the same
    property. Mitigated at the seam by a conservativeness postcondition with a
    brute-force oracle and by an all-pairs verification of guarantee 14 on small
    inputs ("the split pass"); the residual is that every one of these
    mitigations lives in 5b's suite, so 5a ships the rule and 5b is where it is
    enforced. If 5b lands without them, this document's guarantees 14 and 15 are
    worth less than they read.

## Documentation this PR fixes

Per `docs/increments/README.md`, these are fixed in this PR or not recorded.
There is no ledger.

**Six of the eleven did not land, and that is recorded here because the section
asserts otherwise.** `git diff --stat 41054aa~1 e090909 -- testing.md
parallel_refinement.md project_structure.md docs/increments/03-pslg.md
docs/increments/01-predicates.md` prints two files, and checking those two line
by line rather than by filename: **2a, 2b and 11 shipped; 1, 2, 2c, 3, 4, 5, 6, 7
and 8 did not** — nine of eleven. Step 6 of the algorithm still emitted a flat
segment list and step 4 still said "lying on it". They were re-found while
designing 5b and are fixed in **5b's** PR instead; item 9, the `[planned]`
marker, is 5c's, since it describes what CI enforces over a module that ships.
Line numbers quoted below are the pre-fix ones and are stale by design: the text
they pointed at no longer exists. See `docs/increments/05b-noder-driver.md`,
"Documentation this PR fixes".

**`parallel_refinement.md`**

1. Step 6 of the noding algorithm (line 148) — "Output a list of
   `(p0, p1, is_river)` segments plus the deduplicated vertex array". A flat edge
   list cannot drive the CDT, which needs `addOutline`/`addHole` per ring with
   winding (`04-cdt.md:547-556`). Replace with: output a `NodedPslg` — the
   deduplicated node array, the chains with their roles and windings preserved,
   the per-edge `is_river` array and the snap grid.
2. The edge-metadata section says the merge rule is "logical OR" and does not say
   where the bits live. Add: one byte per output edge, index-aligned with the
   flat edge enumeration; not a sparse override set (see 3). Add also that the
   merge is keyed on node-id pairs, not on a collinear-overlap classification —
   ruling 4's measurement is why.
2a. The "manufactures the hard cases" paragraph repeated
   `01-predicates.md`'s claim verbatim. Same defect, same fix (see 11), applied
   where it is quoted; the corrected text is now at lines 41-51. **Fixed in this PR, not
   recorded as debt** — leaving the quotation while fixing the original is
   exactly the ledger failure `docs/increments/README.md` describes.
2b. The dyadic-lever paragraph bundled two properties of a dyadic spacing and
   dismissed both on the cost of the one that needs local coordinates. Unbundled
   (lines 53-74): exact affine scaling
   needs dyadic spacing **alone** and is what the noder's topology depends on;
   exact determinants need local coordinates and stay dismissed. See "Choosing
   the spacing"; the measurement is `kernel-sufficiency-audit.md` §5.3.
2c. Step 4 of the algorithm, "split each input segment at all intersections lying
   on it" — "lying on" is exact incidence and is the wrong test. Replace with
   the hot-pixel form: split at every node whose cell the segment meets. Step 1's
   "pairwise robust intersection" stays as it is; it is where crossings come
   from.

**`docs/increments/03-pslg.md`**

3. Lines 67-69: `is_river` "contributes a sparse override set then, on `NodedPslg`".
   Overturned above — after noding an edge has no single source chain to override.
   Replace with the dense per-edge array, naming this increment.
4. The `NodedPslg` promise at lines 486-493 lists two guarantees (no crossing
   interiors, coordinates on the snap grid). Add the three it is missing and that
   a reader will otherwise assume: no two vertices equal, no zero-length edge,
   and — the one that surprises — **guarantee 9 does not survive**, replaced by
   `node_of_input_vertex`.

**`project_structure.md`**

5. The `noding` section repeats the sparse-override-set claim. Same fix as 3.
6. The directory listing has no `include/terrain/noding/` at all and no
   `core/snap_grid.hpp` or `core/noded_pslg.hpp`. Add them, and mark
   `src/noding/` as not needed rather than *(planned)*: the noder is header-only
   because its driver is a template on `K`.
7. The `noding` section says the broad phase may use the raster's grid. True and
   already caveated, but add the number that makes the independence
   non-negotiable: a snap-spaced broad phase over a 100 km domain is 4e12
   buckets.

**`testing.md`**

8. The `noding` invariant list, **all five bullets** (lines 106-110). Three are
   false, one is true and kept, one is true and moves:
   - "Every input vertex appears in the output vertex set (post-snap)" — false
     after dedup: two input vertices in one cell produce **one** output vertex.
     Replace with: every input vertex's snapped image is an output node, and
     `node_of_input_vertex` maps it there.
   - "Every output vertex is either an input vertex (snapped) or lies on at least
     two input segments" — false for a constructed crossing point, which
     `e412a43` establishes need not lie on *either* segment after snapping.
     Replace with: every output vertex is either the snapped image of an input
     vertex or the snapped image of a constructed crossing of two input
     segments.
   - "Sum of output segment lengths equals sum of input segment lengths, modulo
     snap perturbation. The bound is `√2·h` per segment" — false twice over. A
     collinear overlap merge **deletes** length outright (the road along the
     river is geometrically forgotten, `parallel_refinement.md:176`), and a
     segment split at `m` points accumulates up to `m` displacements, so the
     bound is not per segment. The `√2·h` arithmetic itself is right and worth
     keeping; the claim it is attached to is not. Replace with the one-way
     Hausdorff form: **every output edge derived from input segment `e` lies
     within the closed `h/√2`-neighbourhood of `e`**, which follows from that
     neighbourhood being convex and both endpoints lying in it — and which holds
     per noding round, so state the `k`-round form as `k·h/√2`.
   - "No two output segments intersect in their interior" — this one is true and
     is guarantee 14(a). Keep, and add that it is *verified* by a second pass,
     not established by construction. Add 14(b) beside it, because it is the
     half a reader will assume follows and it does not: **no node's cell meets
     an edge it is not an endpoint of**, in the hot-pixel form, not in the
     exact-incidence form.
   - "`is_river` bit on any output segment is the OR of the bits on the input
     segments that contributed" (line 110) — **true, and it is guarantee 15.**
     Keep it, and add what it is missing: the OR is over contributing *chains*
     and the merge happens at the node-id edge-key dedup, not at a collinear
     overlap classification. Add also that it is verified rather than assumed,
     and verified against the input — not against the dedup's own provenance
     map, which would only restate the dedup — for the same reason 14 is.
   Add the bullets the list lacks: no two output vertices are equal; no output
   edge is zero-length; every output coordinate satisfies
   `grid.snapped(v) == v` bitwise; node ids are independent of input order and
   of thread count.
9. The section keeps its `[planned]` marker through 5a — 5a ships no noder — and
   becomes `[live]` when 5b merges. Say so in the marker rather than leaving the
   transition to be noticed.
10. The test-layout convention block (lines 250-258) already names
    `tests/cpp/unit/test_noding_snap_rounding.cpp` and
    `test_noding_pairwise_intersection.cpp` as examples. They stop being examples
    in this PR; no change is needed **to those two lines**. Two neighbouring
    lines do need one, and an earlier draft of this item missed them by reading
    only the unit block:
    - Line 257 names `prop_noding_no_crossings.cpp`. 5a ships
      `prop_noding_snap_invariants.cpp` instead, and `prop_noding_no_crossings`
      becomes real at 5b as guarantee 14's brute-force check (see "the split
      pass"). Add the 5a file beside it and mark the other as 5b's, so the
      example does not read as a file that exists.
    - Line 177 says `tests/cpp/property/noding_generators.h` "produces random
      sets of polylines with controllable density of intersections". That is 5b's
      generator; 5a's property suite generates points, spacings and segment
      pairs, not polyline sets. Say which increment the header arrives in rather
      than leaving a reader to look for it in this PR.

**`docs/increments/01-predicates.md`**

11. Lines 149-152: "collinear-but-unresolvable triples are the common case in
    snapped breakline data". Measured false in its most important half —
    `world(g) = g·spacing` is not affine at a non-dyadic spacing, so only 207 of
    2976 general-direction grid-collinear triples report `Collinear` at 0.1 m,
    while 38 % of triples fall through the filter to a **definite** sign. An
    exactly collinear triple is the *cheapest* input measured, not the most
    expensive, because the filter short-circuits it. Replace with the definite-
    sign description, keep the surrounding cost ruling — which is unaffected —
    and cite `kernel-sufficiency-audit.md` §2 and §5.1. **This is a defect found
    during this increment, so it is fixed in this PR and there is no ledger
    entry.** Checked: the claim lives in exactly two files, this one and
    `parallel_refinement.md` (2a); `grep -rn "unresolvable\|breakline"
    include/terrain/predicates/` is empty, so the prose in `kernel.hpp` that the
    audit lists as a third site does not in fact carry it.
