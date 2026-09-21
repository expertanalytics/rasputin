# Increment 5d — the corner graze

Status: design, written after increment 5b's red round found the defect it
fixes. 5d is a **predicate change**, small in lines and not small in what it
touches: `include/terrain/noding/intersect.hpp` is a shipped header with its own
mutation round behind it.

`docs/increments/05b-noder-driver.md` risk 16 is the finding and the ruling that
produced this increment; it is not restated here beyond the one paragraph the
next section needs.

**Clipping to the catchment polygon is parked and is unrelated to this
increment.** It was raised alongside and deliberately deferred by the user; it
is named here only so that nobody later reads 5d's scope as covering it.

## The defect, in one paragraph

The split pass splits an edge at **every** node whose closed cell it meets, and
guarantee 14(b) demands exactly that. On a segment passing exactly through a
lattice corner, the two flanking cells are met **at that single point**, so both
become split points; the resolved chain's middle edge then passes through the
same corner and grazes the two cells it did not split at, which 14(b) refuses.
Round 2 splits it there, and the two configurations map into each other. The
chain grows by two positions per round and never reaches a fixpoint. `05b`'s
section "The tie does not converge, and the probe that shows it" has the orbit.

**The split rule and 14(b) are inconsistent on a single-point corner touch.**
That is the defect. It is not a bug in either one.

## The motivating case: ordinary GIS input, nothing constructed

An earlier account of this defect — mine, in `05b` — used a unit-square
configuration and described the trigger as needing a **constructed** point, which
made it read as a rounding artifact of the noder's own intersection points. It is
not one, and the distinction matters because it decides whether preprocessing the
input to insert intersection nodes could dodge the problem. It could not.

A 45° road and two rivers, at `SnapGrid{0.1}`, **every coordinate a round
decimetre**:

- road, a breakline from `(0.0, 0.0)` to `(1.0, 1.0)`;
- river A, ending at `(0.4, 0.5)`, running above the road and never meeting it;
- river B, ending at `(0.5, 0.4)`, running below the road and never meeting it.

```
A. decimetre GIS input at SnapGrid{0.1}
  road vs riverA = Disjoint | road vs riverB = Disjoint  <-- nothing constructed
  road (0,0)-(1,1)                   vs cell of (0.40,0.50) = 1
  road (0,0)-(1,1)                   vs cell of (0.50,0.40) = 1
  middle edge (0.4,0.5)-(0.5,0.4)    vs cell of (0.40,0.40) = 1
  middle edge (0.4,0.5)-(0.5,0.4)    vs cell of (0.50,0.50) = 1
B. minimal reproduction at SnapGrid{1.0}
  diag (0,0)-(1,1)                   vs cell of (1.00,0.00) = 1
  diag (0,0)-(1,1)                   vs cell of (0.00,1.00) = 1
  middle edge (0,1)-(1,0)            vs cell of (0.00,0.00) = 1
  middle edge (0,1)-(1,0)            vs cell of (1.00,1.00) = 1
```

The first line is the part that settles it: `classify<K>` reports **`Disjoint`**
for the road against each river, so the noder constructs **nothing** on this
input — there is no crossing, no intersection point, no `crossing_point<K>` call.
The road is nevertheless split at `(4, 5)` and `(5, 4)`, because both rivers'
end vertices snap to cells whose corner `(0.45, 0.45)` lies exactly on the road.
The resulting middle edge then grazes the cells of `(4, 4)` and `(5, 5)`.

**The trigger is where input vertices sit relative to the lattice, not a point
the noder had to invent.** Round decimetres on a decimetre grid is not an exotic
input; it is what a decimetre-precision cadastral or hydrographic extract looks
like. Block B above is the same configuration reduced to the unit grid, and it is
kept as the minimal reproduction — the smallest input that reaches the defect,
and the one to reach for when debugging.

Both blocks are from one program, re-run for this document:

```cpp
// /tmp/probe5d.cpp -- compile and run:
//   c++ -std=c++20 -Iinclude /tmp/probe5d.cpp build/libterrain_predicates.a -o /tmp/probe5d
#include <cstdio>
#include "terrain/core/snap_grid.hpp"
#include "terrain/noding/intersect.hpp"
#include "terrain/predicates/default_kernel.hpp"
using namespace terrain;
using K = pred::DefaultKernel;
using noding::SegmentRelation;

static const char* rel(SegmentRelation r) {
    switch (r) {
        case SegmentRelation::Disjoint:    return "Disjoint";
        case SegmentRelation::Touching:    return "Touching";
        case SegmentRelation::Crossing:    return "Crossing";
        case SegmentRelation::Overlapping: return "Overlapping";
    }
    return "?";
}
static void meets(const SnapGrid& g, const char* what, Segment2 s, Point2 v) {
    std::printf("  %-34s vs cell of (%.2f,%.2f) = %d\n", what, v.x, v.y,
                (int)noding::segment_meets_cell<K>(g, s, g.snap(v)));
}

int main() {
    // ---- A. decimetre GIS input, SnapGrid{0.1} -------------------------------
    const SnapGrid g{0.1};
    const Segment2 road{{0.0, 0.0}, {1.0, 1.0}};          // 45 deg, y = x
    const Segment2 riverA{{0.1, 0.5}, {0.4, 0.5}};        // ends at (0.4,0.5), stays above y=x
    const Segment2 riverB{{0.5, 0.1}, {0.5, 0.4}};        // ends at (0.5,0.4), stays below y=x
    std::printf("A. decimetre GIS input at SnapGrid{0.1}\n");
    std::printf("  road vs riverA = %s | road vs riverB = %s  <-- nothing constructed\n",
                rel(noding::classify<K>(road, riverA)), rel(noding::classify<K>(road, riverB)));
    meets(g, "road (0,0)-(1,1)", road, Point2{0.4, 0.5});
    meets(g, "road (0,0)-(1,1)", road, Point2{0.5, 0.4});
    const Segment2 mid{g.world(g.snap(Point2{0.4, 0.5})), g.world(g.snap(Point2{0.5, 0.4}))};
    meets(g, "middle edge (0.4,0.5)-(0.5,0.4)", mid, Point2{0.4, 0.4});
    meets(g, "middle edge (0.4,0.5)-(0.5,0.4)", mid, Point2{0.5, 0.5});

    // ---- B. minimal reproduction, SnapGrid{1.0} ------------------------------
    const SnapGrid u{1.0};
    std::printf("B. minimal reproduction at SnapGrid{1.0}\n");
    const Segment2 diag{u.world({0, 0}), u.world({1, 1})};
    meets(u, "diag (0,0)-(1,1)", diag, u.world({1, 0}));
    meets(u, "diag (0,0)-(1,1)", diag, u.world({0, 1}));
    const Segment2 anti{u.world({0, 1}), u.world({1, 0})};
    meets(u, "middle edge (0,1)-(1,0)", anti, u.world({0, 0}));
    meets(u, "middle edge (0,1)-(1,0)", anti, u.world({1, 1}));
    return 0;
}
```

## The fix

**`segment_meets_cell` is replaced, in both the split rule and guarantee 14(b),
by a predicate that excludes a contact at exactly one point.** Both, not one:
changing the split rule alone leaves 14(b) demanding a split the pass no longer
performs, and changing 14(b) alone leaves the pass manufacturing nodes nothing
requires. The inconsistency is between the two, so the change is to the thing
they share.

**Not half-open cells.** A half-open cell is direction-dependent — it resolves
the tie by preferring one side of the lattice — and it reopens the
T-junction-detection question the kernel audit closed. Excluding a single-point
contact is symmetric, is a measure-zero change to the accepted set, and leaves
every genuine T-junction detected: a node whose cell a segment genuinely passes
through meets it in a sub-segment of **positive length**.

### The predicate, characterised

> The segment meets the **open** cell, **or** it meets the **closed** cell in a
> set of **positive length**.

Two disjuncts, and both are load-bearing. The first is what keeps a **zero-length
segment at a cell centre** true — a degenerate chain must still meet its own
cell. The second is what keeps a segment running **along a cell edge** true — it
never enters the open cell, and it must not be excluded. The first formulation
this design tried was "positive length" alone, and measuring it against the
existing fixtures showed it turning the zero-length-at-centre case false, which
would have been a silent regression in the endpoint invariant. **Only a
single-point contact with the closed cell is excluded, and that is the whole
change.**

This is **not** a sign flip in the existing four-corner test. That test answers
"does the segment's line separate the cell", and a single-point contact is not a
question about the line alone — the contact point must also lie within the
segment's own extent. `@developer` is owed the freedom to pick the formulation;
the observable contract is the table below, and the slab-clip form used to
produce it is one correct implementation, not a mandate.

## What it costs: the shipped header, measured

`include/terrain/noding/intersect.hpp` is 5a's, shipped, with a mutation round
and mutants 13–17 aimed at this function. Every currently-passing assertion over
it was run against the new predicate. **Three assertions change, all in
`tests/cpp/unit/test_noding_pairwise_intersection.cpp`, and no others:**

```
test_noding_pairwise_intersection.cpp, SECTION by SECTION:
[crossing the cell's interior]
  seg(-2,0.2,2,0.3)                                    asserted true  | shipped true  | 5d true  
  seg(-0.4,-0.4,0.4,0.45)                              asserted true  | shipped true  | 5d true  
[grazing exactly one corner]   <-- mutant 13's fixture
  seg(0,1,1,0)                                         asserted true  | shipped true  | 5d false   <== CHANGES
  seg(-1,0,0,-1)                                       asserted true  | shipped true  | 5d false   <== CHANGES
[far-away segment]  (mutant 14)
  seg(100,100,200,200)                                 asserted false | shipped false | 5d false 
[long diagonal, line misses]  (mutant 15)
  seg(-10,-8,10,12)                                    asserted false | shipped false | 5d false 
[axis-aligned along a cell edge]
  seg(-5,0.5,5,0.5) along top edge                     asserted true  | shipped true  | 5d true  
  seg(-5,-0.5,5,-0.5) along bottom edge                asserted true  | shipped true  | 5d true  
  seg(0.5,-5,0.5,5) along right edge                   asserted true  | shipped true  | 5d true  
  seg(-5,1.5,5,1.5) one half-cell out                  asserted false | shipped false | 5d false 
[zero-length segments]
  seg(0,0,0,0) cell centre                             asserted true  | shipped true  | 5d true  
  seg(0.4999,...) just inside a corner                 asserted true  | shipped true  | 5d true  
  seg(0.5,0.5,...) exactly on a corner                 asserted true  | shipped true  | 5d false   <== CHANGES
  seg(0.5001,0.5,...) just outside                     asserted false | shipped false | 5d false 

Dominance (prop_noding_snap_invariants): on_segment => meets_cell
  host through world(1,1), p = world(1,1) [a CENTRE]   asserted true  | shipped true  | 5d true  
  snap(0.5,0.5) = (1,1); segment ENDS there and leaves:
  on_segment TRUE at its own endpoint                  asserted true  | shipped true  | 5d false   <== CHANGES
```

- **`SECTION("grazing exactly one corner")` is mutant 13's fixture, and it
  inverts.** Both assertions become `false`. This is the one to be careful of in
  review: a reader meeting it cold will read the diff as a regression, because
  the section's comment currently explains at length why `true` is the right
  answer and calls the alternative "the most innocent-looking diff in 5a". The
  comment is not wrong about 5a; it is describing the behaviour 5d removes on
  purpose. **Rename the section** rather than editing the assertions under the
  old name, so the diff cannot be misread.
- **`SECTION("zero-length segments")` loses one assertion**, the point *exactly*
  on a corner. The cell-centre and just-inside-a-corner cases are unchanged,
  which is the check that the first disjunct is doing its job.
- **Mutants 14 and 15 keep their fixtures unchanged** — the far-away segment and
  the long diagonal whose line misses are `false` under both predicates.

### The second file, which is 5b's and is red today

`tests/cpp/unit/test_noding_node.cpp`'s first lattice-corner case
(`"the unit diagonal grazes both flanking cells and the tie is exact"`) encodes
the old single-point behaviour in **four** `REQUIRE`s: the diagonal against the
cells of `(1, 0)` and `(0, 1)`, and the anti-diagonal against `(0, 0)` and
`(1, 1)`. Both pairs are grazes through the *same* lattice corner, so **all four
invert together** — there is no partial change to make here.

It is worth being exact about its state, because "currently passing" is the wrong
description and would send a reader to run it. The case touches nothing from 5b —
only `SnapGrid`, `Segment2`, `segment_meets_cell<DefaultKernel>` and
`GridPoint::operator<` — but the **file** also exercises `node<K>`, and
`include/terrain/noding/` currently holds only `intersect.hpp` and
`node_set.hpp`. The suite is registered
(`add_terrain_backend_test(test_noding_node unit/test_noding_node.cpp)`) and does
not build. **It is red now, green when 5b lands, and passing by the time 5d
starts** — which is the state that matters, since 5d inverts it. The rule for
checking this rather than the value: `ls include/terrain/noding/` against the
header list in `05b`'s "Files and LOC".

So the full list of assertions 5d inverts is **three in 5a's shipped suite and
four in 5b's**, and the 5b four are in a case deliberately written to depend on
nothing but 5a, precisely so that this increment could flip it without rewriting
it.

### Does mutant 13 survive?

**Its current fixture stops killing it, and there is a candidate replacement
already in the file**: `SECTION("axis-aligned along a cell edge and one half-cell
out")`, which stays `true` under 5d and which a strict-sign-test mutant gets
wrong, because along a cell edge two corners are collinear rather than one.

That is a candidate, not a result. **5d's round must confirm by construction
that each of mutants 13–17 still has a killing fixture under the new predicate**,
and must not assume it from this paragraph — the mutants were written against a
predicate that no longer exists, and a mutant re-pointed at a new fixture by
argument is exactly the claim `.claude/REQUIRED-READING.md` says to run rather
than reason about. This is the single largest piece of work in 5d and it is
`@tester`'s.

### One invariant narrows, and no assertion of it changes

`tests/cpp/property/prop_noding_snap_invariants.cpp` states hot-pixel dominance:
`on_segment<K>(s, p)` implies `segment_meets_cell(grid, s, snap(p))`. **Under 5d
that is false in general and true for every `p` the suite or the noder uses.**

The counterexample is in the table: `llround` rounds half away from zero, so a
point exactly on a cell boundary snaps to the cell **farther** from the origin,
and `(0.5, 0.5)` is therefore the *minimum corner of its own cell* `(1, 1)`. A
segment ending there and leaving touches that cell at one point. `on_segment` is
true — it is an endpoint — and 5d's predicate is false.

It cannot arise where it would matter, and the reason is structural rather than
lucky: **in the split pass every segment is built from `world(g)`, cell centres,
never from raw input coordinates**, and a centre is interior to its own cell by
half a spacing. The suite probes only `p = grid.world(g)` for the same reason.
So the **comment narrows** — dominance holds for `p = world(g)` — and **no
assertion in that file changes**. Narrowing it is not optional: an unqualified
invariant that is false is worse than a qualified one that is true, and the next
reader will believe it.

## Effect on increment 5b, which must not be lost across the seam

`tests/cpp/unit/test_noding_node.cpp` carries two lattice-corner cases, written
during 5b's red round with this increment already in view:

1. **Case 1 asserts 5a's arithmetic alone** — `segment_meets_cell<K>` on the
   diagonal and the anti-diagonal, calling nothing in 5b. **It becomes 5d's
   regression test with all four of its expectations flipped, without being
   rewritten.** That independence was designed in and is the reason the case
   exists separately; see "The second file, which is 5b's and is red today".
2. **Case 2 asserts `NotConverged`** on the full driver. **It is 5d's to change**,
   to `Ok` with the L-shaped resolved chain.

**Case 2 is mutant 5's only *loud* killer, and dropping it rather than changing
it costs that.**
Under the wrong loop condition — "until a round produces no new node" — the tie
converges to `Ok` at the end of round 1, because the flanking nodes already
exist and round 1 manufactures none. It is the only fixture on which the two loop
spellings give different *statuses* rather than different edge counts. Mutant 5
has another killer — the road-along-a-river fixture in the same file — so the
claim is not that coverage vanishes; it is that the only failure a reader cannot
misread does. `05b`'s test-plan section says this in the same narrow form; it is
repeated here because it is precisely the kind of thing a later increment deletes
in good faith.

`node.hpp` and `noded_pslg_builder.hpp` each call the predicate by name and keep
calling it by name, so no production line in 5b moves. **One 5b design claim does
move, and it is the one below.**

### `RingDegenerateAfterSnap`'s demotion rests on this predicate, and 5d must re-run the search that established it

`05b-noder-driver.md` demotes `NodeStatus::RingDegenerateAfterSnap` to a
self-check with no fixture owed. The argument is not about windings at all: a
ring whose snapped winding flipped has a vertex within half a cell of an edge it
is not an endpoint of, so **14(b) requires that edge to be split there**, the
split repeats a node id, and `NonSimpleRing` fires one stage earlier. "Within half
a cell of an edge" is `segment_meets_cell` — **the predicate this increment
narrows.** Narrow it and a ring that is split-and-repeated today may stop being
split, at which point the ring reaches the winding check intact and the demoted
status is reachable again.

So 5d carries two obligations here, and the second is the one easily missed:

1. **Re-run the orientation probe.** `05b-noder-driver.md`'s residual: over 500 000
   split-shaped rings `orientation<K>` never returned `Collinear`, but its reading
   matched the pre-split ring's in only 496 590 — the other 3 410 would report
   `RingDegenerateAfterSnap` where the honest diagnosis is `NonSimpleRing`. One
   line in the ring search reports it.
2. **Re-run the bounded ring search itself**, the whole program, against the *new*
   predicate — not only the new probe. `05b`'s `SURVIVORS 0` and `FLIPPED ==
   grazed` in every block are measurements of the **old** `segment_meets_cell`,
   and they are what the demotion actually rests on. If a narrowed predicate makes
   `SURVIVORS` nonzero, the demotion is refuted and `RingDegenerateAfterSnap`
   becomes a diagnosis owed a fixture and a lever.

Obligation 2 is `.claude/REQUIRED-READING.md`'s rule about deriving the probe set
from the code **as fixed** rather than from the bug as found, in its exact shape:
this increment *narrows* what the predicate accepts, which narrows what gets
split, which widens what can reach the winding check. Re-running only the probe
that 5b wrote against the predicate 5d ships tests the old question with the new
code — the same move as the five UTF-8 inputs in `fdbd532` that verified a fix
which had just widened the glob. The program is pasted inline in `05b`'s "The
bounded ring search" section for precisely this: a design that names a throwaway
it does not contain is not re-runnable by the reader who needs it.

**And the orientation probe belongs here rather than at 5b**, which `@tester` and
`@orchestrator` both concluded independently and this design adopts. Its subject
is a predicate with a scheduled expiry — this increment — and at 5b the 3 410 in
500 000 changes *which name a refusal carries*, not whether the input is refused:
either way the ring is rejected and the caller is told to repair the polygon.
Nobody is owed that distinction at 5b. Here they are, because here the predicate
that produces it changes.

## Degeneracy policy

Unchanged from `05-noder.md` and `05b`, with one row's meaning altered:
`NodeStatus::NotConverged` stops being reachable from a corner graze. It stays
reachable from `05-noder.md` risk 1 — snap rounding creating a genuine crossing
that needs more rounds than the cap — so the status is not demoted and its
"finer spacing" lever becomes unconditionally honest. **`05b`'s status table
carries an exception note for `NotConverged` that 5d deletes**, and that deletion
is part of this increment rather than a later tidy-up.

## Prior art in `legacy/`

**Nothing, and one name collision worth knowing about.**

```sh
grep -rln "graze\|grazing\|corner\|hot pixel\|hot_pixel\|snap_round\|meets_cell" legacy/
```

returns **one** file: `legacy/rasputin/triangulate_dem.h`. The hit is at `:380`,
`// Determine the cell corners`, inside `get_interpolated_value_at_point` — the
corners of a **raster cell**, for bilinear interpolation of a DTM. A different
cell, a different grid and a different question from a snap grid's hot pixel.
The word collides; the concept does not. Anyone grepping the legacy tree for
prior art on this increment will land there first, which is why it is recorded
rather than filtered out.

```sh
grep -rln "noding\|snap\|intersect\|split" legacy/
```

returns the same six files `05b` enumerates — `legacy/bindings.cpp`,
`legacy/rasputin/triangulate_dem.h`, `legacy/rasputin/reader.py`,
`legacy/rasputin/geometry.py`, `legacy/rasputin/gml_repository.py`,
`legacy/tests/test_polygons.py` — and `05b`'s "Prior art" section classifies
every one of them: polygon booleans delegated to CGAL or Shapely, ray–mesh
intersection for shadows, and string splitting. None is a noder. That analysis is
not repeated here.

```sh
grep -rni "predicate\|orient2d\|orientation" legacy/
```

returns CGAL policy headers, `Exact_predicates_inexact_constructions_kernel`, the
`Exact_predicates_tag` on the CDT, and two comments about Shapely-versus-CGAL
ring orientation. The legacy tree has **no geometric predicate of its own** — it
names CGAL's. There is nothing to carry across, and `@migration-expert` is not
spawned for 5d, for the same reason as 5b and 5c.

## Files and LOC

C++ instrument, `grep -vcE '^\s*(//|$)'`, per `05b`'s "Files and LOC".

| File | Contents | Est. LOC |
|---|---|---|
| `include/terrain/noding/intersect.hpp` | the replacement predicate | ~45 |
| `include/terrain/noding/node.hpp`, `include/terrain/noding/noded_pslg_builder.hpp` | call sites, by name | ~4 |

**~50 production lines.** Small, comfortably one PR, and the ceiling argument in
`CLAUDE.md` §2 does not need re-deriving for it. The round is **not** cheap in
proportion: the work is `@tester`'s re-confirmation that mutants 13–17 still die,
over a shipped invariant-critical suite, and that is where the review time goes.

## Risks

1. **A mutant re-pointed at a new fixture by argument rather than by
   construction.** Mutants 13–17 were written against a predicate that 5d
   removes. The mitigation is the ruling above — each must be re-confirmed by
   construction in 5d's round — and the residual is that a mutant quietly loses
   coverage while appearing to keep it.
2. **The inverted fixture reads as a regression.** `SECTION("grazing exactly one
   corner")` currently carries a comment arguing at length for `true`. Mitigated
   by renaming the section rather than editing assertions under the old name, so
   the diff shows an intentional replacement. Residual: a reviewer reading only
   `05-noder.md` or 5a's suite arrives believing `true` is required.
3. **Hot-pixel dominance is narrowed, not preserved.** Mitigated by the
   structural argument that the split pass builds segments from cell centres, and
   by the comment narrowing. Residual: a future caller that passes a raw input
   coordinate rather than `world(g)` re-enters the gap, and nothing in the type
   system stops it.
4. **5d has no production caller until 5c has landed**, inheriting `05b` risk 18
   one increment further. Whether 5d lands before or after 5c is a sequencing
   choice with no correctness content, and this design does not make it.
