# Audit: is the predicates kernel sufficient under the corrected snap doctrine?

Status: **audit, not a design.** Nothing here proposes an increment. Written
against increment 1 as shipped, `parallel_refinement.md` as it stands after
`e412a43`/#67, the vendored `lib/detria/detria.hpp` at the pinned SHA, and the
untracked provisional `docs/increments/05-noder.md`.

The question: increment 1 deferred the interval stage and excluded exact
constructions on a cost model that assumed the filter rarely falls through.
`01-predicates.md:151` already contradicted that assumption and `e412a43`
made the contradiction doctrine — snapping *manufactures* the hard cases. Does
the kernel survive the correction?

## Verdict

**No. The kernel is sufficient, all four deferrals under audit survive, and
none of them survives on the discredited cost model — each is re-justified
below on a measurement.** The single real finding is not a kernel finding: it
is that increment 5a's `classify` asks an exact-incidence question where snap
rounding needs a hot-pixel proximity question, and the fix is four `orient2d`
calls the kernel already provides. Full verdict table at the end.

## Sections

1. Sign-only versus the noder's need for points — *the one that was flagged as
   most worrying, and the one the measurement clears most cleanly*
2. The incircle double-fallback cost
3. `DetriaExact` at the largest coordinates a projected CRS admits
4. Two stages or three
5. What else the corrected snap doctrine invalidates
6. Appendix: the probes, in full, so every table can be re-run

Each section states its conclusion first and its evidence after. Sections were
written to disk as they were settled rather than composed in context, because
the previous attempt at this audit lost everything at the session limit.

---

## 1. Sign-only against the noder's need for points

**Settled: sign-only is sufficient. No exact construction is needed, and the
"construct in double, then snap" pathology is measurably not the one to worry
about. There *is* a gap, it is real, and it lives in increment 5a's algorithm
rather than in the kernel — and it closes with `orient2d` alone.**

### The classical pathology, measured

The worry is that a computed intersection snaps into a cell whose true
intersection lies elsewhere. Measured directly: generate crossing segment pairs
with grid-snapped endpoints, compute the intersection by `05-noder.md`'s exact
formula (`u = cross(t.a - s.a, d2)/den`, clamp, `llround(p/s)`) in double, and
compare the resulting `GridPoint` against the snap of the **exact rational**
intersection computed in `fractions.Fraction`:

```
family                        crossings  non-tie cell mismatch  exact-tie  clamp fired
UTM33  s=0.1  generic             20000                      0       6955            0
WebMerc s=0.1 generic             20000                      0       7062            0
WebMerc s=0.01 near-parallel       5250                      0       1260            0
WebMerc s=2^-4 near-parallel        239                      0          0            0
```

Read that carefully, because the headline number is the one that is *not* a
defect:

- **Non-tie mismatches: zero in ~45,000 crossings.** Whenever the true
  intersection is not sitting exactly on a half-cell boundary, the double
  construction picks the correct cell. The pathology does not show up at any
  spacing or magnitude tried.
- **The ~35 % "exact-tie" column is an artifact of the pipeline's own design,
  not of arithmetic error.** Both segments have grid endpoints, so a large
  fraction of intersections land on an exact half-cell boundary and the tie-break
  decides. `llround` and exact half-away-from-zero rounding of `p/s` then
  disagree by one cell, because `s = 0.1` is not the exact rational 1/10 and
  `p/s` in double lands on the other side. Both answers are legal snap
  roundings; `llround` is deterministic and a pure function of `p`, which is all
  `05-noder.md`'s contract asks for. Exact construction would *change* these
  answers, not correct them.
- **The clamp never fired**, in any family, including the deliberately
  near-parallel ones. The near-parallel regime is also partly self-limiting:
  at `s = 0.1` the near-parallel generator produced *no* surviving crossings at
  all, because two segments closer than a cell snap into coincidence rather than
  into a crossing. `05-noder.md`'s risk 4 is therefore correctly described as
  bounded and accepted; nothing here says it needs exact constructions.

So the first half of the classical worry is closed by measurement. The second
half — snapping creating fresh near-intersections — is real, is what Halperin &
Packer's iterated variant addresses, and is already handled architecturally:
`05-noder.md` guarantee 14 is established by **verification** (a second
broad-phase pass asserting `classify` returns `Disjoint` or `Touching`), with
iteration to a fixpoint and a `NotConverged` status at a cap. Because the
invariant is *checked* rather than assumed, an inexact construction degrades
mesh quality and iteration count; it cannot produce a `NodedPslg` that lies.
That is the structural reason sign-only survives, and it is worth stating as the
reason rather than leaning on the measurement alone.

### The gap that is real

`on_segment<K>` (`include/terrain/core/segment.hpp:68-74`) returns `false`
unless `K::orient2d(s.a, s.b, p) == Collinear` — **exact** incidence. Every
branch of `05-noder.md`'s `classify` that could detect a T-junction or an
overlap goes through that test or through "all four orientations `Collinear`".

Snap rounding does not ask an exact-incidence question. Goodrich-Guibas-
Hershberger-Tanenbaum route each segment through every *hot pixel it passes
through*, not through every pixel whose centre it exactly contains. Those are
different predicates, and §5 shows they differ on the great majority of snapped
input: at `s = 0.1` only 207 of 2976 grid-collinear general-direction triples
come back `Collinear`.

The consequence for `parallel_refinement.md:200`'s watch-list item — "T-junctions
where a polyline endpoint lands on another segment's interior" — is that
`classify` will report `Disjoint` for most of them, the host segment is never
split, and the T-junction survives into the CDT. The verification pass of
guarantee 14 will agree it is `Disjoint`, so nothing catches it. This is not a
wrong mesh; it is a mesh in which a breakline junction silently fails to be a
junction.

### The minimum addition, and it is not in this module

> Does the closed cell of `GridPoint g` — the square
> `[(ix±½)s] × [(iy±½)s]` — intersect segment `s`?

Four `orient2d` calls against the cell's corners (not all the same sign) plus an
interval overlap of the two bounding boxes. Exact, total, division-free,
tolerance-free, **and expressible in the kernel exactly as it stands**. It needs
no exact construction, no new backend entry point and no change to
`ExactPredicates` or `GeometryKernel`.

**That is the answer to the question as posed.** The kernel is sufficient; what
increment 5a is missing is a hot-pixel predicate in `noding/`, built from
`orient2d`. Nothing in `include/terrain/predicates/` needs to change for the
noder to be correct.

One caveat to record for whoever writes it: the cell corners `(ix ± ½) * s` are
themselves constructed, and at a decimal spacing they are constructed
inexactly — which pushes in the same direction as §5's finding about the snap
spacing.

---

## 2. The incircle double-fallback: measured, and it is not a cost problem

**Settled: no. The escape hatch should stay deferred, and the reason is
stronger than "the cost is acceptable" — there is nothing to profile.**

### There is no production caller

`grep -rn incircle include src` returns the predicates module, one *comment* in
`include/terrain/cdt/triangulate.hpp` — obligation 5 of the `CdtBackend`
contract, which `grep -n incircle` over that file is what finds — and nothing
else. The CDT is detria's own triangulator behind increment 4's wrapper;
`04-cdt.md`'s "The property suite runs under `DefaultKernel` only" ruling uses
our `incircle` as a **test oracle**, not in the mesh path. `DefaultKernel::incircle`
is today called only by tests.

An interface addition justified by "profiling demands it"
(`01-predicates.md:179`) cannot be brought forward before there is a profile,
and there cannot be a profile before the first caller that maintains a known
winding — the Lawson flip pass, step 6 of `parallel_refinement.md`. That is the
increment in which the question becomes answerable.

### What the cost actually is, measured

`-O2 -DNDEBUG`, Apple clang, 20k queries × 20 reps per family, ns/call,
linked against `build/libterrain_predicates.a`:

```
family                             orient F orient K  incir F  incir K
well separated                          1.1      2.7      4.0      7.3
grid-collinear, general dir             1.3     11.5     12.0     20.9
grid-collinear, axis aligned            1.1      3.2      1.6      4.3
4 grid-collinear pts (breakline)        1.0      9.4     10.2     27.7
```

`K` is `DefaultKernel`, `F` is `FastKernel`. Grid spacing 0.1 m at UTM33
magnitude; the families are snapped lattice geometry, i.e. exactly the data the
noder emits.

Three things fall out, and two of them correct the increment's prose:

1. **The genuine double-fallback case is the last row** — four grid-collinear
   points, which is a breakline with a fourth point on it, the most ordinary
   input the refiner will ever see. Both the orientation filter *and* the
   incircle filter fall through, and the backend runs twice. It costs
   **27.7 ns against 7.3 ns**, a factor of 3.8 and about 20 ns absolute. That is
   a real slowdown and it is not a cliff. For comparison, one cache miss is
   ~80 ns.
2. **The escape hatch would recover about a third of it.** The orientation half
   of the last row is 9.4 ns; an `incircle_ccw` on `GeometryKernel` for a caller
   holding a known winding removes that and leaves ~18 ns. A 35 % saving on the
   worst family and **zero** on the common path.
3. **`01-predicates.md:151` overstates the worst case in one direction and
   understates the structure.** "Collinear-but-unresolvable triples are the
   common case … that is the expensive path" is only half right: an *exactly*
   collinear triple is the **cheapest** input in the table (row 3, 4.3 ns —
   faster than well-separated). `FilteredKernel::incircle` returns `Cocircular`
   without evaluating a lifted determinant at all, and detria's own filter
   settles a determinant that is exactly zero without entering
   `orient2dadapt`. The expensive family is *near*-collinear with a definite
   sign, which — see §5 — is what snapped data at a decimal spacing actually
   produces, and is not what that sentence describes.

---

## 3. `DetriaExact` at the largest coordinates a projected CRS admits

**Settled: the bounds hold, and they hold for a structural reason, not by luck.
No magnitude admitted by any projected CRS threatens them.** Checked to 1e9 m
and to 1e15 empirically, and read out of the vendored source rather than assumed
from Shewchuk's paper.

### What the source actually says

`lib/detria/detria.hpp:136-178` is a line-for-line port of Shewchuk's
`exactinit`, computing the same constants at `constexpr` time:
`ccwerrboundA = (3 + 16ε)ε`, `iccerrboundA = (10 + 96ε)ε` — the two our
`exact.hpp:54-55` transcribes — plus the B and C bounds that the *adaptive*
stages use internally (`orient2dadapt` at line 391, `incircleadapt` at 521).

Two properties matter for the question and both are visible in that code:

1. **Every bound is relative to a permanent formed from the same coordinates.**
   There is no absolute term anywhere in `ErrorBounds`. A uniform scaling of all
   input coordinates scales determinant and permanent by the same power, so the
   filter's decision is scale-invariant to within the exponent arithmetic.
   Coordinate *magnitude* therefore does not degrade the bound; only the
   *exponent range* could, and that is the second property.
2. **The failure mode Shewchuk's expansions actually have is overflow and
   underflow, not magnitude.** `DETRIA_Split` (line 216) computes
   `splitter * a`, needing `|a| < 2^996`; `incircleadapt` forms fourth-degree
   products of coordinate *differences*. At Web Mercator's ±2.0037e7 the largest
   intermediate is about `(4e7)^4 ≈ 2.6e30`, and even at 1e9 it is `1.6e37` —
   both nowhere near `1.8e308`. A foot-based CRS (~6.6e7 ft) is inside the same
   envelope. There is no projected CRS that gets close.

Note that the integer-overflow asserts near lines 1256/1275/1310 never apply:
`decltype(Point2::x)` is `double`, so `IsInteger` is false and detria takes the
floating-point path (`detria.hpp:1240-1242`), exactly as
`src/predicates/detria_exact.cpp` already records.

### The empirical check

A throwaway probe linked against the real `build/libterrain_predicates.a`,
compared against exact rational arithmetic in Python `fractions.Fraction`
(doubles are exactly representable as rationals, so this is a true oracle):

`orient2d`, triples built as `a`, `a+d`, `a+2d` on an integer lattice — exactly
collinear by construction — then nudged by 0, 1, 2 or 5 ulps so roughly half the
cases have a definite sign that only exact arithmetic can see:

```
mag       n=3000  kernel_wrong  FastKernel_wrong  filter_fallthrough  exactly_collinear
UTM33                        0                 3                1248               1238
WebMercator ±2.0037e7        0                 0                1269               1269
foot-ish 6.6e7               0                 0                1283               1283
1e9                          0                 0                1247               1247
1e15                         0                 0                1260               1260
```

`incircle`, four points drawn from the eight lattice points
`{(±p,±q),(±q,±p)}` on a common circle — exactly cocircular by construction —
with the query point nudged one ulp in or out:

```
mag       n=3000  kernel_wrong  definite(±1) answers  incircle_fallthrough
UTM33                        0                  1508                  2009
WebMercator                  0                  1474                  1579
foot-ish                     0                  1494                  1531
1e9                          0                  1457                  1550
```

Zero disagreements with exact arithmetic in 27,000 predicate calls, of which
~6,000 were definite Inside/Outside answers produced *by the backend after the
filter fell through*. That is a materially larger sample of definite backend
answers than the one test `01-predicates.md` flags as "currently the only thing
in the tree that observes a definite Inside/Outside answer coming out of the
backend", and it agrees.

This also incidentally clears a question nobody had asked: the tree sets no
`-ffp-contract` flag (`CMakeLists.txt` has only `-Wall -Wextra -Wpedantic
-Werror`), and Shewchuk-style expansion arithmetic is in principle sensitive to
FMA contraction of its `Two_Product` tails. On this toolchain, at `-O2`, it is
not being broken by it. That is a measurement on one compiler, not a guarantee
across all of CI — but it is evidence where there was none, and increment 2
already knows FMA contraction is a live concern in this tree.

**Conclusion: this concern is closed.** `project_structure.md` admitting
foot-based CRSs costs the kernel nothing.

---

## 4. Two stages or three

**Settled: two stands, and the corrected snap doctrine strengthens the case
rather than weakening it. The interval stage should not come forward.**

Two independent arguments, either sufficient.

**The backend is already the graded-precision stage.** detria's
`orient2dadapt` (detria.hpp:391) and `incircleadapt` (521) are Shewchuk's
*adaptive* routines: they re-check against `ccwerrboundB` at line 427 and
`ccwerrboundC` at 443, and `iccerrboundB`/`iccerrboundC` at 632/651, returning
as soon as the accumulated expansion settles the sign. Inserting an interval
stage between `FilteredKernel`'s level-A test and that call would add a third
precision level *in front of* a component that already has three. The saving it
could offer — resolving a sign at cheap precision without full expansion
arithmetic — is the saving stages B and C already make, one call deeper.

**On the data that actually falls through, an interval stage resolves nothing.**
The two fall-through families are:

- *exactly degenerate* (determinant exactly 0). An interval containing 0 never
  resolves; the interval stage is pure loss, and the measured cost of this case
  is already 4.3 ns because detria short-circuits it.
- *definite at the level of one or two ulps of the determinant*, which is what
  §5 shows a decimal snap spacing produces. The interval computed from double
  endpoints straddles zero by construction on a quantity that small. The
  interval stage is pure loss here too.

An interval filter earns its keep when fall-throughs are dominated by inputs
that are *moderately* near degenerate — say 10–1000 ulps out. Snapped breakline
data is not that distribution; it is bimodal at exactly-zero and at one-ulp.
The interval stage is the wrong tool for precisely the reason the deferral is
being questioned.

---

## 5. What the corrected snap doctrine invalidates that increment 1 did not catch

**Citations in this section are to the text as it stood at `760dbd9`.** `41054aa`
corrects both sites — `01-predicates.md:151` and `parallel_refinement.md:41-51`
— which is what §5.1 asked for, so the quoted wording is no longer there to find
and the present tense below reads as of that commit. This document is dated
evidence, not a live description of the tree.

### 5.1 "Collinear-but-unresolvable triples are the common case" is false at the recommended spacings

This sentence is now load-bearing in three places — `01-predicates.md:151`,
`parallel_refinement.md:41-51` and `05-noder.md`'s ruling 4 — and it is
measurably wrong in its most important half. (An earlier draft of this audit
counted four, adding `kernel.hpp`'s prose;
`grep -rnE 'unresolvable|breakline' include/terrain/predicates/` is empty, so
the claim never spread into the kernel headers at all.)

Take grid-collinear triples: `a`, `a + k₁·d`, `a + k₂·d` in **index** space,
mapped to world coordinates by `world(g) = g·s`, exactly as `SnapGrid` specifies.
Ask `DefaultKernel::orient2d`:

```
spacing  CRS       general-direction triples reported Collinear   axis-aligned   filter fell through
0.1 m    UTM33                              207 / 2976            1024 / 1024           1517 / 4000
0.1 m    WebMerc                            284 / 3003             997 /  997           1311 / 4000
0.01 m   UTM33                               39 / 3010             990 /  990           1062 / 4000
0.01 m   WebMerc                             47 / 3005             995 /  995           1043 / 4000
2^-4 m   UTM33                             3007 / 3007             993 /  993           4000 / 4000
2^-4 m   WebMerc                           3016 / 3016             984 /  984           4000 / 4000
```

At a decimal spacing, three points that are *exactly* collinear on the grid are
**not collinear in world coordinates** — 93 % of the time at 0.1 m, 98.7 % at
0.01 m. `fl(ix · s)` perturbs each coordinate independently by up to half an ulp,
and `s = 0.1` is not a dyadic rational, so `world` is not an affine map and does
not preserve collinearity. Axis-aligned runs are the exception and are exactly
collinear always, because one coordinate is literally equal across the three
points and the determinant is `0 · something`.

What is actually common in snapped data is therefore **fall-through with a
definite sign** (38 % of triples at `s = 0.1`), not fall-through with a collinear
answer. The kernel handles it correctly and deterministically — §3 verified
that — so this is not a kernel defect. It is a defect in the *description*, and
the description is what two increments are reasoning from.

### 5.2 The consequences land on increment 5a, not on increment 1

Taking 5.1 as given, these claims in the provisional `05-noder.md` do not hold:

- **`classify`'s collinear arm is nearly dead code on real input.** "If all four
  are `Collinear`" fires for axis-aligned pairs and for **at most** ~7 % of
  others — 7 % is the per-*triple* rate measured in §5.1, and the arm needs all
  four of a pair's triples, so the pair rate is lower and was never measured. The
  `Overlapping` relation, and with it `parallel_refinement.md:176`'s rule that a
  road snapped onto a river merges and keeps `is_river = true`, will essentially
  never trigger for a diagonal river. (Increment 7 replaced that one bit with a
  property *set*: the rule now stands at `:168` — the number above is refreshed
  — and says the merged edge is *both* road and river, not that it "keeps
  `is_river = true`". The audit's finding is unaffected: it is
  about how often `Overlapping` fires, not about what is carried when it does.
  The quotation above is left as it stood, this file being an audit dated to
  `e412a43`/#67.) The document names the *opposite* defect —
  reporting `Overlapping` where `Touching` is right — as "the single most likely
  defect in the function"; the measurement says the likely defect is on the other
  side.
- **A collinear-in-grid overlapping pair will be classified as `Crossing` or as
  `Disjoint` depending on which way the ulp-level perturbations fall**, and both
  are reachable. `Crossing` then drives `crossing_point` on a near-parallel pair
  — risk 4's regime — once per pair, per iteration.
- **`05-noder.md`'s ruling 4 quotes `01-predicates.md:151-152` verbatim as
  settled evidence.** It is inherited, not measured, and it is the inherited
  half that is wrong. The ruling's *conclusion* — "the design absorbs that
  because `DefaultKernel` is total and decisive, not because the fallback is
  rare" — survives unchanged and is, if anything, better supported.

None of this changes the kernel. It changes what 5a's suite must fixture.

### 5.3 The dyadic-spacing lever is worth more than `parallel_refinement.md` gives it

The document (lines 53-74) offers dyadic spacing "together with local
coordinates" as the way to get "the exactness property" back. Those are two
different properties and only one of them needs the local coordinates:

1. **`world` is an exact affine scaling, so grid-collinearity survives into
   world coordinates.** Needs dyadic spacing alone. `ix · 2⁻⁴` is exact for
   `|ix| < 2⁵³`, and Web Mercator at 6.25 cm spacing needs `|ix| ≈ 3.2e8 ≈ 2²⁸`.
   Measured above: 3007/3007 and 3016/3016 collinear, at CRS-origin anchoring,
   with no local coordinates anywhere.
2. **The determinants are exact in plain double, so the filter passes.** Needs
   coordinate magnitudes under ~2²⁶, hence local coordinates, and the document's
   2²⁶-not-2⁵³ argument is exactly right about this one.

Property 1 is what the noder's *topology* depends on — `on_segment`, the
collinear arm of `classify`, the hot-pixel cell corners of §1, `Overlapping`.
Property 2 is only a *performance* property, and §3 shows its absence costs
nothing in correctness: the row `2^-4 m` above falls through the filter
100 % of the time and answers `Collinear` 100 % of the time, at the 3.2 ns
measured in §2's third row.

That reframing is not a recommendation — spacing is a declared parameter and
`05-noder.md` risk 2 is right that no value is safe for all input. It is a
correction to the cost/benefit the document attaches to the lever: dyadic
spacing buys the property the noder actually needs, at CRS-origin anchoring, at
Web Mercator magnitudes, for free.

### 5.4 Minor, recorded because they are one-liners

- **Debug builds pay a third exact orientation per fallback.** detria's
  `incircle` asserts its counterclockwise precondition by calling
  `math::orient2d<Robust>` under `#ifndef NDEBUG` (detria.hpp:1236-1238). The
  §2 numbers are `-DNDEBUG`; the asan+ubsan Debug job pays more. Harmless, but
  anyone benchmarking a Debug build will see a number that is not the product's.
- **`01-predicates.md`'s "Deliberately excluded … intersection points are
  constructed in floating point and then snapped" is now *measured* correct**
  rather than assumed. §1's table is the evidence, and it is the piece of
  increment 1's reasoning the corrected doctrine was most expected to break.

---

## Verdict, restated

The kernel is sufficient. Nothing in
`include/terrain/predicates/` needs to change for the noder, and none of
increment 1's four rulings under audit falls:

| Deferral | Status |
|---|---|
| Two stages, not three | **Stands.** The backend is already adaptive (B and C bounds); snapped fall-throughs are bimodal at exactly-zero and one-ulp, where an interval stage resolves nothing. §4 |
| No exact constructions | **Stands, and is now measured.** Zero non-tie cell mismatches in ~45,000 crossings; the clamp never fired. §1 |
| `incircle_ccw` escape hatch deferred | **Stands.** 3.8× on the worst family, ~20 ns absolute, and no production caller exists to profile. §2 |
| `DetriaExact`'s bounds | **Hold**, to 1e9 m and to 1e15, orientation and incircle, against a `Fraction` oracle. §3 |

The one real finding is not a kernel finding: **`05-noder.md`'s `classify` asks
an exact-incidence question where snap rounding needs a hot-pixel proximity
question**, and at a decimal spacing those differ on ~93 % of grid-collinear
input. The fix is four `orient2d` calls against a cell's corners — the kernel
already provides everything it needs.

---

## Appendix: the probes

Throwaway, deliberately not added to the tree. Reproduced in full so the tables
above can be re-run. Build against the real static library:

```
cmake -S . -B build && cmake --build build -j
c++ -std=c++20 -O2 -DNDEBUG -Iinclude -o bench bench.cpp build/libterrain_predicates.a
c++ -std=c++20 -O2          -Iinclude -o probe probe.cpp build/libterrain_predicates.a
```

`probe.cpp` is a stdin filter over hex-encoded doubles (`o` = orient2d,
`i` = incircle) so the Python oracle can drive the real kernel; `bench2.cpp` is
§2's timing harness. The oracle uses `fractions.Fraction`, which represents
every double exactly, so it is a true oracle and not a higher-precision
approximation.

### probe.cpp
```cpp
// Throwaway audit probe. Reads hex-encoded doubles, prints kernel answers.
#include <terrain/predicates/default_kernel.hpp>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>

using namespace terrain;
using namespace terrain::pred;

static double h2d(const std::string& s) {
    std::uint64_t u = std::stoull(s, nullptr, 16);
    double d; std::memcpy(&d, &u, 8); return d;
}

int main() {
    std::string cmd;
    while (std::cin >> cmd) {
        if (cmd == "o") {           // orient2d: 6 hex doubles
            std::string t[6]; for (auto& x : t) std::cin >> x;
            double v[6]; for (int i = 0; i < 6; ++i) v[i] = h2d(t[i]);
            Point2 a{v[0],v[1]}, b{v[2],v[3]}, c{v[4],v[5]};
            std::printf("%d %d\n", static_cast<int>(DefaultKernel::orient2d(a,b,c)),
                                   static_cast<int>(FastKernel::orient2d(a,b,c)));
        } else if (cmd == "i") {    // incircle: 8 hex doubles
            std::string t[8]; for (auto& x : t) std::cin >> x;
            double v[8]; for (int i = 0; i < 8; ++i) v[i] = h2d(t[i]);
            Point2 a{v[0],v[1]}, b{v[2],v[3]}, c{v[4],v[5]}, d{v[6],v[7]};
            std::printf("%d\n", static_cast<int>(DefaultKernel::incircle(a,b,c,d)));
        }
        std::fflush(stdout);
    }
}
```

### bench2.cpp (§2)
```cpp
#include <terrain/predicates/default_kernel.hpp>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>
using namespace terrain; using namespace terrain::pred;
struct Q { Point2 a,b,c,d; };
template <class F> double ns(const std::vector<Q>& v, F f, int reps=20) {
    volatile int sink=0; auto t0=std::chrono::steady_clock::now();
    for (int r=0;r<reps;++r) for (auto& q:v) sink += f(q);
    auto t1=std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double,std::nano>(t1-t0).count()/(double(v.size())*reps);
}
int main(){
    std::mt19937_64 g(1); std::uniform_real_distribution<double> ux(-5e5,5e5), uy(7.0e6,7.9e6), uk(1,1e4);
    const double s=0.1;
    auto snap=[&](double x){ return std::llround(x/s)*s; };
    std::vector<Q> sep, colgen, colaxis, four;
    for (int i=0;i<20000;++i){
        Point2 a{ux(g),uy(g)};
        sep.push_back({a, Point2{a.x+uk(g),a.y+uk(g)}, Point2{a.x-uk(g),a.y+uk(g)}, Point2{a.x+uk(g),a.y-uk(g)}});
        // grid-collinear, general direction: snapped lattice, NOT exactly collinear
        long long ix=std::llround(ux(g)/s), iy=std::llround(uy(g)/s);
        long long dx=(long long)uk(g), dy=(long long)uk(g);
        colgen.push_back({Point2{double(ix)*s,double(iy)*s},
                          Point2{double(ix+dx)*s,double(iy+dy)*s},
                          Point2{double(ix+2*dx)*s,double(iy+2*dy)*s},
                          Point2{snap(ux(g)),snap(uy(g))}});
        // grid-collinear, axis aligned: exactly collinear
        four.push_back({Point2{double(ix)*s,double(iy)*s},
                        Point2{double(ix+dx)*s,double(iy+dy)*s},
                        Point2{double(ix+2*dx)*s,double(iy+2*dy)*s},
                        Point2{double(ix+3*dx)*s,double(iy+3*dy)*s}});
        colaxis.push_back({Point2{double(ix)*s,double(iy)*s},
                           Point2{double(ix+dx)*s,double(iy)*s},
                           Point2{double(ix+2*dx)*s,double(iy)*s},
                           Point2{snap(ux(g)),snap(uy(g))}});
    }
    auto Ok=[](const Q&q){return (int)DefaultKernel::orient2d(q.a,q.b,q.c);};
    auto Ik=[](const Q&q){return (int)DefaultKernel::incircle(q.a,q.b,q.c,q.d);};
    auto Of=[](const Q&q){return (int)FastKernel::orient2d(q.a,q.b,q.c);};
    auto If=[](const Q&q){return (int)FastKernel::incircle(q.a,q.b,q.c,q.d);};
    std::printf("%-34s %8s %8s %8s %8s\n","family","orient F","orient K","incir F","incir K");
    std::printf("%-34s %8.1f %8.1f %8.1f %8.1f\n","well separated",ns(sep,Of),ns(sep,Ok),ns(sep,If),ns(sep,Ik));
    std::printf("%-34s %8.1f %8.1f %8.1f %8.1f\n","grid-collinear, general dir",ns(colgen,Of),ns(colgen,Ok),ns(colgen,If),ns(colgen,Ik));
    std::printf("%-34s %8.1f %8.1f %8.1f %8.1f\n","grid-collinear, axis aligned",ns(colaxis,Of),ns(colaxis,Ok),ns(colaxis,If),ns(colaxis,Ik));
    std::printf("%-34s %8.1f %8.1f %8.1f %8.1f\n","4 grid-collinear pts (breakline)",ns(four,Of),ns(four,Ok),ns(four,If),ns(four,Ik));
}
```

### drive.py — the exact oracle shared by all Python probes
```python
import struct, random, subprocess, sys
from fractions import Fraction as F
PROBE="/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/8810ac89-7bb5-463c-bfdf-d30201db4f47/scratchpad/probe"
def hx(d): return "%016x" % struct.unpack("<Q", struct.pack("<d", d))[0]
EPS = 2.0**-53
OB = (3.0 + 16.0*EPS)*EPS
def filter_passes_orient(a,b,c):
    dl=(b[0]-a[0])*(c[1]-a[1]); dr=(b[1]-a[1])*(c[0]-a[0])
    det=dl-dr; perm=abs(dl)+abs(dr)
    return abs(det) > OB*perm
def exact_orient(a,b,c):
    v=(F(b[0])-F(a[0]))*(F(c[1])-F(a[1]))-(F(b[1])-F(a[1]))*(F(c[0])-F(a[0]))
    return (v>0)-(v<0)
def exact_incircle(a,b,c,d):
    # sign of lifted det, normalized by orientation of abc
    o=exact_orient(a,b,c)
    if o==0: return 0
    if o<0: b,c=c,b
    ad=(F(a[0])-F(d[0]),F(a[1])-F(d[1])); bd=(F(b[0])-F(d[0]),F(b[1])-F(d[1])); cd=(F(c[0])-F(d[0]),F(c[1])-F(d[1]))
    al=ad[0]**2+ad[1]**2; bl=bd[0]**2+bd[1]**2; cl=cd[0]**2+cd[1]**2
    v=al*(bd[0]*cd[1]-cd[0]*bd[1])+bl*(cd[0]*ad[1]-ad[0]*cd[1])+cl*(ad[0]*bd[1]-bd[0]*ad[1])
    return (v>0)-(v<0)

def run(lines):
    p=subprocess.run([PROBE],input="\n".join(lines),capture_output=True,text=True)
    if p.returncode!=0: sys.exit("probe died rc=%s %s"%(p.returncode,p.stderr[:400]))
    return p.stdout.split("\n")
```

### c.py / e.py — §3, exactness at magnitude
```python
exec(open("drive.py").read())
import math
random.seed(3)
mags=[("UTM33",5.0e5,7.9e6),("WebMerc",2.0037e7,2.0037e7),("foot",6.6e7,6.6e7),("1e9",1.0e9,1.0e9),("1e15",1.0e15,1.0e15)]
for name,MX,MY in mags:
    cases=[];lines=[]
    for _ in range(3000):
        a=(float(int(random.uniform(-MX,MX))), float(int(random.uniform(-MY,MY))))
        p,q=random.randint(1,1000),random.randint(-1000,1000)
        b=(a[0]+p,a[1]+q); c=(a[0]+2*p,a[1]+2*q)
        n=random.choice([0,0,1,1,2,5])
        cc=c
        for _k in range(n):
            cc=(math.nextafter(cc[0], cc[0]+random.choice([-1.0,1.0])*1e18), cc[1])
        cases.append((a,b,cc)); lines.append("o "+" ".join(hx(x) for pt in (a,b,cc) for x in pt))
    out=run(lines)
    bad=fast=ft=z=0
    for (a,b,c),o in zip(cases,out):
        k,f=map(int,o.split()); e=exact_orient(a,b,c)
        bad+= k!=e; fast+= f!=e; ft+= not filter_passes_orient(a,b,c); z+= e==0
    print(f"{name:8s} n=3000 kernel_wrong={bad} Fast_wrong={fast} fallthrough={ft} exactly_collinear={z}")

exec(open("drive.py").read())
import math
EPSI=2.0**-53; IB=(10.0+96.0*EPSI)*EPSI
def inc_filter_passes(a,b,c,d):
    ad=(a[0]-d[0],a[1]-d[1]); bd=(b[0]-d[0],b[1]-d[1]); cd=(c[0]-d[0],c[1]-d[1])
    al=ad[0]*ad[0]+ad[1]*ad[1]; bl=bd[0]*bd[0]+bd[1]*bd[1]; cl=cd[0]*cd[0]+cd[1]*cd[1]
    det=al*(bd[0]*cd[1]-cd[0]*bd[1])+bl*(cd[0]*ad[1]-ad[0]*cd[1])+cl*(ad[0]*bd[1]-bd[0]*ad[1])
    perm=(abs(bd[0]*cd[1])+abs(cd[0]*bd[1]))*al+(abs(cd[0]*ad[1])+abs(ad[0]*cd[1]))*bl+(abs(ad[0]*bd[1])+abs(bd[0]*ad[1]))*cl
    return abs(det)>IB*perm
random.seed(9)
for name,MX,MY in [("UTM33",5.0e5,7.9e6),("WebMerc",2.0037e7,2.0037e7),("foot",6.6e7,6.6e7),("1e9",1.0e9,1.0e9)]:
    cases=[];lines=[]
    for _ in range(3000):
        ox,oy=float(int(random.uniform(-MX,MX))),float(int(random.uniform(-MY,MY)))
        p,q=random.randint(1,3000),random.randint(1,3000)
        ring=[(p,q),(-q,p),(-p,-q),(q,-p),(q,p),(-p,q),(-q,-p),(p,-q)]
        random.shuffle(ring)
        sel=ring[:4]
        pts=[(ox+float(u),oy+float(v)) for u,v in sel]
        a,b,c,dq=pts
        mode=random.choice([0,0,1,2])
        if mode: dq=(math.nextafter(dq[0],dq[0]+(1e18 if mode==1 else -1e18)),dq[1])
        cases.append((a,b,c,dq)); lines.append("i "+" ".join(hx(x) for pt in (a,b,c,dq) for x in pt))
    out=run(lines)
    bad=defin=ft=oft=0
    for (a,b,c,dq),o in zip(cases,out):
        k=int(o); e=exact_incircle(a,b,c,dq)
        bad+= k!=e; defin+= e!=0
        oft += not filter_passes_orient(a,b,c)
        if exact_orient(a,b,c)!=0:
            aa,bb,cc=(a,b,c) if exact_orient(a,b,c)>0 else (a,c,b)
            ft += not inc_filter_passes(aa,bb,cc,dq)
    print(f"{name:8s} n=3000 kernel_wrong={bad} definite={defin} orient_fallthrough={oft} incircle_fallthrough={ft}")
```

### f.py — §5.1, does grid-collinearity survive `world()`
```python
exec(open("drive.py").read())
random.seed(17)
def world(ix,s): return ix*s
for sname,s in [("0.1 m",0.1),("0.01 m",0.01),("2^-4 m",2.0**-4)]:
  for name,MX,MY in [("UTM33",5.0e5,7.9e6),("WebMerc",2.0037e7,2.0037e7)]:
    cases=[];lines=[];axis=[]
    for _ in range(4000):
        ix0=int(random.uniform(-MX,MX)/s); iy0=int(random.uniform(-MY,MY)/s)
        ax = random.random()<0.25
        if ax:
            dx,dy=random.randint(1,10000),0
        else:
            dx,dy=random.randint(1,10000),random.randint(1,10000)
        k1=random.randint(1,50); k2=k1+random.randint(1,50)
        a=(world(ix0,s),world(iy0,s))
        b=(world(ix0+k1*dx,s),world(iy0+k1*dy,s))
        c=(world(ix0+k2*dx,s),world(iy0+k2*dy,s))
        cases.append((a,b,c)); axis.append(ax)
        lines.append("o "+" ".join(hx(x) for p in (a,b,c) for x in p))
    out=run(lines)
    coll=0;colla=0;na=0;ft=0
    for (a,b,c),ax,o in zip(cases,axis,out):
        k=int(o.split()[0]); ft += not filter_passes_orient(a,b,c)
        if ax: na+=1; colla += k==0
        else: coll += k==0
    print(f"{sname:7s} {name:8s} grid-collinear triples n=4000: general-direction reported Collinear {coll}/{4000-na}, axis-aligned {colla}/{na}, filter fell through {ft}")
```

### g.py + h.py — §1, snap(double) vs snap(exact)

Two notes for anyone re-running these, both found by `@reviewer` and settled
by running them rather than by reading:

- **`llround` here is not `std::llround`.** `floor(x+0.5)` differs from
  half-away-from-zero exactly at the ties — e.g. at `0.49999999999999994` — and
  the ties are what the `exact-tie` column counts. The headline (zero non-tie
  mismatches) does not depend on it, and ruling 1 of `05-noder.md` forbids
  `rint`/`nearbyint` in the product for a related reason, but a probe published
  so its tables can be re-run should not round differently from the thing under
  test. Treat the tie counts as indicative and the mismatch counts as exact.
- **`crossings = 20000` for the two generic families is `kept`, not `n`.** It
  was queried as a suspected transcription of the loop bound; re-running `h.py`
  reproduces `crossings=20000, non-tie cell mismatch=0` for both. The generic
  construction places `c` and `d` either side of the midpoint of `ab`, so every
  sample passes the exact crossing test and `kept == n` legitimately. The
  near-parallel families, where `den == 0` and near-misses do occur, are the
  ones that come in under `n`.
```python
from fractions import Fraction as F
import random, math
random.seed(23)
def llround(x):                      # NOT std::llround -- see the note below
    return int(math.floor(x+0.5)) if x>=0 else -int(math.floor(-x+0.5))
def exact_snap(fx, s):
    # round half away from zero on exact rational fx/s
    q = fx/F(s)
    n = q.numerator; d = q.denominator
    if n>=0: return (2*n+d)//(2*d)
    else: return -((-2*n+d)//(2*d))
def cross(ax,ay,bx,by): return ax*by-ay*bx
def test(M, s, near_parallel, n=20000):
    mism=0; far=0; clamped_out=0
    for _ in range(n):
        # two segments crossing, endpoints on the grid
        def gp():
            return (llround(random.uniform(-M,M)/s)*s, llround(random.uniform(-M,M)/s)*s)
        a=gp()
        L=random.uniform(10,5000)
        th=random.uniform(0,2*math.pi)
        b=(llround((a[0]+L*math.cos(th))/s)*s, llround((a[1]+L*math.sin(th))/s)*s)
        if near_parallel:
            dth=random.choice([1e-9,1e-7,1e-5])*random.choice([-1,1])
        else:
            dth=random.uniform(0.2,math.pi-0.2)
        mid=((a[0]+b[0])/2,(a[1]+b[1])/2)
        th2=th+dth
        c=(llround((mid[0]-L*math.cos(th2)/2)/s)*s, llround((mid[1]-L*math.sin(th2)/2)/s)*s)
        d=(llround((mid[0]+L*math.cos(th2)/2)/s)*s, llround((mid[1]+L*math.sin(th2)/2)/s)*s)
        d1=(b[0]-a[0],b[1]-a[1]); d2=(d[0]-c[0],d[1]-c[1])
        den=cross(d1[0],d1[1],d2[0],d2[1])
        if den==0: continue
        # exact check that they actually cross
        def eo(p,q,r): 
            v=(F(q[0])-F(p[0]))*(F(r[1])-F(p[1]))-(F(q[1])-F(p[1]))*(F(r[0])-F(p[0])); return (v>0)-(v<0)
        if not (eo(a,b,c)*eo(a,b,d)<0 and eo(c,d,a)*eo(c,d,b)<0): continue
        u=cross(c[0]-a[0],c[1]-a[1],d2[0],d2[1])/den
        px=a[0]+u*d1[0]; py=a[1]+u*d1[1]
        lox=max(min(a[0],b[0]),min(c[0],d[0])); hix=min(max(a[0],b[0]),max(c[0],d[0]))
        loy=max(min(a[1],b[1]),min(c[1],d[1])); hiy=min(max(a[1],b[1]),max(c[1],d[1]))
        if px<lox or px>hix or py<loy or py>hiy: clamped_out+=1
        px=min(max(px,lox),hix); py=min(max(py,loy),hiy)
        # exact
        Fden=F(d1[0])*F(d2[1])-F(d1[1])*F(d2[0])
        Fu=((F(c[0])-F(a[0]))*F(d2[1])-(F(c[1])-F(a[1]))*F(d2[0]))/Fden
        ex=F(a[0])+Fu*F(d1[0]); ey=F(a[1])+Fu*F(d1[1])
        gi=(llround(px/s),llround(py/s)); ge=(exact_snap(ex,s),exact_snap(ey,s))
        if gi!=ge:
            mism+=1
            if abs(gi[0]-ge[0])>1 or abs(gi[1]-ge[1])>1: far+=1
    return mism,far,clamped_out
for M,s,npar,lab in [(5e5,0.1,False,"UTM33 s=0.1 generic"),(5e5,0.1,True,"UTM33 s=0.1 NEAR-PARALLEL"),
                     (2.0037e7,0.1,False,"WebMerc s=0.1 generic"),(2.0037e7,0.1,True,"WebMerc s=0.1 NEAR-PARALLEL"),
                     (2.0037e7,0.01,True,"WebMerc s=0.01 NEAR-PARALLEL"),(2.0037e7,2**-4,True,"WebMerc s=2^-4 NEAR-PARALLEL")]:
    m,f,co=test(M,s,npar)
    print(f"{lab:32s} snap(double)!=snap(exact): {m}   off-by->1 cell: {f}   clamp actually fired: {co}")

exec(open("g.py").read().split("for M,s,npar")[0])
def test2(M,s,near_parallel,n=20000):
    kept=0;mism=0;far=0;co=0;ties=0;worst=0
    for _ in range(n):
        def gp(): return (llround(random.uniform(-M,M)/s)*s, llround(random.uniform(-M,M)/s)*s)
        a=gp(); L=random.uniform(10,5000); th=random.uniform(0,2*math.pi)
        b=(llround((a[0]+L*math.cos(th))/s)*s, llround((a[1]+L*math.sin(th))/s)*s)
        dth=(random.choice([1e-9,1e-7,1e-5])*random.choice([-1,1])) if near_parallel else random.uniform(0.2,math.pi-0.2)
        mid=((a[0]+b[0])/2,(a[1]+b[1])/2); th2=th+dth
        c=(llround((mid[0]-L*math.cos(th2)/2)/s)*s, llround((mid[1]-L*math.sin(th2)/2)/s)*s)
        d=(llround((mid[0]+L*math.cos(th2)/2)/s)*s, llround((mid[1]+L*math.sin(th2)/2)/s)*s)
        d1=(b[0]-a[0],b[1]-a[1]); d2=(d[0]-c[0],d[1]-c[1])
        den=cross(d1[0],d1[1],d2[0],d2[1])
        if den==0: continue
        def eo(p,q,r):
            v=(F(q[0])-F(p[0]))*(F(r[1])-F(p[1]))-(F(q[1])-F(p[1]))*(F(r[0])-F(p[0])); return (v>0)-(v<0)
        if not (eo(a,b,c)*eo(a,b,d)<0 and eo(c,d,a)*eo(c,d,b)<0): continue
        kept+=1
        u=cross(c[0]-a[0],c[1]-a[1],d2[0],d2[1])/den
        px=a[0]+u*d1[0]; py=a[1]+u*d1[1]
        lox=max(min(a[0],b[0]),min(c[0],d[0])); hix=min(max(a[0],b[0]),max(c[0],d[0]))
        loy=max(min(a[1],b[1]),min(c[1],d[1])); hiy=min(max(a[1],b[1]),max(c[1],d[1]))
        if px<lox or px>hix or py<loy or py>hiy: co+=1
        px=min(max(px,lox),hix); py=min(max(py,loy),hiy)
        Fden=F(d1[0])*F(d2[1])-F(d1[1])*F(d2[0])
        Fu=((F(c[0])-F(a[0]))*F(d2[1])-(F(c[1])-F(a[1]))*F(d2[0]))/Fden
        ex=F(a[0])+Fu*F(d1[0]); ey=F(a[1])+Fu*F(d1[1])
        # distance from a half-cell boundary, in cells
        def tie(f):
            q=f/F(s); r=q-((2*q.numerator+q.denominator)//(2*q.denominator) if q>=0 else -((-2*q.numerator+q.denominator)//(2*q.denominator)))
            return abs(abs(r)-F(1,2))
        istie = tie(ex)<F(1,10**6) or tie(ey)<F(1,10**6)
        gi=(llround(px/s),llround(py/s)); ge=(exact_snap(ex,s),exact_snap(ey,s))
        if gi!=ge:
            if istie: ties+=1
            else:
                mism+=1; worst=max(worst,abs(gi[0]-ge[0]),abs(gi[1]-ge[1]))
    return kept,mism,ties,far,co,worst
for M,s,npar,lab in [(5e5,0.1,False,"UTM33 s=0.1 generic"),(5e5,0.1,True,"UTM33 s=0.1 near-parallel"),
                     (2.0037e7,0.1,False,"WebMerc s=0.1 generic"),(2.0037e7,0.1,True,"WebMerc s=0.1 near-parallel"),
                     (2.0037e7,0.01,True,"WebMerc s=0.01 near-parallel"),(2.0037e7,2**-4,True,"WebMerc s=2^-4 near-parallel")]:
    k,m,t,f,co,w=test2(M,s,npar)
    print(f"{lab:28s} crossings={k:6d}  non-tie cell mismatch={m}  exact-tie disagreements={t}  clamp fired={co}  worst cell offset={w}")
```
