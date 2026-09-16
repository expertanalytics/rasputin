# Increment 1 — exact geometric predicates

Status: **shipped**. 1a merged as `d0946bc`, 1b as `aeb667a`.

Written retroactively. The design was settled across several agent
conversations, all of which are gone; this reconstructs it from the merged code,
the PR bodies and the decisions recorded in commit messages. It exists because
increment 1's open items were, until this file, carried nowhere but a chat
transcript.

Ships `include/terrain/predicates/{orientation,exact,kernel,detria_exact,default_kernel}.hpp`,
`src/predicates/detria_exact.cpp`, `lib/detria/`, and
`tools/check_detria_boundary.py`. Namespace `terrain::pred`.

Makes possible: a *sign you can trust*. At UTM33 magnitudes (x ≈ 5e5,
y ≈ 7.9e6) a naive orientation determinant cancels its leading digits and
returns a **wrong** sign — not an imprecise one — exactly when points are near
degenerate and the answer matters. One flipped sign corrupts a triangulation far
from where the mistake happened. Every invariant in `testing.md` for `cdt`,
`noding` and `flip` becomes machine-checkable only once this exists, which is
why it landed before the triangulator rather than after.

## Why it was split

1a shipped the three pure headers. 1b vendored detria and wired the backend.
The split happened because `detria_exact.hpp`, its TU and `default_kernel.hpp`
all depend on a library that was not yet vendored — writing them in 1a would
have landed dead, unbuildable code, and an increment cannot contain a test for
code it does not contain. Nothing downstream consumed `DefaultKernel` at the
time, so nothing was blocked.

## The types

```cpp
enum class Orientation : int { Clockwise = -1, Collinear = 0, CounterClockwise = 1 };
enum class Incircle    : int { Outside   = -1, Cocircular = 0, Inside          = 1 };
```

Underlying values are exactly -1/0/1 so `static_cast<int>` is a usable sign and
`reversed(o) == Orientation{-static_cast<int>(o)}`.

`orientation.hpp` is the vocabulary header and **includes nothing** — not even
`<limits>`. The layering claim is enforced by the preprocessor rather than by
convention. It holds `orientation_of_sign`, `incircle_of_sign`, `reversed`,
`is_left_turn`, `is_collinear`, all `constexpr noexcept`.

`is_left_turn` and `is_collinear` take an `Orientation`, not three points.
Taking points would force this header to depend on `Point2` *and* pick an
arithmetic policy — either hard-coding naive doubles, which makes the header
shared by both kernels quietly endorse the unfiltered one, or depending on a
kernel, which is an upward dependency. `is_left_turn(Collinear) == false`: a
hull walk written from one side keeps collinear points and from the other drops
them, so the degenerate case must not be reported as turning either way.
Callers wanting "left or straight" say so at the call site.

Both sign funnels are written `if (v > 0) … if (v < 0) … else degenerate`, which
gives the signed-zero and NaN behaviour structurally rather than by special
case: `orientation_of_sign(+0.0) == orientation_of_sign(-0.0) == Collinear`, and
`orientation_of_sign(NaN) == Collinear` as the only total choice — documented as
an unreachable state, not a supported one.

**There is no `reversed(Incircle)` overload, and a comment says why.** Its only
plausible use is reintroducing the bug described below.

## The concepts

```cpp
template <typename E> concept ExactPredicates = requires(Point2 a, Point2 b, Point2 c, Point2 d) {
    { E::orient2d(a, b, c)        } -> std::same_as<Orientation>;
    { E::incircle_ccw(a, b, c, d) } -> std::same_as<Incircle>;   // precondition: abc is CCW
};

template <typename K> concept GeometryKernel = requires(...) {
    { K::orient2d(a, b, c)    } -> std::same_as<Orientation>;
    { K::incircle(a, b, c, d) } -> std::same_as<Incircle>;
};
```

Spelled with **qualified static calls**, not the instance form, so the concept
itself requires models to be usable without an instance. Statelessness is
otherwise pinned only by a separate purity test, which future models would not
be covered by.

The exact backend returns **classifications, not sign-carrying doubles**. This
was a revision: the chosen backend returns enums, and the round trip through a
sign was an artifact of assuming a `predicates.c`-shaped API. The `_ccw` suffix
puts the precondition in the name.

`GeometryKernel` is a **semantic** contract. A concept whose models answer the
same question differently is a name for a signature, and everything downstream
templated on it would silently acquire a hidden dependency on which kernel it
was instantiated with, discoverable only as a wrong mesh.

## The kernels

`FastKernel` — naive double determinants, no filter, no fallback. Exists so the
suite can demonstrate the wrong answer that motivates the module.
`orient2d` is `orientation_of_sign(cross(b - a, c - a))`, calling `terrain::cross`
rather than restating the determinant, so the agreement the tests assert is
structural rather than something a future edit can break.

`FilteredKernel<E>` — inline static filter, out-of-line exact fallback only when
the sign is unsafe. Bounds are `orient2d_bound_a = (3 + 16ε)ε` and
`incircle_bound_a = (10 + 96ε)ε` with `ε = epsilon()/2`, derived rather than
literal, and **public only as a transcription check** so a typo in a magic
constant fails a `STATIC_REQUIRE` rather than silently disabling the fallback.
They are not an extension point. The permanent expressions stay private:
exposing them would freeze the expression against future adaptive refinement for
no caller's benefit, and a test asserting the exact threshold would restate the
implementation and pass by construction.

The filter is tested in the two implementation-independent directions — zero
fallbacks on well-separated input, at least one on collinear and on cocircular —
plus mutation testing that kills both a no-op filter and an always-fallback one.

### The incircle rule, and the bug that lived in it

```
incircle(a, b, c, d):
    o = FilteredKernel::orient2d(a, b, c)     // the FILTERED one, not E::orient2d
    Collinear        -> Cocircular            // no backend call, no lifted determinant
    CounterClockwise -> E::incircle_ccw(a, b, c, d)
    Clockwise        -> E::incircle_ccw(a, c, b, d)     // NO reversed()
```

An earlier draft wrapped the clockwise branch in `reversed()`. That inverts a
correct answer. Swapping `b` and `c` makes the triple counterclockwise, which
discharges the precondition; the circle through three points does not depend on
their order, so the reordered call already returns the right answer. The rule was
a leftover from the version where the backend returned a sign-carrying double —
there a row transposition flips the determinant and the caller must flip back,
but once the backend returns a classification that flip moved *inside* the
backend and the caller's compensation became a double negation.

Counterexample: `incircle((5,0), (-5,0), (0,5), (0,0))`. The triple is clockwise,
the circumcircle is centred at the origin with radius 5, and the query point is
the centre — so `Inside`. The old rule returned `Outside`. `@tester` caught it
before any implementation existed; the type-level tell is that `reversed` is
declared only for `Orientation`, so the rule never compiled.

**Normalizing first is not a tuning choice.** The lifted determinant is positive
iff the point is inside the circle *when the triple is counterclockwise*, so its
sign is uninterpretable without the orientation. A filter-first design would
certify a sign as safe and then not know what it meant. The circumcentre
alternative trades the orientation for a division by twice the same determinant —
strictly worse numerically, and still degenerate on collinear input.

Cost is roughly 1.2–1.3× on the filter-passing path, accepted for a total,
precondition-free predicate. It is worse on the *fallback* path: when the
orientation filter fails, the backend runs for orientation and possibly again
for incircle, and collinear-but-unresolvable triples are the common case in
snapped breakline data.

**`FastKernel::incircle` normalizes too**, using its own unfiltered `orient2d` —
never an exact one, since being unfiltered end to end is its defining property.
`FastKernel` is permitted to be *wrong* near degeneracy; it is not permitted to
answer a *different question*.

**Caching or memoizing orientations is prohibited at any scope.** It would put
shared mutable state into the one component in the tree that is provably pure,
and a contended cache line on the hot path. If a caller already knows a
winding, that belongs in the caller's data structure. The non-prohibited escape,
if profiling ever demands it, is a separately named `incircle_ccw` member on
`GeometryKernel` for callers maintaining a known winding — an interface addition
with the precondition in the name, reviewable on its own merits.

## The backend

`DetriaExact` in `detria_exact.hpp` (declaration only) and
`src/predicates/detria_exact.cpp` (the sole TU including detria).
`DefaultKernel = FilteredKernel<DetriaExact>`.

`lib/detria/detria.hpp` vendored at `8aa25f3e0dedf8d37623e7085b69f985e5845ad7`
(2025-03-02), sha256 `d6920a18…`, byte-identical to upstream, no local
modifications.

**MIT taken explicitly**, with `LICENSE-MIT.txt` vendored to record the
election. The library offers WTFPL or MIT at the user's choice and GitHub's own
detector reports the repo as WTFPL, which is not OSI-approved and is rejected by
some corporate license scanners. Recording the election costs one file.

### Why detria and not the alternatives

The deciding criterion was **no `exactinit()` and no file-scope mutable state**.
The classic Shewchuk `predicates.c` computes its error bounds into globals at
runtime: a static-initialisation-order hazard and a data race waiting for
parallel refinement. detria instead has
`static constexpr errorBounds = predicates::calculateErrorBounds<Scalar>()`,
compile-time and templated, deriving the same `(3+16ε)ε` and `(10+96ε)ε` bounds.

Rejected, recorded so the search is not repeated:

- **dengwirda/robust-predicate** — C++ and technically sound, but free only for
  private, research and institutional use, with commercial distribution
  requiring direct arrangement. More restrictive than the LGPL this project is
  leaving.
- **libigl-predicates**, tetgen ports — wrap `predicates.c`, inherit its global
  state.
- **danshapero/predicates** — C, not header-only; repo itself notes "public
  domain is not a license".
- **hporro/robust-predicates** — a Rust crate despite the name.

### The enum mapping

detria returns `math::Orientation{CW=0, CCW=1, Collinear=2}` and
`math::CircleLocation{Inside=0, Outside=1, Cocircular=2}` — not our values and
not our order, so the mapping is an explicit `switch`, **never** a cast. The
suite mutation-tests exactly the failures a cast produces. Call the robust
instantiation; `Robust = false` is killed by the suite.

detria's predicates are templated on any vector with `.x`/`.y`, so
`terrain::Point2` satisfies them directly — no adapter, no copy.

### The SIGTRAP hazard

detria's `math::incircle` (detria.hpp:1233 at the pinned SHA) asserts under
`#ifndef NDEBUG` that its first three points are counterclockwise, and
`detail::detriaAssert` (line 1404) writes to `stderr` then calls
`std::raise(SIGTRAP)` / `__debugbreak()`. It is called directly, not through a
macro, so there is no override short of `-DNDEBUG`.

Our asan+ubsan CI job is a **Debug** build. A normalization regression therefore
kills that job on a signal — no assertion text, no Release reproduction.

This is why `FilteredKernel::incircle` is total and precondition-free, and it is
tested two ways: a **forked probe** that runs the dangerous calls in a child and
reports a signal death as an ordinary test failure with a readable message, and
an **unguarded in-process test** that must simply not die. The unguarded one
carries a `READ BEFORE SIMPLIFYING` comment, because it passes by not crashing
at least as much as by its assertions and its assertions are *supposed* to look
mild.

The assertion was confirmed live rather than compiled out: a standalone probe
built against the sanitizer build's `libterrain_predicates.a`, calling
`DetriaExact::incircle_ccw` with a deliberately clockwise triple, died with exit
133 (signal 5). The green Debug suite is therefore evidence that normalization
works, not that the check is absent.

detria's *other* Debug assertions — the integer-overflow checks around lines
1150, 1256, 1275, 1310 — are **structurally vacuous** for us, not merely
unobserved: the body is behind
`if constexpr (std::is_floating_point_v<Scalar>) { return false; }`, so for
double coordinates only the counterclockwise assertion is live.

### The one-TU boundary, enforced three ways

detria.hpp is included by exactly one translation unit. That claim is what makes
the backend swappable, so it is enforced structurally rather than by convention:

1. **CMake** — `lib/detria` is PRIVATE on `terrain_predicates`, so the include
   directory never reaches the INTERFACE.
2. **`#ifdef DETRIA_HPP_INCLUDED → #error`** in the test suite, verified to fire.
3. **`tools/check_detria_boundary.py`**, wired into the CI governance job and
   negative-tested.

Prose rules erode; this one is checked on every push.

## A coverage trap worth remembering

On integer-coordinate input, the filter reaches the exact backend **only when
the four points are exactly cocircular** — where the answer is `Cocircular` and
an inverted Inside/Outside mapping is invisible. Measured: 93 backend calls
across 2000 lattice quadruples, **zero** with a definite answer. Two
cast-shaped mutations initially survived because of it.

The fix uses the circle through (1,0), (0,1), (-1,0) — exactly the unit circle —
with `nextafter(-1, -2)` and `nextafter(-1, 0)` one ulp outside and inside,
about 1.1e-16 relative, below the filter's bound and known by construction with
no oracle needed. The test `"DefaultKernel reports every circle location the
backend decides"` is currently **the only thing in the tree that observes a
definite Inside/Outside answer coming out of the backend**. Anyone editing it
should know that.

## Deliberately excluded

No segment intersection, no `on_segment`, no point location — those need a
`Segment2`, which is increment 2. No 3D predicates, no `orient3d`: the kernel is
2D forever. No exact *constructions* — signs only; intersection points are
constructed in floating point and then snapped. No interval-arithmetic stage:
two stages, not three, unless profiling demands it. No `operator<=>` on
`Point2`. No pybind11 exposure — predicates are internal.

## Open items carried forward

- **`-fvisibility=hidden` on `terrain_predicates`** when it stops being a static
  library. The *include* boundary is airtight but the *symbol* boundary is not:
  `nm` shows detria's template symbols with external linkage, so anything
  linking the target could declare and call internals, and a second vendored
  library with colliding names would be an ODR problem.
- **Rename `FastKernel`.** It is a full `GeometryKernel` model whose defining
  property is being wrong near degeneracy, and "Fast" reads at a call site like
  a legitimate performance option. `UnfilteredKernel`, or move it behind a
  testing namespace, once real consumers exist.
- **`DetriaExact::incircle_ccw` is a public footgun.** `ExactPredicates`
  requires it public, and with detria behind it any caller skipping
  normalization kills the process rather than returning a wrong answer. Inside
  the module that is fine — `FilteredKernel::incircle` is the only caller. Worth
  revisiting as the CDT lands.
- **`noexcept` asymmetry.** `FilteredKernel`'s predicates cannot be `noexcept`
  because `E::orient2d` may throw; `FastKernel`'s are. Invisible today, matters
  if something downstream wants a `noexcept` predicate in a hot loop.
- **No 128-bit integers anywhere**, here or downstream. GCC 16 rejects
  `__int128` under `-Wpedantic -Werror`, which is this repo's posture and the
  ubuntu CI compiler. First-party exact arithmetic must be limb-based or
  expansion-based.
- **`CLAUDE.md` §4 says `python tools/…`** where some machines have only
  `python3`. Harmless by hand, breaks if shelled out verbatim.
