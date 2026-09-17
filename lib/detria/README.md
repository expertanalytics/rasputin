# detria (vendored)

Third-party, unmodified. Only the exact geometric predicates are used:
`detria::math::orient2d` and `detria::math::incircle`, both instantiated with
`Robust = true`. The triangulator itself is not used and will not be -- this
project builds its own CDT.

## Provenance

* Upstream: <https://github.com/Kimbatt/detria>
* Pinned commit: `8aa25f3e0dedf8d37623e7085b69f985e5845ad7` (2025-03-02)
* Vendored file: `include/detria.hpp` at that exact SHA, fetched from
  `raw.githubusercontent.com`, not from `master`.
* sha256 of `detria.hpp`:
  `d6920a1817cc71d2604762b8a95277c00f0b9603a6ace0981ceb69a7ab4af013`

Re-vendoring means bumping the SHA here, refetching, and rerunning the
predicates suites. Do not fetch from a branch name.

## Licence: MIT, by election

detria is dual-licensed **WTFPL or MIT, at the user's choice**. We elect **MIT**
and vendor `LICENSE-MIT.txt` alongside the header to record that election.

This is not a formality. GitHub's own licence detector reports the repository as
WTFPL, which is not OSI-approved and is rejected outright by some corporate
licence scanners. Recording the election in-tree, next to the file it applies
to, is what keeps a downstream audit from inheriting the WTFPL classification
from upstream metadata we do not control.

## No local patches

`detria.hpp` is byte-identical to upstream and must stay that way. It is
compiled through a `SYSTEM` include directory precisely so that our
`-Wall -Wextra -Wpedantic -Werror` posture does not create pressure to patch it:
it is C++17-era third-party code, and a local fix would have to be re-applied,
by hand, at every version bump.

If a defect needs fixing, fix it upstream or work around it in
`src/predicates/detria_exact.cpp`, which is ours.

## One translation unit

`detria.hpp` is 4500 lines, pulls in `<iostream>`, `<sstream>` and `<csignal>`,
and its Debug assertions raise `SIGTRAP`. The rule is therefore **one
translation unit per backend, and no header under `include/` at all**. Two TUs
include it today:

- `src/predicates/detria_exact.cpp` — the exact predicate backend (increment 1b);
- `src/cdt/detria_backend.cpp` — the CDT backend (increment 4).

A third would need a reason; the point of the rule is that the library is
swappable per backend, not that the count stays at one. That rule is enforced
three ways:

1. `lib/detria` is on the `terrain_predicates` and `terrain_cdt` targets'
   include paths `PRIVATE`, so it never reaches any consumer's interface;
2. `tools/check_detria_boundary.py` parses includes and fails the build if the
   rule is broken (wired into the `governance` CI job);
3. `tests/cpp/unit/test_predicates_detria_exact.cpp` carries an
   `#ifdef DETRIA_HPP_INCLUDED -> #error` guard.

## Why detria and not the alternatives

The search was for a permissively licensed, header-only, thread-safe,
initialisation-free implementation of Shewchuk's adaptive `orient2d` and
`incircle`. This note exists so that search is not repeated.

* **dengwirda/robust-predicate** -- the licence forbids commercial use without a
  separate arrangement. That is *more* restrictive than the LGPL-encumbered CGAL
  this project is in the middle of leaving, so adopting it would defeat the
  migration.
* **libigl's and tetgen's predicate ports** -- both wrap Shewchuk's original
  `predicates.c`, which requires a call to `exactinit()` and stores the derived
  error bounds in mutable globals. Global state initialised lazily is exactly
  what the per-subdomain parallel refinement cannot have, and
  `tests/.../test_predicates_detria_exact.cpp`'s concurrency case exists to keep
  it out.
* **danshapero/predicates** -- C, not header-only, and would add a second
  language and a build step to the core.
* **hporro/robust-predicates** -- Rust; not linkable into this core without an
  FFI layer we have no other reason to have.

detria computes its error bounds as `constexpr`, keeps every expansion on the
stack, and has no global state, so it satisfies all four requirements.
