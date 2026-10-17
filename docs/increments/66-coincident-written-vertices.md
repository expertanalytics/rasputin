# Increment 66: no two vertices written at one point

**Status: draft, not reviewed** (`@architect`, 2026-10-17, on `37bea375`).
After increment 65, whose `flat_as_written` count it uses. Full path: it
changes C++ refine code, so `@perf`'s acceptance run applies. Ola ruled
zero-area slivers bugs (relayed by the main session, 2026-10-17).

## 1. The defect

Fixture W1 of increment 65 (`docs/increments/65-delaunay-as-written.md`,
section 4), outline k = 82, constraint feet **off** (`--no-constraint-feet`;
the CLI's default is on). Measured on master's `_core` (built 2026-10-16)
with the Python bindings, 2026-10-17:

- The outline's edge from (500094.6, 6599920.9) to (500088.6, 6599971.9)
  passes, in exact decimals, through the DEM node (col 9, row 4) at
  (500090.0, 6599960.0); as doubles it misses it by a rounding.
- `refine` (phase 1) inserts that node as vertex 7, off the constraint (no
  output constraint edge has it as an end): phase 1 has no rule for a node
  a hair off a constrained edge.
- `constraint_check_points` puts the edge's crossing of grid row 4 a
  rounding from col 9; `refine_strip` inserts it as vertex 49 on the
  constraint. Both write as (500090.0, 6599960.0), with z -3.0 and
  -2.99999999999589.
- Two output triangles, (49, 7, 11) and (49, 2, 7), have zero area as
  written (increment 65's prototype: `flat_as_written == 2` on k = 82, 0 on
  the other 99 outlines; with feet on, 0 on all 100; why feet on avoids it
  is not yet traced).

## 2. Prior art: legacy and literature

*Literature.* Snap rounding (J. Hershberger, "Stable snap rounding",
*CGTA* 46, 2013; D. H. Greene and F. F. Yao, "Finite-resolution
computational geometry", FOCS 1986) states the discipline: a vertex within
the rounding tolerance of a segment becomes a vertex of that segment. The
repository already applies it in the last stage: L12 of increment 15f (a
point Inside but within `coincidence_radius` of a constrained edge goes in
on it, `near_constraint` in `strip_scan.hpp`). No novelty claimed.

*Legacy.* Not yet checked (`git grep -n <pattern> legacy-archive -- legacy`
is the design's first command).

## 3. Design (draft)

**Choice.** Phase 1 takes L12: in `refine` (`refine.hpp`'s split), a node
Inside a triangle and within `coincidence_radius(g)` of one of its
constrained edges goes in on that edge (`split_edge` at the node, guarded by
`strip_fits` as in `point_loop`); the constraint then passes through the
node, and the strip generator, run after, puts no crossing a rounding from
it. With feet on, k = 82 has no flat triangle; whether the foot rule is the
reason is the design's to check.

Rejected: dropping a strip point that writes onto an existing vertex. Vertex
7 would stay a hair off the constraint and lie on it as written, so a
triangle along the constraint would be flat as written instead.

To settle in the design: the threshold (`coincidence_radius`, as L12, or the
written-point equality of increment 65's step 2); whether `refine_points`
(the geographic path) needs the same for its start; the probe set (W1 both
feet settings, the quick-check cases' `flat_as_written`, which increment 65's
prototype measured at 0 on all seven).

## 4. Tests (draft)

- **T1 (red).** W1, all 100 outlines, feet on and off: no two output
  vertices with equal written points, and `flat_as_written == 0`. Master:
  k = 82 with feet off fails (vertices 7 and 49).
- **T2.** Feet on: every output byte as before this increment on W1
  (the rule changes nothing where the foot rule already put the node on the
  edge).

## 5. Acceptance

`@perf`'s run on the default DEM, both domains; expected hash unchanged
(both run with feet on, and 65's prototype found no flat triangle on
either).

## 6. Size

About 10 net counted lines in `refine.hpp`; to be measured.
