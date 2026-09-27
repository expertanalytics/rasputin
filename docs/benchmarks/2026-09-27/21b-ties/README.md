# 21b tie classification: which incircle tests could QW2's integer path answer?

This is the measurement that `docs/increments/21-parallel-refine.md` §3 QW2 asks
for before 21b is designed in detail. It counts cases and does not time
anything.

**Finding.** In the refine loop, every exact-path incircle test is one of these
ties: a quad with four DEM-node corners, a lattice determinant of exactly 0, and
all of QW2's conditions met. That is 146,962 of 146,962 on the quarter circle
and 154,502 of 154,502 on the tile. The 146,962 is the same count the
2026-09-27 profile gave. The ties are not mostly grid rectangles: on the quarter
circle 38 % are axis-aligned rectangles or squares, and 62 % are other
lattice-cocircular quads. Of **all** incircle tests in the refine loop, QW2's
conditions hold on 99.973 % (quarter circle) and 100 % (tile).

## Method

- **Commit** `d3ee2ce` (branch `increment21b-lattice-incircle`, 21a's head).
  It was measured in a scratch worktree (`git worktree add --detach <wt> d3ee2ce`)
  with `scripts/instrument.patch` applied. That is a local patch and was never
  committed; it is kept here unapplied. `git apply --check
  scripts/instrument.patch` passes on `d3ee2ce`.
- **What the patch records.** `FilteredKernel::incircle_ccw` sets a flag
  saying whether it took the exact (`E::incircle_ccw`) path.
  `detail::must_flip` passes the four `MeshVertex` corners of each incircle
  call to `mesh::prof_ties::record` (new header
  `include/terrain/mesh/prof_ties.hpp`), in the order they went to the kernel
  (`a, b, c` counter-clockwise in the frame, then `d`), together with the
  kernel's answer and that flag. `refine` tags the phase: `p0` is
  `legalise_all`, `p1` the quality pass, `p2` the refine loop. At the end of
  `refine` it appends one JSON line to `$RASPUTIN_TIES_OUT`. For each call it
  records:
  1. the number of corners that are nodes (`MeshVertex::is_node`);
  2. for four nodes, the lattice determinant on `(col, -row)` translated to
     `d`, computed in `__int128` so it is exact at any spread, with its sign
     compared against the kernel's answer;
  3. QW2's conditions, per call: `dx == dy`; every coordinate difference from
     `d` at most 2^14 nodes; the frame exact, meaning `col * dx` and
     `row * dy` are exact for all four corners (checked with `fma(a, b, -a*b)
     == 0`). Once per refine call it also records the sufficient condition
     QW2 states: the significant bits of `dx` plus
     `bit_width(max(rows, cols) - 1)` are at most 53;
  4. for four nodes and a zero determinant, the shape. A quad is an
     axis-aligned square or rectangle when it has two distinct rows and two
     distinct columns, and the size is recorded. It is a rotated square or
     rectangle when its diagonals have equal midpoints and equal lengths. It
     is an isosceles trapezoid when two disjoint chords are parallel: along
     rows or columns, or in another direction. Anything else is "other
     cyclic". A quad with three collinear corners is counted separately.
- **Probes that can fail.** `scripts/classifier_selftest.cpp` checks the shape
  classifier on hand-made quads, one per class, including a circle of radius
  5. `scripts/counter_selftest.cpp` feeds `record` a case for each refusal or
  mismatch counter: an inexact frame at `dx = 0.1`, `dx != dy`, a spread above
  2^14, a nonzero lattice determinant on the exact path, a lattice sign that
  disagrees with the kernel, and three node corners. Every one of those
  counters fired. To build either one:
  `c++ -std=c++20 -Wall -Wextra -Wpedantic -Werror -I<wt>/include <file>`.
- **The instrumented build produces the same output.** The quarter circle's
  binary VTK hashes to `ff705683…`, the value the serial profile got. Hashed
  from `POINTS` onward, the way `bench.py` does it, the ASCII meshes match the
  21a acceptance run (`docs/benchmarks/2026-09-27/21a/run.json`): quarter
  circle `1e531976…`, tile `11741a81…`. Rounds, insertions and flips match as
  well.
- **Build**: Release, configured the way `bench.py`'s `build()` does it
  (`-DCMAKE_BUILD_TYPE=Release -DRASPUTIN_BUILD_PYTHON=ON
  -DRASPUTIN_BUILD_TESTS=OFF`), in `<wt>/build-instr`. The `_core` target was
  built (exit 0) and copied into a fresh `build-instr/pkg/tin_engine/`
  (symlinks to `src_python/tin_engine/*` plus the `.so`). The `.so` sha256 is
  `731a9609…`.
- **Runs.** The inputs are the 1 m benchmark: DEM
  `tests/fixtures/dem_archive/7908_3_10m_z33.tif` (5051 × 5051 nodes,
  `dx = dy = 10`), tolerance 1, all other CLI defaults. The domains are
  `docs/benchmarks/2026-09-26/quarter.geojson` and the whole tile (no
  `--domain`). The command runs through `scripts/prof_driver.py`, copied
  verbatim from `serial-profile/scripts/`:

  ```bash
  RASPUTIN_TIES_OUT=<out>/ties_<domain>_t<N>.jsonl .venv/bin/python scripts/prof_driver.py \
      --pkg <wt>/build-instr/pkg --threads <N> --repeat 1 -- \
      mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 \
      [--domain docs/benchmarks/2026-09-26/quarter.geojson] --out <out>/<domain>_t<N>.vtk --binary
  python scripts/tables.py data/ties_quarter_t1.jsonl data/ties_tile_t1.jsonl > data/tables.md
  ```

  Each domain was run at 1 and 8 threads. The count files are byte-identical
  between the two thread counts (`cmp`), as they should be, since the split
  phase is serial.
- **Machine**: Apple M1 Max, macOS 27.0. **Power: AC**, 100 %
  (`data/pmset-start.txt`). The counts do not depend on power; it is recorded
  because the persona requires it.
- **Timings** printed by the driver are inflated by the instrumentation, which
  does a `std::map` update per call. None is reported.

## Results

`data/tables.md` is produced by `scripts/tables.py` from `data/ties_*_t1.jsonl`.
The rows for the refine loop:

| measure | quarter | tile |
|---|---:|---:|
| incircle calls | 1,615,895 | 1,639,699 |
| exact path | 146,962 (9.095 %) | 154,502 (9.423 %) |
| exact path, answer Cocircular | 146,962 | 154,502 |
| exact path, four node corners | 146,962 | 154,502 |
| exact path, lattice det = 0 | 146,962 | 154,502 |
| exact path, lattice det != 0 | 0 | 0 |
| exact path, all QW2 conditions hold | 146,962 | 154,502 |
| filtered path, four node corners | 1,468,493 | 1,485,197 |
| filtered path, all QW2 conditions hold | 1,468,493 | 1,485,197 |
| **all calls QW2 would answer** | **1,615,455 (99.973 %)** | **1,639,699 (100 %)** |
| calls with fewer than four node corners | 440 | 0 |
| lattice sign differs from the kernel's answer (four nodes) | 0 | 0 |
| max coordinate difference from `d`, four-node calls | 951 nodes | 88 nodes |
| max coordinate difference from `d`, exact-path calls | 146 nodes | 48 nodes |
| frame-exact sufficient condition | 3 + 13 = 16 ≤ 53 | 3 + 13 = 16 ≤ 53 |

The per-call conditions never refused a four-node quad. Every four-node call
had `dx == dy`, a spread at most 2^14 and an exact frame, on both domains and
in every phase.

The shapes of the exact-path quads in the refine loop:

| shape | quarter | tile |
|---|---:|---:|
| other cyclic quad (no parallel sides) | 47,547 (32.4 %) | 47,269 (30.6 %) |
| axis-aligned square | 39,977 (27.2 %) | 39,911 (25.8 %) |
| rotated square | 17,100 (11.6 %) | 17,203 (11.1 %) |
| axis-aligned rectangle, not square | 16,261 (11.1 %) | 16,567 (10.7 %) |
| isosceles trapezoid, parallel sides on rows or columns | 12,697 (8.6 %) | 20,146 (13.0 %) |
| isosceles trapezoid, other direction | 11,169 (7.6 %) | 11,080 (7.2 %) |
| rotated rectangle, not square | 2,211 (1.5 %) | 2,326 (1.5 %) |
| three corners collinear | 0 | 0 |

The most common axis-aligned sizes are 1×1 (37,887 on the quarter circle), 1×2
(13,119), 2×2 (1,845), 1×3 and 2×3. `data/tables.md` has the per-size lists and
the other phases, which also go through `must_flip`:

- On the quarter circle, `legalise_all` makes 533 calls, of which QW2 would
  answer 122. The quality pass makes 7,677, of which QW2 would answer 4,622.
  The remaining 411 and 3,055 calls have an off-node corner (the arc
  vertices).
- On the tile, QW2 would answer all of `legalise_all`'s 48,133 calls, 16,129
  of them exact ties: 15,876 are 40×40 squares and 252 are 10×40 rectangles,
  from the starting mesh. It would also answer all 4,017 of the quality pass's
  calls.

## What is measured and what is inferred

Measured:

- QW2's premise holds. All 146,962 (quarter circle) and 154,502 (tile)
  exact-path ties in the refine loop have four node corners and a zero int64
  lattice determinant on `(col, -row)`, and each of QW2's per-call conditions
  holds for each of them.
- On every four-node call, in every phase and on both domains (3.3 million
  calls), the lattice sign agreed with `FilteredKernel<DetriaExact>`'s answer
  on the frame doubles. This is QW2's bit-identity argument checked on this
  input. It does not prove the argument.
- The share of all refine-loop incircle calls that QW2's conditions would
  admit: 99.973 % and 100 %.
- The ties are lattice-cocircular quads of many shapes. Only 38 % (quarter
  circle) and 37 % (tile) are axis-aligned rectangles.

Not measured, and so inferred or open:

- **The time QW2 saves.** Answering 99.97 % of calls says nothing about the
  cost of an int64 determinant compared with the filtered double determinant
  plus its permanent. Whether QW2 reaches past the 5.0 % exact-path cost into
  the rest of the 13.2 % flip-test cost has to be measured on an
  implementation.
- **Why the ties occur.** One conjecture is that small cocircular quads are
  common because refinement inserts nodes near one another on a grid. That
  has not been checked.
- **Other DEMs.** This covers one DEM, with `dx = dy = 10`. A DEM with
  `dx != dy`, or with a `dx` like 0.1 whose frame is not exact, would get
  nothing from QW2 by design. How common such DEMs are was not surveyed.

## What bears on QW2's design

1. **The general determinant is needed.** An axis-aligned-rectangle shortcut
   would handle only 38 % of the ties, so QW2's general int64 determinant is
   the right shape. The "points on a circle of radius 5" case in the planned
   suite matches what the data shows.
2. **The spread bound never binds here.** The largest spread was 951 nodes
   against a limit of 2^14. By arithmetic, not measurement: for a spread of at
   most 2^12, every intermediate of the plain double determinant on the
   integer `(col, -row)` differences is an integer below 2^53, so it would be
   exact in doubles too. The choice between int64 and doubles can be decided
   by timing, not by correctness.
3. **The conditions are per refine, not per call, on this data.** The per-call
   checks for `dx == dy` and an exact frame never refused anything. The
   once-per-refine decision QW2 already plans is enough, leaving the node check
   and the spread check as the only per-call tests.
4. **`must_flip` is also used outside the refine loop.** `legalise_all`
   (16,129 ties on the tile) and the quality pass also go through it, so
   QW2's call at the top of `must_flip` covers them without extra work.

## Not in the repository

The meshes (15-17 MB binary, 22-23 MB ASCII) and the scratch worktree with its build
stayed in the session scratchpad and are gone with it. To regenerate them,
apply `scripts/instrument.patch` to `d3ee2ce` and follow the Method above.
