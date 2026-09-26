# Increment 18 — the refinement scan walks row spans, not bounding boxes

Status: **implemented, in review.** Red `075113e`; green `62b2926`, `9389bcf`,
`3d395c3`. Output byte-identical to `e0578e5`; scan 15x faster on Ola's
quarter circle at 1 m (see "What landed"). Ola chose C1 (a), C2 (b) simplified to one
source-wide sentinel, C3 (a) and C4 (a) on 2026-09-26 (section "Ruled by Ola").
Written by `@architect` before `@tester`, per `docs/increments/README.md`
step 1, on branch `increment18-row-span-scan` off `increment17-mesh-stats`
(PR #95). Nothing committed.

**Closes.** Ola's idea (2026-09-26): "an iterator over the indices of the DEM
raster inside a triangle, in a row first fashion for better striding
properties". Earlier they called it a "branchless triangle iterator". Nothing
like it exists in this repository's history. After this increment,
`scan` (`include/terrain/refinement/scan.hpp`) visits only the DEM nodes in a
triangle's node set. For each row it computes the column interval once and
walks that row's contiguous memory. Today it tests every node in the
triangle's bounding box.

**Not closed.** The number of carving rounds on an outline-only tile
(`16-domain-polygon.md` R5). Rescanning triangles that had already converged
before a flip. SIMD beyond what the compiler does on its own. A virtual mosaic over several
DEM tiles. R6 makes that additive at the scan, on Ola's constraint, but builds
none of it. See "Not in scope".

## What was measured

### By the main session, with a scratch script

This compares DEM nodes in triangle bounding boxes with the nodes actually
inside the triangles.

| mesh | nodes in boxes | nodes inside | waste |
|---|---|---|---|
| quarter circle, start mesh (534 triangles, all with an off-node arc vertex) | 424 M | 7.1 M | 60× |
| quarter circle, final mesh at 1 m | 105 M | 7.1 M | 15× |
| 14b tile grid, output at 1 m | 59 M | 25.5 M | 2.3× |

`--stats` on the quarter circle at 1 m: `refine` took 3.29 s of 3.39 s, and
`refine: scan (parallel)` took 3.11 s of that.

### By `@architect`, with a scratch prototype

The prototype lived in the session scratchpad and never in the tree. It
includes the real `scan.hpp`, and it links `build-pyext/libterrain_predicates.a`.
It adds `scan_rows`, which is R1 and R2 below with option (a) of C1, for
non-void triangles only. Setup: `-O3`, one thread, Apple clang, and a
5000 × 5000 `Raster<float>` of smooth terrain plus noise.

| triangles | nodes inside | today | row spans | speed-up | results differing |
|---|---|---|---|---|---|
| fan of 534 slivers, off-node centre and arc (4 100 M in boxes) | 18.9 M | 17.5 s | 0.142 s | 123× | 0 |
| the same fan, all vertices rounded to nodes | 18.9 M | 6.07 s | 0.121 s | 50× | 0 |
| fan of 40, nodes (325 M in boxes) | 18.9 M | 0.515 s | 0.062 s | 8.3× | 0 |
| right-triangle grid, stride 20, nodes (2× waste) | 28.3 M | 0.164 s | 0.106 s | 1.5× | 0 |
| right-triangle grid, stride 4, nodes (about 6 nodes per triangle) | 37.4 M | 0.405 s | 0.465 s | **0.87×** | 0 |

"Results differing" compares `ScanResult` per triangle: `max_error` bit for
bit (by `memcmp`), plus `node` and `where`.

A second run used 200 000 random triangles on a 64 × 64 grid. The triangles
mixed nodes, off-node vertices 1e-12 and 1e-9 from a node, arbitrary fractional
vertices, slivers with a vertex 1e-7 off an edge's midpoint, and horizontal
edges. Result: 0 differences.

Two mutants were planted to check whether that probe can fail:

- **Mutant 1:** drop the vertex exclusion at the right end of the interval.
  The probe caught it on only **8 of 199 313** triangles. Comparing the
  argmax alone is weak evidence about the node set, and that is why T1 below
  compares sets.
- **Mutant 2:** change the forward correction from `< 0` to `<= 0`. It
  survived, and it is **equivalent**: the backward step that follows repairs
  it. It is recorded so that the mutation round does not spend time on it.

Three lessons for the design:

1. **The waste is the cost.** On slivers with an off-node vertex, today's
   scan spends 4.3 ns per box node on three exact orientations, almost all of
   them failing.
2. **Row setup has a price.** In the sliver fans the new scan costs about
   6.4 ns per inside node. In the fat fan it costs 3.3 ns. The difference is
   per-row overhead, and it is what makes stride-4 triangles slower than
   today. R2 names that overhead and the acceptance criteria bound it.
3. **The output does not change.** The node set and the per-node error formula
   are unchanged (C1 (a)), so the results are bit-identical.

## Rulings

### R1. The iterator: `for_each_row_span`, in its own header, with no raster

```cpp
// include/terrain/mesh/row_spans.hpp
struct RowSpan {
    std::uint32_t row;
    std::uint32_t c0, c1;       // inclusive, c0 <= c1; empty rows are never reported
    std::int8_t flat_edge;      // k when edge k is horizontal and lies on this row, else -1
};
template <std::invocable<RowSpan> F>
void for_each_row_span(const std::array<MeshVertex, 3>& v, F&& f);  // v counter-clockwise
```

- **It is a visitor, not an iterator class.** All per-triangle state
  (the edge constants and the node-only flag) lives on the caller's stack for
  one call. There is no begin/end pair and no sentinel. A test collects the
  spans into a vector in one line. Ola's "iterator" is the concept; this is its
  zero-overhead spelling. If they want a real range, see C3.
- **It lives in `mesh/`.** It is lattice geometry: it knows `MeshVertex`,
  `orient` and `orient_sign`, and it knows no raster. That lets its suite
  run with no DEM, which the testability pillar requires. `scan.hpp` is its
  only production caller.
- **Rows are visited top row first:** `ceil(min row)` to `floor(max row)`,
  ascending. Row indices grow southward.
- **The contract.** The union of the reported spans, `{(row, c) : c0 ≤ c ≤ c1}`,
  is exactly today's node set. That is every node `p` with all three
  `orient_sign(v[k], v[k+1], p) ≥ 0`, except nodes equal to a vertex. The
  spans come in strictly increasing row order.

**The per-edge bound.** Take edge `k` from `a = v[k]` to `b = v[k+1]`. In the
`(col, −row)` frame, `orient(a, b, (c, r))` is affine in `c`, with slope
`b.row − a.row` (`lattice_mesh.hpp`'s `orient`). Its exact sign is therefore
monotone in `c` along a row. The edge constrains the row as follows:

- slope > 0: a **lower** bound, `c ≥ …`;
- slope < 0: an **upper** bound, `c ≤ …`;
- slope 0: **no column bound**. The sign is the same for every `c`. If it is
  negative, the row is empty. If it is zero, the row lies on the edge, and
  `flat_edge = k`.

Then `c0 = max(ceil(col_min), lower bounds)` and
`c1 = min(floor(col_max), upper bounds)`. All three edges are applied. Where
one does not span the row it is redundant, but it is not wrong, and skipping it
would be one more case for the oracle to cover.

### R1a. Node-only triangles: exact integer floor and ceil

Here `dr = b.row − a.row` and `dc = b.col − a.col` in `int64`, and
`n = dc · (r − a.row)`.

- slope > 0: `c0 = a.col + ceil_div(n, dr)`
- slope < 0: `c1 = a.col + floor_div(n, dr)`

Both divisions are floor-correct for any signs. C++ `/` truncates toward zero,
so `floor_div` adjusts the quotient when the remainder is nonzero and its sign
differs from the divisor's. `ceil_div(a, b) = −floor_div(−a, b)`. Both are
exact, with no estimate and no correction. That is the whole node-only path.

**Overflow bound.** The two factors are bounded by the grid dimensions, not by
`uint32`: `|dc| ≤ cols − 1` and `|r − a.row| ≤ rows − 1`, so
`|n| < rows · cols`. The DEM holds `rows · cols` values in memory. Addressable
memory is at most 2^48 bytes, so `rows · cols < 2^48`, and `n`, `−n` and the
quotients are far below 2^63.

**A documentation defect, fixed in this PR.** `lattice_mesh.hpp` currently
says "uint32 coordinates, so a product of two differences fits int64". That is
false for the full `uint32` range: (2^32 − 1)^2 > 2^63. `orient` is safe for
the same grid-size reason as above. The comment is corrected to say so, per
`README.md`'s rule that a documentation defect is fixed in the PR that finds
it.

### R1b. Triangles with an off-node vertex: a double estimate, corrected exactly

For each non-horizontal edge:

1. **Estimate** `x = a.col + (b.col − a.col) · (r − a.row) / (b.row − a.row)`
   in double. Take `ceil(x)` for a lower bound or `floor(x)` for an upper
   bound. Clamp the result to the column box, extended by one:
   `[ceil(col_min) − 1, floor(col_max) + 1]`. A NaN or infinite estimate
   (which needs `b.row − a.row` to be denormal-small) falls to the box end.
2. **Correct, lower bound.** `s(c) = orient_sign(a, b, MeshVertex{c, r})` is
   exact (`DefaultKernel`), and it is non-decreasing in `c`.
   - While `c ≤ box_hi` and `s(c) < 0`: `++c`.
   - Then, while `c − 1 ≥ box_lo` and `s(c − 1) ≥ 0`: `--c`.
3. **Correct, upper bound.** The mirror image of step 2.

**Why it is exact.** On the box's integers, `{c : s(c) ≥ 0}` is a half-line,
because `s` is monotone. The forward loop leaves `c` at the first member of
that set at or after the estimate, or past the box. The backward loop then
leaves it at the first member overall. That point is the exact bound, so the
result matches three exact tests per node. It uses the same predicate on the
same `(c, −r)` doubles as today (`16-domain-polygon.md` R2's one frame for
inside tests). The node coordinates are small integers, exactly representable,
so both triangles sharing an edge see the same bound.

**Why it terminates.** Each loop moves `c` one step in one direction and stops
at a box end. The worst case is the box width, which is today's per-row cost.

**What it costs.** The estimate's error is a few ulps of the coordinates,
around 1e-12 of a cell at this tile's size. So a normal row pays exactly
**2 exact calls per edge, 6 per row**:

- lower bound: `s(c)` succeeds, then `s(c − 1)` fails;
- upper bound: the mirror image.

Horizontal edges pay 1 call. `DefaultKernel`'s filter answers almost all of
these on the fast path. Today's scan pays 3 per box node. In the off-node fan
that is 4 100 M × 3 calls, against 6 per row.

**Horizontal off-node edges.** When `b.row == a.row` in double, one call
`s(any c)` decides the whole row: empty, full, or `flat_edge`.

### R1c. Degenerate rows

- **The vertex exclusion.** A vertex is an extreme point of the triangle, so
  on its row it is an **end** of the row's interval, never an interior point.
  For each node vertex on row `r`:
  - if `c0` equals its column, `++c0`;
  - if `c1` equals its column, `--c1`.

  An off-node vertex is not a node, so nothing is excluded for it. This is a
  per-row check, not a per-node one.
- **A row that touches only a vertex.**
  - Node vertex: `[vc, vc]` becomes empty after the exclusion.
  - Off-node vertex: the exact bounds cross (`c0 > c1`).

  Both are empty, and empty rows are never reported.
- **A horizontal edge exactly on a row.**
  - Node endpoints: the span is the edge's node range with both endpoints
    removed.
  - Off-node endpoints: the span is the nodes between them.

  `flat_edge = k`.
- **A zero-width row inside a sliver.** `c0 > c1`: not reported.
- **Top and bottom rows.** `ceil` and `floor` of the fractional extents,
  exactly as today. Every vertex is inside the node rectangle, so every span
  is inside the grid. That is an assertion, as today.

### R1d. Tie rule and `where`

- **Ties.** Spans run in ascending row order, each ascending in `c`, and only
  a strictly larger error replaces the best. So ties still go to the smallest
  `(row, col)`. Nothing new is needed.
- **`where` for nodes strictly inside a span:**
  - `Edge{flat_edge}` when `flat_edge ≥ 0`;
  - otherwise `Inside`.

  Proof: suppose a non-horizontal edge had orientation exactly zero at an
  interior `c`. Then `c` would be that edge's bound, and one of `c ± 1` would
  be outside the triangle, which contradicts both being in the span.
- **`where` for the two end nodes.** Classified with today's exact
  zero-tests, in today's priority order (Edge0, then Edge1, then Edge2).
  - It runs only when an end node is recorded, either as the argmax or as the
    carve point. That is a handful of calls per triangle, not per row.
  - A node on two edge lines would be a vertex, so the priority order never
    actually decides anything. It is kept for identity with today's code.

### R2. The inner loop

```text
for each segment of row r over [c0, c1]:        // R6; one segment for one grid
for j in span: z = span[j]; e = |z − plane(c0 + j)|; e = nodata(z) ? 0 : e;
               if (e > best) { best = e; arg = j }        // rarely taken
```

- **No per-node inside test.** There is also no per-node index arithmetic:
  one `row(r)` per row, and a contiguous walk after it. `row(i)`, added in
  increment 12, gets its first caller.
- **NoData is branch-light.**
  - A NaN `z` makes `e` NaN, and `NaN > best` is false. NaN needs nothing.
  - A sentinel value is a select to 0. That is exact, because `best` starts at
    0 and only a strictly larger `e` replaces it.
  - The only remaining branch is the argmax update, which is rarely taken and
    well predicted.

  C2 decides whether the sentinel test can read `z` itself, or must go
  through `is_nodata(cell)`. The row is reached through R6's
  `for_each_row_segment`, which for one grid is exactly one `row(r)` subspan.
- **The plane is C1.** Under (a), today's formula is kept bit for bit:
  - **Node-only triangles:** the three orientations are exact integers, so
    they advance by `+dr_k` per column. That is exact, cheaper than three
    `orient` calls, and bit-identical.
  - **Mixed triangles:** the three double orientations are recomputed per node
    with today's expression. There is no incremental sum, because an
    incremental sum would change the bits.
- **Per-row overhead is the stride-4 regression's cause, and is bounded.**
  - Hoist every per-edge constant out of the row loop: `dr`, `dc`, the
    node-only flag, and the node-vertex rows.
  - Compute `row(r)` once per row.
  - Acceptance (below) requires that a fine node-only grid is not slower than
    today.
  - If hoisting is not enough, the next step is to skip edges whose row range
    excludes `r`. Skipping must be exact, and T1 covers it.
  - There is no fallback to the box walk in production (C4).
- **Branchless is the aim, measured, not a rule.** The acceptance numbers
  below decide. A mask-and-reduce argmax (a vectorised maximum per span, then
  a first-index search) is *not* built in this increment. It would be the next
  step only if the inner loop is shown to dominate after the row setup is
  fixed.

### R3. Carving and the void path use the same iterator

A void triangle's node set is the same set. `scan` passes the same spans to a
different body, which increments `uncovered` and keeps the nearest valid node
to a NoData vertex under the same strict-less, first-wins rule.
`refine.hpp`'s carving reads `ScanResult::node` and is untouched. One
iterator, two bodies, one oracle.

### R4. Parallelism and determinism are unchanged

- `for_each_row_span` is pure. It reads three vertices and holds its state on
  the stack, so `scan` stays pure. The thread model is still one result slot
  per triangle (`14-adaptive-refinement.md` R7), and nothing new is shared.
- Every per-triangle result is bit-identical to today's, by construction
  (R1's contract plus C1 (a)). So the serial split phase sees the same inputs,
  and the output mesh is bit-identical to increment 17's for every thread
  count. That is stronger than "deterministic": it is "unchanged".
- The TSan job has nothing new to see.

### R5. The old scan stays as a test oracle only

- Today's box walk moves verbatim into `tests/cpp/support/scan_oracle.hpp`
  as `scan_bbox`, together with a `bbox_node_set` that returns the membership
  set.
- It is not in `include/`, and nothing in production includes it.
- It is the reference that T1 and T2 compare against. It is frozen: no test
  may change it to agree with new code.

### R6. One DEM row is a sequence of segments, so a later mosaic is additive

**The constraint, from Ola (2026-09-26).** Domain polygons will often span
several DEM files. Stitching the covering tiles into one array does not scale
to large areas. Ola: "(a) does not work for large areas, so perhaps some sort
of indirection when needed". A later increment, a refresh of the parked
increment 15, adds a virtual mosaic: a `RasterSource` over several tiles that
resolves each node to its tile. With a mosaic, one row span can cross a tile
seam. This increment only makes that additive. It builds no multi-tile code.

**The ruling.**

- **The iterator itself is unchanged.** `for_each_row_span` knows no raster
  (R1) and still reports `(row, c0, c1)` in the global lattice. Tiles are
  not its concern.
- **The scan never touches `row(i)` directly.** It walks each span through one
  free function in `raster/`, which hands the body one or more contiguous
  segments in ascending column order:

  ```cpp
  template <typename T>
  struct RowSegment {
      std::span<const T> values;   // contiguous cells first_col .. first_col + size - 1
      std::uint32_t first_col;     // global lattice column of values[0]
  };  // no per-segment sentinel: the source has one nodata() (Ola's C2)
  template <RasterSource R, std::invocable<RowSegment<typename R::value_type>> F>
  void for_each_row_segment(const R& dem, std::size_t row, std::uint32_t c0, std::uint32_t c1, F&& f);
  ```

- **It is zero-cost for today's single grid.** The function is an
  `if constexpr` on a `requires` expression:
  - If `dem.row_segments(row, c0, c1, f)` exists, it forwards to it.
  - Otherwise it makes exactly one call:
    `f({dem.row(row).subspan(c0, c1 − c0 + 1), c0})`; the sentinel is read once from `dem.nodata()` (C2).

  Both `Raster` and `RasterView` take the second branch. That is a single,
  inlined call per row, with no loop over segments and no runtime dispatch.
  The per-node loop inside `f` is identical to R2's. The branch is resolved at
  compile time.

- **What a multi-tile source must provide later, named now:**
  - `row_segments(i, c0, c1, f)`, calling `f` once per tile the column range
    crosses, in ascending `first_col`. The segments must be non-empty, must
    not overlap, and must cover `[c0, c1]` exactly.
  - *Amended by Ola's rulings:* no per-segment sentinel. The source has one
    `nodata()` for the whole dataset (C2, simplified), which increment 15's
    U4 (refuse tiles with differing sentinels) makes sufficient.
  - **Tiles partition the dataset** (Ola: "The tile mosaics partitions the
    area, possibly with overlap ... So there should be no NoData unless out of
    the bounds of the total dataset."). A missing tile inside the dataset's
    extent is a data error the source refuses, not NoData. Only cells outside
    the whole dataset have no value, and a domain reaching there is already
    refused (16 U4). NoData *inside* a tile is still real (the Kartverket
    tile's top row is all −32767) and is carved as today.

- **The concept change, later and not now.** A mosaic cannot return a whole
  row as one span. `RasterSource` will then require "`row(i)` **or**
  `row_segments(i, c0, c1, f)`" instead of `row(i)`. That is a disjunction in
  the concept, or a split into a base concept plus two refinements. Code that
  only goes through `for_each_row_segment` needs no edit. In this increment
  the scan is the only walker, and it takes that route.

- **This increment's concept change.** `nodata()` becomes a requirement, so
  that the single-grid branch can fill `RowSegment::nodata`. `Raster` already
  has it. `RasterView` gains the same one-line accessor. This is C2 (b). Under
  C2 (a) the concept is untouched, `RowSegment` has no `nodata` field yet, and
  the loop calls `is_nodata(cell)`. The mosaic increment then adds the field.

**Why the tie rule and exactness survive segmenting.**

- **The node set does not change.** Segments partition `[c0, c1]` exactly, so
  the cells visited are the same cells, and the node set is still R1's.
- **The visiting order does not change.** Spans come in ascending rows,
  segments in ascending `first_col`, and cells ascending within a segment. So
  the walk is still strictly row-major. The body indexes by global column,
  `first_col + j`, and only a strictly larger error replaces the best. The
  smallest `(row, col)` still wins a tie, and the result is bit-identical
  however a row is cut.
- **T2 pins it with no mosaic.** A test-only `RasterSource` adapter over a
  `Raster` implements `row_segments` by cutting every row at fixed and random
  columns, sharing the source's one sentinel. T2 runs with it
  as well as with the plain raster, and the results must be identical. That
  proves the segmented path in this increment, before any mosaic exists.

**What a virtual mosaic touches beyond the scan.** These are noted here and
not designed; they belong to the mosaic increment.

- **Bilinear z for a start vertex whose cell straddles two tiles.**
  `vertex_z` reads four corners through `value_at`. If the mosaic's
  `value_at` resolves each node to its tile, this works unchanged. That
  requires the tiles to share one lattice, with aligned origin and equal
  `dx`/`dy`. Misaligned tiles need resampling, which is a different design.
- **NoData at tile seams.**
  - One dataset-wide sentinel (above); differing sentinels are refused at
    assembly (increment 15 U4).
  - A missing tile inside the extent is refused, not treated as NoData (Ola,
    above).
  - Overlapping tiles need a rule for which one wins (increment 15 U1, still
    open: data beats NoData, equal values accepted, different values refused).
- **Memory.**
  - A large area cannot be fully resident, so decode must be lazy or windowed.
  - A segment's span must stay valid for the duration of `f`, so a lazy source
    pins the tile row for that call.
  - The parallel scan calls the source from many threads at once
    (`14-adaptive-refinement.md` R7). A lazy source therefore needs
    thread-safe first-touch decode, or tiles decoded before the scan starts.
    Today's sources are pure reads, and that property must not be lost
    silently.
- **Geometry.** The mosaic presents one `RasterGeometry` for the union.
  `geometry()` keeps its meaning, so the lattice frame and every predicate in
  this file are unchanged.

## Ruled by Ola

**Ola chose C1 (a), C2 (b), C3 (a) and C4 (a) on 2026-09-26.** C2 (b) is taken
in its simplified form: `RasterSource` gains one `nodata()` for the whole
source, and segments carry no sentinel of their own, because tiles partition
the dataset and differing sentinels are refused (see R6). The options are kept
so the reasons stay on record.

### C1. The plane along a row: today's formula, or an affine update?

- **(a) Today's formula, bit for bit.** Integer orientations advance exactly;
  mixed triangles recompute today's double expression per node. The output is
  byte-identical to increment 17's, which gives a whole-pipeline golden test
  for free (T3).
- **(b) An affine plane per row,** `plane(c) = α_r + β·c`. This is fewer
  operations and easier for the compiler to vectorise. But its rounding
  differs in the last bits, so a near-tie can pick a different node, and the
  mesh can change. It is still deterministic, but T3 is lost, and the
  equivalence test must allow a tolerance on `max_error` and on argmax ties.

**Recommendation: (a).** The prototype shows that the waste, not the plane
arithmetic, is the cost, and (a) keeps the strongest test. (b) can be its own
measured increment later.

### C2. The NoData sentinel test in the inner loop

- **(a) Keep `dem.is_nodata(cell)` per node.** No interface change. Each node
  pays a second index computation and a second load from the same cache line.
- **(b) The sentinel travels with the segment** (R6). `RasterSource` gains
  `{ r.nodata() } -> std::same_as<const std::optional<value_type>&>`.
  `Raster` already has it, and `RasterView` gains a one-line accessor. The loop
  tests the `z` it already holds: `z != z || (nd && z == *nd)`. About 6 lines.
  There are no fakes of the concept in `tests/cpp` today.

  An earlier draft of this option put `is_nodata_value(v)` on the source
  instead. Ola's mosaic constraint rules that out. Tiles may have different
  sentinels, so a source-wide value test would be wrong at a seam, while a
  per-segment sentinel is right for both one grid and many.

**Recommendation: (b).** It makes the loop branch-light and index-free, and
it is the shape R6 needs anyway. If Ola prefers not to widen the concept in a
performance increment, (a) is correct and costs a little speed. The mosaic
increment would then have to add the per-segment sentinel itself.

### C3. The iterator's form

- **(a) The visitor `for_each_row_span(v, f)`,** as in R1.
- **(b) A forward range, `RowSpans{v}`,** usable in `for (RowSpan s : RowSpans{v})`.
  About 30 more lines (iterator, sentinel, state). The same code runs.

**Recommendation: (a).** The only consumer is one loop, and the visitor is
the smaller interface.

### C4. Small triangles

The stride-4 grid was 0.87× today's speed in the prototype, before any
hoisting.

- **(a) Row spans for every triangle.** Fix the per-row overhead (R2) and
  hold it to the acceptance bar.
- **(b) Keep the box walk below a box-area threshold,** as a second production
  path.

**Recommendation: (a).** (b) keeps two production membership paths, where the
point of R5 is that the box walk is only an oracle. It would also bring a
tuning constant. Revisit only if (a) misses the bar after hoisting.

## Prior art in `legacy/`

```sh
$ grep -rliE 'scanline|scan_line|rasteri[sz]e|span|bounding.?box|bbox' legacy/
legacy/rasputin/gml_repository.py
```

This is a false positive: `legacy/rasputin/gml_repository.py:202` takes a domain polygon's
`bounds` to size a download buffer. The legacy tree handed refinement to CGAL
and never scanned DEM nodes per triangle. Nothing is carried across, and
`@migration-expert` is not needed.

## Files

| file | what |
|---|---|
| `include/terrain/mesh/row_spans.hpp` (new) | `RowSpan`, `for_each_row_span`, `floor_div` / `ceil_div` (R1, R1a to R1c) |
| `include/terrain/refinement/scan.hpp` | the box loop replaced by spans; `where` at the ends; the header comment updated (R1d, R2, R3) |
| `include/terrain/mesh/lattice_mesh.hpp` | the overflow comment corrected (R1a) |
| `include/terrain/raster/row_segments.hpp` (new) | `RowSegment`, `for_each_row_segment` with its `if constexpr` single-grid branch (R6) |
| `include/terrain/raster/raster.hpp`, `view.hpp` | `nodata()` in the concept and on `RasterView` (C2 (b), R6) |
| `tests/cpp/support/scan_oracle.hpp` (new, test-only) | today's scan and node set, frozen (R5) |
| `tests/cpp/CMakeLists.txt` | `test_mesh_row_spans`, `prop_refinement_scan_equivalence` |
| `ROADMAP.md` | row 18 |

## Tests for `@tester`

**Invariant-critical: `test_mesh_row_spans` and T2.** A mutation round is
required on them.

- **T1. Node set = oracle, node for node** (`test_mesh_row_spans`). For each
  triangle, the union of the reported spans must equal `bbox_node_set`
  exactly, as a sorted list of `(row, col)`. The spans must come in strictly
  ascending rows, never be empty, and stay in the grid. Generators, seeded:
  - node-only triangles;
  - mixed triangles;
  - slivers (a vertex 1e-7 and 1e-12 off an edge, and collinear-but-one);
  - huge fans from one off-node centre to an off-node arc (the quarter-circle
    shape);
  - off-node vertices 1e-12 and 1e-9 from a node, on each side;
  - horizontal edges on a row, and horizontal edges between rows;
  - vertical edges;
  - triangles with a vertex on the grid border.

  Run the whole set under `dx ≠ dy` frames as well. The iterator must be
  independent of them, because it never sees world coordinates. The test pins
  that.
- **T1-deg. Named fixtures,** each with its expected span list written out:
  - a row touching only a node vertex, which must give no span;
  - a row touching only an off-node vertex;
  - a horizontal node edge on a row, which must drop both endpoints and set
    `flat_edge`;
  - a single-node triangle interior;
  - a triangle with an empty node set;
  - a negative numerator for `floor_div` / `ceil_div`, and a remainder of
    each sign.
- **T2. Scan equivalence** (`prop_refinement_scan_equivalence`). `scan` and
  `scan_bbox` must give identical `ScanResult`s on random meshes: the error by
  `memcmp`, plus `node`, `where`, `is_void` and `uncovered`. The DEMs must
  include:
  - NaN NoData;
  - sentinel NoData;
  - NoData at a vertex, which makes a void triangle and so exercises R3's
    carve point;
  - ties, from a DEM with repeated values, so that the smallest `(row, col)`
    rule is exercised.

  Per C1 (a), equality is exact.
- **T2-seg. Segmented rows (R6).** Run T2 again through a test-only adapter
  over the same `Raster`. The adapter implements `row_segments` by cutting
  every row at fixed columns and at seeded random ones, including cuts at `c0`,
  at `c1` and at every column. All segments share the source's one `nodata()`
  (Ola's C2, simplified). The results must
  be bit-identical to the plain raster's. This part is invariant-critical.
  Mutants: visiting segments in descending order, and an off-by-one in
  `first_col`.
- **T3. Refine unchanged, end to end.** For the real tile at `--tolerance 1`,
  with and without the quarter-circle domain, the output of `refine` must
  equal a golden digest of vertices, triangles and flags. The digest is
  recorded from increment 17's build **before** the green commit, and the
  commit that records it says so. Run with several thread counts.
- **Mutation targets:**
  - `ceil_div` → truncating `/`;
  - `floor_div` without the sign fix;
  - remove the vertex exclusion at `c0`, and separately at `c1`;
  - treat a horizontal edge with orientation zero as empty;
  - `ceil` → `floor` for the top row;
  - remove the clamp;
  - backward correction `s(c − 1) ≥ 0` → `> 0`;
  - an interior node classified `Inside` on a flat-edge row.

  One mutant is equivalent: the forward correction `< 0` → `<= 0` (see
  "What was measured"). It must not be counted as surviving.

  The prototype showed an argmax-only comparison catching a dropped exclusion
  in 8 of 199 313 triangles. T1's set comparison is what makes the round
  meaningful, so run each mutant against T1 first.

## Performance targets

Measured with the real tool, `rasputin mesh … --stats`, before (this branch's
base) and after, on the same machine. Three runs each, reporting the median
`refine: scan (parallel)` and `refine`.

| run | today's scan | expected after | bar |
|---|---|---|---|
| quarter circle, 1 m | 3.11 s | 0.15 to 0.4 s (8 to 20×) | ≥ 5× on scan |
| quarter circle, 10 m | to measure | similar ratio, smaller absolute | ≥ 3× on scan |
| 14b tile grid, 1 m (2.3× waste, node-only) | to measure | 1.2 to 1.6× | **not slower** |
| outline-only tile (optional; 109 s, 5 054 rounds) | to measure | per-round scan time falls | report only |

**Why the expected range.** The quarter circle wastes 60× at the start and
15× at the end, and every start triangle has an off-node vertex. The
prototype gave 8× at 17× waste with fat triangles, and 50 to 120× on slivers.
The 1 m run spends most of its rounds in between.

**Refine as a whole** should drop from 3.29 s to about 0.4 to 0.7 s. After
that, `split + flip (serial)` becomes the larger share. That is by design and
is the next target, not this increment's.

**The outline-only tile.** The iterator makes each of the 5 054 carving rounds
cheaper, because that tile's triangles are huge fans with large boxes. It does
not change how many rounds there are, and the round count is the cause
recorded in `16-domain-polygon.md` R5. Expect a large cut in scan seconds and
no change in `rounds`.

## LOC estimate

Counted in `CLAUDE.md` §2's unit.

| file | what | est. |
|---|---|---|
| `row_spans.hpp` | struct and signature ~10; div helpers ~10; integer bounds ~20; estimate and correction ~35; horizontal edges and vertex exclusion ~15 | ~90 |
| `scan.hpp` | loop replaced (+45 / −20), `where` at the ends ~10, comment | ~+35 net |
| `lattice_mesh.hpp` | comment only | ~0 |
| `raster.hpp`, `view.hpp`, concept | C2 (b): `nodata()` | ~6 |
| `row_segments.hpp` | `RowSegment`, `for_each_row_segment`, `if constexpr` branch (R6) | ~20 |
| | **total** | **~150** |

Increment 17's measured overrun was +66 %. On that basis the figure is about
250, which is under 700, so there is no split. `scan_oracle.hpp` is test code
and is not counted.

### What landed

Measured by `@developer`: **+168 / −42, net 126** production lines against
~150, under the estimate.

Timings, median of three runs with `--stats`, "before" built from `e0578e5`;
every output file is byte-identical before and after, and T3's digests match:

| run | scan before → after | speed-up | bar | refine | total |
|---|---|---|---|---|---|
| quarter circle, 1 m | 3.120 → 0.207 s | 15× | ≥ 5× | 3.294 → 0.381 s | 3.389 → 0.475 s |
| quarter circle, 10 m | 2.664 → 0.125 s | 21× | ≥ 3× | 2.674 → 0.135 s | 2.745 → 0.205 s |
| tile, no domain, 1 m | 0.100 → 0.076 s | 1.3× | not slower | 0.284 → 0.264 s | 0.431 → 0.414 s |

All three bars met. Settled at green: `where` is classified once per triangle
on the recorded node with the oracle's own exact zero-tests, instead of R1d's
derivation from `flat_edge` (which T1 still tests); a horizontal edge's row
is decided by one exact sign with no predicate call; a NaN estimate in the
mixed path falls to the widened box's low end and the exact correction finds
the bound.

## Acceptance

- T1, T2 and T3 are green. The mutation round is recorded in the green
  commit's message or the PR.
- The performance bars above are met, with the `--stats` tables before and
  after pasted into this file's "What landed".
- All gates in `CLAUDE.md` §4 are green, including TSan. CI is green.

## Not in scope

- **Rescanning converged triangles after flips.** 14b rescans every flipped
  triangle, even one whose node set it has already measured. That is a
  separate saving, with its own correctness argument.
- **The carving-round count on outline-only tiles** (`16-domain-polygon.md`
  R5).
- **SIMD beyond what the compiler does,** including a vectorised
  argmax-with-first-index. It is next only if R2's measurements show that the
  inner loop dominates.
- **An affine plane per row** (C1 (b)).
- **The virtual mosaic itself** (the parked increment 15, refreshed): a
  multi-tile `RasterSource`, the concept disjunction, seams, and lazy decode.
  R6 only shapes the scan so that it can be added.
- **Any change to what `refine` outputs.** R4 requires it to be unchanged.
