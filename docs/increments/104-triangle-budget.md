# Increment 104: the triangle budget

**Status: 104a built and accepted by `@perf` (section 15); code review
round 1's two findings answered, round 2 approved; 104b not started. Design
round 2 (2026-10-25) answers the round-1 design review.** Designed by `@architect`, 2026-10-25, on
master `97759b7b`, worktree `tri-budget`. Two pull requests from this one
file, each with its own red and green commit and review: **104a**, the budget
(Python only; it drives `refine`, so `@perf` times it), then **104b**, the
error summary (a C++ scan pass under `include/terrain/refinement/`, so
`@perf`'s acceptance run applies, `docs/increments/README.md`, *Acceptance*).
Number 104 because 103 is the highest under `docs/increments/` over every
local branch, remote-tracking ref and worktree, checked 2026-10-25
(`git for-each-ref` + `git ls-tree`; 100 to 103 are the private package's).

The showcase is Wednesday 2026-10-28 on Kultsjön (section 8); the plan and
the cut line are section 12.

## 1. What is asked, and why

Today `rasputin mesh --tolerance T` refines until every DEM node is within
T metres of the mesh, and the triangle count is whatever that takes. Ola's
idea (2026-10-17, the ideas file's *A triangle budget*): choose the count
instead, and let rasputin find the smallest error it can reach with it. His
rulings, 2026-10-22, on the study's four questions ("I think a soft and a
hard tolerance is good. Generally, looks fine."): about N triangles and never
more than M (M > N); the tolerance is not given but found and reported, from
the worst triangle (the largest error), with the root-mean-square error
beside it; one tolerance for the whole basin; triangles counted, not
vertices. On 2026-10-25: "yes, design it after the runs" and, for the demo
basin, "Let's go for the smallest": Kultsjön, the smallest catchment of his
northern Sweden delivery of 2026-10-21.

In two lines: `--triangles N` is a tolerance rasputin finds for you. The mesh
it writes is exactly the mesh `--tolerance T` writes at the tolerance it
found, so every guarantee of today's tolerance run holds at that T.

**What the showcase claims, and what it does not.** You ask for a size and
get at most M triangles; the tolerance is found and reported; the mesh is
identical to `--tolerance` at that tolerance. The budget mesh is *not* a
better mesh than the tolerance mesh: both lie on one error-against-count
curve, and a lower largest error always comes with more triangles (on
Kultsjön, 9.79 m costs 2.4 % more triangles than 10 m). That curve is the
demonstration (section 8).

## 2. Prior art: legacy and literature

### Literature

- **Garland and Heckbert 1995**, "Fast Polygonal Approximation of Terrains and
  Height Fields", CMU-CS-95-181 (read in full for the study; sections 2, 4
  and 5.4). Greedy insertion continues "until some specific goal is
  achieved, usually reaching a desired error threshold or exhausting a point
  budget" (2.1); the goal may be the point count, the L∞ error or the L2
  error (4.1). **The budget stop itself is 1995 textbook, not new.** Their
  sequential algorithm inserts one worst point per pass; `refine.hpp` is
  their parallel *constant threshold* variant (5.4: every candidate above the
  tolerance per pass), which they found "much worse" than sequential. The
  study's probes (section 4.2 and 9) found today's rounds spend 2-3 % more
  triangles than a near-worst-first order on a 10 m tile, so this design
  keeps the loop and chooses the stop by a search over T (section 3), outside
  the loop; hmm's `-t` and VTK's count stop act inside theirs.
- **Terra, VTK `vtkGreedyTerrainDecimation`, hmm, delatin** (READMEs and
  class page read for the study): all offer the count stop; hmm's `-e`, `-t`,
  `-p` ("the first one reached is used"); delatin reports both the maximum
  and the RMSD, which 104b mirrors.
- **Agarwal and Suri 1998** (SIAM J. Comput. 27(4), known through arXiv
  1404.3776, not opened): fewest triangles for a given maximum error is
  NP-hard, so the budget form, the same problem turned round, is a heuristic
  here as everywhere.
- **Shoham and Gersho 1988** (bit allocation; recalled, not re-read):
  splitting a budget across parts for a *summed* error needs equal marginal
  gain per part. It is why the stop is on the maximum: for the maximum, one
  tolerance for every piece is the best split (the basin's error is its worst
  piece's), which is what increment 23 already runs. RMS is reported, never
  the stop.
- **Searched and not found** (study, section 3): a greedy-insertion terrain
  method with a budget shared across independently refined pieces. One web
  search, no database search: **no novelty claim** is made, and none is
  needed; the search over T is elementary.

**Departures from the literature, named:** none that drop a guarantee. The
count is reached by a search over the tolerance rather than by stopping the
greedy mid-run (section 3 says why), so the result is "about N" (within a
window the caller sets) rather than exactly N; hmm's `-t` is exact on one
mesh. Nothing in `refine` changes.

### Legacy

`git grep -n -i "Count_stop\|lindstrom_turk\|count_ratio\|-ratio" legacy-archive -- legacy`
(run 2026-10-25) returns `legacy/bindings.cpp` lines 11-12 and 333-342
(CGAL's Lindstrom-Turk edge collapse bound as `lindstrom_turk_by_ratio` with
`Count_ratio_stop_predicate` and `lindstrom_turk_by_size` with
`Count_stop_predicate`), `legacy/rasputin/mesh.py` lines 119-123 (their
callers), and `-ratio` as the CLI's control in `application.py` line 41
(default 0.4), `geo_tiff_reader.py` line 29 and `web_visualize.py` lines 43
and 96. Intent: the old rasputin triangulated every DEM node and
*simplified down* to a share or a count of the triangles; a count was its
only control, and the tolerance stop is the new one. Nothing is carried
across: it is CGAL, runs the other way, and offers no error bound.

## 3. The rule

**Inputs.** `--triangles N` (an integer ≥ 1) and `--max-triangles M`
(an integer > N; default the larger of ⌊1.05 N⌋ and N + 1). The **window** is every count from
2N − M to M inclusive: as far under N as M is over it. Ola's "soft and hard":
N is the aim, M is never exceeded.

**The count that is judged is the written count**: after `refine`, the edge
strip or the final check against the source DEM, and the NoData trim. The
study measured the final check adding 0.3 % to 77 % to refine's count
depending on the input (its section 9.2), and on Kultsjön it adds about 7 %
(2,201 points at 10 m, table in section 7), more than a 5 % window; so a
trial runs the whole chain to the trimmed mesh, and nothing is written until
the chosen trial.

**The floor.** The outline, the feature lines and the start-quality points
make a mesh that no tolerance can shrink: the count at a huge tolerance.
`refine` needs a finite tolerance, so the floor trial runs at
T_FLOOR = 10⁶ m (section 5). Its largest error E_F is the tolerance at and
above which the mesh is the floor, so a floor result reports
`tolerance_m` = E_F, not 10⁶.

**The search** (`budget.search`, pure; section 4.1), given a `trial(T)` that
returns the written count c(T) and the largest error at T. **Every trial of
every step counts toward `MAX_TRIALS` = 12**, and two things end the search
early: a trial in the window, and a **saturated** trial, one whose largest
error is 0 (every valid DEM node is a vertex, so no smaller tolerance can
change the mesh; `refine` stops at tolerance 0 for the same reason,
`refine.hpp`'s header).

1. The floor: F, E_F = trial(T_FLOOR). If F > M, **refuse**: the message
   names F and E_F and says to ask for at least F (nothing is written). If
   2N − M ≤ F, the floor is the answer, reported at tolerance E_F. If E_F = 0
   the floor is saturated: it is the answer, and the record says the window
   was missed (the DEM has no more to give).
2. Coarse to fine, scale-free: T ← E_F / 4; trial. While c(T) < 1.5 F, T is
   not in the window, the trial is not saturated and fewer than twelve
   trials have run: T ← T / 4; trial. If still not in the window (same
   guards), one more trial at T / 4, so that two trials stand well above the
   floor. E_F > 0 here, so every T is positive, and c(T) > F for every
   T < E_F (the node with error E_F is inserted); an unchanged count between
   two quarterings is a plateau the loop steps through, bounded by the cap.
3. Until a trial lands in the window, the last trial is saturated, or twelve
   trials have run: through the
   last two trials (T_a, c_a) and (T_b, c_b), the secant on
   log(c − F) against log T:
   k = ln((c_a − F) / (c_b − F)) / ln(T_a / T_b), with each c − F floored at 1,
   and the next T = T_b · ((N − F) / (c_b − F))^(1/k). The bracket is the
   largest T whose count was over M (T too small) and the smallest T whose
   count was under 2N − M (T too large). When k ≥ 0 (equal counts included),
   T is not finite, or T
   is not strictly inside the bracket, the next T is instead the geometric
   mean of the bracket's two ends; with one end missing, the known end
   divided or multiplied by 4. Trial.
4. The answer is the first trial in the window. Otherwise (twelve trials, or
   a saturated trial outside the window) the answer is the trial with the
   largest count not over M, and the record says the window was missed;
   when that trial's mesh is not the last one kept, it is run again (the
   search keeps no mesh but a window hit's; the rerun writes the same mesh,
   by determinism, and its phases add to the Timings). The cap M holds by
   construction: nothing over it is ever chosen.

Run as a pure script before any test exists (`pure_search.py` in section
7's folder, output `pure-search.log`): on 360 synthetic count functions
(section 10, A6) every search landed in at most 7 trials; the saturated
floor ends in 1 trial, saturation reached mid-search in 3, a 3 % rise lands
in 5, and a gap wider than the window ends missed after exactly 12 with the
largest count under the cap chosen (A8, A9, A11's figures).

Why a search and not an exact stop inside `refine`: a round can nearly double
the count, so stopping mid-run needs a sorted last round in C++; the study
(section 9.4) measured the count moving about 0.1 % per 0.1 % of T, so a 2 %
window holds about twenty reachable counts and the search lands when the
window is reachable. The
written count also moves after `refine`, which an in-loop stop cannot see.
Measured on three inputs (section 7): 5 to 7 trials.

**Determinism.** The rule reads only counts, errors and double arithmetic,
`refine` is deterministic for any thread count, so the trial sequence and the
mesh are functions of the inputs alone; `--tolerance` at the reported
`tolerance_m` writes the same mesh (test A4).

**Monotonicity is observed, not assumed.** The study found the written count
falling at every step of 0.5 % on four inputs (its section 9.2) and the ladders
below agree; the bracket guard in step 3 is what makes a rise harmless (a
secant that leaves the bracket is replaced by a bisection), and test A8
feeds the rule a count function with a rise.

## 4. Blueprint: data flow and boundaries

```
cli.py mesh                         refine_job.py / cli.py            budget.py (pure)
--triangles N --max-triangles M --> Window(asked=N, cap=M) -----------> search(trial, window)
                                    trial(T) = _refined(Once(tile), ...,       |  counts, errors
                                               replace(how, tolerance=T,       |  in; the chosen
                                               checks=store), clock)           v  T out
                                    keeps the Refined of the trial in the window
                                    ---> record: tolerance_m = T found, triangles_asked,
                                         triangles_cap, budget_trials; summary line
104b:  _core.error_summary(view | store, final mesh) --> rms_error_m, error_nodes_measured
       (every budget run, on the chosen trial; a tolerance run with --rms)
```

### 4.1 `budget.py`: the search, pure (104a)

New module `src_python/tin_engine/budget.py`, no numpy, no `_core`, no typer:

- `Window` (frozen dataclass): `asked: int`, `cap: int`; `low` = 2·asked − cap;
  `holds(count)`. Construction refuses cap ≤ asked and asked < 1 (a
  `ValueError` the CLI turns into a `BadParameter`).
- `Trial` (frozen): `tolerance: float`, `triangles: int`, `max_error: float`.
- `search(trial: Callable[[float], Trial], window: Window, *, floor_tolerance: float = FLOOR_TOLERANCE, max_trials: int = MAX_TRIALS) -> Found`
  with `Found` (frozen): `chosen: Trial`, `trials: tuple[Trial, ...]`,
  `missed: bool`; a budget under the floor raises `UnderFloor(floor: Trial)`.
  Section 3's rule, nothing else; `pure_search.py` (section 7's folder) is
  the rule written out once, and `budget.py` is that function with types.

Testable with a synthetic `trial` (a count function of T), no kernel: tests
A6 to A9 and A11.

### 4.2 Driving the trials (104a, `cli.py` and `refine_job.py`)

- `mesh()` gains `--triangles` and `--max-triangles` (section 4.4 for the
  words). `--tolerance` and `--triangles` together are refused ("one of").
  Every check that today says "needs --tolerance" because it needs
  *refinement* (`--domain`, `--start-min-angle`, `--start-quality-gain`,
  `--no-constraint-feet`) accepts either and says "needs --tolerance or
  --triangles". `--tolerance-near` and `--tolerance-slope` keep "needs
  --tolerance" (section 9, left out). `--pieces` keeps "cutting needs
  --tolerance and --domain"; `--piece-budget` is accepted with `--triangles`
  as the budget the piece check below is made against (section 4.5).
- `How` gets the budget: `tolerance: float | None` and
  `budget: Window | None`, exactly one set. `_dem_mesh` with a budget runs
  `budget.search` with a `trial` closure that calls today's `_refined` once
  per T, keeps the `Refined` only when its count is in the window (so at most
  one trial's mesh is alive beside the tile), and reads the trimmed triangle
  count and `values["max_error_m"]`. The accepted `Refined` is what the
  tolerance path returns today; nothing downstream changes. On a missed
  window (section 3, step 4) the closure runs the chosen tolerance once more
  and keeps that mesh.
- **A trial that would need pieces refuses the run.** Today a `--tolerance`
  run with a domain is cut automatically when `partition(cols, rows,
  tolerance, pieces=1, piece_budget)` has more than one cell (`cli.py`, the
  `partition` call after `_open_dem`), because the bytes a node needs grow as
  the tolerance shrinks (`decompose.bytes_per_node`: 19 B above 50 m, 637 B
  at 0). A budget run cannot take that path (4.5), and it does not know its
  tolerances in advance, so **the closure checks `partition` at each trial's
  T before running it**: more than one cell is a `RunError` the CLI turns
  into exit 2, nothing written, with the message "--triangles N needs a
  tolerance of T m here, which would cut the window into nx × ny pieces
  under the piece budget of B bytes; the triangle budget does not run on
  pieces yet: use --tolerance T with --domain, or raise --piece-budget". The
  floor trial at T_FLOOR takes the fewest bytes, so the refusal, when it
  comes, comes after some trials have run; the record of a refused run is not
  written. Without a domain no run is cut today, and the check is not made.
  `check_canvas` guards the window's memory as on every uncut run.
- **The tile stays alive through the trials.** `_refined` takes
  `Once[DemTile]` and drops its reference after the strip or before the final
  check; the closure holds the tile and hands a fresh `Once(tile)` to each
  trial, so `_refined` is unchanged. The peak is the tile plus the check-point
  store plus one trial's mesh, where a tolerance run holds the tile and the
  store at different times; on Kultsjön the tile is 1,575 × 2,153 float32
  (13.6 MB) and the store 3.67 M points × 16 B (59 MB). The piece check
  above is what keeps a budget run's window at the size one piece may have
  under the piece budget.
- **The check-point store is built once.** `How.checks` is an `Iterator[Block]`
  of the source's projected nodes, consumed by `final_check.run`; it cannot
  feed a second trial. `final_check.run` is split: `stored(grid, checks,
  clock) -> CheckPoints` (today's filing and freeze, with its two clock rows)
  and `run(start, store, tolerance, ...)` taking the frozen store;
  `How.checks` becomes `Iterable[Block] | CheckPoints`, and `_refined` files
  blocks into a store when given blocks. The budget closure files once before
  the first trial. The store is read-only after `freeze()` (the header says
  any number of threads may read it), so trials share it.
- **Timings.** Each trial's phases add into today's rows (refine, edge strip,
  check points, final check, trim), so `--stats`'s table sums over the
  trials; the record's `budget_trials` says how many ran and at which
  tolerances, and the Timings note gets one sentence when a budget ran.

### 4.3 The error summary (104b, C++)

RMS is computed nowhere today (`grep -rn -i rms src_python/tin_engine
include` finds no such code). One new pass over the written mesh, once per
run, after the final check, **on every budget run (the chosen trial) and on
a tolerance run that asks with `--rms`** ("also report the root-mean-square
height error over the DEM nodes inside the mesh; one more pass over the
mesh"). Not on every run: the pass costs 3 to 6 % of the whole `rasputin
mesh` time on the quick check's cases with every core, and 22 % at one
thread (measured below), which is outside the quick check's noise band at
one thread, so a run that does not ask does not pay. RMS only: the mean
signed error is left out (section 9).

- `include/terrain/refinement/measure.hpp`: `struct ErrorSummary
  { double max_error; double sum_squares; std::size_t count; }` with
  `rms()`, and two overloads of
  `error_summary(..., unsigned threads)`: one over a `RasterSource` (the
  native path: the DEM's nodes, corner heights as `scan` takes them) and one
  over a `CheckPoints` store with the mesh's z (the reprojected path: the
  source's nodes, as `refine_points` reads them). Both build the lattice mesh
  with `detail::to_lattice` from the mesh as written (vertices in the node
  rectangle; the final outcome's arrays before the trim, as `refine_points`
  takes phase 1's), scan every triangle in parallel with
  `for_each_block`, one partial per slot, and reduce the partials
  serially in slot order, so the sums are bit-identical for any thread count.
  A void triangle (a NoData corner) measures nothing.
- **The walk is `scan`'s and `scan_points`'s, not a copy.** Both get a
  template parameter `class Acc = NoAccumulator`; the `consider` and `visit`
  lambdas call `acc.add(error, owned)` under `if constexpr`, so
  `refine`'s and `refine_points`'s instantiations compile to what they are
  today (bit-identical meshes and no new work in the hot loop; `@perf`
  checks both). The plane expression is therefore the same object `refine`
  evaluates, which is what makes `max_error` of the summary equal the
  recorded `max_error_m` (test B3).
- **Each node once.** A node on an edge shared by two triangles is in both
  closed node sets. The accumulator counts it only in the triangle that owns
  the edge: the lower triangle index of the two, or the only one on a border
  edge; in `scan` the orientation `o[k] == 0` (exact among nodes;
  `orient_sign` for an off-node corner) says which edge, in `scan_points`
  the `orient_sign` already computed. Vertices are never in a node set, and
  a check point coincident with a vertex is skipped as today (it is counted in
  `coincident`, with its own largest difference, which `refined_record`
  already folds into `max_error_m`). So `count` is the nodes strictly inside
  the mesh plus those on interior edges, each once; a NoData node is not a
  point.
- `bindings/core.cpp`: `error_summary(view_or_store, vertices, triangles, z,
  valid, *, threads=0) -> ErrorSummary` (two overloads by the first
  argument's type, as `refine_points` and `refine_strip` take their arrays),
  GIL released; `ErrorSummary` read-only with `max_error`, `rms_error`,
  `count`. `_core.pyi` mirrors it.
- `measure.summary(final, against, threads)` in Python (`refine_job.py`):
  the final outcome's arrays against a `RasterView` (native path, while the
  tile is held) or a `CheckPoints` store (reprojected path). A tolerance run
  with `--rms` calls it from `_refined` after the strip or the final check,
  before the tile reference goes (`How.rms: bool`); the budget closure calls
  it once, after `search` returns, on the chosen trial, with the tile and
  store it holds. The values become `rms_error_m` and
  `error_nodes_measured`; the `--stats` row is "error summary".

**Cost, measured** (`scan_cost.py` in section 7's folder, output
`scan-cost.json` and `scan-cost.log`; AC power, master `97759b7b`'s
`_core`, the quick check's cases from `docs/benchmarks/quick/cases.toml`
and its AC baseline's medians). One full pass is `refine(view, written
mesh, tolerance = 10⁶)` on the native path and `refine_points(store,
written mesh, 10⁶)` on the reprojected one: one round, nothing inserted,
`scan_seconds` read; its largest error equals the record's `max_error_m`
on every case, so the pass measures the population the summary will.
The round-1 premise that a pass is about 2 % of `refine`'s scan was false:
`refine` rescans only the triangles a round touched (`refine.hpp`'s
header), so its whole scan is only two to three full passes.

| case (threads) | run total (s) | `refine` scan total (s) | one full pass (s), median of 5 | pass / run |
|---|---:|---:|---:|---:|
| tile, 25.5 M nodes, 464 k triangles (all cores) | 0.4245 | 0.0568 | 0.0265 | 6.2 % |
| tile (1 thread) | 0.7255 | 0.3678 | 0.1567 | 21.6 % |
| quarter, 9.0 M nodes (all) | 0.3638 | 0.0463 | 0.0160 | 4.4 % |
| velhas, 72,859 check points, reprojected (all) | 0.1248 | final check 0.0074 | 0.0038 | 3.0 % |
| romsdal-slope, 25.5 M nodes, 3.10 M triangles (all) | 3.588 | 0.5649 | 0.1201 | 3.3 % |
| numedalslagen, features, 653 k triangles (all) | 8.846 | 0.3967 | 0.1144 | 1.3 % |
| lagan, features, 12.3 M check points, reprojected (all) | 14.50 | final check 0.2659 | 0.1643 | 1.1 % |

The probe's `refine` legalises the start mesh before its one scan; on the
four cases without feature lines it flipped nothing and its largest error
equals the record's, while on the two feature cases it flipped edges
(33,924 on numedalslagen, with the constraint edges marked or not:
`pass_constrained.py`, `pass-constrained.log`) and so scanned a mesh of the
same size and node population but not the written one (largest error
49.4 m against the record's 10.0 m). Their cost figures stand; their errors
say nothing about the written mesh. `error_summary` scans the mesh as
given and legalises nothing, which is also why no existing binding can
stand in for it.

Kultsjön (reprojected, 3.67 M check points): the final check's first round
scans every point, so one pass is under its whole scan of 0.071 s
(section 7's record): under 2 % of the 4.4 s run, and the budget run pays
it once. Section 11 gates both paths.

### 4.4 The record and the words

`run_record.WORDING` gains, in record order after `max_error_m`:

| key | `--stats` wording | value |
|---|---|---|
| `triangles_asked` | Triangles asked for (about this many) | N |
| `triangles_cap` | Triangles at most | M |
| `budget_trials` | Tolerances tried for the triangle budget, and the triangles each gave | "1e+06 m: 28561; 54 m: 29439; 13.5 m: 49263; 3.37 m: 291983; 9.79 m: 67494" |
| `budget_window_missed` | The triangle budget's window was missed; the nearest count under the cap was written | "yes", only when it was |
| `rms_error_m` (104b) | Root-mean-square height error | over `error_nodes_measured`; only when measured |
| `error_nodes_measured` (104b) | DEM nodes inside the mesh the height errors are over | |

`tolerance_m` keeps its wording; on a budget run it is the tolerance found.
None of the new keys goes into the mesh file's fields (`in_file` false): the
`.vtk`, `.ply` and GeoPackage bytes of a `--tolerance` run do not change,
which is what lets `@perf` check byte-identity (section 11). The stderr
summary (`run_record.summary`) says, on a budget run, before today's
sentence: "67494 triangles, about the 65893 asked for (at most 69187):
tolerance 9.79 m found in 5 trials." and, with 104b, today's sentence ends
"(largest difference 9.7885 m, root-mean-square 3.00 m)" when the RMS was
measured. On a missed window the first sentence instead says "67494
triangles, the most under the cap of 69187 (about 65893 asked for; no
tolerance gave a count in between): tolerance 9.79 m, 12 trials", and when
the chosen trial is saturated "640 triangles: every DEM node is a vertex,
the most this DEM gives (about 1000 asked for)". Plain words, no
labels (the plain-output rule of increment 25).

The MCP server derives its `mesh` tool from the CLI's Typer command
(`mcp/schema.py`, `derive`), so `--triangles` and `--max-triangles` appear
there without a change.

### 4.5 The pieces

A cut run (`--pieces`, `--piece-budget`, increment 23) refuses the budget in
104, as it refuses `--stats`, `--record` and GeoPackage output today ("not
written for a cut run yet"). Reason: pieces run one at a time under the
memory budget and each is dropped before the next is loaded
(`basin_run.cut_run`), so the basin's written count at a trial tolerance is
known only after every piece has run, and a search would re-run the basin
per trial. Ola's ruling (one tolerance for the whole basin) is honoured by
construction when the budget reaches the cut path: the study's step B finds
that tolerance by a search on a sample before the run, and checks the cap
after the pieces are joined; it is a later increment, and the cut path's
`How.tolerance` stays a float. The automatic cut by memory (`partition`
after `_open_dem` in `cli.py`) needs a tolerance to size memory; a budget
run makes that check per trial and refuses rather than cuts (section 4.2),
so a budget run is always one window, guarded by `check_canvas` as every
uncut run is. Kultsjön at 31 m is 1,575 × 2,153 = 3.39 M nodes, under
16 GiB at every tolerance (637 B × 3.39 M = 2.16 GB at tolerance 0), so no
trial of the demo can refuse (section 7).

## 5. Constants

| constant | value | scale it assumes | largest input checked |
|---|---|---|---|
| `FLOOR_TOLERANCE` | 10⁶ m | above any height difference on Earth (under 9,000 m) and any NoData sentinel's magnitude that could leak into an error (32,767) | floor at 10⁴ and 10⁷ m gave the same count on the study's four inputs; Kultsjön at 320 m and 10⁶ m: 28,561 both |
| default cap | max(⌊1.05 N⌋, N + 1) | the study's window; a 2 % cap costs 0 to 2 trials more (the study's section 9.3) | Kultsjön 5 %: 5 trials; 2 %: 6 (section 7) |
| coarse-to-fine start | E_F / 4, quartered | scale-free: the floor's own largest error | Kultsjön E_F = 216 m: 54, 13.5, 3.37 m; velhas E_F = 76.6 m: 19.1, 4.79 m |
| "well above the floor" | 1.5 F | the fit on log(c − F) needs c − F clear of the rounding of F | all three cases below |
| `MAX_TRIALS` | 12 | counts every trial of every step; nearly twice the most seen (7, measured and synthetic) | section 7; the 360 synthetic searches of section 3 (at most 7) |
| saturation | largest error = 0 | exact: a zero largest error means every valid node is on the mesh, so no tolerance changes it | bumpy without `--stride` (640 triangles, error 0 at every T, `search-bumpy.json`) |

## 6. Judged against the architect's criteria

- **Data against execution:** `Window` is data; `search` is a pure function
  of a `trial` callable; the CLI maps options to a `Window` and refusals, and
  the pipeline supplies the callable. No C++ in 104a.
- **State:** the search holds counts and errors, never a mesh; the pipeline
  holds one tile, one store and one accepted mesh. `error_summary` is pure
  and thread-safe (per-slot partials, serial reduce).
- **Dependency gravity:** none added.
- **Async-readiness:** the trials are the same `_refined` calls the worker
  runs today, each releasing the GIL inside `_core`; the search is a loop
  over them and can be driven from a thread as the tolerance run is.
- **Locality, parallelisation, numerics** (the geometry skill's rule): no
  geometric mutation changes; `refine` and `refine_points` are untouched at
  the instruction level; the only new numerics are the secant on logs in
  Python (double, finite by the guard) and sums of squares in double, reduced
  in a fixed order.

## 7. Measured (2026-10-25), and the predictions fixed before any run

Worktree `tri-budget`'s `.venv`, `_core` built 10:25 from master `97759b7b`,
battery (`pmset -g batt`: 100 %, discharging), 4 threads. Scripts and
outputs: `../rasputin_scratch/ideas/research/2026-10-25-triangle-budget-kultsjon/`
(`rms_probe.py`, `search_probe.py`, `ladder-t*.{json,log}`, `rms.json`,
`search-*.json`; the folder is beside the study's results folder and is not
in git).

### Kultsjön, whole catchment, as the delivery meshed its four sections

The command is the delivery's (`mesh-sections/kultsjon-002.json`, `command`)
with the whole catchment as the domain
(`~/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon.geojson`,
1,702.9 km², 868 vertices, EPSG:3006; the delivery cut it on a 50 km grid
into four sections of 394 + 46,988 + 5,197 + 33,958 = 86,537 triangles at
10 m; one window here gives fewer triangles because the seams and their
cut-polygon vertices are gone):

```
rasputin mesh --out kultsjon-t10.gpkg --dem glo30 --tolerance 10.0 \
  --domain .../catchments/kultsjon.geojson --domain-crs EPSG:3006 --out-crs EPSG:3006 \
  --features ~/.cache/rasputin/clc2018/features.gpkg --features-legend corine \
  --features-repair 1.0 --features-tolerance 50.0 --features-outline-snap 5.0 \
  --threads 4 --memory-budget 4G --z-unit m --stats kultsjon-t10.md --record kultsjon-t10.json
```

GLO-30 (the four cached tiles N64-N65, E014-E015, EPSG:4326) resampled to a
31 m grid of 1,575 × 2,153 in EPSG:3006; CORINE 2018 as constraints (16
features kept, 11,362 line vertices after joining); 3,670,772 source nodes
in the check-point store, of which 2,108,499 fall inside the mesh. Total
4.44 s with the GeoPackage and `--stats` (decode 1.61, features clip 1.38,
refine 0.047, check points project 0.27 + store 0.07, final check scan 0.071,
write 0.41); the `.vtk` runs of the ladder took 2.7 to 3.7 s wall.

The ladder (`--tolerance T`, written counts; errors from `rms_probe.py` over
the 2,108,499 source nodes inside the mesh, whose largest error agrees with
the record's `max_error_m` to 1 × 10⁻⁷ on every row, so the probe measures the
population the CLI reports):

| T (m) | triangles | refine inserted | final check inserted | largest error (m) | RMS (m) | mean signed (m) | wall (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10⁶ (floor) | 28,561 | 0 | 0 | 215.97 | 15.66 | −0.45 | 2.80 |
| 320 | 28,561 | 0 | 0 | 215.97 | 15.66 | −0.45 | 2.81 |
| 160 | 28,564 | 1 | 1 | 158.36 | | | 2.74 |
| 80 | 28,742 | 85 | 6 | 79.93 | 13.37 | 0.91 | 2.79 |
| 40 | 30,590 | 968 | 53 | 39.94 | 8.78 | 1.17 | 2.74 |
| 20 | 38,228 | 4,547 | 321 | 19.99 | 5.35 | 0.76 | 2.80 |
| 10.5 | 62,521 | 15,302 | 1,782 | 10.50 | 3.18 | 0.35 | 2.84 |
| **10** | **65,893** | 16,577 | 2,201 | **9.9998** | **3.05** | 0.33 | 4.44 (gpkg) |
| 9.5 | 69,892 | 18,258 | 2,530 | 9.50 | 2.93 | 0.32 | 2.88 |
| 5 | 165,394 | 56,588 | 12,144 | 5.00 | 1.69 | 0.14 | 3.05 |
| 2.5 | 447,064 | 158,169 | 51,797 | 2.50 | 0.93 | 0.05 | 3.68 |

Read: the count never rises with T; the floor is 43 % of the 10 m count,
so the raw count falls only as T^-1.1 around 10 m while c − F falls as
T^-2.0 (the floor-subtracted fit is what the search needs); the largest
error tracks T to 0.05 %; RMS is about 0.43 · T^0.85 between 5 and 20 m.
The mean signed column is the probe's measurement, kept as data for the
sup-norm bias idea (section 13); the product does not report it (section 9).

The prototype search (`search_probe.py`, section 3's rule run through the
whole CLI per trial) on this input:

| N | M | trials (T in m → triangles) | answer |
|---:|---:|---|---|
| 65,893 | 69,187 (5 %) | 10⁶ → 28,561; 53.99 → 29,439; 13.50 → 49,263; 3.375 → 291,983; **9.789 → 67,494** | T 9.7885 m, 67,494 triangles, largest error 9.7885 m, 5 trials |
| 65,893 | 67,210 (2 %) | the same four, then 9.789 → 67,494 (over M), **10.020 → 65,764** | T 10.020 m, 65,764, 6 trials |
| 30,000 | 31,500 | 10⁶ → 28,561 | the floor is in the window: T = 215.97 m, 1 trial |
| 200,000 | 210,000 | 10⁶; 53.99; 13.50; 3.375 → 291,983; **4.265 → 208,315** | T 4.2647 m, 208,315, 5 trials |

### Predictions for the demo, fixed here

P1. `--triangles 65893` (the 10 m mesh's count, default cap 69,187) writes
**67,494 triangles at tolerance 9.79 m in 5 trials**: 2.4 % more
triangles than asked, within the 5 % cap, with the largest error 9.79 m
that this count buys on the curve (the 10 m mesh has 65,893 triangles and
9.9998 m; neither mesh is better, they are two points on one curve). RMS
**3.00 m** (interpolated on the ladder between 9.5 and 10 m; the 10 m
mesh's is 3.05 m, the same curve). Tolerance for the implementation's
rounding of the secant: T within 9.78-9.80 m, count within 67,400-67,600,
RMS within 2.98-3.02 m. With `--max-triangles 67210` (2 %): 65,764
triangles at 10.02 m, largest error 10.02 m, RMS about 3.055 m, 6 trials.

P2. Error against count, both meshes on one curve: the ladder's columns
above are the curve, and the budget run's `budget_trials` adds five points
on it (the trials are `--tolerance` runs). The budget mesh and the tolerance
mesh lie on the same curve by construction; they differ only in T.

P3. Time: a `--tolerance 10` run with GeoPackage and `--stats` takes 4.4 s;
the budget run adds four trials of refine + strip + final-check scan + trim
(about 0.15 s each here; the decode, the features clip, the check-point
projection and the write are done once) and one error summary pass (under
0.08 s, section 4.3), so **under 5.5 s**; the summary
line and `--stats` say 5 trials. The prototype's 12.7 s is the whole CLI
re-run per trial and is not the product's time.

P4. `--triangles 30000` (cap 31,500, window from 28,500) returns the floor
(28,561 triangles, tolerance reported 215.97 m) in one trial;
`--triangles 27000` (cap ⌊1.05 × 27,000⌋ = 28,350, under the floor) is
refused with the floor's count and error in the message, exit 2, nothing
written. (`--triangles 28000` would not refuse: its cap 29,400 holds the
floor.)

### The suite's fixtures (figures for section 10)

`tests/fixtures/velhas` (ANADEM extract, 239 × 288 at 30 m, EPSG:4674,
`--domain catchment.geojson --out-crs EPSG:31983`; the reprojected path, as
Kultsjön): floor 49 triangles, E_F = 76.598 m; written counts 64 m: 54;
16 m: 204; 8 m: 520; 4 m: 1,320; 2 m: 3,178; 1 m: 7,386; 0.5 m: 14,734;
0.25 m: 24,571. Errors over the 12,681 source nodes inside the mesh: 1 m:
largest 0.99789, RMS 0.42122, mean −0.01554; 0.5 m: 0.49997, 0.21482,
−0.01158; 2 m: 1.99671, 0.77732, −0.05450; floor: 76.598, 21.243, +10.037.
Prototype searches: N 7,386 (5 %) → 5 trials (10⁶; 19.149; 4.7874; 1.4224;
**1.03647 → 7,106**); N 14,734 (5 %) → 6 trials, T 0.50386, 14,655;
N 24,571 (5 %) → 6 trials, T 0.27020, 23,388; N 24,571 (2 %) → 7 trials,
T 0.25180, 24,454; N 48 with M 50 → the floor (49) in one trial.

`bumpy` (`cli_helpers.rough_dem(14)`, 17 × 21 nodes of uniform noise 0-50 m,
no domain) with `--stride 8` (the native path, with the edge strip): floor 12
triangles, E_F = 36.938 m; 10 m: 450; 5 m: 546; 2 m: 597; 1 m: 615; 0.5 m:
631; 0 m: 640 (every node). Without `--stride` the start is every node and
every tolerance gives 640 triangles with error 0 (`search-bumpy.json`): the
saturated floor, test A11.

The velhas target grid at the catchment domain is 158 × 130 = 20,540 nodes
(`open_dem` on the fixture, 2026-10-25); `partition` with
`--piece-budget 2054000` (100 B per node) gives one cell at the floor
(19 B per node), at 19.149 m (27 B) and at 4.787 m (72 B), and 2 × 2
cells at 1.4224 m (220 B): the fourth trial of the N 7,386 search refuses
(test A12). The default budget of 16 GiB never cuts it.

The pure rule on synthetic count functions (`pure_search.py`,
`pure-search.log`, 2026-10-25): c(T) = F for T ≥ E_F, else
F + round(A · ((E_F / T)^k − 1)), largest error min(T, E_F); F in
{0, 50, 28,561} with E_F 50 m for F = 0 and 216 m otherwise; k in {0.6,
1.1, 2.0}; A = 3 × 10⁵ / (1000^k − 1) (so c(E_F / 1000) − F is 3 × 10⁵);
20 targets N log-spaced from max(c(E_F / 8), 1.3 F + 1, 20) to
c(E_F / 2000), each with a 5 % and a 2 % cap: 360 searches, every one in
the window, **at most 7 trials** (F = 28,561 and k = 2.0: 38 of 40 take 7;
F = 0: 3 to 6), the second trial E_F / 4 in every search. The same rule on
the three measured ladders, interpolated in log-log on c − F, reproduces
section 7's trial counts within one (Kultsjön N 65,893: 5 and 6 trials;
velhas N 7,386: 5 and 6; N 24,571 at 2 %: 7).

## 8. The showcase (Wednesday 2026-10-28)

Deliverables, made on Tuesday from the merged build by whoever the main
session assigns (`@perf` after its acceptance run is the natural owner; the
commands are section 7's), into a folder of Ola's choosing under
`~/rasputin-output/` (default
`~/rasputin-output/kultsjon-budget-2026-10-27/`):

1. `kultsjon-t10.gpkg` and `kultsjon-n65893.gpkg` for QGIS, with their
   `--stats` Markdown and `--record` JSON; `kultsjon-t10.vtk` and
   `kultsjon-n65893.vtk` for ParaView (a second run each; the runs are
   deterministic).
2. `error-vs-triangles.png`: largest error and RMS against the triangle count
   on log axes, the ladder as the curve, the two meshes marked, the five
   trials as hollow points; made by a scratch script from the records (not
   product code; matplotlib is not a dependency and stays out of the
   package). Three series at most, legend in plain words. The caption says
   what the plot shows: every mesh rasputin writes, by tolerance or by
   budget, is a point on this one curve; the budget picks the point by its
   count, the tolerance by its error.
3. A `README.txt` with the two commands and the summary lines.

The demonstration itself: run `--tolerance 10` and `--triangles 65893` live
(each under 6 s), show the two summary lines, open the plot. Said plainly:
"ask for a size, get at most that many triangles, the tolerance it took is
found and reported, and the mesh is the one `--tolerance` writes at that
tolerance". Not said: that the budget mesh is better; it is not.

## 9. Left out

- The budget on a cut run (section 4.5): the study's step B, a later
  increment.
- The budget with `--tolerance-near` or `--tolerance-slope`: the study
  proposed scaling the whole policy by one factor; it is a semantics question
  (is N metres on the lines still N?) Ola has not ruled on, and the showcase
  does not need it. They keep "needs --tolerance".
- Exact N (the study's step C): dropped by the study (section 9.4).
- The mean signed error and percentiles of the error (the ideas-file entry
  *Mean, RMS and percentiles of the vertical error in --stats*): RMS is
  taken up by 104b; the mean signed error was dropped on review (the
  showcase does not need it, and the bias question has its own idea entry);
  percentiles need a histogram or a sort and are left.
- RMS on every tolerance run without asking: measured at 3 to 6 % of the run
  with every core and 22 % at one thread (section 4.3), so it is `--rms`.
- Warm starts between trials: the study found a warm step slower than a cold
  run and its mesh not the one `--tolerance T` gives (section 9.4).

## 10. The red suite

Fixture figures are section 7's. `A` tests are 104a's, `B` tests 104b's.
Invariant-critical, with mutation targets: A1-A4, A9, A11, A12 and B1-B3
(the window held, the cap on every trial, saturation, the piece check, the
mesh equal to the tolerance run's, each node once).

**A. The budget (104a)**

- A1 `test_cli_triangles.py`: `--triangles 7386` on velhas writes a count in
  [7,017, 7,755], the record's `tolerance_m` is 1.03647 ± 0.0005,
  `max_error_m` ≤ `tolerance_m`, `triangles_asked` 7386, `triangles_cap`
  7755, `budget_trials` names five tolerances; the stderr summary says
  "about the 7386 asked for (at most 7755)" and "5 trials". Mutation
  targets: the window's lower end (2N − M), the cap test, the "first in the
  window" choice.
- A2 the cap: `--max-triangles 7533` (2 %) lands in [7,239, 7,533]; `--max-triangles 7386`
  and `--max-triangles 7000` are refused naming both options; `--max-triangles`
  without `--triangles` is refused.
- A3 the floor: velhas `--triangles 48` writes 49 triangles with
  `tolerance_m` 76.598 ± 0.001 and one trial; `--triangles 40` is exit 2,
  nothing written, the message names 49 triangles and 76.6 m and "at least
  49".
- A4 same mesh as the tolerance run: after A1, `--tolerance <tolerance_m
  from the record>` on velhas writes the same vertices, z and triangles
  (arrays equal; the file fields equal but for the budget keys, which are not
  in the file). Mutation target: any change to `_refined`'s arguments per
  trial.
- A5 the native path: bumpy `--stride 8 --triangles 546` lands in
  [519, 573]; `--triangles 10` is refused naming 12 triangles.
- A6 `test_budget.py` (pure, synthetic `trial`): section 7's synthetic grid
  exactly (c(T) = F for T ≥ E_F, else F + round(A · ((E_F / T)^k − 1)),
  largest error min(T, E_F); F in {0, 50, 28,561}, E_F 50 m for F = 0 and
  216 m otherwise, k in {0.6, 1.1, 2.0}, A = 3 × 10⁵ / (1000^k − 1), 20
  targets log-spaced from max(c(E_F / 8), 1.3 F + 1, 20) to c(E_F / 2000),
  caps 5 % and 2 %): all 360 searches land in the window in **at most 7
  trials** (measured, `pure-search.log`), the first trial is
  `floor_tolerance`, the second E_F / 4, and the sequence is identical on a
  second call.
- A7 refusals and the floor in `search` itself: `UnderFloor` carries the floor
  trial; a window holding the floor returns it after one trial with
  `tolerance` = E_F.
- A8 a rise: c(T) as A6 (F 28,561, E_F 216, k 2.0, A 80.3) with one step of
  +3 % below T = 3 m; for N in {40,000, 65,893, 150,000, 300,000, 600,000}
  at 5 % the search lands in 5 trials (measured), and in every case
  `chosen.triangles ≤ cap`.
- A9 the twelve-trial cap: A8's function without the rise and with every
  count in [60,000, 72,000] mapped to 72,001 (a gap wider than the window),
  N 65,893 and M 69,187: ends `missed` after exactly 12 trials with the
  largest count ≤ M chosen (measured: 59,949). Mutation targets: a trial
  not counted toward the cap; the choice of the largest count under the
  cap.
- A10 options: `--triangles` with `--tolerance` refused ("one of");
  `--triangles` with `--pieces 2` refused with today's cutting message;
  `--triangles` with `--piece-budget` passes the option checks (A12 then refuses at the fourth trial);
  `--domain` with `--triangles` and no `--tolerance` accepted; the MCP
  `mesh` tool's schema lists `triangles` and `max_triangles`
  (`test_mcp_cli.py`'s derivation helper).
- A11 saturation: in `search`, a `trial` returning largest error 0 at every
  T with count 640: window (639, 670) returns it after one trial with
  `tolerance` 0, not missed; window (1000, 1050) ends after one trial,
  `missed`, chosen 640. A function saturating mid-search (A6's form with F 12, E_F 36.938, k 0.9, A 60, as `pure_search.py`;
  F 12, E_F 36.938, counts capped at 640 with error 0 there) and window
  (700, 735): ends `missed` after 3 trials (measured), chosen 640. On the
  CLI: bumpy without `--stride`, `--triangles 639` writes 640 triangles,
  `tolerance_m` 0, one trial; `--triangles 1000` writes 640, the record's
  `budget_window_missed` is "yes", the summary line says every DEM node is
  a vertex. Mutation target: the saturation test (error = 0) in each step.
- A12 pieces refused: velhas `--triangles 7386 --piece-budget 2054000` is
  exit 2 with nothing written; the message names 1.42 m, 2 × 2 pieces,
  2054000 bytes and `--tolerance`; the same with the default piece budget
  lands (A1). Mutation target: the check made after the trial instead of
  before it (the test counts the trials run with a spy on `_refined`: three, not four).

**B. The error summary (104b)**

- B1 `test_core_error_summary.py`, a plane: a 3 × 4 `elevations()`-style
  tile with z = 0.5 x + 0.25 y + 3 and the two-triangle mesh of its corners:
  `max_error` 0, `rms_error` 0, `count` exactly 8 (the twelve nodes less
  the four vertices, those on the diagonal counted once), the same for
  `threads` 1 and 4.
- B2 each node once, by construction: a 5 × 5 tile of uniform noise (seeded),
  the four-triangle fan from the centre: `count` equals the number of nodes
  minus the five vertices (every node lies inside or on an interior edge);
  a `CheckPoints` store with points placed exactly on the shared edges gives
  the same `count` as the number of points minus coincident ones. Mutation
  target: the ownership rule (lower triangle index).
- B3 agreement with the record: velhas `--tolerance 1 --rms` has
  `rms_error_m` 0.42122 ± 0.0001, `error_nodes_measured` 12,681, and
  `max_error_m` equal to the summary's
  largest error (the test reads both from the record and recomputes the RMS
  with an independent barycentric evaluation over the fixture's projected
  nodes, the oracle using its own plane, never the record). `--tolerance 2`:
  0.77732; `--tolerance 0.5`: 0.21482. The equality is pinned on a run
  without feature lines; on a run with them the edge strip's rescan skips
  nodes within its coincidence radius of a strip point (`scan.hpp`'s
  header), so the summary, which skips nothing, may report a larger value
  than `max_error_m` there. On Kultsjön the probe found them equal to
  1 × 10⁻⁷ on every ladder row (section 7); `@perf` records both on the
  demo run (section 11).
- B4 thread independence: velhas `--rms --threads 1` and `--threads 4`
  records have a bit-equal `rms_error_m`.
- B5 the words: the summary line ends with "root-mean-square 0.42 m)" on
  velhas at 1 m with `--rms`, and without `--rms` is today's line and the
  record has no `rms_error_m`; `--stats` rows for the two new keys and the
  "error summary" timing row; the `.vtk` file
  fields of a `--tolerance` run are unchanged from master's (the golden test
  `test_refine_golden.py` passes unchanged).
- B6 the budget path measures the chosen trial without asking: A1's record
  has an `rms_error_m` and the budget run on velhas reports
  `error_nodes_measured` 12,681; its summary line has "root-mean-square".

## 11. Acceptance (`@perf`)

- 104a touches no C++: the 1 m benchmark's mesh hashes equal the baseline's
  (`docs/benchmarks/2026-10-19/23c-4`, the last run with hashes), timings
  within noise. The budget's own timing on Kultsjön, section 7's command
  with `--triangles 65893` against `--tolerance 10`, both to `.gpkg` with
  `--stats`: the budget run is at most **1.5 s slower** (four more trials at
  about 0.15 s each is 0.6 s, the error summary under 0.08 s, noise on a
  4.4 s run about 0.3 s, doubled), and its `--stats` Timings show the
  trials' refine rows summed. Power state recorded; the predictions P1-P4
  checked and written beside the figures.
- 104b touches `scan.hpp` and `refine_points.hpp`: the 1 m benchmark's mesh
  hashes unchanged (the accumulator is a template default; the instantiations
  `refine` uses must not change), and **without `--rms` the whole
  `rasputin mesh` time within noise** on the quick check (`tools/bench_quick.py`;
  its band is the larger of 5 % and the baseline's spread, and 50 ms): no
  new work runs. **With `--rms`**, the "error summary" row of `--stats` is
  the measured pass of section 4.3 within a factor 1.5 (tile all cores
  0.027 s, one thread 0.157 s; romsdal-slope 0.120 s; velhas 0.004 s), and
  the run total grows by that row and no more (noise aside). The budget run
  on Kultsjön pays one pass: under 0.08 s. The
  thread-scaling sweep as the protocol says.

## 12. Size, plan to Tuesday, and the cut line

| part | files | counted lines |
|---|---|---|
| 104a | `budget.py` (new) about 100 (the saturation and cap guards add about 10); `cli.py` two options, the refusals, the trial closure with the piece check, the summary about 85; `final_check.py` the store split about 15; `refine_job.py` `How` about 3; `run_record.py` four keys and the summary sentences about 15 | **about 220** |
| 104b | `measure.hpp` (new) about 80; `scan.hpp` and `refine_points.hpp` accumulators about 35; `bindings/core.cpp` about 40; `_core.pyi` about 20; `cli.py` `--rms`, the call and two values about 20; `refine_job.py` `summary` and `How.rms` about 15; `run_record.py` two keys and the sentence about 8 | **about 220** |

Prose in the same PRs (not counted): `README.md`'s refinement bullet gets
`--triangles N` and `--max-triangles M` (104a) and `--rms` (104b);
`docs/GLOSSARY.md`, *Refinement*:
*triangle budget*, *floor* and *saturated*; `project_structure.md` rows for
`budget.py` and `measure.hpp`; `ROADMAP.md` row 104.

Plan: Sunday evening design review; Monday 104a red, green, review, `@perf`
timing, merge; Monday afternoon 104b red and green (C++), Tuesday morning
`@perf`'s acceptance and review, merge; Tuesday afternoon the showcase
outputs (section 8). **If Monday slips, 104b is cut**: the showcase then
shows the largest error from rasputin and the RMS from the scratch probe
(`rms_probe.py`, labelled as such on the plot), and 104b follows after
Wednesday. 104a is not cut into smaller parts: without the window, the
floor and the same-mesh guarantee there is nothing to show.

## 13. Ideas-file entries taken up or overlapped

- *A triangle budget: spend N triangles where they cut the TIN error most*
  (Ola, 2026-10-17): taken up, one window; the cut path later.
- *Mean, RMS and percentiles of the vertical error in --stats*
  (2026-10-08): RMS and mean taken up by 104b; percentiles left.
- *Unbiased catchment height statistics under the sup-norm* (Ola,
  2026-10-08): overlapped, not taken up; the probe's measurement of the
  bias on Kultsjön at 10 m (+0.33 m, section 7) is a figure for that idea's
  research, and `error_summary`'s accumulator is where a signed sum would go.

## 14. Questions for Ola

1. The default cap is 5 % over N (`--max-triangles` sets another). Default
   taken: 5 %, the study's window; 2 % costs about one trial more.
2. The demo budget is the 10 m mesh's own count, 65,893, so the two meshes
   compare at equal size; the budget mesh then lands at 9.79 m with 2.4 %
   more triangles. Default taken: 65,893; a round number such as 50,000 or
   100,000 is one more run.
3. The root-mean-square error is reported on every budget run and on a
   tolerance run only with `--rms`, because one pass costs 3 to 6 % of the
   run with every core and 22 % at one thread (section 4.3). Default taken:
   `--rms` opt-in; making it unconditional is a one-line change later.

Defaults taken on review's two open questions (round 1): no mean signed
error (104b is RMS only); RMS on every tolerance run only if its measured
cost were within noise, which it is not (above), hence `--rms`.

## 15. As built (104a)

Red `c36cdbe1`, green `d3bf0f04`, test fix `8f992e23`; then the commits
in items 7 to 10, to `8e55deaf`. 247 net production lines against about 220
(`python3 tools/count_loc.py 47308a83 8e55deaf`). What `@developer`
pinned beyond the design, and whether the design accepts it (`@architect`,
2026-10-25):

1. **`How.tolerance` stays a float**, not `float | None`: on a budget run it
   holds `FLOOR_TOLERANCE` until each trial replaces it
   (`replace(how, tolerance=T)`). Accepted. In `_dem_mesh` the value is read
   only as "refining or not" before the trials, each trial replaces it before
   `_refined` reads it, and the cut path refuses the budget (4.5). The design's "exactly one set"
   is lost as a type; nothing reads the placeholder as a tolerance.
2. **The kept trial is the one with the most triangles not over the cap**,
   so a missed window needs no rerun. Accepted, and it is section 3's
   choice, not a new one: every trial before a window hit is over M or
   under 2N − M, so the hit is the most under the cap; on a miss it is
   step 4's choice, and ties keep the first, as `search` does. The rerun
   and its Timings addition (step 4) are gone. Cost: two trial meshes can be
   alive during a trial (4.2 said one); the kept one has at most M
   triangles, so the extra memory is bounded by the cap asked for.
3. **A floor answer reports the floor's largest error** as its tolerance,
   saturated or not (0 when every node is a vertex). Accepted: section 3,
   step 1, and A11.
4. **`--triangles` has `min=1`**: 0 and −5 get click's range error before
   `Window` sees them. Accepted: the same refusal (4.1), in click's words.
5. **`--piece-budget` without a domain is accepted and unused.** Accepted:
   4.2, "without a domain no run is cut today, and the check is not made".
6. **A trial's largest error is the record's**: `max_error_m` raised to the
   difference at DEM nodes on vertices. Accepted: it is the error the
   record reports, and saturation then means 0 there too.
7. **The Timings sentence for a budget run**: missing at `8f992e23`, sent
   back, and in at `656f5247` (red `13ddebfb`). On a budget run the
   `--stats` Timings note says "one run per trial" and adds "The triangle
   budget ran N trials: the refine, edge strip, final check and trim rows
   sum all of them." The check point store rows are left out: the store is
   filed once per run. Accepted: it is 4.2's sentence.
8. **Two extras from Ola beyond the design** (2026-10-25). On a budget run
   the tolerance found and the largest difference print to three
   significant figures, in the summary's "within ... m" sentence and in the
   warning (red `e862b0cc`, green `632cfce5`; red `0eb20d02`, green
   `67564863`); a `--tolerance` run keeps its digits. The refusal of a
   budget that the floor already exceeds is headed "Invalid value for
   --triangles", not "--dem" (`0eb20d02`, `67564863`); the refusal of a
   budget run that would need pieces followed, as code review round 1
   suggested (red `7ea9ea53`, green `8e55deaf`). Both refusals keep their
   messages.
9. **`@perf` ACCEPTED at `9c7f47aa`**, run on `632cfce5`
   (`docs/benchmarks/2026-10-25/104a/README.md`). Every `--tolerance` mesh
   is byte-identical to master's and its timings within noise. On Kultsjön
   `--triangles 65893` found 9.7885 m in 5 trials, 67,494 triangles, and
   costs 0.82 s over `--tolerance 10` (`--stats` total 3.714 s against
   2.896 s; 1.5 s allowed); `--tolerance 9.788516166348431` writes the same
   `.vtk`. P1, P2, P3 and P4 met; P1's RMS waits for 104b. The two wordings
   that run reported are item 8's.
10. **Mutant kill round at `956a261d`** (`@tester`, code review round 1's
   blocking item): 35 mutants of the named targets A1-A4, A9, A11 and A12,
   planted one at a time in a scratch copy. 33 killed, 4 of them only after
   tests were strengthened in that commit; 1 equivalent (a count equal to
   the cap ends the search before the bracket is read); 1 survivor outside
   the named targets (a trial's largest error without the error at DEM
   nodes on vertices), which needs a fixture whose error there exceeds the
   final check's; none has one today.

## Review

- **Design review, round 1, 2026-10-25.** `@reviewer`, head `0f16382557a8af041b12f790b2dd6b55e15aa15e`, range `97759b7b..0f163825`. Verdict: CHANGES REQUESTED. LOC: 0 (design); about 190 (104a) and 220 (104b) estimated. Blocking: (1) section 3 step 2's quartering loop has no trial cap and loops at floor error 0 (`search-bumpy.json`) or when the DEM's most triangles are under 1.5 times the floor: "Every trial counts toward the twelve; when T is 0 or the count stops changing, the search ends missed", with a test; (2) P4: `--triangles 28000` gives cap 29,400, which holds the floor of 28,561, so it returns the floor; use `--triangles 27000` (cap 28,350); (3) A6's "at most 9 trials" was never run and the synthetic floor trial's largest error and targets are unstated; a simulation took 10-11 trials and some searches missed the window; state both and run the pure rule before `@tester`; (4) 104b's cost premise is false: `refine.hpp:19` rescans only changed triangles; on the quick check's tile (0.425 s) a +3 % gate allows 13 ms; time one full scan pass on the quick check's cases; (5) section 4.2's "basin-scale windows are the cut path's, which refuses the budget" is false: a budget run skips the automatic cut by memory (`cli.py:996`); check `partition` at each trial's tolerance and refuse a budget run that would need pieces. Showcase: can claim size asked, at most M triangles, tolerance found and reported, mesh identical to `--tolerance` at that tolerance; cannot claim a better mesh (9.79 m comes with 2.4 % more triangles, on the same curve); the error-against-count curve is the demonstration; rewrite P1's "(−2.1 %)" and the ROADMAP's "RMS 3.00 m against 3.05 m" so neither reads as an improvement. Cut-run refusal acceptable for Kultsjön; Tuesday plan sound. Not pushed; no CI.
- **Design review, round 2, 2026-10-25.** `@reviewer`, head `a6c00bf9ba24ebebfbe4850b27cab639dcd31aaa`, range `d1d9942f..a6c00bf9`. Verdict: APPROVED WITH WORDINGS. LOC: 0 (design); about 220 (104a) and 220 (104b) estimated. Round 1 blockers closed: trial cap and zero-error stop hold; P4 at 27000 (cap 28,350 under the floor); `pure-search.log` matches A6, A8, A9, A11; `scan-cost.log` timings match the quick check baseline; A12 `partition` figures recomputed (one piece down to 4.787 m, 2 x 2 at 1.4224 m); showcase wording no longer reads as an improvement. Five wordings applied in this commit: default cap max(floor(1.05 N), N + 1) (lines 103, 433); "nearly twice the most seen" (436); `--piece-budget` passes the option checks (703); A6 form named (710); trial count by a spy on `_refined` (721). Ideas: P3 "under 5.5 s" is tight; no flip count logged for lagan; `partition` ValueError handling during a trial unstated. Not pushed; no CI.
- **Code review, round 1, 2026-10-25.** `@reviewer`, head `67564863dce370be18c00f5e0ef3456e237275a9`, range `47308a83..67564863`. Verdict: CHANGES REQUESTED. LOC: 246 net (estimate about 220, ceiling 700). Blocking: (1) no mutant kill record for the invariant-critical tests A1-A4, A9, A11, A12; `@tester` runs the kill round for every named target; (2) ROADMAP row 104 and section 15 item 7 are stale after 656f5247 (Timings sentence in), 9c7f47aa (`@perf` ACCEPTED) and the later wording commits; section 15 omits the Timings, rounding and wording commits, Ola's two extras (three significant figures on a budget run; the floor refusal naming `--triangles`) and the acceptance run. Sound: no red-step scaffolding; section 15 items 1-6 match the code; README and glossary match; the acceptance at 632cfce5 still stands; citations gate exits 0. Suggestion: the refusal of a budget run that would need pieces still reads "Invalid value for --dem". Not pushed; no CI.
- **Code review, round 2, 2026-10-25.** `@reviewer`, head `d316fa269b64c1e5405f5026cda10cb692c4bd8d`, range `67564863..d316fa26`. Verdict: APPROVED WITH WORDINGS. LOC: 247 net (this round +1, `cli.py`). Round 1 closed: the kill record in `956a261d` covers every named target in section 10 (A1-A4, A9, A11, A12), 35 mutants; survivors "bracket-cap>=" (equivalent) and "sat-cli-at" (outside the named targets) acceptable; new tests match the record (433 passed); Status line, section 15 items 7-10 and ROADMAP row match. The pieces refusal raises `typer.BadParameter` headed `--triangles` (`cli.py:2477-2483`), no longer relabelled at `cli.py:1122`. `@perf` acceptance `9c7f47aa` stands. One wording applied in this commit (ROADMAP row 104: the floor refusal named, the pieces refusal "followed on review"); Status and row set to round 2 approved. CI after push.
