# Increment 15e: four memory fixes on the reprojected path

Status: **designed by `@architect`, 2026-10-03; implemented.** Red 9879805
(`@tester`), green efb854f (`@developer`), fixes after `@reviewer` round 1
in 18b7a79 and 6d2f40c. `@reviewer` round 2 approved (78df916).
`@perf`'s acceptance run: **ACCEPTED** (2026-10-03,
`docs/benchmarks/2026-10-03/15e-acceptance.md`; summary under "Acceptance
run" below). Pending: a push with Ola's yes, and a green CI. Departures from
the design are under "As built" below.

Ola's order of 2026-10-03 for the basin: the nine level-3 meshes as the
interim result, then this PR, then 23b-23g. Source: fixes 1-4 of
`docs/research/basin-memory-options.md` §2 (a). Fix 5 (P1's third copy while
loading) and the cache-backed, windowed `SourceWindows` are **not** in scope.

Named 15e because all four fixes sit on 15c-2's reprojected path
(`--out-crs`); 15d names the window-decoding design that 23a-1 replaced, and the edge strip is 15f.

Evidence labels are the note's: **[m]** measured, with its record; **[d]**
derived by arithmetic from measured figures; *est.* from the code alone.
Records: `docs/benchmarks/2026-10-02/basin-phases/` (sub-basin 761, per
phase) and `docs/benchmarks/2026-10-02/basin-memory-probe/` (one block).
Line numbers are those of master at `7810cf8`.

## Prior art: legacy and literature

*Literature.* Nothing here is new. Fix 1 is strip-mining: you bound the
working set by processing a fixed amount of data per pass. Fix 4 is a slab
(arena) allocator, which carves fixed-size chunks out of a few large
allocations (Bonwick, "The Slab Allocator", USENIX Summer 1994). It is
combined with per-row chunk lists, the usual answer when bucket sizes are
not known in advance. Fixes 2 and 3 only remove copies and references. **No
novelty claim is made**, so no search was needed.

*Legacy.* Nothing to carry across. The legacy tree meshes a projected DEM
through CGAL. It has no resampling, no check-point store and no
memory-bounded blocking:

```
$ git grep -l -i -E "resampl|shrink_to_fit|\.reserve\(" legacy-archive -- legacy
legacy-archive:legacy/rasputin/triangulate_dem.h
$ git grep -n -i "resampl" legacy-archive -- legacy
(no output)
```

The one hit is ten `reserve` calls on CGAL output vectors
(lines 297 to 793 of the archived header above, at the
`legacy-archive` tag), sized from a known count. They have nothing
to do with any fix here.

## Fix 1: resample blocks sized by node count

**Site.** `src_python/tin_engine/target_grid.py@44fa7f5:155`
(`block_rows: int = 256`), used at `:166` and `:187`.

**Change.** Add a module constant `BLOCK_NODES = 1 << 20` and a pure helper
`rows_per_block(cols: int) -> int`, which returns `max(1, BLOCK_NODES // cols)`.
`resample`'s keyword becomes `block_rows: int | None = None`. `None` means
`rows_per_block(grid.cols)`; an explicit value still overrides it, because
the J3 test (`tests/python/test_target_grid.py:248`) passes one. Values do
not change: J3 says each node depends on its own coordinates alone.

**Why 2^20 (default; Ola may change it).** One block's transient is
133-135 B per node [m, 761 P2 rise]. That makes ~140 MB per block, and
~1.4 GB with 10 threads [d]. The ~10 GB transient on 761 today [d] comes
from 256 × 27,786 = 7.1 M nodes per block. 761 gets 37 rows per block
(524 blocks, against 76 today), and the basin's 40,943 columns get 25. Each
block builds one transformer (`:167`), so there are more transformer
constructions. @perf times this.

**Red test, @tester.**
- `rows_per_block`: a function of `cols` alone, with `rows × cols <=
  BLOCK_NODES`; one row once `cols > BLOCK_NODES`.
- `resample` uses it. Spy on `tg.reprojector`: wrap it so that each
  transformer it returns records the length of the points it is given.
  Then call `resample` with no `block_rows`, on a grid whose
  `cols * 256 > BLOCK_NODES`. Lower `BLOCK_NODES` with `monkeypatch` so the
  grid stays small. Assert every recorded length is `<= BLOCK_NODES` and that
  the number of calls is `ceil(rows / rows_per_block(cols))`.

  Red today: there is no helper, and blocks are 256 rows.

## Fix 2: adopt the canvas instead of copying it

**Site.** `src_python/tin_engine/target_grid.py@44fa7f5:205`,
`return DemTile(meta=meta, array=canvas)`. This is the public constructor's
copy: 2.15 GB on 761 [m, the "P2 end" row].

**Change.** `return DemTile._adopt(meta, canvas)`. `np.empty` at `:163`
already gives a C-contiguous array of shape `(rows, cols)`, which is what
`_adopt` checks. `resample` never hands the canvas out writable, which is
R7's condition. These texts change with it:
- `DemTile._adopt`'s docstring (`src_python/tin_engine/io/models.py@44fa7f5:114-121`): two callers,
  `assemble` and `resample`;
- the comment at `src_python/tin_engine/io/cog.py@44fa7f5:160-161`, which still says M15 pins `assemble`
  as the one caller. `cog`'s own copy stays: it belongs to fix 5;
- 15c D3's "*Built otherwise (15c-2)*" note
  (`docs/increments/15c-geographic-dem.md:394-396`) becomes "as designed".

**Red test, @tester.**
- M15's grep test
  (`test_adopt_is_called_only_from_the_mosaic_and_the_resampler`, formerly
  `tests/python/test_mosaic.py:1620-1626`) now expects
  `["io/models.py", "mosaic.py", "target_grid.py"]`. Red today. This test
  change is the one the note anticipates. It lands in @tester's red commit.
- Identity: `resample(...)`'s `tile.array.base is None` (the canvas
  itself). Today the array is a view of a copy, so the check fails. Also,
  by M15's spy pattern, `DemTile._adopt` is called once and its argument
  `shares_memory` with `tile.array`.

## Fix 3: drop the target tile before phase 2 (D1)

**Site.** Three references keep the target tile alive through
`final_check.run`. 15c D6 records this as "*Built otherwise (15c-2)*"
(`15c-geographic-dem.md:587-590`). The references are:
- `mesh()` holds `opened`, a frozen `DemInput` (`src_python/tin_engine/cli.py@92e5356:829`). It reads
  `opened.tile` at `:839`, `:846` and `:864-865`.
- The call `_dem_mesh(opened.tile, ...)` (`src_python/tin_engine/cli.py@92e5356:845-859`). Its argument
  lives on the caller's frame until the call returns, so a `del` inside the
  callee frees nothing.
- `_dem_mesh`'s parameter `tile` (`src_python/tin_engine/cli.py@92e5356:1398`), live until `:1491`.

**Change** (CLI only; `open_dem` and `DemInput` keep their public shape):
- A private one-shot holder in `cli.py`, `_Once[T]`, whose `take()` returns
  the value and forgets it. A second `take()` is a programming error
  (`AssertionError`).
- `mesh()` reads into locals what it needs from `opened` later: `meta`,
  `plan`, `seams`, `grid`, `source_crs`, `domain`, `label`, `checks`. It
  wraps the tile in `_Once` and runs `del opened` before calling
  `_dem_mesh`. `:864-865` read `dem_run.meta`, which `_DemMesh` already
  returns.
- `_dem_mesh(held: _Once[DemTile], ..., grid, checks)` takes `grid` and
  `checks` in place of `opened`. It calls `tile = held.take()` first, and
  on the tolerance path it runs `del tile` after phase 1's `refine(...)` and
  before `final_check.run`. `meta` (small) stays.
- The `RasterView` made by `to_core(tile)` is a temporary argument of
  `refine`. It is gone when `refine` returns: `bindings/core.cpp` has no
  `keep_alive`, so the outcome does not hold it. The check-point iterator
  holds the *source* mosaic, not the target tile, and that is unchanged.
- 15c D6's "*Built otherwise*" note (`:589-590`) becomes "as designed".

**Red test, @tester.** This is in-process, through Typer's `CliRunner`, on
an existing small `--out-crs --tolerance` fixture from
`tests/python/test_cli_mesh_geographic.py`.
- Wrap `cli.open_dem` so that it records `weakref.ref(result.tile.array)`.
- Wrap `cli.final_check.run` so that on entry it runs `gc.collect()` and
  records whether that ref is dead, then calls the original.
- Assert the run exits 0 and the record is `[True]`.

Red today: `_dem_mesh` still holds the tile. A NumPy array takes a weakref;
whether a Pydantic model does was not checked, which is why the probe uses
the array.

## Fix 4: the check-point store in one arena (C++)

**Site.** `include/terrain/refinement/check_points.hpp`. Today it keeps
`std::vector<std::vector<Packed>> rows_` (`:130`), grown by `push_back`
(`:73`), and `freeze()` ends with `shrink_to_fit` (`:88`). Live data is
16 B per point, but 30 B per point is charged [m, `probe_store.py`, one
run]: the doubling reallocations fragment malloc's small zone, and
`shrink_to_fit` raised the footprint (2.80 → 3.04 GB [m]).

**Choice (default; the alternative is below).** Fixed chunks from large
slabs, with no reallocation and no per-row `malloc`:
- `kChunk = 1024` points (16 KiB), `kSlab = 4096` chunks (64 MiB, one large
  allocation: `std::make_unique_for_overwrite<Packed[]>`, so pages are
  charged only as they are written).
- Each row has a directory of chunk pointers and a count. `add` appends to
  the row's last chunk and takes a new chunk from the current slab when the
  last one is full. Slabs are freed only with the store.
- `freeze()`: for each row, gather into one reusable scratch vector (the
  size of the largest row), sort and unique as today, scatter back into the
  same chunks, and set the count. No `shrink_to_fit`.
- `for_each_in`: the same `lower_bound` on column, over the logical index
  `i → chunk[i >> 10][i & 1023]`. The order, filing and duplicate rules are
  unchanged, as is thread-safety after freeze.
- One new observer, `reserved_points()` (chunks handed out × `kChunk`), for
  the test. The Python binding does not change.

Expected: 16 B per point plus at most one partly filled chunk per row. On
761 that is ≤ 19,377 × 16 KiB = 0.32 GB, ≤ 1.2 B per point, so
**~16.5-17 B per point *est.*, against 30 [m]**: about 3.5 GB less on 761's
270.4 M points [d].

*Alternative, not chosen:* a counting pass, then one exactly sized block.
`checks` is a single-pass iterator, so counting means projecting every
check point twice. That adds the "check points: project" phase again,
8.6-9.5 s of 48 s on 761 at 10 m [m, `runs/761_t10*.stats.md`], about 20 %
of wall time, to save ~1 B per point.

**Red test, @tester** (Catch2, next to
`tests/cpp/unit/test_refinement_check_points.cpp`):
- Bound. Add N points over R rows, unevenly: one row with
  `3 * kChunk + 1` points, others with 1. Then
  `reserved_points() <= N + R * kChunk`, and the bound still holds after
  `freeze()`. A per-row vector cannot give this observer.
- The existing CP1 cases stay green unchanged. They pin order, filing,
  duplicates and frozen behaviour.
- Add one case across chunk boundaries: a row of `2 * kChunk + 7` points,
  added in reverse order, with duplicates on both sides of a boundary. It
  checks the frozen order, `duplicates()`, and `for_each_in` column ranges
  that start and end inside different chunks.

  **This is the increment's invariant-critical suite.** If any suite gets
  mutation testing, it is this one (README, cost constraints), aimed at the
  gather/scatter and the index arithmetic.
- Red: `reserved_points` does not exist, so the new cases do not compile.
  Put them in their own test source (one target or file) so that the rest
  of `ctest` still builds and runs in the red commit.

## Expected effect on 761's peak

| phase | today, `MallocLargeCache=0` | after 1-4 |
|---|---:|---|
| P1 load | 6.81 GB [m] | unchanged (fix 5 is out of scope) |
| P2 resample | 13.97-14.10 GB [m] | mosaic 2.24 + canvas 2.15 + ~1.4 transients ≈ **5.8 GB** [d] |
| P2 end | 6.63 GB [m] | 4.4 GB [d] (no copy) |
| P4-P5 store + phase 2 | under-counted with the variable | the tile no longer held (−2.15 GB [m]); store ~4.5 GB instead of ~8.1 [d, est. 16.5 B against 30 B per point] |

So the peak moves from P2 (~14.1 GB) to P1's load (~6.8 GB [m]), as the
note says. These are live-memory figures. **As shipped** (no variable),
macOS keeps the freed P1 copies charged (+4.3 GB [m]) but may reuse them
for later large allocations. The as-shipped peak cannot be derived from
the code, so @perf measures it. Mesh-bound peaks (769 at 2 m, 10.76 GB [m])
are not addressed here.

## Acceptance by @perf: applies

The diff touches `include/terrain/refinement/` and what drives the mesh
input, so the README rule applies:
- `tools/bench.py`: the 1 m benchmark and the thread sweep against the
  previous merge, power state recorded. The projected 1 m path does not use
  the store or `resample`, so no change is expected. The run checks that.
- **Memory:** re-run `docs/benchmarks/2026-10-02/basin-phases/run_phases.py`
  on sub-basin 761 at 10 m, both as shipped and with `MallocLargeCache=0`.
  Compare the per-phase peaks with the table in that README (17.83 / 13.97
  GB `time -l` peak; P2 end 10.19 / 6.63; P4 11.51 as shipped). For the
  store, run `probe_store.py` again (30 B per point today), because the
  footprint under-counts the store when the variable is set.
- **Time:** "resample" (more blocks, more transformers) and
  "check points: store" (gather/scatter in freeze) against
  `runs/761_t10*.stats.md`.
- Evidence goes in `docs/benchmarks/<date>/`.

Afterwards, re-measuring `b(T)` for 23's piece budget follows (the note,
§2 (b)). That is not part of this PR.

## LOC estimate

Production lines, by `CLAUDE.md` §2's count:

| fix | lines |
|---|---:|
| 1, `target_grid.py` | ~6 |
| 2, `target_grid.py` (+ comments, which do not count) | ~1 |
| 3, `cli.py` (`_Once`, `mesh()` locals, `_dem_mesh` signature) | ~20 |
| 4, `check_points.hpp` (arena, freeze, `for_each_in`, observer), net | ~40 |
| **total** | **~65-70** (+39 % margin: ~95) |

The note's "~40 lines" under-counts fix 4: it priced a `reserve`, not an
arena.

## Defaults taken, and questions for Ola

None of these blocks the work. Each is a default Ola may override:
- *Default:* `BLOCK_NODES = 2^20` (fix 1).
- *Default:* the slab arena, not the counting pass (fix 4). The trade is
  ~1 B per point against ~9 s per 761 run [m].
- *Default:* fix 3 stays in the CLI. `open_dem`'s `DemInput` keeps
  holding the tile for Python API callers. They drop it themselves.

No question needs Ola before @tester starts.

## As built

Where the code departs from the design above (the design text is left as
it was ruled):
- **`CheckPoints` is move-only.** Its copy constructor and copy assignment
  are deleted, because each row's chunk directory holds raw pointers into
  the store's slabs; a copy would point into the original's slabs. The
  defaulted moves keep those pointers valid: the slabs are
  `std::unique_ptr<Packed[]>` and move with the store without being
  reallocated. The binding still returns the store by value from its
  constructor (`bindings/core.cpp`, the `py::init` lambda), which compiles
  against the move.
- **`for_each_in` uses a hand-written lower bound** over the chunked logical
  index, in place of `std::ranges::lower_bound` (which would need a
  random-access iterator written over the chunk list). It finds the same first position, the first
  point with `col >= c0`, then walks forward while `col <= c1`: the same
  order and the same inclusive `[c0, c1]` rule as before. The index is
  written `chunks[i / kChunk][i % kChunk]`, equal to the design's
  `i >> 10`, `i & 1023` for `kChunk = 1024`.
- **`kSlab = 4096` is a public `static constexpr`**, beside `kChunk = 1024`,
  rather than an implementation detail.
- **`_Once.take` raises `RuntimeError` on a second take**, not
  `AssertionError`: it is an explicit check, because `python -O` strips an
  `assert`.

Not a departure: `_dem_mesh` takes `grid` and `checks` in place of
`opened`, as fix 3 specifies.

## Review

### Round 1 (`@reviewer`, 2026-10-03): CHANGES REQUESTED

Range: `7810cf8..efb854f` (design 92e5356 + d81bbcd, red 9879805, green
efb854f). Not pushed, so **no CI exists**; local runs are not CI. The branch
sits on 7810cf8; master 390b516 merges cleanly (`git merge-tree`).

LOC by `CLAUDE.md` §2, production only: **58 net** (`check_points.hpp`
+54/-14, `cli.py` +32/-18, `target_grid.py` +8/-4; `io/models.py` and
`io/cog.py` change only docstrings and comments). 94 lines added gross.
`@developer` counted 59. The estimate was ~65-70, so no split is needed;
the ceiling is 700.

Local evidence: `cmake --build build` exit 0, `ctest` 847/847; `pytest`
3667 passed, 16 skipped, after rebuilding and copying `_core`; mypy, ruff,
format, prohibited-deps and detria gates OK. The four store suites
(arena, CP1, refine_points unit and property) pass under ASan+UBSan.
Fix 3's test fails when `del tile` is replaced by `pass`, and passes again
once it is restored.

Blocking:
1. Red-step scaffolding written in the present tense, now false (it is
   `@tester`'s to fix, because `@developer` does not edit tests):
   `tests/cpp/unit/test_refinement_check_points_arena.cpp:22-24`,
   `tests/python/test_target_grid.py@44fa7f5:399-400` and `tests/python/test_target_grid.py@44fa7f5:514-515`, and
   `tests/python/test_cli_mesh_geographic.py@3e01580:327-328`. Put them in the past
   tense ("went red at 9879805 because ...") or delete them.
2. This file's status line still says "Not started", and none of the
   as-built departures is recorded. Update the status line and add an
   "As built" note covering three things: `CheckPoints` is move-only
   (its copy operations are deleted, because rows hold raw pointers into
   the slabs); `for_each_in` uses a hand-written lower bound over the
   logical index instead of `std::ranges::lower_bound`; and `kSlab` is
   public. Splitting `_dem_mesh`'s `opened` into `grid` and `checks` is
   not a departure: fix 3 specifies it.

Merge also needs a green CI, `@perf`'s acceptance run, and the ROADMAP row.

### Round 2 (`@reviewer`, 2026-10-03): APPROVED

Range: `3e01580..44b72e8` (18b7a79 `@tester`, 6d2f40c `@developer`,
44b72e8 `@architect`). Still not pushed: no PR and no remote branch, so
**no CI exists**, and this approval does not make the branch merge-ready.
Master 390b516 still merges cleanly (`git merge-tree`).

LOC by `CLAUDE.md` §2, production only: the round adds **+1 net**
(`cli.py`: the `assert` becomes `if ...: raise RuntimeError`; the
`check_points.hpp` change is comment-only), so the branch is **59 net**
against the ~65-70 estimate and the 700 ceiling. No split.

Round-1 items:
1. Closed. The four notes are in the past tense ("went red at 9879805
   because ..."), and the `tests/cpp/CMakeLists.txt` comment no longer
   speaks of the red commit. No line the branch adds under `tests/`,
   `include/`, `src_python/` or `bindings/` still describes a red step in
   the present tense (grep of `git diff 7810cf8..HEAD`).
2. Closed. Status line, "As built" and the ROADMAP row added; each
   as-built claim checked against the code (`check_points.hpp:56-64`,
   `:123-127`, `:182`; `bindings/core.cpp`'s `py::init` lambda returns
   `CheckPoints{...}` by value). The 15e citations into files the branch
   edits were re-read as quotations at 7810cf8 and hold; the fix-2 test is
   now cited by name.

Local evidence: mypy, ruff, format, prohibited-deps and detria gates OK;
`test_cli_mesh_geographic.py`, `test_target_grid.py`, `test_mosaic.py`:
209 passed; the arena target rebuilt (exit 0) and passes. No full C++
rebuild: the round's C++ changes are comments.

Not blocking: the `_Once.take` second-take `RuntimeError` has no test. Its
one caller takes once, so nothing reaches it today; a two-line test would
pin it if `_Once` gains callers.

Merge still needs `@perf`'s acceptance run, a push with Ola's yes, and a
green CI.

### Pre-push check (after `@perf`, f288504)

- Citations: every `check_citations.py` at-risk entry outside this file was
  re-read as a quotation against `390b516` and this branch. None is moved by
  this branch. The `cli.py` hunks start at line 49 (a one-line replacement,
  no shift) and line 107; every cited line at or below 106 reads the same on
  both trees. Already stale on master and not this PR's: the `cli.py:NN`
  citations in `16e-multi-features.md`, `test_features.py` (`:39,42`, `:79`,
  `:84`), `tests/python/test_viz_svg.py@390b516:665` and `tests/python/test_cli_mesh_multi_features.py@97eea35:17`.
  `05b-noder-driver.md:1749` and `15-dem-mosaic.md:616` hold.
- f288504 holds: the status line, the ROADMAP row and "Acceptance run" match
  `docs/benchmarks/2026-10-03/15e-acceptance.md`, whose figures match
  `15e-acceptance/tables.md`, `runs/probe_store_isolated_*.out`, the `pmset`
  records (AC, charged) and `bench.py`'s ACCEPTED verdicts. Its claim that
  `7810cf8` and `390b516` differ only outside the code holds
  (`git diff --stat 7810cf8 390b516 -- include src src_python bindings
  tools/bench.py CMakeLists.txt pyproject.toml` is empty). The acceptance run
  is done. Remaining: the push, with Ola's yes, and a green CI.

## Acceptance run

`@perf`, 2026-10-03, AC power, branch at 78df916 against its base 7810cf8,
back to back: **ACCEPTED**. Evidence and method:
`docs/benchmarks/2026-10-03/15e-acceptance.md`.
- `tools/bench.py` 1 m run and thread sweep: no change (the second pair is
  within ±2.6 %), and the meshes and quality are identical.
- Sub-basin 761 at 10 m, `MallocLargeCache=0`: the peak falls from
  14.3-14.5 GB to 6.9 GB. Resample measures 5.3-5.6 GB (sampled) against
  ~5.8 GB derived. The overall peak is set in the store fill, 0.1 GB above
  the load's 6.81 GB. As shipped, the peak falls from 15.4-15.9 GB to
  10.1 GB. All eight meshes are byte-identical.
- The store is charged 16.4 B per point, against 16.5-17 B expected and
  34-36 B for the base in the same probe.
- Time: resample is 28-33 % faster. The freeze is +0.26 to +0.39 s
  (+7 to +11 %) as shipped and unchanged with the variable. The whole run is
  13-16 % faster.

**Round 3 (post-approval delta), 2026-10-03.** Range `70f8403..87df16c`. Verdict: APPROVED. LOC: 0 production lines (one line at `15e-memory-fixes.md:17`, 1 insertion, 1 deletion). The 15d half is true per `23-basin-scale.md` (23a-1 replaces 15d); the 15f half is true on branch `worktree-agent-a0cb49bea16623b07`; naming 15f before it merges is fine (a label, no path or line). No citation into 15e moved. Not pushed; no CI.
