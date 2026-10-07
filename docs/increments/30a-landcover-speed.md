# Increment 30a — the `land cover` phase of `rasputin mesh`, made fast

Status: **approved by code review round 3; Ola said yes to push, open the PR
and enqueue it (2026-10-06).** Round 1's test fix is
`3cc526f`. Designed by `@architect`
2026-10-06 on branch `worktree-landcover-speed`, on `@perf`'s profile commit
`a7154ec` (master `8199f30` plus the profile); design review round 1
APPROVED. Red `6cde0eb`, green `3066d60` (+12 net production lines,
`python3 tools/count_loc.py a7154ec 3066d60`), `@perf`'s acceptance
ACCEPTED at `f38ff1a` (section 9). Master has since gained #195 and #196; neither
touches `landcover.py`, its tests, `ROADMAP.md` or this phase's path through
`cli.py` (`git diff --stat a7154ec origin/master -- src_python/tin_engine/landcover.py
tests/python/test_landcover.py tests/python/landcover_fixtures.py ROADMAP.md`
prints nothing), so the branch was not merged with master for this design.

**What this is.** The first of three pull requests that remove the
bottlenecks `@perf` measured in `rasputin mesh` on two real catchments
(`docs/benchmarks/2026-10-06/bottlenecks/README.md`). Ola approved them on
2026-10-06 in this order: land cover (this, 30a), then the CORINE clip (30b),
then reading the DEM (30c). The ROADMAP row is 30.

**What changes for a user.** Nothing in any output: every triangle's
`land_cover_code` and the four numbers on the `land cover:` stderr line stay
exactly as they are. The phase gets about four times faster: on battery,
2.90 s to 0.71 s on Numedalslågen and 6.71 s to 1.80 s on Skiensvassdraget,
measured on a prototype (section 8). As built, on AC power: 2.85 s to 0.67 s
and 6.54 s to 1.74 s (section 9).

**Lean.** No mutation round (Ola's standing rule for lean rounds). No C++.

## 1. Prior art: legacy and literature

### Literature

The method is increment 16c's and does not change
(`docs/increments/16c-landcover-labels.md`, *R1. The label: components first,
one point per component*): connected components of the triangles across
edges that are not constraints, the spread of Shewchuk's Triangle `-A`
(Shewchuk, LNCS 1148, 1996), one point per component, a point-in-polygon
lookup. What changes is how three steps are computed:

- **Components** stay the array union-find 16c cites (Tarjan, *J. ACM*
  22(2), 1975; the hook-and-shortcut form of Shiloach and Vishkin, *J.
  Algorithms* 3(1), 1982). Only the sort that finds each edge's two
  triangles changes. The union-find's result does not depend on the order in
  which the joined pairs are listed (section 4.1), so the stable sort can
  become the default one.
- **The best triangle per component** becomes a grouped maximum (a segmented
  reduction, in the sense of Blelloch, *Vector Models for Data-Parallel
  Computing*, MIT Press, 1990: one reduction per group, no sort of the whole
  array) instead of a sort of every triangle. Recalled, not reread.
- **Point in polygon.** Both before and after, GEOS answers the
  `intersects` predicate between a point and a polygon. Before, the polygons
  were the tree and each *point* was the query, so GEOS prepared the point
  (which gains nothing) and walked every vertex of the unprepared polygon for
  each candidate pair: 920 M and 1.38 G vertex visits on the two catchments.
  After, each *polygon* is the query, so GEOS prepares the polygon (an index
  over its edges, JTS/GEOS `PreparedPolygon` with an indexed point-in-area
  locator) and answers each point against the index. The tree is Sort-Tile-
  Recursive packing (Leutenegger, Lopez and Edgington, ICDE 1997), as before.
  That a predicate-filtered `STRtree.query` prepares the *query* geometry and
  frees the prepared form afterwards was read in shapely 2.1.2's
  `src/strtree.c` (`evaluate_predicate`: `GEOSPrepare_r` on the query
  geometry, `GEOSPreparedGeom_destroy_r` after) for this design.

16c's literature section said the lookup ran "over prepared polygons". At
master it did not; this PR corrects that line in 16c, since a
documentation defect found during an increment is fixed in its PR.

**Novelty: none claimed.** This is standard engineering on a standard method;
nothing was searched because nothing is claimed.

### Legacy

```sh
$ git grep -l -i "land_cover\|prepared\|STRtree" legacy-archive -- legacy
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/globcov_repository.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/rasputin/land_cover_repository.py
legacy-archive:legacy/rasputin/tin_repository.py
legacy-archive:legacy/rasputin/web_visualize.py
legacy-archive:legacy/rasputin/wfs_repository.py
legacy-archive:legacy/tests/test_gml_repository.py
legacy-archive:legacy/tests/test_land_cover_repository.py
```

Read: `legacy/rasputin/gml_repository.py:184-226` (`land_cover`), the only
lookup. It tests each point against each polygon in a Python loop, unprepared,
with `# TODO: Move to C++ for speed!`; it clipped the polygons to a buffered
domain first. **Nothing is carried.** 16c already took the per-point test as
its test oracle. Clipping the polygons is not carried either: it changes a
polygon's area, and area decides overlaps (16c's Default D2), so it could
change a code (section 10).

## 2. What is slow, and why

From `@perf`'s profile (section 3 of the bottlenecks README; battery, scratch
copies of the code, identical output):

| step | today, s (Numedalslågen / Skiensvassdraget) | why |
|---|---|---|
| `_lookup` | 1.42 / 2.16 | each candidate pair walks every vertex of an unprepared polygon; one CORINE polygon has 370,070 vertices and is a candidate for many points |
| best triangle per component | 0.31 / 1.06 | `np.lexsort` over all 1.3 M / 3.0 M triangles to pick one per component (1,860 / 2,881 components) |
| `regions` | 1.01 / 3.07 | `np.isin` of 3.9 M / 9.1 M side keys (a sort inside), then a stable `argsort` of what is left |

Inputs: 1,287,334 / 3,029,295 triangles; 156,919 / 275,572 constraint edges;
1,611 / 2,628 coded polygons (all valid MultiPolygons, checked for this design
on the captured inputs with `shapely.is_valid`).

16c's R2 set a trigger: "If the phase costs more than the refine it follows
at 1 m, moving `regions` into the core is the next step". At 10 m on these
catchments it does today (2.87 s against refine's 1.16 s on Numedalslågen;
6.69 against 3.15 on Skiensvassdraget). After this PR it does not (section 8).

## 3. The blueprint: what changes where

One file of production code, `src_python/tin_engine/landcover.py`. No change
to `cli.py`, to the `label_triangles` signature or to `CoverLabels`; no C++,
no binding, no new dependency (numpy and shapely, as today).

```
cli._land_cover ──> label_triangles(vertices, triangles, edges, polygons=, margin=)   unchanged signature
                      ├─ regions(tri, edges)            sort once, then drop constrained sides with _member
                      │     └─ _member(values, table)    NEW: membership in a sorted, unique int64 table
                      ├─ incentres and r                 unchanged
                      ├─ _best(ids, r)                   NEW: per component, largest r, ties to lowest index
                      └─ _lookup(points, polygons)       the polygons become the query; a tree of points
```

### 3.1 `regions` (today `src_python/tin_engine/landcover.py@a7154ec:48-53`)

Sort all side keys once with the default `np.argsort`; the owner of sorted
position `i` is `order[i] % T` (the sides are three blocks of `T`, one per
triangle side, so the tile array is not needed); then drop the sides whose
key is a constraint key, by `_member` against `np.unique` of the constraint
keys. The pairing (`keys[1:] == keys[:-1]`) and the union-find loop are
unchanged. Return value unchanged: each triangle's component as its smallest
triangle index.

### 3.2 `_member(values, table) -> bool array` (new)

`True` where a value is in `table`, which is sorted and unique (int64 both).
`np.searchsorted` into the table, clipped to its last index, then an equality
test. An empty table gives all `False` (the clip would otherwise index an
empty array: the case of a mesh with no constraint edges, which
`TestRegions.test_no_constraints_at_all` already exercises through
`regions`). For `mypy --strict`, the return is wrapped so it is not `Any`
(the prototype's one `mypy` finding, `no-any-return`).

### 3.3 `_best(ids, r) -> int64 array` (new; replaces `src_python/tin_engine/landcover.py@a7154ec:95-97`)

Inputs: `ids` as `regions` returns them (so `ids[i] == i` exactly when `i` is
the smallest index of its component), and `r` finite. Output: one triangle
per component, in increasing component id, the one with the largest `r`, ties
to the lowest triangle index. The same answer as today's `lexsort` for finite
`r`. Steps: `np.maximum.at` of `r` into a per-id array that starts at `-inf`;
the triangles whose `r` equals their component's maximum; `np.minimum.at` of
their indices into a per-id array that starts at `T`; read it at the roots,
`np.flatnonzero(ids == np.arange(T))`, which are in increasing order and are
exactly what `np.unique(ids)` returns today.

`r` is finite because the trimmed mesh's x and y are finite (DEM lattice
nodes, input polygon vertices and refinement's points; the z of a trimmed
vertex is never NaN either, but only x and y enter `r`). A NaN `r` is outside
the contract: today's `lexsort` sorts it last, `np.maximum.at` would spread
it. The design pins this as a precondition rather than spending lines on it.

### 3.4 `_lookup` (today `src_python/tin_engine/landcover.py@a7154ec:124-125`)

Build the `STRtree` over the points and query it with the polygons, predicate
`intersects`; the result's two rows swap (`which, found` instead of `found,
which`). The ranking by area then code, the hit count and the first-hit rule
below it are unchanged.

**Why not `@perf`'s first variant** (query the polygon tree by box, then
`shapely.prepare` the candidate polygons and test): `shapely.prepare` stores
the prepared form *on the caller's geometry objects*. `label_triangles` would
then change the `FeatureSet` it was handed (more memory held, and a hidden
side effect in a function 16c calls pure). The reverse query is as fast (0.084
against 0.088 s on Numedalslågen in `raw/lc_numedalslagen.json`), and GEOS
frees its prepared form after each query. Red test P2 (section 7) pins that
the caller's polygons come back unprepared.

## 4. Why the output cannot change

### 4.1 `regions`

- *Which sides are dropped.* `_member(keys, unique(cut keys))` is the same set
  test as `np.isin(keys, cut keys)`, for any int64 keys.
- *Which triangles are joined.* After the constrained sides are dropped, two
  equal neighbouring keys join their two owners. Every key's sides sit in one
  run in the sorted order, whatever the sort; a run of two gives one pair, a
  run of three (an edge with three triangles, which a trimmed mesh does not
  have) gives two pairs that join all three. The set of components is
  therefore the same under any sort.
- *Which label each component gets.* The union-find hooks the larger root
  onto the smaller and stops at a fixed point. A label only falls and is
  always a member of the component, so the component's smallest index is
  never rehooked and ends as the root of all its members. The label is the
  smallest index whatever the order of the pairs. 16c's tests already pin
  this through permuted triangles (`test_ids_are_the_smallest_triangle_index_of_each_component`,
  `test_a_strip_of_1000_triangles_is_one_component`).

### 4.2 `_best`

The maximum per component is one of that component's `r` values, so `r ==
top[ids]` finds exactly the triangles at the maximum (no rounding: equality
of a value with itself), and the minimum index among them is today's
tie-break. Roots in increasing order equal `np.unique(ids)`, the order
today's `best` is in.

### 4.3 `_lookup`

`intersects` is symmetric, so querying points with polygons gives the same
(point, polygon) pairs as querying polygons with points; the order of the
pairs differs, and the code after it sorts by point, then area, then code,
and counts hits, neither of which depends on the input order (16c's D2
tie-break is to the smaller code, so two equal-area polygons with equal codes
give the same code whichever comes first). What GEOS computes differs:
before, a prepared point against a polygon; after, a prepared polygon against
a point. For valid polygons both are the exact point-in-area test (boundary
counts as inside). Checked for this design on hand-made points on a shared
edge, on a shared corner, on a slanted edge, on a hole's boundary and inside a
hole: identical codes and hit counts, base against prototype. Invalid polygons
are the risk (section 10); the captured CORINE sets have none.

## 5. The blueprint against the assessment framework

- *Data and execution:* no configuration involved; `margin` stays a
  parameter from the CLI.
- *State:* pure before and after. The change removes a side effect a
  tempting variant would add (section 3.4).
- *Dependency gravity:* none added.
- *Async-readiness:* unchanged; still one synchronous numpy and GEOS call.

## 6. The byte-identical gate

`docs/increments/30a-probes/landcover_bytes.py`, with its base output
`docs/increments/30a-probes/base_a7154ec.txt`, both committed with this
design. Its docstring says how to run it. Two modes are the gate:

- **`fixtures`** runs the three suites that reach the land-cover code
  (`test_landcover.py`, `test_cli_mesh_landcover.py`,
  `test_cli_mesh_multi_features.py`) in-process and records every call of
  `label_triangles` (from the module or from the CLI) and every direct call of
  `regions`: a hash of the codes and the four counts, keyed by test id and
  call number. 66 calls at the base.
- **`mesh`** runs `rasputin mesh` on Numedalslågen and Skiensvassdraget as
  `@perf`'s profile did and records a hash of the codes, the four counts and a
  hash of the whole `.vtk` written.

The gate: **the base lines are unchanged; the new test lines are not
compared.** Every line of the base file that starts `fixture ` or `mesh `
appears unchanged in the branch's output. A `fixture` line from a test this
increment adds has no base line and is not compared (at acceptance there were
six, from `TestCallerPolygons` and `TestDegenerateMeshes`).
A third mode, `replay`, reruns `label_triangles` on inputs `mesh --save` kept,
for a quick check while developing; it prints the same line without the
`.vtk` hash.

**The base.** Run with `a7154ec`'s own install (`@perf`'s
`worktree-bottlenecks/.venv`, a non-editable install of that commit, with
pytest put on the path from a scratch directory): the `mesh` mode twice and
the `fixtures` mode twice gave identical lines, and a third time after the
probe's save format changed. Counts on both catchments: 0 in no polygon, 0 in
more than one, 0 too narrow.

**It can fail.** With plants applied at run time to the prototype (the base
lines unchanged):

| plant | fixture lines that differ (of 66) | catchment lines that differ (of 2) |
|---|---|---|
| `_member` misses every key divisible by 7 | 49 | 2 (106 and 129 components instead of 1,860 and 2,881) |
| overlaps go to the larger code (D2 broken) | 3 | 0 |
| ties in `_best` go to the highest index | 0 | 0 |

**What the gate cannot see**, and so the red tests pin (section 7):

- the tie rule in `_best`. A component bounded by constraint lines lies in
  one polygon, so any triangle whose incircle is wider than the margin gives
  the same code; the tie only shows for a thin component, and both catchments
  have none;
- the overlap rule on real data: both catchments have no overlaps (CORINE is
  a partition), so D2 is guarded only by the fixtures.

**The prototype passes the gate**: identical `fixture` lines (66), identical
codes and counts on both catchments by `replay`.

## 7. Tests for `@tester` (the red suite)

In `tests/python/test_landcover.py`. Behaviour does not change, so the red
tests are the contracts of the two new helpers; the pins guard what the
rewrite could break and are green today, which is their point (they are
committed with the red ones, before any code). No mutation round.

**Red (fail today: the helpers do not exist).**

- **R1. `_best`, ties to the lowest index.** For `ids` built as `regions`
  builds them (each id the smallest index of its group) and `r` with ties at
  the maximum, the lowest tied index is chosen, with the groups' members
  listed in any order (interleaved groups, the tied pair not adjacent, the
  maximum not first). Include a group where every `r` is 0 (zero-area
  triangles) and single-triangle groups. The result is int64 and in
  increasing id.
- **R2. `_best` against a sort.** On a seeded random case (a few thousand
  triangles, a few dozen groups, `r` drawn from a handful of values so ties
  are common), the result equals the `lexsort` rule written in the test
  (largest `r`, then lowest index, per id). The oracle lives in the test.
- **R3. `_member`.** An empty table gives all `False` (and no error); values
  below the first entry, above the last, equal to the first and last, and
  between entries; a table of one entry; int64 values near `2**62` (the keys
  are `min * N + max`, about `N**2`).

**Pins (green today, must stay green).**

- **P1. Points on boundaries**, through `_lookup` (it exists today) or
  through `label_triangles` with hand-made components: a point on the edge
  two polygons share counts in both and gets the smaller polygon's code; a
  point on a corner shared by three polygons counts in all three; a point on
  a hole's boundary counts as in the holed polygon; a point inside the hole
  in none. Coordinates exactly representable (axis-parallel edges at UTM
  magnitudes, the existing `X0, Y0` offset).
- **P2. The caller's polygons are left as they came**: after
  `label_triangles`, `shapely.is_prepared` is `False` for each polygon passed
  in, and each one's WKB is unchanged.
- **P3. A thin component made of one zero-area triangle** (three collinear
  vertices, all three sides constraints) is labelled from its point and
  counted `thin`.
- **P4. No triangles at all**: `label_triangles` returns an empty int32 array
  and four zero counts; `regions` returns an empty array.

The existing tests (`TestRegions`, `TestDeterminism`, `TestOverlaps`) already
pin ids under permutation and the overlap rule, and stay as they are.

### What `@tester` pinned beyond this section, and the ruling

**Deviation:** these six pins (listed in the module docstring of
`tests/python/test_landcover.py@6cde0eb:31-37`) were ruled after the green
commit `3066d60`, not before it, to save time; the green code satisfies all six.

`@architect`'s ruling, 2026-10-06: **all six stand.**

1. *`_member` returns a bool array.* Stands: section 3.2 says so; the test
   only makes it checkable.
2. *`_member` takes empty `values`.* Stands: `regions` on a mesh with no
   triangles (P4) hands it an empty key array.
3. *`_best` and `_member` take numpy arrays.* Stands: they are private helpers
   fed only by `regions` and `label_triangles`, which hold numpy arrays.
4. *Both are called positionally.* Stands: it leaves the argument names free.
5. *`_lookup(points, polygons)` keeps its signature and returns `(codes,
   hits)`.* Stands: section 3.4 changes only how the pairs are found, and P1
   was allowed to go through `_lookup`.
6. *P4 also holds with vertices present and no triangles.* Stands: it is the
   shape a mesh trimmed to nothing would have, and it costs one parameter.

## 8. Net production lines, and the prototype

Prototyped in a scratch clone at `a7154ec` (removed after measuring):
`python3 tools/count_loc.py a7154ec <prototype commit>` gave **+11 net**
(19 added, 8 removed), all in `landcover.py`. The `mypy` fix of section 3.2
adds none if the return is wrapped in place. Estimate for the PR: **+10 to
+15**. Built: **+12** (20 added, 8 removed, all in `landcover.py`).

Prototype timings of `label_triangles` alone (`replay`, battery, the same
captured inputs, median of three runs):

| | Numedalslågen s | Skiensvassdraget s |
|---|---|---|
| `a7154ec` | 2.90 | 6.71 |
| prototype | 0.71 | 1.80 |

What is left is mostly the one `argsort` in `regions` (0.28 / 0.74 s), then
the union-find (0.12 / 0.38 s, six rounds on both) and the membership test
(0.08 / 0.18 s). Peak memory traced by `tracemalloc` on Skiensvassdraget:
743 MiB at `a7154ec`, 695 MiB for the prototype.

## 9. Acceptance (`@perf`)

`tools/bench.py` and the thread sweep are not required: nothing under
`include/terrain/refinement/` or `include/terrain/mesh/` changes, and labelling
runs after the mesh is finished (16c's acceptance said the same). `@perf`
runs, with the branch's own install and a base install of the merge base,
**the same power state for both, recorded with `pmset -g batt` before and
after each run**:

1. **Byte-identical meshes**: the probe's `mesh` mode on both catchments,
   every line equal to `base_a7154ec.txt` (codes, counts and the `.vtk` hash),
   and its `fixtures` mode: the base lines unchanged, the new test lines not
   compared (section 6).
2. **Time**: `rasputin mesh --stats` on both catchments, three repeats each,
   base and branch alternated, as `docs/benchmarks/2026-10-06/bottlenecks/scripts/stats.sh`
   does; the `land cover` row's median. **Pass:** the branch's median is at
   most a third of the base's on both catchments. (Checked only at these two
   inputs, 1.3 M and 3.0 M triangles; the prototype was about 4.0× and 3.6× faster.)
3. The whole-run total, recorded, not gated.

Evidence under `docs/benchmarks/<date>/30a-landcover/`.

**Result: ACCEPTED** (`@perf`, 2026-10-06, `f38ff1a`; evidence
`docs/benchmarks/2026-10-06/30a-landcover/README.md`). Base `a7154ec` against
branch `3066d60`, AC power for every run, median of 3:

| catchment | land cover, base → branch, s | total, base → branch, s |
|---|---|---|
| Numedalslågen | 2.847 → 0.674 (0.237 of base) | 13.658 → 11.407 |
| Skiensvassdraget | 6.539 → 1.744 (0.267 of base) | 20.458 → 15.679 |

The `.vtk` sha256 is the same in all six runs on each catchment and equal to
the probe's base `vtk=` hash; every base probe line is unchanged on the
branch (68 of 68), and the six new lines come from tests 30a adds.

## 10. Risks

- **Invalid polygons.** GEOS's prepared and unprepared point-in-polygon paths
  agree on valid polygons; on an invalid one (a self-crossing ring) they may
  not. Both CORINE sets are valid (section 2); a GeoJSON from a user need not
  be. Not mitigated beyond the gate: `feature_input` does not check validity
  today, and adding a check would change what is accepted, which is not this
  PR's to do.
- **A rewrite that clips polygons to the domain for speed** (as the legacy
  did) would change areas and so overlap winners. Out of scope here, and
  named so 30b does not reach for it on the land-cover side.
- **numpy's default sort** is not stable and may change algorithm between
  versions; section 4.1 is why the result does not depend on it.
- **Memory at basin scale.** `regions` holds a few arrays of `3T` int64; the
  prototype's peak is a little below today's on Skiensvassdraget (section 8).
  Not measured at São Francisco scale.

## 11. Not in scope

- A C++ flood fill over the triangle adjacency the core already has. It would
  avoid the remaining sort. **An untested idea**, not measured; a possible
  later step if the São Francisco runs show `regions` costing more than
  refine (16c's R2 trigger).
- The CORINE clip (30b) and DEM reading (30c), the next two PRs.
- Any change to the stderr line's wording or the codes' meaning.

## 12. Questions for Ola

1. **After this PR, should moving the component search into C++ wait for the
   São Francisco runs?** After this change, land cover takes about 0.7 s and 1.8 s on the two
   catchments, under refine's 1.2 s and 3.2 s, so 16c's rule for moving
   it into C++ no longer fires. **Default: yes, wait; the C++ flood fill
   stays an untested idea until a basin run shows land cover costing more
   than refine.**
2. **Is the gate on two catchments plus the test fixtures enough, given that
   neither catchment has overlapping polygons or slivers too thin to label?**
   Those two rules are then guarded by the fixtures and the new tests only.
   **Default: yes, enough; no third catchment is added.**

## Review

**Design review, round 1, 2026-10-06.** Range `a7154ec..78b682e`. Verdict: APPROVED. LOC: 0 (design only); estimate +10 to +15. Not pushed; no CI.

**Code review, round 1, 2026-10-06.** Range `a7154ec..46a5e56` (design 78b682e, red 6cde0eb, green 3066d60, @perf f38ff1a ACCEPTED, records 46a5e56). Verdict: CHANGES REQUESTED. LOC: +12 net (20 added, 8 removed, all `src_python/tin_engine/landcover.py`) against +10 to +15. pytest 5221 passed, 17 skipped; mypy, ruff, prohibited-deps, detria boundary, check_citations clean; probe gate re-checked: 68 of 68 base lines unchanged, 6 new, .vtk hashes equal. Blocking: red-step scaffolding at `tests/python/test_landcover.py@46a5e56:464-477` (`helper`'s call-time getattr lookup) and the docstring's 30a paragraph at :31-33; use `landcover._best`/`_member` directly (test-audit R10) and put the paragraph in the past tense. The design record's range and LOC 0 confirmed by count_loc. The profile commit a7154ec goes out in this PR; round 2 reviews `8199f30..` whole. Not pushed; no CI.

**Code review, round 2, 2026-10-06.** Range `8199f30..3cc526f`. Verdict: CHANGES REQUESTED. LOC: +12 (estimate +10 to +15). Round 1's red-step item fixed in `3cc526f`; blocking: six land-cover figures in docs/benchmarks/2026-10-06/bottlenecks/README.md section 3 (two copied into this file's section 2) disagree with raw/lc_*.json, its battery percentages are absent from raw/stats/*_power.txt, and 30a-landcover/README.md's "spread under 3 %" is 3.4 % on Skiensvassdraget's branch land cover; prose only. Merges cleanly with origin/master 9dc3c11. Not pushed; no CI.

**Code review, round 3, 2026-10-06.** Range `3cc526f..5b910de5` (prose only). Verdict: APPROVED. LOC: +12 for `8199f30..5b910de5` (20 added, 8 removed, all in `src_python/tin_engine/landcover.py`; estimate +10 to +15). Round 2's three findings are fixed: the profile table's figures match `raw/lc_*.json`, the battery percentages are gone, and the spread sentence is correct (3.4 % on Skiensvassdraget's branch land cover, 2.2 % or less on the other seven rows). The 58 at-risk citations from `check_citations.py --base 8199f30` were re-read: every cited line is the same at `8199f30` and at the head. The branch merges cleanly with origin/master `9dc3c11`. Not pushed; no CI.
