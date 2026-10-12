# Increment 35 — the land-cover clean-up step, made faster

Status: **merged (#229, 2026-10-11)** from branch `worktree-cleanup-speed` (design `7da1d663`, review
round 1 approved `7218dfb7`; red `c4130194` (`@tester`), green `e805d37f`
(`@developer`, 23 net production lines by `python3 tools/count_loc.py c7dba32e
HEAD`: `crs.py` +14, `feature_input.py` +9), mutation round `7621a99e` (M1-M4
all killed, one count-based test added for M3); section 11 is the as-built
record; `@perf`'s quick check `2cf224a6` ACCEPTED, section 11; code review round 1 APPROVED, Review section; Numedalslagen's 68 % hotspot ruled expected by Ola 2026-10-12, `docs/benchmarks/quick/hotspots.toml`). Master
`c7dba32e` is `@perf`'s profile commit
(`docs/benchmarks/2026-10-11/cleanup-profile/README.md`). Every `file:line`
below is pinned to `c7dba32e`. Every figure was measured on this worktree's
`.venv` with the scripts in `docs/increments/35-probes/` (its README has the
raw output). Ola picked this row on 2026-10-11 ("land-cover clean-up speed-up,
then GeoPackage"). Two questions for Ola, decided on their defaults while he
was away (section 9).

**What this is.** `rasputin mesh --features` with a class map runs 20c-3's
land-cover clean-up: clip to the read region, repair (`coverage_clean`), merge
by class, Ola's outline rule, clip to the domain, the border simplifier (32),
cut where on the outline, then each line through `_add`. On Numedalslagen the
step is 6.08 s of a 9.06 s run (67 %), on Lagan 7.17 s of 19.61 s (37 %), and
before the step Lagan reprojects every CORINE polygon whole, 15.9 million
vertices, 4.08 s (20 % of the run).

**What changes for a user.** Nothing in any output: every `.vtk` byte, the
`features:` stderr line and the run record stay as they are. Lagan's
`features clip` phase goes from 11.7 s to 6.7 s and the run from 20.6 s to
16.3 s (AC power, one run each, section 7); Numedalslagen's phase from 6.3 s
to 6.0 s. Numedalslagen's step stays near 6 s: what is left is GEOS overlays
whose output no cheaper call reproduces (section 4).

**Lean.** Python only, four small changes, about 20 net lines (section 6).
No C++; the simplifier is not touched. Mutation round only on the
invariant-critical tests of section 5.

## 1. Prior art: legacy and literature

### Literature

No novelty is claimed.

- **Threads over independent points (3.1).** A map projection is applied
  point by point; PROJ's `proj_trans_generic` transforms an array in one
  call, and pyproj releases the GIL around it, so chunks of one array on
  several threads give the same bits per point as one call (measured, 15.9
  million points: `same bits: True` at 2, 4 and 8 threads). Each thread uses
  its own `Transformer`, as pyproj's documentation asks ("Transformer objects
  are not thread-safe"; pyproj 3.8, *Advanced examples: multithreading*,
  recalled, not reread).
- **The touching-line test (3.2).** The DE-9IM (Clementini, Di Felice and
  van Oosterom, *A small set of formal topological relationships suitable for
  end-user interaction*, SSD 1993; Egenhofer and Franzosa 1991). The pattern
  `T**F**FF*`, polygon first, is OGC's `covers` (`T*****FF*`) with one more
  `F`: the polygon's boundary meets no interior point of the line. Where it
  holds, the line's point set lies in the polygon and the overlay adds no
  node inside the line, so the result is the input edge; JTS's OverlayNG
  `LineBuilder` emits an unnoded result edge with the input's vertices and
  direction (recalled from JTS's OverlayNG design notes, not reread; the
  argument stops there, the identity itself is measured: 921 of 921 such
  lines on the two catchments, `clean_probe.py`, row C4). 30b's rule
  (`docs/increments/30b-clip-speed.md`, 3.4) is the same filter-and-refine
  pattern with the stricter `contains_properly`; 35 widens the filter, the
  refine step (GEOS's intersection) is unchanged for every other line.
- **One index, many queries (3.4).** Nothing to cite: a tree built once
  instead of once per polygon.

### Legacy

`git grep -n -iE "ThreadPoolExecutor|relate_pattern|Transformer\.from_crs|\.covers\(" legacy-archive -- legacy`
returns `legacy/rasputin/application.py:118`, `legacy/rasputin/avalanche.py:21`,
`legacy/rasputin/geometry.py:237`, `legacy/rasputin/globcov_repository.py:142`,
`legacy/rasputin/gml_repository.py:194`, `legacy/tests/test_gml_repository.py:25`
and `legacy/tests/test_mesh.py:24`: one-shot pyproj transforms, no threads, no
DE-9IM, nothing about cost. Nothing is carried across.

## 2. What is slow (measured, `clean_probe.py`, battery, one run)

Seconds inside the step, Numedalslagen / Lagan; the candidate's column says
whether a cheaper call gives the same bytes.

| sub-step (`feature_input.py`) | N | L | candidate | same bytes |
|---|---|---|---|---|
| reprojection of the whole cover polygons, before the step (`:357`) | 0 | 4.01 | C1 one `shapely.transform` over all polygons: 4.08 s, **no saving alone**; C1' chunks on threads: 3.94 → 0.63 s at 8 | yes, bit for bit |
| clip to the read region (`:420`) | 0.68 | 0.55 | C2 skip where the region contains the polygon properly | **no** (section 4) |
| repair, `coverage_clean` (`:427`) | 1.15 | 0.86 | none: one GEOS call | — |
| merge by class (`:435`) | 0.12 | 0.08 | none | — |
| outline rule (`:640`): `union_all` of rebuilt parts `:680` / `make_valid` of shells `:676` / `_Outline.ring` `:776` / the rest | 0.56 / 0.29 / 0.56 / 0.84 | 0.36 / 0.51 / 0.42 / 0.84 | none exact (section 4) | — |
| clip to the domain (`:444`) | 0.81 | 0.70 | C6 part by part | **no** (section 4) |
| simplifier (C++) | 0.68 | 0.64 | out of scope | — |
| cut where on the outline, `_outline_lines` (`:714`) | 0.27 | 0.71 | C5 one `STRtree` per step: saves 0.06 / 0.22 | yes, by construction |
| `_add`: `intersection(line, dom)` for lines not `_untouched` (`:395`) | 0.15 | 0.62 | C4 the DE-9IM pattern: 71 of 626 / 850 of 875 lines need no overlay | yes, 921 of 921 |
| `_add`: `covers` per line (`:402`) | 0.02 | 0.76 | C3 skip where `_untouched` held | yes, by logic |

## 3. The blueprint

All four changes are Python inside the I/O boundary (`CLAUDE.md` §2); the
C++ core, the simplifier and every model are untouched.

```
crs.reprojector(src, dst) -> apply(xy)                       (3.1 b)
   len(xy) >= 2 * CHUNK_POINTS: chunks on a ThreadPoolExecutor,
   one Transformer per thread (threading.local); else one call

feature_input._Tally._take                                   (3.1 a)
   cover polygons appended unmoved; after the loop, one
   shapely.transform(array of cover polygons, move); _finite per polygon

feature_input._untouched(line, domain)                       (3.2)
   relate_pattern(domain, line, "T**F**FF*") replaces contains_properly

feature_input._Tally._add                                    (3.3)
   covers only for the lines _untouched did not pass

feature_input._clean -> _outline_lines(polygon, outline, tree)  (3.4)
   the outline's segment STRtree built once per _clean, passed in
```

### 3.1 The reprojection in one call, on threads (`feature_input.py:357`, `crs.py:42-52`)

(a) `_take`'s cover branch (`:355-363`) appends `(fid, mask, code,
geometry)` unmoved. After the loop, when `move is not None`, one
`shapely.transform(np.array(geoms, dtype=object), move)` moves every cover
polygon (shapely gathers all their coordinates and calls `move` once); then
`_finite` per moved polygon, raising today's `FeatureError` naming the fid
(`:359-361`), and the `cover` list is rebuilt with the moved geometries. The
time goes to `self.seconds` as today. Measured: one call over all polygons
gives the same bits as per-polygon calls (`C1 ... same bits: True`).

(b) `reprojector`'s `apply` (`crs.py:48-51`): when `len(points) >= 2 *
CHUNK_POINTS`, split with `np.array_split` into `min(WORKERS, len(points) //
CHUNK_POINTS)` chunks and map them over a `ThreadPoolExecutor`, each worker
thread holding its own `Transformer` in a `threading.local`; concatenate.
Smaller arrays take today's single call. Constants: `CHUNK_POINTS =
250_000`, `WORKERS = min(8, os.cpu_count() or 1)`. Scale: Lagan's 15.9
million points, 8 threads, 3.94 → 0.63 s; a 250 000-point chunk is about
60 ms of PROJ work, far above a thread's start cost. Largest input checked:
those 15.9 million points. Below 500 000 points no thread starts, so every
other caller (domain, region, 33's lines) is byte-for-byte today's code
path. Why both halves: without (a), no single polygon reaches two chunks and
Lagan's reprojection stays at 4.05 s (probes README, item 4).

### 3.2 `_untouched` widened to lines whose ends touch the outline (`:244-254`)

Clause 3 of 30b's rule becomes `shapely.relate_pattern(domain, line,
"T**F**FF*")`: the line's interior lies in the domain's interior and its
ends in the closure. Clauses 1, 2 and 4 stay; the order stays. `domain`
first: with GEOS 3.13+ the prepared polygon's RelateNG answers it (measured
0.01 s for 875 lines against 0.76 s with the line first). A line whose
interior touches the outline (along it, or at a domain vertex) still fails,
as today. Where the new clause holds and the old did not, GEOS's
`intersection(line, dom)` returned the line itself in every measured case
(921 of 921, both catchments); the mechanism is in section 1.

### 3.3 `covers` only where it can be false (`:395-402`)

`_add` evaluates `_untouched` once per chain, keeps the flags, and counts
`clipped` as `not all(flag or dom.covers(g) ...)`. A chain in the domain's
closure is covered, so the count is unchanged by logic, not by measurement.
(Both `features:` lines were nevertheless identical in `mesh_check.py`.)

### 3.4 One outline tree per clean-up (`:714-735`)

`_outline_lines(polygon, outline, tree)` takes the `(segments, STRtree)`
pair; `_clean` builds it once before the comprehension at `:448`. Same tree,
same queries, same output.

### 3.5 Against the assessment framework

1. Data and execution: no option added; two module constants.
2. State: `reprojector`'s thread-local `Transformer`s live inside the
   closure; `_untouched`, `_outline_lines` stay pure.
3. Dependencies: the standard library's `concurrent.futures` and `threading`;
   nothing from §2's list.
4. Async: the step runs off the loop today; the pool is created and closed
   inside `apply`, so no thread outlives the call.

## 4. Left out, and why

- **C2, no clip for polygons inside the read region, and C6, the domain
  clip part by part** (0.19 + 0.36 s on Numedalslagen, 0.13 + 0.29 s on
  Lagan). GEOS's overlay returns the same point set but rewrites every ring:
  shell clockwise, starting at the input's last vertex, holes reordered
  (`ring_order.py`: 3066 of 3077 follow one rule, 11 do not), and neither
  `coverage_clean` nor the simplifier washes that order out
  (`stage_check.py`: `same WKB: False` at both stages). Byte-identity would
  mean reproducing OverlayNG's ring construction; not for half a second.
- **The outline rule's GEOS calls** (`union_all` of rebuilt parts,
  `make_valid` of shells; 0.85 / 0.87 s) are overlays with the same
  property. Its Python loop (`_Outline.ring`, 0.25 / 0.13 s of own time)
  would cost tens of lines to vectorise for a tenth of a second.
- **`coverage_clean`** (1.15 / 0.86 s): one GEOS call. The way to drop it is
  Ola's idea "Repair CORINE once per area of interest" (private ideas file):
  a cached repaired extract. A separate row; it changes what a run reads.
- **Clipping in the source CRS before reprojecting**: would cut Lagan's 15.9
  million vertices to about 0.4 million before PROJ sees them, but the
  overlay's ring rewrite (above) changes the vertex order the region clip
  then starts from. Listed, not pursued; 3.1 takes the same time away.
- **Numedalslagen** keeps a step of about 6 s. Every remaining second is a
  GEOS overlay (repair, outline rule, domain clip) or the C++ simplifier.

## 5. Tests for `@tester` (the red suite)

Small fixtures; `tests/python/test_feature_input.py` (30b's `TestUntouched`
table) and `tests/python/test_crs.py`. No real catchment in CI.

1. **Chunked transform, bit for bit** (`reprojector`): 1.2 million points
   drawn uniformly in EPSG:3035 over x 4.3-4.9e6, y 3.5-4.1e6 (southern
   Sweden), to EPSG:3006, against one `Transformer.transform` call;
   `np.array_equal` on both axes. Also 10 points (the single-call route)
   and exactly `2 * CHUNK_POINTS` points. Measured here: 15.9 million
   points equal at 2, 4 and 8 threads.
2. **One transform per source** (`_take`): a GeoJSON source in EPSG:3035
   with three polygons under a `corine` map and a clean-up request;
   `fi.reprojector` monkeypatched to wrap the real one and count calls of
   `apply` over polygon coordinates: one call carrying all the cover
   polygons' vertices. And a polygon with a vertex PROJ cannot map (x =
   1e300) is refused with today's text, naming its fid.
3. **`_untouched` table**: `UNTOUCHED_CASES` gains `ends-on-the-outline`
   (both ends on the `EXTRA` domain's outline, every interior point inside
   it) → True, and `one-end-on-the-outline` → True; `along-the-outline`,
   `along-then-inside`, `domain-vertex-on-the-line`, `touching-a-hole`
   stay False. `test_r2` (GEOS gives the line back where True) runs over
   the new True cases too: this oracle is what makes 3.2 safe.
4. **`clipped` count** (`_add`): a feature of two chains, one untouched and
   one crossing the outline → `clipped` 1; two untouched → 0; a spy that
   `covers` is called once in the first case and never in the second.
5. **One tree per clean-up**: `shapely.STRtree` spied during one
   `open_features` with the CORINE fixture
   (`tests/fixtures/corine/clc2018_7908_3.gpkg`, the `dtm10` fixture's
   domain): built once inside `_outline_lines`' call chain per source, not
   once per polygon; the lines equal to those of `_outline_lines` called
   with a fresh tree.
6. **The set is unchanged**: `open_features` on the CORINE fixture twice in
   one test, once as built, once with `fi._untouched` monkeypatched to 30b's
   rule (kept verbatim in the test) and `fi.reprojector` to a single-call
   transform; every feature's lines' WKB, `clipped`, `outside`, `empty` and
   `cover_vertices` equal.

**Invariant-critical** (mutation round): tests 1, 3 and 4. Mutation
targets: M1 the pattern `T*****FF*` (plain `covers`; killed by
`along-then-inside`, `domain-vertex-on-the-line` and `touching-a-hole`, mutant
True where False is expected; `along-the-outline` does not kill it, both
patterns return False there), M2 `is_simple` dropped (`self-touching-at-a-vertex`), M3
`clipped` counted from the flags alone (test 4, first case), M4 a chunk
boundary off by one in `apply` (test 1: a point dropped or doubled).

## 6. Size

| file | counted lines, net (estimate) |
|---|---|
| `crs.py` (`reprojector`) | about +12 |
| `feature_input.py` (`_take`, `_untouched`, `_add`, `_clean`, `_outline_lines`) | about +8 |
| total | **about +20**, far under 700 |

Check: `python3 tools/count_loc.py c7dba32e <head>`. Short path
(`docs/increments/README.md`): under 50 lines, no C++, nothing the
refine/mesh acceptance covers; `@reviewer` may send it to the full path.

## 7. Acceptance (`@perf`)

`tools/bench_quick.py run` on the branch head against its baseline for the
power state. Every mesh byte-identical to master's on every case: on
`numedalslagen` and `lagan` because the step ran, on `tile`, `quarter`,
`geilo-al-ramp` and `romsdal-slope` by construction (no `--features`, so
`_Tally` is never built; `cases.toml` shows it). Targets, as ratios to the
baseline's phase medians: Lagan `features clip` at most 0.65 (measured
here 11.74 → 6.74 s, 0.57; the clean-up row 7.55 → 5.76 s) and Lagan's
run at most 0.85 (20.64 → 16.28 s); Numedalslagen `features clip` at most
1.0 (6.29 → 6.00 s; the 5 % saving sits inside the quick check's band, so
there the check is identity, not speed). Numbers here are one AC run each;
`@perf`'s medians rule.

## 8. Risks

- PROJ's choice of operation is made once per `Transformer`; a thread's
  own `Transformer` is built from the same CRS pair on the same machine, so
  it picks the same one. A PROJ network or grid download on one thread but
  not another is not possible here (network off, pyproj default).
- `relate_pattern` without GEOS 3.13's prepared RelateNG is correct but
  slow (0.76 s on Lagan with the line first). shapely 2.2 wheels ship GEOS
  3.13 or later; the floor `shapely>=2.2` (`pyproject.toml:38`) holds.
- 3.2 rests on OverlayNG keeping an unnoded line edge as the input gave it,
  measured on 921 lines, not read in GEOS's source; test 3's oracle guards
  it on every GEOS CI runs.

## 9. Questions for Ola

1. **Threads for the reprojection**: up to 8, from the machine's core count,
   independent of `--threads` (which sizes the C++ refinement). *Default:
   yes.* A "no, follow `--threads`" adds the value to `FeatureRequest` and
   about 5 lines.
2. **Numedalslagen's step stays near 6 s.** Every exact change is in; what
   remains is GEOS overlay time or Ola's cached repair idea, which changes
   the inputs a run reads and so is its own row. *Default: ship 35 as
   designed; the cached repair stays an idea.*

## 10. ROADMAP

Row 35 in `ROADMAP.md`.

## 11. As built (recorded by the main session from the handbacks)

**Red** `c4130194` (`@tester`): 15 failing tests in `tests/python/test_crs.py`
and `tests/python/test_feature_input.py`, 275 passing; `TestChunkedReprojector`
(test 1: `crs.CHUNK_POINTS`, `crs.WORKERS` missing; also `2 * CHUNK_POINTS + 1`
points and `WORKERS == min(8, os.cpu_count() or 1)`), `TestUntouched` (test 3:
four new True rows red under `contains_properly`, each in `test_r1` and the
per-row GEOS oracle `test_r2`; table adds `end-on-a-domain-vertex` True and
`ends-on-the-outline-through-a-domain-vertex` False, measured on GEOS 3.14.1),
`TestOneTransformPerSource` (test 2: apply calls `[5, 5, 5]`, expected `[15]`),
`TestClippedCount` (test 4: `covers` called 2 times, expected 1 and 0; calls
`_Tally(BOX, CRS, V)._add(...)` directly, spy on both `BaseGeometry.covers` and
`shapely.covers`), `TestOneOutlineTreePerCleanup` (test 5: 8 trees under
`_clean`, expected 1; quarter-circle domain), `TestTheSetIsUnchanged` (test 6,
the guard for green). No red-step scaffolding. No real-catchment test in CI;
byte-identity is `@perf`'s (section 7).

**Green** `e805d37f` (`@developer`): `src_python/tin_engine/crs.py` (+14) and
`src_python/tin_engine/feature_input.py` (+9) only; 23 net counted lines
against the estimate of about 20. Whole Python suite 7143 passed, 17 skipped;
mypy strict, ruff, the governance gates clean. Output identity at `e805d37f`,
the quick check's two catchment cases (`--out X.vtk --binary`, whole-file
sha256): numedalslagen
`349d0b5632a28af6024210556a0bd0b6f35e5ba99378c18c742d4b666acfdcfe`, lagan
`6c110f9eb1fb6a741bd4a19a954a91b7a7424d43f7ba04c354d25f9a42495fa0`, both equal
to the design's digests (`35-probes/README.md`, item 4); both `features:`
stderr lines identical. Pinned beyond the design: `_outline_lines(polygon,
tree)` drops the unused `outline` argument (the `(segments, STRtree)` pair
carries its segments; `_clean` is the sole caller); the calling thread's
`Transformer` is still built eagerly in `reprojector`, so a bad CRS refuses at
construction as before, worker threads build theirs lazily in a
`threading.local`; `_finite` is checked per cover polygon after the batch, as
the old per-polygon branch did.

**Mutation round** `7621a99e` (`@tester`), mutants planted in a
`tools/scratch_copy.py` copy of `e805d37f`, unmutated run green after each.
The record, verbatim from the handback (the commit message carries it too):

| mutant | planted in | killed by | result |
|---|---|---|---|
| M1 `T*****FF*` | `_untouched` | `test_r1` rows along-then-inside, domain-vertex-on-the-line, touching-a-hole, ends-on-the-outline-through-a-domain-vertex, both closed-ring-touching rows; `test_r2` seeded. along-the-outline does not kill it (as the review said) | killed, 7/36 |
| M2 `is_simple` dropped | `_untouched` | self-touching-at-a-vertex, self-crossing, doubling-back; `test_r2` seeded | killed, 4/36 |
| M3 flags-alone count | `_Tally._add` | test 4 first case, but only by the `covers` spy: the count stays 1 there | killed by spy; new test added (below), re-planted: killed by the count |
| M4 dropped / doubled boundary point | `crs.reprojector.apply` | `test_one_point_above_the_threshold`; plus exactly-two-chunks and 1.2-million rows for the doubled variant | killed, 1/36 and 3/36 |

M4 is local evidence (10 cores; one chunk on a one-core runner). New test for
M3: `TestClippedCount::test_a_covered_chain_the_rule_refused_is_not_clipped`
(a chain along BOX's east side kept whole, `clipped` 0, `covers` asked once);
the design's test 4 cases could not tell M3 apart by the count. Lesson: a
spy-only kill hides a count mutant; check which assertion fails.

**Decided on defaults while Ola was away** (section 9): Q1 reprojection
threads `min(8, cpu_count)` independent of `--threads`, yes; Q2 ship with
Numedalslagen's step near 6 s, the cached-repair idea staying an idea, yes.


**Acceptance** `2cf224a6` (`@perf`, ACCEPTED against section 7;
`docs/benchmarks/2026-10-11/cleanup-speed/README.md`): head `357808b1` against
the branch's base `3e6f87eb` measured back to back with `--tree`, AC power
before and after every run, medians of 3. Mesh SHA-256 equal on all 7 quick
cases. Lagan features clip 11.156 -> 6.342 s (ratio 0.569, limit 0.65), Lagan
total 19.296 -> 14.542 s (0.754); Numedalslagen features clip 6.092 -> 6.000 s
(0.985; clean-up 5.94 s, near 6 s as Q2 expected); other totals -0.7 % to
+0.8 %. Tool verdict `FASTER`. `docs/benchmarks/quick/baseline-ac.json`
refreshed to the head record (the stored one predated increment 32, and the
venv's numpy 2.5.4 differed from its 2.5.3). Hotspots: Lagan features clip
43.6 % of the run (ruled 40.7 %; 57.8 % on master after increment 32);
Numedalslagen features clip 68 %, unruled, GEOS overlay time (section 4).
Decided on defaults while Ola was away: the 68 % is recorded as expected here;
the `hotspots.toml` entry and the Lagan 40.7 % ruling wait for Ola's word.
Ruled by Ola 2026-10-12 ("Rule it expected at 6 s, cached extract is the fix"):
the entry is in `docs/benchmarks/quick/hotspots.toml` at 68.0 %; the fix is the
cached repaired extract per area of interest, when picked. Lagan's 43.6 % is
within 10 points of its 40.7 % ruling, so the tool is quiet; no new entry.
Lesson: a stored quick-check baseline goes stale when a merge changes outputs
without a refresh (32 merged without one).

## Review

Design review round 1, 2026-10-11, @reviewer (7da1d663, prose only, 0 counted lines): APPROVED. Every figure in sections 2, 4 and 7 matches /Users/skavhaug/projects/rasputin/.claude/worktrees/cleanup-speed/docs/increments/35-probes/README.md@7da1d663 and the profile README; every feature_input.py and crs.py citation quotes the pinned code; both mesh digests recorded; check_citations clean. The relate pattern `T**F**FF*` read in words: the line's interior meets the domain's interior and neither its boundary nor its exterior, and no end lies outside; a line touching the outline only at an interior vertex is excluded, which is necessary since GEOS nodes it there. One correction applied in this recording commit: /Users/skavhaug/projects/rasputin/.claude/worktrees/cleanup-speed/docs/increments/35-cleanup-speed.md@7da1d663:234 claimed `along-the-outline` kills M1; it does not (measured, GEOS 3.14.1), M1 is killed by `along-then-inside`, `domain-vertex-on-the-line` and `touching-a-hole`. Nit applied: the comprehension cited at `:450` is at `:448`. Ideas for @tester: add `end-on-a-hole` to test 3 (new rule True, GEOS returns the line unchanged, measured); M4 cannot be killed on a one-core runner, so the mutation round stays local. Lesson: a claimed mutant kill is checked by running the fixture table against the mutant.

Code review round 1, 2026-10-11, @reviewer (3e6f87eb..4de504af, 23 net counted lines by count_loc.py from both 3e6f87eb and c7dba32e, estimate about 20; not pushed, no CI yet): APPROVED. The green diff touches only the two production files and does exactly the four section-3 changes plus the three section-11 pins; the pattern is `relate_pattern(domain, line, "T**F**FF*")` with the domain first (/Users/skavhaug/projects/rasputin/.claude/worktrees/cleanup-speed/src_python/tin_engine/feature_input.py@4de504af:255); chunking cannot change bits (same Transformer per thread, `np.array_split` keeps order, PROJ maps point by point, arrays under 500k take the old single call: /Users/skavhaug/projects/rasputin/.claude/worktrees/cleanup-speed/src_python/tin_engine/crs.py@4de504af:57-72); one threading.local per closure, the pool opened and closed inside `apply`; every mutant M1-M4 names an existing test and fixture row (/Users/skavhaug/projects/rasputin/.claude/worktrees/cleanup-speed/tests/python/test_feature_input.py@4de504af:1274-1295, 1680-1717; /Users/skavhaug/projects/rasputin/.claude/worktrees/cleanup-speed/tests/python/test_crs.py@4de504af:219-237); no red-step scaffolding; the acceptance README's figures equal the raw medians in its base and head records; check_citations clean; the two test files 291 passed. Correction applied in this recording commit: the ROADMAP row's stale "next @tester" phrase. Non-blocking note: `_finite` on cover polygons runs after the loop (/Users/skavhaug/projects/rasputin/.claude/worktrees/cleanup-speed/src_python/tin_engine/feature_input.py@4de504af:381-386), so when a cover polygon and a later line feature are both unmappable, the fid named is the line's.

Ruling-record review round 1, 2026-10-12, @reviewer (`origin/master..bc7639fb`, 0 counted lines): CHANGES REQUESTED, two prose claims, applied in this recording commit: /Users/skavhaug/projects/rasputin/.claude/worktrees/hotspot-ruling/ROADMAP.md@bc7639fb:69 said "code review rounds 1 and 2 APPROVED" for 35, which had one design round and one code round; /Users/skavhaug/projects/rasputin/.claude/worktrees/hotspot-ruling/docs/benchmarks/quick/hotspots.toml@bc7639fb:16 said "about 85 % GEOS calls", the profile README gives the three largest sub-steps 80 % of the step with about 60 % in GEOS proper. Verified: PR merges and dates, every hash, 6.000 s of 8.829 s = 68.0 %, Lagan 43.6 under 40.7 + 10 (`tools/bench_quick.py:166`), TOML loads and the shipped-file tests pass, check_citations exit 0. Suggestion applied: increment 36's Status line said "built ... Next: the push", now merged (#231). Lesson: a percentage in a comment names the README line it sums from.
