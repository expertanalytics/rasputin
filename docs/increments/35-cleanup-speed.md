# Increment 35 — the land-cover clean-up step, made faster

Status: **designed** by `@architect`, 2026-10-11, on branch
`worktree-cleanup-speed` at master `c7dba32e` (`@perf`'s profile commit,
`docs/benchmarks/2026-10-11/cleanup-profile/README.md`). Every `file:line`
below is pinned to `c7dba32e`. Every figure was measured on this worktree's
`.venv` with the scripts in `docs/increments/35-probes/` (its README has the
raw output). Ola picked this row on 2026-10-11 ("land-cover clean-up speed-up,
then GeoPackage"). Two questions for Ola, with defaults (section 9).

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
pair; `_clean` builds it once before the comprehension at `:450`. Same tree,
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
targets: M1 the pattern `T*****FF*` (plain `covers`: `along-the-outline`
must fail it), M2 `is_simple` dropped (`self-touching-at-a-vertex`), M3
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
