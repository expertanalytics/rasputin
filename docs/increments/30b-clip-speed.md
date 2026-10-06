# Increment 30b — the `features clip` phase of `rasputin mesh`, made fast

Status: **designed, not built.** Designed by `@architect` 2026-10-06 on branch
`worktree-clip-speed`, on `@perf`'s profile commit `a7154ec` (master
`8199f30` plus the profile). Master has since gained #195 and #196. Of the
files this PR touches or reads, only `cli.py` changed, and not in
`_open_features` or anything it calls (`git diff a7154ec origin/master --
src_python/tin_engine/cli.py` touches imports, `_off_node`,
`_placement_report` and `catchment`; `git diff --stat a7154ec origin/master --
src_python/tin_engine/feature_input.py src_python/tin_engine/io/geopackage.py
tests/python/test_feature_input.py tests/python/feature_fixtures.py ROADMAP.md`
prints nothing). So the branch was not merged with master for this design.
Amended by `@architect` the same day after design review round 1 (`## Review`):
the whole-feature test now tests the feature's linework (3.2, 4.2), with pin
P5 and three new probe cases.

**What this is.** The second of three pull requests that remove the
bottlenecks `@perf` measured in `rasputin mesh` on two real catchments
(`docs/benchmarks/2026-10-06/bottlenecks/README.md`). Ola approved them on
2026-10-06 in this order: land cover (30a,
`docs/increments/30a-landcover-speed.md` on `worktree-landcover-speed`), then
the CORINE clip (this one, 30b), then reading the DEM (30c). The ROADMAP row
is 30 (section 13 says how this branch records it).

**What changes for a user.** Nothing in any output. Every kept feature's lines
(their vertices and order), its polygon, and the `features:` stderr line's
four numbers (kept, cut at the domain outline, outside, empty) stay exactly as
they are, and so do the noded graph and the mesh. The phase gets about three
to four times faster: on AC power, 5.9 s to 1.9 s on Numedalslågen and 5.7 s
to 1.3 s on Skiensvassdraget, measured on a prototype (section 8).

**Lean.** No mutation round (Ola's standing rule for lean rounds). No C++.

## 1. Prior art: legacy and literature

### Literature

The method is increment 16b's and does not change
(`docs/increments/16b-terrain-polygons.md`, R5 and R6): read the features
whose R-tree box meets the box of a region around the domain, keep each ring's
or line's whole edges that meet the region (the pre-clip), then cut what is
kept with the domain as linework. What changes is that cheap tests answer
first and the exact test runs only where they cannot decide. This is the
**filter-and-refine** pattern of spatial query processing (Orenstein, SIGMOD
1986; Brinkhoff, Kriegel, Schneider and Seeger, *Multi-step processing of
spatial joins*, SIGMOD 1994): a cheap, conservative filter, then the exact
test on what is left. Recalled, not reread. Here the filters are exact rather
than conservative, so nothing is left for a refine step to correct:

- **Whole feature first.** One test of the feature's linework (every ring of a
  polygon, holes included, or the line itself) against the region. If the
  linework misses the region, none of its edges can meet it.
- **Endpoints first.** For each edge, test its two vertices against the region
  (point in polygon). A closed edge whose endpoint lies in the closed region
  meets the region. Only edges with both ends outside get the segment test.
- **Prepared geometry** (JTS `PreparedGeometry`, Martin Davis; GEOS's port)
  answers each test against an index of the region's or the domain's edges
  instead of walking all of them. The facts this design relies on were read in
  the source of the versions installed (shapely 2.1.2, GEOS 3.13.1):
  - shapely's binary predicates use the prepared form of their **first**
    argument only. In `src/ufuncs.c`, `YY_b_p_func` calls the prepared GEOS
    function when `ip1` is prepared and never looks at `ip2`. So today's
    `shapely.intersects(segments, region)` in `pre_clip` runs unprepared,
    although `pre_clip` prepares `region` two lines above. The fix is to swap
    the arguments.
  - `shapely.intersects_xy(region, x, y)` calls `GEOSPreparedIntersectsXY_r`
    (`Ydd_b_p_func` in `src/ufuncs.c`). It prepares an unprepared geometry
    afresh for every point (so `region` must already be prepared), and it
    returns `False` for a NaN coordinate. GEOS's `GEOSPreparedIntersectsXY_r`
    (`capi/geos_ts_c.cpp`, 3.13.1) sets a reused point and calls
    `GEOSPreparedIntersects_r`.
  - GEOS's `PreparedPolygon::intersects` (`src/geom/prep/PreparedPolygon.cpp`,
    3.13.1) first returns `false` when the two **envelopes** (bounding boxes)
    do not meet. A polygon's envelope is its shell's: shapely gives `bounds`
    `(20, 20, 30, 30)` for a shell (20..30)² with a hole (1..2)². So for an
    invalid polygon whose hole lies outside its shell, the hole is outside the
    envelope GEOS tests (design review round 1, B1; section 4.2).
  - Past the envelope test (and a separate path when the region is an
    axis-aligned rectangle, which a buffered hull is not),
    `PreparedPolygonIntersects::intersects`
    (`src/geom/prep/PreparedPolygonIntersects.cpp`, 3.13.1) first locates one
    point of each component of the test geometry with the region's
    point-in-area locator (boundary counts). The components are visited by a
    `GeometryComponentFilter` (`LocationNotMatchingFilter` in
    `PreparedPolygonPredicate.cpp`), which reaches each ring of a polygon and
    each line of a multi-line, and takes its first coordinate. It then looks
    for any segment of the test geometry meeting the region's boundary, and
    for a test geometry of dimension 2 it checks whether the region lies
    inside it. Section 4.2 rests on this order.
- **The R-tree** is SQLite's R*Tree module (Beckmann, Kriegel, Schneider and
  Seeger, SIGMOD 1990), which the GeoPackage R-tree extension uses (OGC
  12-128r18 Annex F.3, as `io/geopackage.py` cites it). From the SQLite
  documentation (`sqlite.org/rtree.html`, sections 3.4 and 6, read for this
  design):
  - A built-in query is a range constraint on the box columns. A query by any
    other shape needs a callback registered through the C API
    (`sqlite3_rtree_query_callback`), which Python's `sqlite3` module does not
    expose.
  - Boxes are stored as 32-bit floats, rounded outward, so an overlap query
    never misses an entry.

  So the read cannot ask the index for "features meeting the hull". It could
  only filter the index's boxes in Python after the query, and section 11 says
  why this PR does not do that either.

**Novelty: none claimed.** This is standard engineering on a standard method,
so nothing was searched.

### Legacy

```sh
$ git grep -l -i "intersection\|clip\|prepared" legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
legacy-archive:legacy/rasputin/geometry.py
legacy-archive:legacy/rasputin/gml_repository.py
legacy-archive:legacy/rasputin/reader.py
legacy-archive:legacy/rasputin/solar_position.h
legacy-archive:legacy/rasputin/triangulate_dem.h
legacy-archive:legacy/tests/test_polygons.py
```

The only clip of land-cover features is
`legacy/rasputin/gml_repository.py:160-178` (`read`). It parses every GML
feature and tests each whole polygon with `polygon.intersects(domain.polygon)`
(unprepared). Each polygon that meets the domain is then clipped as an area,
`polygon.intersection(domain.polygon)`. `legacy/rasputin/geometry.py` (line 243) wraps the same area
intersection, and `legacy/rasputin/reader.py` (line 362) clips raster indices with `np.clip`, which
is unrelated. **Nothing is carried.** The whole-feature test before the
detailed work is the same idea as section 3.2, but here it runs against the
prepared region, after the checks that can refuse a feature. The area clip is
what 16b's R6 replaced with linework.

## 2. What is slow, and why

From `@perf`'s profile (section 2 of the bottlenecks README; battery, scratch
copies of the code, identical output). The CORINE file
(`corine2018_dtm10_utm33.gpkg`) is in EPSG:25833, the DEM's CRS, so every
feature takes the *moved first* branch of `_Tally._take`: no move, then
`pre_clip(geometry, region)` with no widening, then the intersection with
the domain.

| step | today, s (Numedalslågen / Skiensvassdraget) | why |
|---|---|---|
| `pre_clip`, all features | 3.54 / 2.70 | every edge of every feature read (5.4 M / 5.0 M) becomes a `LineString` and is tested against the region, **unprepared** (section 1). 4,635 / 2,409 of the 7,932 / 5,522 features read miss the region altogether |
| intersection, chains the domain covers | 1.14 / 2.32 | GEOS overlay of each chain against the whole domain (14k / 16k vertices), although the chain lies inside it |
| intersection, chains the domain cuts | 0.87 / 0.46 | the real work; stays |
| the rest | 0.22 / 0.13 | of it, `self.domain.polygon.point_on_surface()` is recomputed for every feature left with no line (6,321 / 2,894 of them), about 0.18 s on Numedalslågen in the prototype's profile |

The time follows the number of vertices read, which follows the bounding box
of the hull. Numedalslågen is long and thin, and the box is 5.1 times its
area. The read itself is 0.43 / 0.28 s and stays as it is (section 11).

## 3. The blueprint: what changes where

One file of production code, `src_python/tin_engine/feature_input.py`. No
change to `cli.py`, to any signature (`pre_clip`, `open_features`,
`read_source`, `source_region`), to `FeatureSet` or `TerrainFeature`, to the
GeoPackage reader, and no C++, binding or new dependency (numpy and shapely,
as today).

```
cli._open_features ──> open_features(request, domain, dem_crs)          unchanged
                         └─ _Tally.__init__        + self.inner: the domain's point_on_surface, once
                         └─ _Tally._take           region prepared once; per feature:
                              ├─ refusal checks (geometry type, class map, code, finite)   unchanged, same order
                              ├─ moved-first branch: region.intersects(_linework(moved)) ? pre_clip(...) : []   NEW whole-feature test
                              ├─ per kept chain: _untouched(chain, domain) ? (chain,) : parts of intersection   NEW
                              └─ covering polygon test against self.inner                     (was recomputed per feature)
                         pre_clip(geometry, region, widening=None)          same contract
                              ├─ endpoints: intersects_xy(region, x, y); an edge with an end inside is kept   NEW
                              ├─ the other edges: intersects(region, segments), prepared region first           arguments swapped
                              └─ widening, on the edges still dropped                                         unchanged
                         _runs(xy, keep, closed)                       + returns [] at once when no edge is kept
                         _untouched(line, domain) -> bool                NEW
                         _linework(geometry) -> BaseGeometry             NEW
```

### 3.1 `pre_clip` (today `src_python/tin_engine/feature_input.py@a7154ec:199-209`)

For each ring or line, with `xy` its vertices:

1. `inside = shapely.intersects_xy(region, xy[:, 0], xy[:, 1])`, with `region`
   prepared (it is: `pre_clip` already calls `shapely.prepare(region)`).
2. `keep = inside[:-1] | inside[1:]`: an edge with an end in the region is
   kept.
3. `rest = np.flatnonzero(~keep)`; segments are built **only** for `rest`, and
   `keep[rest] = shapely.intersects(region, segments)`, with the prepared
   region as the first argument. This works when `rest` is empty: the
   prototype's fixtures include rings wholly inside.
4. The widening loop (geographic sources only) runs over the edges in `rest`
   still not kept, exactly the set it runs over today, with the same
   `shapely.distance(segment, region) <= w`.

The docstring and the return value are unchanged.

### 3.2 The whole-feature test (today `src_python/tin_engine/feature_input.py@a7154ec:333-334`)

In the *moved first* branch only:
`kept = list(pre_clip(moved, region)) if region.intersects(_linework(moved)) else []`,
with a new pure helper

- `_linework(geometry) -> BaseGeometry`: for a `Polygon` or `MultiPolygon`,
  `geometry.boundary` (every ring, each hole's included, as lines); for a
  `LineString` or `MultiLineString`, the geometry itself.

**Not `region.intersects(moved)`** (round 1's design): for an invalid polygon
whose hole lies outside its shell, GEOS's envelope test (section 1) drops the
hole, and `pre_clip` keeps that hole's edges today (16b's R1 makes such a
polygon linework). **Not `moved.boundary` for every type** (the review's
suggestion, taken literally): a line's boundary is its two end points, so a
line crossing the region with both ends outside would be dropped, and a
closed line's boundary is empty. Both were planted in the
prototype: the existing suites pass under each, and the probe's new cases
(section 6) catch each.

`region` is prepared once, right after `source_region(...)` in `_take`. The
test sits where `pre_clip` is called today. That is after the geometry-type
check, the class-map check (which can refuse a feature), the code check and
the finite check. So a feature that misses the region is still refused, or
still counted `outside`, exactly as today.

**Not** in the geographic branch (`bound is not None`). There, an edge outside
the region is kept when it lies within its own widening, so a feature that
misses the region can still keep edges. That branch serves a geographic
source over a Transverse Mercator DEM (a GeoJSON in EPSG:4326). Neither
catchment exercises it, and it keeps today's code apart from 3.1.

### 3.3 `_runs` (today `src_python/tin_engine/feature_input.py@a7154ec:218-219`)

After the `keep.all()` return, `if not keep.any(): return []`. Today the loop
returns `[]` in that case too, after walking every vertex in Python.

### 3.4 `_untouched(line, domain) -> bool` (new) and its use (today `src_python/tin_engine/feature_input.py@a7154ec:336-341`)

`True` exactly when all four hold, tested in this order (cheapest and most
selective first):

1. the line has at least two vertices;
2. no two consecutive vertices are equal (`(xy[1:] == xy[:-1]).all(axis=1)`
   is all `False`);
3. `shapely.contains_properly(domain, line)`: the line lies in the domain's
   **interior**, touching neither the outline nor a hole's boundary. `domain`
   comes first, so its prepared form is used (`_Tally.__init__` already
   prepares `domain.polygon`);
4. `shapely.is_simple(line)`: no self-crossing and no self-touching. A closed
   ring that touches itself only at its start is simple.

In `_take`, each kept chain gives `(chain,)` when `_untouched(chain,
self.domain.polygon)` holds, else the parts of `shapely.intersection(chain,
self.domain.polygon)` as today. The filter after it
(`isinstance(piece, LineString) and piece.length > 0`) is unchanged and is
applied to both.

`_untouched` is pure: it reads its arguments and changes neither of them.

### 3.5 `self.inner` (today `src_python/tin_engine/feature_input.py@a7154ec:346`)

`_Tally.__init__` computes `domain.polygon.point_on_surface()` once, after
`shapely.prepare(domain.polygon)`. The covering-polygon test uses it instead
of recomputing it for every feature left with no line. The point depends only
on the domain polygon, which does not change during the run, so it is the same
point.

**Why not `@perf`'s rule** (skip the intersection whenever the prepared domain
`covers` the chain). It is not byte-identical. GEOS's intersection splits a
covered line wherever the line touches the domain's boundary, at a domain
vertex lying on the line, and at the line's own crossings and self-touches. It
also drops repeated vertices and returns nothing for a zero-length line
(section 4.3's table). On the two catchments no chain happened to do any of
these. In the fixtures two do: with the `covers` rule the probe's lines differ
for `test_a_self_intersecting_ring_is_linework_not_refused` and
`test_a_feature_edge_along_the_boundary_is_kept`, **and the existing suite
still passes** (both tests pin lengths, not pieces). The four conditions of
3.4 skip the same chains on the two catchments and cost nothing measurable
beyond `covers` (section 8).

## 4. Why the output cannot change

### 4.1 `pre_clip`: the same `keep` for every edge

- *An edge with an end in the region.* Today it is kept when
  `intersects(segment, region)` is true. A closed segment containing a point of
  the closed region meets the region, and GEOS decides both tests with exact
  orientation predicates. The point test is the region's point-in-area
  locator, with the boundary counting as inside. The segment test is the same
  locator on the segment's first point, then an exact segment-against-boundary
  test (section 1), which also finds a second endpoint lying on the boundary.
  So the two agree.
- *An edge with both ends outside.* It gets the segment test as today, with
  the arguments swapped so that the prepared path runs. `intersects` is
  symmetric, and both paths are exact for a polygon and a segment. The probe
  (section 6) found the same keep on every fixture and on both catchments.
- *The widening.* It runs on the same set of edges with the same arguments.

### 4.2 The whole-feature test: no feature loses an edge it keeps today

Let `L = _linework(moved)`. Its segments are exactly the edges `pre_clip`
tests: every ring of every polygon part (shell and holes) or every line part,
which `pre_clip` walks with `get_parts`, `exterior` and `interiors`. Its
components are those rings and lines, each a connected set. The claim: if the
prepared `region.intersects(L)` is false, no edge of `L` meets the region.

Take an edge `e` that meets the region `R`, in a component `C`. GEOS's steps
(section 1) find it:

1. *Envelopes.* `e` has a point in `R`, and that point lies in `L`'s envelope
   and in `R`'s, so the envelopes meet and the test goes on. This is the step
   round 1's `region.intersects(moved)` failed: there the envelope was the
   shell's, and `e` (a hole's edge) could lie outside it.
2. *Segments.* If `C` meets `R`'s boundary anywhere, the exact
   segment-against-boundary test over all of `L`'s segments finds it.
3. *One point per component.* Otherwise `C` is connected, misses `R`'s
   boundary and has a point (on `e`) in `R`'s interior, so all of `C` lies in
   the interior, and the located first point of `C` is inside.

The step that does not apply: `L` has dimension 1, so the "region inside the
test geometry" check never runs; a line set cannot contain an area.

So when the test is false, every `keep` in `pre_clip` would be false and it
would return `()` (no widening in this branch). The new code gives `kept =
[]`, the same as `list(())`. From there both take the same path: no lines,
then the covering-polygon test with the same `polygon` and `self.inner`
(3.5, the same point). That includes a polygon that surrounds the region
without an edge in it: its linework misses the region, today's `pre_clip`
keeps nothing, and both go on to the covering test
(`test_a_polygon_around_the_domain_with_no_edge_kept_is_dropped_and_counted`).

Round 1's text said GEOS locates "one point per ring". That is true of the
component filter, but it never mattered: the envelope test runs before it,
and a polygon's envelope is its shell's. The reviewer's case (region
`box(0, 0, 10, 10)`, shell (20..30)², hole (1..2)²) was rerun on `a7154ec`'s
install (shapely 2.1.2, GEOS 3.13.1): prepared `region.intersects(g)` is
`False`, `region.intersects(g.boundary)` is `True`, and `pre_clip` keeps the
hole ring. A variant whose shell's envelope overlaps the region's but whose
shell misses it gives `True` even with `region.intersects(g)`, through step 3
on the hole's first point: the envelope is what drops the hole.

### 4.3 `_untouched`: the intersection would have returned the chain itself

GEOS's overlay (OverlayNG, the default since GEOS 3.9) nodes the line against
the domain's boundary and against itself, removes repeated points, and builds
the result from the noded edges. When the line lies in the domain's interior,
is simple and repeats no vertex, there is no node to insert and no point to
remove. One line comes back, in the original direction, starting at its
original first vertex. That account is recalled, not reread in the source. So
the claim rests on these checks, run for this design with shapely 2.1.2 and
GEOS 3.13.1:

- `docs/increments/30b-probes/geos_cases.py`, 14 hand-made cases at UTM
  magnitudes. `covers` is true in all 14. The line comes back unchanged in 3:
  inside, a closed ring inside, a collinear middle vertex. It is split or
  changed in 11: a repeated vertex; a closed ring touching the outline at its
  start or mid-ring; a domain vertex on the line; a line along the outline;
  along, then inside; a self-crossing; a self-touch; doubling back; touching a
  hole; zero length. The prototype's `_untouched` is true in exactly the 3.
- A seeded property run (20,000 random lines near a 400-vertex domain with a
  hole, with lattice-rounded, closed, repeated-vertex and domain-vertex
  variants): `_untouched` was true for 4,549, and the intersection returned
  every one of them unchanged. With `covers` in its place it was true for
  14,789, and 10,079 of those came back changed, so the run can fail. The
  red suite's R2 (section 7) is this run, made small.
- The probe's gate (section 6) on the fixtures and both catchments.

GEOS does not promise the third bullet's behaviour in a specification, so a
GEOS upgrade could change it. That is what R2 and the gate are for (section
10).

### 4.4 The counts

`outside` and `empty` are counted on the same features as today: no feature is
read or skipped differently, and the refusal checks run first (4.2).
`clipped` is still `not all(domain.covers(g) for g in kept)` over the same
`kept`.

## 5. The blueprint against the assessment framework

- *Data and execution:* no configuration involved. `MARGIN` and `DENSIFY` are
  unchanged.
- *State:* `_untouched` is pure. Two side effects exist today and are kept,
  not added: `_Tally.__init__` prepares the caller's `domain.polygon`, and
  `pre_clip` prepares its `region` argument. `region` is now prepared in
  `_take` before the loop; it is `_take`'s own object.
- *Dependency gravity:* none added.
- *Async-readiness:* unchanged. Still blocking GEOS and sqlite3, run by an
  async caller in `asyncio.to_thread` as the module docstring says.

## 6. The byte-identical gate

`docs/increments/30b-probes/clip_bytes.py` and its base output
`docs/increments/30b-probes/base_a7154ec.txt` are committed with this design.
Its docstring says how to run it. Both modes are the gate:

- **`fixtures`** runs the six suites that reach the features code
  (`test_feature_input.py`, `test_cli_mesh_features.py`,
  `test_cli_mesh_multi_features.py`, `test_cli_mesh_landcover.py`,
  `test_cli_mesh_edge_strip.py`, `test_cli_mesh_plain_output.py`) in-process.
  It records every call of `open_features` (from the module or from the CLI)
  and every direct call of `pre_clip`, keyed by test id and call number. A
  call's line holds the four counts, the number of lines and vertices, and a
  hash of every feature's fid, mask, code, lines (WKB, in order) and polygon
  (WKB). Then it runs `open_features` on three hand-made sources no suite at
  the base has (`cases`, keyed `probe::`, added after design review round 1):
  an invalid polygon whose hole lies in the domain and whose shell's envelope
  misses the region's; the same with a shell whose envelope overlaps the
  region's; and a line crossing the domain with both ends outside the region.
  168 calls at the base (165 from the suites, 3 cases); the base keeps one
  feature in each case.
- **`mesh`** runs `rasputin mesh` on Numedalslågen and Skiensvassdraget as
  `@perf`'s profile did, and records the same line for the feature set plus a
  hash of the whole `.vtk` written.

The gate: every line starting `fixture ` or `mesh ` is equal to the base
file's.

**The base.** It was run with `a7154ec`'s own install (`@perf`'s
`worktree-bottlenecks/.venv`, a non-editable install of that commit; pytest,
pytest-asyncio and hypothesis put on the path from a scratch directory). The
`fixtures` mode was run twice and the `mesh` mode twice, with identical lines.
After round 1 added the cases, the base was rerun with the same install:
`fixtures` twice (identical lines; the 165 suite lines equal round 1's), `mesh`
once (both lines equal round 1's), and `base_a7154ec.txt` rewritten from it.
Numedalslågen: 1,611 kept, 271 cut, 6,321 outside, 0 empty.
Skiensvassdraget: 2,628 kept, 216 cut, 2,894 outside, 0 empty.

**It can fail.** These are plants applied to the prototype, with the base
lines unchanged:

| plant | base fixture lines not matched (of 168) | pytest | catchment lines that differ (of 2) |
|---|---|---|---|
| `_untouched` replaced by `domain.covers(line)` (`@perf`'s rule) | 2 | **passes** | 0 |
| `_runs` returns `[]` for a single kept edge | 3 | fails | 0 |
| an edge kept only when **both** ends are inside, no segment test | 23 (one a `probe::` line) | fails | 1 |
| whole-feature test on `moved` (round 1's 3.2) | 1 (`probe::`, shell's envelope misses) | **passes** | 0 |
| whole-feature test on `moved.boundary` for lines too | 1 (`probe::`, the line) | **passes** | 0 |

Rerun after round 1 on the 168-line base, with the plants rewritten in a new
scratch prototype; round 1's third row counted 16 of 165 with its own plant.

**What the gate cannot see**, and so the red tests pin (section 7):

- the `covers` plant on real data: neither catchment has a chain that touches
  the outline, crosses itself or repeats a vertex, so only the fixtures and R1
  and R2 guard 3.4;
- the two whole-feature plants on real data: the CORINE layer has polygons
  only, none with a hole outside its shell, so only the `probe::` cases and
  P5 guard 3.2;
- the order of the refusal checks and the whole-feature test: no fixture has a
  refused feature lying outside the region.

**The prototype passes the gate**: identical `fixture` lines (168) and
identical `mesh` lines (both catchments, `.vtk` hashes included).

## 7. Tests for `@tester` (the red suite)

In `tests/python/test_feature_input.py`. Behaviour does not change, so the red
tests are the contract of the one new helper. The pins guard what the rewrite
could break and are green today, which is their point; they are committed with
the red ones, before any code. No mutation round.

**Red (fail today: `_untouched` does not exist).**

- **R1. `_untouched`, case by case.** These are the 14 cases of
  `geos_cases.py`, written in the test. It is true for: a line inside, a closed
  ring inside, a collinear middle vertex. It is false for: a repeated vertex; a
  closed ring touching the outline at its start, and one touching it mid-ring;
  a domain vertex on the line; a line along the outline; along, then inside; a
  self-crossing; a self-touch at a vertex; doubling back; a line touching a
  hole's boundary; a zero-length line. Also false for a line partly outside and
  a line wholly outside. Use a domain with an extra vertex on one side and a
  domain with a hole, at the existing `X0, Y0` offset. Coordinates must be
  exactly representable.
- **R2. `_untouched` against GEOS.** On a seeded random set (a few hundred
  short lines near a domain with a hole, some rounded to a lattice, some
  closed, some with a repeated vertex, some through a domain vertex): wherever
  `_untouched` is true, `shapely.intersection(line, domain)` gives one
  `LineString` with the same WKB as the line. Assert that it was true for at
  least a tenth of the set, and for at least one closed ring of each
  orientation (counter-clockwise and clockwise) and one open line, so the test
  cannot pass vacuously and covers the shapes a reoriented or restarted ring
  would break. The oracle lives in the test. R2 is the only test that runs in
  CI and would catch a GEOS change to how a covered line comes back
  (section 10).

**Pins (green today, must stay green).**

- **P1. The pieces of a chain the domain covers but touches.** Through
  `open_features` on a GeoJSON source, for: a ring inside touching the outline
  at one vertex; a line running along the outline; a self-crossing ring inside
  (the bowtie of `test_a_self_intersecting_ring_is_linework_not_refused`); a
  ring inside with a repeated vertex. The feature's `lines` equal, piece by
  piece and vertex by vertex, the parts of `shapely.intersection(chain,
  domain)` that are lines of positive length, where `chain` is the source
  ring or line. The test computes that oracle. This is the pin that the
  `covers` plant broke while the suite passed.
- **P2. Refusal before the whole-feature test.** Under a map whose
  `otherwise` is `"refuse"`, a feature far outside the region with a value
  the map lacks is refused (`FeatureError` naming the feature), as one inside
  is.
- **P3. Outside counts with no edge kept.** A source with three features
  that miss the region and one inside gives `outside == 3` and one feature.
- **P4. `pre_clip` decides by the edge, not by its ends.** On the existing
  `REGION`: an edge with both ends outside that crosses the region is kept; an
  edge with one end exactly on the region's boundary and the other outside is
  kept; an edge with both ends outside that passes by a corner without
  touching is dropped. (`test_an_edge_touching_the_region_at_one_point_is_kept`
  already pins the corner touch.)
- **P5. The whole-feature test keeps what `pre_clip` keeps.** Through
  `open_features` on a GeoJSON source and the `BOX` domain, the three
  `probe::` cases of `clip_bytes.py` (section 6): an invalid polygon with
  shell `square(2_000, 2_000, 2_100, 2_100)` and hole `INNER`'s ring; the same
  hole under an L-shaped shell whose envelope overlaps the region's but which
  misses the region; and a line from `at(-500, 150)` to `at(800, 150)`. Each
  gives one feature and `outside == 0`; the two polygons' `lines` are the
  hole's ring, vertex by vertex, and the line's `lines` are its piece inside
  the domain. Green today. Round 1's whole-feature test fails the first, and
  `moved.boundary` for lines fails the third, while the existing suites pass
  under both (section 6).

The existing tests (`TestClip`, `TestPreClip`, `TestPreClipKeepsWholeEdges`,
`TestPreClipThroughOpenFeatures`, `TestLongEdgeWidening`) stay as they are.

## 8. Net production lines, and the prototype

The prototype is the installed `a7154ec` package copied to a scratch directory
and edited there (no C++ build; removed after measuring).
`python3 tools/count_loc.py` between its base and the edit gave **+20 net**
(27 added, 7 removed), all in `feature_input.py`. After round 1, a second
prototype with 3.1 to 3.5 and `_linework`, typed and without the
plant switches, committed in a scratch repository over `a7154ec`'s
`feature_input.py`, gave **+17 net** (24 added, 7 removed); it passes the
gate's `fixtures` mode (168 lines equal). The PR's docstrings and any
spelling-out add a few. Estimate for the PR: **+20 to +30**, unchanged.

Timings of `open_features`' `clip_seconds` (the `features clip` row), on
**AC power** (`pmset -g batt` before and after: AC), base and prototype
alternated, two rounds of three runs each, medians. Round 1 measured the
whole-feature test on `moved`; after round 1 the same prototype was timed
with that test and with the linework test (3.2) alternated in each run, the
`rasputin mesh` run stopped right after `open_features`:

| | Numedalslågen s | Skiensvassdraget s |
|---|---|---|
| `a7154ec`, round 1 | 5.93 / 5.90 | 5.84 / 5.72 |
| prototype, test on `moved`, round 1 | 1.82 / 1.79 | 1.31 / 1.31 |
| `a7154ec`, after round 1 | 5.83 / 5.88 | 5.65 / 5.68 |
| prototype, test on `moved`, after round 1 | 1.81 / 1.81 | 1.31 / 1.30 |
| prototype, linework test (this design) | 1.86 / 1.85 | 1.35 / 1.33 |
| speed-up, this design | 3.2 | 4.2 |

Building the boundary costs about 0.04 s on Numedalslågen and 0.03 s on
Skiensvassdraget (7,932 and 5,522 features read).

A cProfile of a first prototype on Numedalslågen (clip 2.41 s; it had sections 3.1
to 3.4, except that the segment test kept today's argument order, and it
still recomputed the domain's `point_on_surface` per feature) showed per run: the intersection of
the chains the domain cuts, about 0.84 s (2,912 chains, the count `@perf`
measured as cut); `intersects` calls, about 0.49 s; `point_on_surface`, about
0.18 s. Swapping the arguments (3.1, step 3) and `self.inner` (3.5) took it to
round 1's 1.8 s above. What is left is mostly the intersection of the cut chains. `@perf`'s scratch variant "A + both" measured 2.11 / 1.48 s on battery, so
the prototype is at least as fast as that variant, with the stricter rule of
3.4.

## 9. Acceptance (`@perf`)

`tools/bench.py` and the thread sweep are not required: nothing under
`include/terrain/refinement/` or `include/terrain/mesh/` changes. `@perf`
runs, with the branch's own install and a base install of the merge base,
**the same power state for both, recorded with `pmset -g batt` before and
after each run**:

1. **Byte-identical**: the probe's `mesh` mode on both catchments, every line
   equal to `base_a7154ec.txt` (feature set and `.vtk` hash), and its
   `fixtures` mode, every line equal. If 30a has merged into the branch's
   base by then, the `.vtk` hashes still hold, because 30a's own gate is
   byte-identical meshes.
2. **Time**: `rasputin mesh --stats` on both catchments, three repeats each,
   base and branch alternated, as
   `docs/benchmarks/2026-10-06/bottlenecks/scripts/stats.sh` does. Read the
   `features clip` row's median. **Pass:** the branch's median is at most 0.4
   of the base's on both catchments. This was checked only at these two
   inputs (7,932 and 5,522 features read, 5.4 M and 5.0 M vertices); the
   prototype gave 0.32 and 0.24.
3. The `features read` row and the whole-run total, recorded, not gated. The
   read is unchanged, so its row should not move.

Evidence goes under `docs/benchmarks/<date>/30b-clip/`.

## 10. Risks

- **A GEOS upgrade that changes how a covered line comes back.** If GEOS
  started returning, say, a reoriented ring for a line in the interior, the
  skipped chains would differ from what the intersection gives. Only R2
  guards this in CI. P1 does not: its oracle is GEOS's own intersection, so
  it moves with GEOS. The gate does, but it is run by hand, against a base
  recorded with shapely 2.1.2 and GEOS 3.13.1. And nothing holds GEOS still:
  `pyproject.toml` asks for `shapely>=2.0` with no upper bound, so CI takes
  the newest shapely wheel and the GEOS inside it. The mesh would still be
  valid either way; only the byte-identity is at stake.
- **Invalid feature polygons** are linework, not refused (16b R1). Rings that
  cross themselves are not simple, so `_untouched` sends them through the
  intersection as today (pinned by P1). A hole lying outside its shell is
  why the whole-feature test reads the linework (3.2, 4.2; pinned by P5).
- **The geographic branch** keeps its own path (3.2). A future source that
  is geographic and large would not get the whole-feature speed-up. Not
  measured, and not needed for CORINE, which is projected in both copies Ola
  uses (EPSG:25833 here, EPSG:3035 in the full file).
- **Memory.** It falls a little: segments are built only for edges with both
  ends outside. Not measured.

## 11. Not in scope

- **Reading fewer features.** Testing each R-tree hit's box against the hull
  before decoding its geometry, or querying the index with several boxes along
  a long, thin catchment, would skip decoding the features that miss the hull
  (4,635 of 7,932 on Numedalslågen). It would also change the `features:`
  stderr line, because those features would no longer be counted `outside`
  (6,321 today on Numedalslågen, of which 4,635 miss the hull).
  The read is 0.43 / 0.28 s in all, so the saving is a fraction of that. Not
  measured. It is left to Ola as question 1.
- **Clipping the chains the domain cuts** (about 0.84 s left on Numedalslågen)
  any other way than GEOS's intersection. That would change the vertices GEOS
  puts on the outline.
- **Keeping a covered chain whole when it touches the outline** (`@perf`'s
  `covers` rule). It would change the output (3.5), and so it would be a
  behaviour change with its own design, not a speed-up.
- 30a (land cover) and 30c (reading the DEM).

## 12. Questions for Ola

1. **Should a later change be allowed to change the "outside the domain"
   count on the `features:` line, so that features far from the catchment are
   never decoded?** Today that count includes every feature whose box meets
   the search box. On Numedalslågen that is 6,321, of which 4,635 miss the
   hull and are dropped without a test. Skipping them would save
   part of the 0.4 s read. That gain is small and not measured, and the number
   would then mean "near the domain but outside it". **Default: no; the count
   stays as it is, and this is not pursued.**

## 13. ROADMAP

Row 30 is on 30a's branch (`worktree-landcover-speed`) and not yet on master
or here. So this branch does not add or edit row 30. It adds its own line,
`| 30b |`, placed after the GeoPackage output row (`| — |`), the row that row
30 is inserted before on 30a's branch. The two insertions are separated by
that unchanged row, so either branch merges onto the other cleanly. This was
checked with `git merge-tree --write-tree` on the two branch heads: no
conflict. Whichever of 30a and 30b lands second folds its status into row 30
and removes the separate line, in its own PR.

**The fold also rewrites row 30's description of 30b.** On
`worktree-landcover-speed`, row 30 describes 30b as "skip features that miss
the hull, cheap tests before the exact ones, no clip of a line the domain
already covers". This design rejects the first (section 11: the reading is
unchanged and features outside are still counted) and the third (3.5: a line
the domain covers but touches is still clipped). The fold makes it match this
design, in words like: "30b the CORINE clip (point tests before the
segment test, a feature whose lines miss the search region skipped
whole, no clip of a line lying in the domain's interior)". This branch cannot
write in that worktree, so the rewrite happens at whichever merge lands
second, as part of the fold.

## Review

**Design review, round 1, 2026-10-06.** Range `a7154ec..781451d`. Verdict: CHANGES REQUESTED: the whole-feature test (3.2) drops the hole edges of an invalid polygon whose shell misses the region (B1); "one point per ring" in 4.2 is wrong (B2); row 30's 30b description on `worktree-landcover-speed` contradicts this design (B3). LOC: 0 (design only); estimate +20 to +30. Not pushed; no CI.
