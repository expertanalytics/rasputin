# Increment 17 — `rasputin mesh --stats`: sizes, quality and timings as Markdown

Status: **implemented, in review.** Red `5843d10`; green `a5d2674` (C++),
`9d6ba5d` (`stats.py`), `4a1d806` (`cli.py`). Measured 421 added / 73 removed,
net 348 production lines against ~210 (see "What landed"). Written by `@architect` before `@tester`, per
`docs/increments/README.md` step 1, on branch `increment17-mesh-stats` off
`increment16-domain-polygon` (PR #94, merged). The user chose C1 (a) stdout
and C2 (a) count carving on 2026-09-26 (section "Ruled by the user").

**Closes.** The user's request (2026-09-26): "when running the scripts and
generating geometries, we should also have the option of --stats or similar to
show simple statistics like what you have shown me. Sizes, quality metrics,
timings, etc. Output as a markdown." After this increment,

```sh
rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif \
    --domain quarter.geojson --tolerance 1 --out quarter.vtk --stats quarter.md
```

writes the mesh as today and a Markdown report next to it: the sizes, the
triangle quality, the refinement counters and a per-phase time table. It also
gives the next increment, the performance look at NoData carving (5 054 rounds,
109 s on an outline-only tile; `16-domain-polygon.md` R5), a timing table it can
trust without wrapping functions by hand as 14b's M4 and M5 did.

**Not closed.** Per-round series (time and inserts per round), which the carving
investigation will want, but as a throwaway probe of its own (see "Not in
scope"). `--stats` on `draw`. Machine-readable output (JSON).

## Rulings

### R1. Where the stats are computed: a pure module, `tin_engine/stats.py`

A new module with no `_core` import and no typer. It takes numpy arrays, plain
numbers and a phase list, and returns Markdown. It is testable without the
extension, which is what `@tester`'s hand-built meshes need.

Its surface:

- `PhaseClock`: an append-only recorder. `with clock.phase("decode"):` times a
  block; `clock.add(name, seconds)` records a time measured elsewhere (the C++
  sub-phases of R5). A name used twice accumulates, so two files' encodes are
  one "encode" row. `clock.phases()` returns the rows in first-seen order.
- `Quality` (frozen dataclass) and `quality(xy, triangles) -> Quality`. R3.
- `Refinement` (frozen dataclass): tolerance, achieved max error, rounds,
  inserted, flips, uncovered, and optionally carved (C2). Built by `cli.py` from
  `RefineOutcome`, so `stats.py` never sees a `_core` type.
- `Sizes` (frozen dataclass): the rows of R4's "Sizes" table.
- `Report` (frozen dataclass): the command line, `Sizes`, `Quality`,
  `Refinement | None`, the phases, the total, and the stats' own time.
- `render(report) -> str`: the Markdown. Pure; no clock, no I/O.

Frozen dataclasses, not Pydantic: these are results, not configuration, and
`elevation.Trimmed` already sets that precedent.

`cli.py` owns the wiring and nothing else: it creates one `PhaseClock` per
`mesh` call, passes it down, and, only if `--stats` was given, builds the
`Report` and writes `render(report)`. `_engine` gains an optional `clock`
parameter (a fresh, discarded one when absent, so `draw` is untouched).
`_dem_mesh` takes the clock, and its 4-tuple return becomes a small frozen
dataclass that also carries `Sizes` inputs and `Refinement | None`: a fifth and
sixth tuple element is where that signature stops being readable.

### R2. Where the report goes: a file, `--stats PATH`; `-` means stdout

`mesh` prints the written path on stdout, one line per file, and tests and
scripts read it. `--stats` must not change that when absent, and a flag that
bare-writes Markdown to stdout would change it whenever present.

- `--stats PATH` writes the report to `PATH`. It goes through `_destination`, so
  `--out-parent` applies, and it is refused if it resolves to `--out` or
  `--out-edges` (the same overwrite trap `--out-edges` already guards). Its
  path is echoed on stdout after the mesh's, like every other file written.
- `--stats -` writes the report to stdout, after the path line(s). This is
  opt-in, so the stdout contract holds for every caller who did not ask. It is
  what `rasputin mesh ... --stats - | glow` or a quick look in the terminal
  wants.
- No suffix is enforced; `.md` is the recommendation in the help text.
- A refused run (any `BadParameter`) writes no report, as it writes no mesh.
- Today's one-line stderr summary with `--tolerance` is unchanged
  (`test_cli_mesh_refine.py` reads it).

Why not stderr: stderr already carries that summary line and every refusal, so
`2> report.md` would capture refusals as "the report". A file is what the user
asked for ("output as a markdown"), and `-` keeps the terminal case one flag
long. C1 offers stderr instead of stdout for `-`.

### R3. Quality metrics: 2-D, in world x/y, from the output arrays

Computed in numpy from the **written** vertices and triangles, after the trim,
so the numbers describe the file.

- **Minimum angle per triangle, in degrees, in the x/y plane of the world
  coordinates.** 2-D, not on the 3-D surface. Reasons: the Delaunay property
  the pipeline guarantees is a plan-view property (the refinement frame
  `(col·dx, −row·dy)` is world x/y up to a translation, so the angles are the
  same); plan-view angles are what 14b and 16 measured, so the numbers stay
  comparable; and surface angles on steep slopes mix mesh quality with terrain
  steepness. The report's heading says "plan view (x, y)" so nobody reads it as
  3-D.
- Each angle is `atan2(|e × f|, e · f)` of the two edges at the corner. Not
  `acos` of a normalised dot product, which loses all precision below about
  1e-4 rad, exactly where the worst angles live (16's T-real worst is 0.0117°).
  A zero-area triangle reports 0°.
- Reported: median, share under 1°, share under 10°, worst.
- **Vertex degree** = triangles incident on the vertex (`np.bincount` over the
  triangle array), as `16-domain-polygon.md` T-real defines it. Boundary
  vertices count too; that is why the median is 6 and not 6 in the interior
  only. Reported: median, p99, max, count ≥ 12, count ≥ 20.
- Percentiles use `np.percentile(..., method="inverted_cdf")`, so p99 is an
  observed degree (an integer), not an interpolation between two. The median of
  angles uses `np.median`.

### R4. The Markdown layout

Four sections, in this order: a header, **Sizes**, **Quality**, **Refinement**
(only with `--tolerance`), **Timings**. Tables throughout. Integers are plain
(no thousands separator, so a report can be grepped and diffed); seconds to
3 decimals; angles to 2 decimals except "worst", which is `:.3g` so 0.0117°
survives; shares as percent to 2 decimals; errors through `cli._exact` (moved
to `stats.py`), so an achieved error never prints above the tolerance it met.

Rows that do not apply are omitted, not printed empty: a fixture has no DEM
grid, a stride run no domain, a run without `--tolerance` no Refinement section
and no refine rows.

The sample is the quarter-circle run at 1 m. Triangles, rounds, flips, the
achieved error, degree median/p99/max and the angle median, < 1° share and
worst are 16's T-real figures; every other number (vertex and edge counts, the
file size, < 10°, the ≥ 12 and ≥ 20 counts, and all times) is illustrative, not
measured.

```markdown
# rasputin mesh — statistics

`rasputin mesh --dem 7908_3_10m_z33.tif --domain quarter.geojson --tolerance 1 --out quarter.vtk --stats quarter.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 5051 × 5051 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 536 |
| start triangles | 534 |
| output vertices | 214210 |
| output triangles | 427779 |
| constraint edges | 536 |
| vertices without data dropped | 0 |
| quarter.vtk | 31.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.03 % | 2.41 % | 0.0117° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 43 | 118 | 9 |

## Refinement

| tolerance | achieved max error | rounds | inserted | flips | uncovered |
|---|---|---|---|---|---|
| 1 m | 0.999998 m | 40 | 213674 | 452067 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| decode | 0.061 | 1.5 % |
| domain read | 0.004 | 0.1 % |
| start mesh: build | 0.001 | 0.0 % |
| start mesh: node | 0.001 | 0.0 % |
| start mesh: triangulate | 0.001 | 0.0 % |
| start mesh: constraint edges | 0.001 | 0.0 % |
| refine | 0.412 | 10.2 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: scan (parallel) | 0.118 | 2.9 % |
| refine: split + flip (serial) | 0.161 | 4.0 % |
| refine: setup + output | 0.133 | 3.3 % |
| trim | 0.028 | 0.7 % |
| write: encode | 3.471 | 86.1 % |
| write: disk | 0.019 | 0.5 % |
| other | 0.031 | 0.8 % |
| **total** | **4.030** | **100 %** |

Statistics computed in 0.072 s, not included above.
```

The `refine: …` rows are inside `refine` and sum to it; they are not added to
the total a second time (the share column of a sub-row is of the total, and the
table says "sub-rows sum to their parent" in a line under it). "other" is the
total minus the top-level rows, so nothing is hidden: if it grows, a phase is
missing.

### R5. Timings: what each phase covers, and the clock

Python phases use `time.perf_counter_ns()`: monotonic, unaffected by wall-clock
changes, nanosecond API. Each phase is wall time on the calling thread.

| phase | covers exactly |
|---|---|
| decode | open the GeoTIFF, `decode_dem`, close |
| domain read | `read_domain` and `_domain_chains` |
| start mesh: build / node / triangulate | the three calls inside `_engine` |
| start mesh: constraint edges | `_constraint_arrays` on the start mesh |
| sample | `to_core` and `sample` (no-tolerance path only) |
| refine | `to_core` and the `_core.refine` call, GIL released inside |
| trim | `elevation.trim` |
| write: encode | `write_vtk` / `write_ply` producing bytes, all files |
| write: disk | `write_bytes`, all files (mesh files only; the report is not timed) |
| other | total minus the rows above: argument checks, `_destination`, glue |

For a gallery fixture the phases are the start-mesh rows (build, node,
triangulate, constraint edges) and write.

**Trustworthiness, stated in the report and here.** One run, no warm-up: a
first-touch cost (page faults on a fresh allocation, numpy's first call into a
routine) lands in whichever phase pays it. The phase boundaries are the call
boundaries above, so a row can be reproduced by wrapping the same call. The
clock calls themselves cost about 50 ns each, a dozen per run.

**The clock always runs; the metrics do not.** A dozen `perf_counter_ns` calls
cost microseconds, and one code path is simpler than a null clock. Without
`--stats` nothing else happens: no quality pass, no report, no extra output.

### R6. `refine`'s phase split: in, as three `RefineOutcome` fields

In, because the next increment is a `refine` performance investigation and
14b's M5 had to build prototype timers to get exactly this split.

```cpp
double legalise_seconds = 0.0;  // the start mesh's legalise_all
double scan_seconds = 0.0;      // every round's parallel scan, summed
double split_seconds = 0.0;     // every round's serial split + flip phase, summed
```

- Measured with `std::chrono::steady_clock` (monotonic, the standard's clock for
  intervals), **on the calling thread only**: one `now()` before and after
  `legalise_all`, before and after each round's `for_each_chunk`, and before and
  after each round's serial loop. `for_each_chunk` joins every worker before it
  returns (`parallel_util/chunks.hpp`), so the scan interval is the parallel
  region's wall time and no worker ever touches a clock or a shared counter.
  Nothing new is shared between threads, so the TSan job has nothing new to see.
- Cost: four `now()` per round, about 25 ns each on macOS and Linux: 40 rounds
  cost 4 µs; 5 054 rounds cost about 0.5 ms against 109 s. Always on, like
  `rounds` and `flips`, so there is no build flag and no second code path.
- "setup + output" in the report is `refine` minus the three, computed in
  Python: `to_lattice`, the final result sweep, the output assembly, and the
  binding's conversion back to numpy.
- Bound read-only in `bindings/core.cpp`, declared in `_core.pyi`. The output
  mesh is unaffected, so the determinism guarantee (bit-identical for any
  `threads`) is untouched: timings are not part of it and no test may compare
  them.

## Ruled by the user

**The user chose C1 (a) and C2 (a) on 2026-09-26.** The options are kept so the
reasons stay on record.

### C1. Where `--stats -` goes

- (a) **stdout, after the path line(s). Recommended.** Opt-in, pipes cleanly to
  a Markdown viewer, and callers who read stdout for the path are unaffected
  unless they asked for the report.
- (b) stderr. Keeps stdout path-only in every case, but mixes the report with
  the existing summary line and refusals.

### C2. Count carving splits in `RefineOutcome` (`carved`)

- (a) **Yes. Recommended.** One counter in the serial loop (`r.is_void` on a
  split), one binding line, one table column. The carving investigation's first
  question is how many of 5 054 rounds' inserts were carving, and today nothing
  says.
- (b) No; leave it to the investigation's own probe.

## Prior art in `legacy/`

```sh
$ grep -rliE 'min_angle|minimum angle|perf_counter|time\.time|chrono|statistic' legacy/
legacy/bindings.cpp
legacy/rasputin/triangulate_dem.h
legacy/rasputin/avalanche.py
legacy/rasputin/solar_position.h
```

Every hit is a false positive. The C++ ones are `std::chrono` converting
timestamps for solar position and shading (`legacy/bindings.cpp:381, 394`,
`legacy/rasputin/triangulate_dem.h:645`, `legacy/rasputin/solar_position.h`);
`legacy/rasputin/avalanche.py:26` matches `time\.time` inside
`datetime.timedelta`. The legacy tree never reported mesh quality or phase
timings. Nothing is carried across; `@migration-expert` is not needed.

## Files

| file | what |
|---|---|
| `src_python/tin_engine/stats.py` (new) | `PhaseClock`, `Quality`/`quality`, `Refinement`, `Sizes`, `Report`, `render`, `_exact` (moved from `cli.py`). Pure: numpy only; no `_core`, no typer |
| `src_python/tin_engine/cli.py` | `--stats PATH`; its `_destination` and overwrite refusals; the clock through `_engine` and `_dem_mesh`; `_dem_mesh`'s return as a dataclass; building and writing the report |
| `include/terrain/refinement/refine.hpp` | `legalise_seconds`, `scan_seconds`, `split_seconds` (and `carved` under C2 a) |
| `bindings/core.cpp`, `_core.pyi` | the new read-only fields |
| `project_structure.md` | `stats.py` in the tree |
| `ROADMAP.md` | row 17 |

## Tests for `@tester`

**No invariant-critical suite, so no mutation round.** Nothing here changes a
mesh; a wrong statistic is visible on the page.

`stats.py`, without the extension (`tests/python/test_stats.py`):

- **Q1. Known angles.** A hand-built mesh of two triangles: a right isosceles
  (min angle 45°) and a 30-60-90 (min 30°). Median 37.5°, worst 30°, both
  shares 0.
- **Q2. A sliver.** Vertices (0, 0), (1, 0), (0.5, 1e-6): min angle
  `atan(2e-6)` in degrees to 1e-9 relative. This is the case `acos` gets wrong;
  the test must fail against an `acos` implementation (run it against one once
  and say so in the commit).
- **Q3. Known degrees.** A fan of 12 triangles around one centre vertex, and a
  second centre with 20: degree max 20, count ≥ 12 is 2, count ≥ 20 is 1, p99
  an observed integer. Rim vertices have degree 2.
- **Q4. Plan view.** Lifting the vertices to any z (x, y, z input) gives the
  same angles as (x, y): the function reads columns 0 and 1 only.
- **M1. Markdown.** `render` on a fixed `Report` equals a golden string: section
  order, table headers, number formats, omitted rows (no Refinement section when
  `None`, no DEM row for a fixture), "other" = total minus the top-level rows,
  sub-rows not added to the total.
- **P1. `PhaseClock`.** Order is first-seen; a repeated name accumulates; `add`
  records the given seconds. Use a monotonic fake time source injected through
  the constructor, so no test sleeps or asserts on a real duration.

CLI (`tests/python/test_cli_mesh_stats.py`), with the extension:

- **C1. Off by default.** Without `--stats`, stdout is exactly today's path
  line(s) and no `.md` appears in the output directory, for a fixture, a DEM,
  a `--tolerance` and a `--domain` run.
- **C2. `--stats PATH`.** The file exists, stdout is the mesh path then the
  report path, and the mesh file is byte-identical to a run without `--stats`.
- **C3. `--stats -`.** stdout is the path line(s) followed by the report
  (per C1's choice).
- **C4. Sections per run kind.** Fixture: Sizes, Quality, Timings; no
  Refinement, no DEM row. DEM without tolerance: a `sample` row, no Refinement.
  `--tolerance`: Refinement with the same rounds, inserted and flips as the
  stderr line, and the `refine: scan` / `split + flip` / `legalise start` rows.
  `--domain`: the domain row.
- **C5. Refusals.** `--stats` resolving to `--out` or `--out-edges`; outside
  `--out-parent`. A refused mesh run writes no report.
- **C6. Timings are sane, not exact.** Every seconds cell parses as a
  non-negative float; the refine sub-rows sum to at most `refine`. No
  thresholds.

C++ (`test_refinement_refine`, new in this increment): the three seconds fields are
non-negative and sum to no more than a wall clock around the call; with C2 (a),
`carved` is 0 on a DEM without NoData and positive on one with a NoData corner.

## LOC estimate

Counted in `CLAUDE.md` §2's unit.

| file | what | est. |
|---|---|---|
| `stats.py` | clock ~20, dataclasses ~30, quality ~20, render ~50 | ~120 |
| `cli.py` | option and refusals ~15, clock wiring ~25, `_dem_mesh` return ~10, report build and write ~15 | ~65 |
| `refine.hpp` | three fields, six `now()` sites (+ `carved`) | ~15 |
| `bindings/core.cpp`, `_core.pyi` | four read-only fields | ~8 |
| | **total** | **~210** |

On the worst overrun seen so far (+39 %), about 290. Under 700; not split.
Increment 16 overran most in `cli.py` option declarations and refusals, so
that is where to look if this grows.

### What landed

Measured by `@developer` and `@reviewer` (`CLAUDE.md` §2's unit): **+421 / −73,
net 348**, against ~210: 66 % over the estimate and 20 % over the 290
contingency, still under 700.

| file | measured | est. |
|---|---|---|
| `cli.py` | +199 / −72 | ~65 |
| `stats.py` | 189 | ~120 |
| `refine.hpp` | 16 | ~15 |
| `bindings/core.cpp` | +9 / −1 | ~8 (with `.pyi`) |
| `_core.pyi` | 8 | |

The prediction above did not hold. The overrun is not option declarations and
refusals; it is the write section restructured into one encode/target loop
(every destination resolved before any write), the `_DemMesh` dataclass
replacing `_dem_mesh`'s tuple, `_write_report` / `_report_target`, and the
`with clock.phase(...)` re-indentation of existing calls, much of which is
churn rather than new behaviour. `stats.py` is larger mainly in `render`'s
omission rules and formatting.

The real quarter-circle run at 1 m (Ola's example) reported refine at 3.29 s
of 3.39 s, with `refine: scan (parallel)` at 3.11 s; carving 0. This is the
starting point of the performance investigation that follows.

## What was measured

The cost of the quality pass, in the scratchpad (never in the tree): a jittered
490 × 490 grid, 240 100 vertices and 478 242 triangles, the size of the 1 m
runs. The R3 computation (angles by `atan2`, degrees by `bincount`, median,
p99, shares) took **0.072 to 0.089 s** over three runs with the repository's
`.venv` numpy. Against 14b's M4 that is about 15 % of a 0.50 s `--binary` run
and under 5 % of a 1.6 s text run. It is reported on its own line and kept
out of the phase table, so it never distorts the timings it sits beside.

## Acceptance

- The quarter circle at 1 m with `--stats quarter.md` writes a report whose
  Sizes, Quality and Refinement match 16's T-real table, and whose Timings rows
  cover the run with "other" under 5 % of the total.
- A run without `--stats` is byte-identical on stdout and on disk to today's.
- All gates in `CLAUDE.md` §4 green, including the TSan job; CI green.

## Not in scope

- Per-round series (time, inserts, carving per round). The carving
  investigation should build it as a probe; if it earns a place, it becomes a
  `--stats` extension then.
- Surface (3-D) angles, aspect ratios, edge-length histograms.
- JSON output, `draw --stats`, repeated runs or warm-up for timing.
- Any change to what `mesh` writes or prints without `--stats`.
