# tools/bench.py: design

Status: design (@perf, 2026-09-27), reviewed by @reviewer with the green commit. Red suite written (see "Pinned by the red suite"); green: `tools/bench.py` (see "Settled at green"). The rule it serves is
"Acceptance: an increment that touches refine or mesh code" in
`docs/increments/README.md`; its specification is the one-off scripts in
`2026-09-26/` (`bench1m/*.sh`, `run.py`, `quality.py`, `summarize.py`,
`scaling/scale.py`) and that directory's `README.md`.

## CLI

Typer. It is already a runtime dependency and the project's CLI stack, and
`typer.testing.CliRunner` lets the tests call it without a subprocess.
`tools/session_state.py` uses argparse, but it runs under a bare `python3`;
bench.py needs numpy and `tin_engine` anyway, so it runs from the venv.

```
python tools/bench.py run    --label i21 [--tree .] [--dem PATH] [--domain PATH]
                             [--tolerance 1.0] [--threads 1,2,...,20] [--repeats 5]
                             [--mesh-dir DIR] [--no-build] [--baseline DIR]
                             [--out-root docs/benchmarks] [--threshold 5]
                             [--accept-quality]
                             [-- extra mesh args, e.g. --start-min-angle 0]
python tools/bench.py compare NEW_DIR [--baseline DIR] [--out-root DIR]
```

- `run` does both parts of the acceptance run: the 1 m benchmark (the CLI's
  default threads: `cli.py` passes none, so `threads == 0`,
  hardware_concurrency, `include/terrain/parallel_util/chunks.hpp`) and the
  scaling sweep (the refine call forced to each `--threads` value), then
  writes the evidence and the verdict. One command, so the two parts cannot
  come from different builds or power states.
- Defaults are the 2026-09-26 set: DEM `tests/fixtures/dem_archive/7908_3_10m_z33.tif`,
  tolerance 1, repeats 5, threads 1 to 20. Domain: the tile (no `--domain`) and
  `docs/benchmarks/2026-09-26/quarter.geojson`; each `run` measures both unless
  `--domain` names one.
- `--tree` is the source tree measured (default: the repo). Pointing it at a
  worktree of an older commit is how a baseline is re-measured on the other
  power state; `build.sh`'s multi-commit loop becomes one call per tree.
- `compare` re-judges stored evidence without re-running.

## How a run is made

1. **Build** (skipped by `--no-build`, which is then recorded): configure
   `<tree>/build-bench` with `CMAKE_BUILD_TYPE=Release`, build `_core`, and
   assemble `<tree>/build-bench/pkg/tin_engine/`: symlinks to the tree's
   `src_python/tin_engine/*` plus the fresh `.so`. `build-*/` is gitignored, and
   nothing is copied into `src_python/` or the venv, so a bench build can never
   shadow what `pytest` loads. A non-Release cache is refused.
2. **Child process per sample**: `python tools/bench.py _child --pkg DIR
   --threads N -- <mesh args>`. The child drops the editable finder and puts
   `pkg/` first on `sys.path` (run.py's technique), asserts that
   `tin_engine.__file__` and `_core.__file__` are under it, wraps `cli.refine`
   (scale.py's technique: force `threads`, time the call), runs `cli.app()`,
   and prints one JSON line on stderr: `refine_s`, `app_s`, and the result's
   `max_error`, `rounds`, `inserted`, `flips`. The parent times `proc_s`.
3. **Samples**: `--binary` runs to a throwaway path, `--repeats` per
   (domain, thread count), interleaved round-robin over thread counts so a
   thermal or power drift spreads over all of them rather than landing on one.
4. **Quality**: one extra `--ascii` run per domain at the default thread count,
   kept in `--mesh-dir` (default: the scratchpad or `$TMPDIR`, never the repo),
   measured by bench.py's quality function (below).
5. **Power state** read before and after; if they differ the run is marked
   `mixed` and cannot be compared to anything.

## What each run records

`docs/benchmarks/<date>/<label>/`:

- `run.json`, one Pydantic `RunRecord`: tree commit and dirty flag; bench.py's
  blob hash; the child argv; build type, `CMAKE_CXX_FLAGS_RELEASE`, compiler
  and `.so` sha256 (or `no-build`); CPU brand, P/E core counts, memory, macOS,
  Python and numpy; power (`ac` | `battery` | `mixed` | `unknown`, percent,
  raw `pmset -g batt` text); DEM and domain paths and sha256, tolerance, extra
  args; every raw sample, and median/min/max per (domain, threads); quality
  per domain: worst angle, angle median, share under 1 degree, max degree (from
  `tin_engine.stats.quality`, as `--stats` computes them), `max_error <=
  tolerance`, the Delaunay check (checked, ambiguous, violations), and the
  mesh's sha256 from `POINTS` on (the report's "identical mesh" check).
- `raw.tsv`: the samples, in `scaling/raw_*.tsv`'s columns plus power.
- `README.md`, generated: the method (every field above), the 1 m table in
  `REPORT.md`'s layout, the scaling table in the 2026-09-26 README's layout
  with the ceiling (best speed-up over 1 thread, and the thread count where it
  is reached) beside 2.2x, and the verdict. @perf adds the prose by hand below
  a marker line that a rerun preserves.

Meshes stay in `--mesh-dir`; the README says where and how to regenerate.

## Comparison and verdict

A baseline is a `run.json`. `--baseline DIR` names it; otherwise bench.py
takes the newest `docs/benchmarks/*/*/run.json` that is comparable and whose
commit is an ancestor of the new one. Comparable means the same power state
(`ac` with `ac`, `battery` with `battery`), the same machine (CPU brand and
core counts), the same DEM and domain hashes, tolerance and extra mesh args.
Anything else, including `mixed` or `unknown` power on either side, gives
**NO BASELINE** with the mismatching field named: never a cross comparison.

With a baseline, per domain:

- **Time**: median `refine_s` at the default thread count and at each swept
  count. A median more than the threshold (default 5 %, `--threshold`) above
  the baseline's is a regression, reported with its size.
- **Quality**: a tolerance failure, any Delaunay violation, a lower worst
  angle or a higher max degree is a regression, whatever the time, unless
  `--accept-quality` waives the last two (see "Ruled by Ola").
- **Ceiling**: reported, not judged: speed-up at 20 threads over 1, and the
  best, against 2.2x and the baseline's.

Verdict: `ACCEPTED`, `REGRESSION: <domain> <measure> <baseline> -> <new> (+x %)`
for each, or `NO BASELINE: <reason>`. Exit code 0, 1, 2 respectively.

The 2026-09-26 logs are not converted (3 repeats, timings parsed from
stdout: a different method). The first acceptance run measures the previous
increment's merge commit with `--tree <worktree>` back to back with the new
one on the same power state; that pair is the first stored baseline.

## Seams for tests/python/test_bench.py

Loaded by `importlib` from `tools/bench.py`, as `test_session_state.py` does.
Nothing below needs the 1 m DEM or a real build.

- **`Runner` Protocol**: `run(argv, cwd) -> Completed(returncode, stdout,
  stderr, wall_s)`. Every subprocess goes through it: `pmset`, `sysctl`,
  `sw_vers`, `git`, `cmake`, the child. A fake runner returns canned output
  keyed on argv[0], so the orchestration (repeat count, interleaving, build
  refusal on a Debug cache, `--no-build` recorded, power read before and after)
  is tested with no process started. The Typer app gets it from a module-level
  factory the test overrides.
- The signatures in this list are the design's first sketch; where they
  differ, "Pinned by the red suite" and the code win.
- **`parse_pmset(text) -> Power`**, pure. Canned texts: `'AC Power'` with and
  without a battery line (desktop), `'Battery Power'` with a percent,
  `'UPS Power'`, empty and garbage (`unknown`).
- **`parse_child(stderr) -> Sample`**, pure; a missing or malformed JSON line
  is an error naming the run, not a zero. `median_stats`, `ceiling`: pure.
- **`quality(points, triangles, constraint_edges, tolerance, max_error)`**,
  pure over arrays, and **`read_vtk_ascii(path)`**. Hand-built meshes: a quad
  split the non-Delaunay way (1 violation), the same with that diagonal
  constrained (0), a cocircular square (ambiguous, exact test decides 0), a
  sliver (worst angle), a fan (max degree). One ASCII file written by
  `tin_engine.io.vtk_legacy.write_vtk` checks the reader against the writer.
- **`comparable(a, b)`, `find_baseline(root, record)`, `verdict(new, base)`**,
  pure over `RunRecord`s built in the test: AC vs battery is NO BASELINE,
  `mixed` is NO BASELINE, 4 % slower is ACCEPTED, 6 % is REGRESSION with the
  size, equal time with one Delaunay violation is REGRESSION.
- **`write_evidence(record, dir)`** into `tmp_path`: the files exist, the
  JSON round-trips, a rerun keeps the hand-written prose below the marker.
- **One real child run**, skipped when `_core` does not import: a tiny
  synthetic DEM from `tests/python/geotiff_fixtures.py`, threads 1 and 2,
  one repeat, `--no-build` against the installed package. It proves the
  `cli.refine` wrap still matches the real signature, which the fakes cannot.

## Dropped from the 2026-09-26 scripts

- Hard-coded paths (arguments now); `perl`/`bc` timing (`perf_counter`);
  `summarize.py`'s stdout regex (the child's JSON line); `build.sh`'s commit
  list and the `all_*.sh` loops (one `--tree` per call).
- Copying `_core*.so` into the worktree's `src_python/`: the `build-bench/pkg`
  symlink tree, since in the main tree a stray `.so` there is not gitignored.
- The `--stats` probe in `bench.sh` and the `refine_start_stride` wrap:
  `--stats` stays available through the extra mesh args.
- `quality.py`'s degree buckets (>= 12, >= 20) and count under 0.1 degree:
  the acceptance measures are worst angle and max degree. Its Delaunay check
  (float determinant, exact `Fraction` fallback inside the error bound) is kept.

## Size

Estimate, counted as CLAUDE.md section 2 counts: models 60, subprocess and
metadata (runner, pmset, sysctl, git, CMake cache) 70, build and pkg tree 35,
child 40, sampling loop 40, quality and VTK reader 75, comparison and
baseline search 60, evidence writing 70, Typer app 50: about 500 lines, under
the 700 ceiling. Plus `tools/bench.py` added to mypy's `files` in
`pyproject.toml`, so the strict gate covers it.

Measured at green: about 690 lines counted as section 2 counts them (blank
lines in, comments and docstrings out; 688 to 691 depending on whether the two
`# fmt: off/on` lines count), 582 without blank lines, so the estimate was low
by about 80 lines. The margin rests on 16 regions packed by hand under
`# fmt: skip` / `# fmt: off`; formatted normally it is about 780 counted lines
(827 raw). The next change to `tools/bench.py` has about 10 lines of room.

## Ruled by Ola (2026-09-27)

- **Time threshold: 5 % on the median, adjustable per run** (`--threshold`).
  In `REPORT.md` the min-to-max spread of 3 runs, over the median, was 0.3 to
  5.9 %, so single runs are noisier than 5 %. The first real acceptance run
  measures the spread of the median of 5 and reports it next to the threshold.
- **Quality loss is a regression by default.** `--accept-quality` overrides
  a lower worst angle or a higher max degree, never a tolerance failure or a
  Delaunay violation, and every use is written into that run's `README.md`
  and `run.json`, with the verdict still naming the measure it waived.

## Pinned by the red suite (@tester, 2026-09-27)

The seams above name functions but not every signature, field or string the
tests have to assert on. `tests/python/test_bench.py` pins the following; where
this section settles something the design above left open, it says so.

**Module surface.** `tools/bench.py` defines, at module level: `app` (Typer),
`make_runner()` (the factory the tests replace), `Completed`, `Runner`,
`Power`, `Child`, `ChildError(ValueError)`, `Sample`, `TimeStats`, `Ceiling`,
`MeshQuality`, `VtkMesh`, `RunRecord`, `Verdict`, `MARKER`, and the functions
`parse_pmset`, `combine_power`, `parse_child`, `median_stats`, `ceiling`,
`quality`, `read_vtk_ascii`, `comparable`, `find_baseline`, `verdict`,
`write_evidence`. The module must import with nothing but its path (the tests
load it by `importlib`, registered in `sys.modules` as `bench`).

**Power.** `Power(state, percent, raw)`, `state` one of `ac`, `battery`,
`mixed`, `unknown`; `percent` is `int | None`. `parse_pmset`: `'AC Power'` is
`ac` (percent from a battery line if one is present, else `None`),
`'Battery Power'` is `battery`. *Settled here:* `'UPS Power'` is `unknown` — it
is neither of the two states the rule compares. `combine_power(before, after)`
is `before`'s state when the states agree, else `mixed`.

**Child.** Invoked as `python tools/bench.py _child [--pkg DIR] --threads N --
<rasputin argv>`, where the rasputin argv starts with the subcommand (`mesh
--dem ... --out ... --binary`). `--threads 0` means "do not force": the CLI's
own default (the design's "default thread count"); any other N is forced into
the `refine` call. *Settled here:* with `--no-build` the parent passes no
`--pkg` and the child imports the installed `tin_engine`. The child's result
is one stderr line `BENCH ` followed by a JSON object with the keys `refine_s`,
`app_s`, `max_error`, `rounds`, `inserted`, `flips`. `parse_child(stderr, run)
-> Child` raises `ChildError` whose message contains `run` when there is no
such line, more than one, malformed JSON, or a missing key.

**Samples and statistics.** `Sample(domain, threads, repeat, refine_s, app_s,
proc_s, max_error, rounds, inserted, flips)`. `median_stats(samples) ->
list[TimeStats]`, one `TimeStats(domain, threads, n, median, min, max)` per
(domain, threads) over `refine_s`, sorted by that key. `ceiling(medians:
Mapping[int, float]) -> Ceiling(top_threads, at_top, best, best_threads)`:
speed-up over 1 thread at the largest thread count, and the best speed-up and
where it is reached; key 0 (the default) is ignored.

**Domains.** A domain is named `tile` (the whole DEM, no polygon) or by its
file's stem (`quarter`). *Settled here:* `--domain` is repeatable and the
value `tile` names the whole tile; default is `tile` and
`docs/benchmarks/2026-09-26/quarter.geojson`.

**Quality.** `quality(points, triangles, constraint_edges, tolerance,
max_error) -> MeshQuality(worst_angle, angle_median, share_under_1,
max_degree, within_tolerance, delaunay_checked, delaunay_ambiguous,
delaunay_violations, mesh_sha256="")`. Angles in degrees, the share a fraction,
the angle and degree fields as `tin_engine.stats.quality` computes them;
`within_tolerance` is `max_error <= tolerance`; the Delaunay counts are per
interior non-constraint edge, as `quality.py` counted them. `read_vtk_ascii(path)
-> VtkMesh(points (N, 3), triangles (T, 3), edges (E, 2), sha256)`, `sha256`
over the file's bytes from `POINTS` on; a `BINARY` file is a `ValueError`.

**Record.** `RunRecord` fields: `label`, `started` (aware datetime), `tree`
(`commit`, `dirty`), `bench_blob`, `child_argv`, `build` (`no_build`, `type`,
`cxx_flags_release`, `compiler`, `so_sha256`, all but `no_build` optional),
`machine` (`cpu_brand`, `p_cores`, `e_cores`, `memory_bytes`, `macos`,
`python`, `numpy`), `power` (a combined `Power`), `inputs` (`dem`,
`dem_sha256`, `domains`: list of `name`, `path`, `sha256`; `tolerance`,
`extra_args`), `samples`, `stats`, `quality` (domain name to `MeshQuality`),
`accept_quality`, `threshold_pct` (default 5.0), `verdict` (the verdict lines).
The two rulings are carried by the record itself, so `run.json` records them.

**Comparison.** `comparable(a, b) -> str | None`: `None`, or a reason naming
the first mismatching field by its record name (`power`, `cpu_brand`,
`p_cores`, `e_cores`, `dem_sha256`, `domains`, `tolerance`, `extra_args`).
`find_baseline(root, record, is_ancestor) -> tuple[Path, RunRecord] | None`:
the directory and record of the newest (by `started`) `root/*/*/run.json`
strictly older than `record`, comparable, and whose commit
`is_ancestor(candidate_commit, record_commit)` accepts. `verdict(new, base) ->
Verdict(status, lines, exit_code)` reads `threshold_pct` and `accept_quality`
from `new`; `base` may be `None`. Statuses `ACCEPTED`, `REGRESSION`,
`NO BASELINE` with exit codes 0, 1, 2. Line formats:
`NO BASELINE: <reason>`; `REGRESSION: <domain> refine_s[t=<N>] <base> -> <new>
(+x.x %)` for time (strictly more than the threshold); `REGRESSION: <domain>
<measure> <base> -> <new>` for `tolerance`, `delaunay_violations`,
`worst_angle`, `max_degree`; and, for a measure `--accept-quality` waived,
the same line prefixed `WAIVED (--accept-quality): ` instead of `REGRESSION: `.

**CLI.** *Settled here:* `run` takes `--out-root DIR` (default
`docs/benchmarks`), writes `<out-root>/<local date>/<label>/`, and searches
`--out-root` for a baseline; `compare` takes the same option. `--threshold` is
in percent (`--threshold 7.5`). `--threshold` and `--accept-quality` are
`run` options only: they are stored in the record, and `compare` re-judges
with what the new record carries. The verdict's exit code is the command's; a refused build (a
non-Release `CMakeCache.txt` under `<tree>/build-bench`) or a failing child
exits 3 and writes no evidence. Timing runs go in the order domain, then
repeat, then thread count (0 first, then the `--threads` list); the quality
run is one `--ascii` child at `--threads 0` per domain, after its timing runs.
`pmset -g batt` is read exactly twice, before the first child and after the
last. Machine facts come from `sysctl -n machdep.cpu.brand_string`,
`hw.perflevel0.physicalcpu`, `hw.perflevel1.physicalcpu`, `hw.memsize`, and
`sw_vers -productVersion`.

**Evidence.** `write_evidence(record, directory)` writes `run.json`
(`record.model_dump_json()`, round-tripping), `raw.tsv` (no header, one line
per sample: domain, threads, repeat, refine_s, power state), and `README.md`,
whose generated part ends at the line `MARKER`; whatever follows that line in
an existing README is kept on a rerun. The README names `2.2x` beside the
ceiling, every verdict line, and `--accept-quality` when it was used.

## Settled at green (@perf, 2026-09-27)

**Where the Release `.so` lands and how `pkg/` is assembled** (the red suite
pins only the Debug-cache refusal):

- Configure: `cmake -S <tree> -B <tree>/build-bench -DCMAKE_BUILD_TYPE=Release
  -DRASPUTIN_BUILD_PYTHON=ON -DRASPUTIN_BUILD_TESTS=OFF
  -DPython_EXECUTABLE=<the running interpreter> -DPYBIND11_FINDPYTHON=ON`;
  pybind11 is found by `find_package` as `build-pyext` finds it, with no
  `pybind11_DIR` passed. Build: `cmake --build <tree>/build-bench -j --target
  _core`. A failing configure or build exits 3 with no evidence.
- The `.so` lands where `pybind11_add_module` puts it, the build directory's
  root: exactly one `<tree>/build-bench/_core*.so`, or exit 3.
- `pkg/` is rebuilt on every run: `<tree>/build-bench/pkg/tin_engine/` holds a
  symlink to each entry of `<tree>/src_python/tin_engine/` except `__pycache__`
  and any `_core.*`, plus a **copy** of the `.so`, whose sha256 goes into
  `run.json`. The child gets `--pkg <tree>/build-bench/pkg` and refuses to run
  if `tin_engine` or `_core` resolves outside it.
- `build.compiler` is `CMAKE_CXX_COMPILER_ID` and `_VERSION` from
  `build-bench/CMakeFiles/*/CMakeCXXCompiler.cmake`.
- `quality()` angle and degree figures come from the **main tree's**
  `tin_engine.stats`, imported in the parent, whatever `--tree` measures: the
  same measuring code for both sides of a comparison.

