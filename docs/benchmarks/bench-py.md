# tools/bench.py: design

Status: design (@perf, 2026-09-27), not yet reviewed. No code or tests yet. The rule it serves is
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
                             [-- extra mesh args, e.g. --start-min-angle 0]
python tools/bench.py compare NEW_DIR [--baseline DIR]
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
  angle or a higher max degree is a regression, whatever the time.
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

## Open for Ola

- The time threshold. 5 % on the median is proposed. In `REPORT.md` the
  min-to-max spread of 3 runs, over the median, was 0.3 to 5.9 % (worst: the
  quarter circle at head, 0.237 to 0.251 s), so single runs are noisier than
  5 %; the median of 5 is expected to be steadier, not yet measured to be.
- Whether a quality change is always a regression. As designed, a lower worst
  angle or a higher max degree fails the run even when an increment trades
  them on purpose; the alternative is an `--accept-quality` flag whose use the
  README records.
