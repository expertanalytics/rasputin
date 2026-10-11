# Increment 26a quick check: MapBiomas in the catalogue, and class windows (@perf, 2026-10-11)

Branch `worktree-mapbiomas-26a`, head `f51f0a18`, base master `de9b6c03`
(the merge of increment 35). The rule: `docs/increments/perf-quick-check.md`,
section 3, the quick check only. 26a touches no refinement or mesh code, so
no `tools/bench.py` 1 m run and no thread sweep. It does touch the DEM
decode path: `8eec354b` moved the block-decoding loop of `decode_window`
(`src_python/tin_engine/io/cog.py`) into a shared `_decode_into`, which the
new `decode_class_window` also calls; `decode` is the largest phase of the
`geilo-al-ramp` case. So the check is: every case's mesh byte-identical to
the baseline's, every total within its band, and `decode` within its band
where the tool judges it (a phase at 5 % or more of the case's total).

## Verdict: ACCEPTED

`tools/bench_quick.py run` at `f51f0a18`: **`NO CHANGE`, exit 0**
(`run.log`). The tool prints a `MESH CHANGED:` line for any case whose mesh
SHA-256 differs from the baseline's, a `BROKEN:` line for a `max_error` over
the tolerance, and a `SLOWER:`/`FASTER:` line for any measure outside its
band (the case's own spread, at least 5 %, and more than 0.05 s), on the
total and on every phase at 5 % or more of the total; `decode` is judged on
every case. It printed none of them, so:

| gate | result |
|---|---|
| mesh SHA-256 unchanged, all 7 cases | met (no `MESH CHANGED` line) |
| `max_error` within tolerance, all 7 cases | met (no `BROKEN` line) |
| every case's total within its band | met (no `SLOWER`/`FASTER` line) |
| `decode` within its band, `geilo-al-ramp` and the tile cases | met (no `SLOWER`/`FASTER` line) |

Baseline: `docs/benchmarks/quick/baseline-ac.json`, commit `357808b1`, AC
power, Apple M1 Max (8 P + 2 E cores, 32 GB), macOS 27.0, Python 3.14.7,
numpy 2.5.4, shapely 2.2.0 / GEOS 3.14.1, pyproj 3.8.0 / PROJ 9.8.1; saved
by increment 35's acceptance (#229). Its per-case medians are the "head"
column of `docs/benchmarks/2026-10-11/cleanup-speed/README.md`: tile 0.416 s,
tile t=1 0.699 s, quarter 0.355 s, numedalslagen 8.829 s, lagan 14.542 s,
geilo-al-ramp 0.329 s (decode 0.138 s), romsdal-slope 3.578 s. The tool
compared this run's medians with those and found each inside its band; the
run's own medians are not on disk, see "Limits" below.

## Hotspots, as printed (`docs/benchmarks/quick/hotspots.toml` rules only Lagan `features clip`)

```
HOTSPOT: tile refine 44 % of the run (not in hotspots.toml)
HOTSPOT: tile t=1 refine 67 % of the run (not in hotspots.toml)
HOTSPOT: tile t=1 refine: scan (parallel) 50 % of the run (not in hotspots.toml)
HOTSPOT: quarter refine 47 % of the run (not in hotspots.toml)
HOTSPOT: numedalslagen features clip 68 % of the run (not in hotspots.toml)
HOTSPOT: numedalslagen features clip: clean-up 67 % of the run (not in hotspots.toml)
HOTSPOT: geilo-al-ramp decode 42 % of the run (not in hotspots.toml)
HOTSPOT: romsdal-slope refine 56 % of the run (not in hotspots.toml)
```

The same eight lines as the baseline run printed
(`docs/benchmarks/2026-10-11/cleanup-speed/quick/head-run.log`), six at the
same share; tile `refine` 46 -> 44 % and quarter `refine` 48 -> 47 % (both
within those cases' run-to-run noise on a 0.4 s case). geilo-al-ramp `decode`
is 42 % on both, 0.138 of 0.329 s on the baseline. None is new to 26a. Lagan
`features clip` (ruled at 40.7 %) did not print, so it is under 50.7 % of
Lagan's run, as after increment 35.

## Method

* Tree `/Users/skavhaug/projects/rasputin/.claude/worktrees/mapbiomas-26a`
  at `f51f0a18`, clean. `tools/bench_quick.py run` builds Release into
  `build-bench/` and runs the seven cases of `docs/benchmarks/quick/cases.toml`
  from there (3 to 5 timed runs per case, the median judged), so no stale
  extension is measured. `RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data`.
* Power: AC, battery 100 %, `caffeinate` running, before and after
  (`power-state.txt`). Nothing else ran beside the timing run
  (`pgrep` for cmake, clang, pytest, bench and `rasputin mesh` found nothing).
* Full tool output: `run.log` (Python's two `sys.prefix` RuntimeWarnings from
  the Homebrew interpreter are noise, not from the tool).
* Meshes: written by the tool to its own temporary directories and removed;
  regenerate by re-running the command above.

## Limits

* `bench_quick.py` writes a run's record (per-case medians, spreads and mesh
  hashes) only with `--save-baseline`, which replaces
  `docs/benchmarks/quick/baseline-ac.json`. A second run with that flag,
  meant to be copied here and the baseline restored, was refused by the
  session's permission classifier, so this directory holds the tool's
  verdict and hotspot lines, not the medians. A `--record <path>` option that
  saves the record without touching the baseline would close this gap; it
  is tooling (`tools/`) and goes through the TDD loop.
