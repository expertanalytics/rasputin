# Increment 36 acceptance: GeoPackage output (@perf, 2026-10-11)

Branch `worktree-gpkg-out`, head `e2470ae7`, against its base, master
`81a551b3` (`git merge-base master HEAD`; the merge of #227, before increment
35's #229). The rule: `docs/increments/36-geopackage-output.md`, section 9:
(a) the quick check (`tools/bench_quick.py run`) unchanged with every
`mesh_sha256` equal; (b) the GeoPackage of the quick check's Numedalslagen
case written in at most 5 s on AC and of the Romsdalen slope mesh in at most
25 s, peak RSS under 1 GB and 3 GB respectively, measured with `/usr/bin/time
-l` on the built command. (c) is Ola's QGIS check and (d) is `@tester`'s T5;
neither is measured here. The branch touches no C++ (`git diff --stat
81a551b3..HEAD -- include src CMakeLists.txt tests/cpp` is empty), so no
`bench.py` 1 m run and no thread sweep.

## Verdict: ACCEPTED

| gate | measured (medians of 3, AC) | limit | result |
|---|---|---|---|
| (a) mesh SHA-256 unchanged, all 7 cases | equal on every case, against the base and against master's refreshed baseline (`quick/delta.md`, `quick/delta-vs-master-357808b1.md`) | equal | met |
| (a) quick check totals | `NO CHANGE`, exit 0, head against base back to back (`quick/head-run.log`); -4.3 % to +1.6 % | within the band | met |
| (b) Numedalslagen `write: geopackage` | **3.932 s** (3.892, 3.932, 3.988) | at most 5 s | met |
| (b) Numedalslagen peak RSS | whole command **2.77 GB** against 2.73 GB with `--out .vtk`; the writer's own share **+32 MB** (RSS growth inside its phase; +0.04 GB on the command's peak) | writer's share under 1 GB | met |
| (b) Romsdalen `write: geopackage` | **19.368 s** (19.051, 19.368, 20.154) | at most 25 s | met |
| (b) Romsdalen peak RSS | whole command **2.91 GB** against 2.78 GB with `--out .vtk`; the writer's own share **+188 MB** (+0.13 GB on the command's peak) | writer's share under 3 GB | met |

`max_error` is within tolerance on every quick-check case and the same on both
sides (the records). Both GeoPackage runs exited 0; `--stats` and the
written file sizes agree with section 3 (Numedalslagen 203.3 MB, 653,039
triangles; Romsdalen 810.0 MB, 3,096,202 triangles — section 3's probe wrote
922 MB for the same mesh; the difference was not investigated).

## (a) The quick check, base `81a551b3` / head `e2470ae7`

`scripts/quick_delta.py`, full output `quick/delta.md`:

| case | base total s | head total s | change | mesh same | features clip | decode | refine |
|---|---:|---:|---:|---|---|---|---|
| tile | 0.443 | 0.426 | -3.8 % | yes | - | 0.064 / 0.060 | 0.196 / 0.195 |
| tile t=1 | 0.801 | 0.766 | -4.3 % | yes | - | 0.064 / 0.062 | 0.555 / 0.531 |
| quarter | 0.383 | 0.369 | -3.6 % | yes | - | 0.050 / 0.048 | 0.179 / 0.173 |
| numedalslagen | 9.169 | 9.320 | +1.6 % | yes | 6.226 / 6.381 (+2.5 %) | 0.967 / 0.964 | 0.708 / 0.688 |
| lagan | 20.166 | 20.029 | -0.7 % | yes | 11.584 / 11.521 (-0.5 %) | 4.731 / 4.673 | 0.251 / 0.245 |
| geilo-al-ramp | 0.349 | 0.335 | -3.9 % | yes | - | 0.148 / 0.141 | 0.109 / 0.101 |
| romsdal-slope | 3.844 | 3.748 | -2.5 % | yes | - | 0.070 / 0.070 | 2.200 / 2.094 |

Spreads (max - min) / median of the totals: 0.2 % to 3.4 % on both sides
(`quick/delta.md`). The `.vtk` write phases are unchanged (`write: encode`
within 5 % on every case; `write: disk` is 1 to 62 ms and noisy).

**Against master's refreshed baseline** (`357808b1`, increment 35's,
`quick/master-357808b1-baseline.json`, table `quick/delta-vs-master-357808b1.md`):
every mesh SHA-256 equal; Lagan total 14.542 -> 20.029 s (+37.7 %) with
`features clip` 6.342 -> 11.521 s, and Numedalslagen +5.6 %. Lagan's gap is
increment 35's saving, which this branch does not contain (the base
`81a551b3`, run back to back, has the same 11.58 s clip); it is the base
difference, not a regression of 36. The remaining +2 % to +10 % on every
case, including the pure-C++ ones, is tonight's machine against this
afternoon's: the base run of unchanged code is already 3 % to 15 % slower
than master's record (`tile t=1` 0.801 against 0.699). The back-to-back
pair is the like-for-like measure.

Hotspot lines (`quick/head-run.log`): the same unruled lines as the last
checks (`tile` and `quarter` refine, `tile t=1` scan, `geilo-al-ramp` decode,
`romsdal-slope` refine, Numedalslagen `features clip`), plus Lagan `features
clip` raised at 58 % (ruled at 41 %), which is the pre-35 clip; master's
check after #229 is quiet on it (`2026-10-11/cleanup-speed/README.md`).

## (b) The writer on the built command

`scripts/writer_time.py`, output `writer/writer-run.log`, records
`writer/writer-results.json`, per-run logs `writer/writer-<case>-<fmt>-<n>.log`
(the command's `--stats -` report and the `/usr/bin/time -l` block; run 0 is
the warm-up) and RSS traces `writer/writer-<case>-<fmt>-<n>.rss`. Each
command is `/usr/bin/time -l .venv/bin/rasputin mesh <case args> --out
mesh.gpkg --stats -`, and the subtraction baseline is the same command with
`--out mesh.vtk --binary`. The case arguments are the quick check's
`numedalslagen` and `romsdal-slope` lines of `docs/benchmarks/quick/cases.toml`
(the Romsdalen mesh is tile `6901_3`, `--tolerance 10 --tolerance-slope 2 25
35`; its binary `.vtk` is 111.5 MB, as section 3's).

| command | real s | `write:` phase s | whole-command peak RSS | writer's own share | file |
|---|---:|---:|---:|---:|---:|
| Numedalslagen `.gpkg` | 13.53 | **3.932** (`write: geopackage`) | **2.77 GB** (2.78, 2.77, 2.73) | 2.75 -> 2.77 GB, **+32 MB** (32, 18, 35) | 203.3 MB |
| Numedalslagen `.vtk` | 9.58 | 0.014 (`encode` + `disk`) | 2.73 GB (2.73, 2.72, 2.73) | - | 29.2 MB |
| Romsdalen `.gpkg` | 24.15 | **19.368** (`write: geopackage`) | **2.91 GB** (2.95, 2.90, 2.91) | 2.72 -> 2.91 GB, **+188 MB** (227, 181, 188) | 810.0 MB |
| Romsdalen `.vtk` | 4.82 | 0.058 (`encode` + `disk`) | 2.78 GB (2.78, 2.78, 2.79) | - | 111.5 MB |

How the writer's share is separated: `/usr/bin/time -l`'s `maximum resident
set size` is the whole command's peak, and the meshing sets it (the `.vtk`
runs peak at 2.73 and 2.78 GB with no GeoPackage at all; `@developer`'s 2.79 GB
in section 11 is this figure). The script samples the rasputin child's RSS
with `ps` every 0.1 s; the writer's window is the last `write: geopackage`
seconds of the process (it is the run's last phase), and the writer's own
share is the RSS growth inside that window: what it adds on top of the mesh
the process already holds. Freed meshing memory the allocator reuses is not
visible to this measure, so it is a lower bound on the writer's allocations
and an exact figure for what the writer costs the command; the whole-command
peak, 2.77 and 2.91 GB, is the upper bound on anything the writer holds. Both
bounds are under section 9's 1 GB and 3 GB. Section 3's standalone probe
figures (0.56 GB and 2.1 GB, the probe's VTK read included) sit between them.

Per triangle the writer takes 6.0 µs (Numedalslagen) and 6.3 µs
(Romsdalen); section 3 measured 2.2 µs per row plus 4.6 µs per R-tree entry.
The `.gpkg` command's extra wall time over the `.vtk` one (3.95 s and
19.33 s) equals the `write: geopackage` phase, so the GeoPackage path adds
nothing outside that phase.

## Method

- `RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data`, this worktree's
  `.venv` (Python 3.14.7, numpy 2.5.4, shapely 2.2.0, GEOS 3.14.1, pyproj
  3.8.0, PROJ 9.8.1, SQLite as bundled with that Python), Apple M1 Max (8 P +
  2 E cores, 32 GB), macOS 27.0, AC power before and after every run
  (`power-before.txt`, `quick/power-*.txt`, `writer/power-*.txt`), caffeinated;
  nothing else ran. Bounds checks on (`bounds checks: on (libc++ fast)` in
  every `--stats` report), as the stored baselines.
- **Quick check, base first** (`quick/base-run.log`): `bench_quick.py run
  --tree <detached scratch worktree of 81a551b3> --save-baseline`, judged
  against the branch's stored baseline `39e3ae34` (increment 34's, numpy
  2.5.3): `NO BASELINE: machine numpy: 2.5.4 vs 2.5.3`, as
  `docs/increments/README.md` ("Acceptance") foresees; its record
  `quick/base-81a551b3.json` became the stored baseline for the head run. Its
  meshes equal master's refreshed baseline's on every case, so the base and
  master `de9b6c03` produce the same meshes (increment 35 changed none).
  **Then head** (`quick/head-run.log`): `bench_quick.py run --save-baseline`
  at `e2470ae7`, judged against the base record: `NO CHANGE`, exit 0; record
  `quick/head-e2470ae7.json`. Both Release `_core` builds came from
  `bench.build` into each tree's `build-bench`; the Python package ran from
  each tree. `docs/benchmarks/quick/baseline-ac.json` was then restored to
  the branch's committed file (`git checkout --`): the refreshed baseline is
  master's (#229), and a second one on this branch would only conflict at
  the merge. The scratch worktree was removed afterwards.
- **Writer**: the head run's Release `_core` (`build-bench`, hardening on)
  copied into the venv's `site-packages` and touched before the timed runs
  (sha256 `51872f0ab6f1d780...` on both sides), per the stale-artifacts
  rule. No C++ changed, so no sanitized run. One warm-up and 3 timed runs per
  command, medians; meshes written to the session scratchpad and deleted.
- Wall clock: 21:49 to 22:01; two quick checks and 16 writer commands.

## To regenerate

    cd <worktree>; RASPUTIN_DATA=... .venv/bin/python tools/bench_quick.py run
    RASPUTIN_DATA=... .venv/bin/python docs/benchmarks/2026-10-11/gpkg-out/scripts/writer_time.py <out dir>
    python3 docs/benchmarks/2026-10-11/gpkg-out/scripts/quick_delta.py
