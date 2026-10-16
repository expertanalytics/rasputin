# 30d, the speed check of the built fix (`4f13ee90`)

`@perf`, 2026-10-07, as `docs/increments/30d-outline-buffer-speed.md` section 7
asks: one case, one warm-up and three runs, against the master medians in
`README.md` (same machine, same inputs, master not re-run).

## Result

Lagan with the European GeoPackage, `--tolerance 10 --binary`, AC power before
and after every run. Median in bold, then the three runs, in seconds.

| measure | master `f81b20b7` | built `4f13ee90` | change |
|---|---|---|---|
| total | 103.38 | **16.97** (16.83, 16.97, 17.07) | -83.6 % (6.1 times faster) |
| decode | 35.12 | **4.65** (4.65, 4.64, 4.68) | -86.8 % |
| features read | 57.23 | **1.07** (1.07, 1.07, 1.06) | -98.1 % |
| features clip | 6.79 | **6.89** (7.03, 6.85, 6.89) | +1.5 %, inside the 5 % band (unchanged code) |

- **Mesh:** the `.vtk` SHA-256 equals master's on all four runs
  (`7defb85b…`), and the quality section of `--stats` is the same text
  (867,678 triangles, largest difference 9.9998 m against a 10 m tolerance).
- **Numedalslågen** (UTM33 GeoPackage, one run, 7.09 s): SHA-256 `73b0d8fc…`,
  equal to master's. The gate keeps this outline on GEOS's single buffer.

**Largest phase now:** `features clip`, 6.89 s, **40.7 %** of the run. It was
not changed by 30d and was 6.6 % of master's run. It is at or above the 40 %
hotspot line, so it goes to Ola. Why it takes 6.9 s was not measured
(`README.md` notes 15.9 M CORINE vertices handed over from the GeoPackage
against 1.74 M from the GeoJSON window; a likely cause, not a measured one).
Next: `decode` 27.5 %.

## Method

- Code: `git archive 4f13ee90 src_python/tin_engine` in the session
  scratchpad, with the Release `_core` from the main checkout's `.venv`
  (SHA-256 `86765ff1…`, the same file the master baseline ran;
  `include/`, `src/`, `bindings/` and `CMakeLists.txt` are unchanged between
  `f81b20b7` and `4f13ee90`). Every run's stderr names the copy's
  `tin_engine/__init__.py` and `_core` (`raw/runs/built_*_stderr.txt`).
- Script: `scripts/timed_built.sh` through `scripts/run_mesh.sh` and
  `scripts/launch.py` (no trace, no patch). Machine and versions as in
  `README.md`. Nothing else ran.
- Raw: `raw/runs/built_*` (stats, stderr, power) and the hashes appended to
  `raw/runs/vtk_sha256.txt`. Meshes were written to the scratchpad and
  removed; the script regenerates them.
- Wall clock for this check: about 6 minutes of runs (23:16 to 23:22), 7 minutes start to commit.
