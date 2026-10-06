# Audit PR A (`audit-lattice`, DEM grid code) acceptance (@perf, 2026-10-06)

Branch `worktree-audit-lattice`, HEAD `1a84ea9` (code `b3b38d2`), against
`618328b` (PR B's approved head, worktree `audit-crs`). The rule is
`docs/increments/python-audit-pr-a.md`, section `@perf`: every mesh file
byte-identical, on the 1 m set and on one reprojected Velhas run.

## Verdict

**ACCEPTED.**

- **Meshes byte-identical, every run, both gates.** The hash `bench.py`
  records (`mesh_sha256`, the `--ascii` mesh from `POINTS` on), and the
  whole `.vtk` file as well (`shasum -a 256`):

  | set | domain | `mesh_sha256` (base = head, every run) | whole file |
  |---|---|---|---|
  | 1 m | tile | `11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c` | `298905f3…`, 9 of 9 equal |
  | 1 m | quarter | `a60fb597b5b738268a3bd27730989392ee86fb030fe060b78fd3a0bcdbb60a34` | `9c7adff5…`, 9 of 9 equal |
  | Velhas, EPSG:31983 | outline | `3d3ea54fcf16b1fca538d2bfce312b8e1951b7a1c3af23535e0a79ddf48503cd` | `aa4218c8…`, 2 of 2 equal |

  (9 = four base runs, four head runs and the head build run `pra-head-r0`.)
- **Quality identical** (it follows from the hashes). 1 m tile / quarter:
  worst angle 0.6296° / 0.3955°, max degree 74 / 18, `max_error <= 1`,
  0 Delaunay violations. Velhas: worst angle 0.0525°, max degree 14,
  `max_error <= 10`, 0 Delaunay violations of 1,059,205 checked.
- **Timing within noise.** 1 m set, four runs a side pooled: refine
  -2.5 .. +2.6 % over every cell (threads default and 1 to 20, both domains);
  the four base runs alone spread by a median 2.1 % and at most 5.3 % per
  cell. At 1 / 20 threads refine is +0.7 / -0.3 % (tile) and -0.1 / +2.6 %
  (quarter); process time is within -0.5 .. +0.6 % on every cell. Velhas,
  one run a side: refine -2.0 .. +1.3 %, process time -1.2 .. +0.7 %.
  Tables: `tables-1m.md`, `tables-velhas.md`.

## Method

- **Machine and power.** Apple M1 Max. AC for every run (`pmset -g batt`
  before and after each, in each `run.json`); 80 %, not charging, until
  the Velhas head run, which started at 81 % charging. caffeinate running.
- **Builds.** `bench.py run` built Release `_core` into each tree's
  `build-bench`; both hash `6215e9f2…` (no C++ change between the trees).
  Hardening: the default (on), `libc++ fast` in every run.
- **Script.** `tools/bench.py` blob `2b7dca1a…`, unchanged between
  `618328b` and `1a84ea9`, run from this worktree's `.venv` (Python 3.14.7,
  audit-crs's pinned dependencies) for both trees. The head runs record the
  tree dirty: the only change was this untracked evidence directory.
- **1 m set** (`pairs.sh`, log `pairs.log`): DEM
  `tests/fixtures/dem_archive/7908_3_10m_z33.tif`, domains the tile and
  `docs/benchmarks/2026-09-26/quarter.geojson`, tolerance 1, 5 repeats,
  threads default and 1 to 20. Order base, head, head, base; head, base,
  base, head (`pra-base-r1..4`, `pra-head-r1..4`), 02:22-02:42Z.
  `pra-head-r0` is the run that built the head tree, before the batches; it
  is not in the timing pool (its parent process read `tin_engine.stats`
  from the base package, a file with no diff between the trees).
- **Velhas** (`velhas.sh`, log `velhas.log`): the increment's command with
  absolute paths, `--dem
  /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_anadem_window_epsg4674.tif
  --domain .../bho2017_5k_76949_outline_epsg4674.geojson --tolerance 10 --
  --out-crs EPSG:31983`, defaults otherwise. Checked first at `618328b`
  (threads 0 and 20, 1 repeat: it ran, 15 s). Then base, head
  (`vel-base-r1`, `vel-head-r1`), 02:44-02:57Z. The meshes are in
  EPSG:31983 (`crs` field, `POINTS` at x 534,630, y 7,933,196).
- **Tables.** `summarize.py <dir> pra|vel` pools the runs' samples and
  prints the hashes, the per-cell medians and the base-vs-base spread.
- **Meshes** were in the session scratchpad and are not kept; rerun
  `pairs.sh` / `velhas.sh` (paths inside) to regenerate them.
