# 15c-2 after the crash fix: the crash rate, and the whole São Francisco basin (@perf, 2026-10-02)

Branch `worktree-15c-2` at `8b7b0dd`. Apple M1 Max (8 P + 2 E cores, 32 GiB),
macOS 27.0, Python 3.14.7, **AC power** for every run (`pmset -g batt`
before and after, in `runs/crash_loop_pmset_*.txt` and `runs/basin.json`).
`_core` was built Release (`-O3 -DNDEBUG`, AppleClang 21) by `tools/bench.py`'s
`build()` into `build-bench/` from this tree (sha256 `41a22a6168db…`), after
`include/terrain/refinement/refine_points.hpp` last changed. The run used a
venv of its own with `build-bench/pkg` on its path, not the main `.venv`.

## Verdict

- **Crash rate: 0 of 200.** The random crash in the 15c-2 acceptance
  (`../15c-2-acceptance/`: 7 of 47 process starts, 15 %) did not occur
  again. If the rate were still 15 %, the chance of 0 crashes in 200 starts
  would be about 1e-14.
- **Velhas at 5 m:** the independent final check still finds **0 source nodes
  over tolerance** in both CRSs. Its control fails as it should.
- **The whole basin: NOT DONE. It does not fit in memory on this machine,
  even at 50 m.** The 50 m run reached a **33.2 GB peak memory footprint**
  (`time -l`) on the 32 GiB machine. Swap grew by 17 GB, and at 189 s, not
  finished, the run was stopped by the script's swap-growth guard. A rerun
  with the guard raised to 60 GB was refused by the session's permission
  classifier, and no attempt was made to get round that. No basin mesh exists,
  and 20 m, 10 m and 5 m were not run: their grids are the same size as
  50 m's, so they need at least as much memory. That decision is Ola's
  (see *Open*).

## (a) Crash rate on the Velhas piece

`crash_loop.py` ran this command 200 times, alternating between the two CRSs
(100 starts each):

```
rasputin mesh --dem anadem-v1 --cache ../rasputin_data/cache \
  --domain ../rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson \
  --out-crs <CRS> --tolerance 10 --binary --out <scratch>
```

| `--out-crs` | starts | non-zero exit or signal | median wall s |
|---|---:|---:|---:|
| EPSG:31983 | 100 | 0 | 3.78 |
| the suggestion pasted from the refusal (`runs/velhas_out_crs.wkt`) | 100 | 0 | 3.99 |

No new `python*.ips` crash report appeared in
`~/Library/Logs/DiagnosticReports/` during the loop (`runs/crash_loop.json`,
`new_crash_reports: []`). The suggestion is the first line of the refusal
(`runs/velhas_refusal.out`), cut out by `sed` exactly as a shell user would
paste it. It is WKT2 on SIRGAS 2000 (EPSG:4674's datum), transverse Mercator
at lon_0 -44.1, k 0.999972.

**Independent final check at 5 m** (`recheck_5m.py`, which reuses
`../15c-2-acceptance/run_geo.py`'s `independent_check`: matplotlib's linear
interpolator at every valid node of the ANADEM source window inside the
domain; `runs/recheck_5m.json`):

| `--out-crs` | phase 2 source nodes | phase 2 inserted (rounds) | interior over 5 m | strip over 5 m |
|---|---:|---:|---:|---:|
| EPSG:31983 | 18,791,066 | 142,643 (12) | 0 of 13,791,972 | 0 of 27,806 |
| pasted WKT | 18,860,772 | 143,257 (10) | 0 of 13,791,979 | 0 of 27,799 |

Control: the EPSG:31983 mesh moved 15 m east finds 1,048,434 interior nodes
over tolerance (max 47.6 m). The node counts and the insertions match the
15c-2 acceptance's 5 m rows.

## (b) The whole São Francisco basin

Domain: BHO level-2 ottobasin 76, `bho2017_level2_76_raw.geojson`
(EPSG:4674, one polygon, 77,310 vertices, sha256 `8cdfda07…22356a`),
633,624 km² in the suggested CRS. The domain reader accepted it as given, and
it was not simplified. Without `--out-crs` the CLI refuses and suggests
(`runs/basin_refusal.out`) transverse Mercator at lon_0 -42, k 0.997548 on
SIRGAS 2000, with a worst scale error of 0.245 % over the box.
`runs/basin_out_crs.wkt` holds the pasted value.

| attempt | what happened | wall s | peak footprint | max RSS | swap growth |
|---|---|---:|---:|---:|---:|
| 1 | refused: **400 of the 8,700 blocks needed are not in the cache** (the cache held 8,300) | 1.3 | n/a | n/a | n/a |
| 2 | stopped at 189 s by the swap-growth guard (16 GB); not finished | 190 | 33.2 GB | 16.5 GB | 12.7 → 29.7 GB |

Attempt 1: the refusal named the fix, `rasputin fetch anadem-v1 --domain …
--out-crs …`. Run as given, it fetched the 400 blocks: 104 MB in 200
requests, 35.5 s (`runs/fetch_missing_blocks.out`; the manifest before the
fetch is `runs/cache_manifest_before_fetch.json`). The 8,300 blocks already
cached had been fetched for the same outline, but without this `--out-crs`.
The mesh's source region is the image of the target grid's rectangle,
which needs about 5 % more blocks. This was measured once and not traced
further.

Attempt 2: the memory trace (5 s samples) is
`runs/basin/t50_attempt2_stopped.mem.csv`. RSS swung between 6 and 13 GB
while swap grew, so most of the footprint was compressed or swapped out.
**Which phase held the memory was not measured.** The run was stopped before
the `--stats` report was written. The following sizes come from the code's
own planning (`target_grid_for`) and the cache header, without a run:

| array | nodes | bytes |
|---|---:|---:|
| source mosaic (float32), the basin's lon/lat box at ANADEM's 0.000269° | 2.13 G | 8.5 GB |
| resampled target grid (float32), 50,315 × 40,943 at 30 m | 2.06 G | 8.2 GB |
| source nodes inside the domain (the final check keeps those inside the grid's rectangle, so more) | ~0.73 G | ≥ 14.6 GB at 16 B xy + 4 B z |

These three arrays alone come to more than 31 GB, which is thought to be why
the run does not fit. The memory need is set by these sizes, not by the
tolerance, so 20 m and 10 m would need at least as much as 50 m.

## Open

- **ASK OLA:** the whole basin needs more than 32 GiB as the pipeline stands
  (33.2 GB footprint measured at 189 s, not finished). There are three
  options. (1) Allow the run to page into swap with a larger guard; the
  classifier refused that in this session. (2) Mesh the basin in pieces, for
  example by BHO level-3 basin. (3) Make the geographic path stream the source
  mosaic and the check points instead of holding all of them at once. That
  would be a design item for @architect.
- The independent check of a basin mesh was not run, because there is no
  mesh. No basin `.vtk` was written to
  `../rasputin_data/sao_francisco_piece/meshes/`, and its README was not
  changed.

## Reproduce

From the worktree root, with `PY` and `RASPUTIN` the venv's python and its
`rasputin` (any venv with numpy, shapely, pyproj, tifffile, imagecodecs,
typer, pydantic, matplotlib, and `build-bench/pkg` on its path), and
`D=../rasputin_data`:

```bash
E=docs/benchmarks/2026-10-02/basin-anadem
O=$D/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson
B=$D/sao_francisco_piece/bho2017_level2_76_raw.geojson
$RASPUTIN mesh --dem anadem-v1 --cache $D/cache --domain $O --tolerance 10 --out x.vtk 2>&1 \
  | head -1 | sed -e "s/^--out-crs '//" -e "s/'$//" > $E/runs/velhas_out_crs.wkt
$PY $E/crash_loop.py $RASPUTIN $D/cache $O $E/runs/crash_loop.json SCRATCH 100 \
  EPSG:31983 "$(cat $E/runs/velhas_out_crs.wkt)"                       # ~13 min
echo EPSG:31983 > epsg.txt
$PY $E/recheck_5m.py $RASPUTIN $D/cache $O \
  $D/sao_francisco_piece/bho2017_5k_76949_anadem_window_epsg4674.tif SCRATCH \
  $E/runs/recheck_5m.json epsg.txt $E/runs/velhas_out_crs.wkt           # ~2.5 min
$RASPUTIN fetch anadem-v1 --domain $B --out-crs "$(cat $E/runs/basin_out_crs.wkt)" --cache $D/cache
$PY $E/run_basin.py $RASPUTIN $E/runs/basin.json $E/runs/basin t50 -- mesh --dem anadem-v1 \
  --cache $D/cache --domain $B --out-crs "$(cat $E/runs/basin_out_crs.wkt)" --tolerance 50 \
  --ascii --out OUT.vtk --stats $E/runs/basin/t50.stats.md
```

Data credits as in `../basin-piece-anadem/README.md` (ANADEM, CC BY 4.0;
BHO 2017, ANA).
