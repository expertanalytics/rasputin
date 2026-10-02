# One resample block's memory on the basin grid (@architect, 2026-10-02)

Input to `docs/research/basin-memory-options.md` §1. Branch
`worktree-basin-memory` (from `worktree-15c-2` at `c47bff6`); `_core` from
`.claude/worktrees/15c-2/build-bench/` (the probe never calls into it, only
imports it); Python 3.14.7, NumPy from the main `.venv`; Apple M1 Max, 32 GiB.

| input | value |
|---|---|
| domain | `../rasputin_data/sao_francisco_piece/bho2017_level2_76_raw.geojson` |
| target CRS | `docs/benchmarks/2026-10-02/basin-anadem/runs/basin_out_crs.wkt` |
| grid | `target_grid_for(basin, CRS, 30)`: 50,315 rows x 40,943 cols |
| source window | fake EPSG:4674, 0.000269494585236 deg, basin box + 0.5 deg; zero-byte views |
| block | 256 rows (the default `block_rows`) x 40,943 cols, middle row, threads 1 |
| measured by | `tracemalloc` peak (NumPy arrays; pyproj's C buffers not traced) |

Result (`blockprobe.out`): **1.562 GB peak, 149.0 B per node**, 1.9 s.
`resample` runs `os.cpu_count()` (10 here) blocks at once, so up to ~15.6 GB
(arithmetic, not measured).

```bash
mkdir -p /tmp/probe_pkg && cp -r src_python/tin_engine /tmp/probe_pkg/ \
  && cp <a built _core*.so> /tmp/probe_pkg/tin_engine/
python docs/benchmarks/2026-10-02/basin-memory-probe/blockprobe.py /tmp/probe_pkg \
  docs/benchmarks/2026-10-02/basin-anadem/runs/basin_out_crs.wkt \
  ../rasputin_data/sao_francisco_piece/bho2017_level2_76_raw.geojson
```

The 1.562 GB peak includes the block's own output canvas and its copy by
`DemTile(...)` (~8 B per node), so 149 B per node is slightly high for the
temporaries alone.
