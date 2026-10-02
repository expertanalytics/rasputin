# The São Francisco basin at 20 m, one mesh per BHO level-3 sub-basin (@perf, 2026-10-02)

An interim measurement of option (c) in `docs/research/basin-memory-options.md`
(branch `worktree-basin-memory`): mesh each level-3 unit separately, so that
no single run holds the whole basin. Per Ola's B13, the BHO outlines are used
as measurement domains. **The unit edges do not match.** Each unit is meshed
on its own, so neighbouring meshes have different vertices along their shared
border. They share one CRS and line up in space, but together they are not one
conforming mesh.

Branch `worktree-basin-level3`, cut from `worktree-basin-phases` at `55b701f`.
No production code changed.

## Method

- **Machine**: Apple M1 Max (8 P + 2 E cores, 32 GiB), macOS 27.0, Python
  3.14.7. **AC power** for every run (`pmset -g batt` before and after each
  run, in `runs/<unit>.json`). Swap stood at about 17.6 GB in use before the
  runs, left over from earlier work. No run added to it: `swap_max_mb` equals
  `swap_before_mb` in every `runs/<unit>.json`.
- **Software**: rasputin `worktree-15c-2` at `0b42f5b`. `_core` was built
  Release by `tools/bench.py` into that tree's `build-bench/` (sha256
  `41a22a6168db…250485f3`, the same `.so` as basin-phases and as the 761
  mesh). It ran from a venv of its own with a `.pth` to `build-bench/pkg`.
  `rasputin fetch` reads the `rasputin` distribution's version for its
  manifest, so the venv carries a `METADATA`-only `rasputin-0.2.0.dev0`
  dist-info stub (the same stub as basin-phases' venv). 10 threads.
- **Run**: `run_level3.py` handles one unit at a time. First it estimates the
  unit's peak footprint from basin-phases' per-unit figures (the largest of
  load, resample and store, all with `MallocLargeCache=0`) and skips any unit
  whose estimate is above 24 GB. Then it runs `rasputin fetch anadem-v1
  --domain … --out-crs …`, and then `MallocLargeCache=0 /usr/bin/time -l
  rasputin mesh --dem anadem-v1 --cache … --domain … --out-crs … --tolerance
  20`. Swap is polled once a second during the mesh run, and the run is killed
  if swap grows by 3 GB (this never happened). `tabulate.py` builds the table
  below.
- **Domain, DEM, tolerance**: BHO 2017 50k level-3 outlines (`fetch_level3.py`
  and `runs/level3_units.json` in `../basin-phases/`). ANADEM v1 from the
  cache, with all blocks present for every unit (`runs/<unit>.fetch.out`).
  Resampled to a 30 m grid in the output CRS. Tolerance 20 m.
- **Output CRS**: the basin box's suggested `--out-crs`,
  `../basin-phases/runs/basin_out_crs.wkt`. This is a Transverse Mercator on
  SIRGAS 2000, central meridian 42° W, scale factor 0.997548, false easting and
  northing 0. All nine files carry the same `crs` field.
- **Unit 761** was not re-run. Its figures come from the earlier 20 m run with
  the same software, settings and power state: `runs/761.log` (the `time -l`
  output) and `runs/761.stats.md`.
- One run per unit. A test agent was running on the machine at the same time,
  so **the wall times only place the runs and are not timings.** Memory and
  counts are unaffected by it.

## Per unit, 20 m

"Box nodes" is the unit's target grid (rows × cols at 30 m). "Check points"
is the number of source nodes the final check tested. "Final-check inserted"
is phase 2's insertions, with the number of rounds in brackets. The peak is
`/usr/bin/time -l`'s peak memory footprint. Worst angle, "< 10°" and max
degree are plan-view figures from `--stats`.

| unit | area km² | box nodes | check points | triangles | vertices | final-check inserted (rounds) | peak footprint GB | wall s | worst angle | < 10° | max degree | .vtk MB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 761 | 209,316 | 538,409,322 | 270,407,394 | 1,094,297 | 565,650 | 16,715 (11) | 13.88 | 48.9 | 0.00584° | 0.59 % | 17 | 59.0 |
| 762 | 77,174 | 148,938,480 | 104,287,612 | 309,862 | 164,190 | 5,034 (9) | 5.77 | 13.0 | 0.0129° | 1.04 % | 16 | 16.9 |
| 763 | 42,134 | 92,050,400 | 58,985,817 | 451,093 | 232,215 | 9,102 (8) | 3.81 | 9.2 | 0.147° | 0.43 % | 15 | 24.0 |
| 764 | 34,243 | 81,616,394 | 48,440,852 | 115,578 | 63,417 | 1,141 (8) | 4.47 | 7.2 | 0.0221° | 0.68 % | 16 | 6.3 |
| 765 | 37,876 | 114,786,552 | 55,122,267 | 193,742 | 104,086 | 2,680 (8) | 5.65 | 9.2 | 0.551° | 0.62 % | 15 | 10.4 |
| 766 | 30,526 | 63,760,968 | 46,225,333 | 302,211 | 157,263 | 5,530 (8) | 2.82 | 6.6 | 0.0255° | 0.28 % | 14 | 16.1 |
| 767 | 50,881 | 100,694,256 | 71,857,853 | 337,536 | 176,431 | 5,522 (9) | 5.41 | 10.2 | 0.00511° | 0.44 % | 14 | 18.2 |
| 768 | 45,080 | 107,687,349 | 65,762,024 | 425,648 | 220,529 | 8,361 (9) | 4.20 | 10.4 | 0.00787° | 0.41 % | 16 | 22.9 |
| 769 | 106,394 | 195,552,240 | 143,885,192 | 1,760,731 | 891,585 | 34,946 (10) | 5.52 | 25.7 | 0.00994° | 0.21 % | 16 | 94.6 |
| **total** | 633,624 | 1,443,495,961 | 864,974,344 | **4,990,698** | 2,575,366 | 89,031 | max 13.88 | 140.4 | | | | 268.3 |

The basin at 20 m, as nine separate meshes, has **4,990,698 triangles**. The
achieved maximum error was 19.9988-20 m in every unit (`runs/<unit>.stats.md`,
Refinement). No unit was skipped. No run dropped vertices for missing data,
and every run exited 0.

No unit came close to the memory limit. The largest peak was 761's 13.9 GB;
the other eight peaked between 2.8 and 5.8 GB. The basin-phases estimate
(`estimate_gb` in `runs/<unit>.json`) was within 7 % of the measured peak in
all nine units. For 761 it was 0.8 % high, and for the other eight it was
1-7 % low.
In these units the resample transient (135 B × 256 rows × cols × 10 threads)
is the term that sets the peak.

The check points summed over the nine units come to 865 M for 633.6 k km²,
or 1,365 per km². basin-phases measured 1,292 per km² on 761 alone.

## Where the meshes are

They are in `../rasputin_data/sao_francisco_piece/meshes/level3/` (762-769),
with the 761 mesh one level up. They are not in the repository. That folder's
`README.md` is the index: units, CRS, counts, credits and regenerate commands.
Their sha256 values are in `runs/meshes.sha256`.

## Reproduce

From this directory, with `PY` a python whose environment imports the Release
`tin_engine` and has `rasputin` distribution metadata, `D` the
`rasputin_data` folder, and `R` a runs folder:

```
$PY run_level3.py $R $D/sao_francisco_piece/meshes/level3 ../basin-phases/runs/level3_units.json \
  $D/sao_francisco_piece/bho2017_50k_level3 ../basin-phases/runs/basin_out_crs.wkt $D/cache 20 \
  762 763 764 765 766 767 768 769
$PY tabulate.py $R ../basin-phases/runs/level3_units.json $D/sao_francisco_piece/meshes \
  $D/sao_francisco_piece/meshes/level3
```

## Files

- `run_level3.py`: estimate, fetch and mesh per unit → `runs/<unit>.{json,fetch.out,log,stats.md}`.
- `tabulate.py`: the table above → `runs/table.{json,md}`.
- `runs/761.{log,stats.md}`: copied from the earlier 761 run.
- `runs/meshes.sha256`: the nine meshes.
