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
  run, in `runs/t20/<unit>.json`). Swap stood at about 17.6 GB in use before the
  runs, left over from earlier work. No run added to it: `swap_max_mb` equals
  `swap_before_mb` in every `runs/t20/<unit>.json`.
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
  cache, with all blocks present for every unit (`runs/t20/<unit>.fetch.out`).
  Resampled to a 30 m grid in the output CRS. Tolerance 20 m.
- **Output CRS**: the basin box's suggested `--out-crs`,
  `../basin-phases/runs/basin_out_crs.wkt`. This is a Transverse Mercator on
  SIRGAS 2000, central meridian 42° W, scale factor 0.997548, false easting and
  northing 0. All nine files carry the same `crs` field.
- **Unit 761** was not re-run. Its figures come from the earlier 20 m run with
  the same software, settings and power state: `runs/t20/761.log` (the `time -l`
  output) and `runs/t20/761.stats.md`.
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
achieved maximum error was 19.9988-20 m in every unit (`runs/t20/<unit>.stats.md`,
Refinement). No unit was skipped. No run dropped vertices for missing data,
and every run exited 0.

No unit came close to the memory limit. The largest peak was 761's 13.9 GB;
the other eight peaked between 2.8 and 5.8 GB. The basin-phases estimate
(`estimate_gb` in `runs/t20/<unit>.json`) was within 7 % of the measured peak in
all nine units. For 761 it was 0.8 % high, and for the other eight it was
1-7 % low.
In these units the resample transient (135 B × 256 rows × cols × 10 threads)
is the term that sets the peak.

The check points summed over the nine units come to 865 M for 633.6 k km²,
or 1,365 per km². basin-phases measured 1,292 per km² on 761 alone.

## Per unit, 10 m

The same nine units with the same build, CRS, `MallocLargeCache=0` and skip
rule, all nine run this time (761 included). Records in `runs/t10/`. Every
mesh is text VTK. AC power before and after every run, and no run added to swap.

| unit | area km² | box nodes | check points | triangles | vertices | final-check inserted (rounds) | peak footprint GB | wall s | worst angle | < 10° | max degree | .vtk MB | format |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 761 | 209,316 | 538,409,322 | 270,407,394 | 2,940,024 | 1,488,537 | 83,744 (13) | 14.38 | 62.8 | 0.00584° | 0.41 % | 17 | 159.8 | ascii |
| 762 | 77,174 | 148,938,480 | 104,287,612 | 741,469 | 379,995 | 24,387 (11) | 5.67 | 15.4 | 0.0129° | 0.74 % | 16 | 40.1 | ascii |
| 763 | 42,134 | 92,050,400 | 58,985,817 | 1,140,594 | 576,976 | 44,080 (18) | 3.68 | 12.5 | 0.147° | 0.33 % | 15 | 61.0 | ascii |
| 764 | 34,243 | 81,616,394 | 48,440,852 | 275,421 | 143,340 | 5,709 (9) | 4.42 | 8.1 | 0.0221° | 0.46 % | 16 | 14.7 | ascii |
| 765 | 37,876 | 114,786,552 | 55,122,267 | 473,833 | 244,132 | 12,864 (15) | 5.52 | 10.8 | 0.551° | 0.52 % | 16 | 25.4 | ascii |
| 766 | 30,526 | 63,760,968 | 46,225,333 | 787,486 | 399,906 | 27,393 (9) | 2.79 | 9.1 | 0.00586° | 0.29 % | 14 | 42.2 | ascii |
| 767 | 50,881 | 100,694,256 | 71,857,853 | 892,431 | 453,879 | 28,800 (10) | 5.37 | 13.5 | 0.00511° | 0.42 % | 15 | 48.0 | ascii |
| 768 | 45,080 | 107,687,349 | 65,762,024 | 1,145,198 | 580,307 | 42,516 (10) | 4.15 | 14.4 | 0.00787° | 0.40 % | 14 | 61.5 | ascii |
| 769 | 106,394 | 195,552,240 | 143,885,192 | 4,688,908 | 2,355,725 | 180,495 (13) | 5.53 | 40.5 | 0.00994° | 0.20 % | 16 | 260.5 | ascii |
| **total** | 633,624 | 1,443,495,961 | 864,974,344 | **13,085,364** | 6,622,797 | 449,988 | max 14.38 | 187.1 |  |  |  | 713.1 |  |

The basin at 10 m, as nine separate meshes, has **13,085,364 triangles**.
The achieved maximum error was 9.99986-10 m in every unit. The largest peak
was 761's 14.38 GB; the other eight peaked at 2.8-5.7 GB, within 0.5 GB of
their 20 m peaks. Unit 761 at 10 m gave 2,940,024 triangles and 83,744
final-check insertions, the same counts as basin-phases' `761_t10` runs.

## Per unit, 5 m

All nine units again, with the same build, CRS, `MallocLargeCache=0` and skip
rule. Records are in `runs/t5/`. AC power before and after every run, and no
run added to swap. Unit 769 is written as **binary** VTK (`--binary=769`)
because the 10 m figures predicted it would pass 600 MB as text: 4.69 M
triangles × 2.65 (761's 5 m / 10 m ratio) × 54 B per triangle as text is
about 675 MB. As binary it is 449.5 MB. The other eight are text.

| unit | area km² | box nodes | check points | triangles | vertices | final-check inserted (rounds) | peak footprint GB | wall s | worst angle | < 10° | max degree | .vtk MB | format |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 761 | 209,316 | 538,409,322 | 270,407,394 | 7,791,814 | 3,914,535 | 393,389 (13) | 14.27 | 84.4 | 0.00584° | 0.35 % | 18 | 434.3 | ascii |
| 762 | 77,174 | 148,938,480 | 104,287,612 | 1,891,531 | 955,035 | 114,135 (11) | 5.71 | 20.0 | 0.0129° | 0.66 % | 16 | 102.7 | ascii |
| 763 | 42,134 | 92,050,400 | 58,985,817 | 2,948,812 | 1,481,154 | 201,853 (12) | 3.72 | 20.1 | 0.00699° | 0.38 % | 15 | 162.0 | ascii |
| 764 | 34,243 | 81,616,394 | 48,440,852 | 682,531 | 346,900 | 28,249 (10) | 4.48 | 9.5 | 0.0221° | 0.41 % | 16 | 36.6 | ascii |
| 765 | 37,876 | 114,786,552 | 55,122,267 | 1,250,434 | 632,445 | 63,331 (17) | 5.74 | 14.0 | 0.0565° | 0.45 % | 15 | 67.4 | ascii |
| 766 | 30,526 | 63,760,968 | 46,225,333 | 2,088,331 | 1,050,369 | 130,191 (13) | 2.81 | 13.6 | 0.00586° | 0.30 % | 14 | 113.2 | ascii |
| 767 | 50,881 | 100,694,256 | 71,857,853 | 2,378,326 | 1,196,833 | 142,869 (13) | 5.42 | 19.0 | 0.00511° | 0.40 % | 15 | 129.7 | ascii |
| 768 | 45,080 | 107,687,349 | 65,762,024 | 3,078,437 | 1,546,932 | 210,809 (16) | 4.14 | 21.6 | 0.00347° | 0.40 % | 14 | 169.9 | ascii |
| 769 | 106,394 | 195,552,240 | 143,885,192 | 12,466,005 | 6,244,412 | 887,034 (14) | 5.52 | 56.2 | 0.00994° | 0.29 % | 16 | 449.5 | binary |
| **total** | 633,624 | 1,443,495,961 | 864,974,344 | **34,576,221** | 17,368,615 | 2,171,860 | max 14.27 | 258.5 |  |  |  | 1665.4 |  |

The basin at 5 m, as nine separate meshes, has **34,576,221 triangles**. The
achieved maximum error was 4.99999-5 m in every unit. The largest peak was
761's 14.27 GB, and the other eight peaked at 2.8-5.7 GB. Over 20, 10 and 5 m
each unit's peak moved by at most 0.5 GB, so the peak is set before
refinement, as basin-phases found. Unit 761 at 5 m gave 7,791,814 triangles,
the same count as basin-phases' `761_t5_nolargecache` run.

## Where the meshes are

They are in `../rasputin_data/sao_francisco_piece/meshes/level3/`, named
`sub_basin_<unit>_anadem_tol<t>m.vtk`. The 761 mesh at 20 m sits one level up. They are not in the repository. That folder's
`README.md` is the index: units, CRS, counts, credits and regenerate commands.
Their sha256 values are in `runs/t<t>/meshes.sha256`.

## Reproduce

From this directory, with `PY` a python whose environment imports the Release
`tin_engine` and has `rasputin` distribution metadata, `D` the
`rasputin_data` folder, and `R` a runs folder:

```
$PY run_level3.py $R $D/sao_francisco_piece/meshes/level3 ../basin-phases/runs/level3_units.json \
  $D/sao_francisco_piece/bho2017_50k_level3 ../basin-phases/runs/basin_out_crs.wkt $D/cache 20 \
  762 763 764 765 766 767 768 769
$PY tabulate.py $R ../basin-phases/runs/level3_units.json 20 $D/sao_francisco_piece/meshes/level3 \
  $D/sao_francisco_piece/meshes
# 10 m: R=runs/t10, tolerance 10, all nine units (761 762 … 769)
# 5 m: R=runs/t5, tolerance 5, all nine units, plus --binary=<units> (see the 5 m section)
```

## Files

- `run_level3.py`: estimate, fetch and mesh per unit → `runs/t20/<unit>.{json,fetch.out,log,stats.md}`.
- `tabulate.py`: the tables above → `runs/t<tol>/table.{json,md}`.
- `runs/t10/`, `runs/t5/`: the same records at 10 m and 5 m, with `meshes.sha256`.
- `runs/t20/761.{log,stats.md}`: copied from the earlier 761 run.
- `runs/t20/meshes.sha256`: the nine meshes.
