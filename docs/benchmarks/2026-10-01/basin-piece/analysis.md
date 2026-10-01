Basin (BHO level 2, ottobasin 76): 635,194.5 km² geodesic; OAS: 636,920 km².

## The piece: velhas76949, 11,667.6 km² (geodesic)

Medians of 5 timed runs per tolerance; one more ASCII run for quality and the final check.

| tolerance | triangles | vertices | vertices / grid nodes | triangles / km² | refine | process wall | max RSS | threads |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 10,431,955 | 5,220,159 | 40.3% | 894.1 | 6.33 s | 9.67 s | 3.60 GiB | 10 |
| 2 m | 5,736,406 | 2,872,292 | 22.2% | 491.7 | 3.36 s | 5.54 s | 2.18 GiB | 10 |
| 5 m | 1,993,312 | 1,000,480 | 7.7% | 170.8 | 1.15 s | 2.38 s | 1.01 GiB | 10 |
| 10 m | 807,337 | 407,376 | 3.1% | 69.2 | 0.50 s | 1.44 s | 0.79 GiB | 10 |
| 20 m | 314,507 | 160,916 | 1.2% | 27.0 | 0.24 s | 1.03 s | 0.70 GiB | 10 |
| 50 m | 82,517 | 44,914 | 0.3% | 7.1 | 0.10 s | 0.83 s | 0.62 GiB | 10 |

| tolerance | worst angle | max degree | Delaunay violations (checked, ambiguous) | within tolerance | control: grid nodes over tol | source nodes over tol, interior | source nodes over tol, boundary strip |
|---:|---:|---:|---:|---|---:|---:|---:|
| 1 m | 0.113° | 14 | 0 (15,643,752, 2,301,997) | True | 0 of 12,957,257 | 2,894,624 of 12,981,642 (22.30%; max 27.1 m) | 7,892 of 26,294 (30.01%; max 22.1 m) |
| 2 m | 0.113° | 13 | 0 (8,600,521, 849,114) | True | 0 of 12,957,257 | 1,095,911 of 12,981,642 (8.44%; max 28.6 m) | 3,772 of 26,294 (14.35%; max 22.1 m) |
| 5 m | 0.107° | 14 | 0 (2,986,145, 109,047) | True | 0 of 12,957,257 | 184,961 of 12,981,642 (1.42%; max 28.6 m) | 773 of 26,294 (2.94%; max 23.6 m) |
| 10 m | 0.113° | 14 | 0 (1,207,299, 14,356) | True | 0 of 12,957,257 | 39,884 of 12,981,642 (0.31%; max 28.4 m) | 113 of 26,294 (0.43%; max 17.9 m) |
| 20 m | 0.113° | 14 | 0 (468,099, 1,701) | True | 0 of 12,957,257 | 7,908 of 12,981,642 (0.06%; max 32.4 m) | 11 of 26,294 (0.04%; max 24.3 m) |
| 50 m | 0.113° | 14 | 0 (120,121, 91) | True | 0 of 12,957,257 | 858 of 12,981,642 (0.01%; max 55.6 m) | 0 of 26,294 (0.00%; max 38.5 m) |

Power: AC before and after every tolerance block: True.

Max RSS against triangles over the sweep, least squares: 0.55 GiB + 310 B per triangle (max residual 121 MiB).

## The basin, extrapolated from 200 random 10 km boxes (seed 1)

Mean density over the boxes times the basin's area; the interval is the bootstrap 95 % interval of the mean. Memory floor: the piece's per-triangle RSS slope (310 B) times the basin's triangles, plus the 8.6 GiB canvas of Q9. Refine time: the piece's refine time per triangle (10 threads), times the basin's triangles, i.e. linear.

| tolerance | triangles / km², mean (95 %) | median | p10-p90 | basin triangles (95 %) | basin vertices | refine, linear | memory floor |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 373.4 (338.1-410.9) | 315.6 | 101.3-757.0 | 237.2 M (214.8-261.0 M) | 118.6 M | 2.4 min | 77 GiB |
| 2 m | 178.6 (158.4-199.6) | 138.5 | 40.3-391.5 | 113.5 M (100.6-126.8 M) | 56.7 M | 1.1 min | 41 GiB |
| 5 m | 56.6 (49.2-64.6) | 39.6 | 8.4-139.3 | 35.9 M (31.3-41.0 M) | 18.0 M | 0.3 min | 19 GiB |
| 10 m | 21.7 (18.5-24.9) | 13.6 | 2.2-53.4 | 13.8 M (11.8-15.8 M) | 6.9 M | 0.1 min | 13 GiB |
| 20 m | 7.7 (6.5-9.1) | 4.3 | 0.1-19.9 | 4.9 M (4.1-5.8 M) | 2.5 M | 0.1 min | 10 GiB |
| 50 m | 1.6 (1.3-2.0) | 0.7 | 0.0-4.6 | 1.0 M (0.8-1.2 M) | 0.5 M | 0.0 min | 9 GiB |

The final check of Q6, measured: source-DEM nodes outside the tolerance after meshing the resampled grid, **interior only** (a box's 4-corner edges make its boundary strip unrepresentative of a real outline; the piece gives the strip). Basin count = mean per km² times the basin's area.

| tolerance | share of source nodes | per km², mean (95 %) | basin source nodes over tol | as a share of the basin's vertices |
|---:|---:|---:|---:|---:|
| 1 m | 7.71% | 82.6 (72.5-93.2) | 52.44 M | 44% |
| 2 m | 2.50% | 26.8 (22.9-31.1) | 17.03 M | 30% |
| 5 m | 0.39% | 4.1 (3.4-4.9) | 2.62 M | 15% |
| 10 m | 0.08% | 0.8 (0.7-1.0) | 0.53 M | 8% |
| 20 m | 0.02% | 0.2 (0.1-0.2) | 0.12 M | 5% |
| 50 m | 0.00% | 0.0 (0.0-0.0) | 0.01 M | 3% |

Relief inside a box (z max - z min of its grid): median 192 m, p10 55 m, p90 413 m.

## Check of the box method: 30 random 10 km boxes inside the piece against the piece meshed whole

| tolerance | piece, whole catchment, triangles / km² | boxes inside it, mean (95 %) | inside the interval |
|---:|---:|---:|---|
| 1 m | 894.1 | 934.2 (831.3-1,037.5) | True |
| 2 m | 491.7 | 516.5 (445.6-588.9) | True |
| 5 m | 170.8 | 180.8 (151.7-211.3) | True |
| 10 m | 69.2 | 73.7 (60.7-86.9) | True |
| 20 m | 27.0 | 28.7 (23.2-34.4) | True |
| 50 m | 7.1 | 6.5 (4.8-8.2) | True |

Relief inside a box: median 285 m.

## Refine thread scaling on the piece at 1 m

| threads | refine (median) | speed-up | process wall (median) |
|---:|---:|---:|---:|
| 1 | 12.95 s | 1.00 | 14.12 s |
| 2 | 9.12 s | 1.42 | 10.31 s |
| 4 | 7.30 s | 1.77 | 8.49 s |
| 6 | 6.69 s | 1.93 | 7.87 s |
| 8 | 6.36 s | 2.03 | 7.55 s |
| 10 | 6.32 s | 2.05 | 7.50 s |

Power: AC before and after: True.

