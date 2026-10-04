### Bygdin, tolerance 10 m

| item | base c193cb1 | 15f-3 83c7fd2 |
|---|---:|---:|
| process wall, median of 4 / 4 (s) | 2.36 (2.25-2.59) | 2.32 (2.28-2.32) |
| peak RSS, max (GiB) | 0.35 | 0.35 |
| output vertices | 42110 | 42349 |
| output triangles | 83171 | 83590 |
| constraint edges | 8882 | 9081 |
| `line_points_checked` | - | 175246 |
| `line_points_on_nodata` | - | 0 |
| `line_points_inserted` | - | 199 |
| `line_check_dem_nodes_inserted` | - | 40 |
| `line_points_refused` | - | 0 |
| `line_points_refused_max_error_m` | - | 0 |
| `line_max_error_m` | - | 9.98675991396044 |
| `line_points_duplicate` | - | 1 |
| `max_error_m` | 9.999732907715497 | 9.999732907715497 |
| `refinement_rounds` | 17 | 17 |
| `points_inserted` | 17879 | 17879 |
| `points_snapped_to_lines` | 534 | 534 |
| phase: refine (s, median) | 0.049 | 0.050 |
| phase: edge strip: generate (s, median) | 0.000 | 0.006 |
| phase: edge strip: scan (parallel) (s, median) | 0.000 | 0.001 |
| phase: edge strip: split + flip (serial) (s, median) | 0.000 | 0.001 |
| phase: trim (s, median) | 0.003 | 0.003 |
| phase: land cover (s, median) | 0.072 | 0.073 |
| phase: write: encode (s, median) | 0.002 | 0.001 |
| phase: other (s, median) | 0.005 | 0.048 |

| quality: worst_angle | 0.004078799958387991 | 0.00875010258693808 |
| quality: angle_median | 36.86989764584402 | 36.86989764584402 |
| quality: share_under_1 | 0.0006973584542689158 | 0.0005263787534394067 |
| quality: max_degree | 16 | 14 |
| quality: within_tolerance | True | True |
| quality: delaunay_violations | 0 | 0 |
| quality: delaunay_checked | 116398 | 116857 |
| mesh sha256 | `ba93d4220aa0721a` | `c33d2736f1173aa0` |

### Bygdin, tolerance 1 m

| item | base c193cb1 | 15f-3 83c7fd2 |
|---|---:|---:|
| process wall, median of 4 / 4 (s) | 4.32 (4.30-4.33) | 5.81 (5.78-5.85) |
| peak RSS, max (GiB) | 0.82 | 1.11 |
| output vertices | 573758 | 584771 |
| output triangles | 1145434 | 1165056 |
| constraint edges | 14074 | 23086 |
| `line_points_checked` | - | 180436 |
| `line_points_on_nodata` | - | 0 |
| `line_points_inserted` | - | 9012 |
| `line_check_dem_nodes_inserted` | - | 2001 |
| `line_points_refused` | - | 0 |
| `line_points_refused_max_error_m` | - | 0 |
| `line_max_error_m` | - | 0.999842752797349 |
| `line_points_duplicate` | - | 0 |
| `max_error_m` | 1 | 1 |
| `refinement_rounds` | 24 | 24 |
| `points_inserted` | 549527 | 549527 |
| `points_snapped_to_lines` | 5726 | 5726 |
| phase: refine (s, median) | 0.600 | 0.598 |
| phase: edge strip: generate (s, median) | 0.000 | 0.007 |
| phase: edge strip: scan (parallel) (s, median) | 0.000 | 0.008 |
| phase: edge strip: split + flip (serial) (s, median) | 0.000 | 0.026 |
| phase: trim (s, median) | 0.042 | 0.042 |
| phase: land cover (s, median) | 1.295 | 1.316 |
| phase: write: encode (s, median) | 0.017 | 0.017 |
| phase: other (s, median) | 0.005 | 1.420 |

| quality: worst_angle | 0.004078799958387991 | 0.021765796134158345 |
| quality: angle_median | 45.0 | 45.0 |
| quality: share_under_1 | 0.0009961289781864342 | 8.926609536365634e-05 |
| quality: max_degree | 29 | 18 |
| quality: within_tolerance | True | True |
| quality: delaunay_violations | 0 | 3 |
| quality: delaunay_checked | 1705117 | 1726740 |
| mesh sha256 | `1709cbbe8e47fa0c` | `661bf16b8322b405` |

