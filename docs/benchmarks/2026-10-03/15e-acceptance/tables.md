
### `761_t10_base_r1`

`time -l` peak memory footprint 15.36 GB, wall 54.3 s, 205 samples, power: AC Power

| phase | s | footprint at start GB | peak GB | rise over start GB |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 4.6 | 0.06 | 6.84 | 6.79 |
| P2 resample (blocks) | 19.5 | 6.84 | 15.36 | 8.51 |
| P2 end (DemTile copy) | 0.5 | 8.00 | 10.15 (sampled) | 2.15 |
| start mesh | 0.2 | 10.13 | 10.13 (sampled) | 0.00 |
| P3 refine phase 1 | 4.5 | 7.61 | 7.64 (sampled) | 0.03 |
| P4 store fill + freeze | 13.9 | 5.12 | 11.11 (sampled) | 5.99 |
| P5 phase 2 (refine_points) | 9.5 | 11.11 | 11.75 (sampled) | 0.65 |
| P6 trim + write | 0.8 | 7.66 | 7.89 (sampled) | 0.24 |

| measure | value |
|---|---:|
| source mosaic nodes | 560,943,760 |
| target grid nodes | 538,409,322 |
| check points | 270,407,394 |
| phase 1 triangles (2 V - hull) | 2,772,536 |
| final triangles | 2,940,024 |
| P1 rise per mosaic node, B | 12.1 |
| P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols) | 89.4 |
| P2 end rise per grid node, B (copy) | 6.1 |
| P4 rise per check point, B | 22.2 |
| P4 held after freeze per check point, B | 22.2 |
| P3 rise per phase-1 triangle, B | 9.6 |
| P5 rise per final triangle, B | 220.3 |

### `761_t10_base_r2`

`time -l` peak memory footprint 15.91 GB, wall 51.8 s, 197 samples, power: AC Power

| phase | s | footprint at start GB | peak GB | rise over start GB |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 3.4 | 0.06 | 7.41 | 7.35 |
| P2 resample (blocks) | 19.6 | 7.41 | 15.91 | 8.50 |
| P2 end (DemTile copy) | 0.3 | 8.19 | 10.35 (sampled) | 2.15 |
| start mesh | 0.1 | 10.29 | 10.29 (sampled) | 0.00 |
| P3 refine phase 1 | 4.1 | 10.25 | 10.28 (sampled) | 0.03 |
| P4 store fill + freeze | 14.0 | 9.90 | 11.72 (sampled) | 1.81 |
| P5 phase 2 (refine_points) | 8.8 | 11.72 | 12.36 (sampled) | 0.65 |
| P6 trim + write | 0.8 | 7.71 | 7.80 (sampled) | 0.10 |

| measure | value |
|---|---:|
| source mosaic nodes | 560,943,760 |
| target grid nodes | 538,409,322 |
| check points | 270,407,394 |
| phase 1 triangles (2 V - hull) | 2,772,536 |
| final triangles | 2,940,024 |
| P1 rise per mosaic node, B | 13.1 |
| P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols) | 89.2 |
| P2 end rise per grid node, B (copy) | 5.5 |
| P4 rise per check point, B | 6.7 |
| P4 held after freeze per check point, B | 6.7 |
| P3 rise per phase-1 triangle, B | 9.7 |
| P5 rise per final triangle, B | 220.0 |

### `761_t10_new_r1`

`time -l` peak memory footprint 10.13 GB, wall 45.2 s, 172 samples, power: AC Power

| phase | s | footprint at start GB | peak GB | rise over start GB |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 2.3 | 0.06 | 7.41 | 7.35 |
| P2 resample (blocks) | 13.5 | 7.41 | 8.76 | 1.36 |
| P2 end (DemTile copy) | 0.0 | 5.19 | 5.19 (sampled) | 0.00 |
| start mesh | 0.1 | 5.19 | 5.19 (sampled) | 0.01 |
| P3 refine phase 1 | 4.3 | 5.19 | 5.19 (sampled) | 0.00 |
| P4 store fill + freeze | 14.3 | 5.13 | 9.76 | 4.62 |
| P5 phase 2 (refine_points) | 9.4 | 9.76 | 10.13 | 0.38 |
| P6 trim + write | 0.9 | 5.42 | 5.42 (sampled) | 0.00 |

| measure | value |
|---|---:|
| source mosaic nodes | 560,943,760 |
| target grid nodes | 538,409,322 |
| check points | 270,407,394 |
| phase 1 triangles (2 V - hull) | 2,772,536 |
| final triangles | 2,940,024 |
| P1 rise per mosaic node, B | 13.1 |
| P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols) | -77.4 |
| P2 end rise per grid node, B (copy) | -4.1 |
| P4 rise per check point, B | 17.1 |
| P4 held after freeze per check point, B | 17.1 |
| P3 rise per phase-1 triangle, B | 0.0 |
| P5 rise per final triangle, B | 128.1 |

### `761_t10_new_r2`

`time -l` peak memory footprint 10.12 GB, wall 44.0 s, 167 samples, power: AC Power

| phase | s | footprint at start GB | peak GB | rise over start GB |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 2.3 | 0.06 | 7.40 | 7.35 |
| P2 resample (blocks) | 13.0 | 7.40 | 8.82 | 1.42 |
| P2 end (DemTile copy) | 0.0 | 5.18 | 5.18 (sampled) | 0.00 |
| start mesh | 0.1 | 5.18 | 5.19 (sampled) | 0.01 |
| P3 refine phase 1 | 4.1 | 5.19 | 5.19 (sampled) | 0.00 |
| P4 store fill + freeze | 14.0 | 5.13 | 9.75 | 4.62 |
| P5 phase 2 (refine_points) | 9.0 | 9.75 | 10.12 | 0.37 |
| P6 trim + write | 0.9 | 5.41 | 5.41 (sampled) | 0.00 |

| measure | value |
|---|---:|
| source mosaic nodes | 560,943,760 |
| target grid nodes | 538,409,322 |
| check points | 270,407,394 |
| phase 1 triangles (2 V - hull) | 2,772,536 |
| final triangles | 2,940,024 |
| P1 rise per mosaic node, B | 13.1 |
| P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols) | -71.5 |
| P2 end rise per grid node, B (copy) | -4.1 |
| P4 rise per check point, B | 17.1 |
| P4 held after freeze per check point, B | 17.1 |
| P3 rise per phase-1 triangle, B | 0.0 |
| P5 rise per final triangle, B | 125.8 |

### `761_t10_base_nolargecache_r1`

`time -l` peak memory footprint 14.26 GB, wall 59.8 s, 226 samples, power: AC Power

| phase | s | footprint at start GB | peak GB | rise over start GB |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 2.5 | 0.05 | 6.81 | 6.76 |
| P2 resample (blocks) | 25.7 | 2.32 | 14.26 | 11.94 |
| P2 end (DemTile copy) | 0.6 | 4.48 | 6.64 (sampled) | 2.15 |
| start mesh | 0.1 | 4.48 | 4.49 (sampled) | 0.01 |
| P3 refine phase 1 | 4.3 | 4.49 | 4.70 (sampled) | 0.22 |
| P4 store fill + freeze | 15.7 | 4.55 | 6.15 (sampled) | 1.59 |
| P5 phase 2 (refine_points) | 9.2 | 2.55 | 3.11 (sampled) | 0.55 |
| P6 trim + write | 1.1 | 2.63 | 2.82 (sampled) | 0.19 |

| measure | value |
|---|---:|
| source mosaic nodes | 560,943,760 |
| target grid nodes | 538,409,322 |
| check points | 270,407,394 |
| phase 1 triangles (2 V - hull) | 2,772,536 |
| final triangles | 2,940,024 |
| P1 rise per mosaic node, B | 4.0 |
| P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols) | 137.6 |
| P2 end rise per grid node, B (copy) | 8.0 |
| P4 rise per check point, B | 5.9 |
| P4 held after freeze per check point, B | -7.4 |
| P3 rise per phase-1 triangle, B | 78.8 |
| P5 rise per final triangle, B | 188.0 |

### `761_t10_base_nolargecache_r2`

`time -l` peak memory footprint 14.49 GB, wall 57.3 s, 215 samples, power: AC Power

| phase | s | footprint at start GB | peak GB | rise over start GB |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 2.4 | 0.05 | 6.81 | 6.76 |
| P2 resample (blocks) | 24.1 | 2.32 | 14.49 | 12.17 |
| P2 end (DemTile copy) | 0.5 | 4.48 | 6.63 (sampled) | 2.15 |
| start mesh | 0.1 | 4.48 | 4.49 (sampled) | 0.01 |
| P3 refine phase 1 | 4.1 | 4.48 | 4.75 (sampled) | 0.27 |
| P4 store fill + freeze | 15.7 | 4.55 | 6.16 (sampled) | 1.61 |
| P5 phase 2 (refine_points) | 8.7 | 2.57 | 3.12 (sampled) | 0.55 |
| P6 trim + write | 1.0 | 2.64 | 2.82 (sampled) | 0.18 |

| measure | value |
|---|---:|
| source mosaic nodes | 560,943,760 |
| target grid nodes | 538,409,322 |
| check points | 270,407,394 |
| phase 1 triangles (2 V - hull) | 2,772,536 |
| final triangles | 2,940,024 |
| P1 rise per mosaic node, B | 4.0 |
| P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols) | 140.8 |
| P2 end rise per grid node, B (copy) | 8.0 |
| P4 rise per check point, B | 6.0 |
| P4 held after freeze per check point, B | -7.3 |
| P3 rise per phase-1 triangle, B | 97.5 |
| P5 rise per final triangle, B | 188.1 |

### `761_t10_new_nolargecache_r1`

`time -l` peak memory footprint 6.93 GB, wall 51.3 s, 195 samples, power: AC Power

| phase | s | footprint at start GB | peak GB | rise over start GB |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 2.4 | 0.05 | 6.81 | 6.76 |
| P2 resample (blocks) | 18.1 | 2.32 | 5.32 (sampled) | 3.00 |
| P2 end (DemTile copy) | 0.0 | 4.48 | 4.48 (sampled) | 0.00 |
| start mesh | 0.1 | 4.48 | 4.48 (sampled) | 0.00 |
| P3 refine phase 1 | 4.4 | 4.48 | 4.76 (sampled) | 0.28 |
| P4 store fill + freeze | 15.6 | 2.39 | 6.93 | 4.53 |
| P5 phase 2 (refine_points) | 9.2 | 4.65 | 5.20 (sampled) | 0.55 |
| P6 trim + write | 1.0 | 4.72 | 4.72 (sampled) | 0.00 |

| measure | value |
|---|---:|
| source mosaic nodes | 560,943,760 |
| target grid nodes | 538,409,322 |
| check points | 270,407,394 |
| phase 1 triangles (2 V - hull) | 2,772,536 |
| final triangles | 2,940,024 |
| P1 rise per mosaic node, B | 4.0 |
| P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols) | 82.2 |
| P2 end rise per grid node, B (copy) | 4.0 |
| P4 rise per check point, B | 16.8 |
| P4 held after freeze per check point, B | 8.4 |
| P3 rise per phase-1 triangle, B | 102.1 |
| P5 rise per final triangle, B | 187.7 |

### `761_t10_new_nolargecache_r2`

`time -l` peak memory footprint 6.92 GB, wall 50.3 s, 192 samples, power: AC Power

| phase | s | footprint at start GB | peak GB | rise over start GB |
|---|---:|---:|---:|---:|
| P1 decode + assemble | 2.4 | 0.05 | 6.81 | 6.76 |
| P2 resample (blocks) | 18.4 | 2.32 | 5.61 (sampled) | 3.29 |
| P2 end (DemTile copy) | 0.0 | 4.48 | 4.48 (sampled) | 0.00 |
| start mesh | 0.1 | 4.48 | 4.48 (sampled) | 0.00 |
| P3 refine phase 1 | 4.0 | 4.48 | 4.76 (sampled) | 0.28 |
| P4 store fill + freeze | 15.1 | 2.39 | 6.92 | 4.53 |
| P5 phase 2 (refine_points) | 8.9 | 4.65 | 5.20 (sampled) | 0.55 |
| P6 trim + write | 1.0 | 4.72 | 4.72 (sampled) | 0.00 |

| measure | value |
|---|---:|
| source mosaic nodes | 560,943,760 |
| target grid nodes | 538,409,322 |
| check points | 270,407,394 |
| phase 1 triangles (2 V - hull) | 2,772,536 |
| final triangles | 2,940,024 |
| P1 rise per mosaic node, B | 4.0 |
| P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols) | 110.2 |
| P2 end rise per grid node, B (copy) | 4.0 |
| P4 rise per check point, B | 16.8 |
| P4 held after freeze per check point, B | 8.4 |
| P3 rise per phase-1 triangle, B | 102.1 |
| P5 rise per final triangle, B | 187.7 |
