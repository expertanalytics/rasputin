Paired boxes: 200 (same seed, same centres).

| tolerance | piece triangles, GLO-30 | piece triangles, ANADEM | ANADEM / GLO-30 | basin triangles, GLO-30 | basin triangles, ANADEM | ANADEM / GLO-30 | paired box ratio, median (p10-p90) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 10,431,955 | 8,699,303 | 0.834 | 237.2 M | 174.8 M | 0.737 | 0.69 (0.54-0.83) |
| 2 m | 5,736,406 | 4,584,784 | 0.799 | 113.5 M | 84.2 M | 0.742 | 0.72 (0.58-0.83) |
| 5 m | 1,993,312 | 1,576,681 | 0.791 | 35.9 M | 27.6 M | 0.767 | 0.78 (0.65-0.86) |
| 10 m | 807,337 | 649,845 | 0.805 | 13.8 M | 11.0 M | 0.802 | 0.81 (0.72-0.91) |
| 20 m | 314,507 | 264,833 | 0.842 | 4.9 M | 4.2 M | 0.852 | 0.88 (0.77-1.00) |
| 50 m | 82,517 | 74,841 | 0.907 | 1.0 M | 0.9 M | 0.903 | 0.94 (0.79-1.04) |

The final check's excess: source-DEM nodes (each DEM's own) off the mesh by more than the tolerance, interior only. Piece: share of its interior source nodes. Basin: the boxes' mean per km² times the basin's area, and as a share of the basin's vertices.

| tolerance | piece, GLO-30 | piece, ANADEM | basin nodes, GLO-30 | basin nodes, ANADEM | share of vertices, GLO-30 | share of vertices, ANADEM |
|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 22.30% (max 27.1 m) | 16.64% (max 31.7 m) | 52.44 M | 32.96 M | 44% | 38% |
| 2 m | 8.44% (max 28.6 m) | 5.52% (max 31.7 m) | 17.03 M | 9.92 M | 30% | 24% |
| 5 m | 1.42% (max 28.6 m) | 0.85% (max 31.7 m) | 2.62 M | 1.45 M | 15% | 10% |
| 10 m | 0.31% (max 28.4 m) | 0.18% (max 31.7 m) | 0.53 M | 0.30 M | 8% | 5% |
| 20 m | 0.06% (max 32.4 m) | 0.03% (max 27.9 m) | 0.12 M | 0.07 M | 5% | 3% |
| 50 m | 0.01% (max 55.6 m) | 0.00% (max 56.6 m) | 0.01 M | 0.01 M | 3% | 2% |

The boundary strip on the piece (within 30 m of the outline), worst source node:

| tolerance | GLO-30 | ANADEM |
|---:|---:|---:|
| 1 m | 22.1 m (30.0% over) | 16.7 m (24.5% over) |
| 2 m | 22.1 m (14.3% over) | 15.8 m (10.4% over) |
| 5 m | 23.6 m (2.9% over) | 15.8 m (1.7% over) |
| 10 m | 17.9 m (0.4% over) | 14.1 m (0.3% over) |
| 20 m | 24.3 m (0.0% over) | 20.6 m (0.0% over) |
| 50 m | 38.5 m (0.0% over) | 36.0 m (0.0% over) |

Relief inside a box, median: GLO-30 192 m, ANADEM 191 m.
