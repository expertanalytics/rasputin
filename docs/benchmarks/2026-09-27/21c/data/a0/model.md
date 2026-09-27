| domain | option | per-item costs (ns, serial) | steps | syncs | T | split, no sync | split, spawn | split, pool | refine at 8, spawn | refine at 8, pool |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| quarter | A0 | eval 279, commit 434 | 968 | 1448 | 1 | 729.1 | 729.1 | 729.1 |  |  |
| quarter | A0 |  |  |  | 4 | 182.4 | 265.8 | 183.1 |  |  |
| quarter | A0 |  |  |  | 8 | 91.3 | 230.6 | 92.6 | 295 | 157 |
| quarter | A0 |  |  |  | 16 | 45.7 | 273.6 | 75.4 |  |  |
| quarter | A1 | eval 279, commit 434 | 879 | 1341 | 1 | 691.7 | 691.7 | 691.7 |  |  |
| quarter | A1 |  |  |  | 4 | 173.1 | 250.3 | 173.7 |  |  |
| quarter | A1 |  |  |  | 8 | 86.6 | 215.6 | 87.8 | 280 | 152 |
| quarter | A1 |  |  |  | 16 | 43.4 | 254.4 | 70.8 |  |  |
| quarter | C | cbatch 113, ctest 103, cflip 77 | 298 | 447 | 1 | 394.5 | 394.5 | 394.5 |  |  |
| quarter | C |  |  |  | 4 | 98.6 | 124.4 | 98.9 |  |  |
| quarter | C |  |  |  | 8 | 49.3 | 92.3 | 49.7 | 157 | 114 |
| quarter | C |  |  |  | 16 | 24.7 | 95.0 | 33.8 |  |  |
| quarter | A0-once | eval 279, commit 434 | 539 | 1448 | 1 | 303.3 | 303.3 | 303.3 |  |  |
| quarter | A0-once |  |  |  | 4 | 75.9 | 159.3 | 76.6 |  |  |
| quarter | A0-once |  |  |  | 8 | 38.0 | 177.3 | 39.3 | 242 | 104 |
| quarter | A0-once |  |  |  | 16 | 19.1 | 246.9 | 48.7 |  |  |
| quarter | A1-once | eval 279, commit 434 | 467 | 1341 | 1 | 238.7 | 238.7 | 238.7 |  |  |
| quarter | A1-once |  |  |  | 4 | 59.8 | 137.0 | 60.4 |  |  |
| quarter | A1-once |  |  |  | 8 | 29.9 | 158.9 | 31.1 | 223 | 95 |
| quarter | A1-once |  |  |  | 16 | 15.0 | 226.1 | 42.5 |  |  |
| tile | A0 | eval 259, commit 421 | 1215 | 1829 | 1 | 797.3 | 797.3 | 797.3 |  |  |
| tile | A0 |  |  |  | 4 | 199.5 | 304.8 | 200.4 |  |  |
| tile | A0 |  |  |  | 8 | 99.9 | 275.8 | 101.5 | 364 | 190 |
| tile | A0 |  |  |  | 16 | 50.0 | 337.9 | 87.5 |  |  |
| tile | A1 | eval 259, commit 421 | 729 | 1128 | 1 | 637.9 | 637.9 | 637.9 |  |  |
| tile | A1 |  |  |  | 4 | 159.6 | 224.5 | 160.2 |  |  |
| tile | A1 |  |  |  | 8 | 79.9 | 188.4 | 80.9 | 277 | 169 |
| tile | A1 |  |  |  | 16 | 40.0 | 217.5 | 63.1 |  |  |
| tile | C | cbatch 118, ctest 104, cflip 87 | 436 | 654 | 1 | 406.9 | 406.9 | 406.9 |  |  |
| tile | C |  |  |  | 4 | 101.8 | 139.4 | 102.1 |  |  |
| tile | C |  |  |  | 8 | 50.9 | 113.8 | 51.5 | 202 | 140 |
| tile | C |  |  |  | 16 | 25.5 | 128.4 | 38.8 |  |  |
| tile | A0-once | eval 259, commit 421 | 675 | 1829 | 1 | 335.8 | 335.8 | 335.8 |  |  |
| tile | A0-once |  |  |  | 4 | 84.0 | 189.4 | 85.0 |  |  |
| tile | A0-once |  |  |  | 8 | 42.1 | 218.1 | 43.7 | 307 | 132 |
| tile | A0-once |  |  |  | 16 | 21.1 | 309.0 | 58.6 |  |  |
| tile | A1-once | eval 259, commit 421 | 403 | 1128 | 1 | 229.7 | 229.7 | 229.7 |  |  |
| tile | A1-once |  |  |  | 4 | 57.5 | 122.4 | 58.1 |  |  |
| tile | A1-once |  |  |  | 8 | 28.8 | 137.3 | 29.8 | 226 | 118 |
| tile | A1-once |  |  |  | 16 | 14.4 | 192.0 | 37.5 |  |  |

Sync costs used (us, median of runs' medians): T=4: spawn 57.6, for_each_block 58.9, std::barrier 0.50, spin 0.18; T=8: spawn 96.2, for_each_block 85.8, std::barrier 0.90, spin 0.81; T=16: spawn 157.4, for_each_block 156.9, std::barrier 20.46, spin 15.01
Today, measured: split 92.6 ms (quarter) / 92.6 ms (tile) at 1 thread; refine at 8 threads 161.7 / 186.3 ms (21c section 6).
