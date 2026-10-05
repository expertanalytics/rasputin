### Every `uncertain` station, with its causes

Causes from the sensitivity (`causes`): `swing` the area changes over 5 % within the position uncertainty `U` up or down the river; `downstream_unread` the river was not read to `U` below the gauge; `chain_not_draining` a node of the burnt path does not drain into the next, or the counts do not rise; `chain_end_open` the path's lowered end found no way down within the cap. Largest step: the largest rise in area between two path nodes within `U`, and where (m, + downstream). Bypass: in the burnt DEM, the largest count within 30 m of the placed node, against the placed node's own (`bypass.py`).

| station | name | NVE km² | ours km² | NVE's in ours % | causes | largest step km² at m | U m | lake line / lake above P | bypass km² |
|---|---|---:|---:|---:|---|---|---:|---|---:|
| 2.279.0 | Kråkfoss | 434.5 | 0.000 | 0.0 | swing, chain_not_draining, chain_end_open | 4.89 at 14 | 84 |  | 45.96 |
| 2.284.0 | Sælatunga | 457.5 | 460.994 | 99.1 | swing | 460.99 at 0 | 30 |  | 41.06 |
| 2.303.0 | Dombås | 493.9 | 0.003 | 0.0 | swing, downstream_unread, chain_not_draining | 0.00 at -14 | 30 |  | 26.06 |
| 3.22.0 | Høgfoss | 299.8 | 330.409 | 99.2 | chain_not_draining, chain_end_open | 330.41 at 38 | 40 |  | 24.56 |
| 12.70.0 | Etna | 568.5 | 565.452 | 98.8 | swing, downstream_unread, chain_not_draining | 565.45 at 0 | 130 |  | 54.37 |
| 12.178.0 | Eggedal | 310.6 | 0.000 | 0.0 | swing, chain_not_draining | 309.97 at 28 | 32 | lake above P | 34.07 |
| 12.188.0 | Langtjernbekk | 4.7 | 0.001 | 0.0 | swing, chain_not_draining | 0.00 at 24 | 31 | lake above P | 4.71 |
| 19.80.0 | Stigvassåi | 14.5 | 15.197 | 98.9 | chain_not_draining | 0.00 at 0 | 30 |  | 15.20 |
| 19.96.0 | Storgama ovf. | 0.6 | 0.000 | 0.1 | swing, chain_not_draining | 0.00 at 28 | 30 | lake above P | 0.60 |
| 20.2.0 | Austenå | 277.2 | 290.883 | 99.6 | chain_not_draining | 290.88 at 0 | 30 |  | 36.19 |
| 20.11.0 | Tveitdalen | 0.4 | 0.003 | 0.7 | swing, chain_not_draining | 0.42 at 10 | 33 |  | 0.45 |
| 22.16.0 | Myglevatn ndf. | 182.2 | 182.405 | 99.4 | swing, chain_not_draining | 182.42 at 153 | 462 | lake line, lake above P | 33.91 |
| 22.22.0 | Søgne | 203.0 | 204.644 | 96.5 | chain_not_draining, chain_end_open | 204.71 at 310 | 454 |  | 16.36 |
| 38.1.0 | Holmen | 116.7 | 117.488 | 99.2 | swing, downstream_unread | 117.49 at 0 | 30 |  | 28.71 |
| 41.8.0 | Hellaugvatn | 27.5 | 0.124 | 0.4 | swing, chain_not_draining | 0.13 at 91 | 110 |  | 0.13 |
| 48.5.0 | Reinsnosvatn | 120.3 | 0.001 | 0.0 | swing, downstream_unread, chain_not_draining, chain_end_open | 0.01 at -10 | 129 |  | 0.35 |
| 62.15.0 | Kinne | 510.9 | 0.195 | 0.0 | swing, chain_not_draining | 0.19 at 0 | 30 |  | 56.66 |
| 73.21.0 | Frostdalen | 25.8 | 25.521 | 97.2 | chain_not_draining | 0.01 at -14 | 30 |  | 22.11 |
| 73.27.0 | Sula | 30.4 | 0.003 | 0.0 | swing, chain_not_draining | 0.00 at 14 | 30 |  | 19.82 |
| 79.3.0 | Nessedalselv | 30.2 | 0.006 | 0.0 | swing, chain_not_draining | 30.11 at 28 | 30 |  | 27.94 |
| 83.6.0 | Byttevatn | 104.5 | 0.000 | 0.0 | swing, chain_not_draining, chain_end_open | 0.06 at 24 | 31 |  | 5.71 |
| 83.12.0 | Haukedalsvatn ndf. | 205.3 | 206.980 | 99.5 | chain_not_draining | 0.03 at 0 | 30 | lake above P | 36.43 |
| 87.10.0 | Gloppenelv v/Bergheim | 218.5 | 0.092 | 0.0 | swing, downstream_unread, chain_not_draining | 0.08 at 0 | 93 |  | 0.09 |
| 88.11.0 | Strynsvatn | 485.1 | 0.677 | 0.1 | chain_not_draining, chain_end_open | 0.73 at 177 | 222 |  | 0.68 |
| 109.9.0 | Driva v/Risefoss | 744.4 | 0.000 | 0.0 | downstream_unread, chain_not_draining | 0.00 at 0 | 291 |  | 41.00 |
| 112.8.0 | Rinna | 87.9 | 88.228 | 99.2 | chain_not_draining | 88.23 at 0 | 30 |  | 31.15 |
| 122.11.0 | Eggafoss | 655.2 | 655.632 | 99.4 | chain_not_draining | 0.00 at -10 | 30 |  | 63.68 |
| 122.14.0 | Lillebudal bru | 168.1 | 168.951 | 99.3 | chain_not_draining | 0.00 at -14 | 30 |  | 49.85 |
| 122.17.0 | Hugdal bru | 545.9 | 0.013 | 0.0 | swing, chain_not_draining | 546.06 at 24 | 33 |  | 63.88 |
| 124.2.0 | Høggås bru | 494.5 | 643.919 | 99.4 | swing | 643.92 at -10 | 30 |  | 38.42 |
| 139.35.0 | Trangen | 852.3 | 0.001 | 0.0 | swing, downstream_unread, chain_not_draining | 0.00 at 14 | 87 |  | 19.26 |
| 150.1.0 | Sørra | 6.6 | 5.859 | 85.6 | chain_not_draining | 0.01 at 0 | 30 |  | 5.87 |
| 196.11.0 | Lille Rostavatn | 637.3 | 0.002 | 0.0 | swing, chain_not_draining | 0.31 at 119 | 186 |  | 0.06 |
| 200.4.0 | Skogsfjordvatn | 136.0 | 0.000 | 0.0 | swing, chain_not_draining, chain_end_open | 0.00 at 28 | 30 |  | 5.18 |
| 205.6.0 | Didnojokka | 111.0 | 111.171 | 98.4 | swing, downstream_unread, chain_not_draining | 111.17 at -14 | 186 |  | 38.23 |
| 212.10.0 | Masi | 5618.3 | 3.515 | 0.0 | swing, chain_not_draining | 3.72 at 318 | 337 |  | 3.52 |
| 234.18.0 | Polmak nye | 14171.0 | 6.279 | 0.0 | downstream_unread, chain_not_draining | 6.28 at 62 | 158 |  | 6.28 |
| 237.1.0 | Båtsfjord | 23.1 | 22.512 | 94.6 | swing | 22.51 at -20 | 30 |  | 22.52 |
| 311.6.0 | Nybergsund | 4418.1 | 0.001 | 0.0 | chain_not_draining, chain_end_open | 0.32 at -10 | 73 |  | 34.25 |
