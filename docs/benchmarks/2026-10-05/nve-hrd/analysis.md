## Classes

| class | all | river-seeded | lake-seeded | no seed (refused before) |
|---|---:|---:|---:|---:|
| match | 74 | 33 | 41 | 0 |
| close | 5 | 2 | 3 | 0 |
| miss | 6 | 5 | 1 | 0 |
| uncertain | 39 | 39 | 0 | 0 |
| refused | 16 | 12 | 3 | 1 |
| total | 140 | | | |

`match_by`: {'overlap': 74, 'offset': 0}

## By size band (NVE's polygon area, km²)

| band | stations | match | close | miss | uncertain | refused | share uncertain |
|---|---:|---:|---:|---:|---:|---:|---:|
| under 10 | 12 | 5 | 3 | 0 | 4 | 0 | 33 % |
| 10-100 | 44 | 29 | 2 | 3 | 7 | 3 | 17 % |
| 100-1000 | 75 | 38 | 0 | 3 | 25 | 9 | 38 % |
| over 1000 | 9 | 2 | 0 | 0 | 3 | 4 | 60 % |

## By tile count (tiles our fine outline meets)

| tiles | stations | match | close | miss | uncertain | share uncertain |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 81 | 40 | 4 | 4 | 33 | 41 % |
| 2 | 26 | 22 | 0 | 1 | 3 | 12 % |
| 3-4 | 17 | 12 | 1 | 1 | 3 | 18 % |
| 5+ | 0 | 0 | 0 | 0 | 0 |  |

## Scored stations (match, close, miss): percentiles

| measure | min | p10 | p25 | p50 | p75 | p90 | max |
|---|---:|---:|---:|---:|---:|---:|---:|
| area_ratio | 0.005 | 0.989 | 0.996 | 1.001 | 1.009 | 1.023 | 1.730 |
| nve_in_ours | 0.004 | 0.967 | 0.982 | 0.988 | 0.992 | 0.994 | 0.997 |
| ours_in_nve | 0.571 | 0.953 | 0.976 | 0.987 | 0.991 | 0.993 | 0.996 |

Uncertain causes (a station can have several): {'swing': 25, 'downstream_unread': 9, 'chain_not_draining': 35, 'chain_end_open': 8, 'direction': 0}

Refusal causes: {'mixed_grid': 13, 'no_river': 1, 'other': 2}

## Refusals

| station | name | cause | NVE km² | message |
|---|---|---|---:|---|
| 2.142.0 | Knappom | other | 1643.0 | 642878 nodes the request needs are in no tile: x 400260 to 400870, y 6684070 to 6787750 (EPSG:25833) |
| 156.24.0 | Bogvatn | mixed_grid | 36.9 | the request selects tiles on two lattices, 7404_2_10m_z33.tif and 7304_1_10m_z33.tif: 7304_1_10m_z33.tif is 0.5 cell (5 m) east-west off 7404_2_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7404_2_10m_z33.tif, 7304_1_10m_z33.tif) |
| 156.15.0 | Forsbakk | mixed_grid | 56.2 | the request selects tiles on two lattices, 7304_4_10m_z33.tif and 7304_1_10m_z33.tif: 7304_1_10m_z33.tif is 0.5 cell (5 m) east-west off 7304_4_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7304_4_10m_z33.tif, 7304_1_10m_z33.tif) |
| 191.2.0 | Øvrevatn | mixed_grid | 526.8 | the request selects tiles on two lattices, 7505_1_10m_z33.tif and 7606_2_10m_z33.tif: 7606_2_10m_z33.tif is 0.5 cell (5 m) east-west off 7505_1_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7505_1_10m_z33.tif, 7606_2_10m_z33.tif) |
| 203.2.0 | Jægervatn | mixed_grid | 93.7 | the request selects tiles on two lattices, 7706_1_10m_z33.tif and 7707_3_10m_z33.tif: 7707_3_10m_z33.tif is 0.5 cell (5 m) east-west off 7706_1_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7706_1_10m_z33.tif, 7707_3_10m_z33.tif) |
| 206.3.0 | Manndalen bru | mixed_grid | 200.4 | the request selects tiles on two lattices, 7606_1_10m_z33.tif and 7707_3_10m_z33.tif: 7707_3_10m_z33.tif is 0.5 cell (5 m) east-west off 7606_1_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7606_1_10m_z33.tif, 7707_3_10m_z33.tif) |
| 208.2.0 | Oksfjordvatn | mixed_grid | 265.8 | the request selects tiles on two lattices, 7707_4_10m_z33.tif and 7707_1_10m_z33.tif: 7707_1_10m_z33.tif is 0.5 cell (5 m) east-west off 7707_4_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7707_4_10m_z33.tif, 7707_1_10m_z33.tif) |
| 208.3.0 | Svartfossberget | mixed_grid | 1932.4 | the request selects tiles on two lattices, 7707_2_10m_z33.tif and 7707_3_10m_z33.tif: 7707_3_10m_z33.tif is 0.5 cell (5 m) east-west off 7707_2_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7707_2_10m_z33.tif, 7707_3_10m_z33.tif) |
| 209.4.0 | Lillefossen | mixed_grid | 330.6 | the request selects tiles on two lattices, 7707_2_10m_z33.tif and 7707_1_10m_z33.tif: 7707_1_10m_z33.tif is 0.5 cell (5 m) east-west off 7707_2_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7707_2_10m_z33.tif, 7707_1_10m_z33.tif) |
| 212.48.0 | Sagafoss | mixed_grid | 234.2 | the request selects tiles on two lattices, 7707_2_10m_z33.tif and 7707_1_10m_z33.tif: 7707_1_10m_z33.tif is 0.5 cell (5 m) east-west off 7707_2_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7707_2_10m_z33.tif, 7707_1_10m_z33.tif) |
| 212.49.0 | Halsnes | mixed_grid | 145.2 | the request selects tiles on two lattices, 7708_4_10m_z33.tif and 7707_1_10m_z33.tif: 7707_1_10m_z33.tif is 0.5 cell (5 m) east-west off 7708_4_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7708_4_10m_z33.tif, 7707_1_10m_z33.tif) |
| 213.2.0 | Leirbotnvatn | mixed_grid | 135.2 | the request selects tiles on two lattices, 7708_4_10m_z33.tif and 7808_3_10m_z33.tif: 7808_3_10m_z33.tif is 0.5 cell (5 m) east-west off 7708_4_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7708_4_10m_z33.tif, 7808_3_10m_z33.tif) |
| 213.4.0 | Kvalsund | mixed_grid | 124.8 | the request selects tiles on two lattices, 7808_4_10m_z33.tif and 7808_3_10m_z33.tif: 7808_3_10m_z33.tif is 0.5 cell (5 m) east-west off 7808_4_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7808_4_10m_z33.tif, 7808_3_10m_z33.tif) |
| 223.2.0 | Lombola | mixed_grid | 873.9 | the request selects tiles on two lattices, 7708_1_10m_z33.tif and 7808_3_10m_z33.tif: 7808_3_10m_z33.tif is 0.5 cell (5 m) east-west off 7708_1_10m_z33.tif's lattice; tiles are not resampled onto one another, but a --bbox inside one lattice is meshed (tiles: 7708_1_10m_z33.tif, 7808_3_10m_z33.tif) |
| 244.2.0 | Neiden | other | 2944.6 | 1714266 nodes the request needs are in no tile: x 987960 to 999740, y 7735210 to 7749740 (EPSG:25833) |
| 311.4.0 | Femundsenden (Femunden) | no_river | 1788.7 | no mapped river line within 500 m of the station |

## Every station

| station | name | class | by | seed | NVE km² | ours km² | NVE's in ours % | ours in NVE's % | ratio | offset m | causes | tiles | s | peak GB |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| 2.11.0 | Narsjø | match | overlap | lake | 119.4 | 120.26 | 98.2 | 97.5 | 1.007 | 75.7 |  | 1 | 7.9 | 1.12 |
| 2.32.0 | Atnasjø | match | overlap | lake | 463.3 | 460.72 | 98.5 | 99.1 | 0.994 | 87.9 |  | 1 | 11.2 | 1.62 |
| 2.142.0 | Knappom | refused |  | river | 1643.0 |  |  |  |  |  |  |  | 7.1 | 1.95 |
| 2.265.0 | Unsetåa | match | overlap | river | 620.2 | 623.40 | 99.4 | 98.9 | 1.005 | 75.8 |  | 4 | 95.9 | 4.37 |
| 2.279.0 | Kråkfoss | uncertain |  | river | 434.5 | 0.00 | 0.0 | 0.0 | 0.000 | 3057.7 | swing, chain_not_draining, chain_end_open | 1 | 1.2 | 4.38 |
| 2.284.0 | Sælatunga | uncertain |  | river | 457.5 | 460.99 | 99.1 | 98.4 | 1.008 | 99.4 | swing | 1 | 88.2 | 4.38 |
| 2.290.0 | Brustuen | match | overlap | river | 254.0 | 252.76 | 98.9 | 99.4 | 0.995 | 40.8 |  | 2 | 21.4 | 4.38 |
| 2.303.0 | Dombås | uncertain |  | river | 493.9 | 0.00 | 0.0 | 0.0 | 0.000 | 4123.0 | swing, downstream_unread, chain_not_draining | 1 | 1.3 | 2.69 |
| 2.439.0 | Kvarstadseter | match | overlap | river | 375.1 | 374.38 | 99.0 | 99.2 | 0.998 | 58.7 |  | 2 | 21.6 | 2.79 |
| 2.633.0 | Stortorp | match | overlap | river | 87.2 | 87.02 | 97.7 | 97.9 | 0.998 | 50.4 |  | 1 | 6.3 | 1.49 |
| 3.22.0 | Høgfoss | uncertain |  | river | 299.8 | 330.41 | 99.2 | 90.0 | 1.102 | 239.3 | chain_not_draining, chain_end_open | 1 | 94.0 | 3.30 |
| 6.10.0 | Gryta | close |  | river | 7.1 | 7.65 | 97.5 | 90.0 | 1.084 | 61.0 |  | 1 | 1.2 | 3.30 |
| 11.4.0 | Elgtjern | match | overlap | lake | 7.1 | 6.98 | 96.1 | 97.5 | 0.985 | 31.9 |  | 2 | 1.4 | 3.35 |
| 12.70.0 | Etna | uncertain |  | river | 568.5 | 565.45 | 98.8 | 99.3 | 0.995 | 64.5 | swing, downstream_unread, chain_not_draining | 4 | 51.0 | 4.16 |
| 12.171.0 | Hølervatn | match | overlap | lake | 79.4 | 81.26 | 99.1 | 96.8 | 1.023 | 66.7 |  | 1 | 2.5 | 4.17 |
| 12.178.0 | Eggedal | uncertain |  | river | 310.6 | 0.00 | 0.0 | 100.0 | 0.000 | 2520.1 | swing, chain_not_draining | 1 | 15.6 | 4.33 |
| 12.188.0 | Langtjernbekk | uncertain |  | river | 4.7 | 0.00 | 0.0 | 0.0 | 0.000 | 394.2 | swing, chain_not_draining | 1 | 0.5 | 4.33 |
| 12.192.0 | Sundbyfoss | match | overlap | river | 74.8 | 75.59 | 99.0 | 98.0 | 1.011 | 34.8 |  | 2 | 5.5 | 4.42 |
| 12.193.0 | Fiskum | match | overlap | river | 51.5 | 51.68 | 98.8 | 98.5 | 1.003 | 35.4 |  | 2 | 5.0 | 4.52 |
| 12.197.0 | Grunke | match | overlap | river | 184.7 | 184.03 | 98.9 | 99.3 | 0.996 | 47.8 |  | 2 | 19.4 | 4.52 |
| 12.207.0 | Vinde-elv | match | overlap | river | 269.9 | 273.13 | 99.4 | 98.2 | 1.012 | 59.6 |  | 2 | 24.9 | 5.01 |
| 12.215.0 | Storeskar | miss |  | river | 119.7 | 158.07 | 98.7 | 74.7 | 1.321 | 695.9 |  | 1 | 18.5 | 2.34 |
| 15.21.0 | Jondalselv | match | overlap | river | 126.8 | 127.51 | 99.2 | 98.6 | 1.006 | 40.7 |  | 1 | 18.7 | 1.50 |
| 15.49.0 | Halledalsvatn | miss |  | river | 59.4 | 102.71 | 98.8 | 57.1 | 1.730 | 1066.0 |  | 2 | 18.0 | 1.50 |
| 16.66.0 | Grosettjern | close |  | lake | 6.5 | 6.21 | 93.6 | 98.7 | 0.948 | 35.4 |  | 1 | 0.8 | 1.25 |
| 16.122.0 | Grovåi | match | overlap | river | 42.2 | 42.69 | 98.6 | 97.4 | 1.012 | 35.4 |  | 1 | 4.4 | 1.27 |
| 16.127.0 | Viertjern | match | overlap | river | 46.7 | 45.63 | 96.5 | 98.8 | 0.977 | 46.5 |  | 1 | 4.9 | 1.29 |
| 16.194.0 | Kilen | match | overlap | river | 118.1 | 118.92 | 99.0 | 98.3 | 1.007 | 51.2 |  | 3 | 19.7 | 1.59 |
| 18.10.0 | Gjerstad | match | overlap | river | 236.2 | 235.32 | 99.0 | 99.3 | 0.996 | 36.6 |  | 4 | 24.4 | 1.67 |
| 18.11.0 | Tjellingtjernbekk | match | overlap | river | 2.0 | 1.96 | 96.5 | 96.3 | 1.004 | 20.1 |  | 2 | 1.5 | 1.26 |
| 19.79.0 | Gravå | close |  | river | 6.1 | 5.33 | 82.4 | 95.0 | 0.867 | 100.4 |  | 1 | 1.3 | 1.27 |
| 19.80.0 | Stigvassåi | uncertain |  | river | 14.5 | 15.20 | 98.9 | 94.3 | 1.049 | 46.2 | chain_not_draining | 1 | 3.7 | 1.29 |
| 19.82.0 | Rauåna | match | overlap | river | 8.9 | 8.87 | 97.2 | 97.6 | 0.996 | 23.7 |  | 1 | 3.6 | 1.29 |
| 19.96.0 | Storgama ovf. | uncertain |  | river | 0.6 | 0.00 | 0.1 | 75.0 | 0.001 | 143.4 | swing, chain_not_draining | 1 | 0.5 | 1.29 |
| 19.104.0 | Songedalsåi | match | overlap | river | 65.6 | 65.44 | 99.0 | 99.3 | 0.997 | 25.0 |  | 3 | 5.5 | 1.36 |
| 20.2.0 | Austenå | uncertain |  | river | 277.2 | 290.88 | 99.6 | 94.9 | 1.049 | 102.0 | chain_not_draining | 4 | 96.5 | 3.44 |
| 20.11.0 | Tveitdalen | uncertain |  | river | 0.4 | 0.00 | 0.7 | 96.9 | 0.007 | 167.9 | swing, chain_not_draining | 1 | 0.7 | 3.48 |
| 22.16.0 | Myglevatn ndf. | uncertain |  | river | 182.2 | 182.41 | 99.4 | 99.3 | 1.001 | 27.3 | swing, chain_not_draining | 1 | 18.3 | 3.74 |
| 22.22.0 | Søgne | uncertain |  | river | 203.0 | 204.64 | 96.5 | 95.8 | 1.008 | 102.1 | chain_not_draining, chain_end_open | 1 | 67.3 | 3.75 |
| 24.8.0 | Møska (Skolandsvatnet) | match | overlap | lake | 120.9 | 121.24 | 99.1 | 98.9 | 1.003 | 31.9 |  | 1 | 6.4 | 3.13 |
| 24.9.0 | Tingvatn (Lygne) | match | overlap | lake | 272.2 | 271.76 | 99.3 | 99.5 | 0.998 | 24.5 |  | 3 | 10.0 | 3.33 |
| 26.29.0 | Refsvatn | close |  | lake | 53.2 | 62.79 | 99.3 | 84.0 | 1.181 | 201.0 |  | 4 | 2.3 | 3.33 |
| 27.15.0 | Austrumdal (Austrumdalsvatnet) | match | overlap | lake | 60.8 | 60.80 | 98.8 | 98.8 | 0.999 | 28.0 |  | 2 | 2.8 | 2.18 |
| 35.9.0 | Osali (Botnavatnet) | close |  | lake | 22.5 | 23.78 | 98.7 | 93.3 | 1.058 | 67.6 |  | 1 | 1.5 | 2.18 |
| 35.16.0 | Djupadalsvatn | match | overlap | lake | 45.3 | 46.20 | 98.6 | 96.8 | 1.019 | 47.4 |  | 2 | 5.2 | 2.18 |
| 36.13.0 | Grimsvatn | match | overlap | lake | 34.4 | 33.56 | 96.5 | 98.8 | 0.977 | 46.0 |  | 1 | 1.7 | 1.32 |
| 38.1.0 | Holmen | uncertain |  | river | 116.7 | 117.49 | 99.2 | 98.5 | 1.007 | 37.3 | swing, downstream_unread | 2 | 9.9 | 1.50 |
| 41.8.0 | Hellaugvatn | uncertain |  | river | 27.5 | 0.12 | 0.4 | 99.8 | 0.005 | 916.8 | swing, chain_not_draining | 1 | 1.0 | 1.26 |
| 42.2.0 | Djupevad | match | overlap | river | 31.0 | 30.98 | 98.8 | 98.9 | 1.000 | 22.2 |  | 1 | 3.8 | 1.26 |
| 48.1.0 | Sandvenvatn | match | overlap | lake | 469.6 | 467.72 | 98.8 | 99.2 | 0.996 | 63.9 |  | 1 | 11.8 | 1.73 |
| 48.5.0 | Reinsnosvatn | uncertain |  | river | 120.3 | 0.00 | 0.0 | 100.0 | 0.000 | 1897.9 | swing, downstream_unread, chain_not_draining, chain_end_open | 1 | 0.9 | 1.74 |
| 50.1.0 | Hølen | match | overlap | lake | 231.4 | 232.04 | 98.8 | 98.5 | 1.003 | 54.6 |  | 3 | 7.3 | 1.74 |
| 55.4.0 | Røykenes | match | overlap | river | 50.2 | 50.34 | 99.2 | 98.9 | 1.004 | 18.6 |  | 1 | 15.4 | 1.24 |
| 62.18.0 | Svartavatn | match | overlap | lake | 72.4 | 72.31 | 99.1 | 99.2 | 0.999 | 24.6 |  | 1 | 6.1 | 1.28 |
| 62.5.0 | Bulken (Vangsvatnet) | match | overlap | lake | 1091.3 | 1093.63 | 99.6 | 99.4 | 1.002 | 37.0 |  | 4 | 50.5 | 3.47 |
| 62.10.0 | Myrkdalsvatn | match | overlap | lake | 157.6 | 157.82 | 99.2 | 99.0 | 1.001 | 34.7 |  | 1 | 8.3 | 3.76 |
| 62.14.0 | Slondalsvatn | match | overlap | lake | 41.8 | 41.45 | 97.9 | 98.8 | 0.991 | 47.4 |  | 2 | 1.9 | 3.80 |
| 62.15.0 | Kinne | uncertain |  | river | 510.9 | 0.19 | 0.0 | 70.9 | 0.000 | 3229.5 | swing, chain_not_draining | 1 | 79.0 | 4.27 |
| 73.21.0 | Frostdalen | uncertain |  | river | 25.8 | 25.52 | 97.2 | 98.2 | 0.990 | 49.0 | chain_not_draining | 1 | 4.3 | 4.28 |
| 73.27.0 | Sula | uncertain |  | river | 30.4 | 0.00 | 0.0 | 0.0 | 0.000 | 886.0 | swing, chain_not_draining | 1 | 3.4 | 4.28 |
| 76.5.0 | Nigardsbrevatn | match | overlap | lake | 65.2 | 66.82 | 98.1 | 95.8 | 1.024 | 107.6 |  | 1 | 1.9 | 4.28 |
| 77.3.0 | Sogndalsvatn | match | overlap | lake | 111.5 | 110.81 | 98.5 | 99.1 | 0.994 | 47.5 |  | 1 | 2.6 | 4.28 |
| 79.3.0 | Nessedalselv | uncertain |  | river | 30.2 | 0.01 | 0.0 | 96.9 | 0.000 | 1158.4 | swing, chain_not_draining | 1 | 2.9 | 4.28 |
| 81.1.0 | Hersvikvatn (Hagevatnet) | match | overlap | lake | 7.1 | 7.08 | 96.9 | 96.5 | 1.003 | 22.6 |  | 1 | 0.7 | 4.28 |
| 82.4.0 | Nautsundvatn | match | overlap | lake | 218.8 | 218.46 | 99.3 | 99.4 | 0.998 | 25.6 |  | 2 | 8.5 | 4.28 |
| 83.2.0 | Viksvatn (Hestadfjorden) | miss |  | lake | 507.9 | 19.35 | 3.8 | 99.3 | 0.038 | 2960.3 |  | 1 | 1.5 | 3.16 |
| 83.6.0 | Byttevatn | uncertain |  | river | 104.5 | 0.00 | 0.0 | 100.0 | 0.000 | 1609.9 | swing, chain_not_draining, chain_end_open | 1 | 0.8 | 3.16 |
| 83.7.0 | Grønengstølsvatn | match | overlap | lake | 65.5 | 65.62 | 98.9 | 98.7 | 1.002 | 38.1 |  | 2 | 2.3 | 3.16 |
| 83.12.0 | Haukedalsvatn ndf. | uncertain |  | river | 205.3 | 206.98 | 99.5 | 98.7 | 1.008 | 44.7 | chain_not_draining | 2 | 19.8 | 2.13 |
| 84.20.0 | Holsenvatn | match | overlap | lake | 71.2 | 70.52 | 98.3 | 99.2 | 0.991 | 35.0 |  | 1 | 2.6 | 1.35 |
| 85.4.0 | Straumstad (Solheimsvatnet) | match | overlap | lake | 109.7 | 113.88 | 98.7 | 95.1 | 1.038 | 108.3 |  | 2 | 2.7 | 1.35 |
| 86.10.0 | Åvatn (Ommedalsvatnet) | match | overlap | lake | 162.1 | 162.18 | 99.2 | 99.1 | 1.000 | 33.8 |  | 1 | 7.3 | 1.35 |
| 86.12.0 | Skjerdalselv | match | overlap | river | 23.7 | 23.93 | 98.7 | 97.6 | 1.011 | 33.0 |  | 1 | 3.9 | 1.18 |
| 87.10.0 | Gloppenelv v/Bergheim | uncertain |  | river | 218.5 | 0.09 | 0.0 | 82.7 | 0.000 | 2467.5 | swing, downstream_unread, chain_not_draining | 1 | 1.2 | 1.18 |
| 88.4.0 | Lovatn | match | overlap | lake | 234.9 | 236.00 | 99.0 | 98.5 | 1.005 | 80.6 |  | 1 | 4.0 | 1.24 |
| 88.11.0 | Strynsvatn | uncertain |  | river | 485.1 | 0.68 | 0.1 | 100.0 | 0.001 | 3543.0 | chain_not_draining, chain_end_open | 1 | 1.7 | 1.24 |
| 88.30.0 | Nordre Oldevatn | match | overlap | lake | 201.7 | 199.92 | 97.9 | 98.7 | 0.991 | 93.4 |  | 1 | 7.7 | 1.50 |
| 97.1.0 | Fetvatn (Fitjavatnet) | match | overlap | lake | 89.0 | 89.45 | 99.4 | 98.9 | 1.005 | 27.6 |  | 1 | 2.2 | 1.50 |
| 98.4.0 | Øye ndf. | match | overlap | river | 139.2 | 140.84 | 99.5 | 98.4 | 1.012 | 47.7 |  | 2 | 20.4 | 1.28 |
| 101.1.0 | Engsetvatn | match | overlap | river | 39.9 | 40.64 | 98.9 | 97.1 | 1.018 | 54.2 |  | 1 | 3.2 | 1.32 |
| 103.1.0 | Ulvåa v/Storhølen | match | overlap | river | 435.3 | 439.61 | 99.5 | 98.5 | 1.010 | 71.1 |  | 1 | 26.8 | 1.76 |
| 104.23.0 | Vistdal | match | overlap | river | 66.5 | 66.19 | 99.1 | 99.6 | 0.995 | 21.5 |  | 1 | 4.6 | 1.70 |
| 105.1.0 | Osenelv v/Øren | miss |  | river | 138.1 | 0.75 | 0.4 | 77.6 | 0.005 | 2051.5 |  | 1 | 1.1 | 1.70 |
| 109.9.0 | Driva v/Risefoss | uncertain |  | river | 744.4 | 0.00 | 0.0 | 100.0 | 0.000 | 5192.5 | downstream_unread, chain_not_draining | 1 | 1.2 | 1.70 |
| 109.21.0 | Driva v/Svoni | match | overlap | river | 136.0 | 137.54 | 97.2 | 96.1 | 1.012 | 149.1 |  | 1 | 19.5 | 1.70 |
| 112.8.0 | Rinna | uncertain |  | river | 87.9 | 88.23 | 99.2 | 98.8 | 1.004 | 36.4 | chain_not_draining | 1 | 14.9 | 1.36 |
| 121.20.0 | Åmot | match | overlap | river | 282.7 | 283.77 | 99.7 | 99.3 | 1.004 | 24.5 |  | 4 | 22.9 | 1.36 |
| 122.11.0 | Eggafoss | uncertain |  | river | 655.2 | 655.63 | 99.4 | 99.3 | 1.001 | 50.0 | chain_not_draining | 2 | 107.7 | 4.69 |
| 122.14.0 | Lillebudal bru | uncertain |  | river | 168.1 | 168.95 | 99.3 | 98.8 | 1.005 | 44.2 | chain_not_draining | 1 | 19.2 | 4.69 |
| 122.17.0 | Hugdal bru | uncertain |  | river | 545.9 | 0.01 | 0.0 | 64.8 | 0.000 | 3708.7 | swing, chain_not_draining | 1 | 71.8 | 3.87 |
| 124.2.0 | Høggås bru | uncertain |  | river | 494.5 | 643.92 | 99.4 | 76.3 | 1.302 | 927.0 | swing | 4 | 83.4 | 3.73 |
| 127.11.0 | Veravatn | match | overlap | lake | 175.5 | 175.37 | 98.2 | 98.3 | 0.999 | 75.8 |  | 1 | 9.2 | 3.58 |
| 128.9.0 | Leksdalsvatn | match | overlap | lake | 178.4 | 177.72 | 98.9 | 99.3 | 0.996 | 41.6 |  | 1 | 8.0 | 3.58 |
| 133.7.0 | Krinsvatn (Kringsvatnet) | match | overlap | lake | 205.7 | 205.42 | 99.2 | 99.4 | 0.999 | 31.6 |  | 1 | 7.5 | 2.10 |
| 138.1.0 | Øyungen | match | overlap | lake | 238.9 | 239.25 | 99.4 | 99.2 | 1.001 | 30.4 |  | 2 | 10.1 | 2.10 |
| 139.35.0 | Trangen | uncertain |  | river | 852.3 | 0.00 | 0.0 | 100.0 | 0.000 | 3547.8 | swing, downstream_unread, chain_not_draining | 1 | 2.5 | 1.77 |
| 140.2.0 | Salsvatn | match | overlap | lake | 432.1 | 432.26 | 99.6 | 99.6 | 1.000 | 23.6 |  | 2 | 16.5 | 1.57 |
| 148.2.0 | Mevatnet | match | overlap | river | 108.9 | 109.83 | 99.7 | 98.9 | 1.009 | 23.3 |  | 2 | 17.4 | 1.64 |
| 150.1.0 | Sørra | uncertain |  | river | 6.6 | 5.86 | 85.6 | 96.1 | 0.890 | 100.7 | chain_not_draining | 1 | 1.0 | 1.65 |
| 151.13.0 | Øvre Glugvatn | match | overlap | lake | 60.7 | 60.54 | 98.5 | 98.7 | 0.998 | 45.9 |  | 1 | 2.0 | 1.65 |
| 151.15.0 | Nervoll | match | overlap | river | 655.0 | 653.68 | 98.4 | 98.6 | 0.998 | 136.5 |  | 4 | 28.3 | 1.87 |
| 152.4.0 | Fustvatn | match | overlap | lake | 525.8 | 527.60 | 99.5 | 99.1 | 1.003 | 53.5 |  | 2 | 13.8 | 2.02 |
| 153.1.0 | Storvatn | miss |  | river | 49.1 | 2.70 | 5.4 | 98.0 | 0.055 | 1194.5 |  | 1 | 1.3 | 2.02 |
| 156.24.0 | Bogvatn | refused |  | river | 36.9 |  |  |  |  |  |  |  | 0.6 | 2.02 |
| 156.15.0 | Forsbakk | refused |  | river | 56.2 |  |  |  |  |  |  |  | 0.4 | 2.02 |
| 168.3.0 | Lakså bru | match | overlap | river | 26.8 | 26.79 | 99.0 | 99.0 | 0.999 | 20.2 |  | 1 | 3.8 | 2.04 |
| 172.8.0 | Rauvatn | match | overlap | lake | 19.9 | 19.99 | 98.8 | 98.5 | 1.003 | 20.6 |  | 1 | 1.5 | 2.04 |
| 177.4.0 | Sneisvatn | match | overlap | lake | 29.3 | 29.12 | 98.9 | 99.4 | 0.994 | 15.3 |  | 1 | 1.7 | 1.54 |
| 178.1.0 | Langvatn | match | overlap | lake | 18.7 | 18.95 | 99.3 | 97.8 | 1.015 | 23.7 |  | 1 | 1.6 | 1.54 |
| 185.1.0 | Gåslandsvatn | match | overlap | lake | 7.7 | 7.76 | 97.8 | 97.6 | 1.002 | 23.7 |  | 1 | 0.6 | 1.54 |
| 186.2.0 | Ånesvatn | match | overlap | lake | 46.9 | 46.86 | 97.7 | 97.7 | 1.000 | 58.5 |  | 1 | 1.8 | 1.59 |
| 189.3.0 | Tennevikvatn | match | overlap | river | 85.3 | 85.51 | 99.1 | 98.8 | 1.003 | 31.9 |  | 1 | 4.9 | 1.59 |
| 191.2.0 | Øvrevatn | refused |  | lake | 526.8 |  |  |  |  |  |  |  | 8.0 | 1.51 |
| 196.7.0 | Ytre Fiskeløsvatn | match | overlap | lake | 54.5 | 54.26 | 98.4 | 98.9 | 0.995 | 39.8 |  | 2 | 2.3 | 1.51 |
| 196.11.0 | Lille Rostavatn | uncertain |  | river | 637.3 | 0.00 | 0.0 | 100.0 | 0.000 | 3765.8 | swing, chain_not_draining | 1 | 2.0 | 0.96 |
| 200.4.0 | Skogsfjordvatn | uncertain |  | river | 136.0 | 0.00 | 0.0 | 100.0 | 0.000 | 1979.9 | swing, chain_not_draining, chain_end_open | 1 | 0.7 | 1.02 |
| 203.2.0 | Jægervatn | refused |  | lake | 93.7 |  |  |  |  |  |  |  | 0.8 | 1.02 |
| 205.6.0 | Didnojokka | uncertain |  | river | 111.0 | 111.17 | 98.4 | 98.3 | 1.001 | 66.8 | swing, downstream_unread, chain_not_draining | 1 | 3.2 | 1.02 |
| 206.3.0 | Manndalen bru | refused |  | river | 200.4 |  |  |  |  |  |  |  | 1.4 | 1.11 |
| 208.2.0 | Oksfjordvatn | refused |  | river | 265.8 |  |  |  |  |  |  |  | 0.3 | 1.11 |
| 208.3.0 | Svartfossberget | refused |  | river | 1932.4 |  |  |  |  |  |  |  | 0.3 | 1.11 |
| 209.4.0 | Lillefossen | refused |  | river | 330.6 |  |  |  |  |  |  |  | 0.6 | 1.11 |
| 212.10.0 | Masi | uncertain |  | river | 5618.3 | 3.51 | 0.0 | 4.0 | 0.001 | 10205.7 | swing, chain_not_draining | 1 | 7.8 | 1.11 |
| 212.48.0 | Sagafoss | refused |  | river | 234.2 |  |  |  |  |  |  |  | 5.3 | 1.11 |
| 212.49.0 | Halsnes | refused |  | river | 145.2 |  |  |  |  |  |  |  | 0.3 | 0.91 |
| 213.2.0 | Leirbotnvatn | refused |  | lake | 135.2 |  |  |  |  |  |  |  | 0.3 | 0.91 |
| 213.4.0 | Kvalsund | refused |  | river | 124.8 |  |  |  |  |  |  |  | 0.5 | 1.00 |
| 223.2.0 | Lombola | refused |  | river | 873.9 |  |  |  |  |  |  |  | 1.6 | 1.00 |
| 234.13.0 | Veahkkava, Iesjokka | match | overlap | river | 2081.0 | 2077.35 | 99.2 | 99.4 | 0.998 | 104.9 |  | 4 | 144.6 | 5.54 |
| 234.18.0 | Polmak nye | uncertain |  | river | 14171.0 | 6.28 | 0.0 | 97.1 | 0.000 | 13162.9 | downstream_unread, chain_not_draining | 1 | 62.5 | 6.66 |
| 237.1.0 | Båtsfjord | uncertain |  | river | 23.1 | 22.51 | 94.6 | 97.1 | 0.974 | 81.8 | swing | 1 | 4.6 | 6.69 |
| 244.2.0 | Neiden | refused |  | river | 2944.6 |  |  |  |  |  |  |  | 6.9 | 6.81 |
| 246.9.0 | Sametielv | match | overlap | river | 254.7 | 251.47 | 97.8 | 99.0 | 0.987 | 55.8 |  | 1 | 24.2 | 7.12 |
| 247.3.0 | Karpelva | match | overlap | river | 125.3 | 128.17 | 97.7 | 95.5 | 1.023 | 122.8 |  | 4 | 19.4 | 7.12 |
| 307.5.0 | Murusjø | match | overlap | lake | 346.2 | 352.24 | 98.8 | 97.1 | 1.017 | 104.4 |  | 4 | 11.9 | 7.20 |
| 307.7.0 | Landbru | miss |  | river | 61.5 | 78.39 | 98.5 | 77.3 | 1.275 | 407.2 |  | 4 | 5.5 | 7.20 |
| 308.1.0 | Lenglingen | match | overlap | lake | 452.5 | 452.16 | 99.2 | 99.3 | 0.999 | 49.9 |  | 2 | 13.7 | 7.21 |
| 311.4.0 | Femundsenden (Femunden) | refused |  |  | 1788.7 |  |  |  |  |  |  |  | 0.1 | 3.47 |
| 311.6.0 | Nybergsund | uncertain |  | river | 4418.1 | 0.00 | 0.0 | 0.0 | 0.000 | 6746.0 | chain_not_draining, chain_end_open | 1 | 8.6 | 3.47 |
| 313.10.0 | Magnor | match | overlap | river | 357.9 | 356.75 | 98.8 | 99.1 | 0.997 | 54.9 |  | 2 | 25.6 | 1.49 |

## Time and memory

Stations' own seconds: total 1993 s, median 4.9 s, max 145 s (234.13.0).
Sampled peak RSS per station: median 1.82 GB, max 7.21 GB (308.1.0).
Whole process (`/usr/bin/time -l`): 7745503232  maximum resident set size

## Step 6 checks and increment 22's outline guarantees

Stations checked: 124 not refused (79 river-seeded, 45 lake-seeded).

| check | failures | stations |
|---|---:|---|
| placed node's count = catchment nodes (oracle; a failure is a defect) | 0 |  |
| counts strictly rising downstream along the chain, end closed (`monotone`) | 32 | 2.279.0, 3.22.0, 12.70.0, 12.178.0, 12.188.0, 19.80.0, 19.96.0, 20.2.0, 20.11.0, 22.16.0, 22.22.0, 41.8.0, 48.5.0, 73.21.0, 73.27.0, 79.3.0, 83.6.0, 83.12.0, 87.10.0, 88.11.0, 112.8.0, 122.11.0, 122.14.0, 122.17.0, 139.35.0, 150.1.0, 196.11.0, 200.4.0, 205.6.0, 212.10.0, 234.18.0, 311.6.0 |
| each chain node drains into the next (`drains`) | 35 | 2.279.0, 2.303.0, 3.22.0, 12.70.0, 12.178.0, 12.188.0, 19.80.0, 19.96.0, 20.2.0, 20.11.0, 22.16.0, 22.22.0, 41.8.0, 48.5.0, 62.15.0, 73.21.0, 73.27.0, 79.3.0, 83.6.0, 83.12.0, 87.10.0, 88.11.0, 109.9.0, 112.8.0, 122.11.0, 122.14.0, 122.17.0, 139.35.0, 150.1.0, 196.11.0, 200.4.0, 205.6.0, 212.10.0, 234.18.0, 311.6.0 |
| burn lowered no node off the chain (first window) | 0 |  |
| reduced outline simple | 0 |  |
| seed strictly inside the reduced outline | 0 |  |
| area kept (reduced − fine) | 0 over 1e-9 of the area | largest difference 2.62e-06 m² |

Stations whose burn did not hold (monotone or drains failed, or a node off the chain lowered): 35; of the 39 `uncertain`, they explain 35.

River rows with a lake on the reach above P: 13

| station | name | class | NVE's in ours % | ours in NVE's % | causes |
|---|---|---|---:|---:|---|
| 12.178.0 | Eggedal | uncertain | 0.0 | 100.0 | swing, chain_not_draining |
| 12.188.0 | Langtjernbekk | uncertain | 0.0 | 0.0 | swing, chain_not_draining |
| 12.197.0 | Grunke | match | 98.9 | 99.3 |  |
| 15.49.0 | Halledalsvatn | miss | 98.8 | 57.1 |  |
| 16.127.0 | Viertjern | match | 96.5 | 98.8 |  |
| 18.11.0 | Tjellingtjernbekk | match | 96.5 | 96.3 |  |
| 19.96.0 | Storgama ovf. | uncertain | 0.1 | 75.0 | swing, chain_not_draining |
| 22.16.0 | Myglevatn ndf. | uncertain | 99.4 | 99.3 | swing, chain_not_draining |
| 83.12.0 | Haukedalsvatn ndf. | uncertain | 99.5 | 98.7 | chain_not_draining |
| 101.1.0 | Engsetvatn | match | 98.9 | 97.1 |  |
| 148.2.0 | Mevatnet | match | 99.7 | 98.9 |  |
| 168.3.0 | Lakså bru | match | 99.0 | 99.0 |  |
| 189.3.0 | Tennevikvatn | match | 99.1 | 98.8 |  |
