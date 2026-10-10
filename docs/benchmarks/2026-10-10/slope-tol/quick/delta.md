| case | base total s | head total s | change | mesh same | refine | edge strip: scan (parallel) | edge strip: split + flip (serial) | check points: store | final check: scan (parallel) | final check: split + flip (serial) |
|---|---:|---:|---:|---|---|---|---|---|---|---|
| tile | 0.420 | 0.420 | -0.0 % | True | 0.187 / 0.186 (-0.7 %) | 0.001 / 0.001 (+7.0 %) | 0.001 / 0.001 (+9.5 %) | - | - | - |
| tile t=1 | 0.711 | 0.717 | +0.8 % | True | 0.478 / 0.480 (+0.2 %) | 0.001 / 0.001 (+8.4 %) | 0.001 / 0.001 (+10.7 %) | - | - | - |
| quarter | 0.355 | 0.357 | +0.5 % | True | 0.173 / 0.166 (-3.7 %) | 0.001 / 0.002 (+85.6 %) | 0.002 / 0.003 (+11.9 %) | - | - | - |
| numedalslagen | 9.902 | 9.970 | +0.7 % | True | 1.209 / 1.194 (-1.2 %) | 0.012 / 0.014 (+11.9 %) | 0.008 / 0.009 (+8.9 %) | - | - | - |
| lagan | 19.456 | 19.557 | +0.5 % | True | 0.663 / 0.650 (-2.1 %) | - | - | 0.276 / 0.283 (+2.8 %) | 0.301 / 0.302 (+0.4 %) | 0.033 / 0.035 (+4.4 %) |
| geilo-al-ramp | 0.334 | 0.338 | +1.2 % | True | 0.100 / 0.100 (+0.3 %) | 0.001 / 0.001 (+4.2 %) | 0.000 / 0.000 (-19.4 %) | - | - | - |

new case romsdal-slope: total 3.639 s (min 3.622, max 3.641), refine 2.032 s, slope 0.1058 s, max_error 9.9999 of tolerance 10.0
  total: 3.639 s (100 %)
  refine: 2.032 s (56 %)
  refine: split + flip (serial): 1.201 s (33 %)
  refine: scan (parallel): 0.565 s (16 %)
  refine: setup + output: 0.265 s (7 %)
  trim: 0.119 s (3 %)
  slope: 0.106 s (3 %)
  decode: 0.072 s (2 %)
