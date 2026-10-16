| case | base total s | head total s | change | mesh same | features read | features clip | decode | refine |
|---|---:|---:|---:|---|---|---|---|---|
| tile | 0.417 | 0.416 | -0.3 % | True | - | - | 0.061 / 0.060 (-0.9 %) | 0.186 / 0.190 (+2.0 %) |
| tile t=1 | 0.704 | 0.699 | -0.7 % | True | - | - | 0.063 / 0.060 (-3.4 %) | 0.472 / 0.472 (+0.0 %) |
| quarter | 0.352 | 0.355 | +0.8 % | True | - | - | 0.044 / 0.044 (-1.5 %) | 0.171 / 0.172 (+0.7 %) |
| numedalslagen | 8.891 | 8.829 | -0.7 % | True | 0.194 / 0.197 (+1.7 %) | 6.092 / 6.000 (-1.5 %) | 0.921 / 0.921 (+0.0 %) | 0.660 / 0.659 (-0.1 %) |
| lagan | 19.296 | 14.542 | -24.6 % | True | 0.926 / 0.922 (-0.5 %) | 11.156 / 6.342 (-43.1 %) | 4.545 / 4.578 (+0.7 %) | 0.238 / 0.238 (+0.0 %) |
| geilo-al-ramp | 0.328 | 0.329 | +0.1 % | True | - | - | 0.139 / 0.138 (-0.2 %) | 0.100 / 0.099 (-0.3 %) |
| romsdal-slope | 3.582 | 3.578 | -0.1 % | True | - | - | 0.068 / 0.068 (+0.4 %) | 2.000 / 2.009 (+0.4 %) |

Section 7 targets (ratio head median / base median):
  lagan features clip: 11.156 -> 6.342 s, ratio 0.569, limit 0.65: met
  lagan total: 19.296 -> 14.542 s, ratio 0.754, limit 0.85: met
  numedalslagen features clip: 6.092 -> 6.000 s, ratio 0.985, limit 1.0: met

lagan base (3e6f87eb): total 19.296 s; phases at 10 % or more:
  features clip: 11.156 s (57.8 %)
  features clip: clean-up: 7.082 s (36.7 %)
  decode: 4.545 s (23.6 %)

lagan head (357808b1): total 14.542 s; phases at 10 % or more:
  features clip: 6.342 s (43.6 %)
  features clip: clean-up: 5.392 s (37.1 %)
  decode: 4.578 s (31.5 %)

numedalslagen base (3e6f87eb): total 8.891 s; phases at 10 % or more:
  features clip: 6.092 s (68.5 %)
  features clip: clean-up: 6.027 s (67.8 %)
  decode: 0.921 s (10.4 %)

numedalslagen head (357808b1): total 8.829 s; phases at 10 % or more:
  features clip: 6.000 s (68.0 %)
  features clip: clean-up: 5.937 s (67.2 %)
  decode: 0.921 s (10.4 %)
