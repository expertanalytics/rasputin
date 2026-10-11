| case | base total s | head total s | change | mesh same | features read | features clip | decode | refine |
|---|---:|---:|---:|---|---|---|---|---|
| tile | 0.420 | 0.423 | +0.9 % | True | - | - | 0.065 / 0.062 (-4.9 %) | 0.186 / 0.191 (+2.7 %) |
| tile t=1 | 0.717 | 0.735 | +2.6 % | True | - | - | 0.067 / 0.064 (-4.9 %) | 0.480 / 0.497 (+3.7 %) |
| quarter | 0.357 | 0.362 | +1.3 % | True | - | - | 0.048 / 0.048 (-1.9 %) | 0.166 / 0.172 (+3.2 %) |
| numedalslagen | 9.970 | 8.775 | -12.0 % | False | 0.208 / 0.199 (-4.4 %) | 5.242 / 5.993 (+14.3 %) | 1.032 / 0.895 (-13.3 %) | 1.194 / 0.657 (-45.0 %) |
| lagan | 19.557 | 14.570 | -25.5 % | False | 0.963 / 0.924 (-4.1 %) | 9.696 / 6.345 (-34.6 %) | 4.598 / 4.574 (-0.5 %) | 0.650 / 0.236 (-63.6 %) |
| geilo-al-ramp | 0.338 | 0.330 | -2.5 % | True | - | - | 0.147 / 0.139 (-5.2 %) | 0.100 / 0.100 (-0.5 %) |
| romsdal-slope | 3.639 | 3.564 | -2.1 % | True | - | - | 0.072 / 0.066 (-8.0 %) | 2.032 / 2.005 (-1.3 %) |

Section 7 targets (ratio head median / base median):
  lagan features clip: 9.696 -> 6.345 s, ratio 0.654, limit 0.65: NOT MET
  lagan total: 19.557 -> 14.570 s, ratio 0.745, limit 0.85: met
  numedalslagen features clip: 5.242 -> 5.993 s, ratio 1.143, limit 1.0: NOT MET

lagan base (39e3ae34): total 19.557 s; phases at 10 % or more:
  features clip: 9.696 s (49.6 %)
  features clip: clean-up: 5.627 s (28.8 %)
  decode: 4.598 s (23.5 %)

lagan head (357808b1): total 14.570 s; phases at 10 % or more:
  features clip: 6.345 s (43.5 %)
  features clip: clean-up: 5.413 s (37.1 %)
  decode: 4.574 s (31.4 %)

numedalslagen base (39e3ae34): total 9.970 s; phases at 10 % or more:
  features clip: 5.242 s (52.6 %)
  features clip: clean-up: 5.177 s (51.9 %)
  refine: 1.194 s (12.0 %)
  decode: 1.032 s (10.4 %)

numedalslagen head (357808b1): total 8.775 s; phases at 10 % or more:
  features clip: 5.993 s (68.3 %)
  features clip: clean-up: 5.928 s (67.6 %)
  decode: 0.895 s (10.2 %)
