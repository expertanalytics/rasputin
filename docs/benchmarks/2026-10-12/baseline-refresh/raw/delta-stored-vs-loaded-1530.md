stored baseline commit 357808b1 vs run 2 (15:30, loaded machine) commit c317c7e4; medians of the timed runs, seconds

| case | phase | stored | run 2 | change | mesh same |
|---|---|---:|---:|---:|---|
| tile | total | 0.416 | 0.453 | +8.9 % | yes |
| tile | decode | 0.060 | 0.067 | +12.1 % |  |
| tile | start mesh: triangulate | 0.029 | 0.030 | +2.5 % |  |
| tile | refine | 0.190 | 0.201 | +5.9 % |  |
| tile | refine: scan (parallel) | 0.056 | 0.061 | +8.0 % |  |
| tile | refine: split + flip (serial) | 0.106 | 0.113 | +6.3 % |  |
| tile | refine: setup + output | 0.025 | 0.027 | +6.4 % |  |
| tile t=1 | total | 0.699 | 0.789 | +12.8 % | yes |
| tile t=1 | decode | 0.060 | 0.068 | +13.3 % |  |
| tile t=1 | refine | 0.472 | 0.545 | +15.5 % |  |
| tile t=1 | refine: scan (parallel) | 0.352 | 0.402 | +14.3 % |  |
| tile t=1 | refine: split + flip (serial) | 0.094 | 0.113 | +20.5 % |  |
| quarter | total | 0.355 | 0.388 | +9.3 % | yes |
| quarter | decode | 0.044 | 0.053 | +20.1 % |  |
| quarter | refine | 0.172 | 0.181 | +5.4 % |  |
| quarter | refine: scan (parallel) | 0.045 | 0.050 | +9.9 % |  |
| quarter | refine: split + flip (serial) | 0.107 | 0.112 | +4.1 % |  |
| numedalslagen | total | 8.829 | 9.343 | +5.8 % | yes |
| numedalslagen | decode | 0.921 | 0.986 | +7.1 % |  |
| numedalslagen | features clip | 6.000 | 6.342 | +5.7 % |  |
| numedalslagen | features clip: clean-up | 5.937 | 6.275 | +5.7 % |  |
| numedalslagen | refine | 0.659 | 0.719 | +9.1 % |  |
| lagan | total | 14.542 | 15.483 | +6.5 % | yes |
| lagan | decode | 4.578 | 4.783 | +4.5 % |  |
| lagan | features read | 0.922 | 1.094 | +18.7 % |  |
| lagan | features clip | 6.342 | 6.737 | +6.2 % |  |
| lagan | features clip: clean-up | 5.392 | 5.726 | +6.2 % |  |
| geilo-al-ramp | total | 0.329 | 0.364 | +10.9 % | yes |
| geilo-al-ramp | decode | 0.138 | 0.159 | +15.0 % |  |
| geilo-al-ramp | refine | 0.099 | 0.107 | +7.9 % |  |
| geilo-al-ramp | refine: scan (parallel) | 0.083 | 0.089 | +7.3 % |  |
| romsdal-slope | total | 3.578 | 3.876 | +8.3 % | yes |
| romsdal-slope | refine | 2.009 | 2.238 | +11.4 % |  |
| romsdal-slope | refine: scan (parallel) | 0.564 | 0.625 | +10.9 % |  |
| romsdal-slope | refine: split + flip (serial) | 1.184 | 1.337 | +12.9 % |  |
| romsdal-slope | refine: setup + output | 0.257 | 0.276 | +7.5 % |  |
