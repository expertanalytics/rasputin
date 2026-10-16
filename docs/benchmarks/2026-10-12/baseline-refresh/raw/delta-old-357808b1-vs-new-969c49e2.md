# Stored baseline `357808b1` (2026-10-11) vs the refreshed baseline `969c49e2` (run 4, 2026-10-12 20:29, idle, AC)

Every phase at 5 % or more of its case's total, as `bench_quick.verdict`
selects them; band = max(5 %, the old record's (max - min) / median). The
tool itself printed `NO BASELINE: input $RASPUTIN_DATA/sweden_glo30_cache
changed` (the Lagan GLO-30 cache gained tiles), so this table was made by
calling `verdict` on the two records with the old inputs substituted into
the new one (`.venv/bin/python`, `tools/bench_quick.py` at `969c49e2`); the
judge lines it gives are:

    SLOWER: lagan features read 0.92 -> 1.03 s (+12.1 %, band 5.5 %)

No `MESH CHANGED`, no `BROKEN`; every case's mesh SHA-256 and `max_error`
equal the stored ones.

| line | old s | new s | change | band |
|---|---:|---:|---:|---:|
| tile total | 0.416 | 0.430 | +3.3 % | 5.0 % |
| tile decode | 0.060 | 0.070 | +15.8 % | 5.5 % |
| tile start mesh: triangulate | 0.029 | 0.029 | -0.1 % | 5.0 % |
| tile refine | 0.190 | 0.192 | +1.2 % | 5.0 % |
| tile refine: scan (parallel) | 0.056 | 0.057 | +1.5 % | 5.0 % |
| tile refine: split + flip (serial) | 0.106 | 0.108 | +1.3 % | 6.2 % |
| tile refine: setup + output | 0.025 | 0.026 | +4.0 % | 5.0 % |
| tile t=1 total | 0.699 | 0.711 | +1.7 % | 5.0 % |
| tile t=1 decode | 0.060 | 0.068 | +13.3 % | 5.8 % |
| tile t=1 refine | 0.472 | 0.473 | +0.2 % | 5.0 % |
| tile t=1 refine: scan (parallel) | 0.352 | 0.352 | +0.1 % | 5.0 % |
| tile t=1 refine: split + flip (serial) | 0.094 | 0.094 | -0.1 % | 5.0 % |
| quarter total | 0.355 | 0.366 | +3.1 % | 5.0 % |
| quarter decode | 0.044 | 0.051 | +15.9 % | 5.7 % |
| quarter refine | 0.172 | 0.171 | -0.4 % | 5.6 % |
| quarter refine: scan (parallel) | 0.045 | 0.046 | +1.2 % | 5.0 % |
| quarter refine: split + flip (serial) | 0.107 | 0.106 | -1.4 % | 7.8 % |
| numedalslagen total | 8.829 | 8.894 | +0.7 % | 5.0 % |
| numedalslagen decode | 0.921 | 0.950 | +3.2 % | 5.0 % |
| numedalslagen features clip | 6.000 | 6.022 | +0.4 % | 5.0 % |
| numedalslagen features clip: clean-up | 5.937 | 5.956 | +0.3 % | 5.0 % |
| numedalslagen refine | 0.659 | 0.661 | +0.4 % | 5.0 % |
| lagan total | 14.542 | 14.736 | +1.3 % | 5.0 % |
| lagan decode | 4.578 | 4.615 | +0.8 % | 5.0 % |
| lagan features read | 0.922 | 1.033 | +12.1 % | 5.5 % |
| lagan features clip | 6.342 | 6.353 | +0.2 % | 5.0 % |
| lagan features clip: clean-up | 5.392 | 5.397 | +0.1 % | 5.0 % |
| geilo-al-ramp total | 0.329 | 0.344 | +4.7 % | 5.0 % |
| geilo-al-ramp decode | 0.138 | 0.148 | +7.2 % | 5.0 % |
| geilo-al-ramp refine | 0.099 | 0.101 | +1.5 % | 5.0 % |
| geilo-al-ramp refine: scan (parallel) | 0.083 | 0.083 | +1.0 % | 5.0 % |
| romsdal-slope total | 3.578 | 3.624 | +1.3 % | 5.0 % |
| romsdal-slope refine | 2.009 | 2.010 | +0.1 % | 5.0 % |
| romsdal-slope refine: scan (parallel) | 0.564 | 0.567 | +0.6 % | 5.0 % |
| romsdal-slope refine: split + flip (serial) | 1.184 | 1.183 | -0.1 % | 5.0 % |
| romsdal-slope refine: setup + output | 0.257 | 0.262 | +2.1 % | 5.0 % |

The small `decode` phases (tile, tile t=1, quarter, geilo-al-ramp: +7 to
+16 %) are over their percentage band but under the tool's 0.05 s floor
(0.007 to 0.010 s each), so the tool does not flag them; `lagan decode`,
the one large decode, is +0.8 %.
