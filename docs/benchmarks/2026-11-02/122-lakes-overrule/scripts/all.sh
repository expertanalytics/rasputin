#!/bin/bash
R=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/bin/runs.sh
NUM=(mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin_data/numedalslagen_outline_nve.geojson --land-cover /Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg --land-cover-layer corine2018 --land-cover-classes corine --edges all --tolerance 10)
LAG=(mesh --dem glo30 --cache /Users/skavhaug/projects/rasputin_data/sweden_glo30_cache --out-crs EPSG:3006 --domain /Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --land-cover /Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --land-cover-layer U2018_CLC2018_V2020_20u1 --land-cover-crs EPSG:3035 --land-cover-classes corine --edges all --tolerance 10)
BYG=(mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin/.claude/worktrees/lakes-overrule/docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson --tolerance 10 --edges all --threads 8 --prefer lake=/Users/skavhaug/projects/rasputin_scratch/temp-check-20261021/nve-lakes-25833.geojson)
# warm-up of file cache, then alternating repeats for Numedalslagen timing
$R num-warm branch "${NUM[@]}"
for i in 1 2 3; do $R num$i master "${NUM[@]}"; $R num$i branch "${NUM[@]}"; done
$R lag-warm branch "${LAG[@]}"
$R lag master "${LAG[@]}"; $R lag branch "${LAG[@]}"
$R byg-clc master "${BYG[@]}" --land-cover clc2018-clean-25833
$R byg-clc branch "${BYG[@]}" --land-cover clc2018-clean-25833
$R byg-inn master "${BYG[@]}" --land-cover ar50-innlandet-clean-25833
$R byg-inn branch "${BYG[@]}" --land-cover ar50-innlandet-clean-25833
$R byg-inn-partial branch "${BYG[@]}" --land-cover ar50-innlandet-clean-25833 --land-cover-partial
$R byg-innvest master "${BYG[@]}" --land-cover ar50-innlandet-clean-25833 --land-cover ar50-vestland-clean-25833
$R byg-innvest branch "${BYG[@]}" --land-cover ar50-innlandet-clean-25833 --land-cover ar50-vestland-clean-25833
echo ALLDONE
