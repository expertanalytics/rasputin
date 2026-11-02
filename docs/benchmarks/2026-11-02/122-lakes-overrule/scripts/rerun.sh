#!/bin/bash
# The 2026-11-02 rerun (bash, so each pair splits into two words): one untimed
# master and branch run per case, then compare.py on each pair into
# raw/compare-<case>.txt. Needs scripts/runs.sh, compare.py and run_master.py in $S/bin.
set -u
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/lakes-overrule
D=/Users/skavhaug/projects/rasputin_data
RAW=$W/docs/benchmarks/2026-11-02/122-lakes-overrule/raw
R=$S/bin/runs.sh
NUM=(mesh --dem $D/DTM10_UTM33_20260925 --domain $D/numedalslagen_outline_nve.geojson --land-cover $D/corine2018_dtm10_utm33.gpkg --land-cover-layer corine2018 --land-cover-classes corine --edges all --tolerance 10)
LAG=(mesh --dem glo30 --cache $D/sweden_glo30_cache --out-crs EPSG:3006 --domain $D/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --land-cover $D/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --land-cover-layer U2018_CLC2018_V2020_20u1 --land-cover-crs EPSG:3035 --land-cover-classes corine --edges all --tolerance 10)
BYG=(mesh --dem $D/DTM10_UTM33_20260925 --domain $W/docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson --tolerance 10 --edges all --threads 8 --prefer lake=/Users/skavhaug/projects/rasputin_scratch/temp-check-20261021/nve-lakes-25833.geojson)
INN=(--land-cover ar50-innlandet-clean-25833)
for side in master branch; do
  $R num $side "${NUM[@]}"
  $R lag $side "${LAG[@]}"
  $R byg-clc $side "${BYG[@]}" --land-cover clc2018-clean-25833
  $R byg-inn $side "${BYG[@]}" "${INN[@]}"
  $R byg-innvest $side "${BYG[@]}" "${INN[@]}" --land-cover ar50-vestland-clean-25833
done
$R byg-inn-partial branch "${BYG[@]}" "${INN[@]}" --land-cover-partial
cmp() { "$W/.venv/bin/python" "$S/bin/compare.py" "$S/out/$2-master.vtk" "$S/out/$3-branch.vtk" 2>&1 | grep -v RuntimeWarning > "$RAW/compare-$1.txt"; echo "compare $1 status ${PIPESTATUS[0]}"; }
cmp numedalslagen num num
cmp lagan lag lag
cmp bygdin-clc byg-clc byg-clc
cmp bygdin-innlandet-vestland byg-innvest byg-innvest
cmp bygdin-innlandet-partial byg-inn byg-inn-partial
echo ALLDONE
