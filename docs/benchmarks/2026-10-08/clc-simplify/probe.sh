#!/bin/zsh
# Land-cover simplification probe: --features-tolerance (FT) 0, 10, 30, 100 m,
# each at vertical tolerance 50, 20, 10 m with --start-min-angle at the default
# (25) and at 0, plus the "minimal" mesh (--tolerance 1e6 --start-min-angle 0).
# German fused catchment (Isar above Kruen + Walchensee + Loisach above Kochelsee),
# CORINE 2018 window, GLO-30. Same inputs as
# rasputin_scratch/germany/isar-loisach/run_meshes_clc.sh.
# Usage: zsh probe.sh [WORKTREE] [OUTDIR]; meshes go to OUTDIR, the --stats
# files and *_errors.json are copied next to this script.
set -e
HERE=${0:A:h}
W=${1:-/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify}
O=${2:-/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify}
D=/Users/skavhaug/projects/rasputin_data
K=/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach
S=$D/lands/scripts
C=$D/germany_corine/isar_loisach_clc2018_3035.geojson
TILES=($D/germany_glo30/Copernicus_DSM_COG_10_N47_00_E010_00_DEM.tif $D/germany_glo30/Copernicus_DSM_COG_10_N47_00_E011_00_DEM.tif)
export RASPUTIN_DATA=$D
P=/Users/skavhaug/projects/rasputin/.venv/bin/python  # has matplotlib, which the helper needs
mkdir -p $O
pmset -g batt | head -1
for ft in 0 10 30 100; do
  for t in min 50 20 10; do
    for a in 25 0; do
      [[ $t == min && $a == 25 ]] && continue
      if [[ $t == min ]]; then tol=(--tolerance 1e6 --start-min-angle 0); n=ft${ft}_minimal
      else tol=(--tolerance $t --start-min-angle $a); n=ft${ft}_tol${t}_a${a}; fi
      echo "== $n"
      $W/.venv/bin/rasputin mesh --dem glo30 --cache $D/germany_glo30_cache --out-crs EPSG:25832 \
        --domain $K/fused_catchment_glo30.geojson --features $C --features-crs EPSG:3035 \
        --features-map corine --features-tolerance $ft $tol \
        --out $O/$n.vtk --stats $O/${n}_stats.md --record $O/${n}_record.json >/dev/null
      $P $S/mesh_error_stats.py $O/$n.vtk $K/fused_catchment_glo30.geojson $TILES \
        --json $O/${n}_errors.json >/dev/null
      cp $O/${n}_stats.md $O/${n}_errors.json $HERE/
    done
  done
done
