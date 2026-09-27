#!/bin/zsh
# sweep.sh PKG OUT [domain-args...]: phase timings per thread count,
# 3 interleaved passes x 5 calls per process.
R=/Users/skavhaug/projects/rasputin
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad
PKG=$1; OUT=$2; shift 2
pmset -g batt > $OUT.pmset_before
: > $OUT
for pass in 1 2 3; do
  for t in 1 2 3 4 5 6 7 8 10 12 16 20; do
    $R/.venv/bin/python $S/prof_driver.py --pkg $PKG --threads $t --repeat 5 -- mesh \
      --dem $R/tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 "$@" \
      --out $S/sweep.vtk --binary 2>&1 | grep PHASES | sed "s/^PHASES /$pass /" >> $OUT
  done
done
pmset -g batt > $OUT.pmset_after
