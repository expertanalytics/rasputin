#!/bin/zsh
# 30d's speed check (design section 7) on the green step 4f13ee90: Lagan with the
# European GeoPackage, one warm-up and three runs; Numedalslagen (UTM33 GeoPackage)
# once, for its hash. PKG: `git archive 4f13ee90 src_python/tin_engine`, plus the
# Release `_core` from master (C++ unchanged since f81b20b7; SHA-256 86765ff1...).
D=${0:A:h}
for r in w0 r1 r2 r3; do $D/run_mesh.sh built_$r lagan gpkg; done
$D/run_mesh.sh built_r1 numedalslagen gpkg33
echo done
