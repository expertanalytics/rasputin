#!/bin/bash
# build.sh: one worktree + Release _core per labelled commit
set -e
B=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m
R=/Users/skavhaug/projects/rasputin
PB=/Users/skavhaug/.cache/uv/archive-v0/a1_2Hvc_PmgI3tx5/pybind11/share/cmake/pybind11
for lc in i14:935271c i14b:333bcf8 i16:1018901 i17:2760d36 i18:a9a93bc i20:f5489ab i20b:75a5b6c head:08a7711; do
  l=${lc%%:*}; c=${lc##*:}; W=$B/wt/$l
  [ -d $W ] || git -C $R worktree add --detach $W $c >/dev/null 2>&1
  cmake -S $W -B $W/b -DCMAKE_BUILD_TYPE=Release -DRASPUTIN_BUILD_PYTHON=ON -DRASPUTIN_BUILD_TESTS=OFF \
    -Dpybind11_DIR=$PB -DPython_EXECUTABLE=$R/.venv/bin/python -DPYBIND11_FINDPYTHON=ON > $W.cfg.log 2>&1
  cmake --build $W/b -j 10 --target _core > $W.build.log 2>&1
  cp $W/b/_core*.so $W/src_python/tin_engine/
  echo "$l $c ok"
done
