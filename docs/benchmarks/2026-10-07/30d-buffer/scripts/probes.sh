#!/bin/zsh
# isolated.py and candidates.py on each catchment, power state before and after
# each, output to raw/probes/. PKG comes from the environment.
D=${0:A:h}; R=${D:h}/raw/probes
PY=/Users/skavhaug/projects/rasputin/.venv/bin/python
export RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data
mkdir -p $R
for c in lagan ljungan_flasjo numedalslagen; do
  for job in "isolated.py $c" "candidates.py $c pieces" "candidates.py $c region" "candidates.py $c scaling"; do
    o=${${job// /_}/.py/}
    pmset -g batt | head -1 > $R/${o}_power.txt
    $PY $D/${=job} > $R/${o}.jsonl 2> $R/${o}_stderr.txt
    echo "exit $?" >> $R/${o}_stderr.txt
    pmset -g batt | head -1 >> $R/${o}_power.txt
  done
done
echo done
