#!/bin/bash
B=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m
for l in i14 i14b i16 i17 i18 i20 i20b head; do $B/bench.sh tile $l $l; done
$B/bench.sh tile head_nofeet_q0 head --start-min-angle 0 --no-constraint-feet
