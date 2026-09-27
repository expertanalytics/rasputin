UBSan + libc++ extensive hardening runs of the evaluate-once simulation (2026-09-27, 21:18-21:20, battery, awake).
Build: the scratch worktree (sim.patch + sim_a0.patch + sim_a0_eo.patch), RelWithDebInfo with
  -fsanitize=undefined -fno-sanitize-recover=all -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE -fno-omit-frame-pointer -g
Run: bench.py _child --threads 1 under `lldb --batch` (a trap stops in the debugger with a stack, no crash dialog).
ASan was not used: the harness strips DYLD_INSERT_LIBRARIES, so the ASan runtime cannot be preloaded into Python.
Result: A0eot, A0eo, A0eofp on quarter and tile, and A0eot with RASPUTIN_SIM_EOVERIFY on quarter: exit 0, no
  runtime error, rounds/inserted/flips as the Release runs. plant_stalecommit.log: the planted defect traps
  (EXC_BREAKPOINT, a libc++ hardening check) in legalise_around, called from the plant's commit
  (refine.hpp:1229 in the patched tree): a stale plan was split.
