# Increment 24 acceptance: release hardening (@perf, 2026-10-03)

Branch `worktree-agent-a08216a68d5534b18` at `9ce01ee` (increment 24: the
CMake option `RASPUTIN_HARDENING`, default ON; `@reviewer` round 3 APPROVED),
against the previous merge, origin/master's tip `6cdc8cc`. `6cdc8cc` differs
from the branch's merge base `586fbc1` only in two retrospective files
(`git diff --stat 586fbc1 6cdc8cc`), so it is the code base the branch was
made from. The acceptance form is the design's §7
(`docs/increments/24-release-hardening.md`): the branch built OFF must match
the previous merge, and the same branch built ON, measured back to back with
OFF, becomes the first hardened baseline. The ON/OFF difference is the
hardening cost.

## Verdict

**ACCEPTED.**

1. **Nothing else moved.** The OFF build and the previous merge produce
   byte-identical meshes and the same quality and refine counters. Refine
   time agrees within 5 % at every thread count in both clean pairs (below):
   OFF against master is -2.9 % to +1.6 % on refine, with one cell at +3.9 %
   (quarter, default threads). The `_core` binaries differ only because the
   branch adds the `hardening` attribute to the bindings.
2. **Hardening cost, ON against OFF:** refine +4.2 % to +5.9 % at 1 thread
   and +3.2 % to +11.1 % over 2 to 20 threads (most cells +5 % to +8 %). In
   process (`app_s`: read, mesh, write) at the default thread count the cost
   is +1.3 % to +5.5 %, and for the whole child process -1.1 % to +2.9 %.
   This reproduces the study's libc++ FAST figures
   (`release-hardening/README.md`: refine +4.9 / +4.5 % at 1 thread, +4.4 to
   +8.0 % over 2 to 20 threads, `app_s` +4.0 / +4.1 %, process +1.9 /
   +2.6 %). The meshes are identical ON and OFF.
3. **The new hardened baseline** is `24-on-r2/` (AC, `libc++ fast`, `_core`
   `0b3362e5…`). A later `bench.py run` with the default `--hardening on`
   finds it on its own, because it is the newest comparable stored run.
   `24-on/` is the same build, from the disturbed pair.

## Method

- `tools/bench.py run`, blob `2b7dca1aa1184afea53b33c357ab3b836403d2f3` (the
  branch's; it drove every run, including master's through `--tree`).
  Release `_core` built by `bench.build()` into each tree's `build-bench/`.
  Defaults otherwise: DEM `tests/fixtures/dem_archive/7908_3_10m_z33.tif`,
  tolerance 1 m, domains `tile` and `quarter`, threads 0 (the CLI default) and
  1 to 20, 5 repeats interleaved over the thread counts.
- Apple M1 Max (8 P + 2 E, 32 GiB), macOS 27.0, Python 3.14.7, numpy 2.5.3,
  AppleClang 21.0.0. `caffeinate -i` was held during each batch. No other
  agent ran and no other build ran during the measurements.
- Three trees, each built once before its batch:
  - M: a detached worktree of `6cdc8cc`, built with `--hardening off`. That
    tree's CMake ignores the define and warns. Its `_core` (`4b5f875d…`) has
    no `hardening` attribute, so the child reports `none`.
  - OFF: a detached worktree of `9ce01ee`, `--hardening off`, `_core`
    `42ff8b45…`, reports `none`.
  - ON: the branch worktree itself, `--hardening on`, `_core` `0b3362e5…`,
    reports `libc++ fast`.
  - The probe that the objects differ as claimed: `objdump -d` counts 170
    `brk` traps in M's and OFF's `_core` and 321 in ON's (+151; the study
    found +149 on `390b516`). Importing each `pkg/` the way the child does
    (the editable finder removed) gives `_core.hardening` absent, `'none'`
    and `'libc++ fast'`.
- Order per batch: M, OFF, ON, ON2, OFF2, M2. That gives two pairs for each
  comparison, (M, OFF) with (OFF2, M2), and (OFF, ON) with (ON2, OFF2), with
  the order balanced against drift. The driver is
  `24-acceptance/pairs.sh W M OFF S`, and the tables come from
  `24-acceptance/summarize.py [DIR]`.
- The sanitize-first rule does not apply: no C++ was patched for this run.

### Two batches, and which pairs are clean

- **Batch 1** (17:11 to 17:24 UTC, `24-acceptance/batch1/`): the first three
  runs ran undisturbed on AC (charged, 100 %). From 17:17:52 the keyboard and
  trackpad were in use (`batch1/pmset-activity.txt`, from `pmset -g log`).
  At 17:22:25 the machine went to battery, during M2, which `bench.py`
  recorded as `mixed` and refused to judge. Same-binary drift inside the
  batch was +9 to +11 % (ON2 against ON). The batch is kept as evidence. Only
  its first pair (M, OFF, ON) is used.
- **Batch 2** (17:52 to 18:06 UTC, the run directories directly under
  `2026-10-03/`) started automatically once AC had held for 120 s. AC power
  for all six runs: 97 to 98 %, charging, `pmset -g batt` before and after
  each run, in each `run.json`. Its **first pair is disturbed**: M and M2 are
  the same binary, yet M is +6 to +7 % slower at 1 thread, and the AC
  reconnect and keyboard activity (17:50:30 and 17:51:33 UTC) were minutes
  before. The cause was not measured. Its **second pair (ON2, OFF2, M2) is
  clean**: every same-side check agrees.
- The verdict rests on the two clean pairs: batch 1 (M, OFF, ON) and batch 2
  (ON2, OFF2, M2). The full batch-2 tables, both pairs and the pooled column,
  are in `24-acceptance/tables.md`, and batch 1's are in
  `24-acceptance/batch1/tables.md`. In batch 2 the pooled column (10 samples
  a side) keeps OFF against M within 5 % at every cell (-2.5 % to +3.6 %).

## Results (refine_s, change of medians, %)

| comparison | pair | domain | default | 1 thread | 2 to 20 threads |
|---|---|---|---:|---:|---:|
| OFF vs M | batch 1 (M, OFF) | tile | -0.1 | -0.8 | -1.7..+0.2 |
| OFF vs M | batch 1 (M, OFF) | quarter | +0.1 | -0.4 | -2.9..+0.2 |
| OFF vs M | batch 2 (OFF2, M2) | tile | -0.3 | +0.6 | -2.6..+1.6 |
| OFF vs M | batch 2 (OFF2, M2) | quarter | +3.9 | +0.5 | -1.3..+1.4 |
| ON vs OFF | batch 1 (OFF, ON) | tile | +7.2 | +5.2 | +4.3..+7.7 |
| ON vs OFF | batch 1 (OFF, ON) | quarter | +6.2 | +5.9 | +3.2..+11.1 |
| ON vs OFF | batch 2 (ON2, OFF2) | tile | +7.5 | +4.2 | +5.3..+7.8 |
| ON vs OFF | batch 2 (ON2, OFF2) | quarter | +3.3 | +5.1 | +5.4..+7.4 |
| *for reference, disturbed:* OFF vs M | batch 2 (M, OFF) | tile / quarter | +3.6 / -2.0 | +6.7 / -1.1 | -0.4..+10.3 / -1.9..+11.5 |
| *for reference, disturbed:* ON vs OFF | batch 2 (OFF, ON) | tile / quarter | +4.7 / +17.9 | +0.5 / +14.0 | -2.2..+12.9 / -4.2..+14.5 |

The disturbed pair's deviations go both ways, and the same binary moves by
the same amount (M against M2: +6.0 % tile and +7.3 % quarter at 1 thread),
which is why it is set aside rather than read as a code effect. Hardening
cost end to end, clean pairs, default threads: `app_s` +5.5 / +5.4 % (tile)
and +3.8 / +1.3 % (quarter); `proc_s` +2.9 / +2.4 % (tile) and +1.7 /
-1.1 % (quarter), each batch 1 / batch 2.

**Scaling ceiling** (1 thread over 20, best in brackets, from each run's
README), AC: ON2 tile 2.48x (2.51x at 11), quarter 2.46x (2.50x at 10); OFF2
tile 2.52x (2.58x at 12), quarter 2.50x (2.52x at 10); M2 tile 2.49x (2.54x
at 9). The study's AC figures were plain 2.49x and hardened 2.46x (tile,
best). The checks cost a few hundredths of the ceiling, as before.

## Quality (identical in all twelve runs, both batches)

| domain | worst angle | max degree | within tolerance | Delaunay violations | mesh sha256 |
|---|---:|---:|---|---:|---|
| tile | 0.6296° | 74 | yes | 0 of 692,056 | `11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c` |
| quarter | 0.3955° | 18 | yes | 0 of 641,791 | `ccebf96a86c6c5e244e4a0281919de4e866fcfe789b66024a290ac2af33771a1` |

These are the same hashes as the study and the 15f acceptances. The batch-2
ASCII meshes were compared with `cmp` against M's, per domain: all five
identical. The probe can fail: the tile and quarter meshes differ at line 5.
Over all six runs of each batch, every (domain, thread count) cell has a
single (max_error, rounds, inserted, flips) tuple (42 cells), so the checks
and the OFF build changed no decision in refine.

## `bench.py`'s mode handling, checked

- **`run.json` records the mode, read from the loaded module.** The
  `"hardening"` value is `none` for M, OFF and OFF2 and `libc++ fast` for ON
  and ON2, and each README's build line ends with `bounds checks: <mode>`.
  For M the value comes from the absent attribute.
- **`comparable()` refuses across modes.** ON was run with
  `--baseline 24-off`, and in both batches `bench.py` printed
  `NO BASELINE: hardening: libc++ fast vs none` and exited 2 (`pairs.log`).
  `bench.py compare` gives the same answer offline on the study's runs:
  h2-fast against h1-plain gives `hardening: libc++ fast vs none`, and
  g2-gassert against h2-fast gives
  `hardening: libstdc++ assertions vs libc++ fast`, both exit 2. Within a
  mode the comparison goes through: h3-fast against h2-fast is `ACCEPTED`.
  Before the study's runs were re-stamped (next section), h2-fast judged
  against h1-plain printed 39 `REGRESSION` lines
  (`release-hardening/raw/h2-fast/README.md`), so the refusal is the change
  that made the difference.
- `bench.py`'s own verdicts in batch 2: OFF2 `ACCEPTED` against M, M2
  `ACCEPTED` against OFF2, ON2 `ACCEPTED` against ON. OFF against M printed
  7 `REGRESSION` lines, from the disturbed pair. M had no `--baseline`. It
  judged itself against batch 1's M, from a different machine state, and
  printed 7 `REGRESSION` lines. That shows the session drift, not a change.

## The study's hardened runs, re-stamped

The six hardened runs in `release-hardening/raw/` (`h2-fast`, `h3-fast`,
`h6-fast`: libc++ FAST; `g2-gassert`, `g3-gassert`, `g6-gassert`:
`_GLIBCXX_ASSERTIONS`) were made before `RunRecord.hardening` existed. They
therefore loaded as `none`, the mode `bench.py` gives an old `run.json`,
which is wrong for exactly these six. Their `run.json` now carries
`"hardening": "libc++ fast"` or `"libstdc++ assertions"`, the strings
`terrain::stdlib_hardening()` reports for those builds
(`include/terrain/build_info.hpp`). Each was identified by its `_core`
sha256: `bcfda7ba…` for FAST and `f058d5ef…` for the GCC build with
assertions, as listed in the study's Method. The edit was a load, set and
dump through `bench.py`'s own `RunRecord`. A check first confirmed that a
load-and-dump of all twelve files reproduces them byte for byte apart from
the new field, so the diff is one added line per file. The plain runs and
every README are unchanged. Without the stamp, a future unchecked run could
take a hardened study run as its baseline.

## Files and clean-up

- `24-base-6cdc8cc/`, `24-off/`, `24-on/`, `24-on-r2/`, `24-off-r2/`,
  `24-base-6cdc8cc-r2/`: batch 2, in `bench.py`'s format.
- `24-acceptance/pairs.sh`, `summarize.py`, `tables.md`, `pairs.log`
  (batch 2's log with `pmset` per run). `24-acceptance/batch1/`: batch 1's
  six runs, its log, tables and the `pmset -g log` activity extract.
- The ON runs are marked dirty because of the untracked evidence directories
  and the six re-stamped `run.json` files in the branch worktree. No source
  had changed.
- The meshes (about 920 MB of ASCII VTK) were in this session's scratchpad
  and have been deleted, along with the two detached worktrees and their
  `build-bench/`. The branch worktree's `build-bench/` was removed too. To
  regenerate them, create the detached worktrees, build each tree with
  `bench.build(bench.make_runner(), Path(T), "on"|"off")` as in
  `15f-2-acceptance.md`, and run `pairs.sh`. The meshes then land in
  `S/meshes/<label>/`.
