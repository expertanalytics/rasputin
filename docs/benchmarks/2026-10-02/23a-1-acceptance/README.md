# 23a-1 acceptance: windowed source reads (@perf, 2026-10-02)

The increment is "Windowed source reads (23a-1)" in
`docs/increments/23-basin-scale.md`, and its line under "@perf acceptance".
It changes how DTM10 tiles are decoded (`io/cog.decode_window`, block by
block, one extra copy per window) and no C++. Old is master `5a57793`, new
is the branch head `3b729fb`. They were measured back to back, on AC power
throughout (`pmset -g batt` before and after each part, 68-73 % and
charging), on an Apple M1 Max (8 P + 2 E cores, 32 GiB, macOS 27.0).
Both trees built the same `_core`: sha256 `eb5f06e8...c5d5` in both
`run.json` files.

## Claims

1. **The 1 m benchmark's mesh is unchanged.** `tools/bench.py run` (default
   DEM `tests/fixtures/dem_archive/7908_3_10m_z33.tif`, tolerance 1, threads
   0-20, 5 repeats, interleaved). Mesh sha256 from `POINTS` on: quarter
   `11741a81...005c`, tile `ccebf96a...71a1`, the same on both trees. Quality
   is also the same on both: worst angle 0.6296 / 0.3955 deg, max degree
   74 / 18, within tolerance, 0 Delaunay violations. The verdict from
   `bench.py` is **ACCEPTED** (new against old, 5 % threshold on `refine_s`).
   The old run reports `NO BASELINE`, as expected: the stored runs are on
   battery.
2. **Process time did not get worse.** Median `proc_s`, which is the whole
   child process, over all 105 samples per domain: quarter 0.534 s old,
   0.516 s new; tile 0.586 s old, 0.583 s new. At the default thread
   count: quarter 0.533 s old, 0.513 s new; tile 0.584 s on both.
3. **The multi-tile run gives the identical file, a little faster, with a
   lower peak RSS.** `rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925
   --bbox 799750 7849750 900250 7950250 --tolerance 1 --binary`. This is 15a's
   acceptance box, the 2 × 2 block 7908_3, 7908_2, 7808_4, 7808_1. It reads
   nine tiles and 10051 × 10051 nodes, and inserts 5,025,207 vertices.
   There were 5 rounds, alternating old and new. All 10 output files are
   **byte-identical**: the whole file has the same sha256 (`f86bb64d2e9e9dfd...`),
   and the bytes from `POINTS` on hash to `f9f297f0...63e1`.

   | | wall s (median, min-max) | peak RSS MiB (median, min-max) | refine_s | app_s |
   |---|---:|---:|---:|---:|
   | old 5a57793 | 9.00 (8.95-9.02) | 3026 (3013-3077) | 7.255 | 8.690 |
   | new 3b729fb | 8.78 (8.73-8.87) | 2803 (2802-2807) | 7.271 | 8.475 |

   Peak RSS fell by about 223 MiB (7 %), and the new run varies less. The
   extra copy per window does not show in the peak. The cause of the drop
   has not been profiled.
4. **Meshing never used the network.** (a) All 10 multi-tile runs ran under
   `sandbox-exec -p '(version 1)(allow default)(deny network*)'` and exited 0.
   The control: under the same profile, a Python `socket.create_connection`
   fails with `EPERM`, and without the profile it connects. (b) One more
   new-tree run under `socket_audit.py` (a `sys.addaudithook` probe) recorded
   `SOCKET_EVENTS 0`. Its `--control` run, which opens one socket first,
   recorded `SOCKET_EVENTS 1 ['socket.__new__']`, so the probe can fail.

Not measured here: decode throughput per block and per window at threads
1-8 (the increment's @perf line). The brief for this run did not ask for it.

## Files

- `23a1-base-5a57793/`, `23a1-3b729fb/`: `bench.py`'s `run.json`, `raw.tsv`
  and the generated `README.md`. bench.py's blob hash is in `run.json`.
  The absolute paths in them are the run's own.
- `multitile/`: `run_b.sh` (as run; the master worktree was at
  `$CLAUDE_JOB_DIR/tmp/master`), `run.log`, one `*.err` per run (the child's
  `BENCH` JSON line and `/usr/bin/time -l`), and `trial.log` (one timing
  round run before the five). `run_bench.sh` is the benchmark pair as run.
- `socket_audit.py`: the socket probe.

The meshes (363 MB each) are not kept. Rerun to regenerate them.

## Reproduce

From the repository root, with the project venv's `python`:

```bash
git worktree add --detach ../master-5a57793 5a57793
python tools/bench.py run --label 23a1-base-5a57793 --tree ../master-5a57793 --out-root /tmp/bench
python tools/bench.py run --label 23a1-3b729fb --tree . --out-root /tmp/bench \
    --baseline /tmp/bench/<date>/23a1-base-5a57793
# multi-tile, per tree T (both trees' build-bench/pkg exist after the runs above):
/usr/bin/time -l sandbox-exec -p '(version 1)(allow default)(deny network*)' \
  python tools/bench.py _child --pkg T/build-bench/pkg --threads 0 -- \
  mesh --dem ../rasputin_data/DTM10_UTM33_20260925 \
  --bbox 799750 7849750 900250 7950250 --tolerance 1 --out OUT.vtk --binary
python docs/benchmarks/2026-10-02/23a-1-acceptance/socket_audit.py [--control] \
  tools/bench.py _child --pkg build-bench/pkg --threads 0 -- mesh ...   # same args
```

## Verdict

**ACCEPTED** (AC against AC, back to back). The meshes are identical, the
time is within noise or a little lower, and peak RSS is lower.
