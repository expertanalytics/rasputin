# Increment 137, section 8: `--dem-heights` timings and peak memory

Reported, not gated (`docs/increments/137-dem-stats.md` section 8; an
option off by default). Head `4d25f35a` (branch of worktree `dem-stats`)
against public master `98fe80c3`, 2026-11-04.

## Result

With the flag off, head runs as master does, and writes the same bytes.
With it on, the whole fixture tile takes 0.07 s more (the design predicted
+0.06 s) and 10 MB more peak memory. The cut run takes +0.03 s and 2 MB.
Velhas shows no change we can measure. The flag leaves every mesh as it is
without it (same vertices and triangles, in all 6 files).

| case | triangles | master, off | head, off | head, on | on minus off |
|---|---|---|---|---|---|
| tile: wall time | 69,571 | 0.46 s | 0.45 s | 0.52 s | +0.07 s |
| tile: in-process time | | 0.194 s | 0.195 s | 0.264 s | +0.069 s |
| tile: peak resident | | 341 MB | 342 MB | 352 MB | +10 MB |
| cut4: wall time | 31,237 in 4 files | 0.40 s | 0.40 s | 0.43 s | +0.03 s |
| cut4: in-process time | | 0.139 s | 0.138 s | 0.163 s | +0.025 s |
| cut4: peak resident | | 134 MB | 133 MB | 135 MB | +2 MB |
| velhas: wall time | 969 | 0.32 s | 0.32 s | 0.32 s | 0 |
| velhas: in-process time | | 0.065 s | 0.064 s | 0.065 s | +0.001 s |
| velhas: peak resident | | 91 MB | 91 MB | 92 MB | +1 MB |

These are medians of 5 interleaved repeats; every range is in
`samples.jsonl`. Wall time is `/usr/bin/time`'s real time for the whole
process, interpreter start included. In-process time is
`cli.app`'s own time. The extra time is measured, not guessed: it is the new
`dem heights` phase in the mesh's phase clock (median 0.065 s on the tile,
0.023 s on the cut, 0.001 s on Velhas). Encoding the file takes +0.001 s on
the tile. No other phase moves by more than 3 ms.

The design's memory estimate (section 7) for the tile is 27 bytes per
triangle, twice that while the block buffers are joined, plus the 72-byte
chunk of up to 2^18 triangles: about 9 MB. We measured +10 MB. On Velhas
the source tile kept through the final check is too small to show. The
piece budget's estimate is unchanged, as the design says.

| file | mean classes per triangle | bytes off | bytes on |
|---|---|---|---|
| tile | 3.79 | 2,643,736 | 4,534,275 |
| cut4 pieces 0-0, 0-1, 1-0, 1-1 | 1.00, 5.74, 7.72, 7.08 | 1,179,395 in all | 2,451,709 in all |
| velhas | 3.99 | 37,531 | 65,319 |

## Cases

- **tile:** `rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif
  --tolerance 10 --out X.vtk` (the design's fixture command).
- **cut4:** the same DEM with `--domain docs/benchmarks/2026-09-26/quarter.geojson
  --tolerance 10 --pieces 4`. `--pieces` needs a `--domain`, and the design
  names none, so this uses the quarter domain already committed for this
  fixture.
- **velhas:** `--dem tests/fixtures/velhas/anadem_velhas.tif --domain
  tests/fixtures/velhas/catchment.geojson --out-crs EPSG:31983 --tolerance 5`
  (the reprojected path, test 11's command).

Each is run with master and the flag off, head and the flag off, and head and
`--dem-heights`. Output is binary `.vtk`, with the CLI's default threads.

## Method

- Machine: Apple M1 Max, 10 cores, 32 GB. Power: `pmset -g batt` said "Now
  drawing from 'AC Power'" (100 %, charged) before and after the run.
- Builds: both trees are `CMAKE_BUILD_TYPE=Release`, `RASPUTIN_HARDENING=ON`
  (`_core.hardening` "libc++ fast" in every sample), and use the same
  compiler and the same venv's Python 3.14 and pybind11. Head's `_core` was
  rebuilt by `tools/rebuild_core.py` (sha256 `61f1c0ae9d63...`). Master is
  a detached worktree of `98fe80c3` with `_core` built in its own
  `build-pyext`.
- Each sample is one `tools/bench.py _child --pkg <pkg> --threads 0 -- mesh ...`
  process (head's `bench.py`) under `/usr/bin/time -l`. `<pkg>` is a symlink
  tree: each entry of the tree's `src_python/tin_engine/` except
  `__pycache__` and `_core*.so`, plus that tree's built `_core`. The child
  removes the editable finder, so each tree's code is the code that runs.
  The repeats go round by case and variant, and one warm-up pass came first.
- Byte identity: `meshes.json` holds each output's sha256 (a cut run's is
  taken over its piece files and `index.json`). Off on master and off at
  head agree for all three cases. The same file has, per head file, the
  triangle count and whether on and off have the same vertices and
  triangles (`same_mesh`).
- The sanitizer runs for this C++ change (ASan with UBSan on the two touched
  suites) are `@developer`'s, recorded in the increment's As built.
- No refine or mesh code changed (design, header), so the quality checks
  (worst angle, vertex degree, tolerance, Delaunay) were not run. The meshes
  are the same with the flag off or on, and off is byte for byte master's.

## Regenerate

```bash
git worktree add --detach <scratch>/base 98fe80c3
cmake -S <scratch>/base -B <scratch>/base/build-pyext -DCMAKE_BUILD_TYPE=Release \
  -DRASPUTIN_HARDENING=ON -DRASPUTIN_BUILD_PYTHON=ON -DRASPUTIN_BUILD_TESTS=OFF \
  -DPython_EXECUTABLE=<head>/.venv/bin/python -Dpybind11_DIR=<pybind11 cmake dir of that venv>
cmake --build <scratch>/base/build-pyext -j --target _core
# pkg-base, pkg-head: symlink trees as in Method
<head>/.venv/bin/python run_137.py --repo <head> --scratch <scratch> --repeats 5
```

Meshes go to `<scratch>/out` and are not committed.
