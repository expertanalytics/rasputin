# 23a-2 acceptance: `rasputin fetch` (@perf, 2026-10-02)

The increment is "The fetch step and the tile cache" in
`docs/increments/23-basin-scale.md`, and its 23a-2 line under "@perf
acceptance". It adds network and cache-writing code in Python and no C++,
so there is no `bench.py` run here: nothing on `rasputin mesh`'s refine or
mesh path changed. Measured on branch `worktree-23a-2` at `8adb459`,
installed with `pip install ".[codecs]"` into a fresh venv (Python 3.14.7,
tifffile 2026.9.20). Apple M1 Max, 32 GiB, macOS 27.0, on AC power,
battery 100 % and charged (`pmset_*.txt`). Network: the office line to
`opentopography.s3.sdsc.edu` and the GLO-30 bucket. Network figures depend
on that line and on the server, so they are single runs, not medians.

Requests are counted twice: by `rasputin fetch`'s own report (range
requests for blocks) and outside the code under test by `net_audit.py`, a
`sys.addaudithook` probe counting `urllib.Request`, `socket.connect` and
`socket.__new__` events and logging each request's `Range` header.

## Claims

1. **The Velhas piece fetches 160 blocks, 107,977,770 bytes of block data,
   in 14.4 s.** Domain: `bho2017_5k_76949_outline_epsg4674.geojson`
   (EPSG:4674), frame = the source CRS (no `--out-crs`). `--dry-run`'s
   plan and the run agree: 160 blocks, 16 range requests, 107,978,922
   bytes (the 1,152 more than the blocks are the gaps coalesced into
   ranges). The audit counted 19 requests: 3 header prefixes (1, 2 and
   4 MiB; 4 MiB was enough, the design's "8 MiB (measured)" was not
   needed) and the 16 ranges. `--dry-run` wrote nothing (no cache
   directory existed after it).

   | run | blocks fetched | bytes fetched | report requests | audited requests | wall s |
   |---|---:|---:|---:|---:|---:|
   | dry run | 0 (plan: 160) | 0 (plan: 107,978,922) | 0 (plan: 16) | 3 | 7.4 |
   | 1, fresh cache A | 160 | 107,978,922 | 16 | 19 | 14.4 |
   | 2, fresh cache B, SIGKILL at 11 s | 80 present after | - | - | - | 11 (killed) |
   | 3, resume on B | 80 | 54,877,894 | 8 | 9 | 11.2 |
   | 4, rerun on B | 0 | 0 | 0 | 1 | 2.9 |

2. **An interrupted fetch resumes to the same bytes and fetches only what
   was missing.** Run 2 was killed with SIGKILL 11 s in: 80 of 160 blocks
   present, no `.part` file left. Run 3 asked for 8 ranges whose blocks are
   exactly the 80 absent ones (overlap with the present set 0, union all
   160). After it, caches A and B have the same 160 block files, sha256
   equal file by file, every size equal to the header's byte count, and
   the same `header.bin`. Run 4 fetched nothing; its one request is the
   identity re-read of the 4 MiB header prefix (design, "Identity").
   An earlier kill at 7 s (`run2a_killed7s.*`) landed before any block was
   written (0 blocks, 0 `.part`); it is kept as a record, not used.
3. **The cache's blocks are the COG's bytes.** 12 blocks sampled with a
   fixed seed were read again from the COG with `curl -r`, not with
   rasputin's `RangeClient`, and match the cache's files by sha256
   (`verify_cache.out`). The probe can fail: with one bit flipped in one
   block of a copy of cache A, it fails both the A-against-B check and that
   block's COG check (`verify_cache.control.out`).
4. **`NOTICE.txt` equals the catalogue's rendering.** For `anadem-v1` (both
   caches) and `glo30`, the file equals `sources.notice(SOURCES[id])`
   (`check_notice.out`; copies in `NOTICE.*.txt`). With one byte appended
   to a copy, the check fails (`check_notice.control.out`).
5. **No mesh is written from the cache yet, and the refusal is the GeoTIFF
   reader's, after the cache's own checks.** `rasputin mesh --dem anadem-v1
   --cache <B>` on the full cache exits 2 with the reader's GeoKey refusal
   ("ProjectedCSTypeGeoKey (3072) is absent; GeographicTypeGeoKey (2048) =
   4674 ..."), not a cache message (`mesh_full_cache.*`). The cache checks
   that come before it do fire: a `header.bin` with one byte appended gives
   "header.bin is not the one the manifest lists" (`mesh_cacheHdr.out`),
   and a removed `manifest.json` gives "is not in the cache; run: rasputin
   fetch ..." (`mesh_cacheNoMan.out`). The command that refusal prints runs
   as given and finds all 160 blocks present (`refusal_command_dryrun.out`).
   **The block-completeness check is not among the checks the refusal comes
   after:** with one block file removed, the mesh gives the same GeoKey
   refusal (`mesh_hole_cache.*`), because `check(plan)` needs the parsed
   header, and parsing is where the geographic refusal is. This matches
   23a-1's W9 test (manifest and header checks, then the reader), and the
   block check on a geographic source waits for 15c-2.
   GLO-30 behaves the same (`glo30_mesh.out`). The design names no
   projected or GLO-30 mesh case for this acceptance, so none was run.
6. **GLO-30's tile path works against the real bucket.** `--bbox -44.05
   -19.05 -43.95 -18.95` (EPSG:4326) chose 4 tiles, one block each,
   5,904,251 bytes in 4 range requests, 1.8 s wall; the audit counted 9
   requests (the tile list, 4 header prefixes, 4 ranges). A rerun fetched
   nothing, with 5 requests (the tile list and 4 header re-reads).
7. **`rasputin mesh` opened no socket; `rasputin fetch` did.** Every mesh
   run above recorded `urllib.Request=0 socket.connect=0
   socket.__new__=0`. Control: the same mesh run with `--control`, which
   creates one socket first, recorded `socket.__new__=1`
   (`mesh_control.out`), so the hook was live in that process; the fetch
   runs' non-zero counts are the positive case. These mesh runs stop at the
   header refusal, so this covers the mesh path up to there, not a full
   mesh; a full offline mesh from a cache is 23a-1's projected `monkeypatch`
   test and its `socket.socket` check.
8. **The basin box** (`bho2017_level2_76_raw.geojson`, the increment's
   23a-2 line): `--dry-run` plans 8,300 blocks, 4,271,915,858 bytes, 578
   requests, as the design counted (`basin_dryrun.out`). The run fetched
   all 8,300 blocks, 4,271,915,858 bytes, in 578 range requests (581
   audited: 3 header prefixes and 578 ranges), in 266.9 s wall, about
   16 MB/s with the default 8 connections; peak RSS 658 MiB (`/usr/bin/time
   -l`, the probe's Python included). A rerun fetched nothing in 2 s, with
   one request (`basin_run2_noop.*`). The Velhas piece's 160 blocks are
   byte-identical in the basin cache (`velhas_in_basin.out`).

## Files

- `net_audit.py`: the audit-hook wrapper. `verify_cache.py`: claims 2 and
  3, tin_engine-free (tifffile and curl). `check_notice.py`: claim 4.
- `run*.stdout|stderr`, `glo30_*`, `basin_*`: each run's report, progress
  lines, audit counts and requested ranges. `run2_killed.present.txt`: the
  block files present after the kill.
- `manifest.*.json`: the manifests written (Velhas cache B, GLO-30).
- `pmset_*.txt`: power state at the start and around the runs. `dryrun.txt`: the Velhas dry run.

The caches (about 110 MB for Velhas, about 4.3 GB for the basin) were under
the job's scratch directory and are not kept; the commands below recreate
them.

## Reproduce

From the repository root, with `V` a venv holding `pip install ".[codecs]"`
of this tree and `C` an empty scratch directory:

```bash
E=docs/benchmarks/2026-10-02/23a-2-acceptance
D=../rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson
F="anadem-v1 --domain $D --domain-crs EPSG:4674"
$V/bin/rasputin fetch $F --cache $C/A --dry-run
$V/bin/python $E/net_audit.py $V/bin/rasputin fetch $F --cache $C/A
$V/bin/python $E/net_audit.py $V/bin/rasputin fetch $F --cache $C/B & P=$!; sleep 11; kill -9 $P
find $C/B -path '*blocks*' -name '*.bin' | sed "s|$C/B/||" | sort > present.txt
$V/bin/python $E/net_audit.py $V/bin/rasputin fetch $F --cache $C/B 2> resume.stderr
$V/bin/python $E/net_audit.py $V/bin/rasputin fetch $F --cache $C/B       # fetches nothing
$V/bin/python $E/verify_cache.py $C/A $C/B present.txt resume.stderr 12
$V/bin/python $E/check_notice.py $C/A $C/B
$V/bin/python $E/net_audit.py $V/bin/rasputin mesh --dem anadem-v1 --cache $C/B \
  --domain $D --domain-crs EPSG:4674 --tolerance 5 --out $C/v.vtk      # GeoKey refusal, 0 sockets
$V/bin/rasputin fetch glo30 --bbox -44.05 -19.05 -43.95 -18.95 --cache $C/G
$V/bin/rasputin fetch anadem-v1 --domain \
  ../rasputin_data/sao_francisco_piece/bho2017_level2_76_raw.geojson \
  --domain-crs EPSG:4674 --cache $C/basin
```

The kill time depends on the line: at 7 s here no block had landed yet.
