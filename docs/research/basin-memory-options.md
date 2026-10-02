# Fitting the São Francisco basin in memory: where it goes, and the options

Status: research note, 2026-10-02, `@architect`, on `worktree-basin-memory`
(from `worktree-15c-2`), updated the same night with @perf's measurements.
Analysis only: the choice between the options is Ola's.

Each figure is marked:
- **[m]**: measured, with its record;
- **[d]**: derived by arithmetic from measured figures;
- *est.*: worked out from the code alone.

Evidence:
- the uncut basin run, `docs/benchmarks/2026-10-02/basin-anadem/`;
- per phase on level-3 sub-basin 761, `docs/benchmarks/2026-10-02/basin-phases/`;
- the whole basin as nine level-3 sub-basins, `docs/benchmarks/2026-10-02/basin-level3/`;
- one resample block, `docs/benchmarks/2026-10-02/basin-memory-probe/`.

## 1. Where the memory goes

Sub-basin 761 has a 0.538 G-node target grid (27,786 columns) and 270 M
check points. Its peak memory footprint, phase by phase, in GB [m]:

| phase | what holds it (code) | as shipped | `MallocLargeCache=0` |
|---|---|---:|---:|
| P1 load + assemble | the source mosaic, 4 B per node held. While loading, 12 B per node: three copies, one not identified | 7.41 | 6.81 |
| **P2 resample** | mosaic + target canvas + `os.cpu_count()` (10) blocks of 256 rows in flight, **133-135 B per node** of a block (`target_grid.py`, `resample`, `_pool`) | **17.8-18.3** | **14.0-14.1** |
| P2 end | the canvas copied by the public `DemTile(...)` (15c D3, "Built otherwise (15c-2)") | 10.2-10.4 | 6.63 |
| P3 refine | mosaic (held by the check-point iterator) + target tile + mesh | 10.1-10.3 | 4.7-5.0 |
| P4-P5 store, phase 2 | the store: 16 B per point live, **30 B charged** (small-zone fragmentation; `freeze` does not give it back). The mosaic is freed once the iterator is exhausted. The target tile stays (15c D6, "Built otherwise (15c-2)") | 11.1-12.2 | under-counted |

Three points:

- **The peak is set in resampling, at every tolerance**: it moved by less
  than 0.5 GB from 50 m to 5 m [m].
- macOS keeps freed large blocks charged to the process. The environment
  variable `MallocLargeCache=0` returns them [m], with **no code change**.
- The mesh costs less than this note first assumed (2 × 310 B per
  triangle). On 761, phase 2 measured ≥ 170-220 B per final triangle
  [m, a lower bound]. At 2 m, with binary output, unit 769 (43.3 M
  triangles) peaked at 10.76 GB, set by the mesh, and 761 (27.1 M) at 8.81 GB
  late in the run [m]. Less canvas and store, that is **~90-230 B per
  triangle** [d]. The range comes from how much of the store the footprint
  saw: between 16 B per point and the fit's 1.5 B.

**The whole basin uncut** (2.06 G target nodes, 40,943 columns, ~0.82 G check
points), extrapolated from 761 [d]:
- P2 needs **31 GB live, ~47 GB as shipped**;
- P4 needs 29.8 GB live, ~41 GB as shipped.

As shipped it fits 32 GiB at no tolerance. That agrees with the 33.2 GB
measured when the run was stopped in its resampling phase.

## 2. The options, and what each reaches

**(c) Pieces: separate sub-basins, measured.** No production lines: a loop
over `rasputin mesh --domain <unit>` with one shared `--out-crs`, and
`MallocLargeCache=0`. The nine BHO level-3 units give the whole basin [m]:

| tolerance | 20 m | 10 m | 5 m | 2 m |
|---|---:|---:|---:|---:|
| triangles | 5.0 M | 13.1 M | 34.6 M | 120.2 M |
| max peak, all units | 13.9 GB | 14.4 GB | 14.3 GB | 14.3 GB |

**So the basin at 2 m runs on 32 GiB today.** What is missing is **matching
edges**. Each unit splits its border on its own (the edge strip and phase 2,
J9), so neighbouring meshes have hanging vertices, and z along a border can
differ by up to ~2T. The units are BHO outlines, so the result is a
measurement and not a product: B1 drops BHO as geometry and admits
sub-catchments "never as seams"; B13 allows the BHO outlines only as
measurement domains. Risk: the largest unit sets the peak (761 holds a third
of the basin); smaller units would lower it.

**(b) Increment 23's decomposition, 23b-23g** (~1,555 lines estimated,
~2,160 at +39 %; plus 23e). It supplies exactly what (c) lacks:
- seams on lattice lines, with a seam pass both pieces agree on (23b-23c);
- a stitched file with the seams removed (23f-23g);
- pieces in parallel, which takes refine's serial phase off the critical
  path (23d);
- pieces cut by the budget, not by BHO.

Memory is no longer the open question: (c) shows that pieces of ≤ 0.54 G box
nodes fit at ~14 GB. What remains for 23 is conformity and speed.

Risk: **`b(T)` does not match this path.** On 761 the peak is 26 B per box
node at every tolerance [d] (14.1 GB / 0.538 G), set by the resample
transient, which grows with a piece's columns × threads, not with its area.
`b(T)` gives 19 B at 50 m (too low) and 65 B at 5 m (far too high, which
would make pieces needlessly small). It should be re-measured on the
reprojected path after the fixes in (a).

**(a) Stream within one process.** These fixes in `resample`, `DemTile`, the
store and the loader cut the measured costs; each also lowers every piece's
peak under (b) and (c):
1. resample blocks sized by node count (~3 lines): removes most of the
   ~7 GB transient on 761 [m];
2. adopt the canvas instead of copying it (~3 lines; 15a's M15 test reserves
   `_adopt` for `mosaic.py`, so that test changes too);
3. drop the target tile before phase 2, as D1 intends (~15 lines);
4. reserve the store's rows from a counting pass, or allocate them as one
   block (~20 lines with C++): 30 → 16 B per point [m];
5. find and remove P1's third copy (unknown lines; 12 → 4-8 B per source node
   while loading [m]).

Then a cache-backed `SourceWindows`, a windowed `assemble` per window, which
honours 23a-1's decision 1 and J7 (~60-90 lines), so the mosaic never exists.
Items 1-4 plus the stream come to ~100-130 lines.

What remains for the uncut basin [d]:
- phase 1: canvas 8.2 GB + mesh;
- phase 2: store 13.1 GB live + mesh at ~90-230 B per triangle.
- 5 m: **~16-21 GB**, which fits.
- 2 m: 120 M triangles, so 11-28 GB of mesh plus the store, **~24-41 GB**,
  which is marginal to no.

(a) keeps all of 23: D1 names a cached `SourceWindows` as "a change of
provider, not of algorithm", and 23c's per-piece windows need that provider.
Risks: more decoding (an LRU of decoded blocks bounds it); one serial phase;
2 m is marginal uncut.

**(d) Other routes.**
- (d1) `MallocLargeCache=0`, no code: 17.8 → 14.0 GB on 761 [m]. As an
  environment setting it is documentation, not a product fix. With it, the
  footprint under-counts the store [m], so a guard that watches the footprint
  sees less than is live.
- (d2) A ≥ 64 GB machine, 0 lines; legitimate under B14 ("on limited
  resources, you should not expect to reproduce huge catchments").
- (d3) `np.memmap` of the large arrays: it hides the cost rather than
  removing it. Not pursued.

Prior art: out-of-core and streaming Delaunay construction exists (Isenburg,
Liu, Shewchuk, Snoeyink, SIGGRAPH 2006; Agarwal, Arge, Yi, ESA 2005;
TerraStream, Danner et al., ACM GIS 2007; cited from memory, **unverified**).
None of them is a tolerance-refined TIN checked against a reprojected
source, so (a) and (b) are engineering here, not novelty.

## 3. What a 2-5 m tolerance makes necessary at basin scale

At 2-5 m, memory alone no longer forces anything: (c) meshes the whole basin
at 2 m in 14.3 GB [m]. What 2-5 m needs is **matching edges between the
pieces**, and only (b) provides them. (a) is not needed to make the basin
fit as pieces. It makes the uncut basin plausible at 5 m [d], and it lowers
every piece's peak by about half. That raises the piece size 23 can use, so
fewer seams are needed.

| option | prod. lines | basin at 5 m / 2 m on 32 GiB | edges match | keeps 23 | main risk |
|---|---:|---|---|---|---|
| (a) stream + fixes, uncut | ~100-130 (+ P1 copy) | ~16-21 GB / ~24-41 GB, marginal [d] | one mesh | yes, prerequisite of 23c | one serial phase; decoding |
| (b) 23b-23g | ~1,555-2,160 | yes / yes (pieces fit, [m] via (c)) | yes, after 23g | is 23 | time to first conforming mesh; `b(T)` miscalibrated |
| (c) BHO level-3 pieces | 0 | **14.3 GB / 14.3 GB [m]** | no (hanging vertices, ≤ ~2T) | cache and fetch | not a product (B1); largest unit sets peak |
| (d1) `MallocLargeCache=0` | 0 | (with (c): as above) | n/a | yes | footprint under-counts the store |
| (d2) larger machine | 0 | yes / likely [d] | one mesh | yes | cost |

**Recommendation (the decision is Ola's):**
1. **Now:** treat the nine level-3 meshes, run with `MallocLargeCache=0`, as
   the basin's interim 2-5 m result, labelled as non-conforming at unit
   borders.
2. **Next:** (a)'s fixes 1-4 as one small PR (~40 lines). They halve every
   piece's peak whichever route follows, and re-measuring `b(T)` follows from
   them. The streaming `SourceWindows` can wait until 23c needs it.
3. **Then:** 23b-23g, which add matching edges and parallel speed. A
   conforming 2-5 m basin mesh needs them; the memory budget does not.
