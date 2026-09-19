# Automatic Catchments and River Networks from DEM

Design sketch for the feature-extraction stage that produces the constraint inputs (catchment polygons, river/creek polylines) consumed by the noding + CDT pipeline in `parallel_refinement.md`.

## Goal

Given a DEM raster and a seed pixel (or set of pixels), output:

- **Catchment polygon** — the closed boundary of every pixel that drains to the seed.
- **River and creek polylines** — vectorized stream network within the catchment, ordered so creeks and minor rivers are distinguishable from major ones.

Both feed directly into the constraint noding step. Catchment boundaries provide the outer constraint polygon for the mesh; river polylines provide constraints carrying the `river` property.

## Pipeline

1. **DEM conditioning.** Detect sinks (pits) and remove them by filling or breaching so flow routing can complete to a domain outlet.
2. **Flow direction.** Per pixel, encode the direction(s) water leaves the cell.
3. **Flow accumulation.** Per pixel, count the upstream contributing area.
4. **Catchment delineation.** From a seed pixel, walk the flow-direction graph upstream and collect every contributing pixel.
5. **Stream extraction.** Threshold flow accumulation to mask stream pixels; vectorize the mask into polylines; classify by stream order.
6. **Output.** Polygonize the catchment mask to a boundary polyline (simplified upstream of the noder); emit stream polylines with stream-order metadata so callers can pick a creek/river cutoff.

## DEM conditioning

Raw DEMs contain artifact pits — single-cell or multi-cell depressions with no outflow — from sensor noise and resampling. Flow routing cannot escape them, so they must be removed first.

Two families of methods:

- **Filling** (Wang-Liu / priority-flood) — raises pit pixels to the lowest outlet elevation. Simple, preserves general topography, but creates flat areas (lakes) where the original DEM had a depression.
- **Breaching** — carves a path from the pit to a lower neighbor. Better preserves natural landforms but can introduce artifacts in flat regions.

**Default: priority-flood with epsilon (Barnes 2014).** Fills pits but adds a tiny per-step monotonic increment so downstream routing through the filled region is well-defined and no genuine flats are introduced. Parallelizable: tile the DEM, fill tiles independently, reconcile tile boundaries.

**Rationale for fill over breach as default.** Breaching preserves natural landforms better in pure hydrological terms, but a breach far from a natural outlet carves a deep, narrow channel through high terrain — highly visible in a TIN and visually unacceptable. Filling produces small plateaus instead, which are forgivable in visualization. The hybrid Lindsay (2016) "breach if short and shallow, fill otherwise" approach is the right escape hatch for users who care more about hydrological realism than visual fidelity, exposed as a configuration flag.

**Real lakes are handled by the caller, not auto-detected.** Lake polygons must come from external data (land-cover datasets, OSM, national topo). Detection from post-fill flats would have too many false positives (agricultural flats, river plains, real plateaus). Given lake polygons:

- Pixels inside a lake polygon are excluded from D8 flow direction (no internal routing).
- The lake outline is treated as a spillway: flow exits the lake at the lowest pixel on the boundary.
- The lake polygon becomes a CDT constraint polygon in its own right.

## Flow direction

- **D8** (O'Callaghan & Mark 1984): each pixel drains to its single steepest-descent neighbor among the eight surrounding cells. Output is one of 8 directions per pixel.
- **D-infinity** (Tarboton 1997): drainage can split between two adjacent downhill neighbors weighted by aspect angle. More physically realistic, smoother accumulation maps.

**Recommendation:** D8 for this pipeline. D-infinity gives better accumulation rasters but introduces fractional flow that complicates stream vectorization — streams are no longer single-pixel-width paths. Since downstream we need clean polylines as CDT constraints, D8's discrete one-cell-wide streams vectorize cleanly without thinning.

Flow direction is per-pixel — embarrassingly parallel.

### DEM edge convention

The DEM has an arbitrary rectangular boundary, but water flowing toward the edge has nowhere to go in the data. Convention: **boundary pixels are virtual outlets** — effectively `−∞` elevation outside the raster, so any flow reaching the edge terminates there and routing always converges.

This is necessary for correctness (no infinite loops, well-defined sinks) but means catchments touching the DEM boundary are **incomplete** — the true catchment likely extends beyond the data. The catchment delineation step must therefore flag any computed catchment that touches the boundary so the caller knows the result is truncated. The same flag applies to lake polygons that intersect the boundary.

## Flow accumulation

For each pixel, count the number of upstream pixels (or total upstream area) draining through it.

Naive serial: topological sort by elevation, then per-pixel accumulate from upstream. Hard to parallelize directly because of the DAG dependency.

**Recommendation:** the Barnes-Lehman-Mulla parallel algorithm — tile the raster, accumulate within tiles, propagate boundary flows iteratively until convergence. Scales well to many cores.

Output: one float (or int) per pixel.

## Catchment delineation

Given a seed (or pour) pixel, walk the inverse flow graph: collect every pixel whose D8 outflow path eventually reaches the seed.

Algorithm: BFS from the seed in the upstream direction (i.e., for each pixel in the current frontier, find all neighbors whose D8 direction points to it; add them to the next frontier). Mark every visited pixel as in-catchment.

Parallelism: level-synchronous BFS — frontier expansion is parallel-for over the current frontier.

Output: a boolean mask raster (in/out of catchment). For the CDT constraint, polygonize the mask boundary using marching squares or contour tracing, then simplify (Visvalingam-Whyatt at the noding-step resolution).

**Boundary-touch flag.** If any in-catchment pixel sits on the DEM boundary, set a flag on the output. The catchment is truncated to the data extent and the true upstream area extends beyond the DEM. Callers can decide whether to error, warn, or proceed.

## Stream extraction

Streams are pixels with above-threshold flow accumulation. Two parameters:

- **Threshold** for "is a stream" — typically a constant upstream area (e.g., 1 km² for minor streams, more for rivers). Calibratable.
- **Threshold for "is a river" vs. "is a creek"** — separate, higher cutoff applied to stream order (Strahler).

### Vectorization

D8 streams are one pixel wide. To produce polylines:

1. Identify stream pixels (accumulation > threshold).
2. Identify confluence pixels (stream pixels with 2+ upstream stream neighbors).
3. Trace each stream segment from a headwater or confluence down to the next confluence or outlet, walking the D8 direction.
4. Each trace becomes a polyline; confluences are shared endpoints.

The resulting graph is a tree rooted at the outlet. Polylines correspond to tree edges.

### Stream ordering

**Strahler ordering** assigns order 1 to headwater segments; when two segments of order `k` meet, the downstream segment is order `k+1`; when segments of different orders meet, the downstream segment inherits the higher order.

Output: stream order per polyline. The caller picks a cutoff to separate creeks (low order) from rivers (high order) and decides which to emit as constraints.

## Output for the noding step

- **Catchment polygon** — outer boundary of the catchment mask, simplified to the noding-grid resolution. Becomes the surrounding polygon constraint.
- **Internal polygons** (lakes, optionally pre-classified land cover) — fed separately, with the empty property set. Empty is the legal default, not a missing flag: an unclassified constraint is the normal state of the system. A lake that the caller *does* want classified may carry a property like any other chain; nothing about an area feature forbids it.
- **River/creek polylines** — each polyline carries the `river` property. Caller decides the Strahler-order cutoff for inclusion.

A chain's tag is a **set** of properties, not one boolean: the same polyline may
be both `river` and, say, `road` where a track runs along the bank, and the
noder merges contributors by union. Which bit means `river` is a Python-side
vocabulary (`src_python/tin_engine/features.py`); the C++ core merges opaque
bits and names no feature. See `docs/increments/07-edge-properties.md`.

All these go into the noder (see `parallel_refinement.md` § Constraint noding), which resolves any intersections (e.g. a creek crossing a forest boundary) before the CDT.

## Library choice

**RichDEM** (Richard Barnes, MIT) is the obvious reference and possible dependency. It implements:

- Parallel priority-flood pit filling (Barnes 2014).
- D8 and D-infinity flow direction.
- Parallel flow accumulation.
- Watershed delineation.

It's the standard implementation of these algorithms, MIT-licensed, and actively maintained. Either depend on it directly or use it as a reference for in-tree implementations.

Other libraries (not viable on license grounds):

- **WhiteboxTools** — GPL, Rust.
- **TauDEM** — for D-infinity workflows; license needs checking.
- **GRASS GIS** — GPL, heavyweight.
- **PySheds** — BSD-3, Python; useful for prototyping but not as production C++ dependency.

## Parallelization summary

| Step                       | Parallel pattern                                  |
|----------------------------|---------------------------------------------------|
| Pit filling                | Tiled priority-flood (Barnes 2014)                |
| Flow direction             | Per-pixel parallel-for                            |
| Flow accumulation          | Tiled accumulate + iterative boundary propagation |
| Catchment BFS              | Level-synchronous frontier expansion              |
| Stream mask                | Per-pixel threshold                               |
| Stream vectorization       | Parallel per-trace (independent traces)           |
| Strahler ordering          | Post-order traversal of stream tree; parallel via subtree partitioning |
| Polygonize catchment       | Marching squares (per-cell parallel)              |

Memory bandwidth is the typical bottleneck, not compute. Tile sizes should be chosen to fit per-core L2/L3 cache.

## Settled defaults

- **Pit removal:** priority-flood with epsilon. Hybrid Lindsay breach-or-fill exposed as a flag for hydrology-realism users.
- **Lakes:** caller-supplied polygons; no auto-detection.
- **DEM edge:** virtual outlet at boundary; flag any catchment or lake that touches it as truncated.
- **Flow direction:** D8 (chosen for clean stream vectorization).

## Runtime parameters (not design decisions)

- **Stream accumulation threshold** — minimum upstream area to qualify as a stream pixel.
- **River-vs-creek Strahler order cutoff** — caller picks the order at which a stream is classified as a river for `river` tagging. The cutoff survives the widening to property sets unchanged: deciding when a creek becomes a river is a real caller decision about the *data*, not about the carrier.
- **Hybrid pit-removal flag** — opt in to Lindsay (2016) breach-or-fill instead of pure fill.
