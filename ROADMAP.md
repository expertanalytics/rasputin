# Rasputin roadmap

An index, and only an index. What each increment *is* lives in its own record
under `docs/increments/`, and nothing is restated here: a duplicated detail is
how the previous version of this file rotted into a seven-year-stale description
of the CGAL project this one replaced.

**Authoring order is not ship order.** Increment 6 is designed sixth and ships
before `05b`, because the viewer depends on nothing the noder produces and 5b is
the increment that most needs a picture to be checked against. The reasoning is
in `docs/increments/06-cdt-viewer.md`, "The ordering ruling".

| # | What it is | Status | Record |
|---|---|---|---|
| 1 | Exact geometric predicates | shipped (`d0946bc`, `aeb667a`) | `docs/increments/01-predicates.md` |
| 2 | 2D geometry value types | shipped (`1a60ad4`) | `docs/increments/02-core-geometry.md` |
| 3 | The PSLG and its validator | shipped (`d984ff8`) | `docs/increments/03-pslg.md` |
| 4 | The constrained Delaunay triangulation | shipped (`f22ddd3`, `238e511`) | `docs/increments/04-cdt.md` |
| 5a | The noder's numeric floor: snap grid and hot pixels | shipped (`e090909`) | `docs/increments/05-noder.md` |
| 5b | The noder proper: segment splitting and the noded PSLG | not started; record not yet written | `docs/increments/05b-noder-driver.md` |
| 6a | The pybind11 CDT surface | shipped (`311459b`) | `docs/increments/06-cdt-viewer.md` |
| 6b-i | `build_scene`: the mesh-and-PSLG-to-geometry mapping | shipped (`ab4681c`) | `docs/increments/06-cdt-viewer.md` |
| 6b-ii | The SVG renderer, the fixture gallery and `rasputin draw` | in flight | `docs/increments/06-cdt-viewer.md` |

Beyond the current sequence, and deliberately undated because none of it has a
design record yet: refinement and its quality policy
(`parallel_refinement.md`), raster ingestion and the GeoTIFF reader
(`project_structure.md`'s `io/` and `raster` sections), and catchment
extraction (`auto_catchments.md`).

`docs/increments/README.md` is the protocol those records follow, and it — not
this table — is the file that states how work gets done.
