"""Scenes, readers and oracles for increment 23c's end-to-end suites.

Test support only; nothing in `src_python` uses it, and it imports nothing
23c adds, so a missing module fails the tests that reach it, not the
collection of this helper.

`docs/increments/23-basin-scale.md`, "Tests @tester can write red", 23c. The
oracles read the **piece files** a run writes, never the seam records or the
index's own counts (DC2: "The oracle reads the piece files"), and share no
code with `tin_engine`:

- conformity (K4): every seam edge, as its two end points `(x, y, z)` bit for
  bit, is in exactly two piece files;
- the union (DC3): pieces are glued on exact `(x, y)`; every edge is in one
  or two triangles, one only on the domain's boundary; the areas sum to the
  run's noded outline's (taken by a spy on `node`), which lies within the
  noder's snap of the domain; no two triangles overlap;
- constrained Delaunay (DC4, and tester.md section 3D's first oracle): an
  exact incircle (Python `fractions` over the doubles as written) on every
  interior edge of the union that is not a constraint of any piece file,
  in world coordinates as written, which is the frame refine's oracle
  uses up to a translation (`delaunay_oracle` in
  `tests/cpp/property/prop_refinement_constraint_feet.cpp` scales the lattice
  fraction back by the spacing);
- the tolerance (DC5, section 3D's second oracle): every valid DEM node
  inside the domain against the plane of a union triangle holding it,
  recomputed from the written vertices (`geographic_fixtures.source_check`,
  with a 1e-9 slack).

PINNED by these readers, where the design is silent (handback, choices 2-5):
pieces go to the one directory `<out name>.pieces` beside `--out` (found here
by its suffix); `index.json` there holds `crs`, `tolerance`, a `partition`
object (`pieces`, `memory_budget`, `bytes_per_node`, `nx`, `ny`, `dx`, `dy`,
`cols`, `rows`) and a `pieces` list whose entries carry `id` (`[j, i, k]`),
`file` (relative to the directory), `sha256` and `counts`; a piece file names
the seam bit `seam` in its own vocabulary.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry import Polygon

from geographic_fixtures import SourceCheck, source_check
from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff
from test_cli_mesh_dem import invoke
from vtkread import VtkFile, lines_as_array, polygons_as_array, read_vtk

Ring = list[tuple[float, float]]
Key = tuple[float, float]
XYZ = tuple[float, float, float]
PieceId = tuple[int, int, int]

#: The scene's lattice: square 10 m cells, 41 x 41 nodes, micro_tiff's tie point.
H = 10.0
N = 41
UTM33 = "urn:ogc:def:crs:EPSG::25833"
SNAP = 1e-3


def rel(x: float, y: float) -> tuple[float, float]:
    return (TIE_X + x, TIE_Y + y)


#: Off-node in both axes, whole decimetres (so on the noder's 1 mm grid).
SQUARE: Ring = [rel(32.3, -28.7), rel(27.1, -368.9), rel(366.7, -372.3), rel(371.9, -33.1)]


def hills(n: int = N, seed: int = 23, amplitude: float = 20.0, noise: float = 3.0) -> np.ndarray:
    """Smooth hills with noise on top, so a 1 m tolerance inserts many nodes."""
    r, c = np.indices((n, n), dtype=np.float64)
    smooth = amplitude * np.sin(c / 7.0) * np.cos(r / 9.0) + 0.4 * c
    rough = np.random.default_rng(seed).uniform(0.0, noise, (n, n))
    return np.asarray(50.0 + smooth + rough, dtype=np.float32)


def write_dem(path: Path, array: np.ndarray, h: float = H) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(micro_tiff(array, scale=(h, h, 0.0)).getvalue())
    return path


def write_domain(path: Path, outer: Ring, holes: tuple[Ring, ...] = ()) -> Path:
    doc: dict[str, Any] = {
        "type": "Polygon",
        "coordinates": [[*r, r[0]] for r in (outer, *holes)],
        "crs": {"type": "name", "properties": {"name": UTM33}},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc))
    return path


def node_xy(array: np.ndarray, h: float = H) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Every node's world x, y and z, in row-major order."""
    rows, cols = np.indices(array.shape)
    x = TIE_X + cols.ravel() * h
    y = TIE_Y - rows.ravel() * h
    return x.astype(np.float64), y.astype(np.float64), array.astype(np.float64).ravel()


# ------------------------------------------------------------------ the run


@dataclass
class Piece:
    id: PieceId
    path: Path
    vtk: VtkFile
    points: npt.NDArray[np.float64]
    triangles: npt.NDArray[np.int64]
    edges: npt.NDArray[np.int64]
    masks: npt.NDArray[np.int64]
    seam: int  # the seam's bit value in this file's vocabulary


@dataclass
class Cut:
    code: int
    output: str
    out: Path
    dir: Path | None
    index: dict[str, Any] = field(default_factory=dict)
    pieces: list[Piece] = field(default_factory=list)
    noded: Polygon | None = None  # the run's own noded outline, from its `node` call


def bit_named(vtk: VtkFile, name: str) -> int:
    """The bit value a file's own vocabulary gives `name`."""
    bits = [int(b) for b in np.asarray(vtk.field_data["feature_bits"].values)]
    names = list(vtk.field_data["feature_names"].values)
    assert name in names, f"{name!r} is not in the file's vocabulary {names}"
    return 1 << bits[names.index(name)]


def read_piece(piece_id: PieceId, path: Path) -> Piece:
    vtk = read_vtk(path.read_bytes())
    edges = lines_as_array(vtk) if vtk.lines else np.zeros((0, 2), dtype=np.int64)
    masks = np.asarray(vtk.cell_array("feature_mask").values, dtype=np.int64)[: len(edges)]
    return Piece(
        id=piece_id,
        path=path,
        vtk=vtk,
        points=np.asarray(vtk.points, dtype=np.float64),
        triangles=polygons_as_array(vtk),
        edges=edges,
        masks=masks,
        seam=bit_named(vtk, "seam"),
    )


def pieces_dir(out: Path) -> Path | None:
    found = [p for p in out.parent.iterdir() if p.is_dir() and p.name.endswith(".pieces")]
    assert len(found) <= 1, found
    return found[0] if found else None


@contextmanager
def node_spy() -> Iterator[list[Any]]:
    """Every `node` outcome the run makes, wherever `tin_engine` holds `node`
    (the spy K9's test uses); restored on exit."""
    from tin_engine import _core

    original = _core.node
    seen: list[Any] = []

    def spy(*args: Any, **kwargs: Any) -> Any:
        outcome = original(*args, **kwargs)
        seen.append(outcome)
        return outcome

    held = [m for k, m in list(sys.modules.items()) if k.startswith("tin_engine")]
    patched = [m for m in held if getattr(m, "node", None) is original]
    for module in patched:
        module.node = spy
    try:
        yield seen
    finally:
        for module in patched:
            module.node = original


def noded_outline(outcome: Any) -> Polygon:
    """The noded domain: the outer chain as the shell, the hole chains as holes."""
    from tin_engine._core import ChainRole

    pslg = outcome.pslg
    xy = np.asarray(pslg.vertices, dtype=np.float64)
    rings: dict[str, list[npt.NDArray[np.float64]]] = {"outer": [], "hole": []}
    for c, chain in enumerate(pslg.chains):
        if chain.role in (ChainRole.Outer, ChainRole.Hole):
            ring = xy[np.asarray(pslg.indices_of(c), dtype=np.int64)]
            rings["outer" if chain.role == ChainRole.Outer else "hole"].append(ring)
    assert len(rings["outer"]) == 1, "not one outer ring"
    return Polygon(rings["outer"][0], rings["hole"])


def run(where: Path, *args: str, out_name: str = "x.vtk") -> Cut:
    """`rasputin mesh *args --out where/out_name`, in a directory of its own."""
    where.mkdir(parents=True, exist_ok=True)
    out = where / out_name
    with node_spy() as seen:
        code, output = invoke(*args, "--out", str(out))
    done = Cut(code=code, output=output, out=out, dir=pieces_dir(out))
    if len(seen) == 1 and seen[0].ok():
        done.noded = noded_outline(seen[0])
    return done


def cut(where: Path, *args: str, out_name: str = "x.vtk") -> Cut:
    """A run that must succeed and be cut: its index and every piece file read."""
    done = run(where, *args, out_name=out_name)
    assert done.code == 0, done.output
    assert done.dir is not None, f"no <out>.pieces directory beside {done.out}"
    done.index = json.loads((done.dir / "index.json").read_text())
    for entry in done.index["pieces"]:
        piece_id = tuple(int(v) for v in entry["id"])
        assert len(piece_id) == 3, entry
        path = done.dir / entry["file"]
        if path.suffix == ".vtk":
            done.pieces.append(read_piece(piece_id, path))  # type: ignore[arg-type]
    return done


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------ K4, DC2


def xyz(p: npt.NDArray[np.float64]) -> XYZ:
    return (float(p[0]), float(p[1]), float(p[2]))


def seam_edges(piece: Piece) -> list[tuple[XYZ, XYZ]]:
    out = []
    for (a, b), m in zip(piece.edges, piece.masks, strict=True):
        if m & piece.seam:
            pa, pb = xyz(piece.points[a]), xyz(piece.points[b])
            out.append((min(pa, pb), max(pa, pb)))
    return out


def conformity_problems(pieces: list[Piece]) -> list[str]:
    """Every seam edge, bit for bit with z, in exactly two distinct pieces."""
    holders: dict[tuple[XYZ, XYZ], list[PieceId]] = defaultdict(list)
    for piece in pieces:
        for edge in seam_edges(piece):
            holders[edge].append(piece.id)
    problems = [
        f"seam edge {edge} held by {ids}"
        for edge, ids in holders.items()
        if len(ids) != 2 or ids[0] == ids[1]
    ]
    if not holders:
        problems.append("no seam edge in any piece file")
    return problems


def off_lattice_seams(pieces: list[Piece], h: float = H) -> list[tuple[XYZ, XYZ]]:
    """Seam edges not on a lattice column or row (x or y at a whole node)."""
    bad = []
    for piece in pieces:
        for a, b in seam_edges(piece):
            vertical = a[0] == b[0] and float((a[0] - TIE_X) / h).is_integer()
            horizontal = a[1] == b[1] and float((TIE_Y - a[1]) / h).is_integer()
            if not (vertical or horizontal):
                bad.append((a, b))
    return bad


def seam_lines(pieces: list[Piece]) -> tuple[set[float], set[float]]:
    """The x of every vertical seam edge and the y of every horizontal one."""
    xs, ys = set(), set()
    for piece in pieces:
        for a, b in seam_edges(piece):
            if a[0] == b[0]:
                xs.add(a[0])
            if a[1] == b[1]:
                ys.add(a[1])
    return xs, ys


# ------------------------------------------------------------------ the union


@dataclass
class Union:
    points: npt.NDArray[np.float64]  # (N, 3)
    triangles: npt.NDArray[np.int64]  # (T, 3)
    constrained: set[tuple[int, int]]  # every constraint edge of any piece file
    seam_only: set[tuple[int, int]]  # those whose only bit is the seam's
    z_conflicts: list[Key]


def union(pieces: list[Piece]) -> Union:
    """The piece files glued on exact `(x, y)`; the stitcher is 23d's."""
    index: dict[Key, int] = {}
    points: list[XYZ] = []
    conflicts: list[Key] = []

    def vertex(p: npt.NDArray[np.float64]) -> int:
        key = (float(p[0]), float(p[1]))
        if key not in index:
            index[key] = len(points)
            points.append(xyz(p))
        elif np.float64(points[index[key]][2]).tobytes() != np.float64(p[2]).tobytes():
            conflicts.append(key)
        return index[key]

    triangles: list[tuple[int, int, int]] = []
    constrained: set[tuple[int, int]] = set()
    seam_only: set[tuple[int, int]] = set()
    for piece in pieces:
        ids = [vertex(p) for p in piece.points]
        triangles += [(ids[a], ids[b], ids[c]) for a, b, c in piece.triangles]
        for (a, b), m in zip(piece.edges, piece.masks, strict=True):
            e = (min(ids[a], ids[b]), max(ids[a], ids[b]))
            constrained.add(e)
            if m == piece.seam:
                seam_only.add(e)
    return Union(
        points=np.array(points, dtype=np.float64).reshape(-1, 3),
        triangles=np.array(triangles, dtype=np.int64).reshape(-1, 3),
        constrained=constrained,
        seam_only=seam_only,
        z_conflicts=conflicts,
    )


def sides(triangles: npt.NDArray[np.int64]) -> dict[tuple[int, int], list[tuple[int, int]]]:
    """Each undirected edge: (triangle, apex) for every triangle holding it."""
    out: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for t, tri in enumerate(triangles.tolist()):
        for e in range(3):
            a, b, c = tri[e], tri[(e + 1) % 3], tri[(e + 2) % 3]
            out[(min(a, b), max(a, b))].append((t, c))
    return out


#: How far the noder may move a vertex, a crossing included (05-noder.md, 14(b)).
SNAP_MOVE = SNAP * math.sqrt(2) / 2


def noded_problems(noded: Polygon, domain: Polygon) -> list[str]:
    """The noded outline against the domain ("23c-2 stopped at dbbe1c6", 1):
    every vertex within the snap's reach of the domain's boundary, and the
    areas within that reach times the domain's perimeter."""
    problems = []
    boundary = domain.boundary
    rings = [noded.exterior, *noded.interiors]
    for ring in rings:
        for x, y in shapely.get_coordinates(ring):
            d = boundary.distance(shapely.Point(x, y))
            if d > SNAP_MOVE + 1e-9:
                problems.append(f"noded vertex ({x!r}, {y!r}) {d} from the domain")
    if abs(noded.area - domain.area) > SNAP_MOVE * domain.length:
        problems.append(f"noded area {noded.area!r} against the domain's {domain.area!r}")
    return problems


def validity_problems(u: Union, domain: Polygon, noded: Polygon | None) -> list[str]:
    """DC3: a valid union of the piece files over `domain`: its area is the
    run's noded outline's to 1e-9, and that outline is the domain's to the
    noder's snap (`noded_problems`)."""
    if noded is None:
        return ["no noded outline: the run's node call was not seen"]
    problems = noded_problems(noded, domain)
    boundary = domain.boundary
    for (a, b), held in sides(u.triangles).items():
        if len(held) > 2:
            problems.append(f"edge {a}-{b} in {len(held)} triangles")
        elif len(held) == 1:
            mid = shapely.Point(*((u.points[a, :2] + u.points[b, :2]) / 2))
            if boundary.distance(mid) > 2 * SNAP:
                problems.append(f"edge {a}-{b} in one triangle, inside the domain")
    xy = u.points[:, :2]
    a, b, c = xy[u.triangles[:, 0]], xy[u.triangles[:, 1]], xy[u.triangles[:, 2]]
    twice = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    if (twice <= 0).any():
        problems.append(f"{int((twice <= 0).sum())} triangles not counter-clockwise")
    area = float(np.sum(twice) / 2)
    if abs(area - noded.area) > 1e-9 * noded.area:
        problems.append(f"area {area!r} against the noded outline's {noded.area!r}")
    return problems


def overlapping_pairs(u: Union) -> list[tuple[int, int]]:
    """Brute force over pairs whose boxes meet: interiors that intersect."""
    polys = [Polygon(u.points[t, :2]) for t in u.triangles]
    tree = shapely.STRtree(polys)
    left, right = tree.query(polys, predicate="intersects")
    bad = []
    for i, j in zip(left.tolist(), right.tolist(), strict=True):
        if i < j:
            shared = polys[i].intersection(polys[j]).area
            if shared > 1e-9 * min(polys[i].area, polys[j].area):
                bad.append((i, j))
    return bad


# ------------------------------------------------------------------ section 3D


def _f(v: float) -> Fraction:
    return Fraction(float(v))


def orient(a: Any, b: Any, c: Any) -> int:
    det = (_f(b[0]) - _f(a[0])) * (_f(c[1]) - _f(a[1])) - (_f(b[1]) - _f(a[1])) * (
        _f(c[0]) - _f(a[0])
    )
    return (det > 0) - (det < 0)


def incircle(a: Any, b: Any, c: Any, d: Any) -> int:
    """The sign of the incircle determinant of counter-clockwise `a, b, c`, exactly."""
    rows = []
    for p in (a, b, c):
        dx, dy = _f(p[0]) - _f(d[0]), _f(p[1]) - _f(d[1])
        rows.append((dx, dy, dx * dx + dy * dy))
    (a1, a2, a3), (b1, b2, b3), (c1, c2, c3) = rows
    det = a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)
    return (det > 0) - (det < 0)


def delaunay_violations(
    xy: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.int64],
    constrained: set[tuple[int, int]],
) -> list[tuple[tuple[int, int], int]]:
    """Every interior edge not in `constrained` with an apex strictly inside
    the other triangle's circumcircle, decided exactly."""
    bad = []
    for edge, held in sides(triangles).items():
        if len(held) != 2 or edge in constrained:
            continue
        for (t, _), (_, apex) in (held, held[::-1]):
            a, b, c = (xy[i] for i in triangles[t])
            if orient(a, b, c) < 0:
                a, b = b, a
            if incircle(a, b, c, xy[apex]) > 0:
                bad.append((edge, apex))
    return bad


def tolerance_check(
    u: Union, array: np.ndarray, domain: Polygon, tolerance: float, h: float = H
) -> SourceCheck:
    """Every valid DEM node strictly inside `domain` against the union, from
    the written vertices, with refine's 1e-9 slack."""
    x, y, z = node_xy(array, h)
    return source_check(
        u.points,
        u.triangles,
        np.column_stack([x, y]),
        z,
        domain,
        tolerance,
        h,
        slack=1e-9,
    )


def triangle_set(points: npt.NDArray[np.float64], triangles: Any) -> set[frozenset[Key]]:
    return {
        frozenset((float(points[i, 0]), float(points[i, 1])) for i in tri)
        for tri in np.asarray(triangles).tolist()
    }


def problems(done: Cut, array: np.ndarray, domain: Polygon, tolerance: float) -> list[str]:
    """DC2-DC5 over one cut run (both section 3D oracles), as what is wrong."""
    found = conformity_problems(done.pieces)
    found += [f"seam edge off the lattice: {e}" for e in off_lattice_seams(done.pieces)]
    u = union(done.pieces)
    found += [f"two z at {k}" for k in u.z_conflicts]
    found += validity_problems(u, domain, done.noded)
    found += [f"triangles {p} overlap" for p in overlapping_pairs(u)]
    violations = delaunay_violations(u.points[:, :2], u.triangles, u.constrained)
    found += [f"not constrained Delaunay: {v}" for v in violations]
    check = tolerance_check(u, array, domain, tolerance)
    if check.over or check.outside_mesh:
        found.append(f"tolerance: {check}")
    return found


def located_over(
    u: Union, array: np.ndarray, valid: np.ndarray, tolerance: float, h: float = H
) -> tuple[int, int]:
    """Section 3D's tolerance oracle where NoData carves the mesh: every valid
    DEM node that some union triangle holds (closed, 1e-12 slack) against
    that triangle's plane, from the written vertices; every output vertex has
    data, so every triangle has three valid vertices. Returns (located, over)."""
    x, y, z = node_xy(array, h)
    keep = valid.ravel()
    px, py, pz = x[keep], y[keep], z[keep]
    origin = u.points[:, :2].min(axis=0)
    v = u.points[:, :2] - origin
    t = u.triangles
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    za, zb, zc = u.points[t[:, 0], 2], u.points[t[:, 1], 2], u.points[t[:, 2], 2]
    two_a = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    located = over = 0
    for lo in range(0, px.size, 128):
        qx = px[lo : lo + 128, None] - origin[0]
        qy = py[lo : lo + 128, None] - origin[1]

        def weight(p: Any, q: Any, qx: Any = qx, qy: Any = qy) -> Any:
            return (
                (q[:, 0] - p[:, 0]) * (qy - p[:, 1]) - (q[:, 1] - p[:, 1]) * (qx - p[:, 0])
            ) / two_a

        wa, wb, wc = weight(b, c), weight(c, a), weight(a, b)
        holds = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
        error = np.abs(wa * za + wb * zb + wc * zc - pz[lo : lo + 128, None])
        worst = np.where(holds, error, -np.inf).max(axis=1)
        located += int(np.isfinite(worst).sum())
        over += int((worst > tolerance + 1e-9).sum())
    return located, over
