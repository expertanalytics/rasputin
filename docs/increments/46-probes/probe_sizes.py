# ruff: noqa: E501 -- a one-off design probe; its formulas and report lines are long
"""Increment 46 design probe: mesh sizes, grid convergence and the EnSight animation size.

Run after probe_spring.py, with the same WORKDIR (it reuses synth.tif):
    .venv/bin/python docs/increments/46-probes/probe_sizes.py WORKDIR > docs/increments/46-probes/probe_sizes.out
Writes in WORKDIR: flat.tif, mesh_1m.vtk, mesh_0p5m.vtk, mesh_flat_1m.vtk and an
EnSight set ens/ of 92 frames, read back by ParaView's pvpython if installed.
"""

import pathlib
import subprocess
import sys
import time

import numpy as np
import tifffile
from pyproj import CRS, Proj, Transformer

work = pathlib.Path(sys.argv[1])
rasputin = pathlib.Path(sys.executable).parent / "rasputin"

# The terrain of probe_spring.py, without its noise; the noisy grid is checked equal to synth.tif.
n_side, step, x0, ytop = 1201, 10.0, 494000.0, 6806000.0
xs = (np.arange(n_side) * step) - 6000.0
ys = 6000.0 - (np.arange(n_side) * step)
X, Y = np.meshgrid(xs, ys)
r, phi = np.hypot(X, Y), np.arctan2(Y, X)
z = 300 + 1600 * np.exp(-((r / 2500) ** 2)) * (
    1 + 0.15 * (1 - np.exp(-((r / 1500) ** 2))) * np.cos(5 * phi)
)
z += 250 * np.exp(-((X - 3500) ** 2 + (Y + 3000) ** 2) / 900**2)
rng = np.random.default_rng(46)
k = np.hypot(*np.meshgrid(np.fft.fftfreq(n_side, step), np.fft.fftfreq(n_side, step)))
amp = np.zeros_like(k)
keep = (k > 0) & (k <= 1 / 40.0)
amp[keep] = k[keep] ** -1.5
noise = np.real(np.fft.ifft2(amp * np.exp(2j * np.pi * rng.random(k.shape))))
noisy = z + noise * (25.0 / noise.std())
assert np.array_equal(noisy.astype(np.float32), tifffile.imread(work / "synth.tif")), (
    "terrain differs from probe_spring.py"
)
geokeys = (1, 1, 0, 3, 1024, 0, 1, 1, 1025, 0, 1, 2, 3072, 0, 1, 25833)
tifffile.imwrite(
    work / "flat.tif",
    z.astype(np.float32),
    extratags=[
        (33550, "d", 3, (step, step, 0.0), False),
        (33922, "d", 6, (0, 0, 0, x0, ytop, 0), False),
        (34735, "H", len(geokeys), geokeys, False),
    ],
)


def mesh(tif, tol, out):
    t = time.perf_counter()
    subprocess.run(
        [
            str(rasputin),
            "mesh",
            "--dem",
            str(work / tif),
            "--tolerance",
            str(tol),
            "--out",
            str(work / out),
        ],
        check=True,
        capture_output=True,
    )
    wall = time.perf_counter() - t
    data = (work / out).read_bytes()
    i = data.index(b"POINTS ")
    nv = int(data[i : data.index(b"\n", i)].split()[1])
    p0 = data.index(b"\n", i) + 1
    pts = np.frombuffer(data, ">f8", nv * 3, p0).reshape(-1, 3)
    i = data.index(b"POLYGONS ", p0 + nv * 24)
    line = data[i : data.index(b"\n", i)].split()
    tri = np.frombuffer(data, ">i4", int(line[2]), data.index(b"\n", i) + 1).reshape(-1, 4)[:, 1:]
    a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    plan = 0.5 * np.abs(
        (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    )
    print(
        f"{out}: tolerance {tol} m, {len(tri)} triangles, {nv} vertices, mean plan area {plan.mean():.0f} m2, "
        f"file {len(data)} bytes, {wall:.2f} s wall incl. process start and write"
    )
    return pts, tri


for tif, tol, out in (
    ("synth.tif", 2, "mesh_2m.vtk"),
    ("flat.tif", 1, "mesh_flat_1m.vtk"),
    ("synth.tif", 0.5, "mesh_0p5m.vtk"),
):
    mesh(tif, tol, out)
pts, tri = mesh("synth.tif", 1, "mesh_1m.vtk")

proj = Proj(CRS.from_epsg(25833))
to_ll = Transformer.from_crs(25833, 4326, always_xy=True)
for x in (494000.0, 500000.0, 506000.0):
    lon, lat = to_ll.transform(x, 6800000.0)
    print(
        f"gamma at ({x:.0f}, 6800000): lat {lat:.3f} lon {lon:.3f}, {proj.get_factors(lon, lat).meridian_convergence:+.4f} deg"
    )

# EnSight Gold binary: geometry once, two per-element values per frame, 92 frames.
ens = work / "ens"
ens.mkdir(exist_ok=True)
for f in ens.iterdir():
    f.unlink()


def s80(text):
    return text.encode("ascii").ljust(80, b"\0")


t = time.perf_counter()
ox, oy = pts[:, 0].min(), pts[:, 1].min()
geo = [
    s80("C Binary"),
    s80("rasputin snow"),
    s80(f"origin {ox} {oy}"),
    s80("node id off"),
    s80("element id off"),
    s80("part"),
    np.int32(1).tobytes(),
    s80("mesh"),
    s80("coordinates"),
    np.int32(len(pts)).tobytes(),
    (pts[:, 0] - ox).astype("<f4").tobytes(),
    (pts[:, 1] - oy).astype("<f4").tobytes(),
    pts[:, 2].astype("<f4").tobytes(),
    s80("tria3"),
    np.int32(len(tri)).tobytes(),
    (tri + 1).astype("<i4").tobytes(),
]
(ens / "snow.geo").write_bytes(b"".join(geo))
frames = 92
for name in ("swe", "sca"):
    for q in range(frames):
        v = np.full(len(tri), float(q), dtype="<f4")
        (ens / f"{name}{q:04d}.var").write_bytes(
            b"".join([s80(name), s80("part"), np.int32(1).tobytes(), s80("tria3"), v.tobytes()])
        )
times = "\n".join(
    str(24 * q) for q in range(frames)
)  # one per line: VTK reads case lines of at most 256 chars
(ens / "snow.case").write_text(
    "FORMAT\ntype: ensight gold\nGEOMETRY\nmodel: snow.geo\nVARIABLE\n"
    "scalar per element: 1 swe swe****.var\nscalar per element: 1 sca sca****.var\n"
    f"TIME\ntime set: 1\nnumber of steps: {frames}\nfilename start number: 0\n"
    f"filename increment: 1\ntime values:\n{times}\n"
)
wall = time.perf_counter() - t
sizes = {f.name: f.stat().st_size for f in ens.iterdir()}
total = sum(sizes.values())
print(
    f"EnSight: {frames} frames x 2 values written in {wall:.2f} s; geometry {sizes['snow.geo']} bytes; "
    f"one value one frame {sizes['swe0000.var']} bytes (244 + 4 x {len(tri)} = {244 + 4 * len(tri)}); "
    f"total {total} bytes = {total / 1e6:.1f} MB = {total / 2**20:.1f} MiB"
)
vtk_series = (work / "mesh_1m.vtk").stat().st_size + 2 * 4 * len(tri) + 200
print(
    f".vtk series estimate: ({(work / 'mesh_1m.vtk').stat().st_size} + 2 x 4 x {len(tri)}) x {frames} = {vtk_series * frames / 1e9:.2f} GB"
)
pv = pathlib.Path("/Applications/ParaView-6.1.1.app/Contents/bin/pvpython")
if pv.exists():
    script = (
        "from paraview.simple import *\nr = OpenDataFile({!r})\nr.UpdatePipeline()\nts = r.TimestepValues\n"
        "r.UpdatePipeline(ts[50])\nfrom paraview import servermanager as sm\nd = sm.Fetch(r)\n"
        "b = d.GetBlock(0) if d.IsA('vtkMultiBlockDataSet') else d\n"
        "print('pvpython read back:', len(ts), 'time steps', ts[0], ts[1], ts[-1], 'cells', b.GetNumberOfCells(),"
        " 'swe at step 50', b.GetCellData().GetArray('swe').GetValue(0))\n"
    ).format(str(ens / "snow.case"))
    res = subprocess.run([str(pv), "-c", script], capture_output=True, text=True)
    print(res.stdout.strip() or res.stderr.strip()[-500:])
