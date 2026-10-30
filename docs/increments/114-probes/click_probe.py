"""Increment 114: what click says when a five-value --bbox (CRS last) gets four numbers.

    python docs/increments/114-probes/click_probe.py
"""

from typing import Annotated

import typer
from typer.testing import CliRunner

app = typer.Typer()
@app.command()
def mesh(out: Annotated[str, typer.Option("--out")],
         bbox: Annotated[tuple[float, float, float, float, str] | None,
                         typer.Option("--bbox", metavar="XMIN YMIN XMAX YMAX CRS")] = None) -> None:
    print("ok", bbox)
@app.command()
def other() -> None: pass
r = CliRunner()
for argv in (["mesh", "--bbox", "1", "2", "3", "4", "--out", "a.vtk"],
             ["mesh", "--out", "a.vtk", "--bbox", "1", "2", "3", "4"],
             ["mesh", "--out", "a.vtk", "--bbox", "1", "2", "3", "4", "EPSG:4326"]):
    res = r.invoke(app, argv, env={"TERM": "dumb", "NO_COLOR": "1", "TYPER_USE_RICH": "0"})
    print(argv, res.exit_code, [line for line in res.output.splitlines() if "Error" in line])
