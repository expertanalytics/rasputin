"""File formats. Decoding and encoding both live here.

`project_structure.md` marked this package as the home of "all file decoding";
increment 10 makes it the first *encoding* module too, and widens that line.
The rule the package keeps either way is the one that matters: `io/` knows
formats and knows nothing about `_core`, so everything in it is testable with
no compiled extension in the process and no filesystem.

Nothing here opens a file. A path belongs to `cli.py`, which is the only module
that has one (`06-cdt-viewer.md`, "no file is written below `cli.py`").
"""

from __future__ import annotations

from .ply import write_ply
from .vtk_legacy import write_vtk

__all__ = ["write_ply", "write_vtk"]
