"""Master's CORINE legends and palette, snapshotted before they became data.

Increment 41 (`docs/increments/41-class-maps-as-data.md`, section 7): the
expected values that come from master are read from `tests/python/data/`,
written at `f09b9863` (the production code of `ab0ca7eb` is the same) from
`feature_input`'s built-in table of legends, `feature_input._mask`, `feature_input._code` and
`palettes.CORINE_NATURAL`, and from `rasputin palette corine`'s stdout. They
are master's objects, not the shipped files, so a slip in transcribing a
class, a label, a colour or a water flag into the files is caught here.

`legends_f09b9863.json`:
- `palette`: `[code, label, colour]` in `CORINE_NATURAL`'s key order, 0 first.
- `maps[name]`: `name`, `attribute`, `codes` (the code system), `notice`;
  `masks`: `[value, mask]` for T2's values, `null` for a dropped value;
  `every_class`: each level-3 code's mask (as text), `null` when dropped;
  `class_codes`: `[value, code]` from `_code`, `"refused"` where it raised;
  `opened`: `[value, features]`, what `open_features` made of one polygon
  carrying that value, as `[code, mask]` per feature (`[]`: dropped), or
  `"refused"` where it raised `FeatureError`.

`palette_corine_f09b9863.json`: the exact bytes `rasputin palette corine`
wrote to stdout.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

DATA = Path(__file__).parent / "data"
SNAPSHOT: dict[str, Any] = json.loads((DATA / "legends_f09b9863.json").read_text("utf-8"))
PRESET_BYTES = (DATA / "palette_corine_f09b9863.json").read_bytes()
PRESET_NAME: str = SNAPSHOT["preset_title"]

#: Master's `palettes.CORINE_NATURAL`, in its key order: code -> (label, colour).
CORINE_NATURAL: dict[int, tuple[str, str]] = {
    int(code): (label, colour) for code, label, colour in SNAPSHOT["palette"]
}
CORINE_MAPS = ("corine", "corine-water", "clc18_kode")
