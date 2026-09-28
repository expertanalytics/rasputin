"""Repair the legacy CORINE GML so that it is well-formed XML.

    git show 8e30af4:.circleci/rasputin_data/corine/0000_4326_corine2018_4e6064_GML.gml \
        > /tmp/gml-2020.gml
    python tests/fixtures/corine/repair_gml.py /tmp/gml-2020.gml \
        tests/fixtures/corine/0000_4326_corine2018_4e6064_GML.gml

Test support, not production code (`docs/increments/16b-terrain-polygons.md`,
"Test data": Ola, 2026-09-28, "Repear the broken file."). The file as
committed in 2020 (`8e30af4`) has two defects and no others:

1. `</gml:featureMember>` is missing after the `</ogr:sql_statement>` that
   closes feature `sql_statement.40965` (line 1097), so the next 45 members
   nest inside it;
2. the document ends without `</ogr:FeatureCollection>`.

This script adds exactly those two closing tags and changes nothing else. It
refuses any input other than the 2020 file (checked by its git blob hash, and
then by the defect itself), so running it on the repaired file, or on any
other file, is refused and writes nothing. The output is checked to be
well-formed and to hold 200 feature members before it is written.

Standard library only.
"""

from __future__ import annotations

import hashlib
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

#: The git blob hash of the file committed in 8e30af4.
BLOB_2020 = "97ddc98c4262b6b817e96cd91a23720b9b94d2a5"
FEATURE = '    <ogr:sql_statement fid="sql_statement.40965">\n'
CLOSE_FEATURE = "    </ogr:sql_statement>\n"
OPEN_MEMBER = "  <gml:featureMember>\n"
CLOSE_MEMBER = "  </gml:featureMember>\n"
CLOSE_ROOT = "</ogr:FeatureCollection>\n"
MEMBERS = 200


class RepairError(ValueError):
    """The input is not the 2020 file with exactly the known defect."""


def blob_hash(data: bytes) -> str:
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def repair(data: bytes) -> bytes:
    """The 2020 file's bytes -> the repaired file's bytes."""
    if blob_hash(data) != BLOB_2020:
        raise RepairError(f"not the 2020 file: git blob {blob_hash(data)}, expected {BLOB_2020}")
    lines = data.decode("utf-8").splitlines(keepends=True)
    if lines.count(FEATURE) != 1:
        raise RepairError("feature sql_statement.40965 is not in the file exactly once")
    at = lines.index(FEATURE)
    close = lines.index(CLOSE_FEATURE, at)
    if lines[close + 1] != OPEN_MEMBER:
        raise RepairError(f"line {close + 2} does not open the next member unclosed")
    if lines[-1] != CLOSE_MEMBER or CLOSE_ROOT in lines:
        raise RepairError("the document does not end with its root left open")
    repaired = "".join([*lines[: close + 1], CLOSE_MEMBER, *lines[close + 1 :], CLOSE_ROOT])
    root = ET.fromstring(repaired.encode("utf-8"))  # raises ParseError if not well-formed
    members = root.findall("{http://www.opengis.net/gml}featureMember")
    if len(members) != MEMBERS:
        raise RepairError(f"the repaired document has {len(members)} members, not {MEMBERS}")
    return repaired.encode("utf-8")


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    src, dst = Path(argv[1]), Path(argv[2])
    try:
        out = repair(src.read_bytes())
    except RepairError as err:
        print(f"{src}: {err}; nothing written", file=sys.stderr)
        return 1
    dst.write_bytes(out)
    print(f"{dst}: repaired, git blob {blob_hash(out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
