"""Assert that detria.hpp stays inside exactly one translation unit.

The vendored detria header (lib/detria, MIT by election -- see its README) is
4500 lines, pulls in <iostream>, <sstream> and <csignal>, and its Debug
assertions call std::raise(SIGTRAP). It is used for two functions: the exact
`orient2d` and `incircle` predicates behind `terrain::pred::DetriaExact`.

The project's claim that the exact backend is swappable rests entirely on that
header not leaking. If it appears in a header under include/, every consumer of
the predicates module compiles against it, the swap stops being a one-line
change in default_kernel.hpp, and a signal-raising assertion becomes reachable
from translation units that never meant to ask a geometric question.

CMake already makes the rule structural -- lib/detria is PRIVATE on the
terrain_predicates target, so the include would not resolve elsewhere -- and the
unit suite carries an `#ifdef DETRIA_HPP_INCLUDED -> #error` guard. This script
is the third line of defence, and the only one that survives someone "fixing"
the build by adding lib/detria to a broader include path: a structural guard
that a build edit can dissolve is not a guard, and an unenforced rule erodes.

Like tools/check_prohibited_deps.py, this parses #include directives rather than
grepping prose: lib/detria/README.md, CLAUDE.md and the headers themselves all
name detria.hpp precisely in order to explain the rule, and a naive grep would
flag those and get switched off.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# The only translation units permitted to include detria.hpp, relative to ROOT.
# One per backend: the exact predicates, and the CDT wrapper.
PERMITTED = {"src/predicates/detria_exact.cpp", "src/cdt/detria_backend.cpp"}

# Trees that must stay clean apart from PERMITTED. include/ is called out
# separately below because zero is the required count there, not "at most the
# permitted ones" -- no header may include it, ever.
SOURCE_DIRS = ["include", "src", "src_python", "bindings", "tests", "tools"]

CXX_SUFFIXES = {".cpp", ".hpp", ".h", ".cc", ".cxx"}

INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]')

# Matches `detria.hpp` however it is spelled on the include path.
DETRIA_RE = re.compile(r"(^|/)detria\.hpp$", re.IGNORECASE)


def detria_includes(path: Path) -> list[tuple[int, str]]:
    """Every (line number, included path) naming detria.hpp in this file."""
    found = []
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        m = INCLUDE_RE.match(line)
        if m and DETRIA_RE.search(m.group(1)):
            found.append((lineno, m.group(1)))
    return found


def main() -> int:
    findings: list[str] = []
    scanned = 0
    including: set[str] = set()

    for directory in SOURCE_DIRS:
        base = ROOT / directory
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file() or path.suffix not in CXX_SUFFIXES:
                continue
            scanned += 1
            rel = path.relative_to(ROOT).as_posix()
            for lineno, spelling in detria_includes(path):
                including.add(rel)
                if rel.startswith("include/"):
                    findings.append(
                        f"{rel}:{lineno}: includes {spelling} -- no header under include/ "
                        f"may include detria.hpp; declare in the header and define in "
                        f"src/predicates/detria_exact.cpp"
                    )
                elif rel not in PERMITTED:
                    findings.append(
                        f"{rel}:{lineno}: includes {spelling} -- only "
                        f"{', '.join(sorted(PERMITTED))} may do so"
                    )

    # The rule has a floor as well as a ceiling: if the permitted translation
    # unit stops including detria.hpp, either the backend was reimplemented
    # (and this script's PERMITTED set is now stale) or the vendored library is
    # dead weight. Either way somebody should say so deliberately.
    for permitted in sorted(PERMITTED):
        if not (ROOT / permitted).is_file():
            findings.append(f"{permitted}: permitted translation unit does not exist")
        elif permitted not in including:
            findings.append(
                f"{permitted}: does not include detria.hpp, but is the only file "
                f"permitted to -- update PERMITTED in this script, or drop lib/detria"
            )

    if findings:
        print(f"detria boundary check FAILED ({len(findings)}):\n", file=sys.stderr)
        for line in findings:
            print(f"  {line}", file=sys.stderr)
        print(
            "\ndetria.hpp is an implementation detail of one translation unit.\n"
            "See lib/detria/README.md and include/terrain/predicates/detria_exact.hpp.",
            file=sys.stderr,
        )
        return 1

    print(
        f"detria boundary check OK: {scanned} C/C++ files scanned, "
        f"detria.hpp included by {', '.join(sorted(including))} and nothing else."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
