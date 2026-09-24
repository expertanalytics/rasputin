"""Assert the project's defining constraint: no CGAL, no GDAL, no wrappers.

CLAUDE.md section 2 prohibits CGAL, GDAL, OGR, Fiona, Rasterio (it wraps GDAL)
and external date libraries. That prohibition is the reason this project exists,
and until now it was enforced only by agent vigilance -- a reviewer noticing.
Reviewers are fallible and this is trivially machine-checkable.

The check reads declarations, not prose. `pyproject.toml` argues at length about
why Rasterio is excluded, and the raster docs name GDAL to explain what is being
avoided; a naive grep flags those and gets switched off. So this inspects Python
import statements, C/C++ #include directives, and the dependency arrays in
pyproject.toml -- the places where a dependency is actually contracted.

legacy/ is exempt: it is the archived CGAL implementation, kept deliberately as
the porting reference, and legacy/bindings.cpp is the CGAL original.
"""

from __future__ import annotations

import ast
import re
import sys
import tomllib
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).resolve().parent.parent

# Directories that must stay clean. legacy/ is deliberately absent.
SOURCE_DIRS = ["include", "src", "src_python", "bindings", "tests", "tools"]

# Every prohibition carries an AUTHORITY, and the check refuses to run without
# one, because adding a key here otherwise costs one line and no evidence.
# RULED quotes the human who asked, with a date. SPELLING means "the same
# library under another name" and points at the key it derives from, which must
# itself resolve to RULED. A key that is neither does not pass.
RULED = "ruled"      #: a human asked for it; the value quotes them, with a date
SPELLING = "spelling"  #: another name for a RULED key; the value is that key
PENDING = "pending"  #: nobody has ruled; the value says what is unresolved


class Rule(NamedTuple):
    """One prohibited token: why, who said so, and which CLAUDE.md name it is."""

    reason: str
    kind: str
    authority: str
    family: str  #: the name this appears under in CLAUDE.md section 2


_CGAL = "the user, 2026-09-16: 'strict zero-GDAL/CGAL dependencies'"
_GDAL = _CGAL
_RASTERIO = "the user, 2026-09-16: 'Let's drop rasterio first.'"

# Directories that must stay clean. legacy/ is deliberately absent.
SOURCE_DIRS = ["include", "src", "src_python", "bindings", "tests", "tools"]

PROHIBITED = {
    "cgal": Rule(
        "CGAL: GPL-encumbered; the migration replaces it with an MIT-licensed CDT",
        RULED, _CGAL, "CGAL"),
    "gdal": Rule("GDAL: prohibited by CLAUDE.md section 2", RULED, _GDAL, "GDAL"),
    "ogr": Rule("OGR: part of GDAL", SPELLING, "gdal", "OGR"),
    "fiona": Rule("Fiona: wraps GDAL", SPELLING, "gdal", "Fiona"),
    "rasterio": Rule("Rasterio: wraps GDAL", RULED, _RASTERIO, "Rasterio"),
    "osgeo": Rule("osgeo: the GDAL Python bindings", SPELLING, "gdal", "GDAL"),
    "gdalwarper.h": Rule(
        "GDAL: the warper header does not start with a bare 'gdal'",
        SPELLING, "gdal", "GDAL"),
    "cpl": Rule(
        "CPL: GDAL's portability layer (cpl_conv.h, cpl_string.h, ...)",
        SPELLING, "gdal", "GDAL"),
    "date/date.h": Rule(
        "external date library: superseded by C++20 <chrono>",
        PENDING,
        "deleting the vendored lib/date copy was asked for (2026-09-16); "
        "prohibiting the class was not",
        "external date libraries"),
    "date/tz.h": Rule(
        "external date library: superseded by C++20 <chrono>",
        SPELLING, "date/date.h", "external date libraries"),
}

INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]')


# `_` is a boundary alongside `/` and `.` because GDAL's C++ headers are spelled
# that way and nothing else is: gdal_priv.h is the ordinary entry point, and
# ogr_spatialref.h is how a CRS would arrive. Measured before this was widened --
# a planted #include <gdal_priv.h> passed the gate while the find_package(GDAL)
# on the next line was caught, so the include layer was enforcing nothing that
# the build layer was not already catching.
_BOUNDARY = ("/", ".", "_")


def offence(token: str) -> str | None:
    """Reason this import/include is prohibited, or None if it is fine."""
    low = token.lower()
    for bad, rule in PROHIBITED.items():
        if low == bad or any(low.startswith(f"{bad}{sep}") for sep in _BOUNDARY):
            return rule.reason
    return None


def check_python(path: Path) -> list[str]:
    findings = []
    try:
        tree = ast.parse(path.read_text(), str(path))
    except SyntaxError:
        return []  # Syntax is another gate's problem.
    for node in ast.walk(tree):
        names: list[tuple[str, int]] = []
        if isinstance(node, ast.Import):
            names = [(a.name, node.lineno) for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            names = [(node.module, node.lineno)]
        for name, lineno in names:
            if reason := offence(name.split(".")[0]):
                findings.append(f"{path.relative_to(ROOT)}:{lineno}: imports {name} -- {reason}")
    return findings


def check_cxx(path: Path) -> list[str]:
    findings = []
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        m = INCLUDE_RE.match(line)
        if m and (reason := offence(m.group(1))):
            findings.append(f"{path.relative_to(ROOT)}:{lineno}: includes {m.group(1)} -- {reason}")
    return findings


# The build system is where a prohibited dependency actually enters a C++ project:
# an #include only compiles because something put the library on the include path.
# Scanning headers alone would pass a find_package(CGAL) or an apt install
# libgdal-dev, which is the gate's whole purpose defeated one layer down.
def build_files() -> list[Path]:
    """Every CMakeLists and workflow in the tree, excluding the exempt archive.

    Globbed rather than listed: project_structure.md already plans per-module
    CMakeLists under src/, and a hardcoded list reports OK on the file it does
    not know about.
    """
    found = [p for p in ROOT.rglob("CMakeLists.txt")]
    found += [p for p in (ROOT / ".github" / "workflows").glob("*.y*ml")]
    skip = ("legacy", "build", "lib", ".venv")
    return sorted(
        p for p in found
        if not any(part.startswith(skip) or part in skip for part in p.relative_to(ROOT).parts)
    )

BUILD_RE = re.compile(
    r"(find_package|target_link_libraries|link_libraries|find_library|FetchContent_Declare"
    r"|apt-get install|apt install|brew install|vcpkg install|conan install)\b(?P<rest>[^\n]*)",
    re.IGNORECASE,
)


# Build directives name a package, not a header path, so the include-layer keys
# (which are path-shaped: "date/date.h") never match a bare CMake token. These
# are the build-layer spellings.
BUILD_PROHIBITED = {
    "date": Rule("external date library: superseded by C++20 <chrono>",
                 SPELLING, "date/date.h", "external date libraries"),
}


def build_offence(token: str) -> str | None:
    """Reason this build-directive token is prohibited, or None."""
    if reason := offence(token):
        return reason
    rule = BUILD_PROHIBITED.get(token.lower())
    return rule.reason if rule else None


def check_build_file(path: Path) -> list[str]:
    """Prohibited names appearing in dependency-acquiring build directives."""
    findings = []
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        m = BUILD_RE.search(line)
        if not m:
            continue
        for word in re.split(r"[^A-Za-z0-9_.:+-]+", m.group("rest")):
            # Normalise the spellings a package or CMake target actually takes:
            # GDAL::GDAL, libgdal-dev, gdal-devel, CGAL_ROOT. A gate that only
            # matched the bare name would pass `apt install libgdal-dev`, which
            # is precisely how the dependency would arrive.
            bare = word.split("::")[0].split("-")[0].removesuffix("_ROOT")
            if bare.lower().startswith("lib"):
                bare = bare[3:]
            if bare and (reason := build_offence(bare)):
                findings.append(
                    f"{path.relative_to(ROOT)}:{lineno}: build directive names {word} -- {reason}"
                )
    return findings


def check_pyproject() -> list[str]:
    """Declared dependencies only -- the prose around them is allowed to say 'GDAL'."""
    path = ROOT / "pyproject.toml"
    data = tomllib.loads(path.read_text())
    project = data.get("project", {})
    declared: list[str] = list(project.get("dependencies", []))
    for extra in project.get("optional-dependencies", {}).values():
        declared += extra
    declared += data.get("build-system", {}).get("requires", [])

    findings = []
    for spec in declared:
        name = re.split(r"[<>=!~\[; ]", spec.strip(), maxsplit=1)[0]
        if reason := offence(name):
            findings.append(f"pyproject.toml: declares dependency {spec!r} -- {reason}")
    return findings


#: CLAUDE.md section 2's prose list. Parsed rather than hardcoded, because the
#: two drifting apart is how "osgeo" and later "boost" entered the key set
#: without ever appearing in the prose a human reads.
CLAUDE_LIST_RE = re.compile(
    r"\*\*Prohibited Dependencies:\*\*\s*Never introduce\s+(?P<list>.+?)\.\s*\n",
    re.DOTALL,
)


def prose_families() -> set[str]:
    """The dependency names CLAUDE.md section 2 actually spells out."""
    text = (ROOT / "CLAUDE.md").read_text()
    m = CLAUDE_LIST_RE.search(text)
    if not m:
        return set()
    # Strip backticks, parentheticals like "(it wraps GDAL)", and list glue.
    raw = re.sub(r"\([^)]*\)", "", m.group("list")).replace("`", "")
    parts = re.split(r",|\bor\b", raw.replace("\n", " "))
    return {" ".join(part.split()) for part in parts if part.strip()}


def check_authorities() -> list[str]:
    """Every prohibition names who asked for it, and the prose agrees.

    Two failures, kept separate because they have different fixes. A key with no
    resolvable authority is a rule nobody asked for and must be removed or
    ruled on. A family in one list and not the other is drift: the gate quietly
    enforcing something the prose never told a human about, or the prose
    claiming something the gate does not check.
    """
    findings = []
    everything = {**PROHIBITED, **BUILD_PROHIBITED}

    for token, rule in sorted(everything.items()):
        if not rule.authority.strip():
            findings.append(f"key {token!r}: empty authority -- who asked for this?")
        elif rule.kind == SPELLING:
            # Follow the chain to a RULED root, refusing a cycle. A SPELLING
            # pointing at another SPELLING is fine; one pointing at a PENDING
            # inherits the PENDING, which is why date/tz.h is not reported
            # separately from date/date.h.
            seen, cursor = {token}, rule
            while cursor.kind == SPELLING:
                parent = everything.get(cursor.authority)
                if parent is None:
                    findings.append(
                        f"key {token!r}: derives from {cursor.authority!r}, which is not a key")
                    break
                if cursor.authority in seen:
                    findings.append(f"key {token!r}: authority chain is a cycle")
                    break
                seen.add(cursor.authority)
                cursor = parent
            else:
                if cursor.kind == PENDING:
                    findings.append(
                        f"key {token!r}: derives from {cursor.authority[:40]}..., unresolved")
        elif rule.kind == PENDING:
            findings.append(f"key {token!r}: NOT RULED ON -- {rule.authority}")

    prose, gated = prose_families(), {r.family for r in everything.values()}
    if not prose:
        findings.append("CLAUDE.md section 2's prohibited list could not be parsed")
    else:
        for missing in sorted(gated - prose):
            findings.append(f"family {missing!r} is gated but absent from CLAUDE.md section 2")
        for missing in sorted(prose - gated):
            findings.append(f"CLAUDE.md section 2 names {missing!r}, which no key enforces")
    return findings


def main() -> int:
    findings: list[str] = []
    scanned = 0

    for directory in SOURCE_DIRS:
        base = ROOT / directory
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file():
                continue
            if path.suffix == ".py":
                scanned += 1
                findings += check_python(path)
            elif path.suffix in {".cpp", ".hpp", ".h", ".cc", ".cxx"}:
                scanned += 1
                findings += check_cxx(path)

    for path in build_files():
        scanned += 1
        findings += check_build_file(path)

    findings += check_pyproject()

    # Run last so it reports below the dependency findings, but it gates the
    # same exit code: a prohibition nobody authorised is a defect in this file,
    # not a lesser kind of problem than an import that violates one.
    governance = check_authorities()

    if findings or governance:
        total = len(findings) + len(governance)
        print(f"Prohibited dependency check FAILED ({total}):\n", file=sys.stderr)
        for line in findings:
            print(f"  {line}", file=sys.stderr)
        if governance:
            print("\n  -- prohibitions with no author, or out of step with the prose --",
                  file=sys.stderr)
            for line in governance:
                print(f"  {line}", file=sys.stderr)
        print(
            "\nlegacy/ is exempt; the new tree is not. A governance finding above is "
            "not fixed\nby editing this file: it is a question for the user, and the "
            "answer goes in\nCLAUDE.md section 2 and in the key's authority together.",
            file=sys.stderr,
        )
        return 1

    print(f"Prohibited dependency check OK: {scanned} source and build files "
          f"plus pyproject.toml are clean, and all {len(PROHIBITED)} prohibitions "
          "name who asked for them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
