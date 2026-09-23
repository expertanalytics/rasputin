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

ROOT = Path(__file__).resolve().parent.parent

# Directories that must stay clean. legacy/ is deliberately absent.
SOURCE_DIRS = ["include", "src", "src_python", "bindings", "tests", "tools"]

# Import/include tokens -> why they are banned.
PROHIBITED = {
    "cgal": "CGAL: GPL-encumbered; the migration replaces it with an MIT-licensed CDT",
    "gdal": "GDAL: prohibited by CLAUDE.md section 2",
    "ogr": "OGR: part of GDAL",
    "fiona": "Fiona: wraps GDAL",
    "rasterio": "Rasterio: wraps GDAL",
    "osgeo": "osgeo: the GDAL Python bindings",
    "gdalwarper.h": "GDAL: the warper header does not start with a bare 'gdal'",
    "cpl": "CPL: GDAL's portability layer (cpl_conv.h, cpl_string.h, ...)",
    "date/date.h": "external date library: superseded by C++20 <chrono>",
    "date/tz.h": "external date library: superseded by C++20 <chrono>",
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
    for bad, reason in PROHIBITED.items():
        if low == bad or any(low.startswith(f"{bad}{sep}") for sep in _BOUNDARY):
            return reason
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
    "date": "external date library: superseded by C++20 <chrono>",
}


def build_offence(token: str) -> str | None:
    """Reason this build-directive token is prohibited, or None."""
    if reason := offence(token):
        return reason
    return BUILD_PROHIBITED.get(token.lower())


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

    if findings:
        print(f"Prohibited dependency check FAILED ({len(findings)}):\n", file=sys.stderr)
        for line in findings:
            print(f"  {line}", file=sys.stderr)
        print(
            "\nCLAUDE.md section 2 prohibits CGAL, GDAL, OGR, Fiona, Rasterio and\n"
            "external date libraries. legacy/ is exempt; the new tree is not.",
            file=sys.stderr,
        )
        return 1

    print(f"Prohibited dependency check OK: {scanned} source and build files "
          "plus pyproject.toml are clean.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
