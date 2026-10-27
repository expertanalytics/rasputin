# CLAUDE.md

The main session also loads the `dispatch` skill before its first act.

## 1. Personas
Each task goes to its persona in `.claude/agents/`; the file's `description`
says what it does.

## 2. Core Constraints & Technical Mandates
* **Size limit:** under **700 net production lines per pull request**, as
  `python3 tools/count_loc.py <base> <head>` counts them; its docstring
  defines a counted line. Code packed by hand under `# fmt: skip` /
  `# fmt: off` stays readable, and the review says why each new region is
  packed.
* **Rule text pays for itself:** a change that adds words to the files
  `tools/rule_sizes.py` counts cuts as many in the same change, unless Ola
  allows the growth.
* **Prohibited dependencies:** `CGAL`, `GDAL`, `OGR`, `Fiona`, `Rasterio`,
  and external `date` libraries (`tools/check_prohibited_deps.py`).
* **Core stack:** C++20 (Concepts, Pybind11, `std::chrono`); Python 3.12+
  (asyncio, Pydantic V2, Typer, Shapely, PyProj, NumPy, tifffile).
* **I/O boundary:** file decoding is Python's. The C++ core never opens a
  file, sees a path, links a codec or receives a CRS.

## 3. Test-Driven Development (TDD) Protocol
Before acting, read `.claude/REQUIRED-READING.md` and
`docs/increments/README.md`.
* **Step order.** Design (`@architect`), failing tests (`@tester`), code
  (`@developer`), test run, review (`@reviewer`), and `@perf`'s acceptance
  run when refine or mesh code is touched. A performance fix is timed by
  `@perf` before review. Failing tests go back to `@developer`.

## 4. Operational Commands
```bash
pytest
cmake -S . -B build && cmake --build build -j   # -j alone: nproc is Linux-only
ctest --test-dir build
mypy && ruff check . && ruff format --check .
python3 tools/check_code_history.py   # not run by the hook below
```
The C++ targets build with `-Wall -Wextra -Wpedantic -Werror`: any warning
fails. After each commit a hook runs ruff and the other three
`tools/check_*.py` gates, printing their output. CI
(`.github/workflows/main.yaml`) is authoritative: nothing is merge-ready
until `gh pr checks <pr>` is green.
