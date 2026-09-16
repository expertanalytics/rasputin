# CLAUDE.md - Core Agent & Orchestration Guide

This file drives the automation loops for Claude Code (`@orchestrator`) in this repository.

## 1. Multi-Agent Ecosystem
This project is governed by specialized sub-agents. Always defer tasks to the correct persona in `.claude/agents/`:
* `@orchestrator`: Master project driver. Handles human input and chains the TDD loop.
* `@architect`: Enforces declarative structures, component boundaries, and interface decoupling.
* `@migration-expert`: Porting lead. Deconstructs legacy logic into the new target architecture.
* `@tester`: Owns the test suites. Enforces the 85% coverage floor (see `testing.md`) and adversarial geometry fuzzing.
* `@developer`: Writes clean, high-performance C++20 and async Python code.
* `@reviewer`: Final gatekeeper. Audits LOC, code quality, readability, and documentation.

## 2. Core Constraints & Technical Mandates
* **Strict Size Limit:** Production code changes must be strictly under **700 LOC** per interaction.
* **Prohibited Dependencies:** Never introduce `CGAL`, `GDAL`, `OGR`, `Fiona`, `Rasterio` (it wraps GDAL), or external `date` libraries.
* **Core Stack:** Modern C++ (C++20 Concepts, Pybind11, `std::chrono`) + Async Python 3.11+ (Pydantic V2, Typer, Shapely, PyProj, NumPy, tifffile).
* **I/O Boundary:** File decoding is Python's. The C++ core never opens a file, sees a path, or links a codec, and CRS never crosses into it. See the `raster` section of `project_structure.md`.

## 3. Test-Driven Development (TDD) Protocol
Every code alteration or legacy migration step must execute this strict pipeline via `@orchestrator`:
1. `@architect` or `@migration-expert` defines interfaces and types.
2. `@tester` writes failing unit/async test cases *first* (including happy path and edge cases).
3. `@developer` writes the minimal code needed to pass the active tests.
4. `@tester` and `@reviewer` validate results and type consistency before merge readiness.

## 4. Operational Commands

### Python Layer (Async & CLI)
```bash
pytest                 # Run all modern tests (uses pytest-asyncio)
pytest tests/python/   # Run target Python testing suite
```

### C++ Core Layer (C++20 Concepts)
```bash
cmake -S . -B build && cmake --build build -j   # -j alone: nproc is Linux-only
ctest --test-dir build                          # runs test_point, test_raster, test_solar_position
```

### Static gates (Python)
```bash
mypy                   # strict, over src_python/tin_engine
ruff check .           # legacy/ is excluded
```

### Governance gates
```bash
python tools/check_prohibited_deps.py   # section 2, checked against real imports
python tools/check_legacy_imports.py    # legacy/ must stay self-consistent
```

### CI
The gates above also run in GitHub Actions (`.github/workflows/main.yaml`), and
CI is authoritative: local green does not mean the branch is green. A whole
branch once merged with CI red on every commit because the workflow was never
checked, so **verify check status before declaring anything merge-ready**:
```bash
gh pr checks <pr>      # must be green; required for merge on master
```

