---
name: python-development
description: Idiomatic Python 3.12+ for the orchestration layer. Use when writing Pydantic V2 models, Typer CLIs, asyncio pipelines, pytest/pytest-asyncio suites, or satisfying the strict mypy and ruff gates over src_python/tin_engine.
---

# Agent Skill: Idiomatic Python & Async CLI Generation

Enforces cutting-edge, safe, and testable Python 3.12+ architecture for orchestration, DTM generation pipelines, and CLI tools.

## 1. Data Integrity & Validation (Pydantic V2)
* **Strict Modeling:** Every external data input, configuration option, and API boundary must use **Pydantic models** (`pydantic.BaseModel`).
* **Type Safety:** Enforce strict runtime and static type hinting. Leverage `Annotated` types for precise validation criteria (e.g., coordinate boundaries, file path validation).

## 2. Asynchronous & Decoupled Architecture
* **Async by Default:** Design I/O and orchestration loops as `async` functions using `asyncio`. Prepare all execution paths to run inside external loops (e.g., GUI backends or FastAPI services).
* **Non-Blocking C++ Invocation:** When calling long-running Pybind11 meshing methods, run them via `asyncio.to_thread()` to prevent blocking the main Python event loop (assuming the C++ layer releases the GIL).
* **Frontend-Backend Separation:** Keep the business logic completely decoupled from UI or CLI elements. Expose generic, stable APIs that can be wrapped by either a command-line interface or a WebSocket/REST API.

## 3. CLI Engineering & Security
* **Modern CLI:** Use `typer` or `click` to build user-friendly, structured CLI applications for generating TIN models.
* **Security & Input Validation:**
  * Validate all file paths before ingestion; reject symlink exploits and path traversal attacks (use `pathlib.Path.resolve()`).
  * Never use `eval()`, `exec()`, or insecure YAML/Pickle parsing. Use `pydantic` or `json` for serialization.
  * Sanitise input coordinates and parameters at the CLI boundary before passing them to Pybind11.

## 5. Pragmatic Code Design
* **Pragmatic Polymorphism:** Prefer composition over inheritance. Use structural subtyping (`typing.Protocol`) or simple functional abstractions instead of heavy object-oriented hierarchies.

## 6. Docstrings and Comments
Docstrings and comments describe the code as it is, for a reader who has never
seen the project's history.
* A docstring says what the object does: what it takes, what it returns, what
  it refuses, and the reason for anything a careful reader would otherwise
  question. A comment says why a line is as it is, when the code cannot say it.
* They never carry history: no increment numbers or labels (`increment 29`,
  `16c`, `23a-2`), no IDs from a design file (`R5`, `D2`, `Q11`), no review or
  PR names (`audit PR C`, `round 3`), no record of who decided, no commit
  hashes, no account of what the code used to do. Git and `docs/increments/`
  keep that.
* A reason that came from a decision is written as the reason, never as the
  decision's label. A value that was chosen rather than derived says so: "a
  chosen threshold: below 1 mm two tiles are taken to agree".
* A pointer to a design file is kept only where the derivation or measurement
  behind a value is too long to repeat. It names the file and a heading, in
  the form `tools/check_citations.py` verifies (`` `<path>.md`, *<heading>* ``),
  on one line. At most one per module docstring; none in classes, functions,
  comments or help text.
* The same holds for the help text the CLI prints (`help=` strings).
* A link to an outside bug report (an issue in NumPy, PROJ, GitHub) is allowed
  where the code works around that bug, and only with the reason written
  beside it, so the comment still reads if the link dies.
* Tool directives (`# type: ignore[...]`, `# noqa: ...`, `# ruff: noqa`,
  `# mypy: ...`, `# fmt: skip`) are code: keep them exactly as they are.

Checked by `tools/check_code_history.py`, over the paths listed under
`[tool.check_code_history]` in `pyproject.toml`. The reviewer judges the rest:
whether a sentence still describes the code, and history written without a
label ("used to", a bare increment number).
