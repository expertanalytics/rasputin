---
name: python-development
description: Idiomatic Python 3.11+ for the orchestration layer. Use when writing Pydantic V2 models, Typer CLIs, asyncio pipelines, pytest/pytest-asyncio suites, or satisfying the strict mypy and ruff gates over src_python/tin_engine.
---

# Agent Skill: Idiomatic Python & Async CLI Generation

Enforces cutting-edge, safe, and testable Python 3.11+ architecture for orchestration, DTM generation pipelines, and CLI tools.

## 1. Data Integrity & Validation (Pydantic V2)
* **Strict Modeling:** Every external data input, configuration option, and API boundary must use **Pydantic models** (`pydantic.BaseModel`).
* **Settings Management:** Use `pydantic-settings` for application and environment configuration.
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

## 4. Testability & Quality Assurance
* **Test Architecture:** Maintain 100% testable logic. Use `pytest` for all unit and integration/system-level tests.
* **Unit Level:** Mock the Pybind11 geometry engine using standard fixtures to isolate pure Python logic, DTM parsing, and CLI formatting.
* **System Level:** Implement complete integration tests that orchestrate the pipeline from a mock DTM file input to a generated TIN output.
* **Async Testing:** Use `pytest-asyncio` for verifying asynchronous execution loops and non-blocking tasks.

## 5. Pragmatic Code Design
* **Pragmatic Polymorphism:** Prefer composition over inheritance. Use structural subtyping (`typing.Protocol`) or simple functional abstractions instead of heavy object-oriented hierarchies.
* **Change Limit:** Strict maximum of **700 LOC** per pull request (excluding tests). Keep code blocks punchy, clean, and documentation-driven.
