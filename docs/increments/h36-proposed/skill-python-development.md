---
name: python-development
description: Idiomatic Python 3.12+ for the orchestration layer. Use when writing Pydantic V2 models, Typer CLIs, asyncio pipelines, pytest/pytest-asyncio suites, or satisfying the strict mypy and ruff gates over src_python/tin_engine.
---

# Python Orchestration Layer

## 1. Models
External input, configuration and API boundaries are strictly typed
Pydantic V2 models; `Annotated` carries validation bounds.

## 2. Async and decoupled
I/O and orchestration are `async` and run inside an outside event loop; long
Pybind11 calls go through `asyncio.to_thread()`. Business logic knows
nothing of the CLI or UI that wraps it.

## 3. CLI and safety
`typer` for the CLI. Resolve and check every path before reading. No `eval`,
`exec`, pickle or unsafe YAML. Sanitise coordinates and parameters at the
CLI boundary.

## 5. Design
Composition over inheritance; `typing.Protocol` or functions over class
hierarchies.

## 6. Docstrings and Comments
They describe the code as it is, for a reader who never saw its history.
* A docstring says what the object takes, returns and refuses, and why
  anything surprising is so; a comment says why a line is as it is.
* No history: no increment or design labels, review names, deciders, commit
  hashes or past behaviour. A decision is written as its reason; a chosen
  value says it was chosen.
* A design-file pointer only where a derivation is too long to repeat:
  `` `<path>.md`, *<heading>* ``, one line, at most one per module
  docstring.
* The same for CLI `help=` text.
* An outside bug link only beside its workaround, with the reason written.
* Tool directives (`# type: ignore[...]`, `# noqa`, `# fmt: skip`) stay.

`tools/check_code_history.py` checks the labels; the reviewer the rest.
